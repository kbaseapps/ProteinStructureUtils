import hashlib
import logging
import os
from readline import read_history_file
import sys
import re
import gzip
import string
import shutil
import uuid
import errno
import pathlib
import subprocess
from urllib.parse import urlparse
from copy import deepcopy

import json, requests
import aiohttp, asyncio
#import modelseedpy
from requests.exceptions import ConnectionError, HTTPError, RequestException
from python_graphql_client import GraphqlClient

from installed_clients.AbstractHandleClient import AbstractHandle
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace
from installed_clients.baseclient import ServerError as WorkspaceError

from ProteinStructureUtils.Utils.PDBUtil import PDBUtil


class RCSBUtil:

    # Set thresholds to restrict which alignments will be significant or sequence identity matches
    EVALUE_CUTOFF = 1e-5
    IDENTITY_CUTOFF = 0.75
    LOGICAL_AND = 0

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.token = config['KB_AUTH_TOKEN']
        self.user_id = config['USER_ID']
        self.dfu = DataFileUtil(self.callback_url)
        self.hs = AbstractHandle(config['handle-service-url'])
        self.ws_client = Workspace(config['workspace-url'])
        self.shock_url = config['shock-url']
        self.pdb_util = PDBUtil(config)
        self.inchK_cpd_jsonObj = {}

        self.__baseDownloadUrl = 'https://files.rcsb.org/download'
        self.__baseSearchUrl = 'https://search.rcsb.org/rcsbsearch/v2/query'
        self.__baseGraphqlUrl = 'https://data.rcsb.org/graphql'
        self.__graphqlClient = GraphqlClient(endpoint=self.__baseGraphqlUrl)
        self.__graphqlQueryTemplate = """{
            entries(entry_ids:["%s"]) {
                rcsb_id
                exptl {
                    method
                }
                rcsb_primary_citation {
                    title
                    rcsb_authors
                    journal_abbrev
                    year
                    pdbx_database_id_PubMed
                }
                polymer_entities {
                    rcsb_id
                    entity_poly {
                        pdbx_seq_one_letter_code
                        pdbx_strand_id
                    }
                    rcsb_entity_source_organism {
                        ncbi_taxonomy_id
                        ncbi_scientific_name
                    }
                    rcsb_polymer_entity_container_identifiers {
                        reference_sequence_identifiers {
                            database_accession
                            database_name
                        }
                    }
                    rcsb_polymer_entity {
                        rcsb_ec_lineage {
                            id
                        }
                    }
                    uniprots {
                      rcsb_uniprot_protein {
                        name {
                          value
                        }
                        ec {
                          number
                          provenance_code
                        }
                      }
                    }
                }
                nonpolymer_entities {
                  nonpolymer_comp {
                    rcsb_chem_comp_descriptor {
                      InChI
                      InChIKey
                      SMILES
                    }
                  }
                }
            }
        }
        """

    def _validate_rcsb_seqquery_params(self, params):
        """
            _validate_rcsb_seqquery_params:
                validates input params to query_structure_anno
        """
        # check for required parameters
        for p in ['workspace_name', 'sequence_strings']:
            if p not in params:
                raise ValueError(f'Parameter "{p}" is required, but missing!')

        if params.get('evalue_cutoff', None):
            self.EVALUE_CUTOFF = float(params['evalue_cutoff'])

        if params.get('identity_cutoff', None):
            self.IDENTITY_CUTOFF = params['identity_cutoff']

        return params

    def _validate_rcsb_query_params(self, params):
        """
            _validate_rcsb_query_params:
                validates input params to query_structure_info
        """
        # check for required parameters
        for p in ['workspace_name']:
            if p not in params:
                raise ValueError(f'Parameter "{p}" is required, but missing!')

        queriables = {'sequence_strings': 'sequence',
                      'uniprot_ids': 'uniprot_id',
                      'ec_numbers': 'ec_number',
                      'inchis': 'InChI',
                      'smiles': 'SMILES'}

        if params.get('evalue_cutoff', None):
            self.EVALUE_CUTOFF = float(params['evalue_cutoff'])

        if params.get('identity_cutoff', None):
            self.IDENTITY_CUTOFF = params['identity_cutoff']

        if params.get('logical_and', None):
            self.LOGICAL_AND = params['logical_and']

        # check for queriable parameters
        inputJsonObj = {}
        for pk in queriables:
            if params.get(pk, None):
                inputJsonObj[queriables[pk]] = params[pk]

        if not inputJsonObj:
            raise ValueError('At least one (1) queriable parameter should be specified!')

        return inputJsonObj, params.get('workspace_name')

    def _create_seq_ec_uniprot_params(self, searchType, valList):
        """
            _create_seq_ec_uniprot_params - Based on sequence_strings/uniprot_ids/ec_numbers input,
                                            create query parameters in Json object format.
                                            return the seqEcUniprotJsonObject
        """
        seqEcUniprotJsonObject = {}
        if searchType not in ('sequence', 'sequence_strings', 'uniprot_id',
                              'uniprot_ids', 'ec_number', 'ec_numbers'):
            return seqEcUniprotJsonObject
        #
        if searchType == "sequence_strings" or searchType == "sequence":
            if len(valList) == 1:
                seqEcUniprotJsonObject = {
                    "query": {
                        "type": "terminal",
                        "service": "sequence",
                        "parameters": {
                            "evalue_cutoff": self.EVALUE_CUTOFF,
                            "identity_cutoff": self.IDENTITY_CUTOFF,
                            "target": "pdb_protein_sequence",
                            "value": valList[0]
                        }
                    },
                    "return_type": "entry",
                    "request_options": {
                        "return_all_hits": True
                    }
                }
            elif len(valList) > 1:
                nodeList = []
                for value in valList:
                    nodeList.append({"type": "terminal",
                                     "service": "sequence",
                                     "parameters": {
                                        "evalue_cutoff": self.EVALUE_CUTOFF,
                                        "identity_cutoff": self.IDENTITY_CUTOFF,
                                        "target": "pdb_protein_sequence",
                                        "value": value
                                     }})

                seqEcUniprotJsonObject = {
                    "query": {
                         "type": "group",
                         "logical_operator": "or",
                         "nodes": nodeList,
                         "label": "text"
                    },
                    "return_type": "entry",
                    "request_options": {
                        "return_all_hits": True
                    }
                }
            #
        elif searchType == "uniprot_id" or searchType == 'uniprot_ids':
            seqEcUniprotJsonObject = {
                "query": {
                    "type": "group",
                    "logical_operator": "and",
                    "nodes": [
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {
                                "attribute": ("rcsb_polymer_entity_container_identifiers."
                                              "reference_sequence_identifiers.database_accession"),
                                "operator": "in",
                                "negation": False,
                                "value": valList
                            }
                        },
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {
                                "attribute": ("rcsb_polymer_entity_container_identifiers."
                                              "reference_sequence_identifiers.database_name"),
                                "operator": "exact_match",
                                "value": "UniProt",
                                "negation": False
                            }
                        }
                    ],
                    "label": "nested-attribute"
                },
                "return_type": "entry",
                "request_options": {
                    "return_all_hits": True
                }
            }
        elif searchType == "ec_number" or searchType == 'ec_numbers':
            seqEcUniprotJsonObject = {
                "query": {
                    "type": "terminal",
                    "label": "text",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity.rcsb_ec_lineage.id",
                        "operator": "in",
                        "negation": False,
                        "value": valList
                    }
                },
                "return_type": "entry",
                "request_options": {
                    "return_all_hits": True
                }
            }

        return seqEcUniprotJsonObject

    def _create_chem_params(self, searchType, valList):
        """
            _create_chem_params - Based on inchis/smiles input data, define query parameters and
                                  create it in Json object format.
                                  return the chemJsonObject
        """
        chemJsonObject = {}
        if searchType not in ("InChI", "SMILES", "InChIKey", 'inchis', 'smiles'):
            return chemJsonObject
        #
        if len(valList) == 1:
            chemJsonObject = {
                "query": {
                    "type": "terminal",
                    "service": "chemical",
                    "parameters": {
                        "value": valList[0],
                        "type": "descriptor",
                        "descriptor_type": searchType,
                        "match_type": "graph-exact"
                    }
                },
                "return_type": "entry",
                "request_options": {
                    "results_content_type": [
                        "computational",
                        "experimental"
                    ],
                    "return_all_hits": True
                }
            }
        elif len(valList) > 1:
            nodeList = []
            for value in valList:
                nodeList.append({"type": "terminal",
                                 "service": "chemical",
                                 "parameters": {
                                    "value": value,
                                    "type": "descriptor",
                                    "descriptor_type": searchType,
                                    "match_type": "graph-exact"
                                 }})

            chemJsonObject = {
                "query": {
                    "type": "group",
                    "logical_operator": "or",
                    "nodes": nodeList,
                    "label": "text"
                },
                "return_type": "entry",
                "request_options": {
                    "results_content_type": [
                        "computational",
                        "experimental"
                    ],
                    "return_all_hits": True
                }
            }
        #
        return chemJsonObject

    def _queryRCSB(self, jsonQueryObj):
        """ Run search API to get RCSB return as a json object
            The return json object should look like:
            {
               "query_id" : "4aab471d-0134-45d8-b777-11d6e3b5645c",
               "result_type" : "entry",
               "total_count" : 57742,
               "result_set" : [
                  { "identifier" : "10GS",
                    "score" : 1.0
                  },
                  { "identifier" : "11GS",
                    "score" : 1.0
                  },
                  ......
               ]
            }
        """
        try:
            reqH = requests.post(self.__baseSearchUrl, json=jsonQueryObj)
            reqH.raise_for_status()
            return json.loads(reqH.text)
        except (HTTPError, ConnectionError, RequestException) as e:
            logging.info(" _queryRCSB ERROR ".center(30, "-"))
            logging.info(f'Querying RCSB db with {jsonQueryObj} had an Error: {e}')
            return {}
        except Exception as e:
            print(e)
            return {}

    def _readRCSBResult(self, rcsb_retObj):
        """
            _readRCSBResult: parse the RCSB query result and return a PDB ID list & an id_score list
        """
        total_number = 0
        if rcsb_retObj and rcsb_retObj.get('total_count', None):
            total_number = rcsb_retObj["total_count"]
        else:
            return [], {}

        # Read entry ID from "result_set" list
        retIdList = []
        retIdScores = {}
        if "result_set" in rcsb_retObj:
            for obj in rcsb_retObj["result_set"]:
                if "identifier" in obj:
                    rcsb_id = obj['identifier']
                    retIdList.append(rcsb_id)
                    retIdScores[rcsb_id] = obj['score']
            if len(retIdList) > 1:
                retIdList.sort()

        if len(retIdList) != total_number:
            errMsg = (f"The total_count={total_number}, "
                      "but only got {len(retIdList)} id(s) from result_set.")
            raise ValueError(errMsg)

        return retIdList, retIdScores

    def _run_rcsb_search(self, jsonQueryObj={}):
        """ Run search API to return an object that has a PDB ID list like:
                [
                    "117E",
                    "11AS",
                    "12AS",
                  ......
                ]
            and an id_score dictionary like:
                {
                    "117E": "score": 1.0,
                    "11AS": "score": 1.0,
                    "12AS": "score": 1.0,
                  ......
                }
        """
        if not jsonQueryObj:
            return [], {}

        retJsonObj = self._queryRCSB(jsonQueryObj)
        if not retJsonObj:
            return [], {}
        return self._readRCSBResult(retJsonObj)

    def _queryGraphql(self, id_list=[]):
        """
            _queryGraphql: Given a list of rcsb_ids, Query the RCSB GraphQL API
        """
        logging.info(f'Querying GraphQL db for the structures of {len(id_list)} rcsd_ids')
        if not id_list:
            return {'data': None}
        queryString = self.__graphqlQueryTemplate % '", "'.join(id_list)
        try:
            #graphql_ret = self.__graphqlClient.execute(query=queryString)
            # asyncio.run() is available ONLY for python 3.7 or newer
            # graphql_ret = asyncio.run(self.__graphqlClient.execute_async(query=queryString))
            evt_loop = asyncio.get_event_loop()
            graphql_ret = evt_loop.run_until_complete(
                          self.__graphqlClient.execute_async(query=queryString))
            if 'errors' in graphql_ret:
                raise ConnectionError(graphql_ret['errors'][0]['message'])
            return graphql_ret
        except (aiohttp.ClientConnectionError, asyncio.TimeoutError):
            logging.info(" _queryGraphql ERROR ".center(30, "-"))
            logging.info(f'Connecting to RCSB GraphQL host at "{self.__baseGraphqlUrl} errored.')
            return {'data': None}
        except (HTTPError, ConnectionError, RequestException) as e:
            # not raising error to allow continue with other chunks
            logging.info(f'Querying RCSB GraphQL had a Connection Error:************\n {e}.\n'
                         'Or database connection request had no response!')
            return {'data': None}
        except (RuntimeError, TypeError, KeyError, ValueError) as e:
            err_msg = f'Querying RCSB errored with message: {e.message} and data: {e.data}'
            raise ValueError(err_msg)

    def _formatRCSBJson(self, entries):
        """
            _formatRCSBJson:Format rcsb GraphQL returned entries into one Json object per rcsb entry
                            Return an unspecifiedObject as the following:
                    dic[entry['rcsb_id']] = {
                        'method': [entry['exptl'][0]['method']],
                        'primary_citation': entry[prim_cite]<references>,
                        'polymer_entities': [{
                            'id': '1',
                            'one_letter_code_sequence': 'VNIKTNPFKAVSFVESAIKKALDNAGYLIAEI...',
                            'pdbx_strand_id': ['A'],
                            'source_organism': [{'ncbi_taxonomy_id': 10760, 'ncbi_scientific_name': 'Escherichia phage T7'}],
                            'taxonomy': [(10760, 'Escherichia phage T7')],
                            'uniprotID': ['P00969'],
                            'ref_sequence_ids': [{'database_accession': 'P00969', 'database_name': 'UniProt'}],
                            'uniprot_name': ['Hemoglobin subunit beta'],
                            'ec_numbers': ['6', '6.5', '6.5.1', '6.5.1.1'],
                            'uniprot_ec': []
                        },...],
                        'nonpolymer_entities': [{
                            'InChIKey': [...]
                        },...]
                    }
        """
        if not entries:
            return {}

        # short-naming the long rcsb data attribute strings
        src_organism = 'rcsb_entity_source_organism'
        polym_entity = 'rcsb_polymer_entity'
        ec_lineage = 'rcsb_ec_lineage'
        entity_container_ids = 'rcsb_polymer_entity_container_identifiers'
        ref_sequence_ids = 'reference_sequence_identifiers'
        uprot_prot = 'rcsb_uniprot_protein'
        chem_comp_desc = 'rcsb_chem_comp_descriptor'

        dic = {}
        for entry in entries:
            dic[entry['rcsb_id']] = {
                'primary_citation': entry.get('rcsb_primary_citation', ''),
                'polymer_entities': [],
                'nonpolymer_entities': []
            }
            methods = []
            for exp in entry.get('exptl', []):
                if exp.get('method', ''):
                    methods.append(exp['method'])
            dic[entry['rcsb_id']]['method'] = methods

            for pe in entry['polymer_entities']:
                pe_dic = {
                    'id': pe['rcsb_id'].split("_")[1],
                    'one_letter_code_sequence': pe['entity_poly']['pdbx_seq_one_letter_code'],
                    'pdbx_strand_id': pe['entity_poly']['pdbx_strand_id'].split(',')
                }
                if (pe.get(src_organism, None)):
                    pe_dic['source_organism'] = pe[src_organism]
                    pe_dic['taxonomy'] = []
                    for srco in pe[src_organism]:
                        pe_dic['taxonomy'].append(
                          (srco.get('ncbi_taxonomy_id', ''), srco.get('ncbi_scientific_name', '')))
                if (pe.get(entity_container_ids, None) and
                        pe[entity_container_ids].get(ref_sequence_ids, None)):
                    pe_dic['uniprotID'] = []
                    pe_dic['ref_sequence_ids'] = pe[entity_container_ids][ref_sequence_ids]
                    for rsid in pe_dic['ref_sequence_ids']:
                        pe_dic['uniprotID'].append(rsid.get('database_accession', ''))

                if (pe.get(polym_entity, None) and pe[polym_entity].get(ec_lineage, None)):
                    ec_numbers = []
                    for idic in pe[polym_entity][ec_lineage]:
                        if idic['id'] and (idic['id'] not in ec_numbers):
                            ec_numbers.append(idic['id'])

                    if ec_numbers:
                        pe_dic['ec_numbers'] = ec_numbers

                if pe.get('uniprots', None):
                    for unp in pe['uniprots']:
                        pe_dic['uniprot_name'] = []
                        pe_dic['uniprot_ec'] = []
                        if unp.get(uprot_prot, None):
                            uniprot_prot = unp[uprot_prot]
                            if uniprot_prot.get('name', None):
                                pe_dic['uniprot_name'].append(
                                    uniprot_prot['name'].get('value', ''))
                            if uniprot_prot.get('ec', None):
                                pe_dic['uniprot_ec'].extend(uniprot_prot['ec'])

                dic[entry['rcsb_id']]['polymer_entities'].append(pe_dic)

            dic[entry['rcsb_id']]['nonpolymer_entities'] = []
            if entry.get('nonpolymer_entities', None):
                descriptorDic = {}
                for npe in entry['nonpolymer_entities']:
                    if (npe.get('nonpolymer_comp', None) and
                       npe['nonpolymer_comp'].get(chem_comp_desc, None)):
                        # for desType in ('InChI', 'InChIKey', 'SMILES'):
                        if npe['nonpolymer_comp'][chem_comp_desc].get('InChIKey', None):
                            dt_val = npe['nonpolymer_comp'][chem_comp_desc]['InChIKey']
                            if 'InChIKey' in descriptorDic:
                                descriptorDic['InChIKey'].append(dt_val)
                            else:
                                descriptorDic['InChIKey'] = [dt_val]
                dic[entry['rcsb_id']]['nonpolymer_entities'].append(descriptorDic)
        return dic

    def _get_pdb_ids(self, inputJsonObj, logic_and=0):
        """
            _get_pdb_ids: creates rcsb query parameters & run the query to fetch a lost of rcsb_ids
        """
        logging.info(f'Fetching the list of rcsd_ids with query filter {inputJsonObj}')
        try:
            retIdList = []
            id_score_dict = {}
            i = 0
            for searchType, valList in inputJsonObj.items():
                params = {}
                if searchType in ('InChI', 'InChIKey', 'SMILES'):
                    params = self._create_chem_params(searchType, valList)
                elif searchType in ('sequence', 'uniprot_id', 'ec_number'):
                    params = self._create_seq_ec_uniprot_params(searchType, valList)
                if not params:
                    continue

                new_list, new_scoreDict = self._run_rcsb_search(jsonQueryObj=params)
                if not new_list:
                    continue

                if i == 0:
                    retIdList = new_list
                    id_score_dict = new_scoreDict
                    i = 1
                    continue

                if logic_and:
                    # print('AND logic')
                    # intersecting two lists and two dictionaries using & operator
                    retIdList = list(set(retIdList) & set(new_list))
                    id_score_dict = dict(id_score_dict.items() & new_scoreDict.items())
                else:
                    # print('OR logic')
                    retIdList.extend(new_list)
                    id_score_dict.update(new_scoreDict)

            retIdList = sorted(list(set(retIdList)))

            if not retIdList:
                return {}

            outputObj = {}
            outputObj['total_count'] = len(retIdList)
            outputObj['id_list'] = retIdList
            outputObj['id_score_dict'] = id_score_dict
            outputObj['inputJsonObj'] = inputJsonObj
            logging.info(f'Fetched {len(retIdList)} rcsd_ids with query filter {inputJsonObj}')
            return outputObj
        except Exception as e:
            raise e

    def _filter_by_identity(self, seq_idens, evals, iden_cutoff):
        """
            _filter_by_identity: Filter the max sequence identity with the iden_cutoff
        """
        seq_identity = 0.0
        e_val = 0.0
        exact_match = 0

        if seq_idens:
            seq_idens.sort()
            max_iden = seq_idens.pop()
            if max_iden >= iden_cutoff:  # get the good match
                seq_identity = max_iden
                e_val = evals['Eval_' + str(max_iden)]
                exact_match = 1 if max_iden > 0.99 else 0

        return seq_identity, exact_match, e_val

    def _get_inchiK_cpd(self):
        json_file_path = os.path.join(os.path.dirname(__file__), 'inchikey_cpd.json')

        with open(json_file_path) as DATA:
            inchK_cpd_jsonObj = json.load(DATA)
        return inchK_cpd_jsonObj

    def _blastRCSBSequence(self, input_seq, gql_data, evalue_cutoff, identity_cutoff):
        """
            _blastRCSBSequence: Call self.pdb_util._compute_sequence_identity(seq1, seq2, Eval)
                                 and self._filter_by_identity(seq_idens, iden_cutoff) to fetch the
                                 best sequence identity (similarity) matches
        """
        logging.info(f'Blast {input_seq} against rcsb matches...')

        if not gql_data:
            return []

        rcsb_hits = []
        entries = gql_data['data']['entries']
        rcsb_data = self._formatRCSBJson(entries)
        if not self.inchK_cpd_jsonObj:
            self.inchK_cpd_jsonObj = self._get_inchiK_cpd()

        for entry in entries:
            rcsb_id = entry['rcsb_id']
            rcsb_data_entry = rcsb_data[rcsb_id]
            method = rcsb_data_entry.get('method', [])
            components = []
            if rcsb_data_entry.get('nonpolymer_entities', []):
                inchiKs = rcsb_data_entry['nonpolymer_entities'][0].get('InChIKey', [])
                for inK in inchiKs:
                    components.append(self.inchK_cpd_jsonObj.get(inK, {'InChIKey': inK}))
            prim_cite = rcsb_data_entry.get('primary_citation', {})
            if prim_cite:
                references = [prim_cite.get('pdbx_database_id_PubMed', ''),
                              prim_cite.get('title', ''),
                              prim_cite.get('journal_abbrev', ''),
                              prim_cite.get('rcsb_authors', [])[0],  # Show first author only
                              str(prim_cite.get('year', 'TBD'))]
            else:
                references = []

            for pe in rcsb_data_entry['polymer_entities']:
                pe_seq = pe['one_letter_code_sequence']
                seq_idens, exact_mats, evals = self.pdb_util._compute_sequence_identity(
                                                input_seq, pe_seq, evalue_threshold=evalue_cutoff)
                seq_iden, e_mat, e_val = self._filter_by_identity(seq_idens, evals, identity_cutoff)
                # assemble the data per downstream app requirements
                if seq_iden >= identity_cutoff:
                    rcsb_hits.append({
                        'rcsbid': '_'.join([rcsb_id, pe['id']]),
                        'name': pe.get('uniprot_name', []),
                        'pdbx_sequence': pe_seq,
                        'rcsbec': pe.get('ec_numbers', []),
                        'uniprotec': pe.get('uniprot_ec', []),
                        'identity': seq_iden,
                        'uniprotID': pe.get('uniprotID', []),
                        'pdbx_strand_id': pe.get('pdbx_strand_id', []),
                        'taxonomy': pe.get('taxonomy', []),
                        'evalue': e_val,
                        'method': method,
                        'components': components,
                        'references': references
                    })
        return rcsb_hits

    def _get_graphql_data(self, id_list=[]):
        """
            _get_graphql_data - Query the RCSB GraphQL API, fetch data from GraphQL return
        """
        logging.info(f'Fetch GraphQL data for the structures of {len(id_list)} rcsd_ids')
        # Split large id list into multiple 100 id lists to avoid server connection time-out problem
        chunkSize = 100
        idListChunks = [id_list[i:i+chunkSize] for i in range(0, len(id_list), chunkSize)]

        output_obj = {}
        output_obj['total_count'] = len(id_list)
        output_obj['id_list'] = []

        for idChunk in idListChunks:
            retData = self._queryGraphql(id_list=idChunk)
            if retData.get('data', None) and retData['data'].get('entries', None):
                for key, val in self._formatRCSBJson(retData['data']['entries']).items():
                    output_obj[key] = val
                output_obj['id_list'].extend(idChunk)
            else:
                logging.info("_get_graphql_data by chunks ERROR ".center(30, "-"))

        return output_obj

    def _get_graphql_data_with_cutoffs(self, id_list, input_seq, evalue_cutoff, identity_cutoff):
        """
            _get_graphql_data_with_cutoffs - Query the RCSB GraphQL API, fetch data from GraphQL,
                                             then compute using specified Evalue and filter with
                                             sequence similarity threshold; return formatted data
        """
        logging.info(f'Fetch and filter GraphQL data for the structures of {len(id_list)} rcsd_ids')
        # Split large id list into multiple 100 id lists to avoid server connection time-out problem
        chunkSize = 100
        idListChunks = [id_list[i:i+chunkSize] for i in range(0, len(id_list), chunkSize)]

        output_obj = {}
        output_obj['id_list'] = []
        output_obj[input_seq] = []
        for idChunk in idListChunks:
            retData = self._queryGraphql(id_list=idChunk)
            if retData.get('data', None):
                blast_hits = self._blastRCSBSequence(input_seq, retData,
                                                     evalue_cutoff, identity_cutoff)
                output_obj['id_list'].extend(idChunk)
                output_obj[input_seq].extend(blast_hits)
            else:
                logging.info("_get_graphql_data_with_cutoffs by chunks ERROR ".center(30, "-"))

        return output_obj

    def _rcsb_file_download(self, rcsb_id, ext='pdb', zp='gz'):
        """
            _rcsb_file_download: Download the rcsb structure file and save it to the scratch area
                                Return the saved file path
        """
        rcsb_filename = os.path.join(rcsb_id, ext, zp)
        dir_name = os.path.dirname(__file__)
        rcsb_file = os.path.join(dir_name, rcsb_filename)

        if not self.download_dir:
            self.download_dir = os.path.join(self.scratch, str(uuid.uuid4()))
            os.mkdir(self.download_dir)
        rcsb_filepath = os.path.join(self.download_dir, rcsb_filename)

        try:
            resp = requests.get(os.path.join(self.__baseDownloadUrl, rcsb_filename))
            resp.raise_for_status()
            with open(rcsb_file, 'wb') as rcsb_pt:
                rcsb_pt.write(resp.content)
            shutil.copy(rcsb_file, rcsb_filepath)
            logging.info(f'The rcsb file {rcsb_filename} has been downloaded to {rcsb_filepath}.')
            return rcsb_filepath
        except (HTTPError, ConnectionError, RequestException) as e:
            logging.info(" ERROR ".center(30, "-"))
            logging.info(f'RCSB file downloading for {rcsb_id} had an Error: {e}')
            return ''
        except Exception as e:
            print(e)
            return ''

    def _validate_import_rcsb_params(self, params):
        """
            _validate_import_rcsb_params:
                1) validates input params to import_rcsb_structures and remove unsupported rcsb_id's
                2) For the structure ids given by params['rcsb_ids'], return the
                params with rcsb_infos having the attribute 'file_path' defined.
                Now the params has the following data structure:
                    {
                      'rcsb_infos': [{
                          'file_path': file_path,
                          'file_extension': extension,
                          'narrative_id': narrative_id,
                          'genome_name': genome_name,
                          'feature_id': feature_id,
                          'is_model': is_model
                      }, ...],
                      'structures_name': structures_name,
                      'workspace_name': workspace_name
                    }
            Note: If the file download failed, the corresponding structure entry is removed
                  from params and will be skipped from importing.
        """
        # check for required parameters
        for p in ['workspace_name', 'rcsb_infos', 'structures_name']:
            if p not in params:
                raise ValueError(f'Parameter "{p}" is required, but missing!')

        # Only extensions ‘.cif’ or ‘.pdb’ are valid
        accepted_extensions = ['.pdb', '.cif']
        rcsb_infos = params.get('rcsb_infos', None)
        params['skipped_rcsb_ids'] = list()

        rinfos_deepcopy = deepcopy(rcsb_infos)
        for rinfo in rinfos_deepcopy:
            ext = rinfo.get('extension', '')
            if ext:
                if ext.find('.') == -1:
                    ext = '.' + ext

            if ext not in accepted_extensions:
                logging.info(f"File extension {ext} is not supported at this time, "
                             f"therefore structure {rinfo['rcsb_id']} will not be imported.")
                rcsb_infos.remove(rinfo)
                params['skipped_rcsb_ids'].append(rinfo['rcsb_id'])
                continue

        rinfos_deepcopy = deepcopy(rcsb_infos)
        for rinfo in rinfos_deepcopy:
            file_path = self._rcsb_file_download(rinfo['rcsb_id'], rinfo['extension'])
            if file_path:
                rinfo['file_path'] = file_path
            else:
                logging.info(f"File download for structure {rinfo['rcsb_id']} failed, "
                             f"therefore structure {rinfo['rcsb_id']} will not be imported.")
                rcsb_infos.remove(rinfo)
                params['skipped_rcsb_ids'].append(rinfo['rcsb_id'])

        params['rcsb_infos'] = rcsb_infos
        return params

    def _write_struct_info(self, struct_info):
        """
            _write_struct_info: write the query result info to replace the string
                                '<!--replace query result tbody-->' in the jQuery
                                DataTable jQuery in the template file 'pdbids_report_template.html'
        """
        tbody_html = ''
        id_list = struct_info.get('id_list', [])
        for rcsb_id in id_list:
            if not struct_info.get(rcsb_id, None):
                continue
            pdb_struct = struct_info[rcsb_id]
            tbody_html += '<tr>'

            expl_method = '<br>'.join(pdb_struct.get('method', []))
            """
            prim_cite = pdb_struct.get('primary_citation', {})
            if prim_cite:
                prim_cite_str = '<br>'.join([prim_cite.get('title', ''),
                                             ','.join(prim_cite.get('rcsb_authors', [])),
                                             prim_cite.get('journal_abbrev', ''),
                                             str(prim_cite.get('year', 'TBD'))])
            else:
                prim_cite_str = ''
            """
            prot_sequences = []
            pdb_chains = []
            src_organisms = []
            src_dbs = []
            ec_numbers = []
            for poly_entity in pdb_struct['polymer_entities']:
                if poly_entity.get('one_letter_code_sequence', ''):
                    prot_sequences.append(poly_entity['one_letter_code_sequence'])
                if poly_entity.get('pdbx_strand_id', None):
                    pdb_chains.append('[' + ','.join(poly_entity['pdbx_strand_id']) + ']')
                if poly_entity.get('source_organism', None):
                    for srcg in poly_entity['source_organism']:
                        sci_name = srcg.get('ncbi_scientific_name', '')
                        ncbi_num = f"({srcg.get('ncbi_taxonomy_id', '')})"
                        if sci_name and ncbi_num:
                            src_organisms.append(''.join([sci_name, ncbi_num]))
                if poly_entity.get('ref_sequence_ids', None):
                    for poly_en in poly_entity['ref_sequence_ids']:
                        db_acc = f"({poly_en.get('database_accession', '')})"
                        db_nm = poly_en.get('database_name', '')
                        if db_acc and db_nm:
                            src_dbs.append(''.join([db_nm, db_acc]))
                if poly_entity.get('ec_numbers', None):
                    ec_numbers.append('[' + ','.join(poly_entity['ec_numbers']) + ']')

            inchis = []
            inchikeys = []
            smiles = []
            for nonpoly in pdb_struct['nonpolymer_entities']:
                if nonpoly.get('InChI', None):
                    inchis.append('[' + ','.join(nonpoly['InChI']) + ']')
                if nonpoly.get('InChIKey', None):
                    inchikeys.append('[' + ','.join(nonpoly['InChIKey']) + ']')
                if nonpoly.get('SMILES', None):
                    smiles.append('[' + ','.join(nonpoly['SMILES']) + ']')

            tbody_html += (f'\n<td><a href="https://www.rcsb.org/3d-view/{rcsb_id}"'
                           f'style="cursor: pointer;" target="_blank" '
                           f'title="3D Structure Viewer">{rcsb_id}</a></td>')
            tbody_html += f'\n<td>{expl_method} </td>'
            tbody_html += f'\n<td>{"<br>".join(src_organisms)}</td>'
            tbody_html += f'\n<td>{"<br>".join(ec_numbers)} </td>'
            tbody_html += f'\n<td>{"<br>".join(pdb_chains)} </td>'
            tbody_html += f'\n<td>{"<br>".join(src_dbs)}</td>'
            tbody_html += f'\n<td>{"<br>".join(inchis)} </td>'
            tbody_html += f'\n<td>{"<br>".join(inchikeys)} </td>'
            tbody_html += f'\n<td>{"<br>".join(smiles)}</td>'
            #tbody_html += f'\n<td>{"<br>".join(prot_sequences)}</td>'
            #tbody_html += f'\n<td>{prim_cite_str} </td>'
            tbody_html += '\n</tr>'

        return tbody_html

    def _write_queryNresults(self, rcsb_out, output_directory):
        """
            _write_queryNresults: write to the scratch directory the query and its result
                                 (a list of PDB ids and scores) in JSON files
        """
        input_json = rcsb_out.get('inputJsonObj', {})
        id_score_dict = rcsb_out.get('id_score_dict', {})
        dir_name = os.path.dirname(__file__)

        query_json = os.path.join(dir_name, 'query.json')
        with open(query_json, 'w') as query_pt:
            query_pt.write(json.dumps(input_json, indent=4))

        query_path = os.path.join(output_directory, 'rcsb_query.json')
        shutil.copy(query_json, query_path)
        logging.info(f'The JSON query sent to rcsb has been written to {query_path}.')

        idscores_json = os.path.join(dir_name, 'idscores.json')
        with open(idscores_json, 'w') as rcsbids_pt:
            rcsbids_pt.write(json.dumps(id_score_dict, indent=4))

        idscores_path = os.path.join(output_directory, 'rcsb_ids_scores.json')
        shutil.copy(idscores_json, idscores_path)
        logging.info(f'The rcsb_ids and scores JSON files have been written to {idscores_path}.')
        return (query_path, idscores_path)

    def _query_rcsb_report_html(self, struct_info, rcsb_out, output_directory):
        """
            _query_rcsb_report_html: generates the HTML for the query result-a list of PDB ids
        """
        input_json = rcsb_out.get('inputJsonObj', {})
        html_report = list()

        pdblist_html = self._write_struct_info(struct_info)

        dir_name = os.path.dirname(__file__)
        report_template_file = os.path.join(dir_name, 'templates', 'rcsb_query_report_template.html')
        report_html = os.path.join(dir_name, 'query_structures.html')
        qry_details = json.dumps(input_json, indent=4)
        qry_details += f'\n\nWith Evalue_cutoff={str(self.EVALUE_CUTOFF)} and '
        qry_details += f'\nIdentity_cutoff={str(self.IDENTITY_CUTOFF)}'

        with open(report_html, 'w') as report_html_pt:
            with open(report_template_file, 'r') as report_template_pt:
                # Fetch & fill in detailed info into template HTML
                struct_list_html_report = report_template_pt.read()\
                    .replace('<!--replace query result body-->', pdblist_html)\
                    .replace('<!--replace query details-->', qry_details)
                report_html_pt.write(struct_list_html_report)

        structs_html_report_path = os.path.join(output_directory, 'rcsb_query_report.html')
        shutil.copy(report_html, structs_html_report_path)
        logging.info(f'Full rcsb query report has been written to {structs_html_report_path}')

        html_report.append({'path': output_directory,
                            'name': os.path.basename(structs_html_report_path),
                            'description': 'HTML report for RCSB query results'})

        return html_report

    def _generate_query_report(self, workspace_name, struct_info, rcsb_out,
                               output_directory, rcsb_query, rcsb_ids_scores):
        """
            _generate_query_report: generate summary report for the query
        """
        struct_count = struct_info.get('total_count', 0)
        if not struct_count:
            return {'report_name': None, 'report_ref': None}

        output_html_file = self._query_rcsb_report_html(struct_info, rcsb_out, output_directory)
        output_html_file.append({'path': output_directory,
                                 'name': os.path.basename(rcsb_query),
                                 'description': 'Query applied for retrieving RCSB results'})
        output_html_file.append({'path': output_directory,
                                 'name': os.path.basename(rcsb_ids_scores),
                                 'description': 'RCSB results with scores'})

        report_params = {'message': f'Query has resulted in {struct_count} structures in RCSB DB.',
                         'html_links': output_html_file,
                         'direct_html_link_index': 0,
                         'objects_created': [],
                         'workspace_name': workspace_name,
                         'report_object_name': 'query_rcsb_report_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)
        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def query_rcsb_with_file(self, inputJsonFilePath, logic='or'):
        try:
            if not os.access(inputJsonFilePath, os.R_OK):
                return

            with open(inputJsonFilePath) as DATA:
                inputJsonObj = json.load(DATA)

            rcsb_output = self._get_pdb_ids(inputJsonObj)
            logging.info(f'{rcsb_output["total_count"]} RCSB Structures found.')
            return self._get_graphql_data(rcsb_output.get('id_list', []))
        except Exception as e:
            raise e

    def query_structure_info(self, params, create_report=1):
        """
            query_structure_info: with given constraints, query structure info from RCSB database
        """
        # logging.info(f'query_structure_info with params: {params}')
        # fetch the query filters and assemble them into a json object
        inputJsonObj, workspace_name = self._validate_rcsb_query_params(params)

        idlist = []
        struct_info = {}
        returnVal = {}
        returnVal['rcsb_ids'] = []
        returnVal['rcsb_scores'] = {}
        returnVal['report_ref'] = None
        returnVal['report_name'] = None

        # Make output/report directory and copy over uploaded pdb files
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(output_directory)

        rcsb_output = self._get_pdb_ids(inputJsonObj, self.LOGICAL_AND)
        total_count = rcsb_output.get('total_count', 0)
        idlist = rcsb_output.get('id_list', [])
        id_score_dict = rcsb_output.get('id_score_dict', {})
        if total_count and idlist:
            logging.info(f'{total_count} RCSB Structures found.')
            struct_info = self._get_graphql_data(idlist)
            logging.info('Retrieved structure information:')
            logging.info(f"total_count={struct_info.get('total_count', 0)}")
            returnVal['rcsb_ids'] = struct_info.get('id_list', [])
            returnVal['rcsb_scores'] = id_score_dict
            if create_report:
                (query_json, ids_scores) = self._write_queryNresults(rcsb_output, output_directory)
                report_output = self._generate_query_report(workspace_name, struct_info,
                                                            rcsb_output, output_directory,
                                                            query_json, ids_scores)
                returnVal.update(report_output)
            else:
                returnVal['report_ref'] = None
                returnVal['report_name'] = None
        return returnVal

    def query_structure_anno(self, params):
        """
            query_structure_anno: with given constraints, query structure info from RCSB database
        """
        # logging.info(f'query_structure_anno with params: {params}')
        params = self._validate_rcsb_seqquery_params(params)
        returnVal = {}

        # search by sequence one-by-one
        for input_seq in params['sequence_strings']:
            rcsb_output = self._get_pdb_ids({'sequence': [input_seq]}, self.LOGICAL_AND)
            idlist = rcsb_output.get('id_list', [])
            returnVal[input_seq] = []
            if idlist:
                search_ret = self._get_graphql_data_with_cutoffs(idlist, input_seq,
                                                                 self.EVALUE_CUTOFF,
                                                                 self.IDENTITY_CUTOFF)
                logging.info(f'Retrieved structure information for {input_seq}')
                returnVal[input_seq] = search_ret[input_seq]

        return returnVal

    def import_rcsbs(self, params, workspace_name):
        """
            import_rcsbs: uploading the rcsb structures
        """
        pdb_objects = list()
        pdb_infos = list()
        successful_ids = list()
        skipped_ids = params.get('skipped_rcsb_ids', list())

        rcsb_infos = params['rcsb_infos']
        # loop through the list of pdb_file_paths
        for rinfo in rcsb_infos:
            pdb_params = {}
            file_path = rinfo['file_path']
            rid = rinfo['rcsb_id']

            pdb_params['pdb_info'] = rinfo
            pdb_params['input_staging_file_path'] = None
            pdb_params['input_file_path'] = file_path
            pdb_params['input_shock_id'] = None
            pdb_params['workspace_name'] = workspace_name
            pdb_params['structure_name'] = rid
            pdb_params['is_model'] = rinfo['is_model']

            if 'pdb' in rinfo['file_extension']:
                pdb_data, pdb_info = self.pdb_util.import_pdb_file(pdb_params)
                if pdb_data:
                    pdb_objects.append(pdb_data)
                    pdb_infos.append(pdb_info)
                    successful_ids.append(file_path)
                else:
                    skipped_ids.append(rid)
            elif 'cif' in rinfo['file_extension']:
                cif_data, pdb_info = self.pdb_util.import_mmcif_file(pdb_params)
                if cif_data:
                    pdb_objects.append(cif_data)
                    pdb_infos.append(pdb_info)
                    successful_ids.append(file_path)
                else:
                    skipped_ids.append(rid)

        return pdb_objects, pdb_infos, successful_ids, skipped_ids

    def batch_import_rcsbs(self, params):
        """
            batch_import_rcsbs: download a list of rcsb files and then upload them to create a
                                   KBaseStructure.ProteinStructures object
            required params:
                rcsb_ids: a list of rcsb_id's
                structures_name: name of the ProteinStructures object to be generated
                workspace_name: workspace name that the protein structure(s) will be saved
            return:
                structures_ref: return ProteinStructures object reference
                report_name: name of generated report (if any)
                report_ref: report reference (if any)

            1. call _validate_import_rcsb_params to validate input params and download
               the rcsb structure files one by one
            2. upload: call import_rcsbs()
            3. assemble the data for a ProteinStructures and save the data object
            4. call PDBUtil.saveStructures_createReport to generate a report for batch_import_pdbs'.
        """
        params = self._validate_import_rcsb_params(params)

        workspace_name = params.get('workspace_name', '')
        structures_name = params.get('structures_name', '')
        if not isinstance(workspace_name, int):
            workspace_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            workspace_id = workspace_name
        params['workspace_id'] = workspace_id

        pdb_objects = list()
        pdb_infos = list()
        successful_ids = list()
        protein_structures = dict()

        pdb_objects, pdb_infos,
        successful_ids, skipped_ids = self.import_rcsbs(params, workspace_name)

        if not pdb_objects:
            logging.info("No pdb structure was created/saved!")
            return {}

        total_structures = len(pdb_objects)
        protein_structures['protein_structures'] = pdb_objects
        protein_structures['total_structures'] = total_structures
        protein_structures['description'] = (f'Created {total_structures} '
                                             f'structures in {params.get("structures_name")}')

        return self.pdb_util.saveStructures_createReport(structures_name, workspace_id,
                                                         workspace_name, protein_structures,
                                                         pdb_infos, skipped_ids)
