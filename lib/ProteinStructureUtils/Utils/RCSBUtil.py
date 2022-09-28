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

import json, requests
import aiohttp, asyncio
from requests.exceptions import ConnectionError, HTTPError, RequestException
from python_graphql_client import GraphqlClient

from installed_clients.AbstractHandleClient import AbstractHandle
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace


class RCSBUtil:

    # Set thresholds to restrict which alignments will be significant or sequence identity matches
    EVALUE_CUTOFF = 0.1
    IDENTITY_CUTOFF = 0.75

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.token = config['KB_AUTH_TOKEN']
        self.user_id = config['USER_ID']
        self.dfu = DataFileUtil(self.callback_url)
        self.hs = AbstractHandle(config['handle-service-url'])
        self.ws_client = Workspace(config['workspace-url'])
        self.shock_url = config['shock-url']

        self.__baseSearchUrl = "https://search.rcsb.org/rcsbsearch/v2/query"
        self.__baseGraphqlUrl = "https://data.rcsb.org/graphql"
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

    def _validate_rcsb_query_params(self, params):
        """
            _validate_rcsb_query_params:
                validates input params to query_structure_info
        """
        # check for required parameters
        for p in ['workspace_name']:
            if p not in params:
                raise ValueError(f'Parameter "{p}" is required, but missing!')

        # check for queriable parameters
        inputJsonObj = {}
        for p in ['sequence_strings', 'uniprot_ids', 'ec_numbers', 'inchis', 'smiles']:
            if params.get(p, None):
                inputJsonObj[p] = params[p]

        if not inputJsonObj:
            raise ValueError('No queriable parameters found!')

        return inputJsonObj, params.get('workspace_name')

    def _create_seq_ec_uniprot_params(self, searchType, valList):
        """
            _create_seq_ec_uniprot_params - Based on sequence_strings/uniprot_ids/ec_numbers input,
                                            create query parameters in Json object format.
                                            return the seqEcUniprotJsonObject
        """
        seqEcUniprotJsonObject = {}
        if searchType not in ('sequence', 'sequence_strings', 'uniprot_id',
                              'uniprot_ids','ec_number', 'ec_numbers'):
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
            logging.info(" ERROR ".center(30, "-"))
            logging.info(f'Querying RCSB db with {jsonQueryObj} had an Error: {e}')
            return {}
        except Exception as e:
            print(e)
            return {}

    def _readRCSBResult(self, rcsb_retObj):
        """parse the RCSB query result and return a PDB ID list
        """
        total_number = 0
        retIdList = []
        if rcsb_retObj and rcsb_retObj.get('total_count', None):
            total_number = rcsb_retObj["total_count"]
        else:
            return []

        # Read entry ID from "result_set" list
        if "result_set" in rcsb_retObj:
            for obj in rcsb_retObj["result_set"]:
                if "identifier" in obj:
                    retIdList.append(obj["identifier"])
            if len(retIdList) > 1:
                retIdList.sort()

        if len(retIdList) != total_number:
            errMsg = (f"The total_count={total_number}, "
                      "but only got {len(retIdList)} id(s) from result_set.")
            raise ValueError(errMsg)

        return retIdList

    def _run_rcsb_search(self, jsonQueryObj={}):
        """ Run search API to get return PDB ID list like:
                [
                    "117E",
                    "11AS",
                    "12AS",
                    "148L",
                    "1A04",
                    "1A0A",
                    "1A0B",
                    "1A0F",
                  ......
               ]
        """
        if not jsonQueryObj:
            return []
        #
        retJsonObj = self._queryRCSB(jsonQueryObj)
        if not retJsonObj:
            return []
        return self._readRCSBResult(retJsonObj)

    def _queryGraphql(self, id_list=[]):
        """
            _queryGraphql: Given a list of rcsb_ids, Query the RCSB GraphQL API
        """
        logging.info(f'Querying GraphQL db for the structures of {len(id_list)} rcsd_ids')
        if not id_list:
            return {}

        queryString = self.__graphqlQueryTemplate % '", "'.join(id_list)
        try:
            #graphql_ret = self.__graphqlClient.execute(query=queryString)
            # asyncio.run() is available ONLY for python 3.7 or newer
            # graphql_ret = asyncio.run(self.__graphqlClient.execute_async(query=queryString))
            evt_loop = asyncio.get_event_loop()
            graphql_ret = evt_loop.run_until_complete(
                          self.__graphqlClient.execute_async(query=queryString))
            #logging.info(f'Querying GraphQL db returned the structures {graphql_ret}')
            return graphql_ret
        except (aiohttp.ClientConnectionError, asyncio.TimeoutError):
            logging.info('Connecting to RCSB GraphQL host with URL "{self.__baseGraphqlUrl} errored.')
            return {}
        except (HTTPError, ConnectionError, RequestException) as e:
            # not raising error to allow continue with other chunks
            logging.info(f'Querying RCSB GraphQL had a Connection Error:************\n {e}.\n'
                         'Or database connection request had no response!')
            return {}
        except (RuntimeError, TypeError, KeyError, ValueError) as e:
            err_msg = f'Querying RCSB errored with message: {e.message} and data: {e.data}'
            raise ValueError(err_msg)

    def _formatRCSBJson(self, data):
        """
            _formatRCSBJson: Format rcsb GraphQL returned data into Json objects per rcsb structure
        """
        if not data:
            return {}

        dic = {}
        # short-naming the long rcsb data attribute strings
        src_organism = 'rcsb_entity_source_organism'
        prim_cite = 'rcsb_primary_citation'
        polym_entity = 'rcsb_polymer_entity'
        ec_lineage = 'rcsb_ec_lineage'
        entity_container_ids = 'rcsb_polymer_entity_container_identifiers'
        ref_seq_ids = 'reference_sequence_identifiers'
        chem_comp_desc = 'rcsb_chem_comp_descriptor'

        for entry in data['data']['entries']:
            dic[entry['rcsb_id']] = {
                'method': [entry['exptl'][0]['method']],
                'primary_citation': entry[prim_cite],
                'polymer_entities': []
            }
            for pe in entry['polymer_entities']:
                pe_dic = {
                    'id': pe['rcsb_id'].split("_")[1],
                    'one_letter_code_sequence': pe['entity_poly']['pdbx_seq_one_letter_code'],
                    'pdb_chain_ids': pe['entity_poly']['pdbx_strand_id'].split(',')
                }
                if (pe.get(src_organism, None)):
                    pe_dic['source_organism'] = pe[src_organism]
                #
                if (pe.get(entity_container_ids, None) and
                        pe[entity_container_ids].get(ref_seq_ids, None)):
                    pe_dic['identifiers'] = pe[entity_container_ids]
                #
                if (pe.get(polym_entity, None) and pe[polym_entity].get(ec_lineage, None)):
                    ec_numbers = []
                    for idic in pe[polym_entity][ec_lineage]:
                        if idic['id'] and (idic['id'] not in ec_numbers):
                            ec_numbers.append(idic['id'])

                    if ec_numbers:
                        pe_dic['ec_numbers'] = ec_numbers
                #
                dic[entry['rcsb_id']]['polymer_entities'].append(pe_dic)
            #
            dic[entry['rcsb_id']]['nonpolymer_entities'] = []
            if entry.get('nonpolymer_entities', None):
                for nonpolymerDic in entry['nonpolymer_entities']:
                    descriptorDic = {}
                    if (not nonpolymerDic.get('nonpolymer_comp', None) or not
                       nonpolymerDic['nonpolymer_comp'].get('rcsb_chem_comp_descriptor', None)):
                        continue
                    #
                    for desType in ('InChI', 'InChIKey', 'SMILES'):
                        if not nonpolymerDic['nonpolymer_comp'][chem_comp_desc].get(desType, None):
                            continue
                        #
                        dt_val = nonpolymerDic['nonpolymer_comp'][chem_comp_desc][desType]
                        if desType in descriptorDic:
                            descriptorDic[desType].append(dt_val)
                        else:
                            descriptorDic[desType] = [dt_val]
                        #
                    #
                    dic[entry['rcsb_id']]['nonpolymer_entities'].append(descriptorDic)
                #
            #
        #
        return dic

    def _get_pdb_ids(self, inputJsonObj, logic='or'):
        """
            _get_pdb_ids: creates rcsb query parameters & run the query to fetch a lost of rcsb_ids
        """
        try:
            retIdList = []
            i = 0
            for searchType, valList in inputJsonObj.items():
                params = {}
                if searchType in ('InChI', 'InChIKey', 'SMILES', 'inchis', 'smiles'):
                    params = self._create_chem_params(searchType, valList)
                elif searchType in ('sequence', 'sequence_strings', 'uniprot_id', 'ec_number',
                                    'uniprot_ids', 'ec_numbers'):
                    params = self._create_seq_ec_uniprot_params(searchType, valList)
                if not params:
                    continue

                new_list = self._run_rcsb_search(jsonQueryObj=params)
                if i == 0:
                    retIdList = new_list
                    i = 1
                    continue

                if logic == 'and':
                    # print('AND logic')
                    retIdList = list(set(retIdList) & set(new_list))
                else:
                    # print('OR logic')
                    retIdList.extend(new_list)
            #
            retIdList = sorted(list(set(retIdList)))
            #
            outputObj = {}
            outputObj['total_count'] = len(retIdList)
            outputObj['id_list'] = retIdList
            return outputObj
        except Exception as e:
            raise e

    def _get_graphql_data(self, id_list=[]):
        """
            _get_graphql_data - Query the RCSB GraphQL API, fetch data from GraphQL return
        """
        logging.info(f'Fetch GraphQL data for the structures of {len(id_list)} rcsd_ids')
        # Split large id list into multiple 100 id lists to avoid server connection time-out problem
        chunkSize = 100
        idListChunks = [id_list[i:i+chunkSize] for i in range(0, len(id_list), chunkSize)]

        output_obj = {}
        output_obj["total_count"] = len(id_list)
        output_obj["id_list"] = id_list

        for idChunk in idListChunks:
            retData = self._queryGraphql(id_list=idChunk)
            for key, val in self._formatRCSBJson(retData).items():
                output_obj[key] = val

        return output_obj

    def _write_struct_info(self, output_dir, struct_info):
        """
            _write_struct_info: write the query result info to replace the string
                                '<!--replace query result tbody-->' in the jQuery
                                DataTable jQuery in the template file 'pdbids_report_template.html'
        """
        tbody_html = ''
        id_list = struct_info.get('id_list', [])
        for rcsb_id in id_list:
            pdb_struct = struct_info[rcsb_id]
            tbody_html += '<tr>'

            expl_method = ';'.join(pdb_struct.get('method', []))
            prim_cite = pdb_struct.get('primary_citation', {})
            if prim_cite:
                prim_cite_str = ';'.join([prim_cite.get('title', ''),
                                          ','.join(prim_cite.get('rcsb_authors', [])),
                                          prim_cite.get('journal_abbrev', ''),
                                          str(prim_cite.get('year', 'TBD'))])
            else:
                prim_cite_str = ''
            prot_sequences = []
            pdb_chains = []
            src_organisms = []
            src_dbs = []
            ec_numbers = []
            for poly_entity in pdb_struct['polymer_entities']:
                if poly_entity.get('one_letter_code_sequence', None):
                    prot_sequences.extend(poly_entity['one_letter_code_sequence'])
                if poly_entity.get('pdb_chain_ids', None):
                    pdb_chains.extend(poly_entity['pdb_chain_ids'])
                if poly_entity.get('source_organism', None):
                    src_organisms.extend(poly_entity['source_organism'])
                if poly_entity.get('identifiers', None):
                    src_dbs.extend(poly_entity['identifiers'])
                if poly_entity.get('ec_numbers', None):
                    ec_numbers.extend(poly_entity['ec_numbers'])

            inchis = []
            inchikeys = []
            smiles = []
            for nonpoly in pdb_struct['nonpolymer_entities']:
                if nonpoly.get('InChI', None):
                    inchis.extend(nonpoly['InChIKey'])
                if nonpoly.get('InChIKey', None):
                    inchikeys.extend(nonpoly['InChIKey'])
                if nonpoly.get('SMILES', None):
                    smiles.extend(nonpoly['SMILES'])

            tbody_html += (f'\n<td><a href="https://www.rcsb.org/3d-view/{rcsb_id}"'
                           f'style="cursor: pointer;" '
                           f'title="3D Structure Viewer">{rcsb_id}</a></td>')
            tbody_html += f'\n<td>{expl_method} </td>'
            tbody_html += (f'\n<td>{prim_cite_str} </td>')
            tbody_html += f'\n<td>{pdb_chains} </td>'
            tbody_html += f'\n<td>{src_organisms}</td>'
            tbody_html += f'\n<td>{ec_numbers} </td>'
            tbody_html += f'\n<td>{src_dbs}</td>'
            tbody_html += f'\n<td>{inchis} </td>'
            tbody_html += f'\n<td>{inchikeys} </td>'
            tbody_html += f'\n<td>{smiles}</td>'
            tbody_html += f'\n<td>{prot_sequences} </td>'
            tbody_html += '\n</tr>'

        return tbody_html

    def _query_rcsb_report_html(self, struct_info):
        """
            _query_rcsb_report_html: generates the HTML for the query result-a list of PDB ids
        """
        html_report = list()

        # Make report directory and copy over uploaded pdb files
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(output_directory)

        pdblist_html = self._write_struct_info(output_directory, struct_info)

        dir_name = os.path.dirname(__file__)
        report_template_file = os.path.join(dir_name, 'templates', 'rcsb_query_report_template.html')
        report_html = os.path.join(dir_name, 'query_structures.html')

        with open(report_html, 'w') as report_html_pt:
            with open(report_template_file, 'r') as report_template_pt:
                # Fetch & fill in detailed info into template HTML
                struct_list_html_report = report_template_pt.read()\
                    .replace('<!--replace query result body-->', pdblist_html)
                report_html_pt.write(struct_list_html_report)

        structs_html_report_path = os.path.join(output_directory, 'rcsb_query_report.html')
        shutil.copy(report_html, structs_html_report_path)
        logging.info(f'Full rcsb query report has been written to {structs_html_report_path}')

        html_report.append({'path': output_directory,
                            'name': os.path.basename(structs_html_report_path),
                            'description': 'HTML report for RCSB query results'})

        return html_report

    def _generate_query_report(self, workspace_name, struct_info):
        """
            _generate_query_report: generate summary report for the query
        """
        struct_count = struct_info.get('total_count', 0)
        output_html_file = self._query_rcsb_report_html(struct_info)

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

    def querey_structure_info(self, params):
        """
            query_structure_info: with given constraints, query structure info from RCSB database
        """
        logging.info(f'query_structure_info with params: {params}')

        # fetch the query filters and assemble them into a json object
        inputJsonObj, workspace_name = self._validate_rcsb_query_params(params)

        idlist = []
        struct_info = {}
        returnVal = {}
        try:
            rcsb_output = self._get_pdb_ids(inputJsonObj)
            logging.info(f'{rcsb_output["total_count"]} RCSB Structures found.')
            idlist = rcsb_output.get('id_list', [])
            struct_info = self._get_graphql_data(idlist)
        except (RuntimeError, TypeError, KeyError, ValueError) as e:
            logging.info(e)
            return {}
        else:
            print('Retrieved structure information:')
            print(f'total_count={struct_info.get("total_count", 0)}')
            print(f'first structure: {struct_info.get(idlist[0]).get("primary_citation")}')
            returnVal['rcsb_ids'] = idlist
            report_output = self._generate_query_report(workspace_name, struct_info)
            returnVal.update(report_output)
            return returnVal
