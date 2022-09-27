# -*- coding: utf-8 -*-
import os
import shutil
import time
import json
import unittest
from configparser import ConfigParser
#from mock import patch
#from unittest.mock import Mock, patch

from ProteinStructureUtils.ProteinStructureUtilsImpl import ProteinStructureUtils
from ProteinStructureUtils.ProteinStructureUtilsServer import MethodContext
from ProteinStructureUtils.authclient import KBaseAuth as _KBaseAuth

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace

from ProteinStructureUtils.Utils.RCSBUtil import RCSBUtil


class RCSBUtilsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('ProteinStructureUtils'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'ProteinStructureUtils',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = ProteinStructureUtils(cls.cfg)
        cls.cfg['USER_ID'] = cls.ctx['user_id']
        cls.rcsb_util = RCSBUtil(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        cls.dfu = DataFileUtil(cls.callback_url)
        suffix = int(time.time() * 1000)
        cls.wsName = "test_RCSBUtils_" + str(suffix)
        cls.ws_id = cls.wsClient.create_workspace({'workspace': cls.wsName})[0]
        cls.prepareInputJsonData()
        cls.EVALUE_CUTOFF = 0.1
        cls.IDENTITY_CUTOFF = 0.75

    @classmethod
    def prepareInputJsonData(cls):
        file1 = 'input_sequence_uniprot_ec.json'
        file2 = 'input_source_organism.json'
        file3 = 'input_chemical_descriptor.json'

        cls.inputjson_file_path1 = os.path.join(cls.scratch, file1)
        cls.inputjson_file_path2 = os.path.join(cls.scratch, file2)
        cls.inputjson_file_path3 = os.path.join(cls.scratch, file3)
        shutil.copy(os.path.join('data/rcsb_data', file1), cls.inputjson_file_path1)
        shutil.copy(os.path.join('data/rcsb_data', file2), cls.inputjson_file_path2)
        shutil.copy(os.path.join('data/rcsb_data', file3), cls.inputjson_file_path3)

        with open(cls.inputjson_file_path1) as DATA1:
            cls.inputJsonObj1 = json.load(DATA1)
        with open(cls.inputjson_file_path2) as DATA2:
            cls.inputJsonObj2 = json.load(DATA2)
        with open(cls.inputjson_file_path3) as DATA3:
            cls.inputJsonObj3 = json.load(DATA3)

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    @unittest.skip('test_check_input_jsons')
    def test_check_input_jsons(self):
        keys1 = ['sequence', 'ec_number', 'uniprot_id']
        for k1 in keys1:
            self.assertIn(k1, self.inputJsonObj1)

        keys2 = ['Sorghum', 'Miscanthus', 'Camelina', 'Pennycress', 'Poplar', 'Switchgrass',
                 'Caulobacter', 'Clostridial species', 'Cyanobacteria', 'Yeasts', 'E. coli',
                 'Pseudomonas species']
        for k2 in keys2:
            self.assertIn(k2, self.inputJsonObj2)

        keys3 = ['InChI', 'SMILES']
        for k3 in keys3:
            self.assertIn(k3, self.inputJsonObj3)

    # Testing RCSBUtil module functions
    @unittest.skip('test_create_seq_ec_uniprot_params_seq')
    def test_create_seq_ec_uniprot_params_seq(self):
        # test single sequence
        search_type11 = 'sequence'
        valList11 = ["MNLTELKNTPVSELITLGENMGLENLARMRKQDIIF"]
        pm11 = self.rcsb_util._create_seq_ec_uniprot_params(search_type11, valList11)
        expected_pm11 = {
            "query": {
                "type": "terminal",
                "service": "sequence",
                "parameters": {
                    "evalue_cutoff": self.EVALUE_CUTOFF,
                    "identity_cutoff": self.IDENTITY_CUTOFF,
                    "target": "pdb_protein_sequence",
                    "value": valList11[0]
                }
            },
            "return_type": "entry",
            "request_options": {
                "return_all_hits": True
            }
        }
        self.assertEqual(pm11['query']['service'], 'sequence')
        self.assertEqual(pm11['query']['parameters']['value'], valList11[0])
        self.assertEqual(pm11['query']['parameters']['target'], 'pdb_protein_sequence')
        self.assertEqual(pm11['request_options'], expected_pm11['request_options'])
        self.assertDictEqual(pm11, expected_pm11)

        # test mutiple sequences, with search_type as sequence_strings
        search_type12 = 'sequence_strings'
        valList12 = ["MNLTELKNTPVSELITLGENMGLENLARMRKQDIIF",
                     "ATPHINAEMGDFADVVLMPGDPLRAKYIAETFLEDAREVNICTVSDHIRTHEQTTAAER"]
        pm12 = self.rcsb_util._create_seq_ec_uniprot_params(search_type12, valList12)
        self.assertEqual(pm12['query']['type'], 'group')
        self.assertEqual(pm12['query']['logical_operator'], 'or')
        self.assertEqual(pm12['request_options'], expected_pm11['request_options'])
        ret_nodelist = pm12['query']['nodes']
        self.assertEqual(len(ret_nodelist), 2)
        self.assertEqual(ret_nodelist[0]['parameters']['value'], valList12[0])
        self.assertEqual(ret_nodelist[1]['parameters']['value'], valList12[1])

    @unittest.skip('test_create_seq_ec_uniprot_params_uniprot')
    def test_create_seq_ec_uniprot_params_uniprot(self):
        search_type2 = "uniprot_ids"
        valList2 = ["Q01532", "R9RYW2"]
        pm2 = self.rcsb_util._create_seq_ec_uniprot_params(search_type2, valList2)

        expected_pm2 = {
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
                            "value": valList2
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
        self.assertDictEqual(pm2, expected_pm2)

        search_type2b = "uniprot_wrong_key"
        pm_empty = self.rcsb_util._create_seq_ec_uniprot_params(search_type2b, valList2)
        self.assertFalse(pm_empty)

    @unittest.skip('test_create_seq_ec_uniprot_params_ec')
    def test_create_seq_ec_uniprot_params_ec(self):
        search_type3 = "ec_numbers"
        valList3 = ["3.5.2.6", "3.4.13.9"]
        pm3 = self.rcsb_util._create_seq_ec_uniprot_params(search_type3, valList3)

        expected_pm3 = {
            "query": {
                "type": "terminal",
                "label": "text",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity.rcsb_ec_lineage.id",
                    "operator": "in",
                    "negation": False,
                    "value": valList3
                }
            },
            "return_type": "entry",
            "request_options": {
                "return_all_hits": True
            }
        }
        self.assertDictEqual(pm3, expected_pm3)

        search_type3b = "ec_num_wrong_key"
        pm_empty = self.rcsb_util._create_seq_ec_uniprot_params(search_type3b, valList3)
        self.assertFalse(pm_empty)

    @unittest.skip('test_get_pdb_ids_by_sequence_uniprot_ec')
    def test_get_pdb_ids_by_sequence_uniprot_ec(self):
        ret = self.rcsb_util._get_pdb_ids(self.inputJsonObj1)
        k1 = 'total_count'
        k2 = 'id_list'
        self.assertIn(k1, ret, f'Key {k1} not in returned object')
        self.assertIn(k2, ret, f'Key {k2} not in returned object')
        self.assertEqual(len(ret[k2]), ret[k1])
        print(f'RCSB query by sequence/ecnum/uniprotid returned {ret[k1]} ids.')

    @unittest.skip('test_get_pdb_ids_by_source_organism')
    def test_get_pdb_ids_by_source_organism(self):
        ret = self.rcsb_util._get_pdb_ids(self.inputJsonObj2)
        k1 = 'total_count'
        k2 = 'id_list'
        self.assertIn(k1, ret, f'Key {k1} not in returned object')
        self.assertIn(k2, ret, f'Key {k2} not in returned object')
        self.assertEqual(len(ret[k2]), ret[k1])
        print(f'RCSB query by source organism returned {ret[k1]} ids.')

    @unittest.skip('test_get_pdb_ids_by_chem')
    def test_get_pdb_ids_by_chem(self):
        ret = self.rcsb_util._get_pdb_ids(self.inputJsonObj3)
        k1 = 'total_count'
        k2 = 'id_list'
        self.assertIn(k1, ret, f'Key {k1} not in returned object')
        self.assertIn(k2, ret, f'Key {k2} not in returned object')
        self.assertEqual(len(ret[k2]), ret[k1])
        print(f'RCSB query by chem returned {ret[k1]} ids.')

    #@unittest.skip('test_get_graphql_data')
    def test_get_graphql_data(self):
        id_list = ['1A0I', '1A49', '1A5U', '1A82', '1AQ2']
        gql_ret = self.rcsb_util._get_graphql_data(id_list)

        expected_keys = ['total_count', 'id_list']
        expected_keys.extend(id_list)
        #print(gql_ret)
        for k in expected_keys:
            self.assertIn(k, gql_ret)

        self.assertEqual(gql_ret['total_count'], len(id_list))
        self.assertEqual(gql_ret['id_list'], id_list)
