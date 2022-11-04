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
        cls.wsName = "test_RCSBUtils_" + str(int(time.time() * 1000))
        cls.ws_id = cls.wsClient.create_workspace({'workspace': cls.wsName})[0]
        cls.prepareInputJsonData()
        cls.EVALUE_CUTOFF = 1e-5
        cls.IDENTITY_CUTOFF = 0.75

    @classmethod
    def prepareInputJsonData(cls):
        file1 = 'input_sequence_uniprot_ec.json'
        file2 = 'input_source_organism.json'
        file3 = 'input_chemical_descriptor.json'
        file4 = 'gql_ret.json'
        file5 = 'gqldata_itms.json'
        file6 = 'input_chemical_pyridine.json'

        cls.inputjson_file_path1 = os.path.join(cls.scratch, file1)
        cls.inputjson_file_path2 = os.path.join(cls.scratch, file2)
        cls.inputjson_file_path3 = os.path.join(cls.scratch, file3)
        cls.gqlret_file_path = os.path.join(cls.scratch, file4)
        cls.gqldata_itm_path = os.path.join(cls.scratch, file5)
        cls.pyridine_file_path = os.path.join(cls.scratch, file6)
        shutil.copy(os.path.join('data/rcsb_data', file1), cls.inputjson_file_path1)
        shutil.copy(os.path.join('data/rcsb_data', file2), cls.inputjson_file_path2)
        shutil.copy(os.path.join('data/rcsb_data', file3), cls.inputjson_file_path3)
        shutil.copy(os.path.join('data/rcsb_data', file4), cls.gqlret_file_path)
        shutil.copy(os.path.join('data/rcsb_data', file5), cls.gqldata_itm_path)
        shutil.copy(os.path.join('data/rcsb_data', file6), cls.pyridine_file_path)

        with open(cls.inputjson_file_path1) as DATA1:
            cls.inputJsonObj1 = json.load(DATA1)
        with open(cls.inputjson_file_path2) as DATA2:
            cls.inputJsonObj2 = json.load(DATA2)
        with open(cls.inputjson_file_path3) as DATA3:
            cls.inputJsonObj3 = json.load(DATA3)
        with open(cls.gqlret_file_path) as GQL_DATA:
            cls.gql_data = json.load(GQL_DATA)
        with open(cls.gqldata_itm_path) as ITM_DATA:
            cls.gqldata_items = json.load(ITM_DATA)
        with open(cls.pyridine_file_path) as PRD_DATA:
            cls.prd_data = json.load(PRD_DATA)

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

        keys4 = ['rcsb_id', 'exptl', 'rcsb_primary_citation',
                 'polymer_entities', 'nonpolymer_entities']
        for k4 in keys4:
            gql_ret_entries = self.gql_data.get('data').get('entries')
            self.assertIn(k4, gql_ret_entries[0])
            self.assertIn(k4, gql_ret_entries[1])

        keys5 = ['total_count', 'id_list', '1A0I', '1A49']
        for k5 in keys5:
            self.assertIn(k5, self.gqldata_items)

        keys6 = ['InChI', 'SMILES', 'InChIKey']
        for k6 in keys6:
            self.assertIn(k6, self.prd_data)

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

    @unittest.skip('test_create_chem_params_inchi')
    def test_create_chem_params_inchi(self):
        # test single InChI
        search_type1 = "InChI"
        valList11 = ["InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)"]
        chem_pm11 = self.rcsb_util._create_chem_params(search_type1, valList11)

        expected_chem_pm1 = {
            "query": {
                "type": "terminal",
                "service": "chemical",
                "parameters": {
                    "value": valList11[0],
                    "type": "descriptor",
                    "descriptor_type": search_type1,
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
        self.assertDictEqual(chem_pm11, expected_chem_pm1)
        self.assertEqual(chem_pm11['query']['type'], 'terminal')
        self.assertEqual(chem_pm11['query']['service'], 'chemical')
        self.assertEqual(chem_pm11['query']['parameters']['type'], 'descriptor')
        self.assertEqual(chem_pm11['query']['parameters']['value'], valList11[0])
        self.assertEqual(chem_pm11['query']['parameters']['descriptor_type'], search_type1)
        self.assertEqual(chem_pm11['query']['parameters']['match_type'], 'graph-exact')

        # test multiple InChIs
        valList12 = ["InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)",
                    "InChI=1S/C10H16N5O13P3/madeup"]
        chem_pm12 = self.rcsb_util._create_chem_params(search_type1, valList12)

        self.assertEqual(chem_pm12['query']['type'], 'group')
        self.assertEqual(chem_pm12['query']['logical_operator'], 'or')
        self.assertEqual(chem_pm12['request_options'], expected_chem_pm1['request_options'])
        ret_nodelist = chem_pm12['query']['nodes']
        self.assertEqual(len(ret_nodelist), 2)
        self.assertEqual(ret_nodelist[0]['type'], 'terminal')
        self.assertEqual(ret_nodelist[0]['service'], 'chemical')
        self.assertEqual(ret_nodelist[0]['parameters']['value'], valList12[0])
        self.assertEqual(ret_nodelist[0]['parameters']['descriptor_type'], search_type1)
        self.assertEqual(ret_nodelist[0]['parameters']['match_type'], 'graph-exact')
        self.assertEqual(ret_nodelist[1]['parameters']['value'], valList12[1])

    @unittest.skip('test_create_chem_params_smiles')
    def test_create_chem_params_smiles(self):
        # test single SMILES
        search_type2 = "SMILES"
        valList21 = ["c1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO[P@@](=O)(O)O)N"]
        chem_pm21 = self.rcsb_util._create_chem_params(search_type2, valList21)

        expected_chem_pm2 = {
            "query": {
                "type": "terminal",
                "service": "chemical",
                "parameters": {
                    "value": valList21[0],
                    "type": "descriptor",
                    "descriptor_type": search_type2,
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
        self.assertDictEqual(chem_pm21, expected_chem_pm2)
        self.assertEqual(chem_pm21['query']['type'], 'terminal')
        self.assertEqual(chem_pm21['query']['service'], 'chemical')
        self.assertEqual(chem_pm21['query']['parameters']['type'], 'descriptor')
        self.assertEqual(chem_pm21['query']['parameters']['value'], valList21[0])
        self.assertEqual(chem_pm21['query']['parameters']['descriptor_type'], search_type2)
        self.assertEqual(chem_pm21['query']['parameters']['match_type'], 'graph-exact')

        # test multiple SMILES
        valList22 = ["c1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO[P@@](=O)(O)O)N",
                     "c1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@H]4[C@H](O3)CO[P@](=O)(O4)O)O)N"]
        chem_pm22 = self.rcsb_util._create_chem_params(search_type2, valList22)

        self.assertEqual(chem_pm22['query']['type'], 'group')
        self.assertEqual(chem_pm22['query']['logical_operator'], 'or')
        self.assertEqual(chem_pm22['request_options'], expected_chem_pm2['request_options'])
        ret_nodelist = chem_pm22['query']['nodes']
        self.assertEqual(len(ret_nodelist), 2)
        self.assertEqual(ret_nodelist[0]['type'], 'terminal')
        self.assertEqual(ret_nodelist[0]['service'], 'chemical')
        self.assertEqual(ret_nodelist[0]['parameters']['type'], 'descriptor')
        self.assertEqual(ret_nodelist[0]['parameters']['value'], valList22[0])
        self.assertEqual(ret_nodelist[0]['parameters']['descriptor_type'], search_type2)
        self.assertEqual(ret_nodelist[0]['parameters']['match_type'], 'graph-exact')
        self.assertEqual(ret_nodelist[1]['parameters']['value'], valList22[1])

    @unittest.skip('test_run_rcsb_search_seq')
    def test_run_rcsb_search_seq(self):
        srch_type = "sequence"
        val_list1 = [self.inputJsonObj1['sequence'][0]]
        json_qry_obj1 = self.rcsb_util._create_seq_ec_uniprot_params(srch_type, val_list1)
        ret_list1, score_dict1 = self.rcsb_util._run_rcsb_search(json_qry_obj1)
        exp_list1 = ['1A62', '1A63', '1A8V', '1PV4', '1PVO', '1XPO', '1XPR', '1XPU', '2A8V',
                     '2HT1', '3ICE', '5JJI', '5JJK', '5JJL', '6DUQ', '6WA8', '6XAS', '6XAV',
                     '6Z9P', '6Z9Q', '6Z9R', '6Z9S', '6Z9T', '7ADB', '7ADC', '7ADD', '7ADE',
                     '7X2R', '8E3H', '8E5L', '8E5P', '8E6W', '8E70']
        # Use assertGreaterEqual() to make sure the test is True even after the remote db expands
        if ret_list1:
            self.assertGreaterEqual(set(ret_list1), set(exp_list1))
            self.assertCountEqual(score_dict1.keys(), ret_list1)

        val_list2 = [self.inputJsonObj1['sequence'][1]]
        json_qry_obj2 = self.rcsb_util._create_seq_ec_uniprot_params(srch_type, val_list2)
        ret_list2, score_dict2 = self.rcsb_util._run_rcsb_search(json_qry_obj2)
        exp_list2 = ['1A69', '1ECP', '1K9S', '1OTX', '1OTY', '1OU4', '1OUM', '1OV6', '1OVG', '1PK7',
                     '1PK9', '1PKE', '1PR0', '1PR1', '1PR2', '1PR4', '1PR5', '1PR6', '1PW7', '3OCC',
                     '3ONV', '3OOE', '3OOH', '3OPV', '3UT6', '4RJ2', '4TS3', '4TS9', '4TTA', '4TTI',
                     '4TTJ', '5I3C', '5IU6', '6XZ2']
        if ret_list2:
            self.assertGreaterEqual(set(ret_list2), set(exp_list2))
            self.assertCountEqual(score_dict2.keys(), ret_list2)

        # test the logical 'or' sequence search
        val_list1.extend(val_list2)
        json_qry_obj3 = self.rcsb_util._create_seq_ec_uniprot_params(srch_type, val_list1)
        ret_list3, score_dict3 = self.rcsb_util._run_rcsb_search(json_qry_obj3)
        exp_list1.extend(exp_list2)
        if ret_list3:
            # self.assertCountEqual(ret_list3, exp_list1)
            self.assertGreaterEqual(set(ret_list3), set(exp_list1))
            self.assertCountEqual(score_dict3.keys(), ret_list3)

    @unittest.skip('test_run_rcsb_search_seq_ec_uniprot')
    def test_run_rcsb_search_seq_ec_uniprot(self):
        srch_type1 = "sequence"
        val_list1 = self.inputJsonObj1['sequence']
        json_qry_obj1 = self.rcsb_util._create_seq_ec_uniprot_params(srch_type1, val_list1)
        ret_list1, score_dict1 = self.rcsb_util._run_rcsb_search(json_qry_obj1)
        exp_list1 = ['1A62', '1A63', '1A8V', '1PV4', '1PVO', '1XPO', '1XPR', '1XPU', '2A8V', '2HT1',
                     '3ICE', '5JJI', '5JJK', '5JJL', '6DUQ', '6WA8', '6XAS', '6XAV', '6Z9P', '6Z9Q',
                     '6Z9R', '6Z9S', '6Z9T', '7ADB', '7ADC', '7ADD', '7ADE', '7X2R', '8E3H', '8E5L',
                     '8E5P', '8E6W', '8E70', '1A69', '1ECP', '1K9S', '1OTX', '1OTY', '1OU4', '1OUM',
                     '1OV6', '1OVG', '1PK7', '1PK9', '1PKE', '1PR0', '1PR1', '1PR2', '1PR4', '1PR5',
                     '1PR6', '1PW7', '3OCC', '3ONV', '3OOE', '3OOH', '3OPV', '3UT6', '4RJ2', '4TS3',
                     '4TS9', '4TTA', '4TTI', '4TTJ', '5I3C', '5IU6', '6XZ2']
        # Use assertGreaterEqual() to make sure the test is True even after the remote db expands
        if ret_list1:
            self.assertGreaterEqual(set(ret_list1), set(exp_list1))
            self.assertCountEqual(score_dict1.keys(), ret_list1)

        srch_type2 = "uniprot_ids"
        val_list2 = self.inputJsonObj1['uniprot_id']  # ["Q01532", "R9RYW2"]
        json_qry_obj2 = self.rcsb_util._create_seq_ec_uniprot_params(srch_type2, val_list2)
        ret_list2, score_dict2 = self.rcsb_util._run_rcsb_search(json_qry_obj2)
        exp_list2 = ['1A6R', '1GCB', '2DZY', '2DZZ', '2E00', '2E01', '2E02', '2E03',
                     '3GCB', '5FAL', '5FAN']
        if ret_list2:
            # self.assertCountEqual(ret_list2, exp_list2)
            self.assertGreaterEqual(set(ret_list2), set(exp_list2))
            self.assertCountEqual(score_dict2.keys(), ret_list2)

        srch_type3 = "ec_numbers"
        val_list3 = self.inputJsonObj1['ec_number']  # ["3.5.2.6", "3.4.13.9"]
        json_qry_obj3 = self.rcsb_util._create_seq_ec_uniprot_params(srch_type3, val_list3)
        ret_list3, score_dict3 = self.rcsb_util._run_rcsb_search(json_qry_obj3)
        if ret_list3:
            self.assertGreaterEqual(len(ret_list3), 1740)
            self.assertCountEqual(ret_list3, score_dict3.keys())

    @unittest.skip('test_run_rcsb_search_chem')
    def test_run_rcsb_search_chem(self):
        srch_type = "InChI"
        val_list1 = self.inputJsonObj3['InChI']
        json_qry_obj1 = self.rcsb_util._create_chem_params(srch_type, val_list1)
        ret_list1, score_dict1 = self.rcsb_util._run_rcsb_search(json_qry_obj1)

        if ret_list1:
            # Use assertGreaterEqual() to make sure the test is True even after remote dbs expand
            self.assertGreaterEqual(len(ret_list1), 2054)
            self.assertCountEqual(score_dict1.keys(), ret_list1)

        val_list2 = self.inputJsonObj3['SMILES']
        json_qry_obj2 = self.rcsb_util._create_chem_params(srch_type, val_list2)
        ret_list2, score_dict2 = self.rcsb_util._run_rcsb_search(json_qry_obj2)
        if ret_list2:
            self.assertGreaterVyEqual(len(ret_list2), 2185)
            self.assertCountEqual(score_dict2.keys(), ret_list2)

    @unittest.skip('test_get_pdb_ids_by_sequence_uniprot_ec')
    def test_get_pdb_ids_by_sequence_uniprot_ec(self):
        ret = self.rcsb_util._get_pdb_ids(self.inputJsonObj1)
        k1 = 'total_count'
        k2 = 'id_list'
        k3 = 'id_score_dict'
        if ret:
            self.assertIn(k1, ret, f'Key {k1} not in returned object')
            self.assertIn(k2, ret, f'Key {k2} not in returned object')
            self.assertIn(k3, ret, f'Key {k3} not in returned object')
            self.assertEqual(len(ret[k2]), ret[k1])
            print(f'RCSB query by sequence/ecnum/uniprotid returned {ret[k1]} ids.')

    # Not implemented yet, so skipped for now
    @unittest.skip('test_get_pdb_ids_by_source_organism')
    def test_get_pdb_ids_by_source_organism(self):
        ret = self.rcsb_util._get_pdb_ids(self.inputJsonObj2)
        k1 = 'total_count'
        k2 = 'id_list'
        k3 = 'id_score_dict'
        if ret:
            self.assertIn(k1, ret, f'Key {k1} not in returned object')
            self.assertIn(k2, ret, f'Key {k2} not in returned object')
            self.assertIn(k3, ret, f'Key {k3} not in returned object')
            self.assertEqual(len(ret[k2]), ret[k1])
            print(f'RCSB query by source organism returned {ret[k1]} ids.')

    @unittest.skip('test_get_pdb_ids_by_chem')
    def test_get_pdb_ids_by_chem(self):
        ret = self.rcsb_util._get_pdb_ids(self.inputJsonObj3)
        k1 = 'total_count'
        k2 = 'id_list'
        k3 = 'id_score_dict'
        if ret:
            self.assertIn(k1, ret, f'Key {k1} not in returned object')
            self.assertIn(k2, ret, f'Key {k2} not in returned object')
            self.assertIn(k3, ret, f'Key {k3} not in returned object')
            self.assertEqual(len(ret[k2]), ret[k1])
            print(f'RCSB query by chem returned {ret[k1]} ids.')

    @unittest.skip('test_queryGraphql')
    def test_queryGraphql(self):
        id_list = ['1A0I', '1A49']
        gqlData = self.rcsb_util._queryGraphql(id_list).get('data', {})
        gqlqry_entries = gqlData.get('entries', [])

        if gqlqry_entries:
            self.assertEqual(len(gqlqry_entries), len(id_list))

            expected_keys = ['rcsb_id', 'exptl', 'rcsb_primary_citation',
                             'polymer_entities', 'nonpolymer_entities']
            for ent in gqlqry_entries:
                for k in expected_keys:
                    self.assertIn(k, ent)

                if ent.get('polymer_entities', None):
                    poly_en = ent.get('polymer_entities', None)
                    for pe in poly_en:
                        if pe.get('entity_poly', None):
                            self.assertIn('pdbx_seq_one_letter_code', pe['entity_poly'])
                            self.assertIn('pdbx_strand_id', pe['entity_poly'])
                        if pe.get('rcsb_entity_source_organism', None):
                            pe_srcg = pe['rcsb_entity_source_organism']
                            for pes in pe_srcg:
                                self.assertIn('ncbi_taxonomy_id', pes)
                                self.assertIn('ncbi_scientific_name', pes)
                        if pe.get('rcsb_polymer_entity_container_identifiers', None):
                            pe_cont = pe['rcsb_polymer_entity_container_identifiers']
                            if pe_cont.get('reference_sequence_identifiers', None):
                                pe_idfs = pe_cont['reference_sequence_identifiers']
                                for pei in pe_idfs:
                                    self.assertIn('database_accession', pei)
                                    self.assertIn('database_name', pei)

    @unittest.skip('test_formatRCSBJson')
    def test_formatRCSBJson(self):
        formatted_json = self.rcsb_util._formatRCSBJson(self.gql_data['data']['entries'])
        expected_keys = ['1A0I', '1A49']
        self.assertIn(expected_keys[0], formatted_json)
        self.assertIn(expected_keys[1], formatted_json)
        for k in expected_keys:
            struct = formatted_json.get(k, {})
            self.assertIn('method', struct)
            self.assertIn('primary_citation', struct)
            self.assertIn('polymer_entities', struct)
            self.assertIn('nonpolymer_entities', struct)

            if struct.get('polymer_entities', None):
                for pe in struct.get('polymer_entities', None):
                    self.assertIn('id', pe)
                    self.assertIn('uniprot_name', pe)
                    self.assertIn('one_letter_code_sequence', pe)
                    self.assertIn('pdbx_strand_id', pe)
                    self.assertIn('taxonomy', pe)
                    self.assertIn('source_organism', pe)
                    self.assertIn('ec_numbers', pe)
                    self.assertIn('uniprot_ec', pe)
                    self.assertIn('uniprotID', pe)
            if struct.get('nonpolymer_entities', None):
                for npe in struct.get('nonpolymer_entities', None):
                    self.assertIn('pdb_ref_name', npe)
                    self.assertIn('InChI', npe)
                    self.assertIn('InChIKey', npe)
                    self.assertIn('SMILES', npe)

        gql_itms = {}
        for key, val in formatted_json.items():
            gql_itms[key] = val
        self.assertCountEqual(gql_itms['1A0I'], self.gqldata_items['1A0I'])
        self.assertCountEqual(gql_itms['1A49'], self.gqldata_items['1A49'])

    @unittest.skip('test_get_graphql_data')
    def test_get_graphql_data(self):
        id_list = ['1A0I', '1A49', '1A5U', '1A82', '1AQ2']
        gql_itms = self.rcsb_util._get_graphql_data(id_list)
        for k in ['total_count', 'id_list']:
            self.assertIn(k, gql_itms)
        for id in id_list:
            self.assertIn(id, gql_itms['id_list'])
            self.assertIn(id, gql_itms.keys())

        self.assertEqual(gql_itms['total_count'], len(id_list))
        self.assertEqual(gql_itms['id_list'], id_list)
        self.assertCountEqual(gql_itms['1A0I'], self.gqldata_items['1A0I'])
        self.assertCountEqual(gql_itms['1A49'], self.gqldata_items['1A49'])

    @unittest.skip('test_get_cpd_match')
    def test_get_cpd_match(self):
        npes = [
            {"InChI": "InChI=1S/K/q+1",
             "InChIKey": "NPYPAHLBTDXSSS-UHFFFAOYSA-N",
             "SMILES": "[K+]",
             "pdb_ref_name": "MERCURY (II) ION"},
            {"InChI": "InChI=1S/C2H2O4/c3-1(4)2(5)6/h(H,3,4)(H,5,6)/p-2",
             "InChIKey": "MUBZPKHOEPUJKR-UHFFFAOYSA-L",
             "SMILES": "C(=O)(C(=O)[O-])[O-]"},
            {"InChI": "InChI=1S/Mg/q+2",
             "InChIKey": "JLVVSXFLKOJNIY-UHFFFAOYSA-N",
             "SMILES": "[Mg+2]"},
            {"InChI": "InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1",
             "InChIKey": "ZKHQWZAMYRWXGA-KQYNXXCUSA-N",
             "SMILES": "c1nc(c2c(n1)n(cn2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N"},
            {"InChI": "InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-fake",
             "InChIKey": "ZKHQWZAMYRWXGA-KQYNXXFAKE-N",
             "SMILES": "c1nc(c2c(n1)n(cn2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N"}
        ]
        cpds = []
        for npe in npes:
            cpd = self.rcsb_util._get_cpd_match(npe)
            cpds.append(cpd)
        self.assertEqual(cpds[0], 'cpd00205 (K+)')
        self.assertEqual(cpds[1], 'cpd00180 (Oxalate)')
        self.assertEqual(cpds[2], 'cpd00254 (Mg)')
        self.assertEqual(cpds[3], 'cpd00002 (ATP)')
        self.assertEqual(cpds[4], 'InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-fake')

    @unittest.skip('test_blastRCSBSequence')
    def test_blastRCSBSequence(self):
        in_seq1 = 'VNIKTNPFKAVSFVESAIKKALDNAGYLIAEIKYDGVRGNICVDNTANSYWLSRVSKTIPALEHLNGFDVRWKRLLNDDRCFYKDGFMLDGELMVKGVDFNTGSGLLRTKWTDTKNQEFHEELFVEPIRKKDKVPFKLHTGHLHIKLYAILPLHIVESGEDCDVMTLLMQEHVKNMLPLLQEYFPEIEWQAAESYEVYDMVELQQLYEQKRAEGHEGLIVKDPMCIYKRGKKSGWWKMKPENEADGIIQGLVWGTKGLANEGKVIGFEVLLESGRLVNATNISRALMDEFTETVKEATLSQWGFFSPYGIGDNDACTINPYDGWACQISYMEETPDGSLRHPSFVMFR'
        e_val = 1e-3
        iden = 0.75
        exp_keys = ['rcsbid', 'name', 'pdbx_sequence', 'rcsbec', 'uniprotec', 'identity',
                    'uniprotID', 'pdbx_strand_id', 'taxonomy', 'evalue', 'method', 'components',
                    'references']
        gql_itms1 = self.rcsb_util._blastRCSBSequence(in_seq1, self.gql_data, e_val, iden)
        for k in exp_keys:
            self.assertIn(k, gql_itms1[0])

        in_seq2 = 'SKSHSEAGSAFIQTQQLHAAMADTFLEHMCRLDIDSAPITARNTGIICTIGPASRSVETLKEMIKSGMNVARMNFSHGTHEYHAETIKNVRTATESFASDPILYRPVAVALDTKGPEIRTGLIKGSGTAEVELKKGATLKITLDNAYMEKCDENILWLDYKNICKVVDVGSKVYVDDGLISLQVKQKGPDFLVTEVENGGFLGSKKGVNLPGAAVDLPAVSEKDIQDLKFGVEQDVDMVFASFIRKAADVHEVRKILGEKGKNIKIISKIENHEGVRRFDEILEASDGIMVARGDLGIEIPAEKVFLAQKMIIGRCNRAGKPVICATQMLESMIKKPRPTRAEGSDVANAVLDGADCIMLSGETAKGDYPLEAVRMQHLIAREAEAAMFHRKLFEELARSSSHSTDLMEAMAMGSVEASYKCLAAALIVLTESGRSAHQVARYRPRAPIIAVTRNHQTARQAHLYRGIFPVVCKDPVQEAWAEDVDLRVNLAMNVGKARGFFKKGDVVIVLTGWRPGSGFTNTMRVVPVP'
        gql_itms2 = self.rcsb_util._blastRCSBSequence(in_seq2, self.gql_data, e_val, iden)
        for k in exp_keys:
            self.assertIn(k, gql_itms2[0])

        in_seq3 = 'SKSHSEAGSAFIQTQQLHAAMADTFLEHMCRLDIDSAPITARNTGIICTIGPASRSVETLKEMIKSGMNVARMNFSHGTHE'
        gql_itms3 = self.rcsb_util._blastRCSBSequence(in_seq3, self.gql_data, e_val, iden)
        self.assertTrue(gql_itms3)
        for k in exp_keys:
            self.assertIn(k, gql_itms3[0])

    @unittest.skip('test_get_graphql_data_with_cutoffs')
    def test_get_graphql_data_with_cutoffs(self):
        id_list = ['1A0I', '1A49', '1A5U', '1A82', '1AQ2']
        in_seq = 'VNIKTNPFKAVSFVESAIKKALDNAGYLIAEIKYDGVRGNICVDNTANSYWLSRVSKTIPALEHLNGFDVRWKRLLNDDRCFYKDGFMLDGELMVKGVDFNTGSGLLRTKWTDTKNQEFHEELFVEPIRKKDKVPFKLHTGHLHIKLYAILPLHIVESGEDCDVMTLLMQEHVKNMLPLLQEYFPEIEWQAAESYEVYDMVELQQLYEQKRAEGHEGLIVKDPMCIYKRGKKSGWWKMKPENEADGIIQGLVWGTKGLANEGKVIGFEVLLESGRLVNATNISRALMDEFTETVKEATLSQWGFFSPYGIGDNDACTINPYDGWACQISYMEETPDGSLRHPSFVMFR'
        eval = 1e-3
        iden = 0.75
        gql_itms = self.rcsb_util._get_graphql_data_with_cutoffs(id_list, in_seq, eval, iden)
        self.assertIn('id_list', gql_itms)
        self.assertIn(in_seq, gql_itms)

    @unittest.skip('test_write_struct_info')
    def test_write_struct_info(self):
        tbody_html = self.rcsb_util._write_struct_info(self.gqldata_items)
        self.assertEqual(tbody_html.count('</tr>'), 2)
        self.assertIn('1A0I', tbody_html)
        self.assertIn('1A49', tbody_html)

    @unittest.skip('test_query_structure_info')
    def test_query_structure_info(self):
        params1 = {
            'workspace_name': self.wsName,
            'sequence_strings': self.inputJsonObj1['sequence'],
            'uniprot_ids': self.inputJsonObj1['uniprot_id'],
            'ec_numbers': self.inputJsonObj1['ec_number'],
            'inchis': self.inputJsonObj3['InChI'],
            'smiles': self.inputJsonObj3['SMILES'],
            'evalue_cutoff': 1e-10,
            'identity_cutoff': 0.9,
            'logical_and': 0
        }

        struct_ret = self.rcsb_util.query_structure_info(params1)
        if struct_ret:
            rcsb_ids = struct_ret.get('rcsb_ids', [])
            self.assertEqual(len(rcsb_ids), 4045)
            self.assertIn('rcsb_ids', struct_ret)
            self.assertIn('report_name', struct_ret)
            self.assertIn('report_ref', struct_ret)

        params2 = {
            'workspace_name': self.wsName,
            'sequence_strings': self.inputJsonObj1['sequence'],
            'evalue_cutoff': 1e-10,
            'identity_cutoff': 0.9,
            'logical_and': 1
        }
        struct_ret_and = self.rcsb_util.query_structure_info(params2)
        self.assertEqual(len(struct_ret_and.get('rcsb_ids', [])), 67)

    @unittest.skip('test_validate_import_rcsb_params')
    def test_validate_import_rcsb_params(self):
        params1 = {
            'workspace_name': self.wsName,
            'structures_name': 'structs_name',
            'rcsb_infos': [{
                'rcsb_id': '1A0I',
                'extension': '.pdb',
                'narrative_id': 555,
                'genome_name': 'ATCC_49442',
                'feature_id': 'CDS.123',
                'is_model': 1
            }, {
                'rcsb_id': '1fat',
                'extension': 'cif',
                'narrative_id': 555,
                'genome_name': 'ATCC_49442',
                'feature_id': 'CDS.133',
                'is_model': 1
            }, {
                'rcsb_id': '1A49',
                'extension': 'abcd',
                'narrative_id': 555,
                'genome_name': 'ATCC_49442',
                'feature_id': 'CDS.143',
                'is_model': 1
            }]
        }
        self.assertEqual(len(params1['rcsb_infos']), 3)
        params2 = self.rcsb_util._validate_import_rcsb_params(params1)
        self.assertEqual(len(params2['rcsb_infos']), 2)

    #@unittest.skip('test_upload_rcsbs_1')
    def test_upload_rcsbs_1(self):
        id_list = ['6ifs', '6ift', '6ifw', '1fat']
        params = {
            'workspace_name': self.wsName,
            'structures_name': 'structs_name',
            'rcsb_infos': [{
                'rcsb_id': id_list[0],
                'extension': '.pdb',
                'narrative_id': 57196,
                'genome_name': 'Synthetic_bacterium_JCVI_Syn3_genome',
                'feature_id': 'JCVISYN3_0004',
                'is_model': 1
            }, {
                'rcsb_id': id_list[1],
                'extension': '.pdb',
                'narrative_id': 57196,
                'genome_name': 'Synthetic_bacterium_JCVI_Syn3_genome',
                'feature_id': 'JCVISYN3_0004',
                'is_model': 1
            }, {
                'rcsb_id': id_list[2],
                'extension': '.pdb',
                'narrative_id': 57196,
                'genome_name': 'Synthetic_bacterium_JCVI_Syn3_genome',
                'feature_id': 'JCVISYN3_0004',
                'is_model': 1
            }, {
                'rcsb_id': id_list[3],
                'extension': '.cif',
                'narrative_id': 57196,
                'genome_name': 'Synthetic_bacterium_JCVI_Syn3_genome',
                'feature_id': 'JCVISYN3_0001',
                'is_model': 0
            }]
        }
        pdb_objects = list()
        pdb_infos = list()
        successful_ids = list()
        skipped_ids = list()
        pdb_objects, pdb_infos, successful_ids, skipped_ids = self.rcsb_util.upload_rcsbs(
                                                                    params, self.wsName)
        self.assertTrue(pdb_infos)
        self.assertEqual(len(pdb_infos), 3)

    @unittest.skip('test_upload_rcsbs_2')
    def test_upload_rcsbs_2(self):
        id_list = ['1A0I', '1A49', '1A5U', '1A82', '1AQ2']
        params = {
            'workspace_name': self.wsName,
            'structures_name': 'structs_name',
            'rcsb_infos': [{
                'rcsb_id': id_list[0],
                'extension': '.pdb',
                'narrative_id': 63679,
                'genome_name': 'MLuteus_ATCC_49442',
                'feature_id': 'MLuteus_masurca_RAST.CDS.133',
                'is_model': 1
            }, {
                'rcsb_id': id_list[1],
                'extension': '.pdb',
                'narrative_id': 63679,
                'genome_name': 'MLuteus_ATCC_49442',
                'feature_id': 'MLuteus_masurca_RAST.CDS.133',
                'is_model': 1
            }]
        }
        pdb_objects = list()
        pdb_infos = list()
        successful_ids = list()
        skipped_ids = list()
        pdb_objects, pdb_infos, successful_ids, skipped_ids = self.rcsb_util.upload_rcsbs(
                                                                    params, self.wsName)
        self.assertFalse(pdb_infos)

    # Testing self.serviceImpl functions
    @unittest.skip('test_Impl_query_rcsb_structures')
    def test_Impl_query_rcsb_structures(self):
        params = {
            'workspace_name': self.wsName,
            'sequence_strings': self.inputJsonObj1['sequence'],
            'evalue_cutoff': 1e-10,
            'identity_cutoff': 0.9,
            'logical_and': 0
        }
        qry_ret1 = self.serviceImpl.query_rcsb_structures(self.ctx, params)
        if qry_ret1:
            self.assertCountEqual(qry_ret1[0].keys(),
                                  ['rcsb_ids', 'rcsb_scores', 'report_ref', 'report_name'])
            self.assertEqual(len(qry_ret1[0]['rcsb_ids']), 67)

        params['uniprot_ids'] = self.inputJsonObj1['uniprot_id']
        params['ec_numbers'] = self.inputJsonObj1['ec_number']
        params['inchis'] = self.inputJsonObj3['InChI']
        params['smiles'] = self.inputJsonObj3['SMILES']

        qry_ret2 = self.serviceImpl.query_rcsb_structures(self.ctx, params)
        if qry_ret2:
            self.assertCountEqual(qry_ret2[0].keys(),
                                  ['rcsb_ids', 'rcsb_scores', 'report_ref', 'report_name'])
            self.assertEqual(len(qry_ret2[0]['rcsb_ids']), 4045)

    @unittest.skip('test_Impl_query_rcsb_annotations')
    def test_Impl_query_rcsb_annotations(self):
        params = {
            'workspace_name': self.wsName,
            'sequence_strings': self.inputJsonObj1['sequence'],
            'evalue_cutoff': 1e-20,
            'identity_cutoff': 0.999,
            'logical_and': 0
        }
        qry_ret = self.serviceImpl.query_rcsb_annotations(self.ctx, params)
        if qry_ret:
            self.assertCountEqual(qry_ret[0].keys(), params['sequence_strings'])
