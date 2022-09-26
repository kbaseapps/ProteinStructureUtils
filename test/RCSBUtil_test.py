# -*- coding: utf-8 -*-
import os
import shutil
import time
import json
import unittest
from configparser import ConfigParser
#from mock import patch
from unittest.mock import Mock, patch

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

    # Testing RCSBUtil module functions
    #@unittest.skip('test_get_pdb_ids_by_sequence_uniprot_ec')
    def test_get_pdb_ids_by_sequence_uniprot_ec(self):
        ret = self.rcsb_util._get_pdb_ids(self.inputJsonObj1)
        k1 = 'total_count'
        k2 = 'id_list'
        self.assertIn(k1, ret, f'Key {k1} not in returned object')
        self.assertIn(k2, ret, f'Key {k2} not in returned object')
        self.assertEqual(len(ret[k2]), ret[k1])
        print(f'RCSB query returned ret{k1} ids.')

    @unittest.skip('test_get_pdb_ids_by_source_organism')
    def test_get_pdb_ids_by_source_organism(self):
        ret = self.rcsb_util._get_pdb_ids(self.inputJsonObj2)
        k1 = 'total_count'
        k2 = 'id_list'
        self.assertIn(k1, ret, f'Key {k1} not in returned object')
        self.assertIn(k2, ret, f'Key {k2} not in returned object')
        self.assertEqual(len(ret[k2]), ret[k1])
        print(f'RCSB query returned ret{k1} ids.')

    #@unittest.skip('test_get_pdb_ids_by_chem')
    def test_get_pdb_ids_by_chem(self):
        ret = self.rcsb_util._get_pdb_ids(self.inputJsonObj3)
        k1 = 'total_count'
        k2 = 'id_list'
        self.assertIn(k1, ret, f'Key {k1} not in returned object')
        self.assertIn(k2, ret, f'Key {k2} not in returned object')
        self.assertEqual(len(ret[k2]), ret[k1])
        print(f'RCSB query returned ret{k1} ids.')

    def test_get_graphql_data(self):
        id_list = ['1A0I', '1A49', '1A5U']  # , '1A82', '1AQ2']
        gql_ret = self.rcsb_util._get_graphql_data(id_list)

        expected_keys = ['total_count', 'id_list']
        expected_keys.extend(id_list)
        print(expected_keys)
        for k in expected_keys:
            self.assertIn(k, gql_ret)

        self.assertEqual(gql_ret['total_count'], len(id_list))
        self.assertEqual(gql_ret['id_list'], id_list)


