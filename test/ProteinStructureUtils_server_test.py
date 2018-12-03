# -*- coding: utf-8 -*-
import os
import shutil
import time
import unittest
from configparser import ConfigParser

from ProteinStructureUtils.ProteinStructureUtilsImpl import ProteinStructureUtils
from ProteinStructureUtils.ProteinStructureUtilsServer import MethodContext
from ProteinStructureUtils.authclient import KBaseAuth as _KBaseAuth

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


class ProteinStructureUtilsTest(unittest.TestCase):

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
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        cls.dfu = DataFileUtil(cls.callback_url)
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ProteinStructureUtils_" + str(suffix)
        cls.ws_id = cls.wsClient.create_workspace({'workspace': cls.wsName})[0]
        cls.prepareData()

    @classmethod
    def prepareData(cls):
        file = '5o5y.pdb'
        cls.pdb_file_path = os.path.join(cls.scratch, file)
        shutil.copy(os.path.join('data', file), cls.pdb_file_path)
        file_to_shock_params = {
            'file_path': cls.pdb_file_path,
            'pack': 'gzip',
            'make_handle': True,
        }
        shock_id = cls.dfu.file_to_shock(file_to_shock_params)['handle']['hid']
        data = {
            'name': '',
            'num_chains': 0,
            'num_residues': 0,
            'num_atoms': 0,
            'protein': {
                'id': '',
                'sequence': '',
                'md5': ''
            },
            'user_data': '',
            'pdb_handle': shock_id,
        }
        info = cls.dfu.save_objects({
            'id': cls.ws_id,
            'objects': [
                {'type': 'KBaseStructure.ModelProteinStructure',
                 'name': file,
                 'data': data}]
        })[0]
        cls.pdb_ref = f"{info[6]}/{info[0]}/{info[4]}"


    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def test_model_upload(self):
        ret = self.serviceImpl.import_model_pdb_file(
            self.ctx, {
                'input_file_path': self.pdb_file_path,
                'structure_name': 'import_test',
                'workspace_name': self.wsName,
            })[0]
        self.assertCountEqual(ret.keys(), ["structure_obj_ref", "report_ref", "report_name"])

    def test_structure_to_pdb_file(self):
        ret = self.serviceImpl.structure_to_pdb_file(self.ctx, {'input_ref': self.pdb_ref,
                                                                'destination_dir': self.scratch})

    def test_export_structure(self):
        ret = self.serviceImpl.export_pdb(self.ctx, {'input_ref': self.pdb_ref})
