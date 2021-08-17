# -*- coding: utf-8 -*-
import os
import shutil
import time
import unittest
from Bio import PDB
from configparser import ConfigParser

from ProteinStructureUtils.ProteinStructureUtilsImpl import ProteinStructureUtils
from ProteinStructureUtils.ProteinStructureUtilsServer import MethodContext
from ProteinStructureUtils.authclient import KBaseAuth as _KBaseAuth

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace

from ProteinStructureUtils.Utils.PDBUtils import PDBUtil

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
        cls.pdb_util = PDBUtil(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        cls.dfu = DataFileUtil(cls.callback_url)
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ProteinStructureUtils_" + str(suffix)
        cls.ws_id = cls.wsClient.create_workspace({'workspace': cls.wsName})[0]
        cls.prepareData()
        cls.prepareCIFData()

    @classmethod
    def prepareData(cls):
        file = '1nqg.pdb'
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
            'proteins': [{
                'id': '',
                'sequence': '',
                'md5': ''
            }],
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
    def prepareCIFData(cls):
        file_name = '1fat.cif'
        cls.pdb_mmcif_file_path = os.path.join(cls.scratch, file_name)
        shutil.copy(os.path.join('data', file_name), cls.pdb_mmcif_file_path)
        file_to_shock_params = {
            'file_path': cls.pdb_mmcif_file_path,
            'pack': 'gzip',
            'make_handle': True,
        }

        cpd = {'misc': '',
               'molecule': 'adp-dependent glucokinase,adp-dependent glucokinase,adp-dependent glucokinase ',
               'chain': 'a, b',
               'synonym': 'adpgk,adpgk,adpgk',
               'ec_number': '2.7.1.147',
               'ec': '2.7.1.147,2.7.1.147,2.7.1.147',
               'engineered': 'yes'}
        src = {'misc': '',
               'organism_scientific': 'thermococcus litoralis (strain atcc 51850 / dsm5473 / jcm 8560 / ns-c) ',
               'organism_taxid': '523849',
               'gene': 'glka, occ_09701',
               'expression_system': 'escherichia coli bl21(de3)',
               'expression_system_taxid': '469008',
               'expression_system_vector_type': 'plasmid',
               'expression_system_vector': 'pet17'}

        shock_id = cls.dfu.file_to_shock(file_to_shock_params)['handle']['hid']
        data = {
            'name': 'PHYTOHEMAGGLUTININ-L',
            'rcsb_id': 'unknown',
            'deposition_date': '1996-06-12',
            'head': 'LECTIN',
            'release_date': 'unknown',
            'structure_method': 'x-ray diffraction',
            'resolution': 2.8,
            'author': 'unknown',
            'compound': cpd,
            'source': src,
            'num_chains': 0,
            'num_atoms': 0,
            'num_models': 0,
            'num_residues': 0,
            'num_het_atoms': 0,
            'num_water_atoms': 0,
            'num_disordered_atoms': 0,
            'num_disordered_residues': 0,
            'proteins': [{
                'id': '',
                'sequence': '',
                'md5': ''
            }],
            'user_data': '',
            'pdb_handle': shock_id,
            'mmcif_handle': shock_id,
            'xml_handle': shock_id
        }
        info = cls.dfu.save_objects({
            'id': cls.ws_id,
            'objects': [
                {'type': 'KBaseStructure.ExperimentalProteinStructure',
                 'name': file_name,
                 'data': data}]
        })[0]
        cls.pdb_mmCif_ref = f"{info[6]}/{info[0]}/{info[4]}"

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    # To be called when need to download a PDB file
    @classmethod
    def downloadPDB(cls, file_name, fmt='pdb'):
        """ download structures from the PDB database
            e.g., file_name='1fat.cif', format='mmCif'
        """
        pdbl = PDB.PDBList()
        return pdbl.retrieve_pdb_file(file_name, pdir='data', format=fmt)

    # @unittest.skip('test_model_upload')
    def test_model_upload(self):
        ret = self.serviceImpl.import_model_pdb_file(
            self.ctx, {
                'input_file_path': self.pdb_file_path,
                'structure_name': 'import_model_pdb_test',
                'workspace_name': self.wsName,
            })[0]
        self.assertCountEqual(ret.keys(), ["structure_obj_ref", "report_ref", "report_name"])

    # @unittest.skip('test_model_upload2')
    def test_model_upload2(self):
        ret = self.serviceImpl.import_model_pdb_file(
            self.ctx, {
                'input_file_path': self.pdb_file_path,
                'structure_name': 'import_model_pdb_test2',
                'workspace_name': self.wsName,
            })[0]
        self.assertCountEqual(ret.keys(), ["structure_obj_ref", "report_ref", "report_name"])

    # @unittest.skip('test_structure_to_pdb_file')
    def test_structure_to_pdb_file(self):
        ret = self.serviceImpl.structure_to_pdb_file(self.ctx, {'input_ref': self.pdb_ref,
                                                                'destination_dir': self.scratch})
        self.assertEqual(ret[0]['file_path'], os.path.join(self.scratch, '1nqg.pdb'))

    # @unittest.skip('test_export_pdb_structure')
    def test_export_pdb_structure(self):
        ret = self.serviceImpl.export_pdb(self.ctx, {'input_ref': self.pdb_ref})
        self.assertCountEqual(ret[0].keys(), ['shock_id'])

    # @unittest.skip('experiment_upload')
    def test_experiment_upload(self):
        ret = self.serviceImpl.import_experiment_pdb_file(
            self.ctx, {
                'input_file_path': self.pdb_mmcif_file_path,
                'structure_name': 'import_mmcif_test',
                'workspace_name': self.wsName,
            })[0]
        self.assertCountEqual(ret.keys(), ["structure_obj_ref", "report_ref", "report_name"])

    # @unittest.skip('test_structure_to_mmcif_file')
    def test_structure_to_mmcif_file(self):
        ret = self.serviceImpl.structure_to_pdb_file(self.ctx, {'input_ref': self.pdb_mmCif_ref,
                                                                'destination_dir': self.scratch})
        self.assertEqual(ret[0]['file_path'], os.path.join(self.scratch, '1fat.cif'))

    # @unittest.skip('test_export_structure')
    def test_export_mmcif_structure(self):
        ret = self.serviceImpl.export_pdb(self.ctx, {'input_ref': self.pdb_mmCif_ref})
        self.assertCountEqual(ret[0].keys(), ['shock_id'])

    # Testing PDBUtils module functions
    # @unittest.skip('test_model_file_to_data')
    def test_model_file_to_data(self):
        fileName = '1fat.pdb'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)
        (data1, pp_no1) = self.pdb_util._model_file_to_data(pdb_file_path)
        self.assertEqual(pp_no1, 7)
        self.assertEqual(data1['name'], 'phytohemagglutinin-l')
        self.assertEqual(data1['num_chains'], 4)
        self.assertEqual(data1['num_residues'], 928)
        self.assertEqual(data1['num_atoms'], 7248)
        self.assertCountEqual(data1['proteins'][0].keys(), ['id', 'sequence', 'md5'])
        self.assertCountEqual(data1['compound'].keys(), ['misc', 'molecule', 'chain', 'synonym'])
        self.assertCountEqual(data1['source'].keys(),
                              ['misc', 'organism_scientific',
                               'organism_taxid', 'organ', 'other_details'])
        fileName = '5o5y.pdb'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)
        (data2, pp_no2) = self.pdb_util._model_file_to_data(pdb_file_path)

        self.assertEqual(pp_no2, 2)
        self.assertEqual(data2['name'],
            'crystal structure of thermococcus litoralis adp-dependent glucokinase (gk)')
        self.assertEqual(data2['num_chains'], 2)
        self.assertEqual(data2['num_residues'], 903)
        self.assertEqual(data2['num_atoms'], 8013)
        self.assertCountEqual(data2['proteins'][0].keys(), ['id', 'sequence', 'md5'])
        self.assertCountEqual(data2['compound'].keys(), ['misc', 'molecule', 'chain', 'synonym',
                                                         'ec_number', 'ec', 'engineered'])
        self.assertCountEqual(data2['source'].keys(),
                              ['misc', 'organism_scientific', 'expression_system',
                               'gene', 'expression_system_taxid', 'expression_system_vector_type',
                               'organism_taxid', 'expression_system_vector'])

    # @unittest.skip('test_exp_file_to_data')
    def test_exp_file_to_data(self):
        fileName = '1fat.cif'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)
        (data, pp_no) = self.pdb_util._exp_file_to_data(pdb_file_path)
        self.assertEqual(pp_no, 7)
        self.assertEqual(data['name'], 'PHYTOHEMAGGLUTININ-L')
        self.assertEqual(data['head'], 'LECTIN')
        self.assertEqual(data['rcsb_id'], '')
        self.assertEqual(data['deposition_date'], '1996-06-12')
        self.assertEqual(data['release_date'], '')
        self.assertEqual(data['structure_method'], 'X-RAY DIFFRACTION')
        self.assertEqual(data['resolution'], 2.8)
        self.assertCountEqual(data['proteins'][0].keys(), ['id', 'sequence', 'md5'])
        self.assertEqual(data['compound'], {})
        self.assertEqual(data['source'], {})
        self.assertTrue('pdb_handle' in data.keys())

