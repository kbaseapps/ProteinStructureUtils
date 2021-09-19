# -*- coding: utf-8 -*-
import os
import shutil
import time
import unittest
from configparser import ConfigParser
from mock import patch

from ProteinStructureUtils.ProteinStructureUtilsImpl import ProteinStructureUtils
from ProteinStructureUtils.ProteinStructureUtilsServer import MethodContext
from ProteinStructureUtils.authclient import KBaseAuth as _KBaseAuth

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace

from ProteinStructureUtils.Utils.PDBUtils import PDBUtil

from Bio import PDB
from Bio.PDB.Polypeptide import PPBuilder

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
        cls.cfg['USER_ID'] = cls.ctx['user_id']
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

    def mock_download_staging_file(params):
        print('Mocking DataFileUtilClient.download_staging_file')
        print(params)

        staging_file_name = params.get('staging_file_subdir_path')
        file_name = os.path.basename(staging_file_name)
        staging_file_subdir_path = os.path.join('/kb/module/test/data', file_name)
        file_path = os.path.join('/kb/module/work/tmp', file_name)
        shutil.copy(staging_file_subdir_path, file_path)

        return {'copy_file_path': file_path}

    # Testing self.serviceImpl functions
    @unittest.skip('test_model_upload1')
    def test_model_upload1(self):
        ret = self.serviceImpl.import_model_pdb_file(
            self.ctx, {
                'input_file_path': self.pdb_file_path,
                'structure_name': 'import_model_pdb_test1',
                'workspace_name': self.wsName,
            })[0]
        self.assertCountEqual(ret.keys(), ["structure_obj_ref", "report_ref", "report_name"])

    @unittest.skip('test_model_upload2')
    def test_model_upload2(self):
        fileName = '1fat.pdb'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)
        ret = self.serviceImpl.import_model_pdb_file(
            self.ctx, {
                'input_file_path': pdb_file_path,
                'structure_name': 'import_model_pdb_test2',
                'workspace_name': self.wsName,
            })[0]
        self.assertCountEqual(ret.keys(), ["structure_obj_ref", "report_ref", "report_name"])

        # Check saved object with the ref value of ret['structure_obj_ref'] against datatype
        pdb_obj_data = self.wsClient.get_objects2(
            {'objects': [{'ref': ret['structure_obj_ref']}]})['data'][0]['data']
        self.assertEqual(pdb_obj_data['name'], 'phytohemagglutinin-l')
        self.assertEqual(pdb_obj_data['num_chains'], 4)
        self.assertEqual(pdb_obj_data['num_residues'], 928)
        self.assertEqual(pdb_obj_data['num_atoms'], 7248)
        self.assertCountEqual(pdb_obj_data['proteins'][0].keys(),
                              ['id', 'sequence', 'md5', 'model_id', 'chain_id'])
        self.assertCountEqual(pdb_obj_data['compound'].keys(),
                              ['misc', 'molecule', 'chain', 'synonym'])
        self.assertEqual(pdb_obj_data['compound']['molecule'], 'phytohemagglutinin-l')
        self.assertEqual(pdb_obj_data['compound']['chain'], 'a, b, c, d')
        self.assertEqual(pdb_obj_data['compound']['synonym'],
                         'leucoagglutinating phytohemagglutinin, pha-l')
        self.assertCountEqual(pdb_obj_data['source'].keys(),
                              ['misc', 'organism_scientific',
                               'organism_taxid', 'organ', 'other_details'])
        self.assertEqual(pdb_obj_data['source']['organism_scientific'], 'phaseolus vulgaris')
        self.assertEqual(pdb_obj_data['source']['organism_taxid'], '3885')
        self.assertEqual(pdb_obj_data['source']['organ'], 'seed')
        self.assertEqual(pdb_obj_data['source']['other_details'],
                         'purified pha-l was purchased from sigma')

    @unittest.skip('test_model_upload3')
    def test_model_upload3(self):
        fileName = '6ift.pdb'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)
        ret = self.serviceImpl.import_model_pdb_file(
            self.ctx, {
                'input_file_path': pdb_file_path,
                'structure_name': 'import_model_pdb_6ift',
                'workspace_name': self.wsName,
            })[0]
        self.assertCountEqual(ret.keys(), ["structure_obj_ref", "report_ref", "report_name"])

        # Check saved object with the ref value of ret['structure_obj_ref'] against datatype
        pdb_obj_data = self.wsClient.get_objects2(
            {'objects': [{'ref': ret['structure_obj_ref']}]})['data'][0]['data']
        self.assertEqual(pdb_obj_data['name'], 'ksga from bacillus subtilis in complex with sam')
        self.assertEqual(pdb_obj_data['num_chains'], 1)
        self.assertEqual(pdb_obj_data['num_residues'], 292)
        self.assertEqual(pdb_obj_data['num_atoms'], 2508)
        self.assertCountEqual(pdb_obj_data['proteins'][0].keys(),
                              ['id', 'sequence', 'md5', 'model_id', 'chain_id'])
        self.assertTrue(set(pdb_obj_data['compound'].keys()).issuperset(
                        set(['misc', 'molecule', 'chain', 'synonym'])))
        self.assertIn('ec', pdb_obj_data['compound'].keys())
        self.assertIn('ec_number', pdb_obj_data['compound'].keys())
        self.assertIn('engineered', pdb_obj_data['compound'].keys())
        self.assertEqual(pdb_obj_data['compound']['molecule'],
                         'ribosomal rna small subunit methyltransferase a')
        self.assertEqual(pdb_obj_data['compound']['chain'], 'a')
        self.assertEqual(set(pdb_obj_data['compound']['synonym'].split(',')),
                         set(["16s rrna (adenine(1518)-n(6)/adenine(1519)-n(6))-dimethyltransferase",
                              "16s rrna dimethyladenosine transferase",
                              "16s rrna dimethylase",
                              "s-adenosylmethionine-6-n'",
                              "n'-adenosyl(rrna) dimethyltransferase "]))
        self.assertCountEqual(pdb_obj_data['source'].keys(),
                              ['misc', 'organism_scientific', 'strain', 'gene',
                               'organism_taxid', 'expression_system', 'expression_system_taxid'])
        self.assertEqual(pdb_obj_data['source']['organism_scientific'],
                         'bacillus subtilis (strain 168)')
        self.assertEqual(pdb_obj_data['source']['organism_taxid'], '224308')
        self.assertEqual(pdb_obj_data['source']['gene'], 'rsma, ksga, bsu00420')
        self.assertEqual(pdb_obj_data['source']['expression_system'], 'escherichia coli')
        self.assertEqual(pdb_obj_data['source']['expression_system_taxid'], '562')
        self.assertEqual(pdb_obj_data['source']['strain'], '168')

    @unittest.skip('test_structure_to_pdb_file')
    def test_structure_to_pdb_file(self):
        ret = self.serviceImpl.structure_to_pdb_file(self.ctx, {'input_ref': self.pdb_ref,
                                                                'destination_dir': self.scratch})
        self.assertEqual(ret[0]['file_path'], os.path.join(self.scratch, '1nqg.pdb'))

    @unittest.skip('test_export_pdb_structure')
    def test_export_pdb_structure(self):
        ret = self.serviceImpl.export_pdb(self.ctx, {'input_ref': self.pdb_ref})
        self.assertCountEqual(ret[0].keys(), ['shock_id'])

    @unittest.skip('experiment_upload')
    def test_experiment_upload(self):
        ret = self.serviceImpl.import_experiment_pdb_file(
            self.ctx, {
                'input_file_path': self.pdb_mmcif_file_path,
                'structure_name': 'import_mmcif_test',
                'workspace_name': self.wsName,
            })[0]
        self.assertCountEqual(ret.keys(), ["structure_obj_ref", "report_ref", "report_name"])

        # Check saved object with the ref value of ret['structure_obj_ref'] against datatype
        pdb_obj_data = self.wsClient.get_objects2(
            {'objects': [{'ref': ret['structure_obj_ref']}]})['data'][0]['data']
        self.assertEqual(pdb_obj_data['name'], 'PHYTOHEMAGGLUTININ-L')
        self.assertEqual(pdb_obj_data['head'], 'LECTIN')
        self.assertEqual(pdb_obj_data['deposition_date'], '1996-06-12')
        self.assertEqual(pdb_obj_data['resolution'], 2.8)
        self.assertEqual(pdb_obj_data['structure_method'], 'X-RAY DIFFRACTION')
        self.assertEqual(pdb_obj_data['compound'], {})
        self.assertEqual(pdb_obj_data['source'], {})
        self.assertCountEqual(pdb_obj_data['proteins'][0].keys(),
                              ['id', 'sequence', 'md5', 'model_id', 'chain_id'])
        self.assertIn('mmcif_handle', pdb_obj_data.keys())
        self.assertIn('pdb_handle', pdb_obj_data.keys())
        self.assertIn('xml_handle', pdb_obj_data.keys())
        self.assertIn('rcsb_id', pdb_obj_data.keys())
        self.assertIn('release_date', pdb_obj_data.keys())

    @unittest.skip('test_structure_to_mmcif_file')
    def test_structure_to_mmcif_file(self):
        ret = self.serviceImpl.structure_to_pdb_file(self.ctx, {'input_ref': self.pdb_mmCif_ref,
                                                                'destination_dir': self.scratch})
        self.assertEqual(ret[0]['file_path'], os.path.join(self.scratch, '1fat.cif'))

    @unittest.skip('test_export_mmcif_structure')
    def test_export_mmcif_structure(self):
        ret = self.serviceImpl.export_pdb(self.ctx, {'input_ref': self.pdb_mmCif_ref})
        self.assertCountEqual(ret[0].keys(), ['shock_id'])

    # Testing PDBUtils module functions
    @unittest.skip('test_validate_file')
    def test_validate_file(self):
        not_exist_file = 'abc.csv'
        with self.assertRaisesRegex(
                ValueError,
                f'"{not_exist_file}" does not exist!'):
            self.pdb_util._validate_file(not_exist_file)

    @unittest.skip('test_read_file_by_type_csv')
    def test_read_file_by_type_csv(self):
        required_cols = ['Narrative ID', 'Object name (Genome AMA feature set)', 'Feature ID',
                         'PDB molecule', 'PDB filename']

        metafile = 'pdb_metafile_sample1a.csv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        ret_df = self.pdb_util._read_file_by_type(meta_file_path)
        self.assertCountEqual(ret_df.columns, required_cols)
        for col in required_cols:
            self.assertEqual(len(ret_df[col]), 7)
            col_list = ret_df[col].values.tolist()
            self.assertEqual(len(col_list), 7)

        required_cols.append('Is model')
        metafile = 'pdb_metafile_sample1b.csv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        ret_df = self.pdb_util._read_file_by_type(meta_file_path)
        self.assertCountEqual(ret_df.columns, required_cols)
        for col in required_cols:
            self.assertEqual(len(ret_df[col]), 4)
            col_list = ret_df[col].values.tolist()
            self.assertEqual(len(col_list), 4)

        metafile = 'pdb_metafile_sample2.csv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        ret_df = self.pdb_util._read_file_by_type(meta_file_path)
        self.assertCountEqual(ret_df.columns, required_cols)
        for col in required_cols:
            self.assertEqual(len(ret_df[col]), 4)
            col_list = ret_df[col].values.tolist()
            self.assertEqual(len(col_list), 4)

    @unittest.skip('test_read_file_by_type_tsv')
    def test_read_file_by_type_tsv(self):
        required_cols = ['Narrative ID', 'Object name (Genome AMA feature set)', 'Feature ID',
                         'PDB molecule', 'PDB filename']

        metafile = 'pdb_metafile_sample1a.tsv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        ret_df = self.pdb_util._read_file_by_type(meta_file_path)
        self.assertCountEqual(ret_df.columns, required_cols)
        for col in required_cols:
            self.assertEqual(len(ret_df[col]), 7)
            col_list = ret_df[col].values.tolist()
            self.assertEqual(len(col_list), 7)

        required_cols.append('Is model')
        metafile = 'pdb_metafile_sample1c.tsv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        ret_df = self.pdb_util._read_file_by_type(meta_file_path)
        self.assertCountEqual(ret_df.columns, required_cols)
        for col in required_cols:
            self.assertEqual(len(ret_df[col]), 4)
            col_list = ret_df[col].values.tolist()
            self.assertEqual(len(col_list), 4)

        metafile = 'pdb_metafile_sample2.tsv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        ret_df = self.pdb_util._read_file_by_type(meta_file_path)
        self.assertCountEqual(ret_df.columns, required_cols)
        for col in required_cols:
            self.assertEqual(len(ret_df[col]), 4)
            col_list = ret_df[col].values.tolist()
            self.assertEqual(len(col_list), 4)

    @unittest.skip('test_read_file_by_type_xlsx')
    def test_read_file_by_type_xlsx(self):
        required_cols = ['Narrative ID', 'Object name (Genome AMA feature set)', 'Feature ID',
                         'PDB molecule', 'PDB filename']

        metafile = 'pdb_metafile_sample1a.xlsx'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        ret_df = self.pdb_util._read_file_by_type(meta_file_path)
        self.assertCountEqual(ret_df.columns, required_cols)
        for col in required_cols:
            self.assertEqual(len(ret_df[col]), 7)
            col_list = ret_df[col].values.tolist()
            self.assertEqual(len(col_list), 7)

        required_cols.append('Is model')
        metafile = 'pdb_metafile_sample1d.xlsx'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        ret_df = self.pdb_util._read_file_by_type(meta_file_path)
        self.assertCountEqual(ret_df.columns, required_cols)
        for col in required_cols:
            self.assertEqual(len(ret_df[col]), 4)
            col_list = ret_df[col].values.tolist()
            self.assertEqual(len(col_list), 4)

        metafile = 'pdb_metafile_sample2.xlsx'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)

        ret_df = self.pdb_util._read_file_by_type(meta_file_path)
        self.assertCountEqual(ret_df.columns, required_cols)
        for col in required_cols:
            self.assertEqual(len(ret_df[col]), 4)
            col_list = ret_df[col].values.tolist()
            self.assertEqual(len(col_list), 4)

    @unittest.skip('test_incomplete_parse_metadata_csv_files')
    def test_parse_incomplete_metadata_csv_files(self):
        metafile = 'pdb_metafile_sample1a.csv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Trying to parse metadata file '{metafile}'")
        with self.assertRaisesRegex(
                ValueError,
                "Required column 'Is model' is missing!"):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        metafile = 'pdb_metafile_sample1b.csv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Trying to parse metadata file '{metafile}'")
        with self.assertRaisesRegex(
                ValueError,
                "Please fill all the rows in column 'PDB molecule'!"):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        metafile = 'pdb_metafile_sample1c.csv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Trying to parse metadata file '{metafile}'")
        with self.assertRaisesRegex(
                ValueError,
                "Please fill all the rows in column 'Feature ID'!"):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        metafile = 'pdb_metafile_sample1d.csv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Trying to parse metadata file '{metafile}'")
        with self.assertRaisesRegex(
                ValueError,
                "Please fill all the rows in column 'Object name'!"):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

    @unittest.skip('test_incomplete_parse_metadata_tsv_files')
    def test_parse_incomplete_metadata_tsv_files(self):
        metafile = 'pdb_metafile_sample1a.tsv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Trying to parse metadata file '{metafile}'")
        with self.assertRaisesRegex(
                ValueError,
                "Required column 'Is model' is missing!"):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        metafile = 'pdb_metafile_sample1b.tsv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Trying to parse metadata file '{metafile}'")
        with self.assertRaisesRegex(
                ValueError,
                "Please fill all the rows in column 'PDB molecule'!"):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        metafile = 'pdb_metafile_sample1c.tsv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Trying to parse metadata file '{metafile}'")
        with self.assertRaisesRegex(
                ValueError,
                "Please fill all the rows in column 'Feature ID'!"):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        metafile = 'pdb_metafile_sample1d.tsv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Trying to parse metadata file '{metafile}'")
        with self.assertRaisesRegex(
                ValueError,
                "Please fill all the rows in column 'Object name'!"):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

    @unittest.skip('test_incomplete_parse_metadata_xlsx_files')
    def test_parse_incomplete_metadata_xlsx_files(self):
        metafile = 'pdb_metafile_sample1a.xlsx'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Trying to parse metadata file '{metafile}'")
        with self.assertRaisesRegex(
                ValueError,
                "Required column 'Is model' is missing!"):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        metafile = 'pdb_metafile_sample1b.xlsx'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Trying to parse metadata file '{metafile}'")
        with self.assertRaisesRegex(
                ValueError,
                "Please fill all the rows in column 'PDB molecule'!"):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        metafile = 'pdb_metafile_sample1c.xlsx'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Trying to parse metadata file '{metafile}'")
        with self.assertRaisesRegex(
                ValueError,
                "Please fill all the rows in column 'Feature ID'!"):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        metafile = 'pdb_metafile_sample1d.xlsx'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Trying to parse metadata file '{metafile}'")
        with self.assertRaisesRegex(
                ValueError,
                "Please fill all the rows in column 'Object name'!"):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

    @unittest.skip('test_parse_complete_metadata_csv_file')
    def test_parse_complete_metadata_csv_file(self):
        metafile = 'pdb_metafile_sample2.csv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        (pdb_data, genome_objs, feature_ids,
            pdb_mols) = self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        self.assertEqual(len(pdb_data), 4)
        self.assertEqual(len(genome_objs), 4)
        self.assertEqual(len(feature_ids), 4)
        self.assertEqual(len(pdb_mols), 4)
        self.assertCountEqual(
            pdb_data[0].keys(),
            ['file_path', 'is_model', 'genome_object', 'feature_id', 'pdb_molecule'])

    @unittest.skip('test_parse_complete_metadata_tsv_file')
    def test_parse_complete_metadata_tsv_file(self):
        metafile = 'pdb_metafile_sample2.tsv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        (pdb_data, genome_objs, feature_ids,
            pdb_mols) = self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        self.assertEqual(len(pdb_data), 4)
        self.assertEqual(len(genome_objs), 4)
        self.assertEqual(len(feature_ids), 4)
        self.assertEqual(len(pdb_mols), 4)
        self.assertCountEqual(
            pdb_data[0].keys(),
            ['file_path', 'is_model', 'genome_object', 'feature_id', 'pdb_molecule'])

    @unittest.skip('test_parse_complete_metadata_xlsx_file')
    def test_parse_complete_metadata_xlsx_file(self):
        metafile = 'pdb_metafile_sample2.xlsx'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        (pdb_data, genome_objs, feature_ids,
            pdb_mols) = self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        self.assertEqual(len(pdb_data), 4)
        self.assertEqual(len(genome_objs), 4)
        self.assertEqual(len(feature_ids), 4)
        self.assertEqual(len(pdb_mols), 4)
        self.assertCountEqual(
            pdb_data[0].keys(),
            ['file_path', 'is_model', 'genome_object', 'feature_id', 'pdb_molecule'])

    @unittest.skip('test_compute_sequence_identity')
    def test_compute_sequence_identity(self):
        seq1 = 'TGTGACTA'
        seq2 = 'CATGGTCA'
        iden, exact_mat = self.pdb_util._compute_sequence_identity(seq1, seq2)
        self.assertEqual(iden, 0.0)
        self.assertTrue(not exact_mat)

        seq1 = 'SNDIYFNFQRFNETNLILQRDASVSSSGQLRLTNLN'
        seq2 = 'SNDIYFNFQRFNETNLILQRDASVSSSGQLRLTNLN'
        iden, exact_mat = self.pdb_util._compute_sequence_identity(seq1, seq2)
        self.assertEqual(iden, 1.0)
        self.assertTrue(exact_mat)

        seq1 = 'SNDIYFNFQRFNETNLILQRDASVSSSGQLRLTNL'
        seq2 = 'SNDIYFNFQRFNETNLILQRDASVSSSGQLRLTNLN'
        iden, exact_mat = self.pdb_util._compute_sequence_identity(seq1, seq2)
        self.assertEqual(iden, 1.0)
        self.assertTrue(exact_mat)

        seq1 = 'SNDIYFNFQRFNETNLILQRDASVSSSGQLRLTNLN'
        seq2 = 'NDIYFNFQRFNETNLILQRDASVSSSGQLRLTNLN'
        iden, exact_mat = self.pdb_util._compute_sequence_identity(seq1, seq2)
        self.assertEqual(iden, 1.0)
        self.assertTrue(exact_mat)

        seq1 = 'FQTWEEFSRAAEKLYLADPMKVRVVLKYRHVDGNLCIKVTDDLVCLVYRTDQAQDVKKIEKF'
        seq2 = 'FQTWEEFSRAEKLYLADPMKVRVVLRYRHVDGNLCIKVTDDLICLVYRTDQAQDVKKIEKF'
        iden, exact_mat = self.pdb_util._compute_sequence_identity(seq1, seq2)
        expected_iden = 0.967213
        self.assertEqual(iden, expected_iden)
        self.assertTrue(not exact_mat)

    @unittest.skip('test_match_features')  # Note the genome is from an Appdev narrative 57196
    def test_match_features(self):
        fileName = '6ift.pdb'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)

        parser = PDB.PDBParser(PERMISSIVE=1)
        structure = parser.get_structure("test", pdb_file_path)
        model = structure[0]
        protein_data1 = self.pdb_util._get_proteins_by_structure(structure, model.get_id(),
                                                                 pdb_file_path)
        self.assertCountEqual(protein_data1[0].keys(),
                              ['id', 'model_id', 'chain_id', 'sequence', 'md5'])
        # A genome that has matching features to the pdb protein
        narr_id = 57196
        genome_name = 'Synthetic_bacterium_JCVI_Syn3.0_genome'
        feat_id = 'JCVISYN3_0004'
        protein_data1 = self.pdb_util._match_features(narr_id, genome_name, feat_id, protein_data1)
        self.assertCountEqual(protein_data1[0].keys(),
                              ['id', 'sequence', 'md5', 'model_id', 'chain_id',
                               'seq_identity', 'exact_match'])

        # A genome that has NO matching features to the pdb protein
        narr_id = 42297
        genome_name = 'OntSer_GCF_001699635_Feb20b'
        feat_id = 'GCF_001699635.1.CDS.1'
        protein_data2 = self.pdb_util._get_proteins_by_structure(structure, model.get_id(),
                                                                 pdb_file_path)
        protein_data2 = self.pdb_util._match_features(narr_id, genome_name, feat_id, protein_data2)
        self.assertCountEqual(protein_data2[0].keys(),
                              ['id', 'sequence', 'md5', 'model_id', 'chain_id'])

    @unittest.skip('test_model_file_to_data')  # Note the genome is from an Appdev narrative 42297
    def test_model_file_to_data(self):
        fileName = '1fat.pdb'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)
        narr_id = 42297
        genome_name = 'OntSer_GCF_001699635_Feb20b'
        feat_id = 'GCF_001699635.1.CDS.1'

        params = {
                'narrative_id': narr_id,
                'genome_name': genome_name,
                'feature_id': feat_id
        }
        (data1, pp_no1) = self.pdb_util._model_file_to_data(pdb_file_path, params)

        self.assertEqual(pp_no1, 7)
        self.assertEqual(data1['name'], 'phytohemagglutinin-l')
        self.assertEqual(data1['num_chains'], 4)
        self.assertEqual(data1['num_residues'], 928)
        self.assertEqual(data1['num_atoms'], 7248)
        self.assertEqual(len(data1['proteins']), data1['num_chains'])
        self.assertCountEqual(data1['proteins'][0].keys(),
                              ['id', 'sequence', 'md5', 'model_id', 'chain_id'])
        self.assertCountEqual(data1['compound'].keys(), ['misc', 'molecule', 'chain', 'synonym'])
        self.assertCountEqual(data1['source'].keys(),
                              ['misc', 'organism_scientific',
                               'organism_taxid', 'organ', 'other_details'])

        fileName = '5o5y.pdb'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)
        (data2, pp_no2) = self.pdb_util._model_file_to_data(pdb_file_path, params)

        self.assertEqual(pp_no2, 2)
        self.assertEqual(
            data2['name'],
            'crystal structure of thermococcus litoralis adp-dependent glucokinase (gk)')
        self.assertEqual(data2['num_chains'], 2)
        self.assertEqual(data2['num_residues'], 903)
        self.assertEqual(data2['num_atoms'], 8013)
        self.assertEqual(len(data2['proteins']), data2['num_chains'])
        self.assertCountEqual(data2['proteins'][0].keys(),
                              ['id', 'sequence', 'md5', 'model_id', 'chain_id'])
        self.assertCountEqual(data2['compound'].keys(),
                              ['misc', 'molecule', 'chain', 'synonym', 'ec_number', 'ec',
                               'engineered'])
        self.assertCountEqual(data2['source'].keys(),
                              ['misc', 'organism_scientific', 'expression_system',
                               'gene', 'expression_system_taxid', 'expression_system_vector_type',
                               'organism_taxid', 'expression_system_vector'])

        fileName = '6ift.pdb'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)
        narr_id = 57196
        genome_name = 'Synthetic_bacterium_JCVI_Syn3.0_genome'
        feat_id = 'JCVISYN3_0004'
        params = {
                'narrative_id': narr_id,
                'genome_name': genome_name,
                'feature_id': feat_id
        }
        (data3, pp_no3) = self.pdb_util._model_file_to_data(pdb_file_path, params)

        self.assertEqual(pp_no3, 1)
        self.assertEqual(
            data3['name'],
            'ksga from bacillus subtilis in complex with sam')
        self.assertEqual(data3['num_chains'], 1)
        self.assertEqual(data3['num_residues'], 292)
        self.assertEqual(data3['num_atoms'], 2508)
        self.assertEqual(len(data3['proteins']), data3['num_chains'])
        self.assertCountEqual(data3['proteins'][0].keys(),
                              ['id', 'sequence', 'md5', 'model_id', 'chain_id',
                               'seq_identity', 'exact_match'])
        self.assertCountEqual(data3['compound'].keys(),
                              ['misc', 'molecule', 'chain', 'synonym', 'ec_number', 'ec',
                               'engineered'])
        self.assertEqual(data3['compound']['ec_number'], '2.1.1.182')
        self.assertCountEqual(data3['source'].keys(),
                              ['misc', 'organism_scientific', 'expression_system', 'gene',
                               'expression_system_taxid', 'strain', 'organism_taxid'])

    @unittest.skip('test_exp_file_to_data')
    def test_exp_file_to_data(self):
        fileName = '1fat.cif'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)
        narr_id = 42297
        genome_name = 'OntSer_GCF_001699635_Feb20b'
        feat_id = 'GCF_001699635.1.CDS.1'

        params = {
                'narrative_id': narr_id,
                'genome_name': genome_name,
                'feature_id': feat_id
        }
        (data, pp_no) = self.pdb_util._exp_file_to_data(pdb_file_path, params)

        self.assertEqual(pp_no, 7)
        self.assertEqual(data['name'], 'PHYTOHEMAGGLUTININ-L')
        self.assertEqual(data['head'], 'LECTIN')
        self.assertEqual(data['rcsb_id'], '')
        self.assertEqual(data['deposition_date'], '1996-06-12')
        self.assertEqual(data['release_date'], '')
        self.assertEqual(data['structure_method'], 'X-RAY DIFFRACTION')
        self.assertEqual(data['resolution'], 2.8)
        self.assertEqual(data['compound'], {})
        self.assertEqual(data['source'], {})
        self.assertTrue('pdb_handle' in data.keys())
        self.assertEqual(len(data['proteins']), data['num_chains'])
        self.assertCountEqual(data['proteins'][0].keys(),
                              ['id', 'sequence', 'md5', 'model_id', 'chain_id'])
        self.assertEqual(data['num_chains'], 4)
        self.assertEqual(data['num_residues'], 928)
        self.assertEqual(data['num_atoms'], 7248)
        self.assertEqual(len(data['proteins']), data['num_chains'])
        self.assertCountEqual(data['proteins'][0].keys(),
                              ['id', 'sequence', 'md5', 'model_id', 'chain_id'])
        self.assertCountEqual(data['compound'].keys(), [])
        self.assertCountEqual(data['source'].keys(), [])

    @unittest.skip('test_import_model_pdb_file_nopatch')
    def test_import_model_pdb_file_nopatch(self):
        fileName = '1fat.pdb'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)
        params = {
            'input_shock_id': '',
            'input_file_path': fileName,
            'input_staging_file_path': '',
            'structure_name': 'test_pdb_structure_1fat',
            'description': 'for test PDBUtils.import_model_pdb_file',
            'workspace_name': self.wsName
        }
        ret = self.pdb_util.import_model_pdb_file(params, False)
        self.assertIn('structure_obj_ref', ret)

        fileName = '6ift.pdb'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)
        narr_id = 42297
        genome_name = 'OntSer_GCF_001699635_Feb20b'
        feat_id = 'GCF_001699635.1.CDS.1'
        params = {
            'input_shock_id': '',
            'input_file_path': fileName,
            'input_staging_file_path': '',
            'structure_name': 'test_pdb_structure_6ift',
            'narrative_id': narr_id,
            'genome_name': genome_name,
            'feature_id': feat_id,
            'description': 'for test PDBUtils.import_model_pdb_file',
            'workspace_name': self.wsName
        }
        ret = self.pdb_util.import_model_pdb_file(params, False)
        self.assertIn('structure_obj_ref', ret)

    @unittest.skip('test_import_model_pdb_file_patched')
    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    def test_import_model_pdb_file_patched(self, download_staging_file):
        fileName = '1fat.pdb'
        params = {
            'input_shock_id': '',
            'input_file_path': '',
            'input_staging_file_path': fileName,
            'structure_name': 'test_pdb_structure_1fat',
            'description': 'for test PDBUtils.import_model_pdb_file',
            'workspace_name': self.wsName
        }
        ret = self.pdb_util.import_model_pdb_file(params, False)
        self.assertIn('structure_obj_ref', ret)

        fileName = '6ift.pdb'
        narr_id = 42297
        genome_name = 'OntSer_GCF_001699635_Feb20b'
        feat_id = 'GCF_001699635.1.CDS.1'
        params = {
            'input_shock_id': '',
            'input_file_path': '',
            'input_staging_file_path': fileName,
            'structure_name': 'test_pdb_structure_6ift',
            'narrative_id': narr_id,
            'genome_name': genome_name,
            'feature_id': feat_id,
            'description': 'for test PDBUtils.import_model_pdb_file',
            'workspace_name': self.wsName
        }
        ret = self.pdb_util.import_model_pdb_file(params, False)
        self.assertIn('structure_obj_ref', ret)

    @unittest.skip('test_import_experiment_pdb_file_nopatch')
    def test_import_experiment_pdb_file_nopatch(self):
        fileName = '1fat.cif'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)
        params = {
            'input_shock_id': '',
            'input_file_path': pdb_file_path,
            'input_staging_file_path': '',
            'structure_name': 'test_pdb_structure_name',
            'description': 'for test PDBUtils.import_exp_pdb_file',
            'workspace_name': self.wsName
        }
        ret = self.pdb_util.import_experiment_pdb_file(params, False)
        self.assertIn('structure_obj_ref', ret)

    @unittest.skip('test_import_experiment_pdb_file_pathched')
    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    def test_import_experiment_pdb_file_patched(self, download_staging_file):
        fileName = '1fat.cif'
        params = {
            'input_shock_id': '',
            'input_file_path': '',
            'input_staging_file_path': fileName,
            'structure_name': 'test_pdb_structure_name',
            'description': 'for test PDBUtils.import_exp_pdb_file',
            'workspace_name': self.wsName
        }
        ret = self.pdb_util.import_experiment_pdb_file(params, False)
        self.assertIn('structure_obj_ref', ret)

    #@unittest.skip('test_batch_import_pdbs_from_metafile')
    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    def test_batch_import_pdbs_from_metafile(self, download_staging_file):
        metafile = 'pdb_metafile_sample1.csv'
        metafile = os.path.join('/kb/module/test/data', metafile)

        params = {
            'metadata_staging_file_path': metafile,
            'structures_name': 'batch1_test_structures',
            'workspace_name': self.wsName
        }
        ret = self.serviceImpl.batch_import_pdbs_from_metafile(self.ctx, params)
        self.assertCountEqual(ret[0].keys(), ["structures_ref"])  #, "report_ref", "report_name"])

        metafile = 'pdb_metafile_sample2.csv'
        metafile = os.path.join('/kb/module/test/data', metafile)

        params = {
            'metadata_staging_file_path': metafile,
            'structures_name': 'batch2_test_structures',
            'workspace_name': self.wsName
        }
        ret = self.serviceImpl.batch_import_pdbs_from_metafile(self.ctx, params)
        self.assertCountEqual(ret[0].keys(), ["structures_ref"])  #, "report_ref", "report_name"])
