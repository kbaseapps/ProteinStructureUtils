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

from ProteinStructureUtils.Utils.PDBUtil import PDBUtil

from Bio import PDB


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
            'is_model': 1,
        }
        info = cls.dfu.save_objects({
            'id': cls.ws_id,
            'objects': [
                {'type': 'KBaseStructure.ProteinStructure',
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
            'is_model': 0,
            'pdb_handle': shock_id,
            'mmcif_handle': shock_id,
            'xml_handle': shock_id
        }
        info = cls.dfu.save_objects({
            'id': cls.ws_id,
            'objects': [
                {'type': 'KBaseStructure.ProteinStructure',
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

    # Testing PDBUtil module functions
    #@unittest.skip('test_validate_file')
    def test_validate_file(self):
        not_exist_file = 'abc.csv'
        with self.assertRaisesRegex(
                ValueError,
                f'"{not_exist_file}" does not exist!'):
            self.pdb_util._validate_file(not_exist_file)

    #@unittest.skip('test_read_file_by_type_csv')
    def test_read_file_by_type_csv(self):
        required_cols = ['Narrative ID', 'Object name (Genome AMA feature set)',
                         'Feature ID', 'PDB or RCSB ID', 'File Extension', 'Is model', 'From RCSB']

        metafile = 'pdb_metafile_sample1a.csv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        ret_df = self.pdb_util._read_file_by_type(meta_file_path)
        self.assertCountEqual(ret_df.columns, required_cols)
        for col in required_cols:
            self.assertEqual(len(ret_df[col]), 7)
            col_list = ret_df[col].values.tolist()
            self.assertEqual(len(col_list), 7)

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
            self.assertEqual(len(ret_df[col]), 8)
            col_list = ret_df[col].values.tolist()
            self.assertEqual(len(col_list), 8)

    #@unittest.skip('test_read_file_by_type_tsv')
    def test_read_file_by_type_tsv(self):
        required_cols = ['Narrative ID', 'Object name (Genome AMA feature set)',
                         'Feature ID', 'PDB or RCSB ID', 'File Extension', 'Is model', 'From RCSB']

        metafile = 'pdb_metafile_sample1a.tsv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        ret_df = self.pdb_util._read_file_by_type(meta_file_path)
        self.assertCountEqual(ret_df.columns, required_cols)
        for col in required_cols:
            self.assertEqual(len(ret_df[col]), 7)
            col_list = ret_df[col].values.tolist()
            self.assertEqual(len(col_list), 7)

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
            self.assertEqual(len(ret_df[col]), 7)
            col_list = ret_df[col].values.tolist()
            self.assertEqual(len(col_list), 7)

    #@unittest.skip('test_read_file_by_type_xlsx')
    def test_read_file_by_type_xlsx(self):
        required_cols = ['Narrative ID', 'Object name (Genome AMA feature set)',
                         'Feature ID', 'PDB or RCSB ID', 'File Extension', 'Is model', 'From RCSB']

        metafile = 'pdb_metafile_sample1a.xlsx'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        ret_df = self.pdb_util._read_file_by_type(meta_file_path)
        self.assertCountEqual(ret_df.columns, required_cols)
        for col in required_cols:
            self.assertEqual(len(ret_df[col]), 7)
            col_list = ret_df[col].values.tolist()
            self.assertEqual(len(col_list), 7)

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
            self.assertEqual(len(ret_df[col]), 7)
            col_list = ret_df[col].values.tolist()
            self.assertEqual(len(col_list), 7)

    #@unittest.skip('test_incomplete_parse_metadata_csv_files')
    def test_parse_incomplete_metadata_csv_files(self):
        metafile = 'pdb_metafile_sample1a.csv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Parsing metadata file '{metafile}'......")
        with self.assertRaisesRegex(
                ValueError,
                'Please fill all the rows in column: Feature ID!'):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        metafile = 'pdb_metafile_sample1b.csv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Parsing metadata file '{metafile}'......")
        with self.assertRaisesRegex(
                ValueError,
                "Please fill all the rows in column: PDB or RCSB ID!"):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        metafile = 'pdb_metafile_sample1c.csv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Parsing metadata file '{metafile}'......")
        with self.assertRaisesRegex(
                ValueError,
                "Please fill all the rows in column: Feature ID!"):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        metafile = 'pdb_metafile_sample1d.csv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Parsing metadata file '{metafile}'......")
        with self.assertRaisesRegex(
                ValueError,
                'Please fill all the rows in column: Object name!'):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

    #@unittest.skip('test_incomplete_parse_metadata_tsv_files')
    def test_parse_incomplete_metadata_tsv_files(self):
        metafile = 'pdb_metafile_sample1a.tsv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Parsing metadata file '{metafile}'......")
        with self.assertRaisesRegex(
                ValueError,
                'Please fill all the rows in column: Feature ID!'):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        metafile = 'pdb_metafile_sample1b.tsv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Parsing metadata file '{metafile}'......")
        with self.assertRaisesRegex(
                ValueError,
                "Please fill all the rows in column: PDB or RCSB ID!"):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        metafile = 'pdb_metafile_sample1c.tsv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Parsing metadata file '{metafile}'......")
        with self.assertRaisesRegex(
                ValueError,
                "Please fill all the rows in column: Feature ID!"):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        metafile = 'pdb_metafile_sample1d.tsv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Parsing metadata file '{metafile}'......")
        with self.assertRaisesRegex(
                ValueError,
                'Please fill all the rows in column: Object name!'):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

    #@unittest.skip('test_incomplete_parse_metadata_xlsx_files')
    def test_parse_incomplete_metadata_xlsx_files(self):
        metafile = 'pdb_metafile_sample1a.xlsx'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Parsing metadata file '{metafile}'......")
        with self.assertRaisesRegex(
                ValueError,
                'Please fill all the rows in column: Feature ID!'):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        metafile = 'pdb_metafile_sample1b.xlsx'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Parsing metadata file '{metafile}'......")
        with self.assertRaisesRegex(
                ValueError,
                'Required column "PDB or RCSB ID" is missing!'):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        metafile = 'pdb_metafile_sample1c.xlsx'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Parsing metadata file '{metafile}'......")
        with self.assertRaisesRegex(
                ValueError,
                'Required column "PDB or RCSB ID" is missing!'):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        metafile = 'pdb_metafile_sample1d.xlsx'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        print(f"Parsing metadata file '{metafile}'......")
        with self.assertRaisesRegex(
                ValueError,
                'Please fill all the rows in column: Object name!'):
            self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

    #@unittest.skip('test_parse_complete_metadata_csv_file')
    def test_parse_complete_metadata_csv_file(self):
        metafile = 'pdb_metafile_sample2.csv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        (pdb_data, narr_ids, genome_objs,
            feature_ids) = self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        self.assertEqual(len(pdb_data), 8)
        self.assertEqual(len(genome_objs), 8)
        self.assertEqual(len(feature_ids), 8)
        self.assertCountEqual(
            pdb_data[0].keys(),
            ['file_path', 'is_model', 'genome_name', 'structure_name', 'narrative_id',
             'feature_id', 'file_extension', 'from_rcsb'])

    #@unittest.skip('test_parse_complete_metadata_tsv_file')
    def test_parse_complete_metadata_tsv_file(self):
        metafile = 'pdb_metafile_sample2.tsv'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        (pdb_data, narr_ids, genome_objs,
            feature_ids) = self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        self.assertEqual(len(pdb_data), 7)
        self.assertEqual(len(narr_ids), 7)
        self.assertEqual(len(genome_objs), 7)
        self.assertEqual(len(feature_ids), 7)
        self.assertCountEqual(
            pdb_data[0].keys(),
            ['file_path', 'is_model', 'genome_name', 'structure_name', 'narrative_id',
             'feature_id', 'file_extension', 'from_rcsb'])

    #@unittest.skip('test_parse_complete_metadata_xlsx_file')
    def test_parse_complete_metadata_xlsx_file(self):
        metafile = 'pdb_metafile_sample2.xlsx'
        meta_file_path = os.path.join(self.scratch, metafile)
        shutil.copy(os.path.join('data', metafile), meta_file_path)
        (pdb_data, narr_ids, genome_objs,
            feature_ids) = self.pdb_util._parse_metadata_file(meta_file_path, self.ws_id)

        self.assertEqual(len(pdb_data), 7)
        self.assertEqual(len(narr_ids), 7)
        self.assertEqual(len(genome_objs), 7)
        self.assertEqual(len(feature_ids), 7)
        self.assertCountEqual(
            pdb_data[0].keys(),
            ['file_path', 'is_model', 'genome_name', 'structure_name', 'narrative_id',
             'feature_id', 'file_extension', 'from_rcsb'])

    #@unittest.skip('test_structure_to_pdb_file')
    def test_structure_to_pdb_file(self):
        ret = self.pdb_util._structure_to_pdb_file({'input_ref': self.pdb_ref,
                                                    'destination_dir': self.scratch})
        self.assertEqual(ret['file_path'], os.path.join(self.scratch, '1nqg.pdb'))

    #@unittest.skip('test_structure_to_mmcif_file')
    def test_structure_to_mmcif_file(self):
        ret = self.pdb_util._structure_to_pdb_file({'input_ref': self.pdb_mmCif_ref,
                                                    'destination_dir': self.scratch})
        self.assertEqual(ret['file_path'], os.path.join(self.scratch, '1fat.cif'))

    #@unittest.skip('test_export_pdb_structure')
    def test_export_pdb_structure(self):
        ret = self.pdb_util._export_pdb({'input_ref': self.pdb_ref})
        self.assertCountEqual(ret.keys(), ['shock_id'])

    #@unittest.skip('test_export_mmcif_structure')
    def test_export_mmcif_structure(self):
        ret = self.pdb_util._export_pdb({'input_ref': self.pdb_mmCif_ref})
        self.assertCountEqual(ret.keys(), ['shock_id'])

    #@unittest.skip('test_compute_sequence_identity')
    def test_compute_sequence_identity(self):
        seq1 = 'TGTGACTA'
        seq2 = 'CATGGTCA'
        idens, exact_mat, evals = self.pdb_util._compute_sequence_identity(seq1, seq2)
        self.assertEqual(idens, [])
        self.assertTrue(not exact_mat)

        seq1 = 'SNDIYFNFQRFNETNLILQRDASVSSSGQLRLTNLN'
        seq2 = 'SNDIYFNFQRFNETNLILQRDASVSSSGQLRLTNLN'
        idens, exact_mat, evals = self.pdb_util._compute_sequence_identity(seq1, seq2)
        self.assertEqual(idens, [1.0])
        self.assertTrue(exact_mat)

        seq1 = 'SNDIYFNFQRFNETNLILQRDASVSSSGQLRLTNL'
        seq2 = 'SNDIYFNFQRFNETNLILQRDASVSSSGQLRLTNLN'
        idens, exact_mat, evals = self.pdb_util._compute_sequence_identity(seq1, seq2)
        self.assertEqual(idens, [1.0])
        self.assertTrue(exact_mat)

        seq1 = 'SNDIYFNFQRFNETNLILQRDASVSSSGQLRLTNLN'
        seq2 = 'NDIYFNFQRFNETNLILQRDASVSSSGQLRLTNLN'
        idens, exact_mat, evals = self.pdb_util._compute_sequence_identity(seq1, seq2)
        self.assertEqual(idens, [1.0])
        self.assertTrue(exact_mat)

        seq1 = 'FQTWEEFSRAAEKLYLADPMKVRVVLKYRHVDGNLCIKVTDDLVCLVYRTDQAQDVKKIEKF'
        seq2 = 'FQTWEEFSRAEKLYLADPMKVRVVLRYRHVDGNLCIKVTDDLICLVYRTDQAQDVKKIEKF'
        idens, exact_mat, evals = self.pdb_util._compute_sequence_identity(seq1, seq2)
        expected_iden = 0.967213
        self.assertEqual(idens, [expected_iden])
        self.assertTrue(not exact_mat)

        # Change the Evalue and identity threshold
        idens, exact_mat, evals = self.pdb_util._compute_sequence_identity(
                                                    seq1, seq2, evalue_threshold=1e-10)
        expected_iden = 0.967213
        self.assertEqual(idens, [expected_iden])
        self.assertTrue(not exact_mat)

    #@unittest.skip('test_match_features')  # Note the genome is from an Appdev narrative 57196
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
        # A genome that has a matching feature to the pdb protein
        pdb_info = {
            'narrative_id': 57196,
            'genome_name': 'Synthetic_bacterium_JCVI_Syn3_genome',
            'feature_id': 'JCVISYN3_0004'
        }
        params = {'pdb_info': pdb_info}
        protein_data1, params1 = self.pdb_util._match_features(params, protein_data1)
        self.assertCountEqual(protein_data1[0].keys(),
                              ['id', 'sequence', 'md5', 'model_id', 'chain_id',
                               'seq_identity', 'exact_match', 'genome_ref',
                               'feature_id', 'feature_type'])
        self.assertIn('genome_name', params1['pdb_info'])
        self.assertIn('genome_ref', params1['pdb_info'])
        self.assertIn('feature_id', params1['pdb_info'])
        self.assertIn('feature_type', params1['pdb_info'])
        self.assertIn('chain_ids', params1['pdb_info'])
        self.assertIn('model_ids', params1['pdb_info'])
        self.assertIn('sequence_identities', params1['pdb_info'])
        self.assertIn('exact_matches', params1['pdb_info'])

        # A genome that has NO matching feature to the pdb protein
        pdb_info = {
            'narrative_id': 42297,
            'genome_name': 'OntSer_GCF_001699635_Feb20b',
            'feature_id': 'GCF_001699635.1.CDS.1'
        }
        params = {'pdb_info': pdb_info}
        protein_data2 = self.pdb_util._get_proteins_by_structure(structure, model.get_id(),
                                                                 pdb_file_path)
        protein_data2, params2 = self.pdb_util._match_features(params, protein_data2)
        self.assertCountEqual(protein_data2[0].keys(),
                              ['id', 'sequence', 'md5', 'model_id', 'chain_id'])
        self.assertIn('genome_name', params2['pdb_info'])
        self.assertIn('genome_ref', params2['pdb_info'])
        self.assertIn('feature_id', params2['pdb_info'])
        self.assertIn('feature_type', params2['pdb_info'])
        self.assertNotIn('chain_ids', params2['pdb_info'])
        self.assertNotIn('model_ids', params2['pdb_info'])
        self.assertNotIn('sequence_identities', params2['pdb_info'])
        self.assertNotIn('exact_matches', params2['pdb_info'])

    #@unittest.skip('test_match_features_nonexist_genome')
    def test_match_features_nonexist_genome(self):
        fileName = 'MLuteus_AlphaFold_133.pdb'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)

        parser = PDB.PDBParser(PERMISSIVE=1)
        struct = parser.get_structure("ML_alphafold_133", pdb_file_path)
        model = struct[0]
        prot_data = self.pdb_util._get_proteins_by_structure(struct, model.get_id(),
                                                             pdb_file_path)
        self.assertCountEqual(prot_data[0].keys(),
                              ['id', 'model_id', 'chain_id', 'sequence', 'md5'])

        # Try to match with features in a non-existent genome
        pdb_info = {
            'narrative_id': 63679,
            'genome_name': 'MLuteus_ATCC_nonexist',
            'feature_id': 'MLuteus_masurca_RAST.CDS.133'
        }
        params = {'pdb_info': pdb_info}

        # No features matched, nothing is changed by the _match_features(..) call
        prot_data1, params1 = self.pdb_util._match_features(params, prot_data)
        self.assertEqual(params1, params)
        self.assertEqual(prot_data1, prot_data)

    #@unittest.skip('test_pdb_file_to_data')  # Note the genome is from an Appdev narrative 42297
    def test_pdb_file_to_data(self):
        fileName = '1fat.pdb'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)
        narr_id = 42297
        genome_name = 'OntSer_GCF_001699635_Feb20b'
        feat_id = 'GCF_001699635.1.CDS.1'

        pdb_info = {
                'narrative_id': narr_id,
                'genome_name': genome_name,
                'feature_id': feat_id
        }
        params = {'pdb_info': pdb_info}
        (data1, pp_no1, params1) = self.pdb_util._pdb_file_to_data(pdb_file_path, params)
        self.assertFalse(data1)

        fileName = '5o5y.pdb'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)
        (data2, pp_no2, params2) = self.pdb_util._pdb_file_to_data(pdb_file_path, params)
        self.assertFalse(data2)

        fileName = '6ift.pdb'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)
        narr_id = 57196
        genome_name = 'Synthetic_bacterium_JCVI_Syn3_genome'
        feat_id = 'JCVISYN3_0004'
        pdb_info = {
                'narrative_id': narr_id,
                'genome_name': genome_name,
                'feature_id': feat_id
        }
        params = {'pdb_info': pdb_info}
        (data3, pp_no3, params3) = self.pdb_util._pdb_file_to_data(pdb_file_path, params)

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
                               'seq_identity', 'exact_match', 'genome_ref',
                               'feature_id', 'feature_type'])
        self.assertCountEqual(data3['compound'].keys(),
                              ['misc', 'molecule', 'chain', 'synonym', 'ec_number', 'ec',
                               'engineered'])
        self.assertEqual(data3['compound']['ec_number'], '2.1.1.182')
        self.assertCountEqual(data3['source'].keys(),
                              ['misc', 'organism_scientific', 'expression_system', 'gene',
                               'expression_system_taxid', 'strain', 'organism_taxid'])

        self.assertIn('genome_name', params3['pdb_info'])
        self.assertIn('genome_ref', params3['pdb_info'])
        self.assertIn('feature_id', params3['pdb_info'])
        self.assertIn('feature_type', params3['pdb_info'])
        self.assertIn('chain_ids', params3['pdb_info'])
        self.assertIn('model_ids', params3['pdb_info'])
        self.assertIn('sequence_identities', params3['pdb_info'])
        self.assertIn('exact_matches', params3['pdb_info'])

    #@unittest.skip('test_pdb_file_to_data_erred')  # Note the genome is from an Appdev narrative 63679
    def test_pdb_file_to_data_erred(self):
        fileName = 'MLuteus_AlphaFold_133.pdb'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)

        genome_name = 'MLuteus_ATCC_nonexist'
        feat_id = 'MLuteus_masurca_RAST.CDS.133'

        pdb_info = {
                'narrative_id': 63679,
                'genome_name': genome_name,
                'feature_id': feat_id
        }
        params = {'pdb_info': pdb_info}
        (data1, pp_no1, params1) = self.pdb_util._pdb_file_to_data(pdb_file_path, params)
        self.assertFalse(data1)

    #@unittest.skip('test_mmcif_file_to_data')
    def test_mmcif_file_to_data(self):
        fileName = '1fat.cif'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)
        narr_id = 42297
        genome_name = 'OntSer_GCF_001699635_Feb20b'
        feat_id = 'GCF_001699635.1.CDS.1'

        pdb_info = {
                'narrative_id': narr_id,
                'genome_name': genome_name,
                'feature_id': feat_id
        }
        params = {'pdb_info': pdb_info}
        (data, pp_no, params1) = self.pdb_util._mmcif_file_to_data(pdb_file_path, params)
        self.assertFalse(data)
        self.assertEqual(pp_no, 7)
        print(params1)

    #@unittest.skip('test_import_pdb_file_nopatch')
    def test_import_pdb_file_nopatch(self):
        fileName = '1fat.pdb'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)
        params = {
            'input_shock_id': '',
            'input_file_path': pdb_file_path,
            'input_staging_file_path': '',
            'structure_name': 'test_pdb_structure_1fat',
            'description': 'for test PDBUtil.import_pdb_file_nopatch',
            'workspace_name': self.wsName
        }
        ret_data, ret_info = self.pdb_util.import_pdb_file(params)
        self.assertFalse(ret_data)
        self.assertFalse(ret_info)

        fileName = '6ift.pdb'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)
        narr_id = 57196
        genome_name = 'Synthetic_bacterium_JCVI_Syn3_genome'
        feat_id = 'JCVISYN3_0004'
        params = {
            'input_shock_id': '',
            'input_file_path': pdb_file_path,
            'input_staging_file_path': '',
            'structure_name': 'test_pdb_structure_6ift',
            'pdb_info': {'narrative_id': narr_id,
                         'genome_name': genome_name,
                         'feature_id': feat_id},
            'description': 'for test PDBUtil.import_pdb_file_nopatch',
            'workspace_name': self.wsName
        }
        ret_data, ret_info = self.pdb_util.import_pdb_file(params)
        self._check_import_model_pdb_6ift(ret_data, ret_info)
        self.assertEqual(ret_data['user_data'], params['description'])

    #@unittest.skip('test_import_pdb_file_patched')
    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    def test_import_pdb_file_patched(self, download_staging_file):
        fileName = '1fat.pdb'
        params = {
            'input_shock_id': '',
            'input_file_path': '',
            'input_staging_file_path': fileName,
            'structure_name': 'test_pdb_structure_1fat',
            'description': 'for test PDBUtil.import_pdb_file_patch',
            'workspace_name': self.wsName
        }
        ret_data, ret_info = self.pdb_util.import_pdb_file(params)
        self.assertFalse(ret_data)
        self.assertFalse(ret_info)

        fileName = '6ift.pdb'
        narr_id = 57196
        genome_name = 'Synthetic_bacterium_JCVI_Syn3_genome'
        feat_id = 'JCVISYN3_0004'
        params = {
            'input_shock_id': '',
            'input_file_path': '',
            'input_staging_file_path': fileName,
            'structure_name': 'test_pdb_structure_6ift',
            'pdb_info': {'narrative_id': narr_id,
                         'genome_name': genome_name,
                         'feature_id': feat_id},
            'description': 'for test PDBUtil.import_pdb_file_patch',
            'workspace_name': self.wsName
        }
        ret_data, ret_info = self.pdb_util.import_pdb_file(params)
        self._check_import_model_pdb_6ift(ret_data, ret_info)
        self.assertEqual(ret_data['user_data'], params['description'])

    def _check_import_model_pdb_6ift(self, pdb_data, pdb_info):
        self.assertEqual(pdb_data['name'], 'ksga from bacillus subtilis in complex with sam')
        self.assertEqual(pdb_data['num_chains'], 1)
        self.assertEqual(pdb_data['num_residues'], 292)
        self.assertEqual(pdb_data['num_atoms'], 2508)
        self.assertEqual(pdb_data['compound']['molecule'], 'ribosomal rna small subunit methyltransferase a')
        self.assertEqual(pdb_data['compound']['ec_number'], '2.1.1.182')
        self.assertEqual(pdb_data['source']['organism_scientific'], 'bacillus subtilis (strain 168)')
        self.assertEqual(pdb_data['source']['organism_taxid'], '224308')
        self.assertEqual(len(pdb_data['proteins']), 1)

        self.assertEqual(pdb_info['narrative_id'], 57196)
        self.assertEqual(pdb_info['genome_name'], 'Synthetic_bacterium_JCVI_Syn3_genome')
        self.assertEqual(pdb_info['feature_id'], 'JCVISYN3_0004')
        self.assertEqual(pdb_info['genome_ref'], '57196/6/1')
        self.assertEqual(pdb_info['sequence_identities'], '67.95%')
        self.assertEqual(pdb_info['chain_ids'], 'Model 1.Chain A')
        self.assertEqual(pdb_info['feature_type'], 'gene')
        self.assertEqual(pdb_info['model_ids'], '0')
        self.assertEqual(pdb_info['scratch_path'], '/kb/module/work/tmp/6ift.pdb')

    #@unittest.skip('test_import_mmcif_file_nopatch')
    def test_import_mmcif_file_nopatch(self):
        fileName = '1fat.cif'
        pdb_file_path = os.path.join(self.scratch, fileName)
        shutil.copy(os.path.join('data', fileName), pdb_file_path)
        params = {
            'input_shock_id': '',
            'input_file_path': pdb_file_path,
            'input_staging_file_path': '',
            'structure_name': 'test_pdb_structure_name',
            'description': 'for test PDBUtil.import_mmcif_file',
            'workspace_name': self.wsName
        }
        ret_data, ret_info = self.pdb_util.import_mmcif_file(params)
        self.assertEqual(ret_data, {})
        self.assertEqual(ret_info, {})

    #@unittest.skip('test_import_mmcif_file_pathched')
    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    def test_import_mmcif_file_patched(self, download_staging_file):
        fileName = '1fat.cif'
        params = {
            'input_shock_id': '',
            'input_file_path': '',
            'input_staging_file_path': fileName,
            'structure_name': 'test_pdb_structure_name',
            'description': 'for test PDBUtil.import_mmcif_file',
            'workspace_name': self.wsName
        }
        ret_data, ret_info = self.pdb_util.import_mmcif_file(params)
        self.assertEqual(ret_data, {})
        self.assertEqual(ret_info, {})

    # saveStructures_createReport
    #@unittest.skip('test_saveStructures_createReport')
    def test_saveStructures_createReport(self):
        struct_objects = {
            'protein_structures': [{
                'name': 'crystal structure of fdr9',
                'num_chains': 1,
                'num_residues': 393,
                'num_atoms': 6226,
                'compound': {'misc': '', 'molecule': 'putative oxidoreductase',
                             'chain': 'b', 'engineered': 'yes'},
                'source': {'misc': '', 'organism_scientific': 'thermobifida fusca yx',
                           'organism_taxid': '269800', 'gene': 'tfu_1273',
                           'expression_system': 'escherichia coli',
                           'expression_system_taxid': '562',
                           'expression_system_vector_type': 'plasmid',
                           'expression_system_plasmid': 'pet28a'},
                           'proteins': [{
                                'id': '6TUK.pdb.gz', 'model_id': 0, 'chain_id': 'B',
                                'sequence': 'MERIVIVGGGLAASRTCEQLRSRGYEGELVAEPHPPYDPSHORT',
                                'md5': 'a493e65f63ab6d3e4a7c1ff86c75bf0a',
                                'seq_identity': 0.68932,
                                'exact_match': 0,
                                'genome_ref': '63679/38/1',  # for appdev
                                # 'genome_ref': '107138/2/1',  # for prod
                                'feature_id': 'MLuteus_masurca_RAST.CDS.133',
                                'feature_type': 'gene'}],
                           'pdb_handle': 'KBH_209585',
                           'user_data': '',
                           'is_model': 1
            }],
            'total_structures': 1,'total_structures': 1,
            'description': 'ONE structure created'
        }
        pdb_infos = [{
            'rcsb_id': '6TUK',
            'extension': '.pdb',
            'narrative_id': 63679,
            'genome_name': 'MLuteus_ATCC_49442',
            'feature_id': 'MLuteus_masurca_RAST.CDS.133',
            'is_model': 1,
            'file_path': os.path.join('/kb/module/test/data', '6TUK.pdb.gz'),
            'genome_ref': '63679/38/1',  # for appdev
            # 'genome_ref': '107138/2/1',  # for prod
            'feature_type': 'gene',
            'sequence_identities': '68.93%',
            'chain_ids': 'Model 1.Chain B',
            'model_ids': '0',
            'exact_matches': '0',
            'scratch_path': os.path.join('/kb/module/test/data', '6TUK.pdb.gz')
        }]
        failed_ids = ['1A0I']

        structs_name = 'test_saveNreport_struct'
        if not isinstance(self.wsName, int):
            workspace_id = self.dfu.ws_name_to_id(self.wsName)
        else:
            workspace_id = self.wsName

        ret = self.pdb_util.saveStructures_createReport(structs_name, workspace_id, self.wsName,
                                                        struct_objects, pdb_infos, failed_ids,
                                                        rcsb=True)
        self.assertIn('batch_import_rcsb_report', ret['report_name'])

    #@unittest.skip('test_batch_import_pdbs_pathched')
    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    def test_batch_import_pdbs_patched(self, download_staging_file):
        fileName = 'pdb_metafile_sample2.csv'

        params = {
            'metadata_staging_file_path': os.path.join('/kb/module/test/data', fileName),
            'structures_name': 'batch_import_pdbs_structures',
            'workspace_name': self.wsName
        }
        ret = self.pdb_util.batch_import_pdbs(params)
        self.assertCountEqual(ret.keys(), ["structures_ref", "report_ref", "report_name"])
        #structs_ref = ret['structures_ref']
        #print(f'Checking the newly saved object data and info for {structs_ref}\n')

    # Testing self.serviceImpl functions
    ##@unittest.skip('test_Impl_batch_import_pdbs_from_metafile1')
    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    def test_Impl_batch_import_pdbs_from_metafile1(self, download_staging_file):
        metafile = 'pdb_metafile_sample1.csv'
        metafile = os.path.join('/kb/module/test/data', metafile)

        params = {
            'metadata_staging_file_path': metafile,
            'structures_name': 'batch1_test_structures',
            'workspace_name': self.wsName
        }
        ret = self.serviceImpl.batch_import_pdbs_from_metafile(self.ctx, params)
        self.assertCountEqual(ret[0].keys(), ["structures_ref", "report_ref", "report_name"])

    #@unittest.skip('test_Impl_batch_import_pdbs_from_metafile2')
    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    def test_Impl_batch_import_pdbs_from_metafile2(self, download_staging_file):
        metafile = 'pdb_metafile_sample2.csv'
        metafile = os.path.join('/kb/module/test/data', metafile)

        params = {
            'metadata_staging_file_path': metafile,
            'structures_name': 'batch2_test_structures',
            'workspace_name': self.wsName
        }
        ret1 = self.serviceImpl.batch_import_pdbs_from_metafile(self.ctx, params)
        self.assertCountEqual(ret1[0].keys(), ["structures_ref", "report_ref", "report_name"])

    #@unittest.skip('test_Impl_batch_import_pdbs_for_MLuteus_ATCC_alphafolds')
    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    def test_Impl_batch_import_pdbs_for_MLuteus_ATCC_alphafolds(self, download_staging_file):
        metafile = 'pdb_metafile_mluteus_afpdbs.csv'
        metafile = os.path.join('/kb/module/test/data', metafile)

        params = {
            'metadata_staging_file_path': metafile,
            'structures_name': 'mluteus_afpdb_structures',
            'workspace_name': self.wsName
        }
        ret1 = self.serviceImpl.batch_import_pdbs_from_metafile(self.ctx, params)
        self.assertCountEqual(ret1[0].keys(), ["structures_ref", "report_ref", "report_name"])
        parms = {'input_ref': ret1[0]['structures_ref']}
        exp_pdb_shockIDs = self.pdb_util.export_pdb_structures(parms)
        self.assertEqual(len(exp_pdb_shockIDs['shock_ids']), 12)
        self.assertCountEqual(exp_pdb_shockIDs.keys(), ['shock_ids'])

    #@unittest.skip('test_Impl_batch_import_pdbs_for_MLuteus_ATCC_cifs')
    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    def test_Impl_batch_import_pdbs_for_MLuteus_ATCC_cifs(self, download_staging_file):
        metafile = 'pdb_metafile_mluteus_afcifs.csv'
        metafile = os.path.join('/kb/module/test/data', metafile)

        params = {
            'metadata_staging_file_path': metafile,
            'structures_name': 'mluteus_afpdb_structures',
            'workspace_name': self.wsName
        }
        ret = self.serviceImpl.batch_import_pdbs_from_metafile(self.ctx, params)
        self.assertCountEqual(ret[0].keys(), ["structures_ref", "report_ref", "report_name"])
        parms = {'input_ref': ret[0]['structures_ref']}
        exp_pdb_shockIDs = self.pdb_util.export_pdb_structures(parms)
        self.assertEqual(len(exp_pdb_shockIDs['shock_ids']), 2)
        self.assertCountEqual(exp_pdb_shockIDs.keys(), ['shock_ids'])

    #@unittest.skip('test_Impl_batch_import_pdbs_for_MLuteus_ATCC_docked')
    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    def test_Impl_batch_import_pdbs_for_MLuteus_ATCC_docked(self, download_staging_file):
        metafile = 'pdb_metafile_mluteus_docked.csv'
        metafile = os.path.join('/kb/module/test/data', metafile)

        params = {
            'metadata_staging_file_path': metafile,
            'structures_name': 'mluteus_docked_structures',
            'workspace_name': self.wsName
        }
        ret2 = self.serviceImpl.batch_import_pdbs_from_metafile(self.ctx, params)
        self.assertCountEqual(ret2[0].keys(), ["structures_ref", "report_ref", "report_name"])
        parms = {'input_ref': ret2[0]['structures_ref']}
        exp_pdb_shockIDs = self.pdb_util.export_pdb_structures(parms)
        self.assertEqual(len(exp_pdb_shockIDs['shock_ids']), 12)
        self.assertCountEqual(exp_pdb_shockIDs.keys(), ['shock_ids'])

    #@unittest.skip('test_Impl_batch_import_pdbs_mix_rcsb')
    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    def test_Impl_batch_import_pdbs__mix_rcsb(self, download_staging_file):
        metafile = 'pdb_metafile_sample_mix_rcsb.csv'
        metafile = os.path.join('/kb/module/test/data', metafile)

        params = {
            'metadata_staging_file_path': metafile,
            'structures_name': 'mluteus_docked_structures',
            'workspace_name': self.wsName
        }
        ret2 = self.serviceImpl.batch_import_pdbs_from_metafile(self.ctx, params)
        self.assertCountEqual(ret2[0].keys(), ["structures_ref", "report_ref", "report_name"])
        parms = {'input_ref': ret2[0]['structures_ref']}
        exp_pdb_shockIDs = self.pdb_util.export_pdb_structures(parms)
        self.assertEqual(len(exp_pdb_shockIDs['shock_ids']), 3)
        self.assertCountEqual(exp_pdb_shockIDs.keys(), ['shock_ids'])

    # '57196/53/1' is in AppDev
    #@unittest.skip('test_export_pdb_structures')
    def test_export_pdb_structures(self):
        params = {'input_ref': '57196/53/1'}
        ret = self.serviceImpl.export_pdb_structures(self.ctx, params)
        self.assertCountEqual(ret[0].keys(), ['shock_ids'])
        self.assertEqual(len(ret[0]['shock_ids']), 12)

    # '62713/24/1' is in CI, so skipped here by starting the function without 'test_'.
    # ('ci_export_pdb_structures')
    def ci_export_pdb_structures_ci(self):
        params = {'input_ref': '62713/24/1'}
        ret = self.serviceImpl.export_pdb_structures(self.ctx, params)
        self.assertCountEqual(ret[0].keys(), ['shock_ids'])

    def dfu_save_proteinstructure(self):
        """Just for testing dfu saving a well-defined KBaseStructure.ProteinStructures"""

        obj_to_save = {
            'protein_structures': [
              {'name': 'ksga from bacillus subtilis 168', 'num_chains': 2, 'num_residues': 567, 'num_atoms': 4583,
               'compound': {'misc': '', 'molecule': 'ribosomal rna small subunit methyltransferase a', 'chain': 'a, b',
                            'synonym': "16s rrna (adenine(1518)-n(6)/adenine(1519)-n(6))-dimethyltransferase,16s rrna dimethyladenosine transferase,16s rrna dimethylase,s-adenosylmethionine-6-n',n'-adenosyl(rrna) dimethyltransferase ",
                            'ec_number': '2.1.1.182', 'ec': '2.1.1.182', 'engineered': 'yes'},
               'source': {'misc': '', 'organism_scientific': 'bacillus subtilis (strain 168)', 'organism_taxid': '224308', 'strain': '168', 'gene': 'rsma, ksga, bsu00420', 'expression_system': 'escherichia coli', 'expression_system_taxid': '562'},
               'proteins': [
                  {'id': '6ifs.pdb', 'model_id': 0, 'chain_id': 'A', 'sequence': 'KDIATPIRTKEILKKYGFSFKKSLGQNFLIDTNILNRIVDHAEVTEKTGVIEIGPGIGALTEQLAKRAKKVVAFEIDQRLLPILKDTLSPYENVTVIHQDVLKADVKSVIEEQFQDCDEIMVVANLPYYVTTPIIMKLLEEHLPLKGIVVMLQKEVAERMAADPSSKEYGSLSIAVQFYTEAKTVMIVPKTVFVPQPNVDSAVIRLILRDGPAVDVENESFFFQLIKASFAQRRKTLLNNLVNNLPEGKAQKSTIEQVLEETNIDGKRRGESLSIEEFAALSNGLYKALF', 'md5': 'f0eb43b0bf610adb501695053cecc5a6', 'seq_identity': 0.679487, 'exact_match': 0, 'genome_ref': '57196/6/1', 'feature_id': 'JCVISYN3_0004', 'feature_type': 'gene'},
                  {'id': '6ifs.pdb', 'model_id': 0, 'chain_id': 'B', 'sequence': 'ATPIRTKEILKKYNFLIDTNILNRIVDHAEVTEKTGVIEIGPGIGALTEQLAKRAKKVVAFEIDQRLLPILKDTLSPYENVTVIHQDVLKADVKSVIEEQFQDCDEIMVVANLPYYVTTPIIMKLLEEHLPLKGIVVMLQKEVAERMAADPSSKEYGSLSIAVQFYTEAKTVMIVPKTVFVPQPNVDSAVIRLILRDGPAVDVENESFFFQLIKASFAQRRKTLLNNLVNNLPEGKAQKSTIEQVLEETNIDGKRRGESLSIEEFAALSNGLYKALF', 'md5': 'a6c3a37bcc8c3be55fce6bfdae74e3be', 'seq_identity': 0.677419, 'exact_match': 0, 'genome_ref': '57196/6/1', 'feature_id': 'JCVISYN3_0004', 'feature_type': 'gene'}
                ],
                'pdb_handle': 'KBH_182761', 'user_data': '', 'is_model': 1
              },
              {'name': 'ksga from bacillus subtilis in complex with sam', 'num_chains': 1, 'num_residues': 292, 'num_atoms': 2508, 'compound': {'misc': '', 'molecule': 'ribosomal rna small subunit methyltransferase a', 'chain': 'a', 'synonym': "16s rrna (adenine(1518)-n(6)/adenine(1519)-n(6))-dimethyltransferase,16s rrna dimethyladenosine transferase,16s rrna dimethylase,s-adenosylmethionine-6-n',n'-adenosyl(rrna) dimethyltransferase ", 'ec_number': '2.1.1.182', 'ec': '2.1.1.182', 'engineered': 'yes'}, 'source': {'misc': '', 'organism_scientific': 'bacillus subtilis (strain 168)', 'organism_taxid': '224308', 'strain': '168', 'gene': 'rsma, ksga, bsu00420', 'expression_system': 'escherichia coli', 'expression_system_taxid': '562'},
               'proteins': [{'id': '6ift.pdb', 'model_id': 0, 'chain_id': 'A', 'sequence': 'MNKDIATPIRTKEILKKYGFSFKKSLGQNFLIDTNILNRIVDHAEVTEKTGVIEIGPGIGALTEQLAKRAKKVVAFEIDQRLLPILKDTLSPYENVTVIHQDVLKADVKSVIEEQFQDCDEIMVVANLPYYVTTPIIMKLLEEHLPLKGIVVMLQKEVAERMAADPSSKEYGSLSIAVQFYTEAKTVMIVPKTVFVPQPNVDSAVIRLILRDGPAVDVENESFFFQLIKASFAQRRKTLLNNLVNNLPEGKAQKSTIEQVLEETNIDGKRRGESLSIEEFAALSNGLYKALF', 'md5': '076fc7704e83dab1565ae973a11c435d', 'seq_identity': 0.679487, 'exact_match': 0, 'genome_ref': '57196/6/1', 'feature_id': 'JCVISYN3_0004', 'feature_type': 'gene'}],
               'pdb_handle': 'KBH_182762', 'user_data': '', 'is_model': 1
              },
              {'name': 'c-terminal truncated ksga from bacillus subtilis 168', 'num_chains': 2, 'num_residues': 410, 'num_atoms': 3068,'compound': {'misc': '', 'molecule': 'ribosomal rna small subunit methyltransferase a', 'chain': 'a, b', 'fragment': 'c-terminal truncated', 'synonym': "16s rrna (adenine(1518)-n(6)/adenine(1519)-n(6))-dimethyltransferase,16s rrna dimethyladenosine transferase,16s rrna dimethylase,s-adenosylmethionine-6-n',n'-adenosyl(rrna) dimethyltransferase ", 'ec_number': '2.1.1.182', 'ec': '2.1.1.182', 'engineered': 'yes'}, 'source': {'misc': '', 'organism_scientific': 'bacillus subtilis (strain 168)', 'organism_taxid': '224308', 'strain': '168', 'gene': 'rsma, ksga, bsu00420', 'expression_system': 'escherichia coli', 'expression_system_taxid': '562'},
               'proteins': [{'id': '6ifv.pdb', 'model_id': 0, 'chain_id': 'A', 'sequence': 'KDIATPIRTKEILKKYGFSFQNFLIDTNILNRIVDHAEVTEKTGVIEIGPGIGALTEQLAKRAKKVVAFEIDQRLLPILKDTLSPYENVTVIHQDVLKADVKSVIEEQFQDCDEIMVVANLPYYVTTPIIMKLLEEHLPLKGIVVMLQKEVAERMAADPSSKEYGSLSIAVQFYTEAKTVMIVPKTVFVPQPNVDSAVIRLILR', 'md5': '08a6b2ba57f21c334db58633867a40e9', 'seq_identity': 0.694915, 'exact_match': 0, 'genome_ref': '57196/6/1', 'feature_id': 'JCVISYN3_0004', 'feature_type': 'gene'},
                            {'id': '6ifv.pdb', 'model_id': 0, 'chain_id': 'B', 'sequence': 'DIATPIRTKEILKKYGFSFKKSQNFLIDTNILNRIVDHAEVTEKTGVIEIGPGIGALTEQLAKRAKKVVAFEIDQRLLPILKDTLSPYENVTVIHQDVLKADVKSVIEEQFQDCDEIMVVANLPYYVTTPIIMKLLEEHLPLKGIVVMLQKEVAERMAADPSSKEYGSLSIAVQFYTEAKTVMIVPKTVFVPQPNVDSAVIRLILR', 'md5': 'e614ffbf074ae76a9c3fd447adf9e447', 'seq_identity': 0.694915, 'exact_match': 0, 'genome_ref': '57196/6/1', 'feature_id': 'JCVISYN3_0004', 'feature_type': 'gene'}],
               'pdb_handle': 'KBH_182763', 'user_data': '', 'is_model': 1
              },
              {'name': 'crystal structure of chimeric construct of ksga with loop 1 from erm', 'num_chains': 2, 'num_residues': 486, 'num_atoms': 3794, 'compound': {'misc': '', 'molecule': 'ribosomal rna small subunit methyltransferase a', 'chain': 'a, b', 'synonym': "16s rrna (adenine(1518)-n(6)/adenine(1519)-n(6))-dimethyltransferase,16s rrna dimethyladenosine transferase,16s rrna dimethylase,s-adenosylmethionine-6-n',n'-adenosyl(rrna) dimethyltransferase ", 'ec_number': '2.1.1.182', 'ec': '2.1.1.182', 'engineered': 'yes'}, 'source': {'misc': '', 'organism_scientific': 'bacillus subtilis (strain 168)', 'organism_taxid': '224308', 'strain': '168', 'gene': 'rsma, ksga, bsu00420', 'expression_system': 'escherichia coli', 'expression_system_taxid': '562'},
               'proteins': [{'id': '6ifw.pdb', 'model_id': 0, 'chain_id': 'A', 'sequence': 'NFLIDTNILNRIVDHAEVTEKTGVIEIGPGIGALTEQLAKRAKKVVAFEIDQRLLPILKDTLSPYENVTVIHQDVLKADVKSVIEEQFQDCDEIMVVANLPYYVTTPIIMKLLEEHLPLKGIVVMLQKEVAERMAADPSSKEYGSLSIAVQFYTEAKTVMIVPKTVFVPQPNVDSAVIRLILRDGPAVDVENESFFFQLIKASFAQRRKTLLNNLVNNLPEGKAQKSTIEQVLEETNIDGKRRGESLSIEEFAALSNGLYKALF', 'md5': '662e3cf38381c833487b97de1ad4a5f2', 'seq_identity': 0.671053, 'exact_match': 0, 'genome_ref': '57196/6/1', 'feature_id': 'JCVISYN3_0004', 'feature_type': 'gene'},
                            {'id': '6ifw.pdb', 'model_id': 0, 'chain_id': 'B', 'sequence': 'QNFLIDTNILNRIVDHAEVTEKTGVIEIGPGIGALTEQLAKRAKKVVAFEIDQRLLPILKDTLSPYENVTVIHQDVLKADVKSVIEEQFQDCDEIMVVANLPYYVTTPIIMKLLEEHLPLKGIVVMLQKEVAERMAADPSSKEYGSLSIAVQFYTEAKTVMIVPKTVFVPQPNVDSAVIRLILRDGPAVDVENESFFFQLIKASFNNLVNNLSIEEFAALSN', 'md5': 'c7f304f9a4fb589e19b516ea1bf3a36c', 'seq_identity': 0.671756, 'exact_match': 0, 'genome_ref': '57196/6/1', 'feature_id': 'JCVISYN3_0004', 'feature_type': 'gene'}],
               'pdb_handle': 'KBH_182764', 'user_data': '', 'is_model': 1
              }
            ],
            'total_structures': 4,
            'description': 'Created 4 structures in batch2_test_structures'
        }

        try:
            info = self.dfu.save_objects({
                'id': self.ws_id,
                'objects': [
                    {'type': 'KBaseStructure.ProteinStructures',
                     'name': 'structures_test',
                     'data': obj_to_save}]
            })[0]
        except (RuntimeError, TypeError, KeyError, ValueError) as e:
            err_msg = f'DFU.save_objects errored with message: {e.message} and data: {e.data}'
            print(err_msg)
            raise ValueError(err_msg)
        else:
            structs_ref = f"{info[6]}/{info[0]}/{info[4]}"
            self.check_object(structs_ref)
            return structs_ref

    def check_object(self, obj_ref):
        print(f'Checking the given object data and info for {obj_ref}\n')
        obj_data, obj_info = self.pdb_util._dfu_get_objects(obj_ref)
        print(f'\n------------------------Data-------------------------------\n{obj_data}')
        print(f'\n------------------------Info-------------------------------\n{obj_info}')
