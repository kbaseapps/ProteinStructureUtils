import hashlib
import logging
import os
import sys
import shutil
import uuid
import errno
import pandas as pd
import pathlib
import subprocess
from urllib.parse import urlparse

from Bio import PDB
from Bio.PDB.Polypeptide import PPBuilder
# from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from installed_clients.AbstractHandleClient import AbstractHandle
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace
from installed_clients.baseclient import ServerError as WorkspaceError


class PDBUtil:

    E_VALUE_THRESH = 1e-20
    # staging file prefix
    STAGING_GLOBAL_FILE_PREFIX = '/data/bulk/'
    STAGING_USER_FILE_PREFIX = '/staging/'
    PDB_FILE_EXT = ['pdb', 'cif']

    def _validate_import_pdb_file_params(self, params):
        """
            _validate_import_pdb_file_params:
                validates input params to import_model_pdb_file and import_experiment_pdb_file
        """
        # check for required parameters
        for p in ['structure_name', 'workspace_name']:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

        if params.get('input_file_path'):
            file_path = params.get('input_file_path')
        elif params.get('input_shock_id'):
            file_path = self.dfu.shock_to_file(
                {'shock_id': params['input_shock_id'],
                 'file_path': self.scratch}).get('file_path')
        elif params.get('input_staging_file_path'):
            file_path = self.dfu.download_staging_file(
                {'staging_file_subdir_path': params.get('input_staging_file_path')}
                 ).get('copy_file_path')
        else:
            error_msg = "Must supply either a input_shock_id or input_file_path "
            error_msg += "or input_staging_file_path"
            raise ValueError(error_msg)

        return file_path, params.get('workspace_name'), params.get('structure_name')

    def _model_file_to_data(self, file_path, params):
        """
            _model_file_to_data:
                Do the PDB conversion--parse the model pdb file for creating a pdb data object
        """
        logging.info(f'Parsing pdb file {file_path} to a pdb structure with params: {params}')

        parser = PDB.PDBParser(PERMISSIVE=1)
        pdb1 = file_path
        structure = parser.get_structure("test", pdb1)

        pp_no = 0
        ppb = PPBuilder()
        for pp in ppb.build_peptides(structure):
            pp_no += 1

        (compound, source) = self._get_compound_source(structure)
        (num_models, model_ids) = self._get_models_from_structure(structure)
        (num_chains, chain_ids) = self._get_chains_from_structure(structure)
        (num_residues, residue_ids) = self._get_residues_from_structure(structure)
        (num_atoms, atom_ids) = self._get_atoms_from_structure(structure)
        model = structure[0]
        protein_data = self._get_proteins_by_structure(structure, model.get_id(), file_path)
        (protein_data, params) = self._match_features(params, protein_data)

        data = {
            'name': structure.header.get('name', ''),
            'num_chains': num_chains,
            'num_residues': num_residues,
            'num_atoms': num_atoms,
            'compound': compound,
            'source': source,
            'proteins': protein_data
        }

        return data, pp_no, params

    def _exp_file_to_data(self, file_path, params):
        """
            _exp_file_to_data:
                Do the PDB conversion--parse the experiment pdb file for creating a pdb data object
        """
        logging.info(f'Parsing pdb file {file_path} to a pdb structure with params: {params}')

        parser = PDB.MMCIFParser()
        cif = file_path
        structure = parser.get_structure("PHA-L", cif)
        # create a dictionary that maps all mmCIF tags in an mmCIF file to their values
        # mmcif_dict = MMCIF2Dict(cif)  # TypeError: 'module' object is not callable

        pp_no = 0
        ppb = PPBuilder()
        for pp in ppb.build_peptides(structure):
            pp_no += 1

        struc_name = structure.header.get('name', '')
        hd = self._upload_to_shock(file_path)

        (cpd, src) = self._get_compound_source(structure)
        (num_models, model_ids) = self._get_models_from_structure(structure)
        (num_chains, chain_ids) = self._get_chains_from_structure(structure)
        (num_residues, residue_ids) = self._get_residues_from_structure(structure)
        (num_atoms, atom_ids) = self._get_atoms_from_structure(structure)
        protein_data = self._get_proteins_by_structure(structure, model_ids[0], file_path)
        (protein_data, params) = self._match_features(params, protein_data)

        mmcif_data = {
            'name': struc_name,
            'head': structure.header.get('head', ''),
            'rcsb_id': structure.header.get('rcsb_id', ''),
            'deposition_date': structure.header.get('deposition_date', ''),
            'release_date': structure.header.get('release_date', ''),
            'structure_method': structure.header.get('structure_method', ''),
            'resolution': structure.header.get('resolution', 0.0),
            'structure_reference': structure.header.get('structure_reference', []),
            'keywords': structure.header.get('keywords', ''),
            'author': structure.header.get('author', ''),
            'compound': cpd,
            'source': src,
            'num_models': num_models,
            'num_chains': num_chains,
            'num_residues': num_residues,
            'num_atoms': num_atoms,
            'num_het_atoms': structure.header.get('num_het_atoms', 0),
            'num_water_atoms': structure.header.get('num_water_atoms', 0),
            'num_disordered_atoms': structure.header.get('num_disordered_atoms', 0),
            'num_disordered_residues': structure.header.get('num_disordered_residues', 0),
            'pdb_handle': hd,
            'mmcif_handle': hd,
            'xml_handle': hd,
            'proteins': protein_data
        }

        return mmcif_data, pp_no, params

    def _match_features(self, params, protein_data):
        """
            _match_features: match the protein_translation in feature_id with chain sequences in
                             protein_data and compute the seq_identity and determine the exact_match
            example (in appdev):
                    genome_obj = '57196/6/1', genome_name = 'Synthetic_bacterium_JCVI_Syn3.0_genome'
                    feature_id = 'JCVISYN3_0004_CDS_1', feature_type = 'CDS' OR
                    feature_id = 'JCVISYN3_0004', feature_type = 'gene'
        """
        pdb_info = params.get('pdb_info', None)
        if pdb_info:
            feature_type = ''
            feat_prot_seq = ''
            genome_name = pdb_info['genome_name']
            narr_id = pdb_info['narrative_id']
            feature_id = pdb_info['feature_id']

            logging.info(f"Looking up for feature {feature_id} in genome {genome_name}'s features")
            # 1. Get the genome's features and reference
            (gn_ref, gn_features) = self._get_genome_ref_features(narr_id, genome_name)
            pdb_info['genome_ref'] = gn_ref

            # 2. Match the genome features with the specified feature_id to obtain feature sequence
            for feat in gn_features:
                if feat['id'] == feature_id:
                    logging.info(f'Found genome feature match for {feature_id}')
                    prot_trans = feat['protein_translation']
                    feat_prot_seq = prot_trans
                    feature_type = self._get_feature_type(feat)
                    break
            pdb_info['feature_type'] = feature_type

            # 3. Call self._compute_sequence_identity with the feature sequence and the the pdb
            # proteins' translationsto to get the seq_identity and exact_match
            if feat_prot_seq:
                pdb_info['protein_meta'] = []
                for prot in protein_data:
                    seq_idens, seq_mats = self._compute_sequence_identity(feat_prot_seq,
                                                                          prot.get('sequence', ''))
                    if seq_idens:
                        seq_idens.sort()
                        max_iden = seq_idens.pop()
                        prot['seq_identity'] = max_iden
                        prot['exact_match'] = 1 if max_iden > 0.99 else 0
                        prot['genome_ref'] = gn_ref
                        prot['feature_id'] = feature_id
                        prot['feature_type'] = feature_type
                        pdb_info['protein_meta'].append({
                            'model_id': prot['model_id'],
                            'chain_id': prot['chain_id'],
                            'seq_identity': prot['seq_identity'],
                            'exact_match': prot['exact_match']})
            else:
                logging.info(f'Found NO feature in genome that matches with {feature_id}')

        return protein_data, params

    def _compute_sequence_identity(self, seq1, seq2):
        """
            _compute_sequence_identity: Given two input sequences, do a blast identity check and
                                        then compute and return the matching percentage.
        """
        # Create two sequence files
        Seq1 = SeqRecord(Seq(seq1), id="query_seq")
        Seq2 = SeqRecord(Seq(seq2), id="subject_seq")

        blast_dir = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(blast_dir)
        query_seq = os.path.join(blast_dir, 'seq_qry.fasta')
        subject_seq = os.path.join(blast_dir, 'seq_sbj.fasta')
        SeqIO.write(Seq1, query_seq, "fasta")
        SeqIO.write(Seq2, subject_seq, "fasta")

        # on my laptop: blastp_path = '/Users/qzhang/miniconda3/bin/blastp'
        blastp_path = 'blastp'
        output_file_path = os.path.join(blast_dir, 'blast_output.xml')

        # Build the BLASTp command
        blastp_cmd = [blastp_path]
        blastp_cmd.append('-out')
        blastp_cmd.append(output_file_path)
        blastp_cmd.append('-outfmt')
        blastp_cmd.append('5')
        blastp_cmd.append('-query')
        blastp_cmd.append(query_seq)
        blastp_cmd.append('-subject')
        blastp_cmd.append(subject_seq)

        # Run BLASTp and parse the output as XML and then parse the xml file for identity matches
        exact_matches = []
        idens = []
        try:
            p = subprocess.Popen(blastp_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 universal_newlines=True)
            output, errors = p.communicate()
            if not output:
                logging.info(f'BLASTp returned: {p.returncode}')
                logging.info(f'OK> output: {output}')
            if errors:
                e = subprocess.CalledProcessError(p.returncode, blastp_cmd, output=output)
                raise e
        except OSError as e:
            logging.info(f'OSError > {e.errno}')
            logging.info(f'OSError > {e.strerror}')
            logging.info(f'OSError > {e.filename}')
        except subprocess.CalledProcessError as e:
            logging.info(f'CalledError > {e.returncode}')
            logging.info(f'CalledError > {e.output}')
        except:
            logging.info(f'Unexpected error > {sys.exc_info()[0]}')
        else:
            with open(output_file_path) as blast_fhd:
                blast_record = NCBIXML.read(blast_fhd)
                if blast_record:
                    logging.info(f'query: {blast_record.query[:100]}')
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            if hsp.expect < self.E_VALUE_THRESH:
                                logging.info('****Alignment****')
                                logging.info(f'sequence: {alignment.title}')
                                logging.info(f'length: {alignment.length}')
                                logging.info(f'e value: {hsp.expect}')
                                logging.info(f'hsp query: {hsp.query}')
                                logging.info(f'hsp match: {hsp.match}')
                                logging.info(f'hsp subject: {hsp.sbjct}')
                                logging.info(f'hsp identities: {hsp.identities}')
                                logging.info(f'hsp positives: {hsp.positives}')
                                iden = round(hsp.identities/hsp.positives, 6)
                                logging.info(f'identity={iden}')
                                idens.append(iden)
                                if hsp.positives == hsp.identities:
                                    exact_matches.append(alignment.title[:100])
        return idens, exact_matches

    def _get_genome_ref_features(self, narr_id, genome_name):
        """
            _get_genome_ref_features: Get the genome reference and features for the given genome_name
        """
        genome_ref = ''
        genome_features = []
        (genome_info, genome_data) = self._get_object_info_data(narr_id, genome_name)
        genome_ref = '/'.join([str(narr_id), str(genome_info[0]), str(genome_info[4])])
        logging.info(f"Genome {genome_name}'s ref: {genome_ref}")
        genome_features = genome_data['features']
        logging.info(f'There are {len(genome_features)} features in {genome_name}')

        return (genome_ref, genome_features)

    def _get_feature_type(self, feature_obj):
        """
            _get_feature_type: Get the type for the feature object of given feature_obj
        """
        feat_type = feature_obj.get('type', '')
        if not feat_type:
            if feature_obj.get('protein_translation'):
                feat_type = 'gene'
            else:
                feat_type = 'other'

        return feat_type

    def _get_object_info_data(self, narr_id, obj_name):
        """
            _get_object_info_data: Get the object info/data with given obj_name in narrative narr_id
        """
        obj_info = None
        obj_data = None
        if narr_id and obj_name:
            obj_data_res = self.ws_client.get_objects2(
                {'objects': [{'wsid': narr_id, 'name': obj_name}]})['data'][0]
            obj_info = obj_data_res['info']
            obj_data = obj_data_res['data']

        return (obj_info, obj_data)

    def _get_atoms_from_structure(self, pdb_structure):
        """
            _get_atoms_from_structure: Given a pdb_structure object, parse atoms into a list of
                                        atoms and return it
        """
        atom_ids = []
        num_atoms = 0
        my_residues = pdb_structure.get_residues()
        for r_ele in my_residues:
            for a_ele in r_ele.get_atoms():
                num_atoms += 1
                atom_ids.append(a_ele.get_id())

        return (num_atoms, atom_ids)

    def _get_residues_from_structure(self, pdb_structure):
        """
            _get_residues_from_structure: Given a pdb_structure object, parse residues into a list
                                          and return it
        """
        res_ids = []
        num_res = 0
        my_res = pdb_structure.get_residues()
        for r_ele in my_res:
            if PDB.is_aa(r_ele):
                num_res += 1
                res_ids.append(r_ele.get_id())

        return (num_res, res_ids)

    def _get_chains_from_structure(self, pdb_structure):
        """
            _get_chains: Given a pdb_structure object, parse chain ids into a list and return it
        """
        chain_ids = []
        num_chains = 0
        my_chains = pdb_structure.get_chains()
        for c_ele in my_chains:
            if(c_ele):
                num_chains += 1
                chain_ids.append(c_ele.get_id())

        return (num_chains, chain_ids)

    def _get_models_from_structure(self, pdb_structure):
        """
            _get_models_from_structure: Given a pdb_structure object, parse model ids into a list
                                        and return it
        """
        model_ids = []
        num_models = 0
        my_models = pdb_structure.get_models()
        for m_ele in my_models:
            if(m_ele):
                num_models += 1
                model_ids.append(m_ele.get_id())

        return (num_models, model_ids)

    def _get_compound_source(self, structure):
        """
            _get_compound_source: Parse data from given structure for compound and source
        """
        cpd_dict = dict()
        cpd = structure.header.get('compound', {})
        # logging.info(f'Compound:\n {cpd}')
        if cpd and cpd.get('1'):
            cpd_dict = cpd.get('1')

        src_dict = dict()
        src = structure.header.get('source', {})
        # logging.info(f'Source:\n {src}')
        if src and src.get('1'):
            src_dict = src.get('1')

        return (cpd_dict, src_dict)

    def _get_proteins_by_structure(self, pdb_structure, model, file_path):
        """
            _get_proteins_by_structure: Given a pdb_structure, parse the essential protein data
        """
        ppb = PPBuilder()
        protein_data = []

        # Parse for the chain_id and chain sequence
        for c_ele in pdb_structure.get_chains():
            if(c_ele):
                c_ppd_list = []
                for c_ppd in ppb.build_peptides(c_ele):
                    c_pp_seq = str(c_ppd.get_sequence())
                    c_ppd_list.append(c_pp_seq)
                c_seq = ''.join(c_ppd_list)
                protein_data.append({
                    'id': os.path.basename(file_path),
                    'model_id': model,
                    'chain_id': c_ele.get_id(),
                    'sequence': c_seq,
                    'md5': hashlib.md5(c_seq.encode()).hexdigest()
                })

        return protein_data

    def _validate_file(self, file_path):
        """
            _validate_file: Check if file_path is accessable, if yes, return the handle
        """
        try:
            fh = open(file_path, 'r')
        except IOError as e:
            if e.errno == errno.ENOENT:  # No such file or directory
                raise ValueError(f'"{file_path}" does not exist!')
            elif e.errno == errno.EACCES:  # Permission denied
                raise ValueError(f'"{file_path}" cannot be read!')
            else:
                raise ValueError(f'"{e.strerror}" error occurred')
        else:
            fh.close()
            return True

    def _get_pdb_shock_id(self, obj_ref):
        """
            _get_pdb_shock_id: Return the shock id for the PDB file
        """
        obj_data = self.dfu.get_objects({"object_refs": [obj_ref]})['data'][0]['data']
        return self.hs.hids_to_handles([obj_data['pdb_handle']])[0]['id']

    def _upload_to_shock(self, file_path):
        """
            _upload_to_shock: upload target file to shock using DataFileUtil
        """
        logging.info(f'Start uploading file to shock: {file_path}')

        file_to_shock_params = {
            'file_path': file_path,
            'pack': 'gzip',
            'make_handle': True,
        }
        shock_id = self.dfu.file_to_shock(file_to_shock_params)['handle']['hid']

        return shock_id

    def _generate_report_html(self, pdb_name, pdb_path):
        """
            _generate_report_html: generates the HTML for the upload report
        """
        html_report = list()

        # Make report directory and copy over files
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(output_directory)
        result_file_path = os.path.join(output_directory, 'viewer.html')
        new_pdb_path = os.path.join(output_directory, os.path.basename(pdb_path))
        shutil.copy(pdb_path, new_pdb_path)

        # Fill in template HTML
        with open(os.path.join(os.path.dirname(__file__), 'templates', 'viewer_template.html')
                  ) as report_template_file:
            report_template = report_template_file.read()\
                .replace('*PDB_NAME*', pdb_name)\
                .replace('*PDB_PATH*', os.path.basename(pdb_path))

        with open(result_file_path, 'w') as result_file:
            result_file.write(report_template)

        html_report.append({'path': output_directory,
                            'name': os.path.basename(result_file_path),
                            'description': 'HTML report for PDB upload'})

        return html_report

    def _generate_report(self, method_name, pdb_obj_ref, workspace_name,
                         n_poly_pep, pdb_name, pdb_path):
        """
            _generate_report: generate summary report for upload
        """
        output_html_files = self._generate_report_html(pdb_name, pdb_path)

        report_params = {'message': f'You uploaded a PDB file. {n_poly_pep} polypeptides detected.',
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'objects_created': [{'ref': pdb_obj_ref,
                                              'description': 'Imported PDB'}],
                         'workspace_name': workspace_name,
                         'report_object_name': method_name + '_report_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _validate_batch_import_pdbs_params(self, params):
        """
            _validate_batch_import_pdbs_params:
                validates params passed to batch_import_pdbs method
        """
        # check for required parameters
        for p in ['structures_name', 'workspace_name', 'metadata_staging_file_path']:
            if p not in params:
                raise ValueError(f'"{p}" parameter is required, but missing')

        # metadata_staging_file_path must be from the staging area--must have the staging dir prefix
        if params.get('metadata_staging_file_path', None):
            staging_file_path = self.dfu.download_staging_file(
                {'staging_file_subdir_path': params.get('metadata_staging_file_path')}
                 ).get('copy_file_path')
            return (staging_file_path, params['workspace_name'], params['structures_name'])
        else:
            error_msg = "Must supply a 'metadata_staging_file_path'"
            raise ValueError(error_msg)

    def _read_file_by_type(self, file_path):
        """
            _read_file_by_type: read the file given by file_path depending on its type,
                               return a DataFrame object
        """
        logging.info(f'Reading input from file: {file_path}...')

        if not self._validate_file(file_path):
            raise ValueError('Input file is invalid or not found')

        df = None
        file_ext = pathlib.Path(file_path).suffix
        try:  # read the data from file_path depending on its extension
            if 'csv' in file_ext:
                df = pd.read_csv(file_path)
            elif 'tsv' in file_ext:
                df = pd.read_csv(file_path, '\t')
            elif 'xls' in file_ext or 'od' in file_ext:
                # handle xls, xlsx, xlsm, xlsb, odf, ods and odt file extensions
                df = pd.read_excel(file_path, index_col=None, engine='openpyxl')
            else:  # invalid file type
                error_msg = "Invalid input file type, only 'csv/tsv/xlsx' are accepted"
                raise ValueError(error_msg)
            # strip off the leading and trailing whitespaces of the column names
            df.columns = df.columns.str.strip()
        except WorkspaceError as e:
            logging.info(
                f'Reading file {file_path} errored with message: {e.message} and data: {e.data}')
            raise
        return df

    def _parse_metadata_file(self, metadata_file_path, ws_id):
        """
            _parse_metadata_file:
                From metadata_file_path, a spreadsheet file, sort out the model_pdb_file_paths,
            exp_pdb_file_paths and the kbase_meta_data

            return: lists model_pdb_file_paths, exp_pdb_file_paths and dict kbase_meta_data
        """
        logging.info(f'parsing metadata from input file {metadata_file_path}...')

        required_columns = ['Narrative ID', 'Object name (Genome AMA feature set)', 'Feature ID',
                            'PDB filename', 'Is model']

        pdb_file_paths = list()
        narrative_ids = list()
        genome_names = list()
        feature_ids = list()
        PDB_filenames = list()

        # df_meta_data is a Panda DataFrame object
        df_meta_data = self._read_file_by_type(metadata_file_path)
        df_col_list = df_meta_data.columns.values.tolist()

        # check if required columns are read in correctly
        for col in required_columns:
            if col not in df_col_list:
                missing_required = f"Required column '{col}' is missing!"
                raise ValueError(missing_required)

        df_indexes = df_meta_data.columns
        for i in range(len(df_meta_data[df_indexes[0]])):
            narr_id = int(df_meta_data[df_indexes[0]][i])
            if not pd.isna(narr_id):
                narrative_ids.append(narr_id)
            else:
                missing_narr_id = "Please fill all the rows in column 'Narrative ID'!"
                raise ValueError(missing_narr_id)

            obj_name = df_meta_data[df_indexes[1]][i]
            if not pd.isna(obj_name):
                genome_names.append(obj_name)
            else:
                missing_obj_name = "Please fill all the rows in column 'Object name'!"
                raise ValueError(missing_obj_name)

            feat_id = df_meta_data[df_indexes[2]][i]
            if not pd.isna(feat_id):
                feature_ids.append(feat_id)
            else:
                missing_feature_id = f"Please fill all the rows in column '{required_columns[2]}'!"
                raise ValueError(missing_feature_id)

            pdb_fn = df_meta_data[df_indexes[3]][i]  # pdb_fn does not have staging dir prefix
            if not pd.isna(pdb_fn):
                PDB_filenames.append(pdb_fn)
            else:
                missing_pdb_file = f"Please fill all the rows in column '{required_columns[3]}'!"
                raise ValueError(missing_pdb_file)

            f_name = os.path.basename(pdb_fn)
            (struct_name, ext) = os.path.splitext(f_name)
            is_model = df_meta_data[df_indexes[4]][i]
            if not pd.isna(is_model):
                pdb_file_paths.append(
                    {'file_path': self._get_staging_file_path(self.user_id, pdb_fn),
                     'structure_name': struct_name,
                     'is_model': 'y' in is_model or 'Y' in is_model,
                     'narrative_id': narr_id,
                     'genome_name': obj_name,
                     'feature_id': feat_id}
                )
            else:
                missing_pdb_md = f"Please fill all the rows in columns '{required_columns[4]}'!"
                raise ValueError(missing_pdb_md)

        if not pdb_file_paths:
            error_msg = "No PDB file info is provided!"
            raise ValueError(error_msg)

        return (pdb_file_paths, narrative_ids, genome_names, feature_ids)

    def _get_staging_file_path(self, token_user, staging_file_subdir_path):
        """
            _get_staging_file_path: return staging area file path
            directory pattern:
            perfered to return user specific path: /staging/sub_dir/file_name
            If this path is not visible to user, use global bulk path:
                        /data/bulk/user_name/sub_dir/file_name
        """
        user_path = os.path.join(self.STAGING_USER_FILE_PREFIX, staging_file_subdir_path.strip('/'))

        if os.path.exists(user_path):
            self.staging_path_prefix = self.STAGING_USER_FILE_PREFIX
            return user_path
        else:
            self.staging_path_prefix = os.path.join(self.STAGING_GLOBAL_FILE_PREFIX, token_user)
            return os.path.join(self.STAGING_GLOBAL_FILE_PREFIX, token_user,
                                staging_file_subdir_path.strip('/'))

    def _generate_batch_report(self, workspace_name, structs_ref, structs_name,
                               successful_paths, failed_paths):
        """
            _generate_batch_report: generate summary report for upload
        """
        output_html_files = self._generate_batch_report_html(structs_name, successful_paths)

        failed_files = ','.join(failed_paths)
        description = (f'Imported PDBs into a ProteinStructures object {structs_ref}, '
                       f'named {structs_name}')
        if failed_files:
            description += f', with files "{failed_files}" failed to load.'
        report_params = {'message': f'You have uploaded a batch of PDB files into {structs_name}.',
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'objects_created': [{'ref': structs_ref,
                                              'description': description}],
                         'workspace_name': workspace_name,
                         'report_object_name': 'batch_import_pdb_files_report_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _write_pdb_htmls(self, output_dir, succ_pdb_paths):
        """ Write the batch pdb info in a 3 column table of 'cards' into HTML files"""
        pdb_html = '<h2>PDBs successfully uploaded</h2>'
        row_html = ''
        srv_domain = urlparse(self.shock_url).netloc  # parse url to get the domain portion
        srv_base_url = f'https://{srv_domain}'
        logging.info(f'Get the url for building the anchors: {srv_base_url}')

        for i in range(len(succ_pdb_paths)):
            succ_pdb = succ_pdb_paths[i]
            file_path = succ_pdb['file_path']
            file_name = os.path.basename(file_path);
            new_pdb_path = os.path.join(output_dir, os.path.basename(file_path))
            shutil.copy(file_path, new_pdb_path)
            struct_name = succ_pdb['structure_name']
            struct_ref = succ_pdb['pdb_struct_ref']
            genome_name = succ_pdb['genome_name']
            genome_ref = succ_pdb['genome_ref']
            feat_id = succ_pdb['feature_id']
            feat_type = succ_pdb['feature_type']
            print(succ_pdb)

            pdb_chains = []
            seq_idens = []
            for prot in succ_pdb['protein_meta']:
                pdb_chains.append(prot['chain_id'])
                seq_idens.append(prot['seq_identity'])

            col_count = 3
            if i % col_count == 0:  # write col_count columns
                if row_html:  # close the row tags before opening a new row tag
                    row_html += '</div>'
                row_html += '<div class="row">'

            row_html += '<div class="column"><div class="card">'
            row_html += (f'<h3><a href="{srv_base_url}/#dataview/{struct_ref}" target="_blank">'
                         f'{struct_name} Data View</a> or <a href="#" '
                         f'onclick="showModal(\'{struct_name}\',\'{file_name}\');return false;">'
                         f'Structure View</a></h3>'
                         f'<p>Genome: <a href="{srv_base_url}/#dataview/{genome_ref}"'
                         f' target="_blank">{genome_name}</a></p>'
                         f'<p>Feature: {feat_id}, type: {feat_type}</p>')
            if pdb_chains:
                row_html += f'<p>Chains: {pdb_chains}</p>'
            if seq_idens:
                row_html += f'<p>Sequence identity: {seq_idens}</p>'
            row_html += '</div></div>'
        row_html += '</div>'
        pdb_html += row_html
        return pdb_html

    def _generate_batch_report_html(self, prot_structs_name, succ_pdb_paths):
        """
            _generate_batch_report_html: generates the HTML for the upload report
        """
        html_report = list()

        # Make report directory and copy over uploaded pdb files
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(output_directory)

        # Create the template html file for reporting batch-uploaded pdb files
        batch_html_report_path = os.path.join(output_directory, 'batch_pdb_viewer.html')

        pdb_html = self._write_pdb_htmls(output_directory, succ_pdb_paths)

        # Fetch & fill in detailed info into template HTML
        with open(os.path.join(os.path.dirname(__file__), 'templates', 'batch_pdb_template.html')
                  ) as batch_template_html:
            batch_html_report = batch_template_html.read()\
                .replace('<!--replace this content-->', pdb_html)

        with open(batch_html_report_path, 'w') as html_report_file:
            html_report_file.write(batch_html_report)
        print(f'Full batch_html_report has been written to {batch_html_report_path}')

        html_report.append({'path': output_directory,
                            'name': os.path.basename(batch_html_report_path),
                            'description': 'HTML report for PDB upload'})

        return html_report

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.token = config['KB_AUTH_TOKEN']
        self.user_id = config['USER_ID']
        self.dfu = DataFileUtil(self.callback_url)
        self.hs = AbstractHandle(config['handle-service-url'])
        self.ws_client = Workspace(config['workspace-url'])
        self.shock_url = config['shock-url']

    def import_model_pdb_file(self, params, create_report=True):
        """
            import_model_pdb_file: upload an experiment pdb file and convert into a
                                  KBaseStructure.ModelProteinStructure object
        """
        logging.info(f'Importing pdb file to a pdb structure with params: {params}')

        file_path, workspace_name, pdb_name = self._validate_import_pdb_file_params(params)

        if not isinstance(workspace_name, int):
            workspace_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            workspace_id = workspace_name

        (data, n_polypeptides, params) = self._model_file_to_data(file_path, params)
        data['pdb_handle'] = self._upload_to_shock(file_path)
        data['user_data'] = params.get('description', '')
        logging.info(data)

        try:
            info = self.dfu.save_objects({
                'id': workspace_id,
                'objects': [
                    {'type': 'KBaseStructure.ModelProteinStructure',
                     'name': pdb_name,
                     'data': data}]
            })[0]
        except WorkspaceError as e:
            logging.info(f'DFU.save_objects errored with message: {e.message} and data: {e.data}')
            raise
        else:
            obj_ref = f"{info[6]}/{info[0]}/{info[4]}"
            returnVal = {'structure_obj_ref': obj_ref}
            if params.get('pdb_info'):
                returnVal['pdb_info'] = params.get('pdb_info')

            if create_report:
                report_output = self._generate_report('import_model_pdb_file', obj_ref,
                                                      workspace_name, n_polypeptides,
                                                      pdb_name, file_path)
                returnVal.update(report_output)
            return returnVal

    def import_experiment_pdb_file(self, params, create_report=True):
        """
            import_experiment_pdb_file: upload an experiment pdb file and convert into a
                                       KBaseStructure.ExperimentalProteinStructure object
        """

        file_path, workspace_name, mmcif_name = self._validate_import_pdb_file_params(params)

        if not isinstance(workspace_name, int):
            workspace_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            workspace_id = workspace_name

        # Parse the experimental pdb file for an experimental data structure
        (data, n_polypeptides, params) = self._exp_file_to_data(file_path, params)
        logging.info(data)

        try:
            info = self.dfu.save_objects({
                'id': workspace_id,
                'objects': [
                    {'type': 'KBaseStructure.ExperimentalProteinStructure',
                     'name': mmcif_name,
                     'data': data}]
            })[0]
        except WorkspaceError as e:
            logging.info(f'DFU.save_objects errored with message: {e.message} and data: {e.data}')
            raise
        else:
            obj_ref = f"{info[6]}/{info[0]}/{info[4]}"
            returnVal = {'structure_obj_ref': obj_ref}
            if params.get('pdb_info', None):
                returnVal['pdb_info'] = params.get('pdb_info')

            if create_report:
                report_output = self._generate_report('import_experiment_pdb_file', obj_ref,
                                                      workspace_name, n_polypeptides,
                                                      mmcif_name, file_path)
                returnVal.update(report_output)
            return returnVal

    def export_pdb(self, params):
        """
            export_pdb: return the shock_id of the uploaded pdb object
        """
        if "input_ref" not in params:
            raise ValueError("'input_ref' not in supplied params")

        return {'shock_id': self._get_pdb_shock_id(params['input_ref'])}

    def structure_to_pdb_file(self, params):
        """
            structure_to_pdb_file: get the file path for the given pdb object
        """
        if "input_ref" not in params:
            raise ValueError("input_ref not in supplied params")
        if "destination_dir" not in params:
            raise ValueError("destination_dir not in supplied params")

        shock_id = self._get_pdb_shock_id(params['input_ref'])
        file_path = self.dfu.shock_to_file({
            'shock_id': shock_id,
            'file_path': params['destination_dir'],
            'unpack': 'uncompress'
        })['file_path']

        return {'file_path': file_path}

    def batch_import_pdbs(self, params):
        """
            batch_import_pdbs: upload two sets of pdb files and create a
                                   KBaseStructure.ProteinStructures object
            required params:
                metadata_staging_file_path: a metafile from the user's staging area that must be a
                    subdirectory file path in staging area,
                    e.g., /data/bulk/user_name/metadata_staging_file_path
                          staging_file_subdir_path is metadata_staging_file_path
                structures_name: name of the ProteinStructures object to be generated
                workspace_name: workspace name that the protein structure(s) will be saved
            return:
                structures_ref: return ProteinStructures object reference
                report_name: name of generated report (if any)
                report_ref: report reference (if any)

            1. call _validate_batch_import_pdbs_params to validate input params
            2. call _parse_metadata to parse for model_pdb_files, exp_pdb_files and kbase_meta_data
            3. call import_model_pdb_file on each entry in model_pdb_paths, and
               call import_experiment_pdb_file on each entry in exp_pdb_paths
            4. assemble the data for a ProteinStructures and save the data object
            5. call _generate_batch_report to generate a report for batch_import_pdbs' result
        """

        (metadata_file_path, workspace_name,
         structures_name) = self._validate_batch_import_pdbs_params(params)

        if not isinstance(workspace_name, int):
            workspace_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            workspace_id = workspace_name
        params['workspace_id'] = workspace_id

        (pdb_file_paths, narrative_ids, genome_names,
         feature_ids) = self._parse_metadata_file(metadata_file_path, workspace_id)

        model_pdb_objects = list()
        exp_pdb_objects = list()
        successful_files = list()
        failed_files = list()
        protein_structures = dict()
        total_structures = 0

        pdb_params = {}
        # loop through the list of pdb_file_paths
        for pdb_file in pdb_file_paths:
            pdb_params['input_staging_file_path'] = pdb_file['file_path']
            pdb_params['input_file_path'] = None
            pdb_params['input_shock_id'] = None
            pdb_params['workspace_name'] = workspace_name
            pdb_params['structure_name'] = pdb_file['structure_name']
            pdb_params['pdb_info'] = pdb_file

            if pdb_file['is_model']:
                model_pdb_ret = self.import_model_pdb_file(pdb_params, False)
                if model_pdb_ret:
                    ret_strct_ref = model_pdb_ret['structure_obj_ref']
                    model_pdb_objects.append(ret_strct_ref)
                    model_pdb_ret['pdb_info']['pdb_struct_ref'] = ret_strct_ref
                    successful_files.append(model_pdb_ret['pdb_info'])
                    total_structures += 1
                else:
                    failed_files.append(pdb_file)
            else:
                exp_pdb_ret = self.import_experiment_pdb_file(pdb_params, False)
                if exp_pdb_ret:
                    ret_strct_ref = exp_pdb_ret['structure_obj_ref']
                    exp_pdb_objects.append(ret_strct_ref)
                    exp_pdb_ret['pdb_info']['pdb_struct_ref'] = ret_strct_ref
                    successful_files.append(exp_pdb_ret['pdb_info'])
                    total_structures += 1
                else:
                    failed_files.append(pdb_file)

        if model_pdb_objects:
            protein_structures['model_structures'] = model_pdb_objects
        if exp_pdb_objects:
            protein_structures['experimental_structures'] = exp_pdb_objects
        protein_structures['total_structures'] = total_structures
        protein_structures['description'] = (f'Created {total_structures} '
                                             f'structures in {structures_name}')
        logging.info(protein_structures)
        try:
            info = self.dfu.save_objects({
                'id': workspace_id,
                'objects': [
                    {'type': 'KBaseStructure.ProteinStructures',
                     'name': structures_name,
                     'data': protein_structures}]
            })[0]
        except WorkspaceError as e:
            logging.info(f'DFU.save_objects errored with message: {e.message} and data: {e.data}')
            raise
        else:
            structs_ref = f"{info[6]}/{info[0]}/{info[4]}"

        returnVal = {'structures_ref': structs_ref}

        report_output = self._generate_batch_report(workspace_name, structs_ref, structures_name,
                                                    successful_files, failed_files)
        returnVal.update(report_output)
        return returnVal
