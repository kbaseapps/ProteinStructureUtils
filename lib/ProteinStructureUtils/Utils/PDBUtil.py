import hashlib
import logging
import os
import sys
import re
import gzip
import string
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

    # “Expect Value” threshold to restrict which alignments will be significant
    E_VALUE_THRESH = 1e-20

    # BLAST sequence identity threshold to determine which pdb structures will be
    # matched to a KBase genome/feature
    B_IDENTITY_THRESH = 0.6

    def _validate_import_file_params(self, params):
        """
            _validate_import_file_params:
                validates input params to import_pdb_file and import_mmdif_file
        """
        # check for required parameters
        for p in ['structure_name', 'workspace_name']:
            if p not in params:
                raise ValueError(f'Parameter "{p}" is required, but missing!')

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

    def _get_pdb_id(self, pdb_file):
        pdb_id = os.path.basename(pdb_file)
        for ext in [r'\.gz$', r'\.pdb$', r'\.ent$', r'\.cif$', r'^pdb']:
            pdb_id = re.sub(ext, '', pdb_id)
        if pdb_id.startswith('ent') and len(pdb_id) > 4:
            pdb_id = pdb_id[3:]

        return pdb_id

    def _get_pdb_structure(self, pdb_file, pdb_id=None, quiet=True):
        """Set QUIET to False to output warnings like incomplete chains etc."""

        if not pdb_id:
            pdb_id = self._get_pdb_id(pdb_file)
        parser = PDB.PDBParser(PERMISSIVE=1, get_header=True, QUIET=quiet)
        if pdb_file.endswith('.gz'):
            with gzip.open(pdb_file, 'rt') as ifh:
                structure = parser.get_structure(pdb_id, ifh)
        else:
            structure = parser.get_structure(pdb_id, pdb_file)

        # Rename empty chains (i.e. chain.id == ' ')
        model = structure[0]
        chain_ids = {chain.id for chain in model.child_list}
        for chain in model.child_list:
            if chain.id in [' ', 'Z']:
                chain_ids.remove(chain.id)
                chain.id = next(c for c in string.ascii_uppercase if c not in chain_ids)
                chain_ids.add(chain.id)
        model.child_dict = {chain.id: chain for chain in model.child_list}

        return structure

    def _pdb_file_to_data(self, file_path, params):
        """
            _pdb_file_to_data:
                Do the PDB conversion--use PDB.PDBParser to parse the pdb file for
                                       creating a pdb data object
        """
        logging.info(f'Parsing pdb file {file_path} to a ProteinData object with params: {params}')

        pp_no = 0
        pdb_data = {}
        try:
            structure = self._get_pdb_structure(file_path, pdb_id=None, quiet=False)
        except (RuntimeError, TypeError, KeyError, ValueError) as e:
            logging.info(f'PDBParser errored with message: {e.message}')
            raise
        else:
            ppb = PPBuilder()
            for pp in ppb.build_peptides(structure):
                pp_no += 1

            # logging.info(f'Getting pdb structure data for {structure}!')
            (compound, source) = self._get_compound_source(structure)
            (num_models, model_ids) = self._get_models_from_structure(structure)
            (num_chains, chain_ids) = self._get_chains_from_structure(structure)
            (num_residues, residue_ids) = self._get_residues_from_structure(structure)
            (num_atoms, atom_ids) = self._get_atoms_from_structure(structure)
            model = structure[0]
            protein_data = self._get_proteins_by_structure(structure, model.get_id(), file_path)
            (protein_data, params) = self._match_features(params, protein_data)

            pdb_info = params.get('pdb_info', None)
            if pdb_info and pdb_info.get('sequence_identities', None):
                pdb_data = {
                    'name': structure.header.get('name', ''),
                    'num_chains': num_chains,
                    'num_residues': num_residues,
                    'num_atoms': num_atoms,
                    'compound': compound,
                    'source': source,
                    'proteins': protein_data
                }
            else:
                logging.info(f'PDB file {file_path} failed to match KBase genome/features!')
                pdb_data = {}
        finally:
            return pdb_data, pp_no, params

    def _mmcif_file_to_data(self, file_path, params):
        """
            _mmcif_file_to_data:
                Do the PDB conversion--use PDB.MMCIFParser to parse the experiment pdb file for
                                       creating a pdb data object
        """
        logging.info(f'Parsing pdb file {file_path} to a ProteinData object with params: {params}')

        parser = PDB.MMCIFParser()
        cif = file_path
        pp_no = 0
        mmcif_data = None

        try:
            structure = parser.get_structure("PHA-L", cif)
        except (RuntimeError, TypeError, KeyError, ValueError) as e:
            logging.info(f'MMCIFParser errored with message: {e.message}')
            raise
        else:
            ppb = PPBuilder()
            for pp in ppb.build_peptides(structure):
                pp_no += 1

            struc_name = structure.header.get('name', '')
            hd = self._upload_to_shock(file_path)

            # logging.info(f'Getting pdb structure data for {structure}!')
            (cpd, src) = self._get_compound_source(structure)
            (num_models, model_ids) = self._get_models_from_structure(structure)
            (num_chains, chain_ids) = self._get_chains_from_structure(structure)
            (num_residues, residue_ids) = self._get_residues_from_structure(structure)
            (num_atoms, atom_ids) = self._get_atoms_from_structure(structure)
            protein_data = self._get_proteins_by_structure(structure, model_ids[0], file_path)
            (protein_data, params) = self._match_features(params, protein_data)

            pdb_info = params.get('pdb_info', None)
            if pdb_info and pdb_info.get('sequence_identities', None):
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
            else:
                mmcif_data = {}
                logging.info(f'PDB file {file_path} failed to match KBase genome/features!')
        finally:
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
            kb_feature_type = ''
            kb_feature_seq = ''
            genome_name = pdb_info['genome_name']
            narr_id = pdb_info['narrative_id']
            feature_id = pdb_info['feature_id']

            logging.info(f"Looking up for feature {feature_id} in genome {genome_name}'s features")
            # 1. Get the genome's features and reference
            (gn_ref, kb_genome_features) = self._get_genome_ref_features(narr_id, genome_name)
            if not gn_ref:
                logging.info(f"Given genome {genome_name} does not exist in workspace {narr_id}!")
                return protein_data, params

            pdb_info['genome_ref'] = gn_ref
            # 2. Match the genome features with the specified feature_id to obtain feature sequence
            for feat in kb_genome_features:
                if feat['id'] == feature_id:
                    logging.info(f'Found genome feature match for {feature_id}')
                    kb_feature_type = self._get_feature_type(feat)
                    kb_feature_seq = feat.get('protein_translation', '')
                    break

            pdb_info['feature_type'] = kb_feature_type

            # 3. Call self._compute_sequence_identity with the feature sequence and the the pdb
            # proteins' translations to to get the seq_identity and exact_match
            if kb_feature_seq:
                logging.info(f"Finding seq_identity and exact_match for feature {feature_id}"
                             f" in genome {genome_name}'s features...")
                pdb_chain_ids = []
                pdb_model_ids = []
                pdb_seq_idens = []
                pdb_exact_matches = []
                for prot in protein_data:
                    seq_idens, seq_mats = self._compute_sequence_identity(kb_feature_seq,
                                                                          prot.get('sequence', ''))
                    if seq_idens:
                        seq_idens.sort()
                        max_iden = seq_idens.pop()
                        if max_iden >= self.B_IDENTITY_THRESH:  # get the good matches
                            prot['seq_identity'] = max_iden
                            prot['exact_match'] = 1 if max_iden > 0.99 else 0
                            prot['genome_ref'] = gn_ref
                            prot['feature_id'] = feature_id
                            prot['feature_type'] = kb_feature_type
                            pdb_chain_ids.append(f"Model {prot['model_id']+1}.Chain {prot['chain_id']}")
                            pdb_model_ids.append(str(prot['model_id']))
                            pdb_seq_idens.append(f"{round(prot['seq_identity']*100, 2)}%")
                            pdb_exact_matches.append(str(prot['exact_match']))
                if pdb_seq_idens:
                    pdb_info['sequence_identities'] = ','.join(pdb_seq_idens)
                if pdb_chain_ids:
                    pdb_info['chain_ids'] = ','.join(pdb_chain_ids)
                if pdb_model_ids:
                    pdb_info['model_ids'] = ','.join(pdb_model_ids)
                if pdb_exact_matches:
                    pdb_info['exact_matches'] = ','.join(pdb_exact_matches)
            else:
                logging.info(f'Found NO feature in genome that matches with {feature_id}')
        else:
            logging.info('NO KBase genome/feature object info were given for uploading')

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
            _get_genome_ref_features: Get the genome reference and features for genome_name
        """
        genome_ref = ''
        genome_features = []
        (genome_info, genome_data) = self._get_object_info_data(narr_id, genome_name)
        if genome_info and genome_data:
            genome_ref = '/'.join([str(narr_id), str(genome_info[0]), str(genome_info[4])])
            genome_features = genome_data['features']

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
            try:
                obj_data_res = self.ws_client.get_objects2(
                    {'objects': [{'wsid': narr_id, 'name': obj_name}]})['data'][0]
                obj_info = obj_data_res['info']
                obj_data = obj_data_res['data']
            except:
                logging.info(f'No object with name {obj_name} exists in workspace {narr_id}')
                logging.info(f'Unexpected error occurred while getting object for {obj_name}')
                pass

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

    def _dfu_get_objects(self, obj_ref):
        """
            _dfu_get_objects: call dfu.get_objects to return object data and info
        """
        obj = self.dfu.get_objects({"object_refs": [obj_ref]})['data'][0]
        return obj['data'], obj['info']

    def _get_pdb_shock_id(self, obj_ref):
        """
            _get_pdb_shock_id: Return the shock id for the PDB file
        """
        obj_data, obj_info = self._dfu_get_objects(obj_ref)
        return self.hs.hids_to_handles([obj_data['pdb_handle']])[0]['id']

    def _get_struct_shock_id(self, struct_ref, file_ext):
        """
            _get_struct_shock_id: Return the shock id for the structure file of the given extension
        """
        obj_data, obj_info = self._dfu_get_objects(struct_ref)
        shock_id = ''
        if file_ext in ('pdb', 'mmcif', 'xml'):
            handle_type = '_'.join([file_ext, 'handle'])
            shock_id = self.hs.hids_to_handles([obj_data[handle_type]])[0]['id']

        return shock_id

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

    def _validate_batch_import_params(self, params):
        """
            _validate_batch_import_params:
                validates params passed to batch_import_pdbs method
        """
        # check for required parameters
        for p in ['structures_name', 'workspace_name', 'metadata_staging_file_path']:
            if p not in params:
                raise ValueError(f'Parameter "{p}" is required, but missing')

        # metadata_staging_file_path must be from the staging area--must have the staging dir prefix
        if params.get('metadata_staging_file_path', None):
            staging_file_path = self.dfu.download_staging_file(
                {'staging_file_subdir_path': params.get('metadata_staging_file_path')}
                 ).get('copy_file_path')
            return (staging_file_path, params['workspace_name'], params['structures_name'])
        else:
            raise ValueError('Must supply a "metadata_staging_file_path"')

    def _read_file_by_type(self, file_path):
        """
            _read_file_by_type: read the file given by file_path depending on its type,
                               return a DataFrame object
        """
        logging.info(f'Reading input from file: {file_path}...')

        if not self._validate_file(file_path):
            raise ValueError('Input file is invalid or not found!')

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
                error_msg = 'Invalid input file type, only "csv/tsv/xlsx" are accepted!'
                raise ValueError(error_msg)
            # strip off the leading and trailing whitespaces of the column names
            df.columns = df.columns.str.strip()
        except (RuntimeError, TypeError, KeyError, ValueError, WorkspaceError) as e:
            logging.info(
                f'Reading file {file_path} errored with message: {e.message} and data: {e.data}')
            raise
        return df

    def _parse_metadata_file(self, metadata_file_path, ws_id):
        """
            _parse_metadata_file:
                From metadata_file_path, a spreadsheet file, sort out the pdb_file_paths
                and the kbase_meta_data

            return: a list pdb_file_paths containing objects of structure:
                    {'file_path': pdb_fn,
                     'file_extension': ext,
                     'structure_name': struct_name,
                     'narrative_id': narr_id,
                     'genome_name': obj_name,
                     'feature_id': feat_id,
                     'is_model': 1 if 'y' in is_model or 'Y' in is_model else 0},
                    and kbase object lists--narrative_ids, genome_names, feature_ids
        """
        logging.info(f'Parsing metadata from input file {metadata_file_path}...')

        required_columns = ['Narrative ID', 'Object name (Genome AMA feature set)', 'Feature ID',
                            'PDB filename', 'Is model']
        # Only extensions ‘.cif’ or ‘.pdb’ are valid
        accepted_extensions = ['.pdb', '.cif']

        pdb_file_paths = list()
        narrative_ids = list()
        genome_names = list()
        feature_ids = list()

        # df_meta_data is a Panda DataFrame object
        df_meta_data = self._read_file_by_type(metadata_file_path)
        df_columns = df_meta_data.columns
        df_col_list = df_columns.values.tolist()

        # check if required columns are read in correctly
        for col in required_columns:
            if col not in df_col_list:
                raise ValueError(f'Required column "{col}" is missing!')

        for i in range(len(df_meta_data[df_columns[0]])):
            narr_id = int(df_meta_data[df_columns[0]][i])
            if not pd.isna(narr_id):
                narrative_ids.append(narr_id)
            else:
                raise ValueError(f'Please fill all the rows in column: {required_columns[0]}!')

            obj_name = df_meta_data[df_columns[1]][i]
            if not pd.isna(obj_name):
                genome_names.append(obj_name)
            else:
                raise ValueError('Please fill all the rows in column: Object name!')

            feat_id = df_meta_data[df_columns[2]][i]
            if not pd.isna(feat_id):
                feature_ids.append(feat_id)
            else:
                raise ValueError(f'Please fill all the rows in column: {required_columns[2]}!')

            pdb_fn = df_meta_data[df_columns[3]][i]  # pdb_fn does not have staging dir prefix
            if pd.isna(pdb_fn):
                raise ValueError(f'Please fill all the rows in column: {required_columns[3]}!')

            (struct_name, ext) = os.path.splitext(os.path.basename(pdb_fn))
            if ext not in accepted_extensions:
                # raise ValueError('Only files with extensions ".cif" or ".pdb" are accepted.')
                print('Only files with extensions ".cif" or ".pdb" are accepted.')

            is_model = df_meta_data[df_columns[4]][i]
            if not pd.isna(is_model):
                pdb_file_paths.append(
                    {'file_path': pdb_fn,
                     'file_extension': ext,
                     'structure_name': struct_name,
                     'narrative_id': narr_id,
                     'genome_name': obj_name,
                     'feature_id': feat_id,
                     'is_model': 1 if 'y' in is_model or 'Y' in is_model else 0}
                )
            else:
                raise ValueError(f'Please fill all the rows in column: {required_columns[4]}!')

        if not pdb_file_paths:
            raise ValueError('No PDB file info is provided!')

        return (pdb_file_paths, narrative_ids, genome_names, feature_ids)

    def _generate_batch_report(self, workspace_name, structs_ref, structs_name,
                               pdb_infos, failed_pdbs):
        """
            _generate_batch_report: generate summary report for upload
        """

        output_html_files = self._generate_batch_report_html(pdb_infos)

        description = (f'Imported PDBs into a ProteinStructures object "{structs_ref}", '
                       f'named "{structs_name}".')

        if failed_pdbs:
            failed_files = ','.join(failed_pdbs)
            description += f' These files "{failed_files}" failed to load.'

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

    def _config_viewer(self, viewer_nm, div_id):
        """
            _config_viewer: write the mol* viewer configurations
        """
        return (f'<script type="text/javascript">\n'
                f'let {viewer_nm} = new molstar.Viewer("{div_id}", {{\n'
                f'layoutIsExpanded: false,\n'
                f'layoutShowControls: true,\n'
                f'layoutShowRemoteState: false,\n'
                f'layoutShowSequence: true,\n'
                f'layoutShowLog: false,\n'
                f'layoutShowLeftPanel: true,\n'
                f'viewportShowExpand: false,\n'
                f'viewportShowSelectionMode: false,\n'
                f'viewportShowAnimation: true,\n'
                f'collapseLeftPanel: false,\n'
                f'}});\n')

    def _write_viewer_content_single(self, output_dir, succ_pdb_infos):
        """
            _write_viewer_content_single: write the mol* viewer html content to fill in the
                                        subtabcontent and replace the string <!--replace subtabs-->
                                        in the templatefile 'batch_pdb_template.html'
        """
        viewer_content = ''
        pdb_index = 0

        for succ_pdb in succ_pdb_infos:
            file_path = succ_pdb['file_path']
            file_ext = succ_pdb['file_extension'][1:]
            if file_ext == 'cif':
                file_ext = 'mmcif'
            pdb_file_path = succ_pdb['scratch_path']  # this is the scratch path for this pdb file
            base_filename = os.path.basename(file_path)
            new_pdb_path = os.path.join(output_dir, base_filename)
            shutil.copy(pdb_file_path, new_pdb_path)

            struct_nm = succ_pdb['structure_name'].upper()
            viewer_name = 'viewer' + str(pdb_index + 1)
            div_id = 'struct' + str(pdb_index + 1)

            sub_div = (f'<div id="{struct_nm}" class="subtabcontent">\n'
                       f'<h2>{struct_nm}</h2>\n'
                       f'<div id="{div_id}" class="app"></div>\n')

            script_content = self._config_viewer(viewer_name, div_id)
            script_content += (f'{viewer_name}.loadStructureFromUrl("./{base_filename}", '
                               f'"{file_ext}", false, {{ representationParams: '
                               f'{{ theme: {{ globalName: "operator-name" }} }} }});')
            script_content += '\n</script>'

            sub_div += script_content
            sub_div += '\n</div>\n'
            viewer_content += sub_div
            pdb_index += 1

        return viewer_content

    def _write_viewer_content_multi(self, output_dir, succ_pdb_infos):
        """
            _write_viewer_content_multi: write the mol* viewer html content to fill in the
                                         subtabcontent and replace the string <!--replace subtabs-->
                                         in the templatefile 'batch_pdb_template.html'
        """
        viewer_nm = 'viewer_all'
        div_id = 'app_all'
        pre_loads = ''
        viewer_tabs = ('<div class="tab">'
                       '<button id="AllStructures_sub" class="subtablinks" '
                       'onclick="openSubTab(event, this)">ALL STRUCTURES</button>')
        viewer_content = ('<div id="AllStructures" class="subtabcontent">'
                          '<h2>Uploaded Structure(s)</h2><div id="app_all" class="app"></div>')

        script_content = self._config_viewer(viewer_nm, div_id)

        for succ_pdb in succ_pdb_infos:
            struct_nm = succ_pdb['structure_name'].upper()
            file_path = succ_pdb['file_path']
            file_ext = succ_pdb['file_extension'][1:]
            if file_ext == 'cif':
                file_ext = 'mmcif'
            pdb_file_path = succ_pdb['scratch_path']  # this is the scratch path for this pdb file
            new_pdb_path = os.path.join(output_dir, os.path.basename(file_path))
            shutil.copy(pdb_file_path, new_pdb_path)

            viewer_tabs += (f'\n<button id="{struct_nm}_sub" '
                            f'class="subtablinks" onclick="openSubTab(event, this)">'
                            f'{struct_nm}</button>')

            pre_loads += (f'\n{viewer_nm}.loadStructureFromUrl("./{os.path.basename(file_path)}", '
                          f'"{file_ext}", false, {{ representationParams: '
                          f'{{ theme: {{ globalName: "operator-name" }} }} }});')

        viewer_tabs += '\n</div>'

        # insert the structure file for preloading
        script_content += pre_loads
        script_content += '\n</script>'
        viewer_content += script_content
        viewer_content += '\n</div>\n'

        return viewer_tabs, viewer_content

    def _write_structure_info(self, output_dir, succ_pdb_infos):
        """
            _write_structure_info: write the batch uploaded structure info to replace the string
                                   '<!--replace uploaded pdbs tbody-->' in the tboday tag of the
                                   jQuery DataTable in the template file 'batch_pdb_template.html'
        """

        tbody_html = ''
        srv_domain = urlparse(self.shock_url).netloc  # parse url to get the domain portion
        srv_base_url = f'https://{srv_domain}'
        logging.info(f'Get the url for building the anchors: {srv_base_url}')

        for succ_pdb in succ_pdb_infos:
            tbody_html += '<tr>'
            file_path = succ_pdb['file_path']
            pdb_file_path = succ_pdb['scratch_path']  # this is the scratch path for this pdb file
            new_pdb_path = os.path.join(output_dir, os.path.basename(file_path))
            shutil.copy(pdb_file_path, new_pdb_path)

            struct_nm = succ_pdb['structure_name'].upper()
            genome_name = succ_pdb['genome_name']
            genome_ref = succ_pdb['genome_ref']
            feat_id = succ_pdb['feature_id']

            pdb_chains = []
            seq_idens = []
            if succ_pdb.get('chain_ids', None):
                pdb_chains = succ_pdb['chain_ids']
            if succ_pdb.get('sequence_identities', None):
                seq_idens = succ_pdb['sequence_identities']

            tbody_html += (f'\n<td><div class="subtablinks" '
                           f'onclick="openSubTab(event, this, false)" '
                           f'style="cursor:pointer;color:blue;text-decoration:underline;" '
                           f'title="Click to see in mol*">{struct_nm}</div></td>')
            tbody_html += (f'\n<td><a href="{srv_base_url}/#dataview/{genome_ref}"'
                           f' target="_blank">{genome_name}</a></td><td>{feat_id}</td>')
            tbody_html += f'\n<td>{pdb_chains} </td>'
            tbody_html += f'\n<td>{seq_idens}</td>'
            tbody_html += '\n</tr>'

        return tbody_html

    def _generate_batch_report_html(self, succ_pdb_infos):
        """
            _generate_batch_report_html: generates the HTML for the upload report
        """
        html_report = list()

        # Make report directory and copy over uploaded pdb files
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(output_directory)

        pdb_html = self._write_structure_info(output_directory, succ_pdb_infos)
        single_viewer = self._write_viewer_content_single(output_directory, succ_pdb_infos)
        viewer_tabs, multi_viewer = self._write_viewer_content_multi(output_directory,
                                                                     succ_pdb_infos)

        dir_name = os.path.dirname(__file__)
        report_template_file = os.path.join(dir_name, 'templates', 'batch_pdb_template.html')
        report_html = os.path.join(dir_name, 'batch_pdb_viewer.html')

        with open(report_html, 'w') as report_html_pt:
            with open(report_template_file, 'r') as report_template_pt:
                # Fetch & fill in detailed info into template HTML
                batch_html_report = report_template_pt.read()\
                    .replace('<!--replace uploaded pdbs tbody-->', pdb_html)
                batch_html_report = batch_html_report\
                    .replace('<!--replace StructureViewer subtabs-->', viewer_tabs)
                batch_html_report = batch_html_report\
                    .replace('<!--replace subtab content multi-->', multi_viewer)
                batch_html_report = batch_html_report\
                    .replace('<!--replace subtab content single-->', single_viewer)
                report_html_pt.write(batch_html_report)

        molstar_js_file = os.path.join(dir_name, 'templates', 'molstar.js')
        molstar_css_file = os.path.join(dir_name, 'templates', 'molstar.css')
        molstar_ico_file = os.path.join(dir_name, 'templates', 'favicon.ico')
        shutil.copy(molstar_js_file, os.path.join(output_directory, 'molstar.js'))
        shutil.copy(molstar_css_file, os.path.join(output_directory, 'molstar.css'))
        shutil.copy(molstar_ico_file, os.path.join(output_directory, 'favicon.ico'))

        batch_html_report_path = os.path.join(output_directory, 'batch_pdb_report.html')
        shutil.copy(report_html, batch_html_report_path)
        logging.info(f'Full batch_report has been written to {batch_html_report_path}')

        html_report.append({'path': output_directory,
                            'name': os.path.basename(batch_html_report_path),
                            'description': 'HTML report for PDB upload'})

        return html_report

    def _export_pdb(self, params):
        """
            _export_pdb: return the shock_id of the uploaded pdb object
        """
        if "input_ref" not in params:
            raise ValueError('"input_ref" not in supplied params')

        return {'shock_id': self._get_pdb_shock_id(params['input_ref'])}

    def _structure_to_pdb_file(self, params):
        """
            _structure_to_pdb_file: get the file path for the given pdb object
        """
        if "input_ref" not in params:
            raise ValueError('input_ref not in supplied params')
        if "destination_dir" not in params:
            raise ValueError('destination_dir not in supplied params')

        shock_id = self._get_pdb_shock_id(params['input_ref'])
        file_path = self.dfu.shock_to_file({
            'shock_id': shock_id,
            'file_path': params['destination_dir'],
            'unpack': 'uncompress'
        })['file_path']

        return {'file_path': file_path}

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.token = config['KB_AUTH_TOKEN']
        self.user_id = config['USER_ID']
        self.dfu = DataFileUtil(self.callback_url)
        self.hs = AbstractHandle(config['handle-service-url'])
        self.ws_client = Workspace(config['workspace-url'])
        self.shock_url = config['shock-url']

    def import_pdb_file(self, params, create_report=False):
        """
            import_pdb_file: upload a pdb file and convert into a
                            KBaseStructure.ProteinStructure object
        """
        logging.info(f'import_pdb_file to a pdb data structure with params: {params}')

        # file_path is the pdb file's working area path (after dfu.download_staging_file call)
        file_path, workspace_name, pdb_name = self._validate_import_file_params(params)

        (data, n_polypeptides, params) = self._pdb_file_to_data(file_path, params)
        if not data:
            logging.info(f'PDB file {file_path} import with "import_pdb_file" failed!')
            return {}, {}

        data['pdb_handle'] = self._upload_to_shock(file_path)
        data['user_data'] = params.get('description', '')
        data['is_model'] = params.get('is_model', 0)
        # logging.info(f'Protein structure data from a .pdb file:{data}')

        pdb_info = params.get('pdb_info', None)
        if pdb_info:
            pdb_info['scratch_path'] = file_path

        return data, pdb_info

    def import_mmcif_file(self, params, create_report=False):
        """
            import_mmcif_file: upload an mmcif file and convert into a
                              KBaseStructure.ProteinStructure object
        """
        logging.info(f'import_mmcif_file to a pdb structure with params: {params}')

        # file_path is the pdb file's working area path (after dfu.download_staging_file call)
        file_path, workspace_name, mmcif_name = self._validate_import_file_params(params)

        # Parse the experimental pdb file for an experimental data structure
        (data, n_polypeptides, params) = self._mmcif_file_to_data(file_path, params)
        if not data:
            logging.info(f'Import {file_path} with "_mmcif_file_to_data" failed!')
            return {}, {}

        data['mmcif_handle'] = self._upload_to_shock(file_path)
        data['user_data'] = params.get('description', '')
        data['is_model'] = params.get('is_model', 0)
        # logging.info(f'Protein structure data from a .cif file:{data}')

        pdb_info = params.get('pdb_info', None)
        if pdb_info:
            pdb_info['scratch_path'] = file_path
        return data, pdb_info

    def export_pdb_structures(self, params):
        """
            export_structure_handles: return the handles of the ProteinStructures object
        """
        if 'input_ref' not in params:
            raise ValueError('Variable "input_ref" is required!')

        handles = []
        objData, objInfo = self._dfu_get_objects(params['input_ref'])
        for pro_str in objData['protein_structures']:
            if pro_str.get('pdb_handle', None):
                handles.append(pro_str['pdb_handle'])
                continue
            elif pro_str.get('mmcif_handle', None):
                handles.append(pro_str['mmcif_handle'])
                continue
            elif pro_str.get('xml_handle', None):
                handles.append(pro_str['xml_handle'])
                continue

        return {'shock_ids': handles}

    def export_protein_structures(self, params):
        """
            export_pdb_structures: return the shock_ids of the ProteinStructures object
        """
        if 'input_ref' not in params:
            raise ValueError('Variable "input_ref" is required!')

        structsData, structsInfo = self._dfu_get_objects(params['input_ref'])

        export_package_dir = os.path.join(self.scratch, "output")
        if not os.path.isdir(export_package_dir):
            os.mkdir(export_package_dir)

        # TODO!!!!
        # output_file_format = params.get('file_format', 'pdb')
        # prot_structs_name = structsInfo[1]
        # output_file = os.path.join(export_package_dir, '_'.join(prot_structs_name.split()) + ".csv")
        # self._prot_structures_to_output(structsData, output_file, self.token, output_file_format)

        # package it up
        package_details = self.dfu.package_for_download({
            'file_path': export_package_dir,
            'ws_refs': [params['input_ref']]
        })

        return {
            'shock_id': package_details['shock_id'],
            'result_dir': export_package_dir
        }

    def saveStructures_createReport(self, structures_name, workspace_id, workspace_name,
                                    protein_structures, pdb_infos, failed_files):
        """
            saveStructures_createReport: With given inputs, save the ProteinStructures object
                                         create a report and return the final results
        """
        returnVal = {}
        try:
            info = self.dfu.save_objects({
                'id': workspace_id,
                'objects': [
                    {'type': 'KBaseStructure.ProteinStructures',
                     'name': structures_name,
                     'data': protein_structures}]
            })[0]
        except (RuntimeError, TypeError, KeyError, ValueError, WorkspaceError) as e:
            err_msg = f'DFU.save_objects errored with message: {e.message} and data: {e.data}'
            raise ValueError(err_msg)
        else:
            structs_ref = f"{info[6]}/{info[0]}/{info[4]}"
            returnVal = {'structures_ref': structs_ref}
            report_output = self._generate_batch_report(
                        workspace_name, structs_ref, structures_name, pdb_infos, failed_files)
            returnVal.update(report_output)
            logging.info(f'ProteinStructures data structure saved as:\n{structs_ref}')
        finally:
            return returnVal

    def batch_import_pdbs(self, params):
        """
            batch_import_pdbs: upload a list of pdb files and create a
                                   KBaseStructure.ProteinStructures object
            required params:
                metadata_staging_file_path: a metafile from the user's staging area that must be a
                    subdirectory file path in staging area,
                    e.g., /data/bulk/user_name/staging_file_subdir_path
                          where staging_file_subdir_path is metadata_staging_file_path
                structures_name: name of the ProteinStructures object to be generated
                workspace_name: workspace name that the protein structure(s) will be saved
            return:
                structures_ref: return ProteinStructures object reference
                report_name: name of generated report (if any)
                report_ref: report reference (if any)

            1. call _validate_batch_import_params to validate input params
            2. call _parse_metadata to parse for model_pdb_files, exp_pdb_files and kbase_meta_data
            3. call import_pdb_file on each entry in pdb_paths, and/or
               call import_mmcif_file on each entry in pdb_paths
            4. assemble the data for a ProteinStructures for saving the data object
            5. call saveStructures_createReport to save the object and
               generate a report for batch_import_pdbs' result
        """
        (metadata_file_path, workspace_name,
         structures_name) = self._validate_batch_import_params(params)

        if not isinstance(workspace_name, int):
            workspace_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            workspace_id = workspace_name
        params['workspace_id'] = workspace_id

        (pdb_file_paths, narrative_ids, genome_names,
         feature_ids) = self._parse_metadata_file(metadata_file_path, workspace_id)

        pdb_objects = list()
        pdb_infos = list()
        successful_files = list()
        failed_files = list()
        protein_structures = dict()

        # loop through the list of pdb_file_paths
        for pdb in pdb_file_paths:
            pdb_params = {}
            pdb_params['pdb_info'] = pdb
            pdb_params['input_staging_file_path'] = pdb['file_path']
            pdb_params['input_file_path'] = None
            pdb_params['input_shock_id'] = None
            pdb_params['workspace_name'] = workspace_name
            pdb_params['structure_name'] = pdb['structure_name']
            pdb_params['is_model'] = pdb['is_model']

            if pdb['file_extension'] == '.pdb':
                pdb_data, pdb_info = self.import_pdb_file(pdb_params)
                if pdb_data:
                    pdb_objects.append(pdb_data)
                    pdb_infos.append(pdb_info)
                    successful_files.append(pdb['file_path'])
                else:
                    failed_files.append(pdb['file_path'])
            elif pdb['file_extension'] == '.cif':
                cif_data, pdb_info = self.import_mmcif_file(pdb_params)
                if cif_data:
                    pdb_objects.append(cif_data)
                    pdb_infos.append(pdb_info)
                    successful_files.append(pdb['file_path'])
                else:
                    failed_files.append(pdb['file_path'])

        if not pdb_objects:
            logging.info("No pdb structure was created/saved!")
            return {}

        total_structures = len(pdb_objects)
        protein_structures['protein_structures'] = pdb_objects
        protein_structures['total_structures'] = total_structures
        protein_structures['description'] = (f'Created {total_structures} '
                                             f'structures in {structures_name}')

        return self.saveStructures_createReport(structures_name, workspace_id, workspace_name,
                                                protein_structures, pdb_infos, failed_files)
