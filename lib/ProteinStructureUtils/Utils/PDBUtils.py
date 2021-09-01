import hashlib
import logging
import os
import shutil
import uuid
import errno
import pandas as pd
import numpy as np

from Bio import PDB
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from installed_clients.AbstractHandleClient import AbstractHandle
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport


class PDBUtil:

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

    def _model_file_to_data(self, file_path, ref_seq='TGTGACTA'):
        """
            _model_file_to_data:
                Do the PDB conversion--parse the model pdb file for creating a pdb data object
        """
        parser = PDB.PDBParser(PERMISSIVE=1)
        ppb = PPBuilder()
        pdb1 = file_path
        structure = parser.get_structure("test", pdb1)
        model = structure[0]
        chain_no = 0
        res_no = 0
        atom_no = 0
        pp_list = []
        pp_no = 0
        for model in structure:
            for chain in model:
                chain_no += 1
        for residue in model.get_residues():
            if PDB.is_aa(residue):
                res_no += 1
            for atom in residue.get_atoms():
                atom_no += 1

        for pp in ppb.build_peptides(structure):
            pp_no += 1
            my_seq = pp.get_sequence()
            pp_list += str(my_seq)
        seq = ''.join(pp_list)

        (compound, source) = self._get_compound_source(structure)
        seq_iden = self._compute_sequence_identity(ref_seq, seq)

        data = {
            'name': structure.header.get('name', ''),
            'num_chains': chain_no,
            'num_residues': res_no,
            'num_atoms': atom_no,
            'model_id': model.get_id(),
            'identity': seq_iden,
            'exact_match': 1 if seq_iden == 1.0 else 0,
            'compound': compound,
            'source': source,
            'proteins': [{
                'id': os.path.basename(file_path),
                'sequence': seq,
                'md5': hashlib.md5(seq.encode()).hexdigest()
            }],
        }

        return data, pp_no

    def _exp_file_to_data(self, file_path, ref_seq='TGTGACTA'):
        """
            _exp_file_to_data:
                Do the PDB conversion--parse the experiment pdb file for creating a pdb data object
        """
        # TODO: Figure out how to parse the experimental pdb file for an experimental data structure
        parser = PDB.MMCIFParser()
        cif = file_path
        structure = parser.get_structure("PHA-L", cif)
        # create a dictionary that maps all mmCIF tags in an mmCIF file to their values
        # mmcif_dict = MMCIF2Dict(cif)  # TypeError: 'module' object is not callable
        # print("The mmcif dictionary:*****************************************\n\n")
        # E.g., get the solvent content from an mmCIF file:
        # sc = mmcif_dict["_exptl_crystal.density_percent_sol"]
        # print(f'_exptl_crystal.density_percent_sol: {sc}')
        # get the list of the y coordinates of all atoms
        # y_list = mmcif_dict["_atom_site.Cartn_y"]
        # print(f'_atom_site.Cartn_y: {y_list}')

        ppb = PPBuilder()
        pp_list = []
        pp_no = 0
        for pp in ppb.build_peptides(structure):
            pp_no += 1
            my_seq = pp.get_sequence()
            pp_list += str(my_seq)
        seq = ''.join(pp_list)
        struc_name = structure.header.get('name', '')
        hd = self._upload_to_shock(file_path)

        (cpd, src) = self._get_compound_source(structure)
        (num_models, model_ids) = self._get_models(structure)
        (num_chains, chain_ids) = self._get_chains(structure)
        (num_residues, residue_ids) = self._get_residues(structure)
        (num_atoms, atom_ids) = self._get_atoms(structure)
        seq_iden = self._compute_sequence_identity(ref_seq, seq)

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
            'model_id': model_ids[0],
            'identity': seq_iden,
            'exact_match': 1 if seq_iden == 1.0 else 0,
            'pdb_handle': hd,
            'mmcif_handle': hd,
            'xml_handle': hd,
            'proteins': [{
                'id': struc_name,
                'sequence': seq,
                'md5': hashlib.md5(seq.encode()).hexdigest()
            }]
        }

        return mmcif_data, pp_no

    def _compute_sequence_identity(self, seq1, seq2):
        """
            _compute_sequence_identity: Given two input sequences, do a pairwise comparison and then
                                        compute and return the matching percentage.
        """
        AA_Match = 0.0
        for k in range(len(seq1)):
            if(seq1[k] == "-" or seq2[k] == "-"):
                continue

            if(seq1[k] == seq2[k]):
                AA_Match += 0

        seq1 = seq1.replace('-', '')
        seq2 = seq2.replace('-', '')

        ID1 = AA_Match / len(seq1)
        ID2 = AA_Match / len(seq2)
        return float(np.round((ID1 + ID2) / 2.0, 6))

    def _get_atoms(self, pdb_structure):
        """
            _get_atoms: Given a pdb_structure object, parse atoms into a list and return it
        """
        atom_ids = []
        num_atoms = 0
        my_atoms = pdb_structure.get_atoms()
        for a_ele in my_atoms:
            if(a_ele):
                num_atoms += 1
                atom_ids.append(a_ele.get_id())

        return (num_atoms, atom_ids)

    def _get_residues(self, pdb_structure):
        """
            _get_residues: Given a pdb_structure object, parse residues into a list and return it
        """
        res_ids = []
        num_res = 0
        my_res = pdb_structure.get_residues()
        for r_ele in my_res:
            if(r_ele):
                num_res += 1
                res_ids.append(r_ele.get_id())

        return (num_res, res_ids)

    def _get_chains(self, pdb_structure):
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

    def _get_models(self, pdb_structure):
        """
            _get_models: Given a pdb_structure object, parse model ids into a list and return it
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
        # print(f'Compound:\n {cpd}')
        if cpd and cpd.get('1'):
            cpd_dict = cpd.get('1')

        src_dict = dict()
        src = structure.header.get('source', {})
        # print(f'Source:\n {src}')
        if src and src.get('1'):
            src_dict = src.get('1')

        return (cpd_dict, src_dict)

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
                raise ValueError(f'"{e.strerror}" occurred')
        else:
            return fh

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
        logging.info('Start uploading file to shock: {}'.format(file_path))

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

        if params.get('metadata_staging_file_path'):
            file_path = self.dfu.download_staging_file(
                {'staging_file_subdir_path': params.get('metadata_staging_file_path')}
                 ).get('copy_file_path')
        else:
            error_msg = "Must supply a 'metadata_staging_file_path'"
            raise ValueError(error_msg)

        return (params['metadata_staging_file_path'], params['workspace_name'],
                params['structures_name'])

    def _parse_metadata_file(self, metadata_staging_file_path):
        """
            _parse_metadata_file:
                From metadata_file_path, a spreadsheet file, sort out the model_pdb_file_paths,
            exp_pdb_file_paths and the kbase_meta_data

            return: lists model_pdb_file_paths and exp_pdb_file_paths and dict kbase_meta_data
        """
        exp_pdb_file_paths = []
        model_pdb_file_paths = []
        kb_meta_data = dict()

        download_staging_file_params = {
            'staging_file_subdir_path': metadata_staging_file_path
        }
        scratch_file_path = self.dfu.download_staging_file(
                        download_staging_file_params).get('copy_file_path')

        # read the data from scratch_file_path, assuming it is a .tsv file
        tsv_df = pd.read_csv(scratch_file_path, sep='\t')  # tsv_read is a Panda DataFrame object
        # print the first 10 records
        print(tsv_df.head(10))
        # TODO!!!: Handle the file parsing details, depending on the file

        if not model_pdb_file_paths and not exp_pdb_file_paths:
            error_msg = "At least one list of pdb file(s) should be provided!"
            raise ValueError(error_msg)

        return (model_pdb_file_paths, exp_pdb_file_paths, kbase_meta_data)

    def _generate_batch_report(self, structs_ref, workspace_name, structs_name,
                               successful_paths, failed_paths):
        """
            _generate_batch_report: generate summary report for upload
        """
        output_html_files = self._generate_batch_report_html(structs_name,
                                                            successful_paths, failed_paths)

        failed_files = ','.join(failed_paths)
        description = (f'Imported PDBs into a ProteinStructures object {structs_ref}, '
                       f'with files "{failed_files}" failed to load.')
        report_params = {'message': f'You uploaded a batch of PDB files into {structs_name}.',
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

    def _generate_batch_report_html(self, pdb_name, succ_pdb_paths, fail_pdb_paths):
        """
            _generate_batch_report_html: generates the HTML for the upload report
        """
        html_report = list()

        # Make report directory and copy over files
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(output_directory)
        result_file_path = os.path.join(output_directory, 'viewer.html')
        new_pdb_path = os.path.join(output_directory, os.path.basename(pdb_path))
        shutil.copy(pdb_path, new_pdb_path)

        # Fill in template HTML--TODO: Need to create a new template.html for the batch report!!!!!
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

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.token = config['KB_AUTH_TOKEN']
        self.dfu = DataFileUtil(self.callback_url)
        self.hs = AbstractHandle(config['handle-service-url'])

    def import_model_pdb_file(self, params, create_report=True):
        """
            import_model_pdb_file: upload an experiment pdb file and convert into a
                                  KBaseStructure.ModelProteinStructure object
        """

        file_path, workspace_name, pdb_name = self._validate_import_pdb_file_params(params)

        if not isinstance(workspace_name, int):
            workspace_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            workspace_id = workspace_name

        data, n_polypeptides = self._model_file_to_data(file_path)
        data['pdb_handle'] = self._upload_to_shock(file_path)
        data['user_data'] = params.get('description', '')
        logging.info(data)

        info = self.dfu.save_objects({
            'id': workspace_id,
            'objects': [
                {'type': 'KBaseStructure.ModelProteinStructure',
                 'name': pdb_name,
                 'data': data}]
        })[0]
        obj_ref = f"{info[6]}/{info[0]}/{info[4]}"

        returnVal = {'structure_obj_ref': obj_ref}

        if create_report:
            report_output = self._generate_report('import_model_pdb_file', obj_ref, workspace_name,
                                                  n_polypeptides, pdb_name, file_path)
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
        data, n_polypeptides,  = self._exp_file_to_data(file_path)
        logging.info(data)

        info = self.dfu.save_objects({
            'id': workspace_id,
            'objects': [
                {'type': 'KBaseStructure.ExperimentalProteinStructure',
                 'name': mmcif_name,
                 'data': data}]
        })[0]
        obj_ref = f"{info[6]}/{info[0]}/{info[4]}"

        returnVal = {'structure_obj_ref': obj_ref}

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

            1. call _validate_batch_import_pdbs_params to validate input params
            2. call _parse_metadata to parse for model_pdb_files, exp_pdb_files and kbase_meta_data
            3. call import_model_pdb_file on each entry in model_pdb_paths, and
               call import_experiment_pdb_file on each entry in exp_pdb_paths
            4. assemble the data for a ProteinStructures and save the data object
            5. call _generate_batch_report to generate a report for batch_import_pdbs' result

            required params: metadata_staging_file_path, a file path in the staging area
                e.g., /data/bulk/user_name/metadata_file_path
                      staging_file_subdir_path is metadata_staging_file_path
            return:
            structures_ref: return ProteinStructures object reference
            report_name: name of generated report (if any)
            report_ref: report reference (if any)
        """

        (metadata_staging_file_path, workspace_name,
         structures_name) = self._validate_batch_import_pdbs_params(params)

        if not isinstance(workspace_name, int):
            workspace_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            workspace_id = workspace_name

        (model_pdb_files, exp_pdb_files,
         kbase_meta_data) = self._parse_metadata_file(metadata_staging_file_path)

        model_pdb_objects = list()
        exp_pdb_objects = list()
        successful_files = list()
        failed_files = list()
        protein_structures = dict()
        total_structures = 0

        # loop through the lists model_pdb_file_paths and exp_pdb_file_paths
        if model_pdb_files:
            for model_file in model_pdb_files:
                model_pdb_ref = self.import_model_pdb_file(params, False)

                if model_pdb_ref:
                    model_pdb_objects.append(model_pdb_ref)
                    successful_files.append(model_file)
                    total_structures += 1
                else:
                    failed_files.append(model_file)

        if exp_pdb_files:
            for exp_file in exp_pdb_files:
                exp_pdb_ref = self.import_experiment_pdb_file(params, False)

                if exp_pdb_ref:
                    exp_pdb_objects.append(exp_pdb_ref)
                    successful_files.append(exp_file)
                    total_structures += 1
                else:
                    failed_files.append(exp_file)

        if model_pdb_objects:
            protein_structures['model_structures'] = model_pdb_objects
        if exp_pdb_objects:
            protein_structures['experimental_structures'] = exp_pdb_objects
        protein_structures['total_structures'] = total_structures
        protein_structures['description'] = (f'Created {total_structures} '
                                             f'structures in {structures_name}')

        logging.info(protein_structures)

        info = self.dfu.save_objects({
            'id': workspace_id,
            'objects': [
                {'type': 'KBaseStructure.ProteinStructures',
                 'name': structures_name,
                 'data': protein_structures}]
        })[0]
        obj_ref = f"{info[6]}/{info[0]}/{info[4]}"

        returnVal = {'structures_obj_ref': obj_ref}

        report_output = self._generate_batch_report(obj_ref, workspace_name, structures_name,
                                                    successful_files, failed_files)

        returnVal.update(report_output)

        return returnVal
