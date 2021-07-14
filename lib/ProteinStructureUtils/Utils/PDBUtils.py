import hashlib
import logging
import os
import shutil
import uuid

from Bio import PDB
from Bio.PDB.Polypeptide import PPBuilder

from installed_clients.AbstractHandleClient import AbstractHandle
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport


class PDBUtil:

    def _validate_import_pdb_file_params(self, params):
        """
        _validate_import_pdb_file_params:
            validates params passed to import_model_pdb_file and import_experiment_pdb_file methods
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

    def _model_file_to_data(self, file_path):
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

        data = {
            'name': os.path.basename(file_path),
            'num_chains': chain_no,
            'num_residues': res_no,
            'num_atoms': atom_no,
            'protein': {
                'id': os.path.basename(file_path),
                'sequence': seq,
                'md5': hashlib.md5(seq.encode()).hexdigest()
            },
        }

        return data, pp_no

    def _exp_file_to_data(self, file_path):
        """
        _exp_file_to_data:
            Do the PDB conversion--parse the experiment pdb file for creating a pdb data object

        A KBaseStructure.ExperimentalProteinStructure object is defined to have a structure like:
          typedef structure {
            /*Experimental header*/
            string rcsb_id;
            string name;
            string deposition_date;
            string head;
            string release_date;
            string structure_method;
            float resolution;
            string author;
            list<string> compound;
            list<string> source;

            /*Structure metadata*/
            int num_models;
            int num_chains;
            int num_residues;
            int num_atoms;
            int num_het_atoms;
            int num_water_atoms;
            int num_disordered_atoms;
            int num_disordered_residues;

            /*Protein links*/
            list<ProteinData> proteins;

            /*File links*/
            handle_ref pdb_handle;
            handle_ref mmcif_handle;
            handle_ref xml_handle;
        } ExperimentalProteinStructure;
        """
        # TODO: Figure out how to parse the experimental pdb file for an experimental data structure
        parser = PDB.PDBParser(PERMISSIVE=1)
        ppb = PPBuilder()
        pdb2 = file_path
        structure = parser.get_structure("test", pdb2)
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

        data = {
            'name': os.path.basename(file_path),
            'num_chains': chain_no,
            'num_residues': res_no,
            'num_atoms': atom_no,
            'protein': {
                'id': os.path.basename(file_path),
                'sequence': seq,
                'md5': hashlib.md5(seq.encode()).hexdigest()
            },
        }

        return data, pp_no

    def _get_pdb_shock_id(self, obj_ref):
        """Return the shock id for the PDB file"""
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
            _generate_report: generates the HTML for the upload report
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

    def _generate_batch_report(self, structs_ref, workspace_name, structs_name,
                               successful_paths, failed_paths):
        """
        _generate_batch_report: generate summary report for upload
        """
        output_html_files = self._generate_report_html(structs_name, ','.join(successful_paths))

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

        file_path, workspace_name, pdb_name = self._validate_import_pdb_file_params(params)

        if not isinstance(workspace_name, int):
            workspace_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            workspace_id = workspace_name

        # TODO: Figure out how to parse the experimental pdb file for an experimental data structure
        data = dict()
        data, n_polypeptides = self._exp_file_to_data(file_path)
        # data['pdb_handle'] = self._upload_to_shock(file_path)
        # data['user_data'] = params.get('description', '')
        logging.info(data)

        info = self.dfu.save_objects({
            'id': workspace_id,
            'objects': [
                {'type': 'KBaseStructure.ExperimentalProteinStructure',
                 'name': pdb_name,
                 'data': data}]
        })[0]
        obj_ref = f"{info[6]}/{info[0]}/{info[4]}"

        returnVal = {'structure_obj_ref': obj_ref}

        if create_report:
            report_output = self._generate_report('import_experiment_pdb_file', obj_ref,
                                                  workspace_name, n_polypeptides,
                                                  pdb_name, file_path)
            returnVal.update(report_output)

        return returnVal

    def batch_import_pdb_files(self, params):
        """
        batch_import_pdb_files: upload two sets of pdb files and convert into a 
                               KBaseStructure.ProteinStructures object
        """

        workspace_name = params.get('workspace_name', None)
        if not isinstance(workspace_name, int):
            workspace_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            workspace_id = workspace_name

        exp_pdb_files = params.get('exp_pdb_file_paths', None)
        model_pdb_files = params.get('model_pdb_file_paths', None)
        structures_name = params.get('structures_name', 'protein_structures_name')

        model_pdb_objects = list()
        exp_pdb_objects = list()
        successful_files = list()
        failed_files = list()
        protein_structures = dict()
        total_structures = 0

        # loop through the lists model_pdb_file_paths and exp_pdb_file_paths
        if model_pdb_files:
            for model_file in model_pdb_files():
                model_pdb_ref = self.import_model_pdb_file(params, False)

                if model_pdb_ref:
                    model_pdb_objects.append(model_pdb_ref)
                    successful_files.append(model_file)
                    total_structures += 1
                else:
                    failed_files.append(model_file)

        if exp_pdb_files:
            for exp_file in exp_pdb_files():
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

    def export_pdb(self, params):
        if "input_ref" not in params:
            raise ValueError("input_ref not in supplied params")

        return {'shock_id': self._get_pdb_shock_id(params['input_ref'])}

    def structure_to_pdb_file(self, params):
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
