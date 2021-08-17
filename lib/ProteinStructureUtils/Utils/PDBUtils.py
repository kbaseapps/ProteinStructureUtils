import hashlib
import logging
import os
import shutil
import uuid
import errno
import pandas as pd

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
                validates input params to import_model_pdb_file and import_experiment_pdb_file methods
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

        (compound, source) = self._get_compound_source(structure)

        data = {
            'name': structure.header.get('name', ''),
            'num_chains': chain_no,
            'num_residues': res_no,
            'num_atoms': atom_no,
            'compound': compound,
            'source': source,
            'proteins': [{
                'id': os.path.basename(file_path),
                'sequence': seq,
                'md5': hashlib.md5(seq.encode()).hexdigest()
            }],
        }

        return data, pp_no

    def _exp_file_to_data(self, file_path):
        """
            _exp_file_to_data:
                Do the PDB conversion--parse the experiment pdb file for creating a pdb data object
        """
        # TODO: Figure out how to parse the experimental pdb file for an experimental data structure
        parser = PDB.MMCIFParser()
        cif = file_path
        structure = parser.get_structure("PHA-L", cif)
        # structure.header = structure.header
        # print(f'The mmcif structure header for {cif}*************************:\n\n')
        # kwds = structure.header.get('keywords')
        # print(f'Keywords: {kwds}')
        # resl = structure.header.get('resolution')
        # print(f'Resolution: {resl}')

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

        hd = self._upload_to_shock(file_path)

        (cpd, src) = self._get_compound_source(structure)

        mmcif_data = {
            'name': structure.header.get('name', ''),
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
            'num_models': structure.header.get('num_models', 0),
            'num_chains': structure.header.get('num_chains', 0),
            'num_residues': structure.header.get('residues', 0),
            'num_atoms': structure.header.get('num_atoms', 0),
            'num_het_atoms': structure.header.get('num_het_atoms', 0),
            'num_water_atoms': structure.header.get('num_water_atoms', 0),
            'num_disordered_atoms': structure.header.get('num_disordered_atoms', 0),
            'num_disordered_residues': structure.header.get('num_disordered_residues', 0),
            'pdb_handle': hd,
            'mmcif_handle': hd,
            'xml_handle': hd,
            'proteins': [{
                'id': os.path.basename(file_path),
                'sequence': seq,
                'md5': hashlib.md5(seq.encode()).hexdigest()
            }]
        }

        return mmcif_data, pp_no

    def _get_compound_source(self, structure):
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
