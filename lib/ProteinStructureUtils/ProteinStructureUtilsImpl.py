# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os

from ProteinStructureUtils.Utils.PDBUtils import PDBUtil
#END_HEADER


class ProteinStructureUtils:
    '''
    Module Name:
    ProteinStructureUtils

    Module Description:
    A KBase module: ProteinStructureUtils
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.2"
    GIT_URL = "https://github.com/qzzhang/ProteinStructureUtils.git"
    GIT_COMMIT_HASH = "48d0e649c9b29b6036e062b5ab85861f9cebc717"

    #BEGIN_CLASS_HEADER
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.config['SDK_CALLBACK_URL'] = os.environ['SDK_CALLBACK_URL']
        self.config['KB_AUTH_TOKEN'] = os.environ['KB_AUTH_TOKEN']
        self.scratch = config['scratch']
        self.pdb_util = PDBUtil(self.config)
        #END_CONSTRUCTOR
        pass

    def structure_to_pdb_file(self, ctx, params):
        """
        :param params: instance of type "StructureToPDBFileParams" ->
           structure: parameter "input_ref" of type "obj_ref" (An X/Y/Z style
           reference @id ws), parameter "destination_dir" of String
        :returns: instance of type "StructureToPDBFileOutput" -> structure:
           parameter "file_path" of String
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN structure_to_pdb_file
        logging.info('Starting structure_to_pdb_file with params:\n{}'.format(params))
        result = self.pdb_util.structure_to_pdb_file(params)
        #END structure_to_pdb_file

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method structure_to_pdb_file return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def export_pdb(self, ctx, params):
        """
        :param params: instance of type "ExportParams" (Input of the
           export_pdb function obj_ref: generics object reference) ->
           structure: parameter "input_ref" of type "obj_ref" (An X/Y/Z style
           reference @id ws)
        :returns: instance of type "ExportOutput" -> structure: parameter
           "shock_id" of String
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN export_pdb
        logging.info('Starting export_pdb with params:\n{}'.format(params))
        result = self.pdb_util.export_pdb(params)
        #END export_pdb

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method export_pdb return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def import_model_pdb_file(self, ctx, params):
        """
        import_model_pdb_file: import a ProteinStructure from PDB
        :param params: instance of type "ImportPDBParams" (Input of the
           import_model_pdb_file and import_experiment_pdb_file functions
           input_shock_id: file shock id input_file_path: absolute file path
           input_staging_file_path: staging area file path structure_name:
           structure object name workspace_name: workspace name for object to
           be saved to) -> structure: parameter "input_shock_id" of String,
           parameter "input_file_path" of String, parameter
           "input_staging_file_path" of String, parameter "structure_name" of
           String, parameter "description" of String, parameter
           "workspace_name" of type "workspace_name" (workspace name of the
           object)
        :returns: instance of type "ImportPDBOutput" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String,
           parameter "structure_obj_ref" of type "obj_ref" (An X/Y/Z style
           reference @id ws)
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN import_model_pdb_file
        logging.info('Starting import_model_pdb_file with params:\n{}'.format(params))
        result = self.pdb_util.import_model_pdb_file(params)
        #END import_model_pdb_file

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method import_model_pdb_file return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def import_experiment_pdb_file(self, ctx, params):
        """
        import_experiment_pdb_file: import a ProteinStructure from PDB
        :param params: instance of type "ImportPDBParams" (Input of the
           import_model_pdb_file and import_experiment_pdb_file functions
           input_shock_id: file shock id input_file_path: absolute file path
           input_staging_file_path: staging area file path structure_name:
           structure object name workspace_name: workspace name for object to
           be saved to) -> structure: parameter "input_shock_id" of String,
           parameter "input_file_path" of String, parameter
           "input_staging_file_path" of String, parameter "structure_name" of
           String, parameter "description" of String, parameter
           "workspace_name" of type "workspace_name" (workspace name of the
           object)
        :returns: instance of type "ImportPDBOutput" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String,
           parameter "structure_obj_ref" of type "obj_ref" (An X/Y/Z style
           reference @id ws)
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN import_experiment_pdb_file
        logging.info('Starting import_experiment_pdb_file with params:\n{}'.format(params))
        result = self.pdb_util.import_experiment_pdb_file(params)
        #END import_experiment_pdb_file

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method import_experiment_pdb_file return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
