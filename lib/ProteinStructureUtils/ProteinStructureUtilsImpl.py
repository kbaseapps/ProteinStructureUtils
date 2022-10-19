# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os

from ProteinStructureUtils.Utils.PDBUtil import PDBUtil
from ProteinStructureUtils.Utils.RCSBUtil import RCSBUtil
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
    GIT_URL = ""
    GIT_COMMIT_HASH = "835a6e397b8078d13d1076b2665bb35b5531a21f"

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
        #END_CONSTRUCTOR
        pass


    def batch_import_pdbs_from_metafile(self, ctx, params):
        """
        batch_import_pdbs_from_metafile: import a batch of ProteinStructures from PDB files
        :param params: instance of type "BatchPDBImportParams" (Input/Output
           of the batch_import_pdbs_from_metafile structures_name:
           Proteinstructures object name workspace_name: workspace name for
           object to be saved to metadata_staging_file_path: path to a
           spreadsheet file that lists the metadata of PDB files and their
           KBase metadata) -> structure: parameter
           "metadata_staging_file_path" of String, parameter
           "structures_name" of String, parameter "workspace_name" of type
           "workspace_name" (workspace name of the object)
        :returns: instance of type "BatchPDBImportOutput" -> structure:
           parameter "structures_ref" of String, parameter "report_name" of
           String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN batch_import_pdbs_from_metafile
        logging.info(f'Starting batch_import_pdbs_from_metafile with params:\n{params}.')
        self.config['USER_ID'] = ctx['user_id']
        self.pdb_util = PDBUtil(self.config)
        result = self.pdb_util.batch_import_pdbs(params)
        #END batch_import_pdbs_from_metafile

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method batch_import_pdbs_from_metafile return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def import_rcsb_structures(self, ctx, params):
        """
        :param params: instance of type "ImportRCSBParams" (Input/output of
           the import_rcsb_structures function rcsb_ids: a list of rcsb
           structure id's exts: a list of rcsb structure file extensions
           ('pdb' or 'cif') narrative_ids: a list of KBase narrative ids
           genome_names: a list of KBase genome names in the respective
           narratives of narrative_ids feature_ids: a list of KBase feature
           ids in the respective narratives of narrative_ids is_models: a
           list of 0s and/or 1s to indicate the structure is exprimental or
           computational structures_name: Proteinstructures object name
           workspace_name: workspace name for object to be saved to) ->
           structure: parameter "rcsb_ids" of list of String, parameter
           "exts" of list of String, parameter "narrative_ids" of list of
           String, parameter "genome_names" of list of String, parameter
           "feature_ids" of list of String, parameter "is_models" of list of
           type "boolean" (A boolean - 0 for false, 1 for true. @range (0,
           1)), parameter "structures_name" of String, parameter
           "workspace_name" of type "workspace_name" (workspace name of the
           object)
        :returns: instance of type "ImportRCSBStructOutput" -> structure:
           parameter "structures_ref" of String, parameter "report_name" of
           String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN import_rcsb_structures
        logging.info(f'Starting import_rcsb_structures with params:\n{params}.')
        self.config['USER_ID'] = ctx['user_id']
        self.rcsb_util = RCSBUtil(self.config)
        result = self.rcsb_util.batch_import_rcsbs(params)
        #END import_rcsb_structures

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method import_rcsb_structures return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def export_pdb_structures(self, ctx, params):
        """
        :param params: instance of type "ExportParams" (Input/output of the
           export_pdb_structures function input_ref: generics object
           reference) -> structure: parameter "input_ref" of type "obj_ref"
           (An X/Y/Z style reference @id ws)
        :returns: instance of type "ExportStructOutput" -> structure:
           parameter "shock_ids" of list of String
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN export_pdb_structures
        logging.info(f'Starting export_pdb_structures with params:\n{params}.')
        self.config['USER_ID'] = ctx['user_id']
        self.pdb_util = PDBUtil(self.config)
        result = self.pdb_util.export_pdb_structures(params)
        #END export_pdb_structures

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method export_pdb_structures return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def query_rcsb_structures(self, ctx, params):
        """
        :param params: instance of type "QueryRCSBStructsParams"
           (Input/output of the query_rcsb_structures function
           sequence_strings: a list of protein sequences uniprot_ids: a list
           of uniprot ids ec_numbers: a list of ec numbers inchis: a list of
           InChI strings smiles: a list of SMILES strings evalue_cutoff:
           threshold of homology search identity_cutoff: threshold for
           sequence identity match workspace_name: workspace name for objects
           to be saved to @optional sequence_strings uniprot_ids ec_numbers
           inchis smiles evalue_cutoff identity_cutoff) -> structure:
           parameter "sequence_strings" of list of String, parameter
           "uniprot_ids" of list of String, parameter "ec_numbers" of list of
           String, parameter "inchis" of list of String, parameter "smiles"
           of list of String, parameter "evalue_cutoff" of Double, parameter
           "identity_cutoff" of Double, parameter "logical_and" of type
           "boolean" (A boolean - 0 for false, 1 for true. @range (0, 1)),
           parameter "workspace_name" of type "workspace_name" (workspace
           name of the object)
        :returns: instance of type "QueryRCSBStructsOutput" -> structure:
           parameter "rcsb_ids" of list of String, parameter "rcsb_scores" of
           unspecified object, parameter "report_name" of String, parameter
           "report_ref" of String
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN query_rcsb_structures
        logging.info(f'Starting query_rcsb_structures with params:\n{params}.')
        self.config['USER_ID'] = ctx['user_id']
        self.rcsb_util = RCSBUtil(self.config)
        result = self.rcsb_util.querey_structure_info(params)
        #END query_rcsb_structures

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method query_rcsb_structures return value ' +
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
