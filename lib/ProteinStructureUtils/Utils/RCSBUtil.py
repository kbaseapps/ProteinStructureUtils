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


class RCSBUtil:

    # “Sample Value” threshold to restrict which alignments will be significant
    SAMPLE_THRESH = 1e-20

    # BLAST sequence identity threshold to determine which pdb structures will be
    # matched to a KBase genome/feature
    ANOTHER_THRESH = 0.6

    def _validate_rcsb_import_params(self, params):
        """
            _validate_rcsb_import_params:
                validates input params to import_rcsb_pdbs
        """

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.token = config['KB_AUTH_TOKEN']
        self.user_id = config['USER_ID']
        self.dfu = DataFileUtil(self.callback_url)
        self.hs = AbstractHandle(config['handle-service-url'])
        self.ws_client = Workspace(config['workspace-url'])
        self.shock_url = config['shock-url']

    def import_rcsb_pdbs(self, params, create_report=False):
        """
            import_rcsb_pdbs: upload pdb files from RCSB pdb database
        """
        logging.info(f'import_rcsb_pdbs pdb data structures with params: {params}')

        # file_path is the pdb file's working area path (after dfu.download_staging_file call)
        file_path, workspace_name, pdb_name = self._validate_rcsb_import_params(params)
