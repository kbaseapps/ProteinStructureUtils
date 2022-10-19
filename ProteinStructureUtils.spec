/*
A KBase module: ProteinStructureUtils
*/

module ProteinStructureUtils {
  /* A boolean - 0 for false, 1 for true.
    @range (0, 1)
  */
  typedef int boolean;

  /* An X/Y/Z style reference
    @id ws
  */
  typedef string obj_ref;

  /* workspace name of the object */
  typedef string workspace_name;

  /* Input/Output of the batch_import_pdbs_from_metafile
    structures_name: Proteinstructures object name
    workspace_name: workspace name for object to be saved to
    metadata_staging_file_path: path to a spreadsheet file that lists the metadata of PDB files and their KBase metadata
  */
  typedef structure {
      string metadata_staging_file_path;
      string structures_name;
      workspace_name workspace_name;
  } BatchPDBImportParams;

  typedef structure {
      string structures_ref;
      string report_name;
      string report_ref;
  } BatchPDBImportOutput;

  /* batch_import_pdbs_from_metafile: import a batch of ProteinStructures from PDB files*/
  funcdef batch_import_pdbs_from_metafile (BatchPDBImportParams params) returns (BatchPDBImportOutput result) authentication required;

  /* Input/output of the import_rcsb_structures function
    rcsb_ids: a list of rcsb structure id's
    exts: a list of rcsb structure file extensions ('pdb' or 'cif')
    narrative_ids: a list of KBase narrative ids
    genome_names: a list of KBase genome names in the respective narratives of narrative_ids
    feature_ids: a list of KBase feature ids in the respective narratives of narrative_ids
    is_models: a list of 0s and/or 1s to indicate the structure is exprimental or computational
    structures_name: Proteinstructures object name
    workspace_name: workspace name for object to be saved to
  */
  typedef structure {
      list<string> rcsb_ids;
      list<string> exts;
      list<string> narrative_ids;
      list<string> genome_names;
      list<string> feature_ids;
      list<boolean> is_models;
      string structures_name;
      workspace_name workspace_name;
  } ImportRCSBParams;

  typedef structure {
      string structures_ref;
      string report_name;
      string report_ref;
  } ImportRCSBStructOutput;

  funcdef import_rcsb_structures (ImportRCSBParams params) returns (ImportRCSBStructOutput result) authentication required;

  /* Input/output of the export_pdb_structures function
    input_ref: generics object reference
  */
  typedef structure {
      obj_ref input_ref;
  } ExportParams;

  typedef structure {
      list<string> shock_ids;
  } ExportStructOutput;

  funcdef export_pdb_structures (ExportParams params) returns (ExportStructOutput result) authentication required;

  /* Input/output of the query_rcsb_structures function
    sequence_strings: a list of protein sequences
    uniprot_ids: a list of uniprot ids
    ec_numbers: a list of ec numbers
    inchis: a list of InChI strings
    smiles: a list of SMILES strings
    evalue_cutoff: threshold of homology search
    identity_cutoff: threshold for sequence identity match
    workspace_name: workspace name for objects to be saved to
    @optional sequence_strings uniprot_ids ec_numbers inchis smiles evalue_cutoff identity_cutoff
  */
  typedef structure {
      list<string> sequence_strings;
      list<string> uniprot_ids;
      list<string> ec_numbers;
      list<string> inchis;
      list<string> smiles;
      float evalue_cutoff;
      float identity_cutoff;
      boolean logical_and;
      workspace_name workspace_name;
  } QueryRCSBStructsParams;

  typedef structure {
      list<string> rcsb_ids;
      UnspecifiedObject rcsb_scores;
      string report_name;
      string report_ref;
  } QueryRCSBStructsOutput;

  funcdef query_rcsb_structures (QueryRCSBStructsParams params) returns (QueryRCSBStructsOutput result) authentication required;
};
