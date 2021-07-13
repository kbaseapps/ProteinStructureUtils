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

  typedef structure {
      obj_ref input_ref;
      string destination_dir;
   } StructureToPDBFileParams;

  typedef structure {
      string file_path;
  } StructureToPDBFileOutput;

  funcdef structure_to_pdb_file(StructureToPDBFileParams params)
      returns (StructureToPDBFileOutput result) authentication required;

  /* Input of the export_pdb function
    obj_ref: generics object reference
  */
  typedef structure {
      obj_ref obj_ref;
  } ExportParams;

  typedef structure {
      string shock_id;
  } ExportOutput;

  funcdef export_pdb (ExportParams params) returns (ExportOutput result) authentication required;

  /* Input of the import_model_pdb_file and import_model_pdb_file functions
    input_shock_id: file shock id
    input_file_path: absolute file path
    input_staging_file_path: staging area file path
    structure_name: structure object name
    workspace_name: workspace name for object to be saved to

  */
  typedef structure {
      string input_shock_id;
      string input_file_path;
      string input_staging_file_path;
      string structure_name;
      string description;
      workspace_name workspace_name;
  } ImportPDBParams;

  typedef structure {
      string report_name;
      string report_ref;
      obj_ref structure_obj_ref;
  } ImportPDBOutput;

  /* import_model_pdb_file: import a ProteinStructure from PDB*/
  funcdef import_model_pdb_file (ImportPDBParams params) returns (ImportPDBOutput result) authentication required;

  /* import_experiment_pdb_file: import a ProteinStructure from PDB*/
  funcdef import_experiment_pdb_file (ImportPDBParams params) returns (ImportPDBOutput result) authentication required;


  /* Input of the batch_import_pdb_files
    exp_pdb_file_paths: a list experiment pdb files (absolute path)
    model_pdb_file_paths: a list model pdb files (absolute path)
    structures_name: Proteinstructures object name
    workspace_name: workspace name for object to be saved to
  */
  typedef structure {
      list<string> exp_pdb_file_paths;
      list<string> model_pdb_file_paths;
      string structures_name;
      workspace_name workspace_name;
  } BatchPDBImportParams;

  typedef structure {
      string structures_ref;
      string report_name;
      string report_ref;
  } BatchPDBImportOutput;

  /* batch_import_pdb_files: import a batch of ProteinStructures from PDB files*/
  funcdef batch_import_pdb_files (BatchPDBImportParams params) returns (BatchPDBImportOutput result) authentication required;

};
