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

  /* Input/output of the export_pdb_structures function
    obj_ref: generics object reference
  */
  typedef structure {
      obj_ref input_ref;
  } ExportParams;

  typedef structure {
      list<string> shock_ids;
  } ExportStructOutput;

  funcdef export_pdb_structures (ExportParams params) returns (ExportStructOutput result) authentication required;
};
