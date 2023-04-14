
/*
@author chenry jmc jjeffryes tgu2 qzhang
*/
module KBaseStructure {
  typedef int bool;

  /*
    Reference to KBase object
    @id ws KBaseGenomes.Genome KBaseGenomeAnnotations.GenomeAnnotation KBaseGenomes.Feature
  */
  typedef string object_ref;

  /*
    type of a KBase object (e.g., type of a KBaseGenomes.Feature object)
  */
  typedef string object_type;

  /*
    Reference to KBase genome
    @id ws KBaseGenomes.Genome KBaseGenomeAnnotations.GenomeAnnotation
  */
  typedef string genome_ref;

  /*
    Reference to KBaseStructure objects
    @id ws KBaseStructure.ExperimentalProteinStructure KBaseStructure.ModelProteinStructure
  */
  typedef string structure_ref;

  /*
    Reference to KBase metagenome
    @id ws KBaseMetagenomes.AnnotatedMetagenomeAssembly KBaseGenomeAnnotations.Assembly
  */
  typedef string metagenome_ref;

  /*
    Reference to KBaseCollections.FeatureSet
    @id ws KBaseCollections.FeatureSet
  */
  typedef string feature_set_ref;

  /*
    CDS ID
    @id kb
  */
  typedef string cds_id;

  /*
    Molecule ID
    @id external
  */
  typedef string mol_id;

  /*
    Model ID
    @id external
  */
  typedef int mod_id;

  /*
    Uniref ID
    @id external
  */
  typedef string uniref_id;

  /*
    Reference to a file handle in shock
    @id handle
  */
  typedef string handle_ref;

  
  /*
    ProteinData
    mol_id id: ID for the protein
    string sequence: amino acid sequence
    string md5: hash of the amino acid sequence
    uniref_id uniref_id: from uniprot
    genome_ref genome_ref: from a KBase genome
    object_ref feature_ref: from a KBase feature
    object_type feature_type: from a KBase feature
    cds_id cds_id: from a KBase genome
    mod_id model_id: from PDB file
    mol_id chain_id: from PDB file
    float seq_identity: computed by comparing with KBase feature sequence
    bool exact_match: computed according to seq_identity

    @optional id uniref_id cds_id model_id chain_id seq_identity exact_match
    @optional genome_ref metagenome_ref feature_ref feature_set_ref feature_type
  */
  typedef structure {
    mol_id id;
    string sequence;
    string md5;
    uniref_id uniref_id;
    genome_ref genome_ref;
    object_ref feature_ref;
    object_type feature_type;
    cds_id cds_id;
    metagenome_ref metagenome_ref;
    feature_set_ref feature_set_ref;

    /*Parsed from PDB data and computed by comparing to the KBase feature sequence*/
    mod_id model_id;
    mol_id chain_id;
    float seq_identity;
    bool exact_match;
  } ProteinData;

  /*
    PDBInfo
    file_path: PDB structure file path
    file_extension: PDB structure file extenstion
    structure_name: PDB structure name
    narrative_id: KBase narrative id
    genome_name: KBase genome object name
    feature_id: KBase feature id
    is_model: 1 if structure is a computational model 0 otherwise
    from_rcsb': 1 if structure is from RCSB database 0 otherwise
    sequence_identities: sequence identities matched
    chain_ids: protein chain ids
    model_ids: model ids

    rcsb_id: The structure id for RCSB database
    genome_ref: name of a KBase genome object reference
    feature_type: id of a KBase feature object's type, default to 'gene'
    stratch_path: path on the shared folder where the structure file resides, default to file_path
    exact_matches: a string comma seperated '0' and '1' indicating an exact match not
    found ('0') or found ('1') for the structure's proteins with a given KBase genome feature

    @optional rcsb_id genome_ref feature_type scratch_path exact_matches
  */
  typedef structure {
    string file_path;
    string file_extension;
    string structure_name;
    int narrative_id;
    string genome_name;
    string feature_id;
    string sequence_identities;
    string chain_ids;
    string model_ids;
    bool is_model;
    bool from_rcsb;

    string rcsb_id;
    string genome_ref;
    string feature_type;
    string scratch_path;
    string exact_matches;
  } PDBInfo;

  /*
    ProteinStructure - merged from previous ModelProteinStructure and ExperimentalProteinStructure
    compound: a compound dict with keys in ['molecule', 'chain', 'synonym', 'misc', ...]
    source: a source dict with keys in ['organism_scientific', 'organism_taxid', 'other_details', 'organ', 'misc',...]
    @optional compound source
    @optional user_data num_models num_het_atoms num_water_atoms num_disordered_atoms num_disordered_residues
    @optional rcsb_id deposition_date head release_date structure_method resolution author
    @optional mmcif_handle xml_handle
  */
  typedef structure {
    string name;
    string user_data;
    int num_chains;
    int num_residues;
    int num_atoms;

    /*Experimental header from .cif file*/
    string rcsb_id;
    string deposition_date;
    string head;
    string release_date;
    string structure_method;
    float resolution;
    string author;

    /*Structure metadata from .cif file*/
    int num_models;
    int num_het_atoms;
    int num_water_atoms;
    int num_disordered_atoms;
    int num_disordered_residues;

    mapping<string, string> compound;
    mapping<string, string> source;

    /*Protein links*/
    list<ProteinData> proteins;

    /*File links*/
    handle_ref pdb_handle;
    handle_ref mmcif_handle;
    handle_ref xml_handle;

    /*Label provided by the user*/
    bool is_model;
  } ProteinStructure;

      
  /*
    ProteinStructures - using the merged ProteinStructure
    protein_structures: a list of ProteinStructure objects
    pdb_infos: a list of PDBInfo objects
    total_structures: total count of protein structures
    description: description/remarks
    @optional description
  */
  typedef structure {
    list<ProteinStructure> protein_structures;
    list<PDBInfo> pdb_infos;
    int total_structures;
    string description;
  } ProteinStructures;

};
