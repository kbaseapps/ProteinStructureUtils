
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
    @optional  genome_ref metagenome_ref feature_ref feature_set_ref feature_type
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
    ExperimentalProteinStructure
    compound: a compound dict with keys in ['molecule', 'chain', 'synonym', 'misc', ...]
    source: a source dict with keys in ['organism_scientific', 'organism_taxid', 'other_details', 'organ', 'misc',...]
    @optional mmcif_handle xml_handle
    @optional compound source
  */
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

    mapping<string, string> compound;
    mapping<string, string> source;

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

  /*
    ModelProteinStructure
    compound: a compound dict with keys in ['molecule', 'chain', 'synonym', 'misc', ...]
    source: a source dict with keys in ['organism_scientific', 'organism_taxid', 'other_details', 'organ', 'misc',...]
    @optional compound source
  */
  typedef structure {
    string user_data;
    string name;
    int num_chains;
    int num_residues;
    int num_atoms;

    mapping<string, string> compound;
    mapping<string, string> source;

    /*Protein links*/
    list<ProteinData> proteins;

    /*File links*/
    handle_ref pdb_handle;
    handle_ref mmcif_handle;
  } ModelProteinStructure;


  /*
    ProteinStructures
    model_structures: a list of references to ModelProteinStructures
    experimental_structures: a list of references to ExperimentalProteinStructures
    total_structures: total count of protein structures
    description: description/remarks
    @optional experimental_structures description
  */
  typedef structure {
    list<ModelProteinStructure> model_structures;
    list<ExperimentalProteinStructure> experimental_structures;
    int total_structures;
    string description;
  } ProteinStructures;

};
