/************************************************************
This top section contains potential new terms which have NOT BEEN IMPLEMENTED. //\d shows counts of how many times observed to determine if should implement,=.
************************************************************/

terms.header = "header";//1
terms.rows = "rows";//1

terms.number_precision = "number-precision";
terms.table_options.header = terms.header; // 1
terms.table_options.minumum_column_width = "min-column-width";// 1
terms.table_options.align    = "align";// 1
terms.table_options.column_widths = "column-widths";// 1


terms.branch_length_constrain       = "branch length constrain"; // THIS EXISTS BUT UNCLEAR IF CORRECT NAMESPACE
terms.constraint          = "constraint";
terms.fix                 = "fix-me";
terms.likelihood_function = "LF";
terms.model               = "model";



//EVERYBODY
terms.branch_length = "branch length";
terms.parameters    = "parameters";
terms.local         = "local";
terms.global        = "global";
terms.default        = "default"; // WAS: terms.data.default


terms.category = "category";


/* Descriptions */


/* Frequencies */
//terms.model_description.model = "model"; // 
terms.freqs.4x1 = "Nucleotide 4x1 estimator";
terms.freqs.equal = "Equal frequencies";
terms.freqs.CF3x4 = "Corrected 3x4 frequency estimator";    
terms.freqs.20x1 = "Protein 20x1 estimator";     
terms.freqs.predefined = "Based on a training set";    

/* Model description terms */
terms.model_description.alphabet = "alphabet";
terms.model_description.description = "description";
terms.model_description.canonical   = "canonical";
terms.model_description.reversible  =  "reversible";
terms.model_description.efv_estimate_name =  "Equilibrium frequency estimator";
//terms.model_description.parameters = "parameters";
//terms.model_description.parameters.global = "global";
//terms.model_description.parameters.local = "local";
//terms.model_description.parameters.empirical = "empirical";
terms.model_description.empirical = "empirical";
terms.model_description.type = "type";
terms.model_description.get_branch_length = "get-branch-length";
terms.model_description.set_branch_length = "set-branch-length";
terms.model_description.constrain_branch_length = "constrain-branch-length";
terms.model_description.branch_length_string = "branch-length-string";
terms.model_description.branch_length_scaler  = "branch length scaler";

terms.model_description.efv_estimate = "EFV"; // woah buddy
terms.model_description.rate_entry = "rate-entry";
terms.model_description.rate_matrix = "Q";
terms.model_description.efv_matrix = "pi";
terms.model_description.matrix_id = "matrix-id";
terms.model_description.efv_id = "efv-id";
terms.model_description.id = "id";
terms.model_description.data = "data";

terms.model_description.frequency_estimator = "frequency-estimator";
terms.model_description.q_ij = "q_ij";
terms.model_description.time = "time";
terms.model_description.defineQ = "defineQ";
terms.model_description.post_definition = "post-definition";
terms.model_description.bases = "bases"; 
terms.model_description.stop =  "stop";
terms.model_description.translation_table = "translation-table";
terms.model_description.components = "components";


terms.amino_acid     = "amino-acid";
terms.codon          = "codon";
terms.sense_codons   = "sense";
terms.stop_codons    = "stop"; //TODO
//terms.bases          = "bases";
terms.codons         = "codons";


/* Mixture */
terms.mixture                 = "mixture";
terms.mixture_components      = "mixture components";
terms.mixture_weight          = "Mixture weight";
terms.mixture_aux_weight      = "Mixture auxiliary weight";
terms.initial_values          = "Initial values";

/* Rate variation */
terms.rate_variation = "Rate variation";
terms.rate_variation.description = "description";
terms.rate_variation.description = "distribution";
terms.rate_variation.rate_modifier = "rate_modifier";
terms.rate_variation.bins = "Rate variation bins";
terms.rate_variation.gamma_alpha = "Shape parameter for the gamma distribution (alpha)";
terms.rate_variation.gamma_beta = "Shape parameter for the gamma distribution (beta)";
terms.rate_variation.gamma_p_inv = "Estimated proportion of invariant sites";





/* I THHINK THESE ARE MODEL DESCRIPTIONS */

// terms.efv_estimate_name = "Equilibrium frequency estimator"; // ---> terms.model_description.efv_estimate_name
//terms.efv_estimate = "EFV";  //   --> terms.fit.efv_estimate
terms.default_time = "t";

terms.omega_ratio                   = "non-synonymous/synonymous rate ratio"; // in estimators.bf
terms.synonymous_rate               = "synonymous rate";
terms.nonsynonymous_rate            = "non-synonymous rate";



terms.branch_length_constrain       = "branch length constrain";
terms.transition_transversion_ratio = "Transitition/transversion ratio";
terms.transition                    = "transition rate";
terms.transversion                  = "transversion rate";
terms.synonymous_sub_count          = "synonymous substitution count"; //SLAC
terms.nonsynonymous_sub_count       = "nonsynonymous substitution count"; //SLAC


/* Rate variation */
//terms.rate_variation = "Rate variation";
//terms.rate_variation.bins = "Rate variation bins";
//terms.rate_variation.gamma_alpha = "Shape parameter for the gamma distribution (alpha)";
//terms.rate_variation.gamma_beta = "Shape parameter for the gamma distribution (beta)";
//terms.rate_variation.gamma_p_inv = "Estimated proportion of invariant sites";

/* Term functions */
function terms.nucleotideRate(fromC, toC) {
    return "Substitution rate from nucleotide " + fromC + " to nucleotide " + toC;
}

function terms.aminoacidRate(fromA, toA) {
    return "Substitution rate from aminoacid " + fromA + " to aminoacid " + toA;
}

// TODO: Why is this here?
function terms.timeParameter() {
    return "Evolutionary time parameter";
}

function terms.AddCategory (term, categoryID) {
    return term + " for category " + categoryID;
}





/* Generic */

terms.lf.local.constrained = "Local Constrained";
terms.lf.global.constrained = "Global Constrained";





terms.timers = "timers";
terms.timers.timer = "timer";
terms.timers.order = "order";


terms.lower_bound = "LB";
terms.upper_bound = "UB";
terms.range01 = {
    terms.lower_bound: "0",
    terms.upper_bound: "1"
};

terms.range_almost_01 = {
    terms.lower_bound: "1e-8",
    terms.upper_bound: "1"
};

terms.range_gte1 = {
    terms.lower_bound: "1",
    terms.upper_bound: "1e26"
};








/* Models and fits */
//terms.MLE              = "MLE";
//terms.log_likelihood   = "LogL";
//terms.empirical        = "empirical";
//

terms.alternative      = "alternative";
terms.null             = "null";
terms.LRT              = "LRT";
terms.LR               = "LR";
terms.p_value          = "p-value"; //   p_value = "p";


terms.length           = "length";
terms.length_parameter = "length parameter";



// Used in site-level selection methods MEME, FEL, estimators.bf
terms.run_options.retain_lf_object   = "retain-lf-object";
terms.run_options.proportional_branch_length_scaler = "proportional-branch-length-scaler";
terms.run_options.retain_model_object = "retain-model-object";
terms.run_options.model_type          = "model-type";
terms.run_options.partitioned_omega = "partitioned-omega";

/* Codon model terms */
//terms.id      = "id";
terms.omegas  = "omegas";
terms.weights = "weights";
terms.omega   = "omega";
terms.weight  = "weight";
terms.f       = "f";



//terms.data.trees_         = "Trees";


/* Likelihood fit structures */
//terms.fit.global = terms.global;
terms.fit.efv_estimate = "EFV";
terms.fit.evolutionary_time = "Evolutionary time parameter";
terms.fit.MLE     = "MLE";
terms.fit.ID      = "ID";
terms.fit.trees     = "Trees";
terms.fit.log_likelihood = "LogL";
//terms.fit.branch_length = "branch length"; 
//terms.fit.parameters = terms.parameters;
terms.fit.filters    = "Filters";



/* Genetics */
terms.code   = "code";
terms.code.stops = "stops";
terms.code.ordering = "ordering";  // GONE: terms.data.ordering  = "ordering";
terms.code.mapping = "mapping"; // GONE: terms.data.mapping  = "mapping";



terms.three_way = "three-way";
terms.reduced  = "reduced";



/* Data */
terms.data.sites          = "sites";
terms.data.sequences      = "sequences";
terms.data.sequence       = "sequence";
terms.data.file           = "file";
terms.data.id             = "id"; //NOTE: SAME VALUE AS terms.model_description.id
terms.data.cache          = "cache";
terms.data.name           = "name";
terms.data.name_mapping   = "name-mapping";
terms.data.partitions     = "partitions";
terms.data.tree           = "tree";
terms.data.dataset        = "dataset";
terms.data.data_filter    = "datafilter";
terms.data.coverage       = "coverage";
terms.data.filter_string  = "filter-string";
terms.data.is_constant    = "is_constant";
terms.data.pattern_id     = "pattern id";




/* Tree structures */
terms.trees = {};
//terms.trees.branch_length = "branch length";
terms.trees.newick = "string";
terms.trees.newick_with_lengths = "string_with_lengths";
terms.trees.newick_annotated = "annotated_string";
terms.trees.model_map = "model_map";
terms.trees.partitioned = "partitioned";
terms.trees.model_list = "model_list";

terms.tree_attributes.internal = "internal";
terms.tree_attributes.leaf = "leaf";
terms.tree_attributes.test   = "test";
terms.tree_attributes.background  = "background";

terms.trees.name = "name";

// TODO: Unclear if these are actually usable
//terms.trees.node_name = "Name";
//terms.trees.children = "Children";
//terms.trees.parent = "Parent";
//

/* Input/output banner */
terms.io.info         = "info";
terms.io.requirements = "requirements";
terms.io.reference    = "citation";
terms.io.authors      = "authors";
terms.io.contact      = "contact";
terms.io.version      = "version";


/* Attribute terms */
terms.attribute = {};
/******* TO DO: WHICH ONE IS CORRECT. Assuming the string????? *********/
terms.attribute.meta = "attributes"; //terms.attribute.meta = {};
terms.attribute.meta.type      = "attribute type";
terms.attribute.meta.order     = "display order";
terms.attribute.branch_length  = terms.trees.branch_length;
terms.attribute.node_label     = "node label";
terms.attribute.branch_label   = "branch label";
//terms.attribute.test           = "test";
//terms.attribute.background     = "background";

terms.attribute.data           = "data";

/* JSON terms. Many of these are copied from above.*/
terms.json = {};

terms.json.json                     = "json"; 
terms.json.fits                     = "fits";
terms.json.timers                   = "timers";
terms.json.pmid                     = "PMID";
terms.json.trees                    = "trees"; 
terms.json.MLE                      = "MLE";
terms.json.headers                  = "headers";
terms.json.content                  = "content";
terms.json.rate_distributions       = "rate distributions";
terms.json.log_likelihood           = "log likelihood";
terms.json.estimated_parameters     = "estimated parameters";
terms.json.AICc                     = "AIC-c";

terms.json.branch_attributes        = "branch attributes";
terms.json.branch_annotations       = "branch annotations";
terms.json.branch_lengths           = "branch lengths";
terms.json.annotation_tag           = "annotation tag";
terms.json.display_order            = terms.attribute.meta.order;

terms.json.tree.newick              = "newick"; // NEEDED?
terms.json.partitions               = "data partitions";
//terms.json.name_mapping             = terms.name_mapping;
terms.json.tested                   = "tested";    

terms.json.test_results = "test results";
terms.json.p_value = "P-value";
terms.json.input = "input"; 
terms.json.file = "file name";
terms.json.sequences = "number of sequences";
terms.json.partition_count = "partition count";
terms.json.sites = "number of sites";
terms.json.tree_string = "tree"; 
terms.json.tree_length = "tree length";

terms.json.site_log_likelihood = "site log likelihood";
terms.json.evidence_ratios = "evidence ratios";

terms.json.runtime = "runtime";
terms.json.version = "version";
terms.json.convergence_failures = "convergence failures";



