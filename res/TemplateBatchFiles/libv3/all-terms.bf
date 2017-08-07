// MAKE SURE TO SWAP EVERYTHING
//terms.model.id = "id"; --> terms.id
//terms.model.model = "model";      ---> terms.model
//terms.model.parameters = "parameters"; ---> terms.parameters
//terms.model.parameters.global = "global"; ---> terms.global 
//terms.model.parameters.local = "local";   ---> terms.local 
//terms.model.parameters.empirical = "empirical"; ---> terms.model.empirical
//terms.model.translation_table = "translation-table";   ----> terms.translation_table
//terms.model.alphabet = "alphabet";       ----> terms.alphabet
//terms.model.bases = "bases";  --> terms.bases
//terms.model.stop =  "stop"; -->  terms.stop_codons
// terms.efv_estimate_name = "Equilibrium frequency estimator"; // ---> terms.model.efv_estimate_name
//terms.efv_estimate = "EFV";  //   --> terms.fit.efv_estimate


/* Models and fits */
//terms.MLE              = "MLE";   ----> terms.fit.MLE
//terms.log_likelihood   = "LogL";   -----> terms.fit.log_likelihood
//terms.empirical        = "empirical";    ----> terms.model.empirical
//terms.data.id     -------> terms.id
//terms.model.id     -------> terms.id
//terms.attribute.test           = "test";      ------> terms.tree_attributes.test
//terms.attribute.background     = "background";------> terms.tree_attributes.background
//terms.data.ordering  = "ordering";     ----> terms.code.ordering
//terms.data.mapping  = "mapping";     ----> terms.code.mapping
//terms.rate_variation.description = "description"; ------> terms.desription
//terms.model.description = "description"; ------> terms.desription



terms.mpi.job_id = "job_id";
terms.mpi.callback = "callback";
terms.mpi.arguments = "arguments";
terms.mpi.Models = "Models";
terms.mpi.LikelihoodFunctions = "LikelihoodFunctions";
terms.mpi.Headers = "Headers";
terms.mpi.Variables = "Variables";
terms.mpi.Functions = "Functions";












terms.math.count = "Count";
terms.math.mean = "Mean";
terms.math.median = "Median";
terms.math.stddev = "Std.Dev";
terms.math.minimum = "Min";
terms.math.maximum = "Max";
terms.math.2.5 = "2.5%";
terms.math.97.5 = "97.5%";
terms.math.sum = "Sum";
terms.math.variance = "Variance";
terms.math.cov = "COV";
terms.math.skewness = "Skewness";
terms.math.kurtosis = "Kurtosis";
terms.math.square_sum = "Sq. sum";
terms.math.non_negative = "Non-negative";


terms.header = "header";//1
terms.rows = "rows";//1

terms.table_options.header = terms.header;
terms.table_options.minumum_column_width = "min-column-width";
terms.table_options.align    = "align";
terms.table_options.column_widths = "column-widths";

terms.alphabet       = "alphabet";
terms.amino_acid     = "amino-acid";
terms.codon          = "codon"; 
terms.sense_codons   = "sense";
terms.stop_codons    = "stop"; //TODO. This is actually ok.
terms.bases          = "bases"; 
terms.codons         = "codons"; 
terms.translation_table = "translation-table";

terms.synonymous_sub_count          = "synonymous substitution count";
terms.nonsynonymous_sub_count       = "nonsynonymous substitution count";

terms.genetic_code.synonymous             = "synonymous";
terms.genetic_code.nonsynonymous          = "nonsynonymous";
terms.genetic_code.weighting_matrix       = "weighting-matrix";
terms.genetic_code.count_stop_codons      = "count-stop-codons";
terms.genetic_code.EPS                    = "EPS";
terms.genetic_code.EPN                    = "EPN";
terms.genetic_code.OPS                    = "OPS";
terms.genetic_code.OPN                    = "OPN";
terms.genetic_code.NTP                    = "NTP";
terms.genetic_code.SS                     = "SS";
terms.genetic_code.NS                     = "NS";

terms.code   = "code";
terms.code.stops = "stops";
terms.code.ordering = "ordering";
terms.code.mapping = "mapping";



terms.diff.to = "to";
terms.diff.from = "from";
terms.diff.position = "position";



terms.transition_transversion_ratio = "Transitition/transversion ratio";
terms.transition                    = "transition rate";
terms.transversion                  = "transversion rate";




terms.branch_length_constrain       = "branch length constrain"; // THIS EXISTS BUT UNCLEAR IF CORRECT NAMESPACE
terms.constraint          = "constraint";
terms.fix                 = "fix-me";
terms.likelihood_function = "LF";
terms.model               = "model";
terms.number_precision    = "number-precision";
terms.three_way = "three-way";
terms.reduced  = "reduced";
terms.default_time = "t"; 

//EVERYBODY
terms.branch_length = "branch length";
terms.parameters    = "parameters";
terms.local         = "local";
terms.global        = "global";
terms.default        = "default"; // WAS: terms.data.default
terms.id            = "id";
terms.description = "description";

terms.alternative      = "alternative";
terms.null             = "null";
terms.LRT              = "LRT";
//terms.LR               = "LR";  // change all to terms.LRT
terms.p_value          = "p-value"; //
/* Codon model terms */
terms.omegas  = "omegas";
terms.weights = "weights";
terms.freqs   = "f";
terms.omega   = "omega";
terms.weight  = "weight";
terms.category = "category";


/* Frequencies */
terms.freqs.4x1 = "Nucleotide 4x1 estimator";
terms.freqs.equal = "Equal frequencies";
terms.freqs.CF3x4 = "Corrected 3x4 frequency estimator";    
terms.freqs.20x1 = "Protein 20x1 estimator";     
terms.freqs.predefined = "Based on a training set";    

/* Model description terms */

terms.model.branch_length_constrain       = "branch length constrain"; //?

terms.model.length           = "length";
terms.model.length_parameter = "length parameter";
terms.model.length_expression = "length-expression";

terms.model.type = "type";
terms.model.components = "components";
terms.model.description = "description";
terms.model.canonical   = "canonical";
terms.model.reversible  =  "reversible";
terms.model.efv_estimate_name =  "Equilibrium frequency estimator";
terms.model.efv_estimate = "EFV"; 
terms.model.empirical = "empirical";
terms.model.empirical_rates = "empirical-rates";
terms.model.get_branch_length       = "get-branch-length";
terms.model.set_branch_length       = "set-branch-length";
terms.model.constrain_branch_length = "constrain-branch-length";
terms.model.frequency_estimator = "frequency-estimator";
terms.model.q_ij = "q_ij";
terms.model.time = "time";
terms.model.defineQ = "defineQ";
terms.model.post_definition = "post-definition";
terms.model.rate_entry = "rate-entry";
terms.model.rate_matrix = "Q";
terms.model.efv_matrix = "pi";
terms.model.matrix_id = "matrix-id";
terms.model.efv_id = "efv-id";
terms.model.data = "data";


terms.model.branch_length_string    = "branch-length-string";
terms.model.branch_length_scaler    = "branch length scaler";





terms.omega_ratio                   = "non-synonymous/synonymous rate ratio"; // in estimators.bf
terms.synonymous_rate               = "synonymous rate";
terms.nonsynonymous_rate            = "non-synonymous rate";
terms.branch_selection_attributes   = "Branch selection attributes";  

terms.globals_only = "globals only";
terms.empirical_bayes_factor  = "BF";
terms.posterior     = "posterior";

terms.initial_values          = "Initial values";


/* Mixture */
terms.mixture                 = "mixture";
terms.mixture_components      = "mixture components";
terms.mixture_weight          = "Mixture weight";
terms.mixture_aux_weight      = "Mixture auxiliary weight";

/* Rate variation */
terms.rate_variation = "Rate variation";
terms.rate_variation.distribution = "distribution";
terms.rate_variation.options = "options";
terms.rate_variation.rate_modifier = "rate_modifier";

terms.category_parameters = "category parameters";
terms.category.bins = "bins";
terms.category.weights = "weights";
terms.category.represent = "represent";
terms.category.PDF = "PDF";
terms.category.CDF = "CDF";
terms.category.dCDF = "dCDF";


terms.before = "before";
terms.after  = "after";

terms.rate_variation.bins = "Rate variation bins";

terms.rate_variation.Gamma = "Gamma";
terms.rate_variation.GammaI = "Gamma+I";
terms.rate_variation.gamma_alpha = "Shape parameter for the gamma distribution (alpha)";
terms.rate_variation.gamma_beta = "Variance parameter for the gamma distribution (beta)";
terms.rate_variation.gamma_p_inv = "Estimated proportion of invariant sites";




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






terms.run_options.retain_lf_object   = "retain-lf-object";
terms.run_options.proportional_branch_length_scaler = "proportional-branch-length-scaler";
terms.run_options.retain_model_object = "retain-model-object";
terms.run_options.model_type          = "model-type";
terms.run_options.partitioned_omega = "partitioned-omega";


//terms.data.trees_         = "Trees";


/* Likelihood fit structures */
terms.fit.efv_estimate = "EFV";
terms.fit.evolutionary_time = "Evolutionary time parameter";
terms.fit.MLE     = "MLE";
terms.fit.ID      = "ID";
terms.fit.trees   = "Trees";
terms.fit.log_likelihood = "LogL";
terms.fit.filters    = "Filters";


/* Genetics */







/* Data */
terms.data.sites          = "sites";
terms.data.sequences      = "sequences";
terms.data.sequence       = "sequence";
terms.data.file           = "file";
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
//terms.branch_attributes = "branch_attributes";
//terms.trees.name = "name"; // todo: is this real?

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
terms.attribute.branch_length  = terms.branch_length;
terms.attribute.node_label     = "node label";
terms.attribute.branch_label   = "branch label";

terms.attribute.data           = "data";

/* JSON terms. Many of these are copied from above.*/
terms.json = {};

terms.json.json                     = "json"; 
terms.json.fits                     = "fits";
terms.json.timers                   = "timers";
terms.json.PMID                     = "PMID";
terms.json.PMCID                    = "PMCID"; // because we can't all have a pmid.
terms.json.trees                    = "trees"; 
terms.json.MLE                      = "MLE";
terms.json.headers                  = "headers";
terms.json.content                  = "content";
terms.json.rate_distributions       = "rate distributions";
terms.json.log_likelihood           = "log likelihood";
terms.json.estimated_parameters     = "estimated parameters";
terms.json.AICc                     = "AIC-c";
terms.json.model                    = "model";
terms.json.global                   = "Global model fit";
terms.json.branch_attributes        = "branch attributes";
terms.json.branch_annotations       = "branch annotations";
terms.json.branch_lengths           = "branch lengths";
terms.json.annotation_tag           = "annotation tag";
terms.json.display_order            = terms.attribute.meta.order;

terms.json.tree.newick              = "newick"; // NEEDED?
terms.json.partitions               = "data partitions";
terms.json.tested                   = "tested";    
terms.json.uncorrected_pvalue       = "uncorrected p-value";
terms.json.test_results = "test results";
terms.json.input = "input"; 
terms.json.file = "file name";
terms.json.sequences = "number of sequences";
terms.json.partition_count = "partition count";
terms.json.sites = "number of sites";
terms.json.tree_string = "tree"; 
terms.json.tree_length = "tree length";
terms.json.relative_site_rates = "Relative site rate estimates";
terms.json.analysis = "analysis";

terms.json.site_log_likelihood = "site log likelihood";
terms.json.evidence_ratios = "evidence ratios";
terms.json.options  = "options";
terms.json.runtime = "runtime";
terms.json.version = "version";
terms.json.convergence_failures = "convergence failures";



