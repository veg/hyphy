/* General terms without a specific scope necessarily */
terms.local = "local";
terms.global = "global";
terms.category = "category";
terms.rate_entry = "rate-entry";
terms.default_time = "t";


terms.MLE            = "MLE";
terms.log_likelihood = "LogL";
terms.empirical      = "empirical";
terms.parameters     = "parameters";

terms.sites          = "sites";
terms.sequences      = "sequences";
terms.file           = "file";
terms.genetic_code   = "code";
terms.sense_codons   = "sense";
terms.stop_codons    = "stop";

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

/* Terms involved in model definitions to frequencies */
terms.freqs.4x1 = "Nucleotide 4x1 estimator";
terms.freqs.equal = "Equal frequencies";
terms.freqs.CF3x4 = "Corrected 3x4 frequency estimator";




/* Miscellaneous terms */
terms.misc.timer = "timer";




terms.amino_acid = "amino-acid";
terms.codon      = "codon";
terms.synonymous_sub_count = "synonymous substitution count";
terms.nonsynonymous_sub_count = "nonsynonymous substitution count";





terms.rate_matrix = "Q";
terms.efv_matrix = "pi";
terms.efv_estimate_name = "Equilibrium frequency estimator";
terms.efv_estimate = "EFV";
terms.freqs.20x1 = "Protein 20x1 estimator";
terms.freqs.predefined = "Based on a training set";

terms.lf.local.constrained = "Local Constrained";
terms.lf.global.constrained = "Global Constrained";

terms.synonymous_rate = "synonymous rate";
terms.nonsynonymous_rate = "non-synonymous rate";
terms.omega_ratio = "non-synonymous/synonymous rate ratio";

terms.rate_variation = "Rate variation";
terms.rate_variation.bins = "Rate variation bins";
terms.rate_variation.gamma_alpha = "Shape parameter for the gamma distribution (alpha)";
terms.rate_variation.gamma_beta = "Shape parameter for the gamma distribution (beta)";
terms.rate_variation.gamma_p_inv = "Estimated proportion of invariant sites";

terms.branch_length = "branch length";
terms.branch_length_scaler = "branch length scaler";
terms.branch_length_constrain = "branch length constrain";

terms.transition_transversion_ratio = "Transitition/transversion ratio";
terms.transition = "transition rate";
terms.transversion = "transversion rate";

terms.mixture      = "mixture";
terms.mixture_components      = "mixture components";

terms.mixture_weight = "Mixture weight";
terms.mixture_aux_weight = "Mixture auxiliary weight";

terms.initial_values = "Initial values";



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



/* model fits */
terms.model.alternative = "alternative";
terms.model.null = "null";
terms.model.LRT = "LRT";
terms.model.p_value = "p-value";
terms.model.MLE = terms.MLE;




/* Terms specific to tree structures */
terms.trees = {};
terms.trees.newick = "string";
terms.trees.newick_with_lengths = "string_with_lengths";
terms.trees.newick_annotated = "annotated_string";
terms.trees.branch_length = "branch length"; // MUST BE SAME AS THE JSON
terms.trees.model_map = "model_map";
terms.trees.partitioned = "partitioned";
terms.trees.model_list = "model_list";

terms.trees.internal = "internal";
terms.trees.leaf = "leaf";
terms.trees.name = "name";
terms.trees.tree = "tree";

// TODO: Unclear if these are actually usable
terms.trees.node_name = "Name";
terms.trees.children = "Children";
terms.trees.parent = "Parent";
terms.trees.filter_string = "filter-string";






/* Terms specific to alignments, data info */
terms.alignments.code = "code";
terms.alignments.stops = "stops";
//terms.alignments.stops = "stop";

terms.alignments.file           = "file";
terms.alignments.sites          = "sites";
terms.alignments.sequences      = "sequences";
terms.alignments.name_mapping   = "name-mapping";
terms.alignments.partitions     = "partitions";
terms.alignments.sequence       = "sequence"; ///////////////////////////////////////////////////
terms.alignments.dataset        = "dataset";
terms.alignments.data_filter    = "datafilter";
terms.alignments.name           = "name";
terms.alignments.coverage       = "coverage";
terms.alignments.mapping        = "mapping";
terms.alignments.ordering       = "ordering";
terms.alignments.is_constant    = "is_constant";
terms.alignments.default        = "default";

/* Terms specific to I/O */
terms.io.info = "info";
terms.io.requirements = "requirements";
terms.io.authors = "authors";
terms.io.contact = "contact";
terms.io.version = "version";






/*---- SJS -------*/

terms.json_ = "json";
terms.MLE = "MLE";
terms.p_value = "p";
terms.LR = "LR";
terms.constraint = "constraint";
terms.fix = "fix-me";
terms.trees_ = "Trees";
terms.filters = "Filters";
terms.likelihood_function = "LF";
terms.model = "model";
terms.length = "length";
terms.length_parameter = "length parameter";
terms.bases = "bases";
terms.codons = "codons";
terms.id = "id";
terms.omegas = "omegas";
terms.weights = "weights";
terms.omega = "omega";
terms.weight = "weight";
terms.f = "f";
/*----------------*/






/* Terms specific to model construction and fitting */

/* Terms specific to output JSONs */
terms.json = {};
terms.json.attribute = {};

/******* TO DO: WHICH ONE IS CORRECT. Assuming the string????? *********/
//terms.json.attribute.meta = {};
terms.json.attribute.meta           = "attributes";


terms.json.info = terms.io.info;

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
terms.json.parameters               = "parameters";
terms.json.AICc                     = "AIC-c";
terms.json.attribute.branch_length  = "branch length";
terms.json.attribute.node_label     = "node label";
terms.json.attribute.branch_label   = "branch label";
terms.json.attribute.test           = "test";
terms.json.attribute.background     = "background";
terms.json.attribute.meta.type      = "attribute type";
terms.json.attribute.meta.order     = "display order";
terms.json.attribute.data           = "data";
terms.json.branch_attributes        = "branch attributes";
terms.json.display_order            = "display order";
terms.json.tree.newick              = "newick";
terms.json.partitions               = "partitions";
terms.json.name_mapping             = "name-mapping";

terms.json.test_results = "test results";
terms.json.p_value = "p";
terms.json.input = "input"; // All analyses
terms.json.input.filename = "filename";
terms.json.input.sequences = "sequences";
terms.json.input.sites = "sites";
terms.json.tree_string = "tree string"; // Possibly redundant w/ terms.trees.newick_with_lengths ?
terms.json.tree_length = "tree length";

terms.json.site_log_likelihood = "site log likelihood";
terms.json.evidence_ratios = "evidence ratios";


terms.json.runtime = "runtime";
terms.json.version = "version";
terms.json.convergence_failures = "convergence failures";
terms.json.branch_annotations = "branch annotations";
terms.json.branch_lengths = "branch lengths";
terms.json.annotation_tag = "annotation tag";




