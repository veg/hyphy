/*terms.local = "local";
terms.global = "global";
terms.category = "category";

terms.rate_entry = "rate-entry";
terms.default_time = "t";

terms.lower_bound = "LB";
terms.upper_bound = "UB";
terms.MLE   = "MLE";

terms.range01 = {
    "LB": "0",
    "UB": "1"
};

terms.range_almost_01 = {
    "LB": "1e-8",
    "UB": "1"
};

terms.range_gte1 = {
    "LB": "1",
    "UB": "1e26"
};

terms.freqs.4x1 = "Nucleotide 4x1 estimator";
terms.freqs.equal = "Equal frequencies";
terms.freqs.CF3x4 = "Corrected 3x4 frequency estimator";
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
terms.transition = "transision rate";
terms.transversion = "transversion rate";

terms.mixture      = "mixture";
terms.mixture_components      = "mixture components";

terms.mixture_weight = "Mixture weight";
terms.mixture_aux_weight = "Mixture auxiliary weight";

terms.initial_values = "Initial values";



function terms.nucleotideRate(fromC, toC) {
    return "Substitution rate from nucleotide " + fromC + " to nucleotide " + toC;
}

function terms.aminoacidRate(fromA, toA) {
    return "Substitution rate from aminoacid " + fromA + " to aminoacid " + toA;
}

function terms.timeParameter() {
    return "Evolutionary time parameter";
}

function terms.AddCategory (term, categoryID) {
    return term + " for category " + categoryID;
}
*/