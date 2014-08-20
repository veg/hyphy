terms.local 			= "local";
terms.global 			= "global";
terms.default_time 		= "t";
terms.lower_bound 		= "LB";
terms.upper_bound 		= "UB";
terms.freqs.4x1			= "Nucleotide 4x1 estimator";
terms.freqs.equal		= "Equal frequencies";
terms.rate_matrix       = "Q";
terms.efv_matrix        = "pi";
terms.efv_estimate_name	= "Equilibrium frequency estimator";
terms.efv_estimate		= "EFV";

function  terms.nucleotideRate (fromC, toC) {
	return "Substitution rate from nucleotide " + fromC + " to nucleotide " + toC;
}

function  terms.timeParameter () {
	return "Evolutionary time parameter";
}
