terms.local 			= "local";
terms.global 			= "global";
terms.rate_entry        = "rate-entry";
terms.default_time 		= "t";

terms.lower_bound 		= "LB";
terms.upper_bound 		= "UB";
terms.range01    		= {"LB" : "0", 
                           "UB" : "1"};
                           
terms.freqs.4x1			= "Nucleotide 4x1 estimator";
terms.freqs.equal		= "Equal frequencies";
terms.freqs.CF3x4       = "Corrected 3x4 frequency estimator";
terms.rate_matrix       = "Q";
terms.efv_matrix        = "pi";
terms.efv_estimate_name	= "Equilibrium frequency estimator";
terms.efv_estimate		= "EFV";

terms.synonymous_rate    = "synonymous rate";
terms.nonsynonymous_rate = "non-synonymous rate";
terms.omega_ratio        = "non-synonymous/synonymous rate ratio";
 
function  terms.nucleotideRate (fromC, toC) {
	return "Substitution rate from nucleotide " + fromC + " to nucleotide " + toC;
}

function  terms.aminoacidRate (fromA, toA) {
	return "Substitution rate from amino-acid " + fromA + " to amino-acid " + toA;
}

function  terms.timeParameter () {
	return "Evolutionary time parameter";
}
