LoadFunctionLibrary ("../DNA.bf");
LoadFunctionLibrary ("../parameters.bf");
LoadFunctionLibrary ("../frequencies.bf");
LoadFunctionLibrary ("../../UtilityFunctions.bf");

function models.DNA.GTR.modelDescription () {

    return {"alphabet" : models.DNA.alphabet,
    		"description" : "The general time reversible (GTR) model of nucleotide substitution",
    		"canonical" : 1, // is of the r_ij \times \pi_j form
    		"reversible" : 1,
    		terms.efv_estimate_name: terms.freqs.4x1,
    		"parameters" : 	{
    				"global" : {}, 
    				"local" : {}
    			},
    		"get_branch_length" : "",
    		"set_branch_length" : "models.generic.set_branch_length",
    		"constrain_branch_length" : "models.generic.constrain_branch_length",
    		"frequency_estimator" : "frequencies.empirical.nucleotide",
    		"q_ij" : "models.DNA.GTR.generateRate",
    		"time" : "models.DNA.generic.time",
    		"defineQ" : "models.DNA.GTR.defineQ"
    		};
}


function models.DNA.GTR.generateRate (fromChar, toChar) {
	if (fromChar < toChar) {
		return "theta_" + fromChar + toChar;
	}
	return "theta_" + toChar + fromChar;
}

function models.DNA.GTR.defineQ (option, namespace) {
	gtr = models.DNA.GTR.modelDescription ();
	gtr ["type"] = option;
	models.DNA.generic.defineQMatrix (gtr, namespace);
	if (option == terms.global) {
		parameters.setConstraint (((gtr["parameters"])[terms.global])[terms.nucleotideRate ("A","G")], "1", "");
	}
	if (option == terms.local) {
		parameters.setConstraint (((gtr["parameters"])[terms.local])[terms.nucleotideRate ("A","G")], "1", "");
	}
	return gtr;
}
