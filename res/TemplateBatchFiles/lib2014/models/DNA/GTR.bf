LoadFunctionLibrary ("../DNA.bf");
LoadFunctionLibrary ("../parameters.bf");
LoadFunctionLibrary ("../frequencies.bf");
LoadFunctionLibrary ("../../UtilityFunctions.bf");

function models.DNA.GTR.modelDescription (type) {

    return {"alphabet" : models.DNA.alphabet,
    		"description" : "The general time reversible (GTR) model of nucleotide substitution",
    		"canonical" : 1, // is of the r_ij \times \pi_j form
    		"reversible" : 1,
    		terms.efv_estimate_name: terms.freqs.4x1,
    		"parameters" : 	{
    				"global" : {}, 
    				"local" : {}
    			},
    		"type" : type,
     		"get-branch-length" : "",
    		"set-branch-length" : "models.generic.set_branch_length",
    		"constrain-branch-length" : "models.generic.constrain_branch_length",
       		"frequency-estimator" : "frequencies.empirical.nucleotide",
    		"q_ij" : "models.DNA.GTR.generateRate",
    		"time" : "models.DNA.generic.time",
    		"defineQ" : "models.DNA.GTR.defineQ",
    		"post-definition" : "models.generic.post.definition"
    		};
}


function models.DNA.GTR.generateRate (fromChar, toChar, namespace, model_type) {
    models.DNA.GTR.generateRate.p = {};
    models.DNA.GTR.generateRate.p [model_type] = {};
    
	if (fromChar > toChar) {
		 models.DNA.GTR.parameter_name = "theta_" + toChar + fromChar;
	} else {
	    models.DNA.GTR.parameter_name = "theta_" + fromChar + toChar;
	}
	
	if (model_type == terms.global) {
	    models.DNA.GTR.parameter_name = parameters.applyNameSpace (models.DNA.GTR.parameter_name, namespace); 
	}
	
	(models.DNA.GTR.generateRate.p [model_type])[terms.nucleotideRate (fromChar, toChar)] = models.DNA.GTR.parameter_name;
	 models.DNA.GTR.generateRate.p [terms.rate_entry] = models.DNA.GTR.parameter_name;
	
	return models.DNA.GTR.generateRate.p;
}

function models.DNA.GTR.defineQ (gtr, namespace) {

	models.DNA.generic.defineQMatrix (gtr, namespace);
	if (gtr ["type"] == terms.global) {
		parameters.setConstraint (((gtr["parameters"])[terms.global])[terms.nucleotideRate ("A","G")], "1", "");
	} 
	if (gtr ["type"] == terms.local) {
		parameters.setConstraint (((gtr["parameters"])[terms.local])[terms.nucleotideRate ("A","G")], "1", "");
	}
	return gtr;
}
