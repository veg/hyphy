LoadFunctionLibrary("parameters.bf");
LoadFunctionLibrary("libv3/UtilityFunctions.bf");

/** @module rate_variation */

rate_variation.types = {
	utility.getGlobalValue("terms.rate_variation.Gamma"): "rate_variation.types.Gamma.factory",
    utility.getGlobalValue("terms.rate_variation.GammaI"): "rate_variation.types.GammaI.factory"
};


lfunction rate_variation.types.Gamma.factory (options) {
	return {utility.getGlobalValue("terms.rate_variation.distribution") : "rate_variation_define_gamma",
	        utility.getGlobalValue("terms.rate_variation.options") : options, 
			utility.getGlobalValue("terms.rate_variation.rate_modifier"): "rate_variation.modifier_everything"};
}

lfunction rate_variation.types.GammaInv.factory (options) {
	return {utility.getGlobalValue("terms.rate_variation.distribution") : "rate_variation_define_gamma_inv",
	        utility.getGlobalValue("terms.rate_variation.options") : options, 
			utility.getGlobalValue("terms.rate_variation.rate_modifier"): "rate_variation.modifier_everything"};
}

lfunction rate_variation.types.GDD.factory (options) {
	return {utility.getGlobalValue("terms.rate_variation.distribution") : "rate_variation_define_gdd",
	        utility.getGlobalValue("terms.rate_variation.options") : options, 
			utility.getGlobalValue("terms.rate_variation.rate_modifier"): "rate_variation.modifier_everything"};
}

function rate_variation.modifier_everything (q_ij, from, to, namespace, cat_name) {
	if (Abs (q_ij[terms.model.rate_entry])) {
		q_ij[terms.model.rate_entry] = "(" + q_ij[terms.model.rate_entry] + ")*" + cat_name;
	}
	return q_ij;
}


lfunction rate_variation_define_HMM (categories, namespace, globals, definition) {
    
    lambda = parameters.ApplyNameSpace   ("hmm_lambda", namespace);
    parameters.DeclareGlobalWithRanges (lambda, .1/categories, 0, 1/(categories-1));
    globals[utility.getGlobalValue("terms.rate_variation.hmm_lambda")] = lambda;

    hmm_T = parameters.ApplyNameSpace   ("hmm_transition_matrix", namespace);
    hmm_F = parameters.ApplyNameSpace   ("hmm_starting_freqs", namespace);
    hmm_M = parameters.ApplyNameSpace   ("hmm_model", namespace);
    matrix_def = ""; matrix_def * 128;
    
    matrix_def * "`hmm_T` = {`categories`, `categories`}; \n`hmm_F` = {`categories`, 1};\n";
    
    f_def = {categories,1};
    
    cm1 = categories-1;
    
    for (k = 0; k < categories; k+=1) {
        matrix_def * ('`hmm_T`[`k`][`k`] := 1-`cm1`*`lambda`;\n');
        matrix_def * ('`hmm_F`[`k`] = 1/`categories`;\n');
        //f_def[k] = "" + 1/categories;
        for (k2 = k+1; k2 < categories; k2+=1) {
        matrix_def * ('`hmm_T`[`k`][`k2`] := `lambda`;\n');
        matrix_def * ('`hmm_T`[`k2`][`k`] := `lambda`;\n');
        }
    }
    
    matrix_def * "Model `hmm_M` =  (`hmm_T`, `hmm_F`, 0);\n";
    matrix_def * 0;
    utility.ExecuteInGlobalNamespace (matrix_def);
    (definition[utility.getGlobalValue("terms.category.category_parameters")])[utility.getGlobalValue("terms.category.HMM")] = hmm_M;
}

lfunction rate_variation_define_gamma (options, namespace) {

 	if (utility.Has (options, utility.getGlobalValue("terms._namespace"), "String")) {
 	    namespace = options [utility.getGlobalValue("terms._namespace")];
 	}

	alpha = parameters.ApplyNameSpace("rv_gamma_alpha", namespace);
	parameters.DeclareGlobalWithRanges (alpha, 0.5, 0.01, 100);
								 
	if (utility.Has (options, utility.getGlobalValue("terms.rate_variation.bins"), "Number")) {
		categories = options[utility.getGlobalValue("terms.rate_variation.bins")] $ 1;
		assert (categories > 1 && categories < 64, 'Invalid number of `utility.getGlobalValue("terms.rate_variation.bins")` in rate_variation_define_gamma');
	} else {
		categories = 4;
	}
	
	globals = {utility.getGlobalValue("terms.rate_variation.gamma_alpha") : alpha};
	
	definition = {utility.getGlobalValue("terms.id") : parameters.ApplyNameSpace("rv_gamma", namespace),
				  utility.getGlobalValue("terms.description") : "Discretized unit mean Gamma distribution with " + categories + " categories", 
				  utility.getGlobalValue("terms.category.category_parameters") : {
				  		utility.getGlobalValue("terms.category.bins") : ""+categories,
				  		utility.getGlobalValue("terms.category.weights") : "EQUAL",
				  		utility.getGlobalValue("terms.category.represent") : "MEAN",
				  		utility.getGlobalValue("terms.category.PDF"): "GammaDist(_x_,`alpha`,`alpha`)",
				  		utility.getGlobalValue("terms.category.CDF"): "CGammaDist(_x_,`alpha`,`alpha`)",
				  		utility.getGlobalValue("terms.lower_bound") : 0,
				  		utility.getGlobalValue("terms.upper_bound") : 1e25,
				  		utility.getGlobalValue("terms.category.dCDF"): "CGammaDist(_x_,`alpha`+1,`alpha`)"
				  },
				  utility.getGlobalValue("terms.rate_variation.before") : None,
				  utility.getGlobalValue("terms.rate_variation.after") : None,
				  };

	if (utility.Has (options, utility.getGlobalValue("terms.rate_variation.HMM"), "String")) {
	    rate_variation_define_HMM (categories, namespace, globals, definition);
	}
				  
	return {utility.getGlobalValue("terms.global") : globals,
		    utility.getGlobalValue("terms.category") : definition};


}

lfunction rate_variation_define_gdd(options, namespace) {

 	if (utility.Has (options, utility.getGlobalValue("terms._namespace"), "String")) {
 	    namespace = options [utility.getGlobalValue("terms._namespace")];
 	}

	if (utility.Has (options, utility.getGlobalValue("terms.rate_variation.bins"), "Number")) {
		categories = options[utility.getGlobalValue("terms.rate_variation.bins")] $ 1;
		assert (categories > 1 && categories < 32, 'Invalid number of `utility.getGlobalValue("terms.rate_variation.bins")` in rate_variation_define_gamma [2-32 expected]');
	} else {
		categories = 3;
	}
	
	
	globals = {};
		
	rates   = parameters.GenerateSequentialNames (parameters.ApplyNameSpace("rv_gdd_rates", namespace), categories, "_");
	utility.ForEachPair (rates, "_key_", "_value_", '(`&globals`)["GDD rate category " + (1+_key_)] = _value_');
	weights = parameters.GenerateSequentialNames (parameters.ApplyNameSpace("rv_gdd_weights", namespace), categories-1, "_");
	utility.ForEachPair (weights, "_key_", "_value_", '(`&globals`)[utility.getGlobalValue("terms.mixture.mixture_aux_weight") + " for GDD category " + (1+_key_)] = _value_');

	parameters.DeclareGlobalWithRanges (weights, 1/categories, 0, 1);
	parameters.DeclareGlobal (rates, {});
	if (utility.Has (options, utility.getGlobalValue("terms.initial_values"), "String")) {
		if (options[utility.getGlobalValue("terms.initial_values")] == "Randomize") {
			utility.ForEachPair (rates, "_key_", "_value_", '^_key_ = Random (0,1)');
		}
	} else {
	
	}
	
	weight_vector = parameters.helper.stick_breaking(weights, {categories,1}["1/categories"]);
	mean = Join ("+", utility.Map ({1,categories}["_MATRIX_ELEMENT_COLUMN_"], "_index_", "'('+`&rates`[_index_]+')*('+`&weight_vector`[_index_]+')'"));

	normalizer = parameters.ApplyNameSpace("rv_gdd_norm", namespace);
	parameters.SetConstraint (normalizer, mean, utility.getGlobalValue("terms.global"));
	
	definition = {utility.getGlobalValue("terms.id"): parameters.ApplyNameSpace("rv_gdd", namespace),
				  utility.getGlobalValue("terms.description") : "General discrete unit mean variable on " + categories + " categories", 
				  utility.getGlobalValue("terms.category.category_parameters") : {
				  		utility.getGlobalValue("terms.category.bins") : ""+categories,
				  		utility.getGlobalValue("terms.category.weights") : "{{" + Join (",", weight_vector) + "}}",
				  		utility.getGlobalValue("terms.category.represent") : "MEAN",
				  		utility.getGlobalValue("terms.category.CDF"): "{{" + Join (",", utility.Map (rates, "_value_", "'('+_value_+')/`normalizer`'")) + "}}",
				  		utility.getGlobalValue("terms.lower_bound") : 0,
				  		utility.getGlobalValue("terms.upper_bound") : 1e25,
				  },
				  utility.getGlobalValue("terms.rate_variation.before") : None,
				  utility.getGlobalValue("terms.rate_variation.after") : None,
				  };
				
	if (utility.Has (options, utility.getGlobalValue("terms.rate_variation.HMM"), "String")) {
	    rate_variation_define_HMM (categories, namespace, globals, definition);
	}

	return {utility.getGlobalValue("terms.global") : globals,
		    utility.getGlobalValue("terms.category") : definition};


}

lfunction rate_variation_define_gamma_inv (options, namespace) {

 	if (utility.Has (options, utility.getGlobalValue("terms._namespace"), "String")) {
 	    namespace = options [utility.getGlobalValue("terms._namespace")];
 	}

	alpha = parameters.ApplyNameSpace("rv_gamma_inv_alpha", namespace);
	beta  = parameters.ApplyNameSpace("rv_gamma_inv_beta", namespace);
	p     = parameters.ApplyNameSpace("rv_gamma_inv_p", namespace);
	
	parameters.DeclareGlobalWithRanges (alpha, 0.5, 0.01, 100);
	parameters.DeclareGlobalWithRanges (beta, 1, 0.01, 200);
								 
	if (utility.Has (options, utility.getGlobalValue("terms.rate_variation.bins"), "Number")) {
		categories = options[utility.getGlobalValue("terms.rate_variation.bins")] $ 1;
		assert (categories > 1 && categories < 64, 'Invalid number of `utility.getGlobalValue("terms.rate_variation.bins")` in rate_variation_define_gamma_inv');
	} else {
		categories = 4;
	}

	globals = { utility.getGlobalValue("terms.rate_variation.gamma_alpha") : alpha,
							                           utility.getGlobalValue("terms.rate_variation.gamma_beta") : beta,
							                           utility.getGlobalValue("terms.rate_variation.gamma_p_inv") : p};

	parameters.DeclareGlobalWithRanges (p, 1/(categories+1), 0, 1);

	parameters.SetConstraint (beta, "(1-`p`)*`alpha`", "");
	weights = parameters.ApplyNameSpace("rv_gamma_inv_weights", namespace);

	weights = {1, 1+categories};
	weights[0] = p;
	

	for (k = 1; k <= categories; k+=1) {
		weights[k] = "(1-`p`)/" + categories;
	}

	weights = "{{" + Join (",", utility.Map (weights, "_value_", "_value_")) + "}}";
	
	definition = {utility.getGlobalValue("terms.id") : parameters.ApplyNameSpace("rv_gamma_inv", namespace),
				  utility.getGlobalValue("terms.description"): "Discretized (unit mean) invariant class + Gamma distribution with " + categories + " categories", 
				  utility.getGlobalValue("terms.category.category_parameters"): {
				  		utility.getGlobalValue("terms.category.bins"): ""+(categories+1),
				  		utility.getGlobalValue("terms.category.weights"): weights,
				  		utility.getGlobalValue("terms.category.represent"): "MEAN",
				  		utility.getGlobalValue("terms.category.PDF"): "(1-`p`)*GammaDist(_x_,`alpha`,`beta`)*(_x_>0)",
				  		utility.getGlobalValue("terms.category.CDF"): "(1-`p`)*CGammaDist(_x_,`alpha`,`beta`)*(_x_>0)+`p`",
				  		utility.getGlobalValue("terms.lower_bound"): 0,
				  		utility.getGlobalValue("terms.upper_bound"): 1e25,
				  		utility.getGlobalValue("terms.category.dCDF"): "(1-`p`)*CGammaDist(_x_,`alpha`+1,`alpha`)*(`alpha`/`beta`)*(_x_>0)"
                  },
				  utility.getGlobalValue("terms.rate_variation.before"): None,
				  utility.getGlobalValue("terms.rate_variation.after"): None,
				  };
				  

	if (utility.Has (options, utility.getGlobalValue("terms.rate_variation.HMM"), "String")) {
	    rate_variation_define_HMM (categories, namespace, globals, definition);
	}

	return {utility.getGlobalValue("terms.global") : globals,
		    utility.getGlobalValue("terms.category") : definition};


}

lfunction rate_variation.extract_category_information (model) {
    if (utility.Has (model[^"terms.parameters"], ^"terms.category", "AssociativeList")) {
        info = {};
        utility.ForEachPair ((model[^"terms.parameters"])[^"terms.category"], "_k_", "_v_", 
        "
            GetInformation (`&v_info`, ^_k_);
            (`&info`)[_k_] = `&v_info`;
            
        ");
        return info;
    }
    return None;
}

lfunction rate_variation.compute_mean (rv) {
    GetInformation (info, rv);
    return + (info[0][-1]$info[1][-1]);
    
}