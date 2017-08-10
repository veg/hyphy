LoadFunctionLibrary("parameters.bf");


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

lfunction rate_variation_define_gamma (options, namespace) {


	alpha = parameters.ApplyNameSpace("rv_gamma_alpha", namespace);
	parameters.DeclareGlobalWithRanges (alpha, 0.5, 0.01, 100);
								 
	if (utility.Has (options, utility.getGlobalValue("terms.rate_variation.bins"), "Number")) {
		categories = options[utility.getGlobalValue("terms.rate_variation.bins")] $ 1;
		assert (categories > 1 && categories < 64, 'Invalid number of `utility.getGlobalValue("terms.rate_variation.bins")` in rate_variation_define_gamma');
	} else {
		categories = 4;
	}
	
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
				  
	return {utility.getGlobalValue("terms.global") : {"Gamma distribution shape parameter" : alpha},
		    utility.getGlobalValue("terms.category") : definition};


}

lfunction rate_variation_define_gdd(options, namespace) {

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
	utility.ForEachPair (weights, "_key_", "_value_", '(`&globals`)[utility.getGlobalValue("terms.mixture.mixture_aux_weight") + " " + (1+_key_)] = _value_');

	parameters.DeclareGlobalWithRanges (weights, 1/categories, 0, 1);
	parameters.DeclareGlobal (rates, {});
	//console.log (rates);
	//utility.ForEach ({1,categories}["_MATRIX_ELEMENT_COLUMN_"], "_index_", '^(`&rates`[_index_])=_index_+1');
	//console.log (rates);
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
				  utility.getGlobalValue("terms.description") : "Discretized unit mean Gamma distribution with " + categories + " categories", 
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
				  
	return {utility.getGlobalValue("terms.global") : globals,
		    utility.getGlobalValue("terms.category") : definition};


}

lfunction rate_variation_define_gamma_inv (options, namespace) {

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
				  
				  
	return {utility.getGlobalValue("terms.global") : { utility.getGlobalValue("terms.rate_variation.gamma_alpha") : alpha,
							                           utility.getGlobalValue("terms.rate_variation.gamma_beta") : beta,
							                           utility.getGlobalValue("terms.rate_variation.gamma_p_inv") : p},
		    utility.getGlobalValue("terms.category") : definition};


}


