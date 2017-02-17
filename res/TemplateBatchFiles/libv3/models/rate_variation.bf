LoadFunctionLibrary("parameters.bf");
Loa

/** @module rate_variation */

rate_variation.types = {
	"Gamma": "rate_variation.types.Gamma.factory",
    "Gamma+I": "rate_variation.types.GammaI.factory"
};


lfunction rate_variation.types.Gamma.factory (options) {
	return {"distribution" : "rate_variation_define_gamma",
	        "options" : options, 
			"rate_modifier": "rate_variation.modifier_everything"};
}

lfunction rate_variation.types.GammaInv.factory (options) {
	return {"distribution" : "rate_variation_define_gamma_inv",
	        "options" : options, 
			"rate_modifier": "rate_variation.modifier_everything"};
}

lfunction rate_variation.types.GDD.factory (options) {
	return {"distribution" : "rate_variation_define_gdd",
	        "options" : options, 
			"rate_modifier": "rate_variation.modifier_everything"};
}

function rate_variation.modifier_everything (q_ij, from, to, namespace, cat_name) {
	if (Abs (q_ij[terms.rate_entry])) {
		q_ij[terms.rate_entry] = "(" + q_ij[terms.rate_entry] + ")*" + cat_name;
	}
	return q_ij;
}

lfunction rate_variation_define_gamma (options, namespace) {

	alpha = parameters.ApplyNameSpace("rv_gamma_alpha", namespace);
	parameters.DeclareGlobalWithRanges (alpha, 0.5, 0.01, 100);
								 
	if (utility.Has (options, ^"terms.rate_variation.bins", "Number")) {
		categories = options[^"terms.rate_variation.bins"] $ 1;
		assert (categories > 1 && categories < 64, 'Invalid number of `^"terms.rate_variation.bins"` in rate_variation_define_gamma');
	} else {
		categories = 4;
	}
	
	definition = {"id" : parameters.ApplyNameSpace("rv_gamma", namespace),
				  "description" : "Discretized unit mean Gamma distribution with " + categories + " categories", 
				  "category parameters" : {
				  		"bins" : ""+categories,
				  		"weights" : "EQUAL",
				  		"represent" : "MEAN",
				  		"PDF": "GammaDist(_x_,`alpha`,`alpha`)",
				  		"CDF": "CGammaDist(_x_,`alpha`,`alpha`)",
				  		^"terms.lower_bound" : 0,
				  		^"terms.upper_bound" : 1e25,
				  		"dCDF": "CGammaDist(_x_,`alpha`+1,`alpha`)"
				  },
				  "before" : None,
				  "after" : None,
				  };
				  
	return {^"terms.global" : {"Gamma distribution shape parameter" : alpha},
		    ^"terms.category" : definition};


}

lfunction rate_variation_define_gdd(options, namespace) {

	if (utility.Has (options, ^"terms.rate_variation.bins", "Number")) {
		categories = options[^"terms.rate_variation.bins"] $ 1;
		assert (categories > 1 && categories < 32, 'Invalid number of `^"terms.rate_variation.bins"` in rate_variation_define_gamma [2-32 expected]');
	} else {
		categories = 3;
	}
	
	
	globals = {};
		
	rates   = parameters.GenerateSequentialNames (parameters.ApplyNameSpace("rv_gdd_rates", namespace), categories, "_");
	
	utility.ForEachPair (rates, "_key_", "_value_", '(`&globals`)["GDD rate category " + (1+_key_)] = _value_');
	
	weights = parameters.GenerateSequentialNames (parameters.ApplyNameSpace("rv_gdd_weights", namespace), categories-1, "_");
	utility.ForEachPair (weights, "_key_", "_value_", '(`&globals`)[terms.mixture_aux_weight + " " + (1+_key_)] = _value_');

	parameters.DeclareGlobalWithRanges (weights, 1/categories, 0, 1);
	parameters.DeclareGlobal (rates, {});
	//console.log (rates);
	//utility.ForEach ({1,categories}["_MATRIX_ELEMENT_COLUMN_"], "_index_", '^(`&rates`[_index_])=_index_+1');
	//console.log (rates);
	if (utility.Has (options, ^"terms.initial_values", "String")) {
		if (options[^"terms.initial_values"] == "Randomize") {
			utility.ForEachPair (rates, "_key_", "_value_", '^_key_ = Random (0,1)');
		}
	} else {
	
	}
	
	weight_vector = parameters.helper.stick_breaking(weights, {categories,1}["1/categories"]);
	mean = Join ("+", utility.Map ({1,categories}["_MATRIX_ELEMENT_COLUMN_"], "_index_", "'('+`&rates`[_index_]+')*('+`&weight_vector`[_index_]+')'"));

	normalizer = parameters.ApplyNameSpace("rv_gdd_norm", namespace);
	parameters.SetConstraint (normalizer, mean, "global");
	
	definition = {"id" : parameters.ApplyNameSpace("rv_gdd", namespace),
				  "description" : "Discretized unit mean Gamma distribution with " + categories + " categories", 
				  "category parameters" : {
				  		"bins" : ""+categories,
				  		"weights" : "{{" + Join (",", weight_vector) + "}}",
				  		"represent" : "MEAN",
				  		"CDF": "{{" + Join (",", utility.Map (rates, "_value_", "'('+_value_+')/`normalizer`'")) + "}}",
				  		^"terms.lower_bound" : 0,
				  		^"terms.upper_bound" : 1e25,
				  },
				  "before" : None,
				  "after" : None,
				  };
				  
	return {^"terms.global" : globals,
		    ^"terms.category" : definition};


}

lfunction rate_variation_define_gamma_inv (options, namespace) {

	alpha = parameters.ApplyNameSpace("rv_gamma_inv_alpha", namespace);
	beta  = parameters.ApplyNameSpace("rv_gamma_inv_beta", namespace);
	p     = parameters.ApplyNameSpace("rv_gamma_inv_p", namespace);
	
	parameters.DeclareGlobalWithRanges (alpha, 0.5, 0.01, 100);
	parameters.DeclareGlobalWithRanges (beta, 1, 0.01, 200);
								 
	if (utility.Has (options, ^"terms.rate_variation.bins", "Number")) {
		categories = options[^"terms.rate_variation.bins"] $ 1;
		assert (categories > 1 && categories < 64, 'Invalid number of `^"terms.rate_variation.bins"` in rate_variation_define_gamma_inv');
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
	
	definition = {"id" : parameters.ApplyNameSpace("rv_gamma_inv", namespace),
				  "description" : "Discretized (unit mean) invariant class + Gamma distribution with " + categories + " categories", 
				  "category parameters" : {
				  		"bins" : ""+(categories+1),
				  		"weights" : weights,
				  		"represent" : "MEAN",
				  		"PDF": "(1-`p`)*GammaDist(_x_,`alpha`,`beta`)*(_x_>0)",
				  		"CDF": "(1-`p`)*CGammaDist(_x_,`alpha`,`beta`)*(_x_>0)+`p`",
				  		^"terms.lower_bound" : 0,
				  		^"terms.upper_bound" : 1e25,
				  		"dCDF": "(1-`p`)*CGammaDist(_x_,`alpha`+1,`alpha`)*(`alpha`/`beta`)*(_x_>0)"
				  },
				  "before" : None,
				  "after" : None,
				  };
				  
	return {^"terms.global" : {"Gamma distribution shape parameter" : alpha,
							   "Gamma distribution variance parameter" : beta,
							   "Proportion of invariant sites" : p},
		    ^"terms.category" : definition};


}


