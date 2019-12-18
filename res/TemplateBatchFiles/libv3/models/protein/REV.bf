LoadFunctionLibrary("../protein.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");

/** @module models.protein.REV */

/**
 * @name models.protein.REV.ModelDescription
 * @param {String} type
 */

lfunction models.protein.REV.ModelDescription(type) {
    return {
        utility.getGlobalValue("terms.alphabet"): utility.getGlobalValue("models.protein.alphabet"),
        utility.getGlobalValue("terms.description"): "General time reversible model for protein sequences",
        utility.getGlobalValue("terms.model.canonical"): 1, // is of the r_ij \times \pi_j form
        utility.getGlobalValue("terms.model.reversible"): 1,
        utility.getGlobalValue("terms.model.efv_estimate_name"): utility.getGlobalValue("terms.frequencies._20x1"),
        utility.getGlobalValue("terms.parameters"): {
            utility.getGlobalValue("terms.global"): {},
            utility.getGlobalValue("terms.local"): {},
            utility.getGlobalValue("terms.model.empirical"): 0
        },
        utility.getGlobalValue("terms.model.type"): type,
        utility.getGlobalValue("terms.model.get_branch_length"): "",
        utility.getGlobalValue("terms.model.set_branch_length"): "models.generic.SetBranchLength",
        utility.getGlobalValue("terms.model.constrain_branch_length"): "models.generic.ConstrainBranchLength",
        utility.getGlobalValue("terms.model.frequency_estimator"): "frequencies.empirical.protein",
        utility.getGlobalValue("terms.model.q_ij"): "models.protein.REV._GenerateRate",
        utility.getGlobalValue("terms.model.time"): "models.protein.generic.Time",
        utility.getGlobalValue("terms.model.defineQ"): "models.protein.REV._DefineQ",
        utility.getGlobalValue("terms.model.post_definition"): "models.generic.post.definition"
    };
}

lfunction models.protein.REV._GenerateRate(fromChar, toChar, namespace, model_type, model) {
    models.protein.REV._generateRate.p = {};
    models.protein.REV._generateRate.p[model_type] = {};

    if (fromChar > toChar) {
        models.protein.REV.parameter_name = "theta_" + toChar + fromChar;
    } else {
        models.protein.REV.parameter_name = "theta_" + fromChar + toChar;
    }

    if (model_type == utility.getGlobalValue("terms.global")) {
        models.protein.REV.parameter_name = parameters.ApplyNameSpace(models.protein.REV.parameter_name, namespace);
    }

    (models.protein.REV._generateRate.p[model_type])[terms.aminoacidRate(fromChar, toChar)] = models.protein.REV.parameter_name;
    models.protein.REV._generateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = models.protein.REV.parameter_name;

    return models.protein.REV._generateRate.p;

}

/**
 * @name models.protein.REV._DefineQ
 * @param {Dictionary} model definition
 * @param {String} namespace
 */
lfunction models.protein.REV._DefineQ(model_dict, namespace) {
    models.protein.REV.DefineQMatrix (model_dict, namespace);
    if (model_dict[utility.getGlobalValue("terms.model.type")] == utility.getGlobalValue("terms.global")) {
        parameters.SetConstraint(((model_dict[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.global")])[terms.aminoacidRate("I", "L")], "1", "");
    }
    if (model_dict[utility.getGlobalValue("terms.model.type")] == utility.getGlobalValue("terms.local")) {
        parameters.SetConstraint(((model_dict[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.local")])[terms.aminoacidRate("I", "L")], "1", "");
    }


    return model_dict;
}

/**
 * @name models.protein.REV.ModelDescription.withGamma
 * @description Define REV model with four-category gamma rate variation
 */
lfunction models.protein.REV.ModelDescription.withGamma (type) {
	def = models.protein.REV.ModelDescription (type);
	def [utility.getGlobalValue("terms.model.rate_variation")] = rate_variation.types.Gamma.factory ({utility.getGlobalValue("terms.rate_variation.bins") : 4});
	return def;
};


/**
 * @name models.protein.REV.ModelDescription.withGD4
 * @description Define REV model with 4bin General discrete rate variation
 */
lfunction models.protein.REV.ModelDescription.withGDD4 (type) {
	def = models.protein.REV.ModelDescription (type);
	def [utility.getGlobalValue("terms.model.rate_variation")] = rate_variation.types.GDD.factory ({utility.getGlobalValue("terms.rate_variation.bins") : 4});
	return def;
};

/**
 * @name models.protein.REV.DefineQMatrix
 * @param {Dictionary} modelSpec
 * @param {String} namespace
 */
function models.protein.REV.DefineQMatrix (modelSpec, namespace) {

	__alphabet = modelSpec [terms.alphabet];
	assert (Type (__alphabet) == "Matrix" && Columns (__alphabet) == models.protein.dimensions, "Unsupported or missing alphabet '" + __alphabet + "'");

	__modelType = modelSpec[terms.model.type];
	if (Type (__modelType) == "None" || Type (__modelType) == "Number") {
		__modelType = terms.global;
	}
	modelSpec[terms.model.type] = __modelType;
	assert (__modelType == terms.local || __modelType == terms.global, "Unsupported or missing model type '" + __modelType + "'");

	__rate_function = modelSpec [terms.model.q_ij];
	assert (utility.IsFunction (__rate_function), "Missing q_ij callback in model specification");

	__time_function = modelSpec [terms.model.time];
	assert (utility.IsFunction (__time_function), "Missing time callback in model specification");


	__rate_matrix = {models.protein.dimensions,models.protein.dimensions};
	__rate_matrix [0][0] = "";

	__rate_variation = model.generic.get_rate_variation (modelSpec);


	__global_cache = {};

	if (None != __rate_variation) {


		__rp = Call (__rate_variation[terms.rate_variation.distribution], __rate_variation[terms.rate_variation.options], namespace);
		__rate_variation [terms.id] = (__rp[terms.category])[terms.id];
				
		parameters.DeclareCategory   (__rp[terms.category]);
        parameters.helper.copy_definitions (modelSpec[terms.parameters], __rp);
	} 

	for (_rowChar = 0; _rowChar < models.protein.dimensions; _rowChar +=1 ){
		for (_colChar = 0; _colChar < models.protein.dimensions; _colChar += 1) {
		    
		    if (_rowChar == _colChar) {
		        continue;
		    }
			__rp = Call (__rate_function, __alphabet[_rowChar],
															  __alphabet[_colChar],
															   namespace,
															  __modelType,
															  modelSpec);
  

		 	if (None != __rate_variation) {
				__rp = Call (__rate_variation[terms.rate_variation.rate_modifier], 
							 __rp,
							 __alphabet[_rowChar],
							 __alphabet[_colChar],
							 namespace,
							 __rate_variation [terms.id]);
 			}

            if (Abs (__rp[terms.model.rate_entry])) {
                parameters.DeclareGlobal (__rp[terms.global], __global_cache);
                parameters.helper.copy_definitions (modelSpec[terms.parameters], __rp);

                __rate_matrix [_rowChar][_colChar] = __rp[terms.model.rate_entry];
                //__rate_matrix [_colChar][_rowChar] = __rp[terms.model.rate_entry];
                continue;
            }
			__rate_matrix [_rowChar][_colChar] = "";
		}
	}

	__rp = Call (__time_function, __modelType);

	if (Abs (__rp)) {
		((modelSpec[terms.parameters])[terms.local])[terms.timeParameter ()] = __rp;
	    modelSpec [terms.model.rate_matrix] = parameters.AddMultiplicativeTerm (__rate_matrix, __rp, 0);
	} else {
	    modelSpec [terms.model.rate_matrix] = __rate_matrix;
	}

}
