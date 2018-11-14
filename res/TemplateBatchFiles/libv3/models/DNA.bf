LoadFunctionLibrary ("frequencies.bf");
LoadFunctionLibrary ("libv3/all-terms.bf");


/** @module models.DNA */

models.DNA.models = {{"GTR", "General time reversible model"}, 
                     {"HKY85", "Hasegawa Kishino Yano 85 (HKY85) model"},
                     {"JC69", "Jukes-Cantor 69 (JC69) model"}};

models.DNA.generators ={"GTR": "models.DNA.GTR.ModelDescription",
                        "HKY85": "models.DNA.HKY85.ModelDescription",
                        "JC69": "models.DNA.JC69.ModelDescription"};
                                       

models.DNA.alphabet = {{"A","C","G","T"}};
models.DNA.dimensions = 4;

/**
 * @name models.DNA.generic.DefineQMatrix
 * @param {Dictionary} modelSpec
 * @param {String} namespace
 */
function models.DNA.generic.DefineQMatrix (modelSpec, namespace) {
	
	__alphabet = modelSpec [terms.alphabet];
	assert (Type (__alphabet) == "Matrix" && Columns (__alphabet) == 4, "Unsupported or missing alphabet '" + __alphabet + "'");
	
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
	
	__rate_variation = model.generic.get_rate_variation (modelSpec);

	__global_cache = {};

	if (None != __rate_variation) {

    //TODO
		__rp = Call (__rate_variation[terms.rate_variation.distribution], __rate_variation[terms.rate_variation.options], namespace);
		__rate_variation [terms.id] = (__rp[terms.category])[terms.id];
		parameters.DeclareCategory   (__rp[terms.category]);
        parameters.helper.copy_definitions (modelSpec[terms.parameters], __rp);
	} 
	
	__rate_matrix = {4,4};
	__rate_matrix [0][0] = "";
	
	
	for (_rowChar = 0; _rowChar < 4; _rowChar +=1 ){
		for (_colChar = _rowChar + 1; _colChar < 4; _colChar += 1) {
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
                __rate_matrix [_colChar][_rowChar] = __rp[terms.model.rate_entry];
                
            } else {
				__rate_matrix [_rowChar][_colChar] = "";
			}
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

/**
 * @name models.DNA.generic.Time
 * @param option - does nothing
 * @returns default time
 */
function models.DNA.generic.Time (option) {
	return terms.parameters.default_time;
}



