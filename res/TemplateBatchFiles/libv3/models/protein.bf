LoadFunctionLibrary ("terms.bf");
LoadFunctionLibrary ("frequencies.bf");

/** @module models.protein */

models.protein.alphabet = {{"A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"}};

/**
 * @name models.protein.generic.DefineQMatrix
 * @param {Dictionary} modelSpec
 * @param {String} namespace
 */

function models.protein.generic.DefineQMatrix (modelSpec, namespace) {

	__alphabet = modelSpec ["alphabet"];
	assert (Type (__alphabet) == "Matrix" && Columns (__alphabet) == 20, "Unsupported or missing alphabet '" + __alphabet + "'");

	__modelType = modelSpec["type"];
	if (Type (__modelType) == "None" || Type (__modelType) == "Number") {
		__modelType = "global";
	}
	modelSpec["type"] = __modelType;
	assert (__modelType == terms.local || __modelType == terms.global, "Unsupported or missing model type '" + __modelType + "'");

	__rate_function = modelSpec ["q_ij"];
	assert (utility.IsFunction (__rate_function), "Missing q_ij callback in model specification");

	__time_function = modelSpec ["time"];
	assert (utility.IsFunction (__time_function), "Missing time callback in model specification");


	__rate_matrix = {20,20};
	__rate_matrix [0][0] = "";

	__global_cache = {};

	for (_rowChar = 0; _rowChar < 20; _rowChar +=1 ){
		for (_colChar = _rowChar + 1; _colChar < 20; _colChar += 1) {
			__rp = Call (__rate_function, __alphabet[_rowChar],
															__alphabet[_colChar],
															 namespace,
															__modelType);


            if (Abs (__rp[terms.rate_entry])) {
                parameters.DeclareGlobal (__rp[terms.global], __global_cache);
                parameters.helper.copy_definitions (modelSpec["parameters"], __rp);

                __rate_matrix [_rowChar][_colChar] = __rp[terms.rate_entry];
                __rate_matrix [_colChar][_rowChar] = __rp[terms.rate_entry];
                continue;
            }
			__rate_matrix [_rowChar][_colChar] = "";
		}
	}

	__rp = Call (__time_function, __modelType);

	if (Abs (__rp)) {
		((modelSpec["parameters"])[terms.local])[terms.timeParameter ()] = __rp;
	    modelSpec [terms.rate_matrix] = parameters.AddMultiplicativeTerm (__rate_matrix, __rp, 0);
	} else {
	    modelSpec [terms.rate_matrix] = __rate_matrix;
	}

}

/**
 * @name models.protein.generic.Time
 * @param option - does nothing
 * @returns default time
 */
function models.protein.generic.Time (option) {
	return terms.default_time;
}



