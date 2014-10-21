LoadFunctionLibrary ("terms.bf");
LoadFunctionLibrary ("frequencies.bf");

models.DNA.alphabet = {{"A","C","G","T"}};

function models.DNA.generic.defineQMatrix (modelSpec, namespace) {
	
	__alphabet = modelSpec ["alphabet"];
	assert (Type (__alphabet) == "Matrix" && Columns (__alphabet) == 4, "Unsupported or missing alphabet '" + __alphabet + "'");
	
	__modelType = modelSpec["type"];
	if (Type (__modelType) == "None" || Type (__modelType) == "Number") {
		__modelType = "global";
	}
	modelSpec["type"] = __modelType;
	assert (__modelType == terms.local || __modelType == terms.global, "Unsupported or missing model type '" + __modelType + "'");
	
	__rate_function = modelSpec ["q_ij"];
	assert (utility.isFunction (__rate_function), "Missing q_ij callback in model specification");

	__time_function = modelSpec ["time"];
	assert (utility.isFunction (__time_function), "Missing time callback in model specification");
	
	
	__rate_matrix = {4,4};
	__rate_matrix [0][0] = "";
	__global_cache = {};
	
	
	for (_rowChar = 0; _rowChar < 4; _rowChar +=1 ){
		for (_colChar = _rowChar + 1; _colChar < 4; _colChar += 1) {
			__rp = utility.callFunction (__rate_function, {"0": parameters.quote(__alphabet[_rowChar]), 
															   "1": parameters.quote(__alphabet[_colChar]),
															   "2": "namespace",
															   "3": "__modelType"});
															   
            if (Abs (__rp[terms.rate_entry])) {			
                parameters.declareGlobal (__rp[terms.global], __global_cache);
                parameters.helper.copy_definitions (modelSpec["parameters"], __rp);
                            
                __rate_matrix [_rowChar][_colChar] = __rp[terms.rate_entry];
                __rate_matrix [_colChar][_rowChar] = __rp[terms.rate_entry];
                continue;
            } 
			__rate_matrix [_rowChar][_colChar] = "";
		}	
	}
	
	__rp = utility.callFunction (__time_function, {"0" : parameters.quote(__modelType)});
	
	if (Abs (__rp)) {
		((modelSpec["parameters"])[terms.local])[terms.timeParameter ()] = __rp; 
	    modelSpec [terms.rate_matrix] = parameters.addMultiplicativeTerm (__rate_matrix, __rp, 0);
	} else {
	    modelSpec [terms.rate_matrix] = __rate_matrix;
	}
	
}

function models.DNA.generic.time (option) {
	return terms.default_time;
}



