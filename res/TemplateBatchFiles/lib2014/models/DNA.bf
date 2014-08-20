LoadFunctionLibrary ("terms.bf");
LoadFunctionLibrary ("frequencies.bf");

models.DNA.alphabet = "ACGT";


function models.DNA.generic.defineQMatrix (modelSpec, namespace) {
	
	__alphabet = modelSpec ["alphabet"];
	assert (Type (__alphabet) == "String" && Abs (__alphabet) == 4, "Unsupported or missing alphabet '" + __alphabet + "'");
	
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
		for (_colChar = 0; _colChar < 4; _colChar += 1) {
			if (_colChar != _rowChar) {
				__rp = utility.callFunction (__rate_function, {"0": parameters.quote(__alphabet[_rowChar]), 
															   "1": parameters.quote(__alphabet[_colChar])});
				if (__modelType == terms.global) {
					__rp = parameters.applyNameSpace (__rp, namespace);
					parameters.declareGlobal (__rp, __global_cache);
					((modelSpec["parameters"])[terms.global]) [terms.nucleotideRate (__alphabet[_rowChar], __alphabet[_colChar])] = __rp;
				} else {
					((modelSpec["parameters"])[terms.local]) [terms.nucleotideRate (__alphabet[_rowChar], __alphabet[_colChar])] = __rp;				
				}
				
				__rate_matrix [_rowChar][_colChar] = __rp;
			}
		}	
	}
	
	__rp = utility.callFunction (__time_function, {"0" : parameters.quote(__modelType)});
	
	if (Abs (__rp)) {
		((modelSpec["parameters"])[terms.local])[terms.timeParameter ()] = __rp; 
	}
	
	modelSpec [terms.rate_matrix] = parameters.addMultiplicativeTerm (__rate_matrix, __rp);
}

function models.DNA.generic.time (option) {
	return terms.default_time;
}

function models.matchAlphabets (a1, a2) {
	_validStates = {};
	for (_k = 0; _k < Abs (a1); _k += 1) {
		_validStates [a1[_k]] = 1;
	}
	for (_k = 0; _k < Abs (a2); _k += 1) {
		if (_validStates [a2[_k]] == 0) {
			return 0;
		}
	}
	return 1;

}

function models.DNA.generic.attachFilter (model, filter) {
	GetDataInfo (_givenAlphabet, *filter, "CHARACTERS");
	__alphabet = model ["alphabet"];
	assert (Abs (__alphabet) == Abs (_givenAlphabet) && matchAlphabets (_givenAlphabet, __alphabet), "The declared model alphabet '" + __alphabet + "' does not match the `filter` filter: '" + _givenAlphabet + "'");
	
	model ["alphabet"] = _givenAlphabet;
	model ["data"] = filter;
	return model;
}
