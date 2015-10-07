LoadFunctionLibrary ("chooseGeneticCode", {"0" : "Universal"});

function models.codon.map_code (genetic_code) {
	return {"sense" : Columns (ComputeCodonCodeToStringMap (genetic_code)), 
	        "stop"  : Columns (ComputeCodonCodeToStringMapStop (genetic_code)),
	        "translation-table" : defineCodonToAAGivenCode (genetic_code) };
}

function models.codon.generic.defineQMatrix (modelSpec, namespace) {
	
	__base_alphabet = modelSpec ["bases"];
	assert (Type (__base_alphabet) == "Matrix" && Columns (__base_alphabet) == 4, "Unsupported codon bases '" + __base_alphabet + "'");
	
	__alphabet      = modelSpec ["alphabet"];
	// need at least one codon per amino-acid
	
	__dimension     = model.dimension (modelSpec);
	
	__stops         = modelSpec ["stop"];
	__table         = modelSpec ["translation-table"];
		
	assert (Type (__alphabet) == "Matrix" && __dimension > 20 && __dimension + Columns (__stops) && Abs (__table) == 64, "Unsupported or missing alphabet '" + __alphabet + "' (stop codons :'" + __codons + "')");

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
		
	__rate_matrix = {__dimension,__dimension};
	__rate_matrix [0][0] = "";
	__global_cache = {};
	
	
	for (_rowChar = 0; _rowChar < __dimension; _rowChar +=1 ){
		for (_colChar = _rowChar + 1; _colChar < __dimension; _colChar += 1) {
			    
            if (None == models.codon.diff (__alphabet[_rowChar], __alphabet[_colChar])) {
                continue;
            }
        
            __rp = utility.callFunction (__rate_function, {"0": parameters.quote(__alphabet[_rowChar]), 
                                                           "1": parameters.quote(__alphabet[_colChar]),
                                                           "2": "namespace",
                                                           "3": "__modelType",
                                                           "4": "__table"});
            
            
            if (Abs (__rp[terms.rate_entry])) {
                parameters.declareGlobal (__rp[terms.global], __global_cache);
                parameters.helper.copy_definitions (modelSpec["parameters"], __rp);
                __rate_matrix [_rowChar][_colChar] = __rp[terms.rate_entry];
                __rate_matrix [_colChar][_rowChar] = __rp[terms.rate_entry];
            }			
		}	
	}
	
	if (__modelType == "global") {
	    __rp = utility.callFunction (__time_function, {"0" : parameters.quote(__modelType)});
	
        if (Abs (__rp)) {
            ((modelSpec["parameters"])[terms.local])[terms.synonymous_rate] = __rp; 
            modelSpec [terms.rate_matrix] = parameters.addMultiplicativeTerm (__rate_matrix, __rp, 0);
        } else {
            modelSpec [terms.rate_matrix] = __rate_matrix;
        }
    } else {
        modelSpec [terms.rate_matrix] = __rate_matrix;
    }
	
	return modelSpec;
	
}


lfunction models.codon.diff (a,b) {
    r = {"from" : None,
         "to" : None,
         "position" : None};
    
    for (i = 0; i < 3; i += 1) {
        if (a[i] != b[i]) {
            if (r["position"] != None) {
                return None;
            }
            r["from"] = a[i];
            r["to"] = b[i];
            r["position"] = i;
        }
    }
    
    if (r ["position"] == None) {
        return None;
    }
    return r;
}
