LoadFunctionLibrary ("../tasks/genetic_code.bf");
LoadFunctionLibrary ("../UtilityFunctions.bf");

/** @module models.codon */

/**
 * @name models.codon.MapCode
 * @param {String} genetic_code
 * @returns {Dictionary} the sense, stop, and translation-table for the genetic code
 */
function models.codon.MapCode (genetic_code) {

	return {"sense" : utility.Values (genetic_code.ComputeCodonCodeToStringMap (genetic_code)),
	        "stop"  : utility.Values (genetic_code.ComputeCodonCodeToStringMapStop (genetic_code)),
	        "translation-table" : genetic_code.DefineCodonToAAGivenCode (genetic_code) };
}

/**
 * @name models.codon.generic.DefineQMatrix
 * @param modelSpect
 * @param namespace
 */
function models.codon.generic.DefineQMatrix (modelSpec, namespace) {

	__base_alphabet = modelSpec ["bases"];
	assert (Type (__base_alphabet) == "Matrix" && Columns (__base_alphabet) == 4, "Unsupported codon bases '" + __base_alphabet + "'");

	__alphabet      = modelSpec ["alphabet"];
	// need at least one codon per amino-acid

	__dimension     = model.Dimension (modelSpec);

	__stops         = modelSpec ["stop"];
	__table         = modelSpec ["translation-table"];

	assert (Type (__alphabet) == "Matrix" && __dimension > 20 && __dimension + Columns (__stops) && Abs (__table) == 64, "Unsupported or missing alphabet '" + __alphabet + "' (stop codons :'" + __codons + "')");

	__modelType = modelSpec["type"];
	if (Type (__modelType) == "None" || Type (__modelType) == "Number") {
		__modelType = terms.global;
	}
	modelSpec["type"] = __modelType;
	assert (__modelType == terms.local || __modelType == terms.global, "Unsupported or missing model type '" + __modelType + "'");

	__rate_function = modelSpec ["q_ij"];
	assert (utility.IsFunction (__rate_function), "Missing q_ij callback in model specification");

	__time_function = modelSpec ["time"];
	assert (utility.IsFunction (__time_function), "Missing time callback in model specification");

	__rate_matrix = {__dimension,__dimension};
	__rate_matrix [0][0] = "";

	__rate_variation = model.generic.get_rate_variation (modelSpec);

	__global_cache = {};

	if (None != __rate_variation) {

		__rp = Call (__rate_variation["distribution"], __rate_variation["options"], namespace);
		__rate_variation ["id"] = (__rp[terms.category])["id"];
				
		parameters.DeclareCategory   (__rp[terms.category]);
        parameters.helper.copy_definitions (modelSpec["parameters"], __rp);
	} 


	for (_rowChar = 0; _rowChar < __dimension; _rowChar +=1 ){
		for (_colChar = _rowChar + 1; _colChar < __dimension; _colChar += 1) {

            if (None == models.codon.diff (__alphabet[_rowChar], __alphabet[_colChar])) {
                continue;
            }

            __rp = Call (__rate_function, __alphabet[_rowChar],
                                          __alphabet[_colChar],
                                          namespace,
                                          __modelType,
                                          __table);


		 	if (None != __rate_variation) {
				__rp = Call (__rate_variation["rate_modifier"], 
							 __rp,
							 __alphabet[_rowChar],
							 __alphabet[_colChar],
							 namespace,
							 __rate_variation ["id"]);
 			}

            if (Abs (__rp[terms.model_description.rate_entry])) {
                parameters.DeclareGlobal (__rp[terms.global], __global_cache);
                parameters.helper.copy_definitions (modelSpec["parameters"], __rp);
                __rate_matrix [_rowChar][_colChar] = __rp[terms.model_description.rate_entry];
                __rate_matrix [_colChar][_rowChar] = __rp[terms.model_description.rate_entry];
            }
		}
	}

	if (__modelType == terms.global) {
	    __rp = Call (__time_function, __modelType);

        if (Abs (__rp)) {
            ((modelSpec["parameters"])[terms.local])[terms.synonymous_rate] = __rp;
            modelSpec [terms.model_description.rate_matrix] = parameters.AddMultiplicativeTerm (__rate_matrix, __rp, FALSE);
        } else {
            modelSpec [terms.model_description.rate_matrix] = __rate_matrix;
        }
    } else {
        modelSpec [terms.model_description.rate_matrix] = __rate_matrix;
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
