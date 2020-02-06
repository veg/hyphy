LoadFunctionLibrary ("../all-terms.bf");
LoadFunctionLibrary ("frequencies.bf");

/** @module models.protein */



/* Alphabet of amino acids */
models.protein.alphabet = {{"A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"}};


/* Available empirical models */
models.protein.empirical_models = {"LG" : "Generalist empirical model of protein evolution from Le and Gascuel (2008). Ref: https://doi.org/10.1093/molbev/msn067",
                                   "WAG" : "Generalist empirical model of protein evolution from Whelan and Goldman (2001). Ref: https://doi.org/10.1093/oxfordjournals.molbev.a003851",
                                   "JTT" : "Generalist empirical model of protein evolution from Jones, Taylor, and Thornton (1996). Ref: https://doi.org/10.1093/bioinformatics/8.3.275",
                                   "JC69" : "Generalist empirical model of protein evolution with equal exchangeability rates among all amino acids, also known as JC69.",
                                   "mtMet" : "Specialist empirical model of protein evolution for metazoan mitochondrial genomes from Le, Dang, and Le. (2017). Ref: 10.1186/s12862-017-0987-y",
                                   "mtVer" : "Specialist empirical model of protein evolution for vertebrate mitochondrial genomes from Le, Dang, and Le. (2017). Ref: 10.1186/s12862-017-0987-y",
                                   "mtInv" : "Specialist empirical model of protein evolution for invertebrate mitochondrial genomes from Le, Dang, and Le. (2017). Ref: 10.1186/s12862-017-0987-y",
                                   "gcpREV" : "Specialist empirical model of protein evolution for green plant chloroplast genomes from  from Cox and Foster (2013). Ref: https://doi.org/10.1016/j.ympev.2013.03.030",
                                   "HIVBm" : "Specialist empirical model of protein evolution for between-host HIV sequences from Nickle et al. (2007). Ref: https://doi.org/10.1371/journal.pone.0000503",
                                   "HIVWm" : "Specialist empirical model of protein evolution for within-host HIV sequences from Nickle et al. (2007). Ref: https://doi.org/10.1371/journal.pone.0000503"
                                  };


models.protein.dimensions = 20;


/**
 * @name models.protein.generic.Time
 * @param option - does nothing
 * @returns default time
 */
function models.protein.generic.Time (option) {
	return terms.parameters.default_time;
}




/* Function below relocated to protein/empirical.bf and protein/REV.bf, each.
/**
 * @name models.protein.generic.DefineQMatrix
 * @param {Dictionary} modelSpec
 * @param {String} namespace
 */
 /*
function models.protein.generic.DefineQMatrix (modelSpec, namespace) {

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
*/
