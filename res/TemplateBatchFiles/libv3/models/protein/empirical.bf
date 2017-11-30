LoadFunctionLibrary("../protein.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");
LoadFunctionLibrary("../../all-terms.bf");



models.protein.empirical.default_generators = {"LG": "models.protein.LG.ModelDescription",
                                               "WAG": "models.protein.WAG.ModelDescription",
                                               "JTT": "models.protein.JTT.ModelDescription",
                                               "JC69": "models.protein.JC69.ModelDescription",
                                               "mtMAM": "models.protein.mtMAM.ModelDescription",
                                               "mtMet": "models.protein.mtMet.ModelDescription",
                                               "cpREV": "models.protein.cpREV.ModelDescription",
                                               "HIVBm": "models.protein.HIVBm.ModelDescription",
                                               "HIVWm": "models.protein.HIVWm.ModelDescription",                              
                                               "AB"   : "models.protein.AB.ModelDescription"};
                                           
models.protein.empirical.plusF_generators = {"LG": "models.protein.LGF.ModelDescription",
                                             "WAG": "models.protein.WAGF.ModelDescription",
                                             "JTT": "models.protein.JTTF.ModelDescription",
                                             "JC69": "models.protein.JC69F.ModelDescription",
                                             "mtMAM": "models.protein.mtMAMF.ModelDescription",                                           
                                             "mtMet": "models.protein.mtMetF.ModelDescription",
                                             "cpREV": "models.protein.cpREVF.ModelDescription",
                                             "HIVBm": "models.protein.HIVBmF.ModelDescription",
                                             "HIVWm": "models.protein.HIVWmF.ModelDescription",
                                             "AB"   : "models.protein.ABF.ModelDescription"};

models.protein.empirical.mleF_generators = {"LG": "models.protein.LGML.ModelDescription",
                                             "WAG": "models.protein.WAGML.ModelDescription",
                                             "JTT": "models.protein.JTTML.ModelDescription",
                                             "JC69": "models.protein.JC69ML.ModelDescription",
                                             "mtMAM": "models.protein.mtMAMML.ModelDescription",
                                              "mtMet": "models.protein.mtMetML.ModelDescription",
                                             "cpREV": "models.protein.cpREVML.ModelDescription",
                                             "HIVBm": "models.protein.HIVBmML.ModelDescription",
                                             "HIVWm": "models.protein.HIVWmML.ModelDescription",
                                             "AB"   : "models.protein.ABML.ModelDescription"};
/** @module models.protein.empirical */

/**
 * @name models.protein.empirical.ModelDescription
 * @param {String} type
 * @returns {Dictionary} model description
 * @description Create the baseline schema (dictionary) for empirical protein model definitions
 */
lfunction models.protein.empirical.ModelDescription(type) {
    return {
        utility.getGlobalValue("terms.alphabet"): utility.getGlobalValue("models.protein.alphabet"),
        utility.getGlobalValue("terms.description"): "General class of empirical substitution matrices for amino-acids",
        utility.getGlobalValue("terms.model.canonical"): 1, // is of the r_ij \times \pi_j form
        utility.getGlobalValue("terms.model.reversible"): 1,
        utility.getGlobalValue("terms.model.efv_estimate_name"): utility.getGlobalValue("terms.frequencies.predefined"),
        utility.getGlobalValue("terms.parameters"): {
            utility.getGlobalValue("terms.global"): {},
            utility.getGlobalValue("terms.local"): {},
            utility.getGlobalValue("terms.model.empirical"): 0
        },
        utility.getGlobalValue("terms.model.type"): type,
        utility.getGlobalValue("terms.model.get_branch_length"): "",
        utility.getGlobalValue("terms.model.set_branch_length"): "models.generic.SetBranchLength",
        utility.getGlobalValue("terms.model.constrain_branch_length"): "models.generic.constrain_branch_length",
        utility.getGlobalValue("terms.model.frequency_estimator"): "frequencies.empirical.protein",
        utility.getGlobalValue("terms.model.q_ij"): "models.protein.empirical._GenerateRate",
        utility.getGlobalValue("terms.model.time"): "models.protein.generic.Time",
        utility.getGlobalValue("terms.model.defineQ"): "models.protein.empirical._DefineQ",
        utility.getGlobalValue("terms.model.post_definition"): "models.generic.post.definition"
    };
}

/**
 * @name models.protein.Empirical._GenerateRate
 * @description Generates the r_ij component of q_ij := r_ij * time * freq_j
 * @param {Dict} rateDict
 * @param {Number} fromChar
 * @param {Number} toChar
 * @param {String} namespace
 * @param {String} model_type
 * @return list of parameters
 */
lfunction models.protein.empirical._GenerateRate(rateDict, fromChar, toChar, namespace, model_type) {


    models.protein.empirical._GenerateRate.p = {};
    models.protein.empirical._GenerateRate.p  [model_type]       = {};
    if (fromChar < toChar) {
        models.protein.empirical._GenerateRate.p  [utility.getGlobalValue("terms.model.rate_entry")] = "" + (rateDict[fromChar])[toChar];
    } else {
        models.protein.empirical._GenerateRate.p  [utility.getGlobalValue("terms.model.rate_entry")] = "" + (rateDict[toChar])[fromChar];
    }
    return models.protein.empirical._GenerateRate.p;
}

/**
 * @name models.protein.empirical._DefineQ
 * @param {Dictionary} model definition
 * @param {String} namespace
 */
lfunction models.protein.empirical._DefineQ(model_dict, namespace) {

    // Call frequencies here. Will be repeated in model.generic.DefineModel, but we are ok with that.
    frequencies._aux.empirical.singlechar(model_dict, namespace, model_dict[utility.getGlobalValue("terms.model.data")]);

    models.protein.empirical._NormalizeEmpiricalRates(model_dict, namespace); 
    models.protein.empirical.DefineQMatrix (model_dict, namespace);    
    return model_dict;
}




/**
 * @name models.protein.empirical._NormalizeEmpiricalRates
 * @param {Dictionary} model definition
 * @param {String} namespace
 */
lfunction models.protein.empirical._NormalizeEmpiricalRates(model_dict, namespace) {

    
    alphabet  = model_dict[utility.getGlobalValue("terms.alphabet")];
    dim       = utility.Array1D (alphabet);
    raw_rates = model_dict[utility.getGlobalValue("terms.model.empirical_rates")];
    EFV       = model_dict[utility.getGlobalValue("terms.efv_estimate")];


    // Create Q from R, EFV
    Q = {dim,dim};
    rowSums = {dim,1};

    // Fill in nondiagonal
    for (i = 0; i < dim; i +=1 ){
        for (j = i+1; j < dim; j += 1){
            if ( i!=j ){
                 
                
                rate = (raw_rates[alphabet[i]])[alphabet[j]];
                
                Q[i][j] = rate * EFV[j];
                Q[j][i] = rate * EFV[i];

                rowSums[i] += Q[i][j];
                rowSums[j] += Q[j][i];
            }
        }
    }

    // Fill in diagonal
    for (i = 0; i < dim; i += 1){
        Q[i][i] = -1 * rowSums[i];
    }


    // Get normalization factor
    norm = 0;
    for (i = 0; i < dim; i +=1 ){
        norm += Q[i][i] * EFV[i];
    }
    norm = -1*norm;
    
    // perform normalization
    for (i = 0; i < dim; i +=1 ){
        for ( j = 0; j < dim; j += 1){
            Q[i][j] = Q[i][j] / norm;
        }
    }   
    // Now convert it BACK TO hyphy dictionary with frequencies divided out. //
    // ************** This sets the new empirical rates. ************** //
    
    new_empirical_rates  = {};
    for (l1 = 0; l1 < dim - 1; l1 += 1) {
        new_empirical_rates[alphabet[l1]] = {};

        for (l2 = l1 + 1; l2 < dim; l2 += 1) {
            
            if (EFV[l2] == 0){
                nof_rate = 0.;
            }
            else {
                nof_rate =  Q [l1][l2] / EFV[l2];
            }
            (new_empirical_rates[alphabet[l1]])[alphabet[l2]] = nof_rate;
        }
    }
    model_dict[ utility.getGlobalValue("terms.model.empirical_rates")] = new_empirical_rates;

    return model_dict;
}





/**
 * @name models.protein.empirical.DefineQMatrix
 * @param {Dictionary} modelSpec
 * @param {String} namespace
 */
function models.protein.empirical.DefineQMatrix (modelSpec, namespace) {

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

    // ADDED FOR EMPIRICAL MODELS
    __empirical_rates = modelSpec[terms.model.empirical_rates];

    
    
    
	__global_cache = {};

	if (None != __rate_variation) {


		__rp = Call (__rate_variation[terms.rate_variation.distribution], __rate_variation[terms.rate_variation.options], namespace);
		__rate_variation [terms.id] = (__rp[terms.category])[terms.id];
				
		parameters.DeclareCategory   (__rp[terms.category]);
        parameters.helper.copy_definitions (modelSpec[terms.parameters], __rp);
	} 

	for (_rowChar = 0; _rowChar < models.protein.dimensions; _rowChar +=1 ){
		for (_colChar = _rowChar + 1; _colChar < models.protein.dimensions; _colChar += 1) {
		    
		    // NOTE the extra argument for protein models.
			__rp = Call (__rate_function, __empirical_rates, __alphabet[_rowChar],
															  __alphabet[_colChar],
															   namespace,
															  __modelType);
  

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





/**************************************** WAG functions *************************************/



/**
 * @name models.protein.WAG.ModelDescription
 * @description Create the baseline schema (dictionary) for the WAG model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.WAG.ModelDescription(type) {
    models.protein.WAG.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.WAG.ModelDescription.model_definition [terms.model.empirical_rates] = models.protein.WAG.Rij;
    models.protein.WAG.ModelDescription.model_definition [terms.model.frequency_estimator] = "models.protein.WAG.frequencies";
    return models.protein.WAG.ModelDescription.model_definition;
}

/**
 * @name models.protein.WAGF.ModelDescription
 * @description Create the baseline schema (dictionary) for the WAG+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.WAGF.ModelDescription(type) {
    models.protein.WAGF.ModelDescription.model_definition = models.protein.WAG.ModelDescription(type);
    models.protein.WAGF.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.empirical.protein";
    models.protein.WAGF.ModelDescription.model_definition [terms.model.efv_estimate_name] = utility.getGlobalValue("terms.frequencies._20x1");
    return models.protein.WAGF.ModelDescription.model_definition;
}

/**
 * @name models.protein.WAGML.ModelDescription
 * @description Create the baseline schema (dictionary) for the WAG+ML model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.WAGML.ModelDescription(type) {
    models.protein.WAGML.ModelDescription.model_definition = models.protein.WAG.ModelDescription(type);
    models.protein.WAGML.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.ML.protein";
    models.protein.WAGML.ModelDescription.model_definition [terms.model.efv_estimate_name]   =  utility.getGlobalValue("terms.frequencies.MLE");
    return models.protein.WAGML.ModelDescription.model_definition;
}




/**************************************** LG functions *************************************/


/**
 * @name models.protein.LG.ModelDescription
 * @description Create the baseline schema (dictionary) for the LG model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
 function models.protein.LG.ModelDescription(type) {
    models.protein.LG.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.LG.ModelDescription.model_definition [terms.model.empirical_rates] = models.protein.LG.Rij;
    models.protein.LG.ModelDescription.model_definition [terms.model.frequency_estimator] = "models.protein.LG.frequencies";
    return models.protein.LG.ModelDescription.model_definition;
}

/**
 * @name models.protein.LGF.ModelDescription
 * @description Create the baseline schema (dictionary) for the LG+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.LGF.ModelDescription(type) {
    models.protein.LGF.ModelDescription.model_definition = models.protein.LG.ModelDescription(type);
    models.protein.LGF.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.empirical.protein";
    models.protein.LGF.ModelDescription.model_definition [terms.model.efv_estimate_name] = utility.getGlobalValue("terms.frequencies._20x1");
    return models.protein.LGF.ModelDescription.model_definition;
}


/**
 * @name models.protein.LGML.ModelDescription
 * @description Create the baseline schema (dictionary) for the LG+ML model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.LGML.ModelDescription(type) {
    models.protein.LGML.ModelDescription.model_definition = models.protein.LG.ModelDescription(type);
    models.protein.LGML.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.ML.protein";
    models.protein.LGML.ModelDescription.model_definition [terms.model.efv_estimate_name]   =  utility.getGlobalValue("terms.frequencies.MLE");
    return models.protein.LGML.ModelDescription.model_definition;
}


/**************************************** JTT functions *************************************/


/**
 * @name models.protein.JTT.ModelDescription
 * @description Create the baseline schema (dictionary) for the JTT model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
 function models.protein.JTT.ModelDescription(type) {
    models.protein.JTT.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.JTT.ModelDescription.model_definition [terms.model.empirical_rates] = models.protein.JTT.Rij;
    models.protein.JTT.ModelDescription.model_definition [terms.model.frequency_estimator] = "models.protein.JTT.frequencies";
    return models.protein.JTT.ModelDescription.model_definition;
}


/**
 * @name models.protein.JTTF.ModelDescription
 * @description Create the baseline schema (dictionary) for the JTT+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.JTTF.ModelDescription(type) {
    models.protein.JTTF.ModelDescription.model_definition = models.protein.JTT.ModelDescription(type);
    models.protein.JTTF.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.empirical.protein";
    models.protein.JTTF.ModelDescription.model_definition [terms.model.efv_estimate_name] = utility.getGlobalValue("terms.frequencies._20x1");
    return models.protein.JTTF.ModelDescription.model_definition;
}


/**
 * @name models.protein.JTTML.ModelDescription
 * @description Create the baseline schema (dictionary) for the JTT+ML model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.JTTML.ModelDescription(type) {
    models.protein.JTTML.ModelDescription.model_definition = models.protein.JTT.ModelDescription(type);
    models.protein.JTTML.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.ML.protein";
    models.protein.JTTML.ModelDescription.model_definition [terms.model.efv_estimate_name]   =  utility.getGlobalValue("terms.frequencies.MLE");
    return models.protein.JTTML.ModelDescription.model_definition;
}

/**************************************** JC69 functions *************************************/


 /**
 * @name models.protein.JC69.ModelDescription
 * @description Create the baseline schema (dictionary) for the JC69 (equal rates) model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.JC69.ModelDescription(type) {
    models.protein.JC69.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.JC69.ModelDescription.model_definition [terms.model.empirical_rates] = models.protein.JC69.Rij;
    models.protein.JC69.ModelDescription.model_definition [terms.model.frequency_estimator] = "models.protein.JC69.frequencies";
    return models.protein.JC69.ModelDescription.model_definition;
}

/**
 * @name models.protein.JC69F.ModelDescription
 * @description Create the baseline schema (dictionary) for the JC69+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.JC69F.ModelDescription(type) {
    models.protein.JC69F.ModelDescription.model_definition = models.protein.JC69.ModelDescription(type);
    models.protein.JC69F.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.empirical.protein";
    models.protein.JC69F.ModelDescription.model_definition [terms.model.efv_estimate_name] = utility.getGlobalValue("terms.frequencies._20x1");
    return models.protein.JC69F.ModelDescription.model_definition;
}


/**
 * @name models.protein.JC69ML.ModelDescription
 * @description Create the baseline schema (dictionary) for the JC69+ML model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.JC69ML.ModelDescription(type) {
    models.protein.JC69ML.ModelDescription.model_definition = models.protein.JC69.ModelDescription(type);
    models.protein.JC69ML.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.ML.protein";
    models.protein.JC69ML.ModelDescription.model_definition [terms.model.efv_estimate_name]   =  utility.getGlobalValue("terms.frequencies.MLE");
    return models.protein.JC69ML.ModelDescription.model_definition;
}



/**************************************** mtMAM functions *************************************/


 /**
 * @name models.protein.mtMAM.ModelDescription
 * @description Create the baseline schema (dictionary) for the mtMAM model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.mtMAM.ModelDescription(type) {
    models.protein.mtMAM.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.mtMAM.ModelDescription.model_definition [terms.model.empirical_rates] = models.protein.mtMAM.Rij;
    models.protein.mtMAM.ModelDescription.model_definition [terms.model.frequency_estimator] = "models.protein.mtMAM.frequencies";
    return models.protein.mtMAM.ModelDescription.model_definition;
}

/**
 * @name models.protein.mtMAMF.ModelDescription
 * @description Create the baseline schema (dictionary) for the mtMAM+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.mtMAMF.ModelDescription(type) {
    models.protein.mtMAMF.ModelDescription.model_definition = models.protein.mtMAM.ModelDescription(type);
    models.protein.mtMAMF.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.empirical.protein";
    models.protein.mtMAMF.ModelDescription.model_definition [terms.model.efv_estimate_name] = utility.getGlobalValue("terms.frequencies._20x1");
    return models.protein.mtMAMF.ModelDescription.model_definition;
}


/**
 * @name models.protein.mtMAMML.ModelDescription
 * @description Create the baseline schema (dictionary) for the mtMAM+ML model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.mtMAMML.ModelDescription(type) {
    models.protein.mtMAMML.ModelDescription.model_definition = models.protein.mtMAM.ModelDescription(type);
    models.protein.mtMAMML.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.ML.protein";
    models.protein.mtMAMML.ModelDescription.model_definition [terms.model.efv_estimate_name]   =  utility.getGlobalValue("terms.frequencies.MLE");
    return models.protein.mtMAMML.ModelDescription.model_definition;
}



/**************************************** cpREV functions *************************************/

/**
 * @name models.protein.cpREV.ModelDescription
 * @description Create the baseline schema (dictionary) for the cpREV model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.cpREV.ModelDescription(type) {
    models.protein.cpREV.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.cpREV.ModelDescription.model_definition [terms.model.empirical_rates] = models.protein.cpREV.Rij;
    models.protein.cpREV.ModelDescription.model_definition [terms.model.frequency_estimator] = "models.protein.cpREV.frequencies";
    return models.protein.cpREV.ModelDescription.model_definition;
}

/**
 * @name models.protein.cpREVF.ModelDescription
 * @description Create the baseline schema (dictionary) for the cpREV+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.cpREVF.ModelDescription(type) {
    models.protein.cpREVF.ModelDescription.model_definition = models.protein.cpREV.ModelDescription(type);
    models.protein.cpREVF.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.empirical.protein";
    models.protein.cpREVF.ModelDescription.model_definition [terms.model.efv_estimate_name] = utility.getGlobalValue("terms.frequencies._20x1");
    return models.protein.cpREVF.ModelDescription.model_definition;
}

/**
 * @name models.protein.cpREVML.ModelDescription
 * @description Create the baseline schema (dictionary) for the cpREV+ML model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.cpREVML.ModelDescription(type) {
    models.protein.cpREVML.ModelDescription.model_definition = models.protein.cpREV.ModelDescription(type);
    models.protein.cpREVML.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.ML.protein";
    models.protein.cpREVML.ModelDescription.model_definition [terms.model.efv_estimate_name]   =  utility.getGlobalValue("terms.frequencies.MLE");
    return models.protein.cpREVML.ModelDescription.model_definition;
}


/**************************************** HIVBm functions *************************************/


/**
 * @name models.protein.HIVBm.ModelDescription
 * @description Create the baseline schema (dictionary) for the HIVBm model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.HIVBm.ModelDescription(type) {
    models.protein.HIVBm.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.HIVBm.ModelDescription.model_definition [terms.model.empirical_rates] = models.protein.HIVBm.Rij;
    models.protein.HIVBm.ModelDescription.model_definition [terms.model.frequency_estimator] = "models.protein.HIVBm.frequencies";
    return models.protein.HIVBm.ModelDescription.model_definition;
}

/**
 * @name models.protein.HIVBmF.ModelDescription
 * @description Create the baseline schema (dictionary) for the HIVBm+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.HIVBmF.ModelDescription(type) {
    models.protein.HIVBmF.ModelDescription.model_definition = models.protein.HIVBm.ModelDescription(type);
    models.protein.HIVBmF.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.empirical.protein";
    models.protein.HIVBmF.ModelDescription.model_definition [terms.model.efv_estimate_name] = utility.getGlobalValue("terms.frequencies._20x1");
    return models.protein.HIVBmF.ModelDescription.model_definition;
}

/**
 * @name models.protein.HIVBmML.ModelDescription
 * @description Create the baseline schema (dictionary) for the HIVBm+ML model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.HIVBmML.ModelDescription(type) {
    models.protein.HIVBmML.ModelDescription.model_definition = models.protein.HIVBm.ModelDescription(type);
    models.protein.HIVBmML.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.ML.protein";
    models.protein.HIVBmML.ModelDescription.model_definition [terms.model.efv_estimate_name]   =  utility.getGlobalValue("terms.frequencies.MLE");
    return models.protein.HIVBmML.ModelDescription.model_definition;
}



/**************************************** HIVWm functions *************************************/


/**
 * @name models.protein.HIVWm.ModelDescription
 * @description Create the baseline schema (dictionary) for the HIVWm model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.HIVWm.ModelDescription(type) {
    models.protein.HIVWm.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.HIVWm.ModelDescription.model_definition [terms.model.empirical_rates] = models.protein.HIVWm.Rij;
    models.protein.HIVWm.ModelDescription.model_definition [terms.model.frequency_estimator] = "models.protein.HIVWm.frequencies";
    return models.protein.HIVWm.ModelDescription.model_definition;
}

/**
 * @name models.protein.HIVWmF.ModelDescription
 * @description Create the baseline schema (dictionary) for the HIVWm+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.HIVWmF.ModelDescription(type) {
    models.protein.HIVWmF.ModelDescription.model_definition = models.protein.HIVWm.ModelDescription(type);
    models.protein.HIVWmF.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.empirical.protein";
    models.protein.HIVWmF.ModelDescription.model_definition [terms.model.efv_estimate_name] = utility.getGlobalValue("terms.frequencies._20x1");
    return models.protein.HIVWmF.ModelDescription.model_definition;
}

/**
 * @name models.protein.HIVWmML.ModelDescription
 * @description Create the baseline schema (dictionary) for the HIVWm+ML model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.HIVWmML.ModelDescription(type) {
    models.protein.HIVWmML.ModelDescription.model_definition = models.protein.HIVWm.ModelDescription(type);
    models.protein.HIVWmML.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.ML.protein";
    models.protein.HIVWmML.ModelDescription.model_definition [terms.model.efv_estimate_name]   =  utility.getGlobalValue("terms.frequencies.MLE");
    return models.protein.HIVWmML.ModelDescription.model_definition;
}


/**************************************** AB functions *************************************/



/**
 * @name models.protein.AB.ModelDescription
 * @description Create the baseline schema (dictionary) for the AB model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.AB.ModelDescription(type) {
    models.protein.AB.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.AB.ModelDescription.model_definition [terms.model.empirical_rates] = models.protein.AB.Rij;
    models.protein.AB.ModelDescription.model_definition [terms.model.frequency_estimator] = "models.protein.AB.frequencies";
    return models.protein.AB.ModelDescription.model_definition;
}

/**
 * @name models.protein.ABF.ModelDescription
 * @description Create the baseline schema (dictionary) for the AB+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.ABF.ModelDescription(type) {
    models.protein.ABF.ModelDescription.model_definition = models.protein.AB.ModelDescription(type);
    models.protein.ABF.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.empirical.protein";
    models.protein.ABF.ModelDescription.model_definition [terms.model.efv_estimate_name] = utility.getGlobalValue("terms.frequencies._20x1");
    return models.protein.ABF.ModelDescription.model_definition;
}

/**
 * @name models.protein.ABML.ModelDescription
 * @description Create the baseline schema (dictionary) for the AB+ML model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.ABML.ModelDescription(type) {
    models.protein.ABML.ModelDescription.model_definition = models.protein.AB.ModelDescription(type);
    models.protein.ABML.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.ML.protein";
    models.protein.ABML.ModelDescription.model_definition [terms.model.efv_estimate_name]   =  utility.getGlobalValue("terms.frequencies.MLE");
    return models.protein.ABML.ModelDescription.model_definition;
}




/**************************************** mtMet functions *************************************/


/**
 * @name models.protein.mtMet.ModelDescription
 * @description Create the baseline schema (dictionary) for the mtMet model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
 function models.protein.mtMet.ModelDescription(type) {
    models.protein.mtMet.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.mtMet.ModelDescription.model_definition [terms.model.empirical_rates] = models.protein.mtMet.Rij;
    models.protein.mtMet.ModelDescription.model_definition [terms.model.frequency_estimator] = "models.protein.mtMet.frequencies";
    return models.protein.mtMet.ModelDescription.model_definition;
}

/**
 * @name models.protein.mtMetF.ModelDescription
 * @description Create the baseline schema (dictionary) for the mtMet+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.mtMetF.ModelDescription(type) {
    models.protein.mtMetF.ModelDescription.model_definition = models.protein.mtMet.ModelDescription(type);
    models.protein.mtMetF.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.empirical.protein";
    models.protein.mtMetF.ModelDescription.model_definition [terms.model.efv_estimate_name] = utility.getGlobalValue("terms.frequencies._20x1");
    return models.protein.mtMetF.ModelDescription.model_definition;
}


/**
 * @name models.protein.mtMetML.ModelDescription
 * @description Create the baseline schema (dictionary) for the mtMet+ML model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.mtMetML.ModelDescription(type) {
    models.protein.mtMetML.ModelDescription.model_definition = models.protein.mtMet.ModelDescription(type);
    models.protein.mtMetML.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.ML.protein";
    models.protein.mtMetML.ModelDescription.model_definition [terms.model.efv_estimate_name]   =  utility.getGlobalValue("terms.frequencies.MLE");
    return models.protein.mtMetML.ModelDescription.model_definition;
}




/*=============================================================================================*/
/** Below this section are all of the empirical matrices and frequency vectors, including Rij **/



/**
 * @name models.protein.WAG.frequencies
 * @param {Dictionary} Baseline WAG model
 * @returns {Dictionary} Updated WAG model with empirical frequencies
 * @description Define the empirical amino acid frequencies associated with the WAG model of protein evolution
 */
lfunction models.protein.WAG.frequencies (model, namespace, datafilter) {
    model[utility.getGlobalValue("terms.efv_estimate")] =
        {{     0.0866279}
        {     0.0193078}
        {     0.0570451}
        {     0.0580589}
        {     0.0384319}
        {     0.0832518}
        {     0.0244313}
        {      0.048466}
        {     0.0620286}
        {      0.086209}
        {     0.0195027}
        {     0.0390894}
        {     0.0457631}
        {     0.0367281}
        {      0.043972}
        {     0.0695179}
        {     0.0610127}
        {     0.0708956}
        {     0.0143859}
        {     0.0352742}
        };

    model[utility.getGlobalValue("terms.model.efv_estimate_name")] = utility.getGlobalValue("terms.frequencies.predefined");
    (model[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.model.empirical")] = 0;
    return model;
}


/* Define a dictionary of amino-acid exchangeability rates for the WAG model of protein evolution. */ 
models.protein.WAG.Rij = {
 "A":{
   "C":0.02081175,
   "D":0.04424356,
   "E":0.09644885,
   "F":0.008490242,
   "G":0.1237844,
   "H":0.008127018999999999,
   "I":0.009834135000000001,
   "K":0.05899778,
   "L":0.03600239,
   "M":0.0182884,
   "N":0.02091646,
   "P":0.06909218,
   "Q":0.03502343,
   "R":0.02545459,
   "S":0.245933,
   "T":0.1358225,
   "V":0.1492591,
   "W":0.001708106,
   "Y":0.008912199000000001
  },
 "C":{
   "D":0.001813745,
   "E":0.001301055,
   "F":0.01605407,
   "G":0.02679532,
   "H":0.006383892,
   "I":0.008654049,
   "K":0.004819601,
   "L":0.03476936,
   "M":0.007992529,
   "N":0.0108821,
   "P":0.005254569,
   "Q":0.003809101,
   "R":0.02437562,
   "S":0.1027029,
   "T":0.03284827,
   "V":0.0745652,
   "W":0.01082647,
   "Y":0.02013312
  },
 "D":{
   "E":0.3762142,
   "F":0.001884863,
   "G":0.07562952000000001,
   "H":0.02386347,
   "I":0.002005993,
   "K":0.03123852,
   "L":0.007672926,
   "M":0.002123675,
   "N":0.2227414,
   "P":0.02036354,
   "Q":0.02377493,
   "R":0.00679797,
   "S":0.07819566,
   "T":0.02400406,
   "V":0.01133463,
   "W":0.001959249,
   "Y":0.01205807
  },
 "E":{
   "F":0.003272523,
   "G":0.04960369,
   "H":0.01461601,
   "I":0.006480045,
   "K":0.1682461,
   "L":0.01395734,
   "M":0.006450074,
   "N":0.0388587,
   "P":0.03277285,
   "Q":0.2108299,
   "R":0.02026676,
   "S":0.05143238,
   "T":0.0526847,
   "V":0.0438051,
   "W":0.00236373,
   "Y":0.007267292
  },
 "F":{
   "G":0.00436267,
   "H":0.01741975,
   "I":0.05389076,
   "K":0.005783216,
   "L":0.1913755,
   "M":0.02437025,
   "N":0.00394504,
   "P":0.007754001,
   "Q":0.003851615,
   "R":0.004740036,
   "S":0.03983115,
   "T":0.01100758,
   "V":0.04835584,
   "W":0.02309483,
   "Y":0.2389425
  },
 "G":{
   "H":0.006395123,
   "I":0.001548868,
   "K":0.02431859,
   "L":0.005546612,
   "M":0.003563543,
   "N":0.04617598,
   "P":0.01169843,
   "Q":0.0127224,
   "R":0.02698185,
   "S":0.09789926,
   "T":0.01446092,
   "V":0.01393229,
   "W":0.005087841,
   "Y":0.003835501
  },
 "H":{
   "I":0.007029141,
   "K":0.05796705,
   "L":0.04519012,
   "M":0.008272107000000001,
   "N":0.1623063,
   "P":0.03343772,
   "Q":0.1655236,
   "R":0.09862788,
   "S":0.05400277,
   "T":0.0303076,
   "V":0.008806540999999999,
   "W":0.003964322,
   "Y":0.1433978
  },
 "I":{
   "K":0.02108143,
   "L":0.2869017,
   "M":0.08714326,
   "N":0.02273747,
   "P":0.004799484,
   "Q":0.004391122,
   "R":0.008628942000000001,
   "S":0.02330635,
   "T":0.09337141,
   "V":0.5819514,
   "W":0.003208113,
   "Y":0.01555502
  },
 "K":{
   "L":0.02330296,
   "M":0.0191231,
   "N":0.1235674,
   "P":0.02674718,
   "Q":0.1501354,
   "R":0.246964,
   "S":0.07056185,
   "T":0.08881348999999999,
   "V":0.02272611,
   "W":0.002076079,
   "Y":0.004933538
  },
 "L":{
   "M":0.09935387,
   "N":0.005395922,
   "P":0.01997258,
   "Q":0.03351591,
   "R":0.02296714,
   "S":0.02515217,
   "T":0.02091482,
   "V":0.133956,
   "W":0.01004497,
   "Y":0.01475715
  },
 "M":{
   "N":0.008131996000000001,
   "P":0.008228767999999999,
   "Q":0.05956464,
   "R":0.03152741,
   "S":0.03603533,
   "T":0.0970828,
   "V":0.1531609,
   "W":0.007786238,
   "Y":0.01586107
  },
 "N":{
   "P":0.009369552999999999,
   "Q":0.05950219,
   "R":0.02932074,
   "S":0.28996,
   "T":0.1299922,
   "V":0.01460187,
   "W":0.001085813,
   "Y":0.04020457
  },
 "P":{
   "Q":0.03597839,
   "R":0.03135791,
   "S":0.1177049,
   "T":0.05093139,
   "V":0.02342947,
   "W":0.002104766,
   "Y":0.007998193000000001
  },
 "Q":{
   "R":0.140086,
   "S":0.07506641,
   "T":0.05493632,
   "V":0.0224171,
   "W":0.003257243,
   "Y":0.008430004
  },
 "R":{
   "S":0.08931696,
   "T":0.03550112,
   "V":0.01873906,
   "W":0.01757311,
   "Y":0.01412465
  },
 "S":{
   "T":0.2803409,
   "V":0.01731717,
   "W":0.007907568,
   "Y":0.0291351
  },
 "T":{
   "V":0.1032926,
   "W":0.001673848,
   "Y":0.01077852
  },
 "V":{
   "W":0.005516418,
   "Y":0.01165155
  },
 "W":{
   "Y":0.0920111
  }
};






/**
 * @name models.protein.LG.frequencies
 * @param {Dictionary} Baseline LG model
 * @returns {Dictionary} Updated LG model with empirical frequencies
 * @description Define the empirical amino acid frequencies associated with the LG model of protein evolution
 */
lfunction models.protein.LG.frequencies (model, namespace, datafilter) {
    model[utility.getGlobalValue("terms.efv_estimate")] =
 
          {{     0.07906500000000008}
           {     0.012937}
           {     0.053052}
           {     0.071586}
           {     0.042302}
           {     0.057337}
           {     0.022355}
           {     0.062157}
           {     0.0646}
           {     0.099081}
           {     0.022951}
           {     0.041977}
           {     0.04404}
           {     0.040767}
           {     0.055941}
           {     0.061197}
           {     0.053287}
           {     0.069147}
           {     0.012066}
           {     0.034155}
          };

    model[utility.getGlobalValue("terms.model.efv_estimate_name")] = utility.getGlobalValue("terms.frequencies.predefined");
    (model[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.model.empirical")] = 0;
    return model;
}

/* Define a dictionary of amino-acid exchangeability rates for the LG model of protein evolution. */ 
models.protein.LG.Rij = {
 "A":{
   "C":0.03522956,
   "D":0.0229346,
   "E":0.08133689,
   "F":0.01174132,
   "G":0.1296008,
   "H":0.008776704999999999,
   "I":0.01018879,
   "K":0.03791848,
   "L":0.04285406,
   "M":0.02822381,
   "N":0.01271276,
   "P":0.05674114,
   "Q":0.04325807,
   "R":0.02601647,
   "S":0.3164948,
   "T":0.1247291,
   "V":0.1927457,
   "W":0.002385593,
   "Y":0.008181845
  },
 "C":{
   "D":0.003630821,
   "E":0.0002740351,
   "F":0.05115122,
   "G":0.03570948,
   "H":0.01566596,
   "I":0.0218034,
   "K":0.0009375764,
   "L":0.06438966,
   "M":0.02243974,
   "N":0.02428347,
   "P":0.003632027,
   "Q":0.003782507,
   "R":0.0327155,
   "S":0.1864267,
   "T":0.06666287,
   "V":0.1482198,
   "W":0.00884617,
   "Y":0.04355245
  },
 "D":{
   "E":0.41069,
   "F":0.0008060157,
   "G":0.05300146,
   "H":0.02267472,
   "I":0.0007269456,
   "K":0.01999816,
   "L":0.00163422,
   "M":0.0006414941,
   "N":0.2331202,
   "P":0.01900553,
   "Q":0.02334345,
   "R":0.007586211,
   "S":0.08303903999999999,
   "T":0.02482688,
   "V":0.002872194,
   "W":0.0003945694,
   "Y":0.005048546
  },
 "E":{
   "F":0.0008705765,
   "G":0.02188286,
   "H":0.01036699,
   "I":0.003010126,
   "K":0.1277224,
   "L":0.007552471,
   "M":0.004362376,
   "N":0.02487791,
   "P":0.02020781,
   "Q":0.1841385,
   "R":0.02227563,
   "S":0.04097288,
   "T":0.03524391,
   "V":0.01853676,
   "W":0.001027702,
   "Y":0.004485425
  },
 "F":{
   "G":0.00561965,
   "H":0.01668329,
   "I":0.0756681,
   "K":0.001690408,
   "L":0.2810447,
   "M":0.04516806,
   "N":0.004111401,
   "P":0.004551429,
   "Q":0.001599162,
   "R":0.003226682,
   "S":0.02422454,
   "T":0.009619268,
   "V":0.04952661,
   "W":0.03243575,
   "Y":0.2916085
  },
 "G":{
   "H":0.007618063,
   "I":0.0005919608,
   "K":0.02096479,
   "L":0.00479784,
   "M":0.003503711,
   "N":0.06602329999999999,
   "P":0.009489902,
   "Q":0.01195119,
   "R":0.02388046,
   "S":0.116496,
   "T":0.00756921,
   "V":0.005802412,
   "W":0.003544273,
   "Y":0.002043191
  },
 "H":{
   "I":0.007404237,
   "K":0.04927923,
   "L":0.03970833,
   "M":0.01111019,
   "N":0.207085,
   "P":0.02451727,
   "Q":0.2146863,
   "R":0.1485124,
   "S":0.06628340000000001,
   "T":0.03406144,
   "V":0.009003304,
   "W":0.007881539999999999,
   "Y":0.1983005
  },
 "I":{
   "K":0.01124222,
   "L":0.4493203,
   "M":0.1073075,
   "N":0.008794702,
   "P":0.003771706,
   "Q":0.003249348,
   "R":0.007772081,
   "S":0.004291965,
   "T":0.06026516,
   "V":0.04910478,
   "W":0.001473992,
   "Y":0.008688692
  },
 "K":{
   "L":0.01490483,
   "M":0.01648691,
   "N":0.09851189,
   "P":0.01880635,
   "Q":0.1442522,
   "R":0.3871668,
   "S":0.05012591,
   "T":0.06627711,
   "V":0.01401048,
   "W":0.0006587949,
   "Y":0.004929905
  },
 "L":{
   "M":0.1584993,
   "N":0.003142484,
   "P":0.01200011,
   "Q":0.02597806,
   "R":0.01847365,
   "S":0.0122045,
   "T":0.01766063,
   "V":0.1288122,
   "W":0.008179586000000001,
   "Y":0.01119695
  },
 "M":{
   "N":0.01703821,
   "P":0.004810887,
   "Q":0.07459796,
   "R":0.02962982,
   "S":0.0232297,
   "T":0.1177837,
   "V":0.1436375,
   "W":0.009190009000000001,
   "Y":0.01798497
  },
 "N":{
   "P":0.007795161,
   "Q":0.07563193999999999,
   "R":0.04601631,
   "S":0.268368,
   "T":0.116636,
   "V":0.006330976,
   "W":0.0005989957,
   "Y":0.02286955
  },
 "P":{
   "Q":0.02784403,
   "R":0.02035162,
   "S":0.08959077,
   "T":0.03331558,
   "V":0.02243022,
   "W":0.001255797,
   "Y":0.00334857
  },
 "Q":{
   "R":0.1718491,
   "S":0.08193787,
   "T":0.06297003,
   "V":0.01591156,
   "W":0.003117996,
   "Y":0.009615879000000001
  },
 "R":{
   "S":0.05745502,
   "T":0.03375392,
   "V":0.01292756,
   "W":0.007836037000000001,
   "Y":0.01174968
  },
 "S":{
   "T":0.3773224,
   "V":0.00744159,
   "W":0.003285156,
   "Y":0.01496724
  },
 "T":{
   "V":0.1655336,
   "W":0.00185899,
   "Y":0.009186346
  },
 "V":{
   "W":0.002501667,
   "Y":0.009316084000000001
  },
 "W":{
   "Y":0.1177739
  }
};
//==================================================================================================================//





/**
 * @name models.protein.JTT.frequencies
 * @param {Dictionary} Baseline JTT model
 * @returns {Dictionary} Updated JTT model with empirical frequencies
 * @description Define the empirical amino acid frequencies associated with the JTT model of protein evolution
 */
function models.protein.JTT.frequencies (model, namespace, datafilter) {
    model[utility.getGlobalValue("terms.efv_estimate")] =
      {{	0.07686099999999986}
        {	0.020279}
        {	0.051269}
        {	0.06182}
        {	0.04053}
        {	0.074714}
        {	0.022983}
        {	0.052569}
        {	0.059498}
        {	0.091111}
        {	0.023414}
        {	0.042546}
        {	0.050532}
        {	0.041061}
        {	0.051057}
        {	0.068225}
        {	0.058518}
        {	0.066374}
        {	0.014336}
        {   0.032303}
    };

    model[utility.getGlobalValue("terms.model.efv_estimate_name")] = utility.getGlobalValue("terms.frequencies.predefined");
    (model[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.model.empirical")] = 0;
    return model;
}



/* Define a dictionary of amino-acid exchangeability rates for the JTT model of protein evolution.  */ 
models.protein.JTT.Rij = {
 "A":{
   "C":0.01164983,
   "D":0.04242226,
   "E":0.06594220000000001,
   "F":0.005605013,
   "G":0.1300142,
   "H":0.005055569,
   "I":0.01901336,
   "K":0.02198075,
   "L":0.02824504,
   "M":0.01099041,
   "N":0.02373925,
   "P":0.09902242999999999,
   "Q":0.02285967,
   "R":0.02714587,
   "S":0.2651969,
   "T":0.2681624,
   "V":0.1940882,
   "W":0.00120894,
   "Y":0.004506008
  },
 "C":{
   "D":0.005415286,
   "E":0.003332529,
   "F":0.02749291,
   "G":0.04082289,
   "H":0.01666262,
   "I":0.007914733,
   "K":0.002915936,
   "L":0.01499622,
   "M":0.009581053000000001,
   "N":0.01333012,
   "P":0.006248431,
   "Q":0.003749032,
   "R":0.05207011,
   "S":0.1470474,
   "T":0.02749309,
   "V":0.04123968,
   "W":0.01582953,
   "Y":0.06831603999999999
  },
 "D":{
   "E":0.4801284,
   "F":0.001318116,
   "G":0.0950686,
   "H":0.0237263,
   "I":0.00609632,
   "K":0.01680615,
   "L":0.005602049,
   "M":0.004448682,
   "N":0.2361102,
   "P":0.006425849,
   "Q":0.0214193,
   "R":0.007908676,
   "S":0.04020279,
   "T":0.02487944,
   "V":0.02092512,
   "W":0.0008238323,
   "Y":0.014664
  },
 "E":{
   "F":0.001776388,
   "G":0.08335330000000001,
   "H":0.005602518,
   "I":0.005875793,
   "K":0.1030317,
   "L":0.008881953,
   "M":0.004099415,
   "N":0.02459647,
   "P":0.009701838000000001,
   "Q":0.1403343,
   "R":0.01626078,
   "S":0.02131682,
   "T":0.01940362,
   "V":0.03088188,
   "W":0.001639765,
   "Y":0.002049689
  },
 "F":{
   "G":0.003751538,
   "H":0.01042113,
   "I":0.04085083,
   "K":0.00145895,
   "L":0.2278042,
   "M":0.01021273,
   "N":0.003126321,
   "P":0.00750314,
   "Q":0.001875789,
   "R":0.003334736,
   "S":0.0644024,
   "T":0.008128382,
   "V":0.03939149,
   "W":0.007711647,
   "Y":0.1771572
  },
 "G":{
   "H":0.004635577,
   "I":0.002826581,
   "K":0.01605493,
   "L":0.006331483,
   "M":0.003052693,
   "N":0.03290136,
   "P":0.01051474,
   "Q":0.009497159999999999,
   "R":0.06941973,
   "S":0.1278738,
   "T":0.01854212,
   "V":0.03120506,
   "W":0.007801362,
   "Y":0.001695907
  },
 "H":{
   "I":0.00955641,
   "K":0.03124215,
   "L":0.04925195,
   "M":0.007718657,
   "N":0.1712807,
   "P":0.05770555,
   "Q":0.2333939,
   "R":0.1639271,
   "S":0.0507224,
   "T":0.02793385,
   "V":0.008086143,
   "W":0.001837774,
   "Y":0.1889208
  },
 "I":{
   "K":0.01205203,
   "L":0.2127567,
   "M":0.1131285,
   "N":0.0208902,
   "P":0.004981443,
   "Q":0.003213843,
   "R":0.01221257,
   "S":0.02763923,
   "T":0.1494435,
   "V":0.6328057,
   "W":0.001928334,
   "Y":0.009802181
  },
 "K":{
   "L":0.01334602,
   "M":0.01462393,
   "N":0.1076208,
   "P":0.01093234,
   "Q":0.1218169,
   "R":0.333364,
   "S":0.03237125,
   "T":0.05650736,
   "V":0.008234754,
   "W":0.001277824,
   "Y":0.002839562
  },
 "L":{
   "M":0.09030557,
   "N":0.005841096,
   "P":0.05358937,
   "Q":0.0291124,
   "R":0.01900652,
   "S":0.04042405,
   "T":0.01594697,
   "V":0.1169137,
   "W":0.007602722,
   "Y":0.007788057
  },
 "M":{
   "N":0.0140708,
   "P":0.008298109999999999,
   "Q":0.0187608,
   "R":0.02200785,
   "S":0.01948259,
   "T":0.1237496,
   "V":0.2016795,
   "W":0.002886323,
   "Y":0.006133368
  },
 "N":{
   "P":0.006154998,
   "Q":0.03156908,
   "R":0.02303155,
   "S":0.3450795,
   "T":0.1375939,
   "V":0.01092017,
   "W":0.000397107,
   "Y":0.02263447
  },
 "P":{
   "Q":0.06603124,
   "R":0.03627542,
   "S":0.1902389,
   "T":0.06887338,
   "V":0.01404214,
   "W":0.001003017,
   "Y":0.003677695
  },
 "Q":{
   "R":0.1542939,
   "S":0.03744234,
   "T":0.03065318,
   "V":0.01193211,
   "W":0.002468744,
   "Y":0.008229024999999999
  },
 "R":{
   "S":0.06833079,
   "T":0.03805319,
   "V":0.01141599,
   "W":0.01803412,
   "Y":0.007610617
  },
 "S":{
   "T":0.2795782,
   "V":0.02711589,
   "W":0.004457448,
   "Y":0.02030591
  },
 "T":{
   "V":0.0759305,
   "W":0.00115485,
   "Y":0.006495937
  },
 "V":{
   "W":0.003436295,
   "Y":0.005345272
  },
 "W":{
   "Y":0.02415905
  }
};

//==================================================================================================================//



/**
 * @name models.protein.JC69.frequencies
 * @param {Dictionary} Baseline JC69 model
 * @returns {Dictionary} Updated JC69 model with empirical frequencies
 * @description Define the empirical amino acid frequencies associated with the JC69 model of protein evolution
 */
lfunction models.protein.JC69.frequencies (model, namespace, datafilter) {
    model[utility.getGlobalValue("terms.efv_estimate")] =
      {{	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {   0.05}
    };

    model[utility.getGlobalValue("terms.model.efv_estimate_name")] = utility.getGlobalValue("terms.frequencies.predefined");
    (model[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.model.empirical")] = 0;
    return model;
}

/* Define a dictionary of equal amino-acid exchangeability rates for the JC69 model of protein evolution. */ 
models.protein.JC69.Rij = {
    "A":
        {"C":0.05,
        "D":0.05,
        "E":0.05,
        "F":0.05,
        "G":0.05,
        "H":0.05,
        "I":0.05,
        "K":0.05,
        "L":0.05,
        "M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "C":
        {"D":0.05,
        "E":0.05,
        "F":0.05,
        "G":0.05,
        "H":0.05,
        "I":0.05,
        "K":0.05,
        "L":0.05,
        "M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "D":
        {"E":0.05,
        "F":0.05,
        "G":0.05,
        "H":0.05,
        "I":0.05,
        "K":0.05,
        "L":0.05,
        "M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "E":
        {"F":0.05,
        "G":0.05,
        "H":0.05,
        "I":0.05,
        "K":0.05,
        "L":0.05,
        "M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "F":
        {"G":0.05,
        "H":0.05,
        "I":0.05,
        "K":0.05,
        "L":0.05,
        "M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "G":
        {"H":0.05,
        "I":0.05,
        "K":0.05,
        "L":0.05,
        "M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "H":
        {"I":0.05,
        "K":0.05,
        "L":0.05,
        "M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "I":
        {"K":0.05,
        "L":0.05,
        "M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":30.05,
        "W":0.05,
        "Y":0.05},
    "K":
        {"L":0.05,
        "M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "L":
        {"M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "M":
        {"N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "N":
        {"P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "P":
        {"Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "Q":
        {"R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "R":
        {"S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "S":
        {"T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "T":
        {"V":0.05,
        "W":0.05,
        "Y":0.05},
    "V":
        {"W":0.05,
        "Y":0.05},
    "W":
        {"Y":0.05}
};



//==================================================================================================================//



/**
 * @name models.protein.mtMAM.frequencies
 * @param {Dictionary} Baseline mtMAM model
 * @returns {Dictionary} Updated mtMAM model with empirical frequencies
 * @description Define the empirical amino acid frequencies associated with the mtMAM model of protein evolution
 */
lfunction models.protein.mtMAM.frequencies (model, namespace, datafilter) {
    model[utility.getGlobalValue("terms.efv_estimate")] =
       {{	0.0692}
        {   0.0065}
        {	0.0186}
        {	0.0236}
        {	0.0611}
        {	0.0557}
        {	0.0277}
        {	0.0905}
        {	0.0221}
        {	0.1675}
        {	0.0561}
        {	0.04}
        {	0.0536}
        {	0.0238}
        {	0.0184}
        {	0.0725}
        {	0.087}
        {	0.0428}
        {	0.0293}
        {	0.034}
    };

    model[utility.getGlobalValue("terms.model.efv_estimate_name")] = utility.getGlobalValue("terms.frequencies.predefined");
    (model[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.model.empirical")] = 0;
    return model;
}

/* Define a dictionary of equal amino-acid exchangeability rates for the mtMAM model of protein evolution.  */ 
models.protein.mtMAM.Rij = {
    "A": 
          {"C": 0.0,
           "D": 0.11,
           "E": 0.0,
           "F": 0.0,
           "G": 0.78,
           "H": 0.08,
           "I": 0.75,
           "K": 0.0,
           "L": 0.21,
           "M": 0.76,
           "N": 0.02,
           "P": 0.53,
           "Q": 0.0,
           "R": 0.32,
           "S": 3.42,
           "T": 6.81,
           "V": 3.98,
           "W": 0.05,
           "Y": 0.0},
     "C": 
          {"D": 0.0,
           "E": 0.0,
           "F": 0.07,
           "G": 0.0,
           "H": 3.05,
           "I": 0.41,
           "K": 0.0,
           "L": 0.27,
           "M": 0.0,
           "N": 0.0,
           "P": 0.0,
           "Q": 0.0,
           "R": 1.86,
           "S": 3.47,
           "T": 1.14,
           "V": 0.0,
           "W": 0.65,
           "Y": 5.3},
     "D": 
          {"E": 5.69,
           "F": 0.05,
           "G": 0.79,
           "H": 0.11,
           "I": 0.0,
           "K": 0.0,
           "L": 0.0,
           "M": 0.0,
           "N": 8.64,
           "P": 0.02,
           "Q": 0.49,
           "R": 0.0,
           "S": 0.16,
           "T": 0.0,
           "V": 0.1,
           "W": 0.0,
           "Y": 0.0},
     "E": 
          {"F": 0.0,
           "G": 0.22,
           "H": 0.22,
           "I": 0.0,
           "K": 2.15,
           "L": 0.0,
           "M": 0.0,
           "N": 0.0,
           "P": 0.0,
           "Q": 2.74,
           "R": 0.0,
           "S": 0.21,
           "T": 0.04,
           "V": 0.2,
           "W": 0.0,
           "Y": 0.0},
     "F": 
          {"G": 0.0,
           "H": 0.0,
           "I": 0.57,
           "K": 0.0,
           "L": 2.46,
           "M": 0.11,
           "N": 0.06,
           "P": 0.17,
           "Q": 0.0,
           "R": 0.0,
           "S": 0.9,
           "T": 0.08,
           "V": 0.06,
           "W": 0.0,
           "Y": 6.82},
     "G": 
          {"H": 0.0,
           "I": 0.0,
           "K": 0.0,
           "L": 0.0,
           "M": 0.0,
           "N": 0.47,
           "P": 0.0,
           "Q": 0.0,
           "R": 0.18,
           "S": 1.12,
           "T": 0.0,
           "V": 0.05,
           "W": 0.0,
           "Y": 0.01},
     "H": 
          {"I": 0.0,
           "K": 0.0,
           "L": 0.26,
           "M": 0.0,
           "N": 4.58,
           "P": 0.53,
           "Q": 5.5,
           "R": 2.32,
           "S": 0.2,
           "T": 0.01,
           "V": 0.0,
           "W": 0.0,
           "Y": 15.25},
     "I": 
          {"K": 0.06,
           "L": 2.32,
           "M": 3.78,
           "N": 0.19,
           "P": 0.05,
           "Q": 0.0,
           "R": 0.0,
           "S": 0.0,
           "T": 3.6,
           "V": 22.2,
           "W": 0.0,
           "Y": 0.16},
     "K": {"L": 0.04,
           "M": 0.59,
           "N": 4.08,
           "P": 0.18,
           "Q": 2.42,
           "R": 0.5,
           "S": 0.65,
           "T": 0.5,
           "V": 0.0,
           "W": 0.0,
           "Y": 0.67},
     "L": 
          {"M": 6.09,
           "N": 0.0,
           "P": 0.43,
           "Q": 0.2,
           "R": 0.06,
           "S": 0.74,
           "T": 0.34,
           "V": 1.0,
           "W": 0.12,
           "Y": 0.25},
     "M": 
          {"N": 0.21,
           "P": 0.0,
           "Q": 0.22,
           "R": 0.0,
           "S": 0.47,
           "T": 6.91,
           "V": 8.32,
           "W": 0.13,
           "Y": 0.0},
     "N": 
          {"P": 0.33,
           "Q": 0.08,
           "R": 0.04,
           "S": 4.46,
           "T": 1.1,
           "V": 0.0,
           "W": 0.06,
           "Y": 1.56},
     "P": 
          {"Q": 0.51,
           "R": 0.09,
           "S": 2.02,
           "T": 0.78,
           "V": 0.0,
           "W": 0.07,
           "Y": 0.08},
     "Q": 
          {"R": 2.46, 
           "S": 0.3, 
           "T": 0.0, 
           "V": 0.33, 
           "W": 0.0, 
           "Y": 0.54},
     "R": 
          {"S": 0.03,
           "T": 0.0, 
           "V": 0.0, 
           "W": 0.16, 
           "Y": 0.0},
     "S": 
          {"T": 6.14, 
           "V": 0.0, 
           "W": 0.17, 
           "Y": 1.07},
     "T": 
          {"V": 2.37, 
           "W": 0.0, 
           "Y": 0.0},
     "V": 
        {"W": 0.0, 
         "Y": 0.0},
     "W": 
         {"Y": 0.14}
};


lfunction models.protein.HIVWm.frequencies (model, namespace, datafilter) {
    model[utility.getGlobalValue("terms.efv_estimate")] =
      {{0.0377494}
        {0.0240105}
        {0.0342034}
        {0.0618606}
        {0.0422741}
        {0.0838496}
        {0.0156076}
        {0.0983641}
        {0.0641682}
        {0.0577867}
        {0.0158419}
        {0.0891129}
        {0.0458601}
        {0.0437824}
        {0.057321}
        {0.0550846}
        {0.0813774}
        {0.0515639}
        {0.019597}
        {0.0205847}};
    model[utility.getGlobalValue("terms.model.efv_estimate_name")] = utility.getGlobalValue("terms.frequencies.predefined");
    (model[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.model.empirical")] = 0;
    return model;
}
        
models.protein.HIVWm.Rij = {
 "A":{
   "C":0.167653,
   "D":4.43521,
   "E":5.56325,
   "F":0.597923,
   "G":1.8685,
   "H":0.005,
   "I":0.005,
   "K":0.592784,
   "L":0.16024,
   "M":0.005,
   "N":0.617509,
   "P":1.00981,
   "Q":0.005,
   "R":0.0744808,
   "S":8.594200000000001,
   "T":24.1422,
   "V":24.8094,
   "W":0.005,
   "Y":0.005
  },
 "C":{
   "D":0.005,
   "E":0.005,
   "F":0.362959,
   "G":0.0489798,
   "H":0.005,
   "I":0.005,
   "K":0.005,
   "L":0.005,
   "M":0.005,
   "N":0.0604932,
   "P":0.005,
   "Q":0.005,
   "R":2.86364,
   "S":1.12195,
   "T":0.005,
   "V":0.005,
   "W":5.49894,
   "Y":8.34835
  },
 "D":{
   "E":12.1233,
   "F":0.005,
   "G":10.3969,
   "H":2.31779,
   "I":0.145124,
   "K":0.894313,
   "L":0.005,
   "M":0.005,
   "N":29.4087,
   "P":0.005,
   "Q":0.005,
   "R":0.0674539,
   "S":0.427881,
   "T":0.630395,
   "V":2.91786,
   "W":0.005,
   "Y":2.28154
  },
 "E":{
   "F":0.005,
   "G":14.7801,
   "H":0.005,
   "I":0.0390512,
   "K":23.9626,
   "L":0.129839,
   "M":0.005,
   "N":0.201526,
   "P":0.005,
   "Q":3.20656,
   "R":0.0251632,
   "S":0.005,
   "T":0.458743,
   "V":2.19952,
   "W":0.005,
   "Y":0.005
  },
 "F":{
   "G":0.005,
   "H":0.005,
   "I":1.48288,
   "K":0.005,
   "L":7.48781,
   "M":0.005,
   "N":0.005,
   "P":0.0342252,
   "Q":0.005,
   "R":0.005,
   "S":4.27939,
   "T":0.114512,
   "V":2.28,
   "W":0.005,
   "Y":4.12728
  },
 "G":{
   "H":0.005,
   "I":0.005,
   "K":0.279425,
   "L":0.0489798,
   "M":0.0489798,
   "N":0.0604932,
   "P":0.005,
   "Q":0.0604932,
   "R":13.4379,
   "S":6.27966,
   "T":0.0489798,
   "V":2.79622,
   "W":2.8258,
   "Y":0.005
  },
 "H":{
   "I":0.005,
   "K":0.22406,
   "L":1.76382,
   "M":0.005,
   "N":8.59876,
   "P":13.9444,
   "Q":18.5465,
   "R":6.84405,
   "S":0.7251570000000001,
   "T":0.95956,
   "V":0.827479,
   "W":0.005,
   "Y":47.4889
  },
 "I":{
   "K":0.817481,
   "L":9.102460000000001,
   "M":17.3064,
   "N":0.987028,
   "P":0.005,
   "Q":0.0342252,
   "R":1.34069,
   "S":0.7400910000000001,
   "T":9.36345,
   "V":24.8231,
   "W":0.005,
   "Y":0.114512
  },
 "K":{
   "L":0.005,
   "M":4.09564,
   "N":10.6655,
   "P":0.111928,
   "Q":13.0705,
   "R":39.8897,
   "S":0.005,
   "T":4.04802,
   "V":0.128065,
   "W":0.005,
   "Y":0.005
  },
 "L":{
   "M":11.3839,
   "N":0.005,
   "P":9.83095,
   "Q":2.89048,
   "R":0.586757,
   "S":6.14396,
   "T":0.005,
   "V":2.95344,
   "W":1.37031,
   "Y":0.005
  },
 "M":{
   "N":0.201526,
   "P":0.005,
   "Q":0.005,
   "R":3.28652,
   "S":0.392575,
   "T":7.41313,
   "V":14.7683,
   "W":0.005,
   "Y":0.579198
  },
 "N":{
   "P":0.344848,
   "Q":0.342068,
   "R":0.16024,
   "S":14.5699,
   "T":4.54206,
   "V":0.0744808,
   "W":0.005,
   "Y":5.06475
  },
 "P":{
   "Q":3.04502,
   "R":0.404723,
   "S":14.249,
   "T":4.33701,
   "V":0.005,
   "W":0.005,
   "Y":0.005
  },
 "Q":{
   "R":10.6746,
   "S":0.16024,
   "T":0.203091,
   "V":0.005,
   "W":0.0443298,
   "Y":0.005
  },
 "R":{
   "S":8.350239999999999,
   "T":0.928203,
   "V":0.279425,
   "W":5.96564,
   "Y":0.005
  },
 "S":{
   "T":6.34079,
   "V":0.862637,
   "W":1.10156,
   "Y":0.933142
  },
 "T":{
   "V":0.005,
   "W":0.005,
   "Y":0.490608
  },
 "V":{
   "W":0.005,
   "Y":1.35482
  },
 "W":{
   "Y":0.005
  }
};


lfunction models.protein.HIVBm.frequencies (model, namespace, datafilter) {
    model[utility.getGlobalValue("terms.efv_estimate")] =
      {{   0.060490222}
        {   0.020075899}
        {   0.042109048}
        {   0.071567447}
        {   0.028809447}
        {   0.072308239}
        {   0.022293943}
        {   0.069730629}
        {   0.056968211}
        {   0.098851122}
        {   0.019768318}
        {   0.044127815}
        {   0.046025282}
        {   0.053606488}
        {   0.066039665}
        {    0.05060433}
        {   0.053636813}
        {   0.061625237}
        {   0.033011601}
        {   0.028350243}};
    model[utility.getGlobalValue("terms.model.efv_estimate_name")] = utility.getGlobalValue("terms.frequencies.predefined");
    (model[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.model.empirical")] = 0;
    return model;
}
          
models.protein.HIVBm.Rij = {
 "A":{
   "C":0.123758,
   "D":1.45504,
   "E":1.48135,
   "F":0.0141269,
   "G":2.13536,
   "H":0.0847613,
   "I":0.005,
   "K":0.005,
   "L":0.215256,
   "M":0.0186643,
   "N":0.005,
   "P":2.12217,
   "Q":0.0551128,
   "R":0.307507,
   "S":2.46633,
   "T":15.9183,
   "V":7.61428,
   "W":0.005,
   "Y":0.005
  },
 "C":{
   "D":0.005,
   "E":0.005,
   "F":9.29815,
   "G":0.897871,
   "H":0.240073,
   "I":0.005,
   "K":0.005,
   "L":0.129777,
   "M":0.005,
   "N":0.08606419999999999,
   "P":0.005,
   "Q":0.005,
   "R":0.351721,
   "S":4.69314,
   "T":0.739969,
   "V":0.420027,
   "W":2.63277,
   "Y":7.57932
  },
 "D":{
   "E":10.5872,
   "F":0.005,
   "G":2.83806,
   "H":1.9169,
   "I":0.0176792,
   "K":0.005,
   "L":0.008760479999999999,
   "M":0.005,
   "N":17.6612,
   "P":0.0342658,
   "Q":0.005,
   "R":0.005,
   "S":0.52823,
   "T":0.274724,
   "V":1.04793,
   "W":0.005,
   "Y":0.6746529999999999
  },
 "E":{
   "F":0.005,
   "G":3.92775,
   "H":0.11974,
   "I":0.00609079,
   "K":4.61482,
   "L":0.005,
   "M":0.175789,
   "N":0.07926329999999999,
   "P":0.0120226,
   "Q":2.5602,
   "R":0.0749218,
   "S":0.005,
   "T":0.289774,
   "V":1.02847,
   "W":0.005,
   "Y":0.07926329999999999
  },
 "F":{
   "G":0.291561,
   "H":0.145558,
   "I":3.39836,
   "K":0.0342658,
   "L":8.524839999999999,
   "M":0.188025,
   "N":0.005,
   "P":0.005,
   "Q":0.005,
   "R":0.005,
   "S":0.956472,
   "T":0.0141269,
   "V":0.723274,
   "W":0.8293430000000001,
   "Y":15.34
  },
 "G":{
   "H":0.005,
   "I":0.005,
   "K":0.521705,
   "L":0.005,
   "M":0.005,
   "N":0.323401,
   "P":0.005,
   "Q":0.0619137,
   "R":3.65345,
   "S":4.38041,
   "T":0.369615,
   "V":0.953155,
   "W":1.21674,
   "Y":0.005
  },
 "H":{
   "I":0.103111,
   "K":0.005,
   "L":1.74171,
   "M":0.005,
   "N":7.64585,
   "P":2.45318,
   "Q":7.05545,
   "R":9.04044,
   "S":0.382747,
   "T":0.7115939999999999,
   "V":0.005,
   "W":0.06951789999999999,
   "Y":18.6943
  },
 "I":{
   "K":0.322319,
   "L":5.95879,
   "M":11.2065,
   "N":0.680565,
   "P":0.0410593,
   "Q":0.005,
   "R":0.677289,
   "S":1.21803,
   "T":8.612170000000001,
   "V":17.7389,
   "W":0.005,
   "Y":0.148168
  },
 "K":{
   "L":0.0814995,
   "M":1.28246,
   "N":7.90443,
   "P":0.0313862,
   "Q":6.54737,
   "R":20.45,
   "S":0.504111,
   "T":4.67142,
   "V":0.265829,
   "W":0.005,
   "Y":0.005
  },
 "L":{
   "M":5.31961,
   "N":0.005,
   "P":2.07757,
   "Q":1.49456,
   "R":0.701427,
   "S":0.927656,
   "T":0.0437673,
   "V":1.41036,
   "W":0.748843,
   "Y":0.111986
  },
 "M":{
   "N":0.005,
   "P":0.005,
   "Q":0.303676,
   "R":2.51394,
   "S":0.005,
   "T":4.94026,
   "V":6.8532,
   "W":0.089078,
   "Y":0.005
  },
 "N":{
   "P":0.00739578,
   "Q":0.672052,
   "R":0.295543,
   "S":13.1447,
   "T":6.88667,
   "V":0.026656,
   "W":0.005,
   "Y":1.76417
  },
 "P":{
   "Q":4.47211,
   "R":1.28355,
   "S":5.37762,
   "T":2.01417,
   "V":0.005,
   "W":0.0444506,
   "Y":0.0304381
  },
 "Q":{
   "R":3.4215,
   "S":0.116311,
   "T":0.243589,
   "V":0.0209153,
   "W":0.026656,
   "Y":0.113033
  },
 "R":{
   "S":3.4791,
   "T":2.86868,
   "V":0.0812454,
   "W":0.9913380000000001,
   "Y":0.00991826
  },
 "S":{
   "T":8.93107,
   "V":0.0749218,
   "W":0.0248728,
   "Y":0.648024
  },
 "T":{
   "V":0.709226,
   "W":0.005,
   "Y":0.105652
  },
 "V":{
   "W":0.005,
   "Y":0.0410593
  },
 "W":{
   "Y":1.28022
  }
};

lfunction models.protein.cpREV.frequencies (model, namespace, datafilter) {
    model[utility.getGlobalValue("terms.efv_estimate")] =
      {{         0.076}
        {         0.009}
        {         0.037}
        {         0.049}
        {         0.051}
        {         0.084}
        {         0.025}
        {         0.081}
        {          0.05}
        {         0.101}
        {         0.022}
        {         0.041}
        {         0.043}
        {         0.038}
        {         0.062}
        {         0.062}
        {         0.054}
        {         0.066}
        {         0.018}
        {         0.031}};
    model[utility.getGlobalValue("terms.model.efv_estimate_name")] = utility.getGlobalValue("terms.frequencies.predefined");
    (model[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.model.empirical")] = 0;
    return model;
}
        
models.protein.cpREV.Rij = {
 "A":{
   "C":1.2816975,
   "D":0.33527215,
   "E":0.95600458,
   "F":0.13027718,
   "G":1.2740342,
   "H":0.1264455,
   "I":0.27779692,
   "K":0.45213844,
   "L":0.37742065,
   "M":0.35443056,
   "N":0.43489587,
   "P":0.93876201,
   "Q":0.30078701,
   "R":0.20116329,
   "S":4.6746517,
   "T":2.5672267,
   "V":1.8545339,
   "W":0.026821772,
   "Y":0.10728709
  },
 "C":{
   "D":0.019158408,
   "E":0.019158408,
   "F":1.3909005,
   "G":0.58049978,
   "H":0.84488581,
   "I":0.53643544,
   "K":0.091960361,
   "L":0.75867297,
   "M":0.30461869,
   "N":1.0307224,
   "P":0.54601464,
   "Q":0.019158408,
   "R":1.576737,
   "S":4.465825,
   "T":1.1035243,
   "V":1.1341778,
   "W":0.8333907699999999,
   "Y":2.8086227
  },
 "D":{
   "E":7.0713686,
   "F":0.042148499,
   "G":0.8257274,
   "H":0.63414332,
   "I":0.019158408,
   "K":0.78932643,
   "L":0.019158408,
   "M":0.09004452,
   "N":8.4967541,
   "P":0.32569294,
   "Q":0.7663363399999999,
   "R":0.082381156,
   "S":1.1303461,
   "T":0.50961366,
   "V":0.14368806,
   "W":0.034485135,
   "Y":0.53835128
  },
 "E":{
   "F":0.27779692,
   "G":0.72610368,
   "H":0.31036622,
   "I":0.28354444,
   "K":5.0367456,
   "L":0.15709895,
   "M":0.21649002,
   "N":2.0212121,
   "P":0.35443056,
   "Q":5.9812551,
   "R":0.29120781,
   "S":1.0881976,
   "T":0.70694527,
   "V":0.38316817,
   "W":0.12069797,
   "Y":0.2720494
  },
 "F":{
   "G":0.047896021,
   "H":0.24331179,
   "I":0.86979174,
   "K":0.13794054,
   "L":2.4292862,
   "M":0.6264799599999999,
   "N":0.18583656,
   "P":0.082381156,
   "Q":0.019158408,
   "R":0.10153956,
   "S":0.93301449,
   "T":0.28354444,
   "V":0.60732155,
   "W":0.89661351,
   "Y":4.5405428
  },
 "G":{
   "H":0.036400976,
   "I":0.07663363400000001,
   "K":0.50386614,
   "L":0.038316817,
   "M":0.040232658,
   "N":1.2510441,
   "P":0.053643544,
   "Q":0.25480683,
   "R":0.46554933,
   "S":1.323846,
   "T":0.17625736,
   "V":0.17434152,
   "W":0.15709895,
   "Y":0.019158408
  },
 "H":{
   "I":0.055559384,
   "K":0.58433146,
   "L":0.1264455,
   "M":0.019158408,
   "N":2.6917564,
   "P":0.29120781,
   "Q":2.431202,
   "R":1.3698262,
   "S":0.58049978,
   "T":0.061306907,
   "V":0.047896021,
   "W":0.13219302,
   "Y":3.7761223
  },
 "I":{
   "K":0.66096509,
   "L":3.3431423,
   "M":3.39487,
   "N":0.32186126,
   "P":0.22415338,
   "Q":0.17625736,
   "R":0.26055435,
   "S":0.41382162,
   "T":1.9924745,
   "V":9.190288499999999,
   "W":0.080465315,
   "Y":0.17050984
  },
 "K":{
   "L":0.4176533,
   "M":0.36975728,
   "N":4.6554933,
   "P":0.57858393,
   "Q":6.3471807,
   "R":8.586798699999999,
   "S":1.6629499,
   "T":1.7587419,
   "V":0.47704437,
   "W":0.019158408,
   "Y":0.47321269
  },
 "L":{
   "M":2.588301,
   "N":0.21649002,
   "P":0.41956914,
   "Q":0.54793048,
   "R":0.38891569,
   "S":0.98857388,
   "T":0.29887117,
   "V":1.6572023,
   "W":0.30461869,
   "Y":0.36209392
  },
 "M":{
   "N":0.11686629,
   "P":0.19158408,
   "Q":0.38699985,
   "R":0.23948011,
   "S":0.1781732,
   "T":1.2357173,
   "V":0.9100244,
   "W":0.16476231,
   "Y":0.41190578
  },
 "N":{
   "P":0.33144047,
   "Q":1.4713658,
   "R":0.68395518,
   "S":3.9945282,
   "T":2.6687663,
   "V":0.15901479,
   "W":0.07663363400000001,
   "Y":1.444544
  },
 "P":{
   "Q":0.61881659,
   "R":0.16667815,
   "S":2.3028407,
   "T":0.49811862,
   "V":0.23373258,
   "W":0.09387620100000001,
   "Y":0.18583656
  },
 "Q":{
   "R":3.3431423,
   "S":0.75867297,
   "T":0.46171764,
   "V":0.10345541,
   "W":0.10153956,
   "Y":0.74909377
  },
 "R":{
   "S":0.73759872,
   "T":0.60157402,
   "V":0.17625736,
   "W":0.44064339,
   "Y":0.61881659
  },
 "S":{
   "T":4.1209737,
   "V":0.31994542,
   "W":0.13985638,
   "Y":1.0000689
  },
 "T":{
   "V":1.456039,
   "W":0.055559384,
   "Y":0.1360247
  },
 "V":{
   "W":0.019158408,
   "Y":0.22798506
  },
 "W":{
   "Y":0.66288093
  }
};



lfunction models.protein.AB.frequencies (model, namespace, datafilter) {
    model[utility.getGlobalValue("terms.efv_estimate")] =
      {{6.541704e-02}
        {4.708366e-02}
        {3.168984e-02}
        {4.688141e-02}
        {2.150693e-02}
        {4.240711e-02}
        {2.842211e-02}
        {1.005278e-01}
        {9.812606e-03}
        {3.424424e-02}
        {6.222565e-02}
        {4.844488e-02}
        {1.760370e-02}
        {3.478555e-02}
        {3.962469e-02}
        {1.280566e-01}
        {8.199314e-02}
        {3.393045e-02}
        {7.586119e-02}
        {4.948141e-02}};
    model[utility.getGlobalValue("terms.model.efv_estimate_name")] = utility.getGlobalValue("terms.frequencies.predefined");
    (model[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.model.empirical")] = 0;
    return model;
}
        
models.protein.AB.Rij = {
 "A":{
   "C":0.1784266,
   "D":0.09291290000000001,
   "E":1.241095,
   "F":0.008929181,
   "G":0.1992269,
   "H":0.9521821,
   "I":1.851951,
   "K":5.241316,
   "L":0.1140412,
   "M":0.06969101,
   "N":0.07388355000000001,
   "P":0.06299270999999999,
   "Q":0.1130146,
   "R":1.800713,
   "S":0.9988358000000001,
   "T":2.912317,
   "V":0.07939549,
   "W":0.1433528,
   "Y":3.774477
  },
 "C":{
   "D":0.782913,
   "E":0.05795374,
   "F":0.1821885,
   "G":1.923901,
   "H":0.06273863,
   "I":1.0894,
   "K":10.4955,
   "L":0.3245175,
   "M":0.3932002,
   "N":7.54924,
   "P":0.3362326,
   "Q":0.08208677,
   "R":0.3498713,
   "S":1.926435,
   "T":1.135258,
   "V":0.5724286,
   "W":0.1711315,
   "Y":0.1366145
  },
 "D":{
   "E":7.185182,
   "F":1.374268e-06,
   "G":0.08705989,
   "H":0.5038373,
   "I":0.4868901,
   "K":14.05444,
   "L":1.762721,
   "M":0.02769442,
   "N":6.190318,
   "P":0.03972173,
   "Q":0.1955446,
   "R":0.007342554,
   "S":7.348346,
   "T":2.147175,
   "V":0.0007310937,
   "W":2.622763,
   "Y":0.04931206
  },
 "E":{
   "F":0.02340019,
   "G":0.1843856,
   "H":7.426619,
   "I":2.1124,
   "K":11.26995,
   "L":0.03916999,
   "M":0.03020502,
   "N":0.06622772,
   "P":0.03357577,
   "Q":0.1031734,
   "R":0.1509482,
   "S":0.5822988,
   "T":0.1516881,
   "V":0.01423897,
   "W":0.9078338,
   "Y":0.4076074
  },
 "F":{
   "G":1.046446e-08,
   "H":7.519215e-11,
   "I":0.05891123,
   "K":3.963388,
   "L":0.0006594967,
   "M":0.006079219,
   "N":3.722878e-16,
   "P":0.007213178,
   "Q":0.1993818,
   "R":0.004878395,
   "S":0.2639482,
   "T":3.225214e-06,
   "V":0.4440833,
   "W":0.7741612,
   "Y":0.02243512
  },
 "G":{
   "H":3.691671,
   "I":0.0551634,
   "K":8.908434,
   "L":6.712736e-06,
   "M":0.6802781,
   "N":3.030805,
   "P":0.001233336,
   "Q":0.00149661,
   "R":0.7426909,
   "S":0.0005906405,
   "T":0.1202094,
   "V":4.332983e-05,
   "W":0.02737091,
   "Y":0.009047737
  },
 "H":{
   "I":1.38937,
   "K":7.29808,
   "L":0.0001029959,
   "M":0.001283121,
   "N":3.608816,
   "P":0.07659566,
   "Q":0.05288625,
   "R":0.02889815,
   "S":0.06776709,
   "T":0.06016624,
   "V":0.02252612,
   "W":0.1240642,
   "Y":0.5795409
  },
 "I":{
   "K":9.139518000000001,
   "L":0.03560482,
   "M":0.02157936,
   "N":0.055044,
   "P":0.02187264,
   "Q":0.1984772,
   "R":0.07915055999999999,
   "S":0.9984215,
   "T":0.07862767,
   "V":0.1386853,
   "W":0.2295842,
   "Y":0.42282
  },
 "K":{
   "L":4.706586,
   "M":5.879103,
   "N":1.455741,
   "P":2.298295,
   "Q":5.642309,
   "R":10.49496,
   "S":5.439116,
   "T":3.443285,
   "V":7.01389,
   "W":20.55414,
   "Y":6.890244
  },
 "L":{
   "M":1.601123,
   "N":0.5059793,
   "P":10.96748,
   "Q":2.714705,
   "R":0.05016568,
   "S":0.6007607,
   "T":3.087152,
   "V":0.06318748,
   "W":0.2903165,
   "Y":7.926675
  },
 "M":{
   "N":0.02158451,
   "P":5.647985,
   "Q":3.390618,
   "R":1.149931,
   "S":0.1580539,
   "T":0.5702792,
   "V":0.3378544,
   "W":0.152132,
   "Y":3.59531
  },
 "N":{
   "P":1.238634,
   "Q":0.004649035,
   "R":0.009948993999999999,
   "S":0.08688405,
   "T":1.039298,
   "V":0.008024263,
   "W":0.07109973,
   "Y":0.0349344
  },
 "P":{
   "Q":3.94794,
   "R":0.07417279,
   "S":0.01861354,
   "T":1.415612,
   "V":0.1011149,
   "W":0.002246759,
   "Y":4.39672
  },
 "Q":{
   "R":0.3556198,
   "S":0.9813064,
   "T":0.03674486,
   "V":0.2199856,
   "W":7.074464,
   "Y":1.643946
  },
 "R":{
   "S":1.284651,
   "T":0.9057112,
   "V":0.005516074,
   "W":0.1992133,
   "Y":0.2217442
  },
 "S":{
   "T":3.058575,
   "V":0.1385142,
   "W":0.8104751,
   "Y":0.07477041
  },
 "T":{
   "V":0.01412361,
   "W":0.09984255,
   "Y":0.2166054
  },
 "V":{
   "W":0.6121284,
   "Y":0.09663569
  },
 "W":{
   "Y":0.5010635
  }
};




lfunction models.protein.mtMet.frequencies (model, namespace, datafilter) {
    model[utility.getGlobalValue("terms.efv_estimate")] =
      {{0.043793200}
        {0.012957800}
        {0.057001300}
        {0.016899000}
        {0.011330500}
        {0.018018100}
        {0.022538500}
        {0.047050100}
        {0.017183700}
        {0.089779400}
        {0.155226000}
        {0.039913500}
        {0.067444300}
        {0.088448000}
        {0.037528200}
        {0.093752200}
        {0.063579000}
        {0.022671300}
        {0.041568200}
        {0.053317400}
       };
    model[utility.getGlobalValue("terms.model.efv_estimate_name")] = utility.getGlobalValue("terms.frequencies.predefined");
    (model[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.model.empirical")] = 0;
    return model;
}

models.protein.mtMet.Rij = {
'A': {'C': 0.63325584800000001,
       'D': 0.11915685500000001,
       'E': 0.17916388799999999,
       'F': 0.062762255000000003,
       'G': 1.4658622800000001,
       'H': 0.030192130000000001,
       'I': 0.367600449,
       'K': 0.020509507999999999,
       'L': 0.109872766,
       'M': 0.65336399300000003,
       'N': 0.03289392,
       'P': 0.40807705300000002,
       'Q': 0.052454947000000002,
       'R': 0.058078194999999999,
       'S': 2.7716860149999998,
       'T': 6.7308851599999997,
       'V': 2.815163085,
       'W': 0.013623415999999999,
       'Y': 0.014501406999999999},
 'C': {'D': 0.077419373999999999,
       'E': 0.050609064000000002,
       'F': 0.92581086400000001,
       'G': 0.91412559000000004,
       'H': 0.60383390000000003,
       'I': 0.235804245,
       'K': 0.015753762000000001,
       'L': 0.29951899700000001,
       'M': 0.49234014399999998,
       'N': 0.29328103,
       'P': 0.029408410999999999,
       'Q': 0.15259520800000001,
       'R': 0.73981385700000002,
       'S': 3.2830148709999998,
       'T': 0.33866819599999998,
       'V': 1.3945280440000001,
       'W': 1.018410485,
       'Y': 1.967371255},
 'D': {'E': 6.0337889819999999,
       'F': 0.012560742999999999,
       'G': 0.63075329899999999,
       'H': 0.47979111200000002,
       'I': 0.010668856000000001,
       'K': 0.049007455999999998,
       'L': 0.0055291439999999997,
       'M': 0.026109947000000001,
       'N': 4.6584200710000001,
       'P': 0.044609562999999998,
       'Q': 0.13135570199999999,
       'R': 0.049700411999999999,
       'S': 0.36080478100000002,
       'T': 0.102136221,
       'V': 0.084589028999999996,
       'W': 0.040920527999999998,
       'Y': 0.16028995800000001},
 'E': {'F': 0.017716308,
       'G': 0.76885329499999999,
       'H': 0.105414735,
       'I': 0.014004526,
       'K': 1.3792177830000001,
       'L': 0.019157619000000001,
       'M': 0.128410054,
       'N': 0.81224112400000004,
       'P': 0.048786298999999998,
       'Q': 2.2366176229999999,
       'R': 0.080835481000000001,
       'S': 0.36310446600000001,
       'T': 0.13480267100000001,
       'V': 0.227827051,
       'W': 0.086028795000000005,
       'Y': 0.093214721},
 'F': {'G': 0.068139280999999996,
       'H': 0.090353066999999995,
       'I': 0.75090054100000003,
       'K': 0.097125533999999999,
       'L': 1.811101233,
       'M': 0.74842499699999998,
       'N': 0.13875929100000001,
       'P': 0.054271888999999997,
       'Q': 0.026306324999999998,
       'R': 0.0080439580000000004,
       'S': 0.49934990099999998,
       'T': 0.053947742999999999,
       'V': 0.46624344200000001,
       'W': 0.330781928,
       'Y': 3.2090833029999999},
 'G': {'H': 0.025252655999999998,
       'I': 0.013781055,
       'K': 0.13418717499999999,
       'L': 0.027264554,
       'M': 0.14533146599999999,
       'N': 0.54375075699999997,
       'P': 0.005914206,
       'Q': 0.072395535999999996,
       'R': 0.21996712400000001,
       'S': 1.746570145,
       'T': 0.02455829,
       'V': 0.41714895400000002,
       'W': 0.233963371,
       'Y': 0.046746340999999997},
 'H': {'I': 0.017140138999999999,
       'K': 0.13515366300000001,
       'L': 0.11163893699999999,
       'M': 0.032834314000000003,
       'N': 1.7386796440000001,
       'P': 0.51995437499999997,
       'Q': 4.5184508909999996,
       'R': 1.5222568649999999,
       'S': 0.29758608399999997,
       'T': 0.221010609,
       'V': 0.0035110079999999999,
       'W': 0.037480926999999997,
       'Y': 3.9079185509999999},
 'I': {'K': 0.064936611000000005,
       'L': 1.8979743680000001,
       'M': 2.9183532080000001,
       'N': 0.244934765,
       'P': 0.024850021,
       'Q': 0.0088756860000000007,
       'R': 0.012428576,
       'S': 0.096272864,
       'T': 2.4534581430000002,
       'V': 10.953425842,
       'W': 0.028656797000000001,
       'Y': 0.135319461},
 'K': {'L': 0.06132452,
       'M': 0.65931076,
       'N': 2.53039843,
       'P': 0.121234921,
       'Q': 1.8272181860000001,
       'R': 1.057185633,
       'S': 0.69508812799999997,
       'T': 0.39385170400000002,
       'V': 0.055461435000000003,
       'W': 0.073508962999999997,
       'Y': 0.281699174},
 'L': {'M': 3.4255537089999999,
       'N': 0.046318944000000001,
       'P': 0.27026078100000001,
       'Q': 0.25445246700000002,
       'R': 0.058180015000000002,
       'S': 0.31152513100000001,
       'T': 0.253366704,
       'V': 0.95827374300000001,
       'W': 0.25324301300000002,
       'Y': 0.123555332},
 'M': {'N': 0.39982772300000002,
       'P': 0.032714699,
       'Q': 0.237094366,
       'R': 0.013494034,
       'S': 0.45873409599999998,
       'T': 3.0352157260000001,
       'V': 2.5624848949999999,
       'W': 0.167575318,
       'Y': 0.31659903099999998},
 'N': {'P': 0.080313958000000005,
       'Q': 0.83279153299999997,
       'R': 0.14136427500000001,
       'S': 2.6343785139999998,
       'T': 0.96128509299999998,
       'V': 0.051741626999999998,
       'W': 0.049019408,
       'Y': 1.020785491},
 'P': {'Q': 0.84951243499999995,
       'R': 0.15500856599999999,
       'S': 1.231180819,
       'T': 0.73460491000000006,
       'V': 0.054078532999999998,
       'W': 0.029433866,
       'Y': 0.054012182999999998},
 'Q': {'R': 2.6731080889999999,
       'S': 0.38480028399999999,
       'T': 0.274195947,
       'V': 0.027669233000000001,
       'W': 0.12314062000000001,
       'Y': 0.319105788},
 'R': {'S': 0.19737918500000001,
       'T': 0.056079812999999999,
       'V': 0.041063684000000003,
       'W': 0.37081989199999998,
       'Y': 0.12751933200000001},
 'S': {'T': 3.1147429070000001,
       'V': 0.26710946499999999,
       'W': 0.16921202900000001,
       'Y': 0.37418428599999998},
 'T': {'V': 1.5140596740000001, 'W': 0.014378616, 'Y': 0.091031787000000003},
 'V': {'W': 0.093136256000000001, 'Y': 0.069964540000000006},
 'W': {'Y': 0.48104431600000003},
 'Y': {}}