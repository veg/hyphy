LoadFunctionLibrary("../protein.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");
LoadFunctionLibrary("../../all-terms.bf");


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

/*
Load model matrices
*/
LoadFunctionLibrary("matrices/JC69.ibf");
LoadFunctionLibrary("matrices/JTT.ibf");
LoadFunctionLibrary("matrices/LG.ibf");
LoadFunctionLibrary("matrices/WAG.ibf");
LoadFunctionLibrary("matrices/mt.ibf"); // all three mtMet/mtVer/mtInv models
LoadFunctionLibrary("matrices/gcpREV.ibf");
LoadFunctionLibrary("matrices/HIV.ibf");

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


models.protein.empirical.default_generators = {"LG": "models.protein.LG.ModelDescription",
                                               "WAG": "models.protein.WAG.ModelDescription",
                                               "JTT": "models.protein.JTT.ModelDescription",
                                               "JC69": "models.protein.JC69.ModelDescription",
                                               "mtInv": "models.protein.mtInv.ModelDescription",
                                               "mtMet": "models.protein.mtMet.ModelDescription",
                                               "mtVer": "models.protein.mtVer.ModelDescription",
                                               "gcpREV": "models.protein.gcpREV.ModelDescription",
                                               "HIVBm": "models.protein.HIVBm.ModelDescription",
                                               "HIVWm": "models.protein.HIVWm.ModelDescription"};
                                           
models.protein.empirical.plusF_generators = {"LG": "models.protein.LGF.ModelDescription",
                                             "WAG": "models.protein.WAGF.ModelDescription",
                                             "JTT": "models.protein.JTTF.ModelDescription",
                                             "JC69": "models.protein.JC69F.ModelDescription",
                                             "mtMet": "models.protein.mtMetF.ModelDescription",
                                             "mtVer": "models.protein.mtVerF.ModelDescription",
                                             "gcpREV": "models.protein.gcpREVF.ModelDescription",
                                             "HIVBm": "models.protein.HIVBmF.ModelDescription",
                                             "HIVWm": "models.protein.HIVWmF.ModelDescription"};
                                             
models.protein.empirical.mleF_generators = {"LG": "models.protein.LGML.ModelDescription",
                                             "WAG": "models.protein.WAGML.ModelDescription",
                                             "JTT": "models.protein.JTTML.ModelDescription",
                                             "JC69": "models.protein.JC69ML.ModelDescription",
                                             "mtMet": "models.protein.mtMetML.ModelDescription",
                                             "mtVer": "models.protein.mtVerML.ModelDescription",
                                             "gcpREV": "models.protein.gcpREVML.ModelDescription",
                                             "HIVBm": "models.protein.HIVBmML.ModelDescription",
                                             "HIVWm": "models.protein.HIVWmML.ModelDescription"};










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


/**************************************** mtVer functions *************************************/


/**
 * @name models.protein.mtVer.ModelDescription
 * @description Create the baseline schema (dictionary) for the mtVer model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
 function models.protein.mtVer.ModelDescription(type) {
    models.protein.mtVer.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.mtVer.ModelDescription.model_definition [terms.model.empirical_rates] = models.protein.mtVer.Rij;
    models.protein.mtVer.ModelDescription.model_definition [terms.model.frequency_estimator] = "models.protein.mtVer.frequencies";
    return models.protein.mtVer.ModelDescription.model_definition;
}

/**
 * @name models.protein.mtVerF.ModelDescription
 * @description Create the baseline schema (dictionary) for the mtVer+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.mtVerF.ModelDescription(type) {
    models.protein.mtVerF.ModelDescription.model_definition = models.protein.mtVer.ModelDescription(type);
    models.protein.mtVerF.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.empirical.protein";
    models.protein.mtVerF.ModelDescription.model_definition [terms.model.efv_estimate_name] = utility.getGlobalValue("terms.frequencies._20x1");
    return models.protein.mtVerF.ModelDescription.model_definition;
}


/**
 * @name models.protein.mtVerML.ModelDescription
 * @description Create the baseline schema (dictionary) for the mtVer+ML model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.mtVerML.ModelDescription(type) {
    models.protein.mtVerML.ModelDescription.model_definition = models.protein.mtVer.ModelDescription(type);
    models.protein.mtVerML.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.ML.protein";
    models.protein.mtVerML.ModelDescription.model_definition [terms.model.efv_estimate_name]   =  utility.getGlobalValue("terms.frequencies.MLE");
    return models.protein.mtVerML.ModelDescription.model_definition;
}



/**************************************** gcpREV functions *************************************/


/**
 * @name models.protein.gcpREV.ModelDescription
 * @description Create the baseline schema (dictionary) for the gcpREV model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
 function models.protein.gcpREV.ModelDescription(type) {
    models.protein.gcpREV.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.gcpREV.ModelDescription.model_definition [terms.model.empirical_rates] = models.protein.gcpREV.Rij;
    models.protein.gcpREV.ModelDescription.model_definition [terms.model.frequency_estimator] = "models.protein.gcpREV.frequencies";
    return models.protein.gcpREV.ModelDescription.model_definition;
}

/**
 * @name models.protein.gcpREVF.ModelDescription
 * @description Create the baseline schema (dictionary) for the gcpREV+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.gcpREVF.ModelDescription(type) {
    models.protein.gcpREVF.ModelDescription.model_definition = models.protein.gcpREV.ModelDescription(type);
    models.protein.gcpREVF.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.empirical.protein";
    models.protein.gcpREVF.ModelDescription.model_definition [terms.model.efv_estimate_name] = utility.getGlobalValue("terms.frequencies._20x1");
    return models.protein.gcpREVF.ModelDescription.model_definition;
}


/**
 * @name models.protein.gcpREVML.ModelDescription
 * @description Create the baseline schema (dictionary) for the gcpREV+ML model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.gcpREVML.ModelDescription(type) {
    models.protein.gcpREVML.ModelDescription.model_definition = models.protein.gcpREV.ModelDescription(type);
    models.protein.gcpREVML.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.ML.protein";
    models.protein.gcpREVML.ModelDescription.model_definition [terms.model.efv_estimate_name]   =  utility.getGlobalValue("terms.frequencies.MLE");
    return models.protein.gcpREVML.ModelDescription.model_definition;
}





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








 
lfunction models.protein.gcpREV.frequencies (model, namespace, datafilter) {
    model[utility.getGlobalValue("terms.efv_estimate")] =
      {{ 0.07951}
        {0.009051}
        {0.03322}
        {0.049675}
        {0.047731}
        {0.080233}
        {0.02188}
        {0.080496}
        {0.049324}
        {0.107512}
        {0.020776}
        {0.040459}
        {0.039916}
        {0.037505}
        {0.056001}
        {0.07382}
        {0.053615}
        {0.071781}
        {0.016705}
        {0.03079}};
       
    model[utility.getGlobalValue("terms.model.efv_estimate_name")] = utility.getGlobalValue("terms.frequencies.predefined");
    (model[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.model.empirical")] = 0;
    return model;
}

models.protein.gcpREV.Rij = {
'A': {'C': 699.39999999999998,
       'D': 87.700000000000003,
       'E': 288.30000000000001,
       'F': 59.299999999999997,
       'G': 677.89999999999998,
       'H': 32.0,
       'I': 39.5,
       'K': 78.299999999999997,
       'L': 33.799999999999997,
       'M': 197.09999999999999,
       'N': 59.100000000000001,
       'P': 505.80000000000001,
       'Q': 199.5,
       'R': 12.199999999999999,
       'S': 2443.9000000000001,
       'T': 1646.2,
       'V': 1143.8,
       'W': 30.0,
       'Y': 27.899999999999999},
 'C': {'D': 125.3,
       'E': 102.40000000000001,
       'F': 1540.0999999999999,
       'G': 433.60000000000002,
       'H': 574.79999999999995,
       'I': 96.0,
       'K': 89.200000000000003,
       'L': 446.5,
       'M': 288.5,
       'N': 133.69999999999999,
       'P': 124.90000000000001,
       'Q': 195.09999999999999,
       'R': 1867.9000000000001,
       'S': 3029.8000000000002,
       'T': 516.79999999999995,
       'V': 722.20000000000005,
       'W': 562.60000000000002,
       'Y': 1571.2},
 'D': {'E': 3638.6999999999998,
       'F': 15.9,
       'G': 464.80000000000001,
       'H': 527.70000000000005,
       'I': 16.800000000000001,
       'K': 78.299999999999997,
       'L': 8.6999999999999993,
       'M': 33.100000000000001,
       'N': 5083.8000000000002,
       'P': 43.600000000000001,
       'Q': 387.60000000000002,
       'R': 40.600000000000001,
       'S': 142.09999999999999,
       'T': 66.099999999999994,
       'V': 18.5,
       'W': 35.0,
       'Y': 272.10000000000002},
 'E': {'F': 30.699999999999999,
       'G': 393.80000000000001,
       'H': 144.59999999999999,
       'I': 23.399999999999999,
       'K': 1870.9000000000001,
       'L': 26.800000000000001,
       'M': 33.299999999999997,
       'N': 230.90000000000001,
       'P': 39.399999999999999,
       'Q': 2965.3000000000002,
       'R': 51.899999999999999,
       'S': 156.69999999999999,
       'T': 230.0,
       'V': 103.8,
       'W': 57.200000000000003,
       'Y': 79.200000000000003},
 'F': {'G': 10.1,
       'H': 87.200000000000003,
       'I': 354.39999999999998,
       'K': 9.0,
       'L': 1585.8,
       'M': 183.40000000000001,
       'N': 13.4,
       'P': 82.099999999999994,
       'Q': 10.699999999999999,
       'R': 26.100000000000001,
       'S': 790.5,
       'T': 29.800000000000001,
       'V': 180.59999999999999,
       'W': 558.29999999999995,
       'Y': 2360.1999999999998},
 'G': {'H': 35.600000000000001,
       'I': 11.800000000000001,
       'K': 193.40000000000001,
       'L': 4.2999999999999998,
       'M': 11.800000000000001,
       'N': 446.5,
       'P': 11.800000000000001,
       'Q': 40.799999999999997,
       'R': 391.19999999999999,
       'S': 709.0,
       'T': 65.5,
       'V': 72.599999999999994,
       'W': 92.599999999999994,
       'Y': 20.600000000000001},
 'H': {'I': 38.100000000000001,
       'K': 96.900000000000006,
       'L': 46.100000000000001,
       'M': 54.5,
       'N': 2156.8000000000002,
       'P': 185.09999999999999,
       'Q': 1889.5999999999999,
       'R': 1228.7,
       'S': 277.0,
       'T': 77.099999999999994,
       'V': 22.800000000000001,
       'W': 58.5,
       'Y': 4117.5},
 'I': {'K': 125.7,
       'L': 1657.9000000000001,
       'M': 2493.4000000000001,
       'N': 146.0,
       'P': 63.100000000000001,
       'Q': 48.5,
       'R': 105.40000000000001,
       'S': 75.599999999999994,
       'T': 1227.2,
       'V': 6006.5,
       'W': 21.100000000000001,
       'Y': 35.799999999999997},
 'K': {'L': 37.600000000000001,
       'M': 220.09999999999999,
       'N': 1930.8,
       'P': 49.799999999999997,
       'Q': 2571.1999999999998,
       'R': 4666.3999999999996,
       'S': 209.30000000000001,
       'T': 648.5,
       'V': 66.200000000000003,
       'W': 15.9,
       'Y': 84.200000000000003},
 'L': {'M': 1389.5999999999999,
       'N': 11.6,
       'P': 347.69999999999999,
       'Q': 235.40000000000001,
       'R': 124.90000000000001,
       'S': 584.89999999999998,
       'T': 62.899999999999999,
       'V': 661.70000000000005,
       'W': 218.19999999999999,
       'Y': 82.0},
 'M': {'N': 41.399999999999999,
       'P': 19.199999999999999,
       'Q': 341.80000000000001,
       'R': 112.09999999999999,
       'S': 32.899999999999999,
       'T': 1143.5,
       'V': 465.69999999999999,
       'W': 70.5,
       'Y': 41.100000000000001},
 'N': {'P': 18.699999999999999,
       'Q': 827.79999999999995,
       'R': 313.0,
       'S': 2505.4000000000001,
       'T': 1373.3,
       'V': 27.600000000000001,
       'W': 19.899999999999999,
       'Y': 316.19999999999999},
 'P': {'Q': 223.09999999999999,
       'R': 89.400000000000006,
       'S': 1154.5,
       'T': 322.80000000000001,
       'V': 75.5,
       'W': 18.5,
       'Y': 15.800000000000001},
 'Q': {'R': 2321.3000000000002,
       'S': 246.90000000000001,
       'T': 128.5,
       'V': 50.200000000000003,
       'W': 24.600000000000001,
       'Y': 179.09999999999999},
 'R': {'S': 269.60000000000002,
       'T': 173.40000000000001,
       'V': 62.0,
       'W': 257.80000000000001,
       'Y': 137.0},
 'S': {'T': 2042.5,
       'V': 61.200000000000003,
       'W': 102.59999999999999,
       'Y': 501.89999999999998},
 'T': {'V': 583.70000000000005, 'W': 22.0, 'Y': 64.599999999999994},
 'V': {'W': 10.800000000000001, 'Y': 37.299999999999997},
 'W': {'Y': 296.30000000000001},
 'Y': {}};