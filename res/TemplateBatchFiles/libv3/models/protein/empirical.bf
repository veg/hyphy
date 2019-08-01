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
        utility.getGlobalValue("terms.model.constrain_branch_length"): "models.generic.ConstrainBranchLength",
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
lfunction models.protein.empirical._GenerateRate(fromChar, toChar, namespace, model_type, model) {
    models.protein.empirical._GenerateRate.p = {};
    models.protein.empirical._GenerateRate.p  [model_type]       = {};

    if (fromChar < toChar) {
        models.protein.empirical._GenerateRate.p  [utility.getGlobalValue("terms.model.rate_entry")] = "" + ((model[utility.getGlobalValue ("terms.model.empirical_rates")])[fromChar])[toChar];
    } else {
        models.protein.empirical._GenerateRate.p  [utility.getGlobalValue("terms.model.rate_entry")] = "" + ((model[utility.getGlobalValue ("terms.model.empirical_rates")])[toChar])[fromChar];
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

    if (utility.Has (model_dict, utility.getGlobalValue("terms.model.data"), "String")) {
        frequencies._aux.empirical.singlechar(model_dict, namespace, model_dict[utility.getGlobalValue("terms.model.data")]);
        models.protein.empirical._NormalizeEmpiricalRates(model_dict, namespace);
    }
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



