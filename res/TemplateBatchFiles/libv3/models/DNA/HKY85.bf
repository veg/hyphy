LoadFunctionLibrary("../DNA.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");

/**
 * @name models.DNA.HKY85.modelDescription
 * @param {String} type
 */
function models.DNA.HKY85.modelDescription(type) {

    return {
        "alphabet": models.DNA.alphabet,
        "description": "Hasegawa Kishino Yano 85 (HKY85) model of nucleotide substitution",
        "canonical": 1, // is of the r_ij \times \pi_j form
        "reversible": 1,
        terms.efv_estimate_name: terms.freqs.4x1,
        "parameters": {
            "global": {},
            "local": {},
            "empirical": 3
        },
        "type": type,
        "get-branch-length": "",
        "set-branch-length": "models.generic.set_branch_length",
        "constrain-branch-length": "models.generic.constrain_branch_length",
        "frequency-estimator": "frequencies.empirical.nucleotide",
        "q_ij": "models.DNA.HKY85.generateRate",
        "time": "models.DNA.generic.time",
        "defineQ": "models.DNA.HKY85.defineQ",
        "post-definition": "models.generic.post.definition"
    };
}

/**
 * @name models.DNA.HKY85.is_transition
 * @param {Number} fromChar
 * @param {Number} toChar
 */
function models.DNA.HKY85.is_transition(fromChar, toChar) {
    return fromChar == "A" && toChar == "G" || fromChar == "G" && toChar == "A" || fromChar == "C" && toChar == "T" || fromChar == "T" && toChar == "C";
}

/**
 * @name models.DNA.HKY85.generateRate
 * @param {Number} fromChar
 * @param {Number} toChar
 * @param {String} namespace
 * @param {String} model_type
 * @return models.DNA.HKY85.generateRate.p 
 */
function models.DNA.HKY85.generateRate(fromChar, toChar, namespace, model_type) {
    models.DNA.HKY85.generateRate.p = {};
    models.DNA.HKY85.generateRate.p[model_type] = {};


    if (model_type == terms.global) {
        models.DNA.HKY85.parameter_name = parameters.applyNameSpace(models.DNA.HKY85.parameter_name, namespace);
        if (models.DNA.HKY85.is_transition(fromChar, toChar)) {
            models.DNA.HKY85.parameter_name = parameters.applyNameSpace("kappa", namespace);
            (models.DNA.HKY85.generateRate.p[model_type])[terms.transition_transversion_ratio] = models.DNA.HKY85.parameter_name;
        } else {
            models.DNA.HKY85.parameter_name = "1";
        }
    } else {
        if (models.DNA.HKY85.is_transition(fromChar, toChar)) {
            models.DNA.HKY85.parameter_name = "transition";
            (models.DNA.HKY85.generateRate.p[model_type])[terms.transition] = models.DNA.HKY85.parameter_name;
        } else {
            models.DNA.HKY85.parameter_name = "transversion";
            (models.DNA.HKY85.generateRate.p[model_type])[terms.transversion] = models.DNA.HKY85.parameter_name;
        }
    }
    models.DNA.HKY85.generateRate.p[terms.rate_entry] = models.DNA.HKY85.parameter_name;
    return models.DNA.HKY85.generateRate.p;
}

/**
 * @name models.DNA.HKY85.defineQ
 * @param {Dictionary} hky85
 * @param {String} namespace
 */
function models.DNA.HKY85.defineQ(hky85, namespace) {
    models.DNA.generic.defineQMatrix(hky85, namespace);
    return hky85;
}
