LoadFunctionLibrary("../DNA.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");

/** @module models.DNA.HKY85 */

/**
 * @name models.DNA.HKY85.ModelDescription
 * @param {String} type
 */
function models.DNA.HKY85.ModelDescription(type) {

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
        "set-branch-length": "models.generic.SetBranchLength",
        "constrain-branch-length": "models.generic.constrain_branch_length",
        "frequency-estimator": "frequencies.empirical.nucleotide",
        "q_ij": "models.DNA.HKY85._GenerateRate",
        "time": "models.DNA.generic.Time",
        "defineQ": "models.DNA.HKY85._DefineQ",
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
 * @name models.DNA.HKY85._GenerateRate
 * @param {Number} fromChar
 * @param {Number} toChar
 * @param {String} namespace
 * @param {String} model_type
 * @return models.DNA.HKY85._GenerateRate.p 
 */
function models.DNA.HKY85._GenerateRate(fromChar, toChar, namespace, model_type) {
    models.DNA.HKY85._GenerateRate.p = {};
    models.DNA.HKY85._GenerateRate.p[model_type] = {};


    if (model_type == terms.global) {
        models.DNA.HKY85.parameter_name = parameters.ApplyNameSpace(models.DNA.HKY85.parameter_name, namespace);
        if (models.DNA.HKY85.is_transition(fromChar, toChar)) {
            models.DNA.HKY85.parameter_name = parameters.ApplyNameSpace("kappa", namespace);
            (models.DNA.HKY85._GenerateRate.p[model_type])[terms.transition_transversion_ratio] = models.DNA.HKY85.parameter_name;
        } else {
            models.DNA.HKY85.parameter_name = "1";
        }
    } else {
        if (models.DNA.HKY85.is_transition(fromChar, toChar)) {
            models.DNA.HKY85.parameter_name = "transition";
            (models.DNA.HKY85._GenerateRate.p[model_type])[terms.transition] = models.DNA.HKY85.parameter_name;
        } else {
            models.DNA.HKY85.parameter_name = "transversion";
            (models.DNA.HKY85._GenerateRate.p[model_type])[terms.transversion] = models.DNA.HKY85.parameter_name;
        }
    }
    models.DNA.HKY85._GenerateRate.p[terms.rate_entry] = models.DNA.HKY85.parameter_name;
    return models.DNA.HKY85._GenerateRate.p;
}

/**
 * @name models.DNA.HKY85._DefineQ
 * @param {Dictionary} hky85
 * @param {String} namespace
 */
function models.DNA.HKY85._DefineQ(hky85, namespace) {
    models.DNA.generic.DefineQMatrix(hky85, namespace);
    return hky85;
}
