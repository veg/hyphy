LoadFunctionLibrary("../DNA.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");

/** @module models.DNA.HKY85 */

/**
 * @name models.DNA.HKY85.ModelDescription
 * @param {String} type
 */
lfunction models.DNA.HKY85.ModelDescription(type) {

    return {
        utility.getGlobalValue("terms.alphabet"): utility.getGlobalValue("models.DNA.alphabet"),
        utility.getGlobalValue("terms.description"): "Hasegawa Kishino Yano 85 (HKY85) model of nucleotide substitution",
        utility.getGlobalValue("terms.model.canonical"): 1, // is of the r_ij \times \pi_j form
        utility.getGlobalValue("terms.model.reversible"): 1,
        utility.getGlobalValue("terms.model.efv_estimate_name"): utility.getGlobalValue("terms.frequencies._4x1"),
        utility.getGlobalValue("terms.parameters"): {
            utility.getGlobalValue("terms.global"): {},
            utility.getGlobalValue("terms.local"): {},
            utility.getGlobalValue("terms.model.empirical"): 3
        },
        utility.getGlobalValue("terms.model.type"): type,
        utility.getGlobalValue("terms.model.get_branch_length"): "",
        utility.getGlobalValue("terms.model.set_branch_length"): "models.generic.SetBranchLength",
        utility.getGlobalValue("terms.model.constrain_branch_length"): "models.generic.ConstrainBranchLength",
        utility.getGlobalValue("terms.model.frequency_estimator"): "frequencies.empirical.nucleotide",
        utility.getGlobalValue("terms.model.q_ij"): "models.DNA.HKY85._GenerateRate",
        utility.getGlobalValue("terms.model.time"): "models.DNA.generic.Time",
        utility.getGlobalValue("terms.model.defineQ"): "models.DNA.HKY85._DefineQ",
        utility.getGlobalValue("terms.model.post_definition"): "models.generic.post.definition"
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
function models.DNA.HKY85._GenerateRate(fromChar, toChar, namespace, model_type, model) {
    models.DNA.HKY85._GenerateRate.p = {};
    models.DNA.HKY85._GenerateRate.p[model_type] = {};


    if (model_type == terms.global) {
        if (models.DNA.HKY85.is_transition(fromChar, toChar)) {
            models.DNA.HKY85.parameter_name = parameters.ApplyNameSpace("kappa", namespace);
            (models.DNA.HKY85._GenerateRate.p[model_type])[terms.parameters.transition_transversion_ratio] = models.DNA.HKY85.parameter_name;
        } else {
            models.DNA.HKY85.parameter_name = "1";
        }
    } else {
        if (models.DNA.HKY85.is_transition(fromChar, toChar)) {
            models.DNA.HKY85.parameter_name = "transition";
            (models.DNA.HKY85._GenerateRate.p[model_type])[terms.parameters.transition] = models.DNA.HKY85.parameter_name;
        } else {
            models.DNA.HKY85.parameter_name = "transversion";
            (models.DNA.HKY85._GenerateRate.p[model_type])[terms.parameters.transversion] = models.DNA.HKY85.parameter_name;
        }
    }
    models.DNA.HKY85._GenerateRate.p[terms.model.rate_entry] = models.DNA.HKY85.parameter_name;
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
