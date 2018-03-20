LoadFunctionLibrary("../DNA.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");
LoadFunctionLibrary ("libv3/all-terms.bf");

/** @module models.DNA.JC69 */

/**
 * @name models.DNA.JC69.ModelDescription
 * @param {String} type
 */
lfunction models.DNA.JC69.ModelDescription(type) {

    return {
        utility.getGlobalValue("terms.alphabet"): utility.getGlobalValue("models.DNA.alphabet"),
        utility.getGlobalValue("terms.description"): "The JC69 equal-rates model of nucleotide substitution",
        utility.getGlobalValue("terms.model.canonical"): 1, // is of the r_ij \times \pi_j form
        utility.getGlobalValue("terms.model.reversible"): 1,
        utility.getGlobalValue("terms.model.efv_estimate_name"): utility.getGlobalValue("terms.frequencies.equal"),
        utility.getGlobalValue("terms.parameters"): {
            utility.getGlobalValue("terms.global"): {},
            utility.getGlobalValue("terms.local"): {},
            utility.getGlobalValue("terms.model.empirical"): 0
        },
        utility.getGlobalValue("terms.model.type"): type,
        utility.getGlobalValue("terms.model.get_branch_length"): "",
        utility.getGlobalValue("terms.model.set_branch_length"): "models.generic.SetBranchLength",
        utility.getGlobalValue("terms.model.constrain_branch_length"): "models.generic.ConstrainBranchLength",
        utility.getGlobalValue("terms.model.frequency_estimator"): "frequencies.equal",
        utility.getGlobalValue("terms.model.q_ij"): "models.DNA.JC69._generateRate",
        utility.getGlobalValue("terms.model.time"): "models.DNA.generic.Time",
        utility.getGlobalValue("terms.model.defineQ"): "models.DNA.JC69._defineQ",
        utility.getGlobalValue("terms.model.post_definition"): "models.generic.post.definition"
    };
}

/**
 * @name models.DNA.JC69._generateRate
 * @param {Number} fromChar
 * @param {Number} toChar
 * @param {String} namespace
 * @param {String} model_type
 * @returns models.DNA.JC69._generateRate.p
 */
function models.DNA.JC69._generateRate(fromChar, toChar, namespace, model_type, model) {
    models.DNA.JC69._generateRate.p = {};
    models.DNA.JC69._generateRate.p[model_type] = {};

    models.DNA.JC69.parameter_name = "theta";

    if (model_type == terms.global) {
        models.DNA.JC69.parameter_name = parameters.ApplyNameSpace(models.DNA.JC69.parameter_name, namespace);
    }

    (models.DNA.JC69._generateRate.p[model_type])[terms.nucleotideRate(fromChar, toChar)] = models.DNA.JC69.parameter_name;
    models.DNA.JC69._generateRate.p[terms.model.rate_entry] = models.DNA.JC69.parameter_name;

    return models.DNA.JC69._generateRate.p;
}

/**
 * @name models.DNA.JC69._defineQ
 * @param {Dictionary} jc69
 * @param {String} namespace
 */
function models.DNA.JC69._defineQ(jc69, namespace) {

    models.DNA.generic.DefineQMatrix(jc69, namespace);
    return jc69;
}
