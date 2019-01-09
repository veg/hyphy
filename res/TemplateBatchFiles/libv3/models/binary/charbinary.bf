LoadFunctionLibrary("../binary.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");

/** @module models.protein.REV */

/**
 * @name models.protein.REV.ModelDescription
 * @param {String} type
 */

lfunction models.binary.ModelDescription(type) {
    return {
        utility.getGlobalValue("terms.alphabet"): utility.getGlobalValue("models.binary.alphabet"),
        utility.getGlobalValue("terms.description"): "General time reversible model for binary characters",
        utility.getGlobalValue("terms.model.canonical"): 1, // is of the r_ij \times \pi_j form
        utility.getGlobalValue("terms.model.reversible"): 1,
        utility.getGlobalValue("terms.model.efv_estimate_name"): utility.getGlobalValue("terms.frequencies.binary"),
        utility.getGlobalValue("terms.parameters"): {
            utility.getGlobalValue("terms.global"): {},
            utility.getGlobalValue("terms.local"): {},
            utility.getGlobalValue("terms.model.empirical"): 1
        },
        utility.getGlobalValue("terms.model.type"): type,
        utility.getGlobalValue("terms.model.get_branch_length"): "",
        utility.getGlobalValue("terms.model.set_branch_length"): "models.generic.SetBranchLength",
        utility.getGlobalValue("terms.model.constrain_branch_length"): "models.generic.constrain_branch_length",
        utility.getGlobalValue("terms.model.frequency_estimator"): "frequencies.empirical.binary",
        utility.getGlobalValue("terms.model.q_ij"): "models.binary._GenerateRate",
        utility.getGlobalValue("terms.model.time"): "models.binary.generic.Time",
        utility.getGlobalValue("terms.model.defineQ"): "models.binary._DefineQ",
        utility.getGlobalValue("terms.model.post_definition"): "models.generic.post.definition"
    };
}

lfunction models.binary._GenerateRate(fromChar, toChar, namespace, model_type, model) {
    models.binary._generateRate.p = {}; 
    models.binary._generateRate.p[model_type] = {};

    if (fromChar > toChar) {
        models.binary.parameter_name = "mu_" + toChar + fromChar;
    } else {
        models.binary.parameter_name = "mu_" + fromChar + toChar;
    }

    if (model_type == utility.getGlobalValue("terms.global")) {
        models.binary.parameter_name = parameters.ApplyNameSpace(models.binary.parameter_name, namespace);
    }

    (models.binary._generateRate.p[model_type])[terms.binaryRate(fromChar, toChar)] = models.binary.parameter_name;

    models.binary._generateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = models.binary.parameter_name;

    return models.binary._generateRate.p;

}

/**
 * @name models.protein.REV._DefineQ
 * @param {Dictionary} model definition
 * @param {String} namespace
 */
lfunction models.binary._DefineQ(model_dict, namespace) {
    models.binary.generic.DefineQMatrix (model_dict, namespace);
/*
    if (model_dict[utility.getGlobalValue("terms.model.type")] == utility.getGlobalValue("terms.global")) {
        parameters.SetConstraint(  ((model_dict[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.global")])[terms.binaryRate("0", "1")], "1", "" );
    }
    if (model_dict[utility.getGlobalValue("terms.model.type")] == utility.getGlobalValue("terms.local")) {
        parameters.SetConstraint(((model_dict[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.local")])[terms.binaryRate("0", "1")], "1", "");
    }
*/

    return model_dict;
}
