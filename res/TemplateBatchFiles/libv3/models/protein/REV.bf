LoadFunctionLibrary("../protein.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");

/** @module models.protein.REV */

/**
 * @name models.protein.REV.ModelDescription
 * @param {String} type
 */

lfunction models.protein.REV.ModelDescription(type) {
    return {
        utility.getGlobalValue("terms.alphabet"): utility.getGlobalValue("models.protein.alphabet"),
        utility.getGlobalValue("terms.description"): "General time reversible model for protein sequences",
        utility.getGlobalValue("terms.model.canonical"): 1, // is of the r_ij \times \pi_j form
        utility.getGlobalValue("terms.model.reversible"): 1,
        utility.getGlobalValue("terms.model.efv_estimate_name"): utility.getGlobalValue("terms.frequencies._20x1"),
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
        utility.getGlobalValue("terms.model.q_ij"): "models.protein.REV._GenerateRate",
        utility.getGlobalValue("terms.model.time"): "models.protein.generic.Time",
        utility.getGlobalValue("terms.model.defineQ"): "models.protein.REV._DefineQ",
        utility.getGlobalValue("terms.model.post_definition"): "models.generic.post.definition"
    };
}

lfunction models.protein.REV._GenerateRate(fromChar, toChar, namespace, model_type) {
    models.protein.REV._generateRate.p = {};
    models.protein.REV._generateRate.p[model_type] = {};

    if (fromChar > toChar) {
        models.protein.REV.parameter_name = "theta_" + toChar + fromChar;
    } else {
        models.protein.REV.parameter_name = "theta_" + fromChar + toChar;
    }

    if (model_type == utility.getGlobalValue("terms.global")) {
        models.protein.REV.parameter_name = parameters.ApplyNameSpace(models.protein.REV.parameter_name, namespace);
    }

    (models.protein.REV._generateRate.p[model_type])[terms.aminoacidRate(fromChar, toChar)] = models.protein.REV.parameter_name;
    models.protein.REV._generateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = models.protein.REV.parameter_name;

    return models.protein.REV._generateRate.p;

}

/**
 * @name models.protein.REV._DefineQ
 * @param {Dictionary} model definition
 * @param {String} namespace
 */
lfunction models.protein.REV._DefineQ(model_dict, namespace) {
    models.protein.generic.DefineQMatrix (model_dict, namespace);
    if (model_dict[utility.getGlobalValue("terms.model.type")] == utility.getGlobalValue("terms.global")) {
        parameters.SetConstraint(((model_dict[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.global")])[terms.aminoacidRate("I", "L")], "1", "");
    }
    if (model_dict[utility.getGlobalValue("terms.model.type")] == utility.getGlobalValue("terms.local")) {
        parameters.SetConstraint(((model_dict[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.local")])[terms.aminoacidRate("I", "L")], "1", "");
    }


    return model_dict;
}
