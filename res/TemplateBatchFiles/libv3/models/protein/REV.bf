LoadFunctionLibrary("../protein.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");

/** @module models.protein.REV */

/**
 * @name models.protein.REV.ModelDescription
 * @param {String} type
 */

function models.protein.REV.ModelDescription(type) {
    return {
        terms.alphabet: models.protein.alphabet,
        terms.description: "General time reversible model for protein sequences",
        terms.model.canonical: 1, // is of the r_ij \times \pi_j form
        terms.model.reversible: 1,
        terms.model.efv_estimate_name: terms.frequencies._20x1,
        terms.parameters: {
            terms.global: {},
            terms.local: {},
            terms.model.empirical: 0
        },
        terms.model.type: type,
        terms.model.get_branch_length: "",
        terms.model.set_branch_length: "models.generic.SetBranchLength",
        terms.model.constrain_branch_length: "models.generic.constrain_branch_length",
        terms.model.frequency_estimator: "frequencies.empirical.protein",
        terms.model.q_ij: "models.protein.REV._GenerateRate",
        terms.model.time: "models.protein.generic.Time",
        terms.model.defineQ: "models.protein.REV._DefineQ",
        terms.model.post_definition: "models.generic.post.definition"
    };
}

function models.protein.REV._GenerateRate(fromChar, toChar, namespace, model_type) {
    models.protein.REV._generateRate.p = {};
    models.protein.REV._generateRate.p[model_type] = {};

    if (fromChar > toChar) {
        models.protein.REV.parameter_name = "theta_" + toChar + fromChar;
    } else {
        models.protein.REV.parameter_name = "theta_" + fromChar + toChar;
    }

    if (model_type == terms.global) {
        models.protein.REV.parameter_name = parameters.ApplyNameSpace(models.protein.REV.parameter_name, namespace);
    }

    (models.protein.REV._generateRate.p[model_type])[terms.aminoacidRate(fromChar, toChar)] = models.protein.REV.parameter_name;
    models.protein.REV._generateRate.p[terms.model.rate_entry] = models.protein.REV.parameter_name;

    return models.protein.REV._generateRate.p;

}

/**
 * @name models.protein.REV._DefineQ
 * @param {Dictionary} model definition
 * @param {String} namespace
 */
function models.protein.REV._DefineQ(model_dict, namespace) {
    models.protein.generic.DefineQMatrix (model_dict, namespace);
    if (model_dict[terms.model.type] == terms.global) {
        parameters.SetConstraint(((model_dict[terms.parameters])[terms.global])[terms.aminoacidRate("I", "L")], "1", "");
    }
    if (model_dict[terms.model.type] == terms.local) {
        parameters.SetConstraint(((model_dict[terms.parameters])[terms.local])[terms.aminoacidRate("I", "L")], "1", "");
    }


    return model_dict;
}
