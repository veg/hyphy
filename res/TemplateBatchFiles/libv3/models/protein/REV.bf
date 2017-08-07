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
        "alphabet": models.protein.alphabet,
        "description": "General time reversible model for protein sequences",
        "canonical": 1, // is of the r_ij \times \pi_j form
        "reversible": 1,
        terms.efv_estimate_name: terms.freqs.20x1,
        "parameters": {
            "global": {},
            "local": {},
            "empirical": 0
        },
        "type": type,
        "get-branch-length": "",
        "set-branch-length": "models.generic.SetBranchLength",
        "constrain-branch-length": "models.generic.constrain_branch_length",
        "frequency-estimator": "frequencies.empirical.protein",
        "q_ij": "models.protein.REV._GenerateRate",
        "time": "models.protein.generic.Time",
        "defineQ": "models.protein.REV._DefineQ",
        "post-definition": "models.generic.post.definition"
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
    models.protein.REV._generateRate.p[terms.rate_entry] = models.protein.REV.parameter_name;

    return models.protein.REV._generateRate.p;

}

/**
 * @name models.protein.REV._DefineQ
 * @param {Dictionary} model definition
 * @param {String} namespace
 */
function models.protein.REV._DefineQ(model_dict, namespace) {
    models.protein.generic.DefineQMatrix (model_dict, namespace);
    if (model_dict["type"] == terms.global) {
        parameters.SetConstraint(((model_dict["parameters"])[terms.global])[terms.aminoacidRate("I", "L")], "1", "");
    }
    if (model_dict["type"] == terms.local) {
        parameters.SetConstraint(((model_dict["parameters"])[terms.local])[terms.aminoacidRate("I", "L")], "1", "");
    }


    return model_dict;
}
