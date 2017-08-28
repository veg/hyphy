LoadFunctionLibrary("../DNA.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");

/** @module models.DNA.GTR */

/**
 * @name models.DNA.GTR.ModelDescription
 * @param {String} type
 */
function models.DNA.GTR.ModelDescription(type) {

    return {
        "alphabet": models.DNA.alphabet,
        "description": "The general time reversible (GTR) model of nucleotide substitution",
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
        "q_ij": "models.DNA.GTR._generateRate",
        "time": "models.DNA.generic.Time",
        "defineQ": "models.DNA.GTR._defineQ",
        "post-definition": "models.generic.post.definition"
    };
}

/**
 * @name models.DNA.GTR._generateRate
 * @param {Number} fromChar
 * @param {Number} toChar
 * @param {String} namespace
 * @param {String} model_type
 * @returns models.DNA.GTR._generateRate.p
 */
function models.DNA.GTR._generateRate(fromChar, toChar, namespace, model_type) {
    models.DNA.GTR._generateRate.p = {};
    models.DNA.GTR._generateRate.p[model_type] = {};

    if (fromChar > toChar) {
        models.DNA.GTR.parameter_name = "theta_" + toChar + fromChar;
    } else {
        models.DNA.GTR.parameter_name = "theta_" + fromChar + toChar;
    }

    if (model_type == terms.global) {
        models.DNA.GTR.parameter_name = parameters.ApplyNameSpace(models.DNA.GTR.parameter_name, namespace);
    }

    (models.DNA.GTR._generateRate.p[model_type])[terms.nucleotideRate(fromChar, toChar)] = models.DNA.GTR.parameter_name;
    models.DNA.GTR._generateRate.p[terms.rate_entry] = models.DNA.GTR.parameter_name;

    return models.DNA.GTR._generateRate.p;
}

/**
 * @name models.DNA.GTR._defineQ
 * @param {Dictionary} gtr
 * @param {String} namespace
 */
function models.DNA.GTR._defineQ(gtr, namespace) {

    models.DNA.generic.DefineQMatrix(gtr, namespace);
    if (gtr["type"] == terms.global) {
        parameters.SetConstraint(((gtr["parameters"])[terms.global])[terms.nucleotideRate("A", "G")], "1", "");
    }
    if (gtr["type"] == terms.local) {
        parameters.SetConstraint(((gtr["parameters"])[terms.local])[terms.nucleotideRate("A", "G")], "1", "");
    }
    return gtr;
}
