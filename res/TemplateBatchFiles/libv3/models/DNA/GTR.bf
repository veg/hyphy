LoadFunctionLibrary("../DNA.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");
LoadFunctionLibrary ("libv3/all-terms.bf");

/** @module models.DNA.GTR */

/**
 * @name models.DNA.GTR.ModelDescription
 * @param {String} type
 */
function models.DNA.GTR.ModelDescription(type) {

    return {
        terms.alphabet: models.DNA.alphabet,
        terms.description: "The general time reversible (GTR) model of nucleotide substitution",
        terms.model.canonical: 1, // is of the r_ij \times \pi_j form
        terms.model.reversible: 1,
        terms.model.efv_estimate_name: terms.freqs.4x1,
        terms.parameters: {
            terms.global: {},
            terms.local: {},
            terms.model.empirical: 3
        },
        terms.model.type: type,
        terms.model.get_branch_length: "",
        terms.model.set_branch_length: "models.generic.SetBranchLength",
        terms.model.constrain_branch_length: "models.generic.constrain_branch_length",
        terms.model.frequency_estimator: "frequencies.empirical.nucleotide",
        terms.model.q_ij: "models.DNA.GTR._generateRate",
        terms.model.time: "models.DNA.generic.Time",
        terms.model.defineQ: "models.DNA.GTR._defineQ",
        terms.model.post_definition: "models.generic.post.definition"
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
    models.DNA.GTR._generateRate.p[terms.model.rate_entry] = models.DNA.GTR.parameter_name;

    return models.DNA.GTR._generateRate.p;
}

/**
 * @name models.DNA.GTR._defineQ
 * @param {Dictionary} gtr
 * @param {String} namespace
 */
function models.DNA.GTR._defineQ(gtr, namespace) {

    models.DNA.generic.DefineQMatrix(gtr, namespace);
    
    if (gtr[terms.model.type] == terms.global) {
        parameters.SetConstraint(((gtr[terms.parameters])[terms.global])[terms.nucleotideRate("A", "G")], "1", "");
    }
    if (gtr[terms.model.type] == terms.local) {
        parameters.SetConstraint(((gtr[terms.parameters])[terms.local])[terms.nucleotideRate("A", "G")], "1", "");
    }
    return gtr;
}
