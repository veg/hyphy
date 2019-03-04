LoadFunctionLibrary("../DNA.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");
LoadFunctionLibrary ("libv3/all-terms.bf");

/** @module models.DNA.MS.GTR */

/**
 * @name models.DNA.MS.GTR.ModelDescription
 */
lfunction models.DNA.MS.GTR.ModelDescription(ignored) {

    return {
        utility.getGlobalValue("terms.alphabet"): utility.getGlobalValue("models.DNA.alphabet"),
        utility.getGlobalValue("terms.description"): "The general time reversible (GTR) model of nucleotide substitution in the mutation selection form",
        utility.getGlobalValue("terms.model.canonical"):  0, // is of the r_ij \times \pi_j form
        utility.getGlobalValue("terms.model.reversible"): 1,
        utility.getGlobalValue("terms.model.efv_estimate_name"): utility.getGlobalValue("terms.frequencies.run_time"),
        utility.getGlobalValue("terms.parameters"): {
            utility.getGlobalValue("terms.global"): {},
            utility.getGlobalValue("terms.local"): {}
        },
        utility.getGlobalValue("terms.model.type"): utility.getGlobalValue("terms.global"),
        utility.getGlobalValue("terms.model.get_branch_length"): "",
        utility.getGlobalValue("terms.model.set_branch_length"): "models.generic.SetBranchLength",
        utility.getGlobalValue("terms.model.constrain_branch_length"): "models.generic.ConstrainBranchLength",
        utility.getGlobalValue("terms.model.frequency_estimator"): "frequencies.runtime.nucleotide",
        utility.getGlobalValue("terms.model.q_ij"): "models.DNA.MS.GTR._generateRate",
        utility.getGlobalValue("terms.model.time"): "models.DNA.generic.Time",
        utility.getGlobalValue("terms.model.defineQ"): "models.DNA.MS.GTR._defineQ",
        utility.getGlobalValue("terms.model.post_definition"): "models.DNA.MS.GTR.post.definition"//"models.generic.post.definition"
    };
}

/**
 * @name models.DNA.MS.GTR._generateRate
 * @param {Number} fromChar
 * @param {Number} toChar
 * @param {String} namespace
 * @param {String} model_type
 * @returns models.DNA.GTR.MS._generateRate.p
 */
 
lfunction models.DNA.MS.GTR._generateRate(fromChar, toChar, namespace, model_type, model) {

    p = {};
    p[utility.getGlobalValue("terms.global")] = {};
    if (fromChar > toChar) {
        parameter_name = "theta_" + toChar + fromChar;
    } else {
        parameter_name = "theta_" + fromChar + toChar;
    }
    parameter_name = parameters.ApplyNameSpace(parameter_name, namespace);
    (p[utility.getGlobalValue("terms.global")])[terms.nucleotideRate(fromChar, toChar)] = parameter_name;
    p[utility.getGlobalValue("terms.model.rate_entry")] = parameter_name;
    return p;
}

/**
 * @name models.DNA.GTR._defineQ
 * @param {Dictionary} gtr
 * @param {String} namespace
 */
function models.DNA.MS.GTR._defineQ(gtr, namespace) {

    models.DNA.generic.DefineQMatrix(gtr, namespace);
    
    if (gtr[terms.model.type] == terms.global) {
        parameters.SetConstraint(((gtr[terms.parameters])[terms.global])[terms.nucleotideRate("A", "G")], "1", "");
    }

    
    return gtr;
}

/**
 * @name models.DNA.GTR._defineQ
 * @param {Dictionary} gtr
 * @param {String} namespace
 */
function models.DNA.MS.GTR.post.definition (gtr, namespace) {

    console.log (gtr);
    
    return models.generic.post.definition (gtr, namespace);
}
