RequireVersion ("2.5.60");

LoadFunctionLibrary("../codon.bf");
LoadFunctionLibrary("../DNA.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");
LoadFunctionLibrary("MG_REV.bf");
LoadFunctionLibrary("../protein.bf");

/** @module models.codon.GTR */
//----------------------------------------------------------------------------------------------------------------


lfunction model.codon.MG_GTR.prompt_and_define (type, code) {
     console.log (type);
     return models.codon.MG_GTR.ModelDescription(type, code);
}

lfunction models.codon.MG_GTR.ModelDescription(type, code) {


    io.CheckAssertion ("`&type` == '`^'terms.global'`'", "MG_GTR only supports `^'terms.global'` type");

    // piggyback on the standard MG_REV model for most of the code

    mg_base = models.codon.MG_REV.ModelDescription (type, code);
    mg_base[utility.getGlobalValue("terms.description")] = "The Muse-Gaut 94 codon-substitution model coupled with the general time reversible (GTR) model of nucleotide substitution, which defines a separate dN/dS ratio for each pair of amino-acids";
    mg_base[utility.getGlobalValue("terms.model.q_ij")] = "models.codon.MG_GTR._GenerateRate";
    mg_base[utility.getGlobalValue("terms.model.post_definition")] = "models.codon.MG_GTR.post_definition";
    return mg_base;
}

lfunction models.codon.MG_GTR._GenerateRate(fromChar, toChar, namespace, model_type, model) {

    return models.codon.MG_GTR._GenerateRate_generic (fromChar, toChar, namespace, model_type,
        model[utility.getGlobalValue("terms.translation_table")],
        "alpha", utility.getGlobalValue("terms.parameters.synonymous_rate"),
        "omega",  utility.getGlobalValue("terms.parameters.omega_ratio"),
    );
}

/**
 * @name models.codon.MG_GTR._GenerateRate
 * @param {Number} fromChar
 * @param {Number} toChar
 * @param {String} namespace
 * @param {String} model_type
 * @param {Matrix} _tt - translation table
 */


lfunction models.codon.MG_GTR._GenerateRate_generic (fromChar, toChar, namespace, model_type, _tt, alpha, alpha_term, omega, omega_term) {

    _GenerateRate.p = {};
    _GenerateRate.diff = models.codon.diff.complete(fromChar, toChar);
    diff_count = utility.Array1D (_GenerateRate.diff);

    if (diff_count == 1) {

        _GenerateRate.p[model_type] = {};
        _GenerateRate.p[utility.getGlobalValue("terms.global")] = {};

        nuc_rate = "";

        for (i = 0; i < diff_count; i += 1) {
            if ((_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")] > (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")]) {
                nuc_p = "theta_" + (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")] + (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")];
            } else {
                nuc_p = "theta_" + (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")] +(_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")];
            }
            nuc_p = parameters.ApplyNameSpace(nuc_p, namespace);
            (_GenerateRate.p[utility.getGlobalValue("terms.global")])[terms.nucleotideRateReversible((_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")], (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")])] = nuc_p;

            nuc_rate = parameters.AppendMultiplicativeTerm (nuc_rate, nuc_p);
       }

        rate_entry = nuc_rate;

        if (_tt[fromChar] != _tt[toChar]) {
            
            if (_tt[fromChar] < _tt[toChar]) {
                aa_pair = _tt[fromChar] + " and " + _tt[toChar];
                aa_tag = "_" + _tt[fromChar] + _tt[toChar];
            } else {
                aa_pair = _tt[toChar] + " and " + _tt[fromChar];
                aa_tag = "_" + _tt[toChar] + _tt[fromChar];
            }
            
            omega_p = parameters.ApplyNameSpace (omega + aa_tag , namespace);
            (_GenerateRate.p[model_type]) [^"terms.parameters.omega_ratio" + " for " + aa_pair] = omega_p;
                
            rate_entry = parameters.AppendMultiplicativeTerm (rate_entry, omega_p);
        } else {
            _GenerateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate;
        }

        _GenerateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = rate_entry;
    }
    return _GenerateRate.p;
}

lfunction  models.codon.MG_GTR.post_definition (model) {

    for (id; in ; model.GetParameters_RegExp (model,  ^"terms.parameters.omega_ratio")) {
        parameters.SetValue(id, 0.25);
    }

    models.generic.post.definition (model);
}

lfunction models.codon.MG_GTR.set_branch_length(model, value, parameter) {
    return models.codon.MG_REV.set_branch_length(model,value,parameter);
}
