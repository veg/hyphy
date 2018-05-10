LoadFunctionLibrary("../codon.bf");
LoadFunctionLibrary("../DNA.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");

/** @module models.codon.MG_REV */

/**
 * @name models.codon.MG_REV.ModelDescription
 * @param {String} type
 * @param {String} code
 */
lfunction models.codon.MG_REV.ModelDescription(type, code) {
    codons = models.codon.MapCode(code);

    return {
        utility.getGlobalValue("terms.alphabet"): codons[utility.getGlobalValue("terms.sense_codons")],
        utility.getGlobalValue("terms.bases"): utility.getGlobalValue("models.DNA.alphabet"),
        utility.getGlobalValue("terms.stop_codons"): codons[utility.getGlobalValue("terms.stop_codons")],
        utility.getGlobalValue("terms.model.type"): type,
        utility.getGlobalValue("terms.translation_table"): codons[utility.getGlobalValue("terms.translation_table")],
        utility.getGlobalValue("terms.description"): "The Muse-Gaut 94 codon-substitution model coupled with the general time reversible (GTR) model of nucleotide substitution",
        utility.getGlobalValue("terms.model.canonical"): 0, // is NOT of the r_ij \times \pi_j form
        utility.getGlobalValue("terms.model.reversible"): 1,
        utility.getGlobalValue("terms.model.efv_estimate_name"): utility.getGlobalValue("terms.frequencies.CF3x4"),
        utility.getGlobalValue("terms.parameters"): {
            utility.getGlobalValue("terms.global"): {},
            utility.getGlobalValue("terms.local"): {}
        },
        utility.getGlobalValue("terms.model.get_branch_length"): "",
        utility.getGlobalValue("terms.model.set_branch_length"): "models.codon.MG_REV.set_branch_length",
        utility.getGlobalValue("terms.model.constrain_branch_length"): "models.generic.ConstrainBranchLength",
        utility.getGlobalValue("terms.model.frequency_estimator"): "frequencies.empirical.corrected.CF3x4",
        utility.getGlobalValue("terms.model.q_ij"): "models.codon.MG_REV._GenerateRate",
        utility.getGlobalValue("terms.model.time"): "models.DNA.generic.Time",
        utility.getGlobalValue("terms.model.defineQ"): "models.codon.MG_REV._DefineQ",
        utility.getGlobalValue("terms.model.post_definition"): "models.generic.post.definition"
    };
}


lfunction models.codon.MG_REV._GenerateRate(fromChar, toChar, namespace, model_type, model) {
    return models.codon.MG_REV._GenerateRate_generic (fromChar, toChar, namespace, model_type,
    model[utility.getGlobalValue("terms.translation_table")],
    "alpha", utility.getGlobalValue("terms.parameters.synonymous_rate"), "beta", utility.getGlobalValue("terms.parameters.nonsynonymous_rate"), "omega", utility.getGlobalValue("terms.parameters.omega_ratio"));
}

/**
 * @name models.codon.MG_REV._GenerateRate_generic
 * @param {Number} fromChar
 * @param {Number} toChar
 * @param {String} namespace
 * @param {String} model_type
 * @param {Matrix} _tt - translation table
 */


lfunction models.codon.MG_REV._GenerateRate_generic (fromChar, toChar, namespace, model_type, _tt, alpha, alpha_term, beta, beta_term, omega, omega_term) {

    _GenerateRate.p = {};
    _GenerateRate.diff = models.codon.diff(fromChar, toChar);

    if (None != _GenerateRate.diff) {
        _GenerateRate.p[model_type] = {};
        _GenerateRate.p[utility.getGlobalValue("terms.global")] = {};

        if (_GenerateRate.diff[utility.getGlobalValue("terms.diff.from")] > _GenerateRate.diff[utility.getGlobalValue("terms.diff.to")]) {
            nuc_rate = "theta_" + _GenerateRate.diff[utility.getGlobalValue("terms.diff.to")] + _GenerateRate.diff[utility.getGlobalValue("terms.diff.from")];
        } else {
            nuc_rate = "theta_" + _GenerateRate.diff[utility.getGlobalValue("terms.diff.from")] + _GenerateRate.diff[utility.getGlobalValue("terms.diff.to")];
        }

        nuc_rate = parameters.ApplyNameSpace(nuc_rate, namespace);
        (_GenerateRate.p[utility.getGlobalValue("terms.global")])[terms.nucleotideRateReversible(_GenerateRate.diff[utility.getGlobalValue("terms.diff.from")], _GenerateRate.diff[utility.getGlobalValue("terms.diff.to")])] = nuc_rate;

        if (_tt[fromChar] != _tt[toChar]) {
            if (model_type == utility.getGlobalValue("terms.global")) {
                aa_rate = parameters.ApplyNameSpace(omega, namespace);
                (_GenerateRate.p[model_type])[omega_term] = aa_rate;
            } else {
                aa_rate = beta;
                (_GenerateRate.p[model_type])[beta_term] = aa_rate;
            }
            _GenerateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate + "*" + aa_rate;
        } else {
            if (model_type == utility.getGlobalValue("terms.local")) {
                (_GenerateRate.p[model_type])[alpha_term] = alpha;
                _GenerateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate + "*" + alpha;
            } else {
                _GenerateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate;
            }
        }
    }


    return _GenerateRate.p;
}

/**
 * @name models.codon.MG_REV._DefineQ
 * @param {Model} mg_rev
 * @param {String} namespace
 * @returns {Model} updated model
 */
function models.codon.MG_REV._DefineQ(mg_rev, namespace) {
    models.codon.generic.DefineQMatrix(mg_rev, namespace);
    parameters.SetConstraint(((mg_rev[terms.parameters])[terms.global])[terms.nucleotideRate("A", "G")], "1", "");
    return mg_rev;
}

/**
 * @name models.codon.MG_REV.set_branch_length
 * @param {Model} model
 * @param {AssociativeList|Number} value
 * @param {String} parameter
 * @returns {Number} 0
 */
function models.codon.MG_REV.set_branch_length(model, value, parameter) {
    if (model[terms.model.type] == terms.global) {
        return models.generic.SetBranchLength(model, value, parameter);
    }


    models.codon.MG_REV.set_branch_length.beta = model.generic.GetLocalParameter(model, terms.parameters.nonsynonymous_rate);
    models.codon.MG_REV.set_branch_length.alpha = model.generic.GetLocalParameter(model, terms.parameters.synonymous_rate);

    models.codon.MG_REV.set_branch_length.alpha.p = parameter + "." + models.codon.MG_REV.set_branch_length.alpha;
    models.codon.MG_REV.set_branch_length.beta.p = parameter + "." + models.codon.MG_REV.set_branch_length.beta;



    if (Type(value) == "AssociativeList") {
        if (value[terms.model.branch_length_scaler] == terms.model.branch_length_constrain) {
            if (parameters.IsIndependent(models.codon.MG_REV.set_branch_length.alpha.p)) {
                if (Abs(model[terms.parameters.synonymous_rate]) == 0) {
                    bl_string = model[terms.model.branch_length_string];

                    models.codon.MG_REV.set_branch_length.sub = {};
                    models.codon.MG_REV.set_branch_length.sub[models.codon.MG_REV.set_branch_length.alpha] = 1;
                    models.codon.MG_REV.set_branch_length.sub[models.codon.MG_REV.set_branch_length.beta] = 0;
                    model[terms.parameters.synonymous_rate] = Simplify(bl_string, models.codon.MG_REV.set_branch_length.sub);

                    models.codon.MG_REV.set_branch_length.sub[models.codon.MG_REV.set_branch_length.alpha] = 0;
                    models.codon.MG_REV.set_branch_length.sub[models.codon.MG_REV.set_branch_length.beta] = 1;
                    model[terms.parameters.nonsynonymous_rate] = Simplify(bl_string, models.codon.MG_REV.set_branch_length.sub);
                }

                //fprintf (stdout, models.codon.MG_REV.set_branch_length.alpha.p, "\n");

                parameters.SetConstraint(models.codon.MG_REV.set_branch_length.beta.p, "(" + 3 * value[terms.branch_length] + " - " + models.codon.MG_REV.set_branch_length.alpha.p + "*(" + model[terms.parameters.synonymous_rate] + "))/(" + model[terms.parameters.nonsynonymous_rate] + ")", "");
                return 1;

            } else {
                assert(0, "TBA in models.codon.MG_REV.set_branch_length");
            }
        } else {
            models.codon.MG_REV.set_branch_length.lp = 0;
            if (parameters.IsIndependent(models.codon.MG_REV.set_branch_length.alpha.p)) {
                Eval(models.codon.MG_REV.set_branch_length.alpha.p + ":=(" + value[terms.model.branch_length_scaler] + ")*" + value[terms.branch_length]);
                models.codon.MG_REV.set_branch_length.lp += 1;
            }
            if (parameters.IsIndependent(models.codon.MG_REV.set_branch_length.beta.p)) {
                Eval(models.codon.MG_REV.set_branch_length.beta.p + ":=(" + value[terms.model.branch_length_scaler] + ")*" + value[terms.branch_length]);
                models.codon.MG_REV.set_branch_length.lp += 1;
            }
            return models.codon.MG_REV.set_branch_length.lp;
        }
    } else {
        if (parameters.IsIndependent(models.codon.MG_REV.set_branch_length.alpha.p)) {
            if (parameters.IsIndependent(models.codon.MG_REV.set_branch_length.beta.p)) {
                models.codon.MG_REV.set_branch_length.lp = parameters.NormalizeRatio(Eval(models.codon.MG_REV.set_branch_length.beta.p), Eval(models.codon.MG_REV.set_branch_length.alpha.p));
                parameters.SetConstraint(models.codon.MG_REV.set_branch_length.beta, models.codon.MG_REV.set_branch_length.alpha + "*" + models.codon.MG_REV.set_branch_length.lp, "");
                ExecuteCommands("FindRoot (models.codon.MG_REV.set_branch_length.lp,(" + model[terms.model.branch_length_string] + ")-" + 3*value + "," + models.codon.MG_REV.set_branch_length.alpha + ",0,10000)");
                Eval("`models.codon.MG_REV.set_branch_length.alpha.p` =" + models.codon.MG_REV.set_branch_length.lp);
                parameters.RemoveConstraint(models.codon.MG_REV.set_branch_length.beta);
                Eval ("`models.codon.MG_REV.set_branch_length.beta.p` =" + Eval(models.codon.MG_REV.set_branch_length.beta));
            } else {
                parameters.SetConstraint (models.codon.MG_REV.set_branch_length.beta, parameters.GetConstraint (models.codon.MG_REV.set_branch_length.alpha.p),"");
                parameters.SetConstraint (models.codon.MG_REV.set_branch_length.alpha, models.codon.MG_REV.set_branch_length.alpha.p,"");
                ExecuteCommands("FindRoot (models.codon.MG_REV.set_branch_length.lp,(" + model[terms.model.branch_length_string] + ")-" + 3*value + "," + models.codon.MG_REV.set_branch_length.alpha + ",0,10000)");
                Eval("`models.codon.MG_REV.set_branch_length.alpha.p` =" + models.codon.MG_REV.set_branch_length.lp);
                messages.log ("models.codon.MG_REV.set_branch_length: " + models.codon.MG_REV.set_branch_length.alpha.p + "=" + models.codon.MG_REV.set_branch_length.lp);
                parameters.RemoveConstraint (models.codon.MG_REV.set_branch_length.alpha);
                parameters.RemoveConstraint (models.codon.MG_REV.set_branch_length.beta);
            }
        }
    }

    return 0;
}
