LoadFunctionLibrary("MG_REV.bf");

/** @module models.codon.MG_REV */

/**
 * @name models.codon.MSS.ModelDescription
 * @param {String} type
 * @param {String} code
 * @param {Dict} codon_classes : codon => neutral class
  */
  
terms.model.MSS.syn_rate_within  = " within codon class ";
terms.model.MSS.syn_rate_between = " between codon classes ";

  
lfunction models.codon.MSS.ModelDescription(type, code, codon_classes) {
    m = Call ("models.codon.MG_REV.ModelDescription", type, code);
    if (^"fitter.frequency_type" == "F3x4") {
        m[utility.getGlobalValue("terms.model.frequency_estimator")] = "frequencies.empirical.F3x4";
    } else {
        if (^"fitter.frequency_type" == "F1x4") {
            m[utility.getGlobalValue("terms.model.frequency_estimator")] = "frequencies.empirical.F1x4";
        }
    }
    m[utility.getGlobalValue("terms.description")] = "The Muse-Gaut 94 codon-substitution model coupled with the general time reversible (GTR) model of nucleotide substitution, which allows multiple classes of synonymous substitution rates";
    m[utility.getGlobalValue("terms.model.q_ij")] = "models.codon.MSS._GenerateRate";
    m[utility.getGlobalValue("terms.model.mss.codon_classes")] = codon_classes;
    return m;
}

//----------------------------------------------------------------------------------------------------------------

lfunction models.codon.MSS._GenerateRate (fromChar, toChar, namespace, model_type, model) {

    _GenerateRate.p = {};
    _GenerateRate.diff = models.codon.diff.complete(fromChar, toChar);
    diff_count = utility.Array1D (_GenerateRate.diff);

    omega_term = utility.getGlobalValue ("terms.parameters.omega_ratio");
    alpha_term = utility.getGlobalValue ("terms.parameters.synonymous_rate");
    beta_term  = utility.getGlobalValue ("terms.parameters.nonsynonymous_rate");
    omega      = "omega";
    alpha      = "alpha";
    beta       = "beta";

    _tt = model[utility.getGlobalValue("terms.translation_table")];

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

            if (model_type == utility.getGlobalValue("terms.global")) {
                aa_rate = parameters.ApplyNameSpace(omega, namespace);
                (_GenerateRate.p[model_type])[omega_term] = aa_rate;
            } else {
                aa_rate = beta;
                (_GenerateRate.p[model_type])[beta_term] = aa_rate;
            }
            rate_entry += "*" + aa_rate;
        } else {

            class_from = (model[^"terms.model.mss.codon_classes"])[fromChar];
            class_to   = (model[^"terms.model.mss.codon_classes"])[toChar];

            if (class_from == class_to) {
                if (class_from == ^"mss.neutral_reference") {
                    if (model_type == utility.getGlobalValue("terms.local")) {
                        codon_rate = alpha + "_" + class_from;
                        (_GenerateRate.p[model_type])[alpha_term + ^"terms.model.MSS.syn_rate_within" + class_from] = codon_rate;
                        rate_entry += "*" + codon_rate;
                    } else {
                        rate_entry = nuc_rate;
                    }
                } else {
                    if (model_type == utility.getGlobalValue("terms.local")) {
                        codon_rate = alpha + "_" + class_from;
                    } else {
                        codon_rate = parameters.ApplyNameSpace(alpha + "_" + class_from, namespace);
                    }
                    (_GenerateRate.p[model_type])[alpha_term + ^"terms.model.MSS.syn_rate_within" + class_from] = codon_rate;
                    rate_entry += "*" + codon_rate;
                }
            } else {
                if (class_from > class_to) {
                    codon_rate = class_to;
                    class_to = class_from;
                    class_from = codon_rate;
                }
                if (model_type == utility.getGlobalValue("terms.local")) {
                    codon_rate = alpha + "_" + class_from + "_" + class_to;
                } else {
                    codon_rate = parameters.ApplyNameSpace(alpha + "_" + class_from + "_" + class_to, namespace);
                }
                (_GenerateRate.p[model_type])[alpha_term + ^"terms.model.MSS.syn_rate_between" + class_from + " and "  + class_to] = codon_rate;
                rate_entry += "*" + codon_rate;
            }
        }

        _GenerateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = rate_entry;
    }

    return _GenerateRate.p;
}
 
