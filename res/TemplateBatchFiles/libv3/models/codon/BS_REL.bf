LoadFunctionLibrary("../codon.bf");
LoadFunctionLibrary("../DNA.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");
LoadFunctionLibrary("MG_REV.bf");
LoadFunctionLibrary("../../convenience/math.bf");

/** @module models.codon.BS_REL

    Define libv3 style branch-site REL models with N mixture components

*/

/**
 * @name models.codon.BS_REL.ModelDescription
 * @param {String} type
 * @param {String} code
 * @param {Number} components (>=2)
 */
lfunction models.codon.BS_REL.ModelDescription(type, code, components) {

    components = math.Int (components);

    io.CheckAssertion ("`&components` >= 2 && `&components` <= 10", "must have between 2 and 10 components in call to models.codon.BS_REL.ModelDescription");

    codons = models.codon.MapCode(code);

    return {
        "alphabet": codons["sense"],
        "bases": ^ "models.DNA.alphabet",
        "components" : components,
        "stop": codons["stop"],
        "type": type,
        "translation-table": codons["translation-table"],
        "description": "The branch-site mixture of N Muse-Gaut 94 codon-substitution model coupled with the general time reversible (GTR) model of nucleotide substitution",
        "canonical": "EXPLICIT_FORM_MATRIX_EXPONENTIAL", // is NOT of the r_ij \times \pi_j form
        "reversible": TRUE,
        ^ "terms.efv_estimate_name": ^ "terms.freqs.CF3x4",
        "parameters": {
            "global": {},
            "local":  {}
        },
        "get-branch-length": "models.codon.BS_REL.get_branch_length",
        "set-branch-length": "models.codon.BS_REL.set_branch_length",
        "constrain-branch-length": "models.codon.BS_REL.constrain_branch_length",
        "frequency-estimator": "frequencies.empirical.corrected.CF3x4",
        "q_ij": "", // set for individual components in models.codon.BS_REL._DefineQ
        "time": "models.DNA.generic.Time",
        "defineQ": "models.codon.BS_REL._DefineQ",
        "post-definition": "models.generic.post.definition"
    };
}

/**
 * @name models.codon.BS_REL._DefineQ
 * @param {Dict} mg_rev
 * @param {String} namespace
 * @returns {Dict} updated model
 */

lfunction models.codon.BS_REL._DefineQ(bs_rel, namespace) {

    rate_matrices = {};


    bs_rel ["q_ij"] = &rate_generator;
    bs_rel [^'terms.mixture'] = {};

    _aux = parameters.GenerateSequentialNames (namespace + "_mixture_aux", bs_rel["components"] - 1, "_");
    parameters.DeclareGlobal (_aux, None);
    _wts = parameters.helper.stick_breaking (_aux, None);
    mixture = {};


    for (component = 1; component <= bs_rel["components"]; component += 1) {
       key = "component_" + component;
       ExecuteCommands ("
        function rate_generator (fromChar, toChar, namespace, model_type, _tt) {
            return models.codon.MG_REV._GenerateRate_generic (fromChar, toChar, namespace, model_type, _tt,
                'alpha', ^'terms.synonymous_rate',
                'beta_`component`', terms.AddCategory (^'terms.nonsynonymous_rate', component),
                'omega`component`', terms.AddCategory (^'terms.omega_ratio', component));
        }"
       );
       models.codon.generic.DefineQMatrix(bs_rel, namespace);
       rate_matrices [key] = bs_rel[^"terms.rate_matrix"];
       (bs_rel [^'terms.mixture'])[key] = _wts [component-1];

       if ( component < bs_rel["components"]) {
            model.generic.AddGlobal ( bs_rel, _aux[component-1], terms.AddCategory ("Stick-breaking proportion", component ));
            parameters.SetRange (_aux[component-1], ^"terms.range_almost_01");
        }
         mixture [key] = namespace + "_mixture_weight_" + component;
        model.generic.AddGlobal ( bs_rel, mixture [key], terms.AddCategory ("Mixture weight", component));
        parameters.SetConstraint (mixture [key], _wts[component-1], "global ");

    }

    bs_rel[^"terms.rate_matrix"] = rate_matrices;
    bs_rel[^"terms.mixture"] = mixture;


    parameters.SetConstraint(((bs_rel["parameters"])[^'terms.global'])[terms.nucleotideRate("A", "G")], "1", "");
    return bs_rel;
}


/**
 * @name models.codon.MG_REV.set_branch_length
 * @param {Model} model
 * @param {AssociativeList|Number} value
 * @param {String} parameter
 * @returns {Number} 0
 */
function models.codon.BS_REL.set_branch_length(model, value, parameter) {
    if (model["type"] == terms.global) {
        return models.generic.SetBranchLength(model, value, parameter);
    }


    models.codon.MG_REV.set_branch_length.beta = model.generic.GetLocalParameter(model, terms.nonsynonymous_rate);
    models.codon.MG_REV.set_branch_length.alpha = model.generic.GetLocalParameter(model, terms.synonymous_rate);

    models.codon.MG_REV.set_branch_length.alpha.p = parameter + "." + models.codon.MG_REV.set_branch_length.alpha;
    models.codon.MG_REV.set_branch_length.beta.p = parameter + "." + models.codon.MG_REV.set_branch_length.beta;



    if (Type(value) == "AssociativeList") {
        if (value[terms.branch_length_scaler] == terms.branch_length_constrain) {
            if (parameters.IsIndependent(models.codon.MG_REV.set_branch_length.alpha.p)) {
                if (Abs(model[terms.synonymous_rate]) == 0) {
                    bl_string = model["branch-length-string"];

                    models.codon.MG_REV.set_branch_length.sub = {};
                    models.codon.MG_REV.set_branch_length.sub[models.codon.MG_REV.set_branch_length.alpha] = 1;
                    models.codon.MG_REV.set_branch_length.sub[models.codon.MG_REV.set_branch_length.beta] = 0;
                    model[terms.synonymous_rate] = Simplify(bl_string, models.codon.MG_REV.set_branch_length.sub);

                    models.codon.MG_REV.set_branch_length.sub[models.codon.MG_REV.set_branch_length.alpha] = 0;
                    models.codon.MG_REV.set_branch_length.sub[models.codon.MG_REV.set_branch_length.beta] = 1;
                    model[terms.nonsynonymous_rate] = Simplify(bl_string, models.codon.MG_REV.set_branch_length.sub);
                }

                //fprintf (stdout, models.codon.MG_REV.set_branch_length.alpha.p, "\n");

                parameters.SetConstraint(models.codon.MG_REV.set_branch_length.beta.p, "(" + 3 * value[terms.branch_length] + " - " + models.codon.MG_REV.set_branch_length.alpha.p + "*(" + model[terms.synonymous_rate] + "))/(" + model[terms.nonsynonymous_rate] + ")", "");
                return 1;

            } else {
                assert(0, "TBA in models.codon.MG_REV.set_branch_length");
            }
        } else {

            models.codon.MG_REV.set_branch_length.lp = 0;
            if (parameters.IsIndependent(models.codon.MG_REV.set_branch_length.alpha.p)) {
                Eval(models.codon.MG_REV.set_branch_length.alpha.p + ":=(" + value[terms.branch_length_scaler] + ")*" + value[terms.branch_length]);
                models.codon.MG_REV.set_branch_length.lp += 1;
            }
            if (parameters.IsIndependent(models.codon.MG_REV.set_branch_length.beta.p)) {
                Eval(models.codon.MG_REV.set_branch_length.beta.p + ":=(" + value[terms.branch_length_scaler] + ")*" + value[terms.branch_length]);
                models.codon.MG_REV.set_branch_length.lp += 1;
            }
            return models.codon.MG_REV.set_branch_length.lp;
        }
    } else {
        if (parameters.IsIndependent(models.codon.MG_REV.set_branch_length.alpha.p)) {
            if (parameters.IsIndependent(models.codon.MG_REV.set_branch_length.beta.p)) {
                models.codon.MG_REV.set_branch_length.lp = parameters.NormalizeRatio(Eval(models.codon.MG_REV.set_branch_length.beta), Eval(models.codon.MG_REV.set_branch_length.alpha));
                parameters.SetConstraint(models.codon.MG_REV.set_branch_length.beta, models.codon.MG_REV.set_branch_length.alpha + "*" + models.codon.MG_REV.set_branch_length.lp, "");
                ExecuteCommands("FindRoot (models.codon.MG_REV.set_branch_length.lp,(" + model["branch-length-string"] + ")-" + value + "," + models.codon.MG_REV.set_branch_length.alpha + ",0,10000)");
                parameters.RemoveConstraint(models.codon.MG_REV.set_branch_length.beta);
                Eval("`models.codon.MG_REV.set_branch_length.alpha.p` =" + models.codon.MG_REV.set_branch_length.lp);
                Eval("`models.codon.MG_REV.set_branch_length.beta.p` =" + Eval(models.codon.MG_REV.set_branch_length.beta.p));
            } else {
                ExecuteCommands("FindRoot (models.codon.MG_REV.set_branch_length.lp,(" + model["branch-length-string"] + ")-" + value + "," + models.codon.MG_REV.set_branch_length.alpha + ",0,10000)");
                Eval("`models.codon.MG_REV.set_branch_length.alpha.p` =" + models.codon.MG_REV.set_branch_length.lp);
            }
        }
    }

    return 0;
}
