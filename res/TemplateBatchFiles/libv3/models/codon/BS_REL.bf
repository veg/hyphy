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
    io.CheckAssertion ("`&components` >= 1 && `&components` <= 10", "must have between 1 and 10 components in call to models.codon.BS_REL.ModelDescription");
    codons = models.codon.MapCode(code);

    return {
        utility.getGlobalValue("terms.alphabet"): codons[utility.getGlobalValue("terms.sense_codons")],
        utility.getGlobalValue("terms.bases"): utility.getGlobalValue("models.DNA.alphabet"),
        utility.getGlobalValue("terms.model.components"): components,
        utility.getGlobalValue("terms.stop_codons"): codons[utility.getGlobalValue("terms.stop_codons")],
        utility.getGlobalValue("terms.model.type"): type,
        utility.getGlobalValue("terms.translation_table"): codons[utility.getGlobalValue("terms.translation_table")],
        utility.getGlobalValue("terms.description"): "The branch-site mixture of N Muse-Gaut 94 codon-substitution model coupled with the general time reversible (GTR) model of nucleotide substitution",
        utility.getGlobalValue("terms.model.canonical"): "EXPLICIT_FORM_MATRIX_EXPONENTIAL", // is NOT of the r_ij \times \pi_j form
        utility.getGlobalValue("terms.model.reversible"): TRUE,
        utility.getGlobalValue("terms.model.efv_estimate_name"): utility.getGlobalValue("terms.frequencies.CF3x4"),
        utility.getGlobalValue("terms.parameters"): {
            utility.getGlobalValue("terms.global"): {},
            utility.getGlobalValue("terms.local"):  {}
        },
        utility.getGlobalValue("terms.model.get_branch_length"): "models.codon.BS_REL.get_branch_length",
        utility.getGlobalValue("terms.model.set_branch_length"): "models.codon.BS_REL.set_branch_length",
        utility.getGlobalValue("terms.model.constrain_branch_length"): "models.codon.BS_REL.constrain_branch_length",
        utility.getGlobalValue("terms.model.frequency_estimator"): "frequencies.empirical.corrected.CF3x4",
        utility.getGlobalValue("terms.model.q_ij"): "", // set for individual components in models.codon.BS_REL._DefineQ
        utility.getGlobalValue("terms.model.time"): "models.DNA.generic.Time",
        utility.getGlobalValue("terms.model.defineQ"): "models.codon.BS_REL._DefineQ",
        utility.getGlobalValue("terms.model.post_definition"): "models.codon.BS_REL.post_definition"
    };
}



/**
 * @name models.codon.BS_REL_Per_Branch_Mixing.ModelDescription
 * @param {String} type
 * @param {String} code
 * @param {Number} components (>=2)
 */
lfunction models.codon.BS_REL_Per_Branch_Mixing.ModelDescription(type, code, components) {
	template = models.codon.BS_REL.ModelDescription(type, code, components);
	template [utility.getGlobalValue("terms.model.defineQ")] = "models.codon.BS_REL_Per_Branch_Mixing._DefineQ";
	return template;
}

/**
 * @name models.codon.BS_REL_SRV.ModelDescription
 * @param {String} type
 * @param {String} code
 * @param {Matrix} components (alpha, beta) components
 */
lfunction models.codon.BS_REL_SRV.ModelDescription(type, code, components) {

    io.CheckAssertion ('`&type`==terms.global', 'Only ' + ^'terms.global' + ' model type is supported for BS_REL_SRV');

    io.CheckAssertion ('Type (`&components`) == "Matrix" && utility.Array1D (`&components`) == 2', "must have a 2 dimensional matrix of rate counts in models.codon.BS_REL_SRV.ModelDescription");
    components = utility.Map (components, "_value_", "math.Int(_value_)");
    io.CheckAssertion ("Min(`&components`,0) >= 1 && Max(`&components`,0) <= 10", "must have between 1 and 10 components in call to models.codon.BS_REL_SRV.ModelDescription");


 	template = models.codon.BS_REL.ModelDescription(type, code, 3);
	template [utility.getGlobalValue("terms.model.defineQ")] = "models.codon.BS_REL_SRV._DefineQ";
	template [utility.getGlobalValue("terms.model.components")] = components;
	return template;
}

/**
 * @name models.codon.BS_REL.BS_REL_Per_Branch_Mixing._DefineQ
 * @param {Dict} mg_rev
 * @param {String} namespace
 * @returns {Dict} updated model
 */

lfunction models.codon.BS_REL_Per_Branch_Mixing._DefineQ(bs_rel, namespace) {
    rate_matrices = {};

    bs_rel [utility.getGlobalValue("terms.model.q_ij")] = &rate_generator;
    bs_rel [utility.getGlobalValue("terms.mixture.mixture_components")] = {};


    _aux = parameters.GenerateSequentialNames ("bsrel_mixture_aux", bs_rel[utility.getGlobalValue("terms.model.components")] - 1, "_");
    _wts = parameters.helper.stick_breaking (_aux, None);
    mixture = {};


    for (component = 1; component <= bs_rel[utility.getGlobalValue("terms.model.components")]; component += 1) {
       key = "component_" + component;
       ExecuteCommands ("
        function rate_generator (fromChar, toChar, namespace, model_type, model) {
            return models.codon.MG_REV._GenerateRate_generic (fromChar, toChar, namespace, model_type, model[utility.getGlobalValue('terms.translation_table')],
                'alpha', utility.getGlobalValue('terms.parameters.synonymous_rate'),
                'beta_`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.nonsynonymous_rate'), component),
                'omega`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component));
        }"
       );
       models.codon.generic.DefineQMatrix(bs_rel, namespace);
       rate_matrices [key] = bs_rel[utility.getGlobalValue("terms.model.rate_matrix")];
       (bs_rel [utility.getGlobalValue("terms.mixture.mixture_components")])[key] = _wts [component-1];

       if ( component < bs_rel[utility.getGlobalValue("terms.model.components")]) {
            model.generic.AddLocal ( bs_rel, _aux[component-1], terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), component ));
            parameters.SetRange (_aux[component-1], utility.getGlobalValue("terms.range_almost_01"));
        }
    }

    bs_rel[utility.getGlobalValue("terms.model.rate_matrix")] = rate_matrices;
    parameters.SetConstraint(((bs_rel[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.global")])[terms.nucleotideRate("A", "G")], "1", "");
    return bs_rel;
}

/**
 * @name models.codon.BS_REL.ExtractMixtureDistribution
 * @param {Dict} bs_rel
 * @returns {Dict} mixture distribution parameters
 */

lfunction models.codon.BS_REL.ExtractMixtureDistribution (bs_rel) {
    count = bs_rel [utility.getGlobalValue ("terms.model.components")];

    if (Type (count) == "Matrix") {
        count = count[1];
    }

    rates = {count, 1};
    weights = {count-1, 1};

    for (i = 1; i <= count; i+=1) {
        rates [i-1] = ((bs_rel[utility.getGlobalValue ("terms.parameters")])[utility.getGlobalValue ("terms.global")])[terms.AddCategory (utility.getGlobalValue ("terms.parameters.omega_ratio"), i)];
        if (i < count ) {
            weights [i-1] = ((bs_rel[utility.getGlobalValue ("terms.parameters")])[utility.getGlobalValue ("terms.global")])[terms.AddCategory (utility.getGlobalValue ("terms.mixture.mixture_aux_weight"), i )];
        }
    }


    return {"rates" : rates, "weights" : weights };
}

/**
 * @name models.codon.BS_REL.ExtractMixtureDistribution
 * @param {Dict} bs_rel
 * @returns {Dict} mixture distribution parameters
 */

lfunction models.codon.BS_REL.ExtractMixtureDistributionFromFit (bs_rel, fit) {
    count = bs_rel [utility.getGlobalValue ("terms.model.components")];
    rates = {count, 1};
    weights = {count-1, 1};

    for (i = 1; i <= count; i+=1) {
        rates [i-1] = estimators.GetGlobalMLE (fit, terms.AddCategory (utility.getGlobalValue ("terms.parameters.omega_ratio"), i));
        if (i < count ) {
            weights [i-1] = estimators.GetGlobalMLE (fit, terms.AddCategory (utility.getGlobalValue ("terms.mixture.mixture_aux_weight"), i ));
        }
    }

    return {utility.getGlobalValue ("terms.parameters.rates") : rates, utility.getGlobalValue ("terms.parameters.weights") : weights };
}

/**
 * @name models.codon.BS_REL.BS_REL._DefineQ
 * @param {Dict} mg_rev
 * @param {String} namespace
 * @returns {Dict} updated model
 */

lfunction models.codon.BS_REL._DefineQ(bs_rel, namespace) {

    rate_matrices = {};

    bs_rel [utility.getGlobalValue("terms.model.q_ij")] = &rate_generator;
    bs_rel [utility.getGlobalValue("terms.mixture.mixture_components")] = {};

    _aux = parameters.GenerateSequentialNames (namespace + ".bsrel_mixture_aux", bs_rel[utility.getGlobalValue("terms.model.components")] - 1, "_");
    _wts = parameters.helper.stick_breaking (_aux, None);
    mixture = {};


    for (component = 1; component <= bs_rel[utility.getGlobalValue("terms.model.components")]; component += 1) {
       key = "component_" + component;
       ExecuteCommands ("
        function rate_generator (fromChar, toChar, namespace, model_type, model) {
               return models.codon.MG_REV._GenerateRate_generic (fromChar, toChar, namespace, model_type, model[utility.getGlobalValue('terms.translation_table')],
                'alpha', utility.getGlobalValue('terms.parameters.synonymous_rate'),
                'beta_`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.nonsynonymous_rate'), component),
                'omega`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component));
            }"
       );

       if ( component < bs_rel[utility.getGlobalValue("terms.model.components")]) {
            model.generic.AddGlobal ( bs_rel, _aux[component-1], terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), component ));
            parameters.DeclareGlobalWithRanges (_aux[component-1], 0.5, 0, 1);
       }
       models.codon.generic.DefineQMatrix(bs_rel, namespace);
       rate_matrices [key] = bs_rel[utility.getGlobalValue("terms.model.rate_matrix")];
       (bs_rel [^'terms.mixture.mixture_components'])[key] = _wts [component-1];
    }


    bs_rel[utility.getGlobalValue("terms.model.rate_matrix")] = rate_matrices;
    parameters.SetConstraint(((bs_rel[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.global")])[terms.nucleotideRate("A", "G")], "1", "");

    return bs_rel;
}

/**
 * @name models.codon.BS_REL.BS_REL_SRV._DefineQ
 * @param {Dict} mg_rev
 * @param {String} namespace
 * @returns {Dict} updated model
 */

lfunction models.codon.BS_REL_SRV._DefineQ(bs_rel, namespace) {

    rate_matrices = {};

    bs_rel [utility.getGlobalValue("terms.model.q_ij")] = &rate_generator;
    bs_rel [utility.getGlobalValue("terms.model.time")] = &rate_multiplier;

    bs_rel [utility.getGlobalValue("terms.mixture.mixture_components")] = {};

    ns_components  = (bs_rel[(utility.getGlobalValue("terms.model.components"))])[1];
    syn_components = (bs_rel[(utility.getGlobalValue("terms.model.components"))])[0];


    _aux = parameters.GenerateSequentialNames (namespace + ".bsrel_mixture_aux", ns_components - 1, "_");
    _wts = parameters.helper.stick_breaking (_aux, None);

    _aux_srv = parameters.GenerateSequentialNames (namespace + ".bsrel_mixture_aux_srv", syn_components - 1, "_");
    _wts_srv = parameters.helper.stick_breaking (_aux_srv, None);

    _alphas = parameters.GenerateSequentialNames (namespace + ".alpha", syn_components, "_");


    mixture = {};


    for (s_component = 1; s_component <= syn_components; s_component += 1) {

        ExecuteCommands ("
            function rate_multiplier (option) {
                return {
                    ^'terms.global' : {
                        terms.AddCategory (utility.getGlobalValue('terms.parameters.synonymous_rate'), `s_component`) : '`_alphas[s_component-1]`'
                    },
                    ^'terms.local' : {
                        ^'terms.parameters.synonymous_rate' : models.DNA.generic.Time (null)
                    },
                    ^'terms.model.rate_entry' :  '`_alphas[s_component-1]`*' + models.DNA.generic.Time (null)
                };
           }
        ");

        for (component = 1; component <= ns_components; component += 1) {
           key = "component_" + s_component + "_" + component;



           ExecuteCommands ("
            function rate_generator (fromChar, toChar, namespace, model_type, model) {

                   return models.codon.MG_REV._GenerateRate_generic (fromChar, toChar, namespace, model_type, model[utility.getGlobalValue('terms.translation_table')],
                    'alpha_`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.synonymous_rate'), s_component),
                    'beta_`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.nonsynonymous_rate'), component),
                    'omega`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component));

                }"
           );

           if ( component < ns_components) {
                model.generic.AddGlobal ( bs_rel, _aux[component-1], terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), component ));
                parameters.DeclareGlobalWithRanges (_aux[component-1], 0.5, 0, 1);
           }
           if ( s_component < syn_components) {
                model.generic.AddGlobal ( bs_rel, _aux_srv[s_component-1], terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), "SRV " + s_component ));
                parameters.DeclareGlobalWithRanges (_aux_srv[s_component-1], 0.5, 0, 1);
           }
           models.codon.generic.DefineQMatrix(bs_rel, namespace);
           rate_matrices [key] = bs_rel[utility.getGlobalValue("terms.model.rate_matrix")];
           (bs_rel [^'terms.mixture.mixture_components'])[key] = parameters.AppendMultiplicativeTerm (_wts [component-1], _wts_srv[s_component- 1]);
        }
    }

    __rp    = parameters.ConstrainMeanOfSet (_alphas, _wts_srv, 1, namespace);

    parameters.DeclareGlobal (__rp[^'terms.global'], null);
    parameters.helper.copy_definitions (bs_rel[^'terms.parameters'], __rp);


    bs_rel[utility.getGlobalValue("terms.model.rate_matrix")] = rate_matrices;
    parameters.SetConstraint(((bs_rel[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.global")])[terms.nucleotideRate("A", "G")], "1", "");

    return bs_rel;
}


/**
 * @name models.codon.BS_REL.set_branch_length
 * @param {Model} model
 * @param {AssociativeList|Number} value
 * @param {String} parameter
 * @returns {Number} 0
 */

lfunction models.codon.BS_REL.set_branch_length(model, value, parameter) {
    if (Type (value) == "Number") {
        return  models.generic.SetBranchLength (model, value*3, parameter);
    } else {
        if (Type (value) == "AssociativeList") {
            //vcopy = value;
            //vcopy[^"terms.branch_length"] = vcopy[^"terms.branch_length"] * 3;
            return  models.generic.SetBranchLength (model, value, parameter);

        }
    }
    return 0;

}

/**
 * @name models.codon.BS_REL.post_definition
 * @param {Dict} model
*/

function models.codon.BS_REL.post_definition(model) {
    model [terms.model.branch_length_string] = model.BranchLengthExpression (model);
    model [terms.model.branch_length_string_conditional] = {};
    return model;
}

/**
 * @name models.codon.BS_REL.get_branch_length
 * @param {Model} model
 * @param {AssociativeList|Number} value
 * @param {String} parameter
 * @returns {Number} 0
 */

function models.codon.BS_REL.get_branch_length(model, tree, node) {
	parameters.SetLocalModelParameters (model, tree, node);
	parameters.SetCategoryVariables   (model);
    bl = utility.GetEnvVariable ("BRANCH_LENGTH_STENCIL");
    if (Type (bl) == "Matrix") {
        if (utility.Has (model [terms.model.branch_length_string_conditional], bl, "String") == FALSE) {
            (model [terms.model.branch_length_string_conditional])[bl] = model.BranchLengthExpression (model);
        }
        bl = Eval ((model [utility.getGlobalValue("terms.model.branch_length_string_conditional")])[bl]);
    } else {
      bl = Eval (model [utility.getGlobalValue("terms.model.branch_length_string")]);
    }
	return bl;
}
