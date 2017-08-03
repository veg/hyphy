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
        utility.getGlobalValue("terms.alphabet"): codons[utility.getGlobalValue("terms.sense_codons")],
        utility.getGlobalValue("terms.bases"): utility.getGlobalValue("models.DNA.alphabet"),
        utility.getGlobalValue("terms.model.components") : components,
        utility.getGlobalValue("terms.stop_codons"): codons["terms.stop_codons"],
        utility.getGlobalValue("terms.model.type"): type,
        utility.getGlobalValue("terms.translation_table"): codons["terms.translation_table],
        utility.getGlobalValue("terms.description"): "The branch-site mixture of N Muse-Gaut 94 codon-substitution model coupled with the general time reversible (GTR) model of nucleotide substitution",
        utility.getGlobalValue("terms.model.canonical"): "EXPLICIT_FORM_MATRIX_EXPONENTIAL", // is NOT of the r_ij \times \pi_j form
        utility.getGlobalValue("terms.model.reversible"): TRUE,
        utility.getGlobalValue("terms.model.efv_estimate_name"): utility.getGlobalValue("terms.freqs.CF3x4"),
        utility.getGlobalValue("terms.parameters"): {
            utility.getGlobalValue("terms.global"): {},
            utility.getGlobalValue("terms.local"):  {}
        },
        utility.getGlobalValue("terms.model.get_branch_length"): "models.codon.BS_REL.get_branch_length",
        utility.getGlobalValue("terms.model.set_branch_length"): "models.codon.BS_REL.set_branch_length",
        utility.getGlobalValue("terms.model.constrain-branch_length"): "models.codon.BS_REL.constrain_branch_length",
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
 * @name models.codon.BS_REL._DefineQ
 * @param {Dict} mg_rev
 * @param {String} namespace
 * @returns {Dict} updated model
 */

lfunction models.codon.BS_REL_Per_Branch_Mixing._DefineQ(bs_rel, namespace) {

    rate_matrices = {};

    bs_rel [utility.getGlobalValue("terms.model.q_ij")] = &rate_generator;
    bs_rel [utility.getGlobalValue("terms.mixture_components")] = {};

    _aux = parameters.GenerateSequentialNames ("bsrel_mixture_aux", bs_rel[utility.getGlobalValue("terms.model.components")] - 1, "_");
    _wts = parameters.helper.stick_breaking (_aux, None);
    mixture = {};


    for (component = 1; component <= bs_rel[utility.getGlobalValue("terms.model.components")]; component += 1) {
       key = "component_" + component;
       ExecuteCommands ("
        function rate_generator (fromChar, toChar, namespace, model_type, _tt) {
            return models.codon.MG_REV._GenerateRate_generic (fromChar, toChar, namespace, model_type, _tt,
                'alpha', utility.getGlobalValue("terms.synonymous_rate"),
                'beta_`component`', terms.AddCategory (utility.getGlobalValue("terms.nonsynonymous_rate")', component),
                'omega`component`', terms.AddCategory (utility.getGlobalValue("terms.omega_ratio"), component));
        }"
       );
       models.codon.generic.DefineQMatrix(bs_rel, namespace);
       rate_matrices [key] = bs_rel[utility.getGlobalValue("terms.model.rate_matrix")];
       (bs_rel [utility.getGlobalValue("terms.mixture_components")])[key] = _wts [component-1];

       if ( component < bs_rel[utility.getGlobalValue("terms.model.components")]) {
            model.generic.AddLocal ( bs_rel, _aux[component-1], terms.AddCategory (utility.getGlobalValue("terms.mixture_aux_weight"), component ));
            parameters.SetRange (_aux[component-1], utility.getGlobalValue("terms.range_almost_01"));
        }
    	//mixture [key] = "mixture_weight_" + component;
        //model.generic.AddLocal ( bs_rel, mixture [key], terms.AddCategory (^'terms.mixture_weight', component));
        //parameters.SetConstraint (mixture [key], _wts[component-1], "");

    }

    bs_rel[utility.getGlobalValue("terms.model.rate_matrix")] = rate_matrices;
    //bs_rel[^"terms.mixture_components"] = mixture;
    
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
 
function models.codon.BS_REL.set_branch_length(model, value, parameter) {
    assert (FALSE, "models.codon.BS_REL.set_branch_length is not yet implemented");
    return 0;
}

/**
 * @name models.codon.BS_REL.post_definition 
 * @param {Dict} model
*/

function models.codon.BS_REL.post_definition(model) {
    model [terms.model.branch_length_string] = model.BranchLengthExpression (model);
    return model;
}

/**
 * @name models.codon.BS_REL.set_branch_length
 * @param {Model} model
 * @param {AssociativeList|Number} value
 * @param {String} parameter
 * @returns {Number} 0
 */
 
function models.codon.BS_REL.get_branch_length(model, tree, node) {	
	parameters.SetLocalModelParameters (model, tree, node);
	bl = Eval (model [terms.model.branch_length_string]);
	//console.log ("Branch length = " + bl);
	return bl;
}

