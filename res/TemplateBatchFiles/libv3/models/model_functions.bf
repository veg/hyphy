LoadFunctionLibrary ("DNA.bf");
LoadFunctionLibrary ("parameters.bf");
LoadFunctionLibrary ("frequencies.bf");
LoadFunctionLibrary ("../UtilityFunctions.bf");
LoadFunctionLibrary ("../convenience/regexp.bf");

/** @module model */
/**
 * @name model.GetParameters_RegExp
 * @param model {String} - model ID
 * @param re {String} - regular expression
 * @return a dictionary of global model parameters that match a regexp
 */
lfunction model.GetParameters_RegExp(model, re) {

    names = utility.Filter (utility.Keys ((model[utility.getGlobalValue ("terms.parameters")])[ utility.getGlobalValue("terms.global")]),
                            "_parameter_description_",
                            "None != regexp.Find (_parameter_description_, `&re`)");

    result = {};
    count  = utility.Array1D (names);
    for (k = 0; k < count; k += 1) {
        result [names[k]] = ((model[utility.getGlobalValue ("terms.parameters")])[ utility.getGlobalValue("terms.global")])[names[k]];
    }

    return result;
}

lfunction model.GetLocalParameters_RegExp(model, re) {

    names = utility.Filter (utility.Keys ((model[utility.getGlobalValue ("terms.parameters")])[ utility.getGlobalValue("terms.local")]),
                            "_parameter_description_",
                            "None != regexp.Find (_parameter_description_, `&re`)");

    result = {};
    count  = utility.Array1D (names);
    for (k = 0; k < count; k += 1) {
        result [names[k]] = ((model[utility.getGlobalValue ("terms.parameters")])[ utility.getGlobalValue("terms.local")])[names[k]];
    }

    return result;
}


/**
 * @name model.ApplyModelToTree
 * @param id {String}
 * @param tree
 * @param model_list
 * @param rules
 */
function model.ApplyModelToTree (id, tree, model_list, rules) {


	if (Type (rules) == "AssociativeList") {
	    // this has the form
	    // model id : list of branches to apply the model (as a string COLUMN matrix with branch names,
	    // or a dictionary where keys are the branch names)
	    // OR
	    // DEFAULT : model id


	    if (Abs (rules["DEFAULT"])) {
            ExecuteCommands ("UseModel (" + rules["DEFAULT"] + ");
                              Tree `id` = " + tree["string"] + ";
                              ");
	    } else {
            ExecuteCommands ("UseModel (USE_NO_MODEL);
                              Tree `id` = " + tree["string"] + ";
                              ");

	    }


	    /*

	    debug.log (tree);
	    _t = Eval ("Format (`id`,1,1)");
	    debug.log (_t);

	    */

	    model.ApplyModelToTree.ids = Rows (rules);
	    for (model.ApplyModelToTree.k = 0; model.ApplyModelToTree.k < Abs (rules); model.ApplyModelToTree.k += 1) {
	        model.ApplyModelToTree.name = model.ApplyModelToTree.ids[model.ApplyModelToTree.k];
	        //debug.log (model.ApplyModelToTree.name);
	        if ( model.ApplyModelToTree.name != "DEFAULT") {
                model.ApplyModelToTree.list = rules[model.ApplyModelToTree.name];
                if (Type (model.ApplyModelToTree.list) == "AssociativeList") {
                    model.ApplyModelToTree.list = Rows (model.ApplyModelToTree.list);
                }


                if (Type (model_list) == "AssociativeList") {
                    model.ApplyModelToTree.apply_model = model_list[model.ApplyModelToTree.name];
                } else {
                    model.ApplyModelToTree.apply_model = model.ApplyModelToTree.name;
                }

                for (model.ApplyModelToTree.b = 0; model.ApplyModelToTree.b < Columns (model.ApplyModelToTree.list); model.ApplyModelToTree.b += 1) {
                    //debug.log ( "SetParameter (`id`." + model.ApplyModelToTree.list[model.ApplyModelToTree.b] + ",MODEL," + model.ApplyModelToTree.apply_model + ")");
                    ExecuteCommands ("SetParameter (`id`." + model.ApplyModelToTree.list[model.ApplyModelToTree.b] + ",MODEL," + model.ApplyModelToTree.apply_model + ")");
                }
            }
	    }

	} else {
	    // TO DO: REMOVE HARDCODING


		model.ApplyModelToTree.modelID = model_list[model_list ["INDEXORDER"][0]];
		ExecuteCommands ("UseModel (" + model.ApplyModelToTree.modelID[terms.id] + ");
						  Tree `id` = " + tree["string"] + ";
						  ");

	}
}

/**
 * @name model.ApplyToBranch
 * @param model_id
 * @param tree
 * @param branch
 */
function model.ApplyToBranch (model_id, tree, branch) {
	ExecuteCommands ("SetParameter (`tree`.`branch`,MODEL,`model_id`)");
}


/**
 * @name model.define_from_components
 * @param id - {String} Name of of Model
 * @param q - {Matrix} instantaneous transition matrix
 * @param evf - {Matrix} the equilibrium frequencies
 * @param {Number} canonical - matrix exponential
 * @example
 * q = {{*,t,kappa*t,t}
 *  {t,*,t,t*kappa}
 *  {t*kappa,t,*,t}
 *  {t,t*kappa,t,*}};
 * evf =  {{0.4}{0.3}{0.2}{0.1}};
 * model.define_from_components(q,evf,1);
 * @returns nothing, sets a global {Model}
 */
function model.define_from_components (id,q,efv,canonical) {
	ExecuteCommands ("Model `id` = (" + q + "," + efv + "," + canonical + ")");
}

/**
 * @name model.generic.AddLocal
 * @param model_spec
 * @param id
 * @param tag
 */
function model.generic.AddLocal (model_spec, id, tag) {
    ((model_spec[terms.parameters])[terms.local])[tag] = id;
}

/**
 * @name model.generic.AddGlobal
 * @param model_spec
 * @param id
 * @param tag
 */
function model.generic.AddGlobal (model_spec, id, tag) {
    ((model_spec[terms.parameters])[terms.global])[tag] = id;
}


/**
 * @name model.generic.GetLocalParameter
 * @param model_spec
 * @param tag
 */
lfunction model.generic.GetLocalParameter (model_spec, tag) {
   return model.generic.get_a_parameter (model_spec, tag, utility.getGlobalValue("terms.local"));
}

/**
 * @name model.generic.GetGlobalParameter
 * @param model_spec
 * @param tag
 */
lfunction model.generic.GetGlobalParameter (model_spec, tag) {
   return model.generic.get_a_parameter (model_spec, tag, utility.getGlobalValue("terms.global"));
}

/**
 * @name model.generic.get_a_parameter
 * @param model_spec
 * @param tag
 * @param type
 */
lfunction model.generic.get_a_parameter (model_spec, tag, type) {
    v = ((model_spec[utility.getGlobalValue("terms.parameters")])[type])[tag];
    if (Type (v) == "String") {
        return v;
    }
    return None;
}

/**
 * @name  model.generic.get_rate_variation
 * @param model_spec
 * @param tag
 * @param type
 */

lfunction model.generic.get_rate_variation (model_spec) {

	__rate_variation = model_spec[utility.getGlobalValue("terms.model.rate_variation")];

	if (Type (__rate_variation) == "AssociativeList") {
		assert (utility.Has (__rate_variation, utility.getGlobalValue("terms.rate_variation.distribution"), "String"), "Missing required key `distribution` in the rate variation component definition");
		assert (utility.Has (__rate_variation, utility.getGlobalValue("terms.rate_variation.rate_modifier"), "String"), "Missing required key `rate_modifier' in the rate variation component definition");
		assert (utility.IsFunction (__rate_variation[utility.getGlobalValue("terms.rate_variation.distribution")]), "'" + parameters.Quote (__rate_variation[utility.getGlobalValue("terms.rate_variation.distribution")]) + " must be a function ID");
		assert (utility.IsFunction (__rate_variation[utility.getGlobalValue("terms.rate_variation.rate_modifier")]), "'" + parameters.Quote (__rate_variation[utility.getGlobalValue("terms.rate_variation.rate_modifier")]) + " must be a function ID");
		return __rate_variation;
	}


	return None;
}



/**
 * @name model.generic.DefineModel
 * @param model_spec
 * @param id
 * @param arguments
 * @param data_filter
 * @param estimator_type
 */
function model.generic.DefineModel (model_spec, id, arguments, data_filter, estimator_type) {


    // Basic model definition
	model.generic.DefineModel.model = utility.CallFunction (model_spec, arguments);



	// Add data filter information to model description
	if ( None != data_filter) {
	    models.generic.AttachFilter (model.generic.DefineModel.model, data_filter);
	}

    // Set Q field
	model.generic.DefineModel.model = Call (model.generic.DefineModel.model [terms.model.defineQ], model.generic.DefineModel.model, id);


    // Define type of frequency estimator


	if (None != estimator_type) {
		model.generic.DefineModel.model [terms.model.frequency_estimator] = estimator_type;
	}

    // Set EFV field
	model.generic.DefineModel.model = Call (model.generic.DefineModel.model [terms.model.frequency_estimator], model.generic.DefineModel.model,
													    id,
													    data_filter);


    // Define id's for frequencies, Q, and id
	model.generic.DefineModel.model [terms.model.matrix_id] = "`id`_" + terms.model.rate_matrix;
	model.generic.DefineModel.model [terms.model.efv_id] = "`id`_" + terms.model.efv_matrix;
	model.generic.DefineModel.model [terms.id] = id;


	parameters.StringMatrixToFormulas (model.generic.DefineModel.model [terms.model.matrix_id],model.generic.DefineModel.model[terms.model.rate_matrix]);

	if (Type ((model.generic.DefineModel.model[terms.efv_estimate])[0]) == "String") {
	    parameters.StringMatrixToFormulas (model.generic.DefineModel.model [terms.model.efv_id],model.generic.DefineModel.model[terms.efv_estimate]);
	} else {
	    utility.SetEnvVariable (model.generic.DefineModel.model [terms.model.efv_id], model.generic.DefineModel.model[terms.efv_estimate]);
	}

	model.define_from_components (id, 	model.generic.DefineModel.model [terms.model.matrix_id], model.generic.DefineModel.model [terms.model.efv_id], model.generic.DefineModel.model [terms.model.canonical]);

    if (Type (model.generic.DefineModel.model[terms.model.post_definition]) == "String") {
       Call (model.generic.DefineModel.model[terms.model.post_definition], model.generic.DefineModel.model);
    }

	return model.generic.DefineModel.model;
}


/**
 * @name model.generic.DefineMixtureModel
 * @param model_spec
 * @param id
 * @param arguments
 * @param data_filter
 * @param estimator_type
 */
function model.generic.DefineMixtureModel (model_spec, id, arguments, data_filter, estimator_type) {

	model.generic.DefineModel.model = utility.CallFunction (model_spec, arguments);
	if (None != estimator_type) {
	    models.generic.AttachFilter (model.generic.DefineModel.model, data_filter);
	}

    // for mixture models this will define the mixture components as well
	model.generic.DefineModel.model = Call (model.generic.DefineModel.model [terms.model.defineQ], model.generic.DefineModel.model, id);


	if (estimator_type != None) {
		model.generic.DefineModel.model [terms.model.frequency_estimator] = estimator_type;
	}

    // this sets the EFV field
 	Call (model.generic.DefineModel.model [terms.model.frequency_estimator], model.generic.DefineModel.model,
													    id,
													    data_filter);

	model.generic.DefineModel.model [terms.model.matrix_id] = {};

	model.generic.mixture_expr = {};

    utility.ForEachPair (model.generic.DefineModel.model[terms.mixture.mixture_components], "model.generic.key", "model.generic.value",
        '
            (model.generic.DefineModel.model [terms.model.matrix_id])[model.generic.key] = "`id`_" + terms.model.rate_matrix + "_"+ model.generic.key;
            model.generic.mixture_expr + ("Exp(" + (model.generic.DefineModel.model [terms.model.matrix_id])[model.generic.key] + ")*(" +  model.generic.value + ")");
	        parameters.StringMatrixToFormulas ((model.generic.DefineModel.model [terms.model.matrix_id])[model.generic.key],(model.generic.DefineModel.model[terms.model.rate_matrix])[model.generic.key]);
        '
    );

    model.generic.DefineModel.model [terms.mixture]= Join ("+",model.generic.mixture_expr);

	model.generic.DefineModel.model [terms.model.efv_id] = "`id`_" + terms.model.efv_matrix;
	model.generic.DefineModel.model [terms.id] = id;

	utility.SetEnvVariable (model.generic.DefineModel.model [terms.model.efv_id], model.generic.DefineModel.model[terms.efv_estimate]);

	model.define_from_components (id, 	parameters.Quote (model.generic.DefineModel.model [terms.mixture]), model.generic.DefineModel.model [terms.model.efv_id], model.generic.DefineModel.model [terms.model.canonical]);


    if (Type (model.generic.DefineModel.model[terms.model.post_definition]) == "String") {
        Call (model.generic.DefineModel.model[terms.model.post_definition], model.generic.DefineModel.model);
    }

	return model.generic.DefineModel.model;
}

/**
 * @name model.generic.GetLocalParameter
 * @param {Model} model
 * @returns {String}
 */
function models.generic.post.definition  (model) {
    if (Type (model [terms.id]) == "String") {
        ExecuteCommands ("GetString (models.generic.post.definition.bl,"+model [terms.id]+",-1)");
        model [terms.model.branch_length_string] = models.generic.post.definition.bl;
    }
    return model;
}

/**
 * @name models.generic.ConstrainBranchLength
 * @param {Model} model
 * @param {AssociativeList} or {Number} value
 * @param {String} parameter
 * @returns the number of constraints generated (0 or 1)
 */
function models.generic.ConstrainBranchLength (model, value, parameter) {
    if (Type (value) == "Number") {
        if (Abs((model[terms.parameters])[terms.local]) == 1) {
            if (Type (model [terms.model.branch_length_string]) == "String") {

                models.generic.ConstrainBranchLength.expression = model [terms.model.branch_length_string];
                models.generic.ConstrainBranchLength.bl = (Columns ((model[terms.parameters])[terms.local]))[0];
                models.generic.ConstrainBranchLength.bl.p = parameter + "." + models.generic.ConstrainBranchLength.bl;
                models.generic.ConstrainBranchLength.substitution = {models.generic.ConstrainBranchLength.bl : 1};
                utility.Extend (models.generic.ConstrainBranchLength.substitution, parameters.SetCategoryVariables (model));
                models.generic.ConstrainBranchLength.expression = "(" + Simplify (models.generic.ConstrainBranchLength.expression, models.generic.ConstrainBranchLength.substitution) + ")";
                Eval (models.generic.ConstrainBranchLength.bl.p + ":=" + value + "/" + models.generic.ConstrainBranchLength.expression);
                messages.log ("models.generic.ConstrainBranchLength: " + models.generic.ConstrainBranchLength.bl.p + ":=" + value + "/" + models.generic.ConstrainBranchLength.expression);
            } else {
                messages.log ("models.generic.ConstrainBranchLength: model branch length expression not available");
            }

        }
        else {
            messages.log ("models.generic.ConstrainBranchLength: not exactly one local model parameter");
        }
    } else {
         messages.log ("models.generic.ConstrainBranchLength: unsupported value type " + Type (value) + "\n" + value);
    }
    return 0;
}

/**
 * @name models.generic.SetBranchLength
 * @param {Model} model
 * @param {AssociativeList} or {Number} value
 * @param {String} parameter
 * @returns the number of constraints generated (0 or 1)
 */
function models.generic.SetBranchLength (model, value, parameter) {


     if (Abs((model[terms.parameters])[terms.local]) >= 1) {
        if (Type (model [terms.model.branch_length_string]) == "String") {
            models.generic.SetBranchLength.expression = model [terms.model.branch_length_string];

            if (Abs((model[terms.parameters])[terms.local]) > 1) {
                models.generic.SetBranchLength.bl = Call (model [terms.model.time], model[terms.model.type]);
                models.generic.SetBranchLength.bl.p = parameter + "." + models.generic.SetBranchLength.bl;
                if (parameters.IsIndependent (models.generic.SetBranchLength.bl.p) == FALSE) {
                     models.generic.SetBranchLength.bl = utility.First (utility.UniqueValues ((model[terms.parameters])[terms.local]), "_name_",
                            'parameters.IsIndependent (parameter + "." + _name_)');

                     if (None == models.generic.SetBranchLength.bl) {
                        messages.log ("models.generic.SetBranchLength: no independent model parameters");
                        return 0;
                     }
                     models.generic.SetBranchLength.bl.p = parameter + "." + models.generic.SetBranchLength.bl;
                }

                models.generic.SetBranchLength.substitution = {};
                utility.ForEach (utility.UniqueValues ((model[terms.parameters])[terms.local]), "_name_",
                            'if (_name_ != models.generic.SetBranchLength.bl) {
                                models.generic.SetBranchLength.substitution [_name_] = Eval (parameter + "." + _name_);
                             }');

                models.generic.SetBranchLength.expression = Simplify (models.generic.SetBranchLength.expression, models.generic.SetBranchLength.substitution);
            } else {
                models.generic.SetBranchLength.bl = (Columns ((model[terms.parameters])[terms.local]))[0];
                models.generic.SetBranchLength.bl.p = parameter + "." + models.generic.SetBranchLength.bl;
            }

             if (parameters.IsIndependent (models.generic.SetBranchLength.bl.p)) {

                if (Type (value) == "AssociativeList") {
                  	if (Abs (models.generic.SetBranchLength.expression)) {
                    	ExecuteCommands ("FindRoot (models.generic.SetBranchLength.t,(" + models.generic.SetBranchLength.expression + ")-" + value[terms.branch_length] + "," + models.generic.SetBranchLength.bl + ",0,10000)");
                    	Eval (parameter + "." + models.generic.SetBranchLength.bl + ":=(" + value[terms.model.branch_length_scaler] + ")*" + models.generic.SetBranchLength.t);
					    messages.log ("models.generic.SetBranchLength: " + parameter + "." + models.generic.SetBranchLength.bl + ":=(" + value[terms.model.branch_length_scaler] + ")*" + models.generic.SetBranchLength.t);

					} else {
                    	Eval (parameter + "." + models.generic.SetBranchLength.bl + ":=(" + value[terms.model.branch_length_scaler] + ")*" + value[terms.branch_length]);
                    	messages.log ("models.generic.SetBranchLength: " + parameter + "." + models.generic.SetBranchLength.bl + ":=(" + value[terms.model.branch_length_scaler] + ")*" + models.generic.SetBranchLength.t);
					}

                    return 1;

                } else {
                     ExecuteCommands ("FindRoot (models.generic.SetBranchLength.t,(" + models.generic.SetBranchLength.expression + ")-" + value + "," + models.generic.SetBranchLength.bl + ",0,10000)");
                     Eval (parameter + "." + models.generic.SetBranchLength.bl + "=" + models.generic.SetBranchLength.t);
                     messages.log ("models.generic.SetBranchLength: " + parameter + "." + models.generic.SetBranchLength.bl + "=" + models.generic.SetBranchLength.t);
              }
            } else {
                messages.log (models.generic.SetBranchLength.bl.p + " was already constrained in models.generic.SetBranchLength");
            }
        } else {
	        messages.log ("models.generic.SetBranchLength: missing branch-length-string");
        }
    } else {
        messages.log ("models.generic.SetBranchLength: no local model parameters");
    }
    return 0;
}



/**
 * @name models.generic.AttachFilter
 * @param {Model} model
 * @param {DataSetFilter} filter
 * @returns 0
 */
lfunction models.generic.AttachFilter (model, filter) {

    if (Type (filter) != "String") {
        utility.ForEach (filter, "_filter_", "models.generic.AttachFilter (`&model`, _filter_)");
        model[utility.getGlobalValue("terms.model.data")] = filter;
        return model;
    }

	GetDataInfo (givenAlphabet, ^filter, "CHARACTERS");
	alphabet = model [utility.getGlobalValue("terms.alphabet")];

	assert (Columns (alphabet) == Columns (givenAlphabet) && model.MatchAlphabets (givenAlphabet, alphabet), "The declared model alphabet '" +alphabet + "' does not match the `filter` filter: '" + givenAlphabet + "'");

	model [utility.getGlobalValue("terms.alphabet")] = givenAlphabet;
	model [utility.getGlobalValue("terms.model.data")] = filter;
	return model;
}

/**
 * @name model.Dimension
 * @param {Model} model
 * @returns {Matrix} dimensions of model
 */
function model.Dimension (model) {
    if (Type (model[terms.alphabet]) == "Matrix") {
        return Columns (model[terms.alphabet]);
    }
    return None;
}

/**
 * @name model.parameters.Local
 * @param {Model} model
 * @returns {Matrix} local parameters
 */
function model.parameters.Local (model) {
    return (model[terms.parameters])[terms.local];
}


/**
 * @name model.parameters.Global
 * @param {Model} model
 * @returns {Matrix} global parameters
 */
function model.parameters.Global (model) {
    return (model[terms.parameters])[terms.global];
}

/**
 * @name model.MatchAlphabets
 * @param {Matrix} a1 - first alphabet to compare
 * @param {Matrix} a2 - second alphabet to compare
 * @returns {Number} 1 if they are equal, 0 if they are not
 */
lfunction model.MatchAlphabets (a1, a2) {


	_validStates = {};
	for (_k = 0; _k < Columns (a1); _k += 1) {
		_validStates [a1[_k]] = _k;
	}
	for (_k = 0; _k < Columns (a2); _k += 1) {
		if (_validStates [a2[_k]] != _k) {
			return 0;
		}
	}
	return 1;

}

/**
 * @name models.BindGlobalParameters
 * compute the algebraic expression for the branch length
 * @param {Dict} models - list of model objects
 * @param {RegEx} filter - only apply to parameters matching this expression
 * @return {Dict} the set of constraints applied
 */

lfunction models.BindGlobalParameters (models, filter) {


    if (Type (models) == "AssociativeList" && utility.Array1D (models) > 1) {

        reference_set = (((models[0])[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.global")]);
        candidate_set = utility.UniqueValues(utility.Filter (utility.Keys (reference_set), "_key_",
            "regexp.Find (_key_,`&filter`)"
        ));

        constraints_set = {};

        for (k = 1; k < Abs (models); k+=1) {
            parameter_set = (((models[k])[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.global")]);
            if (Type (parameter_set) == "AssociativeList") {
                utility.ForEach (candidate_set, "_p_",
                    "if (`&parameter_set` / _p_) {
                        if (parameters.IsIndependent (`&parameter_set`[_p_])) {
                            parameters.SetConstraint (`&parameter_set`[_p_], `&reference_set`[_p_], '');
                            `&constraints_set` [`&parameter_set`[_p_]] = `&reference_set`[_p_];
                        }
                    }"
                );
            }
        }

        return constraints_set;

    }
    return None;
}

/**
 * @name model.BranchLengthExpression
 * compute the algebraic expression for the branch length
 * @param {Dict} model - first alphabet to compare
 * @sa model.BranchLengthExpressionFromMatrix
 * @returns {String} the branch length expression string
 */

lfunction model.BranchLengthExpression (model) {


	if (Type (model[utility.getGlobalValue("terms.model.rate_matrix")]) == "Matrix") {
		expr = model.BranchLengthExpressionFromMatrix (model[utility.getGlobalValue("terms.model.rate_matrix")], model[utility.getGlobalValue("terms.efv_estimate")], model[utility.getGlobalValue("terms.model.canonical")]);
	} else {
		components = {};
		matrix_count = Abs (model[utility.getGlobalValue("terms.model.rate_matrix")]);
		//keys = utility.Keys (model[utility.getGlobalValue("terms.model.rate_matrix"));
	    keys = utility.Keys (model[utility.getGlobalValue("terms.model.rate_matrix")]);



		for (i = 0; i <  matrix_count; i+=1) {
			expr = model.BranchLengthExpressionFromMatrix ((model[utility.getGlobalValue("terms.model.rate_matrix")])[keys[i]], model[utility.getGlobalValue("terms.efv_estimate")], model[utility.getGlobalValue("terms.model.canonical")]);
			components + ( "(" + expr + ")*(" + (model[utility.getGlobalValue("terms.mixture.mixture_components")])[keys[i]] + ")");
		}
		expr = Join ("+", components);
	}

	return "(" + expr + ")/3";
}

lfunction model.BranchLengthExpressionFromMatrix (q,freqs,is_canonical) {

    by_expr = {};
    dim = Rows (q);
    can_alias = Type (freqs[0]) == "Number";
    stencil = utility.GetEnvVariable ("BRANCH_LENGTH_STENCIL");

    if (Type (stencil) == "Matrix") {
        if (Rows (stencil) != dim) {
          stencil = None;
        }
    }

    if (Type (stencil) != "Matrix") {
        stencil = {dim,dim}["1"];
    }

    if (can_alias) {
        for (i = 0; i < dim; i+=1) {
            for (j = 0; j < dim; j+=1) {
                if (i != j && stencil[i][j]) {
                    expr = q[i][j];
                    if (Abs (expr)) {
                        if (is_canonical == TRUE) {
                            by_expr[expr] += freqs[i]*freqs[j];
                        } else {
                            by_expr[expr] += freqs[i];
                        }
                    }
                }
            }
        }
        expr = {};

        for (_expr_, _wt_; in; by_expr) {
          expr + ("(" + _expr_ + ")*" + _wt_);
        }
          
    } else {
        expr = {};
        for (i = 0; i < dim; i+=1) {
            for (j = 0; j < dim; j+=1) {
                if (i != j && stencil[i][j]) {
                    rate = q[i][j];
                    if (Abs (rate)) {
                        if (is_canonical) {
                            expr + "(`rate`) * (`freqs[i]`)*(`freqs[j]`)";
                        } else {
                            expr +  "(`rate`) * (`freqs[i]`)";
                        }
                    }
                }
            }
        }

    }

    expr = Join ("+", expr);
    return expr;
}
