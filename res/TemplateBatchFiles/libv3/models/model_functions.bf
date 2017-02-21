LoadFunctionLibrary ("DNA.bf");
LoadFunctionLibrary ("parameters.bf");
LoadFunctionLibrary ("frequencies.bf");
LoadFunctionLibrary ("../UtilityFunctions.bf");

/** @module model */

/**
 * @name model.ApplyModelToTree
 * @param id
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

	    model.ApplyModelToTree.ids = Rows (rules);
	    for (model.ApplyModelToTree.k = 0; model.ApplyModelToTree.k < Abs (rules); model.ApplyModelToTree.k += 1) {
	        model.ApplyModelToTree.name = model.ApplyModelToTree.ids[model.ApplyModelToTree.k];
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
                    //fprintf (stdout, "SetParameter (`id`." + model.ApplyModelToTree.list[model.ApplyModelToTree.b] + ",MODEL," + model.ApplyModelToTree.apply_model + ")", "\n");
                    ExecuteCommands ("SetParameter (`id`." + model.ApplyModelToTree.list[model.ApplyModelToTree.b] + ",MODEL," + model.ApplyModelToTree.apply_model + ")");
                }
            }
	    }

	} else {
		model.ApplyModelToTree.modelID = model_list[model_list ["INDEXORDER"][0]];
		ExecuteCommands ("UseModel (" + model.ApplyModelToTree.modelID["id"] + ");
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
 * @name model.generic.AddGlobal
 * @param model_spec
 * @param id
 * @param tag
 */
function model.generic.AddLocal (model_spec, id, tag) {
    ((model_spec["parameters"])[terms.local])[tag] = id;
}

/**
 * @name model.generic.AddGlobal
 * @param model_spec
 * @param id
 * @param tag
 */
function model.generic.AddGlobal (model_spec, id, tag) {
    ((model_spec["parameters"])[terms.global])[tag] = id;
}


/**
 * @name model.generic.GetLocalParameter
 * @param model_spec
 * @param tag
 */
lfunction model.generic.GetLocalParameter (model_spec, tag) {
   return model.generic.get_a_parameter (model_spec, tag, ^"terms.local");
}

/**
 * @name model.generic.GetGlobalParameter
 * @param model_spec
 * @param tag
 */
lfunction model.generic.GetGlobalParameter (model_spec, tag) {
   return model.generic.get_a_parameter (model_spec, tag, ^"terms.global");
}

/**
 * @name model.generic.get_a_parameter
 * @param model_spec
 * @param tag
 * @param type
 */
lfunction model.generic.get_a_parameter (model_spec, tag, type) {
    v = ((model_spec["parameters"])[type])[tag];
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
	
	__rate_variation = model_spec[^"terms.rate_variation"];

	if (Type (__rate_variation) == "AssociativeList") {
		assert (utility.Has (__rate_variation, "distribution", "String"), "Missing required key `distribution' in the rate variation component definition");
		assert (utility.Has (__rate_variation, "rate_modifier", "String"), "Missing required key `rate_modifier' in the rate variation component definition");
		assert (utility.IsFunction (__rate_variation["distribution"]), "'" + parameters.Quote (__rate_variation["distribution"]) + " must be a function ID");
		assert (utility.IsFunction (__rate_variation["rate_modifier"]), "'" + parameters.Quote (__rate_variation["rate_modifier"]) + " must be a function ID");	
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

	model.generic.DefineModel.model = utility.CallFunction (model_spec, arguments);
	models.generic.AttachFilter (model.generic.DefineModel.model, data_filter);

	model.generic.DefineModel.model = Call (model.generic.DefineModel.model ["defineQ"], model.generic.DefineModel.model, id);

	model.generic.DefineModel.model ["matrix-id"] = "`id`_" + terms.rate_matrix;
	model.generic.DefineModel.model ["efv-id"] = "`id`_" + terms.efv_matrix;
	model.generic.DefineModel.model ["id"] = id;

	if (estimator_type != None) {
		model.generic.DefineModel.model ["frequency-estimator"] = estimator_type;
	}

	Call (model.generic.DefineModel.model ["frequency-estimator"], model.generic.DefineModel.model,
													    id,
													    data_filter); // this sets the EFV field


	parameters.StringMatrixToFormulas (model.generic.DefineModel.model ["matrix-id"],model.generic.DefineModel.model[terms.rate_matrix]);
	utility.SetEnvVariable (model.generic.DefineModel.model ["efv-id"], model.generic.DefineModel.model[terms.efv_estimate]);

	model.define_from_components (id, 	model.generic.DefineModel.model ["matrix-id"], model.generic.DefineModel.model ["efv-id"], model.generic.DefineModel.model ["canonical"]);

    if (Type (model.generic.DefineModel.model["post-definition"]) == "String") {
        Call (model.generic.DefineModel.model["post-definition"], model.generic.DefineModel.model);
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
	models.generic.AttachFilter (model.generic.DefineModel.model, data_filter);

	model.generic.DefineModel.model = Call (model.generic.DefineModel.model ["defineQ"], model.generic.DefineModel.model, id);
	    // for mixture models this will define the mixture components as well

	if (estimator_type != None) {
		model.generic.DefineModel.model ["frequency-estimator"] = estimator_type;
	}


 	Call (model.generic.DefineModel.model ["frequency-estimator"], model.generic.DefineModel.model,
													    id,
													    data_filter); // this sets the EFV field

	model.generic.DefineModel.model ["matrix-id"] = {};

	model.generic.mixture_expr = {};
	
    utility.ForEachPair (model.generic.DefineModel.model[terms.mixture_components], "model.generic.key", "model.generic.value",
        '
            (model.generic.DefineModel.model ["matrix-id"])[model.generic.key] = "`id`_" + terms.rate_matrix + "_"+ model.generic.key;
            model.generic.mixture_expr + ("Exp(" + (model.generic.DefineModel.model ["matrix-id"])[model.generic.key] + ")*(" +  model.generic.value + ")");
	        parameters.StringMatrixToFormulas ((model.generic.DefineModel.model ["matrix-id"])[model.generic.key],(model.generic.DefineModel.model[terms.rate_matrix])[model.generic.key]);
        '
    );

    model.generic.DefineModel.model [terms.mixture]= Join ("+",model.generic.mixture_expr);

	model.generic.DefineModel.model ["efv-id"] = "`id`_" + terms.efv_matrix;
	model.generic.DefineModel.model ["id"] = id;

	utility.SetEnvVariable (model.generic.DefineModel.model ["efv-id"], model.generic.DefineModel.model[terms.efv_estimate]);

	model.define_from_components (id, 	parameters.Quote (model.generic.DefineModel.model [terms.mixture]), model.generic.DefineModel.model ["efv-id"], model.generic.DefineModel.model ["canonical"]);


    if (Type (model.generic.DefineModel.model["post-definition"]) == "String") {
        Call (model.generic.DefineModel.model["post-definition"], model.generic.DefineModel.model);
    }

	return model.generic.DefineModel.model;
}

/**
 * @name model.generic.GetLocalParameter
 * @param {Model} model
 * @returns {String}
 */
function models.generic.post.definition  (model) {
    if (Type (model ["id"]) == "String") {
        ExecuteCommands ("GetString (models.generic.post.definition.bl,"+model ["id"]+",-1)");
        model ["branch-length-string"] = models.generic.post.definition.bl;
    } 
    return model;
}

/**
 * @name models.generic.SetBranchLength
 * @param {Model} model
 * @param {AssociativeList} or {Number} value
 * @param {String} parameter
 * @returns the number of constraints generated (0 or 1)
 */
function models.generic.SetBranchLength (model, value, parameter) {

    if (Abs((model["parameters"])["local"]) == 1) {
        if (Type (model ["branch-length-string"]) == "String") {
            models.generic.SetBranchLength.bl = (Columns ((model["parameters"])["local"]))[0];
            models.generic.SetBranchLength.bl.p = parameter + "." + models.generic.SetBranchLength.bl;
             if (parameters.IsIndependent (models.generic.SetBranchLength.bl.p)) {
                if (Type (value) == "AssociativeList") {
                	if (Abs (model ["branch-length-string"])) {
                    	ExecuteCommands ("FindRoot (models.generic.SetBranchLength.t,(" + model ["branch-length-string"] + ")-" + value[terms.branch_length] + "," + models.generic.SetBranchLength.bl + ",0,10000)");
                    	Eval (parameter + "." + models.generic.SetBranchLength.bl + ":=(" + value[terms.branch_length_scaler] + ")*" + models.generic.SetBranchLength.t);
					} else {
                    	Eval (parameter + "." + models.generic.SetBranchLength.bl + ":=(" + value[terms.branch_length_scaler] + ")*" + value[terms.branch_length]);
					}

                    return 1;

                } else {
                    ExecuteCommands ("FindRoot (models.generic.SetBranchLength.t,(" + model ["branch-length-string"] + ")-" + value + "," + models.generic.SetBranchLength.bl + ",0,10000)");
                    Eval (parameter + "." + models.generic.SetBranchLength.bl + "=" + models.generic.SetBranchLength.t);
               }
            } else {
                warning.log (models.generic.SetBranchLength.bl.p + " was already constrained in models.generic.SetBranchLength");
            }
        }
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
        model["data"] = filter;
        return model;
    }

	GetDataInfo (givenAlphabet, ^filter, "CHARACTERS");
	alphabet = model ["alphabet"];

	assert (Columns (alphabet) == Columns (givenAlphabet) && model.MatchAlphabets (givenAlphabet, alphabet), "The declared model alphabet '" +alphabet + "' does not match the `filter` filter: '" + givenAlphabet + "'");

	model ["alphabet"] = givenAlphabet;
	model ["data"] = filter;
	return model;
}

/**
 * @name model.Dimension
 * @param {Model} model
 * @returns {Matrix} dimensions of model
 */
function model.Dimension (model) {
    if (Type (model["alphabet"]) == "Matrix") {
        return Columns (model["alphabet"]);
    }
    return None;
}

/**
 * @name model.parameters.Local
 * @param {Model} model
 * @returns {Matrix} local parameters
 */
function model.parameters.Local (model) {
    return (model["parameters"])["local"];
}


/**
 * @name model.parameters.Global
 * @param {Model} model
 * @returns {Matrix} global parameters
 */
function model.parameters.Global (model) {
    return (model["parameters"])["global"];
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
 * @name model.BranchLengthExpression
 * compute the algebraic expression for the branch length
 * @param {Dict} model - first alphabet to compare
 * @sa model.BranchLengthExpressionFromMatrix
 * @returns {String} the branch length expression string
 */
 
lfunction model.BranchLengthExpression (model) {
	if (Type (model[^'terms.rate_matrix']) == "Matrix") {
		expr = model.BranchLengthExpressionFromMatrix (model[^'terms.rate_matrix'], model[^'terms.efv_estimate'], model['canonical']);
	} else {
		components = {};
		matrix_count = Abs (model[^'terms.rate_matrix']);
		keys = utility.Keys (model[^'terms.rate_matrix']);
		for (i = 0; i <  matrix_count; i+=1) {
			expr = model.BranchLengthExpressionFromMatrix ((model[^'terms.rate_matrix'])[keys[i]], model[^'terms.efv_estimate'], model['canonical']);
			components + ( "(" + expr + ")*(" + (model[^'terms.mixture_components'])[keys[i]] + ")");
		}
		expr = Join ("+", components);
	}
	
	return "(" + expr + ")/3";
}

lfunction model.BranchLengthExpressionFromMatrix (q,freqs,is_canonical) {
	
	by_expr = {};
	dim = Rows (q);
	can_alias = Type (freqs[0]) == "Number";
	
	if (can_alias) {	
		for (i = 0; i < dim; i+=1) {
			for (j = 0; j < dim; j+=1) {
				if (i != j) {
					expr = q[i][j];
					if (Abs (expr)) {
						if (is_canonical) {
							by_expr[expr] += freqs[i]*freqs[j];
						} else {
							by_expr[expr] += freqs[i];
						}
					}
				}
			}
		}
		expr = {};
		utility.ForEachPair (by_expr, "_expr_", "_wt_", 
		'
			`&expr` + ("(" + _expr_ + ")*" + _wt_)
		');
	} else {
		expr = {};
		for (i = 0; i < dim; i+=1) {
			for (j = 0; j < dim; j+=1) {
				if (i != j) {
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
