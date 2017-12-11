LoadFunctionLibrary("../models/model_functions.bf");
LoadFunctionLibrary("../models/DNA/GTR.bf");
LoadFunctionLibrary("../convenience/regexp.bf");
LoadFunctionLibrary("libv3/all-terms.bf");


/**
 * @name estimators.TakeLFStateSnapshot
 * @param {String} lf_id
 * @returns {Dict} parameter -> {"MLE" : value, "constraint" : string (if present)}
 */

lfunction estimators.TakeLFStateSnapshot(lf_id) {
    snapshot = {};
    GetString (info, ^lf_id,-1);

    utility.ForEach (info[utility.getGlobalValue("terms.parameters.global_independent")], "_name_",
                    '`&snapshot`[_name_] = {terms.fit.MLE : Eval (_name_)};');
    utility.ForEach (info[utility.getGlobalValue("terms.parameters.local_independent")], "_name_",
                    '`&snapshot`[_name_] = {terms.fit.MLE : Eval (_name_)};');
    utility.ForEach (info[utility.getGlobalValue("terms.parameters.global_constrained")], "_name_",
                    '`&snapshot`[_name_] = {terms.fit.MLE : Eval (_name_),
                                            terms.constraint : parameters.GetConstraint (_name_)};');
    utility.ForEach (info[utility.getGlobalValue("terms.parameters.local_constrained")], "_name_",
                    '`&snapshot`[_name_] = {terms.fit.MLE : Eval (_name_),
                                            terms.constraint : parameters.GetConstraint (_name_)};');

    return snapshot;
}

lfunction estimators.RestoreLFStateFromSnapshot(lf_id, snapshot) {
    utility.ForEachPair (snapshot, "_name_", "_info_",
    '
        if (_info_ / terms.constraint) {
            parameters.SetConstraint (_name_, _info_ [terms.constraint], "");

        } else {
            parameters.SetValue (_name_, _info_ [terms.fit.MLE]);
        }
    ');
}

/**
 * @name estimators.GetGlobalMLE
 * @param {Dictionary} results
 * @param {String} tag
 * @returns None
 */
lfunction estimators.GetGlobalMLE(results, tag) {
    estimate = (results[ utility.getGlobalValue("terms.global")])[tag];

    if (Type(estimate) == "AssociativeList") {
        return estimate[ utility.getGlobalValue("terms.fit.MLE")];
    }
    return None;
}

/**
 * @name estimators.GetBranchEstimates
 * @param {Dictionary} results
 * @param {Number} partition_index
 * @param {String} node_name
 * @returns None
 */
lfunction estimators.GetBranchEstimates (results, partition_index, node_name) {

    estimate = ((results[ utility.getGlobalValue("terms.branch_length")])[partition_index])[node_name];

    if (Type(estimate) == "AssociativeList") {
        return estimate;
    }
    return None;
}


/**
 * Extract global scope parameter estimates that match a regular expression
 * @name estimators.GetGlobalMLE_RegExp
 * @param {Dictionary} results
 * @param {String} regular expression to match
 * @returns {Dict} parameter description : value (could be empty)
 */
lfunction estimators.GetGlobalMLE_RegExp(results, re) {

    names = utility.Filter (utility.Keys (results[ utility.getGlobalValue("terms.global")]),
                            "_parameter_description_",
                            "None != regexp.Find (_parameter_description_, `&re`)");

    result = {};
    count  = utility.Array1D (names);
    for (k = 0; k < count; k += 1) {
        result [names[k]] = (results[ utility.getGlobalValue("terms.global")])[names[k]];
    }

    return result;
}

/**
 * @name estimators.copyGlobals2
 * @private
 * @param {String} key2
 * @param {String} value2
 * @returns nothing
 */
function estimators.copyGlobals2(key2, value2) {

    if (Type((estimators.ExtractMLEs.results[terms.global])[key2]) == "AssociativeList") {
        key2 = "[`key`] `key2`";
    }

    (estimators.ExtractMLEs.results[terms.global])[key2] = {
        terms.id: value2,
        terms.fit.MLE: Eval(value2)
    };

    if (parameters.IsIndependent(value2) != TRUE) {
        ((estimators.ExtractMLEs.results[terms.global])[key2])[terms.constraint] = parameters.GetConstraint(value2);
    }
}

/**
 * @name estimators.copyGlobals
 * @private
 * @param {String} key
 * @param {String} value
 * @returns nothing
 */
function estimators.copyGlobals(key, value) {
    ((value[terms.parameters])[terms.global])["estimators.copyGlobals2"][""];
}

/**
 * @name estimators.CopyFrequencies
 * @private
 * @param {String} key
 * @param {Dictionary} value
 * @returns nothing
 */
function estimators.CopyFrequencies(model_name, model_description) {
    (estimators.ExtractMLEs.results[terms.efv_estimate])[model_name] = model_description[terms.efv_estimate];
}



function estimators.SetGlobals2(key, value) {
    __init_value = (initial_values[terms.global])[key];
    if (Type(__init_value) == "AssociativeList") {
        if (__init_value[terms.fix]) {
            estimators.ApplyExistingEstimates.df_correction += parameters.IsIndependent(value);
            ExecuteCommands("`value` := " + __init_value[terms.fit.MLE]);
        } else {
            if (parameters.IsIndependent (value)) {
                //fprintf (stdout, "Setting `value` to " + __init_value[terms.fit.MLE] + "\n", parameters.IsIndependent (value), "\n");
                ExecuteCommands("`value` = " + __init_value[terms.fit.MLE]);
            } else {
                messages.log (value + " was already constrained in estimators.SetGlobals2");
            }
        }
    }
}

/**
 * @name estimators.SetGlobals
 * @param {String} key
 * @param {String} value
 * @returns nothing
 */
function estimators.SetGlobals(key, value) {
    ((value[terms.parameters])[terms.global])["estimators.SetGlobals2"][""];
}

/**
 * @name estimators.SetCategory
 * @param {String} key
 * @param {String} value
 * @returns nothing
 */
function estimators.SetCategory (key, value) {
    ((value[terms.parameters])[terms.category])["estimators.SetCategory2"][""];
}

function estimators.SetCategory2(key, value) {
	^key = 1;
}


/**
 * @name estimators.ExtractBranchInformation.copy_local
 * @param {String} key
 * @param {String} value
 */
function estimators.ExtractBranchInformation.copy_local(key, value) {

    estimators.ExtractBranchInformation.copy_local.var_name = estimators.extractBranchLength.parameter_tag + "." + value;

    estimators.extractBranchLength.result[key] = {
        terms.id: value,
        terms.fit.MLE: Eval(estimators.ExtractBranchInformation.copy_local.var_name)
    };

    if (parameters.IsIndependent(estimators.ExtractBranchInformation.copy_local.var_name) != TRUE) {
        (estimators.extractBranchLength.result[key])[terms.constraint] = parameters.GetConstraint(estimators.ExtractBranchInformation.copy_local.var_name);
    }

}

/**
 * @name estimators.ExtractBranchInformation
 * @param {String} type
 * @param {String} node
 * @param {String} model
 * @returns {Dictionary} branch information
 */
function estimators.ExtractBranchInformation(tree, node, model) {
    estimators.extractBranchLength.result = {};

    if (Abs(model[terms.model.get_branch_length])) {
        estimators.extractBranchLength.result[terms.fit.MLE] = Call(model[terms.model.get_branch_length], model, tree, node);
    } else {
        estimators.extractBranchLength.result[terms.fit.MLE] = Eval("BranchLength (`tree`, \"`node`\")");
    }

    estimators.extractBranchLength.parameter_tag = tree + "." + node;
    (model.parameters.Local(model))["estimators.ExtractBranchInformation.copy_local"][""];

    return estimators.extractBranchLength.result;
}

/**
 * @name estimators.applyBranchLength
 * @private
 * @param {String} tree
 * @param {String} node
 * @param {String} model
 * @param {String} length
 */
function estimators.applyBranchLength(tree, node, model, length) {
    return Call(model[terms.model.set_branch_length], model, length, tree + "." + node);
}

/**
 * @name estimators.fixSubsetOfEstimates.helper
 * @private
 * @param {String} key
 * @param {String} value
 */
function estimators.fixSubsetOfEstimates.helper(key, value) {
    if (value / terms.constraint == FALSE) {
        estimators.fixSubsetOfEstimates.fixed + value[terms.id];
        value[terms.fix] = TRUE;
    }
}

/**
 * @name estimators.fixSubsetOfEstimates.helper_condition
 * @private
 * @param {String} key
 */
function estimators.fixSubsetOfEstimates.helper_condition(key) {
    return Type(variables[key]) != Unknown;
}

/**
 * @name estimators.fixSubsetOfEstimates
 * @private
 * @param {String} estimates
 * @param {String} variables
 * @returns nothing
 */
function estimators.fixSubsetOfEstimates(estimates, variables) {
    estimators.fixSubsetOfEstimates.fixed = {};
    (estimates[terms.global])["estimators.fixSubsetOfEstimates.helper"]["estimators.fixSubsetOfEstimates.helper_condition"];
    return estimators.fixSubsetOfEstimates.fixed;
}

/**
 * @name estimators.branch_lengths_in_string.map
 * @private
 * @param {String} id
 * @param {String} value
 */
function estimators.branch_lengths_in_string.map(id, value) {
    estimators.branch_lengths_in_string.lookup[id] = value[terms.fit.MLE];
}

/**
 * @name estimators.branch_lengths_in_string
 * @param {String} tree_id
 * @param {String} lookup
 * @returns {String} branch lengths in tree string
 */
function estimators.branch_lengths_in_string(tree_id, lookup) {
    estimators.branch_lengths_in_string.lookup = {};
    lookup["estimators.branch_lengths_in_string.map"][""];
    utility.ToggleEnvVariable("BRANCH_LENGTH_STENCIL", estimators.branch_lengths_in_string.lookup);
    estimators.branch_lengths_in_string.string = Eval("Format (`tree_id`,1,1)");
    utility.ToggleEnvVariable("BRANCH_LENGTH_STENCIL", None);
    return estimators.branch_lengths_in_string.string;
}

/**
 * @name estimators.ExtractMLEs
 * @param {String} likelihood_function_id
 * @param {String} model_descriptions
 * @returns results
 */
function estimators.ExtractMLEs(likelihood_function_id, model_descriptions) {

    ExecuteCommands("GetString (estimators.ExtractMLEs.lfInfo, `likelihood_function_id`,-1)");

    estimators.ExtractMLEs.results = {};
    estimators.ExtractMLEs.partitions = utility.Array1D(estimators.ExtractMLEs.lfInfo[terms.fit.trees]);

    // copy global variables first

    estimators.ExtractMLEs.results[terms.global] = {};
    model_descriptions["estimators.copyGlobals"][""];
    estimators.ExtractMLEs.results[terms.efv_estimate] = {};
    model_descriptions["estimators.CopyFrequencies"][""];
    estimators.ExtractMLEs.results[terms.branch_length] = {};
    estimators.ExtractMLEs.results[terms.fit.trees] = estimators.ExtractMLEs.lfInfo[terms.fit.trees];


    for (estimators.ExtractMLEs.i = 0; estimators.ExtractMLEs.i < estimators.ExtractMLEs.partitions; estimators.ExtractMLEs.i += 1) {
        _tree_name = (estimators.ExtractMLEs.lfInfo[terms.fit.trees])[estimators.ExtractMLEs.i];


        GetInformation (estimators.ExtractMLEs.map, *_tree_name);
        estimators.ExtractMLEs.branch_names = Rows(estimators.ExtractMLEs.map);
        (estimators.ExtractMLEs.results[terms.branch_length])[estimators.ExtractMLEs.i] = {};

        for (estimators.ExtractMLEs.b = 0; estimators.ExtractMLEs.b < Abs(estimators.ExtractMLEs.map); estimators.ExtractMLEs.b += 1) {
            _branch_name = estimators.ExtractMLEs.branch_names[estimators.ExtractMLEs.b];

            ((estimators.ExtractMLEs.results[terms.branch_length])[estimators.ExtractMLEs.i])[_branch_name] =
            estimators.ExtractBranchInformation(_tree_name, _branch_name, model_descriptions[estimators.ExtractMLEs.map[_branch_name]]);

        }


        (estimators.ExtractMLEs.results[terms.fit.trees])[estimators.ExtractMLEs.i] =
        estimators.branch_lengths_in_string((estimators.ExtractMLEs.results[terms.fit.trees])[estimators.ExtractMLEs.i], (estimators.ExtractMLEs.results[terms.branch_length])[estimators.ExtractMLEs.i]);

    }

    return estimators.ExtractMLEs.results;
}

/**
 * @name estimators.TraverseLocalParameters
 * @param {String} likelihood_function_id
 * @param {Dictionary} model_descriptions
 * @param {String} callback (tree, node, parameter_list)

 */
lfunction estimators.TraverseLocalParameters (likelihood_function_id, model_descriptions, callback) {
    GetString(lf_info, ^ likelihood_function_id, -1);
    partitions = utility.Array1D(lf_info[utility.getGlobalValue ("terms.fit.trees")]);
    result = {};
    for (i = 0; i < partitions; i += 1) {
        tree_name = (lf_info[utility.getGlobalValue ("terms.fit.trees")])[i];
        GetInformation (map, ^tree_name);
        branch_names = Rows (map);
        for (b = 0; b < Abs(map); b += 1) {
            _branch_name = branch_names[b];
            result[_branch_name] = Call (callback, tree_name, _branch_name, (model_descriptions [map[_branch_name]])[utility.getGlobalValue ("terms.parameters")]);
        }
    }
    return result;
}


/**
 * @name estimators.ApplyExistingEstimates
 * @param {String} likelihood_function_id
 * @param {Dictionary} model_descriptions
 * @param {Matrix} initial_values
 * @param branch_length_conditions
 * @returns estimators.ApplyExistingEstimates.df_correction - Abs(estimators.ApplyExistingEstimates.keep_track_of_proportional_scalers);
 */
function estimators.ApplyExistingEstimates(likelihood_function_id, model_descriptions, initial_values, branch_length_conditions) {
    //fprintf (stdout, model_descriptions, "\n", initial_values, "\n");

	/* set all category variable values to one */

    GetString(estimators.ApplyExistingEstimates.lfInfo, ^ likelihood_function_id, -1);
    estimators.ApplyExistingEstimates.results = {};
    estimators.ApplyExistingEstimates.partitions = utility.Array1D(estimators.ApplyExistingEstimates.lfInfo[terms.fit.trees]);


    estimators.ApplyExistingEstimates.df_correction = 0;
    // copy global variables first

    estimators.ApplyExistingEstimates.results[terms.global] = {};
    model_descriptions["estimators.SetCategory"][""];
    // the above line traverses all model descriptions and sets
    // the _value_ of category variables to 1, so that we can
    // compute branch lengths

    model_descriptions["estimators.SetGlobals"][""];


    if (Type(branch_length_conditions) == "String") {
        if (branch_length_conditions == terms.globals_only) {
            return estimators.ApplyExistingEstimates.df_correction;
        }
        assert("0", "Unsupported value for 'branch_length_conditions' in estimators.ApplyExistingEstimates");
        return 0;
    }

    estimators.ApplyExistingEstimates.keep_track_of_proportional_scalers = {};

    for (estimators.ApplyExistingEstimates.i = 0; estimators.ApplyExistingEstimates.i < estimators.ApplyExistingEstimates.partitions; estimators.ApplyExistingEstimates.i += 1) {

        if (Type((initial_values[terms.branch_length])[estimators.ApplyExistingEstimates.i]) == "AssociativeList") { // have branch lengths for this partition

            _tree_name = (estimators.ApplyExistingEstimates.lfInfo[terms.fit.trees])[estimators.ApplyExistingEstimates.i];

            ExecuteCommands("GetInformation (estimators.ApplyExistingEstimates.map, `_tree_name`);");
            estimators.ApplyExistingEstimates.branch_names = Rows(estimators.ApplyExistingEstimates.map);

            for (estimators.ApplyExistingEstimates.b = 0; estimators.ApplyExistingEstimates.b < Abs(estimators.ApplyExistingEstimates.map); estimators.ApplyExistingEstimates.b += 1) {
                _branch_name = estimators.ApplyExistingEstimates.branch_names[estimators.ApplyExistingEstimates.b];

                if ((initial_values[terms.branch_length])[estimators.ApplyExistingEstimates.i] / _branch_name) { // have an entry for this branch name
                   _existing_estimate = ((initial_values[terms.branch_length])[estimators.ApplyExistingEstimates.i])[_branch_name];

                   if (Type(_existing_estimate) == "AssociativeList") {
                       _set_branch_length_to = (((initial_values[terms.branch_length])[estimators.ApplyExistingEstimates.i])[_branch_name])[terms.fit.MLE];
                        if (None != branch_length_conditions) {
                            if (Abs(branch_length_conditions)) {
                                _application_type = branch_length_conditions[estimators.ApplyExistingEstimates.i];

                                if (Type(_application_type) == "String") {
                                    _set_branch_length_to = {};
                                    _set_branch_length_to[terms.branch_length] = _existing_estimate[terms.fit.MLE];
                                    _set_branch_length_to[terms.model.branch_length_scaler] = _application_type;
                                    estimators.ApplyExistingEstimates.keep_track_of_proportional_scalers[_application_type] = 1;
                                }
                            }
                        }

                        estimators.ApplyExistingEstimates.df_correction += estimators.applyBranchLength(_tree_name, _branch_name, model_descriptions[estimators.ApplyExistingEstimates.map[_branch_name]], _set_branch_length_to);
                    } else {
                        if (Type(_existing_estimate) != "Unknown") {
                            warning.log ("Incorrect type for the initial values object of for branch '" + _branch_name + "' : " + _existing_estimate);
                        }
                   }
                }
            }

        } else {
        	if (Type((initial_values[terms.branch_length])[estimators.ApplyExistingEstimates.i]) != "Unknown") {
        		warning.log ("Incorrect type for the initial values object for partition " + estimators.ApplyExistingEstimates.i
        					+ ". " + (initial_values[terms.branch_length])[estimators.ApplyExistingEstimates.i]);
        	}
        }
    } // have branch lengths for this partition


    return estimators.ApplyExistingEstimates.df_correction - Abs(estimators.ApplyExistingEstimates.keep_track_of_proportional_scalers);
}

/**
 * @name estimators._aux.countEmpiricalParameters
 * @private
 * @param {String} id
 * @param {Dictionary} model
 * @returns nothing
 */
function estimators._aux.countEmpiricalParameters(id, model) {
    estimators._aux.parameter_counter += (model[terms.parameters])[terms.model.empirical];
}

/**
 * Fits a LikelihoodFunction
 * @name estimators.FitExistingLF
 * @param model_map
 * @returns LF results
 */
lfunction estimators.FitExistingLF (lf_id, model_objects) {

    utility.ToggleEnvVariable("USE_LAST_RESULTS", TRUE);
    Optimize (mles, ^lf_id);
    utility.ToggleEnvVariable("USE_LAST_RESULTS", None);

    results = estimators.ExtractMLEs( lf_id, model_objects);
    results[utility.getGlobalValue ("terms.fit.log_likelihood")] = mles[1][0];
    results[utility.getGlobalValue ("terms.parameters")] = mles[1][1];

    return results;
}

/**
 * Fits a LikelihoodFunction
 * @name estimators.FitLF
 * @param {Matrix} data_filters_list  - a vector of {DataFilter}s
 * @param {Matrix} tree_list  - a vector of {Tree}s
 * @param model_map
 * @param initial_values
 * @returns LF results
 */
lfunction estimators.FitLF(data_filter, tree, model_map, initial_values, model_objects, run_options) {

    if (Type(data_filter) == "String") {
        return estimators.FitLF ({
            {
                data_filter__
            }
        }, {
            "0": tree
        },
        {
            "0" : model_map
        },
        initial_values, model_objects, run_options);
    }

    components = utility.Array1D(data_filter);


    lf_components = {
        2 * components,
        1
    };


    for (i = 0; i < components; i += 1) {
        lf_components[2 * i] = data_filter[i];
        lf_components[2 * i + 1] = &tree_id + "_" + i;
        model.ApplyModelToTree(lf_components[2*i + 1], tree[i], None, model_map[i]);
    }




    lf_id = &likelihoodFunction;
    utility.ExecuteInGlobalNamespace ("LikelihoodFunction `lf_id` = (`&lf_components`)");


    df = 0;

    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", 1);
            df = estimators.ApplyExistingEstimates("`&likelihoodFunction`", model_objects, initial_values, run_options[utility.getGlobalValue("terms.run_options.proportional_branch_length_scaler")]);
    }

    if (utility.Has (run_options,utility.getGlobalValue("terms.run_options.apply_user_constraints"),"String")) {
        df += Call (run_options[utility.getGlobalValue("terms.run_options.apply_user_constraints")], lf_id, lf_components, data_filter, tree, model_map, initial_values, model_objects);
    }


    //assert (0);

    //utility.SetEnvVariable ("VERBOSITY_LEVEL", 10);

   	Optimize (mles, likelihoodFunction);

    //Export (lf,likelihoodFunction);
    //console.log (lf);

    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", None);
    }

    results = estimators.ExtractMLEs( & likelihoodFunction, model_objects);

    results[utility.getGlobalValue ("terms.fit.log_likelihood")] = mles[1][0];
    results[utility.getGlobalValue ("terms.parameters")] = mles[1][1] + df;

    results[utility.getGlobalValue ("terms.fit.filters")] = {
        1,
        components
    };

    for (i = 0; i < components; i += 1) {
        (results[utility.getGlobalValue ("terms.fit.filters")])[i] = lf_components[2 * i];

    }

    if (run_options[utility.getGlobalValue("terms.run_options.retain_lf_object")]) {
        results[utility.getGlobalValue("terms.likelihood_function")] = & likelihoodFunction;
    } else {
        DeleteObject(likelihoodFunction);
    }

    return results;
}

lfunction estimators.CreateLFObject (context, data_filter, tree, model_template, initial_values, run_options) {
    if (Type(data_filter) == "String") {
        return estimators.FitSingleModel_Ext ({
            {
                data_filter__
            }
        }, {
            "0": tree
        }, model_template, initial_values, run_options)
    }

    components = utility.Array1D(data_filter);

    filters = utility.Map({
        components,
        1
    }["_MATRIX_ELEMENT_ROW_"], "_value_", "''+ '`context`.nuc_data_' + _value_");

    lf_components = {
        2 * components,
        1
    };


    for (i = 0; i < components; i += 1) {
        lf_components[2 * i] = filters[i];
        DataSetFilter ^ (filters[i]) = CreateFilter( ^ (data_filter[i]), 1);
    }

    user_model_id = context + ".user_model";
    utility.ExecuteInGlobalNamespace ("`user_model_id` = 0");
    ^(user_model_id) = model.generic.DefineModel(model_template, context + ".model", {
            "0": "terms.global"
        }, filters, None);


    for (i = 0; i < components; i += 1) {
        lf_components[2 * i + 1] = "`context`.tree_" + i;
        model.ApplyModelToTree(lf_components[2 * i + 1], tree[i], {
            "default": ^(user_model_id)
        }, None);
    }


    lfid = context + ".likelihoodFunction";
    utility.ExecuteInGlobalNamespace ("LikelihoodFunction `lfid` = (`&lf_components`)");

    df = 0;
    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", 1);
            df = estimators.ApplyExistingEstimates(lfid, {
                (^user_model_id)[utility.getGlobalValue ("terms.id")]: ^(user_model_id)
            }, initial_values, run_options[utility.getGlobalValue("terms.run_options.proportional_branch_length_scaler")]);
    }

    return df;
}

/**
 * @name estimators.FitSingleModel_Ext
 * @param {DataFilter} data_filter
 * @param {Tree} tree
 * @param {Dict} model
 * @param {Matrix} initial_values
 * @param {Dict} run_options
 * @returns results
 */

lfunction estimators.FitSingleModel_Ext (data_filter, tree, model_template, initial_values, run_options) {

    this_namespace = (&_);
    this_namespace = this_namespace[0][Abs (this_namespace)-3];

    df = estimators.CreateLFObject (this_namespace, data_filter, tree, model_template, initial_values, run_options);

   	Optimize(mles, likelihoodFunction);

    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", None);
    }

    model_id_to_object = {
        (this_namespace + ".model"): user_model
    };


    results = estimators.ExtractMLEs( & likelihoodFunction, model_id_to_object);


    results[utility.getGlobalValue("terms.fit.log_likelihood")] = mles[1][0];
    results[utility.getGlobalValue("terms.parameters")] = mles[1][1] + (user_model [utility.getGlobalValue("terms.parameters")]) [utility.getGlobalValue("terms.model.empirical")] + df;


    if (option[utility.getGlobalValue("terms.run_options.retain_model_object")]) {
        results[utility.getGlobalValue("terms.model")] = model_id_to_object;
    }

   if (run_options[utility.getGlobalValue("terms.run_options.retain_lf_object")]) {
        results[utility.getGlobalValue("terms.likelihood_function")] = & likelihoodFunction;
    } else {
        DeleteObject(likelihoodFunction);
    }

    return results;
}

/**
 * @name estimators.FitGTR_Ext
 * @param {DataFilter} data_filter
 * @param {Tree} tree
 * @param {Matrix} initial_values
 * @param {Dict} run_options
 * @returns results
 */

lfunction estimators.FitGTR_Ext (data_filter, tree, initial_values, run_options) {
    return estimators.FitSingleModel_Ext (data_filter, tree, "models.DNA.GTR.ModelDescription", initial_values, run_options)
}

/**
 * @name estimators.FitGTR
 * @param {DataFilter} data_filter
 * @param {Tree} tree
 * @param {Matrix} initial_values
 * @returns results
 */
lfunction estimators.FitGTR(data_filter, tree, initial_values) {
    return estimators.FitGTR_Ext(data_filter, tree, initial_values, {});
}

/**
 * @name estimators.FitMGREV.set_partition_omega
 * @private
 * @param {String} key
 * @param {String} value
 */
function estimators.FitMGREV.set_partition_omega(key, value) {
    Eval("estimators.FitMGREV.tree.`key`.`estimators.FitMGREV.alpha` = 0.1");
    ExecuteCommands("estimators.FitMGREV.tree.`key`.`estimators.FitMGREV.beta`:=estimators.FitMGREV.tree.`key`.`estimators.FitMGREV.alpha`*" + estimators.FitMGREV.partitioned_omega.parameters[value]);
}

/**
 * @name estimators.FitMGREV.set_partition_omega
 * @private
 * @param {String} key
 * @param {String} value
 */
function estimators.FitMGREVExtractComponentBranchLengths(codon_data, fit_results) {

    //extract fitted trees with branch lengths scaled on synonymous and non-synonymous
    //substitutions per site

    estimators.FitMGREVExtractComponentBranchLengths.stencils = genetic_code.ComputeBranchLengthStencils(codon_data[terms.code]);


    BRANCH_LENGTH_STENCIL = estimators.FitMGREVExtractComponentBranchLengths.stencils["synonymous"];
    fit_results["synonymous-trees"] = (estimators.ExtractMLEs(fit_results[terms.likelihood_function], fit_results[terms.model]))[terms.fit.trees];

    BRANCH_LENGTH_STENCIL = estimators.FitMGREVExtractComponentBranchLengths.stencils["non-synonymous"];
    fit_results["non-synonymous-trees"] = (estimators.ExtractMLEs(fit_results[terms.likelihood_function], fit_results[terms.model]))[terms.fit.trees];

    BRANCH_LENGTH_STENCIL = None;

    return fit_results;
}

/**
 * @name estimators.FitMGREV
 * @param {DataFilter} codon_data
 * @param {Tree} tree
 * @param {String} genetic_code
 * @param {Dictionary} option
 * @param {Dictionary} initial_values
 * @returns MGREV results
 */
lfunction estimators.FitMGREV(codon_data, tree, genetic_code, option, initial_values) {



    //TODO: Where is data_filter being set?
    if (Type(data_filter) == "String") {
        return estimators.FitMGREV({
                {
                    codon_data__
                }
            }, {
                "0": tree
            },
            genetic_code,
            option,
            initial_values)
    }

    components = utility.Array1D(codon_data);

    lf_components = {
        2 * components,
        1
    };



    for (i = 0; i < components; i += 1) {
        GetDataInfo(fi, ^ (codon_data[i]), "PARAMETERS");
        DataSetFilter * ("filter_" + i) = CreateFilter( ^ (codon_data[i]), 3, '', '', fi["EXCLUSIONS"]);
        // need to do this for global references
        lf_components[2 * i] = "filter_" + i;
    }


    name_space = & model_MGREV;

    mg_rev = model.generic.DefineModel("models.codon.MG_REV.ModelDescription",
        name_space, {
            "0": parameters.Quote(option[utility.getGlobalValue("terms.run_options.model_type")]),
            "1": genetic_code
        },
        codon_data,
        None);




    //utility.ToggleEnvVariable("VERBOSITY_LEVEL", 1);

    df = 0;
    model_assignment = {
        "default": mg_rev
    };
    rules = None;
    model_id_to_object = {
        name_space: mg_rev
    };

    for (i = 0; i < components; i += 1) {
        lf_components[2 * i + 1] = "tree_" + i;
        model.ApplyModelToTree(Eval("&`lf_components[2*i + 1]`"), tree[i], model_assignment, None);
    }


    partition_omega = {};

    if (option[utility.getGlobalValue("terms.run_options.model_type")] == utility.getGlobalValue("terms.local") && Type(option[utility.getGlobalValue("terms.run_options.partitioned_omega")]) == "AssociativeList") {
        /**
            Assumes that option["partitioned-omega"] is a dictionary where each partition has
            an entry (0-index based), which itself is a dictionary of the form: "branch-name" : "branch-set"
        */
        utility.ForEach(option[utility.getGlobalValue("terms.run_options.partitioned_omega")], "_value_", "utility.AddToSet(`&partition_omega`,utility.Values(_value_))");
    }



    if (Abs(partition_omega)) {

        /**
            declare the global ratios for each branch set
            and add them to the model parameter set
        */



        new_globals = {};
        utility.ForEachPair(partition_omega, "_key_", "_value_",
            '`&new_globals` [_key_] = (`&name_space` + ".omega_" + Abs (`&new_globals`)); model.generic.AddGlobal (`&mg_rev`, `&new_globals` [_key_] , (utility.getGlobalValue("terms.parameters.omega_ratio")) + " for *" + _key_ + "*")');
        parameters.DeclareGlobal(new_globals, None);


        /**
            now replicate the local constraint for individual branches
        */

        alpha = model.generic.GetLocalParameter(mg_rev, utility.getGlobalValue("terms.parameters.synonymous_rate"));
        beta = model.generic.GetLocalParameter(mg_rev, utility.getGlobalValue("terms.parameters.nonsynonymous_rate"));
        io.CheckAssertion("None!=`&alpha` && None!=`&beta`", "Could not find expected local synonymous and non-synonymous rate parameters in \`estimators.FitMGREV\`");

        apply_constraint: = component_tree + "." + node_name + "." + beta + ":=" + component_tree + "." + node_name + "." + alpha + "*" + new_globals[branch_map[node_name]];

        for (i = 0; i < components; i += 1) {
            component_tree = lf_components[2 * i + 1];
            ClearConstraints( * component_tree);
            branch_map = (option[utility.getGlobalValue("terms.run_options.partitioned_omega")])[i];
            component_branches = BranchName( * component_tree, -1);
            for (j = 0; j < Columns(component_branches) - 1; j += 1) {
                /**
                    -1 in the upper bound because we don't want to count the root node
                */

                node_name = (component_branches[j]);
                ExecuteCommands(apply_constraint);
            }
        }
    } else {}


    LikelihoodFunction likelihoodFunction = (lf_components);



    //fprintf (stdout, option["proportional-branch-length-scaler"], "\n");

    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", 1);
        df += estimators.ApplyExistingEstimates("`&likelihoodFunction`", model_id_to_object, initial_values, option[utility.getGlobalValue("terms.run_options.proportional_branch_length_scaler")]);
    }

    //Export (lfe, likelihoodFunction);
    //console.log (lfe);

    Optimize(mles, likelihoodFunction);


    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", None);
    }

    results = estimators.ExtractMLEs( & likelihoodFunction, model_id_to_object);


    results[utility.getGlobalValue("terms.fit.log_likelihood")] = mles[1][0];
    results[utility.getGlobalValue("terms.parameters")] = mles[1][1] + (mg_rev [utility.getGlobalValue("terms.parameters")]) [utility.getGlobalValue("terms.model.empirical")] + df;


    if (option[utility.getGlobalValue("terms.run_options.retain_lf_object")]) {
        results[utility.getGlobalValue("terms.likelihood_function")] = & likelihoodFunction;
    } else {
        DeleteObject(likelihoodFunction);
    }

    if (option[utility.getGlobalValue("terms.run_options.retain_model_object")]) {
        results[utility.getGlobalValue("terms.model")] = model_id_to_object;
    }

    return results;
}

/**
 * @name estimators.FitMGREV
 * @description compute the asymptotic (chi^2) p-value for the LRT
 * @param {Number} alternative log likelihood for the alternative (more general model)
 * @param {Number} null log likelihood for the null
 * @param {Number} df degrees of freedom
 * @returns p-value
 */
lfunction estimators.LRT (alternative, null, df) {
    if (alternative > null) {
        return 1-CChi2 (2*(alternative-null), df);
    }
    return 1;
}

/**
 * @name estimators.ComputeLF
 * @description compute the current value of the log likelihood function
 * @param {String} lfid name of the function
 * @returns log likelihood
 */
lfunction estimators.ComputeLF (id) {
	LFCompute (^id,LF_START_COMPUTE);
	LFCompute (^id,logl);
	LFCompute (^id,LF_DONE_COMPUTE);
	return logl;
}
