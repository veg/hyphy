LoadFunctionLibrary("../models/model_functions.bf");
LoadFunctionLibrary("../models/terms.bf");
LoadFunctionLibrary("../models/DNA/GTR.bf");

/**
 * @name estimators.GetGlobalMLE
 * @param {Dictionary} results
 * @param {String} tag
 * @returns None
 */
lfunction estimators.GetGlobalMLE(results, tag) {
    estimate = (results[ ^ "terms.global"])[tag];
    if (Type(estimate) == "AssociativeList") {
        return estimate[ ^ "terms.MLE"];
    }
    return None;
}

/**
 * @name estimators.copyGlobals2
 * @private
 * @param {String} key2
 * @param {String} value2
 * @returns nothing
 */
function estimators.copyGlobals2(key2, value2) {

    if (Type((estimators.ExtractMLEs.results["global"])[key2]) == "AssociativeList") {
        key2 = "[`key`] `key2`";
    }

    (estimators.ExtractMLEs.results["global"])[key2] = {
        "ID": value2,
        "MLE": Eval(value2)
    };

    if (parameters.IsIndependent(value2) != TRUE) {
        ((estimators.ExtractMLEs.results["global"])[key2])["constraint"] = parameters.getConstraint(value2);
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
    ((value["parameters"])["global"])["estimators.copyGlobals2"][""];
}

/**
 * @name estimators.CopyFrequencies
 * @private
 * @param {String} key
 * @param {Dictionary} value
 * @returns nothing
 */
function estimators.CopyFrequencies(model_name, model_decription) {
    (estimators.ExtractMLEs.results[terms.efv_estimate])[model_name] = model_decription[terms.efv_estimate];
}



function estimators.SetGlobals2(key, value) {
    __init_value = (initial_values["global"])[key];
    if (Type(__init_value) == "AssociativeList") {
        if (__init_value["fix-me"]) {
            estimators.ApplyExistingEstimates.df_correction += parameters.IsIndependent(value);
            ExecuteCommands("`value` := " + __init_value["MLE"]);
        } else {
            //fprintf (stdout, "Setting `value` to " + __init_value["MLE"] + "\n");
            ExecuteCommands("`value` = " + __init_value["MLE"]);
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
    ((value["parameters"])["global"])["estimators.SetGlobals2"][""];
}

/**
 * @name estimators.ExtractBranchInformation.copy_local
 * @param {String} key
 * @param {String} value
 */
function estimators.ExtractBranchInformation.copy_local(key, value) {

    estimators.ExtractBranchInformation.copy_local.var_name = estimators.extractBranchLength.parameter_tag + "." + value;

    estimators.extractBranchLength.result[key] = {
        "ID": value,
        "MLE": Eval(estimators.ExtractBranchInformation.copy_local.var_name)
    };

    if (parameters.IsIndependent(estimators.ExtractBranchInformation.copy_local.var_name) != TRUE) {
        (estimators.extractBranchLength.result[key])["constraint"] = parameters.getConstraint(estimators.ExtractBranchInformation.copy_local.var_name);
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

    if (Abs(model["get-branch-length"])) {
        estimators.extractBranchLength.result["MLE"] = Call(model["get-branch-length"], model, tree, node);
    } else {
        estimators.extractBranchLength.result["MLE"] = Eval("BranchLength (`tree`, \"`node`\")");
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
    return Call(model["set-branch-length"], model, length, tree + "." + node);
}

/**
 * @name estimators.fixSubsetOfEstimates.helper
 * @private
 * @param {String} key
 * @param {String} value
 */
function estimators.fixSubsetOfEstimates.helper(key, value) {
    value["fix-me"] = 1;
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
    (estimates["global"])["estimators.fixSubsetOfEstimates.helper"]["estimators.fixSubsetOfEstimates.helper_condition"];
}

/**
 * @name estimators.branch_lengths_in_string.map
 * @private
 * @param {String} id
 * @param {String} value
 */
function estimators.branch_lengths_in_string.map(id, value) {
    estimators.branch_lengths_in_string.lookup[id] = value["MLE"];
}

/**
 * @name estimators.branch_lengths_in_string
 * @param {String} tree_id
 * @param {String} lookup
 * @returns {String} branch lenghts in tree string
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
    estimators.ExtractMLEs.partitions = utility.Array1D(estimators.ExtractMLEs.lfInfo["Trees"]);

    // copy global variables first

    estimators.ExtractMLEs.results["global"] = {};
    model_descriptions["estimators.copyGlobals"][""];
    estimators.ExtractMLEs.results[terms.efv_estimate] = {};
    model_descriptions["estimators.CopyFrequencies"][""];
    estimators.ExtractMLEs.results[terms.json.attribute.branch_length] = {};
    estimators.ExtractMLEs.results["Trees"] = estimators.ExtractMLEs.lfInfo["Trees"];

    for (estimators.ExtractMLEs.i = 0; estimators.ExtractMLEs.i < estimators.ExtractMLEs.partitions; estimators.ExtractMLEs.i += 1) {
        _tree_name = (estimators.ExtractMLEs.lfInfo["Trees"])[estimators.ExtractMLEs.i];

        GetInformation (estimators.ExtractMLEs.map, *_tree_name);
        estimators.ExtractMLEs.branch_names = Rows(estimators.ExtractMLEs.map);
        (estimators.ExtractMLEs.results[terms.json.attribute.branch_length])[estimators.ExtractMLEs.i] = {};

        for (estimators.ExtractMLEs.b = 0; estimators.ExtractMLEs.b < Abs(estimators.ExtractMLEs.map); estimators.ExtractMLEs.b += 1) {
            _branch_name = estimators.ExtractMLEs.branch_names[estimators.ExtractMLEs.b];
            ((estimators.ExtractMLEs.results[terms.json.attribute.branch_length])[estimators.ExtractMLEs.i])[_branch_name] =
            estimators.ExtractBranchInformation(_tree_name, _branch_name, model_descriptions[estimators.ExtractMLEs.map[_branch_name]]);
        }

        (estimators.ExtractMLEs.results["Trees"])[estimators.ExtractMLEs.i] =
        estimators.branch_lengths_in_string((estimators.ExtractMLEs.results["Trees"])[estimators.ExtractMLEs.i], (estimators.ExtractMLEs.results[terms.json.attribute.branch_length])[estimators.ExtractMLEs.i]);
    }

    return estimators.ExtractMLEs.results;
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

    // fprintf (stdout, model_descriptions, "\n", initial_values, "\n");

    GetString(estimators.ApplyExistingEstimates.lfInfo, ^ likelihood_function_id, -1);
    estimators.ApplyExistingEstimates.results = {};
    estimators.ApplyExistingEstimates.partitions = utility.Array1D(estimators.ApplyExistingEstimates.lfInfo["Trees"]);


    estimators.ApplyExistingEstimates.df_correction = 0;
    // copy global variables first

    estimators.ApplyExistingEstimates.results["global"] = {};
    model_descriptions["estimators.SetGlobals"][""];


    if (Type(branch_length_conditions) == "String") {
        if (branch_length_conditions == "globals only") {
            return estimators.ApplyExistingEstimates.df_correction;
        }
        assert("0", "Unsupported value for 'branch_length_conditions' in estimators.ApplyExistingEstimates");
        return 0;
    }

    estimators.ApplyExistingEstimates.keep_track_of_proportional_scalers = {};

    for (estimators.ApplyExistingEstimates.i = 0; estimators.ApplyExistingEstimates.i < estimators.ApplyExistingEstimates.partitions; estimators.ApplyExistingEstimates.i += 1) {

        if (Type((initial_values[terms.json.attribute.branch_length])[estimators.ApplyExistingEstimates.i]) == "AssociativeList") {

            _tree_name = (estimators.ApplyExistingEstimates.lfInfo["Trees"])[estimators.ApplyExistingEstimates.i];

            ExecuteCommands("GetInformation (estimators.ApplyExistingEstimates.map, `_tree_name`);");
            estimators.ApplyExistingEstimates.branch_names = Rows(estimators.ApplyExistingEstimates.map);

            for (estimators.ApplyExistingEstimates.b = 0; estimators.ApplyExistingEstimates.b < Abs(estimators.ApplyExistingEstimates.map); estimators.ApplyExistingEstimates.b += 1) {
                _branch_name = estimators.ApplyExistingEstimates.branch_names[estimators.ApplyExistingEstimates.b];
                _existing_estimate = ((initial_values[terms.json.attribute.branch_length])[estimators.ApplyExistingEstimates.i])[_branch_name];

                if (Type(_existing_estimate) == "AssociativeList") {
                    _set_branch_length_to = _existing_estimate["MLE"];

                    if (None != branch_length_conditions) {
                        if (Abs(branch_length_conditions)) {
                            _application_type = branch_length_conditions[estimators.ApplyExistingEstimates.i];
                            if (Type(_application_type) == "String") {
                                _set_branch_length_to = {};
                                _set_branch_length_to[terms.branch_length] = _existing_estimate["MLE"];
                                _set_branch_length_to[terms.branch_length_scaler] = _application_type;
                                estimators.ApplyExistingEstimates.keep_track_of_proportional_scalers[_application_type] = 1;
                            }
                        }
                    }

                    estimators.ApplyExistingEstimates.df_correction += estimators.applyBranchLength(_tree_name, _branch_name, model_descriptions[estimators.ApplyExistingEstimates.map[_branch_name]], _set_branch_length_to);
                }
            }

        }
    }

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
    estimators._aux.parameter_counter += (model["parameters"])["empirical"];
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
function estimators.FitLF(data_filters_list, tree_list, model_map, initial_values) {
    estimators.FitLF.component_count = utility.Array1D(data_filters_list);

    assert(estimators.FitLF.component_count == utility.Array1D(tree_list),
        "Data filters and tree lists must have the same dimensions in call to estimators.FitLF");


    estimators.FitLF.components = {
        estimators.FitLF.component_count,
        2
    };

    for (estimators.FitLF.i = 0; estimators.FitLF.i < estimators.FitLF.component_count; estimators.FitLF.i += 1) {
        estimators.FitLF.components[estimators.FitLF.i][0] = data_filters_list[estimators.FitLF.i];
        estimators.FitLF.components[estimators.FitLF.i][1] = tree_list[estimators.FitLF.i];
    }

    LikelihoodFunction estimators.FitLF.likelihoodFunction = (estimators.FitLF.components);

    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", 1);
        estimators.ApplyExistingEstimates("estimators.FitLF.likelihoodFunction", model_map, initial_values, None);
    }

    /*Export (boom, estimators.FitLF.likelihoodFunction);
    fprintf (stdout, boom, "\n");
    assert (0);*/
    Optimize(estimators.FitLF.mles, estimators.FitLF.likelihoodFunction);
    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", None);
    }

    estimators.FitLF.results = estimators.ExtractMLEs("estimators.FitLF.likelihoodFunction", model_map);

    estimators._aux.parameter_counter = 0;
    model_map["estimators._aux.countEmpiricalParameters"][""];

    estimators.FitLF.results["LogL"] = estimators.FitLF.mles[1][0];
    estimators.FitLF.results["parameters"] = estimators.FitLF.mles[1][1] + estimators._aux.parameter_counter;
    estimators.FitLF.results["Filters"] = {
        1,
        estimators.FitLF.component_count
    };

    for (estimators.FitLF.i = 0; estimators.FitLF.i < estimators.FitLF.component_count; estimators.FitLF.i += 1) {
        (estimators.FitLF.results["Filters"])[estimators.FitLF.i] = data_filters_list[estimators.FitLF.i];
        //(estimators.FitLF.results["Trees"])[estimators.FitLF.i]   = Eval ("Format("+tree_list[estimators.FitLF.i]+",1,1)");
    }


    DeleteObject(estimators.FitLF.likelihoodFunction);

    return estimators.FitLF.results;
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

    if (Type(data_filter) == "String") {
        return estimators.FitSingleModel_Ext ({
            {
                data_filter__
            }
        }, {
            "0": tree
        }, model, initial_values, run_options)
    }

    components = utility.Array1D(data_filter);

    filters = utility.Map({
        components,
        1
    }["_MATRIX_ELEMENT_ROW_"], "_value_", "''+ '`&nuc_data`_' + _value_");

    lf_components = {
        2 * components,
        1
    };


    for (i = 0; i < components; i += 1) {
        lf_components[2 * i] = filters[i];
        DataSetFilter ^ (filters[i]) = CreateFilter( ^ (data_filter[i]), 1);
    }

    name_space = & user;

    user_model = model.generic.DefineModel(model_template, name_space, {
            "0": "terms.global"
        }, filters, None);

    for (i = 0; i < components; i += 1) {

        lf_components[2 * i + 1] = "tree_" + i;
        model.ApplyModelToTree(Eval("&`lf_components[2*i + 1]`"), tree[i], {
            "default": user_model
        }, None);
    }

    LikelihoodFunction likelihoodFunction = (lf_components);

    df = 0;
    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", 1);
        df = estimators.ApplyExistingEstimates("`&likelihoodFunction`", {
            name_space: user_model
        }, initial_values, run_options["proportional-branch-length-scaler"]);
    }

    Optimize(mles, likelihoodFunction);

    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", None);
    }


    results = estimators.ExtractMLEs( & likelihoodFunction, {
        name_space: user_model
    });

    results["LogL"] = mles[1][0];
    results["parameters"] = mles[1][1] + 3 + df;

    if (run_options["retain-lf-object"]) {
        results["LF"] = & likelihoodFunction;
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

    estimators.FitMGREVExtractComponentBranchLengths.stencils = genetic_code.ComputeBranchLengthStencils(codon_data["code"]);


    BRANCH_LENGTH_STENCIL = estimators.FitMGREVExtractComponentBranchLengths.stencils["synonymous"];
    fit_results["synonymous-trees"] = (estimators.ExtractMLEs(fit_results["LF"], fit_results["model"]))["Trees"];

    BRANCH_LENGTH_STENCIL = estimators.FitMGREVExtractComponentBranchLengths.stencils["non-synonymous"];
    fit_results["non-synonymous-trees"] = (estimators.ExtractMLEs(fit_results["LF"], fit_results["model"]))["Trees"];

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
        lf_components[2 * i] = "filter_" + i;
        // need to do this for global references
    }

    name_space = & model_MGREV;

    mg_rev = model.generic.DefineModel("models.codon.MG_REV.ModelDescription",
        name_space, {
            "0": parameters.Quote(option["model-type"]),
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
    if (option["model-type"] == ^ "terms.local" && Type(option["partitioned-omega"]) == "AssociativeList") {
        /**
            Assumes that option["partitioned-omega"] is a dictionary where each partition has
            an entry (0-index based), which itself is a dictionary of the form: "branch-name" : "branch-set"
        */
        utility.ForEach(option["partitioned-omega"], "_value_", "utility.AddToSet(`&partition_omega`,utility.Values(_value_))");
    }


    if (Abs(partition_omega)) {

        /**
            declare the global ratios for each branch set
            and add them to the model parameter set
        */

        new_globals = {};
        utility.ForEachPair(partition_omega, "_key_", "_value_",
            '`&new_globals` [_key_] = (`&name_space` + ".omega_" + Abs (`&new_globals`)); model.generic.AddGlobal (`&mg_rev`, `&new_globals` [_key_] , (^"terms.omega_ratio") + " for *" + _key_ + "*")');
        parameters.DeclareGlobal(new_globals, None);


        /**
            now replicate the local constraint for individual branches
        */

        alpha = model.generic.GetLocalParameter(mg_rev, ^ "terms.synonymous_rate");
        beta = model.generic.GetLocalParameter(mg_rev, ^ "terms.nonsynonymous_rate");
        io.CheckAssertion("None!=`&alpha` && None!=`&beta`", "Could not find expected local synonymous and non-synonymous rate parameters in \`estimators.FitMGREV\`");


        apply_constraint: = component_tree + "." + node_name + "." + beta + ":=" + component_tree + "." + node_name + "." + alpha + "*" + new_globals[branch_map[node_name && 1]];

        for (i = 0; i < components; i += 1) {
            component_tree = lf_components[2 * i + 1];
            ClearConstraints( * component_tree);
            branch_map = (option["partitioned-omega"])[i];
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
        df += estimators.ApplyExistingEstimates("`&likelihoodFunction`", model_id_to_object, initial_values, option["proportional-branch-length-scaler"]);
    }

    //Export (lfs, likelihoodFunction);
    //fprintf (stdout, lfs, "\n");

    Optimize(mles, likelihoodFunction);

    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", None);
    }

    results = estimators.ExtractMLEs( & likelihoodFunction, model_id_to_object);

    results["LogL"] = mles[1][0];
    results["parameters"] = mles[1][1] + 9 + df; /* 9 frequency parameters */

    if (option["retain-lf-object"]) {
        results["LF"] = & likelihoodFunction;
    } else {
        DeleteObject(likelihoodFunction);
    }

    if (option["retain-model-object"]) {
        results["model"] = model_id_to_object;
    }

    return results;
}
