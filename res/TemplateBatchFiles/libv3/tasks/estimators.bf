LoadFunctionLibrary("../models/model_functions.bf");
LoadFunctionLibrary("../models/terms.bf");
LoadFunctionLibrary("../models/DNA/GTR.bf");

/**
 * @name estimators.getGlobalMLE
 * @param {Dictionary} results
 * @param {String} tag
 * @returns None
 */
lfunction estimators.getGlobalMLE(results, tag) {
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

    if (Type((estimators.extractMLEs.results["global"])[key2]) == "AssociativeList") {
        key2 = "[`key`] `key2`";
    }

    (estimators.extractMLEs.results["global"])[key2] = {
        "ID": value2,
        "MLE": Eval(value2)
    };

    if (parameters.isIndependent(value2) != TRUE) {
        ((estimators.extractMLEs.results["global"])[key2])["constraint"] = parameters.getConstraint(value2);
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

function estimators.setGlobals2(key, value) {
    __init_value = (initial_values["global"])[key];
    if (Type(__init_value) == "AssociativeList") {
        if (__init_value["fix-me"]) {
            estimators.applyExistingEstimates.df_correction += parameters.isIndependent(value);
            ExecuteCommands("`value` := " + __init_value["MLE"]);
        } else {
            ExecuteCommands("`value` = " + __init_value["MLE"]);
        }
    }
}

/**
 * @name estimators.setGlobals
 * @param {String} key
 * @param {String} value
 * @returns nothing
 */
function estimators.setGlobals(key, value) {
    ((value["parameters"])["global"])["estimators.setGlobals2"][""];
}

/**
 * @name estimators.extractBranchInformation.copy_local
 * @param {String} key
 * @param {String} value
 */
function estimators.extractBranchInformation.copy_local(key, value) {

    estimators.extractBranchInformation.copy_local.var_name = estimators.extractBranchLength.parameter_tag + "." + value;

    estimators.extractBranchLength.result[key] = {
        "ID": value,
        "MLE": Eval(estimators.extractBranchInformation.copy_local.var_name)
    };

    if (parameters.isIndependent(estimators.extractBranchInformation.copy_local.var_name) != TRUE) {
        (estimators.extractBranchLength.result[key])["constraint"] = parameters.getConstraint(estimators.extractBranchInformation.copy_local.var_name);
    }

}

/**
 * @name estimators.extractBranchInformation
 * @param {String} type
 * @param {String} node
 * @param {String} model
 * @returns {Dictionary} branch information
 */
function estimators.extractBranchInformation(tree, node, model) {
    estimators.extractBranchLength.result = {};

    if (Abs(model["get-branch-length"])) {
        estimators.extractBranchLength.result["MLE"] = Call(model["get-branch-length"], model, tree, node);
    } else {
        estimators.extractBranchLength.result["MLE"] = Eval("BranchLength (`tree`, \"`node`\")");
    }

    estimators.extractBranchLength.parameter_tag = tree + "." + node;
    (model.parameters.local(model))["estimators.extractBranchInformation.copy_local"][""];

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
    utility.toggleEnvVariable("BRANCH_LENGTH_STENCIL", estimators.branch_lengths_in_string.lookup);
    estimators.branch_lengths_in_string.string = Eval("Format (`tree_id`,1,1)");
    utility.toggleEnvVariable("BRANCH_LENGTH_STENCIL", None);
    return estimators.branch_lengths_in_string.string;
}

/**
 * @name estimators.extractMLEs
 * @param {String} likelihood_function_id
 * @param {String} model_descriptions
 * @returns results
 */
function estimators.extractMLEs(likelihood_function_id, model_descriptions) {

    ExecuteCommands("GetString (estimators.extractMLEs.lfInfo, `likelihood_function_id`,-1)");

    estimators.extractMLEs.results = {};
    estimators.extractMLEs.partitions = utility.array1D(estimators.extractMLEs.lfInfo["Trees"]);

    // copy global variables first

    estimators.extractMLEs.results["global"] = {};
    model_descriptions["estimators.copyGlobals"][""];
    estimators.extractMLEs.results[terms.json.attribute.branch_length] = {};
    estimators.extractMLEs.results["Trees"] = estimators.extractMLEs.lfInfo["Trees"];

    for (estimators.extractMLEs.i = 0; estimators.extractMLEs.i < estimators.extractMLEs.partitions; estimators.extractMLEs.i += 1) {
        _tree_name = (estimators.extractMLEs.lfInfo["Trees"])[estimators.extractMLEs.i];

        ExecuteCommands("GetInformation (estimators.extractMLEs.map, `_tree_name`);");
        estimators.extractMLEs.branch_names = Rows(estimators.extractMLEs.map);
        (estimators.extractMLEs.results[terms.json.attribute.branch_length])[estimators.extractMLEs.i] = {};

        for (estimators.extractMLEs.b = 0; estimators.extractMLEs.b < Abs(estimators.extractMLEs.map); estimators.extractMLEs.b += 1) {
            _branch_name = estimators.extractMLEs.branch_names[estimators.extractMLEs.b];
            ((estimators.extractMLEs.results[terms.json.attribute.branch_length])[estimators.extractMLEs.i])[_branch_name] =
            estimators.extractBranchInformation(_tree_name, _branch_name, model_descriptions[estimators.extractMLEs.map[_branch_name]]);
        }

        (estimators.extractMLEs.results["Trees"])[estimators.extractMLEs.i] =
        estimators.branch_lengths_in_string((estimators.extractMLEs.results["Trees"])[estimators.extractMLEs.i], (estimators.extractMLEs.results[terms.json.attribute.branch_length])[estimators.extractMLEs.i]);
    }

    return estimators.extractMLEs.results;
}

/**
 * @name estimators.applyExistingEstimates
 * @param {String} likelihood_function_id 
 * @param {Dictionary} model_descriptions 
 * @param {Matrix} initial_values 
 * @param branch_length_conditions
 * @returns estimators.applyExistingEstimates.df_correction - Abs(estimators.applyExistingEstimates.keep_track_of_proportional_scalers);
 */
function estimators.applyExistingEstimates(likelihood_function_id, model_descriptions, initial_values, branch_length_conditions) {

    // fprintf (stdout, model_descriptions, "\n", initial_values, "\n");

    GetString(estimators.applyExistingEstimates.lfInfo, ^ likelihood_function_id, -1);
    estimators.applyExistingEstimates.results = {};
    estimators.applyExistingEstimates.partitions = utility.array1D(estimators.applyExistingEstimates.lfInfo["Trees"]);


    estimators.applyExistingEstimates.df_correction = 0;
    // copy global variables first

    estimators.applyExistingEstimates.results["global"] = {};
    model_descriptions["estimators.setGlobals"][""];


    if (Type(branch_length_conditions) == "String") {
        if (branch_length_conditions == "globals only") {
            return estimators.applyExistingEstimates.df_correction;
        }
        assert("0", "Unsupported value for 'branch_length_conditions' in estimators.applyExistingEstimates");
        return 0;
    }

    estimators.applyExistingEstimates.keep_track_of_proportional_scalers = {};

    for (estimators.applyExistingEstimates.i = 0; estimators.applyExistingEstimates.i < estimators.applyExistingEstimates.partitions; estimators.applyExistingEstimates.i += 1) {

        if (Type((initial_values[terms.json.attribute.branch_length])[estimators.applyExistingEstimates.i]) == "AssociativeList") {

            _tree_name = (estimators.applyExistingEstimates.lfInfo["Trees"])[estimators.applyExistingEstimates.i];

            ExecuteCommands("GetInformation (estimators.applyExistingEstimates.map, `_tree_name`);");
            estimators.applyExistingEstimates.branch_names = Rows(estimators.applyExistingEstimates.map);

            for (estimators.applyExistingEstimates.b = 0; estimators.applyExistingEstimates.b < Abs(estimators.applyExistingEstimates.map); estimators.applyExistingEstimates.b += 1) {
                _branch_name = estimators.applyExistingEstimates.branch_names[estimators.applyExistingEstimates.b];
                _existing_estimate = ((initial_values[terms.json.attribute.branch_length])[estimators.applyExistingEstimates.i])[_branch_name];

                if (Type(_existing_estimate) == "AssociativeList") {
                    _set_branch_length_to = _existing_estimate["MLE"];

                    if (None != branch_length_conditions) {
                        if (Abs(branch_length_conditions)) {
                            _application_type = branch_length_conditions[estimators.applyExistingEstimates.i];
                            if (Type(_application_type) == "String") {
                                _set_branch_length_to = {};
                                _set_branch_length_to[terms.branch_length] = _existing_estimate["MLE"];
                                _set_branch_length_to[terms.branch_length_scaler] = _application_type;
                                estimators.applyExistingEstimates.keep_track_of_proportional_scalers[_application_type] = 1;
                            }
                        }
                    }

                    estimators.applyExistingEstimates.df_correction += estimators.applyBranchLength(_tree_name, _branch_name, model_descriptions[estimators.applyExistingEstimates.map[_branch_name]], _set_branch_length_to);
                }
            }

        }
    }

    return estimators.applyExistingEstimates.df_correction - Abs(estimators.applyExistingEstimates.keep_track_of_proportional_scalers);
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
 * @name estimators.fitLF
 * @param {Matrix} data_filters_list  - a vector of {DataFilter}s
 * @param {Matrix} tree_list  - a vector of {Tree}s
 * @param model_map 
 * @param initial_values
 * @returns LF results
 */
function estimators.fitLF(data_filters_list, tree_list, model_map, initial_values) {
    estimators.fitLF.component_count = utility.array1D(data_filters_list);

    assert(estimators.fitLF.component_count == utility.array1D(tree_list),
        "Data filters and tree lists must have the same dimensions in call to estimators.fitLF");


    estimators.fitLF.components = {
        estimators.fitLF.component_count,
        2
    };

    for (estimators.fitLF.i = 0; estimators.fitLF.i < estimators.fitLF.component_count; estimators.fitLF.i += 1) {
        estimators.fitLF.components[estimators.fitLF.i][0] = data_filters_list[estimators.fitLF.i];
        estimators.fitLF.components[estimators.fitLF.i][1] = tree_list[estimators.fitLF.i];
    }

    LikelihoodFunction estimators.fitLF.likelihoodFunction = (estimators.fitLF.components);

    if (Type(initial_values) == "AssociativeList") {
        utility.toggleEnvVariable("USE_LAST_RESULTS", 1);
        estimators.applyExistingEstimates("estimators.fitLF.likelihoodFunction", model_map, initial_values, None);
    }

    /*Export (boom, estimators.fitLF.likelihoodFunction);
    fprintf (stdout, boom, "\n");
    assert (0);*/
    Optimize(estimators.fitLF.mles, estimators.fitLF.likelihoodFunction);
    if (Type(initial_values) == "AssociativeList") {
        utility.toggleEnvVariable("USE_LAST_RESULTS", None);
    }

    estimators.fitLF.results = estimators.extractMLEs("estimators.fitLF.likelihoodFunction", model_map);

    estimators._aux.parameter_counter = 0;
    model_map["estimators._aux.countEmpiricalParameters"][""];

    estimators.fitLF.results["LogL"] = estimators.fitLF.mles[1][0];
    estimators.fitLF.results["parameters"] = estimators.fitLF.mles[1][1] + estimators._aux.parameter_counter;
    estimators.fitLF.results["Filters"] = {
        1,
        estimators.fitLF.component_count
    };

    for (estimators.fitLF.i = 0; estimators.fitLF.i < estimators.fitLF.component_count; estimators.fitLF.i += 1) {
        (estimators.fitLF.results["Filters"])[estimators.fitLF.i] = data_filters_list[estimators.fitLF.i];
        //(estimators.fitLF.results["Trees"])[estimators.fitLF.i]   = Eval ("Format("+tree_list[estimators.fitLF.i]+",1,1)");
    }


    DeleteObject(estimators.fitLF.likelihoodFunction);

    return estimators.fitLF.results;
}

/**
 * @name estimators.fitGTR
 * @param {DataFilter} data_filter
 * @param {Tree} tree
 * @param {Matrix} initial_values
 * @returns results
 */
lfunction estimators.fitGTR(data_filter, tree, initial_values) {

    if (Type(data_filter) == "String") {
        return estimators.fitGTR({
            {
                data_filter__
            }
        }, {
            "0": tree
        }, initial_values)
    }

    components = utility.array1D(data_filter);

    filters = utility.map({
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

    //utility.setEnvVariable ("VERBOSITY_LEVEL", 10);

    name_space = & gtr;

    gtr_model = model.generic.define_model("models.DNA.GTR.modelDescription", name_space, {
        "0": "terms.global"
    }, filters, None);


    for (i = 0; i < components; i += 1) {

        lf_components[2 * i + 1] = "tree_" + i;
        model.applyModelToTree(Eval("&`lf_components[2*i + 1]`"), tree[i], {
            "default": gtr_model
        }, None);
    }

    LikelihoodFunction likelihoodFunction = (lf_components);

    df = 0;
    if (Type(initial_values) == "AssociativeList") {
        utility.toggleEnvVariable("USE_LAST_RESULTS", 1);
        df = estimators.applyExistingEstimates("`&likelihoodFunction`", {
            name_space: gtr_model
        }, initial_values, None);
    }

    Optimize(mles, likelihoodFunction);

    if (Type(initial_values) == "AssociativeList") {
        utility.toggleEnvVariable("USE_LAST_RESULTS", None);
    }


    results = estimators.extractMLEs( & likelihoodFunction, {
        name_space: gtr_model
    });

    results["LogL"] = mles[1][0];
    results["parameters"] = mles[1][1] + 3 + df;

    DeleteObject(likelihoodFunction);

    return results;
}

/**
 * @name estimators.fitMGREV.set_partition_omega
 * @private
 * @param {String} key
 * @param {String} value
 */
function estimators.fitMGREV.set_partition_omega(key, value) {
    Eval("estimators.fitMGREV.tree.`key`.`estimators.fitMGREV.alpha` = 0.1");
    ExecuteCommands("estimators.fitMGREV.tree.`key`.`estimators.fitMGREV.beta`:=estimators.fitMGREV.tree.`key`.`estimators.fitMGREV.alpha`*" + estimators.fitMGREV.partitioned_omega.parameters[value]);
}

/**
 * @name estimators.fitMGREV.set_partition_omega
 * @private
 * @param {String} key
 * @param {String} value
 */
function estimators.fitMGREVExtractComponentBranchLengths(codon_data, fit_results) {

    //extract fitted trees with branch lengths scaled on synonymous and non-synonymous
    //substitutions per site

    estimators.fitMGREVExtractComponentBranchLengths.stencils = genetic_code.ComputeBranchLengthStencils(codon_data["code"]);


    BRANCH_LENGTH_STENCIL = estimators.fitMGREVExtractComponentBranchLengths.stencils["synonymous"];
    fit_results["synonymous-trees"] = (estimators.extractMLEs(fit_results["LF"], fit_results["model"]))["Trees"];

    BRANCH_LENGTH_STENCIL = estimators.fitMGREVExtractComponentBranchLengths.stencils["non-synonymous"];
    fit_results["non-synonymous-trees"] = (estimators.extractMLEs(fit_results["LF"], fit_results["model"]))["Trees"];

    BRANCH_LENGTH_STENCIL = None;

    return fit_results;
}

/**
 * @name estimators.fitMGREV
 * @param {DataFilter} codon_data
 * @param {Tree} tree
 * @param {String} genetic_code
 * @param {Dictionary} option
 * @param {Dictionary} initial_values
 * @returns MGREV results
 */
lfunction estimators.fitMGREV(codon_data, tree, genetic_code, option, initial_values) {

    //TODO: Where is data_filter being set?
    if (Type(data_filter) == "String") {
        return estimators.fitMGREV({
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

    components = utility.array1D(codon_data);

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

    mg_rev = model.generic.define_model("models.codon.MG_REV.modelDescription",
        name_space, {
            "0": parameters.quote(option["model-type"]),
            "1": genetic_code
        },
        codon_data,
        None);

    //utility.toggleEnvVariable("VERBOSITY_LEVEL", 1);

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
        model.applyModelToTree(Eval("&`lf_components[2*i + 1]`"), tree[i], model_assignment, None);
    }


    partition_omega = {};
    if (option["model-type"] == ^ "terms.local" && Type(option["partitioned-omega"]) == "AssociativeList") {
        /**
            Assumes that option["partitioned-omega"] is a dictionary where each partition has
            an entry (0-index based), which itself is a dictionary of the form: "branch-name" : "branch-set"
        */
        utility.forEach(option["partitioned-omega"], "_value_", "utility.addToSet(`&partition_omega`,utility.values(_value_))");
    }


    if (Abs(partition_omega)) {

        /**
            declare the global ratios for each branch set
            and add them to the model parameter set
        */

        new_globals = {};
        utility.forEachPair(partition_omega, "_key_", "_value_",
            '`&new_globals` [_key_] = (`&name_space` + ".omega_" + Abs (`&new_globals`)); model.generic.add_global (`&mg_rev`, `&new_globals` [_key_] , (^"terms.omega_ratio") + " for *" + _key_ + "*")');
        parameters.declareGlobal(new_globals, None);


        /**
            now replicate the local constraint for individual branches
        */

        alpha = model.generic.get_local_parameter(mg_rev, ^ "terms.synonymous_rate");
        beta = model.generic.get_local_parameter(mg_rev, ^ "terms.nonsynonymous_rate");
        io.checkAssertion("None!=`&alpha` && None!=`&beta`", "Could not find expected local synonymous and non-synonymous rate parameters in \`estimators.fitMGREV\`");


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
        utility.toggleEnvVariable("USE_LAST_RESULTS", 1);
        df += estimators.applyExistingEstimates("`&likelihoodFunction`", model_id_to_object, initial_values, option["proportional-branch-length-scaler"]);
    }

    //Export (lfs, likelihoodFunction);
    //fprintf (stdout, lfs, "\n");

    Optimize(mles, likelihoodFunction);

    if (Type(initial_values) == "AssociativeList") {
        utility.toggleEnvVariable("USE_LAST_RESULTS", None);
    }

    results = estimators.extractMLEs( & likelihoodFunction, model_id_to_object);

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
