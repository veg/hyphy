LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("libv3/all-terms.bf");

/** @module parameters */

parameters.infinity = 1e10;

/**
 * Applies a namespace to parameter ids
 * @name parameters.ApplyNameSpace
 * @param {String} id
 * @param {String} namespace
 */
function parameters.ApplyNameSpace(id, namespace) {
    if (Type(namespace) == "String") {
        if (Abs(namespace) > 0) {
            return namespace + "." + id;
        }
    }
    return id;
}

/**
 * @name parameters.UnconstrainParameterSet
 * @param {LikelihoodFunction} lf - the likelihood function to operate on
 * @param {Matrix} set - set of parameters to unconstrain
 * @returns nothing
 */
function parameters.UnconstrainParameterSet(lf, set) {
    ExecuteCommands("GetString(parameters.UnconstrainParameterSet.info, `lf`, -1)");
    if (None == set) {
        set = {
            {
                terms.parameters.global_constrained, terms.parameters.local_constrained
            }
        };
    }
    for (parameters.UnconstrainParameterSet.s = 0; parameters.UnconstrainParameterSet.s < Columns(set); parameters.UnconstrainParameterSet.s += 1) {
        parameters.UnconstrainParameterSet.m = parameters.UnconstrainParameterSet.info[set[parameters.UnconstrainParameterSet.s]];
        for (parameters.UnconstrainParameterSet.i = 0; parameters.UnconstrainParameterSet.i < Columns(parameters.UnconstrainParameterSet.m); parameters.UnconstrainParameterSet.i += 1) {
            Eval(parameters.UnconstrainParameterSet.m[parameters.UnconstrainParameterSet.i] + "=" + Eval(parameters.UnconstrainParameterSet.m[parameters.UnconstrainParameterSet.i]));
        }
    }
}

/**
 * @name parameters.DeclareGlobal
 * @param {String} id
 * @param {Matrix} cache
 * @returns nothing
 */
function parameters.DeclareGlobal(id, cache) {
    if (Type(id) == "String") {
        if (Abs(id)) {
            if (Type(cache) == "AssociativeList") {
                if (Abs(cache[id]) > 0) {
                    return;
                } else {
                    cache[id] = 1;
                }
            }
            ExecuteCommands("global `id` = 1;");
        }
    } else {
        if (Type(id) == "AssociativeList") {
            parameters.DeclareGlobal.var_count = Abs(id);
            parameters.DeclareGlobal.names = Columns(id);
            for (parameters.DeclareGlobal.k = 0; parameters.DeclareGlobal.k < parameters.DeclareGlobal.var_count; parameters.DeclareGlobal.k += 1) {
                parameters.DeclareGlobal(parameters.DeclareGlobal.names[parameters.DeclareGlobal.k], cache);
            }
        } else {
            if (Type(id) == "Matrix") {
                parameters.DeclareGlobal.var_count = Columns(id) * Rows(id);
                for (parameters.DeclareGlobal.k = 0; parameters.DeclareGlobal.k < parameters.DeclareGlobal.var_count; parameters.DeclareGlobal.k += 1) {
                    parameters.DeclareGlobal(id[parameters.DeclareGlobal.k], cache);
                }
            }
        }
    }
}

/**
 * @name parameters.DeclareGlobalWithRanges
 * @param {String} id variable id
 * @param {Number} init initial value (could be None)
 * @param {Number} lb lower bound (could be None)
 * @param {Number} ub upper bound (could be None)
 * @returns nothing
 */
function parameters.DeclareGlobalWithRanges(id, init, lb, ub) {
	if (Type (id) == "String") {
		if (None != init) {
    		ExecuteCommands("global `id` = " + init);
    	} else {
    		ExecuteCommands("global `id`; ");    	
    	}
    	
    	if (None != lb) {
    		ExecuteCommands ("`id` :> " + lb);
    	}	
    	if (None != ub) {
    		ExecuteCommands ("`id` :< " + ub);
    	}	
    } else {
        if (Type(id) == "AssociativeList") {
            parameters.DeclareGlobalWithRanges.var_count = Abs(id);
            parameters.DeclareGlobalWithRanges.names = Columns(id);
            for (parameters.DeclareGlobalWithRanges.k = 0; parameters.DeclareGlobalWithRanges.k < parameters.DeclareGlobalWithRanges.var_count; parameters.DeclareGlobalWithRanges.k += 1) {
                parameters.DeclareGlobalWithRanges(parameters.DeclareGlobalWithRanges.names[parameters.DeclareGlobalWithRanges.k], init, lb, ub);
            }
        } else {
            if (Type(id) == "Matrix") {
                parameters.DeclareGlobalWithRanges.var_count = Columns(id) * Rows(id);
                for (parameters.DeclareGlobalWithRanges.k = 0; parameters.DeclareGlobalWithRanges.k < parameters.DeclareGlobalWithRanges.var_count; parameters.DeclareGlobalWithRanges.k += 1) {
                    parameters.DeclareGlobalWithRanges(id[parameters.DeclareGlobalWithRanges.k], init, lb, ub);
                }
            }
        }
    }
}

function parameters.DeclareCategory.helper (dict, key, default) {
	if (dict / key) {
		return dict[key];
	}
	return default;
}

/**
 * @name parameters.DeclareCategory
 * @param {Dict} def category definition components
 */
function parameters.DeclareCategory (def) {
	 
	 
	 ExecuteCommands ("category " + def[terms.id] + "= (" + 
	 			  Join (",", 
	 			  			utility.Map ({"0": terms.category.bins, "1": terms.category.weights, "2": terms.category.represent, "3": terms.category.PDF, "4": terms.category.CDF, "5": terms.lower_bound, "6": terms.upper_bound, "7": terms.category.dCDF}, 
	 			  						  "_value_",
	 			  						  'parameters.DeclareCategory.helper(def[terms.category.category_parameters], _value_, "")')
	 			  		) + ");");
	 
}


/**
 * @name parameters.NormalizeRatio
 * @param {Number} n
 * @param {Number} d
 * @returns n/d
 */
function parameters.NormalizeRatio(n, d) {
    if (d == 0) {
        if (n == 0) {
            return 1;
        } else {
            return parameters.infinity;
        }
    }
    return n / d;
}

/**
 * Sets value of passed parameter id
 * @name parameters.SetValue
 * @param {String} id - id of parameter to set value to
 * @param {Number} value - value to set
 * @returns nothing
 */
function parameters.SetValue(id, value) {
    Eval("`id` = " + value);
}

/**
 * Returns mean of values
 * @name parameters.Mean
 * @param {Matrix} values - values to return mean of
 * @param {Matrix} weights - weights to multiply values by
 * @param {Number} d - does nothing
 * @returns {Number} mean
 */
lfunction parameters.Mean(values, weights, d) {
    m = 0;
    d = Rows(values) * Columns(values);
    for (i = 0; i < d; i += 1) {
        m += Eval(values[i]) * Eval(weights[i]);
    }
    return m;
}

/**
 * Quotes the argument
 * @name parameters.Quote
 * @param {String} arg - string to be quoted
 * @returns {String} string in quotes
 */
function parameters.Quote(arg) {
    if (Type(arg) == "String") {
        return "\"" + (arg && 2) + "\"";
    }
    return arg;
}

/**
 * @name parameters.AddMultiplicativeTerm
 * @param {Matrix} matrix - matrix to scale
 * @param {Number} term - scalar to multiply matrix by
 * @param {Number} do_empties - if element matrix is empty, fill with term
 * @returns {Matrix} New matrix
 */
lfunction parameters.AddMultiplicativeTerm(matrix, term, do_empties) {

    if (Abs(term) > 0) {
        __N = Rows(matrix);

        for (__r = 0; __r < __N; __r += 1) {
            for (__c = 0; __c < __N; __c += 1) {
                if (__r != __c) {
                    if (Abs(matrix[__r][__c])) {
                        matrix[__r][__c] = "(" + matrix[__r][__c] + ")*(" + term + ")";
                    } else {
                        if (do_empties) {
                            matrix[__r][__c] = term;
                        }
                    }
                }
            }
        }
    }

    return matrix;
}

/**
 * @name parameters.StringMatrixToFormulas
 * @param {String} id - matrix to scale
 * @param {Matrix} matrix - if element matrix is empty, fill with term
 * @returns nothing
 */
function parameters.StringMatrixToFormulas(id, matrix) {
    __N = Rows(matrix);

    ExecuteCommands("`id` = {__N,__N}");

    for (__r = 0; __r < __N; __r += 1) {
        for (__c = 0; __c < __N; __c += 1) {

            if (__r != __c && Abs(matrix[__r][__c])) {
                ExecuteCommands("`id`[__r][__c] := " + matrix[__r][__c]);
            }
        }
    }

}

/**
 * @name parameters.GenerateAttributedNames
 * @param {String} prefix
 * @param {Dictionary} attributes
 * @param {String} delimiter
 */
function parameters.GenerateAttributedNames(prefix, attributes, delimiter) {
    if (delimiter == None) {
        delimiter = "_";
    }
    parameters.generate_names.holder = {};
    for (parameters.generate_names.k = 0; parameters.generate_names.k < Columns(attributes); parameters.generate_names.k += 1) {
        parameters.generate_names.holder + (prefix + delimiter + attributes[parameters.generate_names.k]);
    }
    return parameters.generate_names.holder;
}

/**
 * @name parameters.GenerateSequentialNames
 * @param {String} prefix
 * @param {Number} count
 * @name parameters.GenerateSequentialNames
 * @param {String} prefix
 * @param {Number} count
 * @param {String} delimiter
 * @returns {Matrix} 1 x <count> row vector of generated names
 */
lfunction parameters.GenerateSequentialNames(prefix, count, delimiter) {
    if (delimiter == None) {
        delimiter = "_";
    }
    holder = {};
    for (k = 0; k < count; k += 1) {
        holder + (prefix + delimiter + k);
    }
    return holder;
}

/**
 * @name parameters.SetRange
 * @param id
 * @param ranges
 * @returns nothing
 */
function parameters.SetRange(id, ranges) {
    if (Type(id) == "String") {
        if (Abs(id)) {
            if (Type(ranges) == "AssociativeList") {
                if (Abs(ranges[terms.lower_bound])) {
                    ExecuteCommands("`id` :> " + ranges[terms.lower_bound]);
                }
                if (Abs(ranges[terms.upper_bound])) {
                    ExecuteCommands("`id` :< " + ranges[terms.upper_bound]);
                }
            }
        }
    } else {
        if (Type(id) == "AssociativeList") {
            parameters.SetRange.var_count = Abs(id);
            for (parameters.SetRange.k = 0; parameters.SetRange.k < parameters.SetRange.var_count; parameters.SetRange.k += 1) {
                parameters.SetRange(id[parameters.SetRange.k], ranges);
            }
        }
    }
}

/**
 * Check if parameter is independent
 * @name parameters.IsIndependent
 * @param parameter - id of parameter to check
 * @returns {Bool} TRUE if independent, FALSE otherwise
 */
lfunction parameters.IsIndependent(parameter) {

    //console.log(parameter);

    GetString(info, ^ parameter, -1);
    if (Type(info) == "AssociativeList") {
        return (utility.CheckKey(info, "Local", "Matrix") && utility.CheckKey(info, "Global", "Matrix")) == FALSE;
    }
    return TRUE;
}

lfunction parameters.getConstraint(parameter) {
    GetString(info, ^ parameter, -2);
    return info;
}

/**
 * sets constraint on parameter
 * @name parameters.SetConstraint
 * @param {String} or {AssociativeList} id - id(s) of parameter(s) to set constraint on
 * @param {Number} value - the constraint to set on the parameter
 * @param {String} global_tag - the global namespace of the parameter
 * @returns nothing
 */
function parameters.SetConstraint(id, value, global_tag) {
    if (Type(id) == "String") {
        if (Abs(id)) {
            ExecuteCommands("`global_tag` `id` := " + value);
        }
    } else {
        if (Type(id) == "AssociativeList" && Type(value) == "AssociativeList") {
            parameters.SetConstraint.var_count = Abs(id);
            for (parameters.SetConstraint.k = 0; parameters.SetConstraint.k < parameters.SetConstraint.var_count; parameters.SetConstraint.k += 1) {
                parameters.SetConstraint(id[parameters.SetConstraint.k],
                    value[parameters.SetConstraint.k],
                    global_tag);
            }
        }
    }
}

/**
 * constraint set of parameters
 * @name parameters.ConstrainSets
 * @param {AssociativeList} set1 -
 * @param {AssociativeList} set2 -
 * @returns nothing
 */
function parameters.ConstrainSets(set1, set2) {
    parameters.ConstrainSets.tags = Rows(set1);
    for (parameters.ConstrainSets.k = 0; parameters.ConstrainSets.k < Abs(set1); parameters.ConstrainSets.k += 1) {

        if (Type(set2[parameters.ConstrainSets.tags[parameters.ConstrainSets.k]]) == "String") {
            ExecuteCommands(set2[parameters.ConstrainSets.tags[parameters.ConstrainSets.k]] + ":=" +
                set1[parameters.ConstrainSets.tags[parameters.ConstrainSets.k]]);
        }
    }
}

/**
 * Removes a constraint from a parameter
 * @name parameters.RemoveConstraint
 * @param {String} id - id of parameter to remove constraint from
 * @returns nothing
 */
function parameters.RemoveConstraint(id) {
    if (Type(id) == "String") {
        if (Abs(id)) {
            Eval("`id` = " + Eval(id));
        }
    } else {
        if (Type(id) == "AssociativeList") {
            return parameters.RemoveConstraint(Columns(id));
        }
        if (Type(id) == "Matrix") {
            parameters.RemoveConstraint.var_count = Columns(id) * Rows(id);
            for (parameters.RemoveConstraint.k = 0; parameters.RemoveConstraint.k < parameters.RemoveConstraint.var_count; parameters.RemoveConstraint.k += 1) {
                parameters.RemoveConstraint(id[parameters.RemoveConstraint.k]);
            }
        }
    }
}


/**
 * Copies definitions from target to source
 * @name parameters.helper.copy_definitions
 * @param {Dictionary} target - the target dictionary
 * @param {Dictionary} source - the source element to copy to target
 * @returns nothing
 */
function parameters.helper.copy_definitions(target, source) {
    parameters.helper.copy_definitions.key_iterator = {
        {
            terms.local, terms.global
        }
    };

    for (parameters.helper.copy_definitions.i = 0; parameters.helper.copy_definitions.i < Columns(parameters.helper.copy_definitions.key_iterator); parameters.helper.copy_definitions.i += 1) {
        parameters.helper.copy_definitions.key = parameters.helper.copy_definitions.key_iterator[parameters.helper.copy_definitions.i];
        if (Type(source[parameters.helper.copy_definitions.key]) == "AssociativeList") {
            target[parameters.helper.copy_definitions.key] * source[parameters.helper.copy_definitions.key];
        }
    }
    
    if (utility.Has (source, terms.category, "AssociativeList")) {
    	utility.EnsureKey (target, terms.category);
    	(target[terms.category])[(source[terms.category])[terms.id]] = (source[terms.category])[terms.description];
    }
}

/**
 * @name parameters.helper.stick_breaking
 * @param {AssociativeList} parameters
 * @param {Matrix} initial_values
 * @returns weights
 */
lfunction parameters.helper.stick_breaking(parameters, initial_values) {
    left_over = "";
    weights = {};
    accumulator = 1;


    for (k = 0; k < Abs(parameters); k += 1) {
        if (None != initial_values) {
            vid = parameters[k]; ^ (vid) = initial_values[k] / accumulator;
            accumulator = accumulator * (1 - ^ (vid));
        }
        weights[k] = left_over + parameters[k];
        left_over += "(1-" + parameters[k] + ")*";
    }

    weights[k] = left_over[0][Abs(left_over) - 2];
    return weights;
}

/**
 * Prints matrix to screen
 * @name parameters.helper.dump_matrix
 * @param {Matrix} matrix
 * @returns nothing
 */
lfunction parameters.helper.dump_matrix(matrix) {
    for (i = 0; i < Rows( ^ matrix); i += 1) {
        for (j = 0; j < Columns( ^ matrix); j += 1) {
            ExecuteCommands("GetString (cell, `matrix`, i, j)");
            fprintf(stdout, "`matrix`[", i, "][", j, "] := ", cell, "\n");
        }
    }
    return None;
}

/**
 * Sets tree lengths to initial values
 * @name parameters.helper.tree_lengths_to_initial_values
 * @param dict - a [0 to N-1] dictrionary of tree objects
 * @param type - codon or nucleotide
 * @returns {Dictionary} dictionary of initial branch lengths
 */
lfunction parameters.helper.tree_lengths_to_initial_values(dict, type) {

    components = Abs(dict);

    //result = {"branch lengths" : { "0" : {} } };

    if (type == "codon") { // TODO
        factor = 1;
    } else {
        factor = 1;
    }

    result = {};

    for (i = 0; i < components; i += 1) {
        //((result["branch lengths"])[0])[keys[i]] = {"MLE": factor * dict[keys[i]]};
        this_component = {};
        utility.ForEachPair((dict[i])[ utility.getGlobalValue("terms.branch_length")], "_branch_name_", "_branch_length_", "`&this_component`[_branch_name_] = {utility.getGlobalValue('terms.fit.MLE') : `&factor`*_branch_length_}");
        result[i] = this_component;
    }

    return { utility.getGlobalValue("terms.branch_length"): result
    };
}

/**
 * Profiles likelihood function based on covariance precision level
 * @name parameters.GetProfileCI
 * @param {String} id - covariance parameter id
 * @param {LikelihoodFunction} lf - likelihood function to profile
 * @param {Number} - Covariance precision level
 * @returns {Dictionary} a dictionary containing profiling information
 */
function parameters.GetProfileCI(id, lf, level) {

    utility.ToggleEnvVariable("COVARIANCE_PRECISION", level);
    utility.ToggleEnvVariable("COVARIANCE_PARAMETER", id);
    CovarianceMatrix(parameters.GetProfileCI.mx, * lf);
    utility.ToggleEnvVariable("COVARIANCE_PRECISION", None);
    utility.ToggleEnvVariable("COVARIANCE_PARAMETER", None);

    //TODO: used to be "`terms.lower_bound`". 
    return {
        terms.lower_bound: parameters.GetProfileCI.mx[0],
        terms.fit.MLE: parameters.GetProfileCI.mx[1],
        terms.upper_bound: parameters.GetProfileCI.mx[2]
    };
}

/**
 * Geneate an HBL string needed to define a parameter
 * @name parameters.ExportParameterDefinition
 * @param {String} id - the name of the parameter to export;
 * @returns {String} the string for an HBL definition of the parameter
 * e.g. "global x := 2.3; x :> 1; x :< 3"; note that the definition will
 * be NOT recursive, so if x depends on y and z, then y and z need to be
 * exported separately
 */

lfunction parameters.ExportParameterDefinition (id) {
    GetString (parameter_definition, ^id, -3);
    return parameter_definition;
}

/**
 * Copy local parameters from a node to 'template' variables
 * @name parameters.SetLocalModelParameters
 * @param {Dict} model - model description
 * @param {String} tree id
 * @param {String} node id
 
 * e.g. 
 * 		set alpha = Tree.Node.alpha
 * 		set beta  = Tree.Node.beta
 */

lfunction parameters.SetLocalModelParameters (model, tree, node) {
	node_name = tree + "." + node + ".";
    utility.ForEach ((model[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.local")], "_parameter_", '
    	^_parameter_ = ^(`&node_name` + _parameter_);
    ');
    
}

