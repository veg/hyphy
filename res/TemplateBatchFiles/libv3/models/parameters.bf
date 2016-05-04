LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("terms.bf");

parameters.infinity = 1e10;

/**
 * Applies a namespace to parameter ids
 * @param {String} id
 * @param {String} namespace
 */
function parameters.applyNameSpace (id, namespace) {
	if (Type (namespace) == "String") {
		if (Abs (namespace) > 0) {
			return namespace + "." + id;
		}
	}
	return id;
}

/**
 * parameters.unconstrain_parameter_set 
 * @param {LikelihoodFunction} lf - the likelihood function to operate on
 * @param {Matrix} set - set of parameters to unconstrain
 * @returns nothing
 */
function parameters.unconstrain_parameter_set (lf, set) {
    ExecuteCommands ("GetString(parameters.unconstrain_parameter_set.info, `lf`, -1)");
    if (None == set) {
        set = {{terms.lf.global.constrained,terms.lf.local.constrained}};
    }
    for (parameters.unconstrain_parameter_set.s = 0; parameters.unconstrain_parameter_set.s < Columns (set); parameters.unconstrain_parameter_set.s += 1) {
        parameters.unconstrain_parameter_set.m = parameters.unconstrain_parameter_set.info [set [parameters.unconstrain_parameter_set.s]];
        for (parameters.unconstrain_parameter_set.i = 0; parameters.unconstrain_parameter_set.i < Columns (parameters.unconstrain_parameter_set.m); parameters.unconstrain_parameter_set.i += 1) {
            Eval (parameters.unconstrain_parameter_set.m[parameters.unconstrain_parameter_set.i] + "=" + Eval (parameters.unconstrain_parameter_set.m[parameters.unconstrain_parameter_set.i]));
        }
    }
}

/**
 * parameters.declareGlobal 
 * @param {String} id
 * @param {Matrix} cache
 * @returns nothing
 */
function parameters.declareGlobal (id, cache) {
    if (Type (id) == "String") {
        if (Abs (id)) {
            if (Type (cache) == "AssociativeList") {
                if (Abs (cache[id]) > 0) {
                    return;
                } else {
                    cache[id] = 1;
                }
            }
            ExecuteCommands ("global `id` = 1;");
        }
    } else {
        if (Type (id) == "AssociativeList") {
            parameters.declareGlobal.var_count = Abs (id);
            parameters.declareGlobal.names = Columns (id);
            for (parameters.declareGlobal.k = 0; parameters.declareGlobal.k <  parameters.declareGlobal.var_count; parameters.declareGlobal.k += 1) {
                parameters.declareGlobal (parameters.declareGlobal.names[parameters.declareGlobal.k], cache);
            }
        } else {
            if (Type (id) == "Matrix") {
                 parameters.declareGlobal.var_count = Columns (id) * Rows (id);
                 for (parameters.declareGlobal.k = 0; parameters.declareGlobal.k <  parameters.declareGlobal.var_count; parameters.declareGlobal.k += 1) {
                    parameters.declareGlobal (id[parameters.declareGlobal.k], cache);
                 }
            }
        }
    }
}

/**
 * parameters.normalize_ratio
 * @param {Number} n 
 * @param {Number} d
 * @returns n/d
 */
function parameters.normalize_ratio (n, d) {
    if (d == 0) {
        if (n == 0) {
            return 1;
        } else {
            return parameters.infinity;
        }
    }
    return n/d;
}

/**
 * Sets value of passed parameter id
 * @param {String} id - id of parameter to set value to
 * @param {Number} value - value to set
 * @returns nothing
 */
function parameters.set_value (id, value) {
    Eval ("`id` = " + value);
}

/**
 * Returns mean of values
 * @param {Matrix} values - values to return mean of
 * @param {Matrix} weights - weights to multiply values by
 * @param {Number} d - does nothing
 * @returns {Number} mean
 */
lfunction parameters.mean (values, weights, d) {
    m = 0;
    d = Rows  (values)*Columns (values);
    for (i = 0; i < d; i+=1) {
        m += Eval (values[i]) * Eval(weights[i]);
    }
    return m;
}

/**
 * Quotes the argument
 * @param {String} arg - string to be quoted
 * @returns {String} string in quotes
 */
function parameters.quote (arg) {
    if (Type (arg) == "String") {
	    return "\"" + (arg&&2) + "\"";
    }
    return arg;
}

/**
 * addMultiplicativeTerm 
 * @param {Matrix} matrix - matrix to scale
 * @param {Number} term - scalar to multiply matrix by
 * @param {Number} do_empties - if element matrix is empty, fill with term
 * @returns {Matrix} New matrix
 */
lfunction parameters.addMultiplicativeTerm (matrix, term, do_empties) {

	if (Abs (term) > 0) {
		__N = Rows (matrix);

		for (__r = 0; __r < __N; __r+=1) {
			for (__c = 0; __c < __N; __c+=1) {
				if (__r != __c) {
					if (Abs (matrix[__r][__c])) {
						matrix[__r][__c] = "(" + matrix[__r][__c] + ")*(" + term + ")";
					} else {
					    if (do_empties) {
						    matrix[__r][__c] =  term;
						}
					}
				}
			}
		}
	}

	return matrix;
}

/**
 * stringMatrixToFormulas 
 * @param {String} id - matrix to scale
 * @param {Matrix} matrix - if element matrix is empty, fill with term
 * @returns nothing
 */
function parameters.stringMatrixToFormulas (id, matrix) {
	__N = Rows (matrix);

	ExecuteCommands ("`id` = {__N,__N}");

	for (__r = 0; __r < __N; __r+=1) {
		for (__c = 0; __c < __N; __c+=1) {

			if (__r != __c && Abs (matrix[__r][__c])) {
				ExecuteCommands ("`id`[__r][__c] := " + matrix[__r][__c]);
			}
		}
	}

}

/**
 * parameters.generate_attributed_names 
 * @param {String} prefix 
 * @param {Dictionary} attributes 
 * @param {String} delimiter
 */
function parameters.generate_attributed_names (prefix, attributes, delimiter) {
    if (delimiter == None) {
        delimiter = "_";
    }
    parameters.generate_names.holder = {};
    for (parameters.generate_names.k = 0; parameters.generate_names.k < Columns (attributes); parameters.generate_names.k += 1) {
        parameters.generate_names.holder + (prefix + delimiter + attributes[parameters.generate_names.k]);
    }
    return parameters.generate_names.holder;
}

/**
 * parameters.generate_sequential_names 
 * @param {String} prefix 
 * @param {Number} count 
 * @param {String} delimiter
 * @returns {Matrix} 1 x <count> row vector of generated names
 */
function parameters.generate_sequential_names (prefix, count, delimiter) {
    if (delimiter == None) {
        delimiter = "_";
    }
    parameters.generate_names.holder = {};
    for (parameters.generate_names.k = 0; parameters.generate_names.k < count; parameters.generate_names.k += 1) {
        parameters.generate_names.holder + (prefix + delimiter + parameters.generate_names.k);
    }
    return parameters.generate_names.holder;
}

/**
 * parameters.setRange 
 * @param id
 * @param ranges
 * @returns nothing
 */
function parameters.setRange (id, ranges) {
    if (Type (id) == "String") {
        if (Abs (id)) {
            if (Type (ranges) == "AssociativeList") {
                if (Abs (ranges[terms.lower_bound])) {
                    ExecuteCommands ("`id` :> " + ranges[terms.lower_bound]);
                }
                if (Abs (ranges[terms.upper_bound])) {
                    ExecuteCommands ("`id` :< " + ranges[terms.upper_bound]);
                }
            }
        }
    } else {
        if (Type (id) == "AssociativeList") {
            parameters.setRange.var_count = Abs (id);
            for (parameters.setRange.k = 0; parameters.setRange.k <  parameters.setRange.var_count; parameters.setRange.k += 1) {
                parameters.setRange (id[parameters.setRange.k], ranges);
            }
        }
    }
}

/**
 * Check if parameter is independent
 * @param parameter - id of parameter to check
 * @returns {BOOL} TRUE if independent, FALSE otherwise
 */
lfunction parameters.isIndependent (parameter) {
    GetString (info, ^parameter, -1);
    if (Type (info) == "AssociativeList") {
        return (utility.checkKey (info, "Local", "Matrix") && utility.checkKey (info, "Global", "Matrix")) == FALSE;
    }
    return TRUE;
}

lfunction parameters.getConstraint (parameter) {
    GetString (info, ^parameter, -2);
    return info;
}

/**
 * sets constraint on parameter
 * @param {String} or {AssociativeList} id - id(s) of parameter(s) to set constraint on
 * @param {Number} value - the constraint to set on the parameter
 * @param {String} global_tag - the global namespace of the parameter
 * @returns nothing
 */
function parameters.setConstraint (id, value, global_tag) {
    if (Type (id) == "String") {
        if (Abs (id)) {
            ExecuteCommands ("`global_tag` `id` := " + value);
        }
    } else {
        if (Type (id) == "AssociativeList" && Type (value) == "AssociativeList") {

            parameters.setConstraint.var_count = Abs (id);
            for (parameters.setConstraint.k = 0; parameters.setConstraint.k <  parameters.setConstraint.var_count; parameters.setConstraint.k += 1) {
                parameters.setConstraint (id[parameters.setConstraint.k],
                                          value[parameters.setConstraint.k],
                                          global_tag);
            }
        }
    }
}

/**
 * constraint set of parameters
 * @param {AssociativeList} set1 - 
 * @param {AssociativeList} set2 - 
 * @returns nothing
 */
function parameters.constrainSets (set1, set2) {
    parameters.constrainSets.tags = Rows (set1);
    for (parameters.constrainSets.k = 0; parameters.constrainSets.k < Abs (set1); parameters.constrainSets.k += 1) {

        if (Type (set2[parameters.constrainSets.tags [parameters.constrainSets.k]]) == "String") {
            ExecuteCommands (set2[parameters.constrainSets.tags [parameters.constrainSets.k]] + ":=" +
                             set1[parameters.constrainSets.tags [parameters.constrainSets.k]]);
        }
    }
}

/**
 * Removes a constraint from a parameter
 * @param {String} id - id of parameter to remove constraint from
 * @returns nothing
 */
function parameters.removeConstraint (id) {
    if (Type (id) == "String") {
        if (Abs (id)) {
            Eval ("`id` = " + Eval(id));
        }
    } else {
        if (Type (id) == "AssociativeList") {
            return parameters.removeConstraint (Columns (id));
        }
        if (Type (id) == "Matrix") {
            parameters.removeConstraint.var_count = Columns (id) * Rows (id);
            for (parameters.removeConstraint.k = 0; parameters.removeConstraint.k <  parameters.removeConstraint.var_count; parameters.removeConstraint.k += 1) {
                parameters.removeConstraint (id[parameters.removeConstraint.k]);
            }
        }
    }
}


/**
 * Copies definitions from target to source
 * @param {Dictionary} target - the target dictionary 
 * @param {Dictionary} source - the source element to copy to target
 * @returns nothing
 */
function parameters.helper.copy_definitions (target, source) {
    parameters.helper.copy_definitions.key_iterator = {{terms.local, terms.global}};

    for (parameters.helper.copy_definitions.i = 0;
         parameters.helper.copy_definitions.i < Columns (parameters.helper.copy_definitions.key_iterator);
         parameters.helper.copy_definitions.i += 1) {
         parameters.helper.copy_definitions.key = parameters.helper.copy_definitions.key_iterator[parameters.helper.copy_definitions.i];
         if (Type (source[parameters.helper.copy_definitions.key]) == "AssociativeList") {
            target [parameters.helper.copy_definitions.key] * source [parameters.helper.copy_definitions.key];
         }
    }
}

/**
 * parameters.helper.stick_breaking 
 * @param {AssociativeList} parameters
 * @param {Matrix} initial_values
 * @returns weights
 */
lfunction parameters.helper.stick_breaking (parameters, initial_values) {
    left_over   = "";
    weights     = {};
    accumulator = 1;


    for (k = 0; k < Abs (parameters); k += 1) {
        if (None != initial_values) {
            vid = parameters[k];
            ^(vid) = initial_values[k] / accumulator;
            accumulator = accumulator * (1-^(vid));
        }
        weights [k] = left_over + parameters[k];
        left_over += "(1-" + parameters[k] + ")*";
     }

    weights[k] = left_over[0][Abs (left_over)-2];
    return weights;
}

lfunction parameters.helper.dump_matrix (matrix) {
    for (i = 0; i < Rows (^matrix); i+=1) {
        for (j = 0; j < Columns (^matrix); j+=1) {
            ExecuteCommands ("GetString (cell, `matrix`, i, j)");
            fprintf (stdout, "`matrix`[", i, "][", j, "] := ", cell, "\n");
        }
    }
    return None;
}

/**
 * Sets tree lengths to initial values
 * @param dict - a [0 to N-1] dictrionary of tree objects 
 * @param type - codon or nucleotide
 * @returns {Dictionary} dictionary of initial branch lengths
 */
lfunction parameters.helper.tree_lengths_to_initial_values (dict, type) {

    components = Abs (dict);

    //result = {"branch lengths" : { "0" : {} } };

    if (type == "codon") {
        factor = 1;
    } else {
        factor = 1;
    }

    result  = {};

    for (i = 0; i < components; i += 1) {
        //((result["branch lengths"])[0])[keys[i]] = {"MLE": factor * dict[keys[i]]};
        this_component = {};
        utility.forEachPair ((dict[i])[^"terms.json.attribute.branch_length"], "_branch_name_", "_branch_length_", "`&this_component`[_branch_name_] = {^'terms.json.MLE' : `&factor`*_branch_length_}");
        result[i] = this_component;
    }

    return {^"terms.json.attribute.branch_length" : result} ;
}

/**
 * Profiles likelihood function based on covariance precision level
 * @param {String} id - covariance parameter id
 * @param {LikelihoodFunction} lf - likelihood function to profile
 * @param {Number} - Covariance precision level 
 * @returns {Dictionary} a dictionary containing profiling information
 */
function parameters.getProfileCI (id, lf, level) {

    utility.toggleEnvVariable ("COVARIANCE_PRECISION", level);
    utility.toggleEnvVariable ("COVARIANCE_PARAMETER", id);
    CovarianceMatrix (parameters.getProfileCI.mx, *lf);
    utility.toggleEnvVariable ("COVARIANCE_PRECISION", None);
    utility.toggleEnvVariable ("COVARIANCE_PARAMETER", None);

    return {"`terms.lower_bound`" : parameters.getProfileCI.mx[0], "`terms.MLE`" : parameters.getProfileCI.mx[1],"`terms.upper_bound`" : parameters.getProfileCI.mx[2]};
}


