LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("terms.bf");

parameters.infinity = 1e10;

function parameters.applyNameSpace (id, namespace) {
	if (Type (namespace) == "String") {
		if (Abs (namespace) > 0) {
			return namespace + "." + id;
		}
	}
	return id;
}

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

function parameters.declareGlobal (id, cache) {
    if (Type (id) == "String") {
        if (Abs (id)) {
            if (Type (cache) == "AssociativeList") {
                if (Abs (cache[id]) == 0) {
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
                    parameters.declareGlobal (parameters.declareGlobal.names[parameters.declareGlobal.k], cache);
                 }                       
            }
        }
    }
}

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

function parameters.set_value (id, value) {
    Eval ("`id` = " + value);
}

lfunction parameters.mean (values, weights, d) {
    m = 0;
    for (i = 0; i < 3; i+=1) {
        m += Eval (values[i]) * Eval(weights[i]);
    }
    return m;
}

function parameters.quote (arg) {
	return "\"" + arg + "\"";
}

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

function parameters.isIndependent (id) {
    ExecuteCommands ("GetString (parameters.isIndependent.t, `id`, -1);");
    return Type (parameters.isIndependent.t) != "AssociativeList";
}

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

function parameters.constrainSets (set1, set2) {
    parameters.constrainSets.tags = Rows (set1);
    for (parameters.constrainSets.k = 0; parameters.constrainSets.k < Abs (set1); parameters.constrainSets.k += 1) {
        
        if (Type (set2[parameters.constrainSets.tags [parameters.constrainSets.k]]) == "String") {
            ExecuteCommands (set2[parameters.constrainSets.tags [parameters.constrainSets.k]] + ":=" +
                             set1[parameters.constrainSets.tags [parameters.constrainSets.k]]);
        }
    }
}


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
