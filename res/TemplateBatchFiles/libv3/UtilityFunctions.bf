LoadFunctionLibrary ("IOFunctions.bf");

/**
 * @name utility.AssociativeListToJSON
 * @param associative_list
 */
function utility.AssociativeListToJSON(associative_list) {

    // Replace inf and nan with 1+e9999 and null, respectively

    // SW20150122 TODO: Add recursion
    keys = Rows(associative_list);
    for (_k = 0; _k < Columns(keys); _k = _k+1) {
      current_val = associative_list[keys[_k]];
      if(Type(current_val) == "String") {
			  current_val = current_val^{{"inf"}{"1e999"}}; /* search and replace */
			  current_val = current_val^{{"-nan"}{"null"}}; /* search and replace */
        associative_list[keys[_k]] = current_val;
      }
    }

    return associative_list;
}

/**
 * @name utility.CallFunction
 * @param id
 * @param arguments
 */
function utility.CallFunction (id, arguments) {

    if (Type (id) == "String") {
        if (Type (arguments) == "AssociativeList") {
            return Eval ("`id` (" + Join (",", arguments) + ")");
        }
        return Eval ("`id`()");
    }
    return None;
}

function utility.convertToArgumentString (argument) {
    if (Type (argument) == "String") {
        return '"' + (argument && 2) + '"';
    }
    return argument;
}

/**
 * @name utility.Array1D
 * @param m
 */
lfunction utility.Array1D (m) {
    if (Type (m) == "Matrix") {
        return Rows (m) * Columns (m);
    } else {
        if (Type (m) == "AssociativeList") {
            return Abs (m);
        }
    }
    return None;
}

/**
 * @name utility.IsFunction
 * @param id
 */
function utility.IsFunction (id) {
	if (Type (id) == "String" && Abs (id) > 0) {
		ExecuteCommands ("GetString (__funcInfo, `id`, -1)");
		if (Type (__funcInfo) == "AssociativeList") {
			return Type (__funcInfo["ID"]) == "String" && Type (__funcInfo["Arguments"]) == "Matrix" && Type (__funcInfo["Body"]) == "String";
		}
	}
	return 0;
}

function utility.getGlobalValue (val) {
    return ^val;
}



/**
 * @name utility.ToggleEnvVariable
 * @param var
 * @param value
 */
function utility.ToggleEnvVariable (var, value) {
	if (None != value) {
	    if (Type (utility.ToggleEnvVariable.cache) != "AssociativeList") {
	        utility.ToggleEnvVariable.cache = {};
	    }
		utility.ToggleEnvVariable.cache[var] = Eval (var);
		*var = value;
	} else {
		*var = utility.ToggleEnvVariable.cache[var];
	}
}

/**
 * @name utility.GetEnvVariable
 * @param var
 */
function utility.GetEnvVariable (var) {
	return Eval(var);
}

/**
 * @name utility.SetEnvVariable
 * @param var
 */
function utility.SetEnvVariable (var, value) {
    Eval (var); // this is hack to make sure the variable exists before assigning to it
	^var = value;
}

/**
 * @name utility.CheckCacheFile
 * @param data_info
 */
lfunction utility.CheckCacheFile (data_info) {
    cache_info = {};
    cache_info["file"] = data_info["file"] + ".hyphy_cache";
    if (!(cache_info["file"])) {
        fscanf (cache_info["file"], "Raw", _cache);
        cache_info["cache"] = Eval (_cache);
    } else {
         cache_info["cache"] = {};
    }
    return cache_info;
}

/**
 * Find an element in the array.
 * @name utility.Find
 * @param {Matrix} array
 * @param {String|Number} value
 * @returns -1 if value is not found, otherwise returns the position
 */
lfunction utility.Find (array, value) {
    d = Rows (array) * Columns (array);
    for (i = 0; i < d; i+=1) {
        if (array [i] == value) {
            return i;
        }
    }
    return -1;
}

/**
 * @name utility.SwapKeysAndValues
 * @param dict
 */
lfunction utility.SwapKeysAndValues (dict) {
    swapped_dict = {};
    keys         = Rows (dict);
    items        = Abs (dict);

    for (i = 0; i < items; i+=1) {
        key = keys[i];
        swapped_dict [dict[key]] = key;
    }

    return swapped_dict;
}

/**
 * @name utility.ArrayToDict
 * @param object
 */
lfunction utility.ArrayToDict (object) {
    result = {};
    utility.ForEach(object, "_value_", "(`&result`)[_value_['key']] = _value_['value']");
    return result;
}

/**
 * Converts a matrix into a lookup dict, of the form value -> index
 * @name   utility.MatrixToDict
 * @param  {Matrix} object
 * @return {Dictionary} lookup table
 * @example utility.MatrixToDict ({{"a","b","c"}}) => {"a" : 0, "b" : 1, "c" : 2}
 */

lfunction utility.MatrixToDict (matrix) {
    result = {};
    counter = 0;
    utility.ForEach(matrix, "_value_", "(`&result`)[_value_] = `&counter`; `&counter` += 1;");
    return result;
}

/**
 * @name utility.DictToArray
 * @param object
 */
lfunction utility.DictToArray (object) {
    result = {1,Abs (object)};
    utility.ForEachPair(object, "key", "_value_", "(`&result`)[+key] = _value_");
    return result;
}

/**
 * @name utility.Map
 * @param {AssociativeList|Matrix} object - object to iterate over
 * @param {String} lambda_name - variable name for transform
 * @param {String} transform - function transform
 */
function utility.Map (object, lambda_name, transform) {

    Eval ("`lambda_name` = None");

    if (Type (object) == "AssociativeList") {
        utility.Map.return_object = {};
        utility.Map.keys = Rows (object);
        ^(lambda_name) := object [utility.Map.keys[utility.Map.k]];
        for (utility.Map.k = 0; utility.Map.k < Abs (object); utility.Map.k += 1) {
            utility.Map.return_object [utility.Map.keys[utility.Map.k]] = Eval (transform);
        }
        return utility.Map.return_object;
    }

    if (Type (object) == "Matrix") {
        utility.Map.rows = Rows (object);
        utility.Map.columns = Columns (object);
        utility.Map.return_object = {utility.Map.rows,  utility.Map.columns};

        ^(lambda_name) := object [utility.Map.r][utility.Map.c];
        for (utility.Map.r = 0; utility.Map.r < utility.Map.rows; utility.Map.r += 1) {
            for (utility.Map.c = 0; utility.Map.c < utility.Map.columns; utility.Map.c += 1) {
                utility.Map.temp = Eval (transform);
                assert (Type (utility.Map.temp) == "Number" || Type (utility.Map.temp) == "String", "Unsupported object type in call to utility.Map [Matrix]");
                utility.Map.return_object [utility.Map.r][utility.Map.c] = utility.Map.temp;

            }
        }
        return utility.Map.return_object;
    }

    return None;
}

/**
 * @name utility.MatrixToListOfRows
 * @param {Matrix} object - MxN matrix to convert
 * @param {Matrix} converted 1 x (M*N) Row vector
 */
lfunction utility.MatrixToListOfRows (object) {
    result = {};
    rows = Rows (object);
    cols = Columns (object);

    for (r = 0; r < rows; r += 1) {
        row = {1, cols};
        for (c = 0; c < cols; c+=1) {
            row[c] = object[r][c];
        }
        result + row;
    }
    return result;
}

/**
 * Filters 1xN Matrix or Dictionary
 * @name utility.Filter
 * @param {Dictionary|Matrix} object - matrix to convert
 * @param {String} lambda_name - function to discern whether element is filtered. All elements of iterable object that are false will be removed.
 * @returns filtered object
 * @example
 * _nonnegatives = utility.Filter (_data_vector, "_value_", "_value_ >= 0");
 */
function utility.Filter (object, lambda_name, condition) {

    Eval ("`lambda_name` = None");

    utility.Filter.return_object = {};
    if (Type (object) == "AssociativeList") {
        utility.Filter.keys = Rows (object);
        ^(lambda_name) := object [utility.Filter.keys[utility.Filter.k]];
        for (utility.Filter.k = 0; utility.Filter.k < Abs (object); utility.Filter.k += 1) {
            if (Eval (condition)) {
                utility.Filter.return_object [utility.Filter.keys[utility.Filter.k]] = ^(lambda_name);
            }
        }
        return utility.Filter.return_object;
    }

    if (Type (object) == "Matrix") {
        utility.Filter.rows = Rows (object);
        utility.Filter.columns = Columns (object);
        ^(lambda_name) := object [utility.Filter.r][utility.Filter.c];
        for (utility.Filter.r = 0; utility.Filter.r < utility.Filter.rows; utility.Filter.r += 1) {
            for (utility.Filter.c = 0; utility.Filter.c < utility.Filter.columns; utility.Filter.c += 1) {
                if (Eval (condition)) {
                    utility.Filter.return_object + ^(lambda_name);
                }
            }
        }
        return utility.Filter.return_object;
    }

    return None;
}

/**
 * Find the first element in a matrix or an array that matches the predicate; none otherwise
 * Keys are traversed alphabetically; matrices, by column - then row
 * @name utility.First
 * @param {Dictionary|Matrix} object - matrix to convert
 * @param {String} lambda_name - function to discern whether element is filtered.
 * @returns first matched object or none
 */
function utility.First (object, lambda_name, condition) {


    Eval ("`lambda_name` = None");

    utility.Filter.return_object = {};
    if (Type (object) == "AssociativeList") {
        utility.Filter.keys = Rows (object);
        ^(lambda_name) := object [utility.Filter.keys[utility.Filter.k]];
        for (utility.Filter.k = 0; utility.Filter.k < Abs (object); utility.Filter.k += 1) {


            if (Eval (condition)) {
                return ^(lambda_name);
            }
        }
        return None;
    }

    if (Type (object) == "Matrix") {
        utility.Filter.rows = Rows (object);
        utility.Filter.columns = Columns (object);
        ^(lambda_name) := object [utility.Filter.r][utility.Filter.c];
        for (utility.Filter.r = 0; utility.Filter.r < utility.Filter.rows; utility.Filter.r += 1) {
            for (utility.Filter.c = 0; utility.Filter.c < utility.Filter.columns; utility.Filter.c += 1) {
                if (Eval (condition)) {
                    return ^(lambda_name);
                }
            }
        }
        return None;
    }

    return None;
}

/**
 * @name utility.ForEach
 * @param {Tree|Dictionary|Matrix} object - matrix to convert
 * @param {String} lambda_name
 * @param {String} transform
 * @returns nothing
 */
function utility.ForEach (object, lambda_name, transform) {

    if (Type (object) == "Tree" || Type (object) == "Topology") {
        utility.ForEach (BranchName (object, -1), lambda_name, transform);
        return;
    }

    Eval ("`lambda_name` = None");

    if (Type (object) == "AssociativeList") {
        utility.ForEach.keys = Rows (object);
        ^(lambda_name) := object [utility.ForEach.keys[utility.ForEach.k]];
        for (utility.ForEach.k = 0; utility.ForEach.k < Abs (object); utility.ForEach.k += 1) {
            ExecuteCommands (transform);
        }
        return;
    }

    if (Type (object) == "Matrix") {
        utility.ForEach.rows = Rows (object);
        utility.ForEach.columns = Columns (object);
        utility.ForEach.return_object = {utility.ForEach.rows,  utility.ForEach.columns};
        ^(lambda_name) := object [utility.ForEach.r][utility.ForEach.c];
        for (utility.ForEach.r = 0; utility.ForEach.r < utility.ForEach.rows; utility.ForEach.r += 1) {
            for (utility.ForEach.c = 0; utility.ForEach.c < utility.ForEach.columns; utility.ForEach.c += 1) {
                ExecuteCommands (transform);

            }
        }
    }
}

/**
 * Checks whether key exists in dictionary
 * @name utility.KeyExists
 * @param {Dictionary} dict - dictionary to check
 * @param {String} key - key to check for existence
 * @returns {Bool} TRUE if key exists and is of expected type, otherwise FALSE
 */
function utility.KeyExists(dict, key) {

    keys = utility.Keys(dict);

    if(utility.Find(keys, key) != -1) {
        return TRUE;
    } else {
        return FALSE;
    }
}

/**
 * Checks whether key is a certain type
 * @name utility.CheckKey
 * @param {Dictionary} dict - dictionary to check
 * @param {String} key - key to check
 * @param {String} type - check whether key is "Matrix", "AssociativeList", "String", or "Tree"
 * @returns {Bool} TRUE if key exists and is of expected type, otherwise FALSE
 */
function utility.CheckKey (dict, key, type) {
    if (None != dict) {
        if (Type (dict[key]) == type) {
            return TRUE;
        }
    }
    return FALSE;
}

/**
 * Adds string or list of strings to a dictionary and sets value to 1
 * @name utility.AddToSet
 * @param {AssociativeList} set -
 * @param {String}, {Matrix}, or {AssociativeList} object
 * @returns nothing
 */
function utility.AddToSet (set, object) {

    if (Type(object) == "String") {
        set[object] = 1;
        return;
    }

    if (Type(object) == "AssociativeList" || Type(object) == "Matrix") {
        utility.ForEach (object, "_utility.AddToSet.value_", "set[_utility.AddToSet.value_] = 1");
        return;
    }

    set ["" + object] = 1;


}

/**
 * Set intersection
 * @name utility.Intersect
 * @param {AssociativeList} set - associative list to hold intersection
 * @param {AssociativeList} set1 - First set to intersect
 * @param {AssociativeList} set2 - Second set to intersect
 * @returns nothing
 */
function utility.Intersect(set, set1, set2) {

    keys1 = utility.Keys(set1);
    keys2 = utility.Keys(set2);

    if(Abs(Columns(keys1)) >  Abs(Columns(keys2))) {
        tmp_keys = keys2;
        keys2 = keys1;
        keys1 = tmp_keys;
    }

    for(k=0; k<Abs(Columns(keys1)); k+=1) {
        if(utility.Find(keys2, keys1[k]) != -1) {
            item = keys1[k];
            set["" + item] = 1;
        }
    }

}


/**
 * @name utility.PopulateDict
 * @param {Number} from
 * @param {Number} to
 * @param {Number|AssociativeList|String|Matrix} value
 * @param {String} lambda
 * @returns nothing
 */
function utility.PopulateDict (from, to, value, lambda) {
    utility.PopulateDict.result = {};
    if (Type (lambda) == "String" && Type (value) == "String") {
        Eval ("`lambda` = None");
        ^lambda := utility.PopulateDict.k;
        for (utility.PopulateDict.k = from; utility.PopulateDict.k < to; utility.PopulateDict.k+=1) {
            utility.PopulateDict.result[utility.PopulateDict.k] = Eval (value);
        }
   }
    else {
        for (utility.PopulateDict.k = from; utility.PopulateDict.k < to; utility.PopulateDict.k+=1) {
            utility.PopulateDict.result[utility.PopulateDict.k] = value;
        }
    }
    return utility.PopulateDict.result;
}

/**
 * Iterates over dictionaries
 * @name utility.ForEachPair
 * @param {Dictionary} object - the dictionary to iterate over
 * @param {String} key_name - the variable name for the key in the lambda expression
 * @param {String} value_name - the variable name for the value in the lambda expression
 * @param {String} transform - the lambda expression to use
 * @returns nothing
 * @example
 *    utility.ForEachPair (dict, "_key_", "_selection_",
 *        "fprintf(stdout, _key_);");
 */
function utility.ForEachPair(object, key_name, value_name, transform) {

    Eval ("`key_name` = None");
    Eval ("`value_name` = None");

    if (Type (object) == "AssociativeList") {
        utility.ForEachPair.keys = Rows (object);
        ^(key_name) := utility.ForEachPair.keys[utility.ForEachPair.k];
        ^(value_name) := object [utility.ForEachPair.keys[utility.ForEachPair.k]];
        for (utility.ForEachPair.k = 0; utility.ForEachPair.k < Abs (object); utility.ForEachPair.k += 1) {
            ExecuteCommands (transform);
        }
        return;
    }
    if (Type (object) == "Matrix") {
        utility.ForEachPair.rows = Rows (object);
        utility.ForEachPair.columns = Columns (object);
        utility.ForEachPair.return_object = {utility.ForEachPair.rows,  utility.ForEachPair.columns};
        ^(key_name) = {{utility.ForEachPair.r,utility.ForEachPair.c}};
        ^(value_name) := object [utility.ForEachPair.r][utility.ForEachPair.c];

        for (utility.ForEachPair.r = 0; utility.ForEachPair.r < utility.ForEachPair.rows; utility.ForEachPair.r += 1) {
            for (utility.ForEachPair.c = 0; utility.ForEachPair.c < utility.ForEachPair.columns; utility.ForEachPair.c += 1) {
                ExecuteCommands (transform);

            }
        }
    }

}

/**
 * Returns keys from a dictionary
 * @name utility.Keys
 * @param object - {Dictionary} the object to return keys from
 * @returns {Matrix} List of keys in dictionary
 * @example
 *   keys = utility.Keys(absrel.stats);
 *   keys == {"Count", "Mean", "Median", "Min", "Max", "2.5%", "97.5%", "Sum", "Std.Dev"};
 */
lfunction utility.Keys (object) {
    if (Type (object) == "AssociativeList") {
        return Rows (object);
    }
    return None;
}

/**
 * Returns values from a dictionary. Only returns unique values
 * @name utility.Values
 * @param object - {Dictionary} the object to return keys from
 * @returns {Matrix} List of keys in dictionary
 */
lfunction utility.Values (object) {
    if (Type (object) == "AssociativeList") {
        return Columns (object);
    }
    return None;
}

/**
 * Returns values from a dictionary. Only returns unique values
 * @name utility.UniqueValues
 * @param object - {Dictionary} the object to return keys from
 * @returns {Matrix} List of keys in dictionary
 */
lfunction utility.UniqueValues (object) {
    if (Type (object) == "AssociativeList") {
        return Columns (object);
    }
    if (Type (object) == "Matrix") {
        result = {};
        utility.ForEach (object, "_value_", "`&result`[_value_] += 1");
        return Rows (result);
    }
    return None;
}

/**
 * Ensures a key exists in a dictionary
 * @name utility.EnsureKey
 * @param {Dictionary} dict - the object to return keys from
 * @param {String} key - key to ensure exists
 * @returns nothing
 */
lfunction utility.EnsureKey (dict, key) {

    if (Type (dict[key]) != "AssociativeList") {
        dict[key] = {};
    }
}

/**
 * Ensures a key exists in a dictionary
 * @name utility.Has
 * @param {Dictionary} d - the object to return keys from
 * @param {String/Matrix} key - check to see if key is in the dictionary
 * if the value is a matrix, checks for [key[0]][key[1]][...][key[N]] nested keys
 * type check is performed on the last level
 * @param {String/None} type - if specified will further check if the object has the right type
 * @returns value mapped to the key or None
 */
lfunction utility.Has (d, key, type) {

    if (Type (key) == "String") {
        if (d/key) {
            if (type == None) {
                return TRUE;
            }
            return Type (d[key]) == type;
        }
        return FALSE;
    }

    if (Type (key) == "Matrix") {
        depth = utility.Array1D (key);
        current_check = &d;
        current_key   = key[0];

        for (i = 1; i < depth; i += 1) {
             if (Eval ("`current_check`/'`current_key`'")) {
                current_check = "(" + current_check + ")['`current_key`']";
                if (Eval ("Type(`current_check`)") != "AssociativeList") {
                    return FALSE;
                }
                current_key = key[i];
            } else {
                return FALSE;
            }
        }


        if ( Eval ("`current_check`/'`current_key`'") ) {
            if (type == None) {
                return TRUE;
            }
            return Eval ("Type((`current_check`)['`current_key`'])") == type;
        }
    }
    return FALSE;
}


/**
 * Returns the list of modules loaded with `LoadFunctionLibrary`
 * @returns a string matrix with (absolute) file paths for loaded modules
 */
lfunction utility.GetListOfLoadedModules () {
    GetString (res, LIST_OF_LOADED_LIBRARIES, -1);
    return res;
}


