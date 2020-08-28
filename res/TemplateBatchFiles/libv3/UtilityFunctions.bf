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
    return 0;
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
	return FALSE;
}

function utility.getGlobalValue (val) {
    return ^val;
}

function utility.setGlobalValue (id, val) {
    Eval (id);
    ^id = val;
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
    Eval (var); // this is a hack to make sure the variable exists before assigning to it
	^var = value;
}

/**
 * @name utility.CheckCacheFile
 * @param data_info
 */
lfunction utility.CheckCacheFile (data_info) {
    cache_info = {};
    cache_info[utility.getGlobalValue("terms.data.file")] = data_info[utility.getGlobalValue("terms.data.file")] + ".hyphy_cache";
    if (!(cache_info[utility.getGlobalValue("terms.data.file")])) {
        fscanf (cache_info[utility.getGlobalValue("terms.data.file")], "Raw", _cache);
        cache_info[utility.getGlobalValue("terms.data.cache")] = Eval (_cache);
    } else {
         cache_info[utility.getGlobalValue("terms.data.cache")] = {};
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
    for (i, v; in; array) {
        if (v == value) {
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
 * @name utility.DictToSortedArray
 * @param object
 */
lfunction utility.DictToSortedArray (object) {
    result = {Abs (object),2};
    utility.ForEachPair(object, "_key_", "_value_", "(`&result`)[+_key_][0] = _value_;(`&result`)[+_key_][1] = +_key_;");
    return result % 0;
}


/**
 * @name utility.Range
 * @param elements
 * @param from
 * @param step
 */
lfunction utility.Range (elements, from, step) {
    range = {};
    for (i = 0; i < elements; i+=1) {
        range + (from + i * step);
    }
    return range;
}

/**
 * @name utility.Map
 * @param {AssociativeList|Matrix} object - object to iterate over
 * @param {String} lambda_name - variable name for transform
 * @param {String} transform - function transform
 */
function utility.Map (object, lambda_name, transform) {

    Eval ("`lambda_name` = None");
    utility.Map.return_object = None;

    if (Type (object) == "AssociativeList") {
        ExecuteCommands ("function  utility.Map.CB (`lambda_name`) {return `transform`;}", enclosing_namespace);
        utility.Map.return_object = {};
        utility.Map.keys = Rows (object);
        for (utility.Map.k = 0; utility.Map.k < Abs (object); utility.Map.k += 1) {
            //^(lambda_name) = object [utility.Map.keys[utility.Map.k]];
            //utility.Map.return_object [utility.Map.keys[utility.Map.k]] = Eval (transform);
            utility.Map.return_object [utility.Map.keys[utility.Map.k]] = Call ("utility.Map.CB", object [utility.Map.keys[utility.Map.k]]);
        }
    } else {
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

        }
    }
    Eval ("`lambda_name` = None");

    return utility.Map.return_object;
}

/**
 * @name utility.MapWithKey
 * @param {AssociativeList|Matrix} object - object to iterate over
 * @param {String} lambda_name - variable name for transform
 * @param {String} transform - function transform
 */
function utility.MapWithKey (object, key_name, lambda_name, transform) {

    Eval ("`lambda_name` = None");
    Eval ("`key_name` = None");
    utility.MapWithKey.return_object = None;
    

    if (Type (object) == "AssociativeList") {
        ExecuteCommands ("function  utility.MapWithKey.CB (`key_name`, `lambda_name`) {return `transform`;}", enclosing_namespace);
        utility.MapWithKey.return_object = {};
        utility.MapWithKey.keys = Rows (object);
        for (utility.MapWithKey.k = 0; utility.MapWithKey.k < Abs (object); utility.MapWithKey.k += 1) {
            utility.MapWithKey.key = utility.MapWithKey.keys[utility.MapWithKey.k];
            utility.MapWithKey.return_object [utility.MapWithKey.key] = Call ("utility.MapWithKey.CB", utility.MapWithKey.key, object[utility.MapWithKey.key]);
        }
    } else {
        if (Type (object) == "Matrix") {
            utility.MapWithKey.rows = Rows (object);
            utility.MapWithKey.columns = Columns (object);
            utility.MapWithKey.return_object = {utility.MapWithKey.rows,  utility.MapWithKey.columns};

            ^(lambda_name) := object [utility.MapWithKey.r][utility.MapWithKey.c];
            ^(key_name) := {{utility.MapWithKey.r,utility.MapWithKey.c}};
            for (utility.MapWithKey.r = 0; utility.MapWithKey.r < utility.MapWithKey.rows; utility.MapWithKey.r += 1) {
                for (utility.MapWithKey.c = 0; utility.MapWithKey.c < utility.MapWithKey.columns; utility.MapWithKey.c += 1) {
                    utility.MapWithKey.temp = Eval (transform);
                    assert (Type (utility.MapWithKey.temp) == "Number" || Type (utility.MapWithKey.temp) == "String", "Unsupported object type in call to utility.MapWithKey [Matrix]");
                    utility.MapWithKey.return_object [utility.MapWithKey.r][utility.MapWithKey.c] = utility.MapWithKey.temp;

                }
            }

        }
    }
    Eval ("`lambda_name` = None");
    Eval ("`key_name` = None");

    return utility.MapWithKey.return_object;
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

    utility.Filter.return_object = None;

    Eval ("`lambda_name` = None");

    utility.Filter.return_object = {};
    if (Type (object) == "AssociativeList") {
        utility.Filter.keys = Rows (object);
        for (utility.Filter.k = 0; utility.Filter.k < Abs (object); utility.Filter.k += 1) {
           ^(lambda_name) = object [utility.Filter.keys[utility.Filter.k]];
           if (Eval (condition)) {
                utility.Filter.return_object [utility.Filter.keys[utility.Filter.k]] = ^(lambda_name);
            }
        }
    } else {
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
        }
    }
    Eval ("`lambda_name` = None");
    return utility.Filter.return_object;
}

/**
 * Find the first element in a matrix or an array that matches the predicate; none otherwise
 * Keys are traversed alphabetically; matrices, by column - then row
 * @name utility.First
 * @param {Dictionary|Matrix} object - matrix to convert
 * @param {String} lambda_name - function to discern whether element is filtered.
 * @returns first matched object or none
 */
function utility.First (object_utility_first, lambda_name, condition) {

     utility.First.return_object = None;

     Eval ("`lambda_name` = None");

     if (Type (object_utility_first) == "AssociativeList") {
        utility.First.keys = Rows (object_utility_first);
        for (utility.First.k = 0; utility.First.k < Abs (object_utility_first); utility.First.k += 1) {
            ^(lambda_name) = object_utility_first [utility.First.keys[utility.First.k]];
            if (Eval (condition)) {
                utility.First.return_object = ^(lambda_name);
                break;
            }
        }
    }

    if (Type (object_utility_first) == "Matrix") {
        utility.First.rows = Rows (object_utility_first);
        utility.First.columns = Columns (object_utility_first);
        ^(lambda_name) := object_utility_first [utility.First.r][utility.First.c];

        for (utility.First.r = 0; utility.First.r < utility.First.rows; utility.First.r += 1) {
            for (utility.First.c = 0; utility.First.c < utility.First.columns; utility.First.c += 1) {
                if (Eval (condition)) {
                     utility.First.return_object = ^(lambda_name);
                     break;
                }
            }
        }
    }

    Eval ("`lambda_name` = None");

    return utility.First.return_object;
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
    	// strip out the root node
    	utility.ForEach._aux2 = BranchName (object, -1);
    	utility.ForEach.rows = utility.Array1D (utility.ForEach._aux2) - 1;
    	utility.ForEach._aux = {utility.ForEach.rows, 1};
    	for (utility.ForEach.r = 0; utility.ForEach.r < utility.ForEach.rows; utility.ForEach.r += 1) {
    		utility.ForEach._aux [utility.ForEach.r ] = utility.ForEach._aux2 [utility.ForEach.r ];
    	}

        DeleteObject (utility.ForEach._aux2);
        utility.ForEach (utility.ForEach._aux, lambda_name, transform);
        DeleteObject (utility.ForEach._aux);
        return;
    }

    Eval ("`lambda_name` = None");

    if (Type (object) == "AssociativeList") {


        utility.ForEach.keys = Rows (object);
        utility.ForEach.size = Abs (object);


        for (utility.ForEach.k = 0; utility.ForEach.k < utility.ForEach.size; utility.ForEach.k += 1) {
          ^(lambda_name) = object [utility.ForEach.keys[utility.ForEach.k]];
           ExecuteCommands (transform, enclosing_namespace);
        }
    } else {

        if (Type (object) == "Matrix") {
            utility.ForEach.rows = Rows (object);
            utility.ForEach.columns = Columns (object);

            if (utility.ForEach.rows && utility.ForEach.columns)  {
                 ^(lambda_name) := object [utility.ForEach.r][utility.ForEach.c];
                for (utility.ForEach.r = 0; utility.ForEach.r < utility.ForEach.rows; utility.ForEach.r += 1) {
                    for (utility.ForEach.c = 0; utility.ForEach.c < utility.ForEach.columns; utility.ForEach.c += 1) {
                        ExecuteCommands (transform, enclosing_namespace);
                    }
                }
            }
        }
    }

    Eval ("`lambda_name` = None");
}

/**
 * Checks whether key exists in dictionary
 * @name utility.KeyExists
 * @param {Dictionary} dict - dictionary to check
 * @param {String} key - key to check for existence
 * @returns {Bool} TRUE if key exists and is of expected type, otherwise FALSE
 */
function utility.KeyExists(dict, key) {
    return dict / key;
}

/**
 * Given the lookup string (lookup), maps each character in
 * string to the index of the same character in lookup (>=0)
 * or -1 if no such character exists
 * @name utility.MapStrings
 * @param {String} string - the string to map
 * @param {String} lookup - mapping lookup string
 * @returns {Matrix} the mapping
 * @example utility.MapStrings ("ABCD","DOBAZ") => {
 *        "0":3,
 *        "1":2,
 *        "2":-1,
 *        "3":0
 *       }
 */
lfunction utility.MapStrings (string,lookup) {

// source ID -> target ID (-1 means no correspondence)

	mapping 	  = {};
	targetIndexing = {};
	_d = Abs(lookup);

	for (_i = 0; _i < _d; _i += 1) {
		targetIndexing [lookup[_i]] = _i + 1;
	}

	_d = Abs (string);
	for (_i = 0; _i < _d; _i += 1) {
		mapping [_i] = targetIndexing[string[_i]] - 1;
	}

	return mapping;
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
        Eval ("`lambda` = None");
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

    io.CheckAssertion ("!utility.ForEachPair.warn_non_reentrant", "utility.ForEachPair is non re-entrant");
    utility.ForEachPair.warn_non_reentrant = TRUE;

    ExecuteCommands ("function  utility.ForEachPair.CB (`key_name`, `value_name`) {`transform`}", enclosing_namespace);

    if (Type (object) == "AssociativeList") {
        utility.ForEachPair.keys = Rows (object);
         for (utility.ForEachPair.k = 0; utility.ForEachPair.k < Abs (object); utility.ForEachPair.k += 1) {            //ExecuteCommands (transform, enclosing_namespace);
            utility.ForEachPair.key = utility.ForEachPair.keys[utility.ForEachPair.k];
            Call ("utility.ForEachPair.CB",utility.ForEachPair.key , object [utility.ForEachPair.key]);
        }
    } else {
        if (Type (object) == "Matrix") {
            utility.ForEachPair.rows = Rows (object);
            utility.ForEachPair.columns = Columns (object);

            if (utility.ForEachPair.rows && utility.ForEachPair.columns) {

                //^(key_name) = {{utility.ForEachPair.r,utility.ForEachPair.c}};
                //^(value_name) := object [utility.ForEachPair.r][utility.ForEachPair.c];
                utility.ForEachPair.key = {{utility.ForEachPair.r,utility.ForEachPair.c}};
                for (utility.ForEachPair.r = 0; utility.ForEachPair.r < utility.ForEachPair.rows; utility.ForEachPair.r += 1) {
                    for (utility.ForEachPair.c = 0; utility.ForEachPair.c < utility.ForEachPair.columns; utility.ForEachPair.c += 1) {
                        //ExecuteCommands (transform, enclosing_namespace);
                        Call ("utility.ForEachPair.CB",utility.ForEachPair.key , object [utility.ForEachPair.r][utility.ForEachPair.c]);

                    }
                }
            }
        }
    }

    Eval ("`key_name` = None");
    Eval ("`value_name` = None");
    // reset bindings here to avoid calls on stale formula in subsequent invocations

    utility.ForEachPair.warn_non_reentrant = FALSE;

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
 * Returns values from a dictionary. Not just unique values and the values aren't coerced to strings.
 * @name utility.Values
 * @param object - {Dictionary} the object to return keys from
 * @returns {Matrix} List of keys in dictionary
 */
lfunction utility.Values (object) {
    if (Type (object) == "AssociativeList") {
        keys = utility.Keys(object);
        values = {1,Abs(object)};
        for(i=0; i<Abs(object); i=i+1) {
            values[i] = object[keys[i]];
        }
        return values;
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
 * @returns T/F
 */
lfunction utility.Has (d, key, type) {

    if (Type (d) == "AssociativeList") {
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
    }
    return FALSE;
}

/**
 * Ensures a key exists in a dictionary
 * @name utility.GetByKey
 * @param {Dictionary} d - the object to return keys from
 * @param {String/Matrix} key - check to see if key is in the dictionary
 * if the value is a matrix, checks for [key[0]][key[1]][...][key[N]] nested keys
 * type check is performed on the last level
 * @param {String/None} type - if specified will further check if the object has the right type
 * @returns value mapped to the key or None
 */
lfunction utility.GetByKey (d, key, type) {

    if (Type (d) == "AssociativeList") {
        if (Type (key) == "String") {
            if (d/key) {
                if (type == None) {
                    return d[key];
                }
                if (Type (d[key]) == type) {
                    return d[key];
                }
            }
            return None;
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
                return_value = Eval ("(`current_check`)['`current_key`']");
                if (type == None) {
                    return return_value;
                }
                if (Type (return_value) == type) {
                    return return_value;
                }
             }
        }
    }
    return None;
}

/**
 * Extends a dict with keys from another dict
 * @name utility.Extend
 * @param {Dictionary} d - the object to extend
 * @param {Dictionary} n - the object whose key/value pairs will be added to 'd'
 * @returns d with new keys added
 */
lfunction utility.Extend (d, n) {

    if (Type (d) == "AssociativeList" && Type (n) == "AssociativeList") {
        nkeys = utility.Keys (n);
        size = Abs (n);
        for (i = 0; i < size; i+=1) {
            if (d/nkeys[i] == FALSE) {
                d[nkeys[i]] = n[nkeys[i]];
            }
        }
    }

    return d;
}


/**
 * Returns the list of currently defined variables whose names match a regex
 * @param selector {String} : reg exp selector; show all if None
 * @return {MATRIX} : all the selected variables
 */
lfunction utility.DefinedVariablesMatchingRegex (selector) {
    if (Type (selector) == "String") {
        ExecuteCommands ('GetInformation (`&res`, "`selector`")');
    } else {
        GetInformation (res, ".");
    }
    return res;
}

/**
 * Bin an iterable-type object by value
 * @param  obj {Dict | Matrix}
 * @return Dict value->list of keys
 */
lfunction utility.BinByValue (obj) {
    if (Type (obj) == "AssociativeList") {
        keys  = Rows (obj);
        count = Abs  (obj);
        result = {};
        for (i = 0; i < count; i+=1) {
            str_value = "" + obj[keys[i]];
            if (result / str_value) {
                result [str_value] + keys[i];
            } else {
                result [str_value] = {"0" : keys[i]};
            }
        }
        return result;
    } else {
        if (Type (obj) == "Matrix") {
            count = Rows (obj)*Columns (obj);
            result = {};
            for (i = 0; i < count; i+=1) {
                str_value = "" + obj[i];
                if (result / str_value) {
                    result [str_value] + i;
                } else {
                    result [str_value] = {"0" : i};
                }
            }
        }
        return result;
    }

    return None;
}



/**
 * Returns the list of modules loaded with `LoadFunctionLibrary`
 * param: {String} filter, if provided, will filter the list of library paths via a regexp
 * @returns a string matrix with (absolute) file paths for loaded modules
 */
lfunction utility.GetListOfLoadedModules (filter) {
    GetString (res, LIST_OF_LOADED_LIBRARIES, -1);
    if (None != filter) {
        return utility.UniqueValues (utility.Filter (res, "_path_", "regexp.Find(_path_,`&filter`)"));
    }
    return res;
}

/**
 * Returns the list of likelihood function IDs that are currently available
 * param: {String} filter, if provided, will filter the list of library paths via a regexp
 * @returns a string matrix with (absolute) file paths for loaded modules
 */
lfunction utility.GetListOfLoadedLikelihoodFunctions (filter) {

    lf_count = Rows ("LikelihoodFunction");

    res = {lf_count, 1};

    for (k = 0; k < lf_count; k+=1) {
        GetString (lf_id, LikelihoodFunction, k);
        res[k] = lf_id;
    }

    if (None != filter) {
        return utility.UniqueValues (utility.Filter (res, "_path_", "regexp.Find(_path_,`&filter`)"));
    }
    return res;
}

/**
 * Execute some commands in the global namespace
 * @param commands {String} : the commands to execture
 * @returns None
 */
function utility.ExecuteInGlobalNamespace (commands) {
    ExecuteCommands (commands);
}

/**
 * A product function, which takes a list of k variable IDs (dicts)
 * and returns a dictionary of sets of keys (N1 x N2 x ... Nk)
 * @param {Dictionary} : a list of variable IDs to scan
 * @returns {Dictionary} : {"0" : {{"key_11","key_21",...,"key_k1"}}, "1" : {"key_11","key_21",...,"key_k2"}...}
 * or None if something went wrong
 */
lfunction utility.CatersianProduct (arguments) {
    if (Type (arguments) != "AssociativeList") {
        return None;
    }

    arg_count = Abs (arguments);

    if (arg_count == 0) {
        return None;
    }

    product = {};

    key_counts = {1,arg_count};
    arg_keys   = utility.Keys (arguments);
    key_sets   = {};

    for (k = 0; k < arg_count; k+=1) {
        this_key = arg_keys[k];
        if (Type (arguments[this_key]) != "AssociativeList") {
            return None;
        }
        key_sets   [k] = utility.Keys (arguments[this_key]);
        key_counts [k] = Abs (arguments[this_key]);
        if (key_counts[k] == 0) {
            return None;
        }
    }

    indices = {1,arg_count};
    done    = FALSE;

    while (1) {
        an_item = {};
        for (k = 0; k < arg_count; k+=1) {
            an_item + (key_sets[k])[indices[k]];
        }
        product + an_item;

        for (k = arg_count - 1; k >= 0; k = k - 1) {
            indices [k] += 1;
            if (indices[k] == key_counts[k]) {
                indices[k] = 0;
            } else {
                break;
            }
        }

        if (k < 0) {
            break;
        }

    }

    return product;
}

function _sortStringsAux (theKey, theValue) {
    
 	for (_gb_idx2 = 0; _gb_idx2 < theValue; _gb_idx2 += 1) {
		_gb_sortedStrings [_gb_idx] = theKey;
		_gb_idx += 1;
	}
	return 0;
}

function utility.sortStrings (_theList) {
	_gb_dim = Rows (_theList)*Columns (_theList);
	_toSort = {};
	for (_gb_idx = 0; _gb_idx < _gb_dim; _gb_idx = _gb_idx + 1) {
		_toSort[_theList[_gb_idx]] = _toSort[_theList[_gb_idx]]+1;
	}
	_gb_sortedStrings = {_gb_dim,1};
	_gb_idx = 0;
	_toSort["_sortStringsAux"][""];
	return _gb_sortedStrings;
}

/*
function utility.FinishAndPrintProfile () {
#profile _hyphy_profile_dump;

    stats  			= _hyphy_profile_dump["STATS"];
    _profile_summer = ({1,Rows(stats)}["1"]) * stats;
    _instructions   = _hyphy_profile_dump["INSTRUCTION"];
    _indices	    = _hyphy_profile_dump["INSTRUCTION INDEX"];

    fprintf (stdout, "\nTotal run time (seconds)      : ", Format(_profile_summer[1],15,6),
                     "\nTotal number of steps         : ", Format(_profile_summer[0],15,0), "\n\n");

    to_sort        =  stats["-_MATRIX_ELEMENT_VALUE_*_MATRIX_ELEMENT_COLUMN_+(_MATRIX_ELEMENT_COLUMN_==0)*_MATRIX_ELEMENT_ROW_"] % 1;

    for (k=0; k<Columns(_instructions); k=k+1)
    {
        k2 = to_sort[k][0];
        fprintf (stdout, Format (_indices[k2],6,0), " : ", _instructions[k2], "\n\tCall count: ", stats[k2][0],
                                                       "\n\tTime (seconds): ", stats[k2][1], "\n");
    }
    return 0;
}
*/
