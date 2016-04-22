LoadFunctionLibrary ("IOFunctions.bf");

/**
* utility.associativeListToJSON
* @param associative_list
*/
function utility.associativeListToJSON(associative_list) {

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
* utility.callFunction
* @param id
* @param arguments
*/
function utility.callFunction (id, arguments) {

    if (Type (id) == "String") {
        if (Type (arguments) == "AssociativeList") {
            return Eval ("`id` (" + Join (",", arguments) + ")");
        }
        return Eval ("`id`()");
    }
    return None;
}

/**
* utility.array1D
* @param m
*/
lfunction utility.array1D (m) {
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
* utility.isFunction
* @param id
*/
function utility.isFunction (id) {
	if (Type (id) == "String" && Abs (id) > 0) {
		ExecuteCommands ("GetString (__funcInfo, `id`, -1)");
		if (Type (__funcInfo) == "AssociativeList") {
			return Type (__funcInfo["ID"]) == "String" && Type (__funcInfo["Arguments"]) == "Matrix" && Type (__funcInfo["Body"]) == "String";
		}
	}
	return 0;
}

/**
* utility.toggleEnvVariable
* @param var
* @param value
*/
function utility.toggleEnvVariable (var, value) {
	if (None != value) {
		utilityFunction.toggleEnvVariable.__stash = Eval (var);
		*var = value;
	} else {
		*var = utilityFunction.toggleEnvVariable.__stash;
	}
}

/**
* utility.getEnvVariable
* @param var
*/
function utility.getEnvVariable (var) {
	return Eval(var);
}

/**
* utility.setEnvVariable
* @param var
*/
function utility.setEnvVariable (var, value) {
    Eval (var); // this is hack to make sure the variable exists before assigning to it
	^var = value;
}

/**
* utility.checkCacheFile
* @param data_info
*/
lfunction utility.checkCacheFile (data_info) {
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
* utility.array.find
* @param array
* @param value
*/
lfunction utility.array.find (array, value) {
    d = Rows (array) * Columns (array);
    for (i = 0; i < d; i+=1) {
        if (array [i] == value) {
            return i;
        }
    }
    return -1;
}

/**
* utility.dict.swap_keys_and_values
* @param dict
*/
lfunction utility.dict.swap_keys_and_values (dict) {
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
* utility.array_to_dict
* @param object
*/
lfunction utility.array_to_dict (object) {
    result = {};
    utility.forEach(object, "_value_", "(`&result`)[_value_['key']] = _value_['value']");
    return result;
}

// Only returns unique values
/**
* utility.dict_to_array
* @param object
*/
lfunction utility.dict_to_array (object) {
    result = {1,Abs (object)};
    utility.forEachPair(object, "key", "_value_", "(`&result`)[+key] = _value_");
    return result;
}

/**
* utility.map
* @param object
* @param lambda_name
* @param transform
*/
function utility.map (object, lambda_name, transform) {

    Eval ("`lambda_name` = None");

    if (Type (object) == "AssociativeList") {
        utility.map.return_object = {};
        utility.map.keys = Rows (object);
        ^(lambda_name) := object [utility.map.keys[utility.map.k]];
        for (utility.map.k = 0; utility.map.k < Abs (object); utility.map.k += 1) {
            utility.map.return_object [utility.map.keys[utility.map.k]] = Eval (transform);
        }
        return utility.map.return_object;
    }

    if (Type (object) == "Matrix") {
        utility.map.rows = Rows (object);
        utility.map.columns = Columns (object);
        utility.map.return_object = {utility.map.rows,  utility.map.columns};

        ^(lambda_name) := object [utility.map.r][utility.map.c];
        for (utility.map.r = 0; utility.map.r < utility.map.rows; utility.map.r += 1) {
            for (utility.map.c = 0; utility.map.c < utility.map.columns; utility.map.c += 1) {
                utility.map.temp = Eval (transform);
                assert (Type (utility.map.temp) == "Number" || Type (utility.map.temp) == "String", "Unsupported object type in call to utility.map [Matrix]");
                utility.map.return_object [utility.map.r][utility.map.c] = utility.map.temp;

            }
        }
        return utility.map.return_object;
    }

    return None;
}

lfunction utility.matrix_to_list_of_rows (object) {
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

function utility.filter (object, lambda_name, condition) {

    Eval ("`lambda_name` = None");

    utility.filter.return_object = {};
    if (Type (object) == "AssociativeList") {
        utility.filter.keys = Rows (object);
        ^(lambda_name) := object [utility.filter.keys[utility.filter.k]];
        for (utility.filter.k = 0; utility.filter.k < Abs (object); utility.filter.k += 1) {
            if (Eval (condition)) {
                utility.filter.return_object [utility.filter.keys[utility.filter.k]] = ^(lambda_name);
            }
        }
        return utility.filter.return_object;
    }

    if (Type (object) == "Matrix") {
        utility.filter.rows = Rows (object);
        utility.filter.columns = Columns (object);
        ^(lambda_name) := object [utility.filter.r][utility.filter.c];
        for (utility.filter.r = 0; utility.filter.r < utility.filter.rows; utility.filter.r += 1) {
            for (utility.filter.c = 0; utility.filter.c < utility.filter.columns; utility.filter.c += 1) {
                if (Eval (condition)) {
                    utility.filter.return_object + ^(lambda_name);
                }
            }
        }
        return utility.filter.return_object;
    }

    return None;
}

function utility.forEach (object, lambda_name, transform) {

    if (Type (object) == "Tree" || Type (object) == "Topology") {
        utility.forEach (BranchName (object, -1), lambda_name, transform);
        return;
    }

    Eval ("`lambda_name` = None");

    if (Type (object) == "AssociativeList") {
        utility.forEach.keys = Rows (object);
        ^(lambda_name) := object [utility.forEach.keys[utility.forEach.k]];
        for (utility.forEach.k = 0; utility.forEach.k < Abs (object); utility.forEach.k += 1) {
            ExecuteCommands (transform);
        }
        return;
    }

    if (Type (object) == "Matrix") {
        utility.forEach.rows = Rows (object);
        utility.forEach.columns = Columns (object);
        utility.forEach.return_object = {utility.forEach.rows,  utility.forEach.columns};
        ^(lambda_name) := object [utility.forEach.r][utility.forEach.c];
        for (utility.forEach.r = 0; utility.forEach.r < utility.forEach.rows; utility.forEach.r += 1) {
            for (utility.forEach.c = 0; utility.forEach.c < utility.forEach.columns; utility.forEach.c += 1) {
                ExecuteCommands (transform);

            }
        }
    }
}

function utility.checkKey (dict, key, type) {
    if (None != dict) {
        if (Type (dict[key]) == type) {
            return TRUE;
        }
    }
    return FALSE;
}

function utility.addToSet (set, object) {

    if (Type(object) == "String") {
        set[object] = 1;
        return;
    }
    if (Type(object) == "AssociativeList" || Type(object) == "Matrix") {
        utility.forEach (object, "_utility.addToSet.value_", "set[_utility.addToSet.value_] = 1");
        return;
    }
    set ["" + object] = 1;
}

function utility.populateDict (from, to, value, lambda) {
    utility.populateDict.result = {};
    if (Type (lambda) == "String" && Type (value) == "String") {
        Eval ("`lambda` = None");
        ^lambda = utility.populateDict.k;
        for (utility.populateDict.k = from; utility.populateDict.k < to; utility.populateDict.k+=1) {
            utility.populateDict.result[utility.populateDict.k] = Eval (value);
        }
   }
    else {
        for (utility.populateDict.k = from; utility.populateDict.k < to; utility.populateDict.k+=1) {
            utility.populateDict.result[utility.populateDict.k] = value;
        }
    }
    return utility.populateDict.result;
}

function utility.forEachPair (object, key_name, value_name, transform) {

    Eval ("`key_name` = None");
    Eval ("`value_name` = None");



    if (Type (object) == "AssociativeList") {
        utility.forEachPair.keys = Rows (object);
        ^(key_name) := utility.forEachPair.keys[utility.forEachPair.k];
        ^(value_name) := object [utility.forEachPair.keys[utility.forEachPair.k]];
        for (utility.forEachPair.k = 0; utility.forEachPair.k < Abs (object); utility.forEachPair.k += 1) {
            ExecuteCommands (transform);
        }
        return;
    }
    if (Type (object) == "Matrix") {
        utility.forEachPair.rows = Rows (object);
        utility.forEachPair.columns = Columns (object);
        utility.forEachPair.return_object = {utility.forEachPair.rows,  utility.forEachPair.columns};
        ^(key_name) = {{utility.forEachPair.r,utility.forEachPair.c}};
        ^(value_name) := object [utility.forEachPair.r][utility.forEachPair.c];

        for (utility.forEachPair.r = 0; utility.forEachPair.r < utility.forEachPair.rows; utility.forEachPair.r += 1) {
            for (utility.forEachPair.c = 0; utility.forEachPair.c < utility.forEachPair.columns; utility.forEachPair.c += 1) {
                ExecuteCommands (transform);

            }
        }
    }

}

lfunction utility.keys (object) {
    if (Type (object) == "AssociativeList") {
        return Rows (object);
    }
    return None;
}

// Only returns unique values
lfunction utility.values (object) {
    if (Type (object) == "AssociativeList") {
        return Columns (object);
    }
    return None;
}

lfunction utility.dict.ensure_key (dict, key) {

    if (Type (dict[key]) != "AssociativeList") {
        dict[key] = {};
    }
}


