LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("IOFunctions.bf");

function utility.promptForGeneticCodeAndAlignment (dataset_name, datafilter_name) {
    data_info = io.readCodonDataSet (dataset_name);
    ExecuteCommands ("DataSetFilter `datafilter_name` = CreateFilter (`dataset_name`, 3, , , data_info[\"stop\"])");
    io.checkAssertion ("`datafilter_name`.sites*3==`dataset_name`.sites","The input alignment must not contain stop codons");
    data_info ["sites"] = Eval ("`datafilter_name`.sites");
    data_info ["dataset"] = dataset_name;
    data_info ["datafilter"] = datafilter_name;
    
    return data_info;
}

function utility.replace(search, replace, string) {
    string_search = "*" + search + "*";

    while (string/string_search) {
      inf_instances = string$search;
      new_string = "";

      for(x = 0; x < inf_instances[0]; x=x+1) {
        new_string = new_string + string[x];  
      } 

      new_string = new_string + replace;

      for(x = inf_instances[1] + 1; x < Abs(string); x=x+1) {
        new_string = new_string + string[x];  
      } 

      string = new_string;

    }

    return string;
}


function utility.associativeListToJSON(associative_list) {

    // Convert associative list to string
    string =  " " + associative_list + " ";
    string = utility.replace("inf,", "1e9999,", string);
    string = utility.replace("-nan,", "null,", string);
    string = utility.replace("inf]", "1e9999]", string);
    string = utility.replace("-nan]", "null]", string);

    // Replace inf and nan with 1+e9999 and null, respectively
    return Eval(string);
}

function utility.defineFrequencies (datafilter_name) {
    HarvestFrequencies	          (nuc3, *datafilter_name, 3, 1, 1);
    nucCF						= CF3x4	(nuc3, GeneticCodeExclusions);
    codon3x4					= BuildCodonFrequencies (nucCF);
    return {"raw": nuc3,
            "nucleotide": nucCF,
            "codon": codon3x4};
}

function utility.loadAnnotatedTopology (look_for_newick_tree) {
    tree_string = io.getTreeString(look_for_newick_tree);
    Topology     T = tree_string;
    GetInformation (modelMap, T);
    utility.toggleEnvVariable ("INCLUDE_MODEL_SPECS", 1);
    T.str = "" + T;
    utility.toggleEnvVariable ("INCLUDE_MODEL_SPECS", None);
    
    return {"string"     : Format (T,1,0),
            "annotated_string" : T.str ,
            "model_map"  : modelMap,
            "model_list" :  Columns (modelMap)};
}

function utility.callFunction (id, arguments) {
    if (Type (id) == "String") {
        if (Type (arguments) == "AssociativeList") {
            return Eval ("`id` (" + Join (",", arguments) + ")");	
        }
        return Eval ("`id`()");
    } 
    return None;
}

function utility.isFunction (id) {
	if (Type (id) == "String" && Abs (id) > 0) {
		ExecuteCommands ("GetString (__funcInfo, `id`, -1)");
		if (Type (__funcInfo) == "AssociativeList") {
			return Type (__funcInfo["ID"]) == "String" && Type (__funcInfo["Arguments"]) == "Matrix" && Type (__funcInfo["Body"]) == "String";
		}
	}
	return 0;
}

function utility.toggleEnvVariable (var, value) {
	if (value != None) {
		utilityFunction.toggleEnvVariable.__stash = Eval (var);
		*var = value;
	} else {
		*var = utilityFunction.toggleEnvVariable.__stash;
	}
}

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

lfunction utility.array.find (array, value) {
    d = Rows (array) * Columns (array);
    for (i = 0; i < d; i+=1) {
        if (array [i] == value) {
            return i;
        }
    }
    return -1;
}


