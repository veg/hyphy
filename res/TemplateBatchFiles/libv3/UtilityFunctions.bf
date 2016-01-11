LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("IOFunctions.bf");

function utility.promptForGeneticCodeAndAlignment (dataset_name, datafilter_name) {
    return utility.loadCodonDataFile (dataset_name, io.readCodonDataSet (dataset_name));
}

function utility.loadGeneticCodeAndAlignment (dataset_name, datafilter_name, path) {
    return utility.loadCodonDataFile (dataset_name, io.readCodonDataSetFromPath (dataset_name, path));
}

function utility.loadCodonDataFile (dataset_name, data_info) {
    ExecuteCommands ("DataSetFilter `datafilter_name` = CreateFilter (`dataset_name`, 3, , , data_info[\"stop\"])");
    io.checkAssertion ("`datafilter_name`.sites*3==`dataset_name`.sites","The input alignment must not contain stop codons");
    data_info ["sites"] = Eval ("`datafilter_name`.sites");
    data_info ["dataset"] = dataset_name;
    data_info ["datafilter"] = datafilter_name;

    return data_info;
}

function utility.readNucleotideAlignment (file_name, dataset_name, datafilter_name) {
    data_info = io.readNucleotideDataSet (dataset_name, file_name);
    ExecuteCommands ("DataSetFilter `datafilter_name` = CreateFilter (`dataset_name`,1)");
    data_info ["sites"] = Eval ("`datafilter_name`.sites");
    data_info ["dataset"] = dataset_name;
    data_info ["datafilter"] = datafilter_name;

    return data_info;
}

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

function utility.defineFrequencies (datafilter_name) {
    HarvestFrequencies	          (nuc3, *datafilter_name, 3, 1, 1);
    nucCF						= CF3x4	(nuc3, GeneticCodeExclusions);
    codon3x4					= BuildCodonFrequencies (nucCF);
    return {"raw": nuc3,
            "nucleotide": nucCF,
            "codon": codon3x4};
}

lfunction utility.partition_tree (avl, l) {
    for (k = 0; k < Abs (avl); k+=1) {
        if ((avl[k])["Parent"]) {
            if (Abs ((avl[k])["Children"])) {
                l [(avl[k])["Name"]] = "internal";
            } else {
                l [(avl[k])["Name"]] = "leaf";
            }
        }
    }
}

function utility.loadAnnotatedTopology (look_for_newick_tree) {
    return utility.extractTreeInfo (io.getTreeString(look_for_newick_tree));
 }

function utility.extractTreeInfo (tree_string) {
    Topology     T = tree_string;

    utility.loadAnnotatedTopology.branch_lengths  = BranchLength (T, -1);
    utility.loadAnnotatedTopology.branch_names    = BranchName (T, -1);

    utility.loadAnnotatedTopology.bls            = {};

    for (utility.loadAnnotatedTopology.k = 0; utility.loadAnnotatedTopology.k < Columns (utility.loadAnnotatedTopology.branch_names) - 1; utility.loadAnnotatedTopology.k += 1) {
        if (utility.loadAnnotatedTopology.branch_lengths[utility.loadAnnotatedTopology.k] >= 0.) {
            utility.loadAnnotatedTopology.bls [utility.loadAnnotatedTopology.branch_names[utility.loadAnnotatedTopology.k]] =
               utility.loadAnnotatedTopology.branch_lengths[utility.loadAnnotatedTopology.k];
        }
    }

    GetInformation (modelMap, T);

    leaves_internals    = {};

    utility.partition_tree (T^0, leaves_internals);

    utility.toggleEnvVariable ("INCLUDE_MODEL_SPECS", 1);
    T.str = "" + T;
    utility.toggleEnvVariable ("INCLUDE_MODEL_SPECS", None);

    return {"string"     : Format (T,1,0),
            "branch_lengths" :  utility.loadAnnotatedTopology.bls,
            "annotated_string" : T.str ,
            "model_map"  : modelMap,
            "partitioned" : leaves_internals,
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

function utility.array1D (m) {
    return Rows (m) * Columns (m);
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
	if (None != value) {
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


