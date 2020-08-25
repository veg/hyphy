LoadFunctionLibrary("../IOFunctions.bf");
LoadFunctionLibrary("../all-terms.bf");
LoadFunctionLibrary("../convenience/regexp.bf");
LoadFunctionLibrary("../UtilityFunctions.bf");
LoadFunctionLibrary("TreeTools");
LoadFunctionLibrary("distances.bf");

/** @module trees */

/*
 * Sample tree object
 * trees = {
 *       "string":"((((PIG,COW)Node3,HORSE,CAT)Node2,((RHMONKEY,BABOON)Node9,(HUMAN,CHIMP)Node12)Node8)Node1,RAT,MOUSE)",
 *       "string_with_lengths":"((((PIG:0.147969,COW:0.21343)Node3:0.085099,HORSE:0.165787,CAT:0.264806)Node2:0.058611,((RHMONKEY:0.002015,BABOON:0.003108)Node9:0.022733,(HUMAN:0.004349,CHIMP:0.000799)Node12:0.011873)Node8:0.101856)Node1:0.340802,RAT:0.050958,MOUSE:0.09795)",
 *       "branch length":{
 *         "PIG":0.147969,
 *         "COW":0.21343,
 *         "Node3":0.08509899999999999,
 *         "HORSE":0.165787,
 *         "CAT":0.264806,
 *         "Node2":0.058611,
 *         "RHMONKEY":0.002015,
 *         "BABOON":0.003108,
 *         "Node9":0.022733,
 *         "HUMAN":0.004349,
 *         "CHIMP":0.000799,
 *         "Node12":0.011873,
 *         "Node8":0.101856,
 *         "Node1":0.340802,
 *         "RAT":0.050958,
 *         "MOUSE":0.09795
 *        },
 *       "annotated_string":"((((PIG,COW)Node3,HORSE,CAT)Node2,((RHMONKEY{PR},BABOON{PR})Node9{PR},(HUMAN{PR},CHIMP{PR})Node12{PR})Node8{PR})Node1{PR},RAT,MOUSE);",
 *       "model_map":{
 *         "PIG":"",
 *         "COW":"",
 *         "Node3":"",
 *         "HORSE":"",
 *         "CAT":"",
 *         "Node2":"",
 *         "RHMONKEY":"PR",
 *         "BABOON":"PR",
 *         "Node9":"PR",
 *         "HUMAN":"PR",
 *         "CHIMP":"PR",
 *         "Node12":"PR",
 *         "Node8":"PR",
 *         "Node1":"PR",
 *         "RAT":"",
 *         "MOUSE":""
 *        },
 *       "partitioned":{
 *         "PIG":"leaf",
 *         "COW":"leaf",
 *         "Node3":"internal",
 *         "HORSE":"leaf",
 *         "CAT":"leaf",
 *         "Node2":"internal",
 *         "RHMONKEY":"leaf",
 *         "BABOON":"leaf",
 *         "Node9":"internal",
 *         "HUMAN":"leaf",
 *         "CHIMP":"leaf",
 *         "Node12":"internal",
 *         "Node8":"internal",
 *         "Node1":"internal",
 *         "RAT":"leaf",
 *         "MOUSE":"leaf"
 *        },
 *       "model_list":{
 *        {"", "PR"}
 *        }
 *      }
 *
 */


/**
 * Returns sanitized Newick tree string
 * @name trees.GetTreeString._sanitize
 * @private
 * @param {String} string
 * @returns {String} sanitized string
 */
lfunction trees.GetTreeString._sanitize(string) {

    if (utility.GetEnvVariable("_DO_TREE_REBALANCE_")) {
        string = RerootTree(string, 0);
    }

    if (utility.GetEnvVariable("_KEEP_I_LABELS_")) {
        utility.ToggleEnvVariable("INTERNAL_NODE_PREFIX", "intNode");
    }
    /* this strips bootstrap values; not necessary for 2.4 */
    /* string = string ^ {
        {
            "\\)[0-9]+(\\.[0-9]*)?\:",
            "):"
        }
    }; */

    if (utility.GetEnvVariable("_KEEP_I_LABELS_")) {
        utility.ToggleEnvVariable("INTERNAL_NODE_PREFIX", None);
    }

    return string;
}

/**
 * Looks for a newick tree in an alignment file
 * @name trees.GetTreeString
 * @param {String|Bool} look_for_newick_tree - If a string, sanitizes and returns the string.
 * If TRUE, search the alignment file for a newick tree. If FALSE, the user will be prompted for a nwk tree file.
 * @returns {String} a newick tree string
 */
lfunction trees.GetTreeString(look_for_newick_tree) {


    UseModel(USE_NO_MODEL);

    if (Type(look_for_newick_tree) == "String") {
        treeString = trees.GetTreeString._sanitize(look_for_newick_tree);
    } else {
        if (look_for_newick_tree == FALSE) {
            utility.SetEnvVariable("IS_TREE_PRESENT_IN_DATA", FALSE);
        }

        if (utility.GetEnvVariable("IS_TREE_PRESENT_IN_DATA")) {
            treeString = trees.GetTreeString._sanitize(utility.GetEnvVariable("DATAFILE_TREE"));
            utility.SetEnvVariable("IS_TREE_PRESENT_IN_DATA", TRUE);
        }

        if (!utility.GetEnvVariable("IS_TREE_PRESENT_IN_DATA")) {
            SetDialogPrompt("Please select a tree file for the data:");
            utility.SetEnvVariable ("LAST_FILE_IO_EXCEPTION", None);
            utility.ToggleEnvVariable ("SOFT_FILE_IO_EXCEPTIONS",TRUE);
            fscanf(PROMPT_FOR_FILE, REWIND, "Raw", treeString);

            look_for_newick_tree = utility.getGlobalValue ("LAST_FILE_PATH");

            if (None != utility.GetEnvVariable ("LAST_FILE_IO_EXCEPTION")) {
                if (utility.getGlobalValue ("LAST_RAW_FILE_PROMPT") ==  utility.getGlobalValue ("terms.trees.neighbor_joining")) {
                    datafilter_name = utility.GetEnvVariable(utility.getGlobalValue ("terms.trees.data_for_neighbor_joining"));
                    assert (Type (datafilter_name) == "String", "Expected a string for the datafilter to  build a NJ tree from");
                    treeString = tree.infer.NJ (datafilter_name, None);
                } else {
                    assert (0, utility.GetEnvVariable ("LAST_FILE_IO_EXCEPTION"));
                }
                utility.SetEnvVariable ("LAST_FILE_IO_EXCEPTION", None);
            }

            utility.ToggleEnvVariable ("SOFT_FILE_IO_EXCEPTIONS",None);
            fprintf(stdout, "\n");


            if (regexp.Find(treeString, "^#NEXUS")) {
                ExecuteCommands(treeString);

                if (!utility.GetEnvVariable("IS_TREE_PRESENT_IN_DATA")) {
                    fprintf(stdout, "\n> **This NEXUS file doesn't contain a valid tree block**");
                    return 1;
                }

                nftm = utility.GetEnvVariable("NEXUS_FILE_TREE_MATRIX");


                if (Rows(nftm) > 1) {
                    ChoiceList(treeChoice, "Select a tree", 1, SKIP_NONE, nftm);
                    if (treeChoice < 0) {
                        return 1;
                    }
                    treeString = nftm[treeChoice][1];
                } else {
                    treeString = nftm[0][1];
                }
            } else {
                start = (treeString $ "\\(")[0];
                if (start < 0) {
                    fprintf(stdout, "\n> **This doesn't seem to be a valid Newick string file**. Can't find the opening parenthesis. ``", treeString, "``\n");
                    return 1;
                } else {
                    parenCounter = 1;
                    current = start + 1;
                    while (current < Abs(treeString) && parenCounter) {
                        char = treeString[current];
                        if (char == "(") {
                            parenCounter += 1;
                        } else {
                            if (char == ")") {
                                parenCounter += (-1);
                            }
                        }
                        current += 1;
                    }

                    if (parenCounter) {
                        fprintf(stdout, "\n> ** This doesn't seem to be a valid Newick string file**. Can't match the parentheses. \n``", treeString, "``\n");
                        return 1;
                    }

                    treeString = treeString[start][current - 1];
                }

                treeString = trees.GetTreeString._sanitize(treeString);

            }
        }
    }

    return
    {
        utility.getGlobalValue("terms.data.file"): look_for_newick_tree,
        utility.getGlobalValue("terms.data.tree"): treeString
    };
}


///////// TODO: HOW TO UN-HARDCODE THIS FUNCTION? Regular strategies are not working. //////////////
/**
 * Partitions a tree by assigning nodes to either being internal or leaf
 * @name trees.PartitionTree
 * @param {Dictionary} avl - an AVL representation of the tree to be partitioned
 * @returns nothing
 */
lfunction trees.PartitionTree(avl, l) {
    for (k = 0; k < Abs(avl); k += 1) {
        if ((avl[k])["Parent"]) {
            if (Abs((avl[k])["Children"])) {
                l[(avl[k])["Name"]] = "internal";
            } else {
                l[(avl[k])["Name"]] = "leaf";
            }
        }
    }
}

/**
 * Loads an annatotaed tree topology from a newick tree
 * @name trees.LoadAnnotatedTopology
 * @param {String|Bool} look_for_newick_tree - If a string, sanitizes and
 *  returns the string. If TRUE, search the alignment file for a newick tree. If
 *  FALSE, the user will be prompted for a nwk tree file.
 * @returns {Dictionary} an annotated tree
 */
lfunction trees.LoadAnnotatedTopology(look_for_newick_tree) {
    return trees.ExtractTreeInfo(trees.GetTreeString(look_for_newick_tree));
}

/**
 * @name trees.LoadAnnotatedTopologyAndMap
 * @param dataset_name
 * @returns {Dictionary} an annotated tree
 */
lfunction trees.LoadAnnotatedTopologyAndMap(look_for_newick_tree, mapping) {

    reverse = {};

    for (k,v; in; mapping) {
      reverse[v] = k;
    }

    //utility.ForEach(utility.Keys(mapping), "_key_", "`&reverse`[`&mapping`[_key_]] = _key_");

    io.CheckAssertion("Abs (`&mapping`) == Abs (`&reverse`)", "The mapping between original and normalized tree sequence names must be one to one");
    utility.ToggleEnvVariable("TREE_NODE_NAME_MAPPING", reverse);
    result = trees.ExtractTreeInfo(trees.GetTreeString(look_for_newick_tree));
    utility.ToggleEnvVariable("TREE_NODE_NAME_MAPPING", None);

    return result;
}

/**
 * Loads a tree topology with node name label mappings and annotations from a list of partitions
 * @name trees.LoadAnnotatedTreeTopology.match_partitions
 * @param {Matrix} partitions - a 1xN vector of partition names, typically
 * retrieved from the `partitions` key in the dictionary returned from an
 * alignments parser function.
 * @param {Dictionary} mapping - a mapping of node names to labels
 * @returns {Dictionary} of matched partitions
 * @example
 *  hky85_nucdata_info = alignments.ReadNucleotideAlignment(file_name, "hky85.nuc_data", "hky85.nuc_filter");
 *  name_mapping = {
 *   "HUMAN":"HUMAN",
 *   "CHIMP":"CHIMP",
 *   "BABOON":"BABOON",
 *   "RHMONKEY":"RHMONKEY",
 *   "COW":"COW",
 *   "PIG":"PIG",
 *   "HORSE":"HORSE",
 *   "CAT":"CAT",
 *   "MOUSE":"MOUSE",
 *   "RAT":"RAT"
 *  };
 *  trees.LoadAnnotatedTreeTopology.match_partitions(hky85_nucdata_info[terms.json.partitions], name_mapping);
 *  =>
 *  {
 *   "0":{
 *     "name":"default",
 *     "filter-string":"",
 *     "tree": *    }
 *  }
 */
lfunction trees.LoadAnnotatedTreeTopology.match_partitions(partitions, mapping) {

    partition_count = Rows(partitions);
    partrees = {};

    tree_matrix = utility.GetEnvVariable("NEXUS_FILE_TREE_MATRIX");

    if (Type(tree_matrix) == "Matrix") {
        io.CheckAssertion("Rows(`&tree_matrix`) >= partition_count", "The number of trees in the NEXUS block cannot be smaller than the number of partitions in the file");
        for (i = 0; i < partition_count; i += 1) {
            partrees + {
                utility.getGlobalValue("terms.data.name"): partitions[i][0],
                utility.getGlobalValue("terms.data.filter_string"): partitions[i][1],
                utility.getGlobalValue("terms.data.tree"): trees.LoadAnnotatedTopologyAndMap(tree_matrix[i][1], mapping)
            };
        }
    } else { // no tree matrix; apply the same tree to all partitions

        tree_info = trees.LoadAnnotatedTopologyAndMap(TRUE, mapping);

        for (i = 0; i < partition_count; i += 1) {
            partrees + {
                utility.getGlobalValue("terms.data.name"): partitions[i][0],
                utility.getGlobalValue("terms.data.filter_string"): partitions[i][1],
                utility.getGlobalValue("terms.data.tree"): tree_info
            };
        }
    }

    return partrees;

}

/**
 * @name trees.branch_names
 * @param tree
 * @param respect_case
 * @returns result
 */
lfunction trees.branch_names(tree, respect_case) {
    result = {};
    branch_names = BranchName(tree, -1);
    branch_count = Columns(branch_names) - 1;
    if (respect_case) {
        for (k = 0; k < branch_count; k += 1) {
            result + branch_names[k];
        }
    } else {
        for (k = 0; k < branch_count; k += 1) {
            result + (branch_names[k] && 1);
        }

    }
    return result;
}


/**
 * @name trees.KillZeroBranches
 * Given a tree dict and a list of matching parameter estimates
 * this returns a modified tree dict with all zero-branch length
 * internal branches removed and modifies in in place
 * @param tree - dict of the tree object
 * @param estimates - dict with branch length estimates
 * @param branch_set - if not null, treat as test set and delete branches form it as well
 * @param zero_internal -- store branches that were deleted here
 * @return modified tree
 */

lfunction trees.KillZeroBranches (tree, estimates, branch_set, zero_internal) {

    for (branch, value; in; tree [^"terms.trees.partitioned"]) {
        if (value == ^"terms.tree_attributes.internal") {
            if (estimates / branch) {
                if ((estimates[branch])[^"terms.fit.MLE"] < 1e-10) {
                    zero_internal + branch;
                }
            }
        }
    }

    if (Abs (zero_internal) > 0) { // has zero internal branches
        Topology T = tree[^"terms.trees.newick_annotated"];
        T -  Columns(zero_internal);
        if (None != branch_set) {
            for (branch; in; zero_internal) {
                branch_set - branch;
            }
        }
        return trees.ExtractTreeInfoFromTopology (&T);

    }

    return tree;

}

/**
 * @name trees.RootTree
 * @param {Dict} tree_info
 * @param {String} root on this node (or prompt if empty)
 * @returns a {Dictionary} (same as ExtractTreeInfo) for the rerooted tree
 */

lfunction trees.RootTree(tree_info, root_on) {
    if (Type (root_on) != "String") {
        root_on = io.SelectAnOption (tree_info[^"terms.trees.partitioned"],
                                     "Select a root");
    }


    io.CheckAssertion("`&tree_info`[^'terms.trees.partitioned']/`&root_on`", "Not a valid root choice '" + root_on + "'");

    Topology T = tree_info[^"terms.trees.newick_with_lengths"];

    utility.ToggleEnvVariable("ACCEPT_ROOTED_TREES", TRUE);
    tree_info = trees.ExtractTreeInfo(RerootTree (T, root_on));
    utility.ToggleEnvVariable("ACCEPT_ROOTED_TREES", None);

    return {
        ^"terms.trees.root"   : root_on,
        ^"terms.data.tree"    : tree_info
    };
}

/**
 * @name trees.ExtractTreeInfoFromTopology
 * @param {String} Topology
 * @returns a {Dictionary} of the following tree information :
 * - newick string
 * - newick string with branch lengths
 * - annotated string
 * - model map
 * - internal leaves
 * - list of models
 */
lfunction trees.ExtractTreeInfoFromTopology(topology_object) {


    branch_lengths = BranchLength(^topology_object, -1);
    branch_names   = BranchName(^topology_object, -1);
    branch_count   = utility.Array1D (branch_names) - 1;

    bls = {};

    for (k = 0; k < branch_count; k+=1) {
        if (branch_lengths[k] >= 0.) {
            bls [branch_names[k]] = branch_lengths[k];
        }
    }

    GetInformation(modelMap, ^topology_object);


    leaves_internals = {};
    flat_tree = (^topology_object) ^ 0;
    trees.PartitionTree(flat_tree, leaves_internals);

    utility.ToggleEnvVariable("INCLUDE_MODEL_SPECS", TRUE);
    T.str = "" + ^topology_object;
    utility.ToggleEnvVariable("INCLUDE_MODEL_SPECS", None);

    rooted = utility.Array1D ((flat_tree[(flat_tree[0])["Root"]])["Children"]) == 2;

    flat_tree       = None;
    branch_lengths  = None;
    branch_names    = None;
    branch_count    = None;

    return {
        ^"terms.trees.newick": Format(^topology_object, 1, 0),
        ^"terms.trees.newick_with_lengths": Format(^topology_object, 1, 1),
        ^"terms.branch_length": bls,
        ^"terms.trees.newick_annotated": T.str,
        ^"terms.trees.model_map": modelMap,
        ^"terms.trees.partitioned": leaves_internals,
        ^"terms.trees.model_list": Columns(modelMap),
        ^"terms.trees.rooted" : rooted,
        ^"terms.trees.meta" : Eval(^(topology_object+".__meta")),
        ^"terms.data.file" : file_name

    };
}

/**
 * @name trees.ExtractTreeInfo
 * @param {String} tree_string
 * @returns a {Dictionary} of the following tree information :
 * - newick string
 * - newick string with branch lengths
 * - annotated string
 * - model map
 * - internal leaves
 * - list of models
 */
lfunction trees.ExtractTreeInfo(tree_string) {

    if (Type (tree_string) == "AssociativeList") {
        file_name   = tree_string[^"terms.data.file"];
        tree_string = tree_string[^"terms.data.tree"];
    } else {
        file_name = None;
    }

    Topology T = tree_string;
    return trees.ExtractTreeInfoFromTopology (&T);

}

/**
 * @name trees.HasBranchLengths
 * @param {Dictionary} tree information object (e.g. as returned by LoadAnnotatedTopology)
 * @returns a {Boolean} to indicate whether the tree has a valid branch length array
 */

lfunction trees.HasBranchLengths (tree_info) {
    return utility.Array1D (tree_info [^"terms.trees.partitioned"]) == utility.Array1D (tree_info [^"terms.branch_length"]);
}

/**
 * @name trees.BootstrapSupport
 * @param {Dictionary} tree information object (e.g. as returned by LoadAnnotatedTopology)
 * @returns a {Dictionary} (could be empty or partially filled) with "node name" -> bootstrap support
 */

lfunction trees.BootstrapSupport (tree_info) {
    return utility.Map (utility.Filter (tree_info[^"terms.trees.meta"], "_value_", "(_value_/'bootstrap')"), "_value_",
                        "_value_['bootstrap']");
}

/**
 * Gets total branch count of supplied tree
 * @name trees.GetBranchCount
 * @param {String} tree_string
 * @returns {Number} total branch count
 */
function trees.GetBranchCount(tree_string) {
    Topology T = tree_string;
    return BranchCount(T) + TipCount(T);
}

/**
 * Sorts branch lengths in a descending order
 * @name trees.SortedBranchLengths
 * @param {String} tree_string
 * @returns {Number} total branch count
 */

lfunction trees.SortedBranchLengths(tree_string) {

    tree_count = trees.GetBranchCount(tree_string);
    Tree T = tree_string;

    branch_lengths = BranchLength(T,-1);

    sorted_bls = {};
    sorted_bls = {tree_count, 2}["branch_lengths[_MATRIX_ELEMENT_ROW_]*(_MATRIX_ELEMENT_COLUMN_==0)+_MATRIX_ELEMENT_ROW_*(_MATRIX_ELEMENT_COLUMN_==1)"];
    sorted_bls = sorted_bls%0;
    return sorted_bls;

}

/**
 * Return branch names
 * @name trees.BranchNames
 * @param {String} tree_string
 * @returns {Matrix} 1xN sorted branch names
 */

lfunction trees.BranchNames(tree) {
    tree_string = tree[^"terms.trees.newick"];
    Topology T = tree_string;
    branch_names = BranchName(T, -1);
    return branch_names;
}

/**
 * Compute branch labeling using parsimony
 * @name trees.ParsimonyLabel
 * @param 	{String} tree ID
 * @param 	{Dict} 	 leaf -> label
 					 labels may be missing for some of the leaves to induce partial labeling of the tree
 * @returns {Dict} 	 {"score" : value, "labels" : Internal Branch -> label}
 */

lfunction trees.ParsimonyLabel(tree_id, given_labels) {
   tree_avl = (^tree_id) ^ 0;
   label_values = utility.UniqueValues (given_labels);
   label_count  = utility.Array1D (label_values);
   labels = {};
   scores = {}; // node name -> score of optimal labeling staring at this now given parent state
   optimal_labeling = {}; // node name -> current node labeling which achieves this score
   resulting_labels = {}; // internal nodes -> label

   // pass 1 to fill in the score matrix
   for (k = 0; k < Abs (tree_avl) ; k += 1) {
   	 	node_name = (tree_avl[k])["Name"];
   	 	node_children = (tree_avl[k])["Children"];
   	 	c_count = Abs (node_children);
   	 	if (c_count) { // internal node
   	 		// first check to see if all the children are labeled

   	 		c_names = {c_count , 1};

			for (c = 0; c < c_count; c+=1) {
				c_names[c] = (tree_avl[node_children[c]])["Name"];
   	 			if (scores / c_names[c] == FALSE)  {
   	 				break;
   	 			}

   	 		}
   	 		if (c == c_count) {

   	 			scores [node_name] = {1, label_count};
   	 			labels [node_name] = {1, label_count};


    	 		for (parent_state = 0; parent_state < label_count; parent_state+=1) {
    	 			best_score = 1e100;
    	 			best_state = None;
    	 			for (my_state = 0; my_state < label_count; my_state+=1) {
    	 				my_score = 0;
						for (c = 0; c < c_count; c+=1) {
							my_score += (scores [c_names[c]])[my_state];
						}
						if (my_state != parent_state) {
							my_score += 1;
						}
						if (my_score < best_score) {
							best_score = my_score;
							best_state = my_state;
						}
					}
   	 	 			(scores [node_name])[parent_state] = best_score;
   	 	 			(labels [node_name])[parent_state] = best_state;
   	 	 		}

   	 	 	}
   	 	} else {
   	 		if (utility.Has (given_labels, node_name, "String")) {
   	 			node_label = given_labels[node_name];
   	 			scores [node_name] = {1, label_count};
   	 			labels [node_name] = {1, label_count};

   	 			for (z = 0; z < label_count; z+=1) {
   	 				if (node_label == label_values[z]) {
   	 					break;
   	 				}
   	 			}

    	 			for (i = 0; i < label_count; i+=1) {
   	 				(scores [node_name]) [i] = 1 - (z == i);
   	 				(labels [node_name]) [i] = z;
   	 			}
   	 		}
   	 	}
   }


   // pass 2 to choose the best state for subtree parents


   total_score = 0;

   for (k = 0; k < Abs (tree_avl); k += 1) {
   	 	node_name = (tree_avl[k])["Name"];
   	 	node_children = (tree_avl[k])["Children"];
   	 	c_count = Abs (node_children);
   	 	if (c_count) { // internal node
   	 		parent = (tree_avl[k])["Parent"];
   	 		if (parent > 0) {
   	 			if (utility.Has (scores, (tree_avl[parent])["Name"], "Matrix")) {
   	 				continue;
   	 			}
   	 		}

   	 		if (utility.Has (scores, node_name, "Matrix") == FALSE) {
   	 			continue;
   	 		}

   	 		best_score = 1e100;
   	 		best_label = None;

   	 		for (i = 0; i < label_count; i+=1) {
   	 			my_score = (scores [node_name])[i];

   	 			if (my_score < best_score) {
   	 				best_score = my_score;
   	 				best_label = i;
   	 			}
   	 		}

   	 		total_score += best_score;
   	 		resulting_labels [node_name] = best_label;
   	 	}
 	}

   tree_avl = (^tree_id) ^ 1;
   for (k = 2; k < Abs (tree_avl); k += 1) {
   	 node_name = (tree_avl[k])["Name"];
   	 if (utility.Array1D((tree_avl[k])["Children"])) {
   	 	parent = (tree_avl[k])["Parent"];
   	 	if (utility.Has (resulting_labels, (tree_avl[parent])["Name"], "Number")) {
   	 		parent = resulting_labels[(tree_avl[parent])["Name"]];
   	 		resulting_labels[node_name] = (labels [node_name])[parent];
   	 	}
   	 }
   }

   return {"score" : total_score, "labels" : utility.Map (resulting_labels, "_value_", "`&label_values`[_value_]")};
   // pass 1 to choose the best state for subtree parents
}

/**
 * Compute branch labeling using conjunction, i.e. node N is labeled 'X' iff
 * all of the nodes that are in the subtree rooted at 'N' are also labeled 'N'
 * @name trees.ConjunctionLabel
 * @param 	{String} tree ID
 * @param 	{Dict} 	 leaf -> label
 					 labels may be missing for some of the leaves to induce partial labeling of the tree
 * @returns {Dict} 	 {"labeled" : # of labeled internal  nodes, "labels" : Internal Branch -> label}
 */

lfunction trees.ConjunctionLabel (tree_id, given_labels) {
   tree_avl = (^tree_id) ^ 0;
   label_values = utility.UniqueValues (given_labels);
   label_count  = utility.Array1D (label_values);
   labels = {};
   resulting_labels = {}; // internal nodes -> label
   inodes_labeled = 0;

   // pass 1 to fill in the score matrix
   for (k = 0; k < Abs (tree_avl) ; k += 1) {
   	 	node_name = (tree_avl[k])["Name"];
   	 	node_children = (tree_avl[k])["Children"];
   	 	c_count = utility.Array1D (node_children);

   	 	if (c_count) { // internal node
   	 		// first check to see if all the children are labeled

   	 	   	child_labels = {};

			for (c = 0; c < c_count; c+=1) {
				c_name = (tree_avl[node_children[c]])["Name"];
   	 			if (utility.Has (labels, c_name, "String") == FALSE)  {
   	 				break;
   	 			}
   	 			child_labels [labels[c_name]] = TRUE;

   	 		}
   	 		if (c == c_count) { // all children labeled
    	 	   inodes_labeled += 1;
    	 	   if (utility.Array1D (child_labels) == 1) {
    	 	       labels [node_name] = (utility.Keys (child_labels))[0];
    	 	       resulting_labels[node_name] = labels [node_name];
    	 	   }
   	 	 	}
   	 	} else { // leaf
   	 		if (utility.Has (given_labels, node_name, "String")) {
                labels[node_name] = given_labels[node_name];
   	 		}
   	 	}
   }



   return {"labeled" : inodes_labeled, "labels" : resulting_labels};
   // pass 1 to choose the best state for subtree parents
}

/**
 * Compute branch labeling using conjunction, i.e. node N is labeled 'X' iff
 * SOME of the nodes that are in the subtree rooted at 'N' are also labeled 'N'
 * @name trees.ConjunctionLabel
 * @param 	{String} tree ID
 * @param 	{Dict} 	 leaf -> label
 					 labels may be missing for some of the leaves to induce partial labeling of the tree
 * @returns {Dict} 	 {"labeled" : # of labeled internal  nodes, "labels" : Internal Branch -> label}
 */

lfunction trees.DisjunctionLabel (tree_id, given_labels) {
   tree_avl = (^tree_id) ^ 0;
   label_values = utility.UniqueValues (given_labels);
   label_count  = utility.Array1D (label_values);
   labels = {};
   resulting_labels = {}; // internal nodes -> label
   inodes_labeled = 0;

   // pass 1 to fill in the score matrix
   for (k = 0; k < Abs (tree_avl) ; k += 1) {
   	 	node_name = (tree_avl[k])["Name"];
   	 	node_children = (tree_avl[k])["Children"];
   	 	c_count = utility.Array1D (node_children);

   	 	if (c_count) { // internal node
   	 		// first check to see if all the children are labeled

   	 	   	child_labels = {};

			for (c = 0; c < c_count; c+=1) {
				c_name = (tree_avl[node_children[c]])["Name"];
   	 			if (utility.Has (labels, c_name, "String") == FALSE)  {
   	 				break;
   	 			}
   	 			child_labels [labels[c_name]] = TRUE;

   	 		}
   	 		if (c > 0) { // all children labeled
    	 	   inodes_labeled += 1;
    	 	   if (utility.Array1D (child_labels) == 1) {
    	 	       labels [node_name] = (utility.Keys (child_labels))[0];
    	 	       resulting_labels[node_name] = labels [node_name];
    	 	   }
   	 	 	}
   	 	} else { // leaf
   	 		if (utility.Has (given_labels, node_name, "String")) {
                labels[node_name] = given_labels[node_name];
   	 		}
   	 	}
   }



   return {"labeled" : inodes_labeled, "labels" : resulting_labels};
   // pass 1 to choose the best state for subtree parents
}

/**
 * Annotate a tree string with using user-specified labels
 * @name trees.ParsimonyLabel
 * @param 	{String} tree ID
 * @param 	{Dict/String} 	 node_name -> label OR (if string) (node_name) => name + annotation
 * @param {String} a pair of characters to enclose the label in
 * @param {Dict/String} node_name -> length OR (if string) (node_name) => length annotation
 * @return {String} annotated string
 */

lfunction tree.Annotate (tree_id, labels, chars, doLengths) {
	theAVL = (^tree_id)^0;
	_ost = "";
	_ost * 256;

	lastLevel = 0;
	treeSize  = Abs(theAVL);
	treeInfo  = theAVL[0];
	rootIndex = treeInfo["Root"];
	lastDepth = 0;

	for (nodeIndex = 1; nodeIndex < treeSize; nodeIndex += 1) {
        nodeInfo = theAVL[nodeIndex];
        myDepth = nodeInfo["Depth"];
        if (lastDepth < myDepth) {
            if (lastDepth) {
                _ost * ",";
            }
            for (pidx = lastDepth; pidx < myDepth; pidx += 1) {
                _ost * "(";
            }
        } else {
            if (lastDepth > myDepth) {
                for (pidx = myDepth; pidx < lastDepth; pidx += 1) {
                    _ost * ")";
                }
            } else {
                _ost * ",";
            }
        }


        if (Type (labels) == "String") {
            _ost * Call (labels, nodeInfo["Name"]);
        } else {
            _ost * nodeInfo["Name"];
            if (labels / nodeInfo["Name"]) {
                if (Abs(labels[nodeInfo["Name"]])) {
                    _ost * (chars[0] + labels[nodeInfo["Name"]] + chars[1]);
                }
            }
        }

        if (doLengths ) {

             if (nodeIndex < treeSize - 1) {
                _ost * ":";
                if (Type (doLengths) == "String") {
                    _ost * Call (doLengths, nodeInfo);
                } else {
                    _ost * (""+nodeInfo ["Length"]);

                }
            }
        }
        lastDepth = myDepth;
	}

	_ost * 0;
	return _ost;
}


/**
 * Generate a symmetric binary tree on N leaves (only perfectly symmetric if N is a power of 2)
 * @name trees.GenerateSymmetricTree
 * @param 	{Number} N : number of leavers
 * @param   {Bool} rooted : whether the tree is rooted
 * @param 	{None/String} branch_name   : if a string, then it is assumed to be a function with an integer argument (node index) that generates branch names
                                        default is to use numeric names
 * @param 	{None/String} branch_length : if a string, then it is assumed to be a function with no arguments that generates branch lengths
 * @return  {String} Newick tree string
 */

lfunction tree.GenerateSymmetricTree (N, rooted, branch_name, branch_length) {
    assert (N>=2, "Can't generate trees with fewer than 2 leaves");
    internal_nodes = N-2;
    if (rooted) {
        internal_nodes += 1;
    }

    total_nodes = N + internal_nodes;
    flat_tree  = {total_nodes, 4}["-1"];


    // each row represents the corresponding node (in its' index), [parent, child 1, child 2, child 3] record
    // each node is represented with a [0, total_nodes-1] integer
    // -1 means NULL
    // for example ((1,2),(3,4)) will be represented as
        /*
           {4,-1,-1,-1}
           {4,-1,-1,-1}
           {5,-1,-1,-1}
           {5,-1,-1,-1}
           {6,0,1,-1}
           {6,2,3,-1}
           {-1,4,5,-1}
        */

    current_parent_node = N;
    current_child_node  = 0;


    while (current_parent_node < total_nodes) {
        flat_tree[current_child_node][0] = current_parent_node;
        flat_tree[current_child_node+1][0] = current_parent_node;
        flat_tree[current_parent_node][1] = current_child_node;
        if (current_child_node < total_nodes-2) {
            flat_tree[current_parent_node][2] = current_child_node + 1;
        }
        current_parent_node += 1;
        current_child_node += 2;
        if (current_child_node == N-1) {
            // if the number of leaves in not divisible for 2, we skip the unpaired leaf and attach it to the root
            flat_tree[total_nodes-1][3-(rooted!=FALSE)] = N-1;
            flat_tree[current_child_node][0] = total_nodes-1;
            current_child_node += 1;
        }
    }


    if (current_child_node < total_nodes - 1) { // rooted tree
        /*
        {
            {4, -1, -1, -1}
            {4, -1, -1, -1}
            {5, -1, -1, -1}
            {5, -1, -1, -1}
            {-1, 0, 1, -1}
            {-1, 2, 3, -1}
            }
        */
        if (flat_tree[total_nodes-1][3] < 0) { // doesn't already have the 3rd leaf attached
            flat_tree[total_nodes-1][3] = total_nodes-2;
            flat_tree[total_nodes-2][0] = total_nodes-1;
            flat_tree[flat_tree[total_nodes-1][1]] = total_nodes-1;
            flat_tree[flat_tree[total_nodes-1][2]] = total_nodes-1;
        }
    }

    return tree._NewickFromMatrix (&flat_tree, total_nodes-1, branch_name, branch_length);

}


/**
 * Generate a ladder tree on N leaves
 * @name trees.GenerateSymmetricTree
 * @param 	{Number} N : number of leavers
 * @param   {Bool} rooted : whether the tree is rooted
 * @param 	{None/String} branch_name   : if a string, then it is assumed to be a function with an integer argument (node index) that generates branch names
                                        default is to use numeric names
 * @param 	{None/String} branch_length : if a string, then it is assumed to be a function with no arguments that generates branch lengths
 * @return  {String} Newick tree string
 */

lfunction tree.GenerateLadderTree (N, rooted, branch_name, branch_length) {
    assert (N>=2, "Can't generate trees with fewer than 2 leaves");
    internal_nodes = N-2;
    if (rooted) {
        internal_nodes += 1;
    }

    total_nodes = N + internal_nodes;
    flat_tree  = {total_nodes, 4}["-1"];


    current_parent_node = N;
    current_child_node  = 0;

    // create the first cherry

    flat_tree[current_child_node][0] = current_parent_node;
    flat_tree[current_child_node+1][0] = current_parent_node;
    flat_tree[current_parent_node][1] = current_child_node;
    flat_tree[current_parent_node][2] = current_child_node + 1;
    current_child_node = 2;
    current_parent_node += 1;


    while (current_parent_node < total_nodes) {
        flat_tree[current_child_node][0] = current_parent_node;
        flat_tree[current_parent_node-1][0] = current_parent_node;
        flat_tree[current_parent_node][1] = current_parent_node-1;
        flat_tree[current_parent_node][2] = current_child_node;
        current_parent_node += 1;
        current_child_node += 1;

    }


    if (current_child_node < N) { // unrooted tree
        flat_tree[total_nodes-1][3] = N-1;
        flat_tree[N-1][0] = total_nodes-1;
    }

    return tree._NewickFromMatrix (&flat_tree, total_nodes-1, branch_name, branch_length);

}

/**
 * Generate a RANDOM tree on N leaves
 * @name trees.GenerateRandomTree
 * @param 	{Number} N : number of leavers
 * @param   {Bool} rooted : whether the tree is rooted
 * @param 	{None/String} branch_name   : if a string, then it is assumed to be a function with an integer argument (node index) that generates branch names
                                        default is to use numeric names
 * @param 	{None/String} branch_length : if a string, then it is assumed to be a function with no arguments that generates branch lengths
 * @return  {String} Newick tree string
 */


lfunction tree._GenerateRandomTreeDraw2 (nodes) {
    r = Rows (nodes);
    n = Abs (nodes);
    n1 = Random (0,n) $ 1;
    do {
        n2 = Random (0,n) $ 1;
    } while (n1 == n2);


    nodes - r[n1];
    nodes - r[n2];

    return {"0" : +r[n1], "1" : +r[n2]};
}


lfunction tree.GenerateRandomTree (N, rooted, branch_name, branch_length) {
    assert (N>=2, "Can't generate trees with fewer than 2 leaves");
    internal_nodes = N-2;
    if (rooted) {
        internal_nodes += 1;
    }

    total_nodes = N + internal_nodes;
    flat_tree  = {total_nodes, 4}["-1"];


    available_to_join = {};
    for (k = 0; k < N; k+=1) {
        available_to_join[k] = TRUE;
    }


    current_parent_node = N;
    downto = 1 + (rooted == 0);


    while (Abs (available_to_join) > downto) {
        pair = tree._GenerateRandomTreeDraw2 (available_to_join);
        flat_tree[pair[0]][0] = current_parent_node;
        flat_tree[pair[1]][0] = current_parent_node;
        flat_tree[current_parent_node][1] = pair[0];
        flat_tree[current_parent_node][2] = pair[1];
        available_to_join[current_parent_node] = TRUE;
        current_parent_node += 1;
    }

    if (!rooted) { // attach the last node to the root
        available_to_join - (total_nodes-1);
        leaf_index = +((Rows(available_to_join))[0]);
        flat_tree[leaf_index][0] = total_nodes-1;
        flat_tree[total_nodes-1][3] = leaf_index;

    }

    return tree._NewickFromMatrix (&flat_tree, total_nodes-1, branch_name, branch_length);

}

lfunction tree._NewickFromMatrix (flat_tree, index, branch_name, branch_length) {
    if ((^flat_tree)[index][0] >= 0) { // not a root
        if ((^flat_tree)[index][1] < 0) { // a leaf
            if (branch_name) {
                name = Call (branch_name, index);
            } else {
                name = "" + index;
            }
            if (branch_length) {
                return name + ":" + Call (branch_length);
            }
            return name;
        } else {
            if (branch_length) {
                bl := ":" + Call (branch_length);
            } else {
                bl := "";
            }
            if (branch_name) {
                bn := Call (branch_name, index);
            } else {
                bn := "";
            }

            return "("   +  tree._NewickFromMatrix (flat_tree, (^flat_tree)[index][1], branch_name, branch_length) +
                     "," +  tree._NewickFromMatrix (flat_tree, (^flat_tree)[index][2], branch_name, branch_length) +
                     ")" + bn + bl;
        }
    }  else {
        if ((^flat_tree)[index][3] < 0) { // 2 root children
            return "(" + tree._NewickFromMatrix (flat_tree, (^flat_tree)[index][1], branch_name, branch_length) +
                    "," +  tree._NewickFromMatrix (flat_tree, (^flat_tree)[index][2], branch_name, branch_length) +
                    ")";
        } else {
            return "(" + tree._NewickFromMatrix (flat_tree, (^flat_tree)[index][1], branch_name, branch_length) +
            "," +  tree._NewickFromMatrix (flat_tree, (^flat_tree)[index][2], branch_name, branch_length) +
            "," +  tree._NewickFromMatrix (flat_tree, (^flat_tree)[index][3], branch_name, branch_length) +
            ")";
        }
    }
    return None;
}

/**
 * Infer a neighbor joining tree from a data filter using a specific distance
 * @name trees.infer.NJ
 * @param 	{String} datafilter the ID of an existing datafilter
 * @param 	{null/String/Matrix},
        null   : use the default distance calculation appropriate for the datatype
        String : a callback which takes (filter_id, seq1, seq2) arguments and returns d(seq1,seq2)
        Matrix : a precomputed distance matrix (same order of rows/column as datafilter names)

 * @return  {String} inferred tree string
 */

lfunction tree.infer.NJ (datafilter, distances) {

    flush_distances = FALSE;
    if (None == distances) { // use default distances
        type = alignments.FilterType (datafilter);
        if (type == ^"terms.nucleotide" || type == ^"terms.codon") {
            distances = distances.nucleotide.tn93 (datafilter, null, null);
            flush_distances = TRUE;
        } else {
            if (type == ^"terms.amino_acid" || type == ^"terms.binary") {
                distances = distances.p_distance (datafilter, null);
                flush_distances = TRUE;
            } else {
                io.ReportAnExecutionError ("No default distance available for filter type of `type`");
            }
        }
    }

    N = ^(datafilter + ".species");
    assert (N == Rows (distances), "Incompatible dimensions for the distance matrix and datafilter");

    if (N == 2) {
		d1 = distances[0][1]/2;
		treeNodes = {{0,1,d1__},
					 {1,1,d1__},
					 {2,0,0}};
		cladesInfo = {{2,0}};
	}
	else {
		if (N == 3) {
			d1 = (distances[0][1]+distances[0][2]-distances[1][2])/2;
			d2 = (distances[0][1]-distances[0][2]+distances[1][2])/2;
			d3 = (distances[1][2]+distances[0][2]-distances[0][1])/2;
			treeNodes = {{0,1,d1__},
						 {1,1,d2__},
						 {2,1,d3__}
						 {3,0,0}};

			cladesInfo = {{3,0}};
		}
		else {
			njm = (distances > 0) >= N;

			treeNodes 		= {2*N,3};

           for (i = 0; i < 2*N; i += 1) {
                treeNodes[i][0] = njm[i][0]; // node index; leaves in [0,N), internal nodes N and higher
                treeNodes[i][1] = njm[i][1]; // node depth relative to the root
                treeNodes[i][2] = njm[i][2]; // branch length
            }

            njm = None;
		}
	}

    if (flush_distances) {
        distances = None;
    }

    return tree._matrix2string (treeNodes, N, alignments.GetSequenceNames (datafilter), TRUE);

}

/* ____________________________________________*/

/**
 * Convert a matrix representation of a tree into Newick
 * @name trees._matrix2string
 * @param 	{Matrix} matrix_form : Kx3 matrix;
            each row represents ?
 * @param 	{Number} N : number of leaves
 * @param 	{Matrix/Dict} names : leaf index -> name
 * @param 	{Bool} do_lengths : include branch lengths in the output

 * @return  {String} newick tree string
 */

lfunction tree._matrix2string (matrix_form, N, names, do_lengths) {

	newick = ""; newick * 1024;

	p = 0;                 // tree depth of previous node
	k = 0;                 // current row in matrix_form

	m = matrix_form[0][1]; // current tree depth (0 == root)
	n = matrix_form[0][0]; // current node index

	while (m) {
		if (m>p) { // going down in the tree
			if (p) {
				newick * ",";
			}
			for (j=p; j<m; j += 1 ) {
				newick * "(";
			}
		} else { // going up in the tree
			if (m<p) {
				for (j=m;j<p; j += 1) {
					newick * ")";
				}
			}
			else {
				newick * ",";
			}
		}

		if (n < N) {
		    newick * names [n];
		}

		if (do_lengths) {
			newick * (":"+matrix_form[k][2]);
		}

		k += 1;
		p = m;
		n = matrix_form[k][0];
		m = matrix_form[k][1];
	}

	for (j=m;j<p; j += 1) {
		newick *")";
	}

	newick *0;
	return newick;
}
