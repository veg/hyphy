LoadFunctionLibrary("../IOFunctions.bf");
LoadFunctionLibrary("../all-terms.bf");
LoadFunctionLibrary("../convenience/regexp.bf");
LoadFunctionLibrary("../UtilityFunctions.bf");
LoadFunctionLibrary("TreeTools");

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
    string = string ^ {
        {
            "\\)[0-9]+(\\.[0-9]*)?\:",
            "):"
        }
    };

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
            fprintf(stdout, "\n\n>A tree was found in the data file: ``", utility.GetEnvVariable("DATAFILE_TREE"), "``\n\n>Would you like to use it (y/n)? ");
            fscanf(stdin, "String", response);
            if (response == "n" || response == "N") {
                utility.SetEnvVariable("IS_TREE_PRESENT_IN_DATA", FALSE);
                IS_TREE_PRESENT_IN_DATA = 0;
            } else {
                treeString = trees.GetTreeString._sanitize(utility.GetEnvVariable("DATAFILE_TREE"));
                utility.SetEnvVariable("IS_TREE_PRESENT_IN_DATA", TRUE);
            }
            fprintf(stdout, "\n\n");
        }

        if (!utility.GetEnvVariable("IS_TREE_PRESENT_IN_DATA")) {
            SetDialogPrompt("Please select a tree file for the data:");
            fscanf(PROMPT_FOR_FILE, REWIND, "Raw", treeString);
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

    return treeString;
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
    utility.ForEach(utility.Keys(mapping), "_key_", "`&reverse`[`&mapping`[_key_]] = _key_");

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
    } else { // no tree matrix; allow if there is a single partition

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

    Topology T = tree_string;

    branch_lengths = BranchLength(T, -1);
    branch_names   = BranchName(T, -1);
    branch_count   = utility.Array1D (branch_names) - 1;

    bls = {};

    for (k = 0; k < branch_count; k+=1) {
        if (branch_lengths[k] >= 0.) {
            bls [branch_names[k]] = branch_lengths[k];
        }
    }

    GetInformation(modelMap, T);
    

    leaves_internals = {};
    flat_tree = T ^ 0;
    trees.PartitionTree(flat_tree, leaves_internals);

    utility.ToggleEnvVariable("INCLUDE_MODEL_SPECS", TRUE);
    T.str = "" + T;
    utility.ToggleEnvVariable("INCLUDE_MODEL_SPECS", None);

    rooted = utility.Array1D ((flat_tree[(flat_tree[0])["Root"]])["Children"]) == 2;

    DeleteObject (flat_tree, branch_lengths, branch_names, branch_count);

    return {
        ^"terms.trees.newick": Format(T, 1, 0),
        ^"terms.trees.newick_with_lengths": Format(T, 1, 1),
        ^"terms.branch_length": bls,
        ^"terms.trees.newick_annotated": T.str,
        ^"terms.trees.model_map": modelMap,
        ^"terms.trees.partitioned": leaves_internals,
        ^"terms.trees.model_list": Columns(modelMap),
        ^"terms.trees.rooted" : rooted
    };
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
   label_values = utility.Values (given_labels);
   label_count  = utility.Array1D (label_values);
   labels = {};
   scores = {}; // node name -> score of optimal labeling staring at this now given parent state 
   optimal_labeling = {}; // node name -> current node labeling which achieves this score
   resulting_labels = {}; // internal nodes -> label
   
   // pass 1 to fill in the score matrix
   for (k = 0; k < Abs (tree_avl) ; k += 1) {
   	 	node_name = (tree_avl[k])["Name"];
   	 	node_children = (tree_avl[k])["Children"];
   	 	c_count = utility.Array1D (node_children);
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
   	 	c_count = utility.Array1D (node_children);
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
   label_values = utility.Values (given_labels);
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
 * Annotate a tree string with using user-specified labels  
 * @name trees.ParsimonyLabel
 * @param 	{String} tree ID
 * @param 	{Dict/String} 	 node_name -> label OR (if string) (node_name) => name + annotation
 * @param {String} a pair of characters to enclose the label in 
 * @param {Bool} whether or not to include branch lentghs
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

        if (doLengths) {
            if (nodeIndex < treeSize - 1) {
                _ost * ":";
                _ost * (""+nodeInfo ["Length"]); 
            }
        }
        lastDepth = myDepth;
	}
	
	_ost * 0;
	return _ost;
}


