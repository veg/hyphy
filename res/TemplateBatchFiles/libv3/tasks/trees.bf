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
 * Loads an annotated tree topology from a newick tree
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
function trees.ExtractTreeInfo(tree_string) {

    Topology T = tree_string;

    trees.LoadAnnotatedTopology.branch_lengths = BranchLength(T, -1);
    trees.LoadAnnotatedTopology.branch_names = BranchName(T, -1);

    trees.LoadAnnotatedTopology.bls = {};

    for (trees.LoadAnnotatedTopology.k = 0; trees.LoadAnnotatedTopology.k < Columns(trees.LoadAnnotatedTopology.branch_names) - 1; trees.LoadAnnotatedTopology.k += 1) {
        if (trees.LoadAnnotatedTopology.branch_lengths[trees.LoadAnnotatedTopology.k] >= 0.) {
            trees.LoadAnnotatedTopology.bls[trees.LoadAnnotatedTopology.branch_names[trees.LoadAnnotatedTopology.k]] =
                trees.LoadAnnotatedTopology.branch_lengths[trees.LoadAnnotatedTopology.k];
        }
    }

    GetInformation(modelMap, T);

    leaves_internals = {};

    trees.PartitionTree(T ^ 0, leaves_internals);

    utility.ToggleEnvVariable("INCLUDE_MODEL_SPECS", 1);
    T.str = "" + T;
    utility.ToggleEnvVariable("INCLUDE_MODEL_SPECS", None);


    return {
        terms.trees.newick: Format(T, 1, 0),
        terms.trees.newick_with_lengths: Format(T, 1, 1),
        terms.branch_length: trees.LoadAnnotatedTopology.bls,
        terms.trees.newick_annotated: T.str,
        terms.trees.model_map: modelMap,
        terms.trees.partitioned: leaves_internals,
        terms.trees.model_list: Columns(modelMap)
    };
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
function trees.SortedBranchLengths(tree_string) {

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
function trees.BranchNames(tree) {
    tree_string = tree[terms.trees.newick];
    Topology T = tree_string;
    branch_names = BranchName(T, -1);
    return branch_names;
}

