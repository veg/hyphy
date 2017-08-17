LoadFunctionLibrary("../IOFunctions.bf");
LoadFunctionLibrary("../terms-json.bf");
LoadFunctionLibrary("../convenience/regexp.bf");
LoadFunctionLibrary("../UtilityFunctions.bf");

/** @module trees */

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
 * @param {String} or {Bool} look_for_newick_tree - If a string, sanitizes and returns the string. If TRUE, search the alignment file for a newick tree. If FALSE, the user will be prompted for a nwk tree file. 
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
            fprintf(stdout, "\n> A tree was found in the data file: ``", utility.GetEnvVariable("DATAFILE_TREE"), "``\n>Would you like to use it? ");
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

                nftm = utility.GetEnvVariable("NEXUS_FILE_TREE_MATRIX")

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
 * Returns tree information from extractTreeInfo after calling getTreeString
 * @name trees.LoadAnnotatedTopology
 * @param {String} or {Bool} look_for_newick_tree - If a string, sanitizes and returns the string. If TRUE, search the alignment file for a newick tree. If FALSE, the user will be prompted for a nwk tree file. 
 * @returns {Dictionary} extracted tree information
 */
lfunction trees.LoadAnnotatedTopology(look_for_newick_tree) {
    return trees.ExtractTreeInfo(trees.GetTreeString(look_for_newick_tree));
}

/**
 * @name trees.LoadAnnotatedTopologyAndMap
 * @param dataset_name
 * @returns nothing
 */
lfunction trees.LoadAnnotatedTopologyAndMap(look_for_newick_tree, mapping) {

    reverse = {};
    utility.ForEach(utility.Keys(mapping), "_key_", "`&reverse`[`&mapping`[_key_]] = _key_");

    io.CheckAssertion("Abs (mapping) == Abs (reverse)", "The mapping between original and normalized tree sequence names must be one to one");
    utility.ToggleEnvVariable("TREE_NODE_NAME_MAPPING", reverse);

    result = trees.ExtractTreeInfo(trees.GetTreeString(look_for_newick_tree));
    utility.ToggleEnvVariable("TREE_NODE_NAME_MAPPING", None);
    return result;
}

/**
 * @name trees.LoadAnnotatedTreeTopology.match_partitions
 * @param partitions
 * @param mapping
 * @returns {Dictionary} of matched partitions
 */
lfunction trees.LoadAnnotatedTreeTopology.match_partitions(partitions, mapping) {


    partition_count = Rows(partitions);
    partrees = {};

    tree_matrix = utility.GetEnvVariable("NEXUS_FILE_TREE_MATRIX");

    if (Type(tree_matrix) == "Matrix") {
        io.CheckAssertion("Rows(`&tree_matrix`) >= partition_count", "The number of trees in the NEXUS block cannot be smaller than the number of partitions in the file");
        for (i = 0; i < partition_count; i += 1) {
            partrees + {
                "name": partitions[i][0],
                "filter-string": partitions[i][1],
                "tree": trees.LoadAnnotatedTopologyAndMap(tree_matrix[i][1], mapping)
            };
        }
    } else { // no tree matrix; allow if there is a single partition

        partrees + {
            "name": partitions[0][0],
            "filter-string": partitions[0][1],
            "tree": trees.LoadAnnotatedTopologyAndMap(TRUE, mapping)
        };
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
 * * newick string
 * * newick string with branch lengths
 * * annotated string
 * * model map
 * * internal leaves
 * * list of models
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
        "string": Format(T, 1, 0),
        "string_with_lengths": Format(T, 1, 1),
        terms.json.attribute.branch_length: trees.LoadAnnotatedTopology.bls,
        "annotated_string": T.str,
        "model_map": modelMap,
        "partitioned": leaves_internals,
        "model_list": Columns(modelMap)
    };
}
