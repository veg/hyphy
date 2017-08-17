LoadFunctionLibrary("libv3/all-terms.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/UtilityFunctions.bf");

function test_match_partitions () {

    // Test dataset
    file_name = "`PATH_TO_CURRENT_BF`data/CD2.nex";

    // Read in the nucleotide alignments into hky85.nuc_data and hk85.nuc_filter variables
    hky85_nucdata_info = alignments.ReadNucleotideAlignment(file_name, "hky85.nuc_data", "hky85.nuc_filter");
    nuc_data_default = { "0" : "hky85.nuc_filter"};

    hky85.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions(hky85_nucdata_info[terms.data.partitions], hky85.name_mapping);

    tree =  hky85.partitions_and_trees["0"];

    assert(utility.KeyExists(tree, terms.data.name), "names key not found");
    assert(utility.KeyExists(tree, terms.data.filter_string), "filter-string key not found");
    assert(utility.KeyExists(tree, terms.data.tree), "tree key not found");

}

function test_load_annotated_topology() {

    // Test dataset
    file_name = "`PATH_TO_CURRENT_BF`data/CD2.nex";

    // Read in the nucleotide alignments into hky85.nuc_data and hk85.nuc_filter variables
    hky85_nucdata_info = alignments.ReadNucleotideAlignment(file_name, "hky85.nuc_data", "hky85.nuc_filter");
    nuc_data_default = { "0" : "hky85.nuc_filter"};

    mapping = {};
    
    // Get tree from last file loaded
    tree_matrix = utility.GetEnvVariable("NEXUS_FILE_TREE_MATRIX");
    tree = trees.LoadAnnotatedTopology(tree_matrix[0][1]);

    assert(utility.KeyExists(tree, terms.trees.newick), "newick string key not found");
    assert(utility.KeyExists(tree, terms.trees.newick_with_lengths), "newick with lengths key not found");
    assert(utility.KeyExists(tree, terms.trees.newick_annotated), "newick annotated key not found");
    assert(utility.KeyExists(tree, terms.branch_length), "branch lengths key not found");
    assert(utility.KeyExists(tree, terms.trees.model_map), "model map key not found");
    assert(utility.KeyExists(tree, terms.trees.partitioned),"partitioned key not found");
    assert(utility.KeyExists(tree, terms.trees.model_list), "model_list key not found");

}

function test_get_branch_count() {

    // Test dataset
    file_name = "`PATH_TO_CURRENT_BF`data/CD2.nex";

    // Read in the nucleotide alignments into hky85.nuc_data and hk85.nuc_filter variables
    hky85_nucdata_info = alignments.ReadNucleotideAlignment(file_name, "hky85.nuc_data", "hky85.nuc_filter");
    nuc_data_default = { "0" : "hky85.nuc_filter"};

    mapping = {};
    
    // Get tree from last file loaded
    tree_matrix = utility.GetEnvVariable("NEXUS_FILE_TREE_MATRIX");
    tree = trees.LoadAnnotatedTopology(tree_matrix[0][1]);
    branch_count = trees.GetBranchCount(tree_matrix[0][1]);
    assert(branch_count == 16, "branch count is not as expected");

}

function test_sort_branch_lengths() {

    // Test dataset
    file_name = "`PATH_TO_CURRENT_BF`data/CD2.nex";

    // Read in the nucleotide alignments into hky85.nuc_data and hk85.nuc_filter variables
    hky85_nucdata_info = alignments.ReadNucleotideAlignment(file_name, "hky85.nuc_data", "hky85.nuc_filter");
    nuc_data_default = { "0" : "hky85.nuc_filter"};

    mapping = {};
    
    // Get tree from last file loaded
    tree_matrix = utility.GetEnvVariable("NEXUS_FILE_TREE_MATRIX");
    
    expected = {
                   {0.000799, 10}
                   {0.002015, 6}
                   {0.003108, 7}
                   {0.004349, 9}
                   {0.011873, 11}
                   {0.022733, 8}
                   {0.050958, 14}
                   {0.058611, 5}
                   {0.08509899999999999, 2}
                   {0.09795, 15}
                   {0.101856, 12}
                   {0.147969, 0}
                   {0.165787, 3}
                   {0.21343, 1}
                   {0.264806, 4}
                   {0.340802, 13}
               };
 

    actual = trees.SortedBranchLengths(tree_matrix[0][1]);
    utility.ForEach(expected, "_key_", "assert(expected[utility.ForEach.r] == actual[utility.ForEach.r], \"not equal\n\")");

}

function test_get_branch_names() {

    // Test dataset
    file_name = "`PATH_TO_CURRENT_BF`data/CD2.nex";

    // Read in the nucleotide alignments into hky85.nuc_data and hk85.nuc_filter variables
    hky85_nucdata_info = alignments.ReadNucleotideAlignment(file_name, "hky85.nuc_data", "hky85.nuc_filter");
    nuc_data_default = { "0" : "hky85.nuc_filter"};

    mapping = {};
    
    // Get tree from last file loaded
    tree_matrix = utility.GetEnvVariable("NEXUS_FILE_TREE_MATRIX");
    tree = trees.LoadAnnotatedTopology(tree_matrix[0][1]);
    
    expected = {
        {"PIG", "COW", "Node3", "HORSE", "CAT", "Node2", "RHMONKEY", "BABOON", "Node9", "HUMAN", "CHIMP", "Node12", "Node8", "Node1", "RAT", "MOUSE", "Node0"}
    };

    actual = trees.BranchNames(tree);
    utility.ForEach(expected, "_key_", "assert(expected[utility.ForEach.r] == actual[utility.ForEach.r], \"not equal\n\")");

}

test_match_partitions();
test_load_annotated_topology();
test_get_branch_count();
test_sort_branch_lengths();
test_get_branch_names();
