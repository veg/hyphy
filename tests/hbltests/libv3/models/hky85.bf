LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/all-terms.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/models/DNA/HKY85.bf");

function test_hky85() {

    // Test dataset
    margin_of_error = .05;
    file_name = "`PATH_TO_CURRENT_BF`../data/CD2.nex";

    // Read in the nucleotide alignments into hky85.nuc_data and hk85.nuc_filter variables
    hky85_nucdata_info = alignments.ReadNucleotideAlignment(file_name, "hky85.nuc_data", "hky85.nuc_filter");
    nuc_data_default = { "0" : "hky85.nuc_filter"};

    name_space = & model_HKY85;

    // define the model
    hky85.model = model.generic.DefineModel("models.DNA.HKY85.ModelDescription",
        name_space, {
            "0": terms.local
        },
        nuc_data_default,
        None);

    //fprintf(stdout, hky85.model);
    hky85.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions(hky85_nucdata_info[terms.data.partitions], hky85.name_mapping);
    hky85.trees = utility.Map(hky85.partitions_and_trees, "_partition_", '_partition_[terms.data.tree]');

    lf_components = {1,2};
    lf_components[0] = "hky85.nuc_filter";
    lf_components[1] = "tree_0";

    tree = hky85.trees["0"];

    model_assignment = {
        terms.default: hky85.model
    };

    model.ApplyModelToTree(lf_components[1], 
                           tree, 
                           model_assignment, 
                           None);

    model_id_to_object = {
        name_space: hky85.model
    };

    LikelihoodFunction likelihoodFunction = (lf_components);
    Optimize(mles, likelihoodFunction);

    results = estimators.ExtractMLEs( & likelihoodFunction, model_id_to_object);

    // Ensure results are within limits
    actual_globals_mle = ((results[terms.global])[terms.parameters.transition_transversion_ratio])[terms.fit.MLE];
    expected_globals_mle = 2.592147288219786;
    assert(Abs(expected_globals_mle - actual_globals_mle) <= margin_of_error, "wrong global mle");

    // Test baboon
    baboon = ((results[terms.branch_length])["0"])["BABOON"];
    actual_baboon_mle = baboon[terms.fit.MLE];
    expected_baboon_mle = 0.001668634501409925;
    assert(Abs(expected_baboon_mle - actual_baboon_mle) <= margin_of_error, "wrong baboon mle");

    actual_baboon_time_parameter_id = (baboon[terms.timeParameter()])[terms.id];
    expected_baboon_time_parameter_id = terms.parameters.default_time;
    assert(actual_baboon_time_parameter_id == expected_baboon_time_parameter_id, "wrong evolutionary time parameter id");

    actual_baboon_time_parameter_mle = (baboon[terms.timeParameter()])[terms.fit.MLE];
    expected_baboon_time_parameter_mle = 0.001470878718119435;
    assert(Abs(expected_baboon_time_parameter_mle - actual_baboon_time_parameter_mle) <= margin_of_error, "wrong evolutionary time parameter mle");

}

test_hky85();

