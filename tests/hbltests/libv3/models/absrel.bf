LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/all-terms.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/models/codon/BRANCH_SITE.bf");

function test_absrel() {

    // Test dataset
    margin_of_error = .05;
    file_name = "`PATH_TO_CURRENT_BF`../data/CD2.nex";

    // Read in the nucleotide alignments into absrel.nuc_data and absrel.nuc_filter variables
    absrel.codon_data_info = alignments.SetModelTypeAndReadCodonDataSetFromPath("absrel.codon_data", file_name, "Universal");
    absrel.codon_data_info = alignments.LoadCodonDataFile("absrel.codon_data",  "absrel.codon_filter", absrel.codon_data_info);
    codon_data_default = { "0" : "absrel.codon_filter"};

    name_space = & model_BSREL;

    // define the model
    absrel.model = model.generic.DefineModel("models.codon.BSREL.ModelDescription",
        name_space, {
            "0": terms.local,
            "1": absrel.codon_data_info[terms.code]
        },
        codon_data_default,
        None);

    absrel.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions(
        absrel.codon_data_info[terms.data.partitions], 
        absrel.name_mapping);

    absrel.trees = utility.Map(absrel.partitions_and_trees, "_partition_", '_partition_[terms.data.tree]');
    absrel.filter_specification = alignments.DefineFiltersForPartitions(absrel.partitions_and_trees, "absrel.codon_data" , "absrel.filter.", absrel.codon_data_info);
    absrel.filter_names = utility.Map (absrel.filter_specification, "_partition_", '_partition_[terms.data.name]');

    lf_components = {1,2};
    lf_components[0] = "absrel.filter.default";
    lf_components[1] = "tree_0";

    tree = absrel.trees["0"];

    model_assignment = {
        terms.default: absrel.model
    };

    model.ApplyModelToTree(lf_components[1], 
                           tree, 
                           model_assignment, 
                           None);

    model_id_to_object = {
        name_space: absrel.model
    };

    LikelihoodFunction likelihoodFunction = (lf_components);
    Optimize(mles, likelihoodFunction);

    results = estimators.ExtractMLEs( & likelihoodFunction, model_id_to_object);

    // Ensure results are within limits
    actual_globals_mle = ((results[terms.global])[terms.transition_transversion_ratio])[terms.fit.MLE];
    expected_globals_mle = 2.592147288219786;
    assert(Abs(expected_globals_mle - actual_globals_mle) <= margin_of_error, "wrong global mle");

    // Test baboon
    baboon = ((results[terms.branch_length])["0"])["BABOON"];
    actual_baboon_mle = baboon[terms.fit.MLE];
    expected_baboon_mle = 0.001668634501409925;
    assert(Abs(expected_baboon_mle - actual_baboon_mle) <= margin_of_error, "wrong baboon mle");

    actual_baboon_time_parameter_id = (baboon[terms.timeParameter()])[terms.id];
    expected_baboon_time_parameter_id = terms.default_time;
    assert(actual_baboon_time_parameter_id == expected_baboon_time_parameter_id, "wrong evolutionary time parameter id");

    actual_baboon_time_parameter_mle = (baboon[terms.timeParameter()])[terms.fit.MLE];
    expected_baboon_time_parameter_mle = 0.001470878718119435;
    assert(Abs(expected_baboon_time_parameter_mle - actual_baboon_time_parameter_mle) <= margin_of_error, "wrong evolutionary time parameter mle");

}

function test_absrel_stepup() {

    // Test dataset
    margin_of_error = .05;
    file_name = "`PATH_TO_CURRENT_BF`../data/CD2.nex";

    // Read in the nucleotide alignments into absrel.nuc_data and absrel.nuc_filter variables
    absrel.codon_data_info = alignments.SetModelTypeAndReadCodonDataSetFromPath("absrel.codon_data", file_name, "Universal");
    absrel.codon_data_info = alignments.LoadCodonDataFile("absrel.codon_data",  "absrel.codon_filter", absrel.codon_data_info);
    codon_data_default = { "0" : "absrel.codon_filter"};

    name_space = & model_BSREL;

    // define the model
    absrel.model = model.generic.DefineModel("models.codon.BSREL.ModelDescription",
        name_space, {
            "0": terms.local,
            "1": absrel.codon_data_info[terms.code]
        },
        codon_data_default,
        None);

    model_assignment = {
        terms.default: absrel.model
    };



    absrel.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions(
        absrel.codon_data_info[terms.json.partitions], 
        absrel.name_mapping);

    absrel.trees = utility.Map(absrel.partitions_and_trees, "_partition_", '_partition_[terms.data.tree]');
    absrel.filter_specification = alignments.DefineFiltersForPartitions(absrel.partitions_and_trees, "absrel.codon_data" , "absrel.filter.", absrel.codon_data_info);
    absrel.filter_names = utility.Map (absrel.filter_specification, "_partition_", '_partition_[terms.dat.name]');

    tree = absrel.trees["0"];

    // Options
    use_existing_model_spec = 0;
    tree_info_map = {};

    models.codon.BSREL.StepUp(tree, use_existing_model_spec, model_assignment);

    //// Ensure results are within limits
    //actual_globals_mle = ((results[terms.global])[terms.transition_transversion_ratio])[terms.fit.MLE];
    //expected_globals_mle = 2.592147288219786;
    //assert(Abs(expected_globals_mle - actual_globals_mle) <= margin_of_error, "wrong global mle");

    //// Test baboon
    //baboon = ((results[terms.branch_length])["0"])["BABOON"];
    //actual_baboon_mle = baboon[terms.fit.MLE];
    //expected_baboon_mle = 0.001668634501409925;
    //assert(Abs(expected_baboon_mle - actual_baboon_mle) <= margin_of_error, "wrong baboon mle");

    //actual_baboon_time_parameter_id = (baboon[terms.timeParameter()])[terms.id];
    //expected_baboon_time_parameter_id = terms.default_time;
    //assert(actual_baboon_time_parameter_id == expected_baboon_time_parameter_id, "wrong evolutionary time parameter id");

    //actual_baboon_time_parameter_mle = (baboon[terms.timeParameter()])[terms.fit.MLE];
    //expected_baboon_time_parameter_mle = 0.001470878718119435;
    //assert(Abs(expected_baboon_time_parameter_mle - actual_baboon_time_parameter_mle) <= margin_of_error, "wrong evolutionary time parameter mle");

}

//test_absrel();
test_absrel_stepup();

