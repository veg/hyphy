/**
 * Loads file information into a namespace specified by prefix
 * Sets the following variables:
 * * <prefix>.codon_data.mapping
 * * <prefix>.codon_data.sites
 * * <prefix>.codon_data.species
 * * <prefix>.codon_data.unique_sites
 * * <prefix>.codon_data_info
 * * <prefix>.codon_filter.sequence_map
 * * <prefix>.codon_filter.site_freqs
 * * <prefix>.codon_filter.site_map
 * * <prefix>.codon_filter.sites
 * * <prefix>.codon_filter.species
 * * <prefix>.filter_names
 * * <prefix>.filter_specification
 * * <prefix>.name_mapping
 * * <prefix>.partition_count
 * * <prefix>.partitions_and_trees
 * * <prefix>.prefix
 * * <prefix>.sample_size
 * * <prefix>.selected_branches
 * * <prefix>.trees
 * @param prefix {String} : The namespace to prefix all file information variables with
 * @return nothing, the function sets variables within a namespace
 */

utility.SetEnvVariable ("MARKDOWN_OUTPUT", TRUE);

function load_file (prefix) {

    settings = None;

    if (Type (prefix) == "AssociativeList") {
        settings = prefix[utility.getGlobalValue("terms.settings")];
        prefix = prefix [utility.getGlobalValue("terms.prefix")];
    }


    codon_data_info = alignments.PromptForGeneticCodeAndAlignment(prefix+".codon_data", prefix+".codon_filter");

        /** example output
        {
            "sequences": 13,
            "sites": 897,
            "name-mapping": {
                "AF082576": "AF082576",
    ...
                "AF234767": "AF234767",
                "AF_231119": "AF-231119"
            },
            "partitions": {
                {
                    "SPAN_1", "0-587"
                }
                ...
                {
                    "SPAN_5", "1693-2690"
                }
            },
            "code": {
                {
                    14, 13, 14, 13, 7, 7, 7, 7, 19, 5, 19, 5, 2, 2, 3, 2, 12, 11, 12, 11, 6, 6, 6, 6, 19, 19, 19, 19, 1, 1, 1, 1, 16, 15, 16, 15, 8, 8, 8, 8, 20, 20, 20, 20, 4, 4, 4, 4, 10, 9, 10, 9, 5, 5, 5, 5, 10, 17, 18, 17, 1, 0, 1, 0
                }
            },
            "file": "/Volumes/sergei-raid/hyphy2.2/tests/hbltests/libv3/data/partitioned.nex",
            "stop": "TAA,TAG,TGA",
            "dataset": "slac.codon_data",
            "datafilter": "slac.codon_filter"
        }

        */

    sample_size=codon_data_info[utility.getGlobalValue("terms.data.sites")]*codon_data_info[utility.getGlobalValue("terms.data.sequences")];


    codon_data_info[utility.getGlobalValue("terms.data.sample_size")] = sample_size;
    upper_prefix = prefix && 1; //uppercase the prefix for json name
    codon_data_info[utility.getGlobalValue("terms.json.json")] = codon_data_info[utility.getGlobalValue("terms.data.file")] + "."+upper_prefix+".json";

    name_mapping = codon_data_info[utility.getGlobalValue("terms.data.name_mapping")];

        /**
            will contain "mapped" -> "original" associations with sequence names; or null if no mapping was necessary
        */

    if (None == name_mapping) { /** create a 1-1 mapping if nothing was done */
        name_mapping = {};
        utility.ForEach (alignments.GetSequenceNames (prefix+".codon_data"), "_value_", "`&name_mapping`[_value_] = _value_");
    }

    // check for duplicates
    duplicate_sequences = codon_data_info[^"terms.data.sequences"] - alignments.HasDuplicateSequences (codon_data_info[^"terms.data.datafilter"],-1);
    if (duplicate_sequences > 0) {
      fprintf(stdout, "\n-------\n", io.FormatLongStringToWidth(
      ">[WARNING] This dataset contains " + duplicate_sequences + " duplicate " + io.SingularOrPlural (duplicate_sequences, 'sequence', 'sequences') + ".
      Identical sequences do not contribute any information to the analysis and only slow down computation.
      Please consider removing duplicate or 'nearly' duplicate sequences,
      e.g. using https://github.com/veg/hyphy-analyses/tree/master/remove-duplicates
      prior to running selection analyses", 72),
      "\n-------\n");
    }

    utility.SetEnvVariable(utility.getGlobalValue ("terms.trees.data_for_neighbor_joining"),
                           codon_data_info[utility.getGlobalValue("terms.data.datafilter")]);

    partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (codon_data_info[utility.getGlobalValue("terms.data.partitions")], name_mapping);

    utility.SetEnvVariable(utility.getGlobalValue ("terms.trees.data_for_neighbor_joining"), None);

        /**  this will return a dictionary of partition strings and trees; one set per partition, as in
        {
            "0": {
                "name:" ... ,
                "filter-string": "0-587",
                "tree": {
                    "string": ...
                    "string_with_lengths": ...
                    "branch_lengths": {
                        "AF_231119": 0.00519475,
                        ...
                    },
                    "annotated_string": ... ,
                    "model_map": {
                        "AF_231119": "",
                        ...
                    },
                    "partitioned": {
                        "AF_231119": "leaf",
                     },
                    "model_list": {
                        {
                    }
                }
            },
            ...
        */


    partition_count = Abs (partitions_and_trees);


    // TODO: DE-HARDCODE "filter-string"
    utility.ForEachPair (partitions_and_trees,
                            "_key_",
                            "_value_",
                            '(`&partitions_and_trees`[_key_])[utility.getGlobalValue("terms.data.filter_string")] = selection.io.adjust_partition_string (_value_[utility.getGlobalValue("terms.data.filter_string")], 3*`&codon_data_info`[utility.getGlobalValue("terms.data.sites")])');
        /**
            ensure that all partitions fall on codon boundaries if they are contiguous
        */

    io.ReportProgressMessage ("", ">Loaded a multiple sequence alignment with **" + codon_data_info[utility.getGlobalValue("terms.data.sequences")] + "** sequences, **" + codon_data_info[utility.getGlobalValue("terms.data.sites")] + "** codons, and **" + partition_count + "** partitions from \`" + codon_data_info[utility.getGlobalValue("terms.data.file")] + "\`");

    if (utility.Has (settings, utility.getGlobalValue("terms.settings.branch_selector"), "String")) {
        selected_branches =  Call (settings[utility.getGlobalValue("terms.settings.branch_selector")], partitions_and_trees);
    } else {
        selected_branches = selection.io.defineBranchSets(partitions_and_trees);
    }

    // Place in own attribute called `tested`
     selection.io.json_store_key_value_pair (json, None, utility.getGlobalValue("terms.json.tested"), selected_branches);

        /**  this will return a dictionary of selected branches; one set per partition, like in
        {
            "0": {
                "NODE3": "test",
                "NODE6": "background",
        ...
                "NODE15": "test"
            },

            ...
            "4": {
                "NODE4": "test",
         ...
                 "NODE2": "background"
            }
        }
        */


    /** Input attribute to JSON **/
    json[utility.getGlobalValue("terms.json.input")] = {};
    (json[utility.getGlobalValue("terms.json.input")])[utility.getGlobalValue("terms.json.file")] =  codon_data_info[utility.getGlobalValue("terms.data.file")];
    (json[utility.getGlobalValue("terms.json.input")])[utility.getGlobalValue("terms.json.sequences")] = codon_data_info[utility.getGlobalValue("terms.data.sequences")];
    (json[utility.getGlobalValue("terms.json.input")])[utility.getGlobalValue("terms.json.sites")] = codon_data_info[utility.getGlobalValue("terms.data.sites")];
    (json[utility.getGlobalValue("terms.json.input")])[utility.getGlobalValue("terms.json.partition_count")] = partition_count;

    // The trees should go into input as well and they should be w/ their branch lengths but ONLY if they have any.



     filter_specification = alignments.DefineFiltersForPartitions (partitions_and_trees, "`prefix`.codon_data" , "`prefix`.filter.", codon_data_info);
    /** defines codon filters for each partition, and returns the (codon) sites mapped to each filter
    {
        {
            "0": {
                "name": "slac.filter.SPAN_1",
                "coverage": {
                    {
                        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, ...
                    }
                }
            },
            "1": {
                "name": "...",
                "coverage": "..."
            },
            ...
        }
    }
    */

    store_tree_information ();
 }

function store_tree_information () {
    // Place in own attribute called `tested`

     selection.io.json_store_key_value_pair (json, None, utility.getGlobalValue("terms.json.tested"), selected_branches);

        /**  this will return a dictionary of selected branches; one set per partition, like in
        {
            "0": {
                "NODE3": "test",
                "NODE6": "background",
        ...
                "NODE15": "test"
            },

            ...
            "4": {
                "NODE4": "test",
         ...
                 "NODE2": "background"
            }
        }
        */


    t = (partitions_and_trees["0"])[utility.getGlobalValue("terms.data.tree")];
    abs_branch_lengths = Abs(t[utility.getGlobalValue("terms.branch_length")]);

    if (abs_branch_lengths == 0){
        selection.io.json_store_key_value_pair (json,
                                             utility.getGlobalValue("terms.json.input"), utility.getGlobalValue("terms.json.trees"),
                                             utility.Map (partitions_and_trees, "_pt_", '(_pt_[terms.data.tree])[terms.trees.newick]')
                                             );
    } else {
        selection.io.json_store_key_value_pair (json,
                                             utility.getGlobalValue("terms.json.input"), utility.getGlobalValue("terms.json.trees"),
                                             utility.Map (partitions_and_trees, "_pt_", '(_pt_[terms.data.tree])[terms.trees.newick_with_lengths]')
                                             );
    }


    selection.io.json_store_key_value_pair (json, None, utility.getGlobalValue("terms.json.partitions"),
                                                         filter_specification);
     trees = utility.Map (partitions_and_trees, "_partition_", '_partition_[terms.data.tree]');


     filter_names = utility.Map (filter_specification, "_partition_", '_partition_[terms.data.name]');

     /* Store original name mapping */
     for (partition_index = 0; partition_index < partition_count; partition_index += 1) {

        selection.io.json_store_branch_attribute(json, utility.getGlobalValue ("terms.original_name"), utility.getGlobalValue ("terms.json.node_label"), display_orders[utility.getGlobalValue ("terms.original_name")],
                                         partition_index,
                                         name_mapping);
    }


}






function doGTR (prefix) {

    io.ReportProgressMessageMD (prefix, "nuc-fit", "Obtaining branch lengths and nucleotide substitution biases under the nucleotide GTR model");

    gtr_results = utility.Extend (parameters.helper.tree_lengths_to_initial_values (trees, None),
                                  {
                                    utility.getGlobalValue ("terms.global") : {
                                        terms.nucleotideRate ("A","C") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25} ,
                                        terms.nucleotideRate ("A","T") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25} ,
                                        terms.nucleotideRate ("C","G") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25} ,
                                        terms.nucleotideRate ("G","T") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25}
                                    }
                                 });


    //utility.ToggleEnvVariable("VERBOSITY_LEVEL", 10);

    gtr_results = estimators.FitGTR(filter_names,
                                         trees,
                                         gtr_results);


    KeywordArgument ("kill-zero-lengths", "Automatically delete internal zero-length branches for computational efficiency (will not affect results otherwise)", "Yes");

    kill0 = io.SelectAnOption (
        {
            "Yes":"Automatically delete internal zero-length branches for computational efficiency (will not affect results otherwise)",
            "No":"Keep all branches"
        },
        "The set of properties to use in the model") == "Yes";


    if (kill0) {
        for (index, tree; in; trees) {
            deleted = {};
            if (^(prefix + ".selected_branches") / index) {
                trees[index] = trees.KillZeroBranches (tree, (gtr_results[^"terms.branch_length"])[index], (^(prefix + ".selected_branches"))[index], deleted);
            } else {
                trees[index] = trees.KillZeroBranches (tree, (gtr_results[^"terms.branch_length"])[index], null, deleted);
            }

            if (utility.Array1D (deleted)) {
                io.ReportProgressMessageMD(prefix,  'selector', 'Deleted ' + Abs(deleted) + ' zero-length internal branches: \`' + Join (', ',utility.Values(deleted)) + '\`');
                for (i,v; in; deleted) {
                    (gtr_results[^"terms.branch_length"])[index] - v;
                }
            }
        }
        for (i = 0; i < partition_count; i+=1) {
            (partitions_and_trees[i])[^"terms.data.tree"] = trees[i];
        }
        store_tree_information ();
    }


    io.ReportProgressMessageMD (prefix, "nuc-fit", "* " +
        selection.io.report_fit (gtr_results, 0, 3*(^"`prefix`.sample_size")));



    /* Store nucleotide fit */
    gtr_rates = utility.Map(
                    utility.Map (gtr_results[utility.getGlobalValue("terms.global")], "_value_", '   {terms.fit.MLE : _value_[terms.fit.MLE]}'),
                    "_value_",
                    "_value_[terms.fit.MLE]");

    efv = (gtr_results[utility.getGlobalValue("terms.efv_estimate")])["VALUEINDEXORDER"][0];
    selection.io.json_store_lf_withEFV (json,

                                utility.getGlobalValue ("terms.json.nucleotide_gtr"),
                                gtr_results[utility.getGlobalValue ("terms.fit.log_likelihood")],
                                gtr_results[utility.getGlobalValue ("terms.parameters")] ,
                                3*(^"`prefix`.sample_size"),
                                gtr_rates,
                                efv,
                                display_orders[utility.getGlobalValue ("terms.json.nucleotide_gtr")]);



    /* Store branch lengths */
    for (partition_index = 0; partition_index < Abs(filter_specification); partition_index += 1) {
        selection.io.json_store_branch_attribute(json, utility.getGlobalValue ("terms.json.nucleotide_gtr"), utility.getGlobalValue ("terms.branch_length"), display_orders[terms.json.nucleotide_gtr],
                                         partition_index,
                                         selection.io.extract_branch_info((gtr_results[utility.getGlobalValue ("terms.branch_length")])[partition_index], "selection.io.branch.length"));
        }

}









/**
 * @name doPartitionMG
 * Can only be used after including shared-load-file
 * @return
 */
function doPartitionedMG (prefix, keep_lf) {
    io.ReportProgressMessageMD ("`prefix`", "codon-fit", "Obtaining the global omega estimate based on relative GTR branch lengths and nucleotide substitution biases");


    /**
        Declare per-partition branch length scalers
        slac.scaler_prefix_0
        slac.scaler_prefix_1
        etc
    */
    scaler_variables = utility.PopulateDict (0, partition_count, "`prefix`.scaler_prefix + '_' + _k_", "_k_");

    utility.ForEach (scaler_variables, "_value_", "parameters.DeclareGlobal(_value_, None);parameters.SetValue(_value_, 3);");

    partitioned_mg_results = estimators.FitMGREV(filter_names, trees, codon_data_info [utility.getGlobalValue("terms.code")], {
        utility.getGlobalValue("terms.run_options.model_type"): utility.getGlobalValue("terms.local"),
        utility.getGlobalValue("terms.run_options.proportional_branch_length_scaler"): scaler_variables,
        utility.getGlobalValue("terms.run_options.partitioned_omega"): selected_branches,
        utility.getGlobalValue("terms.run_options.retain_lf_object"): keep_lf
    }, gtr_results);


    io.ReportProgressMessageMD("`prefix`", "codon-fit", "* " + selection.io.report_fit (partitioned_mg_results, 0, (^"`prefix`.codon_data_info")[utility.getGlobalValue ("terms.data.sample_size")]));
    global_dnds = selection.io.extract_global_MLE_re (partitioned_mg_results, "^" + utility.getGlobalValue("terms.parameters.omega_ratio"));
    utility.ForEach (global_dnds, "_value_", 'io.ReportProgressMessageMD ("`prefix`", "codon-fit", "* " + _value_[utility.getGlobalValue("terms.description")] + " = " + Format (_value_[utility.getGlobalValue("terms.fit.MLE")],8,4));');

    if (partition_count > 1) {
       //partition_scalers = selection.io.extract_global_MLE_re (partitioned_mg_results, "^" + utility.getGlobalValue("terms.parameters.omega_ratio"));

    }

    /** extract and report dN/dS estimates */
}
