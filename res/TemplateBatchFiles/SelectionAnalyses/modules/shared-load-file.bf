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


function load_nuc_file (prefix) {

    nuc_data_info = alignments.ReadNucleotideDataSet(prefix+".nuc_data", null);
    sample_size=nuc_data_info[utility.getGlobalValue("terms.data.sites")]*nuc_data_info[utility.getGlobalValue("terms.data.sequences")];


    nuc_data_info[utility.getGlobalValue("terms.data.sample_size")] = sample_size;
    upper_prefix = prefix && 1; //uppercase the prefix for json name
    nuc_data_info[utility.getGlobalValue("terms.json.json")] = nuc_data_info[utility.getGlobalValue("terms.data.file")] + "."+upper_prefix+".json";

    name_mapping = nuc_data_info[utility.getGlobalValue("terms.data.name_mapping")];

        /**
            will contain "mapped" -> "original" associations with sequence names; or null if no mapping was necessary
        */

    if (None == name_mapping) { /** create a 1-1 mapping if nothing was done */
        name_mapping = {};
        utility.ForEach (alignments.GetSequenceNames (prefix+".nuc_data"), "_value_", "`&name_mapping`[_value_] = _value_");
    }

    utility.SetEnvVariable(utility.getGlobalValue ("terms.trees.data_for_neighbor_joining"),
                           nuc_data_info[utility.getGlobalValue("terms.data.datafilter")]);

    partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (nuc_data_info[utility.getGlobalValue("terms.data.partitions")], name_mapping);

    utility.SetEnvVariable(utility.getGlobalValue ("terms.trees.data_for_neighbor_joining"),
                           None);

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

    io.ReportProgressMessage ("", ">Loaded a multiple sequence alignment with **" + nuc_data_info[utility.getGlobalValue("terms.data.sequences")] + "** sequences, **" + nuc_data_info[utility.getGlobalValue("terms.data.sites")] + "** sites, and **" + partition_count + "** partitions from \`" + nuc_data_info[utility.getGlobalValue("terms.data.file")] + "\`");

    

    /** Input attribute to JSON **/
    json[utility.getGlobalValue("terms.json.input")] = {};
    (json[utility.getGlobalValue("terms.json.input")])[utility.getGlobalValue("terms.json.file")] =  nuc_data_info[utility.getGlobalValue("terms.data.file")];
    (json[utility.getGlobalValue("terms.json.input")])[utility.getGlobalValue("terms.json.sequences")] = nuc_data_info[utility.getGlobalValue("terms.data.sequences")];
    (json[utility.getGlobalValue("terms.json.input")])[utility.getGlobalValue("terms.json.sites")] = nuc_data_info[utility.getGlobalValue("terms.data.sites")];
    (json[utility.getGlobalValue("terms.json.input")])[utility.getGlobalValue("terms.json.partition_count")] = partition_count;

    // The trees should go into input as well and they should be w/ their branch lengths but ONLY if they have any.



     filter_specification = alignments.DefineFiltersForPartitions (partitions_and_trees, "`prefix`.nuc_data" , "`prefix`.filter.", nuc_data_info);
   

     store_tree_information ();
 
 }
 
function annotate_codon_info (codon_info, filter_name) {
    sample_size=codon_info[utility.getGlobalValue("terms.data.sites")]*codon_info[utility.getGlobalValue("terms.data.sequences")];


    codon_info[utility.getGlobalValue("terms.data.sample_size")] = sample_size;
    upper_prefix = prefix && 1; //uppercase the prefix for json name
    codon_info[utility.getGlobalValue("terms.json.json")] = codon_info[utility.getGlobalValue("terms.data.file")] + "."+upper_prefix+".json";

    name_mapping = codon_info[utility.getGlobalValue("terms.data.name_mapping")];

        /**
            will contain "mapped" -> "original" associations with sequence names; or null if no mapping was necessary
        */

    if (None == name_mapping) { /** create a 1-1 mapping if nothing was done */
        name_mapping = {};
        for (_value_; in; alignments.GetSequenceNames (filter_name)) {
            name_mapping[_value_] = _value_;
        }
    }
    
    codon_info[^"terms.data.name_mapping"] = name_mapping;

    // check for duplicates
    duplicate_sequences = codon_info[^"terms.data.sequences"] - alignments.HasDuplicateSequences (codon_info[^"terms.data.datafilter"],-1);
    if (duplicate_sequences > 0) {
      fprintf(stdout, "\n-------\n", io.FormatLongStringToWidth(
      ">[WARNING] '`codon_info[utility.getGlobalValue("terms.data.file")]`' contains " + duplicate_sequences + " duplicate " + io.SingularOrPlural (duplicate_sequences, 'sequence', 'sequences') + ".
      Identical sequences do not contribute any information to the analysis and only slow down computation.
      Please consider removing duplicate or 'nearly' duplicate sequences,
      e.g. using https://github.com/veg/hyphy-analyses/tree/master/remove-duplicates
      prior to running selection analyses", 72),
      "\n-------\n");
    }
}

function load_file (prefix) {

    settings = None;
    multiple_files = FALSE;
    blb = 1.0;

    if (Type (prefix) == "AssociativeList") {
        multiple_files  = prefix [utility.getGlobalValue("terms.multiple_files")];
        settings        = prefix [utility.getGlobalValue("terms.settings")];
        blb             = prefix [utility.getGlobalValue("terms.data.blb_subsample")];
        prefix          = prefix [utility.getGlobalValue("terms.prefix")];
    }

    if (multiple_files) {
        SetDialogPrompt ("Supply a list of files to include in the analysis (one per line)");
        codon_data_info = {};
        fscanf (PROMPT_FOR_FILE, "Lines", file_list_path);
        codon_data_info[utility.getGlobalValue("terms.data.file")] = ^"LAST_FILE_PATH";
        file_list = io.validate_a_list_of_files (file_list_path);
        file_count = utility.Array1D (file_list);
        io.CheckAssertion("`&file_count` > 1", "At least one valid filepath is required");
        partitions_and_trees = {};
        codon_data_info [utility.getGlobalValue("terms.data.sequences")] = {1, file_count};
        codon_data_info [utility.getGlobalValue("terms.data.sites")] = {1, file_count};
        codon_data_info [utility.getGlobalValue("terms.data.sample_size")] = 0;
        datasets = {};
        
        for (i = 0; i < file_count; i+=1) {
           filepath = file_list[i];
           datasets[i] = prefix+".codon_data_" + i;
           if (+i == 0) {
                codon_data_info [filepath] = 
                     alignments.LoadCodonDataFile(datasets[i],  prefix+".codon_filter_" + i , alignments.ReadCodonDataSetFromPath(prefix+".codon_data_" + i, filepath));
                (codon_data_info)[^"terms.code"]  = (codon_data_info[filepath])[^"terms.code"];
                (codon_data_info)[^"terms.stop_codons"]  = (codon_data_info[filepath])[^"terms.stop_codons"];
            } else {
                codon_data_info [filepath] = 
                     alignments.LoadCodonDataFile(datasets[i],  prefix+".codon_filter_" + i , alignments.ReadCodonDataSetFromPathGivenCode(prefix+".codon_data_" + i, filepath, (codon_data_info)[^"terms.code"] , (codon_data_info)[^"terms.stop_codons"] ));
            }
            annotate_codon_info ( codon_data_info [filepath], prefix+".codon_filter_" + i);
            partitions_and_trees + (trees.LoadAnnotatedTreeTopology.match_partitions ({{"file_" + i ,""}}, (codon_data_info [filepath])[^"terms.data.name_mapping"]))[0];
            (codon_data_info [utility.getGlobalValue("terms.data.sequences")])[+i] =  (codon_data_info[filepath]) [utility.getGlobalValue("terms.data.sequences")];
            (codon_data_info [utility.getGlobalValue("terms.data.sites")])[+i] =  (codon_data_info[filepath]) [utility.getGlobalValue("terms.data.sites")];
            codon_data_info [utility.getGlobalValue("terms.data.sample_size")] +=  (codon_data_info[filepath]) [utility.getGlobalValue("terms.data.sample_size")];
        }
        

    } else {
        datasets = prefix+".codon_data";
        codon_data_info = alignments.PromptForGeneticCodeAndAlignment(datasets, prefix+".codon_filter");
        if (blb < 1 && blb > 0) {
            console.log (blb);
            codon_data_info [utility.getGlobalValue("terms.data.blb_subsample")] = blb;
        }
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
        annotate_codon_info (codon_data_info, prefix+".codon_filter");
        utility.SetEnvVariable(utility.getGlobalValue ("terms.trees.data_for_neighbor_joining"),
                               codon_data_info[utility.getGlobalValue("terms.data.datafilter")]);
        partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (codon_data_info[utility.getGlobalValue("terms.data.partitions")], name_mapping);
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
                
        for (_key_, _value_; in; partitions_and_trees) {
            (partitions_and_trees[_key_])[utility.getGlobalValue("terms.data.filter_string")] = 
                    selection.io.adjust_partition_string (_value_[utility.getGlobalValue("terms.data.filter_string")], 3*codon_data_info[utility.getGlobalValue("terms.data.sites")]);
        }
     }   

    

    utility.SetEnvVariable(utility.getGlobalValue ("terms.trees.data_for_neighbor_joining"), None);
    partition_count = Abs (partitions_and_trees);

    if (multiple_files) {
        io.ReportProgressMessage ("", ">Loaded **`partition_count`** alignments from from \`" + codon_data_info[utility.getGlobalValue("terms.data.file")] + "\`");
        for (filepath, fileinfo; in; codon_data_info) {
            if (Type (fileinfo) == "AssociativeList") {
                if (utility.Has (fileinfo, utility.getGlobalValue("terms.data.sequences"), "Number")) {
                     io.ReportProgressMessage ("", "**\``filepath`\`** : " +  fileinfo[utility.getGlobalValue("terms.data.sequences")] + "** sequences, **" + fileinfo[utility.getGlobalValue("terms.data.sites")] + "** codons");
                }
            }
        }
    
    } else {
        io.ReportProgressMessage ("", ">Loaded a multiple sequence alignment with **" + codon_data_info[utility.getGlobalValue("terms.data.sequences")] + "** sequences, **" + codon_data_info[utility.getGlobalValue("terms.data.sites")] + "** codons, and **" + partition_count + "** partitions from \`" + codon_data_info[utility.getGlobalValue("terms.data.file")] + "\`");
    }

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

    
    filter_specification = alignments.DefineFiltersForPartitions (partitions_and_trees, datasets , "`prefix`.filter.", codon_data_info);
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
     
     if (None != name_mapping) {
         name_mapping_upper_case = {};
         for (i,n; in; name_mapping) {
             name_mapping_upper_case[i&&1] = n;
         }

         for (partition_index = 0; partition_index < partition_count; partition_index += 1) {
        
            local_name_mapping = {};
            for (l, label; in; (trees[partition_index])[utility.getGlobalValue ("terms.trees.partitioned")]) {
                if (label == ^"terms.tree_attributes.leaf") {
                    if (name_mapping / l) {
                        local_name_mapping [l] = name_mapping [l];
                    } else {
                        local_name_mapping [l] = name_mapping_upper_case [l];
                    }
                }
            }
        
    
            selection.io.json_store_branch_attribute(json, utility.getGlobalValue ("terms.original_name"), utility.getGlobalValue ("terms.json.node_label"), display_orders[utility.getGlobalValue ("terms.original_name")],
                                             partition_index,
                                             local_name_mapping);
        }
    }



}




//------------------------------------------------------------------------------

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
    
    KeywordArgument ("intermediate-fits", "Use/save parameter estimates from 'initial-guess' model fits to a JSON file (default is not to save)", "/dev/null");
    save_intermediate_fits = io.ReadFromOrCreate ("Use/Save parameter estimates from 'initial-guess' model fits", {});
    
    run_gtr = TRUE;
    
    if (None != save_intermediate_fits[^"terms.data.value"]) {
        if (utility.Has (save_intermediate_fits[^"terms.data.value"], "GTR", "AssociativeList")) {
            gtr_results = (save_intermediate_fits[^"terms.data.value"])["GTR"];
            run_gtr = FALSE;
        }        
    } else {
        save_intermediate_fits = None;
    }
    


    if (run_gtr) {
        utility.ToggleEnvVariable ("TREE_DO_NOT_AUTO_RENAME", TRUE);

        gtr_results = estimators.FitGTR(filter_names,
                                             trees,
                                             gtr_results);
                                             
         utility.ToggleEnvVariable ("TREE_DO_NOT_AUTO_RENAME", None);

          if (Type (save_intermediate_fits) == "AssociativeList") {                                             
               (save_intermediate_fits[^"terms.data.value"])["GTR"] = gtr_results;
               io.SpoolJSON (save_intermediate_fits[^"terms.data.value"],save_intermediate_fits[^"terms.data.file"]);     
        }                         
    }

    KeywordArgument ("kill-zero-lengths", "Automatically delete internal zero-length branches for computational efficiency (will not affect results otherwise)", "Yes");

    kill0 = io.SelectAnOption (
        {
            "Yes":"Automatically delete internal zero-length branches for computational efficiency (will not affect results otherwise)",
            "Constrain":"Keep zero-length branches, but constrain their values to 0",
            "No":"Keep all branches"
        },
        "Reduce zero-length branches");

    zero_branch_length_constrain = NULL;
    deleted_by_tree = NULL;
    
    if (kill0 == "Yes") {
        for (index, tree; in; trees) {
            deleted = {};
            
            if (^(prefix + ".selected_branches") / index) {
                deleted_test = utility.Array1D ((^(prefix + ".selected_branches"))[index]);
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
    } else {
        if (kill0 == "Constrain") {
            zero_branch_length_constrain = ^"namespace_tag"+".constrain_zero_branches";
            deleted_by_tree = {};
            for (index, tree; in; trees) {
                deleted = {};
                trees.KillZeroBranches (tree, (gtr_results[^"terms.branch_length"])[index], null, deleted);
 
                if (utility.Array1D (deleted)) {
                    io.ReportProgressMessageMD(prefix,  'selector', 'Marked ' + Abs(deleted) + ' zero-length internal branches to be constrained: \`' + Join (', ',utility.Values(deleted)) + '\`');
                }
                if (Abs (deleted)) { 
                   deleted_dict = {};
                   for (v; in; deleted) {
                    deleted_dict[v] = 1;
                   }   
                   deleted_by_tree[index] = deleted_dict;
                } else {
                    deleted_by_tree[index] = deleted;
                }
            }    
        }
    }


    io.ReportProgressMessageMD (prefix, "nuc-fit", "* " +
        selection.io.report_fit (gtr_results, 0, 3*(^"`prefix`.sample_size")));

    io.ReportProgressMessageMD (prefix, "nuc-fit", "* " +
        selection.io.report_fit_secondary_stats (gtr_results));
    
    

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
    
    gtr_bls_over_10 = 0;
    gtr_bl_sum = 0;
    
    for (partition_index = 0; partition_index < Abs(filter_specification); partition_index += 1) {
        selection.io.json_store_branch_attribute(json, utility.getGlobalValue ("terms.json.nucleotide_gtr"), utility.getGlobalValue ("terms.branch_length"), display_orders[terms.json.nucleotide_gtr],
                                         partition_index,
                                         selection.io.extract_branch_info((gtr_results[utility.getGlobalValue ("terms.branch_length")])[partition_index], "selection.io.branch.length"));
    
        for (bl; in; (gtr_results[utility.getGlobalValue ("terms.branch_length")])[partition_index]) {
            if (bl[^"terms.fit.MLE"] > 10.) {
                gtr_bls_over_100 += 1;
            }
        }
    }
    
    if (gtr_bls_over_10 > 0) {
      fprintf(stdout, "\n-------\n", io.FormatLongStringToWidth(
      ">[WARNING] Some of the branches (" + gtr_bls_over_10 + ") in this alignment appear to be 'infinitely' long. 
      While in *rare* cases this may be biologically reasonable, such overly long branches typically indicate either alignment issues (e.g. codon sequences aligned as nucleotides, improperly trimmed genes, etc),
      or lack of sequence homology (too distantly related species). Consider inspecting your alignment to confirm that it does not contain problematic sequences and that it is of suitable quality.", 72),
      "\n-------\n");
   }

}




/**
 * @name doPartitionedMG
 * Can only be used after including shared-load-file
 * @return
 */



//------------------------------------------------------------------------------

lfunction mg94mh.constrain2H (tree_name, node_name, model_description, partition) {
    if (utility.Has (model_description [utility.getGlobalValue ("terms.local")], utility.getGlobalValue ("terms.parameters.multiple_hit_rate"), "String")) {
        delta = (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.multiple_hit_rate")];
        g_delta = (model_description[utility.getGlobalValue ("terms.global")])[utility.getGlobalValue ("terms.parameters.multiple_hit_rate")];
        parameters.SetConstraint (tree_name + "." + node_name + "." + delta, g_delta, "");
    }
    return tree_name + "." + node_name + "." + delta;
}




//------------------------------------------------------------------------------

function constrain_2H (lf_id, components, data_filter, tree, model_map, initial_values, model_objects) {
    delta = prefix + ".delta_shared";
    parameters.DeclareGlobal(delta, None);
    for (k,i; in; model_objects) {
        model.generic.AddGlobal (i, delta, utility.getGlobalValue ("terms.parameters.multiple_hit_rate"));
    }
    parameter_set = estimators.TraverseLocalParameters (lf_id, model_objects, prefix + ".mg94mh.constrain2H");
    return 0;
}

//------------------------------------------------------------------------------

lfunction mg94mh.constrain3H (tree_name, node_name, model_description, partition) {
    if (utility.Has (model_description [utility.getGlobalValue ("terms.local")], utility.getGlobalValue ("terms.parameters.multiple_hit_rate"), "String")) {
        delta = (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.multiple_hit_rate")];
        g_delta = (model_description[utility.getGlobalValue ("terms.global")])[utility.getGlobalValue ("terms.parameters.multiple_hit_rate")];
        parameters.SetConstraint (tree_name + "." + node_name + "." + delta, g_delta, "");
    }
    if (utility.Has (model_description [utility.getGlobalValue ("terms.local")], utility.getGlobalValue ("terms.parameters.triple_hit_rate"), "String")) {
        psi = (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.triple_hit_rate")];
        g_psi = (model_description[utility.getGlobalValue ("terms.global")])[utility.getGlobalValue ("terms.parameters.triple_hit_rate")];
        parameters.SetConstraint (tree_name + "." + node_name + "." + psi, g_psi, "");
    }
    if (utility.Has (model_description [utility.getGlobalValue ("terms.local")], utility.getGlobalValue ("terms.parameters.triple_hit_rate_syn"), "String")) {
        psi = (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.triple_hit_rate_syn")];
        g_psi = (model_description[utility.getGlobalValue ("terms.global")])[utility.getGlobalValue ("terms.parameters.triple_hit_rate")];
        parameters.SetConstraint (tree_name + "." + node_name + "." + psi, g_psi, "");
    }
    return {"0" : tree_name + "." + node_name + "." + delta, "1" : tree_name + "." + node_name + "." + psi};
}




//------------------------------------------------------------------------------

function constrain_3H (lf_id, components, data_filter, tree, model_map, initial_values, model_objects) {
    delta = prefix + ".delta_shared";
    psi = prefix + ".psi_shared";
    parameters.DeclareGlobal(delta, None);
    parameters.DeclareGlobal(psi, None);

    for (k,i; in; model_objects) {
        model.generic.AddGlobal (i, delta, utility.getGlobalValue ("terms.parameters.multiple_hit_rate"));
        model.generic.AddGlobal (i, psi, utility.getGlobalValue ("terms.parameters.triple_hit_rate"));
    }
    parameter_set = estimators.TraverseLocalParameters (lf_id, model_objects, prefix + ".mg94mh.constrain3H");
    return 0;
}

/**
 * @name doPartitionedMG
 * Can only be used after including shared-load-file
 * @return
 */
function doPartitionedMG (prefix, keep_lf) {
    doPartitionedMGModel (prefix, keep_lf, "models.codon.MG_REV.ModelDescription", None);
}

/**
 * @name doPartitionedMGModel
 * Can only be used after including shared-load-file
 * @return
 */
function doPartitionedMGModel (prefix, keep_lf, model, constraint) {
    io.ReportProgressMessageMD ("`prefix`", "codon-fit", "Obtaining the global omega estimate based on relative GTR branch lengths and nucleotide substitution biases");


    /**
        Declare per-partition branch length scalers
        slac.scaler_prefix_0
        slac.scaler_prefix_1
        etc
    */
    scaler_variables = utility.PopulateDict (0, partition_count, "`prefix`.scaler_prefix + '_' + _k_", "_k_");
    
    for (_value_; in; scaler_variables) {
        parameters.DeclareGlobal(_value_, None);
        parameters.SetValue(_value_, 3);
    }

    
    run_mg94 = TRUE;
    
    if (Type (save_intermediate_fits) == "AssociativeList") {
        if (None != save_intermediate_fits[^"terms.data.value"]) {
            if (utility.Has (save_intermediate_fits[^"terms.data.value"], "MG94", "AssociativeList")) {
                partitioned_mg_results = (save_intermediate_fits[^"terms.data.value"])["MG94"];
                if (keep_lf) {
                    if (utility.Has (save_intermediate_fits[^"terms.data.value"], "MG94-LF", "String")) {
                        ExecuteCommands ((save_intermediate_fits[^"terms.data.value"])["MG94-LF"]);
                        run_mg94 = FALSE;
                    } 
                } else {
                    run_mg94 = FALSE;
                }
            }        
        }
    }
    
    if (run_mg94) {
    
        partitioned_mg_results = estimators.FitCodonModel (filter_names, trees, model, codon_data_info [utility.getGlobalValue("terms.code")], {
            utility.getGlobalValue("terms.run_options.model_type"): utility.getGlobalValue("terms.local"),
            utility.getGlobalValue("terms.run_options.proportional_branch_length_scaler"): scaler_variables,
            utility.getGlobalValue("terms.run_options.partitioned_omega"): selected_branches,
            utility.getGlobalValue("terms.run_options.retain_lf_object"): keep_lf,
            utility.getGlobalValue("terms.run_options.apply_user_constraints") : constraint
        }, gtr_results);
        if (Type (save_intermediate_fits) == "AssociativeList") {
            (save_intermediate_fits[^"terms.data.value"])["MG94"] = partitioned_mg_results;
            if (keep_lf) {
                Export (lfe, ^(partitioned_mg_results[^"terms.likelihood_function"]));
                (save_intermediate_fits[^"terms.data.value"])["MG94-LF"] = lfe;
            }
            io.SpoolJSON (save_intermediate_fits[^"terms.data.value"],save_intermediate_fits[^"terms.data.file"]);      
        }
    }
    

    io.ReportProgressMessageMD("`prefix`", "codon-fit", "* " + selection.io.report_fit (partitioned_mg_results, 0, (^"`prefix`.codon_data_info")[utility.getGlobalValue ("terms.data.sample_size")]));
    io.ReportProgressMessageMD (prefix, "nuc-fit", "* " +
        selection.io.report_fit_secondary_stats (partitioned_mg_results));

    global_dnds = partitioned_mg_results[utility.getGlobalValue("terms.global")];
    //selection.io.extract_global_MLE_re (partitioned_mg_results, ".");
    
    for (_desc_, _value_; in; global_dnds) {
        if (utility.Has (_value_, utility.getGlobalValue("terms.constraint"), "String")) {
            continue;
        }
        io.ReportProgressMessageMD ("`prefix`", "codon-fit", "* " + _desc_ + " = " + Format (_value_[utility.getGlobalValue("terms.fit.MLE")],8,4));
    }
    
    
}

//------------------------------------------------------------------------------


function constrain_zero_branches (lf_id, components, data_filter, tree, model_map, initial_values, model_objects) {
   SetParameter (DEFER_CONSTRAINT_APPLICATION, 1, 0);
   //console.log (">>>IN constrain_zero_branches");
   constrain_zero_branches.parameters = estimators.TraverseLocalParameters (lf_id, model_objects, ^"namespace_tag" + ".constrain_one_branch");
   constrain_zero_branches.parameters = +constrain_zero_branches.parameters;
   SetParameter (DEFER_CONSTRAINT_APPLICATION, 0, 0);
   return constrain_zero_branches.parameters;
}

//------------------------------------------------------------------------------

function constrain_one_branch (tree_name, node_name, model_description,i) {
   //console.log (">>>IN constrain_one_branch") 
   constrain_one_branch.tag = ^"namespace_tag" + ".deleted_by_tree";
    if (Type (^constrain_one_branch.tag) == "AssociativeList") {
        if ((^constrain_one_branch.tag)[i] / node_name) {
            if (utility.Has (model_description [utility.getGlobalValue ("terms.local")], utility.getGlobalValue ("terms.parameters.synonymous_rate"), "String")) {
                constrain_one_branch.t = (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.synonymous_rate")];
                //console.log (tree_name + "." + node_name + "." + constrain_one_branch.t + " => " + Eval (tree_name + "." + node_name + "." + constrain_one_branch.t));
                parameters.SetConstraint (tree_name + "." + node_name + "." + constrain_one_branch.t, "0", "");
                return 1;
            }
        }
    }
    return 0;
 }

