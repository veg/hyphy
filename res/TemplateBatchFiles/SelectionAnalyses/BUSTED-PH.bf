RequireVersion ("2.5.73");


LoadFunctionLibrary("libv3/all-terms.bf"); // must be loaded before CF3x4
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/tasks/genetic_code.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("modules/io_functions.ibf");
LoadFunctionLibrary("modules/selection_lib.ibf");


utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

terms.busted.background_test = "Background selection test results";
terms.busted.comparative_test = "Comparative selection test results";
terms.busted.constrained_null = "Background null";
terms.busted.constrained_same = "Same distribution";
terms.busted.busted_ph        = "BUSTED-PH tests";

bustedph.test_set = "";
bustedph.analysis_description = {

                               terms.io.info : "BUSTED-PH (phenotype) tests if episodic diversifying selection is associated with the set of designated (FG) branches.
The BUSTED model with selected options is fitted to the FG branches (test 1), then the rest of the branches (BG, test 2), and finally to all branches together (test 3). A BUSTED-PH detects selection associated with the phenotype is test 1 IS significant, test 2 is NOT significant, and test 3 IS significant when comparing it to the FG+BG model (i.e., selection regimes differ between the sets of branches)",
                               terms.io.version : "1.0",
                               terms.io.reference : "TBD",
                               terms.io.authors : "Sergei L Kosakovsky Pond",
                               terms.io.contact : "spond@temple.edu",
                               terms.settings: {},
                               terms.io.requirements : "in-frame codon alignment and a phylogenetic tree with at least two sets of labeled branches"
                              };

bustedph.json    = { terms.json.analysis: bustedph.analysis_description };


io.DisplayAnalysisBanner (bustedph.analysis_description);
bustedph.timers = {};

selection.io.startTimer (bustedph.timers, "Overall", 0);

KeywordArgument ("code",      "Which genetic code should be used", "Universal");
KeywordArgument ("alignment", "An in-frame codon alignment in one of the formats supported by HyPhy");
KeywordArgument ("tree",      "A phylogenetic tree (optionally annotated with {})", null, "Please select a tree file for the data:");
KeywordArgument ("branches",  "Branches to test for association (FG)");
GetString (bustedph.kwargs, KWARGS, 0);


namespace bustedph {
    ExecuteAFile ("modules/shared-load-file.bf");
    load_file ({
                utility.getGlobalValue("terms.prefix"): "bustedph", 
                utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "bustedph.select_branches"}
               });
}


bustedph.pass_through = {"--alignment" : bustedph.codon_data_info[terms.data.file],
                         "--branches" : bustedph.test_set};
                         
if (Type ((bustedph.trees[0])[bustedph.codon_data_info[terms.data.file]]) == "String") {
    bustedph.pass_through  ["--tree"] = (bustedph.trees[0])[bustedph.codon_data_info[terms.data.file]];
}

for (i, k; in; bustedph.kwargs) {
    if (bustedph.pass_through / ("--" + i) == FALSE) {
        bustedph.pass_through ["--" + i] = k;
    }
}


console.log ('
===============================================================================================
[BUSTED-PH PHASE 1] Testing for selection on foreground branches by running the BUSTED analysis
===============================================================================================
');

LoadFunctionLibrary ("BUSTED.bf", bustedph.pass_through);

console.log ('
================================================================
[BUSTED-PH PHASE 2] Testing for selection on background branches
================================================================
');

if (busted.error_sink) {
    busted.run_background_test = busted.inferred_background_distribution_raw [busted.rate_classes-1][0] > 1 && busted.inferred_background_distribution_raw [busted.rate_classes-1][1] > 0;
} else {
    busted.run_background_test = busted.inferred_background_distribution [busted.rate_classes-1][0] > 1 && busted.inferred_background_distribution [busted.rate_classes-1][1] > 0;
}

if (!busted.run_test) {
    io.ReportProgressMessageMD ("BUSTED", "Results", "No evidence for episodic diversifying positive selection under the unconstrained model on BACKGROUND BRANCHES, skipping constrained model fitting");
    busted.json [terms.busted.background_test ] = busted.ComputeLRT (0, 0);
} else {
    parameters.RemoveConstraint (model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted.rate_classes)));
    estimators.ApplyExistingEstimates (busted.full_model[terms.likelihood_function], busted.model_object_map,  busted.full_model, {});
    parameters.SetConstraint (model.generic.GetGlobalParameter (busted.background.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted.rate_classes)), terms.parameters.one, terms.global);
    
    busted.background_null_results = estimators.FitExistingLF (busted.full_model[terms.likelihood_function], busted.model_object_map);
    
    io.ReportProgressMessageMD ("BUSTED", "test", "* " + selection.io.report_fit (busted.background_null_results, 9, busted.codon_data_info[terms.data.sample_size]));
    busted.LRT = busted.ComputeLRT (busted.full_model[terms.fit.log_likelihood], busted.background_null_results[terms.fit.log_likelihood]);
    busted.json [terms.busted.background_test] = busted.LRT;

 
    utility.ForEachPair (busted.filter_specification, "_key_", "_value_",
        'selection.io.json_store_branch_attribute(busted.json, terms.busted.constrained_null, terms.branch_length, 1,
                                                 _key_,
                                                 selection.io.extract_branch_info((busted.background_null_results[terms.branch_length])[_key_], "selection.io.branch.length"));');

    
    io.ReportProgressMessageMD("BUSTED", "test", "* For *test* branches under the null (no dN/dS > 1 model), the following rate distribution for branch-site combinations was inferred");
    busted.inferred_test_distribution_raw = parameters.GetStickBreakingDistribution (busted.distribution);
    busted.inferred_test_distribution = busted.inferred_test_distribution_raw  % 0;

    busted.distribution_for_json = {busted.FG : utility.Map (utility.Range (busted.rate_classes, 0, 1),
                                                         "_index_",
                                                         "{terms.json.omega_ratio : busted.inferred_test_distribution_raw [_index_][0],
                                                           terms.json.proportion : busted.inferred_test_distribution_raw [_index_][1]}")};

    busted.report_multi_hit  (busted.null_results, busted.distribution_for_json, "MultiHit", "null-mh",busted.branch_length_string, busted.model_parameters);
    busted.null_distro_raw = parameters.GetStickBreakingDistribution (busted.distribution);
    
    if (busted.error_sink) {
        selection.io.report_dnds_with_sink (busted.null_distro_raw % 0, busted.null_distro_raw[0][0], busted.null_distro_raw[0][1]);
    } else {
        selection.io.report_dnds (busted.null_distro_raw % 0);
    }


    if (busted.has_background) {
        busted.inferred_background_distribution_raw = parameters.GetStickBreakingDistribution (busted.background_distribution);
        io.ReportProgressMessageMD("BUSTED", "test", "* For *background* branches under the null (no dN/dS > 1 model), the following rate distribution for branch-site combinations was inferred");
        selection.io.report_dnds (busted.inferred_background_distribution_raw % 0);
        busted.distribution_for_json [busted.BG] = utility.Map (utility.Range (busted.rate_classes, 0, 1),
                                                             "_index_",
                                                             "{terms.json.omega_ratio : busted.inferred_background_distribution_raw [_index_][0],
                                                               terms.json.proportion : busted.inferred_background_distribution_raw [_index_][1]}");
    }

    busted.report_srv_fit (FALSE);
    
    console.log ("----\n## Branch-site unrestricted statistical test of episodic diversification on background branches [BUSTED]");
    console.log ( "Likelihood ratio test for episodic diversifying positive selection, **p = " + Format ((busted.json [terms.busted.background_test])[terms.p_value], 8, 4) + "**.");


    selection.io.json_store_lf (busted.json,
                            "Constrained model on background",
                            busted.background_null_results[terms.fit.log_likelihood],
                            busted.background_null_results[terms.parameters] + 9 , // +9 comes from CF3x4
                            busted.codon_data_info[terms.data.sample_size],
                            busted.distribution_for_json,
                            busted.display_orders[busted.constrained] + 1);
                            
    
}

console.log ('
=====================================================================================================
[BUSTED-PH PHASE 3] Testing for the equality of omega distributions between foreground and background
=====================================================================================================
');

   estimators.ApplyExistingEstimates (busted.full_model[terms.likelihood_function], busted.model_object_map,  busted.full_model, {});
   

busted.same_df = 0;

for (i,p; in; busted.distribution[terms.parameters.rates]) {
     parameters.SetConstraint ((busted.background_distribution[terms.parameters.rates])[i], p, terms.global);
     busted.same_df += 1;
}

for (i,p; in; busted.distribution[terms.parameters.weights]) {
     parameters.SetConstraint ((busted.background_distribution[terms.parameters.weights])[i], p, terms.global);
     busted.same_df += 1;
}


busted.background_same_results = estimators.FitExistingLF (busted.full_model[terms.likelihood_function], busted.model_object_map);
busted.LRT = math.DoLRT (busted.background_same_results[terms.fit.log_likelihood], busted.full_model[terms.fit.log_likelihood], busted.same_df);
busted.json [terms.busted.comparative_test] = busted.LRT;

io.ReportProgressMessageMD ("BUSTED", "test", "* " + selection.io.report_fit (busted.background_same_results, 9, busted.codon_data_info[terms.data.sample_size]));


utility.ForEachPair (busted.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(busted.json, terms.busted.constrained_same, terms.branch_length, 1,
                                             _key_,
                                             selection.io.extract_branch_info((busted.background_same_results[terms.branch_length])[_key_], "selection.io.branch.length"));');


io.ReportProgressMessageMD("BUSTED", "same", "* For single distributions model, the following rate distribution for branch-site combinations was inferred");
busted.inferred_test_distribution_raw = parameters.GetStickBreakingDistribution (busted.distribution);
busted.inferred_test_distribution = busted.inferred_test_distribution_raw  % 0;

busted.distribution_for_json = {busted.FG : utility.Map (utility.Range (busted.rate_classes, 0, 1),
                                                     "_index_",
                                                     "{terms.json.omega_ratio : busted.inferred_test_distribution_raw [_index_][0],
                                                       terms.json.proportion : busted.inferred_test_distribution_raw [_index_][1]}")};

busted.report_multi_hit  (busted.background_same_results, busted.distribution_for_json, "MultiHit", "null-mh",busted.branch_length_string, busted.model_parameters);
busted.null_distro_raw = parameters.GetStickBreakingDistribution (busted.distribution);

if (busted.error_sink) {
    selection.io.report_dnds_with_sink (busted.null_distro_raw % 0, busted.null_distro_raw[0][0], busted.null_distro_raw[0][1]);
} else {
    selection.io.report_dnds (busted.null_distro_raw % 0);
}

busted.report_srv_fit (FALSE);

console.log ("----\n## Comparison of one vs two dN/dS distributions");
console.log ( "Likelihood ratio test for non-identity of distributions, **p = " + Format ((busted.json [terms.busted.comparative_test])[terms.p_value], 8, 4) + "**.");

selection.io.json_store_lf (busted.json,
                        "Same distributions model",
                        busted.background_same_results[terms.fit.log_likelihood],
                        busted.background_same_results[terms.parameters] + 9 , // +9 comes from CF3x4
                        busted.codon_data_info[terms.data.sample_size],
                        busted.distribution_for_json,
                        busted.display_orders[busted.constrained] + 2);


console.log ('
==================================================
[BUSTED-PH PHASE 4] Performing association testing
==================================================
');

busted.p_value = 0.05;

bustedph.test_results = {
    terms.json.uncorrected_pvalue : {
        'FG' : (busted.json [terms.json.test_results])[terms.p_value],
        'BG' : (busted.json [terms.busted.background_test])[terms.p_value],
        'Comparative' : (busted.json [terms.busted.comparative_test])[terms.p_value]
    },
    'Level' : busted.p_value
};

bustedph.test_results [terms.json.corrected_pvalue ] = math.HolmBonferroniCorrection (bustedph.test_results[terms.json.uncorrected_pvalue]);

console.log ("----\n## Branch-site unrestricted statistical test of episodic diversification and association with phenotype/trait [BUSTED-PH]");
console.log ( "Likelihood ratio test for episodic diversifying positive selection on test branches , **p = " + Format ((bustedph.test_results [terms.json.corrected_pvalue ])['FG'], 8, 4) + "**.");
console.log ( "Likelihood ratio test for episodic diversifying positive selection on background branches , **p = " + Format ((bustedph.test_results [terms.json.corrected_pvalue ])['BG'], 8, 4) + "**.");
console.log ( "Likelihood ratio test for differences in distributions between **test** and **background** , **p = " + Format ((bustedph.test_results [terms.json.corrected_pvalue ])['Comparative'], 8, 4) + "**.");

console.log ("\n\n## Analysis summary (p = " + bustedph.test_results['Level'] + ")");

bustedph.summary = "";

if ((bustedph.test_results [terms.json.corrected_pvalue ])['FG'] <= busted.p_value ) {
    if ((bustedph.test_results [terms.json.corrected_pvalue ])['BG'] > busted.p_value) {
        if ((bustedph.test_results [terms.json.corrected_pvalue ])['Comparative'] <= busted.p_value) {
            bustedph.summary = ("Selection is associated with the phenotype / trait");
        } else {
            bustedph.summary = ("Selection is associated with the phenotype / trait, but there is no significant difference between test and background branches in terms of selective pressure");
        }
    } else {
        console.log ("Selection is acting on the branches with the phenotype / trait, but is **also** acting on background branches.");
        if ((bustedph.test_results [terms.json.corrected_pvalue ])['Comparative'] <= busted.p_value) {
            bustedph.summary = ("There is a significant difference between test and background branches in terms of selective pressure");
        } else {
            bustedph.summary = ("There is no significant difference between test and background branches in terms of selective pressure");
        }
        
    }
} else {
    bustedph.summary = ("There is **no evidence** of episodic diversifying selection on test branches; selection is not associated with phenotype/trait");
}

console.log (bustedph.summary);
bustedph.test_results ["Summary"] = bustedph.summary;
busted.json["BUSTED-PH"] = bustedph.test_results;

(bustedph.json[terms.json.analysis])[terms.settings] = (busted.json[terms.json.analysis])[terms.settings];

busted.json[terms.json.analysis] = bustedph.json[terms.json.analysis];

io.SpoolJSON (busted.json, busted.codon_data_info [terms.json.json]);



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


lfunction bustedph.select_branches(partition_info) {

 
    tree_count = utility.Array1D (partition_info);
    //io.CheckAssertion("`&tree_count` == 1", "BUSTED-PH only works on a single partition dataset");
    available_models = {};

    return_set         = {};
    tree_configuration = {};

    for (i,p; in; partition_info) {
        tree_for_analysis = p[utility.getGlobalValue("terms.data.tree")];
        tree_configuration [i] = {};
        if (+i == 0) {
            for (_value_; in; tree_for_analysis[utility.getGlobalValue("terms.trees.model_map")]) {
                available_models [_value_] = 1;
            }
        } else {
            partition_models = {};
            for (_value_; in; tree_for_analysis[utility.getGlobalValue("terms.trees.model_map")]) {
                partition_models [_value_] = 1;
            }
            for (k, v; in; partition_models) {
                available_models[k] += 1;
            }
        }
    }
    
    list_models = {};
    for (i,p; in; available_models) {
        if (p == tree_count) {
            list_models[i] = 1;
        }
    }
    
      
    list_models   = utility.sortStrings(utility.Keys(list_models)); // get keys
    option_count  = utility.Array1D (list_models);
    can_run_group_mode = (option_count == Abs (available_models));
    
    io.CheckAssertion	("`&option_count` >= 2", "BUSTED-PH requires at least one designated set of branches in the tree.");
    
	

    selectTheseForTesting = {
        option_count, 2
    };

    for (k = 0; k < option_count; k += 1) {
        if (list_models[k] != "") {
            selectTheseForTesting[k][0] = list_models[k];
            selectTheseForTesting[k][1] = "Set " + list_models[k] + " with " + available_models[list_models[k]] + " branches";
        } else {
            selectTheseForTesting[k][0] = "Unlabeled branches";
            selectTheseForTesting[k][1] = "Set of " + available_models[list_models[k]] + " unlabeled branches";
        }
    }

    ChoiceList(testSet, "Choose the set of branches to use as the _foreground_ (FG) set", 1, NO_SKIP, selectTheseForTesting);
    io.CheckAssertion ("`&testSet` >= 0", "User cancelled branch selection; analysis terminating");


    tag_test = selectTheseForTesting [testSet][0];
    if (tag_test == "Unlabeled branches") {
        tag_test = "";
    }
    
    ^"bustedph.test_set" = tag_test;


    for (i,p; in; partition_info) {		
        tree_for_analysis = p[utility.getGlobalValue("terms.data.tree")];	 
        for (_key_,_value_; in; tree_for_analysis[utility.getGlobalValue("terms.trees.model_map")]) {
            if (tag_test == _value_ ) {
                (tree_configuration[i])[_key_] = utility.getGlobalValue('terms.tree_attributes.test');
            } else {
                if (tag_reference == _value_ ) {
                    (tree_configuration[i])[_key_] = utility.getGlobalValue('terms.tree_attributes.background');
                } 
            }
        }
    }


    return tree_configuration;
}
