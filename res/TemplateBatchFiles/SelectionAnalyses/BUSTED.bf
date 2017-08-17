RequireVersion ("2.31");


LoadFunctionLibrary("libv3/all-terms.bf"); // must be loaded before CF3x4
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/models/codon.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/genetic_code.bf");
LoadFunctionLibrary("modules/io_functions.ibf");
LoadFunctionLibrary("modules/selection_lib.ibf");
LoadFunctionLibrary("libv3/models/codon/BS_REL.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);


busted.analysis_description = {terms.io.info : "BUSTED (branch-site unrestricted statistical test of episodic diversification) uses a random effects branch-site model fitted jointly to all or a subset of tree branches in order to test for alignment-wide evidence of episodic diversifying selection. Assuming there is evidence of positive selection (i.e. there is an omega > 1), BUSTED will also perform a quick evidence-ratio style analysis to explore which individual sites may have been subject to selection.",
                           terms.io.version : "1.2",
                           terms.io.reference : "*Gene-wide identification of episodic selection*, Mol Biol Evol. 32(5):1365-71",
                           terms.io.authors : "Sergei L Kosakovsky Pond",
                           terms.io.contact : "spond@temple.edu",
                           terms.io.requirements : "in-frame codon alignment and a phylogenetic tree (optionally annotated with {})"
                          };

io.DisplayAnalysisBanner (busted.analysis_description);

busted.FG = "Test";
busted.BG = "Background";
busted.background = "background";
busted.unconstrained = "unconstrained";
busted.constrained = "constrained";
busted.optimized_null = "optimized null";
busted.MG94 = "MG94xREV with separate rates for test/background branches";

busted.json.background = busted.background;
busted.json.site_logl  = "Site Log Likelihood";
busted.json.evidence_ratios  = "Evidence Ratios";
busted.rate_classes = 3;

busted.json    = { terms.json.analysis: busted.analysis_description,
                   terms.json.input: {},
                   busted.json.background: {},
                   terms.json.fits : {},
                   terms.json.timers : {},
                   busted.json.site_logl : {},
                   busted.json.evidence_ratios: {},
                   busted.json.site_logl : {}
                  };






selection.io.startTimer (busted.json [terms.json.timers], "Overall", 0);

namespace busted {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ("busted");
}


io.ReportProgressMessageMD('BUSTED',  'selector', 'Branches to test for selection in the BUSTED analysis');

utility.ForEachPair (busted.selected_branches, "_partition_", "_selection_",
    "_selection_ = utility.Filter (_selection_, '_value_', '_value_ == terms.tree_attributes.test');
     io.ReportProgressMessageMD('BUSTED',  'selector', '* Selected ' + Abs(_selection_) + ' branches to test in the BUSTED analysis: \\\`' + Join (', ',utility.Keys(_selection_)) + '\\\`')");

// check to see if there are any background branches

busted.has_background = FALSE;

utility.ForEachPair (busted.selected_branches, "_partition_", "_selection_",
    "_selection_ = utility.Filter (_selection_, '_value_', '_value_ != terms.tree_attributes.test');
     if (utility.Array1D (_selection_)) { busted.has_background = TRUE;} ");

busted.json[busted.json.background] =  busted.has_background;

selection.io.startTimer (busted.json [terms.json.timers], "Preliminary model fitting", 1);

namespace busted {
    doGTR ("busted");
}


estimators.fixSubsetOfEstimates(busted.gtr_results, busted.gtr_results[terms.global]);

namespace busted {
    scaler_prefix = "busted.scaler";
    doPartitionedMG ("busted", FALSE);
}


io.ReportProgressMessageMD ("BUSTED", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");

busted.final_partitioned_mg_results = estimators.FitMGREV (busted.filter_names, busted.trees, busted.codon_data_info [terms.code], {
    terms.run_options.model_type: terms.local,
    terms.run_options.partitioned_omega: busted.selected_branches,
}, busted.partitioned_mg_results);


io.ReportProgressMessageMD("BUSTED", "codon-refit", "* " + selection.io.report_fit (busted.final_partitioned_mg_results, 0, busted.codon_data_info[terms.data.sample_size]));
busted.global_dnds = selection.io.extract_global_MLE_re (busted.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio);

utility.ForEach (busted.global_dnds, "_value_", 'io.ReportProgressMessageMD ("BUSTED", "codon-refit", "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));');

selection.io.stopTimer (busted.json [terms.json.timers], "Preliminary model fitting");

utility.ForEachPair (busted.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(busted.json, busted.MG94, terms.branch_length, 2,
                                             _key_,
                                             selection.io.extract_branch_info((busted.final_partitioned_mg_results[terms.branch_length])[_key_], "selection.io.branch.length"));');


busted.test.bsrel_model =  model.generic.DefineMixtureModel("models.codon.BS_REL.ModelDescription",
        "busted.test", {
            "0": parameters.Quote(terms.global),
            "1": busted.codon_data_info[terms.code],
            "2": parameters.Quote (busted.rate_classes) // the number of rate classes
        },
        busted.filter_names,
        None);


busted.background.bsrel_model =  model.generic.DefineMixtureModel("models.codon.BS_REL.ModelDescription",
        "busted.background", {
            "0": parameters.Quote(terms.global),
            "1": busted.codon_data_info[terms.code],
            "2": parameters.Quote (busted.rate_classes) // the number of rate classes
        },
        busted.filter_names,
        None);


models.BindGlobalParameters ({"0" : busted.test.bsrel_model, "1" : busted.background.bsrel_model}, terms.nucleotideRate("[ACGT]","[ACGT]"));

busted.test_guess = busted.DistributionGuess(utility.Map (selection.io.extract_global_MLE_re (busted.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio + ".+test.+"), "_value_",
            "_value_[terms.fit.MLE]"));


busted.distribution = models.codon.BS_REL.ExtractMixtureDistribution(busted.test.bsrel_model);
parameters.SetStickBreakingDistribution (busted.distribution, busted.test_guess);


if (busted.has_background) {
    busted.model_object_map = { "busted.background" : busted.background.bsrel_model,
                                "busted.test" :       busted.test.bsrel_model };
    busted.background_distribution = models.codon.BS_REL.ExtractMixtureDistribution(busted.background.bsrel_model);
    busted.background_guess = busted.DistributionGuess(utility.Map (selection.io.extract_global_MLE_re (busted.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio + ".+background.+"), "_value_",
            "_value_[terms.fit.MLE]"));
    parameters.SetStickBreakingDistribution (busted.background_distribution, busted.background_guess);
} else {
    busted.model_object_map = { "busted.test" :       busted.test.bsrel_model };
}

// set up parameter constraints

for (busted.i = 1; busted.i < busted.rate_classes; busted.i += 1) {
    parameters.SetRange (model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted.i)), terms.range01);
    parameters.SetRange (model.generic.GetGlobalParameter (busted.background.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted.i)), terms.range01);
}

parameters.SetRange (model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted.rate_classes)), terms.range_gte1);

busted.model_map = {};

for (busted.partition_index = 0; busted.partition_index < busted.partition_count; busted.partition_index += 1) {

    selection.io.json_store_branch_attribute(busted.json, terms.original_name, terms.json.node_label, 0,
                                             busted.partition_index,
                                             busted.name_mapping);

    busted.model_map + { "busted.test" : utility.Filter (busted.selected_branches[busted.partition_index], '_value_', '_value_ == terms.tree_attributes.test'),
					     "busted.background" : utility.Filter (busted.selected_branches[busted.partition_index], '_value_', '_value_ != terms.tree_attributes.test')};
}

utility.SetEnvVariable ("ASSUME_REVERSIBLE_MODELS", TRUE);

selection.io.startTimer (busted.json [terms.json.timers], "Unconstrained BUSTED model fitting", 2);

io.ReportProgressMessageMD ("BUSTED", "main", "Performing the full (dN/dS > 1 allowed) branch-site model fit");
busted.full_model =  estimators.FitLF (busted.filter_names, busted.trees, busted.model_map, busted.final_partitioned_mg_results, busted.model_object_map, {"retain-lf-object": TRUE});
io.ReportProgressMessageMD("BUSTED", "main", "* " + selection.io.report_fit (busted.full_model, 9, busted.codon_data_info[terms.data.sample_size]));
io.ReportProgressMessageMD("BUSTED", "main", "* For *test* branches, the following rate distribution for branch-site combinations was inferred");

selection.io.stopTimer (busted.json [terms.json.timers], "Unconstrained BUSTED model fitting");

busted.inferred_test_distribution = parameters.GetStickBreakingDistribution (busted.distribution) % 0;
selection.io.report_dnds (busted.inferred_test_distribution);


busted.distribution_for_json = {busted.FG : utility.Map (utility.Range (busted.rate_classes, 0, 1),
                                                         "_index_",
                                                         "{terms.json.omega_ratio : busted.inferred_test_distribution [_index_][0],
                                                           terms.json.proportion : busted.inferred_test_distribution [_index_][1]}")
                                };


if (busted.has_background) {
    io.ReportProgressMessageMD("BUSTED", "main", "* For *background* branches, the following rate distribution for branch-site combinations was inferred");
    busted.inferred_background_distribution = parameters.GetStickBreakingDistribution (busted.background_distribution) % 0;
    selection.io.report_dnds (busted.inferred_background_distribution);
    busted.distribution_for_json [busted.BG] = utility.Map (utility.Range (busted.rate_classes, 0, 1),
                                                         "_index_",
                                                         "{terms.json.omega_ratio : busted.inferred_background_distribution [_index_][0],
                                                           terms.json.proportion : busted.inferred_background_distribution [_index_][1]}");
}

selection.io.json_store_lf_spool (busted.codon_data_info [terms.json.json], busted.json,
                            "Unconstrained model",
                            busted.full_model[terms.fit.log_likelihood],
                            busted.full_model[terms.parameters] + 9 , // +9 comes from CF3x4
                            busted.codon_data_info[terms.data.sample_size],
                            busted.distribution_for_json);


busted.run_test = busted.inferred_test_distribution [busted.rate_classes-1][0] > 1 && busted.inferred_test_distribution [busted.rate_classes-1][1] > 0;

utility.ForEachPair (busted.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(busted.json, busted.unconstrained, terms.branch_length, 0,
                                             _key_,
                                             selection.io.extract_branch_info((busted.full_model[terms.branch_length])[_key_], "selection.io.branch.length"));');


(busted.json [busted.json.site_logl])[busted.unconstrained] = busted.ComputeSiteLikelihoods (busted.full_model[terms.likelihood_function]);

if (!busted.run_test) {
    io.ReportProgressMessageMD ("BUSTED", "Results", "No evidence for episodic diversifying positive selection under the unconstrained model, skipping constrained model fitting");
    busted.json [terms.json.test_results] = busted.ComputeLRT (0, 0);
} else {
    selection.io.startTimer (busted.json [terms.json.timers], "Constrained BUSTED model fitting", 3);


    io.ReportProgressMessageMD ("BUSTED", "test", "Performing the constrained (dN/dS > 1 not allowed) model fit");
    parameters.SetConstraint (model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted.rate_classes)), terms.parameters.one, terms.global);
    (busted.json [busted.json.site_logl])[busted.constrained] = busted.ComputeSiteLikelihoods (busted.full_model[terms.likelihood_function]);
    busted.null_results = estimators.FitExistingLF (busted.full_model[terms.likelihood_function], busted.model_object_map);
    (busted.json [busted.json.site_logl])[busted.optimized_null] = busted.ComputeSiteLikelihoods (busted.full_model[terms.likelihood_function]);
    io.ReportProgressMessageMD ("BUSTED", "test", "* " + selection.io.report_fit (busted.null_results, 9, busted.codon_data_info[terms.data.sample_size]));
    busted.LRT = busted.ComputeLRT (busted.full_model[terms.fit.log_likelihood], busted.null_results[terms.fit.log_likelihood]);
    busted.json [terms.json.test_results] = busted.LRT;


    selection.io.stopTimer (busted.json [terms.json.timers], "Constrained BUSTED model fitting");

    utility.ForEachPair (busted.filter_specification, "_key_", "_value_",
        'selection.io.json_store_branch_attribute(busted.json, busted.constrained, terms.branch_length, 1,
                                                 _key_,
                                                 selection.io.extract_branch_info((busted.null_results[terms.branch_length])[_key_], "selection.io.branch.length"));');

    io.ReportProgressMessageMD("BUSTED", "test", "* For *test* branches under the null (no dN/dS > 1 model), the following rate distribution for branch-site combinations was inferred");
    busted.inferred_test_distribution = parameters.GetStickBreakingDistribution (busted.distribution) % 0;
    selection.io.report_dnds (parameters.GetStickBreakingDistribution (busted.distribution) % 0);

    busted.distribution_for_json = {busted.FG : utility.Map (utility.Range (busted.rate_classes, 0, 1),
                                                         "_index_",
                                                         "{terms.json.omega_ratio : busted.inferred_test_distribution [_index_][0],
                                                           terms.json.proportion : busted.inferred_test_distribution [_index_][1]}")};

    if (busted.has_background) {
        busted.inferred_background_distribution = parameters.GetStickBreakingDistribution (busted.background_distribution) % 0;
        //selection.io.report_dnds (busted.inferred_background_distribution);
        busted.distribution_for_json [busted.BG] = utility.Map (utility.Range (busted.rate_classes, 0, 1),
                                                             "_index_",
                                                             "{terms.json.omega_ratio : busted.inferred_background_distribution [_index_][0],
                                                               terms.json.proportion : busted.inferred_background_distribution [_index_][1]}");
    }


    selection.io.json_store_lf_spool (busted.codon_data_info [terms.json.json], busted.json,
                            "Constrained model",
                            busted.null_results[terms.fit.log_likelihood],
                            busted.null_results[terms.parameters] + 9 , // +9 comes from CF3x4
                            busted.codon_data_info[terms.data.sample_size],
                            busted.distribution_for_json);

    (busted.json [busted.json.evidence_ratios])[busted.constrained] = busted.EvidenceRatios ( (busted.json [busted.json.site_logl])[busted.unconstrained],  (busted.json [busted.json.site_logl])[busted.constrained]);
    (busted.json [busted.json.evidence_ratios ])[busted.optimized_null] = busted.EvidenceRatios ( (busted.json [busted.json.site_logl])[busted.unconstrained],  (busted.json [busted.json.site_logl])[busted.optimized_null]);
}


console.log ("----\n## Branch-site unrestricted statistical test of episodic diversification [BUSTED]");
console.log ( "Likelihood ratio test for episodic diversifying positive selection, **p = " + Format ((busted.json [terms.json.test_results])[terms.p_value], 8, 4) + "**.");

selection.io.stopTimer (busted.json [terms.json.timers], "Overall");

io.SpoolJSON (busted.json, busted.codon_data_info [terms.json.json]);

return busted.json;

/// HELPERS

//------------------------------------------------------------------------------
lfunction busted.ComputeSiteLikelihoods (id) {
    ConstructCategoryMatrix (sl, ^id, SITE_LOG_LIKELIHOODS);
   return sl;
}
//------------------------------------------------------------------------------

function busted.ComputeLRT (ha, h0) {
    return {terms.LRT : 2*(ha-h0),
            terms.p_value : 1-CChi2 (2*(ha-h0),2)};
}

//------------------------------------------------------------------------------
lfunction busted.EvidenceRatios (ha, h0) {
    return ha["Exp(_MATRIX_ELEMENT_VALUE_-h0[_MATRIX_ELEMENT_COLUMN_])"];
}
//------------------------------------------------------------------------------

lfunction busted.DistributionGuess (mean) {
    /*guess = {{Random (0,0.25),Random (0.5,0.7)}
             {Random (0.25,0.6),Random (0.1,0.2)}
             {Random (2,20),Random (0.01, 0.1)}};
    */

    guess = {{0,0.7}{0.25,0.2}{10,0.1}};

    norm = + guess[-1][1];
    guess_mean = 1/(+(guess [-1][0] $ guess [-1][1]))/norm;
    return guess["_MATRIX_ELEMENT_VALUE_*(guess_mean*(_MATRIX_ELEMENT_COLUMN_==0)+(_MATRIX_ELEMENT_COLUMN_==1)*(1/norm))"];
}


