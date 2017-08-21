RequireVersion ("2.31");



LoadFunctionLibrary("libv3/all-terms.bf"); // must be loaded before CF3x4

LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("CF3x4");
LoadFunctionLibrary("TreeTools");


// namespace 'utility' for convenience functions
LoadFunctionLibrary("libv3/UtilityFunctions.bf");

// namespace 'io' for interactive/datamonkey i/o functions
LoadFunctionLibrary("libv3/IOFunctions.bf");

LoadFunctionLibrary ("libv3/models/codon/MG_REV.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("libv3/tasks/estimators.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("libv3/tasks/alignments.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("libv3/tasks/trees.bf");

LoadFunctionLibrary("modules/io_functions.ibf");
LoadFunctionLibrary("modules/selection_lib.ibf");
LoadFunctionLibrary("libv3/models/codon/BS_REL.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");


utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);
utility.SetEnvVariable ("ASSUME_REVERSIBLE_MODELS", TRUE);

/*------------------------------------------------------------------------------*/

relax.json    = { terms.json.input: {},
                  terms.json.fits : {},
                  terms.json.timers : {},
                  terms.json.test_results : {}
                  };

relax.relaxation_parameter        = "relax.K";
relax.rate_classes     = 3;
relax.MG94 = "MG94xREV with separate rates for branch sets";

terms.relax.k          = "relaxation or intensification parameter";
terms.relax.k_range    = {
        terms.lower_bound: "0",
        terms.upper_bound: "100"
    };

relax.p_threshold = 0.05;

relax.test_branches_name = "Test";
relax.reference_branches_name = "Reference";
relax.unclassified_branches_name = "Unclassified";

/*------------------------------------------------------------------------------*/


relax.analysis_description = {
                               terms.io.info : "RELAX (a random effects test of selection relaxation) uses a random effects branch-site model framework to test whether a set of 'Test' branches evolves under relaxed selection relative to a set of 'Reference' branches (R), as measured by the relaxation parameter (K).",
                               terms.io.version : "2.0",
                               terms.io.reference : "RELAX: Detecting Relaxed Selection in a Phylogenetic Framework (2015). Mol Biol Evol 32 (3): 820-832",
                               terms.io.authors : "Sergei L Kosakovsky Pond, Ben Murrell, Steven Weaver and Temple iGEM / UCSD viral evolution group",
                               terms.io.contact : "spond@temple.edu",
                               terms.io.requirements : "in-frame codon alignment and a phylogenetic tree, with at least two groups of branches defined using the {} notation (one group can be defined as all unlabeled branches)"
                              };

io.DisplayAnalysisBanner ( relax.analysis_description );

selection.io.startTimer (relax.json [terms.json.timers], "Overall", 0);

namespace relax {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ({utility.getGlobalValue("terms.prefix"): "relax", utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "relax.select_branches"}});
}


io.ReportProgressMessageMD('RELAX',  'selector', 'Branch sets for RELAX analysis');

relax.has_unclassified = FALSE;


utility.ForEachPair (relax.selected_branches, "_partition_", "_selection_",
    "_selection_ = utility.Filter (_selection_, '_value_', '_value_ == relax.test_branches_name');
     io.ReportProgressMessageMD('RELAX',  'selector', '* Selected ' + Abs(_selection_) + ' branches as the _test_ set: \\\`' + Join (', ',utility.Keys(_selection_)) + '\\\`')");

utility.ForEachPair (relax.selected_branches, "_partition_", "_selection_",
    "_selection_ = utility.Filter (_selection_, '_value_', '_value_ == relax.reference_branches_name');
     io.ReportProgressMessageMD('RELAX',  'selector', '* Selected ' + Abs(_selection_) + ' branches as the _reference_ set: \\\`' + Join (', ',utility.Keys(_selection_)) + '\\\`')");

utility.ForEachPair (relax.selected_branches, "_partition_", "_selection_",
    "_selection_ = utility.Filter (_selection_, '_value_', '_value_ == relax.unclassified_branches_name');
     relax.has_unclassified = Abs(_selection_) > 0;
     if (relax.has_unclassified) {
        io.ReportProgressMessageMD('RELAX',  'selector', '* ' + Abs(_selection_) + ' branches are in the unclassified (nuisance) set: \\\`' + Join (', ',utility.Keys(_selection_)) + '\\\`')
     }");


relax.model_set = io.SelectAnOption ({
                                        {"All", "[Default] Fit descriptive models AND run the relax test (4 models)"}
                                        {"Minimal", "Run only the RELAX test (2 models)"}
                                    }, "RELAX analysis type");

selection.io.startTimer (relax.json [terms.json.timers], "Preliminary model fitting", 1);


namespace relax {
    doGTR ("relax");
}

estimators.fixSubsetOfEstimates(relax.gtr_results, relax.gtr_results[terms.global]);

namespace relax {
    scaler_prefix = "relax.scaler";
    doPartitionedMG ("relax", FALSE);
}


io.ReportProgressMessageMD ("RELAX", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");

relax.final_partitioned_mg_results = estimators.FitMGREV (relax.filter_names, relax.trees, relax.codon_data_info [terms.code], {
    terms.run_options.model_type: terms.local,
    terms.run_options.partitioned_omega: relax.selected_branches,
}, relax.partitioned_mg_results);


io.ReportProgressMessageMD("RELAX", "codon-refit", "* " + selection.io.report_fit (relax.final_partitioned_mg_results, 0, relax.codon_data_info[terms.data.sample_size]));

relax.global_dnds  = selection.io.extract_global_MLE_re (relax.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio);
relax.report_dnds = {};

utility.ForEach (relax.global_dnds, "_value_", '
    io.ReportProgressMessageMD ("RELAX", "codon-refit", "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));
    relax.report_dnds [(regexp.FindSubexpressions (_value_[terms.description], "^" + terms.parameters.omega_ratio + ".+\\*(.+)\\*$"))[1]] = {"0" : {terms.json.omega_ratio : _value_[terms.fit.MLE], terms.json.proportion : 1}};
');

selection.io.json_store_branch_attribute(relax.json, terms.original_name, terms.json.node_label, 0,
                                         0,
                                         relax.name_mapping);


selection.io.json_store_lf_spool (relax.codon_data_info [terms.json.json], relax.json,
                            relax.MG94,
                            relax.final_partitioned_mg_results[terms.fit.log_likelihood],
                            relax.final_partitioned_mg_results[terms.parameters] ,
                            math.GetIC (relax.final_partitioned_mg_results[terms.fit.log_likelihood], relax.final_partitioned_mg_results[terms.parameters], relax.codon_data_info[terms.data.sample_size]),
                            relax.report_dnds);


selection.io.stopTimer (relax.json [terms.json.timers], "Preliminary model fitting");


if (relax.model_set == "All") { // run all the models

    relax.ge.bsrel_model =  model.generic.DefineMixtureModel("relax.BS_REL.ModelDescription",
            "relax.ge", {
                "0": parameters.Quote(terms.local),
                "1": relax.codon_data_info[terms.code],
                "2": parameters.Quote (relax.rate_classes) // the number of rate classes
            },
            relax.filter_names,
            None);

    for (relax.i = 1; relax.i < relax.rate_classes; relax.i += 1) {
        parameters.SetRange (model.generic.GetGlobalParameter (relax.ge.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)), terms.range_almost_01);
    }
    parameters.SetRange (model.generic.GetGlobalParameter (relax.ge.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.rate_classes)), terms.range_gte1);

    relax.model_object_map = { "relax.ge" :       relax.ge.bsrel_model };

    io.ReportProgressMessageMD ("RELAX", "gd", "Fitting the general descriptive (separate k per branch) model");
    selection.io.startTimer (relax.json [terms.json.timers], "General descriptive model fitting", 2);

    relax.ge_guess = relax.DistributionGuess(utility.Map (selection.io.extract_global_MLE_re (relax.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio + ".+test.+"), "_value_",
            "_value_[terms.fit.MLE]"));

    relax.distribution = models.codon.BS_REL.ExtractMixtureDistribution(relax.ge.bsrel_model);
    parameters.SetStickBreakingDistribution (relax.distribution, relax.ge_guess);

    relax.general_descriptive.fit =  estimators.FitLF (relax.filter_names,
                                    relax.trees,
                                    { "0" : {"DEFAULT" : "relax.ge"}},
                                    relax.final_partitioned_mg_results,
                                    relax.model_object_map,
                                    {
                                        terms.run_options.apply_user_constraints: "relax.init.k"
                                    });



    selection.io.stopTimer (relax.json [terms.json.timers], "General descriptive model fitting");

    io.ReportProgressMessageMD("RELAX", "ge", "* " + selection.io.report_fit (relax.general_descriptive.fit, 9, relax.codon_data_info[terms.data.sample_size]));
    io.ReportProgressMessageMD("RELAX", "ge", "* The following baseline rate distribution for branch-site combinations was inferred");
    relax.inferred_ge_distribution = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistributionFromFit (relax.ge.bsrel_model, relax.general_descriptive.fit)) % 0;
    selection.io.report_dnds (relax.inferred_ge_distribution);
    relax.distribution_for_json = {'Shared' : utility.Map (utility.Range (relax.rate_classes, 0, 1),
                                                         "_index_",
                                                         "{terms.json.omega_ratio : relax.inferred_ge_distribution [_index_][0],
                                                           terms.json.proportion  : relax.inferred_ge_distribution [_index_][1]}")
                                   };
    selection.io.json_store_lf_spool (relax.codon_data_info [terms.json.json], relax.json,
                                "General descriptive",
                                relax.general_descriptive.fit[terms.fit.log_likelihood],
                                relax.general_descriptive.fit[terms.parameters] + 9 , // +9 comes from CF3x4
                                math.GetIC (relax.general_descriptive.fit[terms.fit.log_likelihood], relax.general_descriptive.fit[terms.parameters] + 9, relax.codon_data_info[terms.data.sample_size]),
                                relax.distribution_for_json
                            );

    selection.io.json_store_branch_attribute(relax.json, "General descriptive", terms.branch_length, 1,
                                                 0,
                                                 selection.io.extract_branch_info((relax.general_descriptive.fit[terms.branch_length])[0], "selection.io.branch.length"));

    relax.k_estimates = selection.io.extract_branch_info((relax.general_descriptive.fit[terms.branch_length])[0], "relax.extract.k");

    relax.k_stats = math.GatherDescriptiveStats (utility.Map (utility.Values (relax.k_estimates), "_value_", "0+_value_"));

    io.ReportProgressMessageMD("RELAX", "ge", "* Branch-level `terms.relax.k` distribution has mean " + Format (relax.k_stats[terms.math.mean], 5,2) + ", median " +
                                                 Format (relax.k_stats[terms.math.median], 5,2) + ", and 95% of the weight in " + Format (relax.k_stats[terms.math._2.5], 5,2) + " - " + Format (relax.k_stats[terms.math._97.5], 5,2));


    selection.io.json_store_branch_attribute(relax.json, "k (general descriptive)", terms.json.branch_label, 1,
                                                 0,
                                                 relax.k_estimates);

} else {
    relax.general_descriptive.fit = relax.final_partitioned_mg_results;
}

/* now fit the two main models for RELAX */

io.ReportProgressMessageMD ("RELAX", "alt", "Fitting the alternative model to test K != 1");

selection.io.startTimer (relax.json [terms.json.timers], "RELAX alternative model fitting", 3);

relax.test.bsrel_model =  model.generic.DefineMixtureModel("models.codon.BS_REL.ModelDescription",
        "relax.test", {
            "0": parameters.Quote(terms.global),
            "1": relax.codon_data_info[terms.code],
            "2": parameters.Quote (relax.rate_classes) // the number of rate classes
        },
        relax.filter_names,
        None);



relax.reference.bsrel_model =  model.generic.DefineMixtureModel("models.codon.BS_REL.ModelDescription",
        "relax.reference", {
            "0": parameters.Quote(terms.global),
            "1": relax.codon_data_info[terms.code],
            "2": parameters.Quote (relax.rate_classes) // the number of rate classes
        },
        relax.filter_names,
        None);

relax.bound_weights = models.BindGlobalParameters ({"0" : relax.reference.bsrel_model, "1" : relax.test.bsrel_model}, terms.mixture.mixture_aux_weight + ".+");
models.BindGlobalParameters ({"0" : relax.test.bsrel_model, "1" : relax.reference.bsrel_model}, terms.nucleotideRate("[ACGT]","[ACGT]"));

parameters.DeclareGlobalWithRanges (relax.relaxation_parameter, 1, 0, 50);
model.generic.AddGlobal (relax.test.bsrel_model, relax.relaxation_parameter, terms.relax.k);

for (relax.i = 1; relax.i < relax.rate_classes; relax.i += 1) {
    parameters.SetRange (model.generic.GetGlobalParameter (relax.reference.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)), terms.range01);
    parameters.SetRange (model.generic.GetGlobalParameter (relax.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)), terms.range01);
    parameters.SetConstraint (model.generic.GetGlobalParameter (relax.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)),
                              model.generic.GetGlobalParameter (relax.reference.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)) + "^" + relax.relaxation_parameter,
                              terms.global);
}
parameters.SetRange (model.generic.GetGlobalParameter (relax.reference.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.rate_classes)), terms.range_gte1);
parameters.SetRange (model.generic.GetGlobalParameter (relax.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.rate_classes)), terms.range_gte1);
parameters.SetConstraint (model.generic.GetGlobalParameter (relax.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)),
                          model.generic.GetGlobalParameter (relax.reference.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)) + "^" + relax.relaxation_parameter,
                          terms.global);

relax.model_map = {
                    "relax.test" : utility.Filter (relax.selected_branches[0], '_value_', '_value_ == relax.test_branches_name'),
                    "relax.reference" : utility.Filter (relax.selected_branches[0], '_value_', '_value_ == relax.reference_branches_name')
                  };


// constrain the proportions to be the same

relax.model_object_map = { "relax.reference" : relax.reference.bsrel_model,
                            "relax.test" :       relax.test.bsrel_model };

if (relax.model_set != "All") {
    relax.ge_guess = relax.DistributionGuess(utility.Map (selection.io.extract_global_MLE_re (relax.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio + ".+test.+"), "_value_",
            "_value_[terms.fit.MLE]"));

    relax.distribution = models.codon.BS_REL.ExtractMixtureDistribution(relax.reference.bsrel_model);
    parameters.SetStickBreakingDistribution (relax.distribution, relax.ge_guess);
}

if (relax.has_unclassified) {
    relax.unclassified.bsrel_model =  model.generic.DefineMixtureModel("models.codon.BS_REL.ModelDescription",
        "relax.unclassified", {
            "0": parameters.Quote(terms.global),
            "1": relax.codon_data_info[terms.code],
            "2": parameters.Quote (relax.rate_classes) // the number of rate classes
        },
        relax.filter_names,
        None);

    for (relax.i = 1; relax.i < relax.rate_classes-1; relax.i += 1) {
        parameters.SetRange (model.generic.GetGlobalParameter (relax.unclassified.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)), terms.range01);
    }

    parameters.SetRange (model.generic.GetGlobalParameter (relax.unclassified.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.rate_classes)), terms.range_gte1);

    relax.model_object_map ["relax.unclassified"] = relax.unclassified.bsrel_model;
    relax.model_map ["relax.unclassified"] = utility.Filter (relax.selected_branches[0], '_value_', '_value_ == relax.unclassified_branches_name');
    models.BindGlobalParameters ({"0" : relax.unclassified.bsrel_model, "1" : relax.reference.bsrel_model}, terms.nucleotideRate("[ACGT]","[ACGT]"));
}

relax.alternative_model.fit =  estimators.FitLF (relax.filter_names, relax.trees, { "0" : relax.model_map}, relax.general_descriptive.fit, relax.model_object_map, {terms.run_options.retain_lf_object: TRUE});

io.ReportProgressMessageMD("RELAX", "alt", "* " + selection.io.report_fit (relax.alternative_model.fit, 9, relax.codon_data_info[terms.data.sample_size]));


relax.fitted.K = estimators.GetGlobalMLE (relax.alternative_model.fit,terms.relax.k);
io.ReportProgressMessageMD("RELAX", "alt", "* Relaxation/intensification parameter (K) = " + Format(relax.fitted.K,8,2));
io.ReportProgressMessageMD("RELAX", "alt", "* The following rate distribution was inferred for **test** branches");
relax.inferred_distribution = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistribution (relax.test.bsrel_model)) % 0;
selection.io.report_dnds (relax.inferred_distribution);

io.ReportProgressMessageMD("RELAX", "alt", "* The following rate distribution was inferred for **reference** branches");
relax.inferred_distribution_ref = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistribution (relax.reference.bsrel_model)) % 0;
selection.io.report_dnds (relax.inferred_distribution_ref);

relax.distribution_for_json = {relax.test_branches_name : utility.Map (utility.Range (relax.rate_classes, 0, 1),
                                                     "_index_",
                                                     "{terms.json.omega_ratio : relax.inferred_distribution [_index_][0],
                                                       terms.json.proportion  : relax.inferred_distribution [_index_][1]}"),

                                relax.reference_branches_name : utility.Map (utility.Range (relax.rate_classes, 0, 1),
                                                     "_index_",
                                                     "{terms.json.omega_ratio : relax.inferred_distribution_ref [_index_][0],
                                                       terms.json.proportion  : relax.inferred_distribution_ref [_index_][1]}")
                               };

selection.io.json_store_lf_spool (relax.codon_data_info [terms.json.json], relax.json,
                            "RELAX alternative",
                            relax.alternative_model.fit[terms.fit.log_likelihood],
                            relax.alternative_model.fit[terms.parameters] + 9 , // +9 comes from CF3x4
                            relax.codon_data_info[terms.data.sample_size],
                            relax.distribution_for_json
                        );

selection.io.json_store_branch_attribute(relax.json, "RELAX alternative", terms.branch_length, 2,
                                             0,
                                             selection.io.extract_branch_info((relax.alternative_model.fit[terms.branch_length])[0], "selection.io.branch.length"));

selection.io.stopTimer (relax.json [terms.json.timers], "RELAX alternative model fitting");

// NULL MODEL

selection.io.startTimer (relax.json [terms.json.timers], "RELAX null model fitting", 4);

io.ReportProgressMessageMD ("RELAX", "null", "Fitting the null (K := 1) model");
parameters.SetConstraint (model.generic.GetGlobalParameter (relax.test.bsrel_model , terms.relax.k), terms.parameters.one, terms.global);
relax.null_model.fit = estimators.FitExistingLF (relax.alternative_model.fit[terms.likelihood_function], relax.model_object_map);
io.ReportProgressMessageMD ("RELAX", "null", "* " + selection.io.report_fit (relax.null_model.fit, 9, relax.codon_data_info[terms.data.sample_size]));
relax.LRT = math.DoLRT (relax.null_model.fit[terms.fit.log_likelihood], relax.alternative_model.fit[terms.fit.log_likelihood], 1);


io.ReportProgressMessageMD("RELAX", "null", "* The following rate distribution for test/reference branches was inferred");
relax.inferred_distribution = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistributionFromFit (relax.test.bsrel_model, relax.null_model.fit)) % 0;
selection.io.report_dnds (relax.inferred_distribution);

relax.distribution_for_json = {relax.test_branches_name : utility.Map (utility.Range (relax.rate_classes, 0, 1),
                                                     "_index_",
                                                     "{terms.json.omega_ratio : relax.inferred_distribution [_index_][0],
                                                       terms.json.proportion  : relax.inferred_distribution [_index_][1]}")};

relax.distribution_for_json   [relax.reference_branches_name] =   relax.distribution_for_json   [relax.test_branches_name];

selection.io.json_store_lf_spool (relax.codon_data_info [terms.json.json], relax.json,
                            "RELAX null",
                            relax.null_model.fit[terms.fit.log_likelihood],
                            relax.null_model.fit[terms.parameters] + 9 , // +9 comes from CF3x4
                            relax.codon_data_info[terms.data.sample_size],
                            relax.distribution_for_json
                        );

selection.io.json_store_branch_attribute(relax.json, "RELAX null", terms.branch_length, 3,
                                             0,
                                             selection.io.extract_branch_info((relax.null_model.fit[terms.branch_length])[0], "selection.io.branch.length"));


console.log ("----\n## Test for relaxation (or intensification) of selection [RELAX]");
console.log ( "Likelihood ratio test **p = " + Format (relax.LRT[terms.p_value], 8, 4) + "**.");

if (relax.LRT[terms.p_value] <= relax.p_threshold) {
    if (relax.fitted.K > 1) {
        console.log (">Evidence for relaxation of selection among **test** branches _relative_ to the **reference** branches at P<="+ relax.p_threshold);
    } else {
        console.log (">Evidence for intensification of selection among **test** branches _relative_ to the **reference** branches at P<="+ relax.p_threshold);
    }
} else {
    console.log (">No significant evidence for relaxation (or intensification) of selection among **test** branches _relative_ to the **reference** branches at P<="+ relax.p_threshold);
}

relax.json [terms.json.test_results] = relax.LRT;
(relax.json [terms.json.test_results])[terms.relax.k] = relax.fitted.K;

console.log ("----\n");

selection.io.stopTimer (relax.json [terms.json.timers], "RELAX null model fitting");

if (relax.model_set == "All") {
    selection.io.startTimer (relax.json [terms.json.timers], "RELAX partitioned exploratory", 5);

    io.ReportProgressMessageMD ("RELAX", "pe", "Fitting the partitioned exploratory model (separate distributions for *test* and *reference* branches)");
    parameters.RemoveConstraint (utility.Keys (relax.bound_weights));
    for (relax.i = 1; relax.i < relax.rate_classes; relax.i += 1) {
        parameters.RemoveConstraint (model.generic.GetGlobalParameter (relax.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)));
    }
    relax.pe.fit = estimators.FitExistingLF (relax.alternative_model.fit[terms.likelihood_function], relax.model_object_map);
    io.ReportProgressMessageMD ("RELAX", "pe", "* " + selection.io.report_fit (relax.pe.fit, 9, relax.codon_data_info[terms.data.sample_size]));
    io.ReportProgressMessageMD ("RELAX", "pe", "* The following rate distribution was inferred for *test* branches ");
    relax.test.inferred_distribution = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistribution(relax.test.bsrel_model)) % 0;
    selection.io.report_dnds (relax.test.inferred_distribution);
    io.ReportProgressMessageMD("RELAX", "pe", "* The following rate distribution was inferred for *reference* branches ");
    relax.reference.inferred_distribution =  parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistribution(relax.reference.bsrel_model)) % 0;
    selection.io.report_dnds (relax.reference.inferred_distribution);

    relax.distribution_for_json = {relax.test_branches_name : utility.Map (utility.Range (relax.rate_classes, 0, 1),
                                                         "_index_",
                                                         "{terms.json.omega_ratio : relax.test.inferred_distribution [_index_][0],
                                                           terms.json.proportion  : relax.test.inferred_distribution [_index_][1]}"),

                                    relax.reference_branches_name : utility.Map (utility.Range (relax.rate_classes, 0, 1),
                                                         "_index_",
                                                         "{terms.json.omega_ratio : relax.reference.inferred_distribution [_index_][0],
                                                           terms.json.proportion  : relax.reference.inferred_distribution [_index_][1]}")
                                   };

    selection.io.json_store_lf_spool (relax.codon_data_info [terms.json.json], relax.json,
                                "RELAX partitioned exploratory",
                                relax.pe.fit[terms.fit.log_likelihood],
                                relax.pe.fit[terms.parameters] + 9 , // +9 comes from CF3x4
                                relax.codon_data_info[terms.data.sample_size],
                                relax.distribution_for_json
                            );

    selection.io.json_store_branch_attribute(relax.json, "RELAX partitioned exploratory", terms.branch_length, 4,
                                                 0,
                                                 selection.io.extract_branch_info((relax.pe.fit[terms.branch_length])[0], "selection.io.branch.length"));


    selection.io.stopTimer (relax.json [terms.json.timers], "RELAX partitioned exploratory");
}

selection.io.stopTimer (relax.json [terms.json.timers], "Overall");
io.SpoolJSON (relax.json, relax.codon_data_info [terms.json.json]);

return relax.json;


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

lfunction relax.extract.k(branch_info) {
    return (branch_info[utility.getGlobalValue("terms.relax.k")])[utility.getGlobalValue("terms.fit.MLE")];
}

//------------------------------------------------------------------------------

lfunction relax.set.k (tree_name, node_name, model_description) {
    if (utility.Has (model_description [utility.getGlobalValue ("terms.local")], utility.getGlobalValue ("terms.relax.k"), "String")) {
        k = (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.relax.k")];
        parameters.SetValue (tree_name + "." + node_name + "." + k, 1);
        parameters.SetRange (tree_name + "." + node_name + "." + k, utility.getGlobalValue ("terms.relax.k_range"));
    }
    return tree_name + "." + node_name + "." + k;
}

//------------------------------------------------------------------------------

lfunction relax.init.k (lf_id, components, data_filter, tree, model_map, initial_values, model_objects) {
    parameter_set = estimators.TraverseLocalParameters (lf_id, model_objects, "relax.set.k");
    parameters.SetConstraint (model.generic.GetGlobalParameter (utility.getGlobalValue("relax.ge.bsrel_model") , terms.AddCategory (utility.getGlobalValue("terms.parameters.omega_ratio"),2)), utility.getGlobalValue("terms.parameters.one"), utility.getGlobalValue("terms.global"));
    /*parameters.SetConstraint (model.generic.GetGlobalParameter (utility.getGlobalValue("relax.ge.bsrel_model") , terms.AddCategory (utility.getGlobalValue("terms.parameters.omega_ratio"),utility.getGlobalValue ("relax.rate_classes"))),
                             "1/(" +
                                Join ("*", utility.Map (
                                    utility.Range (utility.getGlobalValue ("relax.rate_classes") - 1, 1, 1),
                                    "_value_",
                                    'model.generic.GetGlobalParameter (relax.ge.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,_value_))'
                                    ))
                             + ")",
                            "global");*/

    return 0;
}

//------------------------------------------------------------------------------

lfunction relax.BS_REL.ModelDescription (type, code, components) {
    model = models.codon.BS_REL.ModelDescription(utility.getGlobalValue ('terms.global'), code, components);
    model [utility.getGlobalValue("terms.model.defineQ")] = "relax.BS_REL._DefineQ";
    return model;
}

//------------------------------------------------------------------------------


lfunction relax.DistributionGuess (mean) {
    guess = {{0.05,0.7}{0.25,0.2}{10,0.1}};

    norm = + guess[-1][1];
    guess_mean = 1/(+(guess [-1][0] $ guess [-1][1]))/norm;
    return guess["_MATRIX_ELEMENT_VALUE_*(guess_mean*(_MATRIX_ELEMENT_COLUMN_==0)+(_MATRIX_ELEMENT_COLUMN_==1)*(1/norm))"];
}

//------------------------------------------------------------------------------

lfunction relax.BS_REL._GenerateRate (fromChar, toChar, namespace, model_type, _tt, alpha, alpha_term, beta, beta_term, omega, omega_term) {

    p = {};
    diff = models.codon.diff(fromChar, toChar);

    if (None != diff) {
        p[model_type] = {};
        p[utility.getGlobalValue("terms.global")] = {};

        if (diff[utility.getGlobalValue("terms.diff.from")] > diff[utility.getGlobalValue("terms.diff.to")]) {
            nuc_rate = "theta_" + diff[utility.getGlobalValue("terms.diff.to")] + diff[utility.getGlobalValue("terms.diff.from")];
        } else {
            nuc_rate = "theta_" + diff[utility.getGlobalValue("terms.diff.from")] + diff[utility.getGlobalValue("terms.diff.to")];
        }
        nuc_rate = parameters.ApplyNameSpace(nuc_rate, namespace);
        (p[utility.getGlobalValue("terms.global")])[terms.nucleotideRate(diff[utility.getGlobalValue("terms.diff.from")], diff[utility.getGlobalValue("terms.diff.to")])] = nuc_rate;

        if (_tt[fromChar] != _tt[toChar]) {
            if (model_type == utility.getGlobalValue("terms.global")) {
                aa_rate = parameters.ApplyNameSpace(omega, namespace);
                (p[model_type])[omega_term] = aa_rate;
                utility.EnsureKey (p, utility.getGlobalValue("terms.local"));
                 (p[utility.getGlobalValue("terms.local")])[utility.getGlobalValue ("terms.relax.k")] = "k";
                 aa_rate += "^k";
            } else {
                aa_rate = beta;
                (p[model_type])[beta_term] = aa_rate;
            }
            p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate + "*" + aa_rate;
        } else {
            if (model_type == utility.getGlobalValue("terms.local")) {
                (p[model_type])[alpha_term] = alpha;
                p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate + "*" + alpha;
            } else {
                p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate;
            }
        }
    }


    return p;
}

//------------------------------------------------------------------------------

lfunction relax.BS_REL._DefineQ(bs_rel, namespace) {
    rate_matrices = {};

    bs_rel [utility.getGlobalValue("terms.model.q_ij")] = &rate_generator;
    bs_rel [utility.getGlobalValue("terms.mixture.mixture_components")] = {};

    _aux = parameters.GenerateSequentialNames (namespace + ".bsrel_mixture_aux", bs_rel[utility.getGlobalValue("terms.model.components")] - 1, "_");
    _wts = parameters.helper.stick_breaking (_aux, None);
    mixture = {};

    for (component = 1; component <= bs_rel[utility.getGlobalValue("terms.model.components")]; component += 1) {
       key = "component_" + component;
       ExecuteCommands ("
        function rate_generator (fromChar, toChar, namespace, model_type, _tt) {
           return relax.BS_REL._GenerateRate (fromChar, toChar, namespace, model_type, _tt,
                'alpha', utility.getGlobalValue('terms.parameters.synonymous_rate'),
                'beta_`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.nonsynonymous_rate'), component),
                'omega`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component));
            }"
       );

       if ( component < bs_rel[utility.getGlobalValue("terms.model.components")]) {
            model.generic.AddGlobal ( bs_rel, _aux[component-1], terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), component ));
            parameters.DeclareGlobalWithRanges (_aux[component-1], 0.5, 0, 1);
       }
       models.codon.generic.DefineQMatrix(bs_rel, namespace);
       rate_matrices [key] = bs_rel[utility.getGlobalValue("terms.model.rate_matrix")];
       (bs_rel [^'terms.mixture.mixture_components'])[key] = _wts [component-1];
    }


    bs_rel[utility.getGlobalValue("terms.model.rate_matrix")] = rate_matrices;
    parameters.SetConstraint(((bs_rel[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.global")])[terms.nucleotideRate("A", "G")], "1", "");
    return bs_rel;
}



//------------------------------------------------------------------------------
lfunction relax.select_branches(partition_info) {

    io.CheckAssertion("utility.Array1D (`&partition_info`) == 1", "RELAX only works on a single partition dataset");
    available_models = {};
    branch_set = {};


    tree_for_analysis = (partition_info[0])[utility.getGlobalValue("terms.data.tree")];
    utility.ForEach (tree_for_analysis[utility.getGlobalValue("terms.trees.model_map")], "_value_", "`&available_models`[_value_] += 1");
    list_models   = utility.Keys   (available_models); // get keys
    branch_counts = utility.Values (available_models);
    option_count  = Abs (available_models);

    io.CheckAssertion("`&option_count` >= 2", "RELAX requires at least one designated set of branches in the tree.");

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

    ChoiceList(testSet, "Choose the set of branches to use as the _test_ set", 1, NO_SKIP, selectTheseForTesting);
    io.CheckAssertion ("`&testSet` >= 0", "User cancelled branch selection; analysis terminating");
    if (option_count > 2) {
        ChoiceList(referenceSet, "Choose the set of branches to use as the _reference_ set", 1, testSet, selectTheseForTesting);
        io.CheckAssertion ("`&referenceSet` >= 0", "User cancelled branch selection; analysis terminating");
    } else {
        referenceSet = 1-testSet;
    }

    return_set = {};

    tree_configuration = {};
    tree_for_analysis = (partition_info[0])[utility.getGlobalValue("terms.data.tree")];

    tag_test = selectTheseForTesting [testSet][0];
    if (tag_test == "Unlabeled branches") {
        tag_test = "";
    }
    tag_reference = selectTheseForTesting [referenceSet][0];
    if (tag_reference == "Unlabeled branches") {
        tag_reference = "";
    }

    utility.ForEachPair (tree_for_analysis[utility.getGlobalValue("terms.trees.model_map")], "_key_", "_value_", "
        if (`&tag_test` == _value_ ) {
            `&tree_configuration`[_key_] = utility.getGlobalValue('relax.test_branches_name');
        } else {
            if (`&tag_reference` == _value_ ) {
                `&tree_configuration`[_key_] = utility.getGlobalValue('relax.reference_branches_name');
            } else {
                `&tree_configuration`[_key_] = utility.getGlobalValue('relax.unclassified_branches_name');
            }
        }
    ");

    return_set + tree_configuration;
    return return_set;
}
