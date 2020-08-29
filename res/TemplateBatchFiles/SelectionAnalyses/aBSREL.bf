RequireVersion ("2.5.8");

LoadFunctionLibrary("libv3/all-terms.bf"); // must be loaded before CF3x4


LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary ("libv3/models/codon/MG_REV.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");

LoadFunctionLibrary("modules/io_functions.ibf");
LoadFunctionLibrary("modules/selection_lib.ibf");
LoadFunctionLibrary("libv3/models/codon/BS_REL.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");
LoadFunctionLibrary("libv3/models/codon/MG_REV_MH.bf");
LoadFunctionLibrary("libv3/models/codon/MG_REV_TRIP.bf");

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);
utility.SetEnvVariable ("ASSUME_REVERSIBLE_MODELS", TRUE);
utility.SetEnvVariable ("USE_MEMORY_SAVING_DATA_STRUCTURES", 1e8);


/*------------------------------------------------------------------------------*/

absrel.max_rate_classes  = 5;
absrel.p_threshold = 0.05;

absrel.MG94                 = "MG94xREV with separate omega for each branch";
absrel.baseline_mg94xrev    = "Baseline MG94xREV";
absrel.baseline_omega_ratio = "Baseline MG94xREV omega ratio";
absrel.full_adaptive_model  = "Full adaptive model";
absrel.full_adaptive_model.ES  = "Full adaptive model (synonymous subs/site)";
absrel.full_adaptive_model.EN  = "Full adaptive model (non-synonymous subs/site)";
absrel.rate_classes         = "Rate classes";
absrel.per_branch_omega     = "Per-branch omega";
absrel.per_branch_delta     = "Per-branch delta";
absrel.per_branch_psi       = "Per-branch psi";

absrel.display_orders = {terms.original_name: -1,
                         terms.json.nucleotide_gtr: 0,
                              absrel.baseline_mg94xrev: 1,
                              absrel.baseline_omega_ratio: 1,
                              absrel.full_adaptive_model: 2,
                              absrel.rate_classes: 2,
                              terms.json.rate_distribution: 3,
                              terms.LRT: 4,
                              terms.json.uncorrected_pvalue: 5,
                              terms.json.corrected_pvalue: 6,
                              terms.parameters.multiple_hit_rate: 7,
                              terms.parameters.triple_hit_rate : 8,
                              absrel.full_adaptive_model.ES: 9,
                              absrel.full_adaptive_model.EN: 10
                             };



/*------------------------------------------------------------------------------*/


absrel.analysis_description = {terms.io.info : "aBSREL (Adaptive branch-site random effects likelihood)
                            uses an adaptive random effects branch-site model framework
                            to test whether each branch has evolved under positive selection,
                            using a procedure which infers an optimal number of rate categories per branch.",
                            terms.io.version : "2.2",
                            terms.io.reference : "Less Is More: An Adaptive Branch-Site Random Effects Model for Efficient Detection of Episodic Diversifying Selection (2015). Mol Biol Evol 32 (5): 1342-1353. v2.2 adds support for multiple-hit models",
                            terms.io.authors : "Sergei L Kosakovsky Pond, Ben Murrell, Steven Weaver and Temple iGEM / UCSD viral evolution group",
                            terms.io.contact : "spond@temple.edu",
                            terms.io.requirements : "in-frame codon alignment and a phylogenetic tree"
                          };


io.DisplayAnalysisBanner ( absrel.analysis_description );

absrel.json    = {
                    terms.json.analysis: absrel.analysis_description,
                    terms.json.input: {},
                    terms.json.fits : {},
                    terms.json.timers : {},
                    terms.json.test_results : {}
                  };


selection.io.startTimer (absrel.json [terms.json.timers], "Overall", 0);


/*------------------------------------------------------------------------------
    Key word arguments
*/

KeywordArgument ("code", "Which genetic code should be used", "Universal");
KeywordArgument ("alignment", "An in-frame codon alignment in one of the formats supported by HyPhy");
KeywordArgument ("tree", "A phylogenetic tree (optionally annotated with {})", null, "Please select a tree file for the data:");
KeywordArgument ("branches",  "Branches to test", "All");
// One additional KeywordArgument ("output") is called below after namespace absrel.

/*------------------------------------------------------------------------------
    Continued Analysis Setup
*/

namespace absrel {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ("absrel");
}

KeywordArgument ("multiple-hits",  "Include support for multiple nucleotide substitutions", "None");

absrel.multi_hit = io.SelectAnOption ({
                                        {"Double", "Include branch-specific rates for double nucleotide substitutions"}
                                        {"Double+Triple", "Include branch-specific rates for double and triple nucleotide substitutions"}
                                        {"None", "[Default] Use standard models which permit only single nucleotide changes to occur instantly"}
                                  }, "Include support for multiple nucleotide substitutions");


KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'ABSREL.json')", absrel.codon_data_info [terms.json.json]);

absrel.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");

io.CheckAssertion("utility.Array1D (absrel.partitions_and_trees) == 1", "aBSREL only works on a single partition dataset");

utility.ForEachPair (absrel.selected_branches, "_partition_", "_selection_",
    "_selection_ = utility.Filter (_selection_, '_value_', '_value_ == terms.tree_attributes.test');
     io.ReportProgressMessageMD('RELAX',  'selector', '* Selected ' + Abs(_selection_) + ' branches for testing: \\\`' + Join (', ',utility.Keys(_selection_)) + '\\\`')");

selection.io.startTimer (absrel.json [terms.json.timers], "Preliminary model fitting", 1);


namespace absrel {
    doGTR ("absrel");
}

selection.io.stopTimer (absrel.json [terms.json.timers], "Preliminary model fitting");
selection.io.startTimer (absrel.json [terms.json.timers], "Baseline model fitting", 2);

estimators.fixSubsetOfEstimates(absrel.gtr_results, absrel.gtr_results[terms.global]);

io.ReportProgressMessageMD ("absrel", "base", "Fitting the baseline model with a single dN/dS class per branch, and no site-to-site variation. ");


// estimators.FitCodonModel (codon_data, tree, "models.codon.MG_REV.ModelDescription", genetic_code, option, initial_values);

absrel.model_generator = "models.codon.MG_REV.ModelDescription";

if (absrel.multi_hit == "Double") {
    absrel.model_generator = "models.codon.MG_REV_MH.ModelDescription";
} else {
    if (absrel.multi_hit == "Double+Triple") {

        lfunction absrel.codon.MG_REV_TRIP._GenerateRate(fromChar, toChar, namespace, model_type, model) {
            return models.codon.MG_REV_TRIP._GenerateRate_generic (fromChar, toChar, namespace, model_type,
            model[utility.getGlobalValue("terms.translation_table")],
                "alpha", utility.getGlobalValue("terms.parameters.synonymous_rate"),
                "beta", utility.getGlobalValue("terms.parameters.nonsynonymous_rate"),
                "omega", utility.getGlobalValue("terms.parameters.omega_ratio"),
                "delta", utility.getGlobalValue("terms.parameters.multiple_hit_rate"),
                "psi", utility.getGlobalValue("terms.parameters.triple_hit_rate"),
                "psi", utility.getGlobalValue("terms.parameters.triple_hit_rate")
                );
        }


        lfunction absrel.model.MG_REV_no_islands (type, code) {
            mdl = models.codon.MG_REV_TRIP.ModelDescription (type,code);
            mdl[utility.getGlobalValue("terms.model.q_ij")] = "absrel.codon.MG_REV_TRIP._GenerateRate";
            return mdl;
        }
        absrel.model_generator = "absrel.model.MG_REV_no_islands";
    }
}

absrel.base.results = estimators.FitCodonModel (absrel.filter_names, absrel.trees, absrel.model_generator, absrel.codon_data_info [terms.code], {
    terms.run_options.model_type: terms.local,
    terms.run_options.retain_lf_object: TRUE,
    terms.run_options.retain_model_object : TRUE
}, absrel.gtr_results);


io.ReportProgressMessageMD("absrel", "base", "* " + selection.io.report_fit (absrel.base.results, 0, absrel.codon_data_info[terms.data.sample_size]));


absrel.baseline.branch_lengths = selection.io.extract_branch_info((absrel.base.results[terms.branch_length])[0], "selection.io.branch.length");
absrel.baseline.omegas = selection.io.extract_branch_info((absrel.base.results[terms.branch_length])[0], "absrel.local.omega");

if (absrel.multi_hit == "Double" || absrel.multi_hit == "Double+Triple") {
    absrel.baseline.deltas = selection.io.extract_branch_info((absrel.base.results[terms.branch_length])[0], "absrel.local.delta");
}


absrel.omega_stats = math.GatherDescriptiveStats (utility.Map (utility.UniqueValues (absrel.baseline.omegas), "_value_", "0+_value_"));

io.ReportProgressMessageMD("absrel", "base", "* Branch-level `terms.parameters.omega_ratio` distribution has median " +
                                             Format (absrel.omega_stats[terms.math.median], 5,2) + ", and 95% of the weight in " + Format (absrel.omega_stats[terms.math._2.5], 5,2) + " - " + Format (absrel.omega_stats[terms.math._97.5], 5,2));

if (absrel.multi_hit == "Double" || absrel.multi_hit == "Double+Triple") {
    absrel.baseline.deltas = selection.io.extract_branch_info((absrel.base.results[terms.branch_length])[0], "absrel.local.delta");
    absrel.delta_stats = math.GatherDescriptiveStats (utility.Map (utility.UniqueValues (absrel.baseline.deltas), "_value_", "0+_value_"));
    io.ReportProgressMessageMD("absrel", "base", "* Branch-level `terms.parameters.multiple_hit_rate` distribution has median " +
                                             Format (absrel.delta_stats[terms.math.median], 5,2) + ", and 95% of the weight in " + Format (absrel.delta_stats[terms.math._2.5], 5,2) + " - " + Format (absrel.delta_stats[terms.math._97.5], 5,2));
}

if (absrel.multi_hit == "Double+Triple") {
    absrel.baseline.psi = selection.io.extract_branch_info((absrel.base.results[terms.branch_length])[0], "absrel.local.psi");
    absrel.psi_stats = math.GatherDescriptiveStats (utility.Map (utility.UniqueValues (absrel.baseline.psi), "_value_", "0+_value_"));
    io.ReportProgressMessageMD("absrel", "base", "* Branch-level `terms.parameters.triple_hit_rate` distribution has median " +
                                             Format (absrel.psi_stats[terms.math.median], 5,2) + ", and 95% of the weight in " + Format (absrel.psi_stats[terms.math._2.5], 5,2) + " - " + Format (absrel.psi_stats[terms.math._97.5], 5,2));
}



selection.io.stopTimer (absrel.json [terms.json.timers], "Baseline model fitting");

// TODO -- there's gotta be a better way to do this
absrel.branch_count = Abs (absrel.baseline.branch_lengths);
absrel.sorted_branch_lengths = {absrel.branch_count, 2};
absrel.bnames = utility.Keys (absrel.baseline.branch_lengths);
utility.ForEachPair (absrel.bnames, "_index_", "_value_",
        '
            absrel.sorted_branch_lengths [_index_[1]][0] = absrel.baseline.branch_lengths[_value_];
            absrel.sorted_branch_lengths [_index_[1]][1] = _index_[1];
        ');
absrel.sorted_branch_lengths  = absrel.sorted_branch_lengths % 0;
absrel.names_sorted_by_length = {absrel.branch_count, 1};

for (absrel.i = absrel.branch_count - 1; absrel.i >= 0;  absrel.i = absrel.i - 1) {
    absrel.names_sorted_by_length [absrel.branch_count - 1 - absrel.i] =  absrel.bnames [absrel.sorted_branch_lengths[absrel.i][1]];
}

absrel.distribution_for_json = {absrel.per_branch_omega :
                                    {terms.math.mean : absrel.omega_stats[terms.math.mean],
                                    terms.math.median : absrel.omega_stats[terms.math.median],
                                    terms.math._2.5 : absrel.omega_stats[terms.math._2.5],
                                    terms.math._97.5 : absrel.omega_stats[terms.math._97.5]}
                               };


if (absrel.multi_hit == "Double" || absrel.multi_hit == "Double+Triple") {
    absrel.distribution_for_json [absrel.per_branch_delta] =
                                   {
                                        terms.math.mean : absrel.delta_stats[terms.math.mean],
                                        terms.math.median : absrel.delta_stats[terms.math.median],
                                        terms.math._2.5 : absrel.delta_stats[terms.math._2.5],
                                        terms.math._97.5 : absrel.delta_stats[terms.math._97.5]
                                   };

    if (absrel.multi_hit == "Double+Triple") {
            absrel.distribution_for_json [absrel.per_branch_psi] =
                                   {
                                        terms.math.mean : absrel.psi_stats[terms.math.mean],
                                        terms.math.median : absrel.psi_stats[terms.math.median],
                                        terms.math._2.5 : absrel.psi_stats[terms.math._2.5],
                                        terms.math._97.5 : absrel.psi_stats[terms.math._97.5]
                                   };

    }
}

//Store MG94 to JSON
selection.io.json_store_lf_withEFV (absrel.json,
                                     absrel.baseline_mg94xrev,
                                     absrel.base.results[terms.fit.log_likelihood],
                                     absrel.base.results[terms.parameters] ,
                                     absrel.codon_data_info[terms.data.sample_size],
                                     absrel.distribution_for_json,
                                     (absrel.base.results[terms.efv_estimate])["VALUEINDEXORDER"][0],
                                     absrel.display_orders[absrel.baseline_mg94xrev]);


selection.io.json_store_branch_attribute(absrel.json, absrel.baseline_mg94xrev, terms.branch_length, absrel.display_orders[absrel.baseline_mg94xrev],
                                                      0,
                                                      absrel.baseline.branch_lengths);

selection.io.json_store_branch_attribute(absrel.json, absrel.baseline_omega_ratio, terms.json.branch_label, absrel.display_orders[absrel.baseline_omega_ratio],
                                                      0,
                                                      absrel.baseline.omegas);



// define BS-REL models with up to N rate classes

absrel.model_defintions = {};

absrel.likelihood_function_id = absrel.base.results [terms.likelihood_function];
absrel.constrain_everything (absrel.likelihood_function_id);
absrel.tree_id = absrel.get_tree_name (absrel.likelihood_function_id);
absrel.model_id = absrel.get_model_id (absrel.likelihood_function_id);
absrel.MG94.model = (absrel.base.results[terms.model])[(utility.Keys (absrel.base.results[terms.model]))[0]];

absrel.temp = model.GetParameters_RegExp (absrel.MG94.model, terms.nucleotideRate("[ACGT]","[ACGT]"));
absrel.temp - terms.nucleotideRate("A","G");
absrel.full_model_parameters = {};
utility.AddToSet (absrel.full_model_parameters, absrel.temp);

selection.io.startTimer (absrel.json [terms.json.timers], "Complexity analysis", 3);

absrel.model_object_map = {
    absrel.MG94.model [terms.id] : absrel.MG94.model
};

absrel.model_defintions [1] = absrel.MG94.model;
for (absrel.i = 2; absrel.i <= absrel.max_rate_classes; absrel.i += 1) {
    absrel.model_defintions [absrel.i] = model.generic.DefineMixtureModel("absrel.BS_REL.ModelDescription",
            "absrel.model." + absrel.i, {
                "0": parameters.Quote(terms.local),
                "1": absrel.codon_data_info[terms.code],
                "2": parameters.Quote (absrel.i) // the number of rate classes
            },
            absrel.filter_names,
            None);

    absrel.model_object_map [(absrel.model_defintions [absrel.i])[terms.id]] = absrel.model_defintions [absrel.i];
    models.BindGlobalParameters ({"1" : absrel.model_defintions [absrel.i], "0" : absrel.MG94.model}, terms.nucleotideRate("[ACGT]","[ACGT]"));
}

io.ReportProgressMessageMD ("absrel", "complexity", "Determining the optimal number of rate classes per branch using a step up procedure");


absrel.current_parameter_count    = absrel.base.results[terms.parameters];
absrel.current_best_score         = math.GetIC (absrel.base.results[terms.fit.log_likelihood], absrel.current_parameter_count, absrel.codon_data_info[terms.data.sample_size]);
absrel.branch.complexity       = {};

utility.ToggleEnvVariable ("USE_LAST_RESULTS", TRUE);


if (absrel.multi_hit == "Double") {
    absrel.complexity_table.settings = {terms.table_options.header : TRUE, terms.table_options.column_widths: {
                "0": 35,
                "1": 10,
                "2": 10,
                "3": 20,
                "4": 12,
                "5": 15,
                "6": 15,
                "7": 15
                },
                terms.number_precision : 2};

    fprintf (stdout, "\n", io.FormatTableRow ({{"Branch", "Length", "Rates", "Max. dN/dS", "2x rate", "Log(L)", "AIC-c", "Best AIC-c so far"}}, absrel.complexity_table.settings));
} else {
    if (absrel.multi_hit == "Double+Triple") {
        absrel.complexity_table.settings = {terms.table_options.header : TRUE, terms.table_options.column_widths: {
                    "0": 35,
                    "1": 10,
                    "2": 10,
                    "3": 20,
                    "4": 12,
                    "5": 12,
                    "6": 15,
                    "7": 15,
                    "8": 15
                    },
                    terms.number_precision : 2};

        fprintf (stdout, "\n", io.FormatTableRow ({{"Branch", "Length", "Rates", "Max. dN/dS", "2x rate", "3x rate", "Log(L)", "AIC-c", "Best AIC-c so far"}}, absrel.complexity_table.settings));
    } else {
        absrel.complexity_table.settings = {terms.table_options.header : TRUE, terms.table_options.column_widths: {
                    "0": 35,
                    "1": 10,
                    "2": 10,
                    "3": 20,
                    "4": 15,
                    "5": 15,
                    "6": 15
                    },
                    terms.number_precision : 2};

        fprintf (stdout, "\n", io.FormatTableRow ({{"Branch", "Length", "Rates", "Max. dN/dS", "Log(L)", "AIC-c", "Best AIC-c so far"}}, absrel.complexity_table.settings));
    }

}

absrel.complexity_table.settings [terms.table_options.header] = FALSE;

utility.SetEnvVariable ("LF_SMOOTHING_SCALER", 0);

for (absrel.branch_id = 0; absrel.branch_id < absrel.branch_count; absrel.branch_id += 1) {

    absrel.current_branch           = absrel.names_sorted_by_length[absrel.branch_id];
    absrel.current_branch_estimates = absrel.GetBranchEstimates (absrel.MG94.model, absrel.tree_id, absrel.current_branch);
    absrel.report.row = {};
    absrel.report.row [0] =  absrel.current_branch;
    absrel.report.row [1] =  absrel.baseline.branch_lengths[absrel.current_branch];

    absrel.current_rate_count = 2;

    while (TRUE) {
        absrel.report.row [2] =  Format(absrel.current_rate_count,0,0);
        model.ApplyToBranch ((absrel.model_defintions [absrel.current_rate_count])[terms.id], absrel.tree_id, absrel.current_branch);
        parameters.SetValues (absrel.current_branch_estimates);


        //absrel.initial_guess =
        absrel.SetBranchConstraints (absrel.model_defintions [absrel.current_rate_count], absrel.tree_id, absrel.current_branch);

        Optimize (absrel.stepup.mles, ^absrel.likelihood_function_id
           , {
               "OPTIMIZATION_METHOD" : "nedler-mead",
               "OPTIMIZATION_PRECISION": 1e-4,
               "OPTIMIZATION_START_GRID" : absrel.ComputeOnAGrid (absrel.PopulateInitialGrid (absrel.model_defintions [absrel.current_rate_count], absrel.tree_id, absrel.current_branch, absrel.current_branch_estimates), absrel.likelihood_function_id)
             }
        );

        absrel.current_test_score = math.GetIC (absrel.stepup.mles[1][0], absrel.current_parameter_count + 2, absrel.codon_data_info[terms.data.sample_size]);

        absrel.provisional_estimates = absrel.GetBranchEstimates(absrel.model_defintions [absrel.current_rate_count], absrel.tree_id, absrel.current_branch);
        absrel.dn_ds.distro = absrel.GetRateDistribution (absrel.provisional_estimates);
        if (absrel.dn_ds.distro[absrel.current_rate_count-1][0] < 1000) {
            absrel.report.row [3] = Format (absrel.dn_ds.distro[absrel.current_rate_count-1][0],5,2) + " (" + Format (absrel.dn_ds.distro[absrel.current_rate_count-1][1]*100,5,2) + "%)";
        } else {
             absrel.report.row [3] = ">1000 (" + Format (absrel.dn_ds.distro[absrel.current_rate_count-1][1]*100,5,2) + "%)";
        }

        absrel.offset = 0;
        if (absrel.multi_hit == "Double" || absrel.multi_hit == "Double+Triple") {
            absrel.report.row [4] = (absrel.provisional_estimates[utility.getGlobalValue ("terms.parameters.multiple_hit_rate")])[utility.getGlobalValue ("terms.fit.MLE")];
            absrel.offset = 1;
        }
        if (absrel.multi_hit == "Double+Triple") {
            absrel.report.row [5] = (absrel.provisional_estimates[utility.getGlobalValue ("terms.parameters.triple_hit_rate")])[utility.getGlobalValue ("terms.fit.MLE")];
            absrel.offset += 1;
        }
        absrel.report.row [4+absrel.offset] = absrel.stepup.mles[1][0];
        absrel.report.row [5+absrel.offset] = absrel.current_test_score;


        if (absrel.current_test_score < absrel.current_best_score) {
            absrel.current_branch_estimates = absrel.provisional_estimates;
            absrel.current_best_score = absrel.current_test_score;
            absrel.report.row [6+absrel.offset] = absrel.current_test_score;
            fprintf (stdout, io.FormatTableRow (absrel.report.row, absrel.complexity_table.settings));

            if (absrel.current_rate_count >= absrel.max_rate_classes) {
                break;
            }
            absrel.current_rate_count      += 1;
            absrel.current_parameter_count += 2;
        } else {
            absrel.report.row [6+absrel.offset] = absrel.current_best_score;
            fprintf (stdout, io.FormatTableRow (absrel.report.row, absrel.complexity_table.settings));
            break;
        }
    }

    if (absrel.current_test_score >= absrel.current_best_score) { // reset the model
        absrel.current_rate_count = absrel.current_rate_count - 1;
        if ( absrel.current_rate_count >= 2) {
            model.ApplyToBranch ((absrel.model_defintions [absrel.current_rate_count])[terms.id], absrel.tree_id, absrel.current_branch);
            absrel.SetBranchConstraints (absrel.model_defintions [absrel.current_rate_count], absrel.tree_id, absrel.current_branch);
        } else {
            model.ApplyToBranch (absrel.MG94.model[terms.id], absrel.tree_id, absrel.current_branch);
        }
        parameters.SetValues (absrel.current_branch_estimates);
    }

    utility.AddToSet (absrel.full_model_parameters,
                     utility.Map (absrel.current_branch_estimates, "_parameter_", "_parameter_[terms.id]"));


    absrel.branch.complexity [absrel.current_branch] = absrel.current_rate_count;
    absrel.constrain_everything (absrel.likelihood_function_id);

}


selection.io.json_store_branch_attribute(absrel.json, absrel.rate_classes, terms.json.branch_label, absrel.display_orders[absrel.rate_classes],
                                                      0,
                                                      absrel.branch.complexity);


io.ReportProgressMessageMD ("absrel", "complexity-summary", "Rate class analyses summary");
utility.ForEachPair (utility.BinByValue (absrel.branch.complexity), "_rates_", "_branches_",
    "io.ReportProgressMessageMD('absrel',  'complexity-summary', '*  ' + Abs(_branches_) + ' branches with **' + _rates_ + '** rate classes')");

selection.io.stopTimer (absrel.json [terms.json.timers], "Complexity analysis");


selection.io.startTimer (absrel.json [terms.json.timers], "Full adaptive model fitting", 4);
io.ReportProgressMessageMD ("absrel", "Full adaptive model", "Improving parameter estimates of the adaptive rate class model");
parameters.RemoveConstraint (utility.Keys (absrel.full_model_parameters));


absrel.full_model.fit = estimators.FitExistingLF (absrel.likelihood_function_id,absrel.model_object_map);

absrel.full_model.mle_set = estimators.TakeLFStateSnapshot (absrel.likelihood_function_id);

io.ReportProgressMessageMD("absrel", "Full adaptive model", "* " + selection.io.report_fit (absrel.full_model.fit, 9, absrel.codon_data_info[terms.data.sample_size]));

selection.io.stopTimer (absrel.json [terms.json.timers], "Full adaptive model fitting");

selection.io.json_store_branch_attribute(absrel.json, absrel.full_adaptive_model, terms.branch_length, absrel.display_orders[absrel.full_adaptive_model],
                                             0,
                                             selection.io.extract_branch_info((absrel.full_model.fit[terms.branch_length])[0], "selection.io.branch.length"));

absrel.branch.rate_distributions = selection.io.extract_branch_info((absrel.full_model.fit[terms.branch_length])[0], "absrel.GetRateDistribution");


selection.io.json_store_branch_attribute(absrel.json, terms.json.rate_distribution, terms.json.branch_label, absrel.display_orders[terms.json.rate_distribution],
                                                      0,
                                                      absrel.branch.rate_distributions);


absrel.full_model.fit [terms.likelihood_function] = absrel.likelihood_function_id;
absrel.full_model.fit [terms.model] = absrel.model_object_map;

absrel.ESEN_trees = estimators.FitMGREVExtractComponentBranchLengths (absrel.codon_data_info , absrel.full_model.fit );

absrel.stree_info  = trees.ExtractTreeInfo ((absrel.ESEN_trees [terms.fit.synonymous_trees])[0]);
absrel.nstree_info = trees.ExtractTreeInfo ((absrel.ESEN_trees [terms.fit.nonsynonymous_trees])[0]);

selection.io.json_store_branch_attribute(absrel.json,
                                        absrel.full_adaptive_model.ES, terms.json.branch_label,absrel.display_orders[absrel.full_adaptive_model.ES],
                                             0,
                                             absrel.stree_info[terms.branch_length]);

selection.io.json_store_branch_attribute(absrel.json, absrel.full_adaptive_model.EN, terms.json.branch_label, absrel.display_orders[absrel.full_adaptive_model.EN],
                                            0,
                                            absrel.nstree_info[terms.branch_length]);




selection.io.json_store_lf (absrel.json,
                            absrel.full_adaptive_model,
                            absrel.full_model.fit[terms.fit.log_likelihood],
                            absrel.full_model.fit[terms.parameters] + 9 ,
                            absrel.codon_data_info[terms.data.sample_size],
                            {},
                            absrel.display_orders[absrel.full_adaptive_model]);

/***
    Testing individual branches for selection
***/

selection.io.startTimer (absrel.json [terms.json.timers], "Testing for selection", 5);
io.ReportProgressMessageMD ("absrel", "testing", "Testing selected branches for selection");

if (absrel.multi_hit == "Double") {
        absrel.testing_table.settings = {terms.table_options.header : TRUE, terms.table_options.column_widths: {
                    "0": 35,
                    "1": 10,
                    "2": 20,
                    "3": 12,
                    "4": 20,
                    "5": 20,
                    },
                    terms.number_precision : 2};

        fprintf (stdout, "\n", io.FormatTableRow ({{"Branch", "Rates", "Max. dN/dS", "2x rate", "Test LRT", "Uncorrected p-value"}}, absrel.testing_table.settings));
} else {
    if (absrel.multi_hit == "Double+Triple") {
        absrel.testing_table.settings = {terms.table_options.header : TRUE, terms.table_options.column_widths: {
                    "0": 35,
                    "1": 10,
                    "2": 20,
                    "3": 12,
                    "4": 12,
                    "5": 20,
                    "6": 20,
                    },
                    terms.number_precision : 2};

        fprintf (stdout, "\n", io.FormatTableRow ({{"Branch", "Rates", "Max. dN/dS", "2x rate","3x rate", "Test LRT", "Uncorrected p-value"}}, absrel.testing_table.settings));
    } else {
        absrel.testing_table.settings = {terms.table_options.header : TRUE, terms.table_options.column_widths: {
                    "0": 35,
                    "1": 10,
                    "2": 20,
                    "3": 20,
                    "4": 20,
                    },
                    terms.number_precision : 2};

        fprintf (stdout, "\n", io.FormatTableRow ({{"Branch", "Rates", "Max. dN/dS", "Test LRT", "Uncorrected p-value"}}, absrel.testing_table.settings));
    }

}



absrel.testing_table.settings [terms.table_options.header] = FALSE;
absrel.branch.p_values                                        = {};
absrel.branch.lrt                                             = {};
absrel.branch.delta                                           = {};
absrel.branch.psi                                             = {};

for (absrel.branch_id = 0; absrel.branch_id < absrel.branch_count; absrel.branch_id += 1) {
    absrel.current_branch           = absrel.names_sorted_by_length[absrel.branch_id];

    absrel.report.row = {};

    absrel.report.row [0] =  absrel.current_branch;
    absrel.report.row [1] =  Format (absrel.branch.complexity[absrel.current_branch], 3, 0);
    absrel.dn_ds.distro   = absrel.branch.rate_distributions [absrel.current_branch];

    if (absrel.dn_ds.distro[absrel.branch.complexity[absrel.current_branch]-1][0] < 1000) {
        absrel.report.row [2] = Format (absrel.dn_ds.distro[absrel.branch.complexity[absrel.current_branch]-1][0],5,2) + " (" + Format (absrel.dn_ds.distro[absrel.branch.complexity[absrel.current_branch]-1][1]*100,5,2) + "%)";
    } else {
        absrel.report.row [2] = ">1000 (" + Format (absrel.dn_ds.distro[absrel.branch.complexity[absrel.current_branch]-1][1]*100,5,2) + "%)";
    }

    absrel.offset = 0;
    absrel.final_estimates = absrel.GetBranchEstimates(absrel.model_defintions [absrel.branch.complexity[absrel.current_branch]], absrel.tree_id, absrel.current_branch);

    if (absrel.multi_hit == "Double" || absrel.multi_hit == "Double+Triple") {
        absrel.report.row [3] = (absrel.final_estimates[utility.getGlobalValue ("terms.parameters.multiple_hit_rate")])[utility.getGlobalValue ("terms.fit.MLE")];
        absrel.branch.delta [absrel.current_branch] =  absrel.report.row [3];
        absrel.offset = 1;
    }
    if (absrel.multi_hit == "Double+Triple") {
        absrel.report.row [4] = (absrel.final_estimates[utility.getGlobalValue ("terms.parameters.triple_hit_rate")])[utility.getGlobalValue ("terms.fit.MLE")];
        absrel.branch.psi [absrel.current_branch] =  absrel.report.row [4];
        absrel.offset += 1;
    }


    if ((absrel.selected_branches[0])[absrel.current_branch] == terms.tree_attributes.test) {
        if (absrel.dn_ds.distro [absrel.branch.complexity[absrel.current_branch]-1][0] > 1) {
            absrel.branch.ConstrainForTesting (absrel.model_defintions [absrel.branch.complexity[absrel.current_branch]], absrel.tree_id, absrel.current_branch);
            Optimize (absrel.null.mles, ^absrel.likelihood_function_id
                //, {"OPTIMIZATION_METHOD" : "nedler-mead", "OPTIMIZATION_PRECISION": 1e-3}
            );
            absrel.branch.test = absrel.ComputeLRT ( absrel.full_model.fit[terms.fit.log_likelihood], absrel.null.mles[1][0]);
            estimators.RestoreLFStateFromSnapshot (absrel.likelihood_function_id, absrel.full_model.mle_set);
        } else {
            absrel.branch.test = {terms.LRT : 0, terms.p_value : 1};
        }
        absrel.branch.p_values [ absrel.current_branch ] = absrel.branch.test [terms.p_value];
        absrel.branch.lrt       [absrel.current_branch] = absrel.branch.test [terms.LRT];
        absrel.report.row [3+absrel.offset] = absrel.branch.test [terms.LRT];
        absrel.report.row [4+absrel.offset] = Format (absrel.branch.test [terms.p_value], 8, 5);
    } else {
        absrel.branch.lrt [absrel.current_branch] = None;
        absrel.branch.p_values [absrel.current_branch] = None;
        absrel.report.row [3+absrel.offset] = "Not selected";
        absrel.report.row [4+absrel.offset] = "for testing";
    }

    fprintf (stdout, io.FormatTableRow (absrel.report.row, absrel.testing_table.settings));
}

if (Abs (absrel.branch.delta)) {
    selection.io.json_store_branch_attribute(absrel.json, terms.parameters.multiple_hit_rate, terms.json.branch_label, absrel.display_orders[terms.parameters.multiple_hit_rate],
                                                      0,
                                                      absrel.branch.delta);

}

if (Abs (absrel.branch.psi)) {
    selection.io.json_store_branch_attribute(absrel.json, terms.parameters.triple_hit_rate, terms.json.branch_label, absrel.display_orders[terms.parameters.triple_hit_rate],
                                                      0,
                                                      absrel.branch.psi);

}

selection.io.json_store_branch_attribute(absrel.json, terms.LRT, terms.json.branch_label, absrel.display_orders[terms.LRT],
                                                      0,
                                                       absrel.branch.lrt);

selection.io.json_store_branch_attribute(absrel.json, terms.json.uncorrected_pvalue, terms.json.branch_label, absrel.display_orders[terms.json.uncorrected_pvalue],
                                                      0,
                                                       absrel.branch.p_values);

absrel.branch.p_values.corrected = math.HolmBonferroniCorrection (absrel.branch.p_values);

selection.io.json_store_branch_attribute (absrel.json, terms.json.corrected_pvalue, terms.json.branch_label,  absrel.display_orders[terms.json.corrected_pvalue],
                                                       0,
                                                       absrel.branch.p_values.corrected);

absrel.test.all      = utility.Filter (absrel.branch.p_values.corrected, "_value_", "None!=_value_");
absrel.test.positive = utility.Filter (absrel.test.all, "_value_", "_value_<=absrel.p_threshold");

selection.io.stopTimer (absrel.json [terms.json.timers], "Testing for selection");



console.log ("----\n### Adaptive branch site random effects likelihood test ");
console.log ( "Likelihood ratio test for episodic diversifying positive selection at Holm-Bonferroni corrected _p = " + Format (absrel.p_threshold, 8, 4) + "_ found **" + Abs(absrel.test.positive) + "** branches under selection among **"+ Abs (absrel.test.all) + "** tested.\n");
utility.ForEachPair (absrel.test.positive, "_name_", "_p_",
            '
            console.log ("* " + _name_ + ", p-value = " + Format (_p_, 8,5));
            ');

absrel.json [terms.json.test_results] = {
                                             terms.json.pvalue_threshold : absrel.p_threshold,
                                             terms.json.tested  : Abs (absrel.test.all),
                                             terms.json.positive : Abs (absrel.test.positive)
                                         };

/***
    Cleanup
***/

selection.io.stopTimer (absrel.json [terms.json.timers], "Overall");
utility.ToggleEnvVariable ("USE_LAST_RESULTS", None);
io.SpoolJSON (absrel.json, absrel.codon_data_info [terms.json.json]);


return absrel.json;

//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------

lfunction absrel.ComputeLRT (ha, h0) {
    lrt  = 2*(ha-h0);
    return {utility.getGlobalValue("terms.LRT") : lrt,
            utility.getGlobalValue("terms.p_value") : (1-0.4*CChi2 (lrt,1)-0.6* CChi2 (lrt,2))*.5};
}


lfunction absrel.GetBranchEstimates (model, tree_id, branch_id) {
    values = {};
    utility.ForEachPair ((model[utility.getGlobalValue ("terms.parameters")])[utility.getGlobalValue ("terms.local")],
                         "_description_",
                         "_id_",
                         "`&values`[_description_] = {
                            terms.fit.MLE : Eval (`&tree_id` + '.' + `&branch_id` + '.' + _id_),
                            terms.id : `&tree_id` + '.' + `&branch_id` + '.' + _id_
                         };");

    return values;
}

//------------------------------------------------------------------------------------------------------------------------

lfunction absrel.GetRateDistribution (local_parameters) {


    offset = (^"absrel.multi_hit" == "Double") + (^"absrel.multi_hit" == "Double+Triple")*2;
    result = None;
    component_count = (Abs (local_parameters)-offset)$2;
    if (component_count > 1) {
        rates = {"rates" : {component_count,1}, "weights" : {component_count,1}};
        for (k = 1; k < component_count; k+=1) {
            (rates["rates"])[k-1]   = (local_parameters[terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), k)])[utility.getGlobalValue ("terms.fit.MLE")];
            (rates["weights"])[k-1]   = (local_parameters[terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), k)])[utility.getGlobalValue ("terms.fit.MLE")];
        }
        (rates["rates"])[component_count-1]   = (local_parameters[terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component_count)])[utility.getGlobalValue ("terms.fit.MLE")];
        result = parameters.GetStickBreakingDistribution (rates) % 0;
    } else {
        result = {{parameters.NormalizeRatio (
                                            (local_parameters[utility.getGlobalValue('terms.parameters.nonsynonymous_rate')])[utility.getGlobalValue ("terms.fit.MLE")],
                                            (local_parameters[utility.getGlobalValue('terms.parameters.synonymous_rate')])[utility.getGlobalValue ("terms.fit.MLE")]
                                          ),
                1}};
    }
    return result;
}

//------------------------------------------------------------------------------------------------------------------------

lfunction absrel.SetBranchConstraints (model, tree_id, branch_id) {
    component_count = model[utility.getGlobalValue ("terms.model.components")];
    local_parameters = (model[utility.getGlobalValue ("terms.parameters")])[utility.getGlobalValue ("terms.local")];
    parameters.SetRange ("`tree_id`.`branch_id`.`local_parameters[utility.getGlobalValue ('terms.parameters.synonymous_rate')]`", {utility.getGlobalValue ("terms.lower_bound") : "0", utility.getGlobalValue ("terms.upper_bound") : "50"});
    for (k = 1; k < component_count; k+=1) {
        omega_k   = terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), k);
        parameters.SetRange ("`tree_id`.`branch_id`.`local_parameters[omega_k]`", utility.getGlobalValue ("terms.range01"));
    }
    omega_k   = terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), k);
    parameters.SetRange ("`tree_id`.`branch_id`.`local_parameters[omega_k]`", utility.getGlobalValue ("terms.range_positive"));
}

//------------------------------------------------------------------------------------------------------------------------

lfunction absrel.branch.ConstrainForTesting (model, tree_id, branch_id) {
    component_count = model[utility.getGlobalValue ("terms.model.components")];
    local_parameters = (model[utility.getGlobalValue ("terms.parameters")])[utility.getGlobalValue ("terms.local")];
     if (component_count > 1) {
        omega_k   = terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component_count);
        parameters.SetConstraint ("`tree_id`.`branch_id`.`local_parameters[omega_k]`", "1", '');
    } else {
        parameters.SetConstraint (
            "`tree_id`.`branch_id`.`local_parameters[utility.getGlobalValue ('terms.parameters.nonsynonymous_rate')]`",
            "`tree_id`.`branch_id`.`local_parameters[utility.getGlobalValue ('terms.parameters.synonymous_rate')]`", '');
    }

}
//------------------------------------------------------------------------------------------------------------------------

lfunction absrel.PopulateInitialGrid (model, tree_id, branch_id, current_estimates) {

    component_count = model[utility.getGlobalValue ("terms.model.components")];
    local_parameters = (model[utility.getGlobalValue ("terms.parameters")])[utility.getGlobalValue ("terms.local")];

    grid = {};


    if (component_count == 2) {
        omega1   = terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), 1);
        omega2   = terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), 2);
        mixture1 = terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), 1 );

        if (^"absrel.multi_hit" == "Double")  {
                doubleH = utility.getGlobalValue ("terms.parameters.multiple_hit_rate");
                grid ["`tree_id`.`branch_id`.`local_parameters[^'terms.parameters.synonymous_rate']`"] = {4,1}["(current_estimates[^'terms.parameters.synonymous_rate'])[^'terms.fit.MLE']*(1+(2-_MATRIX_ELEMENT_ROW_)*0.25)"];
                grid ["`tree_id`.`branch_id`.`local_parameters[omega1]`"]   = {4,1}["_MATRIX_ELEMENT_ROW_ * 0.2"];
                grid ["`tree_id`.`branch_id`.`local_parameters[omega2]`"]   = {5,1}["(1+(_MATRIX_ELEMENT_ROW_-3)^3)*(_MATRIX_ELEMENT_ROW_>=3)+(_MATRIX_ELEMENT_ROW_*0.25+0.25)*(_MATRIX_ELEMENT_ROW_<3)"];
                grid ["`tree_id`.`branch_id`.`local_parameters[mixture1]`"] = {{0.98}{0.90}{0.75}{0.5}};
                grid ["`tree_id`.`branch_id`.`local_parameters[doubleH]`"] =  {{0.01}{0.2}{0.5}{1.0}};

        } else {
            if (^"absrel.multi_hit" == "Double+Triple")  {
                doubleH = utility.getGlobalValue ("terms.parameters.multiple_hit_rate");
                tripleH = utility.getGlobalValue ("terms.parameters.triple_hit_rate");
                grid ["`tree_id`.`branch_id`.`local_parameters[^'terms.parameters.synonymous_rate']`"] = {4,1}["(current_estimates[^'terms.parameters.synonymous_rate'])[^'terms.fit.MLE']*(1+(2-_MATRIX_ELEMENT_ROW_)*0.25)"];
                grid ["`tree_id`.`branch_id`.`local_parameters[omega1]`"]   = {3,1}["_MATRIX_ELEMENT_ROW_ * 0.2"];
                grid ["`tree_id`.`branch_id`.`local_parameters[omega2]`"]   = {4,1}["(1+(_MATRIX_ELEMENT_ROW_-3)^3)*(_MATRIX_ELEMENT_ROW_>=3)+(_MATRIX_ELEMENT_ROW_*0.25+0.25)*(_MATRIX_ELEMENT_ROW_<3)"];
                grid ["`tree_id`.`branch_id`.`local_parameters[mixture1]`"] = {{0.98}{0.90}{0.5}};
                grid ["`tree_id`.`branch_id`.`local_parameters[doubleH]`"] =  {{0.01}{0.2}{1.0}};
                grid ["`tree_id`.`branch_id`.`local_parameters[tripleH]`"] =  {{0.1}{1.0}{5.0}};

            } else {
                grid ["`tree_id`.`branch_id`.`local_parameters[^'terms.parameters.synonymous_rate']`"] = {5,1}["(current_estimates[^'terms.parameters.synonymous_rate'])[^'terms.fit.MLE']*(1+(2-_MATRIX_ELEMENT_ROW_)*0.25)"];
                grid ["`tree_id`.`branch_id`.`local_parameters[omega1]`"]   = {5,1}["_MATRIX_ELEMENT_ROW_ * 0.2"];
                grid ["`tree_id`.`branch_id`.`local_parameters[omega2]`"]   = {7,1}["(1+(_MATRIX_ELEMENT_ROW_-3)^3)*(_MATRIX_ELEMENT_ROW_>=3)+(_MATRIX_ELEMENT_ROW_*0.25+0.25)*(_MATRIX_ELEMENT_ROW_<3)"];
                grid ["`tree_id`.`branch_id`.`local_parameters[mixture1]`"] = {{0.98}{0.95}{0.90}{0.75}{0.5}};
            }
        }
    } else {
        omega_prev = current_estimates [terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component_count - 1)];
        omega_last = terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component_count);
        mixture_last = terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"),  component_count - 1);
        if (omega_prev [utility.getGlobalValue('terms.fit.MLE')] > 1) {
            parameters.SetValue (omega_prev [utility.getGlobalValue('terms.id')], 0.8);
        }
        grid ["`tree_id`.`branch_id`.`local_parameters[omega_last]`"]   = {10,1}["(1+(_MATRIX_ELEMENT_ROW_-5)^3)*(_MATRIX_ELEMENT_ROW_>=5)+(_MATRIX_ELEMENT_ROW_*0.15+0.15)*(_MATRIX_ELEMENT_ROW_<5)"];
        grid ["`tree_id`.`branch_id`.`local_parameters[mixture_last]`"] = {{0.98}{0.95}{0.90}{0.75}{0.5}};
    }

    return grid;
}


//------------------------------------------------------------------------------------------------------------------------

lfunction absrel.ComputeOnAGrid (grid_definition, lfname) {

    parameter_names = utility.Keys (grid_definition);
    parameter_count = Abs (grid_definition);
    grid_dimensions = {};
    total_grid_points = 1;
    explicit.grid = {};

    utility.ForEachPair (grid_definition, "_key_", "_value_",
    '
        `&grid_dimensions`[_key_] = utility.Array1D (_value_);
        `&total_grid_points` = `&total_grid_points` * `&grid_dimensions`[_key_];
    ');


    for (grid_point = 0; grid_point < total_grid_points; grid_point += 1) {
        index = grid_point;

        current_state = grid_dimensions;

        for (p_id = 0; p_id < parameter_count; p_id += 1) {
            p_name = parameter_names[p_id];
            current_state[p_name] = (grid_definition[p_name])[index % grid_dimensions[p_name]];
            index = index $ grid_dimensions[p_name];
        }

        explicit.grid + current_state;

    }


    return explicit.grid;

}

function absrel.SetValues(set) {
    if (Type (set) == "AssociativeList") {
        utility.ForEachPair (set, "_key_", "_value_",
        '
            parameters.SetValue (_key_, _value_);
        ');
    }
}

//----------------------------------------------------
lfunction absrel.get_tree_name (lf_id) {
    GetString (info, ^lf_id, -1);
    return (info["Trees"])[0];
}

lfunction absrel.get_model_id (lf_id) {
    GetString (info, ^lf_id, -1);
    return (info["Models"])[0];
}

function absrel.constrain_everything (lf_id) {
    GetString (absrel.constrain_everything.info, ^lf_id, -1);

    utility.ForEach (absrel.constrain_everything.info [utility.getGlobalValue("terms.parameters.global_independent")], "_value_",
                     "parameters.SetConstraint (_value_, Eval (_value_), terms.global)");
    utility.ForEach (absrel.constrain_everything.info [utility.getGlobalValue("terms.parameters.local_independent")], "_value_",
                     "parameters.SetConstraint (_value_, Eval (_value_), '')");
}

lfunction absrel.local.omega(branch_info) {
    return parameters.NormalizeRatio ((branch_info[utility.getGlobalValue ("terms.parameters.nonsynonymous_rate")])[utility.getGlobalValue("terms.fit.MLE")],
                                      (branch_info[utility.getGlobalValue ("terms.parameters.synonymous_rate")])[utility.getGlobalValue("terms.fit.MLE")]);
}

lfunction absrel.local.delta (branch_info) {
    return (branch_info[utility.getGlobalValue ("terms.parameters.multiple_hit_rate")])[utility.getGlobalValue("terms.fit.MLE")];
}

lfunction absrel.local.psi (branch_info) {
    return (branch_info[utility.getGlobalValue ("terms.parameters.triple_hit_rate")])[utility.getGlobalValue("terms.fit.MLE")];
}

//------------------------------------------------------------------------------

lfunction absrel.BS_REL.ModelDescription (type, code, components) {
    model = models.codon.BS_REL.ModelDescription(type, code, components);
    if (^"absrel.multi_hit" == "Double") {
            model [utility.getGlobalValue("terms.model.defineQ")] = "absrel.codon.BS_REL._DefineQ.DH";

    } else {
        if (^"absrel.multi_hit" == "Double+Triple") {
            model [utility.getGlobalValue("terms.model.defineQ")] = "absrel.codon.BS_REL._DefineQ.TH";
        } else {
            model [utility.getGlobalValue("terms.model.defineQ")] = "absrel.BS_REL._DefineQ";

        }
    }
    return model;
}


//------------------------------------------------------------------------------


lfunction absrel.codon.BS_REL._DefineQ.TH (bs_rel, namespace) {


    rate_matrices = {};

    bs_rel [utility.getGlobalValue("terms.model.q_ij")] = &rate_generator;
    bs_rel [utility.getGlobalValue("terms.mixture.mixture_components")] = {};

    _aux = parameters.GenerateSequentialNames ("bsrel_mixture_aux", bs_rel[utility.getGlobalValue("terms.model.components")] - 1, "_");
    _wts = parameters.helper.stick_breaking (_aux, None);


    for (component = 1; component <= bs_rel[utility.getGlobalValue("terms.model.components")]; component += 1) {
       key = "component_" + component;
       ExecuteCommands ("
        function rate_generator (fromChar, toChar, namespace, model_type, model) {
               return models.codon.MG_REV_TRIP._GenerateRate_generic (fromChar, toChar, namespace, model_type, model[utility.getGlobalValue('terms.translation_table')],
                'alpha', utility.getGlobalValue('terms.parameters.synonymous_rate'),
                'omega`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component),
                'omega`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component),
                'delta', utility.getGlobalValue('terms.parameters.multiple_hit_rate'),
                'psi', utility.getGlobalValue('terms.parameters.triple_hit_rate'),
                'psi', utility.getGlobalValue('terms.parameters.triple_hit_rate')
               );
            }"
       );

       if ( component < bs_rel[utility.getGlobalValue("terms.model.components")]) {
            model.generic.AddLocal ( bs_rel, _aux[component-1], terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), component ));
            parameters.SetRange (_aux[component-1], utility.getGlobalValue("terms.range_almost_01"));
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



lfunction absrel.codon.BS_REL._DefineQ.DH (bs_rel, namespace) {

    rate_matrices = {};

    bs_rel [utility.getGlobalValue("terms.model.q_ij")] = &rate_generator;
    bs_rel [utility.getGlobalValue("terms.mixture.mixture_components")] = {};

    _aux = parameters.GenerateSequentialNames ("bsrel_mixture_aux", bs_rel[utility.getGlobalValue("terms.model.components")] - 1, "_");
    _wts = parameters.helper.stick_breaking (_aux, None);


    for (component = 1; component <= bs_rel[utility.getGlobalValue("terms.model.components")]; component += 1) {
       key = "component_" + component;
       ExecuteCommands ("
        function rate_generator (fromChar, toChar, namespace, model_type, model) {
               return  models.codon.MG_REV_MH._GenerateRate_generic (fromChar, toChar, namespace, model_type, model[utility.getGlobalValue('terms.translation_table')],
                'alpha', utility.getGlobalValue('terms.parameters.synonymous_rate'),
                'omega`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component),
                'omega`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component),
                'delta', utility.getGlobalValue('terms.parameters.multiple_hit_rate')
               );
            }"
       );

       if ( component < bs_rel[utility.getGlobalValue("terms.model.components")]) {
            model.generic.AddLocal ( bs_rel, _aux[component-1], terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), component ));
            parameters.SetRange (_aux[component-1], utility.getGlobalValue("terms.range_almost_01"));
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

lfunction absrel.BS_REL._GenerateRate (fromChar, toChar, namespace, model_type, _tt, alpha, alpha_term, beta, beta_term, omega, omega_term) {

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
            } else {
                aa_rate = omega + "*" + alpha;
                (p[model_type])[omega_term] = omega;
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

lfunction absrel.BS_REL._DefineQ(bs_rel, namespace) {
    rate_matrices = {};

    bs_rel [utility.getGlobalValue("terms.model.q_ij")] = &rate_generator;
    bs_rel [utility.getGlobalValue("terms.mixture.mixture_components")] = {};

    _aux = parameters.GenerateSequentialNames ("bsrel_mixture_aux", bs_rel[utility.getGlobalValue("terms.model.components")] - 1, "_");
    _wts = parameters.helper.stick_breaking (_aux, None);
    mixture = {};

    component_count = bs_rel[utility.getGlobalValue("terms.model.components")];

    for (component = 1; component <= component_count; component += 1) {
       key = "component_" + component;
       ExecuteCommands ("
        function rate_generator (fromChar, toChar, namespace, model_type, model) {
           return absrel.BS_REL._GenerateRate (fromChar, toChar, namespace, model_type, model[utility.getGlobalValue('terms.translation_table')],
                'alpha', utility.getGlobalValue('terms.parameters.synonymous_rate'),
                'beta_`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.nonsynonymous_rate'), component),
                'omega`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component));
            }"
       );

       if ( component < component_count) {
            model.generic.AddLocal ( bs_rel, _aux[component-1], terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), component ));
            parameters.SetRange (_aux[component-1], utility.getGlobalValue("terms.range_almost_01"));
       }

       models.codon.generic.DefineQMatrix(bs_rel, namespace);
       rate_matrices [key] = bs_rel[utility.getGlobalValue("terms.model.rate_matrix")];
       (bs_rel [^'terms.mixture.mixture_components'])[key] = _wts [component-1];
    }


    bs_rel[utility.getGlobalValue("terms.model.rate_matrix")] = rate_matrices;
    parameters.SetConstraint(((bs_rel[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.global")])[terms.nucleotideRate("A", "G")], "1", "");
    return bs_rel;
}
