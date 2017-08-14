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

/*------------------------------------------------------------------------------*/

relax.json    = { terms.json.input: {},
                  terms.json.fits : {},
                  terms.json.timers : {},
                  terms.json.test_results : None
                  };

relax.test            = "RELAX.test";
relax.reference       = "RELAX.reference";
relax.unclassified    = "RELAX.unclassified";
relax.relaxation_parameter        = "RELAX.K";
relax.rate_classes     = 3;
relax.MG94 = "MG94xREV with separate rates for branch sets";

terms.relax.k          = "relaxation coefficient";

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
    load_file ({"prefix": "relax", "settings" : {"branch-selector" : "relax.select_branches"}});
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

ChoiceList  (relax.model_set,"RELAX analysis type",1,NO_SKIP,
            "All", "[Default] Fit descriptive models AND run the relax test (4 models)",
            "Minimal", "Run only the RELAX test (2 models)"
            );

io.CheckAssertion ("`&relax.model_set` >= 0", "User cancelled analysis selection");

selection.io.startTimer (relax.json [terms.json.timers], "Preliminary model fitting", 1);


namespace relax {
    doGTR ("relax");
}

estimators.fixSubsetOfEstimates(relax.gtr_results, relax.gtr_results["global"]);

namespace relax {
    scaler_prefix = "relax.scaler";
    doPartitionedMG ("relax", FALSE);
}


io.ReportProgressMessageMD ("RELAX", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");

relax.final_partitioned_mg_results = estimators.FitMGREV (relax.filter_names, relax.trees, relax.codon_data_info [terms.code], {
    "model-type": terms.local,
    "partitioned-omega": relax.selected_branches,
}, relax.partitioned_mg_results);


selection.io.stopTimer (relax.json [terms.json.timers], "Preliminary model fitting");
io.ReportProgressMessageMD("RELAX", "codon-refit", "* Log(L) = " + Format(relax.final_partitioned_mg_results[terms.fit.log_likelihood],8,2));



relax.global_dnds  = selection.io.extract_global_MLE_re (relax.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio);
relax.report_dnds  = {};



utility.ForEach (relax.global_dnds, "_value_", '
    io.ReportProgressMessageMD ("RELAX", "codon-refit", "* " + _value_["description"] + " = " + Format (_value_[terms.json.MLE],8,4));

    relax.report_dnds [(regexp.FindSubexpressions (_value_["description"], "^" + terms.parameters.omega_ratio + ".+\\*(.+)\\*$"))[1]] = {"0" : {terms.json.omega_ratio : _value_[terms.json.MLE], terms.json.proportion : 1}};
');

selection.io.json_store_lf (relax.json,
                            relax.MG94,
                            relax.final_partitioned_mg_results[terms.fit.log_likelihood],
                            relax.final_partitioned_mg_results[terms.parameters] + 9 , // +9 comes from CF3x4
                            math.GetIC (relax.final_partitioned_mg_results[terms.fit.log_likelihood], relax.final_partitioned_mg_results[terms.parameters] + 9, relax.codon_data_info[terms.data.sample_size]),
                            relax.report_dnds);


selection.io.stopTimer (relax.json [terms.json.timers], "Preliminary model fitting");

relax.test.bsrel_model =  model.generic.DefineMixtureModel("models.codon.BS_REL.ModelDescription",
        "relax.test", {
            "0": parameters.Quote(terms.local),
            "1": relax.codon_data_info["code"],
            "2": parameters.Quote (relax.rate_classes) // the number of rate classes
        },
        relax.filter_names,
        None);

relax.reference.bsrel_model =  model.generic.DefineMixtureModel("models.codon.BS_REL.ModelDescription",
        "relax.reference", {
            "0": parameters.Quote(terms.local),
            "1": relax.codon_data_info["code"],
            "2": parameters.Quote (relax.rate_classes) // the number of rate classes
        },
        relax.filter_names,
        None);


relax.global_distribution = {};

for (relax.i = 0; relax.i < relax.rate_classes; relax.i += 1) {
    relax.omega = "relax.global_omega_" + relax.i;
    parameters.DeclareGlobal (relax.omega, None);
    if (relax.i < relax.rate_classes - 1) {
        parameters.SetRange (relax.omega, terms.range01);
    } else {
        parameters.SetRange (relax.omega, terms.range_gte1);
    }
    relax.global_distribution [model.generic.GetLocalParameter (relax.test.bsrel_model, terms.AddCategory (utility.getGlobalValue('terms.parameters.nonsynonymous_rate'), relax.i + 1))] = relax.omega;
}


//models.BindGlobalParameters ({"0" : relax.test.bsrel_model, "1" : relax.reference.bsrel_model}, terms.mixture.mixture_aux_weight));

relax.model_object_map = { "relax.test" :       relax.test.bsrel_model };

if (relax.has_unclassified) {
    relax.unclassified.bsrel_model =  model.generic.DefineMixtureModel("models.codon.BS_REL.ModelDescription",
        "relax.unclassified", {
            "0": parameters.Quote(terms.global),
            "1": relax.codon_data_info["code"],
            "2": parameters.Quote (relax.rate_classes) // the number of rate classes
        },
        relax.filter_names,
        None);
}


if (relax.model_set == 0) { // run all the models
    selection.io.startTimer (relax.json [terms.json.timers], "General descriptive model fitting", 2);
    relax.general_exploratory =  estimators.FitLF (relax.filter_names,
                                        relax.trees,
                                        { "0" : {"DEFAULT" : "relax.test"}},
                                        relax.final_partitioned_mg_results,
                                        relax.model_object_map,
                                        {
                                            terms.run_options.retain_lf_object: TRUE,
                                            terms.run_options.apply_user_constraints: "relax.set_ge_constaints"
                                        });

    //
} else {

}



lfunction relax.set_ge_constaints(components, data_filter, tree, model_map, initial_values, model_objects) {
    console.log (components);
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

    ChoiceList(testSet, "Choose the set of branches to use the _test_ set", 1, NO_SKIP, selectTheseForTesting);
    io.CheckAssertion ("`&testSet` >= 0", "User cancelled branch selection; analysis terminating");
    if (option_count > 2) {
        ChoiceList(referenceSet, "Choose the set of branches to use the _reference_ set", 1, testSet, selectTheseForTesting);
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



if (RELAX.runModel == 0) {

    io.ReportProgressMessage ("RELAX", "Two-stage fit of the general descriptive model (separate relaxation parameter for each branch)");

    OPTIMIZATION_PRECISION = 0.1;

    Optimize (relax.MLE.general_descriptive, relax.LF);

    OPTIMIZATION_PRECISION = 0.001;
    parameters.UnconstrainParameterSet ("relax.LF", {{terms.parameters.local_constrained}});

    Optimize (relax.MLE.general_descriptive, relax.LF);
    io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.MLE.general_descriptive[1][0]);

    relax.general_descriptive = estimators.ExtractMLEs ("relax.LF", RELAX.model_specification);
    relax.add_scores (relax.general_descriptive, relax.MLE.general_descriptive);

    relax.taskTimerStop (2);

    relax.json_store_lf (relax.json, "General Descriptive",
                         relax.general_descriptive[terms.fit.log_likelihood], relax.general_descriptive[terms.parameters],
                         RELAX.timers[2],
                         relax._aux.extract_branch_info ((relax.general_descriptive[terms.branch_length])[0], "relax.branch.length"),
                         relax._aux.extract_branch_info ((relax.general_descriptive[terms.branch_length])[0], "relax.branch.local_k"),
                         {"All" : relax.getRateDistribution (RELAX.reference.model, 1)},
                         None,
                         "k"
                        );
    relax.json_spool (relax.json, relax.codon_data_info[terms.json.json]);

} else {
    parameters.UnconstrainParameterSet ("relax.LF", {{terms.parameters.local_constrained}});
}





if (RELAX.has_unclassified) {
    parameters.RemoveConstraint (RELAX.unclassified.model [terms.parameters.omegas]);
    parameters.RemoveConstraint (RELAX.unclassified.model [terms.parameters.freqs]);
}

relax.taskTimerStart (3);
io.ReportProgressMessage ("RELAX", "Fitting the RELAX null model");



RELAX.null = relax.define.null ("RELAX.tree", RELAX.reference.model, relax.selected_branches);




Optimize (relax.MLE.null, relax.LF);



io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.MLE.null[1][0]);


relax.null = estimators.ExtractMLEs ("relax.LF", RELAX.model_specification);


relax.add_scores (relax.null, relax.MLE.null);

relax.taskTimerStop (3);


relax.omega_distributions = {};
relax.omega_distributions[RELAX.test_branches_name] = relax.getRateDistribution (RELAX.test.model, 1);
relax.omega_distributions[RELAX.reference_branches_name] =  relax.getRateDistribution (RELAX.reference.model, 1);

if (RELAX.has_unclassified) {
    relax.omega_distributions [RELAX.unclassified_branches_name] = relax.getRateDistribution (RELAX.unclassified.model, 1);
}

//io.SpoolLF ("relax.LF", "/Volumes/home-raid/Desktop/null", None);

relax.json_store_lf (relax.json, "Null",
                     relax.null[terms.fit.log_likelihood], relax.null[terms.parameters],
                     RELAX.timers[3],
                     relax._aux.extract_branch_info ((relax.null[terms.branch_length])[0], "relax.branch.length"),
                     relax._aux.extract_branch_info ((relax.null[terms.branch_length])[0], "relax.branch.local_k"),
                     relax.omega_distributions,
                     1,
                     "k"
                    );

relax.json_spool (relax.json, relax.codon_data_info[terms.json.json]);

io.ReportProgressMessage ("RELAX", "Fitting the RELAX alternative model");

relax.taskTimerStart (4);
parameters.RemoveConstraint (RELAX.relaxation_parameter);
Optimize (relax.MLE.alt, relax.LF);
io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.MLE.alt[1][0] + ". Relaxation parameter K = " + Eval (RELAX.relaxation_parameter));


// Uncomment following lines to obtain full alternative fit information
//LIKELIHOOD_FUNCTION_OUTPUT = 7;
//fprintf (relax.codon_data_info["file"] + ".alternative.fit", CLEAR_FILE, relax.LF);
//LIKELIHOOD_FUNCTION_OUTPUT = 2;


relax.alt = estimators.ExtractMLEs ("relax.LF", RELAX.model_specification);
relax.add_scores (relax.alt, relax.MLE.alt);

relax.relaxation_test = relax.runLRT (relax.alt[terms.fit.log_likelihood], relax.null[terms.fit.log_likelihood]);

io.ReportProgressMessage ("RELAX", "Likelihood ratio test for relaxation on Test branches, p = " + relax.relaxation_test[terms.p_value]);

relax.taskTimerStop  (4);


relax.omega_distributions[RELAX.test_branches_name] = relax.getRateDistribution (RELAX.null, Eval (RELAX.relaxation_parameter));
relax.omega_distributions[RELAX.reference_branches_name] =  relax.getRateDistribution (RELAX.null, 1);

if (RELAX.has_unclassified) {
    relax.omega_distributions [RELAX.unclassified_branches_name] = relax.getRateDistribution (RELAX.unclassified.model, 1);
}

relax.json_store_lf (relax.json, "Alternative",
                     relax.alt[terms.fit.log_likelihood], relax.alt[terms.parameters],
                     RELAX.timers[4],
                     relax._aux.extract_branch_info ((relax.alt[terms.branch_length])[0], "relax.branch.length"),
                     relax._aux.extract_branch_info ((relax.alt[terms.branch_length])[0], "relax.branch.local_k"),
                     relax.omega_distributions,
                     Eval (RELAX.relaxation_parameter),
                     "k"
                    );


if (RELAX.runModel == 0) {

    relax.taskTimerStart  (5);

    parameters.RemoveConstraint (RELAX.test.model [terms.parameters.omegas]);
    parameters.RemoveConstraint (RELAX.test.model [terms.parameters.freqs]);
    parameters.SetConstraint    (RELAX.relaxation_parameter, Eval (RELAX.relaxation_parameter), "");


    io.ReportProgressMessage ("RELAX", "Fitting the RELAX partitioned exploratory model");
    Optimize (relax.MLE.part.expl, relax.LF);
    io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.MLE.part.expl [1][0]);

    // Uncomment following lines to obtain full partitioned descriptive fit information
    //LIKELIHOOD_FUNCTION_OUTPUT = 7;
    //fprintf (relax.codon_data_info["file"] + ".partitioned_descriptive.fit", CLEAR_FILE, relax.LF);
    //LIKELIHOOD_FUNCTION_OUTPUT = 2;

    relax.part.expl = estimators.ExtractMLEs ("relax.LF", RELAX.model_specification);

    relax.omega_distributions[RELAX.test_branches_name]       = relax.getRateDistribution (RELAX.test.model,Eval (RELAX.relaxation_parameter));
    relax.omega_distributions[RELAX.reference_branches_name] =  relax.getRateDistribution (RELAX.reference.model, 1);

    if (RELAX.has_unclassified) {
        relax.omega_distributions [RELAX.unclassified_branches_name] = relax.getRateDistribution (RELAX.unclassified.model, 1);
    }

    relax.add_scores (relax.part.expl, relax.MLE.part.expl);

    relax.taskTimerStop  (5);
    relax.json_store_lf (relax.json, "Partitioned Exploratory",
                         relax.part.expl[terms.fit.log_likelihood], relax.part.expl[terms.parameters],
                         RELAX.timers[5],
                         relax._aux.extract_branch_info ((relax.part.expl[terms.branch_length])[0], "relax.branch.length"),
                         None,
                         relax.omega_distributions,
                         None,
                         ""
                        );
}


relax.taskTimerStop  (0);


(relax.json [terms.json.timers])["Overall"]                  = RELAX.timers[0];
(relax.json [terms.json.timers])["Preliminaries"]            = RELAX.timers[1];
(relax.json [terms.json.timers])["General Descriptive"]      = RELAX.timers[2];
(relax.json [terms.json.timers])["Null"]                     = RELAX.timers[3];
(relax.json [terms.json.timers])["Alternative"]              = RELAX.timers[4];
(relax.json [terms.json.timers])["Partitioned Descriptive"]  = RELAX.timers[5];

//Test results
relax.json[terms.json.test_results] = relax.relaxation_test;

relax.json_spool (relax.json, relax.codon_data_info[terms.json.json]);

return relax.json;

//------------------------------------------------------------------------------
// HELPER FUNCTIONS FROM THIS POINT ON
//------------------------------------------------------------------------------


function relax.branch.length (branch_info) {
    return branch_info[terms.fit.MLE];
}

function relax.branch.omega  (branch_info) {
    return parameters.NormalizeRatio ((branch_info[terms.parameters.nonsynonymous_rate])[terms.fit.MLE], (branch_info[terms.parameters.synonymous_rate])[terms.fit.MLE]);
}

function relax.branch.local_k  (branch_info) {
    return (branch_info[term.RELAX.k])[terms.fit.MLE];
}

function relax._aux.extract_branch_info.callback (key, value) {
    relax._aux.extract_branch_info_result [key] = utility.CallFunction (callback, {"0" : "value"});
}

function relax._aux.extract_branch_info (branch_spec, callback) {
    relax._aux.extract_branch_info_result = {};
    branch_spec ["relax._aux.extract_branch_info.callback"][""];
    return relax._aux.extract_branch_info_result;
}


//------------------------------------------------------------------------------
function relax.getRateDistribution (model_description, K) {
  relax.getRateInformation.rate_classes = Abs (model_description[terms.parameters.omegas]);
  relax.getRateInformation.omega_info = {relax.getRateInformation.rate_classes,2};


  for (relax.getRateInformation.k = 0; relax.getRateInformation.k < relax.getRateInformation.rate_classes; relax.getRateInformation.k += 1) {
    relax.getRateInformation.omega_info[relax.getRateInformation.k][0] = Eval ((model_description[terms.parameters.omegas])[relax.getRateInformation.k])^K;
    relax.getRateInformation.omega_info[relax.getRateInformation.k][1] = Eval ((model_description[terms.parameters.weights])[relax.getRateInformation.k]);
  }
  return relax.getRateInformation.omega_info % 0;
}



//------------------------------------------------------------------------------
function relax.define.null._aux (key, value) {
    //fprintf (stdout, "`tree_id`.`key`.`relax.define.null.local` := `relax.define.null.global`\n");
    ExecuteCommands ("`tree_id`.`key`.`relax.define.null.local` := `relax.define.null.global`");
}

//------------------------------------------------------------------------------
function relax.define.null (tree_id, general_model, partition) {
    RELAX.null = general_model;
    parameters.RemoveConstraint ((general_model[terms.parameters.omegas])[2]);

    relax.define.null.local = ((general_model[terms.parameters])[terms.local])[term.RELAX.k];

    //fprintf (stdout, "\n**", relax.define.null.par , "\n");
    relax.define.null.global = "1";
    (partition[RELAX.reference])["relax.define.null._aux"][""];

    parameters.SetConstraint (RELAX.relaxation_parameter, 1, terms.global);
    relax.define.null.global = RELAX.relaxation_parameter;

    (partition[RELAX.test])["relax.define.null._aux"][""];

    ((general_model[terms.parameters])[terms.global])[RELAX.relaxation_parameter] = RELAX.relaxation_parameter;

    return RELAX.null;
}



//------------------------------------------------------------------------------
function relax.add_scores (desc, mles) {
    if (Type (mles) == "Matrix") {
        desc [terms.fit.log_likelihood] = mles[1][0];
        desc [terms.parameters] = mles[1][1] + 9 + 5 * (RELAX.settings["Estimate GTR"] != 1);
    }
}

//------------------------------------------------------------------------------
function relax.runLRT (ha, h0) {
    return {terms.LRT : 2*(ha-h0),
            terms.p_value : 1-CChi2 (2*(ha-h0),1)};
}


//------------------------------------------------------------------------------
function relax._aux.define_bsrel_model (id,Q,weights,freqs) {
    rate_count = Abs (Q);
    components = {};
    length = "";
    Eval ("`id`_eqf = freqs");

    for (k = 0; k < rate_count; k+=1) {
        components[k] = "Exp(" + Q[k] + ")*" + weights[k];
        ExecuteCommands ("Model busted._aux.define_bsrel_model_bl = (" + Q[k] + ",`id`_eqf,0)");
        GetString (blExp, busted._aux.define_bsrel_model_bl, -1);
        if (k) {
            length += "+";
        }
        length += "(`blExp`)*" + weights[k];
    }

    ExecuteCommands ("Model `id` =(\"" + Join("+",components) + "\",`id`_eqf,EXPLICIT_FORM_MATRIX_EXPONENTIAL);");
    return length;
}



//------------------------------------------------------------------------------

function relax.io.evaluator (key, value) {
    fprintf (stdout, key, "->", Eval (value), "\n");
}

//------------------------------------------------------------------------------
function relax.io.define_a_bsrel_model (id, frequencies, mean_omega, do_local) {

    model_parameters = {terms.parameters : {terms.global : {}, terms.local : {}}, terms.parameters.omegas : {}, terms.parameters.weights : {}, terms.parameters.freqs : {}, terms.model.rate_matrix : {}, terms.model.length : ""};

    model_parameters[terms.parameters.omegas] = parameters.GenerateSequentialNames ("`id`.omega",    3, "_");
    model_parameters[terms.parameters.freqs]      = parameters.GenerateSequentialNames ("`id`.aux_freq", 2, "_");

    parameters.DeclareGlobal    (model_parameters[terms.parameters.freqs], None);
    parameters.SetRange         (model_parameters[terms.parameters.freqs], terms.range01);

    parameters.DeclareGlobal    (model_parameters[terms.parameters.omegas], None);

    model_parameters[terms.parameters.weights] = parameters.helper.stick_breaking (model_parameters[terms.parameters.freqs], {{0.7,0.25,0.05}});

    relax.init_omegas = {{0.05,0.25,4}};
    relax.init_omegas = relax.init_omegas * (1/ parameters.Mean (relax.init_omegas, model_parameters[terms.parameters.weights], Abs (model_parameters[terms.parameters.omegas])));

    parameters.SetRange ((model_parameters[terms.parameters.omegas])[0], terms.range_almost_01);
    parameters.SetValue ((model_parameters[terms.parameters.omegas])[0], relax.init_omegas[0]);

    parameters.SetRange ((model_parameters[terms.parameters.omegas])[1], terms.range_almost_01);
    parameters.SetValue ((model_parameters[terms.parameters.omegas])[1], relax.init_omegas[1]);

    parameters.SetRange ((model_parameters[terms.parameters.omegas])[2], terms.range_gte1);
    parameters.SetValue ((model_parameters[terms.parameters.omegas])[2], relax.init_omegas[2]);

    if (do_local) {
        parameters.SetConstraint ((model_parameters[terms.parameters.omegas])[2], " 1/" + ((model_parameters[terms.parameters.omegas])[0]) + "/" + ((model_parameters[terms.parameters.omegas])[1]) , "");
        relax.io.define_a_bsrel_model_r = {terms.lower_bound : 1e-4, terms.upper_bound : 1};
        parameters.SetRange ((model_parameters[terms.parameters.omegas])[1], relax.io.define_a_bsrel_model_r);
    }


    local_k :< 50;

    relax.nuc = {4,3};
    for (k = 0; k < 4; k+=1) {
        for (k2 = 0; k2 < 3; k2 += 1) {
            relax.nuc [k][k2] = ((frequencies[terms.bases])[models.DNA.alphabet[k]])[k2];
        }
    }


    for (k = 1; k <= 3; k +=1) {
        ((model_parameters[terms.parameters])[terms.global])[relax.define_omega_term (k)] = (model_parameters[terms.parameters.omegas])[k-1];
        if (k < 3) {
            ((model_parameters[terms.parameters])[terms.global])[relax.define_weight_term (k)] = (model_parameters[terms.parameters.freqs])[k-1];
        }

        model_parameters[terms.model.rate_matrix] + ("Q_`id`_" + k);
        if (do_local) {
            PopulateModelMatrix			  ((model_parameters[terms.model.rate_matrix])[k-1],  relax.nuc, terms.parameters.default_time, "Min (1000," + (model_parameters[terms.parameters.omegas])[k-1] +"^local_k)", "");
        } else {
            PopulateModelMatrix			  ((model_parameters[terms.model.rate_matrix])[k-1],  relax.nuc, terms.parameters.default_time, (model_parameters[terms.parameters.omegas])[k-1], "");
        }
    }

    model_parameters[terms.id] = "`id`_model";
    model_parameters[terms.model.length_expression] = relax._aux.define_bsrel_model ("`id`_model", model_parameters[terms.model.rate_matrix], model_parameters[terms.parameters.weights], frequencies[terms.codons]);

    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("A","C")] = "AC";
    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("A","T")] = "AT";
    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("C","G")] = "CG";
    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("C","T")] = "CT";
    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("G","T")] = "GT";

    model_parameters[terms.model.set_branch_length] = "relax.aux.copy_branch_length";

    model_parameters[terms.model.length_parameter] = terms.parameters.default_time;
    ((model_parameters[terms.parameters])[terms.local])[terms.timeParameter ()] = terms.parameters.default_time;
    if (do_local) {
        ((model_parameters[terms.parameters])[terms.local])[term.RELAX.k] = "local_k";
    }
    model_parameters[terms.model.get_branch_length] = "relax.aux.retrieve_branch_length";
    return model_parameters;
}

//------------------------------------------------------------------------------


function relax.aux.retrieve_branch_length (model, tree, node) {
    relax.aux.retrieve_branch_length.locals = Columns ((model_parameters[terms.parameters])[terms.local]);
    for (relax.aux.retrieve_branch_length.i = 0; relax.aux.retrieve_branch_length.i < Columns (relax.aux.retrieve_branch_length.locals); relax.aux.retrieve_branch_length.i += 1) {
        Eval (relax.aux.retrieve_branch_length.locals[relax.aux.retrieve_branch_length.i] + " = `tree`.`node`." + relax.aux.retrieve_branch_length.locals[relax.aux.retrieve_branch_length.i]);
    }
    return Eval (model[terms.model.length_expression]);
}

//------------------------------------------------------------------------------

function relax.aux.copy_branch_length (model, value, parameter) {

    relax.aux.copy_branch_length.t = ((model[terms.parameters])[terms.local])[terms.timeParameter ()];
    relax.aux.copy_branch_length.k = ((model[terms.parameters])[terms.local])[term.RELAX.k];

    if (Abs (RELAX.proportional_constraint)) {
        Eval ("`parameter`.`relax.aux.copy_branch_length.t` := `RELAX.proportional_constraint` * " + value);
    } else {
        Eval ("`parameter`.`relax.aux.copy_branch_length.t` = " + value);
    }

    if (Type (relax.aux.copy_branch_length.k) == "String") {
        Eval ("`parameter`.`relax.aux.copy_branch_length.k` = 1");
    }
}

function relax._aux.io.countBranchSets (key, value) {
    available_models[value] += 1;
    return None;
}

function relax._aux.io.mapBranchSets (key, value) {
    /*if (Abs (value) == 0) {
        value = RELAX.unclassified;
    }*/
    (relax.tree [terms.trees.model_map])[key] = branch_set[value];
    (return_set[branch_set[value]])[key] = 1;
    return None;
}

function relax.handle.unlabeled (label) {
    if (label == "Unlabeled branches") {
        return "";
    } else {
        return label;
    }
}

//------------------------------------------------------------------------------
function relax.io.defineBranchSets (relax.tree) {

    available_models        = {};
    branch_set              = {};


    for (k = 0; k < Columns (relax.tree[terms.trees.model_list]); k += 1) {
        available_models  [(relax.tree[terms.trees.model_list])[k]] = 0;
    }

    (relax.tree[terms.trees.model_map])["relax._aux.io.countBranchSets"][""];


    list_models = Rows (available_models); // get keys
    io.CheckAssertion ("Columns (list_models) > 1", "The tree string must include at least one two sets of branches, at least one of which is annotated using {}");

    option_count = Columns (list_models);

    selectTheseForTesting = {option_count,2};
    for (k = 0; k < Columns (list_models); k+=1) {
        if (list_models[k] != "") {
            selectTheseForTesting [k][0] = list_models[k];
            selectTheseForTesting [k][1] = "Set " + list_models [k] + " with " + available_models[list_models[k]] + " branches";
        } else {
            selectTheseForTesting [k][0] = "Unlabeled branches";
            selectTheseForTesting [k][1] = "Set of " + available_models[list_models [k]] + " unlabeled branches";
        }
    }


    ChoiceList  (fgSet,"Choose the set of branches to test for relaxed selection (T set)",1,NO_SKIP,selectTheseForTesting);

    return_set = {};

    if (fgSet >= 0) {
        branch_set [relax.handle.unlabeled(selectTheseForTesting[fgSet][0])] = RELAX.test;
        return_set [RELAX.test] = {};
        if (option_count > 2) {
            ChoiceList  (bgSet,"Choose the set of reference branches (R set)",1,fgSet,selectTheseForTesting);
            if (bgSet < 0) {
                return {};
            }
            for (k = 0; k < option_count; k+=1) {
                if (k != bgSet && k != fgSet) {
                    branch_set [relax.handle.unlabeled(selectTheseForTesting[k][0])] = RELAX.unclassified;
                    return_set [RELAX.unclassified] = {};
                }
            }
        }
        else {
            bgSet = 1-fgSet;
        }
        branch_set [relax.handle.unlabeled(selectTheseForTesting[bgSet][0])] = RELAX.reference;
        return_set [RELAX.reference] = {};
     }

    (relax.tree[terms.trees.model_map])["relax._aux.io.mapBranchSets"][""];
    relax.tree[terms.trees.model_list] = Columns (relax.tree[terms.trees.model_map]);
    //fprintf (stdout, "\n", relax.tree, "\n", return_set, "\n");

    return return_set;

}


//------------------------------------------------------------------------------------------------------------------------

function relax.taskTimerStart (index) {
    RELAX.timers[index] = Time(1);
}

function relax.taskTimerStop (index) {
    RELAX.timers[index] = Time(1) - RELAX.timers[index];
}

//------------------------------------------------------------------------------------------------------------------------

function relax.getIC (logl,params,samples) {
    return -2*logl + 2*samples/(samples-params-1)*params;
}

//------------------------------------------------------------------------------------------------------------------------

lfunction relax.define_omega_term (cat) {
    return "Omega for category " + cat;
}

lfunction relax.define_weight_term (cat) {
    return "Omega frequency parameter " + cat;
}

//------------------------------------------------------------------------------------------------------------------------

function relax.json_spool (json, file) {
    USE_JSON_FOR_MATRIX = 1;
    fprintf (file, CLEAR_FILE, json);
    USE_JSON_FOR_MATRIX = 0;

}

//------------------------------------------------------------------------------------------------------------------------

// SHIM to the old global namespace for CountSenseCodons

lfunction CountSenseCodons (code) {
    return genetic_code.CountSense (code);
}



//------------------------------------------------------------------------------------------------------------------------

function relax.json_store_lf (json, name, ll, df, time, branch_length, branch_annotation, omega_distribution, K, annotation_tag) {

    (json[terms.json.fits])[name] = {terms.json.log_likelihood     : ll,
                            terms.json.parameters         : df,
                            terms.json.AICc              : relax.getIC (ll, df, relax.sample_size),
                            terms.json.runtime           : time,
                            terms.json.branch_lengths     : branch_length,
                            terms.json.branch_annotations : branch_annotation,
                            terms.json.rate_distributions : omega_distribution,
                            "K" : K,
                            terms.json.annotation_tag : annotation_tag,
                            terms.json.display_order : Abs (json[terms.json.fits])};

 }
