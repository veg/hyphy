RequireVersion ("2.31");


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

// namespace 'alignments' for various alignments related functions
LoadFunctionLibrary("libv3/tasks/alignments.bf");

// namespace 'tree' for various tree related functions
LoadFunctionLibrary("libv3/tasks/trees.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");


LoadFunctionLibrary("BranchSiteTemplate");
/*------------------------------------------------------------------------------
    BranchSiteTemplate Defines

    BuildCodonFrequencies (obsF);
    PopulateModelMatrix (ModelMatrixName&, EFV, synrateP, globalP, nonsynRateP);

------------------------------------------------------------------------------*/

RELAX.settings = {"GTR" : 1,
                  "LocalMG" : 1,
                  "Estimate GTR" : 1};

RELAX.timers  = {6,1};


relax.taskTimerStart (0);

RELAX.json    = {"fits" : {},
                  "timers" : {},
                  "relaxation-test" : None
                  };

RELAX.test            = "RELAX.test";
RELAX.reference       = "RELAX.reference";
RELAX.unclassified    = "RELAX.unclassified";
RELAX.relaxation_parameter        = "RELAX.K";

term.RELAX.k          = "relaxation coefficient";

/*------------------------------------------------------------------------------*/


io.DisplayAnalysisBanner ({"info" : "RELAX (a random effects test of selection relaxation)
                            uses a random effects branch-site model framework
                            to test whether a set of 'Test' branches evolves under relaxed
                            selection relative to a set of 'Reference' branches (R), as measured
                            by the relaxation parameter (K).",
                           "version" : "1.0",
                           "reference" : "RELAX: Detecting Relaxed Selection in a Phylogenetic Framework (2015). Mol Biol Evol 32 (3): 820-832",
                           "authors" : "Sergei L Kosakovsky Pond, Ben Murrell, Steven Weaver and Temple iGEM / UCSD viral evolution group",
                           "contact" : "spond@temple.edu",
                           "requirements" : "in-frame codon alignment and a phylogenetic tree, with at least two groups of branches defined using the {} notation (one group can be defined as all unlabeled branches)"
                          } );

relax.codon_data_info     = alignments.PromptForGeneticCodeAndAlignment ("RELAX.codon_data", "RELAX.codon_filter");
relax.sample_size         = relax.codon_data_info["sites"] * relax.codon_data_info["sequences"];

relax.name_mapping = relax.codon_data_info[utility.getGlobalValue("terms.json.name_mapping")];
    /**
        will contain "mapped" -> "original" associations with sequence names; or null if no mapping was necessary
    */

if (None == name_mapping) { /** create a 1-1 mapping if nothing was done */
    relax.name_mapping = {};
    utility.ForEach (alignments.GetSequenceNames ("RELAX.codon_data"), "_value_", "`&name_mapping`[_value_] = _value_");
}

relax.codon_data_info["json"] = relax.codon_data_info["file"] + ".RELAX.json";
io.ReportProgressMessage ("RELAX", "Loaded an MSA with " + relax.codon_data_info["sequences"] + " sequences and " + relax.codon_data_info["sites"] + " codons from '" + relax.codon_data_info["file"] + "'");

relax.codon_lists = models.codon.MapCode (relax.codon_data_info["code"]);
_Genetic_Code = relax.codon_data_info["code"];
    /*

     hack to make PopulateModelMatrix work

    */
relax.codon_frequencies     = frequencies._aux.CF3x4(frequencies._aux.empirical.collect_data ("RELAX.codon_filter",3,1,1),
models.DNA.alphabet, relax.codon_lists["sense"], relax.codon_lists["stop"]);


relax.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (relax.codon_data_info[utility.getGlobalValue("terms.json.partitions")], relax.name_mapping);


io.CheckAssertion("utility.Array1D (relax.partitions_and_trees) == 1", "RELAX only works on a single partition dataset");

relax.filter_specification = alignments.DefineFiltersForPartitions (relax.partitions_and_trees, "RELAX.codon_data" , "RELAX.codon_filter.", relax.codon_data_info);

relax.trees = utility.Map (relax.partitions_and_trees, "_partition_", '_partition_["tree"]');
relax.filter_names = utility.Map (relax.filter_specification, "_partition_", '_partition_["name"]');

relax.tree = relax.trees[0];

utility.SetEnvVariable ("VERBOSITY_LEVEL", 1);

relax.group_count = utility.Array1D(relax.tree["model_list"]);

assert (relax.group_count >= 2, "This analysis expects a partitioned tree with at least two groups defined (one can be the default group without a model)");

RELAX.groups = {};

utility.ForEach (relax.tree ["model_list"], "_model_", "
    if (_model_ == '') {
        RELAX.groups [RELAX.unclassified] = 'Unclassified branches';
    } else {
        RELAX.groups [_model_] =  _model_ + ' branches';
    }
");

RELAX.reference_set = io.SelectAnOption (RELAX.groups, "Choose the reference set of branches");
RELAX.reference_branches = utility.Filter (relax.tree["model_map"], "_value_", "_value_ == RELAX.reference_set");

io.ReportProgressMessage ("RELAX", "Selected " + Abs (RELAX.reference_branches) + " branches as the test set: " + Join (",", utility.Keys (RELAX.reference_branches)));

RELAX.has_unclassified = RELAX.groups / RELAX.unclassified;
RELAX.branch_to_partiton = {};
utility.ForEachPair (relax.tree ["model_map"], "_key_", "_value_", "if (Abs (_value_)) { RELAX.branch_to_partiton[_key_&&1] = _value_; } else { RELAX.branch_to_partiton[_key_&&1] = RELAX.unclassified;}");

RELAX.json ["partition"] = relax.selected_branches;
RELAX.json ["tree"] = relax.tree ["string"];

RELAX.runModel = io.SelectAnOption ({
                                        {"All", "[Default] Fit descriptive models AND run the relax test (4 models)"}
                                        {"Minimal", "Run only the RELAX test (2 models)"}
                                    }, "Analysis type");

relax.taskTimerStart (1);


if (RELAX.settings["GTR"]) {
    io.ReportProgressMessage ("RELAX", "Obtaining branch lengths under the GTR model");
    relax.gtr_results = estimators.FitGTR     ("RELAX.codon_filter", relax.tree, None);
    io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.gtr_results["LogL"]);
    estimators.fixSubsetOfEstimates (relax.gtr_results, relax.gtr_results["global"]);
} else {
    relax.gtr_results = None;
}

/*if (RELAX.settings["LocalMG"] && RELAX.runModel == "All") {
  io.ReportProgressMessage ("RELAX", "Obtaining  omega and branch length estimates under the local MG94xGTR model");
  relax.local_mg_results  = estimators.FitMGREV     (relax.filter_names, relax.trees, relax.codon_data_info ["code"], {"model-type" : terms.local}, relax.gtr_results);
  io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.local_mg_results["LogL"]);
  estimators.fixSubsetOfEstimates (relax.local_mg_results, relax.local_mg_results["global"]);
} else {
  relax.local_mg_results = relax.gtr_results;
}*/

parameters.DeclareGlobal ("relax.codon_branch_scaler", None);

utility.SetEnvVariable ("VERBOSITY_LEVEL", 1);

io.ReportProgressMessageMD ("RELAX", "mg-rev", "Obtaining omega and branch length estimates under the partitioned MG94xGTR model");
relax.mg_results  = estimators.FitMGREV     (relax.filter_names, relax.trees, relax.codon_data_info ["code"],
                                             {"model-type" : terms.local, "partitioned-omega" : {"0" : RELAX.branch_to_partiton}, "proportional-branch-length-scaler": {"0" : "relax.codon_branch_scaler"}},
                                             relax.local_mg_results);
relax.taskTimerStop (1);

io.ReportProgressMessageMD("RELAX", "mg-rev", "* Log(L) = " + Format(relax.mg_results["LogL"],8,2));
relax.global_dnds = selection.io.extract_global_MLE_re (relax.mg_results , "^" + terms.omega_ratio);

relax.mg_results_rate = {};

utility.ForEach (relax.global_dnds, "_value_",
    '
        io.ReportProgressMessageMD ("fel", "mg-rev", "* " + _value_["description"] + " = " + Format (_value_["MLE"],8,4));
    '
);


relax.json_store_lf (RELAX.json, "Partitioned MG94xREV",
                     relax.mg_results["LogL"], relax.mg_results["parameters"] + 5,
                     RELAX.timers[1],
                     relax._aux.extract_branch_info ((relax.mg_results["branch lengths"])[0], "relax.branch.length"),
                     relax._aux.extract_branch_info ((relax.mg_results["branch lengths"])[0], "relax.branch.omega"),
                     None,
                     None,
                     "&omega;"
                    );
relax.json_spool (RELAX.json, relax.codon_data_info["json"]);

relax.taskTimerStart (2);

RELAX.model_mapping                = {};
RELAX.model_assignment             = {};
RELAX.model_specification          = {};
RELAX.unclassified_model_id        = None;


for (g = 0; g < relax.group_count; g += 1) {
    relax.group_name = (relax.tree ["model_list"])[g];
    relax.branch_set = utility.Filter (relax.tree ["model_map"], "_value_", "_value_ == relax.group_name");

    if (Abs (  relax.group_name ) == 0) {
        relax.group_name = RELAX.unclassified;
    }


    relax.omega_estimate = Eval(utility.Values (selection.io.extract_global_MLE_re (relax.mg_results , relax.group_name + "\*$"))[0]);

    RELAX.model = relax.io.define_a_bsrel_model ("RELAX.`g`",
                                                  relax.codon_frequencies,
                                                  relax.omega_estimate["MLE"],
                                                  1);

    RELAX.model_assignment [RELAX.model["id"]]  =  RELAX.model;
    RELAX.model_specification [relax.group_name] = relax.branch_set;
    RELAX.model_mapping [relax.group_name] = RELAX.model["id"];

    if ((relax.tree ["model_list"])[g] == RELAX.reference_set) {
        RELAX.reference_model = RELAX.model["id"];
    }

    if (relax.group_name == RELAX.unclassified) {
        RELAX.unclassified_model_id = RELAX.model["id"];
    }
}

utility.ForEachPair ( RELAX.model_assignment, "_key_", "_value_", '
        if (_key_ != RELAX.reference_model) {
            parameters.ConstrainSets ((RELAX.model_assignment[RELAX.reference_model]) ["omegas"], (RELAX.model_assignment[_key_]) ["omegas"]);
            parameters.ConstrainSets ((RELAX.model_assignment[RELAX.reference_model]) ["f"], (RELAX.model_assignment[_key_]) ["f"]);
        }
    ');



model.ApplyModelToTree          ("RELAX.tree", relax.tree, RELAX.model_mapping, RELAX.model_specification);

ASSUME_REVERSIBLE_MODELS = 1;
LikelihoodFunction relax.LF = (RELAX.codon_filter, RELAX.tree);

global RELAX.branch_scaler = 4;
RELAX.proportional_constraint = "RELAX.branch_scaler";

if (RELAX.settings["Estimate GTR"] != TRUE) {
    estimators.fixSubsetOfEstimates   (relax.mg_results, relax.mg_results["global"]);
}

estimators.ApplyExistingEstimates ("relax.LF",  RELAX.model_assignment, relax.mg_results, None);

utility.SetEnvVariable ("USE_LAST_RESULTS", 1);

if (RELAX.runModel == "All") {

    io.ReportProgressMessageMD ("RELAX", "GDM", "## Two-stage fit of the general descriptive model (separate relaxation parameter for each branch)");

    OPTIMIZATION_PRECISION = 0.1;

    Optimize (relax.MLE.general_descriptive, relax.LF);

    OPTIMIZATION_PRECISION = 0.001;
    parameters.UnconstrainParameterSet ("relax.LF", {{terms.lf.local.constrained}});

    Optimize (relax.MLE.general_descriptive, relax.LF);
    io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.MLE.general_descriptive[1][0]);

    relax.general_descriptive = estimators.ExtractMLEs ("relax.LF", RELAX.model_specification);
    relax.add_scores (relax.general_descriptive, relax.MLE.general_descriptive);

    relax.taskTimerStop (2);

    relax.json_store_lf (RELAX.json, "General Descriptive",
                         relax.general_descriptive["LogL"], relax.general_descriptive["parameters"],
                         RELAX.timers[2],
                         relax._aux.extract_branch_info ((relax.general_descriptive["branch lengths"])[0], "relax.branch.length"),
                         relax._aux.extract_branch_info ((relax.general_descriptive["branch lengths"])[0], "relax.branch.local_k"),
                         {"All" : relax.getRateDistribution (RELAX.reference.model, 1)},
                         None,
                         "k"
                        );
    relax.json_spool (RELAX.json, relax.codon_data_info["json"]);

} else {
    parameters.UnconstrainParameterSet ("relax.LF", {{terms.lf.local.constrained}});
}

if (RELAX.has_unclassified) {
    parameters.RemoveConstraint ((RELAX.model_assignment[RELAX.unclassified_model_id]) ["omegas"]);
    parameters.RemoveConstraint ((RELAX.model_assignment[RELAX.unclassified_model_id]) ["f"]);
}

relax.taskTimerStart (3);
io.ReportProgressMessage ("RELAX", "Fitting the RELAX null model");

RELAX.null = relax.define.null ("RELAX.tree", RELAX.model_assignment[RELAX.reference_model], RELAX.model_specification);

Optimize (relax.MLE.null, relax.LF);
io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.MLE.null[1][0]);

LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (relax.codon_data_info["file"] + ".null.fit", CLEAR_FILE, relax.LF);
LIKELIHOOD_FUNCTION_OUTPUT = 2;

relax.null = estimators.ExtractMLEs ("relax.LF", RELAX.model_assignment);

relax.add_scores (relax.null, relax.MLE.null);

relax.taskTimerStop (3);

relax.omega_distributions = {};
relax.omega_distributions["Joint"] =  relax.getRateDistribution (RELAX.model_assignment[RELAX.reference_model], 1);
relax.omega_distributions["SRV"] = relax.getRateDistributionSRV (RELAX.model_assignment[RELAX.reference_model]);

if (RELAX.has_unclassified) {
    relax.omega_distributions ["Unclassified"] = relax.getRateDistribution (RELAX.model_assignment[RELAX.unclassified_model_id], 1);
}


//io.SpoolLF ("relax.LF", "/Volumes/home-raid/Desktop/null", None);

relax.json_store_lf (RELAX.json, "Null",
                     relax.null["LogL"], relax.null["parameters"],
                     RELAX.timers[3],
                     relax._aux.extract_branch_info ((relax.null["branch lengths"])[0], "relax.branch.length"),
                     relax._aux.extract_branch_info ((relax.null["branch lengths"])[0], "relax.branch.local_k"),
                     relax.omega_distributions,
                     1,
                     "k"
                    );

relax.json_spool (RELAX.json, relax.codon_data_info["json"]);

io.ReportProgressMessage ("RELAX", "Fitting the RELAX alternative model");

relax.taskTimerStart (4);

utility.ForEach (estimators.GetGlobalMLE_RegExp (relax.null, "^Relaxation parameter"), "_parameter_",
'
  parameters.RemoveConstraint (_parameter_["ID"]);
');


utility.SetEnvVariable ("VERBOSITY_LEVEL", 10);

Optimize (relax.MLE.alt, relax.LF);
relax.alt = estimators.ExtractMLEs ("relax.LF", RELAX.model_assignment);

io.ReportProgressMessageMD ("RELAX",  "AltModel", "Log(L) = " + relax.MLE.alt[1][0]);


utility.ForEachPair (estimators.GetGlobalMLE_RegExp (relax.alt, "^Relaxation parameter"), "_name_", "_parameter_",
'
  io.ReportProgressMessageMD ("RELAX",  "AltModel", "* " + _name_ + " = " + _parameter_ [terms.MLE]);

');


LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (relax.codon_data_info["file"] + ".alternative.fit", CLEAR_FILE, relax.LF);
LIKELIHOOD_FUNCTION_OUTPUT = 2;

relax.add_scores (relax.alt, relax.MLE.alt);

RELAX.json ["relaxation-test"] = relax.runLRT (relax.alt["LogL"], relax.null["LogL"]);

io.ReportProgressMessage ("RELAX", "AltModel", "Likelihood ratio test for relaxation on Test branches, p = " + (RELAX.json ["relaxation-test"])["p"]);

return 0;

relax.taskTimerStop  (4);

relax.omega_distributions["Test"] = relax.getRateDistribution (RELAX.null, Eval (RELAX.relaxation_parameter));
relax.omega_distributions["Reference"] =  relax.getRateDistribution (RELAX.null, 1);

if (RELAX.has_unclassified) {
    relax.omega_distributions ["Unclassified"] = relax.getRateDistribution (RELAX.unclassified.model, 1);
}

relax.json_store_lf (RELAX.json, "Alternative",
                     relax.alt["LogL"], relax.alt["parameters"],
                     RELAX.timers[4],
                     relax._aux.extract_branch_info ((relax.alt["branch lengths"])[0], "relax.branch.length"),
                     relax._aux.extract_branch_info ((relax.alt["branch lengths"])[0], "relax.branch.local_k"),
                     relax.omega_distributions,
                     Eval (RELAX.relaxation_parameter),
                     "k"
                    );


if (RELAX.runModel == "All") {

    relax.taskTimerStart  (5);

    parameters.RemoveConstraint (RELAX.test.model ["omegas"]);
    parameters.RemoveConstraint (RELAX.test.model ["f"]);
    parameters.SetConstraint    (RELAX.relaxation_parameter, Eval (RELAX.relaxation_parameter), "");


    io.ReportProgressMessage ("RELAX", "Fitting the RELAX partitioned exploratory model");
    Optimize (relax.MLE.part.expl, relax.LF);
    io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.MLE.part.expl [1][0]);

    LIKELIHOOD_FUNCTION_OUTPUT = 7;
    fprintf (relax.codon_data_info["file"] + ".partitioned_descriptive.fit", CLEAR_FILE, relax.LF);
    LIKELIHOOD_FUNCTION_OUTPUT = 2;

    relax.part.expl = estimators.ExtractMLEs ("relax.LF", RELAX.model_assignment);

    relax.omega_distributions["Test"]       = relax.getRateDistribution (RELAX.test.model,Eval (RELAX.relaxation_parameter));
    relax.omega_distributions["Reference"] =  relax.getRateDistribution (RELAX.reference.model, 1);

    if (RELAX.has_unclassified) {
        relax.omega_distributions ["Unclassified"] = relax.getRateDistribution (RELAX.unclassified.model, 1);
    }

    relax.add_scores (relax.part.expl, relax.MLE.part.expl);

    relax.taskTimerStop  (5);
    relax.json_store_lf (RELAX.json, "Partitioned Exploratory",
                         relax.part.expl["LogL"], relax.part.expl["parameters"],
                         RELAX.timers[5],
                         relax._aux.extract_branch_info ((relax.part.expl["branch lengths"])[0], "relax.branch.length"),
                         None,
                         relax.omega_distributions,
                         None,
                         ""
                        );
}


relax.taskTimerStop  (0);


(RELAX.json ["timers"])["Overall"]                  = RELAX.timers[0];
(RELAX.json ["timers"])["Preliminaries"]            = RELAX.timers[1];
(RELAX.json ["timers"])["General Descriptive"]      = RELAX.timers[2];
(RELAX.json ["timers"])["Null"]                     = RELAX.timers[3];
(RELAX.json ["timers"])["Alternative"]              = RELAX.timers[4];
(RELAX.json ["timers"])["Partitioned Descriptive"]  = RELAX.timers[5];

relax.json_spool (RELAX.json, relax.codon_data_info["json"]);

return RELAX.json;

//------------------------------------------------------------------------------
// HELPER FUNCTIONS FROM THIS POINT ON
//------------------------------------------------------------------------------

function relax.branch.length (branch_info) {
    return branch_info["MLE"];
}

function relax.branch.omega  (branch_info) {
    return parameters.NormalizeRatio ((branch_info[terms.nonsynonymous_rate])["MLE"], (branch_info[terms.synonymous_rate])["MLE"]);
}

function relax.branch.local_k  (branch_info) {
    return (branch_info[term.RELAX.k])["MLE"];
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
  relax.getRateInformation.rate_classes = Abs (model_description["omegas"]);
  relax.getRateInformation.omega_info = {relax.getRateInformation.rate_classes,2};

  for (relax.getRateInformation.k = 0; relax.getRateInformation.k < relax.getRateInformation.rate_classes; relax.getRateInformation.k += 1) {
    relax.getRateInformation.omega_info[relax.getRateInformation.k][0] = Eval ((model_description["omegas"])[relax.getRateInformation.k])^K;
    relax.getRateInformation.omega_info[relax.getRateInformation.k][1] = Eval ((model_description["weights"])[relax.getRateInformation.k]);
  }
  return relax.getRateInformation.omega_info % 0;
}

//------------------------------------------------------------------------------
function relax.getRateDistributionSRV (model_description) {
  relax.getRateInformation.rate_classes = Abs (model_description["srv_rates"]);
  relax.getRateInformation.omega_info = {relax.getRateInformation.rate_classes,2};

  for (relax.getRateInformation.k = 0; relax.getRateInformation.k < relax.getRateInformation.rate_classes; relax.getRateInformation.k += 1) {
    relax.getRateInformation.omega_info[relax.getRateInformation.k][0] = Eval ((model_description["srv_rates"])[relax.getRateInformation.k]);
    relax.getRateInformation.omega_info[relax.getRateInformation.k][1] = Eval ((model_description["srv_weights"])[relax.getRateInformation.k]);
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

    parameters.RemoveConstraint ((general_model["omegas"])[2]);

    relax.define.null.local = ((general_model["parameters"])["local"])[term.RELAX.k];

    relax.define.null.index = 0;

    utility.ForEachPair (partition, "_part_name_","_branch_list_", '
        relax.define.null.index += 1;
        if (_part_name_ == RELAX.reference_set || _part_name_ == RELAX.unclassified) {
            relax.define.null.global = "1";
            _branch_list_["relax.define.null._aux"][""];
        } else {
            relax.define.null.global = RELAX.relaxation_parameter + "_" + relax.define.null.index;
            (((RELAX.model_assignment[ RELAX.model_mapping[_part_name_]])["parameters"])["global"])["Relaxation parameter for *" + _part_name_ + "*"] = relax.define.null.global;
            parameters.SetConstraint (relax.define.null.global, 1, terms.global);
            _branch_list_["relax.define.null._aux"][""];
        }
    ');

    return general_model;
}



//------------------------------------------------------------------------------
function relax.add_scores (desc, mles) {
    if (Type (mles) == "Matrix") {
        desc ["LogL"] = mles[1][0];
        desc ["parameters"] = mles[1][1] + 9 + 5 * (RELAX.settings["Estimate GTR"] != 1);
    }
}

//------------------------------------------------------------------------------
function relax.runLRT (ha, h0) {
    return {"LR" : 2*(ha-h0),
            "p" : 1-CChi2 (2*(ha-h0),1)};
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
function relax._aux.define_srv_category (rates, weights, cat_name) {

    rate_count           = Abs (weights);
    expected_value  = "(" + weights[0] + ")*(" + rates[0] + ")";

    cat_freqs = "`cat_name`.weights";
    cat_rates = "`cat_name`.rates";
    cat_norm  = "`cat_name`.normalizer";
    ExecuteCommands ("`cat_freqs` = {rate_count,1}");
    ExecuteCommands ("`cat_rates` = {rate_count,1}");
    ExecuteCommands ("`cat_freqs`[0] := " + weights[0]);

    for (k = 1; k < rate_count; k+=1) {
        expected_value += "+(" + weights[k] + ")*(" + rates[k] + ")";
        ExecuteCommands ("`cat_freqs`[k] := " + weights[k]);
    }

    ExecuteCommands ("global `cat_norm` := `expected_value`");

    for (k = 0; k < rate_count; k+=1) {
        ExecuteCommands ("`cat_rates`[k] := (" + rates[k] + ")/`cat_norm`");
    }
    //category c  = ("+resp+", categFreqMatrix , MEAN, ,categRateMatrix, 0, 1e25);\n\n"
    ExecuteCommands ("category `cat_name`  = (rate_count, `cat_freqs`, MEAN, , `cat_rates`, 0,1e25)");
}

//------------------------------------------------------------------------------
function relax.io.define_a_bsrel_model (id, frequencies, mean_omega, do_local) {

    model_parameters = {"parameters" : {"global" : {}, "local" : {}}, "omegas" : {}, "weights" : {}, "f" : {}, "Q" : {}, "length" : ""};

    model_parameters["omegas"] = parameters.GenerateSequentialNames ("`id`.omega",    3, "_");
    model_parameters["f"]      = parameters.GenerateSequentialNames ("`id`.aux_freq", 2, "_");

    parameters.DeclareGlobal    (model_parameters["f"], None);
    parameters.SetRange         (model_parameters["f"], terms.range01);

    parameters.DeclareGlobal    (model_parameters["omegas"], None);

    model_parameters["weights"] = parameters.helper.stick_breaking (model_parameters["f"], {{0.7,0.25,0.05}});

    relax.init_omegas = {{0.05,0.25,4}};
    relax.init_omegas = relax.init_omegas * (1/ parameters.Mean (relax.init_omegas, model_parameters["weights"], Abs (model_parameters["omegas"])));

    parameters.SetRange ((model_parameters["omegas"])[0], terms.range_almost_01);
    parameters.SetValue ((model_parameters["omegas"])[0], relax.init_omegas[0]);

    parameters.SetRange ((model_parameters["omegas"])[1], terms.range_almost_01);
    parameters.SetValue ((model_parameters["omegas"])[1], relax.init_omegas[1]);

    parameters.SetRange ((model_parameters["omegas"])[2], terms.range_gte1);
    parameters.SetValue ((model_parameters["omegas"])[2], relax.init_omegas[2]);

    model_parameters["srv_rates"]    = parameters.GenerateSequentialNames ("relax.srv_rate", 3, "_");
    model_parameters["srv_f"]  = parameters.GenerateSequentialNames ("relax.aux_freq_srv", 2, "_");
    parameters.DeclareGlobal    (model_parameters["srv_f"], None);
    parameters.SetRange         (model_parameters["srv_f"], terms.range01);

    parameters.DeclareGlobal    (model_parameters["srv_rates"], None);
    model_parameters["srv_weights"] = parameters.helper.stick_breaking (model_parameters["srv_f"], {{0.7,0.25,0.05}});

    relax.init_srv = {{0.05,1,4}};

    parameters.SetRange ((model_parameters["srv_rates"])[0], terms.range01);
    parameters.SetValue ((model_parameters["srv_rates"])[0], relax.init_omegas[0]);

    parameters.SetConstraint ((model_parameters["srv_rates"])[1], 1,"");

    parameters.SetRange ((model_parameters["srv_rates"])[2], terms.range_gte1);
    parameters.SetValue ((model_parameters["srv_rates"])[2], relax.init_omegas[2]);

    relax._srv_cat_name = "relax.srv.cat";
    relax._aux.define_srv_category (model_parameters["srv_rates"],  model_parameters["srv_weights"], relax._srv_cat_name);


    if (do_local) {
        parameters.SetConstraint ((model_parameters["omegas"])[2], " 1/" + ((model_parameters["omegas"])[0]) + "/" + ((model_parameters["omegas"])[1]) , "");
        relax.io.define_a_bsrel_model_r = {"LB" : 1e-4, "UB" : 1};
        parameters.SetRange ((model_parameters["omegas"])[1], relax.io.define_a_bsrel_model_r);
    }


    local_k :< 50;

    relax.nuc = {4,3};
    for (k = 0; k < 4; k+=1) {
        for (k2 = 0; k2 < 3; k2 += 1) {
            relax.nuc [k][k2] = ((frequencies["bases"])[models.DNA.alphabet[k]])[k2];
        }
    }


    for (k = 1; k <= 3; k +=1) {
        ((model_parameters["parameters"])["global"])[relax.define_omega_term (k)] = (model_parameters["omegas"])[k-1];
        ((model_parameters["parameters"])["global"])[relax.define_srv_term (k)] = (model_parameters["srv_rates"])[k-1];
        if (k < 3) {
            ((model_parameters["parameters"])["global"])[relax.define_weight_term (k)] = (model_parameters["f"])[k-1];
            ((model_parameters["parameters"])["global"])[relax.define_srv_weight_term (k)] = (model_parameters["srv_f"])[k-1];
        }

        model_parameters["Q"] + ("Q_`id`_" + k);
        if (do_local) {
            PopulateModelMatrix			  ((model_parameters["Q"])[k-1],  relax.nuc, "t*`relax._srv_cat_name`", "Min (1000," + (model_parameters["omegas"])[k-1] +"^local_k)", "");
        } else {
            PopulateModelMatrix			  ((model_parameters["Q"])[k-1],  relax.nuc, "t*`relax._srv_cat_name`", (model_parameters["omegas"])[k-1], "");
        }
    }

    model_parameters["id"] = "`id`_model";
    model_parameters["length-expression"] = relax._aux.define_bsrel_model ("`id`_model", model_parameters["Q"], model_parameters["weights"], frequencies["codons"]);

    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("A","C")] = "AC";
    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("A","T")] = "AT";
    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("C","G")] = "CG";
    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("C","T")] = "CT";
    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("G","T")] = "GT";

    model_parameters["set-branch-length"] = "relax.aux.copy_branch_length";

    model_parameters["length parameter"] = "t";
    ((model_parameters["parameters"])[terms.local])[terms.timeParameter ()] = "t";
    if (do_local) {
        ((model_parameters["parameters"])[terms.local])[term.RELAX.k] = "local_k";
    }
    model_parameters["get-branch-length"] = "relax.aux.retrieve_branch_length";


    return model_parameters;
}

//------------------------------------------------------------------------------


function relax.aux.retrieve_branch_length (model, tree, node) {
    relax.aux.retrieve_branch_length.locals = Columns ((model_parameters["parameters"])[terms.local]);
    for (relax.aux.retrieve_branch_length.i = 0; relax.aux.retrieve_branch_length.i < Columns (relax.aux.retrieve_branch_length.locals); relax.aux.retrieve_branch_length.i += 1) {
        Eval (relax.aux.retrieve_branch_length.locals[relax.aux.retrieve_branch_length.i] + " = `tree`.`node`." + relax.aux.retrieve_branch_length.locals[relax.aux.retrieve_branch_length.i]);
    }
    return Eval (model["length-expression"]);
}

//------------------------------------------------------------------------------

function relax.aux.copy_branch_length (model, value, parameter) {

    relax.aux.copy_branch_length.t = ((model["parameters"])["local"])[terms.timeParameter ()];
    relax.aux.copy_branch_length.k = ((model["parameters"])["local"])[term.RELAX.k];

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
    (relax.tree ["model_map"])[key] = branch_set[value];
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

lfunction relax.define_srv_term (cat) {
    return "Synonymous rate variation for category " + cat;
}

lfunction relax.define_srv_weight_term (cat) {
    return "Synonymous frequency parameter " + cat;
}

//------------------------------------------------------------------------------------------------------------------------

function relax.json_spool (json, file) {
    USE_JSON_FOR_MATRIX = 1;
    fprintf (file, CLEAR_FILE, json);
    USE_JSON_FOR_MATRIX = 0;

}

//------------------------------------------------------------------------------------------------------------------------

function relax.json_store_lf (json, name, ll, df, time, branch_length, branch_annotation, omega_distribution, K, annotation_tag) {

    (json["fits"])[name] = {"log-likelihood"     : ll,
                            "parameters"         : df,
                            "AIC-c"              : relax.getIC (ll, df, relax.sample_size),
                            "runtime"            : time,
                            "branch-lengths"     : branch_length,
                            "branch-annotations" : branch_annotation,
                            "rate-distributions" : omega_distribution,
                            "K" : K,
                            "annotation-tag" : annotation_tag,
                            "display-order" : Abs (json["fits"])};

 }
