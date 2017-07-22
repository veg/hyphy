RequireVersion ("2.31");


//LF_SMOOTHING_SCALER         = 0.1;

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


LoadFunctionLibrary("libv3/models/terms.bf");
LoadFunctionLibrary("libv3/terms-json.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");


LoadFunctionLibrary("BranchSiteTemplate");
/*------------------------------------------------------------------------------
    BranchSiteTemplate Defines

    BuildCodonFrequencies (obsF);
    PopulateModelMatrix (ModelMatrixName&, EFV, synrateP, globalP, nonsynRateP);

------------------------------------------------------------------------------*/

// More sensible for users to see test results first. Can easily update this.
RELAX.display_orders = { "Null": 0,
                         "Alternative": 1,
                         "Partitioned MG94xREV": 2,
                         "General Descriptive": 3,
                         "Partitioned Descriptive": 4};
                        

RELAX.settings = {"GTR" : 1,
                  "LocalMG" : 1,
                  "Estimate GTR" : 1};

RELAX.timers  = {6,1};
RELAX.timers_d = {"Overall":0, 
                  "Preliminaries":1,
                  "General Descriptive":2,
                  "Null":3,
                  "Alternative":4,
                  "Partitioned Descriptive":5};
                  

relax.taskTimerStart (RELAX.timers_d["Overall"]);
relax.taskTimerStart (0);

RELAX.json    = { terms.json.input: {},
                  terms.json.fits : {},
                  terms.json.timers : {},
                  terms.json.test_results : None
                };

RELAX.test            = "RELAX.test";
RELAX.reference       = "RELAX.reference";
RELAX.unclassified    = "RELAX.unclassified";
RELAX.relaxation_parameter = "RELAX.K";

// move to terms?
term.RELAX.k          = "relaxation coefficient";

relax.terms.test = "Test";
relax.terms.reference = "Reference";



/*------------------------------------------------------------------------------*/


io.DisplayAnalysisBanner ({terms.io.info : "RELAX (a random effects test of selection relaxation)
                            uses a random effects branch-site model framework
                            to test whether a set of 'Test' branches evolves under relaxed
                            selection relative to a set of 'Reference' branches (R), as measured
                            by the relaxation parameter (K).",
                           terms.io.version : "1.0",
                           terms.io.reference : "RELAX: Detecting Relaxed Selection in a Phylogenetic Framework (2015). Mol Biol Evol 32 (3): 820-832",
                           terms.io.authors : "Sergei L Kosakovsky Pond, Ben Murrell, Steven Weaver and Temple iGEM / UCSD viral evolution group",
                           terms.io.contact : "spond@temple.edu",
                           terms.io.requirements : "in-frame codon alignment and a phylogenetic tree, with at least two groups of branches defined using the {} notation (one group can be defined as all unlabeled branches)"
                          } );


/*------------------------------------------------------------------------------
                     Input information and setup
------------------------------------------------------------------------------*/
relax.codon_data_info     = alignments.PromptForGeneticCodeAndAlignment ("RELAX.codon_data", "RELAX.codon_filter");
relax.sample_size         = relax.codon_data_info[terms.sites] * relax.codon_data_info[terms.sequences];

relax.name_mapping = relax.codon_data_info[utility.getGlobalValue("terms.json.name_mapping")];
    /**
        will contain "mapped" -> "original" associations with sequence names; or null if no mapping was necessary
    */

// SJS edited, this was never entered due to missing namespace
if (None == relax.name_mapping) { /** create a 1-1 mapping if nothing was done */
    relax.name_mapping = {};
    utility.ForEach (alignments.GetSequenceNames ("RELAX.codon_data"), "_value_", "`&relax.name_mapping`[_value_] = _value_");
}

relax.codon_data_info[terms.json_] = relax.codon_data_info[terms.file] + ".RELAX.json";
io.ReportProgressMessage ("RELAX", "Loaded an MSA with " + relax.codon_data_info[terms.sequences] + " sequences and " + relax.codon_data_info[terms.sites] + " codons from '" + relax.codon_data_info[terms.file] + "'");


relax.codon_lists = models.codon.MapCode (relax.codon_data_info[terms.genetic_code]);

_Genetic_Code = relax.codon_data_info[terms.genetic_code];
    /*

     hack to make PopulateModelMatrix work

    */

relax.codon_frequencies     = frequencies._aux.CF3x4(frequencies._aux.empirical.collect_data ("RELAX.codon_filter",3,1,1),
models.DNA.alphabet, relax.codon_lists[terms.sense_codons], relax.codon_lists[terms.stop_codons]);
 ("RELAX.codon_filter");


relax.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (relax.codon_data_info[utility.getGlobalValue("terms.json.partitions")], relax.name_mapping);


io.CheckAssertion("utility.Array1D (relax.partitions_and_trees) == 1", "RELAX only works on a single partition dataset");

relax.filter_specification = alignments.DefineFiltersForPartitions (relax.partitions_and_trees, "RELAX.codon_data" , "RELAX.codon_filter.", relax.codon_data_info);

relax.trees = utility.Map (relax.partitions_and_trees, "_partition_", '_partition_["tree"]');
relax.filter_names = utility.Map (relax.filter_specification, "_partition_", '_partition_["name"]');

relax.tree = relax.trees[0];

utility.SetEnvVariable ("VERBOSITY_LEVEL", 0);

relax.selected_branches = relax.io.defineBranchSets (relax.tree);

RELAX.has_unclassified = relax.selected_branches / RELAX.unclassified;

RELAX.branch_to_partiton = {};
utility.ForEachPair (relax.selected_branches, "_key_", "_value_", "utility.ForEach (utility.Keys(_value_), '_branch_', 'RELAX.branch_to_partiton[_branch_] = _key_')");


io.ReportProgressMessage ("RELAX", "Selected " + Abs (relax.selected_branches[RELAX.test]) + " branches as the test set: " + Join (",", Rows (relax.selected_branches[RELAX.test])));

ChoiceList  (RELAX.runModel,"Analysis type",1,NO_SKIP,
            "All", "[Default] Fit descriptive models AND run the relax test (4 models)",
            "Minimal", "Run only the RELAX test (2 models)"                              
            );

if (RELAX.runModel < 0) {
    return None;
}



/*------------------------------------------------------------------------------
                                Fit global GTR
------------------------------------------------------------------------------*/

relax.taskTimerStart(RELAX.timers_d["Preliminaries"]);
//relax.taskTimerStart (1);

if (RELAX.settings["GTR"]) {
    io.ReportProgressMessage ("RELAX", "Obtaining branch lengths under the GTR model");
    relax.gtr_results = estimators.FitGTR     ("RELAX.codon_filter", relax.tree, None);
    io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.gtr_results[terms.log_likelihood]);
    estimators.fixSubsetOfEstimates (relax.gtr_results, relax.gtr_results[terms.global]);
} else {
    relax.gtr_results = None;
}


/*------------------------------------------------------------------------------
                                Fit MG94
------------------------------------------------------------------------------*/


if (RELAX.settings["LocalMG"] && RELAX.runModel == 0) {
  io.ReportProgressMessage ("RELAX", "Obtaining  omega and branch length estimates under the local MG94xGTR model");
  relax.local_mg_results  = estimators.FitMGREV     (relax.filter_names, relax.trees, relax.codon_data_info [terms.genetic_code], {"model-type" : terms.local}, relax.gtr_results);
  io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.local_mg_results[terms.log_likelihood]]);
  estimators.fixSubsetOfEstimates (relax.local_mg_results, relax.local_mg_results[terms.global]);
} else {
  relax.local_mg_results = relax.gtr_results;
}

parameters.DeclareGlobal ("relax.codon_branch_scaler", None);



/*------------------------------------------------------------------------------
                                Fit Partitioned MG94xREV
------------------------------------------------------------------------------*/

io.ReportProgressMessage ("RELAX", "Obtaining omega and branch length estimates under the partitioned MG94xGTR model");
relax.mg_results  = estimators.FitMGREV     (relax.filter_names, relax.trees, relax.codon_data_info [terms.genetic_code],
                                             {"model-type" : terms.local, "partitioned-omega" : {"0" : RELAX.branch_to_partiton}, "proportional-branch-length-scaler": {"0" : "relax.codon_branch_scaler"}},
                                             relax.local_mg_results);

relax.taskTimerStop(RELAX.timers_d["Preliminaries"]);
//relax.taskTimerStop (1);

relax.mg_results_rate =
                     {"Reference"   : {{estimators.GetGlobalMLE (relax.mg_results, RELAX.reference),1}},
                      "Test"        : {{estimators.GetGlobalMLE (relax.mg_results, RELAX.test),1}}};



io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.mg_results[terms.log_likelihood]);



relax.json_store_lf (RELAX.json, "Partitioned MG94xREV",
                     relax.mg_results[terms.log_likelihood], relax.mg_results[terms.parameters] + 5,
                     RELAX.timers[1],
                     relax._aux.extract_branch_info ((relax.mg_results[terms.branch_length])[0], "relax.branch.length"),
                     relax._aux.extract_branch_info ((relax.mg_results[terms.branch_length])[0], "relax.branch.omega"),
                     relax._aux.extract_tree_length(relax.mg_results),
                     utility.Keys(utility.MatrixToDict(relax.mg_results[terms.trees_]))[0],
                     relax.mg_results_rate,
                     None,
                     "&omega;",
                    RELAX.display_orders["Partitioned MG94xREV"]);
                    
relax.json_spool (RELAX.json, relax.codon_data_info[terms.json_]);


/*------------------------------------------------------------------------------
                     Setup RELAX preliminaries
------------------------------------------------------------------------------*/

relax.taskTimerStart(RELAX.timers_d["General Descriptive"]);
//relax.taskTimerStart (2);

RELAX.model_assignment             = {};
RELAX.model_specification          = {};

RELAX.reference.model      = relax.io.define_a_bsrel_model (RELAX.reference, relax.codon_frequencies, estimators.GetGlobalMLE (relax.mg_results, RELAX.reference) ,1);
RELAX.model_assignment[RELAX.reference] = RELAX.reference.model["id"];
RELAX.model_specification[RELAX.reference.model["id"]] = RELAX.reference.model;

RELAX.test.model           = relax.io.define_a_bsrel_model (RELAX.test, relax.codon_frequencies, estimators.GetGlobalMLE (relax.mg_results, RELAX.test) ,1);
RELAX.model_assignment[RELAX.test] = RELAX.test.model["id"];
RELAX.model_specification[RELAX.test.model["id"]] = RELAX.test.model;

parameters.ConstrainSets (RELAX.reference.model [terms.omegas], RELAX.test.model [terms.omegas]);
parameters.ConstrainSets (RELAX.reference.model [terms.f], RELAX.test.model [terms.f]);

if (RELAX.has_unclassified) {
    RELAX.unclassified.model = relax.io.define_a_bsrel_model (RELAX.unclassified, relax.codon_frequencies, estimators.GetGlobalMLE (relax.mg_results, RELAX.test) ,1);
    RELAX.model_assignment[RELAX.unclassified] = RELAX.unclassified.model[terms.id];
    RELAX.model_specification[RELAX.unclassified.model[terms.id]] = RELAX.unclassified.model;

    parameters.ConstrainSets (RELAX.reference.model [terms.omegas], RELAX.unclassified.model [terms.omegas]);
    parameters.ConstrainSets (RELAX.reference.model [terms.f], RELAX.unclassified.model [terms.f]);

}

model.ApplyModelToTree          ("RELAX.tree", relax.tree, RELAX.model_assignment, relax.selected_branches);

ASSUME_REVERSIBLE_MODELS = 1;
LikelihoodFunction relax.LF = (RELAX.codon_filter, RELAX.tree);

global RELAX.branch_scaler = 4;
RELAX.proportional_constraint = "RELAX.branch_scaler";

if (RELAX.settings["Estimate GTR"] != 1) {
    estimators.fixSubsetOfEstimates   (relax.mg_results, relax.mg_results[terms.global]);
}

estimators.ApplyExistingEstimates ("relax.LF",  RELAX.model_specification, relax.mg_results, None);


utility.SetEnvVariable ("USE_LAST_RESULTS", 1);



/*------------------------------------------------------------------------------
                                Fit General Descriptive 
------------------------------------------------------------------------------*/
if (RELAX.runModel == 0) {

    io.ReportProgressMessage ("RELAX", "Two-stage fit of the general descriptive model (separate relaxation parameter for each branch)");

    OPTIMIZATION_PRECISION = 0.1;

    Optimize (relax.MLE.general_descriptive, relax.LF);

    OPTIMIZATION_PRECISION = 0.001;
    parameters.UnconstrainParameterSet ("relax.LF", {{terms.lf.local.constrained}});

    Optimize (relax.MLE.general_descriptive, relax.LF);
    io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.MLE.general_descriptive[1][0]);

    relax.general_descriptive = estimators.ExtractMLEs ("relax.LF", RELAX.model_specification);
    relax.add_scores (relax.general_descriptive, relax.MLE.general_descriptive);

    relax.taskTimerStop(RELAX.timers_d["General Descriptive"]);
    //relax.taskTimerStop (2);

    relax.json_store_lf (RELAX.json, "General Descriptive",
                         relax.general_descriptive[terms.log_likelihood], relax.general_descriptive[terms.parameters],
                         RELAX.timers[2],
                         relax._aux.extract_branch_info ((relax.general_descriptive[terms.branch_length])[0], "relax.branch.length"),
                         relax._aux.extract_branch_info ((relax.general_descriptive[terms.branch_length])[0], "relax.branch.local_k"),
                         relax._aux.extract_tree_length(relax.general_descriptive),
                         utility.Keys(utility.MatrixToDict(relax.general_descriptive[terms.trees_]))[0],
                         {"All" : relax.getRateDistribution (RELAX.reference.model, 1)},
                         None,
                         "k",
                         RELAX.display_orders["General Descriptive"]);
    relax.json_spool (RELAX.json, relax.codon_data_info[terms.json_]);

} else {
    parameters.UnconstrainParameterSet ("relax.LF", {{terms.lf.local.constrained}});
}





/*------------------------------------------------------------------------------
                                Fit Null model
------------------------------------------------------------------------------*/


if (RELAX.has_unclassified) {
    parameters.RemoveConstraint (RELAX.unclassified.model [terms.omegas]);
    parameters.RemoveConstraint (RELAX.unclassified.model [terms.f]);
}
relax.taskTimerStart(RELAX.timers_d["Null"]);
//relax.taskTimerStart (3);
io.ReportProgressMessage ("RELAX", "Fitting the RELAX null model");



RELAX.null = relax.define.null ("RELAX.tree", RELAX.reference.model, relax.selected_branches);

Optimize (relax.MLE.null, relax.LF);
io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.MLE.null[1][0]);


relax.null = estimators.ExtractMLEs ("relax.LF", RELAX.model_specification);
relax.add_scores (relax.null, relax.MLE.null);

relax.taskTimerStop(RELAX.timers_d["Null"]);
//relax.taskTimerStop (3);

relax.omega_distributions = {};
relax.omega_distributions[relax.terms.test] = relax.getRateDistribution (RELAX.test.model, 1);
relax.omega_distributions[relax.terms.reference] =  relax.getRateDistribution (RELAX.reference.model, 1);

if (RELAX.has_unclassified) {
    relax.omega_distributions ["Unclassified"] = relax.getRateDistribution (RELAX.unclassified.model, 1);
}


relax.json_store_lf (RELAX.json, "Null",
                     relax.null[terms.log_likelihood], relax.null[terms.parameters],
                     RELAX.timers[3],
                     relax._aux.extract_branch_info ((relax.null[terms.branch_length])[0], "relax.branch.length"),
                     relax._aux.extract_branch_info ((relax.null[terms.branch_length])[0], "relax.branch.local_k"),
                     relax._aux.extract_tree_length(relax.null),
                     utility.Keys(utility.MatrixToDict(relax.null[terms.trees_]))[0],
                     relax.omega_distributions,
                     1,
                     "k",
                    RELAX.display_orders["Null"]);

relax.json_spool (RELAX.json, relax.codon_data_info[terms.json_]);



/*------------------------------------------------------------------------------
                                Fit Alternative model
------------------------------------------------------------------------------*/

io.ReportProgressMessage ("RELAX", "Fitting the RELAX alternative model");

relax.taskTimerStart(RELAX.timers_d["Alternative"]);
//relax.taskTimerStart (4);
parameters.RemoveConstraint (RELAX.relaxation_parameter);
Optimize (relax.MLE.alt, relax.LF);
io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.MLE.alt[1][0] + ". Relaxation parameter K = " + Eval (RELAX.relaxation_parameter));

LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (relax.codon_data_info[terms.file] + ".alternative.fit", CLEAR_FILE, relax.LF);
LIKELIHOOD_FUNCTION_OUTPUT = 2;


relax.alt = estimators.ExtractMLEs ("relax.LF", RELAX.model_specification);
relax.add_scores (relax.alt, relax.MLE.alt);

RELAX.json[terms.json.test_results] = relax.runLRT (relax.alt[terms.log_likelihood], relax.null[terms.log_likelihood]);

io.ReportProgressMessage ("RELAX", "Likelihood ratio test for relaxation on Test branches, p = " + (RELAX.json [terms.json.test_results])[terms.json.p_value]);

relax.taskTimerStop(RELAX.timers_d["Alternative"]);
//relax.taskTimerStop  (4);


relax.omega_distributions[relax.terms.test] = relax.getRateDistribution (RELAX.null, Eval (RELAX.relaxation_parameter));
relax.omega_distributions[relax.terms.reference] =  relax.getRateDistribution (RELAX.null, 1);

if (RELAX.has_unclassified) {
    relax.omega_distributions ["Unclassified"] = relax.getRateDistribution (RELAX.unclassified.model, 1);
}

relax.json_store_lf (RELAX.json, "Alternative",
                     relax.alt[terms.log_likelihood], relax.alt[terms.parameters],
                     RELAX.timers[4],
                     relax._aux.extract_branch_info ((relax.alt[terms.branch_lengths])[0], "relax.branch.length"),
                     relax._aux.extract_branch_info ((relax.alt[terms.branch_length])[0], "relax.branch.local_k"),
                     relax._aux.extract_tree_length(relax.alt),
                     utility.Keys(utility.MatrixToDict(relax.alt[terms.trees_]))[0],
                     relax.omega_distributions,
                     Eval (RELAX.relaxation_parameter),
                     "k",
                    RELAX.display_orders["Alternative"]);



/*------------------------------------------------------------------------------
                                Fit Partitioned Descriptive
------------------------------------------------------------------------------*/



if (RELAX.runModel == 0) {

    relax.taskTimerStart(RELAX.timers_d["Partitioned Descriptive"]);
    //relax.taskTimerStart  (5);

    parameters.RemoveConstraint (RELAX.test.model [terms.relax.omegas]);
    parameters.RemoveConstraint (RELAX.test.model [terms.relax.f]);
    parameters.SetConstraint    (RELAX.relaxation_parameter, Eval (RELAX.relaxation_parameter), "");


    io.ReportProgressMessage ("RELAX", "Fitting the RELAX partitioned descriptive model");
    Optimize (relax.MLE.part.expl, relax.LF);
    io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.MLE.part.expl [1][0]);

    LIKELIHOOD_FUNCTION_OUTPUT = 7;
    fprintf (relax.codon_data_info[terms.file] + ".partitioned_descriptive.fit", CLEAR_FILE, relax.LF);
    LIKELIHOOD_FUNCTION_OUTPUT = 2;

    relax.part.expl = estimators.ExtractMLEs ("relax.LF", RELAX.model_specification);

    relax.omega_distributions[relax.terms.test]       = relax.getRateDistribution (RELAX.test.model,Eval (RELAX.relaxation_parameter));
    relax.omega_distributions[relax.terms.reference] =  relax.getRateDistribution (RELAX.reference.model, 1);

    if (RELAX.has_unclassified) {
        relax.omega_distributions ["Unclassified"] = relax.getRateDistribution (RELAX.unclassified.model, 1);
    }

    relax.add_scores (relax.part.expl, relax.MLE.part.expl);

    relax.taskTimerStop(RELAX.timers_d["Partitioned Descriptive"]);
    //relax.taskTimerStop  (5);
    relax.json_store_lf (RELAX.json, "Partitioned Descriptive",
                         relax.part.expl[terms.log_likelihood], relax.part.expl[terms.parameters],
                         RELAX.timers[5],
                         relax._aux.extract_branch_info ((relax.part.expl[terms.branch_length])[0], "relax.branch.length"),
                         None,
                         relax._aux.extract_tree_length(relax.part.expl),
                         utility.Keys(utility.MatrixToDict(relax.part.expl[terms.trees_]))[0],
                         relax.omega_distributions,
                         None,
                         "",
                        RELAX.display_orders["Partitioned Descriptive"]);
}


relax.taskTimerStop(RELAX.timers_d["Overall"]);
//relax.taskTimerStop  (0); 
 
(RELAX.json [terms.json.timers])["Overall"]                  = RELAX.timers[RELAX.timers_d["Overall"]];
(RELAX.json [terms.json.timers])["Preliminaries"]            = RELAX.timers[RELAX.timers_d["Preliminaries"]];
(RELAX.json [terms.json.timers])["General Descriptive"]      = RELAX.timers[RELAX.timers_d["General Descriptive"]];
(RELAX.json [terms.json.timers])["Null"]                     = RELAX.timers[RELAX.timers_d["Null"]];
(RELAX.json [terms.json.timers])["Alternative"]              = RELAX.timers[RELAX.timers_d["Alternative"]];
(RELAX.json [terms.json.timers])["Partitioned Descriptive"]  = RELAX.timers[RELAX.timers_d["Partitioned Descriptive"]];


/*
Add input information to the JSON
*/
(RELAX.json[terms.json.input])[terms.json.input.filename]  = relax.codon_data_info[terms.file];
(RELAX.json[terms.json.input])[terms.json.input.sequences]  = relax.codon_data_info[terms.sequences];
(RELAX.json[terms.json.input])[terms.json.input.sites]  = relax.codon_data_info[terms.sites];
(RELAX.json[terms.json.input])[terms.json.tree_string] = relax.tree [terms.trees.newick_with_lengths];   //RELAX.json ["tree"] = relax.tree ["string"];

/* 
partition field should be lists rather than the current situation. For consistency, it is also now named partitions (plural!).
    former code: RELAX.json ["partition"] = relax.selected_branches;
*/
RELAX.json[terms.json.partitions] = {relax.terms.test: utility.Keys(relax.selected_branches[RELAX.test]), 
                                     relax.terms.reference: utility.Keys(relax.selected_branches[RELAX.reference])};




relax.json_spool (RELAX.json, relax.codon_data_info[terms.json_]);

return RELAX.json;











//------------------------------------------------------------------------------
// HELPER FUNCTIONS FROM THIS POINT ON
//------------------------------------------------------------------------------


function relax.branch.length (branch_info) {
    return branch_info[terms.MLE];
}

function relax.branch.omega  (branch_info) {
    return parameters.NormalizeRatio ((branch_info[terms.nonsynonymous_rate])[terms.MLE], (branch_info[terms.synonymous_rate])[terms.MLE]);
}

function relax.branch.local_k  (branch_info) {
    return (branch_info[term.RELAX.k])[terms.MLE];
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
  relax.getRateInformation.rate_classes = Abs (model_description[terms.omegas]);
  relax.getRateInformation.omega_info = {relax.getRateInformation.rate_classes,2};

  for (relax.getRateInformation.k = 0; relax.getRateInformation.k < relax.getRateInformation.rate_classes; relax.getRateInformation.k += 1) {
    relax.getRateInformation.omega_info[relax.getRateInformation.k][0] = Eval ((model_description[terms.omegas])[relax.getRateInformation.k])^K;
    relax.getRateInformation.omega_info[relax.getRateInformation.k][1] = Eval ((model_description[terms.weights])[relax.getRateInformation.k]);
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
    parameters.RemoveConstraint ((general_model[terms.omegas])[2]);

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
        desc [terms.log_likelihood] = mles[1][0];
        desc [terms.parameters] = mles[1][1] + 9 + 5 * (RELAX.settings["Estimate GTR"] != 1);
    }
}

//------------------------------------------------------------------------------
function relax.runLRT (ha, h0) {
    return {terms.LR : 2*(ha-h0),
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

    model_parameters = {terms.parameters : {terms.global : {}, terms.local : {}}, terms.omegas : {}, terms.weights : {}, terms.f : {}, terms.rate_matrix : {}, terms.length : ""};

    model_parameters[terms.omegas] = parameters.GenerateSequentialNames ("`id`.omega",    3, "_");
    model_parameters[terms.f]      = parameters.GenerateSequentialNames ("`id`.aux_freq", 2, "_");

    parameters.DeclareGlobal    (model_parameters[terms.f], None);
    parameters.SetRange         (model_parameters[terms.f], terms.range01);

    parameters.DeclareGlobal    (model_parameters[terms.omegas], None);

    model_parameters[terms.weights] = parameters.helper.stick_breaking (model_parameters[terms.f], {{0.7,0.25,0.05}});

    relax.init_omegas = {{0.05,0.25,4}};
    relax.init_omegas = relax.init_omegas * (1/ parameters.Mean (relax.init_omegas, model_parameters[terms.weights], Abs (model_parameters[terms.omegas])));

    parameters.SetRange ((model_parameters[terms.omegas])[0], terms.range_almost_01);
    parameters.SetValue ((model_parameters[terms.omegas])[0], relax.init_omegas[0]);

    parameters.SetRange ((model_parameters[terms.omegas])[1], terms.range_almost_01);
    parameters.SetValue ((model_parameters[terms.omegas])[1], relax.init_omegas[1]);

    parameters.SetRange ((model_parameters[terms.omegas])[2], terms.range_gte1);
    parameters.SetValue ((model_parameters[terms.omegas])[2], relax.init_omegas[2]);

    if (do_local) {
        parameters.SetConstraint ((model_parameters[terms.omegas])[2], " 1/" + ((model_parameters[terms.omegas])[0]) + "/" + ((model_parameters[terms.omegas])[1]) , "");
        relax.io.define_a_bsrel_model_r = {terms.lower_bound : 1e-4, terms.upper_bound : 1};
        parameters.SetRange ((model_parameters[terms.omegas])[1], relax.io.define_a_bsrel_model_r);
    }


    local_k :< 50;

    relax.nuc = {4,3};
    for (k = 0; k < 4; k+=1) {
        for (k2 = 0; k2 < 3; k2 += 1) {
            relax.nuc [k][k2] = ((frequencies["bases"])[models.DNA.alphabet[k]])[k2];
        }
    }


    for (k = 1; k <= 3; k +=1) {
        ((model_parameters[terms.parameters])[terms.global])[relax.define_omega_term (k)] = (model_parameters[terms.omegas])[k-1];
        if (k < 3) {
            ((model_parameters[terms.parameters])[terms.global])[relax.define_weight_term (k)] = (model_parameters[terms.f])[k-1];
        }

        model_parameters[terms.rate_matrix] + ("Q_`id`_" + k);
        if (do_local) {
            PopulateModelMatrix			  ((model_parameters[terms.rate_matrix])[k-1],  relax.nuc, "t", "Min (1000," + (model_parameters[terms.omegas])[k-1] +"^local_k)", "");
        } else {
            PopulateModelMatrix			  ((model_parameters[terms.rate_matrix])[k-1],  relax.nuc, "t", (model_parameters[terms.omegas])[k-1], "");
        }
    }

    model_parameters[terms.id] = "`id`_model";
    model_parameters["length-expression"] = relax._aux.define_bsrel_model ("`id`_model", model_parameters[terms.rate_matrix], model_parameters[terms.weights], frequencies[terms.codons]);

    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("A","C")] = "AC";
    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("A","T")] = "AT";
    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("C","G")] = "CG";
    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("C","T")] = "CT";
    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("G","T")] = "GT";

    model_parameters["set-branch-length"] = "relax.aux.copy_branch_length";

    model_parameters["length parameter"] = "t";
    ((model_parameters[terms.parameters])[terms.local])[terms.timeParameter ()] = "t";
    if (do_local) {
        ((model_parameters[terms.parameters])[terms.local])[term.RELAX.k] = "local_k";
    }
    model_parameters["get-branch-length"] = "relax.aux.retrieve_branch_length";
    return model_parameters;
}

//------------------------------------------------------------------------------


function relax.aux.retrieve_branch_length (model, tree, node) {
    relax.aux.retrieve_branch_length.locals = Columns ((model_parameters[terms.parameters])[terms.local]);
    for (relax.aux.retrieve_branch_length.i = 0; relax.aux.retrieve_branch_length.i < Columns (relax.aux.retrieve_branch_length.locals); relax.aux.retrieve_branch_length.i += 1) {
        Eval (relax.aux.retrieve_branch_length.locals[relax.aux.retrieve_branch_length.i] + " = `tree`.`node`." + relax.aux.retrieve_branch_length.locals[relax.aux.retrieve_branch_length.i]);
    }
    return Eval (model["length-expression"]);
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

function relax.json_store_lf (json, name, ll, df, time, branch_length, branch_annotation, tree_length, tree_string, omega_distribution, K, annotation_tag, display_order) {

    (json[terms.json.fits])[name] = {terms.json.log_likelihood : ll,
                            terms.json.parameters              : df,
                            terms.json.AICc                    : relax.getIC (ll, df, relax.sample_size),
                            terms.json.runtime            : time,
                            terms.json.branch_lengths     : branch_length,
                            terms.json.branch_annotations  : branch_annotation,
                            terms.json.tree_length: tree_length,
                            terms.json.tree_string: tree_string,
                            terms.json.rate_distributions : omega_distribution,
                            "K" : K,
                            terms.json.annotation_tag : annotation_tag,
                            terms.json.display_order : display_order};

 }


//------------------------------------------------------------------------------

// Extract tree length from a fitted result dictionary
lfunction relax._aux.extract_tree_length(results){
    
    mles_of_interest = (results[^"terms.branch_length"])[0];

    tl = 0;
    for (i = 0; i < Abs(mles_of_interest); i+=1){
        t = utility.Keys(mles_of_interest)[i];
        tl += (mles_of_interest[t])[^"terms.MLE"];
    }
    return tl;
}
