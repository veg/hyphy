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
//LoadFunctionLibrary("libv3/terms-json.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");


LoadFunctionLibrary("BranchSiteTemplate");
/*------------------------------------------------------------------------------
    BranchSiteTemplate Defines

    BuildCodonFrequencies (obsF);
    PopulateModelMatrix (ModelMatrixName&, EFV, synrateP, globalP, nonsynRateP);

------------------------------------------------------------------------------*/

relax.analysis_description = {terms.io.info : "RELAX (a random effects test of selection relaxation)
                            uses a random effects branch-site model framework
                            to test whether a set of 'Test' branches evolves under relaxed
                            selection relative to a set of 'Reference' branches (R), as measured
                            by the relaxation parameter (K).",
                           terms.io.version : "1.0",
                           terms.io.reference : "RELAX: Detecting Relaxed Selection in a Phylogenetic Framework (2015). Mol Biol Evol 32 (3): 820-832",
                           terms.io.authors : "Sergei L Kosakovsky Pond, Ben Murrell, Steven Weaver and Temple iGEM / UCSD viral evolution group",
                           terms.io.contact : "spond@temple.edu",
                           terms.io.requirements : "in-frame codon alignment and a phylogenetic tree, with at least two groups of branches defined using the {} notation (one group can be defined as all unlabeled branches)"
                          };

io.DisplayAnalysisBanner ( relax.analysis_description );
                          
                          
// More sensible for users to see test results first. Can easily update this.
relax.display_orders = { "Null": 0,
                         "Alternative": 1,
                         "Partitioned MG94xREV": 2,
                         "General Descriptive": 3,
                         "Partitioned Descriptive": 4};
                        

relax.settings = {"GTR" : 1,
                  "LocalMG" : 1,
                  "Estimate GTR" : 1};

relax.timers  = {6,1};
relax.timers_d = {"Overall":0, 
                  "Preliminaries":1,
                  "General Descriptive":2,
                  "Null":3,
                  "Alternative":4,
                  "Partitioned Descriptive":5};
                  

relax.taskTimerStart (relax.timers_d["Overall"]);
relax.taskTimerStart (0);

relax.json    = { terms.json.analysis: relax.analysis_description,
                  terms.json.PMID: "25540451",
                  terms.json.input: {},
                  terms.json.fits : {},
                  terms.json.timers : {},
                  terms.json.test_results : None
                };

relax.test_branches          = "relax.test";
relax.reference_branches     = "relax.reference";
relax.unclassified_branches  = "relax.unclassified";
relax.k                      = "relax.K";                  // The parameter itself
relax.model.relaxation_coefficient = "relaxation coefficient";   // The term used in model definition

relax.test_distribution           = "Test";
relax.reference_distribution      = "Reference";
relax.unclassified_distribution   = "Unclassified";


/*------------------------------------------------------------------------------*/




/*------------------------------------------------------------------------------
                     Input information and setup
------------------------------------------------------------------------------*/
relax.codon_data_info     = alignments.PromptForGeneticCodeAndAlignment ("relax.codon_data", "relax.codon_filter");
relax.sample_size         = relax.codon_data_info[terms.data.sites] * relax.codon_data_info[terms.data.sequences];

relax.name_mapping = relax.codon_data_info[utility.getGlobalValue("terms.data.name_mapping")];
    /**
        will contain "mapped" -> "original" associations with sequence names; or null if no mapping was necessary
    */

// SJS edited, this was never entered due to missing namespace
if (None == relax.name_mapping) { /** create a 1-1 mapping if nothing was done */
    relax.name_mapping = {};
    utility.ForEach (alignments.GetSequenceNames ("relax.codon_data"), "_value_", "`&relax.name_mapping`[_value_] = _value_");
}

relax.codon_data_info[terms.json.json] = relax.codon_data_info[terms.data.file] + ".RELAX.json";
io.ReportProgressMessage ("RELAX", "Loaded an MSA with " + relax.codon_data_info[terms.data.sequences] + " sequences and " + relax.codon_data_info[terms.data.sites] + " codons from '" + relax.codon_data_info[terms.data.file] + "'");


relax.codon_lists = models.codon.MapCode (relax.codon_data_info[terms.code]);

_Genetic_Code = relax.codon_data_info[terms.code];
    /*

     hack to make PopulateModelMatrix work

    */

relax.codon_frequencies     = frequencies._aux.CF3x4(frequencies._aux.empirical.collect_data ("relax.codon_filter",3,1,1),
models.DNA.alphabet, relax.codon_lists[terms.sense_codons], relax.codon_lists[terms.stop_codons]);
 ("relax.codon_filter");


relax.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (relax.codon_data_info[utility.getGlobalValue("terms.data.partitions")], relax.name_mapping);


io.CheckAssertion("utility.Array1D (relax.partitions_and_trees) == 1", "RELAX only works on a single partition dataset");

relax.filter_specification = alignments.DefineFiltersForPartitions (relax.partitions_and_trees, "relax.codon_data" , "relax.codon_filter.", relax.codon_data_info);

relax.trees = utility.Map (relax.partitions_and_trees, "_partition_", '_partition_[terms.data.tree]');
relax.filter_names = utility.Map (relax.filter_specification, "_partition_", '_partition_[terms.data.name]');

relax.tree = relax.trees[0];


utility.SetEnvVariable ("VERBOSITY_LEVEL", 0);

relax.selected_branches = relax.io.defineBranchSets (relax.tree);

relax.has_unclassified = relax.selected_branches / relax.unclassified_branches;

relax.branch_to_partition = {};
utility.ForEachPair (relax.selected_branches, "_key_", "_value_", "utility.ForEach (utility.Keys(_value_), '_branch_', 'relax.branch_to_partition[_branch_] = _key_')");


io.ReportProgressMessage ("RELAX", "Selected " + Abs (relax.selected_branches[relax.test_branches]) + " branches as the test set: " + Join (",", Rows (relax.selected_branches[relax.test_branches])));

ChoiceList  (relax.runModel,"Analysis type",1,NO_SKIP,
            "All", "[Default] Fit descriptive models AND run the relax test (4 models)",
            "Minimal", "Run only the RELAX test (2 models)"                              
            );

if (relax.runModel < 0) {
    return None;
}


/* Add input and partition information to JSON */
(relax.json[terms.json.input])[terms.json.file]  = relax.codon_data_info[terms.data.file];
(relax.json[terms.json.input])[terms.json.sequences]  = relax.codon_data_info[terms.data.sequences];
(relax.json[terms.json.input])[terms.json.sites]  = relax.codon_data_info[terms.data.sites];
(relax.json[terms.json.input])[terms.json.tree_string] = relax.tree[terms.trees.newick_with_lengths];   //relax.json ["tree"] = relax.tree ["string"];

relax.partitions = {};
utility.ForEach(utility.Keys(relax.selected_branches[relax.test_branches]), "_value_", "`&relax.partitions`[_value_] = relax.test_distribution");
utility.ForEach(utility.Keys(relax.selected_branches[relax.reference_branches]), "_value_", "`&relax.partitions`[_value_] = relax.reference_distribution");
if (relax.has_unclassified){
    utility.ForEach(utility.Keys(relax.selected_branches[relax.unclassified_branches]), "_value_", "`&relax.partitions`[_value_] = relax.unclassified_distribution");
}
relax.json[terms.json.partitions] = relax.partitions;
relax.json_spool (relax.json, relax.codon_data_info[terms.json.json]);


/*------------------------------------------------------------------------------
                                Fit global GTR
------------------------------------------------------------------------------*/

relax.taskTimerStart(relax.timers_d["Preliminaries"]);
//relax.taskTimerStart (1);

if (relax.settings["GTR"]) {
    io.ReportProgressMessage ("RELAX", "Obtaining branch lengths under the GTR model");
    relax.gtr_results = estimators.FitGTR     ("relax.codon_filter", relax.tree, None);
    io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.gtr_results[terms.fit.log_likelihood]);
    estimators.fixSubsetOfEstimates (relax.gtr_results, relax.gtr_results[terms.global]);
} else {
    relax.gtr_results = None;
}


/*------------------------------------------------------------------------------
                                Fit MG94
------------------------------------------------------------------------------*/

if (relax.settings["LocalMG"] && relax.runModel == 0) {

  io.ReportProgressMessage ("RELAX", "Obtaining  omega and branch length estimates under the local MG94xGTR model");
  relax.local_mg_results  = estimators.FitMGREV(relax.filter_names, relax.trees, relax.codon_data_info [terms.code], {terms.run_options.model_type : terms.local}, relax.gtr_results);
  io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.local_mg_results[terms.fit.log_likelihood]);
  estimators.fixSubsetOfEstimates (relax.local_mg_results, relax.local_mg_results[terms.global]);
} 
else {
  relax.local_mg_results = relax.gtr_results;
}

parameters.DeclareGlobal ("relax.codon_branch_scaler", None);



/*------------------------------------------------------------------------------
                                Fit Partitioned MG94xREV
------------------------------------------------------------------------------*/

io.ReportProgressMessage ("RELAX", "Obtaining omega and branch length estimates under the partitioned MG94xGTR model");
relax.mg_results  = estimators.FitMGREV     (relax.filter_names, relax.trees, relax.codon_data_info [terms.code],
                                             {terms.run_options.model_type : terms.local, terms.run_options.partitioned_omega : {"0" : relax.branch_to_partition}, terms.run_options.proportional_branch_length_scaler: {"0" : "relax.codon_branch_scaler"}},
                                             relax.local_mg_results);

relax.taskTimerStop(relax.timers_d["Preliminaries"]);
//relax.taskTimerStop (1);

relax.mg_results_rate =
                     {relax.reference_distribution  : {{estimators.GetGlobalMLE (relax.mg_results, relax.reference_branches),1}},
                      relax.test_distribution        : {{estimators.GetGlobalMLE (relax.mg_results, relax.test_branches),1}}};



io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.mg_results[terms.fit.log_likelihood]);


relax.json_store_lf (relax.json, "Partitioned MG94xREV",
                     relax.mg_results[terms.fit.log_likelihood], relax.mg_results[terms.parameters] + 5,
                     relax.timers[1],
                     relax._aux.extract_branch_info ((relax.mg_results[terms.branch_length])[0], "relax.branch.length"),
                     relax._aux.extract_branch_info ((relax.mg_results[terms.branch_length])[0], "relax.branch.omega"),
                     relax._aux.extract_tree_length(relax.mg_results),
                     utility.Keys(utility.MatrixToDict(relax.mg_results[terms.fit.trees]))[0],
                     relax.mg_results_rate,
                     None,
                     "&omega;",
                    relax.display_orders["Partitioned MG94xREV"]);
                    
relax.json_spool (relax.json, relax.codon_data_info[terms.json.json]);


/*------------------------------------------------------------------------------
                     Setup RELAX preliminaries
------------------------------------------------------------------------------*/

relax.taskTimerStart(relax.timers_d["General Descriptive"]);
//relax.taskTimerStart (2);

relax.model_assignment             = {};
relax.model_specification          = {};

relax.reference.model      = relax.io.define_a_bsrel_model (relax.reference_branches, relax.codon_frequencies, estimators.GetGlobalMLE (relax.mg_results, relax.reference_branches) ,1);
relax.model_assignment[relax.reference_branches] = relax.reference.model[terms.id];
relax.model_specification[relax.reference.model[terms.id]] = relax.reference.model;

relax.test.model           = relax.io.define_a_bsrel_model (relax.test_branches, relax.codon_frequencies, estimators.GetGlobalMLE (relax.mg_results, relax.test_branches) ,1);
relax.model_assignment[relax.test_branches] = relax.test.model[terms.id];
relax.model_specification[relax.test.model[terms.id]] = relax.test.model;

parameters.ConstrainSets (relax.reference.model [terms.omegas], relax.test.model [terms.omegas]);
parameters.ConstrainSets (relax.reference.model [terms.freqs], relax.test.model [terms.freqs]);

if (relax.has_unclassified) {
    relax.unclassified.model = relax.io.define_a_bsrel_model (relax.unclassified_branches, relax.codon_frequencies, estimators.GetGlobalMLE (relax.mg_results, relax.test_branches) ,1);
    relax.model_assignment[relax.unclassified_branches] = relax.unclassified.model[terms.id];
    relax.model_specification[relax.unclassified.model[terms.id]] = relax.unclassified.model;

    parameters.ConstrainSets (relax.reference.model [terms.omegas], relax.unclassified.model [terms.omegas]);
    parameters.ConstrainSets (relax.reference.model [terms.freqs], relax.unclassified.model [terms.freqs]);

}

model.ApplyModelToTree          ("relax.tree", relax.tree, relax.model_assignment, relax.selected_branches);

ASSUME_REVERSIBLE_MODELS = 1;
LikelihoodFunction relax.LF = (relax.codon_filter, relax.tree);

global relax.branch_scaler = 4;
relax.proportional_constraint = "relax.branch_scaler";

if (relax.settings["Estimate GTR"] != 1) {
    estimators.fixSubsetOfEstimates   (relax.mg_results, relax.mg_results[terms.global]);
}

estimators.ApplyExistingEstimates ("relax.LF",  relax.model_specification, relax.mg_results, None);


utility.SetEnvVariable ("USE_LAST_RESULTS", 1);



/*------------------------------------------------------------------------------
                                Fit General Descriptive 
------------------------------------------------------------------------------*/
if (relax.runModel == 0) {

    io.ReportProgressMessage ("RELAX", "Two-stage fit of the general descriptive model (separate relaxation parameter for each branch)");

    OPTIMIZATION_PRECISION = 0.1;

    Optimize (relax.MLE.general_descriptive, relax.LF);

    OPTIMIZATION_PRECISION = 0.001;
    parameters.UnconstrainParameterSet ("relax.LF", {{terms.lf.local.constrained}});

    Optimize (relax.MLE.general_descriptive, relax.LF);
    io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.MLE.general_descriptive[1][0]);

    relax.general_descriptive = estimators.ExtractMLEs ("relax.LF", relax.model_specification);
    relax.add_scores (relax.general_descriptive, relax.MLE.general_descriptive);

    relax.taskTimerStop(relax.timers_d["General Descriptive"]);
    //relax.taskTimerStop (2);

    relax.json_store_lf (relax.json, "General Descriptive",
                         relax.general_descriptive[terms.fit.log_likelihood], relax.general_descriptive[terms.parameters],
                         relax.timers[2],
                         relax._aux.extract_branch_info ((relax.general_descriptive[terms.branch_length])[0], "relax.branch.length"),
                         relax._aux.extract_branch_info ((relax.general_descriptive[terms.branch_length])[0], "relax.branch.local_k"),
                         relax._aux.extract_tree_length(relax.general_descriptive),
                         utility.Keys(utility.MatrixToDict(relax.general_descriptive[terms.fit.trees]))[0],
                         {"All" : relax.getRateDistribution (relax.reference.model, 1)},
                         None,
                         "k",
                         relax.display_orders["General Descriptive"]);
    relax.json_spool (relax.json, relax.codon_data_info[terms.json.json]);

} else {
    parameters.UnconstrainParameterSet ("relax.LF", {{terms.lf.local.constrained}});
}





/*------------------------------------------------------------------------------
                                Fit Null model
------------------------------------------------------------------------------*/


if (relax.has_unclassified) {
    parameters.RemoveConstraint (relax.unclassified.model [terms.omegas]);
    parameters.RemoveConstraint (relax.unclassified.model [terms.freqs]);
}
relax.taskTimerStart(relax.timers_d["Null"]);
//relax.taskTimerStart (3);
io.ReportProgressMessage ("RELAX", "Fitting the RELAX null model");



relax.null.model = relax.define.null ("relax.tree", relax.reference.model, relax.selected_branches);


Optimize (relax.MLE.null, relax.LF);
io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.MLE.null[1][0]);


relax.null = estimators.ExtractMLEs ("relax.LF", relax.model_specification);


relax.add_scores (relax.null, relax.MLE.null);

relax.taskTimerStop(relax.timers_d["Null"]);
//relax.taskTimerStop (3);

relax.omega_distributions = {};
relax.omega_distributions[relax.test_distribution] = relax.getRateDistribution (relax.test.model, 1);
relax.omega_distributions[relax.reference_distribution] =  relax.getRateDistribution (relax.reference.model, 1);

if (relax.has_unclassified) {
    relax.omega_distributions [relax.unclassified_distribution] = relax.getRateDistribution (relax.unclassified.model, 1);
}


relax.json_store_lf (relax.json, "Null",
                     relax.null[terms.fit.log_likelihood], relax.null[terms.parameters],
                     relax.timers[3],
                     relax._aux.extract_branch_info ((relax.null[terms.branch_length])[0], "relax.branch.length"),
                     relax._aux.extract_branch_info ((relax.null[terms.branch_length])[0], "relax.branch.local_k"),
                     relax._aux.extract_tree_length(relax.null),
                     utility.Keys(utility.MatrixToDict(relax.null[terms.fit.trees]))[0],
                     relax.omega_distributions,
                     1,
                     "k",
                    relax.display_orders["Null"]);

relax.json_spool (relax.json, relax.codon_data_info[terms.json.json]);


/*------------------------------------------------------------------------------
                                Fit Alternative model
------------------------------------------------------------------------------*/

io.ReportProgressMessage ("RELAX", "Fitting the RELAX alternative model");

relax.taskTimerStart(relax.timers_d["Alternative"]);
//relax.taskTimerStart (4);
parameters.RemoveConstraint (relax.k);
Optimize (relax.MLE.alt, relax.LF);
io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.MLE.alt[1][0] + ". Relaxation parameter K = " + Eval (relax.k));

// Uncomment following lines to get the LF file
/*
LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (relax.codon_data_info[terms.data.file] + ".alternative.fit", CLEAR_FILE, relax.LF);
LIKELIHOOD_FUNCTION_OUTPUT = 2;
*/

relax.alt = estimators.ExtractMLEs ("relax.LF", relax.model_specification);
relax.add_scores (relax.alt, relax.MLE.alt);


relax.relaxation_test = relax.runLRT (relax.alt[terms.fit.log_likelihood], relax.null[terms.fit.log_likelihood]);

io.ReportProgressMessage ("RELAX", "Likelihood ratio test for relaxation on Test branches, p = " + (relax.relaxation_test)[terms.p_value]);

relax.taskTimerStop(relax.timers_d["Alternative"]);
//relax.taskTimerStop  (4);


relax.omega_distributions[relax.test_distribution] = relax.getRateDistribution (relax.null.model, Eval (relax.k));
relax.omega_distributions[relax.reference_distribution] =  relax.getRateDistribution (relax.null.model, 1);

if (relax.has_unclassified) {
    relax.omega_distributions [relax.unclassified_distribution] = relax.getRateDistribution (relax.unclassified.model, 1);
}

relax.json_store_lf (relax.json, "Alternative",
                     relax.alt[terms.fit.log_likelihood], relax.alt[terms.parameters],
                     relax.timers[4],
                     relax._aux.extract_branch_info ((relax.alt[terms.branch_length])[0], "relax.branch.length"),
                     relax._aux.extract_branch_info ((relax.alt[terms.branch_length])[0], "relax.branch.local_k"),
                     relax._aux.extract_tree_length(relax.alt),
                     utility.Keys(utility.MatrixToDict(relax.alt[terms.fit.trees]))[0],
                     relax.omega_distributions,
                     Eval (relax.k),
                     "k",
                    relax.display_orders["Alternative"]);



/*------------------------------------------------------------------------------
                                Fit Partitioned Descriptive
------------------------------------------------------------------------------*/



if (relax.runModel == 0) {

    relax.taskTimerStart(relax.timers_d["Partitioned Descriptive"]);
    //relax.taskTimerStart  (5);

    parameters.RemoveConstraint (relax.test.model [terms.omegas]);
    parameters.RemoveConstraint (relax.test.model [terms.freqs]);
    parameters.SetConstraint    (relax.k, Eval (relax.k), "");

    io.ReportProgressMessage ("RELAX", "Fitting the RELAX partitioned descriptive model");
    Optimize (relax.MLE.part.expl, relax.LF);
    io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.MLE.part.expl [1][0]);


    // Uncomment following lines to get the LF file
    /*
    LIKELIHOOD_FUNCTION_OUTPUT = 7;
    fprintf (relax.codon_data_info[terms.file] + ".partitioned_descriptive.fit", CLEAR_FILE, relax.LF);
    LIKELIHOOD_FUNCTION_OUTPUT = 2;
    */
    relax.part.expl = estimators.ExtractMLEs ("relax.LF", relax.model_specification);




    relax.omega_distributions[relax.test_distribution]       = relax.getRateDistribution (relax.test.model,Eval (relax.k));
    relax.omega_distributions[relax.reference_distribution] =  relax.getRateDistribution (relax.reference.model, 1);

    if (relax.has_unclassified) {
        relax.omega_distributions [relax.unclassified_distribution] = relax.getRateDistribution (relax.unclassified.model, 1);
    }

    relax.add_scores (relax.part.expl, relax.MLE.part.expl);

    relax.taskTimerStop(relax.timers_d["Partitioned Descriptive"]);
    //relax.taskTimerStop  (5);
    relax.json_store_lf (relax.json, "Partitioned Descriptive",
                         relax.part.expl[terms.fit.log_likelihood], relax.part.expl[terms.parameters],
                         relax.timers[5],
                         relax._aux.extract_branch_info ((relax.part.expl[terms.branch_length])[0], "relax.branch.length"),
                         None,
                         relax._aux.extract_tree_length(relax.part.expl),
                         utility.Keys(utility.MatrixToDict(relax.part.expl[terms.fit.trees]))[0],
                         relax.omega_distributions,
                         None,
                         "",
                        relax.display_orders["Partitioned Descriptive"]);
}


relax.taskTimerStop(relax.timers_d["Overall"]);
//relax.taskTimerStop  (0); 
 
(relax.json [terms.json.timers])["Overall"]                  = relax.timers[relax.timers_d["Overall"]];
(relax.json [terms.json.timers])["Preliminaries"]            = relax.timers[relax.timers_d["Preliminaries"]];
(relax.json [terms.json.timers])["General Descriptive"]      = relax.timers[relax.timers_d["General Descriptive"]];
(relax.json [terms.json.timers])["Null"]                     = relax.timers[relax.timers_d["Null"]];
(relax.json [terms.json.timers])["Alternative"]              = relax.timers[relax.timers_d["Alternative"]];
(relax.json [terms.json.timers])["Partitioned Descriptive"]  = relax.timers[relax.timers_d["Partitioned Descriptive"]];

// Results from RELAX test.
relax.json[terms.json.test_results] = {terms.LRT: relax.relaxation_test[terms.LRT], terms.p_value: relax.relaxation_test[terms.p_value]};



relax.json_spool (relax.json, relax.codon_data_info[terms.json.json]);


return relax.json;











//------------------------------------------------------------------------------
// HELPER FUNCTIONS FROM THIS POINT ON
//------------------------------------------------------------------------------


function relax.branch.length (branch_info) {
    return branch_info[terms.fit.MLE];
}

function relax.branch.omega  (branch_info) {
    return parameters.NormalizeRatio ((branch_info[terms.nonsynonymous_rate])[terms.fit.MLE], (branch_info[terms.synonymous_rate])[terms.fit.MLE]);
}

function relax.branch.local_k  (branch_info) {
    return (branch_info[relax.model.relaxation_coefficient])[terms.fit.MLE];
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
    relax.null.model = general_model;
    parameters.RemoveConstraint ((general_model[terms.omegas])[2]);

    relax.define.null.local = ((general_model[terms.parameters])[terms.local])[relax.model.relaxation_coefficient];

    //fprintf (stdout, "\n**", relax.define.null.par , "\n");
    relax.define.null.global = "1";
    (partition[relax.reference_branches])["relax.define.null._aux"][""];

    parameters.SetConstraint (relax.k, 1, terms.global);
    relax.define.null.global = relax.k;

    (partition[relax.test_branches])["relax.define.null._aux"][""];

    ((general_model[terms.parameters])[terms.global])[relax.k] = relax.k;

    return relax.null.model;
}



//------------------------------------------------------------------------------
function relax.add_scores (desc, mles) {
    if (Type (mles) == "Matrix") {
        desc [terms.fit.log_likelihood] = mles[1][0];
        desc [terms.parameters] = mles[1][1] + 9 + 5 * (relax.settings["Estimate GTR"] != 1);
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

    model_parameters = {terms.parameters : {terms.global : {}, terms.local : {}}, terms.omegas : {}, terms.weights : {}, terms.freqs : {}, terms.model.rate_matrix : {}, terms.model.length : ""};

    model_parameters[terms.omegas] = parameters.GenerateSequentialNames ("`id`.omega",    3, "_");
    model_parameters[terms.freqs]      = parameters.GenerateSequentialNames ("`id`.aux_freq", 2, "_");

    parameters.DeclareGlobal    (model_parameters[terms.freqs], None);
    parameters.SetRange         (model_parameters[terms.freqs], terms.range01);

    parameters.DeclareGlobal    (model_parameters[terms.omegas], None);

    model_parameters[terms.weights] = parameters.helper.stick_breaking (model_parameters[terms.freqs], {{0.7,0.25,0.05}});

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
            relax.nuc [k][k2] = ((frequencies[terms.bases])[models.DNA.alphabet[k]])[k2];
        }
    }


    for (k = 1; k <= 3; k +=1) {
        ((model_parameters[terms.parameters])[terms.global])[relax.define_omega_term (k)] = (model_parameters[terms.omegas])[k-1];
        if (k < 3) {
            ((model_parameters[terms.parameters])[terms.global])[relax.define_weight_term (k)] = (model_parameters[terms.freqs])[k-1];
        }

        model_parameters[terms.model.rate_matrix] + ("Q_`id`_" + k);
        if (do_local) {
            PopulateModelMatrix			  ((model_parameters[terms.model.rate_matrix])[k-1],  relax.nuc, "t", "Min (1000," + (model_parameters[terms.omegas])[k-1] +"^local_k)", "");
        } else {
            PopulateModelMatrix			  ((model_parameters[terms.model.rate_matrix])[k-1],  relax.nuc, "t", (model_parameters[terms.omegas])[k-1], "");
        }
    }

    model_parameters[terms.id] = "`id`_model";
    model_parameters[terms.model.length_expression] = relax._aux.define_bsrel_model ("`id`_model", model_parameters[terms.model.rate_matrix], model_parameters[terms.weights], frequencies[terms.codons]);

    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("A","C")] = "AC";
    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("A","T")] = "AT";
    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("C","G")] = "CG";
    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("C","T")] = "CT";
    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("G","T")] = "GT";

    model_parameters[terms.model.set_branch_length] = "relax.aux.copy_branch_length";

    model_parameters[terms.model.length_parameter] = terms.default_time;
    ((model_parameters[terms.parameters])[terms.local])[terms.timeParameter ()] = terms.default_time;
    if (do_local) {
        ((model_parameters[terms.parameters])[terms.local])[relax.model.relaxation_coefficient] = "local_k";
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
    relax.aux.copy_branch_length.k = ((model[terms.parameters])[terms.local])[relax.model.relaxation_coefficient];

    if (Abs (relax.proportional_constraint)) {
        Eval ("`parameter`.`relax.aux.copy_branch_length.t` := `relax.proportional_constraint` * " + value);
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
        value = relax.unclassified_branches;
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
        branch_set [relax.handle.unlabeled(selectTheseForTesting[fgSet][0])] = relax.test_branches;
        return_set [relax.test_branches] = {};
        if (option_count > 2) {
            ChoiceList  (bgSet,"Choose the set of reference branches (R set)",1,fgSet,selectTheseForTesting);
            if (bgSet < 0) {
                return {};
            }
            for (k = 0; k < option_count; k+=1) {
                if (k != bgSet && k != fgSet) {
                    branch_set [relax.handle.unlabeled(selectTheseForTesting[k][0])] = relax.unclassified_branches;
                    return_set [relax.unclassified_branches] = {};
                }
            }
        }
        else {
            bgSet = 1-fgSet;
        }
        branch_set [relax.handle.unlabeled(selectTheseForTesting[bgSet][0])] = relax.reference_branches;
        return_set [relax.reference_branches] = {};
     }

    (relax.tree[terms.trees.model_map])["relax._aux.io.mapBranchSets"][""];
    relax.tree[terms.trees.model_list] = Columns (relax.tree[terms.trees.model_map]);
    //fprintf (stdout, "\n", relax.tree, "\n", return_set, "\n");

    return return_set;

}


//------------------------------------------------------------------------------------------------------------------------

function relax.taskTimerStart (index) {
    relax.timers[index] = Time(1);
}

function relax.taskTimerStop (index) {
    relax.timers[index] = Time(1) - relax.timers[index];
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
                            terms.parameters              : df,
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
    
    mles_of_interest = (results[utility.getGlobalValue("terms.branch_length")])[0];

    tl = 0;
    for (i = 0; i < Abs(mles_of_interest); i+=1){
        t = utility.Keys(mles_of_interest)[i];
        tl += (mles_of_interest[t])[utility.getGlobalValue("terms.fit.MLE")];
    }
    return tl;
}
