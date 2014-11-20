RequireVersion ("2.220141023");


VERBOSITY_LEVEL				= 0;
LF_SMOOTHING_SCALER         = 0.1;

LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("CF3x4");
LoadFunctionLibrary("TreeTools");


// namespace 'utility' for convenience functions 
LoadFunctionLibrary("lib2014/UtilityFunctions.bf");

// namespace 'io' for interactive/datamonkey i/o functions
LoadFunctionLibrary("lib2014/IOFunctions.bf");

LoadFunctionLibrary ("lib2014/models/codon/MG_REV.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("lib2014/tasks/estimators.bf");


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


io.displayAnalysisBanner ({"info" : "RELAX (a random effects test of selection relaxation)
                            uses a random effects branch-site model framework 
                            to test whether a set of 'Test' branches evolves under relaxed 
                            selection relative to a set of 'Reference' branches (R), as measured
                            by the relaxation parameter (K).",
                           "version" : "1.00",
                           RELAX.reference : "In revision, preprint at xxx",
                           "authors" : "Sergei L Kosakovsky Pond, Ben Murrell, Steven Weaver and the UCSD viral evolution group",
                           "contact" : "spond@ucsd.edu",
                           "requirements" : "in-frame codon alignment and a phylogenetic tree, with at least two groups of branches defined using the {} notation (one group can be defined as all unlabeled branches)"         
                          } );




                  

relax.codon_data_info     = utility.promptForGeneticCodeAndAlignment ("RELAX.codon_data", "RELAX.codon_filter");
relax.sample_size         = relax.codon_data_info["sites"] * relax.codon_data_info["sequences"];


relax.codon_data_info["json"] = relax.codon_data_info["file"] + ".RELAX.json";
io.reportProgressMessage ("RELAX", "Loaded an MSA with " + relax.codon_data_info["sequences"] + " sequences and " + relax.codon_data_info["sites"] + " codons from '" + relax.codon_data_info["file"] + "'");

relax.codon_frequencies     = utility.defineFrequencies ("RELAX.codon_filter");
relax.tree 	  = utility.loadAnnotatedTopology (1);

relax.selected_branches = relax.io.defineBranchSets (relax.tree);
RELAX.has_unclassified = utility.array.find (Rows (relax.selected_branches), RELAX.unclassified) >= 0;

RELAX.json ["partition"] = relax.selected_branches;
RELAX.json ["tree"] = relax.tree ["string"];

io.reportProgressMessage ("RELAX", "Selected " + Abs (relax.selected_branches[RELAX.test]) + " branches as the test set: " + Join (",", Rows (relax.selected_branches[RELAX.test])));


ChoiceList  (RELAX.runModel,"Analysis type",1,NO_SKIP,
            "All", "[Default] Fit descriptive models AND run the relax test (4 models)",
            "Minimal", "Run only the RELAX test (2 models)"
            );

if (RELAX.runModel < 0) {
    return;
}

relax.taskTimerStart (1);

if (RELAX.settings["GTR"]) {
    io.reportProgressMessage ("RELAX", "Obtaining branch lengths under the GTR model");
    relax.gtr_results = estimators.fitGTR     ("RELAX.codon_filter", relax.tree, None);
    io.reportProgressMessage ("RELAX", "Log(L) = " + relax.gtr_results["LogL"]);
    estimators.fixSubsetOfEstimates (relax.gtr_results, relax.gtr_results["global"]);
} else {
    relax.gtr_results = None;
}

if (RELAX.settings["LocalMG"] && RELAX.runModel == 0) {
  io.reportProgressMessage ("RELAX", "Obtaining  omega and branch length estimates under the local MG94xGTR model");
  relax.local_mg_results  = estimators.fitMGREV     (relax.codon_data_info, relax.tree, {"model-type" : terms.local}, relax.gtr_results);
  io.reportProgressMessage ("RELAX", "Log(L) = " + relax.local_mg_results["LogL"]);
  estimators.fixSubsetOfEstimates (relax.local_mg_results, relax.local_mg_results["global"]);
} else {
  relax.local_mg_results = relax.gtr_results;
}

io.reportProgressMessage ("RELAX", "Obtaining omega and branch length estimates under the partitioned MG94xGTR model");
relax.mg_results  = estimators.fitMGREV     (relax.codon_data_info, relax.tree, {"model-type" : terms.local, "partitioned-omega" : 1, "fixed-branch-proportions": 1}, relax.local_mg_results);
relax.taskTimerStop (1);

relax.mg_results_rate =                      
                     {"Reference"   : {{relax.extract_global_MLE (relax.mg_results, RELAX.reference),1}},
                      "Test"        : {{relax.extract_global_MLE (relax.mg_results, RELAX.test),1}}};
                      



io.reportProgressMessage ("RELAX", "Log(L) = " + relax.mg_results["LogL"]);


relax.json_store_lf (RELAX.json, "Partitioned MG94xREV", 
                     relax.mg_results["LogL"], relax.mg_results["parameters"] + 5,
                     RELAX.timers[1],
                     relax._aux.extract_branch_info ((relax.mg_results["branch lengths"])[0], "relax.branch.length"),
                     relax._aux.extract_branch_info ((relax.mg_results["branch lengths"])[0], "relax.branch.omega"),
                     relax.mg_results_rate,
                     None,
                     "&omega;"
                    );
relax.json_spool (RELAX.json, relax.codon_data_info["json"]);

relax.taskTimerStart (2);

RELAX.model_assignment             = {};
RELAX.model_specification          = {};

RELAX.reference.model      = relax.io.define_a_bsrel_model (RELAX.reference, relax.codon_frequencies, relax.extract_global_MLE (relax.mg_results, RELAX.reference) ,1);
RELAX.model_assignment[RELAX.reference] = RELAX.reference.model["id"];
RELAX.model_specification[RELAX.reference.model["id"]] = RELAX.reference.model;

RELAX.test.model           = relax.io.define_a_bsrel_model (RELAX.test, relax.codon_frequencies, relax.extract_global_MLE (relax.mg_results, RELAX.test) ,1);
RELAX.model_assignment[RELAX.test] = RELAX.test.model["id"];
RELAX.model_specification[RELAX.test.model["id"]] = RELAX.test.model;

parameters.constrainSets (RELAX.reference.model ["omegas"], RELAX.test.model ["omegas"]);
parameters.constrainSets (RELAX.reference.model ["f"], RELAX.test.model ["f"]);

if (RELAX.has_unclassified) {
    RELAX.unclassified.model = relax.io.define_a_bsrel_model (RELAX.unclassified, relax.codon_frequencies, relax.extract_global_MLE (relax.mg_results, RELAX.test) ,1);
    RELAX.model_assignment[RELAX.unclassified] = RELAX.unclassified.model["id"];
    RELAX.model_specification[RELAX.unclassified.model["id"]] = RELAX.unclassified.model;
    
    parameters.constrainSets (RELAX.reference.model ["omegas"], RELAX.unclassified.model ["omegas"]);
    parameters.constrainSets (RELAX.reference.model ["f"], RELAX.unclassified.model ["f"]);
    
}

model.applyModelToTree          ("RELAX.tree", relax.tree, RELAX.model_assignment, relax.selected_branches);

ASSUME_REVERSIBLE_MODELS = 1;
LikelihoodFunction relax.LF = (RELAX.codon_filter, RELAX.tree);

global RELAX.branch_scaler = 4;
RELAX.proportional_constraint = "RELAX.branch_scaler";

if (RELAX.settings["Estimate GTR"] != 1) {
    estimators.fixSubsetOfEstimates   (relax.mg_results, relax.mg_results["global"]);
}

estimators.applyExistingEstimates ("relax.LF",  RELAX.model_specification, relax.mg_results);

USE_LAST_RESULTS       = 1;

if (RELAX.runModel == 0) {

    io.reportProgressMessage ("RELAX", "Two-stage fit of the general descriptive model (separate relaxation parameter for each branch)");

    OPTIMIZATION_PRECISION = 0.1;

    Optimize (relax.MLE.general_descriptive, relax.LF);

    OPTIMIZATION_PRECISION = 0.001;
    parameters.unconstrain_parameter_set ("relax.LF", {{terms.lf.local.constrained}});

    Optimize (relax.MLE.general_descriptive, relax.LF);
    io.reportProgressMessage ("RELAX", "Log(L) = " + relax.MLE.general_descriptive[1][0]);
    
    relax.general_descriptive = estimators.extractMLEs ("relax.LF", RELAX.model_specification);   
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
    parameters.unconstrain_parameter_set ("relax.LF", {{terms.lf.local.constrained}});
}



   
//io.spoolLF ("relax.LF", "/Volumes/home-raid/Desktop/explore", None);


if (RELAX.has_unclassified) {
    parameters.removeConstraint (RELAX.unclassified.model ["omegas"]);
    parameters.removeConstraint (RELAX.unclassified.model ["f"]);
}

relax.taskTimerStart (3);
io.reportProgressMessage ("RELAX", "Fitting the RELAX null model");


//io.spoolLF ("relax.LF", "/Volumes/home-raid/Desktop/null", None);

RELAX.null = relax.define.null ("RELAX.tree", RELAX.reference.model, relax.selected_branches);

Optimize (relax.MLE.null, relax.LF);
io.reportProgressMessage ("RELAX", "Log(L) = " + relax.MLE.null[1][0]);


relax.null = estimators.extractMLEs ("relax.LF", RELAX.model_specification);   
relax.add_scores (relax.null, relax.MLE.null);

relax.taskTimerStop (3);

relax.omega_distributions = {};
relax.omega_distributions["Test"] = relax.getRateDistribution (RELAX.test.model, 1);
relax.omega_distributions["Reference"] =  relax.getRateDistribution (RELAX.reference.model, 1);
                             
if (RELAX.has_unclassified) {
    relax.omega_distributions ["Unclassified"] = relax.getRateDistribution (RELAX.unclassified.model, 1);
}

//io.spoolLF ("relax.LF", "/Volumes/home-raid/Desktop/null", None);

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

io.reportProgressMessage ("RELAX", "Fitting the RELAX alternative model");

relax.taskTimerStart (4);
parameters.removeConstraint (RELAX.relaxation_parameter);
Optimize (relax.MLE.alt, relax.LF);
io.reportProgressMessage ("RELAX", "Log(L) = " + relax.MLE.alt[1][0] + ". Relaxation parameter K = " + Eval (RELAX.relaxation_parameter));

LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (relax.codon_data_info["file"] + ".alternative.fit", CLEAR_FILE, relax.LF);
LIKELIHOOD_FUNCTION_OUTPUT = 2;


relax.alt = estimators.extractMLEs ("relax.LF", RELAX.model_specification);   
relax.add_scores (relax.alt, relax.MLE.alt);

RELAX.json ["relaxation-test"] = relax.runLRT (relax.alt["LogL"], relax.null["LogL"]);

io.reportProgressMessage ("RELAX", "Likelihood ratio test for relaxation on Test branches, p = " + (RELAX.json ["relaxation-test"])["p"]);
 
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


if (RELAX.runModel == 0) {

    relax.taskTimerStart  (5);

    parameters.removeConstraint (RELAX.test.model ["omegas"]);
    parameters.removeConstraint (RELAX.test.model ["f"]);
    parameters.setConstraint    (RELAX.relaxation_parameter, Eval (RELAX.relaxation_parameter), "");


    io.reportProgressMessage ("RELAX", "Fitting the RELAX partitioned exploratory model");
    Optimize (relax.MLE.part.expl, relax.LF);
    io.reportProgressMessage ("RELAX", "Log(L) = " + relax.MLE.part.expl [1][0]);

    LIKELIHOOD_FUNCTION_OUTPUT = 7;
    fprintf (relax.codon_data_info["file"] + ".partitioned_descriptive.fit", CLEAR_FILE, relax.LF);
    LIKELIHOOD_FUNCTION_OUTPUT = 2;

    relax.part.expl = estimators.extractMLEs ("relax.LF", RELAX.model_specification);   

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

//------------------------------------------------------------------------------ 
// HELPER FUNCTIONS FROM THIS POINT ON
//------------------------------------------------------------------------------ 

function relax.extract_global_MLE (fit, id) {
    return ((fit["global"])[id])["MLE"];
}

function relax.branch.length (branch_info) {
    return branch_info["MLE"];
}

function relax.branch.omega  (branch_info) {
    return parameters.normalize_ratio ((branch_info[terms.nonsynonymous_rate])["MLE"], (branch_info[terms.synonymous_rate])["MLE"]);
}

function relax.branch.local_k  (branch_info) {
    return (branch_info[term.RELAX.k])["MLE"];
}

function relax._aux.extract_branch_info.callback (key, value) {
    relax._aux.extract_branch_info_result [key] = utility.callFunction (callback, {"0" : "value"});
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
function relax.define.null._aux (key, value) {
    //fprintf (stdout, "`tree_id`.`key`.`relax.define.null.local` := `relax.define.null.global`\n");
    ExecuteCommands ("`tree_id`.`key`.`relax.define.null.local` := `relax.define.null.global`");  
}

//------------------------------------------------------------------------------ 
function relax.define.null (tree_id, general_model, partition) {
    RELAX.null = general_model;
    parameters.removeConstraint ((general_model["omegas"])[2]);
    
    relax.define.null.local = ((general_model["parameters"])["local"])[term.RELAX.k];

    //fprintf (stdout, "\n**", relax.define.null.par , "\n");
    relax.define.null.global = "1";
    (partition[RELAX.reference])["relax.define.null._aux"][""];
    
    parameters.setConstraint (RELAX.relaxation_parameter, 1, terms.global);
    relax.define.null.global = RELAX.relaxation_parameter;

    (partition[RELAX.test])["relax.define.null._aux"][""];
    
    ((general_model["parameters"])["global"])[RELAX.relaxation_parameter] = RELAX.relaxation_parameter;
    
    return RELAX.null;
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
function relax.io.define_a_bsrel_model (id, frequencies, mean_omega, do_local) {

    model_parameters = {"parameters" : {"global" : {}, "local" : {}}, "omegas" : {}, "weights" : {}, "f" : {}, "Q" : {}, "length" : ""};
    
    model_parameters["omegas"] = parameters.generate_sequential_names ("`id`.omega",    3, "_");
    model_parameters["f"]      = parameters.generate_sequential_names ("`id`.aux_freq", 2, "_");
    
    parameters.declareGlobal    (model_parameters["f"], None);
    parameters.setRange         (model_parameters["f"], terms.range01);

    parameters.declareGlobal    (model_parameters["omegas"], None);
    
    
    model_parameters["weights"] = parameters.helper.stick_breaking (model_parameters["f"], {{0.7,0.25,0.05}});
    
    relax.init_omegas = {{0.05,0.25,4}};
    relax.init_omegas = relax.init_omegas * (1/ parameters.mean (relax.init_omegas, model_parameters["weights"], Abs (model_parameters["omegas"])));

    parameters.setRange ((model_parameters["omegas"])[0], terms.range01);
    parameters.set_value ((model_parameters["omegas"])[0], relax.init_omegas[0]);
    
    parameters.setRange ((model_parameters["omegas"])[1], terms.range01);
    parameters.set_value ((model_parameters["omegas"])[1], relax.init_omegas[1]);
    
    parameters.setRange ((model_parameters["omegas"])[2], terms.range_gte1);
    parameters.set_value ((model_parameters["omegas"])[2], relax.init_omegas[2]);
    
    if (do_local) {
        parameters.setConstraint ((model_parameters["omegas"])[2], " 1/" + ((model_parameters["omegas"])[0]) + "/" + ((model_parameters["omegas"])[1]) , "");
        relax.io.define_a_bsrel_model_r = {"LB" : 1e-4, "UB" : 1};        
        parameters.setRange ((model_parameters["omegas"])[1], relax.io.define_a_bsrel_model_r);
    }
    

   local_k :< 50;
      
    for (k = 1; k <= 3; k +=1) {
        ((model_parameters["parameters"])["global"])[relax.define_omega_term (k)] = (model_parameters["omegas"])[k-1];
        if (k < 3) {
            ((model_parameters["parameters"])["global"])[relax.define_weight_term (k)] = (model_parameters["f"])[k-1];
        }
    
        model_parameters["Q"] + ("Q_`id`_" + k);
        if (do_local) {
            PopulateModelMatrix			  ((model_parameters["Q"])[k-1],  frequencies["nucleotide"], "t", "Min (1000," + (model_parameters["omegas"])[k-1] +"^local_k)", "");        
        } else {
            PopulateModelMatrix			  ((model_parameters["Q"])[k-1],  frequencies["nucleotide"], "t", (model_parameters["omegas"])[k-1], "");
        }
    }
    
    model_parameters["id"] = "`id`_model";
    model_parameters["length-expression"] = relax._aux.define_bsrel_model ("`id`_model", model_parameters["Q"], model_parameters["weights"], frequencies["codon"]);
    
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

//------------------------------------------------------------------------------ 
function relax.io.defineBranchSets (relax.tree) {
    
    available_models        = {};
    branch_set              = {};
    
    
    for (k = 0; k < Columns (relax.tree["model_list"]); k += 1) {
        available_models  [(relax.tree["model_list"])[k]] = 0;
    }
    
    (relax.tree["model_map"])["relax._aux.io.countBranchSets"][""];
    
     
    list_models = Rows (available_models); // get keys    
    io.checkAssertion ("Columns (list_models) > 1", "The tree string must include at least one two sets of branches, at least one of which is annotated using {}");
    
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
                    break;
                }
            }
        }
        else {
            
            bgSet = 1-fgSet;
        }
        branch_set [relax.handle.unlabeled(selectTheseForTesting[bgSet][0])] = RELAX.reference;
        return_set [RELAX.reference] = {};
     }
         
    (relax.tree["model_map"])["relax._aux.io.mapBranchSets"][""];
    relax.tree["model_list"] = Columns (relax.tree["model_map"]);
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
