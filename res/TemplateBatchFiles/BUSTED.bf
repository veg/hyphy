RequireVersion ("2.1320140810");

_BUSTED_timers  = {3,1};
busted.taskTimerStart (2);

VERBOSITY_LEVEL				= 0;
LF_SMOOTHING_SCALER         = 0.01;


LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("dSdNTreeTools");
LoadFunctionLibrary("CF3x4");

// namespace 'utility' for convenience functions 
LoadFunctionLibrary("lib2014/UtilityFunctions.bf");

// namespace 'io' for interactive/datamonkey i/o functions
LoadFunctionLibrary("lib2014/IOFunctions.bf");

// namespace 'models.DNA.GTR' for the nucleotide GTR model
LoadFunctionLibrary("lib2014/tasks/estimators.bf");


io.displayAnalysisBanner ({"info" : "BUSTED (branch-site unrestricted statistical test of episodic diversification)
                            uses a random effects branch-site model fitted jointly to all or a subset of tree branches
                            in order to test for alignment-wide evidence of episodic diversifying selection. Assuming
                            there is evidence of positive selection (i.e. there is an omega > 1), BUSTED will also perform
                            a quick evidence-ratio style analysis to explore which individual sites may have been subject to selection.",
                           "version" : "1.00",
                           "reference" : "In preparation, preprint at xxx",
                           "authors" : "Sergei L Kosakovsky Pond and the UCSD VEG group (expand)",
                           "contact" : "spond@ucsd.edu",
                           "requirements" : "in-frame codon alignment and a phylogenetic tree"         
                          } );



/*------------------------------------------------------------------------------ 
    BranchSiteTemplate Defines

    BuildCodonFrequencies (obsF);
    PopulateModelMatrix (ModelMatrixName&, EFV, synrateP, globalP, nonsynRateP);

------------------------------------------------------------------------------*/


LoadFunctionLibrary("BranchSiteTemplate");



_BUSTED_json    = {"fits" : {},
                  "timers" : {},
                  "profiles" : {}
                  };
                  

codon_data_info = utility.promptForGeneticCodeAndAlignment ("codon_data", "codon_filter");

codon_data_info["json"] = codon_data_info["file"] + ".BUSTED.json";

io.reportProgressMessage ("BUSTED", "Loaded an MSA with " + codon_data_info["sequences"] + " sequences and " + codon_data_info["sites"] + " codons from '" + codon_data_info["file"] + "'");

codon_frequencies     = utility.defineFrequencies ("codon_filter");
tree_definition 	  = utility.loadAnnotatedTopology ();

busted.selected_branches = busted.io.selectBranchesToTest (tree_definition);

_BUSTED_json ["test set"] = busted.selected_branches;

io.reportProgressMessage ("BUSTED", "Selected " + Abs (busted.selected_branches) + " branches as the test (foreground) set: " + Join (",", Rows (busted.selected_branches)));

busted.model_definitions = busted.io.define_bsrel_models  ("FG","BG", codon_frequencies);

io.reportProgressMessage ("BUSTED", "Obtaining initial branch lengths under the GTR model");
busted.gtr_results = estimators.fitGTR     ("codon_filter", tree_definition, None);
io.reportProgressMessage ("BUSTED", "Log(L) = " + busted.gtr_results["LogL"]);

model.applyModelToTree ("busted.tree", tree_definition, "", {"DEFAULT" : (busted.model_definitions["BG"])["model"], 
                                                             (busted.model_definitions["FG"])["model"] : Rows (busted.selected_branches)});
                                                             
                                                             
                                                             
busted.taskTimerStart (0);
LikelihoodFunction busted.LF = (codon_filter, busted.tree);

_BUSTED_json["background"] =  busted.hasBackground  ("busted.tree");

estimators.applyExistingEstimates ("busted.LF", {"0":busted.model_definitions}, busted.gtr_results);
global busted.T_scaler = 1;
gtr_lengths = (busted.gtr_results["branch lengths"])[0];
gtr_lengths ["busted._aux.copy_lengths"][""];

io.reportProgressMessage ("BUSTED", "Fitting the unconstrained branch-site model");

USE_LAST_RESULTS = 1;
OPTIMIZATION_PRECISION = 0.1;

busted.computeSiteLikelihoods ("busted.LF");

Optimize (busted.MLE_HA, busted.LF);
gtr_lengths ["busted._aux.free_lengths"][""];

OPTIMIZATION_PRECISION = 0.001;
Optimize (busted.MLE_HA, busted.LF);
busted_positive_class = busted.checkForPS (busted.model_definitions);
io.reportProgressMessage ("BUSTED", "Log(L) = " + busted.MLE_HA[1][0] + ". Unrestricted class omega = " + busted_positive_class["omega"] + " (weight = " + busted_positive_class["weight"] + ")");


busted.sample_size             =codon_data_info["sites"] * codon_data_info["sequences"];
busted.taskTimerStop (0);

busted.json_store_lf                (_BUSTED_json, "Unconstrained model", busted.MLE_HA[1][0], busted.MLE_HA[1][1]+9, busted.getIC (busted.MLE_HA[1][0], busted.MLE_HA[1][1]+9, busted.sample_size) , _BUSTED_timers[0], 0, 0);
// +BranchLength (givenTree,-1), Format (givenTree, 1,1)); 


busted.profiles = {};
busted.profiles ["HA"] = busted.computeSiteLikelihoods ("busted.LF");
(_BUSTED_json ["profiles"])["unconstrained"] = busted.matrix2json (busted.profiles ["HA"]);


if (busted_positive_class["omega"] < 1 || busted_positive_class["weight"] < 1e-8) {
    io.reportProgressMessage ("BUSTED", "No evidence for positive selection under the unconstrained model, skipping constrained model fitting");
    return 0; 
} else {
    busted.taskTimerStart (1);

    io.reportProgressMessage ("BUSTED", "Fitting the branch-site model that disallows omega > 1 among foreground branches");
    busted.constrainTheModel (busted.model_definitions);
    Optimize (busted.MLE_H0, busted.LF);
    io.reportProgressMessage ("BUSTED", "Log(L) = " + busted.MLE_H0[1][0]);
    busted.LRT = busted.runLRT (busted.MLE_HA[1][0], busted.MLE_H0[1][0]);
    
    _BUSTED_json ["test results"] = busted.LRT;
    
    io.reportProgressMessage ("BUSTED", "Likelihood ratio test for episodic positive selection, p = " + busted.LRT["p"]);
    busted.profiles ["H0"] = busted.computeSiteLikelihoods ("busted.LF");
    (_BUSTED_json ["profiles"])["constrained"] = busted.matrix2json (busted.profiles ["H0"]);
    busted.taskTimerStop (1);
    busted.json_store_lf                (_BUSTED_json, "Constrained model", busted.MLE_H0[1][0], busted.MLE_H0[1][1]+9, busted.getIC (busted.MLE_H0[1][0], busted.MLE_H0[1][1]+9, busted.sample_size) , _BUSTED_timers[1], 0, 0);
    (_BUSTED_json ["profiles"])["er"] = busted.evidenceRatios ( busted.profiles ["HA"],  busted.profiles ["H0"]);

}

busted.taskTimerStop (2);

(_BUSTED_json ["timers"])["overall"] = _BUSTED_timers[2];
(_BUSTED_json ["timers"])["unconstrained"] = _BUSTED_timers[0];
(_BUSTED_json ["timers"])["constrained"] = _BUSTED_timers[1];

fprintf (codon_data_info["json"], CLEAR_FILE, _BUSTED_json);


//------------------------------------------------------------------------------ 
// HELPER FUNCTIONS FROM HTHIS POINT ON
//------------------------------------------------------------------------------ 
function busted._aux.copy_lengths (key, value) {
    ExecuteCommands ("busted.tree." + key + ".t := busted.T_scaler * " + value["MLE"]);
}

//------------------------------------------------------------------------------ 
function busted.hasBackground (id) {
   ExecuteCommands ("GetInformation (busted.nodeMap, `id`)");
   busted.nodeMap = Columns(busted.nodeMap);
   return Rows (busted.nodeMap) * Columns (busted.nodeMap) > 1; 
}


//------------------------------------------------------------------------------ 
function busted._aux.free_lengths (key, value) {
    ExecuteCommands ("busted.tree." + key + ".t = busted.tree." + key + ".t");
}


//------------------------------------------------------------------------------ 
function busted.checkForPS (model_parameters) {
   return {"omega" :Eval (((model_parameters["FG"])["omegas"])[2]),
           "weight" : Eval (((model_parameters["FG"])["weights"])[2])};
           
}

//------------------------------------------------------------------------------ 
function busted.constrainTheModel (model_parameters) {
   Eval (((model_parameters["FG"])["omegas"])[2] + ":<1");           
}

//------------------------------------------------------------------------------ 
function busted.computeSiteLikelihoods (id) {
   ConstructCategoryMatrix (_siteLike, *id, SITE_LOG_LIKELIHOODS);
   return _siteLike;         
}



//------------------------------------------------------------------------------ 
function busted.runLRT (ha, h0) {
    return {"LR" : 2*(ha-h0),
            "p" : (1-CChi2 (2*(ha-h0),1))*.5};
}

//------------------------------------------------------------------------------ 
lfunction busted._aux.stick_breaking (parameters) {
    left_over = ""; 
    weights = {};
    for (k = 0; k < Abs (parameters); k += 1) {
        weights [k] = left_over + parameters[k];
        left_over += "(1-" + parameters[k] + ")*";
    }
    weights[k] = left_over[0][Abs (left_over)-2];
    return weights;
}

//------------------------------------------------------------------------------ 
function busted._aux.define_bsrel_model (id,Q,weights,freqs) {
    rate_count = Abs (Q);
    components = {};
    for (k = 0; k < rate_count; k+=1) {
        components[k] = "Exp(" + Q[k] + ")*" + weights[k];
    }
    Eval ("`id`_eqf = freqs");
    ExecuteCommands ("Model `id` =(\"" + Join("+",components) + "\",`id`_eqf,EXPLICIT_FORM_MATRIX_EXPONENTIAL);");
    return id;
}

//------------------------------------------------------------------------------ 
function busted._aux.define_parameter (key, value) {
   ExecuteCommands ("global `value` :< 1;");
   ExecuteCommands ("`value` :> 0;");
   ExecuteCommands ("`value` = " + busted.init_values[key]);
}


//------------------------------------------------------------------------------ 
function busted.io.define_bsrel_models (foreground_id, background_id, frequencies) {

    model_parameters = 
        {"FG": {"omegas" : {}, "weights" : {}, "f" : {}, "Q" : {}},
         "BG": {"omegas" : {}, "weights" : {}, "f" : {}, "Q" : {}},
         "parameters" : {"global" : {}}
          };
    
    for (k = 1; k <= 3; k +=1) {
        tag = "" + k;
        ((model_parameters["FG"])["omegas"]) + "`foreground_id`_omega_`tag`";
        ((model_parameters["BG"])["omegas"]) + "`background_id`_omega_`tag`";
        if (k < 3) {
            ((model_parameters["FG"])["f"]) + "`foreground_id`_f_`tag`";
            ((model_parameters["BG"])["f"]) + "`background_id`_f_`tag`";
        }
        
    }
    
    
    
    ((model_parameters["FG"])["weights"])  = busted._aux.stick_breaking (((model_parameters["FG"])["f"]));
    ((model_parameters["BG"])["weights"])  = busted._aux.stick_breaking (((model_parameters["BG"])["f"]));
     
    busted.init_values = {"0" : 0.1, "1" : 0.5, "2" : 1};

    ((model_parameters["FG"])["omegas"])["busted._aux.define_parameter"][""];
    ((model_parameters["BG"])["omegas"])["busted._aux.define_parameter"][""];
    
    Eval (((model_parameters["FG"])["omegas"])[2] + ":<1e26");
    Eval (((model_parameters["BG"])["omegas"])[2] + ":<1e26");

    busted.init_values = {"0" : 0.8, "1" : 0.75};
    
    ((model_parameters["FG"])["f"])["busted._aux.define_parameter"][""];
    ((model_parameters["BG"])["f"])["busted._aux.define_parameter"][""];
    

    for (k = 1; k <= 3; k +=1) {
        
        ((model_parameters["FG"])["Q"]) + ("Q_`foreground_id`_" + k);
        PopulateModelMatrix			  (((model_parameters["FG"])["Q"])[k-1],  frequencies["nucleotide"], "t",((model_parameters["FG"])["omegas"])[k-1], "");
        ((model_parameters["BG"])["Q"]) + ("Q_`background_id`_" + k);
        PopulateModelMatrix			  (((model_parameters["BG"])["Q"])[k-1],  frequencies["nucleotide"], "t",((model_parameters["BG"])["omegas"])[k-1], "");
    }
    
    (model_parameters["BG"])["model"] = busted._aux.define_bsrel_model ("`background_id`_model", (model_parameters["BG"])["Q"], (model_parameters["BG"])["weights"], frequencies["codon"]);
    (model_parameters["FG"])["model"] = busted._aux.define_bsrel_model ("`foreground_id`_model", (model_parameters["FG"])["Q"], (model_parameters["FG"])["weights"], frequencies["codon"]);
    
    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("A","C")] = {"ID" : "AC"};
    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("A","T")] = {"ID" : "AT"};
    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("C","G")] = {"ID" : "CG"};
    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("C","T")] = {"ID" : "CT"};
    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("G","T")] = {"ID" : "GT"};
    
    return model_parameters;
}


//------------------------------------------------------------------------------ 
lfunction busted.io.selectBranchesToTest (tree_definition) {
    
    UseModel (USE_NO_MODEL);
    ExecuteCommands ("Topology  T = " + tree_definition["string"]);
    tAVL = T^0;

    totalBranchCount = Abs (tAVL) - 2;

    selectedBranches = {};
    selectTheseForTesting = {totalBranchCount + 3, 2};
    selectTheseForTesting [0][0] = "All";  selectTheseForTesting [0][1] = "Test for selection on all branches jointly";
    selectTheseForTesting [1][0] = "Internal";  selectTheseForTesting [1][1] = "Test for selection on all internal branches jointly";
    selectTheseForTesting [2][0] = "Leaves";  selectTheseForTesting [2][1] = "Test for selection on all terminal branches jointly";

    for (k = 0; k < totalBranchCount; k += 1) {
        selectTheseForTesting [k+3][0] = (tAVL[k+1])["Name"];
        selectTheseForTesting [k+3][1] = "Add branch '" + selectTheseForTesting [k+3][0] + "' to the test set";

    }

    ChoiceList  (whichBranchesToTest,"Which branches to test?",0,NO_SKIP,selectTheseForTesting);

    for (k = 0; k < Columns (whichBranchesToTest); k += 1) {
       if (whichBranchesToTest [k] < 3) {
            for (k2 = 1; k2 <=  totalBranchCount; k2 += 1) {
                if (whichBranchesToTest[k] == 0 || whichBranchesToTest[k] == 1 && Abs ((tAVL[k2])["Children"]) > 0 || whichBranchesToTest[k] == 2 && Abs ((tAVL[k2])["Children"]) == 0) {
                    selectedBranches [(tAVL[k2])["Name"]] = 1;
                }
            }
            return selectedBranches;
        } 
        selectedBranches [whichBranchesToTest [k] - 3] = 1;
    }
    
    return selectedBranches;
}

//------------------------------------------------------------------------------------------------------------------------

function busted.evidenceRatios (ha, h0) {
    sites = Rows (ha) * Columns (ha);
    return ha["Exp(_MATRIX_ELEMENT_VALUE_-h0[_MATRIX_ELEMENT_COLUMN_])"];
}


//------------------------------------------------------------------------------------------------------------------------

function busted.taskTimerStart (index) {
    _BUSTED_timers[index] = Time(1);
}

function busted.taskTimerStop (index) {
    _BUSTED_timers[index] = Time(1) - _BUSTED_timers[index];

}

//------------------------------------------------------------------------------------------------------------------------

function busted.getIC (logl,params,samples) {
    return -2*logl + 2*samples/(samples-params-1)*params;
}

//------------------------------------------------------------------------------------------------------------------------

lfunction busted.matrix2json (mx) {
    result_str = ""; result_str * 128;
    result_str * "[";
    for (r = 0; r < Rows (mx); r+=1) {
        if (r) {
            result_str * ",";
        }
        result_str * ("[" + Join(",",mx[r][-1]) + "]");
    }
    result_str * "]";
    result_str * 0;
    return result_str;
}

//------------------------------------------------------------------------------------------------------------------------

lfunction busted.json_store_lf (json, name, ll, df, aicc, time, tree_length, tree_string) {
    (json["fits"])[name] = {"log-likelihood": ll,
                            "parameters": df,
                            "AIC-c" : aicc,
                            "runtime" : time,
                            "tree length" : tree_length,
                            "tree string" : tree_string,
                            "rate distributions" : {}};
}
