RequireVersion ("2.1320140810");

_RELAX_timers  = {5,1};
busted.taskTimerStart (4);

VERBOSITY_LEVEL				= 0;
LF_SMOOTHING_SCALER         = 0.01;


LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("CF3x4");
LoadFunctionLibrary("TreeTools");


// namespace 'utility' for convenience functions 
LoadFunctionLibrary("lib2014/UtilityFunctions.bf");

// namespace 'io' for interactive/datamonkey i/o functions
LoadFunctionLibrary("lib2014/IOFunctions.bf");

// namespace 'models.DNA.GTR' for the nucleotide GTR model
LoadFunctionLibrary("lib2014/tasks/estimators.bf");


io.displayAnalysisBanner ({"info" : "RELAX (a random effects test of selection relaxation)
                            uses a random effects branch-site model framework 
                            to test whether a set of test branches (T) evolves under relaxed 
                            selection relative to a set of reference branches (R), as measured
                            by the relaxation parameter (K).",
                           "version" : "1.00",
                           "reference" : "In revision, preprint at xxx",
                           "authors" : "Sergei L Kosakovsky Pond, Ben Murrell, Steven Weaver and the UCSD viral evolution group",
                           "contact" : "spond@ucsd.edu",
                           "requirements" : "in-frame codon alignment and a phylogenetic tree, with at least two groups of branches defined using the {} notation (one group can be defined as all unlabeled branches)"         
                          } );



/*------------------------------------------------------------------------------ 
    BranchSiteTemplate Defines

    BuildCodonFrequencies (obsF);
    PopulateModelMatrix (ModelMatrixName&, EFV, synrateP, globalP, nonsynRateP);

------------------------------------------------------------------------------*/



LoadFunctionLibrary("BranchSiteTemplate");



_RELAX_json    = {"fits" : {},
                  "timers" : {},
                  "tests" : {}
                  };
                  

codon_data_info = utility.promptForGeneticCodeAndAlignment ("RELAX.codon_data", "RELAX.codon_filter");

LoadFunctionLibrary ("lib2014/models/codon/MG_REV.bf");


codon_data_info["json"] = codon_data_info["file"] + ".RELAX.json";
io.reportProgressMessage ("RELAX", "Loaded an MSA with " + codon_data_info["sequences"] + " sequences and " + codon_data_info["sites"] + " codons from '" + codon_data_info["file"] + "'");

codon_frequencies     = utility.defineFrequencies ("RELAX.codon_filter");
tree_definition 	  = utility.loadAnnotatedTopology ();

relax.selected_branches = relax.io.defineBranchSets (tree_definition);
_RELAX_json ["partition"] = relax.selected_branches;

io.reportProgressMessage ("RELAX", "Selected " + Abs (relax.selected_branches["Test"]) + " branches as the test set: " + Join (",", Rows (relax.selected_branches["Test"])));

io.reportProgressMessage ("RELAX", "Obtaining branch lengths under the GTR model");
relax.gtr_results = estimators.fitGTR     ("RELAX.codon_filter", tree_definition, None);
io.reportProgressMessage ("RELAX", "Log(L) = " + relax.gtr_results["LogL"]);

io.reportProgressMessage ("RELAX", "Obtaining omega and branch length estimates under the MG94xGTR model using GTR as the starting point");
relax.mg_results  = estimators.fitMGREV     (codon_data_info, tree_definition, terms.local, relax.gtr_results);

return 0;

busted.model_definitions = busted.io.define_bsrel_models  ("FG","BG","Unclassified",codon_frequencies);

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

busted.bls = busted.io.evaluate_branch_lengths (busted.model_definitions, "busted.tree", busted.selected_branches);
Optimize (busted.MLE_HA, busted.LF);
gtr_lengths ["busted._aux.free_lengths"][""];

OPTIMIZATION_PRECISION = 0.001;
Optimize (busted.MLE_HA, busted.LF);
busted_positive_class = busted.checkForPS (busted.model_definitions);
io.reportProgressMessage ("BUSTED", "Log(L) = " + busted.MLE_HA[1][0] + ". Unrestricted class omega = " + busted_positive_class["omega"] + " (weight = " + busted_positive_class["weight"] + ")");


busted.sample_size             =codon_data_info["sites"] * codon_data_info["sequences"];
busted.taskTimerStop (0);

busted.bls = busted.io.evaluate_branch_lengths (busted.model_definitions, "busted.tree", busted.selected_branches);
busted.tavl         = busted.tree ^ 0;
busted.renderString = PostOrderAVL2StringDistances (busted.tavl, busted.bls);
UseModel (USE_NO_MODEL);
Tree busted.T = busted.renderString;

busted.json_store_lf                (_BUSTED_json, "Unconstrained model", 
        busted.MLE_HA[1][0], busted.MLE_HA[1][1]+9, 
        busted.getIC (busted.MLE_HA[1][0], busted.MLE_HA[1][1]+9, busted.sample_size) , 
        _BUSTED_timers[0], 
        +BranchLength (busted.T,-1), 
        Format (busted.T, 1,1),
        busted.model_definitions,
        _BUSTED_json["background"]
        );


busted.profiles = {};
(_BUSTED_json ["profiles"])["unconstrained"] = busted.computeSiteLikelihoods ("busted.LF");


if (busted_positive_class["omega"] < 1 || busted_positive_class["weight"] < 1e-8) {
    io.reportProgressMessage ("BUSTED", "No evidence for positive selection under the unconstrained model, skipping constrained model fitting");
    _BUSTED_json ["test results"] = busted.runLRT (0, 0);
} else {
    busted.taskTimerStart (1);

    io.reportProgressMessage ("BUSTED", "Fitting the branch-site model that disallows omega > 1 among foreground branches");
    busted.constrainTheModel (busted.model_definitions);
    (_BUSTED_json ["profiles"])["constrained"] = busted.computeSiteLikelihoods ("busted.LF");;
    Optimize (busted.MLE_H0, busted.LF);
    (_BUSTED_json ["profiles"])["optimized null"] = busted.computeSiteLikelihoods ("busted.LF");;
    io.reportProgressMessage ("BUSTED", "Log(L) = " + busted.MLE_H0[1][0]);
    busted.LRT = busted.runLRT (busted.MLE_HA[1][0], busted.MLE_H0[1][0]);
    
    _BUSTED_json ["test results"] = busted.LRT;
    
    io.reportProgressMessage ("BUSTED", "Likelihood ratio test for episodic positive selection, p = " + busted.LRT["p"]);
     busted.taskTimerStop (1);
    
    busted.bls = busted.io.evaluate_branch_lengths (busted.model_definitions, "busted.tree", busted.selected_branches);
    busted.tavl         = busted.tree ^ 0;
    busted.renderString = PostOrderAVL2StringDistances (busted.tavl, busted.bls);
    UseModel (USE_NO_MODEL);
    Tree busted.T = busted.renderString;

    busted.json_store_lf                (_BUSTED_json, 
                                        "Constrained model", busted.MLE_H0[1][0], 
                                        busted.MLE_H0[1][1]+9, 
                                        busted.getIC (busted.MLE_H0[1][0], busted.MLE_H0[1][1]+9, busted.sample_size) , 
                                        _BUSTED_timers[1],
                                         +BranchLength (busted.T,-1), 
                                        Format (busted.T, 1,1),
                                        busted.model_definitions,
                                        _BUSTED_json["background"]
                                       );
                                       
    (_BUSTED_json ["evidence ratios"])["constrained"] = busted.evidenceRatios ( (_BUSTED_json ["profiles"])["unconstrained"],  (_BUSTED_json ["profiles"])["constrained"]);
    (_BUSTED_json ["evidence ratios"])["optimized null"] = busted.evidenceRatios ( (_BUSTED_json ["profiles"])["unconstrained"],  (_BUSTED_json ["profiles"])["optimized null"]);
}

busted.taskTimerStop (2);

(_BUSTED_json ["timers"])["overall"] = _BUSTED_timers[2];
(_BUSTED_json ["timers"])["unconstrained"] = _BUSTED_timers[0];
(_BUSTED_json ["timers"])["constrained"] = _BUSTED_timers[1];

USE_JSON_FOR_MATRIX = 1;
fprintf (codon_data_info["json"], CLEAR_FILE, _BUSTED_json);
USE_JSON_FOR_MATRIX = 0;

//------------------------------------------------------------------------------ 
// HELPER FUNCTIONS FROM THIS POINT ON
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
function busted.getRateDistribution (model_description, key) {
  busted.getRateInformation.rate_classes = Abs ((model_description[key])["omegas"]);
  busted.getRateInformation.omega_info = {busted.getRateInformation.rate_classes,2};
  
  for (busted.getRateInformation.k = 0; busted.getRateInformation.k < busted.getRateInformation.rate_classes; busted.getRateInformation.k += 1) {
    busted.getRateInformation.omega_info[busted.getRateInformation.k][0] = Eval (((model_description[key])["omegas"])[busted.getRateInformation.k]);
    busted.getRateInformation.omega_info[busted.getRateInformation.k][1] = Eval (((model_description[key])["weights"])[busted.getRateInformation.k]);
  }
  return busted.getRateInformation.omega_info;
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
  ExecuteCommands (((model_parameters["FG"])["omegas"])[2] + ":=1");           
}

//------------------------------------------------------------------------------ 
function busted.computeSiteLikelihoods (id) {
   ConstructCategoryMatrix (_siteLike, *id, SITE_LOG_LIKELIHOODS);
   return _siteLike;         
}



//------------------------------------------------------------------------------ 
function busted.runLRT (ha, h0) {
    return {"LR" : 2*(ha-h0),
            "p" : 1-CChi2 (2*(ha-h0),2)};
}


//------------------------------------------------------------------------------ 
function busted._aux.define_bsrel_model (id,Q,weights,freqs) {
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
function busted._aux.define_parameter (key, value) {
   ExecuteCommands ("global `value` :< 1;");
   ExecuteCommands ("`value` :> 0;");
   ExecuteCommands ("`value` = " + busted.init_values[key]);
}

//------------------------------------------------------------------------------ 
function busted.io.evaluate_branch_lengths (model_parameters, tree_id, fg_set) {
    busted.io.evaluate_branch_lengths.res    = {};
    busted.io.evaluate_branch_lengths.bnames = BranchName (*tree_id, -1);
    for (busted.io.evaluate_branch_lengths.k = 0; 
         busted.io.evaluate_branch_lengths.k < Columns (busted.io.evaluate_branch_lengths.bnames)-1;
         busted.io.evaluate_branch_lengths.k += 1) {
         
         busted.io.evaluate_branch_lengths.lexpr = "";
         
         busted.io.evaluate_branch_lengths.b_name = busted.io.evaluate_branch_lengths.bnames[busted.io.evaluate_branch_lengths.k];
         if (fg_set [busted.io.evaluate_branch_lengths.b_name]) {
            //fprintf (stdout, busted.io.evaluate_branch_lengths.b_name, "=>FG\n");
            busted.io.evaluate_branch_lengths.lexpr = (model_parameters["FG"])["length"];
         } else {
            //fprintf (stdout, busted.io.evaluate_branch_lengths.b_name, "=>BG\n");
            busted.io.evaluate_branch_lengths.lexpr = (model_parameters["BG"])["length"];
         }
         Eval (model_parameters["length parameter"] + " = `tree_id`.`busted.io.evaluate_branch_lengths.b_name`." + model_parameters["length parameter"]);
         busted.io.evaluate_branch_lengths.res [ busted.io.evaluate_branch_lengths.b_name ] =
            Eval (busted.io.evaluate_branch_lengths.lexpr);
    }
    return busted.io.evaluate_branch_lengths.res;
}

//------------------------------------------------------------------------------ 
function busted.io.define_bsrel_models (foreground_id, background_id, frequencies) {

    model_parameters = 
        {"FG": {"omegas" : {}, "weights" : {}, "f" : {}, "Q" : {}, "length" : ""},
         "BG": {"omegas" : {}, "weights" : {}, "f" : {}, "Q" : {}, "length" : ""},
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
    
    ((model_parameters["FG"])["weights"])  = parameters.helper.stick_breaking (((model_parameters["FG"])["f"]), None);
    ((model_parameters["BG"])["weights"])  = parameters.helper.stick_breaking (((model_parameters["BG"])["f"]), None);
     
    busted.init_values = {"0" : 0.1, "1" : 0.5, "2" : 1};

    ((model_parameters["FG"])["omegas"])["busted._aux.define_parameter"][""];
    ((model_parameters["BG"])["omegas"])["busted._aux.define_parameter"][""];
    
    Eval (((model_parameters["FG"])["omegas"])[2] + ":<1e26");
    Eval (((model_parameters["FG"])["omegas"])[2] + ":>1");
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
    
    (model_parameters["BG"])["model"] = "`background_id`_model";
    (model_parameters["BG"])["length"] = busted._aux.define_bsrel_model ("`background_id`_model", (model_parameters["BG"])["Q"], (model_parameters["BG"])["weights"], frequencies["codon"]);
    (model_parameters["FG"])["model"] = "`foreground_id`_model";
    (model_parameters["FG"])["length"] = busted._aux.define_bsrel_model ("`foreground_id`_model", (model_parameters["FG"])["Q"], (model_parameters["FG"])["weights"], frequencies["codon"]);
    
    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("A","C")] = {"ID" : "AC"};
    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("A","T")] = {"ID" : "AT"};
    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("C","G")] = {"ID" : "CG"};
    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("C","T")] = {"ID" : "CT"};
    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("G","T")] = {"ID" : "GT"};
    
    model_parameters["length parameter"] = "t";
    
    return model_parameters;
}

//------------------------------------------------------------------------------ 

function relax._aux.io.countBranchSets (key, value) {
    available_models[value] += 1;
    return None;
}

function relax._aux.io.mapBranchSets (key, value) {
    (return_set[branch_set[value]])[key] = 1;
    return None;
}

//------------------------------------------------------------------------------ 
function relax.io.defineBranchSets (tree_definition) {
    
    available_models        = {};
    branch_set              = {};
    
    
    for (k = 0; k < Columns (tree_definition["model_list"]); k += 1) {
        available_models  [(tree_definition["model_list"])[k]] = 0;
    }
    
    (tree_definition["model_map"])["relax._aux.io.countBranchSets"][""];
    
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
    
    if (fgSet > 0) {
        branch_set [selectTheseForTesting[fgSet][0]] = "Test";
        return_set ["Test"] = {};
        if (option_count > 2) {
            ChoiceList  (bgSet,"Choose the set of reference branches (R set)",1,fgSet,selectTheseForTesting);    
            if (bgSet < 0) {
                return {};
            }    
            for (k = 0; k < option_count; k+=1) {
                if (k != bgSet && k != fgSet) {
                    branch_set [""] = "Ignore";
                    return_set ["Ignore"] = {};
                    break;
                }
            }
        }
        else {
            bgSet = 1-fgSet;
        }
        branch_set [selectTheseForTesting[bgSet][0]] = "Reference";
        return_set ["Reference"] = {};
    }
    
    (tree_definition["model_map"])["relax._aux.io.mapBranchSets"][""];
    
    return return_set;
    
}

//------------------------------------------------------------------------------------------------------------------------

function busted.evidenceRatios (ha, h0) {
    sites = Rows (ha) * Columns (ha);
    return ha["Exp(_MATRIX_ELEMENT_VALUE_-h0[_MATRIX_ELEMENT_COLUMN_])"];
}

//------------------------------------------------------------------------------------------------------------------------

function relax.taskTimerStart (index) {
    _RELAX_timers[index] = Time(1);
}

function relax.taskTimerStop (index) {
    _RELAX_timers[index] = Time(1) - _RELAX_timers[index];

}

//------------------------------------------------------------------------------------------------------------------------

function busted.getIC (logl,params,samples) {
    return -2*logl + 2*samples/(samples-params-1)*params;
}



//------------------------------------------------------------------------------------------------------------------------

lfunction busted.json_store_lf (json, name, ll, df, aicc, time, tree_length, tree_string, defs, has_bg) {
    (json["fits"])[name] = {"log-likelihood": ll,
                            "parameters": df,
                            "AIC-c" : aicc,
                            "runtime" : time,
                            "tree length" : tree_length,
                            "tree string" : tree_string,
                            "rate distributions" : {}};
                            
    (((json["fits"])[name])["rate distributions"])["FG"] = busted.getRateDistribution (defs, "FG");
    if (has_bg) {
        (((json["fits"])[name])["rate distributions"])["BG"] = busted.getRateDistribution (defs, "BG");
    }
}
