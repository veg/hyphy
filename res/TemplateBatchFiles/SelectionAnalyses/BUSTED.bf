RequireVersion ("2.31");


// namespace 'utility' for convenience functions
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


io.DisplayAnalysisBanner ({"info" : "BUSTED (branch-site unrestricted statistical test of episodic diversification)
                            uses a random effects branch-site model fitted jointly to all or a subset of tree branches
                            in order to test for alignment-wide evidence of episodic diversifying selection. Assuming
                            there is evidence of positive selection (i.e. there is an omega > 1), BUSTED will also perform
                            a quick evidence-ratio style analysis to explore which individual sites may have been subject to selection.",
                           "version" : "1.2",
                           "reference" : "*Gene-wide identification of episodic selection*, Mol Biol Evol. 32(5):1365-71",
                           "authors" : "Sergei L Kosakovsky Pond",
                           "contact" : "spond@temple.edu",
                           "requirements" : "in-frame codon alignment and a phylogenetic tree (optionally annotated with {})"
                          } );



utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

busted.timers  = {3,1};
busted.taskTimerStart (2);

busted.json    = { terms.json.fits: {},
                   terms.json.timers: {},
                   "profiles" : {},
                   "evidence ratios": {}
                  };

namespace busted {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ("busted");
}

io.ReportProgressMessageMD('BUSTED',  'selector', 'Branches to test for selection in the BUSTED analysis');

utility.ForEachPair (busted.selected_branches, "_partition_", "_selection_",
    "_selection_ = utility.Filter (_selection_, '_value_', '_value_ == terms.json.attribute.test');
     io.ReportProgressMessageMD('BUSTED',  'selector', 'Selected ' + Abs(_selection_) + ' branches to test in the BUSTED analysis: \\\`' + Join (', ',utility.Keys(_selection_)) + '\\\`')");

// check to see if there are any background branches

busted.has_background = FALSE;

namespace busted {
    doGTR ("busted");
}

busted.scaler_prefix = "busted.scaler";

estimators.fixSubsetOfEstimates(busted.gtr_results, busted.gtr_results["global"]);

namespace busted {
    doPartitionedMG ("busted", FALSE);
}

io.ReportProgressMessageMD ("BUSTED", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");

busted.final_partitioned_mg_results = estimators.FitMGREV (busted.filter_names, busted.trees, busted.codon_data_info ["code"], {
    "model-type": terms.global,
    "partitioned-omega": busted.selected_branches,
    "retain-lf-object": TRUE
}, busted.partitioned_mg_results);


io.ReportProgressMessageMD("BUSTED", "codon-refit", "* Log(L) = " + Format(busted.final_partitioned_mg_results["LogL"],8,2));
busted.global_dnds = selection.io.extract_global_MLE_re (busted.final_partitioned_mg_results, "^" + terms.omega_ratio);
utility.ForEach (busted.global_dnds, "_value_", 'io.ReportProgressMessageMD ("BUSTED", "codon-refit", "* " + _value_["description"] + " = " + Format (_value_["MLE"],8,4));');

busted.test.bsrel_model =  model.generic.DefineMixtureModel("models.codon.BS_REL.ModelDescription",
        "busted.test", {
            "0": parameters.Quote(terms.global),
            "1": busted.codon_data_info["code"],
            "2": parameters.Quote (3) // the number of rate classes
        },
        busted.filter_names,
        None);

busted.background.bsrel_model =  model.generic.DefineMixtureModel("models.codon.BS_REL.ModelDescription",
        "busted.background", {
            "0": parameters.Quote(terms.global),
            "1": busted.codon_data_info["code"],
            "2": parameters.Quote (3) // the number of rate classes
        },
        busted.filter_names,
        None);

// set up parameter constraints

for (busted.i = 1; busted.i < 4; busted.i += 1) {
    parameters.SetRange (model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (terms.omega_ratio,busted.i)), terms.range01);
    parameters.SetRange (model.generic.GetGlobalParameter (busted.background.bsrel_model , terms.AddCategory (terms.omega_ratio,busted.i)), terms.range01);
}

busted.test.omega3  = model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (terms.omega_ratio,3));
parameters.SetRange (busted.test.omega3, terms.range_gte1);

busted.model_map = {};

for (busted.partition_index = 0; busted.partition_index < busted.partition_count; busted.partition_index += 1) {
    busted.model_map + { "busted.test" : utility.Filter (busted.selected_branches[busted.partition_index], '_value_', '_value_ == terms.json.attribute.test'),
					               "busted.background" : utility.Filter (busted.selected_branches[busted.partition_index], '_value_', '_value_ != terms.json.attribute.test')};

}

busted.full_model =  estimators.FitLF (busted.filter_names, busted.trees, busted.model_map, busted.final_partitioned_mg_results);

return 0;

codon_data_info = alignments.PromptForGeneticCodeAndAlignment ("codon_data", "codon_filter");
codon_data_info["json"] = codon_data_info["file"] + ".BUSTED.json";
name_mapping = codon_data_info [utility.getGlobalValue("terms.json.name_mapping")];

if (None == name_mapping) { /** create a 1-1 mapping if nothing was done */
    name_mapping = {};
    utility.ForEach (alignments.GetSequenceNames ("codon_data"), "_value_", "`&name_mapping`[_value_] = _value_");
}

io.ReportProgressMessageMD ("BUSTED", "Data", "Loaded an MSA with " + codon_data_info["sequences"] +
                            " sequences and " + codon_data_info["sites"] + " codons from '" + codon_data_info["file"] + "'");

codon_lists = models.codon.MapCode (codon_data_info["code"]);

_Genetic_Code = codon_data_info["code"];



codon_frequencies     = frequencies._aux.CF3x4(frequencies._aux.empirical.collect_data ("codon_filter",3,1,1),
                        utility.getGlobalValue ("models.DNA.alphabet"), codon_lists["sense"], codon_lists["stop"]);


partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (codon_data_info[utility.getGlobalValue("terms.json.partitions")], name_mapping);

io.CheckAssertion("utility.Array1D (`&partitions_and_trees`) == 1", "BUSTED only works on a single partition dataset");

tree_definition   = utility.Map (partitions_and_trees, "_partition_", '_partition_["tree"]');


busted.selected_branches = busted.io.selectBranchesToTest (tree_definition[0]);
busted.json ["test set"] = Join (",",Rows(busted.selected_branches));

io.ReportProgressMessageMD ("BUSTED", "Data", "Selected " + Abs (busted.selected_branches) + " branches as the test (foreground) set: " + Join (",", Rows (busted.selected_branches)));

busted.model_definitions = busted.io.define_bsrel_models  ("FG","BG", codon_frequencies);
io.ReportProgressMessageMD ("BUSTED", "GTR", "Obtaining initial branch lengths under the GTR model");
busted.gtr_results = estimators.FitGTR     ("codon_filter", tree_definition[0], None);
io.ReportProgressMessageMD ("BUSTED", "GTR", "Log(L) = " + busted.gtr_results["LogL"]);

model.ApplyModelToTree ("busted.tree", tree_definition[0], "", {"DEFAULT" : (busted.model_definitions["BG"])["model"],
                                                             (busted.model_definitions["FG"])["model"] : Rows (busted.selected_branches)});



busted.taskTimerStart (0);
LikelihoodFunction busted.LF = (codon_filter, busted.tree);

busted.json["background"] =  busted.hasBackground  ("busted.tree");

global busted.T_scaler = 4;
BUSTED.proportional_constraint = "busted.T_scaler";

BUSTED.model_specification = {};
BUSTED.model_specification[(busted.model_definitions["FG"])["model"]] = busted.model_definitions;
BUSTED.model_specification[(busted.model_definitions["BG"])["model"]] = busted.model_definitions;

estimators.ApplyExistingEstimates ("busted.LF", BUSTED.model_specification, busted.gtr_results, None);

io.ReportProgressMessageMD ("BUSTED", "BUSTED-alt", "Fitting the unconstrained branch-site model");

utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);
utility.SetEnvVariable ("OPTIMIZATION_PRECISION", 0.1);
utility.SetEnvVariable ("ASSUME_REVERSIBLE_MODELS", TRUE);

busted.bls = busted.io.evaluate_branch_lengths (busted.model_definitions, "busted.tree", busted.selected_branches);

Optimize (busted.MLE_HA, busted.LF);

parameters.UnconstrainParameterSet ("busted.LF", {{terms.lf.local.constrained}});

utility.SetEnvVariable ("USE_LAST_RESULTS", 0.001);
Optimize (busted.MLE_HA, busted.LF);
io.SpoolLF ("busted.LF", codon_data_info["file"], None);
busted_positive_class = busted.checkForPS (busted.model_definitions);
io.ReportProgressMessageMD ("BUSTED", "BUSTED-alt", "Log(L) = " + busted.MLE_HA[1][0] + ". Unrestricted class omega = " + busted_positive_class["omega"] + " (weight = " + busted_positive_class["weight"] + ")");


busted.sample_size             =codon_data_info["sites"] * codon_data_info["sequences"];
busted.taskTimerStop (0);

busted.bls = busted.io.evaluate_branch_lengths (busted.model_definitions, "busted.tree", busted.selected_branches);
busted.tavl         = busted.tree ^ 0;
busted.renderString = PostOrderAVL2StringDistances (busted.tavl, busted.bls);
UseModel (USE_NO_MODEL);
Tree busted.T = busted.renderString;

busted.json_store_lf                (busted.json, "Unconstrained model",
        busted.MLE_HA[1][0], busted.MLE_HA[1][1]+9,
        busted.getIC (busted.MLE_HA[1][0], busted.MLE_HA[1][1]+9, busted.sample_size) ,
        busted.timers[0],
        +BranchLength (busted.T,-1),
        Format (busted.T, 1,1),
        busted.model_definitions,
        busted.json["background"]
        );


busted.profiles = {};
(busted.json ["profiles"])["unconstrained"] = busted.computeSiteLikelihoods ("busted.LF");


if (busted_positive_class["omega"] < 1 || busted_positive_class["weight"] < 1e-8) {
    io.ReportProgressMessageMD ("BUSTED", "Results", "No evidence for positive selection under the unconstrained model, skipping constrained model fitting");
    busted.json ["test results"] = busted.runLRT (0, 0);
} else {
    busted.taskTimerStart (1);

    io.ReportProgressMessageMD ("BUSTED",  "BUSTED-null", "Fitting the branch-site model that disallows omega > 1 among foreground branches");
    busted.constrainTheModel (busted.model_definitions);
    (busted.json ["profiles"])["constrained"] = busted.computeSiteLikelihoods ("busted.LF");;
    Optimize (busted.MLE_H0, busted.LF);
    io.SpoolLF ("busted.LF", codon_data_info["file"], "null");
    (busted.json ["profiles"])["optimized null"] = busted.computeSiteLikelihoods ("busted.LF");;
    io.ReportProgressMessageMD ("BUSTED", "BUSTED-null", "Log(L) = " + busted.MLE_H0[1][0]);
    busted.LRT = busted.runLRT (busted.MLE_HA[1][0], busted.MLE_H0[1][0]);

    busted.json ["test results"] = busted.LRT;

    io.ReportProgressMessageMD ("BUSTED", "BUSTED-results", "Likelihood ratio test for episodic positive selection, p = " + busted.LRT["p"]);
     busted.taskTimerStop (1);

    busted.bls = busted.io.evaluate_branch_lengths (busted.model_definitions, "busted.tree", busted.selected_branches);
    busted.tavl         = busted.tree ^ 0;
    busted.renderString = PostOrderAVL2StringDistances (busted.tavl, busted.bls);
    UseModel (USE_NO_MODEL);
    Tree busted.T = busted.renderString;

    busted.json_store_lf                (busted.json,
                                        "Constrained model", busted.MLE_H0[1][0],
                                        busted.MLE_H0[1][1]+9,
                                        busted.getIC (busted.MLE_H0[1][0], busted.MLE_H0[1][1]+9, busted.sample_size) ,
                                        busted.timers[1],
                                         +BranchLength (busted.T,-1),
                                        Format (busted.T, 1,1),
                                        busted.model_definitions,
                                        busted.json["background"]
                                       );

    (busted.json ["evidence ratios"])["constrained"] = busted.evidenceRatios ( (busted.json ["profiles"])["unconstrained"],  (busted.json ["profiles"])["constrained"]);
    (busted.json ["evidence ratios"])["optimized null"] = busted.evidenceRatios ( (busted.json ["profiles"])["unconstrained"],  (busted.json ["profiles"])["optimized null"]);
}

busted.taskTimerStop (2);

(busted.json ["timers"])["overall"] = busted.timers[2];
(busted.json ["timers"])["unconstrained"] = busted.timers[0];
(busted.json ["timers"])["constrained"] = busted.timers[1];

USE_JSON_FOR_MATRIX = 1;
fprintf (codon_data_info["json"], CLEAR_FILE, busted.json);
USE_JSON_FOR_MATRIX = 0;

return busted.json;

//------------------------------------------------------------------------------
// HELPER FUNCTIONS FROM HTHIS POINT ON
//------------------------------------------------------------------------------


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
         "parameters" : {"global" : {}, "local" : {}}
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
    Eval (((model_parameters["FG"])["omegas"])[2] + ":>1");
    Eval (((model_parameters["BG"])["omegas"])[2] + ":<1e26");

    busted.init_values = {"0" : 0.8, "1" : 0.75};

    ((model_parameters["FG"])["f"])["busted._aux.define_parameter"][""];
    ((model_parameters["BG"])["f"])["busted._aux.define_parameter"][""];

    busted.nuc = {4,3};
    for (k = 0; k < 4; k+=1) {
        for (k2 = 0; k2 < 3; k2 += 1) {
            busted.nuc [k][k2] = ((frequencies["bases"])[utility.getGlobalValue ("models.DNA.alphabet")[k]])[k2];
        }
    }

    for (k = 1; k <= 3; k +=1) {

        ((model_parameters["FG"])["Q"]) + ("Q_`foreground_id`_" + k);
        PopulateModelMatrix			  (((model_parameters["FG"])["Q"])[k-1],  busted.nuc, "t",((model_parameters["FG"])["omegas"])[k-1], "");

        ((model_parameters["BG"])["Q"]) + ("Q_`background_id`_" + k);
        PopulateModelMatrix			  (((model_parameters["BG"])["Q"])[k-1],  busted.nuc, "t",((model_parameters["BG"])["omegas"])[k-1], "");
    }


    (model_parameters["BG"])["model"] = "`background_id`_model";
    (model_parameters["BG"])["length"] = busted._aux.define_bsrel_model ("`background_id`_model", (model_parameters["BG"])["Q"], (model_parameters["BG"])["weights"], frequencies["codons"]);
    (model_parameters["FG"])["model"] = "`foreground_id`_model";
    (model_parameters["FG"])["length"] = busted._aux.define_bsrel_model ("`foreground_id`_model", (model_parameters["FG"])["Q"], (model_parameters["FG"])["weights"], frequencies["codons"]);

    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("A","C")] = "AC";
    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("A","T")] = "AT";
    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("C","G")] = "CG";
    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("C","T")] = "CT";
    ((model_parameters["parameters"])["global"])[terms.nucleotideRate ("G","T")] = "GT";
     model_parameters["set-branch-length"] = "busted.aux.copy_branch_length";

    ((model_parameters["parameters"])[terms.local])[terms.timeParameter ()] = "t";

     model_parameters["length parameter"] = "t";

    return model_parameters;
}

function busted.aux.copy_branch_length (model, value, parameter) {

    busted.aux.copy_branch_length.t = ((model["parameters"])["local"])[terms.timeParameter ()];

    if (Abs (BUSTED.proportional_constraint)) {
        Eval ("`parameter`.`busted.aux.copy_branch_length.t` := `BUSTED.proportional_constraint` * " + value);
    } else {
        Eval ("`parameter`.`busted.aux.copy_branch_length.t` = " + value);
    }



    if (Type (relax.aux.copy_branch_length.k) == "String") {
        Eval ("`parameter`.`relax.aux.copy_branch_length.k` = 1");
    }
}

//------------------------------------------------------------------------------
lfunction busted.io.selectBranchesToTest (tree_definition) {

    extra_models = {};

    for (k = 0; k < Columns (tree_definition["model_list"]); k += 1) {
        model_id = (tree_definition["model_list"])[k];
        if (model_id != "") {
            extra_models  + model_id;
        }
    }


    UseModel (USE_NO_MODEL);
    ExecuteCommands ("Topology  T = " + tree_definition["string"]);
    tAVL = T^0;

    totalBranchCount = Abs (tAVL) - 2;

    selectedBranches = {};
    selectTheseForTesting = {totalBranchCount + 3 + Abs (extra_models), 2};
    selectTheseForTesting [0][0] = "All";  selectTheseForTesting [0][1] = "Test for selection on all branches jointly";
    selectTheseForTesting [1][0] = "Internal";  selectTheseForTesting [1][1] = "Test for selection on all internal branches jointly";
    selectTheseForTesting [2][0] = "Leaves";  selectTheseForTesting [2][1] = "Test for selection on all terminal branches jointly";

    for (k = 0; k < Abs (extra_models); k+=1) {
        selectTheseForTesting [3+k][0] = "Set " + extra_models[k];
        selectTheseForTesting [3+k][1] = "Test for selection on all branches labeled with {" + extra_models[k] + "} jointly";
    }

    for (k = 0; k < totalBranchCount; k += 1) {
        selectTheseForTesting [k+3 + Abs (extra_models)][0] = (tAVL[k+1])["Name"];
        selectTheseForTesting [k+3 + Abs (extra_models)][1] = "Add branch '" + selectTheseForTesting [k+3 + Abs (extra_models)][0] + "' to the test set";

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
        if (whichBranchesToTest [k] < 3 + Abs (extra_models)) {
            model_key = extra_models [whichBranchesToTest [k] - 3];
            for (k2 = 1; k2 <=  totalBranchCount; k2 += 1) {
                bName = (tAVL[k2])["Name"];
                if ((tree_definition["model_map"])[bName] == model_key) {
                    selectedBranches [bName] = 1;
                }
            }
            return selectedBranches;
        }

        selectedBranches [(tAVL[whichBranchesToTest [k] - 2 - Abs (extra_models)])["Name"]] = 1;
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
    busted.timers[index] = Time(1);
}

function busted.taskTimerStop (index) {
    busted.timers[index] = Time(1) - busted.timers[index];

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

