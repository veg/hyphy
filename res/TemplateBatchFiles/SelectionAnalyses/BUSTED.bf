RequireVersion ("2.31");



LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("CF3x4");
LoadFunctionLibrary("TreeTools");


// namespace 'utility' for convenience functions
LoadFunctionLibrary("libv3/UtilityFunctions.bf");

// namespace 'io' for interactive/datamonkey i/o functions
LoadFunctionLibrary("libv3/IOFunctions.bf");

LoadFunctionLibrary("libv3/tasks/estimators.bf");

LoadFunctionLibrary("libv3/tasks/alignments.bf");

LoadFunctionLibrary("libv3/models/codon.bf");

LoadFunctionLibrary("libv3/tasks/trees.bf");

LoadFunctionLibrary("libv3/tasks/genetic_code.bf");

//LoadFunctionLibrary("libv3/models/terms.bf");


busted.analysis_description = {terms.io.info : "BUSTED (branch-site unrestricted statistical test of episodic diversification)
                            uses a random effects branch-site model fitted jointly to all or a subset of tree branches
                            in order to test for alignment-wide evidence of episodic diversifying selection. Assuming
                            there is evidence of positive selection (i.e. there is an omega > 1), BUSTED will also perform
                            a quick evidence-ratio style analysis to explore which individual sites may have been subject to selection.",
                           terms.io.version : "1.2",
                           terms.io.reference : "*Gene-wide identification of episodic selection*, Mol Biol Evol. 32(5):1365-71",
                           terms.io.authors : "Sergei L Kosakovsky Pond",
                           terms.io.contact : "spond@temple.edu",
                           terms.io.requirements : "in-frame codon alignment and a phylogenetic tree (optionally annotated with {})"
                          };
io.DisplayAnalysisBanner (busted.analysis_description);



/*------------------------------------------------------------------------------
    BranchSiteTemplate Defines

    BuildCodonFrequencies (obsF);
    PopulateModelMatrix (ModelMatrixName&, EFV, synrateP, globalP, nonsynRateP);

------------------------------------------------------------------------------*/


skipCodeSelectionStep = TRUE;
LoadFunctionLibrary ("chooseGeneticCode"); 

LoadFunctionLibrary("BranchSiteTemplate");

busted.FG = "FG";
busted.BG = "BG";
busted.BG_FG_map = {0: "background", 1: "foreground"}; // for partitions output
busted.background = "background";
busted.unconstrained = "unconstrained";
busted.constrained = "constrained";
busted.optimized_null = "optimized null";

busted.nrates = 3;



busted.timers  = {3,1};
busted.timers_d = {"Overall":2, 
                  "Unconstrained model":0,
                  "Constrained model":1};

//busted.taskTimerStart (2);
busted.taskTimerStart (busted.timers_d["Overall"]);



busted.json    = { terms.json.analysis: busted.analysis_description,
                   terms.json.input: {},
                   terms.json.partitions: {}, 
                   terms.json.attribute.background: {},
                   terms.json.fits : {},
                   terms.json.timers : {},
                   terms.json.site_log_likelihood : {}, 
                   terms.json.evidence_ratios: {}
                  };


codon_data_info = alignments.PromptForGeneticCodeAndAlignment ("codon_data", "codon_filter");
codon_data_info[terms.json.json] = codon_data_info[terms.file] + ".busted.json";
busted.name_mapping = codon_data_info [utility.getGlobalValue("terms.data.name_mapping")];




if (None == busted.name_mapping) { /** create a 1-1 mapping if nothing was done */
    busted.name_mapping = {};
    utility.ForEach (alignments.GetSequenceNames ("codon_data"), "_value_", "`&busted.name_mapping`[_value_] = _value_");
}


io.ReportProgressMessageMD ("BUSTED", "Data", "Loaded an MSA with " + codon_data_info[terms.data.sequences] +
                            " sequences and " + codon_data_info[terms.data.sites] + " codons from '" + codon_data_info[terms.data.file] + "'");

codon_lists = models.codon.MapCode (codon_data_info[terms.code]);

_Genetic_Code = codon_data_info[terms.code];


codon_frequencies     = frequencies._aux.CF3x4(frequencies._aux.empirical.collect_data ("codon_filter",3,1,1),
                        utility.getGlobalValue ("models.DNA.alphabet"), codon_lists[terms.sense_codons], codon_lists[terms.stop_codons]);


partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (codon_data_info[utility.getGlobalValue("terms.data.partitions")], busted.name_mapping);

io.CheckAssertion("utility.Array1D (`&partitions_and_trees`) == 1", "BUSTED only works on a single partition dataset");

tree_definition   = utility.Map (partitions_and_trees, "_partition_", '_partition_[terms.data.tree]');

busted.selected_branches = busted.io.selectBranchesToTest (tree_definition[0]);


io.ReportProgressMessageMD ("BUSTED", "Data", "Selected " + Abs (busted.selected_branches) + " branches as the test (foreground) set: " + Join (",", Rows (busted.selected_branches)));


busted.model_definitions = busted.io.define_bsrel_models  (busted.FG, busted.BG, codon_frequencies);
io.ReportProgressMessageMD ("BUSTED", "GTR", "Obtaining initial branch lengths under the GTR model");
busted.gtr_results = estimators.FitGTR     ("codon_filter", tree_definition[0], None);
io.ReportProgressMessageMD ("BUSTED", "GTR", "Log(L) = " + busted.gtr_results[terms.fit.log_likelihood]);

model.ApplyModelToTree ("busted.tree", tree_definition[0], "", {"DEFAULT" : (busted.model_definitions[busted.BG])[terms.model],
                                                             (busted.model_definitions[busted.FG])[terms.model] : Rows (busted.selected_branches)});




//busted.taskTimerStart (0);
busted.taskTimerStart(busted.timers_d["Unconstrained model"]);

LikelihoodFunction busted.LF = (codon_filter, busted.tree);

busted.json[terms.json.attribute.background] =  busted.hasBackground  ("busted.tree");

global busted.T_scaler = 4;
busted.proportional_constraint = "busted.T_scaler";

busted.model_specification = {};
busted.model_specification[(busted.model_definitions[busted.FG])[terms.model]] = busted.model_definitions;
busted.model_specification[(busted.model_definitions[busted.BG])[terms.model]] = busted.model_definitions;

estimators.ApplyExistingEstimates ("busted.LF", busted.model_specification, busted.gtr_results, None);

io.ReportProgressMessageMD ("BUSTED", "BUSTED-alt", "Fitting the unconstrained branch-site model");

USE_LAST_RESULTS = 1;
OPTIMIZATION_PRECISION = 0.1;
ASSUME_REVERSIBLE_MODELS = 1;

busted.bls = busted.io.evaluate_branch_lengths (busted.model_definitions, "busted.tree", busted.selected_branches);

Optimize (busted.MLE_HA, busted.LF);

parameters.UnconstrainParameterSet ("busted.LF", {{terms.lf.local.constrained}});

utility.SetEnvVariable ("OPTIMIZATION_PRECISION", 0.001);
Optimize (busted.MLE_HA, busted.LF);
// Uncomment the line below to export the LF fit file, a hyphy nexus file.
//io.SpoolLF ("busted.LF", codon_data_info[terms.data.file], None);
busted_positive_class = busted.checkForPS (busted.model_definitions);
io.ReportProgressMessageMD ("BUSTED", "BUSTED-alt", "Log(L) = " + busted.MLE_HA[1][0] + ". Unrestricted class omega = " + busted_positive_class[terms.omega] + " (weight = " + busted_positive_class[terms.weight] + ")");


busted.sample_size             =codon_data_info[terms.data.sites] * codon_data_info[terms.data.sequences];

//busted.taskTimerStop (0);
busted.taskTimerStop(busted.timers_d["Unconstrained model"]);


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
        busted.json[busted.background]
        );


busted.siteLogL = {};
(busted.json [terms.json.site_log_likelihood])[busted.unconstrained] = busted.computeSiteLikelihoods ("busted.LF");


if (busted_positive_class[terms.omega] < terms.range_almost_01[terms.upper_bound] || busted_positive_class[terms.weight] < terms.range_almost_01[terms.lower_bound]) {
    io.ReportProgressMessageMD ("BUSTED", "Results", "No evidence for positive selection under the unconstrained model, skipping constrained model fitting");
    busted.json [terms.json.test_results] = busted.runLRT (0, 0);
} else {

//    busted.taskTimerStart (1);
    busted.taskTimerStart(busted.timers_d["Constrained model"]);

    io.ReportProgressMessageMD ("BUSTED",  "BUSTED-null", "Fitting the branch-site model that disallows omega > 1 among foreground branches");
    busted.constrainTheModel (busted.model_definitions);
    (busted.json [terms.json.site_log_likelihood])[busted.constrained] = busted.computeSiteLikelihoods ("busted.LF");;
    Optimize (busted.MLE_H0, busted.LF);
    // Uncomment line below to get the LF file.
    //io.SpoolLF ("busted.LF", codon_data_info[terms.data.file], "null");
    (busted.json [terms.json.site_log_likelihood])[busted.optimized_null] = busted.computeSiteLikelihoods ("busted.LF");;
    io.ReportProgressMessageMD ("BUSTED", "BUSTED-null", "Log(L) = " + busted.MLE_H0[1][0]);
    busted.LRT = busted.runLRT (busted.MLE_HA[1][0], busted.MLE_H0[1][0]);

    busted.json [terms.json.test_results] = busted.LRT;

    io.ReportProgressMessageMD ("BUSTED", "BUSTED-results", "Likelihood ratio test for episodic positive selection, p = " + busted.LRT[terms.p_value]);

//    busted.taskTimerStop (1);
    busted.taskTimerStop(busted.timers_d["Constrained model"]);


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
                                        busted.json[busted.background]
                                       );

    (busted.json [terms.json.evidence_ratios])[busted.constrained] = busted.evidenceRatios ( (busted.json [terms.json.site_log_likelihood ])[busted.unconstrained],  (busted.json [terms.json.site_log_likelihood ])[busted.constrained]);
    (busted.json [terms.json.evidence_ratios])[busted.optimized_null] = busted.evidenceRatios ( (busted.json [terms.json.site_log_likelihood ])[busted.unconstrained],  (busted.json [terms.json.site_log_likelihood ])[busted.optimized_null]);
}

//busted.taskTimerStop (2);
busted.taskTimerStop(busted.timers_d["Overall"]);




/* 
    Populate remaining bits of JSON
*/


/*
Add input information to the JSON
*/
(busted.json[terms.json.input])[terms.json.file] = codon_data_info[terms.data.file];
(busted.json[terms.json.input])[terms.json.sequences] = codon_data_info[terms.data.sequences];
(busted.json[terms.json.input])[terms.json.sites] = codon_data_info[terms.data.sites];
(busted.json[terms.json.input])[terms.json.tree_string] = (tree_definition[0])[terms.trees.newick_with_lengths];


/*
 Instead of JSON key test set (see next comment line), we now have a single dictionary of node:<background/foreground>.
former: busted.json ["test set"] = Join (",",Rows(busted.selected_branches));
*/
busted.partitions = {};
all_nodes = utility.Keys( (tree_definition[0])[terms.trees.model_map] );
utility.ForEach (all_nodes, "_value_", "`&busted.partitions`[_value_] = `&busted.BG_FG_map`[utility.Has(busted.selected_branches, _value_, None)]");
busted.json[terms.json.partitions] = busted.partitions; // Improved


(busted.json [terms.json.timers])["Overall"] = busted.timers[busted.timers_d["Overall"]];
(busted.json [terms.json.timers])["Unconstrained model"] = busted.timers[busted.timers_d["Unconstrained model"]];
(busted.json [terms.json.timers])["Constrained model"] = busted.timers[busted.timers_d["Constrained model"]];



USE_JSON_FOR_MATRIX = 1;
fprintf (codon_data_info[terms.json.json], CLEAR_FILE, busted.json);
USE_JSON_FOR_MATRIX = 0;







return busted.json;

//------------------------------------------------------------------------------
// HELPER FUNCTIONS FROM THIS POINT ON
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
function busted.hasBackground (id) {
   ExecuteCommands ("GetInformation (busted.nodeMap, `id`)");
   busted.nodeMap = Columns(busted.nodeMap);
   return Rows (busted.nodeMap) * Columns (busted.nodeMap) > 1;
}

//------------------------------------------------------------------------------
function busted.getRateDistribution (model_description, key) {

  busted.getRateInformation.rate_classes = Abs ((model_description[key])[terms.omegas]);
  busted.getRateInformation.omega_info = {busted.getRateInformation.rate_classes,2};

  for (busted.getRateInformation.k = 0; busted.getRateInformation.k < busted.getRateInformation.rate_classes; busted.getRateInformation.k += 1) {
    busted.getRateInformation.omega_info[busted.getRateInformation.k][0] = Eval (((model_description[key])[terms.omegas])[busted.getRateInformation.k]);
    busted.getRateInformation.omega_info[busted.getRateInformation.k][1] = Eval (((model_description[key])[terms.weights])[busted.getRateInformation.k]);
  }
  return busted.getRateInformation.omega_info;
}


//------------------------------------------------------------------------------
function busted._aux.free_lengths (key, value) {
    ExecuteCommands ("busted.tree." + key + ".t = busted.tree." + key + ".t");
}


//------------------------------------------------------------------------------
function busted.checkForPS (model_parameters) {
   return {terms.omega :Eval (((model_parameters[busted.FG])[terms.omegas])[2]),
           terms.weight: Eval (((model_parameters[busted.FG])[terms.weights])[2])};

}

//------------------------------------------------------------------------------
function busted.constrainTheModel (model_parameters) {
  ExecuteCommands (((model_parameters[busted.FG])[terms.omegas])[2] + ":=1");
}

//------------------------------------------------------------------------------
function busted.computeSiteLikelihoods (id) {
   ConstructCategoryMatrix (_siteLike, *id, SITE_LOG_LIKELIHOODS);
   return _siteLike;
}



//------------------------------------------------------------------------------
function busted.runLRT (ha, h0) {
    return {terms.LR : 2*(ha-h0),
            terms.p_value : 1-CChi2 (2*(ha-h0),2)};
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
            busted.io.evaluate_branch_lengths.lexpr = (model_parameters[busted.FG])[terms.model.length];
         } else {
            //fprintf (stdout, busted.io.evaluate_branch_lengths.b_name, "=>BG\n");
            busted.io.evaluate_branch_lengths.lexpr = (model_parameters[busted.BG])[terms.model.length];
         }
         Eval (model_parameters[terms.model.length_parameter] + " = `tree_id`.`busted.io.evaluate_branch_lengths.b_name`." + model_parameters[terms.model.length_parameter]);
         busted.io.evaluate_branch_lengths.res [ busted.io.evaluate_branch_lengths.b_name ] =
            Eval (busted.io.evaluate_branch_lengths.lexpr);
    }
    return busted.io.evaluate_branch_lengths.res;
}

//------------------------------------------------------------------------------
function busted.io.define_bsrel_models (foreground_id, background_id, frequencies) {

    model_parameters =
        {busted.FG: {terms.omegas : {}, terms.weights : {}, terms.freqs : {}, terms.model.rate_matrix : {}, terms.model.length : ""},
         busted.BG: {terms.omegas : {}, terms.weights : {}, terms.freqs : {}, terms.model.rate_matrix : {}, terms.model.length : ""},
         terms.parameters : {terms.global : {}, terms.local : {}}
          };

    for (k = 1; k <= busted.nrates; k +=1) {
        tag = "" + k;
        ((model_parameters[busted.FG])[terms.omegas]) + "`foreground_id`_omega_`tag`";
        ((model_parameters[busted.BG])[terms.omegas]) + "`background_id`_omega_`tag`";
        if (k < busted.nrates) {
            ((model_parameters[busted.FG])[terms.freqs]) + "`foreground_id`_f_`tag`";
            ((model_parameters[busted.BG])[terms.freqs]) + "`background_id`_f_`tag`";
        }

    }
  

    ((model_parameters[busted.FG])[terms.weights])  = busted._aux.stick_breaking (((model_parameters[busted.FG])[terms.freqs]));
    ((model_parameters[busted.BG])[terms.weights])  = busted._aux.stick_breaking (((model_parameters[busted.BG])[terms.freqs]));

    busted.init_values = {"0" : 0.1, "1" : 0.5, "2" : 1};

    ((model_parameters[busted.FG])[terms.omegas])["busted._aux.define_parameter"][""];
    ((model_parameters[busted.BG])[terms.omegas])["busted._aux.define_parameter"][""];

    Eval (((model_parameters[busted.FG])[terms.omegas])[2] + ":<" + terms.range_gte1[terms.upper_bound]);
    Eval (((model_parameters[busted.FG])[terms.omegas])[2] + ":>" + terms.range_gte1[terms.lower_bound]);
    Eval (((model_parameters[busted.BG])[terms.omegas])[2] + ":<" + terms.range_gte1[terms.upper_bound]);

    busted.init_values = {"0" : 0.8, "1" : 0.75};

    ((model_parameters[busted.FG])[terms.freqs])["busted._aux.define_parameter"][""];
    ((model_parameters[busted.BG])[terms.freqs])["busted._aux.define_parameter"][""];

    busted.nuc = {4,3};
    for (k = 0; k < 4; k+=1) {
        for (k2 = 0; k2 < 3; k2 += 1) {
            busted.nuc [k][k2] = ((frequencies[terms.bases])[utility.getGlobalValue ("models.DNA.alphabet")[k]])[k2];
        }
    }



    for (k = 1; k <= busted.nrates; k +=1) {

        ((model_parameters[busted.FG])[terms.model.rate_matrix]) + ("Q_`foreground_id`_" + k);
        PopulateModelMatrix			  (((model_parameters[busted.FG])[terms.model.rate_matrix])[k-1],  busted.nuc, terms.default_time,((model_parameters[busted.FG])[terms.omegas])[k-1], "");

        ((model_parameters[busted.BG])[terms.model.rate_matrix]) + ("Q_`background_id`_" + k);
        PopulateModelMatrix			  (((model_parameters[busted.BG])[terms.model.rate_matrix])[k-1],  busted.nuc, terms.default_time,((model_parameters[busted.BG])[terms.omegas])[k-1], "");
    }


    (model_parameters[busted.BG])[terms.model] = "`background_id`_model";
    (model_parameters[busted.BG])[terms.model.length] = busted._aux.define_bsrel_model ("`background_id`_model", (model_parameters[busted.BG])[terms.model.rate_matrix], (model_parameters[busted.BG])[terms.weights], frequencies[terms.codons]);
    (model_parameters[busted.FG])[terms.model] = "`foreground_id`_model";
    (model_parameters[busted.FG])[terms.model.length] = busted._aux.define_bsrel_model ("`foreground_id`_model", (model_parameters[busted.FG])[terms.model.rate_matrix], (model_parameters[busted.FG])[terms.weights], frequencies[terms.codons]);

    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("A","C")] = "AC";
    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("A","T")] = "AT";
    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("C","G")] = "CG";
    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("C","T")] = "CT";
    ((model_parameters[terms.parameters])[terms.global])[terms.nucleotideRate ("G","T")] = "GT";
     model_parameters["set-branch-length"] = "busted.aux.copy_branch_length";

    ((model_parameters[terms.parameters])[terms.local])[terms.timeParameter ()] = terms.default_time;

     model_parameters[terms.model.length_parameter] = terms.default_time;

    return model_parameters;
}

function busted.aux.copy_branch_length (model, value, parameter) {

    busted.aux.copy_branch_length.t = ((model[terms.parameters])[terms.local])[terms.timeParameter ()];

    if (Abs (busted.proportional_constraint)) {
        Eval ("`parameter`.`busted.aux.copy_branch_length.t` := `busted.proportional_constraint` * " + value);
    } else {
        Eval ("`parameter`.`busted.aux.copy_branch_length.t` = " + value);
    }



    // This can probably be deleted. Appears to be accidental carry-over from RELAX.bf 
    /*
    if (Type (relax.aux.copy_branch_length.k) == "String") {
        Eval ("`parameter`.`relax.aux.copy_branch_length.k` = 1");
    }
    */
}

//------------------------------------------------------------------------------
lfunction busted.io.selectBranchesToTest (tree_definition) {

    extra_models = {};

    for (k = 0; k < Columns (tree_definition[utility.getGlobalValue("terms.trees.model_list")]); k += 1) {
        model_id = (tree_definition[utility.getGlobalValue("terms.trees.model_list")])[k];
        if (model_id != "") {
            extra_models  + model_id;
        }
    }


    UseModel (USE_NO_MODEL);
    ExecuteCommands ("Topology  T = " + tree_definition[utility.getGlobalValue("terms.trees.newick")]);
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
                if ((tree_definition[terms.trees.model_map])[bName] == model_key) {
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

// TODO: Should this be a function (as it is in RELAX.bf) or lfunction as written here?
function busted.json_store_lf (json, name, ll, df, aicc, time, tree_length, tree_string, defs, has_bg) {
    
    (json[terms.json.fits])[name] = {terms.json.log_likelihood     : ll,
                                     terms.json.parameters         : df,
                                     terms.json.AICc               : aicc,
                                     terms.json.runtime            : time,
                                     terms.json.tree_length        : tree_length,
                                     terms.json.tree_string        : tree_string,
                                     terms.json.rate_distributions : {}
                                     };

    (((json[terms.json.fits])[name])[terms.json.rate_distributions])[busted.FG] = busted.getRateDistribution (defs, busted.FG);
    if (has_bg) {
        (((json[terms.json.fits])[name])[terms.json.rate_distributions])[busted.BG] = busted.getRateDistribution (defs, busted.BG);
    }
}