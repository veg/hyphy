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
busted.rate_classes = 3;

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

utility.ForEachPair (busted.selected_branches, "_partition_", "_selection_",
    "_selection_ = utility.Filter (_selection_, '_value_', '_value_ != terms.json.attribute.test');
     if (utility.Array1D (_selection_)) { busted.has_background = TRUE;} ");

namespace busted {
    doGTR ("busted");
}

busted.scaler_prefix = "busted.scaler";

estimators.fixSubsetOfEstimates(busted.gtr_results, busted.gtr_results["global"]);

namespace busted {
    doPartitionedMG ("busted", FALSE);
}

utility.SetEnvVariable ("VERBOSITY_LEVEL", 0);

io.ReportProgressMessageMD ("BUSTED", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");

busted.final_partitioned_mg_results = estimators.FitMGREV (busted.filter_names, busted.trees, busted.codon_data_info ["code"], {
    "model-type": terms.local,
    "partitioned-omega": busted.selected_branches,
}, busted.partitioned_mg_results);


io.ReportProgressMessageMD("BUSTED", "codon-refit", "* Log(L) = " + Format(busted.final_partitioned_mg_results["LogL"],8,2));
busted.global_dnds = selection.io.extract_global_MLE_re (busted.final_partitioned_mg_results, "^" + terms.omega_ratio);

utility.ForEach (busted.global_dnds, "_value_", 'io.ReportProgressMessageMD ("BUSTED", "codon-refit", "* " + _value_["description"] + " = " + Format (_value_["MLE"],8,4));');

busted.test.bsrel_model =  model.generic.DefineMixtureModel("models.codon.BS_REL.ModelDescription",
        "busted.test", {
            "0": parameters.Quote(terms.global),
            "1": busted.codon_data_info["code"],
            "2": parameters.Quote (busted.rate_classes) // the number of rate classes
        },
        busted.filter_names,
        None);


busted.background.bsrel_model =  model.generic.DefineMixtureModel("models.codon.BS_REL.ModelDescription",
        "busted.background", {
            "0": parameters.Quote(terms.global),
            "1": busted.codon_data_info["code"],
            "2": parameters.Quote (busted.rate_classes) // the number of rate classes
        },
        busted.filter_names,
        None);


models.BindGlobalParameters ({"0" : busted.test.bsrel_model, "1" : busted.background.bsrel_model}, terms.nucleotideRate("[ACGT]","[ACGT]"));

busted.test_guess = busted.DistributionGuess(utility.Map (selection.io.extract_global_MLE_re (busted.final_partitioned_mg_results, "^" + terms.omega_ratio + ".+test.+"), "_value_",
            "_value_[terms.MLE]"));


busted.distribution = models.codon.BS_REL.ExtractMixtureDistribution(busted.test.bsrel_model);
parameters.SetStickBreakingDistribution (busted.distribution, busted.test_guess);


if (busted.has_background) {
    busted.model_object_map = { "busted.background" : busted.background.bsrel_model,
                                "busted.test" :       busted.test.bsrel_model };
    busted.background_distribution = models.codon.BS_REL.ExtractMixtureDistribution(busted.background.bsrel_model);
    busted.background_guess = busted.DistributionGuess(utility.Map (selection.io.extract_global_MLE_re (busted.final_partitioned_mg_results, "^" + terms.omega_ratio + ".+background.+"), "_value_",
            "_value_[terms.MLE]"));
    parameters.SetStickBreakingDistribution (busted.background_distribution, busted.background_guess);
} else {
    busted.model_object_map = { "busted.test" :       busted.test.bsrel_model };
}

// set up parameter constraints

for (busted.i = 1; busted.i < busted.rate_classes; busted.i += 1) {
    parameters.SetRange (model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (terms.omega_ratio,busted.i)), terms.range01);
    parameters.SetRange (model.generic.GetGlobalParameter (busted.background.bsrel_model , terms.AddCategory (terms.omega_ratio,busted.i)), terms.range01);
}

parameters.SetRange (model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (terms.omega_ratio,busted.rate_classes)), terms.range_gte1);

busted.model_map = {};

for (busted.partition_index = 0; busted.partition_index < busted.partition_count; busted.partition_index += 1) {
    busted.model_map + { "busted.test" : utility.Filter (busted.selected_branches[busted.partition_index], '_value_', '_value_ == terms.json.attribute.test'),
					     "busted.background" : utility.Filter (busted.selected_branches[busted.partition_index], '_value_', '_value_ != terms.json.attribute.test')};
}

utility.SetEnvVariable ("ASSUME_REVERSIBLE_MODELS", TRUE);

io.ReportProgressMessageMD ("BUSTED", "main", "Performing the full (dN/dS > 1 allowed) branch-site model fit");
busted.full_model =  estimators.FitLF (busted.filter_names, busted.trees, busted.model_map, busted.final_partitioned_mg_results, busted.model_object_map, {"retain-lf-object": TRUE});
io.ReportProgressMessageMD("BUSTED", "main", "* Log(L) = " + Format(busted.full_model["LogL"],8,2));
io.ReportProgressMessageMD("BUSTED", "main", "* For *test* branches, the following rate distribution for branch-site combinations was inferred");

busted.inferred_test_distribution = parameters.GetStickBreakingDistribution (busted.distribution) % 0;
selection.report_dnds (busted.inferred_test_distribution);

if (busted.has_background) {
    io.ReportProgressMessageMD("BUSTED", "main", "* For *background* branches, the following rate distribution for branch-site combinations was inferred");
    busted.inferred_background_distribution = parameters.GetStickBreakingDistribution (busted.background_distribution) % 0;
    selection.report_dnds (busted.inferred_background_distribution);

}

busted.run_test = busted.inferred_test_distribution [busted.rate_classes-1][0] > 1 && busted.inferred_test_distribution [busted.rate_classes-1][1] > 0;

(busted.json ["profiles"])["unconstrained"] = busted.ComputeSiteLikelihoods (busted.full_model["LF"]);

if (!busted.run_test) {
    io.ReportProgressMessageMD ("BUSTED", "Results", "No evidence for episodic diversifying positive selection under the unconstrained model, skipping constrained model fitting");
    busted.json ["test results"] = busted.LRT (0, 0);
} else {
    io.ReportProgressMessageMD ("BUSTED", "test", "Performing the constrained (dN/dS > 1 not allowed) model fit");
    parameters.SetConstraint (model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (terms.omega_ratio,busted.rate_classes)), "1", terms.global);
    (busted.json ["profiles"])["constrained"] = busted.ComputeSiteLikelihoods (busted.full_model["LF"]);
    Optimize (busted.MLE_H0, ^(busted.full_model["LF"]));
    (busted.json ["profiles"])["optimized null"] = busted.ComputeSiteLikelihoods (busted.full_model["LF"]);
    io.ReportProgressMessageMD ("BUSTED", "test", "Log(L) = " + Format(busted.MLE_H0[1][0],8,2));
    busted.LRT = busted.ComputeLRT (busted.full_model["LogL"], busted.MLE_H0[1][0]);
    busted.json ["test results"] = busted.LRT;
    io.ReportProgressMessageMD("BUSTED", "test", "* For *test* branches under the null (no dN/dS model), the following rate distribution for branch-site combinations was inferred");

    busted.inferred_test_distribution = parameters.GetStickBreakingDistribution (busted.distribution) % 0;
    selection.report_dnds (parameters.GetStickBreakingDistribution (busted.distribution) % 0);

    (busted.json ["evidence ratios"])["constrained"] = busted.EvidenceRatios ( (busted.json ["profiles"])["unconstrained"],  (busted.json ["profiles"])["constrained"]);
    (busted.json ["evidence ratios"])["optimized null"] = busted.EvidenceRatios ( (busted.json ["profiles"])["unconstrained"],  (busted.json ["profiles"])["optimized null"]);
}


console.log ("----\n## Branch-site Unrestricted Test for Episodic Diversification [BUSTED]");
console.log ( "Likelihood ratio test for episodic diversifying positive selection, **p = " + Format ((busted.json ["test results"])["p"], 8, 4) + "**.");


return busted.json;

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
(busted.json ["profiles"])["unconstrained"] = busted.ComputeSiteLikelihoods ("busted.LF");


if (busted_positive_class["omega"] < 1 || busted_positive_class["weight"] < 1e-8) {
    io.ReportProgressMessageMD ("BUSTED", "Results", "No evidence for positive selection under the unconstrained model, skipping constrained model fitting");
    busted.json ["test results"] = busted.ComputeLRT (0, 0);
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


lfunction busted.DistributionGuess (mean) {
    /*guess = {{Random (0,0.25),Random (0.5,0.7)}
             {Random (0.25,0.6),Random (0.1,0.2)}
             {Random (2,20),Random (0.01, 0.1)}};
    */

    guess = {{0,0.7}{0.25,0.2}{10,0.1}};

    norm = + guess[-1][1];
    guess_mean = 1/(+(guess [-1][0] $ guess [-1][1]))/norm;
    return guess["_MATRIX_ELEMENT_VALUE_*(guess_mean*(_MATRIX_ELEMENT_COLUMN_==0)+(_MATRIX_ELEMENT_COLUMN_==1)*(1/norm))"];
}



//------------------------------------------------------------------------------
function busted.ComputeSiteLikelihoods (id) {
   ConstructCategoryMatrix (_siteLike, *id, SITE_LOG_LIKELIHOODS);
   return _siteLike;
}


//------------------------------------------------------------------------------
function busted.ComputeLRT (ha, h0) {
    return {"LR" : 2*(ha-h0),
            "p" : 1-CChi2 (2*(ha-h0),2)};
}

//------------------------------------------------------------------------------------------------------------------------

function busted.EvidenceRatios (ha, h0) {
    sites = Rows (ha) * Columns (ha);
    return ha["Exp(_MATRIX_ELEMENT_VALUE_-h0[_MATRIX_ELEMENT_COLUMN_])"];
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

