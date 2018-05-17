RequireVersion ("2.2.2.0141023");


console.log("ERROR: THIS BATCHFILE WILL NOT RUN PROPERLY WITH UPDATED TERMS BASE. QUITTING.");
exit();


VERBOSITY_LEVEL				= 0;
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


LoadFunctionLibrary("BranchSiteTemplate");
/*------------------------------------------------------------------------------
    BranchSiteTemplate Defines

    BuildCodonFrequencies (obsF);
    PopulateModelMatrix (ModelMatrixName&, EFV, synrateP, globalP, nonsynRateP);

------------------------------------------------------------------------------*/

RELAX.settings = {"GTR" : 1,
                  "LocalMG" : 1,
                  "Estimate GTR" : 1};


RELAX.relaxation_parameter        = "RELAX.K";
term.RELAX.k          = "relaxation coefficient";
RELAX.test            = "RELAX.test";
RELAX.reference       = "RELAX.reference";

/*------------------------------------------------------------------------------*/


io.DisplayAnalysisBanner ({"info" : "RELAX (a random effects test of selection relaxation)
                            uses a random effects branch-site model framework
                            to test whether one alignment ('Test') evolves under relaxed
                            selection relative to a the 'Reference' alignment (R), as measured
                            by the relaxation parameter (K).",
                           "version" : "1.0",
                           "reference" : "RELAX: Detecting Relaxed Selection in a Phylogenetic Framework (2015). Mol Biol Evol 32 (3): 820-832",
                           "authors" : "Sergei L Kosakovsky Pond and Temple iGEM / UCSD viral evolution group",
                           "contact" : "spond@temple.edu",
                           "requirements" : "two in-frame codon alignments and phylogenetic trees (one to be treated as reference, and the other as test)"
                          } );




relax.codon_data_info = {};
relax.tree_info       = {};

SetDialogPrompt ("Load the reference alignment");


relax.codon_data_info     + alignments.PromptForGeneticCodeAndAlignment ("RELAX.codon_data.reference", "RELAX.codon_filter.reference");
io.ReportProgressMessage ("RELAX", "Loaded the reference MSA with " + (relax.codon_data_info[0])["sequences"] + " sequences and " + (relax.codon_data_info[0])["sites"] + " codons from '" + (relax.codon_data_info[0])["file"] + "'");
relax.tree_info + trees.LoadAnnotatedTopology (1);


SetDialogPrompt ("Load the test alignment");
relax.codon_data_info     + alignments.PromptForGeneticCodeAndAlignment ("RELAX.codon_data.test", "RELAX.codon_filter.test");

io.ReportProgressMessage ("RELAX", "Loaded the test MSA with " + (relax.codon_data_info[1])["sequences"] + " sequences and " + (relax.codon_data_info[1])["sites"] + " codons from '" + (relax.codon_data_info[1])["file"] + "'");

relax.sample_size         = (relax.codon_data_info[0])["sites"] * (relax.codon_data_info[0])["sequences"] +
                            (relax.codon_data_info[1])["sites"] * (relax.codon_data_info[1])["sequences"];

relax.tree_info + trees.LoadAnnotatedTopology (1);


relax.codon_lists = models.codon.MapCode ((relax.codon_data_info[0])["code"]);

relax.codon_frequencies_reference =
frequencies._aux.CF3x4(frequencies._aux.empirical.collect_data ("RELAX.codon_filter.reference",3,1,1),
models.DNA.alphabet, relax.codon_lists["sense"], relax.codon_lists["stop"]);

relax.codon_frequencies_test =
frequencies._aux.CF3x4(frequencies._aux.empirical.collect_data ("RELAX.codon_filter.test",3,1,1),
models.DNA.alphabet, relax.codon_lists["sense"], relax.codon_lists["stop"]);


relax.taskTimerStart (1);

relax.codon_filters = {{"RELAX.codon_filter.reference","RELAX.codon_filter.test"}};


if (RELAX.settings["GTR"]) {
    io.ReportProgressMessage ("RELAX", "Obtaining branch lengths under the GTR model");
    relax.gtr_results = estimators.FitGTR     (relax.codon_filters, relax.tree_info, None);
    io.ReportProgressMessage ("RELAX", "Log(L) = " + relax.gtr_results["LogL"]);
    estimators.fixSubsetOfEstimates (relax.gtr_results, relax.gtr_results["global"]);
} else {
    relax.gtr_results = None;
}

RELAX.reference.model      = relax.io.define_a_bsrel_model (RELAX.reference, relax.codon_frequencies_reference, None ,FALSE);
RELAX.test.model           = relax.io.define_a_bsrel_model (RELAX.test, relax.codon_frequencies_test, None , FALSE);

model.ApplyModelToTree              ("RELAX.tree.reference", relax.tree_info[0], {"default": RELAX.reference.model}, None);
model.ApplyModelToTree              ("RELAX.tree.test", relax.tree_info[1], {"default": RELAX.test.model}, None);

((RELAX.test.model["parameters"])["global"])[term.RELAX.k] = RELAX.relaxation_parameter;

parameters.DeclareGlobal (RELAX.relaxation_parameter, {});

relax.set_up_model_rate_constraints (RELAX.reference.model, RELAX.test.model, RELAX.relaxation_parameter);

ASSUME_REVERSIBLE_MODELS = 1;
LikelihoodFunction relax.LF = (RELAX.codon_filter.reference, RELAX.tree.reference, RELAX.codon_filter.test, RELAX.tree.test);

RELAX.model_map = {};
RELAX.model_map [RELAX.reference.model["id"]] = RELAX.reference.model;
RELAX.model_map [RELAX.test.model["id"]] = RELAX.test.model;

RELAX.previous_results = relax.gtr_results;

do {

    estimators.ApplyExistingEstimates   ("relax.LF",
        RELAX.model_map,
        RELAX.previous_results,
        None);




    utility.SetEnvVariable ("VERBOSITY_LEVEL",1);
    utility.SetEnvVariable ("USE_LAST_RESULTS",1);

    // utility.SetEnvVariable ("LF_SMOOTHING_SCALER", 0.1);

    io.ReportProgressMessage ("RELAX", "Fitting the alternative RELAX model");

    Optimize (relax.alt.mle, relax.LF);

    relax.alt = estimators.ExtractMLEs ("relax.LF", RELAX.model_map);
    relax.add_scores (relax.alt, relax.alt.mle);

    fprintf (stdout, "\nReference distribution\n");
    relax.print_distribution (relax.getRateDistribution (RELAX.reference.model, 1));

    fprintf (stdout, "\nTest distribution\n");
    relax.print_distribution (relax.getRateDistribution (RELAX.test.model, 1));


    io.ReportProgressMessage ("RELAX", "Alternative RELAX model fit. log(L) = " + relax.alt["LogL"] + ". Relaxation parameter K = " + Eval (RELAX.relaxation_parameter));
    io.ReportProgressMessage ("RELAX", "Fitting the null RELAX model");
    parameters.SetConstraint (RELAX.relaxation_parameter, "1", "");

    Optimize (relax.null.mle, relax.LF);
    relax.null = estimators.ExtractMLEs ("relax.LF", RELAX.model_map);
    relax.add_scores (relax.null, relax.null.mle);

    fprintf (stdout, "\nReference (=test) distribution\n");
    relax.print_distribution (relax.getRateDistribution (RELAX.reference.model, 1));

    io.ReportProgressMessage ("RELAX", "Null RELAX model fit. log(L) = " + relax.null["LogL"] + ". p-value = " + (relax.runLRT (relax.alt["LogL"], relax.null["LogL"]))["p"]);

    parameters.RemoveConstraint (RELAX.relaxation_parameter);

    RELAX.previous_results = relax.null;

} while (relax.null["LogL"] > relax.alt["LogL"]);


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
        if (k < 3) {
            ((model_parameters["parameters"])["global"])[relax.define_weight_term (k)] = (model_parameters["f"])[k-1];
        }

        model_parameters["Q"] + ("Q_`id`_" + k);
        if (do_local) {
            PopulateModelMatrix			  ((model_parameters["Q"])[k-1],  relax.nuc, "t", "Min (1000," + (model_parameters["omegas"])[k-1] +"^local_k)", "");
        } else {
            PopulateModelMatrix			  ((model_parameters["Q"])[k-1],  relax.nuc, "t", (model_parameters["omegas"])[k-1], "");
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
    relax.aux.copy_branch_length.df = 0;


    if (Abs (RELAX.proportional_constraint)) {
        Eval ("`parameter`.`relax.aux.copy_branch_length.t` := `RELAX.proportional_constraint` * " + value);
        relax.aux.copy_branch_length.df += 1;
    } else {
        Eval ("`parameter`.`relax.aux.copy_branch_length.t` = " + value);
    }

    if (Type (relax.aux.copy_branch_length.k) == "String") {
        Eval ("`parameter`.`relax.aux.copy_branch_length.k` = 1");
    }

    return relax.aux.copy_branch_length.df;
}

//------------------------------------------------------------------------------------------------------------------------

function relax.set_up_model_rate_constraints (reference, test, constrain_parameter) {
    for (relax.set_up_model_rate_constraints.i = 1; relax.set_up_model_rate_constraints.i <= 3; relax.set_up_model_rate_constraints.i += 1) {
        parameters.SetConstraint (((test["parameters"])["global"])[relax.define_omega_term (relax.set_up_model_rate_constraints.i)],
                                  ((reference["parameters"])["global"])[relax.define_omega_term (relax.set_up_model_rate_constraints.i)] + "^" + constrain_parameter,
                                  "");

    }
    for (relax.set_up_model_rate_constraints.i = 1; relax.set_up_model_rate_constraints.i < 3; relax.set_up_model_rate_constraints.i += 1) {
        parameters.SetConstraint (((test["parameters"])["global"])[relax.define_weight_term (relax.set_up_model_rate_constraints.i)],
                                  ((reference["parameters"])["global"])[relax.define_weight_term (relax.set_up_model_rate_constraints.i)],
                                  "");

    }
}


//------------------------------------------------------------------------------------------------------------------------

lfunction relax.define_omega_term (cat) {
    return "Omega for category " + cat;
}

lfunction relax.define_weight_term (cat) {
    return "Omega frequency parameter " + cat;
}

lfunction relax.print_distribution (d) {
    for (i = 0; i < Rows (d); i += 1) {
        fprintf (stdout, "\nomega_", i+1, " = ", d[i][0], " (p = ", d[i][1], ")");
    }
    fprintf (stdout, "\n");
}
