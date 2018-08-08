RequireVersion ("2.3.12");


LoadFunctionLibrary     ("libv3/all-terms.bf");
LoadFunctionLibrary     ("libv3/UtilityFunctions.bf");
LoadFunctionLibrary     ("libv3/IOFunctions.bf");
LoadFunctionLibrary     ("libv3/tasks/estimators.bf");
LoadFunctionLibrary     ("libv3/tasks/alignments.bf");
LoadFunctionLibrary     ("libv3/tasks/ancestral.bf");
LoadFunctionLibrary     ("libv3/tasks/trees.bf");
LoadFunctionLibrary     ("../modules/io_functions.ibf");
LoadFunctionLibrary     ("../modules/selection_lib.ibf");
LoadFunctionLibrary     ("libv3/convenience/math.bf");
LoadFunctionLibrary     ("libv3/models/protein.bf");
LoadFunctionLibrary     ("libv3/models/protein/empirical.bf");
LoadFunctionLibrary     ("libv3/models/protein/REV.bf");
LoadFunctionLibrary     ("libv3/tasks/mpi.bf");
LoadFunctionLibrary     ("libv3/stats.bf");
LoadFunctionLibrary ("libv3/convenience/random.bf");



utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);
utility.SetEnvVariable ("ACCEPT_ROOTED_TREES", TRUE);

namespace terms.fade {
    mode = "mode";
    regimes = "regimes";
    bias = "substitution bias";
    rate = "rate multiplier";
    generator = "generator";
};


fade.parameter.bias = "FADE.bias";
fade.parameter.rate = "FADE.rate";
fade.tree.name      = "FADE.simulated_tree";

fade.settings = {};
fade.alphabet           = "ACDEFGHIKLMNPQRSTVWY";
fade.alphabet.matrix    = {1,20};
fade.simulation.matrix  = {2,20};
fade.evolutionary_modes = {"Null" : "Evolution under baseline model"};

for (r = 0; r < Abs (fade.alphabet); r += 1) {
    fade.evolutionary_modes [fade.alphabet[r]] = "Directional evolution towards `fade.alphabet[r]`";
    fade.alphabet.matrix [r] = fade.alphabet[r];
    fade.simulation.matrix[0][r] = fade.alphabet[r];
}

fade.simulation.matrix[1][0] = "1";

fade.analysis_description = {terms.io.info :
                            "A companion data simulator for FADE",
                            terms.io.version :          "0.1",
                            terms.io.reference :        "TBD",
                            terms.io.authors :          "Sergei L Kosakovsky Pond",
                            terms.io.contact :          "spond@temple.edu",
                            terms.io.requirements :     "A **rooted** phylogenetic tree with branch lengths (optionally annotated with {} to define a branch partition set)"
                          };

io.DisplayAnalysisBanner (fade.analysis_description);


// =========== LOAD DATA AND SET UP CACHES

SetDialogPrompt ("Specify a rooted tree to use for data simulations");
fade.baseline.tree = trees.LoadAnnotatedTopology (TRUE);
assert (trees.HasBranchLengths(fade.baseline.tree), "Input tree MUST have branch lengths");
assert (fade.baseline.tree[terms.trees.rooted], "Input tree MUST be rooted");



fade.replicates = io.PromptUser ("How many replicate datasets be simulated", 100, 1, 10000, true);
fade.sites_class_count = io.PromptUser ("How many types of sites will be simulated", 2, 1, 10000, true);
fade.site_classes      = {};

for (k = 0; k < fade.sites_class_count; k += 1) {
    this_class = {
                    terms.data.sites : io.PromptUser ("How many sites are in class " + (k+1), 100, 1, 10000, true),
                    terms.fade.mode : io.SelectAnOption (fade.evolutionary_modes, "Evolutionary regime for site class " + (k+1)),
                    fade.parameter.rate : io.PromptUser ("Relative overall rate for site class " + (k+1) + " (1 = average)", 1, 0, 1000, false)
                 };

    if (this_class [terms.fade.mode] != "Null") {
        this_class [fade.parameter.bias] = io.PromptUser ("Substitution bias for site class " + (k+1) + " (0 = no bias)", 1, 0, 1000, false);
    }
    fade.site_classes [k] = this_class;
}


fade.selected_branches = (selection.io.defineBranchSets ( {"0" : { terms.data.tree : fade.baseline.tree}} ))[0];

fade.settings [terms.data.tree]    = fade.baseline.tree[terms.trees.newick_annotated]; // SJS changed, should retain labeling in output.s
fade.settings [terms.json.tested]  = fade.selected_branches;
fade.settings [terms.fade.regimes] = fade.site_classes;
fade.settings [terms.replicates]   = fade.replicates;

utility.Extend               (models.protein.empirical_models, {"GTR" : "General time reversible model (189 estimated parameters)."});
fade.baseline_model         = io.SelectAnOption (models.protein.empirical_models, "Baseline substitution model");
fade.generator              = (utility.Extend (models.protein.empirical.plusF_generators , {"GTR" : "models.protein.REV.ModelDescription"}))[fade.baseline_model ];
fade.branch_lengths         = parameters.helper.tree_lengths_to_initial_values ({"0" : fade.baseline.tree}, None);

fade.settings [terms.model]          = fade.baseline_model;
fade.settings [terms.replicates]     = fade.replicates;
lfunction fade.rate.modifier (fromChar, toChar, namespace, model_type, model) {

    baseline = Call (^"fade.baseline_model.rate", fromChar,toChar, namespace, model_type, model);
    utility.EnsureKey (baseline, model_type);
    selection.io.json_store_key_value_pair (baseline, model_type, utility.getGlobalValue("terms.fade.bias"), utility.getGlobalValue("fade.parameter.bias"));
    selection.io.json_store_key_value_pair (baseline, model_type, utility.getGlobalValue("terms.fade.rate"), utility.getGlobalValue("fade.parameter.rate"));
    baseline [utility.getGlobalValue("terms.model.rate_entry")] = parameters.AppendMultiplicativeTerm (baseline [utility.getGlobalValue("terms.model.rate_entry")], utility.getGlobalValue("fade.parameter.rate"));
    if ( Type (model["fade.residue_bias"]) == "String") {
        if (toChar == model["fade.residue_bias"]) {
            baseline [utility.getGlobalValue("terms.model.rate_entry")] =
                parameters.AppendMultiplicativeTerm ( baseline [utility.getGlobalValue("terms.model.rate_entry")],
                                                      "`utility.getGlobalValue("fade.parameter.bias")`/(1-Exp (-`utility.getGlobalValue("fade.parameter.bias")`))");
         } else {
            if (fromChar == model["fade.residue_bias"]) {
                parameters.AppendMultiplicativeTerm ( baseline [utility.getGlobalValue("terms.model.rate_entry")],
                                                      "`utility.getGlobalValue("fade.parameter.bias")`/(Exp (`utility.getGlobalValue("fade.parameter.bias")`-1))");
            }
        }
    }
    return baseline;
}



lfunction fade.biased.model.generator (type, residue) {
    model = Call (^"fade.generator", type);
    utility.setGlobalValue("fade.baseline_model.rate", model[utility.getGlobalValue ("terms.model.q_ij")]);
    model[utility.getGlobalValue ("terms.model.q_ij")] = "fade.rate.modifier";
    model["fade.residue_bias"] = residue;
    model[utility.getGlobalValue ("terms.alphabet")] = utility.getGlobalValue ("fade.alphabet.matrix");
    return model;
}




fade.model.baseline = model.generic.DefineModel("fade.biased.model.generator",
            "fade.baseline_model", {
                "0": "terms.global",
                "1": None
            },
            None,
            "frequencies.equal");

fade.settings [terms.fade.generator] = fade.model.baseline;


fade.model_assignment_with_bias = {
    "fade.baseline_model" : utility.Filter (fade.selected_branches, "_value_", "_value_ == terms.tree_attributes.background"),
    "fade.biased_model" : utility.Filter (fade.selected_branches, "_value_", "_value_ == terms.tree_attributes.test"),
};

fade.model_assignment_without_bias = {
    "fade.baseline_model" : utility.Filter (fade.selected_branches, "_value_", "TRUE"),
};



parameters.DeclareGlobalWithRanges (fade.parameter.rate, 1, 0, 100);
parameters.DeclareGlobalWithRanges (fade.parameter.bias, 1e-10, 1e-10, 100);

fade.replicate_data = {};

fade.simulation_path = io.PromptUserForString ("Save simulation settings to (all simulation files will be saved to the same location with .N.fas extensions)");
fprintf (fade.simulation_path, CLEAR_FILE, fade.settings);

fade.sim_frequencies = fade.model.baseline[terms.efv_estimate];

for (fade.block_id = 0; fade.block_id < Abs (fade.site_classes); fade.block_id += 1) {
    io.ReportProgressBar  ("simulation", "Generating data for selection regime " + (fade.block_id+1) );
    fade.this_block = fade.site_classes[fade.block_id];
    if (utility.Has (  fade.this_block, fade.parameter.bias, "Number")) {
        fade.bias.residue = fade.this_block[terms.fade.mode];
        fade.model.biased = model.generic.DefineModel("fade.biased.model.generator",
                    "fade.biased_model", {
                        "0": "terms.global",
                        "1": parameters.Quote (fade.bias.residue)
                    },
                    None,
                    "frequencies.equal");

        fade.model_id_to_object = {
                    "fade.biased_model": fade.model.biased,
                    "fade.baseline_model": fade.model.baseline
                };
         parameters.SetValue (fade.parameter.bias, fade.this_block[fade.parameter.bias]);
         fade.model_assignment = fade.model_assignment_with_bias;
    } else {
         fade.model_assignment = fade.model_assignment_without_bias;
            fade.model_id_to_object = {
                     "fade.baseline_model": fade.model.baseline
                };
    }

    parameters.SetValue (fade.parameter.bias, 1e-10);
    parameters.SetValue (fade.parameter.rate, 1);
    model.ApplyModelToTree(fade.tree.name, fade.baseline.tree, None, fade.model_assignment);
    estimators.ApplyExistingEstimatesToTree (fade.tree.name, fade.model_id_to_object, (fade.branch_lengths[terms.branch_length])[0], None, {});
    parameters.SetValue (fade.parameter.rate, fade.this_block[fade.parameter.rate]);
    parameters.SetValue (fade.parameter.bias, fade.this_block[fade.parameter.bias]);


    for (fade.replicate_id = 0; fade.replicate_id < fade.replicates; fade.replicate_id += 1) {
        DataSet simulated_block = Simulate (^fade.tree.name, fade.sim_frequencies , fade.simulation.matrix, fade.this_block[terms.data.sites]);
        fade.data_block = alignments.GetAllSequences ("simulated_block");
        io.ReportProgressBar  ("simulation", "Generating data for selection regime " + (fade.block_id+1) + " replicate " + fade.replicate_id + " / " + fade.replicates);
        if (fade.block_id == 0) {
            fade.replicate_data [fade.replicate_id] = fade.data_block;
        } else {
            utility.ForEachPair (fade.data_block, "_id_", "_seq_",
            '
                (fade.replicate_data[fade.replicate_id])[_id_] += _seq_;
            ');
        }
    }


}
io.ClearProgressBar ();

for (fade.replicate_id = 0; fade.replicate_id < fade.replicates; fade.replicate_id += 1) {
    fade.current_file = fade.simulation_path + "." + (fade.replicate_id+1) + ".fas";
    console.log (fade.current_file);
    fprintf (fade.current_file, CLEAR_FILE, KEEP_OPEN);

    utility.ForEachPair ((fade.replicate_data[fade.replicate_id]), "_id_", "_seq_",
    '
        fprintf (fade.current_file, ">", _id_, "\n", _seq_, "\n");
    ');

    fprintf (fade.current_file, "\n", fade.settings [terms.data.tree], ";", CLOSE_FILE);
}
