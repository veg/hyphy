RequireVersion ("2.3.11");


LoadFunctionLibrary("libv3/all-terms.bf"); // must be loaded before CF3x4
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/models/codon.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/genetic_code.bf");
LoadFunctionLibrary("../modules/io_functions.ibf");
LoadFunctionLibrary("../modules/selection_lib.ibf");
LoadFunctionLibrary("libv3/models/codon/BS_REL.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");

busted_sim.rate_classes = 3;

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);


busted_sim.analysis_description = {terms.io.info : "This is a companion simulation module for BUSTED which can be used to generate parametric data from branch-site models.",
                           terms.io.version : "0.01",
                           terms.io.reference : "*Gene-wide identification of episodic selection*, Mol Biol Evol. 32(5):1365-71",
                           terms.io.authors : "Sergei L Kosakovsky Pond",
                           terms.io.contact : "spond@temple.edu",
                           terms.io.requirements : "A phylogenetic tree (optionally annotated with {}); either with branch lengths included, or with an attendant alignment to estimate them from"
                          };

io.DisplayAnalysisBanner (busted_sim.analysis_description);

/* load the nucleotide alignment */


busted_sim.json = {};
// need this for shared-load-file.bf

namespace busted_sim {
    LoadFunctionLibrary ("../modules/shared-load-file.bf");
    load_file ("busted_sim");
}

io.ReportProgressMessageMD('BUSTED',  'selector', 'Focal branch set consists of the following branches');
busted_sim.has_background = FALSE;

utility.ForEachPair (busted_sim.selected_branches, "_partition_", "_selection_",
    "_tested_ = utility.Filter (_selection_, '_value_', '_value_ == terms.tree_attributes.test');
     busted_sim.has_background = busted_sim.has_background || utility.Array1D(_tested_) != utility.Array1D(_selection_);
     io.ReportProgressMessageMD('BUSTED',  'selector', '* Selected ' + Abs(_tested_) + ' branches to test in the BUSTED analysis: \\\`' + Join (', ',utility.Keys(_tested_)) + '\\\`')");

// fit the GTR model to get branch lengths and

namespace busted_sim {
    doGTR ("busted_sim");
}

// define the BUSTED model

busted_sim.test.bsrel_model =  model.generic.DefineMixtureModel("models.codon.BS_REL.ModelDescription",
        "busted_sim.test", {
            "0": parameters.Quote(terms.global),
            "1": busted_sim.codon_data_info[terms.code],
            "2": parameters.Quote (busted_sim.rate_classes) // the number of rate classes
        },
        busted_sim.filter_names,
        None);



busted_sim.distribution = models.codon.BS_REL.ExtractMixtureDistribution(busted_sim.test.bsrel_model);

busted_sim.test.omega = busted.sim.propmt_for_omegas ("Omega distribution for test branches", busted_sim.rate_classes);

parameters.SetStickBreakingDistribution (busted_sim.distribution,
                                         busted_sim.test.omega
                                        );

if (busted_sim.has_background) {
    busted_sim.background.bsrel_model =  model.generic.DefineMixtureModel("models.codon.BS_REL.ModelDescription",
            "busted_sim.background", {
                "0": parameters.Quote(terms.global),
                "1": busted_sim.codon_data_info[terms.code],
                "2": parameters.Quote (busted_sim.rate_classes) // the number of rate classes
            },
            busted_sim.filter_names,
            None);

    models.BindGlobalParameters ({"0" : busted_sim.test.bsrel_model, "1" : busted_sim.background.bsrel_model}, terms.nucleotideRate("[ACGT]","[ACGT]"));
    busted_sim.background_distribution    = models.codon.BS_REL.ExtractMixtureDistribution(busted_sim.background.bsrel_model);
    busted_sim.background_omega = busted.sim.propmt_for_omegas ("Omega distribution for background branches", busted_sim.rate_classes);

    parameters.SetStickBreakingDistribution (
                                             busted_sim.background_distribution,
                                             busted_sim.background_omega
                                            );
    busted_sim.model_object_map = { "busted_sim.background" : busted_sim.background.bsrel_model,
                                    "busted_sim.test" :       busted_sim.test.bsrel_model };

    for (busted_sim.i = 1; busted_sim.i < busted_sim.rate_classes; busted_sim.i += 1) {
        parameters.SetRange (model.generic.GetGlobalParameter (busted_sim.background.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted_sim.i)), terms.range01);
    }
} else {
    busted_sim.model_object_map = { "busted_sim.test" :       busted_sim.test.bsrel_model };
}

// set up parameter constraints

for (busted_sim.i = 1; busted_sim.i < busted_sim.rate_classes; busted_sim.i += 1) {
    parameters.SetRange (model.generic.GetGlobalParameter (busted_sim.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted_sim.i)), terms.range01);
}

parameters.SetRange (model.generic.GetGlobalParameter (busted_sim.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted_sim.rate_classes)), terms.range_gte1);

busted_sim.model_map = {};

for (busted_sim.partition_index = 0; busted_sim.partition_index < busted_sim.partition_count; busted_sim.partition_index += 1) {

    busted_sim.model_map + { "busted_sim.test" : utility.Filter (busted_sim.selected_branches[busted_sim.partition_index], '_value_', '_value_ == terms.tree_attributes.test'),
					         "busted_sim.background" : utility.Filter (busted_sim.selected_branches[busted_sim.partition_index], '_value_', '_value_ != terms.tree_attributes.test')};
}

utility.SetEnvVariable ("ASSUME_REVERSIBLE_MODELS", TRUE);

busted_sim.settings =
    estimators.BuildLFObject ("busted_sim.lf", busted_sim.filter_names, busted_sim.trees, busted_sim.model_map, busted_sim.gtr_results, busted_sim.model_object_map, None)
;

busted_sim.replicates = io.PromptUser ("How many replicates should be generated", 100, 1, 1e6, TRUE);

SetDialogPrompt ("Save simulation settings to (this file will also serve as the base-path for .1, .2, etc simulated alignments)");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, busted_sim.settings);

busted_sim.basepath = LAST_FILE_PATH;

for (busted_sim.i = 0; busted_sim.i < busted_sim.replicates; busted_sim.i += 1) {
    DataSet busted_sim.replicate = SimulateDataSet (busted_sim.lf);
    DataSetFilter busted_sim.replicate.filter = CreateFilter (busted_sim.replicate, 1);

    fprintf (busted_sim.basepath + "." + (busted_sim.i + 1), CLEAR_FILE, busted_sim.replicate.filter);
}



//===========================================================================================================================================

lfunction busted.sim.propmt_for_omegas (prefix, classes) {
    weight_left = 1;
    result = {classes,2};
    for (k = 0; k < classes; k += 1) {
        result [k][0]  = io.PromptUser ("<" + prefix + "> rate for category " + (k + 1), 1, 0, 1 + 9999 * (k == classes - 1), FALSE);

        if (k < classes - 1) {
            result [k][1]  = io.PromptUser ("<" + prefix + "> weight for category " + (k + 1), 1, 0, weight_left, FALSE);
            weight_left = weight_left - result [k][1];
        } else {
            result [k][1] = weight_left;
        }
    }
    return result;
}
