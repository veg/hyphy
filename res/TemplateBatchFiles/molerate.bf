RequireVersion ("2.5.48");


LoadFunctionLibrary("libv3/all-terms.bf"); 
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");

LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary("SelectionAnalyses/modules/selection_lib.ibf");

LoadFunctionLibrary     ("libv3/models/protein.bf");
LoadFunctionLibrary     ("libv3/models/protein/empirical.bf");
LoadFunctionLibrary     ("libv3/models/protein/REV.bf");
LoadFunctionLibrary     ("libv3/models/rate_variation.bf");


utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", FALSE);
utility.SetEnvVariable ("ACCEPT_ROOTED_TREES", TRUE);

rer.analysis_description = {terms.io.info : 
"Perform a relative branch length (RER) test on protein sequences",
                               terms.io.version : "0.0.1",
                               terms.io.authors : "Sergei L Kosakovsky Pond",
                               terms.io.contact : "spond@temple.edu",
                               terms.io.requirements : "a protein alignment, a phylogenetic tree with relative branch lengths, and a list of designated branches to include in the test set"
                              };

io.DisplayAnalysisBanner (fitter.analysis_description);

  
  


namespace terms.rer {
    test         = "test";
    background   = "background";
    proportional = "Proportional";
    proportional_partitoned = "Proportional Partitioned";
    unconstrained_test = "Unconstrained Test";
    unconstrained = "Unconstrained";
    reference     = "Reference";
    strategy      = "labeling strategy";
    model         = "model";
    rv            = "rate variation";
    branch_level  = "branch level analysis"
}

rer.json    = { 
                   terms.json.analysis: rer.analysis_description,
                   terms.json.input: {},
                   terms.json.fits : {},
                   terms.json.timers : {},
                 };
 
selection.io.startTimer (rer.json [terms.json.timers], "Overall", 0);

KeywordArgument ("type",        "The type of data to perform screening on", "protein");
rer.dataType = io.SelectAnOption  ({"nucleotide" : "A nucleotide (DNA/RNA) alignment",
                                      "protein" : "A protein alignment"},
                                      "The type of data to perform screening on");
          
KeywordArgument ("alignment",            "A protein multiple sequence alignment in one of the formats supported by HyPhy (single partition)");
                            
if (rer.dataType == "protein") {
    rer.alignment = alignments.ReadProteinDataSet ("rer.sequences", null);
} else {
    rer.alignment = alignments.ReadNucleotideDataSet ("rer.sequences", null);
}


rer.sample_size = rer.alignment[utility.getGlobalValue("terms.data.sequences")] * rer.alignment[utility.getGlobalValue("terms.data.sites")];

io.ReportProgressMessage ("", ">Loaded a multiple sequence alignment with **" + rer.alignment[utility.getGlobalValue("terms.data.sequences")] + "** sequences, **" + rer.alignment[utility.getGlobalValue("terms.data.sites")] + "** sites from \`" + rer.alignment[utility.getGlobalValue("terms.data.file")] + "\`");

alignments.StoreInJSON (rer.json,rer.alignment,None);

io.CheckAssertion("Rows (rer.alignment[utility.getGlobalValue('terms.data.partitions')]) == 1", "RERConverge only works on a single partition dataset");

KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'MG94.json')",  rer.alignment [terms.data.file] + "molerate.json");

rer.alignment [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");


KeywordArgument ("tree",                 "A phylogenetic tree with branch lengths");
rer.tree =  trees.LoadAnnotatedTopology (FALSE);

io.CheckAssertion("trees.HasBranchLengths (rer.tree)", "The input tree must have fully specified branch lengths");
    
    
rer.original_tree_names = trees.LeafNames (rer.tree);
rer.select_branch_options = {};
rer.branches = {};

for (k; in; rer.original_tree_names) {
    rer.select_branch_options[k] = k;
    rer.branches [k] = terms.rer.background;
}
    
rer.in_tree = rer.array2dict(rer.original_tree_names);  
rer.in_data = rer.array2dict(alignments.GetSequenceNames ("rer.sequences"));

/** 
    compute the intersection of sequences in the tree and the dataset
    this is done in a CASE-INSENSITIVE manner
**/

rer.in_tree_not_in_data = rer.in_tree;
rer.in_tree_not_in_data - rer.in_data;

rer.in_data_not_in_tree = rer.in_data;
rer.in_data_not_in_tree - rer.in_tree;

rer.shared = rer.intersect (rer.in_tree, rer.in_data);

io.CheckAssertion("utility.Array1D (rer.shared)>=3", "There must be 3 or more sequences that appear both in the alignment and in the phylogenetic tree");

Topology T = rer.tree[terms.trees.newick_with_lengths];

if (utility.Array1D (rer.in_tree_not_in_data)) {
    io.ReportProgressMessage ("", ">Trimming the following **" + utility.Array1D (rer.in_tree_not_in_data) + "** tips from the tree because they are missing from the alignment: \`" + Join (", ", Rows(rer.in_tree_not_in_data)) + "\`");

    Topology T = rer.tree[terms.trees.newick_with_lengths];
    rer.to_remove = {};
    for (n; in; T) {
        if (rer.in_tree_not_in_data[n&&1]) {
            rer.to_remove[n] = 1;
        }
    }
    T - Rows (rer.to_remove);
    rer.tree = trees.ExtractTreeInfo (Format (T,1,1));
}

if (utility.Array1D (rer.in_data_not_in_tree)) {
    io.ReportProgressMessage ("", ">Trimming the following **" + utility.Array1D (rer.in_data_not_in_tree) + "** tips from the alignment because they are missing from the tree: \`" + Join (", ", Rows(rer.in_data_not_in_tree)) + "\`");
}

function rer.sequence.filter (id, sequence) {
    return rer.shared [id&&1];
}

DataSetFilter rer.filter = CreateFilter (rer.sequences, 1, "", "rer.sequence.filter");

io.ReportProgressMessage ("", ">Retaining **" + rer.filter.species + "** sequences for the analysis");

KeywordArgument ("branches",             "Designated lineages to test");
KeywordArgument ("labeling-strategy",    "Labeling strategy for internal nodes", "all-descendants");

rer.branch_set = io.MultiSelectOptions (
    rer.select_branch_options,
    "Designated lineages to test",
    -1
);


rer.branches = {};



for (k;in;rer.branch_set) {
    rer.branches [k] = terms.rer.test;
}


rer.branch_count = utility.Array1D (rer.branch_set);

io.CheckAssertion("rer.branch_count > 0 && rer.branch_count < rer.filter.species", "RERConverge requires that both test and background partition be non-empty");
 
 
io.ReportProgressMessage ("", ">Selected **" + rer.branch_count + "** tips as designated lineages: \`" + Join (", ", rer.branch_set) + "\`");

rer.labeling_option = io.SelectAnOption (
    {   
        "all-descendants" : "Label an internal node if an only if ALL of its descendants are labeled" ,
        "some-descendants" : "Label an internal node if an only if SOME of its descendants are labeled",
        "parsimony" : "Label internal nodes using parsimony",
        "none" : "Do not label any internal nodes"  
    },
"Labeling strategy for internal nodes");

(rer.json [terms.json.analysis])[terms.rer.strategy] = rer.labeling_option;

Topology rer.T = rer.tree[terms.trees.newick_with_lengths];

rer.labels = None;

if (rer.labeling_option == "all-descendants") {
    rer.labels = trees.ConjunctionLabel  ("rer.T", rer.branches);
}

if (rer.labeling_option == "some-descendants") {
    rer.labels = trees.DisjunctionLabel  ("rer.T", rer.branches);
} 

if (rer.labeling_option == "parsimony") {
    rer.labels = trees.ParsimonyLabel  ("rer.T", rer.branches);
} 

for (k; in; T) {
    if (rer.branches / k == FALSE) {
        rer.branches [k] = terms.rer.background;
    }
}


rer.tree[terms.trees.model_map] * rer.branches;

if (None != rer.labels) {
    rer.labels = rer.labels["labels"];
    io.ReportProgressMessage ("", ">Additionally labeled **" + utility.Array1D (rer.labels) + "** internal modes as designated lineages: \`" + Join (", ", Rows(rer.labels)) + "\`");
    rer.tree[terms.trees.model_map] *  rer.labels;
    rer.branches * rer.labels;
} 

trees.store_tree_information (rer.json, {"0" : rer.tree}, rer.tree[terms.trees.model_map]);
       
if (rer.dataType == "protein") {
    utility.Extend (models.protein.empirical_models, {"GTR" : "General time reversible model (189 estimated parameters)."});
    KeywordArgument ("model", "The substitution model to use", "WAG");
    rer.substitution_model         = io.SelectAnOption (models.protein.empirical_models, "Baseline substitution model");
    rer.model.generator        = (utility.Extend (models.protein.empirical.plusF_generators , {"GTR" : "models.protein.REV.ModelDescription"}))[rer.substitution_model];
} else {
    rer.substitution_model = "Nucleotide+GTR";
    rer.model.generator = "models.DNA.GTR.ModelDescription"; 
}
       
(rer.json [terms.json.analysis])[terms.rer.model] = rer.substitution_model;
rer.json [terms.json.runtime] = utility.GetVersion (FALSE);

KeywordArgument ("rv",   "Site to site rate variation", "None");
rer.rateVariation = io.SelectAnOption  ({"None"  : "Constant rates",
                                          "Gamma" : "Unit mean gamma distribution discretized into N rates",
                                          "GDD"   : "General discrete distribution on N rates"},
                                          "Site to site rate variation option");



if (rer.rateVariation != "None") {
    KeywordArgument ("rate-classes",   "How many site rate classes to use", "4");
    rer.rateClasses = io.PromptUser(">How many site rate classes to use", 4, 2, 10, TRUE);
    rer.model.generator.base = rer.model.generator;
    if (gard.rateVariation == "Gamma") {
        rer.model.generator = "rer.model.withGamma";
    } else {
        rer.model.generator = "rer.model.withGDD";
    }
}

KeywordArgument ("full-model",   "Fit the full unconstrained model", "Yes");
rer.fit_full_model = io.SelectAnOption  ({"Yes"  : "Fit the model with all branch lengths unconstrained",
                                          "No" : "Do NOT fit the model with all branch lengths unconstrained"},
                                          "Fit the full unconstrained model") == "Yes";



if (rer.branch_count > 1) {
    KeywordArgument ("branch-level-analysis",   "Perform test clade branch-level testing", "No");
    rer.outlier = io.SelectAnOption  ({"No"  : "Do NOT perform branch-level testing",
                                              "Yes" : "Perform  branch-level testing ()"},
                                              "Perform test clade branch-level testing") == "Yes";
} else {
    rer.outlier = FALSE;
}


(rer.json [terms.json.analysis])[terms.rer.rv] = rer.rateVariation;

rer.given_lengths = parameters.helper.tree_lengths_to_initial_values ({"0":rer.tree}, None);

rer.scalers = {
                 terms.rer.background : "rer.global_scaler",
                 terms.rer.test : "rer.global_scaler_test"
              };

for (ignore, param; in; rer.scalers) {
   parameters.DeclareGlobal(param, None);
   parameters.SetValue(param, 1);
}

rer.user_lengths = rer.process_branch_lengths (rer.given_lengths,rer.branches);


selection.io.json_store_branch_attribute(rer.json, terms.rer.reference, terms.branch_length, 0,
                                             0,
                                             selection.io.extract_branch_info((rer.given_lengths[terms.branch_length])[0],
                                              "selection.io.branch.length"));                  


// SET UP THE OUTPUT TABLE
// -----------------------

io.ReportProgressMessageMD ("rer", "fitting", "Fitting different models of branch length constraint to the data");
io.ReportProgressMessageMD ("rer", "fitting", "");


rer.table_output_options = {
        utility.getGlobalValue("terms.table_options.header"): 1,
        utility.getGlobalValue("terms.table_options.minimum_column_width"): 16,
        utility.getGlobalValue("terms.table_options.align"): "left",
        utility.getGlobalValue("terms.table_options.column_widths"): {
            "0": 40,
            "1" :20,
            "2" :12,
            "3" :12,
            "4" :12,
            "5" :12,
            "6" :15
        }
    };
    
rer.header = {
        7,
        1
    };
    
rer.header[0] = "Model";
rer.header[1] = "L(test)";
rer.header[2] = "Rel to ref";
rer.header[3] = "L(background)";
rer.header[4] = "Rel to ref";
rer.header[5] = "Log(L)";
rer.header[6] = "AIC-c";

fprintf(stdout, io.FormatTableRow(rer.header, rer.table_output_options));
rer.table_output_options[utility.getGlobalValue("terms.table_options.header")] = FALSE;

rer.print_row = {
            7,
            1
        };
                
rer.ref_lengths = {
    terms.rer.test : +rer.user_lengths[terms.rer.test],
    terms.rer.background : +rer.user_lengths[terms.rer.background]
};

rer.print_row [0] = terms.rer.reference;
rer.print_row [1] = Format (rer.ref_lengths[terms.rer.test],10,4);
rer.print_row [2] = 1;
rer.print_row [3] = Format (rer.ref_lengths[terms.rer.background],10,4);
rer.print_row [4] = 1;
rer.print_row [5] = "N/A";
rer.print_row [6] = "N/A";

fprintf(stdout, io.FormatTableRow(rer.print_row, rer.table_output_options));

// PROPOTIONAL PARTITIONED
// -----------------------
selection.io.startTimer (rer.json [terms.json.timers], terms.rer.proportional_partitoned, 1);

rer.fit_proportional_partitioned = estimators.FitSingleModel_Ext (
    "rer.filter",
    rer.tree,
    rer.model.generator,
    rer.given_lengths,
    {
        terms.run_options.retain_lf_object : TRUE,
        terms.run_options.apply_user_constraints : "rer.branch_length_processor",
        terms.run_options.retain_model_object : TRUE
    }
);

selection.io.json_store_lf (rer.json,
                            terms.rer.proportional_partitoned ,
                            rer.fit_proportional_partitioned[terms.fit.log_likelihood],
                            rer.fit_proportional_partitioned[terms.parameters],
                            rer.sample_size,
                            utility.Map (rer.fit_proportional_partitioned[terms.global], "_value_", '_value_ [terms.fit.MLE]'), 
                            0);
                            
 
selection.io.json_store_branch_attribute(rer.json, terms.rer.proportional_partitoned, terms.branch_length, 1,
                                             0,
                                             selection.io.extract_branch_info((rer.fit_proportional_partitioned[terms.branch_length])[0],
                                              "selection.io.branch.length"));                  

ref.lf_id = rer.fit_proportional_partitioned[terms.likelihood_function];
rer.lengths_proportional_partitioned = rer.process_branch_lengths (rer.fit_proportional_partitioned, rer.branches);

selection.io.stopTimer (rer.json [terms.json.timers], terms.rer.proportional_partitoned);

rer.empirical_parameters = ((utility.GetFirstDictElement(rer.fit_proportional_partitioned[terms.model]))[terms.parameters])[terms.model.empirical];

rer.print_row [0] = terms.rer.proportional_partitoned;
rer.BL = +rer.lengths_proportional_partitioned[terms.rer.test];
rer.print_row [1] = Format (rer.BL,10,4);
rer.print_row [2] = Format (rer.BL/rer.ref_lengths[terms.rer.test],6,3);
rer.BL = +rer.lengths_proportional_partitioned[terms.rer.background];
rer.print_row [3] = Format (rer.BL,10,4);
rer.print_row [4] = Format (rer.BL/rer.ref_lengths[terms.rer.background],6,3);
rer.print_row [5] = Format (rer.fit_proportional_partitioned[terms.fit.log_likelihood],6,3);
rer.print_row [6] = Format (selection.io.getIC(rer.fit_proportional_partitioned[terms.fit.log_likelihood], rer.fit_proportional_partitioned[terms.parameters], rer.sample_size),6,3);

fprintf(stdout, io.FormatTableRow(rer.print_row, rer.table_output_options));

// PROPOTIONAL 
// -----------

parameters.SetConstraint (rer.scalers[terms.rer.background],rer.scalers[terms.rer.test],"");
^(rer.scalers[terms.rer.test]) = Eval (rer.scalers[terms.rer.test]);

selection.io.startTimer (rer.json [terms.json.timers], terms.rer.proportional, 2);
rer.fit_proportional = estimators.FitExistingLF (ref.lf_id, rer.fit_proportional_partitioned[terms.model]);
rer.lengths_proportional = rer.process_branch_lengths (rer.fit_proportional, rer.branches);

selection.io.json_store_lf (rer.json,
                            terms.rer.proportional ,
                            rer.fit_proportional[terms.fit.log_likelihood],
                            rer.fit_proportional[terms.parameters] + rer.empirical_parameters,
                            rer.sample_size,
                            utility.Map (rer.fit_proportional[terms.global], "_value_", '_value_ [terms.fit.MLE]'), 
                            1);
                            
 
selection.io.json_store_branch_attribute(rer.json, terms.rer.proportional, terms.branch_length, 2,
                                             0,
                                             selection.io.extract_branch_info((rer.fit_proportional[terms.branch_length])[0],
                                              "selection.io.branch.length"));                  

rer.print_row [0] = terms.rer.proportional;
rer.BL = +rer.lengths_proportional[terms.rer.test];
rer.print_row [1] = Format (rer.BL,10,4);
rer.print_row [2] = Format (rer.BL/rer.ref_lengths[terms.rer.test],6,3);
rer.BL = +rer.lengths_proportional[terms.rer.background];
rer.print_row [3] = Format (rer.BL,10,4);
rer.print_row [4] = Format (rer.BL/rer.ref_lengths[terms.rer.background],6,3);
rer.print_row [5] = Format (rer.fit_proportional[terms.fit.log_likelihood],6,3);
rer.print_row [6] = Format (selection.io.getIC(rer.fit_proportional[terms.fit.log_likelihood], rer.fit_proportional[terms.parameters] + rer.empirical_parameters, rer.sample_size),6,3);

fprintf(stdout, io.FormatTableRow(rer.print_row, rer.table_output_options));


selection.io.stopTimer (rer.json [terms.json.timers], terms.rer.proportional);

// UNCONSTRAINED TEST
// ------------------


ref.parameter_set = estimators.TraverseLocalParameters (ref.lf_id, rer.fit_proportional_partitioned[terms.model], "rer.collect_constraints");

if (rer.outlier ) {
    rer.by_test_branch = {};
    for (n,p; in; ref.parameter_set ) {
        if (rer.branches[n] == terms.rer.test) {
            rer.stash_constraint = parameters.GetConstraint (p);
            parameters.RemoveConstraint (p);
            rer.fit_test_unconstrained = estimators.FitExistingLF (ref.lf_id, rer.fit_proportional_partitioned[terms.model]);
            rer.by_test_branch[n] = {terms.fit.log_likelihood : rer.fit_test_unconstrained[terms.fit.log_likelihood]};
            parameters.SetConstraint (p, rer.stash_constraint,"");
        }
    }    
}


for (n,p; in; ref.parameter_set ) {
    if (rer.branches[n] == terms.rer.test) {
        parameters.RemoveConstraint (p);
    }
}

selection.io.startTimer (rer.json [terms.json.timers], terms.rer.unconstrained_test, 3);

rer.fit_test_unconstrained = estimators.FitExistingLF (ref.lf_id, rer.fit_proportional_partitioned[terms.model]);


selection.io.json_store_lf (rer.json,
                            terms.rer.unconstrained_test ,
                            rer.fit_test_unconstrained[terms.fit.log_likelihood],
                            rer.fit_test_unconstrained[terms.parameters] + rer.empirical_parameters,
                            rer.sample_size,
                            utility.Map (rer.fit_test_unconstrained[terms.global], "_value_", '_value_ [terms.fit.MLE]'), 
                            2);
                            
selection.io.json_store_branch_attribute(rer.json, terms.rer.unconstrained_test, terms.branch_length, 3,
                                             0,
                                             selection.io.extract_branch_info((rer.fit_test_unconstrained[terms.branch_length])[0],
                                              "selection.io.branch.length"));                  

rer.lengths_test_unconstrained  = rer.process_branch_lengths (rer.fit_test_unconstrained, rer.branches);

rer.print_row [0] = terms.rer.unconstrained_test;
rer.BL = +rer.lengths_test_unconstrained[terms.rer.test];
rer.print_row [1] = Format (rer.BL,10,4);
rer.print_row [2] = Format (rer.BL/rer.ref_lengths[terms.rer.test],6,3);
rer.BL = +rer.lengths_test_unconstrained[terms.rer.background];
rer.print_row [3] = Format (rer.BL,10,4);
rer.print_row [4] = Format (rer.BL/rer.ref_lengths[terms.rer.background],6,3);
rer.print_row [5] = Format (rer.fit_test_unconstrained[terms.fit.log_likelihood],6,3);
rer.print_row [6] = Format (selection.io.getIC(rer.fit_test_unconstrained[terms.fit.log_likelihood], rer.fit_test_unconstrained[terms.parameters] + rer.empirical_parameters, rer.sample_size),6,3);
fprintf(stdout, io.FormatTableRow(rer.print_row, rer.table_output_options));

if (rer.outlier ) {
   rer.null_ll = rer.fit_proportional[terms.fit.log_likelihood];
   rer.alt_ll = rer.fit_test_unconstrained[terms.fit.log_likelihood];
   rer.outlier_pvalues = {};
   
   rer.i = 0;
   for (n,p; in; rer.by_test_branch) {
        (rer.by_test_branch[n])[terms.Null] = math.DoLRT (rer.null_ll,p[terms.fit.log_likelihood], 1);
        rer.outlier_pvalues[rer.i] = ((rer.by_test_branch[n])[terms.Null])[terms.p_value];
        rer.i += 1;
        (rer.by_test_branch[n])[terms.alternative] = math.DoLRT (p[terms.fit.log_likelihood], rer.alt_ll, rer.branch_count-1);
        rer.outlier_pvalues[rer.i] = ((rer.by_test_branch[n])[terms.alternative])[terms.p_value];
        rer.i += 1;
   }
   rer.outlier_pvalues = math.HolmBonferroniCorrection(rer.outlier_pvalues);
   rer.i = 0;
   for (n,p; in; rer.by_test_branch) {
        ((rer.by_test_branch[n])[terms.Null])[terms.json.corrected_pvalue] = rer.outlier_pvalues[rer.i];
        rer.i += 1;
        ((rer.by_test_branch[n])[terms.alternative])[terms.json.corrected_pvalue] = rer.outlier_pvalues[rer.i];
        rer.i += 1;
  }   
}


selection.io.stopTimer (rer.json [terms.json.timers], terms.rer.unconstrained_test);

if (rer.fit_full_model) {
    // UNCONSTRAINED
    // -------------
    
    for (n,p; in; ref.parameter_set ) {
        if (rer.branches[n] != terms.rer.test) {
            parameters.RemoveConstraint (p);
        }
    }
    
    selection.io.startTimer (rer.json [terms.json.timers], terms.rer.unconstrained, 4);
    
    rer.fit_full_unconstrained = estimators.FitExistingLF (ref.lf_id, rer.fit_proportional_partitioned[terms.model]);
    
    selection.io.stopTimer (rer.json [terms.json.timers], terms.rer.unconstrained);
    
    selection.io.json_store_lf (rer.json,
                                terms.rer.unconstrained ,
                                rer.fit_full_unconstrained[terms.fit.log_likelihood],
                                rer.fit_full_unconstrained[terms.parameters],
                                rer.sample_size,
                                utility.Map (rer.fit_full_unconstrained[terms.global], "_value_", '_value_ [terms.fit.MLE]'), 
                                3);
                                
    selection.io.json_store_branch_attribute(rer.json, terms.rer.unconstrained, terms.branch_length, 4,
                                                 0,
                                                 selection.io.extract_branch_info((rer.fit_full_unconstrained[terms.branch_length])[0],
                                                  "selection.io.branch.length"));       
                                                             
    rer.lengths_unconstrained  = rer.process_branch_lengths (rer.fit_full_unconstrained, rer.branches);
    rer.print_row [0] = terms.rer.unconstrained;
    rer.BL = +rer.lengths_unconstrained[terms.rer.test];
    rer.print_row [1] = Format (rer.BL,10,4);
    rer.print_row [2] = Format (rer.BL/rer.ref_lengths[terms.rer.test],6,3);
    rer.BL = +rer.lengths_unconstrained[terms.rer.background];
    rer.print_row [3] = Format (rer.BL,10,4);
    rer.print_row [4] = Format (rer.BL/rer.ref_lengths[terms.rer.background],6,3);
    rer.print_row [5] = Format (rer.fit_full_unconstrained[terms.fit.log_likelihood],6,3);
    rer.print_row [6] = Format (selection.io.getIC(rer.fit_full_unconstrained[terms.fit.log_likelihood], rer.fit_full_unconstrained[terms.parameters] + rer.empirical_parameters, rer.sample_size),6,3);
    fprintf(stdout, io.FormatTableRow(rer.print_row, rer.table_output_options));
}

// REPORT OUTLIER DETECTION RESULTS, IF DONE
// -----------------------------------------

if (rer.outlier) {
    io.ReportProgressMessageMD ("rer", "outlier", "Single-branch signal");
    io.ReportProgressMessageMD ("rer", "outlier", "> p-values corrected using the Holm-Bonferroni procedure");
    io.ReportProgressMessageMD ("rer", "outlier", "> (*) marks a possible outlier branch which is mostly responsible for rate differences");
    io.ReportProgressMessageMD ("rer", "outlier", "> (#) marks a branch which contributes to the rate difference signal");
    io.ReportProgressMessageMD ("rer", "outlier", "");


    rer.table_output_options = {
        utility.getGlobalValue("terms.table_options.header"): 1,
        utility.getGlobalValue("terms.table_options.minimum_column_width"): 16,
        utility.getGlobalValue("terms.table_options.align"): "left",
        utility.getGlobalValue("terms.table_options.column_widths"): {
            "0": 50,
            "1" :20,
            "2" :20
        }
    };
    
    rer.header = {
            3,
            1
        };
        
    rer.header[0] = "Branch";
    rer.header[1] = "Single branch";
    rer.header[2] = "All but this branch";
    
    fprintf(stdout, io.FormatTableRow(rer.header, rer.table_output_options));
    rer.table_output_options[utility.getGlobalValue("terms.table_options.header")] = FALSE;
    
    rer.print_row = {
            3,
            1
        };
             
    rer.print_row [0] = "";   
    
    rer.cutoff = 0.05;
    
    for (n,p; in; rer.by_test_branch) {
        rer.print_row[0] = n;
        rer.pn = (p[terms.Null])[terms.json.corrected_pvalue];
        rer.pa = (p[terms.alternative])[terms.json.corrected_pvalue];
        if (rer.pn <= rer.cutoff && rer.pa > rer.cutoff) {
            rer.print_row[0] = "(*) " + rer.print_row[0];
        }
        if (rer.pn <= rer.cutoff && rer.pa <= rer.cutoff) {
            rer.print_row[0] = "(#) " + rer.print_row[0];
        }

        rer.print_row[1] = Format ((p[terms.Null])[terms.json.corrected_pvalue], 10, 6);
        rer.print_row[2] = Format ((p[terms.alternative])[terms.json.corrected_pvalue], 10, 6);
        fprintf(stdout, io.FormatTableRow(rer.print_row, rer.table_output_options));
    }   
    
    rer.json [terms.rer.branch_level] = rer.by_test_branch;
}

// COMPUTE AND REPORT TEST STATISTICS
// ----------------------------------

io.ReportProgressMessageMD ("rer", "testing", "LRT test results between pairs of nested models");
io.ReportProgressMessageMD ("rer", "testing", "> p-values corrected using the Holm-Bonferroni procedure (raw values in parentheses)");
io.ReportProgressMessageMD ("rer", "testing", "");

if (rer.fit_full_model) {
    rer.hierarchy = {{terms.rer.proportional,terms.rer.proportional_partitoned,terms.rer.unconstrained_test,terms.rer.unconstrained}};
    rer.fit_results = {
        "0" : {{rer.fit_proportional[terms.fit.log_likelihood],rer.fit_proportional[terms.parameters] + rer.empirical_parameters}},
        "1" : {{rer.fit_proportional_partitioned[terms.fit.log_likelihood],rer.fit_proportional_partitioned[terms.parameters]}},
        "2" : {{rer.fit_test_unconstrained[terms.fit.log_likelihood],rer.fit_test_unconstrained[terms.parameters] + rer.empirical_parameters}},
        "3" : {{rer.fit_full_unconstrained[terms.fit.log_likelihood],rer.fit_full_unconstrained[terms.parameters] + rer.empirical_parameters}}
    };
    rer.model_count = 4;
} else {
    rer.hierarchy = {{terms.rer.proportional,terms.rer.proportional_partitoned,terms.rer.unconstrained_test}};
    rer.fit_results = {
        "0" : {{rer.fit_proportional[terms.fit.log_likelihood],rer.fit_proportional[terms.parameters] + rer.empirical_parameters}},
        "1" : {{rer.fit_proportional_partitioned[terms.fit.log_likelihood],rer.fit_proportional_partitioned[terms.parameters]}},
        "2" : {{rer.fit_test_unconstrained[terms.fit.log_likelihood],rer.fit_test_unconstrained[terms.parameters] + rer.empirical_parameters}}
    };
    rer.model_count = 3;
}


rer.test_matrix = {};
rer.test_lrt = {};


for (i = 0; i < rer.model_count ; i+=1) {
    for (j = i+1; j < rer.model_count ; j+=1) {
        rer.lrt = math.DoLRT ((rer.fit_results[i])[0],(rer.fit_results[j])[0],(rer.fit_results[j])[1]-(rer.fit_results[i])[1]);
        rer.test_matrix[rer.hierarchy[i] + ":" + rer.hierarchy[j]] = rer.lrt[terms.p_value];
        rer.test_lrt[rer.hierarchy[i] + ":" + rer.hierarchy[j]] = rer.lrt[terms.LRT]; 
    }
}   



rer.test_matrix_corrected = math.HolmBonferroniCorrection(rer.test_matrix);
rer.json [terms.json.test_results] = {};

for (k,v; in; rer.test_matrix) {
    (rer.json [terms.json.test_results])[k] = {
        terms.json.uncorrected_pvalue: v,
        terms.LRT: rer.test_lrt [k],
        terms.json.corrected_pvalue : rer.test_matrix_corrected[k]
    };
}


rer.table_output_options = {
        utility.getGlobalValue("terms.table_options.header"): 1,
        utility.getGlobalValue("terms.table_options.minimum_column_width"): 16,
        utility.getGlobalValue("terms.table_options.align"): "left",
        utility.getGlobalValue("terms.table_options.column_widths"): {
            "0": 40,
            "1" :40,
            "2" :40,
            "3" :40
        }
    };
    
rer.header = {
        4,
        1
    };
    
rer.header[0] = "Null";
rer.header[1] = terms.rer.proportional_partitoned;
rer.header[2] = terms.rer.unconstrained_test;
rer.header[3] = terms.rer.unconstrained;

fprintf(stdout, io.FormatTableRow(rer.header, rer.table_output_options));
rer.table_output_options[utility.getGlobalValue("terms.table_options.header")] = FALSE;

rer.print_row = {
            4,
            1
        };
             
rer.print_row [0] = "";   

for (i = 0; i < rer.model_count - 1; i+=1) {
    rer.print_row [0] = rer.hierarchy[i];
    for (j = 1; j <= i; j+=1) {
        rer.print_row [j] = "N/A";
    }
    for (j = i+1; j < rer.model_count ; j+=1) {
        rer.key = rer.hierarchy[i] + ":" + rer.hierarchy[j];
        rer.print_row [j] = Format (rer.test_matrix_corrected[rer.key], 8, 7) + " (" + Format (rer.test_matrix[rer.key], 6, 5) + ")";
    }
    
    fprintf(stdout, io.FormatTableRow(rer.print_row, rer.table_output_options));
}   



// BOOK-KEEPING
// -------------

selection.io.stopTimer (rer.json [terms.json.timers], "Overall");

io.ReportProgressMessageMD ("rer", "writing", "Writing detailed analysis report to \`" + rer.alignment [terms.json.json] + "\'");
io.SpoolJSON (rer.json, rer.alignment [terms.json.json]);

return 0;
                

// SUPPORTING FUNCTIONS
//---------------------

lfunction rer.array2dict (array) {
    dict = {};
    N = 0;
    for (a; in; array) {
        dict[a&&1] = 1;
        N += 1;
    }
    io.CheckAssertion("N == Abs (dict)", "RERConverge performs case-insensitive string matches; it seems that some of the tree labels cannot be resolved in a case-insensitive manner");
    return dict;
}

lfunction rer.intersect (dict1, dict2) {
    dict = {};
    
    for (k,d; in; dict1) {
        if (dict2 / k) {
            dict[k] = 1;
        }
    }
    return dict;
}

lfunction rer.set_branch_length_constraint (tree_name, node_name, model_description, ignore) {
    prms = model_description[^"terms.local"];
    assert (utility.Array1D (prms) == 1, "Node `node_name` was assigned a model with multiple local parameters" );
    prms = prms [terms.timeParameter()];
    assert (Abs (prms) > 0, "Node `node_name` was missing an evolutionary time parameter" );
    prms = "`tree_name`.`node_name`.`prms`";
    return prms;
}

lfunction rer.collect_constraints (tree_name, node_name, model_description, ignore) {
    prms = model_description[^"terms.local"];
    prms = prms [terms.timeParameter()];
    prms = "`tree_name`.`node_name`.`prms`";
    return prms;
}

lfunction rer.process_branch_lengths (results, labels) {
    lengths = {};
    
    for (b,d; in; (results[^"terms.branch_length"])[0]) {
        l = labels[b];
        utility.EnsureKey (lengths, l);
        (lengths[l])[b] = d[^"terms.fit.MLE"];
    }
    return lengths;
}

lfunction rer.branch_length_processor (lf_id, lf_components, data_filter, tree, model_map, initial_values, model_objects) {
    
    parameter_set = estimators.TraverseLocalParameters (lf_id, model_objects, "rer.set_branch_length_constraint");
    model = utility.GetFirstDictElement(model_objects);
                      
    for (node,p; in; parameter_set) {
        bl = (((initial_values[^"terms.branch_length"])["0"])[node])[^"terms.fit.MLE"];
        class_tag = (^"rer.branches")[node];
        io.CheckAssertion ("Abs(`&class_tag`)>0", "Unannotated node `node`");
        class = (^"rer.scalers")[class_tag];
        model.generic.AddGlobal (model, class, ^"terms.model.branch_length_scaler" + " for `class_tag`");
        parameters.SetConstraint (p,parameters.AppendMultiplicativeTerm (class,bl),"");
    }
    
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------
function rer.model.withGamma (options) {
	def = Call  (rer.model.generator.base, options);
	def [utility.getGlobalValue("terms.model.rate_variation")] = rate_variation.types.Gamma.factory (
	    {utility.getGlobalValue("terms.rate_variation.bins") : rer.rateClasses});
	return def;
};

//------------------------------------------------------------------------------------------------------------------------
function rer.model.withGDD  (options) {
	def = Call  (rer.model.generator.base, options);
	def [utility.getGlobalValue("terms.model.rate_variation")] = rate_variation.types.GDD.factory (
	    {utility.getGlobalValue("terms.rate_variation.bins") : rer.rateClasses});
	return def;
};