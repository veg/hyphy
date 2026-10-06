RequireVersion("2.5.99");

LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");

KeywordArgument ("code", "Which genetic code should be used", "Universal", "Choose Genetic Code");
KeywordArgument ("alignment", "An in-frame codon alignment in one of the formats supported by HyPhy", null, "Select a coding sequence alignment file");
KeywordArgument ("use-tree", "Use tree in data file if present", "Yes", "Use tree in data file if present");
KeywordArgument ("tree", "A phylogenetic tree (optionally annotated with {})", null, "Please select a tree file for the data:");

// Initialize analysis description and json results
acd.analysis_description = {
    terms.io.info: "This analysis fits a wide range of codon substitution models to a codon alignment. Supported models include mechanistic options (MG94, GY94) crossed with arbitrary nucleotide substitution models, empirical codon models (ECM, MEC), acceptance bias models (MG94wEX), and property-based models (LCAP). It supports multiple frequency estimators (F1x4, F3x4, CF3x4, F61, positional 9-frequency), branch length options (independent estimation or proportional to input trees), and parameterizations (local or global rates). Rate variation across sites can be modeled using a variety of continuous and discrete distributions (including Gamma, Beta, Log-normal, and General Discrete) optionally coupled with spatial autocorrelation via a Hidden Markov Model (HMM).",
    terms.io.version: "2.0",
    terms.io.authors: "Sergei L Kosakovsky Pond",
    terms.io.contact: "spond@temple.edu",
    terms.io.requirements: "in-frame codon alignment and a phylogenetic tree."
};
io.DisplayAnalysisBanner(acd.analysis_description);

terms.json.BIC = "BIC";
terms.json.AIC = "AIC";
terms.json.global_parameters = "global-parameters";
terms.json.local_parameters = "local-parameters";
terms.json.genetic_code = "genetic-code";

acd.json = {
    terms.json.analysis: acd.analysis_description,
    terms.json.fits: {},
    terms.json.timers: {}
};
json = acd.json;

_Genetic_Code = None;
GeneticCodeExclusions = None;

namespace acd {
    #include "SelectionAnalyses/modules/shared-load-file.bf";
}

namespace acd {
    function load_file (prefix) {
        datasets = prefix+".codon_data";
        codon_data_info = alignments.PromptForGeneticCodeAndAlignment(datasets, prefix+".codon_filter");
        
        // Check use-tree option
        if (utility.GetEnvVariable("IS_TREE_PRESENT_IN_DATA")) {
            use_tree = io.SelectAnOption ({"Yes": "Use tree in data file if present", "No": "Do not use tree in data file"}, "Use tree in data file if present");
            if (use_tree == "No") {
                utility.SetEnvVariable("IS_TREE_PRESENT_IN_DATA", FALSE);
            }
        }
        
        annotate_codon_info (codon_data_info, prefix+".codon_filter");
        utility.SetEnvVariable(utility.getGlobalValue ("terms.trees.data_for_neighbor_joining"),
                               codon_data_info[utility.getGlobalValue("terms.data.datafilter")]);
        partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (codon_data_info[utility.getGlobalValue("terms.data.partitions")], name_mapping);
        
        for (_key_, _value_; in; partitions_and_trees) {
            (partitions_and_trees[_key_])[utility.getGlobalValue("terms.data.filter_string")] = 
                    selection.io.adjust_partition_string (_value_[utility.getGlobalValue("terms.data.filter_string")], 3*codon_data_info[utility.getGlobalValue("terms.data.sites")]);
        }
        
        utility.SetEnvVariable(utility.getGlobalValue ("terms.trees.data_for_neighbor_joining"), None);
        partition_count = Abs (partitions_and_trees);
        
        // Branch selector logic is ignored for this analysis
        selected_branches = {};
        selection.io.json_store_key_value_pair (json, None, utility.getGlobalValue("terms.json.tested"), selected_branches);
        
        filter_specification = alignments.DefineFiltersForPartitions (partitions_and_trees, datasets , "`prefix`.filter.", codon_data_info);
        store_tree_information ();
    }
}

namespace acd {
    load_file ("acd");

    partition_tree_info = partitions_and_trees["0"];
    tree_info = partition_tree_info[^"terms.data.tree"];
    abs_branch_lengths = Abs(tree_info[^"terms.branch_length"]);

    if (abs_branch_lengths == 0) {
        tree_string = tree_info[^"terms.trees.newick"];
    } else {
        tree_string = tree_info[^"terms.trees.newick_with_lengths"];
    }

    // Define global variables for legacy models
    ^"_Genetic_Code" = codon_data_info[^"terms.code"];
    ^"GeneticCodeExclusions" = codon_data_info[^"terms.stop_codons"];

    // Branch lengths choice
    KeywordArgument ("branch-lengths", "Branch lengths scaling", "Estimate", "Branch Lengths");
    branchLengthsOption = io.SelectAnOption ({
        "Estimate": "Estimate branch lengths by ML",
        "Proportional to input tree": "Branch lengths are proportional to those in input tree"
    }, "Branch Lengths");
    // Define filteredData in the global namespace so template model works correctly
    DataSetFilter ^"filteredData" = CreateFilter(codon_data, 3, "", "", ^"GeneticCodeExclusions");
    KeywordArgument ("model", "Model template abbreviation", "MG94");
    SelectTemplateModel(^"filteredData");

    Tree givenTree = tree_string;

    if (branchLengthsOption == "Proportional to input tree") {
        global treeScaler = 1;
        ReplicateConstraint ("this1.?.?:=treeScaler*this2.?.?__", givenTree, givenTree);
    }

    vectorOfFrequencies = ^"vectorOfFrequencies";
    MG94 = ^"MG94";

    LikelihoodFunction lf = (^"filteredData", givenTree);
    Optimize (res, lf);

    // Print Markdown output
    fprintf (stdout, "\n### Analysis Description\n\n");
    fprintf (stdout, analysis_description[^"terms.io.info"], "\n\n");

    fprintf (stdout, "### Model Fit Summary\n\n");
    table_options = {
        ^"terms.table_options.header": TRUE,
        ^"terms.table_options.minimum_column_width": 25,
        ^"terms.table_options.align": "center"
    };
    fprintf (stdout, io.FormatTableRow({"0": "Statistic", "1": "Value"}, table_options));
    table_options[^"terms.table_options.header"] = FALSE;

    log_l = Format(res[1][0], 8, 4);
    params = Format(res[1][1], 0, 0);
    n = codon_data_info[^"terms.data.sample_size"];
    k = res[1][1];
    aic = Format(2 * k - 2 * res[1][0], 8, 4);
    if (n > k + 1) {
        aicc_val = 2 * k - 2 * res[1][0] + 2 * k * (k + 1) / (n - k - 1);
    } else {
        aicc_val = 2 * k - 2 * res[1][0];
    }
    aicc = Format(aicc_val, 8, 4);
    bic = Format(k * Log(n) - 2 * res[1][0], 8, 4);

    fprintf (stdout, io.FormatTableRow({"0": "Log Likelihood", "1": log_l}, table_options));
    fprintf (stdout, io.FormatTableRow({"0": "Estimated Parameters", "1": params}, table_options));
    fprintf (stdout, io.FormatTableRow({"0": "AIC", "1": aic}, table_options));
    fprintf (stdout, io.FormatTableRow({"0": "AIC-c", "1": aicc}, table_options));
    fprintf (stdout, io.FormatTableRow({"0": "BIC", "1": bic}, table_options));

    GetString(lf_info, lf, -1);

    fprintf (stdout, "\n### Global Parameter Estimates\n\n");
    global_table_options = {
        ^"terms.table_options.header": TRUE,
        ^"terms.table_options.minimum_column_width": 16,
        ^"terms.table_options.align": "center"
    };
    fprintf (stdout, io.FormatTableRow({"0": "Parameter", "1": "Estimate"}, global_table_options));
    global_table_options[^"terms.table_options.header"] = FALSE;

    globals = lf_info["Global Independent"];
    global_params = {};
    for (k = 0; k < Columns(globals); k = k + 1) {
        global_params[globals[k]] = Eval(globals[k]);
        val = Format(Eval(globals[k]), 8, 6);
        fprintf (stdout, io.FormatTableRow({"0": globals[k], "1": val}, global_table_options));
    }

    fprintf (stdout, "\n### Branch Lengths\n\n");
    branch_table_options = {
        ^"terms.table_options.header": TRUE,
        ^"terms.table_options.minimum_column_width": 16,
        ^"terms.table_options.align": "center"
    };
    fprintf (stdout, io.FormatTableRow({"0": "Branch", "1": "Length"}, branch_table_options));
    branch_table_options[^"terms.table_options.header"] = FALSE;

    branch_names = BranchName(givenTree, -1);
    for (k = 0; k < Columns(branch_names) - 1; k = k + 1) {
        val = Format(BranchLength(givenTree, branch_names[k]), 8, 6);
        fprintf (stdout, io.FormatTableRow({"0": branch_names[k], "1": val}, branch_table_options));
    }

    local_params = {};
    locals = lf_info["Local Independent"];
    for (k = 0; k < Columns(locals); k = k + 1) {
        local_params[locals[k]] = Eval(locals[k]);
    }

    // Build final JSON
    json = {
        ^"terms.json.analysis": analysis_description,
        ^"terms.json.input": {
            ^"terms.data.file": codon_data_info[^"terms.data.file"],
            ^"terms.data.sequences": codon_data_info[^"terms.data.sequences"],
            ^"terms.data.sites": codon_data_info[^"terms.data.sites"],
            ^"terms.json.genetic_code": codon_data_info[^"terms.id"]
        },
        ^"terms.json.fits": {
            ^"terms.json.model": {
                ^"terms.json.log_likelihood": res[1][0],
                ^"terms.json.parameters": res[1][1],
                ^"terms.json.AIC": 2 * res[1][1] - 2 * res[1][0],
                ^"terms.json.AICc": aicc_val,
                ^"terms.json.BIC": res[1][1] * Log(codon_data_info[^"terms.data.sample_size"]) - 2 * res[1][0],
                ^"terms.json.global_parameters": global_params,
                ^"terms.json.local_parameters": local_params
            }
        },
        ^"terms.json.trees": {
            "0": Format(givenTree, 1, 1)
        }
    };

    KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + '.ACD.json')", codon_data_info[^"terms.data.file"] + ".ACD.json", "Save the resulting JSON to this file");
    output_path = io.PromptUserForFilePath("Save the resulting JSON to this file");

    io.SpoolJSON (json, output_path);

    GetString (sendMeBack_val, lf, -1);
    sendMeBack_val["LogL"] = res[1][0];
    sendMeBack_val["NP"]   = res[1][1];
    sendMeBack = sendMeBack_val;
}

return acd.sendMeBack;
