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
KeywordArgument ("model", "Model template abbreviation", "MG94", "Select a standard model");

// Initialize analysis description and json results
acd.analysis_description = {
    terms.io.info: "This analysis fits standard codon model templates to codon alignment data.",
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

namespace acd {
    #include "SelectionAnalyses/modules/shared-load-file.bf";
}

namespace acd {
    function load_file (prefix) {
        settings = None;
        multiple_files = FALSE;
        blb = 1.0;

        if (Type (prefix) == "AssociativeList") {
            multiple_files  = prefix [utility.getGlobalValue("terms.multiple_files")];
            settings        = prefix [utility.getGlobalValue("terms.settings")];
            blb             = prefix [utility.getGlobalValue("terms.data.blb_subsample")];
            prefix          = prefix [utility.getGlobalValue("terms.prefix")];
        }

        datasets = prefix+".codon_data";
        codon_data_info = alignments.PromptForGeneticCodeAndAlignment(datasets, prefix+".codon_filter");
        
        // Check use-tree option
        if (IS_TREE_PRESENT_IN_DATA) {
            use_tree = io.SelectAnOption ({"Yes": "Use tree in data file if present", "No": "Do not use tree in data file"}, "Use tree in data file if present");
            if (use_tree == "No") {
                IS_TREE_PRESENT_IN_DATA = 0;
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
}

acd_partition_tree_info = acd.partitions_and_trees["0"];
tree_info = acd_partition_tree_info[terms.data.tree];
abs_branch_lengths = Abs(tree_info[terms.branch_length]);

if (abs_branch_lengths == 0) {
    tree_string = tree_info[terms.trees.newick];
} else {
    tree_string = tree_info[terms.trees.newick_with_lengths];
}

// Define global variables for legacy models
_Genetic_Code = acd.codon_data_info[terms.code];
GeneticCodeExclusions = acd.codon_data_info[terms.stop_codons];

// Define filteredData for legacy templates
DataSetFilter filteredData = CreateFilter(acd.codon_data, 3, "", "", GeneticCodeExclusions);

SelectTemplateModel(filteredData);

// Branch lengths choice
KeywordArgument ("branch-lengths", "Branch lengths scaling", "Estimate", "Branch Lengths");
ChoiceList (branchLengths, "Branch Lengths", 1, SKIP_NONE,
                           "Estimate", "Estimate branch lengths by ML",
                           "Proportional to input tree", "Branch lengths are proportional to those in input tree");

if (branchLengths < 0) {
    return;
}

// Define the tree using the Newick string loaded from the library
treeString = tree_string;
Tree givenTree = treeString;

if (branchLengths == 1) {
    global treeScaler = 1;
    ReplicateConstraint ("this1.?.?:=treeScaler*this2.?.?__", givenTree, givenTree);
}

LikelihoodFunction lf = (filteredData, givenTree);
Optimize (res, lf);

// Print Markdown output
fprintf (stdout, "\n### Analysis Description\n\n");
fprintf (stdout, acd.analysis_description[terms.io.info], "\n\n");

fprintf (stdout, "### Model Fit Summary\n\n");
acd.table_options = {
    terms.table_options.header: TRUE,
    terms.table_options.minimum_column_width: 25,
    terms.table_options.align: "center"
};
fprintf (stdout, io.FormatTableRow({"0": "Statistic", "1": "Value"}, acd.table_options));
acd.table_options[terms.table_options.header] = FALSE;

log_l = Format(res[1][0], 8, 4);
params = Format(res[1][1], 0, 0);
n = acd.codon_data_info[terms.data.sample_size];
k = res[1][1];
aic = Format(2 * k - 2 * res[1][0], 8, 4);
if (n > k + 1) {
    aicc_val = 2 * k - 2 * res[1][0] + 2 * k * (k + 1) / (n - k - 1);
} else {
    aicc_val = 2 * k - 2 * res[1][0];
}
aicc = Format(aicc_val, 8, 4);
bic = Format(k * Log(n) - 2 * res[1][0], 8, 4);

fprintf (stdout, io.FormatTableRow({"0": "Log Likelihood", "1": log_l}, acd.table_options));
fprintf (stdout, io.FormatTableRow({"0": "Estimated Parameters", "1": params}, acd.table_options));
fprintf (stdout, io.FormatTableRow({"0": "AIC", "1": aic}, acd.table_options));
fprintf (stdout, io.FormatTableRow({"0": "AIC-c", "1": aicc}, acd.table_options));
fprintf (stdout, io.FormatTableRow({"0": "BIC", "1": bic}, acd.table_options));

GetString(lf_info, lf, -1);

fprintf (stdout, "\n### Global Parameter Estimates\n\n");
acd.global_table_options = {
    terms.table_options.header: TRUE,
    terms.table_options.minimum_column_width: 16,
    terms.table_options.align: "center"
};
fprintf (stdout, io.FormatTableRow({"0": "Parameter", "1": "Estimate"}, acd.global_table_options));
acd.global_table_options[terms.table_options.header] = FALSE;

globals = lf_info["Global Independent"];
global_params = {};
for (k = 0; k < Columns(globals); k = k + 1) {
    global_params[globals[k]] = Eval(globals[k]);
    val = Format(Eval(globals[k]), 8, 6);
    fprintf (stdout, io.FormatTableRow({"0": globals[k], "1": val}, acd.global_table_options));
}

fprintf (stdout, "\n### Branch Lengths\n\n");
acd.branch_table_options = {
    terms.table_options.header: TRUE,
    terms.table_options.minimum_column_width: 16,
    terms.table_options.align: "center"
};
fprintf (stdout, io.FormatTableRow({"0": "Branch", "1": "Length"}, acd.branch_table_options));
acd.branch_table_options[terms.table_options.header] = FALSE;

branch_names = BranchName(givenTree, -1);
for (k = 0; k < Columns(branch_names) - 1; k = k + 1) {
    val = Format(BranchLength(givenTree, branch_names[k]), 8, 6);
    fprintf (stdout, io.FormatTableRow({"0": branch_names[k], "1": val}, acd.branch_table_options));
}

local_params = {};
locals = lf_info["Local Independent"];
for (k = 0; k < Columns(locals); k = k + 1) {
    local_params[locals[k]] = Eval(locals[k]);
}

// Build final JSON
acd.json = {
    terms.json.analysis: acd.analysis_description,
    terms.json.input: {
        terms.data.file: acd.codon_data_info[terms.data.file],
        terms.data.sequences: acd.codon_data_info[terms.data.sequences],
        terms.data.sites: acd.codon_data_info[terms.data.sites],
        terms.json.genetic_code: acd.codon_data_info[terms.id]
    },
    terms.json.fits: {
        terms.json.model: {
            terms.json.log_likelihood: res[1][0],
            terms.json.parameters: res[1][1],
            terms.json.AIC: 2 * res[1][1] - 2 * res[1][0],
            terms.json.AICc: aicc_val,
            terms.json.BIC: res[1][1] * Log(acd.codon_data_info[terms.data.sample_size]) - 2 * res[1][0],
            terms.json.global_parameters: global_params,
            terms.json.local_parameters: local_params
        }
    },
    terms.json.trees: {
        "0": Format(givenTree, 1, 1)
    }
};

KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + '.json')", acd.codon_data_info [terms.json.json], "Save the resulting JSON to this file");
output_path = io.PromptUserForFilePath("Save the resulting JSON to this file");

io.SpoolJSON (acd.json, output_path);

GetString (sendMeBack, lf, -1);
sendMeBack["LogL"] = res[1][0];
sendMeBack["NP"]   = res[1][1];
return sendMeBack;
