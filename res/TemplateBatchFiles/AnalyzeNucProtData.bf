RequireVersion("2.5.99");

LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");
LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");

DeleteObject (inner_run);
DeleteObject (alignment);
DeleteObject (model);
DeleteObject (model_options);
DeleteObject (rate_distribution);
DeleteObject (rate_classes);
DeleteObject (tree);
DeleteObject (display_options);
DeleteObject (output);

ad.inner_run_val = utility.GetEnvVariable("inner_run");
if (Type(ad.inner_run_val) != "String") {
    ad.inner_run_val = "No";
}

// Retrieve keyword arguments list from the process if available
GetString (ad.kwargs_list, KWARGS, 0);
ad.alignment_path = null;
ad.tree_path = null;

if (Type(ad.kwargs_list) == "AssociativeList") {
    if (Type(ad.kwargs_list["alignment"]) == "String") { ad.alignment_path = ad.kwargs_list["alignment"]; }
    if (Type(ad.kwargs_list["tree"]) == "String") { ad.tree_path = ad.kwargs_list["tree"]; }
}

// Set GLOBAL_KWARGS for libv3 in the child run
if (ad.inner_run_val == "Yes") {
    utility.SetEnvVariable("GLOBAL_KWARGS", ad.kwargs_list);
}

// Translate legacy numeric keys in inputOptions to "--keyword" style to trigger user_kwargs in ExecuteAFile
if (Type(inputOptions) == "AssociativeList") {
    if (ad.inner_run_val == "No") {
        mapped_options = {};
        if (Type(inputOptions["0"]) == "String") { mapped_options["--alignment"] = inputOptions["0"]; }
        if (Type(inputOptions["1"]) == "String") { mapped_options["--model"] = inputOptions["1"]; }
        if (Type(inputOptions["2"]) == "String") { mapped_options["--model-options"] = inputOptions["2"]; }
        if (Type(inputOptions["3"]) == "String") { mapped_options["--rate-distribution"] = inputOptions["3"]; }
        if (Type(inputOptions["4"]) == "String") { mapped_options["--rate-classes"] = inputOptions["4"]; }
        if (Type(inputOptions["6"]) == "String") { mapped_options["--tree"] = inputOptions["6"]; }
        if (Type(inputOptions["7"]) == "String") { mapped_options["--display-options"] = inputOptions["7"]; }
        
        // Pass a dummy non-dash key to force redirection isolation in ExecuteAFile
        mapped_options["dummy_redirect"] = "isolate";
        
        path_to_self = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "AnalyzeNucProtData.bf";
        utility.SetEnvVariable("inner_run", "Yes");
        ExecuteAFile (path_to_self, mapped_options);
        utility.SetEnvVariable("inner_run", "No");
        inner_res = utility.GetEnvVariable("inner_res");
        res = inner_res["res"];
        lf = 0;
        return inner_res;
    }
}

ad.analysis_description = {
    terms.io.info: "This analysis fits standard nucleotide or amino-acid model templates to alignment data.",
    terms.io.version: "2.0",
    terms.io.authors: "Sergei L Kosakovsky Pond",
    terms.io.contact: "spond@temple.edu",
    terms.io.requirements: "nucleotide or amino-acid alignment and a phylogenetic tree."
};
io.DisplayAnalysisBanner(ad.analysis_description);

function ad.resolve_path(path) {
    if (Type(path) == "String") {
        if (path[0] != "/" && path[0] != DIRECTORY_SEPARATOR) {
            if (Abs(path) > 1 && path[1] != ":") {
                base_dir = utility.GetEnvVariable ("HYPHY_BASE_DIRECTORY");
                if (Type(base_dir) == "String") {
                    return base_dir + path;
                }
            }
        }
    }
    return path;
}

// Step 1: Load alignment
_Genetic_Code = 0;
ExecuteAFile("simpleBootstrap.bf");

if (Type(ad.alignment_path) == "String") {
    ad.alignment_path = ad.resolve_path(ad.alignment_path);
    DataSet ds = ReadDataFile (ad.alignment_path);
    DataSetFilter filteredData = CreateFilter (ds, 1);
    ad.alignment_path = ad.alignment_path;
} else {
    KeywordArgument ("alignment", "A nucleotide or amino-acid multiple sequence alignment file", null, "Please specify a nucleotide or amino-acid data file:");
    SetDialogPrompt ("Please specify a nucleotide or amino-acid data file:");
    DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
    DataSetFilter filteredData = CreateFilter (ds, 1);
    ad.alignment_path = LAST_FILE_PATH;
}

// Step 2: Select template model (Prompts 1, 2, 3, 4, 5)
KeywordArgument ("model", "Model template abbreviation", "HKY85", "Select a standard model");
SelectTemplateModel(filteredData);
_DO_TREE_REBALANCE_ = 1;

// Step 3: Load tree (Prompt 6)
if (Type(ad.tree_path) == "String") {
    ad.tree_path = ad.resolve_path(ad.tree_path);
    fscanf (ad.tree_path, "Raw", treeString);
    if (_DO_TREE_REBALANCE_) {
        treeString = RerootTree (treeString, 0);
    }
    Tree givenTree = treeString;
} else {
    KeywordArgument ("tree", "A phylogenetic tree", null, "Please select a tree file for the data:");
    SetDialogPrompt ("Please select a tree file for the data:");
    fscanf (PROMPT_FOR_FILE, REWIND, "Raw", treeString);
    if (_DO_TREE_REBALANCE_) {
        treeString = RerootTree (treeString, 0);
    }
    Tree givenTree = treeString;
}

LikelihoodFunction lf = (filteredData, givenTree);
timer = Time(0);
Optimize (res, lf);
timer = Time(0) - timer;

// Step 4: Execute categoryEcho.bf for rate variation display (Prompt 7)
KeywordArgument ("display-options", "Rate class Information Dislpay", "Don't Display", "Rate class Information Dislpay");
ExecuteAFile ("categoryEcho.bf");
GetString	 (lfInfo, lf, -1);

// Print Markdown output
fprintf (stdout, "\n### Analysis Description\n\n");
fprintf (stdout, ad.analysis_description[terms.io.info], "\n\n");

fprintf (stdout, "### Model Fit Summary\n\n");
ad.table_options = {
    terms.table_options.header: TRUE,
    terms.table_options.minimum_column_width: 25,
    terms.table_options.align: "center"
};
fprintf (stdout, io.FormatTableRow({"0": "Statistic", "1": "Value"}, ad.table_options));
ad.table_options[terms.table_options.header] = FALSE;

log_l = Format(res[1][0], 8, 4);
params = Format(res[1][1], 0, 0);
n = ds.species * ds.sites;
k = res[1][1];
aic = Format(2 * k - 2 * res[1][0], 8, 4);
if (n > k + 1) {
    aicc_val = 2 * k - 2 * res[1][0] + 2 * k * (k + 1) / (n - k - 1);
} else {
    aicc_val = 2 * k - 2 * res[1][0];
}
aicc = Format(aicc_val, 8, 4);
bic = Format(k * Log(n) - 2 * res[1][0], 8, 4);

fprintf (stdout, io.FormatTableRow({"0": "Log Likelihood", "1": log_l}, ad.table_options));
fprintf (stdout, io.FormatTableRow({"0": "Estimated Parameters", "1": params}, ad.table_options));
fprintf (stdout, io.FormatTableRow({"0": "AIC", "1": aic}, ad.table_options));
fprintf (stdout, io.FormatTableRow({"0": "AIC-c", "1": aicc}, ad.table_options));
fprintf (stdout, io.FormatTableRow({"0": "BIC", "1": bic}, ad.table_options));

fprintf (stdout, "\n### Global Parameter Estimates\n\n");
ad.global_table_options = {
    terms.table_options.header: TRUE,
    terms.table_options.minimum_column_width: 16,
    terms.table_options.align: "center"
};
fprintf (stdout, io.FormatTableRow({"0": "Parameter", "1": "Estimate"}, ad.global_table_options));
ad.global_table_options[terms.table_options.header] = FALSE;

globals = lfInfo["Global Independent"];
global_params = {};
for (k = 0; k < Columns(globals); k = k + 1) {
    global_params[globals[k]] = Eval(globals[k]);
    val = Format(Eval(globals[k]), 8, 6);
    fprintf (stdout, io.FormatTableRow({"0": globals[k], "1": val}, ad.global_table_options));
}

fprintf (stdout, "\n### Branch Lengths\n\n");
ad.branch_table_options = {
    terms.table_options.header: TRUE,
    terms.table_options.minimum_column_width: 16,
    terms.table_options.align: "center"
};
fprintf (stdout, io.FormatTableRow({"0": "Branch", "1": "Length"}, ad.branch_table_options));
ad.branch_table_options[terms.table_options.header] = FALSE;

branch_names = BranchName(givenTree, -1);
for (k = 0; k < Columns(branch_names) - 1; k = k + 1) {
    val = Format(BranchLength(givenTree, branch_names[k]), 8, 6);
    fprintf (stdout, io.FormatTableRow({"0": branch_names[k], "1": val}, ad.branch_table_options));
}

local_params = {};
locals = lfInfo["Local Independent"];
for (k = 0; k < Columns(locals); k = k + 1) {
    local_params[locals[k]] = Eval(locals[k]);
}

// Build final JSON
ad.json = {
    terms.json.analysis: ad.analysis_description,
    terms.json.input: {
        terms.data.file: ad.alignment_path,
        terms.data.sequences: ds.species,
        terms.data.sites: ds.sites
    },
    terms.json.fits: {
        terms.json.model: {
            terms.json.log_likelihood: res[1][0],
            terms.json.parameters: res[1][1],
            terms.json.AIC: 2 * res[1][1] - 2 * res[1][0],
            terms.json.AICc: aicc_val,
            terms.json.BIC: res[1][1] * Log(n) - 2 * res[1][0],
            terms.json.global_parameters: global_params,
            terms.json.local_parameters: local_params
        }
    },
    terms.json.trees: {
        "0": Format(givenTree, 1, 1)
    }
};

ad.default_json_path = ad.alignment_path + ".AD.json";

ad.global_kwargs = utility.GetEnvVariable("GLOBAL_KWARGS");
output_path = null;
if (Type(ad.global_kwargs) == "AssociativeList") {
    output_path = utility.GetByKey (ad.global_kwargs, "output", "String");
}
if (null == output_path) {
    if (Type(inputOptions) == "AssociativeList") {
        output_path = ad.default_json_path;
    } else {
        KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + '.AD.json')", ad.default_json_path, "Save the resulting JSON to this file");
        output_path = io.PromptUserForFilePath("Save the resulting JSON to this file");
    }
}

output_path = ad.resolve_path(output_path);
io.SpoolJSON (ad.json, output_path);

utility.SetEnvVariable("GLOBAL_KWARGS", None);
utility.SetEnvVariable("inner_res", {"Log(L)": res[1][0], "DF": res[1][1], "Tree": Format (givenTree,1,1), "res": res});
return {"Log(L)": res[1][0], "DF": res[1][1], "Tree": Format (givenTree,1,1)};
