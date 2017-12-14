RequireVersion("2.3");

/*------------------------------------------------------------------------------
    Load library files
*/

LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");
LoadFunctionLibrary("libv3/all-terms.bf");

LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");

LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");

LoadFunctionLibrary("libv3/convenience/math.bf");

LoadFunctionLibrary("modules/io_functions.ibf");

/* Display analysis information */

prime.analysis_description = {
    terms.io.info: "PRIME",
    terms.io.version: "1.00",
    terms.io.reference: "TBD",
    terms.io.authors: "Sergei L Kosakovsky Pond",
    terms.io.contact: "spond@temple.edu",
    terms.io.requirements: "in-frame codon alignment and a phylogenetic tree"
};

io.DisplayAnalysisBanner(prime.analysis_description);

/* Environment Setup */
utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

/* Globals */
prime.site_alpha = "Site relative synonymous rate";
prime.site_beta = "Site relative non-synonymous rate (tested branches)";
prime.site_beta_nuisance = "Site relative non-synonymous rate (untested branches)";

// default cutoff for printing to screen
prime.p_value = 0.1;
prime.scaler_prefix = "PRIME.scaler";

// The dictionary of results to be written to JSON at the end of the run
prime.json = {
    terms.json.analysis: prime.analysis_description,
    terms.json.input: {},
    terms.json.fits: {},
    terms.json.timers: {},
};

selection.io.startTimer (prime.json [terms.json.timers], "Total time", 0);

prime.table_headers = {
                         {"alpha_0", "Synonymous substitution rate at a site"}
                         {"alpha_0_pval", "The rate estimate under the neutral model"}
                         {"alpha_1", "Synonymous substitution rate at a site"}
                         {"alpha_1_pval", "The rate estimate under the neutral model"}
                         {"alpha_2", "Synonymous substitution rate at a site"}
                         {"alpha_2_pval", "The rate estimate under the neutral model"}
                         {"alpha_3", "Synonymous substitution rate at a site"}
                         {"alpha_3_pval", "The rate estimate under the neutral model"}
                         {"alpha_4", "Synonymous substitution rate at a site"}
                         {"alpha_4_pval", "The rate estimate under the neutral model"}
                     };


/**
This table is meant for HTML rendering in the results web-app; can use HTML characters, the second column
is 'pop-over' explanation of terms. This is ONLY saved to the JSON file. For Markdown screen output see
the next set of variables.
*/
prime.table_screen_output  = {{"Codon", "Partition", "alpha", "beta", "LRT", "Selection detected?"}};
prime.table_output_options = {terms.table_options.header : TRUE, terms.table_options.minimum_column_width: 16, terms.table_options.align : "center"};

// Load all the file information
namespace prime {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ("prime");
}


selection.io.startTimer (prime.json [terms.json.timers], "Model fitting",1);

namespace prime {
    doGTR ("prime");
}

estimators.fixSubsetOfEstimates(prime.gtr_results, prime.gtr_results[terms.global]);

io.ReportProgressMessageMD ("prime", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");

namespace prime {
    doPartitionedMG ("prime", FALSE);
}

prime.final_partitioned_mg_results = estimators.FitMGREV(prime.filter_names, prime.trees, prime.codon_data_info [terms.code], {
    terms.run_options.model_type: terms.local,
    terms.run_options.partitioned_omega: prime.selected_branches,
    terms.run_options.retain_lf_object: TRUE
}, prime.partitioned_mg_results);

io.ReportProgressMessageMD("prime", "codon-refit", "* Log(L) = " + Format(prime.final_partitioned_mg_results[terms.fit.log_likelihood],8,2));

prime.global_dnds = selection.io.extract_global_MLE_re (prime.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio);
utility.ForEach (prime.global_dnds, "_value_", 'io.ReportProgressMessageMD ("PRIME", "codon-refit", "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));');

estimators.fixSubsetOfEstimates(prime.final_partitioned_mg_results, prime.final_partitioned_mg_results[terms.global]);

//Store MG94 to JSON
selection.io.json_store_lf_GTR_MG94 (prime.json,
                            terms.json.global_mg94xrev,
                            prime.final_partitioned_mg_results[terms.fit.log_likelihood],
                            prime.final_partitioned_mg_results[terms.parameters],
                            prime.sample_size,
                            utility.ArrayToDict (utility.Map (prime.global_dnds, "_value_", "{'key': _value_[terms.description], 'value' : Eval({{_value_ [terms.fit.MLE],1}})}")),
                            (prime.final_partitioned_mg_results[terms.efv_estimate])["VALUEINDEXORDER"][0],
                            prime.display_orders[terms.json.global_mg94xrev]);

utility.ForEachPair (prime.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(prime.json, terms.json.global_mg94xrev, terms.branch_length, prime.display_orders[terms.json.global_mg94xrev],
                                             _key_,
                                             selection.io.extract_branch_info((prime.final_partitioned_mg_results[terms.branch_length])[_key_], "selection.io.branch.length"));');




selection.io.stopTimer (prime.json [terms.json.timers], "Model fitting");

// define the site-level likelihood function

// TODO: This is where we need to implement the QCAP model

prime.site.mg_rev = model.generic.DefineModel("models.codon.MG_REV.ModelDescription",
        "prime_mg", {
            "0": parameters.Quote(terms.local),
            "1": prime.codon_data_info[terms.code]
        },
        prime.filter_names,
        None);


prime.site_model_mapping = {"prime_mg" : prime.site.mg_rev};

/* set up the local constraint model */

// TODO : will need to update
prime.alpha = model.generic.GetLocalParameter (prime.site.mg_rev, utility.getGlobalValue("terms.parameters.synonymous_rate"));
prime.beta = model.generic.GetLocalParameter (prime.site.mg_rev, utility.getGlobalValue("terms.parameters.nonsynonymous_rate"));
io.CheckAssertion ("None!=prime.alpha && None!=prime.beta", "Could not find expected local synonymous and non-synonymous rate parameters in \`estimators.FitMGREV\`");

selection.io.startTimer (prime.json [terms.json.timers], "prime analysis", 2);

//----------------------------------------------------------------------------------------
function fel.apply_proportional_site_constraint (tree_name, node_name, alpha_parameter, beta_parameter, alpha_factor, beta_factor, branch_length) {

    fel.branch_length = (branch_length[terms.parameters.synonymous_rate])[terms.fit.MLE];

    node_name = tree_name + "." + node_name;

    ExecuteCommands ("
        `node_name`.`alpha_parameter` := (`alpha_factor`) * fel.branch_length__;
        `node_name`.`beta_parameter`  := (`beta_factor`)  * fel.branch_length__;
    ");
}





