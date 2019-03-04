RequireVersion("2.3.12");


/*------------------------------------------------------------------------------
    Load library files
*/

LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");

LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");

LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");

LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary("SelectionAnalyses/modules/selection_lib.ibf");


/*------------------------------------------------------------------------------
    Display analysis information
*/

single_omega.analysis_description = {
    terms.io.info: "This analysis applies the Muse Gaut 94 (MG94) model with the GTR nucleotide rate component 
    to estimate the single dN/dS ratio for the alignment, and report its confidence intervals derived using profile likelihood.
    There is an option to select a subset of branches for estimation, in which case a separate omega is estimated for the selected branches.
    ",
    terms.io.version: "0.1",
    terms.io.reference: "A likelihood approach for comparing synonymous and nonsynonymous nucleotide substitution rates, with application to the chloroplast genome (1994). _Mol Biol Evol_ 11 (5): 715-724",
    terms.io.authors: "Sergei L Kosakovsky Pond",
    terms.io.contact: "spond@temple.edu",
    terms.io.requirements: "in-frame codon alignment and a phylogenetic tree. Multiple partitions with a NEXUS file also supported"
};
io.DisplayAnalysisBanner(single_omega.analysis_description);


/*------------------------------------------------------------------------------
    Environment setup
*/

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

/*------------------------------------------------------------------------------
    Globals
*/


single_omega.json = {
    terms.json.analysis: single_omega.analysis_description,
    terms.json.fits: {},
    terms.json.timers: {}
};
    /**
        The dictionary of results to be written to JSON at the end of the run
    */

single_omega.display_orders =   {terms.original_name: -1,
                        terms.json.nucleotide_gtr: 0,
                        terms.json.global_mg94xrev: 1
                       };


selection.io.startTimer (single_omega.json [terms.json.timers], "Total time", 0);

namespace single_omega {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ("single_omega");
}

namespace single_omega {
    doGTR ("single_omega");
}


single_omega.scaler_prefix = "single_omega.scaler";
estimators.fixSubsetOfEstimates(single_omega.gtr_results, single_omega.gtr_results[terms.global]);

namespace single_omega {
    doPartitionedMG ("single_omega", TRUE);
}

io.ReportProgressMessageMD ("single_omega", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");


single_omega.partitioned_mg_results = estimators.FitMGREV (single_omega.filter_names, single_omega.trees, single_omega.codon_data_info [terms.code], {
    terms.run_options.model_type: terms.local,
    terms.run_options.partitioned_omega: single_omega.selected_branches,
    terms.run_options.retain_lf_object: TRUE
}, single_omega.partitioned_mg_results);


io.ReportProgressMessageMD("single_omega", "codon-refit", "* Log(L) = " + Format(single_omega.partitioned_mg_results[terms.fit.log_likelihood],8,2));

io.ReportProgressMessageMD ("single_omega", "codon-refit", "* Estimating confidence intervals");

single_omega.global_dnds = selection.io.extract_global_MLE_CI_re (single_omega.partitioned_mg_results, "^" + terms.parameters.omega_ratio);


utility.ForEachPair (single_omega.global_dnds, "_tag_", "_value_", '
    io.ReportProgressMessageMD ("single_omega", "codon-refit", "* " + _tag_ + " = " + Format (_value_[terms.fit.MLE],8,4) + " [" 
                                                                                    + Format (_value_[terms.lower_bound],8,4)
                                                                                    + " - " 
                                                                                    + Format (_value_[terms.upper_bound],8,4) 
                                                                                    + " ]");
');


selection.io.json_store_lf_withEFV (single_omega.json,
                            terms.json.global_mg94xrev,
                            single_omega.partitioned_mg_results[terms.fit.log_likelihood],
                            single_omega.partitioned_mg_results[terms.parameters],
                            single_omega.sample_size,
                            utility.ArrayToDict (utility.Map (single_omega.global_dnds, "_value_", "{'key': _value_[terms.description], 'value' : Eval({{_value_ [terms.fit.MLE],1}})}")),
                            (single_omega.partitioned_mg_results[terms.efv_estimate])["VALUEINDEXORDER"][0],
                            single_omega.display_orders[terms.json.global_mg94xrev]);

utility.ForEachPair (single_omega.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(single_omega.json, terms.json.global_mg94xrev, terms.branch_length, single_omega.display_orders[terms.json.global_mg94xrev],
                                             _key_,
                                             selection.io.extract_branch_info((single_omega.partitioned_mg_results[terms.branch_length])[_key_], "selection.io.branch.length"));');


selection.io.json_store_key_value_pair (single_omega.json, terms.json.fits , terms.json.omega_ratio, single_omega.global_dnds);

selection.io.stopTimer (single_omega.json [terms.json.timers], "Total time");
io.SpoolJSON (single_omega.json, single_omega.codon_data_info[terms.json.json]);

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
lfunction selection.io.extract_global_MLE_CI_re (fit, regexp) {
     matches = utility.Filter (utility.Keys(fit[utility.getGlobalValue("terms.global")]), "_tag_", "regexp.Find (_tag_, `&regexp`)");
     matches = utility.Map (matches, "_tag_", "
        {utility.getGlobalValue('terms.description'): _tag_, 'ID' : ((`&fit`[utility.getGlobalValue('terms.global')])[_tag_])[utility.getGlobalValue('terms.id')]}
        "
    );
    
    
    result = {};
    utility.ForEachPair (matches, "_tag_", "_value_", "
        (`&result`) [_value_[utility.getGlobalValue('terms.description')]] = parameters.GetProfileCI(_value_[utility.getGlobalValue('terms.id')], (`&fit`)[utility.getGlobalValue('terms.likelihood_function')], 0.95);
    ");
    
    return result;
}