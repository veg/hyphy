RequireVersion("2.3.5");

LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");
LoadFunctionLibrary("libv3/all-terms.bf");

LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");
LoadFunctionLibrary("libv3/models/rate_variation.bf");

LoadFunctionLibrary("libv3/models/DNA.bf");
LoadFunctionLibrary("libv3/models/DNA/GTR.bf");
LoadFunctionLibrary("libv3/models/DNA/HKY85.bf");
LoadFunctionLibrary("libv3/models/DNA/JC69.bf");
LoadFunctionLibrary("libv3/models/protein.bf");
LoadFunctionLibrary("libv3/models/protein/empirical.bf");

// for JSON storage compatibility
LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");

/*------------------------------------------------------------------------------*/

utility.ToggleEnvVariable ("NORMALIZE_SEQUENCE_NAMES", 1);

leisr.analysis_description = {
    terms.io.info: "LEISR (Likelihood Estimation of Individual Site Rates) infer relative amino-acid or nucleotide rates from a fixed nucleotide or amino-acid alignment and tree. Relative site-specific substitution rates are
    inferred by first optimizing alignment-wide branch lengths, and then inferring a site-specific uniform tree scaler",
    terms.io.version: "0.1alpha",
    terms.io.reference: "Spielman, S.J. and Kosakovsky Pond, S.L. Relative evolutionary rate inference in HyPhy with LEISR. bioRxiv. https://doi.org/10.1101/206011. (2017); Pupko, T., Bell, R. E., Mayrose, I., Glaser, F. & Ben-Tal, N. Rate4Site: an algorithmic tool for the identification of functional regions in proteins by surface mapping of evolutionary determinants within their homologues. Bioinformatics 18, S71–S77 (2002).",
    terms.io.authors: "Sergei L Kosakovsky Pond and Stephanie J Spielman",
    terms.io.contact: "{spond,stephanie.spielman}@temple.edu"
};

io.DisplayAnalysisBanner(leisr.analysis_description);
/*******************************************************************************************************************/


/***************************************** LOAD DATASET **********************************************************/
SetDialogPrompt ("Specify a multiple sequence alignment file");
leisr.alignment_info  = alignments.ReadNucleotideDataSet ("leisr.dataset", None);

leisr.name_mapping = leisr.alignment_info[utility.getGlobalValue("terms.data.name_mapping")];
if (None == leisr.name_mapping) {
    leisr.name_mapping = {};
    utility.ForEach (alignments.GetSequenceNames ("leisr.dataset"), "_value_", "`&leisr.name_mapping`[_value_] = _value_");
}
leisr.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (leisr.alignment_info[utility.getGlobalValue("terms.data.partitions")], leisr.name_mapping);
leisr.partition_count = Abs (leisr.partitions_and_trees);

io.CheckAssertion ("leisr.partition_count==1", "This analysis can only handle a single partition");


io.ReportProgressMessageMD ("relative_rates", "Data", "Input alignment description");
io.ReportProgressMessageMD ("relative_rates", "Data", "Loaded **" +
                            leisr.alignment_info [terms.data.sequences] + "** sequences, **" +
                            leisr.alignment_info [terms.data.sites] + "** sites, and **" + leisr.partition_count + "** partitions from \`" + leisr.alignment_info [terms.data.file] + "\`");
leisr.filter_specification = alignments.DefineFiltersForPartitions (leisr.partitions_and_trees, "leisr.dataset" , "leisr.filter.", leisr.alignment_info);
/*******************************************************************************************************************/


/***************************************** MODEL SELECTION **********************************************************/

leisr.protein_type    = "Protein";
leisr.nucleotide_type = "Nucleotide";
leisr.analysis_type  = io.SelectAnOption ({{leisr.protein_type , "Infer relative rates from a protein (amino-acid) alignment"}, {leisr.nucleotide_type, "Infer relative rates from a nucleotide alignment"}},
                                                    "Select your analysis type:");


if (leisr.analysis_type ==  leisr.protein_type) {
    leisr.baseline_model  = io.SelectAnOption (models.protein.empirical_models, "Select a protein model:");
    leisr.generators = models.protein.empirical.default_generators;
}
else {

    leisr.baseline_model  = io.SelectAnOption (models.DNA.models, "Select a nucleotide model:");
    leisr.generators = models.DNA.generators;
}

leisr.use_rate_variation = io.SelectAnOption( {{"Gamma", "Use a four-category discrete gamma distribution when optimizing branch lengths."},
                                                    {"GDD", "Use a four-category general discrete distribution when optimizing branch lengths."},
                                                    {"No", "Do not consider rate variation when optimizing branch lengths."}
                                                    },
                                                    "Optimize branch lengths with rate variation?");


function leisr.Baseline.ModelDescription(type){
    def = Call( leisr.generators[leisr.baseline_model], type);
    return def;
}

function leisr.Baseline.ModelDescription.withGamma(type){
    def = leisr.Baseline.ModelDescription(type);
	def [terms.model.rate_variation] = rate_variation.types.Gamma.factory ({terms.rate_variation.bins : 4});
    return def;
}
function leisr.Baseline.ModelDescription.withGDD4(type){
    def = leisr.Baseline.ModelDescription(type);
	def [terms.model.rate_variation] = rate_variation.types.GDD.factory ({terms.rate_variation.bins : 4});
    return def;
}


leisr.baseline_model_name = leisr.baseline_model;
if (leisr.analysis_type ==  leisr.protein_type) {
    leisr.baseline_model_name  = leisr.baseline_model_name + "F";
}
if (leisr.use_rate_variation == "Gamma"){
    leisr.baseline_model_name      = leisr.baseline_model_name + "+4Gamma";
    leisr.baseline_model_desc      = "leisr.Baseline.ModelDescription.withGamma";
}
else {
    if (leisr.use_rate_variation == "GDD"){
        leisr.baseline_model_name      = leisr.baseline_model_name + "+4GDD";
        leisr.baseline_model_desc      = "leisr.Baseline.ModelDescription.withGDD4";
    }
    else {
        leisr.baseline_model_name      = leisr.baseline_model_name;
        leisr.baseline_model_desc      = "leisr.Baseline.ModelDescription";
    }
}
/*******************************************************************************************************************/



/***************************************** INFERENCE **********************************************************/


io.ReportProgressMessageMD ("relative_rates", "overall", "Obtaining alignment-wide branch-length estimates");

leisr.trees = utility.Map (leisr.partitions_and_trees, "_value_", "_value_[terms.data.tree]"); // value => value['tree']
leisr.filter_names = utility.Map (leisr.filter_specification, "_value_", "_value_[terms.data.name]"); // value => value['name']
leisr.alignment_wide_MLES = estimators.FitSingleModel_Ext (
                                                          leisr.filter_names,
                                                          leisr.trees,
                                                          leisr.baseline_model_desc,
                                                          None,
                                                          None);



estimators.fixSubsetOfEstimates(leisr.alignment_wide_MLES, leisr.alignment_wide_MLES[terms.global]);

io.ReportProgressMessageMD ("relative_rates", "overall", ">Fitted an alignment-wide model. **Log-L = " + leisr.alignment_wide_MLES [terms.fit.log_likelihood] + "**.");

/**
	Set up the table to display to the screen
*/


leisr.table_screen_output  = {{"Site", "Rel. rate (MLE)", "95% profile likelihood CI"}};
leisr.table_output_options = {terms.table_options.header : TRUE, terms.table_options.minimum_column_width : 16, terms.table_options.align : "center"};
leisr.table_headers = {{"MLE", "Relative rate estimate at a site"}
                       {"Lower", "Lower bound of 95% profile likelihood CI"}
                       {"Upper", "Upper bound of 95% profile likelihood CI"}};
leisr.site_results = {leisr.alignment_info [terms.data.sites], Rows (leisr.table_headers)};      

leisr.site_patterns = alignments.Extract_site_patterns (leisr.filter_names[0]);

// set-up model for site-level fitting in the next couple of lines, where rv turned off
leisr.site_model = model.generic.DefineModel("leisr.Baseline.ModelDescription",
        "relative_rates_site_model_instance", {
            "0": parameters.Quote(terms.global),
        },
        leisr.filter_names[0],
        None);



leisr.site_model_mapping = {"relative_rates_site_model_instance" : leisr.site_model};

// leisr.site_tree is created from the information in  leisr.trees[0]
// and populated with (the default) model
model.ApplyModelToTree( "leisr.site_tree", leisr.trees[0], {terms.default : leisr.site_model}, None);

// create a site filter; this is an ugly hack for the time being
// alignments.serialize_site_filter returns HBL code as string in
// which the function `__make_filter` is defined.
ExecuteCommands (alignments.serialize_site_filter (
								   leisr.filter_names[0],
								   ((leisr.site_patterns[0])[terms.data.sites])[0]));

__make_filter ("leisr.site_filter");

LikelihoodFunction leisr.site_likelihood = (leisr.site_filter, leisr.site_tree);

leisr.site_model_scaler_name = "leisr.site_rate_estimate";

leisr.rate_estimates = {};

/**
	 this will store site estimates, which will then be dumped to JSON
*/

parameters.DeclareGlobal (leisr.site_model_scaler_name, None);

estimators.ApplyExistingEstimates ("leisr.site_likelihood", leisr.site_model_mapping, leisr.alignment_wide_MLES,
									 {"0" : leisr.site_model_scaler_name} // proportional scaler
									);


leisr.queue = mpi.CreateQueue ({terms.mpi.LikelihoodFunctions: {{"leisr.site_likelihood"}},
								    terms.mpi.Models : {{"leisr.site_model"}},
								    terms.mpi.Headers : utility.GetListOfLoadedModules ("libv3/"),
								    terms.mpi.Variables : {{"leisr.site_model_scaler_name"}}
							 });

/* run the main loop over all unique site pattern combinations */
utility.ForEachPair (leisr.site_patterns, "_pattern_", "_pattern_info_",
	'
		mpi.QueueJob (leisr.queue, "leisr.handle_a_site", {"0" : "leisr.site_likelihood",
														"1" : alignments.serialize_site_filter
														   ((leisr.filter_specification[0])[terms.data.name],
														   (_pattern_info_[terms.data.sites])[0]),
														"2" : _pattern_info_,
														"3" : leisr.site_model_mapping
														},
														"leisr.store_results");
	'
);

mpi.QueueComplete (leisr.queue);

leisr.site_rates = utility.Map( utility.Values(utility.Map (leisr.rate_estimates, "_value_", "_value_[terms.fit.MLE]")), "_value_", "0+_value_");
leisr.stats = math.GatherDescriptiveStats(leisr.site_rates);

io.ReportProgressMessageMD ("relative_rates", "Stats", "Rate distribution summary");
io.ReportProgressMessageMD ("relative_rates", "Stats", "* **Mean**: "  + Format (leisr.stats[terms.math.mean], 6, 2));
io.ReportProgressMessageMD ("relative_rates", "Stats", "* **Median**: "  + Format (leisr.stats[terms.math.median], 6, 2));
io.ReportProgressMessageMD ("relative_rates", "Stats", "* **Std.Dev**: "  + Format (leisr.stats[terms.math.stddev], 6, 2));
io.ReportProgressMessageMD ("relative_rates", "Stats", "* **95% Range**: ["  + Format (leisr.stats[terms.math._2.5], 5,2) + "," + Format (leisr.stats[terms.math._97.5], 5,2) + "]");



/************************************* JSON STORAGE ***********************************/

leisr.store_global = utility.Map(
                        utility.Map (leisr.alignment_wide_MLES[utility.getGlobalValue("terms.global")], "_value_", '   {terms.fit.MLE : _value_[terms.fit.MLE]}'),
                        "_value_",
                        "_value_[terms.fit.MLE]");

if (Abs(((leisr.partitions_and_trees["0"])[terms.data.tree])[terms.branch_length]) == 0){
    store_tree = utility.Map (leisr.partitions_and_trees, "_pt_", '(_pt_[terms.data.tree])[terms.trees.newick]');
} else {    
    store_tree = utility.Map (leisr.partitions_and_trees, "_pt_", '(_pt_[terms.data.tree])[terms.trees.newick_with_lengths]');
}
                       
leisr.aicc = leisr.getIC(leisr.alignment_wide_MLES[terms.fit.log_likelihood], leisr.alignment_wide_MLES[terms.parameters], leisr.alignment_info[utility.getGlobalValue("terms.data.sites")] * leisr.alignment_info[utility.getGlobalValue("terms.data.sequences")]);
leisr.json_content = { terms.json.input :
                            {terms.json.file: leisr.alignment_info[terms.data.file],
                                  terms.json.sequences: leisr.alignment_info[terms.data.sequences],
                                  terms.json.sites:     leisr.alignment_info[terms.data.sites],
                                  terms.json.partition_count: leisr.partition_count,
                                  terms.json.trees: store_tree
                            },
                        terms.json.analysis : leisr.analysis_description,
				        terms.json.fits:
				        {
				            leisr.baseline_model_name:
                                {
                                    terms.json.log_likelihood: leisr.alignment_wide_MLES[terms.fit.log_likelihood],
                                    terms.json.parameters: leisr.alignment_wide_MLES[terms.parameters],
                                    terms.json.AICc: leisr.aicc,
                                    terms.json.rate_distribution: leisr.store_global,
                                    terms.efv_estimate: (leisr.alignment_wide_MLES[utility.getGlobalValue("terms.efv_estimate")])["VALUEINDEXORDER"][0]
				                }
				        },
				        terms.json.MLE : {terms.json.headers   : leisr.table_headers,
                                            terms.json.content : {"0":leisr.site_results}
                                        }                   
 	                   };                                        
selection.io.json_store_branch_attribute(leisr.json_content, utility.getGlobalValue ("terms.original_name"), utility.getGlobalValue ("terms.json.node_label"), -1, 0, utility.getGlobalValue ("leisr.name_mapping"));
selection.io.json_store_branch_attribute(leisr.json_content, utility.getGlobalValue ("leisr.baseline_model_name"), utility.getGlobalValue ("terms.branch_length"), 0, 0, 
                                         selection.io.extract_branch_info((leisr.alignment_wide_MLES[terms.branch_length])[0], "selection.io.branch.length")
                                         );

io.SpoolJSON (leisr.json_content, leisr.alignment_info[terms.data.file] + ".LEISR.json");


//----------------------------------------------------------------------------------------
// HANDLERS
//----------------------------------------------------------------------------------------

// fit a rate at a single site

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

lfunction leisr.handle_a_site (lf, filter_data, pattern_info, model_mapping) {

    GetString (lfInfo, ^lf,-1);
    ExecuteCommands (filter_data);

    __make_filter ((lfInfo["Datafilters"])[0]);
    utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);

    if (pattern_info [utility.getGlobalValue("terms.data.is_constant")]) {
    	// the MLE for a constant site is 0;
    	// only the CI is non-trivial
		^(utility.getGlobalValue("leisr.site_model_scaler_name")) = 0;

    } else {

		^(utility.getGlobalValue("leisr.site_model_scaler_name")) = 1;
		Optimize (results, ^lf);
	}
    return parameters.GetProfileCI (utility.getGlobalValue("leisr.site_model_scaler_name"), lf, 0.95);
}


// handle result processing

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------


lfunction leisr.store_results (node, result, arguments) {
    pattern_info    = arguments [2];


	if ((^'leisr.table_output_options')[utility.getGlobalValue("terms.table_options.header")]) {

		io.ReportProgressMessageMD ("relative_rates", "sites", "Site rate estimates and associated 95% profile likelihood estimates\n");

		fprintf (stdout,
			io.FormatTableRow (^'leisr.table_screen_output',^'leisr.table_output_options'));
		(^'leisr.table_output_options')[utility.getGlobalValue("terms.table_options.header")] = FALSE;
	}

	utility.ForEach (pattern_info[utility.getGlobalValue("terms.data.sites")], "_site_index_",
		"
			leisr.rate_estimates [_site_index_+1] = `&result`;
			result_row = {1,3};
			result_row [0] = '' + (_site_index_ + 1);
			result_row [1] = Format((`&result`)[terms.fit.MLE],6,3);
			result_row [2] = Format((`&result`)[terms.lower_bound],6,3) + ' :'  +Format((`&result`)[terms.upper_bound],6,3);
			fprintf (stdout,
				io.FormatTableRow (result_row,leisr.table_output_options));
		
		    // JSON-related
		    //console.log(_site_index_);
		    leisr.site_results[_site_index_][0] = (leisr.rate_estimates[_site_index_+1])[terms.fit.MLE];
		    leisr.site_results[_site_index_][1] = (leisr.rate_estimates[_site_index_+1])[terms.lower_bound];
		    leisr.site_results[_site_index_][2] = (leisr.rate_estimates[_site_index_+1])[terms.upper_bound];
		"
    );

	return rate_statistics;

}



lfunction leisr.getIC(logl, params, samples) {
    return -2 * logl + 2 * samples / (samples - params - 1) * params;
}