RequireVersion("2.3.4");

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

LoadFunctionLibrary("libv3/models/protein/empirical.bf");
LoadFunctionLibrary("libv3/models/protein.bf");

/*------------------------------------------------------------------------------*/

utility.ToggleEnvVariable ("NORMALIZE_SEQUENCE_NAMES", 1);

relative_prot_rates.analysis_description = {
    terms.io.info: "RELprot (RELative protein rates) infers, for a fixed alignment and tree, **relative** site specific substitution rates,
    by first optimizing alignment-wide branch lengths, and then inferring a site-specific uniform tree scaler",
    terms.io.version: "0.1alpha",
    terms.io.reference: "@TBD. Analysis based on Rate4Site method: Pupko, T., Bell, R. E., Mayrose, I., Glaser, F. & Ben-Tal, N. Rate4Site: an algorithmic tool for the identification of functional regions in proteins by surface mapping of evolutionary determinants within their homologues. Bioinformatics 18, S71–S77 (2002).",
    terms.io.authors: "Sergei L Kosakovsky Pond and Stephanie J Spielman",
    terms.io.contact: "{spond,stephanie.spielman}@temple.edu"
};

io.DisplayAnalysisBanner(relative_prot_rates.analysis_description);


console.log("**WARNING**: This analysis will be ported to a new batchfile in a future HyPhy release. Please consider executing the batchfile \`relative_rates_scaler.bf\` instead.");
console.log("");

/***************************************** LOAD DATASET **********************************************************/
SetDialogPrompt ("Specify a protein multiple sequence alignment file");


relative_prot_rates.alignment_info       = alignments.ReadNucleotideDataSet ("relative_prot_rates.dataset", None);

name_mapping = relative_prot_rates.alignment_info[utility.getGlobalValue("terms.data.name_mapping")];
if (None == name_mapping) {  
    name_mapping = {};
    utility.ForEach (alignments.GetSequenceNames ("relative_prot_rates.dataset"), "_value_", "`&name_mapping`[_value_] = _value_");
} 
relative_prot_rates.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (relative_prot_rates.alignment_info[utility.getGlobalValue("terms.data.partitions")], name_mapping);

relative_prot_rates.partition_count      = Abs (relative_prot_rates.partitions_and_trees);
io.CheckAssertion ("relative_prot_rates.partition_count==1", "This analysis can only handle a single partition");


io.ReportProgressMessageMD ("relative_prot_rates", "Data", "Input alignment description");
io.ReportProgressMessageMD ("relative_prot_rates", "Data", "Loaded **" +
                            relative_prot_rates.alignment_info [terms.data.sequences] + "** sequences, **" +
                            relative_prot_rates.alignment_info [terms.data.sites] + "** sites, and **" + relative_prot_rates.partition_count + "** partitions from \`" + relative_prot_rates.alignment_info [terms.data.file] + "\`");


relative_prot_rates.filter_specification = alignments.DefineFiltersForPartitions (relative_prot_rates.partitions_and_trees, "relative_prot_rates.dataset" , "relative_prot_rates.filter.", relative_prot_rates.alignment_info);
/**********************************************************************************************************************************************/




/***************************************** SELECT MODEL, +F **********************************************************/
/** NOTE: there is no rate variation option here ever. This must be discouraged. **/

relative_prot_rates.baseline_model = io.SelectAnOption (models.protein.empirical_models,
                                                     "Select an empirical protein model for fitting rates (we recommend JC69 to avoid matrix-induced biases in inferred rates):");

// "Yes", "No"
relative_prot_rates.plusF          = io.SelectAnOption ({{"Yes", "Use empirical (+F) amino-acid frequencies ."}, {"No", "Use default amino-acid frequencies."}},                 
                                                      "Use a +F model for initial bl optimization? (recommended no for a simple JC69 model)");


// Set up model generator and name as +F or not.
if (relative_prot_rates.plusF == "Yes"){
    relative_prot_rates.full_model_name  =  relative_prot_rates.baseline_model + "+F";
    relative_prot_rates.model_generator = models.protein.empirical.plusF_generators[relative_prot_rates.baseline_model];
}
else {
    relative_prot_rates.full_model_name  =  relative_prot_rates.baseline_model;
    relative_prot_rates.model_generator = models.protein.empirical.default_generators[relative_prot_rates.baseline_model];
}


// Read in data
relative_prot_rates.trees = utility.Map (relative_prot_rates.partitions_and_trees, "_value_", "_value_[terms.data.tree]"); // value => value['tree']
relative_prot_rates.filter_names = utility.Map (relative_prot_rates.filter_specification, "_value_", "_value_[terms.data.name]"); // value => value['name']

/*********************************************************************************************************************************************/




io.ReportProgressMessageMD ("relative_prot_rates", "overall", "Obtaining alignment-wide branch-length estimates");




relative_prot_rates.alignment_wide_MLES = estimators.FitSingleModel_Ext (
                                                          relative_prot_rates.filter_names,
                                                          relative_prot_rates.trees, 
                                                          relative_prot_rates.model_generator,
                                                          None,
                                                          None);
                                                     
estimators.fixSubsetOfEstimates(relative_prot_rates.alignment_wide_MLES, relative_prot_rates.alignment_wide_MLES[terms.global]);

io.ReportProgressMessageMD ("relative_prot_rates", "overall", ">Fitted an alignment-wide model. **Log-L = " + relative_prot_rates.alignment_wide_MLES [terms.fit.log_likelihood] + "**.");


/** 
	Set up the table to display to the screen 
*/


relative_prot_rates.table_screen_output  = {{"Site", "Rel. rate (MLE)", "95% profile likelihood CI"}};
relative_prot_rates.table_output_options = {terms.table_options.header : TRUE, terms.table_options.minimum_column_width : 16, terms.table_options.align : "center"};

relative_prot_rates.site_patterns = alignments.Extract_site_patterns (relative_prot_rates.filter_names[0]);

/**
	relative_prot_rates.site_patterns  maps a unique site pattern (by index) to 
		-- the list of sites that match the pattern
		-- whether or not it is constant (invariable)
	

	{
	 "0":{
	   "sites":{
		 "0":0
		},
	   "is_constant":0
	  },
*/

// set-up model for site-level fitting in the next couple of lines

relative_prot_rates.site_model = model.generic.DefineModel(relative_prot_rates.model_generator,
        "relative_prot_rates_site_model_instance", {
            "0": parameters.Quote(terms.global),
        },
        relative_prot_rates.filter_names[0],
        None);
        

relative_prot_rates.site_model_mapping = {"relative_prot_rates_site_model_instance" : relative_prot_rates.site_model};
        
// relative_prot_rates.site_tree is created from the information in  relative_prot_rates.trees[0]
// and populated with (the default) model        
        
model.ApplyModelToTree( "relative_prot_rates.site_tree", relative_prot_rates.trees[0], {terms.default : relative_prot_rates.site_model}, None);

// create a site filter; this is an ugly hack for the time being
// alignments.serialize_site_filter returns HBL code as string in 
// which the function `__make_filter` is defined.
 
ExecuteCommands (alignments.serialize_site_filter (
								   relative_prot_rates.filter_names[0],
								   ((relative_prot_rates.site_patterns[0])[terms.data.sites])[0]));

__make_filter ("relative_prot_rates.site_filter");

LikelihoodFunction relative_prot_rates.site_likelihood = (relative_prot_rates.site_filter, relative_prot_rates.site_tree);

relative_prot_rates.site_model_scaler_name = "relative_prot_rates.site_rate_estimate";

relative_prot_rates.rate_estimates = {}; 

/**
	 this will store site estimates, which will then be dumped to JSON
*/

parameters.DeclareGlobal (relative_prot_rates.site_model_scaler_name, None);

estimators.ApplyExistingEstimates ("relative_prot_rates.site_likelihood", relative_prot_rates.site_model_mapping, relative_prot_rates.alignment_wide_MLES,
									 {"0" : relative_prot_rates.site_model_scaler_name} // proportional scaler
									);
					

relative_prot_rates.queue = mpi.CreateQueue ({terms.mpi.LikelihoodFunctions: {{"relative_prot_rates.site_likelihood"}},
								    terms.mpi.Models : {{"relative_prot_rates.site_model"}},
								    terms.mpi.Headers : utility.GetListOfLoadedModules ("libv3/"),
								    terms.mpi.Variables : {{"relative_prot_rates.site_model_scaler_name"}}
							 });

/* run the main loop over all unique site pattern combinations */
utility.ForEachPair (relative_prot_rates.site_patterns, "_pattern_", "_pattern_info_",
	'
		mpi.QueueJob (relative_prot_rates.queue, "relative_prot_rates.handle_a_site", {"0" : "relative_prot_rates.site_likelihood",
														"1" : alignments.serialize_site_filter
														   ((relative_prot_rates.filter_specification[0])[terms.data.name],
														   (_pattern_info_[terms.data.sites])[0]),
														"2" : _pattern_info_,
														"3" : relative_prot_rates.site_model_mapping
														},
														"relative_prot_rates.store_results");
	'
);

mpi.QueueComplete (relative_prot_rates.queue);



relative_prot_rates.site_rates = utility.Map( utility.UniqueValues(utility.Map (relative_prot_rates.rate_estimates, "_value_", "_value_[terms.fit.MLE]")), "_value_", "0+_value_");
relative_prot_rates.stats = math.GatherDescriptiveStats(relative_prot_rates.site_rates);

/*
{
 "Count":148,
 "Mean":1.016735774503485,
 "Median":0.9798409545051374,
 "Min":0.4983360276068412,
 "Max":1.729838344017373,
 "2.5%":0.6273256414437858,
 "97.5%":1.434471693559964,
 "Sum":150.4768946265158,
 "Std.Dev":0.2171843066029547,
 "Variance":0.04716902303460622,
 "COV":0.2136093880526776,
 "Skewness":0.492213048250155,
 "Kurtosis":3.242020650192995,
 "Sq. sum":300.9537892530316,
 "Non-negative":148
}
*/

io.ReportProgressMessageMD ("relative_prot_rates", "Stats", "Rate distribution summary");
io.ReportProgressMessageMD ("relative_prot_rates", "Stats", "* **Mean**: "  + Format (relative_prot_rates.stats[terms.math.mean], 6, 2));
io.ReportProgressMessageMD ("relative_prot_rates", "Stats", "* **Median**: "  + Format (relative_prot_rates.stats[terms.math.median], 6, 2));
io.ReportProgressMessageMD ("relative_prot_rates", "Stats", "* **Std.Dev**: "  + Format (relative_prot_rates.stats[terms.math.stddev], 6, 2));
io.ReportProgressMessageMD ("relative_prot_rates", "Stats", "* **95% Range**: ["  + Format (relative_prot_rates.stats[terms.math._2.5], 5,2) + "," + Format (relative_prot_rates.stats[terms.math._97.5], 5,2) + "]");
/* Without formatting:
io.ReportProgressMessageMD ("relative_prot_rates", "Stats", "* [Mean] "  + relative_prot_rates.stats["Mean"]);
io.ReportProgressMessageMD ("relative_prot_rates", "Stats", "* [Median] "  + relative_prot_rates.stats["Median"]);
io.ReportProgressMessageMD ("relative_prot_rates", "Stats", "* [Std.Dev] "  + relative_prot_rates.stats["Std.Dev"]);
io.ReportProgressMessageMD ("relative_prot_rates", "Stats", "* [95% Range] "  + relative_prot_rates.stats["2.5%"] + ":" + relative_prot_rates.stats["97.5%"]);
*/



tree_definition   = utility.Map (relative_prot_rates.partitions_and_trees, "_partition_", '_partition_[terms.data.tree]');

io.SpoolJSON ({ terms.json.input : {terms.json.file: relative_prot_rates.alignment_info[terms.data.file],
                          terms.json.sequences: relative_prot_rates.alignment_info[terms.data.sequences],
                          terms.json.sites:     relative_prot_rates.alignment_info[terms.data.sites],
                          terms.json.tree_string: (tree_definition[0])[terms.trees.newick_with_lengths]},
                terms.json.analysis : relative_prot_rates.analysis_description,       
				terms.json.relative_site_rates : relative_prot_rates.rate_estimates, 
				terms.json.global: {terms.json.model: relative_prot_rates.full_model_name,
				               terms.efv_estimate: (relative_prot_rates.alignment_wide_MLES[utility.getGlobalValue("terms.efv_estimate")])["VALUEINDEXORDER"][0],
				               terms.json.tree_string: (relative_prot_rates.alignment_wide_MLES[terms.fit.trees])[0],
				               terms.json.log_likelihood: relative_prot_rates.alignment_wide_MLES[terms.fit.log_likelihood]}
				},
				relative_prot_rates.alignment_info[terms.data.file] + ".site-rates.json");


//----------------------------------------------------------------------------------------
// HANDLERS
//----------------------------------------------------------------------------------------

// fit a rate at a single site

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

lfunction relative_prot_rates.handle_a_site (lf, filter_data, pattern_info, model_mapping) {
	

    GetString (lfInfo, ^lf,-1);
    ExecuteCommands (filter_data);
    
    __make_filter ((lfInfo["Datafilters"])[0]);
    utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);
    
    if (pattern_info [utility.getGlobalValue("terms.data.is_constant")]) {
    	// the MLE for a constant site is 0; 
    	// only the CI is non-trivial
		^(utility.getGlobalValue("relative_prot_rates.site_model_scaler_name")) = 0;
    
    } else {
    
		^(utility.getGlobalValue("relative_prot_rates.site_model_scaler_name")) = 1;
		Optimize (results, ^lf);
	}
	
    return parameters.GetProfileCI (utility.getGlobalValue("relative_prot_rates.site_model_scaler_name"), lf, 0.95);
}

// handle result processing

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------


lfunction relative_prot_rates.store_results (node, result, arguments) {
    pattern_info    = arguments [2];


	if ((utility.getGlobalValue("relative_prot_rates.table_output_options"))[utility.getGlobalValue("terms.table_options.header")]) {
		
		io.ReportProgressMessageMD ("relative_prot_rates", "sites", "Site rate estimates and associated 95% profile likelihood estimates\n"); 

		fprintf (stdout,
			io.FormatTableRow (utility.getGlobalValue("relative_prot_rates.table_screen_output"),utility.getGlobalValue("relative_prot_rates.table_output_options")));
		(utility.getGlobalValue("relative_prot_rates.table_output_options"))[utility.getGlobalValue("terms.table_options.header")] = FALSE;
	}
	
	
	utility.ForEach (pattern_info[utility.getGlobalValue("terms.data.sites")], "_site_index_",
		"
			relative_prot_rates.rate_estimates [_site_index_+1] = `&result`;
			result_row = {1,3};
			result_row [0] = '' + (_site_index_ + 1);
			result_row [1] = Format((`&result`)[utility.getGlobalValue('terms.fit.MLE')],6,3);
			result_row [2] = Format((`&result`)[utility.getGlobalValue('terms.lower_bound')],6,3) + ' :'  +Format((`&result`)[utility.getGlobalValue('terms.upper_bound')],6,3);
			fprintf (stdout,
				io.FormatTableRow (result_row,relative_prot_rates.table_output_options));
		"
    );



	return rate_statistics;

}


