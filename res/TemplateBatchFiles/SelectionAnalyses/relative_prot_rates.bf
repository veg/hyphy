RequireVersion("2.31");

// namespace 'utility' for convenience functions
LoadFunctionLibrary("libv3/UtilityFunctions.bf");

// namespace 'alignments' for alignment-related functions
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/models/protein/empirical.bf");
LoadFunctionLibrary("libv3/models/protein/REV.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");


/*------------------------------------------------------------------------------*/

utility.ToggleEnvVariable ("NORMALIZE_SEQUENCE_NAMES", 1);

relative_prot_rates.analysis_description = {
    "info": "For a fixed alignment and tree, infer **relative** site specific substitution rates,
    by first optimizing alignment-wide branch lengths, and then inferring a site-specific uniform tree scaler",
    "version": "0.1alpha",
    "reference": "@TBD. Analysis based on Rate4Site method: Pupko, T., Bell, R. E., Mayrose, I., Glaser, F. & Ben-Tal, N. relative_prot_rates: an algorithmic tool for the identification of functional regions in proteins by surface mapping of evolutionary determinants within their homologues. Bioinformatics 18, S71–S77 (2002).",
    "authors": "Sergei L Kosakovsky Pond and Stephanie J Spielman",
    "contact": "{spond,stephanie.spielman}@temple.edu"
};

io.DisplayAnalysisBanner(relative_prot_rates.analysis_description);

SetDialogPrompt ("Specify a protein multiple sequence alignment file");


relative_prot_rates.alignment_info       = alignments.ReadNucleotideDataSet ("relative_prot_rates.dataset", None);
relative_prot_rates.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (relative_prot_rates.alignment_info[utility.getGlobalValue("terms.json.partitions")], None);
relative_prot_rates.partition_count = Abs (relative_prot_rates.partitions_and_trees);

io.CheckAssertion ("relative_prot_rates.partition_count==1", "This analysis can only handle a single partition");


io.ReportProgressMessageMD ("relative_prot_rates", "Data", "Input alignment description");
io.ReportProgressMessageMD ("relative_prot_rates", "Data", "Loaded **" +
                            relative_prot_rates.alignment_info ["sequences"] + "** sequences, **" +
                            relative_prot_rates.alignment_info ["sites"] + "** sites, and **" + relative_prot_rates.partition_count + "** partitions from \`" + relative_prot_rates.alignment_info ["file"] + "\`");

relative_prot_rates.filter_specification = alignments.DefineFiltersForPartitions (relative_prot_rates.partitions_and_trees, "relative_prot_rates.dataset" , "relative_prot_rates.filter.", relative_prot_rates.alignment_info);

relative_prot_rates.model_generator    = "models.protein.WAG.ModelDescription";
// TODO : this needs to be user selectable at run time from the list of available models 

io.ReportProgressMessageMD ("relative_prot_rates", "overall", "Obtaining alignment-wide branch-length estimates");

relative_prot_rates.trees = utility.Map (relative_prot_rates.partitions_and_trees, "_value_", "_value_['tree']"); // value => value['tree']
relative_prot_rates.filter_names = utility.Map (relative_prot_rates.filter_specification, "_value_", "_value_['name']"); // value => value['name']
relative_prot_rates.alignment_wide_MLES = estimators.FitSingleModel_Ext (
                                                          relative_prot_rates.filter_names,
                                                          relative_prot_rates.trees, 
                                                          relative_prot_rates.model_generator,
                                                          None,
                                                          None);
                                                          
estimators.fixSubsetOfEstimates(relative_prot_rates.alignment_wide_MLES, relative_prot_rates.alignment_wide_MLES[terms.global]);


io.ReportProgressMessageMD ("relative_prot_rates", "overall", ">Fitted an alignment-wide model. **Log-L = " + relative_prot_rates.alignment_wide_MLES ['LogL'] + "**.");


/** 
	Set up the table to display to the screen 
*/


relative_prot_rates.table_screen_output  = {{"Site", "Rel. rate (MLE)", "95% profile likelihood CI"}};
relative_prot_rates.table_output_options = {"header" : TRUE, "min-column-width" : 16, "align" : "center"};

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
        
model.ApplyModelToTree( "relative_prot_rates.site_tree", relative_prot_rates.trees[0], {"default" : relative_prot_rates.site_model}, None);

// create a site filter; this is an ugly hack for the time being
// alignments.serialize_site_filter returns HBL code as string in 
// which the function `__make_filter` is defined.
 
ExecuteCommands (alignments.serialize_site_filter (
								   relative_prot_rates.filter_names[0],
								   ((relative_prot_rates.site_patterns[0])["sites"])[0]));

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
					

relative_prot_rates.queue = mpi.CreateQueue ({"LikelihoodFunctions": {{"relative_prot_rates.site_likelihood"}},
								    "Models" : {{"relative_prot_rates.site_model"}},
								    "Headers" : utility.GetListOfLoadedModules (),
								    "Variables" : {{"relative_prot_rates.site_model_scaler_name"}}
							 });

/* run the main loop over all unique site pattern combinations */
utility.ForEachPair (relative_prot_rates.site_patterns, "_pattern_", "_pattern_info_",
	'
		mpi.QueueJob (relative_prot_rates.queue, "relative_prot_rates.handle_a_site", {"0" : "relative_prot_rates.site_likelihood",
														"1" : alignments.serialize_site_filter
														   ((relative_prot_rates.filter_specification[0])["name"],
														   (_pattern_info_["sites"])[0]),
														"2" : _pattern_info_,
														"3" : relative_prot_rates.site_model_mapping
														},
														"relative_prot_rates.store_results");
	'
);

mpi.QueueComplete (relative_prot_rates.queue);

relative_prot_rates.stats = (math.GatherDescriptiveStats(utility.Map (utility.Values(utility.Map (relative_prot_rates.rate_estimates, "_value_", "_value_[terms.MLE]")), "_value_", "0+_value_")));


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
io.ReportProgressMessageMD ("relative_prot_rates", "Stats", "* **Mean**: "  + Format (relative_prot_rates.stats["Mean"], 6, 2));
io.ReportProgressMessageMD ("relative_prot_rates", "Stats", "* **Median**: "  + Format (relative_prot_rates.stats["Median"], 6, 2));
io.ReportProgressMessageMD ("relative_prot_rates", "Stats", "* **Std.Dev**: "  + Format (relative_prot_rates.stats["Std.Dev"], 6, 2));
io.ReportProgressMessageMD ("relative_prot_rates", "Stats", "* **95% Range**: ["  + Format (relative_prot_rates.stats["2.5%"], 5,2) + "," + Format (relative_prot_rates.stats["97.5%"], 5,2) + "]");
/* Without formatting:
io.ReportProgressMessageMD ("relative_prot_rates", "Stats", "* [Mean] "  + relative_prot_rates.stats["Mean"]);
io.ReportProgressMessageMD ("relative_prot_rates", "Stats", "* [Median] "  + relative_prot_rates.stats["Median"]);
io.ReportProgressMessageMD ("relative_prot_rates", "Stats", "* [Std.Dev] "  + relative_prot_rates.stats["Std.Dev"]);
io.ReportProgressMessageMD ("relative_prot_rates", "Stats", "* [95% Range] "  + relative_prot_rates.stats["2.5%"] + ":" + relative_prot_rates.stats["97.5%"]);
*/

							
io.SpoolJSON ({
				'Relative site rate estimates' : relative_prot_rates.rate_estimates, 
				'alignment' : relative_prot_rates.alignment_info["file"],
				'analysis' : relative_prot_rates.analysis_description
				}
				, relative_prot_rates.alignment_info["file"] + ".site-rates.json");

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
    
    if (pattern_info ["is_constant"]) {
    	// the MLE for a constant site is 0; 
    	// only the CI is non-trivial
		^(^"relative_prot_rates.site_model_scaler_name") = 0;
    
    } else {
    
		^(^"relative_prot_rates.site_model_scaler_name") = 1;
		Optimize (results, ^lf);
	}
	
    return parameters.GetProfileCI (^"relative_prot_rates.site_model_scaler_name", lf, 0.95);
}

// handle result processing

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------


lfunction relative_prot_rates.store_results (node, result, arguments) {
    pattern_info    = arguments [2];


	if ((^'relative_prot_rates.table_output_options')["header"]) {
		
		io.ReportProgressMessageMD ("relative_prot_rates", "sites", "Site rate estimates and associated 95% profile likelihood estimates\n"); 

		fprintf (stdout,
			io.FormatTableRow (^'relative_prot_rates.table_screen_output',^'relative_prot_rates.table_output_options'));
		(^'relative_prot_rates.table_output_options')["header"] = FALSE;
	}
	
	
	utility.ForEach (pattern_info["sites"], "_site_index_",
		"
			relative_prot_rates.rate_estimates [_site_index_+1] = `&result`;
			result_row = {1,3};
			result_row [0] = '' + (_site_index_ + 1);
			result_row [1] = Format((`&result`)[terms.MLE],6,3);
			result_row [2] = Format((`&result`)[terms.lower_bound],6,3) + ' :'  +Format((`&result`)[terms.upper_bound],6,3);
			fprintf (stdout,
				io.FormatTableRow (result_row,relative_prot_rates.table_output_options));
		"
    );

	return rate_statistics;

}



