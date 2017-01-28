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

rate4site.analysis_description = {
    "info": "For a fixed alignment and tree, infer **relative** site specific substitution rates,
    by first optimizing alignment-wide branch lengths, and then inferring a site-specific uniform tree scaler",
    "version": "0.1alpha",
    "reference": "@TBD. Analysis based on Rate4Site method: Pupko, T., Bell, R. E., Mayrose, I., Glaser, F. & Ben-Tal, N. Rate4Site: an algorithmic tool for the identification of functional regions in proteins by surface mapping of evolutionary determinants within their homologues. Bioinformatics 18, S71–S77 (2002).",
    "authors": "Sergei L Kosakovsky Pond and Stephanie J Spielman",
    "contact": "{spond,stephanie.spielman}@temple.edu"
};

io.DisplayAnalysisBanner(rate4site.analysis_description);

SetDialogPrompt ("Specify a protein multiple sequence alignment file");


rate4site.alignment_info       = alignments.ReadNucleotideDataSet ("rate4site.dataset", None);
rate4site.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (rate4site.alignment_info[utility.getGlobalValue("terms.json.partitions")], None);
rate4site.partition_count = Abs (rate4site.partitions_and_trees);

io.CheckAssertion ("rate4site.partition_count==1", "This analysis can only handle a single partition");


io.ReportProgressMessageMD ("rate4site", "Data", "Input alignment description");
io.ReportProgressMessageMD ("rate4site", "Data", "Loaded **" +
                            rate4site.alignment_info ["sequences"] + "** sequences, **" +
                            rate4site.alignment_info ["sites"] + "** sites, and **" + rate4site.partition_count + "** partitions from \`" + rate4site.alignment_info ["file"] + "\`");

rate4site.filter_specification = alignments.DefineFiltersForPartitions (rate4site.partitions_and_trees, "rate4site.dataset" , "rate4site.filter.", rate4site.alignment_info);

rate4site.model_generator    = "models.protein.WAG.ModelDescription";
// TODO : this needs to be user selectable at run time from the list of available models 

io.ReportProgressMessageMD ("rate4site", "overall", "Obtaining alignment-wide branch-length estimates");

rate4site.trees = utility.Map (rate4site.partitions_and_trees, "_value_", "_value_['tree']"); // value => value['tree']
rate4site.filter_names = utility.Map (rate4site.filter_specification, "_value_", "_value_['name']"); // value => value['name']
rate4site.alignment_wide_MLES = estimators.FitSingleModel_Ext (
                                                          rate4site.filter_names,
                                                          rate4site.trees, 
                                                          rate4site.model_generator,
                                                          None,
                                                          None);
                                                          
estimators.fixSubsetOfEstimates(rate4site.alignment_wide_MLES, rate4site.alignment_wide_MLES[terms.global]);


io.ReportProgressMessageMD ("rate4site", "overall", ">Fitted an alignment-wide model. **Log-L = " + rate4site.alignment_wide_MLES ['LogL'] + "**.");


/** 
	Set up the table to display to the screen 
*/


rate4site.table_screen_output  = {{"Site", "Rel. rate (MLE)", "95% profile likelihood CI"}};
rate4site.table_output_options = {"header" : TRUE, "min-column-width" : 16, "align" : "center"};

rate4site.site_patterns = alignments.Extract_site_patterns (rate4site.filter_names[0]);

/**
	rate4site.site_patterns  maps a unique site pattern (by index) to 
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

rate4site.site_model = model.generic.DefineModel(rate4site.model_generator,
        "rate4site_site_model_instance", {
            "0": parameters.Quote(terms.global),
        },
        rate4site.filter_names[0],
        None);
        

rate4site.site_model_mapping = {"rate4site_site_model_instance" : rate4site.site_model};
        
// rate4site.site_tree is created from the information in  rate4site.trees[0]
// and populated with (the default) model        
        
model.ApplyModelToTree( "rate4site.site_tree", rate4site.trees[0], {"default" : rate4site.site_model}, None);

// create a site filter; this is an ugly hack for the time being
// alignments.serialize_site_filter returns HBL code as string in 
// which the function `__make_filter` is defined.
 
ExecuteCommands (alignments.serialize_site_filter (
								   rate4site.filter_names[0],
								   ((rate4site.site_patterns[0])["sites"])[0]));

__make_filter ("rate4site.site_filter");

LikelihoodFunction rate4site.site_likelihood = (rate4site.site_filter, rate4site.site_tree);

rate4site.site_model_scaler_name = "rate4site.site_rate_estimate";

rate4site.rate_estimates = {}; 

/**
	 this will store site estimates, which will then be dumped to JSON
*/

parameters.DeclareGlobal (rate4site.site_model_scaler_name, None);

estimators.ApplyExistingEstimates ("rate4site.site_likelihood", rate4site.site_model_mapping, rate4site.alignment_wide_MLES,
									 {"0" : rate4site.site_model_scaler_name} // proportional scaler
									);
					

rate4site.queue = mpi.CreateQueue ({"LikelihoodFunctions": {{"rate4site.site_likelihood"}},
								    "Models" : {{"rate4site.site_model"}},
								    "Headers" : utility.GetListOfLoadedModules (),
								    "Variables" : {{"rate4site.site_model_scaler_name"}}
							 });

/* run the main loop over all unique site pattern combinations */
utility.ForEachPair (rate4site.site_patterns, "_pattern_", "_pattern_info_",
	'
		mpi.QueueJob (rate4site.queue, "rate4site.handle_a_site", {"0" : "rate4site.site_likelihood",
														"1" : alignments.serialize_site_filter
														   ((rate4site.filter_specification[0])["name"],
														   (_pattern_info_["sites"])[0]),
														"2" : _pattern_info_,
														"3" : rate4site.site_model_mapping
														},
														"rate4site.store_results");
	'
);

mpi.QueueComplete (rate4site.queue);

rate4site.stats = (math.GatherDescriptiveStats(utility.Map (utility.Values(utility.Map (rate4site.rate_estimates, "_value_", "_value_[terms.MLE]")), "_value_", "0+_value_")));


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

io.ReportProgressMessageMD ("rate4site", "Stats", "Rate distribution summary");
io.ReportProgressMessageMD ("rate4site", "Stats", "* **Mean**: "  + Format (rate4site.stats["Mean"], 6, 2));
io.ReportProgressMessageMD ("rate4site", "Stats", "* **Median**: "  + Format (rate4site.stats["Median"], 6, 2));
io.ReportProgressMessageMD ("rate4site", "Stats", "* **Std.Dev**: "  + Format (rate4site.stats["Std.Dev"], 6, 2));
io.ReportProgressMessageMD ("rate4site", "Stats", "* **95% Range**: ["  + Format (rate4site.stats["2.5%"], 5,2) + "," + Format (rate4site.stats["97.5%"], 5,2) + "]");
/* Without formatting:
io.ReportProgressMessageMD ("rate4site", "Stats", "* [Mean] "  + rate4site.stats["Mean"]);
io.ReportProgressMessageMD ("rate4site", "Stats", "* [Median] "  + rate4site.stats["Median"]);
io.ReportProgressMessageMD ("rate4site", "Stats", "* [Std.Dev] "  + rate4site.stats["Std.Dev"]);
io.ReportProgressMessageMD ("rate4site", "Stats", "* [95% Range] "  + rate4site.stats["2.5%"] + ":" + rate4site.stats["97.5%"]);
*/

							
io.SpoolJSON ({
				'Relative site rate estimates' : rate4site.rate_estimates, 
				'alignment' : rate4site.alignment_info["file"],
				'analysis' : rate4site.analysis_description
				}
				, rate4site.alignment_info["file"] + ".site-rates.json");

//----------------------------------------------------------------------------------------
// HANDLERS
//----------------------------------------------------------------------------------------

// fit a rate at a single site

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

lfunction rate4site.handle_a_site (lf, filter_data, pattern_info, model_mapping) {
	

    GetString (lfInfo, ^lf,-1);
    ExecuteCommands (filter_data);
    
    __make_filter ((lfInfo["Datafilters"])[0]);
    utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);
    
    if (pattern_info ["is_constant"]) {
    	// the MLE for a constant site is 0; 
    	// only the CI is non-trivial
		^(^"rate4site.site_model_scaler_name") = 0;
    
    } else {
    
		^(^"rate4site.site_model_scaler_name") = 1;
		Optimize (results, ^lf);
	}
	
    return parameters.GetProfileCI (^"rate4site.site_model_scaler_name", lf, 0.95);
}

// handle result processing

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------


lfunction rate4site.store_results (node, result, arguments) {
    pattern_info    = arguments [2];


	if ((^'rate4site.table_output_options')["header"]) {
		
		io.ReportProgressMessageMD ("rate4site", "sites", "Site rate estimates and associated 95% profile likelihood estimates\n"); 

		fprintf (stdout,
			io.FormatTableRow (^'rate4site.table_screen_output',^'rate4site.table_output_options'));
		(^'rate4site.table_output_options')["header"] = FALSE;
	}
	
	
	utility.ForEach (pattern_info["sites"], "_site_index_",
		"
			rate4site.rate_estimates [_site_index_+1] = `&result`;
			result_row = {1,3};
			result_row [0] = '' + (_site_index_ + 1);
			result_row [1] = Format((`&result`)[terms.MLE],6,3);
			result_row [2] = Format((`&result`)[terms.lower_bound],6,3) + ' :'  +Format((`&result`)[terms.upper_bound],6,3);
			fprintf (stdout,
				io.FormatTableRow (result_row,rate4site.table_output_options));
		"
    );

	return rate_statistics;

}



