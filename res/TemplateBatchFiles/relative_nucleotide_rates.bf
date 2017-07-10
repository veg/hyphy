RequireVersion("2.31");
LoadFunctionLibrary("libv3/function-loader.bf");

LoadFunctionLibrary("libv3/models/DNA.bf");
LoadFunctionLibrary("libv3/models/DNA/GTR.bf");


/*------------------------------------------------------------------------------*/

utility.ToggleEnvVariable ("NORMALIZE_SEQUENCE_NAMES", 1);

relative_nuc_rates.analysis_description = {
    "info": "For a fixed **nucleotide** alignment and tree, infer **relative** site specific substitution rates,
    by first optimizing alignment-wide branch lengths, and then inferring a site-specific uniform tree scaler",
    "version": "0.1alpha",
    "reference": "@TBD. Analysis based on Rate4Site method, extended for nucleotides: Pupko, T., Bell, R. E., Mayrose, I., Glaser, F. & Ben-Tal, N. relative_nuc_rates: an algorithmic tool for the identification of functional regions in proteins by surface mapping of evolutionary determinants within their homologues. Bioinformatics 18, S71–S77 (2002).",
    "authors": "Sergei L Kosakovsky Pond and Stephanie J Spielman",
    "contact": "{spond,stephanie.spielman}@temple.edu"
};

io.DisplayAnalysisBanner(relative_nuc_rates.analysis_description);



/***************************************** LOAD DATASET **********************************************************/
SetDialogPrompt ("Specify a nucleotide multiple sequence alignment file");


relative_nuc_rates.alignment_info       = alignments.ReadNucleotideDataSet ("relative_nuc_rates.dataset", None);
relative_nuc_rates.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (relative_nuc_rates.alignment_info[utility.getGlobalValue("terms.json.partitions")], None);
relative_nuc_rates.partition_count = Abs (relative_nuc_rates.partitions_and_trees);

io.CheckAssertion ("relative_nuc_rates.partition_count==1", "This analysis can only handle a single partition");


io.ReportProgressMessageMD ("relative_nuc_rates", "Data", "Input alignment description");
io.ReportProgressMessageMD ("relative_nuc_rates", "Data", "Loaded **" +
                            relative_nuc_rates.alignment_info ["sequences"] + "** sequences, **" +
                            relative_nuc_rates.alignment_info ["sites"] + "** sites, and **" + relative_nuc_rates.partition_count + "** partitions from \`" + relative_nuc_rates.alignment_info ["file"] + "\`");

relative_nuc_rates.filter_specification = alignments.DefineFiltersForPartitions (relative_nuc_rates.partitions_and_trees, "relative_nuc_rates.dataset" , "relative_nuc_rates.filter.", relative_nuc_rates.alignment_info);

relative_nuc_rates.model_generator = "models.DNA.GTR.ModelDescription";





io.ReportProgressMessageMD ("relative_nuc_rates", "overall", "Obtaining alignment-wide branch-length estimates");

relative_nuc_rates.trees = utility.Map (relative_nuc_rates.partitions_and_trees, "_value_", "_value_['tree']"); // value => value['tree']
relative_nuc_rates.filter_names = utility.Map (relative_nuc_rates.filter_specification, "_value_", "_value_['name']"); // value => value['name']
relative_nuc_rates.alignment_wide_MLES = estimators.FitSingleModel_Ext (
                                                          relative_nuc_rates.filter_names,
                                                          relative_nuc_rates.trees, 
                                                          relative_nuc_rates.model_generator,
                                                          None,
                                                          None);
                                                          
estimators.fixSubsetOfEstimates(relative_nuc_rates.alignment_wide_MLES, relative_nuc_rates.alignment_wide_MLES[terms.global]);

io.ReportProgressMessageMD ("relative_nuc_rates", "overall", ">Fitted an alignment-wide model. **Log-L = " + relative_nuc_rates.alignment_wide_MLES ['LogL'] + "**.");


/** 
	Set up the table to display to the screen 
*/


relative_nuc_rates.table_screen_output  = {{"Site", "Rel. rate (MLE)", "95% profile likelihood CI"}};
relative_nuc_rates.table_output_options = {"header" : TRUE, "min-column-width" : 16, "align" : "center"};

relative_nuc_rates.site_patterns = alignments.Extract_site_patterns (relative_nuc_rates.filter_names[0]);

/**
	relative_nuc_rates.site_patterns  maps a unique site pattern (by index) to 
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

relative_nuc_rates.site_model = model.generic.DefineModel(relative_nuc_rates.model_generator,
        "relative_nuc_rates_site_model_instance", {
            "0": parameters.Quote(terms.global),
        },
        relative_nuc_rates.filter_names[0],
        None);
        

relative_nuc_rates.site_model_mapping = {"relative_nuc_rates_site_model_instance" : relative_nuc_rates.site_model};
        
// relative_nuc_rates.site_tree is created from the information in  relative_nuc_rates.trees[0]
// and populated with (the default) model        
        
model.ApplyModelToTree( "relative_nuc_rates.site_tree", relative_nuc_rates.trees[0], {"default" : relative_nuc_rates.site_model}, None);

// create a site filter; this is an ugly hack for the time being
// alignments.serialize_site_filter returns HBL code as string in 
// which the function `__make_filter` is defined.
 
ExecuteCommands (alignments.serialize_site_filter (
								   relative_nuc_rates.filter_names[0],
								   ((relative_nuc_rates.site_patterns[0])["sites"])[0]));

__make_filter ("relative_nuc_rates.site_filter");

LikelihoodFunction relative_nuc_rates.site_likelihood = (relative_nuc_rates.site_filter, relative_nuc_rates.site_tree);

relative_nuc_rates.site_model_scaler_name = "relative_nuc_rates.site_rate_estimate";

relative_nuc_rates.rate_estimates = {}; 

/**
	 this will store site estimates, which will then be dumped to JSON
*/

parameters.DeclareGlobal (relative_nuc_rates.site_model_scaler_name, None);

estimators.ApplyExistingEstimates ("relative_nuc_rates.site_likelihood", relative_nuc_rates.site_model_mapping, relative_nuc_rates.alignment_wide_MLES,
									 {"0" : relative_nuc_rates.site_model_scaler_name} // proportional scaler
									);
					

relative_nuc_rates.queue = mpi.CreateQueue ({"LikelihoodFunctions": {{"relative_nuc_rates.site_likelihood"}},
								    "Models" : {{"relative_nuc_rates.site_model"}},
								    "Headers" : utility.GetListOfLoadedModules (),
								    "Variables" : {{"relative_nuc_rates.site_model_scaler_name"}}
							 });

/* run the main loop over all unique site pattern combinations */
utility.ForEachPair (relative_nuc_rates.site_patterns, "_pattern_", "_pattern_info_",
	'
		mpi.QueueJob (relative_nuc_rates.queue, "relative_nuc_rates.handle_a_site", {"0" : "relative_nuc_rates.site_likelihood",
														"1" : alignments.serialize_site_filter
														   ((relative_nuc_rates.filter_specification[0])["name"],
														   (_pattern_info_["sites"])[0]),
														"2" : _pattern_info_,
														"3" : relative_nuc_rates.site_model_mapping
														},
														"relative_nuc_rates.store_results");
	'
);

mpi.QueueComplete (relative_nuc_rates.queue);

relative_nuc_rates.stats = (math.GatherDescriptiveStats(utility.Map (utility.Values(utility.Map (relative_nuc_rates.rate_estimates, "_value_", "_value_[terms.MLE]")), "_value_", "0+_value_")));


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

io.ReportProgressMessageMD ("relative_nuc_rates", "Stats", "Rate distribution summary");
io.ReportProgressMessageMD ("relative_nuc_rates", "Stats", "* **Mean**: "  + Format (relative_nuc_rates.stats["Mean"], 6, 2));
io.ReportProgressMessageMD ("relative_nuc_rates", "Stats", "* **Median**: "  + Format (relative_nuc_rates.stats["Median"], 6, 2));
io.ReportProgressMessageMD ("relative_nuc_rates", "Stats", "* **Std.Dev**: "  + Format (relative_nuc_rates.stats["Std.Dev"], 6, 2));
io.ReportProgressMessageMD ("relative_nuc_rates", "Stats", "* **95% Range**: ["  + Format (relative_nuc_rates.stats["2.5%"], 5,2) + "," + Format (relative_nuc_rates.stats["97.5%"], 5,2) + "]");
/* Without formatting:
io.ReportProgressMessageMD ("relative_nuc_rates", "Stats", "* [Mean] "  + relative_nuc_rates.stats["Mean"]);
io.ReportProgressMessageMD ("relative_nuc_rates", "Stats", "* [Median] "  + relative_nuc_rates.stats["Median"]);
io.ReportProgressMessageMD ("relative_nuc_rates", "Stats", "* [Std.Dev] "  + relative_nuc_rates.stats["Std.Dev"]);
io.ReportProgressMessageMD ("relative_nuc_rates", "Stats", "* [95% Range] "  + relative_nuc_rates.stats["2.5%"] + ":" + relative_nuc_rates.stats["97.5%"]);
*/



tree_definition   = utility.Map (relative_nuc_rates.partitions_and_trees, "_partition_", '_partition_["tree"]');
io.SpoolJSON ({ "input" : {"filename": relative_nuc_rates.alignment_info["file"],
                          "sequences": relative_nuc_rates.alignment_info["sequences"],
                          "sites":     relative_nuc_rates.alignment_info["sites"],
                          "tree string": (tree_definition[0])["string_with_lengths"]},
                "analysis" : relative_nuc_rates.analysis_description,       
				"Relative site rate estimates" : relative_nuc_rates.rate_estimates, 
				"Global fit": { "tree string": (relative_nuc_rates.alignment_wide_MLES["Trees"])[0],
				        //     "branch lengths": relative_nuc_rates.alignment_wide_MLES["branch length"],
				               "logL": relative_nuc_rates.alignment_wide_MLES["LogL"]}
				},
				relative_nuc_rates.alignment_info["file"] + ".site-rates.json");


//----------------------------------------------------------------------------------------
// HANDLERS
//----------------------------------------------------------------------------------------

// fit a rate at a single site

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

lfunction relative_nuc_rates.handle_a_site (lf, filter_data, pattern_info, model_mapping) {
	

    GetString (lfInfo, ^lf,-1);
    ExecuteCommands (filter_data);
    
    __make_filter ((lfInfo["Datafilters"])[0]);
    utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);
    
    if (pattern_info ["is_constant"]) {
    	// the MLE for a constant site is 0; 
    	// only the CI is non-trivial
		^(^"relative_nuc_rates.site_model_scaler_name") = 0;
    
    } else {
    
		^(^"relative_nuc_rates.site_model_scaler_name") = 1;
		Optimize (results, ^lf);
	}
	
    return parameters.GetProfileCI (^"relative_nuc_rates.site_model_scaler_name", lf, 0.95);
}

// handle result processing

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------


lfunction relative_nuc_rates.store_results (node, result, arguments) {
    pattern_info    = arguments [2];


	if ((^'relative_nuc_rates.table_output_options')["header"]) {
		
		io.ReportProgressMessageMD ("relative_nuc_rates", "sites", "Site rate estimates and associated 95% profile likelihood estimates\n"); 

		fprintf (stdout,
			io.FormatTableRow (^'relative_nuc_rates.table_screen_output',^'relative_nuc_rates.table_output_options'));
		(^'relative_nuc_rates.table_output_options')["header"] = FALSE;
	}
	
	
	utility.ForEach (pattern_info["sites"], "_site_index_",
		"
			relative_nuc_rates.rate_estimates [_site_index_+1] = `&result`;
			result_row = {1,3};
			result_row [0] = '' + (_site_index_ + 1);
			result_row [1] = Format((`&result`)[terms.MLE],6,3);
			result_row [2] = Format((`&result`)[terms.lower_bound],6,3) + ' :'  +Format((`&result`)[terms.upper_bound],6,3);
			fprintf (stdout,
				io.FormatTableRow (result_row,relative_nuc_rates.table_output_options));
		"
    );

	return rate_statistics;

}



