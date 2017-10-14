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


/*------------------------------------------------------------------------------*/

utility.ToggleEnvVariable ("NORMALIZE_SEQUENCE_NAMES", 1);

relative_rates.analysis_description = {
    terms.io.info: "RELrates: Infer relative amino-acid or nucleotide rates from a fixed nucleotide or amino-acid alignment and tree. Relative site-specific substitution rates are
    inferred by first optimizing alignment-wide branch lengths, and then inferring a site-specific uniform tree scaler",
    terms.io.version: "0.1alpha",
    terms.io.reference: "@TBD. Analysis based on Rate4Site method, : Pupko, T., Bell, R. E., Mayrose, I., Glaser, F. & Ben-Tal, N. Rate4Site: an algorithmic tool for the identification of functional regions in proteins by surface mapping of evolutionary determinants within their homologues. Bioinformatics 18, S71–S77 (2002).",
    terms.io.authors: "Sergei L Kosakovsky Pond and Stephanie J Spielman",
    terms.io.contact: "{spond,stephanie.spielman}@temple.edu"
};

io.DisplayAnalysisBanner(relative_rates.analysis_description);

/***************************************** MODEL SELECTION **********************************************************/

relative_rates.protein_type    = "Protein";
relative_rates.nucleotide_type = "Nucleotide";
relative_rates.analysis_type  = io.SelectAnOption ({{relative_rates.protein_type , "Infer relative rates from a protein (amino-acid) alignment"}, {relative_rates.nucleotide_type, "Infer relative rates from a nucleotide alignment"}},
                                                    "Select your analysis type:");


if (relative_rates.analysis_type ==  relative_rates.protein_type) {
    relative_rates.baseline_model  = io.SelectAnOption (models.protein.empirical_models,
                                                    "Select a protein model:");
                                                    
    relative_rates.use_rate_variation = io.SelectAnOption( {{"Gamma", "Use a four-category discrete gamma distribution when optimizing branch lengths."}, 
                                                        {"GDD", "Use a four-category general discrete distribution when optimizing branch lengths."}, 
                                                        {"No", "Do not consider rate variation when optimizing branch lengths."}
                                                        },
                                                        "Optimize branch lengths with rate variation?");
    // "Yes", "No"
    relative_rates.plusF  = io.SelectAnOption ({{"Yes", "Use empirical (+F) amino-acid frequencies ."}, {"No", "Use default amino-acid frequencies."}},                 
                                                        "Use a +F model for initial branch length optimization?");
    if (relative_rates.plusF == "Yes"){
        relative_rates.generators = models.protein.empirical.plusF_generators;
    }
    else {
        relative_rates.generators = models.protein.empirical.default_generators;
    }
}                                                      
else {

    relative_rates.baseline_model  = io.SelectAnOption (models.DNA.models,
                                                    "Select a nucleotide model:");
    relative_rates.use_rate_variation = io.SelectAnOption( {{"Gamma", "Use a four-category discrete gamma distribution when optimizing branch lengths."}, 
                                                        {"GDD", "Use a four-category general discrete distribution when optimizing branch lengths."}, 
                                                        {"No", "Do not consider rate variation when optimizing branch lengths."}
                                                        },
                                                        "Optimize branch lengths with rate variation?");

    relative_rates.generators = models.DNA.generators;
    relative_rates.plusF = "No";
}





function relative_rates.Baseline.ModelDescription(type){
    def = Call( relative_rates.generators[relative_rates.baseline_model], type);
    return def;
}

function relative_rates.Baseline.ModelDescription.withGamma(type){
    def = relative_rates.Baseline.ModelDescription(type);
	def [terms.model.rate_variation] = rate_variation.types.Gamma.factory ({terms.rate_variation.bins : 4});
    return def;
}
function relative_rates.Baseline.ModelDescription.withGDD4(type){
    def = relative_rates.Baseline.ModelDescription(type);
	def [terms.model.rate_variation] = rate_variation.types.GDD.factory ({terms.rate_variation.bins : 4});
    return def;
}


relative_rates.baseline_model_name = relative_rates.baseline_model;
if (relative_rates.plusF == "Yes"){
    relative_rates.baseline_model_name = relative_rates.baseline_model_name + "+F";
}

if (relative_rates.use_rate_variation == "Gamma"){
    relative_rates.baseline_model_name      = relative_rates.baseline_model_name + " with 4 category Gamma rates";
    relative_rates.baseline_model_desc      = "relative_rates.Baseline.ModelDescription.withGamma";
}  
else {
    if (relative_rates.use_rate_variation == "GDD"){
        relative_rates.baseline_model_name      = relative_rates.baseline_model_name + " with 4 category GDD rates";
        relative_rates.baseline_model_desc      = "relative_rates.Baseline.ModelDescription.withGDD4";
    } 
    else {
        relative_rates.baseline_model_name      = relative_rates.baseline_model_name;
        relative_rates.baseline_model_desc      = "relative_rates.Baseline.ModelDescription";
    }
}
/**************************************************************/
   



/*******************************************************************************************************************/


/***************************************** LOAD DATASET **********************************************************/
SetDialogPrompt ("Specify a multiple sequence alignment file");
relative_rates.alignment_info  = alignments.ReadNucleotideDataSet ("relative_rates.dataset", NOne);

name_mapping = relative_rates.alignment_info[utility.getGlobalValue("terms.data.name_mapping")];
if (None == name_mapping) {  
    name_mapping = {};
    utility.ForEach (alignments.GetSequenceNames ("relative_rates.dataset"), "_value_", "`&name_mapping`[_value_] = _value_");
} 
relative_rates.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (relative_rates.alignment_info[utility.getGlobalValue("terms.data.partitions")], name_mapping);
relative_rates.partition_count = Abs (relative_rates.partitions_and_trees);

io.CheckAssertion ("relative_rates.partition_count==1", "This analysis can only handle a single partition");


io.ReportProgressMessageMD ("relative_rates", "Data", "Input alignment description");
io.ReportProgressMessageMD ("relative_rates", "Data", "Loaded **" +
                            relative_rates.alignment_info [terms.data.sequences] + "** sequences, **" +
                            relative_rates.alignment_info [terms.data.sites] + "** sites, and **" + relative_rates.partition_count + "** partitions from \`" + relative_rates.alignment_info [terms.data.file] + "\`");
relative_rates.filter_specification = alignments.DefineFiltersForPartitions (relative_rates.partitions_and_trees, "relative_rates.dataset" , "relative_rates.filter.", relative_rates.alignment_info);
/*******************************************************************************************************************/



/***************************************** INFERENCE **********************************************************/


io.ReportProgressMessageMD ("relative_rates", "overall", "Obtaining alignment-wide branch-length estimates");

relative_rates.trees = utility.Map (relative_rates.partitions_and_trees, "_value_", "_value_[terms.data.tree]"); // value => value['tree']
relative_rates.filter_names = utility.Map (relative_rates.filter_specification, "_value_", "_value_[terms.data.name]"); // value => value['name']
relative_rates.alignment_wide_MLES = estimators.FitSingleModel_Ext (
                                                          relative_rates.filter_names,
                                                          relative_rates.trees, 
                                                          relative_rates.baseline_model_desc,
                                                          None,
                                                          None);
       
       
                                                          
estimators.fixSubsetOfEstimates(relative_rates.alignment_wide_MLES, relative_rates.alignment_wide_MLES[terms.global]);

io.ReportProgressMessageMD ("relative_rates", "overall", ">Fitted an alignment-wide model. **Log-L = " + relative_rates.alignment_wide_MLES [terms.fit.log_likelihood] + "**.");

/** 
	Set up the table to display to the screen 
*/


relative_rates.table_screen_output  = {{"Site", "Rel. rate (MLE)", "95% profile likelihood CI"}};
relative_rates.table_output_options = {terms.table_options.header : TRUE, terms.table_options.minimum_column_width : 16, terms.table_options.align : "center"};

relative_rates.site_patterns = alignments.Extract_site_patterns (relative_rates.filter_names[0]);

// set-up model for site-level fitting in the next couple of lines
relative_rates.site_model = model.generic.DefineModel(relative_rates.baseline_model_desc,
        "relative_rates_site_model_instance", {
            "0": parameters.Quote(terms.global),
        },
        relative_rates.filter_names[0],
        None);
        

relative_rates.site_model_mapping = {"relative_rates_site_model_instance" : relative_rates.site_model};
        
// relative_rates.site_tree is created from the information in  relative_rates.trees[0]
// and populated with (the default) model        
model.ApplyModelToTree( "relative_rates.site_tree", relative_rates.trees[0], {terms.default : relative_rates.site_model}, None);

// create a site filter; this is an ugly hack for the time being
// alignments.serialize_site_filter returns HBL code as string in 
// which the function `__make_filter` is defined.
ExecuteCommands (alignments.serialize_site_filter (
								   relative_rates.filter_names[0],
								   ((relative_rates.site_patterns[0])[terms.data.sites])[0]));

__make_filter ("relative_rates.site_filter");

LikelihoodFunction relative_rates.site_likelihood = (relative_rates.site_filter, relative_rates.site_tree);

relative_rates.site_model_scaler_name = "relative_rates.site_rate_estimate";

relative_rates.rate_estimates = {}; 

/**
	 this will store site estimates, which will then be dumped to JSON
*/

parameters.DeclareGlobal (relative_rates.site_model_scaler_name, None);

estimators.ApplyExistingEstimates ("relative_rates.site_likelihood", relative_rates.site_model_mapping, relative_rates.alignment_wide_MLES,
									 {"0" : relative_rates.site_model_scaler_name} // proportional scaler
									);
					

relative_rates.queue = mpi.CreateQueue ({terms.mpi.LikelihoodFunctions: {{"relative_rates.site_likelihood"}},
								    terms.mpi.Models : {{"relative_rates.site_model"}},
								    terms.mpi.Headers : utility.GetListOfLoadedModules ("libv3/"),
								    terms.mpi.Variables : {{"relative_rates.site_model_scaler_name"}}
							 });

/* run the main loop over all unique site pattern combinations */
utility.ForEachPair (relative_rates.site_patterns, "_pattern_", "_pattern_info_",
	'
		mpi.QueueJob (relative_rates.queue, "relative_rates.handle_a_site", {"0" : "relative_rates.site_likelihood",
														"1" : alignments.serialize_site_filter
														   ((relative_rates.filter_specification[0])[terms.data.name],
														   (_pattern_info_[terms.data.sites])[0]),
														"2" : _pattern_info_,
														"3" : relative_rates.site_model_mapping
														},
														"relative_rates.store_results");
	'
);

mpi.QueueComplete (relative_rates.queue);

relative_rates.site_rates = utility.Map( utility.Values(utility.Map (relative_rates.rate_estimates, "_value_", "_value_[terms.fit.MLE]")), "_value_", "0+_value_");
relative_rates.stats = math.GatherDescriptiveStats(relative_rates.site_rates);

io.ReportProgressMessageMD ("relative_rates", "Stats", "Rate distribution summary");
io.ReportProgressMessageMD ("relative_rates", "Stats", "* **Mean**: "  + Format (relative_rates.stats[terms.math.mean], 6, 2));
io.ReportProgressMessageMD ("relative_rates", "Stats", "* **Median**: "  + Format (relative_rates.stats[terms.math.median], 6, 2));
io.ReportProgressMessageMD ("relative_rates", "Stats", "* **Std.Dev**: "  + Format (relative_rates.stats[terms.math.stddev], 6, 2));
io.ReportProgressMessageMD ("relative_rates", "Stats", "* **95% Range**: ["  + Format (relative_rates.stats[terms.math._2.5], 5,2) + "," + Format (relative_rates.stats[terms.math._97.5], 5,2) + "]");


tree_definition   = utility.Map (relative_rates.partitions_and_trees, "_partition_", '_partition_[terms.data.tree]');
io.SpoolJSON ({ terms.json.input : {terms.json.file: relative_rates.alignment_info[terms.data.file],
                          terms.json.sequences: relative_rates.alignment_info[terms.data.sequences],
                          terms.json.sites:     relative_rates.alignment_info[terms.data.sites],
                          terms.json.tree_string: (tree_definition[0])[terms.trees.newick_with_lengths]},
                terms.json.analysis : relative_rates.analysis_description,       
				terms.json.relative_site_rates : relative_rates.rate_estimates, 
				terms.json.global: {terms.json.model: relative_rates.baseline_model_name,
				               //terms.json.branch_lengths: relative_rates.alignment_wide_MLES[terms.branch_length],
				               terms.json.tree_string: (relative_rates.alignment_wide_MLES[terms.fit.trees])[0],
				               terms.json.log_likelihood: relative_rates.alignment_wide_MLES[terms.fit.log_likelihood]}
				},
				relative_rates.alignment_info[terms.data.file] + ".site-rates.json");


//----------------------------------------------------------------------------------------
// HANDLERS
//----------------------------------------------------------------------------------------

// fit a rate at a single site

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

lfunction relative_rates.handle_a_site (lf, filter_data, pattern_info, model_mapping) {
	

    GetString (lfInfo, ^lf,-1);
    ExecuteCommands (filter_data);
    
    __make_filter ((lfInfo["Datafilters"])[0]);
    utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);
    
    if (pattern_info [utility.getGlobalValue("terms.data.is_constant")]) {
    	// the MLE for a constant site is 0; 
    	// only the CI is non-trivial
		^(utility.getGlobalValue("relative_rates.site_model_scaler_name")) = 0;
    
    } else {
    
		^(utility.getGlobalValue("relative_rates.site_model_scaler_name")) = 1;
		Optimize (results, ^lf);
	}
	
    return parameters.GetProfileCI (utility.getGlobalValue("relative_rates.site_model_scaler_name"), lf, 0.95);
}


// handle result processing

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------


lfunction relative_rates.store_results (node, result, arguments) {
    pattern_info    = arguments [2];


	if ((^'relative_rates.table_output_options')[utility.getGlobalValue("terms.table_options.header")]) {
		
		io.ReportProgressMessageMD ("relative_rates", "sites", "Site rate estimates and associated 95% profile likelihood estimates\n"); 

		fprintf (stdout,
			io.FormatTableRow (^'relative_rates.table_screen_output',^'relative_rates.table_output_options'));
		(^'relative_rates.table_output_options')[utility.getGlobalValue("terms.table_options.header")] = FALSE;
	}
	
	
	utility.ForEach (pattern_info[utility.getGlobalValue("terms.data.sites")], "_site_index_",
		"
			relative_rates.rate_estimates [_site_index_+1] = `&result`;
			result_row = {1,3};
			result_row [0] = '' + (_site_index_ + 1);
			result_row [1] = Format((`&result`)[terms.fit.MLE],6,3);
			result_row [2] = Format((`&result`)[terms.lower_bound],6,3) + ' :'  +Format((`&result`)[terms.upper_bound],6,3);
			fprintf (stdout,
				io.FormatTableRow (result_row,relative_rates.table_output_options));
		"
    );

	return rate_statistics;

}



