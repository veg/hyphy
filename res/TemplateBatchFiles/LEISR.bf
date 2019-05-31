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
    terms.io.info: "LEISR (Likelihood Estimation of Individual Site Rates) infer relative amino-acid or nucleotide rates from a fixed nucleotide or amino-acid alignment and tree, with possibility for partitions. Relative site-specific substitution rates are
    inferred by first optimizing alignment-wide branch lengths, and then inferring a site-specific uniform tree scaler.",
    terms.io.version: "0.4",
    terms.io.reference: "Spielman, S.J. and Kosakovsky Pond, S.L. (2018). Relative evolutionary rate inference in HyPhy with PeerJ 6:e4339. DOI 10.7717/peerj.4339 ; Pupko, T., Bell, R. E., Mayrose, I., Glaser, F. & Ben-Tal, N. Rate4Site: an algorithmic tool for the identification of functional regions in proteins by surface mapping of evolutionary determinants within their homologues. Bioinformatics 18, S71–S77 (2002).",
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
leisr.filter_specification = alignments.DefineFiltersForPartitions (leisr.partitions_and_trees, "leisr.dataset" , "leisr.filter.", leisr.alignment_info);



io.ReportProgressMessageMD ("relative_rates", "Data", "Input alignment description");
io.ReportProgressMessageMD ("relative_rates", "Data", "Loaded **" +
                            leisr.alignment_info [terms.data.sequences] + "** sequences, **" +
                            leisr.alignment_info [terms.data.sites] + "** sites, and **" + leisr.partition_count + "** partitions from \`" + leisr.alignment_info [terms.data.file] + "\`");


/*******************************************************************************************************************/


/***************************************** MODEL SELECTION **********************************************************/

leisr.protein_type    = "Protein";
leisr.nucleotide_type = "Nucleotide";
leisr.analysis_type  = io.SelectAnOption ({{leisr.protein_type , "Infer relative rates from a protein (amino-acid) alignment"}, {leisr.nucleotide_type, "Infer relative rates from a nucleotide alignment"}},
                                                    "Select your analysis type:");


if (leisr.analysis_type ==  leisr.protein_type) {
    leisr.baseline_model  = io.SelectAnOption (models.protein.empirical_models, "Select a protein model:");
    leisr.generators = models.protein.empirical.plusF_generators;
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
                                                          {terms.run_options.retain_lf_object: TRUE});
                                                                  
/* reconstruct ancestral states - */

// SW20180516 - NOT YET SUPPORTED

//leisr.ancestral_data = ancestral.build (leisr.alignment_wide_MLES[terms.likelihood_function], 0, FALSE);
//console.log (leisr.ancestral_data);                                                   
//return 0;

ConstructCategoryMatrix (leisr.site_level_log_likelihoods, ^(leisr.alignment_wide_MLES[terms.likelihood_function]), SITE_LOG_LIKELIHOODS);
DeleteObject (^(leisr.alignment_wide_MLES[terms.likelihood_function]));

estimators.fixSubsetOfEstimates(leisr.alignment_wide_MLES, leisr.alignment_wide_MLES[terms.global]);

io.ReportProgressMessageMD ("relative_rates", "overall", ">Fitted an alignment-wide model. **Log-L = " + leisr.alignment_wide_MLES [terms.fit.log_likelihood] + "**.\n\nTotal tree lengths by partition\n");

utility.ForEachPair (leisr.alignment_wide_MLES[terms.branch_length], "_part_", "_value_", 
'
io.ReportProgressMessageMD ("relative_rates", "overall", "" + (1+_part_) + ". " + Format (+(utility.Map (_value_, "_data_",
    "
        _data_ [terms.fit.MLE]
    "
)),6,3) + " subs/site.");
'
);


/**
	Set up the table to display to the screen
*/


leisr.table_screen_output  = {{"Site", "Partition", "Rel. rate (MLE)", "95% profile likelihood CI", "Log(L) global", "Log(L) local"}};
leisr.table_output_options = {terms.table_options.header : TRUE, 
                              terms.table_options.minimum_column_width : 16, 
                              terms.table_options.align : "center"};

leisr.table_row_report = {{
                                "" + (1+((leisr.filter_specification[leisr.report.partition])[terms.data.coverage])[leisr.report.site]),
                                leisr.report.partition + 1,
                                    Format(leisr.report.row[0],10,6),
                                    Format(leisr.report.row[1],6,2) + " : " + Format(leisr.report.row[2],6,2),
                                    Format (leisr.report.row[3], 6,3),
                                    Format (leisr.report.row[4], 6,3)
                                    }};


leisr.table_headers = {{"MLE", "Relative rate estimate at a site"}
                       {"Lower", "Lower bound of 95% profile likelihood CI"}
                       {"Upper", "Upper bound of 95% profile likelihood CI"}
                       {"LogL global", "Site log likelihood under the global (average rate) model fit"}
                       {"LogL local", "Site log likelihood under the local (site-specific rate) model fit"}};
                       
// set-up model for site-level fitting in the next couple of lines, where rv turned off
leisr.site_model = model.generic.DefineModel("leisr.Baseline.ModelDescription",
        "leisr_site_model_instance", {
            "0": parameters.Quote(terms.global),
        },
        leisr.filter_names[0],
        None);
        
leisr.site_results = {};  
leisr.rate_estimates = {}; // For stats


     
for (leisr.partition_index = 0; leisr.partition_index < leisr.partition_count; leisr.partition_index += 1) {

    leisr.site_patterns = alignments.Extract_site_patterns ((leisr.filter_specification[leisr.partition_index])[utility.getGlobalValue("terms.data.name")]);
    leisr.site_model_mapping = {"leisr_site_model_instance" : leisr.site_model};


    // leisr.site_tree is created from the information in  leisr.trees[leisr.partition_index]
    // and populated with (the default) model
    model.ApplyModelToTree( "leisr.site_tree", leisr.trees[leisr.partition_index], {terms.default : leisr.site_model}, None);

    // create a site filter; this is an ugly hack for the time being
    // alignments.serialize_site_filter returns HBL code as string in
    // which the function `__make_filter` is defined.
    ExecuteCommands (alignments.serialize_site_filter
                                       ((leisr.filter_specification[leisr.partition_index])[utility.getGlobalValue("terms.data.name")],
                                       ((leisr.site_patterns[0])[utility.getGlobalValue("terms.data.sites")])[0],
                   ));
                                       
    __make_filter ("leisr.site_filter");

    LikelihoodFunction leisr.site_likelihood = (leisr.site_filter, leisr.site_tree);

    leisr.site_model_scaler_name = "leisr.site_rate_estimate";
    parameters.DeclareGlobalWithRanges (leisr.site_model_scaler_name, 1, 0, 1e26);
    
    
    /*
    
    SLKP 20171210
    
        leisr.alignment_wide_MLES will have multiple partitions in general, and applying estimates to a site partition 
        will ALWAYS copy values from the 0-index global partition to it. 
         Clearly, this is not the right thing do do. So the solution is to copy the branch lengths from partition 
        i (>0) to the index zero leisr.alignment_wide_MLES
        
        This will invalide the MLE set, but we should be done it with it by now
        If not, it could be deep copied. 
    
    */
    
    if (leisr.partition_index) {
        (leisr.alignment_wide_MLES [terms.branch_length])[0] = (leisr.alignment_wide_MLES [terms.branch_length])[leisr.partition_index];
    }
    
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
                                                               ((leisr.filter_specification[leisr.partition_index])[terms.data.name],
                                                               (_pattern_info_[utility.getGlobalValue("terms.data.sites")])[0]),
                                                            "2": leisr.partition_index,
                                                            "3" : _pattern_info_,
                                                            "4" : leisr.site_model_mapping
                                                            },
                                                            "leisr.store_results");
        '
    );


    mpi.QueueComplete (leisr.queue);

    leisr.partition_matrix = {Abs (leisr.site_results[leisr.partition_index]), 5}; // mle, lower, upper = 3

    utility.ForEachPair (leisr.site_results[leisr.partition_index], "_key_", "_value_",
    '
        for (leisr.index = 0; leisr.index < 5; leisr.index += 1) {
            leisr.partition_matrix [0+_key_][leisr.index] = _value_[leisr.index];
        }
    '
    );

    leisr.site_results[leisr.partition_index] = leisr.partition_matrix;   

}



/* TODO: Update for compatibility with partitioning 
leisr.site_rates = utility.Map( utility.UniqueValues(utility.Map (leisr.site_results, "_value_", "_value_[terms.fit.MLE]")), "_value_", "0+_value_");
leisr.stats = math.GatherDescriptiveStats(leisr.site_rates);

io.ReportProgressMessageMD ("relative_rates", "Stats", "Rate distribution summary");
io.ReportProgressMessageMD ("relative_rates", "Stats", "* **Mean**: "  + Format (leisr.stats[terms.math.mean], 6, 2));
io.ReportProgressMessageMD ("relative_rates", "Stats", "* **Median**: "  + Format (leisr.stats[terms.math.median], 6, 2));
io.ReportProgressMessageMD ("relative_rates", "Stats", "* **Std.Dev**: "  + Format (leisr.stats[terms.math.stddev], 6, 2));
io.ReportProgressMessageMD ("relative_rates", "Stats", "* **95% Range**: ["  + Format (leisr.stats[terms.math._2.5], 5,2) + "," + Format (leisr.stats[terms.math._97.5], 5,2) + "]");
*/


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
                                    terms.json.frequencies: (leisr.alignment_wide_MLES[utility.getGlobalValue("terms.efv_estimate")])["VALUEINDEXORDER"][0]
				                }
				        },
				        terms.json.MLE : {terms.json.headers   : leisr.table_headers,
                                            terms.json.content : leisr.site_results
                                        }                   
 	                   };                                        

for (partition_index = 0; partition_index < leisr.partition_count; partition_index += 1) {
    selection.io.json_store_branch_attribute(leisr.json_content, utility.getGlobalValue ("terms.original_name"), utility.getGlobalValue ("terms.json.node_label"), -1,
                                         partition_index,
                                         leisr.name_mapping);

    selection.io.json_store_branch_attribute(leisr.json_content, utility.getGlobalValue ("leisr.baseline_model_name"), utility.getGlobalValue ("terms.branch_length"), 0,
                                         partition_index,
                                         selection.io.extract_branch_info((leisr.alignment_wide_MLES[utility.getGlobalValue ("terms.branch_length")])[partition_index], "selection.io.branch.length"));

}

selection.io.json_store_key_value_pair (leisr.json_content, None, utility.getGlobalValue("terms.json.partitions"),
                                                         leisr.filter_specification);
                        
                        
                        
io.SpoolJSON (leisr.json_content, leisr.alignment_info[terms.data.file] + ".LEISR.json");


//----------------------------------------------------------------------------------------
// HANDLERS
//----------------------------------------------------------------------------------------

// fit a rate at a single site

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

lfunction leisr.handle_a_site (lf, filter_data, partition_index, pattern_info, model_mapping) {

    GetString (lfInfo, ^lf,-1);
    ExecuteCommands (filter_data);
    __make_filter ((lfInfo["Datafilters"])[0]);
    
    
    
    utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);
    parameters.SetValue (^"leisr.site_model_scaler_name", 1);
    /*
    
     SLKP 20171210
     RESET the global rate value to a default avoid weird optimization issues because of bad starting conditions
     (e.g. if the previous site was saturated)
    */

    
    if (pattern_info [utility.getGlobalValue("terms.data.is_constant")]) {
    	// the MLE for a constant site is 0;
    	// only the CI is non-trivial
		parameters.SetValue (^"leisr.site_model_scaler_name", 0);
		results = {2,1};
		results[1][0] =  estimators.ComputeLF (lf);

    } else {

		^(utility.getGlobalValue("leisr.site_model_scaler_name")) = 1;
		Optimize (results, ^lf);
	}
	
	
    profile = parameters.GetProfileCI (utility.getGlobalValue("leisr.site_model_scaler_name"), lf, 0.95);
    profile[utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];
    
    return profile;
}


// handle result processing

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

function leisr.report.echo (leisr.report.site, leisr.report.partition, leisr.report.row) {
    fprintf (stdout, io.FormatTableRow (leisr.table_row_report,leisr.table_output_options));

}


lfunction leisr.store_results (node, result, arguments) {
    partition_index = arguments [2];
    pattern_info    = arguments [3];


	if ((^'leisr.table_output_options')[utility.getGlobalValue("terms.table_options.header")]) {

		io.ReportProgressMessageMD ("relative_rates", "sites", "Site rate estimates and associated 95% profile likelihood estimates\n");

		fprintf (stdout,
			io.FormatTableRow (^'leisr.table_screen_output',^'leisr.table_output_options'));
		(^'leisr.table_output_options')[utility.getGlobalValue("terms.table_options.header")] = FALSE;
	}
	

    utility.EnsureKey (^"leisr.site_results", partition_index);
	utility.ForEach (pattern_info[utility.getGlobalValue("terms.data.sites")], "_site_index_",
		"
			//leisr.rate_estimates [_site_index_+1] = `&result`;
						
			result_row = {1,5};
			result_row [0] = (`&result`)[terms.fit.MLE];
			result_row [1] = (`&result`)[terms.lower_bound];
			result_row [2] = (`&result`)[terms.upper_bound];
			result_row [3] = leisr.site_level_log_likelihoods[((leisr.filter_specification[`&partition_index`])[terms.data.coverage])[_site_index_]];
 			result_row [4] = (`&result`)[terms.fit.log_likelihood];
 			
 			if (result_row[4] < result_row[3]) {
 			    Export (lfe, ^'leisr.site_likelihood');
 			    console.log (lfe);
 			    console.log ('FUBAR!!!');
 			    console.log (((leisr.filter_specification[`&partition_index`])[terms.data.coverage])[_site_index_]);
 			}
			
            leisr.report.echo (_site_index_, `&partition_index`, result_row);
			
            out_result_row = {1,5};
            out_result_row[0] = (`&result`)[terms.fit.MLE];
            out_result_row[1] = (`&result`)[terms.lower_bound];
            out_result_row[2] = (`&result`)[terms.upper_bound];
            out_result_row[3] = result_row [3];
            out_result_row[4] = (`&result`)[terms.fit.log_likelihood];


		    (leisr.site_results[`&partition_index`])[_site_index_] = out_result_row;

	        /*           		
		    // JSON-related
		    leisr.site_results[_site_index_][0] = (leisr.rate_estimates[_site_index_+1])[terms.fit.MLE];
		    leisr.site_results[_site_index_][1] = (leisr.rate_estimates[_site_index_+1])[terms.lower_bound];
		    leisr.site_results[_site_index_][2] = (leisr.rate_estimates[_site_index_+1])[terms.upper_bound];
		   */
		"
    );

 	//return rate_statistics;
}



lfunction leisr.getIC(logl, params, samples) {
    return -2 * logl + 2 * samples / (samples - params - 1) * params;
}
