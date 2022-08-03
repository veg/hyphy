RequireVersion("2.5.40");

/*------------------------------------------------------------------------------
    Load library files
*/




LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");

LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");

LoadFunctionLibrary("libv3/models/codon/BS_REL.bf");

LoadFunctionLibrary("libv3/convenience/math.bf");

LoadFunctionLibrary("modules/io_functions.ibf");
LoadFunctionLibrary("modules/selection_lib.ibf");


/*------------------------------------------------------------------------------ 
    Display analysis information
*/

meme.analysis_description = {
    terms.io.info: "MEME (Mixed Effects Model of Evolution)
    estimates a site-wise synonymous (&alpha;) and a two-category mixture of non-synonymous
    (&beta;-, with proportion p-, and &beta;+ with proportion [1-p-]) rates, and
    uses a likelihood ratio test to determine if &beta;+ > &alpha; at a site.
    The estimates aggregate information over a proportion of branches at a site,
    so the signal is derived from
    episodic diversification, which is a combination of strength of selection [effect size] and
    the proportion of the tree affected. A subset of branches can be selected
    for testing as well, in which case an additional (nuisance) parameter will be
    inferred -- the non-synonymous rate on branches NOT selected for testing. Multiple partitions within a NEXUS file are also supported
    for recombination - aware analysis. Version 3.0 adds a different format for ancestral state reconstruction, branch-site posterior storage, and site-level heterogeneity testing. 
    ",
    terms.io.version: "3.0",
    terms.io.reference: "Detecting Individual Sites Subject to Episodic Diversifying Selection. _PLoS Genet_ 8(7): e1002764.",
    terms.io.authors: "Sergei L. Kosakovsky Pond, Steven Weaver",
    terms.io.contact: "spond@temple.edu",
    terms.io.requirements: "in-frame codon alignment and a phylogenetic tree"
};


io.DisplayAnalysisBanner(meme.analysis_description);


/*------------------------------------------------------------------------------
    Environment setup
*/

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);


/*------------------------------------------------------------------------------
    Globals
*/

meme.parameter_site_alpha = "Site relative synonymous rate";
meme.parameter_site_omega_minus = "Omega ratio on (tested branches); negative selection or neutral evolution (&omega;- <= 1;)";
meme.parameter_site_beta_minus = "Site relative non-synonymous rate (tested branches); negative selection or neutral evolution (&beta;- <= &alpha;)";
meme.parameter_site_beta_plus = "Site relative non-synonymous rate (tested branches); unconstrained";
meme.parameter_site_mixture_weight = "Beta- category weight";
meme.parameter_site_beta_nuisance = "Site relative non-synonymous rate (untested branches)";
meme.bsER = "Posterior prob omega class by site";


// default cutoff for printing to screen
meme.pvalue = 0.1;
meme.nrate_classes = 2; // there are two rate classes
// The dictionary of results to be written to JSON at the end of the run
meme.json = {
    terms.json.analysis: meme.analysis_description,
    terms.json.fits: {},
    terms.json.timers: {},
    terms.substitutions : {}
};

meme.display_orders =   {terms.original_name: -1,
                        terms.json.nucleotide_gtr: 0,
                        terms.json.global_mg94xrev: 1
                       };

selection.io.startTimer (meme.json [terms.json.timers], "Total time", 0);


/*------------------------------------------------------------------------------
    Key word arguments
*/

KeywordArgument ("code", "Which genetic code should be used", "Universal");
KeywordArgument ("alignment", "An in-frame codon alignment in one of the formats supported by HyPhy");
KeywordArgument ("tree", "A phylogenetic tree (optionally annotated with {})", null, "Please select a tree file for the data:");
KeywordArgument ("branches",  "Branches to test", "All");
KeywordArgument ("pvalue",  "The p-value threshold to use when testing for selection", "0.1");
// One additional KeywordArgument ("output") is called below after namespace meme.


meme.scaler_prefix = "MEME.scaler";

meme.table_headers = {{"&alpha;", "Synonymous substitution rate at a site"}
                     {"&beta;<sup>-</sup>", "Non-synonymous substitution rate at a site for the negative/neutral evolution component"}
                     {"p<sup>-</sup>", "Mixture distribution weight allocated to &beta;<sup>-</sup>; loosely -- the proportion of the tree evolving neutrally or under negative selection"}
                     {"&beta;<sup>+</sup>", "Non-synonymous substitution rate at a site for the positive/neutral evolution component"}
                     {"p<sup>+</sup>", "Mixture distribution weight allocated to &beta;<sup>+</sup>; loosely -- the proportion of the tree evolving neutrally or under positive selection"}
                     {"LRT", "Likelihood ratio test statistic for episodic diversification, i.e., p<sup>+</sup> &gt; 0 <emph>and<emph> &beta;<sup>+</sup> &gt; &alpha;"}
                     {"p-value", "Asymptotic p-value for episodic diversification, i.e., p<sup>+</sup> &gt; 0 <emph>and<emph> &beta;<sup>+</sup> &gt; &alpha;"}
                     {"# branches under selection", "The (very approximate and rough) estimate of how many branches may have been under selection at this site, i.e., had an empirical Bayes factor of 100 or more for the &beta;<sup>+</sup> rate"}
                     {"Total branch length", "The total length of branches contributing to inference at this site, and used to scale dN-dS"}
                     {"MEME LogL", "Site Log-likelihood under the MEME model"}
                     {"FEL LogL", "Site Log-likelihood under the FEL model"},
                     {"Variation p", "Asymptotic p-value for whether or not there is evidence of dN/dS variation across branches"}};


/**
This table is meant for HTML rendering in the results web-app; can use HTML characters, the second column
is 'pop-over' explanation of terms. This is ONLY saved to the JSON file. For Markdown screen output see
the next set of variables.
*/
meme.table_screen_output  = {{"Codon", "Partition", "alpha", "beta+", "p+", "LRT", "Episodic selection detected?", "# branches", "Most common codon substitutions at this site"}};
meme.table_output_options = {terms.table_options.header : TRUE, terms.table_options.minimum_column_width: 12, terms.table_options.align : "center"};


namespace meme {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ("meme");
}

meme.pvalue  = io.PromptUser ("\n>Select the p-value threshold to use when testing for selection",meme.pvalue,0,1,FALSE);

KeywordArgument ("resample",  "[Advanced setting, will result in MUCH SLOWER run time] Perform parametric bootstrap resampling to derive site-level null LRT distributions up to this many replicates per site. Recommended use for small to medium (<30 sequences) datasets", "0");
meme.resample  = io.PromptUser ("\n>[Advanced setting, will result in MUCH SLOWER run time] Perform parametric bootstrap resampling to derive site-level null LRT distributions up to this many replicates per site. Recommended use for small to medium (<30 sequences) datasets",50,0,1000,TRUE);


KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'MEME.json')", meme.codon_data_info [terms.json.json]);
meme.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");

io.ReportProgressMessageMD('MEME',  'selector', 'Branches to include in the MEME analysis');

utility.ForEachPair (meme.selected_branches, "_partition_", "_selection_",
    "_selection_ = utility.Filter (_selection_, '_value_', '_value_ == terms.tree_attributes.test');
     io.ReportProgressMessageMD('MEME',  'selector', 'Selected ' + Abs(_selection_) + ' branches to include in the MEME analysis: \\\`' + Join (', ',utility.Keys(_selection_)) + '\\\`')");


meme.pairwise_counts = genetic_code.ComputePairwiseDifferencesAndExpectedSites(meme.codon_data_info[terms.code], None);

selection.io.startTimer (meme.json [terms.json.timers], "Model fitting",1);

namespace_tag = "meme";

namespace meme {
    doGTR ("meme");
}



estimators.fixSubsetOfEstimates(meme.gtr_results, meme.gtr_results[terms.global]);

// Step 1 - Fit to MG94xREV
namespace meme {
    doPartitionedMG ("meme", FALSE);
}

/*io.ReportProgressMessageMD ("MEME", "codon-refit-prelim", "Improving branch lengths under a full codon model");

estimators.fixSubsetOfEstimates(meme.partitioned_mg_results, meme.partitioned_mg_results[terms.global]);


meme.final_partitioned_mg_results_bl = estimators.FitMGREV (meme.filter_names, meme.trees, meme.codon_data_info [terms.code], {
    terms.run_options.model_type: terms.local,
    terms.run_options.partitioned_omega: meme.selected_branches,
    terms.run_options.retain_lf_object: FALSE,
    terms.run_options.optimization_settings: {
        "OPTIMIZATION_METHOD" : "coordinate-wise",
        "OPTIMIZATION_PRECISION" : Max(meme.partitioned_mg_results[terms.fit.log_likelihood] * 0.00001, 0.1)
    }
}, meme.partitioned_mg_results);


io.ReportProgressMessageMD("MEME", "codon-refit-prelim", "* Log(L) = " + Format(meme.final_partitioned_mg_results_bl[terms.fit.log_likelihood],8,2));
meme.global_dnds = selection.io.extract_global_MLE_re (meme.final_partitioned_mg_results_bl, "^" + terms.parameters.omega_ratio);
utility.ForEach (meme.global_dnds, "_value_", 'io.ReportProgressMessageMD ("MEME", "codon-refit-prelim", "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));');
*/

io.ReportProgressMessageMD ("MEME", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");

KeywordArgument ("full-model", "Perform branch length re-optimization under the full codon model", "Yes");

meme.run_full_mg94 = "Yes" == io.SelectAnOption( {{"Yes", "[Default] Perform branch length re-optimization under the full codon model"}, 
                                         {"No", "Skip branch length re-optimization under the full codon model (faster but less precise)"}},
                                         "Perform branch length re-optimization under the full codon model");


if (meme.run_full_mg94) {

    meme.run_full_mg94 = TRUE;
    
    if (Type (meme.save_intermediate_fits) == "AssociativeList") {
        if (None != meme.save_intermediate_fits[^"terms.data.value"]) {
            if (utility.Has (meme.save_intermediate_fits[^"terms.data.value"], "Full-MG94", "AssociativeList")) {
                meme.final_partitioned_mg_results = (meme.save_intermediate_fits[^"terms.data.value"])["Full-MG94"];
                if (utility.Has (meme.save_intermediate_fits[^"terms.data.value"], "Full-MG94-LF", "String")) {
                    //_PROFILE_NEXUS_LOADS_ = TRUE;
                    ExecuteCommands ((meme.save_intermediate_fits[^"terms.data.value"])["Full-MG94-LF"]);
                    //utility.FinishAndPrintProfile (_NEXUS_PROFILE_DATA_);
                    meme.run_full_mg94 = FALSE;
                
                }
            }        
        }
    }


    if (meme.run_full_mg94) {    
        meme.final_partitioned_mg_results = estimators.FitMGREV (meme.filter_names, meme.trees, meme.codon_data_info [terms.code], {
            terms.run_options.model_type: terms.local,
            terms.run_options.partitioned_omega: meme.selected_branches,
            terms.run_options.retain_lf_object: TRUE,
            terms.run_options.apply_user_constraints: meme.zero_branch_length_constrain,
            terms.run_options.optimization_settings: {
                "OPTIMIZATION_METHOD" : "coordinate-wise"
            }
        }, meme.partitioned_mg_results);

        if (Type (meme.save_intermediate_fits) == "AssociativeList") {
            (meme.save_intermediate_fits[^"terms.data.value"])["Full-MG94"] = meme.final_partitioned_mg_results;        
            Export (lfe, ^meme.final_partitioned_mg_results[^"terms.likelihood_function"]);
            (meme.save_intermediate_fits[^"terms.data.value"])["Full-MG94-LF"] = lfe;
            io.SpoolJSON (meme.save_intermediate_fits[^"terms.data.value"],meme.save_intermediate_fits[^"terms.data.file"]);  
            
        }
    }
} else {
    meme.final_partitioned_mg_results = meme.partitioned_mg_results;
    estimators.RemoveBranchLengthConstraints (meme.final_partitioned_mg_results);
    estimators.RemoveGlobalConstraints (meme.final_partitioned_mg_results);

}
meme.save_intermediate_fits = None;  // clear out memory if used

//meme.final_partitioned_mg_results = meme.partitioned_mg_results;


io.ReportProgressMessageMD("MEME", "codon-refit", "* Log(L) = " + Format(meme.final_partitioned_mg_results[terms.fit.log_likelihood],8,2));
meme.global_dnds = selection.io.extract_global_MLE_re (meme.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio);
utility.ForEach (meme.global_dnds, "_value_", 'io.ReportProgressMessageMD ("MEME", "codon-refit", "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));');


estimators.fixSubsetOfEstimates(meme.final_partitioned_mg_results, meme.final_partitioned_mg_results[terms.global]);



//Store MG94 to JSON
selection.io.json_store_lf_withEFV (meme.json,
                            terms.json.global_mg94xrev,
                            meme.final_partitioned_mg_results[terms.fit.log_likelihood],
                            meme.final_partitioned_mg_results[terms.parameters],
                            meme.sample_size,
                            utility.ArrayToDict (utility.Map (meme.global_dnds, "_value_", "{'key': _value_[terms.description], 'value' : Eval({{_value_ [terms.fit.MLE],1}})}")),
                            (meme.final_partitioned_mg_results[terms.efv_estimate])["VALUEINDEXORDER"][0],
                            meme.display_orders[terms.json.global_mg94xrev]);

utility.ForEachPair (meme.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(meme.json, terms.json.global_mg94xrev, terms.branch_length, meme.display_orders[terms.json.global_mg94xrev],
                                             _key_,
                                             selection.io.extract_branch_info((meme.final_partitioned_mg_results[terms.branch_length])[_key_], "selection.io.branch.length"));');


selection.io.stopTimer (meme.json [terms.json.timers], "Model fitting");

// define the site-level likelihood function

meme.site.background_fel = model.generic.DefineModel("models.codon.MG_REV.ModelDescription",
        "meme.background_fel", {
            "0": parameters.Quote(terms.local),
            "1": meme.codon_data_info[terms.code]
        },
        meme.filter_names,
        None);


meme.alpha = model.generic.GetLocalParameter (meme.site.background_fel, terms.parameters.synonymous_rate);
meme.beta  = model.generic.GetLocalParameter (meme.site.background_fel, terms.parameters.nonsynonymous_rate);

io.CheckAssertion ("None!=meme.alpha && None!=meme.beta", "Could not find expected local synonymous and non-synonymous rate parameters in \`estimators.FitMGREV\`");

meme.site.bsrel =  model.generic.DefineMixtureModel("models.codon.BS_REL_Per_Branch_Mixing.ModelDescription",
        "meme.bsrel", {
            "0": parameters.Quote(terms.local),
            "1": meme.codon_data_info[terms.code],
            "2": parameters.Quote (meme.nrate_classes) // the number of rate classes
        },
        meme.filter_names,
        None);


meme.beta1  = model.generic.GetLocalParameter (meme.site.bsrel , terms.AddCategory (terms.parameters.nonsynonymous_rate,1));
meme.beta2  = model.generic.GetLocalParameter (meme.site.bsrel , terms.AddCategory (terms.parameters.nonsynonymous_rate,2));
meme.branch_mixture = model.generic.GetLocalParameter (meme.site.bsrel , terms.AddCategory (terms.mixture.mixture_aux_weight,1));


io.CheckAssertion ("None!=meme.beta2&&None!=meme.beta1&&None!=meme.branch_mixture", "Could not find expected local rate and mixture parameters for the BS-REL model");

meme.site_model_mapping = {
                           "meme.background_fel" : meme.site.background_fel,
                           "meme.bsrel" : meme.site.bsrel,
                          };


selection.io.startTimer (meme.json [terms.json.timers], "MEME analysis", 2);


// TODO : Why aren't these initially set when the model is first declared?

// FEL parameter declarations
model.generic.AddGlobal (meme.site.background_fel, "meme.site_alpha", meme.parameter_site_alpha);
parameters.DeclareGlobal ("meme.site_alpha", {});

model.generic.AddGlobal (meme.site.background_fel, "meme.site_beta_nuisance", meme.parameter_site_beta_nuisance);
parameters.DeclareGlobal ("meme.site_beta_nuisance", {});


// BSREL mixture model parameter declarations
model.generic.AddGlobal (meme.site.bsrel, "meme.site_alpha", meme.parameter_site_alpha);
model.generic.AddGlobal (meme.site.bsrel, "meme.site_omega_minus", meme.parameter_site_omega_minus);
parameters.DeclareGlobal ("meme.site_omega_minus", {});
parameters.SetRange ("meme.site_omega_minus", terms.range01);

model.generic.AddGlobal (meme.site.bsrel, "meme.site_beta_minus", meme.parameter_site_beta_minus);
parameters.DeclareGlobal ("meme.site_beta_minus", {});
parameters.SetConstraint ("meme.site_beta_minus", "meme.site_alpha * meme.site_omega_minus", "");

model.generic.AddGlobal (meme.site.bsrel, "meme.site_beta_plus", meme.parameter_site_beta_plus);
parameters.DeclareGlobal ("meme.site_beta_plus", {});

model.generic.AddGlobal  (meme.site.bsrel, "meme.site_mixture_weight", meme.parameter_site_mixture_weight);
parameters.DeclareGlobal ("meme.site_mixture_weight", {});
parameters.SetRange ("meme.site_mixture_weight", terms.range_almost_01);

meme.report.count = {{0}};

meme.site.composition.string = "";

meme.report.positive_site = {{"" + (1+((meme.filter_specification[meme.report.partition])[terms.data.coverage])[meme.report.site]),
                                    meme.report.partition + 1,
                                    Format(meme.report.row[0],7,3),
                                    Format(meme.report.row[3],7,3),
                                    Format(meme.report.row[4],7,3),
                                    Format(meme.report.row[5],7,3),
                                    "Yes, p = " + Format(meme.report.row[6],7,4),
                                    Format(meme.report.row[7],0,0),
                                    meme.site.composition.string
}};

meme.site_results = {};
meme.site_LRT = {};

for (meme.partition_index = 0; meme.partition_index < meme.partition_count; meme.partition_index += 1) {
    meme.branch_ebf = {}; 
    meme.report.header_done = FALSE;
    meme.table_output_options[utility.getGlobalValue("terms.table_options.header")] = TRUE;

    meme.model_to_branch_bsrel = { "meme.bsrel" : utility.Filter (meme.selected_branches[meme.partition_index], '_value_', '_value_ == terms.tree_attributes.test'),
                             "meme.background_fel" : utility.Filter (meme.selected_branches[meme.partition_index], '_value_', '_value_ != terms.tree_attributes.test')};


    model.ApplyModelToTree( "meme.site_tree_fel", meme.trees[meme.partition_index], {terms.default : meme.site.background_fel}, None);
    model.ApplyModelToTree( "meme.site_tree_bsrel", meme.trees[meme.partition_index], None, meme.model_to_branch_bsrel);

    meme.site_patterns = alignments.Extract_site_patterns ((meme.filter_specification[meme.partition_index])[utility.getGlobalValue("terms.data.name")]);
    meme.site_count = utility.Array1D ((meme.filter_specification[meme.partition_index])[terms.data.coverage]);
 
    SetParameter (DEFER_CONSTRAINT_APPLICATION, 1, 0);
    
    for (_node_; in; meme.site_tree_fel) {
        _node_class_ = (meme.selected_branches[meme.partition_index])[_node_];
        if (_node_class_ != terms.tree_attributes.test) {
            _beta_scaler = "meme.site_beta_nuisance";
            meme.apply_proportional_site_constraint.fel ("meme.site_tree_bsrel", _node_,
                meme.alpha, meme.beta, "meme.site_alpha", _beta_scaler, (( meme.final_partitioned_mg_results[utility.getGlobalValue("terms.branch_length")])[meme.partition_index])[_node_]);
        } else {
            _beta_scaler = "meme.site_beta_plus";
            meme.apply_proportional_site_constraint.bsrel ("meme.site_tree_bsrel", _node_,
                meme.alpha,  meme.beta1, meme.beta2, meme.branch_mixture, "meme.site_alpha", "meme.site_omega_minus",
                _beta_scaler, "meme.site_mixture_weight", (( meme.final_partitioned_mg_results[utility.getGlobalValue("terms.branch_length")])[meme.partition_index])[_node_]);
        }
        meme.apply_proportional_site_constraint.fel ("meme.site_tree_fel", _node_,
                meme.alpha, meme.beta, "meme.site_alpha", _beta_scaler, (( meme.final_partitioned_mg_results[utility.getGlobalValue("terms.branch_length")])[meme.partition_index])[_node_]);    
    }
    

     SetParameter (DEFER_CONSTRAINT_APPLICATION, 0, 0);



    // create the likelihood function for this site

    ExecuteCommands (alignments.serialize_site_filter
                                       ((meme.filter_specification[meme.partition_index])[utility.getGlobalValue("terms.data.name")],
                                       ((meme.site_patterns[0])[utility.getGlobalValue("terms.data.sites")])[0],
                     ));

    __make_filter ("meme.site_filter");
    LikelihoodFunction meme.site_likelihood = (meme.site_filter, meme.site_tree_fel);



    estimators.ApplyExistingEstimates ("meme.site_likelihood", meme.site_model_mapping, meme.final_partitioned_mg_results,
                                        terms.globals_only);


    __make_filter ("meme.site_filter_bsrel");
    LikelihoodFunction meme.site_likelihood_bsrel = (meme.site_filter_bsrel, meme.site_tree_bsrel);


    estimators.ApplyExistingEstimates ("meme.site_likelihood_bsrel", meme.site_model_mapping, meme.final_partitioned_mg_results,
                                        terms.globals_only);
                                        
                                        
    meme.queue = mpi.CreateQueue ({terms.mpi.LikelihoodFunctions: {{"meme.site_likelihood","meme.site_likelihood_bsrel"}},
                                   terms.mpi.Models : {{"meme.site.background_fel","meme.site.bsrel"}},
                                   terms.mpi.Headers : utility.GetListOfLoadedModules ("libv3/"),
                                   terms.mpi.Variables : {{"meme.selected_branches","meme.branch_mixture","meme.pairwise_counts","meme.codon_data_info","meme.resample"}},
                                   terms.mpi.Functions : {{"meme.compute_branch_EBF"}}
                                 });


    meme.pattern_count_all  = utility.Array1D (meme.site_patterns);
    meme.pattern_count_this = 0;

    /* run the main loop over all unique site pattern combinations */
    utility.ForEachPair (meme.site_patterns, "_pattern_", "_pattern_info_",
        '
            meme.pattern_count_this += 1;
            io.ReportProgressBar("", "Working on site pattern " + (meme.pattern_count_this) + "/" +  meme.pattern_count_all + " in partition " + (1+meme.partition_index));
            if (_pattern_info_[utility.getGlobalValue("terms.data.is_constant")]) {
                meme.store_results (-1,None,{"0" : "meme.site_likelihood",
                                             "1" : "meme.site_likelihood_bsrel",
                                             "2" : None,
                                             "3" : meme.partition_index,
                                             "4" : _pattern_info_,
                                             "5" : meme.site_model_mapping,
                                             "6" : 0
                                     });
            } else {
                mpi.QueueJob (meme.queue, "meme.handle_a_site", {"0" : "meme.site_likelihood",
                                                                 "1" : "meme.site_likelihood_bsrel",
                                                                 "2" : alignments.serialize_site_filter
                                                                   ((meme.filter_specification[meme.partition_index])[utility.getGlobalValue("terms.data.name")],
                                                                   (_pattern_info_[utility.getGlobalValue("terms.data.sites")])[0]),
                                                                 "3" : meme.partition_index,
                                                                 "4" : _pattern_info_,
                                                                 "5" : meme.site_model_mapping,
                                                                 "6" : 0
                                                                    },
                                                                    "meme.store_results");
            }
            pattern_count_all
        '
    );
    
    io.ClearProgressBar();

    mpi.QueueComplete (meme.queue);
    meme.partition_matrix = {Abs (meme.site_results[meme.partition_index]), Rows (meme.table_headers)};
    
    for (_key_,_value_; in; meme.site_results[meme.partition_index]) {
        for (meme.index = 0; meme.index < Rows (meme.table_headers); meme.index += 1) {
            meme.partition_matrix [0+_key_][meme.index] = _value_[meme.index];
        }
    }

    meme.site_results[meme.partition_index] = meme.partition_matrix;
    meme.branch_ebf_matrix = {};
    for (meme.branch,meme.post; in; meme.branch_ebf) {
        meme.post_matrix = {2, meme.site_count};
        for (meme.r,meme.v; in; meme.post) {
            if (Type (meme.v) == "Matrix") {
                meme.post_matrix[0][+meme.r] = meme.v[0];
                meme.post_matrix[1][+meme.r] = meme.v[1];
            }
        }
        meme.branch_ebf_matrix [meme.branch] = meme.post_matrix;
    }
    selection.io.json_store_branch_attribute(meme.json, meme.bsER, terms.json.branch_annotations, -1,
                                     meme.partition_index,
                                     meme.branch_ebf_matrix);

    

}

meme.json [terms.json.MLE ] = {terms.json.headers   : meme.table_headers,
                               terms.json.content : meme.site_results };


io.ReportProgressMessageMD ("MEME", "results", "** Found _" + meme.report.count[0] + "_ sites under episodic diversifying positive selection at p <= " + meme.pvalue + "**");

selection.io.stopTimer (meme.json [terms.json.timers], "Total time");
selection.io.stopTimer (meme.json [terms.json.timers], "MEME analysis");

if (meme.resample)  {
    (meme.json [terms.json.MLE ])[terms.LRT] = meme.site_LRT;
}

GetString (_hpv,HYPHY_VERSION,0);
meme.json[terms.json.runtime] = _hpv;


io.SpoolJSON (meme.json, meme.codon_data_info[terms.json.json]);

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

function meme.apply_proportional_site_constraint.fel (tree_name, node_name, alpha_parameter, beta_parameter, alpha_factor, beta_factor, branch_length) {

    meme.branch_length = (branch_length[terms.parameters.synonymous_rate])[terms.fit.MLE];

    node_name = tree_name + "." + node_name;

    ExecuteCommands ("
        `node_name`.`alpha_parameter` :< 1e10;
        `node_name`.`beta_parameter` :< 1e10;
        `node_name`.`alpha_parameter` := (`alpha_factor`) * meme.branch_length__;
        `node_name`.`beta_parameter`  := (`beta_factor`)  * meme.branch_length__;
    ");
}
//----------------------------------------------------------------------------------------

function meme.apply_proportional_site_constraint.bsrel (tree_name, node_name, alpha_parameter, beta_parameter, beta_parameter2, mixture_parameter, alpha_factor, omega_factor, beta_factor, mixture_global, branch_length) {

    meme.branch_length = (branch_length[terms.parameters.synonymous_rate])[terms.fit.MLE];

    node_name = tree_name + "." + node_name;

    ExecuteCommands ("
        `node_name`.`alpha_parameter` :< 1e10;
        `node_name`.`beta_parameter2` :< 1e10;
        `node_name`.`beta_parameter` :< 1e10;
        `node_name`.`alpha_parameter` := (`alpha_factor`) * meme.branch_length__;
        `node_name`.`beta_parameter`  := (`omega_factor`)  * `node_name`.`alpha_parameter`;
        `node_name`.`beta_parameter2`  := (`beta_factor`)  * meme.branch_length__;
        `node_name`.`mixture_parameter`  := `mixture_global`;
    ");
}

//----------------------------------------------------------------------------------------
lfunction meme.compute_branch_EBF (lf_id, tree_name, branch_name, baseline) {
// TODO: figure out why LFCompute fails if this is run as an `lfunction`

    parameter_name = "`tree_name`.`branch_name`." + ^"meme.branch_mixture";
    ^parameter_name = 1;

    LFCompute (^lf_id,LOGL0);

    utility.ExecuteInGlobalNamespace (parameter_name + ":= meme.site_mixture_weight");

    if (^"meme.site_mixture_weight" != 1 && ^"meme.site_mixture_weight" != 0) {
        _priorOdds = (1-^"meme.site_mixture_weight")/^"meme.site_mixture_weight";
    } else {
        _priorOdds = 0;
    }

    normalizer  = -Max (LOGL0,baseline);


    p1 = Exp(LOGL0+normalizer) * ^"meme.site_mixture_weight";
    p2 = (Exp(baseline+normalizer) - p1);

    _posteriorProb = {{p1,p2}};

    _posteriorProb = _posteriorProb * (1/(+_posteriorProb));
    if ( _priorOdds != 0) {
        eBF = _posteriorProb[1] / (1 - _posteriorProb[1]) / _priorOdds;
    } else {
        eBF = 1;
    }
    return {utility.getGlobalValue("terms.empirical_bayes_factor") : eBF__, 
            utility.getGlobalValue("terms.posterior") : _posteriorProb__};
}

//----------------------------------------------------------------------------------------
lfunction meme.handle_a_site (lf_fel, lf_bsrel, filter_data, partition_index, pattern_info, model_mapping, sim_mode) {

    //console.log (pattern_info);
    //#profile START;
    GetString   (lfInfo, ^lf_fel,-1);   

    //utility.SetEnvVariable ("VERBOSITY_LEVEL", 100);

    //TODO Datafilters hardcode, Trees hardcode.
    ExecuteCommands (filter_data);
    __make_filter ((lfInfo["Datafilters"])[0]);

    GetString (lfInfo, ^lf_bsrel,-1);
    __make_filter ((lfInfo["Datafilters"])[0]);

    bsrel_tree_id = (lfInfo["Trees"])[0];

    utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);
    utility.SetEnvVariable ("ASSUME_REVERSIBLE_MODELS", TRUE);

    ^"meme.site_alpha" = 1;
    ^"meme.site_beta_plus"  = 1;
    ^"meme.site_beta_nuisance"  = 1;

    //console.log ("Optimizing FEL for pattern " + pattern_info);
    //io.SpoolLF (lf_fel, "/tmp/meme.debug" + ^"MPI_NODE_ID", "FEL");
    Optimize (results, ^lf_fel
        , {"OPTIMIZATION_METHOD" : "nedler-mead", OPTIMIZATION_PRECISION: 1e-4}
    );

    fel = estimators.ExtractMLEsOptions (lf_fel, model_mapping, {^"terms.globals_only" : TRUE});
    fel[utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];


    
     if ( ^"meme.site_alpha" <  ^"meme.site_beta_plus") { // FEL returns a +-ve site
        ^"meme.site_mixture_weight" = 0.25;
         ^"meme.site_omega_minus" = 0.25;
     } else { // FEL returns a negative site
        ^"meme.site_mixture_weight" = 0.75;
         if (^"meme.site_alpha" > 0) {
             ^"meme.site_omega_minus" = ^"meme.site_beta_plus" / ^"meme.site_alpha";
        } else {
            ^"meme.site_omega_minus" = 1;
        }
         ^"meme.site_beta_plus" =  ^"meme.site_alpha" * 1.5; 
     }
     
     
    //io.SpoolLF (lf_bsrel, "/tmp/meme.debug", "MEME");
    
    initial_guess_grid = {
                "0" : {
                
                    "meme.site_alpha" : ^"meme.site_alpha",
                    "meme.site_beta_nuisance" : ^"meme.site_beta_nuisance",
                    "meme.site_beta_plus": ^"meme.site_beta_plus",
                    "meme.site_omega_minus": ^"meme.site_omega_minus",
                    "meme.site_mixture_weight": ^"meme.site_mixture_weight"
                   
                },
                "1" : {
                    "meme.site_alpha" : ^"meme.site_alpha",
                    "meme.site_beta_nuisance" : ^"meme.site_beta_nuisance",
                    "meme.site_beta_plus": ^"meme.site_beta_plus" * 2,
                    "meme.site_omega_minus": 0.5,
                    "meme.site_mixture_weight": 0.5                
                },
                "2" : {
                    "meme.site_alpha" : ^"meme.site_alpha",
                    "meme.site_beta_nuisance" : ^"meme.site_beta_nuisance",
                    "meme.site_beta_plus": ^"meme.site_beta_plus" * 4,
                    "meme.site_omega_minus": 0.25,
                    "meme.site_mixture_weight": 0.25                
                },
                "3" : {
                    "meme.site_alpha" : ^"meme.site_alpha",
                    "meme.site_beta_nuisance" : ^"meme.site_beta_nuisance",
                    "meme.site_beta_plus": ^"meme.site_beta_plus",
                    "meme.site_omega_minus": 0.5,
                    "meme.site_mixture_weight": 0.5                
                },
                "4" : {
                    "meme.site_alpha" : ^"meme.site_alpha",
                    "meme.site_beta_nuisance" : ^"meme.site_beta_nuisance",
                    "meme.site_beta_plus": ^"meme.site_beta_plus",
                    "meme.site_omega_minus": 0.75,
                    "meme.site_mixture_weight": 0.8                
                },
                "5" : {
                    "meme.site_alpha" : ^"meme.site_alpha",
                    "meme.site_beta_nuisance" : ^"meme.site_beta_nuisance",
                    "meme.site_beta_plus": ^"meme.site_beta_plus" * 8,
                    "meme.site_omega_minus": 0.5,
                    "meme.site_mixture_weight": 0.8                
                },
                "6" : {
                    "meme.site_alpha" : ^"meme.site_alpha",
                    "meme.site_beta_nuisance" : ^"meme.site_beta_nuisance",
                    "meme.site_beta_plus": ^"meme.site_beta_plus",
                    "meme.site_omega_minus": 0,
                    "meme.site_mixture_weight": 0.01              
                },
                "7" : {
                    "meme.site_alpha" : ^"meme.site_alpha",
                    "meme.site_beta_nuisance" : ^"meme.site_beta_nuisance",
                    "meme.site_beta_plus": ^"meme.site_beta_plus",
                    "meme.site_omega_minus": ^"meme.site_omega_minus",
                    "meme.site_mixture_weight": 1.0              
                }
            
             };
             
    //before_opt = {"alpha" : ^"meme.site_alpha", "other" : initial_guess_grid};
    
    
                      
    ^"meme.site_alpha" = ^"meme.site_alpha";   
    // SLKP 20201028 : without this, there are occasional initialization issues with 
    // the likelihood function below
                      
    Optimize (results, ^lf_bsrel, {
            "OPTIMIZATION_METHOD" : "nedler-mead",
            "OPTIMIZATION_START_GRID" : initial_guess_grid
        });
        
    /*after_opt = {
                    "alpha" : Eval("meme.site_alpha"), 
                    "meme.site_beta_plus": Eval("meme.site_beta_plus"),
                    "meme.site_omega_minus": Eval("meme.site_omega_minus"),
                    "meme.site_mixture_weight": Eval("meme.site_mixture_weight")
                };*/
        
    
    //alternative = estimators.ExtractMLEs (lf_bsrel, model_mapping);
    if (sim_mode) {
        lrt = 2*results[1][0];    
    } else {
        alternative = estimators.ExtractMLEsOptions (lf_bsrel, model_mapping, {^"terms.globals_only" : TRUE});
        alternative [utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];
        ancestral_info = ancestral.build (lf_bsrel,0,FALSE);
        //TODO
        branch_substitution_information = (ancestral.ComputeSubstitutionBySite (ancestral_info,0,None))[^"terms.substitutions"];
        compressed_substitution_info = (ancestral.ComputeCompressedSubstitutions (ancestral_info))[0];
        DeleteObject (ancestral_info);
        branch_ebf       = {};
        branch_posterior = {};
    }

 
    if (^"meme.site_beta_plus" > ^"meme.site_alpha" && ^"meme.site_mixture_weight" < 0.999999) {
        if (!sim_mode) { 
            LFCompute (^lf_bsrel,LF_START_COMPUTE);
            LFCompute (^lf_bsrel,baseline);
     

            for (_node_name_; in; ^bsrel_tree_id) {
                if (((^"meme.selected_branches") [partition_index])[_node_name_]  == utility.getGlobalValue("terms.tree_attributes.test")) {
                    _node_name_res_ = meme.compute_branch_EBF (lf_bsrel, bsrel_tree_id, _node_name_, baseline);
                    branch_ebf[_node_name_] = _node_name_res_[utility.getGlobalValue("terms.empirical_bayes_factor")];
                    branch_posterior[_node_name_] = _node_name_res_[utility.getGlobalValue("terms.posterior")];
                }
            }
        
 
            LFCompute (^lf_bsrel,LF_DONE_COMPUTE);
        }

        ^"meme.site_beta_plus" := ^"meme.site_alpha";
        Optimize (results, ^lf_bsrel, {"OPTIMIZATION_METHOD" : "nedler-mead"});
        //io.SpoolLF (lf_bsrel, "/tmp/meme.debug", "MEME-null");

        //Null = estimators.ExtractMLEs (lf_bsrel, model_mapping);

        if (sim_mode) {
            return lrt - 2*results[1][0];
        } else { 
            Null = estimators.ExtractMLEsOptions (lf_bsrel, model_mapping, {^"terms.globals_only" : TRUE});
            Null [utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];
            if (^"meme.resample") {
                N = ^"meme.resample";
                sims = {};
                GetDataInfo (fi, ^((lfInfo["Datafilters"])[0]), "PARAMETERS");
                for (i = 0; i < N; i+=1) {        
                    DataSet null_sim = SimulateDataSet (^lf_bsrel);
                    DataSetFilter null_filter = CreateFilter (null_sim,3,,,fi["EXCLUSIONS"]);
                    sims + alignments.serialize_site_filter (&null_filter, 0);
                }
                null_LRT = {N,1};
                for (i = 0; i < N; i+=1) {   
                   is = meme.handle_a_site (lf_fel, lf_bsrel, sims[i], partition_index, pattern_info, model_mapping, TRUE);
                   null_LRT[i] = is;
                }
                
                return {"fel" : fel,
                    utility.getGlobalValue("terms.alternative") : alternative,
                    utility.getGlobalValue("terms.posterior") : branch_posterior,
                    utility.getGlobalValue("terms.empirical_bayes_factor") : branch_ebf,
                    utility.getGlobalValue("terms.branch_selection_attributes") : branch_substitution_information, //TODO: keep this attr?
                    utility.getGlobalValue("terms.Null"): Null,
                    utility.getGlobalValue("terms.substitutions") : compressed_substitution_info,
                    utility.getGlobalValue("terms.simulated"): null_LRT
                };
            } 
        }

    } else {
        if (!sim_mode) {
            Null = alternative;
        
            for (_node_name_; in; ^bsrel_tree_id) {
                if (((^"meme.selected_branches") [partition_index])[_node_name_]  == utility.getGlobalValue("terms.tree_attributes.test")) {
                    branch_ebf[_node_name_] = 1.0;
                    branch_posterior[_node_name_] = 0.0;
                }
            }
        } else {
            return 0;
        }        
    }

    //#profile collect;
    //utility.FinishAndPrintProfile (collect);

    return {"fel" : fel,
            utility.getGlobalValue("terms.alternative") : alternative,
            utility.getGlobalValue("terms.posterior") : branch_posterior,
            utility.getGlobalValue("terms.empirical_bayes_factor") : branch_ebf,
            utility.getGlobalValue("terms.branch_selection_attributes") : branch_substitution_information, //TODO: keep this attr?
            utility.getGlobalValue("terms.substitutions") : compressed_substitution_info,
            utility.getGlobalValue("terms.Null"): Null};
}

/* echo to screen calls */

//----------------------------------------------------------------------------------------
function meme.report.echo (meme.report.site, meme.report.partition, meme.report.row) {

    meme.print_row = None;
    if (meme.report.row [6] <= meme.pvalue) {
        meme.print_row = meme.report.positive_site;
        meme.report.count[0] += 1;
    }

     if (None != meme.print_row) {
            io.ClearProgressBar();
            if (!meme.report.header_done) {
                io.ReportProgressMessageMD("MEME", "" + meme.report.partition, "For partition " + (meme.report.partition+1) + " these sites are significant at p <=" + meme.pvalue + "\n");
                fprintf (stdout,
                    io.FormatTableRow (meme.table_screen_output,meme.table_output_options));
                meme.report.header_done = TRUE;
                meme.table_output_options[utility.getGlobalValue("terms.table_options.header")] = FALSE;
            }

            fprintf (stdout,
                io.FormatTableRow (meme.print_row,meme.table_output_options));
        }

}

//----------------------------------------------------------------------------------------

lfunction meme.store_results (node, result, arguments) {

    sub_profile = result[utility.getGlobalValue("terms.branch_selection_attributes")];
    compressed_subs = result[utility.getGlobalValue("terms.substitutions")];
    
/**

{
 "AGA":{
   "AGG":{
     "0":"Node1"
    }
  },
 "AGG":{
   "AAA":{
     "0":"Node3",
     "1":"HORSE"
    },
   "ACG":{
     "0":"Node12"
    },
   "GCA":{
     "0":"CAT"
    }
  }
}

**/

    all_posterior = None;
    has_lrt = FALSE;
    
    if (None != sub_profile) {
        total_sub_count = 0;
    
        sub_by_count = {};
        for (i, v; in; sub_profile) {
            for (j, b; in; v) {
                c = Abs (b);
                total_sub_count += c;
                if ((sub_by_count/c)==0) {
                    sub_by_count [c] = {};
                }
                
                pair_key = {{'','','','>','','',''}};
                for (i2 = 0; i2 < 3; i2+=1) {
                    if (i[i2] != j[i2]) {
                        pair_key[i2] = i[i2];
                        pair_key[i2+4] = j[i2];
                    } else {
                        pair_key[i2] = i[i2] && 0;
                        pair_key[i2+4] = j[i2] && 0;
                   
                    }
                }     
                pair_key = Join ("", pair_key);  
                
                sub_by_count[c] + pair_key;
            }
        }
        
        sorted_subs = {Abs (sub_by_count), 1};
        j = 0;
        for (i, v; in; sub_by_count) {
            sorted_subs[j] = -(+i);
            j += 1; 
        }
        sub_profile = {};
        for (i; in; sorted_subs % 0) {
            sub_profile + ("[" + (-i) + "]" + Join (",",sub_by_count [-i]));
        }
        ^"meme.site.composition.string" = Join ("|", sub_profile);
    }

    partition_index = arguments [3];
    pattern_info    = arguments [4];

    result_row          = { { 0, // alpha 0 
                          0, // beta- 1
                          1, // weight- 2
                          0, // beta + 3
                          0, // weight + 4
                          0, // LRT 5
                          1, // p-value, 6
                          0, // branch count 7
                          0,  // total branch length of tested branches 8
                          0, // Log L | FEL 9 
                          0, // Log L | MEME 10
                          1  // LRT MEME vs FEL
                          } };

      //console.log ( estimators.GetGlobalMLE (result["alternative"], ^"meme.parameter_site_mixture_weight"));

    if (None != result) { // not a constant site

        lrt = {utility.getGlobalValue("terms.LRT") : 2*((result[utility.getGlobalValue("terms.alternative")])[utility.getGlobalValue("terms.fit.log_likelihood")]-(result[utility.getGlobalValue("terms.Null")])[utility.getGlobalValue("terms.fit.log_likelihood")])};
        lrt [utility.getGlobalValue("terms.p_value")] = 2/3-2/3*(0.45*CChi2(lrt[utility.getGlobalValue("terms.LRT")],1)+0.55*CChi2(lrt[utility.getGlobalValue("terms.LRT")],2));
        
        has_lrt = result / ^"terms.simulated";
        
        if (has_lrt) {
           pv = +((result[^"terms.simulated"])["_MATRIX_ELEMENT_VALUE_>=" + lrt [utility.getGlobalValue("terms.LRT")]]);
           lrt [utility.getGlobalValue("terms.p_value")] = (pv+1)/(1+^"meme.resample");
        }     


        result_row [0] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], utility.getGlobalValue("meme.parameter_site_alpha"));
        result_row [1] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], utility.getGlobalValue("meme.parameter_site_omega_minus")) * result_row[0];
        result_row [2] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], utility.getGlobalValue("meme.parameter_site_mixture_weight"));
        result_row [3] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], utility.getGlobalValue("meme.parameter_site_beta_plus"));
        result_row [4] = 1-result_row [2];
        result_row [5] = lrt [utility.getGlobalValue("terms.LRT")];
        result_row [6] = lrt [utility.getGlobalValue("terms.p_value")];
        result_row [9] = (result["fel"])[utility.getGlobalValue("terms.fit.log_likelihood")];
        result_row [10] = (result[utility.getGlobalValue("terms.alternative")])[utility.getGlobalValue("terms.fit.log_likelihood")];
        result_row [11] = 1-CChi2 (2*(result_row [10]-result_row [9]),2);

        all_ebf = result[utility.getGlobalValue("terms.empirical_bayes_factor")];
        all_posterior = result[utility.getGlobalValue("terms.posterior")];

        filtered_ebf = utility.Filter (utility.Filter (all_ebf, "_value_", "None!=_value_"), "_value_", "_value_>=100");

        if(None != filtered_ebf) {
            result_row [7] = utility.Array1D(filtered_ebf);
        } else {
            result_row [7] = 0;
        }

        sum = 0;
        /*
        alternative_lengths = ((result[utility.getGlobalValue("terms.alternative")])[utility.getGlobalValue("terms.branch_length")])[0];
        for (_node_; in; ^"meme.site_tree_fel") {
             _node_class_ = ((^"meme.selected_branches")[partition_index])[_node_];
             if (_node_class_ == utility.getGlobalValue("terms.tree_attributes.test")) {
                sum += ((`&alternative_lengths`)[_node_])[utility.getGlobalValue("terms.json.MLE")];
             }
        }
        */

        result_row [8] = sum;
    } else {
        all_ebf = None;
    }

    utility.EnsureKey (^"meme.site_results", partition_index);
    utility.EnsureKey (((^"meme.json")[^"terms.substitutions"]), partition_index);
    if (has_lrt) {
        utility.EnsureKey (^"meme.site_LRT", partition_index);
    }
    
  
    sites_mapping_to_pattern = pattern_info[utility.getGlobalValue("terms.data.sites")];
    sites_mapping_to_pattern.count = utility.Array1D (sites_mapping_to_pattern);
    

    for (i = 0; i < sites_mapping_to_pattern.count; i+=1) {
        site_index = sites_mapping_to_pattern[i];
        ((^"meme.site_results")[partition_index])[site_index] = result_row;
        (((^"meme.json")[^"terms.substitutions"])[partition_index])[site_index] = compressed_subs;
        if (has_lrt) {
            ((^"meme.site_LRT")[partition_index])[site_index] = result[^"terms.simulated"];
        }
        if (None != all_posterior) {
            for (branch, posterior; in; all_posterior) {
                utility.EnsureKey (^"meme.branch_ebf", branch);
                ((^"meme.branch_ebf")[branch])[site_index] = posterior;
            }
        } 
        meme.report.echo (site_index, partition_index, result_row);
   }



    //assert (0);
}

