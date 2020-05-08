RequireVersion("2.4.0");

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
    for recombination - aware analysis.
    ",
    terms.io.version: "2.1.2",
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

// default cutoff for printing to screen
meme.pvalue = 0.1;
meme.nrate_classes = 2; // there are two rate classes
// The dictionary of results to be written to JSON at the end of the run
meme.json = {
    terms.json.analysis: meme.analysis_description,
    terms.json.fits: {},
    terms.json.timers: {},
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

meme.table_headers = {{"alpha;", "Synonymous substitution rate at a site"}
                     {"&beta;<sup>-</sup>", "Non-synonymous substitution rate at a site for the negative/neutral evolution component"}
                     {"p<sup>-</sup>", "Mixture distribution weight allocated to &beta;<sup>-</sup>; loosely -- the proportion of the tree evolving neutrally or under negative selection"}
                     {"&beta;<sup>+</sup>", "Non-synonymous substitution rate at a site for the positive/neutral evolution component"}
                     {"p<sup>+</sup>", "Mixture distribution weight allocated to &beta;<sup>+</sup>; loosely -- the proportion of the tree evolving neutrally or under positive selection"}
                     {"LRT", "Likelihood ratio test statistic for episodic diversification, i.e., p<sup>+</sup> &gt; 0 <emph>and<emph> &beta;<sup>+</sup> &gt; &alpha;"}
                     {"p-value", "Asymptotic p-value for episodic diversification, i.e., p<sup>+</sup> &gt; 0 <emph>and<emph> &beta;<sup>+</sup> &gt; &alpha;"}
                     {"# branches under selection", "The (very approximate and rough) estimate of how many branches may have been under selection at this site, i.e., had an empirical Bayes factor of 100 or more for the &beta;<sup>+</sup> rate"}
                     {"Total branch length", "The total length of branches contributing to inference at this site, and used to scale dN-dS"}
                     {"MEME LogL", "Site Log-likelihood under the MEME model"}
                     {"FEL LogL", "Site Log-likelihood under the FEL model"}};


/**
This table is meant for HTML rendering in the results web-app; can use HTML characters, the second column
is 'pop-over' explanation of terms. This is ONLY saved to the JSON file. For Markdown screen output see
the next set of variables.
*/
meme.table_screen_output  = {{"Codon", "Partition", "alpha", "beta+", "p+", "LRT", "Episodic selection detected?", "# branches"}};
meme.table_output_options = {terms.table_options.header : TRUE, terms.table_options.minimum_column_width: 12, terms.table_options.align : "center"};


namespace meme {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ("meme");
}

meme.pvalue  = io.PromptUser ("\n>Select the p-value threshold to use when testing for selection",meme.pvalue,0,1,FALSE);

KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'MEME.json')", meme.codon_data_info [terms.json.json]);
meme.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");

io.ReportProgressMessageMD('MEME',  'selector', 'Branches to include in the MEME analysis');

utility.ForEachPair (meme.selected_branches, "_partition_", "_selection_",
    "_selection_ = utility.Filter (_selection_, '_value_', '_value_ == terms.tree_attributes.test');
     io.ReportProgressMessageMD('MEME',  'selector', 'Selected ' + Abs(_selection_) + ' branches to include in the MEME analysis: \\\`' + Join (', ',utility.Keys(_selection_)) + '\\\`')");


meme.pairwise_counts = genetic_code.ComputePairwiseDifferencesAndExpectedSites(meme.codon_data_info[terms.code], None);

selection.io.startTimer (meme.json [terms.json.timers], "Model fitting",1);

namespace meme {
    doGTR ("meme");
}


estimators.fixSubsetOfEstimates(meme.gtr_results, meme.gtr_results[terms.global]);

// Step 1 - Fit to MG94xREV
namespace meme {
    doPartitionedMG ("meme", FALSE);
}


io.ReportProgressMessageMD ("MEME", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");

meme.final_partitioned_mg_results = estimators.FitMGREV (meme.filter_names, meme.trees, meme.codon_data_info [terms.code], {
    terms.run_options.model_type: terms.local,
    terms.run_options.partitioned_omega: meme.selected_branches,
    terms.run_options.retain_lf_object: TRUE
}, meme.partitioned_mg_results);


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

meme.report.positive_site = {{"" + (1+((meme.filter_specification[meme.report.partition])[terms.data.coverage])[meme.report.site]),
                                    meme.report.partition + 1,
                                    Format(meme.report.row[0],7,3),
                                    Format(meme.report.row[3],7,3),
                                    Format(meme.report.row[4],7,3),
                                    Format(meme.report.row[5],7,3),
                                    "Yes, p = " + Format(meme.report.row[6],7,4),
                                    Format(meme.report.row[7],0,0)
}};

meme.site_results = {};

for (meme.partition_index = 0; meme.partition_index < meme.partition_count; meme.partition_index += 1) {
    meme.report.header_done = FALSE;
    meme.table_output_options[utility.getGlobalValue("terms.table_options.header")] = TRUE;

    meme.model_to_branch_bsrel = { "meme.bsrel" : utility.Filter (meme.selected_branches[meme.partition_index], '_value_', '_value_ == terms.tree_attributes.test'),
                             "meme.background_fel" : utility.Filter (meme.selected_branches[meme.partition_index], '_value_', '_value_ != terms.tree_attributes.test')};


    model.ApplyModelToTree( "meme.site_tree_fel", meme.trees[meme.partition_index], {terms.default : meme.site.background_fel}, None);
    model.ApplyModelToTree( "meme.site_tree_bsrel", meme.trees[meme.partition_index], None, meme.model_to_branch_bsrel);

    meme.site_patterns = alignments.Extract_site_patterns ((meme.filter_specification[meme.partition_index])[utility.getGlobalValue("terms.data.name")]);

    utility.ForEach (meme.site_tree_fel, "_node_",
            '_node_class_ = (meme.selected_branches[meme.partition_index])[_node_];
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
        ');




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
                                   terms.mpi.Variables : {{"meme.selected_branches","meme.branch_mixture","meme.pairwise_counts","meme.codon_data_info"}},
                                   terms.mpi.Functions : {{"meme.compute_branch_EBF"}}
                                 });


    /* run the main loop over all unique site pattern combinations */
    utility.ForEachPair (meme.site_patterns, "_pattern_", "_pattern_info_",
        '
            if (_pattern_info_[utility.getGlobalValue("terms.data.is_constant")]) {
                meme.store_results (-1,None,{"0" : "meme.site_likelihood",
                                             "1" : "meme.site_likelihood_bsrel",
                                             "2" : None,
                                             "3" : meme.partition_index,
                                             "4" : _pattern_info_,
                                             "5" : meme.site_model_mapping
                                     });
            } else {
                mpi.QueueJob (meme.queue, "meme.handle_a_site", {"0" : "meme.site_likelihood",
                                                                 "1" : "meme.site_likelihood_bsrel",
                                                                 "2" : alignments.serialize_site_filter
                                                                   ((meme.filter_specification[meme.partition_index])[utility.getGlobalValue("terms.data.name")],
                                                                   (_pattern_info_[utility.getGlobalValue("terms.data.sites")])[0]),
                                                                 "3" : meme.partition_index,
                                                                 "4" : _pattern_info_,
                                                                 "5" : meme.site_model_mapping
                                                                    },
                                                                    "meme.store_results");
            }
        '
    );

    mpi.QueueComplete (meme.queue);
    meme.partition_matrix = {Abs (meme.site_results[meme.partition_index]), Rows (meme.table_headers)};

    utility.ForEachPair (meme.site_results[meme.partition_index], "_key_", "_value_",
    '
        for (meme.index = 0; meme.index < Rows (meme.table_headers); meme.index += 1) {
            meme.partition_matrix [0+_key_][meme.index] = _value_[meme.index];
        }
    '
    );

    meme.site_results[meme.partition_index] = meme.partition_matrix;

}

meme.json [terms.json.MLE ] = {terms.json.headers   : meme.table_headers,
                               terms.json.content : meme.site_results };


io.ReportProgressMessageMD ("MEME", "results", "** Found _" + meme.report.count[0] + "_ sites under episodic diversifying positive selection at p <= " + meme.pvalue + "**");

selection.io.stopTimer (meme.json [terms.json.timers], "Total time");
selection.io.stopTimer (meme.json [terms.json.timers], "MEME analysis");

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
        `node_name`.`alpha_parameter` := (`alpha_factor`) * meme.branch_length__;
        `node_name`.`beta_parameter`  := (`beta_factor`)  * meme.branch_length__;
    ");
}
//----------------------------------------------------------------------------------------

function meme.apply_proportional_site_constraint.bsrel (tree_name, node_name, alpha_parameter, beta_parameter, beta_parameter2, mixture_parameter, alpha_factor, omega_factor, beta_factor, mixture_global, branch_length) {

    meme.branch_length = (branch_length[terms.parameters.synonymous_rate])[terms.fit.MLE];

    node_name = tree_name + "." + node_name;

    ExecuteCommands ("
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
    return {utility.getGlobalValue("terms.empirical_bayes_factor") : eBF__, utility.getGlobalValue("terms.posterior") : _posteriorProb__[1]};
}

//----------------------------------------------------------------------------------------
lfunction meme.handle_a_site (lf_fel, lf_bsrel, filter_data, partition_index, pattern_info, model_mapping) {

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

    fel = estimators.ExtractMLEs (lf_fel, model_mapping);
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
                  
    Optimize (results, ^lf_bsrel, {
            "OPTIMIZATION_METHOD" : "nedler-mead",
            "OPTIMIZATION_START_GRID" : 
             {
                "0" : {
                    "meme.site_beta_plus": ^"meme.site_beta_plus",
                    "meme.site_omega_minus": ^"meme.site_omega_minus",
                    "meme.site_mixture_weight": ^"meme.site_mixture_weight"
                },
                "1" : {
                    "meme.site_beta_plus": ^"meme.site_beta_plus" * 2,
                    "meme.site_omega_minus": 0.5,
                    "meme.site_mixture_weight": 0.5                
                },
                "2" : {
                    "meme.site_beta_plus": ^"meme.site_beta_plus" * 4,
                    "meme.site_omega_minus": 0.25,
                    "meme.site_mixture_weight": 0.25                
                },
                "3" : {
                    "meme.site_beta_plus": ^"meme.site_beta_plus",
                    "meme.site_omega_minus": 0.5,
                    "meme.site_mixture_weight": 0.5                
                },
                "4" : {
                    "meme.site_beta_plus": ^"meme.site_beta_plus",
                    "meme.site_omega_minus": 0.75,
                    "meme.site_mixture_weight": 0.8                
                },
                "5" : {
                    "meme.site_beta_plus": ^"meme.site_beta_plus" * 8,
                    "meme.site_omega_minus": 0.5,
                    "meme.site_mixture_weight": 0.8                
                },
                "6" : {
                    "meme.site_beta_plus": ^"meme.site_beta_plus",
                    "meme.site_omega_minus": 0,
                    "meme.site_mixture_weight": 0.01              
                }
            
             }
        });
        
    
    alternative = estimators.ExtractMLEs (lf_bsrel, model_mapping);
    alternative [utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];


    ancestral_info = ancestral.build (lf_bsrel,0,FALSE);

    //TODO
    branch_substitution_information = selection.substitution_mapper (ancestral_info ["MATRIX"],
                                                      ancestral_info ["TREE_AVL"],
                                                      ancestral_info ["AMBIGS"],
                                                      ^"meme.pairwise_counts",
                                                      ancestral_info ["MAPPING"],
                                                      (^"meme.codon_data_info")[utility.getGlobalValue("terms.code")]);


    DeleteObject (ancestral_info);

    branch_ebf       = {};
    branch_posterior = {};
    
 
    if (^"meme.site_beta_plus" > ^"meme.site_alpha" && ^"meme.site_mixture_weight" < 0.999999) {

        LFCompute (^lf_bsrel,LF_START_COMPUTE);
        LFCompute (^lf_bsrel,baseline);

        utility.ForEach (^bsrel_tree_id, "_node_name_",
        '
            if ((meme.selected_branches [^"`&partition_index`"])[_node_name_]  == utility.getGlobalValue("terms.tree_attributes.test")) {
                _node_name_res_ = meme.compute_branch_EBF (^"`&lf_bsrel`", ^"`&bsrel_tree_id`", _node_name_, ^"`&baseline`");
                (^"`&branch_ebf`")[_node_name_] = _node_name_res_[utility.getGlobalValue("terms.empirical_bayes_factor")];
                (^"`&branch_posterior`")[_node_name_] = _node_name_res_[utility.getGlobalValue("terms.posterior")];
            }
        '
        );

        LFCompute (^lf_bsrel,LF_DONE_COMPUTE);

        ^"meme.site_beta_plus" := ^"meme.site_alpha";
        Optimize (results, ^lf_bsrel);
        //io.SpoolLF (lf_bsrel, "/tmp/meme.debug", "MEME-null");

        Null = estimators.ExtractMLEs (lf_bsrel, model_mapping);
        Null [utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];



    } else {
        Null = alternative;
        utility.ForEach (^bsrel_tree_id, "_node_name_",
        '
            if ((meme.selected_branches [^"`&partition_index`"])[_node_name_]  == utility.getGlobalValue("terms.tree_attributes.test")) {
                (^"`&branch_ebf`")[_node_name_] = 1.0;
                (^"`&branch_posterior`")[_node_name_] = 0.0;
            }
        '
        );
    }


    return {"fel" : fel,
            utility.getGlobalValue("terms.alternative") : alternative,
            utility.getGlobalValue("terms.posterior") : branch_posterior,
            utility.getGlobalValue("terms.empirical_bayes_factor") : branch_ebf,
            utility.getGlobalValue("terms.branch_selection_attributes") : branch_substitution_information, //TODO: keep this attr?
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
                          0  // Log L | MEME 10
                          } };

      //console.log ( estimators.GetGlobalMLE (result["alternative"], ^"meme.parameter_site_mixture_weight"));

    if (None != result) { // not a constant site


        lrt = {utility.getGlobalValue("terms.LRT") : 2*((result[utility.getGlobalValue("terms.alternative")])[utility.getGlobalValue("terms.fit.log_likelihood")]-(result[utility.getGlobalValue("terms.Null")])[utility.getGlobalValue("terms.fit.log_likelihood")])};
        lrt [utility.getGlobalValue("terms.p_value")] = 2/3-2/3*(0.45*CChi2(lrt[utility.getGlobalValue("terms.LRT")],1)+0.55*CChi2(lrt[utility.getGlobalValue("terms.LRT")],2));

        result_row [0] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], utility.getGlobalValue("meme.parameter_site_alpha"));
        result_row [1] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], utility.getGlobalValue("meme.parameter_site_omega_minus")) * result_row[0];
        result_row [2] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], utility.getGlobalValue("meme.parameter_site_mixture_weight"));
        result_row [3] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], utility.getGlobalValue("meme.parameter_site_beta_plus"));
        result_row [4] = 1-result_row [2];
        result_row [5] = lrt [utility.getGlobalValue("terms.LRT")];
        result_row [6] = lrt [utility.getGlobalValue("terms.p_value")];
        result_row [9] = (result["fel"])[utility.getGlobalValue("terms.fit.log_likelihood")];
        result_row [10] = (result[utility.getGlobalValue("terms.alternative")])[utility.getGlobalValue("terms.fit.log_likelihood")];

        all_ebf = result[utility.getGlobalValue("terms.empirical_bayes_factor")];

        filtered_ebf = utility.Filter (utility.Filter (all_ebf, "_value_", "None!=_value_"), "_value_", "_value_>=100");

        if(None != filtered_ebf) {
            result_row [7] = utility.Array1D(filtered_ebf);
        } else {
            result_row [7] = 0;
        }

        sum = 0;
        alternative_lengths = ((result[utility.getGlobalValue("terms.alternative")])[utility.getGlobalValue("terms.branch_length")])[0];

        utility.ForEach (^"meme.site_tree_fel", "_node_",
                '_node_class_ = ((^"meme.selected_branches")[`&partition_index`])[_node_];
                 if (_node_class_ == utility.getGlobalValue("terms.tree_attributes.test")) {
                    `&sum` += ((`&alternative_lengths`)[_node_])[utility.getGlobalValue("terms.json.MLE")];
                 }
            ');

        result_row [8] = sum;
    } else {
        all_ebf = None;
    }

    utility.EnsureKey (^"meme.site_results", partition_index);

    sites_mapping_to_pattern = pattern_info[utility.getGlobalValue("terms.data.sites")];
    sites_mapping_to_pattern.count = utility.Array1D (sites_mapping_to_pattern);

    for (i = 0; i < sites_mapping_to_pattern.count; i+=1) {
        site_index = sites_mapping_to_pattern[i];
        ((^"meme.site_results")[partition_index])[site_index] = result_row;
        meme.report.echo (site_index, partition_index, result_row);
        if (None != all_ebf) {
            direct_index = 1+(((^"meme.filter_specification")[partition_index])[utility.getGlobalValue ("terms.data.coverage")])[site_index];
            selection.io.json_store_branch_attribute(^"meme.json", "EBF site " + direct_index + " (partition " + (1+partition_index) + ")", utility.getGlobalValue ("terms.json.branch_label"), site_index + utility.Array1D (^"meme.display_orders"),
                                                 partition_index,all_ebf);
        }
    }

    /*
    utility.ForEach (pattern_info[utility.getGlobalValue("terms.data.sites")], "_value_",
        '
            (meme.site_results[`&partition_index`])[_value_] = `&result_row`;
            meme.report.echo (_value_, `&partition_index`, `&result_row`);
        '
    );*/


    //assert (0);
}
