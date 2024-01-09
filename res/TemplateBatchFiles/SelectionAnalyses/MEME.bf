RequireVersion("2.5.54");

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
    Version 4 adds support for multiple hits and more than 2 rate classes on omega, as well as site-level imputation option
    ",
    terms.io.version: "4.0",
    terms.io.reference: "Detecting Individual Sites Subject to Episodic Diversifying Selection. _PLoS Genet_ 8(7): e1002764.",
    terms.io.authors: "Sergei L. Kosakovsky Pond, Steven Weaver",
    terms.io.contact: "spond@temple.edu",
    terms.io.requirements: "in-frame codon alignment and a phylogenetic tree",
    terms.settings: {}
};


io.DisplayAnalysisBanner(meme.analysis_description);


/*------------------------------------------------------------------------------
    Environment setup
*/

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);


/*------------------------------------------------------------------------------
    Globals
*/

/*
meme.parameter_site_alpha = "Site relative synonymous rate";
meme.parameter_site_omega_minus = "Omega ratio on (tested branches); negative selection or neutral evolution (&omega;- <= 1;)";
meme.parameter_site_beta_minus = "Site relative non-synonymous rate (tested branches); negative selection or neutral evolution (&beta;- <= &alpha;)";
meme.parameter_site_beta_plus = "Site relative non-synonymous rate (tested branches); unconstrained";
meme.parameter_site_mixture_weight = "Beta- category weight";
meme.parameter_site.delta = "Site 2H rate";
meme.parameter_site.psi = "Site 3H rate";

meme.parameter_site_beta_nuisance = "Site relative non-synonymous rate (untested branches)";
*/

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

terms.meme.bg_param_prefix =  "BG site parameter for ";
terms.meme.fg_param_prefix =  "FG site parameter for ";

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

namespace meme {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ("meme");
}

meme.pvalue  = io.PromptUser ("\n>Select the p-value threshold to use when testing for selection",meme.pvalue,0,1,FALSE);

KeywordArgument ("resample",  "[Advanced setting, will result in MUCH SLOWER run time] Perform parametric bootstrap resampling to derive site-level null LRT distributions up to this many replicates per site. Recommended use for small to medium (<30 sequences) datasets", "0");
meme.resample  = io.PromptUser ("\n>[Advanced setting, will result in MUCH SLOWER run time] Perform parametric bootstrap resampling to derive site-level null LRT distributions up to this many replicates per site. Recommended use for small to medium (<30 sequences) datasets",50,0,1000,TRUE);


KeywordArgument ("rates", "The number omega rate classes to include in the model [2-4]", meme.nrate_classes);
meme.nrate_classes = io.PromptUser ("The number omega rate classes to include in the model [2 is the strongly recommended default]", meme.nrate_classes, 2, 4, TRUE);
selection.io.json_store_setting  (meme.json, "rates", meme.nrate_classes);
meme.report.p_value_index = 2+2*meme.nrate_classes;

KeywordArgument ("multiple-hits",  "Include support for multiple nucleotide substitutions", "None");
meme.multi_hit = io.SelectAnOption ({
                                        {"Double", "Include branch-specific rates for double nucleotide substitutions"}
                                        {"Double+Triple", "Include branch-specific rates for double and triple nucleotide substitutions"}
                                        {"None", "[Default] Use standard models which permit only single nucleotide changes to occur instantly"}
                                  }, "Include support for multiple nucleotide substitutions");

selection.io.json_store_setting  (meme.json, "multihit", meme.multi_hit);

if (meme.multi_hit != "None") {
    KeywordArgument ("site-multihit", "Estimate multiple hit rates for each site", "Estimate");
    meme.multi_hit_option = io.SelectAnOption ({
                                        {"Estimate", "Include branch-specific rates for double nucleotide substitutions"}
                                        {"Global", "Use a plug-in estimate derived from the global model fit"}
                                  }, "Estimate multiple hit rates for each site");
                                  
    selection.io.json_store_setting  (meme.json, "site-multihit", meme.multi_hit);
}

KeywordArgument ("impute-states", "Use site-level model fits to impute likely character states for each sequence", "No");

meme.impute_states = io.SelectAnOption (
            {
                "Yes":"Impute marginal likelihoods for each codon at each sequence and each site",
                "No": "Do not impute marginal likelihoods for each codon at each sequence and each site",
            }, 
            "Impute likely states for sequences") == "Yes";


selection.io.json_store_setting (meme.json, terms.json.imputed_states, meme.impute_states );

KeywordArgument ("precision", "Optimization precision settings for preliminary fits", "standard");
meme.faster = io.SelectAnOption( {{"standard", "[Default] Use standard optimization settings."}, 
                                 {"reduced", "Cruder optimizations settings for faster fitting of preliminary models"}},
                                  "Optimization precision settings for preliminary fits") == "reduced";
selection.io.json_store_setting  (meme.json, "fit-mode", meme.faster);

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

meme.model_name = "models.codon.MG_REV.ModelDescription";
meme.define_constraints = None;
    
if (meme.multi_hit == "Double") {
    meme.model_name = "models.codon.MG_REV_MH.ModelDescription";
    meme.define_constraints = "meme.constrain_2H";
} else {
    if (meme.multi_hit == "Double+Triple") {
        meme.model_name = "models.codon.MG_REV_TRIP.ModelDescription";
        meme.define_constraints = "meme.constrain_3H";
    }
}

// Step 1 - Fit to MG94xREV
namespace meme {
    if (faster) {
        utility.ToggleEnvVariable("OPTIMIZATION_PRECISION", Abs(0.0001*gtr_results[^"terms.fit.log_likelihood"]));
    }
    doPartitionedMGModel ("meme", FALSE, model_name, define_constraints);
    if (faster) {
        utility.ToggleEnvVariable("OPTIMIZATION_PRECISION", None);
    }
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


KeywordArgument ("full-model", "Perform branch length re-optimization under the full codon model", "Yes");

meme.run_full_mg94 = "Yes" == io.SelectAnOption( {{"Yes", "[Default] Perform branch length re-optimization under the full codon model"}, 
                                         {"No", "Skip branch length re-optimization under the full codon model (faster but less precise)"}},
                                         "Perform branch length re-optimization under the full codon model");

meme.site_filter = selection.io.handle_subset_of_sites ();
selection.io.json_store_setting (meme.json, "site-filter", meme.site_filter );


io.ReportProgressMessageMD ("MEME", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");

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
        if (meme.faster) {
            utility.ToggleEnvVariable("OPTIMIZATION_PRECISION", Abs(0.0001*meme.partitioned_mg_results[^"terms.fit.log_likelihood"]));
        }
        meme.prefix = "meme";
        meme.final_partitioned_mg_results = estimators.FitCodonModel (meme.filter_names, meme.trees, meme.model_name, meme.codon_data_info [terms.code], {
            terms.run_options.model_type: terms.local,
            terms.run_options.partitioned_omega: meme.selected_branches,
            terms.run_options.retain_lf_object: TRUE,
            terms.run_options.apply_user_constraints: meme.define_constraints,
            terms.run_options.optimization_settings: {
                "OPTIMIZATION_METHOD" : "coordinate-wise"
            }
        }, meme.partitioned_mg_results);
        

        if (meme.faster) {
            utility.ToggleEnvVariable("OPTIMIZATION_PRECISION", None);
        }
        
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


io.ReportProgressMessageMD ("MEME", "codon-refit", "* Log(L) = " + Format(meme.final_partitioned_mg_results[terms.fit.log_likelihood],8,2));
io.ReportProgressMessageMD ("MEME", "codon-refit", "* " + selection.io.report_fit (meme.final_partitioned_mg_results, 0, (meme.codon_data_info)[utility.getGlobalValue ("terms.data.sample_size")]));
io.ReportProgressMessageMD ("MEME", "codon-refit", "* " +selection.io.report_fit_secondary_stats (meme.final_partitioned_mg_results));
meme.global_dnds = selection.io.extract_global_MLE_re (meme.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio);
meme.initial_dNdS = {};

for(i,k; in; meme.global_dnds) {
    if (None != regexp.Find (k[terms.description],utility.getGlobalValue("terms.tree_attributes.test"))) {
        meme.initial_dNdS['FG'] = k[terms.fit.MLE];
    } else {
        meme.initial_dNdS['BG'] = k[terms.fit.MLE];
    }
}


meme.multi_hit_MLES  = {};

if (meme.multi_hit != "None") {
    meme.global_delta = selection.io.extract_global_MLE_re (meme.final_partitioned_mg_results, "^" + terms.parameters.multiple_hit_rate);
    for (i,k; in; meme.global_delta) {
        meme.global_dnds + k;
        meme.multi_hit_MLES [terms.parameters.multiple_hit_rate] = k [terms.fit.MLE];
    }
    meme.multi_hit_MLES [terms.parameters.triple_hit_rate] = 0;
    if (meme.multi_hit != "Double")  {
        meme.global_delta = selection.io.extract_global_MLE_re (meme.final_partitioned_mg_results, "^" + terms.parameters.triple_hit_rate);
       
        for (i,k; in; meme.global_delta) {
            meme.global_dnds + k;
            meme.multi_hit_MLES [terms.parameters.triple_hit_rate] = k [terms.fit.MLE];
            meme.multi_hit_MLES [terms.parameters.triple_hit_rate_syn] = k [terms.fit.MLE];
        }    
    }
}


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


meme.model_templates = {'BG' : meme.model_name};


if (meme.multi_hit == "Double") {
    meme.model_templates['FG'] = "models.codon.BS_REL_Per_Branch_Mixing.MH.ModelDescription";
} else {
    if (meme.multi_hit == "Double+Triple") {
        meme.model_templates['FG'] = "models.codon.BS_REL_Per_Branch_Mixing.MH.ModelDescription";
    } else {
        meme.model_templates['FG'] = "models.codon.BS_REL_Per_Branch_Mixing.ModelDescription";
    }
}





meme.site.background_fel = model.generic.DefineModel(meme.model_templates["BG"],
        "meme.background_fel", {
            "0": parameters.Quote(terms.local),
            "1": meme.codon_data_info[terms.code]
        },
        meme.filter_names,
        None);



meme.site.bsrel =  model.generic.DefineMixtureModel(meme.model_templates['FG'],
        "meme.bsrel", {
            "0": parameters.Quote(terms.local),
            "1": meme.codon_data_info[terms.code],
            "2": parameters.Quote (meme.nrate_classes) // the number of rate classes
        },
        meme.filter_names,
        None);



/***
    this defines, for each branch in the tree, site-level global scaler parameters,
    and how they correspond to individual local parameters in the tree
    
    'BG' => background / FEL
    'FG' => foreground  / MEME
        'term' => {
            'local' : local rate,
            'scaler' : global scaler,
            'LB' : lower bound or None
            'UB' : upper bound or None
            'constraint' : string or None // which constraint to apply to the local parameter across the tree
        }
        
    
        
**/


lfunction meme.check_parameter_and_store (model, term, target_rate, lb, ub, constraint, target_dict, soft, init_value) {
    local_parameter = model.generic.GetLocalParameter (model, term);
    if (soft && None == local_parameter) {
        return;
    }
    io.CheckAssertion ("None!=`&local_parameter`", "Could not find expected local parameters for `term`");
    target_dict [term] = {
        'local' : local_parameter__,
        'scaler' : target_rate__,
        'LB' : lb__,
        'UB' : ub__,
        'constraint' : constraint__,
        'init' : init_value__
    };
}

meme.scaler_mapping = {
    'FG' : {},
    'FEL-FG' : {},
    'BG' : {}
};

if (meme.multi_hit != "None") {
    
    for (m;in;{
        0 : {{"meme.site.background_fel","BG"}},
        1 : {{"meme.site.bsrel","FG"}},
        2 : {{"meme.site.background_fel","FEL-FG"}}
     } ) {
        
        if (meme.multi_hit != "Double") {
            meme.3H = "meme.site_psi";
        } else {
            meme.3H = 0;
        }
        
        meme.temp = {
            terms.parameters.multiple_hit_rate   : "meme.site_delta", 
            terms.parameters.triple_hit_rate     : meme.3H,
            terms.parameters.triple_hit_rate_syn : meme.3H,
             
        };

        for (k,v; in; meme.temp) {
            if (meme.multi_hit_option == "Global") {
                v = meme.multi_hit_MLES[k];
            } 
            meme.check_parameter_and_store (^(m[0]), k , v, None, None, "'" + v + "'", meme.scaler_mapping[m[1]], TRUE, meme.multi_hit_MLES[k]);
        }
    }    
}


meme.check_parameter_and_store (meme.site.background_fel, terms.parameters.synonymous_rate, "meme.site_alpha", None, 1e4, "'meme.site_alpha*(' + branch_length + ')'", meme.scaler_mapping['BG'], FALSE, 1);
meme.check_parameter_and_store (meme.site.background_fel, terms.parameters.synonymous_rate, "meme.site_alpha", None, 1e4, "'meme.site_alpha*(' + branch_length + ')'", meme.scaler_mapping['FEL-FG'], FALSE, 1);

meme.check_parameter_and_store (meme.site.bsrel,          terms.parameters.synonymous_rate, "meme.site_alpha", None, 1e4, "'meme.site_alpha*(' + branch_length + ')'", meme.scaler_mapping['FG'], FALSE, None);
meme.check_parameter_and_store (meme.site.background_fel, terms.parameters.nonsynonymous_rate, "meme.site_beta_nuisance", None, 1e5, "'meme.site_beta_nuisance*(' + branch_length + ')'", meme.scaler_mapping['BG'], FALSE, meme.initial_dNdS['BG']);

meme.omega_distribution = {
    terms.parameters.rates : {},
    terms.parameters.weights : {}
};

for (i = 1; i <= meme.nrate_classes - 1; i+=1) {
    meme.check_parameter_and_store (meme.site.bsrel,  terms.AddCategory (terms.parameters.nonsynonymous_rate,i) , "meme.site_omega_" + i, 0, 1, "'meme.site_alpha*meme.site_omega_`i`*(' + branch_length + ')'", meme.scaler_mapping['FG'], FALSE, None);
    meme.check_parameter_and_store (meme.site.bsrel,  terms.AddCategory (terms.mixture.mixture_aux_weight,i) , "meme.weight_" + i, 1e-8, 1, "'meme.weight_`i`'", meme.scaler_mapping['FG'], FALSE, None);
    
    meme.omega_distribution [terms.parameters.rates ] + ("meme.site_omega_" + i);
    meme.omega_distribution [terms.parameters.weights ] + ("meme.weight_" + i);
    
}

meme.check_parameter_and_store (meme.site.bsrel,  terms.AddCategory (terms.parameters.nonsynonymous_rate,meme.nrate_classes) , "meme.site_beta_plus", None, 1e5, "'meme.site_beta_plus*(' + branch_length + ')'", meme.scaler_mapping['FG'], FALSE, None);
meme.check_parameter_and_store (meme.site.background_fel,  terms.parameters.nonsynonymous_rate , "meme.site_beta_plus", None, 1e5, "'meme.site_beta_plus*(' + branch_length + ')'", meme.scaler_mapping['FEL-FG'], FALSE, meme.initial_dNdS['FG']);

meme.omega_distribution [terms.parameters.rates ] + "meme.site_beta_plus";

meme.scaler_mapping['OMEGA'] = meme.omega_distribution;
meme.scaler_mapping['OMEGA_DIST'] = parameters.helper.stick_breaking  (meme.omega_distribution[terms.parameters.weights ], None);

meme.global_cache = {};



for (i,v; in; meme.scaler_mapping['BG']) {
    if (Type (v["scaler"]) == "String") {
        model.generic.AddGlobal (meme.site.background_fel, v["scaler"], terms.meme.bg_param_prefix  + i);
        parameters.DeclareGlobalWithRanges (v["scaler"], None, v["LB"], v["UB"]);
    }
}

for (i,v; in; meme.scaler_mapping['FG']) {
    if (Type (v["scaler"]) == "String") {
        model.generic.AddGlobal (meme.site.background_fel, v["scaler"], terms.meme.fg_param_prefix  + i);
        parameters.DeclareGlobalWithRanges (v["scaler"], None, v["LB"], v["UB"]);
    }
}



meme.site_model_mapping = {
                           "meme.background_fel" : meme.site.background_fel,
                           "meme.bsrel" : meme.site.bsrel,
                          };


selection.io.startTimer (meme.json [terms.json.timers], "MEME analysis", 2);




meme.report.count = {{0}};



meme.table_headers = {9 + meme.nrate_classes*2  + (meme.multi_hit == "Double") + 2*(meme.multi_hit == "Double+Triple"),2};

meme.table_headers[0][0] = "&alpha;"; meme.table_headers[0][1] = "Synonymous substitution rate at a site";

for (i=0; i<meme.nrate_classes-1; i+=1) {
    meme.table_headers[1+2*i][0] = "&beta;<sup>" + (i+1) + "</sup>"; meme.table_headers[1+2*i][1] = "Non-synonymous substitution rate at a site for the negative/neutral evolution component " + (i+1);
    meme.table_headers[2+2*i][0] = "p<sup>" + (i+1) + "</sup>"; meme.table_headers[2+2*i][1] = "Mixture distribution weight allocated to negative/neutral evolution component " + (i+1);
}

meme.table_headers[1+2*i][0] = "&beta;<sup>+</sup>"; meme.table_headers[1+2*i][1] = "Non-synonymous substitution rate at a site for the positive selection component";
meme.table_headers[2+2*i][0] = "p<sup>+</sup>"; meme.table_headers[2+2*i][1] = "Mixture distribution weight allocated to the positive selection component";

meme.row_index = 2*meme.nrate_classes+1;
 
meme.table_headers[meme.row_index][0] = "LRT"; meme.table_headers[meme.row_index][1] = "Likelihood ratio test statistic for episodic diversification, i.e., p<sup>+</sup> &gt; 0 <emph>and<emph> &beta;<sup>+</sup> &gt; &alpha;"; meme.row_index+=1;

meme.table_headers[meme.row_index][0] = "p-value"; meme.table_headers[meme.row_index][1] = "Asymptotic p-value for episodic diversification, i.e., p<sup>+</sup> &gt; 0 <emph>and<emph> &beta;<sup>+</sup> &gt; &alpha;"; meme.row_index+=1;

meme.table_headers[meme.row_index][0] = "# branches under selection"; meme.table_headers[meme.row_index][1] = "The (very approximate and rough) estimate of how many branches may have been under selection at this site, i.e., had an empirical Bayes factor of 100 or more for the &beta;<sup>+</sup> rate"; meme.row_index+=1;

meme.table_headers[meme.row_index][0] = "Total branch length"; meme.table_headers[meme.row_index][1] = "The total length of branches contributing to inference at this site, and used to scale dN-dS"; meme.row_index+=1;

meme.table_headers[meme.row_index][0] = "MEME LogL"; meme.table_headers[meme.row_index][1] = "Site Log-likelihood under the MEME model"; meme.row_index+=1;
meme.table_headers[meme.row_index][0] = "FEL LogL"; meme.table_headers[meme.row_index][1] = "Site Log-likelihood under the FEL model"; meme.row_index+=1;
 
meme.table_headers[meme.row_index][0] = "FEL &alpha;"; meme.table_headers[meme.row_index][1] = "Synonymous substitution rate at a site under the FEL model";  meme.row_index+=1;
meme.table_headers[meme.row_index][0] = "FEL &beta;"; meme.table_headers[meme.row_index][1] = "Non-synonymous substitution rate at a site under the FEL model";  meme.row_index+=1;
                       

if (meme.multi_hit == "Double") {
    meme.table_headers[meme.row_index][0] = "&delta;"; meme.table_headers[meme.row_index][0] = "Relative rate estimate for 2-nucleotide substitutions";  meme.row_index+=1;
}

if (meme.multi_hit == "Double+Triple") {
    meme.table_headers[meme.row_index][0] = "&psi;"; meme.table_headers[meme.row_index][0] = "Relative rate estimate for 3-nucleotide substitutions";  meme.row_index+=1;
}


/**
This table is meant for HTML rendering in the results web-app; can use HTML characters, the second column
is 'pop-over' explanation of terms. This is ONLY saved to the JSON file. For Markdown screen output see
the next set of variables.
*/

meme.table_screen_output  = {0 : "Codon", 
                             1 : "Partition", 
                             2:  "alpha", 
                             3:  "non-syn rate (beta) distribution, rates : weights", 
                             4:  "LRT", 
                             5:  "Episodic selection detected?", 
                             6:  "# branches", 
                             7:  "         List of most common codon substitutions at this site          "
                             };
                             
meme.table_output_options = {terms.table_options.header : TRUE, terms.table_options.minimum_column_width: 12, terms.table_options.align : "center"};


meme.site.composition.string = "";
meme.site.distribution_string = "";
meme.site_multihit_string = "";

if (meme.multi_hit != "None") {
    if (meme.multi_hit == "Double") {
        meme.table_screen_output + "Double hit rate";
    } else {
        meme.table_screen_output + "Double/Triple hit rates";
    }
    meme.report.positive_site = {{"" + (1+((meme.filter_specification[meme.report.partition])[terms.data.coverage])[meme.report.site]),
                                    meme.report.partition + 1,
                                    Format(meme.report.row[0],7,3),
                                    meme.site.distribution_string,
                                    Format(meme.report.row[meme.report.p_value_index-1],7,3),
                                    "Yes, p = " + Format(meme.report.row[meme.report.p_value_index],7,4),
                                    Format(meme.report.row[meme.report.p_value_index+1],0,0),
                                    meme.site.composition.string,
                                    meme.site_multihit_string
    }};
} else {
    meme.report.positive_site = {{"" + (1+((meme.filter_specification[meme.report.partition])[terms.data.coverage])[meme.report.site]),
                                    meme.report.partition + 1,
                                    Format(meme.report.row[0],7,3),
                                    meme.site.distribution_string,
                                    Format(meme.report.row[meme.report.p_value_index-1],7,3),
                                    "Yes, p = " + Format(meme.report.row[meme.report.p_value_index],7,4),
                                    Format(meme.report.row[meme.report.p_value_index+1],0,0),
                                    meme.site.composition.string
    }};

}

meme.site_results = {};
meme.imputed_leaf_states = {};
meme.site_LRT = {};

meme.output_file_path = meme.codon_data_info[terms.json.json];

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
    
    utility.ToggleEnvVariable ("LOOP_TREE_ITERATOR_PREORDER", TRUE);
    
    is_root = TRUE;
    
    for (_node_; in; meme.site_tree_fel) {
        if (is_root) {
            is_root = FALSE;
            continue;
        }
        _bl_ = ((meme.final_partitioned_mg_results[terms.branch_length])[meme.partition_index])[_node_];
        _node_class_ = (meme.selected_branches[meme.partition_index])[_node_];
        if (_node_class_ != terms.tree_attributes.test) {
            meme.apply_site_constraints ("meme.site_tree_bsrel",_node_,_bl_,meme.scaler_mapping ['BG']);
            meme.apply_site_constraints ("meme.site_tree_fel",_node_,_bl_,meme.scaler_mapping ['BG']);
        } else {
            meme.apply_site_constraints ("meme.site_tree_bsrel",_node_,_bl_,meme.scaler_mapping ['FG']);     
            meme.apply_site_constraints ("meme.site_tree_fel",_node_,_bl_,meme.scaler_mapping ['FEL-FG']);       
        }
     }
    
     utility.ToggleEnvVariable ("LOOP_TREE_ITERATOR_PREORDER", None);


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
                                   terms.mpi.Variables : {{"meme.selected_branches","meme.impute_states","meme.pairwise_counts","meme.codon_data_info","meme.resample","meme.site_filter","meme.output_file_path"}},
                                   terms.mpi.Functions : {{"meme.compute_branch_EBF","selection.io.sitelist_matches_pattern"}}
                                 });


    meme.pattern_count_all  = utility.Array1D (meme.site_patterns);
    meme.pattern_count_this = 0;

    /* run the main loop over all unique site pattern combinations */
    utility.ForEachPair (meme.site_patterns, "_pattern_", "_pattern_info_",
        '
            meme.pattern_count_this += 1;
            io.ReportProgressBar("", "Working on site pattern " + (meme.pattern_count_this) + "/" +  meme.pattern_count_all + " in partition " + (1+meme.partition_index));
            meme.run_site = selection.io.sitelist_matches_pattern (_pattern_info_[terms.data.sites], meme.site_filter["site-filter"], FALSE);

            if (_pattern_info_[utility.getGlobalValue("terms.data.is_constant")] || (!meme.run_site)) {
                meme.store_results (-1,None,{"0" : "meme.site_likelihood",
                                             "1" : "meme.site_likelihood_bsrel",
                                             "2" : None,
                                             "3" : meme.partition_index,
                                             "4" : _pattern_info_,
                                             "5" : meme.site_model_mapping,
                                             "6" : meme.scaler_mapping,
                                             "7" : 0
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
                                                                 "6" : meme.scaler_mapping,
                                                                 "7" : 0
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

if (meme.impute_states) {
    (meme.json [terms.json.MLE])[terms.json.imputed_states] = meme.imputed_leaf_states;
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
function meme.apply_site_constraints (tree_name, node_name, bl, parameters) {
    node_name = tree_name + "." + node_name;
    branch_length = bl[terms.fit.MLE];
    for (p,i; in; parameters) {
        if (Type (i["constraint"]) == "String") {
            ExecuteCommands ("`node_name`.`i['local']` := " + Eval (i["constraint"]));
        }
    }   
}


//----------------------------------------------------------------------------------------
function meme.apply_2H_constraint (tree_name, node_name, delta_parameter, delta_factor) {
    node_name = tree_name + "." + node_name;
    ExecuteCommands ("
        `node_name`.`delta_parameter` := (`delta_factor`);
    ");
}

//----------------------------------------------------------------------------------------
function meme.apply_3H_constraint (tree_name, node_name, psi_parameter1, psi_parameter2, psi_factor) {
    node_name = tree_name + "." + node_name;
    ExecuteCommands ("
        `node_name`.`psi_parameter1` := (`psi_factor`);
        `node_name`.`psi_parameter2` := (`psi_factor`);
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
lfunction meme.compute_branch_EBF (lf_id, tree_name, branch_name, baseline, scaler_mapping) {
// TODO: figure out why LFCompute fails if this is run as an `lfunction`

    // compute conditional likelihoods for 
    // omega = omega_i
    
    conditionals = {};
    weights      = {};
    locals       = {};
    scaler       = {};
    nrates       = utility.Array1D (scaler_mapping["OMEGA_DIST"]);
    
    for (is,w;in; scaler_mapping["OMEGA_DIST"]) {
        i = +is;
        local = ((scaler_mapping['FG'])[terms.AddCategory (^"terms.mixture.mixture_aux_weight",i+1)])["local"];
        scaler + ((scaler_mapping['FG'])[terms.AddCategory (^"terms.mixture.mixture_aux_weight",i+1)])["scaler"];
        if (None  != local) {
            locals + local;
            weights + utility.EvalInGlobalNamespace (w);
        }
    }
    
    positive_weight = 1 - (+weights);
    
    
    for (i = 0; i <= nrates - 2; i+=1) {
        ^("`tree_name`.`branch_name`." + locals[i]) = 1;
        for (j = 0; j <=  nrates-2; j+=1) {
            if ( i != j) {
                ^("`tree_name`.`branch_name`." + locals[j]) = 0;     
            }   
        }
        LFCompute (^lf_id,LOGL0);
        conditionals[i] = LOGL0;
    }
    


    for (i,l; in; locals) {
        utility.ExecuteInGlobalNamespace ("`tree_name`.`branch_name`." + l + " := " + scaler[i]);  
    }

    _posteriorProb = math.ComputePosteriorsFromLF (conditionals, weights, baseline);
    
    if (positive_weight != 1 && positive_weight!= 0) {
        _priorOdds = positive_weight / (1-positive_weight);
    } else {
        _priorOdds = 0;
    }
    
    //console.log (conditionals);
    //console.log (baseline);
    //console.log (_posteriorProb);

    if ( _priorOdds != 0) {
        eBF = _posteriorProb[nrates-1] / (1 - _posteriorProb[nrates-1]) / _priorOdds;
    } else {
        eBF = 1;
    }
    
    return {utility.getGlobalValue("terms.empirical_bayes_factor") : eBF__, 
            utility.getGlobalValue("terms.posterior") : _posteriorProb__};
}
//----------------------------------------------------------------------------------------
lfunction meme.apply_initial_estimates (spec) {
    for (t,d;in;spec) {
        if (None != d["init"]) {
            if (Type (d["scaler"]) == "String") {
                ^(d["scaler"]) = d["init"];
            }
        }
    }
}


//----------------------------------------------------------------------------------------
lfunction meme.set_ith_omega (spec,i,rate,weight) {
    ^(((spec["OMEGA"])[^"terms.parameters.rates"])[i]) = rate;
    if (None != weight) {
        ^(((spec["OMEGA"])[^"terms.parameters.weights"])[i]) = weight;
    }
}

//----------------------------------------------------------------------------------------
lfunction meme.handle_a_site (lf_fel, lf_bsrel, filter_data, partition_index, pattern_info, model_mapping, scaler_mapping, sim_mode) {

    //console.log (pattern_info);
    //#profile START;
    GetString   (lfInfo, ^lf_fel,-1);   

    //utility.SetEnvVariable ("VERBOSITY_LEVEL", 100);

    //TODO Datafilters hardcode, Trees hardcode.
    ExecuteCommands (filter_data);
    __make_filter ((lfInfo["Datafilters"])[0]);

    GetString (lfInfo, ^lf_bsrel,-1);
    __make_filter ((lfInfo["Datafilters"])[0]);
    
    utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);
    utility.SetEnvVariable ("ASSUME_REVERSIBLE_MODELS", TRUE);

    rate.fel.alpha   = (((scaler_mapping['BG']))[^"terms.parameters.synonymous_rate"])["scaler"];
    rate.fel.beta_fg = (((scaler_mapping['FEL-FG']))[^"terms.parameters.nonsynonymous_rate"])["scaler"];
    rate.fel.beta_bg = (((scaler_mapping['BG']))[^"terms.parameters.nonsynonymous_rate"])["scaler"];
        
    meme.apply_initial_estimates (scaler_mapping["BG"]);
    meme.apply_initial_estimates (scaler_mapping["FEL-FG"]);
    
    //utility.SetEnvVariable ("VERBOSITY_LEVEL", 110);
    
    Optimize (results, ^lf_fel
        , {"OPTIMIZATION_METHOD" : "nedler-mead", OPTIMIZATION_PRECISION: 1e-4}
    );



    fel = estimators.ExtractMLEsOptions (lf_fel, model_mapping, {^"terms.globals_only" : TRUE});
    fel[utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];
 
    rate_count = utility.Array1D ((scaler_mapping["OMEGA"])[^"terms.parameters.rates"]);

    bsrel_tree_id = (lfInfo["Trees"])[0];
    
    character_map = None;    

    if (rate_count == 2) {
        // separate clauses for initial parameter values
        if ( ^rate.fel.alpha <  ^rate.fel.beta_fg) { // FEL estimated beta > alpha
            meme.set_ith_omega (scaler_mapping, 0, 0.25, 0.25);
        } else { // FEL estimated beta <= alpha
            if (^rate.fel.alpha  > 0) {
                omega_rate = ^rate.fel.beta_fg / ^rate.fel.alpha;
            } else {
                omega_rate = 1;
            }
            meme.set_ith_omega (scaler_mapping, 0, omega_rate, 0.75);
            meme.set_ith_omega (scaler_mapping, 1, 1.5*^rate.fel.alpha , None);
         }
        
         rate.meme.mixture_weight1 = ((scaler_mapping["OMEGA"])[^"terms.parameters.weights"])[0];
         rate.meme.omega1          = ((scaler_mapping["OMEGA"])[^"terms.parameters.rates"])[0];

         meme.apply_initial_estimates (scaler_mapping["FG"]);
         
         initial_guess_grid = {
                "0" : {
        
                    rate.fel.alpha : ^rate.fel.alpha,
                    rate.fel.beta_bg : ^rate.fel.beta_bg,
                    rate.fel.beta_fg: ^rate.fel.beta_fg,
                    rate.meme.omega1: ^rate.meme.omega1,
                    rate.meme.mixture_weight1: ^rate.meme.mixture_weight1
           
                },
                "1" : {
                    rate.fel.alpha : ^rate.fel.alpha,
                    rate.fel.beta_bg : ^rate.fel.beta_bg,
                    rate.fel.beta_fg: ^rate.fel.beta_fg * 2,
                    rate.meme.omega1: 0.5,
                    rate.meme.mixture_weight1: 0.5                
                },
                "2" : {
                    rate.fel.alpha : ^rate.fel.alpha,
                    rate.fel.beta_bg : ^rate.fel.beta_bg,
                    rate.fel.beta_fg: ^rate.fel.beta_fg * 4,
                    rate.meme.omega1: 0.25,
                    rate.meme.mixture_weight1: 0.25                
                },
                "3" : {
                    rate.fel.alpha : ^rate.fel.alpha,
                    rate.fel.beta_bg : ^rate.fel.beta_bg,
                    rate.fel.beta_fg: ^rate.fel.beta_fg,
                    rate.meme.omega1: 0.5,
                    rate.meme.mixture_weight1: 0.5                
                },
                "4" : {
                    rate.fel.alpha : ^rate.fel.alpha,
                    rate.fel.beta_bg : ^rate.fel.beta_bg,
                    rate.fel.beta_fg: ^rate.fel.beta_fg,
                    rate.meme.omega1: 0.75,
                    rate.meme.mixture_weight1: 0.8                
                },
                "5" : {
                    rate.fel.alpha : ^rate.fel.alpha,
                    rate.fel.beta_bg : ^rate.fel.beta_bg,
                    rate.fel.beta_fg: ^rate.fel.beta_fg * 8,
                    rate.meme.omega1: 0.5,
                    rate.meme.mixture_weight1: 0.8                
                },
                "6" : {
                    rate.fel.alpha : ^rate.fel.alpha,
                    rate.fel.beta_bg : ^rate.fel.beta_bg,
                    rate.fel.beta_fg: ^rate.fel.beta_fg,
                    rate.meme.omega1: 0,
                    rate.meme.mixture_weight1: 0.01              
                },
                "7" : {
                    rate.fel.alpha : ^rate.fel.alpha,
                    rate.fel.beta_bg : ^rate.fel.beta_bg,
                    rate.fel.beta_fg: ^rate.fel.beta_fg,
                    rate.meme.omega1: 0,
                    rate.meme.mixture_weight1: 0.7         
                }
    
        };
    } else {
        if (rate_count == 3) {
            // separate clauses for initial parameter values
            if ( ^rate.fel.alpha <  ^rate.fel.beta_fg) { // FEL estimated beta > alpha
                meme.set_ith_omega (scaler_mapping, 0, 0.25, 5);
                meme.set_ith_omega (scaler_mapping, 1, 0.5, 0.5);
                meme.set_ith_omega (scaler_mapping, 2, 0.5, None);
            } else { // FEL estimated beta <= alpha
                if (^rate.fel.alpha  > 0) {
                    omega_rate = ^rate.fel.beta_fg / ^rate.fel.alpha;
                } else {
                    omega_rate = 1;
                }
                meme.set_ith_omega (scaler_mapping, 0, omega_rate, 0.5);
                meme.set_ith_omega (scaler_mapping, 1, Min (1, Max (0.5,omega_rate*2)), 0.5);
                meme.set_ith_omega (scaler_mapping, 2, 1.5*^rate.fel.alpha , None);
             }
        
             rate.meme.mixture_weight1 = ((scaler_mapping["OMEGA"])[^"terms.parameters.weights"])[0];
             rate.meme.omega1          = ((scaler_mapping["OMEGA"])[^"terms.parameters.rates"])[0];
             rate.meme.mixture_weight2 = ((scaler_mapping["OMEGA"])[^"terms.parameters.weights"])[1];
             rate.meme.omega2          = ((scaler_mapping["OMEGA"])[^"terms.parameters.rates"])[1];
         
      
             initial_guess_grid = {
                    "0" : {
                        rate.fel.alpha : ^rate.fel.alpha,
                        rate.fel.beta_bg : ^rate.fel.beta_bg,
                        rate.fel.beta_fg: ^rate.fel.beta_fg,
                        rate.meme.omega1: ^rate.meme.omega1,
                        rate.meme.mixture_weight1: ^rate.meme.mixture_weight1,
                        rate.meme.omega2: ^rate.meme.omega2,
                        rate.meme.mixture_weight2: ^rate.meme.mixture_weight2
                    },
                    "1" : {
                        rate.fel.alpha : ^rate.fel.alpha,
                        rate.fel.beta_bg : ^rate.fel.beta_bg,
                        rate.fel.beta_fg: ^rate.fel.beta_fg * 2,
                        rate.meme.omega1: 0.0,
                        rate.meme.mixture_weight1: 0.5,                
                        rate.meme.omega2: 0.5,
                        rate.meme.mixture_weight2: 0.5
                    },
                    "2" : {
                        rate.fel.alpha : ^rate.fel.alpha,
                        rate.fel.beta_bg : ^rate.fel.beta_bg,
                        rate.fel.beta_fg: ^rate.fel.beta_fg * 4,
                        rate.meme.omega1: 0.25,
                        rate.meme.mixture_weight1: 0.25,                
                        rate.meme.omega2: 0.5,
                        rate.meme.mixture_weight2: 0.25
                    },
                    "3" : {
                        rate.fel.alpha : ^rate.fel.alpha,
                        rate.fel.beta_bg : ^rate.fel.beta_bg,
                        rate.fel.beta_fg: ^rate.fel.beta_fg,
                        rate.meme.omega1: 0.5,
                        rate.meme.mixture_weight1: 0.5,                
                        rate.meme.omega2: 1.0,
                        rate.meme.mixture_weight2: 0.5
                    },
                    "4" : {
                        rate.fel.alpha : ^rate.fel.alpha,
                        rate.fel.beta_bg : ^rate.fel.beta_bg,
                        rate.fel.beta_fg: ^rate.fel.beta_fg,
                        rate.meme.omega1: 0.25,
                        rate.meme.mixture_weight1: 0.6,               
                        rate.meme.omega2: 0.75,
                        rate.meme.mixture_weight2: 0.25
                   },
                    "5" : {
                        rate.fel.alpha : ^rate.fel.alpha,
                        rate.fel.beta_bg : ^rate.fel.beta_bg,
                        rate.fel.beta_fg: ^rate.fel.beta_fg * 8,
                        rate.meme.omega1: 0.25,
                        rate.meme.mixture_weight1: 0.6,             
                        rate.meme.omega2: 0.75,
                        rate.meme.mixture_weight2: 0.25
                    },
                    "6" : {
                        rate.fel.alpha : ^rate.fel.alpha,
                        rate.fel.beta_bg : ^rate.fel.beta_bg,
                        rate.fel.beta_fg: ^rate.fel.beta_fg,
                        rate.meme.omega1: 0,
                        rate.meme.mixture_weight1: 0.01,              
                        rate.meme.omega2: 1,
                        rate.meme.mixture_weight2: 0.10
                    },
                    "7" : {
                        rate.fel.alpha : ^rate.fel.alpha,
                        rate.fel.beta_bg : ^rate.fel.beta_bg,
                        rate.fel.beta_fg: ^rate.fel.beta_fg,
                        rate.meme.omega1: 0,
                        rate.meme.mixture_weight1: 0.7,              
                        rate.meme.omega2: 0.5,
                        rate.meme.mixture_weight2: 0.7
                    }
    
            };
        } else {
            // separate clauses for initial parameter values
            if ( ^rate.fel.alpha <  ^rate.fel.beta_fg) { // FEL estimated beta > alpha
                meme.set_ith_omega (scaler_mapping, 0, 0.0, 0.5);
                meme.set_ith_omega (scaler_mapping, 1, 0.25, 0.5);
                meme.set_ith_omega (scaler_mapping, 2, 0.5, 0.5);
                meme.set_ith_omega (scaler_mapping, 3, ^rate.fel.beta_fg, None);
            } else { // FEL estimated beta <= alpha
                if (^rate.fel.alpha  > 0) {
                    omega_rate = ^rate.fel.beta_fg / ^rate.fel.alpha;
                } else {
                    omega_rate = 1;
                }
                meme.set_ith_omega (scaler_mapping, 0, omega_rate, 0.5);
                meme.set_ith_omega (scaler_mapping, 1, Min (1, omega_rate*0.25), 0.5);
                meme.set_ith_omega (scaler_mapping, 2, Min (1, 0.5*omega_rate) , 0.5);
                meme.set_ith_omega (scaler_mapping, 3, 1.5*^rate.fel.alpha , None);
             }
        
             rate.meme.mixture_weight1 = ((scaler_mapping["OMEGA"])[^"terms.parameters.weights"])[0];
             rate.meme.omega1          = ((scaler_mapping["OMEGA"])[^"terms.parameters.rates"])[0];
             rate.meme.mixture_weight2 = ((scaler_mapping["OMEGA"])[^"terms.parameters.weights"])[1];
             rate.meme.omega2          = ((scaler_mapping["OMEGA"])[^"terms.parameters.rates"])[1];
             rate.meme.mixture_weight3 = ((scaler_mapping["OMEGA"])[^"terms.parameters.weights"])[2];
             rate.meme.omega3          = ((scaler_mapping["OMEGA"])[^"terms.parameters.rates"])[2];
             
             initial_guess_grid = {
                    "0" : {
                        rate.fel.alpha : ^rate.fel.alpha,
                        rate.fel.beta_bg : ^rate.fel.beta_bg,
                        rate.fel.beta_fg: ^rate.fel.beta_fg,
                        rate.meme.omega1: ^rate.meme.omega1,
                        rate.meme.mixture_weight1: ^rate.meme.mixture_weight1,
                        rate.meme.omega2: ^rate.meme.omega2,
                        rate.meme.mixture_weight2: ^rate.meme.mixture_weight2,
                        rate.meme.omega3: ^rate.meme.omega3,
                        rate.meme.mixture_weight3: ^rate.meme.mixture_weight3
                    },
                    "1" : {
                        rate.fel.alpha : ^rate.fel.alpha,
                        rate.fel.beta_bg : ^rate.fel.beta_bg,
                        rate.fel.beta_fg: ^rate.fel.beta_fg * 2,
                        rate.meme.omega1: 0.0,
                        rate.meme.mixture_weight1: 0.5,                
                        rate.meme.omega2: 0.25,
                        rate.meme.mixture_weight2: 0.5,
                        rate.meme.omega3: 1,
                        rate.meme.mixture_weight3: 0.5
                   },
                    "2" : {
                        rate.fel.alpha : ^rate.fel.alpha,
                        rate.fel.beta_bg : ^rate.fel.beta_bg,
                        rate.fel.beta_fg: ^rate.fel.beta_fg * 4,
                        rate.meme.omega1: 0.25,
                        rate.meme.mixture_weight1: 0.25,                
                        rate.meme.omega2: 0.5,
                        rate.meme.mixture_weight2: 0.25,
                        rate.meme.omega3: 1,
                        rate.meme.mixture_weight3: 0.5
                    },
                    "3" : {
                        rate.fel.alpha : ^rate.fel.alpha,
                        rate.fel.beta_bg : ^rate.fel.beta_bg,
                        rate.fel.beta_fg: ^rate.fel.beta_fg,
                        rate.meme.omega1: 0.5,
                        rate.meme.mixture_weight1: 0.5,                
                        rate.meme.omega2: 0.75,
                        rate.meme.mixture_weight2: 0.5,
                        rate.meme.omega3: 1,
                        rate.meme.mixture_weight3: 0.5
                    },
                    "4" : {
                        rate.fel.alpha : ^rate.fel.alpha,
                        rate.fel.beta_bg : ^rate.fel.beta_bg,
                        rate.fel.beta_fg: ^rate.fel.beta_fg,
                        rate.meme.omega1: 0.25,
                        rate.meme.mixture_weight1: 0.6,               
                        rate.meme.omega2: 0.75,
                        rate.meme.mixture_weight2: 0.25,
                        rate.meme.omega3: 1,
                        rate.meme.mixture_weight3: 0.25
                   },
                    "5" : {
                        rate.fel.alpha : ^rate.fel.alpha,
                        rate.fel.beta_bg : ^rate.fel.beta_bg,
                        rate.fel.beta_fg: ^rate.fel.beta_fg * 8,
                        rate.meme.omega1: 0.1,
                        rate.meme.mixture_weight1: 0.6,             
                        rate.meme.omega2: 0.5,
                        rate.meme.mixture_weight2: 0.25,
                        rate.meme.omega3: 1,
                        rate.meme.mixture_weight3: 0.5
                    },
                    "6" : {
                        rate.fel.alpha : ^rate.fel.alpha,
                        rate.fel.beta_bg : ^rate.fel.beta_bg,
                        rate.fel.beta_fg: ^rate.fel.beta_fg,
                        rate.meme.omega1: 0,
                        rate.meme.mixture_weight1: 0.01,              
                        rate.meme.omega2: 0.5,
                        rate.meme.mixture_weight2: 0.10,
                        rate.meme.omega3: 1,
                        rate.meme.mixture_weight3: 0.10
                    },
                    "7" : {
                        rate.fel.alpha : ^rate.fel.alpha,
                        rate.fel.beta_bg : ^rate.fel.beta_bg,
                        rate.fel.beta_fg: ^rate.fel.beta_fg,
                        rate.meme.omega1: 0,
                        rate.meme.mixture_weight1: 0.7,              
                        rate.meme.omega2: 0.1,
                        rate.meme.mixture_weight2: 0.7,
                        rate.meme.omega3: 0.25,
                        rate.meme.mixture_weight3: 0.7
                    }
    
            };    
            
      }
    }
    
             
    //before_opt = {"alpha" : ^"meme.site_alpha", "other" : initial_guess_grid};
    
                 
                      
    ^rate.fel.alpha  = ^rate.fel.alpha;   
    // SLKP 20201028 : without this, there are occasional initialization issues with 
    // the likelihood function below
    
                      
    Optimize (results, ^lf_bsrel, {
            "OPTIMIZATION_METHOD" : "nedler-mead",
            "OPTIMIZATION_START_GRID" : initial_guess_grid
            //"OPTIMIZATION_PRECISION" : 1e-5,
        });
 
        
    /*after_opt = {
                    "alpha" : Eval("meme.site_alpha"), 
                    "meme.site_beta_plus": Eval("meme.site_beta_plus"),
                    "meme.site_omega_minus": Eval("meme.site_omega_minus"),
                    "meme.site_mixture_weight": Eval("meme.site_mixture_weight")
                };*/
        
    //console.log (^"LF_INITIAL_GRID_MAXIMUM_VALUE");
    //console.log (^"LF_INITIAL_GRID_MAXIMUM");
        
    if (sim_mode) {
        lrt = 2*results[1][0];    
    } else {
        alternative = estimators.ExtractMLEsOptions (lf_bsrel, model_mapping, {^"terms.globals_only" : TRUE});
        alternative [utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];
        
        
        /*
        Export (lfe, ^lf_bsrel);
        console.log (lfe);
        */
        
        ancestral_info = ancestral.build (lf_bsrel,0,FALSE);
        //TODO
        branch_substitution_information = (ancestral.ComputeSubstitutionBySite (ancestral_info,0,None))[^"terms.substitutions"];
        compressed_substitution_info = (ancestral.ComputeCompressedSubstitutions (ancestral_info))[0];
        DeleteObject (ancestral_info);
        branch_ebf       = {};
        branch_posterior = {};
        
        site_match = selection.io.sitelist_matches_pattern (pattern_info[^"terms.data.sites"], (^"meme.site_filter")["site-save-filter"], TRUE);
    
        if (site_match) {
            Export  (lfe, ^lf_bsrel);
            fprintf (^"meme.output_file_path" + "_site_" + site_match + ".fit", CLEAR_FILE, lfe);
            DeleteObject (lfe);
        }
    }
    
    
    positive_weight_expression = (scaler_mapping["OMEGA_DIST"])[rate_count-1];
    positive_weight_value = utility.EvalInGlobalNamespace (positive_weight_expression);
    
    
    if (!sim_mode) {
        if (^"meme.impute_states") {
            DataSet anc = ReconstructAncestors ( ^lf_bsrel, {{0}}, MARGINAL, DOLEAVES);
            GetString   (names, anc, -1);
            GetDataInfo (codon_chars, ^((lfInfo["Datafilters"])[0]) , "CHARACTERS");
            GetString   (names_obs, ^((lfInfo["Datafilters"])[0]), -1);
            names_obs = utility.MatrixToDict (names_obs);
        
    
            character_map = {};
            for (seq_id, seq_name; in; names) {
                GetDataInfo (site_map, ^((lfInfo["Datafilters"])[0]), names_obs[seq_name], 0);
                character_map [seq_name] = {'observed' :{} , 'imputed' : {}, 'support' : 0};
        
            
                for (char, char_support; in; site_map) {
                    if (char_support > 1e-6) {
                        (((character_map [seq_name]))['observed'])[codon_chars[char]] = char_support;
                    }
                }
            
                for (char, char_support; in; (anc.marginal_support_matrix)[seq_id][-1]) {
                    if (char_support > 1e-6) {
                        (((character_map [seq_name]))['imputed'])[codon_chars[char]] = char_support;
                        if ( (((character_map [seq_name]))['observed'])[codon_chars[char]]) {
                            (character_map [seq_name])['support'] += char_support;
                        }
                    }
                }
            
            }
        }
    }
 
    if (^rate.fel.beta_fg > ^rate.fel.alpha && positive_weight_value > 1e-6) {
        if (!sim_mode) { 
            LFCompute (^lf_bsrel,LF_START_COMPUTE);
            LFCompute (^lf_bsrel,baseline);
     

            for (_node_name_; in; ^bsrel_tree_id) {
                if (((^"meme.selected_branches") [partition_index])[_node_name_]  == utility.getGlobalValue("terms.tree_attributes.test")) {
                    _node_name_res_ = meme.compute_branch_EBF (lf_bsrel, bsrel_tree_id, _node_name_, baseline, scaler_mapping);
                    branch_ebf[_node_name_] = _node_name_res_[utility.getGlobalValue("terms.empirical_bayes_factor")];
                    branch_posterior[_node_name_] = _node_name_res_[utility.getGlobalValue("terms.posterior")];
                }
            }
        
 
            LFCompute (^lf_bsrel,LF_DONE_COMPUTE);
        }

        if (^rate.fel.alpha == 0) {
            ^rate.fel.alpha = 1e-4;
        }

        ^rate.fel.beta_fg := ^rate.fel.alpha;
        Optimize (results, ^lf_bsrel, {"OPTIMIZATION_METHOD" : "nedler-mead"});

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
                   is = meme.handle_a_site (lf_fel, lf_bsrel, sims[i], partition_index, pattern_info, model_mapping, scaler_mapping, TRUE);
                   null_LRT[i] = is;
                }
                
                return {"fel" : fel,
                    utility.getGlobalValue("terms.alternative") : alternative,
                    utility.getGlobalValue("terms.posterior") : branch_posterior,
                    utility.getGlobalValue("terms.empirical_bayes_factor") : branch_ebf,
                    utility.getGlobalValue("terms.branch_selection_attributes") : branch_substitution_information, //TODO: keep this attr?
                    utility.getGlobalValue("terms.Null"): Null,
                    utility.getGlobalValue("terms.substitutions") : compressed_substitution_info,
                    utility.getGlobalValue("terms.simulated"): null_LRT,
                    utility.getGlobalValue("terms.json.imputed_states"): character_map
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
            utility.getGlobalValue("terms.Null"): Null,
            utility.getGlobalValue("terms.json.imputed_states"): character_map};
}

/* echo to screen calls */

//----------------------------------------------------------------------------------------
function meme.report.echo (meme.report.site, meme.report.partition, meme.report.row) {

    meme.print_row = None;
    if (meme.report.row [meme.report.p_value_index] <= meme.pvalue) {
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
    
    ^"meme.site.distribution_string" = "";

    partition_index = arguments [3];
    pattern_info    = arguments [4];
    scalers         = arguments[6];
    n_rates         = utility.Array1D (scalers['OMEGA_DIST']);
    
    result_row = {};
    result_row + 0; // alpha
    result_row + 0; // beta_0
    result_row + 1; // weight 1
    for (i=1;i<n_rates;i+=1) {
        result_row + 0; // beta_i
        result_row + 0; // weight_i
    }
    
    lrt_index = Abs (result_row);
    
    result_row + 0; // LRT
    result_row + 1; // p-value
    result_row + 0; // branch count
    result_row + 0; // total branch length
    result_row + 0; // log L | FEL
    result_row + 0; // log L | MEME
    result_row + 1; // LRT MEME vs FEL
 
    result_row + 0; // FEL alpha
    result_row + 0; // FEL beta
   
    has_delta = scalers['BG'] / ^"terms.parameters.multiple_hit_rate";
    has_psi = scalers['BG'] / ^"terms.parameters.triple_hit_rate";
    
    if (has_delta) {
        has_delta = Abs (result_row);
        result_row + 0;
    }
    
    if (has_psi) {
        has_psi = Abs (result_row);
        result_row + 0;
    }
    
     fel_index = Abs (result_row);
     
    
      //console.log ( estimators.GetGlobalMLE (result["alternative"], ^"meme.parameter_site_mixture_weight"));

    if (None != result) { // not a constant site

        lrt = {utility.getGlobalValue("terms.LRT") : 2*((result[utility.getGlobalValue("terms.alternative")])[utility.getGlobalValue("terms.fit.log_likelihood")]-(result[utility.getGlobalValue("terms.Null")])[utility.getGlobalValue("terms.fit.log_likelihood")])};
        
        // TODO: adjust mixtures for different rate counts
        lrt [utility.getGlobalValue("terms.p_value")] = 2/3-2/3*(0.45*CChi2(lrt[utility.getGlobalValue("terms.LRT")],1)+0.55*CChi2(lrt[utility.getGlobalValue("terms.LRT")],2));
        
        has_lrt = result / ^"terms.simulated";
        
        if (has_lrt) {
           pv = +((result[^"terms.simulated"])["_MATRIX_ELEMENT_VALUE_>=" + lrt [utility.getGlobalValue("terms.LRT")]]);
           lrt [utility.getGlobalValue("terms.p_value")] = (pv+1)/(1+^"meme.resample");
        }     
        

        result_row [0] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], ^"terms.meme.fg_param_prefix" + ^"terms.parameters.synonymous_rate");
        
        rates = "";
        weights = "";
        for (is,w; in; scalers['OMEGA_DIST']) {
           i = +is;
           if (i < n_rates - 1) {
                result_row [1+2*i] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], ^"terms.meme.fg_param_prefix" + terms.AddCategory (^"terms.parameters.nonsynonymous_rate",i+1)) * result_row[0];
           } else {
                result_row [1+2*i] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], ^"terms.meme.fg_param_prefix" + terms.AddCategory (^"terms.parameters.nonsynonymous_rate",i+1));
           }
           result_row[2+2*i] = utility.EvalInGlobalNamespace(w);
           rate_desc = Format (result_row [1+2*i], 5,2) + " (" + Format (result_row [2+2*i]*100, 5,1 ) + "%)";
           if (i > 0) {
                rates += "/" + Format (result_row [1+2*i], 0,2);
                weights += "/" + Format (result_row [2+2*i], 0,2);
           } else {
                rates = Format (result_row [1+2*i], 0,2);
                weights = Format (result_row [2+2*i], 0,2);
           }
           
            ^"meme.site.distribution_string" = rates + " : " + weights;
           
        }
        
        //console.log (result_row);
        
        result_row [lrt_index] = lrt [utility.getGlobalValue("terms.LRT")];
        result_row [lrt_index+1] = lrt [utility.getGlobalValue("terms.p_value")];
        result_row [lrt_index+4] = (result["fel"])[utility.getGlobalValue("terms.fit.log_likelihood")];
        result_row [lrt_index+5] = (result[utility.getGlobalValue("terms.alternative")])[utility.getGlobalValue("terms.fit.log_likelihood")];
        result_row [lrt_index+6] = 1-CChi2 (2*(result_row [lrt_index+5]-result_row [lrt_index+4]),2);

        all_ebf = result[utility.getGlobalValue("terms.empirical_bayes_factor")];
        all_posterior = result[utility.getGlobalValue("terms.posterior")];

        filtered_ebf = utility.Filter (utility.Filter (all_ebf, "_value_", "None!=_value_"), "_value_", "_value_>=100");

        if(None != filtered_ebf) {
            result_row [lrt_index+2] = utility.Array1D(filtered_ebf);
        } else {
            result_row [lrt_index+2] = 0;
        }

        result_row [lrt_index+7]   =  estimators.GetGlobalMLE (result["fel"], ^"terms.meme.bg_param_prefix" + ^"terms.parameters.synonymous_rate");
        result_row [lrt_index+8] =  estimators.GetGlobalMLE (result["fel"], ^"terms.meme.fg_param_prefix" + terms.AddCategory (^"terms.parameters.nonsynonymous_rate",n_rates));
        
 
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

        result_row [lrt_index+3] = sum;
        
        ^"meme.site_multihit_string" = "";
        if (has_delta) {
            result_row[has_delta] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")],^"terms.meme.fg_param_prefix" +  ^"terms.parameters.multiple_hit_rate");
            if (None == result_row[has_delta]) {
                result_row[has_delta] = ((scalers['BG'])[^"terms.parameters.multiple_hit_rate"])["scaler"];
            }
            ^"meme.site_multihit_string" = Format (result_row[has_delta],0,2);
        }

        if (has_psi) {
            result_row[has_psi] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")],^"terms.meme.fg_param_prefix"+  ^"terms.parameters.triple_hit_rate");
            if (None == result_row[has_psi]) {
                result_row[has_delta] = ((scalers['BG'])[^"terms.parameters.triple_hit_rate"])["scaler"];
            }
             ^"meme.site_multihit_string" += "/" + Format (result_row[has_psi],0,2);
        }
        
    } else {
        all_ebf = None;
    }

    utility.EnsureKey (^"meme.site_results", partition_index);
    utility.EnsureKey (((^"meme.json")[^"terms.substitutions"]), partition_index);
    if (has_lrt) {
        utility.EnsureKey (^"meme.site_LRT", partition_index);
    }
    if (^"meme.impute_states") {
        utility.EnsureKey (^"meme.imputed_leaf_states", partition_index);
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
        if (^"meme.impute_states") {
            ((^"meme.imputed_leaf_states")[partition_index])[site_index] = result[^"terms.json.imputed_states"]
        }
        meme.report.echo (site_index, partition_index, result_row);
   }


}

