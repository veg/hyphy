RequireVersion("2.4.0");

/*------------------------------------------------------------------------------
    Load library files
*/

LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");
LoadFunctionLibrary("libv3/all-terms.bf");

LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");

LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");

LoadFunctionLibrary("libv3/convenience/math.bf");

LoadFunctionLibrary("modules/io_functions.ibf");


/*------------------------------------------------------------------------------ Display analysis information
*/

fel.analysis_description = {
    terms.io.info: "FEL (Fixed Effects Likelihood)
    estimates site-wise synonymous (&alpha;) and non-synonymous (&beta;) rates, and
    uses a likelihood ratio test to determine if beta &neq; alpha at a site.
    The estimates aggregate information over all branches,
    so the signal is derived from
    pervasive diversification or conservation. A subset of branches can be selected
    for testing as well, in which case an additional (nuisance) parameter will be
    inferred -- the non-synonymous rate on branches NOT selected for testing.
    Multiple partitions within a NEXUS file are also supported
    for recombination - aware analysis.
    ",
    terms.io.version: "2.1",
    terms.io.reference: "Not So Different After All: A Comparison of Methods for Detecting Amino Acid Sites Under Selection (2005). _Mol Biol Evol_ 22 (5): 1208-1222",
    terms.io.authors: "Sergei L Kosakovsky Pond and Simon DW Frost",
    terms.io.contact: "spond@temple.edu",
    terms.io.requirements: "in-frame codon alignment and a phylogenetic tree"
};
io.DisplayAnalysisBanner(fel.analysis_description);


/*------------------------------------------------------------------------------
    Environment setup
*/

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);


/*------------------------------------------------------------------------------
    Globals
*/


fel.site_alpha = "Site relative synonymous rate";
fel.site_beta = "Site relative non-synonymous rate (tested branches)";
fel.site_beta_nuisance = "Site relative non-synonymous rate (untested branches)";

// default cutoff for printing to screen
fel.p_value = 0.1;
fel.scaler_prefix = "FEL.scaler";

// The dictionary of results to be written to JSON at the end of the run
fel.json = {
    terms.json.analysis: fel.analysis_description,
    terms.json.input: {},
    terms.json.fits: {},
    terms.json.timers: {},
};

fel.display_orders =   {terms.original_name: -1,
                        terms.json.nucleotide_gtr: 0,
                        terms.json.global_mg94xrev: 1
                       };



selection.io.startTimer (fel.json [terms.json.timers], "Total time", 0);


/*------------------------------------------------------------------------------
    Key word arguments
*/
KeywordArgument ("code", "Which genetic code should be used", "Universal");
KeywordArgument ("alignment", "An in-frame codon alignment in one of the formats supported by HyPhy");
KeywordArgument ("tree", "A phylogenetic tree (optionally annotated with {})", null, "Please select a tree file for the data:");
KeywordArgument ("branches",  "Branches to test", "All");
KeywordArgument ("srv", "Include synonymous rate variation in the model", "Yes");
KeywordArgument ("pvalue",  "The p-value threshold to use when testing for selection", "0.1");
KeywordArgument ("ci",  "Compute profile likelihood confidence intervals for each variable site", "No");
// One additional KeywordArgument ("output") is called below after namespace fel.





namespace fel {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ("fel");
}

/* Prompt for one-rate or two-rate analysis */
fel.srv = io.SelectAnOption( {{"Yes", "[Recommended] Consider synonymous rate variation (dS varies across sites)."}, {"No", "Ignore synonymous rate variation (dS := 1 at each site)."}},
                                  "Use synonymous rate variation? Strongly recommended YES for selection inference.");
console.log(fel.srv);
if (fel.srv == "Yes"){
    fel.srv = TRUE
} else {
    fel.srv = FALSE
}
/* Prompt for p value threshold */
fel.pvalue  = io.PromptUser ("\n>Select the p-value threshold to use when testing for selection",0.1,0,1,FALSE);

fel.ci = "Yes" == io.SelectAnOption( {
    {"Yes", "Compute profile likelihood confidence intervals for dN/dS at each site (~2x slower)."}, 
    {"No", "[Default] Do not compute likelihood confidence intervals"}},
    "Compute profile likelihood intervals for site-level dN/dS");

KeywordArgument ("resample",  "[Advanced setting, will result in MUCH SLOWER run time] Perform parametric bootstrap resampling to derive site-level null LRT distributions up to this many replicates per site. Recommended use for small to medium (<30 sequences) datasets", "0");
fel.resample  = io.PromptUser ("\n>[Advanced setting, will result in MUCH SLOWER run time] Perform parametric bootstrap resampling to derive site-level null LRT distributions up to this many replicates per site. Recommended use for small to medium (<30 sequences) datasets",50,0,1000,TRUE);

KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'FEL.json')", fel.codon_data_info [terms.json.json]);
fel.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");


if (^"fel.ci") {
    fel.table_headers = {
                             {"alpha", "Synonymous substitution rate at a site"}
                             {"beta", "Non-synonymous substitution rate at a site"}
                             {"alpha=beta", "The rate estimate under the neutral model"}
                             {"LRT", "Likelihood ratio test statistic for beta = alpha, versus beta &neq; alpha"}
                             {"p-value", "Asymptotic p-value for evidence of selection, i.e. beta &neq; alpha"}
                             {"Total branch length", "The total length of branches contributing to inference at this site, and used to scale dN-dS"}
                             {"dN/dS LB", "95% profile likelihood CI lower bound for dN/dS (if available)"}
                             {"dN/dS MLE", "Point estimate for site dN/dS"}
                             {"dN/dS UB", "95% profile likelihood CI upper bound for dN/dS (if available)"}
                         };
    fel.table_screen_output  = {{"Codon", "Partition", "alpha", "beta", "LRT", "Selection detected?","dN/dS with confidence intervals"}};
    fel.table_output_options = {terms.table_options.header : TRUE, terms.table_options.minimum_column_width: 16, terms.table_options.align : "center"};

} else {
    fel.table_headers = {{"alpha", "Synonymous substitution rate at a site"}
                         {"beta", "Non-synonymous substitution rate at a site"}
                         {"alpha=beta", "The rate estimate under the neutral model"}
                         {"LRT", "Likelihood ratio test statistic for beta = alpha, versus beta &neq; alpha"}
                         {"p-value", "Asymptotic p-value for evidence of selection, i.e. beta &neq; alpha"}
                         {"Total branch length", "The total length of branches contributing to inference at this site, and used to scale dN-dS"}};
    fel.table_screen_output  = {{"Codon", "Partition", "alpha", "beta", "LRT", "Selection detected?"}};
    fel.table_output_options = {terms.table_options.header : TRUE, terms.table_options.minimum_column_width: 16, terms.table_options.align : "center"};

}




fel.partition_count = Abs (fel.filter_specification);
io.ReportProgressMessageMD('FEL',  'selector', 'Branches to include in the FEL analysis');



utility.ForEachPair (fel.selected_branches, "_partition_", "_selection_",
    "_selection_ = utility.Filter (_selection_, '_value_', '_value_ == terms.tree_attributes.test');
     io.ReportProgressMessageMD('FEL',  'selector', 'Selected ' + Abs(_selection_) + ' branches to include in FEL calculations: \\\`' + Join (', ',utility.Keys(_selection_)) + '\\\`')");


selection.io.startTimer (fel.json [terms.json.timers], "Model fitting",1);

namespace_tag = "fel";

KeywordArgument ("precision", "Optimization precision settings for preliminary fits", "standard");
fel.faster = io.SelectAnOption( {{"standard", "[Default] Use standard optimization settings."}, 
                                 {"reduced", "Cruder optimizations settings for faster fitting of preliminary models"}},
                                  "Optimization precision settings for preliminary fits") == "reduced";


namespace fel {
    if (faster) {
        utility.ToggleEnvVariable("OPTIMIZATION_PRECISION", 0.1);
    }
    doGTR ("fel");
    if (faster) {
        utility.ToggleEnvVariable("OPTIMIZATION_PRECISION", None);
    }
}

estimators.fixSubsetOfEstimates(fel.gtr_results, fel.gtr_results[terms.global]);


namespace fel {
    if (faster) {
        utility.ToggleEnvVariable("OPTIMIZATION_PRECISION", Abs(0.0001*gtr_results[^"terms.fit.log_likelihood"]));
    }
    doPartitionedMG ("fel", FALSE);
    if (faster) {
        utility.ToggleEnvVariable("OPTIMIZATION_PRECISION", None);
    }
}

KeywordArgument ("full-model", "Perform branch length re-optimization under the full codon model", "Yes");

fel.run_full_mg94 = "Yes" == io.SelectAnOption( {{"Yes", "[Default] Perform branch length re-optimization under the full codon model"}, 
                                         {"No", "Skip branch length re-optimization under the full codon model (faster but less precise)"}},
                                         "Perform branch length re-optimization under the full codon model");


if (fel.run_full_mg94) {
    io.ReportProgressMessageMD ("fel", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");

    
    if (Type (fel.save_intermediate_fits) == "AssociativeList") {
        if (None != fel.save_intermediate_fits[^"terms.data.value"]) {
            if (utility.Has (fel.save_intermediate_fits[^"terms.data.value"], "Full-MG94", "AssociativeList")) {
                fel.final_partitioned_mg_results = (fel.save_intermediate_fits[^"terms.data.value"])["Full-MG94"];
                if (utility.Has (fel.save_intermediate_fits[^"terms.data.value"], "Full-MG94-LF", "String")) {
                    ExecuteCommands ((fel.save_intermediate_fits[^"terms.data.value"])["Full-MG94-LF"]);
                    fel.run_full_mg94 = FALSE;
                }
            }        
        }
    }
    
    if (fel.run_full_mg94) {    

        if (fel.faster) {
            utility.ToggleEnvVariable("OPTIMIZATION_PRECISION", Abs(0.0001*fel.partitioned_mg_results[^"terms.fit.log_likelihood"]));
        }
        fel.final_partitioned_mg_results = estimators.FitMGREV (fel.filter_names, fel.trees, fel.codon_data_info [terms.code], {
            terms.run_options.model_type: terms.local,
            terms.run_options.partitioned_omega: fel.selected_branches,
            terms.run_options.retain_lf_object: TRUE,
            terms.run_options.apply_user_constraints: fel.zero_branch_length_constrain,
            terms.run_options.optimization_settings: {
                "OPTIMIZATION_METHOD" : "coordinate-wise"
            }
        }, fel.partitioned_mg_results);

        if (faster) {
            utility.ToggleEnvVariable("OPTIMIZATION_PRECISION", None);
        }
    
        if (Type (fel.save_intermediate_fits) == "AssociativeList") {
            (fel.save_intermediate_fits[^"terms.data.value"])["Full-MG94"] = fel.final_partitioned_mg_results;        
            Export (lfe, ^fel.final_partitioned_mg_results[^"terms.likelihood_function"]);
            (fel.save_intermediate_fits[^"terms.data.value"])["Full-MG94-LF"] = lfe;
            io.SpoolJSON (fel.save_intermediate_fits[^"terms.data.value"],fel.save_intermediate_fits[^"terms.data.file"]);      
        }
    }
} else {
     fel.final_partitioned_mg_results = fel.partitioned_mg_results;
     estimators.RemoveBranchLengthConstraints (fel.final_partitioned_mg_results);
     estimators.RemoveGlobalConstraints (fel.final_partitioned_mg_results);
}

fel.save_intermediate_fits = None;


io.ReportProgressMessageMD("fel", "codon-refit", "* Log(L) = " + Format(fel.final_partitioned_mg_results[terms.fit.log_likelihood],8,2));
fel.global_dnds = selection.io.extract_global_MLE_re (fel.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio);
utility.ForEach (fel.global_dnds, "_value_", 'io.ReportProgressMessageMD ("fel", "codon-refit", "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));');



estimators.fixSubsetOfEstimates(fel.final_partitioned_mg_results, fel.final_partitioned_mg_results[terms.global]);

//Store MG94 to JSON
selection.io.json_store_lf_withEFV (fel.json,
                            terms.json.global_mg94xrev,
                            fel.final_partitioned_mg_results[terms.fit.log_likelihood],
                            fel.final_partitioned_mg_results[terms.parameters],
                            fel.sample_size,
                            utility.ArrayToDict (utility.Map (fel.global_dnds, "_value_", "{'key': _value_[terms.description], 'value' : Eval({{_value_ [terms.fit.MLE],1}})}")),
                            (fel.final_partitioned_mg_results[terms.efv_estimate])["VALUEINDEXORDER"][0],
                            fel.display_orders[terms.json.global_mg94xrev]);

utility.ForEachPair (fel.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(fel.json, terms.json.global_mg94xrev, terms.branch_length, fel.display_orders[terms.json.global_mg94xrev],
                                             _key_,
                                             selection.io.extract_branch_info((fel.final_partitioned_mg_results[terms.branch_length])[_key_], "selection.io.branch.length"));');



selection.io.stopTimer (fel.json [terms.json.timers], "Model fitting");

// define the site-level likelihood function

fel.site.mg_rev = model.generic.DefineModel("models.codon.MG_REV.ModelDescription",
        "fel_mg", {
            "0": parameters.Quote(terms.local),
            "1": fel.codon_data_info[terms.code]
        },
        fel.filter_names,
        None);


fel.site_model_mapping = {"fel_mg" : fel.site.mg_rev};

/* set up the local constraint model */

fel.alpha = model.generic.GetLocalParameter (fel.site.mg_rev, utility.getGlobalValue("terms.parameters.synonymous_rate"));
fel.beta = model.generic.GetLocalParameter (fel.site.mg_rev, utility.getGlobalValue("terms.parameters.nonsynonymous_rate"));
io.CheckAssertion ("None!=fel.alpha && None!=fel.beta", "Could not find expected local synonymous and non-synonymous rate parameters in \`estimators.FitMGREV\`");

selection.io.startTimer (fel.json [terms.json.timers], "FEL analysis", 2);

//----------------------------------------------------------------------------------------
function fel.apply_proportional_site_constraint (tree_name, node_name, alpha_parameter, beta_parameter, alpha_factor, beta_factor, branch_length) {

    fel.branch_length = (branch_length[terms.parameters.synonymous_rate])[terms.fit.MLE];
    
    node_name = tree_name + "." + node_name;

    ExecuteCommands ("
        `node_name`.`alpha_parameter` := (`alpha_factor`) * fel.branch_length__;
        `node_name`.`beta_parameter`  := (`beta_factor`)  * fel.branch_length__;
    ");
}
//----------------------------------------------------------------------------------------

fel.scalers = {{"fel.alpha_scaler", "fel.beta_scaler_test", "fel.beta_scaler_nuisance"}};

model.generic.AddGlobal (fel.site.mg_rev, "fel.alpha_scaler", fel.site_alpha);
model.generic.AddGlobal (fel.site.mg_rev, "fel.beta_scaler_test", fel.site_beta);
model.generic.AddGlobal (fel.site.mg_rev, "fel.beta_scaler_nuisance", fel.site_beta_nuisance);
parameters.DeclareGlobal (fel.scalers, {});


//----------------------------------------------------------------------------------------
lfunction fel.handle_a_site (lf, filter_data, partition_index, pattern_info, model_mapping, sim_mode) {
  
    GetString (lfInfo, ^lf,-1);
    ExecuteCommands (filter_data); 
    __make_filter ((lfInfo["Datafilters"])[0]);

 
    utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);

    if (^"fel.srv"){
        ^"fel.alpha_scaler" = 1;
        start.grid = {
            "0" : {
                "fel.beta_scaler_test": 0.1,
                "fel.beta_scaler_nuisance" : 0.1,
                "fel.alpha_scaler" : 0.01
            },
            "1" : {
                "fel.beta_scaler_test": 0.1,
                "fel.beta_scaler_nuisance" : 0.1,
                "fel.alpha_scaler" : 1
            },
            "2" : {
                "fel.beta_scaler_test": 0.5,
                "fel.beta_scaler_nuisance" : 0.5,
                "fel.alpha_scaler" : 1
            },
            "3" : {
                "fel.beta_scaler_test": 1.,
                "fel.beta_scaler_nuisance" : 1.,
                "fel.alpha_scaler" : 1
            },
            "4" : {
                "fel.beta_scaler_test": 5.,
                "fel.beta_scaler_nuisance" : 5.,
                "fel.alpha_scaler" : 1
            },
            "5" : {
                "fel.beta_scaler_test": 0.1,
                "fel.beta_scaler_nuisance" : 0.1,
                "fel.alpha_scaler" : 10.
            }
        };
    } else {
        ^"fel.alpha_scaler" := 1;
        start.grid = {
            "0" : {
                "fel.beta_scaler_test": 0.01,
                "fel.beta_scaler_nuisance" : 0.01 
            },
            "1" : {
                "fel.beta_scaler_test": 0.1,
                "fel.beta_scaler_nuisance" : 0.1 
            },
            "2" : {
                "fel.beta_scaler_test": 0.25,
                "fel.beta_scaler_nuisance" : 1 
            },
            "3" : {
                "fel.beta_scaler_test": 0.5,
                "fel.beta_scaler_nuisance" : 0.5
            },
            "4" : {
                "fel.beta_scaler_test": 1.0,
                "fel.beta_scaler_nuisance" : 1.0
            },
            "5" : {
                "fel.beta_scaler_test": 5.0,
                "fel.beta_scaler_nuisance" : 5.0
            }
        };
    }
    ^"fel.beta_scaler_test"  = 1;
    ^"fel.beta_scaler_nuisance"  = 1;
    

    Optimize (results, ^lf
        , {
            "OPTIMIZATION_METHOD" : "nedler-mead", 
            "OPTIMIZATION_PRECISION" : 1e-5,
            "OPTIMIZATION_START_GRID" : start.grid             
           }
    );
    
    
    if (^"fel.ci") {
        if (!sim_mode) {
            if (^"fel.srv") {
                lf_stash = estimators.TakeLFStateSnapshot (lf);
                global omega_ratio_for_ci;
                if (^"fel.alpha_scaler" == 0.0) { // has no syn_variation
                    ^"fel.alpha_scaler" = 1e-8;
                }
                omega_ratio_for_ci :< 10000;
                omega_ratio_for_ci :>0;
                omega_ratio_for_ci = ^"fel.beta_scaler_test" / ^"fel.alpha_scaler";
            
                ^"fel.beta_scaler_test" := omega_ratio_for_ci * ^"fel.alpha_scaler";
                ci = parameters.GetProfileCI (&omega_ratio_for_ci, lf, 0.95);
                estimators.RestoreLFStateFromSnapshot (lf, lf_stash);
            } else {
                ci = parameters.GetProfileCI ("fel.beta_scaler_test", lf, 0.95);
            }
        }
    } else {
        ci = None;
    }
 
    if (sim_mode) {
        lrt = 2*results[1][0];    
    } else {
        alternative = estimators.ExtractMLEsOptions (lf, model_mapping, {^"terms.globals_only" : TRUE});
        alternative [utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];
    }

    ^"fel.alpha_scaler" = (^"fel.alpha_scaler" + 3*^"fel.beta_scaler_test")/4;
    parameters.SetConstraint ("fel.beta_scaler_test","fel.alpha_scaler", "");

    Optimize (results, ^lf);

    if (sim_mode) {
        return lrt - 2*results[1][0];
    } else {
        Null = estimators.ExtractMLEsOptions (lf, model_mapping, {^"terms.globals_only" : TRUE});
        Null [utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];
    }    


    if (!sim_mode) {
        tree_name = (lfInfo["Trees"])[0];
        sum = (BranchLength (^tree_name, -1)*^"fel.selected_branches_index")[0];
        if (^"fel.resample") {
            N = ^"fel.resample";
            sims = {};
            GetDataInfo (fi, ^((lfInfo["Datafilters"])[0]), "PARAMETERS");
            //Export (lfe, ^lf);
            //console.log (lfe);
            for (i = 0; i < N; i+=1) {        
                DataSet null_sim = SimulateDataSet (^lf);
                DataSetFilter null_filter = CreateFilter (null_sim,3,,,fi["EXCLUSIONS"]);
                sims + alignments.serialize_site_filter (&null_filter, 0);
            }
            null_LRT = {N,1};
            for (i = 0; i < N; i+=1) {   
               is = fel.handle_a_site (lf, sims[i], partition_index, pattern_info, model_mapping, TRUE);
               null_LRT[i] = is;
            }
            return {utility.getGlobalValue("terms.alternative") : alternative, utility.getGlobalValue("terms.Null"): Null, utility.getGlobalValue("terms.simulated"): null_LRT, utility.getGlobalValue("terms.confidence_interval"): ci};
        }
    }
    
 
    return {utility.getGlobalValue("terms.alternative") : alternative, 
            utility.getGlobalValue("terms.Null"): Null,
            utility.getGlobalValue("terms.math.sum") : sum,
            utility.getGlobalValue("terms.confidence_interval"): ci};
}

/* echo to screen calls */

fel.report.counts        = {{0,0}};


if (^"fel.ci") {
    fel.report.positive_site = {{"" + (1+((fel.filter_specification[fel.report.partition])[terms.data.coverage])[fel.report.site]),
                                        fel.report.partition + 1,
                                        Format(fel.report.row[0],10,3),
                                        Format(fel.report.row[1],10,3),
                                        Format(fel.report.row[3],10,3),
                                        "Pos. p = " + Format(fel.report.row[4],6,4),
                                        Format(fel.report.row[7],7,3) + "(" + Format(fel.report.row[6],7,2) + "-" + Format(fel.report.row[8],7,2) + ")"}};

    fel.report.negative_site = {{"" + (1+((fel.filter_specification[fel.report.partition])[terms.data.coverage])[fel.report.site]),
                                        fel.report.partition + 1,
                                        Format(fel.report.row[0],10,3),
                                        Format(fel.report.row[1],10,3),
                                        Format(fel.report.row[3],10,3),
                                        "Neg. p = " + Format(fel.report.row[4],6,4),
                                        Format(fel.report.row[7],7,3) + "(" + Format(fel.report.row[6],7,2) + "-" + Format(fel.report.row[8],7,2) + ")"}};

} else {
    fel.report.positive_site = {{"" + (1+((fel.filter_specification[fel.report.partition])[terms.data.coverage])[fel.report.site]),
                                        fel.report.partition + 1,
                                        Format(fel.report.row[0],10,3),
                                        Format(fel.report.row[1],10,3),
                                        Format(fel.report.row[3],10,3),
                                        "Pos. p = " + Format(fel.report.row[4],6,4)}};

    fel.report.negative_site = {{"" + (1+((fel.filter_specification[fel.report.partition])[terms.data.coverage])[fel.report.site]),
                                        fel.report.partition + 1,
                                        Format(fel.report.row[0],10,3),
                                        Format(fel.report.row[1],10,3),
                                        Format(fel.report.row[3],10,3),
                                        "Neg. p = " + Format(fel.report.row[4],6,4)}};
}

function fel.report.echo (fel.report.site, fel.report.partition, fel.report.row) {

    fel.print_row = None;
    if (fel.report.row [4] < fel.pvalue) {
        if (fel.report.row[0] < fel.report.row[1]) {
            fel.print_row = fel.report.positive_site;
            fel.color = "\033[0;31m";
            fel.report.counts[0] += 1;
        } else {
            fel.print_row = fel.report.negative_site;
            fel.report.counts [1] += 1;
            fel.color = "\033[0;32m";
        }
    }

     if (None != fel.print_row) {
            io.ClearProgressBar();
            if (!fel.report.header_done) {
                io.ReportProgressMessageMD("FEL", "" + fel.report.partition, "For partition " + (fel.report.partition+1) + " these sites are significant at p <=" + fel.pvalue + "\n");
                fprintf (stdout,
                    io.FormatTableRow (fel.table_screen_output,fel.table_output_options));
                fel.report.header_done = TRUE;
                fel.table_output_options[terms.table_options.header] = FALSE;
            }

            fprintf (stdout, fel.color, 
                io.FormatTableRow (fel.print_row,fel.table_output_options),"\033[0m"
                );

        }

}


lfunction fel.store_results (node, result, arguments) {

    partition_index = arguments [2];
    pattern_info    = arguments [3];
    
    if (^"fel.ci") {
        result_row          = { { 0, // alpha
                              0, // beta
                              0, // alpha==beta
                              0, // LRT
                              1, // p-value,
                              0,  // total branch length of tested branches
                              0,   // lower bound dN/dS
                              0,   // dN/dS MLE
                              0    // upper bound dN/dS
                          } };    
    } else {
        result_row          = { { 0, // alpha
                              0, // beta
                              0, // alpha==beta
                              0, // LRT
                              1, // p-value,
                              0,  // total branch length of tested branches
                          } };
    }

    has_lrt = FALSE;

    if (None != result) { // not a constant site
    
        lrt = math.DoLRT ((result[utility.getGlobalValue("terms.Null")])[utility.getGlobalValue("terms.fit.log_likelihood")],
                          (result[utility.getGlobalValue("terms.alternative")])[utility.getGlobalValue("terms.fit.log_likelihood")],
                          1);
                    
        has_lrt = result / ^"terms.simulated";
                  
        if (has_lrt) {
           pv = +((result[^"terms.simulated"])["_MATRIX_ELEMENT_VALUE_>=" + lrt [utility.getGlobalValue("terms.LRT")]]);
           lrt [utility.getGlobalValue("terms.p_value")] = (pv+1)/(1+^"fel.resample");
        }                 
        
        result_row [0] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], ^"fel.site_alpha");
        result_row [1] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], ^"fel.site_beta");
        result_row [2] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.Null")], ^"fel.site_beta");
        result_row [3] = lrt [utility.getGlobalValue("terms.LRT")];
        result_row [4] = lrt [utility.getGlobalValue("terms.p_value")];


        
        /*alternative_lengths = ((result[utility.getGlobalValue("terms.alternative")])[utility.getGlobalValue("terms.branch_length")])[0];
        
        tname = utility.getGlobalValue("terms.tree_attributes.test");
        for (_node_; in; ^"fel.case_respecting_node_names") {
            _node_class_ = ((^"fel.selected_branches")[partition_index])[_node_];
            if (_node_class_ ==  tname) {
                sum += ((alternative_lengths)[_node_])[^"terms.fit.MLE"];
            }
        }*/

    
        result_row [5] = result[^"terms.math.sum"];
        
        if (None != result [^"terms.confidence_interval"]) {
            result_row [6] = (result [^"terms.confidence_interval"])[^"terms.lower_bound"];
            result_row [7] = (result [^"terms.confidence_interval"])[^"terms.fit.MLE"];
            result_row [8] = (result [^"terms.confidence_interval"])[^"terms.upper_bound"];
       }
    }


    utility.EnsureKey (^"fel.site_results", partition_index);
    if (has_lrt) {
        utility.EnsureKey (^"fel.site_LRT", partition_index);
    }
    
    for (_fel_result_; in; pattern_info[utility.getGlobalValue("terms.data.sites")]) {
        ((^"fel.site_results")[partition_index])[_fel_result_] = result_row;
        fel.report.echo (_fel_result_, partition_index, result_row);
        if (has_lrt) {
            ((^"fel.site_LRT")[partition_index])[_fel_result_] = result[^"terms.simulated"];
        }
    }



    //assert (0);
}
//----------------------------------------------------------------------------------------

fel.site_results = {};
fel.site_LRT = {};

for (fel.partition_index = 0; fel.partition_index < fel.partition_count; fel.partition_index += 1) {
    fel.report.header_done = FALSE;
    fel.table_output_options[terms.table_options.header] = TRUE;
    model.ApplyModelToTree( "fel.site_tree", fel.trees[fel.partition_index], {terms.default : fel.site.mg_rev}, None);

    fel.case_respecting_node_names = trees.branch_names (fel.site_tree, TRUE);

    fel.site_patterns = alignments.Extract_site_patterns ((fel.filter_specification[fel.partition_index])[utility.getGlobalValue("terms.data.name")]);

    // apply constraints to the site tree
    // alpha = alpha_scaler * branch_length
    // beta  = beta_scaler_test * branch_length or beta_nuisance_test * branch_length

    SetParameter (DEFER_CONSTRAINT_APPLICATION, 1, 0);
    fel.selected_branches_index = {1,utility.Array1D (fel.case_respecting_node_names)+1};
    fel.i = 0;

    for (_node_; in; fel.case_respecting_node_names) {
        _node_class_ = (fel.selected_branches[fel.partition_index])[_node_];
         if (_node_class_ == terms.tree_attributes.test) {
             fel.selected_branches_index [fel.i] = 1;
            _beta_scaler = fel.scalers[1];
         } else {
            _beta_scaler = fel.scalers[2];
         }
         fel.apply_proportional_site_constraint ("fel.site_tree", _node_, fel.alpha, fel.beta, fel.scalers[0], _beta_scaler, (( fel.final_partitioned_mg_results[terms.branch_length])[fel.partition_index])[_node_]);
         fel.i += 1;
    }
    
 
    /*
    utility.ForEach (fel.case_respecting_node_names, "_node_",
            '_node_class_ = (fel.selected_branches[fel.partition_index])[_node_];
             if (_node_class_ == terms.tree_attributes.test) {
                _beta_scaler = fel.scalers[1];
             } else {
                _beta_scaler = fel.scalers[2];
             }
             fel.apply_proportional_site_constraint ("fel.site_tree", _node_, fel.alpha, fel.beta, fel.scalers[0], _beta_scaler, (( fel.final_partitioned_mg_results[terms.branch_length])[fel.partition_index])[_node_]);
        ');
    */

    SetParameter (DEFER_CONSTRAINT_APPLICATION, 0, 0);


    // create the likelihood function for this site
    ExecuteCommands (alignments.serialize_site_filter
                                       ((fel.filter_specification[fel.partition_index])[utility.getGlobalValue("terms.data.name")],
                                       ((fel.site_patterns[0])[utility.getGlobalValue("terms.data.sites")])[0],
                   ));

    __make_filter ("fel.site_filter");

    LikelihoodFunction fel.site_likelihood = (fel.site_filter, fel.site_tree);



    estimators.ApplyExistingEstimates ("fel.site_likelihood", fel.site_model_mapping, fel.final_partitioned_mg_results,
                                        terms.globals_only);


    fel.queue = mpi.CreateQueue ({"LikelihoodFunctions": {{"fel.site_likelihood"}},
                                   "Models" : {{"fel.site.mg_rev"}},
                                   "Headers" : {{"libv3/all-terms.bf","libv3/tasks/alignments.bf"}},
                                   "Variables" : {{"fel.srv","fel.resample","fel.ci","fel.selected_branches_index"}}
                                 });


    /* run the main loop over all unique site pattern combinations */
    
    fel.pattern_count_all = utility.Array1D (fel.site_patterns);
    fel.pattern_count_this = 0;
    
    utility.ForEachPair (fel.site_patterns, "_pattern_", "_pattern_info_",
        '
           fel.pattern_count_this += 1;
           io.ReportProgressBar("", "Working on site pattern " + (fel.pattern_count_this) + "/" +  fel.pattern_count_all + " in partition " + (1+fel.partition_index));

           if (_pattern_info_[terms.data.is_constant]) {
                fel.store_results (-1,None,{"0" : "fel.site_likelihood",
                                                                "1" : None,
                                                                "2" : fel.partition_index,
                                                                "3" : _pattern_info_,
                                                                "4" : fel.site_model_mapping,
                                                                "5" : 0
                                                                });
            } else {
                mpi.QueueJob (fel.queue, "fel.handle_a_site", {"0" : "fel.site_likelihood",
                                                                "1" : alignments.serialize_site_filter
                                                                   ((fel.filter_specification[fel.partition_index])[terms.data.name],
                                                                   (_pattern_info_[utility.getGlobalValue("terms.data.sites")])[0]),
                                                                "2" : fel.partition_index,
                                                                "3" : _pattern_info_,
                                                                "4" : fel.site_model_mapping,
                                                                "5" : 0
                                                                },
                                                                "fel.store_results");
            }
        '
    );
    
    io.ClearProgressBar();

    mpi.QueueComplete (fel.queue);
    
    fel.partition_matrix = {Abs (fel.site_results[fel.partition_index]), Rows (fel.table_headers)};
    for (_key_, _value_; in; fel.site_results[fel.partition_index]) {
        for (fel.index = 0; fel.index < Rows (fel.table_headers); fel.index += 1) {
            fel.partition_matrix [0+_key_][fel.index] = _value_[fel.index];
        }
    }
    fel.site_results[fel.partition_index] = fel.partition_matrix;

}

fel.json [terms.json.MLE ] = {terms.json.headers   : fel.table_headers,
                               terms.json.content : fel.site_results };

if (fel.resample)  {
    (fel.json [terms.json.MLE ])[terms.LRT] = fel.site_LRT;
}

io.ReportProgressMessageMD ("fel", "results", "** Found _" + fel.report.counts[0] + "_ sites under pervasive positive diversifying and _" + fel.report.counts[1] + "_ sites under negative selection at p <= " + fel.pvalue + "**");

selection.io.stopTimer (fel.json [terms.json.timers], "Total time");
selection.io.stopTimer (fel.json [terms.json.timers], "FEL analysis");

if (fel.resample) {
    fel.json [terms.simulated] =  fel.resample;
}

if (fel.ci) {
    fel.json [terms.confidence_interval] = TRUE;
}

GetString (_hpv,HYPHY_VERSION,0);
fel.json[terms.json.runtime] = _hpv;


io.SpoolJSON (fel.json, fel.codon_data_info[terms.json.json]);


