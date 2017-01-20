RequireVersion("2.3");

/*------------------------------------------------------------------------------
    Load library files
*/

LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");
LoadFunctionLibrary("libv3/terms-json.bf");

LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");

LoadFunctionLibrary("libv3/models/codon/BS_REL.bf");

LoadFunctionLibrary("libv3/convenience/math.bf");

LoadFunctionLibrary("modules/io_functions.ibf");


/*------------------------------------------------------------------------------ Display analysis information
*/

io.DisplayAnalysisBanner({
    "info": "MEME (Mixed Effects Model of Evolution)
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
    "version": "2.00",
    "reference": "Detecting Individual Sites Subject to Episodic Diversifying Selection. _PLoS Genet_ 8(7): e1002764.",
    "authors": "Sergei L. Kosakovsky Pond, Steven Weaver",
    "contact": "spond@temple.edu",
    "requirements": "in-frame codon alignment and a phylogenetic tree"
});



/*------------------------------------------------------------------------------
    Environment setup
*/

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

/*------------------------------------------------------------------------------
    Globals
*/



meme.terms.site_alpha = "Site relative synonymous rate";
meme.terms.site_omega_minus = "Omega ratio on (tested branches); negative selection or neutral evolution (&omega;- <= 1;)";
meme.terms.site_beta_minus = "Site relative non-synonymous rate (tested branches); negative selection or neutral evolution (&beta;- <= &alpha;)";
meme.terms.site_beta_plus = "Site relative non-synonymous rate (tested branches); unconstrained";
meme.terms.site_mixture_weight = "Mixture proportion allocated to &beta;-, i.e. negative selection or neutral evolution";
meme.terms.site_beta_nuisance = "Site relative non-synonymous rate (untested branches)";

meme.pvalue = 0.1;
    /**
        default cutoff for printing to screen
    */

meme.json = {
    terms.json.fits: {},
    terms.json.timers: {},
};
    /**
        The dictionary of results to be written to JSON at the end of the run
    */

selection.io.startTimer (meme.json [terms.json.timers], "Total time", 0);
meme.scaler_prefix = "MEME.scaler";

meme.table_headers = {{"&alpha;", "Synonymous substitution rate at a site"}
                     {"&beta;<sup>-</sup>", "Non-synonymous substitution rate at a site for the negative/neutral evolution component"}
                     {"p<sup>-</sup>", "Mixture distribution weight allocated to &beta;<sup>-</sup>; loosely -- the proportion of the tree evolving neutrally or under negative selection"}
                     {"&beta;<sup>+</sup>", "Non-synonymous substitution rate at a site for the positive/neutral evolution component"}
                     {"p<sup>+</sup>", "Mixture distribution weight allocated to &beta;<sup>+</sup>; loosely -- the proportion of the tree evolving neutrally or under positive selection"}
                     {"LRT", "Likelihood ratio test statistic for episodic diversification, i.e., p<sup>+</sup> &gt; 0 <emph>and<emph> &beta;<sup>+</sup> &gt; &alpha;"}
                     {"p-value", "Asymptotic p-value for for episodic diversification, i.e., p<sup>+</sup> &gt; 0 <emph>and<emph> &beta;<sup>+</sup> &gt; &alpha;"}
                     {"Total branch length", "The total length of branches contributing to inference at this site, and used to scale dN-dS"}};
/**
    This table is meant for HTML rendering in the results web-app; can use HTML characters, the second column
    is 'pop-over' explanation of terms. This is ONLY saved to the JSON file. For Markdown screen output see
    the next set of variables.
*/



meme.table_screen_output  = {{"Codon", "Partition", "alpha", "beta+", "p+", "LRT", "Selection detected?"}};
meme.table_output_options = {"header" : TRUE, "min-column-width" : 16, "align" : "center"};

namespace meme {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ("meme");
}

/**
    TODO: Need to document what variables are available after call to load_file
*/


meme.partition_count = Abs (meme.filter_specification);
meme.pvalue  = io.PromptUser ("\n>Select the p-value used to for perform the test at",meme.pvalue,0,1,FALSE);
io.ReportProgressMessageMD('MEME',  'selector', 'Branches to include in the MEME analysis');

utility.ForEachPair (meme.selected_branches, "_partition_", "_selection_",
    "_selection_ = utility.Filter (_selection_, '_value_', '_value_ == terms.json.attribute.test');
     io.ReportProgressMessageMD('MEME',  'selector', 'Selected ' + Abs(_selection_) + ' branches to include in MEME calculations: \\\`' + Join (', ',utility.Keys(_selection_)) + '\\\`')");


//console.log (utility.DefinedVariablesMatchingRegex ("^meme"));

selection.io.startTimer (meme.json [terms.json.timers], "Model fitting",1);

namespace meme {
    doGTR ("meme");
}


estimators.fixSubsetOfEstimates(meme.gtr_results, meme.gtr_results["global"]);

namespace meme {
    doPartitionedMG ("meme", FALSE);
}


io.ReportProgressMessageMD ("MEME", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");

meme.final_partitioned_mg_results = estimators.FitMGREV (meme.filter_names, meme.trees, meme.codon_data_info ["code"], {
    "model-type": terms.global,
    "partitioned-omega": meme.selected_branches,
    "retain-lf-object": TRUE
}, meme.partitioned_mg_results);

io.ReportProgressMessageMD("MEME", "codon-refit", "* Log(L) = " + Format(meme.final_partitioned_mg_results["LogL"],8,2));
meme.global_dnds = selection.io.extract_global_MLE_re (meme.final_partitioned_mg_results, "^" + terms.omega_ratio);
utility.ForEach (meme.global_dnds, "_value_", 'io.ReportProgressMessageMD ("MEME", "codon-refit", "* " + _value_["description"] + " = " + Format (_value_["MLE"],8,4));');

estimators.fixSubsetOfEstimates(meme.final_partitioned_mg_results, meme.final_partitioned_mg_results["global"]);

selection.io.json_store_lf(
    meme.json,
    "Global MG94xREV",
    meme.final_partitioned_mg_results["LogL"],
    meme.final_partitioned_mg_results["parameters"],
    meme.sample_size,
    utility.ArrayToDict (utility.Map (meme.global_dnds, "_value_", "{'key': _value_['description'], 'value' : Eval({{_value_ ['MLE'],1}})}"))
);

utility.ForEachPair (meme.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(meme.json, "Global MG94xREV model", terms.json.attribute.branch_length, 0,
                                             _key_,
                                             selection.io.extract_branch_info((meme.final_partitioned_mg_results[terms.json.attribute.branch_length])[_key_], "selection.io.branch.length"));');

selection.io.stopTimer (meme.json [terms.json.timers], "Model fitting");



// define the site-level likelihood function

meme.site.background_fel = model.generic.DefineModel("models.codon.MG_REV.ModelDescription",
        "meme.background_fel", {
            "0": parameters.Quote(terms.local),
            "1": meme.codon_data_info["code"]
        },
        meme.filter_names,
        None);

meme.alpha = model.generic.GetLocalParameter (meme.site.background_fel, terms.synonymous_rate);
meme.beta  = model.generic.GetLocalParameter (meme.site.background_fel, terms.nonsynonymous_rate);

io.CheckAssertion ("None!=meme.alpha && None!=meme.beta", "Could not find expected local synonymous and non-synonymous rate parameters in \`estimators.FitMGREV\`");

meme.site.bsrel =  model.generic.DefineMixtureModel("models.codon.BS_REL.ModelDescription",
        "meme.bsrel_negative", {
            "0": parameters.Quote(terms.local),
            "1": meme.codon_data_info["code"],
            "2": parameters.Quote (2)
        },
        meme.filter_names,
        None);


meme.beta1  = model.generic.GetLocalParameter (meme.site.bsrel , terms.AddCategory (terms.nonsynonymous_rate,1));
meme.beta2  = model.generic.GetLocalParameter (meme.site.bsrel , terms.AddCategory (terms.nonsynonymous_rate,2));

io.CheckAssertion ("None!=meme.beta2&&None!=meme.beta1", "Could not find expected local non-synonymous rate parameters for the unrestricted BS-REL class");

meme.site_model_mapping = {
                           "meme.background_fel" : meme.site.background_fel,
                           "meme.bsrel" : meme.site.bsrel,
                         };


selection.io.startTimer (meme.json [terms.json.timers], "MEME analysis", 2);


mneme.global_variables = {
            { "meme.site_alpha",
              "meme.site_omega_minus",
              "meme.site_beta_minus",
              "meme.site_beta_plus",
              "meme.site_beta_nuisance"
             }
        };


model.generic.AddGlobal (meme.site.background_fel, "meme.site_alpha", meme.terms.site_alpha);
parameters.DeclareGlobal ("meme.site_alpha", {});


model.generic.AddGlobal (meme.site.background_fel, "meme.site_beta_nuisance", meme.terms.site_beta_nuisance);
parameters.DeclareGlobal ("meme.site_beta_nuisance", {});

model.generic.AddGlobal (meme.site.bsrel, "meme.site_alpha", meme.terms.site_alpha);
model.generic.AddGlobal (meme.site.bsrel, "meme.site_omega_minus", meme.terms.site_omega_minus);
parameters.DeclareGlobal ("meme.site_omega_minus", {});
parameters.SetRange ("meme.site_omega_minus", terms.range01);

model.generic.AddGlobal (meme.site.bsrel, "meme.site_beta_minus", meme.terms.site_beta_minus);
parameters.DeclareGlobal ("meme.site_beta_minus", {});
parameters.SetConstraint ("meme.site_beta_minus", "meme.site_alpha * meme.site_omega_minus", "");

model.generic.AddGlobal (meme.site.bsrel, "meme.site_beta_plus", meme.terms.site_beta_plus);
parameters.DeclareGlobal ("meme.site_beta_plus", {});
parameters.DeclareGlobal ("meme.site_mixture_weight", {});
parameters.SetRange ("meme.site_mixture_weight", terms.range01);

console.log (meme.site.bsrel);

return 0;


/*
fel.report.counts        = {{0,0}};

fel.report.positive_site = {{"" + (1+((fel.filter_specification[fel.report.partition])["coverage"])[fel.report.site]),
                                    fel.report.partition + 1,
                                    Format(fel.report.row[0],10,3),
                                    Format(fel.report.row[1],10,3),
                                    Format(fel.report.row[3],10,3),
                                    "Pos. p = " + Format(fel.report.row[4],6,4)}};

fel.report.negative_site = {{"" + (1+((fel.filter_specification[fel.report.partition])["coverage"])[fel.report.site]),
                                    fel.report.partition + 1,
                                    Format(fel.report.row[0],10,3),
                                    Format(fel.report.row[1],10,3),
                                    Format(fel.report.row[3],10,3),
                                    "Neg. p = " + Format(fel.report.row[4],6,4)}};

*/


meme.site_results = {};

for (meme.partition_index = 0; meme.partition_index < meme.partition_count; meme.partition_index += 1) {
    meme.report.header_done = FALSE;
    meme.table_output_options["header"] = TRUE;

    model.ApplyModelToTree( "fel.site_tree", fel.trees[fel.partition_index], {"default" : fel.site.mg_rev}, None);

    fel.case_respecting_node_names = trees.branch_names (fel.site_tree, TRUE);

    fel.site_patterns = alignments.Extract_site_patterns ((fel.filter_specification[fel.partition_index])["name"]);

    // apply constraints to the site tree
    // alpha = alpha_scaler * branch_length
    // beta  = beta_scaler_test * branch_length or beta_nuisance_test * branch_length

    utility.ForEach (fel.case_respecting_node_names, "_node_",
            '_node_class_ = (fel.selected_branches[fel.partition_index])[_node_ && 1];
             if (_node_class_ == "test") {
                _beta_scaler = fel.scalers[1];
             } else {
                _beta_scaler = fel.scalers[2];
             }
             fel.apply_proportional_site_constraint ("fel.site_tree", _node_, fel.alpha, fel.beta, fel.scalers[0], _beta_scaler, (( fel.final_partitioned_mg_results[terms.json.attribute.branch_length])[fel.partition_index])[_node_]);
        ');

    // create the likelihood function for this site

    ExecuteCommands (alignments.serialize_site_filter
                                       ((fel.filter_specification[fel.partition_index])["name"],
                                       ((fel.site_patterns[0])["sites"])[0],
                     ));

    __make_filter ("fel.site_filter");

    LikelihoodFunction fel.site_likelihood = (fel.site_filter, fel.site_tree);



    estimators.ApplyExistingEstimates ("fel.site_likelihood", fel.site_model_mapping, fel.final_partitioned_mg_results,
                                        "globals only");

    fel.queue = mpi.CreateQueue ({"LikelihoodFunctions": {{"fel.site_likelihood"}},
                                   "Models" : {{"fel.site.mg_rev"}},
                                   "Headers" : {{"libv3/terms-json.bf"}}
                                 });

    /* run the main loop over all unique site pattern combinations */
    utility.ForEachPair (fel.site_patterns, "_pattern_", "_pattern_info_",
        '
            if (_pattern_info_["is_constant"]) {
                fel.store_results (-1,None,{"0" : "fel.site_likelihood",
                                                                "1" : None,
                                                                "2" : fel.partition_index,
                                                                "3" : _pattern_info_,
                                                                "4" : fel.site_model_mapping
                                                                });
            } else {
                mpi.QueueJob (fel.queue, "fel.handle_a_site", {"0" : "fel.site_likelihood",
                                                                "1" : alignments.serialize_site_filter
                                                                   ((fel.filter_specification[fel.partition_index])["name"],
                                                                   (_pattern_info_["sites"])[0]),
                                                                "2" : fel.partition_index,
                                                                "3" : _pattern_info_,
                                                                "4" : fel.site_model_mapping
                                                                },
                                                                "fel.store_results");
            }
        '
    );

    mpi.QueueComplete (fel.queue);
    fel.partition_matrix = {Abs (fel.site_results[fel.partition_index]), Rows (fel.table_headers)};

    utility.ForEachPair (fel.site_results[fel.partition_index], "_key_", "_value_",
    '
        for (fel.index = 0; fel.index < Rows (fel.table_headers); fel.index += 1) {
            fel.partition_matrix [0+_key_][fel.index] = _value_[fel.index];
        }
    '
    );

    fel.site_results[fel.partition_index] = fel.partition_matrix;


}

fel.json [terms.json.MLE ] = {terms.json.headers   : fel.table_headers,
                               terms.json.content : fel.site_results };


io.ReportProgressMessageMD ("fel", "results", "** Found _" + fel.report.counts[0] + "_ sites under positive and _" + fel.report.counts[1] + "_ sites under negative selection at p <= " + fel.pvalue + "**");

selection.io.stopTimer (fel.json [terms.json.timers], "Total time");
selection.io.stopTimer (fel.json [terms.json.timers], "FEL analysis");

io.SpoolJSON (fel.json, fel.codon_data_info["json"]);

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

function fel.apply_proportional_site_constraint (tree_name, node_name, alpha_parameter, beta_parameter, alpha_factor, beta_factor, branch_length) {

    fel.branch_length = (branch_length[terms.synonymous_rate])[terms.MLE];

    node_name = tree_name + "." + node_name;

    ExecuteCommands ("
        `node_name`.`alpha_parameter` := (`alpha_factor`) * fel.branch_length__;
        `node_name`.`beta_parameter`  := (`beta_factor`)  * fel.branch_length__;
    ");
}

//----------------------------------------------------------------------------------------
lfunction fel.handle_a_site (lf, filter_data, partition_index, pattern_info, model_mapping) {

    GetString (lfInfo, ^lf,-1);
    ExecuteCommands (filter_data);
    __make_filter ((lfInfo["Datafilters"])[0]);

    utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);

    ^"fel.alpha_scaler" = 1;
    ^"fel.beta_scaler_test"  = 1;
    ^"fel.beta_scaler_nuisance"  = 1;

    Optimize (results, ^lf);

    alternative = estimators.ExtractMLEs (lf, model_mapping);
    alternative [utility.getGlobalValue("terms.json.log_likelihood")] = results[1][0];

    ^"fel.alpha_scaler" = (^"fel.alpha_scaler" + 3*^"fel.beta_scaler_test")/4;
    parameters.SetConstraint ("fel.beta_scaler_test","fel.alpha_scaler", "");

    Optimize (results, ^lf);

    null = estimators.ExtractMLEs (lf, model_mapping);
    null [utility.getGlobalValue("terms.json.log_likelihood")] = results[1][0];

    /*
    Export (lfs, ^lf);
    fprintf (MESSAGE_LOG, lfs);

    assert (0);
    */

    return {"alternative" : alternative, "null": null};
}

/* echo to screen calls */

function fel.report.echo (fel.report.site, fel.report.partition, fel.report.row) {
    fel.print_row = None;
    if (fel.report.row [4] < fel.pvalue) {
        if (fel.report.row[0] < fel.report.row[1]) {
            fel.print_row = fel.report.positive_site;
            fel.report.counts[0] += 1;
        } else {
            fel.print_row = fel.report.negative_site;
            fel.report.counts [1] += 1;
        }
    }

     if (None != fel.print_row) {
            if (!fel.report.header_done) {
                io.ReportProgressMessageMD("FEL", "" + fel.report.partition, "For partition " + (fel.report.partition+1) + " these sites are significant at p <=" + fel.pvalue + "\n");
                fprintf (stdout,
                    io.FormatTableRow (fel.table_screen_output,fel.table_output_options));
                fel.report.header_done = TRUE;
                fel.table_output_options["header"] = FALSE;
            }

            fprintf (stdout,
                io.FormatTableRow (fel.print_row,fel.table_output_options));
        }

}


lfunction fel.store_results (node, result, arguments) {

    partition_index = arguments [2];
    pattern_info    = arguments [3];

    result_row          = { { 0, // alpha
                          0, // beta
                          0, // alpha==beta
                          0, // LRT
                          1, // p-value,
                          0  // total branch length of tested branches
                      } };



    if (None != result) { // not a constant site
        lrt = math.DoLRT ((result["null"])[utility.getGlobalValue("terms.json.log_likelihood")],
                          (result["alternative"])[utility.getGlobalValue("terms.json.log_likelihood")],
                          1);
        result_row [0] = estimators.GetGlobalMLE (result["alternative"], ^"fel.site_alpha");
        result_row [1] = estimators.GetGlobalMLE (result["alternative"], ^"fel.site_beta");
        result_row [2] = estimators.GetGlobalMLE (result["null"], ^"fel.site_beta");
        result_row [3] = lrt ["LRT"];
        result_row [4] = lrt ["p-value"];

        sum = 0;
        alternative_lengths = ((result["alternative"])[^"terms.json.attribute.branch_length"])[0];

        utility.ForEach (^"fel.case_respecting_node_names", "_node_",
                '_node_class_ = ((^"fel.selected_branches")[`&partition_index`])[_node_ && 1];
                 if (_node_class_ == "test") {
                    `&sum` += ((`&alternative_lengths`)[_node_])[^"terms.json.MLE"];
                 }
            ');

        result_row [5] = sum;
    }


    utility.EnsureKey (^"fel.site_results", partition_index);

    utility.ForEach (pattern_info["sites"], "_value_",
        '
            (fel.site_results[`&partition_index`])[_value_] = `&result_row`;
            fel.report.echo (_value_, `&partition_index`, `&result_row`);
        '
    );


    //assert (0);
}
