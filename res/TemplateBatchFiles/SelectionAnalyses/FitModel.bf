RequireVersion ("2.5.34");


LoadFunctionLibrary("libv3/all-terms.bf"); 
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");

LoadFunctionLibrary("modules/io_functions.ibf");
LoadFunctionLibrary("modules/selection_lib.ibf");

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

fitter.analysis_description = {terms.io.info : "Fit one of the available codon models to the data, report fit results and model estimates. Model specific options will differ. You can also perform Likelihood Ratio Tests (LRT) on any model parameter by passing --parameter_name value on the command line.",
                               terms.io.version : "0.1",
                               terms.io.authors : "Sergei L Kosakovsky Pond",
                               terms.io.contact : "spond@temple.edu",
                               terms.io.requirements : "in-frame codon alignment and a phylogenetic tree"
                              };

io.DisplayAnalysisBanner (fitter.analysis_description);

  
KeywordArgument ("code",        "Which genetic code should be used", "Universal");  
KeywordArgument ("alignment",   "An in-frame codon alignment in one of the formats supported by HyPhy");
KeywordArgument ("tree",        "A phylogenetic tree", null, "Please select a tree file for the data:");
KeywordArgument ("model",       "Use the following substitution model");  

fitter.json    = { 
                   terms.json.analysis: fitter.analysis_description,
                   terms.json.input: {},
                   terms.json.fits : {},
                   terms.json.timers : {},
                 };
               
fitter.terms.codon_model = "Codon model";
 
fitter.display_orders = {terms.original_name            :  -1,
                         terms.json.nucleotide_gtr      :  0,
                         fitter.terms.codon_model       :  1
                        };


selection.io.startTimer (fitter.json [terms.json.timers], "Overall", 0);

namespace fitter {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ({utility.getGlobalValue("terms.prefix"): "fitter", 
                utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "selection.io.SelectAllBranches"}});
}


fitter.model_name = io.PromptUserForString ("Model to fit");
LoadFunctionLibrary ("libv3/models/codon/`fitter.model_name`.bf");

fitter.model_generator = "model.codon.`fitter.model_name`.`terms.model.argument_collector`";

io.CheckAssertion("utility.HasUserFunction(fitter.model_generator,2)", "This model does not provide the required `terms.model.argument_collector` function");

KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'MG94.json')", fitter.codon_data_info [terms.json.json]);
fitter.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");

namespace fitter {
    doGTR ("fitter");
}


io.ReportProgressMessageMD ("fitter", fitter.model_name,  "Fitting `fitter.model_name`");
selection.io.startTimer (fitter.json [terms.json.timers], fitter.model_name , fitter.display_order [fitter.terms.codon_model ]);

fitter.model_type = terms.global;

fitter.results =  estimators.FitCodonModel (fitter.filter_names, fitter.trees, fitter.model_generator, fitter.codon_data_info [utility.getGlobalValue("terms.code")],
    {
        terms.run_options.model_type: fitter.model_type,
        terms.run_options.retain_lf_object: TRUE,
        terms.run_options.retain_model_object : TRUE
    }, 
    fitter.gtr_results);



io.ReportProgressMessageMD("fitter", fitter.model_name , "* " + selection.io.report_fit (fitter.results, 0, fitter.codon_data_info[terms.data.sample_size]));
io.ReportProgressMessageMD("fitter", fitter.model_name , "* " + "Alignment-wide parameters and their maximum likelihood estimates (MLEs)");
fitter.globals = selection.io.extract_global_MLE_re (fitter.results, '.');

fitter.table_output_options = {
        utility.getGlobalValue("terms.table_options.header"): 1,
        utility.getGlobalValue("terms.table_options.minimum_column_width"): 16,
        utility.getGlobalValue("terms.table_options.align"): "left",
        utility.getGlobalValue("terms.table_options.column_widths"): {
            "0": 120,
            "1" :20
        }
    };
    
fitter.header = {
        2,
        1
    };
    
fitter.header[0] = "Parameter";
fitter.header[1] = "MLE";

fprintf(stdout, io.FormatTableRow(fitter.header, fitter.table_output_options));
fitter.table_output_options[utility.getGlobalValue("terms.table_options.header")] = FALSE;

fitter.print_row = {
            2,
            1
        };
fitter.print_row[0] := fitter.mle[terms.description];
fitter.print_row[1] := Format (fitter.mle[terms.fit.MLE], 16, 6);
for (fitter.mle; in; fitter.globals) {
        fprintf(stdout, io.FormatTableRow(fitter.print_row, fitter.table_output_options));
}

selection.io.json_store_lf (fitter.json,
                            fitter.terms.codon_model ,
                            fitter.results[terms.fit.log_likelihood],
                            fitter.results[terms.parameters],
                            fitter.sample_size,
                            utility.Map (fitter.results[terms.global], "_value_", '_value_ [terms.fit.MLE]'), fitter.display_orders[fitter.terms.codon_model ]);
                  
KeywordArgument ("save-fit", "Save MG94 model fit to this file (default is not to save)", "/dev/null");
io.SpoolLFToPath(fitter.results[terms.likelihood_function], io.PromptUserForFilePath ("Save MG94 model fit to this file ['/dev/null' to skip]"));


// LRT testing

fitter.lrt_testing = {};
fitter.pvalues = {};

function fitter.SetToOne (set) {
    if (set) {
        fitter.SetToOne.stash = Eval (fitter.param.ID);
        parameters.SetConstraint (fitter.param.ID, fitter.use_this_value, "");
    } else {
        parameters.SetValue (fitter.param.ID, fitter.SetToOne.stash);
    }
    return 1;
}

for (fitter.var; in; fitter.globals) {
    fitter.param.ID  = fitter.var[terms.id];
    fitter.param.MLE = fitter.var[terms.fit.MLE];
    fitter.param.desc = fitter.var[terms.description];
    
    ExecuteCommands ('
        KeywordArgument  ("`fitter.param.desc`", "Point hypothesis LRT (set null value, or null to skip) ", "null");
        fitter.use_this_value = io.PromptUserForString ("Point hypothesis LRT (set parameter value for the null, or enter \\\" null\\\" to skip)");
        
    ');
    
    if (fitter.use_this_value != "null") {
        fitter.lrt_testing [fitter.param.desc] = (estimators.ConstrainAndRunLRT (fitter.results[terms.likelihood_function], "fitter.SetToOne"));
        (fitter.lrt_testing [fitter.param.desc])[terms.fit.MLE] =  +fitter.use_this_value;
        io.ReportProgressMessageMD("fitter", "LRT", "\nLikelihood ratio test for _`fitter.param.desc` == `fitter.use_this_value`_, **p = " + Format ((fitter.lrt_testing[fitter.param.desc])[terms.p_value], 8, 4) + "**.");
        fitter.pvalues  [fitter.param.desc] =  (fitter.lrt_testing     [fitter.param.desc])[terms.p_value];

    }
}

fitter.json [terms.json.test_results] = fitter.lrt_testing;

selection.io.stopTimer (fitter.json [terms.json.timers], fitter.model_name);

selection.io.stopTimer (fitter.json [terms.json.timers], "Overall");
io.ReportProgressMessageMD ("fitter", "writing", "Writing detailed analysis report to \`" + fitter.codon_data_info [terms.json.json] + "\'");
io.SpoolJSON (fitter.json, fitter.codon_data_info [terms.json.json]);

return fitter.results;
