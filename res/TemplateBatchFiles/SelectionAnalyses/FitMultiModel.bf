RequireVersion ("2.3.12");


LoadFunctionLibrary("libv3/all-terms.bf"); // must be loaded before CF3x4
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("libv3/models/codon.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/genetic_code.bf");
LoadFunctionLibrary("modules/io_functions.ibf");
LoadFunctionLibrary("modules/selection_lib.ibf");
LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");
LoadFunctionLibrary("libv3/models/codon/MG_REV_MH.bf");
LoadFunctionLibrary("libv3/models/codon/MG_REV_TRIP.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

fitter.analysis_description = {terms.io.info : "Examine whether or not a codon alignment is better fit by models which permit multiple instantaneous substitutions",
                           terms.io.version : "0.1",
                           terms.io.authors : "Sergei L Kosakovsky Pond, Sadie Wisotsky and Alexander Lucacci",
                           terms.io.contact : "spond@temple.edu",
                           terms.io.requirements : "in-frame codon alignment and a phylogenetic tree"
                          };

io.DisplayAnalysisBanner (fitter.analysis_description);

  
namespace fitter.terms {
    MG94 = "Standard MG94";
    MG94x2 = "MG94 with double instantaneous substitutions";
    MG94x3 = "MG94 with double and triple instantaneous substitutions";
    json.site_logl  = "Site Log Likelihood";
    json.evidence_ratios  = "Evidence Ratios";
    json.site_reports = "Site substitutions";
}  
   
fitter.json    = { terms.json.analysis: fitter.analysis_description,
               terms.json.input: {},
               terms.json.fits : {},
               terms.json.timers : {},
               fitter.terms.json.site_logl : {},
               fitter.terms.json.evidence_ratios : {},
               fitter.terms.json.site_reports : {}
              };
 
fitter.display_orders = {terms.original_name: -1,
                         terms.json.nucleotide_gtr: 0,
                         fitter.terms.MG94: 1,
                         fitter.terms.MG94x2: 2,
                         fitter.terms.MG94x3: 3
                        };


selection.io.startTimer (fitter.json [terms.json.timers], "Overall", 0);


namespace fitter {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ({utility.getGlobalValue("terms.prefix"): "fitter", utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "selection.io.SelectAllBranches"}});
}

namespace fitter {
    doGTR ("fitter");
}


//*************** GENERIC FITTER HANDLER ***** 

function fitter.run_model_fit (model_name, model_generator, initial_values) {
    io.ReportProgressMessageMD ("fitter", model_name, "Fitting `model_name`");
    selection.io.startTimer (fitter.json [terms.json.timers], model_name, fitter.display_order [model_name]);

    fitter.results =  estimators.FitCodonModel (fitter.filter_names, fitter.trees, model_generator, fitter.codon_data_info [utility.getGlobalValue("terms.code")],
     {
        utility.getGlobalValue("terms.run_options.model_type"): utility.getGlobalValue("terms.global"),
        utility.getGlobalValue("terms.run_options.retain_lf_object"): TRUE
    }, fitter.gtr_results);


    ConstructCategoryMatrix (fitter.run_model_fit.sl, ^(fitter.results[terms.likelihood_function]), SITE_LOG_LIKELIHOODS);
    (fitter.json [fitter.terms.json.site_logl])[model_name] = fitter.run_model_fit.sl;
 
    io.ReportProgressMessageMD("fitter", model_name, "* " + selection.io.report_fit (fitter.results, 0, fitter.codon_data_info[terms.data.sample_size]));
    fitter.global_dnds = selection.io.extract_global_MLE_re (fitter.results, "^(" + terms.parameters.omega_ratio + "|" + terms.parameters.multiple_hit_rate  + "|" + terms.parameters.triple_hit_rate + ")");

    utility.ForEach (fitter.global_dnds, "_value_", 'io.ReportProgressMessageMD ("fitter", model_name, "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));');

    selection.io.json_store_lf (fitter.json,
                                model_name,
                                fitter.results[terms.fit.log_likelihood],
                                fitter.results[terms.parameters],
                                fitter.sample_size,
                                utility.Map (fitter.results[terms.global], "_value_", '_value_ [terms.fit.MLE]'),
                                fitter.display_orders[model_name]);


    utility.ForEachPair (fitter.filter_specification, "_key_", "_value_",
        'selection.io.json_store_branch_attribute(fitter.json, model_name, terms.branch_length, fitter.display_orders[model_name],
                                                 _key_,
                                                 selection.io.extract_branch_info((fitter.results[terms.branch_length])[_key_], "selection.io.branch.length"));');

    selection.io.stopTimer (fitter.json [terms.json.timers], model_name);
    return fitter.results;
}

//*************** END GENERIC FITTER HANDLER ***** 

fitter.one_hit_results = fitter.run_model_fit (fitter.terms.MG94, "models.codon.MG_REV.ModelDescription", fitter.gtr_results);

utility.Extend (fitter.one_hit_results[terms.global],
                {
                    terms.parameters.multiple_hit_rate : { utility.getGlobalValue ("terms.fit.MLE") : 0.05, terms.fix : FALSE}

                });

fitter.two_hit_results = fitter.run_model_fit (fitter.terms.MG94x2, "models.codon.MG_REV_MH.ModelDescription", fitter.one_hit_results);

utility.Extend (fitter.two_hit_results[terms.global],
                {
                    terms.parameters.triple_hit_rate : { utility.getGlobalValue ("terms.fit.MLE") : 0.05, terms.fix : FALSE}

                });

fitter.three_hit_results = fitter.run_model_fit (fitter.terms.MG94x3, "models.codon.MG_REV_TRIP.ModelDescription", fitter.one_hit_results);

selection.io.stopTimer (fitter.json [terms.json.timers], "Overall");

// PERFORM HYPOTHESIS TESTING AND POPULATE TABLE REPORT

fitter.LRT = {
    "Double-hit vs single-hit" : math.DoLRT (fitter.one_hit_results[terms.fit.log_likelihood], fitter.two_hit_results[terms.fit.log_likelihood], 1),
    "Triple-hit vs single-hit" : math.DoLRT (fitter.one_hit_results[terms.fit.log_likelihood], fitter.three_hit_results[terms.fit.log_likelihood], 2),
    "Triple-hit vs double-hit" : math.DoLRT (fitter.two_hit_results[terms.fit.log_likelihood], fitter.three_hit_results[terms.fit.log_likelihood], 1)
};

fitter.json [terms.json.test_results] = fitter.LRT;

fitter.table_output_options = {
        utility.getGlobalValue("terms.table_options.header"): TRUE,
        utility.getGlobalValue("terms.table_options.minimum_column_width"): 10,
        utility.getGlobalValue("terms.table_options.align"): "center",
        utility.getGlobalValue("terms.table_options.column_widths"): {
            "0": 35,
            "1": 12,
            "2": 12,
            "3": 12,
            "4": 33,
            "5": 12,
            "6": 30
        }
    };
    
fitter.report = {
    {
        "Model", "Log-L", "omega", " 2-hit rate", "p-value", "3-hit rate", "p-value"
    }
};

io.ReportProgressMessageMD ("fitter", "testing", "Summary of rate estimates and significance testing");

fprintf(stdout, io.FormatTableRow(fitter.report , fitter.table_output_options));
fitter.table_output_options[utility.getGlobalValue("terms.table_options.header")] = FALSE;

// single rate model
fprintf(stdout, io.FormatTableRow(
    {
     {
     fitter.terms.MG94, 
     Format (fitter.one_hit_results[terms.fit.log_likelihood], 8, 2), 
     Format (((fitter.one_hit_results[terms.global])[terms.parameters.omega_ratio])[terms.fit.MLE],8,4),
     "N/A",
     "N/A",
     "N/A",
     "N/A"
     }
    },
    fitter.table_output_options)
);

// double-hit rate model
fprintf(stdout, io.FormatTableRow(
    {
     {
     fitter.terms.MG94 + " + 2 hits", 
     Format (fitter.two_hit_results[terms.fit.log_likelihood], 8, 2), 
     Format (((fitter.two_hit_results[terms.global])[terms.parameters.omega_ratio])[terms.fit.MLE],8,4),
     Format (((fitter.two_hit_results[terms.global])[terms.parameters.multiple_hit_rate])[terms.fit.MLE],8,4) ,
     Format ((fitter.LRT["Double-hit vs single-hit"])[terms.p_value], 8, 4) + " (2-hit rate = 0)",
     "N/A",
     "N/A"
     }
    },
    fitter.table_output_options)
);

// triple-hit rate model
fprintf(stdout, io.FormatTableRow(
    {
     {
     fitter.terms.MG94 + " + 2 or 3 hits", 
     Format (fitter.three_hit_results[terms.fit.log_likelihood], 8, 2), 
     Format (((fitter.three_hit_results[terms.global])[terms.parameters.omega_ratio])[terms.fit.MLE],8,4),
     Format (((fitter.three_hit_results[terms.global])[terms.parameters.multiple_hit_rate])[terms.fit.MLE],8,4),
     Format ((fitter.LRT["Triple-hit vs single-hit"])[terms.p_value], 8, 4) + " (2&3-hit rates = 0)",
     Format (((fitter.three_hit_results[terms.global])[terms.parameters.triple_hit_rate])[terms.fit.MLE],8,4),
     Format ((fitter.LRT["Triple-hit vs double-hit"])[terms.p_value], 8, 4) + " (3-hit rate = 0)"
     }
    },
    fitter.table_output_options)
);



(fitter.json [fitter.terms.json.evidence_ratios])["Two-hit"] = 
        fitter.EvidenceRatios ((fitter.json[fitter.terms.json.site_logl])[fitter.terms.MG94x2],
                               (fitter.json[fitter.terms.json.site_logl])[fitter.terms.MG94]);

(fitter.json [fitter.terms.json.evidence_ratios])["Three-hit"] = 
        fitter.EvidenceRatios ((fitter.json[fitter.terms.json.site_logl])[fitter.terms.MG94x3],
                               (fitter.json[fitter.terms.json.site_logl])[fitter.terms.MG94x2]);
                    
fitter.callout = {};

utility.ForEachPair ((fitter.json [fitter.terms.json.evidence_ratios])["Two-hit"], "_index_", "_value_", 
'
    if (Log (_value_) > 2) {
        fitter.callout [_index_[1]] = 1;
    }
');

utility.ForEachPair ((fitter.json [fitter.terms.json.evidence_ratios])["Three-hit"], "_index_", "_value_", 
'
    if (Log (_value_) > 2) {
        fitter.callout [_index_[1]] = 1;
    }
');

if (utility.Array1D (fitter.callout)) {

    // reconstruct ancestors here
    
    fitter.site_reports = {};
    
    fitter.ancestral_cache = ancestral.build (fitter.three_hit_results[terms.likelihood_function], 0, None);
    
    io.ReportProgressMessageMD ("fitter", "sites", "" + utility.Array1D (fitter.callout) + " individual " + io.SingularOrPlural (utility.Array1D (fitter.callout) , "site", "sites") + " which showed sufficiently strong preference for multiple-hit models");
    fitter.table_output_options = {
            utility.getGlobalValue("terms.table_options.header"): TRUE,
            utility.getGlobalValue("terms.table_options.minimum_column_width"): 10,
            utility.getGlobalValue("terms.table_options.align"): "center",
            utility.getGlobalValue("terms.table_options.column_widths"): {
                "0": 10,
                "1": 25,
                "2": 25,
                "3": 60
            }
        };
    
    fitter.report = {
        {
            "Site", "Evidence Ratio (2-hit)", "Evidence Ratio (3-hit)", "Substitutions"
        }
    };
    fprintf(stdout, io.FormatTableRow(fitter.report , fitter.table_output_options));
    
    fitter.table_output_options[utility.getGlobalValue("terms.table_options.header")] = FALSE;
    utility.ForEachPair (fitter.callout, "_site_", "_value_", "
    
         fitter.site_reports [_site_] = (ancestral.ComputeSubstitutionBySite (fitter.ancestral_cache, +_site_, None))[terms.substitutions];
         
         fprintf(stdout, io.FormatTableRow(
            {
                {
                    '' + (+_site_ + 1),
                    Format (((fitter.json [fitter.terms.json.evidence_ratios])['Two-hit'])[+_site_], 10, 4),
                    Format (((fitter.json [fitter.terms.json.evidence_ratios])['Three-hit'])[+_site_], 10, 4),
                    fitter.SubstitutionHistory (fitter.site_reports [_site_])
                }
            }
        , fitter.table_output_options));
    ");
    
    fitter.json [fitter.terms.json.site_reports] = fitter.site_reports;


} else {
    io.ReportProgressMessageMD ("fitter", "sites", "No individual sites showed sufficiently strong preference for multiple-hit models");

}

io.ReportProgressMessageMD ("fitter", "writing", "Writing detailed analysis report to \`" + fitter.codon_data_info [terms.json.json] + "\'");

io.SpoolJSON (fitter.json, fitter.codon_data_info [terms.json.json]);
return fitter.json;

//------------------------------------------------------------------------------
lfunction fitter.EvidenceRatios (ha, h0) {
    return ha["Exp(_MATRIX_ELEMENT_VALUE_-h0[_MATRIX_ELEMENT_COLUMN_])"];
}


//------------------------------------------------------------------------------

lfunction fitter.SubstitutionHistory (subs) {
    result = "";
    result * 128;
    keys = utility.sortStrings (utility.Keys (subs));

    for (i = 0; i < Abs (subs); i+=1) {
        source  = keys[i];
        targets  = subs[source];
        if (i > 0) {
            result * ", ";
        }
        result * (source + "->");
        keys2 =  utility.sortStrings (utility.Keys (targets));
        for (k = 0; k < Abs (targets); k+=1) {
           result * (keys2[k] + "(" + Abs(targets[keys2[k]]) + ")");
        }
    }

    result * 0;
    return result;

}