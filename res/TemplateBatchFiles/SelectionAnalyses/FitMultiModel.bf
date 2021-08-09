RequireVersion ("2.4.0");


LoadFunctionLibrary("libv3/all-terms.bf"); // must be loaded before CF3x4
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("libv3/models/codon.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/genetic_code.bf");
LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary("SelectionAnalyses/modules/selection_lib.ibf");
LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");
LoadFunctionLibrary("libv3/models/codon/MG_REV_MH.bf");
LoadFunctionLibrary("libv3/models/codon/MG_REV_TRIP.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");
LoadFunctionLibrary("libv3/models/rate_variation.bf");
LoadFunctionLibrary("libv3/models/codon/BS_REL.bf");


utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);
//utility.SetEnvVariable ("VERBOSITY_LEVEL", 10);
fitter.rate_classes = 3;

KeywordArgument ("code",      "Which genetic code should be used", "Universal");
    /**
        keyword, description (for inline documentation and help messages), default value
    */
KeywordArgument ("alignment", "An in-frame codon alignment in one of the formats supported by HyPhy");
    /**
        keyword, description (for inline documentation and help messages), no default value,
        meaning that it will be required
    */

KeywordArgument ("tree",      "A phylogenetic tree (optionally annotated with {})", null, "Please select a tree file for the data:");
    /** the use of null as the default argument means that the default expectation is for the
        argument to be missing, i.e. the tree is expected to be in the file
        the fourth, optional argument, can match this keyword with the dialog prompt / choice list title,
        meaning that it can only be consumed when this dialog prompt / choice list is invoked
        This allows handling some branching logic conditionals
    */

KeywordArgument ("rates", "The number omega rate classes to include in the model [1-10, default 3, 1 to turn-off rate variation]", 3);

KeywordArgument ("triple-islands", "Use a separate rate parameter for synonymous triple-hit substitutions", "No");


fitter.analysis_description = {terms.io.info : "Examine whether or not a codon alignment is better fit by models which permit multiple instantaneous substitutions. v0.2 adds a separate rate for codon-island triple-hit rates",
                           terms.io.version : "1.0",
                           terms.io.authors : "Sergei L Kosakovsky Pond, Sadie Wisotsky and Alexander Lucaci",
                           terms.io.contact : "spond@temple.edu",
                           terms.io.reference: "Submitted; preprint at hyphy.org/resources/fmm.pdf",
                           terms.io.requirements : "in-frame codon alignment and a phylogenetic tree"
                          };

io.DisplayAnalysisBanner (fitter.analysis_description);


namespace fitter.terms {
    MG94 = "Standard MG94";
    MG94x2 = "MG94 with double instantaneous substitutions";
    MG94x3 = "MG94 with double and triple instantaneous substitutions";
    MG94x3xS = "MG94 with double and triple instantaneous substitutions [only synonymous islands]";
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
    LoadFunctionLibrary ("SelectionAnalyses/modules/shared-load-file.bf");
    load_file ({utility.getGlobalValue("terms.prefix"): "fitter", utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "selection.io.SelectAllBranches"}});
}

fitter.rate_classes = io.PromptUser ("The number of omega rate classes to include in the model", 3, 1, 10, TRUE);

fitter.do_islands = io.SelectAnOption ({"Yes" : "Use a separate rate parameter for synonymous triple-hit substitutions (e.g. serine islands)",
                                    "No"  : "All triple hits have the same rate multiplier"},
                                    "Synonymous triple-hit substitutions use a separate rate"
                                    ) == "Yes";


namespace fitter {
    doGTR ("fitter");
}

KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'FITTER.json')", fitter.codon_data_info [terms.json.json]);
fitter.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");

KeywordArgument ("save-fit", "Write model fit files (HyPhy NEXUS) to this file path with extensions .MODEL_NAME.bf; default is NOT to save, or 'dev/null'", "/dev/null");
fitter.save_model_path = io.PromptUserForFilePath ("Save model fit files to");


estimators.fixSubsetOfEstimates(fitter.gtr_results, fitter.gtr_results[terms.global]);

namespace fitter {
    scaler_prefix = "fitter.scaler";
    doPartitionedMG ("fitter", FALSE);
}

//Export(lf_serialized, ^(fitter.partioned_mg_results));
//fprintf(PROMPT_FOR_FILE,CLEAR_FILE, lf_serialized, "\n");

//*************** GENERIC FITTER HANDLER *****

function fitter.run_model_fit (model_name, model_generator, initial_values) {
    io.ReportProgressMessageMD ("fitter", model_name, "Fitting `model_name`");
    selection.io.startTimer (fitter.json [terms.json.timers], model_name, fitter.display_order [model_name]);

    fitter.results =  estimators.FitCodonModel (fitter.filter_names, fitter.trees, model_generator, fitter.codon_data_info [utility.getGlobalValue("terms.code")],
     {
        utility.getGlobalValue  ("terms.run_options.model_type"): utility.getGlobalValue("terms.global"),
        utility.getGlobalValue  ("terms.run_options.retain_lf_object"): TRUE,
        utility.getGlobalValue  ("terms.run_options.retain_model_object"): TRUE
    }, initial_values);


    if (fitter.save_model_path != "/dev/null") {
        //^(fitter.results[terms.likelihood_function])
        io.SpoolLF(fitter.results[terms.likelihood_function], fitter.save_model_path, model_name);
    }

    fitter.models = fitter.results [terms.model];
   
    if (^"fitter.rate_classes" > 1) {
        fitter.cat_info = rate_variation.extract_category_information (fitter.models[fitter.models ["INDEXORDER"][0]]);
        fitter.cat_info = Transpose (fitter.cat_info[fitter.cat_info["INDEXORDER"][0]])%0;
    } else {
        fitter.cat_info = {{1,1}};
    }
    

    ConstructCategoryMatrix (fitter.run_model_fit.sl, ^(fitter.results[terms.likelihood_function]), SITE_LOG_LIKELIHOODS);
    (fitter.json [fitter.terms.json.site_logl])[model_name] = fitter.run_model_fit.sl;

    io.ReportProgressMessageMD("fitter", model_name, "* " + selection.io.report_fit (fitter.results, 0, fitter.codon_data_info[terms.data.sample_size]));
    fitter.global_dnds = selection.io.extract_global_MLE_re (fitter.results, "^(" + terms.parameters.omega_ratio + "|" + terms.parameters.multiple_hit_rate  + "|" + terms.parameters.triple_hit_rate + ")");


    utility.ForEach (fitter.global_dnds, "_value_", 'io.ReportProgressMessageMD ("fitter", model_name, "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));');
   
    if (^"fitter.rate_classes" > 1) {
        io.ReportProgressMessageMD("fitter", model_name, "* The following relative rate distribution (mean 1) for site-to-site **non-synonymous** rate variation was inferred");
        selection.io.report_distribution (fitter.cat_info);
    }
    
    fitter.cat_info = {
        terms.rate_variation.distribution : fitter.cat_info,
        terms.parameters: utility.Map (fitter.results[terms.global], "_value_", '_value_ [terms.fit.MLE]')
    };

    selection.io.json_store_lf (fitter.json,
                                model_name,
                                fitter.results[terms.fit.log_likelihood],
                                fitter.results[terms.parameters],
                                fitter.sample_size,
                                fitter.cat_info,
                                fitter.display_orders[model_name]);


    utility.ForEachPair (fitter.filter_specification, "_key_", "_value_",
        'selection.io.json_store_branch_attribute(fitter.json, model_name, terms.branch_length, fitter.display_orders[model_name],
                                                 _key_,
                                                 selection.io.extract_branch_info((fitter.results[terms.branch_length])[_key_], "selection.io.branch.length"));');

    selection.io.stopTimer (fitter.json [terms.json.timers], model_name);
    return fitter.results;
}

//*************** END GENERIC FITTER HANDLER *****

//****** ADDING RATE VARIATION ****

function fitter.modifier_omega (q_ij, from, to, namespace, cat_name) {
	if (Abs (q_ij[terms.model.rate_entry])) {
	    if (utility.Has (q_ij[terms.global], terms.parameters.omega_ratio, "String")) {
		    q_ij[terms.model.rate_entry] = "(" + q_ij[terms.model.rate_entry] + ")*" + cat_name;
		    //console.log (q_ij);
		}
	}
	return q_ij;
}


lfunction MG_REV.model.with.GDD (type, code) {
        def = models.codon.MG_REV.ModelDescription (type, code);
        if (^"fitter.rate_classes" > 1) {
            def [utility.getGlobalValue("terms.model.rate_variation")] = rate_variation.types.GDD.factory ({utility.getGlobalValue("terms.rate_variation.bins") : utility.getGlobalValue("fitter.rate_classes")});
            (def [utility.getGlobalValue("terms.model.rate_variation")])[utility.getGlobalValue("terms.rate_variation.rate_modifier")] = "fitter.modifier_omega";
        }
        return def;
    }

lfunction MG_REV_MH.model.with.GDD (type, code) {
        def = models.codon.MG_REV_MH.ModelDescription (type, code);
        if (^"fitter.rate_classes" > 1) {
            def [utility.getGlobalValue("terms.model.rate_variation")] = rate_variation.types.GDD.factory ({utility.getGlobalValue("terms.rate_variation.bins") : utility.getGlobalValue("fitter.rate_classes")});
            (def [utility.getGlobalValue("terms.model.rate_variation")])[utility.getGlobalValue("terms.rate_variation.rate_modifier")] = "fitter.modifier_omega";
        }
        return def;
    }

lfunction MG_REV_TRIP.model.with.GDD (type, code) {
        def = models.codon.MG_REV_TRIP.ModelDescription (type, code);
        if (^"fitter.rate_classes" > 1) {
            def [utility.getGlobalValue("terms.model.rate_variation")] = rate_variation.types.GDD.factory ({utility.getGlobalValue("terms.rate_variation.bins") : utility.getGlobalValue("fitter.rate_classes")});
            (def [utility.getGlobalValue("terms.model.rate_variation")])[utility.getGlobalValue("terms.rate_variation.rate_modifier")] = "fitter.modifier_omega";
        }
        return def;
    }

//*************** END RATE VARIATION *****


fitter.one_hit_results = fitter.run_model_fit (fitter.terms.MG94, "MG_REV.model.with.GDD", fitter.partitioned_mg_results);

utility.Extend (fitter.one_hit_results[terms.global],
                {
                    terms.parameters.multiple_hit_rate : { utility.getGlobalValue ("terms.fit.MLE") : 0.05, terms.fix : FALSE}

                });

fitter.two_hit_results = fitter.run_model_fit (fitter.terms.MG94x2, "MG_REV_MH.model.with.GDD", fitter.one_hit_results);

utility.Extend (fitter.two_hit_results[terms.global],
                {
                    terms.parameters.triple_hit_rate : { utility.getGlobalValue ("terms.fit.MLE") : 0.05, terms.fix : FALSE},
                    terms.parameters.triple_hit_rate_syn : { utility.getGlobalValue ("terms.fit.MLE") : 0.05, terms.fix : FALSE}

                });


lfunction fitter.MG_REV_TRIP.model.with.GDD (type, code) {
    model = MG_REV_TRIP.model.with.GDD (type, code);
    model[utility.getGlobalValue("terms.model.post_definition")] = "fitter.handle_triple_islands";
    return model;
}

lfunction fitter.MG_REV_TRIP.model.with.GDD.islands (type, code) {
    model = MG_REV_TRIP.model.with.GDD (type, code);
    model[utility.getGlobalValue("terms.model.post_definition")] = "fitter.handle_only_triple_islands";
    return model;
}

lfunction fitter.handle_triple_islands (model) {
    if (utility.getGlobalValue("fitter.do_islands") == FALSE) {
        parameters.SetConstraint (model.generic.GetGlobalParameter (model, utility.getGlobalValue("terms.parameters.triple_hit_rate_syn")),
                                  model.generic.GetGlobalParameter (model, utility.getGlobalValue("terms.parameters.triple_hit_rate")), "");
    }
    return models.generic.post.definition (model);
}

lfunction fitter.handle_only_triple_islands (model) {
    parameters.SetConstraint (model.generic.GetGlobalParameter (model, utility.getGlobalValue("terms.parameters.triple_hit_rate")), "0", "");
    return models.generic.post.definition (model);
}

fitter.three_hit_results = fitter.run_model_fit (fitter.terms.MG94x3, "fitter.MG_REV_TRIP.model.with.GDD", fitter.one_hit_results);

if (fitter.do_islands) {
     fitter.three_hit_island_results = fitter.run_model_fit (fitter.terms.MG94x3xS, "fitter.MG_REV_TRIP.model.with.GDD.islands", fitter.three_hit_results);
}
//

selection.io.stopTimer (fitter.json [terms.json.timers], "Overall");

// PERFORM HYPOTHESIS TESTING AND POPULATE TABLE REPORT

fitter.LRT = {
    "Double-hit vs single-hit" : math.DoLRT (fitter.one_hit_results[terms.fit.log_likelihood], fitter.two_hit_results[terms.fit.log_likelihood], 1),
    "Triple-hit vs single-hit" : math.DoLRT (fitter.one_hit_results[terms.fit.log_likelihood], fitter.three_hit_results[terms.fit.log_likelihood], 2 + fitter.do_islands),
    "Triple-hit vs double-hit" : math.DoLRT (fitter.two_hit_results[terms.fit.log_likelihood], fitter.three_hit_results[terms.fit.log_likelihood], 1 + fitter.do_islands)
};

if (fitter.do_islands) {
    fitter.LRT ["Triple-hit-island vs double-hit"] = math.DoLRT (fitter.two_hit_results[terms.fit.log_likelihood], fitter.three_hit_island_results[terms.fit.log_likelihood],1);
    fitter.LRT ["Triple-hit vs Triple-hit-island"] = math.DoLRT (fitter.three_hit_island_results[terms.fit.log_likelihood], fitter.three_hit_results[terms.fit.log_likelihood],1);
}

fitter.json [terms.json.test_results] = fitter.LRT;

fitter.table_output_options = {
        utility.getGlobalValue("terms.table_options.header"): TRUE,
        utility.getGlobalValue("terms.table_options.minimum_column_width"): 10,
        utility.getGlobalValue("terms.table_options.align"): "center",
        utility.getGlobalValue("terms.table_options.column_widths"): {
            "0": 38,
            "1": 12,
            "2": 12,
            "3": 12,
            "4": 36,
            "5": 12,
            "6": 36
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

if (fitter.do_islands) {



    fprintf(stdout, io.FormatTableRow(
        {
         {
         fitter.terms.MG94 + " + 3 hits (islands)",
         Format (fitter.three_hit_island_results[terms.fit.log_likelihood], 8, 2),
         Format (((fitter.three_hit_island_results[terms.global])[terms.parameters.omega_ratio])[terms.fit.MLE],8,4),
         Format (((fitter.three_hit_island_results[terms.global])[terms.parameters.triple_hit_rate_syn])[terms.fit.MLE],8,4) ,
         Format ((fitter.LRT["Triple-hit-island vs double-hit"])[terms.p_value], 8, 4) + " (3-hit island vs 2-hit)",
         Format (((fitter.three_hit_results[terms.global])[terms.parameters.triple_hit_rate])[terms.fit.MLE],8,4),
         Format ((fitter.LRT["Triple-hit vs Triple-hit-island"])[terms.p_value], 8, 4) + " (3-hit = 0)"
         }
        },
        fitter.table_output_options)
    );
}

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
     Format ((fitter.LRT["Triple-hit vs double-hit"])[terms.p_value], 8, 4) + " (3-hit rate(s) = 0)"
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

if (fitter.do_islands) {
    (fitter.json [fitter.terms.json.evidence_ratios])["Three-hit islands vs 2-hit"] =
        fitter.EvidenceRatios ((fitter.json[fitter.terms.json.site_logl])[fitter.terms.MG94x3xS],
                               (fitter.json[fitter.terms.json.site_logl])[fitter.terms.MG94x2]);

    (fitter.json [fitter.terms.json.evidence_ratios])["Three-hit vs three-hit islands"] =
        fitter.EvidenceRatios ((fitter.json[fitter.terms.json.site_logl])[fitter.terms.MG94x3],
                               (fitter.json[fitter.terms.json.site_logl])[fitter.terms.MG94x3xS]);

}


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

if (fitter.do_islands) {
    utility.ForEachPair ((fitter.json [fitter.terms.json.evidence_ratios])["Three-hit islands vs 2-hit"], "_index_", "_value_",
    '
        if (Log (_value_) > 2) {
            fitter.callout [_index_[1]] = 1;
        }
    ');
    utility.ForEachPair ((fitter.json [fitter.terms.json.evidence_ratios])["Three-hit vs three-hit islands"], "_index_", "_value_",
    '
        if (Log (_value_) > 2) {
            fitter.callout [_index_[1]] = 1;
        }
    ');
}

if (utility.Array1D (fitter.callout)) {
    // reconstruct ancestors here
    
    fitter.site2partition = {};
    fitter.ancestral_caches = {};
    for (index, partition; in; fitter.filter_specification) {
        for (site_index, site; in; partition[terms.data.coverage]) {
            fitter.site2partition [site] = {{index__,site_index__}};
        }
        fitter.ancestral_caches [index] = ancestral.build (fitter.three_hit_results[terms.likelihood_function], +index, None);
    }
    fitter.site_reports = {};
    

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


    if (fitter.do_islands) {
        (fitter.table_output_options[utility.getGlobalValue("terms.table_options.column_widths")]) [3] = 25;
        (fitter.table_output_options[utility.getGlobalValue("terms.table_options.column_widths")]) [4] = 25;
        (fitter.table_output_options[utility.getGlobalValue("terms.table_options.column_widths")]) [5] = 60;
         fitter.report = {
            {
                "Site", "ER (2 vs 1)", "ER (3 vs 2)", "ER (3-island vs 2)", "ER (3-island vs 3)", "Substitutions"
            }
        };
    } else {
        fitter.report = {
            {
                "Site", "ER (2 vs 1)", "ER (3 vs 2)", "Substitutions"
            }
        };
    }
    fprintf(stdout, "\n", io.FormatTableRow(fitter.report , fitter.table_output_options));

    fitter.table_output_options[utility.getGlobalValue("terms.table_options.header")] = FALSE;
    utility.ForEachPair (fitter.callout, "_site_", "_value_", "
    
         _site_map_ = fitter.site2partition [_site_];

         fitter.site_reports [_site_] = (ancestral.ComputeSubstitutionBySite ((fitter.ancestral_caches[_site_map_[0]]), +(_site_map_[1]), None))[terms.substitutions];

             if (fitter.do_islands) {
                 fprintf(stdout, io.FormatTableRow(
                    {
                        {
                            '' + (+_site_ + 1),
                            Format (((fitter.json [fitter.terms.json.evidence_ratios])['Two-hit'])[+_site_], 10, 4),
                            Format (((fitter.json [fitter.terms.json.evidence_ratios])['Three-hit'])[+_site_], 10, 4),
                            Format (((fitter.json [fitter.terms.json.evidence_ratios])['Three-hit islands vs 2-hit'])[+_site_], 10, 4),
                            Format (((fitter.json [fitter.terms.json.evidence_ratios])['Three-hit vs three-hit islands'])[+_site_], 10, 4),
                            fitter.SubstitutionHistory (fitter.site_reports [_site_])
                        }
                    }
                , fitter.table_output_options));

             } else {
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
            }
    ");

    fitter.json [fitter.terms.json.site_reports] = fitter.site_reports;


} else {
    io.ReportProgressMessageMD ("fitter", "sites", "No individual sites showed sufficiently strong preference for multiple-hit models");

}

io.ReportProgressMessageMD ("fitter", "writing", "Writing detailed analysis report to \`" + fitter.codon_data_info [terms.json.json] + "\'");

GetString (_hpv,HYPHY_VERSION,0);
fitter.json[terms.json.runtime] = _hpv;


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
