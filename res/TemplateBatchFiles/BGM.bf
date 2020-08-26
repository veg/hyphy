RequireVersion("2.4.0");

// ---- load library files --------------------------------
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");

LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary ("libv3/tasks/bayesgraph.ibf");

// --- display analysis information -----------------------

bgm.analysis_description = {
    terms.io.info:
"BGM (Bayesian Graphical Model) uses a maximum likelihood ancestral state reconstruction to
map substitution (non-synonymous only for coding data) events to branches in the
phylogeny and then analyzes the joint distribution of the
substitution map using a Bayesian graphical model (network).
Next, a Markov chain Monte Carlo analysis is used to generate
a random sample of network structures from the posterior
distribution given the data.  Each node in the network
represents a codon site in the alignment, and links (edges)
between nodes indicate high posterior support for correlated
substitutions at the two sites over time, which implies coevolution.",
    terms.io.version: "1.2",
    terms.io.reference: "Spidermonkey: rapid detection of co-evolving sites using Bayesian graphical models (2008). _Bioinformatics_ 24(17): 1949-1950",
    terms.io.authors: "Art FY Poon, Fraser I Lewis, Simon DW Frost and Sergei L Kosakovsky Pond",
    terms.io.contact: "apoon42@uwo.ca",
    terms.io.requirements: "in-frame codon alignment and a phylogenetic tree"
};

io.DisplayAnalysisBanner(bgm.analysis_description);

namespace bgm {
    LoadFunctionLibrary ("SelectionAnalyses/modules/shared-load-file.bf");
}

bgm.table_output_options = {terms.table_options.header : TRUE, terms.table_options.minimum_column_width: 12, terms.table_options.align : "center"};

bgm.json = {
    terms.json.analysis: bgm.analysis_description,
    terms.json.fits: {},
    terms.json.timers: {},
};

selection.io.startTimer (bgm.json [terms.json.timers], "Overall", 0);


bgm.data_types = {terms.nucleotide  : "Nucleotide multiple sequence alignment",
                 terms.amino_acid   : "Protein multiple sequence alignment",
                 terms.codon        : "Codon multiple sequence alignment"};

KeywordArgument ("type", "nucleotide, amino-acid or codon", "codon");
bgm.type = io.SelectAnOption (bgm.data_types, "Data type");

SetDialogPrompt ("Specify a `bgm.type` multiple sequence alignment file");

bgm.fit_options = {terms.run_options.retain_lf_object : TRUE};
bgm.reporting_thershold = 0.5;

bgm.run_settings = {
    "steps" : 1e5,
    "burn-in" : 1e4,
    "samples" : 100,
    "max-parents" : 1,
    "min-subs" : 1,
    "data-type" : bgm.type,
    "threshold" : bgm.reporting_thershold
};

KeywordArgument ("code", "Which genetic code should be used", "Universal");
KeywordArgument ("alignment", "An in-frame codon alignment in one of the formats supported by HyPhy");

if (bgm.type == "nucleotide") {
   bgm.alignment_info = alignments.ReadNucleotideDataSet ("bgm.dataset", None);
   bgm.baseline_model = "models.DNA.GTR.ModelDescription";
} else {
    if (bgm.type == "amino-acid") {
        bgm.alignment_info = alignments.ReadProteinDataSet ("bgm.dataset", None);
        LoadFunctionLibrary     ("libv3/models/protein.bf");
        LoadFunctionLibrary     ("libv3/models/protein/empirical.bf");
        LoadFunctionLibrary     ("libv3/models/protein/REV.bf");
        utility.Extend (models.protein.empirical_models, {"GTR" : "General time reversible model (189 estimated parameters)."});
        KeywordArgument ("baseline_model", "Which amino acid substitution model should be used", "LG");
        bgm.run_settings ["model"] = io.SelectAnOption (models.protein.empirical_models, "Baseline substitution model");
        bgm.baseline_model = (utility.Extend (models.protein.empirical.plusF_generators , {"GTR" : "models.protein.REV.ModelDescription"}))[bgm.run_settings ["model"]];
    } else { // codon
        bgm.alignment_info = alignments.PromptForGeneticCodeAndAlignment("bgm.dataset","bgm.codon.filter");
        LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");
        bgm.baseline_model = "models.codon.MG_REV.ModelDescription";
        function bgm.MG94_REV (options) {
            return Call ("models.codon.MG_REV.ModelDescription", options, bgm.alignment_info[terms.code]);
        }

    }
}

bgm.name_mapping = bgm.alignment_info[utility.getGlobalValue("terms.data.name_mapping")];


selection.io.json_store_key_value_pair (bgm.json, terms.json.input, terms.json.file, bgm.alignment_info [terms.data.file]);
selection.io.json_store_key_value_pair (bgm.json, terms.json.input, terms.json.sequences, bgm.alignment_info [terms.data.sequences]);
selection.io.json_store_key_value_pair (bgm.json, terms.json.input, terms.json.sites, bgm.alignment_info [terms.data.sites]);
selection.io.json_store_key_value_pair (bgm.json, terms.json.input, terms.data_type, bgm.type);

bgm.alignment_info[terms.json.json] = bgm.alignment_info[terms.data.file] + ".BGM.json";

bgm.path.base = (bgm.json [terms.json.input])[terms.json.file];


bgm.sample_size = bgm.alignment_info[terms.data.sequences] * bgm.alignment_info[terms.data.sites];
alignments.EnsureMapping ("bgm.dataset", bgm.alignment_info);

KeywordArgument ("tree", "A phylogenetic tree (optionally annotated with {})", null, "Please select a tree file for the data:");
KeywordArgument ("branches",  "Branches to test", "All");
bgm.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (
                                                                             bgm.alignment_info[utility.getGlobalValue("terms.data.partitions")],
                                                                             bgm.alignment_info[utility.getGlobalValue("terms.data.name_mapping")]
                                                                             );


io.CheckAssertion ("Abs (bgm.partitions_and_trees) == 1", "BGM cannot be run on data with multiple site partitions (and trees)");
bgm.filter_specification = alignments.DefineFiltersForPartitions (bgm.partitions_and_trees, "bgm.dataset" , "bgm.filter.", bgm.alignment_info);
bgm.store_tree_information();

io.ReportProgressMessageMD ("BGM", "Data", "Loaded **" +
                            bgm.alignment_info [terms.data.sequences] + "** `bgm.type` sequences, **" +
                            bgm.alignment_info [terms.data.sites] + "** sites, from \`" + bgm.alignment_info [terms.data.file] + "\`");

bgm.initial_values = parameters.helper.tree_lengths_to_initial_values (bgm.trees, None);



console.log ( "\n> BGM will write result file to `bgm.alignment_info[terms.json.json]`\n");

bgm.selected_branches = selection.io.defineBranchSets ( bgm.partitions_and_trees );


KeywordArgument ("steps", "The number of MCMC steps to sample", 100000);
KeywordArgument ("burn-in", "The number of MCMC steps to discard as burn-in", 10000);
KeywordArgument ("samples", "The number of steps to extract from the chain sample", 100);
KeywordArgument ("max-parents", "The maximum number of parents allowed per node", 1);
KeywordArgument ("min-subs", "The minium number of substitutions per site to include it in the analysis", 1);


bgm.run_settings["steps"]      = io.PromptUser("\n>Select the number of MCMC steps to sample", bgm.run_settings["steps"] , 0, 1e9, TRUE);
bgm.run_settings["burn-in"]      = io.PromptUser("\n>Select the number of MCMC steps to discard as burn-in", bgm.run_settings["burn-in"], 0, 1e9, TRUE);
bgm.run_settings["samples"]    = io.PromptUser("\n>Select the number of steps to extract from the chain sample", 100, 0, bgm.run_settings["samples"], TRUE);
bgm.run_settings["max-parents"] = io.PromptUser ("\n>Select the maximum number of parents allowed per node", bgm.run_settings["max-parents"], 1, 3, TRUE);
bgm.run_settings["min-subs"]    = io.PromptUser ("\n>Select the minium number of substitutions per site to include it in the analysis", bgm.run_settings["min-subs"], 1, 1e5, TRUE);


// FIT THE BASELINE MODEL

selection.io.startTimer (bgm.json [terms.json.timers], "Baseline fit", 1);

io.ReportProgressMessageMD("bgm", "phylo", "Performing initial model fit to obtain branch lengths and rate parameters");

if (bgm.type == "nucleotide") {
   bgm.initial_values = utility.Extend (bgm.initial_values,
                                  {
                                    utility.getGlobalValue ("terms.global") : {
                                        terms.nucleotideRate ("A","C") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25} ,
                                        terms.nucleotideRate ("A","T") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25} ,
                                        terms.nucleotideRate ("C","G") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25} ,
                                        terms.nucleotideRate ("G","T") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25}
                                    }
                                 });

 } else {
    if (bgm.type == "codon") {
        bgm.initial_values = utility.Extend (bgm.initial_values,
                                  {
                                    utility.getGlobalValue ("terms.global") : {
                                        terms.nucleotideRate ("A","C") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25} ,
                                        terms.nucleotideRate ("A","T") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25} ,
                                        terms.nucleotideRate ("C","G") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25} ,
                                        terms.nucleotideRate ("G","T") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25} ,
                                        terms.omega_ratio : { utility.getGlobalValue ("terms.fit.MLE") : 0.25}
                                    }
                                 });

        io.ReportProgressMessageMD("bgm", "phylo", "Fitting nucleotide GTR to obtain branch length estimates");

        bgm.initial_values = estimators.FitSingleModel_Ext (
                                                          bgm.filter_names,
                                                          bgm.trees,
                                                          "models.DNA.GTR.ModelDescription" ,
                                                          bgm.initial_values,
                                                          bgm.fit_options
                                                   );

    io.ReportProgressMessageMD ("bgm", "phylo", ">Fitted an alignment-wide GTR model. " + selection.io.report_fit (bgm.initial_values, 0, bgm.sample_size ));

        bgm.filter_names = { "0" : "bgm.codon.filter" };
        bgm.code = bgm.alignment_info [terms.code];
    }
 }

if (bgm.type == "codon") {
    //codon_data, tree, generator, genetic_code, option, initial_values
    bgm.baseline_fit = estimators.FitCodonModel(
                                                     bgm.filter_names,
                                                     bgm.trees,
                                                     bgm.baseline_model ,
                                                     bgm.code,
                                                     bgm.fit_options,
                                                     bgm.initial_values,

                                               );
} else {
    bgm.baseline_fit = estimators.FitSingleModel_Ext (
                                                          bgm.filter_names,
                                                          bgm.trees,
                                                          bgm.baseline_model ,
                                                          bgm.initial_values,
                                                          bgm.fit_options
                                                   );
}


io.ReportProgressMessageMD ("bgm", "phylo", ">Fitted an alignment-wide model. " + selection.io.report_fit (bgm.baseline_fit, 0, bgm.sample_size ) +  "\n\nTotal tree lengths by partition\n");
utility.ForEachPair (bgm.baseline_fit[terms.branch_length], "_part_", "_value_",
'
    io.ReportProgressMessageMD ("bgm", "phylo", "Partition " + (1+_part_) + ". " + Format (+(utility.Map (_value_, "_data_",
    "
        _data_ [terms.fit.MLE]
    "))
    ,6,3) + " subs/site."
    )
'
);

selection.io.json_store_lf(bgm.json, bgm.baseline_model,bgm.baseline_fit[terms.fit.log_likelihood],
                            bgm.baseline_fit[terms.parameters],
                            bgm.sample_size, None, 0);

selection.io.stopTimer (bgm.json [terms.json.timers], "Baseline fit");

utility.ForEachPair (bgm.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(bgm.json, bgm.baseline_model, terms.branch_length, 0,
                                     _key_,
                                     selection.io.extract_branch_info((bgm.baseline_fit[terms.branch_length])[_key_], "selection.io.branch.length"));');



io.ReportProgressMessageMD("bgm", "ancestral", "Performing joint ancestral state reconstruction and mapping substitutions");

bgm.ancestral_cache = ancestral.build (bgm.baseline_fit[terms.likelihood_function], 0, None);
bgm.branch_filter = utility.Filter (bgm.selected_branches[0], "_class_", "_class_ == terms.tree_attributes.test");
DeleteObject (^bgm.baseline_fit[terms.likelihood_function]);

if (bgm.type != "codon") {
   bgm.counts = ancestral.ComputeSubstitutionCounts(
        bgm.ancestral_cache,
        bgm.branch_filter,  // selected branches
        None,  // substitution filter
        "bgm.min_sub_filter"   // site filter (e.g., MinCount)
    );

} else {
    bgm.counts = ancestral.ComputeSubstitutionCounts(
        bgm.ancestral_cache,
        bgm.branch_filter,  // selected branches
        "bgm.nsfilter",  // substitution filter
        "bgm.min_sub_filter"   // site filter (e.g., MinCount)
    );
}

selection.io.startTimer (bgm.json [terms.json.timers], "Network inference", 2);

bgm.site_count = utility.Array1D (bgm.counts["Sites"]);

if (bgm.site_count <= 2) {
    console.log ("###ERROR: NOT ENOUGH SUBSTITUTIONS###");
    console.log ("\n>BGM requires at least three sites to have accumulated sufficient substitutions to run network inference");
} else {
    io.ReportProgressMessageMD("bgm", "inference", "Inferring a BGM on `bgm.site_count` nodes [sites]");

    namespace bgm {

        raw_results = bgm.run (counts, run_settings["burn-in"], run_settings["steps"], run_settings["samples"], run_settings["max-parents"]);
        mcmc_trace = utility.Map ({1, run_settings["samples"]}["_MATRIX_ELEMENT_COLUMN_"], "_index_", "bgm.raw_results[_index_][0]");

    // unpack site results

        table_headers = {
                     {"Site 1", "Index of site 1"}
                     {"Site 2", "Index of site 2"}
                     {"P [Site 1 –> Site 2]", "Probability that site 2 is conditionally dependent on site 1"}
                     {"P [Site 2 –> Site 1]", "Probability that site 1 is conditionally dependent on site 2"}
                     {"P [Site 1 <–> Site 2]", "Probability that sites 1 and 2 are not conditionally independent"}
                     {"Site 1 subs", "Substitution counts inferred for Site 1"}
                     {"Site 2 subs", "Substitution counts inferred for Site 2"}
                     {"Shared subs", "Substitutions shared by both sites"}
                    };

        processed_results = {site_count * (site_count - 1) / 2, 8};
        row_index = 0;

        table_screen_output  = {{"Site 1", "Site 2", "P [Site 1 <-> Site 2]", "Subs (1,2,shared)"}};

        report.coupled_pair = {{Format(processed_results[row_index][0],5,0),
                                Format(processed_results[row_index][1],5,0),
                                Format(processed_results[row_index][4],6,3),
                                "" + processed_results[row_index][5] + ", " + processed_results[row_index][6] + ", " + processed_results[row_index][7]}};

        report.pairs_found = {};
        for (site1 = 0; site1 < site_count; site1 += 1) {
             for (site2 = site1 + 1; site2 < site_count; site2 += 1) {
                processed_results[row_index][0] = 0 + (counts["Sites"])[site1] + 1;
                processed_results[row_index][1] = 0 + (counts["Sites"])[site2] + 1;
                processed_results[row_index][2] = raw_results[site1 * site_count + site2][1];
                processed_results[row_index][3] = raw_results[site2 * site_count + site1][1];
                processed_results[row_index][4] = processed_results[row_index][2] + processed_results[row_index][3];
                processed_results[row_index][5] = +(counts["Counts"])[-1][site1];
                processed_results[row_index][6] = +(counts["Counts"])[-1][site2];
                processed_results[row_index][7] = +((counts["Counts"])[-1][site2]$(counts["Counts"])[-1][site1]);

                if (processed_results[row_index][4] >= reporting_thershold) {
                    if (Abs(report.pairs_found) == 0 && table_output_options[^"terms.table_options.header"]) {
                        fprintf (stdout, "\n", io.FormatTableRow (table_screen_output,table_output_options));
                        table_output_options[^"terms.table_options.header"] = FALSE;
                    }
                    fprintf (stdout, io.FormatTableRow (report.coupled_pair,table_output_options));
                    report.pairs_found + (processed_results[row_index][-1]);
                }

                row_index+=1;
            }
        }

        json [^"terms.fit.MLE"] = {^"terms.json.headers"   : table_headers,
                               ^"terms.json.content" : processed_results };


        console.log ("----\n## BGM analysis summary on `site_count` sites each with at least `run_settings['min-subs']` substitutions. Evidence for conditional dependence was reported at posterior probability of " + reporting_thershold);

        _count_ = Abs(report.pairs_found);

        if (_count_ == 0) {
            console.log ("* No pairs of conditionally independent sites found");
        } else {
            console.log ("* " + _count_ + " " + io.SingularOrPlural (_count_, "pair ", "pairs ") + " of conditionally dependent sites found");
        }


    }
}

selection.io.stopTimer (bgm.json [terms.json.timers], "Network inference");
selection.io.stopTimer (bgm.json [terms.json.timers], "Overall");

bgm.json [terms.settings] = bgm.run_settings;

io.SpoolJSON (bgm.json, bgm.alignment_info[terms.json.json]);


// --- BGM analysis -------------------------------

lfunction bgm.run (_bgm_data, burnin, nsteps, nsamples, max_parents) {
	/* convert data to matrix form */
    nodes = {};
    num_nodes = Abs (_bgm_data["Sites"]);
	for (k = 0; k < num_nodes; k += 1) {
	    /* Arguments:
	        1. node name, must be a string
	        2. maximum number of parents
	        3. prior sample size - always uninformative (count split evenly across levels)
	            - if we were truly Bayesian, we would let the user set informative priors.
	        4. number of levels - always binary in this case (substitution mapped to branch)
	    */
	    node_name = ""+ ((_bgm_data["Sites"])[k] + 1);
	    nodes + bgm.add_discrete_node (node_name, max_parents, 0, 2);
	}

    utility.ToggleEnvVariable ("VERBOSITY_LEVEL",0);

    BayesianGraphicalModel gen_bgm = (nodes);
    bgm.attach_data (&gen_bgm, _bgm_data["Counts"], 0, 0, 0);
    bgm_result = bgm.order_MCMC(&gen_bgm, nsteps, burnin, nsamples);

    utility.ToggleEnvVariable ("VERBOSITY_LEVEL",None);


	return bgm_result;
}


// ==== HELPER FUNCTIONS ====

function bgm.nsfilter(state1, state2, ancestral_data) {

    if (state1 >= 0 && state2 >= 0) {
        if (bgm.code[state1] != bgm.code[state2] &&
            bgm.code[state1] != genetic_code.stop_code &&
            bgm.code[state2] != genetic_code.stop_code) {

            return TRUE;
        }
    }
    return FALSE;
}

function bgm.min_sub_filter (counts) {
    return (+counts) >= bgm.run_settings["min-subs"];
}



