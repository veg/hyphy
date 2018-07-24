RequireVersion("2.3.12");

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
    terms.io.version: "1.0",
    terms.io.reference: "Spidermonkey: rapid detection of co-evolving sites using Bayesian graphical models (2008). _Bioinformatics_ 24(17): 1949-1950",
    terms.io.authors: "Art FY Poon, Fraser I Lewis, Simon DW Frost and Sergei L Kosakovsky Pond",
    terms.io.contact: "apoon42@uwo.ca",
    terms.io.requirements: "in-frame codon alignment and a phylogenetic tree"
};
io.DisplayAnalysisBanner(bgm.analysis_description);


namespace bgm {
    LoadFunctionLibrary ("SelectionAnalyses/modules/shared-load-file.bf");
}

bgm.json = {
    terms.json.analysis: bgm.analysis_description,
    terms.json.fits: {},
    terms.json.timers: {},
};

bgm.data_types = {terms.nucleotide  : "Nucleotide multiple sequence alignment",
                 terms.amino_acid   : "Protein multiple sequence alignment",
                 terms.codon        : "Codon multiple sequence alignment"};

bgm.run_type = io.SelectAnOption (bgm.data_types, "Data type");

SetDialogPrompt ("Specify a `bgm.run_type` multiple sequence alignment file");

bgm.fit_options = {terms.run_options.retain_lf_object : TRUE};


if (bgm.run_type == "nucleotide") {
   bgm.alignment_info = alignments.ReadNucleotideDataSet ("bgm.dataset", None);
   bgm.baseline_model = "models.DNA.GTR.ModelDescription";
} else {
    if (bgm.run_type == "amino-acid") {

    } else { // codon
        LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");
    }
}

bgm.name_mapping = bgm.alignment_info[utility.getGlobalValue("terms.data.name_mapping")];


selection.io.json_store_key_value_pair (bgm.json, terms.json.input, terms.json.file, bgm.alignment_info [terms.data.file]);
selection.io.json_store_key_value_pair (bgm.json, terms.json.input, terms.json.sequences, bgm.alignment_info [terms.data.sequences]);
selection.io.json_store_key_value_pair (bgm.json, terms.json.input, terms.json.sites, bgm.alignment_info [terms.data.sites]);
selection.io.json_store_key_value_pair (bgm.json, terms.json.input, terms.data_type, bgm.run_type);

bgm.alignment_info[terms.json.json] = bgm.alignment_info[terms.data.file] + ".BGM.json";

bgm.path.base = (bgm.json [terms.json.input])[terms.json.file];


bgm.sample_size = bgm.alignment_info[terms.data.sequences] * bgm.alignment_info[terms.data.sites];
alignments.EnsureMapping ("bgm.dataset", bgm.alignment_info);

bgm.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (
                                                                             bgm.alignment_info[utility.getGlobalValue("terms.data.partitions")],
                                                                             bgm.alignment_info[utility.getGlobalValue("terms.data.name_mapping")]
                                                                             );


io.CheckAssertion ("Abs (bgm.partitions_and_trees) == 1", "BGM cannot be run on data with multiple site partitions (and trees)");
bgm.filter_specification = alignments.DefineFiltersForPartitions (bgm.partitions_and_trees, "bgm.dataset" , "bgm.filter.", bgm.alignment_info);
bgm.store_tree_information();

io.ReportProgressMessageMD ("BGM", "Data", "Loaded **" +
                            bgm.alignment_info [terms.data.sequences] + "** `bgm.run_type` sequences, **" +
                            bgm.alignment_info [terms.data.sites] + "** sites, from \`" + bgm.alignment_info [terms.data.file] + "\`");

bgm.initial_values = parameters.helper.tree_lengths_to_initial_values (bgm.trees, None);



console.log ( "\n> BGM will write result file to `bgm.alignment_info[terms.json.json]`\n");

bgm.selected_branches = selection.io.defineBranchSets ( bgm.partitions_and_trees );

bgm.nsteps      = io.PromptUser("\n>Select the number of MCMC steps to sample [default 100000]", 1e5, 0, 1e9, TRUE);
bgm.burnin      = io.PromptUser("\n>Select the number of MCMC steps to discard as burn-in [default 10000]", 1e4, 0, 1e9, TRUE);
bgm.nsamples    = io.PromptUser("\n>Select the number of steps to extract from the chain sample [default 100]", 100, 0, bgm.nsteps, TRUE);
bgm.max_parents = io.PromptUser ("\n>Select the maximum number of parents allowed per node [default 1]", 1, 1, 3, TRUE);
bgm.min_subs    = io.PromptUser ("\n>Select the minium number of substitutions per site to include it in the analysis", 1, 1, 1e5, TRUE);


// FIT THE BASELINE MODEL

io.ReportProgressMessageMD("bgm", "phylo", "Performing initial model fit to obtain branch lengths and rate parameters");

if (bgm.run_type == "nucleotide") {
   bgm.initial_values = utility.Extend (bgm.initial_values,
                                  {
                                    utility.getGlobalValue ("terms.global") : {
                                        terms.nucleotideRate ("A","C") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25} ,
                                        terms.nucleotideRate ("A","T") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25} ,
                                        terms.nucleotideRate ("C","G") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25} ,
                                        terms.nucleotideRate ("G","T") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25}
                                    }
                                 });

 }

bgm.baseline_fit = estimators.FitSingleModel_Ext (
                                                      bgm.filter_names,
                                                      bgm.trees,
                                                      bgm.baseline_model ,
                                                      bgm.initial_values,
                                                      bgm.fit_options
                                               );


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

utility.ForEachPair (bgm.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(bgm.json, bgm.baseline_model, terms.branch_length, 0,
                                     _key_,
                                     selection.io.extract_branch_info((bgm.baseline_fit[terms.branch_length])[_key_], "selection.io.branch.length"));');


io.ReportProgressMessageMD("bgm", "ancestral", "Performing joint ancestral state reconstruction and mapping substitutions");

bgm.ancestral_cache = ancestral.build (bgm.baseline_fit[terms.likelihood_function], 0, None);
bgm.branch_filter = utility.Filter (bgm.selected_branches[0], "_class_", "_class_ == terms.tree_attributes.test");


if (bgm.run_type != "codon") {
   bgm.counts = ancestral.ComputeSubstitutionCounts(
        bgm.ancestral_cache,
        bgm.branch_filter,  // selected branches
        None,  // substitution filter
        "bgm.min_sub_filter"   // site filter (e.g., MinCount)
    );
} else {

}

if (Abs (bgm.counts["Sites"]) <= 2) {
    console.log ("###ERROR: NOT ENOUGH SUBSTITUTIONS###");
    console.log ("\n>BGM requires at least three sites to have accumulated sufficient substitutions to run network inference");
} else {
    bgm.site_count (bgm.counts);

    bgm.raw_results = bgm.run (bgm.counts, bgm.burnin, bgm.nsteps, bgm.nsamples, bgm.max_parents);
    bgm.mcmc_trace = utility.Map ({1, bgm.nsamples}["_MATRIX_ELEMENT_COLUMN_"], "_index_", "bgm.raw_results[_index_][0]");



    //bgm.processed_results = {
}


return 0;

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
	            - if we were truly Bayesian, we would let the user set informative priors..
	        4. number of levels - always binary in this case (substitution mapped to branch)
	    */
	    node_name = ""+ ((_bgm_data["Sites"])[k] + 1);
	    nodes + bgm.add_discrete_node (node_name, max_parents, 0, 2);
	}


    BayesianGraphicalModel gen_bgm = (nodes);
    bgm.attach_data (&gen_bgm, _bgm_data["Counts"], 0, 0, 0);
    bgm_result = bgm.order_MCMC(&gen_bgm, nsteps, burnin, nsamples);

    console.log (Rows(bgm_result));

	return bgm_result;
}


// ==== HELPER FUNCTIONS ====

lfunction bgm.nsfilter(state1, state2, ancestral_data) {
    if (bgm.code[state1] != bgm.code[state2] &&
        bgm.code[state1] != genetic_code.stop_code &&
        bgm.code[state2] != genetic_code.stop_code) {
        return TRUE;
    } else {
        return FALSE;
    }
}

function bgm.min_sub_filter (counts) {
    return (+counts) > bgm.min_subs;
}



