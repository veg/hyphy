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

namespace bgm {
    LoadFunctionLibrary ("libv3/tasks/bayesgraph.ibf");
    /*
    namespace terms {
        namespace settings {
            nsteps
        }
    }
    */

}




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
   bgm.substitution_model_generator = "models.DNA.GTR.ModelDescription";
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
                                                      bgm.substitution_model_generator ,
                                                      bgm.initial_values,
                                                      bgm.fit_options
                                               );

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

} else {

}

console.log (bgm.counts);

return 0;


// --- execution -------------------------

// load and pre-process codon alignment
namespace bgm {
    LoadFunctionLibrary ("SelectionAnalyses/modules/shared-load-file.bf");
    load_file ("bgm");
}


// fit nucleotide general time-reversible model
// we always re-estimate branch lengths?  Can we constrain to scale?
namespace bgm {
    doGTR ("bgm");
}

// what does this do?
estimators.fixSubsetOfEstimates(bgm.gtr_results, bgm.gtr_results[terms.global]);

bgm.user_tree = bgm.trees["0"];

namespace bgm {
    doPartitionedMG("bgm", TRUE);  // keep LF
}


// --- ancestral reconstruction --------------------
bgm.ancestors = ancestral.build (bgm.partitioned_mg_results[terms.likelihood_function], 0, None);


bgm.code = bgm.codon_data_info[utility.getGlobalValue("terms.code")];



bgm.counts = ancestral.ComputeSubstitutionCounts(
    bgm.ancestors,
    None,  // all branches
    "bgm.nsfilter",  // substitution filter
    "bgm.min_sub_filter"   // site filter (e.g., MinCount)
);



// --- BGM analysis -------------------------------

lfunction bgm.run (_bgm_data, burnin, nsteps, nsamples, max_parents) {
	/* convert data to matrix form */
    nodes = {};
    num_nodes = Abs (_bgm_data["Sites"]);
	for (k = 0; k < num_nodes; k = k+1)
	{
	    /* Arguments:
	        1. node name, must be a string
	        2. maximum number of parents
	        3. prior sample size - always uninformative (count split evenly across levels)
	            - if we were truly Bayesian, we would let the user set informative priors..
	        4. number of levels - always binary in this case (substitution mapped to branch)
	    */
	    node_name = ""+ ((_bgm_data["Sites"])[k] + 1);
	    nodes + add_discrete_node (node_name, max_parents, 0, 2);
	}

    BayesianGraphicalModel gen_bgm = (nodes);
    attach_data(&gen_bgm, _bgm_data["Counts"], 0, 0, 0);
    bgm_result = order_MCMC(&gen_bgm, nsteps, burnin, nsamples);
	return bgm_result;
}


bgm.results = bgm.run (bgm.counts, bgm.burnin, bgm.nsteps, bgm.nsamples, bgm.max_parents);

// --- process BGM results -------------------------------

bgm.trace = {1, bgm.nsamples};  // row vector
for (bgm.i = 0; bgm.i < bgm.nsamples; bgm.i += 1) {
    bgm.trace[bgm.i] = bgm.results[bgm.i][0];
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



