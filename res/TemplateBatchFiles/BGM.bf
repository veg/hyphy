RequireVersion("2.3.13");

// ---- load library files --------------------------------
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");

LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");

LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");

LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary("SelectionAnalyses/modules/selection_lib.ibf");

LoadFunctionLibrary("bayesgraph.ibf");





// --- display analysis information -----------------------

bgm.analysis_description = {
    terms.io.info: "BGM (Bayesian Graphical Model) uses a
    maximum likelihood ancestral state reconstruction to
    map non-synonymous substitution events to branches in the
    phylogeny and then analyzes the joint distribution of the
    substitution map using a Bayesian graphical model (network).
    Next, a Markov chain Monte Carlo analysis is used to generate
    a random sample of network structures from the posterior
    distribution given the data.  Each node in the network
    represents a codon site in the alignment, and links (edges)
    between nodes indicate high posterior support for correlated
    substitutions at the two sites over time, which implies
    coevolution.
    ",
    terms.io.version: "1.0",
    terms.io.reference: "Spidermonkey: rapid detection of co-evolving sites using Bayesian graphical models (2008). _Bioinformatics_ 24(17): 1949-1950",
    terms.io.authors: "Art FY Poon, Fraser I Lewis, Simon DW Frost and Sergei LK Pond",
    terms.io.contact: "apoon42@uwo.ca",
    terms.io.requirements: "in-frame codon alignment and a phylogenetic tree"
};
io.DisplayAnalysisBanner(bgm.analysis_description);



// --- enviornment setup -------------------------



// --- globals ----------------------------

bgm.samples = 10;
bgm.pvalue = 0.1;

bgm.json = {
    terms.json.analysis: bgm.analysis_description,
    terms.json.fits: {},
    terms.json.timers: {},
};

bgm.scaler_prefix = "BGM.scaler";

bgm.by_site = "by-site";
bgm.AVERAGED = "AVERAGED";
bgm.RESOLVED = "RESOLVED";

bgm.nsteps = io.PromptUser("\n>Select the number of MCMC steps to sample [default 100000]", 1e5, 0, 1e9, TRUE);
bgm.burnin = io.PromptUser("\n>Select the number of MCMC steps to discard as burn-in [default 10000]", 1e4, 0, 1e9, TRUE);
bgm.nsamples = io.PromptUser("\n>Select the number of steps to extract from the chain sample [default 100]", 100, 0, bgm.nsteps, TRUE);
bgm.max_parents = io.PromptUser ("\n>Select the maximum number of parents allowed per node [default 1]", 1, 1, 3, TRUE);


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

function bgm.nsfilter(state1, state2, ancestral_data) {
    if (bgm.code[state1] != bgm.code[state2] &&
        bgm.code[state1] != genetic_code.stop_code &&
        bgm.code[state2] != genetic_code.stop_code) {
        return 1;
    } else {
        return 0;
    }
}

bgm.counts = ancestral.ComputeSubstitutionCounts(
    bgm.ancestors,
    None,  // all branches
    "bgm.nsfilter",  // substitution filter
    None   // site filter (e.g., MinCount)
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
    attach_data("gen_bgm", _bgm_data["Counts"], 0, 0, 0);
    bgm_result = order_MCMC("gen_bgm", nsteps, burnin, nsamples);
	return bgm_result;
}


bgm.results = bgm.run (bgm.counts, bgm.burnin, bgm.nsteps, bgm.nsamples, bgm.max_parents);


// --- process BGM results -------------------------------

bgm.trace = {1, bgm.nsamples};  // row vector
for (bgm.i = 0; bgm.i < bgm.nsamples; bgm.i += 1) {
    bgm.trace[bgm.i] = bgm.results[bgm.i][0];
}



