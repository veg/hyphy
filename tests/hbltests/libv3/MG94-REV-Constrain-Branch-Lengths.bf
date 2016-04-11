RequireVersion("2.26");


VERBOSITY_LEVEL = 0;

LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("CF3x4");
LoadFunctionLibrary("TreeTools");


// namespace 'utility' for convenience functions

LoadFunctionLibrary("libv3/UtilityFunctions.bf");

// namespace 'io' for interactive/datamonkey i/o functions
LoadFunctionLibrary("libv3/IOFunctions.bf");

LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("libv3/tasks/estimators.bf");

// namespace 'ancestral' for ancestral reconstruction functions
LoadFunctionLibrary("libv3/tasks/ancestral.bf");

// namespace 'stats' for various descriptive stats functions
LoadFunctionLibrary("libv3/stats.bf");


/*------------------------------------------------------------------------------*/


io.displayAnalysisBanner({
    "info": "Load a pre-specified codon file with a tree that has given branch lengths,
    and fit the MG94xREV model (local) while holding branch lengths constant
    ",
    "version": "1.00",
    "reference": "GitHub issue #zz",
    "authors": "Sergei L Kosakovsky Pond",
    "contact": "spond@temple.edu"
});





codon_data_info = utility.loadGeneticCodeAndAlignment("codon_data", "codon_filter",PATH_TO_CURRENT_BF + "../data/CD2.nex");

io.reportProgressMessageMD ("", "data", ">Loaded an MSA with **" + codon_data_info["sequences"] + "** sequences and **" + codon_data_info["sites"] + "** codons from \`" + codon_data_info["file"] + "\`");
codon_frequencies = utility.defineFrequencies("codon_filter");
tree = utility.loadAnnotatedTopology(1);

io.reportProgressMessageMD ("", "fixed", "* Fitting the **fixed branch lengths** model");

fit_options = {
    "model-type": terms.local,
    "proportional-branch-length-scaler": {
    }
};


(fit_options["proportional-branch-length-scaler"])[0] = terms.branch_length_constrain;

constrained_mg_results = estimators.fitMGREV(codon_data_info,
                                             tree,
                                             fit_options ,
                                             parameters.helper.tree_lengths_to_initial_values (tree["branch_lengths"], "codon")
                                             );

io.reportProgressMessageMD ("", "fixed",  "* Log(L) = **" + Format(constrained_mg_results["LogL"],8,2) + "**");
io.reportProgressMessageMD ("", "fixed", "*"  + (constrained_mg_results["Trees"])[0]);

io.reportProgressMessageMD ("", "free", "* Fitting the **free branch length** model");

fit_options - "proportional-branch-length-scaler";

unconstrained_mg_results = estimators.fitMGREV(codon_data_info, tree, fit_options , constrained_mg_results);
io.reportProgressMessageMD ("", "free",  "* Log(L) = **" + Format(unconstrained_mg_results["LogL"],8,2) + "**");
io.reportProgressMessageMD ("", "free", "*" + (unconstrained_mg_results["Trees"])[0]);
