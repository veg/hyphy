RequireVersion("2.26");


VERBOSITY_LEVEL = 0;


LoadFunctionLibrary("libv3/all-terms.bf");

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


console.log("ERROR: DEPRECATED TEST");
exit();

io.DisplayAnalysisBanner({
    terms.io.info: "Load a pre-specified codon file with a tree that has given branch lengths,
    and fit the MG94xREV model (local) while holding branch lengths constant
    ",
    terms.io.version: "1.00",
    terms.io.reference: "GitHub issue #zz",
    terms.io.authors: "Sergei L Kosakovsky Pond",
    terms.io.contact: "spond@temple.edu"
});





codon_data_info = utility.loadGeneticCodeAndAlignment("codon_data", "codon_filter", PATH_TO_CURRENT_BF + "../data/CD2.nex");

io.ReportProgressMessageMD ("", "data", ">Loaded an MSA with **" + codon_data_info[terms.data.sequences] + "** sequences and **" + codon_data_info[terms.data.sites] + "** codons from \`" + codon_data_info[terms.data.file] + "\`");
codon_frequencies = utility.defineFrequencies("codon_filter");
tree = utility.loadAnnotatedTopology(1);

io.ReportProgressMessageMD ("", "fixed", "* Fitting the **fixed branch lengths** model");

fit_options = {
    terms.run_options.model_type: terms.local,
    terms.run_options.proportional_branch_length_scaler: {
    }
};


(fit_options[terms.run_options.proportional_branch_length_scaler])[0] = terms.branch_length_constrain;

constrained_mg_results = estimators.FitMGREV(codon_data_info,
                                             tree,
                                             fit_options ,
                                             parameters.helper.tree_lengths_to_initial_values (tree[terms.branch_lengths], terms.codon)
                                             );

io.ReportProgressMessageMD ("", "fixed",  "* Log(L) = **" + Format(constrained_mg_results[terms.fit.log_likelihood],8,2) + "**");
io.ReportProgressMessageMD ("", "fixed", "*"  + (constrained_mg_results[terms.fit.trees])[0]);

io.ReportProgressMessageMD ("", "free", "* Fitting the **free branch length** model");

fit_options - terms.run_options.proportional_branch_length_scaler;

unconstrained_mg_results = estimators.FitMGREV(codon_data_info, tree, fit_options , constrained_mg_results);
io.ReportProgressMessageMD ("", "free",  "* Log(L) = **" + Format(unconstrained_mg_results[terms.fit.log_likelihood],8,2) + "**");
io.ReportProgressMessageMD ("", "free", "*" + (unconstrained_mg_results[terms.fit.trees])[0]);
