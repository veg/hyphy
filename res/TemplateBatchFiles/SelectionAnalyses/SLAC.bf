RequireVersion("2.25");


VERBOSITY_LEVEL = 0;
//LF_SMOOTHING_SCALER         = 0.1;

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

// namespace 'selection.io' for convenience functions 
LoadFunctionLibrary("modules/io_functions.ibf");

// namespace 'ancestral' for ancestral reconstruction functions 
LoadFunctionLibrary("libv3/tasks/ancestral.bf");


slac.timers = {
    6, 1
};

selection.io.startTimer(slac.timers, 0);

slac.json = {
    "fits": {},
    "timers": {},
};

slac.scaler = "SLAC.scaler";

/*------------------------------------------------------------------------------*/


io.displayAnalysisBanner({
    "info": "SLAC (Single Likelihood Ancestor Counting)
    uses a maximum likelihood ancestral state reconstruction
    and minimum path substitution counting to estimate site - level
    dS and dN,
    and applies a simple binomial - based test to test
    if dS differs drom dN.
    The estimates aggregate information over all branches,
    so the signal is derived from
    pervasive diversification or conservation.A subset of branches can be selected
    for testing as well.
    Multiple partitions within a NEXUS file or multiple files are also supported
    for recombination - aware analysis.
    ",
    "version": "2.00",
    "reference": "Not So Different After All: A Comparison of Methods for Detecting Amino Acid Sites Under Selection (2005). Mol Biol Evol 22 (5): 1208-1222",
    "authors": "Sergei L Kosakovsky Pond and Simon DW Frost",
    "contact": "spond@temple.edu",
    "requirements": "in-frame codon alignment and a phylogenetic tree"
});





slac.codon_data_info = utility.promptForGeneticCodeAndAlignment("slac.codon_data", "slac.codon_filter");
slac.sample_size = slac.codon_data_info["sites"] * slac.codon_data_info["sequences"];
slac.codon_data_info["json"] = slac.codon_data_info["file"] + ".slac.json";


io.reportProgressMessage("SLAC", "Loaded an MSA with " + slac.codon_data_info["sequences"] + " sequences and " + slac.codon_data_info["sites"] + " codons from '" + slac.codon_data_info["file"] + "'");

slac.codon_frequencies = utility.defineFrequencies("slac.codon_filter");
slac.tree = utility.loadAnnotatedTopology(1);

slac.selected_branches = selection.io.defineBranchSets(slac.tree);
/* dict ["selected branch (original case)"] == TRUE */ 

slac.json["tested"] = slac.selected_branches;
slac.json["tree"] = slac.tree["string"];

io.reportProgressMessage("SLAC", "Selected " + Abs(slac.selected_branches) + " branches to test for selection " + Join(",", Rows(slac.selected_branches)));

selection.io.startTimer(slac.timers, 1);


io.reportProgressMessage("SLAC", "Obtaining branch lengths under the nucleotide GTR model");
slac.gtr_results = estimators.fitGTR("slac.codon_filter", slac.tree, None);
io.reportProgressMessage("SLAC", "Log(L) = " + slac.gtr_results["LogL"]);
estimators.fixSubsetOfEstimates(slac.gtr_results, slac.gtr_results["global"]);

io.reportProgressMessage("SLAC", "Obtaining the global omega estimate using relative GTR branch lengths");

parameters.declareGlobal(slac.scaler, None);
parameters.set_value(slac.scaler, 3);

slac.global_mg_results = estimators.fitMGREV(slac.codon_data_info, slac.tree, {
    "model-type": terms.global,
    "proportional-branch-length-scaler": {
        "0": slac.scaler
    },
    "retain-lf-object": TRUE
}, slac.gtr_results);
io.reportProgressMessage("SLAC", "Log(L) = " + slac.global_mg_results["LogL"]);

slac.global_dnds = selection.io.extract_global_MLE(slac.global_mg_results, terms.omega_ratio);
io.reportProgressMessage("SLAC", "Global dN/dS = " + slac.global_dnds);

selection.io.stopTimer(slac.timers, 1);

selection.io.json_store_lf(
    slac.json,
    "Global MG94xREV",
    slac.global_mg_results["LogL"],
    slac.global_mg_results["parameters"] - 1,
    slac.sample_size,
    slac.timers[1],
    selection.io.extract_branch_info((slac.global_mg_results["branch lengths"])[0], "selection.io.branch.length"),
    selection.io.extract_branch_info((slac.global_mg_results["branch lengths"])[0], "selection.io.branch.length"), {
        "All": {
            {
                slac.global_dnds, 1
            }
        }
    },
    "Substitution Counts"
);

io.spoolLF (slac.global_mg_results["LF"], slac.codon_data_info["file"], None);

/*------------------------------------------------------------------------------------*/

lfunction slac.compute_the_counts (matrix, tree, lookup, selected_branches, counts) {
    site_count = Rows (matrix);
    selected_branches_count  = Abs (selected_branches);
    selected_branches_in_avl = {selected_branches_count,1};
    selected_branch_total_length = 0;
    k = 0;
    
    for (i = 1; i < Abs (tree); i+=1) {
        if (selected_branches [(tree[i])["Name"]&&1]) {
            selected_branches_in_avl[k]   = i-1;
            selected_branch_total_length += (tree[i])["Length"];
            k+=1;
        }
    }
    
    /*  columns 
           0 Expected synonymous     sites 
           1 Expected non-synonymous sites 
           2 Observed synonymous subs
           3 Observed non-synonymous subs
           4 Expected ratio
           5 Observed ratio 
           6 p-value
    */
    
    report_resolved = {site_count, 7};
    report_averaged = {site_count, 7};
    
    for (i = 0; i < selected_branches_count; i+=1) {
        this_branch            = selected_branches_in_avl[i];
        if ((tree[this_branch])["Length"]) {
            relative_branch_length = (tree[this_branch])["Length"] / selected_branch_total_length;
            parent_branch          = (tree[this_branch])["Parent"] - 1;
            for (s = 0; s < sites; s += 1) {
                this_state      = matrix[this_branch][s];
                parent_state    = matrix[parent_branch][s];
                
                // parent state can only be resolved or --- (-1)
                
                if (this_state >= 0) { // tip fully resolved
                    if (parent_state >= 0) { // parent fully resolved 
                } 
            }
        }
        
    }
}

/*------------------------------------------------------------------------------------*/


io.reportProgressMessage("SLAC", "Performing joint maximum likelihood ancestral state reconstruction");

slac.counts    = genetic_code.ComputePairwiseDifferencesAndExpectedSites (slac.codon_data_info["code"], {"count-stop-codons" : FALSE});
slac.ancestors = ancestral.build (slac.global_mg_results["LF"], 0, None);

fprintf (stdout, slac.compute_the_counts (slac.ancestors["MATRIX"], slac.ancestors["TREE_AVL"], slac.ancestors["AMBIGS"], slac.selected_branches, slac.counts), "\n");

