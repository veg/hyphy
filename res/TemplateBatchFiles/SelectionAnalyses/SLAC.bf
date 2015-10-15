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
#profile START;

    site_count = Columns (matrix);
    
    selected_branches_count      = Abs (selected_branches);
    selected_branches_in_avl     = {selected_branches_count,1};
    selected_branches_lengths    = {selected_branches_count,1};
    selected_branches_parents    = {selected_branches_count,1};
    selected_branch_total_length = 0;
    k = 0;
    tip_count = 0;
    
    for (i = 1; i < Abs (tree); i+=1) {
        if (selected_branches [(tree[i])["Name"]&&1]) {
            selected_branches_in_avl[k]   = i-1;
            selected_branches_lengths[k] = (tree[i])["Length"];
            selected_branches_parents[k] = (tree[i])["Parent"] - 1;
            
            k+=1;
        }
        if (Abs ((tree[i])["Children"]) == 0) {
            tip_count += 1;        
        }
    }
    
    
    selected_branch_total_length  = +selected_branches_lengths;
    io.checkAssertion ("`&selected_branch_total_length`>0", "SLAC cannot be applied to a branch selection with total zero branch length (i.e. no variation)");
    
    /*  columns 
           0 Expected synonymous     sites 
           1 Expected non-synonymous sites 
           2 Observed synonymous subs
           3 Observed non-synonymous subs
           4 Expected ratio
           5 Observed ratio 
           6 dS
           7 dN
           8 dN-dS 
           9 dN-dS (scaled)
           10 P {S <= observed} // posititve sel
           11 P {S >= observed} // negative sel 
           12 effective branch length
           
    */
    
    column_count    = 13;
    report_resolved = {site_count, column_count};
    report_averaged = {site_count, column_count};
    pairwise_eps = counts ["EPS"];
    state_count  = Rows (pairwise_eps);
    pairwise_epn = counts ["EPN"];
    pairwise_ops = counts ["OPS"];
    pairwise_opn = counts ["OPN"];
    by_site_scaler = {site_count, 1} ["" + selected_branch_total_length];
    column_vector = {state_count,1};
    
    sites_with_ambigs = {};
    
    averaged = {4,1}; // this hack enables deferred execution
    averaged [0] := (+pairwise_eps[-1][parent_state]$resolution)/resolution_count;
    averaged [1] := (+pairwise_epn[-1][parent_state]$resolution)/resolution_count;
    averaged [2] := (+pairwise_ops[-1][parent_state]$resolution)/resolution_count;
    averaged [3] := (+pairwise_opn[-1][parent_state]$resolution)/resolution_count;
    
    fully_resolved = {4,1}; // this hack enables deferred execution
    fully_resolved [0] := pairwise_eps[psi] * relative_branch_length;
    fully_resolved [1] := pairwise_epn[psi] * relative_branch_length;
    fully_resolved [2] := pairwise_ops[psi];
    fully_resolved [3] := pairwise_opn[psi];

      
    for (i = 0; i < selected_branches_count; i+=1) {
        this_branch            = selected_branches_in_avl[i];
        if (selected_branches_lengths[i]) {
            relative_branch_length = selected_branches_lengths[i] / selected_branch_total_length;
            parent_branch          = selected_branches_parents[i];
            
             for (s = 0; s < site_count; s += 1) {
                this_state      = matrix[this_branch][s];
                parent_state    = matrix[parent_branch][s];
                
                // parent state can only be resolved or --- (-1)
                
                if (this_state >= 0) { // child fully resolved (this means that the parent is fully resolved as well
                    psi = this_state*state_count+parent_state;
                    
                    for (k = 0; k < 4; k += 1) {
                        report_averaged[s*column_count + k] += fully_resolved[k];
                        report_resolved[s*column_count + k] += fully_resolved[k];
                    }
                     
                } else {
                    if (this_state == -1) { // the child is fully missing; no counts here
                        if (parent_state != -1) { // but the parent IS resolved
                            by_site_scaler [s] += (-selected_branches_lengths[i]);
                            psi = parent_state*state_count+parent_state;
                            for (k = 0; k < 2; k += 1) {
                                report_averaged[s*column_count + k] += fully_resolved[k];
                                report_resolved[s*column_count + k] += fully_resolved[k];
                            }
                           
                       }
                    } else { // the tip is an ambiguous, but partially resolved character
                             // this implies that the ancestor is fully resolved
                        resolution = lookup [-this_state-2]; // column vector with 1's for possible resolutions
                        resolution_count = + resolution;
                        
                        for (k = 0; k < 4; k += 1) {
                            report_averaged[s*column_count + k] = averaged[k];
                        }
                        
                        extract_site_info = sites_with_ambigs[s];
                        if (Type (extract_site_info) != "Matrix") {
                            extract_site_info    = matrix[{{0,s}}][{{tip_count-1,s}}];
                            extract_site_info    = extract_site_info[extract_site_info["_MATRIX_ELEMENT_VALUE_>=0"]];
                            if (Columns (extract_site_info) > 0) {
                               site_info = column_vector;
                               for (k = 0; k <  Columns (extract_site_info); k+=1) {
                                    site_info[extract_site_info[k]] += 1;
                               }
                               extract_site_info = site_info;
                            } 
                            sites_with_ambigs[s] = extract_site_info;
                            
                        }
                        
                        if (Columns (extract_site_info) > 0) {

                            fprintf (stdout, averaged, "\n", resolution, "\n");
                            resolution_filtered = extract_site_info $ resolution;
                            most_frequent_char = Max (extract_site_info $ resolution,1)[0];
                            if (most_frequent_char) {
                                resolution = resolution_filtered["_MATRIX_ELEMENT_VALUE_==most_frequent_char"];
                                resolution_count = + resolution;
                                //fprintf (stdout, averaged, "\n", resolution, "\nDONE\n");
                            }
                            
                        }
                        for (k = 0; k < 4; k += 1) {
                            report_resolved[s*column_count + k] += averaged[k];
                        }
                    }
                }
            }
        }
    }
    
    //fprintf (stdout, sites_with_ambigs, "\n", by_site_scaler, "\n");
    
    
    for (s = 0; s < site_count; s+=1) {
        k = by_site_scaler[s];
        report_resolved[s*column_count + 12] = k;
        report_averaged[s*column_count + 12] = k;
        
        if (k > 0) {
            k = selected_branch_total_length/k;
            report_resolved[s*column_count + 0] = report_resolved[s*column_count + 0] * k;
            report_resolved[s*column_count + 1] = report_resolved[s*column_count + 1] * k;
            report_averaged[s*column_count + 0] = report_averaged[s*column_count + 0] * k;
            report_averaged[s*column_count + 1] = report_averaged[s*column_count + 1] * k;
        }
    }
    
    #profile _hyphy_profile_dump;

    
/*
stats  			= _hyphy_profile_dump["STATS"];
_profile_summer = ({1,Rows(stats)}["1"]) * stats;
_instructions   = _hyphy_profile_dump["INSTRUCTION"];
_indices	    = _hyphy_profile_dump["INSTRUCTION INDEX"];

fprintf (stdout, "\nTotal run time (seconds)      : ", Format(_profile_summer[1],15,6),
                 "\nTotal number of steps         : ", Format(_profile_summer[0],15,0), "\n\n");
             
to_sort        =  stats["-_MATRIX_ELEMENT_VALUE_*_MATRIX_ELEMENT_COLUMN_+(_MATRIX_ELEMENT_COLUMN_==0)*_MATRIX_ELEMENT_ROW_"] % 1;  
             
for (k=0; k<Columns(_instructions); k=k+1)
{
    k2 = to_sort[k][0];
    fprintf (stdout, Format (_indices[k2],6,0), " : ", _instructions[k2], "\n\tCall count: ", stats[k2][0], 
                                                   "\n\tTime (seconds): ", stats[k2][1], "\n");
}
*/
    return report_averaged;

}

/*------------------------------------------------------------------------------------*/


io.reportProgressMessage("SLAC", "Performing joint maximum likelihood ancestral state reconstruction");

slac.counts    = genetic_code.ComputePairwiseDifferencesAndExpectedSites (slac.codon_data_info["code"], {"count-stop-codons" : FALSE});
slac.ancestors = ancestral.build (slac.global_mg_results["LF"], 0, None);

fprintf (stdout, slac.compute_the_counts (slac.ancestors["MATRIX"], slac.ancestors["TREE_AVL"], slac.ancestors["AMBIGS"], slac.selected_branches, slac.counts), "\n");

