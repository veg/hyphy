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

// namespace 'stats' for various descriptive stats functions 
LoadFunctionLibrary("libv3/stats.bf");

slac.timers = {
    6, 1
};

slac.samples = 100;

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
    "reference": "Not So Different After All: A Comparison of Methods for Detecting Amino Acid Sites Under Selection (2005). _Mol Biol Evol_ 22 (5): 1208-1222",
    "authors": "Sergei L Kosakovsky Pond and Simon DW Frost",
    "contact": "spond@temple.edu",
    "requirements": "in-frame codon alignment and a phylogenetic tree"
});





slac.codon_data_info = utility.promptForGeneticCodeAndAlignment("slac.codon_data", "slac.codon_filter");
slac.sample_size = slac.codon_data_info["sites"] * slac.codon_data_info["sequences"];
slac.codon_data_info["json"] = slac.codon_data_info["file"] + ".slac.json";


io.reportProgressMessage ("", ">Loaded an MSA with **" + slac.codon_data_info["sequences"] + "** sequences and **" + slac.codon_data_info["sites"] + "** codons from \`" + slac.codon_data_info["file"] + "\`");

slac.codon_frequencies = utility.defineFrequencies("slac.codon_filter");
slac.tree = utility.loadAnnotatedTopology(1);

slac.selected_branches = selection.io.defineBranchSets(slac.tree);
/* dict ["selected branch (original case)"] == TRUE */ 

slac.json["tested"] = slac.selected_branches;
slac.json["tree"] = slac.tree["string"];

io.reportProgressMessageMD("SLAC",  "selector", "Selected " + Abs(slac.selected_branches) + " branches to include in SLAC calculations " + Join(", ", Rows(slac.selected_branches)));

selection.io.startTimer(slac.timers, 1);

io.reportProgressMessageMD ("SLAC", "nuc-fit", "Obtaining branch lengths under the nucleotide GTR model");
slac.gtr_results = estimators.fitGTR("slac.codon_filter", slac.tree, None);
io.reportProgressMessageMD ("SLAC", "nuc-fit", "Log(L) = " + slac.gtr_results["LogL"]);
estimators.fixSubsetOfEstimates(slac.gtr_results, slac.gtr_results["global"]);

io.reportProgressMessageMD("SLAC", "codon-fit", "Obtaining the global omega estimate based on relative GTR branch lengths and nucleotide substitution biases");

parameters.declareGlobal(slac.scaler, None);
parameters.set_value(slac.scaler, 3);

slac.global_mg_results = estimators.fitMGREV(slac.codon_data_info, slac.tree, {
    "model-type": terms.global,
    "proportional-branch-length-scaler": {
        "0": slac.scaler
    },
    "retain-lf-object": TRUE
}, slac.gtr_results);
io.reportProgressMessageMD("SLAC", "codon-fit", "Log(L) = " + Format(slac.global_mg_results["LogL"],8,2));

slac.global_dnds = selection.io.extract_global_MLE(slac.global_mg_results, terms.omega_ratio);
io.reportProgressMessageMD ("SLAC", "codon-fit", "Global dN/dS = " + Format (slac.global_dnds,8,2));

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
io.reportProgressMessageMD("SLAC", "anc", "Performing joint maximum likelihood ancestral state reconstruction");

slac.counts    = genetic_code.ComputePairwiseDifferencesAndExpectedSites (slac.codon_data_info["code"], {"count-stop-codons" : FALSE});
slac.ancestors = ancestral.build (slac.global_mg_results["LF"], 0, None);
slac.results   = slac.compute_the_counts (slac.ancestors["MATRIX"], slac.ancestors["TREE_AVL"], slac.ancestors["AMBIGS"], slac.selected_branches, slac.counts);

slac.json ["MLE"] = {"headers" : {{"ES", "Expected synonymous sites"}
                                  {"EN", "Expected non-synonymous sites"}
                                  {"S", "Inferred synonymous substitutions"}
                                  {"N", "Inferred non-synonymous substitutions"}
                                  {"P<sub>S</sub>", "Expected proportion of synonymous sites"}
                                  {"dS", "Inferred synonymous susbsitution rate"}
                                  {"dN", "Inferred non-synonymous susbsitution rate"}
                                  {"dN-dS", "Scaled by the length of the tested branches"}
                                  {"Positive selection", "Binomial probability that S is no greater than the observed value, with P<sub>s</sub> probability of success"}
                                  {"Negative selection", "Binomial probability that S is no less than the observed value, with P<sub>s</sub> probability of success"}
                                  {"Effective branch length", "The total length of branches contributing to inference at this site, and used to scale dN-dS"}},
                     "data" : slac.results };


io.spool_json (slac.json, slac.codon_data_info["json"]);

io.reportProgressMessageMD ("SLAC", "anc", "Generating `slac.samples` ancestral sequence samples to obtain confidence intervals");

slac.sample.results = {};

for (slac.s = 0; slac.s < slac.samples; slac.s+=1) {
    slac.sampled   = ancestral.build (slac.global_mg_results["LF"], 0, {"sample": TRUE});
    slac.sample.results + slac.compute_the_counts (slac.sampled["MATRIX"], slac.sampled["TREE_AVL"], slac.sampled["AMBIGS"], slac.selected_branches, slac.counts);
    io.reportProgressBar("", "\tSample " + (slac.s+1) + "/" + slac.samples);
}

slac.extractor = {slac.samples, 1};

slac.sites = slac.codon_data_info["sites"];
slac.columns = Columns (slac.results["RESOLVED"]);


slac.json ["sample-median"] = {"RESOLVED" : {slac.sites, slac.columns}, "AVERAGED" : {slac.sites, slac.columns}};
slac.json ["sample-2.5"]    = {"RESOLVED" : {slac.sites, slac.columns}, "AVERAGED" : {slac.sites, slac.columns}};
slac.json ["sample-97.5"]   = {"RESOLVED" : {slac.sites, slac.columns}, "AVERAGED" : {slac.sites, slac.columns}};
slac.keys = {{"RESOLVED", "AVERAGED"}};


for (slac.s = 0; slac.s < slac.sites; slac.s += 1) {
    for (slac.c = 0; slac.c <slac.columns; slac.c+=1) {
        for (slac.key = 0; slac.key < Columns(slac.keys); slac.key += 1) {
        
            slac.key_value = slac.keys[slac.key];
            
            slac.col = slac._extract_vector (slac.extractor, slac.sample.results["by-site"], slac.key_value, slac.s,slac.c);
            ((slac.json["sample-median"])[slac.key_value])[slac.s][slac.c] = stats.quantile (slac.col, 0.5);
            ((slac.json["sample-2.5"])[slac.key_value])[slac.s][slac.c] = stats.quantile (slac.col, 0.025);
            ((slac.json["sample-97.5"])[slac.key_value])[slac.s][slac.c] = stats.quantile (slac.col, 0.975);
        }
    }
}

io.spool_json (slac.json, slac.codon_data_info["json"]);


/*___________________________________________________________________________________________________________*/
// HELPER FUNCTIONS
/*___________________________________________________________________________________________________________*/

lfunction slac.extendedBinTail (ebn, ebp, ebx) {
/* 
	returns the LEFT-tail probability for the extended binomial distribution 
	with ebn observations, ebp probability of success and ebx successes 
	i.e, Pr (X <= ebx | enb, ebp)

*/
	if (ebp == 0) {
		return 0;	
	}

	ebr = ebx$1; /* rounded to nearest integer */
	
	currentBinCoeff = (1-ebp)^ebn; /*compute the first binomial coefficient */
	
	binHead = 0;
	
	for (ebk=0; ebk<=ebr; ebk += 1) {
		binHead			+= currentBinCoeff;
		currentBinCoeff = currentBinCoeff * (ebn-ebk) / (ebk+1) * ebp / (1-ebp);
	}
	
	if (ebx <= ebn$1) {
		binHead += currentBinCoeff*(ebx-ebr);
	}
	else {
		binHead += (1-binHead)*(ebx-ebr)/(ebn-ebn$1);	
	}
		
	return binHead;
}


/*------------------------------------------------------------------------------------*/

function slac._extract_vector (to, from, key, row, column) {
    return to ["((from[_MATRIX_ELEMENT_ROW_])[key])[row][column]"]%0;
}

/*------------------------------------------------------------------------------------*/

lfunction slac.compute_the_counts (matrix, tree, lookup, selected_branches, counts) {
//#profile START;

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
           4 Expected S/N ratio
           5 dS
           6 dN
           7 dN-dS (scaled)
           8 P {S <= observed} // posititve sel
           9 P {S >= observed} // negative sel 
           10 effective branch length
           
    */
    
    column_count    = 11;
    
    report_resolved = {site_count, column_count}["0"];
    report_averaged = {site_count, column_count}["0"];
    
    report_resolved_by_branch = {selected_branches_count, column_count}["0"];
    report_averaged_by_branch = {selected_branches_count, column_count}["0"];
    report_branch_names = {selected_branches_count,1};
    
    pairwise_eps = counts ["EPS"];
    state_count  = Rows (pairwise_eps);
    pairwise_epn = counts ["EPN"];
    pairwise_ops = counts ["OPS"];
    pairwise_opn = counts ["OPN"];
    by_site_scaler = {site_count, 1} ["" + selected_branch_total_length];
    column_vector = {state_count,1};
    
    sites_with_ambigs = {};
    
    averaged = {4,1}; // this hack enables deferred execution
    averaged [0] := (+pairwise_eps[-1][parent_state]$resolution)/resolution_count * relative_branch_length;
    averaged [1] := (+pairwise_epn[-1][parent_state]$resolution)/resolution_count * relative_branch_length;
    averaged [2] := (+pairwise_ops[-1][parent_state]$resolution)/resolution_count;
    averaged [3] := (+pairwise_opn[-1][parent_state]$resolution)/resolution_count;
    
    fully_resolved = {4,1}; // this hack enables deferred execution
    fully_resolved [0] := pairwise_eps[psi] * relative_branch_length;
    fully_resolved [1] := pairwise_epn[psi] * relative_branch_length;
    fully_resolved [2] := pairwise_ops[psi];
    fully_resolved [3] := pairwise_opn[psi];
 
    for (i = 0; i < selected_branches_count; i+=1) {
        this_branch            = selected_branches_in_avl[i];
        report_branch_names [i] = (tree[this_branch+1])["Name"];
        if (selected_branches_lengths[i]) {
            
            relative_branch_length = selected_branches_lengths[i] / selected_branch_total_length;
            parent_branch          = selected_branches_parents[i];
            
             for (s = 0; s < site_count; s += 1) {
                this_state      = matrix[this_branch][s];
                parent_state    = matrix[parent_branch][s];
                
                // parent state can only be resolved or --- (-1)
                
                if (this_state >= 0) { // child fully resolved (this means that the parent is fully resolved as well
                    psi = this_state*state_count+parent_state;
                    
                    fr = Eval (fully_resolved);
                    for (k = 0; k < 4; k += 1) {
                        report_averaged[s*column_count + k] += fr[k];
                        report_resolved[s*column_count + k] += fr[k];
                    }
                    
                    relative_branch_length = 1;
                    fr = Eval (fully_resolved);
                    for (k = 0; k < 4; k += 1) {
                        report_resolved_by_branch[i*column_count + k] += fr[k];
                        report_averaged_by_branch[i*column_count + k] += fr[k];
                    }
                     
                } else {
                    if (this_state == -1) { // the child is fully missing; no counts here
                        if (parent_state != -1) { // but the parent IS resolved
                            by_site_scaler [s] += (-selected_branches_lengths[i]);
                            psi = parent_state*state_count+parent_state;
                            
                            
                            fr = Eval (fully_resolved);
                            for (k = 0; k < 2; k += 1) {
                                report_averaged[s*column_count + k] += fr[k];
                                report_resolved[s*column_count + k] += fr[k];
                            }

                            relative_branch_length = 1;
                            fr = Eval (fully_resolved);
                            for (k = 0; k < 2; k += 1) {
                                report_resolved_by_branch[i*column_count + k] += fr[k];
                                report_averaged_by_branch[i*column_count + k] += fr[k];
                            }
                           
                       }
                    } else { // the tip is an ambiguous, but partially resolved character
                             // this implies that the ancestor is fully resolved
                        resolution = lookup [-this_state-2]; // column vector with 1's for possible resolutions
                        resolution_count = + resolution;
                        
                        av = Eval (averaged);
                        for (k = 0; k < 4; k += 1) {
                            report_averaged[s*column_count + k] += av[k];
                            report_averaged_by_branch [i*column_count + k] += av[k];
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

                            resolution_filtered = extract_site_info $ resolution;
                            most_frequent_char = Max (extract_site_info $ resolution,1)[0];
                            if (most_frequent_char) {
                                resolution = resolution_filtered["_MATRIX_ELEMENT_VALUE_==most_frequent_char"];
                                resolution_count = + resolution;
                           }
                            
                        }
                        
                       av = Eval (averaged);
                        for (k = 0; k < 4; k += 1) {
                            report_resolved[s*column_count + k] += av[k];
                            report_resolved_by_branch [i*column_count + k] += av[k];
                        }
                        
                    }
                }
            }
        }
    }
    
    
    receivers = {"0" : "report_resolved", "1": "report_averaged", "2" : "report_resolved_by_branch", "3": "report_averaged_by_branch"};
    for (i = 0; i < Abs (receivers); i+=1) {
        mx = receivers[i];
        upto = Rows (*mx);

        for (s = 0; s < upto; s+=1) {
            k = by_site_scaler[s];
            (*mx)[s*column_count + column_count - 1] = k;
            (*mx)[s*column_count + column_count - 1] = k;
        
            if (k > 0) {
                sc = selected_branch_total_length/k;
                (*mx)[s*column_count + 0] = (*mx)[s*column_count + 0] * sc;
                (*mx)[s*column_count + 1] = (*mx)[s*column_count + 1] * sc;
            }
    
            total_subs = (*mx)[s*column_count + 2] + (*mx)[s*column_count + 3];
            
            if (total_subs) {
                (*mx) [s*column_count + 4] = (*mx) [s*column_count + 0]/((*mx) [s*column_count + 1]+(*mx) [s*column_count + 0]);
                (*mx) [s*column_count + 5] = (*mx) [s*column_count + 2]/(*mx) [s*column_count + 0];
                (*mx) [s*column_count + 6] = (*mx) [s*column_count + 3]/(*mx) [s*column_count + 1];
            
                if (k > 0) {
                    (*mx) [s*column_count + 7] = ((*mx) [s*column_count + 6] - (*mx) [s*column_count + 5])/k;
                }
                
                (*mx) [s*column_count + 8] = slac.extendedBinTail(total_subs,(*mx) [s*column_count + 4],(*mx)[s*column_count + 2]);
                (*mx) [s*column_count + 9] = 1-slac.extendedBinTail(total_subs,(*mx) [s*column_count + 4],Max (0, (*mx)[s*column_count + 2]-1));
            }
            
        }
    }
    
    /*#profile _hyphy_profile_dump;

    

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
}*/



    return {"by-site" :{"AVERAGED" : report_averaged,
                        "RESOLVED" : report_resolved},
                        "by-branch" :{"NAMES": report_branch_names, "AVERAGED" : report_averaged_by_branch,
                        "RESOLVED" : report_resolved_by_branch}
                        };

}

