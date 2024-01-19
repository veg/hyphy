RequireVersion("2.5.58");

/*------------------------------------------------------------------------------
    Load library files
*/

LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");
LoadFunctionLibrary("libv3/all-terms.bf");

LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");

LoadFunctionLibrary("libv3/convenience/math.bf");


/*------------------------------------------------------------------------------ 
    Display analysis information
*/

efilter.analysis_description = {
    terms.io.info: "
The error filter analysis reads a BUSTED-E JSON result file,
identifies sites which may be due to alignment or other error, 
and masks them according to a user-selected strategy, saving the resulting MSA to a new file. 
A machine readable report on the errors filtered (if any)  is also created. 
    ",
    terms.io.version: "0.1",
    terms.io.reference: "TBD",
    terms.io.authors: "Sergei L Kosakovsky Pond",
    terms.io.contact: "spond@temple.edu",
    terms.io.requirements: "A BUSTED-E .json results file (possibly .gz compressed)",
    terms.settings: {}
};
io.DisplayAnalysisBanner(efilter.analysis_description);

KeywordArgument ("json", "Load the BUSTED-E .json result file", None);
efilter.path = io.PromptUserForFilePathRead ("Load the BUSTED-E .json result file");
io.ReportProgressMessageMD ("efilter","loading","Loading .JSON data");

KeywordArgument ("threshold", "Empirical Bayes Factor error threshold for masking sites", 100.);
efilter.threshold = io.PromptUser ("Empirical Bayes Factor threshold for masking sites", 100., 1, 1e100, FALSE);

KeywordArgument ("ratio", "Empirical Bayes Factor for error vs selection ", 20.);
efilter.ratio_threshold = io.PromptUser ("Empirical Bayes Factor for error vs selection", 20., 1, 1e100, FALSE);


KeywordArgument ("site-threshold", "Mask the entire site if more than X% sequences are problematic", .4);
efilter.site_threshold = io.PromptUser ("Mask the entire site if more than X% sequences are problematic", 0.4, 0., 1., FALSE);


efilter.json_out = {
   terms.json.analysis : efilter.analysis_description,
   terms.settings : {
        terms.empirical_bayes_factor : efilter.threshold,
        "BF ratio" : efilter.ratio_threshold,
        "site threshold" : efilter.site_threshold 
   }
};

efilter.json =  io.ParseJSON(efilter.path);

io.CheckAssertion ("utility.GetByKey (efilter.json, {{terms.json.analysis, 'settings', 'error-sink'}},'Number')", "There is no error sink data in the JSON file");


efilter.error_sink_proportion = utility.GetByKey (efilter.json, {{terms.json.fits,"Unconstrained model",terms.json.rate_distribution,"Test","0",terms.json.proportion}},'Number');

io.CheckAssertion ("None!=efilter.error_sink_proportion", "There is no estimate of the error-sink omega proportion");

efilter.rates = utility.GetByKey (efilter.json, {{terms.json.fits,"Unconstrained model",terms.json.rate_distribution,"Test"}},'AssociativeList');
efilter.rates_count = utility.Array1D (efilter.rates);

efilter.fast_omega_proportion = utility.GetByKey (efilter.json, {{terms.json.fits,"Unconstrained model",terms.json.rate_distribution,"Test","" + (efilter.rates_count-1),terms.json.proportion}},'Number');

efilter.input = utility.GetByKey (efilter.json, {{terms.json.input}}, "AssociativeList");
io.CheckAssertion ("None!=efilter.input", "There is no information about the alignment");

io.ReportProgressMessageMD ("efilter","loading","- Original data file path: \`" + efilter.input[terms.json.file] + "\`");
io.ReportProgressMessageMD ("efilter","loading","- Partitions: " + efilter.input[terms.json.partition_count]);
io.ReportProgressMessageMD ("efilter","loading","- Sequences: " + efilter.input[terms.json.sequences]);
io.ReportProgressMessageMD ("efilter","loading","- Sites: " + efilter.input[terms.json.sites]);

efilter.json_out[terms.json.input] = efilter.input;

io.ReportProgressMessageMD ("efilter","loading","- Estimated proportion of the MSA in the error-sink class = " + Format (efilter.error_sink_proportion *100, 8, 3) + "%");

if (efilter.error_sink_proportion == 0) {
    efilter.error_sink_prior_odds = 1e100;
} else {
    efilter.error_sink_prior_odds = efilter.error_sink_proportion / (1-efilter.error_sink_proportion);
}
efilter.error_sink_prior_odds_ratio =  Min (1e25,efilter.error_sink_proportion / efilter.fast_omega_proportion);

efilter.verbose_logging = utility.GetEnvVariable ("VERBOSE_LOGGING");

efilter.site_offset = 0;

utility.ToggleEnvVariable ("LOOP_TREE_ITERATOR_PREORDER", TRUE);



for (p = 0; p < efilter.input[terms.json.partition_count]; p+=1) {
    efilter.sequences    = {};
    efilter.masked_sites = {};
    efilter.branch_data =  utility.GetByKey (efilter.json, {{terms.json.branch_attributes,p}}, "AssociativeList");
    efilter.tested_branches = utility.GetByKey (efilter.json, {{terms.json.tested,p}}, "AssociativeList");
    efilter.tree = utility.GetByKey (efilter.input, {{terms.json.trees,p}}, "String");
    efilter.subs = utility.GetByKey (efilter.json, {{terms.substitutions,p}}, "AssociativeList");
    Topology efilter.T = efilter.tree;
    efilter.sites = utility.Array1D (efilter.subs );
    
    if (p == 0) {
        efilter.leaves = {};
        for (s;in;TipName (efilter.T,-1)) {
            efilter.sequences[s] = "";
            efilter.sequences[s] * efilter.sites;
            efilter.masked_sites[s] = {};  
            efilter.leaves[s] = 1;      
        }
        
    }    
    efilter.descendants = {};
    efilter.leaf_descendants = {};
    
    for (s;in;BranchName (efilter.T,-1)) {
        if (efilter.sequences / s) {
            continue;
        }
        efilter.descendants[s] = {};
        efilter.leaf_descendants[s] = {};
        for (s2, v; in; efilter.T [s]) {
            if (v == 1) {
                (efilter.leaf_descendants[s])[s2] = 1;
            }
            (efilter.descendants[s])[s2] = 1;
        } 
        if (Abs (efilter.leaf_descendants[s]) * 2 > Abs (efilter.leaves)) {
            t = efilter.leaf_descendants[s];
            efilter.leaf_descendants[s] = efilter.leaves;
            efilter.leaf_descendants[s] - t;
        }
    }
    
    
    for (site = 0; site < efilter.sites; site += 1) {
        efilter.subs4site = efilter.subs [site];
        efilter.states = {};
        efilter.masked_already = {};
        efilter.write_out = {};
        for (parent, node; in; efilter.T) {
            if (parent) {
                if (efilter.subs4site/node) {
                    efilter.states [node]  = efilter.subs4site[node];
                } else {
                    efilter.states [node] = efilter.states [parent]
                }
                efilter.write_state = efilter.states [node] ;
                if (efilter.masked_already [node]) {
                    continue;
                }
                if (efilter.branch_data / node) {
                    efilter.site_prob = ((efilter.branch_data[node])["Posterior prob omega class by site"])[0][site];
                    if (None ==  efilter.site_prob) {
                        continue;
                    }
                    efilter.site_BF2 = ((efilter.branch_data[node])["Posterior prob omega class by site"])[efilter.rates_count-1][site];
                    if (efilter.site_prob < 1) {
                        efilter.site_BF = efilter.site_prob/(1-efilter.site_prob)/efilter.error_sink_prior_odds;
                    } else {
                        efilter.site_BF = 1e25;
                    }
                    
                    if (efilter.site_BF2 < 1) {
                        efilter.site_BF2 = efilter.site_prob / efilter.site_BF2 / efilter.error_sink_prior_odds_ratio;
                    } else {
                        efilter.site_BF2 = 1e25;
                    }
                    
                    if (efilter.verbose_logging ) {
                        if (efilter.site_BF >= efilter.threshold || efilter.site_BF2 >= efilter.ratio_threshold) {
                            if (efilter.site_BF >= efilter.threshold && efilter.site_BF2 >= efilter.ratio_threshold) {
                                fprintf (stdout, "FILTERED\n");
                            } else {
                                fprintf (stdout, "KEPT\n");
                            }
                            fprintf (stdout, ">" + (site+1) + "/" + node, " " + efilter.site_BF + ":" +efilter.site_BF2 + "\n");
                        
                        }
                    }
                    
                    
                    if (efilter.site_BF >= efilter.threshold && efilter.site_BF2 >= efilter.ratio_threshold) {
                        if (efilter.masked_sites/node) { // terminal node
                            if (efilter.masked_already [ntm] != 1) {
                                efilter.masked_sites [node] + (site + efilter.site_offset);
                            }
                            efilter.write_out[node] = "---";
                            efilter.masked_already [node] = 1;
                        } else {
                            if (Abs (efilter.leaf_descendants[node]) / Abs (efilter.sequences) >= efilter.site_threshold) {
                                for (ntm, ignore; in; efilter.sequences) {
                                    efilter.write_out[ntm] = "---";
                                    if (efilter.masked_already [ntm] != 1) {
                                        efilter.masked_sites [ntm] + (site + efilter.site_offset);
                                    }  
                                }   
                                if (efilter.verbose_logging) {
                                    console.log ("Masking everything " + site + " " + node + " " + Abs (efilter.leaf_descendants) / Abs (efilter.sequences));               
                                }         
                                break;
                            }
                            for (ntm, ignore; in; efilter.leaf_descendants[node]) {
                                efilter.write_out[ntm] = "---";
                                //console.log ("Masking child " + ntm + " of parent " + node + " at site " + (1+site));   
                                if (efilter.masked_already [ntm] != 1) {
                                    efilter.masked_sites [ntm] + (site + efilter.site_offset);
                                }
                            }
                            efilter.masked_already * efilter.leaf_descendants[node];
                        }
                    }
                }
            } else {
                efilter.states [node] = efilter.subs4site["root"];
            }
            
            
            if (efilter.sequences / node && efilter.masked_already[node] == 0) {
                efilter.write_out[node] = efilter.write_state;
            }
        }
        

        for (s,st; in; efilter.write_out) {
             efilter.sequences[s] * st;
        }
        
    }
    
    efilter.site_offset += efilter.sites;    
}

utility.ToggleEnvVariable ("LOOP_TREE_ITERATOR_PREORDER", None);

KeywordArgument ("output", "Write a filtered MSA to");
efilter.path = io.PromptUserForFilePath ("Write a filtered MSA to");

io.ReportProgressMessageMD ("efilter","results","Filtering results");

efilter.settings = {"header" : TRUE, "column-widths": {
            "0": 50,
            "1": 20,
            "2": 20
            }};

fprintf (stdout, "\n", io.FormatTableRow ({{"Sequence", "Filtered Sites", "Filtered Proportion, %"}}, efilter.settings));
efilter.settings ["header"] = FALSE;
efilter.total_filtered = 0;

for (s;in;TipName (efilter.T,-1)) {
    efilter.row_matrix = {1,3};

    efilter.row_matrix [0] =  s;
    efilter.row_matrix [1] =  utility.Array1D (efilter.masked_sites[s]);
    efilter.row_matrix [2] =  Format (efilter.row_matrix [1] / efilter.site_offset*10, 10, 3);
    
    efilter.total_filtered += efilter.row_matrix [1];
    
    fprintf (stdout, io.FormatTableRow (efilter.row_matrix, efilter.settings));
    efilter.sequences[s] * 0;
    fprintf (efilter.path, ">", s, "\n", efilter.sequences[s], "\n");
}

if (efilter.input[terms.json.partition_count] == 1) {
    fprintf (efilter.path, "\n", Format (efilter.T,1,0));
}

io.ReportProgressMessageMD ("efilter","results","\nMasked a total of **`efilter.total_filtered`** or `Format (efilter.total_filtered/efilter.input[terms.json.sequences]/efilter.input[terms.json.sites]*100,8,3)`% sites");


fprintf (stdout, "\n");

KeywordArgument ("output-json", "Write the resulting JSON to this file");
efilter.json_path = io.PromptUserForFilePath ("Save the resulting JSON file to");

efilter.json_out ["filter"] = efilter.masked_sites ;

io.SpoolJSON (efilter.json_out,efilter.json_path); 





