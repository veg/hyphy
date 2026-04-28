RequireVersion ("2.5.98");

LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");

LoadFunctionLibrary("libv3/tasks/trees.bf");

CS.analysis_description = {
    terms.io.info: 
"This helper analysis processes a BUSTED-PH JSON file to compute post-hot Effective Clade Breadth (ECB() support for selection signal.
The Effective Clade Breadth metric quantifies how many independent phenotypic origins contribute to the overall evidence for positive selection. It calculates the average expected number of selective sites per branch (using empirical Bayes posteriors) for each clade and applies an information-theoretic perplexity measure to the resulting distribution. This identifies whether a gene-phenotype association is driven by a single dominant lineage or represents a broad signal replicated across multiple independent transitions.
    ",
    terms.io.version: "0.0.1",
    terms.io.reference: "",
    terms.io.authors: "Sergei L Kosakovsky Pond",
    terms.io.contact: "spond@temple.edu",
    terms.io.requirements: "a BUSTED-PH JSON results file",
    terms.settings: {}
};

io.DisplayAnalysisBanner (CS.analysis_description);


KeywordArgument ("json", "BUSTED-PH result file (JSON)");

CS.filepath = io.PromptUserForFilePathRead ("BUSTED-PH result file (JSON)");
CS.busted_ph = io.ParseJSON (CS.filepath);

KeywordArgument ("output", "Write the resulting JSON to this file");

CS.output = io.PromptUserForFilePath ("Save the resulting JSON file to");


CS.has_error_sink             = utility.truthy (CS.busted_ph${{"analysis","settings","error-sink"}});
CS.rates = CS.busted_ph${{"fits","Unconstrained model","Rate Distributions","Test"}};

CS.rate_count = utility.Array1D (CS.rates);
CS.positive_rates = {1,CS.rate_count};

for (i,v;in;CS.rates) {
    rc = +i;
    if (v["omega"] > 1) {
        if (rc != 0 || !CS.has_error_sink) {
            CS.positive_rates [rc] = 1;
        }
    }
}

CS.results = {};

CS.table_output_options = {terms.table_options.header : TRUE, 
                            terms.table_options.minimum_column_width: 12, 
                            terms.table_options.align : "left"};

CS.table_header = {"0" : "Clade Root", "1" : "Weight", "2" : "Exp sites/br", "3" : "Branches", "4" : "Tips"};


for (partition, info; in; CS.busted_ph [terms.json.branch_attributes]) {
    if (partition != terms.json.attribute) {
        CS.results[partition] = {};
        CS.branch_posteriors = {};
        CS.tree = CS.busted_ph${{terms.json.input, terms.json.trees,partition}};
        
        Topology CS.T = CS.tree;
        
        // compute by branch posterior expectations
        
        CS.lengths = {};

        for (b,binfo; in; info) {
            CS.bp = binfo["Posterior prob omega class by site"];
            CS.lengths [b] = binfo ["MG94xREV with separate rates for branch sets"];
            if (Type (CS.bp) == "Matrix") {
                (CS.branch_posteriors)[b] = +(CS.positive_rates * CS.bp );
            }
        }
        
        /*   
             identify foreground "clades"
             a clade is a subtree (with possibly one branch)
             where every single branch is labeled as foreground
             this is computed by traversing the tree in post-order
        */
        
        CS.parent_map  = trees.ParentMap ("CS.T");
        CS.parent_tags = {};
        
        for (branch;in;CS.T) {
            CS.parent_tags[branch] = TRUE;
        }
    
        for (branch;in;CS.T) {
            if ((CS.branch_posteriors) / branch) {
                CS.parent_tags[branch] = TRUE;
            } else {
                CS.parent_tags[branch] = FALSE;
                CS.parent_tags[CS.parent_tags[branch]] = FALSE;
            }
        }
        
        CS.clades = {};
        CS.clade_stats = {};
        
        (CS.results[partition])['clades'] = {};
        
                
        for (branch;in;CS.T) {
            if (CS.parent_tags[branch] && CS.parent_tags[CS.parent_map [branch]] == FALSE) {
                
                B_N = 0;
                T_N = 0;
                CS.clades [branch] = (CS.branch_posteriors)[branch];
                
                ((CS.results[partition])['clades'])[branch] = CS.T[branch];
                
                CS.subtree = CS.T[branch];
                for (br, depth; in; CS.subtree) {
                    B_N += 1;
                    if (depth == 1) {
                        T_N += 1;
                    }
                    CS.clades [branch] += (CS.branch_posteriors)[br];
                }
                
                CS.clade_stats [branch] = {"branches" : B_N, "tips" : T_N};
                CS.clades [branch] = CS.clades [branch]/B_N;
            }
        }

        (CS.results[partition])['expected_sites'] = {};
        for (clade, weight; in; CS.clades) {
            ((CS.results[partition])['expected_sites'])[clade] = weight;
        }
        
        (CS.results[partition])['clade_stats'] = CS.clade_stats;
        (CS.results[partition])['weights'] = {};
          
        CS.total         = + CS.clades;
        CS.entropy = 0;
        
        io.ReportProgressMessageMD ("CS", "partition", "Support by independent foreground clade for partition " + partition);
        fprintf (stdout, io.FormatTableRow (CS.table_header, CS.table_output_options));
        CS.table_output_options[terms.table_options.header] = FALSE;

        for (clade, weight; in; CS.clades) {
            CS.weight = weight / CS.total;
            ((CS.results[partition])['weights'])[clade] = CS.weight;
            CS.entropy += Log (CS.weight) * CS.weight;
            
            CS.stats = CS.clade_stats[clade];
            
            CS.table_row = {5, 1};
            CS.table_row [0] = clade;
            CS.table_row [1] = Format(CS.weight, 6, 4);
            CS.table_row [2] = Format(weight, 6, 4);
            CS.table_row [3] = "" + CS.stats["branches"];
            CS.table_row [4] = "" + CS.stats["tips"];
            fprintf (stdout, io.FormatTableRow (CS.table_row, CS.table_output_options));
        }
        
        CS.table_output_options[terms.table_options.header] = TRUE;

        CS.perplexity = Exp(-CS.entropy);
        (CS.results[partition])['perplexity'] = CS.perplexity;
        
        
        function bl (n) {
            return "" + CS.lengths[n["Name"]];
        }
        
        (CS.results[partition])['tree'] = tree.Annotate ("CS.T", {}, '{}', "bl");
        (CS.results[partition])['branch_support'] = CS.branch_posteriors;
        
        fprintf (stdout, "\n**Effective Clade Breadth (Perplexity)**: " + Format(CS.perplexity, 6, 4) + "\n");
    }
}

io.SpoolJSON (CS.results, CS.output);
