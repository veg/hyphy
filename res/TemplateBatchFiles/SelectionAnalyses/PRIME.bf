RequireVersion("2.5.37");

/*------------------------------------------------------------------------------
    Load library files
*/

LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");

LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");

LoadFunctionLibrary("libv3/models/codon/MG_REV_PROPERTIES.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");

LoadFunctionLibrary("modules/io_functions.ibf");
LoadFunctionLibrary("modules/selection_lib.ibf");



debug.checkpoint = utility.GetEnvVariable ("DEBUG_CHECKPOINT");
/*------------------------------------------------------------------------------ 
    Display analysis information
*/

prime.analysis_description = {
    terms.io.info: "PRIME (PRoperty Informed Model of Evolution). Given a set of N amino-acid properties,
    fit a site-level model where non-synonymous rates depend on how much a non-synonymous substitution changes the 
    properties of the residue, beta (X,Y) = Exp (log_omega - lambda_1 * diff_1 (X,Y )- lambda_2 * diff_2 (X,Y) -...).
    When lambda_k > 0, changes in property k are disfavored and when lambda_k < 0 -- they are promoted.
    At each site, N tests are performed
    A subset of branches can be selected
    for testing as well, in which case an additional (nuisance) parameter will be
    inferred -- the non-synonymous rate on branches NOT selected for testing. Multiple partitions within a NEXUS file are also supported
    for recombination - aware analysis.
    ",
    terms.io.version: "0.2.0",
    terms.io.reference: "TBD",
    terms.io.authors: "Sergei L. Kosakovsky Pond",
    terms.io.contact: "spond@temple.edu",
    terms.io.requirements: "in-frame codon alignment and a phylogenetic tree"
};


io.DisplayAnalysisBanner(prime.analysis_description);


/*------------------------------------------------------------------------------
    Environment setup
*/

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);


/*------------------------------------------------------------------------------
    Globals
*/

prime.parameter_site_alpha          = "Site relative synonymous rate";
prime.parameter_site_beta_nuisance  = "Site relative non-synonymous rate (untested branches)";
prime.parameter_property_prefix     = "Site importance for property ";
prime.parameter_site_beta           = "Site log (baseline beta)";
terms.prime_imputed_states          = "Imputed States";

// default cutoff for printing to screen
prime.pvalue = 0.1;
// The dictionary of results to be written to JSON at the end of the run
prime.json = {
    terms.json.analysis: prime.analysis_description,
    terms.json.fits: {},
    terms.json.timers: {}
};

prime.display_orders =   {  
                            terms.original_name: -1,
                            terms.json.nucleotide_gtr: 0,
                            terms.json.global_mg94xrev: 1
                         };

selection.io.startTimer (prime.json [terms.json.timers], "Total time", 0);

/*------------------------------------------------------------------------------
    Key word arguments
*/

KeywordArgument ("code", "Which genetic code should be used", "Universal");
KeywordArgument ("alignment", "An in-frame codon alignment in one of the formats supported by HyPhy");
KeywordArgument ("tree", "A phylogenetic tree (optionally annotated with {})", null, "Please select a tree file for the data:");
KeywordArgument ("branches",  "Branches to test", "All");
KeywordArgument ("pvalue",  "The p-value threshold to use when testing for selection", "0.1");
KeywordArgument ("properties", "Which set of properties to use", "Atchley");

prime.scaler_prefix = "prime.scaler";
 
/**
This table is meant for HTML rendering in the results web-app; can use HTML characters, the second column
is 'pop-over' explanation of terms. This is ONLY saved to the JSON file. For Markdown screen output see
the next set of variables.
*/

prime.site.composition.string = "";


prime.table_screen_output  = {{"Codon", "Partition", "#subs", "#AA", "alpha", "beta", "Property Description", "Importance", "p-value","Most common codon substitutions at this site"}};
prime.table_output_options = {  terms.table_options.header : TRUE, 
                                terms.table_options.minimum_column_width: 12, 
                                terms.table_options.align : "center",
                                terms.table_options.column_widths: {
                                    "0" : 12,
                                    "1" : 12,
                                    "2" : 8,
                                    "3" : 8,
                                    "4" : 12,
                                    "5" : 12,
                                    "6" : 40,
                                    "7" : 12,
                                    "8" : 12,
                                    "9" : 70}
                             };


namespace prime {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ("prime");
}

prime.pvalue  = io.PromptUser ("\n>Select the p-value threshold to use when testing for selection",prime.pvalue,0,1,FALSE);

prime.property_set = io.SelectAnOption (
            {
                "Atchley":"Use the five properties derived from a factor analysis of 500 amino-acid properties [Table 2 in PNAS (2005) 102(18) 6395-6400 doi: 10.1073/pnas.0408677102]",
                "LCAP":"Use the five properties defined in the Conant and Stadler LCAP model [Mol Biol Evol (2009) 26 (5): 1155-1161. doi: 10.1093/molbev/msp031]",
                "Random" : "Random properties (for null hypothesis testing)",
                "Custom":"Load the set of properties from a file"
            }, 
            "The set of properties to use in the model");
           
(prime.json [terms.json.analysis])[terms.model.residue_properties] = prime.property_set;
            
if (prime.property_set == "Custom") {
    KeywordArgument ("property-file", "JSON file which defines amino-acid properties");
    prime.property_set = io.PromptUserForFilePathRead ("JSON file which defines amino-acid properties");
    prime.property_set = io.ParseJSON(prime.property_set);
    console.log (">Loaded a set of `Abs(prime.property_set)` properties");
}
     

KeywordArgument ("impute-states", "Use site-level model fits to impute likely character states for each sequence", "No");

prime.impute_states = io.SelectAnOption (
            {
                "Yes":"Impute marginal likelihoods for each codon at each sequence and each site",
                "No": "Do not impute marginal likelihoods for each codon at each sequence and each site",
            }, 
            "Impute likely states for sequences") == "Yes";


KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'prime.json')", prime.codon_data_info [terms.json.json]);
prime.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");



prime.pairwise_counts = genetic_code.ComputePairwiseDifferencesAndExpectedSites(prime.codon_data_info[terms.code], None);
prime.codon_table = genetic_code.DefineCodonToAAMapping (prime.codon_data_info[terms.code]);

selection.io.startTimer (prime.json [terms.json.timers], "Model fitting",1);

if (Type (debug.checkpoint) != "String") {
    namespace_tag = "prime";
    namespace prime {
        doGTR ("prime");
    }

    io.ReportProgressMessageMD('PRIME',  'selector', 'Branches to include in the PRIME analysis');


    for (_partition_, _selection_; in; prime.selected_branches) {
        _selection_ = utility.Filter (_selection_, '_value_', '_value_ == terms.tree_attributes.test');
        io.ReportProgressMessageMD('PRIME',  'selector', 'Selected ' + Abs(_selection_) + ' branches to include in the PRIME analysis: \`' + Join (', ',utility.Keys(_selection_)) + '\`');
    }

    estimators.fixSubsetOfEstimates(prime.gtr_results, prime.gtr_results[terms.global]);

    // Step 1 - Fit to MG94xREV
    namespace prime {
        doPartitionedMG ("prime", FALSE);
    }


    io.ReportProgressMessageMD ("PRIME", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");

    prime.run_full_mg94 = TRUE;
    
    if (Type (prime.save_intermediate_fits) == "AssociativeList") {
        if (None != prime.save_intermediate_fits[^"terms.data.value"]) {
            if (utility.Has (prime.save_intermediate_fits[^"terms.data.value"], "Full-MG94", "AssociativeList")) {
                prime.final_partitioned_mg_results = (prime.save_intermediate_fits[^"terms.data.value"])["Full-MG94"];
                if (utility.Has (prime.save_intermediate_fits[^"terms.data.value"], "Full-MG94-LF", "String")) {
                    ExecuteCommands ((prime.save_intermediate_fits[^"terms.data.value"])["Full-MG94-LF"]);
                    prime.run_full_mg94 = FALSE;
                }
            }        
        }
    }


    if (prime.run_full_mg94) {   
    
        prime.final_partitioned_mg_results = estimators.FitMGREV (prime.filter_names, prime.trees, prime.codon_data_info [terms.code], {
            terms.run_options.model_type: terms.local,
            terms.run_options.partitioned_omega: prime.selected_branches,
            terms.run_options.retain_lf_object: TRUE,
            terms.run_options.apply_user_constraints: prime.zero_branch_length_constrain,
            terms.run_options.optimization_settings: {
                "OPTIMIZATION_METHOD" : "coordinate-wise"
            }
        }, prime.partitioned_mg_results);
     
        /*prime.final_partitioned_mg_results = estimators.FitMGREV (prime.filter_names, prime.trees, prime.codon_data_info [terms.code], {
            terms.run_options.model_type: terms.local,
            terms.run_options.partitioned_omega: prime.selected_branches,
            terms.run_options.retain_lf_object: TRUE,
            terms.run_options.apply_user_constraints: prime.zero_branch_length_constrain,
            terms.run_options.optimization_settings: {
                "OPTIMIZATION_METHOD" : "coordinate-wise"
            }
        }, prime.partitioned_mg_results);*/

        if (Type (prime.save_intermediate_fits) == "AssociativeList") {
            (prime.save_intermediate_fits[^"terms.data.value"])["Full-MG94"] = prime.final_partitioned_mg_results;        
            Export (lfe, ^prime.final_partitioned_mg_results[^"terms.likelihood_function"]);
            (prime.save_intermediate_fits[^"terms.data.value"])["Full-MG94-LF"] = lfe;
            io.SpoolJSON (prime.save_intermediate_fits[^"terms.data.value"],prime.save_intermediate_fits[^"terms.data.file"]);      
        }
    }


    
    
    debug.spool = utility.GetEnvVariable ("DEBUG_CHECKPOINT_STORE");
    if (Type (debug.spool) == "String" ) {
        io.SpoolJSON (prime.final_partitioned_mg_results, debug.spool);
    }
} else {

    prime.final_partitioned_mg_results = io.ParseJSON (debug.checkpoint);
    io.ReportProgressMessageMD('PRIME',  'selector', 'Branches to include in the PRIME analysis');


    for (_partition_, _selection_; in; prime.selected_branches) {
        _selection_ = utility.Filter (_selection_, '_value_', '_value_ == terms.tree_attributes.test');
        io.ReportProgressMessageMD('PRIME',  'selector', 'Selected ' + Abs(_selection_) + ' branches to include in the PRIME analysis: \`' + Join (', ',utility.Keys(_selection_)) + '\`');
    }   
}

io.ReportProgressMessageMD("PRIME", "codon-refit", "* Log(L) = " + Format(prime.final_partitioned_mg_results[terms.fit.log_likelihood],8,2));
prime.global_dnds = selection.io.extract_global_MLE_re (prime.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio);

for (_value_; in; prime.global_dnds ) {
    io.ReportProgressMessageMD ("PRIME", "codon-refit", "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));
}


estimators.fixSubsetOfEstimates(prime.final_partitioned_mg_results, prime.final_partitioned_mg_results[terms.global]);

//Store MG94 to JSON
selection.io.json_store_lf_withEFV (prime.json,
                            terms.json.global_mg94xrev,
                            prime.final_partitioned_mg_results[terms.fit.log_likelihood],
                            prime.final_partitioned_mg_results[terms.parameters],
                            prime.sample_size,
                            utility.ArrayToDict (utility.Map (prime.global_dnds, "_value_", "{'key': _value_[terms.description], 'value' : Eval({{_value_ [terms.fit.MLE],1}})}")),
                            (prime.final_partitioned_mg_results[terms.efv_estimate])["VALUEINDEXORDER"][0],
                            prime.display_orders[terms.json.global_mg94xrev]);

utility.ForEachPair (prime.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(prime.json, terms.json.global_mg94xrev, terms.branch_length, prime.display_orders[terms.json.global_mg94xrev],
                                             _key_,
                                             selection.io.extract_branch_info((prime.final_partitioned_mg_results[terms.branch_length])[_key_], "selection.io.branch.length"));');


selection.io.stopTimer (prime.json [terms.json.timers], "Model fitting");

// define the site-level likelihood function

prime.site.background_fel = model.generic.DefineModel("models.codon.MG_REV.ModelDescription",
        "prime.background_fel", {
            "0": parameters.Quote(terms.local),
            "1": prime.codon_data_info[terms.code]
        },
        prime.filter_names,
        None);


prime.alpha = model.generic.GetLocalParameter (prime.site.background_fel, terms.parameters.synonymous_rate);
prime.beta  = model.generic.GetLocalParameter (prime.site.background_fel, terms.parameters.nonsynonymous_rate);

io.CheckAssertion ("None!=prime.alpha && None!=prime.beta", "Could not find expected local synonymous and non-synonymous rate parameters in \`estimators.FitMGREV\`");


prime.site.property_model =  model.generic.DefineModel("models.codon.MG_REV_PROPERTIES.ModelDescription",
        "prime.model.prop", {
            "0": parameters.Quote(terms.local),
            "1": prime.codon_data_info[terms.code],
            "2": parameters.Quote (prime.property_set) // property set to use
        },
        prime.filter_names,
        None);
   
      
prime.json [terms.model] = {};     
(prime.json [terms.model] ) [terms.model.residue_properties] = prime.site.property_model [terms.model.residue_properties];   
(prime.json [terms.model] ) [terms.translation_table] = prime.site.property_model [terms.translation_table];
(prime.json [terms.model] ) [terms.alphabet] = prime.site.property_model [terms.alphabet];         




prime.properties               = prime.site.property_model [terms.model.residue_properties];
prime.properties.N             = utility.Array1D (prime.properties);
prime.property_variable_prefix = "prime.site_lambda_";
prime.baseline_non_syn_rate    = "prime.site_beta";
prime.property_variables       = {};
prime.lambda_range = {
        terms.lower_bound: "-15",
        terms.upper_bound: "15"
    };
    
    
    
prime.table_headers = { 11 + 3*prime.properties.N, 2};
    
prime.table_headers[0][0] = "&alpha;"; prime.table_headers[0][1] = "Synonymous substitution rate at a site";
prime.table_headers[1][0] = "&beta;"; prime.table_headers[1][1] = "Property independent non-synonymous rate a site";
prime.table_headers[2][0] = "FEL &alpha;"; prime.table_headers[2][1] = "Synonymous substitution rate at a site (FEL model)";
prime.table_headers[3][0] = "FEL &beta;"; prime.table_headers[3][1] = "Log of property independent non-synonymous rate a site";
prime.table_headers[4][0] = "Total branch length"; prime.table_headers[4][1] = "The total length of branches contributing to inference at this site, and used to scale dN-dS";
prime.table_headers[5][0] = "# subs"; prime.table_headers[5][1] = "The number of amino-acid substitutions that occurred at this site";
prime.table_headers[6][0] = "# aa"; prime.table_headers[6][1] = "The number of unique amino-acids occurring at that site";
prime.table_headers[7][0] = "PRIME LogL"; prime.table_headers[7][1] = "Site Log-likelihood under the PRIME model";
prime.table_headers[8][0] = "FEL LogL"; prime.table_headers[8][1] = "Site Log-likelihood under the FEL model";
prime.table_headers[9][0] = "p-value"; prime.table_headers[9][1] = "Omnibus p-value (any property is important)";
prime.table_headers[10][0] = "q-value"; prime.table_headers[9][1] = "Omnibus q-value (any property is important)";


                         
prime.lambdas = {};                                                                      
prime.property_to_index = {};
io.ReportProgressMessageMD ("PRIME", "property-description", "Using the following set of **" + utility.Array1D (prime.properties) + "** properties");

prime.i = 1;
prime.local_to_property_name = {};
prime.p_value_indices = {"Overall" : 9};
prime.report.count = {"Overall" : 0};


for (key, value; in; prime.properties ) {
    io.ReportProgressMessageMD ("PRIME", "property-description", "* " + key);
    
    prime.p = prime.property_variable_prefix + (prime.site.property_model[terms.model.residue_name_map])[key];
    prime.table_headers [8+prime.i*3][0] = "&lambda;" + prime.i; prime.table_headers [7 + prime.i*3][1] = "Importance for " + key;
    prime.table_headers [9+prime.i*3][0] = "p" + prime.i; prime.table_headers [8 + prime.i*3][1] = "p-value for non-zero effect of " + key;
    prime.table_headers [10+prime.i*3][0] = "LogL" + prime.i; prime.table_headers [9 + prime.i*3][1] = "Log likelihood when there is no effect of " + key;
    
    model.generic.AddGlobal (prime.site.property_model,  prime.p , key);
    parameters.DeclareGlobal ( prime.p, {});
    parameters.SetRange ( prime.p, prime.lambda_range);
    
    parameter.local_lambda = model.generic.GetLocalParameter (prime.site.property_model , terms.propertyImportance(key));
    io.CheckAssertion ("None != parameter.local_lambda", "Could not find expected a local parameter in \`models.codon.MG_REV_PROPERTIES.ModelDescription\`");
    prime.lambdas [parameter.local_lambda] = prime.p;
    prime.local_to_property_name [parameter.local_lambda] = key;
    prime.property_to_index [key] = prime.i-1;
    prime.p_value_indices[key] = 9+prime.i*3;
    prime.report.count[key] = 0;
    prime.i += 1;

}


model.generic.AddGlobal (prime.site.property_model,  prime.baseline_non_syn_rate , prime.parameter_site_beta);
parameters.DeclareGlobal ( prime.baseline_non_syn_rate, {});
//parameters.SetRange ( prime.baseline_non_syn_rate, prime.lambda_range);
model.generic.AddGlobal (prime.site.property_model, "prime.site_alpha", prime.parameter_site_alpha);
parameters.DeclareGlobal ("prime.site_alpha", {});



prime.site_model_mapping = {
                           "prime.background_fel" : prime.site.background_fel,
                           "prime.model.prop" : prime.site.property_model,
                          };


selection.io.startTimer (prime.json [terms.json.timers], "PRIME analysis", 2);


model.generic.AddGlobal (prime.site.background_fel, "prime.site_alpha", prime.parameter_site_alpha);
model.generic.AddGlobal (prime.site.background_fel, "prime.site_beta_nuisance", prime.parameter_site_beta_nuisance);
parameters.DeclareGlobal ("prime.site_beta_nuisance", {});
model.generic.AddGlobal (prime.site.background_fel, "prime.site_beta", prime.parameter_site_beta);


function prime.rate_to_screen (name, value) {
    if (name == "Overall") {
        return "N/A";
    }
    return Format (value, 7, 3);
}

prime.report.significant_site = {{"" + (1+((prime.filter_specification[prime.report.partition])[terms.data.coverage])[prime.report.site]),
                                    prime.report.partition + 1,
                                    Format(prime.report.row[5],6,0),  // #subs
                                    Format(prime.report.row[6],6,0),  // #AA
                                    Format(prime.report.row[0],7,3),  // alpha
                                    Format(prime.report.row[1],7,3),  // beta
                                    prime.property_report_name,  // property name
                                    prime.rate_to_screen ( prime.property_report_name,,prime.report_rate), // property importance 
                                    Format(prime.report.row[prime.print_index],7,3), // property p-value 
                                    Join ("|", prime.site.composition.string[prime.property_report_name])
}};



prime.site_results        = {};
prime.imputed_leaf_states = {};
prime.sub_mapping         = {};
prime.q_values            = {};

for (prime.partition_index = 0; prime.partition_index < prime.partition_count; prime.partition_index += 1) {


    prime.report.header_done = FALSE;
    prime.table_output_options[utility.getGlobalValue("terms.table_options.header")] = TRUE;

    prime.model_to_branch_property = { "prime.model.prop" : utility.Filter (prime.selected_branches[prime.partition_index], '_value_', '_value_ == terms.tree_attributes.test'),
                                       "prime.background_fel" : utility.Filter (prime.selected_branches[prime.partition_index], '_value_', '_value_ != terms.tree_attributes.test')};


    model.ApplyModelToTree( "prime.site_tree_fel", prime.trees[prime.partition_index], {terms.default : prime.site.background_fel}, None);
    model.ApplyModelToTree( "prime.site_tree_property", prime.trees[prime.partition_index], None, prime.model_to_branch_property);

    prime.site_patterns = alignments.Extract_site_patterns ((prime.filter_specification[prime.partition_index])[utility.getGlobalValue("terms.data.name")]);
    

    for (_node_; in; prime.site_tree_fel) {
        _node_class_ = (prime.selected_branches[prime.partition_index])[_node_];
         if (_node_class_ != terms.tree_attributes.test) {
            _beta_scaler = "prime.site_beta_nuisance";
                prime.apply_proportional_site_constraint.fel ("prime.site_tree_property", _node_,
                    prime.alpha, prime.beta, "prime.site_alpha", _beta_scaler, (( prime.final_partitioned_mg_results[utility.getGlobalValue("terms.branch_length")])[prime.partition_index])[_node_]);  
                              
                prime.apply_proportional_site_constraint.fel ("prime.site_tree_fel", _node_,
                 prime.alpha, prime.beta, "prime.site_alpha", _beta_scaler, (( prime.final_partitioned_mg_results[utility.getGlobalValue("terms.branch_length")])[prime.partition_index])[_node_]);
                 
         } else {
                prime.apply_proportional_site_constraint.property ("prime.site_tree_property", _node_,
                    prime.alpha, "prime.site_alpha", prime.beta, prime.baseline_non_syn_rate,
                    prime.lambdas, (( prime.final_partitioned_mg_results[utility.getGlobalValue("terms.branch_length")])[prime.partition_index])[_node_]);
       
                prime.apply_proportional_site_constraint.fel ("prime.site_tree_fel", _node_,
                    prime.alpha, prime.beta, "prime.site_alpha", "prime.site_beta", (( prime.final_partitioned_mg_results[utility.getGlobalValue("terms.branch_length")])[prime.partition_index])[_node_]);
                    
         }
    }
    

    // create the likelihood function for this site

    ExecuteCommands (alignments.serialize_site_filter
                                       ((prime.filter_specification[prime.partition_index])[utility.getGlobalValue("terms.data.name")],
                                       ((prime.site_patterns[0])[utility.getGlobalValue("terms.data.sites")])[0],
                     ));

    __make_filter ("prime.site_filter");
    LikelihoodFunction prime.site_likelihood = (prime.site_filter, prime.site_tree_fel);

    
  

    estimators.ApplyExistingEstimates ("prime.site_likelihood", prime.site_model_mapping, prime.final_partitioned_mg_results,
                                        terms.globals_only);


    __make_filter ("prime.site_filter_property");
    LikelihoodFunction prime.site_likelihood_property = (prime.site_filter_property, prime.site_tree_property);

    
    /*
    Export (lfe, prime.site_likelihood_property);
    console.log (lfe);
    return 0;
    */
    

    estimators.ApplyExistingEstimates ("prime.site_likelihood_property", prime.site_model_mapping, prime.final_partitioned_mg_results,
                                        terms.globals_only);
                                            Export (lf, prime.site_likelihood_property);
                                            
 
 
    if (prime.partition_index == 0) {
        prime.variable_substitutions = {};
        //for (param,mle;in;prime.final_partitioned_mg_results[terms.global]) {
        for (tag,param; in; (prime.site.property_model[terms.parameters])[terms.global]) {
            prime.existing = utility.GetByKey (prime.final_partitioned_mg_results[terms.global], tag, "AssociativeList");
            if (None != prime.existing) {
                prime.variable_substitutions [param] = prime.existing[terms.fit.MLE];
            }
        }
       
        prime.Q = prime.site.property_model [terms.model.rate_matrix];
        prime.Q [0][0] = "";
        for (i,j,e; in; prime.site.property_model [terms.model.rate_matrix]) {
            if (Type (e) == "String" && Abs (e)) {
                prime.Q[+i][+j] = Simplify(e, prime.variable_substitutions);
            } else {
                prime.Q[+i][+j] = "";
            }
        }
        
        
        (prime.json [terms.model]) [terms.model.rate_matrix] = prime.Q;
        (prime.json [terms.model]) [terms.parameters] = prime.local_to_property_name;
        (prime.json [terms.model]) [terms.code.ordering] = prime.property_to_index;
         
    }
    
    


    prime.queue = mpi.CreateQueue ({terms.mpi.LikelihoodFunctions: {{"prime.site_likelihood","prime.site_likelihood_property"}},
                                   terms.mpi.Models : {{"prime.site.background_fel","prime.site.property_model"}},
                                   terms.mpi.Headers : utility.GetListOfLoadedModules ("libv3/"),
                                   terms.mpi.Variables : {{"prime.selected_branches","prime.pairwise_counts","prime.codon_data_info","prime.lambdas","prime.impute_states","terms.prime_imputed_states"}}
                                  });




    /* run the main loop over all unique site pattern combinations */
    
    prime.pattern_count = 1;
    for (_pattern_, _pattern_info_; in; prime.site_patterns) {
          io.ReportProgressBar("", "Working on site pattern " + (prime.pattern_count) + "/" + Abs (prime.site_patterns));
           if (_pattern_info_[utility.getGlobalValue("terms.data.is_constant")]) {
                prime.store_results (-1,None,{"0" : "prime.site_likelihood",
                                             "1" : "prime.site_likelihood_property",
                                             "2" : None,
                                             "3" : prime.partition_index,
                                             "4" : _pattern_info_,
                                             "5" : prime.site_model_mapping
                                     });
            } else {
                mpi.QueueJob (prime.queue, "prime.handle_a_site", {"0" : "prime.site_likelihood",
                                                                 "1" : "prime.site_likelihood_property",
                                                                 "2" : alignments.serialize_site_filter
                                                                   ((prime.filter_specification[prime.partition_index])[utility.getGlobalValue("terms.data.name")],
                                                                   (_pattern_info_[utility.getGlobalValue("terms.data.sites")])[0]),
                                                                 "3" : prime.partition_index,
                                                                 "4" : _pattern_info_,
                                                                 "5" : prime.site_model_mapping
                                                                    },
                                                                    "prime.store_results");
            }
            prime.pattern_count  += 1;
    }  

    mpi.QueueComplete (prime.queue);
    prime.partition_matrix = {Abs (prime.site_results[prime.partition_index]), Rows (prime.table_headers)};
    
    for (_key_, _value_; in; prime.site_results[prime.partition_index]) {
        for (prime.index = 0; prime.index < Rows (prime.table_headers); prime.index += 1) {
            prime.partition_matrix [0+_key_][prime.index] = _value_[prime.index];
            
        }
        prime.q_values ["" + prime.partition_index + "|" + _key_ ] = _value_[prime.p_value_indices["Overall"]];
    }
    
    prime.site_results[prime.partition_index] = prime.partition_matrix;

}

io.ClearProgressBar ();

prime.q_values =  math.HolmBonferroniCorrection (prime.q_values);

//console.log (prime.q_values);

prime.q_index  = prime.p_value_indices["Overall"] + 1;
print.q_count  = 0;
for (prime.partition_index = 0; prime.partition_index < prime.partition_count; prime.partition_index += 1) {
    prime.N = Rows (prime.site_results[prime.partition_index]);
    //console.log ("\n" + prime.N + "\n");
    for (prime.i = 0; prime.i < prime.N; prime.i += 1) {
        (prime.site_results[prime.partition_index])[prime.i][prime.q_index] =  prime.q_values["" + prime.partition_index + "|" + prime.i];
        if ( (prime.site_results[prime.partition_index])[prime.i][prime.q_index] <= prime.pvalue) {
            prime.q_count += 1;
        }
    }
}

prime.json [terms.json.MLE ] = {terms.json.headers   : prime.table_headers,
                                terms.json.content : prime.site_results };

if (prime.impute_states) {
    (prime.json [terms.json.MLE])[terms.prime_imputed_states] = prime.imputed_leaf_states;
}

(prime.json [terms.json.MLE])[terms.substitutions] = prime.sub_mapping;

io.ReportProgressMessageMD ("PRIME", "results", "** Found _" + prime.q_count + "_ sites where the Overall property importance false discovery rate was at q <= " + prime.pvalue + "**");

for (key, value; in;  prime.report.count) {
    io.ReportProgressMessageMD ("PRIME", "results", "** Found _" + value + "_ sites where " + key + " property was important at p <= " + prime.pvalue + "**");
}

selection.io.stopTimer (prime.json [terms.json.timers], "Total time");
selection.io.stopTimer (prime.json [terms.json.timers], "PRIME analysis");

GetString (_hpv,HYPHY_VERSION,0);
prime.json[terms.json.runtime] = _hpv;


io.SpoolJSON (prime.json, prime.codon_data_info[terms.json.json]);

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

function prime.apply_proportional_site_constraint.fel (tree_name, node_name, alpha_parameter, beta_parameter, alpha_factor, beta_factor, branch_length) {

    prime.branch_length = (branch_length[terms.parameters.synonymous_rate])[terms.fit.MLE];

    node_name = tree_name + "." + node_name;
  
    ExecuteCommands ("
        `node_name`.`alpha_parameter` := (`alpha_factor`) * prime.branch_length__;
        `node_name`.`beta_parameter`  := (`beta_factor`)  * prime.branch_length__;
    ");
}
//----------------------------------------------------------------------------------------

function prime.apply_proportional_site_constraint.property (tree_name, node_name, alpha_parameter, alpha_factor, beta_parameter, beta_factor, lambdas, branch_length) {

    prime.branch_length = (branch_length[terms.parameters.synonymous_rate])[terms.fit.MLE];

    node_name = tree_name + "." + node_name;
    
    ExecuteCommands (node_name + "." + alpha_parameter + ":=(" + alpha_factor + ")*" + prime.branch_length__);
    ExecuteCommands (node_name + "." + beta_parameter + ":=(" + beta_factor + ")*" + prime.branch_length__);
    //parameters.SetRange (node_name + "." + beta_parameter , prime.lambda_range);
    
    for (local, glob; in; lambdas) {
        ExecuteCommands (node_name + "." + local + ":=" + glob);
        parameters.SetRange (node_name + "." + local, prime.lambda_range);
    }
}


//----------------------------------------------------------------------------------------
lfunction prime.handle_a_site (lf_fel, lf_prop, filter_data, partition_index, pattern_info, model_mapping) {

    //console.log (pattern_info);
    GetString   (lfInfo, ^lf_fel,-1);   

    //utility.SetEnvVariable ("VERBOSITY_LEVEL", 1);

    //TODO Datafilters hardcode, Trees hardcode.
    ExecuteCommands (filter_data);
    __make_filter ((lfInfo["Datafilters"])[0]);

    GetString (lfInfo, ^lf_prop,-1);
         
    __make_filter ((lfInfo["Datafilters"])[0]);

    bsrel_tree_id = (lfInfo["Trees"])[0];

    utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);
    utility.SetEnvVariable ("ASSUME_REVERSIBLE_MODELS", TRUE);

    ^"prime.site_alpha" = 1;
    ^"prime.site_beta_nuisance"  = 1;
    ^"prime.site_beta"  = 1;
    ^"prime.site_beta" :> 0;

    Optimize (results, ^lf_fel
                , {"OPTIMIZATION_METHOD" : "nedler-mead", "OPTIMIZATION_PRECISION": 1e-4}
    );
    fel = estimators.ExtractMLEs (lf_fel, model_mapping);
    fel[utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];
    fel.ll = results[1][0];
    //console.log ("\nFEL = " + results[1][0]);
    //console.log ("alpha = " + ^"prime.site_alpha");
    //console.log ("beta = " + ^"prime.site_beta");
  
    //Export (lfe, ^lf_prop);
    // fprintf ("/tmp/PRIME-site." + (pattern_info["sites"])[0] + ".bf",CLEAR_FILE,lfe);
    
    if (^"prime.site_beta" > 0) {
  
         // fit the universal alternative
         //console.log (fel);
         ^"prime.site_beta" = Eval (^"prime.site_beta");   
         //parameters.SetConstraint ("prime.site_beta", Eval(^"prime.site_beta"),"");
   
     
     
         start_grid = {};
         propN = utility.Array1D (^"prime.lambdas");
     
         point = {};
         point ["prime.site_alpha"] = ^"prime.site_alpha";
         point ["prime.site_beta"] = ^"prime.site_beta";
         save_beta = ^"prime.site_beta";
         point ["prime.site_beta_nuisance"] = ^"prime.site_beta_nuisance";
     
         for (l; in; ^"prime.lambdas") {
            point [l] = 0;
         }
     
         start_grid + point;
         propNames = utility.Values (^"prime.lambdas");
         ranges = {};
         r_unit = {
            ^"terms.lower_bound" : -1,
            ^"terms.upper_bound" : 1,
         };
     
         for (k,l; in; ^"prime.lambdas") {
            ranges[l] = r_unit;
         }
          
          
         for (sp = 0; sp < 20; sp += 1) {
            for (l; in; ^"prime.lambdas") {
              point [l] = 0.0;
            }
            point [propNames [Random (0, propN) $ 1]] = Random (-1,1);
            start_grid + point;
         }
     

         for (sp; in; estimators.LHC (ranges, 40)) {
            point = {};
        
            point ["prime.site_alpha"] = Random (^"prime.site_alpha" * 0.5, Max(0.05,^"prime.site_alpha") * 2.0);
            point ["prime.site_beta"] = Random (^"prime.site_beta" * 0.5, Max(0.05,^"prime.site_beta") * 2.0);
        
            for (l,v; in; sp) {
              point [l] = v [^"terms.fit.MLE"];
            }
            start_grid + point;
        
            point ["prime.site_alpha"] = ^"prime.site_alpha";
            point ["prime.site_beta"] = ^"prime.site_beta";
        
            start_grid + point;
        
         }
     
     
     
        // Export (lfe, ^lf_prop);
        // fprintf ("/Users/sergei/Desktop/PRIME/site." + (pattern_info["sites"])[0] + ".bf",CLEAR_FILE,lfe);
     
        utility.ToggleEnvVariable ("TOLERATE_CONSTRAINT_VIOLATION", TRUE);       
        //utility.ToggleEnvVariable ("VERBOSITY_LEVEL", 10);       
            
         Optimize (results, ^lf_prop, {
                "OPTIMIZATION_METHOD" : "nedler-mead",
                //"OPTIMIZATION_METHOD" : "gradient-descent",
                "OPTIMIZATION_START_GRID" : start_grid,
                "MAXIMUM_OPTIMIZATION_ITERATIONS" : 1000,
                "OPTIMIZATION_PRECISION": 1e-3
            });
        
        // Export (lfe, ^lf_prop);
        // fprintf ("/tmp/PRIME-site." + (pattern_info["sites"])[0] + ".bf",CLEAR_FILE,lfe);
    
        //console.log ("\n" + ^"LF_INITIAL_GRID_MAXIMUM_VALUE" + "\nGrid best"+  ^"LF_INITIAL_GRID_MAXIMUM" + " / optimized " + results[1][0] + "\n");
        Optimize (results, ^lf_prop);
        //console.log ("\n" +  results[1][0] + "\n");
        //fprintf (stdout, ^lf_prop);  
        altL = results[1][0];
     
        // fit all of the nulls
    
        constrained_models = {};

        alternative = estimators.ExtractMLEsOptions (lf_prop, model_mapping, {});
    
        while (TRUE) {
            done = TRUE;
            for (k,l; in; ^"prime.lambdas") {
                ^l = 0;
                // FORCE UPDATE DEPENDENT VARIABLES
                LFCompute (^lf_prop,LF_START_COMPUTE);
                LFCompute (^lf_prop,results);
                LFCompute (^lf_prop,LF_DONE_COMPUTE);
                ^l := 0;  
                Optimize (results, ^lf_prop,{
                    "OPTIMIZATION_METHOD" : "nedler-mead",
                    "MAXIMUM_OPTIMIZATION_ITERATIONS" : 1000,
                    "OPTIMIZATION_START_GRID" : start_grid,
                    "OPTIMIZATION_PRECISION": 1e-4
                });
                //console.log (k + " => " + (- results[1][0] + altL));
                if (results[1][0] - altL > 1e-2) {
                    done = FALSE;
                    break;
                }
                extra_its = 5;
                while (results[1][0] < fel.ll - 1e-2 && extra_its) {
                    //console.log ("WORSE THAN FEL; REOPTIMIZING...");
                    for (sp; in; estimators.LHC (ranges, 40)) {
                        point = {};
        
                        point ["prime.site_alpha"] = Random (^"prime.site_alpha" * 0.5, ^"prime.site_alpha" * 2.0);
                        point ["prime.site_beta"] = Random (^"prime.site_beta" * 0.5, ^"prime.site_beta" * 2.0);
        
                        for (l,v; in; sp) {
                          point [l] = v [^"terms.fit.MLE"];
                        }
                        start_grid + point;
        
                        point ["prime.site_alpha"] = ^"prime.site_alpha";
                        point ["prime.site_beta"] = ^"prime.site_beta";
        
                        start_grid + point;
        
                     }
                    Optimize (results, ^lf_prop,{
                        "OPTIMIZATION_METHOD" : "nedler-mead",
                        "MAXIMUM_OPTIMIZATION_ITERATIONS" : 1000,
                        "OPTIMIZATION_START_GRID" : start_grid,
                        "OPTIMIZATION_PRECISION": 1e-3
                    });
                    //console.log (k + " => " + (results[1][0] - fel.ll));
                    Optimize (results, ^lf_prop);
                    //console.log (k + " => " + (results[1][0] - fel.ll));
                    extra_its += (-1);
                }
                constrained_models[k] = estimators.ExtractMLEsOptions (lf_prop, model_mapping, {^"terms.globals_only" : TRUE});
                if (extra_its == 0) {
                    (constrained_models[k])[utility.getGlobalValue("terms.fit.log_likelihood")]  = 0;
                } else {
                    (constrained_models[k])[utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];           
                }
                ^l = 0;        
                estimators.ApplyExistingEstimates (lf_prop, model_mapping, alternative, ^"terms.globals_only");
            }
        
            if (done) {
                break;
            } else {
                  ^l = 0;      
                  Optimize (results, ^lf_prop,{
                    "OPTIMIZATION_METHOD" : "nedler-mead",
                    "MAXIMUM_OPTIMIZATION_ITERATIONS" : 1000,
                    "OPTIMIZATION_PRECISION": 1e-4
                  });
                  altL = results[1][0]; 
                  //console.log ("REOPTIMIZED : " + altL);
                  alternative = estimators.ExtractMLEsOptions (lf_prop, model_mapping, {});
                  //Export (lfe, ^lf_prop);
                  //fprintf ("/tmp/PRIME-site." + (pattern_info["sites"])[0] + ".bf",CLEAR_FILE,lfe);
             
            }
        }
    
        //console.log ("\nPRIME = " + altL);
        //console.log ("alpha = " + ^"prime.site_alpha");
        //console.log ("beta = " + ^"prime.site_beta");
 
        character_map = None;    
        if (^"prime.impute_states") {
            DataSet anc = ReconstructAncestors ( ^lf_prop, {{0}}, MARGINAL, DOLEAVES);
            GetString   (names, anc, -1);
            GetDataInfo (codon_chars, ^((lfInfo["Datafilters"])[0]) , "CHARACTERS");
        
            character_map = {};
            for (seq_id, seq_name; in; names) {
                character_map [seq_name] = {};
                for (char, char_support; in; (anc.marginal_support_matrix)[seq_id][-1]) {
                    if (char_support > 1e-6) {
                        (character_map [seq_name])[codon_chars[char]] = char_support;
                    }
                }
            }
        
        }
        
        utility.ToggleEnvVariable ("TOLERATE_CONSTRAINT_VIOLATION", None);       

        ancestral_info = ancestral.build (lf_prop,0,FALSE);
        branch_substitution_information = (ancestral.ComputeSubstitutionBySite (ancestral_info,0,None))[^"terms.substitutions"];
        branch_mapping = ancestral.ComputeCompressedSubstitutionsBySite (ancestral_info,0);
        DeleteObject (ancestral_info);
    
        //console.log (branch_substitution_information);
        //console.log ("END " +(pattern_info["sites"]) + " pattern");
        //console.log ("----------------------------------------------");
    
        alternative [utility.getGlobalValue("terms.fit.log_likelihood")] = altL;
 

        return {
                "fel" : fel,
                utility.getGlobalValue("terms.alternative") : alternative,
                utility.getGlobalValue("terms.Null"): Null,
                utility.getGlobalValue("terms.model.residue_properties") : constrained_models,
                utility.getGlobalValue("terms.substitutions") : branch_mapping,
                utility.getGlobalValue("terms.branch_selection_attributes") : branch_substitution_information,
                utility.getGlobalValue("terms.prime_imputed_states") : character_map
        };
    } else {
        ancestral_info = ancestral.build (lf_fel,0,FALSE);
        branch_substitution_information = (ancestral.ComputeSubstitutionBySite (ancestral_info,0,None))[^"terms.substitutions"];
        branch_mapping = ancestral.ComputeCompressedSubstitutionsBySite (ancestral_info,0);
        DeleteObject (ancestral_info);
     

        return {
                "fel" : fel,
                utility.getGlobalValue("terms.alternative") : None,
                utility.getGlobalValue("terms.Null"): None,
                utility.getGlobalValue("terms.model.residue_properties") : None,
                utility.getGlobalValue("terms.substitutions") : branch_mapping,
                utility.getGlobalValue("terms.branch_selection_attributes") : branch_substitution_information,
                utility.getGlobalValue("terms.prime_imputed_states") : character_map
        };
   
    
    }
}

/* echo to screen calls */

//----------------------------------------------------------------------------------------
function prime.report.echo (prime.report.site, prime.report.partition, prime.report.row) {

    
    for (key, value; in; prime.p_value_indices) {
        prime.print_row = None;
        
        if (prime.report.row [value] <= prime.pvalue) {
            prime.print_row              = prime.report.significant_site;
            prime.print_index            =  value;
            prime.report.count[key]     += 1;
            prime.property_report_name   = key;
            if (key == "Overall") {
                 prime.report_rate = 0;
            } else {
                prime.report_rate = prime.report.row[ prime.print_index-1];
            }
        }

         if (None != prime.print_row) {
            if (!prime.report.header_done) {
                io.ReportProgressMessageMD("PRIME", "" + prime.report.partition, "For partition " + (prime.report.partition+1) + " these sites are significant at p <=" + prime.pvalue + "\n");
                io.ReportProgressMessageMD("PRIME", "" + prime.report.partition, "Properties with positive importance factors are **conserved**, and those with negative -- **changing**\n");
                fprintf (stdout,
                    io.FormatTableRow (prime.table_screen_output,prime.table_output_options));
                prime.report.header_done = TRUE;
                prime.table_output_options[utility.getGlobalValue("terms.table_options.header")] = FALSE;
            }

            io.ClearProgressBar ();
            fprintf (stdout,
                io.FormatTableRow (prime.print_row,prime.table_output_options));
        }
    }

}

//----------------------------------------------------------------------------------------

lfunction prime.store_results (node, result, arguments) {

    sub_profile = result[utility.getGlobalValue("terms.branch_selection_attributes")];
    sub_map = result[utility.getGlobalValue("terms.substitutions")];
    pattern_info    = arguments [4];


    if (None != sub_profile) {   
        total_sub_count = 0;
        sub_by_count = {"Overall" : {}};

        for (p,v; in; ^"prime.properties") {
            sub_by_count[p] = {};     
        }
        
        sub_counts = 0;
        aa_counts = {};
        
        partition_index = arguments [3];

        
        for (i, v; in; sub_profile) {
            for (j, b; in; v) {
                c = Abs (b);
                total_sub_count += c;
                ai = (^"prime.codon_table")[i];
                aj = (^"prime.codon_table")[j];
                
 
                
                if (Abs (ai) && Abs (aj)) {
                    if (ai != aj) {
                        for (idx, br; in; b) {
                            if (((^"prime.selected_branches")[partition_index])[br] == ^"terms.tree_attributes.test") {
                                sub_counts += 1;
                            }
                        }
                    }
                    aa_counts [ai] = 1;
                    aa_counts [aj] = 1;
                }
                
                if (Abs (ai) == 0) {ai = "?";}
                if (Abs (aj) == 0) {aj = "?";}

 
                for (p, counts; in; sub_by_count) {           
                    if ((counts/c)==0) {
                        counts [c] = {};
                    }
            
                    delta = "";
                    if (p != "Overall") {
                        delta = "("+Format (Abs (((^"prime.properties")[p])[ai]-((^"prime.properties")[p])[aj]),0,2) + ")";
                    }
        
                    if (ai != aj) {
                        counts[c] + (((^"prime.codon_table")[i] + ">" + (^"prime.codon_table")[j]) + delta);
                    } else {
                        counts[c] + ai;
                    }
                }
            }
        }
        
        sub_profile = {};
        for (p,count; in; sub_by_count) {

            sorted_subs = {Abs (count), 1};
            j = 0;
            for (i, v; in; count) {
                sorted_subs[j] = -(+i);
                j += 1; 
            }
            sub_profile1 = {};
            for (i; in; sorted_subs % 0) {
                sub_profile1 + ("[" + (-i) + "]" + Join (",",count [-i]));
            }

            sub_profile [p] = sub_profile1;
        }
        ^"prime.site.composition.string" = sub_profile;
    }



    nrows = 11 + 3*^"prime.properties.N";
    result_row          = { nrows, 1 };
    
    for (i; in; ^"prime.p_value_indices") {
        result_row  [i] = 1.;
    }
    result_row  [(^"prime.p_value_indices")["Overall"] + 1]  = 1.; 
    
  
    if (None != result) { // not a constant site

        result_row [0] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], utility.getGlobalValue("prime.parameter_site_alpha"));
        result_row [1] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], utility.getGlobalValue("prime.parameter_site_beta"));
        result_row [5] = sub_counts;
        result_row [6] = utility.Array1D (aa_counts);
        result_row [8] = (result["fel"])[utility.getGlobalValue("terms.fit.log_likelihood")];
       
     
        if (None != result[utility.getGlobalValue("terms.alternative")]) {
            omninbus_ll = (result[utility.getGlobalValue("terms.alternative")])[utility.getGlobalValue("terms.fit.log_likelihood")];
            lrt = {utility.getGlobalValue("terms.LRT") : 2*(omninbus_ll-(result["fel"])[utility.getGlobalValue("terms.fit.log_likelihood")])};
            lrt [utility.getGlobalValue("terms.p_value")] = 1-CChi2(lrt[utility.getGlobalValue("terms.LRT")],^"prime.properties.N");
            result_row [2] = estimators.GetGlobalMLE (result["fel"], utility.getGlobalValue("prime.parameter_site_alpha"));
            result_row [3] = estimators.GetGlobalMLE (result["fel"], utility.getGlobalValue("prime.parameter_site_beta"));
            result_row [7] = (result[utility.getGlobalValue("terms.alternative")])[utility.getGlobalValue("terms.fit.log_likelihood")];
            p_values       = {"Overall" : lrt [utility.getGlobalValue("terms.p_value")]};
            //result_row [5] = lrt [utility.getGlobalValue("terms.p_value")];

            sum = 0;
            alternative_lengths = ((result[utility.getGlobalValue("terms.alternative")])[utility.getGlobalValue("terms.branch_length")])[0];
 
            for (_node_; in; ^"prime.site_tree_fel") {
                _node_class_ = ((^"prime.selected_branches")[partition_index])[_node_];
                if (_node_class_ == utility.getGlobalValue("terms.tree_attributes.test")) {
                        sum += (alternative_lengths[_node_])[utility.getGlobalValue("terms.json.MLE")];
                }
            }
            result_row [4] = sum;
        
        
            for (key, prop; in; result[utility.getGlobalValue ("terms.model.residue_properties")]) {
                property_index = ((^"prime.property_to_index")[(^"prime.local_to_property_name")[key]])*3+11;
                ll             = prop[utility.getGlobalValue("terms.fit.log_likelihood")];
                pv             = 1-CChi2 (2*(omninbus_ll-ll),1);
                rate           = (^"prime.local_to_property_name")[key];
                rate           = ((((result[utility.getGlobalValue("terms.alternative")])[utility.getGlobalValue("terms.global")])[rate])[utility.getGlobalValue("terms.fit.MLE")]);
                //console.log (omninbus_ll);
                result_row [property_index][0] = rate;
                p_values [(^"prime.local_to_property_name")[key]] = pv;
                result_row [property_index+2][0] = ll;
            }
        
       
            p_values = math.HolmBonferroniCorrection(p_values);
            result_row[(^"prime.p_value_indices")["Overall"]] = p_values["Overall"];
            for (key, prop; in; result[utility.getGlobalValue ("terms.model.residue_properties")]) {
                property_index = ((^"prime.p_value_indices")[(^"prime.local_to_property_name")[key]]);
                result_row [property_index] = p_values[(^"prime.local_to_property_name")[key]];
            }
        } else {
            result_row [2] = result_row [0];
            result_row [3] = result_row [1];
            result_row [7] = result_row [8];
        }
        
    }   
    
    //console.log (result_row);

    utility.EnsureKey (^"prime.site_results", partition_index);
    utility.EnsureKey (^"prime.imputed_leaf_states", partition_index);
    utility.EnsureKey (^"prime.sub_mapping", partition_index);

    sites_mapping_to_pattern = pattern_info[utility.getGlobalValue("terms.data.sites")];
    sites_mapping_to_pattern.count = utility.Array1D (sites_mapping_to_pattern);
    
    for (i = 0; i < sites_mapping_to_pattern.count; i+=1) {
        site_index = sites_mapping_to_pattern[i];
        ((^"prime.site_results")[partition_index])[site_index] = result_row;
        ((^"prime.sub_mapping")[partition_index])[site_index] = sub_map;
        ((^"prime.imputed_leaf_states")[partition_index])[site_index] = result[^"terms.prime_imputed_states"];
        prime.report.echo (site_index, partition_index, result_row);
    }
}
