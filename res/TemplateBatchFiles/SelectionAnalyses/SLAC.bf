RequireVersion("2.3");


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

LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");

LoadFunctionLibrary("modules/io_functions.ibf");
LoadFunctionLibrary("modules/selection_lib.ibf");


/*------------------------------------------------------------------------------
    Display analysis information
*/

slac.analysis_description = {
    terms.io.info: "SLAC (Single Likelihood Ancestor Counting)
    uses a maximum likelihood ancestral state reconstruction
    and minimum path substitution counting to estimate site - level
    dS and dN,
    and applies a simple binomial - based test to test
    if dS differs drom dN.
    The estimates aggregate information over all branches,
    so the signal is derived from
    pervasive diversification or conservation. A subset of branches can be selected
    for testing as well.
    Multiple partitions within a NEXUS file are also supported
    for recombination - aware analysis.
    ",
    terms.io.version: "2.00",
    terms.io.reference: "Not So Different After All: A Comparison of Methods for Detecting Amino Acid Sites Under Selection (2005). _Mol Biol Evol_ 22 (5): 1208-1222",
    terms.io.authors: "Sergei L Kosakovsky Pond and Simon DW Frost",
    terms.io.contact: "spond@temple.edu",
    terms.io.requirements: "in-frame codon alignment and a phylogenetic tree"
};
io.DisplayAnalysisBanner(slac.analysis_description);


/*------------------------------------------------------------------------------
    Environment setup
*/

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

/*------------------------------------------------------------------------------
    Globals
*/

slac.samples = 100;
    /**
         how many ancestral reconstruction replicates to run ;
         set to 0 to skip resampling (which is time consuming)
    */

slac.pvalue = 0.1;
    /**
        default cutoff for printing to screen
    */

slac.json = {
    terms.json.analysis: slac.analysis_description,
    terms.json.fits: {},
    terms.json.timers: {}
};
    /**
        The dictionary of results to be written to JSON at the end of the run
    */

slac.display_orders =   {terms.original_name: -1,
                        terms.json.nucleotide_gtr: 0,
                        terms.json.global_mg94xrev: 1
                       };


selection.io.startTimer (slac.json [terms.json.timers], "Total time", 0);

/*------------------------------------------------------------------------------
    Key word arguments
*/

KeywordArgument ("code", "Which genetic code should be used", "Universal");
KeywordArgument ("alignment", "An in-frame codon alignment in one of the formats supported by HyPhy");
KeywordArgument ("tree", "A phylogenetic tree (optionally annotated with {})", null, "Please select a tree file for the data:"); 
KeywordArgument ("branches",  "Branches to test", "All");
KeywordArgument ("samples",  "number of samples to assess ancestral reconstruction uncertainty", "100");
KeywordArgument ("pvalue",  "The p-value threshold to use when testing for selection", "0.1");
// One additional KeywordArgument ("output") is called below after namespace absrel.

slac.scaler_prefix = "SLAC.scaler";

slac.by_site   = "by-site";
slac.by_branch = "by-branch";
slac.AVERAGED  = "AVERAGED";
slac.RESOLVED  = "RESOLVED";
slac.sample_median = "sample-median";
slac.sample_2.5 = "sample-2.5";
slac.sample_97.5 = "sample-97.5";



slac.table_headers = {{"ES", "Expected synonymous sites"}
                      {"EN", "Expected non-synonymous sites"}
                      {"S", "Inferred synonymous substitutions"}
                      {"N", "Inferred non-synonymous substitutions"}
                      {"P[S]", "Expected proportion of synonymous sites"}
                      {"dS", "Inferred synonymous susbsitution rate"}
                      {"dN", "Inferred non-synonymous susbsitution rate"}
                      {"dN-dS", "Scaled by the length of the tested branches"}
                      {"P [dN/dS > 1]", "Binomial probability that S is no greater than the observed value, with P<sub>s</sub> probability of success"}
                      {"P [dN/dS < 1]", "Binomial probability that S is no less than the observed value, with P<sub>s</sub> probability of success"}
                      {"Total branch length", "The total length of branches contributing to inference at this site, and used to scale dN-dS"}};

slac.table_screen_output = {{"Codon", "Partition", "S", "N", "dS", "dN", "Selection detected?"}};
slac.table_output_options =  {terms.table_options.header : TRUE, terms.table_options.minimum_column_width : 16, terms.table_options.align : "center"};


namespace slac {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ("slac");
}

slac.samples = io.PromptUser ("\n>Select the number of samples used to assess ancestral reconstruction uncertainty [select 0 to skip]",100,0,100000,TRUE);
slac.pvalue  = io.PromptUser ("\n>Select the p-value threshold to use when testing for selection",0.1,0,1,FALSE);

KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'SLAC.json')", slac.codon_data_info [terms.json.json]);
slac.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");

io.ReportProgressMessageMD('SLAC',  'selector', 'Branches to include in the SLAC analysis');

utility.ForEachPair (slac.selected_branches, "_partition_", "_selection_",
    "_selection_ = utility.Filter (_selection_, '_value_', '_value_ == terms.tree_attributes.test'); io.ReportProgressMessageMD('SLAC',  'selector', 'Selected ' + Abs(_selection_) + ' branches to include in SLAC calculations: \\\`' + Join (', ',utility.Keys(_selection_)) + '\\\`')");


selection.io.startTimer (slac.json [terms.json.timers], "Model fitting",1 );

namespace slac {
    doGTR ("slac");
}


estimators.fixSubsetOfEstimates(slac.gtr_results, slac.gtr_results[terms.global]);

namespace slac {
    doPartitionedMG ("slac", TRUE);
}

selection.io.stopTimer (slac.json [terms.json.timers], "Model fitting");

slac.global_dnds = selection.io.extract_global_MLE_re (slac.partitioned_mg_results, "^" + terms.parameters.omega_ratio);


//Store MG94 to JSON
selection.io.json_store_lf_withEFV (slac.json,
                            terms.json.global_mg94xrev,
                            slac.partitioned_mg_results[terms.fit.log_likelihood],
                            slac.partitioned_mg_results[terms.parameters],
                            slac.sample_size,
                            utility.ArrayToDict (utility.Map (slac.global_dnds, "_value_", "{'key': _value_[terms.description], 'value' : Eval({{_value_ [terms.fit.MLE],1}})}")),
                            (slac.partitioned_mg_results[terms.efv_estimate])["VALUEINDEXORDER"][0],
                            slac.display_orders[terms.json.global_mg94xrev]);

utility.ForEachPair (slac.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(slac.json, terms.json.global_mg94xrev, terms.branch_length, slac.display_orders[terms.json.global_mg94xrev],
                                             _key_,
                                             selection.io.extract_branch_info((slac.partitioned_mg_results[terms.branch_length])[_key_], "selection.io.branch.length"));');


selection.io.startTimer (slac.json [terms.json.timers], "Primary SLAC analysis", 2);



slac.nucleotide_frequencies    = (slac.gtr_results[terms.efv_estimate])["VALUEINDEXORDER"][0];

/*_EFV_MATRIX0_ = {{1,AC__*pooledFreqs[1],pooledFreqs[2],AT__*pooledFreqs[3]}
				{AC__*pooledFreqs[0],1,CG__*pooledFreqs[2],CT__*pooledFreqs[3]}
				{pooledFreqs[0],CG__*pooledFreqs[3],1,GT__*pooledFreqs[3]}
				{AT__*pooledFreqs[0],CT__*pooledFreqs[3],GT__*pooledFreqs[2],1}};
*/


slac.counting_bias_matrix = {4,4}["1"];

for (slac.i = 0; slac.i < 4; slac.i += 1) {
    for (slac.j = slac.i + 1; slac.j < 4; slac.j += 1) {
        slac.counting_bias_matrix[slac.i][slac.j] = ((slac.partitioned_mg_results[terms.global])[terms.nucleotideRate (models.DNA.alphabet[slac.i], models.DNA.alphabet[slac.j])])[terms.fit.MLE];
        slac.counting_bias_matrix[slac.j][slac.i] = slac.counting_bias_matrix[slac.i][slac.j] * slac.nucleotide_frequencies [slac.i];
        slac.counting_bias_matrix[slac.i][slac.j] = slac.counting_bias_matrix[slac.i][slac.j] * slac.nucleotide_frequencies [slac.j];
    }
}


slac.counting_bias_array = {};
slac.counting_bias_array + slac.counting_bias_matrix;
slac.counting_bias_array + slac.counting_bias_matrix;
slac.counting_bias_array + slac.counting_bias_matrix;

//io.SpoolLF (slac.partitioned_mg_results[terms.likelihood_function], slac.codon_data_info[terms.data.file], "slac");
io.ReportProgressMessageMD("slac", "anc", "Performing joint maximum likelihood ancestral state reconstruction");
slac.counts    = genetic_code.ComputePairwiseDifferencesAndExpectedSites (slac.codon_data_info[terms.code], {terms.genetic_code.count_stop_codons : FALSE, terms.genetic_code.weighting_matrix : slac.counting_bias_array});
slac.results   = {};



slac.report_to_screen = {};
/**
    store pairs of {partition, site, positive selection (if false => negative selection)} to print to screen (after sampling done if necessary)
*/

slac.printed_header = FALSE;

slac.positive_p := Min ((((slac.results[slac.i]) [slac.by_site])[slac.RESOLVED])[slac.site][8], (((slac.results[slac.i]) [slac.by_site])[slac.AVERAGED])[slac.site][8]);
slac.negative_p := Min ((((slac.results[slac.i]) [slac.by_site])[slac.RESOLVED])[slac.site][9], (((slac.results[slac.i]) [slac.by_site])[slac.AVERAGED])[slac.site][9]);

slac.report_positive_site = {{"" + (1+((slac.filter_specification[slac.i])[terms.data.coverage])[slac.site]),
                                    slac.i + 1,
                                    slac.row[2],
                                    slac.row[3],
                                    slac.row[5],
                                    slac.row[6],
                                    "Pos. p = " + slac.row[8]}};

slac.report_negative_site = {{"" + (1+((slac.filter_specification[slac.i])[terms.data.coverage])[slac.site]),
                                    slac.i + 1,
                                    slac.row[2],
                                    slac.row[3],
                                    slac.row[5],
                                    slac.row[6],
                                    "Neg. p = " + slac.row[9]}};



for (slac.i = 0; slac.i < Abs (slac.filter_specification); slac.i += 1) {
    slac.printed_header_sampler = FALSE;
    slac.table_output_options[terms.table_options.header] = TRUE;

    slac.ancestors         = ancestral.build (slac.partitioned_mg_results[terms.likelihood_function], slac.i, None);
    slac.results           [slac.i] = slac.compute_the_counts (slac.ancestors["MATRIX"], slac.ancestors["TREE_AVL"], slac.ancestors["AMBIGS"], slac.selected_branches[slac.i], slac.counts);

    slac.partition_sites   = utility.Array1D ((slac.filter_specification[slac.i])[terms.data.coverage]);

    for (slac.site = 0; slac.site < slac.partition_sites; slac.site += 1) {
        slac.row = utility.Map ((((slac.results[slac.i]) [slac.by_site])[slac.RESOLVED])[slac.site][-1],
                             "_entry_",
                             "Format (0+_entry_, 0, 3)");

        slac.print_row = None;
        if (slac.negative_p <= slac.pvalue) {
            utility.EnsureKey (slac.report_to_screen, slac.i);
            (slac.report_to_screen[slac.i])[slac.site] = 9;
            slac.print_row = slac.report_negative_site;
        } else {
            if (slac.positive_p <= slac.pvalue) {
                utility.EnsureKey (slac.report_to_screen, slac.i);
                (slac.report_to_screen[slac.i])[slac.site] = 8;
                slac.print_row = slac.report_positive_site;
            }
        }


        if (None != slac.print_row) {
            if (!slac.printed_header_sampler) {
                io.ReportProgressMessageMD("slac", "anc" + slac.i, "For partition " + (slac.i+1) + " these sites are significant at p <=" + slac.pvalue + "\n");
                fprintf (stdout,
                    io.FormatTableRow (slac.table_screen_output,slac.table_output_options));
                slac.printed_header_sampler = TRUE;
                slac.table_output_options[terms.table_options.header] = FALSE;
            }
            fprintf (stdout,
                io.FormatTableRow (slac.print_row,slac.table_output_options));
        }
    }

    slac.branch_attributes = selection.substitution_mapper (slac.ancestors["MATRIX"], slac.ancestors["TREE_AVL"], slac.ancestors["AMBIGS"], slac.counts, slac.ancestors ["MAPPING"], slac.codon_data_info[terms.code]);

    /*
    selection.io.json_store_branch_attribute(slac.json, terms.original_name, terms.json.node_label, 0,
                                             slac.i,
                                             slac.name_mapping);
    */

    selection.io.json_store_branch_attribute(slac.json, terms.codon, terms.json.node_label, 0,
                                             slac.i,
                                             slac.branch_attributes[terms.codon]);

    selection.io.json_store_branch_attribute(slac.json, terms.amino_acid, terms.json.node_label, 1,
                                             slac.i,
                                             slac.branch_attributes[terms.amino_acid]);

    selection.io.json_store_branch_attribute(slac.json, terms.synonymous_sub_count, terms.json.branch_label, 0,
                                             slac.i,
                                             slac.branch_attributes[terms.synonymous_sub_count]);

    selection.io.json_store_branch_attribute(slac.json, terms.nonsynonymous_sub_count, terms.json.branch_label, 1,
                                             slac.i,
                                             slac.branch_attributes[terms.nonsynonymous_sub_count]);


}


slac.json [terms.json.MLE ] = {terms.json.headers   : slac.table_headers,
                               terms.json.content : slac.results };


io.SpoolJSON (slac.json, slac.codon_data_info[terms.json.json]);
selection.io.stopTimer (slac.json [terms.json.timers], "Primary SLAC analysis");



lfunction slac.handle_a_sample (lf, partition, branches, counts) {
    //fprintf (stdout, lf, ":", partition, ":", branches, ":", counts, "\n");
    slac.sampled   = ancestral.build (lf, partition, {"sample": TRUE});
    return slac.compute_the_counts (slac.sampled["MATRIX"], slac.sampled["TREE_AVL"], slac.sampled["AMBIGS"], branches, counts);
}

lfunction slac.handle_a_sample_callback (node, result, arguments) {
    (^"slac.sample.results") +  result;
}

if (slac.samples > 0) {
    slac.table_screen_output_samplers = {{"Codon", "Partition", "    S [median, IQR]    ", "    N [median, IQR]    ", "    dS [median, IQR]    ", "    dN [median, IQR]    ", "  p-value [median, IQR]  "}};
    slac.table_output_options_samplers =  {"header" : TRUE, "min-column-width" : 16, "align" : "center"};

    selection.io.startTimer (slac.json [terms.json.timers], "Ancestor sampling analysis", 3);

    io.ReportProgressMessageMD ("slac", "sampling", "Ancestor sampling analysis");
    io.ReportAnalysisStageMD ("Generating `slac.samples` ancestral sequence samples to obtain confidence intervals");

    utility.EnsureKey (slac.json, slac.sample_median);
    utility.EnsureKey (slac.json, slac.sample_2.5);
    utility.EnsureKey (slac.json, slac.sample_97.5);

    for (slac.i = 0; slac.i < Abs (slac.filter_specification); slac.i += 1) {

        slac.printed_header_sampler = FALSE;
        slac.table_output_options_samplers[terms.table_options.header] = TRUE;
        slac.sample.results = {};

        slac.queue = mpi.CreateQueue ({"LikelihoodFunctions": {{slac.partitioned_mg_results[terms.likelihood_function]}},
                                        "Headers" : {{"libv3/all-terms.bf"}},
                                        "Variables" : {{"slac.by_site","slac.AVERAGED","slac.RESOLVED","slac.by_branch"}}
                                        });


        for (slac.s = 0; slac.s < slac.samples; slac.s+=1) {

            io.ReportProgressBar("", "\tSample " + (slac.s+1) + "/" + slac.samples + " for partition " + (1+slac.i));

            mpi.QueueJob (slac.queue, "slac.handle_a_sample", {"0" : slac.partitioned_mg_results[terms.likelihood_function],
                                                                "1" : slac.i,
                                                                "2" : slac.selected_branches[slac.i],
                                                                "3":  slac.counts},
                                                                "slac.handle_a_sample_callback");

        }

        io.ReportProgressBar("", "Done with ancestral sampling              \n");

        mpi.QueueComplete (slac.queue);

        slac.extractor = {slac.samples, 1};
        slac.sites   = utility.Array1D ((slac.filter_specification[slac.i])[terms.data.coverage]);
        slac.columns = Rows (slac.table_headers);


        slac.json [slac.sample_median] + {slac.RESOLVED  : {slac.sites, slac.columns}, slac.AVERAGED  : {slac.sites, slac.columns}};
        slac.json [slac.sample_2.5]    + {slac.RESOLVED  : {slac.sites, slac.columns}, slac.AVERAGED  : {slac.sites, slac.columns}};
        slac.json [slac.sample_97.5]   + {slac.RESOLVED  : {slac.sites, slac.columns}, slac.AVERAGED  : {slac.sites, slac.columns}};
        slac.keys = {{slac.RESOLVED , slac.AVERAGED }};


        for (slac.s = 0; slac.s < slac.sites; slac.s += 1) {
            for (slac.c = 0; slac.c <slac.columns; slac.c+=1) {
                for (slac.key = 0; slac.key < Columns(slac.keys); slac.key += 1) {
                    slac.key_value = slac.keys[slac.key];
                    slac.col = slac._extract_vector (slac.extractor, slac.sample.results, slac.key_value, slac.s,slac.c);

                    (((slac.json[slac.sample_median])[slac.i])[slac.key_value])[slac.s][slac.c] = stats.Quantile (slac.col, 0.5);
                    (((slac.json[slac.sample_2.5])[slac.i])[slac.key_value])[slac.s][slac.c] = stats.Quantile (slac.col, 0.025);
                    (((slac.json[slac.sample_97.5])[slac.i])[slac.key_value])[slac.s][slac.c] = stats.Quantile (slac.col, 0.975);
                }
            }
            if ((slac.report_to_screen[slac.i])[slac.s]) {

                if (!slac.printed_header_sampler) {
                    io.ReportProgressMessageMD("slac", "sampling", "Resampling results for partition " + (slac.i+1) + "\n");
                    fprintf (stdout,
                        io.FormatTableRow (slac.table_screen_output_samplers,slac.table_output_options_samplers));
                    slac.printed_header_sampler = TRUE;
                    slac.table_output_options_samplers[terms.table_options.header] = FALSE;
                }

                slac.row_median = utility.Map ((((slac.json[slac.sample_median])[slac.i])[slac.RESOLVED ])[slac.s][-1],
                             "_entry_",
                             "Format (0+_entry_, 0, 2)");
                slac.row_25 = utility.Map ((((slac.json[slac.sample_2.5])[slac.i])[slac.RESOLVED ])[slac.s][-1],
                             "_entry_",
                             "Format (0+_entry_, 0, 2)");
                slac.row_975 = utility.Map ((((slac.json[slac.sample_97.5])[slac.i])[slac.RESOLVED ])[slac.s][-1],
                             "_entry_",
                             "Format (0+_entry_, 0, 2)");

                slac.pi = (slac.report_to_screen[slac.i])[slac.s];

                slac.print_row = {{"" + (1+((slac.filter_specification[slac.i])[terms.data.coverage])[slac.s]),
                                        slac.i + 1,
                                        slac.row_median[2] + " [" + slac.row_25[2] + "-" + slac.row_975[2] + "]",
                                        slac.row_median[3] + " [" + slac.row_25[3] + "-" + slac.row_975[3] + "]",
                                        slac.row_median[5] + " [" + slac.row_25[5] + "-" + slac.row_975[5] + "]",
                                        slac.row_median[6] + " [" + slac.row_25[6] + "-" + slac.row_975[6] + "]",
                                        slac.row_median[slac.pi] + " [" + slac.row_25[slac.pi] + "-" + slac.row_975[slac.pi] + "]"}};

                fprintf (stdout,
                    io.FormatTableRow (slac.print_row,slac.table_output_options_samplers));
              }
        }
    }

    selection.io.stopTimer (slac.json [terms.json.timers], "Ancestor sampling analysis");
}

selection.io.stopTimer (slac.json [terms.json.timers], "Total time");
io.SpoolJSON (slac.json, slac.codon_data_info[terms.json.json]);

if (Abs (slac.report_to_screen) == 0) {
    io.ReportProgressMessageMD ("slac", "results", "** No sites found to be under positive or negative selection at p <= " + slac.pvalue + "**");
}

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
    return to ["(((from[_MATRIX_ELEMENT_ROW_])[\"by-site\"])[key])[row][column]"]%0;
}



/*------------------------------------------------------------------------------------*/

lfunction slac.compute_the_counts (matrix, tree, lookup, selected_branches, counts) {
//#profile START;

    site_count = Columns (matrix);
    selected_branches = utility.Filter (selected_branches, "_value_", "_value_ == 'test'");



    selected_branches_count      = Abs (selected_branches);
    selected_branches_in_avl     = {selected_branches_count,1};
    selected_branches_lengths    = {selected_branches_count,1};
    selected_branches_parents    = {selected_branches_count,1};
    selected_branch_total_length = 0;
    k = 0;
    tip_count = 0;

    for (i = 1; i < Abs (tree); i+=1) {
        if (selected_branches [(tree[i])["Name"]]) {
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
        
    io.CheckAssertion ("`&selected_branch_total_length`>0", "SLAC cannot be applied to a branch selection with total zero branch length (i.e. no variation)");

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

    pairwise_eps = counts [utility.getGlobalValue("terms.genetic_code.EPS")];
    state_count  = Rows (pairwise_eps);
    pairwise_epn = counts [utility.getGlobalValue("terms.genetic_code.EPN")];
    pairwise_ops = counts [utility.getGlobalValue("terms.genetic_code.OPS")];
    pairwise_opn = counts [utility.getGlobalValue("terms.genetic_code.OPN")];
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

                if (this_state >= 0) { // child fully resolved (this means that the parent is fully resolved as well)
                    psi = this_state*state_count+parent_state;

                    fr = Eval (fully_resolved);
                    for (k = 0; k < 4; k += 1) {
                        report_averaged[s*column_count + k] += fr[k];
                        report_resolved[s*column_count + k] += fr[k];
                    }

                    relative_branch_length = 1;
                    fr = Eval (fully_resolved);
                    relative_branch_length = selected_branches_lengths[i] / selected_branch_total_length;
                    
                    for (k = 0; k < 4; k += 1) {
                        report_resolved_by_branch[i*column_count + k] += fr[k];
                        report_averaged_by_branch[i*column_count + k] += fr[k];
                    }

                } else {
                    if (this_state == -1) { // the child is fully missing; no counts here
                        if (parent_state != -1) { // but the parent IS resolved
                            by_site_scaler [s] += (-selected_branches_lengths[i]);
                            psi = parent_state*state_count+parent_state;

                            /*fr = Eval (fully_resolved);
                            
                            for (k = 0; k < 2; k += 1) {
                                report_averaged[s*column_count + k] += fr[k];
                                report_resolved[s*column_count + k] += fr[k];
                            }*/

                            relative_branch_length = 1;
                            fr = Eval (fully_resolved);
                            relative_branch_length = selected_branches_lengths[i] / selected_branch_total_length;
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
            if (i < 2) {
                k = by_site_scaler[s];
                (*mx)[s*column_count + column_count - 1] = k;

                if (k > 0) {
                    sc = selected_branch_total_length/k;
                    (*mx)[s*column_count + 0] = (*mx)[s*column_count + 0] * sc;
                    (*mx)[s*column_count + 1] = (*mx)[s*column_count + 1] * sc;
                }
            } else {
                (*mx)[s*column_count + column_count - 1] = selected_branches_lengths[s];
            }

            total_subs = (*mx)[s*column_count + 2] + (*mx)[s*column_count + 3];

            (*mx) [s*column_count + 4] = (*mx) [s*column_count]/((*mx) [s*column_count + 1]+(*mx) [s*column_count ]);
            (*mx) [s*column_count + 5] = (*mx) [s*column_count + 2]/(*mx) [s*column_count + 0];
            (*mx) [s*column_count + 6] = (*mx) [s*column_count + 3]/(*mx) [s*column_count + 1];
            if (k > 0) {
                (*mx) [s*column_count + 7] = ((*mx) [s*column_count + 6] - (*mx) [s*column_count + 5])/k;
            }

            if (total_subs) {

                syn_count = (*mx) [s*column_count + 2];
                (*mx) [s*column_count + 8] = slac.extendedBinTail(total_subs,(*mx) [s*column_count + 4],syn_count);
                if (syn_count == 0) {
                    (*mx) [s*column_count + 9] = 1;
                } else {
                    (*mx) [s*column_count + 9] = 1 - slac.extendedBinTail(total_subs,(*mx) [s*column_count + 4],Max (0, syn_count-1));
                }
         } else {
                (*mx) [s*column_count + 8] = 1;
                (*mx) [s*column_count + 9] = 1;

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



    return {utility.getGlobalValue("slac.by_site") :{utility.getGlobalValue("slac.AVERAGED") : report_averaged,
                        utility.getGlobalValue("slac.RESOLVED") : report_resolved},
          utility.getGlobalValue("slac.by_branch") :{"NAMES": report_branch_names,
                        utility.getGlobalValue("slac.AVERAGED") : report_averaged_by_branch,
                        utility.getGlobalValue("slac.RESOLVED") : report_resolved_by_branch}
                        };

}

