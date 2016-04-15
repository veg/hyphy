RequireVersion("2.3");

/*------------------------------------------------------------------------------
    Load library files
*/

LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");
LoadFunctionLibrary("libv3/terms-json.bf");

LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");

LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");

LoadFunctionLibrary("modules/io_functions.ibf");

Export (function_def, slac.compute_the_counts);

fprintf (stdout, function_def, "\n");

return 0;


/*------------------------------------------------------------------------------
    Display analysis information
*/

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
    Multiple partitions within a NEXUS file are also supported
    for recombination - aware analysis.
    ",
    "version": "2.00",
    "reference": "Not So Different After All: A Comparison of Methods for Detecting Amino Acid Sites Under Selection (2005). _Mol Biol Evol_ 22 (5): 1208-1222",
    "authors": "Sergei L Kosakovsky Pond and Simon DW Frost",
    "contact": "spond@temple.edu",
    "requirements": "in-frame codon alignment and a phylogenetic tree"
});



/*------------------------------------------------------------------------------
    Environment setup
*/

utility.setEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

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
    terms.json.fits: {},
    terms.json.timers: {},
};
    /**
        The dictionary of results to be written to JSON at the end of the run
    */

selection.io.startTimer (slac.json [terms.json.timers], "Total time", 0);

slac.scaler_prefix = "SLAC.scaler";

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
slac.table_output_options =  {"header" : TRUE, "min-column-width" : 16, "align" : "center"};

slac.codon_data_info = alignments.promptForGeneticCodeAndAlignment("slac.codon_data", "slac.codon_filter");

    /** example output
    {
        "sequences": 13,
        "sites": 897,
        "name-mapping": {
            "AF082576": "AF082576",
...
            "AF234767": "AF234767",
            "AF_231119": "AF-231119"
        },
        "partitions": {
            {
                "SPAN_1", "0-587"
            }
            ...
            {
                "SPAN_5", "1693-2690"
            }
        },
        "code": {
            {
                14, 13, 14, 13, 7, 7, 7, 7, 19, 5, 19, 5, 2, 2, 3, 2, 12, 11, 12, 11, 6, 6, 6, 6, 19, 19, 19, 19, 1, 1, 1, 1, 16, 15, 16, 15, 8, 8, 8, 8, 20, 20, 20, 20, 4, 4, 4, 4, 10, 9, 10, 9, 5, 5, 5, 5, 10, 17, 18, 17, 1, 0, 1, 0
            }
        },
        "file": "/Volumes/sergei-raid/hyphy2.2/tests/hbltests/libv3/data/partitioned.nex",
        "stop": "TAA,TAG,TGA",
        "dataset": "slac.codon_data",
        "datafilter": "slac.codon_filter"
    }

    */

slac.sample_size = slac.codon_data_info["sites"] * slac.codon_data_info["sequences"];
slac.codon_data_info["json"] = slac.codon_data_info["file"] + ".slac.json";
slac.name_mapping = slac.codon_data_info[terms.json.name_mapping];
    /**
        will contain "mapped" -> "original" associations with sequence names; or null if no mapping was necessary
    */

if (None == slac.name_mapping) { /** create a 1-1 mapping if nothing was done */
    slac.name_mapping = {};
    utility.forEach (alignments.getSequenceNames ("slac.codon_data"), "_value_", "slac.name_mapping[_value_] = _value_");
}

slac.partitions_and_trees = trees.loadAnnotatedTreeTopology.match_partitions (slac.codon_data_info[terms.json.partitions], slac.name_mapping);

    /**  this will return a dictionary of partition strings and trees; one set per partition, as in
    {
        "0": {
            "name:" ... ,
            "filter-string": "0-587",
            "tree": {
                "string": ...
                "string_with_lengths": ...
                "branch_lengths": {
                    "AF_231119": 0.00519475,
                    ...
                },
                "annotated_string": ... ,
                "model_map": {
                    "AF_231119": "",
                    ...
                },
                "partitioned": {
                    "AF_231119": "leaf",
                 },
                "model_list": {
                    {
                }
            }
        },
        ...
    */



slac.partition_count = Abs (slac.partitions_and_trees);

utility.forEachPair (slac.partitions_and_trees, "_key_", "_value_", ' (slac.partitions_and_trees[_key_])["filter-string"] = selection.io.adjust_partition_string (_value_["filter-string"], 3*slac.codon_data_info["sites"])');
    /**
        ensure that all partitions fall on codon boundaries if they are contiguous
    */


io.reportProgressMessage ("", ">Loaded a multiple sequence alignment with **" + slac.codon_data_info["sequences"] + "** sequences, **" + slac.codon_data_info["sites"] + "** codons, and **" + Abs (slac.partitions_and_trees) + "** partitions from \`" + slac.codon_data_info["file"] + "\`");

slac.selected_branches = selection.io.defineBranchSets(slac.partitions_and_trees);
    /**  this will return a dictionary of selected branches; one set per partition, like in
    {
        "0": {
            "NODE3": "test",
            "NODE6": "background",
    ...
            "NODE15": "test"
        },

        ...
        "4": {
            "NODE4": "test",
     ...
             "NODE2": "background"
        }
    }
    */

slac.samples = io.prompt_user ("\n>Select the number of samples used to assess ancestral reconstruction uncertainty [select 0 to skip]",100,0,100000,TRUE);
slac.pvalue  = io.prompt_user ("\n>Select the p-value used to for perform the test at",0.1,0,1,FALSE);

selection.io.json_store_key_value_pair (slac.json, terms.json.trees, terms.json.tree.newick, utility.map (slac.partitions_and_trees, "_pt_", '(_pt_["tree"])["string"]&&1'));
selection.io.json_store_key_value_pair (slac.json, terms.json.trees, "tested", slac.selected_branches);

io.reportProgressMessageMD('SLAC',  'selector', 'Branches to include in the SLAC analysis');

utility.forEachPair (slac.selected_branches, "_partition_", "_selection_",
    "_selection_ = utility.filter (_selection_, '_value_', '_value_ == terms.json.attribute.test'); io.reportProgressMessageMD('SLAC',  'selector', 'Selected ' + Abs(_selection_) + ' branches to include in SLAC calculations: \\\`' + Join (', ',utility.keys(_selection_)) + '\\\`')");


selection.io.startTimer (slac.json [terms.json.timers], "Model fitting", 1);

slac.filter_specification = alignments.defineFiltersForPartitions (slac.partitions_and_trees, "slac.codon_data" , "slac.filter.", slac.codon_data_info);
/** defines codon filters for each partition, and returns the (codon) sites mapped to each filter
{
    {
        "0": {
            "name": "slac.filter.SPAN_1",
            "coverage": {
                {
                    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, ...
                }
            }
        },
        "1": {
            "name": "...",
            "coverage": "..."
        },
        ...
    }
}
*/

selection.io.json_store_key_value_pair (slac.json, None, terms.json.partitions, slac.filter_specification);

io.reportProgressMessageMD ("SLAC", "nuc-fit", "Nucleotide fit");

io.reportAnalysisStageMD ("Obtaining branch lengths and nucleotide rates under the  GTR model");

slac.trees = utility.map (slac.partitions_and_trees, "_partition_", '_partition_["tree"]');
slac.filter_names = utility.map (slac.filter_specification, "_partition_", '_partition_["name"]');

slac.gtr_results = estimators.fitGTR(slac.filter_names,
                                     slac.trees,
                                     parameters.helper.tree_lengths_to_initial_values (slac.trees, None));


io.reportProgressMessageMD ("SLAC", "nuc-fit", "* Log(L) = " + Format (slac.gtr_results["LogL"], 8, 2));
estimators.fixSubsetOfEstimates(slac.gtr_results, slac.gtr_results["global"]);

io.reportProgressMessageMD ("SLAC", "codon-fit", "Obtaining the global omega estimate based on relative GTR branch lengths and nucleotide substitution biases");


slac.scaler_variables = utility.populateDict (0, slac.partition_count, "slac.scaler_prefix + '_' + _k_", "_k_");
utility.forEach (slac.scaler_variables, "_value_", "parameters.declareGlobal(_value_, None);parameters.set_value(_value_, 3);");
/** the previous two lines declare per-partition branch length scalers
    slac.scaler_prefix_0
    slac.scaler_prefix_1
    etc
*/


slac.partitioned_mg_results = estimators.fitMGREV(slac.filter_names, slac.trees, slac.codon_data_info ["code"], {
    "model-type": terms.local,
    "proportional-branch-length-scaler": slac.scaler_variables,
    "partitioned-omega": slac.selected_branches,
    "retain-lf-object": TRUE
}, slac.gtr_results);

io.reportProgressMessageMD("SLAC", "codon-fit", "* Log(L) = " + Format(slac.partitioned_mg_results["LogL"],8,2));
slac.global_dnds = selection.io.extract_global_MLE_re (slac.partitioned_mg_results, "^" + terms.omega_ratio);
utility.forEach (slac.global_dnds, "_value_", 'io.reportProgressMessageMD ("SLAC", "codon-fit", "* " + _value_["description"] + " = " + Format (_value_["MLE"],8,4));');

/** extract and report dN/dS estimates */

selection.io.stopTimer (slac.json [terms.json.timers], "Model fitting");


selection.io.json_store_lf(
    slac.json,
    "Global MG94xREV",
    slac.partitioned_mg_results["LogL"],
    slac.partitioned_mg_results["parameters"],
    slac.sample_size,
    utility.array_to_dict (utility.map (slac.global_dnds, "_value_", "{'key': _value_['description'], 'value' : Eval({{_value_ ['MLE'],1}})}"))
);



utility.forEachPair (slac.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(slac.json, "Global MG94xREV model", terms.json.attribute.branch_length, 0,
                                             _key_,
                                             selection.io.extract_branch_info((slac.partitioned_mg_results[terms.json.attribute.branch_length])[_key_], "selection.io.branch.length"));');


selection.io.startTimer (slac.json [terms.json.timers], "Primary SLAC analysis", 2);


io.spoolLF (slac.partitioned_mg_results["LF"], slac.codon_data_info["file"], "SLAC");
io.reportProgressMessageMD("SLAC", "anc", "Performing joint maximum likelihood ancestral state reconstruction");
slac.counts    = genetic_code.ComputePairwiseDifferencesAndExpectedSites (slac.codon_data_info["code"], {"count-stop-codons" : FALSE});
slac.results   = {};



slac.report_to_screen = {};
/**
    store pairs of {partition, site, positive selection (if false => negative selection)} to print to screen (after sampling done if necessary)
*/

slac.printed_header = FALSE;

slac.positive_p := Min ((((slac.results[slac.i]) ["by-site"])["RESOLVED"])[slac.site][8], (((slac.results[slac.i]) ["by-site"])["AVERAGED"])[slac.site][8]);
slac.negative_p := Min ((((slac.results[slac.i]) ["by-site"])["RESOLVED"])[slac.site][9], (((slac.results[slac.i]) ["by-site"])["AVERAGED"])[slac.site][9]);

slac.report_positive_site = {{"" + (1+((slac.filter_specification[slac.i])["coverage"])[slac.site]),
                                    slac.i + 1,
                                    slac.row[2],
                                    slac.row[3],
                                    slac.row[5],
                                    slac.row[6],
                                    "Pos. p = " + slac.row[8]}};

 slac.report_negative_site = {{"" + (1+((slac.filter_specification[slac.i])["coverage"])[slac.site]),
                                    slac.i + 1,
                                    slac.row[2],
                                    slac.row[3],
                                    slac.row[5],
                                    slac.row[6],
                                    "Neg. p = " + slac.row[9]}};


for (slac.i = 0; slac.i < Abs (slac.filter_specification); slac.i += 1) {
    slac.printed_header_sampler = FALSE;
    slac.table_output_options["header"] = TRUE;

    slac.ancestors         = ancestral.build (slac.partitioned_mg_results["LF"], slac.i, None);
    slac.results           [slac.i] = slac.compute_the_counts (slac.ancestors["MATRIX"], slac.ancestors["TREE_AVL"], slac.ancestors["AMBIGS"], slac.selected_branches[slac.i], slac.counts);

    slac.partition_sites   = utility.array1D ((slac.filter_specification[slac.i])["coverage"]);




    for (slac.site = 0; slac.site < slac.partition_sites; slac.site += 1) {
        slac.row = utility.map ((((slac.results[slac.i]) ["by-site"])["RESOLVED"])[slac.site][-1],
                             "_entry_",
                             "Format (0+_entry_, 0, 3)");

        slac.print_row = None;
        if (slac.negative_p <= slac.pvalue) {
            utility.dict.ensure_key (slac.report_to_screen, slac.i);
            (slac.report_to_screen[slac.i])[slac.site] = 9;
            slac.print_row = slac.report_negative_site;
        } else {
            if (slac.positive_p <= slac.pvalue) {
                utility.dict.ensure_key (slac.report_to_screen, slac.i);
                (slac.report_to_screen[slac.i])[slac.site] = 8;
                slac.print_row = slac.report_positive_site;
            }
        }


        if (None != slac.print_row) {
            if (!slac.printed_header_sampler) {
                io.reportProgressMessageMD("SLAC", "anc" + slac.i, "For partition " + (slac.i+1) + " these sites are significant at p <=" + slac.pvalue + "\n");
                fprintf (stdout,
                    io.format_table_row (slac.table_screen_output,slac.table_output_options));
                slac.printed_header_sampler = TRUE;
                slac.table_output_options["header"] = FALSE;
            }

            fprintf (stdout,
                io.format_table_row (slac.print_row,slac.table_output_options));
        }
    }

    slac.branch_attributes = slac.substituton_mapper (slac.ancestors["MATRIX"], slac.ancestors["TREE_AVL"], slac.ancestors["AMBIGS"], slac.counts, slac.ancestors ["MAPPING"], slac.codon_data_info["code"]);

    selection.io.json_store_branch_attribute(slac.json, "original name", terms.json.attribute.node_label, 0,
                                             slac.i,
                                             slac.name_mapping);

    selection.io.json_store_branch_attribute(slac.json, "codon", terms.json.attribute.node_label, 0,
                                             slac.i,
                                             slac.branch_attributes["codon"]);

    selection.io.json_store_branch_attribute(slac.json, "amino-acid", terms.json.attribute.node_label, 1,
                                             slac.i,
                                             slac.branch_attributes["amino-acid"]);

    selection.io.json_store_branch_attribute(slac.json, "synonymous substitution count", terms.json.attribute.branch_label, 0,
                                             slac.i,
                                             slac.branch_attributes["synonymous substitution count"]);

    selection.io.json_store_branch_attribute(slac.json, "non-synonymous substitution count", terms.json.attribute.branch_label, 1,
                                             slac.i,
                                             slac.branch_attributes["non-synonymous substitution count"]);


}


slac.json [terms.json.MLE ] = {terms.json.headers   : slac.table_headers,
                               terms.json.content : slac.results };


io.spool_json (slac.json, slac.codon_data_info["json"]);
selection.io.stopTimer (slac.json [terms.json.timers], "Primary SLAC analysis");


if (slac.samples > 0) {
    slac.table_screen_output_samplers = {{"Codon", "Partition", "    S [median, IQR]    ", "    N [median, IQR]    ", "    dS [median, IQR]    ", "    dN [median, IQR]    ", "  p-value [median, IQR]  "}};
    slac.table_output_options_samplers =  {"header" : TRUE, "min-column-width" : 16, "align" : "center"};

    selection.io.startTimer (slac.json [terms.json.timers], "Ancestor sampling analysis", 3);

    io.reportProgressMessageMD ("SLAC", "sampling", "Ancestor sampling analysis");
    io.reportAnalysisStageMD ("Generating `slac.samples` ancestral sequence samples to obtain confidence intervals");

    utility.dict.ensure_key (slac.json, "sample-median");
    utility.dict.ensure_key (slac.json, "sample-2.5");
    utility.dict.ensure_key (slac.json, "sample-97.5");

    for (slac.i = 0; slac.i < Abs (slac.filter_specification); slac.i += 1) {

        slac.printed_header_sampler = FALSE;
        slac.table_output_options_samplers["header"] = TRUE;
        slac.sample.results = {};

        for (slac.s = 0; slac.s < slac.samples; slac.s+=1) {
            slac.sampled   = ancestral.build (slac.partitioned_mg_results["LF"], slac.i, {"sample": TRUE});
            slac.sample.results + slac.compute_the_counts (slac.sampled["MATRIX"], slac.sampled["TREE_AVL"], slac.sampled["AMBIGS"], slac.selected_branches[slac.i], slac.counts);
            //io.reportProgressBar("", "\tSample " + (slac.s+1) + "/" + slac.samples + " for partition " + (1+slac.i));
        }

        slac.extractor = {slac.samples, 1};
        slac.sites   = utility.array1D ((slac.filter_specification[slac.i])["coverage"]);
        slac.columns = Rows (slac.table_headers);


        slac.json ["sample-median"] + {"RESOLVED" : {slac.sites, slac.columns}, "AVERAGED" : {slac.sites, slac.columns}};
        slac.json ["sample-2.5"]    + {"RESOLVED" : {slac.sites, slac.columns}, "AVERAGED" : {slac.sites, slac.columns}};
        slac.json ["sample-97.5"]   + {"RESOLVED" : {slac.sites, slac.columns}, "AVERAGED" : {slac.sites, slac.columns}};
        slac.keys = {{"RESOLVED", "AVERAGED"}};


        for (slac.s = 0; slac.s < slac.sites; slac.s += 1) {
            for (slac.c = 0; slac.c <slac.columns; slac.c+=1) {
                for (slac.key = 0; slac.key < Columns(slac.keys); slac.key += 1) {
                    slac.key_value = slac.keys[slac.key];
                    slac.col = slac._extract_vector (slac.extractor, slac.sample.results, slac.key_value, slac.s,slac.c);

                    (((slac.json["sample-median"])[slac.i])[slac.key_value])[slac.s][slac.c] = stats.quantile (slac.col, 0.5);
                    (((slac.json["sample-2.5"])[slac.i])[slac.key_value])[slac.s][slac.c] = stats.quantile (slac.col, 0.025);
                    (((slac.json["sample-97.5"])[slac.i])[slac.key_value])[slac.s][slac.c] = stats.quantile (slac.col, 0.975);
                }
            }
            if ((slac.report_to_screen[slac.i])[slac.s]) {

                if (!slac.printed_header_sampler) {
                    io.reportProgressMessageMD("SLAC", "sampling", "Resampling results for partition " + (slac.i+1) + "\n");
                    fprintf (stdout,
                        io.format_table_row (slac.table_screen_output_samplers,slac.table_output_options_samplers));
                    slac.printed_header_sampler = TRUE;
                    slac.table_output_options_samplers["header"] = FALSE;
                }

                slac.row_median = utility.map ((((slac.json["sample-median"])[slac.i])["RESOLVED"])[slac.s][-1],
                             "_entry_",
                             "Format (0+_entry_, 0, 2)");
                slac.row_25 = utility.map ((((slac.json["sample-2.5"])[slac.i])["RESOLVED"])[slac.s][-1],
                             "_entry_",
                             "Format (0+_entry_, 0, 2)");
                slac.row_975 = utility.map ((((slac.json["sample-97.5"])[slac.i])["RESOLVED"])[slac.s][-1],
                             "_entry_",
                             "Format (0+_entry_, 0, 2)");

                slac.pi = (slac.report_to_screen[slac.i])[slac.s];

                slac.print_row = {{"" + (1+((slac.filter_specification[slac.i])["coverage"])[slac.s]),
                                        slac.i + 1,
                                        slac.row_median[2] + " [" + slac.row_25[2] + "-" + slac.row_975[2] + "]",
                                        slac.row_median[3] + " [" + slac.row_25[3] + "-" + slac.row_975[3] + "]",
                                        slac.row_median[5] + " [" + slac.row_25[5] + "-" + slac.row_975[5] + "]",
                                        slac.row_median[6] + " [" + slac.row_25[6] + "-" + slac.row_975[6] + "]",
                                        slac.row_median[slac.pi] + " [" + slac.row_25[slac.pi] + "-" + slac.row_975[slac.pi] + "]"}};

                fprintf (stdout,
                    io.format_table_row (slac.print_row,slac.table_output_options_samplers));
              }
        }
    }

    selection.io.stopTimer (slac.json [terms.json.timers], "Ancestor sampling analysis");
}

selection.io.stopTimer (slac.json [terms.json.timers], "Total time");
io.spool_json (slac.json, slac.codon_data_info["json"]);

if (Abs (slac.report_to_screen) == 0) {
    io.reportProgressMessageMD ("SLAC", "results", "** No sites found to be under positive or negative selection at p <= " + slac.pvalue + "**");
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

lfunction slac.substituton_mapper (matrix, tree, ambig_lookup, pairwise_counts, code_to_codon, genetic_code) {

    /*
        Return a dictionary with the following entries


            "codon" : {
                "each branch name" : {
                    "codon" :          { an array of strings; codon per site }
                 }
             },
            "amino-acid" : {
                "each branch name" : {
                    { an array of strings; amino-acid per site }
                 }
             },
            "synonymous substitution count" : {
                "each branch name" : {
                    { an array of floats: for each site, # of synonymous substitutions at this branch/site; by convention set to 0 at the root node }
                 }
             },
            "non-synonymous substitution count" : {
                "each branch name" : {
                    { an array of floats: for each site, # of non-synonymous substitutions at this branch/site; by convention set to 0 at the root node }
                 }
             }

    */

    site_count   = Columns (matrix);
    branch_count = Rows (matrix);
    aa_mapping   = genetic_code.DefineCodonToAAMapping (genetic_code);
    integer_mapping = genetic_code.DefineIntegerToAAMapping (genetic_code, TRUE);

    result      = {"codon" : {},
                   "amino-acid" : {},
                   "synonymous substitution count" : {},
                   "non-synonymous substitution count" : {}};

    code_lookup = {"-1" : "-"};

    for (b = 0; b < branch_count; b += 1) {

        bname  = (tree[b+1])["Name"];
        parent = (tree[b+1])["Parent"] - 1;

        branch_info = {"codon" : {1, site_count},
                       "amino-acid" : {1, site_count},
                       "synonymous substitution count" : {1, site_count},
                       "non-synonymous substitution count" : {1, site_count}};


        for (s = 0; s < site_count; s += 1) {
            code        = matrix[b][s];
            parent_code = matrix[parent][s];

            (branch_info["codon"])[s] = code_to_codon[code];

            if (Type(code_lookup [code]) != "String") {
                if (code >= 0) {
                    code_lookup [code] = aa_mapping [code_to_codon[code]];

                } else {
                    collect_residues = {};
                    utility.forEach ( ambig_lookup[-code-2], "_for_each_", "`&collect_residues`[`&integer_mapping`[_for_each_]] = 1");
                    code_lookup [code] = Join ("", utility.keys (collect_residues));
                }
            }

            if (code >= 0) {
                (branch_info["synonymous substitution count"]) [s]     = (pairwise_counts["OPS"])[parent_code][code];
                (branch_info["non-synonymous substitution count"]) [s] = (pairwise_counts["OPN"])[parent_code][code];
            } else {
                if (code != -1) {
                    resolution = (ambig_lookup[-code-2])$(ambig_lookup[-code-2])["_MATRIX_ELEMENT_ROW_"];
                    resolution = resolution[resolution];
                    (branch_info["synonymous substitution count"]) [s] = + utility.map (resolution, "_mapper_", "(`&pairwise_counts`[\"OPS\"])[`&parent_code`][_mapper_]");
                    (branch_info["non-synonymous substitution count"]) [s] = + utility.map (resolution, "_mapper_", "(`&pairwise_counts`[\"OPN\"])[`&parent_code`][_mapper_]");
                }

            }

            (branch_info["amino-acid"])[s] = code_lookup[code];
        }

        utility.forEach (utility.keys (branch_info), "slac.substituton_mapper.key",
                         "(`&result`[slac.substituton_mapper.key])[`&bname`] = `&branch_info`[slac.substituton_mapper.key]");

     }

    return result;

}

/*------------------------------------------------------------------------------------*/

lfunction slac.compute_the_counts (matrix, tree, lookup, selected_branches, counts) {
//#profile START;

    site_count = Columns (matrix);
    selected_branches = utility.filter (selected_branches, "_value_", "_value_ == 'test'");



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


                            fr = Eval (fully_resolved);
                            for (k = 0; k < 2; k += 1) {
                                report_averaged[s*column_count + k] += fr[k];
                                report_resolved[s*column_count + k] += fr[k];
                            }

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



    return {"by-site" :{"AVERAGED" : report_averaged,
                        "RESOLVED" : report_resolved},
          "by-branch" :{"NAMES": report_branch_names,
                        "AVERAGED" : report_averaged_by_branch,
                        "RESOLVED" : report_resolved_by_branch}
                        };

}

