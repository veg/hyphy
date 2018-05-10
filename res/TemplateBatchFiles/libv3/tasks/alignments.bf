LoadFunctionLibrary("../IOFunctions.bf");
LoadFunctionLibrary("../all-terms.bf");


/**
 * Get metadata from existing dataset
 * @name alignments.ReadCodonDataSet
 * @param dataset_name - id of dataset
 * @returns {Dictionary} r - metadata pertaining to the dataset
 */
lfunction alignments.ReadCodonDataSet(dataset_name) {
    return alignments.ReadCodonDataSetFromPath(dataset_name, None);
}

/**
 * Reads dataset from a file path
 * @name alignments.LoadGeneticCode
 * @param {String} code name - name of the genetic code to load, or None to prompt
 * @returns {Dictionary} - metadata pertaining to the genetic code
 */


lfunction alignments.LoadGeneticCode (code) {
     if (Type (code) == "String") {
        ExecuteAFile(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TemplateModels" + DIRECTORY_SEPARATOR + "chooseGeneticCode.def",
            {"0" : code });
     } else {
        ExecuteAFile(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TemplateModels" + DIRECTORY_SEPARATOR + "chooseGeneticCode.def");
     }
     return {
        utility.getGlobalValue("terms.code") : _Genetic_Code,
        utility.getGlobalValue("terms.code.stops") : GeneticCodeExclusions,
        utility.getGlobalValue("terms.code.ordering") : _hyphyAAOrdering,
        utility.getGlobalValue("terms.code.mapping") : defineCodonToAA ()
     };
}
/**
 * Reads dataset from a file path
 * @name alignments.ReadCodonDataSetFromPath
 * @param {String} dataset_name - name of variable you would like to set the dataset to
 * @param {String} path - path to alignment file
 * @returns {Dictionary} r - metadata pertaining to the dataset
 */

lfunction alignments.ReadCodonDataSetFromPath(dataset_name, path) {
     code_info = alignments.LoadGeneticCode (None);
     return alignments.ReadCodonDataSetFromPathGivenCode (dataset_name, path, code_info[utility.getGlobalValue("terms.code")], code_info[utility.getGlobalValue("terms.code.stops")]);
}

/**
 * Reads a codon dataset from a file path given a genetic code
 * @name alignments.ReadCodonDataSetFromPathGivenCode
 * @param {String} dataset_name - name of variable you would like to set the dataset to
 * @param {String} path - path to alignment file
 * @param {Matrix} code - genetic code
 * @param {String} stop_codons - the list of stopcodons
 * @returns {Dictionary} r - metadata pertaining to the dataset
 */
lfunction alignments.ReadCodonDataSetFromPathGivenCode (dataset_name, path, code, stop_codons) {

    if (Type(path) == "String") {
        DataSet ^ dataset_name = ReadDataFile(path);
    } else {
        SetDialogPrompt ("Select a coding sequence alignment file");
        DataSet ^ dataset_name = ReadDataFile(PROMPT_FOR_FILE);
        path = ^ "LAST_FILE_PATH";
    }

    r = alignments.ReadNucleotideDataSet_aux(dataset_name);

    r * {
        utility.getGlobalValue("terms.code"): code,
        utility.getGlobalValue("terms.stop_codons"): stop_codons, // this is indeed stop, not stops
        utility.getGlobalValue("terms.data.file"): path,
        utility.getGlobalValue("terms.data.sequences"): Eval( ^ "`dataset_name`.species")
    };

    return r;
}

/**
 * Helper function to readCodonDataSetFromPath; do not use
 * @name alignments.ReadNucleotideDataSet_aux
 * @private
 * @param dataset_name
 * @returns {Dictionary} r - metadata pertaining to the dataset
 */
lfunction alignments.ReadNucleotideDataSet_aux(dataset_name) {

    partitions = None;
    dfpm = utility.GetEnvVariable("DATA_FILE_PARTITION_MATRIX");

    partitions = {
        {
            utility.getGlobalValue("terms.default"),
            ""
        }
    };


    if (Type(dfpm) == "Matrix") {
        if (Rows(dfpm) > 0) {
            partitions = Transpose(dfpm);
        }
    }

    return {
        utility.getGlobalValue("terms.data.sequences"): Eval("`dataset_name`.species"),
        utility.getGlobalValue("terms.data.sites"): Eval("`dataset_name`.sites"),
        utility.getGlobalValue("terms.data.name_mapping"): Eval("`dataset_name`.mapping"),
        utility.getGlobalValue("terms.data.partitions"): partitions
    };
}

/**
 * Get sequence names from existing dataset
 * @name alignments.GetSequenceNames
 * @param {String} dataset_name - name of dataset to get sequence names from
 * @returns {Matrix} list of sequence names
 */
lfunction alignments.GetSequenceNames(dataset_name) {
    GetString(result, ^ dataset_name, -1);
    return result;
}

/**
 * Get sequence from existing dataset by name
 * @name alignments.GetSequenceNames
 * @param {String} dataset_name - name of dataset to get sequence names from
 * @param {String|None} sequence_name - the name of the sequence to extract or None to set up the initial mapping
 * @returns {String} the corresponding sequence string
 */

lfunction alignments.GetSequenceByName (dataset_name, sequence_name) {
    if (None == sequence_name) {
        cache = utility.MatrixToDict (alignments.GetSequenceNames (dataset_name));
        return None;
    }

    assert (cache/sequence_name, "Invalid sequence name `sequence_name` for data set `dataset_name` in call to alignments.GetSequenceByName");
    GetDataInfo (seq_string, ^dataset_name, cache[sequence_name]);
    return seq_string;
}

/**
 * Get i-th sequence name/value from an alignment
 * @name alignments.GetIthSequence
 * @param {String} dataset_name - name of dataset to get sequence names from
 * @param {String} index - the name of the sequence to extract or None to set up the initial mapping
 * @returns {Dict} {"id" : sequence name, "sequence" : sequence data}
 */

lfunction alignments.GetIthSequence (dataset_name, index) {

    GetString   (seq_id, ^dataset_name, index);
    GetDataInfo (seq_string, ^dataset_name, index);
    return {utility.getGlobalValue("terms.id") : seq_id, utility.getGlobalValue("terms.data.sequence") : seq_string};
}

/**
 * Read Nucleotide dataset from file_name
 * @name alignments.ReadNucleotideDataSet
 * @param dataset_name - the name of the dataset you wish to use
 * @param file_name - path to file
 * @returns {Dictionary} r - metadata pertaining to the dataset
 */
lfunction alignments.ReadNucleotideDataSet(dataset_name, file_name) {
    if (Type(file_name) == "String") {
        DataSet ^dataset_name = ReadDataFile(file_name);
    } else {
        SetDialogPrompt ("Select a sequence alignment file");
        DataSet ^dataset_name = ReadDataFile(PROMPT_FOR_FILE);
        file_name = LAST_FILE_PATH;
    }

    result = alignments.ReadNucleotideDataSet_aux(dataset_name);
    result[utility.getGlobalValue("terms.data.file")] = file_name;
    return result;
}

/**
 * Read a protein dataset from file_name
 * @name alignments.ReadNucleotideDataSet
 * @param dataset_name - the name of the dataset you wish to use
 * @param file_name - path to file
 * @returns {Dictionary} r - metadata pertaining to the dataset
 */
lfunction alignments.ReadProteinDataSet(dataset_name, file_name) {
    result = alignments.ReadNucleotideDataSet (dataset_name, file_name);
    
    /* check that the alignment has protein data */
    DataSetFilter throwaway = CreateFilter (^dataset_name, 1);
    GetDataInfo (alphabet, throwaway, "CHARACTERS");
    DeleteObject (throwaway);
    io.CheckAssertion("alignments.AlphabetType(`&alphabet`)==utility.getGlobalValue ('terms.amino_acid')", 
                      "The input alignment must contain protein data");
                      
    return result;
}

/**
 * Categorize an alphabet for an alignment 
 * @name alignments.AlphabetType
 * @param alphabet - the alphabet vector, e.g. fetched by GetDataInfo (alphabet, ..., "CHARACTERS");
 * @returns {String} one of standard alphabet types or None if unknown
 */
lfunction alignments.AlphabetType (alphabet) {
    alphabet = Join ("", alphabet);
    if (alphabet == "ACGT") {
        return utility.getGlobalValue ("terms.nucleotide");
    } else {
        if (alphabet == "ACDEFGHIKLMNPQRSTVWY") {
            return utility.getGlobalValue ("terms.amino_acid");
        }
    }
    return None;
}   

/**
 * Ensure that name mapping is not None by creating a f(x)=x map if needed
 * @name alignments.EnsureMapping
 * @param dataset_name - the name of the dataset you wish to use
 * @param {Dictionary} r - metadata pertaining to the dataset
 * @param file_name - path to file
 * @returns {Dictionary} r - metadata pertaining to the dataset
 */
 
lfunction alignments.EnsureMapping(dataset_name, data) {
    name_mapping = data[utility.getGlobalValue("terms.data.name_mapping")];
    if (None == name_mapping) { /** create a 1-1 mapping if nothing was done */
        name_mapping = {};
        utility.ForEach (alignments.GetSequenceNames (dataset_name), "_value_", "`&name_mapping`[_value_] = _value_");
        data[utility.getGlobalValue("terms.data.name_mapping")] = name_mapping;
    }    
    return data;
}

/**
 * Read dataset from data
 * @name alignments.ReadNucleotideDataSetString
 * @param {String} dataset_name - the name of the dataset you wish to use
 * @param {String} data - the data you wish to parse
 * @returns {Dictionary} r - metadata pertaining to the dataset
 */
function alignments.ReadNucleotideDataSetString(dataset_name, data) {
    ExecuteCommands("DataSet `dataset_name` = ReadFromString (data);");
    return alignments.ReadNucleotideDataSet_aux(dataset_name);
}

/**
 * Prompts user for genetic code and alignment
 * @name alignments.PromptForGeneticCodeAndAlignment
 * @param {String} dataset_name - the name  of the dataset you wish to use
 * @param {String} datafilter_name - the name  of the dataset filter you wish to use
 * @returns {Dictionary} r - metadata pertaining to the dataset
 */
function alignments.PromptForGeneticCodeAndAlignment(dataset_name, datafilter_name) {
    return alignments.LoadCodonDataFile(dataset_name, datafilter_name, alignments.ReadCodonDataSet(dataset_name));
}

/**
 * Loads genetic code and alignment from file path
 * @name alignments.LoadGeneticCodeAndAlignment
 * @param {String} dataset_name - the name  of the dataset you wish to use
 * @param {String} datafilter_name - the name  of the dataset filter you wish to use
 * @param {String} path - path to file
 * @returns {Dictionary} r - metadata pertaining to the dataset
 */
function alignments.LoadGeneticCodeAndAlignment(dataset_name, datafilter_name, path) {
    return alignments.LoadCodonDataFile(dataset_name, datafilter_name, alignments.ReadCodonDataSetFromPath(dataset_name, path));
}

/**
 * Creates codon dataset filter from data set
 * @name alignments.LoadCodonDataFile
 * @param {String} datafilter_name - the name  of the dataset filter you wish to use
 * @param {String} dataset_name - the name of the existing dataset
 * @param {Dictionary} data_info - DataSet metadata information
 * @returns {Dictionary} updated data_info that includes the number of sites, dataset, and datafilter name
 */
function alignments.LoadCodonDataFile(dataset_name, datafilter_name, data_info) {
    DataSetFilter ^ datafilter_name = CreateFilter( ^ dataset_name, 3, , , data_info[terms.stop_codons]);
    io.CheckAssertion("`datafilter_name`.sites*3==`dataset_name`.sites", "The input alignment must not contain stop codons");
    data_info[terms.data.sites] = ^ "`datafilter_name`.sites";
    data_info[terms.data.dataset] = dataset_name;
    data_info[terms.data.datafilter] = datafilter_name;
    return data_info;
}

/**
 * Creates nucleotide dataset filter from file
 * @name alignments.ReadNucleotideAlignment
 * @param {String} file_name - The location of the alignment file
 * @param {String} datafilter_name - the name  of the dataset filter you wish to use
 * @param {String} dataset_name - the name of the dataset you wish to use
 * @returns {Dictionary} updated data_info that includes the number of sites, dataset, and datafilter name
 */
function alignments.ReadNucleotideAlignment(file_name, dataset_name, datafilter_name) {
    data_info = alignments.ReadNucleotideDataSet(dataset_name, file_name);
    ExecuteCommands("DataSetFilter `datafilter_name` = CreateFilter (`dataset_name`,1)");
    data_info[terms.data.sites] = Eval("`datafilter_name`.sites");
    data_info[terms.data.dataset] = dataset_name;
    data_info[terms.data.datafilter] = datafilter_name;

    return data_info;
}

/**
 * Take an input filter, replace all identical sequences with a single copy
 * Optionally, rename the sequences to indicate copy # by adding ':copies'
 * @name alignments.CompressDuplicateSequences
 * @param {String} filter_in - The name of an existing filter
 * @param {String} filter_out - the name  (to be created) of the filter where the compressed sequences will go)
 * @param {Bool} rename - if true, rename the sequences
 * @returns the number of unique sequences
 */
lfunction alignments.CompressDuplicateSequences (filter_in, filter_out, rename) {
    GetDataInfo (duplicate_info, ^filter_in, -2);
    
    DataSetFilter ^filter_out = CreateFilter (^filter_in, 1, "", Join (",", duplicate_info["UNIQUE_INDICES"]));
    
    if (rename) {
        utility.ForEachPair (duplicate_info["UNIQUE_INDICES"],
                             "_idx_",
                             "_seq_idx_",
                             '
                                GetString (_seq_name_, ^`&filter_in`, _seq_idx_);
                                _seq_name_ += ":" + ((`&duplicate_info`)["UNIQUE_COUNTS"])[_idx_[1]];
                                SetParameter (^`&filter_in`,_seq_idx_,_seq_name_);
                              ');
       
    } 
    
    
    
    return duplicate_info["UNIQUE_SEQUENCES"];
    //DataSetFilter ^filter_out = CreateFilter (filter_in);
}

/**
 * Defines filters for multiple partitions
 * @name alignments.DefineFiltersForPartitions
 * @param {Matrix} partitions - a row vector of dictionaries with partition information containing "name" and "filter-string" attributes
 * @param {DataSet} source_data - the existing dataset to partition
 * @returns {Matrix} filters - a list of newly created dataset filters
 */
lfunction alignments.DefineFiltersForPartitions(partitions, source_data, prefix, data_info) {
    part_count = utility.Array1D(partitions);
    filters = {};
    if (utility.CheckKey(data_info, utility.getGlobalValue("terms.code"), "Matrix")) {
        for (i = 0; i < part_count; i += 1) {
            this_filter = {};
            DataSetFilter test = CreateFilter( ^ source_data, 1, (partitions[i])[utility.getGlobalValue("terms.data.filter_string")]);
            this_filter[utility.getGlobalValue("terms.data.name")] = prefix + (partitions[i])[utility.getGlobalValue("terms.data.name")];
            DataSetFilter ^ (this_filter[utility.getGlobalValue("terms.data.name")]) = CreateFilter( ^ source_data, 3, (partitions[i])[utility.getGlobalValue("terms.data.filter_string")], , data_info[utility.getGlobalValue("terms.stop_codons")]);
            diff = test.sites - 3 * ^ (this_filter[utility.getGlobalValue("terms.data.name")] + ".sites");
            
            //TODO: BELOW, IS THE "names" CORRECT OR SHOULD IT BE "name"????? SJS can't locate another time when the plural is used through libv3.
            io.CheckAssertion("`&diff` == 0", "Partition " + (filters["names"])[i] + " is either has stop codons or is not in frame");
            
            this_filter[utility.getGlobalValue("terms.data.coverage")] = utility.DictToArray(utility.Map(utility.Filter( ^ (this_filter[utility.getGlobalValue("terms.data.name")] + ".site_map"), "_value_", "_value_%3==0"), "_value_", "_value_$3"));
            filters + this_filter;
        }

    } else {
        for (i = 0; i < part_count; i += 1) {
            this_filter = {};
            this_filter[utility.getGlobalValue("terms.data.name")] = prefix + (partitions[i])[utility.getGlobalValue("terms.data.name")];
            DataSetFilter ^ (this_filter[utility.getGlobalValue("terms.data.name")]) = CreateFilter( ^ source_data, 1, (partitions[i])[utility.getGlobalValue("terms.data.filter_string")]);
            this_filter[utility.getGlobalValue("terms.data.coverage")] = ^ (this_filter[utility.getGlobalValue("terms.data.name")] + ".site_map");
            filters + this_filter;

        }

    }
    return filters;
}

/**
 * @name alignments.serialize_site_filter
 * @param {DataFilter} data_filter
 * @param {Number} site_index
 * @returns {String} a string representation of data_filter
 */
lfunction alignments.serialize_site_filter (data_filter, site_index) {
    GetDataInfo (fi, ^data_filter, "PARAMETERS");
    utility.ToggleEnvVariable ("DATA_FILE_PRINT_FORMAT", 6);
    utility.ToggleEnvVariable ("IS_TREE_PRESENT_IN_DATA", FALSE);
    DataSetFilter temp = CreateFilter (^data_filter,fi["ATOM_SIZE"],'' + site_index*fi["ATOM_SIZE"] + '-' + ((site_index+1)*fi["ATOM_SIZE"]-1),'',fi["EXCLUSIONS"]);
    Export (filter_string, temp);
    utility.ToggleEnvVariable ("DATA_FILE_PRINT_FORMAT", None);
    utility.ToggleEnvVariable ("IS_TREE_PRESENT_IN_DATA", None);
    // TODO filter string below
    return '
        lfunction __make_filter (name) {
            DataSet hidden = ReadFromString ("`filter_string`");
            DataSetFilter ^name = CreateFilter (hidden, `""+fi['ATOM_SIZE']`,,,"`fi['EXCLUSIONS']`");
        };
    ';
}

/**
 * @name alignments.TranslateCodonsToAminoAcids
 * Translate a codon sequence to amino-acids using the mapping provided by the
 * genetic code
 * @param {String} sequence - the string to translate
 * @param {Number} offset - start at this position (should be in {0,1,2})
 * @param {Dictionary} code - genetic code description (e.g. returned by alignments.LoadGeneticCode)
 * @returns {String} the amino-acid translation ('?' is used to represent ambiguities; 'X' - stop codons)
 */

lfunction alignments.TranslateCodonsToAminoAcids (sequence, offset, code) {
    l = Abs (sequence);
	translation = "";
	translation * (l/3+1);
	for (k = offset; k < l; k += 3) {
		codon = sequence[k][k+2];
		if (code [utility.getGlobalValue("terms.code.mapping")] / codon) {
		    translation * (code [utility.getGlobalValue("terms.code.mapping")])[codon];
		}
		else {
		    if (codon == "---") {
			    translation * "-";
		    } else {
			    translation * "?";
			}
	    }
	}
	translation * 0;
	return translation;
}

/**
 * @name alignments.TranslateCodonsToAminoAcidsWithAmbigs
 * Translate a codon sequence to amino-acids using the mapping provided by the
 * genetic code
 * @param {String} sequence - the string to translate
 * @param {Number} offset - start at this position (should be in {0,1,2})
 * @param {Dictionary} code - genetic code description (e.g. returned by alignments.LoadGeneticCode)
 * @param {lookup} code - resolution lookup dictionary 
 * @returns {Dict} list of possible amino-acids (as dicts) at this position
 */

lfunction alignments.TranslateCodonsToAminoAcidsWithAmbiguities (sequence, offset, code, lookup) {
    console.log (sequence);
    
    l = Abs (sequence);
	translation = {};
    
    DataSet single_seq                  = ReadFromString (">s\n" + sequence[offset][Abs (sequence)-1]);
    DataSetFilter single_seq_filter     = CreateFilter   (single_seq, 3, "", "");
    
    GetDataInfo (patterns, single_seq_filter);
    GetDataInfo (alphabet, single_seq_filter, "CHARACTERS");
    GetDataInfo (single_seq_data, single_seq_filter, 0);
    code_lookup = code [utility.getGlobalValue("terms.code.ordering")];
    code_table  = code [utility.getGlobalValue("terms.code")];
    
    for (s = 0; s < single_seq_filter.sites; s += 1) {
        codon = single_seq_data[3*s][3*s+2];
        if (lookup / codon) {
            translation[s] = lookup [codon];
        } else {
            GetDataInfo (resolutions, single_seq_filter, 0, patterns[s]);
            resolution_count = +resolutions;
            my_resolution = {};
            if (resolution_count == 1) {
                my_resolution [code_lookup[code_table[(Max (resolutions, 1))[1]]]] = 1;
            } else {
                if (resolution_count == 0 || resolution_count == 64) {
                    my_resolution["-"] = 1;
                } else {
                    for (r = 0; r < 64; r += 1) {
                        if (resolutions[r]) {
                            my_resolution [code_lookup[code_table[r]]] = 1;
                        }
                    }
                }
            }
            lookup [codon] = my_resolution;
            translation[s] = my_resolution;
        }
    }
    
    
	return translation;
}

/**
 * @name alignments.MapAlignmentToReferenceCoordinates
 * Map a query sequence from the aligned coordinates
 * genetic code
 * @param {String} sequence - the string to translate
 * @returns {String} the amino-acid translation ('?' is used to represent ambiguities; 'X' - stop codons)
 */

lfunction alignments.MapAlignmentToReferenceCoordinates (reference, aligned_reference, aligned_qry, offset) {


    realigned         = {};
    mapping           = {};
    coordinates       = {1,Abs(reference)};
    reduced_alignment = {1,Abs(reference)};
    /* this will contain a list of coordinates,
       in terms of the original reference alignment,
       that overlap with non-gap positions in the query
    */

    for (i = 0; i < 3; i+=1) {
        realigned[i] = ""; realigned[i] * Abs (reference);
    }

    reference_coordinate = 0;

    if (offset) {
         for (letter_offset = 0; letter_offset < offset;) {
            realigned[0] * (reference[reference_coordinate]);
            if (reference[reference_coordinate] != '-') {
                letter_offset += 1;
            }
            realigned[1] * "-";
            realigned[2] * "-";
            coordinates[reference_coordinate] = -1;
            reduced_alignment[reference_coordinate] = FALSE;
            reference_coordinate += 1;
        }
    }


    alignment_coordinate = 0;

    while (alignment_coordinate < Abs (aligned_reference) && reference_coordinate < Abs (reference)) {
        coordinates[reference_coordinate] = alignment_coordinate;
        if (reference [reference_coordinate] == "-") { // gap in the reference coordinates; add to all
            reduced_alignment[reference_coordinate] = FALSE;
            reference_coordinate += 1;
            realigned[0] * "-";
            realigned[1] * "-";
            realigned[2] * "-";
        } else {
            if (aligned_reference[alignment_coordinate] == "-") { // insert in the reference
                realigned [0] * "-";
                realigned [1] * "-";
                realigned [2] * aligned_qry[alignment_coordinate];
                alignment_coordinate += 1;
            } else {
                assert ((aligned_reference[alignment_coordinate] &&1) == (reference [reference_coordinate] && 1), "Mismatch between reference and aligned_reference : '`reference [reference_coordinate]`' != '`aligned_reference[alignment_coordinate]`'");
                realigned [0] * reference [reference_coordinate];
                realigned [1] * aligned_reference[alignment_coordinate];
                realigned [2] * aligned_qry[alignment_coordinate];
                reduced_alignment[reference_coordinate] = TRUE;
                alignment_coordinate += 1;
                reference_coordinate += 1;
           }
        }
    }

    while (alignment_coordinate < Abs (aligned_reference)) { // qry is longer than the reference pa
        coordinates[reference_coordinate] = alignment_coordinate;
        realigned [0] * "-";
        realigned [1] * "-";
        realigned [2] * aligned_qry[alignment_coordinate];
        alignment_coordinate += 1;
    }

    while (reference_coordinate < Abs (reference)) { // qry is longer than the reference pa
        coordinates[reference_coordinate] = alignment_coordinate;
        realigned [0] * reference[reference_coordinate];
        realigned [1] * "-";
        realigned [2] * "-";
        reference_coordinate += 1;
    }


    for (i = 0; i < 3; i+=1) {
        realigned[i] * 0;
    }

    return {utility.getGlobalValue("terms.three_way") : Eval (realigned), utility.getGlobalValue("terms.code.mapping") : Eval (coordinates), utility.getGlobalValue("terms.reduced") : Eval (reduced_alignment)};

}

/**
 * @name alignments.Extract_site_patterns
 * @param {DataFilter} data_filter
 * @returns for a data filter, returns a dictionary like this
 *
 *    "pattern id":
 *        "sites" : sites (0-based) mapping to this pattern
 *        "is_constant" : T/F (is the site constant w/matching ambigs)
 *
 *        "0":{
 *           "sites":{
 *             "0":0
 *            },
 *           "is_constant":0
 *          },
 *
 *         "1":{
 *           "sites":{
 *             "0":1
 *            },
 *           "is_constant":1
 *          },
 *...
 *         "34":{
 *           "sites":{
 *             "0":34,
 *             "1":113
 *            },
 *           "is_constant":1
 *          },
 *        ...
 */
lfunction alignments.Extract_site_patterns (data_filter) {
    utility.ToggleEnvVariable ("COUNT_GAPS_IN_FREQUENCIES", FALSE);

    site_info = {};
    GetDataInfo (pattern_list, ^data_filter);
    site_characters = {};
    sequence_count = ^(data_filter + ".species");

    utility.ForEachPair (pattern_list, "_site_index_", "_pattern_",
        '
        utility.EnsureKey (`&site_info`, _pattern_);
        utility.EnsureKey (`&site_info`[_pattern_], utility.getGlobalValue("terms.data.sites"));

        (`&site_info`[_pattern_])[utility.getGlobalValue("terms.data.sites")] + _site_index_[1];

        if (Abs ((`&site_info`[_pattern_])[utility.getGlobalValue("terms.data.sites")]) == 1) {
            // first time we see this site
            GetDataInfo (`&site_characters`, `data_filter`, -1, _pattern_);
            `&site_characters` = utility.Filter (`&site_characters`,
                                                 "_value_",
                                                 "(+_value_>0)");

            (`&site_info`[_pattern_])[utility.getGlobalValue("terms.data.is_constant")] = Abs (`&site_characters`) <= 1;

        }
        '
    );

    utility.ToggleEnvVariable ("COUNT_GAPS_IN_FREQUENCIES", None);

    return site_info;

}
