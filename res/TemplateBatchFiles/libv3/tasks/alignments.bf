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
 * Get i-th sequence name/value from an alignment; retrieve the original sequence name if available
 * @name alignments.GetIthSequenceOriginalName
 * @param {String} dataset_name - name of dataset to get sequence names from
 * @param {String} index - the name of the sequence to extract or None to set up the initial mapping
 * @returns {Dict} {"id" : sequence name, "sequence" : sequence data}
 */

lfunction alignments.GetIthSequenceOriginalName (dataset_name, index) {
    GetString   (seq_id, ^dataset_name, index);
    if (Type (^"`dataset_name`.mapping") == "AssociativeList") {
        seq_id = (^"`dataset_name`.mapping")[seq_id];
    }
    GetDataInfo (seq_string, ^dataset_name, index);
    return {utility.getGlobalValue("terms.id") : seq_id, utility.getGlobalValue("terms.data.sequence") : seq_string};
}

/**
 * Get all sequences as "id" : "sequence" dictionary
 * @name alignments.GetAllSequences
 * @param {String} dataset_name - name of dataset to get sequence names from
 * @returns {Dict} { id -> sequence}
 */

lfunction alignments.GetAllSequences (dataset_name) {

    GetString   (seq_id, ^dataset_name, -1);
    result = {};

    utility.ForEachPair (seq_id, "_index_", "_name_",
    '
         GetDataInfo (`&seq_string`, ^`&dataset_name`, _index_[1]);
        `&result`[_name_] = `&seq_string`;
    ');

    return result;
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
 * Determine what type of data are in the filter
 * @name alignments.FilterType
 * @param {String} filter_name - the name of the dataset filter object
 * @returns {String} one of the standard datatypes, or the alphabet string if it's non-standard
 */

lfunction alignments.FilterType(filter_name) {
    GetDataInfo (parameters, ^filter_name, "PARAMETERS");
    parameters = parameters["ATOM_SIZE"];

    if (parameters == 1) {
        GetDataInfo (alphabet, ^filter_name, "CHARACTERS");
        type = alignments.AlphabetType (alphabet);
        if (null != type) {
            return type;
        }
    } else {
        DataSetFilter throwaway = CreateFilter (^filter_name, 1);
        size = parameters;
        type =  alignments.FilterType(&throwaway);
        DeleteObject (throwaway);
        if (type == ^"terms.nucleotide") {
            if (size == 2) {
                return ^"terms.dinucleotide";
            }
            if (size == 3) {
               return ^"terms.codon";
            }
        }
        GetDataInfo (alphabet, ^filter_name, "CHARACTERS");
    }
    return Join ("", alphabet);
}

/**
 * Read a protein dataset from file_name
 * @name alignments.ReadProteinDataSet
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
        } else {
            if (alphabet == "01") {
                return utility.getGlobalValue ("terms.binary");
            }
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
    if (^"`datafilter_name`.sites"*3 != ^"`dataset_name`.sites") {
        // generate a more diagnostic error here
        for (alignments.LoadCodonDataFile.i = 0; alignments.LoadCodonDataFile.i < ^"`dataset_name`.species"; alignments.LoadCodonDataFile.i += 1) {
            DataSetFilter ^ datafilter_name = CreateFilter( ^ dataset_name, 3,  , "" + alignments.LoadCodonDataFile.i , data_info[terms.stop_codons]);
            if (^"`datafilter_name`.sites"*3 != ^"`dataset_name`.sites") {
                alignments.LoadCodonDataFile.name = alignments.GetIthSequenceOriginalName (dataset_name, alignments.LoadCodonDataFile.i);

                alignments.LoadCodonDataFile.site_map = ^"`datafilter_name`.site_map";

                alignments.LoadCodonDataFile.annotation_string = utility.PopulateDict (0, ^"`dataset_name`.sites",
                                                       '(alignments.LoadCodonDataFile.name[terms.data.sequence])[_idx]',
                                                       '_idx');


                utility.ForEach (alignments.LoadCodonDataFile.site_map, "_value_",
                    '
                        `&alignments.LoadCodonDataFile.annotation_string`[_value_] = `&alignments.LoadCodonDataFile.annotation_string`[_value_] && 0;
                    ');

                console.log ("*** PROBLEM WITH SEQUENCE ' " + alignments.LoadCodonDataFile.name[terms.id] + "' (" + ^"`dataset_name`.sites" + " nt long, stop codons shown in capital letters)\n\n" + Join ("",alignments.LoadCodonDataFile.annotation_string));
                break;
            }
        }
        io.CheckAssertion("`datafilter_name`.sites*3==`dataset_name`.sites", "The input alignment must have the number of sites that is divisible by 3 and must not contain stop codons");
    }
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
 * Take an input filter and check if it has duplicate sequences
 * @name alignments.HasDuplicateSequences
 * @param {String} filter_in - The name of an existing filter
 * @param {Number} check_mode -
      -1 : exact match
      -2 : exact match + gaps match non-gaps");
      -3 : superset filtering
      -4 : partial match filtering
 - @returns the number of *duplicate* sequences
 */
lfunction alignments.HasDuplicateSequences (filter_in, check_mode) {
  GetDataInfo (duplicate_info, ^filter_in, check_mode);
  return duplicate_info["UNIQUE_SEQUENCES"];
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
                                _seq_name_ += "_" + ((`&duplicate_info`)["UNIQUE_COUNTS"])[_idx_[1]];
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

            io.CheckAssertion("`&diff` == 0", "Partition " + this_filter[utility.getGlobalValue("terms.data.name")] + " is either has stop codons or is not in frame");

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
    //console.log (sequence);

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
 * @name alignments.TranslateCodonsToAminoAcidsWithAmbigsAllFrames
 * Translate a codon sequence to amino-acids using the mapping provided by the
 * genetic code in all 3 frames
 * @param {String} sequence - the string to translate
 * @param {Dictionary} code - genetic code description (e.g. returned by alignments.LoadGeneticCode)
 * @param {lookup} code - resolution lookup dictionary
 * @returns {Dict} for each reading frame F in {0, 1, 2} returns

        F -> {
            terms.data.sequence: translated sequence (always choose X if available, otherwise first sense resolution)
            terms.sense_codons : N, // number of sense A/A
            terms.stop_codons : N, // number of stop codons
            terms.terminal_stop : T/F // true if there is a terminal stop codon
        }
 */

 lfunction alignments.TranslateCodonsToAminoAcidsWithAmbigsAllFrames (sequence, code, lookup) {

    result = {};


    for (frame = 0; frame < 3; frame += 1) {
        try_run = alignments.TranslateCodonsToAminoAcidsWithAmbiguities (sequence, frame, code, lookup);

        translation = ""; translation * 128;

        frame_result = {utility.getGlobalValue ("terms.sense_codons") : 0,
                        utility.getGlobalValue ("terms.stop_codons") : 0,
                        utility.getGlobalValue ("terms.terminal_stop") : FALSE
                        };

        upper_bound = Abs (try_run);
        for (i = 0; i < upper_bound; i+=1) {
            if (try_run[i] / "X") { // has_stop
                translation * "X";
                frame_result [^"terms.stop_codons"] += 1;
                if (i == upper_bound - 1) {
                    frame_result [^"terms.terminal_stop"] = TRUE;
                }
            } else {
                translation * (try_run[i])["INDEXORDER"][0];
                frame_result [^"terms.sense_codons"] += 1;
            }
        }


        translation * 0;
        frame_result [utility.getGlobalValue ("terms.data.sequence")] = translation;
        result[frame] = frame_result;
    }

	return result;
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

/**
 * @name alignments.StripGaps
 * Remove gaps from a sequence
 * @param {String} sequence - the input sequence
 * @returns {String} the sequence with all gaps removed
 */

lfunction alignments.StripGaps (sequence) {
    return sequence ^ {{"\-",""}};
}

/**
 * @name alignments.Strip
 * strip non-letters from start and end
 * @param {String} sequence - the input sequence
 * @returns {String} the sequence with all gaps removed
 */

lfunction alignments.Strip (sequence) {
    sequence = sequence ^ {{"[^a-zA-Z]+$",""}};
    return sequence ^ {{"^[^a-zA-Z]+",""}}
}

/**
 * @name alignments.alignment.MapCodonsToAA
 * Map in-frame nucleotides onto a protein alignment string

 * @param {String} codon_sequence - the codon sequence to map
 * @param {String} aa_sequence - the matching aligned a.a. sequence
 * @param {Number} no more than this many mismatches - the codon sequence to map
 * @param {Dict} mapping - code ["terms.code.mapping"]

 * @returns {String} the mapped sequence

 * @example
    GCAAAATCATTAGGGACTATGGAAAACAGA
    -AKSLGTMEN-R

    maps to

    ---GCAAAATCATTAGGGACTATGGAAAAC---AGA

 */

lfunction alignment.MapCodonsToAA(codon_sequence, aa_sequence, this_many_mm, mapping) {

    seqLen = Abs(aa_sequence);
    translString = "";
    translString * (seqLen);
    seqLenN = Abs(codon_sequence);

    aaPos = 0;
    seqPos = 0;
    codon = codon_sequence[seqPos][seqPos + 2];
    currentAA = mapping[codon];
    if (currentAA == 0) {
        currentAA = "X";
    }

    mismatch_count = 0;

    for (aaPos = 0; aaPos < seqLen && seqPos < seqLenN; aaPos += 1) {
        advance = 1;
        copy_codon = 1;

        if (currentAA != 0) {
            if (aa_sequence[aaPos] == "-") {
                //if (currentAA != "X") {
                    translString * "---";
                    advance = 0;
                //}
            } else {
                mismatch_count += (aa_sequence[aaPos] != currentAA && currentAA != "X");
                if (this_many_mm == 1) {
                    if (mismatch_count == this_many_mm) {
                        translString * 0;
                        console.log (translString);
                        console.log ("\n");
                        console.log (codon_sequence);
                        console.log ("\n");
                        console.log (aa_sequence);
                        console.log ("\n");
                    }
                    assert(mismatch_count < this_many_mm, "A mismatch between codon and protein sequences at position " + aaPos + " (codon `seqPos`) : codon '" + codon_sequence[seqPos][seqPos + 2] + "'(`currentAA`) a.a. '`aa_sequence[aaPos]`'");
                } else {
                    if (mismatch_count >= this_many_mm) {
                        translString * 0;
                        return None;
                    }
                }
            }
        } else {
            copy_codon = 0;
        }

        if (advance) {
            if (copy_codon) {
                if (currentAA == "X" && mapping[codon] == "X") {
                    translString * "---";
                } else {
                    translString * codon;
                }
            } else {
                //fprintf (stdout, "Skipping codon ", codon, "\n");
                //aaPos = aaPos - 1;
            }
            seqPos += 3;
            codon = codon_sequence[seqPos][seqPos + 2];
            currentAA = mapping[codon];
            if (currentAA == 0) {
                currentAA = "X";
            }
        }
    }

    for (; aaPos < seqLen; aaPos += 1) {
        translString * "---";
    }


    translString * 0;
    return translString;
}

/**
 * @name alignment.ExportPartitionedNEXUS
 * Export a datafilter with partitions and trees to a file

 * @param {String} filter - filter name
 * @param {Matrix} breakPoints  - locations of breakpoints (Nx1)
 * @param {Matrix} trees  - tree strings for partitions (Nx1)
 * @param {String} file - write the result here
 * @param {Bool}   isCodon - is the filter a codon filter?

 * @returns {String} the mapped sequence

 * @example
    GCAAAATCATTAGGGACTATGGAAAACAGA
    -AKSLGTMEN-R

    maps to

    ---GCAAAATCATTAGGGACTATGGAAAAC---AGA

 */

lfunction alignment.ExportPartitionedNEXUS (filter, breakPoints, trees, file, isCodon) {
    utility.ToggleEnvVariable ("DATA_FILE_PRINT_FORMAT", 4);

    fprintf (file, CLEAR_FILE, KEEP_OPEN, ^filter, "\n");

    breakPointsCount = utility.Array1D (breakPoints);
    partCount      = breakPointsCount + 1;
    currentStart    = 0;

    fprintf (file, "\nBEGIN ASSUMPTIONS;\n");

    for (p = 0; p < partCount; p += 1) {
        lastPartition = p >= breakPointsCount;
        if (!lastPartition) {
            currentEnd = breakPoints[p];
        } else {
            if (isCodon) {
                currentEnd = ^(filter+".sites") * 3;
            } else {
                currentEnd = ^(filter+".sites");
            }
            currentEnd = currentEnd - 1;
        }

        fprintf (file, "\tCHARSET span_", p + 1, " = ", currentStart + 1, "-", currentEnd + 1, ";\n");


        if (!lastPartition) {
            currentStart = breakPoints[p] + 1;
        }
    }

    fprintf (file, "END;\nBEGIN TREES;\n");
    for (p = 0; p < partCount; p += 1) {
        fprintf (file, "\tTREE tree_", p + 1, " = ", trees[p], ";\n");
    }

    fprintf (file, "END;\n");
    utility.ToggleEnvVariable ("DATA_FILE_PRINT_FORMAT", None);
}
