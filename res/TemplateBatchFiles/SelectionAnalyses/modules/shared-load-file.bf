function load_file (prefix) {

    codon_data_info = alignments.promptForGeneticCodeAndAlignment(prefix+".codon_data", prefix+".codon_filter");

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

    sample_size = codon_data_info["sites"] * codon_data_info["sequences"];
    codon_data_info["json"] = codon_data_info["file"] + "."+prefix+".json";
    name_mapping = codon_data_info[utility.getGlobalValue("terms.json.name_mapping")];
        /**
            will contain "mapped" -> "original" associations with sequence names; or null if no mapping was necessary
        */

    if (None == name_mapping) { /** create a 1-1 mapping if nothing was done */
        name_mapping = {};
        utility.forEach (alignments.getSequenceNames (prefix+".codon_data"), "_value_", "`&name_mapping`[_value_] = _value_");
    }

    partitions_and_trees = trees.loadAnnotatedTreeTopology.match_partitions (codon_data_info[utility.getGlobalValue("terms.json.partitions")], name_mapping);

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

    partition_count = Abs (partitions_and_trees);

    utility.forEachPair (partitions_and_trees,
                            "_key_",
                            "_value_",
                            '(`&partitions_and_trees`[_key_])["filter-string"] = selection.io.adjust_partition_string (_value_["filter-string"], 3*`&codon_data_info`["sites"])');
        /**
            ensure that all partitions fall on codon boundaries if they are contiguous
        */


    io.reportProgressMessage ("", ">Loaded a multiple sequence alignment with **" + codon_data_info["sequences"] + "** sequences, **" + codon_data_info["sites"] + "** codons, and **" + partition_count + "** partitions from \`" + codon_data_info["file"] + "\`");

    selected_branches = selection.io.defineBranchSets(partitions_and_trees);
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
}
