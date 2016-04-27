LoadFunctionLibrary ("../IOFunctions.bf");
LoadFunctionLibrary("../terms-json.bf");


lfunction alignments.readCodonDataSet(dataset_name) {
    return alignments.readCodonDataSetFromPath(dataset_name, None);
}

lfunction alignments.readCodonDataSetFromPath(dataset_name, path) {
    ExecuteAFile(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TemplateModels" + DIRECTORY_SEPARATOR + "chooseGeneticCode.def");
    
    if (Type(path) == "String") {
        DataSet ^dataset_name = ReadDataFile (path);
    } else {
        DataSet ^dataset_name = ReadDataFile (PROMPT_FOR_FILE);
        path = ^"LAST_FILE_PATH";
    }
         
    r =  alignments.readNucleotideDataSet_aux (dataset_name);  
    
    r * {
        "code": _Genetic_Code,
        "stop": GeneticCodeExclusions,
        "file": path,
        "sequences": Eval(^"`dataset_name`.species")
    };
    
    
    return r;
}

lfunction alignments.readNucleotideDataSet_aux (dataset_name) {

    partitions = None;
    dfpm = utility.getEnvVariable ("DATA_FILE_PARTITION_MATRIX");
        
    partitions = {{"default",""}};
    
         
    if (Type (dfpm) == "Matrix") {
        if (Rows (dfpm) > 0) {
            partitions = Transpose(dfpm);
        }
    } 

    return {
        "sequences": Eval("`dataset_name`.species"),
        "sites": Eval("`dataset_name`.sites"),
        "name-mapping": Eval("`dataset_name`.mapping"), 
        "partitions" : partitions
    };
}

lfunction alignments.getSequenceNames (dataset_name) {
    GetString (result, ^dataset_name, -1);
    return result;
}

// Creates random unique namespace
lfunction alignments.readNucleotideDataSet(dataset_name, file_name) {
    if (Type(file_name) == "String") {
        // 
        DataSet ^dataset_name = ReadDataFile (^file_name);
    } else {
        DataSet ^dataset_name = ReadDataFile (PROMPT_FOR_FILE);
        file_name = LAST_FILE_PATH;
    }

    result =  alignments.readNucleotideDataSet_aux (dataset_name);
    result["file"] = file_name;
    return result;
}

function alignments.readNucleotideDataSetString(dataset_name, data) {
    ExecuteCommands("DataSet `dataset_name` = ReadFromString (data);");
    return alignments.readNucleotideDataSet_aux (dataset_name);
}

function alignments.promptForGeneticCodeAndAlignment (dataset_name, datafilter_name) {
    return alignments.loadCodonDataFile (dataset_name, datafilter_name, alignments.readCodonDataSet (dataset_name));
}

function alignments.loadGeneticCodeAndAlignment (dataset_name, datafilter_name, path) {
    return alignments.loadCodonDataFile (dataset_name, datafilter_name, alignments.readCodonDataSetFromPath (dataset_name, path));
}

function alignments.loadCodonDataFile (dataset_name, datafilter_name, data_info) {
    DataSetFilter ^datafilter_name = CreateFilter (^dataset_name, 3, , , data_info["stop"]);
    io.checkAssertion ("`datafilter_name`.sites*3==`dataset_name`.sites","The input alignment must not contain stop codons");
    data_info ["sites"] =  ^"`datafilter_name`.sites";
    data_info ["dataset"] = dataset_name;
    data_info ["datafilter"] = datafilter_name;
    return data_info;
}

function alignments.readNucleotideAlignment (file_name, dataset_name, datafilter_name) {
    data_info = alignments.readNucleotideDataSet (dataset_name, file_name);
    ExecuteCommands ("DataSetFilter `datafilter_name` = CreateFilter (`dataset_name`,1)");
    data_info ["sites"] = Eval ("`datafilter_name`.sites");
    data_info ["dataset"] = dataset_name;
    data_info ["datafilter"] = datafilter_name;

    return data_info;
}

lfunction alignments.defineFiltersForPartitions (partitions, source_data, prefix,  data_info) {
    part_count = utility.array1D (partitions);
    filters    = {};
    if (utility.checkKey (data_info, "code", "Matrix")) {
        for (i = 0; i < part_count; i+=1) {
            this_filter = {};
            DataSetFilter test = CreateFilter (^source_data, 1, (partitions[i])["filter-string"]);
            this_filter["name"] = prefix + (partitions[i])["name"];
            DataSetFilter ^(this_filter["name"]) = CreateFilter (^source_data, 3, (partitions[i])["filter-string"],,data_info["stop"]);
            diff = test.sites - 3*^(this_filter["name"]+".sites");
            io.checkAssertion ("`&diff` == 0", "Partition " + (filters["names"])[i] + " is either has stop codons or is not in frame");
            this_filter["coverage"] = utility.dict_to_array (utility.map (utility.filter (^(this_filter["name"]+".site_map"), "_value_", "_value_%3==0"), "_value_", "_value_$3"));
            
            filters + this_filter;
        }
                
    } else {
        for (i = 0; i < part_count; i+=1) {
            this_filter = {};
            this_filter["name"] = prefix + (partitions[i])["name"];
            DataSetFilter ^(this_filter["name"]) = CreateFilter (^source_data, 1, (partitions[i])["filter-string"]);
            this_filter["coverage"] = ^(this_filter["name"]+".site_map");
            filters + this_filter;

        }
   
    }
    return filters;
}
