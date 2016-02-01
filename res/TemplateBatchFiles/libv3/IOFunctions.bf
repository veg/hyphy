LoadFunctionLibrary("ReadDelimitedFiles");

function io.readCodonDataSet(dataset_name) {
    return io.readCodonDataSetFromPath(dataset_name, None);
}

function io.readCodonDataSetFromPath(dataset_name, path) {
    ExecuteAFile(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TemplateModels" + DIRECTORY_SEPARATOR + "chooseGeneticCode.def");

    if (Type(path) == "String") {
        ExecuteCommands("DataSet `dataset_name` = ReadDataFile (`path`);");
    } else {
        ExecuteCommands("DataSet `dataset_name` = ReadDataFile (PROMPT_FOR_FILE);");
        path = LAST_FILE_PATH;
    }
    return {
        "code": _Genetic_Code,
        "stop": GeneticCodeExclusions,
        "file": path,
        "sequences": Eval("`dataset_name`.species")
    };
}

function io.readNucleotideDataSet_aux (dataset_name) {
    return {
        "sequences": Eval("`dataset_name`.species"),
        "sites": Eval("`dataset_name`.sites"),
        "name-mapping": Eval("`dataset_name`.mapping")
    };
}

function io.readNucleotideDataSet(dataset_name, file_name) {
    if (Type(file_name) == "String") {
        ExecuteCommands("DataSet `dataset_name` = ReadDataFile (`file_name`);");
    } else {
        ExecuteCommands("DataSet `dataset_name` = ReadDataFile (PROMPT_FOR_FILE);");
        file_name = LAST_FILE_PATH;
    }

    io.readNucleotideDataSet.result =  io.readNucleotideDataSet_aux (dataset_name);
    io.readNucleotideDataSet.result["file"] = file_name;
    return io.readNucleotideDataSet.result;
}

function io.readNucleotideDataSetString(dataset_name, data) {
    ExecuteCommands("DataSet `dataset_name` = ReadFromString (data);");
    return io.readNucleotideDataSet_aux (dataset_name);
}

function io.getTreeString._sanitize(string) {
    if (_DO_TREE_REBALANCE_) {
        string = RerootTree(string, 0);
    }

    if (_KEEP_I_LABELS_) {
        utility.toggleEnvVariable("INTERNAL_NODE_PREFIX", "intNode");
    }
    string = string ^ {
        {
            "\\)[0-9]+(\\.[0-9]*)?\:", "):"
        }
    };

    if (_KEEP_I_LABELS_) {
        utility.toggleEnvVariable("INTERNAL_NODE_PREFIX", None);
    }
    return string;
}

function io.getTreeString(look_for_newick_tree) {

    UseModel(USE_NO_MODEL);

    if (look_for_newick_tree == 0) {
        IS_TREE_PRESENT_IN_DATA = 0;
    }

    if (IS_TREE_PRESENT_IN_DATA) {
        fprintf(stdout, "\n> A tree was found in the data file: ``", DATAFILE_TREE, "``\n>Would you like to use it? ");
        fscanf(stdin, "String", io.getTreeString.response);
        if (io.getTreeString.response == "n" || io.getTreeString.response == "N") {
            IS_TREE_PRESENT_IN_DATA = 0;
        } else {
            io.getTreeString.treeString = io.getTreeString._sanitize(DATAFILE_TREE);
            IS_TREE_PRESENT_IN_DATA = 1;
        }
        fprintf(stdout, "\n\n");
    }

    if (!IS_TREE_PRESENT_IN_DATA) {
        SetDialogPrompt("Please select a tree file for the data:");
        fscanf(PROMPT_FOR_FILE, REWIND, "Raw", io.getTreeString.treeString);
        fprintf(stdout, "\n");
            
        if (regexp.find(io.getTreeString.treeString, "^#NEXUS")) {
            ExecuteCommands(io.getTreeString.treeString);
            if (IS_TREE_PRESENT_IN_DATA == 0) {
                fprintf(stdout, "\n> **This NEXUS file doesn't contain a valid tree block**");
                return 1;
            }
            if (Rows(NEXUS_FILE_TREE_MATRIX) > 1) {
                ChoiceList(io.getTreeString.treeChoice, "Select a tree", 1, SKIP_NONE, NEXUS_FILE_TREE_MATRIX);
                if (io.getTreeString.treeChoice < 0) {
                    return 1;
                }
                io.getTreeString.treeString = NEXUS_FILE_TREE_MATRIX[io.getTreeString.treeChoice][1];
            } else {
                io.getTreeString.treeString = NEXUS_FILE_TREE_MATRIX[0][1];
            }
        } else {
            io.getTreeString.start = (io.getTreeString.treeString $ "\\(")[0];
            if (io.getTreeString.start < 0) {
                fprintf(stdout, "\n> **This doesn't seem to be a valid Newick string file**. Can't find the opening parenthesis. ``", io.getTreeString, "``\n");
                return 1;
            } else {
                io.getTreeString.parenCounter = 1;
                io.getTreeString.current = io.getTreeString.start + 1;
                while (io.getTreeString.current < Abs(io.getTreeString.treeString) && io.getTreeString.parenCounter) {
                    io.getTreeString.char = io.getTreeString.treeString[io.getTreeString.current];
                    if (io.getTreeString.char == "(") {
                        io.getTreeString.parenCounter += 1;
                    } else {
                        if (io.getTreeString.char == ")") {
                            io.getTreeString.parenCounter += (-1);
                        }
                    }
                    io.getTreeString.current += 1;
                }

                if (io.getTreeString.parenCounter) {
                    fprintf(stdout, "\n> ** This doesn't seem to be a valid Newick string file**. Can't match the parentheses. \n``", io.getTreeString.treeString, "``\n");
                    return 1;
                }

                io.getTreeString.treeString = io.getTreeString.treeString[io.getTreeString.start][io.getTreeString.current - 1];
            }

            io.getTreeString.treeString = io.getTreeString._sanitize(io.getTreeString.treeString);

        }
    }

    return io.getTreeString.treeString;
}

function io.checkAssertion(statement, error_msg) {
    ExecuteCommands("assert (`statement`, error_msg)");
    return None;
}

lfunction io._reportMessageHelper(analysis, text) {
    if (Abs(analysis)) {
        return "[`analysis`] `text`";
    } else {
        return text;
    }
}

lfunction io.spool_json(json, file) {
    utility.toggleEnvVariable("USE_JSON_FOR_MATRIX", 1);
    if (Type(file) == "String") {
        fprintf(file, CLEAR_FILE, json);
    } else {
        fprintf(stdout, "\n", json, "\n");
    }
    utility.toggleEnvVariable("USE_JSON_FOR_MATRIX", None);
}


lfunction io.reportProgressMessage(analysis, text) {
    fprintf(stdout, io._reportMessageHelper(analysis, text), "\n");
}

lfunction io.reportProgressMessageMD(analysis, stage, text) {
    if (Abs(cache) == 0) {
        cache = {};
    }
    advance = TRUE;
    if (Abs(cache[analysis])) {
        if ((cache[analysis])[stage]) {
            advance = FALSE;
        }
        (cache[analysis])[stage] += 1;
    } else {
        cache[analysis] = {};
        (cache[analysis])[stage] = 1;
    }

    if (advance) {
        if (Abs (cache[analysis]) == 1) {
            fprintf (stdout, "\n");
        }
        fprintf(stdout, Abs(cache[analysis]), ". ", text, "\n");
    } else {
        fprintf(stdout, "    ", text, "\n");
    }
}

function io.reportProgressBar(analysis, text) {
    SetParameter(STATUS_BAR_STATUS_STRING, io._reportMessageHelper(analysis, text), 0);
}

function io.validate_a_list_of_files(list) {
    io.validate_a_list_of_files.result = {};
    for (io.validate_a_list_of_files.i = 0; io.validate_a_list_of_files.i < Rows(list) * Columns(list); io.validate_a_list_of_files.i += 1) {
        if (Abs(list[io.validate_a_list_of_files.i])) {
            io.validate_a_list_of_files.fn = list[io.validate_a_list_of_files.i];
            io.checkAssertion("!io.validate_a_list_of_files.fn", "HyPhy cannot open '" + io.validate_a_list_of_files.fn + "' for reading");
            io.validate_a_list_of_files.result + io.validate_a_list_of_files.fn;
        }
    }
    return io.validate_a_list_of_files.result;
}


function io.get_a_list_of_files(filename) {
    if (Type(filename) == "String") {
        if (!filename) { // filename exists
            fscanf(filename, REWIND, "Lines", io.get_a_list_of_files.list);
            return io.validate_a_list_of_files(io.get_a_list_of_files.list);
        }
    }

    io.get_a_list_of_files.result = {};
    io.printAndUnderline("Enter paths to files (blank line to end entry)", "-");
    while (1) {
        fprintf(stdout, "* File ", Abs(io.get_a_list_of_files.result) + 1, " [relative path `PATH_TO_CURRENT_BF`]:");
        io.get_a_list_of_files.current_path = "";
        fscanf(stdin, "String", io.get_a_list_of_files.current_path);
        if (Abs(io.get_a_list_of_files.current_path)) {
            io.checkAssertion("! io.get_a_list_of_files.current_path", "HyPhy cannot open '" + io.get_a_list_of_files.current_path + "' for reading");
        } else {
            break;
        }
        io.get_a_list_of_files.result + io.get_a_list_of_files.current_path;
    }
}

function io.displayAnalysisBanner(analysis_info) {
    if (io.hasStringKey("info", analysis_info)) {
        io.printAndUnderline("Analysis Description", "-");
        fprintf(stdout, io.formatLongStringToWidth(analysis_info["info"], 72), "\n");
    }
    if (io.hasStringKey("requirements", analysis_info)) {
        fprintf(stdout, "\n- __Requirements__: ");
        fprintf(stdout, io.formatLongStringToWidth(analysis_info["requirements"], 72), "\n");
    }
    if (io.hasStringKey("reference", analysis_info)) {
        fprintf(stdout, "\n- __Citation__: ");
        fprintf(stdout, io.formatLongStringToWidth(analysis_info["reference"], 72), "\n");
    }
    if (io.hasStringKey("authors", analysis_info)) {
        fprintf(stdout, "\n- __Written by__: ");
        fprintf(stdout, io.formatLongStringToWidth(analysis_info["authors"], 72), "\n");
    }
    if (io.hasStringKey("contact", analysis_info)) {
        fprintf(stdout, "\n- __Contact Information__: ");
        fprintf(stdout, io.formatLongStringToWidth(analysis_info["contact"], 72), "\n");
    }
    if (io.hasStringKey("version", analysis_info)) {
        fprintf(stdout, "\n- __Analysis Version__: ");
        fprintf(stdout, io.formatLongStringToWidth(analysis_info["version"], 72), "\n");
    }
    fprintf(stdout, "\n");

    return None;
}

function io.hasStringKey(key, dict) {
    return Type(dict[key]) == "String";
}

function io.spoolLF(lf_id, trunk_path, tag) {
    ExecuteCommands("Export (__lf_spool, `lf_id`);");
    if (tag == None) {
        tag = lf_id;
    }
    fprintf(trunk_path + "." + tag + ".bf", CLEAR_FILE, __lf_spool);
}

function io.printAndUnderline(string, char) {
    fprintf(stdout, "\n", string, "\n");
    for (k = 0; k < Abs(string); k += 1) {
        fprintf(stdout, char[0]);
    }
    fprintf(stdout, "\n");
}

function io.formatLongStringToWidth(string, width) {
    words = splitOnRegExp(string, "[\\ \n]+");
    lines = {};

    current_line = "";
    words_in_current_line = 0;
    for (i = 0; i < Abs(words); i += 1) {
        if (words_in_current_line == 0) {
            current_line = words[i];
            words_in_current_line = 1;
        } else {
            if (Abs(current_line) + Abs(words[i]) + 1 <= width) {
                words_in_current_line += 1;
                current_line += " " + words[i];
            } else {
                lines + current_line;
                words_in_current_line = 0;
                current_line = "";
                i = i - 1;
            }
        }
    }

    if (words_in_current_line) {
        lines + current_line;
    }

    return Join("\n", lines);
}