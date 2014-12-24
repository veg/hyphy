LoadFunctionLibrary ("ReadDelimitedFiles");

function io.readCodonDataSet (dataset_name) {
    ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" 
                                      + DIRECTORY_SEPARATOR 
                                      + "TemplateModels" 
                                      + DIRECTORY_SEPARATOR 
                                      + "chooseGeneticCode.def");
    
    ExecuteCommands ("DataSet `dataset_name` = ReadDataFile (PROMPT_FOR_FILE);");
    return {"code": _Genetic_Code, "stop" : GeneticCodeExclusions, "file" : LAST_FILE_PATH, "sequences" : Eval ("`dataset_name`.species")};
}

function io.getTreeString (look_for_newick_tree) {

    UseModel (USE_NO_MODEL);

    if(look_for_newick_tree == 0) {
      IS_TREE_PRESENT_IN_DATA = 0;
    }

    ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" 
                                      + DIRECTORY_SEPARATOR 
                                      + "queryTree.bf");
    return treeString;
}

function io.checkAssertion (statement, error_msg) {
    ExecuteCommands ("assert (`statement`, error_msg)");
    return None;
}

function io.reportProgressMessage (analysis, text) {
    if (Abs (analysis)) {
        fprintf (stdout, "[`analysis`] `text` \n");
    } else {
        fprintf (stdout, text, "\n");
    }
    return None;
}

function io.reportProgressMessage (analysis, text) {
    fprintf (stdout, "[`analysis`] `text` \n");
    return None;
}

function io.displayAnalysisBanner (analysis_info) {
    if (io.hasStringKey ("info", analysis_info)) {
        io.printAndUnderline ("Analysis Description", "=");
        fprintf (stdout, io.formatLongStringToWidth (analysis_info["info"], 72), "\n");
    }
    if (io.hasStringKey ("requirements", analysis_info)) {
        io.printAndUnderline ("Requirements", "=");
        fprintf (stdout, io.formatLongStringToWidth (analysis_info["requirements"], 72), "\n");
    }
    if (io.hasStringKey ("reference", analysis_info)) {
        io.printAndUnderline ("Citation", "=");
        fprintf (stdout, io.formatLongStringToWidth (analysis_info["reference"], 72), "\n");
    }
    if (io.hasStringKey ("authors", analysis_info)) {
        io.printAndUnderline ("Written by", "=");
        fprintf (stdout, io.formatLongStringToWidth (analysis_info["authors"], 72), "\n");
    }   
    if (io.hasStringKey ("contact", analysis_info)) {
        io.printAndUnderline ("Contact Information", "=");
        fprintf (stdout, io.formatLongStringToWidth (analysis_info["contact"], 72), "\n");
    }   
    if (io.hasStringKey ("version", analysis_info)) {
        io.printAndUnderline ("Analysis Version", "=");
        fprintf (stdout, io.formatLongStringToWidth (analysis_info["version"], 72), "\n");
    }   
    fprintf (stdout, "\n");
   
    return None;
}

function io.hasStringKey (key, dict) {
    return Type (dict[key]) == "String";
}

function io.spoolLF (lf_id, trunk_path, tag) {
    ExecuteCommands ("Export (__lf_spool, `lf_id`);");
    if (tag == None) {
        tag = lf_id;
    }
    fprintf (trunk_path + "." + tag + ".bf", CLEAR_FILE, __lf_spool);
}

function io.printAndUnderline (string, char) {
    fprintf (stdout, "\n", string, "\n");
    for (k = 0; k < Abs (string); k+=1) {
        fprintf (stdout, char[0]);
    }
    fprintf (stdout, "\n");
}

function io.formatLongStringToWidth (string, width) {
    words = splitOnRegExp (string, "[\\ \n]+");
    lines = {};
    
    current_line = "";
    words_in_current_line = 0;
    for (i = 0; i < Abs (words); i+=1 ) {
        if (words_in_current_line == 0) {
            current_line = words[i];
            words_in_current_line = 1;
        } else {
            if (Abs (current_line) + Abs (words[i]) + 1 <= width) {
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
    
    return Join ("\n", lines);
}
