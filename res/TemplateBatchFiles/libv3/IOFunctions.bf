LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("convenience/regexp.bf");
LoadFunctionLibrary("libv3/all-terms.bf");



/**
 * @name io.CheckAssertion
 * @param statement
 * @param error_msg
 */
function io.CheckAssertion(statement, error_msg) {
    ExecuteCommands("assert (`statement`, error_msg)");
}

/**
 * @name io.PromptUser
 * @param prompt
 * @param default
 * @param lower_bound
 * @param upper_bound
 * @param is_integer
 */
lfunction io.PromptUser(prompt,
    default, lower_bound, upper_bound, is_integer) {
    value = lower_bound - 1;
    max_tries = 10;

    while (value < lower_bound || value > upper_bound) {
        io.CheckAssertion ("`&max_tries`", "Failed to provide a valid value after the maximum number of propmts");

        fprintf(stdout, prompt, " (permissible range = [", lower_bound, ",", upper_bound, "], default value = ",
            default);
        if (is_integer) {
            fprintf(stdout, ", integer");
        }
        fprintf(stdout, "): ");
        fscanf(stdin, "String", str_val);

        if (Abs(str_val) == 0) {
            value = 0 +
                default;
        } else {
            value = 0 + str_val;
        }
        if (is_integer) {
            value = value $ 1;
        }
        max_tries = max_tries - 1;
    }
    return value;
}

/**
 * @name io.PromptUserString
 * @param prompt
 */
lfunction io.PromptUserForString(prompt) {

    str_val = "";

    while (Abs (str_val) == 0) {
        fprintf(stdout, prompt, " : ");
        fscanf(stdin, "String", str_val);
    }
    return str_val;
}

/**
 * @name io.HandleCacheRequest
 * @param prompt
 */
lfunction io.ReadFromOrCreate (prompt, default_value) {
    SetDialogPrompt (prompt);
    contents = None;
    fscanf (PROMPT_FOR_FILE, CREATE_FILE, "Raw", contents);
        
    if (^"FILE_CREATED") {
        contents = default_value;   
    } else {
        if (None != contents) {
            contents =  Eval (contents);
        } else {
            contents = default_value;
        }
    }
    return {
        ^"terms.data.file" : ^"LAST_FILE_PATH",
        ^"terms.data.value" : contents
    };
}

/**
 * @name io.PromptUserForFilePath
 * @param prompt
 */
lfunction io.PromptUserForFilePath(prompt) {
    SetDialogPrompt (prompt);
    fprintf (PROMPT_FOR_FILE, CLEAR_FILE);
    return  ^"LAST_FILE_PATH";
}

/**
 * @name io.PromptUserForFilePath
 * @param prompt
 */
lfunction io.PromptUserForFilePathRead(prompt) {
    SetDialogPrompt (prompt);
    fscanf (PROMPT_FOR_FILE, REWIND, "String", n);
    return  ^"LAST_FILE_PATH";
}

/**
 * @name  io.LoadFile
 * @param prompt
 */
 
function io.LoadFile(prompt) {
   SetDialogPrompt (prompt);
   ExecuteAFile (PROMPT_FOR_FILE);
   return  ^"LAST_FILE_PATH";
}


/**
 * @name io._reportMessageHelper
 * @param analysis
 * @param text
 * @returns formatted message
 */
lfunction io._reportMessageHelper(analysis, text) {
    if (Abs(analysis)) {
        return "[`analysis`] `text`";
    } else {
        return text;
    }
}

/**
 * @name io.CompressedObjectSpool
 * @param json
 * @param file
 */
 
lfunction io.CompressedObjectSpool (object, file) {
    GetString (have_zip, ZLIB_ENABLED, 0);
    if (have_zip && file != "/dev/null") {
        if (utility.GetEnvVariable ("GZIP_OUTPUT")) {
            utility.ToggleEnvVariable("GZIP_COMPRESSION_LEVEL", 5);
            fprintf(file, CLEAR_FILE, object);
            utility.ToggleEnvVariable("GZIP_COMPRESSION_LEVEL", None);
            return;
        }
    }
    fprintf(file, CLEAR_FILE, object);
 }

/**
 * @name io.SpoolJSON
 * @param json
 * @param file
 */
lfunction io.SpoolJSON(json, file) {
    utility.ToggleEnvVariable("USE_JSON_FOR_MATRIX", 1);
    if (Type(file) == "String" && Abs (file) > 0) {
        io.CompressedObjectSpool (json, file);
    } else {
        fprintf(stdout, "\n", json, "\n");
    }
    utility.ToggleEnvVariable("USE_JSON_FOR_MATRIX", None);
}

/**
 * TODO: Does not support arrays in JSON
 * Parses json from file_path
 * @name io.ParseJSON
 * @param file
 */
lfunction io.ParseJSON(file_path) {
    fscanf(file_path, REWIND, "Raw", test);
    parsed_test = Eval(test);
    DeleteObject (test);
    return parsed_test;
}

/**
 * @name io.ReportProgressMessage
 * @param analysis
 * @param text
 */
lfunction io.ReportProgressMessage(analysis, text) {
    fprintf(stdout, io._reportMessageHelper(analysis, text), "\n");
}

/**
 * @name io.ReportAnalysisStageMD
 * @param stage
 */
lfunction io.ReportAnalysisStageMD(stage) {
    fprintf(stdout, "\n>", stage, "\n\n");
}

/**
 * @name io.ReportProgressMessageMD
 * @param analysis
 * @param stage
 * @param text
 */
lfunction io.ReportProgressMessageMD(analysis, stage, text) {
    if (Abs(cache) == 0) {
        cache = {};
    }
    advance = TRUE;
    utility.EnsureKey(cache, analysis);

    if (utility.Has (cache[analysis],stage,"Number")) {
        advance = FALSE;
    }

    (cache[analysis])[stage] += 1;

    if (advance) {
        if (Abs(cache[analysis]) == 1) {
            fprintf(stdout, "\n");
        }
        fprintf(stdout, "\n### ", text, "\n");
    } else {
        fprintf(stdout, text, "\n");
    }
}

/**
 * Reports stats generated from math.GatherDescriptiveStats
 * @name io.ReportStatsMD
 * @param label - Typically the name of the method used
 * @param stats - The stats dict generated by math.GatherDescriptiveStats
 */
lfunction io.ReportStatsMD(_label, _stats) {

    if (Abs (_label)) {
        console.log (_label + "\n");
    }

    _table_output_options = {
        utility.getGlobalValue("terms.table_options.header"): 1,
        utility.getGlobalValue("terms.table_options.minimum_column_width"): 16,
        utility.getGlobalValue("terms.table_options.align"): "center",
        utility.getGlobalValue("terms.table_options.column_widths"): {
            "0": 16,
            "1": 16,
        }
    };

    _results = {
        Abs(_stats),
        2
    };
    _keys = utility.Keys(_stats);

    for (_k = 0; _k < Abs(_stats); _k = _k + 1) {
        _results[_k][0] = _keys[_k];
        _results[_k][1] = _stats[_keys[_k]];
    }

    _header = {
        2,
        1
    };
    _header[0] = "Metric";
    _header[1] = "Value";
    fprintf(stdout, io.FormatTableRow(_header, _table_output_options));
    _table_output_options[utility.getGlobalValue("terms.table_options.header")] = FALSE;

    for (_k = 0; _k < Abs(_stats); _k = _k + 1) {
        _tmp_matrix = {
            2,
            1
        };
        _tmp_matrix[0] = _results[_k][0];
        _tmp_matrix[1] = Format(_results[_k][1], 3, 3);
        fprintf(stdout, io.FormatTableRow(_tmp_matrix, _table_output_options));
    }

    return 0;
}

/**
 * @name io.ReportProgressBar
 * @param analysis
 * @param text
 */

io.ReportProgressBar._line_width = 0;

lfunction io.ReportProgressBar(analysis, text) {
    str = io._reportMessageHelper(analysis, text);
    SetParameter(STATUS_BAR_STATUS_STRING, str, 0);
    utility.ExecuteInGlobalNamespace ("io.ReportProgressBar._line_width = " + Abs (str));
}

lfunction io.ClearProgressBar() {
    eraser = ""; eraser * (^"io.ReportProgressBar._line_width" + 1);
    for (k = 0; k < ^"io.ReportProgressBar._line_width"; k+=1) {
        eraser * " ";
    }
    eraser * 0;
    SetParameter(STATUS_BAR_STATUS_STRING, eraser, 0);
    SetParameter(STATUS_BAR_STATUS_STRING, "", 0);
}

/**
 * @name io.validate_a_list_of_files
 * @param list
 */
lfunction io.validate_a_list_of_files(list) {
    result = {};
    dim = utility.Array1D (list);
	base_dir = utility.GetEnvVariable ("HYPHY_BASE_DIRECTORY");


    for (i = 0; i < dim; i += 1) {
        if (Abs(list[i])) {
            fn = list[i];
            if (io.FileExists (fn)) {
            	result + fn;
            	continue;
            } else {
            	fn = base_dir + fn;
            	if (io.FileExists (fn)) {
            		result + fn;
            		continue;
				}
            }
            io.CheckAssertion("io.FileExists(`&fn`)", "HyPhy cannot open '" + fn + "' for reading");

        }
    }
    return result;
}

/**
 * @name io.format_object
 * @param object
 * @param options
 */
lfunction io.format_object(object, options) {

    if (Type(object) == "String") {
        return object;
    }
    if (Type(object) == "Number") {
        if (None != options) {
            if (Abs(options[utility.getGlobalValue("terms.number_precision")]) > 0) {
                return Eval("Format (`&object`, 0, " + options[utility.getGlobalValue("terms.number_precision")] + ")");
            }
        }
    }

    return "" + object;
}

/**
 * @name io.FormatTableRow
 * @param row
 * @param options
 */
lfunction io.FormatTableRowDecorators (row, options, prefix, suffix) {

    if (None == options) {
        options = {};
    }

    cells = utility.Map(row, "_value_", "io.format_object(_value_, `&options`)");

    min_width = Max(3, options[utility.getGlobalValue("terms.table_options.minimum_column_width")]);

    underline_chars = {
        {
            "-",
            "-",
            "-"
        }
    };
    dim = utility.Array1D(cells);

    row = "";
    row * 128;
    if (options[utility.getGlobalValue("terms.table_options.header")]) {
        underlines = "";
        underlines * 128;
        widths = {};

        if (utility.Has (options, utility.getGlobalValue("terms.table_options.column_widths"), "AssociativeList")) {
            widths = options[utility.getGlobalValue("terms.table_options.column_widths")]
        }


        if (options[utility.getGlobalValue("terms.table_options.align")] == "center") {
            underline_chars[0] = ':';
            underline_chars[2] = ':';
        } else {
            if (options[utility.getGlobalValue("terms.table_options.align")] == "right") {
                underline_chars[2] = ':';
            }

        }
        for (i = 0; i < dim; i += 1) {
            content_width = Abs(cells[i]);
            cell_width = Max (widths[i], Max(min_width, content_width));
            widths [i] =  cell_width;

            row * "|";
            padding = cell_width - content_width;

            for (k = 0; k < padding$2; k += 1) {
                row * " ";
            }
            row * (prefix + cells[i] + suffix);
            for (k = 0; k < padding - padding$2; k += 1) {
                row * " ";
            }
            underlines * "|";
            underlines * underline_chars[0];
            for (k = 1; k < cell_width - 1; k += 1) {
                underlines * underline_chars[1];
            }
            underlines * underline_chars[2];
        }
        row * "|";
        underlines * "|";
        underlines * 0;
        row * "\n";
        row * underlines;
        row * "\n";
        options[utility.getGlobalValue("terms.table_options.column_widths")] = widths;
    } else {
        for (i = 0; i < dim; i += 1) {
            content_width = Abs(cells[i]);
            cell_width = (options[utility.getGlobalValue("terms.table_options.column_widths")])[i];

            row * "|";
            if (cell_width <= content_width + 3) {
                cells[i] = (cells[i])[0][cell_width - 4] + "...";
                padding = 0;
            } else {
                padding = cell_width - content_width;
            }
            for (k = 0; k < padding$2; k += 1) {
                row * " ";
            }
            row * (prefix + cells[i] + suffix);
            for (k = 0; k < padding - padding$2; k += 1) {
                row * " ";
            }
        }
        row * "|\n";
    }
    row * 0;
    return row;
}

/**
 * @name io.FormatTableRow
 * @param row
 * @param options
 */
lfunction io.FormatTableRow(row, options) {
    return io.FormatTableRowDecorators (row, options, "", "");
}

/**
 * @name io.get_a_list_of_files
 * @param filename
 */
function io.get_a_list_of_files(filename) {
    if (Type(filename) == "String") {
        if (!filename) { // filename exists
            fscanf(filename, REWIND, "Lines", io.get_a_list_of_files.list);
            return io.validate_a_list_of_files(io.get_a_list_of_files.list);
        }
    }

    io.get_a_list_of_files.result = {};
    io.PrintAndUnderline("Enter paths to files (blank line to end entry)", "-");
    while (1) {
        fprintf(stdout, "* File ", Abs(io.get_a_list_of_files.result) + 1, " [relative path `PATH_TO_CURRENT_BF`]:");
        io.get_a_list_of_files.current_path = "";
        fscanf(stdin, "String", io.get_a_list_of_files.current_path);
        if (Abs(io.get_a_list_of_files.current_path)) {
            io.CheckAssertion("! io.get_a_list_of_files.current_path", "HyPhy cannot open '" + io.get_a_list_of_files.current_path + "' for reading");
        } else {
            break;
        }
        io.get_a_list_of_files.result + io.get_a_list_of_files.current_path;
    }
}

/**
 * @name io.DisplayAnalysisBanner
 * @param analysis_info
 */
lfunction io.DisplayAnalysisBanner(analysis_info) {


    if (io.HasStringKey(utility.getGlobalValue("terms.io.info"), analysis_info)) {
        io.PrintAndUnderline("Analysis Description", "-");
        fprintf(stdout, io.FormatLongStringToWidth(analysis_info[utility.getGlobalValue("terms.io.info")], 72), "\n");
    }
    if (io.HasStringKey(utility.getGlobalValue("terms.io.help"), analysis_info)) {
        fprintf(stdout, "\n- __Additional information__: ");
        fprintf(stdout, io.FormatLongStringToWidth(analysis_info[utility.getGlobalValue("terms.io.help")], 72), "\n");
    }
    if (io.HasStringKey(utility.getGlobalValue("terms.io.requirements"), analysis_info)) {
        fprintf(stdout, "\n- __Requirements__: ");
        fprintf(stdout, io.FormatLongStringToWidth(analysis_info[utility.getGlobalValue("terms.io.requirements")], 72), "\n");
    }
    if (io.HasStringKey(utility.getGlobalValue("terms.io.reference"), analysis_info)) {
        fprintf(stdout, "\n- __Citation__: ");
        fprintf(stdout, io.FormatLongStringToWidth(analysis_info[utility.getGlobalValue("terms.io.reference")], 72), "\n");
    }
    if (io.HasStringKey(utility.getGlobalValue("terms.io.authors"), analysis_info)) {
        fprintf(stdout, "\n- __Written by__: ");
        fprintf(stdout, io.FormatLongStringToWidth(analysis_info[utility.getGlobalValue("terms.io.authors")], 72), "\n");
    }
    if (io.HasStringKey(utility.getGlobalValue("terms.io.contact"), analysis_info)) {
        fprintf(stdout, "\n- __Contact Information__: ");
        fprintf(stdout, io.FormatLongStringToWidth(analysis_info[utility.getGlobalValue("terms.io.contact")], 72), "\n");
    }
    if (io.HasStringKey(utility.getGlobalValue("terms.io.version"), analysis_info)) {
        fprintf(stdout, "\n- __Analysis Version__: ");
        fprintf(stdout, io.FormatLongStringToWidth(analysis_info[utility.getGlobalValue("terms.io.version")], 72), "\n");
    }
    fprintf(stdout, "\n");

    return None;
}

/**
 * @name io.HasStringKey
 * @param key
 * @param dict
 */
lfunction io.HasStringKey(key, dict) {
    return Type(dict[key]) == "String";
}

/**
 * @name io.SpoolLF
 * @param lf_id
 * @param trunk_path
 * @param tag
 * @returns nothing
 */
lfunction io.SpoolLF(lf_id, trunk_path, tag) {

    Export(__lf_spool, ^ lf_id);
    if (tag == None) {
        tag = lf_id;
    }
    fprintf(trunk_path + "." + tag + ".bf", CLEAR_FILE, __lf_spool);
}

/**
 * @name io.SpoolLFToPath
 * @param lf_id
 * @param trunk_path
 * @param tag
 * @returns nothing
 */
lfunction io.SpoolLFToPath(lf_id, path) {
    if (path != "/dev/null") {
        Export(__lf_spool, ^ lf_id);
        io.CompressedObjectSpool (__lf_spool, path);
        DeleteObject (__lf_spool);
    }
}

/**
 * @name io.PrintAndUnderline
 * @param string
 * @param char
 * @returns nothing
 */
lfunction io.PrintAndUnderline(string, char) {
    fprintf(stdout, "\n", string, "\n");
    buffer = "";
    buffer * (1 + Abs(string));
    for (k = 0; k < Abs(string); k += 1) {
        buffer * char[0];
    }
    buffer * 0;
    fprintf(stdout, buffer, "\n");
}

/**
 * @name io.FormatLongStringToWidth
 * @param {String} string
 * @param {Number} width
 * @returns {String} formatted string
 */
function io.FormatLongStringToWidth(string, width) {
    words = regexp.Split(string, "[\\ \n]+");
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

/**
 * I am tired of typing fprintf (stdout, ...)
 * @returns nothing
 */
lfunction console.log (arg) {
    fprintf (stdout, arg, "\n");
}

/**
 * I am tired of having dangling console.log debug statements
 * @returns nothing
 */
lfunction debug.log (arg) {
    if (utility.GetEnvVariable ("_DEBUG_MESSAGES_ON_")) {
        fprintf (stdout, arg, "\n");
    }
}

/**
 * I am tired of typing fprintf (MESSAGE_LOG, ...)
 * @returns nothing
 */
lfunction messages.log (arg) {
    fprintf (MESSAGE_LOG, "\n", arg);
}


/**
 * I am tired of typing fprintf (stdout,  ...)
 * @returns nothing
 */
lfunction warning.log (arg) {
    fprintf (stdout, "[**WARNING**] ", arg, "\n");
}

/**
 * Checks if there is a file exists
 * @param {String} path the path to check
 * @returns {Number} TRUE if file exists and is readable; FALSE otherwise
 */

lfunction io.FileExists  (path) {
    return !path;
}

/**
 * Converts a relative path to the absolute path using the current batch file as base
 * @param {String} relative path
 * @returns {String} absolute path
 */

lfunction io.MakeAbsoluteFilePath  (path) {
    if (path [0] != "/") {
        return utility.getGlobalValue ("PATH_TO_CURRENT_BF") + path;
    }
    return path;
}


/**
 * Checks if there is a cache file; creates if empty
 * @param {String} path  the path to the cache file; will be created if it doesn't exist
 * @returns {Dict} the contents of the file
 */
lfunction io.LoadCacheFromFile  (path) {
    if (io.FileExists (path) == TRUE) { // exists
        fscanf (path, REWIND, "Raw", contents);
        contents =  Eval (contents);
        if (Type (contents) == "AssociativeList") {
            return contents;
        }
    } else {
        fprintf (path, CLEAR_FILE, {});
    }
    return {};
}

/**
 * Checks if there is a cache file; creates if empty
 * @param {String} path  the path to the cache file; will be created if it doesn't exist
 * @param {Dict} data  the contents of the cache to save to a file
 * @returns nothing
 */
lfunction io.WriteCacheToFile  (path, data) {
    fprintf (path, CLEAR_FILE, data);
}

/**
 * Format the time in seconds as (DD:)HH:MM:SS ('DD' is only shown if more than HH > 24 hours
 * @param {Number} interval  the time interval in seconds
 * @returns {String} formatted string
 */
lfunction io.FormatTimeInterval  (interval) {
    pieces  = {{interval % 60, interval % 3600 $ 60, interval $ 3600 % 24, interval $ (3600*24)}};

    time_pieces = {};

    for (k = 3; k >= 0; k = k - 1) {
        if (pieces[k] || k < 3) {
            if (pieces[k] < 10) {
                time_pieces + ("0" + pieces[k]);
            } else {
                time_pieces + pieces[k];
            }
        }
    }

    return Join (":", time_pieces);
}

/**
 * Read a delimiter separated file
 * @param {String}  path  path to file; if None, user will be prompted for file
 * @param {String}  separator the separator to use to segment a line into fields
 * @param {Boolean} has_header treat the first line as header
 * @returns {Dict}

  {
      "rows": {
          "0": {
              "0": "KC618402",
              "1": "Human"
          },
          "1": {
              "0": "KC618403",
              "1": "Human"
          },
          "2": {
              "0": "AB248520",
              "1": "Human"
          },
          ...
      }
      "header": {
          "0": "accession",
          "1": "Host"
      }
  }

 */
lfunction io.ReadDelimitedFile  (path, separator, has_header) {
   if (Type (path) == "String") {
        fscanf (path, REWIND, "Lines", data);
   } else {
        fscanf (PROMPT_FOR_FILE, REWIND, "Lines", data);
        path = utility.getGlobalValue("LAST_FILE_PATH");
   }
   result = {utility.getGlobalValue("terms.io.rows") : {}};
   index = 0;
   row_count = utility.Array1D (data);
   if (has_header) {
        result[utility.getGlobalValue("terms.io.header")] = regexp.Split (data[0], separator);
        index = 1;
   }
   for (k = index; k < row_count; k+=1) {
        result [utility.getGlobalValue("terms.io.rows")] + regexp.Split (data[k], separator);
   }
   result[utility.getGlobalValue("terms.json.file")] = path;
   return result;
}

/**
 * Present a selection dialog and return an option
 * @param {Dict/Matrix} options key : value or {{key, value}...}, using the matrix version ensures ordering
 * @param {String} description dialog caption
 * @returns {String/None} selected key or None if selection canceled
 */

lfunction io.SelectAnOption  (options, description) {
    option_set = {utility.Array1D (options),2};
    selection = None;
    if (Rows (option_set) > 0) {
        if (Type (options) == "Matrix") {
            option_set = options;
        } else {
            keys = utility.Keys (options);
            for (k = 0; k < Rows (option_set); k+=1) {
                option_set [k][0] = keys[k];
                option_set [k][1] = options[keys[k]];
            }
        }

        ChoiceList  (selection,description,1,NO_SKIP,option_set);
        
        if (selection >= 0) {
            return option_set[selection][0];
        } else {
            selection = None;
        }
   }
   assert (None != selection, "Selection canceled");
   return None;
}

/**
 * Present a selection dialog and return one or more options (specified via "count")
 * @param {Dict/Matrix} options key : value or {{key, value}...}, using the matrix version ensures ordering
 * @param {String} description dialog caption
 * @returns {Matrix/None} selected key or None if selection canceled
 */

lfunction io.MultiSelectOptions  (options, description, count) {
    option_set = {utility.Array1D (options),2};
    selection = None;
    if (Rows (option_set) > 0) {
        if (Type (options) == "Matrix") {
            option_set = options;
        } else {
            keys = utility.Keys (options);
            for (k = 0; k < Rows (option_set); k+=1) {
                option_set [k][0] = keys[k];
                option_set [k][1] = options[keys[k]];
            }
        }

        ChoiceList  (selection,description,count,NO_SKIP,option_set);
        
        N = utility.Array1D (selection);
        if (N > 0) {
            result = {N,1};
            for (i,k; in; selection) {
                result [i] = option_set[k][0];
            }
            return result;
        } else {
            selection = None;
        }
   }
   assert (None != selection, "Selection canceled");
   return None;
}

/**
 * @name io.ReportExecutionError
 * @param error_msg
 */
function io.ReportAnExecutionError (error_msg) {
    assert (0, "Fatal execution error `error_msg`");
}

/**
 * @name io.ReportWarning
 * @param warn_message
 */
function io.ReportWarning (warn_message) {
    fprintf(stdout, "\n-------\n", 
     io.FormatLongStringToWidth(">[WARNING] " +  warn_message, 72), 
    "\n-------\n");
}

/**
 * @name io.SingularOrPlural
 * @param value
 * @param singular
 * @param plural
 */
lfunction io.SingularOrPlural (value, singular, plural) {
   if (value == 1) {
        return singular;
   }
   return plural;
}

/**
 * @name io.splitFilePath
 * Split file path into a directory, file name, and extension (if any)
 * @param {String} path
 * @return {Matrix} {{dir, filename, ext}}
 */

lfunction io.splitFilePath (path) {
	result = {{"","",""}};
	split     = path $ ("[^\\"+^"DIRECTORY_SEPARATOR"+"]+$");
	if (split[0] != 0 || split[1] != Abs (path)-1) {
		result [0] = path[0][split[0]-1];
		path = path[split[0]][Abs(path)];
	}

	split     = path || "\\.";
	if (split[0] < 0) {
		result[1]  = path;
 	} else {
		result [2] =  path[split[Rows(split)-1]+1][Abs(path)-1];
		result [1]  = path[0][split[Rows(split)-1]-1];
	}
	return result;
}
