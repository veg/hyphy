LoadFunctionLibrary("convenience/regexp.bf");

function io.checkAssertion(statement, error_msg) {
    ExecuteCommands("assert (`statement`, error_msg)");
}

lfunction io.prompt_user (prompt,default,lower_bound,upper_bound,is_integer) {
	value = lower_bound-1;

	while (value < lower_bound || value > upper_bound) {
		fprintf (stdout, prompt, " (permissible range = [", lower_bound, ",", upper_bound, "], default value = ", default);
		if (is_integer) {
			fprintf (stdout, ", integer");
		}
		fprintf (stdout, "): ");
		fscanf  (stdin, "String", str_val);

		if (Abs(str_val) == 0) {
			value = 0+default;
		} else {
		    value = 0+str_val;
		}
		if (is_integer) {
		    value = value $ 1;
		}
	}
	return value;
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

lfunction io.reportAnalysisStageMD(stage) {
    fprintf (stdout, "\n>", stage, "\n\n");
}

lfunction io.reportProgressMessageMD(analysis, stage, text) {
    if (Abs(cache) == 0) {
        cache = {};
    }
    advance = TRUE;
    utility.dict.ensure_key (cache, analysis);

    if ((cache[analysis])[stage]) {
        advance = FALSE;
    }
    (cache[analysis])[stage] += 1;

    if (advance) {
        if (Abs (cache[analysis]) == 1) {
            fprintf (stdout, "\n");
        }
        fprintf(stdout, "\n### ", text, "\n");
    } else {
        fprintf(stdout, text, "\n");
    }
}

lfunction io.reportProgressBar(analysis, text) {
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

lfunction io.format_object (object, options) {
    if (Type (object) == "String") {
        return object;
    }
    if (Type (object) == "Number") {
        if (None != options) {
            if (Abs (options["number-precision"]) > 0) {
                return Eval ("Format (`&object`, 0, " + options["number-precision"] + ")");
            }
        }
    }

    return ""+object;
}


lfunction io.format_table_row (row, options) {

    if (None == options) {
        options = {};
    }

    cells = utility.map (row, "_value_", "io.format_object(_value_, `&options`)");

    min_width = Max (3, options ["min-column-width"]);

    underline_chars = {{"-","-","-"}};
    dim = utility.array1D (cells);

    row = ""; row * 128;
    if (options ["header"]) {
        underlines = ""; underlines * 128;
        widths = {};

        if (options["align"] == "center") {
            underline_chars[0] = ':'; underline_chars[2] = ':';
        } else {
            if (options["align"] == "right") {
                underline_chars[2] = ':';
            }

        }
        for (i = 0; i < dim; i += 1) {
            content_width = Abs (cells[i]);
            cell_width    = Max (min_width, content_width);
            widths + cell_width;
            row * "|";
            padding = cell_width - content_width;

            for (k = 0; k < padding$2; k+=1) {
                row * " ";
            }
            row * cells[i];
            for (k = 0; k < padding - padding$2; k+=1) {
                row * " ";
            }
            underlines * "|";
            underlines * underline_chars[0];
            for (k = 1; k < cell_width - 1; k+=1) {
             underlines * underline_chars[1];
            }
            underlines * underline_chars[2];
        }
        row * "|"; underlines * "|"; underlines * 0;
        row * "\n";
        row * underlines;
        row * "\n";
        options ["column-widths"] = widths;
    } else {
        for (i = 0; i < dim; i += 1) {
            content_width = Abs (cells[i]);
            cell_width    = (options ["column-widths"])[i];

            row * "|";
            if (cell_width <= content_width + 3) {
                cells[i] = (cells[i])[0][cell_width-4] + "...";
                padding = 0;
            } else {
                padding = cell_width - content_width;
            }
            for (k = 0; k < padding$2; k+=1) {
                row * " ";
            }
            row * cells[i];
            for (k = 0; k < padding - padding$2; k+=1) {
                row * " ";
            }
        }
        row * "|\n";
    }
    row * 0;
    return row;
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

lfunction io.displayAnalysisBanner(analysis_info) {
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

lfunction io.hasStringKey(key, dict) {
    return Type(dict[key]) == "String";
}

lfunction io.spoolLF(lf_id, trunk_path, tag) {

    Export (__lf_spool, ^lf_id);
    if (tag == None) {
        tag = lf_id;
    }
    fprintf(trunk_path + "." + tag + ".bf", CLEAR_FILE, __lf_spool);
}

lfunction io.printAndUnderline(string, char) {
    fprintf(stdout, "\n", string, "\n");
    buffer = ""; buffer * (1+Abs (string));
    for (k = 0; k < Abs(string); k += 1) {
        buffer * char[0];
    }
    buffer * 0;
    fprintf(stdout, buffer, "\n");
}

function io.formatLongStringToWidth(string, width) {
    words = regexp.split (string, "[\\ \n]+");
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
