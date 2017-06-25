/*
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon42@uwo.ca)
 Steven Weaver (sweaver@temple.edu)
 
 Module Developers:
 Lance Hepler (nlhepler@gmail.com)
 Martin Smith (martin.audacis@gmail.com)
 
 Significant contributions from:
 Spencer V Muse (muse@stat.ncsu.edu)
 Simon DW Frost (sdf22@cam.ac.uk)
 
 Permission is hereby granted, free of charge, to any person obtaining a
 copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject to
 the following conditions:
 
 The above copyright notice and this permission notice shall be included
 in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

#include "hbl_env.h"

/**
    HBL environment variables and constants go here,
    Getter/setter functions are also here
*/


namespace hy_env {
    
    _List           _hy_env_default_values_aux;
    
    _AVLListXL      _hy_env_default_values (&_hy_env_default_values_aux);
 
    /*********************************************************************************/
    _PMathObj       EnvVariableGetDefault (_String const& name, unsigned long type) {
        _PMathObj default_value = (_PMathObj)_hy_env_default_values.GetDataByKey (&name);
        if (default_value) {
            if (type == HY_ANY_OBJECT || (type | default_value->ObjectClass())) {
                return default_value;
            }
        }
        return nil;
    }

    /*********************************************************************************/
    hyFloat       EnvVariableGetDefaultNumber (_String const& name) {
        _PMathObj default_value = (_PMathObj)EnvVariableGetDefault (name, NUMBER);
        if (default_value) {
            return default_value->Value();
        }
        return HY_INVALID_RETURN_VALUE;
    }

    /*********************************************************************************/
    bool       EnvVariableTrue (_String const& name) {
        _PMathObj default_value = (_PMathObj)EnvVariableGetDefault (name, NUMBER);
        if (default_value) {
            return default_value->Value() != 0.0;
        }
        return false;
    }

    /*********************************************************************************/
    _PMathObj       EnvVariableGet (_String const& name, unsigned long type) {
        _PMathObj current_value = (_PMathObj)FetchObjectFromVariableByType(&name, type);
        if (current_value) {
            if (type == HY_BL_ANY || (type | current_value->ObjectClass())) {
                return current_value;
            }
        }
        return EnvVariableGetDefault (name, type);
    }

    /*********************************************************************************/
    void       EnvVariableSet (_String const& name, _PMathObj value, bool copy) {
        EnvVariableSetNamespace (name, value, nil, copy);
    }

    /*********************************************************************************/
    void       EnvVariableSetNamespace (_String const& name, _PMathObj value, _String* nmspace, bool copy) {
        setParameter(name, value, nmspace, copy);
    }
    /**
        arrange variables alphabetically by C++ ID for easier scanning;
        defaut  values are stored can be obtained by calling
            EnvVariableGetDefault
     
     */
    

    /*********************************************************************************/
    
    void        SetupEnvDefaults (void) {
        _hy_env_default_values.PushPairCopyKey (use_traversal_heuristic,  new HY_CONSTANT_TRUE)
                              .PushPairCopyKey (normalize_sequence_names, new HY_CONSTANT_TRUE)
                              .PushPairCopyKey(message_logging, new HY_CONSTANT_TRUE)
                              .PushPairCopyKey (dataset_save_memory_size, new _Constant (100000.));

    }
    
_String const
    base_directory                                  ("HYPHY_BASE_DIRECTORY"),
        // is set to the base directory for local path names; can be set via a CL argument (BASEPATH)
    blockwise_matrix                                ("BLOCK_LIKELIHOOD"),
        // this _template_ variable is used to define likelihood function evaluator templates
    data_file_partition_matrix                      ("DATA_FILE_PARTITION_MATRIX"),
        // the string of data partitions read from the last valid NEXUS CHARSET block
    data_file_tree                                  ("IS_TREE_PRESENT_IN_DATA"),
        // set to TRUE if the last data load call yielded a Newick trees
    data_file_tree_string                           ("DATAFILE_TREE"),
        // the last tree loaded by via a sequence file read
    dataset_save_memory_size                        ("USE_MEMORY_SAVING_DATA_STRUCTURES"),
        // sets the maximum dimension of a data filter for generating .site_map, .site_freqs, .sequence_map
    directory_separator_char                        ("DIRECTORY_SEPARATOR"),
        // is set to the platform directory separator (e.g. '/')
    error_report_format_expression                  ("ERROR_REPORT_FORMAT_EXPRESSION"),
        // if provided, this expression (assumed string valued), will be used to format the error
        // message, with special placeholder variables (see below) will be replaced with the
        // error related data
    error_report_format_expression_stack            ("_ERROR_CALL_STACK_"),
        // the current HBL call stack, formatted as a list
    error_report_format_expression_stdin           ("_ERROR_CALL_STDIN_"),
        // the current HBL standard input buffer, formatted as a list
    error_report_format_expression_string           ("_ERROR_TEXT_"),
    // the text message explaining the error
    
    false_const                                     ("FALSE"),
        // the FALSE (0.0) constant
    kGetStringFromUser                              ("PROMPT_FOR_STRING"),
        // [LEGACY] placeholder for prompting the user for a string value
    last_file_path                                  ("LAST_FILE_PATH"),
        // is set by various file read/write commands (fscanf, fprintf, dialog prompts)
        // to contain the **absolute** path to the last file interacted with

    lib_directory                                   ("HYPHY_LIB_DIRECTORY"),
        // is set to the library directory for standard library searchers; can be set via a CL argument (LIBPATH)
    matrix_element_column                           ("_MATRIX_ELEMENT_COLUMN_"),
    matrix_element_row                              ("_MATRIX_ELEMENT_ROW_"),
    matrix_element_value                            ("_MATRIX_ELEMENT_VALUE_"),
        // the last three variables are used as _template_ variable for conditional / iterated matrix operations, e.g.,
        // matrix ["_MATRIX_ELEMENT_ROW_+_MATRIX_ELEMENT_COLUMN_"]
    message_logging                                 ("MESSAGE_LOGGING"),
        // if set, then diagnostic messages will be logged
    mpi_node_id                                     ("MPI_NODE_ID"),
        // [MPI only] the ID (0 = master, etc) for teh current process
    mpi_node_count                                  ("MPI_NODE_COUNT"),
        // [MPI only] the count of MPI nodes (master + slaves)
    mpi_last_sent_message                           ("MPI_LAST_SENT_MSG"),
        // [MPI only] the contents of the last message sent by the current node
    nexus_file_tree_matrix                          ("NEXUS_FILE_TREE_MATRIX"),
        // the tree matrix read from the last valid NEXUS TREE block
    normalize_sequence_names                        ("NORMALIZE_SEQUENCE_NAMES"),
        // if set, will trigger automatic renaming of sequence names from files to valid
        // HyPhy IDs, e.g. "awesome monkey!" -> "awesome_monkey_"
        // the mapping will go into dataset_id.mapping
    path_to_current_bf                              ("PATH_TO_CURRENT_BF"),
        // is set to the absolute path for the currently executed batch file (assuming it has one)

    random_seed                                     ("RANDOM_SEED"),
        // the seed used for the Mersenne Twister random number generator
    selection_strings                               ("SELECTION_STRINGS"),
        // populated by a successful call to 'ChoiceList', @see ExecuteCase32
    sitewise_matrix                                 ("SITE_LIKELIHOOD"),
        // this _template_ variable is used to define likelihood function evaluator templates
    status_bar_update_string                        ("STATUS_BAR_STATUS_STRING"),
        // used to set the progress message displayed to the user
    true_const                                      ("TRUE"),
        // the TRUE (1.0) constant
    use_traversal_heuristic                         ("USE_TRAVERSAL_HEURISTIC")
        // TODO (20170413): don't remember what this does; , see @ _DataSetFilter::MatchStartNEnd
        // #DEPRECATE
    ;
    
 
/** default values get set up in hy_global:: */
    
    
}

