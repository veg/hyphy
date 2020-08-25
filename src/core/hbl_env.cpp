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
    HBLObjectRef       EnvVariableGetDefault (_String const& name, unsigned long type) {
        HBLObjectRef default_value = (HBLObjectRef)_hy_env_default_values.GetDataByKey (&name);
        if (default_value) {
            if (type == HY_ANY_OBJECT || (type | default_value->ObjectClass())) {
                return default_value;
            }
        }
        return nil;
    }

    /*********************************************************************************/
    hyFloat       EnvVariableGetDefaultNumber (_String const& name) {
        HBLObjectRef default_value = (HBLObjectRef)EnvVariableGetDefault (name, NUMBER);
        if (default_value) {
            return default_value->Value();
        }
        return HY_INVALID_RETURN_VALUE;
    }

    /*********************************************************************************/
    bool       EnvVariableTrue (_String const& name) {
        HBLObjectRef value = (HBLObjectRef)EnvVariableGet (name, NUMBER);
        if (value) {
            return fabs (value->Value()) > 1e-10;
        }
        return false;
    }

    /*********************************************************************************/
    HBLObjectRef       EnvVariableGet (_String const& name, unsigned long type) {
        HBLObjectRef current_value = (HBLObjectRef)FetchObjectFromVariableByType(&name, type);
        if (current_value) {
            if (type == HY_BL_ANY || (type | current_value->ObjectClass())) {
                return current_value;
            }
        }
        return EnvVariableGetDefault (name, type);
    }

    /*********************************************************************************/
    hyFloat       EnvVariableGetNumber (_String const& name, hyFloat default_value) {
      HBLObjectRef current_value = EnvVariableGet (name, NUMBER);
      return current_value ? current_value -> Value() : default_value;
    }

  /*********************************************************************************/
    void       EnvVariableSet (_String const& name, HBLObjectRef value, bool copy) {
        EnvVariableSetNamespace (name, value, nil, copy);
    }

    /*********************************************************************************/
    void       EnvVariableSetNamespace (_String const& name, HBLObjectRef value, _String* nmspace, bool copy) {
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
                              .PushPairCopyKey (dataset_save_memory_size, new _Constant (100000.))
                              .PushPairCopyKey (harvest_frequencies_gap_options, new HY_CONSTANT_TRUE)
                              .PushPairCopyKey(assertion_behavior, new HY_CONSTANT_FALSE)
                              .PushPairCopyKey (print_float_digits, new _Constant (0.))
                              .PushPairCopyKey (mpi_node_count, new _Constant (1.))
                              .PushPairCopyKey (always_reload_libraries, new HY_CONSTANT_FALSE)
                              .PushPairCopyKey (end_of_file, new HY_CONSTANT_FALSE)
                              .PushPairCopyKey (produce_markdown_output, new HY_CONSTANT_FALSE)
                              .PushPairCopyKey (integration_maximum_iterations, new _Constant (10.))
                              .PushPairCopyKey (integration_precision_factor, new _Constant (1.e-10))
                              .PushPairCopyKey (skip_omissions, new _Constant (HY_CONSTANT_FALSE))
                              .PushPairCopyKey (data_file_print_format, new _Constant (6.0))
                              .PushPairCopyKey (data_file_default_width, new _Constant (50.0))
                              .PushPairCopyKey (data_file_gap_width, new _Constant (10.0))
                              .PushPairCopyKey (accept_branch_lengths, new _Constant (HY_CONSTANT_TRUE))
      ;
    }

_String cli_env_settings;
  
_String const
    accept_branch_lengths                            ("ACCEPT_BRANCH_LENGTHS"),
        // if true (default), then branch lengths from Newick strings will be accepted (whenever possible)
    accept_rooted_trees                             ("ACCEPT_ROOTED_TREES"),
        // if TRUE, do not perform automatic unrooting for topology/tree constructors
    always_reload_libraries                         ("ALWAYS_RELOAD_FUNCTION_LIBRARIES"),
        // if TRUE, reparse and re-execute source code for each call to LoadFunctionLibrary,
        // otherwise load function libraries only once
    assertion_behavior                              ("ASSERTION_BEHAVIOR"),
        // if set to TRUE, then assertions that fail skip to the end of the current script
        // otherwise they terminate the program
    assume_reversible                               ("ASSUME_REVERSIBLE_MODELS"),
        // 0 : check reversibility at run-time
        // 1 : ASSUME reversibility
        // -1 : ASSUME NON-reversibility
  
    automatically_convert_branch_lengths            ("AUTOMATICALLY_CONVERT_BRANCH_LENGTHS"),
        // if TRUE, then HyPhy will attempt to solve BL (t) = C for model parameter t, whenever possible
    base_directory                                  ("HYPHY_BASE_DIRECTORY"),
        // is set to the base directory for local path names; can be set via a CL argument (BASEPATH)
    blockwise_matrix                                ("BLOCK_LIKELIHOOD"),
        // this _template_ variable is used to define likelihood function evaluator templates
    branch_length_stencil                           ("BRANCH_LENGTH_STENCIL"),
    covariance_parameter                            ("COVARIANCE_PARAMETER"),
        // used to control the behavior of CovarianceMatrix
    data_file_default_width                         ("DATA_FILE_DEFAULT_WIDTH"),
      // for file formats with grouped alignment columns (e.g. PHYLIP), determines the width of a column
    data_file_gap_width                             ("DATA_FILE_GAP_WIDTH"),
      // for file formats with grouped alignment columns (e.g. PHYLIP), determines the width of the gap between
      // column groups
    data_file_partition_matrix                      ("DATA_FILE_PARTITION_MATRIX"),
        // the string of data partitions read from the last valid NEXUS CHARSET block
    data_file_print_format                          ("DATA_FILE_PRINT_FORMAT"),
      // determines the file format for datasets and datafilters

    data_file_tree                                  ("IS_TREE_PRESENT_IN_DATA"),
        // set to TRUE if the last data load call yielded a Newick trees
    data_file_tree_string                           ("DATAFILE_TREE"),
        // the last tree loaded by via a sequence file read
    dataset_save_memory_size                        ("USE_MEMORY_SAVING_DATA_STRUCTURES"),
        // sets the maximum dimension of a data filter for generating .site_map, .site_freqs, .sequence_map
    directory_separator_char                        ("DIRECTORY_SEPARATOR"),
        // is set to the platform directory separator (e.g. '/')
    defer_constrain_assignment                        ("DEFER_CONSTRAINT_APPLICATION"),
        // if set to TRUE, then constraint application will be done in a single batch
        // this is helpful when many x := expr statements are strung together to avoid
        // checking the entire namespace for dependancies
    end_of_file                                     ("END_OF_FILE"),
        // set by IO operations, like fscanf to indicate whetehr the end of the input stream has
        // been reached
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
    execution_mode                                  ("HBL_EXECUTION_ERROR_HANDLING"),
    // sets HyPhy exception handling
    // FALSE - bail out on errors
    // TRUE  - return from the current execution list, but keep the program running (e.g. to allow for HBL testing of error handling)

    get_data_info_returns_only_the_index            ("GET_DATA_INFO_RETURNS_ONLY_THE_INDEX"),
    // instead of returing {0,0,1,0} for a 'G' character in GetDataInfo (r, filter, species, pattern)
    // return only the index of 'G', e.g. 2 in this case. -1 is returned for ambigs
    false_const                                     ("FALSE"),
        // the FALSE (0.0) constant
    fprintf_redirect                                ("GLOBAL_FPRINTF_REDIRECT"),
        // if set to a string path, then all stdout will go to the file at the path  (or /dev/null) instead of to stdout

    harvest_frequencies_gap_options                 ("COUNT_GAPS_IN_FREQUENCIES"),
    /* 
        if set to `harvest_frequencies_gap_options` to TRUE, then N-fold ambigs will add 1/N to each character count in HarvestFrequencies,
        otherwise N-folds are ignored in counting
     */
    
    include_model_spec                              ("INCLUDE_MODEL_SPECS"),
    /*
        controls whether or not export / serialization operations (like toStr)
        will include substitution model specifications
     */

    integration_precision_factor                    ("INTEGRATION_PRECISION_FACTOR"),
    integration_maximum_iterations                  ("INTEGRATION_MAX_ITERATES"),
    // used to control integration in _Formula::Integrate
  
    kExpectedNumberOfSubstitutions                  ("EXPECTED_NUMBER_OF_SUBSTITUTIONS"),
        // literal for the expected number of substitions (per unit time)
    kGetStringFromUser                              ("PROMPT_FOR_STRING"),
        // [LEGACY] placeholder for prompting the user for a string value
    kSCFGCorpus                                     ("SCFG_STRING_CORPUS"),
        // set SCFG training corpus
    kStringSuppliedLengths                          ("STRING_SUPPLIED_LENGTHS"),
        // literal for branch lengths from the Newick tree string
    kDevNull                                        ("/dev/null"),
        // literal for branch lengths from the Newick tree string
    last_file_path                                  ("LAST_FILE_PATH"),
        // is set by various file read/write commands (fscanf, fprintf, dialog prompts)
        // to contain the **absolute** path to the last file interacted with
    last_fileio_exception                           ("LAST_FILE_IO_EXCEPTION"),
        // set to the value of the last exception if soft_fileio_exceptions is true,

    last_raw_file_prompt                            ("LAST_RAW_FILE_PROMPT"),
           // the last unprocessed value obtained by PROMPT_FOR_FILE
    last_model_parameter_list                       ("LAST_MODEL_PARAMETER_LIST"),
        // a stand-in for the list of model parameters for the last
        // declared model
    
    lf_convergence_criterion                        ("LF_CONVERGENCE_CRITERION"),
        // if set to a string, provides a callback function ID to LF optimization routines,
        // expected to take two arguments: current log L and a dict with current param values
        // returns T/F for convergence criterion check
    
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
    print_float_digits                              ("PRINT_DIGITS"),
        // controls how many decimal places are generated by fprintf and various number->string conversions
    produce_markdown_output                         ("MARKDOWN_OUTPUT"),
        // controls if certain stdout output is formatted as MarkDown (default is not)
    random_seed                                     ("RANDOM_SEED"),
        // the seed used for the Mersenne Twister random number generator
    selection_strings                               ("SELECTION_STRINGS"),
        // populated by a successful call to 'ChoiceList', @see ExecuteCase32
    sitewise_matrix                                 ("SITE_LIKELIHOOD"),
        // this _template_ variable is used to define likelihood function evaluator templates
    short_mpi_return                                ("SHORT_MPI_RETURN"),
        // controls the return format of optimized functions from MPI slave nodes
    skip_omissions                                  ("SKIP_OMISSIONS"),
        // if set, will cause data filters to _EXCLUDE_ sites with gaps or other N-fold redundancies
    soft_fileio_exceptions                          ("SOFT_FILE_IO_EXCEPTIONS"),
        // if set, read/write errors from fscanf and fprintf will not cause a program termination
        // but rather set the last_fileio_exception variable to the value of the exception
    status_bar_update_string                        ("STATUS_BAR_STATUS_STRING"),
        // used to set the progress message displayed to the user
    tolerate_numerical_errors                       ("TOLERATE_NUMERICAL_ERRORS"),
        // if set, numerical errors that would cause termination are instead trated as warnings
    try_numeric_sequence_match                      ("TRY_NUMERIC_SEQUENCE_MATCH"),
        // try matching sequences by 0 (or 1) based index, if matching by name fails
    true_const                                      ("TRUE"),
        // the TRUE (1.0) constant
    use_last_model                                  ("USE_LAST_MODEL"),
        // a stand-in for the last declared model
    use_traversal_heuristic                         ("USE_TRAVERSAL_HEURISTIC"),
        // TODO (20170413): don't remember what this does; , see @ _DataSetFilter::MatchStartNEnd
        // #DEPRECATE
    verbosity_level_string                           ("VERBOSITY_LEVEL")
        // controls verbosity level during optimization and other long-running operations
    ;
    
 
/** default values get set up in hy_global:: */
    
    
}

