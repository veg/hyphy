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

#include "avllistxl.h"
#include "defines.h"
#include "hy_strings.h"
#include "parser.h"

#ifndef _HY_ENV
#define _HY_ENV

namespace hy_env {

/**
 * @brief Check if the value of the environment variable is "true".
 *
 * @param name The name of the environment variable (in HBL).
 * @param default_true The default value to return if the variable is not
 * defined.
 * @return true if the variable is defined and true-like (not zero), false
 * otherwise.
 */
bool EnvVariableTrue(_String const &name, bool default_true = false);

/**
 * @brief Look up the default numeric value for a given environment variable.
 *
 * @param name The name of the environment variable (in HBL).
 * @return The default value if available and of the correct type, otherwise
 * HY_INVALID_RETURN_VALUE.
 */
hyFloat EnvVariableGetDefaultNumber(_String const &name);

/**
 * @brief Look up the default value for a given environment variable, checking
 * that it is of a particular type.
 *
 * @param name The name of the environment variable (in HBL).
 * @param type The expected type of the variable.
 * @return The default value if available and of the correct type, otherwise
 * nil.
 */
HBLObjectRef EnvVariableGetDefault(_String const &name, unsigned long type);

/**
 * @brief Look up the value for a given environment variable, checking that it
 * is of a particular type.
 *
 * @param name The name of the environment variable (in HBL).
 * @param type The expected type of the variable.
 * @return The current value if set and of the correct type, otherwise the
 * default value if available and of the correct type, otherwise nil.
 */
HBLObjectRef EnvVariableGet(_String const &name, unsigned long type);

/**
 * @brief Look up the numeric value for a given environment variable.
 *
 * @param name The name of the environment variable (in HBL).
 * @param default_value The default value to return if the variable is not
 * defined.
 * @return The current value if set and of the correct type, otherwise the
 * default value if available and of the correct type, otherwise
 * HY_INVALID_RETURN_VALUE.
 */
hyFloat EnvVariableGetNumber(_String const &name,
                             hyFloat default_value = HY_INVALID_RETURN_VALUE);

/**
 * @brief Set the value for an environment variable (global scope).
 *
 * @param name The name of the environment variable (in HBL).
 * @param value The value to set the variable to.
 * @param copy If true, make a deep copy of 'value' prior to setting.
 */
void EnvVariableSet(_String const &name, HBLObjectRef value, bool copy);

/**
 * @brief Set the value for an environment variable (in a given namespace if
 * provided).
 *
 * @param name The name of the environment variable (in HBL).
 * @param value The value to set the variable to.
 * @param nmspace The namespace to set the variable in.
 * @param copy If true, make a deep copy of 'value' prior to setting.
 */
void EnvVariableSetNamespace(_String const &name, HBLObjectRef value,
                             _String *nmspace, bool copy);

/**
 * @brief Populate default value arrays.
 */
void SetupEnvDefaults(void);

extern const _String accept_rooted_trees, accept_branch_lengths,
    automatically_convert_branch_lengths, always_reload_libraries,
    assertion_behavior, assume_reversible, data_file_tree,
    data_file_tree_string, nexus_file_tree_matrix, data_file_partition_matrix,
    use_traversal_heuristic, normalize_sequence_names, dataset_save_memory_size,
    sitewise_matrix, blockwise_matrix, execution_mode, covariance_parameter,
    selection_strings, random_seed, assigned_seed, base_directory,
    lib_directory, directory_separator_char, path_to_current_bf,
    print_float_digits, true_const, false_const, fprintf_redirect,
    harvest_frequencies_gap_options, last_file_path, matrix_element_row,
    matrix_element_column, matrix_element_value, message_logging, mpi_node_id,
    mpi_node_count, mpi_last_sent_message, error_report_format_expression,
    error_report_format_expression_string, error_report_format_expression_stack,
    error_report_format_expression_stdin, status_bar_update_string,
    use_last_model, last_model_parameter_list, kGetStringFromUser,
    get_data_info_returns_only_the_index, defer_constrain_assignment,
    end_of_file, file_created, produce_markdown_output,
    integration_precision_factor, integration_maximum_iterations,
    skip_omissions, data_file_gap_width, data_file_default_width,
    data_file_print_format, branch_length_stencil,
    kExpectedNumberOfSubstitutions, kStringSuppliedLengths, kDevNull,
    soft_fileio_exceptions, last_fileio_exception, last_raw_file_prompt,
    include_model_spec, lf_convergence_criterion, try_numeric_sequence_match,
    short_mpi_return, kSCFGCorpus, verbosity_level_string,
    enforce_constraint_violations, tolerate_numerical_errors,
    tolerate_constraint_violation, number_threads, gzip_compression_level,
    tree_parser_namespace, global_kwargs, strict_alignment_validation_mode;

;

extern _String cli_env_settings;

}; // namespace hy_env

#endif
