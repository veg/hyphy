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

/**
 * @file global_things.h
 * @brief Global functions and variables that don't naturally belong with classes or other headers go here.
 */

#ifndef __HYGLOBALTHINGS__
#define __HYGLOBALTHINGS__

#include "avllistx.h"
#include "defines.h"
#include "hy_types.h"
#include "trie.h"
#include <regex>

#include <stdio.h>

#ifdef __HYPHYMPI__
#include <mpi.h>
#endif

class _Variable;      // forward decl
class _ExecutionList; // forward decl

namespace hy_global {

/**
 Global functions
 */

/**
 * @brief Heap-allocate a specified number of bytes and return the pointer.
 * Optionally zero the memory. If allocation fails, the program is halted via a fatal run-time error.
 *
 * @param bytes The number of bytes to allocate.
 * @param zero If true, zero the memory block.
 * @param alignment The alignment of the memory block.
 * @return A pointer to the new memory block.
 * @see FlagError, MemReallocate
 */
hyPointer MemAllocate(size_t bytes, bool zero = false, size_t alignment = 0);

/**
 * @brief Resize an existing pointer to 'new_bytes' bytes.
 * If allocation fails, the program is halted via a fatal run-time error.
 *
 * @param old_pointer A previously allocated (with MemAllocate) pointer.
 * @param new_bytes The new number of bytes to allocate.
 * @return A pointer to the resized memory block; could be different from old_pointer.
 * @see FlagError, MemReallocate
 */
hyPointer MemReallocate(hyPointer old_pointer, size_t new_bytes);

/**
 * @brief Set up the application environment.
 * This function initializes various global variables, paths, and constants.
 * - Initializes the list of global variables (_hy_application_globals).
 * - Sets up library paths (_hy_standard_library_paths) and extensions (_hy_standard_library_extensions).
 * - Initializes parser arrays (_HBL_Init_Const_Arrays).
 * - Sets up environment variables and language constants (e.g., TRUE, FALSE, HYPHY_BASE_DIRECTORY, HYPHY_LIB_DIRECTORY).
 * - Calls the global setup for hy_env.
 */
void InitializeGlobals(void);

/**
 * @brief Get the directory separator character for the current platform.
 *
 * @return The directory separator character (e.g., '/' for *nix, '\\' for Windows).
 */
char get_platform_directory_char(void);

/* pass-through structure for reading / writing from a file that may or may not
 * be compressed */

/**
 * @brief Open the file located at file_path using mode 'mode'.
 * This function is a wrapper around the standard `fopen` function that also handles compressed files.
 *
 * @param file_path The path of the file to open.
 * @param mode The standard (fopen) open modes.
 * @param error If true, then halt the program if the file could not be opened.
 * @param compress If true, the file will be opened as a compressed file.
 * @return The `hyFile` handle or nil (if file does not exist or could not be open with the requested mode).
 */
hyFile *doFileOpen(const char *file_path, hyFileOpenMode mode,
                   bool error = false, bool compress = false);

/**
 * @brief The omnibus clean-up function that attempts to deallocate all application memory and unwind various objects created.
 * The idea is that upon successful completion, the state of the program is the same as it was after the initial startup.
 *
 * @param all If true, a full purge is undertaken. If false, a partial purge is undertaken, that clears only user functions and attendant stuctures.
 * @todo 20170414: does the function ever get called with the FALSE flag? if not, deprecate I don't recall what the use case for the FALSE flag was
 */
void PurgeAll(bool all = true);

/**
 * @brief This a general initialization function that deals with bookkeeping.
 * - Calls all other initializers.
 * - Defines and sets the global random seed.
 * - Creates and opens messages and errors logs (based on settings).
 * @return true if the startup was successful, false otherwise.
 */
bool GlobalStartup(void);

/**
 * @brief The general clean up function that is called right before HyPhy exits.
 * - Calls other object list destructors.
 * - Shuts down MPI nodes (in MPI mode).
 * - Closes log files.
 * @return true if the shutdown was successful, false otherwise.
 */
bool GlobalShutdown();

//_______________________________________________________________________

/**
 * @brief If the settings request it, write a warning (diagonstic) message to the .log file.
 *
 * @param message The diagnostic message to report.
 */
void ReportWarning(_String const &message);

//_______________________________________________________________________

/**
 * @brief Push a warning message to the console (for high impact warnings).
 *
 * @param message The diagnostic message to report.
 */
void ReportWarningConsole(_String const message);

//_______________________________________________________________________

/**
 * @brief This is a simple wrapper that will either store the error message in the provided pointer (for alternative handling), or go through the standard error handling procedure.
 *
 * @param error_string If not null, then store message here, otherwise use standard error handling.
 * @param message The error message.
 */
void HandleOrStoreApplicationError(_String *error_string,
                                   _String const &message);

/**
 * @brief This is a convenience function to report an error while parsing expressions.
 *
 * @param error_string The error message.
 * @param context The context in which the error occurred.
 */
void HandleErrorWhileParsing(_String const &error_string,
                             _String const &context);

/**
 * @brief When HyPhy encounters a fatal error, this function reports the error to the user and tries to exit gracefully.
 *
 * @param message The error message.
 * @param force_exit If true, the application will exit regardless of the current error handling protocols.
 * @param dump_core If true, the current HyPhy status will be dumped to a file.
 * @param minimal_error_reporting If true, do not report call stacks and such (e.g. during assertion handling).
 */
void HandleApplicationError(_String const &message, bool force_exit = false,
                            bool dump_core = false,
                            bool minimal_error_reporting = false);

/**
 * @brief Wrapper around HandleApplicationError which reports an error message and exits the program.
 *
 * @param message The error message.
 */
[[noreturn]] void HandleApplicationErrorAndExit(_String const &message);

/**
 * @brief When HyPhy encounters an error in a particular expression, it may be useful to report the location of the error in the string.
 * This function extracts the substring of the requested size from the context, padding with ellipses if needed, and returns it.
 *
 * @param context The context in which the error occurred.
 * @param from The starting position of the error.
 * @param size The size of the context to extract.
 * @return The prepared error context.
 */
const _String PrepareErrorContext(_String const &context, long from,
                                  unsigned long size = 32UL);

/**
 * @brief Get a path specification to one of the standard HyPhy directories.
 *
 * @param which_one Which directory to return; supported values are presently limited to HY_HBL_DIRECTORY_TEMPLATE_MODELS.
 * @return A path specification to one of the standard HyPhy directories.
 */
const _String GetStandardDirectory(const unsigned long which_one);

/**
 * @brief Resolve the path name contained in 'path', according to various flag settings.
 *
 * @param path The path to resolve.
 * @param isWrite If true, the path is for writing.
 * @param acceptStringVars If true, string variables are accepted.
 * @param p If not nil, a pointer to a `_String` object.
 * @param assume_platform_specific If true, assume platform-specific path.
 * @param caller The caller of the function.
 * @param relative_to_base If true, the path is relative to the base directory.
 * @param relative_path_passthrough If true, relative paths are passed through.
 * @return true if the path was resolved successfully, false otherwise.
 */
bool ProcessFileName(_String &path, bool isWrite = false,
                     bool acceptStringVars = false, hyPointer = nil,
                     bool assume_platform_specific = false,
                     _ExecutionList *caller = nil,
                     bool relative_to_base = false,
                     bool relative_path_passthrough = false);

/**
 * @brief Logs to the console.
 */
void ConsoleLog(void);

/**
 * @brief Get the full HyPhy version string.
 * @return The full HyPhy version string (e.g., 'HYPHY 2.31alpha20170206beta(MP) for Darwin on x86_64').
 */
const _String GetVersionString(void);

/**
 * @brief Get the current timestamp.
 * @param do_gmt If true, return GMT time as YYYY/[M]M/[D]D H:M.
 * @return The current timestamp (e.g., 'Fri Jun 16 22:29:10 2017' or '2017/6/17 5:29').
 */
const _String GetTimeStamp(bool do_gmt = false);

/**
 Global variables (grouped by type, then alphabetically)
 */

/**
 * @brief Constructs an error message.
 *
 * @param s The string to include in the error message.
 * @return _String* The error message.
 */
_String *ConstructAnErrorMessage(_String const &);

extern _AVLList _hy_application_globals;

extern _List _hy_standard_library_extensions, _hy_standard_library_paths;

extern hyFile *hy_error_log_file, *hy_message_log_file;

extern bool hy_drop_into_debug_mode, hy_need_extra_nl, terminate_execution,
    has_terminal_stdout, has_terminal_stderr, ignore_kw_defaults,
    force_verbosity_from_cli, has_terminal_color;

extern int hy_mpi_node_rank, hy_mpi_node_count;

extern long system_CPU_count, print_digit_specification, verbosity_level;

extern unsigned long matrix_exp_count, taylor_terms_count, squarings_count;

extern hyTreeDefinitionPhase isDefiningATree;

extern _String const kHyPhyVersion, kEmptyString, kPromptForFilePlaceholder,
    kTemporaryFilePlaceholder, kEmptyAssociativeList, kHyphyCiteString,
    kXVariableName, kNVariableName, kNamespaceName,
    kErrorStringIncompatibleOperands, // -101
    kErrorStringBadMatrixDefinition,  // -104
    kErrorStringInvalidMatrixIndex,   // -106
    kErrorStringMemoryFail,           // -108
    kErrorStringDatasetRefIndexError, // -171
    kErrorStringMatrixExportError,    // -200
    kErrorStringNullOperand,          // -666
    kErrorNumerical, kErrorStringUnterminatedMatrix, kNoneToken, kNullToken,
    kNoKWMatch, kEndIteration;

extern _String hy_base_directory, hy_error_log_name, hy_lib_directory,
    hy_messages_log_name, hy_standard_library_directory,
    hy_standard_model_directory, hy_scanf_last_file_path;

extern _Trie availableTemplateFilesAbbreviations;

extern _List availableTemplateFiles, availablePostProcessors;

extern _Variable *hy_x_variable, *hy_n_variable;

extern const hyFloat kMachineEpsilon;

extern std::regex *hy_float_regex, *hy_replicate_constraint_regexp;

} // namespace hy_global

#endif
