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
 Global functions and variables that don't naturally belong with classes
 or other headers go here .
 
 */


#ifndef __HYGLOBALTHINGS__
#define __HYGLOBALTHINGS__

#include "defines.h"
#include "hy_types.h"
#include "avllistx.h"
#include "regex.h"

#include <stdio.h>

#ifdef __HYPHYMPI__
#include <mpi.h>
#endif

class _Variable; // forward decl
class _ExecutionList; // forward decl

namespace hy_global {
  
  /**
   Global functions
   */
  
  
  /**
   Heap-allocate a specified number of bytes and return the pointer.
   Optionally zero the memory. If allocation fails, the program is halted
   via a fatal run-time error
   
   @param bytes the number of bytes to allocate
   @param zero if TRUE, zero the memory block
   
   @see FlagError, MemReallocate
   @return a pointer to the new memory block
   
   */
  hyPointer   MemAllocate (long bytes, bool zero = false, size_t alignment = 0);
  
  /**
   Resize an existing pointer to 'new_bytes' bytes.
   If allocation fails, the program is halted
   via a fatal run-time error
   
   @param bytes the number of bytes to allocate
   @param old_pointer a previously allocated (with MemAllocate) pointer
   
   @see FlagError, MemReallocate
   @return a pointer to the resized; could be different from old_pointer
   
   */
  hyPointer   MemReallocate (hyPointer old_pointer, long new_bytes);
  
  /**
   Set up application environment
   
   -   initialize the list indicating which variables always exist in
   the global namespace  (_hy_application_globals)
   
   -   set up the paths where HyPhy will look for its library files _hy_standard_library_paths (via LoadFunctionLibrary), and which
   standard extensions will be used to look for files
   (_hy_standard_library_extensions)
   
   -   initialize parser arrays (_HBL_Init_Const_Arrays)
   
   -   set up some environment variables and language constants, like
   TRUE, FALSE, HYPHY_BASE_DIRECTORY, HYPHY_LIB_DIRECTORY
   
   -   call global set up for hy_env
   */
  void    InitializeGlobals (void);
  
  /**
   Return directory separator character for the
   build platform (e.g., '/' for *nix)
   
   @return the directory separator char
   */
  char get_platform_directory_char (void);
  
  /**
   Open the file located at file_path using mode 'mode'
   
   @param file_path the path of the file to open
   @param mode standard (fopen) open modes
   @param error if true, then halt the program if the file could not be opened
   
   @return the FILE handle or nil (if file does not exist or could not be open with the requested mode)
   */
  FILE*   doFileOpen                (const char * file_path, const char * mode , bool error = false);
  
  /**
   The omnibus clean-up function that attempts to deallocate all application memory
   and unwind various objects created; the idea is that upon successful completion,
   the state of the program is the same as it was after the initial startup
   
   TODO 20170414: does the function ever get called with the FALSE flag? if not, deprecate
   I don't recall what the use case for the FALSE flag was
   @param all if NOT set, a partial purge is understaken, that clears only user functions and attendant stuctures
   
   
   */
  void    PurgeAll                  (bool   all = true);
  
  /**
   This a general initialization function that deals with bookkeeping like
   * Call all other initializers
   * Define and set the global random seed
   * Create and open messages and errors logs (based on settings)
   
   */
  
  bool    GlobalStartup(void);
  /**
   The general clean up function that is called right before HyPhy exits
   * Call other object list destructors
   * Shut down MPI nodes (in MPI mode)
   * close log files
   
   
   */
  bool    GlobalShutdown();
  
  //_______________________________________________________________________
  
  /** If the settings request it, write a warning (diagonstic) message to the
   .log file
   
   @param message the diagnostic message to report
   */
  void    ReportWarning (_String const & message);
  
  //_______________________________________________________________________
  
  /**
   This is a simple wrapper that will either store the error message in the provided
   pointer (for alternative handling), or go through the standard error handling procedure
   
   @param error_string if not null, then store message here, otherwise use standard error handling
   @message the error message
   
   */
  void      HandleOrStoreApplicationError (_String* error_string, _String const & message);
  
  /**
   This is a convenience function to report an error while parsing expressions (context)
   */
  void        HandleErrorWhileParsing (_String const & error_string, _String const & context);
  
  
  /**
   When HyPhy encounters a fatal error, this function reports the
   error to the user and tries to exit gracefully
   
   If the flag is set, then the error handling protocols for the current execution
   context is ignored, and the application exits
   
   @param message the error message
   @param force_exit force application exit
   */
  void      HandleApplicationError (_String const & message, bool force_exit = false, bool dump_core = false);
  
  /**
      When HyPhy encounters an error in a particular expression, it may be useful
      to report the location of the error in the string.
   
      This function extracts the substring of the requested size from the context,
      padding with ellipses if needed, and returns it
   */
   
   const _String  PrepareErrorContext     (_String const & context, long from, unsigned long size = 32UL);

  /**
   Return a path specification to one of the standard HyPhy directories
   
   
   @param which_one which directory to return; supported values are presently limited to HY_HBL_DIRECTORY_TEMPLATE_MODELS
   @return a path specification to one of the standard HyPhy directories
   */
  const _String      GetStandardDirectory (const unsigned long which_one);
  
  /** resolve the path name contained in 'path', accoding to various flag settings
   
   
   
   * Revision history
   - SLKP 20170616; moved from _String, stripped much of the legacy code, basically only the UNIX branch is needed now
   
   
   */
  
  
  bool  ProcessFileName (_String & path, bool isWrite = false, bool acceptStringVars = false, hyPointer = nil, bool assume_platform_specific = false, _ExecutionList * caller = nil, bool relative_to_base = false, bool relative_path_passthrough = false);
  
  void    ConsoleLog (void);
  
  /** Return the full HyPhy version string
   e.g 'HYPHY 2.31alpha20170206beta(MP) for Darwin on x86_64'
   * Revision history
   - SLKP 20170616; moved from _String
   */
  
  const _String GetVersionString (void);
  
  /** Return the full HyPhy version string
   Example: 'Fri Jun 16 22:29:10 2017' or '2017/6/17 5:29' (if do_gmt is set to true)
   @param do_gmt return GMT time as YYYY/[M]M/[D]D H:M
  
   * Revision history
   - SLKP 20170616; moved from _String
   */
  const _String GetTimeStamp (bool do_gmt = false);
  
  /**
   Global variables (grouped by type, then alphabetically)
   */
  
  extern  _AVLList  _hy_application_globals;
  
  extern  _List     _hy_standard_library_extensions,
  _hy_standard_library_paths;
  
  extern  FILE*     hy_error_log_file,
  *     hy_message_log_file;
  
  
  extern  bool    hy_drop_into_debug_mode,
  hy_need_extra_nl,
  terminate_execution,
  has_terminal_stdout,
  has_terminal_stderr,
  ignore_kw_defaults,
  force_verbosity_from_cli;
  
  extern  int      hy_mpi_node_rank,
  hy_mpi_node_count;
  
  extern  long     system_CPU_count,
                   print_digit_specification,
                   verbosity_level;

    
  extern unsigned long matrix_exp_count,
                       taylor_terms_count,
                       squarings_count;
    
  
  extern   hyTreeDefinitionPhase isDefiningATree;
  
  extern  _String  const kHyPhyVersion,
  kEmptyString,
  kPromptForFilePlaceholder,
  kTemporaryFilePlaceholder,
  kEmptyAssociativeList,
  kHyphyCiteString,
  kXVariableName,
  kNVariableName,
  kErrorStringIncompatibleOperands,        // -101
  kErrorStringBadMatrixDefinition,         // -104
  kErrorStringInvalidMatrixIndex,          // -106
  kErrorStringMemoryFail,                  // -108
  kErrorStringDatasetRefIndexError,        // -171
  kErrorStringMatrixExportError,           // -200
  kErrorStringNullOperand,                  // -666
  kErrorStringUnterminatedMatrix,
  kNoneToken,
  kNullToken,
  kNoKWMatch,
  kEndIteration;
  
  extern  _String  hy_base_directory,
  hy_error_log_name,
  hy_lib_directory,
  hy_messages_log_name,
  hy_standard_library_directory,
  hy_standard_model_directory,
  hy_scanf_last_file_path;
  
  extern _Variable * hy_x_variable, * hy_n_variable;
    
  extern const hyFloat     kMachineEpsilon;
  
  extern regex_t   * hy_float_regex,
                   * hy_replicate_constraint_regexp;
  
}

#endif
