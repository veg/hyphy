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

#include <stdio.h>

#ifdef __HYPHYMPI__
    #include <mpi.h>
#endif

class _Variable; // forward decl

namespace hy_global {
  
  /**
   Global functions
  */
  
    
  hy_pointer   MemAllocate (long bytes, bool zero = false);
    /**
     Heap-allocate a specified number of bytes and return the pointer.
     Optionally zero the memory. If allocation fails, the program is halted
     via a fatal run-time error

     @param bytes the number of bytes to allocate
     @param zero if TRUE, zero the memory block
     
     @see FlagError, MemReallocate
     @return a pointer to the new memory block
     
  */
    
  hy_pointer   MemReallocate (hy_pointer old_pointer, long new_bytes);
    /**
     Resize an existing pointer to 'new_bytes' bytes.
     If allocation fails, the program is halted
     via a fatal run-time error
     
     @param bytes the number of bytes to allocate
     @param old_pointer a previously allocated (with MemAllocate) pointer
     
     @see FlagError, MemReallocate
     @return a pointer to the resized; could be different from old_pointer
     
     */
    
    void    InitializeGlobals (void);
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

    char get_platform_directory_char (void);
    /**
       Return directory separator character for the 
       build platform (e.g., '/' for *nix)
     
       @return the directory separator char
    */
    
    FILE*   doFileOpen                (const char * file_path, const char * mode , bool error = false);
    /** 
        Open the file located at file_path using mode 'mode' 
        
        @param file_path the path of the file to open
        @param mode standard (fopen) open modes
        @param error if true, then halt the program if the file could not be opened
     
        @return the FILE handle or nil (if file does not exist or could not be open with the requested mode)
     */

    void    PurgeAll                  (bool   all = true);
    /** 
        The omnibus clean-up function that attempts to deallocate all application memory
        and unwind various objects created; the idea is that upon successful completion,
        the state of the program is the same as it was after the initial startup
     
        TODO 20170414: does the function ever get called with the FALSE flag? if not, deprecate
        I don't recall what the use case for the FALSE flag was
        @param all if NOT set, a partial purge is understaken, that clears only user functions and attendant stuctures 
     
     
     */
    
    bool    GlobalStartup(void);
    /**
        This a general initialization function that deals with bookkeeping like 
            * Call all other initializers
            * Define and set the global random seed
            * Create and open messages and errors logs (based on settings)
     
     */
    bool    GlobalShutdown();
    /** 
        The general clean up function that is called right before HyPhy exits
            * Call other object list destructors
            * Shut down MPI nodes (in MPI mode)
            * close log files
     
     
     */
    
    //_______________________________________________________________________
    void    ReportWarning (_String const message);
    
    /** If the settings request it, write a warning (diagonstic) message to the 
        .log file 
     
        @param message the diagnostic message to report
     */
     
    //_______________________________________________________________________
    void      HandleOrStoreApplicationError (_String* error_string, _String const message);
    /**
        This is a simple wrapper that will either store the error message in the provided
        pointer (for alternative handling), or go through the standard error handling procedure
     
        @param error_string if not null, then store message here, otherwise use standard error handling
        @message the error message
     
     */
    
    //_______________________________________________________________________
    void        HandleErrorWhileParsing (_String const error_string, _String const context);
    /**
        This is a convenience function to report an error while parsing expressions (context)
     */

    
    void      HandleApplicationError (_String const message, bool force_exit = false);
    /** 
        When HyPhy encounters a fatal error, this function reports the 
        error to the user and tries to exit gracefully
     
        If the flag is set, then the error handling protocols for the current execution
        context is ignored, and the application exits
     
        @param message the error message
        @param force_exit force application exit
     */
    
    const _String      GetStandardDirectory (const unsigned long which_one);
    /**
     Return a path specification to one of the standard HyPhy directories
     
     
     @param which_one which directory to return; supported values are presently limited to HY_HBL_DIRECTORY_TEMPLATE_MODELS
     @return a path specification to one of the standard HyPhy directories
     */

    
    void    ConsoleLog (void);


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
                    terminate_execution;

    extern  int      hy_mpi_node_rank,
                     hy_mpi_node_count;

    extern  long     system_CPU_count;
    
    
    extern  _String  const kEmptyString,
                           kXVariableName,
                           kNVariableName,
                           kErrorStringIncompatibleOperands,        // -101
                           kErrorStringBadMatrixDefinition,         // -104
                           kErrorStringInvalidMatrixIndex,          // -106
                           kErrorStringMemoryFail,                  // -108
                           kErrorStringDatasetRefIndexError,        // -171
                           kErrorStringMatrixExportError,           // -200
                           kErrorStringNullOperand                  // -666
                           ;
    
    extern  _String  hy_base_directory,
                     hy_error_log_name,
                     hy_lib_directory,
                     hy_messages_log_name,
                     hy_standard_library_directory,
                     hy_standard_model_directory;
    
    extern _Variable * hy_x_variable, * hy_n_variable;

}

#endif
