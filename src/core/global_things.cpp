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



#include "global_things.h"
#include "hy_strings.h"
#include "hbl_env.h"
#include "batchlan.h"
#include "mersenne_twister.h"
#include "global_object_lists.h"

#if defined   __UNIX__ 
    #include <unistd.h>
#endif

#include <time.h>

using     namespace hy_env;

extern _SimpleList freeSlots;
extern  _HY_TREE_DEFINITION_PHASE       isDefiningATree;
extern  _String                         hy_scanf_last_file_path;


/*bool        terminate_execution  = false;
 
 
 //____________________________________________________________________________________
 
 
 extern    _SimpleList   freeSlots;
 
 
 */

namespace hy_global {
    
    
    
    bool             hy_drop_into_debug_mode = false,
        /** if set (e.g., via a -d command line flag, errors will drop
            HyPhy into a debug mode, i.e., an interactive console mode
         */
                     hy_need_extra_nl = false,
        /** set to true if it is necessary to print an extra new line character 
            when HyPhy exits */
    
                     terminate_execution = false;
    
    FILE            *hy_error_log_file,
                    *hy_message_log_file;
    
    _String const    kEmptyString,
                     kXVariableName ("hy_x_variable"),
                     kNVariableName ("_n_"),
                     kErrorStringIncompatibleOperands ("Incompatible operands"),
                     kErrorStringBadMatrixDefinition  ("Invalid matrix definition"),
                     kErrorStringInvalidMatrixIndex   ("Invalid matrix index"),
                     kErrorStringMemoryFail           ("Out of memory"),
                     kErrorStringDatasetRefIndexError ("Dataset index reference out of range"),
                     kErrorStringMatrixExportError    ("Export matrix called with a non-polynomial matrix argument"),
                     kErrorStringNullOperand          ("Attempting to operate on an undefined value; this is probably the result of an earlier 'soft' error condition");
    
    _String
                     hy_base_directory,
                     hy_lib_directory,
                     hy_standard_library_directory ("TemplateBatchFiles"),
                     hy_standard_model_directory   ("TemplateModels"),
                     hy_error_log_name             ("errors.log"),
                     hy_messages_log_name          ("messages.log");
    
    
    _List            _hy_application_globals_aux,
                     _hy_standard_library_paths,
                     _hy_standard_library_extensions;
    
    _AVLList         _hy_application_globals (&_hy_application_globals_aux);
    
    _Variable*       hy_x_variable = nil,
             *       hy_n_variable = nil;
    
    int              hy_mpi_node_rank,
        // [MPI only] the MPI rank of the current node (0 = master, 1... = slaves)
                     hy_mpi_node_count;
        // [MPI only]  MPI node count on the system
    
    long             system_CPU_count = 1L;
        /** the number of CPUs for OpenMP */
    

    //____________________________________________________________________________________
    
    hy_pointer MemAllocate (long bytes, bool zero) {
        
        hy_pointer result = (hy_pointer) zero ? calloc (bytes, 1) : malloc (bytes);
        
        if (result == nil) {
            HandleApplicationError (_String ("Failed to allocate '")  & bytes & "' bytes'", true);
        }
        return result;
    }
    
    
    //____________________________________________________________________________________
hy_pointer MemReallocate (hy_pointer old_pointer, long new_size)
    {
        hy_pointer result  = (hy_pointer) realloc (old_pointer, new_size);
        
        if (!result) {
            HandleApplicationError (_String ("Failed to resize memory to '")  & new_size & "' bytes'", true);
        }
        
        return result;
    }
    
    
    //____________________________________________________________________________________
    char get_platform_directory_char (void) {
        #ifdef __MAC__
                return ':';
        #elif defined __WINDOZE__ || defined __MINGW32__
                return '\\';
        #endif
                return '/';
    }

    //____________________________________________________________________________________
    
    FILE *      doFileOpen (const char * fileName, const char * mode, bool error) {
        FILE    *daFile = nil;
        
        if (fileName) {
            daFile = fopen (fileName,mode);
            if (!daFile && error) {
                HandleApplicationError (_String("Could not open file '") & *fileName & "' with mode '" & *mode & "'.");
            }
        }
        return daFile;
    }
    //____________________________________________________________________________________
   

    void    InitializeGlobals (void) {
        
        SetupEnvDefaults();
        
        using namespace hy_env;
        _String const * mark_as_globals [] =
        {
            &data_file_tree,
            &data_file_tree_string,
            &sitewise_matrix,
            &blockwise_matrix,
            &random_seed,
            &selection_strings,
            &base_directory,
            &lib_directory,
            &directory_separator_char,
            &path_to_current_bf,
            &true_const,
            &false_const,
            &matrix_element_value,
            &matrix_element_row,
            &matrix_element_column,
            &nexus_file_tree_matrix,
            &dataset_save_memory_size,
            &last_file_path,
            &mpi_node_id,
            &mpi_last_sent_message,
            &mpi_node_count,
            &error_report_format_expression,
            &error_report_format_expression_stack,
            &error_report_format_expression_stdin,
            &error_report_format_expression_string,
            &kXVariableName,
            &kNVariableName,
            &status_bar_update_string
        };
        
        for (_String const * item : mark_as_globals) {
            _hy_application_globals.Insert (new _String (*item));
        }
         
        
        _String             dd (get_platform_directory_char());
        
        _hy_standard_library_paths.AppendNewInstance      (new _String(hy_lib_directory & hy_standard_library_directory & dd ));
        _hy_standard_library_paths.AppendNewInstance      (new _String(hy_lib_directory & hy_standard_library_directory & dd & hy_standard_model_directory & dd ));
        _hy_standard_library_paths.AppendNewInstance      (new _String(hy_lib_directory & hy_standard_library_directory & dd & "Utility" & dd));
        _hy_standard_library_paths.AppendNewInstance      (new _String(hy_lib_directory & "UserAddIns" & dd));
        _hy_standard_library_paths.AppendNewInstance      (new _String(hy_lib_directory & hy_standard_library_directory & dd & "Distances" & dd));
        
        const char * extensions [] = {"", ".bf", ".ibf", ".def", ".mdl"};
        for (const char * ext : extensions ) {
            _hy_standard_library_extensions < new _String (ext);
        }
        
        //StringToConsole(*((_String*)_hy_standard_library_paths.toStr())); NLToConsole();
        //StringToConsole(*((_String*)_hy_standard_library_extensions.toStr())); NLToConsole();
        
        _HBL_Init_Const_Arrays  ();
        //_Variable* CheckReceptacle (_String const * name, _String const & fID, bool checkValid, bool isGlobal)

        
        CheckReceptacleAndStore(&true_const,  kEmptyString, false, new _Constant (1.), false);
        CheckReceptacleAndStore(&false_const, kEmptyString, false, new _Constant (0.), false);
        hy_x_variable = CheckReceptacle(&kXVariableName, kEmptyString, false, false);
        hy_n_variable = CheckReceptacle(&kNVariableName, kEmptyString, false, false);
        
        EnvVariableSet(directory_separator_char, new _FString (dd, false), false);
        EnvVariableSet(base_directory, new _FString (hy_base_directory, false), false);
        EnvVariableSet(lib_directory,  new _FString (hy_lib_directory, false), false);
        
    }
    //____________________________________________________________________________________
    
    void    PurgeAll (bool all) {
        
            // TODO 20170414: this is ugly; not sure where this datastructure belongs
        
        using namespace hyphy_global_objects;
        
        ClearBFFunctionLists();
        executionStack.Clear();
        loadedLibraryPaths.Clear(true);
        _HY_HBL_Namespaces.Clear();
        if (all) {
            ClearAllGlobals ();
            likeFuncList.Clear();
            likeFuncNamesList.Clear();
            dataSetList.Clear();
            dataSetNamesList.Clear();
            compiledFormulaeParameters.Clear();
            modelNames.Clear();
            KillExplicitModelFormulae ();
            modelMatrixIndices.Clear();
            modelFrequenciesIndices.Clear();
            modelTypeList.Clear();
            listOfCompiledFormulae.Clear();
            variablePtrs.Clear();
            freeSlots.Clear();
            lastMatrixDeclared = -1;
            variableNames.Clear(true);
            _hy_application_globals.Clear(true);
            
            hy_x_variable = nil;
            hy_n_variable = nil;
            pathNames.Clear();
        }
        ::hy_scanf_last_file_path = kEmptyString;
        EnvVariableSet(random_seed, new _Constant (hy_random_seed), false);
        isInFunction        = _HY_NO_FUNCTION;
        isDefiningATree     = kTreeNotBeingDefined;
#ifdef __HYPHYMPI__
        MPI_Comm_size   (MPI_COMM_WORLD, &hy_mpi_node_count);
        EnvVariableSet  (mpi_node_count, new _Constant (hy_mpi_node_count), false);
        EnvVariableSet  (mpi_node_id, new _Constant (hy_mpi_node_rank), false);
#endif
    }

    
    //____________________________________________________________________________________
    bool    GlobalStartup (void) {
        SetupOperationLists     ();
        time_t                  k;
        time                    (&k);
        hy_random_seed          = k + getpid ();
        
        InitializeGlobals();
       
        init_genrand            (hy_random_seed);
        EnvVariableSet(random_seed, new _Constant (hy_random_seed), false);
        
        
#if not defined (__HYPHY_MPI_MESSAGE_LOGGING__) && defined (__HYPHYMPI__)
        if (hy_mpi_node_rank == 0L) {
#endif
            
            
#ifndef __HEADLESS__ // do not create log files for _HEADLESS_
            _String * prefix [2] = {&hy_error_log_name, &hy_messages_log_name};
            FILE ** handle [2] = {&hy_error_log_file, &hy_message_log_file};
            
            for (long file_index = 0; file_index < 2; file_index++) {
                long                    p   = 1L;
                
    #ifndef __HYPHYMPI__
                _String file_name (*prefix[file_index]);
        #if defined  __MINGW32__
                    file_name = hy_base_directory & file_name;
        #endif
    #else // MPI branch
                _String file_name = *prefix[file_index] & ".mpinode" & (long)hy_mpi_node_rank;
    #endif
                
                *handle[file_index] = doFileOpen (file_name.sData,"w+");
                while (*handle[file_index] == nil && p<10) {
                    #ifndef __HYPHYMPI__
                         file_name = *prefix[file_index] & '.' & p;
                    #if defined  __MINGW32__
                        file_name = hy_base_directory & file_name;
                    #endif
                    #else // MPI branch
                        file_name = *prefix[file_index] & ".mpinode" & (long)hy_mpi_node_rank & '.' & p;
                    #endif
                    p++;
                }
                *prefix[file_index] = file_name;
            }
#endif

#if not defined (__HYPHY_MPI_MESSAGE_LOGGING__) && defined (__HYPHYMPI__)
        }
#endif
        
        return hy_error_log_file && hy_message_log_file;
    }
    
    
    //____________________________________________________________________________________
    bool    GlobalShutdown (void) {
        bool no_errors = true;
        
#if defined __UNIX__
        if (hy_need_extra_nl) {
            printf ("\n");
        }
#endif
        
#ifdef  __HYPHYMPI__
        int     flag;
        
        PurgeAll (true);
        ReportWarning ("PurgeAll was successful");
        if (hy_mpi_node_rank == 0) {
            fflush (stdout);
            
            for (long count = 1; count < hy_mpi_node_count; count++) {
                ReportWarning (_String ("Sending shutdown command to node ") & count & '.');
                MPISendString (kEmptyString,count);
            }
        }
#else
        fflush (stdout);
#endif
        
        
        for (AVLListXIteratorKeyValue command_helper : AVLListXIterator (&_HY_HBLCommandHelper)) {
            delete ((_HBLCommandExtras*)command_helper.get_value());
        }
        
        
        PurgeAll(true);
        
        _HY_HBLCommandHelper.Clear();
        _HY_ValidHBLExpressions.Clear();
        listOfCompiledFormulae.Clear();
        
#ifdef  __HYPHYMPI__
        // MPI_Barrier (MPI_COMM_WORLD);
        ReportWarning ("Calling MPI_Finalize");
    #ifdef __USE_ABORT_HACK__
            MPI_Abort(MPI_COMM_WORLD,0);
    #else
            MPI_Finalized(&flag);
            if (!flag)
                MPI_Finalize();
    #endif
        ReportWarning ("Returned from MPI_Finalize");
#endif
        
        
        _String * prefix [2] = {&hy_error_log_name, &hy_messages_log_name};
        char const * messages [] = {"\nCheck %s for execution error details.\n", "\nCheck %s for diagnostic messages.\n"};
        FILE * handle [2] = {hy_error_log_file, hy_message_log_file};
        
        for (long file_index = 0; file_index < 2; file_index++) {
       
            if (handle[file_index]) {
                fflush (handle[file_index]);
                fseek(handle[file_index],0,SEEK_END);
                unsigned long fileSize = ftell(handle[file_index]);
                if (fileSize) {
                    fprintf (stderr, messages[file_index], prefix[file_index]->getStr());
                    if (file_index == 0) { no_errors = false; }
                    fclose (handle[file_index]);
                    
                } else {
                    fclose (handle[file_index]);
                    remove (prefix[file_index]->getStr());
                }
            }
        }
        
        return no_errors;
    }
    
    //____________________________________________________________________________________
    
    void PopulateCurrentCallStack      (_List& calls, _List& stdins) {
    /**
        Fill out the list of execution lists (calls) and standard inputs (stdins)
     
     */
        calls.Clear  ();
        stdins.Clear ();
        
        if (executionStack.lLength) {
            for (long callLevel = executionStack.lLength -1L ; callLevel >=0L ; callLevel --) {
                _ExecutionList * currentLevel = (_ExecutionList*)executionStack (callLevel);
                
                calls.AppendNewInstance(new _String ((_String*)((_ElementaryCommand*)(*currentLevel)(currentLevel->currentCommand?currentLevel->currentCommand-1:0))->toStr()));
                
                if (currentLevel->stdinRedirect) {
                    stdins.AppendNewInstance ((_String*)currentLevel->stdinRedirect->toStr());
                } else {
                    stdins.AppendNewInstance (new _String);
                }
            }
        }
    }
    
    
    //____________________________________________________________________________________
    
    _String* ConstructAnErrorMessage         (_String const& theMessage) {
    /** Prepare the error message using the current formatting specification **/
     
        _String* error_message = new _String (128L,true);
        
        _List    calls,
                 stdins;
        
        PopulateCurrentCallStack	 (calls, stdins);
        
        _FString * error_formatting_expression = (_FString*) EnvVariableGet(error_report_format_expression, STRING);
        
        bool     doDefault = true;
        
        if (error_formatting_expression) {
            _Formula expression;
            _String  expr (*error_formatting_expression->theString),
            errMsgLocal;
            _FormulaParsingContext fpc (&errMsgLocal, nil);
            
            if (Parse    (&expression, expr, fpc, nil) == HY_FORMULA_EXPRESSION) {
                EnvVariableSet(error_report_format_expression_string, new _FString (theMessage, false), false);
                EnvVariableSet(error_report_format_expression_stack,  new _Matrix (calls), false);
                EnvVariableSet(error_report_format_expression_string, new _Matrix (stdins, false), false);
                
                _PMathObj expr = expression.Compute();
                if (!terminate_execution && expr && expr->ObjectClass() == STRING) {
                    (*error_message) << ((_FString*)expr)->theString;
                    doDefault = false;
                }
            }
        }
        
        if (doDefault) {
            (*error_message) << "Error:\n"
                      << theMessage;
            
            if (calls.lLength) {
                (*error_message) << "\n\nFunction call stack\n";
                for (unsigned long k = 0UL; k < calls.lLength; k++) {
                    (*error_message) << (_String((long)k+1) & " : " & (*(_String*)calls(k)) & '\n');
                    _String* redir = (_String*)stdins (k);
                    if (redir->sLength) {
                        (*error_message) << "\tStandard input redirect:\n\t\t"
                                  << redir->Replace ("\n","\n\t\t",true);
                    }
                    (*error_message) << "-------\n";
                    
                }
            }
        }
        
        error_message->Finalize();
        return error_message;
    }
    
    //____________________________________________________________________________________
    const _String    GetStandardDirectory (const unsigned long which_one) {
        _String dirSpacer (get_platform_directory_char());
        
        switch (which_one) {
                
            case HY_HBL_DIRECTORY_TEMPLATE_MODELS:
                return hy_lib_directory & hy_standard_library_directory & dirSpacer & hy_standard_model_directory & dirSpacer;
        }
        
        return kEmptyString;
    }
    
    //____________________________________________________________________________________
    void    ReportWarning (_String const message) {
        
        bool do_logging = EnvVariableTrue(message_logging);
            
#ifdef  __HEADLESS__
        if (globalInterfaceInstance && do_logging >= 0.1) {
            globalInterfaceInstance->PushWarning (&message);
        }
#else
        if ( ! (hy_message_log_file && do_logging )) {
            return;
        }
        
        char   str[] = "\n";
        fwrite (str, 1, 1, hy_message_log_file);
        fwrite (message.getStr(), 1, message.Length(), hy_message_log_file);
        fflush (hy_message_log_file);
#endif
    }
    
    //____________________________________________________________________________________
    void HandleOrStoreApplicationError (_String* error_string, _String const message) {
        if (error_string) {
            *error_string = message;
        } else {
            HandleApplicationError (message);
        }
    }

    //____________________________________________________________________________________
    void HandleErrorWhileParsing (_String const error_string, _String const context) {
        HandleApplicationError (_String ("While parsing:\n") & context & "\n" & error_string);
    }

    
    //____________________________________________________________________________________
    void HandleApplicationError (const _String message, bool force_exit) {

        if (!force_exit && currentExecutionList && currentExecutionList->errorHandlingMode == HY_BL_ERROR_HANDLING_SOFT) {
            currentExecutionList->ReportAnExecutionError(message, true);
            return;
        }
            
    #ifdef  __HEADLESS__
        if (globalInterfaceInstance) {
            globalInterfaceInstance->PushError (&st);
        }
        terminate_execution = true;
    #else
        char  str[] = "\nError:";
        
        if (hy_error_log_file) {
            fwrite (str, 1, 7, hy_error_log_file);
            fwrite (message.getStr(), 1, message.Length(), hy_error_log_file);
            fflush (hy_error_log_file);
        }
            
        if (hy_error_log_file) {
            fprintf (hy_error_log_file, "\n%s", message.sData);
        }
            
        _String errMsg;
    #ifdef __HYPHYMPI__
        errMsg = _String("Received an error state from MPI node ") & (long)hy_mpi_node_rank & '\n' & message;
        
        if (hy_mpi_node_rank > 0) {
            errMsg = ConstructAnErrorMessage (message);
            fprintf (stderr, "HYPHYMPI terminated.\n%s\n\n", errMsg.sData);
            MPI_Abort (MPI_COMM_WORLD,1);
            abort   ();
        } else {
            errMsg = _String ("\nMaster node received an error:") & message;
        }
    #else
            errMsg = message;
    #endif
            errMsg = ConstructAnErrorMessage (errMsg);
            StringToConsole(errMsg);
    #endif 
        
    #if defined __UNIX__
            if (hy_drop_into_debug_mode)
                while (ExpressionCalculator()) ;
    #endif
    #ifdef __HYPHYMPI__
            if (hy_mpi_node_rank==0) {
                MPI_Abort (MPI_COMM_WORLD,1);
            }
    #endif
            GlobalShutdown();
        
            
#ifdef _HY_ABORT_ON_ERROR
            abort ();
#else
            exit(1);
#endif
    }
} // namespace close
