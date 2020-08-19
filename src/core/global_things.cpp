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
#include "hy_string_buffer.h"
#include "hbl_env.h"
#include "batchlan.h"
#include "mersenne_twister.h"
#include "global_object_lists.h"

#if defined   __UNIX__ 
    #include <unistd.h>
    #include <sys/stat.h>
    #include <sys/types.h>
    #if !defined __MINGW32__
        #include <sys/utsname.h>
    #endif
#endif

#include <time.h>
#include <float.h>
#include <signal.h>

//#define __HYPHY_MPI_MESSAGE_LOGGING__

using     namespace hy_env;

extern _SimpleList freeSlots;


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
    
                     terminate_execution = false,
                     has_terminal_stdout = true,
                     has_terminal_stderr = true,
                     ignore_kw_defaults  = false,
                     force_verbosity_from_cli = false;
    
    FILE            *hy_error_log_file = NULL,
                    *hy_message_log_file = NULL;
    
    hyTreeDefinitionPhase
                     isDefiningATree = kTreeNotBeingDefined;
    
    const hyFloat    kMachineEpsilon = 2.*DBL_EPSILON;

    
    _String const    kEmptyString,
                     kPromptForFilePlaceholder        ("PROMPT_FOR_FILE"),
                     kTemporaryFilePlaceholder        ("TEMP_FILE_NAME"),
                     kEmptyAssociativeList ("{}"),
                     kHyPhyCiteString ("\nPlease cite S.L. Kosakovsky Pond, S. D. W. Frost and S.V. Muse. (2005) HyPhy: hypothesis testing using\
 phylogenies. Bioinformatics 21: 676-679 if you use HyPhy in a publication\nIf you are a new HyPhy user, the tutorial located at\
 http://www.hyphy.org/docs/HyphyDocs.pdf may be a good starting point.\n"),
                     kXVariableName ("_x_"),
                     kNVariableName ("_n_"),
                     kErrorStringIncompatibleOperands ("Incompatible operands"),
                     kErrorStringBadMatrixDefinition  ("Invalid matrix definition "),
                     kErrorStringInvalidMatrixIndex   ("Invalid matrix index "),
                     kErrorStringUnterminatedMatrix   ("Unterminated matrix definition "),
                     kErrorStringMemoryFail           ("Out of memory"),
                     kErrorStringDatasetRefIndexError ("Dataset index reference out of range"),
                     kErrorStringMatrixExportError    ("Export matrix called with a non-polynomial matrix argument"),
                     kErrorStringNullOperand          ("Attempting to operate on an undefined value; this is probably the result of an earlier 'soft' error condition"),
                     kHyPhyVersion  = _String ("2.5.17"),
    
                    kNoneToken = "None",
                    kNullToken = "null",
                    kNoKWMatch = "__input_value_not_given__",
                    kEndIteration = "__iterator_end__loop__";
  
    _String
                     hy_base_directory,
                     hy_lib_directory,
                     hy_scanf_last_file_path,
                     hy_standard_library_directory ("TemplateBatchFiles"),
                     hy_standard_model_directory   ("TemplateModels"),
                     hy_error_log_name             ("errors.log"),
                     hy_messages_log_name          ; // 20200222 : DISABLE MESSAGES BY DEFAULT
  
  
    _List            _hy_application_globals_aux,
                     _hy_standard_library_paths,
                     _hy_standard_library_extensions;
  
    _AVLList         _hy_application_globals (&_hy_application_globals_aux);
  
    _Variable*       hy_x_variable = nil,
             *       hy_n_variable = nil;

    unsigned long    matrix_exp_count,
                     taylor_terms_count,
                     squarings_count;
    
    int              hy_mpi_node_rank,
        // [MPI only] the MPI rank of the current node (0 = master, 1... = slaves)
                     hy_mpi_node_count;
        // [MPI only]  MPI node count on the system
  
    long             system_CPU_count = 1L,
        /** the number of CPUs for OpenMP */
                     print_digit_specification = 0L,
                     verbosity_level = 0L;
  
  
    int _reg_exp_err_code = 0;
    regex_t * hy_float_regex                 = _String::PrepRegExp ("^\\ *[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?", _reg_exp_err_code, true),
            * hy_replicate_constraint_regexp = _String::PrepRegExp ("^this([0-9]+)\\.(.+)$", _reg_exp_err_code, true);
  
  

    //____________________________________________________________________________________
  
    hyPointer MemAllocate (long bytes, bool zero, size_t alignment) {
        hyPointer result = nil;
        
//#ifdef _ISOC11_SOURCE
        if (alignment > 0 && bytes >= alignment && bytes % alignment == 0) {
            if (posix_memalign (&result, alignment, bytes) == 0) {
                if (zero) {
                    memset (result, 0, bytes);
                }
                return result;
            }
        }
//#endif
        result = (hyPointer) zero ? calloc (bytes, 1) : malloc (bytes);
        
        if (result == nil) {
            HandleApplicationError (_String ("Failed to allocate '")  & bytes & "' bytes'", true);
        }
        return result;
    }
    
    
    //____________________________________________________________________________________
  
    hyPointer MemReallocate (hyPointer old_pointer, long new_size) {
        hyPointer result  = (hyPointer) realloc (old_pointer, new_size);
        
        if (result == nil) {
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
            &status_bar_update_string,
            &last_model_parameter_list,
            &use_last_model
        };
        
        for (_String const * item : mark_as_globals) {
            _hy_application_globals.Insert (new _String (*item));
        }
         
        
        _String             dd (get_platform_directory_char());
        
        _hy_standard_library_paths.AppendNewInstance      (new _String(hy_lib_directory & hy_standard_library_directory & dd ));
        _hy_standard_library_paths.AppendNewInstance      (new _String(hy_lib_directory & hy_standard_library_directory & dd & hy_standard_model_directory & dd ));
        _hy_standard_library_paths.AppendNewInstance      (new _String(hy_lib_directory & hy_standard_library_directory & dd & "Utility" & dd));
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
        
        executionStack.Clear();
        loadedLibraryPaths.Clear(true);
        _HY_HBL_Namespaces.Clear();
        if (all) {
            ClearAllGlobals ();
            likeFuncList.Clear();
            likeFuncNamesList.Clear();
            dataSetList.Clear();
            dataSetNamesList.Clear();
            ClearBFFunctionLists();
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
            bgmList.Clear();
            scfgList.Clear();
            bgmNamesList.Clear();
            scfgNamesList.Clear();
            
            hy_x_variable = nil;
            hy_n_variable = nil;
            pathNames.Clear();
        }
        hy_scanf_last_file_path = kEmptyString;
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
        
        _Constant::free_slots.Populate ((long)_HY_CONSTANT_PREALLOCATE_SLOTS, (long)_HY_CONSTANT_PREALLOCATE_SLOTS-1, -1L);
        _StringBuffer::free_slots.Populate ((long)_HY_STRING_BUFFER_PREALLOCATE_SLOTS, (long)_HY_STRING_BUFFER_PREALLOCATE_SLOTS-1, -1L);

        
#ifdef __HYPHYMPI__
        hy_env :: EnvVariableSet (hy_env::mpi_node_id, new _Constant (hy_mpi_node_rank), false);
        hy_env :: EnvVariableSet (hy_env::mpi_node_count, new _Constant (hy_mpi_node_count), false);
        has_terminal_stdout = false;
        has_terminal_stderr = false;
 #else
        has_terminal_stdout = isatty (STDOUT_FILENO);
        has_terminal_stderr = isatty (STDERR_FILENO);
#endif

#if not defined (__HYPHY_MPI_MESSAGE_LOGGING__) && defined (__HYPHYMPI__)
        if (hy_mpi_node_rank == 0L) {
            struct stat sb;
            fstat (STDERR_FILENO, &sb);
            has_terminal_stderr = (sb.st_mode & S_IFMT) == S_IFIFO;
#endif
        
            
#ifndef __HEADLESS__ // do not create log files for _HEADLESS_
            _String * prefix [2] = {&hy_error_log_name, &hy_messages_log_name};
            FILE ** handle [2] = {&hy_error_log_file, &hy_message_log_file};
            
            for (long file_index = 0; file_index < 2; file_index++) {
                long                    p   = 1L;
                
                if (prefix[file_index]->empty()) continue;
    #ifndef __HYPHYMPI__
                _String file_name (*prefix[file_index]);
        #if defined  __MINGW32__
                    file_name = hy_base_directory & file_name;
        #endif
    #else // MPI branch
                _String file_name = *prefix[file_index] & ".mpinode" & (long)hy_mpi_node_rank;
    #endif
                
                *handle[file_index] = doFileOpen (file_name.get_str(),"w+");
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
        
        if (hy_env::cli_env_settings.nonempty()) {
            _ExecutionList (hy_env::cli_env_settings).Execute();
        }
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
            //printf ("Deleting %d\n", command_helper.get_index());
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
                    fprintf (stderr, messages[file_index], prefix[file_index]->get_str());
                    if (file_index == 0) { no_errors = false; }
                    fclose (handle[file_index]);
                    
                } else {
                    fclose (handle[file_index]);
                    remove (prefix[file_index]->get_str());
                }
            }
        }
        
        /*if (_Constant::preallocated_buffer) {
            free ((void*)_Constant::preallocated_buffer);
        }
        if (_StringBuffer::preallocated_buffer) {
            free ((void*)_StringBuffer::preallocated_buffer);
        }*/
        return no_errors;
    }
    
    //____________________________________________________________________________________
    
    void PopulateCurrentCallStack      (_List& calls, _List& stdins, _List& contexts, _List& kwargs) {
    /**
        Fill out the list of execution lists (calls) and standard inputs (stdins)
     
     */
        calls.Clear  ();
        stdins.Clear ();
        contexts.Clear ();
        
        if (executionStack.lLength) {
            for (long callLevel = executionStack.lLength -1L ; callLevel >=0L ; callLevel --) {
                _ExecutionList * currentLevel = (_ExecutionList*)executionStack (callLevel);
                
                calls.AppendNewInstance(new _String ((_String*)((_ElementaryCommand*)(*currentLevel)(currentLevel->currentCommand?currentLevel->currentCommand-1:0))->toStr()));
                
                if (currentLevel->nameSpacePrefix) {
                    contexts << currentLevel->nameSpacePrefix->GetName();
                } else {
                    contexts.AppendNewInstance (new _String);
                }
                
                if (currentLevel->stdinRedirect) {
                    stdins.AppendNewInstance ((_String*)currentLevel->stdinRedirect->toStr());
                } else {
                    stdins.AppendNewInstance (new _String);
                }
                
                if (currentLevel->kwargs && currentLevel->kwargs->countitems()) {
                    kwargs.AppendNewInstance ((_String*)currentLevel->kwargs->toStr());
                }else {
                    kwargs.AppendNewInstance (new _String);
                }
            }
        }
    }
    
    
    //____________________________________________________________________________________
    
    _String* ConstructAnErrorMessage         (_String const& theMessage) {
    /** Prepare the error message using the current formatting specification **/
     
        _StringBuffer* error_message = new _StringBuffer (128UL);
        
        _List    calls,
                 stdins,
                 contexts,
                 kwargs;
        
        
        PopulateCurrentCallStack	 (calls, stdins, contexts, kwargs);
        
        _FString * error_formatting_expression = (_FString*) EnvVariableGet(error_report_format_expression, STRING);
        
        bool     doDefault = true;
        
        if (error_formatting_expression) {
            _Formula expression;
            _String  expr (error_formatting_expression->get_str()),
            errMsgLocal;
            _FormulaParsingContext fpc (&errMsgLocal, nil);
            
            if (Parse    (&expression, expr, fpc, nil) == HY_FORMULA_EXPRESSION) {
                EnvVariableSet(error_report_format_expression_string, new _FString (theMessage, false), false);
                EnvVariableSet(error_report_format_expression_stack,  new _Matrix (calls), false);
                EnvVariableSet(error_report_format_expression_string, new _Matrix (stdins, false), false);
                
                HBLObjectRef expr = expression.Compute();
                if (!terminate_execution && expr && expr->ObjectClass() == STRING) {
                    (*error_message) << ((_FString*)expr)->get_str();
                    doDefault = false;
                }
            }
        }
        
        if (doDefault) {
            (*error_message) << "Error:\n" << theMessage;
            
            if (calls.nonempty()) {
                (*error_message) << "\n\nFunction call stack\n";
                for (unsigned long k = 0UL; k < calls.countitems(); k++) {
                    (*error_message) << (_String((long)k+1)) << " :  ";
                    _String * context = (_String*)contexts.GetItem(k);
                    if (context->nonempty()) {
                        (*error_message) << "[namespace = " << *context << "] ";
                    }
                    (*error_message) << (*(_String*)calls(k)) << '\n';
                    
                    if (k == 0) {
                        _String* redir = (_String*)stdins (k);
                        if (redir->nonempty()) {
                            (*error_message) << "\tStandard input redirect:\n\t\t"
                                      << redir->Replace ("\n","\n\t\t",true);
                        }
                        redir = (_String*)kwargs(k);
                        if (redir->nonempty()) {
                            (*error_message) << "\n\tKeyword arguments:\n\t\t"
                            << redir->Replace ("\n","\n\t\t",true);
                        }
                        (*error_message) << "\n";
                    }
                    
                    (*error_message) << "-------\n";
                    
                }
            } else {
                (*error_message) << "\n";
            }
        }
        
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
    void    ReportWarning (_String const & message) {
        
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
      fwrite (message.get_str(), 1, message.length(), hy_message_log_file);
        fflush (hy_message_log_file);
#endif
    }
    
    //____________________________________________________________________________________
    void HandleOrStoreApplicationError (_String* error_string, _String const & message) {
        if (error_string) {
            *error_string = message;
        } else {
            HandleApplicationError (message);
        }
    }

    //____________________________________________________________________________________
    void HandleErrorWhileParsing (_String const & error_string, _String const & context) {
        HandleApplicationError (_String ("While parsing:\n") & context & "\n" & error_string);
    }

  //____________________________________________________________________________________
    const _String  PrepareErrorContext     (_String const & context, long from, unsigned long size) {
      
      _StringBuffer result;
      
      if (from > 0) {
          result << "...";
      }
      
      result << _String (context, from, from + size);
      
      if (from + size > context.length()) {
        result << "...";
      }
      
      
      return result;
    }
  
    //____________________________________________________________________________________
    void HandleApplicationError (const _String & message, bool force_exit, bool dump_core) {

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
            fwrite (message.get_str(), 1, message.length(), hy_error_log_file);
            fflush (hy_error_log_file);
        }
            
        if (hy_error_log_file) {
            fprintf (hy_error_log_file, "\n%s", (const char*)message);
        }
            
        _String errMsg;
    #ifdef __HYPHYMPI__
        errMsg = _String("Received an error state from MPI node ") & (long)hy_mpi_node_rank & '\n' & message;
        
        if (hy_mpi_node_rank > 0) {
            errMsg = ConstructAnErrorMessage (message);
            fprintf (stderr, "HYPHYMPI terminated.\n%s\n\n", errMsg.get_str());
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
        if (dump_core) {
            raise (SIGTERM);
        } else {
            exit(1);
        }
#endif
    }
  
  
  //____________________________________________________________________________________
  const _String GetVersionString (void) {
    _StringBuffer theMessage = _String("HYPHY ")& kHyPhyVersion;
    #ifdef __MP__
      theMessage << "(MP)";
    #endif
    #ifdef __HYPHYMPI__
      theMessage <<  "(MPI)";
    #endif
     theMessage << " for ";
    #ifdef __UNIX__
      #if !defined __HEADLESS_WIN32__ && ! defined __MINGW32__
      struct      utsname      name;
      uname       (&name);
      theMessage << name.sysname << " on " << name.machine;
      #endif
      #if defined __MINGW32__
      theMessage <<  "MinGW ";// " & __MINGW32_VERSION;
      #endif
    #endif
    return theMessage;
  }
  
  //____________________________________________________________________________________
  const _String GetTimeStamp (bool do_gmt) {
    time_t c_time;
    time (&c_time);
    
    if (do_gmt) {
      tm  gmt;
      gmtime_r (&c_time, &gmt);
      
      return _StringBuffer (_String((long)1900+gmt.tm_year)) << '/' << _String (1+(long)gmt.tm_mon) << '/'
            << _String ((long)gmt.tm_mday) << ' ' << _String ((long)gmt.tm_hour) << ':' & _String ((long)gmt.tm_min);
    }
    
    tm    local_time;
    localtime_r (&c_time, &local_time);
    char  time_buffer [128];
    
    asctime_r (&local_time, time_buffer);
    return time_buffer;
  }
  
  //____________________________________________________________________________________
  bool    ProcessFileName (_String & path_name, bool isWrite, bool acceptStringVars, hyPointer theP, bool assume_platform_specific, _ExecutionList * caller, bool relative_to_base, bool relative_path_passthrough) {
      
    static  const _String kRelPathPrefix ("../");
    _String errMsg;
          
    try {
      if (path_name == kPromptForFilePlaceholder || path_name == kEmptyString) {
        // prompt user for file
        if (path_name == kEmptyString) {
#if not defined __MINGW32__ && not defined __WINDOZE__
          char tmpFileName[] = "/tmp/HYPHY-XXXXXX";
          int fileDescriptor = mkstemp(tmpFileName);
          if (fileDescriptor == -1){
            throw _String("Failed to create a temporary file name");
          }
          path_name = tmpFileName;
          CheckReceptacleAndStore(&hy_env::last_file_path,kEmptyString,false, new _FString (path_name, false), false);
          close (fileDescriptor);
          return true;
#else
          throw (tempFString & " is not implemented for this platform");
#endif
        } else {
          if (!isWrite) {
            path_name = ReturnFileDialogInput(&hy_base_directory);
          } else {
            path_name = WriteFileDialogInput (&hy_base_directory);
          }
        }
        ProcessFileName(path_name, false,false,theP,false,caller,true);
        CheckReceptacleAndStore(&hy_env::last_file_path,kEmptyString,false, new _FString (path_name, false), false);
        return true;
      }
      
      if (acceptStringVars) {
        path_name = ProcessLiteralArgument (&path_name,(_VariableContainer*)theP, caller);
        if (caller && caller->IsErrorState()) {
          return false;
        }
        
      } else {
        path_name.StripQuotes();
      }
      
      if (path_name.empty()) {
        CheckReceptacleAndStore(&hy_env::last_file_path,kEmptyString,false, new _FString (path_name, false), false);
        return true;
      }
    }
    
    catch (_String const & errmsg) {
      if (caller) {
        caller->ReportAnExecutionError(errMsg);
      } else {
        HandleApplicationError (errMsg);
      }
      return false;
    }
    
    
#if defined __UNIX__ && !defined __MINGW32__
    //UNIX LINES HERE
    if (path_name.Find('\\') != kNotFound) { // DOS (ASSUME RELATIVE) PATH
      path_name = path_name.Replace ("\\","/",true);
    } else if (path_name.Find(':') != kNotFound) { // Mac (Assume Relative) PATH
      path_name = path_name.Replace ("::",":../", true);
      if (path_name.get_char(0)==':') {
        path_name.Trim(1,-1);
      }
      path_name = path_name.Replace (':','/',true);
    }
    
    if (path_name.get_char(0) != '/') { // relative path
      if (pathNames.nonempty() && !relative_path_passthrough) {
          
        _String*    lastPath = relative_to_base ? & hy_base_directory : (_String*)pathNames(pathNames.lLength-1);

        long        f = (long)lastPath->length ()-2L,
                    k = 0L;
        
        // check the last stored absolute path and reprocess this relative path into an absolute.
        while (path_name.BeginsWith (kRelPathPrefix)) {
          if ( (f = lastPath->FindBackwards('/',0,f)-1) == kNotFound) {
            return true;
          }
          path_name.Trim(3,-1);
          k++;
        }

        if (k==0L) {
          path_name = *lastPath & path_name;
        } else {
          path_name = lastPath->Cut(0,f+1) & path_name;
        }
      }
    }
#endif
    
  #if defined __MINGW32__ // WIN/DOS code'
    
    if (path_name.Find('/')!=-1) { // UNIX PATH
      if (path_name.get_char(0)=='/') {
        path_name.Trim(1,-1);
      }
      path_name = path_name.Replace ("/","\\",true);
    } else {
      if (path_name.Find('\\')==kNotFound) {
        // check to see if this is a relative path
        path_name = path_name.Replace ("::",":..\\", true);
        if ((path_name.get_char (0) ==':')) {
          path_name.Trim(1,-1);
        }
        path_name = path_name.Replace (':','\\',true);
      }
    }
    
    if (path_name.Find(':') == kNotFound && path_name.Find("\\\\",0,1) == kNotFound) { // relative path
      
      if (pathNames.nonempty () && !relative_path_passthrough) {
        _String*    lastPath = relative_to_base ? & hy_base_directory : (_String*)pathNames(pathNames.lLength-1);
        long f = (long)lastPath->length() - 2L, k = 0L;
        // check the last stored absolute path and reprocess this relative path into an absolute.
        while (path_name.BeginsWith("..\\")) {
          f = lastPath->FindBackwards('\\',0,f)-1;
          if (f == kNotFound ) {
            return false;
          }
          path_name.Trim(3,-1);
          k++;
        }
        if (k==0) {
          if (lastPath->get_char (lastPath->length()-1L)!='\\') {
            path_name = *lastPath&'\\'& path_name;
          } else {
            path_name = *lastPath& path_name;
          }
        } else {
          path_name = lastPath->Cut(0,f+1)& path_name;
        }
      }
      
    }
    
    
    _StringBuffer escapedString (path_name.length());
    escapedString.SanitizeAndAppend(path_name);
    path_name = escapedString;
    
  #endif
    
    CheckReceptacleAndStore(&hy_env::last_file_path,kEmptyString,false, new _FString (path_name, false), false);
    return true;
  }

} // namespace close
