/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2002
Sergei L Kosakovsky Pond (sergeilkp@mac.com)
Spencer V Muse (muse@stat.ncsu.edu)

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

#include <stdio.h>
#include <unistd.h>

#include "global_things.h"
#include "function_templates.h"
#include "trie_iterator.h"

using namespace hy_global;

#include "batchlan.h"
#include "calcnode.h"
#include "polynoml.h"

#if defined __MINGW32__
    #include <shlwapi.h>
#else
    #include <termios.h>
    #include <signal.h>
    #include <unistd.h>
    #include <sys/stat.h>
    #define __HYPHY_HANDLE_TERM_SIGNAL__

#endif


const char hy_usage[] =
"usage: hyphy or HYPHYMPI [-h] [--help]"
"[-c] "
"[-d] "
"[-i] "
"[-p] "
"[BASEPATH=directory path] "
"[CPU=integer] "
"[LIBPATH=library path] "
"[USEPATH=library path] "
"[<standard analysis name> or <path to hyphy batch file>] [--keyword value ...] [positional arguments ...]"
"\n";


const char analysis_help_message [] =
"Available analysis command line options\n"
"---------------------------------------\n"
"Use --option VALUE syntax to invoke\n"
"If a [reqired] option is not provided on the command line, the analysis will prompt for its value\n"
"[conditionally required] options may or not be required based on the values of other options\n\n";

const char hy_help_message [] =
"Execute a HyPhy analysis, either interactively, or in batch mode\n"
"optional flags:\n"
"  -h --help                show this help message and exit\n"
"  -c                       calculator mode; causes HyPhy to drop into an expression evaluation until 'exit' is typed\n"
"  -d                       debug mode; causes HyPhy to drop into an expression evaluation mode upon script error\n"
"  -i                       interactive mode; causes HyPhy to always prompt the user for analysis options, even when defaults are available\n"
"  -p                       postprocessor mode; drops HyPhy into an interactive mode where general post-processing scripts can be selected\n"
"                           upon analysis completion\n\n"
"  -m                       write diagnostic messages to messages.log\n"
"optional global arguments:\n"
"  BASEPATH=directory path  defines the base directory for all path operations (default is pwd)\n"
"  CPU=integer              if compiled with OpenMP multithreading support, requests this many threads; HyPhy could use fewer than this\n"
"                           but never more; default is the number of CPU cores (as computed by OpenMP) on the system\n"
"  LIBPATH=directory path   defines the directory where HyPhy library files are located (default installed location is /usr/local/lib/hyphy\n"
"                           or as configured during CMake installation\n"
"  USEPATH=directory path   specifies the optional working and relative path directory (default is BASEPATH)\n\n"
"  ENV=expression           set HBL enviroment variables via explicit statements\n"
"                           for example ENV='DEBUG_MESSAGES=1;WRITE_LOGS=1'\n"
"  batch file to run        if specified, execute this file, otherwise drop into an interactive mode\n"
"  analysis arguments       if batch file is present, all remaining positional arguments are interpreted as inputs to analysis prompts\n\n"
"optional keyword arguments (can appear anywhere); will be consumed by the requested analysis\n"
"  --keyword value          will be passed to the analysis (which uses KeywordArgument directives)\n"
"                           multiple values for the same keywords are treated as an array of values for multiple selectors\n"
"\n"
"usage examples:\n\n"
"Select a standard analysis from the list : \n\thyphy -i \n"
"Run a standard analysis with default options and one required user argument; \n\thyphy busted --alignment path/to/file\n"
"Run a standard analysis with additional keyword arguments \n\thyphy busted --alignment path/to/file --srv No\n"
"See whcih arguments are understood by a standard analysis \n\thyphy busted --help\n"
"Run a custom analysis and pass it some arguments \n\thyphy path/to/hyphy.script argument1 'argument 2' \n"
;

const char hy_available_cli_analyses [] =
"Available standard analyses and their [standard analysis name] are listed below:\n\n"
"        [meme] Test for episodic site-level selection using MEME (Mixed Effects Model of Evolution).\n"
"        [contrast-fel] Use a FEL method to test which sites in a gene may be associated with adaptation to a different environment.\n"
"        [fel] Test for pervasive site-level selection using FEL (Fixed Effects Likelihood).\n"
"        [fubar] Test for pervasive site-level selection using FUBAR (Fast Unconstrained Bayesian AppRoximation for inferring selection).\n"
"        [fade] Test a protein alignment for directional selection towards specific amino acids along a specified set of test branches using FADE (a FUBAR Approach to Directional Evolution).\n"
"        [slac] Test for pervasive site-level selection using SLAC (Single Likelihood Ancestor Counting).\n"
"        [busted] Test for episodic gene-wide selection using BUSTED (Branch-site Unrestricted Statistical Test of Episodic Diversification).\n"
"        [bgm] Apply Bayesian Graphical Model inference to substitution histories at individual sites.\n"
"        [absrel] Test for lineage-specific evolution using the branch-site method aBS-REL (Adaptive Branch-Site Random Effects Likelihood).\n"
"        [relax] Test for relaxation of selection pressure along a specified set of test branches using RELAX (a random effects test of selection relaxation).\n"
"        [gard] Screen an alignment for recombination using GARD (Genetic Algorithm for Recombination Detection).\n\n"
;


#ifdef _MINGW32_MEGA_
  #include <Windows.h>
  HANDLE _HY_MEGA_Pipe = INVALID_HANDLE_VALUE;
#endif


#ifdef  __UNITTEST__
#include "gtest/gtest.h"
#include "ut_strings.h"
#endif

#include "likefunc.h"

#ifndef __HYPHY_NO_CURL__
//#define __HYPHYCURL__
#endif

#ifdef  __HYPHYCURL__
#include <curl/curl.h>
#endif


#ifdef _OPENMP
#include "omp.h"
#endif

//#define _COMPARATIVE_LF_DEBUG_CHECK "/Users/sergei/Desktop/lf.sequence"
//#define _COMPARATIVE_LF_DEBUG_DUMP "/Users/sergei/Desktop/lf.sequence"

#if defined _COMPARATIVE_LF_DEBUG_CHECK
    FILE* _comparative_lf_debug_matrix_content_file = doFileOpen (_COMPARATIVE_LF_DEBUG_CHECK, "r", true);
    _String _comparative_lf_debug_matrix_content (_comparative_lf_debug_matrix_content_file);
    _Matrix* _comparative_lf_debug_matrix = new _Matrix (_comparative_lf_debug_matrix_content, true);
#else
    #if defined _COMPARATIVE_LF_DEBUG_DUMP
    _GrowingVector* _comparative_lf_debug_matrix = new _GrowingVector ();
    #endif
#endif

_List   availableTemplateFiles,
        availablePostProcessors,
        loggedUserInputs;

_Trie   availableTemplateFilesAbbreviations;

_String baseArgDir,
        libArgDir;

const   _String kLoggedFileEntry ("__USER_ENTRY__"),
                kHelpKeyword     ("--help"),
                kVersionKeyword  ("--version"),
                kVerboseKeyword  ("--verbose");

void    ReadInTemplateFiles         (void);
long    DisplayListOfChoices        (void);
void    ProcessConfigStr            (_String const&);
void    ProcessKWStr                (_String const& arg, _String const& arg2, _AssociativeList& kwargs);
void    ReadInPostFiles             (void);
long    DisplayListOfPostChoices    (void);
_String getLibraryPath              (void);





bool    usePostProcessors = false,
        calculatorMode    = false,
        updateMode       = false,
        pipeMode         = false,
        logInputMode   = false,
        displayHelpAndExit = false;
char    prefFileName[] = ".hyphyinit";

#ifdef  __HYPHYMPI__

void            mpiNormalLoop    (int, int, _String &);
void            mpiOptimizerLoop (int, int);
#endif


#ifdef     __HYPHY_HANDLE_TERM_SIGNAL__
    volatile sig_atomic_t hyphy_sigterm_in_progress = 0;
     
    void    hyphy_sigterm_handler (int sig) {
       if (hyphy_sigterm_in_progress)
           raise (sig);
           
        hyphy_sigterm_in_progress = 1;
        if (lockedLFID != -1) {
            ((_LikelihoodFunction*)likeFuncList(lockedLFID))->_TerminateAndDump(_String("HyPhy killed by signal ") & (long)sig);
        } else {
            HandleApplicationError (_String("HyPhy killed by signal ") & (long)sig);

        }
    
       
       signal (sig, SIG_DFL);
       raise  (sig);
    }
     
#endif  

//bool  terminate_execution = false;

//____________________________________________________________________________________


#define _HYPHY_MAX_PATH_LENGTH 8192L

_String getLibraryPath() {
  char    dirSlash = get_platform_directory_char();

#ifdef __MINGW32__
  TCHAR buffer[_HYPHY_MAX_PATH_LENGTH];
  GetModuleFileName(NULL, buffer, _HYPHY_MAX_PATH_LENGTH);
  _String baseDir (buffer);
  baseDir.Trim (0, baseDir.FindBackwards (dirSlash, 0, -1) - 1L);
#else
    char    curWd[_HYPHY_MAX_PATH_LENGTH];
    getcwd (curWd,_HYPHY_MAX_PATH_LENGTH);

    _String baseDir (curWd);
#endif

    if (baseDir.get_char (baseDir.length()-1) != dirSlash) {
        baseDir=baseDir & dirSlash;
    }

#if defined _HYPHY_LIBDIRECTORY_
    _String libDir (_HYPHY_LIBDIRECTORY_);

    if (libDir.get_char (libDir.length()-1) != dirSlash) {
        libDir=libDir & dirSlash;
    }

#else
     pathNames&& &baseDir;
    _String libDir = baseDir;
#endif

    // SW20141119: Check environment libpath and override default path if it exists
    // TODO: Move string to globals in v3
    // TODO: Move function to helpers location in v3
    char* hyphyEnv = getenv("HYPHY_PATH");

    if(hyphyEnv) {
        _String hyphyPath(hyphyEnv);
        if(hyphyPath.nonempty()) {
            libDir = hyphyPath;
        }
    } else {
        _String tryLocal = baseDir & "res" & dirSlash;
        
        struct stat sb;
        
        if (stat((const char*)tryLocal, &sb) == 0 && S_ISDIR(sb.st_mode)) {
            libDir = tryLocal;
        }
        
    }

    return libDir;

}

//__________________________________________________________________________________

void   _helper_clear_screen (void) {
  #ifdef __MINGW32__
    system("cls");
  #else
    printf ("\033[2J\033[H");
  #endif
}

//__________________________________________________________________________________
void    ReadInTemplateFiles(void) {
    _String dir_sep (get_platform_directory_char()),
            fileIndex = *((_String*)pathNames(0)) & hy_standard_library_directory & dir_sep & "files.lst";
  
    FILE* modelList = fopen (fileIndex.get_str(),"r");
    if (!modelList) {
        fileIndex = baseArgDir& hy_standard_library_directory & dir_sep & "files.lst";
        modelList = fopen (fileIndex.get_str(),"r");
        if (!modelList) {
            return;
        }
    } else {
        baseArgDir = *((_String*)pathNames(0));
    }

    _String theData (modelList);
    fclose (modelList);
    
    if (theData.length()) {
        _List extracted_files;
        _ElementaryCommand::ExtractConditions(theData,0,extracted_files);
        extracted_files.ForEach([] (BaseRef item, unsigned long) -> void {
            _String* thisString = (_String*)item;
            _List  *  thisFile  = new _List();
            _ElementaryCommand::ExtractConditions(*thisString,thisString->FirstNonSpaceIndex(),*thisFile,',');
            if (thisFile->countitems() == 3L) {
                thisFile->ForEach ([] (BaseRef item, unsigned long) -> void { ((_String*)item)->StripQuotes();});
                availableTemplateFiles.AppendNewInstance(thisFile);
                bool did_insert = false;
                if (((_String*)thisFile->GetItem(0))->nonempty()) {
                    _String to_insert (((_String*)thisFile->GetItem(0))->ChangeCase(kStringLowerCase));
                    availableTemplateFilesAbbreviations.InsertExtended(to_insert, availableTemplateFiles.countitems()-1, false, &did_insert);
                    if (!did_insert) {
                        HandleApplicationError(_String("Duplicate analysis keyword (not case sensitive) in ") & _String ((_String*)thisFile->toStr()));
                        return;
                    }
                }
            } else {
              DeleteObject (thisFile);
            }
          
        });
      
    }
}

//__________________________________________________________________________________
void    ReadInPostFiles(void) {
    static _String kSeparator ("SEPARATOR");
    //if (!likeFuncList.lLength)
    //  return;
 
    _String dir_sep (get_platform_directory_char());
  
    _String fileIndex = libArgDir & hy_standard_library_directory & dir_sep & "postprocessors.lst";
    FILE*  modelList = fopen (fileIndex.get_str(),"r");
    
    if (!modelList) {
        return;
    }

    _String theData (modelList);
    fclose (modelList);
    
    if (theData.length()) {
        _ElementaryCommand::ExtractConditions(theData,0,availablePostProcessors);
      
        for (unsigned long i = 0; i<availablePostProcessors.countitems(); i++) {
            _String* thisString = (_String*)availablePostProcessors(i);
            _List   thisFile;
            _ElementaryCommand::ExtractConditions(*thisString,thisString->FirstNonSpaceIndex(),thisFile,',');
            if (thisFile.lLength!=3) {
                availablePostProcessors.Delete(i);
                i--;
                continue;
            }
            for (long j = 0; j<3; j++) {
                ((_String*)thisFile(j))->StripQuotes();
            }
            if ( *(_String*)thisFile(0) != kSeparator) {
                fileIndex = *((_String*)pathNames(0)) & hy_standard_library_directory & dir_sep & *(_String*)thisFile(1);
                FILE* dummyFile = fopen (fileIndex,"r");
                if (!dummyFile) {
                    fileIndex =libArgDir& hy_standard_library_directory & dir_sep & *(_String*)thisFile(1);
                    dummyFile = fopen (fileIndex,"r");
                }
                if (dummyFile) {
                    fclose (dummyFile);
                    _String* condition = (_String*)thisFile(2);
                    if (condition->nonempty()) {
                        _Formula condCheck (*condition,nil);
                        HBLObjectRef condCheckRes = condCheck.Compute();
                        if ((!condCheckRes)||(condCheckRes->Value()<.5)) {
                            availablePostProcessors.Delete(i);
                            i--;
                            continue;
                        }
                    }
                    *(_String*)thisFile(1) = fileIndex;
                    availablePostProcessors.Replace(i,&thisFile,true);
                    continue;
                }
            }
            availablePostProcessors.Delete(i);
            i--;
        }

    }
}

//__________________________________________________________________________________
long    DisplayListOfChoices (void) {
    
  
    if (!availableTemplateFiles.lLength) {
        return -1;
    }

    long        choice = -1L;
  
    _SimpleList categoryDelimiters;
    _List       categoryHeadings;

    for (choice = 0; choice< availableTemplateFiles.countitems(); choice++) {
      _String const * this_line = (_String const *)availableTemplateFiles.GetItem (choice, 2);
      
         if ( this_line->char_at (0) == '!') {
            categoryDelimiters<<choice;
            _String * category_heading = new _String (*this_line);
            category_heading->Trim (1,kStringEnd);
            categoryHeadings < category_heading;
        }
    }

    choice = -1;
    if (categoryDelimiters.lLength==0) {
        while (choice == -1) {
            for (choice = 0; choice<availableTemplateFiles.lLength; choice++) {
                printf ("\n\t(%s):%s", ((_String const *)availableTemplateFiles.GetItem (choice, 0))->get_str() , ((_String const *)availableTemplateFiles.GetItem (choice, 1)) -> get_str());
            }
            printf ("\n\n Please type in the abbreviation for the file you want to use (or press ENTER to process custom batch file):");
          
            _String user_input = _String (StringFromConsole()).ChangeCase(kStringUpperCase);
          
            if (user_input.FirstNonSpaceIndex() == kNotFound) {
                return -1;
            }
          
            for (choice = 0; choice<availableTemplateFiles.lLength; choice++) {
                if (user_input == *(_String const *)availableTemplateFiles.GetItem (choice, 0)) {
                    break;
                }
            }
          
            if (choice==availableTemplateFiles.lLength) {
                choice=-1;
            }
        }
    } else {
        long categNumber = -1;
        while (choice==-1) {
            if (categNumber<0) {
                _String   header ("***************** TYPES OF STANDARD ANALYSES *****************"),
                          verString (GetVersionString().get_str());

                if ( verString.length() + 2 < header.length()) {
                    _String    padder (_String (' '), MAX (1, (header.length()-2-verString.length())/2));
                    verString = padder & '/' & verString & "\\" & padder;
                }
                _helper_clear_screen ();
                printf ("%s\n%s\n\n",verString.get_str(), header.get_str());
                for (choice = 0; choice<categoryHeadings.lLength; choice++) {
                    printf ("\n\t(%ld) %s",choice+1,((_String*)categoryHeadings(choice))->get_str());
                }

                printf ("\n\n Please select type of analyses you want to list (or press ENTER to process custom batch file):");


                _String user_input = StringFromConsole();

                if (logInputMode) {
                    loggedUserInputs && & user_input;
                }

                if (user_input.FirstNonSpaceIndex() == kNotFound) {
                    return -1;
                }

                choice = user_input.to_long();

                if ( choice>0 && choice<=categoryHeadings.lLength) {
                    categNumber = choice-1;
                }
            } else {
                _helper_clear_screen ();
                printf ("***************** FILES IN '%s' ***************** \n\n",((_String*)categoryHeadings(categNumber))->get_str());
                long start = categoryDelimiters.list_data[categNumber]+1,
                     end = categNumber==categoryDelimiters.lLength-1?availableTemplateFiles.lLength:categoryDelimiters.list_data[categNumber+1];

                for (choice = start; choice<end; choice++) {
                    printf ("\n\t(%ld) %s",choice-start+1,((_String const *)availableTemplateFiles.GetItem (choice, 1))->get_str());
                }

                printf ("\n\n Please select the analysis you would like to perform (or press ENTER to return to the list of analysis types):");

                _String user_input = StringFromConsole();
              
                if (logInputMode) {
                    loggedUserInputs && & user_input;
                }

                if (user_input.FirstNonSpaceIndex() == kNotFound) {
                    categNumber = -1;
                } else {
                    choice = user_input.to_long();
                    if ( choice>0 && choice<=end-start) {
                        return start+choice-1;
                    }
                }

            }
            choice = -1;
        }
    }
    return choice;
}

//__________________________________________________________________________________
long    DisplayListOfPostChoices (void) {
    long choice = -1;
  
    if (availablePostProcessors.lLength) {
        _helper_clear_screen ();
        printf ("\n\t Available Result Processing Tools\n\t ---------------------------------\n\n");
        while (choice == -1) {
            for (choice = 0; choice<availablePostProcessors.lLength; choice++) {
                printf ("\n\t(%ld):%s",choice+1,
                        ((_String*)(*(_List*)availablePostProcessors(choice))(0))->get_str());
            }
            printf ("\n\n Please type in the abbreviation for the tool you want to use (or press q to exit):");
            _String user_input = _String(StringFromConsole()).ChangeCase (kStringUpperCase);
          
            if (logInputMode) {
                loggedUserInputs && & user_input;
            }
          
            if (user_input.empty () || user_input == 'Q') {
                return -1;
            }
            choice = user_input.to_long();

            if (choice<=0 || choice>availablePostProcessors.lLength) {
                choice = -1;
            }
        }
    }
    return choice;
}


void    DisplayHelpMessage (void) {
   printf ("%s\n%s", hy_usage, hy_help_message);
   printf ("Available standard keyword analyses (located in %s)\n", getLibraryPath().get_str());
   TrieIterator options (&availableTemplateFilesAbbreviations);
   for (TrieIteratorKeyValue ti : options) {
       printf ("\t%s \t%s\n", ti.get_key().get_str(), ((_String*)(availableTemplateFiles.GetItem(availableTemplateFilesAbbreviations.GetValue(ti.get_value()),1)))->get_str());
   }
   printf ("\n");
   //fprintf( stderr, "%s\n%s\n%s", hy_usage, hy_help_message, hy_available_cli_analyses );
   exit (0);
}
            
void    ProcessConfigStr (_String const & conf) {
    for (unsigned long i=1UL; i<conf.length(); i++) {
        switch (char c = conf.char_at (i)) {
            case 'h':
            case 'H': {
                displayHelpAndExit = true;
            }
            break;

            case 'p':
          case 'P': {
              usePostProcessors = true;
              break;
          }
          case 'c':
          case 'C': {
              calculatorMode = true;
              break;
          }
          case 'd':
          case 'D': {
              hy_drop_into_debug_mode = true;
              break;
          }
                
         case 'i' :
         case 'I' : {
             ignore_kw_defaults = true;
             break;
         }
                
          case 'u':
          case 'U': {
              updateMode = true;
              break;
          }
          case 'l':
          case 'L': {
              logInputMode = true;
              break;
          }

          case 'm':
          case 'M': {
            hy_messages_log_name = "messages.log";
            break;
          }

        default: {
              ReportWarning (_String ("Option" ) & _String (c).Enquote() & " is not valid command line option and will be ignored");
          }
        }
    }
}

//__________________________________________________________________________________

void    ProcessKWStr (_String const & conf, _String const & conf2, _AssociativeList & kwargs) {
    if (conf.length() == 2) {
        HandleApplicationError ("Cannot have an empty (--) keyword as a command line argument");
    } else {
        // check to see if the kw already exists and append subsequent ones
        _String key = conf.Cut (2, kStringEnd);
        HBLObjectRef existing_value = kwargs.GetByKey(key);
        if (existing_value) {
            if (existing_value->ObjectClass() == ASSOCIATIVE_LIST) {
                (*(_AssociativeList*)existing_value) < _associative_list_key_value {nil, new _FString (conf2)};
            } else {
                _AssociativeList * replacement_list = new _AssociativeList;
                (*replacement_list) << _associative_list_key_value {nil, existing_value} < _associative_list_key_value {nil, new _FString (conf2)};
                kwargs.MStore(key, replacement_list, false);
            }
        } else {
            kwargs.MStore(key, new _FString (conf2), false);
        }
    }
}


//__________________________________________________________________________________

void hyphyBreak (int signo) {
    //terminate_execution = false;
    printf ("\nInterrupt received %d. HYPHY will break into calculator mode at the earliest possibility...\n", signo);
}

//__________________________________________________________________________________
void    SetStatusBarValue           (long,hyFloat,hyFloat){

}
//__________________________________________________________________________________
void    SetStatusLine               (_String s) {
#ifdef  _MINGW32_MEGA_
    if (_HY_MEGA_Pipe != INVALID_HANDLE_VALUE) {
        DWORD bytesWritten = 0;
        if (WriteFile (_HY_MEGA_Pipe,(LPCVOID)s.sData,s.length(),&bytesWritten,NULL) == FALSE || bytesWritten != s.length()) {
            _String errMsg ("Failed to write the entire status update to a named MEGA pipe");
            StringToConsole (errMsg);
        }
        FlushFileBuffers(_HY_MEGA_Pipe);
    } else {
        StringToConsole (s);
    }
#endif

}

//__________________________________________________________________________________
void    SetStatusLineUser   (_String const s) {
  if ( has_terminal_stderr ) { // only print to terminal devices
    setvbuf(stderr, NULL, _IONBF, 0);
    BufferToConsole("\33[2K\r", stderr);
    StringToConsole(s, stderr);
    hy_need_extra_nl = true;
  }
}

//__________________________________________________________________________________

#ifndef __UNITTEST__
int main (int argc, char* argv[]) {
    
 
#ifdef _COMPARATIVE_LF_DEBUG_DUMP
    FILE * comparative_lf_debug_matrix_content_file = doFileOpen (_COMPARATIVE_LF_DEBUG_DUMP, "w");
#endif

#ifdef  __HYPHYMPI__
    MPI_Init       (&argc, &argv);
    int            rank,
    size;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    hy_mpi_node_rank  = rank;
    hy_mpi_node_count = size;
    
    
    /*if (rank > 0) {
        volatile int i = 0;
        char hostname[256];
        gethostname(hostname, sizeof(hostname));
        printf("PID %d on %s ready for attach\n", getpid(), hostname);
        fflush(stdout);
        while (0 == i)
          sleep(5);
        printf("PID %d continuing\n", getpid());
    }*/
  
#endif

#ifdef __HYPHY_HANDLE_TERM_SIGNAL__
    if (signal (SIGTERM, hyphy_sigterm_handler) == SIG_IGN)
         signal (SIGTERM, SIG_IGN);
     if (signal (SIGINT, hyphy_sigterm_handler) == SIG_IGN)
         signal (SIGINT, SIG_IGN);
#endif
    //printf ("%e\n", mapParameterToInverval (1.322753E-23, 0x2, false));
    //exit (0);

    /*long read = 0L;
    hyFloat value = 0.0;
    
    printf ("%ld\n", sscanf (" 0.1e2 beavis", "%lf%n", &value, &read));
    */
    
    char    curWd[4096],
            dirSlash = get_platform_directory_char();
    
    getcwd (curWd,4096);
  
    _String baseDir (curWd);
  

    if (baseDir.get_char (baseDir.length()-1) != dirSlash) {
        baseDir=baseDir & dirSlash;
    }
  
    _String libDir = getLibraryPath();
    

#ifdef __MINGW32__
     baseDir = libDir;
#endif
  
    pathNames&& &libDir;
    _String argFile;
  

    hy_lib_directory  = libDir;
    libArgDir     = hy_lib_directory;
    hy_base_directory = baseDir;
    baseArgDir    = hy_base_directory;
    
    
#ifdef _OPENMP
    system_CPU_count = omp_get_max_threads();
#endif

    _List positional_arguments;
    _AssociativeList kwargs;
  
    const _String path_consts [] = {"BASEPATH=", "LIBPATH=", "USEPATH=", "CPU=", "ENV="};

    for (unsigned long i=1UL; i<argc; i++) {
      _String thisArg (argv[i]);

      if (thisArg.get_char(0)=='-') { // -[LETTER] arguments
          if (thisArg.get_char (1) == '-') {
              if (thisArg == kHelpKeyword) {
                  displayHelpAndExit = true;
                  continue;
              }
              if (thisArg == kVersionKeyword) {
                  StringToConsole(GetVersionString()); NLToConsole();
                  exit (0);
              }
              if (thisArg == kVerboseKeyword) {
                  force_verbosity_from_cli = true;
                  continue;
              }
              if (i + 1 < argc) {
                  _String payload (argv[++i]);
                  ProcessKWStr (thisArg, payload, kwargs);
              } else {
                  HandleApplicationError ("A keyword (--) command line argument must be followed by a value");
              }
          } else {
              ProcessConfigStr (thisArg);
          }
      } else if (thisArg.BeginsWith (path_consts[0])) { // BASEPATH
          baseArgDir = thisArg.Cut(path_consts[0].length(),kStringEnd);
          if (baseArgDir.length()) {
              if (baseArgDir (-1L) != dirSlash) {
                  baseArgDir = baseArgDir & dirSlash;
              }
              hy_base_directory = baseArgDir;
              pathNames.Delete    (0);
              pathNames&&         &hy_base_directory;
         }
      } else if (thisArg.BeginsWith (path_consts[1])) { // LIBPATH
          libArgDir = thisArg.Cut(path_consts[1].length(),kStringEnd);
          if (libArgDir.length()) {
              if (libArgDir (-1L) != dirSlash) {
                  libArgDir = libArgDir & dirSlash;
              }
              hy_lib_directory = libArgDir;
              pathNames.Delete    (0);
              pathNames&&         &hy_lib_directory;
         }
      } else if (thisArg.BeginsWith (path_consts[2])) {
          baseDir                      = thisArg.Cut(path_consts[2].length(),kStringEnd);
          hy_error_log_name            = baseDir & hy_error_log_name;
          hy_messages_log_name         = baseDir & hy_messages_log_name;
          pathNames.Delete    (0);
          pathNames&&         &baseDir;
      } else if (thisArg.BeginsWith (path_consts[3])) {
          system_CPU_count  = Maximum (1L, thisArg.Cut(path_consts[3].length(),kStringEnd).to_long());
      } else if (thisArg.BeginsWith (path_consts[4])) {
          hy_env::cli_env_settings = thisArg.Cut(path_consts[4].length(),kStringEnd);
      }
      else
      //argFile = thisArg;
      positional_arguments && &thisArg;
    }
    
    GlobalStartup();
    ReadInTemplateFiles();
    
    if (positional_arguments.empty () && displayHelpAndExit) {
        // --help returns the message below
        DisplayHelpMessage();
        GlobalShutdown();
        return 0;
    }
    
    //ObjectToConsole(&availableTemplateFilesAbbreviations);
    
    if (positional_arguments.Count()) {
        argFile = *(_String*) positional_arguments.GetItem(0UL);
        _String locase (argFile.ChangeCase(kStringLowerCase));
        long is_shortcut = availableTemplateFilesAbbreviations.FindKey(locase);
        
        if (is_shortcut != kNotFound) {
            argFile = baseArgDir & hy_standard_library_directory & dirSlash & *(_String*) availableTemplateFiles.GetItem(availableTemplateFilesAbbreviations.GetValue(is_shortcut), 2);
        }
        // check to see if this is an analysis key
    }
  
    _ExecutionList ex;
    
#ifdef __HYPHYMPI__
    if (hy_mpi_node_rank == 0L) 
#endif
    ex.SetKWArgs (&kwargs);
    
    
    if (calculatorMode) {
        if (argFile.length()) {
          PushFilePath  (argFile);
          ReadBatchFile (argFile,ex);
          ex.Execute();
        }
        printf ("\nHYPHY is running in calculator mode. Type 'exit' when you are finished.\n");
        while (ExpressionCalculator()) ;
        GlobalShutdown();
        return 0;
    }

    if (pipeMode) {
        _String bfIn (stdin);
        _ExecutionList exIn (bfIn);
        exIn.Execute();
        GlobalShutdown();
        return 0;
    }

    
    // try to read the preferences
    _String     prefFile (curWd);
    prefFile = prefFile & '/' & prefFileName;
    FILE     * testPrefFile = fopen (prefFile.get_str(),"r");
    if (!testPrefFile) {
        prefFile = baseArgDir & prefFileName;
        testPrefFile = fopen (prefFile.get_str(),"r");
    }
    if (testPrefFile) {
        fclose(testPrefFile);
        ReadBatchFile (prefFile,ex);
        ex.Execute();
        ex.Clear();
    }
    //printf ("Node %d before mpiParallelOptimizer\n", rank);
#ifdef __HYPHYMPI__
    if (rank>0) {
        //if (mpiParallelOptimizer || mpiPartitionOptimizer)
        //  mpiOptimizerLoop (rank, size);
        //else
        _String defaultBaseDirectory = *(_String*)pathNames(0);
        mpiNormalLoop (rank, size, defaultBaseDirectory);
        /*argFile = "SHUTDOWN_CONFIRM";
        MPISendString (argFile, senderID);*/
    } else {
#endif
        if (!argFile.length()) {
            long selection = -2;
            if (!updateMode) {
                selection = DisplayListOfChoices();
            }

            if (selection == -1) {
                dialogPrompt = "Batch file to run:";
                _String fStr (ReturnDialogInput (true));
                if (logInputMode) {
                    _String tts = kLoggedFileEntry&fStr;
                    loggedUserInputs && & tts;
                }

                PushFilePath (fStr);
                ReadBatchFile (fStr,ex);
            } else {
                _String templ;

                if (selection >= 0) {
                    templ = baseArgDir & hy_standard_library_directory & dirSlash;
                } else {
                    templ = baseArgDir & hy_standard_library_directory & dirSlash & "WebUpdate.bf";
                }

                if (selection >= 0) {
                    templ= templ&*(_String*)(*(_List*)availableTemplateFiles(selection))(2);
                }

                PushFilePath (templ);
                ReadBatchFile (templ,ex);
            }
        } else {
            
            if (positional_arguments.Count () > 1) {
                ex.stdinRedirectAux = new _List;
                ex.stdinRedirect = new _AVLListXL (ex.stdinRedirectAux);
                for (unsigned long i = 1UL; i < positional_arguments.countitems (); i++) {
                    char buf[256];
                    snprintf(buf, 255, "%012ld", i);
                    ex.stdinRedirect->Insert (new _String(buf), (long)positional_arguments.GetItem (i), true);
                }
            }

            
#ifndef __MINGW32__
            if (argFile.char_at (0) != '/') {
                argFile       = hy_base_directory & argFile;
            }
#else
            if (argFile.char_at (1) != ':') { // not an absolute path
                argFile       = hy_base_directory & argFile;
            }
#endif
            PushFilePath  (argFile);
            
            // if this is a nexus file, it will be executed here
            ReadBatchFile (argFile,ex);
        }

        if (displayHelpAndExit) {
#ifdef __HYPHYMPI__
            if (hy_mpi_node_rank == 0L) {
#endif
            NLToConsole();
            BufferToConsole(analysis_help_message);
            StringToConsole(ex.GenerateHelpMessage());
            NLToConsole();
#ifdef __HYPHYMPI__
            }
#endif
            PurgeAll                    (true);
            GlobalShutdown              ();
            return 0;
        }


         ex.Execute();

 
        if (usePostProcessors && (!updateMode)) {
            ReadInPostFiles();
            printf ("\n\n**********Continue with result processing (y/n)?");
            _String c_str (StringFromConsole());

            if (logInputMode) {
                loggedUserInputs && & c_str;
            }

            if (c_str.get_char(0) !='n' && c_str.get_char(0)!='N' ) {
                long choice = DisplayListOfPostChoices();
                while (choice != -1) {
                    _ExecutionList postEx;
                    argFile = *(_String*)(*(_List*)availablePostProcessors(choice-1))(1);
                    PushFilePath (argFile);
                    ReadBatchFile (argFile, postEx);
                    postEx.Execute();
                    PopFilePath ();
                    printf ("\n\n**********Continue with result processing (y/n)?");

                    c_str = StringFromConsole();
                    if (logInputMode) {
                        loggedUserInputs && & c_str;
                    }

                    if (c_str.get_char(0)=='n' || c_str.get_char(0)=='N' ) {
                        break;
                    }

                    choice = DisplayListOfPostChoices();
                }
            }
        }
#ifdef __HYPHYMPI__
    }
    ReportWarning               (_String ("Node ") & (long)rank & " is shutting down\n");
#endif


#ifdef _MINGW32_MEGA_
    if (_HY_MEGA_Pipe != INVALID_HANDLE_VALUE) {
        CloseHandle (_HY_MEGA_Pipe);
    }
#endif
  
#if defined _COMPARATIVE_LF_DEBUG_CHECK
    fclose (_comparative_lf_debug_matrix_content_file);
#endif
#if defined _COMPARATIVE_LF_DEBUG_DUMP
    _comparative_lf_debug_matrix->Trim();
    _comparative_lf_debug_matrix->toFileStr(comparative_lf_debug_matrix_content_file);
    fclose (comparative_lf_debug_matrix_content_file);
#endif


    PurgeAll                    (true);
    ex.ClearExecutionList();
    

    GlobalShutdown              ();
    


#ifdef __MINGW32__
  fflush (stdout);
  system("PAUSE");
#endif
  
#ifdef __HYPHYMPI__
    if (rank == 0) {
        printf ("\n\n");
    }
#endif

}

#endif

#ifdef  __UNITTEST__
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
#endif
