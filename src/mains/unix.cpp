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
#include <stdlib.h>
#include "batchlan.h"
#include "calcnode.h"
#include <unistd.h>
#include "polynoml.h"

#if !defined __MINGW32__

    #include <termios.h>
    #include <signal.h>
    #define __HYPHY_HANDLE_TERM_SIGNAL__

#endif


#ifdef _MINGW32_MEGA_
#include <Windows.h>
HANDLE _HY_MEGA_Pipe = INVALID_HANDLE_VALUE;
#endif

#ifdef  __UNITTEST__
#include "gtest/gtest.h"
#include "ut_strings.h"
#endif

#if defined   __MP2__ || defined __MP__
#include <pthread.h>
#endif
#include "likefunc.h"

#ifndef __HYPHY_NO_CURL__
//#define __HYPHYCURL__
#endif

#ifdef  __HYPHYCURL__
#include <curl/curl.h>
#endif

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

#ifdef _OPENMP
#include "omp.h"
#endif

_List   availableTemplateFiles,
        availablePostProcessors,
        loggedUserInputs;

_String baseArgDir,
        libArgDir,
        loggedFileEntry ("__USER_ENTRY__");

void    ReadInTemplateFiles         (void);
long    DisplayListOfChoices        (void);
void    ProcessConfigStr            (_String&);
void    ReadInPostFiles             (void);
long    DisplayListOfPostChoices    (void);
_String getLibraryPath              (void);


extern  long
systemCPUCount;

extern  _String         VerbosityLevelString,
        errorFileName,
        messageFileName,
        baseDirectory,
        libDirectory,
        shortMPIReturn,
        dialogPrompt;

long    mainArgCount = 0;

bool    usePostProcessors = false,
        calculatorMode    = false,
        updateMode       = false,
        pipeMode         = false,
        dropIntoDebugMode = false,
        logInputMode   = false,
        needExtraNL    = false;

char    prefFileName[] = ".hyphyinit";

#ifdef  __HYPHYMPI__
extern int          _hy_mpi_node_rank;

void            mpiNormalLoop    (int, int, _String &);
void            mpiOptimizerLoop (int, int);
#endif


#ifdef     __HYPHY_HANDLE_TERM_SIGNAL__
    volatile sig_atomic_t hyphy_sigterm_in_progress = 0;
     
    void    hyphy_sigterm_handler (int sig)
    {
       if (hyphy_sigterm_in_progress)
           raise (sig);
           
       hyphy_sigterm_in_progress = 1;
     
       WarnError (_String("HyPhy killed by signal ") & (long)sig);
       
       signal (sig, SIG_DFL);
       raise  (sig);
    }
     
#endif  

//bool  terminateExecution = false;

//____________________________________________________________________________________

_String getLibraryPath() {

    char    curWd[4096],
            dirSlash = GetPlatformDirectoryChar();
    getcwd (curWd,4096);

    _String baseDir (curWd);

    if (baseDir.getChar (baseDir.sLength-1) != dirSlash) {
        baseDir=baseDir & dirSlash;
    }

#if defined _HYPHY_LIBDIRECTORY_
    _String libDir (_HYPHY_LIBDIRECTORY_);

    if (libDir.getChar (libDir.sLength-1) != dirSlash) {
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
      if(hyphyPath.sLength != 0) {
        libDir = hyphyPath;
      }
    }

    return libDir;

}

//__________________________________________________________________________________
void    ReadInTemplateFiles(void)
{
    _String fileIndex;
    fileIndex = *((_String*)pathNames(0)) &"TemplateBatchFiles/files.lst";
    FILE* modelList = fopen (fileIndex.getStr(),"r");
    if (!modelList) {
        fileIndex = baseArgDir&"TemplateBatchFiles/files.lst";
        modelList = fopen (fileIndex.getStr(),"r");
        if (!modelList) {
            return;
        }
    } else {
        baseArgDir = *((_String*)pathNames(0));
    }

    _String theData (modelList);
    fclose (modelList);
    
    if (theData.sLength) {
        _ElementaryCommand::ExtractConditions(theData,0,availableTemplateFiles);
        for (long i = 0; i<availableTemplateFiles.countitems(); i++) {
            _String* thisString = (_String*)availableTemplateFiles(i);
            _List   thisFile;
            _ElementaryCommand::ExtractConditions(*thisString,thisString->FirstNonSpaceIndex(),thisFile,',');
            if (thisFile.lLength!=3) {
                availableTemplateFiles.Delete(i);
                i--;
                continue;
            }
            for (long j = 0; j<3; j++) {
                ((_String*)thisFile(j))->StripQuotes();
            }
            availableTemplateFiles.Replace(i,&thisFile,true);
        }

    }
}

//__________________________________________________________________________________
void    ReadInPostFiles(void)
{
    //if (!likeFuncList.lLength)
    //  return;

    _String fileIndex;
    FILE* modelList = fopen (fileIndex.getStr(),"r");
    fileIndex = libArgDir &"TemplateBatchFiles/postprocessors.lst";
    modelList = fopen (fileIndex.getStr(),"r");
    
    if (modelList == NULL) {
        return;
    }

    _String theData (modelList);
    fclose (modelList);
    
    if (theData.sLength) {
        _ElementaryCommand::ExtractConditions(theData,0,availablePostProcessors);
        for (long i = 0; i<availablePostProcessors.countitems(); i++) {
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
            if (*(_String*)thisFile(0)!=_String("SEPARATOR")) {
                fileIndex = *((_String*)pathNames(0)) &"TemplateBatchFiles/" & *(_String*)thisFile(1);
                //printf ("%s\n", fileIndex.sData);
                FILE* dummyFile = fopen (fileIndex,"r");
                if (!dummyFile) {
                    fileIndex =libArgDir&"TemplateBatchFiles/"& *(_String*)thisFile(1);
                    dummyFile = fopen (fileIndex,"r");
                }
                if (dummyFile) {
                    fclose (dummyFile);
                    _String* condition = (_String*)thisFile(2);
                    if (condition->sLength) {
                        _Formula condCheck (*condition,nil);
                        _PMathObj condCheckRes = condCheck.Compute();
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
long    DisplayListOfChoices (void)
{
    ReadInTemplateFiles();

    if (!availableTemplateFiles.lLength) {
        return -1;
    }

    long        choice = -1;
    char        buffer[2048];
    _String     fileAbbr,
                *thisLine;
    _SimpleList categoryDelimiters;
    _List       categoryHeadings;

    for (choice = 0; choice< availableTemplateFiles.lLength; choice++) {
        thisLine = (_String*)(*(_List*)availableTemplateFiles(choice))(2);
        if (thisLine->sData[0]=='!') {
            categoryDelimiters<<choice;
            fileAbbr = *thisLine;
            fileAbbr.Trim (1,-1);
            categoryHeadings && &fileAbbr;
        }
    }

    choice = -1;
    if (categoryDelimiters.lLength==0) {
        while (choice == -1) {
            for (choice = 0; choice<availableTemplateFiles.lLength; choice++) {
                printf ("\n\t(%s):%s",((_String*)(*(_List*)availableTemplateFiles(choice))(0))->getStr(),
                        ((_String*)(*(_List*)availableTemplateFiles(choice))(1))->getStr());
            }
            printf ("\n\n Please type in the abbreviation for the file you want to use (or press ENTER to process custom batch file):");
            fgets (buffer,2048,stdin);
            fgets (buffer,2048,stdin);
            fileAbbr = buffer;
            if (fileAbbr.FirstNonSpaceIndex()<0) {
                return -1;
            }
            fileAbbr.UpCase();
            for (choice = 0; choice<availableTemplateFiles.lLength; choice++) {
                if (fileAbbr.Equal((_String*)(*(_List*)availableTemplateFiles(choice))(0))) {
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
                          verString (GetVersionString().getStr());

                if (verString.sLength<header.sLength-2) {
                    _String padder (128,true);
                    long    poop = (header.sLength-2-verString.sLength)/2;
                    if (!poop) {
                        poop = 1;
                    }
                    for (choice=0; choice<poop; choice++) {
                        padder << ' ';
                    }
                    padder.Finalize();
                    verString = padder & '/' & verString & "\\" & padder;
                }

                printf ("\n\033[2J\033[H%s\n%s\n\n",verString.getStr(), header.getStr());
                for (choice = 0; choice<categoryHeadings.lLength; choice++) {
                    printf ("\n\t(%ld) %s",choice+1,((_String*)categoryHeadings(choice))->getStr());
                }

                printf ("\n\n Please select type of analyses you want to list (or press ENTER to process custom batch file):");


                fgets (buffer,2048,stdin);
                fileAbbr = buffer;

                if (logInputMode) {
                    loggedUserInputs && & fileAbbr;
                }

                if (fileAbbr.FirstNonSpaceIndex()<0) {
                    return -1;
                }

                choice = fileAbbr.toNum();

                if ( choice>0 && choice<=categoryHeadings.lLength) {
                    categNumber = choice-1;
                }
            } else {
                printf ("\n\033[2J\033[H ***************** FILES IN '%s' ***************** \n\n",((_String*)categoryHeadings(categNumber))->getStr());
                long start = categoryDelimiters.lData[categNumber]+1,
                     end = categNumber==categoryDelimiters.lLength-1?availableTemplateFiles.lLength:categoryDelimiters.lData[categNumber+1];

                for (choice = start; choice<end; choice++) {
                    printf ("\n\t(%ld) %s",choice-start+1,((_String*)(*(_List*)availableTemplateFiles(choice))(1))->getStr());
                }

                printf ("\n\n Please select the file you want to use (or press ENTER to return to the list of analysis types):");

                fileAbbr = *StringFromConsole ();

                if (logInputMode) {
                    loggedUserInputs && & fileAbbr;
                }

                if (fileAbbr.FirstNonSpaceIndex()<0) {
                    categNumber = -1;
                } else {
                    choice = fileAbbr.toNum();
                    if ((choice>0 && choice<=end-start)) {
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
long    DisplayListOfPostChoices (void)
{
    long choice = -1;
    if (availablePostProcessors.lLength) {
        _String fileAbbr;
        printf ("\033[2J\033[H\n\t Available Result Processing Tools\n\t ---------------------------------\n\n");
        while (choice == -1) {
            for (choice = 0; choice<availablePostProcessors.lLength; choice++) {
                printf ("\n\t(%ld):%s",choice+1,
                        ((_String*)(*(_List*)availablePostProcessors(choice))(0))->getStr());
            }
            printf ("\n\n Please type in the abbreviation for the tool you want to use (or press q to exit):");
            fileAbbr = *StringFromConsole();
            fileAbbr.UpCase();
            if (logInputMode) {
                loggedUserInputs && & fileAbbr;
            }
            if (!fileAbbr.sLength||((fileAbbr.sLength==1)&&(fileAbbr.sData[0]=='Q'))) {
                return -1;
            }
            choice = fileAbbr.toNum();

            if (choice<=0 || choice>availablePostProcessors.lLength) {
                choice = -1;
            }
        }
    }
    return choice;
}


//__________________________________________________________________________________
void    ProcessConfigStr (_String& conf)
{
    _String errMsg;
    for (long i=1; i<conf.sLength; i++) {
        switch (conf.sData[i]) {
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
            dropIntoDebugMode = true;
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
        //case 'i':
        //case 'I':
        //{
        //pipeMode = true;
        //break;
        //}
        default: {
            errMsg = "Option ";
            errMsg = errMsg & conf.sData[i] & " is not valid and is ignored";
            ReportWarning (errMsg);
        }
        }
    }
}


//__________________________________________________________________________________

void hyphyBreak (int signo)
{
    //terminateExecution = false;
    printf ("\nInterrupt received %d. HYPHY will break into calculator mode at the earliest possibility...\n", signo);
}

//__________________________________________________________________________________
void    SetStatusBarValue           (long,_Parameter,_Parameter)
{

}
//__________________________________________________________________________________
void    SetStatusLine               (_String s)
{
#ifdef  _MINGW32_MEGA_
    if (_HY_MEGA_Pipe != INVALID_HANDLE_VALUE) {
        DWORD bytesWritten = 0;
        if (WriteFile (_HY_MEGA_Pipe,(LPCVOID)s.sData,s.sLength,&bytesWritten,NULL) == FALSE || bytesWritten != s.sLength) {
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
void    SetStatusLineUser   (_String s)
{
    setvbuf(stdout, NULL, _IONBF, 0);
    BufferToConsole("\33[2K\r");
    StringToConsole(s);
    needExtraNL = true;
}

//__________________________________________________________________________________

#ifndef __UNITTEST__
int main (int argc, char* argv[])
{
    mainArgCount = argc - 1;


#ifdef  __HYPHYMPI__
    int            rank,
                   size;

    MPI_Init       (&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    setParameter  (mpiNodeID,    (_Parameter)rank);
    setParameter    (mpiNodeCount, (_Parameter)size);
    _hy_mpi_node_rank = rank;

    if (rank == 0) {
        mpiNodesThatCantSwitch.Populate (size,1,0);
        /* {
              char hostname[256];
              gethostname(hostname, sizeof(hostname));
              printf("PID %d on %s ready for attach\n", getpid(), hostname);
              fflush(stdout);
              //getchar ();
          } */
#endif


        //for (long k=0; k<NSIG; k++)
        //{
        //  signal(k, &hyphyBreak);
        //}

#ifdef  __HYPHYMPI__
    }
#endif

#ifdef __HYPHY_HANDLE_TERM_SIGNAL__
    if (signal (SIGTERM, hyphy_sigterm_handler) == SIG_IGN)
         signal (SIGTERM, SIG_IGN);
     if (signal (SIGINT, hyphy_sigterm_handler) == SIG_IGN)
         signal (SIGINT, SIG_IGN);
         
#endif
    char    curWd[4096],
            dirSlash = GetPlatformDirectoryChar();
    getcwd (curWd,4096);

    _String baseDir (curWd);

    if (baseDir.getChar (baseDir.sLength-1) != dirSlash) {
        baseDir=baseDir & dirSlash;
    }
    
    _String libDir = getLibraryPath();
    pathNames&& &libDir;
    _String argFile;

    libDirectory  = libDir;
    libArgDir     = libDirectory;
    baseDirectory = baseDir;
    baseArgDir    = baseDirectory;

  
#ifdef _OPENMP
    systemCPUCount = omp_get_max_threads();
#endif

#ifdef _MINGW32_MEGA_
    {
        char pid[16];
        snprintf (pid,16,"%u", GetCurrentProcessId());

        _String pipeName = _String("\\\\.\\pipe\\MEGAPipe") & pid;
        printf ("Pipe name = %s\n", pipeName.sData);
        if ((_HY_MEGA_Pipe = CreateFile(pipeName.sData, GENERIC_WRITE, 0, NULL, OPEN_EXISTING, 0, NULL)) == INVALID_HANDLE_VALUE) {
            char* lpMsgBuf;
            FormatMessage(
                FORMAT_MESSAGE_ALLOCATE_BUFFER |
                FORMAT_MESSAGE_FROM_SYSTEM |
                FORMAT_MESSAGE_IGNORE_INSERTS,
                NULL,
                GetLastError(),
                MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                (LPTSTR) &lpMsgBuf,
                0, NULL );
            FlagError (_String("Failed to create a pipe named '") & pipeName & "' to send data from HyPhy to MEGA. Error: "&lpMsgBuf);
        }
    }
#endif

    for (long i=1; i<argc; i++) {
        _String thisArg (argv[i]);
        if (thisArg.sData[0]=='-') {
            ProcessConfigStr (thisArg);
        } else if (thisArg.beginswith ("BASEPATH=")) {
            baseArgDir = thisArg.Cut(9,-1);
            if (baseArgDir.sLength) {
                if (baseArgDir.sData[baseArgDir.sLength-1]!=dirSlash) {
                    baseArgDir = baseArgDir&dirSlash;
                }
                baseDirectory = baseArgDir;
                pathNames.Delete    (0);
                pathNames&&         &baseDirectory;
           }
        } else if (thisArg.beginswith ("LIBPATH=")) {
            libArgDir = thisArg.Cut(8,-1);
            if (libArgDir.sLength) {
                if (libArgDir.sData[libArgDir.sLength-1] != dirSlash) {
                    libArgDir = libArgDir & dirSlash;
                }
                libDirectory = libArgDir;
                pathNames.Delete    (0);
                pathNames&&         &libDirectory;
           }
        } else if (thisArg.beginswith ("USEPATH=")) {
            baseDir             = thisArg.Cut(8,-1);
            errorFileName       = baseDir & errorFileName;
            messageFileName     = baseDir & messageFileName;
            pathNames.Delete    (0);
            pathNames&&         &baseDir;
        } else
#ifdef __MP__
            if (thisArg.beginswith ("CPU=")) {
                _String cpus = thisArg.Cut(4,-1);
                systemCPUCount = cpus.toNum();
                if (systemCPUCount<1) {
                    systemCPUCount = 1;
                }
                pthread_setconcurrency (systemCPUCount+1);
            } else
#endif
                argFile = thisArg;
    }

    GlobalStartup();
  
  /*_String exp2 ("1+x*t+y*t"),
          exp1 ("y*t");
  
  _Formula f1 (exp1),
           f2 (exp2);
  
  _MathObject *pl = ((_Polynomial*)f1.ConstructPolynomial ())->Add (f2.ConstructPolynomial());
  
  fprintf (stdout, "%s\n", _String ((_String*)pl->toStr()).getStr());*/
           

    _ExecutionList ex;
  
  
  if (calculatorMode) {
        printf ("\nHYPHY is running in calculator mode. Type 'exit' when you are finished.\n");
        while (ExpressionCalculator()) ;
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
    FILE     * testPrefFile = fopen (prefFile.sData,"r");
    if (!testPrefFile) {
        prefFile = baseArgDir & prefFileName;
        testPrefFile = fopen (prefFile.sData,"r");
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
        if (!argFile.sLength) {
            long selection = -2;
            if (!updateMode) {
                selection = DisplayListOfChoices();
            }

            if (selection == -1) {
                dialogPrompt = "Batch file to run:";
                _String fStr (ReturnDialogInput (true));
                if (logInputMode) {
                    _String tts = loggedFileEntry&fStr;
                    loggedUserInputs && & tts;
                }

                PushFilePath (fStr);
                ReadBatchFile (fStr,ex);
            } else {
                _String templ;

                if (selection >= 0) {
                    templ = baseArgDir &"TemplateBatchFiles" & dirSlash;
                } else {
                    templ = baseArgDir & "TemplateBatchFiles" & dirSlash & "WebUpdate.bf";
                }

                if (selection >= 0) {
                    templ= templ&*(_String*)(*(_List*)availableTemplateFiles(selection))(2);
                }

                PushFilePath (templ);
                ReadBatchFile (templ,ex);
            }
        } else {
#ifndef __MINGW32__
            if (argFile.sData[0] != '/') {
                argFile       = baseDirectory & argFile;
            }
#else
            if (argFile.sData[1] != ':') { // not an absolute path
                argFile       = baseDirectory & argFile;
            }
#endif
            PushFilePath  (argFile);
            ReadBatchFile (argFile,ex);
        }

        ex.Execute();

        if (usePostProcessors && (!updateMode)) {
            ReadInPostFiles();
            printf ("\n\n**********Continue with result processing (y/n)?");
            _String c_str (StringFromConsole());

            if (logInputMode) {
                loggedUserInputs && & c_str;
            }

            if (c_str.getChar(0) !='n' && c_str.getChar(0)!='N' ) {
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

                    if (c_str.getChar(0)=='n' || c_str.getChar(0)=='N' ) {
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

    PurgeAll                    (true);
    GlobalShutdown              ();

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
