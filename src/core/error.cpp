/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (spond@ucsd.edu)
 Art FY Poon    (apoon@cfenet.ubc.ca)
 Steven Weaver (sweaver@ucsd.edu)
 
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

#include <stdlib.h>
#include <stdio.h>
#include "errorfns.h"
#include "hy_strings.h"



#if !defined __UNIX__ || defined __HEADLESS__ || defined __MAC__ || defined __HYPHY_GTK__
#include "preferences.h"
#endif

#ifdef __MAC__
#include <Dialogs.h>
#include "HYUtils.h"
#include "HYConsoleWindow.h"
#include "HYDialogs.h"
#endif

#ifdef __WINDOZE__
void WinErrorBox(_String&, bool);
#include "HYDialogs.h"
#endif

#ifdef __HYPHY_GTK__
#include <gtk/gtk.h>
#include "HYConsoleWindow.h"
#include "HYDialogs.h"
#endif

#ifdef   __UNIX__
#if !defined __MINGW32__
#include <sys/utsname.h>
#endif

#ifndef __HYPHY_GTK__
extern  bool dropIntoDebugMode;
#endif
#else
void    SaveConsole (void);
#endif

#include "batchlan.h"

//_____________________________________________________________

_String  DecodeError                    (long);
_String* ConstructAnErrorMessage        (_String&);


//_____________________________________________________________

int      gError                         = _HYNOERROR;

bool     isFixable                      = true,
         skipWarningMessages            = false;

_String  errorReportFormatExpression            ("ERROR_REPORT_FORMAT_EXPRESSION"),
         errorReportFormatExpressionStr            ("_ERROR_TEXT_"),
         errorReportFormatExpressionStack      ("_ERROR_CALL_STACK_"),
         errorReportFormatExpressionStdin      ("_ERROR_STDIN_");


//_____________________________________________________________

bool gStatus(void)
{
    return gError == _HYNOERROR;
}

//_____________________________________________________________

_String DecodeError (long errCode)
{
    switch (errCode) {
    case -101:
        return "Incompatible Operands";
        break;
    case -102:
        return "Operation Undefined for Type";
        break;
    case -103:
        return "Incompatible Matrix Dimensions";
        break;
    case -104:
        return "Bad Matrix Definition";
        break;
    case -105:
        return "Matrix Index Out of Range";
        break;
    case -106:
        return "Bad Matrix Index";
        break;
    case -108:
        return "Memory Full";
        break;
    case -109:
        return "Syntax Error";
        break;
    case -110:
        return "Runtime Expression Error";
        break;
    case -111:
        return "Non-polynomial expression encountered in polynomial calculation";
        break;
    case -171:
        return "Dataset index reference out of range";
        break;
    case -200:
        return "Export Matrix Called With a Non-polynomial Matrix Argument";
        break;
    case -666:
        return "Attempting to operate on an undefined value; this is probably a result of an earlier 'soft' error condition";
        break;
    default:
        return "Unclassified Error";
    }
}

//_____________________________________________________________
bool isError (long errCode)
{
//  if (errCode == 0)
//      acknError (false);
    if (errCode!=_HYNOERROR)  {
        if (errCode<0) {
            gError = errCode;
            isFixable = TRUE;
            return FALSE;
        }
    }
    gError = _HYNOERROR;
    return TRUE;
}

//_____________________________________________________________

void    warnError (long errCode)
{
    if (errCode == -108) {
        warnError (DecodeError (errCode)&_String(" Exiting..."));
    } else {
        WarnError (DecodeError (errCode)&_String(" Exiting..."));
    }
}

//_____________________________________________________________

void    flagError (long errCode)
{
    warnError (DecodeError (errCode));
}


//_____________________________________________________________

void    warnError (const char* theError)
{
    FlagError (theError);
}

//_____________________________________________________________

void    acknError (const char* theError)
{
    WarnError (theError);
}

//_____________________________________________________________

void*   checkPointer (void* p)
{
    if (p) {
        return p;
    }

    warnError(-108);
    return nil;
}

//_______________________________________________________________________
void    ReportWarning (_String st)
{
    checkParameter          (MessageLogging, messageLogFlag, 1.0);

#ifdef  __HEADLESS__
    if (globalInterfaceInstance && messageLogFlag >= 0.1) {
        globalInterfaceInstance->PushWarning (&st);
    }
#else
    if ( !globalMessageFile || messageLogFlag<.1 ) {
        return;
    }

    char   str[] = "\n";
    fwrite (str, 1, 1, globalMessageFile);
    fwrite (st.getStr(), 1, st.Length(), globalMessageFile);
    fflush (globalMessageFile);
#endif
}


//_______________________________________________________________________
void    FlagError (_String st)
{
#ifdef  __HEADLESS__
    if (globalInterfaceInstance) {
        globalInterfaceInstance->PushError (&st);
    }

    terminateExecution = true;
#else
    char  str[] = "\nError:";


    if (globalErrorFile) {
        fwrite (str, 1, 7, globalErrorFile);
        fwrite (st.getStr(), 1, st.Length(), globalErrorFile);
        fflush(globalErrorFile);
    }


#ifdef __HYPHYMPI__
    int     rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

#if !defined __MAC__ && !defined __WINDOZE__
    _String errMsg;
#ifdef __HYPHYMPI__
    errMsg = _String("Received an error state from MPI node ") & (long)rank & '\n' & st;

    if (rank > 0) {
        errMsg = ConstructAnErrorMessage (st);
        MPISendString (errMsg,0,true);
    } else {
        errMsg = _String ("\nMaster node received an error:") & st ;
    }
#else
    errMsg = st;
#endif

    errMsg = ConstructAnErrorMessage (errMsg);

#ifdef  _MINGW32_MEGA_
    SetStatusLine  (errMsg);
#else
    _SimpleList color (255,2,0,0);
    StringToConsole(errMsg, &color);
#endif
#endif

#ifdef __HYPHYMPI__
    if (rank==0) {
        MPI_Abort (MPI_COMM_WORLD,1);
    }
#endif


#if defined __UNIX__ && !defined __HYPHYQT_ &&! defined __HYPHY_GTK__
    if (dropIntoDebugMode)
        while (ExpressionCalculator()) ;
#endif

#ifdef _HY_ABORT_ON_ERROR
    abort ();
#else
    PurgeAll(true);
    exit(1);
#endif
#endif
}

//_______________________________________________________________________
void    WarnErrorWhileParsing (_String st, _String& context)
{
    WarnError (_String ("While parsing:\n") & context & "\n" & st);
}


//_______________________________________________________________________
void WarnError (_String st)
{
    if (currentExecutionList && currentExecutionList->errorHandlingMode == HY_BL_ERROR_HANDLING_SOFT) {
        currentExecutionList->ReportAnExecutionError(st, true);
        return;
    }

#ifdef  __HEADLESS__
    if (globalInterfaceInstance) {
        globalInterfaceInstance->PushError (&st);
    }
    terminateExecution = true;
#else
    char  str[] = "\nError:";

    if (globalErrorFile) {
        fwrite (str, 1, 7, globalErrorFile);
        fwrite (st.getStr(), 1, st.Length(), globalErrorFile);
        fflush (globalErrorFile);
    }

#ifdef __HYPHYMPI__
    int     rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    if (globalMessageFile) {
        fprintf (globalMessageFile, "\n%s", st.sData);
    }

    _String errMsg;
    #ifdef __HYPHYMPI__
        errMsg = _String("Received an error state from MPI node ") & (long)rank & '\n' & st;

        if (rank > 0) {
            errMsg = ConstructAnErrorMessage (st);
            fprintf (stderr, "HYPHYMPI terminated.\n%s\n\n", errMsg.sData);
            MPI_Abort (MPI_COMM_WORLD,1);
            abort   ();
        } else {
            errMsg = _String ("\nMaster node received an error:") & st;
        }
    #else
        errMsg = st;
    #endif



    errMsg = ConstructAnErrorMessage (errMsg);
    

    #ifdef  _MINGW32_MEGA_
        SetStatusLine  (errMsg);
    #else
        _SimpleList color (255,2,0,0);
        StringToConsole(errMsg, &color);
#ifdef __HYPHYQT__
    return;
#endif
#endif


#if defined __UNIX__ && !defined __HYPHY_GTK__ && !defined __HYPHY_QT__
    if (dropIntoDebugMode)
        while (ExpressionCalculator()) ;
#endif
#ifdef __HYPHYMPI__
    if (rank==0) {
        MPI_Abort (MPI_COMM_WORLD,1);
    }
#endif
    //GlobalShutdown();

#ifdef _HY_ABORT_ON_ERROR
    abort ();
#else
    #if defined __MAC__ or defined __WINDOZE__ or defined __HYPHY_GTK__
      terminateExecution = true;
      ProblemReport (errMsg);
    #else
      exit(1);
    #endif
  
#endif

#endif
}

//____________________________________________________________________________________

_String* ConstructAnErrorMessage         (_String& theMessage)
{
    _String* errMsg = new _String (128L,true);

    _List    calls,
             stdins;

    ReturnCurrentCallStack (calls, stdins);

    _FString * errorFormattingExpression = (_FString*)FetchObjectFromVariableByType (&errorReportFormatExpression, STRING);
    bool     doDefault = true;

    if (errorFormattingExpression) {
        _Formula expression;
        _String  expr (*errorFormattingExpression->theString),
                 errMsgLocal;
        _FormulaParsingContext fpc (&errMsgLocal, nil);
        
        if (Parse    (&expression, expr, fpc, nil) == HY_FORMULA_EXPRESSION) {
            CheckReceptacleAndStore(&errorReportFormatExpressionStr, empty, false, new _FString (theMessage, false), false);
            CheckReceptacleAndStore(&errorReportFormatExpressionStack, empty, false, new _Matrix (calls), false);
            CheckReceptacleAndStore(&errorReportFormatExpressionStdin, empty, false, new _Matrix (stdins, false), false);
            _PMathObj expr = expression.Compute();
            if (!terminateExecution && expr && expr->ObjectClass() == STRING) {
                (*errMsg) << ((_FString*)expr)->theString;
                doDefault = false;
            }
        }
    }

    if (doDefault) {
        (*errMsg) << "Error:\n";
        (*errMsg) << theMessage;
        if (calls.lLength) {
            (*errMsg) << "\n\nFunction call stack\n";
            for (unsigned long k = 0; k < calls.lLength; k++) {
                (*errMsg) << (_String((long)k+1) & " : " & (*(_String*)calls(k)) & '\n');
                _String* redir = (_String*)stdins (k);
                if (redir->sLength) {
                    (*errMsg) << "\tStandard input redirect:\n\t\t";
                    (*errMsg) << redir->Replace ("\n","\n\t\t",true);
                }
                (*errMsg) << "-------\n";

            }
        }
    }

    errMsg->Finalize();
    return errMsg;
}

//____________________________________________________________________________________

void ReturnCurrentCallStack      (_List& calls, _List& stdins)
{
    calls.Clear  ();
    stdins.Clear ();

    if (executionStack.lLength) {
        for (long callLevel = executionStack.lLength -1 ; callLevel >=0 ; callLevel --) {
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

//_____________________________________________________________
//EOF





