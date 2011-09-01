/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2009
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon              (apoon@cfenet.ubc.ca)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#include <stdlib.h>
#include <stdio.h>
#include "errorfns.h"
#include "hy_strings.h"



#if !defined __UNIX__ || defined __HEADLESS__
#include "preferences.h"
#endif

#ifdef __MAC__
#include <Dialogs.h>
#include "HYUtils.h"
#endif

#ifdef __WINDOZE__
void WinErrorBox(_String&, bool);
#endif

#ifdef __HYPHY_GTK__
#include <gtk/gtk.h>
#include "HYConsoleWindow.h"
#endif

#ifdef   __UNIX__
#if !defined __MINGW32__
#include <sys/utsname.h>
#endif
extern  bool dropIntoDebugMode;
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
    StringToConsole(errMsg);
#endif
#endif

#ifdef __MAC__
    Str255            err;
    StringToStr255   (st,err);
    ParamText        (err,NULL,NULL,NULL);
    Alert            (128, (ModalFilterUPP)NULL);
    WritePreferences ();
    SaveConsole      ();
#endif

#ifdef __WINDOZE__
    if (st.sLength>255) {
        st = st.Cut(0,255);
    }
    WritePreferences();
    WinErrorBox(st,false);
#endif

#ifdef __HYPHYMPI__
    if (rank==0) {
        MPI_Abort (MPI_COMM_WORLD,1);
    }
#endif

#ifdef __UNIX__
    if (dropIntoDebugMode)
        while (ExpressionCalculator()) ;
#endif
    //GlobalShutdown();
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
void    WarnError (_String st)
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
        fflush (globalErrorFile);
    }

#ifdef __HYPHYMPI__
    int     rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (globalMessageFile) {
        fprintf (globalMessageFile, "\n%s", st.sData);
    }
#if !defined __MAC__ && !defined __WINDOZE__
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
    StringToConsole(errMsg);
#endif
#endif

#ifdef __MAC__
    if (!skipWarningMessages) {
        Str255 err;
        err[0] = st.sLength>255?255:st.sLength;
        memcpy (err+1,st.getStr(),st.sLength>255?255:st.sLength);
        ParamText (err,NULL,NULL,NULL);
        char alertCode;
#ifndef __OLDMAC__
#ifdef TARGET_API_MAC_CARBON
        alertCode = Alert (129, (ModalFilterUPP)NULL);
#else
        alertCode = Alert (129, (RoutineDescriptor*)NULL);
#endif
#else
        alertCode = Alert (129, NULL);
#endif
        terminateExecution = true;
        if (alertCode == 2) {
            skipWarningMessages = true;
        } else if (alertCode == 3) {
            WritePreferences();
            SaveConsole();
            //GlobalShutdown();
            // graceless exit; no need to clean stuff up
            exit(1);
        }
    }
    return;
#endif
#ifdef __HYPHY_GTK__
    if (!skipWarningMessages) {
        GtkWidget *dialog = gtk_message_dialog_new (
                                hyphyConsoleWindow?GTK_WINDOW(gtk_widget_get_ancestor(hyphyConsoleWindow->theWindow,GTK_TYPE_WINDOW)):NULL,
                                GTK_DIALOG_MODAL,
                                GTK_MESSAGE_WARNING,
                                GTK_BUTTONS_NONE,
                                "The following error occurred:\n %s",
                                st.sData);

        gtk_dialog_add_button (GTK_DIALOG(dialog),"Skip Further Messages",2);
        gtk_dialog_add_button (GTK_DIALOG(dialog),"Quit",3);
        gtk_dialog_add_button (GTK_DIALOG(dialog),"OK",1);
        char alertCode = gtk_dialog_run (GTK_DIALOG (dialog));
        gtk_widget_destroy (dialog);

        terminateExecution = true;
        if (alertCode == 2) {
            skipWarningMessages = true;
        } else if (alertCode == 3) {
            WritePreferences();
            SaveConsole     ();
            GlobalShutdown  ();
            exit            (1);
        }
    }
    return;
#endif
#ifdef __WINDOZE__
    if (!skipWarningMessages) {
        if (st.sLength>255) {
            st = st.Cut(0,255);
        }
        WinErrorBox     (st, true);
        terminateExecution = true;
    }
    return;
#endif
#ifdef __UNIX__
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
    PurgeAll(true);
    exit(1);
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
        _String  expr (*errorFormattingExpression->theString);
        long     varRef = -1;
        if (Parse    (&expression, expr, varRef, nil, nil, false) == HY_FORMULA_EXPRESSION) {
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
            for (long k = 0; k < calls.lLength; k++) {
                (*errMsg) << (_String(k+1) & " : " & (*(_String*)calls(k)) & '\n');
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





