/*

HyPhy - Hypothesis Testing Using Phylogenies.

This file implements shared function used by
'mains'

Written by SL Kosakovsky Pond
June 11, 2007

Copyright (C) 1997-2007
Primary Development:
  Sergei L Kosakovsky Pond (sergeilkp@mac.com)
Significant contributions from:
  Spencer V Muse (muse@stat.ncsu.edu)
  Simon DW Frost (sdfrost@ucsd.edu)
  Art FY Poon    (apoon@biomail.ucsd.edu)

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

#include "HYSharedMain.h"
#include "HYUtils.h"
#include "HYConsoleWindow.h"
#include "HYDialogs.h"

bool             isSuspended        = false,
                 hasTemplates      = false,
                 highLevelQuit        = false,
                 isRerunAvailable     = false,
                 updateTimer        = false,
                 addToRecent      = true,
                 echoPaused        = false,
                 calculatorMode     = false;

_ExecutionList   ex;

_List            availableTemplateFiles,
                 availablePostProcessors;

_SimpleList      windowPtrs,
                 windowObjects,
                 treeIDReferences,
                 windowObjectRefs;

_String*         argFileName = nil;

//_________________________________________________________________________

void    PrepareToExecuteBatchFile  (void)
{
#ifdef __MAC__
    DisableMenuItem (GetMenuHandle(129),0);
    EnableMenuItem  (GetMenuHandle(131),0);
    EnableMenuItem  (GetMenuHandle(131),1);
    EnableMenuItem  (GetMenuHandle(131),2);
    DisableMenuItem (GetMenuHandle(131),4);
    DisableMenuItem (GetMenuHandle(131),6);
    DisableMenuItem (GetMenuHandle(131),7);
    DisableMenuItem (GetMenuHandle(131),8);
    InvalMenuBar();
#endif

#ifdef __WINDOZE__
    HMENU       hMenu = GetMenu ((HWND)hyphyConsoleWindow->GetOSWindowData()),
                sMenu;
    if (hMenu) {
        EnableMenuItem(hMenu, 0 ,MF_GRAYED|MF_BYPOSITION);
    }
    sMenu = GetSubMenu (hMenu,2);
    if (sMenu) {
        EnableMenuItem(sMenu, 0 ,MF_ENABLED|MF_BYPOSITION);
        EnableMenuItem(sMenu, 1 ,MF_ENABLED|MF_BYPOSITION);
        EnableMenuItem(sMenu, 3 ,MF_GRAYED|MF_BYPOSITION);
        EnableMenuItem(sMenu, 5 ,MF_GRAYED|MF_BYPOSITION);
        EnableMenuItem(sMenu, 6 ,MF_GRAYED|MF_BYPOSITION);
        EnableMenuItem(sMenu, 7 ,MF_GRAYED|MF_BYPOSITION);
    }

#endif

#ifdef __HYPHY_GTK__
    ToggleAnalysisMenu (true);
#endif
}

//_________________________________________________________________________

void    DoneWithExecutionOfBatchFile (bool doPost)
{
#ifdef __MAC__
    if (doPost) {
        EnableMenuItem (GetMenuHandle(131),7);
        EnableMenuItem (GetMenuHandle(201),0);
        for (long i=0; i<availablePostProcessors.lLength; i++) {
            _String* condition = (_String*)(*(_List*)availablePostProcessors(i))(2);
            EnableMenuItem (GetMenuHandle (201),i+1);
            if (condition->sLength) {
                _Formula condCheck (*condition,nil);
                _PMathObj condCheckRes = condCheck.Compute();
                if ((!condCheckRes)||(condCheckRes->Value()<.5)) {
                    DisableMenuItem (GetMenuHandle (201),i+1);
                }
            }
        }
    }

    EnableMenuItem (GetMenuHandle(129),0);
    EnableMenuItem (GetMenuHandle(150),5);
    EnableMenuItem (GetMenuHandle(131),0);
    DisableMenuItem (GetMenuHandle(131),1);
    DisableMenuItem (GetMenuHandle(131),2);
    EnableMenuItem (GetMenuHandle(131),4);
    if (hasTemplates) {
        EnableMenuItem (GetMenuHandle(131),6);
    }
    EnableMenuItem (GetMenuHandle(131),8);
    InvalMenuBar();
#endif

#ifdef __WINDOZE__
    HMENU       hMenu = GetMenu ((HWND)hyphyConsoleWindow->GetOSWindowData()),
                sMenu;

    if (hMenu) {
        EnableMenuItem(hMenu, 0, MF_ENABLED|MF_BYPOSITION);
    }

    sMenu = GetSubMenu (hMenu,2);
    if (sMenu) {
        EnableMenuItem(sMenu, 0 ,MF_GRAYED|MF_BYPOSITION);
        EnableMenuItem(sMenu, 1 ,MF_GRAYED|MF_BYPOSITION);
        EnableMenuItem(sMenu, 3 ,MF_ENABLED|MF_BYPOSITION);
        EnableMenuItem(sMenu, 5 ,MF_ENABLED|MF_BYPOSITION);
        EnableMenuItem(sMenu, 6 ,MF_ENABLED|MF_BYPOSITION);
        EnableMenuItem(sMenu, 7 ,MF_ENABLED|MF_BYPOSITION);
    }
    if (doPost) {
        HMENU       sMenu = GetSubMenu(GetSubMenu (GetMenu((HWND)hyphyConsoleWindow->GetOSWindowData()),2),6);
        for (long i=0; i<availablePostProcessors.lLength; i++) {
            EnableMenuItem(sMenu, i ,MF_ENABLED|MF_BYPOSITION);
            _String* condition = (_String*)(*(_List*)availablePostProcessors(i))(2);
            if (condition->sLength) {
                _Formula condCheck (*condition,nil);
                _PMathObj condCheckRes = condCheck.Compute();
                if ((!condCheckRes)||(condCheckRes->Value()<.5)) {
                    EnableMenuItem(sMenu, i ,MF_GRAYED|MF_BYPOSITION);
                }
            }
        }
    }
    DrawMenuBar((HWND)hyphyConsoleWindow->GetOSWindowData());
    SendMessage ((HWND)hyphyConsoleWindow->GetOSWindowData(), WM_PAINT, NULL, NULL);

#endif

#ifdef __HYPHY_GTK__
    ToggleAnalysisMenu (false);
#endif

}

//_________________________________________________________________________

bool    ExecuteBatchFile (void)
{
    PrepareToExecuteBatchFile ();
    terminateExecution      = false;
    skipWarningMessages     = false;

    _String justTheName     (*argFileName,argFileName->FindBackwards (GetPlatformDirectoryChar(),0,-1)+1,-1);

    AddStringToRecentMenu   (justTheName, *argFileName);
    SetStatusLine           (justTheName,"Loading","00:00:00", -1);

    StartBarTimer           ();
    ApplyPreferences        ();

    ex.Clear                ();
    ReadBatchFile           (*argFileName,ex);
    ex.Execute();

    StopBarTimer            ();
    setParameter            (VerbosityLevelString, 0.0, nil);
    //SetStatusLine             (justTheName, "Finished", empty, -1, HY_SL_FILE|HY_SL_TASK);
    BufferToConsole         ("\n");

    ReportAnalysisAsFinished (empty);

    isRerunAvailable        = true;
    terminateExecution      = false;

    ex.ResetFormulae();
    DoneWithExecutionOfBatchFile ();
    return true;
}

//_________________________________________________________________________
void    RunTemplate (long idx)
{
    PurgeAll        (windowPtrs.lLength==0);
    _String         pathName = libDirectory&"TemplateBatchFiles"&GetPlatformDirectoryChar();
    pathNames&&     &pathName;
    pathName    = pathName&*(_String*)(*(_List*)availableTemplateFiles(idx))(2);

    if (!argFileName) {
        argFileName = new _String(pathName);
    } else {
        *argFileName = pathName;
    }

    ExecuteBatchFile();
}

//_________________________________________________________________________
bool    OpenBatchFile (bool openOrNot, _String* dL)
{
    PurgeAll            (windowPtrs.lLength==1);
    if (openOrNot)
        if (!PopUpFileDialog(" Please select a batch file to run:",dL)) {
            return false;
        }

    if (!argFileName || argFileName->sLength == 0) {
        return false;
    }

    _String       pathName (*argFileName);
    PushFilePath  (pathName);
    return true;
}

//_________________________________________________________________________
void    ExecuteAPostProcessor (_String justTheName)
{
    PrepareToExecuteBatchFile ();
    _ExecutionList postEx;

    _String        postFile    = libDirectory & "TemplateBatchFiles" & GetPlatformDirectoryChar() & justTheName,
                   pathName (postFile);

#ifdef __MAC__
    volumeName  = baseDirectory.Cut (0, baseDirectory.Find(':'));
#endif

    SetStatusLine ( justTheName, "Loading", "00:00:00", -1);

    PushFilePath      (pathName);
    ReadBatchFile     (postFile, postEx);
    StartBarTimer     ();
    terminateExecution = false;
    postEx.Execute    ();
    SetStatusLine     (justTheName, "Finished", empty, -1, HY_SL_FILE|HY_SL_TASK);
    PopFilePath       ();
    terminateExecution = false;
    StopBarTimer      ();
    DoneWithExecutionOfBatchFile ();
}

//_________________________________________________________________________
void    updateTimerF (_String& rec, long time_diff)
{
    rec.FormatTimeString (time_diff);
}

//__________________________________________________________________________________
void    ReadInTemplateFiles(void)
{
    _String     fileIndex = libDirectory&"TemplateBatchFiles"&GetPlatformDirectoryChar()&"files.lst";
    FILE      * modelList = doFileOpen (fileIndex.getStr(), "r");

    if (!modelList) {
        return;
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


    // try reading post processing files

    fileIndex = libDirectory&"TemplateBatchFiles"&GetPlatformDirectoryChar()&"postprocessors.lst";
    if (! (modelList = doFileOpen (fileIndex.getStr(),"r"))) {
        return;
    }

    _String postData (modelList);
    fclose  (modelList);

    if (postData.sLength) {
        _ElementaryCommand::ExtractConditions(postData,0,availablePostProcessors);
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

            availablePostProcessors.Replace(i,&thisFile,true);
        }

        for (long counter=0; counter<availablePostProcessors.countitems(); counter++) {
            _String * postItem = (_String*)(*(_List*)availablePostProcessors(counter))(0);
            if (postItem->Equal(&menuSeparator)) {
                continue;
            }
            postItem = (_String*)(*(_List*)availablePostProcessors(counter))(1);
            _String tryFileName = libDirectory & "TemplateBatchFiles"& GetPlatformDirectoryChar() & (*postItem);
            FILE    * tryFile = doFileOpen (tryFileName.getStr(), "r");
            if (tryFile) {
                fclose(tryFile);
            } else {
                availablePostProcessors.Delete(counter);
                counter--;
            }
        }
    }
}

//____________________________________________________________________________________________
long  SelectATemplate (void)
{
    _SimpleList std,
                vc (availableTemplateFiles.lLength,0,1),
                selection;

    std<<2;
    std<<1;

    return HandleHierListSelection (availableTemplateFiles, std, vc, "Select a standard analysis to run",selection,1);
}

//____________________________________________________________________________________________
_String     MatrixExpCounter             (void)
{
    _Matrix B(61,61,true,true);
    long i;
    for (i=0; i<740; i++) {
        B[((_Parameter)genrand_int32()/RAND_MAX_32)*(61*61-1)] = ((_Parameter)genrand_int32())/RAND_MAX_32;
    }

    clock_t startTime = clock(),
            otherTime;

    for (i=0; i<100000; i++) {
        B.Exponentiate();
        //C *= B;
        otherTime = clock()-startTime;
        if (otherTime>CLOCKS_PER_SEC && i>=5) {
            break;
        }
    }
    if (B.GetHDim()) {
        _Parameter result = otherTime/(double)CLOCKS_PER_SEC;
        result = i*1.2/result;
        return _String("61x61 Random Sparse Matrix Exps/sec : ")&_String(result);
    }
    return empty;
}

//____________________________________________________________________________________________
void SpoolFile (void)
{
    FILE *f = doFileOpen (argFileName->getStr(), "r");
    if (f) {
        _String lateralus (f);
        fclose (f);
        _String fS = _String ("\n-------------------File ") & *argFileName & " ------------------------\n";
        StringToConsole (fS);
        StringToConsole (lateralus);
    }
}

//________________________________________________________

void        RunStandardAnalyses (void)
{
    long menuChoice=SelectATemplate();
    if (menuChoice >= 0) {
        RunTemplate(menuChoice);
    }
}

//_________________________________________________________________________

void    ReportAnalysisAsFinished  (_String text, bool fgTheConsole)
{
#if defined __MAC__ || defined __WINDOZE__ || defined __HYPHY_GTK__
    SetStatusLine (empty,text.sLength?text:_String("Finished"),empty,-HY_SL_DONE,HY_SL_PERCENT|HY_SL_TASK|HY_SL_PERCENT|HY_SL_FORCE|HY_SL_DONE);
    if (fgTheConsole) {
        hyphyConsoleWindow->BringToFront();
    }
#endif
}



//EOF
