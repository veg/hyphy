/*

HyPhy - Hypothesis Testing Using Phylogenies.

This file implements shared function used by
'mains'

Written by SL Kosakovsky Pond
June 11, 2007

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

#include "HYSharedMain.h"

#ifndef __HYPHYQT__
    #include "HYUtils.h"
    #include "HYConsoleWindow.h"
    #include "HYDialogs.h"
#else
    #include "hyphymain.h"
    #include "hyphy_qt_helpers.h"
    _String   menuSeparator       ("SEPARATOR");
#endif

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
#ifdef __HYPHY_GTK__
    ToggleAnalysisMenu (true);
#endif
}

//_________________________________________________________________________

void    DoneWithExecutionOfBatchFile (bool doPost)
{
#ifdef __HYPHY_GTK__
    ToggleAnalysisMenu (false);
#endif

#ifdef __HYPHYQT__
    ClearStatusLine();
    _hyPrimaryConsoleWindow->DisplayPrompt();
#endif

}

//_________________________________________________________________________

bool    ExecuteBatchFile (void)
{
    PrepareToExecuteBatchFile ();
    terminateExecution      = false;
    skipWarningMessages     = false;

    _String justTheName     (*argFileName,argFileName->FindBackwards (GetPlatformDirectoryChar(),0,-1)+1,-1);

    _hyPrimaryConsoleWindow->AddStringToRecentMenu   (justTheName, *argFileName);
    _hyPrimaryConsoleWindow->SetStatusLine           (justTheName,"Loading","00:00:00", -1);
    _hyPrimaryConsoleWindow->StartBarTimer           ();
    
    _String info = _String("Running '") & *argFileName & "'\n";
    StringToConsole(info);

    
    #ifndef __HYPHYQT__
    ApplyPreferences        ();
    #endif 
    
    ex.Clear                ();
    ReadBatchFile           (*argFileName,ex);
    ex.Execute();

    _hyPrimaryConsoleWindow->StopBarTimer            ();
    
    setParameter            (VerbosityLevelString, 0.0, nil);
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
    
    #if !defined __HYPHYQT__
    
    if (openOrNot)
        if (!PopUpFileDialog(" Please select a batch file to run:",dL)) {
            return false;
        }
    
    #else
        if (!argFileName) {
            argFileName = new _String (_hyQTFileDialog("Select an HBL file to run",empty,false));
        } else {
            *argFileName = _hyQTFileDialog("Select an HBL file to run",empty,false);
        }            
    #endif 

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

    _hyPrimaryConsoleWindow->SetStatusLine ( justTheName, "Loading", "00:00:00", -1);

    PushFilePath      (pathName);
    ReadBatchFile     (postFile, postEx);
    _hyPrimaryConsoleWindow->StartBarTimer     ();
    terminateExecution = false;
    postEx.Execute    ();
    _hyPrimaryConsoleWindow->SetStatusLine     (justTheName, "Finished", empty, -1, HY_SL_FILE|HY_SL_TASK);
    PopFilePath       ();
    terminateExecution = false;
    _hyPrimaryConsoleWindow->StopBarTimer      ();
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

    return -1;
    //return HandleHierListSelection (availableTemplateFiles, std, vc, "Select a standard analysis to run",selection,1);
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
