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

#include "HYConsoleWindow.h"
#include "HYTextBox.h"
#include "HYButtonBar.h"
#include "HYLabel.h"
#include "HYUtils.h"
#include "HYEventTypes.h"
#include "HYCanvas.h"
#include "HYDialogs.h"
#include "HYUtils.h"
#include "batchlan.h"

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

#define  HY_CONSOLE_OUT_ROW  0
#define  HY_CONSOLE_IN_ROW   1

_HYColor consoleIBarColor = {168,173,176};

_String  startEchoing  ("Begin saving to file"),
         stopEchoing   ("Stop saving to file"),
         pauseEchoing  ("Pause saving to file"),
         resumeEchoing ("Resume saving to file"),
         cState        ("File"),
         cTask         ("Status"),
         cInput        ("Waiting on user input"),
         suspendedPref ("(Suspended)"),
         newLN         ("\n"),
         userHookinDir ("UserAddins");



extern   _String            dialogPrompt,
         baseDirectory;

_HYConsoleWindow*  hyphyConsoleWindow = nil,
                * _hyPrimaryConsoleWindow = nil;

_List    userHookins;

#ifdef  __WINDOZE__
extern bool hyphyExiting;
#endif

#ifdef __HYPHYMPI__
extern int  _hy_mpi_node_rank;
#endif


//__________________________________________________________

_String WriteFileDialogInput    (void);
void    LoadUserHookins         (void);
void    FlushConsoleBuffer      (void);


//__________________________________________________________
_HYConsoleWindow::_HYConsoleWindow (_String name):_HYTWindow (name, 2)
{

    if (!userHookins.lLength) {
        LoadUserHookins ();
    }

    _HYRect         canvasSettings = {10,10,0xffff,0xffff,HY_COMPONENT_NO_SCROLL};

    _HYTextBox*     ow      = new _HYTextBox (canvasSettings, GetOSWindowData(), true);

    canvasSettings.top      = 90;
    canvasSettings.bottom   = 90;

    _HYTextBox*     iw      = new _HYTextBox (canvasSettings, GetOSWindowData());

    canvasSettings.left     = 25;
    canvasSettings.right    = 25;

    _HYButtonBar*   bb      = new _HYButtonBar (canvasSettings,GetOSWindowData());

    bb->SetMessageRecipient (this);
    iw->SetMessageRecipient (this);
    ow->SetMessageRecipient (this);
    iw->boxFlags |= HY_TB_WRAP;

    AddObject (ow,false);      // 0
    AddObject (iw,false);      // 1
    AddObject (bb,false);      // 2

    SetTableDimensions (2,2);

    SetCell   (HY_CONSOLE_OUT_ROW,0,ow);
    SetCell   (HY_CONSOLE_OUT_ROW,1,ow);

    SetCell   (HY_CONSOLE_IN_ROW,0,iw);
    SetCell   (HY_CONSOLE_IN_ROW,1,bb);

    iw->SetBackColor (consoleIBarColor);
    bb->SetBackColor (consoleIBarColor);
    bb->SetButtonLayoutW (1);

    _HYFont         defFont;

#ifdef __MAC__
    defFont.face = "Monaco";
    defFont.size = 9;
#endif
#ifdef __WINDOZE__
    defFont.face = "MS Sans Serif";
    defFont.size = 12;
#endif
#ifdef __HYPHY_GTK__
    defFont.face = _HY_SANS_FONT;
    defFont.size = 10;
#endif

    defFont.style = HY_FONT_PLAIN;
    ow->SetFont (defFont);

    ow->margins = (_HYRect) {
        0,0,0,0,0
    };
    ow->SetText ("");
    iw->SetAlignFlags (HY_ALIGN_LEFT);
    iw->SetText ("Input");

    bb->SetAlignFlags (HY_ALIGN_LEFT);
    bb->SetButtonDim(16);

    _String s = "User Actions";
    bb->AddButton (ProcureIconResource(6006),&s);
    s = "Look-up Command";
    bb->AddButton (ProcureIconResource(6011),&s);
    s = "File echo";
    bb->AddButton (ProcureIconResource(6042),&s);
    s = "Web Update";
    bb->AddButton (ProcureIconResource(7050),&s);

    bb->EnableButton   (0,userHookins.lLength);
    bb->EnableButton   (1,false);
    bb->EnableButton   (2,true);
    bb->MarkAsPullDown (0,true);
    bb->MarkAsPullDown (1,true);
    bb->MarkAsPullDown (2,true);
    bb->EnableButton   (3,true);
    iw->SetMargins     ((_HYRect) {
        5,5,5,5,0
    });

    SetWindowRectangle (0,0,400,600);

    //keyboardFocusChain << 0;
    keyboardFocusChain << 1;

    echoFileRef     = nil;
    echoStatus      = 0;
    percentDone     = -1;
    timer           = "00:00:00";
    editOptions     = 0;
    inputStatus     = 0;

}

//__________________________________________________________
_HYConsoleWindow::~_HYConsoleWindow()
{
    if (echoFileRef) {
        fclose (echoFileRef);
    }
}

//__________________________________________________________

bool    _HYConsoleWindow::ConfirmClose (void)
{
    return false;
}

//__________________________________________________________

bool    _HYConsoleWindow::ProcessEvent (_HYEvent* e)
{
    bool    done = false;

    _String firstArg;
    long    i,f,k;

    if (e->EventClass().Equal(&_hyButtonPushEvent)) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }
        if (i==2) {
            firstArg        = e->EventCode().Cut (f+1,-1);
            k               = firstArg.toNum();
            _HYButtonBar * bb = (_HYButtonBar*)GetObject (2);

            switch  (k) {
            case 0:
            case 1: {

                _List menuChoices;
                if (k==0)
                    for (long k=0; k<userHookins.lLength; k++) {
                        _String *dName =  (_String*)userHookins(k),
                                 aName (*dName, dName->FindBackwards (
#ifdef __MAC__
                                            ':'
#else
#ifdef __WINDOZE__
                                            '\\'
#else
                                            '/'
#endif
#endif
                                            ,0,-1)+1,-1);

                        menuChoices && & aName;
                    }
                else {
                    firstArg = "Search Commands";
                    menuChoices && & firstArg;
                    firstArg = "Search Descriptions";
                    menuChoices && & firstArg;
                    firstArg = "Search Notes";
                    menuChoices && & firstArg;
                    firstArg = "Export to LaTeX by name";
                    menuChoices && & firstArg;
                    firstArg = menuSeparator;
                    menuChoices && & firstArg;
                    firstArg = "Open in a window";
                    menuChoices && & firstArg;
                }

                int h,v;
                bb->GetButtonLoc(k,h,v,true);
                _String userAction  = HandlePullDown (menuChoices,h,v,0),
                        justTheName = userAction;

                bb->_UnpushButton();
                i = menuChoices.Find (&userAction);
                if (i>=0) {
                    _ExecutionList uxl;
                    if (k==0) {
                        userAction = ((_String*)userHookins(i))->getStr();
                    } else {
                        userAction = baseDirectory&"Help"&baseDirectory.sData[baseDirectory.sLength-1]&"Commands"&baseDirectory.sData[baseDirectory.sLength-1]&"query.bf";
                    }

                    h = PushFilePath  (userAction);
                    ReadBatchFile (userAction, uxl);

                    SetStatusLine ( justTheName, "Loading", "00:00:00", -1);
                    if (k==1) {
                        _String addin;
                        switch (i) {
                        case 0:
                            addin = "Command";
                            break;
                        case 1:
                            addin = "Description";
                            break;
                        case 2:
                            addin = "Notes";
                            break;
                        case 3:
                            addin = "Export";
                            break;
                        case 5:
                            addin = "Window";
                            break;
                        }
                        firstArg = ((_HYTextBox*)GetObject (1))->GetText();
                        _FString qf = _FString (addin,false),
                                 qt = _FString (firstArg,false);

                        firstArg = "QUERY_FIELD";
                        setParameter (firstArg,&qf);
                        firstArg = "QUERY_TERM";
                        setParameter (firstArg,&qt);
                    }

                    uxl.Execute();
                    terminateExecution = false;
                    if (h) {
                        PopFilePath   ();
                    }
                    /*}
                    else
                    {
                        userAction = userAction & " could not be found. Perhaps the file was recently moved or deleted.";
                        ProblemReport (userAction, (Ptr)this);
                    }*/
                    SetStatusLine     (justTheName, "Finished", empty, -1, HY_SL_FILE|HY_SL_TASK);
                }
            }
            break;

            case 2:
                DoEcho();
                break;

            case 3: {
                _String webUpdate = libDirectory&"TemplateBatchFiles"&GetPlatformDirectoryChar()&"WebUpdate.bf";
                _ExecutionList wbl;
                k = PushFilePath  (webUpdate);
                ReadBatchFile (webUpdate, wbl);
                wbl.Execute();
                terminateExecution = false;
                if (k) {
                    PopFilePath   ();
                }

            }
            }

            done = true;
        }
    } else if (e->EventClass().Equal(&_hyTextEditChange)) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }
        if (i==0) { // out box
            _UpdateEditMenu();
            done = true;
        } else if (i==1) { // in box
            firstArg        = e->EventCode().Cut (f+1,-1);
            k               = firstArg.toNum();
            _HYTextBox * ib = (_HYTextBox*)GetObject (1);

            if (k==2) {
                if (inputStatus == 1) {
                    inputStatus = 2;
                }
            } else {
                if (k==3)
                    // down arrow
                {
                    if (inputLocation<recentInputs.lLength-1) {
                        ib->SetText (*(_String*)recentInputs(++inputLocation), false);
                    }
                } else if (k==4)
                    // up arrow
                {
                    if (inputLocation>0) {
                        ib->SetText (*(_String*)recentInputs(--inputLocation), false);
                    }
                } else {
                    if (inputStatus == 1) {
                        inputLocation = recentInputs.lLength;
                    }

                    _HYButtonBar * bb = (_HYButtonBar*)GetObject (2);
                    bool    onOff = ! ib->_IsEmpty();
                    bb->EnableButton   (1,onOff);
                }
            }
            done = true;
        }
    }

    if (done) {
        DeleteObject (e);
        return true;
    }
    return _HYTWindow::ProcessEvent(e);
}

//__________________________________________________________

void    _HYConsoleWindow::PrintString   (_String& s)
{
    if (s.sLength) {
        if (echoFileRef && echoStatus == 1) {
            fprintf (echoFileRef,"%s",s.sData);
        }
        //#ifndef __MAC__
        forceUpdateForScrolling = true;
        //#endif
        ((_HYTextBox *)components(0))->InsertText (s,true,true);
        //#ifndef __MAC__
        forceUpdateForScrolling = false;
        //#endif
#ifdef __WINDOZE__
        yieldCPUTime ();
#endif
    }
}


//__________________________________________________________

_String* _HYConsoleWindow::ReadString    (void)
{
    BringToFront ();
    _String * res = nil;

    _HYTextBox * inbox = ((_HYTextBox*)GetObject(1));

    inbox->SetText        (empty,false);

    ProcessEvent (generateKeyboardFocusEvent(inbox->GetID()));

    inputStatus   = 1;
    inputLocation = recentInputs.lLength;

    _PaintStatusBar();

    inbox->boxFlags |= HY_TB_ARROWS;

#ifndef __WINDOZE__
    while ((inputStatus != 2)&&(!terminateExecution)) {
        handleGUI();
    }
#else
    while ((inputStatus != 2)&&(!terminateExecution)&&(!hyphyExiting)) {
        handleGUI();
    }

    if (hyphyExiting) {
        return new _String;
    }
#endif

    inbox->boxFlags -= HY_TB_ARROWS;

    inbox->StoreText (res);

    inputStatus = 0;

    if (res->sLength) {
        if (recentInputs.lLength > 100) {
            recentInputs.Delete (0);
        }
        recentInputs && res;
    }

    _PaintStatusBar();

    if (echoFileRef && echoStatus == 1) {
        fprintf (echoFileRef,"\n>Input:%s",res->sData);
    }

    return res;
}

//__________________________________________________________

void _HYConsoleWindow::DoFind        (bool prompt)
{
    _String promptS ("Look for:");
    if (prompt)
        if (!EnterStringDialog (searchTerm, promptS, (Ptr)this)) {
            return;
        }

    if (searchTerm.sLength) {
        ((_HYTextBox*)GetObject(0))->_DoFind(searchTerm);
    }
}

//__________________________________________________________

void _HYConsoleWindow::DoSave        (void)
{
}



//__________________________________________________________

void _HYConsoleWindow::DoEcho        (void)
{
    _List           menuChoices;

    if (!echoFileRef) {
        menuChoices && & startEchoing;
    } else {
        menuChoices && & stopEchoing;

        if (echoStatus == 0) { // paused
            menuChoices && & resumeEchoing;
        } else {
            menuChoices && & pauseEchoing;
        }
    }

    int h,v;

    _HYButtonBar * bb = (_HYButtonBar*)GetObject (2);

    bb->GetButtonLoc(2,h,v,true);

    _String choice = HandlePullDown (menuChoices,h,v,0);

    if (choice.Equal (&startEchoing)) {
        _String saveDP (dialogPrompt),
                fileName;

        dialogPrompt = "Echo HyPhy console to:";
        fileName = WriteFileDialogInput ();
        if (fileName.sLength) {
            echoFileRef = doFileOpen (fileName.sData,"w");
            if (!echoFileRef) {
                fileName = _String ("File '")&fileName& "' couldn't be opened for writing.";
                ProblemReport (fileName);
            }
        }
        dialogPrompt = saveDP;
        echoStatus   = 1;
    } else if (choice.Equal (&stopEchoing)) {
        fclose (echoFileRef);
        echoFileRef = nil;
        echoStatus  = 0;
    } else if (choice.Equal (&pauseEchoing)) {
        echoStatus = 0;
    } else if (choice.Equal (&resumeEchoing)) {
        echoStatus = 1;
    }

    bb->_UnpushButton();
    _PaintStatusBar();
}

//_________________________________________________________________________

#define CONSOLE_BUFFERING  1024L

_String * consoleBuffer         = new _String (CONSOLE_BUFFERING, true);
long      consoleBufferAdded    = 0;

//_________________________________________________________________________

void    FlushConsoleBuffer (void)
{
    if (consoleBufferAdded) {
        consoleBuffer->Finalize();
        hyphyConsoleWindow->PrintString (*consoleBuffer);
        consoleBufferAdded = 0;
        DeleteObject (consoleBuffer);
        consoleBuffer     = new _String (CONSOLE_BUFFERING,true);
        checkPointer (consoleBuffer);
    }
}

//_________________________________________________________________________

void    StringToConsole (_String & s, _SimpleList* )
{

    if (s.sLength) {
        if (s.sLength + consoleBufferAdded >= CONSOLE_BUFFERING) {
            FlushConsoleBuffer ();
            hyphyConsoleWindow->PrintString (s);
        } else
            for (long k=0; k<s.sLength; k++, consoleBufferAdded++) {
                char c=s.sData[k];

                if ( c=='\n' || c=='\r' || consoleBufferAdded>=CONSOLE_BUFFERING) {
                    *consoleBuffer << c;
                    consoleBufferAdded++;
                    FlushConsoleBuffer ();
                } else {
                    *consoleBuffer << c;
                }
            }

        //hyphyConsoleWindow->PrintString (s);
    }
}

//_________________________________________________________________________

void    BufferToConsole (const char* buffer, _SimpleList* dummy)
{
    _String s (buffer);
    StringToConsole (s, dummy);
}

//_________________________________________________________________________

void    SaveConsole (void)
{

    _String     sDP = dialogPrompt;

    dialogPrompt = "Save HyPhy console to:";

    _String     fileName(WriteFileDialogInput ());

    if (fileName.sLength == 0) {
        return;
    }

    FILE *      f = doFileOpen (fileName.getStr(),"w");
    if (!f) {
        fileName = _String ("Could not open '") & fileName & "' for writing.";
        ProblemReport (fileName, (Ptr)hyphyConsoleWindow);
    } else {
        fileName = GetVersionString () & " console saved on " & GetTimeStamp () & "\n\n";
        fwrite (fileName.sData,1,fileName.sLength, f);
        _String *cS;
        ((_HYTextBox*)hyphyConsoleWindow->GetObject (0))->StoreText (cS);
#ifdef __WINDOZE__
        for (long k=0; k<cS->sLength; k++) {
            char c = cS->sData[k];
            if (c=='\r') {
                fprintf (f,"\n");
            } else {
                fputc (c,f);
            }
        }
#else
        fprintf (f,"%s",cS->sData);
#endif
        fclose (f);
        DeleteObject (cS);
    }

    dialogPrompt = sDP;
}

//_________________________________________________________________________

void    NLToConsole (void)
{
    StringToConsole (newLN);
}

//_________________________________________________________________________

_String*    StringFromConsole (bool echo)
{
    FlushConsoleBuffer();
    _String * res = hyphyConsoleWindow->ReadString ();
    if (echo) {
        hyphyConsoleWindow->PrintString (*res);
        hyphyConsoleWindow->PrintString (newLN);
    }
    return res;
}

//_________________________________________________________________________

void    SetStatusLine (_String arg)
{
#ifdef __HYPHYMPI__
    if (_hy_mpi_node_rank > 0) {
        return;
    }
#endif

    hyphyConsoleWindow->action = arg;
    hyphyConsoleWindow->_PaintStatusBar();

}

//_________________________________________________________________________

void    SetStatusBarValue (long l,_Parameter max, _Parameter rate)
{
#ifdef __HYPHYMPI__
    if (_hy_mpi_node_rank > 0) {
        return;
    }
#endif

    hyphyConsoleWindow->percentDone = l;
    if (max<=0.0) {
        hyphyConsoleWindow->action   = _String ("LF Optimization. Value=") & _String (max) &", "&_String (rate) & " evals/sec.";
    }
    hyphyConsoleWindow->_PaintStatusBar();
}

//_________________________________________________________________________

void    SetStatusLine (_String arg, _String arg2, _String arg3)
{
#ifdef __HYPHYMPI__
    if (_hy_mpi_node_rank > 0) {
        return;
    }
#endif

    hyphyConsoleWindow->action = arg2;
    hyphyConsoleWindow->fileName = arg;
    hyphyConsoleWindow->timer = arg3;
    hyphyConsoleWindow->_PaintStatusBar();
}

//_________________________________________________________________________

void    SetStatusLine (_String arg, _String arg2, _String arg3, long l)
{
#ifdef __HYPHYMPI__
    if (_hy_mpi_node_rank > 0) {
        return;
    }
#endif

    hyphyConsoleWindow->action = arg2;
    hyphyConsoleWindow->fileName = arg;
    hyphyConsoleWindow->timer = arg3;
    hyphyConsoleWindow->percentDone = l;
    hyphyConsoleWindow->_PaintStatusBar();
}

//_________________________________________________________________________

void    SetStatusLine (_String arg, _String arg2, _String arg3, long l, char c)
{
#ifdef __HYPHYMPI__
    if (_hy_mpi_node_rank > 0) {
        return;
    }
#endif

    if (c & HY_SL_FILE) {
        hyphyConsoleWindow->fileName = arg;
    }
    if (c & HY_SL_TASK) {
        hyphyConsoleWindow->action = arg2;
    }
    if (c & HY_SL_TIMER) {
        hyphyConsoleWindow->timer = arg3;
    }
    if (c & HY_SL_PERCENT) {
        hyphyConsoleWindow->percentDone = l;
    }
    if (c & HY_SL_SUSPEND) {
        hyphyConsoleWindow->action = suspendedPref & hyphyConsoleWindow->action;
    }
    if (c & HY_SL_RESUME) {
        hyphyConsoleWindow->action.Trim (suspendedPref.sLength,-1);
    }

    hyphyConsoleWindow->_PaintStatusBar(nil,(c & HY_SL_FORCE));
}

//_________________________________________________________________________

void    LoadUserHookins (void)
{
    _String baseDir = libDirectory & userHookinDir;
    ScanDirectoryForFileNames (baseDir,userHookins,false);
}


// EOF
