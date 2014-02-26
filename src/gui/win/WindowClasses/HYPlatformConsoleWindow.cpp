/*
    Win 32 Portions of the console window class

    Sergei L. Kosakovsky Pond, December 2003.
*/

#include "HYConsoleWindow.h"
#include "HYTextBox.h"
#include "HYDialogs.h"
#include "HYUtils.h"
#include "HYDataPanel.h"
#include "HYEventTypes.h"
#include "HYSharedMain.h"
#include "preferences.h"

#include "time.h"
#include "parser.h"

extern   HFONT   statusBarBold,
         statusBarPlain;

extern   bool    hyphyExiting;

//extern     _String selectedFileName,
//               savedFileName;


extern  time_t   timerStart,
        lastTimer;

extern  HMENU    recentFilesMenu;

LARGE_INTEGER    lastMeasure = {0,0};

void             RunStandardAnalyses (void);
void             displayAbout        (bool);
void             getUserFont         (void);
void             UpdateTimer         (void);

//__________________________________________________________________

void _HYConsoleWindow::_SetMenuBar(void)
{
    if (!menuKeys) {
        menuKeys = LoadAccelerators (ProgramInstance,"HYPHY");
    }
    //_HYWindow::_SetMenuBar();
}

//__________________________________________________________________

void _HYConsoleWindow::_UpdateEditMenu (void)
{

    _HYTextBox* txb = (_HYTextBox*)GetObject (0);

    HWND        te = txb->te;

    CHARRANGE   cr;

    SendMessage (te,EM_EXGETSEL,0,(LPARAM)&cr);

    bool        haveSelection = cr.cpMax-cr.cpMin,
                canPaste      = SendMessage (te,EM_CANPASTE,CF_TEXT,0),
                canUndo        = SendMessage (te,EM_GETUNDONAME,0,0),
                canRedo         = SendMessage (te,EM_CANREDO,0,0);

    if ((((bool)editOptions&HY_CONSOLE_CAN_COPY)!=haveSelection)||
            (((bool)editOptions&HY_CONSOLE_CAN_PASTE)!=canPaste)||
            (((bool)editOptions&HY_CONSOLE_CAN_UNDO)!=canUndo)||
            (((bool)editOptions&HY_CONSOLE_CAN_UNDO)!=canRedo)) {
        HMENU            windowMenu = GetMenu (theWindow),
                         editMenu   = GetSubMenu(windowMenu,1);

        editOptions = 0;

        if (haveSelection) {
            editOptions |= HY_CONSOLE_CAN_COPY;
        }

        EnableMenuItem (editMenu,3,MF_BYPOSITION|(haveSelection?MF_ENABLED:MF_GRAYED));
        EnableMenuItem (editMenu,4,MF_BYPOSITION|(haveSelection?MF_ENABLED:MF_GRAYED));

        if (canPaste) {
            editOptions |= HY_CONSOLE_CAN_PASTE;
        }

        EnableMenuItem (editMenu,5,MF_BYPOSITION|(canPaste?MF_ENABLED:MF_GRAYED));

        if (canUndo) {
            editOptions |= HY_CONSOLE_CAN_UNDO;
        }

        EnableMenuItem (editMenu,0,MF_BYPOSITION|(canUndo?MF_ENABLED:MF_GRAYED));

        if (canRedo) {
            editOptions |= HY_CONSOLE_CAN_REDO;
        }

        EnableMenuItem (editMenu,1,MF_BYPOSITION|(canRedo?MF_ENABLED:MF_GRAYED));
    }
}

//__________________________________________________________________

void _HYConsoleWindow::_UnsetMenuBar(void)
{
    _HYWindow::_UnsetMenuBar();
}

//__________________________________________________________________

bool        _HYConsoleWindow::_ProcessMenuSelection (long msel)
{
    switch (msel) {
    case 13:
        SaveConsole ();
        return true;

    case 14:
        _DoPrint ();
        return true;

    case 16:
        ((_HYTextBox*)GetObject(0))->_DoUndo(true);
        return true;

    case 17:
        ((_HYTextBox*)GetObject(0))->_DoCut(true);
        return true;

    case 18:
        ((_HYTextBox*)GetObject(0))->_DoCopy(true);
        return true;

    case 20:
        ((_HYTextBox*)GetObject(0))->_DoSelectAll(true);
        return true;

    case 21:
        SetStatusLine ("Canceling");
        terminateExecution = true;
        return true;

    case 22:
        if (!isSuspended) {
            isSuspended = true;
            SetStatusLine (empty,empty,empty,-1,HY_SL_SUSPEND);
            HMENU       sMenu = GetSubMenu (GetMenu((HWND)GetOSWindowData()),2);
            if (sMenu) {
                EnableMenuItem(sMenu, 0 ,MF_GRAYED|MF_BYPOSITION);
            }
            ModifyMenu (sMenu,1,MF_BYPOSITION,22,"&Resume\tCtrl+;");
            updateTimer = false;
        } else {
            isSuspended = false;
            SetStatusLine (empty,empty,empty,-1,HY_SL_RESUME);
            HMENU       sMenu = GetSubMenu (GetMenu((HWND)GetOSWindowData()),2);
            if (sMenu) {
                EnableMenuItem(sMenu, 0 ,MF_ENABLED|MF_BYPOSITION);
            }

            ModifyMenu (sMenu,1,MF_BYPOSITION,22,"&Suspend Execution\tCtrl+;");
            time_t      tt;
            timerStart += time(&tt)-lastTimer;
            updateTimer = true;
        }
        return true;

    case 23:
        ShowMessagesLog();
        return true;

    case 24:
        RunStandardAnalyses();
        return true;

    case 25:
        displayAbout(false);
        return true;

    case 26:
        ((_HYTextBox*)GetObject(0))->_DoClear (true,true);
        return true;

    case 27:
        HandlePreferences (globalPreferencesList,"HYPHY Preferences");
        return true;

    case 28:
        getUserFont();
        return true;

    case 29:
        WinExec ("hh HYPHY HELP.chm",SW_SHOWNORMAL);
        return true;

    case 30:
        if (OpenBatchFile (false)) {
            ExecuteBatchFile();
        }
        return true;

    case 31:
        ((_HYTextBox*)GetObject(0))->_DoRedo(true);
        return true;

    case 15:
        postWindowCloseEvent (GetID());
        return true;

    case 60: // expression calculator
        if (calculatorMode) {
            _HYTextBox         *ib = (_HYTextBox*)hyphyConsoleWindow->GetObject(1);
            ib->SetText ("exit");
            hyphyConsoleWindow->ProcessEvent (generateTextEditChangeEvent(ib->GetID(),2));
            calculatorMode         = false;
            //ib->SetText (empty);
        } else {
            calculatorMode = true;
            while(!hyphyExiting && calculatorMode&&ExpressionCalculator()) {}
            calculatorMode = false;
        }
        return true;

    case 61: { // execute selection
        ExecuteSelection();
        return true;
    }

    case 70: // New Tree
        NewTreeWindow(-1);
        return true;

    case 71: // New Model
        NewModel(nil);
        return true;

    case 72: // New Chart
        NewChartWindow();
        return true;

    case 73: // New Genetic Code
        NewGeneticCodeTable(0);
        return true;

    case 74: // New Database
        NewDatabaseFile (nil);
        return true;

    case 80: // Open Batch File
        if (OpenBatchFile()) {
            ExecuteBatchFile ();
        }
        return true;

    case 81: // Open Data File
        OpenDataFile();
        return true;

    case 82: // Open Tree
        OpenTreeFile();
        return true;

    case 83: // Open Text
        OpenTextFile();
        return true;

    case 84: // Open Table
        OpenTable ();
        return true;

    case 85: // Open Database
        OpenDatabaseFile (nil);
        return true;

    case 200:
        hyphyConsoleWindow->_BringWindowToFront();
        return true;

    case 201:
        ShowObjectInspector ();
        return true;

    default: {
        msel -= 1000;
        if (msel<availablePostProcessors.lLength) {
            ExecuteAPostProcessor (*(_String*)(*(_List*)availablePostProcessors(msel))(1));
            return 0;
        }
        msel-=1000;
        if (msel<(long)recentPaths.lLength) {
            if (msel == -1) {
                for (long mi = 0; mi < recentPaths.lLength; mi++) {
                    DeleteMenu (recentFilesMenu,2,MF_BYPOSITION);
                }

                recentPaths.Clear();
                recentFiles.Clear();
                DrawMenuBar ((HWND)GetOSWindowData());
            } else {
                if (argFileName) {
                    *argFileName = *(_String*)recentPaths(msel);
                } else {
                    argFileName = new _String (*(_String*)recentPaths(msel));
                }
                if (OpenBatchFile(false)) {
                    ExecuteBatchFile ();
                }
            }
            return true;
        }

        return true;
    }
    }

    return _HYTWindow::_ProcessMenuSelection(msel);
}

//__________________________________________________________________

bool _HYConsoleWindow::_ProcessOSEvent (Ptr vEvent)
{
    _HYWindowsUIMessage*    theEvent = (_HYWindowsUIMessage*)vEvent;
    if (theEvent->iMsg == WM_SYSCOMMAND) {
        if (theEvent->wParam == SC_CLOSE) {
            hyphyExiting = true;
            return true;
        }
    } else {
        if (theEvent->iMsg == UPDATE_TIMER) {
            UpdateTimer ();
            return true;
        }
    }
    return _HYTWindow::_ProcessOSEvent (vEvent);
}

//__________________________________________________________________

bool _HYConsoleWindow::_Close (Ptr )
{
    hyphyExiting = true;
    return false;
}

//__________________________________________________________________
void    _HYConsoleWindow::_DoPrint          (void)
{
    _HYTextBox* ob = (_HYTextBox*)components(0);

    DOCINFO                 di = {sizeof(DOCINFO), "HYPHY.out", NULL };
    PRINTDLG                pd;
    BOOL                    SuccessFlag;

    pd.lStructSize         = sizeof(PRINTDLG);
    pd.hwndOwner           = theWindow;
    pd.hDevMode            = NULL;
    pd.hDevNames           = NULL;
    pd.hDC                 = NULL;
    pd.Flags               = PD_COLLATE | PD_RETURNDC | PD_NOSELECTION;
    pd.nFromPage           = 1;
    pd.nToPage             = 0xffff;
    pd.nMinPage            = 1;
    pd.nMaxPage            = 0xffff;
    pd.nCopies             = 1;
    pd.hInstance           = NULL;
    pd.lCustData           = 0L;
    pd.lpfnPrintHook       = NULL;
    pd.lpfnSetupHook       = NULL;
    pd.lpPrintTemplateName = NULL;
    pd.lpSetupTemplateName = NULL;
    pd.hPrintTemplate      = NULL;
    pd.hSetupTemplate      = NULL;

    if (!PrintDlg(&pd)) {
        return;
    }

    if (pd.hDC == NULL) {
        pd.hDC = GetPrinterDeviceContext(theWindow);
    }

    EnableWindow(theWindow, FALSE);

    SuccessFlag   = TRUE;
    UserAbortFlag = FALSE;

    PrintDialogHandle = CreateDialog(GetModuleHandle(NULL), (LPCTSTR)"PrintDlgBox", theWindow,
                                     PrintDialogProc);
    SetDlgItemText(PrintDialogHandle, IDD_FNAME, "Table Printing...");

    SetAbortProc(pd.hDC, AbortProc);

    if (StartDoc(pd.hDC, &di) > 0) {

        long        fromPage = pd.nMinPage,
                    toPage   = pd.nMaxPage;

        _Parameter  CFW = 1440./GetDeviceCaps(pd.hDC, LOGPIXELSX),
                    CFH = 1440./GetDeviceCaps(pd.hDC, LOGPIXELSY);



        FORMATRANGE  fr;

        fr.rc.left  = GetDeviceCaps(pd.hDC, PHYSICALOFFSETX)    * CFW;
        fr.rc.top   = GetDeviceCaps(pd.hDC, PHYSICALOFFSETY)    * CFH;
        fr.rc.right  = fr.rc.left + GetDeviceCaps(pd.hDC, HORZRES) * CFW;
        fr.rc.bottom = fr.rc.top + GetDeviceCaps(pd.hDC, VERTRES) * CFH;

        fr.rcPage.left  = fr.rcPage.top = 0;
        fr.rcPage.right = GetDeviceCaps(pd.hDC, PHYSICALWIDTH) * CFW;
        fr.rcPage.bottom = GetDeviceCaps(pd.hDC, PHYSICALHEIGHT) * CFH;
        fr.chrg.cpMin = 0;

        GETTEXTLENGTHEX gtl = {GTL_NUMCHARS,CP_ACP};
        fr.chrg.cpMax = SendMessage (ob->te,EM_GETTEXTLENGTHEX,(WPARAM)&gtl,0)-1;

        while (fr.chrg.cpMin < fr.chrg.cpMax) {
            if (StartPage (pd.hDC) <= 0) {
                SuccessFlag = FALSE;
                break;
            }

            SetMapMode     (pd.hDC, MM_TEXT);

            fr.hdc          = pd.hDC;
            fr.hdcTarget    = pd.hDC;

            fr.chrg.cpMin = SendMessage(ob->te,EM_FORMATRANGE,1,(LPARAM)&fr);
            SendMessage(ob->te,EM_DISPLAYBAND,1,(LPARAM)&fr.rcPage);
            SendMessage(ob->te,EM_FORMATRANGE,0,nil);

            if (EndPage (pd.hDC) <= 0) {
                SuccessFlag = FALSE;
            }
        }
    } else {
        SuccessFlag = FALSE;
    }


    if (SuccessFlag) {
        SuccessFlag = (EndDoc(pd.hDC)>0);
    }

    if (!UserAbortFlag) {
        EnableWindow(theWindow, TRUE);
        DestroyWindow(PrintDialogHandle);
    }

    DeleteDC (pd.hDC);

    if (!SuccessFlag && !UserAbortFlag) {
        _String errMsg = _String("Failed to print console. Windows Error:") & (long)GetLastError();
        ProblemReport (errMsg,nil);
    }
}

//__________________________________________________________________

void _HYConsoleWindow::_PaintStatusBar(Ptr hdp, bool force)
{
    _Parameter      vL;
    checkParameter (VerbosityLevelString, vL, 0.0);
    if (vL<-0.5 && !force) {
        LARGE_INTEGER curMeasure;
        QueryPerformanceCounter (&curMeasure);

        _Parameter      timeDiff   = (curMeasure.QuadPart-lastMeasure.QuadPart)/1000000;

        if (timeDiff < -vL) {
            return;
        }
        lastMeasure = curMeasure;
    }

    HDC  osc = CreateCompatibleDC (NULL),
         hdc;

    if (hdp) {
        hdc = (HDC)hdp;
    } else {
        hdc = GetDC (theWindow);
    }


    RECT wRC,
         w2RC;

    GetClientRect (theWindow, &wRC);
    wRC.top    = 0;
    wRC.bottom = HY_SCROLLER_WIDTH;

    HBITMAP oBM = CreateCompatibleBitmap (hdc,wRC.right,HY_SCROLLER_WIDTH),
            oldBM;

    //printf ("%x %x\n", osc, oBM);

    if (oBM) {
        oldBM = (HBITMAP)SelectObject (osc,oBM);
        HBRUSH bkBrush = CreateSolidBrush (GetBkColor(hdc));

        FillRect     (osc, &wRC,windowStatusBarBrush);
        MoveToEx     (osc,0,wRC.bottom-HY_SCROLLER_WIDTH+1,NULL);
        SelectObject (osc, (HPEN)GetStockObject(BLACK_PEN));
        LineTo       (osc,wRC.right,wRC.bottom-HY_SCROLLER_WIDTH+1);

        SetBkMode    (osc, TRANSPARENT);
        SetTextAlign (osc, TA_LEFT|TA_BOTTOM|TA_NOUPDATECP);
        SetTextColor (osc,RGB(0,0,0));
        SelectObject (osc,statusBarPlain);
        TextOut      (osc,33,wRC.bottom-1, fileName.getStr(),fileName.sLength);
        if (inputStatus == 1) {
            TextOut      (osc,193,wRC.bottom-1,cInput.getStr(),cInput.sLength);
        } else {
            TextOut      (osc,193,wRC.bottom-1,action.getStr(),action.sLength);
        }

        w2RC = wRC;

        w2RC.right = 30;
        w2RC.top=w2RC.bottom-HY_SCROLLER_WIDTH+2;

        FillRect(osc, &w2RC, (HBRUSH)GetStockObject (DKGRAY_BRUSH));
        w2RC.left  = 150;
        w2RC.right = 190;
        FillRect(osc, &w2RC, (HBRUSH)GetStockObject (DKGRAY_BRUSH));
        w2RC.right = wRC.right;
        w2RC.left = w2RC.right - 50;
        FillRect(osc, &w2RC, (HBRUSH)GetStockObject (DKGRAY_BRUSH));
        SelectObject (osc,statusBarBold);
        SetTextColor (osc,RGB(255,255,255));
        TextOut (osc,3,wRC.bottom-1,cState.sData,cState.sLength);
        TextOut (osc,151,wRC.bottom-1,cTask.sData,cTask.sLength);
        SelectObject (osc,statusBarPlain);

        TextOut (osc,wRC.right-48,wRC.bottom-1,timer.getStr(),timer.sLength);

        if (percentDone>=0 || percentDone == -HY_SL_DONE) {
            HBRUSH blackBrush = CreateSolidBrush (RGB(80,80,80)), orangeBrush = CreateSolidBrush (RGB(255,153,102));
            w2RC.right = w2RC.left-5;
            w2RC.left-=75;
            w2RC.top++;
            w2RC.bottom--;

            w2RC.right = w2RC.left+(w2RC.right-w2RC.left)*(percentDone>=0?percentDone:100)/100;
            FillRect (osc,&w2RC,orangeBrush);
            w2RC.right = w2RC.left+70;
            FrameRect (osc,&w2RC,blackBrush);
            _String pLine;
            if (percentDone>=0) {
                pLine = _String(percentDone)&"%";
            } else {
                pLine = "DONE";
            }

            SetTextColor (osc,RGB(0,0,0));
            TextOut (osc,w2RC.left+28,wRC.bottom-1,pLine.getStr(),pLine.sLength);

            DeleteObject (blackBrush);
            DeleteObject (orangeBrush);
        }

        GetClientRect (theWindow,&w2RC);

        BitBlt  (hdc,0,w2RC.bottom-HY_SCROLLER_WIDTH,wRC.right,HY_SCROLLER_WIDTH,osc, 0,0,SRCCOPY);
        SelectObject (osc,oldBM);
        DeleteObject (bkBrush);
        DeleteObject (oBM);
    }
    DeleteDC     (osc);


    if (!hdp) {
        ReleaseDC (theWindow,hdc);
    }
    
    yieldCPUTime();
}

//EOF