#include "batchlan.h"
#include "preferences.h"
#include "HYSharedMain.h"
#include "HYUtils.h"
#include "HYDialogs.h"
#include "HYChartWindow.h"
#include "HYModelWindow.h"
#include "HYConsoleWindow.h"
#include "HYCanvas.h"
#include "HYTextBox.h"
#include "HYDataPanel.h"
#include "HYEventTypes.h"

#include <io.h>
#include <windows.h>
#include <winuser.h>
#include <wingdi.h>
#include <commdlg.h>
#include <time.h>
#include <stdio.h>

#ifndef _WIN32_IE
#define _WIN32_IE 0x0400
#endif

#include <commctrl.h>

#define  RECENT_FILE_ITEMS       10

//____________________________________________________________________________________________
/* EXTERNAL GLOBAL VARIABLES */
//____________________________________________________________________________________________


extern           _Parameter         printDigits;
extern           _String            dialogPrompt;

//____________________________________________________________________________________________
/* GLOBAL VARIABLES */
//____________________________________________________________________________________________

HMENU                   recentFilesMenu,
                        resultsMenu         = nil;

char                    AppName[] = "HYPHY";
char                    WinName[] = "HPGUIW";


_String                 lastWinPathUsed,
                        consoleWindowTitle ("HYPHY Console");


bool                    onlyStatusBar       = false,
                        hyphyExiting        = false;

time_t                  timerStart          = 0,
                        lastTimer             = 0;

_SimpleList             *currentListBoxSel,
                        *currentListBoxAvChoices,
                        *currentListBoxRec;

long                    currentListBoxSelLength,
                        barPercentage = -1;


HFONT                   statusBarBold,
                        statusBarPlain;

HBRUSH                  windowStatusBarBrush,
                        dlgWhite = CreateSolidBrush (RGB (255,255,255));

HDC                     offScreenContext;

HINSTANCE               ProgramInstance;

BOOL                    UserAbortFlag;
HWND                    PrintDialogHandle;


//__________________________________________________________________________________
// RESOURCE IMPORT DEFS
//__________________________________________________________________________________

HCURSOR                 hSizeCursor,
                        pickUpCursor,
                        dropOffCursor;

HBITMAP                 tablePDMenuIcon,
                        redButtonIcon,
                        yellowButtonIcon,
                        orangeButtonIcon,
                        greenButtonIcon;

HDC                     otherDC;
HBITMAP                 oDCBM;


//____________________________________________________________________________________________
/* FUNCTION PROTOTYPES */
//____________________________________________________________________________________________


LRESULT         CALLBACK HYWndProc               (HWND , UINT , WPARAM  , LPARAM );
LRESULT         CALLBACK HYDlgProc               (HWND , UINT , WPARAM  , LPARAM );
VOID            CALLBACK GlobalQueueTimer        (HWND , UINT , UINT_PTR, DWORD );
BOOL            CALLBACK startupBoxHandlerFunction
(HWND , UINT , WPARAM , LPARAM );

BOOL            CALLBACK aboutBoxHandlerFunction (HWND , UINT , WPARAM , LPARAM );
BOOL            CALLBACK PrintDialogProc         (HWND, UINT, WPARAM, LPARAM);
BOOL            CALLBACK AbortProc               (HDC, int);
HDC             GetPrinterDeviceContext          (HWND);

void            LoadHYPHYResources               (void);
void            UpdateTimer                      (void);
bool            displayFirstDialog               (long&);

void            WinInitSet                       (void);
void            WinErrorBox                      (_String&, bool = false);
bool            runFile                          (void);
void            displayAbout                     (bool);
void            centerDialogBox                  (HWND);
void            getUserFont                      (void);

//__________________________________________________________________________________

BOOL CALLBACK PrintDialogProc(HWND DialogHandle, UINT msg, WPARAM /* wParam */, LPARAM /*lParam */)
{
    switch(msg) {
    case WM_INITDIALOG:
        EnableMenuItem(GetSystemMenu(DialogHandle, FALSE), SC_CLOSE, MF_GRAYED);
        return TRUE;

    case WM_COMMAND:
        UserAbortFlag = TRUE;
        EnableWindow(GetParent(DialogHandle), TRUE);
        DestroyWindow(DialogHandle);
        PrintDialogHandle = 0;
        return TRUE;
    }
    return FALSE;
}

//__________________________________________________________________________________

void        UpdateTimer (void)
{
    time_t tt,
           time_diff = time(&tt)-lastTimer;

    if (time_diff>=1) { // update dispay
        lastTimer = tt;

        _String     r;
        updateTimerF(r,tt-timerStart);
        SetStatusLine (empty,empty,r,0,HY_SL_TIMER);
    }
}


//__________________________________________________________________________________

HDC GetPrinterDeviceContext(HWND WindowHandle)
{
    PRINTER_INFO_4  * PrinterInfo4;             /* mm 990507 */
    PRINTER_INFO_5  * PrinterInfo5;             /* mm 990507 */
    DWORD           BytesCopied, StructuresCopied;
    char            Msg[]  = "Default printer not found";
    HDC             hdc;
    /* begin mm 990507 */
    if (GetVersion() & 0x80000000) { /*  Windows 98  */

        EnumPrinters(PRINTER_ENUM_DEFAULT, NULL, 5, NULL, 0, &BytesCopied,  &StructuresCopied);

        PrinterInfo5 = (PRINTER_INFO_5 *)malloc(BytesCopied);

        if (EnumPrinters(PRINTER_ENUM_DEFAULT, NULL, 5, (PBYTE) PrinterInfo5,
                         BytesCopied, &BytesCopied, &StructuresCopied)) {
            hdc = CreateDC(NULL, PrinterInfo5->pPrinterName, NULL, NULL);
            free(PrinterInfo5);
        } else {

            MessageBox(WindowHandle, Msg, "HYPHY",
                       MB_OK | MB_ICONEXCLAMATION);
            return (0);
        }

    } else {
        /*  Windows NT || Windows 95 */

        EnumPrinters(PRINTER_ENUM_LOCAL, NULL, 4, NULL, 0, &BytesCopied, &StructuresCopied);

        PrinterInfo4 = (PRINTER_INFO_4 *)malloc(BytesCopied);

        if (EnumPrinters(PRINTER_ENUM_LOCAL, NULL, 4, (PBYTE) PrinterInfo4,
                         BytesCopied, &BytesCopied, &StructuresCopied)) {
            hdc = CreateDC(NULL, PrinterInfo4->pPrinterName, NULL, NULL);
            free(PrinterInfo4);

        } else {

            MessageBox(WindowHandle, Msg, "HYPHYKernel",
                       MB_OK | MB_ICONEXCLAMATION);
            return (0);
        }
    }
    /* end mm 990507 */

    return(hdc);
}


//__________________________________________________________________________________

BOOL CALLBACK AbortProc(HDC /* hPrinterDC */, int /* iCode */)
{
    MSG msg;

    while (!UserAbortFlag && PeekMessage (&msg, NULL, 0, 0, PM_REMOVE)) {
        if (!PrintDialogHandle || !IsDialogMessage(PrintDialogHandle, &msg)) {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }
    return !UserAbortFlag;
}

//_________________________________________________________________________
void  getUserFont (void)
{
    LOGFONT fInfo;
    CHOOSEFONT cInfo;
    cInfo.lStructSize = sizeof (CHOOSEFONT);
    cInfo.hwndOwner = (HWND)hyphyConsoleWindow->GetOSWindowData();
    cInfo.lpLogFont = &fInfo;
    cInfo.hDC = NULL;
    cInfo.Flags = CF_BOTH|CF_LIMITSIZE|CF_FIXEDPITCHONLY;
    cInfo.nSizeMin = 8;
    cInfo.nSizeMax = 14;
    cInfo.lCustData = 0;
    cInfo.hInstance = NULL;
    cInfo.lpszStyle = (LPTSTR)NULL;
    cInfo.nFontType = SCREEN_FONTTYPE;

    if (ChooseFont(&cInfo)) {
        _HYFont cF;
        cF.face     = fInfo.lfFaceName;
        cF.size     = cInfo.iPointSize/10;
        cF.style    = HY_FONT_PLAIN;

        ((_HYTextBox*)hyphyConsoleWindow->GetObject(0))->SetFont(cF);
        ((_HYTextBox*)hyphyConsoleWindow->GetObject(1))->SetFont(cF);
    }
}

//_________________________________________________________________________
BOOL CALLBACK startupBoxHandlerFunction (HWND DialogHandle, UINT iMsg, WPARAM wParam, LPARAM )
{
    switch (iMsg) {
    case WM_COMMAND: {
        if ((HIWORD (wParam)==STN_CLICKED)||(HIWORD (wParam)==BN_CLICKED)) {
            if (LOWORD(wParam)) {
                HWND   checkBoxHandle = GetDlgItem (DialogHandle,0);
                char   isChecked = (SendMessage (checkBoxHandle, BM_GETCHECK, 0,0)==BST_CHECKED);
                EndDialog (DialogHandle, LOWORD (wParam)+5*isChecked);
                return 1;
            }
        }
        break;
    }

    case WM_CLOSE:
    case WM_DESTROY: {
        EndDialog (DialogHandle,0);
        return 1;
    }
    case WM_INITDIALOG: {
        _String theMessage = GetVersionString();
        SetDlgItemText (DialogHandle, 5, theMessage.sData);
        centerDialogBox (DialogHandle);
        return 1;
        break;
    }
    case WM_ERASEBKGND: {
        RECT cr;
        GetClientRect (DialogHandle, &cr);
        FillRect ((HDC)wParam, &cr, dlgWhite);
        return 1;
    }
    case WM_CTLCOLORBTN:
    case WM_CTLCOLORSTATIC:
        return (LRESULT)dlgWhite;
    }
    return 0;
}

//_________________________________________________________________________

bool            runFile (void)
{
    return argFileName && argFileName->sLength;
}

//_________________________________________________________________________
void    WinInitSet (void)
{
    GlobalStartup();

    ReadInTemplateFiles     ();
    ReadModelTemplates      ();
    ReadGeneticCodes        ();

    {
        SYSTEM_INFO      sysInfo;
        GetSystemInfo (&sysInfo);

        systemCPUCount = sysInfo.dwNumberOfProcessors;

        if (systemCPUCount < 1) {
            systemCPUCount = 1;
        }

        _String procType;

        switch (sysInfo.wProcessorArchitecture) {
        case PROCESSOR_ARCHITECTURE_INTEL:
            procType = " Intel x86 architecture ";
            break;
        case PROCESSOR_ARCHITECTURE_IA64:
            procType = " Intel 64-bit IA (Itanium) architecture ";
            break;

        case PROCESSOR_ARCHITECTURE_AMD64:
            procType = " AMD64 architecture ";
            break;
        }

        if (systemCPUCount == 1) {
            procType = _String ("One ") & procType & " processor detected\n";
            StringToConsole (procType);
        } else {
            procType = _String ('\n') & _String ((long)systemCPUCount) & procType & " processors detected.\n You can change the number of processors utilized by HyPhy via the HyPhy settings dialog.\n";
            StringToConsole (procType);
        }
    }

    StringToConsole (hyphyCiteString);

    if (availablePostProcessors.lLength) {
        resultsMenu = CreatePopupMenu();
        if (!resultsMenu) {
            return;
        }

        for (long counter=0; counter<availablePostProcessors.countitems(); counter++) {
            _String * postItem = (_String*)(*(_List*)availablePostProcessors(counter))(0);
            if (postItem->Equal(&menuSeparator)) {
                InsertMenu (resultsMenu,-1,MF_BYPOSITION|MF_SEPARATOR,0,(LPSTR)0);
            } else {
                InsertMenu (resultsMenu,-1,MF_BYPOSITION,1000+counter,(LPSTR)postItem->getStr());
            }
        }

        HMENU               sMenu = GetSubMenu (GetMenu((HWND)hyphyConsoleWindow->GetOSWindowData()),2);
        ModifyMenu          (sMenu,6,MF_BYPOSITION|MF_POPUP,(UINT)resultsMenu,"&Results");
        EnableMenuItem      (sMenu,6,MF_BYPOSITION|MF_GRAYED);
        recentFilesMenu  = CreatePopupMenu();
        checkPointer (recentFilesMenu);
        InsertMenu (recentFilesMenu,-1,MF_BYPOSITION,1999,"&Clear Menu");
        InsertMenu (recentFilesMenu,-1,MF_BYPOSITION|MF_SEPARATOR,0,(LPSTR)0);
        ModifyMenu (GetSubMenu(GetSubMenu (GetMenu((HWND)hyphyConsoleWindow->GetOSWindowData()),0),1),4,MF_BYPOSITION|MF_POPUP,(UINT)recentFilesMenu,"Open &Recent");
        DrawMenuBar((HWND)hyphyConsoleWindow->GetOSWindowData());
    }

}

//__________________________________________________________________________________

void    LoadHYPHYResources (void)
{
    tablePDMenuIcon  = LoadBitmap (GetModuleHandle(NULL),MAKEINTRESOURCE(4020));
    pickUpCursor     = LoadCursor (GetModuleHandle(NULL),MAKEINTRESOURCE(128));
    dropOffCursor    = LoadCursor (GetModuleHandle(NULL),MAKEINTRESOURCE(129));
    hSizeCursor      = LoadCursor (GetModuleHandle(NULL),MAKEINTRESOURCE(132));

    redButtonIcon    = LoadBitmap (GetModuleHandle(NULL),MAKEINTRESOURCE(4000));
    yellowButtonIcon = LoadBitmap (GetModuleHandle(NULL),MAKEINTRESOURCE(4001));
    greenButtonIcon  = LoadBitmap (GetModuleHandle(NULL),MAKEINTRESOURCE(4002));
    orangeButtonIcon = LoadBitmap (GetModuleHandle(NULL),MAKEINTRESOURCE(4003));


    otherDC         = CreateCompatibleDC (NULL);
    oDCBM           = CreateCompatibleBitmap (otherDC,1,1);
    DeleteObject    (SelectObject (otherDC, oDCBM));
}

//__________________________________________________________________

VOID CALLBACK GlobalQueueTimer (HWND , UINT , UINT_PTR idEvent, DWORD )
{
    if (idEvent == GLOBAL_QUEUE_TIMER) {
        if (GlobalGUIEventQueue.lLength) {
            HandleGlobalQueueEvent ();
        }
        if (updateTimer) {
            SendMessage ((HWND)hyphyConsoleWindow->GetOSWindowData(), UPDATE_TIMER, 0, 0);
        }
    }
}


//_________________________________________________________________________
bool displayFirstDialog (long& choice)
{
    long res = DialogBoxParam(ProgramInstance, MAKEINTRESOURCE(107), (HWND)hyphyConsoleWindow->GetOSWindowData(),  startupBoxHandlerFunction, nil);
    choice = res%5;
    return !(res/5);
}

//_________________________________________________________________________

int MessageLoop         (bool peek, bool moreThanOne)
{
    MSG msg;
    int loopMore = peek?PeekMessage(&msg, NULL, 0, 0, PM_REMOVE):GetMessage (&msg, NULL, 0, 0);
    while (loopMore) {
        _HYWindow* me = nil;
        HWND       daWindow = msg.hwnd;
        long f = windowPtrs.Find((long)daWindow);

        if (f<0) {
            daWindow = GetParent (msg.hwnd);
            f = windowPtrs.Find ((long)daWindow);
        }

        if (f>=0) {
            me = (_HYWindow*)(windowObjectRefs.lData[f]);
            if ( me->menuKeys!=NULL && TranslateAccelerator (daWindow, me->menuKeys, &msg) && moreThanOne) {
                loopMore = peek?PeekMessage(&msg, NULL, 0, 0, PM_REMOVE):GetMessage (&msg, NULL, 0, 0);
                continue;
            }

        }
        TranslateMessage(&msg);
        DispatchMessage(&msg);

        if (hyphyExiting) {
            return 0;
        }
        if (moreThanOne) {
            loopMore = peek?PeekMessage(&msg, NULL, 0, 0, PM_REMOVE):GetMessage (&msg, NULL, 0, 0);
        } else {
            return msg.wParam;
        }
    }
    return msg.wParam ;
}

//_________________________________________________________________________

void     WinErrorBox (_String& errText, bool warn)
{
    if (warn) {
        errText = errText & " Current Task has been terminated. Would you like to see the remaining error messages, if there are any?";
        int retVal = MessageBox(hyphyConsoleWindow?(HWND)hyphyConsoleWindow->GetOSWindowData():nil,errText.getStr(), AppName, MB_TASKMODAL | MB_ICONWARNING | MB_YESNO);
        if (retVal==IDNO) {
            skipWarningMessages = true;
        }
    } else {
        MessageBox(hyphyConsoleWindow?(HWND)hyphyConsoleWindow->GetOSWindowData():nil,errText.getStr(), AppName, MB_TASKMODAL | MB_ICONERROR);
    }
}

//_________________________________________________________________________
BOOL CALLBACK aboutBoxHandlerFunction (HWND DialogHandle, UINT iMsg, WPARAM wParam, LPARAM lParam)
{
    long k;
    static _List  theStrings;
    static long   currentTopPosition = -1, currentFirstLine = 0, currentDrawnLines;
    static bool   splash = false;

    //static _String theText;
    static  HFONT plainFont,boldFont;
    switch (iMsg) {
    case WM_COMMAND: {
        switch (HIWORD (wParam)) {
        case BN_CLICKED: {
            if (LOWORD(wParam)==1) { // selection
                EndDialog (DialogHandle,0);
                return 1;
            }
            break;
        }
        }
        break;
    }
    case WM_CLOSE:
    case WM_DESTROY: {
        //if (currentTopPosition!=-1)
        KillTimer (DialogHandle,5);
        EndDialog (DialogHandle,-1);
        DeleteObject (plainFont);
        DeleteObject (boldFont);
        return 1;
    }
    case WM_PAINT: {
        PAINTSTRUCT      thePS;
        HWND theBox = GetDlgItem (DialogHandle,129);
        HDC theContext = BeginPaint (theBox,&thePS);
        RECT bRC;
        SIZE stringSize;
        GetClientRect (theBox,&bRC);
        bool paintLast = false;

        long v = bRC.bottom-currentTopPosition, h = (bRC.left+bRC.right)/2;
        SetTextAlign (theContext,TA_NOUPDATECP |TA_TOP |TA_CENTER);
        k = currentFirstLine;
        SetTextColor (theContext, GetSysColor (COLOR_WINDOWTEXT));
        SetBkMode (theContext, TRANSPARENT);
        //SetBkColor (theContext, GetSysColor (COLOR_WINDOW));
        if (currentTopPosition == -1) {
            v = 1;
//              FillRect (theContext, &bRC, (HBRUSH)GetStockObject(WHITE_BRUSH));
            while((v<bRC.bottom)&&(k<theStrings.lLength)) {
                _String * thisLine = (_String*)theStrings.lData[k];
                GetTextExtentPoint32(theContext,thisLine->sData, thisLine->sLength,&stringSize);
                if (thisLine->sData[0]=='_') {
                    SelectObject (theContext, boldFont);
                    TextOut(theContext,h,v,thisLine->sData+1, thisLine->sLength-1);
                } else {
                    SelectObject (theContext, plainFont);
                    TextOut(theContext,h,v,thisLine->sData, thisLine->sLength);
                }
                v+=stringSize.cy+1;
                k++;
            }
        } else {
//              if (currentTopPosition == 15)
//                  FillRect (theContext, &bRC, (HBRUSH)GetStockObject(WHITE_BRUSH));

            if (currentFirstLine<theStrings.lLength) {
                ScrollWindow (theBox,0,-1,NULL,NULL);
                _String * thisLine;;
                while((v<bRC.bottom)&&(k<currentFirstLine+currentDrawnLines)&&(k<theStrings.lLength)) {
                    thisLine = (_String*)theStrings.lData[k];
                    GetTextExtentPoint32(theContext,thisLine->sData, thisLine->sLength,&stringSize);
                    v+=stringSize.cy+1;
                    k++;
                }
                if (k<theStrings.lLength) {
                    thisLine = (_String*)theStrings.lData[k];
                    if (v+stringSize.cy+2<bRC.bottom) { // enough room for a new line
                        currentDrawnLines++;
                        if (thisLine->sData[0]=='_') {
                            SelectObject (theContext, boldFont);
                            TextOut(theContext,h,v,thisLine->sData+1, thisLine->sLength-1);
                        } else {
                            SelectObject (theContext, plainFont);
                            TextOut(theContext,h,v,thisLine->sData, thisLine->sLength);
                        }
                    }
                } else {
                    currentFirstLine++;
                }
            } else {
                currentTopPosition = -1;
            }
        }

        EndPaint (theBox,&thePS);
        return 0;

    }
    case WM_TIMER: {
        //color -= 5;
        if (currentTopPosition == -1) {
            //KillTimer (DialogHandle,5);
            currentTopPosition = 15;
            currentFirstLine = 0;
            currentDrawnLines = 1;
            if (splash) {
                SendMessage (DialogHandle, WM_CLOSE, 0, 0);
            }
            return 1;
        }
        currentTopPosition++;
        HWND theBox = GetDlgItem (DialogHandle,129);
        RECT bRC;
        GetClientRect (theBox,&bRC);
        InvalidateRect (theBox, &bRC, FALSE);
        SendMessage (DialogHandle, WM_PAINT,0,0);
        return 1;

    }
    case WM_INITDIALOG: {
        char    buffer[256];
        _String theText;
        theStrings.Clear();
        splash = (bool)lParam;

        plainFont = CreateFont (8,0,0,0,400,FALSE,FALSE,FALSE,ANSI_CHARSET,OUT_DEFAULT_PRECIS,
                                CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,DEFAULT_PITCH|FF_SWISS,"MS Sans Serif"),
                    boldFont =  CreateFont (8,0,0,0,700,FALSE,TRUE,FALSE,ANSI_CHARSET,OUT_DEFAULT_PRECIS,
                                            CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,DEFAULT_PITCH|FF_SWISS,"MS Sans Serif");


        theStrings && " ";
        theText = GetVersionString();
        theStrings && & theText;
        k = 1;
        while (LoadString(ProgramInstance,k,buffer,255)) {
            if (splash&& buffer[0]==' ') {
                break;
            }
            theStrings && buffer;
            k++;
        }

        theStrings && " ";
        theStrings && "_System Tidbits";
        /*long      cpuCount = 1;

        PROCESSOR_POWER_INFORMATION * cpuInfo = new PROCESSOR_POWER_INFORMATION[cpuCount];
        checkPointer (cpuInfo);

        CallNtPowerInformation (ProcessorInformation,NULL,0,cpuInfo,sizeof (PROCESSOR_POWER_INFORMATION)*cpuCount);

        theText = _String ("CPU Count: ") & cpuCount;
        theStrings && & theText;

        theText = _String ("CPU Speed: ") & (long) cpuInfo[0].CurrentMHz & " MHz";
        theStrings && & theText;*/

        MEMORYSTATUS        memStat;
        GlobalMemoryStatus (&memStat);

        theText = _String("Total Memory:")&_String ((long)(memStat.dwTotalPhys>>10)) & " K";
        theStrings && & theText;

        theText = _String("Available Memory:")&_String ((long)(memStat.dwAvailPhys>>10)) & " K";
        theStrings && & theText;

        theText = MatrixExpCounter ();
        if (theText.sLength) {
            theStrings&&(&theText);
        }
        k += 5;
        currentTopPosition = 15;
        currentFirstLine = 0;
        currentDrawnLines = 1;
        SetTimer (DialogHandle, 5, splash?5:35, (TIMERPROC)aboutBoxHandlerFunction);
        centerDialogBox (DialogHandle);
        return 1;
        break;
    }
    }
    return 0;
}

//_________________________________________________________________________
void  displayAbout (bool splash)
{
    DialogBoxParam(  ProgramInstance, MAKEINTRESOURCE(103), (HWND)hyphyConsoleWindow->GetOSWindowData(),  aboutBoxHandlerFunction, (LPARAM)splash);
}

//_________________________________________________________________________
void  centerDialogBox (HWND hwndDlg)
{
    HWND hwndOwner;
    RECT rc, rcOwner,rcDlg;
    if ((hwndOwner = GetParent(hwndDlg)) == NULL) {
        hwndOwner = GetDesktopWindow();
    }

    GetWindowRect(hwndOwner, &rcOwner);
    GetWindowRect(hwndDlg, &rcDlg);
    CopyRect(&rc, &rcOwner);
    OffsetRect(&rcDlg, -rcDlg.left, -rcDlg.top);
    OffsetRect(&rc, -rc.left, -rc.top);
    OffsetRect(&rc, -rcDlg.right, -rcDlg.bottom);
    SetWindowPos(hwndDlg,HWND_TOP,
                 rcOwner.left + (rc.right / 2),
                 rcOwner.top + (rc.bottom / 2),
                 0, 0,
                 SWP_NOSIZE);
}

//_________________________________________________________________________

bool    handleGUI (bool )
{
    Sleep (10);
    MessageLoop ();
    return true;
}

//_________________________________________________________________________

int WINAPI WinMain (HINSTANCE thisInstance, HINSTANCE , LPSTR , int)
{
    WNDCLASSEX  windowClass;

    HICON       progIcon        = LoadIcon(thisInstance, MAKEINTRESOURCE(101));
    HCURSOR     arrowCursor     = LoadCursor(NULL, IDC_ARROW);

    ProgramInstance           = thisInstance;
    windowClass.cbSize        = sizeof(windowClass);/* Number of bytes in structure  */
    windowClass.cbClsExtra    = 0;              /* No extra class data            */
    windowClass.cbWndExtra    = 0;              /* No extra window data           */
    windowClass.hInstance     = thisInstance;   /* Current instance is the owner  */
    windowClass.hIcon         = progIcon;       /* Set icon                       */
    windowClass.hCursor       = arrowCursor;    /* Set cursor                     */
    windowClass.hIconSm       = progIcon;       /* Set icon                       */
    windowClass.lpszClassName = WinName;
    windowClass.lpszMenuName  = nil;
    windowClass.lpfnWndProc   = HYWndProc;
    windowClass.hbrBackground = NULL;
    windowClass.style         = CS_HREDRAW | CS_VREDRAW | CS_DBLCLKS ;
    windowClass.hCursor       = nil;

    bool regStatus            = RegisterClassEx(&windowClass);

    windowClass.lpszClassName = AppName;
    windowClass.lpszMenuName  = AppName;

    if (regStatus) {
        regStatus = RegisterClassEx(&windowClass);
    }

    if (!regStatus) {
        LPVOID lpMsgBuf;
        _String errMsg = _String ("Failed to register a window class: ") &windowClass.lpszClassName & " with error: ";

        if (!FormatMessage( FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
                            NULL,GetLastError(),MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPTSTR) &lpMsgBuf,0,NULL )) {
            errMsg = errMsg&"unknown system error";
        } else {
            errMsg = errMsg & (char*)lpMsgBuf;
        }
        MessageBox(nil, errMsg.getStr() , GetVersionString(), MB_OK | MB_ICONEXCLAMATION);
    }


    char    dlgName[] = "HPGUIDLG";

    windowClass.lpszMenuName  = empty.sData;
    windowClass.lpszClassName = dlgName;
    windowClass.style        |= CS_NOCLOSE;
    windowClass.lpszMenuName  = nil;
    windowClass.lpfnWndProc   = HYDlgProc;
    windowClass.hbrBackground = (HBRUSH)(COLOR_MENU+1);
    windowClass.cbWndExtra    = DLGWINDOWEXTRA;

    if (!RegisterClassEx(&windowClass)) {
        _String errMsg = _String ("Failed to register a window class: ") &dlgName& " with error " & (long)GetLastError();
        MessageBox(nil, errMsg.getStr() , GetVersionString(), MB_OK | MB_ICONEXCLAMATION);
    }

    HBITMAP        sbPatternBitMap = LoadBitmap (ProgramInstance,MAKEINTRESOURCE(3000));
    windowStatusBarBrush = CreatePatternBrush (sbPatternBitMap);

    if (!LoadLibrary ("Riched20.dll")) {
        MessageBox(nil, "HyPhy reqires Rich Edit 2.0 or better (included in Riched20.dll)", GetVersionString(), MB_OK | MB_ICONEXCLAMATION);
        exit (0);
    }

    //char curWd[4096];

    //_getcwd (curWd,4096);
    
    HMODULE hModule = GetModuleHandleW(NULL);
    CHAR path[MAX_PATH];
    GetModuleFileName(hModule, path, MAX_PATH);

    _String baseDir = path;
    
    baseDir = baseDir.Cut(0,baseDir.FindBackwards ("\\",0,-1));
    
    if (baseDir.sData[baseDir.sLength-1]!='\\') {
        baseDir = baseDir & '\\';
    }

    baseDirectory = baseDir;
    libDirectory  = baseDir;
    pathNames&&     &baseDir;

    statusBarPlain = CreateFont (9,0,0,0,400,FALSE,FALSE,FALSE,ANSI_CHARSET,OUT_DEFAULT_PRECIS,
                                 CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,DEFAULT_PITCH|FF_SWISS,"MS Sans Serif");
    statusBarBold =  CreateFont (9,0,0,0,700,FALSE,FALSE,FALSE,ANSI_CHARSET,OUT_DEFAULT_PRECIS,
                                 CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,DEFAULT_PITCH|FF_SWISS,"MS Sans Serif");

    InitCommonControls();
    INITCOMMONCONTROLSEX icces;
    icces.dwICC = ICC_WIN95_CLASSES;
    icces.dwSize = sizeof (INITCOMMONCONTROLSEX);
    InitCommonControlsEx(&icces);


    LoadHYPHYResources ();

    hyphyConsoleWindow = new _HYConsoleWindow (consoleWindowTitle);

    UINT_PTR gqTimer = SetTimer ((HWND)hyphyConsoleWindow->GetOSWindowData(), GLOBAL_QUEUE_TIMER , 100,(TIMERPROC) GlobalQueueTimer);
    WinInitSet();
    ReadChartProcessors ();
    ReadPreferences     ();
    MoveConsoleWindow   (consolePositionRectangle);
    hyphyConsoleWindow->BringToFront();
    

    if (showDialogAtStartup) {
        long                  choice = 0;

        showDialogAtStartup = displayFirstDialog (choice);

        if (!showDialogAtStartup) {
            _String wStr ("You can later change startup dialog settings in 'Preferences'");
            ProblemReport (wStr);
        }

        if (choice == 2) {
            RunStandardAnalyses();
        } else {
            _String dL (baseDirectory);
            if (choice == 3) {
                dL = dL & "data\\";
                OpenDataFile (&dL);
            } else if (choice == 4) {
                dL = dL & "Saves\\";
                if (OpenBatchFile (true,&dL)) {
                    ExecuteBatchFile ();
                }
            }
        }

        SetShowDialogAtStartup (showDialogAtStartup);
    }


    while (!hyphyExiting) {
        //if (runFile())
        //  ExecuteBatchFile ();
        //else
        MessageLoop(false);
    }

    // -- end main -- //

    WritePreferences    ();
    GlobalShutdown      ();
    KillTimer ((HWND)hyphyConsoleWindow->GetOSWindowData(),gqTimer);
    //exit(0);
    DeleteObject (statusBarPlain);
    DeleteObject (statusBarBold);
    DeleteDC (offScreenContext);
    return   0;

}
