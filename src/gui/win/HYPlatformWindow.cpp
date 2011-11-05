/*
    A Windows window object - a window/title/size-box/scroll-bars handler.

    Sergei L. Kosakovsky Pond, June 2000.
*/

#include "HYWindow.h"
#include "HYGWindow.h"
#include "HYPlatformComponent.h"
#include "HYEventTypes.h"
#include "HYTableWindow.h"
#include "HYTreePanel.h"
#include "HYUtils.h"
#include "HYConsoleWindow.h"


extern  _SimpleList             windowPtrs,
        windowObjects,
        windowObjectRefs;

extern  _String                 consoleWindowTitle;

/* printing stuff */
extern                          FILE*       dbgfile;

BOOL CALLBACK                   PrintDialogProc(HWND, UINT, WPARAM, LPARAM);
BOOL CALLBACK                   AbortProc(HDC, int);
HDC                             GetPrinterDeviceContext(HWND);

extern  BOOL                    UserAbortFlag;
extern  HWND                    PrintDialogHandle;

extern                          HFONT           statusBarPlain;

//__________________________________________________________________

LRESULT coreProcessor (HWND WindowHand, UINT iMsg, WPARAM wParam, LPARAM lParam)
{
    _HYWindow* me = nil;

    long       f  = windowPtrs.Find((long)WindowHand);

    if (f>=0) {
        me = (_HYWindow*)(windowObjects.lData[f]);
    }

    switch (iMsg) {
    case WM_COMMAND:
        if (!lParam) {
            if (me->_ProcessMenuSelection(LOWORD(wParam))) {
                return true;
            }
        }

        _HYWindowsUIMessage thisM;
        thisM.iMsg = iMsg;
        thisM.wParam = wParam;
        thisM.lParam = lParam;
        if (me) {
            if (((_HYPlatformWindow*)me)->_ProcessOSEvent ((Ptr)&thisM)) {
                return true;
            }
        }
        return false;

    case WM_CREATE:

        return true ;

    case WM_SIZE:

        me->Grow(nil);
        return true;

    case WM_PAINT: {
        if (me) {
            me->Update(nil);
        }

        return true;
    }

    case WM_CLOSE: {
        //printf ("WM_CLOSE\n");
        /*f = windowObjects.Find((long)(me));
        if (f>=0)
        {
            windowObjectRefs.Delete(f);
            windowObjects.Delete(f);
            windowPtrs.Delete(f);
        }           */
        if (((_HYWindow*)windowObjectRefs(f))->Close(nil)) {
            DestroyWindow (WindowHand);
            return -2;
        }
        return false;
    }

    case WM_QUERYENDSESSION :
        return true;

    case WM_DESTROY: {
        //printf ("WM_DESTOY\n");
        //printf ("_Close handler %d\n",f);
        if (f>=0) {
            //DestroyWindow ((HWND)windowPtrs(f));
            windowObjects.Delete(f);
            windowPtrs.Delete(f);
            windowObjectRefs.Delete(f);
        }
        return true;
    }

    default: {

        //if ((iMsg == WM_ACTIVATE)&&(LOWORD(wParam) != WA_INACTIVE))
        //me->_SetMenuBar ();

        _HYWindowsUIMessage thisM;
        thisM.iMsg = iMsg;
        thisM.wParam = wParam;
        thisM.lParam = lParam;
        thisM.res    = nil;
        if (me) {
            if (((_HYPlatformWindow*)me)->_ProcessOSEvent ((Ptr)&thisM))
                if (thisM.res) {
                    return thisM.res;
                } else {
                    return 1L;
                }
        }
    }

    }
    return false;
}

//__________________________________________________________________

LRESULT CALLBACK HYWndProc (HWND WindowHand, UINT iMsg, WPARAM wParam, LPARAM lParam)
{
    LRESULT res = coreProcessor (WindowHand, iMsg, wParam, lParam);
    if (res) {
        return res==-2?0:res;
    }

    /*if (iMsg == WM_ACTIVATE)
        if (LOWORD (wParam) != WA_INACTIVE)
        {
            long f = windowPtrs.Find (lParam);
            if (f>=0)
            {
                _HYWindow * deact = (_HYWindow*)windowObjectRefs (f);
                if (deact->flags & HY_WINDOW_DLOG)
                {
                    if ((deact->parentWindowPtr==nil)||(deact->parentWindowPtr==(Ptr)WindowHand))
                    {
                        SetActiveWindow ((HWND)windowPtrs(f));
                        ReportWarning ("SetActiveWindow");
                        return 0;
                    }
                }
            }
        }*/
    return DefWindowProc(WindowHand, iMsg, wParam, lParam);
}

//__________________________________________________________________

LRESULT CALLBACK HYDlgProc (HWND WindowHand, UINT iMsg, WPARAM wParam, LPARAM lParam)
{

    if (iMsg == WM_CREATE) {
        HWND pw = (HWND)((LPCREATESTRUCT)lParam)->lpCreateParams;
        if (pw) {
            EnableWindow (pw,false);
        } else {
            EnableWindow ((HWND)hyphyConsoleWindow->GetOSWindowData(),false);
            for (long k = 0; k<windowPtrs.lLength; k++) {
                EnableWindow ((HWND)windowPtrs(k),false);
            }
        }
    } else {
        if ((iMsg == WM_CLOSE)||(iMsg == WM_DESTROY)) {
            long f = windowPtrs.Find ((long)WindowHand);
            if (f>=0) {
                _HYWindow * pWindow = (_HYWindow*)windowObjectRefs (f);
                if (pWindow->parentWindowPtr) {
                    EnableWindow ((HWND)pWindow->parentWindowPtr,true);
                } else {
                    EnableWindow ((HWND)hyphyConsoleWindow->GetOSWindowData(),true);
                    for (long k = 0; k<windowPtrs.lLength; k++) {
                        EnableWindow ((HWND)windowPtrs(k),true);
                    }
                }
                SetActiveWindow (pWindow->oldFrontWindow);
            }
        }
    }

    LRESULT res = coreProcessor (WindowHand, iMsg, wParam, lParam);

    if (res!=0) {
        return res==-2?0L:res;
    }

    return DefWindowProc(WindowHand, iMsg, wParam, lParam);
}

//__________________________________________________________________

_HYPlatformWindow::_HYPlatformWindow(unsigned char windowFlag,_String windowTitle,bool windowVisibility,Ptr par)
{
    DWORD dwStyle   = WS_CAPTION;
    DWORD exStyle   = 0;
    parentWindowPtr = par;
    oldFrontWindow  = nil;

    bool  isDlg = (windowFlag & HY_WINDOW_DLOG)||par,
          isConsole = (windowTitle.Equal (&consoleWindowTitle));

    if (isDlg) {
        windowMenu      = CreateMenu();
        exStyle         = WS_EX_DLGMODALFRAME;
        oldFrontWindow  = GetActiveWindow();
    } else {
        if (isConsole) {
            windowMenu = nil;
            dwStyle |= WS_SYSMENU;
        } else {
            windowMenu = CreateMenu();
            dwStyle |= WS_SYSMENU;
            HMENU        fileHandle = CreateMenu ();
            checkPointer (fileHandle);
            InsertMenu   (windowMenu, 0, MF_BYPOSITION|MF_POPUP, (UINT) fileHandle , "&File");
        }
    }

    if (windowFlag&HY_WINDOW_SIZE) {
        if (windowFlag&HY_WINDOW_ZOOM) {
            dwStyle = WS_OVERLAPPEDWINDOW;
        } else {
            dwStyle |= WS_MINIMIZEBOX|WS_THICKFRAME;
        }

        if (isConsole) {
            dwStyle |= WS_MAXIMIZEBOX;
        }

    } else {
        if (isDlg)
            ;
        else {
            dwStyle |= WS_MINIMIZEBOX|WS_MAXIMIZEBOX;
        }
    }

    menuKeys = NULL;

    HWND        pWnd = NULL;

    if (par) {
        pWnd = (((_HYWindow*)par)->theWindow);
        parentWindowPtr = (Ptr)pWnd;
    }

    theWindow = CreateWindowEx(exStyle,isConsole?"HYPHY":(isDlg?"HPGUIDLG":"HPGUIW"),           /* window class name              */
                               windowTitle.sData,       /* window caption                 */
                               dwStyle | (((windowFlag&HY_WINDOW_SCROLL)&&(windowFlag&HY_WINDOW_SIZE))?(WS_HSCROLL|WS_VSCROLL):0)| WS_CLIPCHILDREN ,
                               CW_USEDEFAULT,          /* initial x position             */
                               CW_USEDEFAULT,          /* initial y position             */
                               300,                     /* initial x size                 */
                               250,                     /* initial y size                 */
                               pWnd,                   /* parent window handle           */
                               windowMenu,             /* window menu handle             */
                               ProgramInstance,         /* program instance handle        */
                               par?pWnd:nil);                   /* creation parameters            */


    if ((!theWindow)||(!(windowMenu)&&(!isConsole))) {
        _String errMsg ("Could not create a window structure. Windows Error:");
        errMsg = errMsg&_String ((long)GetLastError());
        FlagError (errMsg);
    }


    if ((windowFlag&HY_WINDOW_SCROLL)&&(windowFlag&HY_WINDOW_SIZE)) {
        SetScrollRange (theWindow,SB_HORZ,0,MAX_CONTROL_VALUE,FALSE);
        SetScrollRange (theWindow,SB_VERT,0,MAX_CONTROL_VALUE,FALSE);
    }


    windowPtrs<<(long)theWindow;
    windowObjects<<(long)(_HYWindow*)this;



    flags = windowFlag;

}

//__________________________________________________________________

_HYPlatformWindow::~_HYPlatformWindow(void)
{
    if (windowMenu) {
        DestroyMenu (windowMenu);
    }
    if (menuKeys) {
        DestroyAcceleratorTable (menuKeys);
    }
    DestroyWindow (theWindow);
}

//__________________________________________________________________

void _HYPlatformWindow::_SetTitle(_String windowTitle)
{
    SetWindowText (theWindow,windowTitle);
}

//__________________________________________________________________

void _HYPlatformWindow::_Show(void)
{
    ShowWindow (theWindow,SW_SHOWNORMAL);
}

//__________________________________________________________________

void _HYPlatformWindow::_Hide(void)
{
    ShowWindow (theWindow,SW_HIDE);
}

//__________________________________________________________________

_HYRect _HYPlatformWindow::_GetWindowRect(void)
{
    RECT res;

    ::GetWindowRect (theWindow, &res);

    return (_HYRect) {
        res.top,res.left,res.bottom,res.right,0
    };
}

//__________________________________________________________________

long _HYPlatformWindow::_Grow(Ptr theData)
{
    RECT myDims;
    GetClientRect (theWindow,&myDims);
    _HYWindow* theParent = (_HYWindow*)this;
    theParent->right = theParent->left+myDims.right-myDims.left;
    theParent->bottom = theParent->top+myDims.bottom-myDims.top;
    int windowWidth = theParent->right-theParent->left,
        windowHeight = theParent->bottom-theParent->top;
    if (flags&HY_WINDOW_SIZE) {
        _HYWindow* theParent = (_HYWindow*)this;
        if (windowWidth >= theParent->contentWidth) {
            EnableScrollBar (theWindow,SB_HORZ, ESB_DISABLE_BOTH);
        } else {
            EnableScrollBar (theWindow,SB_HORZ, ESB_ENABLE_BOTH);
        }
        if (windowHeight >= theParent->contentHeight) {
            EnableScrollBar (theWindow,SB_VERT, ESB_DISABLE_BOTH);
        } else {
            EnableScrollBar (theWindow,SB_VERT, ESB_ENABLE_BOTH);
        }
    }
    return 0;
}

//__________________________________________________________________

_String& _HYPlatformWindow::_GetTitle(void)
{
    return      ((_HYWindow*)this)->GetTitle();
}

//__________________________________________________________________

void    _HYPlatformWindow::_SetWindowBackColor (_HYColor newColor)
{
    PAINTSTRUCT      thePS;
    HDC theContext = BeginPaint (theWindow,&thePS);

    SetBkColor  (theContext, HYColor2ColorRef (newColor, theContext));

    EndPaint (theWindow,&thePS);
}

//__________________________________________________________________

void _HYPlatformWindow::_Move(Ptr theData)
{
}

//__________________________________________________________________

void _HYPlatformWindow::_SetPosition(int l,int t)
{
    SetWindowPos (theWindow,NULL,l,t,0,0,SWP_NOZORDER|SWP_NOSIZE);
}

//__________________________________________________________________

void _HYPlatformWindow::_BringWindowToFront(void)
{
    //RedrawWindow (theWindow, NULL, NULL, RDW_ALLCHILDREN);
    //ShowWindow (theWindow, SW_SHOWNORMAL);
    _Activate ();
}

//__________________________________________________________________

bool _HYPlatformWindow::_Close(Ptr theData)
{
    _HYWindow* theParent = (_HYWindow*)this;
    bool        doit = true;

    if (theParent->ConfirmClose()) {
        //long f = windowObjectRefs.Find((long)theParent);
        //printf ("_Close handler %d\n",f);
        /*if (f>=0)
        {
            DestroyWindow ((HWND)windowPtrs(f));
            windowObjects.Delete(f);
            windowPtrs.Delete(f);
            windowObjectRefs.Delete(f);
        }*/
    } else {
        doit = false;
    }

    return doit;
}

//__________________________________________________________________

void _HYPlatformWindow::_Activate(void)
{
    _SetMenuBar();
    _Show();
    SetFocus(theWindow);
}

//__________________________________________________________________

void _HYPlatformWindow::_Deactivate(void)
{
    //_Hide();
}

//__________________________________________________________________

void _HYPlatformWindow::_AddStandardAccels(void)
{
    accels       << (FCONTROL|FVIRTKEY);
    accels       << 'Z';
    accels       << HY_WINDOW_MENU_ID_EDIT;

    accels       << (FCONTROL|FVIRTKEY);
    accels       << 'C';
    accels       << HY_WINDOW_MENU_ID_EDIT+1;

    accels       << (FCONTROL|FVIRTKEY);
    accels       << 'X';
    accels       << HY_WINDOW_MENU_ID_EDIT+2;

    accels       << (FCONTROL|FVIRTKEY);
    accels       << 'V';
    accels       << HY_WINDOW_MENU_ID_EDIT+3;

    accels       << (FCONTROL|FVIRTKEY);
    accels       << 'A';
    accels       << HY_WINDOW_MENU_ID_EDIT+5;

    accels       << (FCONTROL|FVIRTKEY);
    accels       << 'W';
    accels       << HY_WINDOW_MENU_ID_FILE;

    accels       << (FCONTROL|FVIRTKEY);
    accels       << 'S';
    accels       << HY_WINDOW_MENU_ID_FILE+1;

    accels       << (FCONTROL|FVIRTKEY);
    accels       << 'P';
    accels       << HY_WINDOW_MENU_ID_FILE+2;

    accels       << (FCONTROL|FVIRTKEY);
    accels       << '0';
    accels       << HY_WINDOW_MENU_ID_FILE-1;

    accels       << (FCONTROL|FVIRTKEY);
    accels       << 'I';
    accels       << HY_WINDOW_MENU_ID_FILE-2;
}

//__________________________________________________________________

void _HYPlatformWindow::_BuildAccelTable (bool force)
{
    if (force&&menuKeys) {
        DestroyAcceleratorTable (menuKeys);
        menuKeys = nil;
    }

    if ((menuKeys == NULL)&&(accels.lLength>=3)) {
        LPACCEL myAccels = (LPACCEL) new ACCEL[accels.lLength/3];
        checkPointer (myAccels);

        for (long k=0; k<accels.lLength/3; k++) {
            myAccels[k].fVirt = accels.lData[k*3];
            myAccels[k].key = accels.lData[k*3+1];
            myAccels[k].cmd = accels.lData[k*3+2];
        }

        menuKeys = CreateAcceleratorTable (myAccels, accels.lLength/3);
        delete (myAccels);
    }
}

//__________________________________________________________________

void _HYPlatformWindow::_SetMenuBar(void)
{
    HMENU        fileHandle = GetSubMenu(GetMenu (theWindow),0);

    if (windowMenu&&(GetMenuItemCount (fileHandle)==0)) {

        _HYWindow* theParent = (_HYWindow*)this;

        HMENU        editHandle = CreateMenu ();

        checkPointer (fileHandle);
        checkPointer (editHandle);

        InsertMenu   (fileHandle, 0xFFFFFFFF, MF_BYPOSITION, HY_WINDOW_MENU_ID_FILE+1, "&Save\tCtrl-S");
        InsertMenu   (fileHandle, 0xFFFFFFFF, MF_BYPOSITION, HY_WINDOW_MENU_ID_FILE+2, "&Print\tCtrl-P");
        InsertMenu   (fileHandle, 0xFFFFFFFF, MF_BYPOSITION, HY_WINDOW_MENU_ID_FILE,   "&Close\tCtrl-W");
        InsertMenu   (fileHandle, 0xFFFFFFFF, MF_BYPOSITION, HY_WINDOW_MENU_ID_FILE-1, "S&witch to Console\tCtrl-0");
        InsertMenu   (fileHandle, 0xFFFFFFFF, MF_BYPOSITION, HY_WINDOW_MENU_ID_FILE-2, "Object &Inspector\tCtrl-I");

        if (!theParent->IsSaveEnabled()) {
            EnableMenuItem (fileHandle, 0, MF_BYPOSITION|MF_GRAYED);
        }

        if (!theParent->IsPrintEnabled()) {
            EnableMenuItem (fileHandle, 1, MF_BYPOSITION|MF_GRAYED);
        }

        InsertMenu   (editHandle, 0xFFFFFFFF, MF_BYPOSITION|MF_GRAYED, HY_WINDOW_MENU_ID_EDIT, "&Undo\tCtrl-Z");

        InsertMenu   (editHandle, 0xFFFFFFFF, MF_BYPOSITION|MF_SEPARATOR, 0, NULL);

        InsertMenu   (editHandle, 0xFFFFFFFF, MF_BYPOSITION|MF_GRAYED, HY_WINDOW_MENU_ID_EDIT+1, "&Copy\tCtrl-C");

        InsertMenu   (editHandle, 0xFFFFFFFF, MF_BYPOSITION|MF_GRAYED, HY_WINDOW_MENU_ID_EDIT+2, "C&ut\tCtrl-X");

        InsertMenu   (editHandle, 0xFFFFFFFF, MF_BYPOSITION|MF_GRAYED, HY_WINDOW_MENU_ID_EDIT+3, "&Paste\tCtrl-V");

        InsertMenu   (editHandle, 0xFFFFFFFF, MF_BYPOSITION|MF_GRAYED, HY_WINDOW_MENU_ID_EDIT+4, "C&lear");

        InsertMenu   (editHandle, 0xFFFFFFFF, MF_BYPOSITION|MF_GRAYED, HY_WINDOW_MENU_ID_EDIT+5, "Select &All\tCtrl-A");

        //InsertMenu     (windowMenu, 0, MF_BYPOSITION|MF_POPUP, (UINT) fileHandle , "&File");
        InsertMenu   (windowMenu, 1, MF_BYPOSITION|MF_POPUP, (UINT) editHandle , "&Edit");
    }

    _AddStandardAccels();
    _BuildAccelTable  ();

    accels.Clear();
}

//__________________________________________________________________

void _HYPlatformWindow::_UnsetMenuBar(void)
{
}

//__________________________________________________________________

void _HYPlatformWindow::_Paint(Ptr)
{
    PAINTSTRUCT      thePS;
    HDC theContext = BeginPaint (theWindow,&thePS);
    EndPaint (theWindow,&thePS);
}

//__________________________________________________________________

void _HYPlatformWindow::_Update(Ptr)
{
}

//__________________________________________________________________

void _HYPlatformWindow::_SetWindowRectangle(int top, int left, int bottom, int right, bool ss)
{

    int windowWidth  = right-left,
        windowHeight = bottom-top;

    if (flags&HY_WINDOW_SIZE) {
        _HYWindow* theParent = (_HYWindow*)this;

        if (windowWidth >= theParent->contentWidth) {
            EnableScrollBar (theWindow,SB_HORZ, ESB_DISABLE_BOTH);
        } else {
            EnableScrollBar (theWindow,SB_HORZ, ESB_ENABLE_BOTH);
        }
        if (windowHeight >= theParent->contentHeight) {
            EnableScrollBar (theWindow,SB_VERT, ESB_DISABLE_BOTH);
        } else {
            EnableScrollBar (theWindow,SB_VERT, ESB_ENABLE_BOTH);
        }
    }

    if (ss) {
        RECT            smallR,
                        bigR;

        GetWindowRect   (theWindow,&bigR);
        GetClientRect   (theWindow,&smallR);
        SetWindowPos    (theWindow,NULL,0,0,
                         windowWidth + bigR.right  - bigR.left - smallR.right,
                         windowHeight+ bigR.bottom - bigR.top  - smallR.bottom,
                         SWP_NOZORDER|SWP_NOMOVE);

        InvalidateRect  (theWindow, nil, false);
    }
}

//__________________________________________________________________

bool _HYPlatformWindow::_ProcessOSEvent (Ptr vEvent)
{
    _HYWindowsUIMessage * theEvent = (_HYWindowsUIMessage *)vEvent;
    _HYWindow* theParent = (_HYWindow*)this;
    long       invisPixels;
    long       currentValue, oldValue;
    switch (theEvent->iMsg) {
    case WM_HSCROLL:
    case WM_VSCROLL: {
        if (theEvent->iMsg == WM_HSCROLL) {
            invisPixels = theParent->contentWidth-(theParent->right-theParent->left);
            currentValue = GetScrollPos (theWindow,SB_HORZ);
        } else {
            invisPixels = theParent->contentHeight-(theParent->bottom-theParent->top);
            currentValue = GetScrollPos (theWindow,SB_VERT);
        }
        oldValue = currentValue;
        invisPixels = (_Parameter)MAX_CONTROL_VALUE/invisPixels;
        if (!invisPixels) {
            invisPixels = 1;
        }

        switch (LOWORD(theEvent->wParam)) {
        case SB_LINEUP:
            currentValue-=invisPixels;
            break;
        case SB_LINEDOWN:
            currentValue+=invisPixels;
            break;
        case SB_PAGEUP:
            currentValue-=10*invisPixels;
            break;
        case SB_PAGEDOWN:
            currentValue+=10*invisPixels;
            break;
        case SB_THUMBPOSITION:
        case SB_THUMBTRACK:
            currentValue = HIWORD(theEvent->wParam);
            break;
        case SB_TOP:
            currentValue = 0;
            break;
        case  SB_BOTTOM:
            currentValue = MAX_CONTROL_VALUE;
            break;
        default:
            return true;
        }

        if (currentValue<0) {
            currentValue = 0;
        }
        if (currentValue>MAX_CONTROL_VALUE) {
            currentValue = MAX_CONTROL_VALUE;
        }
        if (currentValue!=oldValue) {
            if (theEvent->iMsg == WM_HSCROLL) {
                SetScrollPos (theWindow,SB_HORZ, currentValue, TRUE);
                theParent->ProcessEvent (generateScrollEvent (currentValue-oldValue,0));
            } else {
                SetScrollPos (theWindow,SB_VERT, currentValue, TRUE);
                theParent->ProcessEvent (generateScrollEvent (0,currentValue-oldValue));
            }
        }
    }
    return true;
    break;

    case WM_GETMINMAXINFO:

        MINMAXINFO* windowInfo = (MINMAXINFO*)theEvent->lParam;
        RECT   bigR, smallR;
        GetClientRect (theWindow,&smallR);
        GetWindowRect (theWindow,&bigR);
        windowInfo->ptMinTrackSize.x = 50;
        windowInfo->ptMinTrackSize.y = 50;
        windowInfo->ptMaxTrackSize.x = theParent->contentWidth+bigR.right-bigR.left-(smallR.right);
        windowInfo->ptMaxTrackSize.y = theParent->contentHeight+bigR.bottom-bigR.top-smallR.bottom;
        return true;
        break;
    }
    return false;
}

//__________________________________________________________________

bool        _HYWindow::_ProcessMenuSelection (long msel)
{
    switch (msel) {
    case HY_WINDOW_MENU_ID_FILE-2: { // close window
        ShowObjectInspector();
        return true;
    }
    case HY_WINDOW_MENU_ID_FILE-1: { // close window
        SetForegroundWindow ((HWND)hyphyConsoleWindow->GetOSWindowData());
        return true;
    }
    case HY_WINDOW_MENU_ID_FILE: { // close window
        Close(nil);
        return true;
    }
    }
    return false;
}

//__________________________________________________________________

void    _HYPlatformWindow::_SetContentSize (int w, int h)
{
    _HYWindow* theParent = (_HYWindow*)this;
    int     windowW = theParent->right-theParent->left,
            windowH =  theParent->bottom-theParent->top;
    if (w<=windowW) {
        EnableScrollBar (theWindow,SB_HORZ,ESB_DISABLE_BOTH);
    } else {
        EnableScrollBar (theWindow,SB_HORZ,ESB_ENABLE_BOTH);
        SetScrollPos (theWindow,SB_HORZ,0,TRUE);
    }
    if (h<=windowH) {
        EnableScrollBar (theWindow,SB_VERT,ESB_DISABLE_BOTH);
    } else {
        EnableScrollBar (theWindow,SB_VERT,ESB_ENABLE_BOTH);
        SetScrollPos (theWindow,SB_VERT,0,TRUE);
    }

}

//__________________________________________________________________


void    _HYPlatformWindow::_VisibleContents (int& t,int& l,int& b,int& r)
{
    _HYWindow* theParent = (_HYWindow*)this;
    short v;
    short windowH=theParent->bottom-theParent->top,
          windowW=theParent->right-theParent->left;
    if (flags&HY_WINDOW_SIZE) {
        if (windowW>theParent->contentWidth) {
            l = theParent->left;
            r = theParent->right;
        } else {
            v = GetScrollPos (theWindow,SB_HORZ);
            l = theParent->left+(theParent->contentWidth-windowW)*v/(_Parameter)MAX_CONTROL_VALUE;
            r = l+windowW;
        }
        if (windowH>theParent->contentHeight) {
            t = theParent->top;
            b = theParent->bottom;
        } else {
            v = GetScrollPos (theWindow,SB_VERT);
            t = theParent->top+(theParent->contentHeight-windowH)*v/(_Parameter)MAX_CONTROL_VALUE;
            b = t+windowH;
        }
    } else {
        t = theParent->top;
        l = theParent->left;
        b = theParent->bottom;
        r = theParent->right;
    }

}


//__________________________________________________________________

RECT    _HYPlatformWindow::newVRect (void)
{
    return (RECT) {
        0,0,0,0
    };
}

//__________________________________________________________________

RECT    _HYPlatformWindow::newHRect (void)
{
    return (RECT) {
        0,0,0,0
    };
}

//__________________________________________________________________

void    _HYPlatformWindow::_SelectWindow (void)
{
}

//__________________________________________________________________
void        _PaintTheCircle (HBITMAP aPic, HWND theWindow, HDC theDC)
{
    RECT r;
    GetClientRect (theWindow, &r);

    r.right = r.left+15;
    r.left  += 3;
    r.bottom --;
    r.top = r.bottom-12;

    DrawTransparentBitmap (theDC, aPic, r.left, r.top, 12, 12, RGB(255,255,255));

    /*BITMAP theBM;
    theBM.bmBits = nil;



    GetObject (aPic, sizeof (BITMAP), &theBM);
    if (otherDC)
    {
        if (r.right-r.left<=0)
            r.right  = r.left + theBM.bmWidth;

        if (r.bottom-r.top<=0)
            r.bottom  = r.top + theBM.bmHeight;

        SelectObject (otherDC, aPic);

        TransparentBlt   (theDC, r.left, r.top, r.right-r.left+1, r.bottom-r.top+1,
                            otherDC, 0, 0, theBM.bmWidth, theBM.bmHeight, RGB(255,255,255));

        SelectObject (otherDC, oDCBM);
    }*/
}

//EOF