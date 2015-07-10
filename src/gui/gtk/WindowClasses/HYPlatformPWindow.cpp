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


/*
    TBI parts of HYPWindow for GTK
*/

#include "HYGWindow.h"
#include "HYCanvas.h"
#include "errorfns.h"
#include "HYDialogs.h"

#define  HY_PWINDOW_WIN32_MENU_BASE  9000


//__________________________________________________________________
// _HYPlatformPWindow
//__________________________________________________________________

_HYPlatformPWindow::_HYPlatformPWindow          (void)
{
    // TBI
}

//__________________________________________________________________

_HYPlatformPWindow::~_HYPlatformPWindow         (void)
{
    // TBI
}

//__________________________________________________________________

void    _HYPlatformPWindow::_StartPicture   (void)
{
    // TBI
}

//__________________________________________________________________

void    _HYPlatformPWindow::_EndPicture     (void)
{
    // TBI
}

//__________________________________________________________________

void    _HYPlatformPWindow::_DrawPicture    (_HYRect r)
{
    //TBI
}

//__________________________________________________________________
// _HYPWindow
//__________________________________________________________________

void        _HYPWindow::_PrintPWindow(void)
{
    // TBI
}

//__________________________________________________________________

void _HYPWindow::_SetWindowRectangle(int top, int left, int bottom, int right, bool ss)
{
    if (theWindow && ss && GTK_WIDGET_MAPPED (theWindow)) {
        long    menuHeight = 0;
        if (windowMB) {
            menuHeight = windowMB->allocation.x;
            if (menuHeight <= 0) {
                GtkRequisition sizeReq;
                gtk_widget_size_request (windowMB,&sizeReq);
                menuHeight = sizeReq.height;
            }
        }

        GdkGeometry windowG;
        windowG.min_width = windowG.max_width  = right-left;
        windowG.max_height = bottom-top+menuHeight;
        gtk_window_set_geometry_hints (GTK_WINDOW(theWindow), NULL, &windowG, (GdkWindowHints)(GDK_HINT_MAX_SIZE));

        //SetContentSize    (right-left,bottom-top);
        _HYPlatformWindow::_SetWindowRectangle (top,left,bottom,right, ss);
    }
}

//__________________________________________________________________

void    _HYPWindow::_Paint (Ptr p)
{
    // TBI
    _HYGWindow::_Paint (p);
}

//__________________________________________________________________

void    _HYPWindow::_Update (Ptr p)
{
    _Paint (p);
}

//__________________________________________________________________

long _HYPWindow::_Grow(Ptr theData)
{
    /*_HYPlatformWindow::_Grow (theData);
    RECT myDims;
    GetClientRect (theWindow,&myDims);

    SetWindowRectangle (0,0,myDims.bottom,myDims.right,false);

    GetClientRect  (theWindow,&myDims);
    InvalidateRect (theWindow,&myDims,false);*/
    return 0;
}

//__________________________________________________________________

bool        _HYPWindow::_ProcessMenuSelection (long msel)
{
    return _HYWindow::_ProcessMenuSelection(msel);

    /*bool      done = false;

    switch (msel)
    {
        case HY_WINDOW_MENU_ID_FILE+2:
            _PrintPWindow();
            done = true;
            break;

        case HY_PWINDOW_WIN32_MENU_BASE:
            Zoom (1.1);
            done = true;
            break;

        case HY_PWINDOW_WIN32_MENU_BASE+1:
            Zoom (.9);
            done = true;
            break;

        case HY_PWINDOW_WIN32_MENU_BASE+2:
            OriginalSize();
            done = true;
            break;
    }

    if (!done)
        return _HYGWindow::_ProcessMenuSelection (msel);

    return done;*/
}

//__________________________________________________________________

void _HYPWindow::_SetMenuBar(void)
{
    /*
        _HYGWindow::_SetMenuBar();

        HMENU            windowMenu = GetMenu (theWindow),
                         chartMenu =  GetSubMenu(windowMenu,2);

        if (!chartMenu)
        {
            chartMenu   = CreateMenu();

            InsertMenu   (chartMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_PWINDOW_WIN32_MENU_BASE  , "&Enlarge\tCtrl-1");
            InsertMenu   (chartMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_PWINDOW_WIN32_MENU_BASE+1, "&Shrink\tCtrl-2");
            InsertMenu   (chartMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_PWINDOW_WIN32_MENU_BASE+2, "&Original Size\tCtrl-3");

            InsertMenu   (windowMenu, 2, MF_BYPOSITION|MF_POPUP, (UINT) chartMenu , "&Image");


            accels       << (FCONTROL|FVIRTKEY);
            accels       << '1';
            accels       << HY_PWINDOW_WIN32_MENU_BASE;

            accels       << (FCONTROL|FVIRTKEY);
            accels       << '2';
            accels       << HY_PWINDOW_WIN32_MENU_BASE+1;

            accels       << (FCONTROL|FVIRTKEY);
            accels       << '3';
            accels       << HY_PWINDOW_WIN32_MENU_BASE+2;

            _AddStandardAccels();
            _BuildAccelTable  (true);
            accels.Clear();
        }

        DrawMenuBar(theWindow);*/

}

//__________________________________________________________________

void _HYPWindow::_UnsetMenuBar(void)
{

}

//__________________________________________________________________

bool _HYPWindow::_ProcessOSEvent (Ptr vEvent)
{
    /*_HYWindowsUIMessage * theEvent = (_HYWindowsUIMessage *)vEvent;

    switch (theEvent->iMsg)
    {
        case WM_GETMINMAXINFO:

            MINMAXINFO* windowInfo = (MINMAXINFO*)theEvent->lParam;

            windowInfo->ptMinTrackSize.x = 10;
            windowInfo->ptMinTrackSize.y = 10;
            windowInfo->ptMaxSize.x=windowInfo->ptMaxTrackSize.x = 0x7777;
            windowInfo->ptMaxSize.y=windowInfo->ptMaxTrackSize.y = 0x7777;

            return false;
        break;
    }*/

    return _HYPlatformWindow::_ProcessOSEvent (vEvent);
}

//__________________________________________________________________
//EOF
