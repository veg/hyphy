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
    GTK+ glue for the container window

    Sergei L. Kosakovsky Pond, October-November 2004
*/

#include "HYTableWindow.h"
#include "HYEventTypes.h"
#include "HYUtils.h"
#include "errorfns.h"
#include "HYGraphicPane.h"
#include <gdk/gdkkeysyms.h>


//__________________________________________________________________
extern      _SimpleList     windowPtrs,
            windowObjects;

GdkColor    _BLACK_ = {0,0,0,0};
GdkPixmap   *stripedFill        = nil;
GdkGC       *stripedFillGC      = nil;
PangoLayout *statusBarLayout    = nil;
PangoFontDescription
*statusBarFontDesc  = nil;

extern      _String             baseDirectory;

_HYFont     statusBarFont       = {_HY_SANS_FONT,10, HY_FONT_PLAIN};
//__________________________________________________________________

gboolean IdleWindowTimer (Ptr* userData)
{
    _HYTWindow * myTW = (_HYTWindow*)userData;
    for (long k=0; k<myTW->components.lLength; k++)
        if (myTW->cells.Find(k)>=0) {
            _HYComponent* tC = (_HYComponent*)myTW->components(k);
            tC->IdleHandler();
        }
    myTW->_HandleIdleEvent ();
    return true;
}



//__________________________________________________________________

/*  pascal OSStatus scrollWheelHandler (EventHandlerCallRef , EventRef theEvent, void* userData)
    {
        EventParamType                  actType;
        EventMouseWheelAxis             axis;
        GetEventParameter (theEvent,  kEventParamMouseWheelAxis, typeMouseWheelAxis, &actType,sizeof(EventMouseWheelAxis),nil,&axis);

        if (axis == kEventMouseWheelAxisY)
        {
            Point              mouseLocation;
            GetEventParameter (theEvent, kEventParamMouseLocation, typeQDPoint, &actType,sizeof(Point),nil,&mouseLocation);

            long               mouseWheelDelta;
            GetEventParameter (theEvent, kEventParamMouseWheelDelta, typeLongInteger, &actType,sizeof(Point),nil,&mouseWheelDelta);

            UInt32 modifiers;
            GetEventParameter (theEvent, kEventParamKeyModifiers, typeUInt32, NULL, sizeof(modifiers), NULL, &modifiers);

            if (modifiers & optionKey)
                 mouseWheelDelta *= 10;

            _HYTWindow* thisWindow = (_HYTWindow*)windowObjectRefs(windowPtrs.Find((long)userData));
            GrafPtr savedPort;
            GetPort(&savedPort);
            #ifdef OPAQUE_TOOLBOX_STRUCTS
                SetPort(GetWindowPort(thisWindow->theWindow));
            #else
                SetPort(thisWindow->theWindow);
            #endif
            GlobalToLocal (&mouseLocation);
            int c = thisWindow->FindClickedCell(mouseLocation.h,mouseLocation.v);
            thisWindow->DoMouseWheel (c, mouseWheelDelta);
            return noErr;
        }

        return eventNotHandledErr;
    }*/

//__________________________________________________________________


bool        _HYTWindow::_ProcessMenuSelection (long msel)
{
    // TBI
    /*long          menuChoice = msel&0x0000ffff;
    bool        done = false;

    switch (msel/0xffff)
    {
        case 129:
            if (menuChoice == 7) // print setup
            {
                OSStatus theStatus;
                Boolean isAccepted;

                theStatus = PMBegin();
                if (theStatus != noErr)
                    return false;

                if (InitPrint())
                    theStatus = PMPageSetupDialog(gPageFormat, &isAccepted);

                if (theStatus == noErr)
                {
                    if (gFlattenedFormat != NULL)
                    {
                        DisposeHandle(gFlattenedFormat);
                        gFlattenedFormat = NULL;
                    }

                    theStatus = PMFlattenPageFormat(gPageFormat, &gFlattenedFormat);
                }

                if (gPageFormat != kPMNoPageFormat)
                {
                    theStatus = PMDisposePageFormat(gPageFormat);
                    gPageFormat = kPMNoPageFormat;
                }

                theStatus = PMEnd();
                return true;
            }

        case 130:
        {
            if ((menuChoice == 4)||(menuChoice==5))
            {
                done = true;
                HandleCopyPaste(menuChoice-4);
                break;
            }
        }
    }

    HiliteMenu(0);
    InvalMenuBar();*/

    return _HYWindow::_ProcessMenuSelection(msel);
}

//__________________________________________________________________

long _HYTWindow::_Grow(Ptr theData)
{
    return  0;
}

//__________________________________________________________________

void _HYTWindow::_PaintStatusBar(Ptr,bool)
{
    if (GTK_WIDGET_MAPPED (theWindow)) {
        if (!stripedFillGC) {
            SetUpStatusBarStuff (windowContent);
        }

        GdkRectangle    statusBarR = {windowContent->allocation.x,
                                      windowContent->allocation.y+windowContent->allocation.height-HY_SCROLLER_WIDTH,
                                      windowContent->allocation.width,
                                      HY_SCROLLER_WIDTH
                                     };

        gdk_gc_set_fill (stripedFillGC,GDK_TILED);
        gdk_draw_rectangle (GDK_DRAWABLE (theWindow->window), stripedFillGC, true, statusBarR.x, statusBarR.y,
                            statusBarR.width, statusBarR.height);

        gdk_gc_set_fill (stripedFillGC,GDK_SOLID);
        gdk_draw_line(GDK_DRAWABLE (theWindow->window), stripedFillGC,statusBarR.x,statusBarR.y,statusBarR.x+statusBarR.width,statusBarR.y);
        if (statusBar.sLength && statusBarR.width > statusBarFont.size) {
            statusBarR.width  -= statusBarFont.size;
            //statusBarR.height  = statusBarFont.size;
            statusBarR.x += statusBarFont.size/2 + (flags&HY_WINDOW_STATUS_BAR_LIGHT_LEFT)?20:0;
            statusBarR.y += HY_SCROLLER_WIDTH-4*statusBarFont.size/3;

            pango_layout_set_text (statusBarLayout, statusBar.sData, statusBar.sLength);

            gdk_gc_set_clip_rectangle  (stripedFillGC, &statusBarR);
            gdk_draw_layout (GDK_DRAWABLE (theWindow->window), stripedFillGC, statusBarR.x, statusBarR.y, statusBarLayout);
            gdk_gc_set_clip_rectangle  (stripedFillGC, nil);
        }

    }

    // TBI
    /*Rect clearRect = newHRect();
    FillCRect (&clearRect, statusBarFill);
    if (hScroll)
    {
        clearRect = newVRect();
        EraseRect (&clearRect);
    }
    else
    {
        RGBColor saveColor;
        GetForeColor (&saveColor);
        RGBForeColor (&_BLACK_);
        MoveTo (clearRect.left,clearRect.top);
        LineTo (clearRect.right,clearRect.top);
        RGBForeColor (&saveColor);
        if (statusBar.sLength)
        {
            #ifdef OPAQUE_TOOLBOX_STRUCTS
                Rect    destRect;
                GetWindowBounds (theWindow,kWindowGlobalPortRgn,&destRect);
                OffsetRect (&destRect,-destRect.left,-destRect.top);
                if (flags&HY_WINDOW_STATUS_BAR_LIGHT_LEFT)
                    MoveTo (20,destRect.bottom-4);
                else
                    MoveTo (5,destRect.bottom-4);
            #else
                if (flags&HY_WINDOW_STATUS_BAR_LIGHT_LEFT)
                    MoveTo (20,theWindow->portRect.bottom-4);
                else
                    MoveTo (5,theWindow->portRect.bottom-4);
            #endif
            TextFont(0);
            TextFace(0);
            TextSize(9);
            DrawText(statusBar.sData,0,statusBar.sLength);
        }
    }*/
}

//__________________________________________________________________

void _HYTWindow::_Paint(Ptr p)
{
    _HYRect         relRect;
    _SimpleList     alreadyDone (cells.lLength);

    GdkEventExpose * expEvent = (GdkEventExpose*)p;

    for (int k=0; k<cells.lLength; k++) {
        long i = cells.lData[k];

        if (alreadyDone.lData[i] == 0) {
            relRect.left = componentL.lData[i];
            relRect.right = componentR.lData[i];
            relRect.top = componentT.lData[i];
            relRect.bottom = componentB.lData[i];

            if (expEvent) {
                GdkRectangle paintRect = HYRect2GDKRect(relRect);

                paintRect.x += windowContent->allocation.x;
                paintRect.y += windowContent->allocation.y;

                GdkOverlapType otv = gdk_region_rect_in(expEvent->region,&paintRect);

                if (otv != GDK_OVERLAP_RECTANGLE_OUT) {
                    ((_HYComponent*)components(i))->Update((Ptr)&relRect);
                }
            } else {
                ((_HYComponent*)components(i))->Update((Ptr)&relRect);
            }

            alreadyDone.lData[i] = 1;
        }
    }

    if (!hScroll && (flags&HY_WINDOW_SIZE)) {
        _PaintStatusBar();
    }
}

//__________________________________________________________________

void _HYTWindow::_Update(Ptr p)
{
    //_Paint(p);
}

//__________________________________________________________________

void _HYTWindow::_Activate(void)
{
    for (int i=0; i<components.lLength; i++)
        if (cells.Find(i)>=0) {
            ((_HYComponent*)components(i))->Activate();
        }


    if (!theTimer) {
        theTimer = g_timeout_add  (500,(GSourceFunc)IdleWindowTimer,(gpointer)this);
    }

    _HYPlatformWindow::_Activate();

}

//__________________________________________________________________

void _HYTWindow::_Activate2(void)
{
}

//__________________________________________________________________

void _HYTWindow::_Deactivate2(void)
{
}

//__________________________________________________________________

void _HYTWindow::_Zoom(bool inOut)
{
    if ((((savedLoc.right-savedLoc.left!=right-left)||(savedLoc.bottom-savedLoc.top!=bottom-left)))&&(!inOut)) {
        SetPosition (savedLoc.left,savedLoc.top);
        SetWindowRectangle (0,0,savedLoc.bottom-savedLoc.top, savedLoc.right-savedLoc.left);
    } else {
        _HYRect sd = GetScreenDimensions();
        SetPosition (2,2);
        SetWindowRectangle (0,0, sd.bottom-5, sd.right-5);
    }
}

//__________________________________________________________________

void _HYTWindow::_Deactivate(void)
{
    for (int i=0; i<components.lLength; i++)
        if (cells.Find(i)>=0) {
            ((_HYComponent*)components(i))->Deactivate();
        }

    if (theTimer) {
        g_source_remove (theTimer);
        theTimer = 0;
    }
    _HYPlatformWindow::_Deactivate();
}

//__________________________________________________________________
void        _HYTWindow::_SetStatusBar(_String& text)
{
    statusBar = text;
    GdkRectangle    titleBarRect = {theWindow->allocation.x, theWindow->allocation.y+theWindow->allocation.height - HY_SCROLLER_WIDTH,
                                    theWindow->allocation.width, HY_SCROLLER_WIDTH
                                   };

    if (theWindow->window) {
        gdk_window_invalidate_rect (theWindow->window, &titleBarRect, false);
    }
}

//__________________________________________________________________

bool _HYTWindow::_ProcessOSEvent (Ptr vEvent)
{
    _HYWindow* theParent = (_HYWindow*)this;

    _HY_GTK_UI_Message * theMessage = (_HY_GTK_UI_Message*)vEvent;


    gdouble   xc,
              yc;

    if (theMessage->theEvent->type == GDK_KEY_PRESS) {
        GdkEventKey * kpe = (GdkEventKey*)theMessage->theEvent;

        if (keyboardFocusChain.lLength) {
            int keyCode = kpe->keyval;
            if (keyCode == GDK_Tab) { // tab
                bool    backwards = kpe->state & GDK_SHIFT_MASK;

                if (keyboardFocus==-1) {
                    keyCode = keyboardFocusChain.lData[backwards?keyboardFocusChain.lLength-1:0];
                } else if (keyboardFocusChain.lLength>1) {
                    keyCode = keyboardFocusChain.Find (keyboardFocus);
                    keyCode += backwards?(-1):1;
                    if (keyCode<0) {
                        keyCode = keyboardFocusChain.lLength-1;
                    } else if (keyCode >= keyboardFocusChain.lLength) {
                        keyCode = 0;
                    }
                    keyCode = keyboardFocusChain.lData[keyCode];
                } else {
                    keyCode = -1;
                }

                if (keyCode>=0) {
                    ProcessEvent (generateKeyboardFocusEvent (((_HYComponent*)components(keyCode))->GetID()));
                }
                return true;
            }
        }

        for (long k=0; k<components.lLength; k++)
            if (cells.Find(k)>=0) {
                _HYComponent* tC = (_HYComponent*)components(k);
                if (tC->UnfocusedKeyboardInput())
                    if (tC->_ProcessOSEvent (vEvent)) {
                        return true;
                    }
            }

        if ((keyboardFocus>=0)&&(keyboardFocus<components.lLength)) {
            if (((_HYComponent*)components(keyboardFocus))->_ProcessOSEvent (vEvent)) {
                return true;
            }
        }

        return false;
    }


    if (gdk_event_get_coords (theMessage->theEvent,&xc,&yc)) {
        if (trackMouseComponent) {
            return ((_HYComponent*)trackMouseComponent)->_ProcessOSEvent (vEvent);
        }

        xc-=windowContent->allocation.x;
        yc-=windowContent->allocation.y;

        int c = FindClickedCell(xc,yc);

        if (theMessage->theEvent->type == GDK_MOTION_NOTIFY || GDK_LEAVE_NOTIFY == theMessage->theEvent->type) {
            if (c<0) {
                if (lastMouseComponent>=0) {
                    ((_HYComponent*)components(lastMouseComponent))->_ComponentMouseExit();
                }

                lastMouseComponent = -1;
                return false;
            } else {
                if (lastMouseComponent>=0 && c!=lastMouseComponent) {
                    ((_HYComponent*)components(lastMouseComponent))->_ComponentMouseExit();
                }

                //else
                //  SetCursor (LoadCursor (nil, IDC_ARROW));

                lastMouseComponent = c;
                return ((_HYComponent*)components(c))->_ProcessOSEvent (vEvent);
            }
        } else {
            if (c>=0 && (theMessage->processingResult = ((_HYComponent*)components(c))->_ProcessOSEvent (vEvent))) {
                return true;
            }
        }
    }

    return _HYPlatformWindow::_ProcessOSEvent (vEvent);
}

//__________________________________________________________________

void _HYPlatformTWindow::_SetWindowRectangle(int top, int left, int bottom, int right, bool ss)
{
    _HYTWindow* pW = ((_HYTWindow*)this);
    if (ss) {
        _HYRect dr = pW->MinMaxWindowDimensions ();
        GdkGeometry windowG;
        windowG.max_width  = dr.right;
        long    menuHeight = 0;
        if (pW->windowMB) {
            menuHeight = pW->windowMB->allocation.x;
            if (menuHeight <= 0) {
                GtkRequisition sizeReq;
                gtk_widget_size_request (pW->windowMB,&sizeReq);
                menuHeight = sizeReq.height;
            }
        }
        windowG.max_height = dr.bottom+menuHeight;
        windowG.min_width  = dr.left;
        windowG.min_height = dr.top+menuHeight;

        if (windowG.min_height > windowG.max_height) {
            windowG.min_height = windowG.max_height;
        }
        if (windowG.min_width > windowG.max_width) {
            windowG.min_width = windowG.max_width;
        }
        gtk_window_set_geometry_hints (GTK_WINDOW(pW->theWindow), NULL, &windowG, (GdkWindowHints)(GDK_HINT_MAX_SIZE|GDK_HINT_MIN_SIZE));
    }
    //pW->_HYPlatformWindow::_SetWindowRectangle (top,left,bottom,right,ss);
}

//__________________________________________________________________

void        _HYTWindow::_SetCopyString (_String* str)
{
    // TBI
    //ZeroScrap();
    //PutScrap (str->sLength,'TEXT',str->sData);
}

//__________________________________________________________________

_String*        _HYTWindow::_GetPasteString (void)
{
    _String *res = nil;
    // TBI
    /*Handle    scrapHandle = NewHandle (0);
    long    rc;

    #ifdef TARGET_API_MAC_CARBON
        ScrapRef theScrapRef;
        GetCurrentScrap(&theScrapRef);
        GetScrapFlavorSize(theScrapRef, 'TEXT', &rc);
    #else
        long    scrapOffset;
        rc = GetScrap( scrapHandle, 'TEXT', &scrapOffset );
    #endif
    if ( rc >= 0 )
    {
        SetHandleSize( scrapHandle, rc+1 );
        HLock  (scrapHandle);
        #ifdef TARGET_API_MAC_CARBON
            GetScrapFlavorData(theScrapRef, 'TEXT', &rc, *scrapHandle);
        #endif
        (*scrapHandle)[rc] = 0;
        HUnlock (scrapHandle);
        res = new _String (*scrapHandle);
    }
    else
        res = new _String;

    DisposeHandle (scrapHandle);*/
    return res;
}

//__________________________________________________________________

_HYPlatformTWindow::_HYPlatformTWindow(Ptr)
{
    theTimer            = 0;
    trackMouseComponent = nil;
}

//__________________________________________________________________

_HYPlatformTWindow::~_HYPlatformTWindow(void)
{
    if (theTimer) {
        g_source_remove (theTimer);
        theTimer = 0;
    }
}
//EOF
