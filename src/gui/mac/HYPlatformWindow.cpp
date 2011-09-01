/*
    A Mac OS window object - a window/title/size-box/scroll-bars handler.

    Sergei L. Kosakovsky Pond, Spring 2000 - Fall 2001.
*/

#include "HYWindow.h"
#include "HYGWindow.h"
#include "HYEventTypes.h"
#include "HYConsoleWindow.h"
#include "Controls.h"
#include "ControlDefinitions.h"
#include "HYTableWindow.h"
#include "HYUtils.h"
#include "Scrap.h"

//__________________________________________________________________
extern      _SimpleList windowPtrs,
            windowObjects;

extern      _String     objectInspectorTitle;

long                    smallScrollStep = 5,
                        scrollStepCounter = 0,
                        lastScrollControlValue = 0;

_HYGuiObject*           scrollingWindow;
bool                    hScrollingAction = false;


pascal      void        scrollAction                    (ControlHandle,ControlPartCode);


#ifdef TARGET_API_MAC_CARBON
extern  ModalFilterUPP   myFilterProc;
#else
extern  UniversalProcPtr myFilterProc;
#endif

Point                    lastScrollPoint;



//__________________________________________________________________
pascal void  scrollAction (ControlHandle theControl,ControlPartCode ctlPart)
{
    long   cv  = GetControl32BitValue (theControl),
           cv2 = cv;

    scrollStepCounter++;

    switch (ctlPart) {
    case kControlDownButtonPart:
        HiliteControl (theControl,kControlDownButtonPart);
        cv2 = cv+smallScrollStep;
        if (scrollStepCounter>100) {
            cv2 += 99*smallScrollStep;
        } else if (scrollStepCounter>10) {
            cv2 += 9*smallScrollStep;
        }
        if (cv2>MAX_CONTROL_VALUE) {
            cv2 = MAX_CONTROL_VALUE;
        }
        break;
    case kControlUpButtonPart:
        HiliteControl (theControl,kControlUpButtonPart);
        cv2 = cv-smallScrollStep;
        if (scrollStepCounter>100) {
            cv2 -= 99*smallScrollStep;
        } else if (scrollStepCounter>10) {
            cv2 -= 9*smallScrollStep;
        }
        if (cv2<0) {
            cv2 = 0;
        }
        break;
    case kControlPageUpPart:
        HiliteControl (theControl,kControlPageUpPart);
        cv2 = cv-100*smallScrollStep;
        if (cv2<0) {
            cv2 = 0;
        }
        break;
    case kControlPageDownPart:
        HiliteControl (theControl,kControlPageDownPart);
        cv2 = cv+100*smallScrollStep;
        if (cv2>MAX_CONTROL_VALUE) {
            cv2 = MAX_CONTROL_VALUE;
        }
        break;
    default:
        cv2 = cv;
        cv  = lastScrollControlValue;
        break;
    }

    if (cv!=cv2) {
        SetControl32BitValue (theControl,cv2);
        if (hScrollingAction) {
            scrollingWindow->ProcessEvent (generateScrollEvent (cv2-cv,0));
        } else {
            scrollingWindow->ProcessEvent (generateScrollEvent (0,cv2-cv));
        }
    }

    lastScrollControlValue = GetControl32BitValue (theControl);
}

#ifdef TARGET_API_MAC_CARBON


//__________________________________________________________________

pascal OSStatus wSizeHandler (EventHandlerCallRef, EventRef theEvent, void* userData)
{

    EventParamType  actType;
    UInt32          attr = 0;
    GetEventParameter (theEvent, kEventParamAttributes, typeUInt32, &actType,sizeof(UInt32),nil,&attr);

    if (attr&kWindowBoundsChangeSizeChanged) {
        Rect       newSize;
        GetEventParameter (theEvent, kEventParamCurrentBounds, typeQDRectangle, &actType,sizeof(Rect),nil,&newSize);

        _HYWindow* thisWindow = (_HYWindow*)windowObjectRefs(windowPtrs.Find((long)userData));
        forceUpdateForScrolling = true;
        //printf ("%d %d %d %d\n", newSize.top,newSize.left,newSize.bottom+1,newSize.right+1);
        thisWindow->SetWindowRectangle (newSize.top,newSize.left,newSize.bottom+1,newSize.right+1,false);
        thisWindow->Update(nil);
        forceUpdateForScrolling = false;
        return noErr;
    }

    return eventNotHandledErr;
}


#endif



//__________________________________________________________________
void        _PaintTheCircle (CIconHandle c, WindowPtr theWindow)
{
#ifdef OPAQUE_TOOLBOX_STRUCTS
    Rect r;
    GetWindowBounds (theWindow,kWindowGlobalPortRgn,&r);
    OffsetRect (&r,-r.left,-r.top);
#else
    Rect r = theWindow->portRect;
#endif
    r.right = r.left+15;
    r.left  += 3;
    r.bottom --;
    r.top = r.bottom-12;
    PlotCIcon (&r,c);
}

//__________________________________________________________________

_HYPlatformWindow::_HYPlatformWindow(unsigned char windowFlag,_String windowTitle,bool windowVisibility, Ptr cRef)
{
    Rect defRect = {(50+25*windowObjects.lLength)%500,(50+25*windowObjects.lLength)%400,300,350};
    Str255 theTitle;
    StringToStr255 (windowTitle,theTitle);
    long   proc;
    sheet1st = true;
    if (cRef) {
        if (aquaInterfaceOn) {
            if (IsWindowVisible(((_HYWindow*)cRef)->theWindow)) {
                containerRef = cRef;
            } else {
                containerRef = nil;
            }
        } else {
            windowFlag -= HY_WINDOW_SHEET;
            windowFlag += HY_WINDOW_DLOG;
            containerRef = nil;
        }
    } else {
        containerRef = nil;
    }

#ifdef TARGET_API_MAC_CARBON
    if (windowFlag&HY_WINDOW_SHEET) {
        proc = kWindowSheetProc;
    } else
#endif
        if (windowFlag&HY_WINDOW_SIZE) {
            if (windowFlag&HY_WINDOW_ZOOM) {
                proc = zoomDocProc;
            } else if (windowFlag & HY_WINDOW_FLUSHED)
#ifdef TARGET_API_MAC_CARBON
                if (aquaInterfaceOn) {
                    proc =  zoomDocProc;
                } else {
                    proc =  kWindowFloatGrowProc;
                }
#else
                proc =  kWindowFloatGrowProc;
#endif
            else {
                proc = documentProc;
            }
        } else {
            if (windowFlag & HY_WINDOW_DLOG) {
                proc =  kWindowMovableModalDialogProc;
            } else {
                proc = noGrowDocProc;
            }
        }

    savedLoc = defRect;

    flags = windowFlag;

    if (windowFlag & HY_WINDOW_DLOG) {
#ifdef OPAQUE_TOOLBOX_STRUCTS
        DialogPtr temp = NewColorDialog (nil,&defRect,theTitle,windowVisibility,
                                         proc,(WindowPtr)-1,
                                         windowFlag&HY_WINDOW_CLOSE ,0, nil);
        if (temp) {
            theWindow = GetDialogWindow (temp);
        } else {
            theWindow = nil;
        }
#else
        theWindow = NewColorDialog (nil,&defRect,theTitle,windowVisibility,
                                    proc,(WindowPtr)-1,
                                    windowFlag&HY_WINDOW_CLOSE ,0, nil);
#endif
        if (theWindow) {
            GrafPtr savePort;
            GetPort (&savePort);
#ifdef OPAQUE_TOOLBOX_STRUCTS
            SetPort (GetWindowPort(theWindow));
#else
            SetPort (theWindow);
#endif
            SetThemeWindowBackground (theWindow,kThemeBrushDialogBackgroundActive,false);
            SetPort (savePort);
        }
    } else {
#ifdef TARGET_API_MAC_CARBON
        if (cRef) {
            CreateNewWindow (kSheetWindowClass, kWindowNoAttributes,
                             &defRect, &theWindow);
            //if (theWindow)
            //theWindow = GetDialogWindow (theWindow);
            if (theWindow) {
                GrafPtr savePort;
                GetPort (&savePort);
#ifdef OPAQUE_TOOLBOX_STRUCTS
                SetPort (GetWindowPort(theWindow));
#else
                SetPort (theWindow);
#endif
                SetThemeWindowBackground (theWindow,kThemeBrushDialogBackgroundActive,false);
                SetPort (savePort);
            }
        } else
#endif
            theWindow = NewCWindow (nil,&defRect,theTitle,windowVisibility,
                                    proc,(WindowPtr)-1,
                                    windowFlag&HY_WINDOW_CLOSE ,0);
    }

    if (theWindow) {
        if ((windowFlag&HY_WINDOW_SCROLL)&&(windowFlag&HY_WINDOW_SIZE)) {
            Rect cSize = newHRect();
            hScroll = NewControl (theWindow,&cSize,"\p",true,0,0,MAX_CONTROL_VALUE,kControlScrollBarLiveProc,0);
            cSize = newVRect();
            vScroll = NewControl (theWindow,&cSize,"\p",true,0,0,MAX_CONTROL_VALUE,kControlScrollBarLiveProc,0);
            if (!(hScroll&&vScroll)) {
                theWindow = nil;
            } else {
                SetControl32BitMinimum (hScroll,0);
                SetControl32BitMaximum (hScroll,MAX_CONTROL_VALUE);
                SetControl32BitMinimum (vScroll,0);
                SetControl32BitMaximum (vScroll,MAX_CONTROL_VALUE);
            }
        } else {
            hScroll = vScroll = nil;
            windowFlag &= 123;
        }
    }
    if (!theWindow) {
        _String errMsg = "Could not allocate memory for a window structure";
        FlagError (errMsg);
    }
    windowPtrs<<(long)theWindow;
    windowObjects<<(long)this;

    if (((flags&HY_WINDOW_NOLIST)||(flags&HY_WINDOW_DLOG)||(flags&HY_WINDOW_SHEET))==0) {
        MenuHandle wMenu = GetMenuHandle(132);
        AppendMenu (wMenu,"\pa");
        if (theTitle[0]>250) {
            theTitle[0] = 250;
        }

        SetMenuItemText (wMenu, CountMenuItems(wMenu), theTitle);
    }

#ifdef TARGET_API_MAC_CARBON
    if (windowFlag&HY_WINDOW_SIZE) {
        ChangeWindowAttributes (theWindow,kWindowLiveResizeAttribute, kWindowNoAttributes);
        sizeHandler = NewEventHandlerUPP ((EventHandlerProcPtr)wSizeHandler);
        checkPointer  ((Ptr)sizeHandler);
        EventTypeSpec windowSizeEvent;
        windowSizeEvent.eventClass= kEventClassWindow;
        windowSizeEvent.eventKind = kEventWindowBoundsChanged;
        InstallWindowEventHandler (theWindow,sizeHandler,1,&windowSizeEvent,(Ptr)theWindow,NULL);
    } else {
        sizeHandler = nil;
    }
#endif
}

//__________________________________________________________________

_HYPlatformWindow::~_HYPlatformWindow(void)
{
    if (hScroll) {
        DisposeControl (hScroll);
    }
    if (vScroll) {
        DisposeControl (vScroll);
    }
    if (theWindow) {
#ifdef TARGET_API_MAC_CARBON
        if (containerRef) {
            HideSheetWindow (theWindow);
        }
#endif
        if (flags&HY_WINDOW_DLOG)
#ifdef OPAQUE_TOOLBOX_STRUCTS
            DisposeDialog (GetDialogFromWindow(theWindow));
#else
            DisposeDialog (theWindow);
#endif
        else {
            DisposeWindow (theWindow);
        }
    }

#ifdef TARGET_API_MAC_CARBON
    if (sizeHandler) {
        DisposeEventHandlerUPP (sizeHandler);
    }
#endif
}

//__________________________________________________________________

void _HYPlatformWindow::_SetTitle(_String windowTitle)
{
    Str255 theTitle;
    StringToStr255 (windowTitle,theTitle);
    SetWTitle (theWindow,theTitle);
    long f = windowPtrs.Find((long)theWindow);
    if (f>=0) {
        long f2 = FindWindowByName (objectInspectorTitle);
        if((f2>=0)&&(f>=f2)) {
            f--;
        }
        SetMenuItemText (GetMenuHandle(132),f+6,theTitle);
    }
}

//__________________________________________________________________

void _HYPlatformWindow::_Show(void)
{
#ifdef TARGET_API_MAC_CARBON
    if (containerRef&&sheet1st) {
        ShowSheetWindow (theWindow, ((_HYWindow*)containerRef)->theWindow);
        sheet1st = false;
    } else
#endif

        ShowHide (theWindow,true);
}

//__________________________________________________________________

_HYRect _HYPlatformWindow::_GetWindowRect(void)
{
    _HYRect res;

#ifdef OPAQUE_TOOLBOX_STRUCTS
    Rect    wr;
    GetWindowBounds (theWindow,kWindowGlobalPortRgn,&wr);
#else
    Rect    wr = (*((WindowPeek)theWindow)->contRgn)->rgnBBox;
#endif

    res.left    = wr.left;
    res.top     = wr.top;

    res.right     = wr.right;
    res.bottom    = wr.bottom;

    res.width     = 0;

    return res;
}


//__________________________________________________________________

void _HYPlatformWindow::_Hide(void)
{
    ShowHide (theWindow,false);
}

//__________________________________________________________________

long _HYPlatformWindow::_Grow(Ptr theData)
{
    EventRecord* theEvent = (EventRecord*)theData;
    Rect         sizeRect;
    sizeRect.top = 50;
    sizeRect.left = 50;
    _HYWindow* theParent = (_HYWindow*)this;
    sizeRect.bottom = theParent->contentHeight+15;
    sizeRect.right  = theParent->contentWidth+15;

    return       GrowWindow  (theWindow,theEvent->where,&sizeRect);

}

//__________________________________________________________________

_String& _HYPlatformWindow::_GetTitle(void)
{
    return      ((_HYWindow*)this)->GetTitle();
}

//__________________________________________________________________

void _HYPlatformWindow::_Move(Ptr theData)
{
    EventRecord* theEvent = (EventRecord*)theData;
    Rect         allOfIt = {0,0,0x7fff,0x7fff};
    DragWindow  (theWindow,theEvent->where,&allOfIt);
}

//__________________________________________________________________

void _HYPlatformWindow::_SetPosition(int l,int t)
{
    MoveWindow (theWindow,l,t,false);
    long deltah = savedLoc.left-l,
         deltav = savedLoc.top-t;
    savedLoc.left = l;
    savedLoc.right -= deltah;
    savedLoc.top = t;
    savedLoc.bottom -= deltav;
}

//__________________________________________________________________

bool _HYPlatformWindow::_Close(Ptr theData)
{
    _HYWindow* theParent = (_HYWindow*)this;
    bool        doit = true;

    if (theData) {
        EventRecord* theEvent = (EventRecord*)theData;
        doit = theData?TrackGoAway  (theWindow,theEvent->where):true;
    }
    if (doit) {
        if (theParent->ConfirmClose()) {
            long f = windowObjects.Find((long)this);
            if (f>=0) {
                windowObjects.Delete(f);
                windowPtrs.Delete(f);
                windowObjectRefs.Delete(f);
                if ((flags&HY_WINDOW_NOLIST)==0) {
                    MenuHandle  windowM = GetMenuHandle(132);
                    long k = CountMenuItems(windowM);
                    for (; k>=6; k--) {
                        Str255 buffer;
                        _String conv;
                        GetMenuItemText (windowM,k,buffer);
                        Str255ToStr (conv,buffer);

                        _String * myTitle = &theParent->GetTitle();
                        if (conv.Equal (myTitle)) {
                            DeleteMenuItem (windowM,k);
                            break;
                        }
                    }
                }
            }
            //theWindow = nil;
            if (FrontWindow() == theWindow) {
                _UnsetMenuBar ();
            }
        } else {
            doit = false;
        }
    }
    return doit;
}

//__________________________________________________________________

void _Close2(Ptr p)
{
    WindowPtr theWindow = (WindowPtr)p;
#ifdef TARGET_API_MAC_CARBON
    DisposeWindow (theWindow);
#else
    CloseWindow (theWindow);
#endif
}
//__________________________________________________________________
void _HYPlatformWindow::drawGrowIcon (void)
{
    GrafPtr savedPort;
    GetPort(&savedPort);
#ifdef OPAQUE_TOOLBOX_STRUCTS
    SetPort(GetWindowPort(theWindow));
#else
    SetPort(theWindow);
#endif
    RgnHandle   oldClip;
    oldClip = NewRgn();
    Rect r = newSRect();
    GetClip (oldClip);
    ClipRect (&r);
    DrawGrowIcon (theWindow);
    SetClip (oldClip);
    DisposeRgn (oldClip);
    SetPort(savedPort);
}
//__________________________________________________________________

void _HYPlatformWindow::_Activate(void)
{
    //_String    deAC = _String("Activate ") & (long)((_HYWindow*)this)->GetID() & '\n';
    //StringToConsole (deAC);
#ifdef TARGET_API_MAC_CARBON
    if (containerRef&&sheet1st) {
        sheet1st = false;
        ShowSheetWindow (theWindow, ((_HYWindow*)containerRef)->theWindow);
    }

    else
#endif
    {
        ShowWindow   (theWindow);
        SelectWindow (theWindow);
        _SetMenuBar();
    }

    if (flags&HY_WINDOW_SIZE) {
        drawGrowIcon();
    }
    if (hScroll) {
        ShowControl (hScroll);
        ShowControl (vScroll);
        _HYWindow* theParent = (_HYWindow*)this;
        int     windowW = theParent->right-theParent->left-13,
                windowH =  theParent->bottom-theParent->top-13;
        if (theParent->contentWidth>windowW) {
            HiliteControl (hScroll,0);
        }
        if (theParent->contentHeight>windowH) {
            HiliteControl (vScroll,0);
        }
        DrawControls(theWindow);
    }
}

//__________________________________________________________________

void _HYPlatformWindow::_BringWindowToFront(void)
{
#ifdef TARGET_API_MAC_CARBON
    if (containerRef&&sheet1st) {
        ShowSheetWindow (theWindow, ((_HYWindow*)containerRef)->theWindow);
        sheet1st = false;
    }

    else
#endif
    {
        ShowWindow   (theWindow);
        SelectWindow (theWindow);
    }
}

//__________________________________________________________________

void _HYPlatformWindow::_Deactivate(void)
{
    if (!(flags&HY_WINDOW_SHEET)) {

        //_String    deAC = _String("Deactivate ") & (long)((_HYWindow*)this)->GetID() & '\n';
        //StringToConsole (deAC);
        if (flags&HY_WINDOW_SIZE) {
            drawGrowIcon ();
        }
        if (hScroll) {
            HiliteControl (hScroll, kControlInactivePart);
            HiliteControl (vScroll, kControlInactivePart);
            DrawControls(theWindow);
        }
        _UnsetMenuBar();
#ifdef TARGET_API_MAC_CARBON
        Cursor arrow;
        SetCursor(GetQDGlobalsArrow(&arrow));
#else
        SetCursor (&qd.arrow);
#endif
    }
}



//__________________________________________________________________

void _HYPlatformWindow::_SetMenuBar(void)
{
    //_String mOut = _String("_SetMenuBar in ") & (long)((_HYWindow*)this)->GetID() & "\n";
    //StringToConsole (mOut);
    if (flags& HY_WINDOW_DLOG) {
        for (long j = 129; j<133; j++) {
            MenuHandle mh = GetMenuHandle (j);
            DisableMenuItem (mh,0);
        }
    } else {
        MenuHandle  t = GetMenuHandle (130);
        DisableMenuItem (t,1);
        DisableMenuItem (t,3);
        DisableMenuItem (t,4);
        DisableMenuItem (t,5);
        DisableMenuItem (t,6);
        DisableMenuItem (t,8);
        DisableMenuItem (t,9);

        EnableMenuItem  (GetMenuHandle (129),5);
        _HYWindow* theParent = (_HYWindow*)this;
        if (!theParent->IsSaveEnabled()) {
            DisableMenuItem (GetMenuHandle (129),4);
        }
        if (!theParent->IsPrintEnabled()) {
            DisableMenuItem (GetMenuHandle (129),8);
        }
    }

    InvalMenuBar();
}

//__________________________________________________________________

void _HYPlatformWindow::_UnsetMenuBar(void)
{
    //_String mOut = _String("_UnsetMenuBar in ") & (long)((_HYWindow*)this)->GetID() & "\n";
    //StringToConsole (mOut);
    if (flags& HY_WINDOW_DLOG) {
        for (long j = 129; j<133; j++) {
            MenuHandle mh = GetMenuHandle (j);
            EnableMenuItem (mh,0);
        }
    } else {
        MenuHandle  t = GetMenuHandle (129);
        EnableMenuItem (t,1);
        EnableMenuItem (t,2);
        EnableMenuItem (t,4);
        EnableMenuItem (t,8);
        DisableMenuItem(t,5);
        t = GetMenuHandle (130);
        EnableMenuItem (t,8);
        EnableMenuItem (t,9);
    }
    InvalMenuBar();
}

//__________________________________________________________________

void _HYPlatformWindow::_Paint(Ptr)
{
    if (flags&HY_WINDOW_SIZE) {
        DrawGrowIcon (theWindow);
    }
    if (hScroll) {
        ShowControl (hScroll);
        ShowControl (vScroll);
        DrawControls(theWindow);
    }
}

//__________________________________________________________________

void _HYPlatformWindow::_Update(Ptr)
{
    GrafPtr savedPort;
    GetPort(&savedPort);
#ifdef OPAQUE_TOOLBOX_STRUCTS
    SetPort(GetWindowPort(theWindow));
    Rect    portRect;
    GetWindowBounds (theWindow,kWindowGlobalPortRgn,&portRect);
    OffsetRect (&portRect,-portRect.left,-portRect.top);
    EraseRect (&portRect);
#else
    SetPort(theWindow);
    EraseRect (&theWindow->portRect);
#endif
    BeginUpdate(theWindow);
    EndUpdate(theWindow);
    if (flags&HY_WINDOW_SIZE) {
        DrawGrowIcon (theWindow);
    }
    if (hScroll) {
        DrawControls (theWindow);
    }
    SetPort(savedPort);
}

//__________________________________________________________________

void _HYPlatformWindow::_SetWindowRectangle(int top, int left, int bottom, int right, bool ss)
{
    if (ss) {
        SizeWindow (theWindow, right-left, bottom-top, false);
    }
    GrafPtr   savePort;
    GetPort   (&savePort);
#ifdef OPAQUE_TOOLBOX_STRUCTS
    SetPort(GetWindowPort(theWindow));
    Rect    portRect;
    GetWindowBounds (theWindow,kWindowGlobalPortRgn,&portRect);
    OffsetRect (&portRect,-portRect.left,-portRect.top);
    InvalWindowRect (theWindow,&portRect);
#else
    SetPort(theWindow);
    InvalRect (&theWindow->portRect);
#endif
    SetPort   (savePort);
    long      newSize;
    if (hScroll) {
        int windowWidth = right-left-13,
            windowHeight = bottom-top-13;

        _HYWindow* theParent = (_HYWindow*)this;

        if (windowWidth > theParent->contentWidth) {
            HiliteControl (hScroll,255);
        } else {
            HiliteControl (hScroll,0);
            newSize = MAX_CONTROL_VALUE*(_Parameter)windowWidth/(theParent->contentWidth-windowWidth);
            if (newSize>0x6fffffff) {
                newSize = 0x6fffffff;
            }
            SetControlViewSize (hScroll,(long)newSize);
        }

        if (windowHeight > theParent->contentHeight) {
            HiliteControl (vScroll,255);
        } else {
            HiliteControl (vScroll,0);
            newSize = MAX_CONTROL_VALUE*(_Parameter)windowHeight/(theParent->contentHeight-windowHeight);
            if (newSize>0x6fffffff) {
                newSize = 0x6fffffff;
            }
            SetControlViewSize (vScroll,(long)newSize);
        }

        Rect newRect = newVRect();
        MoveControl (vScroll,newRect.left,newRect.top);
        SizeControl (vScroll,newRect.right-newRect.left,newRect.bottom-newRect.top);
        newRect     = newHRect();
        MoveControl (hScroll,newRect.left,newRect.top);
        SizeControl (hScroll,newRect.right-newRect.left,newRect.bottom-newRect.top);
    }
    savedLoc.bottom=savedLoc.top+bottom-top;
    savedLoc.right=savedLoc.left+right-left;
}

//__________________________________________________________________

bool _HYPlatformWindow::_ProcessOSEvent (Ptr vEvent)
{
    EventRecord*    theEvent = (EventRecord*)vEvent;
    WindowPtr       dummy;
    _HYWindow*      theParent = (_HYWindow*)this;
    switch (theEvent->what) {
    case updateEvt:
        if ((WindowPtr)theEvent->message==theWindow) {
            theParent->Update(nil);
            return true;
        }
        return false;
    case activateEvt: {
        if ((WindowPtr)theEvent->message==theWindow) {
            if (theEvent->modifiers&activeFlag) {
                theParent->Activate();
            } else {
                theParent->Deactivate();
            }
        }
        return true;
    }
    case mouseDown: {
        long evtType = FindWindow (theEvent->where,&dummy);
        switch (evtType) {
        case inDrag:
            theParent->Move(vEvent);
            return true;
        case inGrow:
            theParent->Grow(vEvent);
            return true;
        case inGoAway:
            theParent->Close(vEvent);
            return true;
        case inZoomIn:
        case inZoomOut:
            if (TrackBox(theWindow,theEvent->where,evtType)) {
                theParent->_Zoom(evtType==inZoomIn);
            }
            return true;

        case inContent:
            if (FrontWindow()!=theWindow) {
                SelectWindow (theParent->theWindow);
                //theParent->Activate();
                return true;
            } else {
                GrafPtr savedPort;
                GetPort(&savedPort);
#ifdef OPAQUE_TOOLBOX_STRUCTS
                SetPort(GetWindowPort(theWindow));
#else
                SetPort(theWindow);
#endif
                Point localClick = theEvent->where;
                GlobalToLocal (&localClick);
                ControlHandle whichC;
                short f = FindControl (localClick,theWindow,&whichC);
                scrollingWindow = theParent;
                if (f) {
                    // set scroll step
                    short   invisPixels;
                    if ((whichC!=hScroll)&&(whichC!=vScroll)) {
                        return false;
                    }
                    if (whichC==hScroll) {
                        invisPixels = theParent->contentWidth-(theParent->right-theParent->left-13);
                        hScrollingAction = true;
                    } else {
                        invisPixels = theParent->contentHeight-(theParent->bottom-theParent->top-13);
                        hScrollingAction = false;
                    }
                    scrollStepCounter = 0;
                    smallScrollStep = (double)MAX_CONTROL_VALUE/invisPixels;
                    if (!smallScrollStep) {
                        smallScrollStep = 1;
                    }
#ifdef TARGET_API_MAC_CARBON
                    ControlActionUPP myActionProc;
                    myActionProc =   NewControlActionUPP(scrollAction);
#else
                    UniversalProcPtr myActionProc;
                    myActionProc = NewRoutineDescriptor((ProcPtr)scrollAction,
                                                        uppControlActionProcInfo,
                                                        GetCurrentISA());
#endif

                    lastScrollControlValue = GetControl32BitValue (whichC);

                    switch (f) {
                    case kControlIndicatorPart: {
                        long    cv = GetControl32BitValue(whichC),cv2;
                        TrackControl (whichC,localClick,nil);
                        cv2 = GetControl32BitValue(whichC);
                        if (cv!=cv2)
                            if (whichC==hScroll) {
                                scrollingWindow->ProcessEvent (generateScrollEvent (cv2-cv,0));
                            } else {
                                scrollingWindow->ProcessEvent (generateScrollEvent (0,cv2-cv));
                            }

                        break;
                    }
                    default:
                        forceUpdateForScrolling = true;
#ifndef __OLDMAC__
                        TrackControl (whichC,localClick,myActionProc);
#endif
#ifdef __OLDMAC__
                        TrackControl (whichC,localClick,scrollAction);
#endif
                        forceUpdateForScrolling = false;

                    }
#ifdef TARGET_API_MAC_CARBON
                    DisposeControlActionUPP(myActionProc);
#endif
                    return true;
                }
                SetPort(savedPort);
            }

        default:
            return false;
        }
    }
    }
    return false;
}
//__________________________________________________________________

void    _HYPlatformWindow::_SetContentSize (int w, int h)
{
    if (!hScroll) {
        return;
    }
    _HYWindow* theParent = (_HYWindow*)this;
    int     windowW = theParent->right-theParent->left-15,
            windowH =  theParent->bottom-theParent->top-15;
    if (w<=windowW) {
        HiliteControl (hScroll,255);
    } else {
        SetControl32BitValue (hScroll, 0);
        HiliteControl (hScroll,0);
    }
    if (h<=windowH) {
        HiliteControl (vScroll,255);
    } else {
        SetControl32BitValue (vScroll, 0);
        HiliteControl (hScroll,0);
    }

    DrawControls (theWindow);
}

//__________________________________________________________________


void    _HYPlatformWindow::_VisibleContents (int& t,int& l,int& b,int& r)
{
    _HYWindow*
    theParent = (_HYWindow*)this;

    long  windowH   = theParent->bottom - theParent->top  - HY_SCROLLER_WIDTH,
          windowW   = theParent->right  - theParent->left - HY_SCROLLER_WIDTH;

    _Parameter v;

    if (hScroll) {
        if (windowW>theParent->contentWidth) {
            l = theParent->left;
            r = theParent->right;
        } else {
            v = GetControl32BitValue (hScroll);
            l = theParent->left+(theParent->contentWidth-windowW)*v/(double)MAX_CONTROL_VALUE;
            r = l+windowW;
        }
        if (windowH>theParent->contentHeight) {
            t = theParent->top;
            b = theParent->bottom;
        } else {
            v = GetControl32BitValue (vScroll);
            t = theParent->top+(theParent->contentHeight-windowH)*v/(double)MAX_CONTROL_VALUE;
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

void    _HYPlatformWindow::_SetWindowBackColor (_HYColor newColor)
{
    GrafPtr         savePort;
    GetPort         (&savePort);
#ifdef OPAQUE_TOOLBOX_STRUCTS
    SetPort(GetWindowPort(theWindow));
#else
    SetPort(theWindow);
#endif
    RGBColor        c;

    c.red           = newColor.R*256;
    c.blue          = newColor.B*256;
    c.green         = newColor.G*256;

    RGBBackColor    (&c);
    SetPort         (savePort);
}

//__________________________________________________________________

bool    _HYPlatformWindow::_IsHScroll (ControlHandle ch)
{
    return ch == hScroll;
}

//__________________________________________________________________

Rect    _HYPlatformWindow::newVRect (void)
{
#ifdef OPAQUE_TOOLBOX_STRUCTS
    Rect all;
    GetWindowBounds (theWindow,kWindowGlobalPortRgn,&all);
    OffsetRect (&all,-all.left,-all.top);
#else
    Rect all = theWindow->portRect;
#endif
    all.left = all.right-HY_SCROLLER_WIDTH;
    all.top--;
    all.right++;
    all.bottom-=(HY_SCROLLER_WIDTH-1);
    return all;
}

//__________________________________________________________________

Rect    _HYPlatformWindow::newHRect (void)
{
#ifdef OPAQUE_TOOLBOX_STRUCTS
    Rect all;
    GetWindowBounds (theWindow,kWindowGlobalPortRgn,&all);
    OffsetRect (&all,-all.left,-all.top);
#else
    Rect all = theWindow->portRect;
#endif
    all.right -= (HY_SCROLLER_WIDTH-1);
    all.left--;
    all.top = all.bottom-HY_SCROLLER_WIDTH;
    all.bottom++;
    return all;
}

//__________________________________________________________________

Rect    _HYPlatformWindow::newSRect (void)
{
#ifdef OPAQUE_TOOLBOX_STRUCTS
    Rect all;
    GetWindowBounds (theWindow,kWindowGlobalPortRgn,&all);
    OffsetRect (&all,-all.left,-all.top);
#else
    Rect all = theWindow->portRect;
#endif
    all.left = all.right-HY_SCROLLER_WIDTH;
    all.top = all.bottom-HY_SCROLLER_WIDTH;
    return all;
}

//__________________________________________________________________

bool        _HYWindow::_ProcessMenuSelection (long msel)
{
    long        menuChoice = msel&0x0000ffff;

    switch (msel/0xffff) {
    case 129: { // file menu
        if (menuChoice==5) { // close
            Close(nil);
            HiliteMenu(0);
            InvalMenuBar();
            return true;
        }
        if (menuChoice==7)
            // page setup
        {
            OSStatus theStatus;
            Boolean isAccepted;

            PMPrintSession     hyPrintSession;

            //theStatus = PMBegin();
            theStatus = PMCreateSession (&hyPrintSession);
            if (theStatus != noErr) {
                return true;
            }

            if (InitPrint(hyPrintSession)) {
                theStatus = PMSessionPageSetupDialog(hyPrintSession, gPageFormat, &isAccepted);
            }

            if (theStatus == noErr) {
                if (gFlattenedFormat != NULL) {
                    DisposeHandle(gFlattenedFormat);
                    gFlattenedFormat = NULL;
                }

                theStatus = PMFlattenPageFormat(gPageFormat, &gFlattenedFormat);
            }

            if (gPageFormat != kPMNoPageFormat) {
                theStatus = PMRelease (gPageFormat);
                gPageFormat = kPMNoPageFormat;
            }

            theStatus = PMRelease(hyPrintSession);
        }
    }
    }
    return false;
}

//__________________________________________________________________

void _HYPlatformWindow::_SelectWindow (void)
{
    SelectWindow (theWindow);
}

//EOF