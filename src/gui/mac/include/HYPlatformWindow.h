/*
    A Mac OS window object - a window/title/size-box/scroll-bars handler.

    Sergei L. Kosakovsky Pond, May 2000.
*/

#ifndef _HYPWINDOW_
#define _HYPWINDOW_

#define  MAX_CONTROL_VALUE 100000000
#include "hy_strings.h"

#include "Quickdraw.h"
#include "Controls.h"
#include "HYBaseGUI.h"

#ifdef   TARGET_API_MAC_CARBON
#include "CarbonEvents.h"
#endif

#include "Carbon.h"

//__________________________________________________________________

class _HYPlatformWindow
{

public:

    _HYPlatformWindow   (unsigned char,_String,bool,Ptr = nil);
    // flags, title, visibility

    virtual ~_HYPlatformWindow  (void);

    void        _SetTitle               (_String);
    void        _Show                   (void);
    void        _Hide                   (void);
    virtual long        _Grow                   (Ptr);
    virtual bool        _Close                  (Ptr);
    virtual void        _Move                   (Ptr);
    virtual void        _SetPosition            (int,int);
    void        _SelectWindow           (void);
    virtual void        _SetWindowRectangle     (int,int,int,int,bool=true);
    virtual void        _SetContentSize         (int,int);
    virtual void        _Paint                  (Ptr);
    virtual void        _Update                 (Ptr);
    virtual void        _Activate               (void);
    virtual void        _BringWindowToFront     (void);
    virtual void        _Deactivate             (void);
    virtual bool        _ProcessOSEvent         (Ptr);
    virtual void        _VisibleContents        (int&,int&,int&,int&);
    virtual _HYRect     _GetWindowRect          (void);
    bool        _IsHScroll              (ControlHandle);
    virtual bool        _ProcessMenuSelection   (long) {
        return false;
    }
    virtual void        _SetMenuBar             (void);
    virtual void        _UnsetMenuBar           (void);
    virtual _String&    _GetTitle               (void);
    void        _SetWindowBackColor     (_HYColor);

    Rect        newVRect                (void);
    Rect        newHRect                (void);
    Rect        newSRect                (void);
    void        drawGrowIcon            (void);
    virtual Ptr         _GetOSWindowData (void) {
        return (Ptr)theWindow;
    }


    WindowPtr                           theWindow;
    unsigned char                       flags;
    ControlHandle                       vScroll,
                                        hScroll;
#ifdef  TARGET_API_MAC_CARBON
    EventHandlerUPP                 sizeHandler;
#endif
    Rect                                savedLoc;
    Ptr                                 containerRef;
    bool                                sheet1st;
};

//__________________________________________________________________

class _HYPlatformTWindow
{

public:

    _HYPlatformTWindow  (Ptr);
    // flags, title, visibility

    virtual ~_HYPlatformTWindow (void);

    virtual void        _SetWindowRectangle     (int,int,int,int,bool=true);

#ifdef      TARGET_API_MAC_CARBON
    EventLoopTimerUPP timerUPP;
    EventLoopTimerRef theTimer;
    EventHandlerUPP   scrollWheelH;
#endif
};

//__________________________________________________________________

class _HYPlatformPWindow
{

public:

    _HYPlatformPWindow          (void);

    virtual ~_HYPlatformPWindow         (void);

    virtual void _StartPicture          (void);
    virtual void _EndPicture            (void);
    virtual void _DrawPicture           (_HYRect);

private:

    PicHandle            savedPic;
    RgnHandle            savedClip;

};

//__________________________________________________________________

extern   bool   forceUpdateForScrolling,
         aquaInterfaceOn;

void            _Close2 (Ptr);

extern  Boolean             InitPrint (PMPrintSession);

extern  CIconHandle         redButtonIcon,
        yellowButtonIcon,
        greenButtonIcon,
        orangeButtonIcon;

#ifdef TARGET_API_MAC_CARBON
#ifndef PM_USE_SESSION_APIS
#define PM_USE_SESSION_APIS 0
#endif

#include <PMApplication.h>
extern PMPageFormat     gPageFormat;
extern PMPrintSettings  gPrintSettings;
extern Handle           gFlattenedFormat;
extern Handle           gFlattenedSettings;
pascal OSStatus         wSizeHandler (EventHandlerCallRef, EventRef, void*);
#else
#include <Printing.h>
extern THPrint prRecHdl;
#endif

void                        _PaintTheCircle                 (CIconHandle,WindowPtr);

#endif

//EOF
