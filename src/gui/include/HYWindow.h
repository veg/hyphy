/*
    A general window object - a window/title/size-box/scroll-bars handler.

    Sergei L. Kosakovsky Pond, May 2000.
*/

#ifndef _HYWINDOW_
#define _HYWINDOW_
//#pragma once
#include "HYPlatformWindow.h"
#include "HYBaseGUI.h"

#define  HY_WINDOW_SIZE     1
#define  HY_WINDOW_CLOSE    2
#define  HY_WINDOW_SCROLL   4
#define  HY_WINDOW_NOLIST   8
#define  HY_WINDOW_ZOOM     16
#define  HY_WINDOW_DLOG     32
#define  HY_WINDOW_FLUSHED  64
#define  HY_WINDOW_SHEET    128

#define  HY_WINDOW_KIND_BASE 0

//__________________________________________________________________

class _HYWindow: public _HYGuiObject, public _HYPlatformWindow
{

private:

    _String     windowTitle;

public:

    _HYWindow(unsigned char,_String,bool,Ptr = nil);
    // flags, title, visibility

    virtual ~_HYWindow();

//      BaseRef     makeDynamic();
//      void        Duplicate (BaseRef);

    virtual void        SetTitle            (_String);
    virtual _String&    GetTitle            (void) {
        return windowTitle;
    }

    virtual bool        ProcessEvent        (_HYEvent* e) {
        delete e;
        return false;
    }
    virtual bool        ProcessGEvent       (_HYEvent* e);

    void        Show                (void);
    void        Hide                (void);
    virtual void        Grow                (Ptr);
    virtual void        SetWindowRectangle  (int,int,int,int,bool=true);
    virtual bool        Close               (Ptr);
    virtual void        _Zoom               (bool)
    {}
    virtual bool        _ProcessMenuSelection (long);
    virtual bool        ConfirmClose        (void) {
        return true;
    }
    virtual void        Move                (Ptr);
    virtual void        SetPosition         (int, int);
    virtual void        Paint               (Ptr p) {
        _Paint(p);
    }
    virtual void        Update              (Ptr p) {
        _Update(p);
    }
    virtual void        Activate            (void);
    virtual void        BringToFront        (void) {
        _BringWindowToFront();
    }
    virtual void        Deactivate          (void);
    void        SelectThisWindow    (void);
    virtual void        SetContentSize      (int,int);
    virtual void        SetFont             (_HYFont&);

    virtual _HYRect     GetWindowRect       (void) {
        return _GetWindowRect();
    }

    virtual void        VisibleContents     (int&t,int&l,int&b,int&r) {
        _VisibleContents    (t,l,b,r);
    }

    virtual bool        IsSaveEnabled       (void) {
        return false;
    }
    virtual bool        IsPrintEnabled      (void) {
        return false;
    }

    virtual Ptr         GetOSWindowData (void) {
        return _GetOSWindowData();
    }
    virtual char        WindowKind      (void) {
        return HY_WINDOW_KIND_BASE;
    }

    int         top,
                left,
                bottom,
                right,
                contentHeight,
                contentWidth;
};

//__________________________________________________________________

extern  unsigned long GUIObjectGlobalCounter;
extern  _List    GlobalGUIEventQueue;

#endif

//EOF
