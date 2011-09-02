/*
    A window with a picture.

    Sergei L. Kosakovsky Pond, May 2000.

    Revised: August   2001
             November 2001 - added "picture" window.
*/

#ifndef _HYGWINDOW_
#define _HYGWINDOW_
//#pragma once
#include "HYWindow.h"
#include "HYGraphicPane.h"

//__________________________________________________________________

class _HYGWindow: public _HYWindow, public _HYGraphicPane
{

public:

    _HYGWindow(_String,int,int,int,bool,bool scale = false);
    // title, pic h, pic w, pic d, visible

    virtual ~_HYGWindow() {};

    virtual bool        ProcessEvent            (_HYEvent* e);
    virtual bool        _ProcessMenuSelection   (long);

    virtual void        _SetMenuBar             (void);
    virtual void        _UnsetMenuBar           (void);

// overload  platform specific _paint and _update

    virtual void        _Paint(Ptr);
    virtual void        _Update(Ptr);
    void        _SaveGWindow (void);
    void        _PrintGWindow(void);

    virtual bool        IsSaveEnabled     (void) {
        return true;
    }
    virtual bool        IsPrintEnabled    (void) {
        return true;
    }

};

//__________________________________________________________________

class _HYPWindow: public _HYGWindow, public _HYPlatformPWindow
{

public:

    _HYPWindow(_String,int,int,int,bool);
    // title, pic h, pic w, pic d, visible

    virtual ~_HYPWindow() {};

    virtual void        SetPaneSize             (int,int,int);
    virtual void        StartDraw               (void);
    virtual void        EndDraw                 (void);
    virtual void        Zoom                    (_Parameter);
    virtual void        OriginalSize            (void);


    virtual bool        _ProcessMenuSelection   (long);
    virtual bool        _ProcessOSEvent         (Ptr);
    virtual void        _SetWindowRectangle     (int,int,int,int,bool=true);
    virtual void        _Paint(Ptr);
    virtual void        _Update(Ptr);
    virtual long        _Grow (Ptr);
    void        _PrintPWindow           (void);

    virtual void        _SetMenuBar             (void);
    virtual void        _UnsetMenuBar           (void);

    bool        resizing;
    long        oh,
                ow;

};
#endif

//EOF
