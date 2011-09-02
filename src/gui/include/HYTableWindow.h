/*
    A general composite window object

    Sergei L. Kosakovsky Pond, May 2000.
*/

#ifndef _HYTWINDOW_
#define _HYTWINDOW_
//#pragma once
#include "HYComponent.h"
#include "HYWindow.h"

#define  HY_WINDOW_STATUS_BAR_LIGHT_LEFT 64

//__________________________________________________________________

class _HYTWindow: public _HYWindow, public _HYPlatformTWindow
{

public:

    _HYTWindow(_String, char list = 1, bool dlog = false, Ptr sheetParent = nil);

    virtual ~_HYTWindow() {};

    virtual bool         ProcessEvent (_HYEvent*);

    virtual void         SetTableDimensions     (int, int);
    virtual void         Activate               (void);
    int          Rows                   (void) {
        return rows;
    }
    int          Columns                (void) {
        return columns;
    }
    virtual void         SetCell                (int,int,_HYGuiObject*);
    virtual void         AddObject              (_HYGuiObject*, bool = true,long = -1, long = -1);
    virtual void         AddKeyboardChainObject (_HYGuiObject*);
    virtual void         Paint                  (Ptr);
    virtual void         Update                 (Ptr);
    virtual void         SetWindowRectangle     (int,int,int,int,bool=true);
    virtual _HYRect      MinMaxWindowDimensions (void);
    virtual _HYGuiObject*GetCellObject          (int,int);
    virtual _HYGuiObject*GetObject              (int);
    virtual int          FindClickedCell        (int,int);
    virtual void         UpdateComponentInfo    (void);
    virtual void         SetStatusBar           (_String&);
    virtual void         DoMouseWheel           (long, long);

    virtual void         _Zoom                  (bool);
    virtual void         _SetStatusBar          (_String&);
    virtual long         MatchComponentID       (_String&);
    virtual void         HandleCopyPaste        (bool);
    virtual void         _SetCopyString         (_String*);
    virtual _String*     _GetPasteString        (void);

    virtual long         _Grow                  (Ptr);
    virtual void         _Paint                 (Ptr);
    virtual void         _PaintStatusBar        (Ptr = nil,bool=false);
    virtual void         _Update (Ptr);

    virtual void         _Activate              (void);
    virtual void         _Deactivate            (void);

    virtual void         _Activate2             (void);
    virtual void         _Deactivate2           (void);

    virtual void         _HandleIdleEvent       (void)
    {}


    virtual bool         _ProcessOSEvent        (Ptr);
    virtual bool         _ProcessMenuSelection  (long);
    void         RecomputeCellRects     (void);


    _List       components;
    _SimpleList cells,
                componentT,
                componentL,
                componentB,
                componentR,
                keyboardFocusChain;

    int         rows,
                columns,
                keyboardFocus,
                lastMouseComponent;

    _HYRect     dim;

    _String     statusBar;

};

//__________________________________________________________________

#endif

//EOF
