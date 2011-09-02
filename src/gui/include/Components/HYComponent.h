/*
    A general composite window component object

    Sergei L. Kosakovsky Pond, May 2000.
*/

#ifndef _HYCOMPONENT_
#define _HYCOMPONENT_
//#pragma once
#include "HYBaseGUI.h"
#include "HYPlatformComponent.h"

#define  HY_COMPONENT_NO_SCROLL  0
#define  HY_COMPONENT_H_SCROLL   1
#define  HY_COMPONENT_V_SCROLL   2
#define  HY_COMPONENT_BORDER_L   4
#define  HY_COMPONENT_BORDER_R   8
#define  HY_COMPONENT_BORDER_T   16
#define  HY_COMPONENT_BORDER_B   32
#define  HY_COMPONENT_BORDER     60
#define  HY_COMPONENT_BORDER_REL 64
#define  HY_COMPONENT_WELL       128
#define  HY_COMPONENT_DRAW_FOCUS 256
#define  HY_COMPONENT_TRANSP_BG  512

#define  HY_ALIGN_LEFT   1
#define  HY_ALIGN_RIGHT  2
#define  HY_ALIGN_TOP    4
#define  HY_ALIGN_BOTTOM 8

//__________________________________________________________________

class _HYComponent: public _HYGuiObject, public _HYPlatformComponent
{

public:

    _HYComponent(void);
    _HYComponent(_HYRect,Ptr);
    // settings, window data

    virtual ~_HYComponent();

    virtual BaseRef     makeDynamic();
    virtual void        Duplicate (BaseRef);

    virtual bool        ProcessEvent (_HYEvent*);
    // by default do nothing

    // a few tool functions

    virtual int         GetMinH (void);
    virtual int         GetMaxH (void);
    virtual int         GetMaxW (void);
    virtual int         GetMinW (void);
    virtual int         GetMaxLH(void);
    virtual int         GetMaxLW(void);
    virtual int         GetHSize(void);
    virtual int         GetVSize(void);
    virtual bool        IsHElastic (void);
    virtual bool        IsVElastic (void);
    virtual void        SetDimensions (_HYRect,_HYRect);
    virtual void        SetOrigin (int,int);
    virtual void        SetVisibleSize (_HYRect);
    virtual _HYRect     VisibleContents(Ptr);
    virtual int         Settings (void) {
        return settings.width;
    }
    virtual bool        HasHScroll (void) {
        return settings.width&HY_COMPONENT_H_SCROLL;
    }
    virtual bool        HasVScroll (void) {
        return settings.width&HY_COMPONENT_V_SCROLL;
    }
    virtual bool        UnfocusedKeyboardInput
    (void) {
        return false;
    }
    virtual void        ComponentMouseExit (void) {
        _ComponentMouseExit ();
    }
    virtual void        FocusComponent  (void) {}
    virtual void        UnfocusComponent(void) {}
    virtual bool        CanCopy         (void) {
        return false;
    }
    virtual BaseRef     CanPaste        (_String&) {
        return nil;
    }

    virtual _String*    HandleCopy      (void) {
        return new _String;
    }
    virtual void        HandlePaste     (BaseRef)
    {}


    // paint functions

    virtual void        Paint           (Ptr);
    virtual void        Update          (Ptr);
    virtual void        Activate        (void) {
        _Activate();
    }
    virtual void        Deactivate      (void) {
        UnfocusComponent();
        _Deactivate();
    }
    virtual void        IdleHandler     (void)
    {}


    virtual void        SetMessageRecipient (_HYGuiObject* o) {
        messageRecipient = o;
    }

    _HYRect     settings;
    // left,top - minimum dimensions
    // right, bottom - maximum dimensions
    // width - scrollable flags
    int         hOrigin, vOrigin, hSize, vSize;
    // display origin for non-self maintaining scrollers
    bool        needUpdate;
    // have I been changed?
    _HYGuiObject*
    messageRecipient;
};

//__________________________________________________________________

#endif

//EOF
