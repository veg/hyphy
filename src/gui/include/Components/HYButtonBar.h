/*
    A general button bar object

    Sergei L. Kosakovsky Pond, June 2000.
*/

#ifndef _HYBUTTONBAR_
#define _HYBUTTONBAR_
//#pragma once

#include "HYComponent.h"
#include "HYPlatformButtonBar.h"

#define  HY_BUTTONBAR_BORDER 2

//__________________________________________________________________

class _HYButtonBar: public _HYComponent, public _HYPlatformButtonBar
{

public:

    _HYButtonBar(_HYRect,Ptr);

    virtual ~_HYButtonBar(void);

    virtual void        SetBackColor (_HYColor);
    virtual _HYColor&   GetBackColor (void);
    virtual void        AddButton    (Ptr,_String* toolTip);

    virtual void        ReplaceButton(long,Ptr,_String* = nil);
    virtual Ptr         GetButtonIcon(long);
    virtual void        DeleteButton (long);
    virtual void        EnableButton (long,bool);

    virtual void        SendButtonPush      (int);
    virtual bool        ProcessEvent (_HYEvent* e);
    virtual void        Paint (Ptr p) {
        _HYPlatformButtonBar::_Paint(p);
    }
    virtual void        Update(Ptr p) {
        _HYPlatformButtonBar::_Update(p);
    }
    virtual void        SetDimensions(_HYRect r, _HYRect rel) {
        _HYComponent::SetDimensions(r,rel);
        _HYPlatformButtonBar::_SetDimensions(r,rel);
    }
    virtual void        SetVisibleSize(_HYRect);
    virtual bool        _ProcessOSEvent (Ptr p);
    virtual _HYRect     _SuggestDimensions (void);
    virtual void        SetAlignFlags (unsigned char f) {
        alignFlags = f;
    }
    virtual unsigned char
    GetAlignFlags (void) {
        return alignFlags;
    }
    virtual long        ButtonCount   (void) {
        return buttons.lLength;
    }
    virtual void        SetButtonLayoutW(int w) {
        barW = w;
    }
    virtual int         BarWidth (void) {
        return barW;
    }
    virtual void        SetButtonDim  (int d) {
        buttonDim = d;
    }
    virtual int         GetButtonDim  (void) {
        return buttonDim;
    }
    virtual void        GetButtonLoc  (int, int&, int&, bool);
    virtual void        SetToolTip    (long,_String*);
    void        MarkAsPullDown(long,bool);

    virtual void        _Activate           (void);
    virtual void        _Deactivate         (void);
    virtual void        _ComponentMouseExit (void);
    void        _DisplayToolTip     (void);

    friend  class       _HYPlatformButtonBar;

protected:

    _HYColor        backColor;
    _SimpleList     buttons, enabledButtons, pullDownButtons;
    _List           toolTips;
    unsigned char   alignFlags;
    int             barW,
                    buttonDim;
};

//__________________________________________________________________

#endif

//EOF
