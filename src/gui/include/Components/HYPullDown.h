/*
    A general pull down menu object.

    Sergei L. Kosakovsky Pond, May 2000.
*/

#ifndef _HYPULLDOWNMENU_
#define _HYPULLDOWNMENU_
//#pragma once


#define  HY_PULLDOWN_NO_MARK      0
#define  HY_PULLDOWN_BULLET_MARK  1
#define  HY_PULLDOWN_CHECK_MARK   2
#define  HY_PULLDOWN_DIAMOND_MARK 4

#include "HYPlatformPullDown.h"
#include "HYComponent.h"

//__________________________________________________________________

class _HYPullDown: public _HYComponent, public _HYPlatformPullDown
{

public:

    _HYPullDown(_HYRect,Ptr);
    // flags, title, visibility

    virtual ~_HYPullDown(void);

    virtual BaseRef     makeDynamic();

    virtual void        Activate     (void);
    virtual void        Deactivate   (void);

    virtual void        SetBackColor (_HYColor);
    virtual _HYColor&   GetBackColor (void);
    virtual void        AddMenuItem  (_String, long);
    virtual void        SetMenuItem  (_String, long);
    virtual void        DeleteMenuItem  (long);
    virtual void        DeleteAllItems  (void);
    virtual _String*    GetMenuItem  (long);
    virtual long        GetSelection (void);
    virtual void        SendSelectionChange (void);
    virtual void        EnableMenu   (bool);
    virtual void        EnableItem   (long, bool);
    virtual bool        IsEnabled    (void);
    virtual void        MarkItem     (long item, char kind) {
        _MarkItem (item,kind);
    }
    virtual char        ItemMark     (long item)            {
        return _ItemMark (item);
    }
    virtual bool        ProcessEvent (_HYEvent* e);
    virtual void        Paint (Ptr p) {
        _HYPlatformPullDown::_Paint(p);
    }
    virtual void        Update(Ptr p) {
        _HYPlatformPullDown::_Update(p);
    }
    virtual void        SetDimensions(_HYRect r, _HYRect rel) {
        _HYComponent::SetDimensions(r,rel);
        _HYPlatformPullDown::_SetDimensions(r,rel);
    }
    virtual void        SetVisibleSize(_HYRect);
    virtual void        SetAlignFlags (unsigned char f) {
        alignFlags = f;
    }
    virtual long        FindMenuItem  (_String& s) {
        return menuSelections.Find(&s);
    }
    virtual unsigned char
    GetAlignFlags (void) {
        return alignFlags;
    }
    virtual void        ChangeSelection (long, bool eventSend = true);
    virtual long        MenuItemCount (void) {
        return menuSelections.lLength;
    }

    virtual bool        _ProcessOSEvent             (Ptr p);
    virtual _HYRect     _SuggestDimensions          (void);
    void        _SetMenuItemTextStyle       (long, char);

private:

    _HYColor        backColor;
    _List           menuSelections;
    bool            enabledFlag;
    unsigned char   alignFlags;
};

//__________________________________________________________________

#endif

//EOF
