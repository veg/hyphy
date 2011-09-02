/*
    A general list object

    Sergei L. Kosakovsky Pond, August 2000.
*/

#ifndef _HYLIST_
#define _HYLIST_
//#pragma once
#include "HYComponent.h"
#include "HYPlatformList.h"

//__________________________________________________________________

class _HYList: public _HYComponent, public _HYPlatformList
{

public:

    _HYList(_HYRect,Ptr);

    virtual ~_HYList(void);

    virtual void        InsertItem   (_String&, long);
    virtual _String*    GetItem      (long);
    virtual void        SetItem      (_String&, long);
    virtual void        DeleteItem   (long);
    virtual long        ItemCount    (void) {
        return items.lLength;
    }


    virtual void        ToggleMultSelection (bool);
    virtual bool        HaveMultSelection   (void) {
        return multipleSelections;
    }
    virtual _SimpleList&GetSelection        (void) {
        return selection;
    }
    virtual void        SendSelectionChange ();
    virtual void        SendDblClickEvent   (long);
    virtual void        KillSelection       (void);
    virtual bool        ProcessEvent (_HYEvent* e);
    virtual void        Paint (Ptr p) {
        _HYPlatformList::_Paint(p);
        _HYPlatformComponent::_Paint(p);
    }
    virtual void        Update(Ptr p) {
        _HYPlatformList::_Update(p);
        _HYPlatformComponent::_Update(p);
    }
    virtual void        SetSelection  (_SimpleList&);
    virtual void        SetDimensions(_HYRect r, _HYRect rel) {
        _HYComponent::SetDimensions(r,rel);
        _HYPlatformList::_SetDimensions(r,rel);
    }
    virtual void        SetVisibleSize(_HYRect);
    virtual bool        _ProcessOSEvent (Ptr p);
    virtual void        _EraseRect(_HYRect);
    virtual void        SetFont (_HYFont&);
    virtual _HYFont&    GetFont (void) {
        return listFont;
    }
    virtual void        _ToggleDrawing (bool);
    virtual void        _Activate       (void);
    virtual void        _Deactivate     (void);

protected:

    _List       items;
    _SimpleList selection;
    bool        multipleSelections;
    _HYFont     listFont;
};

//__________________________________________________________________

#endif

//EOF
