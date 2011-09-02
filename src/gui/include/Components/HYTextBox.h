/*
    A text edit box object.

    Sergei L. Kosakovsky Pond, November 2001.
*/

#ifndef _HYTEXTBOX_
#define _HYTEXTBOX_

#include "HYComponent.h"
#include "HYPlatformTextbox.h"

#define  HY_TB_ENABLED 0x01
#define  HY_TB_FOCUSED 0x02
#define  HY_TB_WRAP    0x04
#define  HY_TB_BIGBOX  0x08
#define  HY_TB_ARROWS  0x10

//__________________________________________________________________

class _HYTextBox: public _HYComponent, public _HYPlatformTextBox
{

public:

    _HYTextBox          (_HYRect,Ptr, bool = false);

    virtual                 ~_HYTextBox         (void);

    virtual void            SetBackColor        (_HYColor);
    virtual void            SetForeColor        (_HYColor);
    virtual void            SetBackTColor       (_HYColor);
    virtual _HYColor&       GetBackColor        (void);
    virtual _HYColor&       GetBackTColor       (void);
    virtual _HYColor&       GetForeColor        (void);
    virtual _String         GetText             (void);
    virtual void            StoreText           (_String*&, bool = false);
    virtual void            SetFont             (_HYFont&);
    virtual _HYFont&        GetFont             (void);
    virtual void            SetText             (const _String&, bool update = false);
    virtual void            InsertText          (const _String&, bool update = false, bool append = false);
    virtual void            SetSelection        (long s,long e) {
        _SetSelection(s,e);
    }

    virtual void            Paint               (Ptr p) {
        _HYPlatformTextBox::_Paint(p);
    }
    virtual void            Update              (Ptr p) {
        _HYPlatformTextBox::_Update(p);
    }

    virtual void            SetDimensions       (_HYRect r, _HYRect rel) {
        _HYComponent::SetDimensions(r,rel);
        _HYPlatformTextBox::_SetDimensions(r,rel);
    }

    virtual void            EnableTextEdit      (bool);
    virtual void            SetMargins          (_HYRect newMargins);

    virtual void            SetVisibleSize      (_HYRect);
    virtual void            SetAlignFlags       (unsigned char f) {
        alignFlags = f;
        _SetAlignFlags (f);
    }
    virtual unsigned char
    GetAlignFlags       (void) {
        return alignFlags;
    }

    virtual void            FocusComponent      (void);
    virtual void            UnfocusComponent    (void);
    virtual void            IdleHandler         (void) {
        _IdleHandler();
    }

    virtual bool            _ProcessOSEvent     (Ptr);
    virtual void            _Activate           (void);
    virtual void            _Deactivate         (void);
    virtual void            _SetSelection       (long,long);
    virtual void            _IdleHandler        (void);
    virtual void            _MarkForUpdate      (void);
    virtual bool            _IsEmpty            (void);

    virtual void            _DoUndo             (bool = false);
    virtual void            _DoRedo             (bool = false);
    virtual void            _DoCopy             (bool = false);
    virtual void            _DoCut              (bool = false);
    virtual void            _DoPaste            (bool = false);
    virtual void            _DoSelectAll        (bool = false);
    virtual void            _DoClear            (bool,bool = false);
    virtual void            _DoFind             (_String&);


    _HYFont         editBoxFont;
    _HYColor        backColor,
                    foreColor,
                    backTextColor;

    _HYRect         margins;

    unsigned char   alignFlags;

    int             boxFlags;

    /*bool          isEnabled,
                    isFocused,
                    wrapToView,
                    boxType;*/
};


//__________________________________________________________________

#endif

//EOF
