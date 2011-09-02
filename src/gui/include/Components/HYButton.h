/*
    A general label object.

    Sergei L. Kosakovsky Pond, October 2001.
*/

#ifndef _HYBUTTON_
#define _HYBUTTON_

#define  HY_BUTTON_OK     1
#define  HY_BUTTON_CANCEL 2

#include "HYComponent.h"
#include "HYPlatformButton.h"

//__________________________________________________________________

class _HYButton: public _HYComponent, public _HYPlatformButton
{

public:

    _HYButton(_HYRect,Ptr);

    virtual ~_HYButton (void);

    virtual void        SetBackColor        (_HYColor);
    virtual _HYColor&   GetBackColor        (void);
    virtual _String&    GetText             (void);
    virtual void        SetFont             (_HYFont&);
    virtual _HYFont&    GetFont             (void);
    virtual void        SetText             (_String);
    virtual void        Paint               (Ptr p) {
        _HYPlatformButton::_Paint(p);
    }
    virtual void        Update              (Ptr p) {
        _HYPlatformButton::_Update(p);
    }
    virtual void        SetDimensions       (_HYRect r, _HYRect rel) {
        _HYComponent::SetDimensions(r,rel);
        _HYPlatformButton::_SetDimensions(r,rel);
    }

    virtual void        EnableButton        (bool);
    virtual void        SetButtonKind       (unsigned char);
    virtual void        SetVisibleSize      (_HYRect);
    virtual void        SetAlignFlags       (unsigned char f) {
        alignFlags = f;
    }
    virtual unsigned char
    GetAlignFlags       (void) {
        return alignFlags;
    }

    virtual bool        UnfocusedKeyboardInput
    (void) {
        return isEnabled&&buttonKind;
    }

    virtual _HYRect     _SuggestDimensions  (void);
    virtual bool        _ProcessOSEvent     (Ptr);
    virtual void        _Activate           (void);
    virtual void        _Deactivate         (void);



    _String         buttonText;
    _HYColor        backColor;
    _HYFont         buttonFont;
    unsigned char   alignFlags,
             buttonKind;
    bool            isEnabled;
};

//__________________________________________________________________

#endif

//EOF
