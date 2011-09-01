/*
    A checkbox with optional static label object for GTK glue.

    Sergei L. Kosakovsky Pond, November 3rd (this is the end...) 2004
*/

#ifndef _HYPLCHECKBOX_
#define _HYPLCHECKBOX_

#define HY_DEFAULT_CHECK_SPACE 16

#include "HYPlatformComponent.h"

//__________________________________________________________________

class _HYPlatformCheckbox
{

public:

    _HYPlatformCheckbox(bool);
    virtual ~_HYPlatformCheckbox(void);

    virtual void            _SetVisibleSize     (_HYRect);
    virtual void            _SetState           (bool);

    virtual void            _Paint (Ptr p);
    virtual void            _Update(Ptr p);
    virtual void            _Enable(bool);
    virtual long            _GetBoxDimension (void);

    void            _PaintMe            (void);
    GtkWidget       *checkboxControl,
                    *dummyControl;
    _HYRect         checkboxRect;
    bool            isRadio;
};

#endif