/*

    A checkbox with optional static label object for MacOS.



    Sergei L. Kosakovsky Pond, May 2000-December 2002.

*/



#ifndef _HYPLCHECKBOX_

#define _HYPLCHECKBOX_



#ifdef   TARGET_API_MAC_CARBON

#include "CarbonEvents.h"

#endif



#include "HYPlatformComponent.h"



//__________________________________________________________________



class _HYPlatformCheckbox

{



public:



    _HYPlatformCheckbox(bool);

    virtual ~_HYPlatformCheckbox(void);



    virtual void        _SetVisibleSize     (_HYRect);

    virtual void        _SetState           (bool);



    virtual void        _Paint (Ptr p);

    virtual void        _Update(Ptr p);

    virtual void        _Enable(bool);



    void        _PaintMe            (void);

    ControlHandle   checkboxControl;

    Rect            checkboxRect;

    bool            isRadio;

};



#endif

