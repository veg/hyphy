/*

    A button object for MacOS.



    Sergei L. Kosakovsky Pond, May 2000.

*/



#ifndef _HYPBUTTON_

#define _HYPBUTTON_



#ifdef   TARGET_API_MAC_CARBON

#include "CarbonEvents.h"

#endif



#include "HYPlatformComponent.h"



//__________________________________________________________________



class _HYPlatformButton

{



public:



    _HYPlatformButton(void);



    virtual ~_HYPlatformButton(void);



    virtual void            _SetBackColor    (_HYColor&);

    virtual void            _SetDimensions   (_HYRect,_HYRect);

    virtual void            _SetVisibleSize  (_HYRect);

    virtual void            _SetFont         (_HYFont&);

    void            _SetText         (void);

    void            _SetButtonKind   (unsigned char);

    void            _EnableButton    (bool);

    void            _ApplyFont       (void);

    void            _PaintMe         (void);



    virtual void            _Paint (Ptr p);

    virtual void            _Update(Ptr p);



    PixPatHandle    backFill;

    ControlHandle   buttonControl;

    Rect            buttonRect;

    short           fontID;

};



#endif

