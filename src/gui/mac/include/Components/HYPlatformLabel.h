/*

    A static label object for MacOS.



    Sergei L. Kosakovsky Pond, May 2000-December 2002.

*/



#ifndef _HYPLABEL_

#define _HYPLABEL_



#ifdef   TARGET_API_MAC_CARBON

#include "CarbonEvents.h"

#endif



#include "HYPlatformComponent.h"



//__________________________________________________________________



class _HYPlatformLabel

{



public:



    _HYPlatformLabel(void);

    // flags, title, visibility



    virtual ~_HYPlatformLabel(void);



    virtual void        _SetBackColor    (_HYColor&);

    virtual void        _SetForeColor    (_HYColor&);

    virtual void        _SetDimensions   (_HYRect,_HYRect);

    virtual void        _SetVisibleSize  (_HYRect);

    virtual void        _SetFont         (_HYFont&);

    virtual void        _SetText         (void) {}



    virtual void        _Paint (Ptr p);

    virtual void        _Update(Ptr p);



    PixPatHandle backFill;

    RGBColor     fc;

    long         fontID;

    Rect         labelRect;

};



#endif

