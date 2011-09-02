/*

    At toolbar menu object for MacOS.



    Sergei L. Kosakovsky Pond, May 2000.

*/



#ifndef _HYPBUTTONBAR_

#define _HYPBUTTONBAR_



#include "HYPlatformComponent.h"



#ifdef TARGET_API_MAC_CARBON

#include "CarbonEvents.h"

pascal void ButtonBarTimer (EventLoopTimerRef theTimer,void* userData);

#endif



//__________________________________________________________________



class _HYPlatformButtonBar

{



public:



    _HYPlatformButtonBar(void);



    virtual ~_HYPlatformButtonBar(void);



    virtual void        _SetBackColor    (_HYColor&);

    virtual void        _SetDimensions   (_HYRect,_HYRect);

    virtual void        _SetVisibleSize  (_HYRect);

    _HYRect     _GetButtonRect   (bool conv = false);



    virtual void        _Paint (Ptr p);

    virtual void        _Update(Ptr p);

    void        _DisposeButtons (void);

    void        _DisposeButton  (long);

    void        _MarkButtonForUpdate(int);

    void        _UnpushButton   (void);

    int         _FindClickedButton (int,int);





    PixPatHandle backFill;

    Rect         buttonRect,

                 toolTipBounds;

    int          pushed;

    int          saveMousePosH, saveMousePosV;

    unsigned long

    lastSave;

#ifdef      TARGET_API_MAC_CARBON

    EventLoopTimerUPP timerUPP;

    EventLoopTimerRef theTimer;

#endif

};



//__________________________________________________________________







#endif

