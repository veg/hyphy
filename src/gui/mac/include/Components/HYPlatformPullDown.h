/*

    A pull down menu object for MacOS.



    Sergei L. Kosakovsky Pond, May 2000.

*/



#ifndef _HYPPULLDOWNMENU_

#define _HYPPULLDOWNMENU_

#ifdef   TARGET_API_MAC_CARBON

#include "CarbonEvents.h"

#endif



#include "HYPlatformComponent.h"



//__________________________________________________________________



extern   _String    menuSeparator;



//__________________________________________________________________



class _HYPlatformPullDown

{



public:



    _HYPlatformPullDown(void);

    // flags, title, visibility



    virtual ~_HYPlatformPullDown(void);



    virtual void        _AddMenuItem     (_String&, long);

    virtual void        _SetMenuItem     (_String&, long);

    virtual void        _SetBackColor    (_HYColor&);

    virtual void        _Duplicate       (Ptr);

    virtual void        _DeleteMenuItem  (long);

    virtual long        _GetSelection    (void);

    virtual void        _SetDimensions   (_HYRect,_HYRect);

    virtual void        _SetVisibleSize  (_HYRect);

    virtual void        _EnableItem      (long, bool);

    virtual void        _MarkItem        (long, char);

    virtual char        _ItemMark        (long);

    virtual void        _EnableMenu      (bool) {}





    virtual void        _Paint (Ptr p);

    virtual void        _Update(Ptr p);



    MenuHandle    myMenu;

    PixPatHandle  backFill;

    long          myID, selection;

    Rect          menuRect;

};



//__________________________________________________________________



extern  RGBColor    buttonBorder1,

        buttonBorder2;





#endif

