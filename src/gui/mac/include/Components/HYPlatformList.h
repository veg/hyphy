/*

    A text input box object for MacOS.



    Sergei L. Kosakovsky Pond, May 2000-December 2002.

*/



#ifndef _HYPLLIST_

#define _HYPLLIST_



#ifdef   TARGET_API_MAC_CARBON

#include "CarbonEvents.h"

#endif



#include "HYPlatformComponent.h"

#include "Lists.h"



//__________________________________________________________________



class _HYPlatformList

{



public:



    _HYPlatformList(void);



    virtual ~_HYPlatformList(void);



    virtual void        _SetDimensions   (_HYRect,_HYRect);

    virtual void        _SetVisibleSize  (_HYRect);



    virtual void        _Paint (Ptr p);

    virtual void        _Update(Ptr p);



    virtual void        _InsertItem         (_String&, long);

    virtual void        _SetItem            (_String&, long);

    virtual void        _DeleteItem         (long);

    virtual void        _CheckSelection     (void);

    virtual void        _ToggleMultSelection(bool);

    void        _KillSelection  (void);

    void        _SetSelection   (_SimpleList&);

    void        _SetFont        (void);



    ListHandle      listData;

    short           fontID;

};



#endif

