/*

    A text input box object for MacOS.



    Sergei L. Kosakovsky Pond, May 2000-December 2002.

*/



#ifndef _HYPLTABLE_

#define _HYPLTABLE_



#ifdef   TARGET_API_MAC_CARBON

#include "CarbonEvents.h"

#endif



//__________________________________________________________________



#include "HYPlatformComponent.h"



#define     HY_TABLE_SIZE_CURSOR  0x01

#define     HY_TABLE_DRAG_CURSOR  0x02

#define     HY_TABLE_EDIT_CURSOR  0x04



//__________________________________________________________________



class _HYPlatformTable

{



public:



    _HYPlatformTable        (void);

    virtual             ~_HYPlatformTable       (void);



    void            _SetFont                (void);

    void            _SetBackColor           (_HYColor&);

    void            _SetBackColor2          (_HYColor&);



    void            _CreateTextBox          (_HYRect&,_String&);

    _String         _RetrieveTextValue      (void);

    void            _KillTextBox            (void);

    bool            _HasTextBox             (void) {

        return      editBox;

    }



    Rect            _GetVisibleRowRect      (long);



    void            _HiliteRowForDrag       (long,long);



    short           fontID;

    PixPatHandle    backPattern,

                    backPattern2;

    char            cursorState;



    TEHandle        editBox;

    Rect            textBoxRect;

};



#endif

