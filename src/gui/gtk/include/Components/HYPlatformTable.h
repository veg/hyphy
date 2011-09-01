/*
    A table object glue for GTK+

    Sergei L. Kosakovsky Pond, March 2005.
*/

#ifndef _HYPLTABLE_
#define _HYPLTABLE_


//__________________________________________________________________

#include "HYPlatformComponent.h"

#define     HY_TABLE_SIZE_CURSOR   0x01
#define     HY_TABLE_DRAG_CURSOR   0x02
#define     HY_TABLE_EDIT_CURSOR   0x04

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

    GdkRectangle    _GetVisibleRowRect      (long);

    void            _HiliteRowForDrag       (long,long);
    void            _ResetCursorState       (void);

    void            _FrameRect              (GdkRectangle&);
    bool            _CheckGC                (void);

    GdkColor        backPattern,
                    backPattern2;

    PangoFontDescription
    *tableFont,
    *tableFontB,
    *tableFontI,
    *tableFontBI;

    char            cursorState;

    GtkWidget*      editBox;
    GdkRectangle    textBoxRect,
                    limits;

    long            activeColumn,
                    activeColumn2;

    GdkGC*          theContext;

};

#endif