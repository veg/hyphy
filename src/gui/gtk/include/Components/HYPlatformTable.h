/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (spond@ucsd.edu)
 Art FY Poon    (apoon@cfenet.ubc.ca)
 Steven Weaver (sweaver@ucsd.edu)
 
 Module Developers:
 Lance Hepler (nlhepler@gmail.com)
 Martin Smith (martin.audacis@gmail.com)
 
 Significant contributions from:
 Spencer V Muse (muse@stat.ncsu.edu)
 Simon DW Frost (sdf22@cam.ac.uk)
 
 Permission is hereby granted, free of charge, to any person obtaining a
 copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject to
 the following conditions:
 
 The above copyright notice and this permission notice shall be included
 in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
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
