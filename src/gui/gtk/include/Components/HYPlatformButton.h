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

#ifndef _HYPBUTTON_
#define _HYPBUTTON_

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

    GdkRectangle            buttonRect;
    GdkColor                bgColor;
    GtkWidget*              buttonControl;
    PangoFontDescription*   buttonFontDesc;
    _HYRect                 lastButtonDimension;
};

#endif
