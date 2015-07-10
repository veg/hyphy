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

#ifndef _HYPLTEXTBOX_
#define _HYPLTEXTBOX_

#include "HYPlatformComponent.h"

//__________________________________________________________________

class _HYPlatformTextBox
{

public:

    _HYPlatformTextBox (void);

    virtual ~_HYPlatformTextBox(void);

    virtual void            _SetBackColor    (_HYColor&);
    virtual void            _SetBackTColor   (_HYColor&);
    virtual void            _SetForeColor    (_HYColor&);
    virtual void            _SetDimensions   (_HYRect,_HYRect);
    virtual void            _SetVisibleSize  (_HYRect);
    virtual void            _SetFont         (_HYFont&);
    void            _SetText         (const _String&);
    void            _InsertText      (const _String&, bool);
    _String         _GetText         (void);
    void            _StoreText       (_String*&, bool);
    void            _CreateTE        (void);
    void            _EnableTextBox   (bool);
    void            _SetMargins      (_HYRect&);
    void            _SetAlignFlags   (unsigned char);

    virtual void            _Paint (Ptr p);
    virtual void            _Update(Ptr p);

    virtual void            _FocusComponent  (void);
    virtual void            _UnfocusComponent(void);
    virtual bool            _NeedMultiLines  (void);

    GdkColor        backFill,
                    backTFill,
                    textColor;

    GtkWidget       *te,
                    *scrollWindow;

    GdkRectangle    textBoxRect;
    PangoFontDescription
    *pLabelFont;

    bool            isSingleLine;
};

//__________________________________________________________________

//_String                   retrieveEditControlText (HWND);
//void                  retrieveEditControlText (HWND, _String*&);

#endif
