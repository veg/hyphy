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

#ifndef _HYPGRAPHICPANE_
#define _HYPGRAPHICPANE_
//#pragma once

#include "HYBaseGUI.h"
#include <gtk/gtk.h>

//__________________________________________________________________

class  _HYPlatformGraphicPane
{

public:

    _HYPlatformGraphicPane(int, int,int);
    // initial size

    virtual ~_HYPlatformGraphicPane();

    virtual void        _SetPaneSize  (int,int,int);
    virtual void        _DrawLine     (_HYRect);
    virtual void        _DrawHatchedLine(_HYRect);
    // from, to , width
    virtual void        _DisplayText  (_String,int,int,bool);
    // text, where, left-right or top-bottom
    virtual void        _DisplayText  (_String&,_HYRect&, char);

    virtual void        _DisplayChar  (char,int,int);
    // text, where
    virtual void        _DrawRect    (_HYRect);
    virtual void        _FillRect    (_HYRect);
    virtual void        _EraseRect   (_HYRect);
    virtual void        _DrawOval    (_HYRect);
    virtual void        _FillOval    (_HYRect);
    virtual void        _EraseOval   (_HYRect);
    virtual void        _DrawArc     (_HYRect,int,int);
    virtual void        _FillArc     (_HYRect,int,int);
    virtual void        _EraseArc    (_HYRect,int,int);
    virtual void        _SetColor    (_HYColor);
    virtual void        _SetBColor   (_HYColor);
    virtual void        _SetFont     (_HYFont);
    virtual void        _SetFontSize (long);
    virtual void        _StartDraw   (void);
    virtual void        _EndDraw     (void);
    virtual void        _SetPort     (Ptr);
    virtual void        _SlidePane   (int dv, int dh);
    virtual void        _SlideRect   (_HYRect&r ,int dv, int dh);
    virtual void        _InvertRect  (_HYRect&);
    virtual void        _CopyToClipboard
    (void);

    virtual void        _DrawPicRes  (_HYRect&, long);
    virtual void        _SavePicture (_String);

    virtual Ptr         _DefinePolygon
    (_SimpleList&);
    virtual void        _KillPolygon (Ptr);
    virtual void        _DrawPolygon (Ptr,long = 1);
    virtual void        _FillPolygon (Ptr);
    virtual void        _ErasePolygon(Ptr);
    virtual void        _SetDialogBG (void);
    virtual void        _BuildCharGlyphs (void);
    virtual void        _ResetCharGlyphs (void);

    GdkPixmap    *thePane;

    GdkGC*                      theContext;
    PangoLayout                 *textLayout;
    PangoFontDescription        *theFont;

    GdkColor                    saveFG,
                                saveBG,
                                fillColor;

    GList*                      charCachePangoItems;
    void                        *cachedCharacterGlyphs [61];
};

GdkRectangle    HYRect2GDKRect              (const _HYRect&);
_HYRect         GdkRect2HYRect              (const GdkRectangle& hr);
void            HYFont2PangoFontDesc        (const _HYFont&, PangoFontDescription*);

extern          _List       exportFormats;

void            findGraphicsExporterComponents (_List&);


extern                      PangoContext* screenPContext;
#endif

//EOF
