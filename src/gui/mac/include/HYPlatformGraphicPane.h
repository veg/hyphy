/*
    A painting canvas with double buffer. MacOS.

    Sergei L. Kosakovsky Pond, May 2000-May 2002.
*/

#ifndef _HYPGRAPHICPANE_
#define _HYPGRAPHICPANE_
//#pragma once

#include "HYBaseGUI.h"
#include "QDOffscreen.h"
#include "QuickDraw.h"

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

    GWorldPtr    thePane;

protected:

    CGrafPtr     savedPort;
    GDHandle     savedDevice;
    PixPatHandle fillColor;
    RGBColor     saveFG,
                 saveBG;
};

Rect    HYRect2Rect (_HYRect&);

extern  _List       graphicsFormats;
extern  _SimpleList qtGrexComponents;

void        InitializeQTExporters          (void);

#endif

//EOF
