/*
    A painting canvas glue for GTK

    Sergei L. Kosakovsky Pond, October 2004.
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