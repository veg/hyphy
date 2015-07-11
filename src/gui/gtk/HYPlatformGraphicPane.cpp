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


/*
    A painting canvas with double buffer glue for GTK+

    Sergei L. Kosakovsky Pond, October 2004.
*/

#include "HYGraphicPane.h"

#include "HYUtils.h"
#include "HYPlatformWindow.h"
#include "HYComponent.h"

#include <gdk/gdk.h>
#include <gdk-pixbuf/gdk-pixbuf-io.h>

_String     savePicPrompt               ("Save Picture As:"),
            savePicAs                  ("File Format:");

_List       exportFormats;

//__________________________________________________________________

void add_if_writable (GdkPixbufFormat *data, _List *list)
{
    if (gdk_pixbuf_format_is_writable (data)) {
        _String tStr (data->name);
        (*list) && & tStr ;
    }
}

//__________________________________________________________________

void        findGraphicsExporterComponents (_List& storeIn)
{
    GSList *formats = gdk_pixbuf_get_formats ();
    g_slist_foreach (formats, (GFunc)add_if_writable, &storeIn);
    g_slist_free (formats);
}

//__________________________________________________________________

GdkRectangle HYRect2GDKRect (const _HYRect& hr)
{
    return (GdkRectangle) {
        hr.left, hr.top, hr.right-hr.left+1, hr.bottom-hr.top+1
    };
}

//__________________________________________________________________

_HYRect GdkRect2HYRect (const GdkRectangle& hr)
{
    return (_HYRect) {
        hr.y, hr.x, hr.y+hr.height-1, hr.x+hr.width-1,0
    };
}

//__________________________________________________________________

void HYFont2PangoFontDesc (const _HYFont& f, PangoFontDescription* theFont)
{
    pango_font_description_set_family (theFont, f.face.sData);
    pango_font_description_set_style  (theFont, (f.style & HY_FONT_ITALIC) ? PANGO_STYLE_ITALIC : PANGO_STYLE_NORMAL);
    pango_font_description_set_weight (theFont, (f.style & HY_FONT_BOLD) ? PANGO_WEIGHT_BOLD : PANGO_WEIGHT_NORMAL);
    pango_font_description_set_size   (theFont, f.size*PANGO_SCALE);
    //pango_font_description_set_absolute_size (theFont, f.size*PANGO_SCALE);
}

//__________________________________________________________________

_HYPlatformGraphicPane::_HYPlatformGraphicPane(int h, int w, int d)
{
    fillColor = (GdkColor) {
        0,0,0,0
    };
    //printf ("Allocating pixmap\n");
    thePane = gdk_pixmap_new (NULL, w, h, 24/*d<24?d:24*/);
    //printf ("Allocating context\n");
    theContext   = gdk_gc_new (thePane);
    //printf ("Allocating context\n");
    textLayout   = pango_layout_new (screenPContext);
    theFont      = pango_font_description_new ();
    //printf ("Setting colormaps\n");
    gdk_drawable_set_colormap (thePane, gdk_colormap_get_system ()); // ?
    gdk_gc_set_colormap (theContext, gdk_colormap_get_system ()); // ?

    charCachePangoItems   = NULL;
    for (long k=0; k<=60; k=k+1) {
        cachedCharacterGlyphs[k] = NULL;
    }
}

//__________________________________________________________________

_HYPlatformGraphicPane::~_HYPlatformGraphicPane(void)
{
    g_object_unref (thePane);
    g_object_unref (theContext);
    g_object_unref (textLayout);
    pango_font_description_free (theFont);
    _ResetCharGlyphs ();
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_ResetCharGlyphs(void)
{
    if (charCachePangoItems) {
        g_list_free (charCachePangoItems);
        charCachePangoItems = NULL;
    }
    for (long k=0; k<=60; k=k+1)
        if (cachedCharacterGlyphs[k]) {
            pango_glyph_string_free ((PangoGlyphString* )cachedCharacterGlyphs[k]);
            cachedCharacterGlyphs[k] = NULL;
        }
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_BuildCharGlyphs(void)
{
    char c[] = "ACGT";
    PangoAttrList*  natl = pango_attr_list_new          ();
    PangoAttribute* pafd = pango_attr_font_desc_new    (theFont);
    pango_attr_list_insert (natl, pafd);
    charCachePangoItems = pango_itemize (screenPContext, c, 0, 4, natl,NULL);
    pango_attr_list_unref (natl);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_SetPaneSize  (int h,int w, int d)
{
    g_object_unref (thePane);
    g_object_unref (theContext);
    thePane = gdk_pixmap_new (NULL, w, h, 24/*d<24?d:24*/);
    theContext   = gdk_gc_new (thePane);
    gdk_drawable_set_colormap (thePane, gdk_colormap_get_system ()); // ?
    gdk_gc_set_colormap (theContext, gdk_colormap_get_system ()); // ?
}


//__________________________________________________________________

void _HYPlatformGraphicPane::_DrawLine (_HYRect lineDesc)
{
    gdk_gc_set_line_attributes (theContext, lineDesc.width, GDK_LINE_SOLID, GDK_CAP_PROJECTING, GDK_JOIN_MITER);
    //if (lineDesc.top == lineDesc.bottom)
    //  gdk_draw_line (thePane, theContext, lineDesc.left, lineDesc.top, lineDesc.right+1, lineDesc.bottom);
    //else
    //  if (lineDesc.left == lineDesc.right)
    //      gdk_draw_line (thePane, theContext, lineDesc.left-1, lineDesc.top, lineDesc.right-1, lineDesc.bottom);
    //else
    gdk_draw_line (thePane, theContext, lineDesc.left, lineDesc.top, lineDesc.right, lineDesc.bottom);
}



//__________________________________________________________________

void _HYPlatformGraphicPane::_DrawHatchedLine (_HYRect lineDesc)
{
    gdk_gc_set_line_attributes (theContext, lineDesc.width, GDK_LINE_ON_OFF_DASH, GDK_CAP_NOT_LAST, GDK_JOIN_MITER);
    gdk_draw_line (thePane, theContext, lineDesc.left, lineDesc.top, lineDesc.right, lineDesc.bottom);
}


//__________________________________________________________________
void _HYPlatformGraphicPane::_DisplayText    (_String theText,int t, int l, bool dir)
{
    _HYGraphicPane* theParent = (_HYGraphicPane*)this;
    pango_layout_set_width(textLayout, -1);
    pango_layout_set_text (textLayout, theText.sData, theText.sLength);
    PangoLayoutLine* aLine = pango_layout_get_line      (textLayout,0);
    //gdk_draw_layout (thePane, theContext, l+1, t-theParent->font.size, textLayout);
    gdk_draw_layout_line (thePane, theContext, l+1,t, aLine);
}



//__________________________________________________________________
void _HYPlatformGraphicPane::_DisplayText    (_String& theText,_HYRect& r, char align)
{
    _HYGraphicPane* theParent = (_HYGraphicPane*)this;
    _EraseRect (r);
    pango_layout_set_width(textLayout, PANGO_SCALE*(r.right-r.left+1));
    pango_layout_set_text (textLayout, theText.sData, theText.sLength);
    GdkRectangle clip_rect = HYRect2GDKRect (r);
    gdk_gc_set_clip_rectangle  (theContext, &clip_rect);
    pango_layout_set_alignment (textLayout, align==HY_ALIGN_LEFT?PANGO_ALIGN_LEFT:(align==HY_ALIGN_RIGHT?PANGO_ALIGN_RIGHT:PANGO_ALIGN_CENTER));
    gdk_draw_layout (thePane, theContext, r.left, r.top, textLayout);
    gdk_gc_set_clip_rectangle  (theContext, NULL);
}

//__________________________________________________________________
void _HYPlatformGraphicPane::_DisplayChar  (char c,int t, int l)
{
    if (c>=40 && c<=100) {
        if (charCachePangoItems == NULL) {
            _BuildCharGlyphs();
        }
        PangoItem * pitem = (PangoItem*) g_list_first (charCachePangoItems)->data;
        if (!cachedCharacterGlyphs[c-40]) {
            PangoGlyphString * newString = pango_glyph_string_new ();
            pango_shape (&c, 1, &pitem->analysis, newString);
            cachedCharacterGlyphs[c-40] = newString;
        }
        gdk_draw_glyphs (thePane, theContext, pitem->analysis.font,l,t,(PangoGlyphString* )cachedCharacterGlyphs[c-40]);
    } else {
        _DisplayText (c, t, l, false);
    }
}

//__________________________________________________________________
void _HYPlatformGraphicPane::_SlidePane  (int dv, int dh)
{
    _HYGraphicPane* theParent = (_HYGraphicPane*)this;
    GdkPixbuf* theImage = gdk_pixbuf_get_from_drawable (NULL, thePane, NULL, 0, 0, 0, 0, -1, -1);
    gdk_draw_pixbuf (thePane, NULL, theImage, 0,0, dh, dv, -1, -1, GDK_RGB_DITHER_NONE, 0, 0);
    g_object_unref (theImage);
}


//__________________________________________________________________
void _HYPlatformGraphicPane::_SlideRect (_HYRect& rct, int dv, int dh)
{
    _HYGraphicPane* theParent = (_HYGraphicPane*)this;

    GdkRectangle clipped_rect   = HYRect2GDKRect (rct),
                 full_rect        = {0,0,theParent->w, theParent->h},
                 clip_rect;

    if (gdk_rectangle_intersect(&clipped_rect, &full_rect, &clip_rect)) {
        gdk_gc_set_clip_rectangle  (theContext, &clip_rect);
        GdkPixbuf* theImage = gdk_pixbuf_get_from_drawable (NULL, thePane, NULL, rct.left, rct.top, 0, 0, clip_rect.width, clip_rect.height);
        gdk_draw_pixbuf (thePane, theContext, theImage, 0,0, clip_rect.x+dh, clip_rect.y+dv, -1,-1, GDK_RGB_DITHER_NONE, 0, 0);
        g_object_unref (theImage);
        gdk_gc_set_clip_rectangle  (theContext, NULL);
    }
}

//__________________________________________________________________
void _HYPlatformGraphicPane::_InvertRect (_HYRect& rct)
{
    GdkRectangle clip_rect = HYRect2GDKRect (rct);
    GdkPixbuf* theImage = gdk_pixbuf_get_from_drawable (NULL, thePane, NULL, rct.left, rct.top, 0, 0, clip_rect.width, clip_rect.height);
    gdk_gc_set_function (theContext,GDK_INVERT);
    gdk_draw_pixbuf (thePane, theContext, theImage, 0,0, clip_rect.x, clip_rect.y, clip_rect.width, clip_rect.height, GDK_RGB_DITHER_NORMAL, 0, 0);
    gdk_gc_set_function (theContext,GDK_COPY);
    g_object_unref (theImage);
}


//__________________________________________________________________

void _HYPlatformGraphicPane::_DrawRect (_HYRect rectDesc)
{
    gdk_gc_set_line_attributes (theContext, rectDesc.width, GDK_LINE_SOLID, GDK_CAP_NOT_LAST, GDK_JOIN_MITER);
    gdk_draw_rectangle         (thePane, theContext, false, rectDesc.left, rectDesc.top,
                                rectDesc.right-rectDesc.left,rectDesc.bottom-rectDesc.top);
}


//__________________________________________________________________

void _HYPlatformGraphicPane::_FillRect (_HYRect rectDesc)
{
    gdk_gc_set_line_attributes (theContext, rectDesc.width, GDK_LINE_SOLID, GDK_CAP_NOT_LAST, GDK_JOIN_MITER);
    gdk_draw_rectangle         (thePane, theContext, true, rectDesc.left, rectDesc.top,
                                rectDesc.right-rectDesc.left+1,rectDesc.bottom-rectDesc.top+1);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_EraseRect (_HYRect rectDesc)
{
    gdk_gc_set_foreground (theContext,&saveBG);
    _FillRect (rectDesc);
    gdk_gc_set_foreground (theContext,&saveFG);
}


//__________________________________________________________________

void _HYPlatformGraphicPane::_DrawOval (_HYRect rectDesc)
{
    gdk_gc_set_line_attributes (theContext, rectDesc.width, GDK_LINE_SOLID, GDK_CAP_NOT_LAST, GDK_JOIN_MITER);
    gdk_draw_arc   (thePane, theContext, false, rectDesc.left, rectDesc.top, rectDesc.right-rectDesc.left+1,
                    rectDesc.bottom-rectDesc.top+1, 0, 64*360);
}


//__________________________________________________________________

void _HYPlatformGraphicPane::_FillOval (_HYRect rectDesc)
{
    gdk_gc_set_line_attributes (theContext, rectDesc.width, GDK_LINE_SOLID, GDK_CAP_NOT_LAST, GDK_JOIN_MITER);
    gdk_draw_arc   (thePane, theContext, true, rectDesc.left, rectDesc.top, rectDesc.right-rectDesc.left+1,
                    rectDesc.bottom-rectDesc.top+1, 0, 64*360);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_EraseOval (_HYRect rectDesc)
{
    gdk_gc_set_foreground (theContext,&saveBG);
    _FillOval (rectDesc);
    gdk_gc_set_foreground (theContext,&saveFG);
}


//__________________________________________________________________

void _HYPlatformGraphicPane::_DrawArc (_HYRect rectDesc, int s, int f)
{
    gdk_gc_set_line_attributes (theContext, rectDesc.width, GDK_LINE_SOLID, GDK_CAP_NOT_LAST, GDK_JOIN_MITER);
    gdk_draw_arc   (thePane, theContext, false, rectDesc.left, rectDesc.top, rectDesc.right-rectDesc.left+1,
                    rectDesc.bottom-rectDesc.top+1, (90-s)*64, -64*f);
}
//__________________________________________________________________

void _HYPlatformGraphicPane::_FillArc (_HYRect rectDesc, int s, int f)
{
    gdk_gc_set_line_attributes (theContext, rectDesc.width, GDK_LINE_SOLID, GDK_CAP_NOT_LAST, GDK_JOIN_MITER);
    gdk_draw_arc   (thePane, theContext, true, rectDesc.left, rectDesc.top, rectDesc.right-rectDesc.left+1,
                    rectDesc.bottom-rectDesc.top+1, (90-s)*64, -64*f);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_EraseArc (_HYRect rectDesc, int s, int f)
{
    gdk_gc_set_foreground (theContext,&saveBG);
    _FillArc (rectDesc,s,f);
    gdk_gc_set_foreground (theContext,&saveFG);
}


//__________________________________________________________________

void _HYPlatformGraphicPane::_SetColor  (_HYColor c)
{
    saveFG = HYColorToGDKColor(c);
    gdk_gc_set_foreground (theContext, &saveFG);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_SetBColor  (_HYColor c)
{
    saveBG = HYColorToGDKColor(c);
    gdk_gc_set_background (theContext, &saveBG);
}


//__________________________________________________________________

void _HYPlatformGraphicPane::_SetFont (_HYFont f)
{
    HYFont2PangoFontDesc(f,theFont);
    pango_layout_set_font_description (textLayout, theFont ); // ref ?
    _ResetCharGlyphs ();
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_SetFontSize (long s)
{
    pango_font_description_set_size   (theFont, s*PANGO_SCALE);
    pango_layout_set_font_description (textLayout, theFont ); // ref ?
    _ResetCharGlyphs ();
}



//__________________________________________________________________


void _HYPlatformGraphicPane::_DrawPicRes (_HYRect& r, long id)
{

    GdkPixbuf* theImage = (GdkPixbuf*)ProcureIconResource(id);
    if (theImage) {

        if (r.right-r.left<=0) {
            r.right = r.left + gdk_pixbuf_get_width (theImage);
        }
        if (r.bottom-r.top<=0) {
            r.bottom = r.top + gdk_pixbuf_get_height (theImage);
        }

        GdkRectangle clip_rect = HYRect2GDKRect (r);
        if (clip_rect.height != gdk_pixbuf_get_height (theImage) || clip_rect.width != gdk_pixbuf_get_width (theImage)) {
            GdkPixbuf*  scaledImage = gdk_pixbuf_scale_simple  (theImage, clip_rect.width, clip_rect.height, GDK_INTERP_BILINEAR);
            checkPointer (scaledImage);
            gdk_draw_pixbuf (thePane, NULL, scaledImage, 0,0, clip_rect.x, clip_rect.y, clip_rect.width, clip_rect.height, GDK_RGB_DITHER_NORMAL, 0, 0);
            g_object_unref (scaledImage);
        } else {
            //theImage = gdk_pixbuf_add_alpha(theImage,true,255,255,255);
            // add white transparency

            /*GdkPixbuf * copy = gdk_pixbuf_composite_color_simple (theImage, clip_rect.width, clip_rect.height,GDK_INTERP_BILINEAR,
                                                                            128,1,0x00ffffff,0x00888888);

            gdk_draw_pixbuf (thePane, NULL, copy, 0,0, clip_rect.x, clip_rect.y, clip_rect.width, clip_rect.height, GDK_RGB_DITHER_NORMAL, 0, 0);
            g_object_unref (copy);*/
            // inactive button

            /*GdkPixbuf * copy = gdk_pixbuf_copy(theImage);
            gdk_pixbuf_fill (copy, 0x000000ff);
            gdk_pixbuf_composite (theImage, copy, 0, 0, clip_rect.width, clip_rect.height, 0., 0., 1., 1.,GDK_INTERP_BILINEAR, 100);
            gdk_draw_pixbuf (thePane, NULL, copy, 0,0, clip_rect.x, clip_rect.y, clip_rect.width, clip_rect.height, GDK_RGB_DITHER_NORMAL, 0, 0);
            g_object_unref (copy);*/
            // pressed-in
            gdk_draw_pixbuf (thePane, NULL, theImage, 0,0, clip_rect.x, clip_rect.y, clip_rect.width, clip_rect.height, GDK_RGB_DITHER_NORMAL, 0, 0);
        }
    }
}


//__________________________________________________________________

void _HYPlatformGraphicPane::_SetDialogBG (void)
{
    ((_HYGraphicPane*)this)->SetBColor(GetDialogBackgroundColor());
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_StartDraw  (void)
{
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_EndDraw    (void)
{
}


//__________________________________________________________________

void _HYPlatformGraphicPane::_SetPort    (Ptr)
{
}



//__________________________________________________________________

void _HYPlatformGraphicPane::_CopyToClipboard   (void)
{
    // TBI
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_SavePicture   (_String prompt)
{
    if (exportFormats.lLength==0) {
        findGraphicsExporterComponents (exportFormats);
    }
    _String s1 ("snapshot"),
            s2 ("File Format"),
            filePath;

    long menuChoice = SaveFileWithPopUp (filePath, prompt,s1, s2,exportFormats);

    if (menuChoice >= 0) {
        GdkPixbuf* theImage = gdk_pixbuf_get_from_drawable (NULL, thePane, NULL, 0, 0, 0, 0, -1, -1);
        gdk_pixbuf_save (theImage, filePath.sData,((_String*)exportFormats (menuChoice))->sData,NULL,NULL);
        g_object_unref (theImage);
    }
}

//__________________________________________________________________

struct  _HYGDK_Polygon {

    GdkPoint * thePoints;
    gint       pointCount;

};

//__________________________________________________________________

Ptr _HYPlatformGraphicPane::_DefinePolygon  (_SimpleList& points)
{
    if ((points.lLength>=6)&&(points.lLength%2==0)) {
        _HYGDK_Polygon * thePoly = new _HYGDK_Polygon;
        thePoly->pointCount =   points.lLength/2;
        thePoly->thePoints = new GdkPoint [thePoly->pointCount];
        for (long k=0; k<points.lLength; k+=2) {
            thePoly->thePoints[k/2].x = points.lData[k];
            thePoly->thePoints[k/2].y = points.lData[k+1];
        }
        return (Ptr)thePoly;
    }
    return nil;
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_KillPolygon   (Ptr rgn)
{
    if (rgn) {
        _HYGDK_Polygon * thePoly = (_HYGDK_Polygon *)rgn;
        delete thePoly->thePoints;
        delete thePoly;
    }
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_DrawPolygon (Ptr rgn, long width)
{
    if (rgn) {
        gdk_gc_set_line_attributes (theContext, width, GDK_LINE_SOLID, GDK_CAP_NOT_LAST, GDK_JOIN_MITER);
        _HYGDK_Polygon * thePoly = (_HYGDK_Polygon *)rgn;
        gdk_draw_polygon (thePane, theContext, false, thePoly->thePoints, thePoly->pointCount);
    }
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_FillPolygon (Ptr rgn)
{
    if (rgn) {
        _HYGDK_Polygon * thePoly = (_HYGDK_Polygon *)rgn;
        gdk_draw_polygon (thePane, theContext, true, thePoly->thePoints, thePoly->pointCount);
    }
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_ErasePolygon (Ptr rgn)
{
    if (rgn) {
        gdk_gc_set_foreground (theContext,&saveBG);
        _FillPolygon (rgn);
        gdk_gc_set_foreground (theContext,&saveFG);
    }
}



//EOF
