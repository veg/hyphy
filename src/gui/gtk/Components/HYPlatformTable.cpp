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

#include "errorfns.h"
#include "HYTableComponent.h"
#include "HYUtils.h"
#include "HYEventTypes.h"
#include "HYWindow.h"
#include "HYGraphicPane.h"
#include "HYTextBox.h"
#include "HYTableWindow.h"
#include "HYDialogs.h"

#include <gdk-pixbuf/gdk-pixbuf.h>
#include <gdk/gdkkeysyms.h>

//__________________________________________________________________

extern  GdkCursor               *hSizeCursor,
        *pickUpCursor,
        *dropOffCursor;

GdkColor    menuLine1           = HYColorToGDKColor ((_HYColor)
{
    0xA0,0xA0,0xA0
}),
                                  menuLine2         = HYColorToGDKColor ((_HYColor)
{
    0x04,0x04,0x04
});

GdkColor  _BLACKBRUSH_          = HYColorToGDKColor ((_HYColor)
{
    0,0,0
});

extern     GdkPixbuf*           tablePDMenuIcon;

#define    HY_TABLE_PREEDIT_CURSOR  0x08

//__________________________________________________________________

gboolean   hy_table_key_interceptor (GtkWidget *te, GdkEventKey *kp, gpointer tbxp)
{
    _HYTable * theParent = (_HYTable*)tbxp;

    if (kp->keyval==GDK_KP_Enter || kp->keyval==GDK_Return) {
        theParent->EditBoxHandler (-1,theParent->rel);
        return true;
    }

    return FALSE;
}


//__________________________________________________________________

_HYPlatformTable::_HYPlatformTable(void)
{
    backPattern  = _BLACKBRUSH_;
    backPattern2 = _BLACKBRUSH_;
    cursorState  = false;
    editBox      = nil;
    tableFont    = nil;
    tableFontB   = nil;
    tableFontI   = nil;
    tableFontBI  = nil;
    activeColumn = -1;
    activeColumn2= -1;
    theContext   = nil;
}

//__________________________________________________________________

_HYPlatformTable::~_HYPlatformTable(void)
{
    if (tableFont) {
        pango_font_description_free (tableFont);
    }
    if (tableFontB) {
        pango_font_description_free (tableFontB);
    }
    if (tableFontI) {
        pango_font_description_free (tableFontI);
    }
    if (tableFontBI) {
        pango_font_description_free (tableFontBI);
    }
    if (theContext) {
        g_object_unref (theContext);
    }
}

//__________________________________________________________________

void        _HYPlatformTable::_SetFont (void)
{
    _HYTable* parent = (_HYTable*)this;

    if (tableFont) {
        pango_font_description_free (tableFont);
    }
    if (tableFontB) {
        pango_font_description_free (tableFontB);
    }
    if (tableFontI) {
        pango_font_description_free (tableFontI);
    }
    if (tableFontBI) {
        pango_font_description_free (tableFontBI);
    }

    _HYFont       tf = parent->textFont;
    tf.style = HY_FONT_PLAIN;
    tableFont   = pango_font_description_new();
    HYFont2PangoFontDesc(tf,tableFont);
    tf.style = HY_FONT_BOLD;
    tableFontB  = pango_font_description_new();
    HYFont2PangoFontDesc(tf,tableFontB);
    tf.style = HY_FONT_ITALIC;
    tableFontI  = pango_font_description_new();
    HYFont2PangoFontDesc(tf,tableFontI);
    tf.style = HY_FONT_ITALIC | HY_FONT_BOLD;
    tableFontBI = pango_font_description_new();
    HYFont2PangoFontDesc(tf,tableFontBI);
}

//__________________________________________________________________

bool        _HYPlatformTable::_CheckGC (void)
{
    if (!theContext) {
        _HYTable* theParent = (_HYTable*) this;

        if (theParent->parentWindow->window) {
            theContext           = gdk_gc_new (theParent->parentWindow->window);
        }

        return theContext != nil;

    }
    return true;
}

//__________________________________________________________________
void        _HYTable::_HScrollTable (long h)
{
    if (h && _CheckGC()) {
        _HYTable*       theParent = (_HYTable*) this;
        long            vsShift = ((settings.width&HY_COMPONENT_H_SCROLL)?HY_SCROLLER_WIDTH:0),
                        hsShift = ((settings.width&HY_COMPONENT_V_SCROLL)?HY_SCROLLER_WIDTH:0);

        EditBoxHandler  (-1,rel);

        _HYRect         paintRect = rel;

        //if (abs(h) > (rel.right-rel.left)/2)
        _Paint((Ptr)&paintRect);
        /*else
        {
            _HYRect       scrollRect;
            scrollRect.top      = rel.top;
            scrollRect.bottom   = rel.bottom-vsShift;
            scrollRect.right    = rel.right-hsShift;
            scrollRect.left     = rel.left;
            paintRect.top       = scrollRect.top;
            paintRect.bottom    = rel.bottom;

            GdkRectangle clip_rect  = HYRect2GDKRect (scrollRect);
            clip_rect.x += theParent->parentWindow->allocation.x;
            clip_rect.y += theParent->parentWindow->allocation.y;

            gdk_gc_set_clip_rectangle  (theContext, &clip_rect);
            GdkPixbuf* theImage = gdk_pixbuf_get_from_drawable (NULL, theParent->parentWindow->window, NULL,
                                                                clip_rect.x,
                                                                clip_rect.y,
                                                                0, 0, clip_rect.width-1, clip_rect.height-1);

            if (h>0)
            {
                gdk_draw_pixbuf (theParent->parentWindow->window, theContext, theImage, 0,0,
                                clip_rect.x-h,
                                clip_rect.y,
                                -1,-1, GDK_RGB_DITHER_NONE, 0, 0);

                paintRect.right =   scrollRect.right+hsShift;
                paintRect.left  =   scrollRect.right-h;

                hOrigin         +=  paintRect.left-rel.left;
                Paint((Ptr)&paintRect);
                hOrigin         -=  paintRect.left-rel.left;

            }
            else
            {
                gdk_draw_pixbuf (theParent->parentWindow->window, theContext, theImage, 0,0,
                                clip_rect.x-h,
                                clip_rect.y,
                                -1,-1, GDK_RGB_DITHER_NONE, 0, 0);

                paintRect.left  = scrollRect.left;
                paintRect.right = paintRect.left-h+hsShift;
                Paint   ((Ptr)&paintRect);
            }

            g_object_unref (theImage);
            gdk_gc_set_clip_rectangle  (theContext, NULL);
        }*/
    }
}


//__________________________________________________________________
void        _HYTable::_VScrollTable (long v)
{
    if (v && _CheckGC()) {
        long        vsShift = ((settings.width&HY_COMPONENT_H_SCROLL)?HY_SCROLLER_WIDTH:0),
                    hsShift = ((settings.width&HY_COMPONENT_V_SCROLL)?HY_SCROLLER_WIDTH:0);

        EditBoxHandler  (-1,rel);

        _HYTable*       theParent = (_HYTable*) this;
        _HYRect         paintRect = rel;

        //if (abs(v)>(rel.bottom-rel.top)/2)
        _Paint((Ptr)&paintRect);
        /*else
        {
            _HYRect                   scrollRect;

            scrollRect.left         = rel.left;
            scrollRect.right        = rel.right-hsShift;
            scrollRect.top          = rel.top;
            scrollRect.bottom       = rel.bottom-vsShift;
            paintRect.left          = scrollRect.left;
            paintRect.right         = rel.right;

            if (v<0)
                scrollRect.top--;
            GdkRectangle clip_rect  = HYRect2GDKRect (scrollRect);
            clip_rect.x += theParent->parentWindow->allocation.x;
            clip_rect.y += theParent->parentWindow->allocation.y;

            gdk_gc_set_clip_rectangle  (theContext, &clip_rect);
            GdkPixbuf* theImage = gdk_pixbuf_get_from_drawable (NULL, theParent->parentWindow->window, NULL,
                                                                clip_rect.x,
                                                                clip_rect.y,
                                                                0, 0, clip_rect.width - 1, clip_rect.height - 1);

            if (v>0)
            {
                gdk_draw_pixbuf (theParent->parentWindow->window, theContext, theImage, 0,0,
                                clip_rect.x,
                                clip_rect.y-v,
                                -1,-1, GDK_RGB_DITHER_NONE, 0, 0);
                paintRect.top = rel.bottom-vsShift-v-1;
                paintRect.bottom = rel.bottom-1;
                vOrigin+=paintRect.top-rel.top;
                Paint((Ptr)&paintRect);
                vOrigin-=paintRect.top-rel.top;
            }
            else
            {
                gdk_draw_pixbuf (theParent->parentWindow->window, theContext, theImage, 0,0,
                                clip_rect.x,
                                clip_rect.y-v,
                                -1,-1, GDK_RGB_DITHER_NONE, 0, 0);
                paintRect.top = rel.top;
                paintRect.bottom = rel.top-v+vsShift+2;
                paintRect.bottom = rel.bottom;
                Paint((Ptr)&paintRect);
            }

            g_object_unref (theImage);
            gdk_gc_set_clip_rectangle  (theContext, NULL);
        }*/
    }
}

//__________________________________________________________________

void        _HYTable::_ComponentMouseExit (void)
{
    if (cursorState) {
        if ((cursorState == HY_TABLE_DRAG_CURSOR)&&(activeColumn>=0)&&(activeColumn2>=0)) {
            _HiliteRowForDrag (activeColumn2,activeColumn);
        }

        _ResetCursorState ();
    }
}

//__________________________________________________________________

void        _HYPlatformTable::_SetBackColor (_HYColor& c)
{
    backPattern = HYColorToGDKColor(c);
}

//__________________________________________________________________

void        _HYPlatformTable::_SetBackColor2 (_HYColor& c)
{
    backPattern2 = HYColorToGDKColor(c);
}


//__________________________________________________________________

void        _HYPlatformTable::_KillTextBox (void)
{
    if (editBox) {
        gtk_widget_destroy(editBox);
        editBox = nil;
        _ResetCursorState ();
    }
}

//__________________________________________________________________

_String     _HYPlatformTable::_RetrieveTextValue (void)
{
    if (editBox) {
        return _String((char*)gtk_entry_get_text(GTK_ENTRY(editBox)));
    }

    return empty;
}

//__________________________________________________________________

void        _HYPlatformTable::_CreateTextBox (_HYRect& tBox,_String& textIn)
{
    textBoxRect = HYRect2GDKRect(tBox);

    _HYTable * theParent = (_HYTable*)this;

    editBox = gtk_entry_new ();
    checkPointer   (editBox);

    gtk_container_add(GTK_CONTAINER(theParent->parentWindow),editBox);

    g_signal_connect (G_OBJECT (editBox), "key-press-event", G_CALLBACK (hy_table_key_interceptor), theParent);

    if (tableFont) {
        gtk_widget_modify_font (editBox, tableFont);
    }

    gtk_fixed_move (GTK_FIXED(theParent->parentWindow), editBox, textBoxRect.x, textBoxRect.y);
    gtk_widget_set_size_request(editBox,textBoxRect.width,textBoxRect.height);
    gtk_entry_set_has_frame (GTK_ENTRY(editBox),true);
    gtk_entry_set_text(GTK_ENTRY(editBox), textIn.sData);
    gtk_widget_show(editBox);
    gtk_widget_activate(editBox);
    gtk_widget_grab_focus (editBox);
}



//__________________________________________________________________

GdkRectangle    _HYPlatformTable::_GetVisibleRowRect (long h)
{
    _HYTable*   parent = (_HYTable*)this;
    _HYRect res;

    long        w = (parent->settings.width&HY_COMPONENT_H_SCROLL)?HY_SCROLLER_WIDTH:0;
    res.left   = parent->rel.left;
    res.right  = (parent->settings.width&HY_COMPONENT_H_SCROLL)?parent->rel.right-HY_SCROLLER_WIDTH:parent->rel.right;

    res.bottom = parent->verticalSpaces.lData[h]-parent->vOrigin+parent->rel.top;
    if (res.bottom>parent->rel.bottom-w) {
        res.bottom=parent->rel.bottom-w;
    }

    if (h) {
        res.top = parent->verticalSpaces.lData[h-1]-parent->vOrigin+parent->rel.top;
    } else {
        res.top = parent->rel.top-parent->vOrigin;
    }

    return HYRect2GDKRect(res);
}


//__________________________________________________________________

void        _HYPlatformTable::_HiliteRowForDrag (long row, long old)
{
    if (_CheckGC()) {
        _HYTable * theParent = (_HYTable*)this;
        GdkRectangle        cellRect = _GetVisibleRowRect (row);

        if (row>=old) {
            cellRect.y = cellRect.y+cellRect.height-2;
        }

        cellRect.height = 2;

        GdkPixbuf* theImage = gdk_pixbuf_get_from_drawable (NULL, theParent->parentWindow->window, NULL, cellRect.x+theParent->parentWindow->allocation.x,
                              cellRect.y+theParent->parentWindow->allocation.y, 0, 0,
                              cellRect.width, cellRect.height);
        gdk_gc_set_function (theContext,GDK_INVERT);
        gdk_draw_pixbuf (theParent->parentWindow->window, theContext, theImage, 0,0,
                         cellRect.x+theParent->parentWindow->allocation.x,
                         cellRect.y+theParent->parentWindow->allocation.y,
                         cellRect.width, cellRect.height, GDK_RGB_DITHER_NORMAL, 0, 0);

        gdk_gc_set_function (theContext,GDK_COPY);
        g_object_unref (theImage);
    }
}

//__________________________________________________________________

void        _HYTable::_MarkCellsForUpdate (_SimpleList& cells)
{
    if (_CheckGC()) {
        long hs,hf,vs,vf,t,t2,k;
        GetDisplayRange (&rel, hs, hf, vs, vf);

        GdkRectangle  clipRect;

        clipRect.x     = rel.left;
        clipRect.width = rel.right - rel.left + 1;
        if (settings.width&HY_COMPONENT_V_SCROLL) {
            clipRect.width -= HY_SCROLLER_WIDTH;
        }

        clipRect.y      = rel.top;
        clipRect.height = rel.bottom - rel.top + 1;

        if (settings.width&HY_COMPONENT_H_SCROLL) {
            clipRect.height -= HY_SCROLLER_WIDTH;
        }

        for (k=0; k<cells.lLength; k++) {
            t2 = cells.lData[k]/horizontalSpaces.lLength;
            t  = cells.lData[k]%horizontalSpaces.lLength;

            if (t>=hs && t<=hf && t2<=vf && t2>=vs) {
                GdkRectangle invalRect;

                if (t) {
                    invalRect.x = horizontalSpaces.lData[t-1];
                } else {
                    invalRect.x = 0;
                }

                invalRect.width = horizontalSpaces.lData[t] - invalRect.x + 1;

                if (t2) {
                    invalRect.y = verticalSpaces.lData[t2-1];
                } else {
                    invalRect.y = 0;
                }

                invalRect.height = verticalSpaces.lData[t2] - invalRect.y + 1;
                invalRect.x += rel.left-hOrigin;
                invalRect.y += rel.top -vOrigin;

                GdkRectangle           tempRect;
                if (gdk_rectangle_intersect (&clipRect, &invalRect,&tempRect)) {
                    _HYTable * theParent = (_HYTable*) this;
                    tempRect.x += theParent->parentWindow->allocation.x;
                    tempRect.y += theParent->parentWindow->allocation.y;
                    gdk_window_invalidate_rect (theParent->parentWindow->window, &tempRect, false);
                }
            }
        }
    }
}

//__________________________________________________________________

void        _HYTable::_MarkColumnForUpdate (long index)
{
    if (_CheckGC()) {
        long hs,hf,vs,vf;
        GetDisplayRange (&rel, hs, hf, vs, vf);

        if (index>=hs && index<=hf) {
            GdkRectangle    clipRect;

            clipRect.x     = rel.left;
            clipRect.width = rel.right - rel.left + 1;
            if (settings.width&HY_COMPONENT_V_SCROLL) {
                clipRect.width -= HY_SCROLLER_WIDTH;
            }

            clipRect.y      = rel.top;
            clipRect.height = rel.bottom - rel.top + 1;

            if (settings.width&HY_COMPONENT_H_SCROLL) {
                clipRect.height -= HY_SCROLLER_WIDTH;
            }

            GdkRectangle    invalRect;

            invalRect.height = clipRect.height;
            invalRect.y      = clipRect.y;
            invalRect.x      = index?horizontalSpaces.lData[index-1]:0;
            invalRect.width  = horizontalSpaces.lData[index] - invalRect.x + 1;

            invalRect.x     += rel.left-hOrigin;

            GdkRectangle           tempRect;
            if (gdk_rectangle_intersect (&clipRect, &invalRect,&tempRect)) {
                _HYTable * theParent = (_HYTable*) this;
                tempRect.x += theParent->parentWindow->allocation.x;
                tempRect.y += theParent->parentWindow->allocation.y;
                gdk_window_invalidate_rect (theParent->parentWindow->window, &tempRect, false);
            }
        }
    }
}

//__________________________________________________________________

void        _HYTable::_MarkRowForUpdate (long index)
{
    if (_CheckGC()) {
        long hs,hf,vs,vf;
        GetDisplayRange (&rel, hs, hf, vs, vf);

        if ((index>=vs)&&(index<=vf)) {
            GdkRectangle  clipRect;

            clipRect.x      = rel.left;
            clipRect.width  = rel.right - rel.left + 1;
            if (settings.width&HY_COMPONENT_V_SCROLL) {
                clipRect.width -= HY_SCROLLER_WIDTH;
            }

            clipRect.y      = rel.top;
            clipRect.height = rel.bottom - rel.top + 1;
            if (settings.width&HY_COMPONENT_H_SCROLL) {
                clipRect.height -= HY_SCROLLER_WIDTH;
            }

            GdkRectangle      invalRect;
            invalRect.width     = clipRect.width;
            invalRect.x         = clipRect.x;
            invalRect.y         = index?verticalSpaces.lData[index-1]:0;
            invalRect.height    = verticalSpaces.lData[index]-invalRect.y+1;

            invalRect.y += rel.top-vOrigin;

            GdkRectangle           tempRect;
            if (gdk_rectangle_intersect (&clipRect, &invalRect,&tempRect)) {
                _HYTable * theParent = (_HYTable*) this;
                tempRect.x += theParent->parentWindow->allocation.x;
                tempRect.y += theParent->parentWindow->allocation.y;
                gdk_window_invalidate_rect (theParent->parentWindow->window, &tempRect, false);
            }
        }
    }
}


//__________________________________________________________________

void        _HYTable::_MarkCellForUpdate (long index)
{
    _SimpleList     dummy (index);
    _MarkCellsForUpdate (dummy);
}

//__________________________________________________________________

void        _HYTable::_IdleHandler (void)
{
}


//__________________________________________________________________

void        _HYTable::_FocusComponent (void)
{
    if (!GTK_WIDGET_HAS_FOCUS (parentWindow)) {
        gtk_widget_grab_focus (parentWindow);
    }
}


//__________________________________________________________________
long        _HYTable::_HandlePullDown (_List& data, long h, long v, long currentS)
{
    if (data.lLength) {
        return HandlePullDownWithFont (data,h,v,currentS,textFont.face,textFont.size);
    }

    return -1;
}


//__________________________________________________________________

void        _HYTable::_ScrollVPixels (long offset)
{
    long     voff = ((settings.width&HY_COMPONENT_H_SCROLL)?HY_SCROLLER_WIDTH:0);
    offset = offset/(_Parameter)(GetMaxH()- rel.bottom+rel.top+1-voff)*MAX_CONTROL_VALUE;
    ProcessEvent (generateScrollEvent(0,offset));
    _SetVScrollerPos((double)MAX_CONTROL_VALUE*vOrigin/(verticalSpaces.lData[verticalSpaces.lLength-1]-vSize+voff));
}

//__________________________________________________________________

void        _HYTable::_ScrollHPixels (long offset)
{
    long     hoff = ((settings.width&HY_COMPONENT_V_SCROLL)?HY_SCROLLER_WIDTH:0);
    offset = offset/(_Parameter)(GetMaxW()- rel.right+rel.bottom+1-hoff)*MAX_CONTROL_VALUE;
    ProcessEvent (generateScrollEvent(offset,0));
    _SetHScrollerPos((double)MAX_CONTROL_VALUE*hOrigin/(horizontalSpaces.lData[horizontalSpaces.lLength-1]-hSize+hoff));
}


//__________________________________________________________________

void        _HYTable::_ScrollRowIntoView (long index)
{
    if ((index>=0)&&(index<verticalSpaces.lLength)) {
        long hs, hf, vs, vf;
        GetDisplayRange (&rel,hs,hf,vs,vf);
        if ((index>vf)||(index<vs)) {
            _ScrollVPixels ((index?verticalSpaces.lData[index-1]:0)-vOrigin);
        }
    }
}

//__________________________________________________________________

void        _HYPlatformTable::_ResetCursorState (void)
{
    if (_CheckGC()) {
        gdk_window_set_cursor (((_HYTable*)this)->parentWindow->window, NULL);
    }
    cursorState = 0;
}


//__________________________________________________________________

void        _HYPlatformTable::_FrameRect        (GdkRectangle& theRect)
{
    if (_CheckGC()) {
        GdkGCValues saveGCValues;
        gdk_gc_get_values (theContext, &saveGCValues);

        gdk_gc_set_line_attributes (theContext, 2, GDK_LINE_ON_OFF_DASH,GDK_CAP_NOT_LAST, GDK_JOIN_MITER);
        gdk_gc_set_function        (theContext,GDK_INVERT);
        gdk_gc_set_foreground      (theContext, &_BLACKBRUSH_);
        gdk_gc_set_clip_rectangle  (theContext, NULL);

        gdk_draw_rectangle (((_HYTable*)this)->parentWindow->window,theContext,false,
                            theRect.x+((_HYTable*)this)->parentWindow->allocation.x,
                            theRect.y+((_HYTable*)this)->parentWindow->allocation.y,
                            theRect.width, theRect.height);

        gdk_gc_set_values (theContext, &saveGCValues,
                           (GdkGCValuesMask)(GDK_GC_FOREGROUND|GDK_GC_FUNCTION|GDK_GC_LINE_WIDTH|GDK_GC_LINE_STYLE|GDK_GC_CAP_STYLE|GDK_GC_JOIN_STYLE));
    }
}


//__________________________________________________________________
void        _HYTable::_Paint (Ptr p)
{
    if (!_CheckGC()) {
        return;
    }

    _HYRect         *relRect    = (_HYRect*)p;

    bool            isPrinting = relRect->right<0;

    if (isPrinting) {
        relRect->right = -relRect->right;
    }

    GdkColor        whiteC          = HYColorToGDKColor ((_HYColor) {
        0xff,0xff,0xff
    }),
    fillColor       = HYColorToGDKColor ((_HYColor) {
        102,204,255
    }),
                                              fillTColor        = HYColorToGDKColor (textColor),
                                              tableBkColor    = HYColorToGDKColor (backColor),
                                              tableBkColor2   = HYColorToGDKColor (backColor2),
                                              textColorGDK   = HYColorToGDKColor (textColor);

    GdkGCValues     savedDCValues;
    GdkDrawable*    window2Paint = GDK_DRAWABLE      (parentWindow->window);
    gdk_gc_get_values (theContext, &savedDCValues);

    PangoLayout     *textLayout   = pango_layout_new (screenPContext);
    pango_layout_set_width (textLayout, -1);

    //HBRUSH            themeFill   = CreateSolidBrush (fillColor);
    //checkPointer  (themeFill);

    //HRGN          saveRgn     = CreateRectRgn    (0,0,1,1);
    //checkPointer  (saveRgn);

    long            hs, // starting column
                    hf, // ending column
                    vs, // starting row
                    vf, // ending row
                    k,  // loop index
                    t,  // aux variable
                    t2,
                    st; // a few more auxs

    bool            chop,
                    chopv;


    GetDisplayRange (relRect, hs, hf, vs, vf);

    long            vsShift = ((settings.width&HY_COMPONENT_H_SCROLL)?HY_SCROLLER_WIDTH:0),
                    hsShift = ((settings.width&HY_COMPONENT_V_SCROLL)?HY_SCROLLER_WIDTH:0);


    _HYRect         saveRel;
    GdkRectangle    bRect;

    GdkGC*          offScreenPtr = nil,
                    *           dc           = theContext;

    GdkPixmap*      offBitmap    = nil;

    long            alX = parentWindow->allocation.x,
                    alY = parentWindow->allocation.y;

    //bRect = HYRect2GDKRect(*relRect);
    //gdk_gc_set_foreground(dc,&fillColor);
    //gdk_draw_rectangle(window2Paint,dc,true,bRect.x+alX,bRect.y+alY,bRect.width,bRect.width);
    //return;

    if ((vf-vs>2)&&(!isPrinting)) {
        offScreenPtr = gdk_gc_new (window2Paint);
        if (offScreenPtr) {
            bRect.x          = bRect.y = 0;
            bRect.width      = relRect->right-relRect->left-hsShift + 1;
            bRect.height     = relRect->bottom-relRect->top-vsShift + 1;
            offBitmap        = gdk_pixmap_new (window2Paint, bRect.width, bRect.height,-1);

            if (!offBitmap) {
                g_object_unref (offScreenPtr);
                offScreenPtr = nil;
            } else {
                dc = offScreenPtr;
                window2Paint = GDK_DRAWABLE (offBitmap);
                saveRel = *relRect;
                relRect->bottom -= relRect->top;
                relRect->top = 0;
                relRect->right  -= relRect->left;
                relRect->left = 0;
                alX = 0;
                alY = 0;
            }
        }
    }
    //printf ("Paint row %d-%d; parent offset %d %d\n",vs,vf, alX, alY);


    //SetTextAlign   (dc, TA_BASELINE);
    //SetBkMode    (dc, TRANSPARENT);

    GdkRectangle    clipRect = {
        relRect->left,relRect->top,
        relRect->right-hsShift - relRect->left + 1,
        relRect->bottom-vsShift - relRect->top + 1
    },
    anotherRect,
    clipRect2,
    clipRect3;

    //POINT       mapPts[2];
    //mapPts[0] = (POINT){clipRect.left,clipRect.top};
    //mapPts[1] = (POINT){clipRect.right+1,clipRect.bottom+1};
    //LPtoDP (dc, mapPts,2);

    anotherRect  = clipRect;

    anotherRect.x += alX;
    anotherRect.y += alY;
    //anotherRect.width  ++;
    //anotherRect.height ++;

    GdkRegion     *tempRgn = gdk_region_rectangle (&anotherRect);

    gdk_gc_set_clip_region(dc, tempRgn);

    //::SetTextColor (dc, RGB(textColor.R, textColor.G, textColor.B));

    anotherRect = clipRect;
    t = relRect->top-vOrigin;

    //saveDCPen = (HPEN)SelectObject (dc, whitePen);
    gdk_gc_set_line_attributes(dc,1,GDK_LINE_SOLID,GDK_CAP_NOT_LAST,GDK_JOIN_MITER);

    for (k=vs; k<=vf; k++) {
        anotherRect.y       = k?verticalSpaces.lData[k-1]+t:relRect->top;
        anotherRect.height  = verticalSpaces.lData[k]+t-anotherRect.y + 1;

        if (cellTypes.lData[k*horizontalSpaces.lLength]&HY_TABLE_BEVELED) {
            gdk_gc_set_foreground(dc,&backPattern2);
            gdk_draw_rectangle   (window2Paint,dc,true,anotherRect.x+alX,anotherRect.y+alY,anotherRect.width
                                  ,anotherRect.height);
            if ( k==vs || k<vf ) {
                gdk_gc_set_foreground(dc,&menuLine2);
                //MoveToEx (dc,anotherRect.left,anotherRect.bottom-1,nil);
                //LineTo (dc,anotherRect.right,anotherRect.bottom-1);
                /*gdk_draw_line (window2Paint,dc,anotherRect.x+alX,
                                               anotherRect.y+alY+anotherRect.height-2,
                                               anotherRect.x+alX+anotherRect.width-1,
                                               anotherRect.y+alY+anotherRect.height-2);*/
                gdk_gc_set_foreground(dc,&menuLine1);
                //MoveToEx (dc,anotherRect.left,anotherRect.bottom-2,nil);
                //LineTo (dc,anotherRect.right,anotherRect.bottom-2);
                gdk_draw_line (window2Paint,dc,anotherRect.x+alX,
                               anotherRect.y+alY+anotherRect.height-3,
                               anotherRect.x+alX+anotherRect.width-1,
                               anotherRect.y+alY+anotherRect.height-3);
                //gdk_gc_set_foreground(dc,&whiteC);
            }
        } else {
            gdk_gc_set_foreground(dc,&backPattern);
            gdk_draw_rectangle   (window2Paint,dc,true,anotherRect.x+alX,anotherRect.y+alY,anotherRect.width
                                  ,anotherRect.height);
            if ( k==vs || k<vf ) {
                gdk_gc_set_foreground(dc,&whiteC);
                //MoveToEx (dc,anotherRect.left,anotherRect.bottom-1,nil);
                //LineTo (dc,anotherRect.right,anotherRect.bottom-1);
                gdk_draw_line (window2Paint,dc,anotherRect.x+alX,
                               anotherRect.y+alY+anotherRect.height-2,
                               anotherRect.x+alX+anotherRect.width-1,
                               anotherRect.y+alY+anotherRect.height-2);
            }
        }
    }


    st = 0;
    if (hf<horizontalSpaces.lLength-1) {
        st = hf;
    } else {
        st = hf-1;
    }
    t = relRect->left-hOrigin-2;

    if ((selectionType & HY_TABLE_NO_COLS_LINES) == 0) {
        gdk_gc_set_foreground(dc,&menuLine1);
        for (k=hs; k<=st; k++) {
            t2 = t+horizontalSpaces.lData[k];
            gdk_draw_line (window2Paint,dc,t2+alX,
                           relRect->top+alY,
                           t2+alX,
                           relRect->bottom+alY+1);
            //MoveToEx (dc,t2,relRect->top,nil);
            //LineTo (dc,t2,relRect->bottom);
        }
        gdk_gc_set_foreground(dc,&menuLine2);
        t++;
        for (k=hs; k<=st; k++) {
            t2 = t+horizontalSpaces.lData[k];
            //MoveToEx (dc,t2,relRect->top,nil);
            //LineTo (dc,t2,relRect->bottom);
            gdk_draw_line (window2Paint,dc,t2+alX,
                           relRect->top+alY,
                           t2+alX,
                           relRect->bottom+alY+1);
        }
    }

    st = -1;
    bool  highlightColorOn = false;

    for (k=vs; k<=vf; k++) {
        anotherRect.y       = relRect->top-vOrigin+1;
        anotherRect.height  = verticalSpaces.lData[k];

        if (k) {
            anotherRect.y      += verticalSpaces.lData[k-1];
            anotherRect.height -= verticalSpaces.lData[k-1];
        }

        long    t3,
                st2,
                w = anotherRect.height-1,
                w2,
                shift = (w-textFont.size)/2-1;

        if (anotherRect.y+anotherRect.height -1 > relRect->bottom-vsShift) {
            anotherRect.height = relRect->bottom-vsShift - anotherRect.y + 1;
            chopv = false;
        } else {
            chopv = true;
        }

        for (t2=hs; t2<=hf; t2++) {
            t3 = k*horizontalSpaces.lLength+t2;
            if (t3 == editCellID) {
                continue;
            }

            anotherRect.x     = relRect->left-hOrigin+1;
            anotherRect.width = horizontalSpaces.lData[t2];

            if (t2) {
                anotherRect.x     += horizontalSpaces.lData[t2-1];
                anotherRect.width -= horizontalSpaces.lData[t2-1];
            }

            clipRect2 = anotherRect;
            w2        = anotherRect.width-1;

            chop = true;

            if (anotherRect.x + anotherRect.width - 1 > relRect->right-hsShift) {
                anotherRect.width = relRect->right-hsShift - anotherRect.x + 1;
                chop = false;
            } else {
                chop = true;
            }

            gdk_region_destroy (tempRgn);

            if ((t2==hs)||(k==vs)) {
                gdk_rectangle_intersect(&anotherRect,&clipRect,&clipRect3);
            } else {
                clipRect3 = anotherRect;
            }

            clipRect3.height ++;
            clipRect3.width ++;
            clipRect3.x += alX;
            clipRect3.y += alY;
            tempRgn = gdk_region_rectangle(&clipRect3);
            gdk_gc_set_clip_region (dc, tempRgn);


            if (cellTypes.lData[t3]&HY_TABLE_SELECTED) {
                gdk_gc_set_foreground(dc,&fillColor);
                gdk_draw_rectangle(window2Paint,dc,true,    anotherRect.x+alX,
                                   anotherRect.y+alY,
                                   anotherRect.width-(chop?2:0),
                                   anotherRect.height-chopv);
            }

            if (cellTypes.lData[t3]&HY_TABLE_ICON) {
                if (!isPrinting) {
                    _SimpleList     * cellList = (_SimpleList*)cellData.lData[t3];

                    if (w2-4>cellList->lData[1]) {
                        clipRect2.x += (w2-cellList->lData[1])/2;
                    }

                    clipRect2.width = cellList->lData[1]+1;

                    if (w-2>cellList->lData[2]) {
                        clipRect2.y += (w-cellList->lData[2])/2;
                    }

                    clipRect2.height = cellList->lData[2]+1;

                    if (cellList->lLength==3) {
                        GdkPixbuf* aPic = (GdkPixbuf*)cellList->lData[0];
                        if (aPic) {

                            if (clipRect2.height != gdk_pixbuf_get_height (aPic) || clipRect2.width != gdk_pixbuf_get_width (aPic)) {
                                GdkPixbuf*  scaledImage = gdk_pixbuf_scale_simple  (aPic, clipRect2.width, clipRect2.height, GDK_INTERP_BILINEAR);
                                checkPointer (scaledImage);
                                gdk_draw_pixbuf (window2Paint, dc, scaledImage, 0,0, clipRect2.x+alX, clipRect2.y+alY, clipRect2.width, clipRect2.height, GDK_RGB_DITHER_NORMAL, 0, 0);
                                g_object_unref (scaledImage);
                            } else {
                                gdk_draw_pixbuf (window2Paint, dc, aPic, 0,0, clipRect2.x+alX, clipRect2.y+alY, clipRect2.width, clipRect2.height, GDK_RGB_DITHER_NORMAL, 0, 0);
                            }
                        }
                    } else {
                        if ((cellList->lData[3]==HY_TABLE_COLOR_BOX)||(cellList->lData[3]==HY_TABLE_COLOR_CIRCLE)) {
                            _HYColor    c   = LongToHYColor    (cellList->lData[0]);

                            if (cellList->lData[3]==HY_TABLE_COLOR_BOX) {
                                GdkColor  clr = HYColorToGDKColor(c);
                                gdk_gc_set_foreground(dc,&clr);
                                gdk_draw_rectangle(window2Paint,dc,true,clipRect2.x+alX,clipRect2.y+alY,clipRect2.width,clipRect.height);
                            } else {
                                gdk_gc_set_line_attributes (dc, 1, GDK_LINE_SOLID, GDK_CAP_NOT_LAST, GDK_JOIN_MITER);



                                GdkColor                circColor = _BLACKBRUSH_;
                                gdk_gc_set_foreground   (dc,&circColor);
                                gdk_draw_arc   (window2Paint, dc, true, clipRect2.x-2+alX,
                                                clipRect2.y-2+alY,
                                                clipRect2.width+1,
                                                clipRect2.height+1,
                                                0, 64*360);


                                circColor = HYColorToGDKColor((_HYColor) {
                                    c.R/2,c.G/2,c.B/2
                                });

                                gdk_gc_set_foreground   (dc,&circColor);
                                gdk_draw_arc   (window2Paint, dc, true, clipRect2.x-1+alX,
                                                clipRect2.y-1+alY,
                                                clipRect2.width,
                                                clipRect2.height,
                                                0, 64*360);

                                circColor = HYColorToGDKColor((_HYColor) {
                                    c.R/1.25,c.G/1.25,c.B/1.25
                                });

                                gdk_gc_set_foreground   (dc,&circColor);
                                gdk_draw_arc   (window2Paint, dc, true, clipRect2.x+alX,
                                                clipRect2.y+alY,
                                                clipRect2.width-2,
                                                clipRect2.height-2,
                                                0, 64*360);

                                circColor = HYColorToGDKColor(c);

                                gdk_gc_set_foreground   (dc,&circColor);
                                gdk_draw_arc   (window2Paint, dc, true, clipRect2.x+1+alX,
                                                clipRect2.y+1+alY,
                                                clipRect2.width-4,
                                                clipRect2.height-4,
                                                0, 64*360);
                            }
                        }
                    }
                }
            } else { // text
                st2 = cellTypes.lData[t3]&HY_TABLE_STYLEMASK;
                if (st!=st2) {
                    PangoFontDescription*   setFont = tableFont;
                    st = st2;
                    if (st&HY_TABLE_BOLD) {
                        if (st&HY_TABLE_ITALIC) {
                            setFont = tableFontBI;
                        } else {
                            setFont = tableFontB;
                        }
                    } else if (st&HY_TABLE_ITALIC) {
                        setFont = tableFontI;
                    }

                    pango_layout_set_font_description(textLayout,setFont);
                }


                _String  *thisCell = (_String*)cellData.lData[t3];

                if (cellTypes.lData[t3]&HY_TABLE_SELECTED) {
                    gdk_gc_set_foreground(dc,&fillTColor);
                } else {
                    gdk_gc_set_foreground(dc,&textColorGDK);
                }

                pango_layout_set_text(textLayout,thisCell->sData,thisCell->sLength);
                gdk_draw_layout(window2Paint,dc,anotherRect.x+textFont.size/3+alX,anotherRect.y+shift+alY,textLayout);
                //TextOut (dc,, anotherRect.top+shift+textFont.size, );
                if (cellTypes.lData[t3]&HY_TABLE_PULLDOWN) {
                    if (!isPrinting) {
                        clipRect2.x     = clipRect2.x + clipRect2.width - 4 - tPDMw;
                        clipRect2.width = tPDMw;

                        if (w-2>tPDMh) {
                            t3 = (w-tPDMh)/2;
                            clipRect2.y += t3;
                        }
                        clipRect2.height = tPDMh;
                        gdk_draw_pixbuf (window2Paint, dc, tablePDMenuIcon, 0,0, clipRect2.x+alX, clipRect2.y+alY, clipRect2.width, clipRect2.height, GDK_RGB_DITHER_NORMAL, 0, 0);
                    }
                }
            }
        }
    }


    gdk_gc_set_clip_region  (dc, NULL);

    if (offScreenPtr) {
        *relRect    = saveRel;
        clipRect.x += relRect->left + parentWindow->allocation.x;
        clipRect.y += relRect->top  + parentWindow->allocation.y;

        gdk_draw_drawable (GDK_DRAWABLE(parentWindow->window), theContext, offBitmap, 0, 0, clipRect.x, clipRect.y, -1, -1);
        g_object_unref   (offScreenPtr);
        g_object_unref    (offBitmap);
    }

    //gdk_gc_set_values (theContext, &savedDCValues, 0xff);
    gdk_region_destroy(tempRgn);
    g_object_unref (textLayout);
    if (editBox) {
        gtk_widget_queue_draw(editBox);
    }

    _HYPlatformComponent::_Paint(p);
}


//__________________________________________________________________

bool        _HYTable::_ProcessOSEvent (Ptr vEvent)
{
    static  int     lastH = 0,
                    lastV = 0;

    long            k,
                    h,
                    v,
                    vsShift = ((settings.width&HY_COMPONENT_H_SCROLL)?HY_SCROLLER_WIDTH:0),
                    hsShift = ((settings.width&HY_COMPONENT_V_SCROLL)?HY_SCROLLER_WIDTH:0);

    _HY_GTK_UI_Message *theMessage = (_HY_GTK_UI_Message*)vEvent;

    gdouble   xc,
              yc;

    gdk_event_get_coords (theMessage->theEvent,&xc,&yc);

    switch (theMessage->theEvent->type) {
    case GDK_BUTTON_PRESS:
    case GDK_2BUTTON_PRESS: {
        GdkEventButton * bevent = (GdkEventButton*)theMessage->theEvent;

        lastH = xc-parentWindow->allocation.x;
        lastV = yc-parentWindow->allocation.y;


        if ((selectionType&HY_TABLE_FOCUSABLE)&&messageRecipient&&((selectionType&HY_TABLE_IS_FOCUSED)==0)) {
            messageRecipient->ProcessEvent(generateKeyboardFocusEvent (GetID()));
        }

        if ( cursorState == HY_TABLE_SIZE_CURSOR && theMessage->theEvent->type != GDK_2BUTTON_PRESS ) {
            // do drag here
            EditBoxHandler (-1,rel);

            //SetCapture      (parentWindow);
            gdk_pointer_grab (parentWindow->window,false,
                              (GdkEventMask)(GDK_POINTER_MOTION_MASK | GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK),
                              parentWindow->window, NULL, bevent->time);

            limits.x      = lastH;
            limits.y      = rel.top;
            limits.height = rel.bottom-vsShift-limits.y + 1;
            limits.width  = rel.right -hsShift;

            lastH += hOrigin;

            for (activeColumn = 0; activeColumn<horizontalSpaces.lLength-1; activeColumn++)
                if (horizontalSpaces.lData[activeColumn]>lastH-2-rel.left) {
                    break;
                }

            if (activeColumn) {
                limits.x = rel.left+horizontalSpaces.lData[activeColumn-1]+3-hOrigin;
            } else {
                limits.x = rel.left+3;
            }

            long dragRes = horizontalSpaces.lData[horizontalSpaces.lLength-1]-
                           rel.right+rel.left-hOrigin+hsShift;

            if (dragRes<lastH-limits.x) {
                limits.x = lastH-dragRes;
            }

            lastH        -= hOrigin;
            limits.width += 1-limits.x;

            return true;

        } else {
            if (((k=FindClickedTableCell (lastH-rel.left,lastV-rel.top,h,v))>-1)&&(lastV<rel.bottom-vsShift)&&(lastH<rel.right-hsShift)) {
                if (theMessage->theEvent->type == GDK_2BUTTON_PRESS) {
                    if (cursorState == HY_TABLE_SIZE_CURSOR) {
                        _ResetCursorState ();
                        gdk_pointer_ungrab (((GdkEventButton*)theMessage->theEvent)->time);
                    } else if ((cursorState == HY_TABLE_DRAG_CURSOR || HY_TABLE_PREEDIT_CURSOR == cursorState)&& messageRecipient) {
                        gdk_pointer_ungrab (((GdkEventButton*)theMessage->theEvent)->time);
                        ((_HYTWindow*)messageRecipient)->trackMouseComponent = nil;
                    }

                    if (cellTypes.lData[k]&HY_TABLE_EDIT_TEXT) {
                        EditBoxHandler (k,rel);
                    } else if (messageRecipient) {
                        messageRecipient->ProcessEvent (generateTableDblClickEvent(GetID()));
                    }
                    break;
                }
                ModifySelection (h,v,bevent->state & GDK_SHIFT_MASK, bevent->state & GDK_CONTROL_MASK, true);
            }

            if (k==-2)
                // process pull-down
            {
                if (messageRecipient)
                    messageRecipient->ProcessEvent (generateTablePullDownEvent(GetID(),v*horizontalSpaces.lLength+h,
                                                    (((long)xc)<<16)+(long)yc));
                break;
            }


            if (cursorState == HY_TABLE_DRAG_CURSOR && theMessage->theEvent->type != GDK_2BUTTON_PRESS) {
                gdk_window_set_cursor (parentWindow->window, dropOffCursor);
                FindClickedTableCell(lastH-rel.left,lastV-rel.top,k,activeColumn);
                activeColumn2 = -1;
                if (messageRecipient) {
                    gdk_pointer_grab (parentWindow->window,false,
                                      (GdkEventMask)(GDK_POINTER_MOTION_MASK | GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK),
                                      NULL, NULL, bevent->time);
                    ((_HYTWindow*)messageRecipient)->trackMouseComponent = (Ptr)this;
                }
                return true;
            }


            /*  handle row drag here */
            if ((selectionType&HY_TABLE_SINGLE_SELECTION)==0 &&
                    (selectionType&HY_TABLE_NODRAG_SELECTION)==0 &&
                    lastH<rel.right-hsShift &&
                    lastV<rel.bottom-vsShift ) {
                cursorState = HY_TABLE_PREEDIT_CURSOR;
                limits.x     = lastH;
                limits.width = 0;
                limits.y     = lastV;
                limits.height = 0;

                textBoxRect.x      = -1;
                textBoxRect.width  = lastH;
                textBoxRect.height = lastV;

                if (messageRecipient) {
                    gdk_pointer_grab (parentWindow->window,false,
                                      (GdkEventMask)(GDK_POINTER_MOTION_MASK | GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK),
                                      NULL, NULL, bevent->time);

                    ((_HYTWindow*)messageRecipient)->trackMouseComponent = (Ptr)this;
                }

                return true;
            }
        }
    }
    break;

    case GDK_MOTION_NOTIFY: {
        GdkEventMotion * mevent = (GdkEventMotion*)theMessage->theEvent;

        GdkPoint        downWhere = (GdkPoint) {
            xc-parentWindow->allocation.x,yc-parentWindow->allocation.y
        };

        if (mevent->state & GDK_BUTTON1_MASK) { // left button down
            if (cursorState == HY_TABLE_SIZE_CURSOR) {
                GdkPoint   currentPoint = (GdkPoint) {
                    xc-parentWindow->allocation.x,yc-parentWindow->allocation.y
                };

                if (currentPoint.x >= limits.x && currentPoint.x < limits.x + limits.width &&
                        currentPoint.y >= limits.y && currentPoint.y < limits.y + limits.height ) {
                    if (currentPoint.x-lastH) {
                        SetColumnSpacing (activeColumn,currentPoint.x-lastH,true);
                        if (messageRecipient)
                            messageRecipient->ProcessEvent (generateTableResizeCEvent(GetID(),
                                                            activeColumn,currentPoint.x-lastH));
                        lastH = currentPoint.x;
                    }
                }
            } else if (cursorState == HY_TABLE_DRAG_CURSOR && activeColumn>=0) {
                long h,v,k;

                if ((downWhere.y>rel.bottom-vsShift)||(downWhere.y<rel.top)) {
                    if (activeColumn2>=0) {
                        _HiliteRowForDrag (activeColumn2,activeColumn);
                        activeColumn2 = -1;
                    }
                    h = verticalSpaces.lData[verticalSpaces.lLength-1]/verticalSpaces.lLength;
                    _ScrollVPixels ((downWhere.y<rel.top)?-h:h);
                    break;
                }

                if (downWhere.x>rel.right-hsShift) {
                    downWhere.x=rel.right-hsShift;
                }

                if ( (lastH!=downWhere.x)|| (lastV!=downWhere.y)) {
                    k = FindClickedTableCell(downWhere.x-rel.left,downWhere.y-rel.top,h,v);
                    if ((v!=activeColumn2)&&(k>-1)) {
                        if (activeColumn2>=0) {
                            _HiliteRowForDrag (activeColumn2,activeColumn);
                        }
                        if ((v!=activeColumn)&&(!(cellTypes.lData[k]&HY_TABLE_CANTSELECT))) {
                            _HiliteRowForDrag (v,activeColumn);
                            activeColumn2 = v;
                        } else {
                            activeColumn2 = -1;
                        }
                    }
                    lastH = downWhere.x;
                    lastV = downWhere.y;
                }
            } else if (cursorState == HY_TABLE_EDIT_CURSOR || cursorState == HY_TABLE_PREEDIT_CURSOR) {
                cursorState = HY_TABLE_EDIT_CURSOR;

                GdkRectangle      clippingRect = HYRect2GDKRect(rel),
                                  paintRect;

                while (1) {
                    if (downWhere.y>rel.bottom-vsShift) {
                        if (rel.bottom-rel.top-vsShift+vOrigin <  verticalSpaces.lData[verticalSpaces.lLength-1]-1) {
                            long h = verticalSpaces.lData[verticalSpaces.lLength-1]/verticalSpaces.lLength;

                            gdk_rectangle_intersect(&limits,&clippingRect,&paintRect);
                            _FrameRect (paintRect);
                            _ScrollVPixels (h);

                            limits.y     -= h;
                            limits.height += h;

                            gdk_rectangle_intersect(&limits,&clippingRect,&paintRect);
                            _FrameRect    (paintRect);

                            lastV               -= h;
                            textBoxRect.height  -= h;
                            break;
                        }
                        downWhere.y=rel.bottom-vsShift;
                    }
                    if (downWhere.x>rel.right-hsShift) {
                        if (rel.right-rel.left-hsShift+hOrigin <  horizontalSpaces.lData[horizontalSpaces.lLength-1]-1) {
                            long h = horizontalSpaces.lData[horizontalSpaces.lLength-1]/horizontalSpaces.lLength;

                            gdk_rectangle_intersect(&limits,&clippingRect,&paintRect);
                            _FrameRect    (paintRect);

                            _ScrollHPixels (h);

                            limits.x     -= h;
                            limits.width += h;

                            gdk_rectangle_intersect(&limits,&clippingRect,&paintRect);
                            _FrameRect    (paintRect);

                            lastH             -= h;
                            textBoxRect.width -= h;
                            break;
                        }
                        downWhere.x=rel.right-hsShift;
                    }
                    if (downWhere.y<rel.top) {
                        if (vOrigin>0) {
                            long h = verticalSpaces.lData[verticalSpaces.lLength-1]/verticalSpaces.lLength;

                            gdk_rectangle_intersect(&limits,&clippingRect,&paintRect);
                            _FrameRect    (paintRect);

                            _ScrollVPixels (-h);

                            limits.y      += h;
                            limits.height -= h;

                            gdk_rectangle_intersect(&limits,&clippingRect,&paintRect);
                            _FrameRect    (paintRect);

                            lastV              += h;
                            textBoxRect.height += h;
                            break;
                        }
                        downWhere.y=rel.top;
                    }
                    if (downWhere.x<rel.left) {
                        if (hOrigin>0) {
                            long h = horizontalSpaces.lData[horizontalSpaces.lLength-1]/horizontalSpaces.lLength;

                            gdk_rectangle_intersect(&limits,&clippingRect,&paintRect);
                            _FrameRect    (paintRect);
                            _ScrollHPixels (-h);

                            limits.x     += h;
                            limits.width -= h;

                            gdk_rectangle_intersect(&limits,&clippingRect,&paintRect);
                            _FrameRect    (paintRect);

                            lastH             += h;
                            textBoxRect.width += h;
                            break;
                        }
                        downWhere.x=rel.left;
                    }

                    if ( lastH!=downWhere.x || lastV!=downWhere.y) {
                        if (textBoxRect.x>=0) {
                            gdk_rectangle_intersect(&limits,&clippingRect,&paintRect);
                            _FrameRect    (paintRect);
                        }

                        if (downWhere.x > textBoxRect.width) {
                            limits.width = downWhere.x - limits.x + 1;
                        } else {
                            limits.x = downWhere.x;
                            limits.width =  textBoxRect.width - limits.x + 1;
                        }
                        if (downWhere.y >  textBoxRect.height) {
                            limits.height = downWhere.y - limits.y + 1;
                        } else {
                            limits.y = downWhere.y;
                            limits.height = textBoxRect.height - limits.y + 1;
                        }

                        gdk_rectangle_intersect(&limits,&clippingRect,&paintRect);
                        _FrameRect    (paintRect);

                        lastH = downWhere.x;
                        lastV = downWhere.y;

                        textBoxRect.x = 1;
                    }
                    break;
                }
                return true;
            }

        } else { // treat as plain move
            bool  ch = (!(selectionType&HY_TABLE_DONT_SIZE))&&(CheckForHSizeLocation(downWhere.x-rel.left))&&(downWhere.y<rel.bottom-vsShift);

            if (ch&&(cursorState!=HY_TABLE_SIZE_CURSOR)) {
                cursorState = HY_TABLE_SIZE_CURSOR;
                gdk_window_set_cursor(parentWindow->window,hSizeCursor);
            } else if ((!ch)&&(cursorState==HY_TABLE_SIZE_CURSOR)) {
                gdk_window_set_cursor(parentWindow->window,NULL);
                cursorState = 0;
            }

            if (selectionType & HY_TABLE_SEL_ROWS) {
                if (!ch) {
                    k = FindClickedTableCell(downWhere.x-rel.left,downWhere.y-rel.top,h,v);
                    if (k>=0) {
                        if ((cursorState != HY_TABLE_DRAG_CURSOR)&&
                                ((selectionType&HY_TABLE_NODRAG_SELECTION)==0)) {
                            if (cellTypes.lData[k]&HY_TABLE_SELECTED) {
                                if (IsRowSelectionSimple()) {
                                    gdk_window_set_cursor(parentWindow->window,pickUpCursor);
                                    cursorState = HY_TABLE_DRAG_CURSOR;
                                    activeColumn = -1;
                                }
                            }
                        } else {
                            if (((selectionType&HY_TABLE_NODRAG_SELECTION)==0)&&(!(cellTypes.lData[k]&HY_TABLE_SELECTED))) {
                                gdk_window_set_cursor(parentWindow->window,NULL);
                                cursorState = 0;
                            }
                        }
                    }
                }
            }
        }
        return true;
    }
    break;

    case GDK_KEY_PRESS: {
        GdkEventKey * kp = (GdkEventKey*)theMessage->theEvent;
        bool          ctlDown =  kp->state & GDK_CONTROL_MASK;

        switch (kp->keyval) {
        case GDK_Up:
        case GDK_KP_Up:
            HandleKeyMove (0,ctlDown);
            return true;
        case GDK_Down:
        case GDK_KP_Down:
            HandleKeyMove (1,ctlDown);
            return true;
        case GDK_Left:
        case GDK_KP_Left:
            HandleKeyMove (2,ctlDown);
            return true;
        case GDK_Right:
        case GDK_KP_Right:
            HandleKeyMove (3,ctlDown);
            return true;

        }
        break;
    }
    case GDK_BUTTON_RELEASE: {
        lastH = xc - parentWindow->allocation.x;
        lastV = yc - parentWindow->allocation.y;

        GdkEventButton * be = (GdkEventButton*)theMessage->theEvent;

        bool    isInComponent = lastH>=rel.left && lastV>=rel.top && lastH<rel.right && lastV<rel.bottom;

        if (cursorState == HY_TABLE_SIZE_CURSOR) {
            _ResetCursorState ();
            gdk_pointer_ungrab (be->time);
        } else if (cursorState == HY_TABLE_DRAG_CURSOR) {
            if (activeColumn >= 0) {
                if (messageRecipient) {
                    gdk_pointer_ungrab (be->time);
                    ((_HYTWindow*)messageRecipient)->trackMouseComponent = nil;
                }
                if (activeColumn2>=0) {
                    _HiliteRowForDrag (activeColumn2,activeColumn);
                }

                if (activeColumn!=activeColumn2) {
                    EditBoxHandler (-1,rel);
                    if (isInComponent) {
                        DragRow (activeColumn,activeColumn2);
                    }
                }

                activeColumn = -1;
            }
            _ResetCursorState ();
        } else if (cursorState == HY_TABLE_EDIT_CURSOR || cursorState == HY_TABLE_PREEDIT_CURSOR) {
            if (messageRecipient) {
                gdk_pointer_ungrab (be->time);
                ((_HYTWindow*)messageRecipient)->trackMouseComponent = nil;
            }

            if (cursorState == HY_TABLE_EDIT_CURSOR) {
                GdkRectangle      clippingRect = HYRect2GDKRect (rel),
                                  paintRect;

                gdk_rectangle_intersect(&limits,&clippingRect,&paintRect);
                _FrameRect    (paintRect);

                _HYRect     outlineHRect;
                outlineHRect.left   = limits.x;
                outlineHRect.right  = limits.x+limits.width-1+hsShift;
                outlineHRect.top    = limits.y;
                outlineHRect.bottom = limits.y+limits.height-1+vsShift;

                long    hs,hf,vs,vf;

                hOrigin += limits.x-rel.left;
                vOrigin += limits.y-rel.top;
                GetDisplayRange (&outlineHRect,hs,hf,vs,vf);
                hOrigin -= limits.x-rel.left;
                vOrigin -= limits.y-rel.top;
                ExpungeSelection();
                if ((hf>=hs)||(vs>=vf)) {
                    _SimpleList sel;
                    if (selectionType&HY_TABLE_SEL_ROWS) {
                        sel.RequestSpace (vf-vs+1);
                        for (h=vs; h<=vf; h++) {
                            sel<<h;
                        }
                        SetRowSelection (sel);
                    } else if (selectionType&HY_TABLE_SEL_COLS) {
                        sel.RequestSpace (hf-hs+1);
                        for (v=hs; v<=hf; h++) {
                            sel<<v;
                        }
                        SetColumnSelection(sel);

                    } else {
                        sel.RequestSpace ((vf-vs+1)*(hf-hs+1));
                        for (h=hs; h<=hf; h++)
                            for (v=vs; v<=vf; v++) {
                                sel << v*horizontalSpaces.lLength + h;
                            }
                        SetSelection (sel,true);
                        _MarkCellsForUpdate (sel);
                    }
                }
            }
            _ResetCursorState ();
        }

        break;
    }
    }


    return _HYPlatformComponent::_ProcessOSEvent (vEvent);
}




//__________________________________________________________________

void        _HYTable::_PrintTable (_SimpleList& columns, _SimpleList& rows, _HYTable* ch)
{
    _String TBD ("This feature has not yet been implemented in the GTK+ port of HyPhy.");
    ProblemReport (TBD);
}


//EOF
