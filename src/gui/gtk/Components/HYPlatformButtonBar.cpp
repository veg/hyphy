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

#include "HYPlatformGraphicPane.h"

#include "errorfns.h"
#include "HYButtonBar.h"
#include "HYUtils.h"
#include "HYEventTypes.h"

#include "HYWindow.h"
#include "HYGraphicPane.h"

#include "string.h"

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif
//__________________________________________________________________

GdkColor    buttonBorder1 = HYColorToGDKColor((_HYColor)
{
    0,0,0
}),
                            buttonBorder2 = HYColorToGDKColor((_HYColor)
{
    0x3f,0x3f,0x3f
});

_HYFont     _hyttDisplayFont = {_HY_SANS_FONT,9,HY_FONT_PLAIN};

//__________________________________________________________________

gboolean TooltipPopupTimer (Ptr* userData)
{
    // TBI

    _HYButtonBar * myBB = (_HYButtonBar*)userData;

    gint            x,y;
    GdkModifierType gmt;
    GdkWindow       *pw = gtk_widget_get_parent_window(myBB->parentWindow);
    if (pw) {
        gdk_window_get_pointer (pw, &x, &y, &gmt);

        if ( x==myBB->saveMousePosH && y==myBB->saveMousePosV && myBB->toolTipBounds.x == 0) {
            myBB->_DisplayToolTip();
        }

        myBB->saveMousePosH = x;
        myBB->saveMousePosV = y;
    }
    return true;
}

//__________________________________________________________________

_HYPlatformButtonBar::_HYPlatformButtonBar(void)
{
    backFill           = HYColorToGDKColor ((_HYColor) {
        255,255,255
    });
    pushed             = -1;
    saveMousePosH      = -1;
    saveMousePosV      = -1;
    lastMouseDown      = -1;
    toolTipBounds.x    = 0;
    theTimer           = 0;
    bbGC               = nil;
}

//__________________________________________________________________

_HYPlatformButtonBar::~_HYPlatformButtonBar(void)
{
    if (theTimer) {
        g_source_remove (theTimer);
        theTimer = 0;
    }
    if (bbGC) {
        g_object_unref (bbGC);
    }
}

//__________________________________________________________________

void    _HYPlatformButtonBar::_DisposeButtons(void)
{
}

//__________________________________________________________________

void    _HYPlatformButtonBar::_DisposeButton(long k)
{
}


//__________________________________________________________________

void        _HYPlatformButtonBar::_SetBackColor (_HYColor& c)
{
    backFill = HYColorToGDKColor(c);
}

//__________________________________________________________________
void        _HYPlatformButtonBar::_SetVisibleSize (_HYRect rel)
{
    _HYButtonBar* theParent = (_HYButtonBar*) this;

    buttonRect.x = rel.left;
    buttonRect.y = rel.top;

    _HYRect s = theParent->_SuggestDimensions();

    buttonRect.width  = s.right;
    buttonRect.height = s.bottom;

    AlignRectangle (rel, buttonRect, theParent->GetAlignFlags());
}


//__________________________________________________________________

void        _HYButtonBar::_Activate (void)
{
    if (!activationFlag)
        for (long k=0; k<enabledButtons.lLength; k++) {
            _MarkButtonForUpdate (enabledButtons.lData[k]);
        }

    if (!theTimer) {
        theTimer = g_timeout_add  (1000,(GSourceFunc)TooltipPopupTimer,(gpointer)this);
    }

    _HYPlatformComponent::_Activate();
}

//__________________________________________________________________
void        _HYPlatformButtonBar::_SmiteTooltip (void)
{
    _HYButtonBar * p = (_HYButtonBar*)this;
    if (toolTipBounds.x) {
        if (p->parentWindow->window) {
            toolTipBounds.width++;
            toolTipBounds.height++;
            gdk_window_invalidate_rect (p->parentWindow->window, &toolTipBounds, false);
        }
        toolTipBounds.x = 0;
    }
}

//__________________________________________________________________

void        _HYButtonBar::_Deactivate (void)
{
    if (activationFlag) {
        for (long k=0; k<enabledButtons.lLength; k++) {
            _MarkButtonForUpdate (enabledButtons.lData[k]);
        }
        _SmiteTooltip();
    }

    if (theTimer) {
        g_source_remove (theTimer);
        theTimer = 0;
    }

    _HYPlatformComponent::_Deactivate();
}


//__________________________________________________________________

void        _HYButtonBar::_ComponentMouseExit (void)
{
    _UnpushButton();
    _SmiteTooltip();

    /*if (theTimer)
    {
        g_source_remove (theTimer);
        theTimer = 0;
    }*/
}

//__________________________________________________________________
void        _HYPlatformButtonBar::_MarkButtonForUpdate (int i)
{
    _HYButtonBar* theParent = (_HYButtonBar*)this;

    if ( i>=0 && i<theParent->ButtonCount()) {
        int hR = i%theParent->BarWidth(),
            vR = i/theParent->BarWidth(),
            step = 2*HY_BUTTONBAR_BORDER+theParent->GetButtonDim();

        GdkRectangle          invRect;
        invRect.x           = buttonRect.x + hR*step + theParent->parentWindow->allocation.x;
        invRect.y           = buttonRect.y + vR*step + theParent->parentWindow->allocation.y;
        invRect.width       = invRect.height = step;

        if (theParent->parentWindow->window) {
            gdk_window_invalidate_rect (theParent->parentWindow->window, &invRect, false);
        }
    }
}


//__________________________________________________________________
void        _HYPlatformButtonBar::_Paint (Ptr p)
{
    _HYButtonBar * theParent = (_HYButtonBar*)this;

    if (theParent->parentWindow->window) {
        GdkRectangle    cRect    = HYRect2GDKRect(*(_HYRect*)p),
                        iRect;

        if (!bbGC) {
            bbGC                 = gdk_gc_new (theParent->parentWindow->window);
            gdk_gc_set_line_attributes (bbGC, 1, GDK_LINE_SOLID, GDK_CAP_NOT_LAST, GDK_JOIN_MITER);
        }

        cRect.x += theParent->parentWindow->allocation.x;
        cRect.y += theParent->parentWindow->allocation.y;

        if (!(theParent->settings.width&HY_COMPONENT_TRANSP_BG)) {
            gdk_gc_set_foreground(bbGC,&backFill);
            gdk_draw_rectangle(theParent->parentWindow->window,bbGC,true,cRect.x, cRect.y, cRect.width, cRect.height);
        }


        gdk_gc_set_clip_rectangle(bbGC,&cRect);

        cRect.x  = buttonRect.x;
        cRect.y  = buttonRect.y;

        int   step   = theParent->GetButtonDim()+2*HY_BUTTONBAR_BORDER;

        cRect.width  = cRect.height = step;

        gdk_gc_set_foreground(bbGC,&buttonBorder1);
        for (long i=0; i<theParent->ButtonCount(); i++) {
            if (i&& i%theParent->BarWidth()==0) {
                cRect.x  =  buttonRect.x;
                cRect.y  += step;
            }

            iRect    = cRect;
            iRect.x += theParent->parentWindow->allocation.x + HY_BUTTONBAR_BORDER;
            iRect.y += theParent->parentWindow->allocation.y + HY_BUTTONBAR_BORDER;

            iRect.width  -= 2*HY_BUTTONBAR_BORDER;
            iRect.height -= 2*HY_BUTTONBAR_BORDER;

            GdkPixbuf * buttonToPlot = nil;

            if (theParent->activationFlag) {
                if (i==pushed) {
                    buttonToPlot = gdk_pixbuf_composite_color_simple ((GdkPixbuf*)theParent->buttons.lData[i], iRect.width, iRect.height,GDK_INTERP_BILINEAR,128,1024,0x00000000,0x00000000);
                } else {
                    if (theParent->enabledButtons.Find(i)>=0) {
                        buttonToPlot = gdk_pixbuf_composite_color_simple ((GdkPixbuf*)theParent->buttons.lData[i], iRect.width, iRect.height,GDK_INTERP_BILINEAR,255,1024,0,0);
                    } else {
                        buttonToPlot = gdk_pixbuf_composite_color_simple ((GdkPixbuf*)theParent->buttons.lData[i], iRect.width, iRect.height,GDK_INTERP_BILINEAR,128,1,0x00ffffff,0x00888888);
                    }
                }
            } else {
                buttonToPlot = gdk_pixbuf_composite_color_simple ((GdkPixbuf*)theParent->buttons.lData[i], iRect.width, iRect.height,GDK_INTERP_BILINEAR,128,1,0x00ffffff,0x00888888);
            }
            gdk_draw_pixbuf (theParent->parentWindow->window,bbGC, buttonToPlot, 0,0, iRect.x, iRect.y, iRect.width, iRect.height, GDK_RGB_DITHER_NORMAL, 0, 0);
            g_object_unref (buttonToPlot);
            iRect.x --;
            iRect.y --;
            iRect.width += 2;
            iRect.height += 2;

            gdk_draw_rectangle(theParent->parentWindow->window,bbGC,false,iRect.x, iRect.y, iRect.width, iRect.height);
            cRect.x +=step;
        }
        gdk_gc_set_clip_rectangle(bbGC,nil);

    }

    (*theParent)._HYPlatformComponent::_Paint(p);
}


//__________________________________________________________________
_HYRect _HYPlatformButtonBar::_GetButtonRect (bool conv)
{
    _HYRect         res = GdkRect2HYRect (buttonRect);
    _HYButtonBar* theParent = (_HYButtonBar*) this;

    if (conv) {
        res.left   += theParent->parentWindow->parent->allocation.x;
        res.right  += theParent->parentWindow->parent->allocation.x;
        res.top    += theParent->parentWindow->parent->allocation.y;
        res.bottom += theParent->parentWindow->parent->allocation.y;
    }
    return res;
}


//__________________________________________________________________
void        _HYPlatformButtonBar::_Update (Ptr p)
{
    _Paint (p);
}

//__________________________________________________________________
void        _HYPlatformButtonBar::_UnpushButton (void)
{
    if (pushed>=0) {
        _MarkButtonForUpdate(pushed);
        pushed = -1;
        lastMouseDown = -1;
    }
}

//__________________________________________________________________
void        _HYPlatformButtonBar::_SetDimensions (_HYRect r, _HYRect rel)
{
    _HYButtonBar* theParent = (_HYButtonBar*) this;
    theParent->_HYPlatformComponent::_SetDimensions (r,rel);
    _SetVisibleSize (rel);
}

//__________________________________________________________________
int         _HYPlatformButtonBar::_FindClickedButton (int v, int h)
{
    _HYButtonBar * parent = (_HYButtonBar*)this;

    v -= parent->parentWindow->allocation.y;
    h -= parent->parentWindow->allocation.x;

    if ( h>=buttonRect.x && h<buttonRect.x + buttonRect.width && v>=buttonRect.y && v<buttonRect.y + buttonRect.height ) {
        v -= buttonRect.y;
        h -= buttonRect.x;

        int step = 2*HY_BUTTONBAR_BORDER+parent->buttonDim,
            hR = h/step,
            vR = v/step;

        v  -= vR*step;
        h  -= hR*step;

        if (v>HY_BUTTONBAR_BORDER && v<step-HY_BUTTONBAR_BORDER && h>HY_BUTTONBAR_BORDER && h<step-HY_BUTTONBAR_BORDER) {
            return hR+vR*parent->barW;
        }
    }
    return -1;

}

//__________________________________________________________________

bool _HYButtonBar::_ProcessOSEvent (Ptr vEvent)
{
    if (buttons.lLength) {
        _HY_GTK_UI_Message *theMessage = (_HY_GTK_UI_Message*)vEvent;

        gdouble   xc,
                  yc;

        if (gdk_event_get_coords (theMessage->theEvent,&xc,&yc)) {
            switch (theMessage->theEvent->type) {
            case GDK_BUTTON_PRESS: {
                if (((GdkEventButton*)theMessage->theEvent)->button != 1) {
                    return false;
                }
                _SmiteTooltip();
                int h = _FindClickedButton (yc,xc);
                if ( h>=0 && enabledButtons.Find(h) >=0 ) {
                    if (pullDownButtons.Find (h)>=0) {
                        lastMouseDown = -1;
                        SendButtonPush(h);
                    } else {
                        pushed  = lastMouseDown = h;
                        _MarkButtonForUpdate (h);
                    }
                }
                return true;
            }

            case GDK_BUTTON_RELEASE: {
                if (lastMouseDown>=0 || pushed >= 0) {
                    int h = _FindClickedButton (yc,xc);
                    if (h==pushed && pushed >= 0) {
                        SendButtonPush (h);
                    }
                    _MarkButtonForUpdate (pushed);
                    _MarkButtonForUpdate (lastMouseDown);
                    pushed        = -1;
                    lastMouseDown = -1;
                }
                return true;
            }

            case GDK_MOTION_NOTIFY : {
                _SmiteTooltip();
                if (lastMouseDown>=0) {
                    int h = _FindClickedButton (yc,xc);
                    if (h!=lastMouseDown) {
                        if (pushed>=0) {
                            pushed = -1;
                            _MarkButtonForUpdate (lastMouseDown);
                        }
                    } else {
                        if (pushed==-1) {
                            pushed = lastMouseDown;
                            _MarkButtonForUpdate (lastMouseDown);
                        }
                    }
                }
                return true;
            }
            }
        }
    }
    return _HYPlatformComponent::_ProcessOSEvent (vEvent);
}


//__________________________________________________________________
void        _HYButtonBar::_DisplayToolTip      (void)
{
    // TBI

    int h = _FindClickedButton (saveMousePosV,saveMousePosH);

    if ( h>=0 && h<toolTips.lLength) {
        _String* toolTip = (_String*)toolTips(h);
        if (toolTip->sLength) {
            GdkColor      toolTipColor = HYColorToGDKColor((_HYColor) {
                0xFF,0xCC,0x66
            });

            int        bL,
                       bT;

            GetButtonLoc (h,bL,bT,false);

            PangoLayout          * ttText = pango_layout_new (screenPContext);
            PangoFontDescription * ttFont = pango_font_description_new();
            HYFont2PangoFontDesc(_hyttDisplayFont,ttFont);

            pango_layout_set_font_description (ttText, ttFont);
            pango_layout_set_width(ttText, -1);
            pango_layout_set_text (ttText, toolTip->sData, toolTip->sLength);

            PangoRectangle extents;
            pango_layout_get_pixel_extents (ttText,&extents,nil);

            toolTipBounds.y      = bT - 16;
            toolTipBounds.height = 16;

            if (toolTipBounds.y<0) {
                toolTipBounds.y = bT+buttonDim+1;
            }

            h = extents.width + 4;

            toolTipBounds.x = bL+(buttonDim-h-1)/2;

            if (toolTipBounds.x<=0) {
                toolTipBounds.x = 1;
            }

            toolTipBounds.width = h+2;

            h = toolTipBounds.x+toolTipBounds.width-parentWindow->allocation.width;

            if (h>0) {
                if (h>=toolTipBounds.x) {
                    h = toolTipBounds.x-1;
                }
                toolTipBounds.x -= h;
            }

            toolTipBounds.x += parentWindow->allocation.x;
            toolTipBounds.y += parentWindow->allocation.y;

            gdk_gc_set_foreground(bbGC,&toolTipColor);
            gdk_draw_rectangle(parentWindow->window,bbGC,true,toolTipBounds.x, toolTipBounds.y, toolTipBounds.width, toolTipBounds.height);

            toolTipColor = HYColorToGDKColor((_HYColor) {
                0,0,0
            });
            gdk_gc_set_foreground(bbGC,&toolTipColor);
            gdk_draw_rectangle(parentWindow->window,bbGC,false,toolTipBounds.x, toolTipBounds.y, toolTipBounds.width, toolTipBounds.height);

            gdk_gc_set_clip_rectangle  (bbGC, &toolTipBounds);
            gdk_draw_layout (parentWindow->window, bbGC, toolTipBounds.x+3, toolTipBounds.y+3, ttText);
            gdk_gc_set_clip_rectangle  (bbGC, nil);

            pango_font_description_free (ttFont);
            g_object_unref (ttText);

        }
    }
}
//EOF
