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

#include "HYLabel.h"
#include "HYUtils.h"
#include "HYGraphicPane.h"
#include "HYWindow.h"

//__________________________________________________________________

_HYPlatformLabel::_HYPlatformLabel(void)
{
    labelLayout          = pango_layout_new           (screenPContext);
    labelFontDesc        = pango_font_description_new ();
    if (((_HYLabel*)this)->parentWindow->window) {
        labelGC              = gdk_gc_new                 (((_HYLabel*)this)->parentWindow->window);
    } else {
        labelGC              = nil;
    }

    pango_layout_set_width (labelLayout, -1);
    pango_layout_set_alignment(labelLayout,PANGO_ALIGN_CENTER);
}

//__________________________________________________________________

_HYPlatformLabel::~_HYPlatformLabel(void)
{
    g_object_unref              (labelLayout);
    pango_font_description_free (labelFontDesc);
    if (labelGC) {
        g_object_unref (labelGC);
    }
}

//__________________________________________________________________

void        _HYPlatformLabel::_SetBackColor (_HYColor& c)
{
    if (labelGC) {
        GdkColor    newBG = HYColorToGDKColor(c);
        gdk_gc_set_background(labelGC,&newBG);
    }
}

//__________________________________________________________________

void        _HYPlatformLabel::_SetForeColor (_HYColor& c)
{
    if (labelGC) {
        GdkColor    newFG = HYColorToGDKColor(c);
        gdk_gc_set_foreground(labelGC,&newFG);
    }
}

//__________________________________________________________________

void        _HYPlatformLabel::_SetFont (_HYFont& f)
{
    HYFont2PangoFontDesc(f,labelFontDesc);
    pango_layout_set_font_description (labelLayout, labelFontDesc ); // ref ?
}

//__________________________________________________________________
void        _HYPlatformLabel::_Update (Ptr p)
{
    _Paint (p);
}

//__________________________________________________________________
void        _HYPlatformLabel::_SetDimensions (_HYRect r, _HYRect rel)
{
    _HYLabel* theParent = (_HYLabel*) this;
    theParent->_HYPlatformComponent::_SetDimensions (r,rel);
    _SetVisibleSize (rel);
}


//__________________________________________________________________
void        _HYPlatformLabel::_SetVisibleSize (_HYRect rel)
{
    _HYLabel* theParent = (_HYLabel*) this;
    _HYRect s = theParent->_SuggestDimensions();

    labelRect.x         =   rel.left;
    labelRect.y         =   rel.top;
    labelRect.width     =   s.right<rel.right-rel.left?s.right:rel.right-rel.left+1;
    labelRect.height    =   s.bottom<rel.bottom-rel.top?s.bottom:rel.bottom-rel.top+1;

    AlignRectangle (rel, labelRect, theParent->GetAlignFlags());

    //printf ("Label layout %d %d %d %d\n", labelRect.x, labelRect.y, labelRect.width, labelRect.height);
    //pango_layout_set_width(labelLayout,labelRect.width);
}

//__________________________________________________________________
_HYRect _HYLabel::_SuggestDimensions (void)
{
    _HYRect res = {10,100,10,100,HY_COMPONENT_NO_SCROLL};
    PangoRectangle extents;
    pango_layout_get_pixel_extents (labelLayout,nil,&extents);
    res.top = res.bottom = extents.height+4;
    res.left = res.right = extents.width + 5;
    return res;
}

//__________________________________________________________________
void    _HYPlatformLabel::_SetText (void)
{
    _HYLabel *      theParent = (_HYLabel*)this;
    pango_layout_set_text(labelLayout,theParent->GetText().sData,theParent->GetText().sLength);
}


//__________________________________________________________________
void        _HYPlatformLabel::_Paint (Ptr p)
{
    _HYLabel *      theParent = (_HYLabel*)this;
    if (!labelGC) {
        _SetBackColor       (theParent->GetBackColor());
        _SetForeColor       (theParent->GetForeColor());
        labelGC              = gdk_gc_new                 (((_HYLabel*)this)->parentWindow->window);
    }


    GdkRectangle    cRect   = HYRect2GDKRect(*(_HYRect*)p);
    GdkColor        forC    = HYColorToGDKColor(theParent->GetForeColor()),
                    aColor;
    cRect.x += theParent->parentWindow->allocation.x;
    cRect.y += theParent->parentWindow->allocation.y;

    if (!(theParent->settings.width&HY_COMPONENT_TRANSP_BG)) {
        aColor = HYColorToGDKColor(theParent->GetBackColor());
        gdk_gc_set_foreground(labelGC,&aColor);
        gdk_draw_rectangle(theParent->parentWindow->window,labelGC,true,cRect.x, cRect.y, cRect.width, cRect.height);
        gdk_gc_set_foreground(labelGC,&forC);
    }

    if (theParent->settings.width&HY_COMPONENT_WELL) {
        cRect.x++;
        cRect.y++;
        cRect.width  -= 2;
        cRect.height -= 2;
    }

    gdk_gc_set_clip_rectangle(labelGC,&cRect);

    if (theParent->HasShadow()) {
        _HYColor  tempColor = theParent->GetForeColor();

        tempColor.R /= 4;
        tempColor.G /= 4;
        tempColor.B /= 4;

        GdkColor  tc = HYColorToGDKColor(tempColor);

        gdk_gc_set_foreground(labelGC,&tc);
        gdk_draw_layout(theParent->parentWindow->window,labelGC,theParent->parentWindow->allocation.x+labelRect.x+1,
                        theParent->parentWindow->allocation.y+labelRect.y-1,
                        labelLayout);
        tempColor = theParent->GetForeColor();
        tempColor.R /= 2;
        tempColor.G /= 2;
        tempColor.B /= 2;

        tc = HYColorToGDKColor(tempColor);
        gdk_gc_set_foreground(labelGC,&tc);

        gdk_draw_layout(theParent->parentWindow->window,labelGC,theParent->parentWindow->allocation.x+labelRect.x-1,
                        theParent->parentWindow->allocation.y+labelRect.y+1,
                        labelLayout);
        gdk_gc_set_foreground(labelGC,&forC);
    }

    gdk_draw_layout(theParent->parentWindow->window,labelGC,theParent->parentWindow->allocation.x+labelRect.x,
                    theParent->parentWindow->allocation.y+labelRect.y,
                    labelLayout);
    gdk_gc_set_clip_rectangle(labelGC,nil);

    (*theParent)._HYPlatformComponent::_Paint(p);
}

// EOF
