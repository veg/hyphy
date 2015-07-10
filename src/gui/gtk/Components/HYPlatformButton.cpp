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
#include "HYButton.h"
#include "HYUtils.h"
#include "HYEventTypes.h"
#include "HYWindow.h"
#include "HYGraphicPane.h"
#include <gdk/gdkkeysyms.h>

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

//__________________________________________________________________

void        hy_clicked_button_callback      (GtkButton *button, gpointer user_data)
{
    _HYButton * theb = (_HYButton*)user_data;
    if (theb->messageRecipient) {
        theb->messageRecipient->ProcessEvent (generateButtonPushEvent (theb->GetID(),0));
    }
}

//__________________________________________________________________

_HYPlatformButton::_HYPlatformButton    (void)
{
    _HYButton *     parent = (_HYButton*)this;
    buttonControl          = nil;
    buttonRect             = (GdkRectangle) {
        0,0,100,100
    };
    buttonFontDesc         = pango_font_description_new ();
}

//__________________________________________________________________

_HYPlatformButton::~_HYPlatformButton   (void)
{
    pango_font_description_free (buttonFontDesc);
}

//__________________________________________________________________

void    _HYPlatformButton::_SetBackColor (_HYColor& c)
{
    bgColor = HYColorToGDKColor(c);
}


//__________________________________________________________________

void        _HYPlatformButton::_SetFont (_HYFont& f)
{
    //_HYFont f2 = f;
    //f2.size = f.size*fontConversionFactor;
    HYFont2PangoFontDesc(f,buttonFontDesc);
    if (buttonControl) {
        _ApplyFont ();
    }
}


//__________________________________________________________________

void        _HYPlatformButton::_ApplyFont(void)
{
    gtk_widget_modify_font (gtk_bin_get_child(GTK_BIN(buttonControl)), buttonFontDesc);
}

//__________________________________________________________________
void        _HYPlatformButton::_Update (Ptr p)
{
    _Paint (p);
}

//__________________________________________________________________
void        _HYPlatformButton::_SetDimensions (_HYRect r, _HYRect rel)
{
    _HYButton* theParent = (_HYButton *) this;
    theParent->_HYPlatformComponent::_SetDimensions (r,rel);
    _SetVisibleSize (rel);
}

//__________________________________________________________________
void        _HYPlatformButton::_EnableButton (bool e)
{
    if (buttonControl) {
        gtk_widget_set_sensitive(buttonControl,e);
    }
}


//__________________________________________________________________
void        _HYPlatformButton::_SetButtonKind (unsigned char k)
{
    if (buttonControl)
        if (k==HY_BUTTON_OK) {
            GTK_WIDGET_SET_FLAGS (buttonControl,  GTK_CAN_DEFAULT|GTK_HAS_DEFAULT);
        } else {
            GTK_WIDGET_UNSET_FLAGS (buttonControl,  GTK_CAN_DEFAULT|GTK_HAS_DEFAULT);
        }
}


//__________________________________________________________________
void        _HYPlatformButton::_SetVisibleSize (_HYRect rel)
{
    _HYButton * theParent = (_HYButton*) this;

    buttonRect.x    =   rel.left;
    buttonRect.y    =   rel.top;

    _HYRect s = theParent->_SuggestDimensions();

    buttonRect.width =  s.right;
    buttonRect.height = s.bottom;

    AlignRectangle (rel, buttonRect, theParent->GetAlignFlags());

    if (buttonControl) {
        gtk_fixed_move (GTK_FIXED(theParent->parentWindow), buttonControl, buttonRect.x, buttonRect.y);
        gtk_widget_set_size_request(buttonControl,buttonRect.width,buttonRect.height);
        //SizeControl (buttonControl,buttonRect.right-buttonRect.left+1,buttonRect.bottom-buttonRect.top+1);
        //MoveControl (buttonControl,buttonRect.left,buttonRect.top);
    }
}


//__________________________________________________________________
void        _HYPlatformButton::_Paint (Ptr p)
{
    _HYButton * theParent = (_HYButton*)this;

    GdkRectangle        cRect = HYRect2GDKRect(*(_HYRect*)p);

    if (!(theParent->settings.width&HY_COMPONENT_TRANSP_BG)) {
        if (theParent->parentWindow->window) {
            GdkGC *buttonGC              = gdk_gc_new (theParent->parentWindow->window);
            gdk_gc_set_foreground(buttonGC,&bgColor);
            GdkRegion * r1 = gdk_region_rectangle(&cRect),
                        * r2 = gdk_region_rectangle (&buttonRect);

            gdk_region_subtract    (r1,r2);
            gdk_region_offset      (r1,theParent->parentWindow->allocation.x,theParent->parentWindow->allocation.y);
            gdk_gc_set_clip_region (buttonGC,r1);
            gdk_draw_rectangle(theParent->parentWindow->window,buttonGC,true,cRect.x+theParent->parentWindow->allocation.x,
                               cRect.y+theParent->parentWindow->allocation.y, cRect.width, cRect.height);
            gdk_region_destroy(r1);
            gdk_region_destroy(r2);
            g_object_unref (buttonGC);
        }
    }

    (*theParent)._HYPlatformComponent::_Paint(p);
}


//__________________________________________________________________
void        _HYPlatformButton::_PaintMe (void)
{
}


//__________________________________________________________________

_HYRect _HYButton::_SuggestDimensions (void)
{
    _HYRect res = {10,100,10,100,HY_COMPONENT_NO_SCROLL};

    if (buttonControl) {
        PangoLayout *buttonLayout            = pango_layout_new           (screenPContext);
        pango_layout_set_width (buttonLayout, -1);
        pango_layout_set_font_description(buttonLayout,buttonFontDesc);
        pango_layout_set_text(buttonLayout,buttonText.sData,buttonText.sLength);
        PangoRectangle extents,
                       logical_ext;

        pango_layout_get_pixel_extents (buttonLayout,&extents,&logical_ext);
        //pango_layout_get_extents     (buttonLayout,&logical_ext,NULL);

        res.right   = (res.left = logical_ext.width  + 20);
        res.bottom  = (res.top  = logical_ext.height + 9);
        g_object_unref              (buttonLayout);
    }

    return res;
}


//__________________________________________________________________

void        _HYButton::_Activate (void)
{
}

//__________________________________________________________________

void        _HYButton::_Deactivate (void)
{
}

//__________________________________________________________________

void        _HYPlatformButton::_SetText (void)
{
    _HYButton *parent = (_HYButton*)this;

    if (buttonControl) {
        gtk_button_set_label (GTK_BUTTON(buttonControl),parent->GetText().sData);
    } else {
        buttonControl = gtk_button_new_with_label (parent->GetText().sData);
        gtk_container_add(GTK_CONTAINER(parent->parentWindow),buttonControl);
        gtk_widget_set_app_paintable(buttonControl,true);
        gtk_widget_show (buttonControl);
        g_signal_connect (G_OBJECT (buttonControl), "clicked", G_CALLBACK (hy_clicked_button_callback), parent);
        _ApplyFont ();
    }
}


//__________________________________________________________________

bool _HYButton::_ProcessOSEvent (Ptr vEvent)
{
    _HY_GTK_UI_Message *theMessage = (_HY_GTK_UI_Message*)vEvent;

    if(buttonControl&&isEnabled) {
        switch (theMessage->theEvent->type) {
        case GDK_KEY_PRESS: {
            GdkEventKey * kpe = (GdkEventKey*)theMessage->theEvent;

            gint    keyCode = kpe->keyval;

            bool    good    = false;

            if (buttonKind == HY_BUTTON_OK) {
                good = (keyCode==GDK_Return || keyCode==GDK_KP_Enter);
            } else if (buttonKind == HY_BUTTON_CANCEL) {
                good = keyCode==GDK_Escape;
            }

            if (good) {
                gtk_button_clicked (GTK_BUTTON(buttonControl));
                return true;
            }
        }
        }
    }

    return _HYPlatformComponent::_ProcessOSEvent (vEvent);
}
