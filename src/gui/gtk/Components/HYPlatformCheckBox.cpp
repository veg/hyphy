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
#include "HYWindow.h"
#include "HYGraphicPane.h"

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif
//__________________________________________________________________

void        hy_toggled_button_callback      (GtkToggleButton *togglebutton, gpointer user_data)
{
    _HYCheckbox * thecb = (_HYCheckbox*)user_data;
    if (thecb->isEnabled && gtk_toggle_button_get_active(togglebutton) != thecb->GetState())
        if (!(thecb->isRadio&&thecb->checkState)) {
            thecb->SetState (!thecb->checkState, true);
        }
}

//__________________________________________________________________

void        _HYCheckbox::_Activate (void)
{
    if ((!activationFlag) && (isEnabled&&checkboxControl)) {
        gtk_widget_set_sensitive(checkboxControl,true);
    }

    _HYPlatformComponent::_Activate();
}

//__________________________________________________________________

void        _HYCheckbox::_Deactivate (void)
{
    if (activationFlag&&isEnabled&&checkboxControl) {
        gtk_widget_set_sensitive(checkboxControl,false);
    }

    _HYPlatformComponent::_Deactivate();
}

//__________________________________________________________________

bool _HYCheckbox::_ProcessOSEvent (Ptr vEvent)
{
    return _HYPlatformComponent::_ProcessOSEvent (vEvent);
}

//__________________________________________________________________

_HYRect     _HYCheckbox::_SuggestDimensions (void)
{
    _HYRect res  = _HYLabel::_SuggestDimensions();
    res.left    += checkSpacing*2 + _GetBoxDimension();
    res.right    = res.left;
    return         res;
}

//__________________________________________________________________
void        _HYPlatformCheckbox::_SetVisibleSize (_HYRect rel)
{
    _HYCheckbox * theParent = (_HYCheckbox*) this;

    theParent->labelRect.x= rel.left;
    theParent->labelRect.y = rel.top;

    _HYRect s = theParent->_SuggestDimensions();

    theParent->labelRect.width =  s.right;
    theParent->labelRect.height = s.bottom;

    AlignRectangle (rel, theParent->labelRect, theParent->GetAlignFlags());

    checkboxRect.left   = theParent->labelRect.x + theParent->checkSpacing;
    checkboxRect.bottom = theParent->labelRect.y + theParent->labelRect.height;

    long             boxD = _GetBoxDimension();

    checkboxRect.right = checkboxRect.left  + boxD;
    checkboxRect.top   = checkboxRect.bottom >= boxD?checkboxRect.bottom - boxD:0;

    theParent->labelRect.x += boxD+2*theParent->checkSpacing;

    if (!checkboxControl) {
        if (isRadio) {
            checkboxControl = gtk_radio_button_new(nil);
            dummyControl    = gtk_radio_button_new (gtk_radio_button_get_group(GTK_RADIO_BUTTON(checkboxControl)));
        } else {
            checkboxControl = gtk_check_button_new();
        }

        g_signal_connect (G_OBJECT (checkboxControl), "toggled", G_CALLBACK (hy_toggled_button_callback), theParent);

        if (theParent->GetState() != gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkboxControl))) {
            theParent->_SetState (theParent->GetState());
        }

        gtk_container_add(GTK_CONTAINER(theParent->parentWindow),checkboxControl);
        gtk_widget_set_app_paintable(checkboxControl,true);
        gtk_widget_show (checkboxControl);
        g_signal_connect (G_OBJECT (checkboxControl), "event", G_CALLBACK (hyphy_event_component_callback), theParent);
    }

    gtk_fixed_move (GTK_FIXED(theParent->parentWindow), checkboxControl, checkboxRect.left,checkboxRect.top);
    boxD = MIN(checkboxControl->allocation.x,checkboxControl->allocation.y);
    long bh = MIN(checkboxRect.bottom-checkboxRect.top+1,boxD),
         bw = MIN(checkboxRect.right-checkboxRect.left+1,boxD);

    if (bh<boxD||bw<boxD) {
        gtk_widget_set_size_request(checkboxControl, bw, bh);
    }
}

//__________________________________________________________________

void        _HYPlatformCheckbox::_Enable (bool e)
{
    if (checkboxControl) {
        gtk_widget_set_sensitive(checkboxControl,e);
    }
}

//__________________________________________________________________

long        _HYPlatformCheckbox::_GetBoxDimension (void)
{
    if (checkboxControl) {
        GtkRequisition sizeReq;
        gtk_widget_size_request (checkboxControl, &sizeReq);
        return MAX(sizeReq.width, sizeReq.height);
    }
    return HY_DEFAULT_CHECK_SPACE;
}

//__________________________________________________________________
void        _HYPlatformCheckbox::_Update (Ptr p)
{
    _Paint (p);
}

//__________________________________________________________________
void        _HYPlatformCheckbox::_Paint (Ptr p)
{
    _HYCheckbox *theParent = (_HYCheckbox*) this;
    theParent->_HYPlatformLabel::_Paint(p);
}

//__________________________________________________________________
void        _HYPlatformCheckbox::_PaintMe (void)
{
    _HYCheckbox * theParent = (_HYCheckbox*)this;
    /*bool      drawBk = !(theParent->settings.width&HY_COMPONENT_TRANSP_BG);
    if (theParent->labelGC && drawBk)
    {

        GdkColor saveColor = HYColorToGDKColor(theParent->GetForeColor()),
                 newColor  = HYColorToGDKColor(theParent->GetBackColor());

        gdk_gc_set_foreground(theParent->labelGC, &newColor);
        gdk_draw_rectangle(theParent->parentWindow->window,theParent->labelGC,true,checkboxRect.left,checkboxRect.top,_GetBoxDimension(),_GetBoxDimension());
        gdk_gc_set_foreground(theParent->labelGC, &saveColor);
    }*/
    GdkRectangle    cbr = HYRect2GDKRect(checkboxRect);
    gtk_widget_draw (checkboxControl,&cbr);
}

//__________________________________________________________________
_HYPlatformCheckbox::_HYPlatformCheckbox (bool isR)
{
    checkboxControl = nil;
    dummyControl    = nil;
    isRadio = isR;
    checkboxRect = (_HYRect) {
        0,0,100,100,0
    };
}

//__________________________________________________________________
_HYPlatformCheckbox::~_HYPlatformCheckbox (void)
{
}

//__________________________________________________________________
void    _HYPlatformCheckbox::_SetState (bool v)
{
    if (checkboxControl)
        if (isRadio) {
            gtk_toggle_button_set_active (v?GTK_TOGGLE_BUTTON(checkboxControl):GTK_TOGGLE_BUTTON(dummyControl),true);
        } else {
            gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON(checkboxControl),v);
        }
}

//EOF
