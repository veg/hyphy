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
    A general composite window component object, GTK+ glue

    Sergei L. Kosakovsky Pond, October 2004.
*/

#include "HYComponent.h"
#include "HYPlatformComponent.h"
#include "HYWindow.h"
#include "HYEventTypes.h"
#include "HYTableWindow.h"
#include "HYCanvas.h"
#include "HYUtils.h"

//__________________________________________________________________________________


//extern    RGBColor    menuLine1,
//                      menuLine2; -- have to fix thouse

// GTK Callbacks

bool        forceUpdateForScrolling = false;

//__________________________________________________________________

gboolean hyphy_event_component_callback(GtkWidget *widget, GdkEvent* theEvent, gpointer   data)
{
    _HYComponent* pC = (_HYComponent*) data;
    if (pC->messageRecipient) {
        _HY_GTK_UI_Message theMessage = {FALSE,theEvent};
        ((_HYTWindow*)pC->messageRecipient)->_ProcessOSEvent((Ptr)&theMessage);
        return theMessage.processingResult;
    }
    return FALSE;
}

//__________________________________________________________________________________

void            AlignRectangle (_HYRect& rel , GdkRectangle& target , unsigned char alFlags)
{
    if (alFlags&HY_ALIGN_RIGHT) {
        target.x  = rel.right + 1 - target.width;
    } else if (!(alFlags&HY_ALIGN_LEFT)) {
        target.x = rel.left + (rel.right-rel.left + 1 - target.width)/2;
    }

    if (alFlags&HY_ALIGN_BOTTOM) {
        target.y = rel.bottom - target.height + 1;
    } else if (!(alFlags&HY_ALIGN_TOP)) {
        target.y = rel.top + (rel.bottom-rel.top + 1 - target.height)/2;
    }
}

//__________________________________________________________________

void h_scroll_bar_callback_component (GtkRange *widget, gpointer data)
{
    GtkAdjustment* ta = gtk_range_get_adjustment (widget);
    double newV  = ta->value*(ta->upper/(ta->upper-ta->page_size));
    _HYComponent * parent_obj = (_HYComponent*)data;
    long diff = newV-parent_obj->lastHScroll;
    if (diff) {
        parent_obj->lastHScroll = newV;
        parent_obj->ProcessEvent (generateScrollEvent(diff,0));
    }
}

//__________________________________________________________________

void v_scroll_bar_callback_component (GtkRange *widget, gpointer data)
{
    GtkAdjustment* ta = gtk_range_get_adjustment (widget);
    double newV  = ta->value*(ta->upper/(ta->upper-ta->page_size));
    _HYComponent * parent_obj = (_HYComponent*)data;
    long diff = (newV-parent_obj->lastVScroll);
    if (diff) {
        parent_obj->lastVScroll = newV;
        parent_obj->ProcessEvent (generateScrollEvent(0,diff));
    }
}

//__________________________________________________________________

_HYPlatformComponent::_HYPlatformComponent(void)
{
    vScroll = hScroll = parentWindow = nil;
    activationFlag = false;
}

//__________________________________________________________________

_HYPlatformComponent::_HYPlatformComponent(_HYRect s,Ptr w)
{
    activationFlag = false;
    bool        memError = false;
    vScroll     = hScroll = nil;
    parentWindow = (GtkWidget*)w;
    if (s.width&HY_COMPONENT_V_SCROLL) {
        GtkObject
        *v_adj = gtk_adjustment_new(0,0,MAX_CONTROL_VALUE,MAX_CONTROL_VALUE/100.0,MAX_CONTROL_VALUE/10.0,MAX_CONTROL_VALUE/5.0);

        vScroll = gtk_vscrollbar_new  ((GtkAdjustment*)v_adj),
        gtk_container_add (GTK_CONTAINER (parentWindow), vScroll);
        gtk_widget_show(vScroll);

        g_signal_connect (G_OBJECT (vScroll), "value-changed",G_CALLBACK (v_scroll_bar_callback_component), (gpointer)((_HYComponent*)this));
    }
    if (s.width&HY_COMPONENT_H_SCROLL) {
        GtkObject
        *h_adj = gtk_adjustment_new(0,0,MAX_CONTROL_VALUE,MAX_CONTROL_VALUE/100.0,MAX_CONTROL_VALUE/10.0,MAX_CONTROL_VALUE/5.0);

        hScroll = gtk_hscrollbar_new  ((GtkAdjustment*)h_adj),
        gtk_container_add (GTK_CONTAINER (parentWindow), hScroll);
        gtk_widget_show(hScroll);

        g_signal_connect (G_OBJECT (hScroll), "value-changed",G_CALLBACK (h_scroll_bar_callback_component), (gpointer)((_HYComponent*)this));
    }
    lastHScroll = 0;
    lastVScroll = 0;
}


//__________________________________________________________________

void _HYPlatformComponent::_CleanUp(void)
{
    if (hScroll) {
        //gtk_widget_destroy(hScroll);
        hScroll = nil;
    }
    if (vScroll) {
        //gtk_widget_destroy(vScroll);
        vScroll = nil;
    }
}

//__________________________________________________________________
long        _HYPlatformComponent::_GetScrollerPos (GtkWidget *theSB)
{
    if (theSB) {
        GtkAdjustment* old = gtk_range_get_adjustment (GTK_RANGE(theSB));
        return old->upper-old->page_size?old->value*(old->upper/(old->upper-old->page_size)):0.0;
    } else {
        return 0;
    }
}


//__________________________________________________________________
long        _HYPlatformComponent::_GetHScrollerPos (void)
{
    return _GetScrollerPos(hScroll);
}
//__________________________________________________________________
long        _HYPlatformComponent::_GetVScrollerPos (void)
{
    return _GetScrollerPos(vScroll);
}

//__________________________________________________________________
void        _HYPlatformComponent::_SetScrollerPos (GtkWidget *theSB, long nv)
{
    if (theSB) {
        if (nv<0) {
            nv = 0;
        } else if (nv>MAX_CONTROL_VALUE) {
            nv = MAX_CONTROL_VALUE;
        }

        GtkAdjustment* old = gtk_range_get_adjustment (GTK_RANGE(theSB));
        //printf ("Setting adjustment value at %d (%d)\n",nv,lastVScroll);
        //gtk_range_set_value (GTK_RANGE(theSB),old->upper-old->page_size?nv*((old->upper-old->page_size)/old->upper):0.0);
        gtk_adjustment_set_value (old, old->upper-old->page_size?nv*((old->upper-old->page_size)/old->upper):0.0);
    } 
}

//__________________________________________________________________
void        _HYPlatformComponent::_SetHScrollerPos (long nv)
{
    lastHScroll = nv;
    _SetScrollerPos (hScroll, nv);
}

//__________________________________________________________________
void        _HYPlatformComponent::_SetVScrollerPos (long nv)
{
    lastVScroll = nv;
    _SetScrollerPos (vScroll, nv);
}


//__________________________________________________________________
void _HYPlatformComponent::Duplicate (BaseRef s)
{
    // TBI - depending on what calls use it - for now should barf

    _String errMsg ("_HYPlatformComponent::Duplicate is not yet implemented");
    FlagError (errMsg);

    //_HYComponent* theS = (_HYComponent*)s;
    //_CleanUp();
    //hScroll = theS->hScroll;
    //vScroll = theS->vScroll;
}


//__________________________________________________________________
void _HYPlatformComponent::_SetDimensions (_HYRect,_HYRect r)
{
    _SetVisibleSize (r);
}

//__________________________________________________________________
void _HYPlatformComponent::_MarkForUpdate (void)
{
    GdkWindow * pW = gtk_widget_get_parent_window (parentWindow);
    if (pW) {
        GdkRectangle inv = {parentWindow->allocation.x+rel.left, parentWindow->allocation.y+rel.top,
                            rel.right-rel.left+1, rel.bottom-rel.top+1
                           }; // ? offset by 1
        gdk_window_invalidate_rect (pW, &inv, false);
        if (forceUpdateForScrolling) {
            gdk_window_process_updates (pW, false);
        }
    }
}
//__________________________________________________________________
void _HYPlatformComponent::_MarkContentsForUpdate (void)
{
    GdkWindow * pW = gtk_widget_get_parent_window (parentWindow);
    if (pW) {
        GdkRectangle inv = {parentWindow->allocation.x+rel.left, parentWindow->allocation.y+rel.top,
                            rel.right-rel.left+1, rel.bottom-rel.top+1
                           }; // ? offset by 1
        if (hScroll) {
            inv.height -= hScroll->allocation.height;
        }
        if (vScroll) {
            inv.width  -= vScroll->allocation.height;
        }
        gdk_window_invalidate_rect (pW, &inv, false);
    }
}

//__________________________________________________________________

void        _HYPlatformComponent::_SetVisibleSize (_HYRect r)
{
    _HYComponent * theParent = (_HYComponent*)this;

    long           t,
                   v,
                   tt;

    _Parameter  newSize;

    GtkFixed *  wPane = GTK_FIXED(parentWindow);

    if (hScroll&&!vScroll) { // only horizontal scroll bar
        tt = hScroll->allocation.height;
        if (tt <= 1) {
            tt = HY_SCROLLER_WIDTH;
        }

        gtk_fixed_move (wPane, hScroll, r.left, r.bottom - tt);
        gtk_widget_set_size_request(hScroll,r.right-r.left+1,tt);
        t = theParent->GetMaxW();
        v = r.right-r.left+1;

        if (t>v) {
            AdjustScroller (hScroll, v, t);
        } else {
            AdjustScroller (hScroll, t, t);
        }

        gtk_widget_set_sensitive(hScroll,t>v);
    }
    if (vScroll&&!hScroll) { // only vertical scroll bar
        tt = vScroll->allocation.width;
        if (tt <= 1) {
            tt = HY_SCROLLER_WIDTH;
        }

        gtk_fixed_move (wPane, vScroll, r.right-tt, r.top);
        gtk_widget_set_size_request(vScroll,tt,r.bottom-r.top+1);
        t = theParent->GetMaxH();
        v = r.bottom-r.top+1;
        //printf ("Vertical Scroll size %d %d\n", t, v);
        if (t>v) {
            AdjustScroller (vScroll, v, t);
        } else {
            AdjustScroller (vScroll, t, t);
        }
        gtk_widget_set_sensitive(vScroll,t>v);
    }

    if (vScroll&&hScroll) {
        tt = hScroll->allocation.height;
        if (tt <= 1) {
            tt = HY_SCROLLER_WIDTH;
        }

        long ttv = vScroll->allocation.width;
        if (ttv <= 1) {
            ttv = HY_SCROLLER_WIDTH;
        }

        gtk_fixed_move                  (wPane, hScroll, r.left, r.bottom - tt + 1);
        gtk_widget_set_size_request     (hScroll,r.right-r.left-tt+2,tt);

        t = theParent->GetMaxW();
        v = r.right-r.left+1+ttv;

        if (t>v) {
            AdjustScroller (hScroll, v, t);
        } else {
            AdjustScroller (hScroll, t, t);
        }

        gtk_widget_set_sensitive(hScroll,t>v);
        gtk_fixed_move              (wPane, vScroll, r.right-ttv+1, r.top);
        gtk_widget_set_size_request (vScroll,ttv,r.bottom-r.top+2);

        t = theParent->GetMaxH();
        v = r.bottom-r.top+1+tt;
        if (t>v) {
            AdjustScroller (vScroll, v, t);
        } else {
            AdjustScroller (vScroll, t, t);
        }

        gtk_widget_set_sensitive(vScroll,t>v);
    }
    rel = r;
}

//__________________________________________________________________

void        _HYPlatformComponent::_Paint (Ptr)
{
    _HYComponent* parent = (_HYComponent*)this;
    if (parent->settings.width&(HY_COMPONENT_BORDER|HY_COMPONENT_WELL)) {
        GdkWindow* pWindow     = gtk_widget_get_parent_window (parentWindow);
        GdkGC    * tempGC      = gdk_gc_new (pWindow);

        GdkColor lineColor = HYColorToGDKColor ((_HYColor) {
            0x40,0x40,0x40
        });
        gdk_gc_set_line_attributes (tempGC, 1, GDK_LINE_SOLID, GDK_CAP_NOT_LAST, GDK_JOIN_MITER);

        long    hplus = parent->parentWindow->allocation.x,
                vplus = parent->parentWindow->allocation.y;

        if (parent->settings.width&HY_COMPONENT_BORDER) {
            gdk_gc_set_foreground      (tempGC, &lineColor);

            if (parent->settings.width&HY_COMPONENT_BORDER_T) {
                gdk_draw_line (pWindow, tempGC, rel.left+hplus,rel.top+vplus, rel.right-1+hplus,rel.top+vplus);
            }

            if (parent->settings.width&HY_COMPONENT_BORDER_R) {
                gdk_draw_line (pWindow, tempGC, rel.right-1+hplus,rel.top+vplus, rel.right-1+hplus,rel.bottom-1+vplus);
            }

            if (parent->settings.width&HY_COMPONENT_BORDER_B) {
                gdk_draw_line (pWindow, tempGC, rel.right-1+hplus,rel.bottom-1+vplus, rel.left+hplus,rel.bottom-1+vplus);
            }

            if (parent->settings.width&HY_COMPONENT_BORDER_L) {
                gdk_draw_line (pWindow, tempGC, rel.left+hplus,rel.bottom-1+vplus, rel.left+hplus,rel.top+vplus);
            }
        }

        if (parent->settings.width&HY_COMPONENT_WELL) {
            // TBI (the actual sunken thing)
            gdk_draw_rectangle(pWindow,tempGC,false,rel.left+2+hplus,rel.top+2+vplus,rel.right-rel.left-4+hplus,rel.bottom-rel.top-4+vplus);
        }
        g_object_unref (tempGC);
    }
}

//__________________________________________________________________

void      _HYPlatformComponent::_Update (Ptr p)
{
    _Paint (p);
}


//__________________________________________________________________

void        _HYPlatformComponent::_Activate (void)
{
    activationFlag = true;
}

//__________________________________________________________________

void        _HYPlatformComponent::_Deactivate (void)
{
    //printf ("De-act comp\n");
    activationFlag = false;
}


//__________________________________________________________________

bool _HYPlatformComponent::_ProcessOSEvent (Ptr vEvent)
{
    // TBI?
    return false;
}


//__________________________________________________________________


_HYRect _HYPlatformComponent::_VisibleContents (Ptr p)
{
    _HYComponent* theParent = (_HYComponent*)this;

    _HYRect     * r = (_HYRect*)p,
                  res;

    _Parameter    v;

    short windowW   =   r->right-r->left,
          windowH =   r->bottom-r->top;

    if ((!hScroll)&&(!vScroll)) {
        res.left    = theParent->hOrigin;
        res.right   = res.left+windowW;
        res.top     = theParent->vOrigin;
        res.bottom  = res.top+windowH;
        return res;
    }

    if (hScroll) {
        windowH-=hScroll->allocation.height>1?hScroll->allocation.height:HY_SCROLLER_WIDTH;
    }
    if (vScroll) {
        windowW-=vScroll->allocation.width>1?vScroll->allocation.width:HY_SCROLLER_WIDTH;
    }

    if (hScroll) {
        if (windowW>theParent->GetMaxW()) {
            res.left = 0;
            res.right = theParent->GetMaxW();
        } else {
            v = _GetHScrollerPos();
            res.left = (theParent->GetMaxW()-windowW)*v/(double)MAX_CONTROL_VALUE;
            res.right = res.left+windowW;
        }
    } else {
        res.left = 0;
        res.right = windowW;
    }
    if (vScroll) {
        if (windowH>theParent->GetMaxH()) {
            res.top = 0;
            res.bottom = theParent->GetMaxH();
        } else {
            v = _GetVScrollerPos ();
            res.top = (theParent->GetMaxH()-windowH)*v/(double)MAX_CONTROL_VALUE;
            res.bottom = res.top+windowH;
        }
    } else {
        res.top = 0;
        res.bottom = windowH;
    }
    return res;
}

//__________________________________________________________________

void    _HYCanvas::_Paint (Ptr p)
{
    _HYRect * destR = (_HYRect*)p;

    GdkRectangle       destRect,
                       srcRect;

    destRect.width  = destR->right;
    destRect.height = destR->bottom;

    if (HasHScroll()) {
        destRect.height-=hScroll->allocation.height;
    }
    if (HasVScroll()) {
        destRect.width-=vScroll->allocation.width;
    }

    destRect.width  -= destR->left-1;
    destRect.height -= destR->top-1;
    destRect.x  = destR->left+parentWindow->allocation.x;
    destRect.y  = destR->top+parentWindow->allocation.y;

    srcRect = HYRect2GDKRect(_VisibleContents (p));

//  printf ("Canvas paint %d %d %d %d\n", w, h, destRect.width, destRect.height);

    gdk_draw_drawable (GDK_DRAWABLE(parentWindow->window), theContext, thePane, srcRect.x, srcRect.y,
                       destRect.x, destRect.y, destRect.width, destRect.height);

    _HYPlatformComponent::_Paint(p);
}

//__________________________________________________________________

void    _HYCanvas::_Update (Ptr p)
{
    _Paint(p);
}


//__________________________________________________________________

bool _HYCanvas::_ProcessOSEvent (Ptr vEvent)
{
    _HY_GTK_UI_Message *theMessage = (_HY_GTK_UI_Message*)vEvent;

    gdouble   xc,
              yc;

    gdk_event_get_coords (theMessage->theEvent,&xc,&yc);

    switch (theMessage->theEvent->type) {
    case GDK_BUTTON_PRESS: {
        GdkEventButton * be = (GdkEventButton*)theMessage->theEvent;

        if (be->button > 1) {
            if (exportFormats.lLength==0) {
                findGraphicsExporterComponents (exportFormats);
            }

            _String s1 ("Save as a picture"),
                    s2 ("Save canvas as:"),
                    filePath;
            _List   menuOptions;
            menuOptions && & s1;
            long    menuChoice;
            s1 = HandlePullDown (menuOptions,xc,yc,0);
            menuChoice  = menuOptions.Find (&s1);
            if (menuChoice==0) {
                _SavePicture (s2);
            }

            return true;
        } else {
            if (messageRecipient && doMouseClicks) {
                messageRecipient->ProcessEvent(generateContextPopUpEvent (GetID(),
                                               xc-parentWindow->allocation.x-rel.left,
                                               yc-parentWindow->allocation.y-rel.top));
                return true;
            }
        }
    }
    }
    return false;
}

