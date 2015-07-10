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
    A  window object - GTK glue.

    Sergei L. Kosakovsky Pond, Fall 2004.
*/

#include "HYWindow.h"
#include "HYEventTypes.h"
#include "HYUtils.h"
#include "HYPlatformGraphicPane.h"
#include "HYConsoleWindow.h"
#include <gtk/gtk.h>
#include <gdk-pixbuf/gdk-pixbuf.h>


//__________________________________________________________________
extern      _SimpleList windowPtrs,
            windowObjects;

extern      _String     objectInspectorTitle;

// GTK callbacks

//__________________________________________________________________

void hyphy_menu_item_callback(gpointer  data, guint menuItem, GtkWidget *widget)
{
    _HYWindow * parentWindow = (_HYWindow*)(data);
    parentWindow->_ProcessMenuSelection (menuItem);
}

//__________________________________________________________________

static gboolean window_event_callback( GtkWidget *widget, GdkEvent* theEvent, gpointer   data )
{
    if (theEvent) {
        _HY_GTK_UI_Message theMessage = {FALSE,theEvent};
        ((_HYWindow*)data)->_ProcessOSEvent((Ptr)&theMessage);
        return theMessage.processingResult;
    }
    return false;
}

//__________________________________________________________________

static gboolean window_expose_callback( GtkWidget *widget,  GdkEventExpose* theEvent, gpointer   data )
{
    ((_HYWindow*)data)->_PaintHook ((Ptr)theEvent);
    return FALSE;
}

//__________________________________________________________________

static gboolean window_noexpose_callback( GtkWidget *widget,  GdkEventNoExpose* theEvent, gpointer   data )
{
    ((_HYWindow*)data)->_Paint ((Ptr)theEvent);
    return TRUE;
}



//__________________________________________________________________

static gboolean window_expose_callback_null ( GtkWidget *,  GdkEventExpose*, gpointer  )
{
    return TRUE;
}

//__________________________________________________________________

static void activate_window_callback( GtkWidget *widget, gpointer   data )
{
    _HYWindow * parentWindow = (_HYWindow*)(data);
    parentWindow->Activate();
}

//__________________________________________________________________

void h_scroll_bar_callback_window (GtkRange *widget, gpointer data)
{
    GtkAdjustment* ta = gtk_range_get_adjustment (widget);
    double newV  = ta->value*(ta->upper/(ta->upper-ta->page_size));
    _HYWindow * parent_obj = (_HYWindow*)data;
    parent_obj->ProcessEvent (generateScrollEvent(newV-parent_obj->last_H_Position,0));
    parent_obj->last_H_Position = newV;
}

//__________________________________________________________________

void v_scroll_bar_callback_window (GtkRange *widget, gpointer data)
{
    GtkAdjustment* ta = gtk_range_get_adjustment (widget);
    double newV  = ta->value*(ta->upper/(ta->upper-ta->page_size));
    _HYWindow * parent_obj = (_HYWindow*)data;
    parent_obj->ProcessEvent (generateScrollEvent(0,newV-parent_obj->last_V_Position));
    parent_obj->last_V_Position = newV;
}

//__________________________________________________________________

static void window_resize_callback (GtkWidget *window, GtkAllocation* allocation, gpointer data)
{
    _HYWindow * parent_obj = (_HYWindow*)data;
    //parent_obj->printf ("Size-allocate: %x %d %d %d %d\n",window, allocation->width, allocation->height,parent_obj->windowContent->allocation.width,parent_obj->windowContent->allocation.height);
    if (parent_obj->last_V_Position<0.) {
        parent_obj->last_V_Position = 0.0;
        parent_obj->SetContentSize (parent_obj->contentWidth, parent_obj->contentHeight);
    }
    if (allocation->width != parent_obj->lastWW || allocation->height != parent_obj->lastWH) {
        parent_obj->lastWW = allocation->width;
        parent_obj->lastWH = allocation->height;
        parent_obj->SetWindowRectangle (0,0,allocation->height, allocation->width,false);
    }
}

//__________________________________________________________________

void AdjustScroller (GtkWidget* scrollbar, long viewport, long size)
{
    GtkAdjustment* old = gtk_range_get_adjustment (GTK_RANGE(scrollbar));

    double         newPageSize    = viewport*MAX_CONTROL_VALUE/size,
                   scaledValue =  old->upper-old->page_size?old->value*(old->upper/(old->upper-old->page_size)):0.0;

    scaledValue *= (old->upper-newPageSize)/old->upper;
    GtkObject* newAdj = gtk_adjustment_new (scaledValue,0,MAX_CONTROL_VALUE,MAX_CONTROL_VALUE/100.0,MAX_CONTROL_VALUE/10.0,newPageSize);
    gtk_range_set_adjustment (GTK_RANGE(scrollbar),GTK_ADJUSTMENT(newAdj));
}

//__________________________________________________________________
GdkColor    HYColorToGDKColor   (_HYColor hc)
{
    GdkColor gc;

    if (hc.R+hc.B+(long)hc.G==765) {
        gc.pixel = 0x00ffffff;
        gc.red   =
            gc.green =
                gc.blue  = 0xffff;
    } else {
        gc.pixel = hc.R*65536L+hc.G*256L+hc.B;
        gc.red   = hc.R*256;
        gc.green = hc.G*256;
        gc.blue  = hc.B*256;
    }

    return gc;
}


//__________________________________________________________________
void        _PaintTheCircle (GdkPixbuf * theIcon, GtkWidget * theWindow)
{
    if (GTK_WIDGET_MAPPED (theWindow) && GDK_IS_PIXBUF (theIcon))
        gdk_draw_pixbuf (theWindow->window, nil, theIcon, 0,0,
                         theWindow->allocation.x + 3,
                         theWindow->allocation.y+theWindow->allocation.height - 13,
                         -1,-1, GDK_RGB_DITHER_NONE, 0, 0);
}


//__________________________________________________________________

_HYPlatformWindow::_HYPlatformWindow(unsigned char windowFlag,_String windowTitle,bool windowVisibility, Ptr theParent)
{
    lastWW = lastWH = 0;
    last_H_Position = 0.0;
    last_V_Position = -1.0;

    if (windowFlag & HY_WINDOW_SHEET) {
        windowFlag |= HY_WINDOW_DLOG;
    }

    theWindow = gtk_window_new (GTK_WINDOW_TOPLEVEL);

    // debug color set here
    /*_HYColor red = {255,0,0};
    GdkColor redGDK = HYColorToGDKColor(red);
    gtk_widget_modify_bg (theWindow, GTK_STATE_NORMAL, &redGDK);*/

    gtk_window_set_title (GTK_WINDOW (theWindow), windowTitle.sData);
    gtk_container_set_border_width (GTK_CONTAINER (theWindow), 0);

    gtk_window_set_default_size (GTK_WINDOW(theWindow),300,250);
    gtk_window_set_resizable    (GTK_WINDOW(theWindow),windowFlag & HY_WINDOW_SIZE);
    gtk_widget_set_events       (theWindow,GDK_ALL_EVENTS_MASK-GDK_POINTER_MOTION_HINT_MASK);
    //gtk_widget_set_events     (theWindow,GDK_POINTER_MOTION_MASK);

    g_signal_connect (G_OBJECT (theWindow), "event",            G_CALLBACK (window_event_callback),     (_HYWindow*)this);
    g_signal_connect (G_OBJECT (theWindow), "activate-default", G_CALLBACK (activate_window_callback),  (_HYWindow*)this);
    //g_signal_connect (G_OBJECT (theWindow), "expose-event",       G_CALLBACK (window_expose_callback),    (_HYWindow*)this);
    // a(and other menu user UI events)

    if (theParent) {
        gtk_window_set_transient_for (GTK_WINDOW(theWindow), GTK_WINDOW((GtkWidget*)((_HYWindow*)theParent)->theWindow));
        gtk_window_set_destroy_with_parent (GTK_WINDOW(theWindow),true);
    }

    windowContent = gtk_fixed_new ();
    g_signal_connect (G_OBJECT (windowContent), "expose-event", G_CALLBACK (window_expose_callback),    (_HYWindow*)this);
    //g_signal_connect (G_OBJECT (windowContent), "window-state-event", G_CALLBACK (null_window_state), nil);
    gtk_widget_set_app_paintable(windowContent,TRUE);
    gtk_container_set_resize_mode (GTK_CONTAINER(windowContent),GTK_RESIZE_IMMEDIATE);

    menu_items = nil;

    if (windowFlag & HY_WINDOW_DLOG) {
        windowMB = nil;
    } else {
        windowMB = theWindow;
        ((_HYWindow*)this)->_SetMenuBar  ();
    }

    gtk_widget_show (windowContent);

    if ((windowFlag&HY_WINDOW_SCROLL)&&(windowFlag&HY_WINDOW_SIZE)) {
        vbox = gtk_vbox_new (FALSE, 0);
        hbox = gtk_hbox_new (FALSE, 0);

        GtkObject
        *v_adj = gtk_adjustment_new(0,0,MAX_CONTROL_VALUE,MAX_CONTROL_VALUE/100.0,MAX_CONTROL_VALUE/10.0,MAX_CONTROL_VALUE/5.0),
         *h_adj = gtk_adjustment_new(0,0,MAX_CONTROL_VALUE,MAX_CONTROL_VALUE/100.0,MAX_CONTROL_VALUE/10.0,MAX_CONTROL_VALUE/5.0);

        vScroll = gtk_vscrollbar_new  ((GtkAdjustment*)v_adj),
        hScroll = gtk_hscrollbar_new  ((GtkAdjustment*)h_adj);

        /*GtkWidget * tScroll = gtk_vscrollbar_new ((GtkAdjustment*)t_adj);
        gtk_widget_set_size_request(tScroll,20,100);

        gtk_fixed_put (GTK_FIXED(windowContent), tScroll, 30,40);
        gtk_widget_show (tScroll);*/

        gtk_widget_show (hScroll);
        gtk_widget_show (vScroll);
        gtk_widget_show (hbox);
        gtk_widget_show (vbox);

        g_signal_connect (G_OBJECT (vScroll), "value-changed",G_CALLBACK (v_scroll_bar_callback_window), (gpointer)((_HYWindow*)this));
        g_signal_connect (G_OBJECT (hScroll), "value-changed",G_CALLBACK (h_scroll_bar_callback_window), (gpointer)((_HYWindow*)this));

        gtk_box_pack_start (GTK_BOX (hbox), windowContent, TRUE, TRUE, 0);
        gtk_box_pack_end   (GTK_BOX (hbox), vScroll, FALSE, FALSE, 0);

        if (windowMB) {
            gtk_box_pack_start (GTK_BOX (vbox), windowMB, FALSE, FALSE, 0);
        }
        gtk_box_pack_start (GTK_BOX (vbox), hbox, TRUE, TRUE, 0);
        gtk_box_pack_end (GTK_BOX (vbox), hScroll, FALSE, FALSE, 0);
        gtk_container_add (GTK_CONTAINER (theWindow), vbox);
    } else {
        if (windowMB) {
            hScroll = vScroll = hbox = nil;
            vbox = gtk_vbox_new (FALSE, 0);
            gtk_box_pack_start   (GTK_BOX (vbox), windowMB, FALSE, FALSE, 0);
            gtk_box_pack_end (GTK_BOX (vbox), windowContent, TRUE, TRUE, 0);
            gtk_widget_show (vbox);
            gtk_container_add (GTK_CONTAINER (theWindow), vbox);
        } else {
            hScroll = vScroll = hbox = vbox = nil;
            windowFlag &= 123;
            gtk_container_add (GTK_CONTAINER (theWindow), windowContent);
        }
    }

    windowPtrs   << (long)theWindow;
    windowObjects<< (long)this;

    gtk_window_set_type_hint (GTK_WINDOW(theWindow),windowFlag & HY_WINDOW_DLOG ? GDK_WINDOW_TYPE_HINT_DIALOG : GDK_WINDOW_TYPE_HINT_NORMAL);
    gtk_window_set_modal        (GTK_WINDOW(theWindow),windowFlag & HY_WINDOW_DLOG);
    if (windowVisibility) {
        gtk_widget_show(theWindow);
    }

    flags = windowFlag;

    g_signal_connect (G_OBJECT (theWindow), "size-allocate",    G_CALLBACK (window_resize_callback),    (_HYWindow*)this);
}

//__________________________________________________________________

_HYPlatformWindow::~_HYPlatformWindow(void)
{
    if (menu_items) {
        gtk_object_sink (GTK_OBJECT(menu_items));
    }
}

//__________________________________________________________________

void _HYPlatformWindow::_SetTitle(const _String & windowTitle)
{
    gtk_window_set_title (GTK_WINDOW (theWindow), windowTitle.sData);
}

//__________________________________________________________________

void _HYPlatformWindow::_Show(void)
{
    bool firstShow = !GTK_WIDGET_MAPPED (theWindow);
    gtk_window_present(GTK_WINDOW (theWindow));
    if (firstShow) {
        _HYWindow *pw = ((_HYWindow*)this);
        pw->SetWindowRectangle (pw->top,pw->left,pw->bottom,pw->right,true);
    }
}

//__________________________________________________________________

_HYRect _HYPlatformWindow::_GetWindowRect(void)
{
    _HYRect res ;

    gint    x,y;
    gtk_window_get_position(GTK_WINDOW (theWindow),&x,&y);
    res.left = x;
    res.top  = y;
    gtk_window_get_size(GTK_WINDOW (theWindow),&x,&y);
    res.right = x+res.left-1;
    res.bottom = y+res.top-1;

    return res;
}


//__________________________________________________________________

void _HYPlatformWindow::_Hide(void)
{
    gtk_widget_hide (GTK_WIDGET(theWindow));
}

//__________________________________________________________________

_String& _HYPlatformWindow::_GetTitle(void)
{
    return      ((_HYWindow*)this)->GetTitle();
}

//__________________________________________________________________

void _HYPlatformWindow::_SetPosition(int l,int t)
{
    gtk_window_move(GTK_WINDOW (theWindow),l,t);

    long deltah = savedLoc.left-l,
         deltav = savedLoc.top-t;
    savedLoc.left = l;
    savedLoc.right -= deltah;
    savedLoc.top = t;
    savedLoc.bottom -= deltav;
}

//__________________________________________________________________

bool _HYPlatformWindow::_Close(Ptr theData)
{
    _HYWindow* theParent = (_HYWindow*)this;
    bool       doit = theParent->ConfirmClose();

    if (doit) {
        long f = windowObjects.Find((long)this);
        if (f>=0) {
            windowObjects.Delete(f);
            windowPtrs.Delete(f);
            windowObjectRefs.Delete(f);
        }
        if (!theData) {
            gtk_object_destroy (GTK_OBJECT(theWindow));
        }

        if (windowPtrs.lLength == 0) {
            gtk_main_quit();
        }
    }

    return doit;
}


//__________________________________________________________________

void _HYPlatformWindow::_Activate(void)
{
    if (!GTK_WIDGET_REALIZED (theWindow)) {
        _Show();
    }

    _SetMenuBar ();
}

//__________________________________________________________________

void _HYPlatformWindow::_BringWindowToFront(void)
{
    _Show();
}

//__________________________________________________________________

void _HYPlatformWindow::_Deactivate(void)
{
    _UnsetMenuBar();
}

static GtkItemFactoryEntry hyphy_window_menu_file[] = {
    { "/_File",         NULL,         NULL,           0, "<Branch>" },
    { "/File/_Save",     "<control>S", (GtkItemFactoryCallback)hyphy_menu_item_callback,    HY_WINDOW_MENU_ID_FILE+1, "<StockItem>", GTK_STOCK_SAVE },
    { "/File/_Print",    "<control>P", (GtkItemFactoryCallback)hyphy_menu_item_callback,    HY_WINDOW_MENU_ID_FILE+2, "<StockItem>", GTK_STOCK_OPEN },
    { "/File/_Close",    "<control>W", (GtkItemFactoryCallback)hyphy_menu_item_callback,    HY_WINDOW_MENU_ID_FILE, "<StockItem>", GTK_STOCK_CLOSE },
    { "/File/sep1",     NULL,         NULL,           0,                    "<Separator>" },
    { "/File/S_witch to console", "<control>0",         (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_WINDOW_MENU_ID_FILE-1, "<Item>" },
    { "/File/Object _inspector", "<control>I",         (GtkItemFactoryCallback)hyphy_menu_item_callback,  HY_WINDOW_MENU_ID_FILE-2, "<Item>" }
};
static GtkItemFactoryEntry hyphy_window_menu_edit[] = {
    { "/_Edit",         NULL,         NULL,           0, "<Branch>" },
    { "/Edit/_Undo",     "<control>Z", (GtkItemFactoryCallback)hyphy_menu_item_callback,    HY_WINDOW_MENU_ID_EDIT, "<StockItem>", GTK_STOCK_UNDO },
    { "/Edit/sep1",     NULL,         NULL,           0,                    "<Separator>" },
    { "/Edit/_Copy",    "<control>C", (GtkItemFactoryCallback)hyphy_menu_item_callback,    HY_WINDOW_MENU_ID_EDIT+1, "<StockItem>", GTK_STOCK_COPY },
    { "/Edit/_Cut",    "<control>X", (GtkItemFactoryCallback)hyphy_menu_item_callback,    HY_WINDOW_MENU_ID_EDIT+2, "<StockItem>", GTK_STOCK_CUT },
    { "/Edit/_Paste", "<control>V",         (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_WINDOW_MENU_ID_EDIT+3, "<StockItem>",GTK_STOCK_PASTE },
    { "/Edit/_Clear", NULL ,         (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_WINDOW_MENU_ID_EDIT+4, "<StockItem>",GTK_STOCK_CLEAR },
    { "/Edit/Select _All", "<control>A",         (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_WINDOW_MENU_ID_EDIT+5, "<Item>" },
};

//__________________________________________________________________

bool _HYPlatformWindow::_CleanDefaultMenu (void)
{
    if (windowMB && menu_items) {
        GtkWidget * oiw = gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/File/Object Inspector");
        if (oiw) {
            gtk_item_factory_delete_entries (menu_items, sizeof (hyphy_window_menu_file) / sizeof (hyphy_window_menu_file[0]), hyphy_window_menu_file);
            gtk_item_factory_delete_entries (menu_items, sizeof (hyphy_window_menu_edit) / sizeof (hyphy_window_menu_edit[0]), hyphy_window_menu_edit);
            return true;
        }

    }
    return false;
}


//__________________________________________________________________

void _HYPlatformWindow::_SetMenuBar(void)
{
    _HYWindow* theParent = (_HYWindow*)this;
    if (windowMB && menu_items == nil) {

        GtkAccelGroup * accel_group = gtk_accel_group_new();
        menu_items = gtk_item_factory_new (GTK_TYPE_MENU_BAR, "<HY_WINDOW>",accel_group);
        gtk_item_factory_create_items (menu_items,  sizeof (hyphy_window_menu_file) / sizeof (hyphy_window_menu_file[0]), hyphy_window_menu_file, theParent);
        gtk_item_factory_create_items (menu_items,  sizeof (hyphy_window_menu_edit) / sizeof (hyphy_window_menu_edit[0]), hyphy_window_menu_edit, theParent);
        gtk_window_add_accel_group (GTK_WINDOW (theWindow), accel_group);
        windowMB   = gtk_item_factory_get_widget (menu_items, "<HY_WINDOW>");
        gtk_widget_show (windowMB);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget (menu_items,"<HY_WINDOW>/Edit/Undo"),  FALSE);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget (menu_items,"<HY_WINDOW>/Edit/Copy"),  FALSE);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget (menu_items,"<HY_WINDOW>/Edit/Cut"),  FALSE);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget (menu_items,"<HY_WINDOW>/Edit/Paste"),  FALSE);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget (menu_items,"<HY_WINDOW>/Edit/Clear"),  FALSE);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget (menu_items,"<HY_WINDOW>/Edit/Select All"),  FALSE);

    } else if (menu_items) {
        gtk_widget_set_sensitive(gtk_item_factory_get_widget (menu_items,"<HY_WINDOW>/File/Save"),  theParent->IsSaveEnabled());
        gtk_widget_set_sensitive(gtk_item_factory_get_widget (menu_items,"<HY_WINDOW>/File/Print"), theParent->IsPrintEnabled());
    }
}

//__________________________________________________________________

void _HYPlatformWindow::_UnsetMenuBar(void)
{

}

//__________________________________________________________________

void _HYPlatformWindow::_Paint(Ptr)
{

}

//__________________________________________________________________

void _HYPlatformWindow::_PaintHook (Ptr p)
{
    ((_HYWindow*)this)->Paint(p);
}


//__________________________________________________________________

void _HYPlatformWindow::_Update(Ptr)
{

}

//__________________________________________________________________

void _HYPlatformWindow::_SetWindowRectangle(int top, int left, int bottom, int right, bool ss)
{
    if (ss) {
        _HYWindow *pW = (_HYWindow*)this;
        long    menuHeight = 0;
        if (pW->windowMB) {
            menuHeight = pW->windowMB->allocation.x;
            if (menuHeight <= 0) {
                GtkRequisition sizeReq;
                gtk_widget_size_request (pW->windowMB,&sizeReq);
                menuHeight = sizeReq.height;
            }
        }
        gtk_window_resize (GTK_WINDOW(theWindow), right-left+1,bottom-top+1+menuHeight);
        //if (left || top)
        //  gtk_window_move(GTK_WINDOW(theWindow),left,top);
    } else {
        long      newSize;
        if (hScroll && vScroll) {
            _HYWindow* theParent = (_HYWindow*)this;
            int windowWidth  = windowContent->allocation.width,
                windowHeight = windowContent->allocation.height;

            gtk_widget_set_sensitive(hScroll,windowWidth < theParent->contentWidth);
            AdjustScroller (hScroll, windowWidth, theParent->contentWidth);

            gtk_widget_set_sensitive(vScroll,windowHeight < theParent->contentHeight);
            AdjustScroller (vScroll, windowHeight, theParent->contentHeight);
        }
        savedLoc.bottom=savedLoc.top+bottom-top;
        savedLoc.right=savedLoc.left+right-left;
    }
}

//__________________________________________________________________

bool _HYPlatformWindow::_ProcessOSEvent (Ptr vEvent)
{
    if (vEvent) {
        _HYWindow          *theParent = (_HYWindow*)this;
        _HY_GTK_UI_Message *theMessage = (_HY_GTK_UI_Message*)vEvent;

        switch (theMessage->theEvent->type) {
        case GDK_DELETE:
            theMessage->processingResult = (theParent->Close(nil)==false);
            return true;

        case GDK_FOCUS_CHANGE:
            if (((GdkEventFocus*)(theMessage->theEvent))->in) {
                theParent->Activate();
            }
            return true;

        case GDK_BUTTON_PRESS:
            switch (((GdkEventButton*)(theMessage->theEvent))->button) {
            case 1:
                //printf ("Left button\n");
                break;
            case 2:
                //printf ("Middle button\n");
                break;
            case 3:
                //printf ("Right button\n");
                /*{
                    _List       testMenu;
                    _String     item ("George");
                    testMenu && & item;
                    item = "Bush";
                    testMenu && & item;
                    item = "SEPARATOR";
                    testMenu && & item;
                    item = "Sucks Ass";
                    testMenu && & item;
                    printf ("Menu Selection %s\n", HandlePullDown(testMenu,100,100,0).sData);
                }*/

                break;
            }
            return false;

            //case GDK_CONFIGURE:
            //{
            //  GdkEventConfigure * ci = (GdkEventConfigure*)theMessage->theEvent;
            //  gtk_container_resize_children (GTK_CONTAINER(theWindow));
            //  _SetWindowRectangle (0,0,ci->height, ci->width,false);
            //  return true;
            //}
        }
    }
    return false;
}
//__________________________________________________________________

void    _HYPlatformWindow::_SetContentSize (int w, int h)
{
    if (hScroll == nil) {
        return;
    }

    _HYWindow* theParent = (_HYWindow*)this;

    int     windowW =  windowContent->allocation.width,
            windowH =  windowContent->allocation.height;

    GdkGeometry windowG;
    windowG.max_width  = w+vScroll->allocation.width;
    windowG.max_height = h+hScroll->allocation.height+(windowMB?windowMB->allocation.height:0);
    gtk_window_set_geometry_hints (GTK_WINDOW(theWindow), NULL, &windowG, (GdkWindowHints)(GDK_HINT_MAX_SIZE));

    //printf ("Window Size Allocation %d %d %d\n", windowG.max_width, windowG.max_height, windowMB->allocation.height);


    gtk_widget_set_sensitive(hScroll,w>windowW);
    AdjustScroller (hScroll, windowW, w);

    gtk_widget_set_sensitive(vScroll,h>windowH);
    AdjustScroller (vScroll, windowH, h);

}

//__________________________________________________________________


void    _HYPlatformWindow::_VisibleContents (int& t,int& l,int& b,int& r)
{
    _HYWindow*   theParent = (_HYWindow*)this;

    if (GTK_WIDGET_MAPPED (windowContent)) {
        if (hScroll) {
            _Parameter v;

            long  windowH   = windowContent->allocation.height,
                  windowW   = windowContent->allocation.width;

            if (windowW>theParent->contentWidth) {
                l = theParent->left;
                r = theParent->right;
            } else {
                GtkAdjustment* old = gtk_range_get_adjustment (GTK_RANGE(hScroll));
                //printf ("GtkAdjustment %g %g %g\n", old->upper, old->page_size, old->value);
                if (old->upper-old->page_size) {
                    v = old->value*(old->upper/(old->upper-old->page_size));
                    l = theParent->left+(theParent->contentWidth-windowW)*v/(double)MAX_CONTROL_VALUE;
                } else {
                    l = theParent->left;
                }
                r = l+windowW;
            }
            if (windowH>theParent->contentHeight) {
                t = theParent->top;
                b = theParent->bottom;
            } else {
                GtkAdjustment* old = gtk_range_get_adjustment (GTK_RANGE(vScroll));
                if (old->upper-old->page_size) {
                    v = old->value*(old->upper/(old->upper-old->page_size));
                    t = theParent->top+(theParent->contentHeight-windowH)*v/(double)MAX_CONTROL_VALUE;
                } else {
                    t = theParent->top;
                }
                b = t+windowH;
            }
            //printf ("Visible area %d %d %d %d %d %d\n", windowW, windowH, l, r, t, b);
        } else {
            t = theParent->top;
            l = theParent->left;
            b = windowContent->allocation.height+t;
            r = windowContent->allocation.width+l;
        }
    } else {
        t = theParent->top;
        l = theParent->left;
        b = theParent->bottom;
        r = theParent->right;
    }
}


//__________________________________________________________________

void    _HYPlatformWindow::_SetWindowBackColor (_HYColor newColor)
{
    /*Fl_Color new_color = fl_rgb_color (newColor.R, newColor.G, newColor.B);
    color (new_color);
    redraw();*/
}

//__________________________________________________________________

bool    _HYPlatformWindow::_IsHScroll (GtkWidget* ch)
{
    return ch == hScroll;
}


//__________________________________________________________________

bool        _HYWindow::_ProcessMenuSelection (long msel)
{
    switch (msel) {
    case HY_WINDOW_MENU_ID_FILE-2: { // show object inspector
        ShowObjectInspector();
        return true;
    }
    case HY_WINDOW_MENU_ID_FILE-1: { // switch to console
        if (hyphyConsoleWindow) {
            hyphyConsoleWindow->BringToFront();
        }
        return true;
    }
    case HY_WINDOW_MENU_ID_FILE: { // close window
        Close (nil);
        return true;
    }
    }
    return false;
}

//__________________________________________________________________

void _HYPlatformWindow::_SelectWindow (void)
{
    _Show();
}

//EOF
