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
    GTK+ Portions of the chart window class

    Sergei L. Kosakovsky Pond, March 2005.
*/

#include "HYChartWindow.h"
#include "HYCanvas.h"
#include "HYUtils.h"
#include "HYPullDown.h"
#include "HYDialogs.h"

#include "math.h"

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

extern   _Parameter                pi_const;

#define  HY_CHART_WIN32_MENU_BASE   6000
#define  HY_CHARTD_WIN32_MENU_BASE  27000

static GtkItemFactoryEntry hyphy_char_window_menu[] = {
    { "/File/Save _Graphic", "<control><alt>S", (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_WINDOW_MENU_ID_FILE+3, "<Item>"},
    { "/File/Save _Table", "<control><shift>S", (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_WINDOW_MENU_ID_FILE+4, "<Item>"},
    { "/File/Pri_nt Table", "<control><shift>S", (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_WINDOW_MENU_ID_FILE+5, "<Item>"},
    { "/_Chart",            NULL,         NULL,           0,                    "<Branch>" },
    { "/Chart/Chart _Name", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_CHART_WIN32_MENU_BASE+4, "<Item>"},
    { "/Chart/Chart _Options", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_CHART_WIN32_MENU_BASE, "<Item>"},
    { "/Chart/_Fonts", NULL,         NULL,           0,                 "<Branch>" },
    { "/Chart/Fonts/_Tickmark Font", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_CHART_WIN32_MENU_BASE+1, "<Item>"},
    { "/Chart/Fonts/_Legend Font", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_CHART_WIN32_MENU_BASE+2, "<Item>"},
    { "/Chart/Fonts/_Axis Label Font", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_CHART_WIN32_MENU_BASE+3, "<Item>"},
    { "/Chart/sep1",            NULL,         NULL,           0,        "<Separator>" },
    { "/Chart/_Data Processing",            NULL,         NULL,           0,                    "<Branch>" }
};

static GtkItemFactoryEntry hyphy_distro_window_menu[] = {
    { "/Cate_gories", NULL, NULL, 0, "<Branch>"},
    { "/Categories/Define new _variable", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_CHARTD_WIN32_MENU_BASE, "<Item>"},
    { "/Categories/_Delete variable", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_CHARTD_WIN32_MENU_BASE+1, "<Item>"},
    { "/Categories/_Conditional Distribution", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_CHARTD_WIN32_MENU_BASE+2, "<Item>"},
    { "/Categories/sep1",           NULL,         NULL,           0,        "<Separator>" }
};

//__________________________________________________________________

void _HYChartWindow::_SetMenuBar(void)
{
    if (menu_items && !gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Chart")) {
        gtk_item_factory_create_items (menu_items,  sizeof (hyphy_char_window_menu) / sizeof (hyphy_char_window_menu[0]), hyphy_char_window_menu, this);

        GtkMenu *fileMenu = GTK_MENU(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/File"));
        gtk_menu_reorder_child(fileMenu,gtk_item_factory_get_widget_by_action(menu_items,HY_WINDOW_MENU_ID_FILE+3),1);
        gtk_menu_reorder_child(fileMenu,gtk_item_factory_get_widget_by_action(menu_items,HY_WINDOW_MENU_ID_FILE+4),2);
        gtk_menu_reorder_child(fileMenu,gtk_item_factory_get_widget_by_action(menu_items,HY_WINDOW_MENU_ID_FILE+5),4);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget_by_action(menu_items,HY_WINDOW_MENU_ID_EDIT+1),true);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget_by_action(menu_items,HY_WINDOW_MENU_ID_EDIT+3),true);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget_by_action(menu_items,HY_WINDOW_MENU_ID_EDIT+5),true);

        /*InsertMenu    (saveMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_WINDOW_MENU_ID_FILE+1, "Save &Chart\tCtrl-S");
        InsertMenu      (saveMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_WINDOW_MENU_ID_FILE+3, "Save &Graphic");
        InsertMenu      (saveMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_WINDOW_MENU_ID_FILE+4, "Save &Table");

        InsertMenu      (printMenu, 0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_WINDOW_MENU_ID_FILE+2, "Print &Graphic\tCtrl-P");
        InsertMenu      (printMenu, 0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_WINDOW_MENU_ID_FILE+5, "Print &Data");*/



        for (long k=0; k<chartProcessors.lLength; k++) {
            _String *thisItem = (_String*)chartProcessors (k),
                     chopped = thisItem->Cut (thisItem->FindBackwards ('/',0,-1)+1,-1),
                     type = "<Item>";

            GtkItemFactoryEntry aProcEntry = {NULL,NULL,(GtkItemFactoryCallback)hyphy_menu_item_callback,HY_CHART_WIN32_MENU_BASE+5+k,type.sData};
            chopped = _String("/Chart/Data Processing/")&chopped;
            aProcEntry.path = chopped.sData;

            gtk_item_factory_create_items (menu_items,  1, &aProcEntry, this);
        }

        //ModifyMenu     (chartMenu, 0, MF_BYPOSITION|MF_POPUP, (UINT) saveMenu , "&Save");
        //ModifyMenu     (chartMenu, 1, MF_BYPOSITION|MF_POPUP, (UINT) printMenu , "&Print");

    }
}

//__________________________________________________________________

void _HYChartWindow::_UnsetMenuBar(void)
{

}

//__________________________________________________________________

void        _HYChartWindow::_PrintChart(void)
{
    _String printTBI ("Chart printing has not yet been implemented\n");
    ReportWarning (printTBI);
}


//__________________________________________________________________


bool        _HYChartWindow::_ProcessMenuSelection (long msel)
{
    switch (msel) {
    case HY_CHART_WIN32_MENU_BASE: { // chart menu
        HandleChartOptions ();
        return true;
    }
    case HY_WINDOW_MENU_ID_FILE+1: // save menu
    case HY_WINDOW_MENU_ID_FILE+3: // save menu
    case HY_WINDOW_MENU_ID_FILE+4: { // save menu
        DoSave ((msel==HY_WINDOW_MENU_ID_FILE-1)?0:msel-HY_WINDOW_MENU_ID_FILE-2);
        return true;
    }
    case HY_WINDOW_MENU_ID_FILE+2: // print menu
    case HY_WINDOW_MENU_ID_FILE+5: { // print menu
        DoPrint ((msel==HY_WINDOW_MENU_ID_FILE+2)?0:-1);
        return true;
    }
    case HY_CHART_WIN32_MENU_BASE+1: // font menu
    case HY_CHART_WIN32_MENU_BASE+2: // font menu
    case HY_CHART_WIN32_MENU_BASE+3: { // font menu
        DoChangeFont (msel-HY_CHART_WIN32_MENU_BASE-1);
        return true;
    }
    case HY_CHART_WIN32_MENU_BASE+4: { // chart name
        RenameChartWindow ();
        return true;
    }
    default: { // proc menu
        if (msel>=HY_CHART_WIN32_MENU_BASE+5) {
            ExecuteProcessor (msel-HY_CHART_WIN32_MENU_BASE-5);
            return true;
        }
    }
    }

    return _HYTWindow::_ProcessMenuSelection(msel);
}

//__________________________________________________________________

bool _HYChartWindow::_ProcessOSEvent (Ptr vEvent)
{
    static long   lastH = -1,
                  lastV = -1;

    if (!_HYTWindow::_ProcessOSEvent (vEvent)) {
        if (components.lLength == 0) {
            return false;
        }

        _HYPullDown *p1 = (_HYPullDown*)GetObject (4);

        if (p1&&(p1->GetSelection()>=8)&&(ySeries.lLength)) {

            _HY_GTK_UI_Message *theMessage = (_HY_GTK_UI_Message*)vEvent;

            gdouble   xc,
                      yc;

            gdk_event_get_coords (theMessage->theEvent,&xc,&yc);
            switch (theMessage->theEvent->type) {
            case GDK_BUTTON_PRESS: {
                if (((GdkEventButton*)theMessage->theEvent)->button != 1) {
                    return false;
                }

                lastH = xc-theWindow->allocation.x-windowContent->allocation.x;
                lastV = yc-theWindow->allocation.y-windowContent->allocation.y;

                if (FindClickedCell(lastH, lastV)!=0) { // the chart
                    lastH = -1;
                    lastV = -1;
                } else {
                    gdk_pointer_grab (theWindow->window,false,
                                      (GdkEventMask)(GDK_POINTER_MOTION_MASK | GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK),
                                      theWindow->window, NULL, ((GdkEventButton*)theMessage->theEvent)->time);
                    return      true;
                }
                break;
            }

            case GDK_BUTTON_RELEASE: {
                if (lastH>=0) {
                    gdk_pointer_ungrab (((GdkEventButton*)theMessage->theEvent)->time);
                    lastH = -1;
                    lastV = -1;
                    return  true;
                }
                break;
            }

            case GDK_MOTION_NOTIFY: {
                if (lastH>=0) {
                    long        newH = xc-theWindow->allocation.x-windowContent->allocation.x,
                                newV = yc-theWindow->allocation.y-windowContent->allocation.y;

                    bool        redraw = false;

                    _Parameter  stepper = pi_const/180.;

                    if (abs(newH-lastH)>abs(newV-lastV)) {
                        stepper *= 1+log (fabs(newH-lastH))/log(2.0);
                        if (newH-lastH<0) {
                            if (xyAngle>0.0) {
                                xyAngle -= stepper;
                                if (xyAngle<0) {
                                    xyAngle = 0;
                                }
                                redraw = true;
                            }
                        } else if (xyAngle<pi_const/2) {
                            xyAngle += stepper;
                            if (xyAngle>pi_const/2) {
                                xyAngle = pi_const/2;
                            }
                            redraw = true;
                        }
                    } else {
                        if (newV==lastV) {
                            return false;
                        }
                        stepper *= 1+log (fabs(newV-lastV))/log(2.0);
                        if (newV-lastV>0) {
                            if (zAngle<pi_const/2) {
                                zAngle += stepper;
                                if (zAngle>pi_const/2) {
                                    zAngle = pi_const/2;
                                }
                                redraw = true;
                            }
                        } else if (zAngle>0.0) {
                            zAngle -= stepper;
                            if (zAngle<0) {
                                zAngle = 0;
                            }
                            redraw = true;
                        }

                    }

                    if (redraw) {
                        ComputeProjectionSettings();
                        projectionMatrix = ComputeProjectionMatrix   ();
                        forceUpdateForScrolling = true;
                        DrawChart();
                        forceUpdateForScrolling = false;
                    }

                    lastH = newH;
                    lastV = newV;
                }
                break;
            }
            }
        }
        return false;
    }
    return true;
}

//__________________________________________________________________

void _HYChartWindow::_CopyChart (void)
{
    //_HYStretchCanvas    *sc = (_HYStretchCanvas*)GetObject (0);
    //PlaceBitmapInClipboard (sc->paneBitMap, theWindow);
}

//__________________________________________________________________

void _HYDistributionChartWindow::_SetMenuBar(void)
{
    if (menu_items && !gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Categories")) {
        _HYChartWindow::_SetMenuBar();
        gtk_item_factory_create_items (menu_items,  sizeof (hyphy_distro_window_menu) / sizeof (hyphy_distro_window_menu[0]), hyphy_distro_window_menu, this);

        if (distribProcessors.lLength > 0) {
            for (long k=0; k<distribProcessors.lLength; k++) {
                _String *thisItem = (_String*)distribProcessors (k),
                         chopped = thisItem->Cut (thisItem->FindBackwards ('/',0,-1)+1,-1),
                         type = "<Item>";

                GtkItemFactoryEntry aProcEntry = {NULL,NULL,(GtkItemFactoryCallback)hyphy_menu_item_callback,HY_CHARTD_WIN32_MENU_BASE+3+k,type.sData};
                chopped = _String("/Categories/")&chopped;
                aProcEntry.path = chopped.sData;

                gtk_item_factory_create_items (menu_items,  1, &aProcEntry, this);
            }
        }
    }
}

//__________________________________________________________________

void _HYDistributionChartWindow::_UnsetMenuBar(void)
{
    _HYChartWindow::_UnsetMenuBar();
}

//__________________________________________________________________

bool _HYDistributionChartWindow::_ProcessMenuSelection (long msel)
{
    switch (msel) {
    case HY_CHARTD_WIN32_MENU_BASE: {
        AddVariable ();
        return true;
    }
    case HY_CHARTD_WIN32_MENU_BASE+1: {
        RemoveVariable ();
        return true;
    }
    case HY_CHARTD_WIN32_MENU_BASE+2: {
        ShowMarginals ();
        return true;
    }
    default: {
        if (msel>=HY_CHARTD_WIN32_MENU_BASE+3) {
            HandleCatPostProcessor (msel-HY_CHARTD_WIN32_MENU_BASE-3);
            return true;
        }
    }
    }

    return _HYChartWindow::_ProcessMenuSelection(msel);
}



//EOF
