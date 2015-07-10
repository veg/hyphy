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
    GTK Portions of the console window class

    Sergei L. Kosakovsky Pond, Spring 2005.
*/

#include "HYConsoleWindow.h"
#include "HYTextBox.h"
#include "HYDialogs.h"
#include "HYUtils.h"
#include "HYEventTypes.h"
#include "time.h"
#include "batchlan.h"

#include "HYGraphicPane.h"
#include "HYDataPanel.h"
#include "HYSharedMain.h"

GdkColor        _DARKGREYBRUSH_    = HYColorToGDKColor((_HYColor)
{
    100,100,100
});
GdkColor        _WHITEBRUSH_       = HYColorToGDKColor((_HYColor)
{
    255,255,255
});

GdkGC*          _hyConsoleWindowGC = nil;
PangoLayout*    _hyConsoleLayout   = nil;

PangoFontDescription    *statusBarBold,
                        *statusBarNormal;

clock_t                 lastMeasure = 0;

#define                 RECENT_FILE_ITEMS 10


//__________________________________________________________________


extern  clock_t  timerStart,
        lastTimer;

void             displayAbout        (bool);
void             getUserFont         (void);

//__________________________________________________________________

static GtkItemFactoryEntry hyphy_console_window_menu_file[] = {
    { "/File/_New",             NULL,         NULL,         0, "<Branch>" },
    { "/File/New/New _Tree",      "<control>N",       (GtkItemFactoryCallback)hyphy_menu_item_callback,         70, "<Item>" },
    { "/File/New/New _Model",     "<control><mod>N",      (GtkItemFactoryCallback)hyphy_menu_item_callback,         71, "<Item>" },
    { "/File/New/New _Chart",     "<control><shift>N",  (GtkItemFactoryCallback)hyphy_menu_item_callback,           72, "<Item>" },
    {
        "/File/New/New _Genetic Code",
        NULL,         (GtkItemFactoryCallback)hyphy_menu_item_callback,         73, "<Item>"
    },
    {
        "/File/New/New S_QLite Database",
        NULL,         (GtkItemFactoryCallback)hyphy_menu_item_callback,         74, "<Item>"
    },

    { "/File/_Open",                NULL,         NULL,         0, "<Branch>" },
    { "/File/Open/Open _Batch File",      "<control>O",       (GtkItemFactoryCallback)hyphy_menu_item_callback,         80, "<Item>" },
    { "/File/Open/Open _Data File",     "<control><alt>O",    (GtkItemFactoryCallback)hyphy_menu_item_callback,         81, "<Item>" },
    { "/File/Open/Open _Tree File",     "<control><shift>O",    (GtkItemFactoryCallback)hyphy_menu_item_callback,           82, "<Item>" },
    { "/File/Open/sep1",            NULL,         NULL,           0,                    "<Separator>" },
    { "/File/Open/Open _Recent",            NULL,         NULL,           0,                    "<Branch>" },
    { "/File/Open/Open Recent/_Clear menu",      NULL,        (GtkItemFactoryCallback)hyphy_menu_item_callback,         1999, "<Item>" },
    { "/File/Open/Open Recent/sep1",            NULL,         NULL,           0,                    "<Separator>" },
    { "/File/Open/sep2",            NULL,         NULL,           0,                    "<Separator>" },
    { "/File/Open/Open Te_xt File",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,          83, "<Item>" },
    { "/File/Open/Open Ta_ble",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,              84, "<Item>" },
    { "/File/Open/Open S_QLite Database",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,                85, "<Item>" },
    { "/File/sep1",         NULL,         NULL,           0,                    "<Separator>" },
    { "/File/_Quit",     "<control>Q",    (GtkItemFactoryCallback)hyphy_menu_item_callback,             15, "<StockItem>",GTK_STOCK_QUIT },
    { "/Edit/sep_pref",     NULL,    NULL,              0, "<Separator>" },
    { "/Edit/Pre_ferences",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,              27, "<Item>" },
    { "/_Analysis",         NULL,         NULL,           0,                    "<Branch>" } ,
    { "/Analysis/_Cancel execution",     "<control>.",    (GtkItemFactoryCallback)hyphy_menu_item_callback,             21, "<Item>"},
    { "/Analysis/_Suspend execution",     "<control>;",    (GtkItemFactoryCallback)hyphy_menu_item_callback,                22, "<Item>"},
    { "/Analysis/sep1",         NULL,         NULL,           0,                    "<Separator>" },
    { "/Analysis/_View messages.log",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,                23, "<Item>"},
    { "/Analysis/sep2",         NULL,         NULL,           0,                    "<Separator>" },
    { "/Analysis/_Standard analyses...",     "<control>E",    (GtkItemFactoryCallback)hyphy_menu_item_callback,             24, "<Item>"},
    { "/Analysis/Results",          NULL,         NULL,           0,                    "<Branch>" } ,
    { "/Analysis/_Rerun last analysis",     "<control>R",    (GtkItemFactoryCallback)hyphy_menu_item_callback,              30, "<Item>"},
    { "/Analysis/sep3",         NULL,         NULL,           0,                    "<Separator>" },
    { "/Analysis/Co_mpute expression",     "<control>K",    (GtkItemFactoryCallback)hyphy_menu_item_callback,               60, "<Item>"},
    { "/Analysis/_Execute selection",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,                61, "<Item>"}
};

//__________________________________________________________________

void _HYConsoleWindow::_SetMenuBar(void)
{
    if (menu_items && !gtk_item_factory_get_item(menu_items,"<HY_WINDOW>/Analysis")) {
        gtk_item_factory_delete_item (menu_items, "<HY_WINDOW>/File/Open");
        gtk_item_factory_delete_item (menu_items, "<HY_WINDOW>/File/Close");
        gtk_item_factory_delete_item (menu_items, "<HY_WINDOW>/File/Switch to console");
        gtk_item_factory_create_items (menu_items,  sizeof (hyphy_console_window_menu_file) / sizeof (hyphy_console_window_menu_file[0]), hyphy_console_window_menu_file, this);

        GtkMenu * fileMenu = GTK_MENU (gtk_item_factory_get_widget (menu_items, "<HY_WINDOW>/File"));
        gtk_menu_reorder_child (fileMenu, gtk_item_factory_get_item (menu_items, "<HY_WINDOW>/File/New"), 0);
        gtk_menu_reorder_child (fileMenu, gtk_item_factory_get_item (menu_items, "<HY_WINDOW>/File/Open"), 1);
        gtk_menu_reorder_child (fileMenu, gtk_item_factory_get_item (menu_items, "<HY_WINDOW>/File/sep1"), 2);

        gtk_widget_set_sensitive(gtk_item_factory_get_widget (menu_items,"<HY_WINDOW>/Edit/Clear"),  true);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget (menu_items,"<HY_WINDOW>/Edit/Select All"),  true);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget_by_action(menu_items,21),  false);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget_by_action(menu_items,22),  false);
        gtk_widget_set_sensitive(gtk_item_factory_get_item (menu_items,"<HY_WINDOW>/Analysis/Results"),  false);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget_by_action(menu_items,30),  false);

        for (long counter=0; counter<availablePostProcessors.countitems(); counter++) {
            _String * postItem = (_String*)(*(_List*)availablePostProcessors(counter))(0);
            if (!postItem->Equal (&menuSeparator)) {
                postItem = (_String*)(*(_List*)availablePostProcessors(counter))(1);
                _String tryFileName (*postItem);
                FILE    * tryFile   = doFileOpen (tryFileName.getStr(), "r");
                if (tryFile) {
                    fclose(tryFile);
                } else {
                    availablePostProcessors.Delete(counter);
                    counter--;
                }
            }
        }
        if (availablePostProcessors.lLength) {
            long sepCount = 0;
            for (long counter=0; counter<availablePostProcessors.countitems(); counter++) {
                _String * postItem = (_String*)(*(_List*)availablePostProcessors(counter))(0),
                          chopped  = _String("/Analysis/Results/")&*postItem;

                if (postItem->Equal (&menuSeparator)) {
                    GtkItemFactoryEntry aProcEntry = {NULL,NULL,NULL,0,"<Separator>"};
                    chopped = chopped & sepCount++;
                    aProcEntry.path = chopped.sData;
                    gtk_item_factory_create_items (menu_items,  1, &aProcEntry, this);
                } else {
                    GtkItemFactoryEntry aProcEntry = {NULL,NULL,(GtkItemFactoryCallback)hyphy_menu_item_callback,1000+counter,"<Item>"};
                    aProcEntry.path = chopped.sData;
                    gtk_item_factory_create_items (menu_items,  1, &aProcEntry, this);
                }
            }

        }
    }
}


//__________________________________________________________________

void _HYConsoleWindow::_UpdateEditMenu (void)
{
    /*
    _HYTextBox* txb = (_HYTextBox*)GetObject (0);

    HWND        te = txb->te;

    CHARRANGE   cr;

    SendMessage (te,EM_EXGETSEL,0,(LPARAM)&cr);

    bool        haveSelection = cr.cpMax-cr.cpMin,
                canPaste      = SendMessage (te,EM_CANPASTE,CF_TEXT,0),
                canUndo       = SendMessage (te,EM_GETUNDONAME,0,0),
                canRedo       = SendMessage (te,EM_CANREDO,0,0);

    if ((((bool)editOptions&HY_CONSOLE_CAN_COPY)!=haveSelection)||
        (((bool)editOptions&HY_CONSOLE_CAN_PASTE)!=canPaste)||
        (((bool)editOptions&HY_CONSOLE_CAN_UNDO)!=canUndo)||
        (((bool)editOptions&HY_CONSOLE_CAN_UNDO)!=canRedo))
    {
        HMENU            windowMenu = GetMenu (theWindow),
                         editMenu   = GetSubMenu(windowMenu,1);

        editOptions = 0;

        if (haveSelection)
            editOptions |= HY_CONSOLE_CAN_COPY;

        EnableMenuItem (editMenu,3,MF_BYPOSITION|(haveSelection?MF_ENABLED:MF_GRAYED));
        EnableMenuItem (editMenu,4,MF_BYPOSITION|(haveSelection?MF_ENABLED:MF_GRAYED));

        if (canPaste)
            editOptions |= HY_CONSOLE_CAN_PASTE;

        EnableMenuItem (editMenu,5,MF_BYPOSITION|(canPaste?MF_ENABLED:MF_GRAYED));

        if (canUndo)
            editOptions |= HY_CONSOLE_CAN_UNDO;

        EnableMenuItem (editMenu,0,MF_BYPOSITION|(canUndo?MF_ENABLED:MF_GRAYED));

        if (canRedo)
            editOptions |= HY_CONSOLE_CAN_REDO;

        EnableMenuItem (editMenu,1,MF_BYPOSITION|(canRedo?MF_ENABLED:MF_GRAYED));
    }*/
}

//__________________________________________________________________

void _HYConsoleWindow::_UnsetMenuBar(void)
{
    _HYWindow::_UnsetMenuBar();
}

//__________________________________________________________________

bool        _HYConsoleWindow::_ProcessMenuSelection (long msel)
{
    switch (msel) {
    case HY_WINDOW_MENU_ID_FILE+1:
        SaveConsole ();
        return true;

    case 14:
        _DoPrint ();
        return true;

    case 16:
        ((_HYTextBox*)GetObject(0))->_DoUndo(true);
        return true;

    case 17:
        ((_HYTextBox*)GetObject(0))->_DoCut(true);
        return true;

    case 18:
        ((_HYTextBox*)GetObject(0))->_DoCopy(true);
        return true;

    case HY_WINDOW_MENU_ID_EDIT+5:
        ((_HYTextBox*)GetObject(0))->_DoSelectAll(true);
        return true;

    case 21:
        SetStatusLine ("Canceling");
        terminateExecution = true;
        return true;

    case 22: {
        GtkWidget * suspendItem = gtk_item_factory_get_widget_by_action(menu_items,22);
        if (!isSuspended) {
            isSuspended = true;
            SetStatusLine (empty,empty,empty,-1,HY_SL_SUSPEND);
            gtk_label_set_text (GTK_LABEL (gtk_bin_get_child(GTK_BIN(suspendItem))), "Resume");
            updateTimer = false;
            while (isSuspended) {
                gtk_main_iteration_do(true);
            }
        } else {
            isSuspended = false;
            SetStatusLine (empty,empty,empty,-1,HY_SL_RESUME);
            gtk_label_set_text (GTK_LABEL (gtk_bin_get_child(GTK_BIN(suspendItem))), "Suspend execution");
            timerStart += clock()-lastTimer;
            updateTimer = true;
        }
        return true;
    }
    case 23:
        ShowMessagesLog();
        return true;

    case 24:
        RunStandardAnalyses();
        return true;

    case 25:
        displayAbout(false);
        return true;

    case HY_WINDOW_MENU_ID_EDIT+4:
        ((_HYTextBox*)GetObject(0))->_DoClear (true,true);
        return true;

    case 27:
        HandlePreferences (globalPreferencesList,"HYPHY Preferences");
        return true;

    case 29:
        //WinExec ("hh HYPHY HELP.chm",SW_SHOWNORMAL);
        return true;

    case 30:
        if (OpenBatchFile (false)) {
            ExecuteBatchFile();
            PopFilePath();
        }
        return true;

    case 31:
        ((_HYTextBox*)GetObject(0))->_DoRedo(true);
        return true;

    case 15:
        postWindowCloseEvent (GetID());
        gtk_main_quit ();
        return true;

    case 60: // expression calculator
        if (calculatorMode) {
            _HYTextBox         *ib = (_HYTextBox*)hyphyConsoleWindow->GetObject(1);
            ib->SetText ("exit");
            hyphyConsoleWindow->ProcessEvent (generateTextEditChangeEvent(ib->GetID(),2));
            calculatorMode         = false;
            //ib->SetText (empty);
        } else {
            calculatorMode = true;
            while(calculatorMode&&ExpressionCalculator()) {}
            calculatorMode = false;
        }
        return true;

    case 61: { // execute selection
        ExecuteSelection();
        return true;
    }

    case 70: // New Tree
        NewTreeWindow(-1);
        return true;

    case 71: // New Model
        NewModel(nil);
        return true;

    case 72: // New Chart
        NewChartWindow();
        return true;

    case 73: // New Genetic Code
        NewGeneticCodeTable(0);
        return true;

    case 74: // New Database
        NewDatabaseFile(0);
        return true;

    case 80: // Open Batch File
        if (OpenBatchFile()) {
            ExecuteBatchFile ();
        }
        return true;

    case 81: // Open Data File
        OpenDataFile();
        return true;

    case 82: // Open Tree
        OpenTreeFile();
        return true;

    case 83: // Open Text
        OpenTextFile();
        return true;

    case 84: // Open Table
        OpenTable ();
        return true;

    case 85: // Open SQLite database
        OpenDatabaseFile (nil);
        return true;

    case HY_WINDOW_MENU_ID_FILE-2:
        ShowObjectInspector ();
        return true;

    default: {
        msel -= 1000;
        if (msel<availablePostProcessors.lLength) {
            ExecuteAPostProcessor (*(_String*)(*(_List*)availablePostProcessors(msel))(1));
            return 0;
        }

        msel-=1000;
        if (msel<(long)recentPaths.lLength) {
            if (msel == -1) {
                for (long mi=0; mi<recentFiles.lLength; mi++) {
                    GtkWidget * recFile = gtk_item_factory_get_widget_by_action(hyphyConsoleWindow->menu_items,2000+mi);
                    if (recFile) {
                        gtk_item_factory_delete_item(hyphyConsoleWindow->menu_items,gtk_item_factory_path_from_widget(recFile));
                    }
                }
                recentPaths.Clear();
                recentFiles.Clear();
            } else {
                if (argFileName) {
                    *argFileName = *(_String*)recentPaths(msel);
                } else {
                    argFileName = new _String (*(_String*)recentPaths(msel));
                }
                if (OpenBatchFile(false)) {
                    ExecuteBatchFile ();
                }
            }
            return true;
        }
        return true;
    }
    }

    return _HYTWindow::_ProcessMenuSelection(msel);
}


//__________________________________________________________________

bool _HYConsoleWindow::_ProcessOSEvent (Ptr vEvent)
{
    /*_HYWindowsUIMessage*  theEvent = (_HYWindowsUIMessage*)vEvent;
    if (theEvent->iMsg == WM_SYSCOMMAND)
    {
        if (theEvent->wParam == SC_CLOSE)
        {
            hyphyExiting = true;
            return true;
        }
    }
    else
    {
        if (theEvent->iMsg == UPDATE_TIMER)
        {
            UpdateTimer ();
            return true;
        }
    }*/
    return _HYTWindow::_ProcessOSEvent (vEvent);
}

//__________________________________________________________________

bool _HYConsoleWindow::_Close (Ptr )
{
    _String closingNow ("Please confirm that you want to close HyPhy, terminating all unfinished tasks.");
    if (ProceedPrompt (closingNow, (Ptr)this)) {
        gtk_main_quit();
    }
    return false;
}

//__________________________________________________________________
void    _HYConsoleWindow::_DoPrint          (void)
{
    _String NIY ("Printing has not yet been implemented\n");
    ProblemReport (NIY);
}

//__________________________________________________________________

void _HYConsoleWindow::_PaintStatusBar(Ptr,bool force)
{
    if (GTK_WIDGET_MAPPED (theWindow)) {
        _Parameter      vL;
        checkParameter (VerbosityLevelString, vL, 0.0);

        if (vL<-0.5 && !force) {
            clock_t curMeasure = clock();
            _Parameter diff = 1.0/CLOCKS_PER_SEC*(curMeasure-lastMeasure);
            if (diff < -vL) {
                return;
            }
            lastMeasure = curMeasure;
        }

        if (!stripedFillGC) {
            SetUpStatusBarStuff (theWindow);
        }

        GdkRectangle wRC  = {0,0,theWindow->allocation.width,HY_SCROLLER_WIDTH},
                     w2RC;

        if (!_hyConsoleWindowGC) {
            _hyConsoleWindowGC   = gdk_gc_new     (theWindow->window);
            gdk_gc_set_tile (_hyConsoleWindowGC, stripedFill);
            gdk_gc_set_line_attributes(_hyConsoleWindowGC,1,GDK_LINE_SOLID,GDK_CAP_NOT_LAST,GDK_JOIN_MITER);
            _hyConsoleLayout = pango_layout_new (screenPContext);
            pango_layout_set_width (_hyConsoleLayout, -1);
            statusBarBold = pango_font_description_new ();
            statusBarNormal = pango_font_description_new ();
            _HYFont         consoleFont = {_HY_MONO_FONT,9,HY_FONT_PLAIN};
            HYFont2PangoFontDesc (consoleFont,statusBarNormal);
            consoleFont.style = HY_FONT_BOLD;
            HYFont2PangoFontDesc (consoleFont,statusBarBold);
        }

        GdkPixmap * offBitmap        = gdk_pixmap_new (theWindow->window, wRC.width, wRC.height,-1);

        if (offBitmap) {
            gdk_gc_set_fill         (_hyConsoleWindowGC,GDK_TILED);
            gdk_draw_rectangle      (offBitmap,_hyConsoleWindowGC,true,-1,-1,wRC.width+2,wRC.height+2);
            gdk_gc_set_fill         (_hyConsoleWindowGC,GDK_SOLID);
            gdk_gc_set_foreground   (_hyConsoleWindowGC,&_BLACKBRUSH_);
            gdk_draw_line(offBitmap,_hyConsoleWindowGC,0,0,wRC.width,0);

            pango_layout_set_font_description(_hyConsoleLayout,statusBarNormal);
            pango_layout_set_text(_hyConsoleLayout,fileName.getStr(),fileName.sLength);
            gdk_draw_layout(offBitmap,_hyConsoleWindowGC,33,wRC.height-13,_hyConsoleLayout);

            if (inputStatus == 1) {
                pango_layout_set_text(_hyConsoleLayout,cInput.getStr(),cInput.sLength);
            } else {
                pango_layout_set_text(_hyConsoleLayout,action.getStr(),action.sLength);
            }

            gdk_draw_layout(offBitmap,_hyConsoleWindowGC,193,wRC.height-13,_hyConsoleLayout);

            gdk_gc_set_foreground   (_hyConsoleWindowGC,&_DARKGREYBRUSH_);

            gdk_draw_rectangle      (offBitmap,_hyConsoleWindowGC,true,0,1,30,wRC.height-1);
            gdk_draw_rectangle      (offBitmap,_hyConsoleWindowGC,true,150,1,40,wRC.height-1);
            gdk_draw_rectangle      (offBitmap,_hyConsoleWindowGC,true,wRC.width-55,1,55,wRC.height-1);

            gdk_gc_set_foreground   (_hyConsoleWindowGC,&_WHITEBRUSH_);
            pango_layout_set_font_description(_hyConsoleLayout,statusBarBold);

            pango_layout_set_text(_hyConsoleLayout,cState.getStr(),cState.sLength);
            gdk_draw_layout(offBitmap,_hyConsoleWindowGC,3,wRC.height-13,_hyConsoleLayout);
            pango_layout_set_text(_hyConsoleLayout,cTask.getStr(),cTask.sLength);
            gdk_draw_layout(offBitmap,_hyConsoleWindowGC,151,wRC.height-13,_hyConsoleLayout);

            pango_layout_set_font_description(_hyConsoleLayout,statusBarNormal);

            pango_layout_set_text(_hyConsoleLayout,timer.getStr(),timer.sLength);
            gdk_draw_layout(offBitmap,_hyConsoleWindowGC,wRC.width-53,wRC.height-13,_hyConsoleLayout);

            if (percentDone>0 || percentDone == -HY_SL_DONE) {
                GdkColor blackBrush = HYColorToGDKColor((_HYColor) {
                    80,80,80
                }),
                orangeBrush = HYColorToGDKColor((_HYColor) {
                    255,153,102
                });

                gdk_gc_set_foreground   (_hyConsoleWindowGC,&orangeBrush);
                gdk_draw_rectangle      (offBitmap,_hyConsoleWindowGC,true,wRC.width-135,wRC.height-14,(percentDone>=0?percentDone:100.)*0.75,12);

                gdk_gc_set_foreground   (_hyConsoleWindowGC,&blackBrush);
                gdk_draw_rectangle      (offBitmap,_hyConsoleWindowGC,false,wRC.width-135,wRC.height-14,75,12);

                //gdk_gc_set_foreground (_hyConsoleWindowGC,&_WHITEBRUSH_);
                _String pLine;
                if (percentDone>=0) {
                    pLine = _String(percentDone)&"%";
                } else {
                    pLine = "DONE";
                }

                pango_layout_set_text(_hyConsoleLayout,pLine.getStr(),pLine.sLength);
                gdk_draw_layout(offBitmap,_hyConsoleWindowGC,wRC.width-107,wRC.height-13,_hyConsoleLayout);
            }
        }

        gdk_draw_drawable (theWindow->window, _hyConsoleWindowGC, offBitmap, 0, 0, theWindow->allocation.x, theWindow->allocation.y+theWindow->allocation.height-wRC.height, -1, -1);
        g_object_unref (offBitmap);
    }
    yieldCPUTime();
}

