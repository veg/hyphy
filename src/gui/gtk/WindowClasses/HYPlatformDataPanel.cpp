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
    GTK Portions of the data panel class

    Sergei L. Kosakovsky Pond, Spring 2005.
*/

#include "HYTreePanel.h"
#include "HYDataPanel.h"
#include "HYUtils.h"

#include "likefunc.h"
#include <gdk/gdkkeysyms.h>
#include <gdk/gdkx.h>

#define     HY_DATA_WIN32_MENU_BASE    8000
#define     HY_DATALF_WIN32_MENU_BASE  8500

//__________________________________________________________________

static GtkItemFactoryEntry hyphy_data_window_menu[] = {
    { "/File/Save _as",     "<control>F",    (GtkItemFactoryCallback)hyphy_menu_item_callback,              HY_WINDOW_MENU_ID_FILE+3, "<Item>" },
    { "/Edit/_Find",     "<control>F",    (GtkItemFactoryCallback)hyphy_menu_item_callback,             HY_WINDOW_MENU_ID_EDIT+6, "<Item>" },
    { "/Edit/_Search and Replace",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,               HY_WINDOW_MENU_ID_EDIT+7, "<Item>" },
    { "/_Data",         NULL,         NULL,           0,                    "<Branch>" },
    { "/Data/_Partition->Selection", "<control>1", (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE, "<Item>"},
    { "/Data/_Selection->Partition", "<control>2", (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+1, "<Item>"},
    { "/Data/_Invert selection", "<control>3", (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+2, "<Item>"},
    { "/Data/sep1",         NULL,         NULL,           0,        "<Separator>" },
    { "/Data/_Block width",         NULL,         NULL,           0,        "<Branch>" },
    { "/Data/Block width/_9", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+50, "<CheckItem>"},
    { "/Data/Block width/_10", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+51, "<CheckItem>"},
    { "/Data/_Repeat characters",           NULL,         NULL,           0,        "<Branch>" },
    { "/Data/Repeat characters/Display _actual character", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+60, "<CheckItem>"},
    { "/Data/Repeat characters/_Display '.'", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+61, "<CheckItem>"},
    { "/Data/_Name Display",            NULL,         NULL,           0,        "<Branch>" },
    { "/Data/Name Display/_None", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+70, "<CheckItem>"},
    { "/Data/Name Display/_First 10 characters", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+71, "<CheckItem>"},
    { "/Data/Name Display/F_ull names", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+72, "<CheckItem>"},
    { "/Data/Name Display/sep1",            NULL,         NULL,           0,        "<Separator>" },
    { "/Data/Name Display/_Alphabetize names", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+73, "<Item>"},
    { "/Data/Name Display/_Revert to file order", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+74, "<Item>"},
    { "/Data/Name Display/_Clean up sequence names", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+75, "<Item>"},
    { "/Data/_Omitted sequences",           NULL,         NULL,           0,        "<Branch>" },
    { "/Data/Omitted sequences/_Restore all", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+10000, "<Item>"},
    { "/Data/_Omitted sequences/sep1",          NULL,         NULL,           0,        "<Separator>" },
    { "/Data/_Additional information",          NULL,         NULL,           0,        "<Branch>" },
    { "/Data/Additional information/_Consensus sequence", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+80, "<CheckItem>"},
    { "/Data/Additional information/_Rate class", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+81, "<CheckItem>"},
    { "/Data/Additional information/_Aminoacid translation", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+82, "<CheckItem>"},
    { "/Data/Additional information/Re_ference sequence", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+83, "<CheckItem>"},
    { "/Data/sep2",         NULL,         NULL,           0,        "<Separator>" },
    { "/Data/Part_ition properties", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+3, "<Item>"},
    { "/Data/Inp_ut partition", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+4, "<Item>"},
    { "/Data/sep3",         NULL,         NULL,           0,        "<Separator>" },
    { "/Data/_Simulation",          NULL,         NULL,           0,        "<Branch>" },
    { "/Data/Simulation/_Simulate 1",NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+90, "<Item>"},
    { "/Data/Simulation/Simulate 1 to _file",NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+91, "<Item>"},
    { "/Data/Simulation/Simulate _many",NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+92, "<Item>"},
    { "/Data/An_cestors", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+5, "<Item>"},
    { "/Data/sep4",         NULL,         NULL,           0,        "<Separator>" },
    { "/Data/_Font options", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATA_WIN32_MENU_BASE+6, "<Item>"},
    { "/Data/sep5",         NULL,         NULL,           0,        "<Separator>" },
    { "/Data/_Data Processing",         NULL,         NULL,           0,                    "<Branch>" },
    { "/_Likelihood",           NULL,         NULL,           0,                    "<Branch>" },
    { "/Likelihood/Build", "<control>L", (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATALF_WIN32_MENU_BASE, "<Item>"},
    { "/Likelihood/_Display",           NULL,         NULL,           0,                    "<Branch>" },
    { "/Likelihood/Display/Log-Lkhd _only", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATALF_WIN32_MENU_BASE+50, "<CheckItem>"},
    { "/Likelihood/Display/Log-Lkhd with &parameter values", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATALF_WIN32_MENU_BASE+51, "<CheckItem>"},
    { "/Likelihood/Display/Log-Lkhd with &concise trees", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATALF_WIN32_MENU_BASE+52, "<CheckItem>"},
    { "/Likelihood/Display/_Parameter Listing", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATALF_WIN32_MENU_BASE+53, "<CheckItem>"},
    { "/Likelihood/Display/Log-Lkhd with complete &trees", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATALF_WIN32_MENU_BASE+54, "<CheckItem>"},
    { "/Likelihood/_Optimize", "<control>T", (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATALF_WIN32_MENU_BASE+1, "<Item>"},
    { "/Likelihood/sep4",           NULL,         NULL,           0,        "<Separator>" },
    { "/Likelihood/_Show parameters", "<control>H", (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATALF_WIN32_MENU_BASE+2, "<Item>"},
    { "/Likelihood/_General bootstrap", "<control>B", (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATALF_WIN32_MENU_BASE+3, "<Item>"}
};

static GtkItemFactoryEntry hyphy_data_window_menu2[] = {
    { "/Likelihood/_Inference",         NULL,         NULL,           0,                    "<Branch>" },
    { "/Likelihood/Inference/_Infer topology", "<control>L", (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATALF_WIN32_MENU_BASE+20, "<Item>"},
    { "/Likelihood/Inference/Infer topology with _constraints", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_DATALF_WIN32_MENU_BASE+21, "<Item>"}
};

//__________________________________________________________________

void _HYDataPanel::_SetMenuBar(void)
{
    _HYWindow::_SetMenuBar();

    if (menu_items && !gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Data")) {
        _HYSequencePane* sp = (_HYSequencePane*)GetObject (0);

        gtk_item_factory_create_items (menu_items,  sizeof (hyphy_data_window_menu) / sizeof (hyphy_data_window_menu[0]), hyphy_data_window_menu, this);
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(menu_items,HY_DATA_WIN32_MENU_BASE),false);
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(menu_items,HY_DATA_WIN32_MENU_BASE+1),false);
        gtk_widget_set_sensitive(gtk_item_factory_get_item(menu_items,"<HY_WINDOW>/Data/Omitted sequences"),false);

        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(menu_items,HY_DATA_WIN32_MENU_BASE+5),false);
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(menu_items,HY_DATA_WIN32_MENU_BASE+81),false);
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(menu_items,HY_DATA_WIN32_MENU_BASE+82),false);

        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(menu_items,HY_DATALF_WIN32_MENU_BASE+1),false);
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(menu_items,HY_DATALF_WIN32_MENU_BASE+2),false);
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(menu_items,HY_DATALF_WIN32_MENU_BASE+3),false);

        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(menu_items,HY_WINDOW_MENU_ID_EDIT+1),true);
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(menu_items,HY_WINDOW_MENU_ID_EDIT+5),true);

        gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(gtk_item_factory_get_item_by_action(menu_items,HY_DATA_WIN32_MENU_BASE+50+(sp->blockWidth==10))),true);
        gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(gtk_item_factory_get_item_by_action(menu_items,HY_DATA_WIN32_MENU_BASE+60+sp->showDots)),true);

        if (omittedSeqs.lLength) {
            _OmitSelectedSpecies(omittedSeqs);
        }

        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(menu_items,HY_DATA_WIN32_MENU_BASE+82),dataType&HY_DATAPANEL_NUCDATA);

        if (sp->nameDisplayFlags&HY_SEQUENCE_PANE_NAMES_ALL) {
            gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(gtk_item_factory_get_item_by_action(menu_items,HY_DATA_WIN32_MENU_BASE+72)),true);
        } else if (sp->nameDisplayFlags&HY_SEQUENCE_PANE_NAMES_SHORT) {
            gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(gtk_item_factory_get_item_by_action(menu_items,HY_DATA_WIN32_MENU_BASE+71)),true);
        } else {
            gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(gtk_item_factory_get_item_by_action(menu_items,HY_DATA_WIN32_MENU_BASE+70)),true);
        }

        GtkMenu *fileMenu = GTK_MENU(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/File"));
        gtk_menu_reorder_child(fileMenu,gtk_item_factory_get_widget_by_action(menu_items,HY_WINDOW_MENU_ID_FILE+3),2);

        for (long k=0; k<dataPanelProcessors.lLength; k++) {
            _String *thisItem = (_String*)dataPanelProcessors (k),
                     chopped = thisItem->Cut (thisItem->FindBackwards ('/',0,-1)+1,-1),
                     type = "<Item>";

            GtkItemFactoryEntry aProcEntry = {NULL,NULL,(GtkItemFactoryCallback)hyphy_menu_item_callback,HY_DATA_WIN32_MENU_BASE+100+k,type.sData};
            chopped = _String("/Data/Data Processing/")&chopped;
            aProcEntry.path = chopped.sData;

            gtk_item_factory_create_items (menu_items,  1, &aProcEntry, this);
        }

        _UpdateLFMenu();
        _VerifyInferMenu    ();
    }
}
//__________________________________________________________________


bool        _HYDataPanel::_ProcessMenuSelection (long msel)
{
    if (_HYWindow::_ProcessMenuSelection(msel)) {
        return true;
    }

    _HYSequencePane* sp = (_HYSequencePane*)GetObject (0);
    _HYSequencePane* sp2 =(_HYSequencePane*)GetObject (4);
    _String     prompt;
    bool        done = false;

    switch (msel) {
    case HY_WINDOW_MENU_ID_FILE+2: {
        _PrintData();
        done = true;
        break;
    }

    case HY_WINDOW_MENU_ID_EDIT+1: { // Copy
        _CopySelectionToClipboard   ();
        done = true;
        break;
    }

    case HY_WINDOW_MENU_ID_EDIT+5: { // Select All
        sp->SelectAll(true);
        done = true;
        break;
    }

    case HY_WINDOW_MENU_ID_EDIT+6: { // Find Function
        FindFunction();
        done = true;
        break;
    }

    case HY_WINDOW_MENU_ID_EDIT+7: { // Find Function
        HandleSearchAndReplace();
        done = true;
        break;
    }

    case HY_DATA_WIN32_MENU_BASE: { // Partition->Selection
        SelectPartition();
        done = true;
        break;
    }

    case HY_DATA_WIN32_MENU_BASE+1: { // Selection->Partition
        if (sp->selection.lLength) {
            CreatePartition (sp->selection,1,true);
        }
        done = true;
        break;
    }

    case HY_DATA_WIN32_MENU_BASE+2: { // Invert Selection
        InvertSelection();
        done = true;
        break;
    }

    case HY_DATA_WIN32_MENU_BASE+3: { // Parition props
        PartitionPropsMenu ();
        done = true;
        break;
    }

    case HY_DATA_WIN32_MENU_BASE+4: { // Input part
        InputPartitionString ();
        done = true;
        break;
    }

    case HY_DATA_WIN32_MENU_BASE+5: { // Ancestors
        SimulateDataSet (0,true);
        done = true;
        break;
    }

    case HY_DATA_WIN32_MENU_BASE+6: { // Ancestors
        HandleFontChange();
        done = true;
        break;
    }

    case HY_DATALF_WIN32_MENU_BASE: { // Build LF
        if (gtk_item_factory_get_item (menu_items,"<HY_WINDOW>/Likelihood/Inference")) {
            InferTopologies ();
        } else {
            BuildLikelihoodFunction();
        }
        done = true;
        break;
    }

    case HY_DATALF_WIN32_MENU_BASE+1: { // Optimize LF
        OptimizeLikelihoodFunction();
        done = true;
        break;
    }

    case HY_DATALF_WIN32_MENU_BASE+2: { // Show Parameters
        DisplayParameterTable ();
        done = true;
        break;
    }

    case HY_DATALF_WIN32_MENU_BASE+3: { // General bootstrap
        //Uncomment when rdy
        OpenGeneralBSWindow ();
        done = true;
        break;
    }

    case HY_DATA_WIN32_MENU_BASE+50: // Block size
    case HY_DATA_WIN32_MENU_BASE+51: {
        GtkCheckMenuItem *checkItem = GTK_CHECK_MENU_ITEM(gtk_item_factory_get_item_by_action(menu_items,msel));

        if (gtk_check_menu_item_get_active(checkItem)) {
            long             newBlockSize;
            bool             is9 = (msel==HY_DATA_WIN32_MENU_BASE+50);

            if (is9) {
                newBlockSize = 9;
            } else {
                newBlockSize = 10;
            }

            if (sp->blockWidth!=newBlockSize) {
                sp->blockWidth = newBlockSize;
                sp2->blockWidth = newBlockSize;
                sp->BuildPane();
                sp->_MarkForUpdate();
                sp2->BuildPane();
                sp2->_MarkForUpdate();
                gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(gtk_item_factory_get_item_by_action(menu_items,is9?msel+1:msel-1)),false);
            }
        }
        done = true;
        break;
    }

    case HY_DATA_WIN32_MENU_BASE+60: // Repeating character
    case HY_DATA_WIN32_MENU_BASE+61: {
        GtkCheckMenuItem *checkItem = GTK_CHECK_MENU_ITEM(gtk_item_factory_get_item_by_action(menu_items,msel));

        if (gtk_check_menu_item_get_active(checkItem)) {
            bool       newDisplay;
            if (msel==HY_DATA_WIN32_MENU_BASE+60) {
                newDisplay = false;
            } else {
                newDisplay = true;
            }

            if (sp->showDots!=newDisplay) {
                sp->showDots = newDisplay;
                sp->BuildPane();
                sp->_MarkForUpdate();
                gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(gtk_item_factory_get_item_by_action(menu_items,newDisplay?msel-1:msel+1)),false);
            }
        }
        done = true;
        break;
    }

    case HY_DATA_WIN32_MENU_BASE+70: // Sequence names
    case HY_DATA_WIN32_MENU_BASE+71:
    case HY_DATA_WIN32_MENU_BASE+72:
    case HY_DATA_WIN32_MENU_BASE+73:
    case HY_DATA_WIN32_MENU_BASE+74:
    case HY_DATA_WIN32_MENU_BASE+75: {
        if (msel<=HY_DATA_WIN32_MENU_BASE+72) {
            GtkCheckMenuItem *checkItem = GTK_CHECK_MENU_ITEM(gtk_item_factory_get_item_by_action(menu_items,msel));
            if (gtk_check_menu_item_get_active(checkItem)) {
                unsigned char newDisplay;
                switch (msel) {
                case HY_DATA_WIN32_MENU_BASE+70:
                    newDisplay = HY_SEQUENCE_PANE_NAMES_NONE;
                    break;
                case HY_DATA_WIN32_MENU_BASE+71:
                    newDisplay = HY_SEQUENCE_PANE_NAMES_SHORT;
                    break;
                case HY_DATA_WIN32_MENU_BASE+72:
                    newDisplay = HY_SEQUENCE_PANE_NAMES_ALL;

                }

                if ((sp->nameDisplayFlags&HY_SEQUENCE_PANE_NAMES_MASK)!=newDisplay) {
                    sp->SetNameDisplayMode(newDisplay,true);
                    sp2->SetNameDisplayMode(newDisplay,true);
                    BuildThermometer();
                    BuildMarksPane();
                    for (long k = HY_DATA_WIN32_MENU_BASE+70; k<HY_DATA_WIN32_MENU_BASE+73; k++)
                        if (k!=msel) {
                            gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(gtk_item_factory_get_item_by_action(menu_items,k)),false);
                        }

                }
            }
        } else {
            if (msel==HY_DATA_WIN32_MENU_BASE+73) {
                sp->AlphabetizeSpecies();
            } else if (msel==HY_DATA_WIN32_MENU_BASE+74) {
                sp->RevertFileOrder();
            } else {
                sp->CleanUpSequenceNames();
            }
        }
        done = true;
        break;
    }

    case HY_DATA_WIN32_MENU_BASE+80: // status lines
    case HY_DATA_WIN32_MENU_BASE+81: // status lines
    case HY_DATA_WIN32_MENU_BASE+82: // status lines
    case HY_DATA_WIN32_MENU_BASE+83: { // status lines
        AdjustStatusLine (msel-HY_DATA_WIN32_MENU_BASE-80);
        done = true;
        break;
    }

    case HY_DATALF_WIN32_MENU_BASE+50: // likelihood display
    case HY_DATALF_WIN32_MENU_BASE+51: // likelihood display
    case HY_DATALF_WIN32_MENU_BASE+52: // likelihood display
    case HY_DATALF_WIN32_MENU_BASE+53: // likelihood display
    case HY_DATALF_WIN32_MENU_BASE+54: { // likelihood display
        ComputeLikelihoodFunction (msel-HY_DATALF_WIN32_MENU_BASE-50);
        done = true;
        break;
    }

    case HY_DATA_WIN32_MENU_BASE+90: // simulate data set
    case HY_DATA_WIN32_MENU_BASE+91: // simulate data set
    case HY_DATA_WIN32_MENU_BASE+92: { // simulate data set
        SimulateDataSet (msel-HY_DATA_WIN32_MENU_BASE-90);
        done = true;
        break;
    }

    case HY_WINDOW_MENU_ID_FILE+1: // save/save as
    case HY_WINDOW_MENU_ID_FILE+3: { // save/save as
        SaveDataPanel (msel==HY_WINDOW_MENU_ID_FILE+1);
        done = true;
        break;
    }

    case HY_DATALF_WIN32_MENU_BASE+20: // infer
    case HY_DATALF_WIN32_MENU_BASE+21: { // infer
        InferTopologies (msel==HY_DATALF_WIN32_MENU_BASE+21);
        _VerifyInferMenu ();
        done = true;
        break;
    }

    default: {
        if (msel>=HY_DATA_WIN32_MENU_BASE+10000) {
            RestoreOmittedSequence(msel-HY_DATA_WIN32_MENU_BASE-10003);
            done = true;
            break;
        } else if (msel>=HY_DATA_WIN32_MENU_BASE+100) {
            ExecuteProcessor(msel-HY_DATA_WIN32_MENU_BASE-100);
            done = true;
            break;
        }
    }
    }
    return done;
}

//__________________________________________________________________

void        _HYDataPanel::_PaintThermRect(bool update)
{
    navRect = ComputeNavRect();
    _HYCanvas* theCanvas = (_HYCanvas*)GetObject (1);
    GdkRectangle       r;

    r.x      = navRect.left+theCanvas->rel.left+thermRect.left+1+windowContent->allocation.x+1;
    r.width  = navRect.right-navRect.left-3;
    r.y      = navRect.top+theCanvas->rel.top+thermRect.top+1+windowContent->allocation.y+1;
    r.height = navRect.bottom-navRect.top-3;

    GdkColor   newBr  = HYColorToGDKColor((_HYColor) {
        255,151,51
    });
    GdkGC      *gc    = gdk_gc_new (theWindow->window);

    gdk_gc_set_foreground(gc,&newBr);
    gdk_gc_set_line_attributes(gc,2,GDK_LINE_SOLID,GDK_CAP_NOT_LAST,GDK_JOIN_MITER);
    gdk_draw_rectangle(theWindow->window,gc,false,r.x,r.y,r.width,r.height);

    if (update) {
        _HYRect rect = {componentT.lData[1],componentL.lData[1],componentB.lData[1],componentR.lData[1],0};

        r.x--;
        r.y--;
        r.height+=2;
        r.width +=2;

        GdkRegion * rg1 = gdk_region_rectangle(&r),
                    * rg2 = gdk_region_rectangle(&r);

        gdk_region_shrink(rg1,0,0);
        gdk_region_shrink(rg2,2,2);

        gdk_region_subtract(rg1,rg2);
        r = HYRect2GDKRect(rect);
        r.x += windowContent->allocation.x;
        r.y += windowContent->allocation.y;
        gdk_region_destroy (rg2);
        rg2 = gdk_region_rectangle(&r);
        gdk_region_subtract (rg2,rg1);

        gdk_gc_set_clip_region(theCanvas->theContext,rg2);

        theCanvas->_Paint((char*)&rect);
        gdk_gc_set_clip_region(theCanvas->theContext,NULL);
        gdk_region_destroy (rg1);
        gdk_region_destroy (rg2);
    }

    g_object_unref (gc);
    _PaintLFStatus ();
}

//__________________________________________________________________

void        _HYDataPanel::_PaintLFStatus(void)
{
    if (lfID<0) {
        _SimpleList goodP;
        bool    paintOrange = GenerateGoodPartitions (goodP);

        if (goodP.lLength) {
            _PaintTheCircle (paintOrange?orangeButtonIcon:yellowButtonIcon,theWindow);
        } else {
            _PaintTheCircle (redButtonIcon,theWindow);
        }
    } else {
        _PaintTheCircle (greenButtonIcon,theWindow);
    }

}

//__________________________________________________________________

void        _HYDataPanel::_PrintData(void)
{
    _String ptbi ("DataPanel printing has not yet been implemented");
    ReportWarning (ptbi);
}

//__________________________________________________________________

void _HYDataPanel::_VerifyInferMenu(void)
{
    _SimpleList    gp;

    GtkWidget      *inferSubMenu = gtk_item_factory_get_widget (menu_items, "<HY_WINDOW>/Likelihood/Inference");

    if (GenerateGoodPartitions(gp)) {
        if (!inferSubMenu) {
            gtk_item_factory_delete_item (menu_items, "<HY_WINDOW>/Likelihood/Build");
            gtk_item_factory_create_items (menu_items,  sizeof (hyphy_data_window_menu2) / sizeof (hyphy_data_window_menu2[0]), hyphy_data_window_menu2, this);
            inferSubMenu = gtk_item_factory_get_item(menu_items,"<HY_WINDOW>/Likelihood/Inference");
            gtk_menu_reorder_child(GTK_MENU(gtk_item_factory_get_widget(menu_items, "<HY_WINDOW>/Likelihood")),
                                   inferSubMenu,0);
        }
    } else {
        if (inferSubMenu) {
            gtk_item_factory_delete_item (menu_items, "<HY_WINDOW>/Likelihood/Inference");
            gtk_item_factory_create_items (menu_items,  1, &hyphy_data_window_menu[41], this);
            gtk_menu_reorder_child(GTK_MENU(gtk_item_factory_get_widget(menu_items, "<HY_WINDOW>/Likelihood")),
                                   gtk_item_factory_get_item(menu_items, "<HY_WINDOW>/Likelihood/Build"),0);
        }
    }
}

//__________________________________________________________________

void _HYDataPanel::_UpdateLFMenu (void)
{
    bool onOff = lfID>=0;
    if (gtk_item_factory_get_item(menu_items,"<HY_WINDOW>/Likelihood")) {
        gtk_widget_set_sensitive(gtk_item_factory_get_item (menu_items,"<HY_WINDOW>/Likelihood/Display"), onOff);
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(menu_items,HY_DATALF_WIN32_MENU_BASE+1),onOff);
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(menu_items,HY_DATALF_WIN32_MENU_BASE+2),onOff);
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(menu_items,HY_DATALF_WIN32_MENU_BASE+3),onOff);
        gtk_widget_set_sensitive(gtk_item_factory_get_item (menu_items,"<HY_WINDOW>/Data/Simulation"), onOff);
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(menu_items,HY_DATA_WIN32_MENU_BASE+5),onOff);
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(menu_items,HY_DATA_WIN32_MENU_BASE+81),
                                 lfID>=0 && ((_LikelihoodFunction*)likeFuncList (lfID))->GetCategoryVars().lLength);
    }
}

//__________________________________________________________________

void _HYDataPanel::_UpdateSelectionChoices (bool toggle)
{
    gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(menu_items,HY_DATA_WIN32_MENU_BASE+1),toggle);
}

//__________________________________________________________________

void _HYDataPanel::_CopySelectionToClipboard (void)
{
    _HYSequencePane*    sp = (_HYSequencePane*)GetObject(0);
    _String             cbStr (128L,true);

    if (sp->selection.lLength) {
        for (long m=0; m<sp->speciesIndex.lLength; m++) {
            long idx = sp->speciesIndex.lData[m];
            for (long k=0; k<sp->selection.lLength; k++) {
                cbStr << ((_String*)(sp->columnStrings(sp->selection.lData[k])))->sData[idx];
                if (k&&((k+1)%sp->blockWidth==0)) {
                    cbStr << ' ';
                }
            }
            cbStr << '\r';
            cbStr << '\n';
        }
    } else if (sp->vselection.lLength)
        for (long m=0; m<sp->vselection.lLength; m++) {
            cbStr << (_String*)(sp->rowHeaders(sp->speciesIndex(sp->vselection.lData[m])));
            cbStr << '\r';
            cbStr << '\n';
        }

    cbStr.Finalize();


}

//__________________________________________________________________

void _HYDataPanel::_OmitSelectedSpecies (_SimpleList& idx)
{
    GtkMenu * omittedSpecies = GTK_MENU(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Data/Omitted sequences"));
    if (omittedSpecies) {
        _HYSequencePane*    sp = (_HYSequencePane*)GetObject(0);
        long                idxShift = 10001+g_list_length (GTK_MENU_SHELL(omittedSpecies)->children);

        for (long k=0; k<idx.lLength; k++) {
            _String*        thisSpec = (_String*)sp->rowHeaders(idx.lData[k]),
                            entryPath = _String("/Data/Omitted sequences/") & *thisSpec;

            GtkItemFactoryEntry anItem = {entryPath.sData, NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback,  HY_DATA_WIN32_MENU_BASE+idxShift+k, "<Item>"};
            gtk_item_factory_create_item(menu_items,&anItem,this,1);
        }
        gtk_widget_set_sensitive (gtk_item_factory_get_item(menu_items,"<HY_WINDOW>/Data/Omitted sequences"),true);
    }
}

//__________________________________________________________________

void _HYDataPanel::_RestoreOmittedSequence (long index)
{
    if (index>=0) {
        GtkMenu * omittedSpecies = GTK_MENU(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Data/Omitted sequences"));

        long        mic = g_list_length (GTK_MENU_SHELL(omittedSpecies)->children);

        _List       savedPaths;
        for (long k=10003+index+1; k<10001+mic; k++) {
            _String aPath (gtk_item_factory_path_from_widget (gtk_item_factory_get_item_by_action (menu_items, HY_DATA_WIN32_MENU_BASE+k)));
            gtk_item_factory_delete_item (menu_items, aPath.sData);
            aPath.Trim (aPath.Find ('/'),-1);
            savedPaths && & aPath;
        }

        gtk_item_factory_delete_item (menu_items,
                                      gtk_item_factory_path_from_widget (gtk_item_factory_get_item_by_action (menu_items, HY_DATA_WIN32_MENU_BASE+10003+index)));

        for (long k2 = 0; k2 < savedPaths.lLength; k2++) {
            GtkItemFactoryEntry anItem = {((_String*)savedPaths(k2))->sData,
                                          NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback,  HY_DATA_WIN32_MENU_BASE+10003+index+k2, "<Item>"
                                         };
            gtk_item_factory_create_item(menu_items,&anItem,this,1);
        }

        if (mic==3) {
            gtk_widget_set_sensitive (gtk_item_factory_get_item(menu_items,"<HY_WINDOW>/Data/Omitted sequences"),false);
        }
    } else {
        for (long k=0; k< omittedSeqs.lLength; k++)
            gtk_item_factory_delete_item (menu_items,
                                          gtk_item_factory_path_from_widget (gtk_item_factory_get_item_by_action (menu_items, HY_DATA_WIN32_MENU_BASE+10003+k)));

        gtk_widget_set_sensitive (gtk_item_factory_get_item(menu_items,"<HY_WINDOW>/Data/Omitted sequences"),false);
    }
}

//__________________________________________________________________

void _HYDataPanel::_UpdatePartitionOperations (_SimpleList* sl)
{
    bool enFlag = sl->lData[0];
    gtk_widget_set_sensitive (gtk_item_factory_get_item_by_action (menu_items, HY_DATA_WIN32_MENU_BASE), enFlag);
    gtk_widget_set_sensitive (gtk_item_factory_get_item_by_action (menu_items, HY_DATA_WIN32_MENU_BASE+3), enFlag);
}

//__________________________________________________________________

void _HYDataPanel::_UnsetMenuBar(void)
{
}

//__________________________________________________________________
bool _HYDataPanel::_ProcessOSEvent (Ptr vEvent)
{
    static  long        clickH   = 0;
    static  bool        isNavBar = false;

    if (!_HYTWindow::_ProcessOSEvent (vEvent)) {
        _HY_GTK_UI_Message *theMessage = (_HY_GTK_UI_Message*)vEvent;

        if (theMessage->theEvent->type==GDK_BUTTON_PRESS || theMessage->theEvent->type==GDK_MOTION_NOTIFY
                || theMessage->theEvent->type==GDK_2BUTTON_PRESS || theMessage->theEvent->type==GDK_BUTTON_RELEASE) {
            double xc,
                   yc;
            gdk_event_get_coords(theMessage->theEvent,&xc,&yc);

            long ch = xc - windowContent->allocation.x,
                 cv = yc - windowContent->allocation.y;

            long  c  = FindClickedCell(ch, cv);

            if (c<0) {
                return false;
            }

            _HYComponent* thisComponent = (_HYComponent*)components(c);

            if (c==1 || trackMouseComponent == (Ptr)components(1)) { // navBar
                if (theMessage->theEvent->type==GDK_BUTTON_PRESS) {
                    trackMouseComponent = (Ptr)thisComponent;
                }

                ch = ch-componentL.lData[1]-thermRect.left;
                cv = cv-componentT.lData[1]-thermRect.top;

                if (theMessage->theEvent->type==GDK_2BUTTON_PRESS) {
                    NavBarDblClick (ch);
                    return true;
                }
                if (navRect.Contains(ch,cv)||isNavBar) {
                    if (theMessage->theEvent->type==GDK_MOTION_NOTIFY && (((GdkEventMotion*)theMessage->theEvent)->state&GDK_BUTTON1_MASK)) {
                        // skip events more that 0.1 seconds old
                        /*GdkModifierType keyDown;
                        double td1, td2;
                        gdk_display_get_pointer (gdk_display_get_default(),NULL,&td1,&td2,&keyDown);
                        if (keyDown & GDK_BUTTON1_MASK) && */
                        isNavBar = true;
                        guint32 serverTime = gdk_x11_get_server_time (theWindow->window);
                        if (serverTime-gdk_event_get_time(theMessage->theEvent) < 100) {
                            forceUpdateForScrolling = true;
                            SetNavRectCenter (xc-windowContent->allocation.x-componentL.lData[1]-thermRect.left+clickH,0);
                            //printf ("Nav scroll %g, %g, %g\n", xc, yc,0.001*(serverTime-gdk_event_get_time(theMessage->theEvent)));
                            forceUpdateForScrolling = false;
                        }
                    } else {
                        if (theMessage->theEvent->type==GDK_BUTTON_PRESS && (((GdkEventButton*)theMessage->theEvent)->button == 1)) {
                            clickH = (navRect.right+navRect.left)/2-ch;
                        }
                        isNavBar = false;
                    }
                } else {
                    if (theMessage->theEvent->type==GDK_BUTTON_PRESS) {
                        SetNavRectCenter (ch,cv);
                    }
                }

                return true;
            } else if (theMessage->theEvent->type==GDK_BUTTON_PRESS && (((GdkEventButton*)theMessage->theEvent)->button > 0) && c==4) {
                _HYSequencePane* sp2 = (_HYSequencePane*)components (4);
                sp2->ProcessContextualPopUp (ch,cv);
                return true;
            }
        } else if (theMessage->theEvent->type==GDK_KEY_PRESS) {
            GdkEventKey * kp = (GdkEventKey*)theMessage->theEvent;
            if ( kp->keyval==GDK_Left || kp->keyval==GDK_Right ) { // left/right arrow
                _HYSequencePane* sp = (_HYSequencePane*) GetObject (0);
                if (kp->keyval==GDK_Left && sp->startColumn) {
                    sp->HScrollPane (-1);
                } else if (kp->keyval==GDK_Right && sp->endColumn<sp->columnStrings.lLength) {
                    sp->HScrollPane (1);
                }
                return true;
            }
        } else if (theMessage->theEvent->type==GDK_BUTTON_RELEASE) {
            trackMouseComponent = nil;
            isNavBar = false;
            return   true;
        }
    } else {
        return true;
    }

    return false;
}



//EOF
