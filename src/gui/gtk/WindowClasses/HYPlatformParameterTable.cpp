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
    GTK+ Portions of the parameter table class

    Sergei L. Kosakovsky Pond, Spring 2005.
*/

#include "HYParameterTable.h"
#include "HYUtils.h"

#define     HY_PT_WIN32_MENU_BASE  9000

static GtkItemFactoryEntry hyphy_parameter_table_window_menu[] = {
    { "/_Likelihood",                       NULL,         NULL,           0,    "<Branch>" },
    { "/Likelihood/_View Options",          NULL,         NULL,           0,    "<Branch>" },
    { "/Likelihood/View Options/_Local parameters",           NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,                HY_PT_WIN32_MENU_BASE+10, "<CheckItem>" },
    { "/Likelihood/View Options/_Global parameters",          NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,                HY_PT_WIN32_MENU_BASE+11, "<CheckItem>" },
    { "/Likelihood/View Options/_Constrained parameters",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,                HY_PT_WIN32_MENU_BASE+12, "<CheckItem>" },
    { "/Likelihood/View Options/_Rate classes",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,              HY_PT_WIN32_MENU_BASE+13, "<CheckItem>" },
    { "/Likelihood/View Options/_Trees",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,             HY_PT_WIN32_MENU_BASE+14, "<CheckItem>" },
    { "/Likelihood/sep1",                                     NULL,    NULL,           0,   "<Separator>" },
    { "/Likelihood/_Recalculate LF",     "<control>U",    (GtkItemFactoryCallback)hyphy_menu_item_callback,             HY_PT_WIN32_MENU_BASE, "<Item>" },
    { "/Likelihood/_Optimize LF",     "<control>T",    (GtkItemFactoryCallback)hyphy_menu_item_callback,                HY_PT_WIN32_MENU_BASE+1, "<Item>" },
    { "/Likelihood/sep2",                                     NULL,    NULL,           0,   "<Separator>" },
    { "/Likelihood/_Enter command",     NULL,       (GtkItemFactoryCallback)hyphy_menu_item_callback,               HY_PT_WIN32_MENU_BASE+2, "<Item>" },
    { "/Likelihood/Remove _unused parameters",     NULL,        (GtkItemFactoryCallback)hyphy_menu_item_callback,               HY_PT_WIN32_MENU_BASE+3, "<Item>" },
    { "/Likelihood/sep3",                                     NULL,    NULL,           0,   "<Separator>" },
    { "/Likelihood/Covariance, sampler and C_I",     NULL,      (GtkItemFactoryCallback)hyphy_menu_item_callback,               HY_PT_WIN32_MENU_BASE+4, "<Item>" },
    { "/Likelihood/Likelihood pro_file plot",     NULL,     (GtkItemFactoryCallback)hyphy_menu_item_callback,               HY_PT_WIN32_MENU_BASE+5, "<Item>" },
    { "/Likelihood/sep4",                                     NULL,    NULL,           0,   "<Separator>" },
    { "/Likelihood/Cate_gories Processor",     NULL,        (GtkItemFactoryCallback)hyphy_menu_item_callback,               HY_PT_WIN32_MENU_BASE+6, "<Item>" },
    { "/Likelihood/sep5",                                     NULL,    NULL,           0,   "<Separator>" },
    { "/Likelihood/Sele_ct parameters",     NULL,       (GtkItemFactoryCallback)hyphy_menu_item_callback,               HY_PT_WIN32_MENU_BASE+7, "<Item>" },
    { "/Likelihood/_Open selection in table",     NULL,     (GtkItemFactoryCallback)hyphy_menu_item_callback,               HY_PT_WIN32_MENU_BASE+8, "<Item>" }
};

//__________________________________________________________________

void _HYParameterTable::_SetMenuBar(void)
{
    if (menu_items && !gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Likelihood")) {
        gtk_item_factory_create_items (menu_items,  sizeof (hyphy_parameter_table_window_menu) / sizeof (hyphy_parameter_table_window_menu[0]),
                                       hyphy_parameter_table_window_menu, this);
        _UpdateViewMenu();
        gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_WINDOW_MENU_ID_EDIT),undoCommands.lLength);
    }
}

//__________________________________________________________________


bool        _HYParameterTable::_ProcessMenuSelection (long msel)
{

    _HYTable*   table = (_HYTable*)GetCellObject(HY_PARAMETER_TABLE_TABLE_ROW,0);
    bool        res = false;

    switch (msel) {
    case HY_WINDOW_MENU_ID_FILE+1: { // save
        DoSave ();
        res = true;
        break;
    }
    case HY_WINDOW_MENU_ID_FILE+2: { // print
        _SimpleList columns,
                    sel;
        columns << 0;
        columns << 1;
        columns << 2;
        columns << 3;
        table->GetSelection (sel);
        char    resp = 3;
        if (sel.lLength) {
            _String pr ("Would you like to print only the selected cells? (Click \"No\" to print the entire table).");
            resp = YesNoCancelPrompt (pr);
        }
        if (resp == 3) {
            table->_PrintTable(columns,(_HYTable*)GetCellObject(HY_PARAMETER_TABLE_TABLE_ROW-1,0));
        } else if (resp == 1) {
            _SimpleList rows;
            for (long k = 0; k < sel.lLength; k+=4) {
                rows << sel.lData[k]/4;
            }
            table->_PrintTable(columns,rows,(_HYTable*)GetCellObject(HY_PARAMETER_TABLE_TABLE_ROW-1,0));
        }
        res = true;
        break;
    }

    case HY_WINDOW_MENU_ID_EDIT: { // undo
        UndoCommand();
        _UpdateUndoMenu (nil,nil);
        res = true;
        break;
    }

    case HY_WINDOW_MENU_ID_EDIT+5: { // undo
        SelectAll();
        res = true;
        break;
    }

    case HY_PT_WIN32_MENU_BASE+10:
    case HY_PT_WIN32_MENU_BASE+11:
    case HY_PT_WIN32_MENU_BASE+12:
    case HY_PT_WIN32_MENU_BASE+13: {
        res = true;
        char   toggleFlag;
        GtkWidget * checkMI = gtk_item_factory_get_widget_by_action(menu_items,msel);
        msel -= HY_PT_WIN32_MENU_BASE+10;
        switch (msel) {
        case 0:
            toggleFlag = HY_PARAMETER_TABLE_VIEW_LOCAL;
            break;
        case 1:
            toggleFlag = HY_PARAMETER_TABLE_VIEW_GLOBAL;
            break;
        case 2:
            toggleFlag = HY_PARAMETER_TABLE_VIEW_CONSTRAINED;
            break;
        case 3:
            toggleFlag = HY_PARAMETER_TABLE_VIEW_CATEGORY;
            break;

        }
        if (gtk_check_menu_item_get_active (GTK_CHECK_MENU_ITEM(checkMI)) != (bool)viewOptions&toggleFlag) {
            if (viewOptions&toggleFlag) {
                if (viewOptions-toggleFlag) {
                    viewOptions-=toggleFlag;
                } else {
                    break;
                }
            } else {
                viewOptions+=toggleFlag;
            }

            gtk_check_menu_item_set_active (GTK_CHECK_MENU_ITEM(checkMI), viewOptions&toggleFlag);
            ConstructTheTable();
            SetWindowRectangle (top,left,bottom,right);
        }
        break;
    }

    case HY_PT_WIN32_MENU_BASE+1:
        OptimizeLikelihoodFunction();
        res = true;
        break;

    case HY_PT_WIN32_MENU_BASE+2:
        DoEnterConstraint ();
        res = true;
        break;

    case HY_PT_WIN32_MENU_BASE+3:
        DoCleanUp   ();
        res = true;
        break;

    case HY_PT_WIN32_MENU_BASE+4:
        HandleVarianceEstimates ();
        res = true;
        break;

    case HY_PT_WIN32_MENU_BASE+5:
        HandleProfilePlot   ();
        res = true;
        break;

    case HY_PT_WIN32_MENU_BASE+6:
        HandleCategories    ();
        res = true;
        break;

    case HY_PT_WIN32_MENU_BASE+7:
        HandleSelectParameters  ();
        res = true;
        break;

    case HY_PT_WIN32_MENU_BASE+8:
        HandleOpenInChart   ();
        res = true;
        break;
    }

    if (!res) {
        res = _HYTWindow::_ProcessMenuSelection(msel);
    }
    return res;
}

//__________________________________________________________________

void _HYParameterTable::_UpdateViewMenu(void)
{
    GtkWidget * vlocal  = gtk_item_factory_get_widget_by_action(menu_items,HY_PT_WIN32_MENU_BASE+10),
                * vglobal = gtk_item_factory_get_widget_by_action(menu_items,HY_PT_WIN32_MENU_BASE+11),
                  * vconstr = gtk_item_factory_get_widget_by_action(menu_items,HY_PT_WIN32_MENU_BASE+12),
                    * vcateg  = gtk_item_factory_get_widget_by_action(menu_items,HY_PT_WIN32_MENU_BASE+13),
                      * vtree   = gtk_item_factory_get_widget_by_action(menu_items,HY_PT_WIN32_MENU_BASE+14);

    gtk_widget_set_sensitive (vglobal,avViewOptions&HY_PARAMETER_TABLE_VIEW_GLOBAL);
    gtk_widget_set_sensitive (vconstr,avViewOptions&HY_PARAMETER_TABLE_VIEW_CONSTRAINED);
    gtk_widget_set_sensitive (vcateg,avViewOptions&HY_PARAMETER_TABLE_VIEW_CATEGORY);
    gtk_widget_set_sensitive (vtree,avViewOptions&HY_PARAMETER_TABLE_VIEW_TREES);

    if (gtk_check_menu_item_get_active (GTK_CHECK_MENU_ITEM (vlocal)) != (viewOptions&HY_PARAMETER_TABLE_VIEW_LOCAL)) {
        gtk_check_menu_item_set_active (GTK_CHECK_MENU_ITEM (vlocal), viewOptions&HY_PARAMETER_TABLE_VIEW_LOCAL);
    }
    if (gtk_check_menu_item_get_active (GTK_CHECK_MENU_ITEM (vglobal)) != (viewOptions&HY_PARAMETER_TABLE_VIEW_GLOBAL)) {
        gtk_check_menu_item_set_active (GTK_CHECK_MENU_ITEM (vglobal), viewOptions&HY_PARAMETER_TABLE_VIEW_GLOBAL);
    }
    if (gtk_check_menu_item_get_active (GTK_CHECK_MENU_ITEM (vconstr)) != (viewOptions&HY_PARAMETER_TABLE_VIEW_CONSTRAINED)) {
        gtk_check_menu_item_set_active (GTK_CHECK_MENU_ITEM (vconstr), viewOptions&HY_PARAMETER_TABLE_VIEW_CONSTRAINED);
    }
    if (gtk_check_menu_item_get_active (GTK_CHECK_MENU_ITEM (vcateg)) != (viewOptions&HY_PARAMETER_TABLE_VIEW_CATEGORY)) {
        gtk_check_menu_item_set_active (GTK_CHECK_MENU_ITEM (vcateg), viewOptions&HY_PARAMETER_TABLE_VIEW_CATEGORY);
    }
    if (gtk_check_menu_item_get_active (GTK_CHECK_MENU_ITEM (vtree)) != (viewOptions&HY_PARAMETER_TABLE_VIEW_TREES)) {
        gtk_check_menu_item_set_active (GTK_CHECK_MENU_ITEM (vtree), viewOptions&HY_PARAMETER_TABLE_VIEW_TREES);
    }
}


//__________________________________________________________________

void _HYParameterTable::_UnsetMenuBar(void)
{
}

//__________________________________________________________________

void _HYParameterTable::_UpdateUndoMenu(_String* command, _String* desc)
{
    GtkWidget * undoItem = gtk_item_factory_get_widget_by_action(menu_items,HY_WINDOW_MENU_ID_EDIT);
    if (command&&desc) {
        undoCommands        &&  command;
        undoDescriptions    &&  desc;
        gtk_widget_set_sensitive (undoItem,true);
        gtk_label_set_text (GTK_LABEL (gtk_bin_get_child(GTK_BIN(undoItem))), desc->sData);
    } else {
        if (undoDescriptions.lLength==0) {
            gtk_widget_set_sensitive (undoItem,false);
            gtk_label_set_text (GTK_LABEL (gtk_bin_get_child(GTK_BIN(undoItem))), "Can't _undo");
        } else {
            _String       temp =  *(_String*)undoDescriptions(undoDescriptions.lLength-1);
            gtk_widget_set_sensitive (undoItem,true);
            gtk_label_set_text (GTK_LABEL (gtk_bin_get_child(GTK_BIN(undoItem))), temp.sData);
        }
    }
}

//EOF
