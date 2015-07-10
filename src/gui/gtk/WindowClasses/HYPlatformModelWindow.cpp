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
    GTK+ Portions of the model window

    Sergei L. Kosakovsky Pond, Spring 2005.
*/

#include "HYModelWindow.h"
#include "HYUtils.h"

#define     HY_MDL_WIN32_MENU_BASE  9200

static GtkItemFactoryEntry hyphy_parameter_model_window_menu[] = {
    { "/_Model",                        NULL,         NULL,           0,    "<Branch>" },
    { "/Model/_Model name",             NULL,         (GtkItemFactoryCallback)hyphy_menu_item_callback,           HY_MDL_WIN32_MENU_BASE,   "<Item>" },
    { "/Model/_Rate variation",         NULL,         NULL,           0, "<Branch>" },
    { "/File/_Save",                    NULL,         NULL,           0, "<Branch>" },
    { "/File/Save/Save",                "<control>S",         (GtkItemFactoryCallback)hyphy_menu_item_callback,           HY_WINDOW_MENU_ID_FILE+1, "<Item>" },
    { "/File/Save/Save as...",          NULL,         (GtkItemFactoryCallback)hyphy_menu_item_callback,           HY_WINDOW_MENU_ID_FILE+3, "<Item>" }
};

//__________________________________________________________________

void _HYModelWindow::_SetMenuBar(void)
{
    if (menu_items && !gtk_item_factory_get_item(menu_items,"<HY_WINDOW>/Model")) {
        gtk_item_factory_delete_item (menu_items, "<HY_WINDOW>/File/Save");
        gtk_item_factory_create_items (menu_items,  sizeof (hyphy_parameter_model_window_menu) / sizeof (hyphy_parameter_model_window_menu[0]),
                                       hyphy_parameter_model_window_menu, this);

        GtkMenu * fileMenu = GTK_MENU (gtk_item_factory_get_widget (menu_items, "<HY_WINDOW>/File"));
        gtk_menu_reorder_child (fileMenu, gtk_item_factory_get_item (menu_items, "<HY_WINDOW>/File/Save"), 0);

        for (long k=0; k<rateOptions.lLength; k++) {
            _String *thisItem = (_String*)rateOptions (k),
                     chopped,
                     type = "<CheckItem>";

            GtkItemFactoryEntry aProcEntry = {NULL,NULL,(GtkItemFactoryCallback)hyphy_menu_item_callback,HY_MDL_WIN32_MENU_BASE+100+k,type.sData};
            chopped = _String("/Model/Rate variation/")&*thisItem;
            aProcEntry.path = chopped.sData;

            gtk_item_factory_create_items (menu_items,  1, &aProcEntry, this);
        }
        gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(gtk_item_factory_get_widget_by_action(menu_items,HY_MDL_WIN32_MENU_BASE+100+rateChoice)),
                                       true);

        gtk_widget_set_sensitive(gtk_item_factory_get_widget (menu_items,"<HY_WINDOW>/Edit/Copy"),  true);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget (menu_items,"<HY_WINDOW>/Edit/Select All"),  true);
    }
}


//__________________________________________________________________

void _HYModelWindow::_UnsetMenuBar(void)
{
}

//__________________________________________________________________


bool        _HYModelWindow::_ProcessMenuSelection (long msel)
{
    if (_HYTWindow::_ProcessMenuSelection(msel)) {
        return true;
    }

    switch (msel) {
    case HY_WINDOW_MENU_ID_FILE+1: // save menu
    case HY_WINDOW_MENU_ID_FILE+3: { // save as menu
        DoSave (msel-HY_WINDOW_MENU_ID_FILE-2);
        return true;
    }

    case HY_WINDOW_MENU_ID_FILE+2: { // print
        _HYTable* t = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,4);
        t->_PrintTable((_HYTable*)GetCellObject (MODEL_MATRIX_ROW-1,4));
        return true;
    }

    case HY_WINDOW_MENU_ID_EDIT+1: { // copy
        DoCopyCell ();
        return true;
    }

    case HY_WINDOW_MENU_ID_EDIT+3: { // paste
        DoPasteToCells();
        return true;
    }

    case HY_WINDOW_MENU_ID_EDIT+5: { // select all
        DoSelectAll();
        return true;
    }

    case HY_MDL_WIN32_MENU_BASE: { // model menu
        DoEditModelName ();
        return true;
    }

    default: { // rate menu
        msel -= HY_MDL_WIN32_MENU_BASE+100;
        if (msel >=0 && gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(gtk_item_factory_get_widget_by_action(menu_items,HY_MDL_WIN32_MENU_BASE+100+msel))))
            if (msel!=rateChoice) {
                gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(gtk_item_factory_get_widget_by_action(menu_items,HY_MDL_WIN32_MENU_BASE+100+rateChoice)),false);
                rateChoice = msel;
                taint = true;
            }
        return true;
    }
    }


    return false;
}

//__________________________________________________________________

void _HYModelWindow::_UpdateEditMenu (bool c, bool p)
{
    gtk_widget_set_sensitive(gtk_item_factory_get_widget (menu_items,"<HY_WINDOW>/Edit/Copy"),  c);
    gtk_widget_set_sensitive(gtk_item_factory_get_widget (menu_items,"<HY_WINDOW>/Edit/Paste"),  p);
}

//__________________________________________________________________

bool _HYModelWindow::_CheckClipboard (void)
{
    /*clipboardString = empty;
    HANDLE  scrapHandle = GetClipboardData (CF_TEXT);

    if (scrapHandle)
    {
        _String cText ((char*)scrapHandle);
        skipWarningMessages = true;
        _Formula f (cText,nil,false);
        skipWarningMessages = false;
        if (f.GetList().lLength)
        {
            clipboardString = cText;
            SyncEditBox ();
        }
    }
    return clipboardString.sLength;*/
    return false;
}

//__________________________________________________________________

void _HYModelWindow::_SetClipboard (void)
{
    //if (clipboardString.sLength)
    //  PlaceStringInClipboard (clipboardString, (Ptr)theWindow);
}

//EOF
