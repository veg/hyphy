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


#include "HYUtils.h"
#include "HYWindow.h"
#include "hy_strings.h"
#include "batchlan.h"
#include <HYGraphicPane.h>
#include "HYConsoleWindow.h"
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#include <gdk/gdk.h>
#include <gtk/gtk.h>
#include <gdk-pixbuf/gdk-pixbuf.h>


long                lastPopupMenuSelection      = -1;
_List               lastPopupMenuItemStrings;
GtkItemFactoryEntry *lastPopupFactory = nil;

_SimpleList         loadXPM_Index;
_AVLListX           loadedXPMs (&loadXPM_Index);

extern              PangoContext*               screenPContext;

extern              bool                        updateTimer;

extern              clock_t                     timerStart,
                    lastTimer;

//________________________________________________________
_HYRect     GetScreenDimensions (void)
{
    GdkScreen * theScreen = gdk_screen_get_default ();
    return  (_HYRect) {
        0,0,gdk_screen_get_height (theScreen),gdk_screen_get_width (theScreen),0
    };
}

//________________________________________________________
static void hyphy_popup_menu_callback(gpointer  data, guint menuItem, GtkWidget *widget)
{
    lastPopupMenuSelection = menuItem;
    gtk_main_quit();
}

//________________________________________________________
static void hyphy_popup_menu_close(GtkWidget *widget, gpointer userData)
{
    //lastPopupMenuSelection = menuItem;
    //printf ("Close popup\n");
    if (lastPopupMenuSelection < 0) {
        //gtk_widget_destroy (widget);
        gtk_main_quit();
    }
}

//__________________________________________________________________________________

void    ListToPopUpMenu (_List& menuOptions)
{
    lastPopupFactory = new GtkItemFactoryEntry [menuOptions.lLength];
    checkPointer (lastPopupFactory);

    long         sepCounter = 0;

    for (long counter=0; counter<menuOptions.lLength; counter++) {
        _String *postItem = (_String*)(menuOptions(counter)),
                 * pathID   = new _String("/Popup/"),
        * itemType = new _String;

        lastPopupFactory[counter].accelerator    = NULL;
        lastPopupFactory[counter].callback_action = counter;
        lastPopupFactory[counter].callback       = G_CALLBACK(hyphy_popup_menu_callback);
        if (*postItem==_String("SEPARATOR")) {
            *pathID = *pathID & "sep" & (sepCounter++);
            *itemType = "<Separator>";
        } else {
            *pathID = *pathID & *postItem;
            *itemType = "<Item>";
        }
        lastPopupFactory[counter].path = pathID->sData;
        lastPopupFactory[counter].item_type = itemType->sData;
        lastPopupMenuItemStrings << pathID;
        lastPopupMenuItemStrings << itemType;

        DeleteObject (pathID);
        DeleteObject (itemType);
    }
}

//__________________________________________________________________________________

_String HandlePullDown (_List& menuOptions, long l, long t,long startPos)
{
    lastPopupMenuSelection = -1;
    lastPopupMenuItemStrings.Clear();
    if (lastPopupFactory) {
        delete lastPopupFactory;
    }

    if (menuOptions.lLength) {
        ListToPopUpMenu (menuOptions);
        GtkWidget *menu;

        int button                  = 0;
        guint32     event_time      = gtk_get_current_event_time ();

        GtkItemFactory * menu_items = gtk_item_factory_new (GTK_TYPE_MENU_BAR, "<HY_POPUP>", NULL);
        gtk_item_factory_create_items (menu_items,  menuOptions.lLength, lastPopupFactory, menu);

        menu = gtk_item_factory_get_widget (menu_items, "<HY_POPUP>/Popup");
        g_signal_connect (menu, "selection-done", (GCallback)hyphy_popup_menu_close, menu_items);

        lastPopupMenuSelection = -1;

        gtk_menu_popup (GTK_MENU (menu), NULL, NULL, NULL, NULL,  button, event_time);
        gtk_main ();
        gtk_widget_destroy (menu);
        g_object_unref(menu_items);
        if (lastPopupMenuSelection>=0) {
            return *(_String*)menuOptions (lastPopupMenuSelection);
        }

    }

    return empty;
}

//________________________________________________________

long HandlePullDownWithFont (_List& menuOptions, long l, long t,long startPos,_String fName,long fSize)
{
    _String selRes = HandlePullDown (menuOptions, l, t, startPos);
    return menuOptions.Find (&selRes);
}

//________________________________________________________
void        CenterWindow (_HYGuiObject* g)
{
    _HYWindow* w = (_HYWindow*)g;

    _HYRect   screen = GetScreenDimensions();

    long      cleft = 0, ctop = 0;

    if (screen.right>w->right) {
        cleft = (screen.right-w->right)/2;
    }
    if (screen.bottom>w->bottom) {
        ctop = (screen.bottom-w->bottom)/2;
    }

    w->_SetPosition (cleft,ctop);
}

//_________________________________________________________________________

void    MoveConsoleWindow (_HYRect& newLoc)
{
    if (hyphyConsoleWindow) {
        hyphyConsoleWindow->SetPosition (newLoc.left, newLoc.top);
        hyphyConsoleWindow->SetWindowRectangle (newLoc.top, newLoc.left,newLoc.bottom,newLoc.right,true);
    }
}

//________________________________________________________

void    StartBarTimer(void)
{
    timerStart = clock();
    lastTimer = timerStart;
    updateTimer = true;
}

//________________________________________________________

void    StopBarTimer(void)
{
    updateTimer = false;
}


//________________________________________________________
void    ToggleAnalysisMenu (bool running)
{
    if (hyphyConsoleWindow)
        //if (running)
    {
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(hyphyConsoleWindow->menu_items, 21), running);
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(hyphyConsoleWindow->menu_items, 22), running);
        //gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(hyphyConsoleWindow->menu_items, 23), false);
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(hyphyConsoleWindow->menu_items, 24), !running);
        gtk_widget_set_sensitive(gtk_item_factory_get_item(hyphyConsoleWindow->menu_items, "<HY_WINDOW>/Analysis/Results"), !running);
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(hyphyConsoleWindow->menu_items, 30), !running);
    }
    /*else
    {
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(hyphyConsoleWindow->menu_items, 21), false);
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(hyphyConsoleWindow->menu_items, 22), false);
        //gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(hyphyConsoleWindow->menu_items, 23), true);
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(hyphyConsoleWindow->menu_items, 24), true);
        gtk_widget_set_sensitive(gtk_item_factory_get_item(hyphyConsoleWindow->menu_items, "<HY_WINDOW>/Analysis/Results"), true);
        gtk_widget_set_sensitive(gtk_item_factory_get_item_by_action(hyphyConsoleWindow->menu_items, 30), true);
        //SetStatusBarValue (-1,1,0);
        //SetStatusLine ("Idle");
    }*/
}

//________________________________________________________
Ptr     ProcureIconResource (long id)
{
    long cid = loadedXPMs.Find ((BaseRef)id);
    if (cid < 0) {
        _String picFileName = libDirectory & _String("GTKResources/") & id & ".png";
        GError * pixBufError = nil;
        GdkPixbuf*  thePixMap = gdk_pixbuf_new_from_file  (picFileName.sData, &pixBufError);
        if (!thePixMap) {
            picFileName = _String("Failed to load a pixbuf from file. GTK Error ") & pixBufError -> message;
            WarnError (picFileName);
            g_error_free (pixBufError);
            return nil;
        }
        cid = loadedXPMs.Insert ((BaseRef)id, (long)thePixMap);
    }

    return (Ptr)loadedXPMs.GetXtra(cid);
}

//________________________________________________________
long    GetMaxCharWidth (_HYFont& f)
{
    _String dumb ('W');
    return GetVisibleStringWidth (dumb,f);
    // TBI
}

//________________________________________________________
long    GetVisibleStringWidth (_String& s, _HYFont& f)
{
    static PangoLayout* textLayout       = pango_layout_new (screenPContext);
    static PangoFontDescription * fd     = pango_font_description_new ();
    static _HYFont                         stashedFont;

    if (s.sLength) {
        if (stashedFont.size!=f.size || stashedFont.style != f.style || stashedFont.face != f.face) {
            HYFont2PangoFontDesc(f,fd);
            pango_layout_set_width (textLayout,-1);
            pango_layout_set_font_description(textLayout,fd);
            stashedFont = f;
        }
        pango_layout_set_text (textLayout, s.sData,s.sLength);
        //PangoRectangle charPos;
        //pango_layout_index_to_pos (textLayout,s.sLength-1,&charPos);
        //return PANGO_PIXELS(charPos.x+charPos.width);
        PangoRectangle extents,
                       logical_ext;

        pango_layout_get_pixel_extents (textLayout,&extents,&logical_ext);
        return logical_ext.width;
    } else {
        return 0;
    }

}

//________________________________________________________

bool handleGUI (bool yield)
{
    return gtk_main_iteration_do (!yield);
}

//_________________________________________________________________________
_HYColor    GetDialogBackgroundColor (void)
{
    if (hyphyConsoleWindow) {
        GtkStyle* wstyle = gtk_widget_get_style (hyphyConsoleWindow->theWindow);
        if (wstyle) {
            GdkColor  bgc = wstyle->bg[0];
            return (_HYColor) {
                bgc.red/256, bgc.green/256, bgc.blue/256
            };
        }
    }

    return (_HYColor) {
        255,255,255
    };
}

//_________________________________________________________________________

void    GenerateFontList (_List& fonts)
{
    fonts.Clear();
    PangoFontFamily **families;
    PangoFontFamily *match_family = NULL;
    gint n_families, i;
    pango_context_list_families ((screenPContext),&families, &n_families);

    for (i=0; i<n_families; i++) {
        fonts.AppendNewInstance(new _String(pango_font_family_get_name (families[i])));
    }
    fonts.Sort();
    g_free (families);
}

//__________________________________________________________________________________
char    ScanDirectoryForFileNames (_String& source, _List& rec, bool recurse)
{
    DIR * dirPntr = opendir (source.sData);
    if (dirPntr)
        /* source directory exists */
    {
        struct dirent * curEntry = nil;

        while (curEntry = readdir (dirPntr))
            // index thru the items
        {
            _String childDir (curEntry->d_name);
            if (childDir.sData[0] != '.') { // invisible file
                childDir = source & '/' & childDir;
                struct stat fileInfo;
                if (stat(childDir.sData, &fileInfo) == 0)
                    //if (curEntry->d_type == DT_DIR)
                {
                    if (S_ISDIR(fileInfo.st_mode))
                        // a directory
                    {
                        if (recurse) {
                            ScanDirectoryForFileNames (childDir,rec,true);
                        }
                    } else
                        //if (curEntry->d_type == DT_REG)
                        if (S_ISREG(fileInfo.st_mode)) {
                            rec && & childDir;
                        }
                }
            }
        }
        closedir (dirPntr);
    }
    return '/';
}

//_________________________________________________________________________

bool hyPromptForADirectory = false;

//

void hyphy_store_filename (GtkWidget *widget, gpointer user_data)
{
    GtkWidget *file_selector = GTK_WIDGET (user_data);
    const gchar *selected_filename;
    selected_filename = gtk_file_selection_get_filename (GTK_FILE_SELECTION (file_selector));
    *argFileName = selected_filename;
    DIR * dirPntr = opendir (argFileName->sData);
    if (dirPntr) {
        closedir (dirPntr);
        if (hyPromptForADirectory) {
            gtk_widget_destroy (file_selector);
            gtk_main_quit();
        } else {
            *argFileName = empty;
        }
    } else {
        if (hyPromptForADirectory) {
            *argFileName = empty;
        } else {
            gtk_widget_destroy (file_selector);
            gtk_main_quit();
        }
    }
}



//_________________________________________________________________________

void hyphy_kill_dialog_and_exit (GtkWidget *widget, gpointer user_data)
{
    GtkWidget *file_selector = GTK_WIDGET (user_data);
    gtk_widget_destroy (file_selector);
    gtk_main_quit();
}

//_________________________________________________________________________

void hyphy_grab_option (GtkWidget *widget, gpointer user_data)
{
    GtkOptionMenu *file_selector = GTK_OPTION_MENU (widget);
    *((long*)user_data) = gtk_option_menu_get_history(file_selector);
}

//_________________________________________________________________________

void hyphy_exit_loop (GtkWidget *widget, gpointer user_data)
{
    gtk_main_quit();
}

//_________________________________________________________________________
bool    PopUpFileDialog(_String ps, _String* defaultLocation)
{
    static _String   lastOpenFilePath;

    if (!argFileName) {
        argFileName = new _String;
    } else {
        *argFileName = empty;
    }

    GtkWidget * fileSelector = gtk_file_selection_new (ps.sData);
    if (defaultLocation) {
        gtk_file_selection_set_filename (GTK_FILE_SELECTION (fileSelector), defaultLocation->sData);
    } else if (lastOpenFilePath.sLength) {
        gtk_file_selection_set_filename (GTK_FILE_SELECTION (fileSelector), lastOpenFilePath.sData);
    }

    if (defFileNameValue.sLength) {
        gtk_entry_set_text(GTK_ENTRY(GTK_FILE_SELECTION (fileSelector)->selection_entry), defFileNameValue.sData);
    }

    gtk_file_selection_hide_fileop_buttons(GTK_FILE_SELECTION (fileSelector));
    gtk_file_selection_set_select_multiple (GTK_FILE_SELECTION (fileSelector),false);

    g_signal_connect (GTK_FILE_SELECTION (fileSelector)->ok_button,
                      "clicked",
                      G_CALLBACK (hyphy_store_filename),
                      fileSelector);


    g_signal_connect (GTK_FILE_SELECTION (fileSelector)->cancel_button,
                      "clicked",
                      G_CALLBACK (hyphy_kill_dialog_and_exit),
                      fileSelector);

    g_signal_connect_swapped (fileSelector,
                              "close",
                              G_CALLBACK (hyphy_exit_loop),
                              fileSelector);

    gtk_widget_show (fileSelector);
    gtk_window_set_modal (GTK_WINDOW(fileSelector),true);
    gtk_main ();
    printf ("\nPopUpFileDialog:%s\n", argFileName->sData);
    if (argFileName->sLength) {
        lastOpenFilePath = argFileName->Cut(0,argFileName->FindBackwards('/',0,-1));
    }

    return argFileName->sLength;
}

//________________________________________________________

long    SaveFileWithPopUp (_String& fileName, _String& prompt, _String& defFileName, _String& listLabel, _List& menuOptions)
{
    static _String   lastOpenFilePath;
    long      optionChoice = 0;

    if (!argFileName) {
        argFileName = new _String;
    } else {
        *argFileName = empty;
    }

    GtkWidget * fileSelector = gtk_file_selection_new (prompt.sData);
    _String fName = lastOpenFilePath & defFileName;
    gtk_file_selection_set_filename (GTK_FILE_SELECTION (fileSelector), fName.sData);

    gtk_file_selection_set_select_multiple (GTK_FILE_SELECTION (fileSelector),false);

    GtkWidget * fileFormatOptions = gtk_menu_new ();
    for (long k = 0; k<menuOptions.lLength; k++) {
        GtkWidget * optionWidget      = gtk_menu_item_new_with_label (((_String*)menuOptions(k))->sData);
        gtk_menu_shell_append (GTK_MENU_SHELL (fileFormatOptions), optionWidget);
        gtk_widget_show (optionWidget);
    }
    GtkWidget * fpd = gtk_option_menu_new ();
    gtk_option_menu_set_menu (GTK_OPTION_MENU (fpd), fileFormatOptions);
    gtk_option_menu_set_history(GTK_OPTION_MENU (fpd), 0);
    gtk_widget_show (fpd);
    GtkWidget* secLabel = gtk_label_new (listLabel.sData);
    GtkWidget* newBox   = gtk_hbox_new (true,0);
    gtk_box_pack_start(GTK_BOX(newBox),secLabel,false,false,0);
    gtk_box_pack_start(GTK_BOX(newBox),fpd,true,false,0);
    gtk_widget_show (secLabel);
    gtk_widget_show (newBox);

    gtk_box_pack_end(GTK_BOX(GTK_DIALOG(fileSelector)->vbox),newBox,true,true,0);
    //gtk_button_box_set_child_secondary (GTK_BUTTON_BOX(GTK_DIALOG(fileSelector)->action_area),newBox,true);

    g_signal_connect (GTK_FILE_SELECTION (fileSelector)->ok_button,
                      "clicked",
                      G_CALLBACK (hyphy_store_filename),
                      fileSelector);

    g_signal_connect (fpd,
                      "changed",
                      G_CALLBACK (hyphy_grab_option),
                      &optionChoice);


    g_signal_connect (GTK_FILE_SELECTION (fileSelector)->cancel_button,
                      "clicked",
                      G_CALLBACK (hyphy_kill_dialog_and_exit),
                      fileSelector);

    g_signal_connect_swapped (fileSelector,
                              "close",
                              G_CALLBACK (hyphy_exit_loop),
                              fileSelector);

    gtk_widget_show (fileSelector);
    gtk_window_set_modal (GTK_WINDOW(fileSelector),true);
    gtk_main ();
    fileName = *argFileName;
    if (argFileName->sLength) {
        lastOpenFilePath = argFileName->Cut(0,argFileName->FindBackwards('/',0,-1));
    }

    if (argFileName->sLength) {
        return optionChoice;
    }
    return -1;
}

//_________________________________________________________________________

void hyphy_store_color (GtkWidget *widget, gpointer user_data)
{
    _HYColor *color_ref = (_HYColor*)user_data;
    GdkColor newColor;
    gtk_color_selection_get_current_color (GTK_COLOR_SELECTION(widget),&newColor);
    *color_ref = (_HYColor) {
        newColor.red/256,newColor.green/256,newColor.blue/256
    };
}

//________________________________________________________

_HYColor        SelectAColor (_HYColor& currentColor, _String& prompt)
{
    _HYColor newColor = currentColor;

    GtkWidget * fileSelector = gtk_color_selection_dialog_new(prompt.sData);
    GdkColor    startColor = HYColorToGDKColor(currentColor);

    gtk_color_selection_set_current_color (GTK_COLOR_SELECTION(GTK_COLOR_SELECTION_DIALOG(fileSelector)->colorsel),
                                           &startColor);

    g_signal_connect (GTK_COLOR_SELECTION_DIALOG (fileSelector)->ok_button,
                      "clicked",
                      G_CALLBACK (hyphy_kill_dialog_and_exit),
                      fileSelector);


    g_signal_connect (GTK_COLOR_SELECTION_DIALOG (fileSelector)->colorsel,
                      "color_changed",
                      G_CALLBACK (hyphy_store_color),
                      &newColor);

    g_signal_connect (GTK_COLOR_SELECTION_DIALOG (fileSelector)->cancel_button,
                      "clicked",
                      G_CALLBACK (hyphy_kill_dialog_and_exit),
                      fileSelector);

    g_signal_connect_swapped (GTK_COLOR_SELECTION_DIALOG(fileSelector),
                              "close",
                              G_CALLBACK (hyphy_exit_loop),
                              fileSelector);

    gtk_widget_show (fileSelector);
    gtk_window_set_modal (GTK_WINDOW(fileSelector),true);
    gtk_main ();
    return newColor;
}

//________________________________________________________
char        YesNoCancelPrompt (_String& prompt)
{
    GtkWidget * theDialog = gtk_message_dialog_new (NULL, GTK_DIALOG_MODAL, GTK_MESSAGE_QUESTION, GTK_BUTTONS_NONE,
                            "%s",prompt.sData);
    gtk_dialog_add_buttons (GTK_DIALOG(theDialog), "Cancel", 2, "Yes", 1, "No", 0, NULL);
    char res =  gtk_dialog_run (GTK_DIALOG(theDialog));
    gtk_widget_destroy (theDialog);
    return res;
}

//_________________________________________________________________________
void    DelayNMs (long ms)
{
    g_usleep(ms*1000);
}

//_________________________________________________________________________
void    PositionWindow          (_HYGuiObject* twp, _String* args)
{
    _List * argL = args->Tokenize (",");
    _HYWindow*   tw = (_HYWindow*)twp;
    if (argL->lLength>=4) {
        long R[5],
             k;

        for (k=0; k<4; k++) {
            R[k] = ((_String*)(*argL)(k))->toNum();
        }
        if (argL->lLength>4) {
            R[4] = ((_String*)(*argL)(4))->toNum();
        } else {
            R[4] = 0;
        }

        _HYRect   wR = GetScreenDimensions  (),
                  wiR;
        long      W[4] = {wR.left,wR.top, wR.right, wR.bottom};
        for (k=0; k<4; k++)
            if (R[k]<0) {
                R[k] = W[k] + ((k<2)?-1:1)*R[k];
            }

        wiR.left    = R[0];
        wiR.right   = R[2];
        wiR.top     = R[1];
        wiR.bottom  = R[3];

        if (wiR.left>=wiR.right) {
            wiR.right = 1+wiR.left;
        }
        if (wiR.top>=wiR.bottom) {
            wiR.bottom = 1+wiR.top;
        }

        tw->SetPosition        (wiR.left,wiR.top);
        tw->SetWindowRectangle (0,0,wiR.bottom-wiR.top,wiR.right-wiR.left);

        if (R[4]>0) {
            wiR.top         = wiR.bottom+2;
            wiR.bottom      = wR.bottom - 2;
            wiR.left        = 5;
            wiR.right       = wR.right - 2;
            MoveConsoleWindow (wiR);
        }

    }

    DeleteObject (argL);
}

//__________________________________________________________________________________
void    PlaceStringInClipboard (_String& res,Ptr )
{
    gtk_clipboard_set_text (gtk_clipboard_get(GDK_SELECTION_CLIPBOARD),res.sData,res.sLength);
}

//_________________________________________________________________________

_String ChooseAFolder       (_String& prompt)
{
    hyPromptForADirectory = true;
    PopUpFileDialog (prompt, NULL);
    hyPromptForADirectory = false;
    return *argFileName & '/';
    //return empty;
}
