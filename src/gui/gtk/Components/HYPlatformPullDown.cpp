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
#include "HYPullDown.h"
#include "HYUtils.h"
#include "HYEventTypes.h"
#include "HYWindow.h"
#include "HYGraphicPane.h"


_HYColor                    _hyGTKMenuBackground    = {255,255,255};
_HYFont                     defaultPullDownFont     =  {_HY_SANS_FONT,10,HY_FONT_PLAIN};
PangoFontDescription        *defaultPullDownFontPD   = NULL;


//__________________________________________________________________

void        hy_pulldown_selection_callback      (GtkEditable *button, gpointer user_data)
{
    /*_HYPullDown * thePD = (_HYPullDown*)user_data;
    char * val = gtk_editable_get_chars (GTK_EDITABLE(button),0,-1);
    _String  valS = (val);
    long aSel  = thePD->FindMenuItem(valS);
    g_free (val);
    if (aSel>=0 && thePD->selection!=aSel)
    {
        thePD->selection = aSel;
        if (thePD->messageRecipient)
            thePD->messageRecipient->ProcessEvent (generateMenuSelChangeEvent (thePD->GetID(),thePD->selection));
    }*/
}

//__________________________________________________________________

gboolean hy_pulldown_selection_start_callback_event (GtkWidget *widget, GdkEvent* theEvent, gpointer   data )
{
    if (theEvent->type == GDK_BUTTON_PRESS && ((GdkEventButton*)theEvent)->button == 1) {
        _HYPullDown * thePD = (_HYPullDown*)data;
        if (thePD->messageRecipient) {
            thePD->messageRecipient->ProcessEvent (generateMenuOpenEvent(thePD->GetID()));
        }
    }
    return false;
}

//__________________________________________________________________

gboolean hy_pulldown_unmap_event (GtkWidget *widget,gpointer   user_data )
{
    /*if (theEvent->type == GDK_BUTTON_PRESS && ((GdkEventButton*)theEvent)->button == 1)
    {
        _HYPullDown * thePD = (_HYPullDown*)data;
        if (thePD->messageRecipient)
            thePD->messageRecipient->ProcessEvent (generateMenuOpenEvent(thePD->GetID()));
    }*/
    _HYPullDown * thePD = (_HYPullDown*)user_data;
    char * val = gtk_editable_get_chars (GTK_EDITABLE(GTK_COMBO(thePD->theMenu)->entry),0,-1);
    _String  valS = (val);
    long aSel  = thePD->FindMenuItem(valS);
    g_free (val);
    if (aSel>=0 /*&& thePD->selection!=aSel*/) {
        thePD->selection = aSel;
        if (thePD->messageRecipient) {
            thePD->messageRecipient->ProcessEvent (generateMenuSelChangeEvent (thePD->GetID(),thePD->selection));
        }
    }
    return false;
}

//__________________________________________________________________

_HYPlatformPullDown::_HYPlatformPullDown(void)
{
    _HYPullDown * parent = (_HYPullDown*)this;
    theMenu              = gtk_combo_new ();
    backFill             = HYColorToGDKColor((_HYColor) {
        255,255,255
    });
    gtk_combo_set_value_in_list (GTK_COMBO(theMenu),true,false);
    gtk_container_add(GTK_CONTAINER(parent->parentWindow),theMenu);
    gtk_entry_set_editable(GTK_ENTRY(GTK_COMBO(theMenu)->entry),false);
    GList* children = gtk_container_get_children(GTK_CONTAINER(theMenu));
    g_signal_connect (GTK_WIDGET(children->next->data),"event",(GCallback)hy_pulldown_selection_start_callback_event,(_HYPullDown*)this);
    g_signal_connect (GTK_COMBO(theMenu)->entry,"changed",(GCallback)hy_pulldown_selection_callback,(_HYPullDown*)this);
    g_signal_connect (GTK_COMBO(theMenu)->popwin,"hide",(GCallback)hy_pulldown_unmap_event,(_HYPullDown*)this);
    g_list_free (children);
    //gtk_container_set_resize_mode(GTK_CONTAINER(theMenu), GTK_RESIZE_IMMEDIATE);
    selection      = 0;
    cbSelection    = -1;
    if (!defaultPullDownFontPD) {
        defaultPullDownFontPD = pango_font_description_new();
        defaultPullDownFont.size = defaultPullDownFont.size*fontConversionFactor;
        HYFont2PangoFontDesc (defaultPullDownFont,defaultPullDownFontPD);
    }
    gtk_widget_modify_font (theMenu, defaultPullDownFontPD);
    gtk_widget_modify_font (GTK_COMBO(theMenu)->entry, defaultPullDownFontPD);
    gtk_widget_show (theMenu);
    //g_signal_connect (GTK_COMBO(theMenu)->entry,"changed",hy_pulldown_selection_callback,(_HYPullDown*)this);

}

//__________________________________________________________________

_HYPlatformPullDown::~_HYPlatformPullDown(void)
{
}

//__________________________________________________________________
void        _HYPlatformPullDown::_RefreshComboBox  ()
{
    _HYPullDown * theParent = (_HYPullDown*)this;

    if (theMenu && ( cbSelection!=selection || !theParent->IsEnabled())) {
        long mic = theParent->MenuItemCount();
        if (mic) {
            if (selection >= mic) {
                selection = mic-1;
            }

            _String * mItem = theParent->GetMenuItem (selection);

            if (mItem->Equal(&menuSeparator)) {
                mItem = &empty;
            }

            gtk_entry_set_text (GTK_ENTRY(GTK_COMBO(theMenu)->entry), mItem->sData);
            cbSelection = selection;
        }
    }
}

//__________________________________________________________________
void        _HYPlatformPullDown::_AddMenuItem   (_String& newItem, long index)
{
    if (theMenu) {
        GList * singleItem   = g_list_alloc();
        singleItem->prev = singleItem->next = nil;
        if (newItem.Equal(&menuSeparator)) {
            GtkWidget * itemContents = gtk_separator_menu_item_new ();
            gtk_widget_show (itemContents);
            singleItem->data = gtk_list_item_new();
            gtk_container_add((GtkContainer*)singleItem->data,itemContents);
            gtk_combo_set_item_string (GTK_COMBO (theMenu), GTK_ITEM (singleItem->data), "");
            gtk_widget_set_sensitive((GtkWidget*)singleItem->data,false);
            widgetList.InsertElement ((BaseRef)itemContents,2*index,false,false);
        } else {
            _String inItem = newItem;
            if (newItem.beginswith ("(")) {
                inItem.Trim (1,-1);
            }

            GtkWidget * itemContents = gtk_menu_item_new_with_label (inItem.sData);
            gtk_widget_show (itemContents);
            singleItem->data = gtk_list_item_new();
            gtk_container_add((GtkContainer*)singleItem->data,itemContents);
            gtk_combo_set_item_string (GTK_COMBO (theMenu), GTK_ITEM (singleItem->data), inItem.sData);
            gtk_widget_set_sensitive((GtkWidget*)singleItem->data,inItem.sLength == newItem.sLength);
            widgetList.InsertElement ((BaseRef)itemContents,2*index,false,false);
        }
        widgetList.InsertElement ((BaseRef)singleItem->data,2*index,false,false);
        GdkColor convColor = HYColorToGDKColor(_hyGTKMenuBackground);
        gtk_widget_modify_bg ((GtkWidget*)singleItem->data, GTK_STATE_INSENSITIVE, 
        	&convColor);
        gtk_widget_show ((GtkWidget*)singleItem->data);
        if (index<0) {
            gtk_list_append_items (GTK_LIST (GTK_COMBO (theMenu)->list), singleItem);
        } else {
            gtk_list_insert_items (GTK_LIST (GTK_COMBO (theMenu)->list), singleItem, index);
        }

        /*printf ("\nAdding menu item %s at %d\n", newItem.sData, index);
        for (long k = 0; k<widgetList.lLength; k++)
            printf ("%d %s\n", k, GTK_OBJECT_TYPE_NAME (GTK_WIDGET (widgetList(k))));*/
    }

    if (((_HYPullDown*)this)->MenuItemCount()==1||selection==index) {
        cbSelection = -1;
        _RefreshComboBox();
    }
}


//__________________________________________________________________
void        _HYPlatformPullDown::_SetMenuItem   (_String& newItem, long index)
{
    _DeleteMenuItem (index);
    _AddMenuItem (newItem,index);
    if (index==selection) {
        cbSelection = -1;
        _RefreshComboBox();
    }
}


//__________________________________________________________________
void        _HYPlatformPullDown::_MarkItem      (long index, char mark)
{
    if (theMenu && index*2 < widgetList.lLength) {
        GtkWidget * theImage = nil;
        switch (mark) {
        case HY_PULLDOWN_CHECK_MARK:
            theImage = gtk_image_new_from_stock (GTK_STOCK_APPLY,GTK_ICON_SIZE_MENU);
            break;

        case HY_PULLDOWN_BULLET_MARK:
            theImage = gtk_image_new_from_stock (GTK_STOCK_YES,GTK_ICON_SIZE_MENU);
            break;

        case HY_PULLDOWN_DIAMOND_MARK:
            theImage = gtk_image_new_from_stock (GTK_STOCK_ADD,GTK_ICON_SIZE_MENU);
            break;
        }

        GtkWidget * menuLabel = gtk_bin_get_child (GTK_BIN(widgetList(2*index+1)));
        if (theImage) {
            gtk_widget_show (theImage);
        }

        if (GTK_IS_LABEL(menuLabel)) {
            if (theImage) {
                GtkContainer * cW = GTK_CONTAINER (widgetList(2*index+1));
                GtkWidget * hBox = gtk_hbox_new (FALSE,5);
                gtk_widget_show (hBox);
                g_object_ref(menuLabel);
                gtk_container_remove (cW, menuLabel);
                gtk_box_pack_start(GTK_BOX(hBox), theImage, FALSE, FALSE, 0);
                gtk_box_pack_start(GTK_BOX(hBox), menuLabel, FALSE, FALSE, 0);
                g_object_unref(menuLabel);
                gtk_container_add (cW, hBox);
            }
        } else {
            GList     * children      = gtk_container_get_children (GTK_CONTAINER(menuLabel));
            GtkWidget * restoreWidget = GTK_WIDGET(children->next->data);
            g_object_ref (restoreWidget);
            if (theImage) {
                GtkWidget * replaceWidget = GTK_WIDGET(children->data);
                gtk_container_remove(GTK_CONTAINER(menuLabel),replaceWidget);
                gtk_container_remove(GTK_CONTAINER(menuLabel),restoreWidget);
                gtk_box_pack_start(GTK_BOX(menuLabel), theImage, FALSE, FALSE, 0);
                gtk_box_pack_start(GTK_BOX(menuLabel), restoreWidget, FALSE, FALSE, 0);
            } else {
                GtkWidget * parentWidget =  menuLabel->parent;
                gtk_container_remove(GTK_CONTAINER(parentWidget),menuLabel);
                gtk_container_add (GTK_CONTAINER(parentWidget), restoreWidget);
            }
            gtk_widget_show (restoreWidget);
            g_object_unref (restoreWidget);
            g_list_free (children);
        }
    }
}


//__________________________________________________________________
char        _HYPlatformPullDown::_ItemMark      (long index)
{
    // TBI
}

//__________________________________________________________________
void        _HYPlatformPullDown::_DeleteMenuItem  (long index)
{
    if (theMenu) {
        /*printf ("\nDeleting menu item at %d\n", index);
        for (long k = 0; k<widgetList.lLength; k+=2)
            printf ("%d %s\n", k, GTK_OBJECT_TYPE_NAME (GTK_WIDGET (widgetList(k))));*/
        widgetList.Delete (2*index);
        widgetList.Delete (2*index);
        gtk_list_clear_items(GTK_LIST (GTK_COMBO (theMenu)->list),index,index+1);
        if (selection == index) {
            cbSelection = -1;
            //if (selection)
            //  selection--;
            _RefreshComboBox ();
        }
    }
}

//__________________________________________________________________

void        _HYPlatformPullDown::_SetBackColor (_HYColor& c)
{
    backFill = HYColorToGDKColor(c);
}

//__________________________________________________________________
long        _HYPlatformPullDown::_GetSelection (void)
{
    return selection;
}


//__________________________________________________________________
void        _HYPlatformPullDown::_Duplicate (Ptr p)
{
    _HYPullDown * theSource = (_HYPullDown*) p;
    theMenu                 = theSource->theMenu;
    selection               = theSource->selection;
    backFill                = theSource->backFill;
    theSource->theMenu      = nil;
}

//__________________________________________________________________
void        _HYPlatformPullDown::_Update (Ptr p)
{
    _Paint (p);
}

//__________________________________________________________________
void        _HYPlatformPullDown::_SetDimensions (_HYRect r, _HYRect rel)
{
    _HYPullDown* theParent = (_HYPullDown*) this;
    theParent->_HYPlatformComponent::_SetDimensions (r,rel);
    _SetVisibleSize (rel);
}

//__________________________________________________________________
void        _HYPlatformPullDown::_SetVisibleSize (_HYRect rel)
{
    _HYPullDown* theParent = (_HYPullDown*) this;

    menuRect = HYRect2GDKRect(rel);
    menuRect.x     += 3;
    menuRect.width -= 6;

    if (theMenu) {
        if (menuRect.height>20) {
            menuRect.height = 20;
        }

        long            naturalWidth = theParent->_SuggestDimensions ().right;

        if (menuRect.width > naturalWidth + 25) {
            menuRect.width = naturalWidth + 25;
        }
    }

    AlignRectangle (rel, menuRect, theParent->GetAlignFlags());
    gtk_fixed_move (GTK_FIXED(theParent->parentWindow), theMenu, rel.left, rel.top);
    gtk_widget_set_size_request(theMenu,menuRect.width,menuRect.height);
}


//__________________________________________________________________
void        _HYPlatformPullDown::_EnableItem (long theItem, bool toggle)
{
    gtk_widget_set_sensitive(GTK_WIDGET(widgetList(2*theItem)), toggle);
}


//__________________________________________________________________
void        _HYPlatformPullDown::_Paint (Ptr p)
{
    _HYPullDown * theParent = (_HYPullDown*)this;

    GdkRectangle        cRect = HYRect2GDKRect(*(_HYRect*)p);
    _RefreshComboBox();
    if (!(theParent->settings.width&HY_COMPONENT_TRANSP_BG)) {
        if (theParent->parentWindow->window) {
            GdkGC *buttonGC              = gdk_gc_new (theParent->parentWindow->window);
            gdk_gc_set_foreground(buttonGC,&backFill);
            gdk_draw_rectangle(theParent->parentWindow->window,buttonGC,true,cRect.x+theParent->parentWindow->allocation.x,
                               cRect.y+theParent->parentWindow->allocation.y, cRect.width, cRect.height);
            g_object_unref (buttonGC);
        }
    }


    (*theParent)._HYPlatformComponent::_Paint(p);
}

//__________________________________________________________________
_HYRect _HYPullDown::_SuggestDimensions (void)
{
    _HYRect res = {25,100,25,100,HY_COMPONENT_NO_SCROLL};

    if (theMenu) {
        long            naturalWidth = 0;
        PangoLayout *   fontMeasurer = nil;
        for (long k=0; k<widgetList.lLength; k+=2) {
            GtkWidget * dropWidget = GTK_WIDGET (widgetList(k));
            if (!fontMeasurer) {
                fontMeasurer = pango_layout_new (gtk_widget_get_pango_context (dropWidget));
            }

            pango_layout_set_text(fontMeasurer,GetMenuItem(k/2)->sData ,GetMenuItem(k/2)->sLength);
            PangoRectangle extents,
                           log_ext;
            pango_layout_get_pixel_extents (fontMeasurer,&extents,nil);
            if (extents.width > naturalWidth) {
                naturalWidth = extents.width;
            }
        }
        if (fontMeasurer) {
            g_object_unref (fontMeasurer);
        }

        res.right = res.left = naturalWidth+25;
    }

    return res;
}

//__________________________________________________________________

void    _HYPullDown::_SetMenuItemTextStyle (long ID, char style)
{
    // TBI
}

//__________________________________________________________________

bool _HYPullDown::_ProcessOSEvent (Ptr vEvent)
{
    // TBI
    return _HYPlatformComponent::_ProcessOSEvent (vEvent);
}


//EOF
