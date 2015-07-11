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
    Tree  Panel Object  for GTK

    Sergei L. Kosakovsky Pond, March 2005
*/

#include "HYTreePanel.h"
#include "HYUtils.h"
#include <gdk/gdkkeysyms.h>

_String     saveForTreesPrompt ("Save Tree As:");

GdkColor    navPen = HYColorToGDKColor((_HYColor)
{
    255,151,51
});

#define     HY_TREE_GTK_MENU_BASE  7000


//__________________________________________________________________

static GtkItemFactoryEntry hyphy_tree_window_menu[] = {
    { "/Edit/sep2",         NULL,         NULL,           0,    "<Separator>" },
    { "/Edit/Search and Repla_ce",     "<control>F",    (GtkItemFactoryCallback)hyphy_menu_item_callback,               HY_WINDOW_MENU_ID_EDIT+6, "<Item>" },
    { "/Edit/Search and Re_place in Selection",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,          HY_WINDOW_MENU_ID_EDIT+7, "<Item>" },
    { "/_Tree",         NULL,         NULL,           0,                    "<Branch>" },
    { "/Tree/Tip _Labels", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_TREE_GTK_MENU_BASE, "<CheckItem>"},
    { "/Tree/I_nternal Labels", NULL, (GtkItemFactoryCallback)hyphy_menu_item_callback, HY_TREE_GTK_MENU_BASE+1, "<CheckItem>"},
    { "/Tree/sep2",         NULL,         NULL,           0,        "<Separator>" },
    { "/Tree/S_wap Subtrees",     "<control>1",    (GtkItemFactoryCallback)hyphy_menu_item_callback,                HY_TREE_GTK_MENU_BASE+2, "<Item>" },
    { "/Tree/_Collapse Branch",     "<control>2",    (GtkItemFactoryCallback)hyphy_menu_item_callback,              HY_TREE_GTK_MENU_BASE+3, "<Item>" },
    { "/Tree/_Join",     "<control>3",    (GtkItemFactoryCallback)hyphy_menu_item_callback,             HY_TREE_GTK_MENU_BASE+4, "<Item>" },
    { "/Tree/_Graft A Tip",     "<control>4",    (GtkItemFactoryCallback)hyphy_menu_item_callback,              HY_TREE_GTK_MENU_BASE+5, "<Item>" },
    { "/Tree/_Reroot",     "<control>5",    (GtkItemFactoryCallback)hyphy_menu_item_callback,               HY_TREE_GTK_MENU_BASE+6, "<Item>" },
    { "/Tree/_Flip Tip Ordering",     "<control>6",    (GtkItemFactoryCallback)hyphy_menu_item_callback,                HY_TREE_GTK_MENU_BASE+7, "<Item>" },
    { "/Tree/sep3",         NULL,         NULL,           0,        "<Separator>" },
    { "/Tree/Select Branches",          NULL,         NULL,           0,                    "<Branch>" },
    { "/Tree/Select Branches/Select _Entire Subtree",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,                HY_TREE_GTK_MENU_BASE+8, "<Item>" },
    { "/Tree/Select Branches/Select _Incomplete Branchhes",     "<control>I",    (GtkItemFactoryCallback)hyphy_menu_item_callback,              HY_TREE_GTK_MENU_BASE+9, "<Item>" },
    { "/Tree/Select Branches/Select Bran_ches Without Models",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,               HY_TREE_GTK_MENU_BASE+10, "<Item>" },
    { "/Tree/Select Branches/Select Branches By _Name",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,              HY_TREE_GTK_MENU_BASE+17, "<Item>" },
    { "/Tree/Select Branches/Select Branches By _Length",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,                HY_TREE_GTK_MENU_BASE+18, "<Item>" },
    { "/Tree/Select Branches/_Invert Selection",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,             HY_TREE_GTK_MENU_BASE+30, "<Item>" },
    { "/Tree/Select Branches/_Grow Selection",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,               HY_TREE_GTK_MENU_BASE+31, "<Item>" },
    { "/Tree/Select Branches/_Map Selection to Datapanel",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,               HY_TREE_GTK_MENU_BASE+32, "<Item>" },
    { "/Tree/Select Branches/_Find selection in another tree",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,               HY_TREE_GTK_MENU_BASE+33, "<Item>" },
    { "/Tree/Edit Prope_rties",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,              HY_TREE_GTK_MENU_BASE+11, "<Item>" },
    { "/Tree/sep4",         NULL,         NULL,           0,        "<Separator>" },
    { "/Tree/Optimi_ze Again",     "<control>T",    (GtkItemFactoryCallback)hyphy_menu_item_callback,               HY_TREE_GTK_MENU_BASE+12, "<Item>" },
    { "/Tree/Sh_ow Parameters in Table",     "<control>H",    (GtkItemFactoryCallback)hyphy_menu_item_callback,             HY_TREE_GTK_MENU_BASE+13, "<Item>" },
    { "/Tree/sep5",         NULL,         NULL,           0,        "<Separator>" },
    { "/Tree/Tr_ee Display Options...",     NULL,    (GtkItemFactoryCallback)hyphy_menu_item_callback,              HY_TREE_GTK_MENU_BASE+14, "<Item>" },
    { "/Tree/Branch Labels",            NULL,         NULL,           0,                    "<Branch>" },
    { "/Tree/Branch Labels/_Above Branches",     "<control>8",    (GtkItemFactoryCallback)hyphy_menu_item_callback,             HY_TREE_GTK_MENU_BASE+24, "<Item>" },
    { "/Tree/Branch Labels/_Below Branches",     "<control>9",    (GtkItemFactoryCallback)hyphy_menu_item_callback,             HY_TREE_GTK_MENU_BASE+25, "<Item>" },
    { "/Tree/sep6",         NULL,         NULL,           0,        "<Separator>" },
    { "/Tree/Show Rate Matri_x",            NULL,         (GtkItemFactoryCallback)hyphy_menu_item_callback,           HY_TREE_GTK_MENU_BASE+15,                 "<Item>" },
    { "/Tree/S_how Transition Matrix",          NULL,         (GtkItemFactoryCallback)hyphy_menu_item_callback,           HY_TREE_GTK_MENU_BASE+16,                 "<Item>" },
    { "/Tree/sep7",         NULL,         NULL,           0,        "<Separator>" },
    { "/Tree/Pairwise Distan_ces",          NULL,         (GtkItemFactoryCallback)hyphy_menu_item_callback,           HY_TREE_GTK_MENU_BASE+19,                 "<Item>" },
    { "/Tree/Branch Length Distributi_on",          NULL,         (GtkItemFactoryCallback)hyphy_menu_item_callback,           HY_TREE_GTK_MENU_BASE+20,                 "<Item>" },
    { "/Tree/sep8",         NULL,         NULL,           0,        "<Separator>" },
    { "/Tree/Tree Compar_ison",         NULL,         NULL,           0,                    "<Branch>" },
    { "/Tree/Tree Comparison/Test For _Equality",           NULL,         (GtkItemFactoryCallback)hyphy_menu_item_callback,           HY_TREE_GTK_MENU_BASE+35,                 "<Item>" },
    { "/Tree/Tree Comparison/Find _Subtree In Another Tree",            NULL,         (GtkItemFactoryCallback)hyphy_menu_item_callback,           HY_TREE_GTK_MENU_BASE+36,                 "<Item>" },
    { "/Tree/Tree Comparison/Find _Maximal Common Subtree",         NULL,         (GtkItemFactoryCallback)hyphy_menu_item_callback,           HY_TREE_GTK_MENU_BASE+37,                 "<Item>" },
    { "/Tree/Tree Comparison/Find Maximal Common _Forest",          NULL,         (GtkItemFactoryCallback)hyphy_menu_item_callback,           HY_TREE_GTK_MENU_BASE+38,                 "<Item>" },
    { "/Tree/Tree Comparison/_Match Tree To Pattern",           NULL,         (GtkItemFactoryCallback)hyphy_menu_item_callback,           HY_TREE_GTK_MENU_BASE+39,                 "<Item>" },
    { "/Tree/Match Leaves To Se_quence Data",           NULL,         (GtkItemFactoryCallback)hyphy_menu_item_callback,           HY_TREE_GTK_MENU_BASE+40,                 "<Item>" },
    { "/Tree/sep9",         NULL,         NULL,           0,        "<Separator>" },
    { "/Tree/_Additional Tools",            NULL,         NULL,           0,                    "<Branch>" }
};


//__________________________________________________________________

bool _HYTreePanel::_ProcessOSEvent (Ptr vEvent)
{
    static  char draggin = 0;

    _HY_GTK_UI_Message *theMessage = (_HY_GTK_UI_Message*)vEvent;

    if (theMessage->theEvent->type == GDK_BUTTON_PRESS) {
        GdkEventButton * be = (GdkEventButton*)theMessage->theEvent;
        if (be->button > 1) {
            if (FindClickedCell(be->x-windowContent->allocation.x,be->y-windowContent->allocation.y) == 0) {
                HandleContextPopup(be->x-windowContent->allocation.x,be->y-windowContent->allocation.y);
                return true;
            }
        }
    }

    if(!_HYTWindow::_ProcessOSEvent (vEvent)) {
        if (theMessage->theEvent->type == GDK_BUTTON_PRESS || theMessage->theEvent->type == GDK_2BUTTON_PRESS) {
            GdkEventButton * be = (GdkEventButton*)theMessage->theEvent;

            int ch = be->x-windowContent->allocation.x,
                cv = be->y-windowContent->allocation.y,
                c;

            c = FindClickedCell(ch,cv);
            if (c<0) {
                return false;
            }

            _HYComponent* thisComponent = (_HYComponent*)components(c);

            if (c==1) { // navBar
                ch -= componentL.lData[1];
                cv -= componentT.lData[1];
                if (navRect.Contains(ch,cv)) {
                    draggin = 1;
                } else {
                    SetNavRectCenter (ch,cv);
                }

                return true;
            } else if (c==0) {
                ch -= componentL.lData[0];
                cv -= componentT.lData[0];
                char shiftFlag = 0;

                if (be->state & GDK_SHIFT_MASK) {
                    shiftFlag |= 0x01;
                }
                if (be->state & GDK_CONTROL_MASK) {
                    shiftFlag |= 0x02;
                }

                if (shiftFlag > 2) {
                    draggin = 2;
                    //theEvent->iMsg = WM_MOUSEMOVE;
                } else {
                    if (IsVertical()) {
                        c = ch;
                        ch = cv+thisComponent->vOrigin;
                        cv = c+thisComponent->hOrigin;
                    } else {
                        ch+=thisComponent->hOrigin;
                        cv+=thisComponent->vOrigin;
                    }

                    if (theMessage->theEvent->type == GDK_2BUTTON_PRESS && currentSelection.lLength) {
                        InvokeNodeEditor();
                        return true;
                    }

                    if(FindSelection (ch,cv,shiftFlag)) {
                        _UpdateOperationsMenu();
                        RenderTree();
                    }
                }
            }
        }

        if (theMessage->theEvent->type == GDK_BUTTON_RELEASE) {
            draggin = 0;
            return true;
        }

        if (theMessage->theEvent->type == GDK_MOTION_NOTIFY) {
            GdkEventMotion * me = (GdkEventMotion*)theMessage->theEvent;

            if (draggin == 1) {
                int ch = me->x-componentL.lData[1]-windowContent->allocation.x,
                    cv = me->y-componentT.lData[1]-windowContent->allocation.y;

                SetNavRectCenter (ch,cv);
            } else {
                if (draggin == 2) {
                    _HYCanvas   *theTree    = (_HYCanvas*)GetObject (0);

                    int ch = me->x-windowContent->allocation.x,
                        cv = me->y-windowContent->allocation.y;

                    if (IsVertical()) {
                        cv = cv-componentL.lData[0]+theTree->hOrigin;
                        ch = ch-componentT.lData[0]+theTree->vOrigin;
                    } else {
                        ch = ch-componentL.lData[0]+theTree->hOrigin;
                        cv = cv-componentT.lData[0]+theTree->vOrigin;
                    }
                    if ((ch<theTree->_HYComponent::GetMaxW())&&(cv<theTree->_HYComponent::GetMaxH())) {
                        FishEyeProjection (ch,theTree->_HYComponent::GetMaxH()-cv,theTree->_HYComponent::GetMaxW(),
                                           theTree->_HYComponent::GetMaxH(),coordTree);
                        treeFlags |= HY_TREEPANEL_PROJECTION;
                        RenderTree(false);
                    }
                }
            }
            return true;
        } else {
            if (theMessage->theEvent->type == GDK_KEY_PRESS) {
                GdkEventKey * kp = (GdkEventKey*)theMessage->theEvent;

                if (kp->keyval == GDK_Delete || kp->keyval == GDK_BackSpace) {
                    DeleteCurrentSelection();
                    _UpdateOperationsMenu();
                    return true;
                } else if (kp->keyval == GDK_Return || kp->keyval == GDK_KP_Enter) {
                    InvokeNodeEditor ();
                    _UpdateOperationsMenu();
                    return true;
                }
            }
        }
        return false;
    }
    return true;
}
//__________________________________________________________________


bool        _HYTreePanel::_ProcessMenuSelection (long msel)
{
    if (_HYWindow::_ProcessMenuSelection(msel)) {
        return true;
    }

    switch (msel) {

    case HY_WINDOW_MENU_ID_FILE+1: { // Save tree
        _String           filePath,
                          dtreeName = treeName,
                          ext,
                          fileFormat ("File Format:");

        bool              good = false;

        _List             treeFormats;

        long pidOptions = GetUniversalSaveOptions (treeFormats);
        if (exportFormats.lLength == 0) {
            findGraphicsExporterComponents(exportFormats);
        }
        treeFormats << exportFormats;

        long              menuChoice = SaveFileWithPopUp (filePath, saveForTreesPrompt,dtreeName,
                                       fileFormat,
                                       treeFormats);
        if (menuChoice>=0) {
            if (menuChoice >= pidOptions) {
                _HYCanvas   *theTree = (_HYCanvas*)GetObject (0);

                bool    ruler = (!IsVertical())&&(scaleVariable.sLength);

                GdkPixbuf* theImage = gdk_pixbuf_get_from_drawable (NULL, theTree->thePane, NULL, 0, 0, 0, 0, -1, -1);
                if (ruler) {
                    _HYCanvas   *theRuler = (_HYCanvas*)GetObject (2);
                    GdkPixbuf* theImage2 = gdk_pixbuf_get_from_drawable (NULL, theRuler->thePane, NULL, 0, 0, 0, 0, -1, -1);
                    long        w  = gdk_pixbuf_get_width(theImage),
                                h1 = gdk_pixbuf_get_height (theImage),
                                h2 = gdk_pixbuf_get_height(theImage2);
                    GdkPixbuf* composite = gdk_pixbuf_new (GDK_COLORSPACE_RGB,false,gdk_pixbuf_get_bits_per_sample(theImage),w,h1+h2);
                    gdk_pixbuf_copy_area (theImage,0,0,w,h1,composite,0,h2);
                    g_object_unref (theImage);
                    gdk_pixbuf_copy_area (theImage2,0,0,w,h2,composite,0,0);
                    g_object_unref (theImage2);
                    theImage = composite;
                }
                gdk_pixbuf_save (theImage, filePath.sData,((_String*)treeFormats (menuChoice))->sData,NULL,NULL);
                g_object_unref (theImage);
            } else {
                HandleTreeSave (menuChoice, filePath);
            }
        }
        return true;
    }
    break;
    case HY_WINDOW_MENU_ID_FILE+2: {
        _PrintTree();
        break;
    }

    case HY_WINDOW_MENU_ID_EDIT:   // Undo
        UndoLastOperation();
        break;

    case HY_WINDOW_MENU_ID_EDIT+1: { // COPY
        if (treeFlags&HY_TREEPANEL_CLIPBOARD_READY) {
            CutSelectionToClipboard (false);
        } else {
            /*_HYCanvas   *theTree = (_HYCanvas*)GetObject (0);
            bool ruler = (!IsVertical())&&(scaleVariable.sLength);
            if (ruler)
            {
                _HYCanvas   *theRuler = (_HYCanvas*)GetObject (2);

            }
            else
            {
                PlaceBitmapInClipboard ((HBITMAP)GetCurrentObject (theTree->thePane,OBJ_BITMAP),theWindow);
            }*/
        }
        break;
    }

    case HY_WINDOW_MENU_ID_EDIT+2:   // Cut
        CutSelectionToClipboard ();
        break;

    case HY_WINDOW_MENU_ID_EDIT+3:   // Paste
        PasteClipboardTree();
        break;

    case HY_WINDOW_MENU_ID_EDIT+4:   // Delete
        DeleteCurrentSelection();
        break;

    case HY_WINDOW_MENU_ID_EDIT+5:   // Select All
        SelectAllBranches();
        break;

    case HY_WINDOW_MENU_ID_EDIT+6:   // S & R
        HandleSearchAndReplace(false);
        break;

    case HY_WINDOW_MENU_ID_EDIT+7:   // S & R
        HandleSearchAndReplace(true);
        break;

    case HY_TREE_GTK_MENU_BASE: {
        GtkCheckMenuItem * menuItem = GTK_CHECK_MENU_ITEM(gtk_item_factory_get_item_by_action(menu_items,HY_TREE_GTK_MENU_BASE));
        if ( ((bool)(treeFlags&HY_TREEPANEL_TIP_LABELS)) != gtk_check_menu_item_get_active (menuItem)) {
            unsigned short newF;
            if (treeFlags&HY_TREEPANEL_TIP_LABELS) {
                newF = treeFlags - HY_TREEPANEL_TIP_LABELS;
            } else {
                newF = treeFlags + HY_TREEPANEL_TIP_LABELS;
            }
            SetFlags (newF);
        }
        break;
    }

    case HY_TREE_GTK_MENU_BASE+1: {
        GtkCheckMenuItem * menuItem = GTK_CHECK_MENU_ITEM(gtk_item_factory_get_item_by_action(menu_items,HY_TREE_GTK_MENU_BASE+1));
        if ( ((bool)(treeFlags&HY_TREEPANEL_INT_LABELS)) != gtk_check_menu_item_get_active (menuItem)) {
            unsigned short newF;
            if (treeFlags&HY_TREEPANEL_INT_LABELS) {
                newF = treeFlags - HY_TREEPANEL_INT_LABELS;
            } else {
                newF = treeFlags + HY_TREEPANEL_INT_LABELS;
            }
            SetFlags (newF);
        }
        break;
    }

    case HY_TREE_GTK_MENU_BASE+2: {
        SwapSelectedSubTrees ();
        break;
    }

    case HY_TREE_GTK_MENU_BASE+3: {
        CollapseSelectedBranch ();
        break;
    }

    case HY_TREE_GTK_MENU_BASE+4: {
        JoinSelectedBranches ();
        break;
    }

    case HY_TREE_GTK_MENU_BASE+5: {
        GraftATip ();
        break;
    }

    case HY_TREE_GTK_MENU_BASE+6: {
        RerootTree ();
        break;
    }

    case HY_TREE_GTK_MENU_BASE+7: {
        FlipSelectedBranches ();
        break;
    }

    case HY_TREE_GTK_MENU_BASE+8:
    case HY_TREE_GTK_MENU_BASE+9:
    case HY_TREE_GTK_MENU_BASE+10: {
        HandleSelection (msel-HY_TREE_GTK_MENU_BASE-8);
        break;
    }

    case HY_TREE_GTK_MENU_BASE+17:
    case HY_TREE_GTK_MENU_BASE+18: {
        HandleSelection (msel-HY_TREE_GTK_MENU_BASE-14);
        break;
    }

    case HY_TREE_GTK_MENU_BASE+30:
    case HY_TREE_GTK_MENU_BASE+31:
    case HY_TREE_GTK_MENU_BASE+32:
    case HY_TREE_GTK_MENU_BASE+33: {
        HandleSelection (msel-HY_TREE_GTK_MENU_BASE-25);
        break;
    }

    case HY_TREE_GTK_MENU_BASE+11: {
        InvokeNodeEditor ();
        break;
    }

    case HY_TREE_GTK_MENU_BASE+12: {
        RecalculateLikelihood ();
        break;
    }

    case HY_TREE_GTK_MENU_BASE+13: {
        DisplayParameterTable();
        break;
    }

    case HY_TREE_GTK_MENU_BASE+14: {
        HandleViewOptions ();
        break;
    }

    case HY_TREE_GTK_MENU_BASE+15:
    case HY_TREE_GTK_MENU_BASE+16: {
        ShowModelMatrix (msel-HY_TREE_GTK_MENU_BASE-15);
        break;
    }

    case HY_TREE_GTK_MENU_BASE+19:
    case HY_TREE_GTK_MENU_BASE+20: {
        GenerateDistanceTable (msel-HY_TREE_GTK_MENU_BASE-19);
        break;
    }

    case HY_TREE_GTK_MENU_BASE+24:
    case HY_TREE_GTK_MENU_BASE+25: {
        HandleLabels(msel-HY_TREE_GTK_MENU_BASE-24);
        break;
    }

    case HY_TREE_GTK_MENU_BASE+35:
    case HY_TREE_GTK_MENU_BASE+36:
    case HY_TREE_GTK_MENU_BASE+37:
    case HY_TREE_GTK_MENU_BASE+38:
    case HY_TREE_GTK_MENU_BASE+39: {
        HandleComparison (msel-HY_TREE_GTK_MENU_BASE-35);
        _UpdateOperationsMenu();
        break;
    }

    case HY_TREE_GTK_MENU_BASE+40: {
        MatchToDataSet();
        break;
    }

    default: { // proc menu
        if (msel>=HY_TREE_GTK_MENU_BASE+1000) {
            ExecuteProcessor (msel-HY_TREE_GTK_MENU_BASE-1000);
            _UpdateOperationsMenu();
            return true;
        }
    }
    }

    if (((msel>=HY_TREE_GTK_MENU_BASE)&&(msel<HY_TREE_GTK_MENU_BASE+12))||
            ((msel>=HY_WINDOW_MENU_ID_EDIT)&&(msel<HY_WINDOW_MENU_ID_EDIT+6))) {
        _UpdateOperationsMenu();
    }

    return true;
}



//__________________________________________________________________

void        _HYTreePanel::_PaintNavRect(void)
{
    if (GTK_WIDGET_MAPPED (theWindow)) {
        navRect = ComputeNavRect();
        GdkGC* navDC = gdk_gc_new (theWindow->window);
        gdk_gc_set_line_attributes(navDC,2,GDK_LINE_SOLID,GDK_CAP_NOT_LAST,GDK_JOIN_MITER);
        gdk_gc_set_foreground(navDC,&navPen);

        _HYCanvas* theCanvas = (_HYCanvas*)GetCellObject (0,0);

        _HYRect r;

        r.left   = navRect.left+theCanvas->rel.left+HY_TREEPANEL_NAVSPACING;
        r.right  = navRect.right+theCanvas->rel.left;
        r.top    = navRect.top+theCanvas->rel.top+HY_TREEPANEL_NAVSPACING;
        r.bottom = navRect.bottom+theCanvas->rel.top;

        gdk_draw_rectangle(theWindow->window,navDC,false,
                           r.left + windowContent->allocation.x,
                           r.top  + windowContent->allocation.y,
                           r.right  - r.left + 1,
                           r.bottom - r.top + 1);

        _PaintLFStatus (nil);
        g_object_unref (navDC);
    }
}

//__________________________________________________________________

void        _HYTreePanel::_PrintTree(long hPages, long vPages)
{
    _String TBI ("Tree Printing Has Not Yet Been Implemented\n");
    ReportWarning (TBI);
}

//__________________________________________________________________

void _HYTreePanel::_SetMenuBar(void)
{
    if (menu_items && !gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Edit/Search and Replace")) {
        gtk_item_factory_create_items (menu_items,  sizeof (hyphy_tree_window_menu) / sizeof (hyphy_tree_window_menu[0]), hyphy_tree_window_menu, this);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Tree/Tree Comparison/Find Subtree In Another Tree"),false);

        gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Tree/Tip Labels")),
                                       treeFlags&HY_TREEPANEL_TIP_LABELS);
        gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Tree/Internal Labels")),
                                       treeFlags&HY_TREEPANEL_INT_LABELS);

        gtk_widget_set_sensitive(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Tree/Swap Subtrees"),false);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Tree/Collapse Branch"),false);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Tree/Join"),false);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Tree/Graft A Tip"),false);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Tree/Reroot"),false);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Tree/Flip Tip Ordering"),false);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Tree/Select Branches/Select Entire Subtree"),false);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Tree/Select Branches/Grow Selection"),false);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Tree/Edit Properties"),false);

        gtk_widget_set_sensitive(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Edit/Select All"),true);

        if (likeFuncID < 0) {
            gtk_widget_set_sensitive(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Tree/Optimize Again"),false);
            gtk_widget_set_sensitive(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Tree/Show Parameters in Table"),false);
        }

        gtk_widget_set_sensitive(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Tree/Show Rate Matrix"),false);
        gtk_widget_set_sensitive(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Tree/Show Transition Matrix"),false);

        if (treeProcessors.lLength == 0) {
            gtk_widget_set_sensitive(gtk_item_factory_get_widget(menu_items,"<HY_WINDOW>/Tree/Additional Tools"),false);
        } else {
            for (long k=0; k<treeProcessors.lLength; k++) {
                _String *thisItem = (_String*)treeProcessors (k),
                         chopped = thisItem->Cut (thisItem->FindBackwards ('/',0,-1)+1,-1),
                         type = "<Item>";

                GtkItemFactoryEntry aProcEntry = {NULL,NULL,(GtkItemFactoryCallback)hyphy_menu_item_callback,HY_TREE_GTK_MENU_BASE+1000+k,type.sData};
                chopped = _String("/Tree/Additional Tools/")&chopped;
                aProcEntry.path = chopped.sData;

                gtk_item_factory_create_items (menu_items,  1, &aProcEntry, this);
                //InsertMenu        (procMenu, 0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_TREE_GTK_MENU_BASE+1000+k, chopped.sData);
            }
        }
    }
}

//__________________________________________________________________

void _HYTreePanel::_UnsetMenuBar(void)
{

}

//__________________________________________________________________

void _HYTreePanel::_UpdateOperationsMenu (void)
{
    node<nodeCoord>* node1, *node2, *t;

    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+2),false);
    //EnableMenuItem(treeMenu, 3 ,MF_GRAYED|MF_BYPOSITION);
    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+3),false);
    //EnableMenuItem(treeMenu, 4 ,MF_GRAYED|MF_BYPOSITION);
    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+4),false);
    //EnableMenuItem(treeMenu, 5 ,MF_GRAYED|MF_BYPOSITION);
    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+5),false);
    //EnableMenuItem(treeMenu, 6 ,MF_GRAYED|MF_BYPOSITION);
    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+6),false);
    //EnableMenuItem(treeMenu, 7 ,MF_GRAYED|MF_BYPOSITION);
    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+7),false);
    //EnableMenuItem(treeMenu, 8 ,MF_GRAYED|MF_BYPOSITION);
    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+8),false);
    //EnableMenuItem(selectMenu,0 ,MF_GRAYED|MF_BYPOSITION);
    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+30),false);
    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+32),false);
    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+33),true);
    //EnableMenuItem(selectMenu,5 ,MF_GRAYED|MF_BYPOSITION);
    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+11),false);
    //EnableMenuItem(treeMenu,11 ,MF_GRAYED|MF_BYPOSITION);
    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+15),false);
    //EnableMenuItem(treeMenu,19 ,MF_GRAYED|MF_BYPOSITION);
    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+16),false);
    //EnableMenuItem(treeMenu,20 ,MF_GRAYED|MF_BYPOSITION);

    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_WINDOW_MENU_ID_EDIT+4),false);
    //EnableMenuItem(editMenu, 8 ,MF_GRAYED|MF_BYPOSITION);
    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_WINDOW_MENU_ID_EDIT+7),false);
    //EnableMenuItem(editMenu, 11,MF_GRAYED|MF_BYPOSITION);
    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+36),false);
    //EnableMenuItem(compMenu, 1 ,MF_GRAYED|MF_BYPOSITION);

    bool  good = true;
    long  k,j;
    if (currentSelection.lLength==2) {
        node1 = (node<nodeCoord>*)currentSelection(0);
        node2 = (node<nodeCoord>*)currentSelection(1);
        t = node1->parent;
        while (t) {
            if (t==node2) {
                good = false;
                break;
            }
            t = t->parent;
        }
        if (good) {
            t = node2->parent;
            while (t) {
                if (t==node1) {
                    good = false;
                    break;
                }
                t = t->parent;
            }
        }
        if (good) {
            gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+2),true);
        }
    }
    if (currentSelection.lLength) {
        gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+5),true);
        //EnableMenuItem(treeMenu, 6 ,MF_ENABLED|MF_BYPOSITION);
        gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+11),true);
        //EnableMenuItem(treeMenu, 11 ,MF_ENABLED|MF_BYPOSITION);
        gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+30),true);
        gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+31),true);
        gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+32),true);
        gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+33),true);
        //EnableMenuItem(selectMenu,5 ,MF_ENABLED|MF_BYPOSITION);
        gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_WINDOW_MENU_ID_EDIT+7),true);
        //EnableMenuItem(editMenu, 11,MF_ENABLED|MF_BYPOSITION);
        for (k=0; k<currentSelection.lLength; k++) {
            node1 = (node<nodeCoord>*)currentSelection(k);
            if (node1->get_num_nodes()) {
                gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+3),true);
                //EnableMenuItem(treeMenu, 4 ,MF_ENABLED|MF_BYPOSITION);
                break;
            }
        }
        for (k=0; k<currentSelection.lLength; k++) {
            node1 = (node<nodeCoord>*)currentSelection(k);
            if (!node1->get_num_nodes()) {
                gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+4),true);
                //EnableMenuItem(treeMenu, 5 ,MF_ENABLED|MF_BYPOSITION);
                break;
            }
        }
        if (currentSelection.lLength>=2) {
            node1 = (node<nodeCoord>*)currentSelection(0);
            t = node1->parent;
            if (t&&(t->get_num_nodes()>currentSelection.lLength)) {
                for (k=1; k<currentSelection.lLength; k++) {
                    node1 = (node<nodeCoord>*)currentSelection(k);
                    if (node1->parent!=t) {
                        break;
                    }
                }
                if (k==currentSelection.lLength) {
                    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+4),true);
                    //EnableMenuItem(treeMenu, 5 ,MF_ENABLED|MF_BYPOSITION);
                }
            }
        } else {
            node1 = (node<nodeCoord>*)currentSelection(0);
            if (node1->parent) {
                GtkWidget * rerootItem = gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+6);
                gtk_widget_set_sensitive (rerootItem,true);
                gtk_label_set_text (GTK_LABEL (gtk_bin_get_child(GTK_BIN(rerootItem))), "Reroot");
                //EnableMenuItem(treeMenu, 7 ,MF_ENABLED|MF_BYPOSITION);
                //ModifyMenu    (treeMenu, 7, MF_BYPOSITION, HY_TREE_GTK_MENU_BASE+6, "Reroot");
            }
            if (node1->get_num_nodes()>0) {
                gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+7),true);
                //EnableMenuItem(treeMenu, 8 ,MF_ENABLED|MF_BYPOSITION);
                gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+8),true);
                //EnableMenuItem(selectMenu,0 ,MF_ENABLED|MF_BYPOSITION);
            }
            if (node1->in_object.varRef>=0) {
                _CalcNode* thisCNode = (_CalcNode*)LocateVar(node1->in_object.varRef);

                if (thisCNode&&(thisCNode->GetModelIndex()>=0)) {
                    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+15),true);
                    //EnableMenuItem(treeMenu,19 ,MF_ENABLED|MF_BYPOSITION);
                    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+16),true);
                    //EnableMenuItem(treeMenu,20 ,MF_ENABLED|MF_BYPOSITION);
                }
            }
        }
    } else {
        _TheTree *me = LocateMyTreeVariable();
        if (me) {
            GtkWidget * rerootItem = gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+6);
            gtk_widget_set_sensitive (rerootItem,true);
            if (me->RootedFlag()==UNROOTED) {
                gtk_label_set_text (GTK_LABEL (gtk_bin_get_child(GTK_BIN(rerootItem))), "Balance");
            } else {
                gtk_label_set_text (GTK_LABEL (gtk_bin_get_child(GTK_BIN(rerootItem))), "Unroot");
            }
        }
    }

    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+12),likeFuncID>=0);
    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+13),likeFuncID>=0);


    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_WINDOW_MENU_ID_EDIT+3),treePanelClipboardRoot&&(currentSelection.lLength==1));
    gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_WINDOW_MENU_ID_EDIT+2),false);

    // check if can cut/paste
    t = nil;

    for (k=0; k<currentSelection.lLength; k++) {
        node1 = (node<nodeCoord>*)currentSelection.lData[k];
        if (node1->parent) {
            if (currentSelection.Find((long)node1->parent)<0) {
                if (t) {
                    break;
                } else {
                    t = node1;
                }
            }
            for (j=0; j<node1->nodes.length; j++)
                if (currentSelection.Find((long)node1->nodes.data[j])<0) {
                    break;
                }

            if (j<node1->nodes.length) {
                break;
            }
        } else {
            if (t) {
                break;
            } else {
                t = node1;
            }
        }
    }
    selectionTop = nil;
    treeFlags &= 0xFF7F;
    if (t&&(t->parent!=coordTree)&&(t->parent)) {
        if (k==currentSelection.lLength) {
            gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_WINDOW_MENU_ID_EDIT+2),true);
            //EnableMenuItem(editMenu, 3 ,MF_ENABLED|MF_BYPOSITION);
            gtk_widget_set_sensitive (gtk_item_factory_get_widget_by_action(menu_items,HY_TREE_GTK_MENU_BASE+36),true);
            //EnableMenuItem(compMenu, 1 ,MF_ENABLED|MF_BYPOSITION);
            selectionTop = t;
            treeFlags |= HY_TREEPANEL_CLIPBOARD_READY;
        }
    }
    _String undoMessage;
    GtkWidget * undoItem = gtk_item_factory_get_widget_by_action(menu_items,HY_WINDOW_MENU_ID_EDIT);
    //gtk_widget_set_sensitive (rerootItem,true);
    //EnableMenuItem(editMenu, 0 ,MF_ENABLED|MF_BYPOSITION);
    switch (undoCode) {
    case 1:
        undoMessage = "Undo Swap";
        break;
    case 2:
        undoMessage = "Undo Flip";
        break;
    case 3:
        undoMessage = "Undo Collapse";
        break;
    case 4:
        undoMessage = "Undo Delete";
        break;
    case 5:
        undoMessage = "Undo Join";
        break;
    case 6:
        undoMessage = "Undo Cut";
        break;
    case 7:
        undoMessage = "Undo Graft";
        break;
    case 8:
        undoMessage = "Undo Paste";
        break;
    case 9:
        undoMessage = "Undo Subtree Move";
        break;
    default:
        undoMessage = "Can't Undo";
    }
    gtk_widget_set_sensitive (undoItem,!(undoMessage==_String("Can't Undo")));
    gtk_label_set_text (GTK_LABEL (gtk_bin_get_child(GTK_BIN(undoItem))), undoMessage.sData);
}
//__________________________________________________________________

void _HYTreePanel::_HandleIdleEvent (void)
{
    /*#ifdef TARGET_API_MAC_CARBON
    Point    curMouse;
    GetGlobalMouse (&curMouse);

    unsigned long t;
    GetDateTime(&t);


    if ((abs(curMouse.h-saveMouseH)<=3)
      &&(abs(curMouse.v-saveMouseV)<=3)
      &&(t-lastSave>.5))

    {
        if (!HasToolTip())
        {
            GrafPtr curPort;
            GetPort (&curPort);
            SetPort (GetWindowPort (theWindow));
            _DisplayBranchFloater();
            SetPort (curPort);
        }

        lastSave   = t;
    }

    saveMouseH = curMouse.h;
    saveMouseV = curMouse.v;
    #endif*/
}

//__________________________________________________________________

void        _HYTreePanel::_PaintLFStatus(Ptr p)
{
    if (likeFuncID<0) {
        _PaintTheCircle (redButtonIcon,theWindow);
    } else {
        if (dubiousNodes.lLength) {
            _PaintTheCircle (yellowButtonIcon,theWindow);
        } else {
            _PaintTheCircle (greenButtonIcon,theWindow);
        }
    }
}

//EOF
