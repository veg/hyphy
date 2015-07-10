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
#include "HYTextBox.h"
#include "HYUtils.h"
#include "HYEventTypes.h"
#include "HYGraphicPane.h"
#include "HYPlatformWindow.h"
#include <gdk/gdkkeysyms.h>

void       ApplyTagToAll                (GtkTextBuffer*, GtkTextTag*);
gboolean   hy_textbox_key_interceptor   (GtkWidget *, GdkEventKey*, gpointer);

_HYFont                     defaultTBFont       =  {_HY_SANS_FONT,10,HY_FONT_PLAIN};

//__________________________________________________________________

gboolean   hy_textbox_change_notifier (GtkTextBuffer *te,gpointer tbxp)
{
    _HYTextBox * tbx = (_HYTextBox*)tbxp;

    if (tbx->messageRecipient) {
        tbx->messageRecipient->ProcessEvent (generateTextEditChangeEvent (tbx->GetID(),1));
    }

    return false;
}


//__________________________________________________________________

gboolean   hy_textbox_selchange_notifier (GtkTextBuffer *tbuf,GtkTextIter *arg1,GtkTextMark *arg2,gpointer tbxp)
{
    _HYTextBox * tbx = (_HYTextBox*)tbxp;

    if (tbx->messageRecipient) {
        char * markName = (char*)gtk_text_mark_get_name(arg2);
        if (markName) {
            _String markname (markName);
            if (markname == _String("selection_bound")) { //TBI - will miss many changes
                //printf ("Changed Selection \n");
                tbx->messageRecipient->ProcessEvent (generateTextEditChangeEvent (tbx->GetID(),0));
            }
        }
    }

    return false;
}

//__________________________________________________________________

gboolean   hy_textbox_key_interceptor (GtkWidget *te, GdkEventKey *kp, gpointer tbxp)
{
    _HYTextBox * tbx = (_HYTextBox*)tbxp;

    if (tbx->messageRecipient && (tbx->boxFlags & HY_TB_ARROWS)) {
        if (kp->keyval==GDK_Down || kp->keyval == GDK_KP_Down) {
            tbx->messageRecipient->ProcessEvent (generateTextEditChangeEvent (tbx->GetID(),3));
            return true;
        } else if (kp->keyval==GDK_Up || kp->keyval == GDK_KP_Up) {
            tbx->messageRecipient->ProcessEvent (generateTextEditChangeEvent (tbx->GetID(),4));
            return true;
        }
    }

    if ((tbx->boxFlags & HY_TB_BIGBOX) == 0) {
        if (kp->keyval==GDK_KP_Enter || kp->keyval==GDK_Return || kp->keyval==GDK_Escape) {
            if (tbx->messageRecipient) {
                tbx->messageRecipient->ProcessEvent (generateTextEditChangeEvent (tbx->GetID(),2));
            }

            return true;
        } else if (kp->keyval==GDK_Tab) {
            return true;
        }
    }

    return FALSE;
}

//__________________________________________________________________

void ApplyTagToAll (GtkTextBuffer* tbuf, GtkTextTag* tt)
{
    GtkTextIter               startIT,
                              endIT;

    gtk_text_buffer_get_iter_at_offset (tbuf, &startIT, 0);
    gtk_text_buffer_get_iter_at_offset (tbuf, &endIT, gtk_text_buffer_get_char_count(tbuf));
    gtk_text_buffer_apply_tag (tbuf,tt,&startIT,&endIT);
}

//__________________________________________________________________

_HYPlatformTextBox::_HYPlatformTextBox  (void)
{
    backFill   = backFill = HYColorToGDKColor((_HYColor) {
        255,255,255
    });
    pLabelFont = pango_font_description_new();

    textBoxRect.x       = textBoxRect.y     = 0;
    textBoxRect.width   = textBoxRect.height = 100;
    textColor           = HYColorToGDKColor((_HYColor) {
        0,0,0
    });
    te                  = nil;

    _HYTextBox * theParent = (_HYTextBox*)this;
    isSingleLine        = true;

    te                  = gtk_text_view_new  ();
    gtk_text_view_set_wrap_mode (GTK_TEXT_VIEW (te), GTK_WRAP_NONE);
    scrollWindow        = gtk_scrolled_window_new (NULL,NULL);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW(scrollWindow),GTK_POLICY_NEVER,GTK_POLICY_NEVER);

    gtk_container_add(GTK_CONTAINER(scrollWindow),te);
    gtk_container_add (GTK_CONTAINER(theParent->parentWindow), scrollWindow);
    gtk_widget_set_app_paintable(te,true);
    gtk_widget_show   (te);
    gtk_widget_show   (scrollWindow);
    isSingleLine = true;
    g_signal_connect (G_OBJECT (te), "key-press-event", G_CALLBACK (hy_textbox_key_interceptor), theParent);
    g_signal_connect (G_OBJECT (gtk_text_view_get_buffer(GTK_TEXT_VIEW(te))), "changed", G_CALLBACK (hy_textbox_change_notifier), theParent);
    g_signal_connect (G_OBJECT (gtk_text_view_get_buffer(GTK_TEXT_VIEW(te))), "mark-set", G_CALLBACK (hy_textbox_selchange_notifier), theParent);
}

//__________________________________________________________________

_HYPlatformTextBox::~_HYPlatformTextBox (void)
{
    pango_font_description_free(pLabelFont);
}

//__________________________________________________________________

void    _HYPlatformTextBox::_SetBackColor (_HYColor& c)
{
    backFill = HYColorToGDKColor(c);
}

//__________________________________________________________________

void    _HYPlatformTextBox::_SetBackTColor (_HYColor& c)
{
    backTFill = HYColorToGDKColor(c);
    GtkTextBuffer * tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(te));
    gboolean        trueVal = TRUE;
    GtkTextTag * ct = gtk_text_tag_table_lookup(gtk_text_buffer_get_tag_table(tbuf),"hyphy_texbox_bgtag");
    if (!ct) {
        ct = gtk_text_buffer_create_tag (tbuf,"hyphy_texbox_bgtag","background-gdk",&backTFill,"background-set",&trueVal,NULL);
        ApplyTagToAll (tbuf,ct);
    } else {
        g_object_set (ct, "background-gdk",&backTFill,"background-set",&trueVal,NULL);
    }
}


//__________________________________________________________________

void    _HYPlatformTextBox::_SetForeColor (_HYColor& c)
{
    textColor = HYColorToGDKColor(c);
    GtkTextBuffer * tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(te));
    gboolean        trueVal = TRUE;
    GtkTextTag * ct = gtk_text_tag_table_lookup(gtk_text_buffer_get_tag_table(tbuf),"hyphy_texbox_fgtag");
    if (!ct) {
        ct = gtk_text_buffer_create_tag (tbuf,"hyphy_texbox_fgtag","foreground-gdk",&textColor,"foreground-set",&trueVal,NULL);
        ApplyTagToAll (tbuf,ct);
    } else {
        g_object_set (ct, "foreground-gdk",&textColor,"foreground-set",&trueVal,NULL);
    }
}


//__________________________________________________________________

void    _HYPlatformTextBox::_SetFont (_HYFont& f)
{
    HYFont2PangoFontDesc(f,pLabelFont);
    gtk_widget_modify_font (te, pLabelFont);
}

//__________________________________________________________________
void        _HYPlatformTextBox::_Update (Ptr p)
{
    _Paint (p);
}


//__________________________________________________________________
void        _HYPlatformTextBox::_SetDimensions (_HYRect r, _HYRect rel)
{
    _HYTextBox* theParent = (_HYTextBox *) this;
    theParent->_HYPlatformComponent::_SetDimensions (r,rel);
    _SetVisibleSize (rel);
}


//__________________________________________________________________
void        _HYPlatformTextBox::_SetVisibleSize (_HYRect rel)
{
    _HYTextBox *theParent = (_HYTextBox*) this;

    textBoxRect.y       = rel.top   + theParent->margins.top;
    textBoxRect.x       = rel.left  + theParent->margins.left;
    textBoxRect.height  = rel.bottom- theParent->margins.bottom - textBoxRect.y + 1;
    textBoxRect.width   = rel.right - theParent->margins.right  - textBoxRect.x + 1;

    if (te) {
        bool          newSL = !_NeedMultiLines ();
        if (newSL!=isSingleLine) {
            if (((_HYTextBox*)this)->boxFlags & HY_TB_BIGBOX) {
                if (!newSL) {
                    gtk_text_view_set_wrap_mode (GTK_TEXT_VIEW (te),GTK_WRAP_CHAR);
                    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW(scrollWindow),GTK_POLICY_NEVER,GTK_POLICY_ALWAYS);
                }
            } else {
                gtk_text_view_set_wrap_mode (GTK_TEXT_VIEW (te), newSL?GTK_WRAP_NONE:GTK_WRAP_CHAR);
            }
            isSingleLine = newSL;
        }
        gtk_fixed_move (GTK_FIXED(theParent->parentWindow), scrollWindow, textBoxRect.x, textBoxRect.y);
        gtk_widget_set_size_request(scrollWindow,textBoxRect.width,textBoxRect.height);
    }
}


//__________________________________________________________________
void        _HYPlatformTextBox::_SetAlignFlags (unsigned char f)
{
    if (te) {
        gtk_text_view_set_justification (GTK_TEXT_VIEW(te), f==HY_ALIGN_LEFT?GTK_JUSTIFY_LEFT:(f==HY_ALIGN_RIGHT?GTK_JUSTIFY_RIGHT:GTK_JUSTIFY_CENTER));
    }
}


//__________________________________________________________________
_String     _HYPlatformTextBox::_GetText (void)
{
    GtkTextIter               startIT,
                              endIT;

    GtkTextBuffer             *tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(te));
    gtk_text_buffer_get_iter_at_offset (tbuf, &startIT, 0);
    gtk_text_buffer_get_iter_at_offset (tbuf, &endIT, gtk_text_buffer_get_char_count(tbuf));
    gsize                     br, bw;
    gchar                     *bufText  = gtk_text_buffer_get_text (tbuf,&startIT,&endIT,true),
                               *asciiText = g_locale_from_utf8 (bufText, -1, &br, &bw, NULL);

    _String                  *res = new _String ((char*)asciiText);
    g_free                   (asciiText);
    g_free                   (bufText);
    return res;
}


//__________________________________________________________________
void    _HYPlatformTextBox::_StoreText (_String*& rec, bool selOnly)
{
    GtkTextIter               startIT,
                              endIT;

    GtkTextBuffer             *tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(te));
    if (selOnly) {
        if (!gtk_text_buffer_get_selection_bounds (tbuf, &startIT, &endIT)) {
            rec = new _String;
            return;
        }
    } else {
        gtk_text_buffer_get_iter_at_offset (tbuf, &startIT, 0);
        gtk_text_buffer_get_iter_at_offset (tbuf, &endIT, gtk_text_buffer_get_char_count(tbuf));
    }
    gsize                     br, bw;
    gchar                     *bufText  = gtk_text_buffer_get_text (tbuf,&startIT,&endIT,true),
                               *asciiText = g_locale_from_utf8 (bufText, -1, &br, &bw, NULL);

    rec                      = new _String ((char*)asciiText);
    g_free                   (asciiText);
    g_free                   (bufText);
}


//__________________________________________________________________
bool        _HYPlatformTextBox::_NeedMultiLines (void)
{
    if (((_HYTextBox*)this)->boxFlags & HY_TB_BIGBOX) {
        return 1;
    } else {
        return (textBoxRect.height>((_HYTextBox*)this)->editBoxFont.size * 2.5);
    }
}


//__________________________________________________________________
void        _HYPlatformTextBox::_SetText (const _String& editBoxText)
{
    _HYTextBox * theParent = (_HYTextBox*)this;

    _HYGuiObject* stashRec = theParent->messageRecipient;
    theParent->SetMessageRecipient (nil);
    if (editBoxText.sData) {
        gtk_text_buffer_set_text (gtk_text_view_get_buffer(GTK_TEXT_VIEW(te)),editBoxText.sData,editBoxText.sLength);
    } else {
        gtk_text_buffer_set_text (gtk_text_view_get_buffer(GTK_TEXT_VIEW(te)),empty.sData,editBoxText.sLength);
    }

    theParent->SetMessageRecipient (stashRec);
}


//__________________________________________________________________
void        _HYPlatformTextBox::_InsertText (const _String& editBoxText, bool append)
{
    if (!te) {
        _SetText (editBoxText);
    } else {
        GtkTextIter               startIT,
                                  endIT;
        GtkTextBuffer             *tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(te));

        if (append) {
            gtk_text_buffer_get_selection_bounds (tbuf,&startIT,&endIT);
            gtk_text_iter_set_offset (&endIT,gtk_text_buffer_get_char_count(tbuf));
        } else {
            if (gtk_text_buffer_get_selection_bounds (tbuf,&startIT,&endIT)) {
                gtk_text_buffer_delete (tbuf, &startIT, &endIT);
            } else {
                gtk_text_buffer_insert_at_cursor (tbuf, editBoxText.sData, editBoxText.sLength);
                return;
            }
        }
        gtk_text_buffer_insert (tbuf, &endIT, editBoxText.sData, editBoxText.sLength);
        if (append) {
            gtk_text_iter_set_offset (&endIT,gtk_text_buffer_get_char_count(tbuf));
            GtkTextMark * endMark = gtk_text_buffer_create_mark(tbuf,NULL,&endIT,false);
            gtk_text_view_scroll_mark_onscreen (GTK_TEXT_VIEW(te), endMark);
            gtk_text_buffer_delete_mark (tbuf, endMark);
        }
    }
}

//__________________________________________________________________

void        _HYPlatformTextBox::_EnableTextBox (bool e)
{
    if (te) {
        gtk_widget_set_sensitive(te,e);
    }
}

//__________________________________________________________________

void        _HYPlatformTextBox::_FocusComponent (void)
{
    if (te) {
        gtk_widget_grab_focus (te);
    }
}

//__________________________________________________________________

void        _HYPlatformTextBox::_UnfocusComponent (void)
{
}

//__________________________________________________________________

void        _HYTextBox::_Activate (void)
{
    _HYPlatformComponent::_Activate();
    gtk_widget_show(te);
    //TBI
}

//__________________________________________________________________

void        _HYTextBox::_IdleHandler (void)
{
}

//__________________________________________________________________

void        _HYTextBox::_Deactivate (void)
{
    if (te) {
        gtk_widget_hide(te);
    }
}

//__________________________________________________________________

void    _HYTextBox::_SetSelection (long s, long e)
{
    if (te && s>=0 && e>=s) {
        GtkTextIter               startIT;
        GtkTextBuffer             *tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(te));

        gtk_text_buffer_get_iter_at_offset (tbuf, &startIT, s);
        gtk_text_buffer_move_mark_by_name (tbuf, "selection_bound" ,&startIT);
        long        cCount = gtk_text_buffer_get_char_count (tbuf);
        if (e<cCount) {
            e=e+1;
        } else {
            e = cCount;
        }
        gtk_text_buffer_get_iter_at_offset (tbuf, &startIT, e);
        gtk_text_buffer_move_mark_by_name (tbuf, "insert" ,&startIT);
    }
}



//__________________________________________________________________

void    _HYTextBox::_MarkForUpdate (void)
{
    if (te) {
        gtk_widget_queue_draw_area (te,0,0,textBoxRect.width,textBoxRect.height);
    }

    _HYPlatformComponent::_MarkForUpdate();
}


//__________________________________________________________________

bool    _HYTextBox::_IsEmpty (void)
{
    if (te) {
        return !gtk_text_buffer_get_char_count(gtk_text_view_get_buffer(GTK_TEXT_VIEW(te)));
    }
    return true;
}


//__________________________________________________________________

void _HYTextBox::_DoCut (bool m)
{
    if (te) {
        gtk_text_buffer_cut_clipboard (gtk_text_view_get_buffer(GTK_TEXT_VIEW(te)), gtk_clipboard_get(GDK_SELECTION_CLIPBOARD), true);
        if (m&&messageRecipient) {
            messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),1));
        }
    }
}

//__________________________________________________________________

void _HYTextBox::_DoCopy (bool)
{
    if (te) {
        gtk_text_buffer_copy_clipboard (gtk_text_view_get_buffer(GTK_TEXT_VIEW(te)), gtk_clipboard_get(GDK_SELECTION_CLIPBOARD));
    }
}

//__________________________________________________________________

void _HYTextBox::_DoPaste (bool m)
{
    if (te) {
        gtk_text_buffer_paste_clipboard (gtk_text_view_get_buffer(GTK_TEXT_VIEW(te)), gtk_clipboard_get(GDK_SELECTION_CLIPBOARD),NULL,true);
        if (m&&messageRecipient) {
            messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),1));
        }
    }
}



//__________________________________________________________________

void _HYTextBox::_DoUndo (bool m)
{
    //TBI
    //if (m&&messageRecipient)
    //  messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),1));
}

//__________________________________________________________________

void _HYTextBox::_DoRedo (bool m)
{
    // TBI
    //if (m&&messageRecipient)
    //  messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),1));
}

//__________________________________________________________________

void _HYTextBox::_DoSelectAll (bool m)
{
    ((_HYTextBox*)this)->SetSelection (0,0x7fffffff);
    if (m&&messageRecipient) {
        messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),0));
    }
}


//__________________________________________________________________

void _HYTextBox::_DoClear (bool doAll, bool m)
{
    if (doAll) {
        _DoSelectAll (false);
    }
    _InsertText (empty, false);
    if (m&&messageRecipient) {
        messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),0));
    }
}


//__________________________________________________________________

void _HYTextBox::_DoFind (_String & st)
{
}


//__________________________________________________________________
void        _HYPlatformTextBox::_Paint (Ptr p)
{
    _HYTextBox * theParent = (_HYTextBox*)this;

    GdkRectangle        cRect = HYRect2GDKRect(*(_HYRect*)p);

    if (!(theParent->settings.width&HY_COMPONENT_TRANSP_BG) || (theParent->boxFlags & HY_TB_BIGBOX) == 0) {
        if (theParent->parentWindow->window) {
            GdkGC *textGC                = gdk_gc_new (theParent->parentWindow->window);
            GdkRegion * r1 = gdk_region_rectangle(&cRect),
                        * r2 = gdk_region_rectangle (&textBoxRect);

            gdk_region_subtract    (r1,r2);
            gdk_region_offset      (r1,theParent->parentWindow->allocation.x,theParent->parentWindow->allocation.y);
            gdk_gc_set_clip_region (textGC,r1);
            if ((theParent->settings.width&HY_COMPONENT_TRANSP_BG) == 0) {
                gdk_gc_set_foreground(textGC,&backFill);
                gdk_draw_rectangle(theParent->parentWindow->window,textGC,true,cRect.x+theParent->parentWindow->allocation.x,
                                   cRect.y+theParent->parentWindow->allocation.y, cRect.width, cRect.height);
            }
            if ((theParent->boxFlags & HY_TB_BIGBOX) == 0) {
                gtk_draw_shadow (gtk_widget_get_style (te), theParent->parentWindow->window, GTK_STATE_NORMAL,
                                 GTK_SHADOW_ETCHED_IN,textBoxRect.x-1+theParent->parentWindow->allocation.x,
                                 textBoxRect.y-1+theParent->parentWindow->allocation.y, textBoxRect.width+2, textBoxRect.height+2);

            }
            gdk_region_destroy(r1);
            gdk_region_destroy(r2);
            g_object_unref (textGC);
        }
    }

    (*theParent)._HYPlatformComponent::_Paint(p);
}

//__________________________________________________________________

bool _HYTextBox::_ProcessOSEvent (Ptr vEvent)
{
    return _HYPlatformComponent::_ProcessOSEvent (vEvent);
}

//__________________________________________________________________
void        _HYPlatformTextBox::_CreateTE (void)
{
    te                  = nil;

    if (((_HYTextBox*)this)->boxFlags & HY_TB_BIGBOX) {
        te                  = gtk_text_view_new  ();
        gtk_text_view_set_wrap_mode (GTK_TEXT_VIEW (te), GTK_WRAP_CHAR);
        _HYTextBox * theParent = (_HYTextBox*)this;
        if (te) {
            gtk_widget_set_app_paintable(te,true);
            gtk_widget_show   (te);
            gtk_container_add (GTK_CONTAINER(theParent->parentWindow), te);
        }
    }
}
//EOF
