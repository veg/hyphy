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
    Graphics Window Object  glue for GTK+

    Sergei L. Kosakovsky Pond, October 2004
*/

#include "HYGWindow.h"
#include "HYEventTypes.h"
#include <gtk/gtk.h>
#include <gdk/gdk.h>
//#include "HYUtils.h"

//__________________________________________________________________

void    _HYGWindow::_Paint (Ptr)
{
    int t,l,b,r;
    _VisibleContents (t,l,b,r);
    //printf ("GPaint on %d %d %d %d\n", t, l, b, r);
    gdk_draw_drawable (GDK_DRAWABLE(windowContent->window), theContext, thePane, l, t,
                       windowContent->allocation.x, windowContent->allocation.y, r-l+1, b-t+1);
}

//__________________________________________________________________

void    _HYGWindow::_Update (Ptr)
{
    _Paint (nil);
}

//__________________________________________________________________


bool        _HYGWindow::_ProcessMenuSelection (long msel)
{
    bool        done = false;
    switch (msel) {
    case HY_WINDOW_MENU_ID_EDIT+1: {
        _CopyToClipboard ();
        done = true;
    }
    break;

    case HY_WINDOW_MENU_ID_FILE+1: {
        _SavePicture (GetTitle());
        done = true;
    }
    break;
    }

    return _HYWindow::_ProcessMenuSelection(msel);
}

//__________________________________________________________________

void _HYGWindow::_SetMenuBar(void)
{
    _HYWindow::_SetMenuBar();
    gtk_widget_set_sensitive(gtk_item_factory_get_widget (menu_items,"<HY_WINDOW>/Edit/Copy"),  TRUE);
}

//__________________________________________________________________

void _HYGWindow::_UnsetMenuBar(void)
{
}


//EOF
