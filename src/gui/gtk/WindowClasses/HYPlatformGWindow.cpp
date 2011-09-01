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

void	_HYGWindow::_Paint (Ptr)
{
	int t,l,b,r;
	_VisibleContents (t,l,b,r);
	//printf ("GPaint on %d %d %d %d\n", t, l, b, r);
	gdk_draw_drawable (GDK_DRAWABLE(windowContent->window), theContext, thePane, l, t,
					   windowContent->allocation.x, windowContent->allocation.y, r-l+1, b-t+1);
}

//__________________________________________________________________

void	_HYGWindow::_Update (Ptr)
{
	_Paint (nil);
}

//__________________________________________________________________


bool 		_HYGWindow::_ProcessMenuSelection (long msel)
{
	bool		done = false;
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