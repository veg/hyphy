/*
	A general composite window component object, FLTK glue
	
	Sergei L. Kosakovsky Pond, October 2004.
*/

#ifndef _HYPCOMPONENT_
#define _HYPCOMPONENT_
//#pragma once
#include "HYBaseGUI.h"

#include <gtk/gtk.h>

//__________________________________________________________________

class _HYPlatformComponent{

	public:

		_HYPlatformComponent(void);
		_HYPlatformComponent(_HYRect,Ptr);
		// settings
		
virtual	~_HYPlatformComponent() {};
		
virtual	void		Duplicate (BaseRef);
		
virtual	void		_CleanUp   (void);
		
virtual	void		_SetDimensions (_HYRect,_HYRect);
virtual void		_SetVisibleSize(_HYRect);

virtual	void		_Paint (Ptr);
virtual void		_Update (Ptr);
virtual bool		_ProcessOSEvent (Ptr);
virtual	_HYRect		_VisibleContents(Ptr);
virtual	void		_MarkForUpdate (void);
virtual	void		_MarkContentsForUpdate (void);
virtual long		_GetScrollerPos	 (GtkWidget*);
virtual	long		_GetHScrollerPos (void);
virtual	long		_GetVScrollerPos (void);
virtual	void		_SetScrollerPos  (GtkWidget*,long);
virtual	void		_SetHScrollerPos (long);
virtual	void		_SetVScrollerPos (long);
virtual void		_Activate (void);
virtual void		_Deactivate(void);
virtual void		_ComponentMouseExit (void) {}

				

		GtkWidget	  *parentWindow,
					  *vScroll,
					  *hScroll;

		_HYRect		  rel;
		
		bool		  activationFlag;
		long		  lastHScroll,
					  lastVScroll;
};

//extern  bool 	forceUpdateForScrolling;

void			AlignRectangle (_HYRect& rel , GdkRectangle& target , unsigned char alFlags);


#ifdef TARGET_API_MAC_CARBON

extern bool aquaInterfaceOn;

#endif

//__________________________________________________________________

#endif

//EOF