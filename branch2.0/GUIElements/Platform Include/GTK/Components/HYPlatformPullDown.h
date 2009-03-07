/*
	A pull down menu object glue for GTK+.
	
	Sergei L. Kosakovsky Pond, November 2004.
*/

#ifndef _HYPPULLDOWNMENU_
#define _HYPPULLDOWNMENU_

#include "HYPlatformComponent.h"

//__________________________________________________________________

extern   _String	menuSeparator;

//__________________________________________________________________

class _HYPlatformPullDown {

	public:

		_HYPlatformPullDown(void);
		// flags, title, visibility
		
virtual	~_HYPlatformPullDown(void);
				
virtual void		_AddMenuItem  	 (_String&, long);
virtual void		_SetMenuItem  	 (_String&, long);
virtual void		_SetBackColor  	 (_HYColor&);
virtual void		_Duplicate		 (Ptr);
virtual void		_DeleteMenuItem  (long);
virtual long		_GetSelection 	 (void);
virtual	void		_SetDimensions   (_HYRect,_HYRect);
virtual	void		_SetVisibleSize  (_HYRect);
virtual	void		_EnableItem  	 (long, bool);
virtual void		_MarkItem		 (long, char);
virtual char		_ItemMark		 (long);
virtual void		_EnableMenu		 (bool) {}

		
virtual	void		_Paint (Ptr p);
virtual	void		_Update(Ptr p);
virtual	void		_RefreshComboBox  (void);

		GtkWidget*		theMenu;
						
		GdkColor		backFill;
		long			selection,
						cbSelection;
		GdkRectangle	menuRect;
		_SimpleList		widgetList;
};

//__________________________________________________________________

extern	GdkColor	buttonBorder1,
					buttonBorder2;


#endif