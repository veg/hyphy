/*
	At toolbar menu object glue for GTK+
	
	Sergei L. Kosakovsky Pond, November 2004.
*/

#ifndef _HYPBUTTONBAR_
#define _HYPBUTTONBAR_

#include "HYPlatformComponent.h"

//__________________________________________________________________

class _HYPlatformButtonBar {

	public:

		_HYPlatformButtonBar(void);
		
virtual	~_HYPlatformButtonBar(void);
				
virtual void		_SetBackColor  	 (_HYColor&);
virtual	void		_SetDimensions   (_HYRect,_HYRect);
virtual	void		_SetVisibleSize  (_HYRect);
virtual	void		_SmiteTooltip	 (void);
		_HYRect		_GetButtonRect   (bool conv = false);
		
virtual	void		_Paint (Ptr p);
virtual	void		_Update(Ptr p);
		void		_DisposeButtons (void);
		void		_DisposeButton	(long);
		void		_MarkButtonForUpdate(int);
		void		_UnpushButton	(void);
		int			_FindClickedButton (int,int);


		GdkColor		backFill;
		GdkRectangle	buttonRect,
						toolTipBounds;
		
		GdkGC			*bbGC;
					 
		int				pushed,
						saveMousePosH, 
						saveMousePosV,
						lastMouseDown;
						
		guint			theTimer;

					 
		//unsigned long	lastSave;
		
};

//__________________________________________________________________



#endif