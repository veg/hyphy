/*
	A static label object glue for GTK.
	
	Sergei L. Kosakovsky Pond, November 2nd, 2004.
*/

#ifndef _HYPLABEL_
#define _HYPLABEL_

#include "HYPlatformComponent.h"

//__________________________________________________________________

class _HYPlatformLabel {

	public:

		_HYPlatformLabel(void);
		// flags, title, visibility
		
virtual	~_HYPlatformLabel(void);
				
virtual void		_SetBackColor  	 (_HYColor&);
virtual void		_SetForeColor  	 (_HYColor&);
virtual	void		_SetDimensions   (_HYRect,_HYRect);
virtual	void		_SetVisibleSize  (_HYRect);
virtual	void		_SetFont		 (_HYFont&);
virtual	void		_SetText		 (void);
		
virtual	void		_Paint (Ptr p);
virtual	void		_Update(Ptr p);

		GdkGC					*labelGC;
		PangoLayout				*labelLayout;
		PangoFontDescription	*labelFontDesc;
		
		GdkRectangle			labelRect;
};

#endif