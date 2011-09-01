/*
	A button object glue for GTK.

	Sergei L. Kosakovsky Pond, November 2004.
*/

#ifndef _HYPBUTTON_
#define _HYPBUTTON_

#include "HYPlatformComponent.h"

//__________________________________________________________________

class _HYPlatformButton
{

public:

	_HYPlatformButton(void);

	virtual	~_HYPlatformButton(void);

	virtual void			_SetBackColor  	 (_HYColor&);
	virtual	void			_SetDimensions   (_HYRect,_HYRect);
	virtual	void			_SetVisibleSize  (_HYRect);
	virtual	void			_SetFont		 (_HYFont&);
	void			_SetText		 (void);
	void			_SetButtonKind	 (unsigned char);
	void			_EnableButton	 (bool);
	void			_ApplyFont		 (void);
	void			_PaintMe		 (void);

	virtual	void			_Paint (Ptr p);
	virtual	void			_Update(Ptr p);

	GdkRectangle			buttonRect;
	GdkColor				bgColor;
	GtkWidget*				buttonControl;
	PangoFontDescription*   buttonFontDesc;
	_HYRect					lastButtonDimension;
};

#endif