/*
	A text input box object glue for GTK+

	Sergei L. Kosakovsky Pond, November 2004
*/

#ifndef _HYPLTEXTBOX_
#define _HYPLTEXTBOX_

#include "HYPlatformComponent.h"

//__________________________________________________________________

class _HYPlatformTextBox
{

public:

	_HYPlatformTextBox (void);

	virtual	~_HYPlatformTextBox(void);

	virtual void			_SetBackColor  	 (_HYColor&);
	virtual void			_SetBackTColor   (_HYColor&);
	virtual void			_SetForeColor  	 (_HYColor&);
	virtual	void			_SetDimensions   (_HYRect,_HYRect);
	virtual	void			_SetVisibleSize  (_HYRect);
	virtual	void			_SetFont		 (_HYFont&);
	void			_SetText		 (const _String&);
	void			_InsertText		 (const _String&, bool);
	_String			_GetText		 (void);
	void			_StoreText		 (_String*&, bool);
	void			_CreateTE		 (void);
	void			_EnableTextBox	 (bool);
	void			_SetMargins		 (_HYRect&);
	void			_SetAlignFlags	 (unsigned char);

	virtual	void			_Paint (Ptr p);
	virtual	void			_Update(Ptr p);

	virtual	void			_FocusComponent  (void);
	virtual	void			_UnfocusComponent(void);
	virtual bool			_NeedMultiLines  (void);

	GdkColor		backFill,
					backTFill,
					textColor;

	GtkWidget		*te,
					*scrollWindow;

	GdkRectangle	textBoxRect;
	PangoFontDescription
	*pLabelFont;

	bool			isSingleLine;
};

//__________________________________________________________________

//_String					retrieveEditControlText (HWND);
//void					retrieveEditControlText (HWND, _String*&);

#endif