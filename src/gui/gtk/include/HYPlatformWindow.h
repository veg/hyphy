/*
	GTK window object glue - a window/title/size-box/scroll-bars handler.
	
	Sergei L. Kosakovsky Pond, October 2004.
*/

#ifndef _HYPWINDOW_
#define _HYPWINDOW_

#include <gtk/gtk.h>

#define  MAX_CONTROL_VALUE 100000000.
#define	 HY_WINDOW_MENU_ID_FILE 5000
#define	 HY_WINDOW_MENU_ID_EDIT 5050
#define  HY_SCROLLER_WIDTH		16

#include "hy_strings.h"
#include "HYBaseGUI.h"

//__________________________________________________________________

class _HYPlatformWindow {

	public:

		_HYPlatformWindow	(unsigned char,_String,bool,Ptr = nil);
		// flags, title, visibility
		
virtual ~_HYPlatformWindow	(void);
		
		void		_SetTitle 				(const _String&);
		void		_Show 					(void);
		void		_Hide 					(void);
virtual	bool		_Close					(Ptr);
virtual void		_SetPosition 			(int,int);
		void		_SelectWindow 			(void);
virtual void		_SetWindowRectangle 	(int,int,int,int,bool=true);
virtual void		_SetContentSize	    	(int,int);
virtual void		_Paint 					(Ptr);
virtual void		_PaintHook 				(Ptr);
virtual void		_Update 				(Ptr);
virtual void		_Activate 				(void);
virtual void		_BringWindowToFront 	(void);
virtual void		_Deactivate 			(void);
virtual bool		_ProcessOSEvent 		(Ptr);
virtual long		_Grow					(Ptr) {return 0;}
virtual void		_Move					(Ptr) {}
virtual void		_VisibleContents 		(int&,int&,int&,int&);
virtual _HYRect		_GetWindowRect			(void);
		bool		_IsHScroll 		 		(GtkWidget*);
virtual	bool		_ProcessMenuSelection 	(long) 
												{return false;}
virtual	void		_SetMenuBar 			(void);
virtual	void		_UnsetMenuBar 			(void);
virtual _String&	_GetTitle 				(void);
		void		_SetWindowBackColor		(_HYColor);
		bool		_CleanDefaultMenu		(void);

virtual Ptr			_GetOSWindowData (void) {return (Ptr)(windowContent);}



		GtkWidget   					    *vScroll,
											*hScroll,
											*windowMB,
											*theWindow,
											*windowContent,
											*hbox,
											*vbox;
											
		GtkItemFactory						*menu_items;
											
		double								last_H_Position,
											last_V_Position;
											
		_HYRect								savedLoc;
		Ptr									containerRef;
		unsigned char 						flags;
		long								lastWW,
											lastWH;
};

//__________________________________________________________________

class _HYPlatformTWindow  {

	public:

		_HYPlatformTWindow	(Ptr);
		// flags, title, visibility
		
virtual ~_HYPlatformTWindow	(void);
virtual void		_SetWindowRectangle 	(int,int,int,int,bool=true);

		guint			theTimer;
		Ptr				trackMouseComponent;
					
		
};

//__________________________________________________________________

class _HYPlatformPWindow  {

	public:
	
		_HYPlatformPWindow			(void);

virtual ~_HYPlatformPWindow			(void);

virtual	void _StartPicture			(void);
virtual void _EndPicture    		(void);
virtual void _DrawPicture			(_HYRect);

	private:

		//PicHandle 			 savedPic;
		//RgnHandle			 savedClip;

};

//__________________________________________________________________

struct _HY_GTK_UI_Message {
	gboolean   processingResult;
	GdkEvent * theEvent;
};

//__________________________________________________________________

void	 AdjustScroller 	(GtkWidget*, long, long);

extern	 bool				forceUpdateForScrolling;
extern	 bool			 	InitPrint (void);

extern	 GdkPixmap  		*stripedFill;

extern	 GdkPixbuf			*redButtonIcon,
						 	*yellowButtonIcon,
						 	*greenButtonIcon,
						 	*orangeButtonIcon;
							
extern	 GdkColor			_BLACKBRUSH_,
							_DARKGREYBRUSH_,
							_WHITEBRUSH_;
							
extern   GdkGC				*stripedFillGC;
extern	 PangoLayout		*statusBarLayout;
extern   PangoFontDescription		
							*statusBarFontDesc;
							
extern   _HYFont			statusBarFont;

void						_PaintTheCircle 				(GdkPixbuf*,GtkWidget*);

void						h_scroll_bar_callback_window	(GtkRange*, gpointer);
void						v_scroll_bar_callback_window	(GtkRange*, gpointer);
void						hyphy_menu_item_callback		 (gpointer, guint, GtkWidget *);

GdkColor					HYColorToGDKColor				(_HYColor);

void						SetUpStatusBarStuff				(GtkWidget *);


#endif

//EOF