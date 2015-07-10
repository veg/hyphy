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

#ifndef _HYPWINDOW_
#define _HYPWINDOW_

#include <gtk/gtk.h>

#define  MAX_CONTROL_VALUE 100000000.
#define  HY_WINDOW_MENU_ID_FILE 5000
#define  HY_WINDOW_MENU_ID_EDIT 5050
#define  HY_SCROLLER_WIDTH      16

#include "hy_strings.h"
#include "HYBaseGUI.h"

//__________________________________________________________________

class _HYPlatformWindow
{

public:

    _HYPlatformWindow   (unsigned char,_String,bool,Ptr = nil);
    // flags, title, visibility

    virtual ~_HYPlatformWindow  (void);

    void        _SetTitle               (const _String&);
    void        _Show                   (void);
    void        _Hide                   (void);
    virtual bool        _Close                  (Ptr);
    virtual void        _SetPosition            (int,int);
    void        _SelectWindow           (void);
    virtual void        _SetWindowRectangle     (int,int,int,int,bool=true);
    virtual void        _SetContentSize         (int,int);
    virtual void        _Paint                  (Ptr);
    virtual void        _PaintHook              (Ptr);
    virtual void        _Update                 (Ptr);
    virtual void        _Activate               (void);
    virtual void        _BringWindowToFront     (void);
    virtual void        _Deactivate             (void);
    virtual bool        _ProcessOSEvent         (Ptr);
    virtual long        _Grow                   (Ptr) {
        return 0;
    }
    virtual void        _Move                   (Ptr) {}
    virtual void        _VisibleContents        (int&,int&,int&,int&);
    virtual _HYRect     _GetWindowRect          (void);
    bool        _IsHScroll              (GtkWidget*);
    virtual bool        _ProcessMenuSelection   (long) {
        return false;
    }
    virtual void        _SetMenuBar             (void);
    virtual void        _UnsetMenuBar           (void);
    virtual _String&    _GetTitle               (void);
    void        _SetWindowBackColor     (_HYColor);
    bool        _CleanDefaultMenu       (void);

    virtual Ptr         _GetOSWindowData (void) {
        return (Ptr)(windowContent);
    }



    GtkWidget                           *vScroll,
                                        *hScroll,
                                        *windowMB,
                                        *theWindow,
                                        *windowContent,
                                        *hbox,
                                        *vbox;

    GtkItemFactory                      *menu_items;

    double                              last_H_Position,
                                        last_V_Position;

    _HYRect                             savedLoc;
    Ptr                                 containerRef;
    unsigned char                       flags;
    long                                lastWW,
                                        lastWH;
};

//__________________________________________________________________

class _HYPlatformTWindow
{

public:

    _HYPlatformTWindow  (Ptr);
    // flags, title, visibility

    virtual ~_HYPlatformTWindow (void);
    virtual void        _SetWindowRectangle     (int,int,int,int,bool=true);

    guint           theTimer;
    Ptr             trackMouseComponent;


};

//__________________________________________________________________

class _HYPlatformPWindow
{

public:

    _HYPlatformPWindow          (void);

    virtual ~_HYPlatformPWindow         (void);

    virtual void _StartPicture          (void);
    virtual void _EndPicture            (void);
    virtual void _DrawPicture           (_HYRect);

private:

    //PicHandle              savedPic;
    //RgnHandle          savedClip;

};

//__________________________________________________________________

struct _HY_GTK_UI_Message {
    gboolean   processingResult;
    GdkEvent * theEvent;
};

//__________________________________________________________________

void     AdjustScroller     (GtkWidget*, long, long);

extern   bool               forceUpdateForScrolling;
extern   bool               InitPrint (void);

extern   GdkPixmap          *stripedFill;

extern   GdkPixbuf          *redButtonIcon,
         *yellowButtonIcon,
         *greenButtonIcon,
         *orangeButtonIcon;

extern   GdkColor           _BLACKBRUSH_,
         _DARKGREYBRUSH_,
         _WHITEBRUSH_;

extern   GdkGC              *stripedFillGC;
extern   PangoLayout        *statusBarLayout;
extern   PangoFontDescription
*statusBarFontDesc;

extern   _HYFont            statusBarFont;

void                        _PaintTheCircle                 (GdkPixbuf*,GtkWidget*);

void                        h_scroll_bar_callback_window    (GtkRange*, gpointer);
void                        v_scroll_bar_callback_window    (GtkRange*, gpointer);
void                        hyphy_menu_item_callback         (gpointer, guint, GtkWidget *);

GdkColor                    HYColorToGDKColor               (_HYColor);

void                        SetUpStatusBarStuff             (GtkWidget *);


#endif

//EOF
