/*
    A pull down menu object for MacOS.

    Sergei L. Kosakovsky Pond, May 2000.
*/

#ifndef _HYPPULLDOWNMENU_
#define _HYPPULLDOWNMENU_
//#pragma once
#include "HYBaseGUI.h"
#include "Menus.h"
#include "Lists.h"
#include "TextEdit.h"
#include "strings.h"
#ifdef   TARGET_API_MAC_CARBON
#include "CarbonEvents.h"
#endif



//__________________________________________________________________

class _HYPlatformCheckbox
{

public:

    _HYPlatformCheckbox(bool);
    virtual ~_HYPlatformCheckbox(void);

    virtual void        _SetVisibleSize     (_HYRect);
    virtual void        _SetState           (bool);

    virtual void        _Paint (Ptr p);
    virtual void        _Update(Ptr p);
    virtual void        _Enable(bool);

    ControlHandle   checkboxControl;
    Rect            checkboxRect;
    bool            isRadio;
};

//__________________________________________________________________

class _HYPlatformButtonBar
{

public:

    _HYPlatformButtonBar(void);

    virtual ~_HYPlatformButtonBar(void);

    virtual void        _SetBackColor    (_HYColor&);
    virtual void        _SetDimensions   (_HYRect,_HYRect);
    virtual void        _SetVisibleSize  (_HYRect);
    _HYRect     _GetButtonRect   (bool conv = false);

    virtual void        _Paint (Ptr p);
    virtual void        _Update(Ptr p);
    void        _DisposeButtons (void);
    void        _DisposeButton  (long);
    void        _MarkButtonForUpdate(int);
    void        _UnpushButton   (void);
    int         _FindClickedButton (int,int);


    PixPatHandle backFill;
    Rect         buttonRect,
                 toolTipBounds;
    int          pushed;
    int          saveMousePosH, saveMousePosV;
    unsigned long
    lastSave;
#ifdef      TARGET_API_MAC_CARBON
    EventLoopTimerUPP timerUPP;
    EventLoopTimerRef theTimer;
#endif
};

//__________________________________________________________________

class _HYPlatformList
{

public:

    _HYPlatformList(void);

    virtual ~_HYPlatformList(void);

    virtual void        _SetDimensions   (_HYRect,_HYRect);
    virtual void        _SetVisibleSize  (_HYRect);

    virtual void        _Paint (Ptr p);
    virtual void        _Update(Ptr p);

    virtual void        _InsertItem         (_String&, long);
    virtual void        _SetItem            (_String&, long);
    virtual void        _DeleteItem         (long);
    virtual void        _CheckSelection     (void);
    virtual void        _ToggleMultSelection(bool);
    void        _KillSelection  (void);
    void        _SetSelection   (_SimpleList&);
    void        _SetFont        (void);

    ListHandle      listData;
    short           fontID;
};

//__________________________________________________________________

#define     HY_TABLE_SIZE_CURSOR  0x01
#define     HY_TABLE_DRAG_CURSOR  0x02
#define     HY_TABLE_EDIT_CURSOR  0x04

//__________________________________________________________________

class _HYPlatformTable
{

public:

    _HYPlatformTable        (void);
    virtual             ~_HYPlatformTable       (void);

    void            _SetFont                (void);
    void            _SetBackColor           (_HYColor&);
    void            _SetBackColor2          (_HYColor&);

    void            _CreateTextBox          (_HYRect&,_String&);
    _String         _RetrieveTextValue      (void);
    void            _KillTextBox            (void);
    bool            _HasTextBox             (void) {
        return      editBox;
    }

    Rect            _GetVisibleRowRect      (long);

    void            _HiliteRowForDrag       (long,long);

    short           fontID;
    PixPatHandle    backPattern,
                    backPattern2;
    char            cursorState;

    TEHandle        editBox;
    Rect            textBoxRect;
};



//__________________________________________________________________

class _HYPlatformTextBox
{

public:

    _HYPlatformTextBox (void);

    virtual ~_HYPlatformTextBox(void);

    virtual void            _SetBackColor    (_HYColor&);
    virtual void            _SetBackTColor   (_HYColor&);
    virtual void            _SetForeColor    (_HYColor&);
    virtual void            _SetDimensions   (_HYRect,_HYRect);
    virtual void            _SetVisibleSize  (_HYRect);
    virtual void            _SetFont         (_HYFont&);
    void            _SetText         (_String&);
    void            _InsertText      (_String&);
    _String         _GetText         (void);
    void            _CreateTE        (void);
    void            _EnableTextBox   (bool);
    void            _SetMargins      (_HYRect&);
    void            _SetAlignFlags   (unsigned char);

    virtual void            _Paint (Ptr p);
    virtual void            _Update(Ptr p);

    virtual void            _FocusComponent  (void);
    virtual void            _UnfocusComponent(void);

    PixPatHandle    backFill,
                    backTFill;
    TEHandle        te;
    Rect            textBoxRect;
};
#endif

#ifdef TARGET_API_MAC_CARBON
pascal void ButtonBarTimer (EventLoopTimerRef theTimer,void* userData);
#endif


//EOF
