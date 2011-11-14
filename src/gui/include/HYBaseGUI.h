/*
    Some general GUI object definitions.
    Basic GUI object,
    Basic Event.

    Global event queue, and object ID counter.

    Sergei L. Kosakovsky Pond, May 2000.
*/

#ifndef _HYBASEGUI_
#define _HYBASEGUI_
//#pragma once
#include "hy_strings.h"
#include "list.h"

//__________________________________________________________________

class _HYEvent: public BaseObj
{

public:

    _HYEvent();

    _HYEvent(_HYEvent&);

    _HYEvent(_String, _String);

    virtual ~_HYEvent(void);

    virtual BaseRef  toStr (void);

    virtual BaseRef  makeDynamic (void);

    virtual void     Duplicate (BaseObj* ref);

    _String& EventClass (void);

    _String& EventCode  (void);

private:

    _String  eventClass, eventCode;
};

//__________________________________________________________________

class _HYGuiObject:public BaseObj
{

public:

    _HYGuiObject (void);

    virtual ~_HYGuiObject (void);

    unsigned long GetID (void);

    bool     MatchID (unsigned long);

    virtual bool     ProcessEvent (_HYEvent* e)     {
        DeleteObject(e);
        return false;
    }
    virtual bool     ProcessGEvent (_HYEvent*)      {
        return false;
    }

    virtual void     ProcessContextualPopUp (long, long) {}

private:

    unsigned long ID;
};



//__________________________________________________________________

typedef struct _HYRect {
    long        top,
                left,
                bottom,
                right,
                width;

    bool        Contains (int, int);
    long        Width    ();
    long        Height   ();

} _HYRect;

//__________________________________________________________________

typedef struct _HYColor {
    unsigned char R,G,B;

    _String  HTMLColor (void);
    bool     operator == (_HYColor&);
} _HYColor;

//__________________________________________________________________

typedef struct  _HYFont {
    _String face;
    int     size;
    unsigned char style;
} _HYFont;

//__________________________________________________________________


long     HYColorToLong      (_HYColor);
_HYColor LongToHYColor      (long);

_HYRect  makeHYRect         (int,int,int,int,int);

extern  _List               GlobalGUIEventQueue;
extern  _SimpleList         windowObjectRefs;
extern  _String             ReportThisError;
extern  unsigned long       GUIObjectGlobalCounter;

#endif

//EOF
