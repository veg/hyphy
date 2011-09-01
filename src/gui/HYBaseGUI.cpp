/*
	Some general GUI object functionality.
	Basic GUI object,
	Basic Event.

	Global event queue, and object ID counter.

	Sergei L. Kosakovsky Pond, May 2000.
*/

#include "HYBaseGUI.h"
#include "baseobj.h"

#ifdef 	  __HYPHYDMALLOC__
#include "dmalloc.h"
#endif
//__________________________________________________________________

unsigned long GUIObjectGlobalCounter = 0;
_List	 GlobalGUIEventQueue;
_String	 ReportThisError
(" Please report this error using a bug report form at http://peppercat.stat.ncsu.edu.");

//______________________________________________________________
//  Function Definitions for _HYGuiObject
//______________________________________________________________

_HYGuiObject::_HYGuiObject (void)
{
	ID = ++GUIObjectGlobalCounter;
}

//__________________________________________________________________

_HYGuiObject::~_HYGuiObject (void)
{
}

//__________________________________________________________________

bool _HYGuiObject::MatchID (unsigned long testID)
{
	return ID==testID;
}

//__________________________________________________________________

unsigned long _HYGuiObject::GetID (void)
{
	return ID;
}


//______________________________________________________________
//  Function Definitions for _HYEvent
//______________________________________________________________

_HYEvent::_HYEvent() {}

//______________________________________________________________

_HYEvent::_HYEvent(_HYEvent& he)
{
	Duplicate (&he);
}

//______________________________________________________________

_HYEvent::_HYEvent(_String eCl, _String eCo)
{
	eventClass.Duplicate(&eCl);
	eventCode.Duplicate(&eCo);
}

//______________________________________________________________

_HYEvent::~_HYEvent(void) {}

//______________________________________________________________

BaseRef _HYEvent::toStr (void)
{
	_String result ("Event class = ");
	result = result & eventClass & " code = " & eventCode;
	return result.makeDynamic();
}

//______________________________________________________________

BaseObj* _HYEvent::makeDynamic (void)
{
	return new _HYEvent (*this);
}

//______________________________________________________________

void	_HYEvent::Duplicate (BaseObj* ref)
{
	_HYEvent* he = (_HYEvent*)ref;
	eventClass.Duplicate (&he->eventClass);
	eventCode.Duplicate(&he->eventCode);
}

//______________________________________________________________

_String& _HYEvent::EventClass (void)
{
	return eventClass;
}

//______________________________________________________________

_String& _HYEvent::EventCode (void)
{
	return eventCode;
}

//______________________________________________________________

_HYRect	 makeHYRect (int t,int l,int b,int r,int w)
{
	_HYRect res;
	res.top = t;
	res.bottom = b;
	res.left = l;
	res.right = r;
	res.width = w;
	return res;

}

//______________________________________________________________

bool	_HYRect::Contains (int h, int v)
{
	return ((h>=left)&&(h<=right)&&(v>=top)&&(v<=bottom));
}

//______________________________________________________________

long	_HYRect::Width (void)
{
	return right-left+1;
}

//______________________________________________________________

long	_HYRect::Height (void)
{
	return bottom-top+1;
}



//______________________________________________________________

_String	_HYColor::HTMLColor (void)
{
	char out[7];

	sprintf (out, "%x%x%x", R, G, B);

	return _String (out);
}

//______________________________________________________________

bool	_HYColor::operator ==  (_HYColor& c)
{
	return ((c.R==R)&&(c.G==G)&&(c.B==B));
}

//______________________________________________________________

long	HYColorToLong (_HYColor c)
{
	return (long)c.R*0x00010000+(long)c.G*0x00000100+(long)c.B;
}

//______________________________________________________________

_HYColor
LongToHYColor (long c)
{
	_HYColor hc;
	hc.B = c%256;
	c=c>>8;
	hc.G = c%256;
	c=c>>8;
	hc.R = c;
	return hc;
}





//EOF