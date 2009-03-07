/*
	Pull-down menu component
	
	Sergei L. Kosakovsky Pond, May 2000 - November 2001.
*/

#include "HYEventTypes.h"
#include "HYPullDown.h"

#ifdef 	  __HYPHYDMALLOC__
	#include "dmalloc.h"
#endif

_String			menuSeparator 		("SEPARATOR");

//__________________________________________________________________

_HYPullDown::_HYPullDown(_HYRect r,Ptr p):_HYComponent (r,p)
{
	backColor.R = backColor.G = backColor.B = 255;
	enabledFlag = true;
	alignFlags = 0;
}

//__________________________________________________________________
			
_HYPullDown::~_HYPullDown()
{
}
		
//__________________________________________________________________
BaseRef		_HYPullDown::makeDynamic()
{
	/*_HYPullDown* res = new _HYPullDown();
	memcpy ((Ptr)res,(Ptr)this,sizeof (_HYPullDown));
	res->menuSelections.Duplicate (&menuSelections);
	res->_Duplicate ((Ptr)this);
	return res;*/
	return nil;
}
		
//__________________________________________________________________
void		_HYPullDown::SetBackColor (_HYColor c)
{
	if ((c.R!=backColor.R)||(c.G!=backColor.G)||(c.B!=backColor.B))
	{
		backColor = c;
		_SetBackColor (c);
		_MarkForUpdate();
	}
}


//__________________________________________________________________

_HYColor&		_HYPullDown::GetBackColor (void)
{
	return backColor;
}

//__________________________________________________________________

void			_HYPullDown::AddMenuItem  (_String newItem, long loc)
{
	menuSelections.InsertElement (&newItem,loc,true);
	_AddMenuItem (newItem,loc);
}

//__________________________________________________________________

void			_HYPullDown::SetMenuItem  (_String newItem, long loc)
{
	if ((loc>=0)&&(loc<menuSelections.lLength))
	{
		menuSelections.Replace (loc,&newItem,true);
		_SetMenuItem (newItem,loc);
	}
}

//__________________________________________________________________

void			_HYPullDown::DeleteMenuItem  (long loc)
{
	if ((loc>=0)&&(loc<menuSelections.lLength))
	{
		menuSelections.Delete (loc);
		_DeleteMenuItem (loc);
	}
}

//__________________________________________________________________

void			_HYPullDown::DeleteAllItems  (void)
{
	for (long k=menuSelections.lLength-1; k>=0; k--)
		_DeleteMenuItem (k);
		
	menuSelections.Clear ();
}

//__________________________________________________________________
_String*		_HYPullDown::GetMenuItem  (long loc)
{
	if ((loc>=0)&&(loc<menuSelections.lLength))
	{
		return (_String*)menuSelections (loc);
	}
	return nil;
}

//__________________________________________________________________
long		_HYPullDown::GetSelection (void)
{
	return _GetSelection ();
}

//__________________________________________________________________

void			_HYPullDown::Activate  (void)
{
	_HYComponent::Activate();
	if (IsEnabled())
		_MarkForUpdate();
}

//__________________________________________________________________

void			_HYPullDown::Deactivate  (void)
{
	_HYComponent::Deactivate();
	if (IsEnabled())
		_MarkForUpdate();
}



//__________________________________________________________________
void		_HYPullDown::SendSelectionChange (void)
{
	if (messageRecipient)
	{
		messageRecipient->ProcessEvent (generateMenuSelChangeEvent (GetID(),GetSelection()));
	}
}

//__________________________________________________________________
void		_HYPullDown::EnableMenu 	 (bool flag)
{
	if (flag!=enabledFlag)
	{
		enabledFlag = flag;
		_EnableMenu    (flag);
		_MarkForUpdate ();
	}
}

//__________________________________________________________________
bool		_HYPullDown::IsEnabled	 (void)
{
	return enabledFlag;
}

//__________________________________________________________________
void		_HYPullDown::SetVisibleSize	 (_HYRect rel)
{
	_HYComponent::SetVisibleSize (rel);
	_HYPlatformPullDown::_SetVisibleSize (rel);
}

//__________________________________________________________________
void		_HYPullDown::ChangeSelection	 (long newSel, bool eventSend)
{
	if ((newSel>=0)&&(newSel<menuSelections.lLength))
	{
		selection = newSel;
		if (eventSend)
			SendSelectionChange();
		_MarkForUpdate();
	}
}

//__________________________________________________________________
void		_HYPullDown::EnableItem	 (long theItem, bool toggle)
{
	if ((theItem>=0)&&(theItem<menuSelections.lLength))
	{
		_EnableItem (theItem,toggle);
	}
}

//__________________________________________________________________
bool		_HYPullDown::ProcessEvent (_HYEvent* e) 
{
	DeleteObject (e); 
	return false;
}