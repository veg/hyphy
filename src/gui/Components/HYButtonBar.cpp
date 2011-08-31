/*
	Iconic button bar (toolbar) component
	
	Sergei L. Kosakovsky Pond, May 2000 - December 2002.
*/

#include "HYEventTypes.h"
#include "HYButtonBar.h"
#include "HYGraphicPane.h"

#ifdef 	  __HYPHYDMALLOC__
	#include "dmalloc.h"
#endif

//__________________________________________________________________

_HYButtonBar::_HYButtonBar(_HYRect r,Ptr p):_HYComponent (r,p)
{
	backColor.R = backColor.G = backColor.B = 255;
	barW = 100;
	buttonDim = 32;
	alignFlags = 0;
}

//__________________________________________________________________
			
_HYButtonBar::~_HYButtonBar()
{
	_DisposeButtons();
}

//__________________________________________________________________
void		_HYButtonBar::SetBackColor (_HYColor c)
{
	if ((c.R!=backColor.R)||(c.G!=backColor.G)||(c.B!=backColor.B))
	{
		backColor = c;
		_SetBackColor (c);
		_MarkForUpdate();
	}
}


//__________________________________________________________________

_HYColor&		_HYButtonBar::GetBackColor (void)
{
	return backColor;
}		

//__________________________________________________________________

void		_HYButtonBar::AddButton (Ptr p, _String* toolTip)
{
	if (p)
	{
		enabledButtons<<buttons.lLength;
		buttons<<(long)p;
		if (toolTip)
			toolTips && toolTip;
		else
			toolTips && & empty;
		_MarkForUpdate();
	}
}	

//__________________________________________________________________

void		_HYButtonBar::SetToolTip (long index, _String* toolTip)
{
	if ((index>=0)&&(index<buttons.lLength))
		toolTips.Replace (index,toolTip,true);
}	

//__________________________________________________________________

void		_HYButtonBar::ReplaceButton (long index,Ptr p, _String * tt)
{
	if (p)
	{
		if ((index>=0)&&(index<=buttons.lLength))
		{
			_DisposeButton (index);
			buttons.lData[index] = (long)p;
			if (tt)
				toolTips.Replace (index, tt, true);
			_MarkButtonForUpdate (index);
		}
	}
}	

//__________________________________________________________________

Ptr		_HYButtonBar::GetButtonIcon (long index)
{
	if ((index>=0)&&(index<=buttons.lLength))
		return (Ptr)buttons.lData[index];
		
	return nil;
}	
//__________________________________________________________________

void		_HYButtonBar::DeleteButton (long index)
{	
	long f = enabledButtons.BinaryFind(index);
	if (f>=0)
		enabledButtons.Delete(f);
	f = pullDownButtons.Find (index);
	if (f>=0)
		pullDownButtons.Delete (f);
	buttons.Delete(index);
	toolTips.Delete (index);
	_MarkForUpdate();
}	

//__________________________________________________________________

void		_HYButtonBar::MarkAsPullDown (long index, bool add)
{	
	long f = pullDownButtons.Find (index);
	if ((f>=0)&&(!add))
		pullDownButtons.Delete (f);
	if (add&&(f<=0))
		pullDownButtons << index;
}

//__________________________________________________________________

void		_HYButtonBar::EnableButton (long index, bool onOff)
{	
	long f = enabledButtons.Find(index);
	if (onOff)
	{
		if (f<0)
		{
			_MarkButtonForUpdate(index);
			enabledButtons.InsertElement ((BaseRef)index,-1,false,false);
		}
	}
	else
	{
		if (f>=0)
		{
			enabledButtons.Delete(f);
			_MarkButtonForUpdate(index);
		}
	}
}	


//__________________________________________________________________

void		_HYButtonBar::GetButtonLoc (int b, int& h, int& v, bool c)
{
	_HYRect    cRect = _GetButtonRect(c);
 	
	int	  step = GetButtonDim()+2*HY_BUTTONBAR_BORDER;

 	h = cRect.left+ (b%barW)*step;
 	v = cRect.top + (b/barW)*step;
}

//__________________________________________________________________
void		_HYButtonBar::SendButtonPush (int bID)
{
	if (messageRecipient)
		messageRecipient->ProcessEvent (generateButtonPushEvent (GetID(),bID));
}

//__________________________________________________________________
void		_HYButtonBar::SetVisibleSize	 (_HYRect rel)
{
	_HYComponent::SetVisibleSize (rel);
	_HYPlatformButtonBar::_SetVisibleSize (rel);
}	

//__________________________________________________________________
bool		_HYButtonBar::ProcessEvent (_HYEvent* e) 
{
	DeleteObject (e); 
	return false;
}

//__________________________________________________________________
 _HYRect	_HYButtonBar::_SuggestDimensions (void)
 {
  	_HYRect res = {10,10,10,10,HY_COMPONENT_NO_SCROLL};
 	int w = ButtonCount(), 
 		h=1;
 	
 	if (w>BarWidth())
 	{
 		w = BarWidth();
 		h = ButtonCount()/w+(ButtonCount()%w>0);
 	}
 	
 	res.right =  w*(buttonDim+2*HY_BUTTONBAR_BORDER);
  	res.bottom = h*(buttonDim+2*HY_BUTTONBAR_BORDER);
	return res;
 }