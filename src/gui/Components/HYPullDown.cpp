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

#include "HYEventTypes.h"
#include "HYPullDown.h"

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

_String         menuSeparator       ("SEPARATOR");

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
BaseRef     _HYPullDown::makeDynamic()
{
    /*_HYPullDown* res = new _HYPullDown();
    memcpy ((Ptr)res,(Ptr)this,sizeof (_HYPullDown));
    res->menuSelections.Duplicate (&menuSelections);
    res->_Duplicate ((Ptr)this);
    return res;*/
    return nil;
}

//__________________________________________________________________
void        _HYPullDown::SetBackColor (_HYColor c)
{
    if ((c.R!=backColor.R)||(c.G!=backColor.G)||(c.B!=backColor.B)) {
        backColor = c;
        _SetBackColor (c);
        _MarkForUpdate();
    }
}


//__________________________________________________________________

_HYColor&       _HYPullDown::GetBackColor (void)
{
    return backColor;
}

//__________________________________________________________________

void            _HYPullDown::AddMenuItem  (_String newItem, long loc)
{
    menuSelections.InsertElement (&newItem,loc,true);
    _AddMenuItem (newItem,loc);
}

//__________________________________________________________________

void            _HYPullDown::SetMenuItem  (_String newItem, long loc)
{
    if ((loc>=0)&&(loc<menuSelections.lLength)) {
        menuSelections.Replace (loc,&newItem,true);
        _SetMenuItem (newItem,loc);
    }
}

//__________________________________________________________________

void            _HYPullDown::DeleteMenuItem  (long loc)
{
    if ((loc>=0)&&(loc<menuSelections.lLength)) {
        menuSelections.Delete (loc);
        _DeleteMenuItem (loc);
    }
}

//__________________________________________________________________

void            _HYPullDown::DeleteAllItems  (void)
{
    for (long k=menuSelections.lLength-1; k>=0; k--) {
        _DeleteMenuItem (k);
    }

    menuSelections.Clear ();
}

//__________________________________________________________________
_String*        _HYPullDown::GetMenuItem  (long loc)
{
    if ((loc>=0)&&(loc<menuSelections.lLength)) {
        return (_String*)menuSelections (loc);
    }
    return nil;
}

//__________________________________________________________________
long        _HYPullDown::GetSelection (void)
{
    return _GetSelection ();
}

//__________________________________________________________________

void            _HYPullDown::Activate  (void)
{
    _HYComponent::Activate();
    if (IsEnabled()) {
        _MarkForUpdate();
    }
}

//__________________________________________________________________

void            _HYPullDown::Deactivate  (void)
{
    _HYComponent::Deactivate();
    if (IsEnabled()) {
        _MarkForUpdate();
    }
}



//__________________________________________________________________
void        _HYPullDown::SendSelectionChange (void)
{
    if (messageRecipient) {
        messageRecipient->ProcessEvent (generateMenuSelChangeEvent (GetID(),GetSelection()));
    }
}

//__________________________________________________________________
void        _HYPullDown::EnableMenu      (bool flag)
{
    if (flag!=enabledFlag) {
        enabledFlag = flag;
        _EnableMenu    (flag);
        _MarkForUpdate ();
    }
}

//__________________________________________________________________
bool        _HYPullDown::IsEnabled   (void)
{
    return enabledFlag;
}

//__________________________________________________________________
void        _HYPullDown::SetVisibleSize  (_HYRect rel)
{
    _HYComponent::SetVisibleSize (rel);
    _HYPlatformPullDown::_SetVisibleSize (rel);
}

//__________________________________________________________________
void        _HYPullDown::ChangeSelection     (long newSel, bool eventSend)
{
    if ((newSel>=0)&&(newSel<menuSelections.lLength)) {
        selection = newSel;
        if (eventSend) {
            SendSelectionChange();
        }
        _MarkForUpdate();
    }
}

//__________________________________________________________________
void        _HYPullDown::EnableItem  (long theItem, bool toggle)
{
    if ((theItem>=0)&&(theItem<menuSelections.lLength)) {
        _EnableItem (theItem,toggle);
    }
}

//__________________________________________________________________
bool        _HYPullDown::ProcessEvent (_HYEvent* e)
{
    DeleteObject (e);
    return false;
}
