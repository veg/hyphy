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

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

_String  _hyScrollingEvent       ("Scroll");
_String  _hyMenuSelChangeEvent   ("MenuChange");
_String  _hyMenuOpenEvent        ("MenuOpen");
_String  _hyButtonPushEvent      ("ButtonPush");
_String  _hyRebuildSCanvasEvent  ("CanvasStretch");
_String  _hyListChangeEvent      ("ListChange");
_String  _hyListDblClickEvent    ("ListDClick");
_String  _hyTableDblClickEvent   ("TableDClick");
_String  _hyKeyboardFocusEvent   ("KeyboardFocus");
_String  _hyTablePullDownEvent   ("TablePD");
_String  _hyTableEditCellEvent   ("TableEdit");
_String  _hyTableChangeSelEvent  ("TableChangeSel");
_String  _hyTableResizeCEvent    ("TableResize");
_String  _hyGlobalLFKillEvent    ("KillLF");
_String  _hyGlobalLFSpawnEvent   ("SpawnLF");
_String  _hyGlobalDFKillEvent    ("KillDF");
_String  _hyGlobalDSKillEvent    ("KillDS");
_String  _hyGlobalTreeKillEvent  ("KillTree");
_String  _hyGlobalChangeEvent    ("Change");
_String  _hyContextPopUp         ("CPopUp");
_String  _hyTextEditChange       ("TEChange");
_String  _hyGlobalCloseWindow    ("CW");
_String  _hyGlobalChangeLF       ("CLF");
_String  _hyGlobalChangeLFParams ("CLFP");
_String  _hyGlobalSetTreePanelSelection
("STPS");

//____________________________________________________

_HYEvent*   generateTablePullDownEvent (long objID, long index, long location)
{
    _String eventClass (_hyTablePullDownEvent),
            eventData  (objID);
    eventData = eventData & ','& _String (index) & ',' & _String (location);
    return  new _HYEvent (eventClass,eventData);
}

//____________________________________________________

_HYEvent*   generateTableChangeSelEvent (long objID)
{
    _String eventClass (_hyTableChangeSelEvent),
            eventData  (objID);
    return  new _HYEvent (eventClass,eventData);
}

//____________________________________________________

_HYEvent*   generateTableEditCellEvent (long objID, long index)
{
    _String eventClass (_hyTableEditCellEvent),
            eventData  (objID);
    eventData = eventData & ','& _String (index);
    return  new _HYEvent (eventClass,eventData);
}

//____________________________________________________

_HYEvent*   generateTableResizeCEvent (long objID, long index, long shift)
{
    _String eventClass (_hyTableResizeCEvent),
            eventData  (objID);
    eventData = eventData & ','& _String (index) & ','& _String (shift);
    return  new _HYEvent (eventClass,eventData);
}

//____________________________________________________

_HYEvent*   generateScrollEvent (int h,int v,long objID)
{
    _String eventClass (_hyScrollingEvent),
            eventData ((long)h);
    eventData = eventData & ','& _String ((long)v);
    if (objID>=0) {
        eventData = _String (objID)&','&eventData;
    }

    return  new _HYEvent (eventClass,eventData);
}

//____________________________________________________

_HYEvent*   generateMenuSelChangeEvent (long objID, long newSel)
{
    _String eventClass (_hyMenuSelChangeEvent),
            eventData (objID);
    eventData = eventData & ','& _String (newSel);
    return  new _HYEvent (eventClass,eventData);
}

//____________________________________________________

_HYEvent*   generateMenuOpenEvent (long objID)
{
    _String eventClass (_hyMenuOpenEvent),
            eventData (objID);
    return  new _HYEvent (eventClass,eventData);
}

//____________________________________________________

_HYEvent*   generateListDblClickEvent (long objID, long newSel)
{
    _String eventClass (_hyListDblClickEvent),
            eventData (objID);
    eventData = eventData & ','& _String (newSel);
    return  new _HYEvent (eventClass,eventData);
}

//____________________________________________________

_HYEvent*   generateTableDblClickEvent (long objID)
{
    _String eventClass (_hyTableDblClickEvent),
            eventData (objID);
    return  new _HYEvent (eventClass,eventData);
}

//____________________________________________________

_HYEvent*   generateRebuildSCanvas (long objID)
{
    _String eventClass (_hyRebuildSCanvasEvent),
            eventData (objID);
    return  new _HYEvent (eventClass,eventData);
}

//____________________________________________________

_HYEvent*   generateButtonPushEvent (long objID, long newSel)
{
    _String eventClass (_hyButtonPushEvent),
            eventData (objID);
    eventData = eventData & ','& _String (newSel);
    return  new _HYEvent (eventClass,eventData);
}

//____________________________________________________

_HYEvent*   generateListChangeEvent (long objID)
{
    _String eventClass (_hyListChangeEvent),
            eventData (objID);
    return  new _HYEvent (eventClass,eventData);
}

//____________________________________________________

_HYEvent*   generateKeyboardFocusEvent (long objID)
{
    _String eventClass (_hyKeyboardFocusEvent),
            eventData (objID);
    return  new _HYEvent (eventClass,eventData);
}

//____________________________________________________

void    postLFKillEvent (long objID,long lfID)
{
    _String eventClass (_hyGlobalLFKillEvent),
            eventData (objID);
    eventData = eventData & ','& _String (lfID);
    _HYEvent*
    killEvent = new _HYEvent (eventClass,eventData);

    GlobalGUIEventQueue << killEvent;
}


//____________________________________________________

void    postLFSpawnEvent (long objID,long lfID)
{
    _String eventClass (_hyGlobalLFSpawnEvent),
            eventData (objID);
    eventData = eventData & ','& _String (lfID);
    _HYEvent*
    spawnEvent = new _HYEvent (eventClass,eventData);

    GlobalGUIEventQueue << spawnEvent;
}

//____________________________________________________

void    postDSKillEvent (long objID,long lfID)
{
    _String eventClass (_hyGlobalDSKillEvent),
            eventData (objID);
    eventData = eventData & ','& _String (lfID);
    _HYEvent*
    killEvent = new _HYEvent (eventClass,eventData);

    GlobalGUIEventQueue << killEvent;
}

//____________________________________________________

void    postTreeKillEvent (long objID,long lfID)
{
    _String eventClass (_hyGlobalTreeKillEvent),
            eventData (objID);
    eventData = eventData & ','& _String (lfID);
    _HYEvent*
    killEvent = new _HYEvent (eventClass,eventData);

    GlobalGUIEventQueue << killEvent;
}

//____________________________________________________

void    postDFKillEvent (long objID,long lfID)
{
    _String eventClass (_hyGlobalDFKillEvent),
            eventData (objID);
    eventData = eventData & ','& _String (lfID);
    _HYEvent*
    killEvent = new _HYEvent (eventClass,eventData);

    GlobalGUIEventQueue << killEvent;
}

//____________________________________________________

void    postObjectChangelEvent (long senderID,long objectKind, long objectID)
{
    _String eventClass (_hyGlobalChangeEvent),
            eventData (senderID);
    eventData = eventData & ','& _String ((objectKind<<28)+objectID);
    _HYEvent*
    killEvent = new _HYEvent (eventClass,eventData);

    GlobalGUIEventQueue << killEvent;
}

//____________________________________________________

_HYEvent*   generateContextPopUpEvent (long objID, long locH, long locV)
{
    _String eventClass (_hyContextPopUp),
            eventData  (objID);
    eventData = eventData & ','& _String (locH) & ',' & _String (locV);
    return  new _HYEvent (eventClass,eventData);
}

//____________________________________________________

_HYEvent*   generateTextEditChangeEvent (long objID, long mouseOrKey)
{
    _String eventClass (_hyTextEditChange),
            eventData  (objID);
    eventData = eventData & ','& _String (mouseOrKey);
    return  new _HYEvent (eventClass,eventData);
}

//____________________________________________________

void    postWindowCloseEvent (long senderID)
{
    _String eventClass (_hyGlobalCloseWindow),
            eventData (senderID);
    _HYEvent*
    closeEvent = new _HYEvent (eventClass,eventData);

    GlobalGUIEventQueue << closeEvent;
}

//____________________________________________________

void    postChangeLFEvent (long senderID, long lfID)
{
    _String eventClass (_hyGlobalChangeLF),
            eventData (senderID);

    eventData = eventData & ',' & lfID;
    _HYEvent*
    modEvent = new _HYEvent (eventClass,eventData);

    GlobalGUIEventQueue << modEvent;
}

//____________________________________________________

void    postChangeLFParamsEvent (long senderID, long lfID)
{
    _String eventClass (_hyGlobalChangeLFParams),
            eventData (senderID);

    eventData = eventData & ',' & lfID;
    _HYEvent*
    modEvent = new _HYEvent (eventClass,eventData);

    GlobalGUIEventQueue << modEvent;
}

//____________________________________________________

void    postSetTreePanelSelection (long treeID, _String* sel)
{
    _String eventClass (_hyGlobalSetTreePanelSelection),
            eventData;

    eventData = _String(treeID) & ',' & *sel;
    _HYEvent*
    modEvent = new _HYEvent (eventClass,eventData);

    GlobalGUIEventQueue << modEvent;
}

//____________________________________________________
//  Update/Notification structures
//____________________________________________________

//_List     notificationPairs;

//____________________________________________________

//EOF