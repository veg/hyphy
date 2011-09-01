/*
    Event Types and Generators

    Sergei L. Kosakovsky Pond, May-August 2000.
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