/*
    Event Types and Generators

    Sergei L. Kosakovsky Pond, May 2000.
*/

#ifndef _HYEVENTTYPES_
#define _HYEVENTTYPES_

//#pragma once
#include "HYBaseGUI.h"


_HYEvent*   generateScrollEvent         (int,int,long objID = -1);
_HYEvent*   generateMenuSelChangeEvent  (long,long);
_HYEvent*   generateMenuOpenEvent       (long);
_HYEvent*   generateButtonPushEvent     (long,long);
_HYEvent*   generateRebuildSCanvas      (long);
_HYEvent*   generateListChangeEvent     (long);
_HYEvent*   generateListDblClickEvent   (long,long);
_HYEvent*   generateTableDblClickEvent  (long);
_HYEvent*   generateKeyboardFocusEvent  (long);
_HYEvent*   generateTablePullDownEvent  (long,long,long);
_HYEvent*   generateTableEditCellEvent  (long,long);
_HYEvent*   generateTableResizeCEvent   (long,long,long);
_HYEvent*   generateTableChangeSelEvent (long);
_HYEvent*   generateContextPopUpEvent   (long,long,long);
_HYEvent*   generateTextEditChangeEvent (long,long);

extern _String _hyScrollingEvent,
       _hyMenuSelChangeEvent,
       _hyMenuOpenEvent,
       _hyButtonPushEvent,
       _hyRebuildSCanvasEvent,
       _hyListChangeEvent,
       _hyTableDblClickEvent,
       _hyKeyboardFocusEvent,
       _hyListDblClickEvent,
       _hyTablePullDownEvent,
       _hyTableEditCellEvent,
       _hyTableChangeSelEvent,
       _hyTableResizeCEvent,
       _hyContextPopUp,
       _hyTextEditChange,
       _hyGlobalLFKillEvent,
       _hyGlobalLFSpawnEvent,
       _hyGlobalDFKillEvent,
       _hyGlobalTreeKillEvent,
       _hyGlobalDSKillEvent,
       _hyGlobalDFAlterEvent,
       _hyGlobalCloseWindow,
       _hyGlobalChangeLF,
       _hyGlobalChangeLFParams,
       _hyGlobalSetTreePanelSelection;


void    postLFKillEvent             (long,long);
void    postLFSpawnEvent            (long,long);
void    postDFKillEvent             (long,long);
void    postDSKillEvent             (long,long);
void    postTreeKillEvent           (long,long);
void    postObjectChangelEvent      (long,long,long);
void    postWindowCloseEvent        (long);
void    postChangeLFEvent           (long,long);
void    postChangeLFParamsEvent     (long,long);
void    postSetTreePanelSelection   (long,_String*);


void    HandleGlobalQueueEvent      (void);


/*extern _String _hyNotifyNewTree,
               _hyNotifyChangeTree


void        RegisterForNotification   (_HYGuiObject*,_String&,_String&);
// request that the object be notified of events of a certain type
void        UnregisterFromNotification(_HYGuiObject*);
// remove  the object from notification queues
void        PostNotification          (_HYGuiObject*,_String&,_String&);


void
*/
#endif

//EOF
