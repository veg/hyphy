/*
    A sequence display component.

    Sergei L. Kosakovsky Pond, August 2000.
*/

#ifndef _HYSEQPANEL_
#define _HYSEQPANEL_
//#pragma once
#include "HYCanvas.h"

#define  HY_SEQUENCE_PANE_CHAR_SPACING 4
#define  HY_SEQUENCE_PANE_TOP_MARGIN   15
#define  HY_SEQUENCE_PANE_NAMES_NONE   0
#define  HY_SEQUENCE_PANE_NAMES_SHORT  1
#define  HY_SEQUENCE_PANE_NAMES_ALL    2
#define  HY_SEQUENCE_PANE_NAMES_MASK   0x03
#define  HY_SEQUENCE_PANE_SHORT_WIDTH  10

//__________________________________________________________________

class _HYSequencePane: public _HYCanvas
{

public:

    _HYSequencePane(_HYRect,Ptr,int,int);

    virtual ~_HYSequencePane() {};

    virtual void        _Paint  (Ptr);
    virtual void        BuildPane   (bool setport = true);
    virtual void        BuildHeaders(void);

    virtual bool        ProcessEvent (_HYEvent* e);

    virtual void        InsertColumn (_String*, long, bool update = false);
    virtual void        RemoveColumn (long, bool update = false);
    virtual void        SetCharColor (unsigned char, _HYColor, bool update = false);
    virtual int         GetMaxW (void);
    virtual int         GetMaxH (void);
    virtual int         GetMinH (void);
    virtual void        SetFont                 (_HYFont);
    virtual void        SetBackColor            (_HYColor);
    virtual void        SetHeaderColor          (_HYColor);
    virtual void        SetHighliteColor        (_HYColor);
    virtual void        SetVisibleSize          (_HYRect);

    virtual void        HScrollPane             (long dx);
    virtual long        RowCount                (void);
    virtual void        SetHeaders              (_List*,bool);
    virtual void        SetBlockWidth           (long,bool);
    virtual void        AddColumnToSelection    (long);
    virtual void        AddSpeciesToSelection   (long);
    virtual void        SelectAll               (bool);
    virtual void        SelectRange             (_SimpleList&, bool = false);
    virtual void        ProcessSelectionChange  (long, long, bool, bool, bool drag = false);
    virtual void        ProcessVSelectionChange (long, long, bool, bool, bool drag = false, bool editName = false);
    virtual void        SendSelectionChange     (bool vert = false);
    virtual void        ChangeHeaders           (void);
    virtual void        MoveSpecies             (long, long);
    virtual void        AlphabetizeSpecies      (void);
    virtual void        RevertFileOrder         (void);
    virtual void        ProcessContextualPopUp  (long,long);
    virtual void        SetSequenceOrder        (_SimpleList&);
    virtual void        SelectSequenceNames     (_List&, bool = true);


    virtual bool        _ProcessOSEvent         (Ptr);
    virtual void        SetNameDisplayMode      (unsigned char, bool update = false);
    virtual long        GetSlotHeight           ();
    virtual void        SetActiveOrPassive      (bool ap) {
        active = ap;
    }
    virtual void        SetNumberDisplay        (bool nd) {
        numbers = nd;
    }
    virtual void        ScrollToHSelection      (void);
    virtual void        CleanUpSequenceNames    (void);
    virtual void        EditSequenceName        (long);
    virtual void        BatchRenameSequences    (_List&,_List&);

    // data members

    _List           columnStrings,
                    rowHeaders;

    _SimpleList     hiliteColumns,
                    characterColors,
                    selection,
                    vselection,
                    speciesIndex;

    long            startColumn,
                    endColumn,
                    startRow,
                    endRow,
                    headerWidth,
                    charWidth,
                    blockWidth,
                    recentClick,
                    shortHeaderWidth,
                    fullHeaderWidth,
                    dragScrollCounter;

    _HYColor        backColor,
                    headerColor,
                    highlightColor;

    unsigned char   colorMap[256], nameDisplayFlags;

    bool            showDots,
                    active,
                    numbers,
                    invertColors;

};

//__________________________________________________________________

#endif

//EOF
