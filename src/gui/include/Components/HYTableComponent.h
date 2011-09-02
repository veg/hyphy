/*
    A table display component.

    Sergei L. Kosakovsky Pond, January 2001.
*/

#ifndef _HYTABLE_
#define _HYTABLE_
//#pragma once

#include "HYComponent.h"
#include "HYPlatformTable.h"
#include "batchlan.h"

#define  HY_TABLE_STATIC_TEXT       0x0001
#define  HY_TABLE_EDIT_TEXT         0x0002
#define  HY_TABLE_ICON              0x0004
#define  HY_TABLE_BEVELED           0x0008
#define  HY_TABLE_SELECTED          0x0010
#define  HY_TABLE_BOLD              0x0020
#define  HY_TABLE_ITALIC            0x0040
#define  HY_TABLE_PULLDOWN          0x0080
#define  HY_TABLE_CANTSELECT        0x0100

#define  HY_TABLE_STYLEMASK         0x60
#define  HY_TABLE_DESELECT          0xFFEF

#define  HY_TABLE_SEL_ROWS          0x01
#define  HY_TABLE_SEL_COLS          0x02
#define  HY_TABLE_DONT_SIZE         0x04
#define  HY_TABLE_DONT_GROW_VERT    0x08
#define  HY_TABLE_DONT_GROW_HORIZ   0x10
#define  HY_TABLE_NO_COLS_LINES     0x20
#define  HY_TABLE_SINGLE_SELECTION  0x40
#define  HY_TABLE_FOCUSABLE         0x80
#define  HY_TABLE_IS_FOCUSED        0x100
#define  HY_TABLE_NODRAG_SELECTION  0x200
#define  HY_TABLE_HORZ_STRETCH      0x400
#define  HY_TABLE_VERT_STRETCH      0x800

#define  HY_TABLE_COLOR_BOX         0x01
#define  HY_TABLE_COLOR_CIRCLE      0x02

//__________________________________________________________________

class _HYTable: public _HYComponent, public _HYPlatformTable
{

public:
    _HYTable        (_HYRect,Ptr, long,long,long,long,long);
    // rel rect, window data, rows, columns, h step, v step, default cell type
    virtual             ~_HYTable       ();

    virtual void        _Paint          (Ptr);
    virtual bool        _ProcessOSEvent (Ptr);
    void        _HScrollTable   (long);
    void        _VScrollTable   (long);
    void        _MarkCellsForUpdate
    (_SimpleList&);
    void        _MarkCellForUpdate
    (long);
    void        _MarkColumnForUpdate
    (long);
    void        _MarkRowForUpdate
    (long);
    long        _HandlePullDown (_List&,long,long,long);
    void        _ComponentMouseExit (void);

    virtual bool        ProcessEvent    (_HYEvent*);
    virtual int         GetMaxW         (void);
    virtual int         GetMaxH         (void);
    virtual int         GetMaxLW        (void);
    virtual int         GetMaxLH        (void);
    virtual void        SetVisibleSize  (_HYRect);
    virtual void        SetTableFromMx  (_Matrix*, long, long, long);

    virtual void        SetFont         (_HYFont&);
    virtual void        SetBackColor    (_HYColor);
    virtual void        SetBackColor2   (_HYColor);
    virtual void        SetTextColor    (_HYColor);
    virtual void        SetRowSelection (const _SimpleList&);
    virtual void        SetColumnSelection (const _SimpleList&);
    virtual void        _ScrollRowIntoView
    (long);
    virtual void        ClearSelection  (bool  standAlone = true);

    virtual BaseRef     GetCellData     (long,long);
    virtual void        SetCellData     (BaseRef,long,long,long,bool);
    virtual void        SetTableSize    (long,long,long,long,long);
    virtual void        AddRow          (long,long,long);
    virtual  void       RequestSpace    (long,long);
    virtual void        DeleteRow       (long);
    virtual void        AddColumn       (long,long,long);
    virtual void        DeleteColumn    (long);
    virtual void        ClearTable      (bool all = false);
    virtual void        GetDisplayRange (_HYRect*,long&,long&,long&,long&);

    void        EnforceWidth    (long,long,bool = false);
    void        EnforceHeight   (long,long,bool = false);
    void        GetSelection    (_SimpleList&);
    void        SetSelection    (_SimpleList&, bool update = false);
    void        GetRowSelection (_SimpleList&, long = 0);
    void        ScanRowSelection(_SimpleList&);
    void        ScanColumnSelection
    (_SimpleList&);
    long        GetFirstRowSelection
    (void);
    bool        IsRowSelectionSimple
    (void);
    void        GetColumnSelection
    (_SimpleList&);
    void        SetRowOrder     (_SimpleList&);

    void        SetRowSpacing   (long,long,bool);
    void        SetColumnSpacing(long,long,bool);
    void        AutoFitColumn   (long,bool,bool increaseOnly = false);
    void        AutoFitWidth    (void);
    void        AutoFitWidth    (_HYTable&, long = 0);

    bool        ScrollToRow     (long);
    bool        ScrollToColumn  (long);

    void        InvertSelection (void);

    long        GetRowSpacing   (long);
    long        GetColumnSpacing(long);

    bool        CheckForHSizeLocation
    (long);
    long        FindClickedTableCell
    (long,long,long&,long&);

    void        DragRow         (long, long);

    virtual void        ModifySelection (long,long,bool,bool,bool message);
    void        ExpungeSelection(void);

    void        SaveTable       (_HYTable*, _HYTable*, long, FILE*,_String&,_SimpleList&,_SimpleList&);
    void        SaveTable       (_HYTable*, _HYTable*, long, FILE*,_String&);
    void        ExportCell      (FILE*,long,long,long);
    long        GetCellWidth    (long,long);

    static  void        GetTableFormats (_List&);
    long        FindString      (_String*, long startat = 0);

    void        EditBoxHandler  (long,_HYRect&);
    virtual void        UnfocusComponent(void);
    virtual void        FocusComponent  (void);
    virtual void        _FocusComponent (void);
    virtual void        IdleHandler     (void) {
        _IdleHandler();
    }

    virtual void        _IdleHandler    (void);
    long        _HandlePullDown (_List&,long,long);

    void        _ScrollVPixels  (long);
    void        _ScrollHPixels  (long);

    void        _PrintTable     (_SimpleList&, _SimpleList&,_HYTable* cHeaders = nil);
    void        _PrintTable     (_HYTable* cHeaders = nil);
    void        _PrintTable     (_SimpleList&,_HYTable* cHeaders = nil);

    virtual void        HandleKeyMove   (char, bool);

    virtual bool        CanCopy         (void);
    virtual BaseRef     CanPaste        (_String&);

    virtual _String*    HandleCopy      (void);
    virtual void        HandlePaste     (BaseRef);

//virtual void      _Activate           (void);
//virtual void      _Deactivate         (void);

    // data components

    _List           cellData;
    _SimpleList     cellTypes,
                    horizontalSpaces,
                    verticalSpaces;
    _HYColor        backColor,
                    backColor2,
                    textColor;

    _HYFont         textFont;
    short           editCellID,
                    selectionType;
    _String*        undoString;

    long            undoIndex,
                    undoIndex2,
                    stretchWidth,
                    stretchHeight;
};

//__________________________________________________________________

class _HYHList: public _HYTable
{

public:
    _HYHList        (_HYRect, Ptr, _List&, bool = true);
    // rel rect, window data, list data, single or multiple
    virtual             ~_HYHList       () {}
    virtual void        ModifySelection (long,long,bool,bool,bool);
    virtual long        IsARubrik       (long);
    virtual bool        IsRubrikOpen    (long);
    virtual long        RubrikIndex     (long);
    virtual long        AbsoluteIndex   (long,long);
    virtual void        FitToHeight     (long);
    virtual bool        HasPadding      (void);
    virtual void        HandleKeyMove   (char, bool);
    virtual long        FindString      (_String*, long startat = 0);
    virtual void        AddRubrik       (_String&, _List&, long);
    virtual void        DeleteRubrik    (long);
    virtual long        IsSingleRubrik  (_SimpleList&);

    _List           listData;
    _SimpleList     rubrikIndex;
    bool            single;
};

//__________________________________________________________________

#define   tPDMh         12
#define   tPDMw         13

extern    _HYColor      tableDefaultBk,
          tableDefaultBk2;

extern    _SimpleList*  openArrow,
          *    closedArrow;


#endif

//EOF
