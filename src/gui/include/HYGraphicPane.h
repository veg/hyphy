/*
    A painting canvas with double buffer.

    Sergei L. Kosakovsky Pond, May 2000.
*/

#ifndef _HYGRAPHICPANE_
#define _HYGRAPHICPANE_
//#pragma once
#include "HYBaseGUI.h"
#include "HYPlatformGraphicPane.h"

#define  HY_FONT_PLAIN  0
#define  HY_FONT_BOLD   1
#define  HY_FONT_ITALIC 2

//__________________________________________________________________

class _HYGraphicPane: public _HYPlatformGraphicPane
{

public:

    _HYGraphicPane(int, int, int);
    // initial size, pixel depth

    virtual ~_HYGraphicPane();


    virtual void        SetPaneSize  (int,int,int);
    virtual void        EraseAll     (void);
    virtual void        DrawLine     (_HYRect);
    virtual void        DrawHatchedLine(_HYRect);
    // from, to , width
    virtual void        DisplayText  (_String,int,int,bool);
    // text, where, left-right or top-bottom
    virtual void        DisplayText  (_String, _HYRect, char);
    // text wrapped in a rect
    // with alignment

    virtual void        DisplayChar  (char,int,int);
    // text, where
    virtual void        DrawRect     (_HYRect);
    virtual void        DrawPicRes   (_HYRect& r,long i) {
        _DrawPicRes (r,i);
    }

    virtual void        FillRect     (_HYRect);
    virtual void        InvertRect   (_HYRect);
    virtual void        EraseRect    (_HYRect);
    virtual void        DrawOval     (_HYRect);
    virtual void        FillOval     (_HYRect);
    virtual void        EraseOval    (_HYRect);
    virtual void        DrawPolygon  (_SimpleList&, long = 1);
    virtual void        FillPolygon  (_SimpleList&);
    virtual void        ErasePolygon (_SimpleList&);
    virtual void        DrawArc      (_HYRect,int,int);
    virtual void        FillArc      (_HYRect,int,int);
    virtual void        EraseArc     (_HYRect,int,int);
    virtual void        SetColor     (_HYColor);
    virtual void        SetBColor    (_HYColor);
    virtual void        SetFont      (_HYFont);
    virtual void        SetFontSize  (long);
    virtual void        SetDialogBG  (void) {
        _SetDialogBG();
    }
    virtual void        DrawFrame    (_HYColor c);
    virtual _HYColor&   GetColor     (void) {
        return color;
    }
    virtual _HYColor&   GetBColor    (void) {
        return bColor;
    }
    virtual _HYFont&    GetFont      (void) {
        return font;
    }
    virtual void        StartDraw    (void);
    virtual void        EndDraw      (void);
    virtual void        DrawInfoBox  (_HYRect, _String);
    virtual void        DrawTable    (_List&, _List&, char, _HYRect, _HYRect);
    // strings,
    // placements obtained by ComputeTableCellPlacements,
    // alignment flags
    // _HYRect has the margins with width specifying the width of table lines
    // top - vOffset, left - hOffset, right - table with, bottom - table height
    virtual void        DrawHashes   (_HYRect, _SimpleList&, _List&, int, int);
    // HYRect has
    // left and right or top and bottom : the line to hash
    // width has the offset (to the left or to the bottom) of labels vs the axis
    // int is the length of the hashes
    // 2nd int is the width of the hashes


    int         w,h,depth;
    _HYColor    color,
                bColor;
    _HYFont     font;
};

extern  _HYColor hyDefaultFrameColor;
void    RotateRect90 (_HYRect&);
double  DistanceBetweenPoints     (double,double,double,double);
long    AngleBetweenPoints        (double,double,double,double);
long    ComputeTableCellPlacement (_List&, _List&, _HYRect, _HYFont&);
//      _HYRect contains the margins for each cell (left, right, etc)
//      the 2nd _List will contain _SimpleLists per row with 3 entries per cell:
//      top, left, horizontal extra space (used for alignment)
//      returns table width

long    ComputeHashMarkPlacement  (_HYRect, _Parameter, _Parameter, _Parameter&, _Parameter&, _SimpleList&, _List&, int, _HYFont&);
//      _HYRect contains the endpoints of the axis (if left or right != 0, then horizontal)
//      otherwise - vertical. _SimpleList will contain hashmark offsets
//      _Parameter& will contain the step;
//      _Parameter& #2 will contain the first (smallest hash)
//      _List& contains the labels
//      int is used to specify the minimal number of hashes
//      returns the width/height of the labels
#endif

//EOF
