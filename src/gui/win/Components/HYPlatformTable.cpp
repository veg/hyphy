/*
    Table component for Mac OS API

    Sergei L. Kosakovsky Pond, May 2000-December 2002
*/

#include "errorfns.h"
#include "HYTableComponent.h"
#include "HYUtils.h"
#include "HYEventTypes.h"
#include "HYWindow.h"
#include "HYGraphicPane.h"
#include "HYTextBox.h"
#include "HYTableWindow.h"
#include "HYDialogs.h"

//__________________________________________________________________

BOOL CALLBACK                   PrintDialogProc(HWND, UINT, WPARAM, LPARAM);
BOOL CALLBACK                   AbortProc(HDC, int);
HDC                             GetPrinterDeviceContext(HWND);

extern  BOOL                    UserAbortFlag;
extern  HWND                    PrintDialogHandle;


extern  HCURSOR                 hSizeCursor,
        pickUpCursor,
        dropOffCursor;

HPEN    menuLine1               = CreatePen (PS_SOLID,1,RGB(0xA0,0xA0,0xA0)),
        menuLine2                = CreatePen (PS_SOLID,1,RGB(0x04,0x04,0x04));

HBRUSH  _BLACKBRUSH_            = CreateSolidBrush (RGB(0,0,0));

extern  HBITMAP                 tablePDMenuIcon;


//__________________________________________________________________

static LRESULT  CALLBACK tableEditSubclassHandler (HWND WindowHand, UINT iMsg, WPARAM wParam, LPARAM lParam)
{
    _HYTable * theParent = (_HYTable*)GetWindowLongPtr (WindowHand, GWLP_USERDATA);

    switch (iMsg) {
    case WM_CHAR: {
        int     keyCode = wParam;
        if (keyCode==VK_RETURN) {
            theParent->EditBoxHandler (-1,theParent->rel);
            return true;
        } else if (keyCode==VK_TAB) {
            SendMessage (theParent->parentWindow, WM_CHAR, wParam, lParam);
            return true;
        }

        break;
    }
    }

    return CallWindowProc ((WNDPROC)theParent->defautlTextHandler, WindowHand, iMsg, wParam, lParam);
}


//__________________________________________________________________

_HYPlatformTable::_HYPlatformTable(void)
{
    backPattern  = nil;
    backPattern2 = nil;
    cursorState  = false;
    editBox      = nil;
    tableFont    = nil;
    activeColumn = -1;
    activeColumn2= -1;
}

//__________________________________________________________________

_HYPlatformTable::~_HYPlatformTable(void)
{
    if (backPattern) {
        DeleteObject (backPattern);
    }
    if (backPattern2) {
        DeleteObject (backPattern2);
    }
    if (tableFont) {
        DeleteObject (tableFont);
    }
    if (tableFontB) {
        DeleteObject (tableFontB);
    }
    if (tableFontI) {
        DeleteObject (tableFontI);
    }
    if (tableFontBI) {
        DeleteObject (tableFontBI);
    }
}

//__________________________________________________________________

void        _HYPlatformTable::_SetFont (void)
{
    _HYTable* parent = (_HYTable*)this;

    if (tableFont) {
        DeleteObject (tableFont);
    }
    if (tableFontB) {
        DeleteObject (tableFontB);
    }
    if (tableFontI) {
        DeleteObject (tableFontI);
    }
    if (tableFontBI) {
        DeleteObject (tableFontBI);
    }

    tableFont = CreateFont (parent->textFont.size,0,0,0,FW_NORMAL,FALSE,FALSE,FALSE,ANSI_CHARSET,OUT_DEFAULT_PRECIS,
                            CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,DEFAULT_PITCH|FF_DONTCARE,parent->textFont.face.sData);

    checkPointer (tableFont);

    tableFontB = CreateFont (parent->textFont.size,0,0,0,FW_BOLD,FALSE,FALSE,FALSE,ANSI_CHARSET,OUT_DEFAULT_PRECIS,
                             CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,DEFAULT_PITCH|FF_DONTCARE,parent->textFont.face.sData);

    checkPointer (tableFontB);

    tableFontI = CreateFont (parent->textFont.size,0,0,0,FW_NORMAL,TRUE,FALSE,FALSE,ANSI_CHARSET,OUT_DEFAULT_PRECIS,
                             CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,DEFAULT_PITCH|FF_DONTCARE,parent->textFont.face.sData);

    checkPointer (tableFontI);

    tableFontBI = CreateFont (parent->textFont.size,0,0,0,FW_BOLD,true,FALSE,FALSE,ANSI_CHARSET,OUT_DEFAULT_PRECIS,
                              CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,DEFAULT_PITCH|FF_DONTCARE,parent->textFont.face.sData);

    checkPointer (tableFontBI);
}

//__________________________________________________________________
void        _HYTable::_HScrollTable (long h)
{
    if (h) {
        long        vsShift = ((settings.width&HY_COMPONENT_H_SCROLL)?HY_SCROLLER_WIDTH:0),
                    hsShift = ((settings.width&HY_COMPONENT_V_SCROLL)?HY_SCROLLER_WIDTH:0);

        EditBoxHandler  (-1,rel);

        _HYRect         paintRect = rel;
        paintRect.width = (long)GetDC (parentWindow);

        if (abs(h)>(rel.right-rel.left)/2) {
            _Paint((Ptr)&paintRect);
        } else {
            RECT            scrollRect;
            scrollRect.top = rel.top;
            scrollRect.bottom = rel.bottom-vsShift;
            scrollRect.right = rel.right-hsShift;
            scrollRect.left  = rel.left;
            paintRect.top = scrollRect.top;
            paintRect.bottom = rel.bottom;
            if (h>0) {
                ScrollWindowEx (parentWindow, -h, 0, &scrollRect, &scrollRect, nil, nil, 0); //?
                paintRect.right = scrollRect.right+hsShift;
                paintRect.left = scrollRect.right-h;
                hOrigin+=paintRect.left-rel.left;
                Paint((Ptr)&paintRect);
                hOrigin-=paintRect.left-rel.left;
            } else {
                ScrollWindowEx (parentWindow, -h, 0, &scrollRect, &scrollRect, nil, nil, 0); //?
                paintRect.left = scrollRect.left;
                paintRect.right = paintRect.left-h+hsShift;
                Paint((Ptr)&paintRect);
            }
        }

        ReleaseDC (parentWindow, (HDC)paintRect.width);
    }
}

//__________________________________________________________________
void        _HYTable::_VScrollTable (long v)
{
    if (v) {
        long        vsShift = ((settings.width&HY_COMPONENT_H_SCROLL)?HY_SCROLLER_WIDTH:0),
                    hsShift = ((settings.width&HY_COMPONENT_V_SCROLL)?HY_SCROLLER_WIDTH:0);

        EditBoxHandler  (-1,rel);

        _HYRect         paintRect = rel;
        paintRect.width = (long)GetDC (parentWindow);

        if (abs(v)>(rel.bottom-rel.top)/2) {
            _Paint((Ptr)&paintRect);
        } else {
            RECT      scrollRect;

            scrollRect.left         = rel.left;
            scrollRect.right        = rel.right-hsShift;
            scrollRect.top          = rel.top;
            scrollRect.bottom       = rel.bottom-vsShift;
            paintRect.left          = scrollRect.left;
            paintRect.right         = rel.right;
            if (v>0) {
                ScrollWindowEx (parentWindow, 0, -v, &scrollRect, &scrollRect, nil, nil, 0); //?
                paintRect.top = rel.bottom-vsShift-v-1;
                paintRect.bottom = rel.bottom-1;
                vOrigin+=paintRect.top-rel.top;
                Paint((Ptr)&paintRect);
                vOrigin-=paintRect.top-rel.top;
            } else {
                scrollRect.top --;
                ScrollWindowEx (parentWindow, 0, -v, &scrollRect, &scrollRect, nil, nil, 0); //?
                paintRect.top = rel.top;
                paintRect.bottom = rel.top-v+vsShift+2;
                paintRect.bottom = rel.bottom;
                Paint((Ptr)&paintRect);
            }
        }
        ReleaseDC (parentWindow, (HDC)paintRect.width);
    }
}

//__________________________________________________________________

void        _HYTable::_ComponentMouseExit (void)
{
    if (cursorState) {
        if ((cursorState == HY_TABLE_DRAG_CURSOR)&&(activeColumn>=0)&&(activeColumn2>=0)) {
            _HiliteRowForDrag (activeColumn2,activeColumn);
        }

        _ResetCursorState ();
    }
}

//__________________________________________________________________

void        _HYPlatformTable::_SetBackColor (_HYColor& c)
{
    DeleteObject (backPattern);
    backPattern = CreateSolidBrush (RGB(c.R,c.G,c.B));
    checkPointer (backPattern);
}

//__________________________________________________________________

void        _HYPlatformTable::_SetBackColor2 (_HYColor& c)
{
    DeleteObject (backPattern2);
    backPattern2 = CreateSolidBrush (RGB(c.R,c.G,c.B));
    checkPointer (backPattern2);
}

//__________________________________________________________________

void        _HYPlatformTable::_KillTextBox (void)
{
    if (editBox) {
        DestroyWindow (editBox);
        editBox = nil;
        _ResetCursorState ();
    }
}

//__________________________________________________________________

_String     _HYPlatformTable::_RetrieveTextValue (void)
{
    return retrieveEditControlText(editBox);
}


//__________________________________________________________________

void        _HYPlatformTable::_CreateTextBox (_HYRect& tBox,_String& textIn)
{
    textBoxRect = HYRect2Rect(tBox);

    _HYTable * theParent = (_HYTable*)this;

    editBox = CreateWindow ("EDIT",textIn.sData,WS_VISIBLE|WS_CHILD|ES_LEFT| ES_AUTOHSCROLL|WS_BORDER,
                            textBoxRect.left,textBoxRect.top,textBoxRect.right-textBoxRect.left,textBoxRect.bottom-textBoxRect.top
                            ,theParent->parentWindow,NULL,ProgramInstance, NULL);

    checkPointer   (editBox);

    SetWindowLongPtr (editBox,GWLP_USERDATA, (LONG_PTR)theParent);
    defautlTextHandler  = SetWindowLongPtr (editBox,GWLP_WNDPROC,(LONG_PTR)tableEditSubclassHandler);

    if (tableFont) {
        SendMessage (editBox, WM_SETFONT, (WPARAM)tableFont, 1);
    }

    SendMessage (editBox, EM_SETMARGINS, EC_LEFTMARGIN, 3);
    SendMessage (editBox, EM_SETMARGINS, EC_RIGHTMARGIN, 3);
    SendMessage (editBox, EM_SETSEL, 0, 50000);
    SetFocus    (editBox);

}


//__________________________________________________________________

RECT        _HYPlatformTable::_GetVisibleRowRect (long h)
{
    _HYTable*   parent = (_HYTable*)this;
    RECT        res;

    long        w = (parent->settings.width&HY_COMPONENT_H_SCROLL)?HY_SCROLLER_WIDTH:0;
    res.left   = parent->rel.left;
    res.right  = (parent->settings.width&HY_COMPONENT_H_SCROLL)?parent->rel.right-HY_SCROLLER_WIDTH:parent->rel.right;

    res.bottom = parent->verticalSpaces.lData[h]-parent->vOrigin+parent->rel.top;
    if (res.bottom>parent->rel.bottom-w) {
        res.bottom=parent->rel.bottom-w;
    }

    if (h) {
        res.top = parent->verticalSpaces.lData[h-1]-parent->vOrigin+parent->rel.top;
    } else {
        res.top = parent->rel.top-parent->vOrigin;
    }

    return res;
}


//__________________________________________________________________

void        _HYPlatformTable::_HiliteRowForDrag (long row, long old)
{
    RECT        cellRect = _GetVisibleRowRect (row);

    if (row<old) {
        cellRect.bottom = cellRect.top+2;
    } else {
        cellRect.top = cellRect.bottom-2;
    }

    HDC         winDC = GetDC (((_HYTable*)this)->parentWindow);

    InvertRect  (winDC, &cellRect);

    ReleaseDC   (((_HYTable*)this)->parentWindow,winDC);
}

//__________________________________________________________________

void        _HYTable::_MarkCellsForUpdate (_SimpleList& cells)
{
    long hs,hf,vs,vf,t,t2,k;
    GetDisplayRange (&rel, hs, hf, vs, vf);
    RECT  clipRect;

    clipRect.left = rel.left;
    clipRect.right = rel.right;
    if (settings.width&HY_COMPONENT_V_SCROLL) {
        clipRect.right -= HY_SCROLLER_WIDTH;
    }
    clipRect.top = rel.top;
    clipRect.bottom = rel.bottom;
    if (settings.width&HY_COMPONENT_H_SCROLL) {
        clipRect.bottom -= HY_SCROLLER_WIDTH;
    }

    for (k=0; k<cells.lLength; k++) {
        t2 = cells.lData[k]/horizontalSpaces.lLength;
        t = cells.lData[k]%horizontalSpaces.lLength;

        if ((t>=hs)&&(t<=hf)&&(t2<=vf)&&(t2>=vs)) {
            //printf ("%d %d\n",t,t2);
            RECT invalRect;
            if (t) {
                invalRect.left = horizontalSpaces.lData[t-1];
            } else {
                invalRect.left = 0;
            }
            invalRect.right = horizontalSpaces.lData[t];
            if (t2) {
                invalRect.top = verticalSpaces.lData[t2-1];
            } else {
                invalRect.top = 0;
            }
            invalRect.bottom = verticalSpaces.lData[t2];

            OffsetRect (&invalRect, rel.left-hOrigin, rel.top-vOrigin);

            RECT           tempRect;
            IntersectRect (&tempRect, &clipRect, &invalRect); //?
            if (!IsRectEmpty (&tempRect)) {
                InvalidateRect (parentWindow,&tempRect,false);
            }
        }
    }
}

//__________________________________________________________________

void        _HYTable::_MarkColumnForUpdate (long index)
{
    long hs,hf,vs,vf;
    GetDisplayRange (&rel, hs, hf, vs, vf);

    if ((index>=hs)&&(index<=hf)) {
        RECT  clipRect;

        clipRect.left = rel.left;
        clipRect.right = rel.right;

        if (settings.width&HY_COMPONENT_V_SCROLL) {
            clipRect.right -= HY_SCROLLER_WIDTH;
        }
        clipRect.top = rel.top;
        clipRect.bottom = rel.bottom;
        if (settings.width&HY_COMPONENT_H_SCROLL) {
            clipRect.bottom -= HY_SCROLLER_WIDTH;
        }

        RECT  invalRect;

        invalRect.bottom = clipRect.bottom;
        invalRect.top = clipRect.top;
        invalRect.left = index?horizontalSpaces.lData[index-1]:0;
        invalRect.right = horizontalSpaces.lData[index];

        OffsetRect (&invalRect, rel.left-hOrigin, 0);

        RECT           tempRect;
        IntersectRect (&tempRect, &clipRect, &invalRect); //?
        if (!IsRectEmpty (&tempRect)) {
            InvalidateRect (parentWindow,&tempRect,false);
        }
    }
}


//__________________________________________________________________

void        _HYTable::_MarkRowForUpdate (long index)
{
    long hs,hf,vs,vf;
    GetDisplayRange (&rel, hs, hf, vs, vf);

    if ((index>=vs)&&(index<=vf)) {
        RECT  clipRect;

        clipRect.left = rel.left;
        clipRect.right = rel.right;
        if (settings.width&HY_COMPONENT_V_SCROLL) {
            clipRect.right -= HY_SCROLLER_WIDTH;
        }
        clipRect.top = rel.top;
        clipRect.bottom = rel.bottom;
        if (settings.width&HY_COMPONENT_H_SCROLL) {
            clipRect.bottom -= HY_SCROLLER_WIDTH;
        }

        RECT  invalRect;
        invalRect.right = clipRect.right;
        invalRect.left = clipRect.left;
        invalRect.top = index?verticalSpaces.lData[index-1]:0;
        invalRect.bottom = verticalSpaces.lData[index];
        OffsetRect (&invalRect, 0, rel.top-vOrigin);

        RECT           tempRect;
        IntersectRect (&tempRect, &clipRect, &invalRect); //?
        if (!IsRectEmpty (&tempRect)) {
            InvalidateRect (parentWindow,&tempRect,false);
        }
    }
}


//__________________________________________________________________

void        _HYTable::_MarkCellForUpdate (long index)
{
    _SimpleList     dummy (index);
    _MarkCellsForUpdate (dummy);
}

//__________________________________________________________________

void        _HYTable::_IdleHandler (void)
{
}

//__________________________________________________________________

void        _HYTable::_FocusComponent (void)
{
    if (GetFocus()!=parentWindow) {
        //printf ("Table SetFocus\n");
        SetFocus (parentWindow);
    }
}


//__________________________________________________________________
long        _HYTable::_HandlePullDown (_List& data, long h, long v, long currentS)
{
    if (data.lLength) {
        return HandlePullDownWithFont (data,h,v,currentS,textFont.face,textFont.size);
    }
    return -1;
}


//__________________________________________________________________

void        _HYTable::_ScrollVPixels (long offset)
{
    long     voff = ((settings.width&HY_COMPONENT_H_SCROLL)?HY_SCROLLER_WIDTH:0);
    offset = offset/(_Parameter)(GetMaxH()- rel.bottom+rel.top+1-voff)*MAX_CONTROL_VALUE;
    ProcessEvent (generateScrollEvent(0,offset));
    _SetVScrollerPos((double)MAX_CONTROL_VALUE*vOrigin/(verticalSpaces.lData[verticalSpaces.lLength-1]-vSize+voff));
}

//__________________________________________________________________

void        _HYTable::_ScrollHPixels (long offset)
{
    long     hoff = ((settings.width&HY_COMPONENT_V_SCROLL)?HY_SCROLLER_WIDTH:0);
    offset = offset/(_Parameter)(GetMaxW()- rel.right+rel.bottom+1-hoff)*MAX_CONTROL_VALUE;
    ProcessEvent (generateScrollEvent(offset,0));
    _SetHScrollerPos((double)MAX_CONTROL_VALUE*hOrigin/(horizontalSpaces.lData[horizontalSpaces.lLength-1]-hSize+hoff));
}


//__________________________________________________________________

void        _HYTable::_ScrollRowIntoView (long index)
{
    if ((index>=0)&&(index<verticalSpaces.lLength)) {
        long hs, hf, vs, vf;
        GetDisplayRange (&rel,hs,hf,vs,vf);
        if ((index>vf)||(index<vs)) {
            _ScrollVPixels ((index?verticalSpaces.lData[index-1]:0)-vOrigin);
        }
    }
}

//__________________________________________________________________
void        _HYTable::_Paint (Ptr p)
{

    _HYRect         *relRect    = (_HYRect*)p;

    bool            isPrinting = relRect->right<0;

    if (isPrinting) {
        relRect->right = -relRect->right;
    }

    HDC             dc = (HDC)relRect->width;

    /*RECT          all = HYRect2Rect (*relRect);

    HBRUSH          redFill   = CreateSolidBrush (RGB(255,0,0)); // checked
    FillRect        (dc, &all, redFill);
    DeleteObject    (redFill);

    return;
    */
    COLORREF        whiteC          = RGB(0xff,0xff,0xff),
                    fillColor         = GetSysColor(COLOR_HIGHLIGHT),
                    fillTColor        = GetSysColor(COLOR_HIGHLIGHTTEXT),
                    tableBkColor    = RGB(backColor.R, backColor.G, backColor.B),
                    tableBkColor2   = RGB(backColor2.R, backColor2.G, backColor2.B),
                    saveTextColor;

    HPEN            saveDCPen,
                    whitePen = CreatePen (PS_SOLID, 1, whiteC); // checked


    HFONT           saveDCFont;


    HBRUSH          themeFill   = CreateSolidBrush (fillColor); // checked
    checkPointer    (themeFill);

    HRGN            saveRgn     = CreateRectRgn    (0,0,1,1);   // checked
    checkPointer    (saveRgn);

    long            hs, // starting column
                    hf, // ending column
                    vs, // starting row
                    vf, // ending row
                    k,  // loop index
                    t,  // aux variable
                    t2,
                    st; // a few more auxs

    bool            chop,
                    chopv;

    //printf            ("%d %d %d \n", GetRValue (fillColor), GetGValue (fillColor), GetBValue (fillColor));


    GetDisplayRange (relRect, hs, hf, vs, vf);

    long            vsShift = ((settings.width&HY_COMPONENT_H_SCROLL)?HY_SCROLLER_WIDTH:0),
                    hsShift = ((settings.width&HY_COMPONENT_V_SCROLL)?HY_SCROLLER_WIDTH:0);


    _HYRect         saveRel;
    RECT            bRect;

    HDC             offScreenPtr = nil;  // checked

    HBITMAP         offBitmap    = nil;

    if ((vf-vs>2)&&(!isPrinting)) {
        offScreenPtr = CreateCompatibleDC (dc);
        if (offScreenPtr) {
            bRect.left           = bRect.top = 0;
            bRect.right          = relRect->right-relRect->left-hsShift;
            bRect.bottom         = relRect->bottom-relRect->top-vsShift;
            offBitmap            = CreateCompatibleBitmap (dc, bRect.right, bRect.bottom);

            if (!offBitmap) {
                DeleteDC (offScreenPtr);
                offScreenPtr = nil;
            } else {
                DeleteObject (SelectObject (offScreenPtr, offBitmap));
                dc = offScreenPtr;

                saveRel = *relRect;
                relRect->bottom -= relRect->top;
                relRect->top = 0;
                relRect->right  -= relRect->left;
                relRect->left = 0;
            }
        }
    }

    st =            GetClipRgn (dc, saveRgn);

    if (st == 0) {
        DeleteObject (saveRgn);
        saveRgn = nil;
    }

    int             saveBackMode  = GetBkMode    (dc);
    UINT            saveTextAlign = GetTextAlign (dc);
    COLORREF        saveBkColor   = GetBkColor   (dc);

    SetTextAlign   (dc, TA_BASELINE);
    SetBkMode      (dc, TRANSPARENT);

    RECT clipRect = {relRect->left,relRect->top,
                     relRect->right-hsShift,
                     relRect->bottom-vsShift
                    },

         anotherRect,
         clipRect2,
         clipRect3;

    POINT
    mapPts[2];

    mapPts[0] = (POINT) {
        clipRect.left,clipRect.top
    };
    mapPts[1] = (POINT) {
        clipRect.right+1,clipRect.bottom+1
    };

    LPtoDP (dc, mapPts,2);

    HRGN          tempRgn = CreateRectRgn (mapPts[0].x, mapPts[0].y, mapPts[1].x, mapPts[1].y);
    checkPointer  (tempRgn);

    SelectClipRgn (dc, tempRgn);

    saveTextColor = GetTextColor (dc);

    ::SetTextColor (dc, RGB(textColor.R, textColor.G, textColor.B));

    anotherRect = clipRect;

    t = relRect->top-vOrigin;

    saveDCPen = (HPEN)SelectObject (dc, whitePen);

    for (k=vs; k<=vf; k++) {
        anotherRect.top = k?verticalSpaces.lData[k-1]+t:relRect->top;
        anotherRect.bottom = verticalSpaces.lData[k]+t;
        if (cellTypes.lData[k*horizontalSpaces.lLength]&HY_TABLE_BEVELED) {
            FillRect (dc,&anotherRect,backPattern2);
            if ((k==vs)||(k<vf)) {
                SelectObject (dc, menuLine2);
                MoveToEx (dc,anotherRect.left,anotherRect.bottom-1,nil);
                LineTo (dc,anotherRect.right,anotherRect.bottom-1);
                SelectObject (dc, menuLine1);
                MoveToEx (dc,anotherRect.left,anotherRect.bottom-2,nil);
                LineTo (dc,anotherRect.right,anotherRect.bottom-2);
                SelectObject (dc, whitePen);
            }
        } else {
            FillRect (dc,&anotherRect,backPattern);
            if ((k==vs)||(k<vf)) {
                MoveToEx (dc,anotherRect.left,anotherRect.bottom-1,nil);
                LineTo (dc,anotherRect.right,anotherRect.bottom-1);
            }
        }
    }

    saveTextColor = GetTextColor (dc);
    ::SetTextColor (dc, RGB (textColor.R, textColor.G, textColor.B));

    saveDCFont    = (HFONT)SelectObject (dc, tableFont);
    st = 0;

    if (hf<horizontalSpaces.lLength-1) {
        st = hf;
    } else {
        st = hf-1;
    }
    t = relRect->left-hOrigin-2;

    if ((selectionType & HY_TABLE_NO_COLS_LINES) == 0) {
        SelectObject (dc,menuLine1);
        for (k=hs; k<=st; k++) {
            t2 = t+horizontalSpaces.lData[k];
            MoveToEx (dc,t2,relRect->top,nil);
            LineTo (dc,t2,relRect->bottom);
        }
        SelectObject (dc,menuLine2);
        t++;
        for (k=hs; k<=st; k++) {
            t2 = t+horizontalSpaces.lData[k];
            MoveToEx (dc,t2,relRect->top,nil);
            LineTo (dc,t2,relRect->bottom);
        }
    }

    bool  highlightColorOn = false;

    for (k=vs; k<=vf; k++) {
        anotherRect.top = relRect->top-vOrigin+1;
        anotherRect.bottom = anotherRect.top+verticalSpaces.lData[k]-1;
        if (k) {
            anotherRect.top+=verticalSpaces.lData[k-1];
        }

        long    t3,
                st2,
                w = anotherRect.bottom-anotherRect.top,
                w2,
                shift = (w-textFont.size)/2-1;

        if (anotherRect.bottom>relRect->bottom-vsShift) {
            anotherRect.bottom = relRect->bottom-vsShift;
            chopv = false;
        } else {
            chopv = true;
        }

        for (t2=hs; t2<=hf; t2++) {
            t3 = k*horizontalSpaces.lLength+t2;
            if (t3 == editCellID) {
                continue;
            }
            anotherRect.left  = relRect->left-hOrigin+1;
            anotherRect.right = anotherRect.left+horizontalSpaces.lData[t2]-1;

            if (t2) {
                anotherRect.left+=horizontalSpaces.lData[t2-1];
            }

            clipRect2 = anotherRect;
            w2        = anotherRect.right-anotherRect.left;

            chop = true;

            if (anotherRect.right>relRect->right-hsShift) {
                anotherRect.right=relRect->right-hsShift;
                chop = false;
            } else {
                chop = true;
            }

            if ((t2==hs)||(k==vs)) {
                IntersectRect (&clipRect3,&anotherRect,&clipRect);
                DeleteObject (tempRgn);

                mapPts[0] = (POINT) {
                    clipRect3.left,clipRect3.top
                };
                mapPts[1] = (POINT) {
                    clipRect3.right+1,clipRect3.bottom+1
                };

                LPtoDP (dc, mapPts,2);

                tempRgn = CreateRectRgn (mapPts[0].x, mapPts[0].y, mapPts[1].x, mapPts[1].y);
                checkPointer  (tempRgn);
                SelectClipRgn (dc, tempRgn);
            } else {
                DeleteObject (tempRgn);
                mapPts[0] = (POINT) {
                    anotherRect.left,anotherRect.top
                };
                mapPts[1] = (POINT) {
                    anotherRect.right+1,anotherRect.bottom+1
                };

                LPtoDP (dc, mapPts,2);

                tempRgn = CreateRectRgn (mapPts[0].x, mapPts[0].y, mapPts[1].x, mapPts[1].y);
                checkPointer  (tempRgn);
                SelectClipRgn (dc, tempRgn);
            }

            if (cellTypes.lData[t3]&HY_TABLE_SELECTED) {
                if (chopv) {
                    anotherRect.bottom--;
                }
                if (chop) {
                    anotherRect.right-=2;
                    FillRect (dc,&anotherRect,themeFill);
                    anotherRect.right+=2;
                } else {
                    FillRect (dc,&anotherRect,themeFill);
                }

                if (chopv) {
                    anotherRect.bottom++;
                }
            }


            if (cellTypes.lData[t3]&HY_TABLE_ICON) {
                if (!isPrinting) {
                    _SimpleList     * cellList = (_SimpleList*)cellData.lData[t3];

                    if (w2-4>cellList->lData[1]) {
                        t3 = (w2-cellList->lData[1])/2;
                        clipRect2.left+=t3;
                    }
                    clipRect2.right = clipRect2.left+cellList->lData[1];
                    if (w-2>cellList->lData[2]) {
                        t3 = (w-cellList->lData[2])/2;
                        clipRect2.top+=t3;
                    }
                    clipRect2.bottom=clipRect2.top+cellList->lData[2];

                    if (cellList->lLength==3) {

                        HBITMAP aPic = (HBITMAP)cellList->lData[0];
                        if (aPic) {
                            DrawTransparentBitmap (dc, aPic, clipRect2.left, clipRect2.top, clipRect2.right-clipRect2.left+1, clipRect2.bottom-clipRect2.top+1,RGB(255,255,255));
                        }
                    } else {
                        if ((cellList->lData[3]==HY_TABLE_COLOR_BOX)||(cellList->lData[3]==HY_TABLE_COLOR_CIRCLE)) {
                            _HYColor    c   = LongToHYColor    (cellList->lData[0]);
                            HBRUSH      clr;


                            if (cellList->lData[3]==HY_TABLE_COLOR_BOX) {
                                clr =  CreateSolidBrush (RGB(c.R,c.G,c.B));
                                checkPointer (clr);
                                FillRect  (dc, &clipRect2,clr);
                                DeleteObject (clr);
                            } else {
                                COLORREF        trColor = RGB(0,0,0);
                                HPEN            trPen   = CreatePen (PS_NULL, 1, trColor),
                                                saveTPen;

                                checkPointer    (trPen);

                                saveTPen     = (HPEN)SelectObject (dc, trPen);
                                RECT            circRect = clipRect2;
                                InflateRect    (&circRect,1,1);
                                Ellipse        (dc,circRect.left, circRect.top, circRect.right, circRect.bottom);


                                clr =  CreateSolidBrush (RGB(c.R/2,c.G/2,c.B/2));
                                checkPointer (clr);

                                HBRUSH          saveBRUSH = (HBRUSH)SelectObject (dc, clr);

                                InflateRect    (&circRect,-1,-1);
                                trColor      = RGB(c.R/2,c.G/2,c.B/2);
                                //trPen      = CreatePen (PS_SOLID, 1, trColor);
                                //checkPointer   (trPen);
                                //DeleteObject   (SelectObject   (dc, trPen));
                                Ellipse        (dc,circRect.left, circRect.top, circRect.right, circRect.bottom);

                                InflateRect    (&circRect,-1,-1);
                                trColor        = RGB(c.R/1.25,c.G/1.25,c.B/1.25);

                                clr =  CreateSolidBrush (RGB(c.R/1.25,c.G/1.25,c.B/1.25));
                                checkPointer (clr);
                                DeleteObject (SelectObject (dc, clr));

                                //trPen = CreatePen (PS_SOLID, 1, trColor);
                                //checkPointer    (trPen);
                                //DeleteObject  (SelectObject (dc, trPen));

                                Ellipse        (dc,circRect.left, circRect.top, circRect.right, circRect.bottom);

                                clr =  CreateSolidBrush (RGB(c.R,c.G,c.B));
                                checkPointer (clr);
                                DeleteObject (SelectObject (dc, clr));

                                InflateRect    (&circRect,-1,-1);
                                //trColor      = RGB(c.R,c.G,c.B);
                                //trPen = CreatePen (PS_SOLID, 1,trColor);
                                //checkPointer    (trPen);

                                DeleteObject   (SelectObject (dc, trPen));
                                Ellipse        (dc,circRect.left, circRect.top, circRect.right, circRect.bottom);
                                DeleteObject   (SelectObject   (dc, saveTPen));
                                DeleteObject   (SelectObject   (dc, saveBRUSH));
                            }

                        }
                    }
                }
            } else { // text
                st2 = cellTypes.lData[t3]&HY_TABLE_STYLEMASK;

                if (st!=st2) {
                    HFONT      setFont = tableFont;
                    st = st2;
                    if (st&HY_TABLE_BOLD) {
                        if (st&HY_TABLE_ITALIC) {
                            setFont = tableFontBI;
                        } else {
                            setFont = tableFontB;
                        }
                    } else if (st&HY_TABLE_ITALIC) {
                        setFont = tableFontI;
                    }

                    SelectObject (dc, setFont);
                }


                _String  *thisCell = (_String*)cellData.lData[t3];

                if (cellTypes.lData[t3]&HY_TABLE_SELECTED) {
                    if (!highlightColorOn) {
                        ::SetTextColor (dc, fillTColor);
                        highlightColorOn = true;
                    }
                } else {
                    if (highlightColorOn) {
                        ::SetTextColor (dc, RGB(textColor.R, textColor.G, textColor.B));
                        highlightColorOn = false;
                    }
                }
                TextOut (dc,anotherRect.left+textFont.size/3, anotherRect.top+shift+textFont.size, thisCell->sData,thisCell->sLength);

                if (cellTypes.lData[t3]&HY_TABLE_PULLDOWN) {
                    if (!isPrinting) {
                        clipRect2.right-=4;
                        clipRect2.left=clipRect2.right-tPDMw;
                        if (w-2>tPDMh) {
                            t3 = (w-tPDMh)/2;
                            clipRect2.top+=t3;
                        }
                        clipRect2.bottom=clipRect2.top+tPDMh;

                        BITMAP theBM;
                        GetObject (tablePDMenuIcon, sizeof (BITMAP), &theBM);

                        if (otherDC) {
                            SelectObject         (otherDC, tablePDMenuIcon);

                            StretchBlt           (dc, clipRect2.left, clipRect2.top, clipRect2.right-clipRect2.left+1, clipRect2.bottom-clipRect2.top+1,
                                                  otherDC, 0, 0, theBM.bmWidth, theBM.bmHeight, SRCCOPY);
                            SelectObject         (otherDC,oDCBM);
                        } else {
                            _String errMsg = _String ("Failed to make CompatibleDC in _HYTable::_Paint");
                            ReportWarning (errMsg);
                        }
                    }
                }
            }
        }
    }


    SetBkMode      (dc, saveBackMode);
    SetTextAlign   (dc, saveTextAlign);
    SetBkColor     (dc, saveBkColor);

    SelectClipRgn  (dc, saveRgn);

    if (saveRgn) {
        DeleteObject (saveRgn);
    }

    ::SetTextColor (dc, saveTextColor);

    SelectObject   (dc, saveDCPen);
    SelectObject   (dc, saveDCFont);
    DeleteObject   (whitePen);

    if (themeFill) {
        DeleteObject (themeFill);
    }

    if (tempRgn) {
        DeleteObject (tempRgn);
    }

    if (offScreenPtr) {
        *relRect    = saveRel;
        OffsetRect   (&clipRect,relRect->left, relRect->top);

        dc = (HDC)relRect->width;

        BitBlt   (dc,clipRect.left,clipRect.top, clipRect.right-clipRect.left+1,clipRect.bottom-clipRect.top+1, offScreenPtr, 0,0, SRCCOPY);
        DeleteDC     (offScreenPtr);
        DeleteObject (offBitmap);
    }

    if (editBox) {
        UpdateWindow (editBox);
    }

    _HYPlatformComponent::_Paint(p);
}


//__________________________________________________________________

bool        _HYTable::_ProcessOSEvent (Ptr vEvent)
{
    _HYWindowsUIMessage*    theEvent = (_HYWindowsUIMessage*)vEvent;

    static  int     lastH = 0,
                    lastV = 0;

    long            k,
                    h,
                    v,
                    vsShift = ((settings.width&HY_COMPONENT_H_SCROLL)?HY_SCROLLER_WIDTH:0),
                    hsShift = ((settings.width&HY_COMPONENT_V_SCROLL)?HY_SCROLLER_WIDTH:0);

    switch (theEvent->iMsg) {
    case WM_LBUTTONUP: {
        lastH = (short)LOWORD (theEvent->lParam);
        lastV = (short)HIWORD (theEvent->lParam);

        bool    isInComponent = ((lastH>=rel.left)&&(lastV>=rel.top)&&(lastH<rel.right)&&(lastV<rel.bottom));

        if (cursorState == HY_TABLE_SIZE_CURSOR) {
            _ResetCursorState ();
            //ReleaseCapture      ();
        } else if (cursorState == HY_TABLE_DRAG_CURSOR) {
            if (activeColumn >= 0) {
                if (messageRecipient) {
                    ReleaseCapture ();
                    ((_HYTWindow*)messageRecipient)->trackMouseComponent = nil;
                }
                if (activeColumn2>=0) {
                    _HiliteRowForDrag (activeColumn2,activeColumn);
                }

                if (activeColumn!=activeColumn2) {
                    EditBoxHandler (-1,rel);
                    if (isInComponent) {
                        DragRow (activeColumn,activeColumn2);
                    }
                }

                activeColumn = -1;
            }
            _ResetCursorState ();
        } else if (cursorState == HY_TABLE_EDIT_CURSOR) {
            if (messageRecipient) {
                ReleaseCapture ();
                ((_HYTWindow*)messageRecipient)->trackMouseComponent = nil;
            }

            RECT      clippingRect = HYRect2Rect (rel),
                      paintRect;

            IntersectRect (&paintRect,&limits,&clippingRect);
            _FrameRect (paintRect);

            _HYRect     outlineHRect;
            outlineHRect.left   = limits.left;
            outlineHRect.right  = limits.right+hsShift;
            outlineHRect.top    = limits.top;
            outlineHRect.bottom = limits.bottom+vsShift;
            long    hs,hf,vs,vf;
            hOrigin += limits.left-rel.left;
            vOrigin += limits.top-rel.top;
            GetDisplayRange (&outlineHRect,hs,hf,vs,vf);
            hOrigin -= limits.left-rel.left;
            vOrigin -= limits.top-rel.top;
            ExpungeSelection();
            if ((hf>=hs)||(vs>=vf)) {
                _SimpleList sel;
                if (selectionType&HY_TABLE_SEL_ROWS) {
                    sel.RequestSpace (vf-vs+1);
                    for (h=vs; h<=vf; h++) {
                        sel<<h;
                    }
                    SetRowSelection (sel);
                } else if (selectionType&HY_TABLE_SEL_COLS) {
                    sel.RequestSpace (hf-hs+1);
                    for (v=hs; v<=hf; h++) {
                        sel<<v;
                    }
                    SetColumnSelection(sel);

                } else {
                    sel.RequestSpace ((vf-vs+1)*(hf-hs+1));
                    for (h=hs; h<=hf; h++)
                        for (v=vs; v<=vf; v++) {
                            sel << v*horizontalSpaces.lLength + h;
                        }
                    SetSelection (sel,true);
                    _MarkCellsForUpdate (sel);
                }
            }

            if (messageRecipient) {
                ReleaseCapture ();
                ((_HYTWindow*)messageRecipient)->trackMouseComponent = nil;
            }

            _ResetCursorState ();
        }

        break;
    }

    case WM_LBUTTONDOWN:
    case WM_LBUTTONDBLCLK: {
        lastH = (short)LOWORD (theEvent->lParam);
        lastV = (short)HIWORD (theEvent->lParam);

        POINT   downWhere = (POINT) {
            lastH, lastV
        };

        ClientToScreen (parentWindow, &downWhere);

        if ((selectionType&HY_TABLE_FOCUSABLE)&&messageRecipient&&((selectionType&HY_TABLE_IS_FOCUSED)==0)) {
            messageRecipient->ProcessEvent(generateKeyboardFocusEvent (GetID()));
        }

        if ((cursorState == HY_TABLE_SIZE_CURSOR)&&(theEvent->iMsg!= WM_LBUTTONDBLCLK)) {
            // do drag here
            EditBoxHandler (-1,rel);
            if (!DragDetect (parentWindow, downWhere)) {
                _ResetCursorState ();
                return    true;
            }

            //SetCapture      (parentWindow);

            limits.left   = lastH;
            limits.top    = rel.top;
            limits.bottom = rel.bottom-vsShift;

            limits.right = rel.right-hsShift;

            lastH += hOrigin;

            for (activeColumn = 0; activeColumn<horizontalSpaces.lLength-1; activeColumn++)
                if (horizontalSpaces.lData[activeColumn]>lastH-2-rel.left) {
                    break;
                }

            if (activeColumn) {
                limits.left = rel.left+horizontalSpaces.lData[activeColumn-1]+3-hOrigin;
            } else {
                limits.left = rel.left+3;
            }

            long dragRes = horizontalSpaces.lData[horizontalSpaces.lLength-1]-
                           rel.right+rel.left-hOrigin+hsShift;

            if (dragRes<lastH-limits.left) {
                limits.left = lastH-dragRes;
            }

            lastH-=hOrigin;

            return true;

        } else {
            if (((k=FindClickedTableCell (lastH-rel.left,lastV-rel.top,h,v))>-1)&&
                    (lastV<rel.bottom-vsShift)&&(lastH<rel.right-hsShift)) {
                if (theEvent->iMsg== WM_LBUTTONDBLCLK) {
                    if (cellTypes.lData[k]&HY_TABLE_EDIT_TEXT) {
                        EditBoxHandler (k,rel);
                    } else if (messageRecipient) {
                        messageRecipient->ProcessEvent (generateTableDblClickEvent(GetID()));
                    }
                    break;
                }
                ModifySelection (h,v,GetAsyncKeyState (VK_SHIFT) & 0x8000, GetAsyncKeyState (VK_CONTROL) & 0x8000, true);
            }

            if (k==-2)
                // process pull-down
            {
                if (messageRecipient) {
                    POINT loc = {lastH, lastV};
                    ClientToScreen (parentWindow, &loc);
                    messageRecipient->ProcessEvent (generateTablePullDownEvent(GetID(),v*horizontalSpaces.lLength+h,
                                                    (((long)loc.x)<<16)+loc.y));
                }
                break;
            }


            if ((cursorState == HY_TABLE_DRAG_CURSOR)&&(theEvent->iMsg!= WM_LBUTTONDBLCLK)) {
                SetCursor (dropOffCursor);
                FindClickedTableCell(lastH-rel.left,lastV-rel.top,k,activeColumn);
                activeColumn2 = -1;
                if (messageRecipient) {
                    SetCapture (parentWindow);
                    ((_HYTWindow*)messageRecipient)->trackMouseComponent = this;
                }

                return true;
            }

            /*  handle row drag here */

            if (((selectionType&HY_TABLE_SINGLE_SELECTION)==0)&&
                    ((selectionType&HY_TABLE_NODRAG_SELECTION)==0)&&
                    (lastH<rel.right-hsShift)&&
                    (lastV<rel.bottom-vsShift)&&
                    DragDetect (parentWindow, downWhere)) {
                cursorState = HY_TABLE_EDIT_CURSOR;
                limits.left = limits.right  = lastH;
                limits.top  = limits.bottom = lastV;

                textBoxRect.left   = -1;
                textBoxRect.right  = lastH;
                textBoxRect.bottom = lastV;

                if (messageRecipient) {
                    SetCapture (parentWindow);
                    ((_HYTWindow*)messageRecipient)->trackMouseComponent = this;
                }

                return true;
            }

        }

    }
    break;

    case WM_KEYDOWN: {

        bool    ctlDown = (GetAsyncKeyState (VK_CONTROL) & 0x8000);

        switch (theEvent->wParam) {
        case VK_UP: // up
            HandleKeyMove (0,ctlDown);
            return true;
        case VK_DOWN: // down
            HandleKeyMove (1,ctlDown);
            return true;
        case VK_LEFT: // left
            HandleKeyMove (2,ctlDown);
            return true;
        case VK_RIGHT: // right
            HandleKeyMove (3,ctlDown);
            return true;
        }

        break;
    }

    case WM_MOUSEMOVE: {
        POINT   downWhere = (POINT) {
            (short)LOWORD (theEvent->lParam), (short)HIWORD (theEvent->lParam)
        };

        if (theEvent->wParam & MK_LBUTTON) { // left button down
            if (cursorState == HY_TABLE_SIZE_CURSOR) {
                POINT   currentPoint = (POINT) {
                    LOWORD (theEvent->lParam),HIWORD (theEvent->lParam)
                };
                if (PtInRect (&limits, currentPoint)) {
                    if (currentPoint.x-lastH) {
                        SetColumnSpacing (activeColumn,currentPoint.x-lastH,true);
                        if (messageRecipient)
                            messageRecipient->ProcessEvent (generateTableResizeCEvent(GetID(),
                                                            activeColumn,currentPoint.x-lastH));
                        lastH = currentPoint.x;
                    }

                }
            } else if ((cursorState == HY_TABLE_DRAG_CURSOR)&&(activeColumn>=0)) {
                long h,v,k;

                if ((downWhere.y>rel.bottom-vsShift)||(downWhere.y<rel.top)) {
                    if (activeColumn2>=0) {
                        _HiliteRowForDrag (activeColumn2,activeColumn);
                        activeColumn2 = -1;
                    }
                    h = verticalSpaces.lData[verticalSpaces.lLength-1]/verticalSpaces.lLength;
                    _ScrollVPixels ((downWhere.y<rel.top)?-h:h);
                    break;
                }

                if (downWhere.x>rel.right-hsShift) {
                    downWhere.x=rel.right-hsShift;
                }

                if ( (lastH!=downWhere.x)|| (lastV!=downWhere.y)) {
                    k = FindClickedTableCell(downWhere.x-rel.left,downWhere.y-rel.top,h,v);
                    if ((v!=activeColumn2)&&(k>-1)) {
                        if (activeColumn2>=0) {
                            _HiliteRowForDrag (activeColumn2,activeColumn);
                        }
                        if ((v!=activeColumn)&&(!(cellTypes.lData[k]&HY_TABLE_CANTSELECT))) {
                            _HiliteRowForDrag (v,activeColumn);
                            activeColumn2 = v;
                        } else {
                            activeColumn2 = -1;
                        }
                    }
                    lastH = downWhere.x;
                    lastV = downWhere.y;
                }
            } else if (cursorState == HY_TABLE_EDIT_CURSOR) {
                RECT      clippingRect = HYRect2Rect (rel),
                          paintRect;

                while (1) {
                    if (downWhere.y>rel.bottom-vsShift) {
                        if (rel.bottom-rel.top-vsShift+vOrigin <  verticalSpaces.lData[verticalSpaces.lLength-1]-1) {
                            long h = verticalSpaces.lData[verticalSpaces.lLength-1]/verticalSpaces.lLength;

                            IntersectRect (&paintRect,&limits,&clippingRect);
                            _FrameRect (paintRect);

                            _ScrollVPixels (h);

                            limits.top -= h;
                            IntersectRect (&paintRect,&limits,&clippingRect);
                            _FrameRect    (paintRect);

                            lastV -= h;
                            textBoxRect.bottom -= h;
                            break;
                        }
                        downWhere.y=rel.bottom-vsShift;
                    }
                    if (downWhere.x>rel.right-hsShift) {
                        if (rel.right-rel.left-hsShift+hOrigin <  horizontalSpaces.lData[horizontalSpaces.lLength-1]-1) {
                            long h = horizontalSpaces.lData[horizontalSpaces.lLength-1]/horizontalSpaces.lLength;

                            IntersectRect (&paintRect,&limits,&clippingRect);
                            _FrameRect    (paintRect);

                            _ScrollHPixels (h);

                            limits.left -= h;

                            IntersectRect (&paintRect,&limits,&clippingRect);
                            _FrameRect    (paintRect);

                            lastH -= h;
                            textBoxRect.right -= h;
                            break;
                        }
                        downWhere.x=rel.right-hsShift;
                    }
                    if (downWhere.y<rel.top) {
                        if (vOrigin>0) {
                            long h = verticalSpaces.lData[verticalSpaces.lLength-1]/verticalSpaces.lLength;

                            IntersectRect (&paintRect,&limits,&clippingRect);
                            _FrameRect    (paintRect);

                            _ScrollVPixels (-h);

                            limits.top += h;
                            IntersectRect (&paintRect,&limits,&clippingRect);
                            _FrameRect    (paintRect);

                            lastV += h;
                            textBoxRect.bottom += h;
                            break;
                        }
                        downWhere.y=rel.top;
                    }
                    if (downWhere.x<rel.left) {
                        if (hOrigin>0) {
                            long h = horizontalSpaces.lData[horizontalSpaces.lLength-1]/horizontalSpaces.lLength;

                            IntersectRect (&paintRect,&limits,&clippingRect);
                            _FrameRect    (paintRect);

                            _ScrollHPixels (-h);

                            limits.left += h;

                            IntersectRect (&paintRect,&limits,&clippingRect);
                            _FrameRect    (paintRect);

                            lastH += h;
                            textBoxRect.right += h;
                            break;
                        }
                        downWhere.x=rel.left;
                    }

                    if ( (lastH!=downWhere.x)|| (lastV!=downWhere.y)) {
                        if (textBoxRect.left>=0) {
                            IntersectRect (&paintRect,&limits,&clippingRect);
                            _FrameRect    (paintRect);
                        }

                        if (downWhere.x > textBoxRect.right) {
                            limits.right = downWhere.x;
                        } else {
                            limits.right =  textBoxRect.right;
                            limits.left = downWhere.x;
                        }
                        if (downWhere.y >  textBoxRect.bottom) {
                            limits.bottom = downWhere.y;
                        } else {
                            limits.bottom = textBoxRect.bottom;
                            limits.top = downWhere.y;
                        }

                        IntersectRect (&paintRect,&limits,&clippingRect);
                        _FrameRect    (paintRect);

                        lastH = downWhere.x;
                        lastV = downWhere.y;

                        textBoxRect.left = 1;
                    }
                    break;
                }
                return true;
            }

        } else { // treat as plain move
            bool  ch = (!(selectionType&HY_TABLE_DONT_SIZE))&&(CheckForHSizeLocation(downWhere.x-rel.left))&&(downWhere.y<rel.bottom-vsShift);

            if (ch&&(cursorState!=HY_TABLE_SIZE_CURSOR)) {
                cursorState = HY_TABLE_SIZE_CURSOR;
                SetCursor(hSizeCursor);
            } else if ((!ch)&&(cursorState==HY_TABLE_SIZE_CURSOR)) {
                SetCursor(LoadCursor (nil, IDC_ARROW));
                cursorState = 0;
            }

            if (selectionType & HY_TABLE_SEL_ROWS) {
                if (!ch) {
                    k = FindClickedTableCell(downWhere.x-rel.left,downWhere.y-rel.top,h,v);
                    if (k>=0) {
                        if ((cursorState != HY_TABLE_DRAG_CURSOR)&&
                                ((selectionType&HY_TABLE_NODRAG_SELECTION)==0)) {
                            if (cellTypes.lData[k]&HY_TABLE_SELECTED) {
                                if (IsRowSelectionSimple()) {
                                    SetCursor (pickUpCursor);
                                    cursorState = HY_TABLE_DRAG_CURSOR;
                                    activeColumn = -1;
                                }
                            }
                        } else {
                            if (((selectionType&HY_TABLE_NODRAG_SELECTION)==0)&&(!(cellTypes.lData[k]&HY_TABLE_SELECTED))) {
                                SetCursor(LoadCursor (nil, IDC_ARROW));
                                cursorState = 0;
                            }
                        }
                    }
                }
            }
        }
        return true;
    }
    }
    return _HYPlatformComponent::_ProcessOSEvent (vEvent);
}




//__________________________________________________________________

void        _HYTable::_PrintTable (_SimpleList& columns, _SimpleList& rows, _HYTable* ch)
{

    if ((columns.lLength == 0)||(rows.lLength == 0)) {
        return;
    }

    DOCINFO                 di = {sizeof(DOCINFO), "HYPHY.out", NULL };
    PRINTDLG                pd;
    BOOL                    SuccessFlag;

    pd.lStructSize         = sizeof(PRINTDLG);
    pd.hwndOwner           = parentWindow;
    pd.hDevMode            = NULL;
    pd.hDevNames           = NULL;
    pd.hDC                 = NULL;
    pd.Flags               = PD_COLLATE | PD_RETURNDC | PD_NOSELECTION;
    pd.nFromPage           = 1;
    pd.nToPage             = 0xffff;
    pd.nMinPage            = 1;
    pd.nMaxPage            = 0xffff;
    pd.nCopies             = 1;
    pd.hInstance           = NULL;
    pd.lCustData           = 0L;
    pd.lpfnPrintHook       = NULL;
    pd.lpfnSetupHook       = NULL;
    pd.lpPrintTemplateName = NULL;
    pd.lpSetupTemplateName = NULL;
    pd.hPrintTemplate      = NULL;
    pd.hSetupTemplate      = NULL;

    if (!PrintDlg(&pd)) {
        return;
    }

    if (pd.hDC == NULL) {
        pd.hDC = GetPrinterDeviceContext(parentWindow);
    }


    EnableWindow(parentWindow, FALSE);

    SuccessFlag   = TRUE;
    UserAbortFlag = FALSE;

    PrintDialogHandle = CreateDialog(GetModuleHandle(NULL), (LPCTSTR)"PrintDlgBox", parentWindow,
                                     PrintDialogProc);
    SetDlgItemText(PrintDialogHandle, IDD_FNAME, "Table Printing...");

    SetAbortProc(pd.hDC, AbortProc);

    if (StartDoc(pd.hDC, &di) > 0) {
        HDC         windowDC = GetDC (parentWindow);

        long        printW = GetDeviceCaps(pd.hDC, HORZRES),
                    printH = GetDeviceCaps(pd.hDC, VERTRES),

                    hRes = GetDeviceCaps(pd.hDC, LOGPIXELSX),
                    vRes = GetDeviceCaps(pd.hDC, LOGPIXELSY),

                    screenHRes = GetDeviceCaps(windowDC, LOGPIXELSX),
                    screenVRes = GetDeviceCaps(windowDC, LOGPIXELSY),

                    tW = 0,
                    tH = 0,
                    cC = 0,
                    cE = 0,
                    cP = 1,
                    t,
                    i,
                    pH = ch?ch->GetRowSpacing (0):0,
                    fromPage = pd.nMinPage,
                    toPage   = pd.nMaxPage;


        if (pd.Flags & PD_PAGENUMS) {
            fromPage = pd.nFromPage;
            toPage   = pd.nToPage;
        }
        ReleaseDC   (parentWindow, windowDC);
        hRes = printW*((_Parameter)screenHRes/hRes);
        vRes = printH*((_Parameter)screenVRes/vRes);
        screenHRes = printW;
        screenVRes = printH;

        printW = hRes;
        printH = vRes;
        if (ch) {
            printH -= pH;
            if (printH <= 0) {
                _String errMsg ("Table header is too tall to fit on the page.");
                WarnError (errMsg);
                terminateExecution = false;
                EndDoc(pd.hDC);
                DeleteDC (pd.hDC);
                return;
            }
        }

        _HYRect  relDim = {0,0,0,0,HY_COMPONENT_NO_SCROLL};

        if (ch)
            for (i=0; i<columns.lLength; i++) {
                tW += ch->GetColumnSpacing (columns.lData[i]);
            }
        else
            for (i=0; i<columns.lLength; i++) {
                tW += GetColumnSpacing (columns.lData[i]);
            }

        relDim.left    = relDim.right   = tW > printW ? printW : tW;

        while ((cP < fromPage)&&SuccessFlag&&(!UserAbortFlag)) {
            t = 0;
            while (cC < rows.lLength) {
                i = GetRowSpacing (rows.lData[cC]);
                if ( t+i > printH) {
                    cP ++;
                    if (i > printH) {
                        cC++;
                    }
                    break;
                } else {
                    t += i;
                    cC ++;
                }
            }
        }

        cE = cC;

        for (long pageCount = fromPage; pageCount<=toPage && (cC < rows.lLength); pageCount++) {
            t = 0;
            while (cE < rows.lLength) {
                i = GetRowSpacing (rows.lData[cE]);
                if ( t+i > printH) {
                    if (i > printH) {
                        t = printH;
                        cE++;
                    }
                    break;
                } else {
                    t += i;
                    cE ++;
                }
            }

            relDim.top = relDim.bottom = t+pH;

            _HYTable        *thisPage = new _HYTable (relDim,(Ptr)parentWindow,
                    cE-cC+(ch?1:0),columns.lLength,20,20,HY_TABLE_STATIC_TEXT);

            checkPointer    (thisPage);


            if (StartPage (pd.hDC) <= 0) {
                SuccessFlag = FALSE;
                break;
            }

            SetMapMode  (pd.hDC, MM_ISOTROPIC);
            SetWindowExtEx (pd.hDC, hRes, vRes,nil);
            SetViewportExtEx (pd.hDC, screenHRes, screenVRes, nil);

            thisPage->SetFont       (textFont);
            thisPage->SetBackColor  (backColor);
            thisPage->SetBackColor2 (backColor2);
            thisPage->SetTextColor  (textColor);

            if (ch)
                for (i=0; i<columns.lLength; i++) {
                    thisPage->SetColumnSpacing (i,ch->GetColumnSpacing (columns.lData[i])-20, false);
                }
            else
                for (i=0; i<columns.lLength; i++) {
                    thisPage->SetColumnSpacing (i,GetColumnSpacing (columns.lData[i])-20, false);
                }

            t = 0;
            if (ch) {
                thisPage->SetRowSpacing (0,pH-20,false);
                for (i=0; i<columns.lLength; i++) {
                    BaseRef cellData = ch->GetCellData(columns.lData[i],0);
                    cellData->nInstances++;
                    thisPage->SetCellData(cellData,0,i,ch->cellTypes.lData[columns.lData[i]]&HY_TABLE_DESELECT,false);
                }
                t = 1;
            }

            for (cP = cC; cP < cE; cP++,t++) {
                long     rI = rows.lData[cP];
                thisPage->SetRowSpacing (t,GetRowSpacing(rI)-20,false);
                for (i=0; i<columns.lLength; i++) {
                    BaseRef cellData = GetCellData(columns.lData[i],rI);
                    cellData->nInstances++;
                    thisPage->SetCellData(cellData,t,i,cellTypes.lData[rI*horizontalSpaces.lLength+columns.lData[i]]&HY_TABLE_DESELECT,false);
                }
            }

            _HYRect relDim2 = relDim;
            relDim2.left = relDim2.top = 1;
            relDim2.right = -relDim2.right-1;
            relDim2.bottom ++;

            relDim2.width = (long)pd.hDC;

            thisPage->_Paint ((Ptr)&relDim2); // may need to disable double buffering

            RECT wrect = HYRect2Rect (relDim2);
            
            FrameRect (pd.hDC, &wrect, _BLACKBRUSH_);

            if (EndPage (pd.hDC) <= 0) {
                SuccessFlag = FALSE;
            }

            cC = cE;
        }
    } else {
        SuccessFlag = FALSE;
    }

    if (SuccessFlag) {
        SuccessFlag = (EndDoc(pd.hDC)>0);
    }

    if (!UserAbortFlag) {
        EnableWindow(parentWindow, TRUE);
        DestroyWindow(PrintDialogHandle);
    }

    DeleteDC (pd.hDC);

    if (!SuccessFlag && !UserAbortFlag) {
        _String errMsg = _String("Failed to print the table. Windows Error:") & (long)GetLastError();
        ProblemReport (errMsg,nil);
    }
}

//__________________________________________________________________

void        _HYPlatformTable::_ResetCursorState (void)
{
    SetCursor (LoadCursor (nil, IDC_ARROW));
    cursorState = 0;
}

//__________________________________________________________________

void        _HYPlatformTable::_FrameRect        (RECT& theRect)
{
    HPEN      aPen = CreatePen (PS_DOT, 1, RGB(0,0,0));
    HDC       dc   = GetDC (((_HYTable*)this)->parentWindow);

    aPen           = (HPEN)SelectObject (dc, aPen);

    int       saveBkMode = GetBkMode (dc),
              saveROP2   = GetROP2   (dc);

    SetBkMode (dc, TRANSPARENT);
    SetROP2   (dc, R2_NOTXORPEN);


    Rectangle (dc, theRect.left, theRect.top, theRect.right, theRect.bottom);

    SetBkMode (dc, saveBkMode);
    SetROP2   (dc, saveROP2);

    DeleteObject (SelectObject (dc, aPen));

    ReleaseDC (((_HYTable*)this)->parentWindow, dc);
}


//EOF