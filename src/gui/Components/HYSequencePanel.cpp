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

#include "HYComponent.h"
#include "HYEventTypes.h"
#include "HYCanvas.h"
#include "HYSequencePanel.h"
#include "HYUtils.h"
#include "HYPlatformWindow.h"
#include "math.h"
#include "HYDialogs.h"



#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif
//__________________________________________________________________

_HYSequencePane::_HYSequencePane(_HYRect s, Ptr p, int h, int w):_HYCanvas (s,p,h,w,32)
{
    startColumn = endColumn = startRow = endRow = headerWidth = 0;
    for (long i=0; i<256; i++) {
        colorMap[i] = 0;
    }
    _HYFont displayFont;
#ifdef __MAC__
    displayFont.face = "Monaco";
    displayFont.size = 10;
#endif
#ifdef __WINDOZE__
    displayFont.face = "Courier";
    displayFont.size = 14;
#endif
#ifdef __HYPHY_GTK__
    displayFont.face = _HY_MONO_FONT;
    displayFont.size = 12;
#endif
    displayFont.style = HY_FONT_PLAIN;
    StartDraw();
    _HYGraphicPane::SetFont(displayFont);
    EndDraw();
    charWidth    = GetMaxCharWidth (displayFont)+1;
    invertColors = false;
    backColor = (_HYColor) {
        0xDE,0xDE,0xDE
    };

    headerColor = (_HYColor) {
        180,167,150
    };

    highlightColor = (_HYColor) {
        75,75,75
    };

    characterColors<<0;
    headerWidth = 0;
    blockWidth = 10;
    recentClick = -1;
    showDots = false;
    nameDisplayFlags = HY_SEQUENCE_PANE_NAMES_ALL;
    shortHeaderWidth = 0;
    fullHeaderWidth = 0;
    active = true;
    numbers = true;
    dragScrollCounter = 0;
}

//__________________________________________________________________
void _HYSequencePane::InsertColumn(_String* newColumn, long newIndex, bool update)
{
    if (newColumn->sLength) {
        if (columnStrings.lLength) {
            if (((_String*)columnStrings(0))->sLength!=newColumn->sLength) {
                return;
            }
        }
        columnStrings.InsertElement (newColumn,newIndex,false);
    }
    if (update) {
        BuildPane();
    }
}
//__________________________________________________________________
void _HYSequencePane::RemoveColumn(long killIndex, bool update)
{
    if ((killIndex>=0)&&(killIndex<columnStrings.lLength)) {
        long k;
        if ((k=hiliteColumns.BinaryFind(killIndex))>=0) {
            hiliteColumns.Delete (k);
        }
        columnStrings.Delete(killIndex);
        if (update) {
            BuildPane();
        }
    }
}
//__________________________________________________________________
void _HYSequencePane::SetCharColor (unsigned char c, _HYColor newColor, bool update)
{
    long k      = HYColorToLong (newColor),
         f,
         count   =0;

    for (f=0; f<256; f++)
        if (colorMap[f]==colorMap[c]) {
            count++;
        }

    f = characterColors.Find(k);
    if (count>1) {
        if (f<0) {
            f = characterColors.lLength;
            characterColors<<k;
        }
    } else {
        if (f<0) {
            characterColors.lData[colorMap[c]] = k;
            f = colorMap[c];
        } else {
            count = colorMap[c];
            if (k!=characterColors[count]) {
                characterColors.Delete(count);
                for (k=0; k<256; k++) {
                    if (colorMap[k]>count) {
                        colorMap[k]--;
                    }
                }
            }
        }
    }
    colorMap[c] = f;
    if (update) {
        BuildPane();
    }
}

//__________________________________________________________________
void    _HYSequencePane::SetFont(_HYFont newFont)
{
    charWidth = GetMaxCharWidth (newFont)+1;
    StartDraw ();
    _HYGraphicPane::SetFont (newFont);
    EndDraw ();
    SetHeaders (nil,false);
}

//__________________________________________________________________
void    _HYSequencePane::SetBackColor (_HYColor newColor)
{
    backColor = newColor;
    BuildPane();
    _MarkForUpdate();
}

//__________________________________________________________________
void    _HYSequencePane::SetHeaderColor (_HYColor newColor)
{
    headerColor = newColor;
    BuildPane();
    _MarkForUpdate();
}

//__________________________________________________________________
void    _HYSequencePane::SetHighliteColor (_HYColor newColor)
{
    highlightColor = newColor;
    BuildPane();
    _MarkForUpdate();
}

//__________________________________________________________________
void    _HYSequencePane::BuildHeaders (void)
{
    long h,v,slotHeight = GetSlotHeight(),k;
    long visWidth = _HYCanvas::GetMaxW()-headerWidth-5,
         visHeight = _HYCanvas::GetMaxH()-5,
         selectionIndex = 0;

    if (settings.width&HY_COMPONENT_V_SCROLL) {
        visWidth-=HY_SCROLLER_WIDTH;
    }
    if (settings.width&HY_COMPONENT_H_SCROLL) {
        visHeight-=HY_SCROLLER_WIDTH;
    }

    if (headerWidth) {
        _HYRect  r = {0,0,0,0,1};
        _HYColor    blackC = {0,0,0},
                    selectColor = highlightColor;

        SetColor(headerColor);
        r.right =  headerWidth;
        if (numbers) {
            r.bottom = (GetSlotHeight()+1)-1;
            FillRect(r);
            r.top = r.bottom+2;
        } else {
            r.top = 0;
        }
        r.bottom = visHeight;
        FillRect(r);

        if (numbers) {
            v = slotHeight+(GetSlotHeight()+1);
        } else {
            v = slotHeight;
        }

        if (vselection.lLength) {
            SetColor (selectColor);
            for (h=startRow; h<endRow; h++) {
                if (selectionIndex<vselection.lLength) {
                    k = vselection.lData[selectionIndex];
                    while ((k<h)&&(selectionIndex<vselection.lLength)) {
                        selectionIndex++;
                        k = vselection.lData[selectionIndex];
                    }
                    if (h==k) {
                        _HYRect invR;
                        invR.left = 0;
                        invR.right = headerWidth;
                        invR.top = (h-startRow)*slotHeight+(GetSlotHeight()+1)+2;
                        invR.bottom = invR.top+slotHeight;
                        FillRect (invR);
                        selectionIndex++;
                    }
                } else {
                    break;
                }
            }
        }
        SetColor (blackC);
        for (h=startRow; h<endRow; h++)
            if (nameDisplayFlags&HY_SEQUENCE_PANE_NAMES_SHORT) {
                for (h=startRow; h<endRow; h++,v+=slotHeight) {
                    _String *thisString = (_String*)rowHeaders(speciesIndex.lData[h]);
                    if (thisString->sLength<=HY_SEQUENCE_PANE_SHORT_WIDTH) {
                        DisplayText (*thisString,v,HY_SEQUENCE_PANE_CHAR_SPACING/2,true);
                    } else {
                        _String chopped (*thisString);
                        chopped.Trim(0,HY_SEQUENCE_PANE_SHORT_WIDTH-1);
                        chopped.sData[HY_SEQUENCE_PANE_SHORT_WIDTH-1]='-';
                        DisplayText (chopped,v,HY_SEQUENCE_PANE_CHAR_SPACING/2,true);
                    }
                }
            } else
                for (h=startRow; h<endRow; h++,v+=slotHeight) {
                    DisplayText (*(_String*)rowHeaders(speciesIndex.lData[h]),v,HY_SEQUENCE_PANE_CHAR_SPACING/2,true);
                }
    }
}

//long buildCounter = 0;
//__________________________________________________________________
void    _HYSequencePane::BuildPane (bool setport)
{
    if (setport) {
        StartDraw();
        EraseAll();
    }
    if (columnStrings.lLength) {
        long visWidth  = _HYCanvas::GetMaxW()-headerWidth-5,
             visHeight = _HYCanvas::GetMaxH()-5,
             rowCount = speciesIndex.lLength,
             slotHeight = GetSlotHeight(),h,v,lastColor = -1, selectionIndex = 0,
             c;
        if (settings.width&HY_COMPONENT_V_SCROLL) {
            visWidth-=HY_SCROLLER_WIDTH;
        }
        if (settings.width&HY_COMPONENT_H_SCROLL) {
            visHeight-=HY_SCROLLER_WIDTH;
        }
        if (startColumn<0) {
            startColumn = 0;
        }
        if (startRow<0) {
            startRow = 0;
        }
        if (setport) {
            endColumn = startColumn+visWidth/charWidth;
            h = (endColumn-startColumn)/5;
            endColumn = startColumn+(visWidth-h)/charWidth;
            endRow = startRow+visHeight/slotHeight;
        }
        if (endColumn>columnStrings.lLength) {
            endColumn = columnStrings.lLength;
        }
        if (endRow>rowCount) {
            endRow = rowCount;
        }

        _HYRect     r = {0,0,visHeight,visWidth+headerWidth,1},
                    backCharRect;

        _HYColor    blackC      = {0,0,0},
                    charColor   = {0,0,0},
                    whiteC      = {255,255,255},
                    selectColor = highlightColor;

        if (setport) {
            SetColor    (backColor);
        } else {
            selectColor = (_HYColor) {
                160,160,160
            };
            SetColor      (whiteC);
        }

        FillRect    (r);

        SetColor (blackC);
        BuildHeaders();
        if (numbers) {
            r.left = 0;
            r.right = visWidth+headerWidth;
            r.bottom = (GetSlotHeight()+1);
            SetColor (headerColor);
            FillRect(r);
            if (setport) {
                SetColor (blackC);
                r.top = r.bottom;
                DrawLine    (r);
            }
        }

        if (setport) {
            r.left = r.right = headerWidth;
            r.top = 0;
            r.bottom = visHeight;
            DrawLine    (r);
        }
        r.left = r.right = headerWidth;
        r.top = 0;
        r.bottom = visHeight;
        h = HY_SEQUENCE_PANE_CHAR_SPACING/2+headerWidth;
        visWidth+=headerWidth;
        for (c = startColumn; c < endColumn; c++,h+=charWidth) {
            if (h+charWidth>visWidth) {
                endColumn = c;
                break;
            }

            bool isColumnSelected = false;

            if (selectionIndex<selection.lLength) {
                rowCount = selection.lData[selectionIndex];
                while ((rowCount<c)&&(selectionIndex<selection.lLength)) {
                    selectionIndex++;
                    rowCount = selection.lData[selectionIndex];
                }

                if (c==rowCount) {
                    if (!invertColors) {
                        _HYRect invR;
                        invR.top    = (GetSlotHeight()+1)+1;
                        invR.bottom = visHeight-1;
                        invR.right  = h+charWidth;
                        invR.left   = h;

                        if (c&&(c%blockWidth==0)) {
                            invR.right+=2;
                        } else if (c+1==endColumn) {
                            invR.right+=3;
                        }

                        charColor = GetColor();
                        SetColor (selectColor);
                        FillRect (invR);
                        SetColor (charColor);

                    }
                    selectionIndex++;
                    isColumnSelected = true;
                }
            }

            if (c && (c%blockWidth==0) && (setport|| c>startColumn) ) {
                SetColor (blackC);
                if (numbers) {
                    _String number (c);
#ifdef __MAC__
                    rowCount = h-2-GetVisibleStringWidth(number);
#else
                    rowCount = h-2-GetVisibleStringWidth(number, font);
#endif
                    if (rowCount>headerWidth) {
                        DisplayText (number,slotHeight-3,rowCount,true);
                    }
                }
                if (setport) {
                    r.left = r.right = h;
                    DrawLine (r);
                }
                SetColor (charColor);
                h+=2;
            }

            if (numbers) {
                v = slotHeight+(GetSlotHeight()+1);
            } else {
                v = slotHeight;
            }

            _String       *thisString = (_String*)columnStrings(c);

            unsigned char topChar = showDots?thisString->sData[speciesIndex.lData[0]]:0;

            //if (!setport)
            //SetColor (blackC);

            backCharRect = (_HYRect) {
                v-slotHeight+1,h-1,v+1,h+charWidth-1,0
            };
            if ((c+1)%blockWidth == 0) {
                backCharRect.right++;
            }

            for (long r = startRow; r < endRow; r++,v+=slotHeight) {
                unsigned char thisC = thisString->sData[speciesIndex.lData[r]];
                long     myColor    = colorMap[thisC];

                if (r && thisC==topChar) {
                    thisC = '.';
                }

                //if (setport)
                {
                    if (invertColors && !isColumnSelected) {
                        charColor = LongToHYColor (characterColors.lData[myColor]);
                        if ((long)charColor.R+charColor.G+charColor.B > 100) {
                            SetColor (charColor);
                            FillRect    (backCharRect);
                        }
                        backCharRect.top    += slotHeight;
                        backCharRect.bottom += slotHeight;
                        SetColor    (blackC);
                    } else if (colorMap[thisC]!=lastColor) {
                        lastColor = colorMap[thisC];
                        charColor = LongToHYColor(characterColors.lData[lastColor]);
                        SetColor (charColor);
                    }
                }
                DisplayChar (thisC,v,h);
            }
        }

        if ((!setport)&&(c%blockWidth == 0)) {
            SetColor (blackC);
            if (numbers) {
                _String number (c);
#ifdef __MAC__
                rowCount = h-2-GetVisibleStringWidth(number);
#else
                rowCount = h-2-GetVisibleStringWidth(number, font);
#endif
                if (rowCount>headerWidth) {
                    DisplayText (number,slotHeight-3,rowCount,true);
                }
            }
        }
        SetColor (blackC);
        r.left = 0;
        r.right = visWidth;


        if (settings.width&HY_COMPONENT_BORDER_B) {
            r.top = r.bottom = visHeight-1;
            if (setport) {
                DrawLine(r);
            }
        }

    }
    if (setport) {
        EndDraw();
    }
}

//__________________________________________________________________
void    _HYSequencePane::HScrollPane (long dx)
{
#ifndef __HYPHY_GTK__
    if (abs(dx)>=endColumn-startColumn-5)
#endif
    {
        startColumn+=dx;
        BuildPane();
        _MarkForUpdate();
        return;
    }

    long visWidth = _HYCanvas::GetMaxW()-headerWidth-5,
         visHeight = _HYCanvas::GetMaxH()-5,
         slotHeight = GetSlotHeight(),h,v,lastColor = -1,adx,
         loopStart,loopEnd,selectionIndex = 0;

    if (settings.width&HY_COMPONENT_V_SCROLL) {
        visWidth-=HY_SCROLLER_WIDTH;
    }
    if (settings.width&HY_COMPONENT_H_SCROLL) {
        visHeight-=HY_SCROLLER_WIDTH;
    }

    if (dx>0) {
        if (endColumn+dx>columnStrings.lLength) {
            dx = columnStrings.lLength-endColumn;
        }
        adx = dx;
        v = (startColumn+dx)/blockWidth-startColumn/blockWidth;
        if (startColumn&&(startColumn%blockWidth==0)) {
            v++;
        }
        if ((startColumn+dx)&&((startColumn+dx)%blockWidth==0)) {
            v--;
        }
    } else {
        if  (startColumn<-dx) {
            dx=-startColumn;
        }
        adx = -dx;
        v = startColumn/blockWidth-(startColumn+dx)/blockWidth;
        if (startColumn&&(startColumn%blockWidth==0)) {
            v--;
        }
        if ((startColumn+dx)&&((startColumn+dx)%blockWidth==0)) {
            v++;
        }
    }
    if (!dx) {
        return;
    }

    StartDraw();

    h = adx*charWidth+2*v;
    _HYRect     r = {0,0,visHeight,0,1},
                backCharRect;

    visWidth+=headerWidth;

    startColumn+=dx;
    endColumn+=dx;

    if (dx>0) {
        r.left = headerWidth+1;
        r.right = visWidth;
        _SlideRect (r,0,-h);
    } else {
        r.left = headerWidth+1;
        r.right = visWidth;
        _SlideRect (r,0,h);
    }


    if (dx>0) {
        r.right = visWidth;
        r.left = r.right-h;
    } else {
        r.left = headerWidth+1;
        r.right = headerWidth+h+1;
    }
    SetColor    (backColor);
    FillRect    (r);
    _HYColor    blackC = {0,0,0},
                charColor = {0,0,0},
                selectColor = highlightColor;

    r.top = 0;
    if (numbers) {
        r.bottom = (GetSlotHeight()+1);
        SetColor(headerColor);
        FillRect(r);
        r.top = r.bottom;
    } else {
        r.top = r.bottom = 0;
    }
    SetColor    (blackC);

    if (settings.width&HY_COMPONENT_BORDER_T) {
        DrawLine    (r);
    }

    if (settings.width&HY_COMPONENT_BORDER_B) {
        r.top = r.bottom = visHeight-1;
        DrawLine (r);
    }
    r.top = 0;
    if (dx<0) {
        if (numbers) {
            SetColor (headerColor);
            r.bottom = (GetSlotHeight()+1);
            r.left = 0;
            r.right = headerWidth;
            FillRect (r);
        }
        h = HY_SEQUENCE_PANE_CHAR_SPACING/2+headerWidth;
    } else {
        h = HY_SEQUENCE_PANE_CHAR_SPACING/2+headerWidth;
        for (v=startColumn; v<endColumn-dx; v++,h+=charWidth)
            if (v&&(v%blockWidth)==0) {
                h+=2;
            }
    }
    r.bottom = visHeight;
    loopStart = dx>0?endColumn-dx:startColumn;
    loopEnd =  dx>0?endColumn:startColumn-dx;

    if (loopStart<0) {
        loopStart = 0;
    }
    if (loopEnd<0) {
        loopEnd = 0;
    }
    if (loopStart>columnStrings.lLength) {
        loopStart = columnStrings.lLength;
    }
    if (loopEnd>columnStrings.lLength) {
        loopEnd = columnStrings.lLength;
    }

    for (long c = loopStart; c < loopEnd; c++,h+=charWidth) {
        bool isColumnSelected = false;
        if (selectionIndex<selection.lLength) {
            visWidth = selection.lData[selectionIndex];
            while ((visWidth<c)&&(selectionIndex<selection.lLength)) {
                selectionIndex++;
                visWidth = selection.lData[selectionIndex];
            }
            if (c==visWidth) {
                if (!invertColors) {
                    _HYRect invR;
                    invR.top = (GetSlotHeight()+1)+1;
                    invR.bottom = visHeight-1;
                    invR.right = h+charWidth;
                    invR.left = h-2;
                    if (c&&(c%blockWidth==0)) {
                        invR.right+=2;
                    }
                    /*else
                        if ((c+1==loopEnd)&&(dx>0))
                            invR.right+=3;
                    */
                    charColor = GetColor();
                    SetColor (selectColor);
                    FillRect (invR);
                    SetColor (charColor);
                }
                selectionIndex++;
                isColumnSelected = true;
            }
        }
        if (c&&(c%blockWidth==0)) {
            SetColor (blackC);
            if (numbers) {
                _String number (c);

#ifdef __MAC__
                visWidth = h-2-GetVisibleStringWidth(number);
#else
                visWidth = h-2-GetVisibleStringWidth(number, font);
#endif

                if (visWidth>headerWidth) {
                    DisplayText (number,slotHeight-3,visWidth,true);
                }
            }
            r.left = r.right = h;
            DrawLine (r);
            SetColor (charColor);
            h+=2;
        }
        if (numbers) {
            v = slotHeight+(GetSlotHeight()+1);
        } else {
            v = slotHeight;
        }

        _String *thisString = (_String*)columnStrings(c);
        unsigned char topChar = showDots?thisString->sData[speciesIndex.lData[0]]:0;

        backCharRect = (_HYRect) {
            v-slotHeight+1,h-1,v+1,h+charWidth-1,0
        };
        if ((c+1)%blockWidth == 0) {
            backCharRect.right++;
        }

        for (long r = startRow; r < endRow; r++,v+=slotHeight) {
            /*unsigned char thisC = thisString->sData[speciesIndex.lData[r]];
            if (r&&(thisC==topChar))
                thisC = '.';
            if (colorMap[thisC]!=lastColor)
            {
                lastColor = colorMap[thisC];
                charColor = LongToHYColor(characterColors.lData[lastColor]);
                SetColor (charColor);
            }
            DisplayChar (thisC,v,h);*/
            unsigned char thisC = thisString->sData[speciesIndex.lData[r]];
            long     myColor    = colorMap[thisC];

            if (r && thisC==topChar) {
                thisC = '.';
            }

            if (invertColors && !isColumnSelected) {
                charColor = LongToHYColor (characterColors.lData[myColor]);
                if ((long)charColor.R+charColor.G+charColor.B > 100) {
                    SetColor (charColor);
                    FillRect    (backCharRect);
                }
                backCharRect.top += slotHeight;
                backCharRect.bottom += slotHeight;
                SetColor    (blackC);
            } else if (colorMap[thisC]!=lastColor) {
                lastColor = colorMap[thisC];
                charColor = LongToHYColor(characterColors.lData[lastColor]);
                SetColor (charColor);
            }
            DisplayChar (thisC,v,h);
        }
    }

    if (numbers&&(dx<0)) {
        loopEnd--;
        v = loopEnd%blockWidth;
        visWidth = blockWidth-log(double(loopEnd+blockWidth-1))/log((double)blockWidth);

        if (adx<blockWidth) {
            loopStart = (startColumn-dx)%blockWidth;
            if ((startColumn-dx)/blockWidth!=startColumn/blockWidth)
                // scrolled thru a block divider
            {
                if (loopStart==0) {
                    h += (blockWidth-v-1)*charWidth;
                } else {
                    h = 0;
                }
            } else {
                if (loopStart>visWidth) {
                    h += (blockWidth-v-1)*charWidth;
                } else {
                    h = 0;
                }
            }
        } else {
            loopStart = (startColumn-dx)%blockWidth;

            if (loopStart)
                if (loopStart>=visWidth) {
                    h -= (loopStart-blockWidth)*charWidth;
                } else {
                    h = 0;
                }

        }

        if (h) {
            _String number (v?((loopEnd/blockWidth)+1)*blockWidth:loopEnd);
            SetColor (blackC);
#ifdef __MAC__
            visWidth = h-2-GetVisibleStringWidth(number);
#else
            visWidth = h-2-GetVisibleStringWidth(number, font);
#endif
            while ((visWidth < headerWidth)&&(number.sLength)) {
                number.Trim (1,-1);
#ifdef __MAC__
                visWidth = h-2-GetVisibleStringWidth(number);
#else
                visWidth = h-2-GetVisibleStringWidth(number, font);
#endif
            }
            if (number.sLength) {
                DisplayText (number,slotHeight-3,visWidth,true);
            }
        }
    }

    EndDraw();
    if (messageRecipient) {
        messageRecipient->ProcessEvent(generateScrollEvent(0,0));
    }
    _MarkForUpdate();
}
//__________________________________________________________________
void    _HYSequencePane::ScrollToHSelection (void)
{
    if (selection.lLength) {
        long        dx = selection.lData[0]-startColumn;
        if ((dx<0)&&(-dx>startColumn)) {
            dx = -startColumn;
        } else if ((dx>0)&&(endColumn+dx>columnStrings.lLength)) {
            dx = columnStrings.lLength-endColumn;
        }

        HScrollPane (dx);
    }
}


//__________________________________________________________________
bool    _HYSequencePane::ProcessEvent(_HYEvent* e)
{
    if (e->EventClass() == _hyScrollingEvent) {
        long h,v,k,t;
        _String firstArg = e->EventCode().Cut (0,(v=e->EventCode().Find(','))-1);
        h = firstArg.toNum();
        firstArg = e->EventCode().Cut (v+1,-1);
        DeleteObject(e);
        v = firstArg.toNum();
        if (h||v) {
            if (h) {
                if (settings.width&HY_COMPONENT_H_SCROLL) {
                    k = columnStrings.lLength-endColumn+startColumn+1;
                    t = _GetHScrollerPos()*(_Parameter)k/MAX_CONTROL_VALUE;
                    if (startColumn == t) {
                        if (h>0) {
                            t++;
                        } else {
                            t--;
                        }
                        _SetHScrollerPos((MAX_CONTROL_VALUE*(_Parameter)t)/k);
                    }
                    if (!v) {
                        HScrollPane (t-startColumn);
                        return true;
                    }
                    startColumn = t;
                } else {
                    HScrollPane (h);
                }
                //startColumn += v;

            } else {
                if (settings.width&HY_COMPONENT_V_SCROLL) {
                    k = RowCount()-endRow+startRow+1;
                    t = (_GetVScrollerPos()*(_Parameter)k) /MAX_CONTROL_VALUE;
                    if (startRow != t) {
                        startRow = t;
                    } else {
                        if (v>0) {
                            startRow++;
                        } else {
                            startRow--;
                        }
                        _SetVScrollerPos((MAX_CONTROL_VALUE*(_Parameter)startRow)/k);
                    }
                } else {
                    startRow += v;
                }
                BuildPane();
            }
            if (messageRecipient) {
                messageRecipient->ProcessEvent (generateScrollEvent(0,0));
            }

            _MarkForUpdate();
        }
        return true;
    }
    DeleteObject (e);
    return false;
}

//__________________________________________________________________
int         _HYSequencePane::GetMaxW (void)
{
    int res = charWidth*(columnStrings.lLength)+headerWidth;
    res+=(columnStrings.lLength/blockWidth)*2;
    if (settings.width&HY_COMPONENT_V_SCROLL) {
        res+=2*HY_SCROLLER_WIDTH;
    }
    if (res<settings.right) {
        res = settings.right;
    }
    return res;
}

//__________________________________________________________________
int         _HYSequencePane::GetMaxH (void)
{
    int res = (GetFont().size+HY_SEQUENCE_PANE_CHAR_SPACING)*(RowCount())+3;
    if (settings.width&HY_COMPONENT_H_SCROLL) {
        res+=2*HY_SCROLLER_WIDTH;
    } else {
        res+=2;
    }
    if (numbers) {
        res += (GetSlotHeight()+1);
    }

    return res;
}

//__________________________________________________________________
int         _HYSequencePane::GetMinH (void)
{
    int res = GetFont().size+HY_SEQUENCE_PANE_CHAR_SPACING*2;

    if (numbers) {
        res += (GetSlotHeight()+1)+HY_SEQUENCE_PANE_CHAR_SPACING+2;
    }

    return res;
}
//__________________________________________________________________
long            _HYSequencePane::RowCount (void)
{
    int res = 0;
    if (columnStrings.lLength) {
        res = speciesIndex.lLength;
    }
    return res;
}
//__________________________________________________________________
long            _HYSequencePane::GetSlotHeight (void)
{
    return GetFont().size+HY_SEQUENCE_PANE_CHAR_SPACING;
}

//__________________________________________________________________

void        _HYSequencePane::SetVisibleSize (_HYRect r)
{
    int ht = r.bottom-r.top+5,
        wd = r.right-r.left+5;

    SetPaneSize (ht,wd,32);
    hOrigin = vOrigin = 0;
    settings.right = hSize = wd;
    settings.bottom = vSize = ht;
    needUpdate = true;
//  check if resizing also requires scrolling

    long visWidth  = _HYCanvas::GetMaxW()-headerWidth-5,
         visHeight = _HYCanvas::GetMaxH()-5,
         rowCount = speciesIndex.lLength;

    if (settings.width&HY_COMPONENT_V_SCROLL) {
        visWidth-=HY_SCROLLER_WIDTH;
    }
    if (settings.width&HY_COMPONENT_H_SCROLL) {
        visHeight-=HY_SCROLLER_WIDTH;
    }

    endColumn = startColumn+visWidth/charWidth;
    ht = (endColumn-startColumn)/5;
    endColumn = startColumn+(visWidth-ht)/charWidth;
    endRow = startRow+visHeight/GetSlotHeight();

    if (endColumn>columnStrings.lLength) {
        startColumn -= endColumn-columnStrings.lLength-1;
        if (startColumn<0) {
            startColumn = 0;
        }
        endColumn = columnStrings.lLength;
    }

    if (endRow>rowCount) {
        startRow -= endRow-rowCount-1;
        if (startRow<0) {
            startRow = 0;
        }
        endRow = rowCount;
    }

    if (settings.width&HY_COMPONENT_V_SCROLL) {
        long k = RowCount()-endRow+startRow+1;
        _SetVScrollerPos((MAX_CONTROL_VALUE*(_Parameter)startRow)/k);
    }

    BuildPane();
    r.bottom ++;
    _HYComponent::SetVisibleSize(r);
    rel.bottom--;
}


//__________________________________________________________________

void        _HYSequencePane::SetBlockWidth (long bw, bool update)
{
    if (bw!=blockWidth) {
        blockWidth = bw;
        if (update) {
            BuildPane();
        }
    }
}

//__________________________________________________________________

void        _HYSequencePane::SetHeaders (_List* newH, bool update)
{
    if (newH) {
        rowHeaders.Clear();
        rowHeaders.Duplicate (newH);
        speciesIndex.Clear();
    }

    headerWidth = 0;
    shortHeaderWidth = 0;
    StartDraw();
    for (long k=0; k<rowHeaders.lLength; k++) {
        _String * thisString = (_String*)rowHeaders(k);
#ifdef __MAC__
        long lW = GetVisibleStringWidth(*thisString);
#else
        long lW = GetVisibleStringWidth(*thisString, font);
#endif
        if (lW>headerWidth) {
            headerWidth = lW;
        }
        if (thisString->sLength<=HY_SEQUENCE_PANE_SHORT_WIDTH) {
            if (lW>shortHeaderWidth) {
                shortHeaderWidth = lW;
            }
        } else {
            _String chopped (*thisString);
            chopped.Trim(0,HY_SEQUENCE_PANE_SHORT_WIDTH-1);
#ifdef __MAC__
            long sW = GetVisibleStringWidth (chopped);
#else
            long sW = GetVisibleStringWidth (chopped, font);
#endif
            if (sW>shortHeaderWidth) {
                shortHeaderWidth = sW;
            }
        }
        if (newH) {
            speciesIndex  << k;
        }
    }
    headerWidth         += HY_SEQUENCE_PANE_CHAR_SPACING;
    shortHeaderWidth    += HY_SEQUENCE_PANE_CHAR_SPACING;
    EndDraw();
    fullHeaderWidth = headerWidth;

    if (update) {
        BuildPane();
    }
}

//__________________________________________________________________

void        _HYSequencePane::ChangeHeaders (void)
{
    BuildPane();
    _MarkForUpdate();
    SendSelectionChange(true);
}

//__________________________________________________________________

void        _HYSequencePane::SetNameDisplayMode (unsigned char newMode, bool update)
{
    if (nameDisplayFlags!=newMode) {
        nameDisplayFlags = newMode;
        if (newMode&HY_SEQUENCE_PANE_NAMES_SHORT) {
            headerWidth = shortHeaderWidth;
        } else if (newMode&HY_SEQUENCE_PANE_NAMES_ALL) {
            headerWidth = fullHeaderWidth;
        } else {
            headerWidth = 0;
        }

        if (update) {
            BuildPane();
            _MarkForUpdate();
        }
    }
}

//__________________________________________________________________
void        _HYSequencePane::AddColumnToSelection (long c)
{
    if ((c>0)&&(c<columnStrings.lLength)) {
        if (selection.lLength) {
            long f = selection.BinaryFind (c);
            if (f<0) {
                selection.InsertElement ((BaseRef)c,-f-2,false,false);
            }
        } else {
            selection<<c;
        }
    }
}

//__________________________________________________________________
void        _HYSequencePane::AddSpeciesToSelection (long c)
{
    if (c>0 && c<speciesIndex.lLength)
        if (vselection.lLength) {
            long f = vselection.BinaryFind (c);
            if (f<0) {
                vselection.InsertElement ((BaseRef)c,-f-2,false,false);
            }
        } else {
            vselection<<c;
        }
}


//__________________________________________________________________
void        _HYSequencePane::SelectSequenceNames (_List& list, bool send)
{
    selection.Clear();
    vselection.Clear();

    for (long k=0; k<rowHeaders.lLength; k++) {
        _String * aSeq = (_String*)rowHeaders(speciesIndex.lData[k]);
        if (list.BinaryFind (aSeq) >= 0) {
            vselection << k;
        }
    }

    BuildPane();
    _MarkForUpdate();
    if (send) {
        SendSelectionChange (true);
    }
}


//__________________________________________________________________
void        _HYSequencePane::SelectAll (bool flag)
{
    if (flag) {
        if (selection.lLength!=columnStrings.lLength) {
            selection.Clear();
            vselection.Clear();
            selection.RequestSpace (columnStrings.lLength);
            for (long k=0; k<columnStrings.lLength; k++) {
                selection<<k;
            }
            BuildPane();
            _MarkForUpdate();
            SendSelectionChange();
        }
    } else {
        if (selection.lLength) {
            selection.Clear();
            BuildPane();
            _MarkForUpdate();
            SendSelectionChange();
        }
    }
}

//__________________________________________________________________
void        _HYSequencePane::SelectRange (_SimpleList& range, bool vert)
{
    if (vert) {
        selection.Clear();
        vselection.Clear();
        vselection.Duplicate(&range);
        vselection.Sort();
    } else {
        vselection.Clear();
        selection.Clear();
        selection.Duplicate(&range);
        selection.Sort();
    }
    BuildPane();
    _MarkForUpdate();
}

//__________________________________________________________________
void        _HYSequencePane::MoveSpecies (long oldIndex, long newIndex)
{
    if ((newIndex!=oldIndex)&&(newIndex>=-1)) {
        long c = speciesIndex.lData[oldIndex];
        speciesIndex.InsertElement ((BaseRef)c,newIndex+1,false,false);
        if (oldIndex>newIndex) {
            oldIndex++;
            newIndex++;
        }
        vselection.lData[0] = newIndex;
        speciesIndex.Delete(oldIndex);
        BuildPane();
        _MarkForUpdate();
    }
}

//__________________________________________________________________
void        _HYSequencePane::SetSequenceOrder (_SimpleList& no)
{
    speciesIndex.Duplicate (&no);
    BuildPane();
    _MarkForUpdate();
}

//__________________________________________________________________
void        _HYSequencePane::ProcessSelectionChange (long h, long v, bool shift, bool command, bool drag)
{
    if (vselection.lLength) {
        vselection.Clear();
        SendSelectionChange(true);
    }
    if (drag) {
        long k = columnStrings.lLength-endColumn+startColumn+1,
             scrollAmount = 0;

        if ( h<headerWidth && v>GetSlotHeight()+1)
            // scroll in columns from the left
        {
            if (selection.lData[0]>startColumn) {
                h = selection.lData[0];
                for (v=h-1; v>=startColumn; v--) {
                    selection.InsertElement((BaseRef)v,0,false,false);
                }
                BuildPane();
                if (startColumn == 0) {
                    _MarkForUpdate();
                }
            }
            if (startColumn>0) {
                long dsc = abs(dragScrollCounter);

                if (dsc < 25) {
                    scrollAmount = 1;
                } else if (dsc < 50) {
                    scrollAmount = 2;
                } else if (dsc < 150) {
                    scrollAmount = 3;
                } else {
                    scrollAmount = 5;
                }

                if (scrollAmount>startColumn) {
                    scrollAmount = startColumn;
                }

                if (selection.lData[0] >= startColumn) {
                    for (dsc = startColumn-1; dsc>=startColumn-scrollAmount; dsc--) {
                        selection.InsertElement((BaseRef)(dsc),0,false,false);
                    }
                } else {
                    if (selection.lLength) {
                        if (selection.lLength < scrollAmount) {
                            scrollAmount = selection.lLength;
                        }

                        if (selection.lData[selection.lLength] >= startColumn) {
                            while (selection.lLength && (selection.lData[selection.lLength] >= startColumn)) {
                                selection.Delete (selection.lLength-1);
                            }
                            BuildPane();
                        }

                        for (dsc = 0; (dsc<scrollAmount)&&(selection.lLength>1); dsc++) {
                            selection.Delete (selection.lLength-1);
                        }
                    }
                }


                HScrollPane(-scrollAmount);

                dragScrollCounter -= scrollAmount;

                _SetHScrollerPos(((double)MAX_CONTROL_VALUE*startColumn)/k);
                SendSelectionChange();

                return;
            }
        } else {
            if (h>=_HYCanvas::GetMaxW()-5 && v>GetSlotHeight()+1)
                // scroll in columns from the right
            {
                if (selection.lData[selection.lLength-1]<endColumn-1) {
                    for (v=selection.lData[selection.lLength-1]+1; v<endColumn; v++) {
                        selection<<v;
                    }
                    BuildPane();
                    if (endColumn==columnStrings.lLength) {
                        _MarkForUpdate();
                    }
                }
                if (endColumn<columnStrings.lLength) {
                    long dsc = abs(dragScrollCounter);

                    if (dsc < 25) {
                        scrollAmount = 1;
                    } else if (dsc < 50) {
                        scrollAmount = 2;
                    } else if (dsc < 150) {
                        scrollAmount = 3;
                    } else {
                        scrollAmount = 5;
                    }

                    if (scrollAmount>columnStrings.lLength-endColumn) {
                        scrollAmount = columnStrings.lLength-endColumn;
                    }


                    if ((!selection.lLength)||(selection.lData[selection.lLength-1]<=endColumn))
                        for (dsc = endColumn; dsc<endColumn+scrollAmount; dsc++) {
                            selection << dsc;
                        }
                    else {
                        if (selection.lLength) {
                            if (selection.lLength < scrollAmount) {
                                scrollAmount = selection.lLength;
                            }

                            if (selection.lData[0] < endColumn) {
                                while (selection.lLength && (selection.lData[0] < endColumn)) {
                                    selection.Delete (0);
                                }
                                BuildPane();
                            }

                            for (dsc = 0; (dsc<scrollAmount)&&(selection.lLength>1); dsc++) {
                                selection.Delete (0);
                            }
                        }
                    }


                    HScrollPane(scrollAmount);

                    dragScrollCounter+= scrollAmount;

                    _SetHScrollerPos(((double)MAX_CONTROL_VALUE*startColumn)/k);
                    SendSelectionChange();
                    return;
                }
            } else {
                dragScrollCounter = 0;
            }
        }
    } else {
        dragScrollCounter = 0;
    }

    if ( h>headerWidth && v>GetSlotHeight()+1 ) {
        v=startRow+(v-(GetSlotHeight()+1))/(GetFont().size+HY_SEQUENCE_PANE_CHAR_SPACING);

        long    k,
                p = headerWidth+charWidth+HY_SEQUENCE_PANE_CHAR_SPACING/2;

        for (k=startColumn; k<endColumn; k++,p+=charWidth)
            // find which column was clicked
        {
            if (k&&(k%blockWidth==0)) {
                p+=2;
            }
            if (p>=h) {
                break;
            }
        }
        h = k;
        if (shift) {
            if (recentClick!=-1) {
                if (h==recentClick) {
                    return;
                }

                _SimpleList   * saveSelection = nil;

                if ( !command || drag) {
                    saveSelection = (_SimpleList*)selection.makeDynamic();
                    selection.Clear();
                }

                if (h>recentClick) {
                    k=recentClick;
                    p=h;
                } else {
                    p=recentClick;
                    k=h;
                }
                //if (!drag)
                {
                    if (!command||!selection.lLength||drag)
                        for (v=k; v<=p; v++) {
                            selection<<v;
                        }
                    else {
                        for (v=k; v<=p; v++) {
                            AddColumnToSelection(v);
                        }
                    }
                }
                if (saveSelection) {
                    bool doUpdate = false;
                    if (!saveSelection->Equal(selection)) {
                        doUpdate = true;
                    }
                    DeleteObject (saveSelection);
                    if (!doUpdate) {
                        return;
                    }
                }
                BuildPane();
                _MarkForUpdate();
                SendSelectionChange();
                return;
            }
        } else {
            if (command) {
                v = selection.BinaryFind (h);
                if (v>=0) {
                    selection.Delete(v);
                } else {
                    selection.InsertElement ((BaseRef)h,-v-2,false,false);
                }
                BuildPane();
                _MarkForUpdate();
                SendSelectionChange();
                return;
            }
        }
        //if (selection.BinaryFind(h)>0) return;
        selection.Clear();
        selection<<h;
        recentClick = h;
        BuildPane();
        _MarkForUpdate();
    } else if ((h<headerWidth)&&(v<(GetSlotHeight()+1))) {
        if ((!drag)&&selection.lLength) {
            selection.Clear();
            vselection.Clear();
            BuildPane();
            _MarkForUpdate();
        }
    }
    SendSelectionChange();
}
//__________________________________________________________________
void        _HYSequencePane::AlphabetizeSpecies (void)
{
    _List currentActiveNames;
    long  k;

    for (k=0; k<speciesIndex.lLength; k++) {
        currentActiveNames << rowHeaders (speciesIndex.lData[k]);
    }

    SortLists (&currentActiveNames,&speciesIndex);

    BuildPane ();
    _MarkForUpdate();
}

//__________________________________________________________________
void        _HYSequencePane::RevertFileOrder (void)
{
    speciesIndex.Sort();
    BuildPane ();
    _MarkForUpdate();
}

//__________________________________________________________________
void        _HYSequencePane::CleanUpSequenceNames (void)
{
    bool    doSomething = false;

    _List       namesl;
    _AVLList    names (&namesl);

    for    (long k=0; k<speciesIndex.lLength; k++) {
        _String * thisString = (_String*)rowHeaders (speciesIndex.lData[k]);
        if (!thisString->IsValidIdentifier(false)) {
            BufferToConsole ("Changed ");
            StringToConsole(*thisString);
            thisString->ConvertToAnIdent(false);
            BufferToConsole (" to ");
            StringToConsole(*thisString);
            NLToConsole();
            doSomething = true;
        }

        _String * testString = new _String (*thisString);

        if (!testString) {
            checkPointer (testString);
        }

        long    tryThisSuffix = 2;

        while (names.Find (testString)>=0) {
            *testString = *thisString & '_' & tryThisSuffix;
            tryThisSuffix++;
        }

        if (tryThisSuffix>2) {
            BufferToConsole ("Changed ");
            StringToConsole(*thisString);
            BufferToConsole (" to ");
            StringToConsole(*testString);
            BufferToConsole (" to avoid duplicate identifiers\n");
            doSomething  = true;
            thisString->CopyDynamicString (testString,true);
        } else {
            DeleteObject (testString);
        }

        names.Insert(thisString);
        thisString->nInstances++;
    }

    if (doSomething) {
        SetHeaders (nil,true);
        _MarkForUpdate();
    }
}

//__________________________________________________________________
void        _HYSequencePane::BatchRenameSequences (_List& oldNames,_List& newNames)
{
    bool touched = false;
    for (long k=0; k<oldNames.lLength; k++) {
        long nID = rowHeaders.Find (oldNames(k));
        if (nID >= 0 && rowHeaders.Find (newNames(k)) < 0) {
            rowHeaders (nID)->Duplicate (newNames(k));
            touched = true;
        }
    }
    if (touched) {
        SetHeaders (nil,true);
        _MarkForUpdate();
    }
}


//__________________________________________________________________
void        _HYSequencePane::EditSequenceName (long k)
{
    _String prompt ("Edit Sequence Name"),
            *present = ((_String*)rowHeaders (speciesIndex.lData[k])),
             edited = *present;

    if (EnterStringDialog ( edited, prompt,(Ptr)messageRecipient)) {
        if (!edited.Equal (present)) {
            if (rowHeaders.Find (&edited)>=0) {
                prompt = _String("Another sequence is already named ") & edited & ". Please choose another name.";
                ProblemReport (prompt, (Ptr)messageRecipient);
                return;
            }
            present->Duplicate (&edited);
            SetHeaders (nil,true);
            _MarkForUpdate();
        }
    }
}

//__________________________________________________________________
void        _HYSequencePane::ProcessVSelectionChange (long h, long v, bool shift, bool command, bool drag, bool editName)
{
    if (selection.lLength) {
        selection.Clear();
        BuildPane();
        _MarkForUpdate();
        SendSelectionChange();
    }

    if (drag) {
        return;
    }

    if ((h<headerWidth)&&(v>(GetSlotHeight()+1))) {
        long k, p;
        v = (v-(GetSlotHeight()+1))/(GetFont().size+HY_SEQUENCE_PANE_CHAR_SPACING)+startRow;
        if (v>=endRow) {
            return;
        }

        if (editName) {
            EditSequenceName (v);
            return;
        }

        if (shift) {
            if (recentClick!=-1) {
                if (v==recentClick) {
                    return;
                }
                if ((!command)||drag) {
                    vselection.Clear();
                }
                if (v>recentClick) {
                    k=recentClick;
                    p=v;
                } else {
                    p=recentClick;
                    k=v;
                }
                if (!command||!vselection.lLength||drag)
                    for (h=k; h<=p; h++) {
                        vselection<<h;
                    }
                else {
                    for (h=k; h<=p; h++) {
                        AddSpeciesToSelection(h);
                    }
                }
                StartDraw();
                BuildHeaders();
                EndDraw();
                _MarkForUpdate();
                SendSelectionChange();
                return;
            }
        } else {
            if (command) {
                h = vselection.BinaryFind (v);
                if (h>=0) {
                    vselection.Delete(h);
                } else {
                    vselection.InsertElement ((BaseRef)v,-h-2,false,false);
                }
                StartDraw();
                BuildHeaders();
                EndDraw();
                _MarkForUpdate();
                SendSelectionChange(true);
                return;
            }
        }
        if (vselection.BinaryFind(v)>0) {
            return;
        }
        vselection.Clear();
        vselection<<v;
        recentClick = v;
        StartDraw();
        BuildHeaders();
        EndDraw();
        _MarkForUpdate();
    }

    SendSelectionChange(true);
}



//__________________________________________________________________

void        _HYSequencePane::SendSelectionChange (bool vert)
{
    if (messageRecipient) {
        messageRecipient->ProcessEvent (generateMenuSelChangeEvent(GetID(),vert));
    }
}

//__________________________________________________________________

void        _HYSequencePane::ProcessContextualPopUp (long l, long t)
{
    if (messageRecipient) {
        messageRecipient->ProcessEvent (generateContextPopUpEvent(GetID(),l,t));
    }
}




//EOF
