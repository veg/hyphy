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
#include "HYGraphicPane.h"
#include "HYUtils.h"
#include "HYComponent.h"
#include "errorfns.h"
#include "math.h"

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

_HYColor hyDefaultFrameColor = {100,100,100};

//__________________________________________________________________

_HYGraphicPane::_HYGraphicPane(int ht, int wd, int d):
    _HYPlatformGraphicPane (ht,wd,d)

{
    w = wd;
    h = ht;
    color.R=color.G=color.B=0;
    bColor.R=bColor.G=bColor.B=255;
    StartDraw();
    SetColor  (color);
    SetBColor (bColor);
#ifdef __WINDOZE__
    font.face = "Verdana";
#else
#ifdef __MAC__
    font.face = "Monaco";
#else
    font.face = _HY_SANS_FONT;
#endif
#endif
    font.size = 10;
    font.style = HY_FONT_PLAIN;
    SetFont (font);
    EraseAll();
    EndDraw();
    depth = d;
}


//__________________________________________________________________

_HYGraphicPane::~_HYGraphicPane(void)
{
}

//__________________________________________________________________

void _HYGraphicPane::SetPaneSize  (int ht,int wd, int d)
{
    if ((h==ht)&&(w==wd)&&(d==depth)) {
        return;
    }
    h = ht;
    w = wd;
    _SetPaneSize (h,w,d);
    StartDraw();
    EraseAll();
    /*font.face = "Monaco";
    font.size = 10;
    font.style = HY_FONT_PLAIN;*/
    SetColor (color);
    SetFont (font);
    EndDraw();
    depth = d;
}


//__________________________________________________________________

void _HYGraphicPane::DrawLine (_HYRect lineDesc)
{
    _DrawLine (lineDesc);
}

//__________________________________________________________________

void _HYGraphicPane::DrawHatchedLine (_HYRect lineDesc)
{
    _DrawHatchedLine (lineDesc);
}

//__________________________________________________________________
void _HYGraphicPane::DisplayText     (_String theText,int t, int l, bool dir)
{
    _DisplayText (theText,t,l,dir);
}

//__________________________________________________________________
void _HYGraphicPane::DisplayText     (_String theText,_HYRect r, char align)
{
    _DisplayText (theText,r,align);
}

//__________________________________________________________________
void _HYGraphicPane::DisplayChar     (char c,int t, int l)
{
    _DisplayChar (c,t,l);
}

//__________________________________________________________________

void _HYGraphicPane::DrawRect (_HYRect rectDesc)
{
    _DrawRect (rectDesc);
}
//__________________________________________________________________

void _HYGraphicPane::FillRect (_HYRect rectDesc)
{
    _FillRect (rectDesc);
}
//__________________________________________________________________

void _HYGraphicPane::InvertRect (_HYRect rectDesc)
{
    _InvertRect (rectDesc);
}

//__________________________________________________________________

void _HYGraphicPane::EraseRect (_HYRect rectDesc)
{
    _EraseRect (rectDesc);
}

//__________________________________________________________________

void _HYGraphicPane::DrawOval (_HYRect rectDesc)
{
    _DrawOval (rectDesc);
}
//__________________________________________________________________

void _HYGraphicPane::FillOval (_HYRect rectDesc)
{
    _FillOval (rectDesc);
}

//__________________________________________________________________

void _HYGraphicPane::EraseOval (_HYRect rectDesc)
{
    _EraseOval (rectDesc);
}

//__________________________________________________________________

void _HYGraphicPane::DrawPolygon (_SimpleList& r, long width)
{
    Ptr p = _DefinePolygon (r);
    _DrawPolygon (p, width);
    _KillPolygon (p);
}
//__________________________________________________________________

void _HYGraphicPane::FillPolygon (_SimpleList& r)
{
    Ptr p = _DefinePolygon (r);
    _FillPolygon (p);
    _KillPolygon (p);
}
//__________________________________________________________________

void _HYGraphicPane::ErasePolygon (_SimpleList& r)
{
    Ptr p = _DefinePolygon (r);
    _ErasePolygon (p);
    _KillPolygon (p);
}

//__________________________________________________________________

void _HYGraphicPane::DrawArc (_HYRect rectDesc, int s, int f)
{
    _DrawArc (rectDesc,s,f);
}
//__________________________________________________________________

void _HYGraphicPane::FillArc (_HYRect rectDesc, int s, int f)
{
    _FillArc (rectDesc,s,f);
}

//__________________________________________________________________

void _HYGraphicPane::EraseArc (_HYRect rectDesc, int s, int f)
{
    _EraseArc (rectDesc,s,f);
}

//__________________________________________________________________

void _HYGraphicPane::EraseAll (void)
{
    _HYRect all = {0,0,h,w,0};
    EraseRect (all);
}

//__________________________________________________________________

void _HYGraphicPane::SetColor  (_HYColor c)
{
    color = c;
    _SetColor (c);
}

//__________________________________________________________________

void _HYGraphicPane::SetBColor  (_HYColor c)
{
    bColor = c;
    _SetBColor (c);
}

//__________________________________________________________________

void _HYGraphicPane::DrawFrame  (_HYColor c)
{
    SetColor(c);
    _HYRect frame;
    frame.left=frame.top=0;
    frame.bottom = h;
    frame.right = w;
    frame.width = 2;
    DrawRect (frame);
    _HYColor darker;
    darker.R = c.R/2;
    darker.G = c.G/2;
    darker.B = c.B/2;
    frame.left=frame.top=2;
    frame.bottom = h-2;
    frame.right = w-2;
    frame.width = 1;
    SetColor(darker);
    DrawRect (frame);
}


//__________________________________________________________________

void _HYGraphicPane::SetFont (_HYFont f)
{
    font = f;
    _SetFont (f);
}

//__________________________________________________________________

void _HYGraphicPane::SetFontSize (long s)
{
    font.size = s;
    _SetFontSize (s);
}

//__________________________________________________________________

void _HYGraphicPane::StartDraw   (void)
{
    _StartDraw();
}

//__________________________________________________________________

void _HYGraphicPane::EndDraw     (void)
{
    _EndDraw();
}

//__________________________________________________________________

void _HYGraphicPane::DrawTable   (_List& theStrings, _List& p, char align, _HYRect margins, _HYRect m)
{
    long i1, i2;
    for (i1=0; i1<theStrings.lLength; i1++) {
        _List       *thisRow = (_List*)theStrings(i1);
        _SimpleList *rw      = (_SimpleList*)p(i1);

        for (i2=0; i2<thisRow->lLength; i2++) {
            long w = rw->lData[3*i2+1];
            switch (align) {
            case HY_ALIGN_LEFT:
                break;
            case HY_ALIGN_RIGHT:
                w += rw->lData[3*i2+2];
                break;
            default:
                w += rw->lData[3*i2+2]/2;
                break;
            }
            DisplayText (*(_String*)(*thisRow)(i2),m.top+rw->lData[3*i2],m.left+w,true);
            if (margins.width&&(i1==0)) {
                _HYRect betweenColumns = {m.top, m.left+rw->lData[3*i2+1]-margins.left,
                                          m.top+m.bottom, m.left+rw->lData[3*i2+1]-margins.left, margins.width
                                         };
                DrawLine (betweenColumns);
            }
        }
        if (margins.width) {
            if (i1==0) {
                _HYRect betweenColumns = {m.top, m.left+m.right,
                                          m.top+m.bottom, m.left+m.right, margins.width
                                         };
                DrawLine (betweenColumns);
            }
            _HYRect betweenRows = {rw->lData[0]+margins.bottom+m.top,m.left,
                                   rw->lData[0]+margins.bottom+m.top, m.right+m.left,
                                   margins.width
                                  };
            DrawLine (betweenRows);
        }
    }
    if (margins.width) {
        _HYRect betweenRows = {m.top,m.left,
                               m.top, m.right+m.left,
                               margins.width
                              };
        DrawLine (betweenRows);
    }
}

//__________________________________________________________________

void _HYGraphicPane::DrawHashes  (_HYRect b, _SimpleList& off, _List& labels, int hl, int hw)
{
    long    k;
    if (b.right-b.left > b.bottom-b.top) {
        for (k=0; k<labels.lLength; k++) {
            DisplayText (*(_String*)labels(k),b.bottom+b.width+font.size,b.left+off.lData[2*k]-off.lData[2*k+1],true);
            if (hl&&hw) {
                _HYRect     tick = {b.top,b.left+off.lData[2*k]-hw/2,b.top+hl,b.left+off.lData[2*k]-hw/2,hw};
                DrawLine    (tick);
            }
        }
    } else {
        for (k=0; k<labels.lLength; k++) {
            DisplayText (*(_String*)labels(k),b.bottom-off.lData[2*k],b.left-off.lData[2*k+1]-b.width,true);
            if (hl&hw) {
                _HYRect     tick = {b.bottom-off.lData[2*k]-hw/2,b.left-hw,b.bottom-off.lData[2*k]-hw/2,b.left,hw};
                DrawLine    (tick);
            }
        }
    }
}

//__________________________________________________________________
void RotateRect90 (_HYRect& theRect)
{
    int t = theRect.bottom;
    theRect.bottom = theRect.right;
    theRect.right =  t;
    t = theRect.left;
    theRect.left = theRect.top;
    theRect.top = t;
}

//__________________________________________________________________
double DistanceBetweenPoints (double x1, double y1, double x2, double y2)
{
    double d1 = x1-x2,
           d2 = y1-y2;

    return sqrt (d1*d1+d2*d2);
}

//__________________________________________________________________
long AngleBetweenPoints (double x1, double y1, double x2, double y2)
{
    double xspan = x2-x1,
           yspan = y1-y2;

    if (xspan)
        return round(180.0*(0.5-atan(yspan/xspan)/pi_const)
                     + ((xspan< 0.)?180.:0.0));

    return yspan>0?0:180;
}

//__________________________________________________________________
long     ComputeTableCellPlacement (_List& theStrings, _List& thePlacements, _HYRect margins, _HYFont& f)
{
    _SimpleList     columnWidths;
    _List           stringWidths;
    long            i1,i2,w;

    for (i1=0; i1<theStrings.lLength; i1++) {
        _List       *thisRow = (_List*)theStrings(i1);
        _SimpleList *rowWidths = new _SimpleList;

        checkPointer (rowWidths);

        for (i2=0; i2<thisRow->lLength; i2++) {
            w = GetVisibleStringWidth (*(_String*)((*thisRow)(i2)),f);

            if (columnWidths.lLength<=i2) {
                columnWidths << w;
            } else if (columnWidths.lData[i2]<w) {
                columnWidths.lData[i2] = w;
            }

            (*rowWidths) << w;
        }
        stringWidths << rowWidths;
        DeleteObject (rowWidths);
    }

    thePlacements.Clear();

    w = f.size+margins.top;

    for (i1=0; i1<theStrings.lLength; i1++) {
        _SimpleList *thisRow       = (_SimpleList*)stringWidths(i1),
                     *rowPlacements = new _SimpleList;

        checkPointer (rowPlacements);

        long         cumW = margins.left;

        for (i2=0; i2<thisRow->lLength; i2++) {
            (*rowPlacements) << w;
            (*rowPlacements) << cumW;
            (*rowPlacements) << columnWidths.lData[i2]-thisRow->lData[i2];
            cumW += columnWidths.lData[i2]+margins.right+margins.left;
        }

        thePlacements << rowPlacements;
        w += f.size+margins.top+margins.bottom;
        DeleteObject (rowPlacements);
    }

    w = 0;
    for (i1=0; i1<columnWidths.lLength; i1++) {
        w += columnWidths.lData[i1];
        w += margins.left;
        w += margins.right;
    }

    return w;
}

//__________________________________________________________________
long        ComputeHashMarkPlacement  (_HYRect line, _Parameter min, _Parameter max, _Parameter& tickStep,
                                       _Parameter&minLabel, _SimpleList& offsets, _List& labels,
                                       int count, _HYFont& f)
{
    if (max>min) {
        offsets.Clear();
        labels.Clear();

        long            w = line.right-line.left,
                        h = line.bottom-line.top,
                        t = 0;
        _Parameter      scalingFactor = (w>h)?w/(max-min):h/(max-min);


        tickStep        = log(max-min)/log(10.);
        tickStep        = pow (10.,floor (tickStep))*2.;
        minLabel        = min;

        t = count-1;
        while (t<count) {
            tickStep/=2.;
            if (tickStep>1.) {
                tickStep = floor (tickStep);
            }
            if (min) {
                minLabel = ceil(min/tickStep)*tickStep;
            }
            t = ceil((max-minLabel)/tickStep);
            if (min) {
                minLabel = ceil(min/tickStep)*tickStep;
            }
        }

        if (h>w) {
            _Parameter tracer = minLabel;
            long res  = 0,
                 last = -100000,
                 k,
                 sw;

            w = h;

            for (h=0; h<=t; h++,tracer+=tickStep) {
                _String buffer (tracer);
                sw = GetVisibleStringWidth (buffer,f);
                k  = (tracer-min)*scalingFactor;
                if (k>w) {
                    break;
                }
                if (k-last>f.size) {
                    offsets << k;
                    offsets << sw;
                    labels  && & buffer;
                    last = k;
                }
                if (sw>res) {
                    res = sw;
                }
            }
            return res;
        } else {
            _Parameter tracer = minLabel;
            long last = -100000,
                 k,
                 sw;

            for (h=0; h<=t; h++,tracer+=tickStep) {
                _String buffer (tracer);
                sw = GetVisibleStringWidth (buffer,f);
                k  = (tracer-min)*scalingFactor;
                if (k>w) {
                    break;
                }
                if (k-last>sw) {
                    offsets << k;
                    offsets << sw;
                    labels  && & buffer;
                    last = k;
                }
            }
        }

        return f.size;

    }
    return 0;

}

//__________________________________________________________________
void        _HYGraphicPane::DrawInfoBox (_HYRect theBox, _String boxName)
{
    long   stringWidth,
           avWidth = theBox.right-theBox.left-20;

    theBox.width = 1;

    DrawRect (theBox);
    if (boxName.sLength) {
        stringWidth = GetVisibleStringWidth (boxName, font);
        while (stringWidth>avWidth && boxName.sLength) {
            boxName.Trim (0, boxName.sLength-2);
            stringWidth = GetVisibleStringWidth (boxName, font);
        }
        // measure string length
        _HYRect     textRect = theBox;
        textRect.left += 8;
        textRect.right = textRect.left+stringWidth+4;
        textRect.top-=6;
        textRect.bottom = textRect.top+9;
        EraseRect   (textRect);
        DisplayText (boxName, textRect.bottom, textRect.left+2, true);
    }
}


//EOF