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

#include "HYWindow.h"
#include "HYGWindow.h"
#include "HYTableWindow.h"
#include "HYEventTypes.h"

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

//__________________________________________________________________

_HYWindow::_HYWindow(unsigned char windowFlag,_String wTitle,bool windowVisibility, Ptr p):
    _HYPlatformWindow (windowFlag,wTitle,windowVisibility, p)
{
    top = left = 0;
    right = 300;
    bottom = 250;
    if (windowFlag&HY_WINDOW_SCROLL) {
        contentWidth = 300;
        contentHeight = 250;
    }
    windowTitle = wTitle;
#ifdef __HYPHY_FLTK__
    _SetTitle (windowTitle);
#endif
    if (windowVisibility) {
        Activate();
        Paint(nil);
    }
    windowObjectRefs << (long)this;
}

//__________________________________________________________________

_HYWindow::~_HYWindow(void)
{
}

//__________________________________________________________________

void _HYWindow::SetTitle(_String windowT)
{
    windowTitle = windowT;
    _SetTitle (windowTitle);
}

//__________________________________________________________________

void _HYWindow::Show(void)
{
    _Show ();
    Update(nil);
}

//__________________________________________________________________

void _HYWindow::Hide(void)
{
    _Hide ();
}

//__________________________________________________________

bool    _HYWindow::ProcessGEvent (_HYEvent* e)
{
    long k;
    if (e->EventClass()==_hyGlobalCloseWindow) {
        k = e->EventCode().toNum();
        if (MatchID (k)) {
            Close (nil);
            return true;
        }
    }
    return false;
}
//__________________________________________________________________

void _HYWindow::Grow(Ptr theData)
{
    unsigned long newSize = _Grow (theData);
    if (newSize) {
        right = left+(newSize&0x0000ffff);
        bottom = top+(((newSize&0xffff0000)>>16)&0xffff);
        SetWindowRectangle (top,left,bottom,right);
        Update(nil);
    }
}

//__________________________________________________________________

void _HYWindow::Move(Ptr theData)
{
    _Move (theData);
    _Update (nil);
}

//__________________________________________________________________

void _HYWindow::SetPosition(int left, int top)
{
    _SetPosition (left,top);
}

//__________________________________________________________________

void _HYWindow::SelectThisWindow(void)
{
    _SelectWindow ();
}

//__________________________________________________________________

void _HYWindow::SetFont(_HYFont&)
{
}

//__________________________________________________________________

bool _HYWindow::Close(Ptr theData)
{
    if (!_Close(theData)) {
        return false;
    }

    //Ptr p = GetOSWindowData();
    DeleteObject (this);
    //_Close2 (p);
    return true;
}

//__________________________________________________________________

void _HYWindow::Activate(void)
{
    _Activate ();
    Update(nil);
}

//__________________________________________________________________

void _HYWindow::Deactivate(void)
{
    _Deactivate ();
}

//__________________________________________________________________

void _HYWindow::SetWindowRectangle(int t, int l, int b, int r, bool ss)
{
    _SetWindowRectangle (t,l,b,r,ss);
    top = t;
    left = l;
    right = r;
    bottom = b;
}

//__________________________________________________________________

void _HYWindow::SetContentSize(int w, int h)
{
    contentWidth = w;
    contentHeight = h;
    _SetContentSize (w,h);
}

//__________________________________________________________________

_HYGWindow::_HYGWindow(_String windowTitle,int ih, int iw, int id, bool vis,bool scale):
    _HYWindow         (scale?(HY_WINDOW_CLOSE|HY_WINDOW_SIZE|HY_WINDOW_FLUSHED):(HY_WINDOW_CLOSE|HY_WINDOW_SIZE|HY_WINDOW_SCROLL),windowTitle,vis),
    _HYGraphicPane    (ih,iw,id)
{
    SetContentSize (iw,ih);
}

//__________________________________________________________________

bool    _HYGWindow::ProcessEvent (_HYEvent* e)
{
    if (e->EventClass() == _hyScrollingEvent) {
        Update (nil);
        DeleteObject (e);
        return true;
    }
    return      _HYWindow::ProcessEvent (e);
}


//__________________________________________________________________

_HYTWindow::_HYTWindow (_String title, char list, bool dlog, Ptr sheetParent):
    _HYWindow ((list==2)?HY_WINDOW_SIZE:(list?(HY_WINDOW_CLOSE|HY_WINDOW_SIZE|HY_WINDOW_ZOOM):(dlog?(sheetParent?HY_WINDOW_SHEET:HY_WINDOW_DLOG):(HY_WINDOW_CLOSE|HY_WINDOW_SIZE|HY_WINDOW_NOLIST|HY_WINDOW_ZOOM))),title,false, sheetParent),
    _HYPlatformTWindow ((Ptr)this)
{
    rows = columns = 0;
    keyboardFocus = -1;
    lastMouseComponent = -1;
}

//__________________________________________________________________

void    _HYTWindow::SetTableDimensions (int r, int c)
{
    if (!cells.lLength) {
        rows = r;
        columns = c;
        for (long i=0; i<r; i++)
            for (long j=0; j<c; j++) {
                cells<<0;
            }
    } else {
        if ((rows!=r)||(columns!=c)) {
            cells.Clear();
            SetTableDimensions(r,c);
        }
    }
}

//__________________________________________________________________

void    _HYTWindow::SetCell (int r, int c, _HYGuiObject* d)
{
    long f = components._SimpleList::Find ((long)d);
    if (f>=0) {
        cells.lData[r*columns+c] = f;
    }
}

//__________________________________________________________________

void    _HYTWindow::AddObject (_HYGuiObject* d, bool dup, long r, long c)
{
    if (r >= 0 && c >= 0) {
        cells.lData[r*columns+c] = components.lLength;
    }

    if (dup) {
        components << d;
    } else {
        components.AppendNewInstance(d);
    }

    componentL<<0;
    componentR<<0;
    componentT<<0;
    componentB<<0;

}

//__________________________________________________________________

void    _HYTWindow::AddKeyboardChainObject (_HYGuiObject* d)
{
    long f = components._SimpleList::Find ((long)d);
    if (f>=0) {
        keyboardFocusChain << f;
    }
}

//__________________________________________________________

bool    _HYTWindow::ProcessEvent (_HYEvent* e)
{
    long i;
    if (e->EventClass()==_hyKeyboardFocusEvent) {
        i = MatchComponentID (e->EventCode());

        if (keyboardFocus>=0) {
            ((_HYComponent*)components(keyboardFocus))->UnfocusComponent();
        }

        if (i<components.lLength) {
            keyboardFocus=i;
            ((_HYComponent*)components(keyboardFocus))->FocusComponent();
        } else {
            keyboardFocus=-1;
        }
    }
    DeleteObject (e);
    return false;
}


//__________________________________________________________________

_HYGuiObject*   _HYTWindow::GetCellObject (int r, int c)
{
    return (_HYGuiObject*)components(cells[r*columns+c]);
}

//__________________________________________________________________

_HYGuiObject*   _HYTWindow::GetObject (int index)
{
    return (_HYGuiObject*)components(index);
}

//__________________________________________________________________

_HYRect         _HYTWindow::MinMaxWindowDimensions (void)
{
    _HYRect     res = {0,0,30000,30000,0};

    long    i,
            j,
            t1,
            t2,
            f,
            g,
            minCont = rows;

    for (i=0; i<rows; i++) {
        t1 = t2 = 0;
        for (j=0; j<columns;) {
            g = cells(i*columns+j);
            do {
                f = g;
                j++;
                if (j==columns) {
                    break;
                }
                g = cells(i*columns+j);
            } while (f==g);
            _HYComponent* thisCell = (_HYComponent*)components(f);
            t1 += thisCell->GetMaxLW();
            t2 += thisCell->GetMinW();
        }
        if (t1<res.right) {
            res.right = t1;
        }
        if (t2>res.left) {
            res.left = t2;
        }
    }

    for (j=0; j<columns; j++) {
        t1 = t2 = 0;
        long thisCont = rows;
        for (i=0; i<rows;) {
            g = cells(i*columns+j);
            do {
                f = g;
                i++;
                if (i==rows) {
                    break;
                }
                g = cells(i*columns+j);
                if (f==g) {
                    thisCont --;
                } else {
                    break;
                }
            } while (1);

            _HYComponent* thisCell = (_HYComponent*)components(f);
            t1 += thisCell->GetMaxLH();
            t2 += thisCell->GetMinH();


        }
        if (t1<res.bottom) {
            res.bottom = t1;
        }
        if (t2>res.top) {
            res.top = t2;
        }

        if (thisCont<minCont) {
            minCont = thisCont;
        }
    }

    if (res.left) {
        res.left--;
    }

    res.right--;

    if (res.top) {
        res.top--;
    }

    res.bottom--;

    if (minCont>3) {
        res.bottom-=minCont-3;
    }

    if (!(flags&HY_WINDOW_DLOG)) {
        res.bottom += HY_SCROLLER_WIDTH-1;
        res.top    += HY_SCROLLER_WIDTH-1;
    }

    contentHeight = res.bottom;
    contentWidth  = res.right;
    return res;
}

//__________________________________________________________________

void    _HYTWindow::RecomputeCellRects (void)
{
    int t,l,b,r,i,j,f,g,c,j2,p,p2,mw;
    bool scroller;
    _HYComponent* thisCell;
    VisibleContents (t,l,b,r);
    //r-=l+HY_SCROLLER_WIDTH;
    r-=l;
    b-=t;

    _Parameter lf=(r-dim.left)/columns; // "extra" width per cell
    _Parameter tf=(b-dim.top)/rows;     // "extra" height per cell

    if (!(flags&HY_WINDOW_DLOG)) {
        b-=HY_SCROLLER_WIDTH;
    }


    _SimpleList     alreadyDone;
    for (i=0; i<rows; i++) {
        c = 0;
        p = -1;
        for (j=0; j<columns;) {
            j2 = j;
            f = cells.lData[i*columns+j];
            thisCell = (_HYComponent*)components.lData[f];
            if (alreadyDone.Find(f)>=0) {
                c = componentR.lData[f];
                j++;
                scroller = false;
                continue;
            }
            do {
                g = f;

                j++;
                if (j==columns) {
                    break;
                }
                f = cells.lData[i*columns+j];
            } while (f==g);
            _HYComponent* thisCell;
            thisCell = (_HYComponent*)components.lData[g];
            componentL.lData[g] = c;
            f = lf*(j-j2)+thisCell->GetMinW()-1;
            mw = thisCell->GetMaxLW()-1;
            if (f>=mw) {
                f=mw;
            } else {
                p=g;
                p2=j-1;
            }
            c+=f;
            componentR.lData[g] = c;
            scroller = true;
        }
        if ((c<r)&&(p>=0)) {
            int diff = r-c;
            componentR.lData[p]+=diff;
            _SimpleList localDone;
            for (j=p2+1; j<columns; j++) {
                f = cells.lData[i*columns+j];
                if (localDone.Find(f)>=0) {
                    continue;
                } else {
                    localDone<<f;
                }
                componentL.lData[f]+=diff;
                componentR.lData[f]+=diff;
            }
        }
    }
    alreadyDone.Clear();
    for (j=0; j<columns; j++) {
        c = 0;
        p=-1;
        _SimpleList stretch, stretchI,stretchW;
        p2 = alreadyDone.lLength;
        for (i=0; i<rows;) {
            j2 = i;
            f = cells.lData[i*columns+j];
            thisCell = (_HYComponent*)components.lData[f];
            if (alreadyDone.Find(f)>=0) {
                c += componentB.lData[f]-componentT.lData[f];
                i++;
                scroller = false;
                continue;
            }
            do {
                g = f;
                i++;
                if (i==rows) {
                    break;
                }
                f = cells.lData[i*columns+j];
            } while (f==g);
            _HYComponent* thisCell = (_HYComponent*)components.lData[g];
            componentT.lData[g] = c;
            f = tf*(i-j2)+thisCell->GetMinH();
            mw = thisCell->GetMaxLH()-1;
            if (f>mw) {
                f=mw;
            } else {
                stretch<<g;
                stretchI<<i-1;
                stretchW<<mw-f;
                //p=g;
                //p2=i-1;
            }
            c+=f;
            componentB.lData[g] = c;
            scroller = true;
            alreadyDone<<g;
            //c++;
        }
        while ((c<b)&&(stretch.lLength)) {
            int diff = b-c;
            if (diff>stretchW.lData[stretch.lLength-1]) {
                diff = stretchW.lData[stretch.lLength-1];
            }

            componentB.lData[stretch.lData[stretch.lLength-1]]+=diff;
            _SimpleList localDone;
            for (i=stretchI.lData[stretch.lLength-1]+1; i<rows; i++) {
                f = cells.lData[i*columns+j];
                p = alreadyDone.Find(f);
                if ((localDone.Find(f)>=0)||(j&&(p>=0)&&(p<p2)))
                    //if (localDone.Find(f)>=0)
                {
                    continue;
                } else {
                    localDone<<f;
                }
                componentT.lData[f]+=diff;
                componentB.lData[f]+=diff;
            }

            c+=diff;
            stretchI.Delete(stretch.lLength-1);
            stretchW.Delete(stretch.lLength-1);
            stretch.Delete(stretch.lLength-1);
        }
    }
}

//__________________________________________________________________

int     _HYTWindow::FindClickedCell (int h, int v)
{
    for (int i=0; i<components.lLength; i++) {
        if (cells.Find (i)>=0)
            if ((h>=componentL.lData[i])&&(h<=componentR.lData[i])
                    &&(v>=componentT.lData[i])&&(v<=componentB.lData[i])) {
                return i;
            }
    }
    return -1;
}

//__________________________________________________________________

void        _HYTWindow::DoMouseWheel (long c, long delta)
{
    if (c>=0) {
        _HYComponent * thisC = (_HYComponent*) components (c);

        if (thisC->HasVScroll()) {
            long    invisPixels     = thisC->GetMaxH()-(thisC->GetVSize()),
                    smallScrollStep = (double)10.*MAX_CONTROL_VALUE/invisPixels,
                    cv               = thisC->_GetVScrollerPos (),
                    cv2;

            if (!smallScrollStep) {
                smallScrollStep = 1;
            }

            cv2 = cv - delta*smallScrollStep;

            if (cv2<0) {
                cv2 = 0;
            } else if (cv2>MAX_CONTROL_VALUE) {
                cv2 = MAX_CONTROL_VALUE;
            }

            thisC->_SetVScrollerPos (cv2);

            thisC->ProcessEvent (generateScrollEvent(0,cv2-cv));
        }
    }
}


//__________________________________________________________________
void        _HYTWindow::Paint (Ptr p)
{
    _Paint(p);
}

//__________________________________________________________________
void        _HYTWindow::Update (Ptr p)
{
    _Update(p);
}

//__________________________________________________________________
void        _HYTWindow::UpdateComponentInfo (void)
{
    SetWindowRectangle(top,left,bottom,right);
}

//__________________________________________________________________
void        _HYTWindow::Activate(void)
{
    _HYWindow::Activate();
}

//__________________________________________________________________
void        _HYTWindow::SetStatusBar(_String& text)
{
    statusBar = text;
    _SetStatusBar (text);
}

//__________________________________________________________________
long        _HYTWindow::MatchComponentID (_String& text)
{
    long k = text.toNum(),
         i;
    for (i=0; i<components.lLength; i++)
        if (((_HYGuiObject*)components(i))->MatchID(k)) {
            return i;
        }
    return -1;
}

//__________________________________________________________________
void        _HYTWindow::HandleCopyPaste (bool paste)
{
    int         componentID = keyboardFocus;

    if          (paste) {
        BaseRef      pastingReference = nil;
        _String*     clipData = _GetPasteString();

        if (clipData->sLength)
            if (componentID < 0) {
                for (componentID=0; componentID<components.lLength; componentID++)
                    if (cells.Find(componentID)>=0)
                        if ((pastingReference =     ((_HYComponent*)components(componentID))->CanPaste (*clipData))) {
                            break;
                        }

            } else {
                pastingReference =  ((_HYComponent*)components(componentID))->CanPaste (*clipData);
            }

        DeleteObject (clipData);
        if (pastingReference) {
            ((_HYComponent*)components(componentID))->HandlePaste (pastingReference);
            DeleteObject (pastingReference);
        }
    } else {
        _String*     clipData = nil;

        if (componentID < 0) {
            for (componentID=0; componentID<components.lLength; componentID++)
                if (cells.Find(componentID)>=0)
                    if (((_HYComponent*)components(componentID))->CanCopy()) {
                        break;
                    }

        } else if (!((_HYComponent*)components(componentID))->CanCopy()) {
            componentID = components.lLength;
        }

        if (componentID < components.lLength) {
            clipData = ((_HYComponent*)components(componentID))->HandleCopy ();
            if (clipData) {
                _SetCopyString (clipData);
                DeleteObject   (clipData);
            }
        }

    }

}

//__________________________________________________________________

void        _HYTWindow::SetWindowRectangle(int t, int l, int b, int r, bool ss)
{
    dim = MinMaxWindowDimensions();
    if (b-t>=dim.bottom) {
        b=t+dim.bottom-1;
    }

    if (b-t<=dim.top) {
        b=t+dim.top-1;
    }

    if (r-l>=dim.right) {
        r=l+dim.right-1;
    }

    if (r-l<=dim.left) {
        r=l+dim.left-1;
    }

    _HYWindow::SetWindowRectangle (t,l,b,r,ss);
    _HYPlatformTWindow::_SetWindowRectangle (t,l,b,r,ss);
    RecomputeCellRects();
    _HYRect  relRect;
    for (int i=0; i<components.lLength; i++) {
        if (cells.Find(i)>=0) {
            relRect.left = componentL.lData[i];
            relRect.right = componentR.lData[i];
            relRect.top = componentT.lData[i];
            relRect.bottom = componentB.lData[i];
            ((_HYComponent*)components(i))->SetVisibleSize(relRect);
        }
    }
}


//__________________________________________________________________
// _HYPWindow
//__________________________________________________________________

_HYPWindow::_HYPWindow(_String windowTitle,int ih, int iw, int id, bool vis):_HYGWindow (windowTitle, ih,iw,id,vis,true)
{
    resizing = false;
    oh       = ih;
    ow       = iw;
}

//__________________________________________________________________

void    _HYPWindow::StartDraw(void)
{
    _HYGraphicPane::StartDraw();
    if (!resizing) {
        _StartPicture();
    }
}

//__________________________________________________________________

void    _HYPWindow::EndDraw(void)
{
    if (!resizing) {
        _EndPicture();
    }

    _HYGraphicPane::EndDraw();
}

//__________________________________________________________________

void    _HYPWindow::SetPaneSize (int h, int w, int d)
{
    resizing = true;
    _HYGraphicPane::SetPaneSize (h,w,d);
    _HYRect        drRect = {0,0,h,w,d};
    _HYGraphicPane::StartDraw();
    _DrawPicture  (drRect);
    _HYGraphicPane::EndDraw();
    resizing = false;
}

//__________________________________________________________________

void    _HYPWindow::Zoom (_Parameter factor)
{
    long newW = contentWidth  * factor,
         newH = contentHeight * factor;

    SetWindowRectangle (0,0,newH, newW);
    SetContentSize     (newW,newH);
}

//__________________________________________________________________

void    _HYPWindow::OriginalSize (void)
{
    SetWindowRectangle (0,0,oh,ow);
    SetContentSize     (ow,oh);
}

//__________________________________________________________________________________

void    HandleGlobalQueueEvent (void)
{
    long mSel,
         menuChoice;

    static  bool handling = false;

    _String evCode (((_HYEvent*)GlobalGUIEventQueue.lData[0])->EventCode());

    mSel = evCode.Find(',')-1;
    if (mSel<-1) {
        mSel = -1;
    }

    mSel = evCode.Cut (0,mSel).toNum();
    menuChoice = -1;
    if (!handling) {
        handling = true;
        for (long f=0; f<windowObjectRefs.countitems(); f++)
            if (((_HYGuiObject*)windowObjectRefs.lData[f])->MatchID(mSel)) {
                menuChoice = f;
            } else {
                ((_HYGuiObject*)windowObjectRefs.lData[f])->ProcessGEvent((_HYEvent*)GlobalGUIEventQueue.lData[0]);
            }

        if (menuChoice>=0) {
            ((_HYGuiObject*)windowObjectRefs.lData[menuChoice])->ProcessGEvent((_HYEvent*)GlobalGUIEventQueue.lData[0]);
        }

        GlobalGUIEventQueue.Delete(0);
        handling = false;
    }
}

//EOF