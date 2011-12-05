/*
    A general composite window component object, Windows specifics

    Sergei L. Kosakovsky Pond, June 2000.
*/

#include "HYComponent.h"
#include "HYPlatformComponent.h"
#include "HYPlatformGraphicPane.h"
#include "HYWindow.h"
#include "HYEventTypes.h"
#include "HYCanvas.h"
#include "HYUtils.h"

#include "Windows.h"

extern  long           smallScrollStep;
extern  _HYGuiObject * scrollingWindow;
extern  bool           hScrollingAction;

bool    forceUpdateForScrolling = false;


//__________________________________________________________________________________

void            AlignRectangle (_HYRect& rel , RECT& target , unsigned char alFlags)
{
    long     temp;
    if (alFlags&HY_ALIGN_RIGHT) {
        temp = target.right-target.left;
        target.right = rel.right;
        target.left = target.right - temp;
    } else if (!(alFlags&HY_ALIGN_LEFT)) {
        temp = (rel.right-rel.left-target.right+target.left)/2;
        target.left+=temp;
        target.right+=temp;
    }

    if (alFlags&HY_ALIGN_BOTTOM) {
        temp = target.bottom-target.top;
        target.bottom = rel.bottom;
        target.top = target.bottom - temp;
    } else if (!(alFlags&HY_ALIGN_TOP)) {
        temp = (rel.bottom-rel.top-target.bottom+target.top)/2;
        target.top+=temp;
        target.bottom+=temp;
    }
}

//__________________________________________________________________

_HYPlatformComponent::_HYPlatformComponent(void)
{
    vScroll = hScroll = nil;
    scrollState = 0;
}

//__________________________________________________________________

_HYPlatformComponent::_HYPlatformComponent(_HYRect s,Ptr w)
{
    bool    memError = false;
    vScroll = hScroll = nil;
    parentWindow = (HWND)w;
    if (s.width&HY_COMPONENT_H_SCROLL) {
        hScroll = CreateWindow ("SCROLLBAR","",SBS_HORZ|WS_CHILD|WS_VISIBLE,0,0,100,HY_SCROLLER_WIDTH,parentWindow,NULL,GetModuleHandle (NULL), NULL);
        if (!hScroll) {
            memError = true;
        } else {
            SetScrollRange (hScroll,SB_CTL,0,MAX_CONTROL_VALUE,FALSE);
        }
    }
    if (s.width&HY_COMPONENT_V_SCROLL) {
        vScroll = CreateWindow ("SCROLLBAR","",SBS_VERT|WS_CHILD|WS_VISIBLE,0,0,HY_SCROLLER_WIDTH,100,parentWindow,NULL,GetModuleHandle (NULL), NULL);
        if (!vScroll) {
            memError = true;
        } else {
            SetScrollRange (vScroll,SB_CTL,0,MAX_CONTROL_VALUE,FALSE);
        }
    }
    if (memError) {
        _String errMsg = "Could not allocate memory for a window component structure.";
        FlagError (errMsg);
    }
}

//__________________________________________________________________

void _HYPlatformComponent::_CleanUp(void)
{
}

//__________________________________________________________________
long        _HYPlatformComponent::_GetHScrollerPos (void)
{
    if (hScroll) {
        SCROLLINFO  sI;
        sI.cbSize = sizeof (SCROLLINFO);
        sI.fMask  = SIF_POS|SIF_PAGE;

        GetScrollInfo (hScroll, SB_CTL, &sI);

        return        (sI.nPos*(_Parameter)MAX_CONTROL_VALUE)/((_Parameter)MAX_CONTROL_VALUE-(_Parameter)sI.nPage);
    } else {
        return 0;
    }
}
//__________________________________________________________________
long        _HYPlatformComponent::_GetVScrollerPos (void)
{
    if (vScroll) {
        SCROLLINFO  sI;
        sI.cbSize = sizeof (SCROLLINFO);
        sI.fMask  = SIF_POS|SIF_PAGE;

        GetScrollInfo (vScroll, SB_CTL, &sI);

        return (sI.nPos*(_Parameter)MAX_CONTROL_VALUE)/((_Parameter)MAX_CONTROL_VALUE-sI.nPage);
        //printf ("%d %d %d\n", sI.nPos, sI.nPage, (long)res);

    } else {
        return 0;
    }
}

//__________________________________________________________________
void        _HYPlatformComponent::_SetHScrollerPos (long nv)
{
    if (hScroll&&(scrollState & HY_COMPONENT_ENABLE_H_SCROLL)) {
        if (nv<0) {
            nv = 0;
        } else if (nv>MAX_CONTROL_VALUE) {
            nv = MAX_CONTROL_VALUE;
        }

        SCROLLINFO  sI;
        sI.cbSize = sizeof (SCROLLINFO);
        sI.fMask  = SIF_POS|SIF_PAGE;

        GetScrollInfo (hScroll, SB_CTL, &sI);

        sI.nPos = (nv*((_Parameter)MAX_CONTROL_VALUE-(_Parameter)sI.nPage))/(_Parameter)MAX_CONTROL_VALUE;

        //printf (">%d %d %d\n", nv, sI.nPos,sI.nPage);

        sI.fMask  = SIF_POS;
        SetScrollInfo (hScroll, SB_CTL, &sI, true);

    }
}

//__________________________________________________________________
void        _HYPlatformComponent::_SetVScrollerPos (long nv)
{
    if (vScroll&&(scrollState & HY_COMPONENT_ENABLE_V_SCROLL)) {
        if (nv<0) {
            nv = 0;
        } else if (nv>MAX_CONTROL_VALUE) {
            nv = MAX_CONTROL_VALUE;
        }

        SCROLLINFO  sI;
        sI.cbSize = sizeof (SCROLLINFO);
        sI.fMask  = SIF_POS|SIF_PAGE;

        GetScrollInfo (vScroll, SB_CTL, &sI);

        sI.nPos = (nv*((_Parameter)MAX_CONTROL_VALUE-(_Parameter)sI.nPage))/(_Parameter)MAX_CONTROL_VALUE;

        sI.fMask  = SIF_POS;
        SetScrollInfo (vScroll, SB_CTL, &sI, true);
        //printf ("%d %d\n", nv, sI.nPos);

    }
}

//__________________________________________________________________
void _HYPlatformComponent::Duplicate (BaseRef s)
{
    _HYComponent* theS = (_HYComponent*)s;
    _CleanUp();
    hScroll = theS->hScroll;
    vScroll = theS->vScroll;
}

//__________________________________________________________________
void _HYPlatformComponent::_SetDimensions (_HYRect d,_HYRect r)
{
    _SetVisibleSize (r);
}

//__________________________________________________________________
void _HYPlatformComponent::_MarkForUpdate (void)
{
    RECT inv;
    inv.left = rel.left;
    inv.right = rel.right;
    inv.bottom = rel.bottom;
    inv.top = rel.top;

    if (forceUpdateForScrolling) {
        _HYRect invR = {inv.top, inv.left, inv.bottom, inv.right, (long)GetDC (parentWindow)};
        ((_HYComponent*)this)->Paint ((Ptr)&invR);
        ReleaseDC (parentWindow, (HDC)invR.width);
    } else {
        InvalidateRect (parentWindow,&inv, ((_HYComponent*)this)->settings.width&HY_COMPONENT_TRANSP_BG);
    }
}

//__________________________________________________________________

void _HYPlatformComponent::_MarkContentsForUpdate (void)
{
    RECT inv;
    inv.left = rel.left;
    inv.right = rel.right;
    inv.bottom = rel.bottom;
    inv.top = rel.top;
    if (hScroll) {
        inv.bottom -= HY_SCROLLER_WIDTH;
    }
    if (vScroll) {
        inv.right -= HY_SCROLLER_WIDTH;
    }
    InvalidateRect (parentWindow,&inv, false);
}

//__________________________________________________________________

void        _HYPlatformComponent::_SetVisibleSize (_HYRect r)
{
    _HYComponent * theParent = (_HYComponent*)this;
    SCROLLINFO     si;
    long           t,v;

    _Parameter     newSize;

    scrollState = 0;

    si.cbSize = sizeof (SCROLLINFO);

    if (hScroll&&!vScroll) { // only horizontal scroll bar
        SetWindowPos (hScroll,NULL,r.left,r.bottom-HY_SCROLLER_WIDTH,r.right-r.left+1,HY_SCROLLER_WIDTH,SWP_NOZORDER);

        t = theParent->GetMaxW();
        v = r.right-r.left+1;

        if (t>v) {
            si.fMask = SIF_ALL;
            GetScrollInfo (hScroll, SB_CTL, & si);

            newSize = MAX_CONTROL_VALUE*(_Parameter)v/t;
            si.nPage = newSize;

            SetScrollInfo  (hScroll,SB_CTL,&si,true);
            EnableScrollBar(hScroll,SB_CTL,ESB_ENABLE_BOTH);
            ShowScrollBar  (hScroll,SB_CTL,true);
            scrollState |= HY_COMPONENT_ENABLE_H_SCROLL;

        } else {
            EnableScrollBar(hScroll,SB_CTL,ESB_DISABLE_BOTH);
            ShowScrollBar  (hScroll,SB_CTL,false);
        }

    }
    if (vScroll&&!hScroll) { // only vertical scroll bar
        SetWindowPos (vScroll,NULL,r.right-HY_SCROLLER_WIDTH,r.top,HY_SCROLLER_WIDTH,r.bottom-r.top+1,SWP_NOZORDER);

        t = theParent->GetMaxH();
        v = r.bottom-r.top+1;

        if (t>v) {
            si.fMask = SIF_ALL;
            GetScrollInfo (vScroll, SB_CTL, & si);

            newSize = MAX_CONTROL_VALUE*(_Parameter)v/t;
            si.nPage = newSize;

            SetScrollInfo  (vScroll,SB_CTL,&si,false);
            EnableScrollBar(vScroll,SB_CTL,ESB_ENABLE_BOTH);
            ShowScrollBar  (hScroll,SB_CTL,true);
            scrollState |= HY_COMPONENT_ENABLE_V_SCROLL;
        } else {
            EnableScrollBar(vScroll,SB_CTL,ESB_DISABLE_BOTH);
            ShowScrollBar  (hScroll,SB_CTL,false);
        }
    }
    if (vScroll&&hScroll) {
        SetWindowPos (hScroll,NULL,r.left,r.bottom-HY_SCROLLER_WIDTH,r.right-r.left-HY_SCROLLER_WIDTH,HY_SCROLLER_WIDTH,SWP_NOZORDER|SWP_NOACTIVATE);

        t = theParent->GetMaxW();
        v = r.right-r.left+1+HY_SCROLLER_WIDTH;

        if (t>v) {
            si.fMask = SIF_ALL;
            GetScrollInfo (hScroll, SB_CTL, & si);

            newSize = MAX_CONTROL_VALUE*(_Parameter)v/t;
            si.nPage = newSize;

            SetScrollInfo  (hScroll,SB_CTL,&si,true);
            EnableScrollBar(hScroll,SB_CTL,ESB_ENABLE_BOTH);
            ShowScrollBar  (hScroll,SB_CTL,true);
            scrollState |= HY_COMPONENT_ENABLE_H_SCROLL;
        } else {

            EnableScrollBar(hScroll,SB_CTL,ESB_DISABLE_BOTH);
            ShowScrollBar  (hScroll,SB_CTL,false);
        }

        SetWindowPos (vScroll,NULL,r.right-HY_SCROLLER_WIDTH,r.top,HY_SCROLLER_WIDTH,r.bottom-r.top+1,SWP_NOZORDER);

        t = theParent->GetMaxH();
        v = r.bottom-r.top+1+HY_SCROLLER_WIDTH;

        if (t>v) {
            si.fMask = SIF_ALL;
            GetScrollInfo (vScroll, SB_CTL, & si);

            newSize = MAX_CONTROL_VALUE*(_Parameter)v/t;
            si.nPage = newSize;

            SetScrollInfo  (vScroll,SB_CTL,&si,false);
            EnableScrollBar(vScroll,SB_CTL,ESB_ENABLE_BOTH);
            ShowScrollBar  (hScroll,SB_CTL,true);
            scrollState |= HY_COMPONENT_ENABLE_V_SCROLL;
        } else {
            EnableScrollBar(vScroll,SB_CTL,ESB_DISABLE_BOTH);
            ShowScrollBar  (vScroll,SB_CTL,false);
        }
    }
    rel = r;
}
//__________________________________________________________________

void        _HYPlatformComponent::_Paint (Ptr p)
{
    if (hScroll) {
        ShowScrollBar(hScroll,SB_CTL,TRUE);
    }
    if (vScroll) {
        ShowScrollBar(vScroll,SB_CTL,TRUE);
    }
    _HYComponent* parent = (_HYComponent*)this;

    if (parent->settings.width&HY_COMPONENT_BORDER) {
        HDC      theDC = (HDC)((_HYRect*)p)->width;

        HPEN     blackPen = CreatePen (PS_SOLID, 1, RGB (0,0,0)),
                 oldPen   = (HPEN)SelectObject (theDC, blackPen);

        MoveToEx (theDC,rel.left,rel.top,nil);

        if (parent->settings.width&HY_COMPONENT_BORDER_T) {
            LineTo (theDC,rel.right,rel.top);
        }

        MoveToEx (theDC,rel.right-1,rel.top,nil);

        if (parent->settings.width&HY_COMPONENT_BORDER_R) {
            LineTo (theDC,rel.right-1,rel.bottom);
        }

        MoveToEx (theDC,rel.right,rel.bottom-1,nil);

        if (parent->settings.width&HY_COMPONENT_BORDER_B) {
            LineTo (theDC,rel.left-1,rel.bottom-1);
        }

        MoveToEx (theDC,rel.left,rel.bottom,nil);

        if (parent->settings.width&HY_COMPONENT_BORDER_L) {
            LineTo (theDC,rel.left,rel.top-1);
        }

        SelectObject (theDC, oldPen);
        DeleteObject (blackPen);

        /*UINT  bFlags = 0;

        if (parent->settings.width&HY_COMPONENT_BORDER_T)
            bFlags |= BF_TOP;

        if (parent->settings.width&HY_COMPONENT_BORDER_R)
            bFlags |= BF_RIGHT;

        if (parent->settings.width&HY_COMPONENT_BORDER_B)
            bFlags |= BF_BOTTOM;

        if (parent->settings.width&HY_COMPONENT_BORDER_L)
            bFlags |= BF_LEFT;

        RECT     rr    =  HYRect2Rect (rel);
        DrawEdge (theDC,&rr,EDGE_RAISED, bFlags);*/

    }

    if (parent->settings.width&HY_COMPONENT_WELL) {
        HDC      theDC = (HDC)((_HYRect*)p)->width;

        RECT     rr    =  HYRect2Rect (rel);
        InflateRect (&rr,-2,-2);
        DrawEdge (theDC,&rr,EDGE_SUNKEN, BF_RECT);

    }
}

//__________________________________________________________________

void      _HYPlatformComponent::_Update (Ptr p)
{
    _Paint (p);
}

//__________________________________________________________________

void        _HYPlatformComponent::_Activate (void)
{
    _HYComponent*   theParent = (_HYComponent*)this;
    if (hScroll) {
        if (theParent->GetMaxW()>theParent->rel.right-theParent->rel.left+1+HY_SCROLLER_WIDTH) {
            scrollState |= HY_COMPONENT_ENABLE_H_SCROLL;
            EnableScrollBar(hScroll,SB_CTL,ESB_ENABLE_BOTH);
        }
    }
    if (vScroll) {
        if (theParent->GetMaxH()>theParent->rel.bottom-theParent->rel.top+1+HY_SCROLLER_WIDTH) {
            scrollState |= HY_COMPONENT_ENABLE_V_SCROLL;
            EnableScrollBar(vScroll,SB_CTL,ESB_ENABLE_BOTH);
        }
    }
    activationFlag = true;
}

//__________________________________________________________________

void        _HYPlatformComponent::_Deactivate (void)
{
    if (hScroll) {
        EnableScrollBar(hScroll,SB_CTL,ESB_DISABLE_BOTH);
    }
    if (vScroll) {
        EnableScrollBar(vScroll,SB_CTL,ESB_DISABLE_BOTH);
    }

    scrollState = 0;

    activationFlag = false;
}

//__________________________________________________________________

bool _HYPlatformComponent::_ProcessOSEvent (Ptr vEvent)
{
    _HYWindowsUIMessage * theEvent = (_HYWindowsUIMessage *)vEvent;
    switch (theEvent->iMsg) {
    case WM_HSCROLL:
    case WM_VSCROLL: {
        if((((HWND)(theEvent->lParam))==hScroll)||(((HWND)(theEvent->lParam))==vScroll)) {
            _HYComponent* theParent = (_HYComponent*)this;

            long       invisPixels,
                       currentValue,
                       oldValue;

            if (theEvent->iMsg == WM_HSCROLL) {
                invisPixels = theParent->GetMaxW()-(theParent->GetHSize());
                currentValue = _GetHScrollerPos();
            } else {
                invisPixels = theParent->GetMaxH()-(theParent->GetVSize());
                currentValue = _GetVScrollerPos();
            }


            oldValue = currentValue;
            invisPixels = (_Parameter)MAX_CONTROL_VALUE/invisPixels;
            if (!invisPixels) {
                invisPixels = 1;
            }

            switch (LOWORD(theEvent->wParam)) {
            case SB_LINEUP:
                currentValue-=invisPixels;
                break;
            case SB_LINEDOWN:
                currentValue+=invisPixels;
                break;
            case SB_PAGEUP:
                currentValue-=100*invisPixels;
                break;
            case SB_PAGEDOWN:
                currentValue+=100*invisPixels;
                break;
            case SB_THUMBPOSITION:
            case SB_THUMBTRACK: {

                SCROLLINFO  sI;
                sI.cbSize = sizeof (SCROLLINFO);
                sI.fMask  = SIF_PAGE|SIF_TRACKPOS;
                if (theEvent->iMsg == WM_HSCROLL) {
                    GetScrollInfo (hScroll, SB_CTL, &sI);
                } else {
                    GetScrollInfo (vScroll, SB_CTL, &sI);
                }

                currentValue =  sI.nTrackPos;

                currentValue = (currentValue*(_Parameter)MAX_CONTROL_VALUE)/(MAX_CONTROL_VALUE-(_Parameter)sI.nPage);
                break;
            }
            case SB_TOP:
                currentValue = 0;
                break;
            case  SB_BOTTOM:
                currentValue = MAX_CONTROL_VALUE;
                break;
            default:
                return true;
            }

            if (currentValue<0) {
                currentValue = 0;
            }
            if (currentValue>MAX_CONTROL_VALUE) {
                currentValue = MAX_CONTROL_VALUE;
            }

            if (currentValue!=oldValue) {
                if (theEvent->iMsg == WM_HSCROLL) {
                    _SetHScrollerPos (currentValue);
                    theParent->ProcessEvent (generateScrollEvent (currentValue-oldValue,0));
                } else {
                    _SetVScrollerPos (currentValue);
                    theParent->ProcessEvent (generateScrollEvent (0,currentValue-oldValue));
                }
            }
            return true;
        }
        break;
    }
    }
    return false;
}

//__________________________________________________________________


_HYRect _HYPlatformComponent::_VisibleContents (Ptr p)
{
    _HYComponent* theParent = (_HYComponent*)this;
    _HYRect     * r = (_HYRect*)p,
                  res;

    short v;
    short windowW=r->right-r->left,
          windowH=r->bottom-r->top;
    if ((!hScroll)&&(!vScroll)) {
        res.left = theParent->hOrigin;
        res.right = res.left+windowW;
        res.top = theParent->vOrigin;
        res.bottom = res.top+windowH;
        return res;
    }
    if (hScroll) {
        windowH-=HY_SCROLLER_WIDTH;
    }
    if (vScroll) {
        windowW-=HY_SCROLLER_WIDTH;
    }
    if (hScroll) {
        if (windowW>theParent->GetMaxW()) {
            res.left = 0;
            res.right = theParent->GetMaxW();
        } else {
            v = GetScrollPos (hScroll,SB_CTL);
            res.left = (theParent->GetMaxW()-windowW)*v/(_Parameter)MAX_CONTROL_VALUE;
            res.right = res.left+windowW;
        }
    } else {
        res.left = 0;
        res.right = windowW;
    }
    if (vScroll) {
        if (windowH>theParent->GetMaxH()) {
            res.top = 0;
            res.bottom = theParent->GetMaxH();
        } else {
            v = GetScrollPos (vScroll,SB_CTL);
            res.top = (theParent->GetMaxH()-windowH)*v/(_Parameter)MAX_CONTROL_VALUE;
            res.bottom = res.top+windowH;
        }
    } else {
        res.top = 0;
        res.bottom = windowH;
    }
    return res;
}

//__________________________________________________________________

void    _HYCanvas::_Paint (Ptr p)
{
    _HYPlatformComponent::_Paint(p);
    _HYRect * destR = (_HYRect*)p;
    HDC pDC = (HDC)destR->width;
    RECT srcRect,destRect;
    destRect.right = destR->right;
    destRect.bottom = destR->bottom;
    if (HasHScroll()) {
        destRect.bottom-=HY_SCROLLER_WIDTH;
    }
    if (HasVScroll()) {
        destRect.right-=HY_SCROLLER_WIDTH;
    }
    destRect.left = destR->left;
    destRect.top = destR->top;
    _HYRect srcR = _VisibleContents (p);
    srcRect.right = srcR.right;
    srcRect.left = srcR.left;
    srcRect.top = srcR.top;
    srcRect.bottom = srcR.bottom;
//  printf ("\n%d %d %d %d %d %d %d %d\n",srcRect.top, srcRect.left, srcRect.bottom,srcRect.right,
//                                        destRect.top, destRect.left, destRect.bottom,destRect.right);
    BitBlt (pDC,destRect.left,destRect.top,srcRect.right-srcRect.left,srcRect.bottom-srcRect.top,
            thePane,srcRect.left,srcRect.top,SRCCOPY);
}

//__________________________________________________________________

void    _HYCanvas::_Update (Ptr p)
{
    _Paint(p);
}

//__________________________________________________________________

bool _HYCanvas::_ProcessOSEvent (Ptr vEvent)
{
    _HYWindowsUIMessage*    theEvent = (_HYWindowsUIMessage*)vEvent;

    switch (theEvent->iMsg) {
    case WM_RBUTTONDOWN: {
        _String s ("Save as picture");
        _List   l;
        l && & s;

        POINT loc = {(short)LOWORD (theEvent->lParam),(short)HIWORD (theEvent->lParam)};
        ClientToScreen (parentWindow, &loc);

        _String pulldownValue = HandlePullDown (l,loc.x,loc.y,0);
        switch (l.Find(&pulldownValue)) {
        case 0:
            s = "Save picture";
            _SavePicture (s);
            return true;
        }
        break;
    }
    case WM_LBUTTONDOWN: {
        if ((messageRecipient)&&(doMouseClicks)) {
            messageRecipient->ProcessEvent(generateContextPopUpEvent (GetID(),LOWORD(theEvent->lParam)-rel.left,
                                           HIWORD(theEvent->lParam)-rel.top));
            return true;
        }
        break;
    }
    }
    return false;
}


//EOF