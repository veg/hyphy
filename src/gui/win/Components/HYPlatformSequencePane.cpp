/*
    Sequence Panel for Win32 API

    Sergei L. Kosakovsky Pond, May 2000-January 2003
*/

#include "errorfns.h"
#include "HYSequencePanel.h"
#include "HYUtils.h"
#include "HYEventTypes.h"
#include "HYPlatformWindow.h"
#include "HYTableWindow.h"


//__________________________________________________________________

void    _HYSequencePane::_Paint (Ptr p)
{
    long        saveBorder = settings.width & HY_COMPONENT_BORDER;
    settings.width -= saveBorder;
    _HYPlatformComponent::_Paint(p);
    settings.width += saveBorder;
    _HYRect*    destR = (_HYRect*)p;

    RECT        srcRect,
                destRect;

    HDC pDC = (HDC)destR->width;

    destRect.right      = destR->right;
    destRect.bottom     = destR->bottom;
    if (HasHScroll()) {
        destRect.bottom-= HY_SCROLLER_WIDTH;
    }

    if (HasVScroll()) {
        destRect.right -= HY_SCROLLER_WIDTH;
    }

    destRect.left       = destR->left;
    destRect.top        = destR->top;
    _HYRect srcR        = _VisibleContents (p);
    srcRect.right       = srcR.right-srcR.left;
    srcRect.left        = 0;
    srcRect.bottom      = srcR.bottom-srcR.top;
    srcRect.top         = 0;

    BitBlt (pDC,destRect.left,destRect.top,srcRect.right-srcRect.left,srcRect.bottom-srcRect.top,
            thePane,srcRect.left,srcRect.top,SRCCOPY);
}

//__________________________________________________________________

bool _HYSequencePane::_ProcessOSEvent (Ptr vEvent)
{
    static  bool    amScrolling = false,
                    vertical;

    static  POINT   localPt;

    static  long    originalStart,
            originalSpan,
            lastClick,
            firstClick;

    if (_HYPlatformComponent::_ProcessOSEvent (vEvent)) {
        return true;
    }
    if (!active) {
        return false;
    }

    short lastH, lastV;
    POINT globalPt;

    _HYWindowsUIMessage*    theEvent = (_HYWindowsUIMessage*)vEvent;

    switch (theEvent->iMsg) {
    case WM_RBUTTONDOWN:
    case WM_LBUTTONDOWN:
    case WM_LBUTTONDBLCLK:
        {

            lastH = (short)LOWORD (theEvent->lParam),
            lastV = (short)HIWORD (theEvent->lParam);

            globalPt = (POINT) {
                lastH, lastV
            };

            localPt  = (POINT) {
                lastH - rel.left , lastV - rel.top
            };

            vertical = (localPt.x<headerWidth)&&(localPt.y>=(GetSlotHeight()+1));


            if ((theEvent->iMsg == WM_LBUTTONDOWN)||(theEvent->iMsg == WM_LBUTTONDBLCLK)) {
                forceUpdateForScrolling = true;
                if (vertical)
                    ProcessVSelectionChange (localPt.x,localPt.y,GetAsyncKeyState (VK_SHIFT) & 0x8000,
                                             GetAsyncKeyState (VK_CONTROL) & 0x8000, false, theEvent->iMsg == WM_LBUTTONDBLCLK);
                else
                    ProcessSelectionChange  (localPt.x,localPt.y,GetAsyncKeyState (VK_SHIFT) & 0x8000,
                                             GetAsyncKeyState (VK_CONTROL) & 0x8000);
                forceUpdateForScrolling = false;

                ClientToScreen (parentWindow, &globalPt);

                if (DragDetect (parentWindow, globalPt)) {
                    if (messageRecipient) {
                        SetCapture (parentWindow);
                        ((_HYTWindow*)messageRecipient)->trackMouseComponent = this;
                    }
                    amScrolling = true;
                    if (vertical) {
                        originalStart = startRow,
                        originalSpan  = endRow-startRow;
                        lastClick = -2;
                        firstClick = (localPt.y-(GetSlotHeight()+1))/GetSlotHeight();
                    }
                }

                return true;
            }

            if ((theEvent->iMsg == WM_RBUTTONDOWN)&&(vertical&&vselection.lLength)||((!vertical)&&selection.lLength)) {
                ClientToScreen (parentWindow, &globalPt);
                ProcessContextualPopUp (globalPt.x, globalPt.y);
                return true;
            }
        }
        break;

    case WM_LBUTTONUP:
        if (amScrolling) {
            amScrolling = false;
            if (messageRecipient) {
                ReleaseCapture ();
                ((_HYTWindow*)messageRecipient)->trackMouseComponent = nil;
            }
            if  (vertical) {
                RECT invalRect = {rel.left,rel.top+(GetSlotHeight()+1)+1,rel.left+headerWidth,rel.bottom-HY_SCROLLER_WIDTH};
                InvalidateRect (parentWindow,&invalRect,false);
                if ((localPt.x<headerWidth)&&(localPt.x>0)&&(lastClick>-2)) {
                    MoveSpecies (firstClick+originalStart,lastClick+startRow);
                }
            }
        }
        return true;
        break;

    case WM_MOUSEMOVE:
        if ((theEvent->wParam & MK_LBUTTON)&&(amScrolling)) {
            POINT mousePt = {((short)LOWORD (theEvent->lParam))-rel.left,
                             ((short)HIWORD (theEvent->lParam))-rel.top
                            };
            if (vertical) {

                long  wHeight = rel.bottom-rel.top-HY_SCROLLER_WIDTH,
                      slotHeight = GetSlotHeight();


                forceUpdateForScrolling = true;
                if ((mousePt.y<(GetSlotHeight()+1))||(localPt.y!=mousePt.y)||(mousePt.y>wHeight)) {
                    localPt = mousePt;
                    if (mousePt.y>wHeight) {
                        // scroll down
                        if ((endRow<=speciesIndex.lLength)&&(vselection.lData[0]!=speciesIndex.lLength-1)) {
                            if (endRow-startRow<originalSpan) {
                                break;
                            }
                            startRow++;
                            endRow++;
                            _SetVScrollerPos(((double)MAX_CONTROL_VALUE*startRow)/
                                             (speciesIndex.lLength-endRow+startRow+1));
                            BuildPane();
                            _MarkForUpdate();
                            lastClick = -2;
                        }
                        break;
                    } else {
                        mousePt.y-=(GetSlotHeight()+1);
                        if (mousePt.y<=slotHeight) {
                            if (mousePt.y>=0) {
                                if (mousePt.y<slotHeight/2) {
                                    mousePt.y = -1;
                                } else {
                                    mousePt.y = 0;
                                }
                            } else {
                                // scroll up
                                if (startRow>0) {
                                    startRow--;
                                    endRow--;
                                    _SetVScrollerPos(((double)MAX_CONTROL_VALUE*startRow)/(speciesIndex.lLength-endRow+startRow+1));
                                    BuildPane();
                                    _MarkForUpdate();
                                    lastClick = -2;
                                }
                                break;
                            }
                        } else {
                            mousePt.y=(mousePt.y-(GetSlotHeight()+1))/slotHeight;
                        }
                    }

                    if ((mousePt.y<-1)||(mousePt.y>=(endRow-startRow))) {
                        break;
                    }
                    if (mousePt.y!=lastClick) {
                        HDC winDC    = GetDC   (parentWindow);
                        int saveROP2 = GetROP2 (winDC);
                        SetROP2 (winDC,R2_NOT);
                        if (lastClick>=-1) {
                            lastClick = (GetSlotHeight()+1)+slotHeight*(lastClick+1)+rel.top+1;
                            MoveToEx(winDC,rel.left+1,lastClick,nil);
                            LineTo  (winDC,rel.left+headerWidth-1,lastClick);
                        }
                        lastClick = mousePt.y;
                        if (lastClick+startRow!=firstClick+originalStart) {
                            mousePt.y = (GetSlotHeight()+1)+slotHeight*(lastClick+1)+rel.top+1;
                            MoveToEx (winDC,rel.left+1,mousePt.y,nil);
                            LineTo   (winDC,rel.left+headerWidth-1,mousePt.y);
                        }
                        SetROP2 (winDC,saveROP2);
                        ReleaseDC (parentWindow,winDC);
                    }
                }
                forceUpdateForScrolling = false;
                return true;
            } else {
                if (((mousePt.x<headerWidth)&&(startColumn>0))||(localPt.x!=mousePt.x)||(mousePt.x>_HYCanvas::GetMaxW()-5)) {
                    forceUpdateForScrolling = true;
                    ProcessSelectionChange (mousePt.x,mousePt.y,true,true,true);
                    forceUpdateForScrolling = false;
                    localPt = mousePt;
                }
                return true;
            }
        }
        break;
    }

    return false;
}


//EOF