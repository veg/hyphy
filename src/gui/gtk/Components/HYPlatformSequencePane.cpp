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

#include "errorfns.h"
#include "HYSequencePanel.h"
#include "HYUtils.h"
#include "HYEventTypes.h"
#include "HYPlatformWindow.h"
#include "HYTableWindow.h"
#include <gdk/gdkx.h>


//__________________________________________________________________

void    _HYSequencePane::_Paint (Ptr p)
{
    _HYRect*    destR = (_HYRect*)p;

    _HYRect     srcRect,
                destRect = *destR;

    //printf ("Sequence Paint Called %d %d %d %d\n", destR->left, destR->top, destR->right, destR->bottom);

    if (HasHScroll()) {
        destRect.bottom-= HY_SCROLLER_WIDTH;
    }

    if (HasVScroll()) {
        destRect.right -= HY_SCROLLER_WIDTH;
    }

    //destRect.left         = destR->left;
    //destRect.top      = destR->top;
    _HYRect srcR        = _VisibleContents (p);

    gdk_draw_drawable (GDK_DRAWABLE(parentWindow->window), theContext, thePane, 0, 0,
                       parentWindow->allocation.x+destRect.left, parentWindow->allocation.y+destRect.top, destRect.Width(), destRect.Height());

    long        saveBorder = settings.width & HY_COMPONENT_BORDER;
    settings.width -= saveBorder;
    _HYPlatformComponent::_Paint(p);
    settings.width += saveBorder;
}

//__________________________________________________________________

bool _HYSequencePane::_ProcessOSEvent (Ptr vEvent)
{
    static          bool    amScrolling = false,
                            vertical;

    static          long   localPt_x,
                    localPt_y,
                    originalStart,
                    originalSpan,
                    lastClick,
                    firstClick;


    if (_HYPlatformComponent::_ProcessOSEvent (vEvent)) {
        return true;
    }

    _HY_GTK_UI_Message *theMessage = (_HY_GTK_UI_Message*)vEvent;
    if (active) {
        gdouble   xc,
                  yc;

        if (gdk_event_get_coords (theMessage->theEvent,&xc,&yc)) {
            switch (theMessage->theEvent->type) {
            case GDK_BUTTON_PRESS:
            case GDK_2BUTTON_PRESS: {
                GdkEventButton * bevent = (GdkEventButton*)theMessage->theEvent;

                long  globalPt_x = xc,
                      globalPt_y = yc;

                localPt_x = globalPt_x-rel.left-parentWindow->allocation.x;
                localPt_y = globalPt_y-rel.top-parentWindow->allocation.y;

                vertical = (localPt_x<headerWidth)&&(localPt_y>=(GetSlotHeight()+1));

                if (((GdkEventButton*)bevent)->button == 1) {
                    if (vertical) {
                        ProcessVSelectionChange (localPt_x,localPt_y,bevent->state & GDK_SHIFT_MASK,bevent->state & GDK_CONTROL_MASK, false, bevent->type == GDK_2BUTTON_PRESS);
                    } else {
                        ProcessSelectionChange  (localPt_x,localPt_y,bevent->state & GDK_SHIFT_MASK,bevent->state & GDK_CONTROL_MASK);
                    }
                } else {
                    if ((((GdkEventButton*)bevent)->button == 2 || ((GdkEventButton*)bevent)->button == 3 )&& (vertical&&vselection.lLength || !vertical &&selection.lLength)) {
                        ProcessContextualPopUp (globalPt_x, globalPt_y);
                        return true;
                    }
                }
            }
            break;
            case GDK_BUTTON_RELEASE:
            case GDK_LEAVE_NOTIFY: {
                if (amScrolling) {
                    if (messageRecipient) {
                        ((_HYTWindow*)messageRecipient)->trackMouseComponent = (Ptr)nil;
                    }

                    /*gdk_pointer_ungrab (((GdkEventButton*)theMessage->theEvent)->time);*/

                    if  (vertical) {
                        _HYRect invalRectH = {parentWindow->allocation.x+rel.left,
                                              parentWindow->allocation.y+rel.top+(GetSlotHeight()+1)+1,rel.left+headerWidth,rel.bottom-HY_SCROLLER_WIDTH
                                             };
                        GdkRectangle irect = HYRect2GDKRect(invalRectH);
                        irect.x+=parentWindow->allocation.x;
                        irect.y+=parentWindow->allocation.y;
                        gdk_window_invalidate_rect (parentWindow->window, &irect, false);
                        if ( localPt_x<headerWidth && localPt_x>0 && lastClick>-2) {
                            MoveSpecies (firstClick+originalStart,lastClick+startRow);
                        }
                    }
                    amScrolling = false;
                }
                break;
            }

            case GDK_MOTION_NOTIFY: {
                GdkEventMotion * motEvent = (GdkEventMotion*)theMessage->theEvent;
                if (motEvent->state & GDK_BUTTON1_MASK) {
                    if (amScrolling) {
                        gint mousePt_x = motEvent->x - rel.left - parentWindow->allocation.x,
                             mousePt_y = motEvent->y - rel.top  - parentWindow->allocation.y;

                        if (vertical) { // vertical scrolling
                            long  wHeight = rel.bottom-rel.top-HY_SCROLLER_WIDTH,
                                  slotHeight = GetSlotHeight();

                            if ( mousePt_y <  GetSlotHeight()+1 || localPt_y != mousePt_y || mousePt_y>wHeight ) {
                                localPt_x = mousePt_x;
                                localPt_y = mousePt_y;

                                if (mousePt_y>wHeight) {
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
                                        forceUpdateForScrolling = true;
                                        _MarkForUpdate();
                                        forceUpdateForScrolling = false;
                                        lastClick = -2;
                                    }
                                    break;
                                } else {
                                    mousePt_y-=(GetSlotHeight()+1);
                                    if (mousePt_y<=slotHeight) {
                                        if (mousePt_y>=0) {
                                            if (mousePt_y<slotHeight/2) {
                                                mousePt_y = -1;
                                            } else {
                                                mousePt_y = 0;
                                            }
                                        } else {
                                            // scroll up
                                            if (startRow>0) {
                                                startRow--;
                                                endRow--;
                                                _SetVScrollerPos(((double)MAX_CONTROL_VALUE*startRow)/(speciesIndex.lLength-endRow+startRow+1));
                                                BuildPane();
                                                forceUpdateForScrolling = true;
                                                _MarkForUpdate();
                                                forceUpdateForScrolling = false;
                                                lastClick = -2;
                                            }
                                            break;
                                        }
                                    } else {
                                        mousePt_y=(mousePt_y-(GetSlotHeight()+1))/slotHeight;
                                    }
                                }

                                if ( mousePt_y<-1 || mousePt_y>= endRow-startRow ) {
                                    break;
                                }

                                if (mousePt_y!=lastClick) {
                                    GdkDrawable * tempDr = GDK_DRAWABLE(parentWindow->window);
                                    GdkGC * tempGC = gdk_gc_new (tempDr);

                                    GdkColor black = HYColorToGDKColor((_HYColor) {
                                        0,0,0
                                    });
                                    gdk_gc_set_foreground (theContext, &black);

                                    gdk_gc_set_function (tempGC, GDK_INVERT);
                                    gdk_gc_set_line_attributes (tempGC, 2, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_MITER);
                                    if (lastClick>=-1) {
                                        lastClick = (GetSlotHeight()+1)+slotHeight*(lastClick+1)+rel.top+1;
                                        gdk_draw_line  (tempDr,tempGC,parentWindow->allocation.x+rel.left+1,parentWindow->allocation.y+lastClick,
                                                        parentWindow->allocation.x+rel.left+headerWidth-1,parentWindow->allocation.y+lastClick);

                                    }

                                    lastClick = mousePt_y;

                                    if (lastClick+startRow != firstClick+originalStart) {
                                        mousePt_y = (GetSlotHeight()+1)+slotHeight*(lastClick+1)+rel.top+1;
                                        gdk_draw_line  (tempDr,tempGC,parentWindow->allocation.x+rel.left+1,parentWindow->allocation.y+mousePt_y,
                                                        parentWindow->allocation.x+rel.left+headerWidth-1,parentWindow->allocation.y+mousePt_y);
                                    }
                                    g_object_unref (tempGC);
                                }
                            }
                            return true;
                        } else { // horizontal scrolling
                            long    rightWindowBound = _HYCanvas::GetMaxW()-HY_SCROLLER_WIDTH;
                            guint32 serverTime = gdk_x11_get_server_time (parentWindow->window);
                            if ( mousePt_x<headerWidth && startColumn>0 || localPt_x!=mousePt_x || mousePt_x> rightWindowBound) {
                                forceUpdateForScrolling = true;
                                if (mousePt_x<headerWidth && startColumn>0) {
                                    gint      wx, wy;
                                    gdk_window_get_origin (parentWindow->window,&wx,&wy);
                                    wx += parentWindow->allocation.x+rel.left;
                                    wy += parentWindow->allocation.y+rel.top;
                                    do {
                                        guint32 serverTime2 = gdk_x11_get_server_time (parentWindow->window);
                                        if (serverTime2-serverTime < 100) {
                                            ProcessSelectionChange (mousePt_x,mousePt_y,true,true,true);
                                        }
                                        GdkModifierType keyDown;
                                        gdk_display_get_pointer (gdk_display_get_default(),NULL,&mousePt_x,&mousePt_y,&keyDown);
                                        mousePt_x -= wx;
                                        mousePt_y -= wy;
                                        serverTime = serverTime2;
                                        if ((keyDown & GDK_BUTTON1_MASK)==0) {
                                            break;
                                        }

                                        gtk_main_iteration_do(true);
                                    } while (mousePt_x<headerWidth && startColumn>0);
                                } else {
                                    if (mousePt_x> rightWindowBound) {
                                        gint      wx, wy;
                                        gdk_window_get_origin (parentWindow->window,&wx,&wy);
                                        wx += parentWindow->allocation.x;
                                        wy += parentWindow->allocation.y;
                                        do {
                                            guint32 serverTime2 = gdk_x11_get_server_time (parentWindow->window);
                                            if (serverTime2-serverTime < 100) {
                                                ProcessSelectionChange (mousePt_x+HY_SCROLLER_WIDTH,mousePt_y,true,true,true);
                                            }
                                            GdkModifierType keyDown;
                                            gdk_display_get_pointer (gdk_display_get_default(),NULL,&mousePt_x,&mousePt_y,&keyDown);
                                            mousePt_x -= wx;
                                            mousePt_y -= wy;
                                            serverTime = serverTime2;
                                            if ((keyDown & GDK_BUTTON1_MASK)==0) {
                                                break;
                                            }

                                            gtk_main_iteration_do(true);
                                        } while (mousePt_x> rightWindowBound && endColumn < columnStrings.lLength);
                                    } else if (serverTime-gdk_event_get_time(theMessage->theEvent) < 100) {
                                        ProcessSelectionChange (mousePt_x,mousePt_y,true,true,true);
                                    }
                                }

                                forceUpdateForScrolling = false;
                                localPt_x = mousePt_x;
                                localPt_y = mousePt_y;
                            }
                            return true;
                        }
                    } else {

                        /*gdk_pointer_grab (parentWindow->window,false,
                                  GDK_POINTER_MOTION_MASK | GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK,
                                  NULL, NULL, (((GdkEventMotion*)theMessage->theEvent))->time);*/

                        if (messageRecipient) {
                            ((_HYTWindow*)messageRecipient)->trackMouseComponent = (Ptr)((_HYComponent*)this);
                        }

                        amScrolling = true;
                        originalStart = startRow,
                        originalSpan  = endRow-startRow;
                        lastClick = -2;
                        firstClick = (localPt_y-(GetSlotHeight()+1))/GetSlotHeight();
                        return true;
                    }
                }
                break;
            }
            }
        }
    }


    return false;
}


//EOF
