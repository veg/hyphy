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

#ifndef _HYPCOMPONENT_
#define _HYPCOMPONENT_
//#pragma once
#include "HYBaseGUI.h"

#include <gtk/gtk.h>
#define  _HY_SANS_FONT  "HYPHY_Sans"
#define  _HY_MONO_FONT  "HYPHY_Mono"
#define  _HY_SERF_FONT  "HYPHY_Serif"

//__________________________________________________________________

class _HYPlatformComponent
{

public:

    _HYPlatformComponent(void);
    _HYPlatformComponent(_HYRect,Ptr);
    // settings

    virtual ~_HYPlatformComponent() {};

    virtual void        Duplicate (BaseRef);

    virtual void        _CleanUp   (void);

    virtual void        _SetDimensions (_HYRect,_HYRect);
    virtual void        _SetVisibleSize(_HYRect);

    virtual void        _Paint (Ptr);
    virtual void        _Update (Ptr);
    virtual bool        _ProcessOSEvent (Ptr);
    virtual _HYRect     _VisibleContents(Ptr);
    virtual void        _MarkForUpdate (void);
    virtual void        _MarkContentsForUpdate (void);
    virtual long        _GetScrollerPos  (GtkWidget*);
    virtual long        _GetHScrollerPos (void);
    virtual long        _GetVScrollerPos (void);
    virtual void        _SetScrollerPos  (GtkWidget*,long);
    virtual void        _SetHScrollerPos (long);
    virtual void        _SetVScrollerPos (long);
    virtual void        _Activate (void);
    virtual void        _Deactivate(void);
    virtual void        _ComponentMouseExit (void) {}



    GtkWidget     *parentWindow,
                  *vScroll,
                  *hScroll;

    _HYRect       rel;

    bool          activationFlag;

    long          lastHScroll,
                  lastVScroll;
};

void            AlignRectangle                          (_HYRect&, GdkRectangle&, unsigned char);
gboolean        hyphy_event_component_callback          (GtkWidget *, GdkEvent*, gpointer);


extern          bool    forceUpdateForScrolling;
extern          double  fontConversionFactor;

//__________________________________________________________________

#endif

//EOF
