/*
    A canvas component.
    Added stretchable canvas

    Sergei L. Kosakovsky Pond, May-August 2000.
*/

#ifndef _HYCANVAS_
#define _HYCANVAS_
//#pragma once
#include "HYComponent.h"
#include "HYGraphicPane.h"

#define  HY_SCANVAS_HORIZONTAL  1
#define  HY_SCANVAS_VERTICAL    2

//__________________________________________________________________

class _HYCanvas: public _HYGraphicPane, public _HYComponent
{

public:

    _HYCanvas(_HYRect,Ptr,int,int,int);

    virtual ~_HYCanvas() {};

    virtual bool        ProcessEvent        (_HYEvent*);
    virtual void        _Paint              (Ptr);
    virtual void        _Update             (Ptr);
    virtual bool        _ProcessOSEvent     (Ptr);
    virtual void        SetCanvasSize       (int,int,int);
    virtual _HYRect     GetCanvasSize       (void);
    virtual void        CrossfadeText       (_String& t1, _String& t2, _HYRect&, long d, long n, char, char);
    virtual void        SetMouseClick       (bool mc) {
        doMouseClicks = mc;
    }
    // text currently displayed is morphed into
    // new text, in n steps with d ms per step
    // + alignment


    bool        doMouseClicks;

};

//__________________________________________________________________

class _HYStretchCanvas: public _HYCanvas
{

public:

    _HYStretchCanvas(_HYRect,Ptr,int,int,int,char flags);

    virtual         ~_HYStretchCanvas() {};


    virtual _HYRect GetCanvasSize (void);
    virtual void    SetVisibleSize (_HYRect);
    virtual int     GetMaxH (void);
    virtual int     GetMaxW (void);


    char        canvasFlags,
                canvasDepth;
};

//__________________________________________________________________

#endif

//EOF
