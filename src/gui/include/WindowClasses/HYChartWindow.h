/*
    Chart Table

    Sergei L. Kosakovsky Pond, December 2001.
*/

#ifndef _HYCHARTWINDOW_
#define _HYCHARTWINDOW_
//#pragma once
#include "HYTableWindow.h"
#include "HYTableComponent.h"
#include "category.h"

#define  HY_CHART_ICON_BASE                         6300
#define  HY_CHART_WINDOW_MENU_ID                    8115
#define  HY_CHART_WINDOW_HMENU_ID                   155

#define  HY_CHART_MENU_ROW                          0
#define  HY_CHART_CHART_ROW                         1
#define  HY_CHART_HEADER_ROW                        2
#define  HY_CHART_TABLE_ROW                         3

#define  HY_CHART_TABLE_HEIGHT                      150

#define  HY_CHART_PC_TILTED                         0x01
#define  HY_CHART_PC_LABEL_TOP                      0x02
#define  HY_CHART_PC_LABEL_BOTTOM                   0x04
#define  HY_CHART_PC_TITLE_STRING                   0x08

#define  HY_CHART_COLOR_COUNT                       10

#define  HY_CHART_BC_LOG_X                          0x01
#define  HY_CHART_BC_LOG_Y                          0x02
#define  HY_CHART_BC_TITLE_STRING                   0x04

#define  HY_CHART_LEGEND_NONE                       0
#define  HY_CHART_LEGEND_TOP_LEFT                   1
#define  HY_CHART_LEGEND_TOP_MID                    2
#define  HY_CHART_LEGEND_TOP_RIGHT                  3
#define  HY_CHART_LEGEND_BOT_LEFT                   4
#define  HY_CHART_LEGEND_BOT_MID                    5
#define  HY_CHART_LEGEND_BOT_RIGHT                  6
#define  HY_CHART_LEGEND_MID_LEFT                   7
#define  HY_CHART_LEGEND_MID_RIGHT                  8

//__________________________________________________________________

class _HYChartWindow: public _HYTWindow
{

public:

    _HYChartWindow          (_String,_List&,_Matrix&,_HYWindow* pw = nil);
    virtual             ~_HYChartWindow         ();

    virtual bool        ProcessEvent            (_HYEvent*);

    virtual void        _SetMenuBar             (void);
    virtual void        _UnsetMenuBar           (void);
    virtual bool        _ProcessMenuSelection   (long);
    virtual bool        IsSaveEnabled           (void) {
        return true;
    }
    virtual bool        IsPrintEnabled          (void) {
        return true;
    }
    virtual void        SetChartType            (_String,_String,_String,bool = true);
    void                SetTable                (_List&,_Matrix&);
    void                SetLabels               (_String,_String,_String,long,_String, long, long, long, _Parameter = 0., _Parameter = 0.);
    virtual bool        _ProcessOSEvent         (Ptr);
    virtual void        SetFont                 (_HYFont&);
    virtual void        HandleCopyPaste         (bool);
    void                ExecuteProcessor        (long);
    long                Get3DScaling            (bool);

    void                SetProjection           (_Parameter, _Parameter, _Parameter);
    void                SetFonts                (_List*);
    void                SetColors               (_List*);
    void                ToggleSuspend           (bool);

    void                DrawChart               (_HYRect* printRect = nil);
    _SimpleList&        GetColors               (void) {
        return theColors;
    }
    void                _CopyChart              (void);
    void                SetRowHeaders           (_List&);
    void                RenameChartWindow       (void);
    bool                SetUserBounds           (_Parameter = 0.0, _Parameter = 0.0, bool = false);
protected:
    virtual void        DoSave                  (long, _String* distrib = nil);
private:

    void                DoPrint                 (long);
    void                DoChangeFont            (long);
    void                HandleCellEditEvent     (long);
    void                HandleChartOptions      (void);
    bool                NeedToRedrawChart       (long);
    void                SetYPullDown            (void);
    void                GenerateBarChart        (_SimpleList&, _HYRect plotRect, char options = 0);
    void                Generate3DChart         (_SimpleList&, _HYRect plotRect, char options = 0);
//  void                GenerateStackedBarChart (_SimpleList&, _HYRect plotRect);
    void                GeneratePieChart        (long,_HYRect);
    void                SetDataMatrix           (_Matrix&);
    void                SetColumnHeaders        (_List&);
    void                SetPullDowns            (_List&);
    void                HandleChartTypeChange   (long);
    void                _PrintChart             (void);
    _Matrix*            ComputeProjectionMatrix (void);
    void                DrawAParallelogram      (_Matrix&,char,_Parameter,_HYRect&);
    void                DrawAPlane              (_Matrix&,char,_Parameter,_HYRect&);
    void                DrawALine               (_Matrix&,char,_Parameter,_HYRect&, long, _String* = nil, long = 0, long = 0);
    void                PaintAFace              (_Matrix&,long,long,long,long,_HYRect&,_Parameter = 1.0, char = 0);
    _Matrix*            ComputeProjection       (_Matrix&);
    void                ComputeProjectionSettings
    (void);
    void                DrawLegend              (_HYRect&,long,bool=false);


    _HYWindow*          parentWindow;

    _HYFont             labelFont1,
                        labelFont2,
                        labelFont3,
                        headerFont;

    _HYColor            backColor1,
                        backColor2;

    _SimpleList         ySeries,
                        theColors;

    _Parameter          projectionSettings[7],
                        oR,
                        xyAngle,
                        zAngle,
                        xShift,
                        yShift,
                        xScale,
                        yScale,
                        userMin,
                        userMax;

    _String             xLabel,
                        yLabel,
                        zLabel;

    _Formula            overlayPlot;

    char                showLegend,
                        whichFont;
    bool                suspendDraw;

    long                xAxis3DScale,
                        yAxis3DScale,
                        surfaceDivs;

    _Matrix*            projectionMatrix;

};

//__________________________________________________________________

class _HYDistributionChartWindow: public _HYChartWindow
{

public:
    _HYDistributionChartWindow          (_String,_List&,_Matrix&,_List&_,_HYWindow* pw = nil);
    virtual             ~_HYDistributionChartWindow         (void) {};

    virtual bool        ProcessEvent            (_HYEvent*);

    virtual void        _SetMenuBar             (void);
    virtual void        _UnsetMenuBar           (void);
    virtual bool        _ProcessMenuSelection   (long);
    virtual void        AddVariable             (_String* = nil);
    virtual void        SetAtoms                (_Matrix&,_List&);

private:

    _List*              BuildAtoms              (_Matrix&,_List&);
    virtual void        DoSave                  (long, _String* distrib = nil);
    virtual void        RemoveVariable          (void);
    virtual void        HandleCatPostProcessor  (long);
    void                ShowMarginals           (void);

    virtual void        ProduceDistribution     (_SimpleList&, _SimpleList&, _List&, _Formula&, long, _Parameter);
    virtual _Matrix*    ComputeConditionals     (_SimpleList&, _SimpleList*);
    virtual _Matrix*    ComputeConditionals     (long);
    virtual _Matrix*    MakeCDF                 (long);

    _List               atoms,
                        atomNames,

                        derived,
                        derivedDependencies,
                        derivedNames,
                        derivedRemaps;

    _String             savedExpression;

    _Matrix             marginals;

    _SimpleList         atomSizes,
                        atomMultiples;

};

//__________________________________________________________________

class       _HYChartOptionDialog: public _HYTWindow
{

public:

    _HYChartOptionDialog     (_List*,long*,bool*, _SimpleList*,_HYChartWindow*);
    virtual         ~_HYChartOptionDialog    (void) {}

    virtual bool    ProcessEvent             (_HYEvent*);
    void    DrawColorCanvas          (void);

    _List*          args;
    long *          msel;
    bool *          res;
    _HYChartWindow* parWin;
    _SimpleList*    colorsPicked;

};

bool                    ReadDataFromFile        (_String,char,_Matrix&, _List&);

extern      _List       chartProcessors,
            distribProcessors;

extern      void        ReadChartProcessors     (bool = false);

#endif

//EOF
