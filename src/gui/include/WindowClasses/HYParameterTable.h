/*
    Parameter Display Table

    Sergei L. Kosakovsky Pond, January-July 2001.
*/

#ifndef _HYPARAMETERTABLE_
#define _HYPARAMETERTABLE_
//#pragma once
#include "HYTableWindow.h"
#include "HYTableComponent.h"
#include "calcnode.h"
#include "category.h"

#define  HY_PARAMETER_TABLE_ICON_BASE               6016
#define  HY_PARAMETER_TABLE_BUTTON_ICON_BASE        6030
#define  HY_PARAMETER_BOOTSTRAP_BUTTON_ICON_BASE    6040
#define  HY_PARAMETER_TABLE_VIEW_LOCAL              0x01
#define  HY_PARAMETER_TABLE_VIEW_GLOBAL             0x02
#define  HY_PARAMETER_TABLE_VIEW_CONSTRAINED        0x04
#define  HY_PARAMETER_TABLE_VIEW_CATEGORY           0x08
#define  HY_PARAMETER_TABLE_VIEW_ALPHABETICAL       0x10
#define  HY_PARAMETER_TABLE_VIEW_TREES              0x20
#define  HY_PARAMETER_TABLE_SNAPSHOT                0x80
#define  HY_PARAMETER_TABLE_VIEW_ALL                0xFF

#define  HY_PARAMETER_TABLE_MENU_ID                 7995
#define  HY_PARAMETER_TABLE_HMENU_ID                172

#define  HY_PARAMETER_TABLE_TABLE_ROW               2
#define  HY_PARAMETER_TABLE_BUTTON_ROW              0

#define  HY_WINDOW_KIND_PARAMETERTABLE              2


//__________________________________________________________________

class _HYParameterTable: public _HYTWindow
{

public:

    _HYParameterTable       (_String,long);
    virtual             ~_HYParameterTable      ();

    virtual bool        ProcessEvent            (_HYEvent*);
    virtual bool        ProcessGEvent           (_HYEvent*);

    virtual void        _SetMenuBar             (void);
    virtual void        _UnsetMenuBar           (void);
    //virtual bool      _ProcessOSEvent         (Ptr);
    virtual bool        _ProcessMenuSelection   (long);
    void        _UpdateViewMenu         (void);

    BaseRef             RetrieveIthRowVar       (long);
    void                OpenCategoryPlotWindow  (_CategoryVariable*);
    void                SetRowValues            (long, _Variable*,_HYTable*);
    void                OptimizeLikelihoodFunction
    (void);
    void                RefreshTheTable         (void);
    void                ToggleConstrainedView   (char flag = HY_PARAMETER_TABLE_VIEW_CONSTRAINED);
    void                TaintTheLF              (void);
    void                SetSelectedRows         (_List&);
    void                SelectAll               (void);
    void                UndoCommand             (void);
    virtual bool        IsSaveEnabled     (void) {
        return true;
    }
    virtual bool        IsPrintEnabled    (void) {
        return true;
    }
    virtual char        WindowKind        (void) {
        return HY_WINDOW_KIND_PARAMETERTABLE;
    }

private:

    void                UpdateLogLikelihood     (bool postEvent = true);
    void                ConstructTheTable       (void);
    char                GetViewOptions          (void);
    void                HandleCellPullDown      (long,long);
    void                HandleHeaderPullDown    (long,long);
    void                HandleCellEditEvent     (long);
    void                HandleSaveLF            (void);
    void                HandleRestoreLF         (_String*);
    void                HandleSetHypothesis     (_String*,bool);
    void                HandleDeleteLF          (_String*);
    void                HandleLRT               (void);
    void                HandleSelectParameters  (void);
    void                HandleOpenInChart       (void);
    void                HandleBootstrap         (bool);
    void                PrepareLFMenu           (void);
    void                UpdateSelectionDependentButtons
    (void);
    void                RefreshParameterValues  (void);
    void                DoEqualConstraint       (void);
    void                DoProportionalConstraint(void);
    void                DoProportionalSubtrees  (void);
    void                DoBoundsChange          (void);
    void                DoMolecularClock        (void);
    void                DoSave                  (void);
    void                DoOpenTreeWindow        (void);
    void                DoClearConstraints      (void);
    void                DoEnterConstraint       (void);
    void                DoCleanUp               (void);
    void                UpdateKthRow            (long, bool mark = true);
    void                StretchColumnToFit      (long);
    _CalcNode*          GrabRowsCalcNode        (_String*);
    void                _UpdateUndoMenu         (_String*, _String*);
    void                TakeLFSnapshot          (void);
    void                DoneLFSnapshot          (void) {
        /*viewOptions -= HY_PARAMETER_TABLE_SNAPSHOT;*/
    }
    void                VerifyGlobalVariables   (void);
    long                GatherListOfGlobals     (_List&);
    void                HandleVarianceEstimates (void);
    void                HandleProfilePlot       (void);
    void                HandleCategories        (void);

    _HYGuiObject*       RetrieveParentDataPanel (void);

    // data members
    long                lfID,
                        lastParameterCount;

    _Parameter          lastLFValue;

    char                viewOptions,
                        avViewOptions;

    _List               undoCommands,
                        undoDescriptions;

    _String             lastLFState,
                        stashCustomCommand;

};

//__________________________________________________________________

class _HYBootstrapWindow: public _HYTWindow
{

public:

    _HYBootstrapWindow      (_String,_Parameter);
    virtual             ~_HYBootstrapWindow     (void) {}
    virtual  bool       ProcessEvent            (_HYEvent*);

    virtual  bool       ConfirmClose            (void);
    void                AddIterate              (_Parameter);
    bool                HasStopped              (void) {
        return done;
    }
    void                SetStopped              (bool);
    void                OpenHistogramWindow     (bool, long);
    void                DoSave                  (void);
    virtual bool        IsSaveEnabled           (void) {
        return true;
    }
    virtual bool        IsPrintEnabled          (void) {
        return true;
    }
    virtual bool        _ProcessMenuSelection   (long);

    bool                requestQuit;


protected:

    long                iterateCount,
                        iterateCountP;

    bool                done;

    _Parameter          nullLRT;
};

//__________________________________________________________________

class _HYGeneralBootstrapWindow: public _HYBootstrapWindow
{

public:

    _HYGeneralBootstrapWindow       (_String,long = -1);
    virtual             ~_HYGeneralBootstrapWindow      (void)
    {}

    virtual  bool       ProcessEvent                    (_HYEvent*);
    void                SetNullLF                       (long, bool = false);
    void                HandleBootstrap                 (char);

    virtual  bool       ProcessGEvent                   (_HYEvent*);
    virtual  bool       ConfirmClose                    (void);

    //virtual bool      _ProcessMenuSelection           (long);
    //virtual void      _SetMenuBar                     (long);
    //virtual void      _UnsetMenuBar                   (long);

private:

    void                GetNullStateString              (_String&, long);
    void                GetAltStateString               (_String&, long);
    void                ResetResults                    (void);
    void                BuildCompatibleLF               (bool set = false);
    void                BuildAvailableLF                (void);
    void                UpdateBSButtons                 (void);
    void                UpdateLRValue                   (void);
    long                lfNullID,
                        lfAltID;

    _SimpleList         lfID0,
                        dpID0,
                        lfIDA,
                        dpIDA;

    bool                hasGoIcon;

};

#endif

//EOF
