/*
    Model Display/Edit Window

    Sergei L. Kosakovsky Pond,
        August-September 2001.

*/

#ifndef _HYMODELWINDOW_
#define _HYMODELWINDOW_

#include "HYTableWindow.h"
#include "HYTableComponent.h"
#include "parser.h"

#define  HY_MODELWINDOW_ICON_BASE        7000
#define  HY_MODEL_WINDOW_MENU_ID         8111

#define  HY_MODEL_TYPE_NUC               0
#define  HY_MODEL_TYPE_AA                1
#define  HY_MODEL_TYPE_DINUC             2
#define  HY_MODEL_TYPE_CODON             3

#define  MODEL_PARAMETER_ROW     0
#define  MODEL_MATRIX_HEADER_ROW 1
#define  MODEL_MATRIX_ROW        2
#define  MODEL_CLASS_ROW         3
#define  MODEL_BUTTON_ROW        4

//class     _HYModelParameterDialog;

//__________________________________________________________________

class _HYModelWindow: public _HYTWindow
{

public:

    _HYModelWindow          (_String,_List*,char );
    virtual             ~_HYModelWindow         (void)
    {}

    virtual bool        ProcessEvent            (_HYEvent*);
    virtual void        Activate                (void);
    virtual void        Deactivate              (void);

    virtual bool        ConfirmClose            (void);

    long                FindVarID               (_String&);
    _String&            RetrieveVarID           (long);
    char                VarKind                 (long);
    void                InsertVarID             (_String&, char);
    bool                DeleteVarID             (long, bool clean = false);
    void                RenameVarID             (_String&, long, char);
    long                MenuChoiceToID          (long);
    void                DoFitColumns            (void);
    void                DoCopyCell              (void);
    void                DoPasteToCells          (void);
    void                DoClearCells            (void);
    void                DoSelectCells           (void);
    void                DoMultiplyCells         (void);
    void                DoSelectAll             (void);
    void                DoSelectAllEmpty        (void);
    void                DoMatchCells            (_String*);
    void                DoSave                  (char);
    bool                DoEditModelName         (void);
    void                SymmetrizeSelection     (_SimpleList&, bool pad = false);
    void                UpdateClassMenus        (void);
    void                ProcessClassMenu        (void);
    void                BuildTemplates          (_SimpleList*,char);
    void                GrabEFVs                (char);
    void                GrabRateVariation       (void);
    bool                CopyMatrix              (_String&);
    bool                CopySimpleMatrix        (_Matrix*);
    void                SetEFVVector            (_String&);
    void                HandleEFVEditEvent      (long);
    void                UpdateEFVLastEntry      (void);
    void                SetEFVChoice            (long);
    void                SyncEditBox             (void);

    bool                AddFormulaParametersToList
    (_Formula&,bool=false);

    void                ProcessTemplateSelection(long);
    void                DoScrollToSelection     (_SimpleList&);
    long                ParameterEditBox        (long);
    virtual void        _SetMenuBar             (void);
    virtual void        _UnsetMenuBar           (void);
    virtual bool        _ProcessMenuSelection   (long);
    bool                _CheckClipboard         (void);
    void                _SetClipboard           (void);
    void                SetTaint                (bool t) {
        taint = t;
    }

    void                CopyEFVMatrix           (_Matrix*);
    virtual bool        IsSaveEnabled     (void) {
        return true;
    }
    virtual bool        IsPrintEnabled    (void) {
        return true;
    }

    friend              class _HYModelParameterDialog;

private:

    void                UpdateParameterName     (long, long);
    bool                RemoveParameterName     (long);
    void                HandleCellEditEvent     (long);
    void                SetParameterMenuState   (long);
    void                PrepareParameterMenu    (void);
    void                UpdateButtonState       (void);
    void                _UpdateEditMenu         (bool,bool);
    void                CheckDependencies       (void);

    _List               globalVariables,
                        localVariables,
                        categoryVariables,
                        stateLabels,
                        subClasses,
                        matrixTemplates,
                        EFVOptions,
                        rateOptions;

    _SimpleList         genCode;

    _String             clipboardString;

    long                lastParameterChoice,
                        EFVChoice,
                        rateChoice;

    bool                taint;
    char                type;
};


//__________________________________________________________________

class _HYModelWindowDialog: public _HYTWindow
{

public:

    _HYModelWindowDialog     (long*,long*,_String*,long*);
    virtual     ~_HYModelWindowDialog    (void) {}

    virtual bool ProcessEvent            (_HYEvent*);

    void SetModelChoice          (_String*);
private:

    void FindModelsOfType        (char);
    void SetModelOptions         (_String*);
    void GetGeneticCodes         (void);
    bool SwitchModelOptions      (char);
    void SwitchModelSelection    (long);

    long        *mDim,
                *gCode,
                *pl3;

    _String     *sp1;
};

//__________________________________________________________________

class _HYModelParameterDialog: public _HYTWindow
{

public:

    _HYModelParameterDialog      (_HYModelWindow*, long, long*);
    virtual     ~_HYModelParameterDialog     (void)
    {}

    virtual bool ProcessEvent                    (_HYEvent*);

    long*                            result,
                                     inID,
                                     orb,
                                     rb;

    _HYModelWindow*                  parentModel;
};

//__________________________________________________________________


_String*                DefStringForGCode       (_SimpleList*);
bool                    OpenModelFromFile       (_String&);
bool                    OpenModelFromMatrix     (_String&,bool = true);
bool                    OpenModelFromCalcNode   (_String&,bool);
void                    SpawnStandardLabels     (char,_List&);
_HYModelWindow          *SpawnNewModel          (long, long, _String&, long);

#endif

//EOF
