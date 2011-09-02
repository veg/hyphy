/*
    A tree window - tree canvas, scrolling pane and a button bar

    Sergei L. Kosakovsky Pond, May 2000.
*/

#ifndef _HYTREEPANEL_
#define _HYTREEPANEL_
//#pragma once
#include "HYTableWindow.h"
#include "HYCanvas.h"
#include "HYDialogs.h"
#include "calcnode.h"
#include "classes.h"

#define  HY_TREEPANEL_ALIGNED               0
#define  HY_TREEPANEL_SCALED                1
#define  HY_TREEPANEL_SQUARE                2
#define  HY_TREEPANEL_STRAIGHT              4
#define  HY_TREEPANEL_ARCS                  8
#define  HY_TREEPANEL_CIRCULAR              1024
#define  HY_TREEPANEL_LABEL1                2048
#define  HY_TREEPANEL_LABEL2                4096
#define  HY_TREEPANEL_LAYOUT_MASK           0xfbf1
#define  HY_TREEPANEL_SCALING_MASK          0xfffd
#define  HY_TREEPANEL_TIP_LABELS            16
#define  HY_TREEPANEL_INT_LABELS            32
#define  HY_TREEPANEL_VERTICAL              64
#define  HY_TREEPANEL_CLIPBOARD_READY       128
#define  HY_TREEPANEL_PROJECTION            256
#define  HY_TREEPANEL_FISHEYECURSOR         512
#define  HY_TREEPANEL_SCALE_TO_WINDOW       8192
#define  HY_TREEPANEL_MAXSIZE               4096
#define  HY_TREEPANEL_NAVSIZE               115
#define  HY_TREEPANEL_DEFSIZE               350
#define  HY_TREEPANEL_MIN_VSPACE            10
#define  HY_TREEPANEL_MAX_VSPACE            500
#define  HY_TREEPANEL_MIN_HSPACE            10
#define  HY_TREEPANEL_SPACE_STEP            5
#define  HY_TREEPANEL_MAX_HSPACE            500
#define  HY_TREEPANEL_VSPACE                30
#define  HY_TREEPANEL_HSPACE                60
#define  HY_TREEPANEL_MARGIN                10
#define  HY_TREEPANEL_LABELSPACING          75
#define  HY_TREEPANEL_NAVSPACING            2
#define  HY_TREEPANEL_TEXTHALFSIZE          4
#define  HY_TREEPANEL_TEXTSIZE              9
#define  HY_TREEPANEL_RULER_COLLAPSED       5
#define  HY_TREEPANEL_RULER_EXPANDED        15
#define  HY_TREEPANEL_ICON_ID               5000
#define  HY_TREEPANEL_MENU_ID               7990
#define  HY_TREEPANEL_HMENU_ID              190
#define  HY_CLIPBOARD_CUT                   1
#define  HY_CLIPBOARD_COPY                  2

#define  HY_WINDOW_KIND_TREE                1



//__________________________________________________________________

class _HYTreePanel: public _HYTWindow
{

    friend class _HYLabelDialog;
    friend class _HYTreeViewDialog;

public:

    _HYTreePanel(_String&,_String&);
    // tree name, tree string (or tree ident)

    virtual ~_HYTreePanel();

    virtual bool        ProcessEvent                    (_HYEvent*);
    virtual bool        ProcessGEvent                   (_HYEvent*);
    bool        BuildTree                       (bool);
    void        RenderTree                      (bool map = true, bool printing = false);
    void        RenderRuler                     (_Parameter hsc = 1.0, bool printing = false, int hshift = 0, int vshift = 0);
    void        RenderNavTree                   (void);
    void        SetTreeString                   (_String&);
    void        SetVariableReference            (_String&);
    void        SetFlags                        (unsigned int);
    unsigned int        GetFlags                        (void) {
        return treeFlags;
    }
    void        SetScaleVariable                (_String&);
    void        SetNavRectCenter                (int,int);
    virtual void        Update                          (Ptr);
    virtual void        Paint                           (Ptr);
    long        FindSelection                   (long,long,char);
    long        FindSelectedBranch              (long,long);
    _String     GetTreeString                   (_AVLListXL* = nil);
    void        SwapSelectedSubTrees            (void);
    void        CollapseSelectedBranch          (void);
    void        JoinSelectedBranches            (void);
    void        GraftATip                       (void);
    void        MoveSubTree                     (void);
    void        RerootTree                      (void);
    void        DeleteCurrentSelection          (void);
    virtual bool        _ProcessMenuSelection           (long);
    virtual void        _SetMenuBar                     (void);
    virtual void        _UnsetMenuBar                   (void);
    _TheTree*   LocateMyTreeVariable            (void);
    void        SelectAllBranches               (void);
    void        SelectIncompleteBranches        (void);
    void        SelectBranchesWithoutModel      (void);
    void        SelectBranchesByName            (void);
    void        SelectBranchesByLength          (void);
    void        DisplayLikelihoodValue          (void);
    void        CutSelectionToClipboard         (bool cut = true);
    void        PasteClipboardTree              (void);
    void        FlipSelectedBranches            (void);
    void        SelectEntireSubtree             (void);
    void        UpdateScalingVariablesList      (void);
    void        UndoLastOperation               (void);
    void        SelectRangeAndScroll            (_List&, bool = false);
    void        ScrollToSelection               (_Parameter, _Parameter);
    void        DisplayParameterTable           (void);
    void        SetFont                         (_HYFont&);
    virtual bool        IsSaveEnabled                   (void) {
        return true;
    }
    virtual bool        IsPrintEnabled                  (void) {
        return true;
    }
    void        ToggleScaleOption               (void);
    void        ShowModelMatrix                 (bool);
    void        RenderTree2                     (void);
    virtual char        WindowKind                      (void) {
        return HY_WINDOW_KIND_TREE;
    }
    void        _DisplayBranchFloater           (void);
    bool        HasToolTip                      (void) {
        return toolTipBounds.left>0;
    }
    virtual void        _HandleIdleEvent                (void);
    void        HandleSelection                 (char);
    void        HandleComparison                (char);
    void        GenerateDistanceTable           (char);
    void        InvertSelection                 (void);
    void        MatchToDataSet                  (void);
    void        GrowSelection                   (node<nodeCoord>* = nil, bool = false);
    void        HandleContextPopup              (long, long);

private:

    friend      class _HYNodeInfoDialog;

    void        ExecuteProcessor                (long);
    void        HandleTreeSave                  (long,_String);
    void        HandleSearchAndReplace          (bool);
    void        SelectSubTree                   (node<nodeCoord>*, bool);
    void        Relabel                         (_TheTree*);
    void        InvokeNodeEditor                (void);
    void        PaintSquareBranchLabels         (_HYCanvas*,node<nodeCoord>*,bool);
    long        PaintSquare                     (_HYCanvas*,node<nodeCoord>*);
    long        PaintVSquare                    (_HYCanvas*,node<nodeCoord>*);
    void        PaintStraight                   (_HYCanvas*,node<nodeCoord>*);
    void        PaintRadial                     (_HYCanvas*,node<nodeCoord>*);
    void        PaintArcs                       (_HYCanvas*,node<nodeCoord>*);
    long        PaintNSquare                    (_HYCanvas*,node<nodeCoord>*,_Parameter, _Parameter);
    long        PaintNVSquare                   (_HYCanvas*,node<nodeCoord>*,_Parameter, _Parameter);
    void        PaintNStraight                  (_HYCanvas*,node<nodeCoord>* ,_Parameter, _Parameter);
    void        PaintNRadial                    (_HYCanvas*,node<nodeCoord>* ,_Parameter, _Parameter);
    void        PaintNArcs                      (_HYCanvas*,node<nodeCoord>*,_Parameter, _Parameter);
    void        Convert2ScreenCoordinates       (_Parameter,_Parameter,_Parameter,node<nodeCoord>*);
    void        ShiftScreenCoordinates          (_Parameter,_Parameter,node<nodeCoord>*);
    void        StraightenEdges                 (node<nodeCoord>*,_Parameter,_Parameter,int);
    _Parameter  PreStraightenEdges              (node<nodeCoord>*);
    void        SpaceNodes                      (node<nodeCoord>*,_Parameter&);
    _HYRect     ComputeNavRect                  (void);
    bool        SetHSpace                       (int,bool force = false);
    bool        SetVSpace                       (int,bool force = false, bool render = true);
    bool        IsVertical                      (void) {
        return treeFlags&HY_TREEPANEL_VERTICAL;
    }
    void        SetVertical                     (bool);
    void        ResizeTreeCanvas                (int, int);
    bool        CheckIfNeedUnscaled             (void);
    bool        KillLikeFunc                    (char);
    void        PreserveSelection               (_SimpleList&);
    void        RestoreSelection                (_SimpleList&);
    void        FishEyeProjection               (long,long,long,long,node<nodeCoord>*);
    long        CircularLayoutPass1             (node<nodeCoord>*);
    void        CircularLayoutPass2             (node<nodeCoord>*,_Parameter,_Parameter,_Parameter,_Parameter,_Parameter,long*);
    void        CircularLayoutPass3             (node<nodeCoord>*,_Parameter,_Parameter,_Parameter,_Parameter);
    void        HandleViewOptions               (void);
    void        HandleLabels                    (bool);
    void        HighlightSelectionInDataPanel   (void);
    void        ShowSelectionInAnotherTree      (void);

    void        _PaintNavRect                   (void);
    void        _PaintLFStatus                  (Ptr = nil);
    virtual bool        _ProcessOSEvent                 (Ptr);
    void        _PrintTree                      (long = -1, long = -1);
    void        _UpdateOperationsMenu           (void);
    long        FindUnusedSuffix                (_String);
    void        UpdateLnLikelihood              (void);
    void        RecalculateLikelihood           (void);
    void        ClearClipboardContents          (void);
    bool        ShrinkDataFilter                (bool doit = false);
    void        PrepareUndoData                 (char);
    void        FlushUndoData                   (void);
    void        HandleFisheyeButton             (void);
    void        FitToWindow                     (bool force = false);
    long        GetUniversalSaveOptions         (_List&);
    void        dumpCoordTree                   (void);
    void        DrawStraightEdge                (_HYCanvas*, _HYRect&, node<nodeCoord>*);

    node<nodeCoord>* coordTree,
         * selectionTop;

    node<long>     * undoTree;

    _String          treeName,
                     scaleVariable,
                     undoLFName,
                     branchVar1,
                     branchVar2;

    char             undoCode,
                     topAlign,
                     bottomAlign,
                     labelDigits1,
                     labelDigits2,
                     branchWidth;

    int              vSpace,
                     hSpace,
                     tips,
                     textSpace,
                     saveMouseH,
                     saveMouseV,
                     windowTextMarginH,
                     windowTextMarginV;

    long             likeFuncID;

    unsigned long    lastSave;

    unsigned int     treeFlags;

    _HYRect          navRect,
                     toolTipBounds;

    _Parameter       treeLength,
                     hScale,
                     vScale,
                     distortion,
                     arcStart,
                     arcEnd;

    _SimpleList      currentSelection,
                     dubiousNodes,
                     undoNodeList;

    _HYFont          treeLabelFont,
                     branchLabel1,
                     branchLabel2;

    Ptr              undoLFPointer;
};


//__________________________________________________________________

class _HYLabelDialog: public _HYFontDialog
{

public:

    _HYLabelDialog          (_HYTreePanel*,_List&,bool);
    virtual     ~_HYLabelDialog         (void)
    {}

    virtual bool ProcessEvent            (_HYEvent*);

private:

    bool which;

};

//__________________________________________________________________

class _HYTreeViewDialog: public _HYFontDialog
{

public:

    _HYTreeViewDialog           (_HYTreePanel*,bool);
    virtual     ~_HYTreeViewDialog              (void)
    {}

    virtual bool ProcessEvent               (_HYEvent*);

private:

    void            DrawArcPreview              (void);

    bool            which,
                    tb1,
                    tb2,
                    tb3;

    _Parameter      sAng,
                    fAng;

    char            lw;

};

//__________________________________________________________________

class _HYNodeInfoDialog: public _HYTWindow
{

public:

    _HYNodeInfoDialog   (_String&, long*, bool, _TheTree*, _SimpleList*, bool*, bool*, bool*,_HYTreePanel*);
    virtual     ~_HYNodeInfoDialog  (void)
    {}


    virtual bool ProcessEvent       (_HYEvent*);

private:

    long        AddParameterRow         (_String&,_String&,_String&,long);
    void        DeleteParameterRow      (long);
    void        BuildModelParameters    (void);
    void        UpdateParameterValues   (bool);
    bool        SetNewParameterValue    (long, _String*);
    void        FitToHeight             (long);
    bool        SetNodeModel            (_CalcNode*,long,bool = true);
    void        ToggleFormulaBox        (bool,_String);
    void        SetNewUserFormula       (void);
    void        SetNewParameterName     (_String, long);

    long            *modelSelection;
    bool            isSingle,
                    *updateSCV,
                    *incFlag,
                    *result,
                    validTE,
                    formulaBox,
                    taintFla;

    _TheTree*       treeRef;
    _SimpleList*    tSel;
    _String         nodeName,
                    *iNodeName;

    _List           allowableDataNames,
                    sharedUEs;

    _Matrix         parameterValues;
    _SimpleList     parameterFlags,
                    parameterVRefs,
                    nodeVRefs;

    _HYTreePanel*   parentRef;
};

//__________________________________________________________________

extern  node<long>* treePanelClipboardRoot;
extern  char        treePanelClipboardMode;
extern  _List       treeProcessors;

void    RenameNodesInClipboard  (node<long>*);
void    RenameUndoNodes         (node<long>*,_String&);
void    ReadTreeProcessors      (void);

void    KillNodesInClipboard (node<long>*);
node<long>*
PatchNodesFromClipboard (node<long>*, _TheTree*, bool&);
node<long>*
CopyNodesToClipboard (node<long>*, _TheTree*, bool&);

void    CopyNodeParameters          (node<nodeCoord>*, node<nodeCoord>*);

typedef _String (listDescPopulator) (_HYWindow*);

long     SelectOpenWindowObjects    (long, _String, listDescPopulator*,Ptr);
_String  treeSelectorPopulator      (_HYWindow *);
_String  dsSelectorPopulator        (_HYWindow *);



#endif

//EOF
