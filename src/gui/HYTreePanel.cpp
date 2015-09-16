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

#include "HYTreePanel.h"
#include "HYCanvas.h"
#include "HYPullDown.h"
#include "HYEventTypes.h"
#include "calcnode.h"
#include "batchlan.h"
#include "likefunc.h"
#include "HYLabel.h"
#include "HYButtonBar.h"
#include "math.h"
#include "HYUtils.h"

#include "HYTextBox.h"
#include "HYButton.h"
#include "HYChartWindow.h"

#include "HYDataPanel.h"
#include "HYParameterTable.h"
#include "HYModelWindow.h"

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

#define     TREE_PANEL_PIXEL_DEPTH  32

_HYColor    navFill   = {212,212,212},
            navSelect = {255,0,255},
            black = {0,0,0};

bool        promptDefModel     = false,
            promptShrinkFilter = false,
            promptKillLF       = false,
            makeUndoTree     = false,
            displayNotUndoable = false,
            brLengthIntOnly    = false,
            updateInDataPanel  = false;

node<long>* treePanelClipboardRoot = nil;
char        treePanelClipboardMode = HY_CLIPBOARD_CUT;

_HYTreePanel* feedbackTreePanel = nil;

extern      _SimpleList       windowObjects;
extern      _String           none;

_String     eSubsScale          ("E[Substitutions]"),
            assValScale           ("Assigned Values"),
            notUndoableWarning    ("The change of branch model will result in loss of all current parameter values and expressions and the likelihood function."),
            donotWarnAgain        ("Do not warn me again"),
            treeStringExport  ("_TREE_STRING_FOR_PROCESSING_"),
            treeStringSelect    ("_TREE_PANEL_SELECTION_"),
            selectByNameSave,
            selectByBLSave,
            sandrString1,
            sandrString2;

_List       treeProcessors;


//__________________________________________________________
_HYTreePanel::_HYTreePanel(_String& title,_String& argument):_HYTWindow (_String ("Tree ")&title)
{
    if (treeProcessors.lLength==0) {
        ReadTreeProcessors ();
    }

    coordTree = nil;
    likeFuncID = -1;
    undoCode = -1;
    branchWidth = 2;

    flags |= HY_WINDOW_STATUS_BAR_LIGHT_LEFT;
    _HYRect         canvasSettings = {50,50,HY_TREEPANEL_DEFSIZE,HY_TREEPANEL_DEFSIZE,HY_COMPONENT_NO_SCROLL};
    _HYCanvas*      c = new _HYCanvas (canvasSettings,GetOSWindowData(),HY_TREEPANEL_DEFSIZE,HY_TREEPANEL_DEFSIZE,TREE_PANEL_PIXEL_DEPTH);
    canvasSettings.top = canvasSettings.bottom = HY_TREEPANEL_RULER_EXPANDED;
    _HYCanvas*      r = new _HYCanvas (canvasSettings,GetOSWindowData(),HY_TREEPANEL_RULER_EXPANDED,HY_TREEPANEL_DEFSIZE,TREE_PANEL_PIXEL_DEPTH);
    canvasSettings.right=canvasSettings.bottom=canvasSettings.left=canvasSettings.top=HY_TREEPANEL_NAVSIZE;
    _HYCanvas*      c1 = new _HYCanvas (canvasSettings,GetOSWindowData(),HY_TREEPANEL_NAVSIZE,HY_TREEPANEL_NAVSIZE,32);
    canvasSettings.left     = 30;
    canvasSettings.top      = 30;
    canvasSettings.bottom   = 30;
    canvasSettings.right    = HY_TREEPANEL_MAXSIZE;
    _HYLabel*       l1      = new _HYLabel    (canvasSettings,GetOSWindowData());
    _HYLabel*       l2      = new _HYLabel    (canvasSettings,GetOSWindowData());
    canvasSettings.left     = 50;
    _HYPullDown*    p1      = new _HYPullDown (canvasSettings,GetOSWindowData());
    _HYPullDown*    p2      = new _HYPullDown (canvasSettings,GetOSWindowData());
    p1->SetMessageRecipient (this);
    p1->SetAlignFlags (HY_ALIGN_LEFT);
    p2->SetMessageRecipient (this);
    p2->SetAlignFlags (HY_ALIGN_LEFT);
    canvasSettings.bottom   = 50;
    canvasSettings.top      = 50;
    _HYLabel*       l3      = new _HYLabel (canvasSettings,GetOSWindowData());
    _HYButtonBar*   bb      = new _HYButtonBar (canvasSettings, GetOSWindowData());
    bb->SetMessageRecipient (this);
    _String tt ("Expand Tree Horizontally");
    bb->AddButton (ProcureIconResource(HY_TREEPANEL_ICON_ID),&tt);
    tt = "Expand Tree Vertically";
    bb->AddButton (ProcureIconResource(HY_TREEPANEL_ICON_ID+1),&tt);
    tt = "Contract Tree Horizontally";
    bb->AddButton (ProcureIconResource(HY_TREEPANEL_ICON_ID+2),&tt);
    tt = "Contract Tree Vertically";
    bb->AddButton (ProcureIconResource(HY_TREEPANEL_ICON_ID+3),&tt);
    tt = "Toggle Horizontal/Vertical View";
    bb->AddButton (ProcureIconResource(HY_TREEPANEL_ICON_ID+4),&tt);
    tt = "Scale to Fit Window";
    bb->AddButton (ProcureIconResource(HY_TREEPANEL_ICON_ID+6),&tt);
    tt = "Fisheye View";
    bb->AddButton (ProcureIconResource(HY_TREEPANEL_ICON_ID+7),&tt);
    bb->MarkAsPullDown (6,true);
    bb->SetAlignFlags (HY_ALIGN_LEFT);
    l3->SetAlignFlags (HY_ALIGN_LEFT);
    l1->SetAlignFlags (HY_ALIGN_LEFT);
    l2->SetAlignFlags (HY_ALIGN_LEFT);
    canvasSettings.top = canvasSettings.bottom = HY_TREEPANEL_NAVSIZE-107;
    canvasSettings.width = HY_COMPONENT_BORDER_B;
    _HYLabel*       l4 = new _HYLabel (canvasSettings,GetOSWindowData());
    bb->SetButtonDim (20);
    bb->SetButtonLayoutW (4);
    AddObject (c,false);
    AddObject (c1,false);
    AddObject (r,false);
    AddObject (p1,false);
    AddObject (l1,false);
    AddObject (p2,false);
    AddObject (l2,false);
    AddObject (bb,false);
    AddObject (l3,false);
    AddObject (l4,false);
    SetTableDimensions (6,3);
    SetCell (0,0,c1);
    SetCell (0,1,l1);
    SetCell (0,2,p1);
    SetCell (1,0,c1);
    SetCell (1,1,l2);
    SetCell (1,2,p2);
    SetCell (2,0,c1);
    SetCell (2,1,l3);
    SetCell (2,2,bb);
    SetCell (3,0,c1);
    SetCell (3,1,l4);
    SetCell (3,2,l4);
    SetCell (4,0,r);
    SetCell (4,1,r);
    SetCell (4,2,r);
    SetCell (5,0,c);
    SetCell (5,1,c);
    SetCell (5,2,c);

    p1->AddMenuItem(_String("Rectangular"),-1);
    /*p1->AddMenuItem(_String("Straight"),-1);
    p1->AddMenuItem(_String("Smooth Edges"),-1);
    p1->AddMenuItem(_String("Radial"),-1);*/
    p2->AddMenuItem(_String("Unscaled"),-1);

    p2->AddMenuItem(eSubsScale,-1);
    p2->AddMenuItem(assValScale,-1);
    l1->SetText (_String("Tree Style"));
    l2->SetText (_String("Branch Scaling"));
    l3->SetText (_String("Zoom/Rotate"));
    l1->SetShadow(true);
    l2->SetShadow(true);
    l3->SetShadow(true);
    canvasSettings.top = 0;
    canvasSettings.left = 0;
    canvasSettings.bottom = HY_TREEPANEL_NAVSIZE;
    canvasSettings.right = 500;
    canvasSettings.width = HY_COMPONENT_NO_SCROLL;
    //p1->SetDimensions(canvasSettings,canvasSettings);
    c1->StartDraw();
    _HYColor col = navFill;
    c1->SetColor (col);
    canvasSettings.right = HY_TREEPANEL_NAVSIZE-1;
    canvasSettings.width = 1;
    c1->FillRect (canvasSettings);
    col.R = col.B = col.G = 0;
    c1->SetColor (col);
    c1->DrawRect (canvasSettings);
    c1->EndDraw();
    col.R = 102;
    col.G = 90;
    col.B = 88;
    canvasSettings.top = canvasSettings.left = 0;
    canvasSettings.bottom = canvasSettings.right = HY_TREEPANEL_MAXSIZE;
    canvasSettings.width = 1;
    p1->SetBackColor (col);
    p2->SetBackColor (col);
    _HYFont  labelFont;

#ifdef __WINDOZE__
    labelFont.size = 10;
    labelFont.face      = "MS Sans Serif";
    treeLabelFont.face  = labelFont.face;
    branchLabel1.face   = labelFont.face;
    branchLabel2.face   = labelFont.face;
#endif

#ifdef __MAC__
    labelFont.size      = 12;
    labelFont.face      = "Geneva";
    treeLabelFont.face  = labelFont.face;
    branchLabel1.face   = "Times";
    branchLabel2.face   = branchLabel1.face;
#endif

#ifdef __HYPHY_GTK__
    labelFont.size      = 10;
    labelFont.face      = _HY_SANS_FONT;
    treeLabelFont.face  = labelFont.face;
    branchLabel1.face   = _HY_SERF_FONT;
    branchLabel2.face   = branchLabel1.face;
#endif

    treeLabelFont.size = 9;
    treeLabelFont.style = HY_FONT_PLAIN;
    branchLabel1.size = 9;
    branchLabel1.style = HY_FONT_PLAIN;
    branchLabel2.size = 9;
    branchLabel2.style = HY_FONT_BOLD;

    topAlign    = 0;
    bottomAlign = 0;

    branchVar1 = empty;
    branchVar2 = empty;

    labelDigits1 = 6;
    labelDigits2 = 6;

    windowTextMarginH = HY_TREEPANEL_MARGIN;
    windowTextMarginV = HY_TREEPANEL_MARGIN;

    labelFont.style = HY_FONT_PLAIN;
    l1->SetFont (labelFont);
    l1->SetBackColor (col);
    l2->SetFont (labelFont);
    l2->SetBackColor (col);
    l3->SetFont (labelFont);
    l3->SetBackColor (col);
    l4->SetBackColor (col);
    bb->SetBackColor (col);
    canvasSettings = l2->_SuggestDimensions();
    canvasSettings.width = HY_COMPONENT_NO_SCROLL;
    canvasSettings.top = canvasSettings.bottom = 30;
    canvasSettings.right+=10;
    canvasSettings.left+=10;
    col.R = 254;
    col.G = 242;
    col.B = 208;
    l1->SetForeColor  (col);
    l2->SetForeColor  (col);
    l3->SetForeColor  (col);
    l1->SetDimensions (canvasSettings,canvasSettings);
    l2->SetDimensions (canvasSettings,canvasSettings);
    canvasSettings.top = canvasSettings.bottom = 50;
    l3->SetDimensions (canvasSettings,canvasSettings);
    SetWindowRectangle (0,0,HY_TREEPANEL_DEFSIZE,HY_TREEPANEL_DEFSIZE);

    treeName = title;
    hSpace = HY_TREEPANEL_HSPACE;
    vSpace = HY_TREEPANEL_VSPACE;
    textSpace = HY_TREEPANEL_LABELSPACING;
    treeFlags = HY_TREEPANEL_SCALED|HY_TREEPANEL_SQUARE|HY_TREEPANEL_TIP_LABELS;
    if (argument.sLength) {
        if (argument.sData[0]=='(') {
            SetTreeString (argument);
        } else {
            SetVariableReference (argument);
        }
    }
    distortion         = 1.5;
    arcStart           = 0;
    arcEnd             = 6.284;
    saveMouseH         = 0;
    saveMouseV         = 0;
    toolTipBounds.left = 0;
    lastSave           = 0;

}

//__________________________________________________________
void _HYTreePanel::dumpCoordTree (void)
{
    if (coordTree) {
        coordTree->delete_tree();
        delete coordTree;
        coordTree = nil;
    }
}

//__________________________________________________________
_HYTreePanel::~_HYTreePanel()
{
    dumpCoordTree();
    FlushUndoData();
}

//__________________________________________________________

bool    _HYTreePanel::ProcessGEvent (_HYEvent* e)
{
    _String firstArg;
    long    k,f;

    bool    done = false;

    if (e->EventClass()==_hyGlobalLFKillEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        if (k!=GetID()) {
            firstArg        = e->EventCode().Cut (f+1,-1);
            if (likeFuncID==firstArg.toNum()) {
                likeFuncID = -1;
                firstArg = statusBar.Cut(0,statusBar.FindBackwards (". Ln-likelihood",0,-1)-1);
                SetStatusBar (firstArg);
                _PaintLFStatus();
                if (scaleVariable.sLength&&(!scaleVariable.Equal(&assValScale))) {
                    scaleVariable = assValScale;
                    UpdateScalingVariablesList ();
                    BuildTree (true);
                    RenderTree(true);
                }
                done = true;
            }
        }
    } else if (e->EventClass()==_hyGlobalTreeKillEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        if (k!=GetID()) {
            firstArg = e->EventCode().Cut (f+1,-1);
            k = LocateVarByName(treeName);
            if (k>=0)
#ifndef USE_AVL_NAMES
                if (variableReindex.lData[k] == firstArg.toNum())
#else
                if (variableNames.GetXtra (k) == firstArg.toNum())
#endif
                {
                    postWindowCloseEvent (GetID());
                    done = true;
                }
        }
    } else if (e->EventClass()==_hyGlobalChangeLF) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        if (k!=GetID()) {
            firstArg        = e->EventCode().Cut (f+1,-1);
            if (likeFuncID==firstArg.toNum()) {
                if (scaleVariable.sLength&&(!scaleVariable.Equal(&assValScale))) {
                    UpdateScalingVariablesList ();
                    BuildTree (true);
                    RenderTree(true);
                }
                done = true;
            }
        }
    } else if (e->EventClass()==_hyGlobalLFSpawnEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        if (k!=GetID()) {
            if (likeFuncID == -1) {
                k       = e->EventCode().Cut (f+1,-1).toNum();

                _LikelihoodFunction *lf = (_LikelihoodFunction*)likeFuncList (k);
                if (lf->DependOnTree (treeName)>=0) {
                    likeFuncID = k;
                    UpdateScalingVariablesList ();
                    scaleVariable = expectedNumberOfSubs;
                    BuildTree (true);
                    RenderTree(true);
                }
            }
        }
    } else if (e->EventClass()==_hyGlobalSetTreePanelSelection) {
        _TheTree * me = LocateMyTreeVariable ();
        if (me) {
            firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
            k = firstArg.toNum();
            if (k==me->GetAVariable()) {
                _List * nodeNames = e->EventCode().Tokenize(",");
                if (nodeNames->lLength>1) {
                    nodeNames->Delete(0);
                    SelectRangeAndScroll(*nodeNames);
                    done = true;
                }
                DeleteObject (nodeNames);
            }
        }
    }

    if (!done) {
        return _HYWindow::ProcessGEvent (e);
    }
    return true;
}

//__________________________________________________________

void  _HYTreePanel::HandleContextPopup (long h, long v)
{
    _List       menuOptions;
    _String     menuChoice ("Copy as Picture");

    menuOptions && & menuChoice;

    menuChoice = "Copy as Newick string";
    menuOptions && & menuChoice;

    if (currentSelection.lLength) {
        menuChoice = "Count Nodes in Selection";
        menuOptions && & menuChoice;
        menuChoice = "Copy Names in Selection";
        menuOptions && & menuChoice;
    }

    menuChoice = HandlePullDown (menuOptions, h,v ,-1);

    switch (menuOptions.FindObject(&menuChoice)) {
    case 0: {
        _HYCanvas* theCanvas = (_HYCanvas*)GetObject (0);
        theCanvas->_CopyToClipboard ();
        break;
    }
    case 1: {
        menuChoice = GetTreeString();
        PlaceStringInClipboard (menuChoice, GetOSWindowData());
        break;
    }
    case 2: {
        long             lCount = 0;
        long             iCount = 0;
        for (long k=0; k<currentSelection.lLength; k++) {
            node<nodeCoord>* nn = (node<nodeCoord>*)currentSelection (k);
            if (nn->get_num_nodes()) {
                iCount++;
            } else {
                lCount++;
            }
        }
#ifdef __MAC__
        _String newLine ("\r");
#else
        _String newLine ("\n\r");
#endif
        _String msg = _String("Current selection includes") & newLine &_String(lCount) & " leaves and" & newLine& _String(iCount) & " internal nodes.";
        ProblemReport (msg, (Ptr)this);
        break;
    }
    case 3: {
        _String   res (currentSelection.lLength*16,true);

        if (treeFlags&HY_TREEPANEL_TIP_LABELS)
            for (long k=0; k<currentSelection.lLength; k++) {
                node<nodeCoord>* nn = (node<nodeCoord>*)currentSelection (k);
                if (!nn->get_num_nodes()) {
                    res << nn->in_object.branchName;
                    res << '\n';
                }
            }

        if (treeFlags&HY_TREEPANEL_INT_LABELS)
            for (long k=0; k<currentSelection.lLength; k++) {
                node<nodeCoord>* nn = (node<nodeCoord>*)currentSelection (k);
                if (nn->get_num_nodes()) {
                    res << nn->in_object.branchName;
                    res << '\n';
                }
            }
        res.Finalize();
        PlaceStringInClipboard (res, GetOSWindowData());
        break;
    }
    }
}

//__________________________________________________________

void  _HYTreePanel::DisplayParameterTable (void)
{
    if (likeFuncID>=0) {
        _String     windowName (*(_String*)likeFuncNamesList (likeFuncID));
        if (windowName.endswith ("_LF")) {
            windowName.Trim (0,windowName.sLength-4);
        }
        windowName = _String ("Likelihood parameters for ") & windowName;
        long    k = FindWindowByName (windowName);
        _HYParameterTable*    thePT;
        if (k>=0) {
            _HYPlatformWindow* thisWindow = (_HYPlatformWindow*)windowObjects(k);
            thisWindow->_Activate();
            thePT = (_HYParameterTable*)windowObjectRefs(k);
        } else {
            thePT = new _HYParameterTable (windowName,likeFuncID);
            thePT->BringToFront();
        }
        if (currentSelection.lLength) {
            _List   nodeNames;
            for (k=0; k<currentSelection.lLength; k++) {
                node<nodeCoord>* node1 = (node<nodeCoord>*)currentSelection(k);
                _String          nodeName = treeName&'.'&node1->in_object.branchName&'.';
                nodeNames.BinaryInsert(&nodeName);
            }
            thePT->SetSelectedRows (nodeNames);
        }
    }
}

//__________________________________________________________

void  _HYTreePanel::HandleSelection (char menuChoice)
{
    switch (menuChoice) {
    case 0: { // Select subtree
        SelectEntireSubtree();
        break;
    }

    case 1: { // Show incomplete
        SelectIncompleteBranches();
        break;
    }

    case 2: { // Show w/o models
        SelectBranchesWithoutModel();
        break;
    }

    case 3: { // branches by name
        SelectBranchesByName();
        break;
    }

    case 4: { // branches by name
        SelectBranchesByLength();
        break;
    }

    case 5: { // invert selection
        InvertSelection();
        break;
    }

    case 6: { // grow selection
        GrowSelection();
        break;
    }

    case 7: { // find selection in data panel
        HighlightSelectionInDataPanel();
        break;
    }

    case 8: { // find selection in data panel
        ShowSelectionInAnotherTree ();
        break;
    }
    }
}

//__________________________________________________________

void  _HYTreePanel::HandleSearchAndReplace (bool selectionOnly)
{
    _TheTree* meTree = LocateMyTreeVariable();

    if (meTree) {
        _String p1 ("Search for:"),
                p2 ("Replace with:");

        if (EnterString2Dialog (sandrString1, sandrString2, p1, p2, (Ptr)this)) {
            _List   oldNames,
                    newNames;

            long    cutAt = meTree->GetName()->sLength+1;

            node<nodeCoord>* currentNd = NodeTraverser(coordTree);
            while (currentNd->parent) {
                if (currentNd->in_object.branchName.sLength) {
                    if ((selectionOnly && (currentNd->in_object.flags|=HY_BRANCH_SELECT)) || (!selectionOnly)) {
                        _String nodeName  = LocateVar(currentNd->in_object.varRef)->GetName()->Cut(cutAt,-1),
                                repString = nodeName.Replace (sandrString1, sandrString2, true);

                        if ((!repString.Equal (&nodeName)) && repString.IsValidIdentifier()) {
                            oldNames && & nodeName;
                            newNames && & repString;
                        }
                    }
                }
                currentNd = NodeTraverser((node<nodeCoord>*)nil)    ;
            }


            if (oldNames.lLength) {
                _HYDataPanel* thisDP = nil;

                for (long dpID=0; dpID<windowObjects.lLength; dpID++) {
                    _HYWindow * thisWindow = (_HYWindow*)windowObjectRefs(dpID);

                    if (thisWindow->WindowKind() == HY_WINDOW_KIND_DATAPANEL && ((_HYDataPanel*)thisWindow)->DependOnTree (meTree->GetAVariable())) {
                        thisDP = (_HYDataPanel*) thisWindow;
                        break;
                    }
                }

                if (thisDP) {
                    if (updateInDataPanel == false) {
                        _String rid = _String("This tree is attached to data panel ") & thisDP->GetTitle() & ". Would you like to change sequence names in that window as well?";
                        if (!ProceedPromptWithCheck (rid,donotWarnAgain,updateInDataPanel,(Ptr)this)) {
                            thisDP = nil;
                        }
                    }
                }

                if (thisDP) {
                    ((_HYSequencePane*)thisDP->GetObject (0))->BatchRenameSequences (oldNames,newNames);
                }

                BufferToConsole ("\n");
                _String prefix = *meTree->GetName() & '.';

                for (long k=0; k<oldNames.lLength; k++) {
                    p1 =   *(_String*)oldNames(k);
                    p2 =   *(_String*)newNames(k);

                    StringToConsole (p1);
                    BufferToConsole (" => ");
                    StringToConsole (p2);
                    BufferToConsole ("\n");

                    p1 = prefix & p1;
                    p2 = prefix & p2;

                    RenameVariable (&p1,&p2);
                }

                PrepareUndoData (-1);
                _SimpleList     saveSel;
                PreserveSelection (saveSel);
                BuildTree(true);
                RestoreSelection (saveSel);
                if (likeFuncID>=0) {
                    postChangeLFEvent(GetID(),likeFuncID);
                }
                RenderTree ();
            }
        }
    }
}

//__________________________________________________________

void  _HYTreePanel::HandleComparison (char menuChoice)
{
    _TheTree* meTree = LocateMyTreeVariable(),
              * me = nil;

    if (meTree) {
        _List           choices;

        _SimpleList     all,
                        selectors;

        selectors << 0;
        selectors << 1;

#ifndef USE_AVL_NAMES
        for (long vi = 0; vi < variableNames.lLength; vi++)
#else
        _SimpleList tcache;
        long        ti,
                    vi = variableNames.Traverser (tcache, ti, variableNames.GetRoot());

        for (; vi >=0 ; vi = variableNames.Traverser (tcache, ti))
#endif
        {
            _Variable* vv = (_Variable*)FetchVar(vi);
            if ((vv->ObjectClass () == TREE) && (vv!=meTree)) {
                _List   aChoice;
                aChoice << vv->GetName();
                aChoice << vv->GetName();
                choices && & aChoice;
                all << all.lLength;
            }
        }

        if (choices.lLength) {
            long choice = HandleListSelection (choices,selectors, all, "Target for comparison", all,1);
            if (choice>=0) {
                me = (_TheTree*)FetchVar(LocateVarByName(*(_String*)(*(_List*)choices(choice))(0)));
            }
        } else {
            _String treeError ("No other tree variables are presently resident in memory. Nothing to compare to!");
            ProblemReport (treeError);
        }
        if (me) {
            switch (menuChoice) {
            case 0: { // compare equal
                _String res = meTree->CompareTrees (me);
                ProblemReport (res, (Ptr)this);
                break;
            }

            case 1: { // find subtree
                node <long>* n1 = nil;
                _CalcNode* travNode = meTree->DepthWiseTraversal (true);
                while (travNode) {
                    if (travNode->theIndex==selectionTop->in_object.varRef) {
                        n1 = &meTree->GetCurrentNode();
                        break;
                    }
                    travNode = meTree->DepthWiseTraversal();
                }
                if (n1) {
                    _String res = me->CompareSubTrees (meTree,n1);
                    if (!res.beginswith("No")) {
                        _String cres (res,res.FirstSpaceIndex(0,-1,-1)+1,res.sLength-2);
                        postSetTreePanelSelection (me->GetAVariable(), &cres);
                    }
                    ProblemReport (res, (Ptr)this);
                }
                break;
            }

            case 2: { // Show max common subtree
                long    sizeVar = 0;
                _String res = meTree->FindMaxCommonSubTree (me, sizeVar, nil);

                if (!res.sLength) {
                    res = "No common subtrees";
                    ProblemReport (res, (Ptr)this);
                } else {
                    _List sel;
                    sel && & res;
                    SelectRangeAndScroll (sel);
                    SelectEntireSubtree  ();
                    res = _String ("Maximal subtree found and selected. The subtree contains ") & sizeVar & " leaves.";
                    ProblemReport (res, (Ptr)this);
                }

                break;
            }

            case 3: { // max forest
                long        sizeVar = 0;
                _List       forest;

                _String res = meTree->FindMaxCommonSubTree (me, sizeVar, &forest);

                if (forest.lLength==0) {
                    res = "No common subtrees";
                    ProblemReport (res, (Ptr)this);
                } else {
                    forest.Sort();
                    SelectRangeAndScroll (forest);
                    res = _String ("Maximal forest found and selected. The forest contains ") & sizeVar & " leaves.";
                    GrowSelection();
                    ProblemReport (res, (Ptr)this);
                }

                break;
            }
            case 4: { // match to pattern
                _String res = me->MatchTreePattern (meTree);
                ProblemReport (res, (Ptr)this);
                break;
            }
            }
        }
    }
}

//__________________________________________________________

void  _HYTreePanel::SelectBranchesByName   (void)
{
    _String             res,
                        bName  ("_branch_name"),
                        prompt ("Template(in terms of _branch_name):");

    if (EnterStringDialog   (selectByNameSave,prompt,(Ptr)this)) {
        res = selectByNameSave;
        _Variable*       bv = CheckReceptacle (&bName, empty, false);

        node<nodeCoord>* currentNd = NodeTraverser(coordTree),
                         * parent = nil;

        _FString         bn (empty);
        bv->SetValue    (&bn);
        _Formula        tempF (res);

        if (terminateExecution) {
            terminateExecution = false;
            return;
        }

        _PMathObj       fv = tempF.Compute();

        if (terminateExecution||!fv||!(fv->ObjectClass()==NUMBER)) {
            terminateExecution = false;
            return;
        }

        currentSelection.Clear();

        _Parameter  ch = 0,
                    cv = 0;

        while (currentNd) {
            parent = currentNd->parent;
            if (!parent) {
                break;
            }

            currentNd->in_object.flags&=HY_BRANCH_DESELECT;

            if (currentNd->in_object.varRef >= 0) {
                _String*   thisName = LocateVar(currentNd->in_object.varRef)->GetName();
                bName = thisName->Cut (thisName->FindBackwards (".",0,-1)+1,-1);

                bn.theString->Duplicate(&bName);
                bv->SetValue (&bn);
                fv = tempF.Compute();
                if (fv&&(fv->Value()>0.0)) {
                    currentNd->in_object.flags|=HY_BRANCH_SELECT;
                    currentSelection<<(long)currentNd;
                    ch += currentNd->in_object.h;
                    cv += currentNd->in_object.v;
                }
            }
            currentNd = NodeTraverser((node <nodeCoord>*)nil)   ;
        }

        RenderTree();
        ScrollToSelection (ch, cv);
    }
}

//__________________________________________________________

_String  dsSelectorPopulator (_HYWindow *thisWindow)
{
    _DataSet *ds = (_DataSet*)dataSetList (((_HYDataPanel*)thisWindow)->GetDSID());
    return (_String ("Data Set with ") & ds->NoOfSpecies() & " and " & ds->NoOfColumns() & " alignment columns.");
}

//__________________________________________________________

_String  treeSelectorPopulator (_HYWindow *thisWindow)
{
    _TheTree * meTree = ((_HYTreePanel*)thisWindow)->LocateMyTreeVariable();

    _String res;
    if (meTree) {
        _PMathObj tc = meTree->TipCount(),
                  bc = meTree->BranchCount ();
        res = (_String ("A tree with ") & tc->Value() & " tips and" & bc->Value() & " internal branches.");
        DeleteObject (tc);
        DeleteObject (bc);
    }
    return res;
}


//__________________________________________________________

long SelectOpenWindowObjects (long objectKind, _String objectDesc, listDescPopulator* popFunction, Ptr caller)
{
    _SimpleList  allowedChoices;

    _List           choices;

    _SimpleList     all,
                    selectors;

    selectors << 0;
    selectors << 1;

    for (long k=0; k<windowObjects.lLength; k++) {
        _HYWindow* thisWindow = (_HYWindow*)windowObjectRefs(k);
        if (caller == (Ptr)thisWindow) {
            continue;
        }

        if (thisWindow->WindowKind() == objectKind) {
            _String  *dsName = &thisWindow->GetTitle();
            if (dsName->sLength) {
                _List   aChoice;
                _String tt = (*popFunction)(thisWindow);
                if (tt.sLength) {
                    aChoice << dsName;
                    aChoice && & tt;
                    choices && & aChoice;
                    allowedChoices << k;
                    all << all.lLength;
                }
            }
        }
    }

    if (allowedChoices.lLength == 0) {
        _String       errMsg ("No open ");
        errMsg = errMsg & objectDesc & "s were found.";
        ProblemReport (errMsg, caller);
        return -1;
    }

    long choice = HandleListSelection (choices,selectors, all, _String("Choose a ") & objectDesc, all,1);
    if (choice >= 0) {
        choice = allowedChoices.lData[choice];
    }

    return choice;

}

//__________________________________________________________

void  _HYTreePanel::HighlightSelectionInDataPanel   (void)
{
    long choice = SelectOpenWindowObjects (HY_WINDOW_KIND_DATAPANEL, "data panel", dsSelectorPopulator, (Ptr)this);

    if (choice>=0) {
        _HYDataPanel * theDP = (_HYDataPanel*)windowObjectRefs(choice);

        _List          selectedNodes;

        for (long k=0; k<currentSelection.lLength; k++) {
            node<nodeCoord>* nn = (node<nodeCoord>*)currentSelection (k);
            if (nn->in_object.varRef) {
                _String upn = nn->in_object.branchName;
                upn.UpCase();
                selectedNodes && & upn;
            }
        }

        selectedNodes.Sort();

        ((_HYSequencePane*)theDP->GetObject (0))->SelectSequenceNames(selectedNodes);
    }
}

//__________________________________________________________

void  _HYTreePanel::ShowSelectionInAnotherTree   (void)
{
    long choice = SelectOpenWindowObjects (HY_WINDOW_KIND_TREE, "tree panel", treeSelectorPopulator, (Ptr)this);

    if (choice>=0) {
        _HYTreePanel * tp        = ((_HYTreePanel*)windowObjectRefs(choice));
        _TheTree     * otherTree = tp->LocateMyTreeVariable();

        _String        prefix = *otherTree->GetName() & '.';
        _List          selectedNodes;
        for (long k=0; k<currentSelection.lLength; k++) {
            node<nodeCoord>* nn = (node<nodeCoord>*)currentSelection (k);
            if (nn->in_object.varRef && nn->get_num_nodes() == 0) {
                _String upn = prefix&nn->in_object.branchName;
                upn.UpCase();
                selectedNodes && & upn;

            }
        }

        selectedNodes.Sort();
        tp->SelectRangeAndScroll (selectedNodes, true);

        _String errMsg = _String("Selected ") & (long)tp->currentSelection.lLength & '/' & (long)selectedNodes.lLength & " tips";
        ProblemReport (errMsg, (Ptr)this);
    }
}


//__________________________________________________________

void  _HYTreePanel::MatchToDataSet (void)
{
    _TheTree* meTree = LocateMyTreeVariable();

    if (meTree) {
        _Constant* tc = (_Constant*) meTree->TipCount ();
        long       countLeaves = tc->Value();
        DeleteObject (tc);

        _SimpleList  allowedChoices;

        _List           choices;

        _SimpleList     all,
                        selectors;

        selectors << 0;
        selectors << 1;

        for (long k=0; k<dataSetNamesList.lLength; k++) {
            _String* dsName = (_String*)dataSetNamesList(k);
            if (dsName->sLength) {
                _DataSet* ds = (_DataSet*)dataSetList (k);
                if (ds->NoOfSpecies() >= countLeaves) {
                    _List   aChoice;
                    aChoice << dsName;
                    _String tt = (_String ("Data Set with ") & ds->NoOfSpecies() & " and " & ds->NoOfColumns() & " alignment columns.");
                    aChoice && & tt;
                    choices && & aChoice;
                    allowedChoices << k;
                    all << all.lLength;
                }
            }
        }

        if (allowedChoices.lLength == 0) {
            _String       errMsg ("No sufficiently large data sets were found in memory.");
            ProblemReport (errMsg, (Ptr)this);
            return;
        }

        long choice = HandleListSelection (choices,selectors, all, "Target for comparison", all,1);
        if (choice>=0) {
            _DataSet* ds = (_DataSet*)dataSetList (allowedChoices.lData[choice]);

            _SimpleList     tipMatches;
            _CalcNode*      travNode = meTree->StepWiseTraversal(true);
            _List           tips;
            long            j,k;

            while (travNode) {
                if (meTree->IsCurrentNodeATip()) {
                    _String tipName (*travNode->GetName(),travNode->GetName()->FindBackwards('.',0,-1)+1,-1);
                    tips&& &tipName;
                }
                travNode = meTree->StepWiseTraversal(false);
            }

            _List  * sortedDSNames = (_List*)ds->GetNames().makeDynamic();
            sortedDSNames->Sort();

            for (j=0; j<tips.lLength; j++)
                if ((k = sortedDSNames->BinaryFindObject(tips(j))) < 0) {
                    break;
                } else {
                    tipMatches<<k;
                }

            DeleteObject (sortedDSNames);

            if (j<tips.lLength) {
                long sj = j;
                for (j=0; j<tips.lLength; j++) {
                    _String *thisName = (_String*)tips(j);
                    k = atoi (thisName->sData);
                    _String tryAgain (k);
                    if (tryAgain.Equal(thisName) && k<=tips.lLength) {
                        tipMatches<<k;
                    } else {
                        break;
                    }
                }

                if (j==tips.lLength) {
                    if (tipMatches.Find(0)==-1)
                        for (j=0; j<tips.lLength; j++) {
                            tipMatches.lData[j]--;
                        }
                } else {
                    j=sj;
                }


                if (j != tips.lLength) {
                    _String  errMsg ("Match failed. First offending leaf is labeled '");
                    errMsg = errMsg & *(_String*)tips (j) & "'.";
                    ProblemReport (errMsg, (Ptr)this);

                    _List    unmatchedNode;
                    errMsg = *meTree->GetName () & '.' & *(_String*)tips (j);
                    unmatchedNode && & errMsg;
                    SelectRangeAndScroll (unmatchedNode);
                } else {
                    _String  errMsg ("Matched by index. Would you like to relabel the tree with names from the data? Be careful with this option if the tree is currently a part of a likelihood function, because some batch language instructions may fail then.");
                    if (ProceedPrompt (errMsg, (Ptr)this)) {
                        _CalcNode*      travNode = meTree->StepWiseTraversal(true);
                        j = 0;
                        while (travNode) {
                            if (meTree->IsCurrentNodeATip()) {
                                _String newName = *meTree->GetName() & '.' & *(_String*)(ds->GetNames())(tipMatches.lData[j++]);
                                RenameVariable (travNode->GetName(), &newName);

                            }
                            travNode = meTree->StepWiseTraversal(false);
                        }

                        BuildTree (true);
                        RenderTree();
                    }

                }
            } else {
                _String       errMsg ("Successfully matched.");
                ProblemReport (errMsg, (Ptr)this);
            }
        }
    }
}

//__________________________________________________________

void  _HYTreePanel::SelectBranchesByLength (void)
{
    _String             res,
                        bName  ("_branch_length"),
                        prompt ("Template(in terms of _branch_length):"),
                        cprt   ("Restrict to internal branches");


    if (EnterStringDialogWithCheckbox   (selectByBLSave,prompt,cprt,brLengthIntOnly,(Ptr)this)) {
        res = selectByBLSave;
        _Variable*       bv = CheckReceptacle (&bName, empty, false);

        node<nodeCoord>* currentNd = NodeTraverser(coordTree),
                         * parent = nil;

        _Constant        bn (0.0);
        bv->SetValue    (&bn);
        _Formula        tempF (res);

        if (terminateExecution) {
            terminateExecution = false;
            return;
        }

        _PMathObj       fv = tempF.Compute();

        if (terminateExecution||!fv||!(fv->ObjectClass()==NUMBER)) {
            terminateExecution = false;
            return;
        }

        currentSelection.Clear();

        _Parameter  ch = 0,
                    cv = 0;

        while (currentNd) {
            parent = currentNd->parent;
            if (!parent) {
                break;
            }

            currentNd->in_object.flags&=HY_BRANCH_DESELECT;

            if (currentNd->in_object.varRef >= 0) {
                if (!brLengthIntOnly || currentNd->get_num_nodes()) {
                    bn.SetValue(currentNd->in_object.bL);
                    bv->SetValue (&bn);
                    fv = tempF.Compute();
                    if (fv&&(fv->Value()>0.0)) {
                        currentNd->in_object.flags|=HY_BRANCH_SELECT;
                        currentSelection<<(long)currentNd;
                        ch += currentNd->in_object.h;
                        cv += currentNd->in_object.v;

                        StringToConsole(currentNd->in_object.branchName);
                        BufferToConsole("\n");
                    }
                }
            }
            currentNd = NodeTraverser((node <nodeCoord>*)nil)   ;
        }

        RenderTree();
        ScrollToSelection (ch, cv);
    }

}


//__________________________________________________________

bool    _HYTreePanel::ProcessEvent (_HYEvent* e)
{
    long    f,
            i,
            k;

    _String firstArg;
    if (e->EventClass()==_hyMenuSelChangeEvent)
        // menu change event here
    {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);
        k = (e->EventCode().Cut (f+1,-1)).toNum();
        if (i==3) { //layout mode
            unsigned int nflags = treeFlags&HY_TREEPANEL_LAYOUT_MASK;
            switch (k) {
            case 0:
                nflags|=HY_TREEPANEL_SQUARE;
                break;
            case 1:
                nflags|=HY_TREEPANEL_STRAIGHT;
                break;
            case 2:
                nflags|=HY_TREEPANEL_ARCS;
                break;
            case 3:
                nflags|=HY_TREEPANEL_CIRCULAR;
                break;
            }
            SetFlags (nflags);
        } else if (i==5) { // scaling var
            _HYPullDown*  selVar = (_HYPullDown*)GetObject(5);
            // k == 0 -> unscaled
            // k > 0 -> scaled
            treeFlags    = treeFlags&HY_TREEPANEL_SCALING_MASK;
            treeFlags   |= (k?HY_TREEPANEL_SCALED:HY_TREEPANEL_ALIGNED);
            switch (k) {
            case 0:
                scaleVariable = empty;
                break;
            case 1:
                scaleVariable = expectedNumberOfSubs;
                break;
            case 2:
                scaleVariable = stringSuppliedLengths;
                break;
            default:
                scaleVariable = *selVar->GetMenuItem(k);
            }
            if (BuildTree(true)) {
                RenderTree((treeFlags&HY_TREEPANEL_PROJECTION)==0);
            }
        }
        DeleteObject (e);
        return true;
    } else {
        if (e->EventClass()==_hyButtonPushEvent) {
            firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
            i = MatchComponentID (firstArg);
            k = (e->EventCode().Cut (f+1,-1)).toNum();
            if (i==7) { // zoom button bar
                _HYButtonBar* zoom = (_HYButtonBar*)GetObject (7);
                if (k==5) { // fit to window
                    FitToWindow();
                } else if (k==6) { // fisheye
                    HandleFisheyeButton ();
                } else if (k==4) { // change orientation
                    zoom->ReplaceButton(4,ProcureIconResource(HY_TREEPANEL_ICON_ID+5-IsVertical()));
                    SetVertical (!IsVertical());
                } else // expand or contract
                    if (IsVertical())
                        switch (k) {
                        case 0: // expand horizontally
                            f = vSpace+HY_TREEPANEL_SPACE_STEP;
                            zoom->EnableButton (0,f<HY_TREEPANEL_MAX_VSPACE);
                            zoom->EnableButton (2,vSpace<=HY_TREEPANEL_MIN_VSPACE);
                            SetVSpace(f);
                            break;
                        case 1: // expand vertically
                            f = hSpace+HY_TREEPANEL_SPACE_STEP;
                            if (f>=HY_TREEPANEL_MAX_HSPACE) {
                                zoom->EnableButton (1,false);
                            }
                            if (hSpace<=HY_TREEPANEL_MIN_HSPACE) {
                                zoom->EnableButton (3,true);
                            }
                            SetHSpace(f);
                            break;
                        case 2: // contract horizontally
                            f = vSpace-HY_TREEPANEL_SPACE_STEP;
                            if (f<=HY_TREEPANEL_MIN_VSPACE) {
                                zoom->EnableButton (2,false);
                            }
                            if (vSpace>=HY_TREEPANEL_MAX_VSPACE) {
                                zoom->EnableButton (0,true);
                            }
                            SetVSpace(f);
                            break;
                        case 3: // contract vertically
                            f = hSpace-HY_TREEPANEL_SPACE_STEP;
                            if (f<=HY_TREEPANEL_MIN_HSPACE) {
                                zoom->EnableButton (3,false);
                            }
                            if (hSpace>=HY_TREEPANEL_MAX_HSPACE) {
                                zoom->EnableButton (1,true);
                            }
                            SetHSpace(f);
                            break;
                        }
                    else
                        switch (k) {
                        case 0: // expand horizontally
                            f = hSpace+HY_TREEPANEL_SPACE_STEP;
                            if (f>=HY_TREEPANEL_MAX_HSPACE) {
                                zoom->EnableButton (0,false);
                            }
                            if (hSpace<=HY_TREEPANEL_MIN_HSPACE) {
                                zoom->EnableButton (2,true);
                            }
                            SetHSpace(f);
                            break;
                        case 1: // expand vertically
                            f = vSpace+HY_TREEPANEL_SPACE_STEP;
                            if (f>=HY_TREEPANEL_MAX_VSPACE) {
                                zoom->EnableButton (1,false);
                            }
                            if (vSpace<=HY_TREEPANEL_MIN_VSPACE) {
                                zoom->EnableButton (3,true);
                            }
                            SetVSpace(f);
                            break;
                        case 2: // contract horizontally
                            f = hSpace-HY_TREEPANEL_SPACE_STEP;
                            if (f<=HY_TREEPANEL_MIN_HSPACE) {
                                zoom->EnableButton (2,false);
                            }
                            if (hSpace>=HY_TREEPANEL_MAX_HSPACE) {
                                zoom->EnableButton (0,true);
                            }
                            SetHSpace(f);
                            break;
                        case 3: // contract vertically
                            f = vSpace-HY_TREEPANEL_SPACE_STEP;
                            if (f<=HY_TREEPANEL_MIN_VSPACE) {
                                zoom->EnableButton (3,false);
                            }
                            if (vSpace>=HY_TREEPANEL_MAX_VSPACE) {
                                zoom->EnableButton (1,true);
                            }
                            SetVSpace(f);
                            break;
                        }
            }
        } else if (e->EventClass()==_hyRebuildSCanvasEvent) {
            k = e->EventCode().toNum();
            for (i=0; i<components.lLength; i++)
                if (((_HYGuiObject*)components(i))->MatchID(k)) {
                    break;
                }
            if (i==2) {
                FitToWindow(true);
            }
        }

    }
    DeleteObject (e);
    return false;
}

//__________________________________________________________

void    _HYTreePanel::SetTreeString (_String& theString)
{
    // need to make the tree variable, eh
    _String treeStringPrep = _String ("Tree ") & treeName & " = " & theString;
    ExecuteBLString (treeStringPrep, nil);
    UpdateScalingVariablesList ();
}

//__________________________________________________________

void    _HYTreePanel::ClearClipboardContents (void)
{
    if (treePanelClipboardRoot) {
        KillNodesInClipboard (treePanelClipboardRoot);
        treePanelClipboardRoot = nil;
    }
}

//__________________________________________________________

void    _HYTreePanel::SetVariableReference (_String& varName)
{
    _TheTree * myTree = (_TheTree *)FetchObjectFromVariableByType (&varName, TREE);

    _HYPullDown* scalingVars = (_HYPullDown*)GetObject (5);
    long f = scalingVars->MenuItemCount();

    if (myTree) {
        treeName = varName;
        for (long k=3; k<f; k++) {
            scalingVars->DeleteMenuItem(3);
        }

        _SimpleList newScalingVariables;
        (myTree)->FindScalingVariables(newScalingVariables);
        scaleVariable = empty;
        likeFuncID = -1;

        for (long f=0; f<likeFuncList.lLength; f++) {
            if (((_String*)likeFuncNamesList(f))->sLength) {
                _LikelihoodFunction* theLF = (_LikelihoodFunction*)likeFuncList (f);
                if (theLF->DependOnTree (varName)>=0) {
                    if (likeFuncID>=0) {
                        likeFuncID=-1;
                        break;
                    } else {
                        likeFuncID = f;
                    }
                }
            }
        }
        UpdateScalingVariablesList ();
    } else
        for (long k=3; k<f; k++) {
            scalingVars->DeleteMenuItem(k);
        }
}

//__________________________________________________________

void    _HYTreePanel::UpdateScalingVariablesList (void)
{
    _TheTree *t = LocateMyTreeVariable ();
    if (t) {
        long k,f, oldSel;

        _List           uE;
        _SimpleList     nodeList, modelSCV;
        bool            esubs, blen;

        _HYPullDown* scalingVars = (_HYPullDown*)GetObject (5);

        f = scalingVars->MenuItemCount();

        oldSel = scalingVars->GetSelection();

        for (k=3; k<f; k++) {
            scalingVars->DeleteMenuItem(3);
        }

        esubs = t->FindScalingVariables(modelSCV);
        blen = t->HaveStringBranchLengths();
        scalingVars->EnableItem(1,esubs);
        scalingVars->EnableItem(2,blen);
        if (!esubs) {
            likeFuncID = -1;
        }

        _CalcNode       *travNode = t->StepWiseTraversal (true);

        travNode = t->StepWiseTraversal (false);
        while (travNode) {
            nodeList<<travNode->GetAVariable();
            travNode = t->StepWiseTraversal (false);
        }
        CompileListOfUserExpressions (nodeList,uE,true);

        for (k=0; k<uE.lLength; k++) {
            _String *userVar = (_String*)uE(k);
            if (userVar->sData[0] == '!') {
                userVar->Trim (1,-1);
            }
        }

        uE.Sort();

        for (k=1; k<uE.countitems(); k++) {
            if (((_String*)uE(k-1))->Equal((_String*)uE(k))) {
                uE.Delete(k);
                k--;
            }
        }

        for (k=0; k<modelSCV.lLength; k++) {
            _String modelVar = *LocateVar (modelSCV.lData[k])->GetName();
            if (uE.FindObject(&modelVar)<0) {
                uE.InsertElement (&modelVar,k,true);
            }
        }

        for (k=0; k<uE.lLength; k++) {
            scalingVars->AddMenuItem (*(_String*)uE(k),-1);
        }

        if (oldSel&&(oldSel<scalingVars->MenuItemCount())) {
            if (oldSel>2) {
                scalingVars->ChangeSelection (oldSel);
                return;
            }
            if (oldSel==2) {
                esubs = false;
            }
        }
        if (esubs) {
            scalingVars->ChangeSelection (1);
        } else if (blen) {
            scalingVars->ChangeSelection (2);
        } else {
            scalingVars->ChangeSelection (0);
        }

    }
}

//__________________________________________________________

bool    _HYTreePanel::SetHSpace (int newSpace, bool force)
{
    if ((newSpace!=hSpace)||force) {
        _HYCanvas* theCanvas = (_HYCanvas*)GetObject (0),
                   * theRuler  = (_HYCanvas*)GetObject (2);
        _HYRect    myDims = theCanvas->GetCanvasSize();

        _Parameter t;


        t = ((_Parameter)newSpace)/hSpace;
        hScale *= t;
        if (treeFlags&HY_TREEPANEL_CIRCULAR) {
            ShiftScreenCoordinates (-windowTextMarginH-((treeFlags&HY_TREEPANEL_TIP_LABELS)?textSpace:0),
                                    -windowTextMarginV,coordTree);
            Convert2ScreenCoordinates (t,1,0,coordTree);
            if (treeFlags&HY_TREEPANEL_TIP_LABELS) {
                ShiftScreenCoordinates (textSpace,0,coordTree);
            }
        } else {
            ShiftScreenCoordinates (-windowTextMarginH,-windowTextMarginV,coordTree);
            Convert2ScreenCoordinates (t,1,0,coordTree);
            //ShiftScreenCoordinates (0,-HY_TREEPANEL_MARGIN,coordTree);
        }
        hSpace = newSpace;

        bool       isv = IsVertical();
        t =        coordTree->tree_depth()-1.;
        int        visTreeWidth = t*hSpace+t,
                   visWidth = isv?theCanvas->hOrigin+componentB.lData[0]-componentT.lData[0]:
                              theCanvas->hOrigin+componentR.lData[0]-componentL.lData[0],tmp;

        visTreeWidth += 2*windowTextMarginH
                        +((treeFlags&HY_TREEPANEL_TIP_LABELS)?textSpace:0);

        if (visTreeWidth > HY_TREEPANEL_MAXSIZE) {
            visTreeWidth = HY_TREEPANEL_MAXSIZE;
        }
        bool change = false;

        tmp = isv?myDims.bottom:myDims.right;

        if (visTreeWidth>tmp) {
            change = true;
        } else {
            if (visWidth<tmp) {
                change = true;
                if (visTreeWidth<visWidth) {
                    visTreeWidth = visWidth;
                }
            }
        }
        if (change) {
            if (isv) {
                theCanvas->SetCanvasSize (visTreeWidth,myDims.right,TREE_PANEL_PIXEL_DEPTH);
            } else {
                theCanvas->SetCanvasSize (myDims.bottom,visTreeWidth,TREE_PANEL_PIXEL_DEPTH);
            }
            myDims.left = componentL.lData[0];
            myDims.top = componentT.lData[0];
            myDims.right = componentR.lData[0];
            myDims.bottom = componentB.lData[0];
            theCanvas->SetVisibleSize (myDims);

            if (!isv) {
                theRuler->SetCanvasSize(HY_TREEPANEL_RULER_EXPANDED,visTreeWidth,TREE_PANEL_PIXEL_DEPTH);

                myDims.left = componentL.lData[2];
                myDims.top = componentT.lData[2];
                myDims.right = componentR.lData[2];
                myDims.bottom = componentB.lData[2];
                theRuler->SetVisibleSize (myDims);
            }

#ifdef __HYPHY_GTK__
            UpdateComponentInfo ();
#else
            dim = MinMaxWindowDimensions();
#endif
        }
        RenderTree();
        //RenderRuler();
        return true;
    }
    return false;
}

//__________________________________________________________

bool    _HYTreePanel::SetVSpace (int newSpace, bool force, bool render)
{
    if ((newSpace!=vSpace)||force) {
        _Parameter t = ((_Parameter)newSpace)/vSpace;
        vScale *= t;
        ShiftScreenCoordinates (0,-windowTextMarginV,coordTree);
        Convert2ScreenCoordinates (1,t,0,coordTree);
        ShiftScreenCoordinates (-windowTextMarginH,0,coordTree);
        vSpace = newSpace;
        _HYCanvas* theCanvas = (_HYCanvas*)GetObject (0),
                   * theRuler  = (_HYCanvas*)GetObject (2);
        _HYRect    myDims = theCanvas->GetCanvasSize();
        bool       isv = IsVertical();
        int visTreeWidth = tips*vSpace+2*windowTextMarginV,tmp,
            visWidth = isv?componentR.lData[0]+theCanvas->vOrigin-componentL.lData[0]:componentB.lData[0]+theCanvas->vOrigin-componentT.lData[0];
        if (visTreeWidth > HY_TREEPANEL_MAXSIZE) {
            visTreeWidth = HY_TREEPANEL_MAXSIZE;
        }
        bool change = false;
        tmp = isv?myDims.right:myDims.bottom;
        if (visTreeWidth>tmp) {
            change = true;
        } else {
            if (visWidth<tmp) {
                change = true;
                if (visTreeWidth<visWidth) {
                    visTreeWidth = visWidth;
                }
            }
        }
        if (change) {
            if (isv) {
                theCanvas->SetCanvasSize (myDims.bottom,visTreeWidth,TREE_PANEL_PIXEL_DEPTH);
                theRuler->SetCanvasSize (HY_TREEPANEL_RULER_EXPANDED,visTreeWidth,TREE_PANEL_PIXEL_DEPTH);
                myDims.left = componentL.lData[2];
                myDims.top = componentT.lData[2];
                myDims.right = componentR.lData[2];
                myDims.bottom = componentB.lData[2];
                theRuler->SetVisibleSize (myDims);
            } else {
                theCanvas->SetCanvasSize (visTreeWidth,myDims.right,TREE_PANEL_PIXEL_DEPTH);
            }
            myDims.left = componentL.lData[0];
            myDims.top = componentT.lData[0];
            myDims.right = componentR.lData[0];
            myDims.bottom = componentB.lData[0];
            theCanvas->SetVisibleSize (myDims);

#ifdef __HYPHY_GTK__
            UpdateComponentInfo ();
#else
            dim = MinMaxWindowDimensions();
#endif
        }
        if (render) {
            RenderTree();
            return true;
        }
    }
    return false;
}
//__________________________________________________________

long    _HYTreePanel::GetUniversalSaveOptions (_List&l)
{
    _String option ("Newick String");
    l&& & option;
    option = "NEXUS";
    l&& & option;
    option = "NEXUS Numeric";
    l&& & option;
    option = "PostScript file [cladogram]";
    l&& & option;
    option = "PostScript file [radial]";
    l&& & option;
    return 5;
}

//__________________________________________________________

void    _HYTreePanel::FitToWindow (bool force)
{
    _TheTree* me = LocateMyTreeVariable();
    if (me) {
        _PMathObj   leafCount = me->TipCount();
        long newHSpace,
             newVSpace;

        if (treeFlags&HY_TREEPANEL_CIRCULAR) {
            _Parameter maxHCoord = 0,
                       maxVCoord = 0;

            node<nodeCoord>* currentNd = NodeTraverser(coordTree);
            while (currentNd) {
                if (currentNd->in_object.h>maxHCoord) {
                    maxHCoord = currentNd->in_object.h;
                }

                if (currentNd->in_object.v>maxVCoord) {
                    maxVCoord = currentNd->in_object.v;
                }

                currentNd = NodeTraverser((node<nodeCoord>*)NULL);
            }
            /*if (IsVertical())
            {
                newHSpace = maxVCoord;
                maxVCoord = maxHCoord;
                maxHCoord = newHSpace;
            }*/

            if (IsVertical()) {
                newHSpace = (componentB.lData[0]-componentT.lData[0]-(2*windowTextMarginV
                             +((treeFlags&HY_TREEPANEL_TIP_LABELS)?textSpace:0)))/maxHCoord*hSpace;
                newVSpace = (componentR.lData[0]-componentL.lData[0]-2*windowTextMarginH)/maxVCoord*vSpace;
            } else {
                newHSpace = (componentR.lData[0]-componentL.lData[0]-(2*windowTextMarginH
                             +((treeFlags&HY_TREEPANEL_TIP_LABELS)?textSpace:0)))/maxHCoord*hSpace;
                newVSpace = (componentB.lData[0]-componentT.lData[0]-2*windowTextMarginV)
                            /maxVCoord*vSpace;
            }
        } else {
            if (IsVertical()) {
                newHSpace = (componentB.lData[0]-componentT.lData[0]-(2*windowTextMarginV
                             +((treeFlags&HY_TREEPANEL_TIP_LABELS)?textSpace:0)))/coordTree->tree_depth();
                newVSpace = (componentR.lData[0]-componentL.lData[0]-2*windowTextMarginH)/
                            leafCount->Value();
            } else {
                newHSpace = (componentR.lData[0]-componentL.lData[0]-(2*windowTextMarginH
                             +((treeFlags&HY_TREEPANEL_TIP_LABELS)?textSpace:0)))/coordTree->tree_depth();
                newVSpace = (componentB.lData[0]-componentT.lData[0]-2*windowTextMarginV)/
                            leafCount->Value();
            }
        }
        DeleteObject (leafCount);


        if (force || (newHSpace!=hSpace) || (newVSpace!=vSpace)) {
            ((_HYCanvas*)GetObject (0))->hOrigin = 0;
            ((_HYCanvas*)GetObject (0))->vOrigin = 0;
            if (treeFlags & HY_TREEPANEL_SCALE_TO_WINDOW) {
                ((_HYStretchCanvas*)GetObject (2))->SetMessageRecipient (nil);
            }
            SetVSpace (newVSpace>0?newVSpace:1,false,!force);
            SetHSpace (newHSpace>0?newHSpace:1,force);
            if (treeFlags & HY_TREEPANEL_SCALE_TO_WINDOW) {
                ((_HYStretchCanvas*)GetObject (2))->SetMessageRecipient (this);
            }
        }
    }
}

//__________________________________________________________

void    _HYTreePanel::SetFlags (unsigned int f)
{
    if (treeFlags!=f) {
        if ((treeFlags & HY_TREEPANEL_SCALED) != (f & HY_TREEPANEL_SCALED)) {
            treeFlags = f;
            if (!BuildTree(false)) {
                return;
            }
        } else {
            if ((treeFlags&HY_TREEPANEL_CIRCULAR)&&(!(f&HY_TREEPANEL_CIRCULAR))) {
                _TheTree* me = LocateMyTreeVariable();
                long      tipCount;
                node<nodeCoord>* newTree;
                if (!scaleVariable.sLength) {
                    newTree = me->AlignedTipsMapping(true);
                } else {
                    char     mapMode;
                    _String  scalerString = me->DetermineBranchLengthMappingMode(&scaleVariable, mapMode);
                    newTree = me->ScaledBranchMapping(nil,&scalerString,0,tipCount,mapMode);
                }
                CopyNodeParameters (newTree,coordTree);
                coordTree->delete_tree();
                delete  (coordTree);
                coordTree = newTree;
                Convert2ScreenCoordinates (hScale,vScale,-coordTree->in_object.h,coordTree);
            } else {
                if ((f&HY_TREEPANEL_CIRCULAR)&&(!(treeFlags&HY_TREEPANEL_CIRCULAR))) {
                    _HYCanvas* theCanvas = (_HYCanvas*)GetObject (0);
                    CircularLayoutPass1 (coordTree);
                    CircularLayoutPass2 (coordTree,theCanvas->_HYComponent::GetMaxW()/2,
                                         theCanvas->_HYComponent::GetMaxH()/2,arcStart,arcEnd,0,nil);
                }
            }
            treeFlags = f;
        }
        RenderTree();
    }
}

//__________________________________________________________

void    _HYTreePanel::SetScaleVariable (_String& varName)
{
    _HYPullDown *sm = (_HYPullDown*)GetObject (5);
    if (sm->FindMenuItem (varName)<0) {
        return;
    }

    if (varName!=scaleVariable) {
        scaleVariable = varName;
        if (treeFlags&HY_TREEPANEL_SCALED && BuildTree(true)) {
            RenderTree();
        }
    }
}

//__________________________________________________________

void    _HYTreePanel::RenderTree2 (void)
{
    _HYCanvas* theCanvas = (_HYCanvas*)GetObject (0);
    if (treeFlags&HY_TREEPANEL_ARCS) {
        PaintArcs(theCanvas, coordTree);
    } else if (treeFlags&HY_TREEPANEL_STRAIGHT) {
        PaintStraight(theCanvas, coordTree);
    } else if (treeFlags&HY_TREEPANEL_CIRCULAR) {
        PaintRadial (theCanvas, coordTree);
    } else {
        if (IsVertical()) {
            PaintVSquare(theCanvas, coordTree);
        } else {
            PaintSquare(theCanvas, coordTree);
            if (treeFlags&HY_TREEPANEL_LABEL1) {
                theCanvas->SetFont(branchLabel1);
                PaintSquareBranchLabels (theCanvas,coordTree,true);
            }
            if (treeFlags&HY_TREEPANEL_LABEL2) {
                theCanvas->SetFont(branchLabel2);
                PaintSquareBranchLabels (theCanvas,coordTree,false);
            }
        }
    }
}

//__________________________________________________________

void    _HYTreePanel::RenderTree (bool map, bool printing)
{
    _HYCanvas* theCanvas = (_HYCanvas*)GetObject (0);

    if (map) {
        if (treeFlags&HY_TREEPANEL_PROJECTION) {
            treeFlags -= HY_TREEPANEL_PROJECTION;
            node<nodeCoord>* thisNode = NodeTraverser(coordTree);
            while (thisNode) {
                thisNode->in_object.h = thisNode->in_object.auxL;
                thisNode->in_object.v = thisNode->in_object.auxD;
                thisNode = NodeTraverser((node<nodeCoord>*)NULL);
            }
        }
        if (treeFlags&HY_TREEPANEL_STRAIGHT) {
            PreStraightenEdges (coordTree);
            StraightenEdges(coordTree,0,windowTextMarginV,-1);
        } else {
            if (!(treeFlags&HY_TREEPANEL_CIRCULAR)) {
                _Parameter offset = windowTextMarginV;
                SpaceNodes (coordTree,offset);
            }
        }
    }
    if (!printing) {
        theCanvas->StartDraw();
    }

    theCanvas->SetFont(treeLabelFont);
    theCanvas->EraseAll();

    RenderTree2();

    if (!printing) {
        theCanvas->EndDraw();
        theCanvas->_MarkForUpdate();
        RenderNavTree ();
    }
    RenderRuler(1.0,printing);
}

//__________________________________________________________

void    _HYTreePanel::RenderRuler (_Parameter hsc, bool printing, int hshift, int vshift)
{
    _HYCanvas* theRuler = (_HYCanvas*)GetObject (2);
    if (!printing) {
        theRuler->StartDraw();
        theRuler->EraseAll();
    }
    if ((scaleVariable.sLength)&&(!IsVertical())) {
        long       rulerWidth = treeLength*hScale*hsc, scaleStep, hashScale,k,leftOver,oldRulerW;
        _Parameter tickScale, flScale;

        _HYRect  r = {HY_TREEPANEL_RULER_EXPANDED-4+vshift,windowTextMarginH+hshift,
                      HY_TREEPANEL_RULER_EXPANDED-4+vshift,windowTextMarginV+rulerWidth+hshift,1
                     };

        hashScale = floor(log10(treeLength));
        tickScale = exp(log(10.0)*hashScale);
        hashScale = floor(treeLength/tickScale);
        if (hashScale<2) {
            tickScale/=10.0;
            hashScale = floor(treeLength/tickScale);
        }

        _HYFont  rulerFont;
#ifdef __WINDOZE__
        rulerFont.face = "LucidaConsole";
        rulerFont.size = 10;
#else
#ifdef __HYPHY_GTK__
        rulerFont.face = _HY_SANS_FONT;
        rulerFont.size = 9;
#else
        rulerFont.face = "Helvetica";
        rulerFont.size = 9;
#endif
#endif

        rulerFont.style = HY_FONT_PLAIN;
        theRuler->SetFont (rulerFont);

        theRuler->DrawLine(r);
        r.right = r.left;
        r.top   -=  8;
        theRuler->DrawLine(r);
        r.right +=  rulerWidth;
        r.left+=rulerWidth;
        theRuler->DrawLine(r);

        r.left  = windowTextMarginH+hshift;
        r.right = windowTextMarginH+hshift;
        char    buff[255];
        snprintf (buff, sizeof(buff),"%.4g",treeLength);
        theRuler->DisplayText (_String (buff), HY_TREEPANEL_RULER_EXPANDED-5+vshift,windowTextMarginH+rulerWidth+3+hshift,true);
        oldRulerW = rulerWidth;
        _Parameter t = hashScale*tickScale/treeLength;
        rulerWidth*=t;
        scaleStep = rulerWidth/hashScale;
        leftOver = rulerWidth-scaleStep*hashScale;
        flScale = tickScale;
        tickScale = ((_Parameter)leftOver)/hashScale;
        t = 0.0;
        for (k=1; k<=hashScale; k++) {
            r.left+=scaleStep;
            r.right+=scaleStep;
            t+=tickScale;
            if (t>=1.0) {
                r.left++;
                r.right++;
                t-=1.0;
            }
            theRuler->DrawLine(r);
            _String l (flScale*k);
            if (GetVisibleStringWidth(l,rulerFont)<MIN(scaleStep-8,windowTextMarginH+oldRulerW+hshift-r.right-8)) {
                theRuler->DisplayText (l, r.bottom-2,r.left+4,true);
            }
        }
    }
    if (!printing) {
        theRuler->EndDraw();
        theRuler->_MarkForUpdate();
    }
}


//__________________________________________________________

void    _HYTreePanel::RenderNavTree (void)
{
    _HYCanvas* theCanvas = (_HYCanvas*)GetObject (1),
               * bigCanvas = (_HYCanvas*)GetObject (0);
    _Parameter hScale, vScale;
    hScale = (HY_TREEPANEL_NAVSIZE-4)/(_Parameter)(bigCanvas->GetCanvasSize().right);
    vScale = (HY_TREEPANEL_NAVSIZE-4)/(_Parameter)(bigCanvas->GetCanvasSize().bottom);
    theCanvas->StartDraw();
    theCanvas->EraseAll();
    _HYRect r;
    r.top = 0;
    r.left = 0;
    r.bottom = HY_TREEPANEL_NAVSIZE-1;
    r.right = HY_TREEPANEL_NAVSIZE-1;
    r.width = 1;
    theCanvas->SetColor (navFill);
    theCanvas->FillRect (r);
    theCanvas->SetColor (black);
    theCanvas->DrawRect (r);
    //theCanvas->DrawFrame (hyDefaultFrameColor);
    if (treeFlags&HY_TREEPANEL_ARCS) {
        IsVertical()?PaintNArcs(theCanvas, coordTree, vScale, hScale):PaintNArcs(theCanvas, coordTree, hScale, vScale);
    } else if (treeFlags&(HY_TREEPANEL_STRAIGHT|HY_TREEPANEL_CIRCULAR)) {
        IsVertical()?PaintNStraight(theCanvas, coordTree, vScale, hScale):PaintNStraight(theCanvas, coordTree, hScale, vScale);
    } else {
        IsVertical()?PaintNVSquare(theCanvas, coordTree, vScale, hScale):PaintNSquare(theCanvas, coordTree, hScale, vScale);
    }
    theCanvas->EndDraw();
#ifndef __MAC__
    theCanvas->_MarkForUpdate();
#endif
}

//__________________________________________________________

void    _HYTreePanel::Update (Ptr p)
{
    _HYTWindow::Update(p);
    if ((treeFlags&HY_TREEPANEL_SCALE_TO_WINDOW) == 0) {
#ifdef __MAC__
        _HYCanvas* theCanvas = (_HYCanvas*)GetObject (1);
        forceUpdateForScrolling=true;
        theCanvas->_MarkForUpdate();
        forceUpdateForScrolling=false;
#endif
        _PaintNavRect();
    } else {
        _PaintLFStatus();
    }
}

//__________________________________________________________

void    _HYTreePanel::Paint (Ptr p)
{
    _HYTWindow::Paint(p);
    if ((treeFlags&HY_TREEPANEL_SCALE_TO_WINDOW) == 0) {
#ifdef __MAC__
        _HYCanvas* theCanvas = (_HYCanvas*)GetObject (1);
        forceUpdateForScrolling=true;
        theCanvas->_MarkForUpdate();
        forceUpdateForScrolling=false;
#endif
        _PaintNavRect();
    } else {
        _PaintLFStatus();
    }
}

//__________________________________________________________

void    _HYTreePanel::SetVertical (bool toggle)
{
    _HYCanvas* theCanvas = (_HYCanvas*)GetObject (0);
    int        treeHeight = 2*windowTextMarginV+tips*vSpace,
               treeWidth = 2*windowTextMarginH+textSpace+coordTree->tree_depth()*hSpace;
    if (toggle) {
        if (!(treeFlags&HY_TREEPANEL_VERTICAL)) {
            treeFlags|=HY_TREEPANEL_VERTICAL;
            theCanvas->SetOrigin(0,0);
            ResizeTreeCanvas (treeWidth-textSpace+windowTextMarginH,treeHeight+textSpace);
            RenderTree();
            RenderRuler();
        }
    } else if (treeFlags&HY_TREEPANEL_VERTICAL) {
        treeFlags-=HY_TREEPANEL_VERTICAL;
        theCanvas->SetOrigin(0,0);
        ResizeTreeCanvas (treeHeight,treeWidth);
        RenderTree();
        RenderRuler();
    }
}

//__________________________________________________________
void    _HYTreePanel::ResizeTreeCanvas (int newH, int newW)
{

    _HYCanvas* theCanvas = (_HYCanvas*)GetObject (0);
    _HYRect    myDims = theCanvas->GetCanvasSize();
    if (newH > HY_TREEPANEL_MAXSIZE) {
        newH = HY_TREEPANEL_MAXSIZE;
    }
    if (newW > HY_TREEPANEL_MAXSIZE) {
        newW = HY_TREEPANEL_MAXSIZE;
    }
    int visH = componentB.lData[0]+theCanvas->vOrigin-componentT.lData[0],
        visW = theCanvas->hOrigin+componentR.lData[0]-componentL.lData[0];

    bool change = false;
    if (newH>myDims.bottom) {
        change = true;
    } else {
        if (newH<myDims.bottom) {
            change = true;
            if (newH<visH) {
                newH = visH;
            }
        }
    }
    if (newW>myDims.right) {
        change = true;
    } else {
        if (newW<myDims.right) {
            change = true;
            if (newW<visW) {
                newW = visW;
            }
        }
    }
    if (change) {
        theCanvas->SetCanvasSize (newH,newW,TREE_PANEL_PIXEL_DEPTH);
        myDims.left = componentL.lData[0];
        myDims.top = componentT.lData[0];
        myDims.right = componentR.lData[0];
        myDims.bottom = componentB.lData[0];
        theCanvas->SetVisibleSize (myDims);

        dim = MinMaxWindowDimensions();
    }
}

//__________________________________________________________

bool    _HYTreePanel::BuildTree (bool saveOrigin)
{
    _TheTree* me = LocateMyTreeVariable();
    if      (!me) {
        return false;
    }

    dumpCoordTree         ();
    currentSelection.Clear();
    hScale = vScale = 1.0;

    long            f,
                    tipCount;

    _Parameter      treeWidth,
                    treeHeight,
                    visTreeWidth,
                    visTreeHeight;

    _String         newStatusString;
    bool            scaling = false;


    if ((treeFlags&HY_TREEPANEL_SCALED) && scaleVariable.sLength) {
        scaling = true;
    }

    if (!scaling) {
        coordTree = me->AlignedTipsMapping(true);
    } else {
        char     mapMode;
        _String  scalerString = me->DetermineBranchLengthMappingMode(&scaleVariable, mapMode);
        coordTree = me->ScaledBranchMapping(nil,&scalerString,0,tipCount,mapMode);
    }

    treeWidth       = -coordTree->in_object.h;
    tipCount        = coordTree->tree_depth();
    newStatusString = newStatusString&"Depth = "&tipCount;

    _HYCanvas* theCanvas = (_HYCanvas*)GetObject (0),
               * theRuler  = (_HYCanvas*)GetObject (2);
    textSpace = 0;
    node<nodeCoord>* currentNd = NodeTraverser(coordTree);
    while (currentNd) {
        if (currentNd->in_object.varRef>=0) {
            _String* intLabel = LocateVar(currentNd->in_object.varRef)->GetName();
            currentNd->in_object.branchName = intLabel->Cut(intLabel->FindBackwards('.',0,-1)+1,-1);
            f = GetVisibleStringWidth (currentNd->in_object.branchName,treeLabelFont);
            currentNd->in_object.textWidth = f;
            if ((currentNd->get_num_nodes()==0)&&(f>textSpace)) {
                textSpace = f;
            }
        } else {
            currentNd->in_object.textWidth = 0;
            currentNd->in_object.branchName = empty;
        }
        currentNd->in_object.flags = 0;
        currentNd = NodeTraverser((node<nodeCoord>*)NULL);
    }

    visTreeWidth = (tipCount-1)*hSpace;
    hScale = visTreeWidth/treeWidth;
    visTreeWidth += 2*windowTextMarginH
                    +((treeFlags&HY_TREEPANEL_TIP_LABELS)?textSpace:0);
    if (visTreeWidth>HY_TREEPANEL_MAXSIZE) {
        visTreeWidth=HY_TREEPANEL_MAXSIZE;
        hScale = (visTreeWidth-2*windowTextMarginH
                  -((treeFlags&HY_TREEPANEL_TIP_LABELS)?textSpace:0))/treeWidth;
    }

    currentNd = coordTree;

    tipCount = currentNd->get_num_nodes();

    while (tipCount) {
        currentNd = currentNd->go_down(1);
        tipCount = currentNd->get_num_nodes();
    }

    treeHeight = currentNd->in_object.v;
    currentNd = coordTree;
    tipCount = currentNd->get_num_nodes();

    while (tipCount) {
        currentNd = currentNd->go_down(tipCount);
        tipCount = currentNd->get_num_nodes();
    }

    treeHeight = currentNd->in_object.v - treeHeight;
    _PMathObj   leafCount = me->TipCount();
    newStatusString = newStatusString&". Leaf Count = "&_String((long)leafCount->Value());
    tips = leafCount->Value();
    if (likeFuncID>=0) {
        _LikelihoodFunction* theLF = (_LikelihoodFunction*)likeFuncList(likeFuncID);
        theLF->PrepareToCompute();
        newStatusString = newStatusString&". Ln-likelihood = "&_String(theLF->Compute());
        if (!isInOptimize) {
            theLF->DoneComputing();
        }
    }
    if (scaling) {
        treeLength = treeWidth;
    } else {
        treeLength = 0;
    }
    visTreeHeight = (tips-1)*vSpace;
    vScale = visTreeHeight/treeHeight;
    visTreeHeight+=2*windowTextMarginV;
    if (visTreeHeight>HY_TREEPANEL_MAXSIZE) {
        visTreeHeight=HY_TREEPANEL_MAXSIZE;
        vScale = (visTreeHeight-2*windowTextMarginV)/treeHeight;
    }
    Convert2ScreenCoordinates (hScale,vScale,-coordTree->in_object.h,coordTree);
    _HYRect    myDims = theCanvas->GetCanvasSize();
    scaling = false;
    if (visTreeHeight!=myDims.bottom) {
        scaling = true;
        myDims.bottom = visTreeHeight;
        if (myDims.bottom<HY_TREEPANEL_NAVSIZE) {
            myDims.bottom = HY_TREEPANEL_NAVSIZE;
        }
    }
    if (visTreeWidth!=myDims.right) {
        scaling = true;
        myDims.right = visTreeWidth;
        if (myDims.right<HY_TREEPANEL_DEFSIZE) {
            myDims.right = HY_TREEPANEL_DEFSIZE;
        }
    }
    long  saveH = theCanvas->hOrigin,
          saveV = theCanvas->vOrigin;
    Relabel (me);
    if (scaling) {
        if (IsVertical()) {
            theCanvas->SetCanvasSize (myDims.right,MAX(myDims.bottom,HY_TREEPANEL_DEFSIZE),TREE_PANEL_PIXEL_DEPTH);
            theRuler->SetCanvasSize(HY_TREEPANEL_RULER_EXPANDED,MAX(myDims.bottom,HY_TREEPANEL_DEFSIZE),TREE_PANEL_PIXEL_DEPTH);
        } else {
            theCanvas->SetCanvasSize (myDims.bottom,myDims.right,TREE_PANEL_PIXEL_DEPTH);
            theRuler->SetCanvasSize(HY_TREEPANEL_RULER_EXPANDED,myDims.right,TREE_PANEL_PIXEL_DEPTH);
        }
        UpdateComponentInfo();
        if (saveOrigin) {
            theCanvas->SetOrigin (saveV,saveH);
            if (!IsVertical()) {
                theRuler->SetOrigin (0,saveH);
            }
        }
    }

    if (treeFlags&HY_TREEPANEL_SCALE_TO_WINDOW) {
        FitToWindow();
    }

    if (treeFlags&HY_TREEPANEL_CIRCULAR) {
        CircularLayoutPass1 (coordTree);
        CircularLayoutPass2 (coordTree,theCanvas->_HYComponent::GetMaxW()/2,
                             theCanvas->_HYComponent::GetMaxH()/2,arcStart,arcEnd,0,nil);
    }


    DeleteObject (leafCount);
    SetStatusBar (newStatusString);
    return true;
}

//__________________________________________________________

void  _HYTreePanel::Convert2ScreenCoordinates (_Parameter hScale, _Parameter vScale, _Parameter offset, node<nodeCoord>* cNode)
{
    cNode->in_object.h = (cNode->in_object.h+offset)*hScale+windowTextMarginH;
    cNode->in_object.v = cNode->in_object.v*vScale+windowTextMarginV;
    for (int i=cNode->get_num_nodes(); i; i--) {
        Convert2ScreenCoordinates (hScale,vScale,offset,cNode->go_down(i));
    }
}

//__________________________________________________________

void  _HYTreePanel::ShiftScreenCoordinates (_Parameter hShift, _Parameter vShift,node<nodeCoord>* cNode)
{
    cNode->in_object.h += hShift;
    cNode->in_object.v += vShift;
    for (int i=cNode->get_num_nodes(); i; i--) {
        ShiftScreenCoordinates (hShift,vShift,cNode->go_down(i));
    }
}

//__________________________________________________________

void  _HYTreePanel::StraightenEdges (node<nodeCoord>* cNode,_Parameter slope, _Parameter offset, int m)
{
    node<nodeCoord>* theParent = cNode->get_parent(), *c;
    int n = cNode->get_num_nodes(),j,k=0;
    _Parameter t;
    for (j=1; j<=n/2; j++) {
        k+=cNode->go_down(j)->in_object.auxL;
    }

    if (!theParent) {
        cNode->in_object.v = k*vSpace+offset;
    }

    k = 0;
    for (j=1; j<=n/2; j++) {
        c = cNode->go_down(j);
        k+=c->in_object.auxL;
        if (j!=m) {
            t = -(cNode->in_object.v-offset)/(c->in_object.auxD-cNode->in_object.h);
        } else {
            t = slope;
        }
        c->in_object.v = cNode->in_object.v+t*(c->in_object.h-cNode->in_object.h);
        StraightenEdges (c,t,offset,j);
        offset+=c->in_object.auxL*vSpace;
    }
    k=0;
    for (j=n/2+1; j<=n; j++) {
        c = cNode->go_down(j);
        k+=c->in_object.auxL;
        if (j!=m) {
            t = ((k-1)*vSpace+offset-cNode->in_object.v)/(c->in_object.auxD-cNode->in_object.h);
        } else {
            t = slope;
        }
        c->in_object.v = cNode->in_object.v+t*(c->in_object.h-cNode->in_object.h);
        StraightenEdges (c,t,offset,j);
        offset+=c->in_object.auxL*vSpace;
    }
}

//__________________________________________________________

long  _HYTreePanel::CircularLayoutPass1 (node<nodeCoord>* cNode)
/* this function will populate the auxL member of nodeCoord
   with the total number of leaves at or below this node
*/
{
    int n = cNode->get_num_nodes();
    if (n) {
        cNode->in_object.auxL = 0;
        for (; n; n--) {
            cNode->in_object.auxL += CircularLayoutPass1 (cNode->go_down(n));
        }
    } else {
        cNode->in_object.auxL = 1;
    }
    return cNode->in_object.auxL;
}

//__________________________________________________________

void  _HYTreePanel::CircularLayoutPass2 (node<nodeCoord>* cNode, _Parameter cH, _Parameter cV, _Parameter arc1, _Parameter arc2, _Parameter d, long* minMax)
{
    long n  =   cNode->get_num_nodes(),
         k ;
    if (n) { // an internal node
        _Parameter step     = (arc2-arc1)/cNode->in_object.auxL; // angle per child
        node<nodeCoord>*      tNode;
        if (cNode->parent) {
            _Parameter dn = d-cNode->parent->in_object.h+cNode->in_object.h;
            for (k=1; k<=n; k++) {
                tNode = cNode->go_down(k);
                arc2 = arc1 + step*tNode->in_object.auxL;
                CircularLayoutPass2 (tNode,cH,cV,arc1,arc2,dn,minMax);
                arc1 = arc2;
            }
        } else {
            minMax = new long[4];
            minMax[0] = cH-2;
            minMax[1] = cH+2;
            minMax[2] = cV-2;
            minMax[3] = cV+2;
            for (k=1; k<=n; k++) {
                tNode = cNode->go_down(k);
                arc2 = arc1 + step*tNode->in_object.auxL;
                CircularLayoutPass2 (tNode,cH,cV,arc1,arc2,0,minMax);
                arc1 = arc2;
            }
        }
        arc1 = arc2 - step*cNode->in_object.auxL;
    }
    if (cNode->parent) {
        arc1 = (arc1+arc2)*0.5;
        d += cNode->in_object.h - cNode->parent->in_object.h;
        d *= .25;
        cNode->in_object.h = cH+d*cos (arc1);
        cNode->in_object.v = cV-d*sin (arc1);
        minMax[0] = MIN(minMax[0],cNode->in_object.h);
        minMax[1] = MAX(minMax[1],cNode->in_object.h);
        minMax[2] = MIN(minMax[2],cNode->in_object.v);
        minMax[3] = MAX(minMax[3],cNode->in_object.v);
        cNode->in_object.h -= cH;
        cNode->in_object.v -= cV;
    } else { // this is the root node
        if (treeFlags&HY_TREEPANEL_TIP_LABELS) {
            arc1 = (_Parameter)(cH*2-windowTextMarginH-2*textSpace)/(minMax[1]-minMax[0]);
        } else {
            arc1 = (_Parameter)(cH*2-windowTextMarginH)/(minMax[1]-minMax[0]);
        }
        arc2 = (_Parameter)(cV*2-2*windowTextMarginV-treeLabelFont.size)/(minMax[3]-minMax[2]);
        cH = 2.*cH*(cH-minMax[0])/(minMax[1]-minMax[0]);
        cV = 2.*cV*(cV-minMax[2])/(minMax[3]-minMax[2]);
        cNode->in_object.h = 0;
        cNode->in_object.v = 0;
        CircularLayoutPass3 (cNode,cH,cV,arc1,arc2);
        delete minMax;
    }
    cNode->in_object.auxL = 0;
}

//__________________________________________________________

void  _HYTreePanel::CircularLayoutPass3 (node<nodeCoord>* cNode, _Parameter cH, _Parameter cV, _Parameter hs, _Parameter vs)
{
    int n = cNode->get_num_nodes();
    if (n) {
        for (; n; n--) {
            CircularLayoutPass3 (cNode->go_down(n),cH,cV,hs,vs);
        }
    }
    cNode->in_object.h = cH+(cNode->in_object.h)*hs;
    cNode->in_object.v = cV+(cNode->in_object.v)*vs;
}
//__________________________________________________________

long  _HYTreePanel::FindSelection (long h, long v, char flag)
{
    long res = 0, addWidth;
    bool shiftDown   = (flag&0x01)>0,
         optionDown  = (flag&0x02)>0,
         commandDown = ((flag&0x04)>0)&&selectionTop;
    if (commandDown)
        if (!selectionTop->parent->parent) {
            commandDown = false;
        }
    node<nodeCoord>* currentNd = NodeTraverser(coordTree), *parent = nil;

    if (!(shiftDown||commandDown)) {
        currentSelection.Clear();
    }
    while (currentNd) {
        addWidth = 0;
        parent = currentNd->parent;
        if (!parent) {
            break;
        }
        if (treeFlags&HY_TREEPANEL_TIP_LABELS) {
            if (currentNd->get_num_nodes()==0) {
                addWidth = currentNd->in_object.textWidth+(treeLabelFont.size/2);
            }
        }
        if (treeFlags&HY_TREEPANEL_INT_LABELS) {
            if (currentNd->get_num_nodes()!=0) {
                addWidth = currentNd->in_object.textWidth+(treeLabelFont.size/2);
            }
        }

        if (treeFlags&HY_TREEPANEL_CIRCULAR) {
            if (((h<currentNd->in_object.h)&&(h>parent->in_object.h))||
                    ((h>currentNd->in_object.h)&&(h<parent->in_object.h))) {
                _Parameter lineSlope = ((_Parameter)currentNd->in_object.v-parent->in_object.v)/
                                       ((_Parameter)currentNd->in_object.h-parent->in_object.h);
                if (abs((long)(v-parent->in_object.v-lineSlope*(h-parent->in_object.h)))<=3) {
                    break;
                }
            }

        } else if ((h<currentNd->in_object.h+addWidth)&&(h>parent->in_object.h)) {
            if (treeFlags&HY_TREEPANEL_STRAIGHT) {
                _Parameter lineSlope = ((_Parameter)currentNd->in_object.v-parent->in_object.v)/
                                       ((_Parameter)currentNd->in_object.h-parent->in_object.h);
                if (abs((long)(v-parent->in_object.v-lineSlope*(h-parent->in_object.h)))<=5) {
                    break;
                }
            } else if (treeFlags&HY_TREEPANEL_ARCS) ;
            else {
                if (abs((long)(v-currentNd->in_object.v))<=5) {
                    break;
                }
            }
        }
        if (!shiftDown) {
            if (currentNd->in_object.flags&HY_BRANCH_SELECT) {
                currentNd->in_object.flags&=HY_BRANCH_DESELECT;
                res = 1;
            }
        }
        currentNd = NodeTraverser((node <nodeCoord>*)nil);
    }
    if (currentNd&&parent&&currentNd->in_object.branchName.sLength) {
        if (commandDown) {
            if (currentNd->in_object.flags&HY_BRANCH_SELECT) {
                _String errMsg ("Can't move a subtree to its own branch! Try again.");
                ProblemReport (errMsg,(Ptr)this);
                return 0;
            } else {
                currentSelection<<(long)currentNd;
                MoveSubTree();
                return 1;
            }
        }
        if (optionDown) {
            node<nodeCoord>* saveNode = currentNd;
            if (!shiftDown) {
                currentNd = NodeTraverser((node <nodeCoord>*)nil);
                while (currentNd) {
                    currentNd->in_object.flags&=HY_BRANCH_DESELECT;
                    currentNd = NodeTraverser((node <nodeCoord>*)nil);
                }
            }
            SelectSubTree (saveNode,shiftDown);
            res = 1;
        } else {
            if (currentNd->in_object.flags&HY_BRANCH_SELECT) {
                if (shiftDown) {
                    currentNd->in_object.flags&=HY_BRANCH_DESELECT;
                    currentSelection.Delete (currentSelection.Find((long)currentNd));
                } else {
                    currentSelection<<(long)currentNd;
                }
            } else {
                currentNd->in_object.flags|=HY_BRANCH_SELECT;
                currentSelection<<(long)currentNd;
            }
            res = 1;
            if (!shiftDown) {
                currentNd = NodeTraverser((node <nodeCoord>*)nil);
                while (currentNd) {
                    currentNd->in_object.flags&=HY_BRANCH_DESELECT;
                    currentNd = NodeTraverser((node <nodeCoord>*)nil);
                }
            }
        }
    }
    return res;
}

//__________________________________________________________

long  _HYTreePanel::FindSelectedBranch (long h, long v)
{
    if (coordTree) {
        long            addWidth;

        node<nodeCoord>* currentNd = NodeTraverser(coordTree),
                         * parent = nil;

        while (currentNd) {
            addWidth = 0;
            parent = currentNd->parent;
            if (!parent) {
                break;
            }
            if (treeFlags&HY_TREEPANEL_TIP_LABELS) {
                if (currentNd->get_num_nodes()==0) {
                    addWidth = currentNd->in_object.textWidth+(treeLabelFont.size/2);
                }
            }
            if (treeFlags&HY_TREEPANEL_INT_LABELS) {
                if (currentNd->get_num_nodes()!=0) {
                    addWidth = currentNd->in_object.textWidth+(treeLabelFont.size/2);
                }
            }

            if (treeFlags&HY_TREEPANEL_CIRCULAR) {
                if (((h<currentNd->in_object.h)&&(h>parent->in_object.h))||
                        ((h>currentNd->in_object.h)&&(h<parent->in_object.h))) {
                    _Parameter lineSlope = ((_Parameter)currentNd->in_object.v-parent->in_object.v)/
                                           ((_Parameter)currentNd->in_object.h-parent->in_object.h);
                    if (abs((long)(v-parent->in_object.v-lineSlope*(h-parent->in_object.h)))<=3) {
                        break;
                    }
                }

            } else if ((h<currentNd->in_object.h+addWidth)&&(h>parent->in_object.h)) {
                if (treeFlags&HY_TREEPANEL_STRAIGHT) {
                    _Parameter lineSlope = ((_Parameter)currentNd->in_object.v-parent->in_object.v)/
                                           ((_Parameter)currentNd->in_object.h-parent->in_object.h);
                    if (abs((long)(v-parent->in_object.v-lineSlope*(h-parent->in_object.h)))<=5) {
                        break;
                    }
                } else if (treeFlags&HY_TREEPANEL_ARCS) ;
                else {
                    if (abs((long)(v-currentNd->in_object.v))<=5) {
                        break;
                    }
                }
            }
            currentNd = NodeTraverser((node <nodeCoord>*)nil);
        }
        if (parent) {
            return (long)currentNd;
        }
    }
    return 0;
}

//__________________________________________________________

void  _HYTreePanel::SelectSubTree (node<nodeCoord>*thisNode, bool shiftDown)
{
    if (!thisNode->in_object.flags&HY_BRANCH_SELECT) {
        thisNode->in_object.flags|=HY_BRANCH_SELECT;
        currentSelection<<(long)thisNode;
    } else if (shiftDown) {
        thisNode->in_object.flags&=HY_BRANCH_DESELECT;
        currentSelection.Delete (currentSelection.Find((long)thisNode));
    }
    for (long k=0; k<thisNode->nodes.length; k++) {
        SelectSubTree (thisNode->nodes.data[k],shiftDown);
    }
}

//__________________________________________________________

void  _HYTreePanel::SelectRangeAndScroll (_List& nodeNames, bool upcaseMe)
{
    if (nodeNames.lLength) {
        currentSelection.Clear();

        node<nodeCoord>* currentNd = NodeTraverser(coordTree);
        _Parameter  ch = 0,
                    cv = 0;
        while (currentNd) {
            currentNd->in_object.flags &= HY_BRANCH_DESELECT;
            if (currentNd->in_object.varRef>=0) {
                bool addMe = false;
                _String* nodeName = LocateVar(currentNd->in_object.varRef)->GetName();
                if (upcaseMe) {
                    _String upcased (*nodeName);
                    upcased.UpCase();
                    addMe = nodeNames.BinaryFindObject (&upcased)>=0;
                } else {
                    addMe = nodeNames.BinaryFindObject (nodeName)>=0;
                }

                if (addMe) {
                    ch += currentNd->in_object.h;
                    cv += currentNd->in_object.v;
                    currentNd->in_object.flags |= HY_BRANCH_SELECT;
                    currentSelection << (long)currentNd;
                }
            }
            currentNd = NodeTraverser((node<nodeCoord>*)NULL);
        }

        _UpdateOperationsMenu();
        RenderTree();
        ScrollToSelection (ch, cv);
    }
}

//__________________________________________________________


void  _HYTreePanel::ScrollToSelection (_Parameter ch, _Parameter cv)
{
    if (currentSelection.lLength&&((treeFlags&HY_TREEPANEL_SCALE_TO_WINDOW)==0)) {
        _HYCanvas * bigCanvas = (_HYCanvas*)GetObject (0);
        ch /= currentSelection.lLength;
        cv /= currentSelection.lLength;
        ch /= bigCanvas->GetMaxW();
        cv /= bigCanvas->GetMaxH();
        SetNavRectCenter (ch*HY_TREEPANEL_NAVSIZE,cv*HY_TREEPANEL_NAVSIZE);
    }
}

//__________________________________________________________

void  _HYTreePanel::SpaceNodes (node<nodeCoord>* cNode,_Parameter& offset)
{
    int n = cNode->get_num_nodes(),j;
    if (n) {
        _Parameter t = 0.0;
        for (j=1; j<=n; j++) {
            SpaceNodes (cNode->go_down(j),offset);
            t+=cNode->go_down(j)->in_object.v;
        }
        cNode->in_object.v = t/n;
    } else {
        cNode->in_object.v = offset;
        offset+=vSpace;
    }

}

//__________________________________________________________

_Parameter  _HYTreePanel::PreStraightenEdges (node<nodeCoord>* cNode)
{
    int i = cNode->get_num_nodes();
    if (i) {
        _Parameter t = 0.0;
        cNode->in_object.auxL = 0;
        cNode->in_object.auxD = 0.0;
        for (int j=1; j<=i; j++) {
            t = PreStraightenEdges (cNode->go_down(j));
            if (cNode->in_object.auxD<t) {
                cNode->in_object.auxD = t;
            }
            cNode->in_object.auxL+=cNode->go_down(j)->in_object.auxL;
        }
        return cNode->in_object.auxD;
    } else {
        cNode->in_object.auxL = 1;
        return cNode->in_object.auxD = cNode->in_object.h;
    }
}

//__________________________________________________________

long                lastLeafValue;
node  <nodeCoord>*  lastLeafPointer;

//__________________________________________________________

long  _HYTreePanel::PaintSquare (_HYCanvas* theCanvas, node<nodeCoord>* thisNode)
{
    _HYRect nodeCircle;
    node<nodeCoord>* parent = thisNode->get_parent();
    int     mh, mv;
    mh = thisNode->in_object.h;
    mv = thisNode->in_object.v;
    nodeCircle.width = branchWidth;
    if (parent) {
        nodeCircle.top   = nodeCircle.bottom = mv;
        nodeCircle.right = mh;
        nodeCircle.left = parent->in_object.h;
        DrawStraightEdge (theCanvas,nodeCircle,thisNode);
    } else {
        lastLeafValue = -1;
        lastLeafPointer = nil;
    }

    int n,t,b;
    if ((n=thisNode->get_num_nodes())) {
        b = PaintSquare (theCanvas,thisNode->go_down(n));
        t = b;
        n--;
        for (; n>1; n--) {
            PaintSquare (theCanvas,thisNode->go_down(n));
        }
        if (n) {
            t = PaintSquare (theCanvas,thisNode->go_down(1));
        }
        if (t!=b) {
            nodeCircle.left = nodeCircle.right = mh;
            nodeCircle.top = t;
            nodeCircle.bottom = mv;
            theCanvas->DrawLine (nodeCircle);
            nodeCircle.top = mv;
            nodeCircle.bottom = b;
            theCanvas->DrawLine (nodeCircle);
        }
        if (parent&&(treeFlags&HY_TREEPANEL_INT_LABELS)) {
            if (thisNode->in_object.branchName.sLength>0) {
                n = b-t;
                if (n>treeLabelFont.size+1) {
                    theCanvas->DisplayText (thisNode->in_object.branchName, mv+(treeLabelFont.size/2),mh+(treeLabelFont.size/2),true);
                } else {
                    if (n>6) {
                        theCanvas->SetFontSize(n-1);
                        theCanvas->DisplayText (thisNode->in_object.branchName, mv+(treeLabelFont.size/2),mh+(treeLabelFont.size/2),true);
                        theCanvas->SetFontSize(treeLabelFont.size);
                    }
                }
            }
        }
    } else {
        if (treeFlags&HY_TREEPANEL_TIP_LABELS) {
            if (thisNode->in_object.branchName.sLength>0) {
                if (lastLeafValue < 0) {
                    lastLeafPointer = thisNode;
                } else {
                    n = lastLeafValue-mv;
                    if (n>treeLabelFont.size+1) {
                        theCanvas->DisplayText (thisNode->in_object.branchName, mv+(treeLabelFont.size/2),mh+(treeLabelFont.size/2),true);
                        if (lastLeafPointer)
                            theCanvas->DisplayText (lastLeafPointer->in_object.branchName,
                                                    lastLeafPointer->in_object.v+(treeLabelFont.size/2),
                                                    lastLeafPointer->in_object.h+(treeLabelFont.size/2),true);
                    } else {
                        if (n>6) {
                            theCanvas->SetFontSize(n-1);
                            theCanvas->DisplayText (thisNode->in_object.branchName, mv+(treeLabelFont.size/2),mh+(treeLabelFont.size/2),true);
                            if (lastLeafPointer)
                                theCanvas->DisplayText (lastLeafPointer->in_object.branchName,
                                                        lastLeafPointer->in_object.v+(treeLabelFont.size/2),
                                                        lastLeafPointer->in_object.h+(treeLabelFont.size/2),true);
                            theCanvas->SetFontSize(treeLabelFont.size);
                        }
                    }
                    lastLeafPointer = nil;
                }
            }
            lastLeafValue = mv;
        }
    }
    return mv;
}

//__________________________________________________________

void  _HYTreePanel::PaintSquareBranchLabels (_HYCanvas* theCanvas, node<nodeCoord>* thisNode, bool above)
{
    _HYRect nodeCircle;
    node<nodeCoord>* parent = thisNode->get_parent();
    int     mh, mv;
    mh = thisNode->in_object.h;
    mv = thisNode->in_object.v;
    nodeCircle.width = branchWidth;
    if (parent) {
        _Parameter value = above?thisNode->in_object.label1:thisNode->in_object.label2;

        char conv [255],
             conf [32];

        //if (value>1.0)
        snprintf (conf, sizeof(conf),"%%.%dg",above?labelDigits1:labelDigits2);
        //else
        //  snprintf (conf, sizeof(conf),"%%.%dg",(above?labelDigits1:labelDigits2)+1);

        snprintf (conv, sizeof(conv),conf,value);
        _String label (conv);

        nodeCircle.top = nodeCircle.bottom = mv;
        nodeCircle.right = mh;
        nodeCircle.left = parent->in_object.h;

#ifdef __MAC__
        long lw = GetVisibleStringWidth (label);
#else
        long lw = GetVisibleStringWidth (label,treeLabelFont);
#endif

        int firstDig = 0,
            cDig = label.sLength-1;
        while ((label.sData[firstDig]=='0')||(label.sData[firstDig]=='.')) {
            firstDig++;
        }

        if (value==0.0) {
            firstDig = 1;
            cDig = 1;
        }

        while ((lw>nodeCircle.right-nodeCircle.left-2*(branchWidth+1))&&(cDig>=firstDig)) {
            if (value>1.0) {
                snprintf (conf, sizeof(conf),"%%.%dg",cDig-firstDig);
            } else {
                snprintf (conf, sizeof(conf),"%%.%dg",cDig-firstDig+1);
            }
            snprintf (conv, sizeof(conv),conf,value);
            label = conv;
            label.Trim(0,cDig);
#ifdef __MAC__
            lw = GetVisibleStringWidth (label);
#else
            lw = GetVisibleStringWidth (label,treeLabelFont);
#endif
            cDig--;
        }
        char al = above?topAlign:bottomAlign;
        if (al==HY_ALIGN_LEFT) {
            mh = nodeCircle.left+branchWidth+1;
        } else if (al==HY_ALIGN_RIGHT) {
            mh = nodeCircle.right-lw-branchWidth-1;
        } else {
            mh = nodeCircle.left+(nodeCircle.right-nodeCircle.left-lw)/2;
        }

      if ((lw<=nodeCircle.right-nodeCircle.left-2*(branchWidth+1))&&(cDig>=firstDig)) {
            if (above) {
                theCanvas->DisplayText (label,mv-2,mh,true);
            } else {
                theCanvas->DisplayText (label,mv+branchLabel2.size+2,mh,true);
            }
      }
    }

    int n;
    if ((n=thisNode->get_num_nodes())) {
        for (; n>=1; n--) {
            PaintSquareBranchLabels (theCanvas,thisNode->go_down(n),above);
        }
    }
}

//__________________________________________________________

long  _HYTreePanel::PaintVSquare (_HYCanvas* theCanvas, node<nodeCoord>* thisNode)
{
    _HYRect nodeCircle;
    node<nodeCoord>* parent = thisNode->get_parent();
    int     mh, mv;
    mh = thisNode->in_object.h;
    mv = thisNode->in_object.v;
    nodeCircle.width = branchWidth;
    if (parent) {
        nodeCircle.left = nodeCircle.right = mv;
        nodeCircle.bottom = mh;
        nodeCircle.top = parent->in_object.h;
        DrawStraightEdge (theCanvas,nodeCircle,thisNode);
    }
    int n,t,b;
    if ((n=thisNode->get_num_nodes())) {
        if (parent&&(treeFlags&HY_TREEPANEL_INT_LABELS)) {
            if (thisNode->in_object.branchName.sLength>0) {
                theCanvas->DisplayText (thisNode->in_object.branchName, mh-3,mv+3,false);
            }
        }
        b = PaintVSquare (theCanvas,thisNode->go_down(n));
        t = b;
        n--;
        for (; n>1; n--) {
            PaintVSquare (theCanvas,thisNode->go_down(n));
        }
        if (n) {
            t = PaintVSquare (theCanvas,thisNode->go_down(1));
        }
        if (t!=b) {
            nodeCircle.top = nodeCircle.bottom = mh;
            nodeCircle.left = t;
            nodeCircle.right = mv;
            theCanvas->DrawLine (nodeCircle);
            nodeCircle.left = mv;
            nodeCircle.right = b;
            theCanvas->DrawLine (nodeCircle);
        }
    } else {
        if (treeFlags&HY_TREEPANEL_TIP_LABELS) {
            if (thisNode->in_object.branchName.sLength>0) {
                theCanvas->DisplayText (thisNode->in_object.branchName, mh+3*(treeLabelFont.size/2),mv,false);
            }
        }
    }
    return mv;
}
//__________________________________________________________

long  _HYTreePanel::PaintNSquare (_HYCanvas* theCanvas, node<nodeCoord>* thisNode, _Parameter hScale, _Parameter vScale)
{
    _HYRect nodeCircle;
    node<nodeCoord>* parent = thisNode->get_parent();
    nodeCircle.width = 1;
    if (parent) {
        nodeCircle.top = nodeCircle.bottom = thisNode->in_object.v*vScale+HY_TREEPANEL_NAVSPACING;
        nodeCircle.right = thisNode->in_object.h*hScale+HY_TREEPANEL_NAVSPACING;
        nodeCircle.left = parent->in_object.h*hScale+HY_TREEPANEL_NAVSPACING;
        theCanvas->DrawLine (nodeCircle);
    }
    int n,t,b;
    if ((n=thisNode->get_num_nodes())) {
        b = PaintNSquare (theCanvas,thisNode->go_down(n),hScale,vScale);
        t = b;
        n--;
        for (; n>1; n--) {
            PaintNSquare (theCanvas,thisNode->go_down(n),hScale,vScale);
        }
        if (n) {
            t = PaintNSquare (theCanvas,thisNode->go_down(1),hScale,vScale);
        }
        if (t!=b) {
            nodeCircle.left = nodeCircle.right = thisNode->in_object.h*hScale+HY_TREEPANEL_NAVSPACING;
            nodeCircle.top = t;
            nodeCircle.bottom = thisNode->in_object.v*vScale+HY_TREEPANEL_NAVSPACING;
            theCanvas->DrawLine (nodeCircle);
            nodeCircle.top = thisNode->in_object.v*vScale+HY_TREEPANEL_NAVSPACING;
            nodeCircle.bottom = b;
            theCanvas->DrawLine (nodeCircle);
        }
    }
    return thisNode->in_object.v*vScale+HY_TREEPANEL_NAVSPACING;
}

//__________________________________________________________

long  _HYTreePanel::PaintNVSquare (_HYCanvas* theCanvas, node<nodeCoord>* thisNode, _Parameter hScale, _Parameter vScale)
{
    _HYRect nodeCircle;
    node<nodeCoord>* parent = thisNode->get_parent();
    nodeCircle.width = 1;
    if (parent) {
        nodeCircle.left = nodeCircle.right = thisNode->in_object.v*vScale+HY_TREEPANEL_NAVSPACING;
        nodeCircle.bottom = thisNode->in_object.h*hScale+HY_TREEPANEL_NAVSPACING;
        nodeCircle.top = parent->in_object.h*hScale+HY_TREEPANEL_NAVSPACING;
        theCanvas->DrawLine (nodeCircle);
    }
    int n,t,b;
    if ((n=thisNode->get_num_nodes())) {
        b = PaintNVSquare (theCanvas,thisNode->go_down(n),hScale,vScale);
        t = b;
        n--;
        for (; n>1; n--) {
            PaintNVSquare (theCanvas,thisNode->go_down(n),hScale,vScale);
        }
        if (n) {
            t = PaintNVSquare (theCanvas,thisNode->go_down(1),hScale,vScale);
        }
        if (t!=b) {
            nodeCircle.top = nodeCircle.bottom = thisNode->in_object.h*hScale+HY_TREEPANEL_NAVSPACING;
            nodeCircle.left = t;
            nodeCircle.right = thisNode->in_object.v*vScale+HY_TREEPANEL_NAVSPACING;
            theCanvas->DrawLine (nodeCircle);
            nodeCircle.left = thisNode->in_object.v*vScale+HY_TREEPANEL_NAVSPACING;
            nodeCircle.right = b;
            theCanvas->DrawLine (nodeCircle);
        }
    }
    return thisNode->in_object.v*vScale+HY_TREEPANEL_NAVSPACING;
}

//__________________________________________________________

void  _HYTreePanel::PaintStraight (_HYCanvas* theCanvas, node<nodeCoord>* thisNode)
{
    _HYRect nodeCircle;
    node<nodeCoord>* parent = thisNode->get_parent();
    nodeCircle.width = branchWidth;
    bool    isv = IsVertical();
    long    n;
    if (parent) {
        nodeCircle.bottom = thisNode->in_object.v;
        nodeCircle.right = thisNode->in_object.h;
        nodeCircle.left = parent->in_object.h;
        nodeCircle.top = parent->in_object.v;
        if (isv) {
            RotateRect90(nodeCircle);
        }
        DrawStraightEdge(theCanvas,nodeCircle,thisNode);
    }
    if ((n=thisNode->get_num_nodes())) {
        if (parent&&(treeFlags&HY_TREEPANEL_INT_LABELS)) {
            if (thisNode->in_object.branchName.sLength>0) {

                if (isv) {
                    long nn = thisNode->in_object.v-thisNode->in_object.textWidth-(treeLabelFont.size/2);
                    if (nn<0) {
                        nn = 0;
                    }
                    theCanvas->DisplayText (thisNode->in_object.branchName, thisNode->in_object.h,nn,true);
                } else {
                    theCanvas->DisplayText (thisNode->in_object.branchName, thisNode->in_object.v-(treeLabelFont.size/2),thisNode->in_object.h-(treeLabelFont.size/2)-thisNode->in_object.textWidth,true);
                }
            }
        }
        for (; n; n--) {
            PaintStraight (theCanvas,thisNode->go_down(n));
        }
    } else {
        if (treeFlags&HY_TREEPANEL_TIP_LABELS) {
            if (thisNode->in_object.branchName.sLength>0) {
                if (isv) {
                    _Parameter branchSlope = (thisNode->in_object.h-thisNode->parent->in_object.h)
                                             /(_Parameter)(thisNode->parent->in_object.v-thisNode->in_object.v+1);
                    //n = -1+2*(thisNode->in_object.h>thisNode->parent->in_object.h);
                    if ((branchSlope < 1.)&&(branchSlope > -1.))
                        if (thisNode->parent->in_object.v>thisNode->in_object.v)
                            theCanvas->DisplayText (thisNode->in_object.branchName,
                                                    thisNode->in_object.h+(treeLabelFont.size/2),
                                                    thisNode->in_object.v+(treeLabelFont.size/2),
                                                    true);
                        else
                            theCanvas->DisplayText (thisNode->in_object.branchName,
                                                    thisNode->in_object.h+(treeLabelFont.size/2),
                                                    thisNode->in_object.v-(treeLabelFont.size/2)-thisNode->in_object.textWidth,
                                                    true);
                    else if (thisNode->parent->in_object.h>thisNode->in_object.h)
                        theCanvas->DisplayText (thisNode->in_object.branchName,
                                                thisNode->in_object.h-(treeLabelFont.size/2),
                                                thisNode->in_object.v-thisNode->in_object.textWidth/2,
                                                true);
                    else
                        theCanvas->DisplayText (thisNode->in_object.branchName,
                                                thisNode->in_object.h+3*(treeLabelFont.size/2),
                                                thisNode->in_object.v-thisNode->in_object.textWidth/2,
                                                true);
                } else {
                    _Parameter branchSlope = (thisNode->in_object.v-thisNode->parent->in_object.v)
                                             /(_Parameter)(thisNode->in_object.h-thisNode->parent->in_object.h+1);
                    //n = -1+2*(thisNode->in_object.h>thisNode->parent->in_object.h);
                    if ((branchSlope < 1.)&&(branchSlope > -1.))
                        if (thisNode->parent->in_object.h<=thisNode->in_object.h)
                            theCanvas->DisplayText (thisNode->in_object.branchName,
                                                    thisNode->in_object.v+(treeLabelFont.size/2),
                                                    thisNode->in_object.h+(treeLabelFont.size/2),
                                                    true);
                        else
                            theCanvas->DisplayText (thisNode->in_object.branchName,
                                                    thisNode->in_object.v+(treeLabelFont.size/2),
                                                    thisNode->in_object.h-(treeLabelFont.size/2)-thisNode->in_object.textWidth,
                                                    true);
                    else if (thisNode->parent->in_object.v>thisNode->in_object.v)
                        theCanvas->DisplayText (thisNode->in_object.branchName,
                                                thisNode->in_object.v-(treeLabelFont.size/2),
                                                thisNode->in_object.h-thisNode->in_object.textWidth/2,
                                                true);
                    else
                        theCanvas->DisplayText (thisNode->in_object.branchName,
                                                thisNode->in_object.v+3*(treeLabelFont.size/2),
                                                thisNode->in_object.h-thisNode->in_object.textWidth/2,
                                                true);
                }
            }
        }
    }
}

//__________________________________________________________

void  _HYTreePanel::DrawStraightEdge (_HYCanvas* theCanvas, _HYRect& nodeCircle, node<nodeCoord>* thisNode)
{
    if (thisNode->in_object.flags&HY_BRANCH_SELECT) {
        nodeCircle.width ++;
        theCanvas->DrawHatchedLine (nodeCircle);
        nodeCircle.width --;
    } else {
        theCanvas->DrawLine (nodeCircle);
    }
}


//__________________________________________________________

void  _HYTreePanel::PaintRadial (_HYCanvas* theCanvas, node<nodeCoord>* thisNode)
/* the idea is as follows:
    for every node draw a radial segment (back to the root of the tree) proportional to the length of the branch
    for descendants of the root, simply connect back to the root
    for every internal node, except root; draw an arc from the left-most to the right-most child
*/
{
    _HYRect                   nodeCircle;
    node<nodeCoord>* parent = thisNode->get_parent();
    nodeCircle.width        = branchWidth;
    long                      n = thisNode->get_num_nodes();

    double                    myR, // radial distance from the root
                              rootH = coordTree->in_object.h,
                              rootV = coordTree->in_object.v;
    if (parent) {
        // compute radial distance from this node to the root

        myR         = DistanceBetweenPoints (rootH,rootV,thisNode->in_object.h,thisNode->in_object.v);
        double  pR  = DistanceBetweenPoints (rootH,rootV,parent->in_object.h,parent->in_object.v);


        if (myR-pR >= branchWidth) { // draw a connector
            double  myAngle = (90.-AngleBetweenPoints (rootH,rootV,thisNode->in_object.h,thisNode->in_object.v))*pi_const/180.,
                    cosMA   = cos(myAngle);

            myAngle = sin(myAngle);

            nodeCircle.left   = rootH+cosMA*pR;
            nodeCircle.right  = rootH+cosMA*myR;

            nodeCircle.top    = rootV-myAngle*pR;
            nodeCircle.bottom = rootV-myAngle*myR;
            DrawStraightEdge (theCanvas,nodeCircle,thisNode);
        }

        //theCanvas->DisplayText (thisNode->in_object.branchName,thisNode->in_object.v,thisNode->in_object.h,true);

        if (n>1) // has > 1 children
            // draw the arc; angles measured in degrees clockwise from standard 90 degrees (12 o'clock)
        {
            parent             = thisNode->go_down (1);
            long               angle1 = AngleBetweenPoints (rootH,rootV,parent->in_object.h,parent->in_object.v);
            parent             = thisNode->go_down (n);
            long               angle2 = AngleBetweenPoints (rootH,rootV,parent->in_object.h,parent->in_object.v);

            nodeCircle.left    = rootH - myR;
            nodeCircle.right   = rootH + myR;
            nodeCircle.top     = rootV - myR;
            nodeCircle.bottom  = rootV + myR;
            angle1             = (angle1-angle2>0)?(angle1-angle2):(angle1-angle2+360);

            if (angle1 <= 10) { // draw a straight line

            } else { // draw arc
                theCanvas->DrawArc (nodeCircle,angle2,angle1);
            }
            //theCanvas->SetColor (black);
        }

    }
    // do nothing for the root node

    // paint children
    for (; n; n--) {
        PaintRadial (theCanvas,thisNode->go_down(n));
    }
}

//__________________________________________________________

void  _HYTreePanel::PaintNRadial (_HYCanvas* theCanvas, node<nodeCoord>* thisNode, _Parameter hScale, _Parameter vScale)
{
    PaintNStraight (theCanvas, thisNode, hScale, vScale);
}


//__________________________________________________________

void  _HYTreePanel::PaintNStraight (_HYCanvas* theCanvas, node<nodeCoord>* thisNode, _Parameter hScale, _Parameter vScale)
{
    _HYRect nodeCircle;
    node<nodeCoord>* parent = thisNode->get_parent();
    _String nameString;
    nodeCircle.width = 1;
    if (parent) {
        nodeCircle.bottom = thisNode->in_object.v*vScale+HY_TREEPANEL_NAVSPACING;
        nodeCircle.right = thisNode->in_object.h*hScale+HY_TREEPANEL_NAVSPACING;
        nodeCircle.left = parent->in_object.h*hScale+HY_TREEPANEL_NAVSPACING;
        nodeCircle.top = parent->in_object.v*vScale+HY_TREEPANEL_NAVSPACING;
        if (treeFlags&HY_TREEPANEL_VERTICAL) {
            RotateRect90(nodeCircle);
        }
        theCanvas->DrawLine (nodeCircle);
    }
    int n;
    if ((n=thisNode->get_num_nodes()))
        for (; n; n--) {
            PaintNStraight (theCanvas,thisNode->go_down(n),hScale,vScale);
        }
}

//__________________________________________________________

void  _HYTreePanel::PaintArcs (_HYCanvas* theCanvas, node<nodeCoord>* thisNode)
{
    _HYRect nodeCircle, circleRect;
    node<nodeCoord>* parent = thisNode->get_parent();
    int     h,w,s,f;
    bool    isv = IsVertical();
    nodeCircle.width = branchWidth;
    if (parent) {
        nodeCircle.right = thisNode->in_object.h;
        nodeCircle.left = parent->in_object.h;
        nodeCircle.top = parent->in_object.v;
        nodeCircle.bottom = thisNode->in_object.v;
        w = nodeCircle.right-nodeCircle.left;
        nodeCircle.right += w;
        if (nodeCircle.top<nodeCircle.bottom) {
            h=nodeCircle.bottom-nodeCircle.top;
            nodeCircle.top-=h;
            s=f=-90;
            if (isv) {
                s=90;
            }
        } else {
            nodeCircle.bottom=nodeCircle.top;
            nodeCircle.top=thisNode->in_object.v;
            h=nodeCircle.bottom-nodeCircle.top;
            nodeCircle.bottom+=h;
            s=270;
            f=90;
        }
        if (h<w) {
            circleRect = nodeCircle;
            circleRect.right = circleRect.left+2*h;
            if (isv) {
                RotateRect90(circleRect);
            }
            theCanvas->DrawArc(circleRect,s,f);
            if (isv) {
                RotateRect90(circleRect);
            }
            circleRect.bottom = circleRect.top = thisNode->in_object.v-1;
            circleRect.right = circleRect.left+w;
            circleRect.left += h;
            if (isv) {
                RotateRect90(circleRect);
            }
            theCanvas->DrawLine(circleRect);
        } else {
            if (isv) {
                RotateRect90(nodeCircle);
            }
            theCanvas->DrawArc(nodeCircle,s,f);
        }

    }
    if ((w=thisNode->get_num_nodes())) {
        if (parent&&(treeFlags&HY_TREEPANEL_INT_LABELS)) {
            if (thisNode->in_object.branchName.sLength) {
                if (isv) {
                    theCanvas->DisplayText (thisNode->in_object.branchName, thisNode->in_object.h+(treeLabelFont.size/2),thisNode->in_object.h+(treeLabelFont.size/2),true);
                } else {
                    theCanvas->DisplayText (thisNode->in_object.branchName, thisNode->in_object.v+(treeLabelFont.size/2),thisNode->in_object.h+(treeLabelFont.size/2),true);
                }
            }
        }
        for (; w; w--) {
            PaintArcs (theCanvas,thisNode->go_down(w));
        }
    } else {
        if (treeFlags&HY_TREEPANEL_TIP_LABELS) {
            if (thisNode->in_object.branchName.sLength) {
                if (isv) {
                    w = thisNode->in_object.v-thisNode->in_object.textWidth/2;
                    if (w<0) {
                        w = 0;
                    }
                    theCanvas->DisplayText (thisNode->in_object.branchName, thisNode->in_object.h+3*(treeLabelFont.size/2),w,true);
                } else {
                    theCanvas->DisplayText (thisNode->in_object.branchName, thisNode->in_object.v+(treeLabelFont.size/2),thisNode->in_object.h+(treeLabelFont.size/2),true);
                }
            }
        }
    }
}

//__________________________________________________________

void _HYTreePanel::PaintNArcs (_HYCanvas* theCanvas, node<nodeCoord>* thisNode, _Parameter hScale, _Parameter vScale)
{
    _HYRect nodeCircle, circleRect;
    node<nodeCoord>* parent = thisNode->get_parent();
    _String nameString;
    nodeCircle.width = 1;
    int     h,w,s,f;
    bool    isv = IsVertical();
    if (parent) {
        nodeCircle.right = thisNode->in_object.h*hScale+HY_TREEPANEL_NAVSPACING;
        nodeCircle.left = parent->in_object.h*hScale+HY_TREEPANEL_NAVSPACING;
        w = nodeCircle.right-nodeCircle.left;
        nodeCircle.right += w;
        nodeCircle.top = parent->in_object.v*vScale+HY_TREEPANEL_NAVSPACING;
        nodeCircle.bottom = thisNode->in_object.v*vScale+HY_TREEPANEL_NAVSPACING;
        if (nodeCircle.top<nodeCircle.bottom) {
            h=nodeCircle.bottom-nodeCircle.top;
            nodeCircle.top-=h;
            s=f=-90;
            if (isv) {
                s = 90;
            }
        } else {
            nodeCircle.bottom=nodeCircle.top;
            nodeCircle.top=thisNode->in_object.v*vScale+HY_TREEPANEL_NAVSPACING;
            h=nodeCircle.bottom-nodeCircle.top;
            nodeCircle.bottom+=h;
            s=270;
            f=90;
        }
        if (h<w) {
            circleRect = nodeCircle;
            circleRect.right = circleRect.left+2*h;
            if (isv) {
                RotateRect90(circleRect);
            }
            theCanvas->DrawArc(circleRect,s,f);
            if (isv) {
                RotateRect90(circleRect);
            }
            circleRect.bottom = circleRect.top = thisNode->in_object.v*vScale+HY_TREEPANEL_NAVSPACING;
            circleRect.right = circleRect.left+w;
            circleRect.left += h;
            if (isv) {
                RotateRect90(circleRect);
            }
            theCanvas->DrawLine(circleRect);
        } else {
            if (isv) {
                RotateRect90(nodeCircle);
            }
            theCanvas->DrawArc(nodeCircle,s,f);
        }
    }
    if ((w=thisNode->get_num_nodes())) {
        for (; w; w--) {
            PaintNArcs (theCanvas,thisNode->go_down(w),hScale,vScale);
        }
    }

}

//__________________________________________________________

_HYRect  _HYTreePanel::ComputeNavRect (void)
{
    _HYRect res;
    _HYCanvas* bigCanvas = (_HYCanvas*)GetObject (0);
    res.width = 1;
    int     nS = HY_TREEPANEL_NAVSIZE-HY_TREEPANEL_NAVSPACING;

    _Parameter visProp = bigCanvas->GetHSize()/(_Parameter)bigCanvas->_HYComponent::GetMaxW();
    res.left = bigCanvas->hOrigin*nS/(_Parameter)bigCanvas->_HYComponent::GetMaxW();
    res.right = res.left+visProp*nS;
    visProp = bigCanvas->GetVSize()/(_Parameter)bigCanvas->_HYComponent::GetMaxH();
    res.top = bigCanvas->vOrigin*nS/(_Parameter)bigCanvas->_HYComponent::GetMaxH();
    res.bottom = res.top+visProp*nS;
    return res;
}

//__________________________________________________________

void  _HYTreePanel::SetNavRectCenter (int h, int v)
{
    _HYCanvas * bigCanvas = (_HYCanvas*)GetObject (0);
    _HYCanvas * theRuler = (_HYCanvas*)GetObject (2);
    _HYCanvas * theNP = (_HYCanvas*)GetObject (1);

    int     nS = HY_TREEPANEL_NAVSIZE-HY_TREEPANEL_NAVSPACING-1,t;
    t = navRect.right-navRect.left;
    h-=t/2;
    if (h<0) {
        h=0;
    }
    if (h+t>nS) {
        h=nS-t;
    }
    t = navRect.bottom-navRect.top;
    v-=t/2;
    if (v<0) {
        v=0;
    }
    if (v+t>nS) {
        v=nS-t;
    }
    h = (h*bigCanvas->GetMaxW())/nS;
    v = (v*bigCanvas->GetMaxH())/nS;
    if (v<0) {
        v = 0;
    }
    if (h<0) {
        h = 0;
    }
    if ((h!=bigCanvas->hOrigin)||(v!=bigCanvas->vOrigin)) {
        bigCanvas->SetOrigin (v,h);
        theRuler->SetOrigin (0,h);
        theNP->_MarkForUpdate();
        Update(nil);
    }
}

//__________________________________________________________

void  _HYTreePanel::SetFont (_HYFont& newFont)
{
    treeLabelFont = newFont;
    if ((treeFlags&HY_TREEPANEL_TIP_LABELS)||(treeFlags&HY_TREEPANEL_INT_LABELS)||(treeFlags&HY_TREEPANEL_LABEL1)||(treeFlags&HY_TREEPANEL_LABEL2)) {
        node<nodeCoord>* currentNd = NodeTraverser(coordTree);
        textSpace = 0;
        while (currentNd) {
            if (currentNd->in_object.branchName.sLength) {
                currentNd->in_object.textWidth = GetVisibleStringWidth (currentNd->in_object.branchName,treeLabelFont);
            }

            if (currentNd->in_object.varRef>=0) {
                long f = GetVisibleStringWidth (currentNd->in_object.branchName,treeLabelFont);
                currentNd->in_object.textWidth = f;
                if ((currentNd->get_num_nodes()==0)&&(f>textSpace)) {
                    textSpace = f;
                }
            } else {
                currentNd->in_object.textWidth = 0;
                currentNd->in_object.branchName = empty;
            }
            currentNd = NodeTraverser((node<nodeCoord>*)NULL);
        }

        if (treeFlags & HY_TREEPANEL_TIP_LABELS) {
            SetHSpace (hSpace,true);
        }

        if (treeLabelFont.size/2+5>windowTextMarginV) {
            windowTextMarginV = treeLabelFont.size/2+5;
            SetVSpace (vSpace,true);
            return;
        } else {
            if (treeLabelFont.size/2+5<windowTextMarginV) {
                int newVSpace = treeLabelFont.size/2+5;
                if (((newVSpace>=branchLabel1.size+3)||(!(treeFlags&HY_TREEPANEL_LABEL1)))&&
                        ((newVSpace>=branchLabel2.size+3)||(!(treeFlags&HY_TREEPANEL_LABEL2)))) {
                    windowTextMarginV = newVSpace;
                    SetVSpace (vSpace,true);
                    return;
                }
            }
        }
        RenderTree();
    }
}

//__________________________________________________________

_String     _HYTreePanel::GetTreeString (_AVLListXL * remap)
{
    _Variable* theTree = FetchVar(LocateVarByName (treeName));
    if (theTree&&theTree->ObjectClass()==TREE) {
        _TheTree * myTree = (_TheTree*)theTree;
        myTree->StepWiseTraversal (true);
        _String  res ((unsigned long)16,true);
        long flag = -1;
        if (scaleVariable.sLength) {
            if (scaleVariable == expectedNumberOfSubs) {
                flag = -3;
            } else if (scaleVariable == stringSuppliedLengths) {
                flag = -2;
            } else
#ifndef USE_AVL_NAMES
                flag = variableReindex.lData[LocateVarByName(scaleVariable)];
#else
            {
                _Variable* scaleVar = CheckReceptacle  (&scaleVariable,empty,false);
                //flag = variableNames.GetXtra (LocateVarByName(scaleVariable));
                flag = scaleVar->GetAVariable();
            }
#endif
        }
        if (myTree->rooted == UNROOTED || myTree->theRoot->get_num_nodes() != 3) {
            myTree->SubTreeString (res,treeFlags&HY_TREEPANEL_INT_LABELS,flag, remap);
        } else {
            if (myTree->rooted == ROOTED_LEFT) {
                res << "((";
                myTree->currentNode = myTree->theRoot->go_down(1);
                myTree->SubTreeString (res,treeFlags&HY_TREEPANEL_INT_LABELS,flag, remap);
                res << ",";
                myTree->currentNode = myTree->theRoot->go_down(2);
                myTree->SubTreeString (res,treeFlags&HY_TREEPANEL_INT_LABELS,flag, remap);
                res << "),";
                myTree->currentNode = myTree->theRoot->go_down(3);
                myTree->SubTreeString (res,treeFlags&HY_TREEPANEL_INT_LABELS,flag, remap);
                res << ")";
            } else {
                res << "(";
                myTree->currentNode = myTree->theRoot->go_down(1);
                myTree->SubTreeString (res,treeFlags&HY_TREEPANEL_INT_LABELS,flag, remap);
                res << ",(";
                myTree->currentNode = myTree->theRoot->go_down(2);
                myTree->SubTreeString (res,treeFlags&HY_TREEPANEL_INT_LABELS,flag, remap);
                res << ",";
                myTree->currentNode = myTree->theRoot->go_down(3);
                myTree->SubTreeString (res,treeFlags&HY_TREEPANEL_INT_LABELS,flag, remap);
                res << "))";
            }
        }
        res.Finalize();
        return res;
    }
    return empty;
}

//__________________________________________________________

void    _HYTreePanel::ShowModelMatrix (bool expMx)
{
    node<nodeCoord>* node1 = (node<nodeCoord>*)currentSelection(0);
    if (node1->in_object.varRef>=0) {
        _CalcNode * thisCNode = (_CalcNode*)LocateVar (node1->in_object.varRef);
        OpenModelFromCalcNode (*thisCNode->GetName(),expMx);
    }
}
//__________________________________________________________

void    _HYTreePanel::SwapSelectedSubTrees (void)
{
    _String errMsg ("Internal Tree Error in SwapSelectedSubTrees");
    _TheTree *me = LocateMyTreeVariable();
    if (!me) {
        return;
    }
    node<nodeCoord>* node1 = (node<nodeCoord>*)currentSelection(0),
                     * node2 = (node<nodeCoord>*)currentSelection(1),
                       *t1 = node1->parent, *t2 = node2->parent;
    long i,j;

    for (i=0; i<t1->nodes.length; i++) {
        if (t1->nodes.data[i] == node1) {
            break;
        }
    }
    for (j=0; j<t2->nodes.length; j++) {
        if (t2->nodes.data[j] == node2) {
            break;
        }
    }
    if ((i==t1->nodes.length)||(j==t2->nodes.length)) {
        WarnError (errMsg);
        return;
    }
    node1->parent->nodes.data[i] = node2;
    node2->parent->nodes.data[j] = node1;
    node1->parent = t2;
    node2->parent = t1;
    node<long> *n1 = nil,*n2 = nil,*p1 =nil,*p2 = nil;

    _CalcNode* travNode = me->DepthWiseTraversal (true);
    while ((travNode)&&(!(n1&&n2))) {
        if (travNode->theIndex==node1->in_object.varRef) {
            n1 = &me->GetCurrentNode();
        } else if (travNode->theIndex==node2->in_object.varRef) {
            n2 = &me->GetCurrentNode();
        }
        travNode = me->DepthWiseTraversal();
    }
    p1 = n1->parent;
    p2 = n2->parent;
    if (!(n1&&n2&&p1&&p2)) {
        WarnError (errMsg);
        return;
    }
    for (i=0; i<p1->nodes.length; i++) {
        if (p1->nodes.data[i] == n1) {
            break;
        }
    }
    for (j=0; j<p2->nodes.length; j++) {
        if (p2->nodes.data[j] == n2) {
            break;
        }
    }
    if ((i==p1->nodes.length)||(j==p2->nodes.length)) {
        WarnError (errMsg);
        return;
    }
    n1->parent->nodes.data[i] = n2;
    n2->parent->nodes.data[j] = n1;
    n1->parent = p2;
    n2->parent = p1;
    if (t1!=t2) {
        BuildTree(true);
        node1 = NodeTraverser (coordTree);
        while (node1) {
            if ((node1->in_object.varRef==n1->in_object)||
                    (node1->in_object.varRef==n2->in_object)) {
                currentSelection << (long)(node1);
                node1->in_object.flags|=HY_BRANCH_SELECT;
            }
            node1 = NodeTraverser((node <nodeCoord>*)nil);
        }
    }
    PrepareUndoData (1);

    if (likeFuncID>=0) {
        _LikelihoodFunction* theLF = (_LikelihoodFunction*)likeFuncList(likeFuncID);
        theLF->MapTreeTipsToData (theLF->DependOnTree(*me->GetName()),true);
        theLF->VoidOldResults();
        me->SetUp();
        DisplayLikelihoodValue();
    }
    RenderTree((treeFlags&HY_TREEPANEL_PROJECTION)==0);
}

//__________________________________________________________

void    _HYTreePanel::FlipSelectedBranches (void)
{
    _String errMsg ("Internal Tree Error in FlipSelectedBranches");
    _TheTree *me = LocateMyTreeVariable();
    if (!me) {
        return;
    }
    node<nodeCoord>* node1 = (node<nodeCoord>*)currentSelection(0);
    long i,k;
    _CalcNode* travNode = me->DepthWiseTraversal (true);
    PrepareUndoData (2);

    node<long>*      n1 = nil;

    while (travNode) {
        if (travNode->theIndex==node1->in_object.varRef) {
            n1 = &me->GetCurrentNode();
            break;
        }
        travNode = me->DepthWiseTraversal();
    }

    if (!n1) {
        WarnError (errMsg);
        return;
    }

    i = 0;

    k = node1->nodes.length-1;

    while (i<k) {
        node<long>*     swapL = n1->nodes.data[i];
        n1->nodes.data[i] = n1->nodes.data[k];
        n1->nodes.data[k] = swapL;

        node<nodeCoord>*    swapNC = node1->nodes.data[i];
        node1->nodes.data[i] = node1->nodes.data[k];
        node1->nodes.data[k] = swapNC;
        i++;
        k--;
    }

    if (likeFuncID>=0) {
        _LikelihoodFunction* theLF = (_LikelihoodFunction*)likeFuncList(likeFuncID);
        theLF->MapTreeTipsToData (theLF->DependOnTree(*me->GetName()),true);
        theLF->VoidOldResults();
        me->SetUp();
        DisplayLikelihoodValue();
    }
    RenderTree((treeFlags&HY_TREEPANEL_PROJECTION)==0);
}

//__________________________________________________________

void    _HYTreePanel::JoinSelectedBranches (void)
{
    _String errMsg ("Internal Tree Error in CollapseSelectedBranch");
    _TheTree *me = LocateMyTreeVariable();
    if (!me) {
        return;
    }
    if (!CheckIfNeedUnscaled()) {
        return;
    }
    node<nodeCoord>* node1=(node<nodeCoord>*)currentSelection(0);
    node<long> *n1 = nil,*p1 =nil,*newp;
    PrepareUndoData (5);
    long k;
    _CalcNode* travNode = me->DepthWiseTraversal (true);
    while (travNode) {
        if (travNode->theIndex==node1->in_object.varRef) {
            n1 = &me->GetCurrentNode();
            break;
        }
        travNode = me->DepthWiseTraversal();
    }
    if (n1) {
        p1 = n1->parent;
    }
    if (!(n1&&p1)) {
        WarnError (errMsg);
        return;
    }
    newp = new node<long>;
    checkPointer ((Ptr)newp);
    for (k=0; k<currentSelection.lLength; k++) {
        node1 = (node<nodeCoord>*)currentSelection.lData[k];
        for (long j = 0; j<p1->nodes.length; j++) {
            if (p1->nodes.data[j]) {
                if (p1->nodes.data[j]->in_object==node1->in_object.varRef) {
                    newp->add_node(*p1->nodes.data[j]);
                    p1->nodes.data[j] = nil;
                    break;
                }
            }
        }
    }
    long shift = 0, l = p1->nodes.length;
    newp->set_parent (*p1);
    for (k=0; k<l; k++) {
        if (shift) {
            if (p1->nodes.data[k]) {
                p1->nodes.data[k-shift+1]=p1->nodes.data[k];
            } else {
                shift++;
            }
        } else {
            if (!p1->nodes.data[k]) {
                p1->nodes.data[k]=newp;
                shift = 1;
            }
        }
    }
    p1->nodes.length-=shift-1;
    shift = travNode->GetModelIndex();
    _String nodeName (empty);
    _String dummy ("1.0"), dummy2 ("");
    if (shift != HY_NO_MODEL) {
        nodeName = *(_String*)modelNames(shift);
        me->FinalizeNode(newp,FindUnusedSuffix ("Node"),dummy2,nodeName,dummy);
        _CalcNode* targetNode = (_CalcNode*)LocateVar (newp->in_object);
        targetNode->SetCodeBase (me->GetCodeBase());
        targetNode->CopyMatrixParameters(travNode);
        dubiousNodes<<targetNode->GetAVariable();
        UpdateLnLikelihood();
    } else {
        me->FinalizeNode(newp,FindUnusedSuffix ("Node"),nodeName,dummy2,dummy);
    }
    if (likeFuncID>=0) {
        _LikelihoodFunction* theLF = (_LikelihoodFunction*)likeFuncList(likeFuncID);
        theLF->MapTreeTipsToData (theLF->DependOnTree(*me->GetName()),true);
    }
    BuildTree(true);

    RenderTree();
}

//__________________________________________________________

void    _HYTreePanel::GraftATip (void)
{
    _String errMsg ("Internal Tree Error in GraftATip");
    _TheTree *me = LocateMyTreeVariable();
    if (!me) {
        return;
    }
    if (!CheckIfNeedUnscaled()) {
        return;
    }
    if (!KillLikeFunc(7)) {
        return;
    }

    node<long> *n1 = nil,*p1 =nil,*newt,*newp;
    _SimpleList matchingIDs,matchingNodes;
    long k,f;
    for (k=0; k<currentSelection.lLength; k++) {
        matchingIDs<<((node<nodeCoord>*)currentSelection(k))->in_object.varRef;
    }
    _CalcNode* travNode = me->DepthWiseTraversal (true);
    while (travNode) {
        f = matchingIDs.Find(travNode->theIndex);
        if (f>=0) {
            matchingNodes<<(long)&me->GetCurrentNode();
            matchingIDs.Delete(f);
        }
        travNode = me->DepthWiseTraversal();
    }
    for (k=0; k<matchingNodes.lLength; k++) {
        n1 = (node<long>*)matchingNodes.lData[k];
        p1 = n1->parent;
        if (!p1) {
            warnError (errMsg);
            return;
        }
        newp = new node<long>;
        newt = new node<long>;
        checkPointer ((Ptr)newp);
        checkPointer ((Ptr)newt);
        newp->set_parent(*p1);
        newp->add_node(*newt);
        newp->add_node(*n1);
        for (long j = 0; j<p1->nodes.length; j++) {
            if (p1->nodes.data[j]->in_object==n1->in_object) {
                p1->nodes.data[j] = newp;
                break;
            }
        }
        travNode = (_CalcNode*)LocateVar (n1->in_object);
        _String dummy ("1.0"), dummy2, dummy3;
        me->FinalizeNode(newp,FindUnusedSuffix ("Node"),dummy2,dummy3,dummy);
        _String nodeName ("Species");
        nodeName = nodeName&_String (FindUnusedSuffix(nodeName));
        dummy2 = empty;
        me->FinalizeNode(newt,0,nodeName,dummy2,dummy);
        ((_VariableContainer*)LocateVar(newp->in_object))->CopyMatrixParameters (travNode);
        ((_VariableContainer*)LocateVar(newt->in_object))->CopyMatrixParameters (travNode);
        travNode = (_CalcNode*)LocateVar (newp->in_object);
        travNode->SetCodeBase (me->GetCodeBase());
        travNode = (_CalcNode*)LocateVar (newt->in_object);
        travNode->SetCodeBase (me->GetCodeBase());
    }

    BuildTree(true);

    RenderTree();
}

//__________________________________________________________

void    _HYTreePanel::RecalculateLikelihood (void)
{
    if (likeFuncID>=0) {
        _LikelihoodFunction* lf = (_LikelihoodFunction*)likeFuncList (likeFuncID);
        feedbackTreePanel = this;
        ToggleAnalysisMenu (true);
        StartBarTimer();
        _Matrix* res = lf->Optimize();
        StopBarTimer();
        ToggleAnalysisMenu (false);
        terminateExecution = false;
        feedbackTreePanel = nil;
        DeleteObject (res);
        dubiousNodes.Clear();
        FlushUndoData();
        BuildTree(true);
        RenderTree();
        postChangeLFEvent (GetID(),likeFuncID);
    }
}
//__________________________________________________________

void    _HYTreePanel::RerootTree (void)
{
    _String errMsg ("Internal Tree Error in RerootTree.");
    _TheTree *me = LocateMyTreeVariable();
    if (!me) {
        return;
    }
    if (currentSelection.lLength) {
        node<nodeCoord>* node1;
        node<long> *n1 = nil,*p1 =nil,*t1, *t2, *t3;
        node1 = (node<nodeCoord>*)currentSelection(0);
        long        i,j;
        _CalcNode* travNode = me->DepthWiseTraversal (true);
        FlushUndoData();
        unsigned char newRFlag;
        while (travNode) {
            if (travNode->theIndex==node1->in_object.varRef) {
                n1 = &me->GetCurrentNode();
                break;
            }
            travNode = me->DepthWiseTraversal();
        }
        if (n1) {
            p1 = n1->parent;
        }
        if (!(n1&&p1)) {
            WarnError (errMsg);
            return;
        }

        for (i=0; i<p1->nodes.length; i++) {
            if (p1->nodes.data[i]==n1) {
                break;
            }
        }
        t1 = p1;
        t2 = t1->parent;
        if (t2) {
            newRFlag = ROOTED_RIGHT;
            node <long>* newRoot = new node<long>;
            checkPointer ((Ptr)newRoot);
            newRoot->add_node (*n1);
            for (j=0; j<p1->nodes.length; j++) {
                if (i==j) {
                    continue;
                }
                newRoot->add_node(*(p1->nodes.data[j]));
            }

            newRoot->add_node(*p1);
            while(t2) {
                t1->nodes.length = 0;
                t3 = t2->parent;
                for (i=0; i<t2->nodes.length; i++) {
                    if (t2->nodes.data[i]==t1) {
                        if (t3) {
                            t1->add_node(*t2);
                        }
                    } else {
                        t1->add_node(*(t2->nodes.data[i]));
                    }
                }
                t1 = t2;
                t2 = t3;
            }
            t3 = &me->GetRoot();
            newRoot->in_object = t3->in_object;
            delete (t3);
            me->SetRoot (newRoot);
        } else {
            newRFlag = (i>=1)?ROOTED_LEFT:ROOTED_RIGHT;
        }

        me->RootedFlag() = newRFlag;

        currentSelection.Clear();
        if (likeFuncID>=0) {
            _LikelihoodFunction* theLF = (_LikelihoodFunction*)likeFuncList(likeFuncID);
            theLF->MapTreeTipsToData (theLF->DependOnTree(*me->GetName()),true);
        }
    } else {
        if (me->RootedFlag()!=UNROOTED) {
            me->RootedFlag() = UNROOTED;
        } else {
            long       leafCount = 0;
            _Parameter ratioAB   = 1.e300,
                       t;

            node       <nodeCoord>*   newRoot = nil,
                                      *   trav    = NodeTraverser (coordTree);

            while (trav) {
                if (trav->nodes.length==0) {
                    trav->in_object.auxL = 1;
                    leafCount ++;
                } else {
                    trav->in_object.auxL = 0;
                    for (long k=0; k<trav->nodes.length; k++) {
                        trav->in_object.auxL += trav->nodes.data[k]->in_object.auxL;
                    }
                }
                trav = NodeTraverser((node <nodeCoord>*)nil);
            }

            trav    = NodeTraverser (coordTree);
            while (trav) {
                t = fabs(leafCount/(_Parameter)trav->in_object.auxL-2);
                if (t<ratioAB) {
                    ratioAB = t;
                    newRoot = trav;
                }
                trav = NodeTraverser((node <nodeCoord>*)nil);
            }
            if (newRoot->parent) {
                currentSelection << (long)newRoot;
                RerootTree();
            }
            return;
        }
    }
    //UpdateLnLikelihood ();
    me->SetUp();
    BuildTree(true);
    RenderTree();
}

//__________________________________________________________

void    _HYTreePanel::CollapseSelectedBranch (void)
{
    _String errMsg ("Internal Tree Error in CollapseSelectedBranch");
    _TheTree *me = LocateMyTreeVariable();
    if (!me) {
        return;
    }
    node<nodeCoord>* node1;
    long i,j;
    PrepareUndoData (3);
    for (long k=0; k<currentSelection.lLength; k++) {
        node1 = (node<nodeCoord>*)currentSelection(k);
        if (node1->get_num_nodes()==0) {
            continue;
        }

        node<long> *n1 = nil,*p1 =nil;

        _CalcNode* travNode = me->DepthWiseTraversal (true);
        while (travNode) {
            if (travNode->theIndex==node1->in_object.varRef) {
                n1 = &me->GetCurrentNode();
                break;
            }

            travNode = me->DepthWiseTraversal();
        }
        if (n1) {
            p1 = n1->parent;
        }
        if (!(n1&&p1)) {
            WarnError (errMsg);
            return;
        }
        for (i=0; i<p1->nodes.length; i++) {
            if (p1->nodes.data[i] == n1) {
                break;
            }
        }
        for (j=i+1; j<p1->nodes.length; j++) {
            p1->nodes.data[j-1]=p1->nodes.data[j];
        }
        p1->nodes.length--;
        if (i==0) {
            for (i=n1->nodes.length; i; i--) {
                p1->prepend_node (*n1->go_down(i));
            }
        } else {
            for (i=1; i<=n1->nodes.length; i++) {
                p1->add_node (*n1->go_down(i));
            }
        }
        delete (n1);
        DeleteVariable (*travNode->GetName(),true);
    }

    me->SetUp();
    UpdateLnLikelihood ();
    BuildTree(true);

    RenderTree();
}
//__________________________________________________________

void    _HYTreePanel::UpdateLnLikelihood (void)
{
    if (likeFuncID>=0) {
        _LikelihoodFunction* theLF = (_LikelihoodFunction*)likeFuncList (likeFuncID);
        theLF->RescanAllVariables();
        theLF->HasBeenOptimized () = false;
        postChangeLFParamsEvent(GetID(),likeFuncID);
    }
}

//__________________________________________________________

void    _HYTreePanel::Relabel (_TheTree* mt)
{
    if ((treeFlags & HY_TREEPANEL_LABEL1)&&(branchVar1.sLength)) {
        mt->AssignLabelsToBranches (coordTree,&branchVar1,false);
    }
    if ((treeFlags & HY_TREEPANEL_LABEL2)&&(branchVar2.sLength)) {
        mt->AssignLabelsToBranches (coordTree,&branchVar2,true);
    }
}

//__________________________________________________________

void    _HYTreePanel::DisplayLikelihoodValue (void)
{
    if (likeFuncID>=0) {
        _LikelihoodFunction* theLF = (_LikelihoodFunction*)likeFuncList(likeFuncID);
        statusBar.Trim(0,statusBar.FindBackwards(". Ln-likelihood",0,-1)-1);
        theLF->PrepareToCompute();
        _String newStatusBar = statusBar&". Ln-likelihood = "&_String(theLF->Compute());
        theLF->DoneComputing();
        SetStatusBar(newStatusBar);
    }
}

//__________________________________________________________

void    _HYTreePanel::DeleteCurrentSelection (void)
{
    _String errMsg ("Internal Tree Error in DeleteCurrentSelection");
    _String* lastName;
    _TheTree *me = LocateMyTreeVariable();
    if (!me) {
        return;
    }
    if (!ShrinkDataFilter()) {
        return;
    }
    node<nodeCoord>* node1;
    node<long>* n1, *p1, *p2;
    PrepareUndoData (4);
    long i,j,k;
    _SimpleList      selectionID;
    for (k=0; k<currentSelection.lLength; k++) {
        node1 = (node<nodeCoord>*) currentSelection(k);
        selectionID<<node1->in_object.varRef;
    }
    _CalcNode* travNode = me->DepthWiseTraversal (true);
    while (travNode) {
        k = selectionID.Find (travNode->theIndex);
        if (k>=0) {
            n1 = &me->GetCurrentNode();
            selectionID.Delete(k);
            p1 = n1->parent;
            if (!(p1&&n1)) {
                warnError (errMsg);
                return;
            }
            if (n1->nodes.length==0) { // removing a tip
                if (p1->nodes.length>2) {
                    lastName = travNode->GetName();
                    travNode = me->DepthWiseTraversal();
                    for (i=0; i<p1->nodes.length; i++) {
                        if (p1->nodes.data[i]==n1) {
                            break;
                        }
                    }
                    for (j=i+1; j<p1->nodes.length; j++) {
                        p1->nodes.data[j-1]=p1->nodes.data[j];
                    }
                    p1->nodes.length--;
                    DeleteVariable(*lastName,true);
                    delete (n1);
                    continue;
                } else {
                    p2 = p1->parent;
                    if (p2) {
                        for (i=0; i<p2->nodes.length; i++) {
                            if (p2->nodes.data[i]==p1) {
                                break;
                            }
                        }
                        lastName = travNode->GetName();
                        travNode = me->DepthWiseTraversal();
                        DeleteVariable(*lastName,true);
                        if (p1->nodes.data[0]==n1) {
                            p2->nodes.data[i]=p1->nodes.data[1];
                            p1->nodes.data[1]->parent = p2;
                        } else {
                            lastName = travNode->GetName();
                            travNode = me->DepthWiseTraversal();
                            DeleteVariable(*lastName,true);
                            p2->nodes.data[i]=p1->nodes.data[0];
                            p1->nodes.data[0]->parent = p2;
                        }
                        delete (n1);
                        delete (p1);
                        continue;
                    }
                }
            }
        }
        travNode = me->DepthWiseTraversal ();
    }

    currentSelection.Clear();
    //likeFuncID = -1;
    ShrinkDataFilter (true);
    if (likeFuncID>=0) {
        UpdateLnLikelihood();
        me->SetUp();
    }
    BuildTree(true);

    RenderTree();
}

//__________________________________________________________

void    _HYTreePanel::CutSelectionToClipboard (bool cut)
{
    _String errMsg ("Internal Tree Error in CutSelectionToClipboard");
    _TheTree *me = LocateMyTreeVariable();
    long     i,k;
    node<long>* ct = nil;
    if (!me) {
        return;
    }
    if (cut)
        if (!ShrinkDataFilter()) {
            return;
        }
    PrepareUndoData (6);
    if (selectionTop) {
        _CalcNode* travNode = me->DepthWiseTraversal (true), *killedNode;
        while (travNode) {
            if (travNode->GetAVariable()==selectionTop->in_object.varRef) {
                ct = &me->GetCurrentNode();
                break;
            }
            travNode = me->DepthWiseTraversal();
        }

        if (ct) {
            node<long> *p1, *p2;

            p1 = ct->parent;
            if (p1) {
                if (cut) {
                    for (i=0; i<p1->nodes.length; i++) {
                        if (p1->nodes.data[i]==ct) {
                            break;
                        }
                    }
                    if (p1->nodes.length>2) {
                        for (k=i+1; k<p1->nodes.length; k++) {
                            p1->nodes.data[k-1]=p1->nodes.data[k];
                        }
                        p1->nodes.length--;
                    } else {
                        p2 = p1->parent;
                        if (!p2) {
                            WarnError (errMsg);
                            return;
                        }
                        for (k=0; k<p2->nodes.length; k++) {
                            if (p2->nodes.data[k]==p1) {
                                break;
                            }
                        }
                        i = !i;
                        p2->nodes.data[k] = p1->go_down(i+1);
                        p2->nodes.data[k]->parent = p2;
                        killedNode = (_CalcNode*)LocateVar(p1->in_object);
                        travNode = (_CalcNode*)LocateVar(p2->nodes.data[k]->in_object);

                        if (killedNode->GetModelIndex()==travNode->GetModelIndex()) {
                            _SimpleList l1, l2;
                            {
                                _AVLList la1 (&l1),
                                         la2 (&l2);

                                travNode->ScanForVariables (la1,la1);
                                killedNode->ScanForVariables (la2,la2);
                                la1.ReorderList();
                                la2.ReorderList();
                            }

                            for (k=0; k<l1.lLength; k++) {
                                if (k==l2.lLength) {
                                    break;
                                }
                                _Variable *thisVar = LocateVar(l1.lData[k]);
                                _Constant nv (thisVar->Compute()->Value()+
                                              LocateVar(l2.lData[k])->Compute()->Value());
                                thisVar->SetValue (&nv);
                            }

                        }
                        DeleteVariable (*killedNode->GetName(),true);
                        delete(p1);
                    }
                }
                ClearClipboardContents();
                bool    hasModel = true;
                if (cut) {
                    treeFlags &= 0xFF7F;
                    treePanelClipboardRoot = ct;
                    RenameNodesInClipboard (treePanelClipboardRoot);
                    currentSelection.Clear();
                    selectionTop = nil;
                    ShrinkDataFilter(true);
                    if (likeFuncID>=0) {
                        UpdateLnLikelihood();
                        me->SetUp();
                    }
                    BuildTree(true);
                    RenderTree();
                } else {
                    treePanelClipboardRoot = CopyNodesToClipboard (ct,me,hasModel);
                }
                return;
            }
        }
    }
    WarnError (errMsg);
}

//__________________________________________________________

void    RenameNodesInClipboard  (node<long>* thisNode)
{
    for (long i=0; i<thisNode->nodes.length; i++) {
        RenameNodesInClipboard (thisNode->nodes.data[i]);
    }
    _String oldName (*LocateVar(thisNode->in_object)->GetName()), newName (oldName);
    newName.Trim (newName.Find('.')+1,-1);
    newName = _String("_CLIPBOARDTREENODES_.")&newName;
    RenameVariable(&oldName,&newName);
}

//__________________________________________________________

void    RenameUndoNodes (node<long>* thisNode,_String& newTreeName)
{
    for (long i=0; i<thisNode->nodes.length; i++) {
        RenameUndoNodes (thisNode->nodes.data[i], newTreeName);
    }
    _String oldName (*LocateVar(thisNode->in_object)->GetName()), newName (oldName);
    newName.Trim (newName.Find('.'),-1);
    newName = newTreeName&newName;
    RenameVariable(&oldName,&newName);
}

//__________________________________________________________

void    KillNodesInClipboard (node<long>* thisNode)
{
    for (long i=0; i<thisNode->nodes.length; i++) {
        KillNodesInClipboard (thisNode->nodes.data[i]);
    }
    _String oldName (*LocateVar(thisNode->in_object)->GetName());
    DeleteVariable (oldName,true);
    delete (thisNode);
}

//__________________________________________________________

node<long>* PatchNodesFromClipboard (node<long>* thisNode, _TheTree* me, bool& hasModel)
{
    node<long>* newNode = new node<long>;
    checkPointer (newNode);

    for (long i=0; i<thisNode->nodes.length; i++) {
        newNode->add_node(*PatchNodesFromClipboard (thisNode->nodes.data[i],me,hasModel));
    }


    _CalcNode *sourceNode = (_CalcNode*)LocateVar(thisNode->in_object);

    long      modelID = sourceNode->GetModelIndex(),
              k=2;

    _String   nodeName (*sourceNode->GetName()),
              tryName;

    nodeName.Trim(nodeName.Find('.')+1,-1);
    nodeName = *me->GetName()&'.'&nodeName;
    tryName  = nodeName;

    while   (LocateVarByName(tryName)>=0) {
        tryName = nodeName&'_'&k;
        k++;
    }
    if (modelID != HY_NO_MODEL) {
        nodeName = *(_String*)modelNames(modelID);
    } else {
        nodeName = empty;
        hasModel = false;
    }

    modelID = lastMatrixDeclared;
    lastMatrixDeclared = -1;

    _String dummy ("1.0");
    me->FinalizeNode(newNode,0,tryName,nodeName,dummy);

    lastMatrixDeclared = modelID;

    _CalcNode* thisCNode = (_CalcNode*)LocateVar (newNode->in_object);
    thisCNode->CopyMatrixParameters (sourceNode);
    thisCNode->SetCodeBase (me->GetCodeBase());
    return newNode;
}

//__________________________________________________________

node<long>* CopyNodesToClipboard (node<long>* thisNode, _TheTree* me, bool& hasModel)
{
    node<long>* newNode = new node<long>;
    checkPointer (newNode);
    for (long i=0; i<thisNode->nodes.length; i++) {
        newNode->add_node(*CopyNodesToClipboard (thisNode->nodes.data[i],me,hasModel));
    }

    _CalcNode *sourceNode = (_CalcNode*)LocateVar(thisNode->in_object);
    long      modelID = sourceNode->GetModelIndex();
    _String   nodeName (*sourceNode->GetName()), tryName;
    nodeName.Trim(me->GetName()->sLength+1,-1);
    if (makeUndoTree) {
        nodeName = _String("_UNDOTREENODES_4_")&*me->GetName()&'.'&nodeName;
    } else {
        nodeName = _String("_CLIPBOARDTREENODES_")&nodeName;
    }
    tryName = nodeName;
    if (modelID>=0) {
        nodeName = *(_String*)modelNames(modelID);
    } else {
        nodeName = empty;
        hasModel = false;
    }
    modelID = lastMatrixDeclared;
    lastMatrixDeclared = -1;
    _String dummy ("1.0");
    me->FinalizeNode(newNode,0,tryName,nodeName,dummy);
    lastMatrixDeclared = modelID;
    _CalcNode* thisCNode = (_CalcNode*)LocateVar (newNode->in_object);
    thisCNode->CopyMatrixParameters (sourceNode);
    thisCNode->SetCodeBase (me->GetCodeBase());
    return newNode;
}

//__________________________________________________________

void    CopyNodeParameters (node<nodeCoord>* destination, node<nodeCoord>* source)
{
    int n = source->get_num_nodes();
    if (n!=destination->get_num_nodes()) {
        return;
    }
    for (; n; n--) {
        CopyNodeParameters (destination->go_down(n),source->go_down(n));
    }
    destination->in_object.label1 = source->in_object.label1;
    destination->in_object.label2 = source->in_object.label2;
    destination->in_object.branchName = source->in_object.branchName;
    destination->in_object.textWidth = source->in_object.textWidth;
    destination->in_object.flags = source->in_object.flags;

}

//__________________________________________________________

void    _HYTreePanel::PasteClipboardTree (void)
{
    _String errMsg ("Internal Tree Error in PasteClipboardTree");
    _TheTree *me = LocateMyTreeVariable();
    long     i,k;
    if (!me) {
        return;
    }
    if (!CheckIfNeedUnscaled()) {
        return;
    }
    if (!KillLikeFunc(8)) {
        return;
    }

    node<long> *p1, *ct = nil, *newp;

    _CalcNode* travNode = me->DepthWiseTraversal (true);

    node<nodeCoord>*  node1 = (node<nodeCoord>*)currentSelection(0);

    while (travNode) {
        if (travNode->GetAVariable()==node1->in_object.varRef) {
            ct = &me->GetCurrentNode();
            break;
        }
        travNode = me->DepthWiseTraversal();
    }

    if (!ct) {
        WarnError (errMsg);
        return;
    }
    p1 = ct->parent;
    for (i=0; i<p1->nodes.length; i++) {
        if (p1->nodes.data[i]==ct) {
            break;
        }
    }
    newp = new node<long>;
    checkPointer ((Ptr)newp);
    p1->nodes.data[i] = newp;
    newp->parent = p1;
    bool    hasModel = true;
    node<long>* cliptop = PatchNodesFromClipboard (treePanelClipboardRoot,me,hasModel);
    newp->add_node(* cliptop);
    newp->add_node(*ct);

    k = travNode->GetModelIndex();
    _String nodeName (empty);
    if (k>=0) {
        nodeName = *(_String*)modelNames(k);
        _String dummy ("1.0"),dummy2;
        me->FinalizeNode(newp,FindUnusedSuffix ("Node"),dummy2,nodeName,dummy);
        _CalcNode* targetNode = (_CalcNode*)LocateVar (newp->in_object);
        targetNode->SetCodeBase (me->GetCodeBase());
        targetNode->CopyMatrixParameters(travNode);

        _SimpleList varList;
        {
            _AVLList vla (&varList);
            targetNode->ScanForVariables(vla,vla);
            travNode->ScanForVariables(vla,vla);
        }
        for (k=0; k<varList.lLength; k++) {
            _Constant myValue(LocateVar(varList.lData[k])->Compute()->Value()/2.0);
            LocateVar(varList.lData[k])->SetValue (&myValue);
        }
        dubiousNodes<<targetNode->GetAVariable();
        UpdateLnLikelihood();
    } else {
        _String dummy ("1.0"), dummy2;
        me->FinalizeNode(newp,FindUnusedSuffix ("Node"),dummy2,nodeName,dummy);
    }

    if (likeFuncID>=0) {
        _LikelihoodFunction* theLF = (_LikelihoodFunction*)likeFuncList(likeFuncID);
        theLF->MapTreeTipsToData (theLF->DependOnTree(*me->GetName()),true);
    }
    BuildTree(true);

    RenderTree();
    return;
}

//__________________________________________________________

void    _HYTreePanel::MoveSubTree (void)
{
    _String errMsg ("Internal Tree Error in MoveSubTree");
    bool    done = true,
            savePDM = promptDefModel;
    long     i,k;

    _TheTree *me = LocateMyTreeVariable();

    if (!me) {
        return;
    }

    promptDefModel = true;
    if (!CheckIfNeedUnscaled()) {
        return;
    }
    promptDefModel = savePDM;

    if (selectionTop) {
        PrepareUndoData (9);
        _CalcNode* travNode = me->DepthWiseTraversal (true), *killedNode;
        node<long>* ct = nil;

        while (travNode) {
            if (travNode->GetAVariable()==selectionTop->in_object.varRef) {
                ct = &me->GetCurrentNode();
                break;
            }
            travNode = me->DepthWiseTraversal();
        }

        if (ct) {
            node<long> *p1, *p2, *ct2 = nil;

            p1 = ct->parent;
            if (p1) {
                for (i=0; i<p1->nodes.length; i++) {
                    if (p1->nodes.data[i]==ct) {
                        break;
                    }
                }
                if (p1->nodes.length>2) {
                    for (k=i+1; k<p1->nodes.length; k++) {
                        p1->nodes.data[k-1]=p1->nodes.data[k];
                    }
                    p1->nodes.length--;
                } else {
                    p2 = p1->parent;
                    if (!p2) {
                        WarnError (errMsg);
                        return;
                    }
                    for (k=0; k<p2->nodes.length; k++) {
                        if (p2->nodes.data[k]==p1) {
                            break;
                        }
                    }
                    i = !i;
                    p2->nodes.data[k] = p1->go_down(i+1);
                    p2->nodes.data[k]->parent = p2;
                    killedNode = (_CalcNode*)LocateVar(p1->in_object);
                    travNode = (_CalcNode*)LocateVar(p2->nodes.data[k]->in_object);

                    if (killedNode->GetModelIndex()==travNode->GetModelIndex()) {
                        _SimpleList l1, l2;
                        {
                            _AVLList la1 (&l1),
                                     la2 (&l2);

                            travNode->ScanForVariables (la1,la1);
                            killedNode->ScanForVariables (la2,la2);
                            la1.ReorderList();
                            la2.ReorderList();
                        }
                        for (k=0; k<l1.lLength; k++) {
                            if (k==l2.lLength) {
                                break;
                            }
                            _Variable *thisVar = LocateVar(l1.lData[k]);
                            _Constant nv (thisVar->Compute()->Value()+
                                          LocateVar(l2.lData[k])->Compute()->Value());
                            thisVar->SetValue (&nv);
                        }

                    }
                    DeleteVariable (*killedNode->GetName(),true);
                    delete(p1);

                    // begin paste
                    travNode = me->DepthWiseTraversal (true);

                    node<nodeCoord>*node1 = (node<nodeCoord>*)currentSelection(currentSelection.lLength-1);

                    while (travNode) {
                        if (travNode->GetAVariable()==node1->in_object.varRef) {
                            ct2 = &me->GetCurrentNode();
                            break;
                        }
                        travNode = me->DepthWiseTraversal();
                    }

                    if (!ct2) {
                        WarnError (errMsg);
                        return;
                    }
                    p1 = ct2->parent;
                    for (i=0; i<p1->nodes.length; i++) {
                        if (p1->nodes.data[i]==ct2) {
                            break;
                        }
                    }
                    node <long>* newp = new node<long>;
                    checkPointer ((Ptr)newp);
                    p1->nodes.data[i] = newp;
                    newp->parent = p1;
                    newp->add_node(*ct);
                    newp->add_node(*ct2);

                    k = travNode->GetModelIndex();
                    _String nodeName (empty);
                    if (k != HY_NO_MODEL) {
                        nodeName = *(_String*)modelNames(k);
                        _String dummy ("1.0"), dummy2;
                        me->FinalizeNode(newp,FindUnusedSuffix ("Node"),dummy2,nodeName,dummy);
                        _CalcNode* targetNode = (_CalcNode*)LocateVar (newp->in_object);
                        targetNode->SetCodeBase (me->GetCodeBase());
                        targetNode->CopyMatrixParameters(travNode);
                        _SimpleList varList;
                        {
                            _AVLList vla (&varList);
                            targetNode->ScanForVariables(vla,vla);
                            travNode->ScanForVariables(vla,vla);
                        }
                        //targetNode->ScanForVariables(varList,varList);
                        //travNode->ScanForVariables(varList,varList);
                        for (k=0; k<varList.lLength; k++) {
                            _Constant myValue(LocateVar(varList.lData[k])->Compute()->Value()/2.0);
                            LocateVar(varList.lData[k])->SetValue (&myValue);
                        }
                        dubiousNodes<<targetNode->GetAVariable();
                    } else {
                        _String dummy ("1.0"), dummy2;
                        me->FinalizeNode(newp,FindUnusedSuffix ("Node"),dummy2,nodeName,dummy);
                        _CalcNode* targetNode = (_CalcNode*)LocateVar (newp->in_object);
                        targetNode->CopyMatrixParameters(travNode);
                    }
                }
            } else {
                done = false;
            }
        } else {
            done = false;
        }
    } else {
        done = false;
    }

    if (!done) {
        WarnError (errMsg);
    }

    if (likeFuncID>=0) {
        _LikelihoodFunction* theLF = (_LikelihoodFunction*)likeFuncList(likeFuncID);
        theLF->MapTreeTipsToData (theLF->DependOnTree(*me->GetName()),true);
        theLF->VoidOldResults();
        me->SetUp();
        UpdateLnLikelihood();
    }

    BuildTree(true);

    RenderTree();
    return;
}
//__________________________________________________________

bool    _HYTreePanel::CheckIfNeedUnscaled (void)
{
    _HYPullDown *varMenu = (_HYPullDown*)GetObject (5);
    if (scaleVariable.sLength) {
        if (!promptDefModel) {
            _String msg ("Newly created branches will inherit models/parameter values from the selected node(s)."),
                    msg2("Do not warn again");
            if (!ProceedPromptWithCheck (msg,msg2,promptDefModel, (Ptr)this)) {
                return false;
            }
        }
    } else {
        if (varMenu->MenuItemCount()>3) {
            for (long k=varMenu->MenuItemCount(); k>=3; k--) {
                varMenu->DeleteMenuItem(3);
            }
            varMenu->EnableItem (1,false);
            varMenu->EnableItem (2,false);
        }
    }
    return true;
}

//__________________________________________________________

bool    _HYTreePanel::ShrinkDataFilter (bool doit)
{
    if (likeFuncID>=0) {
        if (!RequestLFDeleteOrAlter (likeFuncID)) {
            return false;
        }

        if (!doit) {
            if (!promptShrinkFilter) {
                _String msg ("I will reduce the number of sequences in the datafilter attached to the tree to compensate for deletion of tips."),
                        msg2("Do not warn again");
                if (!ProceedPromptWithCheck (msg,msg2,promptShrinkFilter)) {
                    return false;
                }
            }
        } else {
            _LikelihoodFunction* theLF = (_LikelihoodFunction*)likeFuncList (likeFuncID);
            if (!theLF->UpdateFilterSize(theLF->DependOnTree(treeName))) {
                _String msg3 ("Could not update data partition. Should I detach this tree from the data?");
                if (ProceedPrompt(msg3,(Ptr)this)) {
                    likeFuncID = -1;
                    return true;
                }
                return false;
            }
        }
    }
    return true;
}

//__________________________________________________________

bool    _HYTreePanel::KillLikeFunc (char code)
{
    if (likeFuncID>=0) {
        if (!RequestLFDeleteOrAlter (likeFuncID)) {
            return false;
        }

        if (!promptKillLF) {
            _String msg ("This operation will destroy the likelihood function currently associated with the tree."),
                    msg2("Do not warn again");
            if (!ProceedPromptWithCheck (msg,msg2,promptKillLF)) {
                return false;
            }
        }
        PrepareUndoData (code);
        postLFKillEvent (GetID(),likeFuncID);
        KillLFRecord (likeFuncID,false);
        likeFuncID = -1;
        //DeleteObject (theLF);
    } else {
        PrepareUndoData (code);
    }
    return true;
}
//__________________________________________________________
_TheTree*  _HYTreePanel::LocateMyTreeVariable()
{

    long f = LocateVarByName (treeName);
    _TheTree* me = nil;
    if (f>=0) {
        _Variable* v = FetchVar(f);
        if (v->ObjectClass()==TREE) {
            me = (_TheTree*)v;
        } else {
            f = -1;
        }
    }
    if (f<0) {
        _String errMsg ("The tree variable ");
        errMsg = errMsg & treeName & " is no longer valid. This tree can still be printed or saved, but can't be rescaled";
        _HYPullDown*  theMenu = (_HYPullDown*)GetObject(3);
        for (f = 0; f<theMenu->MenuItemCount(); f++) {
            theMenu->EnableItem(f,false);
        }
        theMenu = (_HYPullDown*)GetObject(5);
        for (f = 0; f<theMenu->MenuItemCount(); f++) {
            theMenu->EnableItem(f,false);
        }
        WarnError (errMsg);
        return nil;
    }
    return me;
}

//__________________________________________________________
long  _HYTreePanel::FindUnusedSuffix(_String suffix)
{
    _String prefix (treeName), testName;
    prefix = treeName&'.'&suffix;
    long    p = 1;
    while(1) {
        testName = prefix&_String (p);
        if (LocateVarByName(testName)<0) {
            break;
        }
        p++;
    }
    return p;
}

//__________________________________________________________
void _HYTreePanel::SelectAllBranches(void)
{
    node<nodeCoord>* currentNd = NodeTraverser(coordTree), *parent = nil;
    currentSelection.Clear();
    while (currentNd) {
        parent = currentNd->parent;
        if (!parent) {
            break;
        }
        if (currentNd->in_object.branchName.sLength) {
            currentNd->in_object.flags|=HY_BRANCH_SELECT;
            currentSelection<<(long)currentNd;
        }
        currentNd = NodeTraverser((node <nodeCoord>*)nil)   ;
    }
    _UpdateOperationsMenu();
    RenderTree();
}

//__________________________________________________________
void _HYTreePanel::SelectEntireSubtree(void)
{
    node<nodeCoord>* rtNode=(node<nodeCoord>*)currentSelection.lData[0], *currentNd = NodeTraverser(rtNode);
    //currentSelection.Clear();
    while ((currentNd)&&(currentNd!=rtNode)) {
        if (currentNd->in_object.branchName.sLength) {
            currentNd->in_object.flags|=HY_BRANCH_SELECT;
            currentSelection<<(long)currentNd;
        }
        currentNd = NodeTraverser((node <nodeCoord>*)nil)   ;
    }
    _UpdateOperationsMenu();
    RenderTree();
}

//__________________________________________________________
void _HYTreePanel::InvertSelection()
{
    node<nodeCoord> *currentNd = NodeTraverser(coordTree);

    currentSelection.Clear();

    while (currentNd) {
        if (currentNd->in_object.branchName.sLength) {
            if (currentNd->in_object.flags&HY_BRANCH_SELECT) {
                currentNd->in_object.flags-=HY_BRANCH_SELECT;
            } else {
                currentNd->in_object.flags|=HY_BRANCH_SELECT;
                currentSelection<<(long)currentNd;
            }
        }
        currentNd = NodeTraverser((node <nodeCoord>*)nil)   ;
    }
    _UpdateOperationsMenu();
    RenderTree();
}

//__________________________________________________________
void _HYTreePanel::GrowSelection(node<nodeCoord>* currentRoot, bool force)
{
    node<nodeCoord>* startAt = currentRoot;
    if (currentRoot == nil) {
        currentSelection.Clear();
        startAt = coordTree;
    }

    if (startAt->in_object.flags&HY_BRANCH_SELECT) {
        force = true;
    }

    if (force) {
        startAt->in_object.flags|=HY_BRANCH_SELECT;
        currentSelection << (long)startAt;
    }

    for (long kk=startAt->get_num_nodes(); kk>=1; kk--) {
        GrowSelection (startAt->go_down (kk), force);
    }

    if (currentRoot == nil) {
        _UpdateOperationsMenu();
        RenderTree();
    }
}

//__________________________________________________________
void _HYTreePanel::SelectIncompleteBranches(void)
{
    node<nodeCoord>* currentNd = NodeTraverser(coordTree), *parent = nil;
    currentSelection.Clear();
    dubiousNodes.Sort();
    while (currentNd) {
        parent = currentNd->parent;
        if (!parent) {
            break;
        }
        currentNd->in_object.flags&=HY_BRANCH_DESELECT;
        if (dubiousNodes.BinaryFind(currentNd->in_object.varRef)>=0) {
            currentNd->in_object.flags|=HY_BRANCH_SELECT;
            currentSelection<<(long)currentNd;
        }
        currentNd = NodeTraverser((node <nodeCoord>*)nil)   ;
    }
    _UpdateOperationsMenu();
    if (currentSelection.lLength==0) {
        _String msg ("I couldn't find any branches marked as incomplete (That's good!).");
        ProblemReport (msg,(Ptr)this);
    }
    RenderTree();
}

//__________________________________________________________
void _HYTreePanel::SelectBranchesWithoutModel(void)
{
    node<nodeCoord>* currentNd = NodeTraverser(coordTree), *parent = nil;
    currentSelection.Clear();
    while (currentNd) {
        parent = currentNd->parent;
        if (!parent) {
            break;
        }
        if (currentNd->in_object.varRef>=0) {
            currentNd->in_object.flags&=HY_BRANCH_DESELECT;
            if (((_VariableContainer*)LocateVar(currentNd->in_object.varRef))->GetModelIndex()== HY_NO_MODEL) {
                currentNd->in_object.flags|=HY_BRANCH_SELECT;
                currentSelection<<(long)currentNd;
            }
        }
        currentNd = NodeTraverser((node <nodeCoord>*)nil)   ;
    }
    _UpdateOperationsMenu();
    if (currentSelection.lLength==0) {
        _String msg ("All  branches have models attached to them.");
        ProblemReport (msg,(Ptr)this);
    }
    RenderTree();
}
//__________________________________________________________
void _HYTreePanel::PreserveSelection (_SimpleList& s)
{
    for (long k=0; k<currentSelection.lLength; k++) {
        s<<((node<nodeCoord>*)currentSelection(k))->in_object.varRef;
    }
    s.Sort();
}

//__________________________________________________________
void _HYTreePanel::RestoreSelection (_SimpleList& s)
{
    if (s.lLength) {
        node<nodeCoord>* thisNode = NodeTraverser (coordTree);
        while (thisNode) {
            if (s.BinaryFind(thisNode->in_object.varRef)>=0) {
                thisNode->in_object.flags |= HY_BRANCH_SELECT;
                currentSelection<<(long)thisNode;
            }
            thisNode = NodeTraverser((node <nodeCoord>*)nil);
        }
    }
}
//__________________________________________________________
void _HYTreePanel::InvokeNodeEditor (void)
{
    node<nodeCoord>* thisNode = (node<nodeCoord>*) currentSelection.lData[0];

    _TheTree* me = LocateMyTreeVariable();
    long      f = -1L;
    if (!me) {
        return;
    }

    _String     nodeName (treeName);
    nodeName =  nodeName&'.'&thisNode->in_object.branchName;

    _CalcNode*  thisCNode = (_CalcNode*)LocateVar (thisNode->in_object.varRef);

    _String     startingName (*thisCNode->GetName());

    long        nodeModel = thisCNode->GetTheModelID();

    bool        doSCV = false,
                incFlag;

    if (currentSelection.lLength==1UL) {
        f = dubiousNodes.Find(thisNode->in_object.varRef);
    }
    incFlag = (f>=0);

    bool    doBuild;
    _HYNodeInfoDialog * ndd = new _HYNodeInfoDialog (nodeName, &nodeModel, currentSelection.lLength==1,LocateMyTreeVariable(),&currentSelection,&doSCV,&incFlag,&doBuild,this);
    ndd->Activate();
    while (windowObjectRefs.Find ((long)ndd)>=0) {
        handleGUI();
    }

    bool    doDraw = false,
            doLFR  = false;
    if (currentSelection.lLength==1) {
        thisCNode = (_CalcNode*)FetchVar(LocateVarByName(startingName));
        if (!thisCNode) {
            nodeName = "Internal Tree Edit Error.";
            FlagError (nodeName);
        }
        if (nodeName.Find('.')<0) {
            startingName.Trim (startingName.Find('.')+1,-1);
        }
        if (!nodeName.Equal (&startingName))
            // node renamed
        {
            doDraw = (treeFlags&HY_TREEPANEL_TIP_LABELS)||(treeFlags&HY_TREEPANEL_INT_LABELS);
            _String*  nn = thisCNode->GetName();
            thisNode->in_object.branchName = nodeName;
            nodeName = nn->Cut(0,nn->Find('.'))&nodeName;
            RenameVariable (thisCNode->GetName(),&nodeName);
            doLFR = true;
            postChangeLFParamsEvent(GetID(),likeFuncID);
        }
        if (thisCNode->GetAVariable()!=thisNode->in_object.varRef) {
            node<long>* cNode = DepthWiseStepTraverser(&me->GetRoot());
            while (cNode) {
                if (cNode->in_object==thisNode->in_object.varRef) {
                    cNode->in_object = thisCNode->GetAVariable();
                    break;
                }
                cNode = DepthWiseStepTraverser ((node<long>*)nil);
            }
            thisNode->in_object.varRef = thisCNode->GetAVariable();
        }
    }

    if ((doSCV)&&(!doLFR)) {
        doLFR = true;
        if (likeFuncID>=0) {
            ((_LikelihoodFunction*)likeFuncList(likeFuncID))->RescanAllVariables();
            postChangeLFParamsEvent(GetID(),likeFuncID);
        }
    }

    if (incFlag) {
        if (currentSelection.lLength==1) {
            if (f<0) {
                dubiousNodes<<thisCNode->GetAVariable();
            }
        }
    } else {
        if ((f>=0)&&(currentSelection.lLength==1)) {
            dubiousNodes.Delete (f);
        }
    }

    if (!doBuild) {
        doBuild = me->HasChanged();
        if (doBuild&&(likeFuncID==-1)) {
            me->MarkDone();
        }

    } else {
        SetVariableReference (treeName);
        if (!doLFR) {
            postChangeLFEvent(GetID(),likeFuncID);
        }
        return;
    }

    if (doBuild) {
        PrepareUndoData (-1);
        _SimpleList     saveSel;
        PreserveSelection (saveSel);
        BuildTree(true);
        RestoreSelection (saveSel);
        if (!doLFR) {
            postChangeLFEvent(GetID(),likeFuncID);
        }
        doDraw = true;
    } else if (doSCV) {
        UpdateScalingVariablesList();
    }
    if (doDraw) {
        RenderTree ();
    }

}

//__________________________________________________________

void    _HYTreePanel::FlushUndoData (void)
{
    if (undoCode>0) {
        undoNodeList.Clear();
        if ((undoCode>2)&&(undoTree)) {
            KillNodesInClipboard (undoTree);
            undoTree = nil;
        }
        if ((undoCode>6)&&(undoCode!=9)) {
            if (undoLFPointer) {
                DeleteObject ((BaseRef)undoLFPointer);
                undoLFPointer = nil;
            }
        }
    }
    undoCode = -1;
}

//__________________________________________________________

void    _HYTreePanel::PrepareUndoData (char code)
// codes:
// 1 - swap
// 2 - flip
// 3 - collapse
// 4 - Delete Tips
// 5 - Cut
// 6 - Graft
// 7 - Paste
{
    FlushUndoData();
    if ((code==1)||(code==2)) {
        undoCode = code;
        undoNodeList.Duplicate (&currentSelection);
    } else if (code>2) {
        undoCode = code;
        makeUndoTree = true;
        _TheTree* me = LocateMyTreeVariable();
        bool      dummy;
        undoTree = CopyNodesToClipboard(&me->GetRoot(),me,dummy);
        makeUndoTree = false;
        if ((code>6)&&(code!=9)) {
            undoLFPointer = nil;
            if (likeFuncID>=0) {
                undoLFName = *(_String*)likeFuncNamesList (likeFuncID);
                undoLFPointer = (Ptr)likeFuncList (likeFuncID);
                ((BaseRef)undoLFPointer)->nInstances++;
            }
        }
    }
}

//__________________________________________________________

void    _HYTreePanel::UndoLastOperation (void)
{
    switch (undoCode) {
    case 1: {
        currentSelection.Clear();
        currentSelection.Duplicate(&undoNodeList);
        SwapSelectedSubTrees();
        break;
    }
    case 2: {
        currentSelection.Clear();
        currentSelection.Duplicate(&undoNodeList);
        FlipSelectedBranches();
        break;
    }
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9: {
        _TheTree* me = LocateMyTreeVariable();
        KillNodesInClipboard (&me->GetRoot());
        me->SetRoot (undoTree);
        RenameUndoNodes (undoTree,*me->GetName());
        undoTree = nil;
        if ((undoCode>6)&&(undoCode!=9)) {
            if (undoLFPointer) {
                likeFuncID = likeFuncList.lLength;
                likeFuncNamesList && & undoLFName;
                likeFuncList << (_LikelihoodFunction*) undoLFPointer;
                undoLFPointer = nil;
            }
        } else if ((undoCode == 9)&&(likeFuncID>=0)) {
            _LikelihoodFunction* theLF = (_LikelihoodFunction*)likeFuncList(likeFuncID);
            theLF->MapTreeTipsToData (theLF->DependOnTree(*me->GetName()),true);
        }
        UpdateLnLikelihood();
        if ((undoCode%2==0)&&(likeFuncID>=0)) {
            ShrinkDataFilter (true);
        }

        BuildTree(true);
        RenderTree();
    }
    }
    FlushUndoData();
    _UpdateOperationsMenu();
}

//__________________________________________________________

void  _HYTreePanel::FishEyeProjection (long fx, long fy, long hSize, long vSize, node<nodeCoord>* thisNode)
{
    if (treeFlags&HY_TREEPANEL_PROJECTION) {
        thisNode->in_object.h = thisNode->in_object.auxL;
        thisNode->in_object.v = thisNode->in_object.auxD;
    } else {
        thisNode->in_object.auxL = thisNode->in_object.h;
        thisNode->in_object.auxD = thisNode->in_object.v;
    }

    bool       side = thisNode->in_object.h<fx;

    _Parameter dMax = side?fx:hSize-fx,
               gArg = fabs(thisNode->in_object.h-fx)/dMax;

    gArg = ((distortion+1.)*gArg/(distortion*gArg+1.));

    if (side) {
        thisNode->in_object.h = fx-gArg*dMax;
    } else {
        thisNode->in_object.h = fx+gArg*dMax;
    }

    thisNode->in_object.v =  vSize-thisNode->in_object.v;
    side = thisNode->in_object.v<fy;
    dMax = side?fy:vSize-fy;
    gArg = fabs(thisNode->in_object.v-fy)/dMax;
    gArg     = (distortion+1.)*gArg/(distortion*gArg+1.);

    if (side) {
        thisNode->in_object.v = vSize-(fy-gArg*dMax);
    } else {
        thisNode->in_object.v = vSize-(gArg*dMax+fy);
    }

    for (long k = thisNode->get_num_nodes(); k; k--) {
        FishEyeProjection (fx,fy,hSize, vSize,thisNode->go_down(k));
    }

    /*if ((treeFlags&HY_TREEPANEL_STRAIGHT)&&(thisNode->parent == nil))
    {
        PreStraightenEdges (coordTree);
        StraightenEdges(coordTree,0,HY_TREEPANEL_MARGIN,-1);
    }*/
}

//__________________________________________________________

void  _HYTreePanel::HandleFisheyeButton (void)
{
    _HYButtonBar * bb = (_HYButtonBar*)GetObject (7);
    int          h,v;

    bb->GetButtonLoc(5,h,v,true);

    _String menuChoice ("Set Distortion Parameter"),
            prStr      ("Select Distortion Parameter (>=0):");

    _List   menuAll;

    menuAll && & menuChoice;

    if (HY_TREEPANEL_PROJECTION&treeFlags) {
        menuChoice = "Restore Normal View";
        menuAll && & menuChoice;
    }

    menuChoice = HandlePullDown (menuAll,h,v,0);
    bb->_UnpushButton();

    if (menuChoice.sLength) {
        h = menuAll.FindObject (&menuChoice);
        if (h==0) {
            menuChoice = distortion;
            if (EnterStringDialog (menuChoice, prStr, (Ptr)this)) {
                distortion = menuChoice.toNum();
            }
        } else {
            RenderTree ();
        }
    }
}

//__________________________________________________________

void  _HYTreePanel::HandleViewOptions (void)
{
    _HYTreeViewDialog * fD = new _HYTreeViewDialog (this,treeFlags&HY_TREEPANEL_CIRCULAR);
    fD->Activate();
}

//__________________________________________________________

void  _HYTreePanel::HandleLabels (bool below)
{
    _List options;

    options && & none;

    _HYPullDown * p1 = (_HYPullDown*) GetObject (5);

    _TheTree *t = LocateMyTreeVariable ();
    if (t) {
        long            mI = p1->MenuItemCount();

        if (mI>3) {
            options && & eSubsScale;
        }

        if (t->HaveStringBranchLengths()) {
            options && & assValScale;
        }

        for (long k=3; k<mI; k++) {
            options << p1->GetMenuItem (k);
        }

        _HYLabelDialog * ld = new _HYLabelDialog (this, options, below);
        ld->Activate();
    }
}


//__________________________________________________________

void  _HYTreePanel::ToggleScaleOption (void)
{
    bool onOff = treeFlags & HY_TREEPANEL_SCALE_TO_WINDOW;
    _HYButtonBar* zoom = (_HYButtonBar*)GetObject (7);
    for (long k = 0; k<6; k++) {
        zoom->EnableButton (k,onOff);
    }

    _HYCanvas* theCanvas = (_HYCanvas*)GetObject (0),
               * theRuler  = (_HYCanvas*)GetObject (2);

    _HYCanvas* newCanvas,
               * newRuler;

    _HYRect         canvasSettings =  {50,50,theCanvas->rel.bottom-theCanvas->rel.top,
                                       theCanvas->rel.right-theCanvas->rel.left,
                                       HY_COMPONENT_NO_SCROLL
                                      },

                    rulerSettings = canvasSettings;

    rulerSettings.top = rulerSettings.bottom = HY_TREEPANEL_RULER_EXPANDED;
    if (onOff) {
        newCanvas = new _HYCanvas (canvasSettings, GetOSWindowData(), canvasSettings.bottom, canvasSettings.right, theCanvas->depth);
        newRuler  = new _HYCanvas (rulerSettings, GetOSWindowData(), rulerSettings.bottom, rulerSettings.right, theRuler->depth);
        treeFlags -= HY_TREEPANEL_SCALE_TO_WINDOW;
        zoom->settings.top = zoom->settings.bottom = 50;
        zoom->SetButtonLayoutW(4);
        SetTableDimensions (6,3);
        SetCell (0,0,GetObject (1));
        SetCell (0,1,GetObject (4));
        SetCell (0,2,GetObject (3));
        SetCell (1,0,GetObject (1));
        SetCell (1,1,GetObject (6));
        SetCell (1,2,GetObject (5));
        SetCell (2,0,GetObject (1));
        SetCell (2,1,GetObject (8));
        SetCell (2,2,GetObject (7));
        SetCell (3,0,GetObject (1));
        SetCell (3,1,GetObject (9));
        SetCell (3,2,GetObject (9));
        SetCell (4,0,GetObject (2));
        SetCell (4,1,GetObject (2));
        SetCell (4,2,GetObject (2));
        SetCell (5,0,GetObject (0));
        SetCell (5,1,GetObject (0));
        SetCell (5,2,GetObject (0));
    } else {
        newCanvas = new _HYStretchCanvas (canvasSettings, GetOSWindowData(),canvasSettings.bottom, canvasSettings.right, theCanvas->depth, HY_SCANVAS_HORIZONTAL|HY_SCANVAS_VERTICAL);
        newRuler  = new _HYStretchCanvas (rulerSettings, GetOSWindowData(), rulerSettings.bottom, rulerSettings.right , theRuler->depth,HY_SCANVAS_HORIZONTAL);
        treeFlags += HY_TREEPANEL_SCALE_TO_WINDOW;
        zoom->settings.top = zoom->settings.bottom = 30;
        zoom->SetButtonLayoutW(8);
        SetTableDimensions (3,3);
        SetCell (0,0,GetObject (3));
        SetCell (0,1,GetObject (5));
        SetCell (0,2,GetObject (7));
        SetCell (1,0,GetObject (2));
        SetCell (1,1,GetObject (2));
        SetCell (1,2,GetObject (2));
        SetCell (2,0,GetObject (0));
        SetCell (2,1,GetObject (0));
        SetCell (2,2,GetObject (0));
    }
    newCanvas->SetDimensions (canvasSettings,theCanvas->rel);
    newRuler->SetDimensions (rulerSettings,theRuler->rel);

    checkPointer    (newCanvas);
    checkPointer    (newRuler);

    if (!onOff) {
        ((_HYStretchCanvas*)newRuler)->SetMessageRecipient (this);
    }

    components.Replace (0,newCanvas,false);
    components.Replace (2,newRuler,false);

    FitToWindow(true);
    dim = MinMaxWindowDimensions ();
    UpdateComponentInfo ();
    /*RenderTree ();
    RenderRuler();*/
}

//__________________________________________________________
// _HYLabelDialog
//__________________________________________________________

_HYLabelDialog::_HYLabelDialog (_HYTreePanel* rec, _List& options,bool below):_HYFontDialog (below?rec->branchLabel2:rec->branchLabel1,rec)
{

    //long          index,
    //              cellWidth;

    _HYColor        bgc     = GetDialogBackgroundColor();

    _HYRect         canvasSettings = {30,100,30,100,HY_COMPONENT_NO_SCROLL|HY_COMPONENT_TRANSP_BG};

    _HYLabel*       l1      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l2      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l3      = new _HYLabel (canvasSettings, GetOSWindowData());


    canvasSettings.left = canvasSettings.right = 200;

    _HYPullDown*    p1      = new _HYPullDown (canvasSettings, GetOSWindowData());
    _HYPullDown*    p2      = new _HYPullDown (canvasSettings, GetOSWindowData());
    _HYPullDown*    p3      = new _HYPullDown (canvasSettings, GetOSWindowData());

    p1->SetMessageRecipient (this);
    p2->SetMessageRecipient (this);
    p3->SetMessageRecipient (this);

    AddObject (l1);
    AddObject (l2);
    AddObject (l3);

    AddObject (p1);
    AddObject (p2);
    AddObject (p3);

    _SimpleList saveCells (cells);

    long k;

    SetTableDimensions (8,3);

    SetCell   (0,0,l1);
    SetCell   (0,1,p1);
    SetCell   (0,2,p1);

    SetCell   (1,0,l2);
    SetCell   (1,1,p2);
    SetCell   (1,2,p2);

    SetCell   (2,0,l3);
    SetCell   (2,1,p3);
    SetCell   (2,2,p3);

    for (k=0; k<15; k++) {
        cells.lData[k+9] = saveCells.lData[k];
    }

    _HYFont  labelFont;
    l1->SetBackColor (bgc);
    l2->SetBackColor (bgc);
    l3->SetBackColor (bgc);

    p1->SetBackColor (bgc);
    p2->SetBackColor (bgc);
    p3->SetBackColor (bgc);

#ifdef __WINDOZE__
    labelFont.face  = "Arial";
    labelFont.size  = 14;
#else
    labelFont.face  = "System Font";
    labelFont.size  = 12;
#endif
    labelFont.style = HY_FONT_PLAIN;

    l1->SetFont (labelFont);
    l2->SetFont (labelFont);
    l3->SetFont (labelFont);

    l1->SetAlignFlags (HY_ALIGN_LEFT);
    l2->SetAlignFlags (HY_ALIGN_LEFT);
    l3->SetAlignFlags (HY_ALIGN_LEFT);

    p1->SetAlignFlags (HY_ALIGN_LEFT);
    p2->SetAlignFlags (HY_ALIGN_LEFT);
    p3->SetAlignFlags (HY_ALIGN_LEFT);

    l1->SetText       ("Label with:");
    l2->SetText       ("Aligh at:");
    l3->SetText       ("Max. Digits:");

    p2->AddMenuItem   ("Left",-1);
    p2->AddMenuItem   ("Center",-1);
    p2->AddMenuItem   ("Right",-1);

    char flag = below?rec->bottomAlign:rec->topAlign;
    if (flag&HY_ALIGN_LEFT) {
        p2->ChangeSelection (0,false);
    } else if (flag&HY_ALIGN_RIGHT) {
        p2->ChangeSelection (2,false);
    } else {
        p2->ChangeSelection (1,false);
    }


    for (k=0; k<options.lLength; k++) {
        p1->AddMenuItem (*(_String*)options(k),-1);
    }

    _String * scv = &(below?rec->branchVar2:rec->branchVar1);

    k = options.FindObject (scv);

    if (k==-1) {
        long f;
        k=0;
        if (scv->Equal (&expectedNumberOfSubs) && ((f=options.FindObject (&eSubsScale))>=0)) {
            k = f;
        } else if (scv->Equal (&stringSuppliedLengths) && ((f=options.FindObject (&assValScale))>=0)) {
            k = f;
        }
    }

    p1->ChangeSelection (k,false);

    for (k=1; k<16; k++) {
        p3->AddMenuItem (_String (k),-1);
    }

    p3->ChangeSelection ((below?rec->labelDigits2:rec->labelDigits1)-1,false);

    DeleteObject (l1);
    DeleteObject (l2);
    DeleteObject (l3);
    DeleteObject (p1);
    DeleteObject (p2);
    DeleteObject (p3);

    if (below) {
        SetTitle ("Below Branch Labels");
    } else {
        SetTitle ("Above Branch Labels");
    }

    which = below;

    _HYRect     dim = MinMaxWindowDimensions();

    SetWindowRectangle (0,0,dim.bottom,dim.right);
    CenterWindow       (this);
}

//__________________________________________________________

bool    _HYLabelDialog::ProcessEvent (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    long        i,f;

    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);

        if (i==7) { // OK
            _HYTreePanel* tp = (_HYTreePanel*) mr;
            bool          needRedraw = false,
                          labelsOn = which?tp->treeFlags&HY_TREEPANEL_LABEL2:tp->treeFlags&HY_TREEPANEL_LABEL1,
                          needResize = false;

            _HYFont      *oldFont = &(which?tp->branchLabel2:tp->branchLabel1);

            if ((!myFont.face.Equal(&oldFont->face))||(myFont.size!=oldFont->size)||(myFont.style!=oldFont->style)) {
                *oldFont = myFont;
                if (labelsOn)
                    if (myFont.size-3!=tp->windowTextMarginV)
                        if (myFont.size>=(which?tp->branchLabel1.size:tp->branchLabel2.size))
                            if (myFont.size+3>=tp->treeLabelFont.size/2+5) {
                                tp->windowTextMarginV = myFont.size+3;
                                needResize = true;
                            }

                needRedraw = true;
            }

            _HYPullDown * p1 = (_HYPullDown*)components (12),
                          * p2 = (_HYPullDown*)components (13),
                            * p3 = (_HYPullDown*)components (14);

            _String      newScaleVar,
                         *oldScaleVar;

            oldScaleVar = (_String*)p1->GetMenuItem (p1->GetSelection());
            if (oldScaleVar->Equal(&eSubsScale)) {
                newScaleVar = expectedNumberOfSubs;
            } else if (oldScaleVar->Equal(&assValScale)) {
                newScaleVar = stringSuppliedLengths;
            } else if (!oldScaleVar->Equal(&none)) {
                newScaleVar = *oldScaleVar;
            }

            oldScaleVar = &(which?tp->branchVar2:tp->branchVar1);

            if (!newScaleVar.Equal(oldScaleVar)) {
                needRedraw = true;
                if (newScaleVar.sLength==0) {
                    tp->treeFlags -= (which ? HY_TREEPANEL_LABEL2 : HY_TREEPANEL_LABEL1);
                } else {
                    if (oldScaleVar->sLength == 0) {
                        tp->treeFlags += (which ? HY_TREEPANEL_LABEL2 : HY_TREEPANEL_LABEL1);
                        labelsOn = true;
                        if (myFont.size-3!=tp->windowTextMarginV)
                            if (myFont.size>=(which?tp->branchLabel1.size:tp->branchLabel2.size))
                                if (myFont.size+3>=tp->treeLabelFont.size/2+5) {
                                    tp->windowTextMarginV = myFont.size+3;
                                    needResize = true;
                                }
                    }
                    _TheTree * mt = tp->LocateMyTreeVariable();
                    if (mt) {
                        mt->AssignLabelsToBranches (tp->coordTree,&newScaleVar,which);
                    }
                }
                *oldScaleVar = newScaleVar;
            }

            i = p3->GetSelection ();
            char   * oldDigits = &(which?tp->labelDigits2:tp->labelDigits1);

            if (i+1!=*oldDigits) {
                needRedraw = true;
                *oldDigits = i+1;
            }

            i = p2->GetSelection ();
            oldDigits = &(which?tp->bottomAlign:tp->topAlign);

            if (i==0) {
                i = HY_ALIGN_LEFT;
            } else if (i==1) {
                i = 0;
            } else {
                i = HY_ALIGN_RIGHT;
            }
            if (i!=*oldDigits) {
                needRedraw = true;
                *oldDigits = i;
            }

            if (needRedraw&&labelsOn&&(!tp->IsVertical())&&(!((tp->treeFlags&HY_TREEPANEL_STRAIGHT)||(tp->treeFlags&HY_TREEPANEL_ARCS)||(tp->treeFlags&HY_TREEPANEL_CIRCULAR)))) {
                if (needResize) {
                    tp->SetVSpace (tp->vSpace,true);
                } else {
                    tp->RenderTree();
                }
            }
            postWindowCloseEvent (GetID());
            done = true;
        }
    }

    if (done) {
        DeleteObject (e);
        return true;
    }
    return _HYFontDialog::ProcessEvent(e);
}


//__________________________________________________________
// _HYTreeViewDialog
//__________________________________________________________

_HYTreeViewDialog::_HYTreeViewDialog (_HYTreePanel* rec, bool radial):_HYFontDialog (rec->treeLabelFont,rec)
{

    //long          index,
    //              cellWidth;

    _HYColor        bgc     = GetDialogBackgroundColor();
    _HYRect         canvasSettings = {30,300,30,300,HY_COMPONENT_NO_SCROLL|HY_COMPONENT_TRANSP_BG};

    _HYCheckbox*    c1      = new _HYCheckbox (canvasSettings, GetOSWindowData());

    canvasSettings.left = canvasSettings.right = 200;
    _HYTextBox*     t1      = new _HYTextBox (canvasSettings, GetOSWindowData());
    canvasSettings.left = canvasSettings.right = 120;
    _HYTextBox*     t2      = new _HYTextBox (canvasSettings, GetOSWindowData());
    _HYTextBox*     t3      = new _HYTextBox (canvasSettings, GetOSWindowData());

    canvasSettings.left = canvasSettings.right = 100;
    _HYLabel*       l1      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l2      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l3      = new _HYLabel (canvasSettings, GetOSWindowData());

    canvasSettings.left = canvasSettings.right = 81;
    canvasSettings.top = canvasSettings.bottom = 59;

    _HYCanvas*      cc      = new _HYCanvas(canvasSettings, GetOSWindowData(), 59, 81, 32);

    cc->StartDraw();
    cc->SetDialogBG();
    cc->EndDraw();

    c1->SetMessageRecipient (this);
    t1->SetMessageRecipient (this);
    t2->SetMessageRecipient (this);
    t3->SetMessageRecipient (this);

    AddObject (c1); // 9
    AddObject (l1); // 10
    AddObject (t1); // 11
    AddObject (t2); // 12
    AddObject (t3); // 13
    AddObject (cc); // 14
    AddObject (l2); // 15
    AddObject (l3); // 16

    AddKeyboardChainObject (t1);
    AddKeyboardChainObject (t2);
    AddKeyboardChainObject (t3);

    _SimpleList saveCells (cells);

    long k;

    SetTableDimensions (9,3);

    for (k=0; k<12; k++) {
        cells.lData[k] = saveCells.lData[k];
    }

    SetCell (4,0,l1);
    SetCell (4,1,t1);
    SetCell (4,2,t1);

    SetCell (5,0,l2);
    SetCell (5,1,t2);
    SetCell (5,2,cc);

    SetCell (6,0,l3);
    SetCell (6,1,t3);
    SetCell (6,2,cc);

    SetCell (7,0,c1);
    SetCell (7,1,c1);
    SetCell (7,2,c1);

    for (k=24; k<27; k++) {
        cells.lData[k] = saveCells.lData[k-12];
    }


    _HYFont  labelFont;
    c1->SetBackColor (bgc);
    t1->SetBackColor (bgc);
    t2->SetBackColor (bgc);
    t3->SetBackColor (bgc);
    l1->SetBackColor (bgc);
    l2->SetBackColor (bgc);
    l3->SetBackColor (bgc);

#ifdef __WINDOZE__
    labelFont.face  = "Arial";
    labelFont.size  = 14;
#else
    labelFont.face  = "System Font";
    labelFont.size  = 12;
#endif
    labelFont.style = HY_FONT_PLAIN;

    c1->SetFont (labelFont);
    c1->SetAlignFlags (HY_ALIGN_LEFT);
    c1->SetState (rec->treeFlags & HY_TREEPANEL_SCALE_TO_WINDOW);

    t1->SetFont (labelFont);
    t1->SetAlignFlags (HY_ALIGN_LEFT);
    lw    = ((_HYTreePanel*)rec)->branchWidth;
    t1->SetText (_String ((long)lw));


    t2->SetFont (labelFont);
    t2->SetAlignFlags (HY_ALIGN_LEFT);
    sAng      = ((_HYTreePanel*)rec)->arcStart;
    t2->SetText (_String (sAng));

    t3->SetFont (labelFont);
    t3->SetAlignFlags (HY_ALIGN_LEFT);
    fAng      = ((_HYTreePanel*)rec)->arcEnd;
    t3->SetText (_String (fAng));

    l1->SetFont (labelFont);
    l1->SetAlignFlags (HY_ALIGN_LEFT);
    l1->SetText ("Line width:");

    l2->SetFont (labelFont);
    l2->SetAlignFlags (HY_ALIGN_LEFT);
    l2->SetText ("Start Angle:");

    l3->SetFont (labelFont);
    l3->SetAlignFlags (HY_ALIGN_LEFT);
    l3->SetText ("End Angle:");

    c1->SetText       ("Scale tree by resizing the window");

    DeleteObject (c1);
    DeleteObject (t1);
    DeleteObject (l1);
    DeleteObject (t2);
    DeleteObject (l2);
    DeleteObject (t3);
    DeleteObject (l3);
    DeleteObject (cc);


    tb1 = false;
    tb2 = false;
    tb3 = false;

    DrawArcPreview ();

    SetTitle ("Tree Display Options");

    which = radial;

    _HYRect     dim = MinMaxWindowDimensions();

    SetWindowRectangle (0,0,dim.bottom,dim.right);
    CenterWindow       (this);
}

//__________________________________________________________

bool    _HYTreeViewDialog::ProcessEvent (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    long        i,f;

    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);
        if (i==7) { // OK
            done = false;
            bool   rd = false;
            if (lw != ((_HYTreePanel*)mr)->branchWidth) {
                ((_HYTreePanel*)mr)->branchWidth = lw;
                done = true;
            }

            if ((sAng!=((_HYTreePanel*)mr)->arcStart)||(fAng!=((_HYTreePanel*)mr)->arcEnd)) {
                ((_HYTreePanel*)mr)->arcStart = sAng;
                ((_HYTreePanel*)mr)->arcEnd   = fAng;
                done = which;
                if (done) {
                    ((_HYTreePanel*)mr)->BuildTree (true);
                }
            }

            if ((!myFont.face.Equal(&firstFont.face))||(myFont.size!=firstFont.size)||(myFont.style!=firstFont.style)) {
                mr->SetFont (myFont);
                rd = true;
            }

            if (((_HYCheckbox*)GetCellObject (7,0))->GetState() != ((bool)((((_HYTreePanel*)mr)->treeFlags & HY_TREEPANEL_SCALE_TO_WINDOW)))) {
                ((_HYTreePanel*)mr)->ToggleScaleOption ();
                rd = true;
            }

            if (done && (!rd)) {
                ((_HYTreePanel*)mr)->RenderTree();
            }

            postWindowCloseEvent (GetID());
            done = true;
        }
    } else if (e->EventClass()==_hyTextEditChange) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);
        _HYTextBox* tb = nil,
                    * tbb;
        _HYButton*  kb = (_HYButton*)components (7);
        _Parameter  t;
        bool        disable = false;
        tb = (_HYTextBox*)components (i);
        if (i==11) {
            f = tb->GetText().toNum();
            if ((f<=0)||(f>=127)) {
                disable = true;
                tb1 = true;
            } else {
                if (tb1) {
                    tb1 = false;
                    tb->SetForeColor ((_HYColor) {
                        0,0,0
                    });
                }
                lw = f;
            }
        }
        if ((i==12)||(i==13)) {
            t = tb->GetText().toNum();
            if (i==12) {
                if ((t<0.0)||(t>6.284)||(t>=fAng)) {
                    disable = true;
                    tb2 = true;
                } else {
                    if (tb2) {
                        tb2 = false;
                        tb->SetForeColor ((_HYColor) {
                            0,0,0
                        });
                    }
                    sAng = t;
                    if (tb3) {
                        tbb = (_HYTextBox*)components (13);
                        t = tbb->GetText().toNum();
                        if ((t>sAng)&&(t<=6.284)) {
                            tb3 = false;
                            tbb->SetForeColor ((_HYColor) {
                                0,0,0
                            });
                            fAng = t;
                        }
                    }
                }
            } else {
                if ((t<0.0)||(t>6.284)||(t<=sAng)) {
                    disable = true;
                    tb3 = true;
                } else {
                    if (tb3) {
                        tb3 = false;
                        tb->SetForeColor ((_HYColor) {
                            0,0,0
                        });
                    }
                    fAng = t;
                    if (tb2) {
                        tbb = (_HYTextBox*)components (12);
                        t = tbb->GetText().toNum();
                        if ((t<fAng)&&(t>=0)) {
                            tb2 = false;
                            tbb->SetForeColor ((_HYColor) {
                                0,0,0
                            });
                            sAng = t;
                        }
                    }
                }

            }
            DrawArcPreview();
        }


        if (disable ) {
            kb->EnableButton (false);
            tb->SetForeColor ((_HYColor) {
                255,0,0
            });
        } else {
            if (!(tb1||tb2||tb3)) {
                kb->EnableButton(true);
            }
        }
    }

    if (done) {
        DeleteObject (e);
        return true;
    }
    return _HYFontDialog::ProcessEvent(e);
}

//__________________________________________________________

void    _HYTreeViewDialog::DrawArcPreview (void)
{
    _HYCanvas * cc = (_HYCanvas*) components (14);
    cc->StartDraw();
    //cc->SetColor (GetDialogBackgroundColor());
    _HYRect      drawer = {0,0,cc->h, cc->w,1};
    //cc->FillRect (drawer);
    cc->EraseAll();
    if (!(tb2||tb3)) {
        drawer.left  = 15;
        drawer.right = 65;
        drawer.top   = 5;
        drawer.bottom = 55;

        cc->SetColor ((_HYColor) {
            120,120,120
        });
        cc->FillArc  (drawer,sAng*180./3.1415926-90,(fAng-sAng)*180./3.1415926);
        drawer.width = 2;
        cc->SetColor ((_HYColor) {
            0,0,0
        });
        cc->DrawArc  (drawer,sAng*180./3.1415926-90,(fAng-sAng)*180./3.1415926);
    }
    cc->EndDraw();
    cc->_MarkForUpdate();
}

//__________________________________________________________

void    _HYTreePanel::HandleTreeSave (long c, _String filePath)
{
    FILE*       theFile = doFileOpen (filePath.sData,"w");
    if (theFile) {
        switch (c) {
        case 1: {
            _String res = GetTreeString();
            fprintf (theFile,"#NEXUS\nBEGIN TREES;\n\tTree %s = %s;\nEND;",GetTitle().Cut(5,-1).getStr(),res.getStr());
        }
        break;
        case 2:
        case 3:
        case 4:
            if (_TheTree    * me    = LocateMyTreeVariable()) {
                if (c==2) {
                    _List       nodeNames;
                    _AVLListXL  nodeRemap (&nodeNames);
                    _CalcNode*  cNode  = me->LeafWiseTraversal(true);
                    long        nodeCounter = 1;
                    fprintf (theFile,"#NEXUS\nBEGIN TREES;\n\tTRANSLATE");

                    while (cNode) {
                        if (nodeCounter > 1) {
                            fprintf (theFile,",");
                        }
                        _String * rmp = new _String (nodeCounter++),
                        *nname = new _String;
                        checkPointer (nname);
                        checkPointer (rmp);
                        me->GetNodeName (&me->GetCurrentNode(),*nname,false);
                        nodeRemap.Insert (nname,(long)rmp);
                        fprintf (theFile,"\n\t\t%s %s",rmp->sData,nname->sData);
                        DeleteObject (rmp);
                        cNode  = me->LeafWiseTraversal(false);
                    }
                    _String res = GetTreeString(&nodeRemap);
                    fprintf (theFile,";\n\tTree %s = %s;\nEND;",GetTitle().Cut(5,-1).getStr(),res.getStr());
                } else {
                    _Matrix   plot_dim (1,2,false,true);
                    _HYCanvas* theCanvas = (_HYCanvas*)GetObject (0);
                    plot_dim.theData[0] = theCanvas->w;
                    plot_dim.theData[1] = theCanvas->h;
                    _FString  plot_sc  (new _String(((treeFlags&HY_TREEPANEL_SCALED)&&scaleVariable.sLength)?scaleVariable:empty));

                    _AssociativeList*   viewOptions = new _AssociativeList ();
                    viewOptions->MStore (treeOutputLayout,new _Constant(c==4), false);
                    setParameter (treeOutputAVL,viewOptions,false);

                    _FString   *res = (_FString*)me->PlainTreeString (&plot_sc,&plot_dim);
                    fwrite (res->theString->sData,1,res->theString->sLength,theFile);
                    DeleteObject (res);
                }
            }
            break;
        default: {
            _String res = GetTreeString();
            fwrite (res.sData,1,res.sLength,theFile);
        }
        }

        fclose (theFile);
    } else {
        filePath = _String ("Couldn't open file '") & filePath & "' for writing.";
        ProblemReport (filePath,(Ptr)this);
    }
}

//__________________________________________________________

void    _HYTreePanel::GenerateDistanceTable (char opt)
{
    _TheTree* me = LocateMyTreeVariable();
    if (me) {
        long         treeNameLength = me->GetName()->sLength+1;

        _String      reportName;
        _List        vNames;
        _Matrix*     dataMatrix;

        if (opt==1)
            // branch length distribution
        {
            _List        iBranchNames,
                         leafNames,
                         iBranchLengths,
                         leafLengths;

            _SimpleList  sortingOrder,
                         sortingOrder2;

            if (currentSelection.lLength == 0) {
                node<nodeCoord>* meNode = NodeTraverser(coordTree);
                while (meNode->parent) {
                    if (meNode->in_object.varRef>=0) {
                          if (meNode->get_num_nodes()) {
                              sortingOrder    << (long)meNode;
                              iBranchLengths.AppendNewInstance(new _Constant(meNode->in_object.bL));
                          } else {
                              sortingOrder2   << (long)meNode;
                              leafLengths.AppendNewInstance(new _Constant (meNode->in_object.bL));
                          }
                    }

                    meNode = NodeTraverser((node<nodeCoord>*)NULL);
                }
            } else {
                for (long k=0; k<currentSelection.lLength; k++) {
                    node <nodeCoord> * meNode = (node <nodeCoord>*)currentSelection(k);
                    if (meNode->get_num_nodes()) {
                        sortingOrder    << (long)meNode;
                        iBranchLengths.AppendNewInstance(new  _Constant (meNode->in_object.bL));
                    } else {
                        sortingOrder2   << (long)meNode;
                        leafLengths.AppendNewInstance(new _Constant (meNode->in_object.bL));
                    }
                }
            }

            for (long k=0; k<sortingOrder.lLength; k++) {
                node <nodeCoord> * meNode = (node <nodeCoord>*)sortingOrder(k);
                reportName = LocateVar(meNode->in_object.varRef)->GetName()->Cut(treeNameLength,-1);
                iBranchNames && & reportName;
                sortingOrder.lData[k] = k;
            }

            for (long k=0; k<sortingOrder2.lLength; k++) {
                node <nodeCoord> * meNode = (node <nodeCoord>*)sortingOrder2(k);
                reportName = LocateVar(meNode->in_object.varRef)->GetName()->Cut(treeNameLength,-1);
                leafNames && & reportName;
                sortingOrder2.lData[k] = k;
            }

            SortLists (&leafNames,   &sortingOrder2);
            SortLists (&iBranchNames,&sortingOrder);

            reportName = _String ("Branch Length Distribution for ") & *me->GetName();
            //vNames    << leafNames;
            //vNames    << iBranchNames;

            dataMatrix = new _Matrix (leafNames.lLength+iBranchNames.lLength, 1, false, true);
            checkPointer (dataMatrix);

            _String       rowLabels (128L, true);

            rowLabels << "Branch";

            for (long k=0; k<sortingOrder2.lLength; k++) {
                dataMatrix->theData[k] = ((_Constant*)leafLengths(sortingOrder2.lData[k]))->Value();
                rowLabels << ';';
                rowLabels << (_String*)leafNames(k);
            }

            for (long k=0; k<sortingOrder.lLength; k++) {
                dataMatrix->theData[k+sortingOrder2.lLength] = ((_Constant*)iBranchLengths(sortingOrder.lData[k]))->Value();
                rowLabels << ';';
                rowLabels << (_String*)iBranchNames(k);
            }

            rowLabels.Finalize();

            vNames.AppendNewInstance(new _String ("Length"));
            vNames && & rowLabels;
        } else { // pairwise distances
            node<nodeCoord>* traversalTop;

            if (currentSelection.lLength) {
                if (selectionTop&&(currentSelection.lLength>1)) {
                    traversalTop = selectionTop;
                } else {
                    reportName = "HyPhy can only compute pairwise distances on a complete subtree (or the entire tree)";
                    ProblemReport (reportName);
                    return;
                }
            } else {
                traversalTop = coordTree;
            }

            _SimpleList      indexer,
                             indexer2;

            _List            leafNames;

            long             leafCounter  = 0,
                             leafCounter2 = 0;

            node<nodeCoord>* meNode = NodeTraverser(traversalTop);

            while (meNode) {
                if (meNode->get_num_nodes() == 0) {
                    _String   tstr (*LocateVar(meNode->in_object.varRef)->GetName(),treeNameLength,-1);
                    leafNames && & tstr;
                    indexer<< leafCounter++;
                }
                meNode->in_object.auxL = 0;
                meNode = NodeTraverser((node<nodeCoord>*)NULL);
            }

            _Matrix * unsortedDM = nil;
            checkPointer (unsortedDM = new _Matrix (leafCounter, leafCounter, false, true));

            meNode = NodeTraverser(traversalTop);
            while (meNode != traversalTop) {
                node<nodeCoord>* parNode = meNode->parent;
                _SimpleList    * pInfo   = (_SimpleList*)parNode->in_object.auxL;
                if (!pInfo) {
                    checkPointer (pInfo = (_SimpleList*)(parNode->in_object.auxL = (long)new _SimpleList (leafCounter,0,0)));
                }


                if (meNode->get_num_nodes() == 0) {
                    pInfo->lData[leafCounter2] = 1;
                    for (long nid = 0; nid < leafCounter; nid++)
                        if (nid!=leafCounter2) {
                            unsortedDM->theData[nid*leafCounter + leafCounter2] += meNode->in_object.bL;
                            unsortedDM->theData[leafCounter2*leafCounter + nid] += meNode->in_object.bL;
                            /*char buffer [255];
                            snprintf (buffer, sizeof(buffer), "Adding %g %d:%d:%d (%g,%g) to %s %s path from leaf %s\n", meNode->in_object.bL,nid,leafCounter2,leafCounter,unsortedDM->theData[nid*leafCounter + leafCounter2],unsortedDM->theData[leafCounter2*leafCounter + nid],
                                                                                             ((_String*)leafNames(nid))->sData,((_String*)leafNames(leafCounter2))->sData,((_String*)leafNames(leafCounter2))->sData);
                            BufferToConsole (buffer);*/
                        }

                    leafCounter2 ++;
                } else {
                    _SimpleList inL,
                                outL,
                                *myL = (_SimpleList*)meNode->in_object.auxL;


                    for (long nid = 0; nid < leafCounter; nid++)
                        if (myL->lData[nid]) {
                            inL << nid;
                            pInfo->lData[nid] = 1;
                        } else {
                            outL << nid;
                        }


                    DeleteObject (myL);
                    meNode->in_object.auxL = 0;

                    if (meNode->in_object.bL>=0.0) {
                        for (long nid = 0; nid < inL.lLength; nid++) {
                            long idx1 = inL.lData[nid];
                            for (long nid2 = 0; nid2 < outL.lLength; nid2++) {
                                long idx2 = outL.lData[nid2];
                                unsortedDM->theData[idx1*leafCounter + idx2] += meNode->in_object.bL;
                                unsortedDM->theData[idx2*leafCounter + idx1] += meNode->in_object.bL;
                                /*char buffer [255];
                                snprintf (buffer, sizeof(buffer), "Adding %g %d:%d:%d  (%g,%g) to %s %s path from leaf %s\n", meNode->in_object.bL,idx1,idx2,leafCounter,
                                                                                                  unsortedDM->theData[idx1*leafCounter + idx2],
                                                                                                  unsortedDM->theData[idx2*leafCounter + idx1],
                                                                                                  ((_String*)leafNames(idx1))->sData,((_String*)leafNames(idx2))->sData,meNode->in_object.branchName.sData);
                                BufferToConsole (buffer);*/
                            }
                        }
                    }
                }
                meNode = NodeTraverser((node<nodeCoord>*)NULL);
            }
            DeleteObject ((BaseRef)meNode->in_object.auxL);
            meNode->in_object.auxL = 0;

            /*_List          pairwiseDistances;
                             // list of of lists
                             // each list has 3 constants
                             // 1: 1st leaf ID in leafID
                             // 2: 2nd leaf ID in leafID
                             // 3: the length of the path between leaves from entries 1 and 2

            while (meNode)
            {
                if (meNode->get_num_nodes())
                {
                    _List * descInfo = new _List;
                    checkPointer (descInfo);

                    for (long k=1; k<=meNode->get_num_nodes(); k++)
                    {
                        node<nodeCoord>* meChild = meNode->go_down(k);
                        if (meChild->get_num_nodes())
                        {
                            _List * childInfo = (_List*)meChild->in_object.auxL;

                            for (long jj = 0; jj < childInfo->lLength; jj++)
                            {
                                _List    * ccInfo      = (_List*)(*childInfo)(jj);
                                _Constant* pathToChild = (_Constant*)(*ccInfo)(1),
                                         * thisLeafID  = (_Constant*)(*ccInfo)(0);

                                _Parameter pathToMe = pathToChild->Value();

                                for (long ii = 0; ii < descInfo->lLength; ii++)
                                {
                                    _List distanceEntry,
                                         *previousEntry = (_List*)(*descInfo)(ii);

                                    distanceEntry << thisLeafID;
                                    distanceEntry << (*previousEntry) (0);
                                    distanceEntry && & _Constant (((_Constant*)(*previousEntry)(1))->Value()+pathToMe);

                                    pairwiseDistances && & distanceEntry;
                                }
                                pathToChild->SetValue(pathToMe+meChild->in_object.bL);
                            }

                            (*descInfo) << *childInfo;
                            DeleteObject  (childInfo);
                        }
                        else
                        {
                            for (long jj = 0; jj < descInfo->lLength; jj++)
                            {
                                _List distanceEntry,
                                     *previousEntry = (_List*)(*descInfo)(jj);

                                distanceEntry << (*previousEntry)(0);
                                distanceEntry && & _Constant (leafCounter);
                                distanceEntry && & _Constant (((_Constant*)(*previousEntry)(1))->Value()+meChild->in_object.bL);

                                pairwiseDistances && & distanceEntry;
                            }

                            _List * childInfo = new _List;
                            checkPointer (childInfo);
                            (*childInfo) && & _Constant (leafCounter);
                            (*childInfo) && & _Constant (meChild->in_object.bL);
                            (*descInfo) << childInfo;
                            DeleteObject (childInfo);
                            _String   tstr (*LocateVar(meChild->in_object.varRef)->GetName(),treeNameLength,-1);
                            leafNames && & tstr;
                            indexer<< leafCounter++;
                        }
                    }
                    meNode->in_object.auxL = (long)descInfo;

                }

                if (meNode == traversalTop)
                    break;

                meNode = NodeTraverser((node<nodeCoord>*)NULL);
            }

            DeleteObject ((_List*)meNode->in_object.auxL);*/

            indexer2.Duplicate (&indexer);
            SortLists    (&leafNames, &indexer2);
            SortLists    (&indexer2, &indexer);

            dataMatrix      = new _Matrix (leafNames.lLength, leafNames.lLength, false, true);
            checkPointer    (dataMatrix);

            /*for (long idx = 0; idx < pairwiseDistances.lLength; idx++)
            {
                _List * distanceEntry = (_List*) pairwiseDistances (idx);
                long    rIdx = indexer.lData[(long) ((_Constant*)(*distanceEntry)(0))->Value()],
                        cIdx = indexer.lData[(long) ((_Constant*)(*distanceEntry)(1))->Value()];

                _Parameter value = ((_Constant*)(*distanceEntry)(2))->Value();

                dataMatrix->theData[rIdx*leafNames.lLength+cIdx] = value;
                dataMatrix->theData[cIdx*leafNames.lLength+rIdx] = value;
            }*/

            for (long idx = 0; idx < leafCounter; idx ++) {
                long midx = indexer.lData[idx];
                for (long idx2 = idx + 1; idx2 < leafCounter; idx2++) {
                    long midx2 = indexer.lData[idx2];
                    dataMatrix->theData [midx2*leafCounter+midx] =
                        dataMatrix->theData [midx*leafCounter+midx2] =
                            unsortedDM->theData [idx*leafCounter+idx2];
                }
            }

            DeleteObject (unsortedDM);

            _String       rowLabels (128L, true);

            rowLabels << "Leaf";

            for (long k=0; k<leafNames.lLength; k++) {
                rowLabels << ';';
                rowLabels << (_String*)leafNames(k);
                vNames    << (_String*)leafNames(k);
            }

            rowLabels.Finalize();
            vNames && & rowLabels;
            reportName = _String ("Pairwise distances for ") & *me->GetName();

        }
        _HYChartWindow * reportChart = (_HYChartWindow*) FindWindowByNameAndOpen (reportName);
        if (reportChart) {
            reportChart->SetTable(vNames, *dataMatrix);
        } else {
            reportChart = new _HYChartWindow (reportName, vNames, *dataMatrix, nil);
        }
        reportChart->BringToFront();
        DeleteObject (dataMatrix);
    }
}

//__________________________________________________________________

_String         parameterLocal          ("Local"),
                parameterGlobal           ("Global"),
                parameterUE               ("Expression"),
                parameterConstrained  ("Constrained"),
                parameterMV               ("Multiple Values"),
                parameterAV               ("Assigned Value"),
                parameterBL               ("Branch Length"),
                parameterP                ("Parameter"),
                parameterUND          ("To be defined");

//__________________________________________________________________
_HYNodeInfoDialog::_HYNodeInfoDialog     (_String& n, long* ms, bool is, _TheTree* tr, _SimpleList* ts, bool* us, bool* ifl, bool* res,_HYTreePanel* pref)
    :_HYTWindow (n,false,true,(Ptr)pref)
{
    modelSelection = ms;
    isSingle       = is;
    treeRef        = tr;
    updateSCV      = us;
    incFlag        = ifl;
    tSel           = ts;
    result         = res;
    iNodeName      = &n;
    parentRef      = pref;

    _HYColor       bgRGB = GetDialogBackgroundColor ();

    _HYRect        canvasSettings = {30,100,30,100,HY_COMPONENT_NO_SCROLL|HY_COMPONENT_TRANSP_BG};

    _HYLabel*      l1 = new _HYLabel   (canvasSettings, GetOSWindowData());
    checkPointer   (l1);

    canvasSettings.width = HY_COMPONENT_BORDER_T|HY_COMPONENT_TRANSP_BG;
    _HYLabel*      l3 = new _HYLabel   (canvasSettings, GetOSWindowData());
    checkPointer   (l3);
    canvasSettings.width = HY_COMPONENT_NO_SCROLL|HY_COMPONENT_TRANSP_BG;

    if (!isSingle) {
        canvasSettings.width = HY_COMPONENT_BORDER_B|HY_COMPONENT_TRANSP_BG;
        nodeName = "Selected Branches";
    } else {
        nodeName = n.Cut(n.FindBackwards(".",0,-1)+1,-1);
    }

    _HYLabel*       l2 = new _HYLabel   (canvasSettings, GetOSWindowData());
    checkPointer   (l2);

    canvasSettings.left = canvasSettings.right = 200;
    _HYPullDown*    p1 = new _HYPullDown(canvasSettings, GetOSWindowData());
    checkPointer   (p1);


    canvasSettings.width = HY_COMPONENT_BORDER_T|HY_COMPONENT_TRANSP_BG;
    _HYTextBox*     tb2 = new _HYTextBox (canvasSettings, GetOSWindowData());
    checkPointer    (tb2);
    _HYLabel*       l5  = new _HYLabel (canvasSettings, GetOSWindowData());
    checkPointer    (l5);

    canvasSettings.width = HY_COMPONENT_NO_SCROLL|HY_COMPONENT_TRANSP_BG;

    canvasSettings.left = canvasSettings.right = 180;

    _HYTextBox*     tb = new _HYTextBox (canvasSettings, GetOSWindowData());
    checkPointer   (tb);

    canvasSettings.left = canvasSettings.right  = 20;
    canvasSettings.top  = canvasSettings.bottom = 30;

    _HYButtonBar*   bb = new _HYButtonBar (canvasSettings, GetOSWindowData());
    checkPointer   (bb);

    canvasSettings.top  = canvasSettings.bottom = 35;
    _HYLabel*       l4 = new _HYLabel (canvasSettings, GetOSWindowData());
    checkPointer   (l4);

    canvasSettings.left  = canvasSettings.right  = 300;
    canvasSettings.top   = canvasSettings.bottom = 180;
    canvasSettings.width = HY_COMPONENT_V_SCROLL|HY_COMPONENT_TRANSP_BG;

    _HYTable*       pl = new _HYTable    (canvasSettings, GetOSWindowData(), 1,3, 100, 180, HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT);
    checkPointer   (pl);

    pl->selectionType = HY_TABLE_SEL_ROWS|HY_TABLE_FOCUSABLE|HY_TABLE_DONT_SIZE|HY_TABLE_NODRAG_SELECTION|HY_TABLE_SINGLE_SELECTION;
    pl->SetColumnSpacing (0,15,false);
    pl->SetColumnSpacing (1,-15,false);

    canvasSettings.top   = canvasSettings.bottom = 20;
    canvasSettings.width = HY_COMPONENT_NO_SCROLL|HY_COMPONENT_TRANSP_BG;
    _HYTable*       ph   = new _HYTable    (canvasSettings, GetOSWindowData(), 1,3, 100, 20,
                                            HY_TABLE_BEVELED|HY_TABLE_STATIC_TEXT);
    checkPointer   (ph);
    ph->SetColumnSpacing (0,15,false);
    ph->SetColumnSpacing (1,-15,false);

    canvasSettings.width = HY_COMPONENT_BORDER_B|HY_COMPONENT_TRANSP_BG;
    canvasSettings.top   = canvasSettings.bottom = 25;
    _HYCheckbox*    cb   = new _HYCheckbox (canvasSettings, GetOSWindowData());
    checkPointer   (cb);

    canvasSettings.top    = canvasSettings.bottom = 35;
    canvasSettings.left   = canvasSettings.right  = 100;

    canvasSettings.width = HY_COMPONENT_NO_SCROLL|HY_COMPONENT_TRANSP_BG;

    _HYButtonBar*   bb2  = new _HYButtonBar (canvasSettings, GetOSWindowData());
    checkPointer   (bb2);

    canvasSettings.left   = canvasSettings.right  = 100;
    _HYButton*      obtn   = new _HYButton  (canvasSettings, GetOSWindowData());
    checkPointer    (obtn);

    canvasSettings.left   = canvasSettings.right  = 80;
    _HYButton*      cbtn   = new _HYButton  (canvasSettings, GetOSWindowData());
    checkPointer    (cbtn);


    bb->SetMessageRecipient (this);
    bb2->SetMessageRecipient (this);
    pl->SetMessageRecipient (this);
    //p1->SetMessageRecipient (this);
    ph->SetMessageRecipient (this);
    cb->SetMessageRecipient (this);
    cbtn->SetMessageRecipient (this);
    obtn->SetMessageRecipient (this);
    tb->SetMessageRecipient (this);
    tb2->SetMessageRecipient (this);

    SetTitle (_String("Branch Info for \"")&nodeName&"\"");

    AddObject   (pl);  //0
    AddObject   (ph);  //1
    AddObject   (bb);  //2
    AddObject   (bb2); //3
    AddObject   (l1);  //4
    AddObject   (cb);  //5
    AddObject   (obtn);//6
    AddObject   (cbtn);//7
    AddObject   (tb);  //8
    AddObject   (l2);  //9
    AddObject   (p1);  //10
    AddObject   (l4);  //11
    AddObject   (l3);  //12
    AddObject   (l5);  //13
    AddObject   (tb2); //14

    if (isSingle) {
        SetTableDimensions (6,4);
        SetCell        (0,0,l1);
        SetCell        (0,1,tb);
        SetCell        (0,2,tb);
        SetCell        (0,3,bb);
        SetCell        (1,0,cb);
        SetCell        (1,1,cb);
        SetCell        (1,2,cb);
        SetCell        (1,3,cb);
        SetCell        (2,0,ph);
        SetCell        (2,1,ph);
        SetCell        (2,2,ph);
        SetCell        (2,3,ph);
        SetCell        (3,0,pl);
        SetCell        (3,1,pl);
        SetCell        (3,2,pl);
        SetCell        (3,3,pl);
        SetCell        (4,0,l3);
        SetCell        (4,1,l5);
        SetCell        (4,2,l5);
        SetCell        (4,3,l5);
        SetCell        (5,0,bb2);
        SetCell        (5,1,obtn);
        SetCell        (5,2,cbtn);
        SetCell        (5,3,l4);
    } else {
        SetTableDimensions (4,4);
        SetCell        (0,0,ph);
        SetCell        (0,1,ph);
        SetCell        (0,2,ph);
        SetCell        (0,3,ph);
        SetCell        (1,0,pl);
        SetCell        (1,1,pl);
        SetCell        (1,2,pl);
        SetCell        (1,3,pl);
        SetCell        (2,0,l3);
        SetCell        (2,1,l5);
        SetCell        (2,2,l5);
        SetCell        (2,3,l5);
        SetCell        (3,0,bb2);
        SetCell        (3,1,obtn);
        SetCell        (3,2,cbtn);
        SetCell        (3,3,l4);
#ifdef __HYPHY_GTK__
        tb->Deactivate();
#endif
    }

    /*  _HYTable *      hl      = new _HYTable (canvasSettings, GetOSWindowData(),validChoices->lLength,1,300,20,
                                                                HY_TABLE_STATIC_TEXT);*/
    l1->SetBackColor (bgRGB);
    l2->SetBackColor (bgRGB);
    l3->SetBackColor (bgRGB);
    l4->SetBackColor (bgRGB);
    l5->SetBackColor (bgRGB);
    p1->SetBackColor (bgRGB);
    tb->SetBackColor (bgRGB);
    tb2->SetBackColor (bgRGB);
    bb->SetBackColor (bgRGB);
    bb2->SetBackColor (bgRGB);
    cb->SetBackColor (bgRGB);
    obtn->SetBackColor (bgRGB);
    cbtn->SetBackColor (bgRGB);


    _HYFont         labelFont;
#ifdef __WINDOZE__
    labelFont.face  = "Arial";
    labelFont.size  = 14;
#else
    labelFont.face  = "System Font";
    labelFont.size  = 12;
#endif
    labelFont.style = HY_FONT_PLAIN;

    obtn->SetFont (labelFont);
    cbtn->SetFont (labelFont);
    l1->SetFont   (labelFont);
    l1->SetAlignFlags (HY_ALIGN_LEFT);
    l1->SetText   ("Branch ID:");
    l2->SetFont   (labelFont);
    l2->SetAlignFlags (HY_ALIGN_LEFT);
    l2->SetText   ("Model:");
    l3->SetFont   (labelFont);
    l3->SetAlignFlags (HY_ALIGN_LEFT);
    cb->SetFont   (labelFont);
    cb->SetAlignFlags (HY_ALIGN_LEFT);
    cb->SetText   (_String("Mark \"")&nodeName&"\" as incomplete");
    tb->EnableTextEdit (true);
    tb->SetText        (nodeName);
    p1->AddMenuItem     ("None",-1);
    p1->SetAlignFlags   (HY_ALIGN_LEFT);
    bb2->SetAlignFlags  (HY_ALIGN_LEFT);

    long counter;

#ifdef __HYPHY_GTK__
    gtk_widget_hide(p1->theMenu);
#endif

    /*for (counter = 0; counter < modelNames.lLength; counter++)
        p1->AddMenuItem (*(_String*)modelNames(counter),-1);*/

    if (*modelSelection<0) {
        *modelSelection = 0;
    } else {
        (*modelSelection)++;
    }

    //p1->ChangeSelection(*modelSelection,false);

    obtn->SetAlignFlags (HY_ALIGN_RIGHT);

    if (isSingle) {
        _String tryName,
                namePrefix = *treeRef->GetName() & '.';

        for (counter = 0; counter < dataSetList.lLength; counter++) {
            if (counter&&(allowableDataNames.lLength)) {
                allowableDataNames<< &menuSeparator;
            }

            _DataSet* thisSet = (_DataSet*)dataSetList (counter);
            _String*  spName;
            _List     thisList;

            for (long whatW = 0; whatW < thisSet->NoOfSpecies(); whatW++) {
                spName = (_String*)(thisSet->GetNames()(whatW));
                tryName = namePrefix&*spName;
                if (LocateVarByName(tryName)<0) {
                    long iType = thisList.BinaryFindObject (spName);
                    if (iType<0) {
                        if (thisList.lLength) {
                            iType = -iType-2;
                            if (*spName>*(_String*)thisList(iType)) {
                                iType++;
                            }
                            thisList.InsertElement (spName,iType,false);
                        } else {
                            thisList<<spName;
                        }
                    }
                }
            }
            allowableDataNames<<thisList;
        }
    }


#ifdef __WINDOZE__
    labelFont.face = "MS Sans Serif";
    labelFont.size = 12;
#else
    labelFont.face = "Geneva";
    labelFont.size = 10;
#endif
    tb->SetFont (labelFont);
    tb2->SetFont (labelFont);

    obtn->SetText         ("  OK  ");
    obtn->SetButtonKind   (HY_BUTTON_OK);
    cbtn->SetText         (" Cancel ");
    cbtn->SetButtonKind   (HY_BUTTON_CANCEL);

#ifdef __WINDOZE__
    labelFont.face = "Arial";
    labelFont.size = 14;
#else
#ifdef __MAC__
    labelFont.face = "Times";
    labelFont.size = 12;
#else
    labelFont.face = _HY_SANS_FONT;
    labelFont.size = 12;
#endif
#endif

    pl->SetFont (labelFont);
    ph->SetFont (labelFont);

    _String     str ("Parameter ID");

    ph->SetCellData (&str,0,0,HY_TABLE_STATIC_TEXT|HY_TABLE_BOLD|HY_TABLE_BEVELED|HY_TABLE_CANTSELECT,true);
    str = "Type";
    ph->SetCellData (&str,0,1,HY_TABLE_STATIC_TEXT|HY_TABLE_BOLD|HY_TABLE_BEVELED|HY_TABLE_CANTSELECT,true);
    str = "Value";
    ph->SetCellData (&str,0,2,HY_TABLE_STATIC_TEXT|HY_TABLE_BOLD|HY_TABLE_BEVELED|HY_TABLE_CANTSELECT,true);

    str = "Available sequence names from the data filter";
    bb->AddButton (ProcureIconResource(4011),&str);
    bb->MarkAsPullDown (0,true);

    str = "Add User Expression";
    bb2->AddButton (ProcureIconResource(6001),&str);
    str = "Kill User Expression";
    bb2->AddButton (ProcureIconResource(6003),&str);
    str = "Check and Set Formula";
    bb2->AddButton (ProcureIconResource(7020),&str);

    bb->SetButtonDim  (13);
    bb2->SetButtonDim (16);
    bb2->EnableButton (1,false);
    bb2->EnableButton (2,false);
    bb->EnableButton  (0, allowableDataNames.lLength);
    cb->SetState      (*ifl);

    _HYRect     dim = MinMaxWindowDimensions();

    SetWindowRectangle  (0,0,dim.bottom,dim.right);
    CenterWindow        (this);
    BuildModelParameters();

    validTE          =  true;
    formulaBox       =  false;

    tb2->SetVisibleSize (l5->rel);
    tb2->Deactivate();
    componentL.lData[14] = componentL.lData[13];
    componentR.lData[14] = componentR.lData[13];
    componentT.lData[14] = componentT.lData[13];
    componentB.lData[14] = componentB.lData[13];

    if (isSingle) {
        keyboardFocusChain << 8;
    }

    keyboardFocusChain << 0;
    taintFla   = false;

    DeleteObject    (pl);  //0
    DeleteObject    (ph);  //1
    DeleteObject    (bb);  //2
    DeleteObject    (bb2); //3
    DeleteObject    (l1);  //4
    DeleteObject    (cb);  //5
    DeleteObject    (obtn);//6
    DeleteObject    (cbtn);//7
    DeleteObject    (tb);  //8
    DeleteObject    (l2);  //9
    DeleteObject    (p1);  //10
    DeleteObject    (l4);  //11
    DeleteObject    (l3);  //12
    DeleteObject    (l5);  //13
    DeleteObject    (tb2); //14
}

//__________________________________________________________________

void    _HYNodeInfoDialog::ToggleFormulaBox     (bool onOff, _String value)
{
    _HYLabel  *  l3  = (_HYLabel*)GetObject (12),
                 *  l5  = (_HYLabel*)GetObject (13);

    _HYTextBox*  tb2 = (_HYTextBox*)GetObject(14);

    long         row = isSingle?4:2;

    if (onOff) {
        if (!formulaBox) {
            SetCell (row,1,tb2);
            SetCell (row,2,tb2);
            SetCell (row,3,tb2);
            l3->SetText ("Formula:");
            tb2->SetText (value,false);
            tb2->SetSelection (0,30000);
            tb2->Activate();
            tb2->EnableTextEdit(true);
            keyboardFocusChain << 14;
        } else {
            tb2->SetText (value,true);
            tb2->SetSelection (0,30000);
        }
    } else {
        if (formulaBox) {
            SetCell (row,1,l5);
            SetCell (row,2,l5);
            SetCell (row,3,l5);
            tb2->Deactivate();
            l5->_MarkForUpdate();
            l3->SetText (empty);
            keyboardFocusChain.Delete (keyboardFocusChain.lLength-1);
        }
    }
    formulaBox = onOff;
    taintFla   = false;
}
//__________________________________________________________________

void    _HYNodeInfoDialog::BuildModelParameters (void)
{
    _CalcNode       *thisCNode  = (_CalcNode*)FetchVar (LocateVarByName(*iNodeName));

    if (!thisCNode) {
        return ;
    }

    _HYTable        *pl         = (_HYTable*)GetObject (0);
    _HYButtonBar    *bb2        = (_HYButtonBar*)GetObject (3);
    long            counter;

    if (isSingle) {
        _SimpleList       indVars,
                          depVars,
                          globalVars;

        _Variable*        thisVar;
        _String           theMessage;

        while (((_String*)pl->GetCellData(0,0))->sLength) {
            DeleteParameterRow (0);
        }

        {
            _AVLList ia (&indVars),
                     da (&depVars),
                     ga (&globalVars);

            thisCNode->ScanForVariables (ia,da);
            thisCNode->ScanForDVariables (da,ia);
            thisCNode->ScanForGVariables (ga,da);

            ia.ReorderList();
            da.ReorderList();
            ga.ReorderList();
        }
        parameterFlags.Clear();
        parameterVRefs.Clear();

        long choice = indVars.lLength + depVars.lLength + globalVars.lLength;
        if (choice) {
            parameterValues.Clear();
            parameterVRefs.Clear();
            CreateMatrix (&parameterValues,choice,1,false,true,false);

            for (counter = 0; counter < globalVars.lLength; counter++) {
                thisVar    = LocateVar(globalVars.lData[counter]);
                theMessage = *thisVar->GetName();
                AddParameterRow (theMessage,parameterGlobal,empty,-1);
                parameterFlags<<1;
                parameterVRefs<<thisVar->theIndex;
            }
            for (counter = 0; counter < indVars.lLength; counter++) {
                thisVar = LocateVar(indVars.lData[counter]);
                theMessage = thisVar->GetName()->ShortenVarID (*thisCNode->GetName());
                AddParameterRow (theMessage,parameterLocal,empty,-1);
                parameterFlags<<0;
                parameterVRefs<<thisVar->theIndex;
            }
            for (counter = 0; counter < depVars.lLength; counter++) {
                thisVar = LocateVar(depVars.lData[counter]);
                if (thisVar->IsGlobal()) {
                    theMessage = *thisVar->GetName();
                    AddParameterRow (theMessage,parameterGlobal,empty,-1);
                    parameterFlags<<1;
                    parameterVRefs<<thisVar->theIndex;
                } else {
                    theMessage = thisVar->GetName()->ShortenVarID (*thisCNode->GetName());

                    if (thisCNode->IsModelVar(counter)) {
                        AddParameterRow (theMessage,parameterConstrained,empty,-1);
                        parameterFlags<<2;
                    } else {
                        AddParameterRow (theMessage,parameterUE,empty,-1);
                        parameterFlags<<3;
                    }
                    parameterVRefs<<thisVar->theIndex;
                }
            }
            UpdateParameterValues (true);
            bb2->EnableButton (0,true);
        } else {
            bb2->EnableButton (0,true);
            CreateMatrix (&parameterValues,1,1,false,true,false);
            AddParameterRow (parameterBL, parameterAV, empty,-1);
            parameterFlags<<1;
            parameterVRefs<<thisCNode->theIndex;
            UpdateParameterValues (true);
        }
    } else {
        for (counter = 0; counter < tSel->lLength; counter++) {
            nodeVRefs << ((node<nodeCoord>*)tSel->lData[counter])->in_object.varRef;
        }

        CompileListOfUserExpressions (nodeVRefs,sharedUEs,true);

        for (counter = 0; counter < sharedUEs.lLength; counter++) {
            _String theMessage = *(_String*)sharedUEs(counter);
            if (theMessage.sData[0]=='!') {
                ((_String*)sharedUEs(counter))->Trim(1,-1);
                AddParameterRow (*(_String*)sharedUEs(counter),parameterUE, parameterMV,-1);

                for (long iType = 0; iType<tSel->lLength; iType++) {
                    _String tBuffer = *LocateVar(nodeVRefs.lData[iType])->GetName()&'.'&*(_String*)sharedUEs(counter);
                    parameterFlags<<3;
#ifndef USE_AVL_NAMES
                    parameterVRefs<<variableReindex.lData[LocateVarByName(tBuffer)];
#else
                    parameterVRefs<<variableNames.GetXtra(LocateVarByName(tBuffer));
#endif
                }
            } else {
                AddParameterRow (*(_String*)sharedUEs(counter),parameterP, parameterMV,-1);

                for (long iType = 0; iType<tSel->lLength; iType++) {
                    _String tBuffer = *LocateVar(nodeVRefs.lData[iType])->GetName()&'.'&*(_String*)sharedUEs(counter);
                    parameterFlags<<2;
#ifndef USE_AVL_NAMES
                    parameterVRefs<<variableReindex.lData[LocateVarByName(tBuffer)];
#else
                    parameterVRefs<<variableNames.GetXtra(LocateVarByName(tBuffer));
#endif
                }
            }
        }
    }
}

//________________________________________________________

void    _HYNodeInfoDialog::UpdateParameterValues (bool force)
{
    _HYTable        *pl         = (_HYTable*)GetObject (0);
    bool redraw = force;
    for (long i=0; i<parameterVRefs.lLength; i++) {
        _Variable* thisVar = LocateVar (parameterVRefs.lData[i]);
        if (force||(thisVar->HasChanged())) {
            redraw = true;
            _PMathObj  varVal = thisVar->Compute();
            _String varValue (varVal->Value());

            pl->SetCellData (&varValue,i,2,pl->cellTypes[3*i+2],true);

            if (force&&(parameterFlags.lData[i]<3)) {
                parameterValues.theData[i] = varVal->Value();
            }
        }
    }
    if (redraw) {
        pl->_MarkColumnForUpdate(2);
    }
}

//__________________________________________________________________

void    _HYNodeInfoDialog::FitToHeight (long height)
{
    _HYTable       *pl      = (_HYTable*)GetObject (0);
    bool           hasPadding = (pl->cellTypes.lData[pl->cellTypes.lLength-1] & HY_TABLE_CANTSELECT);
    long           h = pl->verticalSpaces.lData[pl->verticalSpaces.lLength-(hasPadding?2:1)];

    if (hasPadding) {
        if (h<height) {
            pl->SetRowSpacing (pl->verticalSpaces.lLength-1,height-h-pl->GetRowSpacing (pl->verticalSpaces.lLength-1),false);
        } else {
            pl->DeleteRow (pl->verticalSpaces.lLength-1);
        }
    } else if (h<height) {
        pl->AddRow (-1,height-h,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT);
    }

    pl->SetVisibleSize (pl->rel);
}

//__________________________________________________________________

void    _HYNodeInfoDialog::DeleteParameterRow (long index)
{
    _HYTable       *pl      = (_HYTable*)GetObject (0);
    pl->DeleteRow (index);
    FitToHeight (pl->rel.bottom-pl->rel.top+1);
}

//__________________________________________________________________

bool    _HYNodeInfoDialog::SetNewParameterValue (long index, _String* val)
{
    _String      newValStr (*val);
    _Parameter   newVal = 0.0;
    newVal=ProcessNumericArgument(&newValStr,nil);
    bool         good   = false;
    if (numericalParameterSuccessFlag) {
        _Variable* thisVar = LocateVar (parameterVRefs.lData[index]);
        if (thisVar) {
            if ((newVal>=thisVar->GetLowerBound())&&(newVal<=thisVar->GetUpperBound())) {
                _Constant newV (newVal);
                thisVar->SetValue (&newV);
                good = true;
            } else {
                _String errMsg ("Valid values for this parameter must be in the interval [");
                errMsg = errMsg & thisVar->GetLowerBound() &','& thisVar->GetUpperBound() &"].";
                ProblemReport (errMsg,(Ptr)this);
            }
        }
    }
    return good;
}

//__________________________________________________________________

long    _HYNodeInfoDialog::AddParameterRow (_String& pname, _String& ptype, _String& pvalue, long where)
{
    _HYTable       *pl      = (_HYTable*)GetObject (0);
    bool           hasPadding = (pl->cellTypes.lData[pl->cellTypes.lLength-1] & HY_TABLE_CANTSELECT);
    if (where<0) {
        where = pl->verticalSpaces.lLength;
        if (hasPadding) {
            where--;
        }
    }

    long      ct1 = HY_TABLE_STATIC_TEXT,
              ct2 = HY_TABLE_STATIC_TEXT,
              ct3 = HY_TABLE_EDIT_TEXT;

    if (ptype.Equal (&parameterGlobal)) {
        ct1 |= HY_TABLE_BOLD;
        ct2 |= HY_TABLE_BOLD;
        ct3 |= HY_TABLE_BOLD;
    } else if (ptype.Equal (&parameterConstrained)||ptype.Equal (&parameterUE)) {
        if (ptype.Equal(&parameterUE)) {
            ct1 = HY_TABLE_EDIT_TEXT;
        }

        ct1 |= HY_TABLE_ITALIC;
        ct2 |= HY_TABLE_ITALIC;
        ct3 = HY_TABLE_ITALIC|HY_TABLE_STATIC_TEXT;
    } else if (pvalue.Equal (&parameterMV)) {
        ct3 = HY_TABLE_ITALIC|HY_TABLE_STATIC_TEXT;
    }


    pl->AddRow (where,20,HY_TABLE_STATIC_TEXT);
    pl->SetCellData (&pname,where,0,ct1,true);
    pl->SetCellData (&ptype,where,1,ct2,true);
    pl->SetCellData (&pvalue,where,2,ct3,true);
    FitToHeight (pl->rel.bottom-pl->rel.top+1);
    return where;
}

//________________________________________________________

bool    _HYNodeInfoDialog::SetNodeModel (_CalcNode* thisCNode, long modelID, bool isFirst)
{
  if (isFirst && !displayNotUndoable) {
    if (ProceedPromptWithCheck (notUndoableWarning,donotWarnAgain,displayNotUndoable)) {
      _SimpleList depVars;
      for (long iType=0; iType<tSel->lLength; iType++) {
        LocateVar (((node<nodeCoord>*)tSel->lData[iType])->in_object.varRef)->CompileListOfDependents(depVars);
      }
      
      if (depVars.lLength) {
        _String depVarList ("The following variables depend on the parameters you are about to delete: ");
        for (long i=0; i<depVars.lLength; i++) {
          if (i) {
            depVarList = depVarList & ", " ;
          }
          depVarList = depVarList & *LocateVar(depVars.lData[i])->GetName();
        }
        depVarList = depVarList & ((char)13)&"You are about to remove those dependencies.";
        if (!ProceedPrompt (depVarList,(Ptr)this)) {
          return false;
        }
      }
      *modelSelection = modelID+1;
    } else {
      return false;
    }
  }
  
  _String treeName (*thisCNode->GetName());
  treeName.Trim (0,treeName.Find('.')-1);
  
  if (modelID<0) {
    if (isFirst&&(!parentRef->KillLikeFunc (-1))) {
      return false;
    }
    
    DeleteVariable (*thisCNode->GetName(),false);
  } else {
    _String theName (*thisCNode->GetName()),
    modelName (*(_String*)modelNames(modelID)),
    dummy;
    
    theName.Trim (theName.Find('.')+1,-1);
    
    _CalcNode* travNode = treeRef->DepthWiseTraversal (true);
    
    while (travNode!=thisCNode) {
      travNode = treeRef->DepthWiseTraversal ();
    }
    
    DeleteVariable (*thisCNode->GetName(),true);
    treeRef->FinalizeNode (&treeRef->GetCurrentNode(),0,theName,modelName,dummy);
    thisCNode = (_CalcNode*)LocateVar (treeRef->GetCurrentNode().in_object);
    _SimpleList newVars;
    {
      _AVLList  nal (&newVars);
      thisCNode->ScanForVariables (nal,nal);
      nal.ReorderList();
    }
    thisCNode->SetCodeBase(thisCNode->GetModelDimension());
    _Constant   newVal (.25);
    for (long i=0; i<newVars.lLength; i++) {
      LocateVar(newVars.lData[i])->SetValue (&newVal);
    }
    
    for (long i = 0; i<likeFuncList.lLength; i++) {
      if (((_String*)likeFuncNamesList(i))->sLength) {
        _LikelihoodFunction * thisLF = (_LikelihoodFunction*)likeFuncList(i);
        if (thisLF->DependOnTree(treeName)>=0) {
          thisLF->RescanAllVariables();
        }
      }
    }
  }
  
  return true;
}

//__________________________________________________________________
bool    _HYNodeInfoDialog::ProcessEvent  (_HYEvent* e)
{
    bool            done = false;
    _String         firstArg;

    long            i,
                    f,
                    g,
                    k;

    _HYButtonBar    * bb  = (_HYButtonBar*)GetObject (2),
                      * bb2 = (_HYButtonBar*)GetObject (3);

    _HYTextBox      * tb  = (_HYTextBox*)  GetObject (8);
    _HYButton       * bok = (_HYButton*)   GetObject (6);
    _HYTable        * pl  = (_HYTable*)    GetObject (0);
    _HYCheckbox     * chb = (_HYCheckbox*) GetObject (5);

    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        g =  e->EventCode().Cut(f+1,-1).toNum();
        for (i=0; i<components.lLength; i++)
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }

        if (i==6) { // OK
            *result = true;
            pl->Deactivate();
            if (taintFla) {
                SetNewUserFormula();
            }

            if (isSingle) {
                *iNodeName = tb->GetText();
            }

            *incFlag = chb->GetState ();

            postWindowCloseEvent (GetID());
        } else if (i==7) { // Cancel
            *result = false;
            if (isSingle) {
                for (long counter = 0; counter < parameterVRefs.lLength; counter++) {
                    _Variable *thisVar = LocateVar (parameterVRefs.lData[counter]);
                    if (thisVar->IsIndependent()||thisVar->IsContainer()) {
                        if (!CheckEqual(thisVar->Value(),parameterValues.theData[counter])) {
                            _Constant newVal (parameterValues.theData[counter]);
                            thisVar->SetValue (&newVal);
                        }
                    }
                }
            }
            postWindowCloseEvent (GetID());
        } else if ((i==2)&&(g==0)) {
            int h,v;
            bb->GetButtonLoc (0,h,v,true);
            _String*tbs;
            tb->StoreText (tbs);
            firstArg = HandlePullDown (allowableDataNames,h,v,allowableDataNames.FindObject(tbs));
            DeleteObject (tbs);
            if (firstArg.sLength) {
                tb->SetText (firstArg);
                ProcessEvent(generateKeyboardFocusEvent  (tb->GetID()));
                ProcessEvent(generateTextEditChangeEvent (tb->GetID(),1));
            }
        }

        done = true;
    } else if (e->EventClass()==_hyTableChangeSelEvent) {
        k = e->EventCode().toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }

        if (i==0) {
            _SimpleList rs;
            pl->GetRowSelection (rs);
            if (rs.lLength<=1) {
                if (rs.lLength) {
                    done = parameterUE.Equal((_String*)pl->GetCellData(1,rs.lData[0]));
                } else {
                    done = false;
                }
                bb2->EnableButton (1,done);
                bb2->EnableButton (2,done);
                if (done) {
                    _Variable* thisVar = LocateVar (parameterVRefs.lData[rs.lData[0]*tSel->lLength]);
                    _String* flaString = thisVar->GetFormulaString(),
                             infoBoxString,
                             theMessage;

                    _CalcNode* thisCNode;

                    if (flaString->sLength==0) {
                        if (parameterFlags.lData[rs.lData[0]*tSel->lLength]>2) {
                            infoBoxString = "Enter your expression";
                        } else {
                            infoBoxString = empty;
                        }
                    } else {
                        if (isSingle) {
                            thisCNode  = (_CalcNode*)FetchVar (LocateVarByName(*iNodeName));
                        } else {
                            thisCNode=(_CalcNode*)LocateVar (nodeVRefs.lData[0]);
                        }

                        theMessage = *thisCNode->GetName()&'.';
                        if (parameterFlags.lData[rs.lData[0]*tSel->lLength]>2) {
                            infoBoxString = *flaString;
                        } else {
                            infoBoxString = *thisVar->GetName()&_String(":=")&*flaString;
                        }
                        infoBoxString = infoBoxString.Replace(theMessage,"",true);
                        theMessage = theMessage.Cut (0,theMessage.Find('.'));
                        infoBoxString = infoBoxString.Replace(theMessage,"",true);
                    }
                    ToggleFormulaBox (true, infoBoxString);

                } else {
                    ToggleFormulaBox (false, empty);
                }
            }
        }
        done = true;
    } else if (e->EventClass()==_hyTextEditChange) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        g =  e->EventCode().Cut(f+1,-1).toNum();
        for (i=0; i<components.lLength; i++)
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }

        if (g) {
            if (i==8) {
                _String     curVal (tb->GetText()),
                            temp = *treeRef->GetName() & '.' & curVal;

                bool isGood = curVal.sLength&&temp.IsValidIdentifier ();
                if (isGood) {
                    if (!temp.Equal(iNodeName)) {
                        isGood = (LocateVarByName (temp)<0);
                    }
                }

                if (isGood&&(!validTE)) {
                    tb->SetForeColor ((_HYColor) {
                        0,0,0
                    });
                    bok->EnableButton (true);
                } else if (!isGood&&(validTE)) {
                    tb->SetForeColor ((_HYColor) {
                        255,0,0
                    });
                    bok->EnableButton (false);
                }

                if (isGood) {
                    nodeName = curVal;
                }

                validTE = isGood;
            } else if (i==14) {
                taintFla = true;
            }
        }

        done = true;
    } else if (e->EventClass()==_hyTableResizeCEvent ) {
        f = e->EventCode().Find(',',0,-1);
        k = e->EventCode().Find(',',f+1,-1);
        long g = e->EventCode().Cut (0,f-1).toNum();
        f = e->EventCode().Cut (f+1,k-1).toNum();
        k = e->EventCode().Cut (k+1,-1).toNum();

        for (i=0; i<components.lLength; i++)
            if (((_HYGuiObject*)components(i))->MatchID(g)) {
                break;
            }

        if (i==1) {
            pl->SetColumnSpacing (f,k,true);
        }

        done = true;
    } else if (e->EventClass()==_hyTableEditCellEvent ) {
        f = e->EventCode().Find(',',0,-1);
        k = e->EventCode().Cut (0,f-1).toNum();
        f = e->EventCode().Cut (f+1,-1).toNum();

        for (i=0; i<components.lLength; i++)
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }

        if (i==0) {
            g = f/pl->horizontalSpaces.lLength;
            f%=pl->horizontalSpaces.lLength;
            if (f==2) {
                UpdateParameterValues (!SetNewParameterValue (g,(_String*)pl->GetCellData (2,g)));
            } else if (f==0) {
                SetNewParameterName (*(_String*)pl->GetCellData (0,g),g);
            }

        }

        done = true;
    } else if (e->EventClass()==_hyMenuSelChangeEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k =  firstArg.toNum();
        g =  e->EventCode().Cut(f+1,-1).toNum();
        for (i=0; i<components.lLength; i++)
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }

        if ((i==10)&&(g!=*modelSelection))
            for (i=0; i<tSel->lLength; i++) {
                _CalcNode* tempNode = (_CalcNode*)LocateVar
                                      (((node<nodeCoord>*)tSel->lData[i])->in_object.varRef);
                if (!SetNodeModel (tempNode,g-1,i==0)) {
                    ((_HYPullDown*)GetObject(10))->ChangeSelection (*modelSelection,false);
                    break;
                }
            }
        done = true;
    }
    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k =  firstArg.toNum();
        g =  e->EventCode().Cut(f+1,-1).toNum();
        for (i=0; i<components.lLength; i++)
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }

        if (i==3) {
            *updateSCV = true;

            _CalcNode* thisCNode;

            if (g==0) { // add user expression
                _String theMessage ("user_param"),
                        tBuffer;

                f = 0;

                if (isSingle) {
                    thisCNode  = (_CalcNode*)FetchVar (LocateVarByName(*iNodeName));
                    f = thisCNode->CheckAndAddUserExpression (theMessage);
                    parameterVRefs << f;
                    parameterFlags << 3;
                } else {
                    for (g = 0; g<tSel->lLength; g++) {
                        k = ((_CalcNode*)LocateVar (((node<nodeCoord>*)tSel->lData[g])->in_object.varRef))->CheckAndAddUserExpression(theMessage,-1);
                        if (k>f) {
                            f=k;
                        }
                    }

                    for (g = 0; g<tSel->lLength; g++) {
                        thisCNode=(_CalcNode*)LocateVar (((node<nodeCoord>*)tSel->lData[g])->in_object.varRef);
                        tBuffer = theMessage;
                        k = thisCNode->CheckAndAddUserExpression(tBuffer,f);
                        parameterVRefs << k;
                        parameterFlags << 3;
                    }

                    if (tSel->lLength>1) {
                        sharedUEs && &tBuffer;
                    }

                    theMessage = tBuffer;
                }

                k = AddParameterRow (theMessage,parameterUE,parameterUND,-1);
                pl->_MarkForUpdate();
                _SimpleList nS;
                nS << k;
                pl->SetRowSelection (nS);
            } else {
                if (g==1) {
                    _SimpleList rs2;
                    pl->GetRowSelection (rs2);

                    f = rs2.lData[0] * tSel->lLength;

                    if (isSingle) {
                        thisCNode  = (_CalcNode*)FetchVar (LocateVarByName(*iNodeName));
                        thisCNode->KillUserExpression (parameterVRefs.lData[rs2.lData[0]]);
                    } else {
                        for (k = 0; k<tSel->lLength; k++) {
                            thisCNode=(_CalcNode*)LocateVar (nodeVRefs.lData[k]);
                            thisCNode->KillUserExpression (parameterVRefs.lData[f+k]);
                        }
                        sharedUEs.Delete (rs2.lData[0]);
                    }

                    for (k=0; k<tSel->lLength; k++) {
                        parameterVRefs.Delete (f);
                        parameterFlags.Delete (f);
                    }

                    DeleteParameterRow (rs2.lData[0]);
                    _SimpleList         nRS;
                    if (rs2.lData[0]>0) {
                        nRS << rs2.lData[0]-1;
                    } else {
                        nRS << 0;
                    }
                    pl->SetRowSelection (nRS);
                    pl->_MarkForUpdate();
                } else {
                    SetNewUserFormula();
                }

            }
        }
        done = true;
    }

    if (done) {
        DeleteObject (e);
        return true;
    }
    return _HYTWindow::ProcessEvent(e);
}

//________________________________________________________

void    _HYNodeInfoDialog::SetNewUserFormula (void)
{
    _HYTextBox  * flaBox = (_HYTextBox*) GetObject (14);
    _HYTable    * pl     = (_HYTable*)   GetObject (0);

    _String     newFla   = flaBox->GetText(),
                newValStr ;
    _SimpleList rSel;

    pl->GetRowSelection (rSel);
    if (rSel.lLength!=1) {
        return;
    }

    long        k = rSel.lData[0],
                j;

    _CalcNode*  thisCNode;

    if (nodeVRefs.lLength) {
        _Variable* thisVar;
        for (j=0; j<nodeVRefs.lLength; j++) {
            thisCNode = (_CalcNode*)LocateVar(nodeVRefs.lData[j]);
            _Formula   newUF (newFla,thisCNode);

            if (!newUF.IsEmpty()) {
                thisVar = LocateVar (parameterVRefs.lData[k*tSel->lLength+j]);
                if (thisVar) {
                    thisVar->SetFormula(newUF);
                }
            }
        }
        newValStr = parameterMV;
    } else {
        thisCNode  = (_CalcNode*)FetchVar (LocateVarByName(*iNodeName));

        _Formula   newUF (newFla,thisCNode);
        if (newUF.IsEmpty()) {
            return;
        }

        _Variable* thisVar = LocateVar (parameterVRefs.lData[k]);
        if (thisVar) {
            thisVar->SetFormula(newUF);
            newValStr = _String (thisVar->Compute()->Value());

            _String*  newFlaString = thisVar->GetFormulaString();
            flaBox->SetText (*newFlaString,false);
            DeleteObject (newFlaString);
        }
    }

    j = k*pl->horizontalSpaces.lLength + 2;
    pl->SetCellData (&newValStr,k,2,pl->cellTypes.lData[j],true);
    pl->_MarkCellForUpdate (j);
}

//________________________________________________________

void    _HYNodeInfoDialog ::SetNewParameterName (_String newValStr, long row)
{
    _HYTable*       pl = (_HYTable*)GetObject (0);
    _CalcNode*      thisCNode;
    long            k = row*tSel->lLength;
    _String         oldName;
    _Variable*      thisVar;

    bool       restore = true;

    if (isSingle) {
        thisCNode  = (_CalcNode*)FetchVar (LocateVarByName(*iNodeName));
    } else {
        thisCNode = (_CalcNode*)LocateVar(nodeVRefs.lData[0]);
    }

    thisVar = LocateVar (parameterVRefs.lData[k]);
    oldName = *(thisVar->GetName());

    if (newValStr.IsValidIdentifier()) {
        if (isSingle) {
            if (thisVar) {
                newValStr = *thisCNode->GetName()&'.'&newValStr;
                if (!thisVar->GetName()->Equal(&newValStr)) {
                    oldName = *thisVar->GetName();
                    if (LocateVarByName(newValStr)>=0) {
                        newValStr = newValStr & " is already in use. Select another name.";
                        ProblemReport (newValStr,(Ptr)this);
                    } else {
                        RenameVariable (thisVar->GetName(),&newValStr);
                        restore = false;
                    }
                } else {
                    restore = false;
                }
            }
        } else {
            long j = 0;
            for (j=0; j<nodeVRefs.lLength; j++) {
                thisVar = LocateVar (parameterVRefs.lData[k+j]);
                if (thisVar) {
                    _String newValStr2 = *LocateVar(nodeVRefs.lData[j])->GetName()&'.'&newValStr;
                    if (!thisVar->GetName()->Equal(&newValStr2)) {
                        if (LocateVarByName(newValStr2)>=0) {
                            newValStr2 = newValStr2 & " is already in use. Select another name.";
                            ProblemReport (newValStr2,(Ptr)this);
                            break;
                        }
                    } else {
                        restore = false;
                        break;
                    }
                }
            }
            if (j==nodeVRefs.lLength) {
                for (j=0; j<nodeVRefs.lLength; j++) {
                    _Variable* thisVar = LocateVar (parameterVRefs.lData[k+j]);
                    if (thisVar) {
                        newValStr = *LocateVar(nodeVRefs.lData[j])->GetName()&'.'&newValStr;
                        RenameVariable (thisVar->GetName(),&newValStr);
                        newValStr.Trim (newValStr.FindBackwards('.',0,-1)+1,-1);
                    }
                }
                restore = false;
            }
        }
    } else {
        newValStr = newValStr & " is an invalid identifier. Can't change parameter name to that!";
        ProblemReport (newValStr,(Ptr)this);
    }

    if (restore) {
        oldName.Trim (oldName.FindBackwards('.',0,-1)+1,-1);
        pl->SetCellData (&oldName,k,0,pl->cellTypes[row*pl->horizontalSpaces.lLength],true);
        pl->_MarkCellForUpdate (k);
    } else {
        *updateSCV = true;
    }
}

//__________________________________________________________

void    _HYTreePanel::ExecuteProcessor (long cID)
{
    if (cID<treeProcessors.lLength) {
        _FString treeStringV;
        *treeStringV.theString = GetTreeString();

        setParameter (treeStringExport, &treeStringV);

        FILE * thisFile = doFileOpen (((_String*)treeProcessors(cID))->sData,"rb");
        if (thisFile) {
            _String buffer (thisFile);
            fclose (thisFile);

            long   g = batchLanguageFunctionNames.lLength;

            _String tpp = *((_String*)treeProcessors(cID));
            PushFilePath (tpp);
            _ExecutionList   thisList;
            terminateExecution = false;
            thisList.BuildList (buffer);
            thisList.ExecuteAndClean(g);
            terminateExecution = false;

            PopFilePath ();
        }
    }
}

//__________________________________________________________

void ReadTreeProcessors (void)
{
    _String     pathToModelTemplates;
    _List       receptacle;

    pathToModelTemplates = libDirectory&"TreeAddIns";

    ScanDirectoryForFileNames (pathToModelTemplates,receptacle,false);

    for (long k=0; k<receptacle.lLength; k++) {
        FILE * thisFile = doFileOpen (((_String*)receptacle(k))->sData,"rb");
        if (thisFile) {
            _String buffer (thisFile);
            fclose (thisFile);
            if (buffer.sLength) {
                treeProcessors << receptacle(k);
            }
        }
    }
}

//EOF
