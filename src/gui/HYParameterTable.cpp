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

#include "HYDataPanel.h"
#include "HYTreePanel.h"
#include "HYParameterTable.h"
#include "HYTableComponent.h"
#include "HYLabel.h"
#include "HYButtonBar.h"
#include "HYUtils.h"
#include "HYEventTypes.h"
#include "parser.h"
#include "HYSharedMain.h"
#include "HYGWindow.h"
#include "HYChartWindow.h"

#include "math.h"

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

#ifdef  __MAC__
#define HY_PARAMETER_TABLE_ROW_HEIGHT 18
#else
#define HY_PARAMETER_TABLE_ROW_HEIGHT 20
#endif

extern  _SimpleList windowObjects;
extern  _HYColor    tableDefaultBk2,
        labelBColor,
        labelFColor;

extern  _String     blReplicate,
        donotWarnAgain,
        optimizationPrecision,
        useLastResults,
        none;

_SimpleList         categoryTableIcon,
                    globalTableIcon,
                    ghostGlobalTableIcon,
                    localTableIcon,
                    constrainedTableIcon,
                    constrainedGlobalTableIcon,
                    typeTableIcon,
                    treeTableIcon;

bool                autoUpdateTableLF = true,
                    doUndo            = true,
                    doTableUpdate   = true;

bool                warnAboutNonOptimized     = false,
                    warnAboutBootstrapReplace = false;

_String             allVariablesOption  ("All Parameters"),
                    undoEditValue       ("Undo Edit Value"),
                    undoSetConstraint   ("Undo Set Constarint"),
                    undoKillConstraint  ("Undo Kill Constraint"),
                    undoClearConstraints("Undo Clear Constraints"),
                    undoEqualConstraints("Undo Equal Constraint"),
                    undoMolecularClock  ("Undo Molecular Clock"),
                    undoProportionalSubtrees
                    ("Undo Proportional Constraint"),
                    undoCustomConstraints
                    ("Undo Custom Constraints"),
                    saveIteratesPrompt  ("Save Bootstrap Iterates to"),
                    bsSavePrompt        ("Save simulated data and MLEs to files"),
                    histChartBins       ("10");

extern              _String             covariancePrecision,
                    covarianceParameterList;

//__________________________________________________________
_HYParameterTable::_HYParameterTable (_String name, long lfIDRef):_HYTWindow (name, true)
{
    _HYRect         canvasSettings =
    {25,145,25,145,HY_COMPONENT_NO_SCROLL|HY_COMPONENT_BORDER_B|HY_COMPONENT_BORDER_REL};
    _HYButtonBar*   bb      = new _HYButtonBar (canvasSettings, GetOSWindowData());

    canvasSettings.left     = 70;
    canvasSettings.right    = 10000;

    _HYPullDown*    p1      = new _HYPullDown (canvasSettings,GetOSWindowData());

    canvasSettings.right    = 10000;
    canvasSettings.top      = 50;
    canvasSettings.left     = 50;
    canvasSettings.bottom   = 10000;
    canvasSettings.width    = HY_COMPONENT_H_SCROLL|HY_COMPONENT_V_SCROLL;

    _HYTable*   table       = new _HYTable (canvasSettings,GetOSWindowData(),0,4,100,HY_PARAMETER_TABLE_ROW_HEIGHT,HY_TABLE_STATIC_TEXT);

    canvasSettings.top      = canvasSettings.bottom     = HY_PARAMETER_TABLE_ROW_HEIGHT;
    canvasSettings.width    = HY_COMPONENT_NO_SCROLL | HY_COMPONENT_BORDER_B;

    _HYTable*   tableHead   = new _HYTable (canvasSettings,GetOSWindowData(),1,4,100,HY_PARAMETER_TABLE_ROW_HEIGHT,HY_TABLE_STATIC_TEXT);


    table->SetMessageRecipient (this);
    tableHead->SetMessageRecipient (this);
    bb->SetMessageRecipient (this);
    p1->SetMessageRecipient (this);

    AddObject (table);
    AddObject (bb);
    AddObject (tableHead);
    AddObject (p1);
    SetTableDimensions (3,2);
    SetCell   (HY_PARAMETER_TABLE_TABLE_ROW,0,table);
    SetCell   (HY_PARAMETER_TABLE_TABLE_ROW,1,table);
    SetCell   (HY_PARAMETER_TABLE_TABLE_ROW-1,0,tableHead);
    SetCell   (HY_PARAMETER_TABLE_TABLE_ROW-1,1,tableHead);
    SetCell   (HY_PARAMETER_TABLE_BUTTON_ROW,0,bb);
    SetCell   (HY_PARAMETER_TABLE_BUTTON_ROW,1,p1);

    _HYFont  labelFont;
    bb->SetBackColor (labelBColor);
    p1->SetBackColor (labelBColor);

#ifdef __MAC__
    labelFont.face  = "Geneva";
    labelFont.size  = 10;
#endif

#ifdef __WINDOZE__
    labelFont.face  = "MS Sans Serif";
    labelFont.size  = 12;
#endif

#ifdef __HYPHY_GTK__
    labelFont.face  = _HY_SANS_FONT;
    labelFont.size  = 12;
#endif


    labelFont.style = HY_FONT_PLAIN;

    table->SetFont(labelFont);
    tableHead->SetFont(labelFont);


    _String toolTipText ("Locate selected parameters in tree");
    bb->AddButton (ProcureIconResource(HY_PARAMETER_TABLE_BUTTON_ICON_BASE),&toolTipText);
    toolTipText = "Constrain selected parameters to be equal";
    bb->AddButton (ProcureIconResource(HY_PARAMETER_TABLE_BUTTON_ICON_BASE+1),&toolTipText);
    toolTipText = "Constrain 2 parameters to be proportional";
    bb->AddButton (ProcureIconResource(HY_PARAMETER_TABLE_BUTTON_ICON_BASE+2),&toolTipText);
    toolTipText = "Set bounds on parameter values";
    bb->AddButton (ProcureIconResource(HY_PARAMETER_TABLE_BUTTON_ICON_BASE+3),&toolTipText);
    toolTipText = "Impose molecular clock starting at the selected node";
    bb->AddButton (ProcureIconResource(HY_PARAMETER_TABLE_BUTTON_ICON_BASE+4),&toolTipText);
    toolTipText = "Relative ratio constraint";
    bb->AddButton (ProcureIconResource(HY_PARAMETER_TABLE_BUTTON_ICON_BASE+5),&toolTipText);
    toolTipText = "Clear constraints on selected variables";
    bb->AddButton (ProcureIconResource(HY_PARAMETER_TABLE_BUTTON_ICON_BASE+6),&toolTipText);

    bb->EnableButton(0,false);
    bb->EnableButton(1,false);
    bb->MarkAsPullDown (1,true);
    bb->EnableButton(2,false);
    bb->MarkAsPullDown (2,true);
    bb->EnableButton(3,false);
    bb->MarkAsPullDown (3,true);
    bb->EnableButton(4,false);
    bb->MarkAsPullDown (4,true);
    bb->EnableButton(5,false);
    bb->MarkAsPullDown (5,true);
    bb->EnableButton(6,false);
    bb->SetAlignFlags (HY_ALIGN_LEFT);

    bb->SetButtonDim(16);

    // set up the table

    lfID =lfIDRef;

    avViewOptions = viewOptions = GetViewOptions();
    table->selectionType = HY_TABLE_SEL_ROWS|HY_TABLE_FOCUSABLE;

    p1->SetAlignFlags (HY_ALIGN_LEFT);
    p1->AddMenuItem ("Current LF",-1);
    p1->AddMenuItem (menuSeparator,-1);
    p1->AddMenuItem ("Save LF state",-1);
    p1->AddMenuItem ("Delete LF state",-1);
    p1->AddMenuItem ("Select as null",-1);
    p1->AddMenuItem ("Select as alternative",-1);
    p1->AddMenuItem (menuSeparator,-1);
    p1->AddMenuItem ("LRT",-1);
    p1->AddMenuItem ("Parametric Bootstrap",-1);
    p1->AddMenuItem ("Nonparametric Bootstrap",-1);

    PrepareLFMenu   ();

    _String cellValue ("Parameter ID");
    tableHead->SetCellData (&cellValue,0,1,HY_TABLE_BOLD|HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_CANTSELECT|HY_TABLE_PULLDOWN,true);
    cellValue = "";
    tableHead->SetCellData (&cellValue,0,0,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_CANTSELECT|HY_TABLE_BOLD|HY_TABLE_PULLDOWN,true);
    cellValue = "Value";
    tableHead->SetCellData (&cellValue,0,2,HY_TABLE_BOLD|HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_CANTSELECT,true);
    cellValue = "Constraint";
    tableHead->SetCellData (&cellValue,0,3,HY_TABLE_BOLD|HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_CANTSELECT,true);

    ConstructTheTable ();

    _HYRect  screenRect = GetScreenDimensions();
    SetWindowRectangle (0,0,screenRect.bottom-80,screenRect.right-80);
    SetPosition (70,70);
    lastParameterCount = -1;

    lastLFState = "Current LF";

    SetWindowRectangle (0,0,screenRect.bottom-200,screenRect.right-20);
    SetPosition (5,40);
    DeleteObject (table);
    DeleteObject (bb);
    DeleteObject (tableHead);
    DeleteObject (p1);

}

//__________________________________________________________
_HYParameterTable::~_HYParameterTable()
{
}

//__________________________________________________________

bool    _HYParameterTable::ProcessGEvent (_HYEvent* e)
{
    _String firstArg,
            secondArg;
    long    k,f;

    bool    done = false;

    if (e->EventClass()==_hyGlobalLFKillEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        if (k!=GetID()) {
            firstArg        = e->EventCode().Cut (f+1,-1);
            if (lfID==firstArg.toNum()) {
                postWindowCloseEvent (GetID());
                done = true;
            }
        }
    } else if ((e->EventClass()==_hyGlobalChangeLF)||(e->EventClass()==_hyGlobalChangeLFParams)) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        if (k!=GetID()) {
            firstArg        = e->EventCode().Cut (f+1,-1);
            if (lfID==firstArg.toNum()) {
                //_HYTable*    theTable  = (_HYTable*)  GetCellObject    (HY_PARAMETER_TABLE_TABLE_ROW,0);

                if (e->EventClass()==_hyGlobalChangeLF) {
                    RefreshParameterValues();
                } else {
                    avViewOptions = viewOptions = GetViewOptions();
                    ConstructTheTable();
                    SetWindowRectangle (0,0,bottom,right);
                }

                done = true;
            }
        }
    }


    if (!done) {
        return _HYWindow::ProcessGEvent (e);
    }
    return true;
}

//__________________________________________________________

bool    _HYParameterTable::ProcessEvent (_HYEvent* e)
{
    _String firstArg,
            secondArg,
            thirdArg;
    long    k,i,f;

    bool    done = false;

    if (e->EventClass()==_hyTablePullDownEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }
        k               = e->EventCode().Find(',',f+1,-1);
        firstArg        = e->EventCode().Cut (f+1,k-1);
        secondArg       = e->EventCode().Cut (k+1,-1);
        if (i==0) {
            HandleCellPullDown (firstArg.toNum(),secondArg.toNum());
            done = true;
        } else if (i==2) {
            HandleHeaderPullDown (firstArg.toNum(),secondArg.toNum());
            done = true;
        }
    } else if (e->EventClass()==_hyTableEditCellEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }
        if (i==0) {
            firstArg        = e->EventCode().Cut (f+1,-1);
            TakeLFSnapshot();
            HandleCellEditEvent (firstArg.toNum());
            DoneLFSnapshot();
            done = true;
        }
    } else if (e->EventClass()==_hyTableChangeSelEvent) {
        k = e->EventCode().toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }
        if (i==0) {
            UpdateSelectionDependentButtons();
        }
    } else if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }
        if (i==1) {
            firstArg        = e->EventCode().Cut (f+1,-1);
            k               = firstArg.toNum();
            switch  (k) {
            case 0:
                DoOpenTreeWindow();
                break;

            case 1: // set equal
                TakeLFSnapshot();
                DoEqualConstraint ();
                DoneLFSnapshot();
                break;

            case 2: // set proportional
                TakeLFSnapshot();
                DoProportionalConstraint ();
                DoneLFSnapshot();
                break;

            case 3: // bounds
                TakeLFSnapshot();
                DoBoundsChange ();
                DoneLFSnapshot();
                break;

            case 4: // do mol clock
                TakeLFSnapshot();
                DoMolecularClock ();
                DoneLFSnapshot();
                break;

            case 5: // do clear constraints
                TakeLFSnapshot();
                DoProportionalSubtrees ();
                DoneLFSnapshot();
                break;

            case 6: // do clear constraints
                TakeLFSnapshot();
                DoClearConstraints ();
                DoneLFSnapshot();
                break;
            }
            done = true;
        }
    } else if (e->EventClass()==_hyTableDblClickEvent) {
        k = e->EventCode().toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }
        if (i==0) {
            DoOpenTreeWindow();
            done = true;
        }
    } else if (e->EventClass()==_hyScrollingEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        long g = e->EventCode().Find(',',f+1);
        if (g>=0) {
            k = firstArg.toNum();
            for (i=0; i<components.lLength; i++) {
                if (((_HYGuiObject*)components(i))->MatchID(k)) {
                    break;
                }
            }
            if (i==0) {
                firstArg = e->EventCode().Cut (f+1,g-1);
                k = firstArg.toNum();
                if (k) {
                    _HYTable*      theTable = (_HYTable*)     GetCellObject  (HY_PARAMETER_TABLE_TABLE_ROW-1,0);
                    theTable->ProcessEvent (generateScrollEvent (k,0));
                }
                done = true;
            }
        }
    } else if (e->EventClass()==_hyTableResizeCEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }
        if ((i==0)||(i==2)) {
            firstArg = e->EventCode().Cut (f+1,-1);
            f = firstArg.Find(',');
            k = firstArg.Cut(f+1,-1).toNum(); // shift
            f = firstArg.Cut(0,f-1).toNum();  // column
            _HYTable*      table;
            if (i==0) {
                table = (_HYTable*)     GetCellObject    (HY_PARAMETER_TABLE_TABLE_ROW-1,0);
            } else {
                table = (_HYTable*)     GetCellObject    (HY_PARAMETER_TABLE_TABLE_ROW,0);
            }
            table->SetColumnSpacing (f,k,true);
#ifdef __HYPHY_GTK__
            UpdateComponentInfo ();
#else
            dim = MinMaxWindowDimensions();
#endif
            done = true;
        }
    } else if (e->EventClass()==_hyMenuOpenEvent) {
        k = e->EventCode().toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }
        if (i==3) {
            PrepareLFMenu();
            done = true;
        }
    } else if (e->EventClass()==_hyMenuSelChangeEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }
        if (i==3) {
            k = e->EventCode().Cut (f+1,-1).toNum();

            _HYPullDown * p1 = (_HYPullDown*) GetCellObject (0,1);

            f = p1->MenuItemCount ();

            if (k == f-8) {
                HandleSaveLF();
            } else if (k&&(k<f-9)) {
                HandleRestoreLF (p1->GetMenuItem (k));
            } else if ((k == f-6)||(k==f-5)) {
                HandleSetHypothesis (&lastLFState, k==f-5);
            } else if (k == f-7) {
                HandleDeleteLF(&lastLFState);
            } else if (k == f-3) {
                HandleLRT ();
            } else if (k >= f-2) {
                HandleBootstrap (k-f+2);
            }

            lastLFState = *p1->GetMenuItem(p1->GetSelection());
        }
    }
    if (done) {
        DeleteObject (e);
        return true;
    }
    return _HYTWindow::ProcessEvent(e);
}

//__________________________________________________________

void    _HYParameterTable::SetSelectedRows (_List& theNames)
{
    _HYTable*      theTable = (_HYTable*)     GetCellObject  (HY_PARAMETER_TABLE_TABLE_ROW,0);
    long           k,f;
    _SimpleList    rowsToSelect,
                   patternIndex;
    _List          patterns;

    for (k=1; k<theTable->verticalSpaces.lLength; k++) {
        _String currentRow = *(_String*)theTable->GetCellData(1,k);
        f = currentRow.FindBackwards ('.',0,-1);
        if (f>=0) {
            currentRow.Trim (0,f);
            patterns && & currentRow;
            patternIndex << k;
        }
    }
    for (k=0; k<patterns.lLength; k++) {
        if (theNames.BinaryFind(patterns(k))>=0) {
            rowsToSelect << patternIndex.lData[k];
        }
    }
    theTable->SetRowSelection (rowsToSelect);
    if (rowsToSelect.lLength) {
        theTable->_ScrollRowIntoView(rowsToSelect.lData[0]);
    }
}

//__________________________________________________________

void    _HYParameterTable::SelectAll (void)
{
    _HYTable*      theTable = (_HYTable*)     GetCellObject  (HY_PARAMETER_TABLE_TABLE_ROW,0);
    long           k;
    _SimpleList    rowsToSelect;

    for (k=1; k<theTable->verticalSpaces.lLength; k++) {
        rowsToSelect << k;
    }
    theTable->SetRowSelection (rowsToSelect);
}
//__________________________________________________________

long    _HYParameterTable::GatherListOfGlobals (_List& glList)
{
    _HYTable*      theTable = (_HYTable*)     GetCellObject  (HY_PARAMETER_TABLE_TABLE_ROW,0);
    long           k,
                   res = -1;

    for (k=0; k<theTable->verticalSpaces.lLength; k++) {
        _Variable* theV = (_Variable*)RetrieveIthRowVar (k);
        if (theV&&(!theV->IsCategory())&&(theV->ObjectClass()!=TREE)) {
            if (theV->IsGlobal()) {
                if (theV->IsIndependent()) {
                    res = k;
                }

                glList << theV->GetName();
            }
        }
    }

    return res+1;
}

//__________________________________________________________

void    _HYParameterTable::SetRowValues (long rowCount, _Variable* tv, _HYTable* theTable)
{
    if (tv->ObjectClass()==TREE) {
        theTable->SetCellData (&treeTableIcon,rowCount,0,HY_TABLE_ICON,true);
        theTable->SetCellData (tv->GetName(),rowCount,1,HY_TABLE_STATIC_TEXT,true);
        theTable->SetCellData (&empty,rowCount,2,HY_TABLE_STATIC_TEXT,true);
        theTable->SetCellData (&empty,rowCount,3,HY_TABLE_STATIC_TEXT,true);
    } else {
        if (tv->IsCategory()) {
            theTable->SetCellData (&categoryTableIcon,rowCount,0,HY_TABLE_ICON,true);
            theTable->SetCellData (tv->GetName(),rowCount,1,HY_TABLE_STATIC_TEXT,true);
            _String   cellValue = _String("1(")&_String (((_CategoryVariable*)tv)->GetNumberOfIntervals())&')';
            theTable->SetCellData (&cellValue,rowCount,2,HY_TABLE_STATIC_TEXT,true);
            ((_CategoryVariable*)tv)->Refresh();
            cellValue = _String (((_CategoryVariable*)tv)->GetIntervalValue(0));
            theTable->SetCellData (&cellValue,rowCount,3,HY_TABLE_STATIC_TEXT|HY_TABLE_PULLDOWN,true);
        } else if (tv->IsIndependent()) {
            if (tv->IsGlobal()) {
                theTable->SetCellData (&globalTableIcon,rowCount,0,HY_TABLE_ICON,true);
            } else {
                theTable->SetCellData (&localTableIcon,rowCount,0,HY_TABLE_ICON,true);
            }
            theTable->SetCellData (tv->GetName(),rowCount,1,HY_TABLE_STATIC_TEXT,true);
            _String cellValue (tv->Value());
            theTable->SetCellData (&cellValue,rowCount,2,HY_TABLE_EDIT_TEXT,true);
            theTable->SetCellData (&empty,rowCount,3,HY_TABLE_EDIT_TEXT,true);
        } else {
            if (tv->IsGlobal()) {
                theTable->SetCellData (&constrainedGlobalTableIcon,rowCount,0,HY_TABLE_ICON,true);
            } else {
                theTable->SetCellData (&constrainedTableIcon,rowCount,0,HY_TABLE_ICON,true);
            }
            theTable->SetCellData (tv->GetName(),rowCount,1,HY_TABLE_STATIC_TEXT,true);
            _String cellValue (tv->Compute()->Value());
            theTable->SetCellData (&cellValue,rowCount,2,HY_TABLE_STATIC_TEXT,true);
            _String*  theF = tv->GetFormulaString();
            theTable->SetCellData (theF,rowCount,3,HY_TABLE_EDIT_TEXT,false);
        }
    }
}

//__________________________________________________________

void    _HYParameterTable::ConstructTheTable (void)
{
    _HYTable*      theTable  = (_HYTable*)  GetCellObject    (HY_PARAMETER_TABLE_TABLE_ROW,0),
                   *       hTable    = (_HYTable*)  GetCellObject    (HY_PARAMETER_TABLE_TABLE_ROW-1,0);

    long           k,
                   rowCount = 0;

    if (globalTableIcon.lLength==0) {
        globalTableIcon<<(long)ProcureIconResource (HY_PARAMETER_TABLE_ICON_BASE);
        globalTableIcon<<16;
        globalTableIcon<<16;
        localTableIcon<<(long)ProcureIconResource (HY_PARAMETER_TABLE_ICON_BASE+1);
        localTableIcon<<16;
        localTableIcon<<16;
        constrainedTableIcon<<(long)ProcureIconResource (HY_PARAMETER_TABLE_ICON_BASE+2);
        constrainedTableIcon<<16;
        constrainedTableIcon<<16;
        categoryTableIcon<<(long)ProcureIconResource (HY_PARAMETER_TABLE_ICON_BASE+3);
        categoryTableIcon<<16;
        categoryTableIcon<<16;
        treeTableIcon<<(long)ProcureIconResource (HY_PARAMETER_TABLE_ICON_BASE+4);
        treeTableIcon<<16;
        treeTableIcon<<16;
        ghostGlobalTableIcon<<(long)ProcureIconResource (HY_PARAMETER_TABLE_ICON_BASE+5);
        ghostGlobalTableIcon<<16;
        ghostGlobalTableIcon<<16;
        constrainedGlobalTableIcon<<(long)ProcureIconResource (HY_PARAMETER_TABLE_ICON_BASE+6);
        constrainedGlobalTableIcon<<16;
        constrainedGlobalTableIcon<<16;
    }

    theTable->ClearTable();

    //theTable->_SetVScrollerPos (0);
    //theTable->_SetHScrollerPos (0);

    if (lfID>=0) {
        _LikelihoodFunction* lf = (_LikelihoodFunction*)likeFuncList (lfID);
        _Variable*           tv;

        if (viewOptions&HY_PARAMETER_TABLE_VIEW_TREES) {
            k = lf->GetTheTrees().lLength;
            for (long j=0; j<k; j++) {
                tv = LocateVar (lf->GetTheTrees().lData[j]);
                theTable->AddRow (-1,HY_PARAMETER_TABLE_ROW_HEIGHT,HY_TABLE_STATIC_TEXT);
                SetRowValues (rowCount, tv, theTable);
                rowCount++;
            }
        }

        _List       varList1, varList2;

        for (k=0; k<lf->GetIndependentVars().lLength; k++) {
            tv = LocateVar (lf->GetIndependentVars().lData[k]);
            if (tv->IsGlobal()) {
                if (viewOptions&HY_PARAMETER_TABLE_VIEW_GLOBAL) {
                    varList1.BinaryInsert (tv->GetName());
                } else {
                    continue;
                }
            } else {
                if (viewOptions&HY_PARAMETER_TABLE_VIEW_LOCAL) {
                    varList2.BinaryInsert (tv->GetName());
                } else {
                    continue;
                }
            }
        }

        for (k=0; k<varList1.lLength; k++) {
            theTable->AddRow (-1,HY_PARAMETER_TABLE_ROW_HEIGHT,HY_TABLE_STATIC_TEXT);
            tv = FetchVar (LocateVarByName (*(_String*)varList1(k)));
            SetRowValues (rowCount, tv, theTable);
            rowCount++;
        }

        for (k=0; k<varList2.lLength; k++) {
            theTable->AddRow (-1,HY_PARAMETER_TABLE_ROW_HEIGHT,HY_TABLE_STATIC_TEXT);
            tv = FetchVar (LocateVarByName (*(_String*)varList2(k)));
            SetRowValues (rowCount, tv, theTable);
            rowCount++;
        }

        varList1.Clear();

        if (viewOptions&HY_PARAMETER_TABLE_VIEW_CONSTRAINED) {
            for (k=0 ; k<lf->GetDependentVars().lLength; k++) {
                tv = LocateVar (lf->GetDependentVars().lData[k]);
                varList1.BinaryInsert (tv->GetName());
            }
            for (k=0; k<varList1.lLength; k++) {
                theTable->AddRow (-1,HY_PARAMETER_TABLE_ROW_HEIGHT,HY_TABLE_STATIC_TEXT);
                tv = FetchVar (LocateVarByName (*(_String*)varList1(k)));
                SetRowValues (rowCount, tv, theTable);
                rowCount++;
            }
            varList1.Clear();
        }
        if (viewOptions&HY_PARAMETER_TABLE_VIEW_CATEGORY) {
            for (k=0 ; k<lf->GetCategoryVars().lLength; k++) {
                tv = LocateVar (lf->GetCategoryVars().lData[k]);
                varList1.BinaryInsert (tv->GetName());
            }
            for (k=0; k<varList1.lLength; k++) {
                theTable->AddRow (-1,HY_PARAMETER_TABLE_ROW_HEIGHT,HY_TABLE_STATIC_TEXT);
                tv = FetchVar (LocateVarByName (*(_String*)varList1(k)));
                SetRowValues (rowCount, tv, theTable);
                rowCount++;
            }
        }
        theTable->AutoFitWidth(*hTable);
        UpdateLogLikelihood ();
    }
}


//__________________________________________________________

void    _HYParameterTable::RefreshTheTable (void)
{
    if (lfID>=0) {
        _HYTable*      table    = (_HYTable*) GetCellObject  (HY_PARAMETER_TABLE_TABLE_ROW,0);
        long idx = 0;
        bool update;

        for (long k=0; k<table->verticalSpaces.lLength; k++,idx+=table->horizontalSpaces.lLength) {
            _Variable* thisV = (_Variable*)RetrieveIthRowVar(k);
            update = false;
            if (thisV->ObjectClass()!=TREE) {
                if (!thisV->IsCategory()) {
                    if (thisV->IsIndependent()) {
                        if (thisV->IsGlobal()) {
                            if ((globalTableIcon.lData[0]!=(long)((_SimpleList*)table->GetCellData(0,k))->lData[0])
                                    &&(ghostGlobalTableIcon.lData[0]!=(long)((_SimpleList*)table->GetCellData(0,k))->lData[0])) {
                                table->SetCellData(&globalTableIcon,k,0,table->cellTypes[idx],true);
                                table->SetCellData (&empty,k,3,table->cellTypes[idx+3],true);
                            }
                        } else {
                            if (localTableIcon.lData[0]!=(long)((_SimpleList*)table->GetCellData(0,k))->lData[0]) {
                                table->SetCellData(&localTableIcon,k,0,table->cellTypes[idx],true);
                                table->SetCellData (&empty,k,3,table->cellTypes[idx+3],true);
                            }
                        }
                    } else {
                        _SimpleList * targetIcon;

                        if (thisV->IsGlobal()) {
                            targetIcon = &constrainedGlobalTableIcon;
                        } else {
                            targetIcon = &constrainedTableIcon;
                        }

                        //if ((targetIcon->lData[0]!=(long)((_SimpleList*)table->GetCellData(0,k))->lData[0])
                        //&&(ghostGlobalTableIcon.lData[0]!=(long)((_SimpleList*)table->GetCellData(0,k))->lData[0]))
                        {
                            table->SetCellData(targetIcon,k,0,table->cellTypes[idx],true);
                            _String* fla = thisV->GetFormulaString();
                            table->SetCellData(fla,k,3,table->cellTypes[idx+3],false);
                            if (thisV->varFormula && thisV->varFormula->IsAConstant ()) {
                                update = true;
                            }
                        }
                    }
                    if (update||(thisV->HasChanged())) {
                        _String newVal (thisV->Compute()->Value());
                        table->SetCellData(&newVal,k,2,table->cellTypes[idx+2],true);
                    }
                } else {
                    UpdateKthRow (k,false);
                }
            }
        }
        UpdateLogLikelihood ();
    }
}

//__________________________________________________________

void    _HYParameterTable::RefreshParameterValues (void)
{
    if (lfID>=0) {
        _HYTable*      table    = (_HYTable*) GetCellObject  (HY_PARAMETER_TABLE_TABLE_ROW,0);
        for (long k=0, idx = 0; k<table->verticalSpaces.lLength; k++,idx+=table->horizontalSpaces.lLength) {
            _Variable* thisV = (_Variable*)RetrieveIthRowVar(k);
            if (thisV->ObjectClass()!=TREE) {
                if (!thisV->IsCategory()) {
                    _String newVal (thisV->Compute()->Value());
                    table->SetCellData(&newVal,k,2,table->cellTypes[idx+2],true);
                }
            }
        }
        table->_MarkColumnForUpdate (2);
        UpdateLogLikelihood ();
    }
}

//__________________________________________________________

void    _HYParameterTable::VerifyGlobalVariables (void)
{
    if (lfID>=0) {
        _HYTable*      table    = (_HYTable*) GetCellObject  (HY_PARAMETER_TABLE_TABLE_ROW,0);
        long idx = 0,
             lastGlobal = table->verticalSpaces.lLength;

        bool update;

        _LikelihoodFunction* lf = (_LikelihoodFunction*)likeFuncList (lfID);

        lf->RescanAllVariables ();

        _SimpleList gv,
                    tc,
                    found;

        lf->GetGlobalVars(gv);
        gv.Sort();

        for (long k=0; k<table->verticalSpaces.lLength; k++,idx+=table->horizontalSpaces.lLength) {
            _Variable* thisV = (_Variable*)RetrieveIthRowVar(k);
            update = false;
            long     selected = table->cellTypes[idx]&HY_TABLE_SELECTED;

            if ((thisV->ObjectClass()!=TREE)&&(!thisV->IsCategory())&&(thisV->IsGlobal())) {
                long f = gv.BinaryFind (thisV->GetAVariable());
                if (f<0) {
                    if (ghostGlobalTableIcon.lData[0]!=(long)((_SimpleList*)table->GetCellData(0,k))->lData[0]) {
                        table->SetCellData (&ghostGlobalTableIcon,k,0,table->cellTypes[idx],true);
                        table->SetCellData (&empty,k,2,HY_TABLE_STATIC_TEXT|HY_TABLE_ITALIC|selected,true);
                        table->SetCellData (&empty,k,3,HY_TABLE_STATIC_TEXT|selected,true);
                        table->cellTypes[idx+1] = HY_TABLE_STATIC_TEXT|HY_TABLE_ITALIC|selected;
                        tc << idx;
                        tc << idx+1;
                        tc << idx+2;
                        tc << idx+3;
                        UpdateKthRow (k);
                    }
                } else {
                    found << gv.lData[f];
                    if (thisV->IsIndependent())
                        if (globalTableIcon.lData[0]!=(long)((_SimpleList*)table->GetCellData(0,k))->lData[0]) {
                            table->SetCellData (&globalTableIcon,k,0,table->cellTypes[idx],true);
                            table->SetCellData (&empty,k,2,HY_TABLE_EDIT_TEXT|selected,true);
                            table->SetCellData (&empty,k,3,HY_TABLE_EDIT_TEXT|selected,true);
                            table->cellTypes[idx+1] = HY_TABLE_STATIC_TEXT|selected;
                            tc << idx;
                            tc << idx+1;
                            tc << idx+2;
                            tc << idx+3;
                            UpdateKthRow (k);
                        }
                    table->cellTypes[idx+3] = HY_TABLE_EDIT_TEXT|selected;
                    lastGlobal = k;
                }
            } else {
                if (thisV->ObjectClass()==TREE) {
                    lastGlobal = k;
                }
            }
        }
        if (found.lLength<gv.lLength) {
            found.Sort();
            _SimpleList diff;
            diff.Subtract (gv,found);
            for (idx=0; idx<diff.lLength; idx++,lastGlobal++) {
                table->AddRow (lastGlobal,HY_PARAMETER_TABLE_ROW_HEIGHT,HY_TABLE_STATIC_TEXT);
                SetRowValues  (lastGlobal,LocateVar (diff.lData[idx]),table);
            }
            table->_MarkForUpdate();
            tc.Clear();
        }

        if (tc.lLength) {
            table->_MarkCellsForUpdate (tc);
        }
    }
}

//__________________________________________________________

void    _HYParameterTable::UpdateLogLikelihood (bool postEvent)
{
    if (lfID>=0) {
        _String newStatLine;
        if (lfID == lockedLFID) {
            newStatLine = _String ("This likelihood function is currently being optimized...");
            SetStatusBar (newStatLine);
            return;
        }

        _LikelihoodFunction* lf = (_LikelihoodFunction*)likeFuncList (lfID);
        lf->PrepareToCompute ();
        _Parameter lfval    = lf->Compute();
        long       parCount = lf->GetIndependentVars().lLength;
        lf->DoneComputing ();

        newStatLine = _String ("Log Likelihood = ")&_String (lfval) &", parameter count = "&
                      _String(parCount) & ", AIC = " & _String(2.*(parCount-lfval)) & '.' ;

        SetStatusBar (newStatLine);
        if (postEvent) {
            postChangeLFEvent (GetID(),lfID);
        }

        _HYPullDown* p1 = (_HYPullDown*)GetCellObject (HY_PARAMETER_TABLE_BUTTON_ROW,1);

        if (p1->GetSelection()) {
            p1->EnableItem (0,true);
            p1->ChangeSelection (0, false);
        }
    }
}

//__________________________________________________________

void    _HYParameterTable::TaintTheLF ()
{
    if (lfID>=0) {
        _LikelihoodFunction* lf = (_LikelihoodFunction*)likeFuncList (lfID);
        for (long k=lf->GetIndependentVars().lLength-1; k>=0; k--) {
            lf->SetIthIndependent (k,lf->GetIthIndependent(k));
        }
    }
}

//__________________________________________________________

void    _HYParameterTable::TakeLFSnapshot ()
{
    /*if (lfID>=0)
    {
        if (!(viewOptions & HY_PARAMETER_TABLE_SNAPSHOT))
        {
            _LikelihoodFunction* lf = (_LikelihoodFunction*)likeFuncList (lfID);
            lastParameterCount = lf->GetIndependentVars().lLength;
            lastLFValue        = lf->Compute();
            viewOptions       |= HY_PARAMETER_TABLE_SNAPSHOT;
        }
    }*/
}

//__________________________________________________________

void    _HYParameterTable::PrepareLFMenu (void)
{
    _HYPullDown* p1 = (_HYPullDown*)GetCellObject (0,1);
    for (; !menuSeparator.Equal (p1->GetMenuItem(1));) {
        p1->DeleteMenuItem (1);
    }

    long k = 0;

    _HYDataPanel*parent = (_HYDataPanel*)RetrieveParentDataPanel ();
    long n = -1,
         a = -1,
         s = p1->GetSelection ();
    if (parent) {
        for (k=0; k<parent->savedLFStates.lLength; k++) {
            p1->AddMenuItem (*(_String*)parent->savedLFNames (k),k+1);
        }
        n = parent->GetHypothesis (false)+1;
        a = parent->GetHypothesis (true)+1;
        s = p1->GetSelection ();

    }

    k = p1->MenuItemCount();

    bool flag = false;

    if ((n>=1)&&(a>=1)) {
        flag = true;
    }

    p1->EnableItem (k-1,flag);
    p1->EnableItem (k-2,flag);
    p1->EnableItem (k-3,flag);

    flag = (s>0) && (n!=s);
    p1->EnableItem (k-6,flag);

    flag = (s>0) && (a!=s);
    p1->EnableItem (k-5,flag);

    p1->EnableItem (k-7,s>0);
}

//__________________________________________________________

void    _HYParameterTable::HandleSaveLF (void)
{
    _HYDataPanel*parent = (_HYDataPanel*)RetrieveParentDataPanel ();
    _HYPullDown* p1 = (_HYPullDown*)GetCellObject (0,1);
    if (parent) {
        _LikelihoodFunction* lf = (_LikelihoodFunction*)likeFuncList (lfID);
        bool    opt = lf->HasBeenOptimized();
        if (!(opt||warnAboutNonOptimized)) {
            _String warnNO ("The parameters values you are about to save are NOT maximum likelihood estimates. LRT and bootstraps may be meaningless.");
            if (!ProceedPromptWithCheck (warnNO,donotWarnAgain,warnAboutNonOptimized, (Ptr)this)) {
                p1->ChangeSelection (p1->FindMenuItem (lastLFState),false);
                return;
            }
        }
        _String newLF ("New LF Snapshot"),
                prompt("Type in a name for the saved LF state:");
        if (EnterStringDialog (newLF, prompt, (Ptr)this)) {
            if (!opt) {
                newLF = _String("[-]")&newLF;
            }
            _HYPullDown* p1 = (_HYPullDown*) GetCellObject (0,1);
            long k = parent->FindLFState (newLF);

            while (k>=0) {
                prompt = "A likelihood function state with that name already exists. Would you like to replace it?";
                if (!ProceedPrompt (prompt,(Ptr)this)) {
                    p1->ChangeSelection (p1->FindMenuItem (lastLFState),false);
                    return;
                }
                /*prompt = "Type in a name for the saved LF state:";
                if (!EnterStringDialog (newLF, prompt, (Ptr)this))
                {
                    p1->ChangeSelection (p1->FindMenuItem (lastLFState),false);
                    return;
                }*/
                if (!opt) {
                    newLF = _String("[-]")&newLF;
                }
            }
            _String * snapshot = parent->LFSnapshot();
            if (snapshot) {
                if (k<0) {
                    parent->savedLFNames.InsertElement (&newLF,k,true);
                    parent->savedLFStates.InsertElement (snapshot,k,false);
                    k = parent->savedLFNames.lLength;
                } else {
                    parent->savedLFStates.Replace (k,snapshot,false);
                    k++;
                }
                PrepareLFMenu ();
                p1->EnableItem (0,false);
                p1->ChangeSelection (k,false);
            } else {
                p1->ChangeSelection (p1->FindMenuItem (lastLFState),false);
            }

        } else {
            p1->ChangeSelection (p1->FindMenuItem (lastLFState),false);
        }

    } else {
        p1->ChangeSelection(0,false);
        _String warnNODS ("Sorry, but you may only save likelihood function states if the likelihood function was created using the data panel interface.");
        ProblemReport (warnNODS);
    }
}

//__________________________________________________________

void    _HYParameterTable::HandleLRT (void)
{
    _HYPullDown*    p1  = (_HYPullDown*) GetCellObject   (HY_PARAMETER_TABLE_BUTTON_ROW,1);
    //p1->EnableItem  (0,false);
    p1->ChangeSelection (p1->FindMenuItem (lastLFState),false);

    _HYDataPanel*parent = (_HYDataPanel*)RetrieveParentDataPanel ();
    if (parent) {
        _String * currentSnapshot = parent->LFSnapshot();

        _Parameter  H_0,
                    H_A,
                    P = 1.0;

        long        DF_0,
                    DF_A;

        _LikelihoodFunction* lf = (_LikelihoodFunction*)likeFuncList (lfID);
        if (!parent->LFRestore(parent->GetHypothesis(false))) {
            return;
        }

        lf->RescanAllVariables();

        lf->PrepareToCompute ();
        H_0 = lf->Compute();
        DF_0= lf->GetIndependentVars().lLength;
        lf->DoneComputing ();

        if (!parent->LFRestore(parent->GetHypothesis(true))) {
            return;
        }

        lf->RescanAllVariables();

        lf->PrepareToCompute ();
        H_A = lf->Compute();
        DF_A= lf->GetIndependentVars().lLength;
        lf->DoneComputing ();

        _ExecutionList exl (*currentSnapshot);
        exl.Execute();

        if (DF_A < DF_0) {
            _String errMsg ("The null model should contain fewer parameters than the alternative model. Unable to obtain the asymptotic P-Value.");
            ProblemReport (errMsg,(Ptr)this);
            P = -1.;
        } else if (H_A < H_0) {
            _String errMsg ("The null hypothesis should have worse likelihood than that of the alternative hypothesis. Unable to obtain the asymptotic P-Value.");
            ProblemReport (errMsg,(Ptr)this);
            P = -1.;
        }


        char            cb[128];
        snprintf (cb, sizeof(cb),"\nLikelihood Ratio Test\n\n\t 2*LR = %g\n\t DF = %ld\n\t P-Value = ",2.*(H_A-H_0),DF_A-DF_0);
        if (P>0.0) {
            _Constant       c1 (2.*(H_A-H_0)),
                            c2 (DF_A-DF_0), *c3 = (_Constant*)c1.CChi2(&c2);

            BufferToConsole (cb);
            snprintf (cb, sizeof(cb),"%g \n ",1.-c3->Value());
            BufferToConsole (cb);
            DeleteObject    (c3);
        } else {
            BufferToConsole     ("Unable to compute \n ");
        }


        DeleteObject (currentSnapshot);
    }
}

//__________________________________________________________

void    _HYParameterTable::HandleBootstrap (bool np)
{
    _HYPullDown*    p1  = (_HYPullDown*) GetCellObject   (HY_PARAMETER_TABLE_BUTTON_ROW,1);
    p1->ChangeSelection (p1->FindMenuItem (lastLFState),false);

    long            iters = -1;
    bool            saveIterates = true;

    _HYDataPanel*parent = (_HYDataPanel*)RetrieveParentDataPanel ();
    if (parent) {
        _String  bsWindow   ("Bootstrap Window for "),
                 folderSave;

        bsWindow = bsWindow & *(_String*)likeFuncNamesList (lfID);

        _HYBootstrapWindow * bsw = nil;

        long              windID = FindWindowByName (bsWindow);

        if (windID>=0) {
            bsw = (_HYBootstrapWindow*)windowObjectRefs (windID);

            _String warnMessage;
            warnMessage = _String ("The Window '") & bsWindow & "' is already open. Should I append more iterates to the existing window? (Click 'No' to start a fresh bootstrap).";
            char res = YesNoCancelPrompt (warnMessage);

            if (res==2) {
                return;
            }

            if (res==3) {
                bsw->Close(nil);
                bsw = nil;
            }
        }

        _String  promptString ("Enter the number of iterations:"),
                 valueString ("100");
        while (1) {
            if (!EnterStringDialogWithCheckbox (valueString, promptString, bsSavePrompt, saveIterates, (Ptr)this)) {
                return;
            }
            iters = valueString.toNum();
            if (iters>0) {
                break;
            } else {
                valueString = "100";
            }
        }

        if (saveIterates) {
            folderSave = ChooseAFolder (saveIteratesPrompt);
            if (folderSave.sLength==0) {
                return;
            }
        }

        parent->SetLockState (true);
        _String  * currentSnapshot = parent->LFSnapshot();
        checkPointer (currentSnapshot);

        // compute LRT
        _Parameter  nullLRT;
        //long      simP = 0;

        _LikelihoodFunction* lf = (_LikelihoodFunction*)likeFuncList (lfID);
        parent->LFRestore(parent->GetHypothesis(false));
        lf->RescanAllVariables();
        lf->PrepareToCompute ();
        nullLRT = -lf->Compute();
        lf->DoneComputing ();
        parent->LFRestore(parent->GetHypothesis(true));
        lf->RescanAllVariables();
        lf->PrepareToCompute ();
        nullLRT += lf->Compute();
        nullLRT *= 2;
        lf->DoneComputing ();

        _DataSet*originalData = (_DataSet*)dataSetList (parent->dataSetID);
        _String *originalName = (_String*)dataSetNamesList (parent->dataSetID);

        _String  dataSetName (*originalName);

        dataSetName = dataSetName & "_sim";
        FindUnusedObjectName (empty,dataSetName,dataSetNamesList);

        _DataSet *ds = nil;
        long     sID = -1;
        if (!np) {
            ds = new _DataSet();
            sID = AddDataSetToList (dataSetName,ds);
        }
        StartBarTimer ();

        //(_HYBootstrapWindow*)FindWindowByNameAndOpen (bsWindow);


        if (!bsw) {
            bsw = new _HYBootstrapWindow (bsWindow,nullLRT);
            checkPointer (bsw);
            bsw->BringToFront();
        } else {
            bsw->SetStopped (false);
        }

        for (long k = 0; k < iters; k++) {
            _String    dataFileName = folderSave& "simulatedData" & (k+1) & ".seq";

            _Parameter simLRT;

            _List       dummy;
            _SimpleList dummy2;

            parent->LFRestore(parent->GetHypothesis(false));
            lf->RescanAllVariables();
            if (np) {
                parent->SpawnLikelihoodFunctionNP (dummy);
            } else {
                ds->Clear();
                lf->Simulate (*ds,dummy);
                parent->SpawnLikelihoodFunction (ds,&dataSetName, dummy, dummy2);
            }


            valueString = _String("Obtaining MLEs for the null hypothesis. Iterate #")& k;
            SetStatusLine (valueString);

            _Matrix* res = lf->Optimize();
            simLRT = -(*res)(1,0);
            DeleteObject (res);

            if (saveIterates) {
                _String     saveFile     = folderSave & "nullMLEs" & (k+1) & ".bf";
                if (np) {
                    _DataSet* dsNP = parent->GenerateOrderedDataSet ();
                    parent->SaveDataPanel (false, &saveFile, &dataFileName, false, dsNP, true);
                    DeleteObject (dsNP);
                } else {
                    parent->SaveDataPanel (false, &saveFile, &dataFileName, false, ds);
                }

            }

            parent->LFRestore(parent->GetHypothesis(true));
            lf->RescanAllVariables();

            valueString = _String("Obtaining MLEs for the alternative hypothesis. Iterate #")& k;
            SetStatusLine (valueString);
            res = lf->Optimize();
            simLRT += (*res)(1,0);
            DeleteObject (res);

            if (saveIterates) {
                _String     saveFile     = folderSave & "altMLEs" & (k+1) & ".bf";
                parent->SaveDataPanel (false, &saveFile, &dataFileName, false, nil, np);
            }

            simLRT *= 2.;

            if (np) {
                parent->SpawnLikelihoodFunctionNP (dummy);
            } else {
                parent->SpawnLikelihoodFunction (originalData,originalName, dummy, dummy2);
            }

            bsw->AddIterate (simLRT);

            if (bsw->requestQuit) {
                break;
            }
        }


        if (!np) {
            KillDataSetRecord (sID);
            //DeleteObject (ds);
        }

        _ExecutionList exl (*currentSnapshot);
        exl.Execute();
        DeleteObject (currentSnapshot);
        parent->SetLockState (false);
        bsw->SetStopped (true);
        StopBarTimer();
        bsw->requestQuit = false;
    }
}

//__________________________________________________________

void    _HYParameterTable::HandleRestoreLF (_String* lfName)
{
    _HYDataPanel*parent = (_HYDataPanel*)RetrieveParentDataPanel ();
    if (parent) {
        long k = parent->savedLFNames.Find (lfName);
        if (k>=0) {
            if (!parent->LFRestore (k)) {
                return;
            }

            VerifyGlobalVariables();
            RefreshTheTable ();


            _HYTable*    table = (_HYTable*) GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);
            StretchColumnToFit (2);
            StretchColumnToFit (3);
            table->_MarkColumnForUpdate (0);

            undoCommands.Clear();
            undoDescriptions.Clear();
            _UpdateUndoMenu (nil,nil);

            PrepareLFMenu ();
            _HYPullDown*    p1  = (_HYPullDown*) GetCellObject   (HY_PARAMETER_TABLE_BUTTON_ROW,1);
            p1->EnableItem  (0,false);
            p1->ChangeSelection (k+1,false);

        }
    }
}

//__________________________________________________________

void    _HYParameterTable::HandleDeleteLF (_String* lfName)
{
    _HYDataPanel*parent = (_HYDataPanel*)RetrieveParentDataPanel ();
    if (parent) {
        long k = parent->savedLFNames.Find (lfName);
        if (k>=0) {
            parent->savedLFNames.Delete (k);
            parent->savedLFStates.Delete (k);
            parent->tainted = true;

            VerifyGlobalVariables();
            RefreshTheTable ();

            PrepareLFMenu ();

            _HYPullDown*    p1  = (_HYPullDown*) GetCellObject   (HY_PARAMETER_TABLE_BUTTON_ROW,1);
            p1->EnableItem  (0,true);
            p1->ChangeSelection (0,false);
        }
    }
}

//__________________________________________________________

void    _HYParameterTable::HandleSetHypothesis (_String* sel, bool alt)
{
    _HYDataPanel*parent = (_HYDataPanel*)RetrieveParentDataPanel ();
    if (parent) {
        parent->SetHypothesis (sel,alt);
        PrepareLFMenu ();
        _HYPullDown*    p1  = (_HYPullDown*) GetCellObject   (HY_PARAMETER_TABLE_BUTTON_ROW,1);
        p1->ChangeSelection (parent->GetHypothesis (alt)+1);
    }
}


//__________________________________________________________

char    _HYParameterTable::GetViewOptions (void)
{
    char res  = 0;
    if (lfID>=0) {
        _LikelihoodFunction* lf = (_LikelihoodFunction*)likeFuncList (lfID);
        res = HY_PARAMETER_TABLE_VIEW_LOCAL|HY_PARAMETER_TABLE_VIEW_TREES;
        if (lf->CountObjects (1)) {
            res |= HY_PARAMETER_TABLE_VIEW_GLOBAL;
        }
        if (lf->CountObjects (3)) {
            res |= HY_PARAMETER_TABLE_VIEW_CONSTRAINED;
        }
        if (lf->GetCategoryVars().lLength) {
            res |= HY_PARAMETER_TABLE_VIEW_CATEGORY;
        }
    }
    return res;
}

//__________________________________________________________

BaseRef _HYParameterTable::RetrieveIthRowVar (long index)
{
    BaseRef     result = nil;
    _HYTable*   table  = (_HYTable*)GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);
    if (index<table->verticalSpaces.lLength) {
        _String* varName = (_String*) table->GetCellData(1,index);
        result  = FetchVar (LocateVarByName (*varName));
    }

    return      result;
}

//__________________________________________________________

_CalcNode* _HYParameterTable::GrabRowsCalcNode (_String* varName)
{
    long k = varName->FindBackwards('.',0,-1);
    if (k>=0) {
        _String    nodeName = varName->Cut(0,k-1);
        return     (_CalcNode*)FetchVar(LocateVarByName (nodeName));
    }
    return      nil;
}

//__________________________________________________________

void    _HYParameterTable::UpdateKthRow (long index, bool mark)
{
    _Variable* thisV = (_Variable*)RetrieveIthRowVar (index);
    if (thisV) {
        long       idx = -1;
        _HYTable*   table  = (_HYTable*)GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);
        if (!(thisV->IsCategory()||(thisV->ObjectClass()==TREE))) {
            _String     newV (thisV->Compute()->Value());
            idx = index*table->horizontalSpaces.lLength+2;
            table->SetCellData (&newV,index,2,table->cellTypes.lData[idx],true);
        } else if (thisV->IsCategory()) {
            long      catIdx = 0;
            _String   *cellValue = (_String*)table->GetCellData (2,index);
            catIdx = cellValue->Cut (0,cellValue->Find ('(')-1).toNum();
            idx = index*table->horizontalSpaces.lLength+3;
            ((_CategoryVariable*)thisV)->Refresh();
            _String   cv = _String (((_CategoryVariable*)thisV)->GetIntervalValue(catIdx-1));
            table->SetCellData (&cv,index,3,HY_TABLE_STATIC_TEXT|HY_TABLE_PULLDOWN,true);
        } else {
            return;
        }

        if (mark&&(idx>=0)) {
            table->_MarkCellForUpdate (idx);
        }
    }
}

//__________________________________________________________

void    _HYParameterTable::HandleCellEditEvent (long index)
{

    _HYTable*       table = (_HYTable*)GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);

    long h = index%table->horizontalSpaces.lLength,
         v = index/table->horizontalSpaces.lLength,
         k;

    _Variable * thisV = (_Variable*) RetrieveIthRowVar (v);
    _String   * currentValue,
              warningString;
    _SimpleList taintCells;

    _Parameter  p1, p2, b1, b2;

    bool        runUpdate = false,
                ind = thisV->IsIndependent();

    currentValue = (_String*)table->GetCellData (h,v);

    switch (h) {
    case 2: { // parameter value column
        if (thisV->IsIndependent()) {
            _Formula f (*currentValue);
            p1 = thisV->Compute()->Value();
            _PMathObj  fv = nil;

            if (!f.IsEmpty()) {
                fv = f.Compute();
                if (fv) {
                    p2 = fv->Value();
                    b1 = thisV->GetLowerBound();
                    b2 = thisV->GetUpperBound();
                }
            }

            if (!fv) {
                warningString = _String("Failed to evaluate the formula as entered. The value will be reset to ")& p1 &'.';
            } else {
                if ((p2>=b1)&&(p2<=b2)) {
                    _Constant   c (p2);
                    warningString = *thisV->GetName() & '=' & _String(thisV->Compute()->Value()) & ';';
                    if (doUndo) {
                        _UpdateUndoMenu (&warningString, &undoEditValue);
                    } else {
                        *(_String*)(undoCommands (undoCommands.lLength-1)) << warningString;
                    }

                    thisV->SetValue(&c);
                    warningString = p2;
                    table->SetCellData(&warningString,v,2,table->cellTypes.lData[index],true);
                    taintCells << index;
                    warningString = empty;
                    runUpdate = true;
                } else {
                    warningString = _String ("The value ")& *currentValue & " is out of bounds ["
                                    & b1 & ',' & b2 &"] for the variable "& *thisV->GetName() &
                                    ". The value will be reset to "& p1 &'.';
                }
            }
            if (warningString.sLength) {
                _String         rv (p1);
                table->SetCellData(&rv,v,2,table->cellTypes.lData[index],true);
                taintCells << index;
                ProblemReport (warningString,(Ptr)this);
            }
        }
        break;
    }

    case 3: { // constraint column
        if (ind||((!ind)&&(currentValue->sLength)))
            // indep=>dep or dep=>dep
        {
            _Formula newF (*currentValue,nil);
            if (!newF.IsEmpty()) {
                if (newF.CheckFForDependence (thisV->GetAVariable(),true)) {
                    warningString = _String ("The constraint '")& *thisV->GetName() &":="& *currentValue & "' contains self-referencing and is thus disallowed.";
                    if (ind) {
                        table->SetCellData(&empty,v,3,table->cellTypes.lData[index],true);
                    } else {
                        table->SetCellData(thisV->GetFormulaString(),v,3,table->cellTypes.lData[index],false);
                    }
                    taintCells << index;
                    ProblemReport (warningString,(Ptr)this);
                    break;
                }

                warningString = *thisV->GetName() & '=' & _String(thisV->Compute()->Value()) & ';';
                if (doUndo) {
                    _UpdateUndoMenu (&warningString, &undoSetConstraint);
                } else {
                    (*(_String*)undoCommands (undoCommands.lLength-1)) <<  warningString;
                }

                thisV->SetFormula (newF);
                if (ind) {
                    if (thisV->IsGlobal()) {
                        table->SetCellData(&constrainedGlobalTableIcon,v,0,table->cellTypes.lData[index-3],true);
                    } else {
                        table->SetCellData(&constrainedTableIcon,v,0,table->cellTypes.lData[index-3],true);
                    }
                    taintCells << index-3;
                    if (table->cellTypes.lData[index-1]&HY_TABLE_EDIT_TEXT) {
                        table->cellTypes.lData[index-1] += HY_TABLE_STATIC_TEXT;
                        table->cellTypes.lData[index-1] -= HY_TABLE_EDIT_TEXT;
                    }
                }
                warningString = thisV->Compute()->Value();
                table->SetCellData(&warningString,v,2,table->cellTypes.lData[index-1],true);
                taintCells << index-1;
                runUpdate = true;
                if (doTableUpdate) {
                    StretchColumnToFit (3);
                    ToggleConstrainedView();
                }

            } else {
                table->SetCellData(&empty,v,3,table->cellTypes.lData[index],true);
                taintCells << index;
            }
        } else
            // dep=>ind
        {
            if ((!ind)&&(currentValue->sLength==0)) {
                _Constant cv (thisV->Compute()->Value());
                currentValue = (_String*)thisV->GetFormulaString();
                warningString = *thisV->GetName() & ":=" & *currentValue & ';';
                DeleteObject (currentValue);
                if (doUndo) {
                    _UpdateUndoMenu (&warningString, &undoKillConstraint);
                } else {
                    (*(_String*)undoCommands (undoCommands.lLength-1)) << warningString;
                }
                thisV->SetValue (&cv);
                table->SetCellData(thisV->IsGlobal()?&globalTableIcon:&localTableIcon,
                                   v,0,table->cellTypes.lData[index-3],true);
                taintCells << index-3;
                if (table->cellTypes.lData[index-1]&HY_TABLE_STATIC_TEXT) {
                    table->cellTypes.lData[index-1] -= HY_TABLE_STATIC_TEXT;
                    table->cellTypes.lData[index-1] += HY_TABLE_EDIT_TEXT;
                }
                taintCells << index-1;
                runUpdate = true;
                if (avViewOptions&HY_PARAMETER_TABLE_VIEW_CONSTRAINED) {
                    avViewOptions = GetViewOptions();
                    if (!(avViewOptions&HY_PARAMETER_TABLE_VIEW_CONSTRAINED))
                        if (viewOptions&HY_PARAMETER_TABLE_VIEW_CONSTRAINED) {
                            viewOptions -= HY_PARAMETER_TABLE_VIEW_CONSTRAINED;
                        }

                    _UpdateViewMenu();
                }
            }

        }
        if (doTableUpdate) {
            VerifyGlobalVariables();
        }
        break;
    }
    }

    if (runUpdate) {
        if (!thisV->IsIndependent()) {
            _SimpleList indVars;
            {
                _AVLList ia (&indVars);
                if (thisV->varFormula) {
                    thisV->varFormula->ScanFForVariables(ia,false,true);
                }
                ia.ReorderList ();
            }
            h = -1;
            for (k=0; k<indVars.lLength; k++) {
                thisV = LocateVar (indVars.lData[k]);
                if (thisV->IsIndependent()) {
                    h = k;
                    if (!thisV->IsGlobal()) {
                        break;
                    }
                }
            }
            if (h>=0) {
                thisV = LocateVar (indVars.lData[h]);
                _Constant c (thisV->Compute()->Value());
                thisV->SetValue (&c);
            } else {
                _LikelihoodFunction* lf = (_LikelihoodFunction*)likeFuncList (lfID);
                if (lf->CountObjects (1)+lf->CountObjects (2)) {
                    lf->SetIthIndependent (0,lf->GetIthIndependent (0));
                }
            }

        }
        for (k = 1; k<table->verticalSpaces.lLength; k++) {
            _Variable * vr = (_Variable*) RetrieveIthRowVar (k);
            if ((!vr->IsIndependent())&&vr->HasChanged()) {
                p1 = vr->Compute()->Value();
                warningString = p1;
                index = k*table->horizontalSpaces.lLength+2;
                table->SetCellData (&warningString,k,2,table->cellTypes.lData[index],true);
                taintCells << index;
            } else {
                if (vr->IsCategory())
                    if (vr->HasChanged()) {
                        ((_CategoryVariable*)vr)->Refresh (true);
                        UpdateKthRow (k,true);
                    }
                /*else
                {
                    _SimpleList vars;
                    vr->ScanForVariables (vars,true);
                    if (vars.Find (thisV->GetAVariable())>=0)
                    {
                        ((_CategoryVariable*)vr)->Refresh (true);
                        UpdateKthRow (k,true);
                    }
                }*/
            }
        }

        if (autoUpdateTableLF) {
            UpdateLogLikelihood();
        }
    }

    if (taintCells.lLength) {
        table->_MarkCellsForUpdate(taintCells);
    }
}

//__________________________________________________________

void    _HYParameterTable::HandleCellPullDown (long index, long location)
{

    _HYTable*       table = (_HYTable*)GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);

    long h = index%table->horizontalSpaces.lLength,
         v = index/table->horizontalSpaces.lLength,
         k;

    _Variable * thisV = (_Variable*) RetrieveIthRowVar (v);
    if (thisV&&thisV->IsCategory()) {
        _CategoryVariable* thisCV = (_CategoryVariable*) thisV;
        _List              menuItems;
        _String            kthOption,
                           *cellValue;

        for (k = 1; k <= thisCV->GetNumberOfIntervals(); k++) {
            kthOption = _String ("Rate Class ") & k;
            menuItems && & kthOption;
        }
        if (!thisCV->GetDensity().IsEmpty()) {
            menuItems && & menuSeparator;
            kthOption = "Density Plot";
            menuItems && & kthOption;
        }

        cellValue = (_String*) table->GetCellData (2,v);

        k = (cellValue->Cut(0,cellValue->Find ('(')-1)).toNum();

        h = table->_HandlePullDown  (menuItems,(location&0xffff0000)>>16,location&0x0000ffff,k);

        if ((h!=-1)&&(h!=(k-1))) {
            if (h<thisCV->GetNumberOfIntervals()) {
                kthOption = _String ((long)(h+1))&'('&thisCV->GetNumberOfIntervals()&')';
                table->SetCellData (&kthOption,v,2,HY_TABLE_STATIC_TEXT,true);
                kthOption = thisCV->GetIntervalValue(h);
                table->SetCellData (&kthOption,v,3,HY_TABLE_STATIC_TEXT|HY_TABLE_PULLDOWN,true);
                _SimpleList cells;
                cells << v*table->horizontalSpaces.lLength+2;
                cells << v*table->horizontalSpaces.lLength+3;
                table->_MarkCellsForUpdate(cells);
                StretchColumnToFit (3);
            } else if (h==thisCV->GetNumberOfIntervals()+1) { // produce a plot
                OpenCategoryPlotWindow(thisCV);
            }
        }
    }
}

//__________________________________________________________

void    _HYParameterTable::HandleHeaderPullDown (long index, long location)
{

    _HYTable*       table = (_HYTable*)GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW-1,0);
    _HYTable*       table2 = (_HYTable*)GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);

    long h = index%table->horizontalSpaces.lLength,
         v = index/table->horizontalSpaces.lLength;

    //_Variable * thisV = (_Variable*) RetrieveIthRowVar (v);
    if ((v==0)&&(h==1)) {
        _List              menuItems;
        _String            kthOption ("Sort in Ascending Order");
        menuItems && & kthOption;
        kthOption = "Sort in Descending Order";
        menuItems && & kthOption;

        h = table->_HandlePullDown  (menuItems,(location&0xffff0000)>>16,location&0x0000ffff,0);

        if (h>=0) {
            menuItems.Clear();
            _SimpleList  indices;
            for (v=0; v<table2->verticalSpaces.lLength; v++) {
                menuItems << table2->GetCellData(1,v);
                indices << v;
            }
            SortLists (&menuItems,&indices);
            if (h) {
                indices.Flip();
            }
            table2->SetRowOrder (indices);
        }
    } else if ((v==0)&&(h==0)) {
        _List              menuItems;
        _String            kthOption ("Sort by type");
        menuItems && & kthOption;

        h = table->_HandlePullDown  (menuItems,(location&0xffff0000)>>16,location&0x0000ffff,0);

        if (h>=0) {
            _SimpleList  indices,
                         types;

            for (v=0; v<table2->verticalSpaces.lLength; v++) {
                _Variable * theV = (_Variable*)RetrieveIthRowVar(v);
                if (theV->ObjectClass()==TREE) {
                    types << 0;
                } else {
                    if (theV->IsGlobal()) {
                        if (theV->IsIndependent()) {
                            types<<1;
                        } else {
                            types<<2;
                        }
                    } else if (theV->IsIndependent()) {
                        types<<3;
                    } else {
                        types<<4;
                    }
                }
                indices << v;
            }
            SortLists (&types,&indices);
            table2->SetRowOrder (indices);
        }
    }

}

//__________________________________________________________

void    _HYParameterTable::UpdateSelectionDependentButtons (void)
{

    _HYTable*       table = (_HYTable*)GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);
    _HYButtonBar*   bb    = (_HYButtonBar*)GetCellObject (HY_PARAMETER_TABLE_BUTTON_ROW,0);

    bool            b1 = false,
                    b2 = false,
                    b3 = false,
                    b4 = false,
                    b5 = false,
                    b6 = false,
                    b7 = false;

    _SimpleList     rows;
    _Variable*      thisV,
                    *      thisV2 = nil;
    long            counterL,
                    counterV,
                    counterC,
                    counterT,
                    k;

    table->GetRowSelection (rows);

    if (rows.lLength>=1) {
        thisV = (_Variable*)RetrieveIthRowVar (rows.lData[0]);

        if ((thisV->ObjectClass()==TREE)&&(rows.lLength == 1)) {
            b5 = true;
        } else if (thisV->IsIndependent()) {
            b4 = true;
            _CalcNode* cNode = GrabRowsCalcNode (thisV->GetName());
            if (cNode) {
                node<long>*nodeD    = cNode->LocateMeInTree();
                if (nodeD&&nodeD->get_num_nodes()) {
                    b5 = true;
                }
            }

            if (b5)
                for (long k=1; k<rows.lLength; k++) {
                    thisV = (_Variable*)RetrieveIthRowVar (rows.lData[k]);
                    _CalcNode* cNode2 = GrabRowsCalcNode (thisV->GetName());
                    if (cNode2 != cNode) {
                        b5 = false;
                        break;
                    }
                }
        }
    }

    if (rows.lLength) {
        counterL = 0;
        counterV = 0;
        counterC = 0;
        counterT = 0;
        for (k=0; k<rows.lLength; k++) {
            thisV = (_Variable*)RetrieveIthRowVar (rows.lData[k]);
            if (thisV->ObjectClass()==TREE) {
                counterT++;
            } else if (!thisV->IsCategory()) {
                counterV++;
                if (!thisV->IsGlobal()) {
                    counterL++;
                }
                if (!thisV->IsIndependent()) {
                    counterC++;
                }
            }
        }
        if (counterV >= 1) {
            b2 = true;
        }

        if (counterV == 2) {
            b3 = true;
        }

        if (counterL) {
            b1 = true;
        }

        if (counterC) {
            b7 = true;
        }

        if (rows.lLength==2) {
            if (counterT==0) {
                thisV =  (_Variable*)GrabRowsCalcNode ((_String*)table->GetCellData(1,rows.lData[0]));
                thisV2 = (_Variable*)GrabRowsCalcNode ((_String*)table->GetCellData(1,rows.lData[1]));
                if (!(thisV&&thisV2)) {
                    thisV2 = nil;
                } else if (!((_CalcNode*)thisV)->MatchSubtree ((_CalcNode*)thisV2)) {
                    thisV2 = nil;
                }
            } else if (counterT==2) {
                thisV = (_Variable*)RetrieveIthRowVar (rows.lData[0]);
                thisV2 = (_Variable*)RetrieveIthRowVar (rows.lData[1]);
            }
            if (thisV2) {
                b6 = true;
            }
        }
    }

    bb->EnableButton (0,b1);
    bb->EnableButton (1,b2);
    bb->EnableButton (2,b3);
    bb->EnableButton (3,b4);
    bb->EnableButton (4,b5);
    bb->EnableButton (5,b6);
    bb->EnableButton (6,b7);

}

//__________________________________________________________

void    _HYParameterTable::DoEqualConstraint (void)
{

    _HYTable*       table = (_HYTable*)GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);
    _HYButtonBar*   bb    = (_HYButtonBar*)GetCellObject (HY_PARAMETER_TABLE_BUTTON_ROW,0);

    int  h,v;
    long k;

    _List           menuOptions;
    _SimpleList     rows, tainted;
    _String         menuChoice;

    table->GetRowSelection (rows);

    if (rows.lLength>1) {
        for (h=0; h<rows.lLength; h++) {
            menuChoice = _String (":=") & *(_String*)table->GetCellData(1,rows.lData[h]);
            menuOptions && & menuChoice;
        }
        menuOptions && & menuSeparator;
    }

    menuChoice = "Equal to own value";
    menuOptions && & menuChoice;
    menuChoice = "Equal to expression";
    menuOptions && & menuChoice;


    bb->GetButtonLoc(1,h,v,true);
    menuChoice = HandlePullDown (menuOptions,h,v,0);
    bb->_UnpushButton();

    v = menuOptions.Find (&menuChoice);
    if (v<0) {
        return;
    }

    doUndo            = false;
    autoUpdateTableLF = false;
    doTableUpdate     = false;

    undoDescriptions && & undoEqualConstraints;
    undoCommands << new _String (128L, true);

    if ((v<rows.lLength)&&(rows.lLength>1)) {
        menuChoice = *(_String*)table->GetCellData(1,rows.lData[v]);
        for (h=0; h<rows.lLength; h++) {
            if (v==h) {
                continue;
            }
            k = table->horizontalSpaces.lLength*rows.lData[h] + 3;
            table->SetCellData (&menuChoice,rows.lData[h],3,table->cellTypes[k],true);
            ProcessEvent (generateTableEditCellEvent(table->GetID(),k));
            tainted << k;
        }
    } else {
        if (v == menuOptions.lLength-1) {
            _String prompt = "Expression for the constraint";
            if (!EnterStringDialog (menuChoice, prompt, (Ptr)this)) {
                ((_String*)undoCommands(undoCommands.lLength-1))->Finalize();
                autoUpdateTableLF = true;
                doUndo            = true;
                doTableUpdate     = true;
                undoDescriptions.Delete (undoDescriptions.lLength-1);
                undoCommands.Delete (undoCommands.lLength-1);
                return;
            }

            for (h=0; h<rows.lLength; h++) {
                k = table->horizontalSpaces.lLength*rows.lData[h] + 3;
                table->SetCellData (&menuChoice,rows.lData[h],3,table->cellTypes[k],true);
                ProcessEvent (generateTableEditCellEvent(table->GetID(),k));
                tainted << k;
            }
        } else
            for (h=0; h<rows.lLength; h++) {
                k = table->horizontalSpaces.lLength*rows.lData[h] + 3;
                menuChoice = *(_String*)table->GetCellData(2,rows.lData[h]);
                table->SetCellData (&menuChoice,rows.lData[h],3,table->cellTypes[k],true);
                ProcessEvent (generateTableEditCellEvent(table->GetID(),k));
                tainted << k;
            }
    }

    ((_String*)undoCommands(undoCommands.lLength-1))->Finalize();
    autoUpdateTableLF = true;
    doUndo            = true;
    doTableUpdate     = true;

    StretchColumnToFit   (3);
    ToggleConstrainedView( );
    VerifyGlobalVariables( );

    table->_MarkCellsForUpdate (tainted);
    _UpdateUndoMenu (nil,nil);
    UpdateLogLikelihood();
    UpdateSelectionDependentButtons();
}

//__________________________________________________________

void    _HYParameterTable::DoClearConstraints (void)
{

    _HYTable*       table = (_HYTable*)GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);
    _SimpleList     rows,
                    tCells;
    long            i,j;
    _Variable*      thisV;

    table->GetRowSelection (rows);

    doUndo            = false;
    autoUpdateTableLF = false;
    doTableUpdate     = false;
    undoDescriptions && & undoClearConstraints;
    undoCommands << new _String (128L, true);
    for (i=0; i<rows.lLength; i++) {
        thisV = (_Variable*)RetrieveIthRowVar(rows.lData[i]);
        if (thisV->ObjectClass()==TREE) {
            continue;
        }
        if (!thisV->IsIndependent()) {
            j = rows.lData[i]*table->horizontalSpaces.lLength + 3;
            table->SetCellData(&empty,rows.lData[i],3,table->cellTypes[j],true);
            ProcessEvent(generateTableEditCellEvent(table->GetID(),j));
            tCells<<j;
        }
    }
    autoUpdateTableLF = true;
    doUndo            = true;
    doTableUpdate     = true;



    ((_String*)undoCommands(undoCommands.lLength-1))->Finalize();
    _UpdateUndoMenu (nil,nil);
    table->_MarkCellsForUpdate(tCells);

    StretchColumnToFit   (3);
    ToggleConstrainedView( );
    VerifyGlobalVariables( );

    UpdateSelectionDependentButtons();
    UpdateLogLikelihood ();
}

//__________________________________________________________

void    _HYParameterTable::DoProportionalConstraint (void)
{

    _HYTable*       table = (_HYTable*)    GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);
    _HYButtonBar*   bb    = (_HYButtonBar*)GetCellObject (HY_PARAMETER_TABLE_BUTTON_ROW,0);

    int  h,v;
    long k;

    _List           menuOptions,
                    globalVars;
    _SimpleList     rows;
    _String         menuChoice,
                    *var1,
                    *var2;

    table->GetRowSelection (rows);

    if ((lfID>=0)&&(rows.lLength==2)) {

        long    insertAt = GatherListOfGlobals (globalVars);

        var1 = (_String*)table->GetCellData(1,rows.lData[0]);
        var2 = (_String*)table->GetCellData(1,rows.lData[1]);

        for (h=0; h<globalVars.lLength; h++) {
            menuChoice = *var1 & ":=" & *(_String*)globalVars(h) & '*' & *var2;
            menuOptions && & menuChoice;
            menuChoice = *var2 & ":=" & *(_String*)globalVars(h) & '*' & *var1;
            menuOptions && & menuChoice;
        }

        if (globalVars.lLength) {
            menuOptions && & menuSeparator;
        }

        menuChoice = *var1 & ":= {New Ratio} *" & *var2;
        menuOptions && & menuChoice;
        menuChoice = *var2 & ":= {New Ratio} *" & *var1;
        menuOptions && & menuChoice;

        menuOptions && & menuSeparator;

        menuChoice = *var1 & ":= {Constant} *" & *var2;
        menuOptions && & menuChoice;
        menuChoice = *var2 & ":= {Constant} *" & *var1;
        menuOptions && & menuChoice;

        bb->GetButtonLoc(2,h,v,true);
        menuChoice = HandlePullDown (menuOptions,h,v,0);
        bb->_UnpushButton();

        v = menuOptions.Find (&menuChoice);

        if (v<0) {
            return;
        }

        if (v>=2*globalVars.lLength+(globalVars.lLength!=0)) {
            _String          newID,
                             prStr;
            if (v<=2*globalVars.lLength+1+(globalVars.lLength!=0)) {
                prStr = "Identifier of the new ratio variable:";
                if (!EnterStringDialog (newID,prStr, (Ptr)this, hyIDValidator)) {
                    return;
                } else {
                    if ((h=globalVars.Find (&newID))>=0) {
                        menuChoice = menuChoice.Replace ("{New Ratio}", newID, true);
                    } else {
                        if (!newID.IsValidIdentifier()) {
                            return;
                        }
                        _Variable newVar (newID,true);
                        LocateVar(newVar.GetAVariable())->SetNumericValue(1.);
                        menuChoice = menuChoice.Replace ("{New Ratio}", *newVar.GetName(), true);
                        table->AddRow (insertAt,HY_PARAMETER_TABLE_ROW_HEIGHT,HY_TABLE_STATIC_TEXT);
                        for (long kk=0; kk<rows.lLength; kk++)
                            if (rows.lData[kk]>=insertAt) {
                                rows.lData[kk]++;
                            }

                        SetRowValues (insertAt,LocateVar(newVar.GetAVariable()),table);
                        table->_MarkForUpdate();
                    }
                }
            } else {
                prStr = "Specify the constant of proportionality:";
                if (!EnterStringDialog (newID,prStr, (Ptr)this)) {
                    return;
                } else {
                    newID = _String(newID.toNum());
                    menuChoice = menuChoice.Replace ("{Constant}", newID, true);
                }
            }
        }

        h = menuChoice.beginswith (*var1)?0:1;
        menuChoice = menuChoice.Cut (menuChoice.Find ('=')+1,-1);
        k = table->horizontalSpaces.lLength*rows.lData[h] + 3;
        table->SetCellData (&menuChoice,rows.lData[h],3,table->cellTypes[k],true);
        ProcessEvent (generateTableEditCellEvent(table->GetID(),k));
        ToggleConstrainedView ();
        ToggleConstrainedView (HY_PARAMETER_TABLE_VIEW_GLOBAL);
        table->_MarkCellForUpdate (k);
        UpdateSelectionDependentButtons();
    }
}

//__________________________________________________________

void    _HYParameterTable::DoCleanUp (void)
{
    if (lfID>=0) {

        _LikelihoodFunction*  theLF =
            (_LikelihoodFunction*)likeFuncList (lfID);

        theLF->RescanAllVariables();
        _SimpleList     allUsedVars (theLF->GetIndependentVars()),
                        rowsToDelete;

        allUsedVars << theLF->GetDependentVars();

        allUsedVars.Sort();

        _HYTable*       table = (_HYTable*)    GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);

        long            k;

        for (k=0; k<table->verticalSpaces.lLength; k++) {
            _Variable* thisV = (_Variable*)RetrieveIthRowVar(k);
            if (thisV->ObjectClass()!=TREE)
                if (!thisV->IsCategory ())
                    if (allUsedVars.BinaryFind (thisV->GetAVariable())<0) {
                        rowsToDelete << k;
                    }
        }

        if (rowsToDelete.lLength) {

            undoCommands.Clear();
            undoDescriptions.Clear();
            _UpdateUndoMenu (nil,nil);

            for (k=rowsToDelete.lLength-1; k>=0; k--) {
                table->DeleteRow (rowsToDelete.lData[k]);
            }

            UpdateComponentInfo();
            table->SetVisibleSize (table->rel);
            table->_MarkForUpdate();


            UpdateSelectionDependentButtons();
        }
        UpdateLogLikelihood ();
    }
}

//__________________________________________________________

void    _HYParameterTable::DoProportionalSubtrees (void)
{

    _HYTable*       table = (_HYTable*)    GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);
    _HYButtonBar*   bb    = (_HYButtonBar*)GetCellObject (HY_PARAMETER_TABLE_BUTTON_ROW,0);

    int             h,v;
    long            k;

    _CalcNode       *thisCN,
                    *thisCN2;

    _List           menuOptions,
                    menuOptions2,
                    undoValues;

    _SimpleList     rows;

    _String         menuChoice;

    _Variable       *thisTree,
                    *thisTree2;

    table->GetRowSelection (rows);

    if ((lfID>=0)&&(rows.lLength==2)) {
        bb->_UnpushButton();
        thisCN = GrabRowsCalcNode(((_Variable*)RetrieveIthRowVar (rows.lData[0]))->GetName());
        if (thisCN) {
            h = thisCN->GetName()->Find ('.');
            if (h<0) {
                return;
            }
            menuChoice = thisCN->GetName()->Cut (0,h-1);
            thisTree = FetchVar (LocateVarByName (menuChoice));
            thisCN2 = GrabRowsCalcNode(((_Variable*)RetrieveIthRowVar (rows.lData[1]))->GetName());
            h = thisCN2->GetName()->Find ('.');
            if (h<0) {
                return;
            }
            menuChoice = thisCN2->GetName()->Cut (0,h-1);
            thisTree2 = FetchVar (LocateVarByName (menuChoice));
        } else {
            menuChoice = *(_String*)table->GetCellData (1,rows.lData[0]);
            thisTree = FetchVar (LocateVarByName (menuChoice));
            menuChoice = *(_String*)table->GetCellData (1,rows.lData[1]);
            thisTree2 = FetchVar (LocateVarByName (menuChoice));
            thisCN2 = nil;
        }

        if (!(thisTree&&(thisTree->ObjectClass()==TREE)&&
                thisTree2&&(thisTree2->ObjectClass()==TREE))) {
            return;
        }

        ((_TheTree*)thisTree) ->ScanSubtreeVars (menuOptions, 1,thisCN);
        ((_TheTree*)thisTree2)->ScanSubtreeVars (menuOptions2,1,thisCN2);

        for (k=menuOptions.lLength-1; k>=0; k--)
            if (menuOptions2.Find (menuOptions(k))<0) {
                menuOptions.Delete (k);
            }

        if (menuOptions.lLength>1) {
            menuOptions && & menuSeparator;
            menuOptions && & allVariablesOption;
        } else if (menuOptions.lLength==0) {
            menuChoice = "There are no parameters eligible for the relative ratio constraint in the selected subtrees.";
            ProblemReport (menuChoice,(Ptr)this);
            return;
        }


        bb->GetButtonLoc(5,h,v,true);
        menuChoice = HandlePullDown (menuOptions,h,v,0);
        bb->_UnpushButton();

        v = menuOptions.Find (&menuChoice);

        if (v<0) {
            return;
        }

        _String prStr = "Identifier of the ratio variable:",
                newID;

        if (!EnterStringDialog (newID,prStr, (Ptr)this)) {
            return;
        } else {
            if (!newID.IsValidIdentifier()) {
                return;
            }

            _List globalVars;
            long insertAt = GatherListOfGlobals (globalVars);
            if (globalVars.Find(&newID)<0) {
                _Variable newVar (newID,true);
                LocateVar(newVar.GetAVariable())->SetNumericValue(1.);
                table->AddRow (insertAt,HY_PARAMETER_TABLE_ROW_HEIGHT,HY_TABLE_STATIC_TEXT);
                SetRowValues (insertAt,LocateVar(newVar.GetAVariable()),table);
            }
        }

        _LikelihoodFunction *lf = (_LikelihoodFunction*)likeFuncList (lfID);
        rows.Clear();
        rows.Duplicate (&lf->GetIndependentVars());
        for (v=0; v<rows.lLength; v++) {
            _Variable* thisV = LocateVar (rows.lData[v]);
            prStr = thisV->Compute()->Value();
            undoValues && & prStr;
        }
        if (thisCN) {
            if (menuChoice.Equal (&allVariablesOption)) {
                menuChoice = blReplicate&"\"this1.?:="&newID&"*this2.?\","&*thisCN->GetName()&','&*thisCN2->GetName()&");";
            } else {
                menuChoice = blReplicate&"\"this1.?."&menuChoice&":="&newID&"*this2.?."&menuChoice&"\","&*thisCN->GetName()&','&*thisCN2->GetName()&");";
            }
        } else {
            if (menuChoice.Equal (&allVariablesOption)) {
                menuChoice = blReplicate&"\"this1.?:="&newID&"*this2.?\","&*thisTree->GetName()&','&*thisTree2->GetName()&");";
            } else {
                menuChoice = blReplicate&"\"this1.?."&menuChoice&":="&newID&"*this2.?."&menuChoice&"\","&*thisTree->GetName()&','&*thisTree2->GetName()&");";
            }
        }
        _ExecutionList ex (menuChoice);
        ex.Execute();
        _String     undoS (16,true);
        for (v=0; v<rows.lLength; v++) {
            _Variable* thisV = LocateVar (rows.lData[v]);
            if (!thisV->IsIndependent()) {
                undoS << thisV->GetName();
                undoS << '=';
                undoS << (_String*)undoValues (v);
                undoS << ';';
                undoS << '\n';
            }
        }
        undoS.Finalize ();
        _UpdateUndoMenu (&undoS,&undoProportionalSubtrees);
        TaintTheLF();
        RefreshTheTable();
        VerifyGlobalVariables();
        StretchColumnToFit(3);
        table->_MarkForUpdate();
        UpdateSelectionDependentButtons();
    }
}

//__________________________________________________________

void    _HYParameterTable::DoBoundsChange (void)
{

    _HYTable*       table = (_HYTable*)    GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);
    _HYButtonBar*   bb    = (_HYButtonBar*)GetCellObject (HY_PARAMETER_TABLE_BUTTON_ROW,0);

    int  h,v;

    _List           menuOptions;
    _SimpleList     rows;
    _String         menuChoice,
                    prStr,
                    newID;
    _Variable       *thisV;
    _Parameter      newB;

    table->GetRowSelection (rows);

    if ((lfID>=0)&&(rows.lLength==1)) {
        thisV = (_Variable*)RetrieveIthRowVar (rows.lData[0]);
        menuChoice = _String("Lower Bound [")& thisV->GetLowerBound()&']';
        menuOptions && & menuChoice;
        menuChoice = _String("Upper Bound [")& thisV->GetUpperBound()&']';
        menuOptions && & menuChoice;

        bb->GetButtonLoc(3,h,v,true);
        menuChoice = HandlePullDown (menuOptions,h,v,0);
        bb->_UnpushButton();

        v = menuOptions.Find (&menuChoice);

        if (v<0) {
            return;
        }

        if (v==0) {
            prStr = "Specify the new lower bound:";
            newID = thisV->GetLowerBound();
        } else {
            prStr = "Specify the new upper bound:";
            newID = thisV->GetUpperBound();
        }

        if (!EnterStringDialog (newID,prStr, (Ptr)this)) {
            return;
        }

        newB = newID.toNum();
        h = 0;

        if (v==0) {
            if (newB>thisV->GetUpperBound()) {
                prStr = "The lower bound cannot exceed the upper bound!";
                ProblemReport(prStr,(Ptr)this);
                return;
            }
            if (newB>thisV->Compute()->Value()) {
                h = 1;
            }
            thisV->SetBounds (newB,thisV->GetUpperBound());
        } else {
            if (newB<thisV->GetLowerBound()) {
                prStr = "The upper bound cannot be less than the lower bound!";
                ProblemReport(prStr,(Ptr)this);
                return;
            }
            if (newB<thisV->Compute()->Value()) {
                h = 1;
            }
            thisV->SetBounds (thisV->GetLowerBound(),newB);
        }

        if (h) {
            newID (newB);
            h = rows.lData[0]*table->horizontalSpaces.lLength+2;
            table->SetCellData (&newID,rows.lData[0],2,table->cellTypes.lData[h],true);
            ProcessEvent (generateTableEditCellEvent(table->GetID(),h));
        }

    }
}

//__________________________________________________________

void    _HYParameterTable::DoMolecularClock (void)
{

    _HYTable*       table = (_HYTable*)    GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);
    _HYButtonBar*   bb    = (_HYButtonBar*)GetCellObject (HY_PARAMETER_TABLE_BUTTON_ROW,0);

    int  h,v;

    _List           menuOptions,
                    undoValues;
    _SimpleList     rows;
    _String         menuChoice,
                    prStr,
                    newID;
    _CalcNode       *thisCN;
    _Variable       *thisTree;

    table->GetRowSelection (rows);

    if ((lfID>=0)&&(rows.lLength>=1)) {
        bb->_UnpushButton();
        thisCN = GrabRowsCalcNode(((_Variable*)RetrieveIthRowVar (rows.lData[0]))->GetName());
        if (thisCN) {
            h = thisCN->GetName()->Find ('.');
            if (h<0) {
                return;
            }
            menuChoice = thisCN->GetName()->Cut (0,h-1);
            thisTree = FetchVar (LocateVarByName (menuChoice));
        } else {
            menuChoice = *(_String*)table->GetCellData (1,rows.lData[0]);
            thisTree = FetchVar (LocateVarByName (menuChoice));
        }

        if (!(thisTree&&(thisTree->ObjectClass()==TREE))) {
            return;
        }
        ((_TheTree*)thisTree)->ScanSubtreeVars (menuOptions,1,thisCN);

        if (menuOptions.lLength>1) {
            menuOptions && & menuSeparator;
            menuOptions && & allVariablesOption;
        } else if (menuOptions.lLength==0) {
            menuChoice = "There are no parameters eligible for the molecular clock constraint in the selected subtree.";
            ProblemReport (menuChoice,(Ptr)this);
            return;
        }

        bb->GetButtonLoc(4,h,v,true);
        menuChoice = HandlePullDown (menuOptions,h,v,0);

        v = menuOptions.Find (&menuChoice);
        if (v<0) {
            return;
        }

        if (v<menuOptions.lLength-2) {
            menuChoice = *(_String*)menuOptions(v);
            menuOptions.Clear();
            menuOptions && & menuChoice;
        } else {
            if (menuOptions.lLength>1) {
                menuOptions.Delete (menuOptions.lLength-1);
                menuOptions.Delete (menuOptions.lLength-1);
            }
        }

        menuOptions.InsertElement(&empty,0,true);
        if (thisCN) {
            menuChoice = *thisCN->GetName();
            menuChoice.Trim (menuChoice.Find('.')+1,-1);
        } else {
            menuChoice = empty;
        }

        _LikelihoodFunction *lf = (_LikelihoodFunction*)likeFuncList (lfID);
        rows.Clear();
        rows.Duplicate (&lf->GetIndependentVars());
        for (v=0; v<rows.lLength; v++) {
            _Variable* thisV = LocateVar (rows.lData[v]);
            prStr = thisV->Compute()->Value();
            undoValues && & prStr;
        }

        ((_TheTree*)thisTree)->MolecularClock(menuChoice,menuOptions);

        _String     undoS (16,true);
        for (v=0; v<rows.lLength; v++) {
            _Variable* thisV = LocateVar (rows.lData[v]);
            if (!thisV->IsIndependent()) {
                undoS << thisV->GetName();
                undoS << '=';
                undoS << (_String*)undoValues (v);
                undoS << ';';
                undoS << '\n';
            }
        }
        undoS.Finalize ();

        _UpdateUndoMenu (&undoS,&undoMolecularClock);

        ToggleConstrainedView();
        TaintTheLF           ();
        RefreshTheTable();
        UpdateSelectionDependentButtons();
        StretchColumnToFit(3);
        table->_MarkContentsForUpdate();
    }
}

//__________________________________________________________

void    _HYParameterTable::ToggleConstrainedView (char flag)
{
    if (!(avViewOptions&flag)) {
        avViewOptions = GetViewOptions();
        if (avViewOptions&flag) {
            viewOptions |= flag;
        }
        _UpdateViewMenu();
    }
}


//__________________________________________________________

void    _HYParameterTable::DoEnterConstraint (void)
{
    _String prStr = "Enter a constraint/command",
            newCons (stashCustomCommand);

    if (!EnterStringDialog (newCons,prStr, (Ptr)this)) {
        return;
    }

    stashCustomCommand = newCons;

    _HYDataPanel*parent = (_HYDataPanel*)RetrieveParentDataPanel ();

    _String * undoSnapshot = nil;

    if (parent) {
        undoSnapshot = parent->LFSnapshot();
    }

    _ExecutionList  ex      (newCons);
    ex.Execute              ();
    _HYTable*               table = (_HYTable*)    GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);
    ToggleConstrainedView(HY_PARAMETER_TABLE_VIEW_CONSTRAINED);

    if (lfID>=0) {
        //_LikelihoodFunction* lf = (_LikelihoodFunction*)likeFuncList (lfID);
        //lf->RescanAllVariables();
        VerifyGlobalVariables();
    }

    avViewOptions = GetViewOptions();

    if (avViewOptions&HY_PARAMETER_TABLE_VIEW_GLOBAL) {
        viewOptions |= HY_PARAMETER_TABLE_VIEW_GLOBAL;
    }
    if (avViewOptions&HY_PARAMETER_TABLE_VIEW_CATEGORY) {
        viewOptions |= HY_PARAMETER_TABLE_VIEW_CATEGORY;
    }

    _UpdateViewMenu();



    ConstructTheTable ();
    StretchColumnToFit (3);
    table->_MarkForUpdate   ();
    if (undoSnapshot) {
        _UpdateUndoMenu (undoSnapshot, &undoCustomConstraints);
        DeleteObject (undoSnapshot);
    } else {
        undoCommands.Clear();
        undoDescriptions.Clear();
        _UpdateUndoMenu (nil,nil);
    }
}


//__________________________________________________________

void    _HYParameterTable::DoSave (void)
{

    _HYTable*       table = (_HYTable*)    GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0),
                    *       head  = (_HYTable*)    GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW-1,0);

    _SimpleList columns,
                sel,
                rows;
    columns << 1;
    columns << 2;
    columns << 3;
    table->GetSelection (sel);

    char    resp = 3;
    if (sel.lLength) {
        _String pr ("Would you like to save only the selected cells? (Click \"No\" to save the entire table).");
        resp = YesNoCancelPrompt (pr);
    }

    if (resp == 1)
        for (long k = 0; k < sel.lLength; k+=4) {
            rows << sel.lData[k]/4;
        }
    else if (resp == 3)
        for (long k = 0; k < table->verticalSpaces.lLength; k++) {
            rows << k;
        }
    else {
        return;
    }

    _List           menuOptions;
    _String         menuChoice,
                    filePath,
                    str1 ("Save table as"),
                    str2;

    _HYTable::GetTableFormats(menuOptions);

    str2 = *(_String*)likeFuncNamesList (lfID) & " Parameter Table";

    long menuSel = SaveFileWithPopUp (filePath, str1 ,str2 , empty, menuOptions);

    if (menuSel>=0) {
        FILE*   outFile = doFileOpen (filePath.sData,"w");
        if (!outFile) {
            menuChoice = filePath & " could not be opened for writing.";
            ProblemReport (menuChoice,(Ptr)this);
            return;
        }
        table->SaveTable (head,nil, menuSel,outFile,GetTitle(),columns, rows);
        fclose (outFile);
    }
}
//__________________________________________________________
void    _HYParameterTable::HandleCategories(void)
{
    if (lfID>=0) {
        _LikelihoodFunction* lf = (_LikelihoodFunction*)likeFuncList  (lfID);
        _SimpleList* cv         = &lf->GetCategoryVars ();
        _SimpleList allParts (lf->CountObjects(0),0,1);
        _Matrix*    marginals   = lf->ConstructCategoryMatrix (allParts, _hyphyLFConstructCategoryMatrixConditionals, true);
        if (cv->lLength) {
            _String     windowTitle = _String ("Category Display for ") & *(_String*)likeFuncNamesList (lfID),
                        tryMe (windowTitle);

            long  k = 2;
            while (FindWindowByName (tryMe)>=0) {
                tryMe = windowTitle & '_' & k++;
            }

            _List       colHeaders;
            for (k=1; k<=marginals->GetHDim(); k++) {
                windowTitle = _String ("Class ") & k;
                colHeaders && & windowTitle;
            }
            _List       varInfo;

            for (long vv = 0; vv < cv->lLength; vv++) {
                _List   anEntry;
                _CategoryVariable * cvar = (_CategoryVariable*) LocateVar (cv->lData[vv]);
                anEntry << cvar->GetName   ();
                anEntry << cvar->GetWeights();
                anEntry << cvar->GetValues ();
                varInfo && & anEntry;
            }

            marginals->Transpose();
            _HYDistributionChartWindow* dcw = new _HYDistributionChartWindow (tryMe,colHeaders,*marginals,varInfo,nil);
            checkPointer (dcw);
            dcw->BringToFront();
        } else {
            //_String status = _String ("Likelihood Function ") & *(_String*)likeFuncNamesList (lfID) & " does not include any category variables. Opening a list of site-by-site log-likelihoods.";
            //ProblemReport (status, (Ptr)this);

            _String     windowTitle = _String ("Site-by-site likelihoods for ") & *(_String*)likeFuncNamesList (lfID),
                        tryMe (windowTitle);

            long  k = 2;
            while (FindWindowByName (tryMe)>=0) {
                tryMe = windowTitle & '_' & k++;
            }

            _SimpleList allParts (lf->CountObjects(0),0,1);
            _Matrix*    marginals           = lf->ConstructCategoryMatrix (allParts, _hyphyLFConstructCategoryMatrixSiteProbabilities, true);
            _List       colHeaders;
            colHeaders.AppendNewInstance (new _String ("Log-likelihood"));

            marginals->Transpose();
            _HYChartWindow* ncw = (_HYChartWindow*)checkPointer(new _HYChartWindow (tryMe,colHeaders,*marginals));
            ncw->BringToFront();
        }
        DeleteObject (marginals);
    }
}
//__________________________________________________________
void    _HYParameterTable::HandleProfilePlot (void)
{
    if (lfID>=0) {
        _LikelihoodFunction* lf = (_LikelihoodFunction*)likeFuncList  (lfID);
        _HYTable*       table   = (_HYTable*)           GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);

        // first build the list of independent parameters from the
        // likelihood function and sort it

        _SimpleList     rows;

        table->GetRowSelection (rows);

        if (rows.lLength == 1) {
            long  f = LocateVarByName (*(_String*)table->GetCellData(1,rows.lData[0]));

            if (f>=0) {
#ifndef USE_AVL_NAMES
                f = variableReindex.lData[f];
#else
                f = variableNames.GetXtra (f);
#endif

                _Variable* varObj  = LocateVar (f);
                _String*   varName = varObj->GetName();

                f = lf->GetIndependentVars().Find(f);

                if (f>=0) {
                    long       stepCount  = 10;
                    bool       computeNormalApproximation = true;

                    _Parameter MLE        = lf->GetIthIndependent (f),
                               leftBound  = varObj->GetLowerBound(),
                               rightBound = varObj->GetUpperBound(),
                               h          = (1.e-7)*MLE,
                               step,
                               MLEVal,
                               sigma     = 1.;

                    if (!HandlePLDialog (leftBound,rightBound,stepCount,computeNormalApproximation,MLE,(Ptr)this)) {
                        return;
                    }


                    if (h<1e-7) {
                        h = 1e-7;
                    }

                    step      = (rightBound-leftBound)/(stepCount-1);
                    lf->PrepareToCompute();
                    MLEVal    = lf->Compute();

                    if (computeNormalApproximation) {
                        lf->SetIthIndependent (f,MLE+h);
                        sigma = lf->Compute () - (MLEVal + MLEVal);
                        lf->SetIthIndependent (f,MLE-h);
                        sigma = (sigma+lf->Compute())/(h*h);
                    }

                    h = leftBound;

                    _Matrix      profilePlot (stepCount+1,2, false, true);

                    profilePlot.Store (0,0,h);
                    lf->SetIthIndependent (f,h);
                    profilePlot.Store (0,1,lf->Compute()-MLEVal);

                    for (long k = 1; k <= stepCount; k++) {
                        h += step;
                        if ((h>MLE)&&(h-step<MLE)) {
                            profilePlot.Store (k,0,MLE);
                            profilePlot.Store (k++,1,0);

                        }
                        lf->SetIthIndependent (f,h);
                        profilePlot.Store (k,0,h);
                        profilePlot.Store (k,1,lf->Compute()-MLEVal);
                    }

                    _List        labels;
                    _String reportName = *varName;
                    labels && & reportName;
                    reportName = "Log[L]";
                    labels && & reportName;

                    reportName = _String ("Likelihood Profile Plot for ") &
                                 *varName &" in " & *(_String*)likeFuncNamesList (lfID);

                    _HYChartWindow *reportChart = (_HYChartWindow*) FindWindowByNameAndOpen (reportName);
                    if (reportChart) {
                        reportChart->SetTable(labels, profilePlot);
                    } else {
                        reportChart = new _HYChartWindow (reportName, labels, profilePlot, nil);
                    }

                    reportChart->SetChartType("Line Plot",*varName,"Log[L]", false);
                    if (computeNormalApproximation) {
                        reportName = _String("0.5*(")&sigma&")*(_x_-"&MLE&")^2";
                    } else {
                        reportName = empty;
                    }

                    reportChart->SetLabels   (*varName,empty,"L_Max-L",HY_CHART_LEGEND_NONE,reportName,-1,-1,-1);

                    lf->SetIthIndependent (f,MLE);
                    lf->DoneComputing();

                    reportChart->BringToFront();
                    return;
                }
            }
        }

        _String             errStr   ("Profile Plot needs exactly one selected independent parameter.");
        ProblemReport       (errStr, (Ptr)this);

    }
}

long        _hypt_lastSelectParamChoice   = 0;
bool        _hypt_lastSelectParamChoiceCS = false;
_String     _hypt_lastSelectParamChoicePT;

//__________________________________________________________
void    _HYParameterTable::HandleSelectParameters (void)
{

    _String prompt  ("Select Parameters Matching RegExp:"),
            cprompt ("Case Sensitive"),
            lEntry;

    _List   choices;

    lEntry = "Parameter Name";
    choices && & lEntry;
    lEntry = "Parameter Value";
    choices && & lEntry;
    lEntry = "Parameter Constraint";
    choices && & lEntry;

    long    msel;

    if (!EnterStringDialogWithPulldown (_hypt_lastSelectParamChoicePT,prompt,cprompt,msel,choices,nil,_hypt_lastSelectParamChoiceCS,_hypt_lastSelectParamChoice,(Ptr)this)) {
        return;
    }

    _HYTable*       table   = (_HYTable*)           GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);
    _SimpleList     newRowSelection;
    _List           stringsToMatch;

    msel++; // row here

    for (long k=0; k<table->verticalSpaces.lLength; k++) {
        stringsToMatch << table->GetCellData (msel,k);
    }

    int errCode;
    Ptr rxp = PrepRegExp (&_hypt_lastSelectParamChoicePT, errCode,_hypt_lastSelectParamChoiceCS);
    if (!rxp) {
        lEntry = GetRegExpError (errCode);
        ProblemReport (lEntry);
        return;
    }

    for (long k=0; k<stringsToMatch.lLength; k++) {
        _SimpleList matchedPairs;
        ((_String*)stringsToMatch(k))->RegExpMatch (rxp,matchedPairs);
        if (matchedPairs.lLength) {
            newRowSelection << k;
        }
    }
    table->SetRowSelection (newRowSelection);
    FlushRegExp (rxp);

}

//__________________________________________________________
void    _HYParameterTable::HandleOpenInChart (void)
{
    _String      aString ("Value");

    _HYTable*       table = (_HYTable*)  GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);
    _SimpleList     rowSel;
    table->GetRowSelection (rowSel);

    if (rowSel.lLength == 0) {
        aString = "No parameters have been selected";
        ProblemReport (aString,(Ptr)this);
        return;
    }

    _List        colHeaders;
    _Matrix      res (rowSel.lLength, 1, false, true);
    _String      parmNames (128L,true);
    parmNames << "Name";

    colHeaders && & aString;
    for (long k=0; k<rowSel.lLength; k++) {
        res.theData[k] = ((_String*)table->GetCellData (2,rowSel.lData[k]))->toNum();
        parmNames << ';';
        parmNames << (_String*)table->GetCellData (1,rowSel.lData[k]);
    }

    parmNames.Finalize();
    colHeaders && & parmNames;

    aString = _String("Selected Parameter Values for ") &
              *(_String*) (likeFuncNamesList.lData[lfID]);

    long f = FindWindowByName (aString);
    _HYChartWindow * nc;
    if (f>=0) {
        nc = (_HYChartWindow*)windowObjectRefs (f);
        nc->SetTable (colHeaders,res);
    } else {
        nc = new _HYChartWindow (aString, colHeaders, res, nil);
        checkPointer (nc);
    }
    nc->SetChartType ("Bar Chart","Index","Value", true);
    nc->BringToFront();
}

//__________________________________________________________
void    _HYParameterTable::HandleVarianceEstimates (void)
{
    if (lfID>=0) {
        char        options;
        _Parameter  cLevel,
                    cLevelR;

        _LikelihoodFunction* lf = (_LikelihoodFunction*)likeFuncList(lfID);
        _HYTable*       table = (_HYTable*)  GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);

        // first build the list of independent parameters from the
        // likelihood function and sort it

        _List           independents;
        _SimpleList     indIndices,
                        * lfVars = &lf->GetIndependentVars(),
                          filteredList,
                          rows;

        _String         errStr ("Please select some independent variables to work with.");

        table->GetRowSelection (rows);

        if (rows.lLength == 0) {
            ProblemReport (errStr,(Ptr)this);
            return;
        }

        for (long k=0; k<lfVars->lLength; k++) {
            independents << LocateVar(lfVars->lData[k])->GetName();
            indIndices   << k;
        }

        SortLists (&independents, &indIndices);

        for (long k2=0; k2<rows.lLength; k2++) {
            long    f = independents.BinaryFind (table->GetCellData(1,rows.lData[k2]));
            if (f>=0) {
                filteredList << indIndices.lData[f];
            }
        }

        if (filteredList.lLength == 0) {
            ProblemReport (errStr, (Ptr)this);
            return;
        }

        if (HandleCIDialog (options, cLevel, cLevelR, (Ptr)this)) {
            if (options < 6) {
                _Matrix   * ci    = nil,
                            * covMx = nil;

                _String reportName,
                        rowLabels (128L, true);

                if (options<2) {
                    stashParameter (covariancePrecision, 1.+options,true);
                } else if (options == 3) {
                    stashParameter (covariancePrecision, cLevel,true);
                } else {
                    stashParameter (covariancePrecision, -cLevel,true);
                }


                covMx = (_Matrix*)lf->CovarianceMatrix (&filteredList);

                independents.Clear();

                for (long k2=0; k2<filteredList.lLength; k2++) {
                    independents << LocateVar(lfVars->lData[filteredList.lData[k2]])->GetName();
                }


                rowLabels << "Parameter";

                for (long k=0; k<independents.lLength; k++) {
                    rowLabels << ';';
                    rowLabels << (_String*)independents(k);
                }

                rowLabels.Finalize();
                independents && & rowLabels;

                if (options < 2) {
                    reportName = _String ("Covariance Matrix For ") & *(_String*)likeFuncNamesList(lfID);

                    _HYChartWindow * reportChart = (_HYChartWindow*) FindWindowByNameAndOpen (reportName);
                    if (reportChart) {
                        reportChart->SetTable(independents, *covMx);
                    } else {
                        reportChart = new _HYChartWindow (reportName, independents, *covMx, nil);
                    }

                    reportChart->BringToFront();
                }

                independents.Clear();
                reportName = "Left Bound";
                independents && & reportName;
                reportName = "MLE";
                independents && & reportName;
                reportName = "Right Bound";
                independents && & reportName;
                independents && & rowLabels;


                if (options<2) {
                    reportName = _String ("Asymptotic Normal 95% CI For ") & *(_String*)likeFuncNamesList(lfID);
                    ci = new _Matrix (covMx->GetHDim(),3,false,true);
                    checkPointer (ci);

                    for (long k=0; k<filteredList.lLength; k++) {
                        _Variable* aParam = LocateVar(lfVars->lData[filteredList.lData[k]]);
                        _Parameter val = aParam->Value(),
                                   lb  = aParam->GetLowerBound(),
                                   ub  = aParam->GetUpperBound(),
                                   size= 1.96*sqrt ((*covMx)(k,k));

                        ci->Store (k,0,val-lb>size?val-size:lb);
                        ci->Store (k,1,val);
                        ci->Store (k,2,ub-val>size?val+size:ub);
                    }
                } else {
                    if (options == 3) {
                        reportName = _String ("Likelihood Profile CI For ") & *(_String*)likeFuncNamesList(lfID) & " [chi2 level " & cLevel & ']';
                    } else {
                        reportName = _String ("Likelihood Profile CI For ") & *(_String*)likeFuncNamesList(lfID) & " [support  " & cLevel & ']';
                    }

                    ci = covMx;
                }

                _HYChartWindow *reportChart = (_HYChartWindow*) FindWindowByNameAndOpen (reportName);
                if (reportChart) {
                    reportChart->SetTable(independents, *ci);
                } else {
                    reportChart = new _HYChartWindow (reportName, independents, *ci, nil);
                }

                reportChart->SetChartType("Line Plot","Index","Left Bound;MLE;Right Bound", false);
                reportChart->SetLabels   ("Parameter Index",empty,"Values",HY_CHART_LEGEND_TOP_RIGHT,empty,-1,-1,-1);
                reportChart->BringToFront();

                stashParameter (covariancePrecision, 0,false);
                if (ci!=covMx) {
                    DeleteObject (ci);
                }

                if (covMx) {
                    DeleteObject (covMx);
                }
            } else {
                /* try to read the plug-in - should be in ChartAddIns/SIR/sampler.bf */
                errStr = libDirectory & "ChartAddIns" & GetPlatformDirectoryChar() &"Samplers" & GetPlatformDirectoryChar() & (options==6?"sir.bf":"lhc.bf");

                FILE *    procFile = doFileOpen (errStr.sData,"rb");
                if (!procFile) {
                    errStr = "Your distribution is missing the $HYPHY/ChartAddIns/SIR/sampler.bf file needed for this function";
                    ProblemReport (errStr, (Ptr)this);
                } else {
                    _String   procCode (procFile),
                              prefix   (128L, true);

                    fclose (procFile);

                    prefix << "SAMPLE_N = ";
                    prefix << _String(cLevel);
                    prefix << ";\nSAMPLE_M = ";
                    prefix << _String(cLevelR);
                    prefix << ";\nLF_NAME = \"";
                    prefix << (_String*)likeFuncNamesList (lfID);
                    prefix << "\";\n";
                    prefix << covarianceParameterList;
                    prefix << " = {};\n";

                    for (long k2=0; k2<filteredList.lLength; k2++) {
                        prefix << covarianceParameterList;
                        prefix << "[\"";
                        prefix << LocateVar(lfVars->lData[filteredList.lData[k2]])->GetName();
                        prefix << "\"]=1;\n";
                    }

                    prefix.Finalize();
                    procCode = prefix & procCode;

                    PushFilePath (errStr);
                    _ExecutionList sir (procCode);
                    PopFilePath ();
                    sir.Execute();
                }
            }


        }
    }
}

//__________________________________________________________

void    _HYParameterTable::DoOpenTreeWindow (void)
{
    _HYTable*       table = (_HYTable*)    GetCellObject (HY_PARAMETER_TABLE_TABLE_ROW,0);
    _HYButtonBar*   bb    = (_HYButtonBar*)GetCellObject (HY_PARAMETER_TABLE_BUTTON_ROW,0);
    _SimpleList     rows;
    _String         nodeName,
                    *cellName;
    _List           selectNodes;
    long            k,
                    f;
    long            treeVar = 0;

    table->GetRowSelection (rows);

    bb->_UnpushButton();
    for (k=0; k<rows.lLength; k++) {
        cellName = (_String*)table->GetCellData (1,rows.lData[k]);
        f = cellName->FindBackwards ('.',0,-1);
        if (f>=0) {
            nodeName = cellName->Cut (0,f-1);
            f = selectNodes.BinaryInsert (&nodeName);
        } else if (((_Variable*)RetrieveIthRowVar(rows.lData[k]))->ObjectClass()==TREE) {
            treeVar = k+1;
        }
    }

    _String treeTitle;
    if (selectNodes.lLength||treeVar) {
        if (selectNodes.lLength) {
            cellName = (_String*)selectNodes(0);
            f = cellName->Find ('.');
            if (f>=0) {
                nodeName = cellName->Cut (0,f-1);
            }
        } else {
            if (!treeVar) {
                return;
            }
            nodeName = *((_Variable*)RetrieveIthRowVar(rows.lData[treeVar-1]))->GetName();
        }
        treeTitle = _String ("Tree ") & nodeName;
    } else {
        return;
    }

    _HYTreePanel* thisPanel = (_HYTreePanel*)FindWindowByNameAndOpen (treeTitle);
    if (!thisPanel) {
        thisPanel = new _HYTreePanel (nodeName,nodeName);
        checkPointer (thisPanel);
        thisPanel->_Zoom (true);
        thisPanel->BringToFront();
    } else {
        if (selectNodes.lLength == 0) {
            _String prompt ("A window for this tree has already been opened. Spawn another instance of this tree plot?");
            if (ProceedPrompt (prompt,(Ptr)thisPanel)) {
                thisPanel = new _HYTreePanel (nodeName,nodeName);
                checkPointer (thisPanel);
                thisPanel->_Zoom (true);
                thisPanel->BringToFront();
            } else {
                return;
            }
        }
    }
    thisPanel->SelectRangeAndScroll (selectNodes);
}

//__________________________________________________________

void    _HYParameterTable::OpenCategoryPlotWindow (_CategoryVariable* thisCV)
{
    _String windowTitle;
    windowTitle = _String ("Density Plot for ")&*thisCV->GetName();
    if (!FindWindowByNameAndOpen (windowTitle)) {
        long            width  = 500,
                        height = 400,
                        i = thisCV->GetNumberOfIntervals(),
                        leftMargin = 10;

        _Matrix*        intervalEnds = thisCV->GetIntervalEnds();
        _HYColor        boxColor = {255,255,255};
        _HYRect         set;

        _Parameter      xMin = thisCV->GetMinX(),
                        xMax = thisCV->GetMaxX(),
                        xTick,
                        step,
                        vstep;

        char            buffer[255];

        if (xMax>1e10) {
            xMax = thisCV->GetIntervalValue(i-1)*1.5;
        }
        if (xMin<1e-25) {
            xMin = 0;
        }

        step = (xMax-xMin)/(width-20);

        _Parameter      *valueArray = new _Parameter [width-20],
        max = 0.0,
        tracer = xMin+step,
        x,
        y = 0,
        z = 0;
        long            counter,
                        tW = 0,
                        tH = 0;

        checkPointer    (valueArray);

        for (counter = 1; counter < width-20; counter++) {
            _Constant c (tracer);
            _x_->SetValue (&c);
            valueArray[counter] = thisCV->GetDensity().Compute()->Value();
            if (valueArray[counter]>max) {
                max = valueArray[counter];
            }
            tracer+=step;
        }

        vstep = (height-20)/max;

        _HYFont         labelFont;
        labelFont.face  = "Times";
        labelFont.size  = 10;
        labelFont.style = HY_FONT_PLAIN;

        _List           classLabels,
                        cellPlacement;

        set.top         = 0;
        cellPlacement.AppendNewInstance (new _String ("Class"));
        cellPlacement.AppendNewInstance (new _String ("Rate"));
        classLabels     && & cellPlacement;

        cellPlacement.Clear();

        for (counter = 0; counter < i; counter++) {
            _List       *thisRow = new _List;
            checkPointer (thisRow);

            x = thisCV->GetIntervalValue(counter);
            snprintf (buffer, sizeof(buffer),"%.6g",x);

            thisRow->AppendNewInstance(new _String (counter));
            thisRow->AppendNewInstance(new _String (buffer));

            tracer = thisCV->GetIntervalWeight(counter);
            y += x*tracer;
            z += x*x*tracer;

            classLabels << thisRow;
        }

        _HYRect     margins = {3,3,3,3,1};

        tW = ComputeTableCellPlacement (classLabels, cellPlacement, margins, labelFont);

        /*snprintf (buffer, sizeof(buffer),"E[X] = %8g",y);
        set.left = GetVisibleStringWidth (*((_String*)classLabels (counter)),labelFont);
        if (set.left > set.top)
            set.top = set.left;
        classLabels.InsertElement (new _String (buffer),0,false);
        snprintf (buffer, sizeof(buffer),"Var[X] = %8g",z-y*y);
        set.left = GetVisibleStringWidth (*((_String*)classLabels (counter)),labelFont);
        if (set.left > set.top)
            set.top = set.left;
        classLabels.InsertElement (new _String (buffer),1,false);*/

        tH = (i+3)*(labelFont.size+6)+10;

        for (counter = width-21; counter>=1; counter--)
            if (height-15-valueArray[counter]*vstep<tH) {
                break;
            }

        counter = width-20-counter;
        if (counter>tW) {
            counter = tW;
        }


        set.right  = set.left = leftMargin;
        set.top    = 10;
        set.bottom = height - 15;
        set.width = 2;

        _SimpleList     tickInfo;
        _List           labelInfo;

        set.width = ComputeHashMarkPlacement (set,max/1e50,max,xTick,tracer,tickInfo,labelInfo,0,labelFont);
        if (set.width+2>10) {
            leftMargin = set.width+2;
        }

        _HYPWindow*     graph = new _HYPWindow(windowTitle,height,width+tW-counter+leftMargin-10,32,false);

        graph->StartDraw();
        graph->_HYGraphicPane::SetFont  (labelFont);
        set.left = 0;
        set.top = 0;
        set.bottom = height;
        set.right = width;
        graph->SetColor (boxColor);
        graph->FillRect (set);

        set.bottom = height-17;
        set.width  = 1;
        graph->SetColor ((_HYColor) {
            0,0,0
        });

        for (counter = 0; counter < i-1; counter++) {
            set.left = set.right = leftMargin+2+(intervalEnds->theData[counter]-xMin)/step;
            set.top = height-15-valueArray[set.left-leftMargin-1]*vstep;
            graph->DrawHatchedLine (set);
        }

        set.top = height - 17;
        set.bottom = height - 13;
        boxColor.R = boxColor.B = 0;
        boxColor.G = 100;
        graph->SetColor (boxColor);
        set.width  = 2;

        for (counter = 0; counter < i; counter++) {
            set.left = set.right = leftMargin+2+(thisCV->GetIntervalValue(counter)-xMin)/step;
            set.top = height-15-valueArray[set.left-leftMargin-1]*vstep;
            graph->DrawLine (set);
        }

        set.left = leftMargin+2;
        set.top = height-15-valueArray[1]*vstep;
        set.width = 2;

        boxColor.B = 50;
        boxColor.G = 50;
        boxColor.R = 255;
        graph->SetColor (boxColor);

        set.right = leftMargin+2;

        for (counter = 2; counter < width-21; counter++) {
            set.right++;
            set.bottom = height-15-valueArray[counter]*vstep;
            //if ((abs(set.top-set.bottom)>4)||(counter==width-22))
            {
                graph->DrawLine (set);
                set.left = set.right;
                set.top = set.bottom;
            }
        }
        set.right++;
        set.bottom = height-15-valueArray[counter]*vstep;
        graph->DrawLine (set);

        delete       valueArray;

        /*set.top = 20;
        set.left = graph->w-5-lW;


        for (counter = 0; counter < i; counter++)
        {
            graph->DisplayText (*(_String*)classLabels(counter),set.top, set.left,true);
            set.top += labelFont.size*1.3;
        }*/
        set.right  = tW;
        set.bottom = tH-10-2*(labelFont.size+6);
        set.left   = graph->w-5-tW;
        set.top    = 5;

        graph->SetColor ((_HYColor) {
            0,0,0
        });
        graph->DrawTable (classLabels,cellPlacement,HY_ALIGN_LEFT, margins,set);
        snprintf (buffer, sizeof(buffer),"E[X]   = %8.3g", y);
        set.top += set.bottom+3+labelFont.size;
        graph->DisplayText (buffer,set.top,set.left+margins.left, true);
        snprintf (buffer, sizeof(buffer),"Var[X] = %8.3g", z-y*y);
        graph->DisplayText (buffer,set.top+6+labelFont.size,set.left+margins.left, true);
        set.right  = set.left = leftMargin;
        set.top    = 10;
        set.bottom = height - 15;
        set.width  = 2;

        graph->DrawHashes (set, tickInfo, labelInfo, 5, 1);
        graph->DrawLine   (set);

        set.right = width - 10;
        set.top = set.bottom;

        ComputeHashMarkPlacement (set,xMin,xMax,xTick,tracer,tickInfo,labelInfo,2,labelFont);
        graph->DrawHashes        (set, tickInfo, labelInfo, 5, 1);
        graph->DrawLine (set);
        graph->EndDraw();
        graph->SetWindowRectangle (0,0,height,graph->w);
        graph->BringToFront();
    }
}

//__________________________________________________________

void    _HYParameterTable::OptimizeLikelihoodFunction (void)
{
    if ((lfID>=0)&&(!isInOptimize)) {
        terminateExecution = false;
        _LikelihoodFunction *lf = (_LikelihoodFunction*)likeFuncList (lfID);
        _Matrix * res;
        if (lf->GetIndependentVars().lLength<10) {
            _Parameter prec,
                       usePrev;
            checkParameter (optimizationPrecision,prec,0.001);
            if (prec>0.001) {
                stashParameter (optimizationPrecision,0.001,true);
            }

            checkParameter (useLastResults,usePrev,0.0);
            if (CheckEqual(usePrev,0.0)) {
                setParameter (useLastResults, 1.0);
            }

            ToggleAnalysisMenu (true);
            StartBarTimer ();
            res = lf->Optimize();
            StopBarTimer  ();
            ToggleAnalysisMenu (false);
            terminateExecution = false;
            if (prec>0.001) {
                stashParameter (optimizationPrecision,0.001,false);
            }

            setParameter (useLastResults, usePrev);
        } else {
            ToggleAnalysisMenu (true);
            StartBarTimer ();
            res = lf->Optimize();
            StopBarTimer  ();
            ToggleAnalysisMenu (false);
            terminateExecution = false;
        }

        DeleteObject (res);
        _HYTable*table = (_HYTable*)GetCellObject(HY_PARAMETER_TABLE_TABLE_ROW,0);
        for (long k=1; k<table->verticalSpaces.lLength; k++) {
            UpdateKthRow (k,false);
        }
        UpdateLogLikelihood ();
        StretchColumnToFit (2);

        ReportAnalysisAsFinished (empty);
    }
}

//__________________________________________________________

void    _HYParameterTable::StretchColumnToFit (long index)
{
    _HYTable*table = (_HYTable*)GetCellObject(HY_PARAMETER_TABLE_TABLE_ROW,0);
    _HYTable*table2= (_HYTable*)GetCellObject(HY_PARAMETER_TABLE_TABLE_ROW-1,0);
    table->AutoFitColumn (index,true,true);
    if (index<3) {
        table2->SetColumnSpacing(index,table->GetColumnSpacing(index)-table2->GetColumnSpacing(index),true);
    } else {
        table2->SetColumnSpacing(index,table->GetColumnSpacing(index)+HY_SCROLLER_WIDTH-table2->GetColumnSpacing(index),true);
    }
    for (; index<4; index++) {
        table->_MarkColumnForUpdate (index);
        table2->_MarkColumnForUpdate (index);
    }
    table->SetVisibleSize(table->rel);
#ifdef __HYPHY_GTK__
    UpdateComponentInfo ();
#else
    dim = MinMaxWindowDimensions();
#endif
}

//__________________________________________________________

void    _HYParameterTable::UndoCommand (void)
{
    long k;
    if ((k = undoCommands.lLength)) {
        k--;
        TakeLFSnapshot();
        _ExecutionList ex (*(_String*)undoCommands(k));
        ex.Execute();
        terminateExecution = false;
        undoCommands.Delete(k);
        undoDescriptions.Delete (k);
        _HYTable*table = (_HYTable*)GetCellObject(HY_PARAMETER_TABLE_TABLE_ROW,0);
        VerifyGlobalVariables();
        RefreshTheTable();
        table->_MarkContentsForUpdate();
        ToggleConstrainedView ();
        ToggleConstrainedView (HY_PARAMETER_TABLE_VIEW_GLOBAL);
        StretchColumnToFit (3);
        UpdateSelectionDependentButtons ();
        DoneLFSnapshot();
    }
}

//__________________________________________________________

_HYGuiObject*   _HYParameterTable::RetrieveParentDataPanel (void)
{

    for (long i=0; i<windowObjects.lLength; i++) {
        _HYWindow* thisWindow = (_HYWindow*)windowObjectRefs(i);
        if (thisWindow->WindowKind () == HY_WINDOW_KIND_DATAPANEL)
            if (((_HYDataPanel*)thisWindow)->GetLFID() == lfID) {
                return (_HYGuiObject*)windowObjectRefs.lData[i];
            }
    }
    return nil;
}

//__________________________________________________________
// HYBootstrapWindow
//__________________________________________________________


_String     gbGoTooltip     ("Start Bootstrapping"),
            gbStopTooltip    ("Stop  Bootstrapping");

_HYBootstrapWindow::_HYBootstrapWindow (_String name, _Parameter nLRT):_HYTWindow (name, true)
{

    done            = false;
    requestQuit     = false;

    nullLRT         = nLRT;
    iterateCount    = 0;
    iterateCountP   = 0;

    _HYRect         canvasSettings = {25,50,25,50,HY_COMPONENT_BORDER_T};

    _HYButtonBar*   bb      = new _HYButtonBar (canvasSettings, GetOSWindowData());

    canvasSettings.left = canvasSettings.right = 250;

    _HYLabel*       l3      = new _HYLabel (canvasSettings, GetOSWindowData());

    canvasSettings.left = canvasSettings.right = 300;
    canvasSettings.top  = canvasSettings.bottom = 21;
    canvasSettings.width = HY_COMPONENT_BORDER_B;

    _HYTable*   table2      = new _HYTable (canvasSettings,GetOSWindowData(),1,2,100,20,HY_TABLE_STATIC_TEXT);

    canvasSettings.top      = 60;
    canvasSettings.bottom   = 100000;
    canvasSettings.width    = HY_COMPONENT_V_SCROLL;

    _HYTable*   table       = new _HYTable (canvasSettings,GetOSWindowData(),5,2,100,18,HY_TABLE_STATIC_TEXT);

    table->selectionType = HY_TABLE_SEL_ROWS|HY_TABLE_DONT_SIZE;
    table->SetColumnSpacing (0,-60,false);
    table->SetColumnSpacing (1,160,false);
    table2->selectionType = HY_TABLE_SEL_ROWS|HY_TABLE_DONT_SIZE;
    table2->SetColumnSpacing (0,-60,false);
    table2->SetColumnSpacing (1,160,false);
    _String cellValue ("#");
    table2->SetCellData (&cellValue,0,0,HY_TABLE_BOLD|HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_CANTSELECT,true);
    cellValue = "Simulated LR";
    table2->SetCellData (&cellValue,0,1,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_CANTSELECT|HY_TABLE_BOLD,true);


    table->SetMessageRecipient (this);
    bb->SetMessageRecipient (this);

    AddObject (l3);     // 0
    AddObject (table);  // 1
    AddObject (table2); // 2
    AddObject (bb);     // 3

    SetTableDimensions (3,2);
    SetCell   (1,0,table2);
    SetCell   (1,1,table2);
    SetCell   (2,0,table);
    SetCell   (2,1,table);
    SetCell   (0,0,bb);
    SetCell   (0,1,l3);

    _HYFont  labelFont;
    bb->SetBackColor (labelBColor);
    l3->SetBackColor (labelBColor);
    l3->SetForeColor (labelFColor);

#ifdef __MAC__
    labelFont.face  = "Geneva";
    labelFont.size  = 10;
#endif

#ifdef __WINDOZE__
    labelFont.face  = "MS Sans Serif";
    labelFont.size  = 12;
#endif

#ifdef __HYPHY_GTK__
    labelFont.face  = _HY_SANS_FONT;
    labelFont.size  = 12;
#endif

    labelFont.style = HY_FONT_PLAIN;
    table->SetFont(labelFont);
    table2->SetFont(labelFont);

#ifdef __MAC__
    labelFont.size  = 12;
#endif

#ifdef __WINDOZE__
    labelFont.size  = 14;
#endif

#ifdef __HYPHY_GTK__
    labelFont.size  = 12;
#endif


    l3->SetFont (labelFont);

    l3->SetText (_String("P-Value: 0"));
    l3->SetShadow(true);
    l3->SetAlignFlags (HY_ALIGN_LEFT);
    bb->SetAlignFlags (HY_ALIGN_LEFT);

    bb->AddButton (ProcureIconResource(HY_PARAMETER_BOOTSTRAP_BUTTON_ICON_BASE),&gbStopTooltip );
    cellValue = "Show LRT Histogram";
    bb->AddButton (ProcureIconResource(HY_PARAMETER_BOOTSTRAP_BUTTON_ICON_BASE+1),&cellValue );

    bb->EnableButton(0,true);
    bb->EnableButton(1,false);
    bb->MarkAsPullDown (1,true);
    bb->SetButtonDim(16);

    // set up the table

    _HYRect  screenRect = GetScreenDimensions();
    SetWindowRectangle (0,0,screenRect.bottom-10,screenRect.right-20);
    SetPosition (5,40);

    DeleteObject (l3);
    DeleteObject (bb);
    DeleteObject (table);
    DeleteObject (table2);

}

//__________________________________________________________

bool    _HYBootstrapWindow::ProcessEvent (_HYEvent* e)
{
    bool eDone = false;
    _String firstArg;
    long    i,k,f;
    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }
        if (i==3) {
            firstArg        = e->EventCode().Cut (f+1,-1);
            k               = firstArg.toNum();
            _HYButtonBar*   bb = (_HYButtonBar*)GetObject (3);
            switch  (k) {
            case 0:
                requestQuit = true;
                break;

            case 1: { // set equal
                int h,v;
                _List menuOptions;

                firstArg = "Histogram [Default]";
                menuOptions && & firstArg;
                firstArg = "Histogram [Custom]";
                menuOptions && & firstArg;
                firstArg = "Sample CDF [Default]";
                menuOptions && & firstArg;
                firstArg = "Sample CDF [Custom]";
                menuOptions && & firstArg;

                bb->GetButtonLoc(1,h,v,true);
                firstArg = HandlePullDown (menuOptions,h,v,0);
                bb->_UnpushButton();

                v = menuOptions.Find (&firstArg);

                if (v>=0) {
                    if (v%2 == 0) {
                        OpenHistogramWindow (v>1,-1);
                    } else {
                        firstArg = "Number of bins:";
                        if (EnterStringDialog (histChartBins, firstArg, (Ptr)this)) {
                            OpenHistogramWindow (v>1,histChartBins.toNum());
                        }
                    }
                }

                break;
            }

            }
            bb->_UnpushButton();
            eDone = true;
        }
    }
    if (eDone) {
        DeleteObject (e);
        return true;
    }
    return _HYTWindow::ProcessEvent(e);
}


//__________________________________________________________

void    _HYBootstrapWindow::SetStopped(bool d)
{
    if (done!=d) {
        _HYButtonBar * bb = (_HYButtonBar*)GetObject (3);
        bb->EnableButton (0,!d);
        bb->EnableButton (1,d);
        done = d;
    }
}

//__________________________________________________________

bool    _HYBootstrapWindow::ConfirmClose (void)
{
    if (!done || requestQuit) {
        _String m ("This window contains an active bootstrap simulation. Please stop it before closing the window.");
        ProblemReport (m,(Ptr)this);
        return false;
    }
    return true;
}

//__________________________________________________________

void    _HYBootstrapWindow::AddIterate (_Parameter lrt)
{
    //if (!done)
    {
        _HYTable* table  = (_HYTable*)GetObject (1);
        _HYTable* tableh = (_HYTable*)GetObject (2);
        _HYLabel* l      = (_HYLabel*)GetObject (0);
        long      type;
        if (lrt>nullLRT) {
            type = HY_TABLE_STATIC_TEXT|HY_TABLE_BOLD;
            iterateCountP++;
        } else {
            type = HY_TABLE_STATIC_TEXT;
        }

        if (iterateCount == 0) {
            _String hValue = _String("Simlated LR (test:") & nullLRT & ')';
            tableh->SetCellData (&hValue,0,1,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_CANTSELECT|HY_TABLE_BOLD, true);
            tableh->_MarkCellForUpdate (1);
        }

        if (iterateCount>=table->verticalSpaces.lLength) {
            table->AddRow     (-1,18,type);
        }

        _String cellValue (iterateCount+1);
        table->SetCellData(&cellValue,iterateCount,0,type,true);
        cellValue = lrt;
        table->SetCellData(&cellValue,iterateCount,1,type,true);
        table->_MarkCellForUpdate (iterateCount*2);
        table->_MarkCellForUpdate (iterateCount*2+1);
        table->ScrollToRow        (iterateCount);
        iterateCount++;
        cellValue = _String("P-Value: ")&(iterateCountP/(_Parameter)iterateCount);
        l->SetText (cellValue);
        table->SetVisibleSize(table->rel);
#ifdef __HYPHY_GTK__
        UpdateComponentInfo ();
#else
        dim = MinMaxWindowDimensions();
#endif
        handleGUI(true);
    }
}

//__________________________________________________________

void    _HYBootstrapWindow::OpenHistogramWindow (bool cdf, long binCount)
{
    _String windowTitle;
    if (cdf) {
        windowTitle = _String ("CDF: ")& GetTitle();
    } else {
        windowTitle = _String ("Histogram: ")& GetTitle();
    }
    if (!FindWindowByNameAndOpen (windowTitle)) {
        _HYTable* table = (_HYTable*)GetObject (1);
        long    width  = 550,
                height = 450,
                k,
                maxBin = 0;

        _Parameter minValue = nullLRT,
                   maxValue = nullLRT,
                   step,
                   vstep,
                   sum = 0.,
                   sum2 = 0.;

        if (binCount < 0) {
            binCount = iterateCount/5;
            if (binCount < 5) {
                binCount = 5;
            }
        }

        if (binCount>iterateCount) {
            binCount = iterateCount;
        }

        _Matrix values (1,iterateCount,false,true);
        _Matrix bins   (1,binCount,false,true);

        for (k=0; k<iterateCount; k++) {
            _Parameter tt = values.theData[k] = ((_String*)table->GetCellData (1,k))->toNum();
            if (tt < minValue) {
                minValue = tt;
            } else if (tt>maxValue) {
                maxValue = tt;
            }
            sum += tt;
            sum2 += tt*tt;
        }
        step = (maxValue-minValue)/binCount;
        for (k=0; k<iterateCount; k++) {
            _Parameter tt = values.theData[k];
            maxBin = (tt-minValue)/step;
            if (maxBin==binCount) {
                maxBin--;
            }
            bins.theData[maxBin]+=1.;
        }

        if (cdf) {
            for (k=1; k< binCount; k++) {
                bins.theData[k] += bins.theData[k-1];
            }
            maxBin = iterateCount;
        } else {
            maxBin = 0;
            for (k=0; k<binCount; k++)
                if(bins.theData[k]>maxBin) {
                    maxBin = bins.theData[k];
                }
        }

        _HYColor boxColor = {200,200,200},
                 barColor = {200,50,50},
                 barColor2= {50,200,50};

        _HYRect  box =  GetScreenDimensions();

        width  = box.right-100;
        height = box.bottom-100;

        /*if (width>550)
            width = 550;
        height = box.bottom-3*HY_SCROLLER_WIDTH;
        if (height>450)
            height = 450;*/

        _HYPWindow*    graph = new _HYPWindow(windowTitle,height,width,32,false);
        _HYFont        graphFont;

        graphFont.face = "Times";
        graphFont.size = 12;
        graphFont.style = HY_FONT_PLAIN;


        graph->StartDraw();
        graph->_HYGraphicPane::SetFont  (graphFont);

        box.left = 0;
        box.top = 0;
        box.bottom = height;
        box.right = width;
        graph->SetColor (boxColor);
        graph->FillRect (box);
        box.left = 10;
        box.right = width - 10;
        box.top = box.bottom = height - 15;
        box.width = 2;
        boxColor.R = boxColor.G = boxColor.B = 0;
        graph->SetColor (boxColor);
        graph->DrawLine (box);
        box.right = box.left;
        box.top = 15;
        graph->DrawLine (box);

        //box.width = 1;

        step = (width-20.)/binCount;
        vstep = (height-30.)/maxBin;

        long mult = 1,
             nValue = 11 + (nullLRT-minValue)/(maxValue-minValue)* (width-20);

        graph->SetColor (boxColor);
        graph->DrawRect (box);

        if (maxBin>20) {
            mult = 10;
        }
        if (maxBin>200) {
            mult = 100;
        }
        if (maxBin>2000) {
            mult = 1000;
        }

        k = mult;
        boxColor.R = boxColor.G = boxColor.B = 255;
        box.left = 11;
        box.right = width-10;

        graph->SetColor (boxColor);
        while (k<=maxBin) {
            box.bottom = height-14-vstep*k;
            box.top = box.bottom - 1;
            graph->DrawRect (box);
            k += mult;
        }

        boxColor.R = boxColor.G = boxColor.B = 0;
        box.bottom = height-13;
        for (k=0; k<binCount; k++) {
            box.left  = 11+step*k;
            box.right = box.left+step+2;
            box.top = box.bottom - vstep*bins.theData[k];
            if (box.bottom>box.top) {
                if (box.right < nValue) {
                    graph->SetColor (barColor2);
                    graph->FillRect (box);
                } else {
                    if (box.left > nValue) {
                        graph->SetColor (barColor);
                        graph->FillRect (box);
                    } else {
                        box.right = nValue;
                        graph->SetColor (barColor2);
                        graph->FillRect (box);
                        box.left = nValue;
                        box.right = 13+step*(k+1);
                        graph->SetColor (barColor);
                        graph->FillRect (box);
                        box.left  = 11+step*k;
                    }
                }
                graph->SetColor (boxColor);
                graph->DrawRect (box);
            }
        }

        /*windowTitle = minValue;
        graph->DisplayText (windowTitle,height-2, 10, true);

        windowTitle = maxValue;
        graph->DisplayText (windowTitle,height-2, width-GetVisibleStringWidth(windowTitle)-10,true);*/

        _HYRect    lr;
        lr.top   = lr.bottom = height-15;
        lr.right = width-10;
        lr.left  = 10;

        boxColor.R = boxColor.G = boxColor.B = 0;
        graph->SetColor (boxColor);

        _List       labels;
        _SimpleList labelInfo;

        ComputeHashMarkPlacement (lr,minValue,maxValue,step,vstep,labelInfo, labels,5,graphFont);

        lr.width = 2;

        graph->DrawHashes (lr,labelInfo,labels,5,1);

        windowTitle = _String (maxBin) & "    P-Value: " & (_Parameter)iterateCountP/iterateCount;

        graph->DisplayText (windowTitle,13,5,true);

        windowTitle = _String ("Mean: ")& sum/iterateCount & (", Variance: ") & (sum2-sum*sum/iterateCount)/(iterateCount-1);
#ifdef __MAC__
        graph->DisplayText (windowTitle,13, width-GetVisibleStringWidth(windowTitle)-10,true);
#else
        graph->DisplayText (windowTitle,13, width-GetVisibleStringWidth(windowTitle, graphFont)-10,true);
#endif

        box.left = nValue-1;
        box.right = nValue+1;
        box.bottom = height-5;
        box.top = box.bottom+10;
        graph->SetColor (barColor2);
        graph->FillRect (box);

        graph->EndDraw();
        graph->SetWindowRectangle (0,0,height,width);
        CenterWindow(graph);
        graph->BringToFront();
    }
}

//__________________________________________________________

void    _HYBootstrapWindow::DoSave (void)
{

    _HYTable*       table = (_HYTable*)    GetObject (1),
                    *       head  = (_HYTable*)    GetObject (2);

    _List           menuOptions;
    _String         menuChoice,
                    filePath,
                    str1 ("Save bootstrap results as:"),
                    str2;

    _HYTable::GetTableFormats(menuOptions);

    str2 = GetTitle();

    long menuSel = SaveFileWithPopUp (filePath, str1 ,str2 , empty, menuOptions);

    if (menuSel>=0) {
        FILE*   outFile = doFileOpen (filePath.sData,"w");
        if (!outFile) {
            menuChoice = filePath & " could not be opened for writing.";
            ProblemReport (menuChoice,(Ptr)this);
            return;
        }
        table->SaveTable (head,nil, menuSel,outFile,GetTitle());
        fclose (outFile);
    }
}

//__________________________________________________________
// HYGeneralBootstrapWindow
//__________________________________________________________


_HYGeneralBootstrapWindow::_HYGeneralBootstrapWindow (_String name, long lID):_HYBootstrapWindow (name, 0.0)
{

    lfNullID        = -1;
    lfAltID         = -1;
    hasGoIcon       = false;

    _HYRect         canvasSettings = {25,50,25,50,HY_COMPONENT_BORDER_T};
    _HYLabel*       l1      = new _HYLabel (canvasSettings, GetOSWindowData());
    checkPointer    (l1);
    canvasSettings.width = HY_COMPONENT_NO_SCROLL;
    _HYLabel*       l2      = new _HYLabel (canvasSettings, GetOSWindowData());
    checkPointer    (l2);

    canvasSettings.left = canvasSettings.right = 250;
    canvasSettings.width = HY_COMPONENT_BORDER_T;

    _HYPullDown *   p1      = new _HYPullDown (canvasSettings, GetOSWindowData());
    checkPointer   (p1);
    canvasSettings.width = HY_COMPONENT_NO_SCROLL;
    _HYPullDown *   p2      = new _HYPullDown (canvasSettings, GetOSWindowData());
    checkPointer   (p2);


    p1->SetMessageRecipient (this);
    p2->SetMessageRecipient (this);

    AddObject (p1);     // 4
    AddObject (p2);     // 5
    AddObject (l1);     // 6
    AddObject (l2);     // 7

    _HYGuiObject* bsComp [4] = {GetObject(0),GetObject(1),GetObject(2),GetObject(3)};

    SetTableDimensions (5,2);

    SetCell   (0,0,bsComp[3]);
    SetCell   (0,1,bsComp[0]);

    SetCell   (3,0,l1);
    SetCell   (3,1,p1);

    SetCell   (4,0,l2);
    SetCell   (4,1,p2);

    SetCell   (1,0,bsComp[2]);
    SetCell   (1,1,bsComp[2]);
    SetCell   (2,0,bsComp[1]);
    SetCell   (2,1,bsComp[1]);

    _HYFont  labelFont;
    p1->SetBackColor (labelBColor);
    p2->SetBackColor (labelBColor);
    l1->SetBackColor (labelBColor);
    l2->SetBackColor (labelBColor);
    l1->SetForeColor (labelFColor);
    l2->SetForeColor (labelFColor);

    labelFont.face = "Geneva";
    labelFont.size = 12;
    labelFont.style = HY_FONT_PLAIN;

    l1->SetFont (labelFont);
    l2->SetFont (labelFont);

    l1->SetText (_String("Null"));
    l2->SetText (_String("Alt."));
    l1->SetAlignFlags (HY_ALIGN_LEFT);
    l2->SetAlignFlags (HY_ALIGN_LEFT);
    p1->SetAlignFlags (HY_ALIGN_LEFT);
    p2->SetAlignFlags (HY_ALIGN_LEFT);

    p1->AddMenuItem (none, -1);
    p2->AddMenuItem (none, -1);

    done = true;

    SetNullLF         (lID,true);
    BuildCompatibleLF ();
    UpdateBSButtons();

    _HYRect  screenRect = GetScreenDimensions();
    SetWindowRectangle (0,0,screenRect.bottom-10,screenRect.right-20);
    SetPosition (5,40);


    DeleteObject (l1);
    DeleteObject (l2);
    DeleteObject (p1);
    DeleteObject (p2);

}

//__________________________________________________________

void    _HYGeneralBootstrapWindow::UpdateBSButtons      (void)
{
    _HYButtonBar *bb =  (_HYButtonBar*) GetObject (3);

    if (done) {
        if (!hasGoIcon) {
            bb->ReplaceButton (0, ProcureIconResource(HY_PARAMETER_BOOTSTRAP_BUTTON_ICON_BASE+2), &gbGoTooltip);
            bb->MarkAsPullDown (0,true);
            hasGoIcon = true;
        }
        bb->EnableButton (0,(lfNullID>=0)&&(lfAltID>=0));
    } else {
        if (hasGoIcon) {
            bb->ReplaceButton  (0,ProcureIconResource(HY_PARAMETER_BOOTSTRAP_BUTTON_ICON_BASE), &gbStopTooltip);
            bb->MarkAsPullDown (0,false);
            hasGoIcon = false;
        }
    }
}

//__________________________________________________________

void    _HYGeneralBootstrapWindow::UpdateLRValue        (void)
{
    /*_HYTable *th =  (_HYTable*) GetObject (2);

    _String  cellValue ("Simulated LR");

    if ((lfNullID>=0)&&(lfAltID>=0))
    {
    }*/

    //TBA?

}

//__________________________________________________________

void    _HYGeneralBootstrapWindow::SetNullLF        (long lID, bool checkFlag)
{
    _HYPullDown *p1 =  (_HYPullDown*) GetObject (4),
                 *p2    =  (_HYPullDown*) GetObject (5);

    long        f;

    if (checkFlag) {
        BuildAvailableLF();
        f = lfID0.Find (lID);
        p1->ChangeSelection (f,false);
    }

    if ((lID>=0)&&(lfNullID!=lID)) {
        lfNullID = lID;

        BuildCompatibleLF ();

        f = lfIDA.Find (lfAltID);
        if (f>=0) {
            p2->ChangeSelection (f,false);
        } else {
            lfAltID = -1;
            p2->ChangeSelection (0,false);
        }
        p2->EnableMenu (lfNullID>=0);
        UpdateBSButtons ();
        return;
    } else if (lID < 0) {
        lfNullID = -1;
        p1->ChangeSelection (0,false);
        p2->EnableMenu (false);
    }

    p2->ChangeSelection (0,false);
    UpdateBSButtons();
}

//__________________________________________________________

bool    _HYGeneralBootstrapWindow::ConfirmClose     (void)
{
    return _HYBootstrapWindow::ConfirmClose();
}

//__________________________________________________________

void    _HYGeneralBootstrapWindow::BuildCompatibleLF (bool)
{
    _HYPullDown * p1 = (_HYPullDown*) GetObject (5);
    while (p1->MenuItemCount()>1) {
        p1->DeleteMenuItem (1);
    }

    lfIDA.Clear();
    lfIDA << -1;
    dpIDA.Clear();
    dpIDA << -1;

    if (lfNullID>=0) {
        _LikelihoodFunction * nLF = (_LikelihoodFunction*) likeFuncList (lfNullID);

        long                nullSeqs  = nLF->SequenceCount (0),
                            nullSites = nLF->SiteCount     (),
                            k;

        _HYDataPanel*       nullData = nil;

        for (k=0; k<windowObjectRefs.lLength; k++) {
            _HYWindow * kthWindow = (_HYWindow*) windowObjectRefs (k);

            if (kthWindow->WindowKind () == HY_WINDOW_KIND_DATAPANEL) {
                nullData = (_HYDataPanel*) kthWindow;
                if (nullData->GetLFID()==lfNullID) {
                    break;
                }
            }
        }

        if (nullData&&((nullData->dataType==HY_DATAPANEL_NUCDATA)||(nullData->dataType==HY_DATAPANEL_PROTDATA)))
            for (k=0; k<windowObjectRefs.lLength; k++) {
                _HYWindow * kthWindow = (_HYWindow*) windowObjectRefs (k);

                if (kthWindow->WindowKind () == HY_WINDOW_KIND_DATAPANEL) {
                    _HYDataPanel * aDP = (_HYDataPanel*) kthWindow;
                    if (aDP->dataType == nullData->dataType) {
                        long         lID = aDP->GetLFID();
                        if ((lID>=0)&&(lID!=lfNullID)) {
                            _LikelihoodFunction * aLF = (_LikelihoodFunction*) likeFuncList (lID);

                            if ((aLF->SequenceCount(0) == nullSeqs) && (aLF->SiteCount() == nullSites)) {
                                _String dsName = *aDP->dataSetName;
                                p1->AddMenuItem (dsName,-1);
                                lfIDA << lID;
                                dpIDA << k;
                                for (long j=0; j<aDP->savedLFNames.lLength; j++) {
                                    _String dsSName = dsName & '[' & aDP->GetLFStateName (j) & ']';
                                    p1->AddMenuItem (dsSName, -1);
                                    lfIDA << lID;
                                    dpIDA << k;
                                }
                            }
                        }
                    }
                }
            }
    }

    p1->SetVisibleSize (p1->rel);

}

//__________________________________________________________

void    _HYGeneralBootstrapWindow::BuildAvailableLF  (void)
{
    _HYPullDown * p1 = (_HYPullDown*) GetObject (4);
    while (p1->MenuItemCount()>1) {
        p1->DeleteMenuItem (1);
    }

    lfID0.Clear();
    lfID0 << -1;
    dpID0.Clear();
    dpID0 << -1;

    for (long k=0; k<windowObjectRefs.lLength; k++) {
        _HYWindow * kthWindow = (_HYWindow*) windowObjectRefs (k);

        if (kthWindow->WindowKind () == HY_WINDOW_KIND_DATAPANEL) {
            _HYDataPanel * aDP = (_HYDataPanel*) kthWindow;
            if (aDP->GetLFID()>=0) {
                _String dsName = *aDP->dataSetName;
                p1->AddMenuItem (dsName,-1);
                lfID0 << aDP->GetLFID();
                dpID0 << k;
                for (long j=0; j<aDP->savedLFNames.lLength; j++) {
                    _String dsSName = dsName & '[' & aDP->GetLFStateName (j) & ']';
                    p1->AddMenuItem (dsSName, -1);
                    lfID0 << aDP->GetLFID();
                    dpID0 << k;
                }
            }
        }
    }

    p1->SetVisibleSize (p1->rel);

}

//__________________________________________________________

void    _HYGeneralBootstrapWindow::ResetResults (void)
{
    iterateCount  = 0;
    iterateCountP = 0;

    _HYTable* table = (_HYTable*)GetObject (1);
    _HYLabel* l     = (_HYLabel*)GetObject (0);

    for (long k=0; k<table->verticalSpaces.lLength; k++) {
        table->SetCellData(&empty,k,0,HY_TABLE_STATIC_TEXT,true);
        table->SetCellData(&empty,k,1,HY_TABLE_STATIC_TEXT,true);
    }

    _String cellValue ("P-Value: N/A");
    l->SetText (cellValue);
    table->SetVisibleSize(table->rel);
    table->_MarkForUpdate();
    handleGUI(true);
}

//__________________________________________________________

void    _HYGeneralBootstrapWindow::GetNullStateString (_String& s, long idx)
{
    long sidx = idx;
    while (sidx && (dpID0.lData[sidx-1] == dpID0.lData[sidx])) {
        sidx --;
    }

    sidx = idx - sidx;

    _HYDataPanel   *dp = (_HYDataPanel*)windowObjectRefs (dpID0.lData[idx]);

    if (sidx) {
        s = *dp->GetLFStateString (sidx-1);
    } else {
        _String  * cs = dp->LFSnapshot();
        s = *cs;
        DeleteObject (cs);
    }
}

//__________________________________________________________

void    _HYGeneralBootstrapWindow::GetAltStateString (_String& s, long idx)
{
    long sidx = idx;
    while (sidx && (dpIDA.lData[sidx-1] == dpIDA.lData[sidx])) {
        sidx --;
    }

    sidx = idx - sidx;

    _HYDataPanel   *dp = (_HYDataPanel*)windowObjectRefs (dpIDA.lData[idx]);

    if (sidx) {
        s = *dp->GetLFStateString (sidx-1);
    } else {
        _String  * cs = dp->LFSnapshot();
        s = *cs;
        DeleteObject (cs);
    }
}

//__________________________________________________________

void    _HYGeneralBootstrapWindow::HandleBootstrap (char np)
{
    _HYPullDown * p1  = (_HYPullDown*) GetObject (4),
                  * p2  = (_HYPullDown*) GetObject (5);

    long        iters = -1,
                sel1  = p1->GetSelection(),
                sel2  = p2->GetSelection(),
                k;

    bool        saveIterates = false;

    // prepare the remap

    _LikelihoodFunction* lfN = (_LikelihoodFunction*)likeFuncList (lfID0.lData[sel1]),
                         * lfA = (_LikelihoodFunction*)likeFuncList (lfIDA.lData[sel2]);


    _List                 seqsN,
                          seqsA;

    _SimpleList           smapN,
                          smapA,
                          remap1,
                          remap2,
                          remap4,
                          *nullFilters = &lfN->GetTheFilters ();

    for (k=0; k < nullFilters->lLength; k++) {
        _DataSetFilter* thisDF = (_DataSetFilter*)dataSetFilterList (nullFilters->lData[k]);
        if (k==0) {
            for (long k2 = 0; k2 < thisDF->theNodeMap.lLength; k2++) {
                seqsN << thisDF->GetSequenceName(k2);
                if (np) {
                    smapN << thisDF->theNodeMap.lData[k2];
                } else {
                    smapN << k2;
                }
            }
        }
        for (long j=0; j<thisDF->theOriginalOrder.lLength; j++) {
            remap1 << thisDF->theOriginalOrder.lData[j];
            remap2 << remap2.lLength;
        }
    }

    SortLists (&remap1, &remap2);

    nullFilters = &lfA->GetTheFilters ();

    unsigned long         maxOO = 0;


    for (k=0; k < nullFilters->lLength; k++) {
        _DataSetFilter* thisDF = (_DataSetFilter*)dataSetFilterList (nullFilters->lData[k]);
        if (k==0) {
            for (long k2 = 0; k2 < thisDF->theNodeMap.lLength; k2++) {
                seqsA << thisDF->GetSequenceName(k2);
                smapA << k2;
            }
        }
        for (long j=0; j<thisDF->theOriginalOrder.lLength; j++) {
            long locVal = thisDF->theOriginalOrder.lData[j];
            remap4 << locVal;
            if (maxOO <= locVal) {
                maxOO = locVal+1;
            }
        }
    }

    remap4.Sort();

    if (!remap1.Equal (remap4)) {
        _String incompatiblePartitions ("The likelihood functions you have selected do not span the same collection of sites. Bootstrap can't proceed.");
        ProblemReport (incompatiblePartitions, (Ptr)this);
        return;
    }

    SortLists (&seqsN, &smapN);
    SortLists (&seqsA, &smapA);

    if (!seqsN.Equal (seqsA)) {
        _String incompatiblePartitionS ("The likelihood functions you have selected do not span the same set of sequences. Bootstrap can't proceed.");
        ProblemReport (incompatiblePartitionS, (Ptr)this);
        return;
    }

    _SimpleList remap (maxOO),
                sremap (smapN.lLength);
    //*lmap = &((_DataSetFilter*)dataSetFilterList (nullFilters->lData[0]))->theNodeMap;

    for (k=0; k<remap4.lLength; k++) {
        remap.lData[remap4.lData[k]] = remap2.lData[k];
    }

    remap.lLength  = maxOO;

    sremap.lLength = smapN.lLength;

    for (k=0; k<smapN.lLength; k++) {
        sremap.lData[smapA.lData[k]] = smapN.lData[k];
    }

    _HYDataPanel   *nullPanel = (_HYDataPanel*)windowObjectRefs (dpID0.lData[sel1]),
                    *altPanel  = (_HYDataPanel*)windowObjectRefs (dpIDA.lData[sel2]);

    if (iterateCount>0) {
        _String warnMessage;
        warnMessage = _String ("Existing bootstrap results are about to be purged.");
        if (!ProceedPromptWithCheck (warnMessage, donotWarnAgain, warnAboutBootstrapReplace, (Ptr)this)) {
            return;
        }
        ResetResults();
    }

    _String  promptString ("Enter the number of iterations:"),
             valueString ("100"),
             folderSave;
    while (1) {
        if (!EnterStringDialogWithCheckbox (valueString, promptString, bsSavePrompt, saveIterates, (Ptr)this)) {
            return;
        }
        iters = valueString.toNum();
        if (iters>0) {
            break;
        } else {
            valueString = "100";
        }
    }

    if (saveIterates) {
        folderSave = ChooseAFolder (saveIteratesPrompt);
        if (folderSave.sLength==0) {
            return;
        }
    }

    nullPanel->SetLockState (true);
    altPanel->SetLockState (true);

    done = false;
    UpdateBSButtons();
    p1->EnableMenu (false);
    p2->EnableMenu (false);

    _String  * currentSnapshotN = nullPanel->LFSnapshot(),
               * currentSnapshotA = altPanel->LFSnapshot(),
                 nullState,
                 altState;


    GetAltStateString  (altState, sel2);
    GetNullStateString (nullState, sel1);

    _ExecutionList
    nullCommands (nullState),
                 altCommands  (altState);

    checkPointer (currentSnapshotN);
    checkPointer (currentSnapshotA);

    // compute LRT

    StartBarTimer ();


    nullCommands.Execute ();

    lfN->RescanAllVariables();
    lfN->PrepareToCompute ();

    nullLRT = -lfN->Compute();

    _Parameter saveNLL = nullLRT;

    lfN->DoneComputing ();

    altCommands.Execute();

    lfA->RescanAllVariables();
    lfA->PrepareToCompute ();
    nullLRT += lfA->Compute();
    nullLRT *= 2;
    lfA->DoneComputing ();

    _DataSet*originalData = (_DataSet*)dataSetList (nullPanel->dataSetID);
    _String *originalName = (_String*)dataSetNamesList (nullPanel->dataSetID);

    _String  dataSetName (*originalName);

    dataSetName = dataSetName & "_sim";
    FindUnusedObjectName (empty,dataSetName,dataSetNamesList);

    _DataSet *ds = nil;
    long     sID = -1;

    if (!np) {
        ds = new _DataSet();
        sID = AddDataSetToList (dataSetName,ds);
    }


    SetTitle  (GetTitle () & " (Running...)");

    for (long k = 0; k < iters; k++) {
        _String dataFileName = folderSave& "simulatedData" & (k+1) & ".seq";

        _Parameter  simLRT;

        _List       dummyN,
                    dummyA;

        _SimpleList dummy2N,
                    dummy2A;

        _DataSet*   npDS = nil;

        nullCommands.Execute ();
        lfN->RescanAllVariables();


        if (np) {
            nullPanel->SpawnLikelihoodFunctionNP (dummyN, np>1);
            npDS = nullPanel->GenerateOrderedDataSet ();
            sID  = AddDataSetToList (dataSetName,npDS);
            altPanel->SpawnLikelihoodFunction (npDS,&dataSetName, dummyA, dummy2A, &remap, &sremap);
        } else {
            ds->Clear();
            lfN->Simulate (*ds,dummyN);
            nullPanel->SpawnLikelihoodFunction (ds,&dataSetName, dummyN, dummy2N);
            altPanel->SpawnLikelihoodFunction (ds,&dataSetName, dummyA, dummy2A, &remap, &sremap);
        }


        if (np!=2) {
            valueString = _String("Obtaining MLEs for the null hypothesis. Iterate #")& k;
            SetStatusLine (valueString);
            _Matrix* res = lfN->Optimize();
            simLRT = -(*res)(1,0);
            DeleteObject (res);

            if (saveIterates) {
                _String     saveFile     = folderSave & "nullMLEs" & (k+1) & ".bf";
                if (np) {
                    nullPanel->SaveDataPanel (false, &saveFile, &dataFileName, false, npDS, true);
                } else {
                    nullPanel->SaveDataPanel (false, &saveFile, &dataFileName, false, ds);
                }
            }
        } else {
            simLRT = saveNLL;
        }

        altCommands.Execute ();
        lfA->RescanAllVariables();


        valueString = _String("Obtaining MLEs for the alternative hypothesis. Iterate #")& k;
        SetStatusLine (valueString);
        _Matrix * res = lfA->Optimize();
        simLRT += (*res)(1,0);
        DeleteObject (res);

        if (saveIterates) {
            _String     saveFile     = folderSave & "altMLEs" & (k+1) & ".bf";
            if (np) {
                altPanel->SaveDataPanel (false, &saveFile, &dataFileName, false, nil);
            } else {
                altPanel->SaveDataPanel (false, &saveFile, &dataFileName, false, nil);
            }
        }

        simLRT *= 2.;

        if (np) {
            nullPanel->SpawnLikelihoodFunctionNP (dummyN);
            altPanel->SpawnLikelihoodFunction (originalData,originalName, dummyA, dummy2A);
            KillDataSetRecord (sID);
            //DeleteObject (npDS);
        } else {
            nullPanel->SpawnLikelihoodFunction (originalData,originalName, dummyN, dummy2N);
            altPanel->SpawnLikelihoodFunction (originalData,originalName, dummyA, dummy2A);
        }

        AddIterate (simLRT);

        if (requestQuit) {
            break;
        }
    }

    StopBarTimer ();

    if (!np) {
        KillDataSetRecord (sID);
        //DeleteObject (ds);
    }

    _ExecutionList exl  (*currentSnapshotN),
                   exl2 (*currentSnapshotA);
    exl.Execute();
    exl2.Execute();
    DeleteObject (currentSnapshotN);
    DeleteObject (currentSnapshotA);
    nullPanel->SetLockState (false);
    altPanel->SetLockState (false);
    SetStopped (true);
    p1->EnableMenu (true);
    p2->EnableMenu (true);
    UpdateBSButtons();
    SetTitle  (GetTitle ().Cut (0, GetTitle ().sLength - 14));
    ReportAnalysisAsFinished("Bootstrap finished");
    requestQuit = false;
}

//__________________________________________________________

bool    _HYGeneralBootstrapWindow::ProcessEvent (_HYEvent* e)
{
    bool eDone = false;
    _String firstArg;
    long    i,f;

    if (e->EventClass()==_hyMenuOpenEvent) {
        i = MatchComponentID (e->EventCode());

        if (i==4) {
            BuildAvailableLF ();
            eDone = true;
        } else if (i==5) {
            BuildCompatibleLF();
            eDone = true;
        }
    } else if (e->EventClass()==_hyMenuSelChangeEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);

        if (i==4) {
            SetNullLF (lfID0 (e->EventCode().Cut(f+1,-1).toNum()));
            eDone = true;
        } else if (i==5) {
            lfAltID = lfIDA (e->EventCode().Cut(f+1,-1).toNum());
            UpdateBSButtons ();
            eDone = true;
        }
    } else if (e->EventClass() == _hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);

        if (i==3) {
            i =  e->EventCode().Cut (f+1,-1).toNum();
            if ((i==0)&&(HasStopped())) {
                int  h,v;
                _HYButtonBar* bb = (_HYButtonBar*)GetObject (3);
                _List        menuOptions;
                _String      menuChoice ("Parametric Bootstrap");

                menuOptions && & menuChoice;
                menuChoice = "Nonparametric Bootstrap";
                menuOptions && & menuChoice;
                menuChoice = "Permutation Test";
                menuOptions && & menuChoice;
                menuChoice = "Full Permutation Test";
                menuOptions && & menuChoice;
                bb->GetButtonLoc(0,h,v,true);
                menuChoice = HandlePullDown (menuOptions,h,v,0);
                bb->_UnpushButton();
                i = menuOptions.Find (&menuChoice);
                if (i>=0) {
                    HandleBootstrap (i);
                }

                eDone = true;
            }
        }
    }

    if (eDone) {
        DeleteObject (e);
        return true;
    }
    return _HYBootstrapWindow::ProcessEvent(e);
}

//__________________________________________________________

bool    _HYGeneralBootstrapWindow::ProcessGEvent (_HYEvent* e)
{
    _String firstArg,
            secondArg;
    long    k,f;

    bool    eDone = false;

    if (e->EventClass()==_hyGlobalLFKillEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        if (k!=GetID()) {
            f   = e->EventCode().Cut (f+1,-1).toNum();

            if ((lfNullID==f)||(lfAltID==f)) {
                _HYPullDown * p1 = (_HYPullDown*)GetObject (4),
                              * p2 = (_HYPullDown*)GetObject (5);

                if (done) {
                    if (lfAltID == f) {
                        p2->ChangeSelection (0,true);
                    } else {
                        p1->ChangeSelection (0,true);
                    }
                } else {
                    firstArg = "Internal Error: deleted active LF function!! Sorry about that.";
                    FlagError (firstArg);
                }
                eDone = true;
            }
        }
    }

    if (!eDone) {
        return _HYWindow::ProcessGEvent (e);
    }
    return true;
}

// EOF
