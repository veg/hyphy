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
#include "HYModelWindow.h"
#include "HYObjectInspector.h"
#include "HYLabel.h"
#include "HYTableComponent.h"
#include "HYButtonBar.h"
#include "HYUtils.h"
#include "HYEventTypes.h"
#include "HYParameterTable.h"
#include "likefunc.h"

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

_HYColor    boxFColor = {254,242,208},
            boxBColor = {102,90,88};


_String     objectInspectorTitle      ("Object Inspector"),
            objectInspectorKillPrompt (" is a part of one or more likelihood functions. Deleting it will cause those functions to be purged.");

extern      _List   dataSetNamesList,
            dataSetList,
            likeFuncNamesList,
            likeFuncList,
            modelTemplates;

extern      _SimpleList
windowObjects;

extern      _String
donotWarnAgain;

bool        warnTree      = false,
            warnDS      = false,
            warnModelOpen = false;

void        NewModel    (_String*);

_HYRect     objectInspectorRect     = {0,0,0,0,-1};
long        objectInspectorObjClass = 0;


//__________________________________________________________

void    TreeDependencies (_SimpleList& lfs, long treeID)
{
    _Variable* thisTree = LocateVar (treeID);

    for (long k=0; k< likeFuncNamesList.lLength; k++)
        if (((_String*)likeFuncNamesList(k))->sLength)
            if (((_LikelihoodFunction*)likeFuncList(k))->DependOnTree (*thisTree->GetName())>=0) {
                lfs << k;
            }
}

//__________________________________________________________

void    ModelDependencies (_SimpleList& lfs, long modelID)
{
    for (long k=0; k< likeFuncNamesList.lLength; k++)
        if (((_String*)likeFuncNamesList(k))->sLength)
            if (((_LikelihoodFunction*)likeFuncList(k))->DependOnModel (*(_String*)modelNames(modelID))>=0) {
                lfs << k;
            }
}

//__________________________________________________________

void    DSDependencies (_SimpleList& lfs, long dsID)
{

    for (long k=0; k< likeFuncNamesList.lLength; k++)
        if (((_String*)likeFuncNamesList(k))->sLength)
            if (((_LikelihoodFunction*)likeFuncList(k))->DependOnDS (dsID)>=0) {
                lfs << k;
            }
}

//__________________________________________________________
_HYObjectInspector::_HYObjectInspector(void):_HYTWindow (objectInspectorTitle, false)
{
    _HYRect         canvasSettings = {30,50,30,50,HY_COMPONENT_NO_SCROLL};
    _HYLabel*       l1 = new _HYLabel (canvasSettings,GetOSWindowData());
    canvasSettings.left  = 110;
    canvasSettings.right = 10000;
    _HYPullDown*    p1 = new _HYPullDown (canvasSettings,GetOSWindowData());
    p1->EnableMenu(true);
    canvasSettings.left = canvasSettings.right = 90;
    _HYButtonBar*   b1 = new _HYButtonBar(canvasSettings,GetOSWindowData());

    canvasSettings.top = 20;
    canvasSettings.left = 250;
    canvasSettings.bottom = 10000;
    canvasSettings.width = HY_COMPONENT_V_SCROLL;
    _HYTable*   ot  = new _HYTable (canvasSettings,GetOSWindowData(),1,2,125,100,HY_TABLE_STATIC_TEXT);
    canvasSettings.width = HY_COMPONENT_NO_SCROLL|HY_COMPONENT_BORDER_B|HY_COMPONENT_BORDER_T;
    _HYTable*   ot2 = new _HYTable (canvasSettings,GetOSWindowData(),1,2,125,20,HY_TABLE_STATIC_TEXT);


    ot->SetMessageRecipient (this);
    ot2->SetMessageRecipient (this);
    b1->SetMessageRecipient (this);
    p1->SetMessageRecipient (this);

    _String toolTipText ("Open window with the object");
    b1->AddButton (ProcureIconResource(HY_DATAPANEL_ICON_ID+7),&toolTipText);
    toolTipText = "Delete selected object";
    b1->AddButton (ProcureIconResource(HY_DATAPANEL_ICON_ID+3),&toolTipText);
    toolTipText = "Read object from file";
    b1->AddButton (ProcureIconResource(HY_DATAPANEL_ICON_ID+8),&toolTipText);
    toolTipText = "Create new object";
    b1->AddButton (ProcureIconResource(HY_DATAPANEL_ICON_ID+9),&toolTipText);

    b1->EnableButton(0,false);
    b1->EnableButton(1,false);
    b1->EnableButton(2,false);
    b1->EnableButton(3,false);
    b1->SetButtonDim(16);

    AddObject (l1);
    AddObject (ot);
    AddObject (ot2);
    AddObject (p1);
    AddObject (b1);
    SetTableDimensions (3,3);
    SetCell (0,0,b1);
    SetCell (0,1,p1);
    SetCell (0,2,p1);
    SetCell (1,0,ot2);
    SetCell (1,1,ot2);
    SetCell (1,2,ot2);
    SetCell (2,0,ot);
    SetCell (2,1,ot);
    SetCell (2,2,ot);

    _HYFont  labelFont;
#ifdef __WINDOZE__
    labelFont.face = "MS Sans Serif";
    labelFont.size = 10;
#else
#ifdef __HYPHY_GTK__
    labelFont.face = _HY_SANS_FONT;
    labelFont.size = 12;
#else
    labelFont.face = "Geneva";
    labelFont.size = 14;
#endif
#endif
    labelFont.style = HY_FONT_PLAIN;
    l1->SetFont (labelFont);
#ifdef __MAC__
    labelFont.face = "Times";
    labelFont.size = 12;
#endif
    ot->SetFont (labelFont);
    ot2->SetFont (labelFont);
    l1->SetBackColor (boxBColor);
    p1->SetBackColor (boxBColor);
    b1->SetBackColor (boxBColor);
    l1->SetForeColor (boxFColor);
    l1->SetAlignFlags (HY_ALIGN_LEFT);
    l1->SetShadow(true);
    p1->SetAlignFlags (HY_ALIGN_LEFT);
    b1->SetAlignFlags (HY_ALIGN_LEFT);
    ot->selectionType = HY_TABLE_SEL_ROWS|HY_TABLE_DONT_SIZE|HY_TABLE_FOCUSABLE ;
    ot2->selectionType = HY_TABLE_SEL_ROWS|HY_TABLE_DONT_SIZE;

    _String cellValue ("Object ID");
    ot2->SetCellData (&cellValue,0,0,HY_TABLE_BOLD|HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_CANTSELECT|HY_TABLE_PULLDOWN,true);
    cellValue = "Object Info";
    ot2->SetCellData (&cellValue,0,1,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_CANTSELECT|HY_TABLE_BOLD,true);

    l1->SetText (_String(" Class"));
    p1->AddMenuItem ("Trees",-1);
    p1->AddMenuItem ("Data Sets",-1);
    p1->AddMenuItem ("Models",-1);
    p1->AddMenuItem ("Likelihood Functions",-1);
    p1->ChangeSelection (objectInspectorObjClass,false);

    if (objectInspectorRect.width<0) {
        _HYRect  screenRect = GetScreenDimensions();
        SetWindowRectangle (0,0,screenRect.bottom-50,screenRect.right-20);
        SetPosition (screenRect.right-260,60);
        objectInspectorRect = GetWindowRect ();
    }

    firstTime = true;

    //p1->EnableItem  (2,false);

    DeleteObject (l1);
    DeleteObject (p1);
    DeleteObject (ot);
    DeleteObject (b1);

    lastKillID = 0;
}

//__________________________________________________________
_HYObjectInspector::~_HYObjectInspector()
{
    objectInspectorObjClass = ((_HYPullDown*)GetCellObject (0,2))->GetSelection();
}

//__________________________________________________________
bool    _HYObjectInspector::ConfirmClose()
{
    objectInspectorRect = GetWindowRect ();
    return true;
}

//__________________________________________________________

bool    _HYObjectInspector::ProcessEvent (_HYEvent* e)
{
    _String firstArg;
    long    k,i,f;
    bool    done = false;
    if (e->EventClass()==_hyMenuSelChangeEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }
        firstArg = e->EventCode().Cut (f+1,-1);
        k = firstArg.toNum();
        BuildListOfObjects (k);
        done = true;

    } else {
        if (e->EventClass()==_hyTableChangeSelEvent) {
            k = e->EventCode().toNum();
            for (i=0; i<components.lLength; i++)
                if (((_HYGuiObject*)components(i))->MatchID(k)) {
                    break;
                }
            if (i==1) {
                UpdateButtonsAndInfo ();
            }
            done = true;
        } else {
            if (e->EventClass()==_hyButtonPushEvent) {
                firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
                k = firstArg.toNum();
                for (i=0; i<components.lLength; i++)
                    if (((_HYGuiObject*)components(i))->MatchID(k)) {
                        break;
                    }

                firstArg = e->EventCode().Cut (f+1,-1);
                k = firstArg.toNum();
                if (i==4) { // button bar
                    _HYButtonBar* bb = (_HYButtonBar*)GetCellObject (0,0);
                    switch (k) {
                    case 0:
                        OpenObjectWindow();
                        break;
                    case 1:
                        KillObject ();
                        break;
                    case 2:
                        OpenObject ();
                        break;
                    case 3:
                        NewObject ();
                        break;
                    }
                    bb->_UnpushButton();
                    UpdateButtonsAndInfo();
                }

                done = true;
            } else {
                if (e->EventClass()==_hyTableDblClickEvent) {
                    k = e->EventCode().toNum();
                    for (i=0; i<components.lLength; i++)
                        if (((_HYGuiObject*)components(i))->MatchID(k)) {
                            break;
                        }
                    if (i==1) {
                        OpenObjectWindow();
                    }
                    done = true;
                } else if (e->EventClass()==_hyTablePullDownEvent) {
                    firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
                    k = firstArg.toNum();
                    for (i=0; i<components.lLength; i++) {
                        if (((_HYGuiObject*)components(i))->MatchID(k)) {
                            break;
                        }
                    }
                    k               = e->EventCode().Find(',',f+1,-1);
                    firstArg        = e->EventCode().Cut (k+1,-1);
                    if (i==2) {
                        SortObjectsByName (firstArg.toNum());
                    }

                    done = true;
                }

            }
        }
    }
    if (done) {
        DeleteObject (e);
        return true;
    }
    return _HYTWindow::ProcessEvent(e);
}

//__________________________________________________________

bool    _HYObjectInspector::ProcessGEvent (_HYEvent* e)
{
    _String firstArg;
    long    k,f;
    bool    done = false;

    lastKillID++;

    _HYPullDown*    p1 = (_HYPullDown*)GetCellObject (0,2);

    if (e->EventClass()==_hyGlobalTreeKillEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        if (MatchID(k)) {
            k = e->EventCode().Cut (f+1,-1).toNum();
            _Variable* theV = LocateVar(k);
            if (theV->ObjectClass()==TREE) {
                _SimpleList lfIDs;
                TreeDependencies (lfIDs, theV->GetAVariable());
                for (k=0; k<lfIDs.lLength; k++) {
                    postLFKillEvent (GetID(), lfIDs.lData[k]);
                }
                if (lfIDs.lLength == 0) {
                    DeleteVariable (*theV->GetName(),true);
                } else if (lastKillID<=2) {
                    postTreeKillEvent (GetID(), theV->GetAVariable());
                }
                BuildListOfObjects (p1->GetSelection());
            }
        }
        done = true;
    } else if (e->EventClass()==_hyGlobalLFKillEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        if (MatchID(k)) {
            k = e->EventCode().Cut (f+1,-1).toNum();
            KillLFRecord (k);
        }
        done = true;
    } else if (e->EventClass()==_hyGlobalDSKillEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        if (MatchID(k)) {
            k = e->EventCode().Cut (f+1,-1).toNum();
            _SimpleList lfIDs;
            DSDependencies (lfIDs, k);
            for (f=0; f<lfIDs.lLength; f++) {
                postLFKillEvent (GetID(), lfIDs.lData[f]);
            }
            if (lfIDs.lLength == 0) {
                KillDataSetRecord (k);
            } else if (lastKillID<=2) {
                postDSKillEvent (GetID(), k);
            }
            BuildListOfObjects (p1->GetSelection());
            done = true;
        }
    }

    lastKillID --;
    if (!done) {
        return _HYWindow::ProcessGEvent (e);
    }
    return true;
}

//__________________________________________________________

void    _HYObjectInspector::Update (Ptr p)
{
    _HYTWindow::Update(p);
}

//__________________________________________________________

void    _HYObjectInspector::Paint (Ptr p)
{
    _HYTWindow::Paint(p);
}

//__________________________________________________________

void    _HYObjectInspector::Activate ()
{
    _HYPullDown * p1 = (_HYPullDown*)GetCellObject (0,2);
    _HYTable*   oList = (_HYTable*)GetCellObject (HY_OBJECT_INSPECTOR_TABLE_ROW,0);
    _SimpleList saveSelection;
    oList->GetSelection (saveSelection);
    if ((objectInspectorRect.width>=0)&&firstTime) {
        Hide();
#ifdef __HYPHY_GTK__
        UpdateComponentInfo ();
#else
        dim = MinMaxWindowDimensions();
#endif
        SetWindowRectangle (0,0,dim.bottom,dim.right);

    }
    BuildListOfObjects (p1->GetSelection());
    oList->SetSelection (saveSelection);
    if ((objectInspectorRect.width>=0)&&firstTime) {
#ifdef __HYPHY_GTK__
        UpdateComponentInfo ();
#else
        dim = MinMaxWindowDimensions();
#endif
        SetWindowRectangle (0,0,objectInspectorRect.bottom-objectInspectorRect.top,objectInspectorRect.right-objectInspectorRect.left);
        SetPosition (objectInspectorRect.left,objectInspectorRect.top);
        objectInspectorRect = GetWindowRect ();
        Show();
    }
    _HYTWindow::Activate ();
    firstTime = false;
}

//__________________________________________________________

void    _HYObjectInspector::BuildListOfObjects (long index)
{
    _HYTable*       oList = (_HYTable*)GetCellObject (HY_OBJECT_INSPECTOR_TABLE_ROW,0);
    _HYTable*       oHead = (_HYTable*)GetCellObject (HY_OBJECT_INSPECTOR_TABLE_ROW-1,0);
    _HYButtonBar*   b1    = (_HYButtonBar*)GetCellObject (0,0);
    _HYPullDown*    p1    = (_HYPullDown*)GetCellObject (0,2);

    long            k,
                    rowCount = 0;

    oList->ClearTable();

    _String newStat;

    switch (index) {
    case 0: { // trees
        for (k=0; k<variablePtrs.lLength; k++) {
            _Variable * thisVar = LocateVar(k);
            if (thisVar&&(thisVar->ObjectClass()==TREE)) {
                oList->AddRow (-1,20,HY_TABLE_STATIC_TEXT);
                _SimpleList deps;
                TreeDependencies    (deps,thisVar->GetAVariable());
                if (deps.lLength) {
                    oList->SetCellData (thisVar->GetName(),rowCount,0,HY_TABLE_STATIC_TEXT|HY_TABLE_ITALIC|HY_TABLE_BOLD,true);
                } else {
                    oList->SetCellData (thisVar->GetName(),rowCount,0,HY_TABLE_STATIC_TEXT|HY_TABLE_ITALIC,true);
                }
                _TheTree * thisTree = (_TheTree*)thisVar;
                _PMathObj tc = thisTree->TipCount();
                newStat = _String ((long)tc->Value()) & " leaves ";
                DeleteObject (tc);
                tc = thisTree->BranchCount();
                newStat = newStat& ", " & _String ((long)tc->Value()+1) & " internal nodes.";
                DeleteObject (tc);
                oList->SetCellData (&newStat,rowCount,1,HY_TABLE_STATIC_TEXT,true);
                rowCount++;
            }
        }
        b1->EnableButton (2,true);
        b1->EnableButton (3,true);
        break;
    }
    case 1: { // data sets
        for (k=0; k<dataSetNamesList.lLength; k++) {
            _String* thisDS = (_String*)dataSetNamesList(k);
            if (thisDS->sLength) {
                oList->AddRow (-1,20,HY_TABLE_STATIC_TEXT);
                _SimpleList deps;
                DSDependencies  (deps,k);
                if (deps.lLength) {
                    oList->SetCellData (thisDS,rowCount,0,HY_TABLE_STATIC_TEXT|HY_TABLE_ITALIC|HY_TABLE_BOLD,true);
                } else {
                    oList->SetCellData (thisDS,rowCount,0,HY_TABLE_STATIC_TEXT|HY_TABLE_ITALIC,true);
                }
                _DataSet*  theDS = (_DataSet*)dataSetList(dataSetNamesList.Find(thisDS));
                if (theDS->GetTT()->IsStandardNucleotide()) {
                    newStat = "Nucleotide data";
                } else if (theDS->GetTT()->IsStandardAA()) {
                    newStat = "Aminoacid data";
                } else if (theDS->GetTT()->IsStandardBinary()) {
                    newStat = "Binary data";
                } else {
                    newStat = "Custom data";
                }

                newStat = newStat & ". "& _String (theDS->NoOfColumns()) &" sites ("&_String (theDS->NoOfUniqueColumns())&" distinct patterns), "&_String(theDS->NoOfSpecies())& " species.";
                oList->SetCellData (&newStat,rowCount,1,HY_TABLE_STATIC_TEXT,true);
                rowCount++;
            }
        }
        b1->EnableButton (2,true);
        b1->EnableButton (3,false);
        break;
    }
    case 2: { // Models
        for (k=0; k<modelTemplates.lLength; k++) {
            _List       * thisModel       = (_List*) modelTemplates(k);
            _SimpleList * modelParameters = (_SimpleList*) (*thisModel)(1);

            oList->AddRow (-1,20,HY_TABLE_STATIC_TEXT);

            if (modelParameters->lData[0] & HY_DATAPANEL_MODEL_MODELS) {
                oList->SetCellData ((*thisModel)(0),rowCount,0,HY_TABLE_STATIC_TEXT,true);
            } else {
                oList->SetCellData ((*thisModel)(0),rowCount,0,HY_TABLE_STATIC_TEXT|HY_TABLE_ITALIC,true);
            }

            if (modelParameters->lData[1]==4) {
                newStat = "Nucleotide model";
            } else if (modelParameters->lData[1]==16) {
                newStat = "Dinucleotide Model";
            } else if (modelParameters->lData[1]==20) {
                newStat = "Aminoacid Model";
            } else {
                newStat = "Codon Model";
            }

            newStat = newStat & ".";
            oList->SetCellData (&newStat,rowCount,1,HY_TABLE_STATIC_TEXT,true);
            rowCount++;
        }

        for (k=0; k<modelNames.lLength; k++) {
            _String * modelID = (_String*)modelNames(k);
            if (modelID->sLength) {
                oList->AddRow (-1,20,HY_TABLE_STATIC_TEXT);

                _SimpleList deps;
                ModelDependencies   (deps,k);

                if (deps.lLength) {
                    oList->SetCellData (modelID,rowCount,0,HY_TABLE_STATIC_TEXT|HY_TABLE_BOLD,true);
                } else {
                    oList->SetCellData (modelID,rowCount,0,HY_TABLE_STATIC_TEXT,true);
                }

                _Matrix * modelMatrix = (_Matrix*)LocateVar (modelMatrixIndices.lData[k])->GetValue();

                long    mDim = modelMatrix->GetHDim();

                if (mDim==4) {
                    newStat = "Nucleotide model";
                } else if (mDim==16) {
                    newStat = "Dinucleotide Model";
                } else if (mDim==20) {
                    newStat = "Aminoacid Model";
                } else {
                    newStat = "Codon Model";
                }

                newStat = newStat & ".";
                oList->SetCellData (&newStat,rowCount,1,HY_TABLE_STATIC_TEXT,true);
                rowCount++;
            }
        }

        b1->EnableButton (2,true);
        b1->EnableButton (3,true);
        break;
    }
    case 3: { // likelihood functions
        for (k=0; k<likeFuncNamesList.lLength; k++) {
            _String *lfName = (_String*)likeFuncNamesList (k);

            if (lfName->sLength) {
                oList->AddRow (-1,20,HY_TABLE_STATIC_TEXT);
                oList->SetCellData (lfName,rowCount,0,HY_TABLE_STATIC_TEXT|HY_TABLE_BOLD,true);
                _LikelihoodFunction *lf = (_LikelihoodFunction*)likeFuncList (k);
                newStat = _String ("Likelihood Function on ") & (long)lf->GetTheTrees().lLength & " partitions/trees.";

                oList->SetCellData (&newStat,rowCount,1,HY_TABLE_STATIC_TEXT,true);
                rowCount++;
            }
        }
        b1->EnableButton (2,true);
        b1->EnableButton (3,false);
        break;
    }
    }
    if (rowCount==0) {
        newStat = _String("No ")& *p1->GetMenuItem(index) &" found.";
        oList->AddRow (-1,20,HY_TABLE_STATIC_TEXT);
        oList->SetCellData (&newStat,rowCount,0,HY_TABLE_STATIC_TEXT|HY_TABLE_ITALIC|HY_TABLE_CANTSELECT,true);
        oList->SetCellData (&empty,rowCount,1,HY_TABLE_STATIC_TEXT|HY_TABLE_ITALIC|HY_TABLE_CANTSELECT,true);
        rowCount ++;
    } else if (rowCount>1) {
        newStat = _String (rowCount) & ' ' & *p1->GetMenuItem(index) &" found.";
    } else {
        newStat = _String("One ") & p1->GetMenuItem(index)->Cut(0,p1->GetMenuItem(index)->sLength-2) & " found.";
    }

    if (rowCount<5) {
        oList->AddRow (-1,(5-rowCount)*20,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT);
    }

    SetStatusBar (newStat);
    oList -> AutoFitWidth (*oHead);
    rowCount = oHead->GetMaxW();
    if (rowCount<oHead->rel.right-oHead->rel.left) {
        oList->SetColumnSpacing (1,-rowCount+oHead->rel.right-oHead->rel.left+1,false);
        oHead->SetColumnSpacing (1,-rowCount+oHead->rel.right-oHead->rel.left+1,false);
    }
    rowCount = oList->GetMaxH();
    if (rowCount<=oList->rel.bottom-oList->rel.top) {
        oList->SetRowSpacing (oList->verticalSpaces.lLength-1,oList->rel.bottom-oList->rel.top+1-rowCount,false);
    }

    oList->SetVisibleSize(oList->rel);
    oHead->SetVisibleSize(oHead->rel);
    oList->_MarkForUpdate ();
    oHead->_MarkForUpdate ();
#ifdef __HYPHY_GTK__
    UpdateComponentInfo ();
#else
    dim = MinMaxWindowDimensions();
#endif
}

//__________________________________________________________

void    _HYObjectInspector::UpdateButtonsAndInfo (void)
{
    _HYTable*       dl = (_HYTable*)GetCellObject (HY_OBJECT_INSPECTOR_TABLE_ROW,0);
    _HYButtonBar*   b1 = (_HYButtonBar*)GetCellObject (0,0);

    _SimpleList     sel;
    dl->GetSelection (sel);
    if (sel.lLength) {
        b1->EnableButton (0,true);
        b1->EnableButton (1,true);
    } else {
        b1->EnableButton (0,false);
        b1->EnableButton (1,false);
    }
}

//__________________________________________________________

void    _HYObjectInspector::OpenObjectWindow (void)
{
    _HYTable*       dl = (_HYTable*)GetCellObject (HY_OBJECT_INSPECTOR_TABLE_ROW,0);
    _HYButtonBar*   b1 = (_HYButtonBar*)GetCellObject (0,0);
    _HYPullDown*    p1 = (_HYPullDown*)GetCellObject (0,2);

    _SimpleList     sel;
    dl->GetSelection (sel);
    if (sel.lLength) {
        b1->EnableButton (0,true);
        for (long k = 0; k<sel.lLength; k+=2) {
            _String          windowTitle (*(_String*)dl->GetCellData (0,sel.lData[k]/2));
            long f;
            switch (p1->GetSelection()) {
            case 0: {
                _String winTitle (windowTitle);
                windowTitle = _String ("Tree ")&windowTitle;
                f = FindWindowByName (windowTitle);
                if (f<0) {
                    _HYTreePanel* newTreePanel = new _HYTreePanel (winTitle,winTitle);
                    newTreePanel->_Zoom(true);
                    newTreePanel->BringToFront();
                    //newTreePanel->Show();
                } else {
                    _HYPlatformWindow* thisWindow = (_HYPlatformWindow*)windowObjects(f);
                    thisWindow->_Activate();
                }
                break;
            }
            case 1: {
                _String winTitle (windowTitle);
                windowTitle = _String ("DataSet ")&windowTitle;
                f = FindWindowByName (windowTitle);
                if (f<0) {
                    _HYDataPanel* newDataPanel = new _HYDataPanel (winTitle,winTitle);
                    newDataPanel->BringToFront();
                } else {
                    _HYPlatformWindow* thisWindow = (_HYPlatformWindow*)windowObjects(f);
                    thisWindow->_Activate();
                }
                break;
            }
            case 2: {
                _String winTitle (windowTitle);
                windowTitle = _String ("Model ")&windowTitle;
                f = FindWindowByName (windowTitle);
                if (f<0) {
                    if (dl->cellTypes.lData[sel.lData[k]] & HY_TABLE_ITALIC) {
                        if (!warnModelOpen) {
                            windowTitle = _String ("Model \"") & winTitle & "\" is a template model, and it can't be opened. You can open a new model based on it, though.";
                            if (!warnModelOpen) {
                                if (!ProceedPromptWithCheck (windowTitle, donotWarnAgain, warnModelOpen)) {
                                    break;
                                }
                            }
                            NewObject ();
                            break;

                        }
                    } else {
                        if (sel.lData[k]/2<modelTemplates.lLength) {
                            _List* modelSpec = FindModelTemplate (&winTitle);
                            if (modelSpec) {
                                OpenModelFromFile(*(_String*)(*modelSpec)(2));
                            }
                        } else {
                            OpenModelFromMatrix (winTitle);
                        }
                    }
                } else {
                    _HYWindow* thisWindow = (_HYWindow*)windowObjectRefs(f);
                    thisWindow->BringToFront();
                }
                break;
            }
            case 3: {
                _String winTitle (windowTitle);
                windowTitle = _String ("Likelihood Function ")&windowTitle;
                f = FindWindowByName (windowTitle);
                if (f<0) {
                    f = likeFuncNamesList.Find(&winTitle);
                    if (f>=0) {
                        _HYParameterTable* newParameterTable = new _HYParameterTable (windowTitle,f);
                        newParameterTable->_Zoom(true);
                        newParameterTable->BringToFront();
                    }
                } else {
                    _HYPlatformWindow* thisWindow = (_HYPlatformWindow*)windowObjects(f);
                    thisWindow->_Activate();
                }
                break;
            }
            }
        }
    }
}

//__________________________________________________________

void    _HYObjectInspector::KillObject (void)
{
    _HYTable*       dl = (_HYTable*)GetCellObject (HY_OBJECT_INSPECTOR_TABLE_ROW,0);
    _HYButtonBar*   b1 = (_HYButtonBar*)GetCellObject (0,0);
    _HYPullDown*    p1 = (_HYPullDown*)GetCellObject (0,2);

    _SimpleList     sel;
    dl->GetSelection (sel);
    if (sel.lLength) {
        b1->EnableButton (0,true);
        long f;
        switch (p1->GetSelection()) {
        case 0: {
            for (f = 0; f<sel.lLength; f+=2) {
                long treeID = LocateVarByName (*(_String*)dl->GetCellData (0,sel.lData[f]/2));
                if (treeID>=0)
#ifndef USE_AVL_NAMES
                    treeID = variableReindex.lData[treeID];
#else
                    treeID = variableNames.GetXtra(treeID);
#endif
                if (dl->cellTypes.lData[sel.lData[f]] & HY_TABLE_BOLD) {
                    _String prompt ("Tree ");
                    prompt = prompt & *LocateVar (treeID)->GetName() & objectInspectorKillPrompt;
                    if (warnTree||(!ProceedPromptWithCheck (prompt,donotWarnAgain,warnTree))) {
                        continue;
                    }
                }
                postTreeKillEvent (GetID(), treeID);
            }
            break;
        }
        case 1: {
            for (f = 0; f<sel.lLength; f+=2) {
                long dsID = dataSetNamesList.Find ((_String*)dl->GetCellData (0,sel.lData[f]/2));
                if (dsID>=0) {
                    if (dl->cellTypes.lData[sel.lData[f]] & HY_TABLE_BOLD) {
                        _String prompt ("Dataset ");
                        prompt = prompt & *(_String*)dataSetNamesList(dsID) & objectInspectorKillPrompt;
                        if (!((!warnTree)&&ProceedPromptWithCheck (prompt,donotWarnAgain,warnDS))) {
                            continue;
                        }
                    }
                    postDSKillEvent (GetID(), dsID);
                }
            }
            break;
        }
        case 2: {
            for (f = 0; f<sel.lLength; f+=2) {
                if (sel.lData[f]/2<modelTemplates.lLength) {
                    _List* modelT = FindModelTemplate((_String*)dl->GetCellData (0,sel.lData[f]/2));
                    if (modelT) { // template model
                        _String*fileLocation = (_String*)(*modelT)(2),
                                errMsg (*(_String*)dl->GetCellData (0,sel.lData[f]/2));

                        if (dl->cellTypes.lData[sel.lData[f]] & HY_TABLE_ITALIC) {
                            errMsg = errMsg & " is a template model. If you wish to delete it, please do so by removing the appropriate file from the \"Substitution Models\" directory.";
                            ProblemReport (errMsg);
                        } else {
                            errMsg = errMsg & " is about to be deleted. This action is NOT undoable.";
                            if (ProceedPrompt (errMsg,(Ptr)this)) {
                                if (remove (fileLocation->sData)) {
                                    errMsg = "File delete operation failed.";
                                    ProblemReport (errMsg);
                                } else {
                                    errMsg = *(_String*)dl->GetCellData (0,sel.lData[f]/2);
                                    for (long k=0; k<modelTemplates.lLength; k++) {
                                        _List* thisList = (_List*)modelTemplates(k);
                                        if (errMsg.Equal((_String*)(*thisList)(0))) {
                                            modelTemplates.Delete (k);
                                        }
                                    }
                                    BuildListOfObjects (2);
                                }
                            }
                        }
                    }
                }
            }
            break;
        }
        }
    }
}

//__________________________________________________________

void    _HYObjectInspector::OpenObject (void)
{
    _HYPullDown*    p1 = (_HYPullDown*)GetCellObject (0,2);
    switch (p1->GetSelection()) {
    case 0: {
        OpenTreeFile();
        break;
    }
    case 1: {
        OpenDataFile();
        break;
    }
    case 2: {
        OpenModelFile();
        break;
    }
    default:
        return;
    }
    BuildListOfObjects (p1->GetSelection());
}

//__________________________________________________________

void    _HYObjectInspector::NewObject (void)
{
    _HYPullDown*    p1 = (_HYPullDown*)GetCellObject (0,2);
    _HYTable*       dl = (_HYTable*)GetCellObject (HY_OBJECT_INSPECTOR_TABLE_ROW,0);
    switch (p1->GetSelection()) {
    case 0: {
        NewTreeWindow (-1);
        break;
    }
    case 2: {
        _SimpleList sel;
        dl->GetSelection (sel);
        _String*    md = nil;
        if (sel.lLength==2) {
            md = (_String*)dl->GetCellData(0,sel.lData[0]/2);
        }

        NewModel(md);
        break;
    }
    default:
        return;
    }
    BuildListOfObjects (p1->GetSelection());
}

//__________________________________________________________

void    _HYObjectInspector::SortObjectsByName (long location)
{

    _HYTable*       table = (_HYTable*)GetCellObject (HY_OBJECT_INSPECTOR_TABLE_ROW,0);

    _List              menuItems;

    menuItems.AppendNewInstance(new _String ("Ascending"));
    menuItems.AppendNewInstance(new _String ("Descending"));

    long    h = table->_HandlePullDown  (menuItems,(location&0xffff0000)>>16,location&0x0000ffff,0),
            k;

    if (h>=0) {
        menuItems.Clear();
        _SimpleList index;
        for (k=0; k<table->verticalSpaces.lLength; k++) {
            index<<k;
            menuItems << table->GetCellData (0,k);
        }
        bool    hasPadding = ((_String*)menuItems(menuItems.lLength-1))->sLength==0;

        if (hasPadding) {
            menuItems.Delete (menuItems.lLength-1);
            index.Delete (index.lLength-1);
        }
        SortLists (&menuItems, &index);
        if (h==1) {
            index.Flip();
        }
        if (hasPadding) {
            index<<index.lLength;
        }

        table->SetRowOrder (index);
    }
}


//EOF