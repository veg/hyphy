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


#include "HYModelWindow.h"
#include "HYButtonBar.h"
#include "HYLabel.h"
#include "HYTableComponent.h"
#include "HYGraphicPane.h"
#include "HYUtils.h"
#include "HYEventTypes.h"
#include "HYDataPanel.h"
#include "HYButton.h"
#include "HYDialogs.h"
#include "HYTextBox.h"

#include "parser.h"
#include "likefunc.h"
#include "batchlan.h"

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

extern  _SimpleList windowObjects;
extern  _HYColor    tableDefaultBk2,
        labelBColor,
        labelFColor;

extern  _String     blReplicate,
        donotWarnAgain,
        modelName,
        modelOptions,
        modelDimension,
        modelFunction,
        globalPrefix,
        categoryPrefix,
        nucDataType,
        dinucDataType,
        codonDataType,
        proteinDataType,
        none,
        parameterOption[3],
        modelEFVVector,
        multiplyByFrequencies;


_String             noneDefined        ("None Defined"),
                    defineNewParameter ("New Parameter"),
                    editParameter       ("Edit Parameter"),
                    deleteParameter     ("Delete Parameter"),
                    updateParemeterList("Remove Unused Parameters"),
                    moveToExpression   ("Set Expression To Parameter"),
                    addToExpression     ("Append Parameter To Expression"),
                    spawnParameters    ("The formula contains parameters not defined previously. Should I add them to the list (as locals) and proceed?"),
                    glStart         ("(__Global__"),
                    catStart            ("(__Category__"),
                    emptySelection      ("Empty selection"),
                    noneApplicable      ("None Applicable"),
                    modelGenCode    ("ModelGeneticCode"),
                    buildTemps          ("CreateTemplates"),
                    buildClasses    ("CreateSubClasses"),
                    nameMatrixT         ("TemplateNames"),
                    nameMatrixS         ("SubClassNames"),
                    resultName          ("MatrixTemplate"),
                    embeddedInModel    ("User Defined Vector"),
                    gammaVariation      ("Unit Mean Gamma"),
                    modelMatrixName    ("ModelMatrixName"),
                    EFVModifier         ("EFVModifierMatrix"),
                    EFVCodeMarker       ("/** MATRIX MODIFICATION CODE **/"),
                    ModelMatrixDimension
                    ("ModelMatrixDimension"),
                    titlePrefix         ("Model "),
                    untitled        ("untitled"),
                    ModelEFVType    ("Model_EFV_Type"),
                    modelDataType       ("Model_Data_Type"),
                    rateClassCount     ("rateClassCount"),
                    invalidEFV          ("Can't add up to 1"),
                    exprNotValidated   ("Your expression was not validated."),
                    deleteParameters   ("Deleting this parameter will erase all matrix elements which depend on it.");

_HYColor midGrey        = {220,220,220},
         darkGrey      = {160,160,160},
         editBoxBad      = {128,0,0},
         editBoxGood = {0,128,0};

bool     warnNewParameters = false;

long     modelSpawnerCount = 1;

_String  directorySep (
#ifdef __MAC__
    ':'
#else
#ifdef __WINDOZE__
    '\\'
#else
    '/'
#endif
#endif
);

//__________________________________________________________
_HYModelWindow::_HYModelWindow (_String name, _List* labels, char t):_HYTWindow (titlePrefix&name, true)
{

    long            index,
                    cellWidth;

    stateLabels.Duplicate (labels);

    _HYRect         canvasSettings = {26,82,26,82,HY_COMPONENT_NO_SCROLL};

    _HYLabel*       l1    = new _HYLabel (canvasSettings, GetOSWindowData());


    canvasSettings.right  = 40;
    canvasSettings.left   = 40;
    _HYLabel*       l5    = new _HYLabel (canvasSettings, GetOSWindowData());

    canvasSettings.width  = HY_COMPONENT_BORDER_T;
    _HYLabel*       l3    = new _HYLabel (canvasSettings, GetOSWindowData()),
    *       l4    = new _HYLabel (canvasSettings, GetOSWindowData());

    canvasSettings.right  = 130;
    canvasSettings.left   = 130;
    _HYPullDown*    p2    = new _HYPullDown (canvasSettings, GetOSWindowData());

    canvasSettings.right  = 10000;
    _HYPullDown*    p3    = new _HYPullDown (canvasSettings, GetOSWindowData());

    canvasSettings.width  = HY_COMPONENT_NO_SCROLL;

    _HYPullDown*    p4    = new _HYPullDown (canvasSettings, GetOSWindowData());
    canvasSettings.left   = 92;

    _HYPullDown*    p1    = new _HYPullDown (canvasSettings, GetOSWindowData());

    canvasSettings.left   = canvasSettings.right  = 170;
    canvasSettings.top    = canvasSettings.bottom = 26;
    canvasSettings.width  = HY_COMPONENT_BORDER_T;

    _HYButtonBar*   bb    = new _HYButtonBar (canvasSettings, GetOSWindowData());

    canvasSettings.left   = 120;
    canvasSettings.right  = 120;
    _HYButtonBar*   bb2   = new _HYButtonBar (canvasSettings, GetOSWindowData());

    canvasSettings.left   = 100;
    canvasSettings.right  = 10000;
    _HYTextBox*   tl      = new _HYTextBox (canvasSettings, GetOSWindowData());

    canvasSettings.left   = canvasSettings.right  = 26;
    _HYButtonBar*   bb3   = new _HYButtonBar (canvasSettings, GetOSWindowData());


    canvasSettings.right  = 10000;
    canvasSettings.left   = 60;
    canvasSettings.top    = 20;
    canvasSettings.bottom = 20;

    canvasSettings.width  = HY_COMPONENT_BORDER_B|HY_COMPONENT_BORDER_T;

    _HYTable*   columnHeaders   = new _HYTable (canvasSettings,GetOSWindowData(),1,labels->lLength+1,25,20,HY_TABLE_STATIC_TEXT);

    for (index=0; index < labels->lLength; index++) {
        columnHeaders->SetCellData (stateLabels(index),0,index,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_BOLD,false);
        stateLabels(index)->nInstances++;
    }

    columnHeaders->SetCellData (&empty,0,index,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_BOLD|HY_TABLE_CANTSELECT,true);

    canvasSettings.left   = canvasSettings.right = 20;
    canvasSettings.top    = 50;
    canvasSettings.bottom = 10000;

    canvasSettings.width = HY_COMPONENT_BORDER_R;

    _HYTable*   rowHeaders  = new _HYTable (canvasSettings,GetOSWindowData(),labels->lLength+1,1,25,20,HY_TABLE_STATIC_TEXT);
    //rowHeaders->selectionType = HY_TABLE_DONT_SIZE;

    for (index=0; index < labels->lLength; index++) {
        rowHeaders->SetCellData (stateLabels(index),0,index,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_BOLD,false);
        stateLabels(index)->nInstances++;
    }

    rowHeaders->SetCellData (&empty,0,index,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_BOLD|HY_TABLE_CANTSELECT,true);
    rowHeaders->SetRowSpacing (index,HY_SCROLLER_WIDTH-20,false);
    rowHeaders->AutoFitColumn(0,false,true);

    rowHeaders->settings.left = rowHeaders->settings.right = canvasSettings.left = canvasSettings.right =
                                    rowHeaders->GetColumnSpacing (0);

    columnHeaders->SetColumnSpacing (index,HY_SCROLLER_WIDTH-25,false);

    cellWidth               = 400/labels->lLength+1;
    canvasSettings.top      = canvasSettings.bottom     = 20;
    canvasSettings.width    = HY_COMPONENT_BORDER_T|HY_COMPONENT_BORDER_B|HY_COMPONENT_BORDER_R;

    _HYLabel*       l2      = new _HYLabel (canvasSettings, GetOSWindowData());

    canvasSettings.top      = canvasSettings.left   = 50;
    canvasSettings.bottom   = canvasSettings.right  = 10000;

    canvasSettings.width    = HY_COMPONENT_V_SCROLL|HY_COMPONENT_H_SCROLL;

    _HYTable*   table       = new _HYTable (canvasSettings,GetOSWindowData(),labels->lLength,labels->lLength,cellWidth,20,HY_TABLE_EDIT_TEXT);

    canvasSettings.right    = canvasSettings.left = 169-rowHeaders->GetRowSpacing(0);
    canvasSettings.width    = HY_COMPONENT_NO_SCROLL|HY_COMPONENT_BORDER_R;
    canvasSettings.top      = 50;
    canvasSettings.bottom   = 10000;
    _HYTable*   EFVColumn   = new _HYTable (canvasSettings,GetOSWindowData(),labels->lLength+1,1,169-rowHeaders->GetRowSpacing(0),20,HY_TABLE_EDIT_TEXT);
    EFVColumn->SetCellData (&empty,labels->lLength-1,0,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT,true);
    EFVColumn->SetCellData (&empty,labels->lLength,0,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_BOLD|HY_TABLE_CANTSELECT,true);

    canvasSettings.top      = canvasSettings.bottom = 20;
    canvasSettings.width    = HY_COMPONENT_NO_SCROLL|(HY_COMPONENT_BORDER-HY_COMPONENT_BORDER_L);
    _HYLabel*   EFVLabel    = new _HYLabel (canvasSettings,GetOSWindowData());

    EFVColumn->SetRowSpacing (labels->lLength,HY_SCROLLER_WIDTH-20,false);
    EFVColumn->selectionType = HY_TABLE_DONT_SIZE;
    EFVColumn->SetMessageRecipient (this);

    _String s (" * ");
    table->SetCellData (&s,0,0,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT|HY_TABLE_BOLD,true);
    for (index=1; index < labels->lLength; index++) {
        table->SetCellData (table->GetCellData(0,0),index,index,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT|HY_TABLE_BOLD,true);
    }

    table->AutoFitWidth(*columnHeaders);
    table->EnforceWidth(400,stateLabels.lLength-1);
    columnHeaders->EnforceWidth(401,stateLabels.lLength-1);

    table->selectionType      = HY_TABLE_DONT_SIZE;
    rowHeaders->selectionType = HY_TABLE_DONT_SIZE;

    columnHeaders->SetColumnSpacing (index,HY_SCROLLER_WIDTH-columnHeaders->GetColumnSpacing(index),false);

    table->SetMessageRecipient (this);
    table->SetBackColor (midGrey);
    EFVColumn->SetBackColor (midGrey);
    rowHeaders->SetMessageRecipient (this);
    columnHeaders->SetMessageRecipient (this);
    bb->SetMessageRecipient (this);
    bb2->SetMessageRecipient (this);
    bb3->SetMessageRecipient (this);
    p1->SetMessageRecipient (this);
    p2->SetMessageRecipient (this);
    p4->SetMessageRecipient (this);
    tl->SetMessageRecipient (this);

    AddObject (table);                          // 0
    AddObject (rowHeaders);                     // 1
    AddObject (columnHeaders);                  // 2
    AddObject (bb);                             // 3
    AddObject (p1);                             // 4
    AddObject (p2);                             // 5
    AddObject (p3);                             // 6
    AddObject (bb2);                            // 7
    AddObject (l1);                             // 8
    AddObject (l2);                             // 9
    AddObject (l3);                             // 10
    AddObject (l4);                             // 11
    AddObject (EFVColumn);                      // 12
    AddObject (EFVLabel);                       // 13
    AddObject (l5);                             // 14
    AddObject (p4);                             // 15
    AddObject (tl);                             // 16
    AddObject (bb3);                            // 17


    SetTableDimensions (5,7);
    SetCell   (MODEL_PARAMETER_ROW,0,l1);
    SetCell   (MODEL_PARAMETER_ROW,1,l1);
    SetCell   (MODEL_PARAMETER_ROW,2,l1); // width 82
    SetCell   (MODEL_PARAMETER_ROW,3,p1);
    SetCell   (MODEL_PARAMETER_ROW,4,l5);
    SetCell   (MODEL_PARAMETER_ROW,5,p4);
    SetCell   (MODEL_PARAMETER_ROW,6,p4);

    SetCell   (MODEL_CLASS_ROW,0,l3); // 40
    SetCell   (MODEL_CLASS_ROW,1,l3); // 40
    SetCell   (MODEL_CLASS_ROW,2,p2); // 130
    SetCell   (MODEL_CLASS_ROW,3,p2);
    SetCell   (MODEL_CLASS_ROW,4,l4); // 40
    SetCell   (MODEL_CLASS_ROW,5,p3);
    SetCell   (MODEL_CLASS_ROW,6,p3);

    SetCell   (MODEL_BUTTON_ROW,0,bb); // width 170
    SetCell   (MODEL_BUTTON_ROW,1,bb);
    SetCell   (MODEL_BUTTON_ROW,2,bb);
    SetCell   (MODEL_BUTTON_ROW,3,bb);
    SetCell   (MODEL_BUTTON_ROW,4,bb2);
    SetCell   (MODEL_BUTTON_ROW,5,tl);
    SetCell   (MODEL_BUTTON_ROW,6,bb3);

    SetCell   (MODEL_MATRIX_HEADER_ROW,0,l2);
    SetCell   (MODEL_MATRIX_HEADER_ROW,1,EFVLabel);
    SetCell   (MODEL_MATRIX_HEADER_ROW,2,EFVLabel);
    SetCell   (MODEL_MATRIX_HEADER_ROW,3,EFVLabel);
    SetCell   (MODEL_MATRIX_HEADER_ROW,4,columnHeaders);
    SetCell   (MODEL_MATRIX_HEADER_ROW,5,columnHeaders);
    SetCell   (MODEL_MATRIX_HEADER_ROW,6,columnHeaders);

    SetCell   (MODEL_MATRIX_ROW,0,rowHeaders);
    SetCell   (MODEL_MATRIX_ROW,1,EFVColumn);
    SetCell   (MODEL_MATRIX_ROW,2,EFVColumn);
    SetCell   (MODEL_MATRIX_ROW,3,EFVColumn);
    SetCell   (MODEL_MATRIX_ROW,4,table);
    SetCell   (MODEL_MATRIX_ROW,5,table);
    SetCell   (MODEL_MATRIX_ROW,6,table);

    _HYFont  labelFont;
    bb->SetBackColor (labelBColor);
    bb2->SetBackColor(labelBColor);
    bb3->SetBackColor(labelBColor);
    tl->SetBackColor (labelBColor);
    l1->SetBackColor (labelBColor);
    l1->SetForeColor (labelFColor);
    l5->SetBackColor (labelBColor);
    l5->SetForeColor (labelFColor);
    l2->SetBackColor (darkGrey);
    l2->SetForeColor (labelFColor);
    p1->SetBackColor (labelBColor);
    p4->SetBackColor (labelBColor);
    l3->SetBackColor (darkGrey);
    l4->SetBackColor (darkGrey);
    EFVLabel->SetBackColor (tableDefaultBk2);
//  l3->SetForeColor (labelBColor);
//  l4->SetForeColor (labelBColor);
    p2->SetBackColor (darkGrey);
    p3->SetBackColor (darkGrey);

#ifdef __MAC__
    labelFont.face = "Geneva";
    labelFont.size = 10;
#endif
#ifdef __WINDOZE__
    labelFont.face = "MS Sans Serif";
    labelFont.size = 12;
#endif
#ifdef __HYPHY_GTK__
    labelFont.face = _HY_SANS_FONT;
    labelFont.size = 10;
#endif

    labelFont.style = HY_FONT_PLAIN;

    table->SetFont(labelFont);
    rowHeaders->SetFont(labelFont);
    columnHeaders->SetFont(labelFont);
    labelFont.style = HY_FONT_BOLD;
    EFVLabel->SetFont (labelFont);

    labelFont.style = HY_FONT_PLAIN;
#ifdef __MAC__
    labelFont.size = 12;
#endif
#ifdef __WINDOZE__
    labelFont.size = 14;
#endif
#ifdef __HYPHY_GTK__
    labelFont.size = 11;
#endif
    l1->SetFont (labelFont);
    l3->SetFont (labelFont);
    l4->SetFont (labelFont);
    l5->SetFont (labelFont);
    EFVLabel->SetAlignFlags (HY_ALIGN_LEFT);

    l1->SetText (_String("Parameters:"));
    l1->SetShadow(true);
    l1->SetAlignFlags (HY_ALIGN_LEFT);

    l5->SetText (_String("Eq.Fr:"));
    l5->SetShadow(true);
    l5->SetAlignFlags (HY_ALIGN_LEFT);

    l3->SetText (_String("Class"));
    l4->SetText (_String("To "));
    l3->SetAlignFlags (HY_ALIGN_LEFT);
    l4->SetAlignFlags (HY_ALIGN_RIGHT);

    EFVLabel->SetText ("Equilibruim Freqs.");

    bb->SetAlignFlags (HY_ALIGN_LEFT);
    bb->SetButtonDim(16);
    bb2->SetAlignFlags (HY_ALIGN_LEFT);
    bb2->SetButtonDim(16);
    bb3->SetAlignFlags (HY_ALIGN_LEFT);
    bb3->SetButtonDim(16);

    p1->AddMenuItem(noneDefined,-1);
    p1->AddMenuItem(menuSeparator,-1);
    p1->AddMenuItem(_String(defineNewParameter),-1);
    p1->EnableItem (0,false);
    p1->SetAlignFlags (HY_ALIGN_LEFT);

    p2->AddMenuItem(noneDefined,-1);
    p2->EnableItem (0,false);
    p2->SetAlignFlags (HY_ALIGN_LEFT);
    p3->AddMenuItem(noneApplicable,-1);
    p3->EnableItem (0,false);
    p3->SetAlignFlags (HY_ALIGN_LEFT);
    p4->SetAlignFlags (HY_ALIGN_LEFT);

    s = "Autofit column widths";
    bb->AddButton (ProcureIconResource(HY_MODELWINDOW_ICON_BASE),&s);
    s = "Copy Cell to Clipboard";
    bb->AddButton (ProcureIconResource(HY_MODELWINDOW_ICON_BASE+1),&s);
    s = "Paste Clipboard to Selection";
    bb->AddButton (ProcureIconResource(HY_MODELWINDOW_ICON_BASE+2),&s);
    s = "Clear Cells in Selection";
    bb->AddButton (ProcureIconResource(HY_MODELWINDOW_ICON_BASE+3),&s);
    s = "Select Cells";
    bb->AddButton (ProcureIconResource(HY_MODELWINDOW_ICON_BASE+4),&s);
    s = "Add Multiplicative Factor to Selection";
    bb->AddButton (ProcureIconResource(HY_MODELWINDOW_ICON_BASE+5),&s);
    bb->EnableButton (1,false);
    bb->EnableButton (2,false);
    bb->EnableButton (3,false);
    bb->EnableButton (5,false);
    bb->MarkAsPullDown (4,true);
    bb->MarkAsPullDown (5,true);

    s = "Replace Current Selection";
    bb2->AddButton (ProcureIconResource(HY_MODELWINDOW_ICON_BASE+10),&s);
    s = "Union w/ Current Selection";
    bb2->AddButton (ProcureIconResource(HY_MODELWINDOW_ICON_BASE+11),&s);
    s = "Intersection w/ Current Selection";
    bb2->AddButton (ProcureIconResource(HY_MODELWINDOW_ICON_BASE+12),&s);
    s = "XOR w/ Current Selection";
    bb2->AddButton (ProcureIconResource(HY_MODELWINDOW_ICON_BASE+13),&s);
    s = "Subtract from  Current Selection";
    bb2->AddButton (ProcureIconResource(HY_MODELWINDOW_ICON_BASE+14),&s);
    bb2->EnableButton (0,false);
    bb2->EnableButton (1,false);
    bb2->EnableButton (2,false);
    bb2->EnableButton (3,false);
    bb2->EnableButton (4,false);

    s = "Validate Expression";
    bb3->AddButton (ProcureIconResource(HY_MODELWINDOW_ICON_BASE+20),&s);

    tl->SetText       ("Type Expression Here");
    tl->SetForeColor  (editBoxBad);

    // set up the table

    lastParameterChoice = -1;

    SetStatusBar (emptySelection);


    EFVOptions && & embeddedInModel;
    EFVChoice = 0;

    rateOptions && & gammaVariation;
    rateChoice = 0;

    type = t;

    GrabEFVs      (type);
    GrabRateVariation ();

    _String      efs (1./stateLabels.lLength);

    for (index = 0; index < stateLabels.lLength-1; index++) {
        EFVColumn->SetCellData (&efs, index, 0, HY_TABLE_EDIT_TEXT, true);
    }

    EFVColumn->SetCellData (&efs,labels->lLength-1,0,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT,true);


    _HYRect  screenRect = GetScreenDimensions();
    SetWindowRectangle (0,0,screenRect.bottom-50,screenRect.right-20);
    DoFitColumns ();
    //SetPosition (5,40);
    CenterWindow(this);

    screenRect = MinMaxWindowDimensions();

    DeleteObject (l1);
    DeleteObject (l2);
    DeleteObject (l3);
    DeleteObject (l4);
    DeleteObject (l5);
    DeleteObject (p1);
    DeleteObject (p2);
    DeleteObject (p3);
    DeleteObject (p4);
    DeleteObject (bb);
    DeleteObject (bb2);
    DeleteObject (bb3);
    DeleteObject (table);
    DeleteObject (rowHeaders);
    DeleteObject (columnHeaders);
    DeleteObject (EFVLabel);
    DeleteObject (EFVColumn);
    DeleteObject (tl);

}

//__________________________________________________________

bool    _HYModelWindow::ProcessEvent (_HYEvent* e)
{
    bool    done = false;
    _String firstArg;
    long    i,f,g,k;
    if (e->EventClass()==_hyScrollingEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        g = e->EventCode().Find(',',f+1);
        if (g>=0) {
            k = firstArg.toNum();
            for (i=0; i<components.lLength; i++) {
                if (((_HYGuiObject*)components(i))->MatchID(k)) {
                    break;
                }
            }
            _HYTable*      theTable  = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,4);
            _HYTable*      EFVColumn = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,1);

            if (i==0) {
                firstArg = e->EventCode().Cut (f+1,g-1);
                k = firstArg.toNum();
                if (k) {
                    theTable = (_HYTable*)  GetCellObject (MODEL_MATRIX_HEADER_ROW,4);
                    theTable->SetMessageRecipient (nil);
                    theTable->ProcessEvent (generateScrollEvent (k,0));
                } else {
                    firstArg = e->EventCode().Cut (g+1,-1);
                    k = firstArg.toNum();
                    if (k) {
                        theTable = (_HYTable*)  GetCellObject (MODEL_MATRIX_ROW,0);
                        theTable->SetMessageRecipient (nil);
                        theTable->ProcessEvent (generateScrollEvent (0,k));
                        if (EFVChoice==0) {
                            EFVColumn->ProcessEvent (generateScrollEvent (0,k));
                        }
                    }
                }
            } else {
                theTable->SetMessageRecipient (nil);
                if (i==1) { // row headers
                    firstArg = e->EventCode().Cut (g+1,-1);
                    k = firstArg.toNum();
                    theTable->ProcessEvent (generateScrollEvent (0,k));
                    if (EFVChoice==0) {
                        EFVColumn->ProcessEvent (generateScrollEvent (0,k));
                    }
                } else if (i==2) {
                    firstArg = e->EventCode().Cut (f+1,g-1);
                    k = firstArg.toNum();
                    theTable->ProcessEvent (generateScrollEvent (k,0));
                }
            }
            theTable->SetMessageRecipient (this);
            done = true;
        }
    } else if (e->EventClass()==_hyMenuOpenEvent) {
        k = e->EventCode().toNum();
        for (i=0; i<components.lLength; i++)
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }

        if (i==4) {
            PrepareParameterMenu();
            done = true;
        }
    } else if (e->EventClass()==_hyMenuSelChangeEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++)
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        if (i==4) {
            k = e->EventCode().Cut (f+1,-1).toNum();
            _HYPullDown * p1 = (_HYPullDown*) GetCellObject (MODEL_PARAMETER_ROW,3);
            f = p1->MenuItemCount ();
            if (k == f-7) {
                lastParameterChoice = ParameterEditBox (-1);
            } else if (k==f-6) {
                lastParameterChoice = ParameterEditBox (lastParameterChoice);
            } else if (k==f-5) {
                if (DeleteVarID (lastParameterChoice,true)) {
                    lastParameterChoice = -1;
                }
            } else {
                if (k==f-1) {
                    lastParameterChoice = -1;
                    CheckDependencies();
                } else if ((k==f-2)||(k==f-3)) {
                    _HYTextBox * tl = (_HYTextBox*) GetCellObject (MODEL_BUTTON_ROW,5);
                    _String    curValue = tl->GetText(),
                               pName = RetrieveVarID(lastParameterChoice);
                    if (k==f-2) {
                        curValue = curValue&pName;
                    } else {
                        curValue = pName;
                    }

                    tl->SetText (curValue);
                }
            }

            if (k>=f-7) {
                PrepareParameterMenu  ();
                SetParameterMenuState (lastParameterChoice);
                if ((k==f-5)||(k==f-1)) {
                    PrepareParameterMenu();
                }
            } else {
                lastParameterChoice = MenuChoiceToID (k);
            }

            done = true;
        } else if (i==5) {
            ProcessClassMenu();
            done = true;
        } else if (i==15) {
            _HYPullDown * p4 = (_HYPullDown*) GetCellObject (MODEL_PARAMETER_ROW,5);
            if (p4->GetSelection() != EFVChoice) {
                SetEFVChoice (p4->GetSelection());
            }
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
        firstArg        = e->EventCode().Cut (f+1,-1);
        if (i==0) {
            HandleCellEditEvent (firstArg.toNum());
        } else if (i==12) {
            HandleEFVEditEvent (firstArg.toNum());
        }

        done = true;
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
                table = (_HYTable*)     GetCellObject (MODEL_MATRIX_HEADER_ROW,4);
            } else {
                table = (_HYTable*)     GetCellObject (MODEL_MATRIX_ROW,4);
            }
            table->SetColumnSpacing (f,k,true);
            dim = MinMaxWindowDimensions();
            done = true;
        }
    } else if (e->EventClass()==_hyTableChangeSelEvent) {
        k = e->EventCode().toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }

        done = true;

        _HYTable* t1 = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,4),
                  * t2 = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,0),
                    * t3 = (_HYTable*)GetCellObject (MODEL_MATRIX_HEADER_ROW,4),
                      * t4 = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,1);

        if (i&&(i<3)) {
            _SimpleList tSel, tSel2;
            t2->GetSelection (tSel);
            t3->GetSelection (tSel2);
            if (tSel.lLength||tSel2.lLength) {
                if (i==1) {
                    t3->ClearSelection();
                    t2->SetMessageRecipient (nil);
                    t1->SetRowSelection(tSel);
                    t2->SetMessageRecipient (this);
                    UpdateButtonState();
                } else {
                    t2->ClearSelection();
                    t2->SetMessageRecipient (nil);
                    t1->SetColumnSelection(tSel2);
                    t2->SetMessageRecipient (this);
                    UpdateButtonState();
                }
            }
        } else {
            if (i==0) {
                if (t2->messageRecipient) {
                    t2->ClearSelection();
                    t3->ClearSelection();
                }
                _SimpleList newSelection;
                t1->GetSelection (newSelection);
                UpdateButtonState();

                if (newSelection.lLength==0) {
                    SetStatusBar (emptySelection);
                } else {
                    _String sb;
                    k = newSelection.lData[0]/stateLabels.lLength;
                    f = newSelection.lData[0]%stateLabels.lLength;

                    sb = (_String)('(') & k
                         & ','
                         & f
                         & "): "
                         & *(_String*)stateLabels(k)
                         & ">>"
                         & *(_String*)stateLabels(f);

                    if (newSelection.lLength>1) {
                        long k2 = newSelection.lData[newSelection.lLength-1]/stateLabels.lLength,
                             f2 = newSelection.lData[newSelection.lLength-1]%stateLabels.lLength;

                        sb = sb& "...(" & k2
                             & ','
                             & f2
                             & "): "
                             & *(_String*)stateLabels(k2)
                             & ">>"
                             & *(_String*)stateLabels(f2);


                        if (newSelection.lLength != (f2-f+1)*(k2-k+1)) {
                            sb = _String ("Part of ") & sb;
                        }
                    }
                    sb = sb & " -- " & _String((long)newSelection.lLength) & " cells (" &
                         100.*(_Parameter)newSelection.lLength/(_Parameter)(stateLabels.lLength*stateLabels.lLength)
                         & "%)";

                    SetStatusBar (sb);
                }
            }
        }

        if (EFVChoice==0) {
            if (i<3) {
                t4->SetMessageRecipient (nil);
                t4->ClearSelection ();
                t4->SetMessageRecipient (this);
            } else {
                if (i==12) {
                    t1->SetMessageRecipient (nil);
                    t1->ClearSelection ();
                    t1->SetMessageRecipient (this);
                    t2->SetMessageRecipient (nil);
                    t2->ClearSelection ();
                    t2->SetMessageRecipient (this);
                    t3->SetMessageRecipient (nil);
                    t3->ClearSelection ();
                    t3->SetMessageRecipient (this);
                    UpdateButtonState();
                }
            }
        }
    } else if (e->EventClass()==_hyButtonPushEvent) {
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
            switch  (k) {
            case 0:
                DoFitColumns();
                break;

            case 1: // copy
                DoCopyCell  ();
                break;

            case 2: // paste
                DoPasteToCells ();
                break;

            case 3: // clear
                DoClearCells ();
                break;

            case 4:
                DoSelectCells ();
                break;

            case 5:
                DoMultiplyCells ();
                break;

            }
            ((_HYButtonBar*)GetCellObject (MODEL_BUTTON_ROW,0))->_UnpushButton();
            done = true;
        } else if (i==7) {
            firstArg        = e->EventCode().Cut (f+1,-1);
            k               = firstArg.toNum();
            ProcessTemplateSelection (k);
            done = true;
        } else if (i==17) {
            _HYTextBox*    tl = (_HYTextBox*)  GetCellObject (MODEL_BUTTON_ROW,5);
            _HYButtonBar* bb3 = (_HYButtonBar*)GetCellObject (MODEL_BUTTON_ROW,6);

            _String cText (tl->GetText());
            _Formula f (cText,nil,false);
            if (f.GetList().lLength) {
                clipboardString = cText;
                SyncEditBox ();
                UpdateButtonState();
            } else {
                ProblemReport (exprNotValidated,(Ptr)this);
            }
            bb3->_UnpushButton ();
            done = true;
        }
    } else if (e->EventClass()==_hyTextEditChange) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++) {
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }
        }
        if (i==16) {
            firstArg        = e->EventCode().Cut (f+1,-1);
            k               = firstArg.toNum();
            if (k==1) {
                if (clipboardString.sLength) {
                    clipboardString = empty;
                    _HYTextBox*    tl = (_HYTextBox*)  GetCellObject (MODEL_BUTTON_ROW,5);
                    _HYButtonBar*  bb3 = (_HYButtonBar*)GetCellObject (MODEL_BUTTON_ROW,6);
                    tl->SetForeColor (editBoxBad);
                    bb3->EnableButton (0,true);
                    UpdateButtonState();
                }
            } else if (k==2) {
                _HYTextBox*    tl = (_HYTextBox*)  GetCellObject (MODEL_BUTTON_ROW,5);
                if (clipboardString.sLength==0) {
                    _HYButtonBar*  bb3 = (_HYButtonBar*)GetCellObject (MODEL_BUTTON_ROW,6);
                    ProcessEvent   (generateButtonPushEvent (bb3->GetID(),0));
                } else {
                    tl->UnfocusComponent ();
                    keyboardFocus = -1;
                }
            }
            done = true;
        }
    }

    if (done) {
        DeleteObject (e);
        return true;
    }
    return _HYTWindow::ProcessEvent(e);
}

//__________________________________________________________

void    _HYModelWindow::PrepareParameterMenu (void)
{
    _HYPullDown * p1 = (_HYPullDown*) GetCellObject (MODEL_PARAMETER_ROW,3);
    p1->DeleteAllItems();
    if (globalVariables.lLength+localVariables.lLength+categoryVariables.lLength) {
        long k;
        for (k=0; k<localVariables.lLength; k++) {
            p1->AddMenuItem (*(_String*)localVariables(k),-1);
        }
        if (globalVariables.lLength+categoryVariables.lLength) {
            p1->AddMenuItem (glStart,-1);
        }
        for (k=0; k<globalVariables.lLength; k++) {
            p1->AddMenuItem (*(_String*)globalVariables(k),-1);
        }
        if (categoryVariables.lLength) {
            p1->AddMenuItem (catStart,-1);
        }
        for (k=0; k<categoryVariables.lLength; k++) {
            p1->AddMenuItem (*(_String*)categoryVariables(k),-1);
        }
    } else {
        p1->AddMenuItem (noneDefined,-1);
    }

    p1->AddMenuItem (menuSeparator,-1);
    p1->AddMenuItem (defineNewParameter,-1);
    p1->AddMenuItem (editParameter,-1);
    p1->AddMenuItem (deleteParameter,-1);
    p1->AddMenuItem (menuSeparator,-1);
    p1->AddMenuItem (moveToExpression,-1);
    p1->AddMenuItem (addToExpression,-1);
    p1->AddMenuItem (updateParemeterList,-1);
    if (lastParameterChoice<0) {
        p1->EnableItem (p1->MenuItemCount()-5,false);
        p1->EnableItem (p1->MenuItemCount()-6,false);
        p1->EnableItem (p1->MenuItemCount()-2,false);
        p1->EnableItem (p1->MenuItemCount()-3,false);
    }

    if (globalVariables.lLength+localVariables.lLength+categoryVariables.lLength==0) {
        p1->EnableItem (p1->MenuItemCount()-1,false);
    }
}

//__________________________________________________________

long    _HYModelWindow::FindVarID (_String& s)
{
    long f = localVariables.Find(&s);

    if (f>=0) {
        return f;
    }

    f = globalVariables.Find (&s);
    if (f>=0) {
        return f+localVariables.lLength;
    }


    f = categoryVariables.Find (&s);
    if (f>=0) {
        return f+localVariables.lLength+globalVariables.lLength;
    }

    return -1;
}

//__________________________________________________________

_String&    _HYModelWindow::RetrieveVarID (long ID)
{
    if (ID<localVariables.lLength) {
        return *(_String*)localVariables(ID);
    }

    if (ID<localVariables.lLength+globalVariables.lLength) {
        return *(_String*)globalVariables(ID-localVariables.lLength);
    }

    return *(_String*)categoryVariables(ID-localVariables.lLength-globalVariables.lLength);

}

//__________________________________________________________

char    _HYModelWindow::VarKind (long ID)
{
    if (ID<localVariables.lLength) {
        return 0;
    }

    if (ID<localVariables.lLength+globalVariables.lLength) {
        return 1;
    }

    return 2;

}

//__________________________________________________________

long    _HYModelWindow::MenuChoiceToID (long mc)
{
    if (mc<localVariables.lLength) {
        return mc;
    }

    if (mc<localVariables.lLength+globalVariables.lLength+1) {
        return mc-1;
    }

    return mc-2;

}

//__________________________________________________________

void    _HYModelWindow::InsertVarID (_String& vID, char type)
{
    _Variable v (vID);
    _List* toIns = nil;

    switch (type) {
    case 0:
        toIns = &localVariables;
        break;
    case 1:
        toIns = &globalVariables;
        break;
    case 2:
        toIns = &categoryVariables;
        break;
    }

    if (toIns) {
        toIns->BinaryInsert (&vID);
    }
}

//__________________________________________________________

void    _HYModelWindow::RenameVarID (_String& vID, long ID, char type)
{
    _Variable v (vID);
#ifndef USE_AVL_NAMES
    UpdateParameterName (variableReindex.lData[LocateVarByName (vID)],
                         variableReindex.lData[LocateVarByName (RetrieveVarID (ID))]);
#else
    UpdateParameterName (variableNames.GetXtra(LocateVarByName (vID)),
                         variableNames.GetXtra(LocateVarByName (RetrieveVarID (ID))));
#endif
    DeleteVarID (ID);
    InsertVarID (vID,type);
}

//__________________________________________________________

bool    _HYModelWindow::DeleteVarID (long ID, bool clean)
{
    char type = VarKind (ID);

    if (clean)
#ifndef USE_AVL_NAMES
        if (!RemoveParameterName (variableReindex.lData[LocateVarByName (RetrieveVarID (ID))]))
#else
        if (!RemoveParameterName (variableNames.GetXtra(LocateVarByName (RetrieveVarID (ID)))))
#endif
            return false;

    switch (type) {
    case 0:
        localVariables.Delete (ID);
        break;
    case 1:
        globalVariables.Delete (ID-localVariables.lLength);
        break;
    case 2:
        categoryVariables.Delete (ID-localVariables.lLength-globalVariables.lLength);
        break;
    }
    return true;
}

//__________________________________________________________

void    _HYModelWindow::SetParameterMenuState (long ID)
{
    _HYPullDown * p1 = (_HYPullDown*) GetCellObject (MODEL_PARAMETER_ROW,3);

    if (ID>=0) {
        p1->ChangeSelection (ID+VarKind(ID),false);
    } else {
        if (localVariables.lLength+globalVariables.lLength+categoryVariables.lLength) {
            lastParameterChoice = 0;
        } else {
            lastParameterChoice = -1;
        }
        p1->ChangeSelection (0,false);
    }
}

//__________________________________________________________

void    _HYModelWindow::SetEFVVector (_String& name)
{
    long f = EFVOptions.Find (&name);
    if (f>=0) {
        _HYPullDown* efv = (_HYPullDown*)GetCellObject (MODEL_PARAMETER_ROW,5);
        efv->ChangeSelection (f,true);
        //SetEFVChoice(f);
    }
}

//__________________________________________________________

void    _HYModelWindow::CopyEFVMatrix (_Matrix* m)
{
    _HYTable * EFVs = (_HYTable*)components(12);
    bool    good = false;
    long    k;

    if ((m->GetHDim()*m->GetVDim()==stateLabels.lLength)&&(!m->IsVariable())) {
        _Parameter partSum = 0.0, v;
        for (k=0; k<stateLabels.lLength-1; k++) {
            v = m->theData[k];
            partSum+=v;
            if ((v<0.0)||(partSum>1.0)) {
                break;
            }
            _String converter (v);
            EFVs->SetCellData (&converter,k,0,HY_TABLE_EDIT_TEXT,true);
        }
        if (k==stateLabels.lLength-1) {
            good = true;
        }
    }

    if (good) {
        UpdateEFVLastEntry ();
        EFVs->_MarkForUpdate();
    } else {
        _Matrix dummy (1,stateLabels.lLength,false,true);
        for (k=0; k<stateLabels.lLength; k++) {
            dummy.theData[k] = 1./stateLabels.lLength;
        }

        CopyEFVMatrix (&dummy);
    }
}

//__________________________________________________________

void    _HYModelWindow::HandleCellEditEvent (long index)
{

    _HYTable*       table = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,4);

    long v = index%table->horizontalSpaces.lLength,
         h = index/table->horizontalSpaces.lLength;

    bool u1 = false,
         u2 = true;

    _String* cellEntry = (_String*)table->GetCellData (v,h);

    if (cellEntry->sLength) {
        _Formula f (*cellEntry);
        if (f.IsEmpty()) {
            cellEntry = (_String*)table->GetCellData(h,v);
            table->SetCellData (cellEntry,h,v,HY_TABLE_EDIT_TEXT,false);
            u1 = true;
            u2 = false;
        } else {
            if (!AddFormulaParametersToList(f)) {
                cellEntry = (_String*)table->GetCellData(h,v);
                table->SetCellData (table->GetCellData(h,v),h,v,HY_TABLE_EDIT_TEXT,false);
                u1 = true;
                u2 = false;
            }

            if (u2) {
                cellEntry = (_String*)f.toStr();
                table->SetCellData (cellEntry,h,v,HY_TABLE_EDIT_TEXT,false);
                table->SetCellData (cellEntry,v,h,HY_TABLE_EDIT_TEXT,false);
                u1 = true;
            }
        }
    } else {
        table->SetCellData (cellEntry,v,h,HY_TABLE_EDIT_TEXT,true);
    }

    cellEntry->nInstances++;
    if (u1) {
        table->_MarkCellForUpdate (index);
    }
    if (u2) {
        table->_MarkCellForUpdate (v*table->horizontalSpaces.lLength+h);
    }

    taint = u1|u2;

}

//__________________________________________________________

void    _HYModelWindow::HandleEFVEditEvent (long index)
{

    _HYTable*       table     = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,1);
    _String*        cellEntry = (_String*)table->GetCellData (index,0);
    _String         newValue  (0.0),
                    errMsg ("Frequency entry failed to evaluate to a number in [0,1].");

    if (cellEntry->sLength) {
        _Formula f (*cellEntry);
        if (!f.IsEmpty()) {
            _PMathObj fv = f.Compute();
            if (fv&&(fv->ObjectClass()==NUMBER)) {
                _Parameter freq = fv->Value();
                if ((freq>=0.0)&&(freq<=1.0)) {
                    newValue = freq;
                    errMsg = empty;
                }
            }
        }
    }

    if (errMsg.sLength) {
        ProblemReport (errMsg,(Ptr)this);
    }

    table->SetCellData (&newValue,0,index,HY_TABLE_EDIT_TEXT,true);
    UpdateEFVLastEntry ();
}

//__________________________________________________________

void    _HYModelWindow::UpdateEFVLastEntry (void)
{

    _HYTable*       table     = (_HYTable*)components (12);
    _String*        cellEntry;
    _Parameter      sum = 0.0;

    for (long   k=0; k<stateLabels.lLength-1; k++) {
        cellEntry = (_String*)table->GetCellData (k,0);
        sum += cellEntry->toNum();
    }

    if (sum<=1.0) {
        _String lastCell (1.-sum);
        table->SetCellData (&lastCell,0,stateLabels.lLength-1,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT,true);
    } else {
        table->SetCellData (&invalidEFV,0,stateLabels.lLength-1,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT,true);
    }

    table->_MarkCellForUpdate (stateLabels.lLength-1);
    taint = true;
}

//__________________________________________________________

void    _HYModelWindow::SetEFVChoice (long c)
{
    if (c!=EFVChoice) {
        bool update = false;
        if (c==0) {
            _HYGuiObject  *EFVLabel  = (_HYGuiObject*)components(13),
                           *EFVColumn = (_HYGuiObject*)components(12);

            SetCell   (MODEL_MATRIX_HEADER_ROW,1,EFVLabel);
            SetCell   (MODEL_MATRIX_HEADER_ROW,2,EFVLabel);
            SetCell   (MODEL_MATRIX_HEADER_ROW,3,EFVLabel);
            SetCell   (MODEL_MATRIX_ROW,1,EFVColumn);
            SetCell   (MODEL_MATRIX_ROW,2,EFVColumn);
            SetCell   (MODEL_MATRIX_ROW,3,EFVColumn);

            update = true;
        } else {
            if (EFVChoice == 0) {
                _HYGuiObject  *EFVLabel  = (_HYGuiObject*)components(2),
                               *EFVColumn = (_HYGuiObject*)components(0);

                _HYTable      *col = (_HYTable*)components(12);

                SetCell   (MODEL_MATRIX_HEADER_ROW,1,EFVLabel);
                SetCell   (MODEL_MATRIX_HEADER_ROW,2,EFVLabel);
                SetCell   (MODEL_MATRIX_HEADER_ROW,3,EFVLabel);
                SetCell   (MODEL_MATRIX_ROW,1,EFVColumn);
                SetCell   (MODEL_MATRIX_ROW,2,EFVColumn);
                SetCell   (MODEL_MATRIX_ROW,3,EFVColumn);


                col->SetMessageRecipient (nil);
                col->ClearSelection();
                col->SetMessageRecipient (this);

                update = true;
            }
        }
        if (update) {
            SetWindowRectangle (top,left,bottom,right);
        }

        EFVChoice = c;
    }
}

//__________________________________________________________

void    _HYModelWindow::UpdateButtonState ()
{
    _HYButtonBar * bb = (_HYButtonBar*) GetCellObject (MODEL_BUTTON_ROW,0);
    _HYTable     * t  = (_HYTable*)     GetCellObject (MODEL_MATRIX_ROW,4);

    _SimpleList  s;

    bool         b1 = false,
                 b2 = false,
                 b3;

    t->GetSelection (s);
    b3 = s.lLength;

    if (clipboardString.sLength)
        if (s.lLength) {
            b2 = true;
        }

    if (s.lLength==1) {
        _String* sv = (_String*)t->GetCellData(s.lData[0]%t->horizontalSpaces.lLength,
                                               s.lData[0]/t->horizontalSpaces.lLength);

        if (sv->sLength) {
            b1 = true;
        }
    }

    bb->EnableButton (1,b1);
    bb->EnableButton (2,b2);
    bb->EnableButton (3,b3);
    bb->EnableButton (5,b3);
    _UpdateEditMenu  (b1,b2);
}

//__________________________________________________________

void    _HYModelWindow::DoCopyCell ()
{
    _HYTable     * t  = (_HYTable*)     GetCellObject (MODEL_MATRIX_ROW,4);

    _SimpleList  s;
    t->GetSelection (s);
    clipboardString = *(_String*)t->GetCellData(s.lData[0]%t->horizontalSpaces.lLength,
                      s.lData[0]/t->horizontalSpaces.lLength);
    SyncEditBox     ();
}

//__________________________________________________________

bool    _HYModelWindow::DoEditModelName ()
{
    _String curName (GetTitle()),
            prompt ("Enter a new model name:"),
            curTitle;

    curName.Trim (titlePrefix.sLength,-1);
    curTitle = curName;

    if (EnterStringDialog (curName,prompt, (Ptr)this)) {
        if (!curName.Equal (&curTitle)) {
            if (!curName.IsValidIdentifier()) {
                curName = _String ('"') & curName & "\" is not a valid identifier.";
                ProblemReport (curName,(Ptr)this);
                return false;
            }
            if (FindModelTemplate (&curName)) {
                curName = _String ('"') & curName & "\" is already being used.";
                ProblemReport (curName,(Ptr)this);
                return false;
            }
            SetTitle (titlePrefix&curName);
        }
        return true;
    }
    return false;
}

//__________________________________________________________

void    _HYModelWindow::DoPasteToCells ()
{
    _Formula f (clipboardString);
    _HYTable*t = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,4);

    if (AddFormulaParametersToList(f)) {
        _SimpleList sel, sel2;
        t->GetSelection (sel);

        SymmetrizeSelection (sel);

        _String * cStr = new _String (clipboardString);
        checkPointer (cStr);

        for (long k=0; k< sel.lLength; k++) {
            long h = sel.lData[k]/stateLabels.lLength,
                 v = sel.lData[k]%stateLabels.lLength;

            t->SetCellData(cStr,h,v,HY_TABLE_EDIT_TEXT|HY_TABLE_SELECTED,false);
            t->SetCellData(cStr,v,h,HY_TABLE_EDIT_TEXT|HY_TABLE_SELECTED,false);

            sel2 << v*stateLabels.lLength + h;
        }
        cStr->nInstances += 2*sel.lLength;
        DeleteObject (cStr);
        t->_MarkCellsForUpdate (sel);
        t->_MarkCellsForUpdate (sel2);
        taint = sel.lLength;
    }

}

//__________________________________________________________

void    _HYModelWindow::UpdateParameterName (long newID, long oldID)
{
    _HYTable*t = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,4);

    _List        processedStrings,
                 newStrings;

    _String*     e = new _String (empty);

    _SimpleList  flagged;

    for (long k = 0; k<stateLabels.lLength; k++)
        for (long j=k+1; j<stateLabels.lLength; j++) {
            _String * cStr = (_String*)t->GetCellData (j,k),
                      * pStr;
            if (cStr->sLength) {
                long    f = processedStrings.BinaryFind (cStr);
                if (f>=0) {
                    pStr = (_String*)newStrings(f);
                    if (pStr->sLength) {
                        flagged << k*stateLabels.lLength+j;
                        flagged << j*stateLabels.lLength+k;
                        t->SetCellData (pStr,k,j,t->cellTypes.lData[k*stateLabels.lLength+j],false);
                        t->SetCellData (pStr,j,k,t->cellTypes.lData[j*stateLabels.lLength+k],false);
                        pStr->nInstances+=2;
                    }
                } else {
                    bool t = false;
                    _Formula fla (*cStr);
                    _String  *newS;
                    for (f = 0; f<fla.GetList().lLength; f++) {
                        _Operation* o = (_Operation*)fla.GetList()(f);
                        if (o->GetAVariable()==oldID) {
                            o->SetAVariable (newID);
                            t = true;
                        }
                    }
                    if (t) {
                        newS = (_String*)fla.toStr();
                        newS -> nInstances = 0;
                        j--;
                    } else {
                        newS = e;
                    }
                    newStrings.InsertElement(newS,processedStrings.BinaryInsert (cStr),false);
                }
            }
        }
    t->_MarkCellsForUpdate (flagged);
    taint = flagged.lLength;
}

//__________________________________________________________

bool    warnDeleteParameter = false;

//__________________________________________________________

bool    _HYModelWindow::RemoveParameterName (long oldID)
{
    _HYTable*t = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,4);

    _List        processedStrings,
                 newStrings;

    _SimpleList  flagged;

    for (long k = 0; k<stateLabels.lLength; k++)
        for (long j=k+1; j<stateLabels.lLength; j++) {
            _String * cStr = (_String*)t->GetCellData (j,k),
                      * pStr;
            if (cStr->sLength) {
                long    f = processedStrings.BinaryFind (cStr);
                if (f>=0) {
                    pStr = (_String*)newStrings(f);
                    if (!pStr->Equal((_String*)processedStrings(f))) {
                        if ((flagged.lLength==0)&&(!warnDeleteParameter))
                            if (!ProceedPromptWithCheck (deleteParameters, donotWarnAgain,warnDeleteParameter, (Ptr)this)) {
                                return false;
                            }

                        flagged << k*stateLabels.lLength+j;
                        flagged << j*stateLabels.lLength+k;
                        t->SetCellData (pStr,k,j,t->cellTypes.lData[k*stateLabels.lLength+j],false);
                        t->SetCellData (pStr,j,k,t->cellTypes.lData[j*stateLabels.lLength+k],false);
                        pStr->nInstances+=2;
                    }
                } else {
                    _Formula fla (*cStr);
                    _String  newS;
                    for (f = 0; f<fla.GetList().lLength; f++) {
                        _Operation* o = (_Operation*)fla.GetList()(f);
                        if (o->GetAVariable()==oldID) {
                            break;
                        }
                    }
                    if (f<fla.GetList().lLength) {
                        newS = empty;
                        j--;
                    } else {
                        newS = *cStr;
                    }
                    newStrings.InsertElement(&newS,processedStrings.BinaryInsert (cStr),true);
                }
            }
        }
    t->_MarkCellsForUpdate (flagged);
    taint = flagged.lLength;
    return true;
}

//__________________________________________________________

void    _HYModelWindow::DoMultiplyCells ()
{
    _HYTable*       t = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,4);
    _HYButtonBar*  bb = (_HYButtonBar*)GetCellObject (MODEL_BUTTON_ROW,0);

    _List           menuChoices;

    long k = 0;

    _String         s ("User expression");
    menuChoices && & s;
    if (localVariables.lLength+globalVariables.lLength+categoryVariables.lLength) {
        menuChoices && & menuSeparator;
    }

    for (k=0; k<localVariables.lLength; k++) {
        menuChoices << localVariables(k);
    }
    if (localVariables.lLength&&(globalVariables.lLength+categoryVariables.lLength)) {
        menuChoices && & menuSeparator;
    }
    for (k=0; k<globalVariables.lLength; k++) {
        menuChoices << globalVariables(k);
    }
    if (categoryVariables.lLength&&(globalVariables.lLength+localVariables.lLength)) {
        menuChoices && & menuSeparator;
    }
    for (k=0; k<categoryVariables.lLength; k++) {
        menuChoices << categoryVariables(k);
    }

    int h,v;

    bb->GetButtonLoc  (4,h,v,true);
    bb->_UnpushButton ();

    s = HandlePullDown (menuChoices,h,v,0);
    h = menuChoices.Find (&s);

    if (h>=0) {
        _String * t1 = (_String*)menuChoices (h);
        if (h==0) {
            _String prompt ("Your expression:");
            s = "1";
            if (!EnterStringDialog(s,prompt, (Ptr)this)) {
                return;
            }
            _Formula fl (s);
            if (fl.IsEmpty()) {
                s = "You entered an invalid formula";
                WarnError (s);
                return;
            }
            _PMathObj  val = fl.Compute();
            if (!val || (val->ObjectClass()!=NUMBER)) {
                s = "Your expression did not evaluate to a number";
                WarnError (s);
                return;
            }
            fl.SimplifyConstants();
            if (!AddFormulaParametersToList (fl)) {
                return;
            }
            t1 = (_String*)fl.toStr();

            if (fl.GetList().lLength > 1) {
                *t1 = _String ('(') & *t1 & ')';
            }
        }

        _SimpleList sel,
                    sel2;

        t->GetSelection (sel);

        SymmetrizeSelection (sel);

        _List       doneStrings,
                    doneStrings2;

        for (long k=0; k< sel.lLength; k++) {
            long h = sel.lData[k]/stateLabels.lLength,
                 v = sel.lData[k]%stateLabels.lLength;

            _String* cStr = (_String*)t->GetCellData(v,h)->makeDynamic();

            if (cStr->sLength) {
                long ff = doneStrings.BinaryFind (cStr);
                if (ff>=0) {
                    DeleteObject (cStr);
                    cStr = (_String*)doneStrings2 (ff);
                    cStr->nInstances++;
                } else {
                    _Formula fla (*cStr);
                    _String  *cStr2 = (_String*)fla.toStr();

                    *cStr2 = _String ('(') & *cStr2 & ")*" & *t1;
                    _Formula fla2(*cStr2);
                    fla2.SimplifyConstants ();
                    DeleteObject (cStr2);

                    cStr2 = (_String*)fla2.toStr();
                    doneStrings.BinaryInsert (cStr);
                    ff = doneStrings.BinaryFind (cStr);
                    doneStrings2.InsertElement (cStr2,ff,false);
                    cStr = cStr2;
                    cStr->nInstances++;
                }
            } else {
                *cStr = *t1;
                cStr->Trim (1,cStr->sLength-2);
            }

            t->SetCellData(cStr,h,v,HY_TABLE_EDIT_TEXT|HY_TABLE_SELECTED,false);
            t->SetCellData(cStr,v,h,HY_TABLE_EDIT_TEXT|HY_TABLE_SELECTED,false);

            sel2 << v*stateLabels.lLength + h;
            cStr->nInstances++;
        }

        if (h==0) {
            DeleteObject (t1);
        }

        t->_MarkCellsForUpdate (sel);
        t->_MarkCellsForUpdate (sel2);
        taint = sel.lLength;
    }
}

//__________________________________________________________

void    _HYModelWindow::DoClearCells ()
{
    _HYTable*t = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,4);

    _SimpleList sel, sel2;
    t->GetSelection (sel);

    for (long k=0; k< sel.lLength; k++) {
        long h = sel.lData[k]/stateLabels.lLength,
             v = sel.lData[k]%stateLabels.lLength;

        t->SetCellData(&empty,h,v,HY_TABLE_EDIT_TEXT|HY_TABLE_SELECTED,true);
        t->SetCellData(&empty,v,h,HY_TABLE_EDIT_TEXT|HY_TABLE_SELECTED,true);

        sel2 << v*stateLabels.lLength + h;
    }

    t->_MarkCellsForUpdate (sel);
    t->_MarkCellsForUpdate (sel2);
    taint = sel.lLength;
}

//__________________________________________________________

bool    _HYModelWindow::CopyMatrix (_String& modelMatrix)
{
    long       k = LocateVarByName (modelMatrix);
    if (k>=0) {
        _Variable* theMx = FetchVar(k);
        if (theMx->ObjectClass()==MATRIX) {
            _Matrix * mx = (_Matrix*)theMx->GetValue();
            return CopySimpleMatrix (mx);
        }
    }
    return false;
}

//__________________________________________________________

bool    _HYModelWindow::CopySimpleMatrix (_Matrix* mx)
{
    long       k;
    if ((mx->GetHDim()==stateLabels.lLength)&&
            (mx->GetVDim()==stateLabels.lLength)) {
        _HYTable* t = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,4);
        if (mx->MatrixType()==2) {
            _List   alreadyDone;
            bool    warnNewParametersSave = warnNewParameters;

            warnNewParameters = true;

            for (k=0; k<stateLabels.lLength; k++)
                for (long j=k+1; j<stateLabels.lLength; j++) {
                    _Formula* ff = mx->GetFormula (k,j);
                    if (ff) {
                        ff->SimplifyConstants();
                        _String* fs = (_String*)ff->toStr();
                        long f = alreadyDone.BinaryFind (fs);
                        if (f<0) {
                            AddFormulaParametersToList (*ff,true);
                            alreadyDone.BinaryInsert (fs);
                            f = alreadyDone.BinaryFind (fs);
                        }
                        DeleteObject (fs);
                        fs = (_String*)alreadyDone(f);
                        t->SetCellData(fs,k,j,HY_TABLE_EDIT_TEXT,false);
                        t->SetCellData(fs,j,k,HY_TABLE_EDIT_TEXT,false);
                        fs->nInstances+=2;
                    } else {
                        t->SetCellData(&empty,k,j,HY_TABLE_EDIT_TEXT,true);
                        t->SetCellData(&empty,j,k,HY_TABLE_EDIT_TEXT,true);
                    }
                }

            warnNewParameters = warnNewParametersSave;
            DoFitColumns();
            return true;
        } else {
            if (mx->MatrixType()==1) {
                for (k=0; k<stateLabels.lLength; k++)
                    for (long j=0; j<stateLabels.lLength; j++) {
                        if (j==k) {
                            continue;
                        }

                        _Parameter cellValue  = (*mx)(k,j);
                        _String    cellString (cellValue);

                        t->SetCellData(&cellString,k,j,HY_TABLE_EDIT_TEXT,true);
                    }
                DoFitColumns();
                return true;
            }
        }
    }
    return false;
}

//__________________________________________________________

void    _HYModelWindow::DoFitColumns ()
{

    _HYTable     * table  = (_HYTable*) GetCellObject (MODEL_MATRIX_ROW,4),
                   * columnHeaders = (_HYTable*) GetCellObject (MODEL_MATRIX_HEADER_ROW,4);

    table->AutoFitWidth(*columnHeaders);
    table->EnforceWidth         (table->rel.right-table->rel.left-HY_SCROLLER_WIDTH+1,stateLabels.lLength-1);
    columnHeaders->EnforceWidth (table->rel.right-table->rel.left-HY_SCROLLER_WIDTH+1,stateLabels.lLength-1);
    columnHeaders->SetColumnSpacing (stateLabels.lLength,HY_SCROLLER_WIDTH-columnHeaders->GetColumnSpacing (stateLabels.lLength),false);
    table->SetVisibleSize(table->rel);
    table->_MarkForUpdate();
    columnHeaders->_MarkForUpdate();
    dim = MinMaxWindowDimensions();
}

//__________________________________________________________

bool    _HYModelWindow::AddFormulaParametersToList (_Formula& f, bool respectClass)
{
    _SimpleList  newVars, classes;
    _List        undefs;

    {
        _AVLList na (&newVars);
        f.ScanFForVariables (na,true,false,false);
        na.ReorderList();
    }

    for (long k=0; k<newVars.lLength; k++) {
        _Variable * thisV = LocateVar (newVars.lData[k]);
        _String* vName = thisV->GetName();
        if (FindVarID (*vName) == -1) {
            undefs << vName;
            if (thisV->IsCategory()) {
                classes<<2;
            } else if (thisV->IsGlobal()) {
                classes<<1;
            } else {
                classes<<0;
            }
        }
    }

    if (undefs.lLength) {
        if (!warnNewParameters)
            if (!ProceedPromptWithCheck (spawnParameters, donotWarnAgain, warnNewParameters, (Ptr)this)) {
                return false;
            }

        if (respectClass)
            for (long k=0; k<undefs.lLength; k++) {
                InsertVarID (*(_String*)undefs(k),classes.lData[k]);
            }
        else
            for (long k=0; k<undefs.lLength; k++) {
                InsertVarID (*(_String*)undefs(k),0);
            }

        PrepareParameterMenu ();
        SetParameterMenuState (lastParameterChoice);
    }
    return true;
}

//__________________________________________________________

void    _HYModelWindow::DoSelectCells ()
{
    _HYTable*       t = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,4);
    _HYButtonBar*  bb = (_HYButtonBar*)GetCellObject (MODEL_BUTTON_ROW,0);

    _SimpleList     sel;
    _List           menuChoices;
    _String         s ("All Cells"),
                    *match;

    t->GetSelection (sel);

    menuChoices && & s;
    s = "All Empty Cells";
    menuChoices && & s;
    s = "Invert Selection";
    menuChoices && & s;
    s = "Scroll To Selection";
    menuChoices && & s;

    if (sel.lLength==1) {
        match = (_String*)t->GetCellData (sel.lData[0]%stateLabels.lLength,sel.lData[0]/stateLabels.lLength);
        if (match->sLength) {
            s = _String("Cells equal to \"") & *match & '"';
            menuChoices && &s;
        }
    }

    if (sel.lLength) {
        s = _String("Symmetrize Selection");
        menuChoices && &s;
    }

    int h,v;

    bb->GetButtonLoc  (3,h,v,true);
    bb->_UnpushButton ();
    s = HandlePullDown (menuChoices,h,v,0);

    h = menuChoices.Find (&s);

    switch (h) {
    case 0:
        DoSelectAll();
        break;
    case 1:
        DoSelectAllEmpty();
        break;
    case 2:
        t->InvertSelection();
        t->_MarkForUpdate();
        break;
    case 3:
        DoScrollToSelection(sel);
        break;
    default: {
        if (h>=0) {
            _String s2 ("Symmetrize Selection");
            if (s.Equal (&s2)) {
                SymmetrizeSelection (sel,true);
                t->SetSelection     (sel,true);
                t->_MarkForUpdate   ();
            } else {
                DoMatchCells(match);
            }
        }
    }
    }
}
//__________________________________________________________

void    _HYModelWindow::DoScrollToSelection (_SimpleList& sel)
{
    if (sel.lLength) {
        _HYTable*       t = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,4);
        t->ScrollToRow    (sel.lData[0]/t->horizontalSpaces.lLength);
        t->ScrollToColumn (sel.lData[0]%t->horizontalSpaces.lLength);
        t->_MarkForUpdate();
    }
}
//__________________________________________________________

void    _HYModelWindow::DoSelectAll ()
{
    _HYTable*       t = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,4);

    _SimpleList     sel (stateLabels.lLength*stateLabels.lLength);

    for (long k=0; k<stateLabels.lLength*stateLabels.lLength; k++) {
        sel << k;
    }

    t->ClearSelection (false);
    t->SetSelection (sel,true);
    t->_MarkForUpdate();
}

//__________________________________________________________

void    _HYModelWindow::DoSave (char mode)
{
    _String pathToModelTemplates,
            ending,
            *varName,
            EFVBit,
            rateBit;

    long    modelType;

    if (GetTitle().Cut(titlePrefix.sLength,-1).beginswith(untitled))
        if (!DoEditModelName ()) {
            return;
        }

    FILE    * F;
    _Matrix * mods = nil;

    _HYTable* t = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,4);

    long    k,j,m;


    if (EFVChoice != 0) {

        pathToModelTemplates = libDirectory&"SubstitutionClasses"&directorySep;

        switch (type) {
        case 0:
            ending = "NucEFV";
            break;
        case 1:
            ending = "AAEFV";
            break;
        case 2:
            ending = "DinucEFV";
            break;
        case 3:
            ending = "CodonEFV";
            break;
        }
        pathToModelTemplates = pathToModelTemplates & ending & directorySep & *(_String*)EFVOptions(EFVChoice);

        F = doFileOpen (pathToModelTemplates.sData,"rb");

        if (!F) {
            ending = _String ("Can't open the batch code for \"") & *(_String*)EFVOptions(EFVChoice) & "\". That file was recently moved or deleted. Bad idea.";
            EFVOptions.Delete (EFVChoice);
            EFVChoice = 0;
            WarnError (ending);
            return;
        } else {
            _String buffer (F);
            fclose (F);
            buffer = buffer.Replace (ModelMatrixDimension,_String ((long)stateLabels.lLength),true);
            long f = buffer.Find (EFVCodeMarker);
            if (f>=0) {
                EFVBit = buffer.Cut (0,f-1);
                buffer.Trim (f+EFVCodeMarker.sLength,-1);
                f = batchLanguageFunctionNames.lLength;
                if (type == HY_MODEL_TYPE_CODON) {
                    _String * gcs = DefStringForGCode (&genCode);
                    _ExecutionList ex1 (*gcs);
                    ex1.Execute();
                    DeleteObject (gcs);
                }

                _ExecutionList   thisList;
                terminateExecution = false;
                thisList.BuildList (buffer);
                thisList.Execute();
                if (terminateExecution==false) {
                    k = LocateVarByName (EFVModifier);
                    if (k>=0) {
                        _Variable* thisVar = FetchVar (k);
                        if (thisVar->ObjectClass() == MATRIX) {
                            mods = (_Matrix*)thisVar->GetValue();
                            if (mods->MatrixType()!=2) {
                                mods = nil;
                            }
                        }
                    }
                } else {
                    terminateExecution = false;
                }

                while (f<batchLanguageFunctionNames.lLength) {
                    batchLanguageFunctionNames.Delete (f);
                    batchLanguageFunctionParameters.Delete (f);
                    batchLanguageFunctions.Delete(f);
                    batchLanguageFunctionParameterLists.Delete(f);
                    batchLanguageFunctionClassification.Delete(f);

                }
            } else {
                EFVBit = buffer;
            }
        }
    } else {
        _HYTable* EFV = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,1);
        _String *cellValue = (_String*)EFV->GetCellData (0,stateLabels.lLength-1);
        if (cellValue->Equal(&invalidEFV)) {
            ending = "This model can't be saved because user equilibrium frequencies don't add up to 1.";
            ProblemReport (ending,(Ptr)this);
            return;
        }
        _String uEFV (256L, true);
        uEFV << &modelEFVVector;
        uEFV << "={{";
        for (k=0; k<stateLabels.lLength; k++) {
            cellValue = (_String*)EFV->GetCellData (0,k);
            uEFV << cellValue;
            if (k<stateLabels.lLength-1) {
                uEFV << ',';
            }
        }
        uEFV << "}};\n";
        uEFV << &multiplyByFrequencies;
        uEFV << "=1;\n";
        uEFV.Finalize();
        EFVBit = uEFV;
    }

    CheckDependencies ();

    if (localVariables.lLength+globalVariables.lLength+categoryVariables.lLength==0) {
        ending = "A model without estimable parameters isn't very interesting. Please define some parameters before saving.";
        ProblemReport (ending,(Ptr)this);
        return;
    }

    if (categoryVariables.lLength) {
        if (rateChoice) {
            pathToModelTemplates = libDirectory&"SubstitutionClasses"&directorySep&"Heterogeneity"&directorySep& *(_String*)rateOptions(rateChoice);

            F = doFileOpen (pathToModelTemplates.sData,"rb");

            if (!F) {
                ending = _String ("Can't open the batch code for \"") & *(_String*)rateOptions(rateChoice) & "\". That file was recently moved or deleted. Bad idea.";
                rateOptions.Delete (rateChoice);
                rateChoice = 0;
                WarnError (ending);
                return;
            } else {
                _String buffer (F);
                fclose (F);
                rateBit = buffer;
            }
        } else {
            rateBit = "global shapeParameter = .5;\nshapeParameter:>0.01;shapeParameter:<100;\ncategory     categoryPlaceholder = (rateClassCount, EQUAL, MEAN, GammaDist(_x_,shapeParameter,shapeParameter), CGammaDist(_x_,shapeParameter,shapeParameter), 0 ,1e25,CGammaDist(_x_,shapeParameter+1,shapeParameter));";
        }
    }


    _String modelOut (1024L, true),
            modelID  (GetTitle().Cut(titlePrefix.sLength,-1));

    modelOut << &modelName;
    modelOut << '=';
    modelOut << '"';
    modelOut << &modelID;
    modelOut << '"';
    modelOut << ';';
    modelOut << '\n';

    modelOut << &modelDimension;
    modelOut << '=';
    modelOut << _String ((long)stateLabels.lLength);
    modelOut << ';';
    modelOut << '\n';

    modelOut << &modelOptions;
    modelOut << '=';
    if (categoryVariables.lLength) {
        modelType = HY_DATAPANEL_MODEL_MODELS+HY_DATAPANEL_MODEL_GLOBALG ;
    } else if (globalVariables.lLength) {
        modelType = HY_DATAPANEL_MODEL_MODELS+HY_DATAPANEL_MODEL_GLOBAL;
    } else {
        modelType = HY_DATAPANEL_MODEL_MODELS;
    }

    modelOut.AppendNewInstance(new _String (modelType)) ;
    modelOut << ";\n";

    modelOut << &ModelEFVType;
    modelOut << "=\"";
    modelOut << (_String*)EFVOptions (EFVChoice);
    modelOut << "\";\n";

    modelOut << &modelDataType;
    modelOut << '=';
    modelOut.AppendNewInstance(new _String ((long)type));
    modelOut << ";\n";

    if (type == HY_MODEL_TYPE_CODON) {
        varName = DefStringForGCode (&genCode);
        modelOut << varName;
        DeleteObject (varName);
    }

    modelOut << "\n\n\n";
    modelOut << "function ";
    modelOut << modelFunction;
    modelOut << "(";
    modelOut << &modelMatrixName;
    modelOut << "&,EFV)\n{\n\n\t";

    modelOut << & modelMatrixName;
    modelOut << "={";
    ending = (long)stateLabels.lLength;
    modelOut << &ending;
    modelOut << ',';
    modelOut << &ending;
    modelOut << "};\n";

    _List     vnames,
              pvnames,
              binames;

    for (k=0; k<globalVariables.lLength; k++) {
        modelOut << "\tglobal ";
        varName = (_String*)globalVariables(k);
        if (!varName->beginswith (globalPrefix)) {
            _String vn;
            vn = globalPrefix & '_' & *varName;
            vnames << varName;
            pvnames && & vn;
            modelOut << &vn;
        } else {
            modelOut << varName;
        }
        modelOut << ";\n";
    }

    for (k=0; k<categoryVariables.lLength; k++) {
        varName = (_String*)categoryVariables(k);
        if (!varName->beginswith (categoryPrefix)) {
            _String vn;
            vn = categoryPrefix & '_' & *varName;
            vnames << varName;
            pvnames && & vn;
            ending = vn;
        } else {
            ending = *varName;
        }
        _String cDef (rateBit.Replace ("categoryPlaceholder",ending,true));
        if (k) {
            cDef = cDef.Replace ("shapeParameter",_String("shapeParameter")&_String(k),true);
        }
        modelOut << &cDef;
    }

    binames && & vnames;
    binames && & pvnames;

    vnames.Clear();
    pvnames.Clear();

    if (mods)
        for (k=0; k<stateLabels.lLength; k++) {
            _Constant r (k);
            for (j=0; j<stateLabels.lLength; j++) {
                if (j==k) {
                    continue;
                }
                varName = (_String*)t->GetCellData (j,k);
                if (varName->sLength) {
                    _Constant c (j);
                    _PMathObj mod = mods->MAccess (&r,&c);

                    m = vnames.BinaryFind (varName);
                    if (m<0) {
                        _Formula fl (*varName);
                        _String  *ss = (_String*) fl.toStr (&binames);
                        vnames.BinaryInsert (varName);
                        m = vnames.BinaryFind (varName);
                        pvnames.InsertElement (ss,m,false);
                    }
                    varName = (_String*)pvnames (m);

                    _String index (k),
                            index2(j);

                    modelOut << '\t';
                    modelOut << &modelMatrixName;
                    modelOut << '[';
                    modelOut << & index;
                    modelOut << "][";
                    modelOut << & index2;
                    modelOut << "]:=";
                    modelOut << varName;
                    if (mod->ObjectClass ()==STRING) {
                        ending = _String( "*(") & *((_FString*)mod)->theString & ')';
                        modelOut << &ending;
                    }
                    modelOut << ";\n";
                    DeleteObject (mod);
                }
            }
        }
    else
        for (k=0; k<stateLabels.lLength; k++)
            for (j=k+1; j<stateLabels.lLength; j++) {
                varName = (_String*)t->GetCellData (j,k);
                if (varName->sLength) {
                    m = vnames.BinaryFind (varName);
                    if (m<0) {
                        _Formula fl (*varName);
                        _String  *ss = (_String*) fl.toStr (&binames);
                        vnames.BinaryInsert (varName);
                        m = vnames.BinaryFind (varName);
                        pvnames.InsertElement (ss,m,false);
                    }
                    varName = (_String*)pvnames (m);
                    _String index (k),
                            index2(j);

                    modelOut << '\t';
                    modelOut << &modelMatrixName;
                    modelOut << '[';
                    modelOut << & index;
                    modelOut << "][";
                    modelOut << & index2;
                    modelOut << "]:=";
                    modelOut << varName;
                    modelOut << ";\n\t";
                    modelOut << &modelMatrixName;
                    modelOut << '[';
                    modelOut << & index2;
                    modelOut << "][";
                    modelOut << & index;
                    modelOut << "]:=";
                    modelOut << varName;
                    modelOut << ";\n";
                }
            }

    ending = "\n\treturn 0;\n}\n\n";
    modelOut << &ending;

    modelOut << &EFVBit;

    modelOut.Finalize();

    if (mode==-1) {
        pathToModelTemplates = libDirectory&"SubstitutionModels"&directorySep&"User"&directorySep;

        switch (type) {
        case 0:
            ending = "Nucleotide";
            break;
        case 1:
            ending = "AA";
            break;
        case 2:
            ending = "Dinucleotide";
            break;
        case 3:
            ending = "Codon";
            break;
        }
        pathToModelTemplates = pathToModelTemplates & ending & directorySep & modelID;
    } else {
        pathToModelTemplates = WriteFileDialogInput ();
    }

    if (pathToModelTemplates.sLength) {
        F = doFileOpen (pathToModelTemplates.sData,"w");

        if (!F) {
            pathToModelTemplates = _String("Can't open \"") & pathToModelTemplates & "\" for writing. Models must be saved in \"Substituion Models/User/\" to be accessible to HyPhy GUI.";
            WarnError (pathToModelTemplates);
            return;
        } else {
            fwrite (modelOut.sData,modelOut.sLength,1,F);
        }

        _List*  modelList = FindModelTemplate (&modelID);

        if (modelList) {
            modelList->Replace (2,&pathToModelTemplates,true);
            _SimpleList * parms = (_SimpleList*)(*modelList)(1);
            parms->lData[1] = stateLabels.lLength;
            parms->lData[0] = modelType;
        } else {
            _List newModel;
            newModel && & modelID;
            _SimpleList parms;
            parms << modelType;
            parms << stateLabels.lLength;
            newModel && & parms;
            newModel && & pathToModelTemplates;
            modelTemplates && & newModel;
        }

        taint = false;
        fclose (F);
    }

    if (mods) {
        DeleteVariable (EFVModifier);
    }
}

//__________________________________________________________

void    _HYModelWindow::DoSelectAllEmpty ()
{
    _HYTable*       t = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,4);

    _SimpleList     sel (stateLabels.lLength*stateLabels.lLength);

    for (long k=0; k<stateLabels.lLength-1; k++)
        for (long n=k+1; n<stateLabels.lLength; n++) {
            _String * v = (_String*)t->GetCellData(n,k);
            if (v->sLength==0) {
                sel << k*stateLabels.lLength+n;
                sel << n*stateLabels.lLength+k;
            }
        }

    t->ClearSelection (false);
    t->SetSelection (sel,true);
    t->_MarkForUpdate();
}

//__________________________________________________________

void    _HYModelWindow::DoMatchCells (_String* m)
{
    _HYTable*       t = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,4);

    _SimpleList     sel (stateLabels.lLength*stateLabels.lLength);

    for (long k=0; k<stateLabels.lLength-1; k++)
        for (long n=k+1; n<stateLabels.lLength; n++) {
            _String * v = (_String*)t->GetCellData(n,k);
            if (v->Equal (m)) {
                sel << k*stateLabels.lLength+n;
                sel << n*stateLabels.lLength+k;
            }
        }

    t->ClearSelection (false);
    t->SetSelection   (sel,true);
    t->_MarkForUpdate ();
}

//__________________________________________________________

void    _HYModelWindow::Activate (void)
{
    //_UpdateEditMenu   ();
    _HYTWindow::Activate ();
    _CheckClipboard();
    UpdateButtonState ();
}

//__________________________________________________________

void    _HYModelWindow::Deactivate (void)
{
    _HYTWindow::Deactivate ();
    _SetClipboard();
}

//__________________________________________________________

void    _HYModelWindow::SyncEditBox (void)
{
    _HYTextBox   *tl = (_HYTextBox*)GetCellObject (MODEL_BUTTON_ROW,5);
    _HYButtonBar *bb3= (_HYButtonBar*)GetCellObject (MODEL_BUTTON_ROW,6);
    tl->SetForeColor    (editBoxGood);
    tl->SetText         (clipboardString);
    tl->UnfocusComponent();
    bb3->EnableButton   (0,false);
}

//__________________________________________________________

void    _HYModelWindow::SymmetrizeSelection (_SimpleList& sel, bool pad)
{
    if (sel.lLength) {
        _SimpleList symSel;
        long k,h,v;
        for (k=0; k< sel.lLength; k++) {
            h = sel.lData[k]/stateLabels.lLength;
            v = sel.lData[k]%stateLabels.lLength;

            if (h<v) {
                symSel << sel.lData[k];
            } else {
                symSel << v*stateLabels.lLength+h;
            }
        }

        symSel.Sort();

        sel.Clear();

        sel << symSel.lData[0];

        if (pad) {
            h = symSel.lData[0]/stateLabels.lLength;
            v = symSel.lData[0]%stateLabels.lLength;
            sel<<v*stateLabels.lLength+h;
            for (k=1; k<symSel.lLength; k++) {
                if (symSel.lData[k]>symSel.lData[k-1]) {
                    h = symSel.lData[k]/stateLabels.lLength;
                    v = symSel.lData[k]%stateLabels.lLength;
                    sel<<symSel.lData[k];
                    sel<<v*stateLabels.lLength+h;
                }
            }
        } else
            for (k=1; k<symSel.lLength; k++) {
                if (symSel.lData[k]>symSel.lData[k-1]) {
                    sel<<symSel.lData[k];
                }
            }
    }
}



//__________________________________________________________

void    _HYModelWindow::UpdateClassMenus (void)
{
    _HYPullDown * p1 = (_HYPullDown*)GetCellObject (MODEL_CLASS_ROW,2);

    p1->DeleteAllItems();

    long k = 0, p;

    for (; k<subClasses.lLength; k++) {
        _List * sC = (_List*)subClasses(k);
        p1->AddMenuItem (*(_String*)(*sC)(0),-1);
    }

    if (k&&matrixTemplates.lLength) {
        p1->AddMenuItem (menuSeparator,-1);
    }

    p = p1->MenuItemCount();

    for (k=0; k<matrixTemplates.lLength; k++) {
        _List * sC = (_List*)matrixTemplates(k);
        p1->AddMenuItem (*(_String*)(*sC)(0),-1);
        p1->_SetMenuItemTextStyle (k+p, HY_FONT_ITALIC);
    }

    p+=k;

    if (!p) {
        p1->AddMenuItem (noneDefined,-1);
        //p1->EnableMenu (false);
    }

    p1->ChangeSelection (0,true);
}

//__________________________________________________________

void    _HYModelWindow::ProcessClassMenu (void)
{
    bool          f = false;

    _HYPullDown * p1 = (_HYPullDown*)GetCellObject (MODEL_CLASS_ROW,2),
                  * p2 = (_HYPullDown*)GetCellObject (MODEL_CLASS_ROW,5);

    _HYButtonBar* bb2 = (_HYButtonBar*)GetCellObject (MODEL_BUTTON_ROW,4);

    p2->DeleteAllItems ();

    if (subClasses.lLength + matrixTemplates.lLength) {
        if (p1->GetSelection() < subClasses.lLength)
            for (long k = 0; k<subClasses.lLength; k++) {
                _List * sC = (_List*)subClasses(k);
                p2->AddMenuItem (*(_String*)(*sC)(0),-1);
            }
        else {
            p2->AddMenuItem (noneApplicable,-1);
            p2->EnableItem  (0,false);
        }

        p2->ChangeSelection (0,true);
        f = true;
    } else {
        p2->AddMenuItem (noneApplicable,-1);
        p2->EnableItem  (0,false);
    }

    p2->ChangeSelection (0,false);

    bb2->EnableButton (0,f);
    bb2->EnableButton (1,f);
    bb2->EnableButton (2,f);
    bb2->EnableButton (3,f);
    bb2->EnableButton (4,f);
}

//__________________________________________________________

void    _HYModelWindow::ProcessTemplateSelection (long k)
{

    _HYPullDown * p1 = (_HYPullDown*)GetCellObject (MODEL_CLASS_ROW,2),
                  * p2 = (_HYPullDown*)GetCellObject (MODEL_CLASS_ROW,5);

    _HYTable    * table = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,4);

    long        s1 = p1->GetSelection(),
                s2 = p2->GetSelection();

    _SimpleList gSel;
    if (s1<subClasses.lLength) {
        _SimpleList *tSel1 =  (_SimpleList*) (*(_List*)subClasses(s1))(1),
                     *tSel2 =  (_SimpleList*) (*(_List*)subClasses(s2))(1);

        s1 = table->horizontalSpaces.lLength;

        for (long i=0; i<tSel1->lLength; i++) {
            s2 = tSel1->lData[i];
            for (long j=0; j<tSel2->lLength; j++) {
                gSel << s1*s2+tSel2->lData[j];
            }
        }
    } else {
        if (subClasses.lLength) {
            s1 -= subClasses.lLength+1;
        }

        _SimpleList* tSel = (_SimpleList*) (*(_List*)matrixTemplates(s1))(1);
        gSel.Duplicate (tSel);
    }


    if (k==0) {
        table->ClearSelection();
        table->SetSelection (gSel,true);
    } else {
        _SimpleList sel, fSel;
        table->GetSelection (sel);
        switch (k) {
        case 1:
            fSel.Union (sel,gSel);
            break;
        case 2:
            fSel.Intersect (sel,gSel);
            break;
        case 3:
            fSel.XOR (sel,gSel);
            break;
        case 4:
            fSel.Subtract (sel,gSel);
            break;
        }
        table->ClearSelection();
        table->SetSelection (fSel,true);

    }

    table->_MarkForUpdate();
}

//__________________________________________________________

void    _HYModelWindow::BuildTemplates (_SimpleList* geneticCode, char type)
{
    _List   receptacle;
    _String pathToModelTemplates;
    long    k;

    pathToModelTemplates = libDirectory&"SubstitutionClasses";

    _String ending;

    switch (type) {
    case HY_MODEL_TYPE_NUC:
        ending = "nuc.bf";
        break;
    case HY_MODEL_TYPE_AA:
        ending = "aa.bf";
        break;
    case HY_MODEL_TYPE_DINUC:
        ending = "dinuc.bf";
        break;
    case HY_MODEL_TYPE_CODON:
        ending = "codon.bf";
        break;
    }

    ending = directorySep & ending;

    ScanDirectoryForFileNames (pathToModelTemplates,receptacle,false);
    for (k=0; k<receptacle.lLength; k++)
        if (((_String*)receptacle(k))->endswith (ending)) {
            break;
        }

    if (k<receptacle.lLength) {
        FILE * thisFile = doFileOpen (((_String*)receptacle(k))->sData,"rb");
        if (thisFile) {
            _String buffer (thisFile);
            fclose (thisFile);
            long   g = batchLanguageFunctionNames.lLength;
            _ExecutionList   thisList;
            terminateExecution = false;
            thisList.BuildList (buffer);
            thisList.Execute();
            if (!terminateExecution) {
                long     tempFunc  = batchLanguageFunctionNames.Find (&buildTemps),
                         classFunc = batchLanguageFunctionNames.Find (&buildClasses),
                         tempNames = LocateVarByName (nameMatrixT),
                         classNames= LocateVarByName (nameMatrixS);

                if (((classNames>=0)||(tempNames>=0))&&(((classNames>=0)&&(classFunc>=0))||
                                                        ((tempNames>=0)&&(tempFunc>=0)))) {
                    if (geneticCode) {
                        genCode.Duplicate (geneticCode);
                        _String buffer2 (256,true);
                        if (type == HY_MODEL_TYPE_CODON) {
                            _String * gcs = DefStringForGCode (geneticCode);
                            _ExecutionList ex1 (*gcs);
                            ex1.Execute();
                            DeleteObject (gcs);
                        }
                    }

                    _Variable * v1,
                              * res;

                    _Matrix   * setNames,
                              * setDescription;

                    if ((tempNames>=0)&&(tempFunc>=0)) {
                        v1 = FetchVar (tempNames);
                        if (v1->ObjectClass() == MATRIX) {
                            setNames = (_Matrix*) v1->GetValue();
                            for (k=0; k<setNames->GetHDim()*setNames->GetVDim(); k++) {
                                buffer  = _String ("out=")&buildTemps&'('&k&");";
                                _ExecutionList  spawn (buffer);
                                spawn.Execute();
                                long  p = LocateVarByName (resultName);
                                if (p>=0) {
                                    res = FetchVar (p);
                                    if (res->ObjectClass() == MATRIX) {
                                        _SimpleList indices;
                                        setDescription = (_Matrix*)res->GetValue();

                                        for (p=0; p<setDescription->GetHDim()*setDescription->GetVDim(); p++) {
                                            long ii = setDescription->theData[p];
                                            if ((p>0)&&(ii==0)) {
                                                break;
                                            }
                                            indices << ii;
                                        }

                                        _List tPair;

                                        indices.Sort();
                                        tPair << ((_FString*)(((_Formula**)setNames->theData)[k])->Compute())->theString;
                                        tPair && & indices;

                                        matrixTemplates && & tPair;
                                    }
                                }
                            }
                        }
                    }
                    if ((classNames>=0)&&(classFunc>=0)) {
                        v1 = FetchVar (LocateVarByName (nameMatrixS));
                        if (v1->ObjectClass() == MATRIX) {
                            setNames = (_Matrix*) v1->GetValue();
                            for (k=0; k<setNames->GetHDim()*setNames->GetVDim(); k++) {
                                buffer  = _String ("out=")&buildClasses&'('&k&");";
                                _ExecutionList  spawn (buffer);
                                spawn.Execute();
                                long  p = LocateVarByName (resultName);
                                if (p>=0) {
                                    res = FetchVar (p);
                                    if (res->ObjectClass() == MATRIX) {
                                        _SimpleList indices;
                                        setDescription = (_Matrix*)res->GetValue();

                                        for (p=0; p<setDescription->GetHDim()*setDescription->GetVDim(); p++) {
                                            long ii = setDescription->theData[p];
                                            if ((p>0)&&(ii==0)) {
                                                break;
                                            }
                                            indices << ii;
                                        }

                                        _List tPair;

                                        indices.Sort();
                                        tPair << ((_FString*)(((_Formula**)setNames->theData)[k])->Compute())->theString;
                                        tPair && & indices;

                                        subClasses && & tPair;
                                    }
                                }
                            }
                        }
                    }
                }
                while (g<batchLanguageFunctionNames.lLength) {
                    batchLanguageFunctionNames.Delete (g);
                    batchLanguageFunctionParameters.Delete (g);
                    batchLanguageFunctions.Delete(g);
                    batchLanguageFunctionParameterLists.Delete(g);
                    batchLanguageFunctionClassification.Delete(g);
                }

            }
        }
    }

    terminateExecution = false;
    DeleteVariable (nameMatrixT);
    DeleteVariable (nameMatrixS);
    DeleteVariable (modelGenCode);

    UpdateClassMenus();
}

//__________________________________________________________

void    _HYModelWindow::GrabEFVs (char type)
{
    _List   receptacle;
    _String pathToModelTemplates = libDirectory&"SubstitutionClasses"&directorySep;
    long    k;

    _String ending;

    switch (type) {
    case 0:
        ending = "NucEFV";
        break;
    case 1:
        ending = "AAEFV";
        break;
    case 2:
        ending = "DinucEFV";
        break;
    case 3:
        ending = "CodonEFV";
        break;
    }

    pathToModelTemplates =  pathToModelTemplates & ending;
    ScanDirectoryForFileNames (pathToModelTemplates,receptacle,true);

    for (k=0; k<receptacle.lLength; k++) {
        pathToModelTemplates = *(_String*) receptacle (k);
        pathToModelTemplates.Trim (pathToModelTemplates.FindBackwards (directorySep,0,-1)+1,-1);

        EFVOptions && & pathToModelTemplates;
    }

    _HYPullDown * efp = (_HYPullDown*)GetCellObject (MODEL_PARAMETER_ROW,5);
    efp->DeleteAllItems();
    for (k=0; k<EFVOptions.lLength; k++) {
        efp->AddMenuItem (*(_String*)EFVOptions(k),-1);
    }

    efp->_MarkForUpdate();
}


//__________________________________________________________

void    _HYModelWindow::GrabRateVariation (void)
{
    _List   receptacle;
    _String pathToModelTemplates;
    long    k;

    pathToModelTemplates =   libDirectory&"SubstitutionClasses" & directorySep & "Heterogeneity";
    ScanDirectoryForFileNames (pathToModelTemplates,receptacle,true);

    for (k=0; k<receptacle.lLength; k++) {
        pathToModelTemplates = *(_String*) receptacle (k);
        pathToModelTemplates.Trim (pathToModelTemplates.FindBackwards (directorySep,0,-1)+1,-1);
        rateOptions && & pathToModelTemplates;
    }
}

//__________________________________________________________

long    _HYModelWindow::ParameterEditBox(long inID)
{
    long res = lastParameterChoice;
    _HYModelParameterDialog* sd = new _HYModelParameterDialog (this, inID, &res);
    sd->Activate();
    while (windowObjectRefs.Find ((long)sd)>=0) {
        handleGUI();
    }
    return res;
}

//__________________________________________________________

void    _HYModelWindow::CheckDependencies (void)
{
    _SimpleList declared;
    _List       alreadyDone;
    _HYTable*   tt = (_HYTable*)GetCellObject (MODEL_MATRIX_ROW,4);

    long        k;


    {
        _AVLList    da (&declared);

        for (k=0; k<stateLabels.lLength; k++)
            for (long j=k+1; j<stateLabels.lLength; j++) {
                _String * thisCell = (_String*)tt->GetCellData (k,j);
                long    f = alreadyDone.BinaryFind (thisCell);
                if (f<0) {
                    _Formula ff (*thisCell);
                    ff.ScanFForVariables (da, true, false,false);
                    alreadyDone.BinaryInsert (thisCell);
                }
            }

        da.ReorderList();
    }

    for (k=localVariables.lLength-1; k>=0; k--) {
        long t = LocateVarByName (*(_String*)localVariables (k));
        if (t>=0) {
#ifndef USE_AVL_NAMES
            t = variableReindex.lData[t];
#else
            t = variableNames.GetXtra(t);
#endif
            if (declared.BinaryFind(t)<0) {
                localVariables.Delete (k);
            }
        }
    }
    for (k=globalVariables.lLength-1; k>=0; k--) {
        long t = LocateVarByName (*(_String*)globalVariables (k));
        if (t>=0) {
#ifndef USE_AVL_NAMES
            t = variableReindex.lData[t];
#else
            t = variableNames.GetXtra(t);
#endif
            if (declared.BinaryFind(t)<0) {
                globalVariables.Delete (k);
            }
        }
    }
    for (k=categoryVariables.lLength-1; k>=0; k--) {
        long t = LocateVarByName (*(_String*)categoryVariables (k));
        if (t>=0) {
#ifndef USE_AVL_NAMES
            t = variableReindex.lData[t];
#else
            t = variableNames.GetXtra(t);
#endif
            if (declared.BinaryFind(t)<0) {
                categoryVariables.Delete (k);
            }
        }
    }

}

//__________________________________________________________

_HYModelWindowDialog::_HYModelWindowDialog (long* ln1, long* ln2, _String* str1, long* ln3):_HYTWindow ("New Model Window", false, true)
{

    //long          index,
    //              cellWidth;

    mDim  = ln1;
    gCode = ln2;
    pl3   = ln3;
    sp1   = str1;


    _HYRect         canvasSettings = {30,100,30,100,HY_COMPONENT_NO_SCROLL|HY_COMPONENT_TRANSP_BG};

    _HYLabel*       l1      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l2      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l3      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l4      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l5      = new _HYLabel (canvasSettings, GetOSWindowData());

    canvasSettings.left = canvasSettings.right = 200;

    _HYPullDown*    p1      = new _HYPullDown (canvasSettings, GetOSWindowData());
    _HYPullDown*    p2      = new _HYPullDown (canvasSettings, GetOSWindowData());
    _HYPullDown*    p3      = new _HYPullDown (canvasSettings, GetOSWindowData());
    _HYPullDown*    p4      = new _HYPullDown (canvasSettings, GetOSWindowData());

    canvasSettings.left = canvasSettings.right = 120;
    _HYButton*      b1      = new _HYButton (canvasSettings, GetOSWindowData());

    canvasSettings.left = canvasSettings.right = 81;
    _HYButton*      b2      = new _HYButton (canvasSettings, GetOSWindowData());

    _HYColor        bgc     = GetDialogBackgroundColor();

    p1->SetMessageRecipient (this);
    p2->SetMessageRecipient (this);
    p3->SetMessageRecipient (this);
    p4->SetMessageRecipient (this);
    b1->SetMessageRecipient (this);
    b2->SetMessageRecipient (this);

    AddObject (l1);
    AddObject (l2);
    AddObject (l3);
    AddObject (l4);
    AddObject (l5);

    AddObject (p1);
    AddObject (p2);
    AddObject (p3);
    AddObject (p4);

    AddObject (b1);
    AddObject (b2);


    SetTableDimensions (5,3);
    SetCell   (0,0,l1);
    SetCell   (1,0,l2);
    SetCell   (2,0,l3);
    SetCell   (3,0,l4);
    SetCell   (4,0,l5);

    SetCell   (0,1,p1);
    SetCell   (0,2,p1);
    SetCell   (1,1,p2);
    SetCell   (1,2,p2);
    SetCell   (2,1,p3);
    SetCell   (2,2,p3);
    SetCell   (3,1,p4);
    SetCell   (3,2,p4);

    SetCell   (4,1,b1);
    SetCell   (4,2,b2);


    _HYFont  labelFont;
    l1->SetBackColor (bgc);
    l2->SetBackColor (bgc);
    l3->SetBackColor (bgc);
    l4->SetBackColor (bgc);
    l5->SetBackColor (bgc);

    p1->SetBackColor (bgc);
    p2->SetBackColor (bgc);
    p3->SetBackColor (bgc);
    p4->SetBackColor (bgc);

    b1->SetBackColor (bgc);
    b2->SetBackColor (bgc);

#ifdef __MAC__
    labelFont.face  = "System Font";
    labelFont.size  = 12;
#endif
#ifdef __WINDOZE__
    labelFont.face = "MS Sans Serif";
    labelFont.size = 12;
#endif
#ifdef __HYPHY_GTK__
    labelFont.face = _HY_SANS_FONT;
    labelFont.size = 10;
#endif

    labelFont.style = HY_FONT_PLAIN;

    l1->SetFont (labelFont);
    l2->SetFont (labelFont);
    l3->SetFont (labelFont);
    l4->SetFont (labelFont);

    b1->SetFont (labelFont);
    b2->SetFont (labelFont);

    l1->SetAlignFlags (HY_ALIGN_LEFT);
    l2->SetAlignFlags (HY_ALIGN_LEFT);
    l3->SetAlignFlags (HY_ALIGN_LEFT);
    l4->SetAlignFlags (HY_ALIGN_LEFT);

    b1->SetAlignFlags (HY_ALIGN_RIGHT);

    p1->SetAlignFlags (HY_ALIGN_LEFT);
    p2->SetAlignFlags (HY_ALIGN_LEFT);
    p3->SetAlignFlags (HY_ALIGN_LEFT);
    p4->SetAlignFlags (HY_ALIGN_LEFT);

    l1->SetText       ("Data type:");
    l2->SetText       ("Genetic Code:");
    l3->SetText       ("Copy from:");
    l4->SetText       ("Using mode:");

    b1->SetText       ("  OK  ");
    b1->SetButtonKind (HY_BUTTON_OK);
    b2->SetText       (" Cancel ");
    b2->SetButtonKind (HY_BUTTON_CANCEL);

    p1->AddMenuItem   (nucDataType,-1);
    p1->AddMenuItem   (dinucDataType,-1);
    p1->AddMenuItem   (proteinDataType,-1);
    p1->AddMenuItem   (codonDataType,-1);

    for (long k=0; k<geneticCodes.lLength; k++) {
        _List * thisCode = (_List*)geneticCodes(k);
        p2->AddMenuItem (*(_String*)(*thisCode)(0), -1);
    }

    p3->AddMenuItem   ("A rather long one",-1);

    p4->AddMenuItem   (parameterOption[0],-1);
    p4->AddMenuItem   (parameterOption[1],-1);
    p4->AddMenuItem   (parameterOption[2],-1);

    p2->EnableMenu    (false);
    p4->EnableMenu    (false);

    DeleteObject (l1);
    DeleteObject (l2);
    DeleteObject (l3);
    DeleteObject (l4);
    DeleteObject (l5);
    DeleteObject (p1);
    DeleteObject (p2);
    DeleteObject (p3);
    DeleteObject (p4);
    DeleteObject (b1);
    DeleteObject (b2);

    _HYRect     dim = MinMaxWindowDimensions();

    *gCode = -1;
    *mDim   = 4;

    SetWindowRectangle (0,0,dim.bottom,dim.right);
    CenterWindow       (this);
    SwitchModelOptions   (0);
    SwitchModelSelection (0);
}

//__________________________________________________________

bool    _HYModelWindowDialog::ProcessEvent (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    _HYPullDown*pd;
    long        i,f,g,k;

    if (e->EventClass()==_hyMenuSelChangeEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        g = e->EventCode().Cut (f+1,-1).toNum();

        for (i=0; i<components.lLength; i++)
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }

        if (i==5) {
            pd = (_HYPullDown*)GetCellObject (1,1);
            pd->EnableMenu(g==3);
            if (SwitchModelOptions   (g)) {
                pd = (_HYPullDown*)GetCellObject (2,1);
                pd->ChangeSelection (0);
                pd->_MarkForUpdate();
            }
        } else if (i==6) {
            SwitchModelOptions (3);
        } else if (i==7) {
            SwitchModelSelection (g);
        }

        done = true;
    } else if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++)
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }

        if (i==9) { // OK
            _HYPullDown*     pd = (_HYPullDown*)GetCellObject (2,1);
            *sp1 =* pd->GetMenuItem (pd->GetSelection());
            pd = (_HYPullDown*)GetCellObject (3,1);
            *pl3 = pd->GetSelection();
            postWindowCloseEvent (GetID());
            //_HYModelWindow* newModel = SpawnNewModel();
            //newModel->BringToFront();
        } else if (i==10) { // Cancel
            *mDim = -1;
            postWindowCloseEvent (GetID());
        }

        done = true;
    }

    if (done) {
        DeleteObject (e);
        return true;
    }
    return _HYTWindow::ProcessEvent(e);
}

//__________________________________________________________

bool    _HYModelWindowDialog::SwitchModelOptions (char type)
{
    long modelDim1 = 0;
    _HYPullDown * pd;
    switch (type) {
    case 0:
        *mDim  = 4;
        break;
    case 1:
        *mDim = 16;
        break;
    case 2:
        *mDim = 20;
        break;
    case 3: {
        pd = (_HYPullDown*) GetCellObject (1,1);
        *gCode = pd->GetSelection();
        modelDim1 = 64;
        *mDim = 0;
        _SimpleList * code = (_SimpleList*)(*(_List*)geneticCodes (*gCode))(2);
        for (long k=0; k<code->lLength; k++)
            if (code->lData[k]!=10) {
                (*mDim)++;
            }
    }
    }

    pd = (_HYPullDown*) GetCellObject (2,1);
    pd->DeleteAllItems ();
    pd->AddMenuItem    (none,-1);

    for (long k=0; k<modelTemplates.lLength; k++) {
        _List * thisModel = (_List*)modelTemplates (k);
        _SimpleList * modO= (_SimpleList*) (*thisModel) (1);
        if ((modO->lData[1]==*mDim)||(modO->lData[1]==modelDim1)) {
            pd->AddMenuItem (*(_String*)(*thisModel)(0),-1);
        }
    }

    modelDim1 = pd->MenuItemCount();
    //pd = (_HYPullDown*) GetCellObject (MODEL_CLASS_ROW,2);
    //pd->EnableMenu (modelDim1>1);
    return  modelDim1;
}

//__________________________________________________________

void    _HYModelWindowDialog::SwitchModelSelection (long index)
{
    _HYPullDown * pd    = (_HYPullDown*) GetCellObject (2,1);

    _List * thisModel   = FindModelTemplate (pd->GetMenuItem(index));

    pd = (_HYPullDown*)GetCellObject (3,1);

    if (thisModel) {
        _SimpleList * modO  = (_SimpleList*) (*thisModel)   (1);
        long        options = modO->lData[0];

        pd->EnableMenu (true);
        pd->EnableItem (1,options&HY_DATAPANEL_MODEL_GLOBAL);
        pd->EnableItem (2,options&HY_DATAPANEL_MODEL_GLOBALG);

        pd->ChangeSelection (0,false);
    } else {
        pd->EnableMenu (false);
    }
}

//__________________________________________________________

void    _HYModelWindowDialog::SetModelChoice (_String* model)
{
    _List*  modelT = FindModelTemplate (model);
    if (modelT) {
        _HYPullDown* pd = (_HYPullDown*)GetCellObject (0,1);
        _SimpleList * mo = (_SimpleList*)(*modelT)(1);
        long        k = 3;
        switch (mo->lData[1]) {
        case 4:
            k=0;
            break;
        case 16:
            k=1;
            break;
        case 20:
            k=2;
            break;
        }
        pd->ChangeSelection (k,true);
        pd = (_HYPullDown*)GetCellObject (2,1);
        k = pd->FindMenuItem(*(_String*)(*modelT)(0));
        if (k>=0) {
            pd->ChangeSelection (k,true);
        }
    }
}


//__________________________________________________________

/*_HYModelWindow*   _HYModelWindowDialog::SpawnNewModel (void)
{
    _HYPullDown*     pd = (_HYPullDown*)GetCellObject (2,1);
    _String          pds = *pd->GetMenuItem (pd->GetSelection());
    pd = (_HYPullDown*)GetCellObject (3,1);
    return ::SpawnNewModel (mDim, gCode,pds,pd->GetSelection());
}*/

//__________________________________________________________
bool    _HYModelWindow::ConfirmClose (void)
{
    if (taint) {
        _String warnMessage ("Would you like to save '");
        warnMessage = warnMessage & GetTitle() & "' before it is closed?";
        char    res = YesNoCancelPrompt(warnMessage);
        if (res!=2) {
            if (res==1) {
                DoSave(-1);
                return taint;
            }
        } else {
            return false;
        }
    }

    return true;
}

//__________________________________________________________

_String*                DefStringForGCode (_SimpleList* geneticCode)
{
    _String * res = new _String (256L, true),
    buffer;

    checkPointer (res);

    (*res) << &modelGenCode;
    (*res) << '=';
    (*res) << '{';
    (*res) << '{';
    for (long k=0; k<geneticCode->lLength-1; k++) {
        buffer = (*geneticCode)(k);
        (*res) << &buffer;
        (*res) << ',';
    }
    buffer = (*geneticCode)(geneticCode->lLength-1);
    (*res) << &buffer;
    (*res) << '}';
    (*res) << '}';
    (*res) << ';';
    (*res).Finalize();
    return res;
}

//__________________________________________________________

bool                OpenModelFromFile (_String& fileName)
{
    FILE    * F    = doFileOpen (fileName.sData,"rb");
    bool    result = false;
    _String errMsg;
    if (!F) {
        errMsg  = _String ("Can't open file \"") & fileName & "\" for reading.";
    } else {
        _String fileData (F);
        fclose (F);

        long        g = LocateVarByName (rateClassCount);

        if (g<0) {
            _String   four ("4");
            _Variable rcv (rateClassCount,four,false);
        }

        g = batchLanguageFunctionNames.lLength;
        _String * modelDef = new _String (fileData);

        DeleteVariable (modelDataType);
        DeleteVariable (modelFunction);
        DeleteVariable (modelName);
        DeleteVariable (ModelEFVType);
        DeleteVariable (modelEFVVector);
        DeleteVariable (modelDataType);
        DeleteVariable (modelGenCode);

        _ExecutionList    exList (*modelDef);
        exList.Execute();
        DeleteObject (modelDef);
        if (terminateExecution) {
            terminateExecution = false;
        } else {

            long        i1 = LocateVarByName (modelDataType),
                        i2 = batchLanguageFunctionNames.Find (&modelFunction),
                        k;

            if (i2>=0) {
                if (i1>=0) {
                    _Variable * modelSL = FetchVar (i1);
                    if (modelSL->ObjectClass ()==NUMBER) {
                        i1 = modelSL->Compute()->Value();
                        if ((i1==HY_MODEL_TYPE_NUC)||(i1==HY_MODEL_TYPE_DINUC)||(i1==HY_MODEL_TYPE_AA)||(i1==HY_MODEL_TYPE_CODON)) {
                            _String     varID ("DummyEFV"),
                                        EFV,
                                        varID2 ("ModelMatrixTemplate");

                            _SimpleList * fileGenCode = nil;

                            FindUnusedObjectName (EFV,varID,variableNames,true);
                            FindUnusedObjectName (EFV,varID2,variableNames,true);

                            if (i1 == HY_MODEL_TYPE_CODON) {
                                i2 = LocateVarByName (modelGenCode);
                                if (i2<0) {
                                    errMsg = "Codon model definitions must also include their genetic code (by defining \"ModelGeneticCode\").";
                                } else {
                                    modelSL = FetchVar (i2);
                                    if (modelSL->ObjectClass() == MATRIX) {
                                        _Matrix * codeMatrix = (_Matrix*)modelSL->GetValue();
                                        if ((codeMatrix->GetHDim()!=1)||(codeMatrix->GetVDim()!=64)) {
                                            errMsg = "Codon model definitions must also include their genetic code (64 component vector).";
                                        } else {
                                            fileGenCode = new _SimpleList;
                                            checkPointer (fileGenCode);
                                            _Matrix * evalMatrix = (_Matrix*)codeMatrix->ComputeNumeric();
                                            for (i2 = 0; i2 < 64; i2++) {
                                                k = evalMatrix->theData[i2];
                                                if ((k<0)||(k>20)) {
                                                    break;
                                                }
                                                (*fileGenCode) << k;
                                            }
                                            if (i2<64) {
                                                errMsg = "Genetic Code definition vector entries must be integers between 0 and 20.";
                                            }

                                        }
                                    }
                                }
                            }

                            if (errMsg.sLength==0) {
                                EFV = varID & "={";

                                for (long k2=0; k2<64; k2++) {
                                    EFV = EFV & "{1,1,1}\n";
                                }

                                EFV = EFV & "};\n dummy=" & modelFunction & "(\"" & varID2 & "\"," & varID & ");";

                                _ExecutionList xl (EFV);
                                xl.Execute();
                                if (terminateExecution) {
                                    terminateExecution = false;
                                } else {
                                    errMsg = untitled;
                                    k = LocateVarByName (modelName);
                                    if (k>=0) {
                                        modelSL = FetchVar (k);
                                        if (modelSL->ObjectClass ()==STRING) {
                                            errMsg = *((_FString*)modelSL->GetValue())->theString;
                                            if (!errMsg.IsValidIdentifier()) {
                                                errMsg = untitled;
                                            }
                                        }
                                    }

                                    if (errMsg.Equal (&untitled)) {
                                        if (modelSpawnerCount>1) {
                                            errMsg = untitled & modelSpawnerCount;
                                        }
                                        modelSpawnerCount++;
                                    }

                                    _List            stateLabels;

                                    SpawnStandardLabels (i1,stateLabels);

                                    if (fileGenCode)
                                        for (i2=63; i2>=0; i2--)
                                            if (fileGenCode->lData[i2]==10) {
                                                stateLabels.Delete (i2);
                                            }

                                    _HYModelWindow * res = new _HYModelWindow (errMsg,&stateLabels,i1);
                                    res->BuildTemplates (fileGenCode,i1);

                                    if (!res->CopyMatrix (varID2)) {
                                        errMsg = "The file doesn't seem to contain a valid model definition function.";
                                        res->SetTaint (false);
                                        postWindowCloseEvent (res->GetID());
                                    } else {
                                        k = LocateVarByName (ModelEFVType);
                                        if (k>=0) {
                                            _Variable* thisVar = FetchVar(k);
                                            if (thisVar->ObjectClass()==STRING) {
                                                modelDef = ((_FString*)thisVar->Compute())->theString;
                                                res->SetEFVVector (*modelDef);
                                            }
                                        }
                                        k = LocateVarByName (modelEFVVector);
                                        if (k>=0) {
                                            _Variable* thisVar = FetchVar(k);
                                            if (thisVar->ObjectClass()==MATRIX) {
                                                _Matrix * mx = (_Matrix*)thisVar->GetValue();
                                                res->CopyEFVMatrix ((_Matrix*)mx->ComputeNumeric());
                                            }
                                        }
                                        errMsg = empty;
                                        res->BringToFront();
                                    }
                                }
                            }
                            DeleteVariable (varID);
                            DeleteVariable (varID2);
                            DeleteVariable (ModelEFVType);
                            DeleteObject (fileGenCode);
                        }
                    } else {
                        errMsg = modelDataType & " is not one of recognized (0,1,2,3) values.";
                    }

                } else {
                    errMsg = modelDataType & " constant must be defined in the model file.";
                }
            } else {
                errMsg = modelFunction & " function must be defined in the model file.";
            }

            DeleteVariable (modelDataType);
            DeleteVariable (modelGenCode);

            if (errMsg.sLength) {
                ProblemReport (errMsg);
            } else {
                result = true;
            }

        }
        while (g<batchLanguageFunctionNames.countitems()) {
            batchLanguageFunctionNames.Delete (g);
            batchLanguageFunctionParameters.Delete (g);
            batchLanguageFunctions.Delete(g);
            batchLanguageFunctionParameterLists.Delete(g);
            batchLanguageFunctionClassification.Delete(g);
        }
    }
    return result;
}

//__________________________________________________________

void                SpawnStandardLabels (char type, _List& labels)
{
    labels.Clear();
    _String stateLabels,
            nucs (dnaOneCharCodes);

    long   k,j,l;

    switch (type) {
    case HY_MODEL_TYPE_NUC:
        labels.AppendNewInstance(new _String('A'));
        labels.AppendNewInstance(new _String('C'));
        labels.AppendNewInstance(new _String('G'));
        labels.AppendNewInstance(new _String('T'));
        break;
    case HY_MODEL_TYPE_DINUC: {
        for (k=0; k<4; k++)
            for (j=0; j<4; j++) {
                stateLabels = _String(nucs.sData[k])&nucs.sData[j];
                labels && & stateLabels;
            }
        break;
    }
    case HY_MODEL_TYPE_AA: {
        for (k=0; k<20; k++) {
            labels.AppendNewInstance (new _String (aminoAcidOneCharCodes.sData[k]));
        }
        break;
    }
    case HY_MODEL_TYPE_CODON: {
        _String stateLabels;

        for (k=0; k<4; k++)
            for (j=0; j<4; j++)
                for (l=0; l<4; l++) {
                    stateLabels = _String(nucs.sData[k])&nucs.sData[j]&nucs.sData[l];
                    labels && & stateLabels;
                }
        break;
    }
    }
}

//__________________________________________________________
bool                OpenModelFromMatrix (_String& matrixName, bool setMx)
{
    long f = modelNames.Find(&matrixName),
         i1,
         i2,
         k;

    if (f>=0) {
        _Matrix * modelMatrix = (_Matrix*)LocateVar (modelMatrixIndices.lData[f])->GetValue();

        switch (modelMatrix->GetHDim()) {
        case 4:
            i1= HY_MODEL_TYPE_NUC;
            break;
        case 16:
            i1= HY_MODEL_TYPE_DINUC;
            break;
        case 20:
            i1= HY_MODEL_TYPE_AA;
            break;
        default:
            i1= HY_MODEL_TYPE_CODON;

        }

        _List            stateLabels;
        _String          errMsg;

        SpawnStandardLabels (i1,stateLabels);
        _SimpleList * code = nil;

        if (i1==HY_MODEL_TYPE_CODON) {
            _SimpleList matchedTables;
            for (k=0; k<geneticCodes.lLength; k++)
                if (DimensionOfGenCode(k) == modelMatrix->GetHDim()) {
                    matchedTables << k;
                }
            if (matchedTables.lLength) {
                if (matchedTables.lLength>1) {
                    _SimpleList picks,
                                sel;
                    picks << 0;
                    picks << 1;
                    i2 = HandleListSelection (geneticCodes,picks,matchedTables,"Pick a genetic code for model display",sel,1);
                    if (i2<0) {
                        return false;
                    } else {
                        i2 = sel.lData[0];
                    }
                } else {
                    i2 = matchedTables.lData[0];
                }
            } else {
                errMsg = _String ("Can't display codon model \"") & matrixName & "\", b/c its dimensions don't match any genetic code known to HyPhy.";
                ProblemReport (errMsg);
                return false;
            }
            code = ((_SimpleList*)(*(_List*)(geneticCodes(i2)))(2));

            for (i2=63; i2>=0; i2--)
                if (code->lData[i2]==10) {
                    stateLabels.Delete (i2);
                }

        }


        if (modelSpawnerCount>1) {
            errMsg = untitled & modelSpawnerCount;
        } else {
            errMsg = untitled;
        }
        modelSpawnerCount++;

        _HYModelWindow * res = new _HYModelWindow (errMsg,&stateLabels,i1);
        res->BuildTemplates (code,i1);

        if (setMx)
            if (!res->CopySimpleMatrix (modelMatrix)) {
                errMsg = "Failed to copy model matrix. Oops.";
                res->SetTaint (false);
                postWindowCloseEvent (res->GetID());
                ProblemReport (errMsg);
                return false;
            } else {
                i2 = modelFrequenciesIndices.lData[f];
                if (i2<0) {
                    i2 = -i2-1;
                }
                res->CopyEFVMatrix ((_Matrix*)((_Matrix*)LocateVar(i2)->GetValue())->ComputeNumeric());
                res->BringToFront();
            }
        return true;
    }

    return false;
}

//__________________________________________________________
bool                OpenModelFromCalcNode (_String& nodeName, bool expMx = false)
{
    long f = LocateVarByName (nodeName),
         i1,
         i2;

    if (f>=0) {
        _Variable* thisNode = FetchVar(f);
        if (thisNode&&(thisNode->ObjectClass()==TREE_NODE)) {
            _CalcNode* thisCNode = (_CalcNode*)thisNode;

            i1 = thisCNode->GetModelIndex();
            if (i1 != HY_NO_MODEL) {
                if (OpenModelFromMatrix(*(_String*)modelNames(i1),false)) {
                    _String errMsg;

                    if (modelSpawnerCount>2) {
                        errMsg = untitled & (modelSpawnerCount-1);
                    } else {
                        errMsg = untitled;
                    }

                    errMsg = titlePrefix & errMsg;

                    _HYModelWindow* res = (_HYModelWindow*)windowObjectRefs(FindWindowByName (errMsg));

                    _Matrix       * modelMatrix = (_Matrix*)thisCNode->ComputeModelMatrix (false)->MultByFreqs (i1);
                    if (expMx) {
                        modelMatrix = (_Matrix*) modelMatrix->Exponentiate();
                    }

                    if (res->CopySimpleMatrix (modelMatrix)) {
                        i2 = modelFrequenciesIndices.lData[i1];
                        if (i2<0) {
                            i2 = -i2-1;
                        }
                        res->CopyEFVMatrix ((_Matrix*)((_Matrix*)LocateVar(i2)->GetValue())->ComputeNumeric());
                        res->BringToFront();
                        res->SetTaint (false);
                        _String tryMe;
                        if (expMx) {
                            tryMe = "Transition_Matrix_For_";
                        } else {
                            tryMe = "Rate_Matrix_For_";
                        }

                        tryMe = tryMe & nodeName;
                        errMsg = tryMe;
                        i2 = 1;
                        while (FindWindowByName(errMsg)>=0) {
                            i2 ++;
                            errMsg = tryMe & i2;
                        }
                        res->SetTitle (errMsg);
                        if (expMx) {
                            DeleteObject (modelMatrix);
                        }
                        return true;
                    }
                    if (expMx) {
                        DeleteObject (modelMatrix);
                    }
                }
            }
        }
    }

    return false;
}

//__________________________________________________________________
_HYModelParameterDialog::_HYModelParameterDialog     (_HYModelWindow* tp, long in, long* r):_HYTWindow ("Parameter Properties",false,true, (Ptr)tp)
{

    result      = r;
    inID        = in;
    parentModel = tp;

    _HYRect         canvasSettings = {30,100,30,100,HY_COMPONENT_TRANSP_BG};

    _HYLabel*       l1 = new _HYLabel (canvasSettings, GetOSWindowData());
    checkPointer    (l1);

    _HYLabel*       l2 = new _HYLabel (canvasSettings, GetOSWindowData());
    checkPointer    (l2);

    canvasSettings.top = canvasSettings.bottom = 60;

    _HYLabel*       l3 = new _HYLabel (canvasSettings, GetOSWindowData());
    checkPointer    (l3);

    canvasSettings.left = canvasSettings.right = 200;
    canvasSettings.top = canvasSettings.bottom = 30;

    _HYTextBox*     tb = new _HYTextBox (canvasSettings, GetOSWindowData());
    checkPointer    (tb);

    _HYCheckbox*    rb1 = new _HYCheckbox (canvasSettings, GetOSWindowData(),true);
    checkPointer    (rb1);
    _HYCheckbox*    rb2 = new _HYCheckbox (canvasSettings, GetOSWindowData(),true);
    checkPointer    (rb2);
    _HYCheckbox*    rb3 = new _HYCheckbox (canvasSettings, GetOSWindowData(),true);
    checkPointer    (rb3);


    canvasSettings.bottom = canvasSettings.top   = 36;
    canvasSettings.left   = canvasSettings.right = 210;

    _HYButton*      b1      = new _HYButton (canvasSettings, GetOSWindowData());
    canvasSettings.left   = canvasSettings.right = 90;
    _HYButton*      b2      = new _HYButton (canvasSettings, GetOSWindowData());

    checkPointer (b1);
    checkPointer (b2);

    tb->SetMessageRecipient (this);
    rb1->SetMessageRecipient (this);
    rb2->SetMessageRecipient (this);
    rb3->SetMessageRecipient (this);
    b1->SetMessageRecipient (this);
    b2->SetMessageRecipient (this);

    AddObject (l1);// 0
    AddObject (l2);// 1
    AddObject (rb1);// 2
    AddObject (rb2);// 3
    AddObject (rb3);// 4
    AddObject (tb);// 5
    AddObject (b1);// 6
    AddObject (b2);// 7
    AddObject (l3);// 7

    SetTableDimensions (5,2);

    SetCell (0,0,l1);
    SetCell (0,1,tb);
    SetCell (1,0,l2);
    SetCell (1,1,rb1);
    SetCell (2,0,l3);
    SetCell (2,1,rb2);
    SetCell (3,0,l3);
    SetCell (3,1,rb3);
    SetCell (4,0,b1);
    SetCell (4,1,b2);

    _HYFont         labelFont;
#ifdef __MAC__
    labelFont.face = "System Font";
    labelFont.size = 12;
#endif
#ifdef __WINDOZE__
    labelFont.face = "MS Sans Serif";
    labelFont.size = 12;
#endif
#ifdef __HYPHY_GTK__
    labelFont.face = _HY_SANS_FONT;
    labelFont.size = 10;
#endif

    labelFont.style = HY_FONT_PLAIN;

    b1->SetFont (labelFont);
    b2->SetFont (labelFont);
    l1->SetFont (labelFont);
    l2->SetFont (labelFont);
    rb1->SetFont (labelFont);
    rb2->SetFont (labelFont);
    rb3->SetFont (labelFont);
    tb->SetText (" ");
    tb->SetFont (labelFont);
    tb->SetText (empty);

    b1->SetText (" OK ");
    b2->SetText (" Cancel ");

    b1->SetAlignFlags (HY_ALIGN_RIGHT);
    l1->SetAlignFlags (HY_ALIGN_LEFT);
    l2->SetAlignFlags (HY_ALIGN_LEFT);
    rb1->SetAlignFlags (HY_ALIGN_LEFT);
    rb2->SetAlignFlags (HY_ALIGN_LEFT);
    rb3->SetAlignFlags (HY_ALIGN_LEFT);

    l1->SetText ("Parameter ID:");
    l2->SetText ("Parameter Class:");
    rb1->SetText ("Local");
    rb2->SetText ("Global");
    rb3->SetText ("Category");

    keyboardFocusChain << 5;
    ProcessEvent (generateKeyboardFocusEvent (tb->GetID()));

    orb = 2;

    if (inID>=0) {
        _String newSID = parentModel->RetrieveVarID (inID);
        tb->SetText  (newSID);
        tb->SetForeColor (editBoxGood);
        orb = 2+parentModel->VarKind (inID);
        b1->EnableButton (true);
    } else {
        b1->EnableButton (false);
    }


    ((_HYCheckbox*)components (orb))->SetState (true,false);
    rb = orb;

    b1->SetButtonKind (HY_BUTTON_OK);
    b2->SetButtonKind (HY_BUTTON_CANCEL);

    _HYRect             dim = MinMaxWindowDimensions();
    SetWindowRectangle  (0,0,dim.bottom,dim.right);
    CenterWindow        (this);

    DeleteObject (b1);
    DeleteObject (b2);
    DeleteObject (l1);
    DeleteObject (l2);
    DeleteObject (l3);
    DeleteObject (rb1);
    DeleteObject (rb2);
    DeleteObject (rb3);
    DeleteObject (tb);
}

//__________________________________________________________

bool    _HYModelParameterDialog::ProcessEvent (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    long        i,f;

    if (e->EventClass()==_hyTextEditChange) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);

        if (i==5) {
            _HYTextBox * tb = (_HYTextBox*)components (i);
            bool    isGood = false;

            _String newSID = tb->GetText();

            if (newSID.IsValidIdentifier()) {
                f = parentModel->FindVarID (newSID);
                if (f>=0) {
                    if (inID == f) {
                        isGood = true;
                    }
                } else {
                    isGood = true;
                }
            }

            tb->SetForeColor (isGood?editBoxGood:editBoxBad);

            ((_HYButton*)components(6))->EnableButton (isGood);

        }
        done = true;
    } else if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);

        if (i==6) { // OK
            _HYTextBox * tb = (_HYTextBox*)components (5);
            _String newSID = tb->GetText();
            *result = parentModel->lastParameterChoice;
            if (newSID.sLength) {
                if (inID!=-1) {
                    if ((newSID.Equal(&parentModel->RetrieveVarID (inID)))&&(orb==rb)) {
                        *result = inID;
                    } else {
                        parentModel->RenameVarID (newSID,inID,rb-2);
                    }
                } else {
                    parentModel->InsertVarID (newSID, rb-2);
                    *result = parentModel->FindVarID (newSID);
                }
            }
            //*result = 0;
            postWindowCloseEvent (GetID());
        } else if (i==7) {
            *result = parentModel->lastParameterChoice;
            postWindowCloseEvent (GetID());
        } else if ((i>=2)&&(i<5)) {
            if (i!=rb) {
                ((_HYCheckbox*)components (rb))->SetState (false,false);
                rb = i;
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

//__________________________________________________________

_HYModelWindow* SpawnNewModel (long mDim, long gCode, _String& s1, long s2)
{
    _List   labels;

    _String stateLabels,
            nucs (dnaOneCharCodes);

    _SimpleList
    *genC = nil;

    char    type;
    long    k,j,l;

    switch (mDim) {
    case 4:
        labels.AppendNewInstance (new _String('A'));
        labels.AppendNewInstance (new _String('C'));
        labels.AppendNewInstance (new _String('G'));
        labels.AppendNewInstance (new _String('T'));
        type = HY_MODEL_TYPE_NUC;
        break;
    case 16: {
        for (k=0; k<4; k++)
            for (j=0; j<4; j++) {
                stateLabels = _String(nucs.sData[k])&nucs.sData[j];
                labels && & stateLabels;
            }
        type = HY_MODEL_TYPE_DINUC;
        break;
    }
    case 20: {
        for (k=0; k<20; k++) {
            labels.AppendNewInstance (new _String (aminoAcidOneCharCodes.sData[k]));
        }
        type = HY_MODEL_TYPE_AA;
        break;
    }
    default: {
        genC = (_SimpleList*)(*(_List*)geneticCodes(gCode))(2);
        for (k=0; k<4; k++)
            for (j=0; j<4; j++)
                for (l=0; l<4; l++) {
                    if (genC->lData[k*16+j*4+l]!=10) {
                        stateLabels = _String(nucs.sData[k])&nucs.sData[j]&nucs.sData[l];
                        labels && & stateLabels;
                    }
                }
        type = HY_MODEL_TYPE_CODON;

    }
    }

    //_HYPullDown*       pd = (_HYPullDown*)GetCellObject (2,1);


    _HYModelWindow * res = new _HYModelWindow (untitled&(modelSpawnerCount==1?empty:modelSpawnerCount), &labels, type);
    modelSpawnerCount++;
    res->BuildTemplates (genC, type);
    if (!s1.Equal (&none)) {
        _List* theModel = FindModelTemplate(&s1);
        if (theModel) {
            FILE * thisModel = doFileOpen (((_String*)(*theModel)(2))->getStr(),"rb");
            if (thisModel) {
                long   g = batchLanguageFunctionNames.lLength;
                _String * modelDef = new _String (thisModel);
                fclose (thisModel);
                //pd = (_HYPullDown*)GetCellObject (3,1);
                *modelDef = _String ("modelType=")&_String (s2)&";\nrateClassCount=4;\n"&*modelDef;
                if (type==HY_MODEL_TYPE_CODON) {
                    _String * gcs = DefStringForGCode (genC);
                    *modelDef = *gcs&*modelDef;
                    DeleteObject (gcs);
                }
                _String varID ("DummyEFV"),
                        EFV,
                        varID2 ("ModelMatrixTemplate");

                FindUnusedObjectName (EFV,varID,variableNames,true);
                FindUnusedObjectName (EFV,varID2,variableNames,true);

                EFV = varID & "={";

                for (long k2=0; k2<64; k2++) {
                    EFV = EFV & "{1,1,1}\n";
                }

                EFV = EFV & "};\n dummy=" & modelFunction & "(\"" &
                      varID2 & "\"," & varID & ");";

                *modelDef = *modelDef & '\n' & EFV;

                DeleteVariable (ModelEFVType);
                DeleteVariable (modelEFVVector);

                _ExecutionList ex (*modelDef);
                ex.Execute();
                DeleteObject   (modelDef);

                res->CopyMatrix (varID2);

                k = LocateVarByName (ModelEFVType);
                if (k>=0) {
                    _Variable* thisVar = FetchVar(k);
                    if (thisVar->ObjectClass()==STRING) {
                        modelDef = ((_FString*)thisVar->Compute())->theString;
                        res->SetEFVVector (*modelDef);
                    }
                }

                k = LocateVarByName (modelEFVVector);
                if (k>=0) {
                    _Variable* thisVar = FetchVar(k);
                    if (thisVar->ObjectClass()==MATRIX) {
                        _Matrix * mx = (_Matrix*)thisVar->GetValue();
                        res->CopyEFVMatrix ((_Matrix*)mx->ComputeNumeric());
                    }
                }


                DeleteVariable (varID);
                DeleteVariable (varID2);
                DeleteVariable (ModelEFVType);

                while (g<batchLanguageFunctionNames.lLength) {
                    batchLanguageFunctionNames.Delete (g);
                    batchLanguageFunctionParameters.Delete (g);
                    batchLanguageFunctions.Delete(g);
                    batchLanguageFunctionParameterLists.Delete(g);
                    batchLanguageFunctionClassification.Delete(g);
                }
            }
        }
    }
    return res;
}

// EOF