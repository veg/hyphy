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
#include "HYDialogs.h"
#include "HYTableWindow.h"
#include "HYButton.h"
#include "HYLabel.h"
#include "HYCanvas.h"
#include "HYUtils.h"
#include "HYEventTypes.h"
#include "HYTextBox.h"
#include "HYButtonBar.h"
#include "HYPullDown.h"
#include "HYTableComponent.h"
#include "HYTreePanel.h"
#include "HYChartWindow.h"
#include "HYConsoleWindow.h"
#include "HYDataPanel.h"
#include "HYModelWindow.h"
#include "HYDBWindow.h"
#include "HYSharedMain.h"

#include "batchlan.h"
#include "string.h"
#include "likefunc.h"
#include "ctype.h"
#include "math.h"

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

_HYRect  hierListTextBoxBounds  = {25,20,145,280,0},
         hierListTextBoxBounds2 = {15,10,145,165,0};

_List    cachedDialogTitles,
         cachedDialogSelections;

void     SetDefaultDialogFont   (_HYFont&);


extern  _String         none;

long            combSizeMemory = 3;
_SimpleList     combSelectionMemory;
_List           wiseCrackButtons;

//__________________________________________________________

void            SetDefaultDialogFont (_HYFont& labelFont)
{
#ifdef __MAC__
    labelFont.face  = "System Font";
    labelFont.size  = 12;
#else
#ifdef __WINDOZE__
    labelFont.face  = "Arial";
    labelFont.size  = 14;
#else
    labelFont.face  = _HY_SANS_FONT;
    labelFont.size  = 10;
#endif
#endif
    labelFont.style = HY_FONT_PLAIN;
}

//__________________________________________________________
// Begin font dialog
//__________________________________________________________

_SimpleList     customFontSizes;

// stores custom font sizes between calls

//__________________________________________________________

_HYFontDialog::_HYFontDialog (_HYFont& startWith,_HYWindow* rec):_HYTWindow ("HyPhy Font Selector", false, true, (Ptr)rec)
{

    _HYRect         canvasSettings = {30,100,30,100,HY_COMPONENT_NO_SCROLL|HY_COMPONENT_TRANSP_BG};

    _HYLabel*       l1      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l2      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l3      = new _HYLabel (canvasSettings, GetOSWindowData());


    canvasSettings.left = canvasSettings.right = 200;

    _HYPullDown*    p1      = new _HYPullDown (canvasSettings, GetOSWindowData());
    _HYPullDown*    p2      = new _HYPullDown (canvasSettings, GetOSWindowData());
    _HYPullDown*    p3      = new _HYPullDown (canvasSettings, GetOSWindowData());

    canvasSettings.left = canvasSettings.right = 220;
    _HYButton*      b1      = new _HYButton (canvasSettings, GetOSWindowData());

    canvasSettings.left = canvasSettings.right = 81;
    _HYButton*      b2      = new _HYButton (canvasSettings, GetOSWindowData());

    canvasSettings.left = canvasSettings.right = 300;
    canvasSettings.top = canvasSettings.bottom = 50;
    canvasSettings.width = HY_COMPONENT_WELL|HY_COMPONENT_TRANSP_BG;

    _HYLabel*       l4      = new _HYLabel (canvasSettings, GetOSWindowData());

    p1->SetMessageRecipient (this);
    p2->SetMessageRecipient (this);
    p3->SetMessageRecipient (this);
    b1->SetMessageRecipient (this);
    b2->SetMessageRecipient (this);

    AddObject (l1);
    AddObject (l2);
    AddObject (l3);
    AddObject (l4);

    AddObject (p1);
    AddObject (p2);
    AddObject (p3);

    AddObject (b1);
    AddObject (b2);


    SetTableDimensions (5,3);
    SetCell   (0,0,l1);
    SetCell   (0,1,p1);
    SetCell   (0,2,p1);

    SetCell   (1,0,l2);
    SetCell   (1,1,p2);
    SetCell   (1,2,p2);

    SetCell   (2,0,l3);
    SetCell   (2,1,p3);
    SetCell   (2,2,p3);

    SetCell   (3,0,l4);
    SetCell   (3,1,l4);
    SetCell   (3,2,l4);

    SetCell   (4,0,b1);
    SetCell   (4,1,b1);
    SetCell   (4,2,b2);


    _HYFont  labelFont;
    SetDefaultDialogFont (labelFont);

    l1->SetFont (labelFont);
    l2->SetFont (labelFont);
    l3->SetFont (labelFont);

    b1->SetFont (labelFont);
    b2->SetFont (labelFont);

    myFont = startWith;
    l4->SetFont (myFont);

    l1->SetAlignFlags (HY_ALIGN_LEFT);
    l2->SetAlignFlags (HY_ALIGN_LEFT);
    l3->SetAlignFlags (HY_ALIGN_LEFT);

    b1->SetAlignFlags (HY_ALIGN_RIGHT);

    p1->SetAlignFlags (HY_ALIGN_LEFT);
    p2->SetAlignFlags (HY_ALIGN_LEFT);
    p3->SetAlignFlags (HY_ALIGN_LEFT);

    l1->SetText       ("Font:");
    l2->SetText       ("Size:");
    l3->SetText       ("Style:");

    l4->SetText       ("The answer is 42.");

    b1->SetText       ("  OK  ");
    b1->SetButtonKind (HY_BUTTON_OK);
    b2->SetText       (" Cancel ");
    b2->SetButtonKind (HY_BUTTON_CANCEL);

    _List fonts;

    GenerateFontList (fonts);

    long k;

    for (k=0; k<fonts.lLength; k++) {
        p1->AddMenuItem (*(_String*)fonts(k),-1);
    }

    if (customFontSizes.lLength==0) {
        customFontSizes << 4;
        customFontSizes << 7;
        customFontSizes << 8;
        customFontSizes << 9;
        customFontSizes << 10;
        customFontSizes << 11;
        customFontSizes << 12;
        customFontSizes << 14;
        customFontSizes << 18;
    }

    k = customFontSizes.Find (myFont.size);
    if (k<0) {
        customFontSizes << myFont.size;
    }

    customFontSizes.Sort();

    for (k=0; k< customFontSizes.lLength; k++) {
        p2->AddMenuItem (_String(customFontSizes.lData[k]),-1);
    }

    p2->AddMenuItem   (menuSeparator,-1);
    p2->AddMenuItem   ("Other...",-1);

    p3->AddMenuItem   ("Plain",-1);
    p3->AddMenuItem   ("Bold",-1);
    p3->AddMenuItem   ("Italic",-1);
    p3->AddMenuItem   ("Bold Italic",-1);


    k = p1->FindMenuItem (myFont.face);
    if (k>=0) {
        p1->ChangeSelection (p1->FindMenuItem (myFont.face),false);
    } else {
        p1->AddMenuItem (myFont.face,0);
    }


    k = customFontSizes.Find (myFont.size);
    p2->ChangeSelection (k,false);

    p3->ChangeSelection (myFont.style, false);

    DeleteObject (l1);
    DeleteObject (l2);
    DeleteObject (l3);
    DeleteObject (l4);
    DeleteObject (p1);
    DeleteObject (p2);
    DeleteObject (p3);
    DeleteObject (b1);
    DeleteObject (b2);

    _HYRect     dim = MinMaxWindowDimensions();

    SetWindowRectangle  (0,0,dim.bottom,dim.right);
    CenterWindow        (this);
    firstFont = myFont;
    mr = rec;
}

//__________________________________________________________

bool    _HYFontDialog::ProcessEvent (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    long        i,f,g,k;

    if (e->EventClass()==_hyMenuSelChangeEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        g = e->EventCode().Cut (f+1,-1).toNum();

        for (i=0; i<components.lLength; i++)
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }

        if (i==4) {
            myFont.face = *((_HYPullDown*)components(i))->GetMenuItem(g);
        } else if (i==5) {
            _HYPullDown * p2 = (_HYPullDown*)components(i);
            if (g!=p2->MenuItemCount()-1) {
                myFont.size = p2->GetMenuItem(g)->toNum();
            } else {
                _String prompt ("Desired Font Size (2-255):"),
                        value  ((long)myFont.size);
                if (EnterStringDialog (value,prompt)) {
                    g = value.toNum();
                    if ((g>1)&&(g<255)) {
                        k=customFontSizes.Find(g);
                        if (k==-1) {
                            customFontSizes<<g;
                            p2->AddMenuItem (_String(g),p2->MenuItemCount()-2);
                        }
                        myFont.size = g;
                    }
                }
                p2->ChangeSelection (customFontSizes.Find (myFont.size),false);
            }
        } else if (i==6) {
            myFont.style = g;
        }

        ((_HYLabel*)components (3))->SetFont (myFont);
        done = true;
    } else if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++)
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }

        if (i==7) { // OK
            if ((!myFont.face.Equal(&firstFont.face))||(myFont.size!=firstFont.size)||(myFont.style!=firstFont.style)) {
                mr->SetFont (myFont);
            }
            postWindowCloseEvent (GetID());
        } else if (i==8) { // Cancel
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

//__________________________________________________________________
_HYListSelectDialog::_HYListSelectDialog     (_List* d, _SimpleList* c, _SimpleList* vc, _String n, _SimpleList* s, long* r)
    :_HYTWindow (n,false,true)
{
    data            = d;
    choices         = c;
    validChoices    = vc;
    selections      = s;

    result          = r;

    _List           hData,
                    rData,
                    iData,
                    diData;


    long            counter;

    keyboardFocusChain << 0;

    for (counter = 0; counter< (*validChoices).lLength; counter++) {
        _String* tStr = (_String*)(*(_List*)((*data)(validChoices->lData[counter])))(choices->lData[0]);
        if (tStr->sData[0]=='!') {
            if (rData.lLength) {
                rData && & iData;
                hData && & rData;
                dData && & diData;
            }
            rData.Clear();
            iData.Clear();
            diData.Clear();
            diData << (*(_List*)(*data)((*validChoices)(counter)))((*choices)(1));
            _String tStr2 (*tStr,1,-1);
            rData && & tStr2;
        } else {
            iData  << tStr;
            diData << (*(_List*)(*data)((*validChoices)(counter)))((*choices)(1));
        }
    }

    rData && & iData;
    hData && & rData;
    dData && & diData;

    _HYRect         canvasSettings = {200,300,200,300,HY_COMPONENT_V_SCROLL|HY_COMPONENT_TRANSP_BG};

    _HYHList*       hl      = new _HYHList (canvasSettings, GetOSWindowData(),hData);

    hl->settings.bottom = MAX(200,hl->GetMaxH());
    hl->selectionType |= HY_TABLE_SINGLE_SELECTION|HY_TABLE_FOCUSABLE;

    canvasSettings.top   = canvasSettings.bottom = 150;

    canvasSettings.width = HY_COMPONENT_NO_SCROLL|HY_COMPONENT_TRANSP_BG;
    _HYCanvas*      cn      = new _HYCanvas (canvasSettings,GetOSWindowData(),150,300,32);

    canvasSettings.bottom = canvasSettings.top   = 40;
    canvasSettings.left   = canvasSettings.right = 220;

    _HYButton*      b1      = new _HYButton (canvasSettings, GetOSWindowData());
    canvasSettings.left   = canvasSettings.right = 80;
    _HYButton*      b2      = new _HYButton (canvasSettings, GetOSWindowData());

    canvasSettings.left  = canvasSettings.right  = 300;
    canvasSettings.top   = canvasSettings.bottom = 10;
    canvasSettings.width = HY_COMPONENT_BORDER_B|HY_COMPONENT_TRANSP_BG;

    _HYLabel*       l1      = new _HYLabel  (canvasSettings, GetOSWindowData());


    hl->SetMessageRecipient (this);
    b1->SetMessageRecipient (this);
    b2->SetMessageRecipient (this);

    AddObject (hl);
    AddObject (cn);
    AddObject (b1);
    AddObject (b2);
    AddObject (l1);

    SetTableDimensions (3,2);
    SetCell   (0,0,l1);
    SetCell   (0,1,l1);

    SetCell   (0,0,hl);
    SetCell   (0,1,hl);

    SetCell   (1,0,cn);
    SetCell   (1,1,cn);

    SetCell   (2,0,b1);
    SetCell   (2,1,b2);

    //b1->SetBackColor (bgRGB);
    //b2->SetBackColor (bgRGB);
    //l1->SetBackColor (bgRGB);


    _HYFont         labelFont;
    SetDefaultDialogFont (labelFont);

    b1->SetFont (labelFont);
    b2->SetFont (labelFont);

    cn->StartDraw();
    //cn->SetBColor (bgRGB);
    cn->SetDialogBG();
    cn->EraseAll  ();
    cn->SetColor  ((_HYColor) {
        0,0,0
    });

#ifndef __MAC__
    labelFont.face  = "Verdana";
    labelFont.size  = 12;
#else
    labelFont.face = "Geneva";
    labelFont.size = 10;
#endif

    labelFont.style = HY_FONT_BOLD;
    canvasSettings.width  = 1;
    canvasSettings.top = canvasSettings.bottom = canvasSettings.left = 0;
    canvasSettings.right  = 300;
    cn->DrawLine (canvasSettings);
    cn->SetFont (labelFont);
    canvasSettings.left = 5;
    canvasSettings.top  = 15;
    canvasSettings.bottom = 148;
    canvasSettings.right  = 293;
    cn->DrawInfoBox (canvasSettings,"Item Description");
    labelFont.style = HY_FONT_PLAIN;
    cn->SetFont (labelFont);
    cn->EndDraw();

    b1->SetAlignFlags (HY_ALIGN_RIGHT);

    b1->SetText       ("  OK  ");
    b1->SetButtonKind (HY_BUTTON_OK);
    b2->SetText       (" Cancel ");
    b2->SetButtonKind (HY_BUTTON_CANCEL);
    b1->EnableButton  (false);

    DeleteObject (b1);
    DeleteObject (b2);
    DeleteObject (l1);
    DeleteObject (hl);
    DeleteObject (cn);

    _HYRect     dim         = MinMaxWindowDimensions(),
                screenDim     = GetScreenDimensions   ();

    SetWindowRectangle  (0,0,MIN(screenDim.bottom-200,dim.bottom),dim.right);
    CenterWindow        (this);
    *result = -1;

}


//__________________________________________________________________
bool    _HYListSelectDialog::ProcessEvent    (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    long        i,f,g,k;

    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++)
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }

        if (i==2) { // OK
            _HYHList * theList = (_HYHList*)components  (0);

            _SimpleList s;
            theList->GetSelection (s);
            if (s.lLength == 1) {
                k = theList->RubrikIndex (s.lData[0]/2);
                f = k>>16;
                g = k&0xffff;
                k = 0;
                for (f=f-1; f>=0; f--) {
                    k+=((_List*)dData(f))->lLength;
                }
                *result = k+g;
                StoreDialogSelection (&GetTitle(),(_String*)theList->GetCellData (s.lData[0]%2,s.lData[0]/2));
            }
            postWindowCloseEvent (GetID());
        } else if (i==3) { // Cancel
            postWindowCloseEvent (GetID());
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
            _HYHList * theList = (_HYHList*)components  (0);
            _HYCanvas* cn      = (_HYCanvas*)components (1);
            _HYButton* b1      = (_HYButton*)components (2);

            _SimpleList s;
            theList->GetSelection (s);
            cn->StartDraw();
            if (s.lLength!=1) {
                cn->EraseRect(hierListTextBoxBounds);
                b1->EnableButton (false);
            } else {
                k = theList->RubrikIndex (s.lData[0]/2);
                f = k>>16;
                g = k&0xffff;
                b1->EnableButton (g);
                _String * thisString = (_String*)(*(_List*)(dData(f)))(g);
                //cn->CrossfadeText (lastString,*thisString,hierListTextBoxBounds,6,30,HY_ALIGN_LEFT,2);
                cn->DisplayText (*thisString,hierListTextBoxBounds,HY_ALIGN_LEFT);
                lastString = *thisString;
            }
            cn->EndDraw();
            cn->_MarkForUpdate();
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
            _HYHList * theList = (_HYHList*)components  (0);

            _SimpleList s;
            theList->GetSelection (s);
            k = theList->RubrikIndex (s.lData[0]/2);
            f = k>>16;
            g = k&0xffff;
            if (g) {
                k = 0;
                for (f=f-1; f>=0; f--) {
                    k+=((_List*)dData(f))->lLength;
                }
                *result = k+g;
                postWindowCloseEvent (GetID());
            } else {
                theList->HandleKeyMove (theList->IsRubrikOpen(f)?2:3,false);
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

//____________________________________________________________________________________________

void        _HYListSelectDialog::SetInitialSelection (void)
{
    long       index = 0;
    _HYHList * theList = (_HYHList*)components  (0);
    _String* sel   = RetrieveCachedSelection (&GetTitle());

    if (sel) {
        index = theList->FindString (sel);
        if (index<0) {
            index = 0;
        } else {
            long rubrikIndex = (index >> 16);
            theList->ModifySelection (0,rubrikIndex,false,false,false);
            index = rubrikIndex + (index&0xffff) + 1;
        }
    }

    _SimpleList  ns (2*index+1);
    theList->ClearSelection      (false);
    theList->SetSelection        (ns,true);
    theList->_MarkCellsForUpdate (ns);
    theList->ScrollToRow         (index);
    ProcessEvent (generateKeyboardFocusEvent (theList->GetID()));
}

//__________________________________________________________________
_HYSimpleListSelectDialog::_HYSimpleListSelectDialog     (_List* d, _SimpleList* c, _SimpleList* vc, _String n, _SimpleList* s, long ns, long* r, Ptr ptr)
    :_HYTWindow (n,false,true,ptr)
{
    choices         = c;
    validChoices    = vc;
    selections      = s;

    result          = r;
    reqSel          = ns;

    long            counter;

    _HYRect         canvasSettings = {150,300,150,300,HY_COMPONENT_V_SCROLL|HY_COMPONENT_TRANSP_BG};

    _HYTable *      hl      = new _HYTable (canvasSettings, GetOSWindowData(),validChoices->lLength,1,301-HY_SCROLLER_WIDTH,20,
                                            HY_TABLE_STATIC_TEXT);

    for (counter = 0; counter< (*validChoices).lLength; counter++) {
        hl->SetCellData ((*(_List*)((*d)(validChoices->lData[counter])))(choices->lData[0]),counter,0,HY_TABLE_STATIC_TEXT,true);
        _String* tStr = (_String*)(*(_List*)((*d)(validChoices->lData[counter])))(choices->lData[1]);
        dData << tStr;
    }

    counter = hl->GetMaxH();
    if (counter<150) {
        hl->AddRow (-1,150-counter,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT);
    }

    hl->selectionType |= HY_TABLE_FOCUSABLE|HY_TABLE_NODRAG_SELECTION;

    if (reqSel==1) {
        hl->selectionType |= HY_TABLE_SINGLE_SELECTION;
    }

    /*canvasSettings.top   = canvasSettings.bottom = 150;
    canvasSettings.left  = canvasSettings.right  = 2;*/

    keyboardFocusChain << 0;

    canvasSettings.width = HY_COMPONENT_NO_SCROLL|HY_COMPONENT_TRANSP_BG;
    _HYCanvas*      cn      = new _HYCanvas (canvasSettings,GetOSWindowData(),150,300,32);

    canvasSettings.bottom = canvasSettings.top   = 40;
    canvasSettings.left   = canvasSettings.right = 80;

    _HYButton*      b1      = new _HYButton (canvasSettings, GetOSWindowData());
    canvasSettings.left   = canvasSettings.right = 80;
    _HYButton*      b2      = new _HYButton (canvasSettings, GetOSWindowData());

    canvasSettings.left  = canvasSettings.right  = 140;

    _HYLabel*       l1      = new _HYLabel  (canvasSettings, GetOSWindowData());


    hl->SetMessageRecipient (this);
    b1->SetMessageRecipient (this);
    b2->SetMessageRecipient (this);

    AddObject (hl);
    AddObject (cn);
    AddObject (b1);
    AddObject (b2);
    AddObject (l1);

    SetTableDimensions (3,3);

    SetCell   (0,0,hl);
    SetCell   (0,1,hl);
    SetCell   (0,2,hl);

    SetCell   (1,0,cn);
    SetCell   (1,1,cn);
    SetCell   (1,2,cn);

    SetCell   (2,0,l1);
    SetCell   (2,1,b1);
    SetCell   (2,2,b2);

    //b1->SetBackColor (bgRGB);
    //b2->SetBackColor (bgRGB);
    //l1->SetBackColor (bgRGB);
    l1->SetForeColor ((_HYColor) {
        0,0,0
    });


    _HYFont         labelFont;
    SetDefaultDialogFont (labelFont);

    b1->SetFont (labelFont);
    b2->SetFont (labelFont);

    cn->StartDraw();
    cn->SetDialogBG();
    cn->EraseAll  ();
    cn->SetColor  ((_HYColor) {
        0,0,0
    });

#ifdef __WINDOZE__
    labelFont.face  = "Verdana";
    labelFont.size  = 12;
#else
#ifdef __MAC__
    labelFont.face = "Geneva";
    labelFont.size = 10;
#else
    labelFont.face = _HY_SANS_FONT;
    labelFont.size = 12;
#endif
#endif

    labelFont.style = HY_FONT_BOLD;
    canvasSettings.width  = 1;
    canvasSettings.top = canvasSettings.bottom = canvasSettings.left = 0;
    canvasSettings.right  = 300;
    cn->DrawLine (canvasSettings);
    cn->SetFont (labelFont);
    canvasSettings.left = 5;
    canvasSettings.top  = 15;
    canvasSettings.bottom = 148;
    canvasSettings.right  = 293;
    cn->DrawInfoBox (canvasSettings,"Item Description");
    labelFont.style = HY_FONT_PLAIN;
    cn->SetFont (labelFont);
    cn->EndDraw();
    l1->SetFont       (labelFont);
    b1->SetAlignFlags (HY_ALIGN_RIGHT);
    _String labelText;
    if (reqSel) {
        labelText = reqSel;
    } else {
        labelText = ">0";
    }
    labelText = labelText &" required (0 chosen).";

    l1->SetText (labelText);

    b1->SetText       ("  OK  ");
    b1->SetButtonKind (HY_BUTTON_OK);
    b2->SetText       (" Cancel ");
    b2->SetButtonKind (HY_BUTTON_CANCEL);
    b1->EnableButton  (false);

    DeleteObject (b1);
    DeleteObject (b2);
    DeleteObject (l1);
    DeleteObject (hl);
    DeleteObject (cn);

    _HYRect     dim  = MinMaxWindowDimensions(),
                sDim = GetScreenDimensions();


    if (dim.bottom>sDim.bottom-sDim.top-90) {
        dim.bottom = sDim.bottom-sDim.top-90;
    }

    if (dim.right>sDim.right-sDim.left-20) {
        dim.right=sDim.right-sDim.left-20;
    }

    SetWindowRectangle  (0,0,dim.bottom,dim.right);
    CenterWindow        (this);
    *result = -1;

}


//__________________________________________________________________
bool    _HYSimpleListSelectDialog::ProcessEvent  (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    long        i,f,k;

    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++)
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }

        if (i==2) { // OK
            _HYHList * theList = (_HYHList*)components  (0);

            selections->Clear();
            theList->GetSelection (*selections);
            //for (k=0; k<selections->lLength; k++)
            //  selections->lData[k] = validChoices->lData[selections->lData[k]];

            if (reqSel == 1) {
                *result = selections->lData[0];
                StoreDialogSelection (&GetTitle(),(_String*)theList->GetCellData (0,*result));
            } else {
                *result = selections->lLength;
            }


            postWindowCloseEvent (GetID());
        } else if (i==3) { // Cancel
            postWindowCloseEvent (GetID());
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
            _HYTable * theList = (_HYHList*)components  (0);
            _HYCanvas* cn      = (_HYCanvas*)components (1);
            _HYButton* b1      = (_HYButton*)components (2);
            _HYLabel * l1      = (_HYLabel*)components (4);

            selections->Clear();
            theList->GetSelection (*selections);

            _String thisItem = l1->GetText ();

            thisItem = thisItem.Cut (0, thisItem.Find ('(')) & _String ((long)selections->lLength) & " chosen).";
            l1->SetText (thisItem);

            cn->StartDraw();
            if (selections->lLength!=1) {
                if (selections->lLength>1) {
                    _String multipleItems (128L, true);
                    for (k=0; (k<selections->lLength)&&(k<8); k++) {
                        thisItem = *(_String*)dData(selections->lData[k]);
                        if (thisItem.sLength>45) {
                            thisItem = thisItem.Cut (0,44) & "...";
                        }
                        thisItem = thisItem.Replace ("\r"," ",true);
                        thisItem = thisItem.Replace ("\n"," ",true);
                        multipleItems << thisItem;
#ifdef __MAC__
                        multipleItems << '\r';
#else
#ifdef __WINDOZE__
                        multipleItems << "\n\r";
#else
                        multipleItems << '\n';
#endif
#endif
                    }
                    if (selections->lLength>8) {
                        multipleItems << "...";
                    }

                    multipleItems.Finalize();
                    if (!multipleItems.Equal(&lastString)) {
                        cn->EraseRect(hierListTextBoxBounds);
                        //cn->CrossfadeText (lastString,multipleItems,hierListTextBoxBounds,8,30,HY_ALIGN_LEFT,2);
                        cn->DisplayText (multipleItems,hierListTextBoxBounds,HY_ALIGN_LEFT);
                        lastString = multipleItems;
                    }
                } else {
                    cn->EraseRect(hierListTextBoxBounds);
                    lastString = empty;
                }

            } else {
                _String * thisString = (_String*)dData(selections->lData[0]);
                //cn->CrossfadeText (lastString,*thisString,hierListTextBoxBounds,6,30,HY_ALIGN_LEFT,2);
                cn->DisplayText (*thisString,hierListTextBoxBounds,HY_ALIGN_LEFT);
                lastString = *thisString;
            }
            if ((selections->lLength == reqSel)||((reqSel==0)&&(selections->lLength))) {
                b1->EnableButton (true);
            } else {
                b1->EnableButton (false);
            }

            cn->EndDraw();
            cn->_MarkForUpdate();
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
            _HYButton* b1      = (_HYButton*)components (2);
            if (b1->isEnabled) {
                ProcessEvent (generateButtonPushEvent (b1->GetID(),1));
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

//____________________________________________________________________________________________

void        _HYSimpleListSelectDialog::SetInitialSelection (void)
{
    long     index = 0;
    _HYTable * theList = (_HYHList*)components  (0);
    _String* sel   = RetrieveCachedSelection (&GetTitle());
    if (sel) {
        index = theList->FindString (sel);
        if (index<0) {
            index = 0;
        }
    }
    _SimpleList  ns         (index);
    theList->ClearSelection (false);
    theList->SetSelection   (ns,true);
    theList->_MarkCellsForUpdate (ns);
    theList->ScrollToRow         (index);
    ProcessEvent (generateKeyboardFocusEvent (theList->GetID()));
}

//__________________________________________________________________
_HYPreferencesDialog::_HYPreferencesDialog   (_List& thePreferences, _String n, bool, bool* r):_HYTWindow (n,false,true)
{
    _List               hData,
                        rData,
                        iData,
                        diData;

    _SimpleList*        preferencesCodes;

    _List               *cellNames,
                        *cellDescriptions,
                        *cellValues;

    prefList         =  &thePreferences;
    result           =  r;


    preferencesCodes = (_SimpleList*) thePreferences (0);
    cellNames        = (_List*)       thePreferences (1);
    cellDescriptions = (_List*)       thePreferences (2);
    cellValues       = (_List*)       thePreferences (4);

    long            counter;

    for (counter = 0; counter< preferencesCodes->lLength; counter++) {
        if (preferencesCodes->lData[counter]) {
            if (rData.lLength) {
                rData && & iData;
                hData && & rData;
                dData && & diData;
            }
            rData.Clear();
            iData.Clear();
            diData.Clear();
            diData << (*cellDescriptions) (counter);
            rData  << (*cellNames) (counter);
        } else {
            iData  << (*cellNames) (counter);
            diData << (*cellDescriptions) (counter);
        }

        backupValues && (*cellValues) (counter);
    }

    rData && & iData;
    hData && & rData;
    dData && & diData;

    _HYRect         canvasSettings = {250,225,250,225,HY_COMPONENT_V_SCROLL|HY_COMPONENT_TRANSP_BG};

    _HYHList*       hl      = new _HYHList (canvasSettings, GetOSWindowData(),hData);

    hl->selectionType |= HY_TABLE_SINGLE_SELECTION|HY_TABLE_FOCUSABLE;

    canvasSettings.top    = canvasSettings.bottom = 150;
    canvasSettings.left   = canvasSettings.right = 175;

    canvasSettings.width  = HY_COMPONENT_NO_SCROLL|HY_COMPONENT_TRANSP_BG;
    _HYCanvas*      cn    = new _HYCanvas (canvasSettings,GetOSWindowData(),150,175,32);

    canvasSettings.bottom = canvasSettings.top   = 40;
    canvasSettings.left   = canvasSettings.right = 95;

    _HYButton*      b1      = new _HYButton (canvasSettings, GetOSWindowData());
    canvasSettings.left   = canvasSettings.right = 80;
    _HYButton*      b2      = new _HYButton (canvasSettings, GetOSWindowData());

    canvasSettings.left   = canvasSettings.right  = 175;
    canvasSettings.top    = canvasSettings.bottom = 47;

    _HYTextBox*     tb    = new _HYTextBox  (canvasSettings, GetOSWindowData());
    _HYPullDown*    pd    = new _HYPullDown (canvasSettings, GetOSWindowData());
    _HYLabel*       l3    = new _HYLabel    (canvasSettings, GetOSWindowData());

    canvasSettings.top    = canvasSettings.bottom = 20;

    _HYLabel*       l6    = new _HYLabel    (canvasSettings, GetOSWindowData());

    canvasSettings.left  = canvasSettings.right  = 225;
    canvasSettings.top   = canvasSettings.bottom = 10;
    canvasSettings.width = HY_COMPONENT_BORDER_B|HY_COMPONENT_TRANSP_BG;

    _HYLabel*       l1      = new _HYLabel  (canvasSettings, GetOSWindowData());
    canvasSettings.width = HY_COMPONENT_BORDER_T|HY_COMPONENT_TRANSP_BG;
    _HYLabel*       l2      = new _HYLabel  (canvasSettings, GetOSWindowData());

    canvasSettings.left  = canvasSettings.right  = 175;
    canvasSettings.top   = canvasSettings.bottom = 10;
    canvasSettings.width = HY_COMPONENT_NO_SCROLL|HY_COMPONENT_TRANSP_BG;

    _HYLabel*       l4      = new _HYLabel  (canvasSettings, GetOSWindowData());
    _HYLabel*       l5      = new _HYLabel  (canvasSettings, GetOSWindowData());


    hl->SetMessageRecipient (this);
    b1->SetMessageRecipient (this);
    b2->SetMessageRecipient (this);
    tb->SetMessageRecipient (this);
    pd->SetMessageRecipient (this);

    AddObject (hl);// 0
    AddObject (cn);// 1
    AddObject (b1);// 2
    AddObject (b2);// 3
    AddObject (l1);// 4
    AddObject (l2);// 5
    AddObject (l3);// 6
    AddObject (l4);// 7
    AddObject (l5);// 8
    AddObject (l6);// 9
    AddObject (tb);// 10
    AddObject (pd);// 11

    SetTableDimensions (6,3);
    SetCell   (0,0,l1);
    SetCell   (0,1,l4);
    SetCell   (0,2,l4);

    SetCell   (1,0,hl);
    SetCell   (1,1,cn);
    SetCell   (1,2,cn);

    SetCell   (2,0,hl);
    SetCell   (2,1,l6);
    SetCell   (2,2,l6);

    SetCell   (3,0,hl);
    SetCell   (3,1,l3);
    SetCell   (3,2,l3);

    SetCell   (4,0,hl);
    SetCell   (4,1,b1);
    SetCell   (4,2,b2);

    SetCell   (5,0,l2);
    SetCell   (5,1,l5);
    SetCell   (5,2,l5);

    /*b1->SetBackColor (bgRGB);
    b2->SetBackColor (bgRGB);
    l1->SetBackColor (bgRGB);
    l2->SetBackColor (bgRGB);
    l3->SetBackColor (bgRGB);
    l4->SetBackColor (bgRGB);
    l5->SetBackColor (bgRGB);
    l6->SetBackColor (bgRGB);
    tb->SetBackColor (bgRGB);
    pd->SetBackColor (bgRGB);*/


    _HYFont         labelFont;
    SetDefaultDialogFont (labelFont);

    b1->SetFont (labelFont);
    b2->SetFont (labelFont);

#ifdef __MAC__
    labelFont.face = "Geneva";
    labelFont.size = 10;
#endif

    labelFont.style = HY_FONT_PLAIN;
    tb->SetFont (labelFont);

    labelFont.style = HY_FONT_BOLD;

    l6->SetFont (labelFont);
    l6->SetAlignFlags (HY_ALIGN_LEFT);
    pd->SetAlignFlags (HY_ALIGN_LEFT);

    cn->StartDraw ();
    //cn->SetBColor (bgRGB);
    cn->SetDialogBG();
    cn->EraseAll  ();
    cn->SetColor  ((_HYColor) {
        0,0,0
    });

    canvasSettings.width  = 1;
    cn->SetFont (labelFont);
    canvasSettings.left   = 5;
    canvasSettings.top    = 5;
    canvasSettings.bottom = 148;
    canvasSettings.right  = 167;
    cn->DrawInfoBox (canvasSettings,"Item Description");
    labelFont.style = HY_FONT_PLAIN;
    cn->SetFont (labelFont);
    cn->EndDraw();

    b1->SetAlignFlags (HY_ALIGN_RIGHT);

    b1->SetText       ("  OK  ");
    b1->SetButtonKind (HY_BUTTON_OK);
    b2->SetText       (" Cancel ");
    b2->SetButtonKind (HY_BUTTON_CANCEL);



    _HYRect     dim = MinMaxWindowDimensions();

    SetWindowRectangle  (0,0,dim.bottom,dim.right);
    CenterWindow        (this);

    tb->SetVisibleSize (l3->rel);
    pd->SetVisibleSize (l3->rel);
    componentL.lData[10] = componentL.lData[6];
    componentL.lData[11] = componentL.lData[6];
    componentR.lData[10] = componentR.lData[6];
    componentR.lData[11] = componentR.lData[6];
    componentT.lData[10] = componentT.lData[6];
    componentT.lData[11] = componentT.lData[6];
    componentB.lData[10] = componentB.lData[6];
    componentB.lData[11] = componentB.lData[6];

    keyboardFocusChain << 0;
    ProcessEvent (generateKeyboardFocusEvent (hl->GetID()));

#ifndef __MAC__
    _HYRect dummy = {30000,30000,30100,30100,0};
    pd->SetVisibleSize(dummy);
    tb->SetVisibleSize(dummy);
    tb->SetText (" ");
#endif

    DeleteObject (b1);
    DeleteObject (b2);
    DeleteObject (l1);
    DeleteObject (hl);
    DeleteObject (cn);
    DeleteObject (l2);
    DeleteObject (l3);
    DeleteObject (l4);
    DeleteObject (l5);
    DeleteObject (l6);
    DeleteObject (tb);
    DeleteObject (pd);

}

//____________________________________________________________________________________________

void        _HYPreferencesDialog::SetInitialSelection (void)
{
    long       index = 0;
    _HYHList * theList = (_HYHList*)components  (0);
    _String* sel   = RetrieveCachedSelection (&GetTitle());

    if (sel) {
        index = theList->FindString (sel);
        if (index<0) {
            index = 0;
        } else {
            long rubrikIndex = (index >> 16);
            theList->ModifySelection (0,rubrikIndex,false,false,false);
            index = rubrikIndex + (index&0xffff) + 1;
        }
    }

    _SimpleList  ns (2*index+1);
    theList->ClearSelection      (false);
    theList->SetSelection        (ns,true);
    theList->_MarkCellsForUpdate (ns);
    theList->ScrollToRow         (index);
    ProcessEvent (generateKeyboardFocusEvent (theList->GetID()));
}


//__________________________________________________________________
bool    _HYPreferencesDialog::ProcessEvent   (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    long        i,f,g,k;

    _HYHList * theList = (_HYHList*)components  (0);
    _HYCanvas* cn      = (_HYCanvas*)components (1);

    _SimpleList s;
    theList->GetSelection (s);

    _List       *values = (_List*)          (*prefList)(4);
    _SimpleList *iTypes = (_SimpleList*)    (*prefList)(3);
    _List       *choices = (_List*)         (*prefList)(5);

    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++)
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }

        _SimpleList s;
        theList->GetSelection(s);
        if (s.lLength==1) {
            StoreDialogSelection (&GetTitle(),(_String*)theList->GetCellData (s.lData[0]%2,s.lData[0]/2));
        }

        if (i==2) { // OK
            if (s.lLength == 1) {
                k = theList->RubrikIndex (s.lData[0]/2);
                f = k>>16;
                g = k&0xffff;
                k = 0;
                for (f=f-1; f>=0; f--) {
                    k+=((_List*)dData(f))->lLength;
                }
            }
            *result = true;
            postWindowCloseEvent (GetID());
        } else if (i==3) { // Cancel
            *result = false;
            for (k=0; k<backupValues.lLength; k++) {
                values->Duplicate (&backupValues);
            }

            if (setFont) {
                SetPreferences();
            }
            postWindowCloseEvent (GetID());
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
            cn->StartDraw();
            if (s.lLength!=1) {
                _HYLabel *l3 = (_HYLabel*)GetObject(6);
                cn->EraseRect(hierListTextBoxBounds2);
                if (GetCellObject (3,2)!=l3) {
#ifndef __MAC__
                    _HYRect dummy = {30000,30000,30100,30100,0},
                            srel  = ((_HYComponent*)GetCellObject (3,2))->rel;
                    l3->SetVisibleSize(srel);
                    ((_HYComponent*)GetCellObject (3,2))->SetVisibleSize (dummy);
#endif
                    SetCell (3,1,l3);
                    SetCell (3,2,l3);
                }
            } else {
                k = theList->RubrikIndex (s.lData[0]/2);
                f = k>>16;
                g = k&0xffff;
                _String * thisString = (_String*)(*(_List*)(dData(f)))(g);

                _HYLabel * vlabel = (_HYLabel*)GetObject (9);

                if (g) {
                    vlabel->SetText ("  Value:");
                    k = theList->AbsoluteIndex (f,g);

                    if (iTypes->lData[k] == PREFITEM_TEXTBOX) {
                        _HYTextBox* tb = (_HYTextBox*)GetObject(10);
                        if (GetCellObject (3,2)!=tb) {
#ifndef __MAC__
                            _HYRect dummy = {30000,30000,30100,30100,0},
                                    srel  = ((_HYComponent*)GetCellObject (3,2))->rel;
                            tb->SetVisibleSize(srel);
                            ((_HYComponent*)GetCellObject (3,2))->SetVisibleSize (dummy);
#endif
                            SetCell (3,1,tb);
                            SetCell (3,2,tb);
                        }
                        tb->EnableTextEdit(true);
                        tb->SetText (*(_String*)(*values)(k),true);
                        tb->SetSelection (0,30000);
                        tb->Activate();
                        keyboardFocusChain << 10;
                    } else {
                        _HYPullDown* pd = (_HYPullDown*)GetObject(11);
                        if (GetCellObject (3,2)!=pd) {
#ifndef __MAC__
                            _HYRect dummy = {30000,30000,30100,30100,0},
                                    srel  = ((_HYComponent*)GetCellObject (3,2))->rel;
                            pd->SetVisibleSize(srel);
                            ((_HYComponent*)GetCellObject (3,2))->SetVisibleSize (dummy);
#endif
                            SetCell (3,1,pd);
                            SetCell (3,2,pd);
                        }
                        pd->DeleteAllItems();
                        _List * iC = (_List*)(*choices)(k);

                        for (f = 0; f<iC->lLength; f++) {
                            pd->AddMenuItem (*((_String*)(*iC)(f)),-1);
                        }

                        _String * cV = (_String*)(*values)(k);

                        f = pd->FindMenuItem (*cV);
                        pd->SetVisibleSize  (pd->rel);
                        pd->ChangeSelection (f,false);
                        pd->Activate();
                        keyboardFocusChain.Delete (1);
                    }
                } else {
                    vlabel->SetText (empty);
                    _HYLabel* l3 = (_HYLabel*)GetObject(6);
                    if (GetCellObject (3,2)!=l3) {
#ifndef __MAC__
                        _HYRect dummy = {30000,30000,30100,30100,0},
                                srel  = ((_HYComponent*)GetCellObject (3,2))->rel;
                        l3->SetVisibleSize(srel);
                        ((_HYComponent*)GetCellObject (3,2))->SetVisibleSize (dummy);
#endif
                        SetCell (3,1,l3);
                        SetCell (3,2,l3);
                        l3->_MarkForUpdate();
                        keyboardFocusChain.Delete (1);
                    }
                }


                cn->EraseRect (hierListTextBoxBounds2);
                cn->DisplayText (*thisString,hierListTextBoxBounds2,HY_ALIGN_LEFT);
            }
            cn->EndDraw();
            cn->_MarkForUpdate();
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
            k = theList->RubrikIndex (s.lData[0]/2);
            f = k>>16;
            g = k&0xffff;
            theList->HandleKeyMove (theList->IsRubrikOpen(f)?2:3,false);
            done = true;
        }
    } else if (e->EventClass()==_hyMenuSelChangeEvent) {
        f = e->EventCode().Find(',');
        firstArg = e->EventCode().Cut(0,f-1);
        i = MatchComponentID (firstArg);
        if (i==11) {
            f = e->EventCode().Cut(f+1,-1).toNum();
            k = theList->RubrikIndex (s.lData[0]/2);
            g = k&0xffff;
            i = k>>16;
            k = theList->AbsoluteIndex (i,g);
            _List * iChoices = (_List*)(*choices)(k);
            values->Replace (k, (_String*)(*iChoices)(f),true);
            if ((i==0)&&setFont) {
                SetPreferences();
            }
        }
    } else if (e->EventClass()==_hyTextEditChange ) {
        f = e->EventCode().Find(',');
        firstArg = e->EventCode().Cut(0,f-1);
        i = MatchComponentID (firstArg);
        if (i==10) {
            k = theList->RubrikIndex (s.lData[0]/2);
            g = k&0xffff;
            k = k>>16;
            k = theList->AbsoluteIndex (k,g);
            _String * tbs;
            ((_HYTextBox*)GetObject(10))->StoreText (tbs);
            values->Replace (k,tbs,false);
        }
    }

    if (done) {
        DeleteObject (e);
        return true;
    }
    return _HYTWindow::ProcessEvent(e);
}

//____________________________________________________________________________________________

long  HandleHierListSelection (_List& data, _SimpleList& choices, _SimpleList& validChoices, _String titleInfo, _SimpleList& selections, long )
{
    long res = 0;
    _HYListSelectDialog* sd = new _HYListSelectDialog
    (&data, &choices, &validChoices,titleInfo ,&selections,&res);
    sd->SetInitialSelection ();
    sd->BringToFront();
    while (windowObjectRefs.Find ((long)sd)>=0) {
        handleGUI();
    }
    return res;
}


//____________________________________________________________________________________________
long  HandleListSelection (_List& data, _SimpleList& choices, _SimpleList& validChoices, _String titleInfo, _SimpleList& selections, long fixedLength, Ptr prt)
{
    long res = -1;
    if (data.lLength < 1 || validChoices.lLength < 1) {
        _String errMsg ("An empty list of choices was passed to 'HandleListSelection'");
        ProblemReport (errMsg);
    } else {
        _HYSimpleListSelectDialog* sd = new _HYSimpleListSelectDialog
        (&data, &choices, &validChoices,titleInfo ,&selections,fixedLength,&res, prt);
        sd->SetInitialSelection ();
        sd->BringToFront();
        while (windowObjectRefs.Find ((long)sd)>=0) {
            handleGUI();
        }
    }
    return res;
}

//____________________________________________________________________________________________
long  HandleListSelection (_List& data, _String titleInfo, Ptr prt)
{
    _SimpleList validChoices,
                choices,
                sels;

    _List       menuData;

    validChoices << 0;
    validChoices << 1;

    for (long k=0; k<data.lLength; k+=2) {
        _List aChoice;

        aChoice << data (k);
        aChoice << data (k+1);

        menuData && & aChoice;

        choices << choices.lLength;
    }

    return  HandleListSelection (menuData, validChoices, choices, titleInfo, sels, 1, prt);
}


//____________________________________________________________________________________________

bool        EnterStringDialog (_String& res, _String& prompt, Ptr parent, _hyStringValidatorType validator)
{
    bool resB = false;
    _HYTextDialog* sd = new _HYTextDialog (&res,prompt,&resB, parent, validator);
    sd->BringToFront();
    while (windowObjectRefs.Find ((long)sd)>=0) {
        handleGUI();
    }
    return resB;
}

//____________________________________________________________________________________________

bool        EnterString2Dialog (_String& res, _String& res2, _String& prompt, _String& prompt2, Ptr parent)
{
    bool resB = false;
    _HYTextDialog2* sd = new _HYTextDialog2 (&res,&res2,prompt,prompt2,&resB, parent);
    sd->BringToFront();
    while (windowObjectRefs.Find ((long)sd)>=0) {
        handleGUI();
    }
    return resB;
}

//____________________________________________________________________________________________

bool        EnterStringDialogWithCheckbox (_String& res, _String& prompt, _String& cprompt, bool& resC, Ptr parent)
{
    bool resB = false;
    _HYTextDialogWithCheckbox* sd = new _HYTextDialogWithCheckbox (&res,prompt,cprompt,&resB,&resC, parent);
    sd->BringToFront();
    while (windowObjectRefs.Find ((long)sd)>=0) {
        handleGUI();
    }
    return resB;
}

//____________________________________________________________________________________________

bool        EnterStringDialogWithPulldown (_String& res, _String& prompt, _String& cprompt, long& msel, _List& pulldown, _List * actions, bool& resC, long initChoice, Ptr parent)
{
    bool resB = false;
    _HYTextDialogWithPulldown* sd = new _HYTextDialogWithPulldown (&res,prompt,cprompt, pulldown, actions, &resB,&resC,&msel,initChoice,parent);
    sd->BringToFront();
    while (windowObjectRefs.Find ((long)sd)>=0) {
        handleGUI();
    }
    return resB;
}


//____________________________________________________________________________________________
bool  HandlePreferences (_List& data, _String title, bool doFonts)
{
    bool res = false;
    _HYPreferencesDialog* sd = new _HYPreferencesDialog (data, title, doFonts, &res);
    sd->SetInitialSelection ();
    sd->BringToFront();
    while (windowObjectRefs.Find ((long)sd)>=0) {
        handleGUI();
    }
    return res;
}

//____________________________________________________________________________________________
bool  ProceedPrompt (_String& pr, Ptr par)
{
    bool res = false;

    _String b1 (" OK "),
            b2 (" Cancel ");

    _HYFont def;
    SetDefaultDialogFont (def);

    _HYProceedPromptBox * sd = new _HYProceedPromptBox (pr,b1,b2,def,133,&res,par);
    sd->BringToFront();
    while (windowObjectRefs.Find ((long)sd)>=0) {
        handleGUI();
    }
    return res;
}

//____________________________________________________________________________________________
void  ProblemReport (_String& pr, Ptr par)
{
    bool res = false;

    _String b1,
            b2;

    if (wiseCrackButtons.lLength == 0) {
        b1 = " Oh Well ";
        wiseCrackButtons && & b1;
        b1 = " Too bad ";
        wiseCrackButtons && & b1;
        b1 = " C'est la Vie ";
        wiseCrackButtons && & b1;
        b1 = " Argh! ";
        wiseCrackButtons && & b1;
        b1 = " Begone! ";
        wiseCrackButtons && & b1;
        b1 = " As if I had a choice! ";
        wiseCrackButtons && & b1;
        b1 = " Really? ";
        wiseCrackButtons && & b1;
        b1 = " Not again! ";
        wiseCrackButtons && & b1;
        b1 = " Why, oh why? ";
        wiseCrackButtons && & b1;
        b1 = " Blast! ";
        wiseCrackButtons && & b1;
        b1 = " Rats! ";
        wiseCrackButtons && & b1;
        b1 = " Blimey! ";
        wiseCrackButtons && & b1;
        b1 = " This is the end... ";
        wiseCrackButtons && & b1;
        b1 = " Whatever, dude ";
        wiseCrackButtons && & b1;
        b1 = " What have I done to deserve this? ";
        wiseCrackButtons && & b1;
        b1 = " Hasta la vista! ";
        wiseCrackButtons && & b1;
        b1 = " Famous last words ";
        wiseCrackButtons && & b1;
    }

    _Constant       messageListLength (wiseCrackButtons.lLength);
    _PMathObj       rn = _Constant (0.0).Random (&messageListLength);
    b1 = *(_String*)wiseCrackButtons (rn->Value());
    DeleteObject (rn);

    _HYFont def;
    SetDefaultDialogFont (def);

    _HYProceedPromptBox * sd = new _HYProceedPromptBox (pr,b1,b2,def,130,&res,par);
    sd->BringToFront();
    while (windowObjectRefs.Find ((long)sd)>=0) {
        handleGUI();
    }
}

//____________________________________________________________________________________________
bool  ProceedPromptWithCheck (_String& pr, _String& cpr, bool& cs, Ptr par)
{
    bool res = false;

    _String b1 (" OK "),
            b2 (" Cancel ");

    _HYFont def;
    SetDefaultDialogFont (def);

    _HYProceedPromptBoxWCheck * sd = new _HYProceedPromptBoxWCheck (pr,b1,b2,cpr,def,133,&res,&cs,par);
    sd->BringToFront();
    while (windowObjectRefs.Find ((long)sd)>=0) {
        handleGUI();
    }
    return res;
}

//__________________________________________________________________
_String*        RetrieveCachedSelection (_String* title)
{
    long        f = cachedDialogTitles.BinaryFind (title);

    if (f>=0) {
        _String * sel = (_String*)cachedDialogSelections(f);
        if (sel->sLength) {
            return sel;
        }
    }
    return nil;
}

//__________________________________________________________________
void            StoreDialogSelection    (_String* title,_String* selection)
{
    long        f = cachedDialogTitles.BinaryFind (title);

    if (f<0) {
        f = cachedDialogTitles.BinaryInsert (title);
        cachedDialogSelections.InsertElement (selection,f,true);
    } else {
        cachedDialogSelections.Replace (f,selection,true);
    }
}


//____________________________________________________________________________________________

_String  NewTreeWindow (long sourceDF)
{
    // set up default preferences
    _List   theList,dsNames;
    _String optionList ("No dataset"), comma(",");
    long    k;

    AddItemToPreferences (1|8,-1,"Tree Identifier","Choose an identifier for the new tree variable.","",nil,theList,false);
    AddItemToPreferences (0,PREFITEM_TEXTBOX,"Identifier","This must be a valid unique HYPHY identifier (begins with a letter or an underscore, containes letters, numbers or underscores).","New_Tree",nil, theList,false);
    AddItemToPreferences (1|8,-1,"Tips and Labels","Specifies whether the new tree will be created based on the sequences of an existing data set, or simply with a given number of tips. Also selects the style for the new tree.","",nil,theList,false);
    AddItemToPreferences (0,PREFITEM_TEXTBOX,"Number of Tips","Use this option to create a tree with a desired number of tips, without linking to an existing data file.","10",nil,theList,false);
    if (sourceDF<0) {
        dsNames && &optionList;
        for (k=0; k<dataSetNamesList.lLength; k++) {
            dsNames<<dataSetNamesList(k);
        }
    } else {
        dsNames<<dataSetFilterNamesList(sourceDF);
    }
    AddItemToPreferences (0,PREFITEM_POPUP,"Source Dataset","Select a dataset from the list (you need to open a data set for it to appear in the list).\"No dataset\" indicates that a tree should with a given number of tips will be created.","No dataset",&dsNames,theList,false);
    AddItemToPreferences (1|8,-1,"Tree Options","Selects starting topology for the new tree, and (optionally) defines assigned branch lengths","",nil,theList,false);
    optionList = "Star,Ladder Left,Ladder Right,Balanced";
    AddItemToPreferences (0,PREFITEM_POPUP,"Tree Style","Selects a topology for the new tree. It can be edited later.","Star",optionList.Tokenize (comma),theList,true);
    AddItemToPreferences (0,PREFITEM_TEXTBOX,"Default Length","Assign default length to all branches (the value -1 is equivalent to no assigned branch lengths).","-1",nil,theList,false);
    AddItemToPreferences (0,PREFITEM_TEXTBOX,"Newick String","Enter a Newick string here to define tree topology. If anything is entered in this field, all other options will be ignored","",nil,theList,false);
    if (HandlePreferences (theList, "New Tree Setup", false)) {
        _List   *treeSettings = (_List*)theList.lData[4],
                 *settingsValues = (_List*)theList.lData[5];

        _String *treeName = (_String*)treeSettings->lData[1],
                 *setName = (_String*)treeSettings->lData[4],
                  *newickString = (_String*)treeSettings->lData[8];

        if (!treeName->IsValidIdentifier()) {
            optionList = _String('"')& *treeName & "\" appears to be an invalid identifier. Please try again.";
            ProblemReport (optionList);
            return empty;
        }
        if (LocateVarByName (*treeName)>=0) {
            optionList = _String('"')& *treeName & "\" already exists. Kill and replace?";
            if (!ProceedPrompt (optionList)) {
                return empty;
            }
            DeleteVariable (*treeName,true);
        }
        long    tipCount = ((_String*)treeSettings->lData[3])->toNum();
        dsNames.Clear();
        k = ((_List*)(*settingsValues)(4))->Find(setName);
        if ((k>=1)||(sourceDF>=0)) { // use data set
            if (sourceDF<0) {
                dsNames.Duplicate (&((_DataSet*)dataSetList(k-1))->GetNames());
            } else {
                _DataSetFilter* sDF = (_DataSetFilter*)dataSetFilterList(sourceDF);
                _DataSet*   sDS = sDF->GetData();
                for (k=0; k<sDF->theNodeMap.lLength; k++) {
                    dsNames<< (sDS->GetNames())(sDF->theNodeMap.lData[k]);
                }
            }
        } else {
            if (tipCount<=1) {
                _String errMsg = _String("Sorry, but I can't create a tree with ")&tipCount&" tips for fairly obvious reasons.";
                ProblemReport (errMsg);
                return empty;
            }

            for (long f=1; f<=tipCount; f++) {
                _String thisSpec ("Species_");
                thisSpec = thisSpec & f;
                dsNames&& &thisSpec;
            }
        }
        tipCount = dsNames.lLength;
        _String treeString ((unsigned long)32,true);
        setName = (_String*)treeSettings->lData[6];
        k = ((_List*)(*settingsValues)(6))->Find(setName);
        setName = (_String*)treeSettings->lData[7];
        _Parameter defL = setName->toNum();

        if (newickString->sLength==0) {
            if (k<=0) {
                treeString<<'(';
                for (k=0; k<tipCount; k++) {
                    if (k) {
                        treeString << ',';
                    }
                    treeString << (_String*)dsNames(k);
                }
                treeString<<')';
            } else if (k==1) { // ladder left
                treeString<<'(';
                treeString << (_String*)dsNames(0);
                treeString<<',';
                treeString << (_String*)dsNames(1);
                for (k=2; k<tipCount-1; k++) {
                    treeString<<',';
                    treeString<<'(';
                    treeString << (_String*)dsNames(k);
                }
                if (tipCount>2) {
                    treeString<<',';
                    treeString << (_String*)dsNames(tipCount-1);
                    for (k=2; k<tipCount-1; k++) {
                        treeString<<')';
                    }
                }
                treeString<<')';

            } else if (k==2) { // ladder right
                treeString<<'(';
                if (tipCount>2) {
                    for (k=2; k<tipCount-1; k++) {
                        treeString<<'(';
                    }
                }
                treeString << (_String*)dsNames(0);
                treeString<<',';
                treeString << (_String*)dsNames(1);
                for (k=2; k<tipCount-1; k++) {
                    treeString<<')';
                    treeString<<',';
                    treeString << (_String*)dsNames(k);
                }
                if (tipCount>2) {
                    treeString<<',';
                    treeString << (_String*)dsNames(k);
                }
                treeString<<')';
            } else if (k==3) { // balanced
                BuildBalancedTree (0,dsNames.lLength-1,treeString,dsNames);
            }
            treeString.Finalize();
        } else {
            treeString.Finalize();
            treeString = *newickString;
        }
        if (treeString.sLength) {
            _HYTreePanel * newTree = new _HYTreePanel (*treeName,treeString);
            if (defL>-1.0) {
                _TheTree* me = newTree->LocateMyTreeVariable();
                _Constant newVal (defL);
                if (me) {
                    _CalcNode* travNode = me->DepthWiseTraversal(TRUE);
                    while (!me->IsCurrentNodeTheRoot()) {
                        travNode->SetValue (&newVal);
                        travNode = me->DepthWiseTraversal();
                    }
                    newTree->UpdateScalingVariablesList();
                }
            }
            newTree->BringToFront();
            newTree->Show();
        }
        return *treeName;
    }
    return empty;
}

//____________________________________________________________________________________________

void  BuildBalancedTree (long start, long end, _String& result, _List& tipNames)
{
    long rangeCount = end-start+1;
    if (rangeCount==1) {
        result<<(_String*)tipNames(start);
    } else if (rangeCount==2) {
        result<<'(';
        result<<(_String*)tipNames(start);
        result<<',';
        result<<(_String*)tipNames(end);
        result<<')';
    } else {
        result<<'(';
        if (rangeCount%2) {
            BuildBalancedTree (start,start+rangeCount/2,result,tipNames);
            result<<',';
            BuildBalancedTree (start+rangeCount/2+1,end,result,tipNames);
        } else {
            BuildBalancedTree (start,start+rangeCount/2-1,result,tipNames);
            result<<',';
            BuildBalancedTree (start+rangeCount/2,end,result,tipNames);
        }
        result<<')';
    }
}


//__________________________________________________________________
_HYCombDialog::_HYCombDialog     (_SimpleList* storage, bool* res, Ptr windowRf):_HYTWindow ("Define the comb",false,true, windowRf)
{

    combResult   = storage;
    combResult->Clear();
    result       = res;

    if (combSelectionMemory.lLength==0) {
        combSelectionMemory << 0;
        combSelectionMemory << 0;
        combSelectionMemory << 0;
    }

    _HYRect         canvasSettings = {30,80,30,80,HY_COMPONENT_TRANSP_BG};

    _HYLabel*       l1 = new _HYLabel (canvasSettings, GetOSWindowData());
    checkPointer    (l1);

    canvasSettings.left   = canvasSettings.right  = 270;
    _HYPullDown*    p1    = new _HYPullDown (canvasSettings, GetOSWindowData());
    checkPointer    (p1);

    long            k;

    _HYCheckbox*    cb[10];

    canvasSettings.left = canvasSettings.right = 35;
    for (k=0; k<10; k++) {
        cb[k] = new _HYCheckbox (canvasSettings, GetOSWindowData());
        cb[k]->SetMessageRecipient (this);
        checkPointer (cb[k]);
        AddObject (cb[k]);
    }

    canvasSettings.bottom = canvasSettings.top   = 40;
    canvasSettings.left   = canvasSettings.right = 260;

    _HYButton*      b1      = new _HYButton (canvasSettings, GetOSWindowData());
    canvasSettings.left   = canvasSettings.right = 90;
    _HYButton*      b2      = new _HYButton (canvasSettings, GetOSWindowData());

    checkPointer (b1);
    checkPointer (b2);

    p1->SetMessageRecipient (this);
    b1->SetMessageRecipient (this);
    b2->SetMessageRecipient (this);

    AddObject (b1);// 10
    AddObject (b2);// 11
    AddObject (l1);// 12
    AddObject (p1);// 13

    SetTableDimensions (3,10);

    for (k=0; k<10; k++) {
        SetCell (1,k,cb[k]);
    }

    for (k=0; k<3; k++) {
        SetCell (0,k,l1);
    }
    for (; k<10; k++) {
        SetCell (0,k,p1);
    }

    for (k=0; k<8; k++) {
        SetCell (2,k,b1);
    }
    for (; k<10; k++) {
        SetCell (2,k,b2);
    }

    _HYFont         labelFont;
    SetDefaultDialogFont (labelFont);

    b1->SetFont (labelFont);
    b2->SetFont (labelFont);
    l1->SetFont (labelFont);

    b1->SetText (" OK ");
    b2->SetText (" Cancel ");

    b1->SetAlignFlags (HY_ALIGN_RIGHT);
    l1->SetAlignFlags (HY_ALIGN_LEFT);
    p1->SetAlignFlags (HY_ALIGN_LEFT);
    l1->SetText ("Comb Size:");

    for (k=0; k<10; k++) {
        cb[k]->SetFont (labelFont);
        cb[k]->SetText (_String (k+1));
        cb[k]->SetSpacing (1);
    }

    for (k=2; k<=10; k++) {
        p1->AddMenuItem (_String (k), -1);
    }

    combResult->Duplicate (&combSelectionMemory);
    p1->ChangeSelection (combSizeMemory-2, true);

    for (k=combSelectionMemory.lLength; k<10; k++) {
        cb[k]->Enable (false);
    }

    for (k=0; k<combSelectionMemory.lLength; k++) {
        cb[k]->SetState (combSelectionMemory.lData[k]);
    }

    _HYRect     dim = MinMaxWindowDimensions();

    b1->SetButtonKind (HY_BUTTON_OK);
    b2->SetButtonKind (HY_BUTTON_CANCEL);
    b1->EnableButton ((combResult->Find(0)>=0)&&(combResult->Find(1)>=0));

    SetWindowRectangle  (0,0,dim.bottom,dim.right);
    CenterWindow        (this);

    DeleteObject (b1);
    DeleteObject (b2);
    DeleteObject (l1);
    DeleteObject (p1);

    for (k=0; k<10; k++) {
        DeleteObject (cb[k]);
    }

}

//__________________________________________________________

bool    _HYCombDialog::ProcessEvent (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    long        i,f,g,k;

    if (e->EventClass()==_hyMenuSelChangeEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        g = e->EventCode().Cut (f+1,-1).toNum();

        for (i=0; i<components.lLength; i++)
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }

        if (i==13) {
            _HYPullDown * p1 = (_HYPullDown*)components(i);
            k = p1->GetSelection()+2;
            if (k<combResult->lLength) {
                for (f=combResult->lLength-1; f>=k; f--) {
                    combResult->Delete (k);
                    ((_HYCheckbox*)components(f))->Enable (false);
                }
            } else if (k>combResult->lLength) {
                for (f=combResult->lLength; f<k; f++) {
                    (*combResult)<<0;
                    ((_HYCheckbox*)components(f))->Enable (true);
                }
            }
        }
        done = true;
    } else if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        for (i=0; i<components.lLength; i++)
            if (((_HYGuiObject*)components(i))->MatchID(k)) {
                break;
            }

        if (i==10) { // OK
            *result = true;
            combSelectionMemory.Duplicate (combResult);
            combSizeMemory = combResult->lLength;
            postWindowCloseEvent (GetID());
        } else if (i==11) {
            *result = false;
            postWindowCloseEvent (GetID());
        } else if (i<10) {
            combResult->lData[i] = !combResult->lData[i];
            ((_HYButton*)components (10))->EnableButton ((combResult->Find(0)>=0)&&(combResult->Find(1)>=0));
        }

        done = true;
    }

    if (done) {
        DeleteObject (e);
        return true;
    }
    return _HYTWindow::ProcessEvent(e);
}


//__________________________________________________________________

_HYRect partitionDialogDescBox  = {15,5,95,295,0},
        partitionDialogTextBox  = {25,20,85,285,0},
        partitionDialogColorBox = {5,5,25,30,1};


//__________________________________________________________________
_HYPartitionDialog::_HYPartitionDialog   (_String& partDesc, _String* pn, _HYColor* pc, bool* dir, char* rf, long* cr, long df, bool* res, char codon, Ptr parentWindow):_HYTWindow ("Partition Properties",false,true, parentWindow)
{

    partName = pn;
    color    = pc;

    direction = dir;
    codeRef   = cr;

    dfID = df;

    readFrame = rf;

    result    = res;

    _HYRect         canvasSettings = {100,300,100,300,HY_COMPONENT_TRANSP_BG};

    _HYCanvas*      cn = new _HYCanvas (canvasSettings, GetOSWindowData(),100,300,32);
    checkPointer   (cn);

    if (codon) {
        canvasSettings.left = canvasSettings.right  = 100;
    } else {
        canvasSettings.left = canvasSettings.right  = 60;
    }

    canvasSettings.top  = canvasSettings.bottom = 30;

    _HYLabel*       l1 = new _HYLabel (canvasSettings, GetOSWindowData());
    checkPointer    (l1);
    _HYLabel*       l2 = new _HYLabel (canvasSettings, GetOSWindowData());
    checkPointer    (l2);

    _HYLabel*       l3,
                    *       l4,
                    *       l5;

    if (codon) {
        l3 = new _HYLabel (canvasSettings, GetOSWindowData());
        checkPointer    (l3);
        l4 = new _HYLabel (canvasSettings, GetOSWindowData());
        checkPointer    (l3);
        if (codon>1) {
            l5 = new _HYLabel (canvasSettings, GetOSWindowData());
            checkPointer    (l5);
        }
    }

    if (codon) {
        canvasSettings.left = canvasSettings.right  = 200;
    } else {
        canvasSettings.left = canvasSettings.right  = 240;
    }

    _HYTextBox*     tb    = new _HYTextBox (canvasSettings, GetOSWindowData());
    checkPointer    (tb);
    _HYCanvas *     cn2   = new _HYCanvas (canvasSettings, GetOSWindowData(),30,175,32);
    checkPointer    (cn2);

    _HYPullDown*    p1,
                    *    p2;

    _HYCheckbox*    cb1,
                    *    cb2;

    if (codon) {
        p1 = new _HYPullDown (canvasSettings, GetOSWindowData());
        checkPointer    (p1);
        if (codon>1) {
            p2 = new _HYPullDown (canvasSettings, GetOSWindowData());
            checkPointer    (p2);
        }
        canvasSettings.left   = canvasSettings.right  = canvasSettings.right/2;
        cb1 = new _HYCheckbox (canvasSettings, GetOSWindowData(),true);
        checkPointer (cb1);
        cb2 = new _HYCheckbox (canvasSettings, GetOSWindowData(),true);
        checkPointer (cb2);
    }


    long            k;


    canvasSettings.bottom = canvasSettings.top   = 40;
    canvasSettings.left   = canvasSettings.right = 215;

    _HYButton*      b1      = new _HYButton (canvasSettings, GetOSWindowData());
    canvasSettings.left   = canvasSettings.right = 85;
    _HYButton*      b2      = new _HYButton (canvasSettings, GetOSWindowData());

    checkPointer (b1);
    checkPointer (b2);



    tb-> SetMessageRecipient (this);
    cn2->SetMessageRecipient (this);
    b1->SetMessageRecipient (this);
    b2->SetMessageRecipient (this);


    AddObject (b1);  // 0
    AddObject (b2);  // 1
    AddObject (cn);  // 2
    AddObject (l1);  // 3
    AddObject (l2);  // 4
    AddObject (tb);  // 5
    AddObject (cn2); // 6

    if (codon) {
        p1->SetMessageRecipient (this);
        cb1->SetMessageRecipient (this);
        cb2->SetMessageRecipient (this);
        AddObject (l3); // 7
        AddObject (l4); // 8
        AddObject (p1); // 9
        AddObject (cb1);// 10
        AddObject (cb2);// 11
        SetTableDimensions (6,3);

        if (codon>1) {
            AddObject (p2);// 12
            AddObject (l5);// 13
            SetTableDimensions (7,3);
        }

        SetCell (0,0,l1);
        SetCell (0,1,tb);
        SetCell (0,2,tb);
        SetCell (1,0,l2);
        SetCell (1,1,cn2);
        SetCell (1,2,cn2);
        SetCell (2,0,l3);
        SetCell (2,1,p1);
        SetCell (2,2,p1);
        SetCell (3,0,l4);
        SetCell (3,1,cb1);
        SetCell (3,2,cb2);

        k = 4;

        if (codon>1) {
            SetCell (4,0,l5);
            SetCell (4,1,p2);
            SetCell (4,2,p2);
            k = 5;
        }

        SetCell (k,0,cn);
        SetCell (k,1,cn);
        SetCell (k,2,cn);
        k++;
        SetCell (k,0,b1);
        SetCell (k,1,b1);
        SetCell (k,2,b2);

    } else {
        SetTableDimensions (4,2);
        SetCell (0,0,l1);
        SetCell (0,1,tb);
        SetCell (1,0,l2);
        SetCell (1,1,cn2);
        SetCell (2,0,cn);
        SetCell (2,1,cn);
        SetCell (3,0,b1);
        SetCell (3,1,b2);
    }

    _HYFont         labelFont;
    SetDefaultDialogFont(labelFont);

    b1->SetFont (labelFont);
    b2->SetFont (labelFont);
    l1->SetFont (labelFont);
    l2->SetFont (labelFont);

    b1->SetText (" OK ");
    b2->SetText (" Cancel ");

    b1->SetAlignFlags (HY_ALIGN_RIGHT);
    l1->SetAlignFlags (HY_ALIGN_LEFT);
    l2->SetAlignFlags (HY_ALIGN_LEFT);
    l1->SetText ("Name:");
    l2->SetText ("Color:");

    if (codon) {
        l3->SetFont (labelFont);
        l4->SetFont (labelFont);
        cb1->SetFont    (labelFont);
        cb2->SetFont    (labelFont);
        l3->SetAlignFlags (HY_ALIGN_LEFT);
        l4->SetAlignFlags (HY_ALIGN_LEFT);
        cb1->SetAlignFlags (HY_ALIGN_LEFT);
        cb2->SetAlignFlags (HY_ALIGN_LEFT);
        p1->SetAlignFlags (HY_ALIGN_LEFT);
        l3->SetText ("Reading Frame:");
        l4->SetText ("Orientation:");
        cb1->SetText ("Normal");
        cb2->SetText ("Reversed");
        cb1->SetState(true,false);

        p1->AddMenuItem (_String("No offset"),-1);
        p1->AddMenuItem (_String("Offset by 1"),-1);

        if (codon>1) {
            p1->AddMenuItem (_String("Offset by 2"),-1);
            p2->SetAlignFlags (HY_ALIGN_LEFT);
            l5->SetFont (labelFont);
            l5->SetAlignFlags (HY_ALIGN_LEFT);
            l5->SetText ("Genetic Code:");

            for (k=0; k<geneticCodes.lLength; k++) {
                p2->AddMenuItem (*(_String*)(*((_List*)geneticCodes(k)))(0),-1);
            }
            p2->ChangeSelection (*codeRef,false);

        }
    }

    b1->SetButtonKind (HY_BUTTON_OK);
    b2->SetButtonKind (HY_BUTTON_CANCEL);

    labelFont.face  = "Times";
    cn->StartDraw();
    cn->SetDialogBG ();
    cn->EraseAll();
    labelFont.style = HY_FONT_BOLD;
    cn->SetFont         (labelFont);
    cn->DrawInfoBox     (partitionDialogDescBox,"Partition Info");
    labelFont.style = HY_FONT_PLAIN;
    cn->SetFont         (labelFont);
    cn->DisplayText     (partDesc,partitionDialogTextBox,HY_ALIGN_LEFT);
    cn->EndDraw();

    tb->SetFont (labelFont);
    tb->SetAlignFlags (HY_ALIGN_LEFT);
    tb->SetText (*(_String*)dataSetFilterNamesList(dfID));

    labelFont.style = HY_FONT_ITALIC;
    cn2->StartDraw();
    cn2->SetDialogBG();
    cn2->EraseAll();
    cn2->SetColor (*color);
    cn2->FillRect (partitionDialogColorBox);
    cn2->SetColor ((_HYColor) {
        0,0,0
    });
    cn2->DrawRect (partitionDialogColorBox);
    cn2->SetFont  (labelFont);
    cn2->DisplayText ("(click to change)",partitionDialogColorBox.bottom, partitionDialogColorBox.right+10,true);
    cn2->EndDraw();
    cn2->SetMouseClick (true);

    _HYRect     dim = MinMaxWindowDimensions();

    SetWindowRectangle  (0,0,dim.bottom,dim.right);
    CenterWindow        (this);

    keyboardFocusChain << 5;
    ProcessEvent (generateKeyboardFocusEvent (tb->GetID()));

    DeleteObject (b1);
    DeleteObject (b2);
    DeleteObject (l1);
    DeleteObject (l2);
    DeleteObject (cn);
    DeleteObject (cn2);
    DeleteObject (tb);

    if (codon) {
        DeleteObject (l3);
        DeleteObject (l4);
        DeleteObject (p1);
        DeleteObject (cb1);
        DeleteObject (cb2);
        if (codon>1) {
            DeleteObject (l5);
            DeleteObject (p2);
        }
    }

}

//__________________________________________________________

bool    _HYPartitionDialog::ProcessEvent (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    long        i,f,g;

    if (e->EventClass()==_hyMenuSelChangeEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        g = e->EventCode().Cut (f+1,-1).toNum();

        i = MatchComponentID (firstArg);

        done = true;
    } else if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);

        if (i==0) { // OK
            *result = true;
            if (components.lLength>11) {
                _HYCheckbox * cb = (_HYCheckbox*)components (11);
                *direction = cb->GetState();
                *readFrame = ((_HYPullDown*)components(9))->GetSelection();
                if (components.lLength>13) {
                    *codeRef = ((_HYPullDown*)components(12))->GetSelection();
                }
            }
            postWindowCloseEvent (GetID());
        } else if (i==1) {
            *result = false;
            postWindowCloseEvent (GetID());
        } else if ((i==10)||(i==11)) {
            _HYCheckbox * cb1 = (_HYCheckbox*)components (10),
                          * cb2 = (_HYCheckbox*)components (11);

            if (i==10) {
                cb2->SetState(!cb1->GetState(),false);
            } else {
                cb1->SetState(!cb2->GetState(),false);
            }
        }


        done = true;
    } else if (e->EventClass()==_hyTextEditChange) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);

        if (i==5) {
            *partName = ((_HYTextBox*)components(5))->GetText();

            done = partName->IsValidIdentifier ();
            if (done) {
                f = dataSetFilterNamesList.Find (partName);
                if ((f>=0)&&(f!=dfID)) {
                    done = false;
                }
            }

            ((_HYButton*)components(0))->EnableButton (done);
        }

        done = true;
    } else if (e->EventClass()==_hyContextPopUp) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);

        if (i==6) {
            g = e->EventCode().Find(',',f+1,-1);
            f = e->EventCode().Cut(f+1,g-1).toNum();
            g = e->EventCode().Cut(g+1,-1).toNum();
            if (partitionDialogColorBox.Contains(f,g)) {
                _String  prompt ("New Partition Color");
                _HYColor newColor = SelectAColor (*color,prompt);
                if (!(newColor==*color)) {
                    *color = newColor;
                    _HYCanvas * cn2 = (_HYCanvas*)components (6);
                    cn2->StartDraw();
                    cn2->SetColor (*color);
                    cn2->FillRect (partitionDialogColorBox);
                    cn2->SetColor ((_HYColor) {
                        0,0,0
                    });
                    cn2->DrawRect (partitionDialogColorBox);
                    cn2->EndDraw();
                    cn2->_MarkForUpdate();
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

//__________________________________________________________________
_HYTextDialog::_HYTextDialog     (_String* res, _String& prompt, bool* resB, Ptr parent, _hyStringValidatorType textValidator):_HYTWindow ("Text Dialog",false,true, parent)
{
    result               = resB;
    textOut              = res;
    validator            = textValidator;
    lastValidationState  = textValidator?(*textValidator) (res):true;

    _HYRect         canvasSettings = {30,270,30,270,HY_COMPONENT_TRANSP_BG};
    _HYLabel*       l1 = (_HYLabel*)checkPointer(new _HYLabel (canvasSettings, GetOSWindowData())); // label

    canvasSettings.top   = canvasSettings.bottom  = 170;
    _HYTextBox*     tb   = (_HYTextBox*)checkPointer(new _HYTextBox (canvasSettings, GetOSWindowData()));
    tb->boxFlags |= HY_TB_WRAP;

    canvasSettings.left   = canvasSettings.right  = 15;
    canvasSettings.top    = 200;
    canvasSettings.bottom = 10000;

    _HYLabel*       l2 = (_HYLabel*)checkPointer(new _HYLabel (canvasSettings, GetOSWindowData())); // left spacer
    _HYLabel*       l3 = (_HYLabel*)checkPointer(new _HYLabel (canvasSettings, GetOSWindowData())); // right spacer

    canvasSettings.bottom = canvasSettings.top   = 40;
    canvasSettings.left   = canvasSettings.right = 210;

    _HYButton*      b1      = (_HYButton*)checkPointer(new _HYButton (canvasSettings, GetOSWindowData()));
    canvasSettings.left   = canvasSettings.right = 90;
    _HYButton*      b2      = (_HYButton*)checkPointer(new _HYButton (canvasSettings, GetOSWindowData()));

    tb->SetMessageRecipient (this);
    b1->SetMessageRecipient (this);
    b2->SetMessageRecipient (this);

    SetTableDimensions (3,3);

    AddObject (b1, false, 2, 0); // 0
    AddObject (b2, false, 2, 2); // 1
    AddObject (l1, false, 0, 1); // 2
    AddObject (l2, false, 0, 0); // 3
    AddObject (l3, false, 0, 2); // 4
    AddObject (tb, false, 1, 1); // 5


    SetCell (1,0,l2);
    SetCell (1,2,l3);
    SetCell (2,1,b1);

    _HYFont         labelFont;
    SetDefaultDialogFont(labelFont);


    b1->SetFont (labelFont);
    b2->SetFont (labelFont);
    l1->SetFont (labelFont);
    l2->SetFont (labelFont);
    l3->SetFont (labelFont);

    tb->SetText (*res);
    tb->SetFont (labelFont);

    b1->SetText (" OK ");
    b2->SetText (" Cancel ");
    b1->EnableButton (lastValidationState);

    b1->SetAlignFlags (HY_ALIGN_RIGHT);
    l1->SetAlignFlags (HY_ALIGN_LEFT);
    tb->SetAlignFlags (HY_ALIGN_LEFT);
    l1->SetText (prompt);

    b1->SetButtonKind (HY_BUTTON_OK);
    b2->SetButtonKind (HY_BUTTON_CANCEL);

    keyboardFocusChain << 5;
    ProcessEvent(generateKeyboardFocusEvent (tb->GetID()));

    _HYRect     dim = MinMaxWindowDimensions();

    SetWindowRectangle  (0,0,dim.bottom,dim.right);
    CenterWindow        (this);
}

//__________________________________________________________

bool    _HYTextDialog::ProcessEvent (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    long        i,f,k;

    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        i = MatchComponentID (firstArg);

        if (i==0) { // OK
            *result = true;
            *textOut = ((_HYTextBox*)GetObject (5))->GetText();
            postWindowCloseEvent (GetID());
        } else if (i==1) {
            *result = false;
            postWindowCloseEvent (GetID());
        }

        done = true;
    } else {
        if (validator) {
            if (e->EventClass() == _hyTextEditChange) {
                firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
                i = MatchComponentID (firstArg);

                if (i == 5) {
                    f =  e->EventCode().Cut(f+1,-1).toNum();
                    _HYTextBox*  myTB = (_HYTextBox*)components(i);
                    _String tt = myTB->GetText();
                    bool         good = (*validator) (&tt);

                    if (good != lastValidationState) {
                        ((_HYButton*)components(0))->EnableButton (good);
                        myTB->SetForeColor ((_HYColor) {
                            good?0:127,good?127:0,0
                        });
                        lastValidationState = good;
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

//__________________________________________________________________
_HYTextDialog2::_HYTextDialog2   (_String* res, _String* res2, _String& prompt, _String& prompt2, bool* resB, Ptr parent):_HYTextDialog (res,prompt,resB,parent)
{
    textOut2     = res2;

    _HYRect         canvasSettings = {30,270,30,270,HY_COMPONENT_TRANSP_BG};

    _HYLabel*       l4 = new _HYLabel (canvasSettings, GetOSWindowData());
    checkPointer    (l4);

    canvasSettings.top   = canvasSettings.bottom  = 170;
    _HYTextBox*     tb2 = new _HYTextBox (canvasSettings, GetOSWindowData());
    checkPointer    (tb2);
    tb2->boxFlags |= HY_TB_WRAP;


    tb2->SetMessageRecipient (this);

    SetTableDimensions (5,3);

    AddObject (l4); // 6
    AddObject (tb2);// 7

    SetCell (0,0,(_HYGuiObject*)components(3));
    SetCell (0,1,(_HYGuiObject*)components(2));
    SetCell (0,2,(_HYGuiObject*)components(4));

    SetCell (1,0,(_HYGuiObject*)components(3));
    SetCell (1,1,(_HYGuiObject*)components(5));
    SetCell (1,2,(_HYGuiObject*)components(4));

    SetCell (2,0,(_HYGuiObject*)components(3));
    SetCell (2,1,(_HYGuiObject*)l4);
    SetCell (2,2,(_HYGuiObject*)components(4));

    SetCell (3,0,(_HYGuiObject*)components(3));
    SetCell (3,1,(_HYGuiObject*)tb2);
    SetCell (3,2,(_HYGuiObject*)components(4));

    SetCell (4,0,(_HYGuiObject*)components(0));
    SetCell (4,1,(_HYGuiObject*)components(0));
    SetCell (4,2,(_HYGuiObject*)components(1));

    _HYFont         labelFont;
    SetDefaultDialogFont(labelFont);

    l4->SetFont (labelFont);
    l4->SetText (prompt2);

    tb2->SetText    (*res2);
    tb2->SetFont    (labelFont);

    l4->SetAlignFlags (HY_ALIGN_LEFT);
    tb2->SetAlignFlags (HY_ALIGN_LEFT);

    keyboardFocusChain << 7;

    _HYRect     dim = MinMaxWindowDimensions();

    SetWindowRectangle  (0,0,dim.bottom,dim.right);
    CenterWindow        (this);

    DeleteObject (l4);
    DeleteObject (tb2);

}

//__________________________________________________________

bool    _HYTextDialog2::ProcessEvent (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    long        i,f,k;

    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        i = MatchComponentID (firstArg);

        if (i==0) { // OK
            *result = true;
            *textOut  = ((_HYTextBox*)GetObject (5))->GetText();
            *textOut2 = ((_HYTextBox*)GetObject (7))->GetText();
            postWindowCloseEvent (GetID());
        }
    }

    if (done) {
        DeleteObject (e);
        return true;
    }
    return _HYTextDialog::ProcessEvent(e);
}


//__________________________________________________________________
_HYTextDialogWithPulldown::_HYTextDialogWithPulldown (_String* res,_String& prompt, _String& cPrompt, _List& menuOptions, _List* actions, bool *resB, bool* resC, long* resP, long initSel, Ptr parent):
    _HYTextDialogWithCheckbox    (res, prompt, cPrompt,resB, resC, parent)
{
    autoFill        = actions;
    menuSelection   = resP;

    _HYRect         canvasSettings = {30,270,30,270,HY_COMPONENT_TRANSP_BG};

    _HYPullDown*    pd = new _HYPullDown (canvasSettings, GetOSWindowData ());
    checkPointer    (pd);

    _HYGuiObject*   cache[7] = {GetObject(0), GetObject(1), GetObject (2),
                                GetObject(3), GetObject(4), GetObject (5), GetObject (6)
                               };

    SetTableDimensions (5,3);

    AddObject (pd);// 7

    SetCell (0,0,cache[3]);
    SetCell (0,1,cache[2]);
    SetCell (0,2,cache[4]);

    SetCell (1,0,cache[3]);
    SetCell (1,1,cache[5]);
    SetCell (1,2,cache[4]);

    SetCell (2,0,cache[3]);
    SetCell (2,1,pd);
    SetCell (2,2,cache[4]);

    SetCell (3,0,cache[3]);
    SetCell (3,1,cache[6]);
    SetCell (3,2,cache[4]);

    SetCell (4,0,cache[0]);
    SetCell (4,1,cache[0]);
    SetCell (4,2,cache[1]);

    ((_HYComponent*)(cache[3]))->settings.top+=30;
    ((_HYComponent*)(cache[3]))->settings.bottom+=30;

    ((_HYComponent*)(cache[4]))->settings.top+=30;
    ((_HYComponent*)(cache[4]))->settings.bottom+=30;

    pd->SetMessageRecipient (this);

    _HYFont         labelFont;

    for (long k=0; k<menuOptions.lLength; k++) {
        pd->AddMenuItem (*(_String*)menuOptions(k), -1);
    }

    pd->SetAlignFlags (HY_ALIGN_LEFT);
    pd->ChangeSelection (initSel, false);


    _HYRect     dim = MinMaxWindowDimensions();

    SetWindowRectangle  (0,0,dim.bottom,dim.right);
    CenterWindow        (this);

    DeleteObject (pd);
}

//__________________________________________________________

bool    _HYTextDialogWithPulldown::ProcessEvent (_HYEvent* e)
{
    _String     firstArg;
    long        i,f,k;

    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);

        if (i==0) { // OK
            *checkState     = ((_HYCheckbox*)GetObject(6))->GetState();
            *menuSelection  = ((_HYPullDown*)GetObject(7))->GetSelection();
        }
    } else {
        if (autoFill&&(e->EventClass()==_hyMenuSelChangeEvent)) {
            firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
            i = MatchComponentID (firstArg);


            if (i==7) { // OK
                k = e->EventCode().Cut(f+1,-1).toNum();
                if (autoFill->lLength>k) {
                    ((_HYTextBox*)GetObject(5))->InsertText (*(_String*)(*autoFill)(k),true);
                }
            }
        }
    }

    return _HYTextDialog::ProcessEvent(e);
}

//__________________________________________________________________
_HYTextDialogWithCheckbox::_HYTextDialogWithCheckbox     (_String* res, _String& prompt, _String& cPrompt, bool* resB, bool *resC, Ptr parent):
    _HYTextDialog                (res,prompt,resB, parent)
{
    checkState      = resC;

    _HYRect         canvasSettings = {30,270,30,270,HY_COMPONENT_TRANSP_BG};

    _HYCheckbox*    cb = new _HYCheckbox (canvasSettings, GetOSWindowData ());
    checkPointer    (cb);

    _HYGuiObject*   cache[6] = {GetObject(0), GetObject(1), GetObject (2),
                                GetObject(3), GetObject(4), GetObject (5)
                               };

    SetTableDimensions (4,3);

    AddObject (cb);// 6

    SetCell (0,0,cache[3]);
    SetCell (0,1,cache[2]);
    SetCell (0,2,cache[4]);

    SetCell (1,0,cache[3]);
    SetCell (1,1,cache[5]);
    SetCell (1,2,cache[4]);

    SetCell (2,0,cache[3]);
    SetCell (2,1,cb);
    SetCell (2,2,cache[4]);

    SetCell (3,0,cache[0]);
    SetCell (3,1,cache[0]);
    SetCell (3,2,cache[1]);

    ((_HYComponent*)(cache[3]))->settings.top+=30;
    ((_HYComponent*)(cache[3]))->settings.bottom+=30;

    ((_HYComponent*)(cache[4]))->settings.top+=30;
    ((_HYComponent*)(cache[4]))->settings.bottom+=30;

    cb->SetMessageRecipient (this);

    _HYFont         labelFont;
    SetDefaultDialogFont(labelFont);


    cb->SetFont (labelFont);
    cb->SetAlignFlags (HY_ALIGN_LEFT);
    cb->SetText (cPrompt);


    _HYRect     dim = MinMaxWindowDimensions();

    SetWindowRectangle  (0,0,dim.bottom,dim.right);
    CenterWindow        (this);

    cb->SetState (*resC,false);

    DeleteObject (cb);

}

//__________________________________________________________

bool    _HYTextDialogWithCheckbox::ProcessEvent (_HYEvent* e)
{
    _String     firstArg;
    long        i,f,k;

    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        i = MatchComponentID (firstArg);

        if (i==0) { // OK
            *checkState = ((_HYCheckbox*)GetObject(6))->GetState();
        }
    }

    return _HYTextDialog::ProcessEvent(e);
}

//__________________________________________________________________

_HYProceedPromptBox::_HYProceedPromptBox     (_String& prompt, _String& sb1, _String& sb2,_HYFont& tf, long pid, bool* resB, Ptr parent):
    _HYTWindow ("HYPHY Prompt",false,true, parent)
{
    result      = resB;

    _HYRect         canvasSettings = {150,400,150,400,HY_COMPONENT_TRANSP_BG};

    _HYCanvas*      cvs = new _HYCanvas (canvasSettings, GetOSWindowData(),150,400,24);
    checkPointer    (cvs);


    _HYFont         labelFont;
    labelFont.face  = "SystemFont";
    labelFont.size  = 14;
    labelFont.style = HY_FONT_PLAIN;

    canvasSettings.left = GetVisibleStringWidth (sb2, labelFont)+35;

    if (canvasSettings.left>325) {
        canvasSettings.left = 325;
    }

    canvasSettings.right  = canvasSettings.left;
    canvasSettings.bottom = canvasSettings.top   = 40;

    _HYButton*      b2      = nil;

    if (sb2.sLength) {
        b2 = new _HYButton (canvasSettings, GetOSWindowData());
        checkPointer (b2);
        b2->SetMessageRecipient (this);
    } else {
        canvasSettings.left = 20;
    }

    canvasSettings.left     = canvasSettings.right = 400 - canvasSettings.left;
    _HYButton*      b1      = new _HYButton (canvasSettings, GetOSWindowData());

    checkPointer (b1);

    b1->SetMessageRecipient (this);

    SetTableDimensions (2,2);

    AddObject (cvs);// 0
    AddObject (b1); // 1

    SetCell (0,0,cvs);
    SetCell (0,1,cvs);

    SetCell (1,0,b1);
    if (b2) {
        AddObject (b2); // 2
        SetCell (1,1,b2);
    } else {
        canvasSettings.left = canvasSettings.right = 20;
        _HYLabel* l1 = new _HYLabel (canvasSettings, GetOSWindowData());
        AddObject (l1);
        SetCell (1,1,l1);
        DeleteObject (l1);
    }


    cvs->StartDraw   ();
    cvs->SetDialogBG ();
    cvs->EraseAll    ();

    canvasSettings.left  = canvasSettings.right  = 10;
    canvasSettings.top   = canvasSettings.bottom = 10;

    cvs->DrawPicRes  (canvasSettings, pid);

    canvasSettings.left     = canvasSettings.right + 10;
    canvasSettings.right    = 390;
    canvasSettings.top      = 15;
    canvasSettings.bottom   = 145;

    cvs->SetFont     (tf);
    cvs->DisplayText (prompt, canvasSettings, HY_ALIGN_LEFT);

    cvs->EndDraw     ();


    b1->SetFont (labelFont);
    b1->SetText (sb1);
    b1->SetAlignFlags (HY_ALIGN_RIGHT);
    b1->SetButtonKind (HY_BUTTON_OK);

    if (b2) {
        b2->SetText (sb2);
        b2->SetFont (labelFont);
        b2->SetButtonKind (HY_BUTTON_CANCEL);
        DeleteObject (b2);
    }

    _HYRect     dim = MinMaxWindowDimensions();

    SetWindowRectangle  (0,0,dim.bottom,dim.right);
    CenterWindow        (this);

    DeleteObject (b1);
    DeleteObject (cvs);
}

//__________________________________________________________

bool    _HYProceedPromptBox::ProcessEvent (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    long        i,f,k;

    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        i = MatchComponentID (firstArg);

        if (i==1) { // OK
            *result = true;
            postWindowCloseEvent (GetID());
        } else if (i==2) {
            *result = false;
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

//__________________________________________________________________

_HYProceedPromptBoxWCheck::_HYProceedPromptBoxWCheck     (_String& prompt, _String& sb1, _String& sb2, _String& checkPrompt,_HYFont& tf, long pid, bool* resB, bool* checkB, Ptr parent):
    _HYProceedPromptBox (prompt, sb1, sb2, tf, pid, resB, parent)
{
    checkState  = checkB;

    _HYRect         canvasSettings = {tf.size + 5,400,tf.size + 5,400,HY_COMPONENT_TRANSP_BG};

    _HYFont         labelFont;
    SetDefaultDialogFont(labelFont);

    long            w = GetVisibleStringWidth (checkPrompt, tf) + 30;

    if ((tf.size < 36)&&(w+ 120 + GetVisibleStringWidth (sb1, labelFont) + GetVisibleStringWidth (sb2, labelFont)<=400)) {
        canvasSettings.left = canvasSettings.right  = w;
        canvasSettings.top  = canvasSettings.bottom = 40;
    }


    _HYCheckbox*    cb = new _HYCheckbox (canvasSettings, GetOSWindowData());
    checkPointer    (cb);


    AddObject (cb); // 3
    cb->SetMessageRecipient (this);

    _HYGuiObject * refs[3] = {(_HYGuiObject*)components(0),(_HYGuiObject*)components(1),(_HYGuiObject*)components(2)};

    if (canvasSettings.left==400) {
        SetTableDimensions (3,2);

        SetCell (0,0,refs[0]);
        SetCell (0,1,refs[0]);

        SetCell (1,0,cb);
        SetCell (1,1,cb);

        SetCell (2,0,refs[1]);
        SetCell (2,1,refs[2]);
    } else {
        SetTableDimensions (2,3);

        SetCell (0,0,refs[0]);
        SetCell (0,1,refs[0]);
        SetCell (0,2,refs[0]);

        SetCell (1,0,cb);
        SetCell (1,1,refs[1]);
        SetCell (1,2,refs[2]);

        _HYButton*      b = (_HYButton*)GetObject (1);
        b->settings.left = b->settings.right = b->settings.left - canvasSettings.left;
    }


    cb->SetFont       (tf);
    cb->SetState      (*checkState);
    cb->SetText       (checkPrompt);
    cb->SetAlignFlags (HY_ALIGN_LEFT);

    _HYRect     dim = MinMaxWindowDimensions();

    SetWindowRectangle  (0,0,dim.bottom,dim.right);
    CenterWindow        (this);

    DeleteObject (cb);
}

//__________________________________________________________

bool    _HYProceedPromptBoxWCheck::ProcessEvent (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    long        i,f,k;

    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        k = firstArg.toNum();
        i = MatchComponentID (firstArg);

        if (i==1) { // OK
            *result     = true;
            *checkState = ((_HYCheckbox*)GetObject (3))->GetState();
            postWindowCloseEvent (GetID());
        } else if (i==2) {
            *result = false;
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

//__________________________________________________________________

_HYInferenceConstraints::_HYInferenceConstraints
(_String& prompt ,  _List& l1, _List& l2, _List& l3, _List& l4,
 _List& L1, _List& L2, _List& L3, _List& L4,
 _HYFont& , bool* resB, _String* ress, Ptr parent):
    _HYTextDialog (ress, prompt, resB, parent)
{
    _HYRect         canvasSettings = {30,67,30,67,HY_COMPONENT_TRANSP_BG};

    _HYFont         labelFont;
    SetDefaultDialogFont(labelFont);

    il1 = &L1;
    il2 = &L2;
    il3 = &L3;
    il4 = &L4;

    _HYLabel       *lb1 = new _HYLabel (canvasSettings, GetOSWindowData ());
    checkPointer   (lb1);
    _HYLabel       *lb2 = new _HYLabel (canvasSettings, GetOSWindowData ());
    checkPointer   (lb2);
    _HYLabel       *lb3 = new _HYLabel (canvasSettings, GetOSWindowData ());
    checkPointer   (lb3);
    _HYLabel       *lb4 = new _HYLabel (canvasSettings, GetOSWindowData ());
    checkPointer   (lb4);

    canvasSettings.left = canvasSettings.right = 178;

    _HYPullDown    *pd1 = new _HYPullDown (canvasSettings, GetOSWindowData  ());
    checkPointer   (pd1);
    _HYPullDown    *pd2 = new _HYPullDown (canvasSettings, GetOSWindowData  ());
    checkPointer   (pd2);
    _HYPullDown    *pd3 = new _HYPullDown (canvasSettings, GetOSWindowData  ());
    checkPointer   (pd3);
    _HYPullDown    *pd4 = new _HYPullDown (canvasSettings, GetOSWindowData  ());
    checkPointer   (pd4);

    canvasSettings.left = canvasSettings.right = 25;

    _HYButtonBar   *btn1 = new _HYButtonBar (canvasSettings, GetOSWindowData    ());
    checkPointer   (btn1);
    _HYButtonBar   *btn2 = new _HYButtonBar (canvasSettings, GetOSWindowData    ());
    checkPointer   (btn2);
    _HYButtonBar   *btn3 = new _HYButtonBar (canvasSettings, GetOSWindowData    ());
    checkPointer   (btn3);
    _HYButtonBar   *btn4 = new _HYButtonBar (canvasSettings, GetOSWindowData    ());
    checkPointer   (btn4);

    AddObject      (lb1); // 6
    AddObject      (lb2); // 7
    AddObject      (lb3); // 8
    AddObject      (lb4); // 9

    AddObject      (pd1); // 10
    AddObject      (pd2); // 11
    AddObject      (pd3); // 12
    AddObject      (pd4); // 13

    AddObject      (btn1);// 14
    AddObject      (btn2);// 15
    AddObject      (btn3);// 16
    AddObject      (btn4);// 17

    _HYGuiObject * refs[6] = {(_HYGuiObject*)components(0),(_HYGuiObject*)components(1),(_HYGuiObject*)components(2),
                              (_HYGuiObject*)components(3),(_HYGuiObject*)components(4),(_HYGuiObject*)components(5)
                             };
    SetTableDimensions (7,7);

    SetCell (0,0,refs[3]);
    SetCell (0,1,refs[2]);
    SetCell (0,2,refs[2]);
    SetCell (0,3,refs[2]);
    SetCell (0,4,refs[2]);
    SetCell (0,5,refs[2]);
    SetCell (0,6,refs[4]);

    SetCell (1,0,refs[3]);
    SetCell (1,1,refs[5]);
    SetCell (1,2,refs[5]);
    SetCell (1,3,refs[5]);
    SetCell (1,4,refs[5]);
    SetCell (1,5,refs[5]);
    SetCell (1,6,refs[4]);

    SetCell (2,0,refs[3]);
    SetCell (2,1,lb1);
    SetCell (2,2,pd1);
    SetCell (2,3,pd1);
    SetCell (2,4,pd1);
    SetCell (2,5,btn1);
    SetCell (2,6,refs[4]);

    SetCell (3,0,refs[3]);
    SetCell (3,1,lb2);
    SetCell (3,2,pd2);
    SetCell (3,3,pd2);
    SetCell (3,4,pd2);
    SetCell (3,5,btn2);
    SetCell (3,6,refs[4]);

    SetCell (4,0,refs[3]);
    SetCell (4,1,lb3);
    SetCell (4,2,pd3);
    SetCell (4,3,pd3);
    SetCell (4,4,pd3);
    SetCell (4,5,btn3);
    SetCell (4,6,refs[4]);

    SetCell (5,0,refs[3]);
    SetCell (5,1,lb4);
    SetCell (5,2,pd4);
    SetCell (5,3,pd4);
    SetCell (5,4,pd4);
    SetCell (5,5,btn4);
    SetCell (5,6,refs[4]);



    SetCell (6,0,refs[3]);
    SetCell (6,1,refs[0]);
    SetCell (6,2,refs[0]);
    SetCell (6,3,refs[1]);
    SetCell (6,4,refs[1]);
    SetCell (6,5,refs[1]);
    SetCell (6,6,refs[4]);

    ((_HYComponent*)refs[0])->settings.left -= 30;
    ((_HYComponent*)refs[0])->settings.right -= 30;

    pd1->SetAlignFlags (HY_ALIGN_LEFT);
    pd2->SetAlignFlags (HY_ALIGN_LEFT);
    pd3->SetAlignFlags (HY_ALIGN_LEFT);
    pd4->SetAlignFlags (HY_ALIGN_LEFT);

    lb1->SetAlignFlags (HY_ALIGN_LEFT);
    lb2->SetAlignFlags (HY_ALIGN_LEFT);
    lb3->SetAlignFlags (HY_ALIGN_LEFT);
    lb4->SetAlignFlags (HY_ALIGN_LEFT);

    lb1->SetFont       (labelFont);
    lb2->SetFont       (labelFont);
    lb3->SetFont       (labelFont);
    lb4->SetFont       (labelFont);

    lb1->SetText ("Locals");
    lb2->SetText ("Globals");
    lb3->SetText ("Trees");
    lb4->SetText ("Templates");

    _String         btnTT ("Paste to expression");

    btn1->SetMessageRecipient (this);
    btn2->SetMessageRecipient (this);
    btn3->SetMessageRecipient (this);
    btn4->SetMessageRecipient (this);

    btn1->SetButtonDim  (16);
    btn2->SetButtonDim  (16);
    btn3->SetButtonDim  (16);
    btn4->SetButtonDim  (16);

    btn1->AddButton (ProcureIconResource(HY_MODELWINDOW_ICON_BASE+2),&btnTT);
    btn2->AddButton (ProcureIconResource(HY_MODELWINDOW_ICON_BASE+2),&btnTT);
    btn3->AddButton (ProcureIconResource(HY_MODELWINDOW_ICON_BASE+2),&btnTT);
    btn4->AddButton (ProcureIconResource(HY_MODELWINDOW_ICON_BASE+2),&btnTT);

    long k;

    for (k=0; k<l1.lLength; k++) {
        pd1->AddMenuItem (*(_String*)l1(k),-1);
    }

    for (k=0; k<l2.lLength; k++) {
        pd2->AddMenuItem (*(_String*)l2(k),-1);
    }

    for (k=0; k<l3.lLength; k++) {
        pd3->AddMenuItem (*(_String*)l3(k),-1);
    }

    for (k=0; k<l4.lLength; k++) {
        pd4->AddMenuItem (*(_String*)l4(k),-1);
    }

    if (l1.lLength == 0) {
        pd1->AddMenuItem (none,0);
        pd1->EnableMenu  (false);
        btn1->EnableButton (0,false);
    }

    if (l2.lLength == 0) {
        pd2->AddMenuItem (none,0);
        pd2->EnableMenu  (false);
        btn2->EnableButton (0,false);
    }

    if (l3.lLength == 0) {
        pd3->AddMenuItem (none,0);
        pd3->EnableMenu  (false);
        btn3->EnableButton (0,false);
    }

    if (l4.lLength == 0) {
        pd4->AddMenuItem (none,0);
        pd4->EnableMenu  (false);
        btn4->EnableButton (0,false);
    }

    pd1->ChangeSelection (0,false);
    pd2->ChangeSelection (0,false);
    pd3->ChangeSelection (0,false);
    pd4->ChangeSelection (0,false);

    _HYRect     dim = MinMaxWindowDimensions();

    SetWindowRectangle  (0,0,dim.bottom,dim.right);
    CenterWindow        (this);

    DeleteObject (pd1);
    DeleteObject (pd2);
    DeleteObject (pd3);
    DeleteObject (pd4);

    DeleteObject (btn1);
    DeleteObject (btn2);
    DeleteObject (btn3);
    DeleteObject (btn4);

    DeleteObject (lb1);
    DeleteObject (lb2);
    DeleteObject (lb3);
    DeleteObject (lb4);
}

//__________________________________________________________

bool    _HYInferenceConstraints::ProcessEvent (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    long        i,f;

    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);

        if (i==0) { // OK
            *result = true;
            *textOut = ((_HYTextBox*)GetObject (5))->GetText();
            postWindowCloseEvent (GetID());
        } else if (i==1) {
            *result = false;
            postWindowCloseEvent (GetID());
        } else {
            _List           * theList = nil;
            _HYPullDown     * thePD   = nil;
            _HYTextBox      * tb = (_HYTextBox*)GetObject (5);
            switch (i) {
            case 14:
                theList = il1;
                thePD   = (_HYPullDown*)GetObject (10);
                break;
            case 15:
                theList = il2;
                thePD   = (_HYPullDown*)GetObject (11);
                break;
            case 16:
                theList = il3;
                thePD   = (_HYPullDown*)GetObject (12);
                break;
            case 17:
                theList = il4;
                thePD   = (_HYPullDown*)GetObject (13);
                break;
            }

            if (theList) {
                tb->InsertText (*(_String*)((*theList)(thePD->GetSelection())));
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

char        ciOption    = 0;
_Parameter  ciConfLevel = 0.95,
            ciResampler = 100;

//__________________________________________________________

_HYCIDialog::_HYCIDialog (bool* res, char* opt, _Parameter* val, _Parameter * val2, Ptr parWindow):_HYTWindow ("Confidence Interval Options", false, true, parWindow)
{


    result          = res;
    option          = opt;
    sigValue        = val;
    sigValue2       = val2;

    _HYRect         canvasSettings = {30,150,30,150,HY_COMPONENT_NO_SCROLL|HY_COMPONENT_TRANSP_BG};

    _HYLabel*       l1      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l2      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l3      = new _HYLabel (canvasSettings, GetOSWindowData());


    canvasSettings.left = canvasSettings.right = 211;

    _HYTextBox  *   tb1     = new _HYTextBox  (canvasSettings, GetOSWindowData());
    _HYTextBox  *   tb2     = new _HYTextBox  (canvasSettings, GetOSWindowData());
    _HYPullDown *   p1      = new _HYPullDown (canvasSettings, GetOSWindowData());
    canvasSettings.left = canvasSettings.right = 280;
    _HYButton*      b1      = new _HYButton (canvasSettings, GetOSWindowData());
    canvasSettings.left = canvasSettings.right = 81;
    _HYButton*      b2      = new _HYButton (canvasSettings, GetOSWindowData());

    tb1->SetMessageRecipient (this);
    tb2->SetMessageRecipient (this);
    b1 ->SetMessageRecipient (this);
    b2 ->SetMessageRecipient (this);
    p1 ->SetMessageRecipient (this);

    AddObject (tb1); // 0

    AddObject (b1);  // 1
    AddObject (b2);  // 2

    AddObject (l1);  // 3
    AddObject (l2);  // 4

    AddObject (p1);  // 5

    AddObject (l3);   // 6
    AddObject (tb2);  // 7


    keyboardFocusChain << 0;
    keyboardFocusChain << 7;

    SetTableDimensions (4,2);

    SetCell   (0,0,l1);
    SetCell   (0,1,p1);

    SetCell   (1,0,l2);
    SetCell   (1,1,tb1);

    SetCell   (2,0,l3);
    SetCell   (2,1,tb2);

    SetCell   (3,0,b1);
    SetCell   (3,1,b2);


    _HYFont  labelFont;
    SetDefaultDialogFont(labelFont);

    l1->SetFont (labelFont);
    l2->SetFont (labelFont);
    l3->SetFont (labelFont);

    b1->SetFont (labelFont);
    b2->SetFont (labelFont);

    l1->SetText (" Estimation Method");

    tb1->SetText(ciConfLevel);
    tb1->SetFont (labelFont);

    tb2->SetText(ciResampler);
    tb2->SetFont (labelFont);

    l1 ->SetAlignFlags (HY_ALIGN_LEFT);
    l2 ->SetAlignFlags (HY_ALIGN_LEFT);
    l3 ->SetAlignFlags (HY_ALIGN_LEFT);
    p1 ->SetAlignFlags (HY_ALIGN_LEFT);
    tb1->SetAlignFlags (HY_ALIGN_LEFT);
    tb2->SetAlignFlags (HY_ALIGN_LEFT);

    b1->SetAlignFlags (HY_ALIGN_RIGHT);

    b1->SetText       ("  OK  ");
    b1->SetButtonKind (HY_BUTTON_OK);
    b2->SetText       (" Cancel ");
    b2->SetButtonKind (HY_BUTTON_CANCEL);

    p1->AddMenuItem   ("Asymptotic Normal [crude]",-1);
    p1->AddMenuItem   ("Asymptotic Normal [finer]",-1);
    p1->AddMenuItem   (menuSeparator,-1);
    p1->AddMenuItem   ("Likelihood Profile [chi2]",-1);
    p1->AddMenuItem   ("Likelihood Profile [custom]",-1);
    p1->AddMenuItem   (menuSeparator,-1);
    p1->AddMenuItem   ("SIR Sampler",-1);
    p1->AddMenuItem   ("Latin Hypercube Sampler",-1);


    _HYRect     dim = MinMaxWindowDimensions();

    lastState = 2;
    *option   = -1;

    tb1->SetText        (ciConfLevel,false);
    tb2->SetText        (ciResampler,false);
    p1->ChangeSelection (ciOption,true);

    SetWindowRectangle  (0,0,dim.bottom,dim.right);
    CenterWindow        (this);

    DeleteObject (l1);
    DeleteObject (l2);
    DeleteObject (l3);
    DeleteObject (tb1);
    DeleteObject (tb2);
    DeleteObject (p1);
    DeleteObject (b1);
    DeleteObject (b2);
}

//__________________________________________________________

bool    _HYCIDialog::ProcessEvent (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    long        i,f;

    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);

        if (i==1) { // OK
            *option     = ciOption = ((_HYPullDown*)components (5))->GetSelection();
            *result     = true;
            ciConfLevel = *sigValue;
            if (*option > 5) {
                ciResampler = *sigValue2;
            }
            postWindowCloseEvent (GetID());
        } else if (i==2) { // Cancel
            *result = false;
            postWindowCloseEvent (GetID());
        }

        done = true;
    } else if (e->EventClass()==_hyTextEditChange) {
        *sigValue = (((_HYTextBox*)components (0))->GetText()).toNum();
        *sigValue2 = (((_HYTextBox*)components (7))->GetText()).toNum();

        bool         good = true;

        switch (((_HYPullDown*)components (5))->GetSelection()) {
        case 3:
            good = *sigValue>0.0 && *sigValue<1.0;
            break;
        case 4:
            good = *sigValue>0.0;
            break;
        case 6:
            good =  *sigValue>=1 && *sigValue2<=*sigValue;
            break;
        }

        if (good!=lastState) {
            lastState = good;
            ((_HYButton*)components(1))->EnableButton (good);
            ((_HYTextBox*)components (0))->SetForeColor ((_HYColor) {
                good?0:127,good?127:0,0
            });
            ((_HYTextBox*)components (7))->SetForeColor ((_HYColor) {
                good?0:127,good?127:0,0
            });
        }

        done = true;
    } else if (e->EventClass()==_hyMenuSelChangeEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        f        = e->EventCode().Cut (f+1,-1).toNum();
        i        = MatchComponentID (firstArg);

        if (i==5) {
            if (*option != f) {
                _HYLabel    *l2     = (_HYLabel*)components(4),
                             *l3     = (_HYLabel*)components(6);

                _HYTextBox  *tb1    = (_HYTextBox*)components(0),
                             *tb2    = (_HYTextBox*)components(7);

                if (f>=6) {
                    if (*option < 6)
                        l3->SetForeColor((_HYColor) {
                        0,0,0
                    });

                    l3->SetText(" Resample M");
                    tb2->SetText(tb2->GetText(), true);
                    lastState = -2;
                } else {
                    l3->SetText(" Not applicable");
                    l3->SetForeColor((_HYColor) {
                        127,127,127
                    });
                    l3->_MarkForUpdate();
                    tb2->SetForeColor   ((_HYColor) {
                        255,255,255
                    });
                }

                if (f<2) {
                    if ((*option > 2)||(*option < 0)) {
                        l2->SetText(" Not applicable");
                        l2->SetForeColor((_HYColor) {
                            127,127,127
                        });
                        l2->_MarkForUpdate();
                        tb1->SetForeColor   ((_HYColor) {
                            255,255,255
                        });
                        ((_HYButton*)components(1))->EnableButton (true);
                        lastState = -2;
                    }
                } else {
                    if (*option < 2)
                        l2->SetForeColor((_HYColor) {
                        0,0,0
                    });

                    if (f==3) {
                        l2->SetText(" Chi^2 significance (0-1)");
                    } else if (f>=6) {
                        l2->SetText(" Sample N");
                    } else {
                        l2->SetText(" Log(L) difference (>0)");
                    }

                    l2->_MarkForUpdate();
                    tb1->SetText(tb1->GetText(), true);
                }

                *option = f;

                tb1->EnableTextEdit (f>2);
                tb2->EnableTextEdit (f>5);
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

bool         lProfileNormPlot   = true;
long         lProfileIntervals  = 16;
_Parameter   lProfileLB         = 0.;
_Parameter   lProfileUB         = 1.;

//__________________________________________________________

_HYLProfDialog::_HYLProfDialog (bool*res,_Parameter*lb,_Parameter*ub,long*st,bool*np,_Parameter ml,_Parameter lbv,_Parameter ubv,Ptr parWindow):
    _HYTWindow ("Likelihood Profile Options", false, true, parWindow)
{
    result          = res;
    left            = lb;
    right           = ub;
    mle             = ml;
    intervals       = st;
    doNorm          = np;
    vlb             = lbv;
    vub             = ubv;


    _HYRect         canvasSettings = {30,230,30,230,HY_COMPONENT_NO_SCROLL|HY_COMPONENT_TRANSP_BG};

    _HYLabel*       l1      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l2      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l3      = new _HYLabel (canvasSettings, GetOSWindowData());


    canvasSettings.left = canvasSettings.right = 111;

    _HYTextBox  *   tb1     = new _HYTextBox  (canvasSettings, GetOSWindowData());
    _HYTextBox  *   tb2     = new _HYTextBox  (canvasSettings, GetOSWindowData());
    _HYTextBox  *   tb3     = new _HYTextBox  (canvasSettings, GetOSWindowData());

    canvasSettings.left = canvasSettings.right = 260;
    _HYButton*      b1      = new _HYButton (canvasSettings, GetOSWindowData());
    canvasSettings.left = canvasSettings.right = 81;
    _HYButton*      b2      = new _HYButton (canvasSettings, GetOSWindowData());

    canvasSettings.left = canvasSettings.right = 340;

    _HYCheckbox*    cb      = new _HYCheckbox (canvasSettings, GetOSWindowData());

    tb1->SetMessageRecipient (this);
    tb2->SetMessageRecipient (this);
    tb3->SetMessageRecipient (this);
    cb ->SetMessageRecipient (this);
    b1 ->SetMessageRecipient (this);
    b2 ->SetMessageRecipient (this);

    AddObject (tb1); // 0
    AddObject (tb2); // 1
    AddObject (tb3); // 2

    AddObject (b1);  // 3
    AddObject (b2);  // 4

    AddObject (l1);  // 5
    AddObject (l2);  // 6
    AddObject (l3);  // 7

    AddObject (cb);  // 8


    keyboardFocusChain << 0;
    keyboardFocusChain << 1;
    keyboardFocusChain << 2;

    SetTableDimensions (5,2);

    SetCell   (0,0,l1);
    SetCell   (0,1,tb1);

    SetCell   (1,0,l2);
    SetCell   (1,1,tb2);

    SetCell   (2,0,l3);
    SetCell   (2,1,tb3);

    SetCell   (3,0,cb);
    SetCell   (3,1,cb);

    SetCell   (4,0,b1);
    SetCell   (4,1,b2);


    _HYFont  labelFont;
    SetDefaultDialogFont(labelFont);

    l1->SetFont (labelFont);
    l2->SetFont (labelFont);
    l3->SetFont (labelFont);
    cb->SetFont (labelFont);

    b1->SetFont (labelFont);
    b2->SetFont (labelFont);

    l1->SetText (_String(" Lower Bound [") & vlb & ',' & mle & ']' );
    l2->SetText (_String(" Upper Bound [") & mle & ',' & vub & ']' );
    l3->SetText (" Plot Points");
    cb->SetText ("Show Quadratic Log-Likelihood Approximation");

    tb1->SetFont (labelFont);
    tb2->SetFont (labelFont);
    tb3->SetFont (labelFont);

    l1 ->SetAlignFlags (HY_ALIGN_LEFT);
    l2 ->SetAlignFlags (HY_ALIGN_LEFT);
    l3 ->SetAlignFlags (HY_ALIGN_LEFT);
    cb ->SetAlignFlags (HY_ALIGN_LEFT);

    tb1->SetAlignFlags (HY_ALIGN_LEFT);
    tb2->SetAlignFlags (HY_ALIGN_LEFT);
    tb3->SetAlignFlags (HY_ALIGN_LEFT);

    b1->SetAlignFlags (HY_ALIGN_RIGHT);
    b1->SetText       ("  OK  ");
    b1->SetButtonKind (HY_BUTTON_OK);
    b1->EnableButton  (false);
    b2->SetText       (" Cancel ");
    b2->SetButtonKind (HY_BUTTON_CANCEL);

    cb->SetState      (lProfileNormPlot,  false);
    tb3->SetText      (lProfileIntervals, false);
    tb1->SetText      (lProfileLB, true);
    tb2->SetText      (lProfileUB, true);

    lastState[0] = -1;
    lastState[1] = -1;
    lastState[2] = -1;

    _HYRect             dim = MinMaxWindowDimensions();

    SetWindowRectangle  (0,0,dim.bottom,dim.right);
    CenterWindow        (this);

    *intervals = lProfileIntervals;

    DeleteObject (l1);
    DeleteObject (l2);
    DeleteObject (l3);
    DeleteObject (tb1);
    DeleteObject (tb2);
    DeleteObject (tb3);
    DeleteObject (cb);
    DeleteObject (b1);
    DeleteObject (b2);
}

//__________________________________________________________

bool    _HYLProfDialog::ProcessEvent (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    long        i,f;

    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);

        if (i==3) { // OK
            *result          = true;
            lProfileNormPlot = *doNorm = ((_HYCheckbox*)components (8))->GetState();
            lProfileIntervals= *intervals;
            lProfileLB       = *left;
            lProfileUB       = *right;
            postWindowCloseEvent (GetID());
        } else if (i==4) { // Cancel
            *result = false;
            postWindowCloseEvent (GetID());
        }

        done = true;
    } else if (e->EventClass()==_hyTextEditChange) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);
        f =  e->EventCode().Cut(f+1,-1).toNum();

        if (f) {
            bool         good;
            _HYTextBox*  myTB = (_HYTextBox*)components(i);

            switch (i) {
            case 0:
                *left = myTB->GetText().toNum();
                if (*left>=vlb && *left<=mle && *left<*right) {
                    good = true;
                } else {
                    good = false;
                }
                break;
            case 1:
                *right = myTB->GetText().toNum();
                if (*right>=mle && *right<=vub && *left<*right) {
                    good = true;
                } else {
                    good = false;
                }
                break;
            case 2: {
                _Parameter temp = myTB->GetText().toNum();
                if (temp>2 && CheckEqual(ceil(temp),temp)) {
                    *intervals = temp;
                    good = true;
                } else {
                    good = false;
                }
                break;
            }
            }

            if (good!=lastState[i]) {
                ((_HYButton*)components(3))->EnableButton (good);
                myTB->SetForeColor ((_HYColor) {
                    good?0:127,good?127:0,0
                });
                lastState[i] = good;
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

long    defNewChartRows     = 4,
        defNewChartColumns   = 4;

//__________________________________________________________

_HYNewChart::_HYNewChart (long* rows, long* columns, bool* res, _String* chartTitle):_HYTWindow ("New Chart Window Setup", false, true)
{

    //long          index,
    //              cellWidth;

    result   = res;
    resR     = rows;
    resC     = columns;
    title    = chartTitle;

    _HYRect         canvasSettings = {30,150,30,150,HY_COMPONENT_NO_SCROLL|HY_COMPONENT_TRANSP_BG};

    _HYLabel*       l1      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l2      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l3      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l4      = new _HYLabel (canvasSettings, GetOSWindowData());


    canvasSettings.left = canvasSettings.right = 151;

    _HYTextBox  *   tb1     = new _HYTextBox  (canvasSettings, GetOSWindowData());
    _HYTextBox  *   tb2     = new _HYTextBox  (canvasSettings, GetOSWindowData());
    _HYTextBox  *   tb3     = new _HYTextBox  (canvasSettings, GetOSWindowData());
    _HYPullDown *   p1      = new _HYPullDown (canvasSettings, GetOSWindowData());
    canvasSettings.left = canvasSettings.right = 220;
    _HYButton*      b1      = new _HYButton (canvasSettings, GetOSWindowData());
    canvasSettings.left = canvasSettings.right = 81;
    _HYButton*      b2      = new _HYButton (canvasSettings, GetOSWindowData());

    tb1->SetMessageRecipient (this);
    tb2->SetMessageRecipient (this);
    tb3->SetMessageRecipient (this);
    b1 ->SetMessageRecipient (this);
    b2 ->SetMessageRecipient (this);
    p1 ->SetMessageRecipient (this);

    AddObject (tb1); // 0
    AddObject (tb2); // 1
    AddObject (tb3); // 2

    AddObject (b1);  // 3
    AddObject (b2);  // 4

    AddObject (l1);  // 5
    AddObject (l2);  // 6
    AddObject (l3);  // 7
    AddObject (l4);  // 8

    AddObject (p1);  // 9


    keyboardFocusChain << 0;
    keyboardFocusChain << 1;
    keyboardFocusChain << 2;

    SetTableDimensions (5,2);

    SetCell   (0,0,l1);
    SetCell   (0,1,tb1);

    SetCell   (1,0,l2);
    SetCell   (1,1,tb2);

    SetCell   (2,0,l3);
    SetCell   (2,1,tb3);

    SetCell   (3,0,l4);
    SetCell   (3,1,p1);

    SetCell   (4,0,b1);
    SetCell   (4,1,b2);


    _HYFont  labelFont;
    SetDefaultDialogFont(labelFont);

    l1->SetFont (labelFont);
    l2->SetFont (labelFont);
    l3->SetFont (labelFont);
    l4->SetFont (labelFont);

    b1->SetFont (labelFont);
    b2->SetFont (labelFont);

    l1->SetText (" Series (Columns)");
    l2->SetText (" Data Points (Rows)");
    l3->SetText (" Window Title");
    l4->SetText (" Spawn from Matrix");

    tb1->SetText(defNewChartColumns);
    tb2->SetText(defNewChartRows);
    tb3->SetText(*title);

    tb1->SetFont (labelFont);
    tb2->SetFont (labelFont);
    tb3->SetFont (labelFont);

    l1 ->SetAlignFlags (HY_ALIGN_LEFT);
    l2 ->SetAlignFlags (HY_ALIGN_LEFT);
    l3 ->SetAlignFlags (HY_ALIGN_LEFT);
    l4 ->SetAlignFlags (HY_ALIGN_LEFT);
    p1 ->SetAlignFlags (HY_ALIGN_LEFT);
    tb1->SetAlignFlags (HY_ALIGN_LEFT);
    tb2->SetAlignFlags (HY_ALIGN_LEFT);
    tb3->SetAlignFlags (HY_ALIGN_LEFT);

    b1->SetAlignFlags (HY_ALIGN_RIGHT);

    b1->SetText       ("  OK  ");
    b1->SetButtonKind (HY_BUTTON_OK);
    b2->SetText       (" Cancel ");
    b2->SetButtonKind (HY_BUTTON_CANCEL);

    p1->AddMenuItem   (none,-1);

    _List   mxNames;

#ifndef USE_AVL_NAMES
    for (long k=0; k<variableNames.lLength; k++)
#else
    _SimpleList tcache;
    long        vi,
                k = variableNames.Traverser (tcache,vi,variableNames.GetRoot());

    for (; k>=0; k = variableNames.Traverser (tcache,vi))
#endif
    {
        _Variable * thisV = FetchVar (k);
        if (thisV->ObjectClass () == MATRIX) {
            mxNames && thisV->GetName();
        }
    }

    if (mxNames.lLength) {
        p1->AddMenuItem (menuSeparator,-1);
        mxNames.Sort();
        for (long kk=0; kk<mxNames.lLength; kk++) {
            p1->AddMenuItem (*(_String*)mxNames(kk),-1);
        }
        p1->EnableMenu(true);
    } else {
        p1->EnableMenu(false);
    }


    ProcessEvent (generateKeyboardFocusEvent (tb1->GetID()));

    DeleteObject (l1);
    DeleteObject (l2);
    DeleteObject (l3);
    DeleteObject (l4);
    DeleteObject (tb1);
    DeleteObject (tb2);
    DeleteObject (tb3);
    DeleteObject (b1);
    DeleteObject (b2);
    DeleteObject (p1);

    _HYRect     dim = MinMaxWindowDimensions();

    SetWindowRectangle  (0,0,dim.bottom,dim.right);
    CenterWindow        (this);
}

//__________________________________________________________

bool    _HYNewChart::ProcessEvent (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    long        i,f;

    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);

        if (i==3) { // OK
            f = ((_HYPullDown*)components (9))->GetSelection();
            if (f) {
                *resR = -LocateVarByName (*(_String*)(((_HYPullDown*)components (9))->GetMenuItem(f)));
                if (*resR>0) {
                    *result = true;
                    postWindowCloseEvent (GetID());
                }
            } else {
                *resC = defNewChartColumns= ((_HYTextBox*)components (0))->GetText().toNum();
                *resR = defNewChartRows   = ((_HYTextBox*)components (1))->GetText().toNum();
            }
            *title= ((_HYTextBox*)components (2))->GetText();
            *result = true;
            postWindowCloseEvent (GetID());
        } else if (i==4) { // Cancel
            *result = false;
            postWindowCloseEvent (GetID());
        }

        done = true;
    } else if (e->EventClass()==_hyTextEditChange) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i        = MatchComponentID (firstArg);

        if ((i==0)||(i==1)||(i==2)) {
            _HYButton* okButton     = (_HYButton*)components  (3);
            _String    textField    = ((_HYTextBox*)components (i))->GetText();
            if (i<2) {
                _Parameter fieldValue   = textField.toNum();

                if ((fieldValue>0.0)&&(fieldValue == (long) (fieldValue))) {
                    done = true;
                }
            } else {
                done = (textField.FirstNonSpaceIndex (0,-1,1)!=-1);
            }

            okButton->EnableButton (done);
        }

        done = true;
    } else if (e->EventClass()==_hyMenuSelChangeEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        f        = e->EventCode().Cut (f+1,-1).toNum();
        i        = MatchComponentID (firstArg);

        if (i==9) {
            if (f==0) {
                ((_HYTextBox*)components (0))->EnableTextEdit (true);
                ((_HYTextBox*)components (1))->EnableTextEdit (true);
            } else {
                firstArg = *(_String*)((_HYPullDown*)components (9))->GetMenuItem(f);
                ((_HYTextBox*)components (2))->SetText(firstArg,true);

                _Matrix *   nt =  (_Matrix*)FetchVar(LocateVarByName(firstArg))->GetValue();

                firstArg = nt->GetVDim();
                ((_HYTextBox*)components (0))->EnableTextEdit (false);
                ((_HYTextBox*)components (0))->SetText        (firstArg,true);
                firstArg = nt->GetHDim();
                ((_HYTextBox*)components (1))->EnableTextEdit (false);
                ((_HYTextBox*)components (1))->SetText        (firstArg,true);
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

long    defTreeHPages       = 1,
        defTreeVPages       = 1;

//__________________________________________________________

_HYTreePrintPrefs::_HYTreePrintPrefs (long* columns, long* rows, bool* res, Ptr parent):_HYTWindow ("Tree Printing Prefs", false, true, parent)
{
    result   = res;
    resH     = rows;
    resW     = columns;

    _HYRect         canvasSettings = {30,150,30,150,HY_COMPONENT_NO_SCROLL|HY_COMPONENT_TRANSP_BG};

    _HYLabel*       l1      = new _HYLabel (canvasSettings, GetOSWindowData());
    _HYLabel*       l2      = new _HYLabel (canvasSettings, GetOSWindowData());

    canvasSettings.left = canvasSettings.right = 151;

    _HYTextBox  *   tb1     = new _HYTextBox  (canvasSettings, GetOSWindowData());
    _HYTextBox  *   tb2     = new _HYTextBox  (canvasSettings, GetOSWindowData());

    canvasSettings.left = canvasSettings.right = 220;
    _HYButton*      b1      = new _HYButton (canvasSettings, GetOSWindowData());
    canvasSettings.left = canvasSettings.right = 81;
    _HYButton*      b2      = new _HYButton (canvasSettings, GetOSWindowData());

    tb1->SetMessageRecipient (this);
    tb2->SetMessageRecipient (this);
    b1 ->SetMessageRecipient (this);
    b2 ->SetMessageRecipient (this);

    AddObject (tb1); // 0
    AddObject (tb2); // 1

    AddObject (b1);  // 2
    AddObject (b2);  // 3

    AddObject (l1);  // 4
    AddObject (l2);  // 5


    keyboardFocusChain << 0;
    keyboardFocusChain << 1;

    SetTableDimensions (3,2);

    SetCell   (0,0,l1);
    SetCell   (0,1,tb1);

    SetCell   (1,0,l2);
    SetCell   (1,1,tb2);

    SetCell   (2,0,b1);
    SetCell   (2,1,b2);


    _HYFont  labelFont;
    SetDefaultDialogFont(labelFont);

    l1->SetFont (labelFont);
    l2->SetFont (labelFont);

    b1->SetFont (labelFont);
    b2->SetFont (labelFont);

    l1->SetText (" Pages Wide (0 = Auto)");
    l2->SetText (" Pages Height (0 = Auto)");

    tb1->SetText(defTreeHPages);
    tb2->SetText(defTreeVPages);

    tb1->SetFont (labelFont);
    tb2->SetFont (labelFont);

    l1 ->SetAlignFlags (HY_ALIGN_LEFT);
    l2 ->SetAlignFlags (HY_ALIGN_LEFT);
    tb1->SetAlignFlags (HY_ALIGN_LEFT);
    tb2->SetAlignFlags (HY_ALIGN_LEFT);

    b1->SetAlignFlags (HY_ALIGN_RIGHT);

    b1->SetText       ("  OK  ");
    b1->SetButtonKind (HY_BUTTON_OK);
    b2->SetText       (" Cancel ");
    b2->SetButtonKind (HY_BUTTON_CANCEL);

    ProcessEvent (generateKeyboardFocusEvent (tb1->GetID()));

    DeleteObject (l1);
    DeleteObject (l2);
    DeleteObject (tb1);
    DeleteObject (tb2);
    DeleteObject (b1);
    DeleteObject (b2);

    _HYRect     dim = MinMaxWindowDimensions();

    SetWindowRectangle  (0,0,dim.bottom,dim.right);
    CenterWindow        (this);
}

//__________________________________________________________

bool    _HYTreePrintPrefs::ProcessEvent (_HYEvent* e)
{
    bool        done = false;
    _String     firstArg;
    long        i,f;

    if (e->EventClass()==_hyButtonPushEvent) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i = MatchComponentID (firstArg);

        if (i==2) { // OK
            *resW = defTreeHPages   = ((_HYTextBox*)components (0))->GetText().toNum();
            *resH = defTreeVPages   = ((_HYTextBox*)components (1))->GetText().toNum();
            *result = true;
            postWindowCloseEvent (GetID());
        } else if (i==3) { // Cancel
            *result = false;
            postWindowCloseEvent (GetID());
        }

        done = true;
    } else if (e->EventClass()==_hyTextEditChange) {
        firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
        i        = MatchComponentID (firstArg);

        if ((i==0)||(i==1)) {
            _HYButton* okButton     = (_HYButton*)components  (3);
            _String    textField    = ((_HYTextBox*)components (i))->GetText();

            _Parameter fieldValue   = textField.toNum();
            if ((fieldValue>=0.0)&&(fieldValue == (long) (fieldValue))) {
                done = true;
            }

            okButton->EnableButton (done);
        }

        done = true;
    }

    if (done) {
        DeleteObject (e);
        return true;
    }
    return _HYTWindow::ProcessEvent(e);
}

//_________________________________________________________________________
bool    OpenTreeFile (void)
{
    static  long openTreeCounter = 1;
    if (!PopUpFileDialog("HY-PHY Open: Choose a tree file")) {
        return false;
    }

    SetStatusLine ("Reading Tree");
    char c;

    FILE *f = doFileOpen (argFileName->getStr(), "r", true);
    if (!f) {
        return false;
    }

    while (!feof(f))
        // guess what type of file this is
    {
        if (!isspace(c = fgetc (f))) {
            if (c!='(') {
                if (c=='#') {
                    char checkNxs [6] = {toupper(fgetc (f)),toupper(fgetc (f)),toupper(fgetc (f)),toupper(fgetc (f)),toupper(fgetc (f)),'\0'},
                                        nxs[6] = "NEXUS";

                    if (strcmp (nxs,checkNxs) == 0)
                        // NEXUS
                    {
                        c = 3;
                        rewind (f);
                    }
                } else { // HBL
                    c = 2;
                }
            } else { // Newick
                c = 1;
            }

            break;
        }
    }

    if (feof(f)) {
        c=0;    // fail
    }

    long stashLMD = lastMatrixDeclared;
    lastMatrixDeclared = -1;

    if (c==2) {
        fclose (f);
        if (OpenBatchFile (false)) {
            ExecuteBatchFile();
        }
        lastMatrixDeclared = stashLMD;
        return true;
    } else if (c==1) {
        unsigned long cfp = (unsigned long)ftell(f),
                      ssi = _String::storageIncrement;

        fseek (f,0,SEEK_END);
        _String::storageIncrement = ftell(f)-cfp+1;
        _String treeString (16L,true);
        fseek (f,cfp,SEEK_SET);


        treeString<<'(';
        long levelCounter = 1;
        while (!feof(f) && levelCounter ) {
            c = fgetc (f);
            if (c=='(') {
                levelCounter++;
            } else if (c==')') {
                levelCounter--;
            }
            treeString<<c;
        }
        treeString.Finalize();
        _String::storageIncrement = ssi;

        if (levelCounter == 0) {
            _String treeName = _String("Tree_") & openTreeCounter++ ;

            _HYTreePanel* newTreePanel = new _HYTreePanel (treeName, treeString);
            newTreePanel->_Zoom (true);
            newTreePanel->BringToFront();
        } else {
            WarnError ("Imbalanced parentheses in the Newick string file");
        }
    } else if (c==3) {
        DeleteObject (ReadDataSetFile (f));
        if (!terminateExecution) {
            _String readTrees (1024L, true);
            readTrees << "UseModel(USE_NO_MODEL);\nfor (i=0; i<Rows(";
            readTrees << nexusFileTreeMatrix;
            readTrees << "); i=i+1)\n\t{\n\tdefineTreeString=\"Tree \"+";
            readTrees << nexusFileTreeMatrix;
            readTrees << "[i][0] +\" = \" +";
            readTrees << nexusFileTreeMatrix;
            readTrees << "[i][1]+\";\";\n\tExecuteCommands (defineTreeString);OpenWindow (TREEWINDOW,{{";
            readTrees << nexusFileTreeMatrix;
            readTrees <<"[i][0]}},\"(SCREEN_WIDTH-50)/2;(SCREEN_HEIGHT-50)/2;(20+i*20)%(SCREEN_WIDTH$2);(45+i*20)%(SCREEN_HEIGHT$2)\");\n\t}";
            readTrees.Finalize();
            _ExecutionList el (readTrees);
            el.Execute();
        }
        terminateExecution = false;
    }

    fclose (f);
    SetStatusLine ("Idle");
    lastMatrixDeclared = stashLMD;
    return true;
}

//_________________________________________________________________________
bool    OpenTextFile (void)
{
    if (!PopUpFileDialog("HY-PHY Open: Choose a text file")) {
        return false;
    }

    FILE *f = doFileOpen (argFileName->getStr(), "r");
    SetStatusLine ("Reading File");
    if (!f) {
        _String errMsg ("Could not read the file:");
        errMsg = errMsg & *argFileName;
        WarnError (errMsg);
        return false;
    }
    BufferToConsole ("\n\n");

    char  * buffer = MemAllocate (65536);
    while (!feof(f)) {
        buffer  [fread(buffer,1,65535,f)] = 0;
        BufferToConsole (buffer);
        handleGUI();
    }
    free (buffer);
    BufferToConsole ("\n\n");
    fclose (f);
    SetStatusLine ("Idle");
    return true;
}


//_________________________________________________________________________
bool    OpenDataFile (_String* defLoc)
{
    if (!PopUpFileDialog("HY-PHY Open: Choose a sequence file",defLoc)) {
        return false;
    }

    char    c            = GetPlatformDirectoryChar();
    _String dataFileName = *argFileName,
            tryName;

    dataFileName.Trim (dataFileName.FindBackwards(c,0,-1)+1,-1);
    c = '.';
    dataFileName.Trim (0,dataFileName.FindBackwards(c,0,-1)-1);
    dataFileName.ConvertToAnIdent();

    if (dataFileName.sLength==0) {
        dataFileName = "DataSet_";
    }
    FindUnusedObjectName (dataFileName,tryName,dataSetNamesList);
    dataFileName = tryName;
    _String BLCommand ("DataSet ");
    BLCommand = BLCommand & dataFileName & " = ReadDataFile (\""& * argFileName
                &"\");";
    _ExecutionList      ex (BLCommand);
    ex.Execute();
    SetStatusLine ("Idle");
    long k = dataSetNamesList.Find(&dataFileName);
    if (k>=0) {
        _DataSet* theDS = (_DataSet*)dataSetList (k);
        if (theDS->NoOfSpecies()>0 && theDS->NoOfColumns()>0) {
            // test the new data panel GUI
            _String dfHandler (128L, true);
            dfHandler << "if (Rows(";
            dfHandler << dataFilePartitionMatrix;
            dfHandler << ")==2){\n\tfor (counter = 0; counter <Columns(";
            dfHandler << dataFilePartitionMatrix;
            dfHandler << "); counter = counter+1){\n\tcommandString = \"DataSetFilter \"+";
            dfHandler << dataFilePartitionMatrix;
            dfHandler << "[0][counter] + \" = CreateFilter (";
            dfHandler << dataFileName;
            dfHandler << ",1,\\\"\"+";
            dfHandler << dataFilePartitionMatrix;
            dfHandler << "[1][counter] + \"\\\")\";\n\tExecuteCommands (commandString);\n}\n}";

            dfHandler.Finalize();

            ExecuteBLString (dfHandler,nil);

            _HYDataPanel* myTWindow = new _HYDataPanel(dataFileName, dataFileName);
            myTWindow->SetFilePath (*argFileName);
            myTWindow->BringToFront();

            _FString * v  = (_FString*)FetchObjectFromVariableByType (&dataFileTreeString, STRING);
            _Matrix  * m  = (_Matrix*)FetchObjectFromVariableByType (&nexusFileTreeMatrix,MATRIX);

            _String  treeBuilder (128L,1);

            if (v) {
                tryName = dataFileName & "_tree";
                FindUnusedObjectName (dataFileName,tryName,variableNames);
                treeBuilder << (_String ("Tree ")&tryName&'='&dataFileTreeString&';');
            }
            if (m) {
                _List treeStrings;
                m->FillInList (treeStrings);
                for (long s = 0; s < treeStrings.lLength; s+=2) {
                    tryName = dataFileName & "_" & *(_String*)treeStrings(s) ;
                    FindUnusedObjectName (dataFileName,tryName,variableNames);
                    treeBuilder << (_String ("Tree ")&tryName&'='&*(_String*)treeStrings(s+1) & ';');
                }
            }
            treeBuilder.Finalize();
            ExecuteBLString (treeBuilder,nil);

            return true;
            // execute the handler for opening data-set defined partitions
        }
        dataSetNamesList.Delete (k);
        dataSetList.Delete (k);
    }
    dataFileName = _String('"')&dataFileName&"\" does not appear to be a valid sequence data file";
    WarnError (dataFileName);
    terminateExecution = false;
    return true;
}

//____________________________________________________________________________________

void    NewModel (_String* model)
{
    long        mDim,
                gCode,
                pl1;

    _String     s1;
    _HYModelWindowDialog* newModel = new _HYModelWindowDialog (&mDim,&gCode,&s1,&pl1);
    if (model) {
        newModel->SetModelChoice (model);
    }

    newModel->BringToFront();
    while (windowObjectRefs.Find ((long)newModel)>=0) {
        handleGUI();
    }

    if (mDim>0) {
        _HYModelWindow* newMdl = SpawnNewModel (mDim, gCode, s1, pl1);
        if (newMdl) {
            newMdl->BringToFront();
        }
    }
}

//____________________________________________________________________________________

bool    TreePrintSetup (long& pW, long& pH, Ptr parent)
{
    bool        resB        = false;

    _HYTreePrintPrefs* tps = new _HYTreePrintPrefs (&pW,&pH,&resB,parent);

    tps->BringToFront();
    while (windowObjectRefs.Find ((long)tps)>=0) {
        handleGUI();
    }

    return      resB;
}


//____________________________________________________________________________________

bool    HandleCIDialog (char& opt, _Parameter& sig, _Parameter & sig2, Ptr parent)
{
    bool        resB        = false;

    _HYCIDialog* cid = new _HYCIDialog (&resB,&opt,&sig,&sig2,parent);

    cid->BringToFront();
    while (windowObjectRefs.Find ((long)cid)>=0) {
        handleGUI();
    }

    return      resB;
}

//____________________________________________________________________________________

bool    HandlePLDialog (_Parameter& lb, _Parameter& ub, long& ivals, bool& donorm, _Parameter mle, Ptr parent)
{
    bool        resB        = false;

    _HYLProfDialog* lpd = new _HYLProfDialog (&resB,&lb,&ub,&ivals,&donorm,mle,lb,ub,parent);

    lpd->BringToFront();
    while (windowObjectRefs.Find ((long)lpd)>=0) {
        handleGUI();
    }

    return      resB;
}

//____________________________________________________________________________________

void    NewChartWindow (void)
{
    bool        resB        = false;
    _String     windowTitle ("New Chart");
    long r,
         c;

    _HYNewChart* nc = new _HYNewChart (&r,&c,&resB,&windowTitle);

    nc->BringToFront();
    while (windowObjectRefs.Find ((long)nc)>=0) {
        handleGUI();
    }

    if  (resB) {
        _List            columnHeaders;
        _HYChartWindow * newChart = nil;
        _Matrix        * vMx;

        if (r<=0) {
            _Variable* me = FetchVar (-r);
            if ((!me) || (me->ObjectClass() != MATRIX)) {
                return;
            }

            vMx =  (_Matrix*)me->Compute();
            c = vMx->GetVDim();
        }

        for   (long k=1; k<=c; k++) {
            _String thisCol = _String ("Series ") & k;
            columnHeaders && & thisCol;
        }

        if (r<=0) {
            newChart = new _HYChartWindow (windowTitle, columnHeaders, *vMx, nil);
        } else {
            _Matrix data (r,c,false,true);
            newChart = new _HYChartWindow (windowTitle, columnHeaders, data, nil);
        }
        checkPointer (newChart);
        newChart->BringToFront();
    }
}


//_________________________________________________________________________
bool    OpenTable (void)
{
    static bool        resC = false;
    static char        sepC = ',';

    if (!PopUpFileDialog("HY-PHY Open: Choose a table file")) {
        return false;
    }

    SetStatusLine  ("Reading Data Table");

    _Matrix        data;
    _List          names;

    if (!resC) {

        _String        res,
                       prompt ("Enter field separator"),
                       cprompt("Don't show this dialog again");

        _List          pulldown,
                       charmeaning;

        res = "Comma";
        pulldown && & res;
        res = "Tab";
        pulldown && & res;
        res = ',';
        charmeaning && & res;
        res = '\t';
        charmeaning && & res;
        res = ',';

        long           msel;

        if (!EnterStringDialogWithPulldown(res,prompt,cprompt,msel,pulldown,&charmeaning,resC,-1,nil)) {
            SetStatusLine  ("Idle");
            return false;
        }

        sepC = res.getChar(0);

        if (ReadDataFromFile (*argFileName, sepC, data, names)) {
            _HYChartWindow * cw = new _HYChartWindow (argFileName->Cut (argFileName->FindBackwards (":",0,-1)+1,-1), names, data, nil);
            checkPointer (cw);
            cw->BringToFront();
        }
    }

    SetStatusLine  ("Idle");
    return true;
}

//_________________________________________________________________________
bool    OpenModelFile (_String* defLoc)
{
    if (!PopUpFileDialog("HY-PHY Open: Choose a model file",defLoc)) {
        return false;
    }

    return OpenModelFromFile (*argFileName);
}

//_________________________________________________________________________
bool    OpenDatabaseFile (_String* defLoc)
{
    if (!PopUpFileDialog("HY-PHY Open: Choose an SQLite database file",defLoc)) {
        return false;
    }

    _HYDBWindow * testDBWindow = new _HYDBWindow (_String("SQLite file: ")&*argFileName,argFileName);
    testDBWindow->BringToFront();
    return true;
}

//_________________________________________________________________________
bool    NewDatabaseFile (_String* defLoc)
{
    _String        str1 ("Save new database file to:"),
                   str2 (defLoc?*defLoc:empty),
                   filePath;

    _List          menuOptions;
    filePath = "SQLite file";
    menuOptions && & filePath;
    filePath = empty;


    long menuSel = SaveFileWithPopUp (filePath, str1 ,str2 , empty, menuOptions);

    if (menuSel>=0) {
        FILE*   outFile = doFileOpen (filePath.sData,"w");
        if (!outFile) {
            str1 = filePath & " could not be opened for writing.";
            ProblemReport (str1);
            return false;
        }
        fclose (outFile);


        _HYDBWindow * testDBWindow = new _HYDBWindow (_String("SQLite file: ")&filePath,&filePath);
        testDBWindow->BringToFront();
        return true;
    }
    return false;
}

//_________________________________________________________________________
void    ShowMessagesLog (void)
{
    if (globalMessageFile) {
        BufferToConsole ("\n\n_________________________<Log of Messages>______________________________\n\n");
        _String mF (globalMessageFile);
        rewind (globalMessageFile);
        StringToConsole(mF);
        BufferToConsole ("\n\n_________________________<End of Messages>______________________________\n\n");
        fseek (globalMessageFile,0,SEEK_END);
    } else {
        BufferToConsole ("\nLog file was not open at startup (probably because HYPHY startup volume is read-only\n");
    }
}

//_________________________________________________________________________

void    ExecuteSelection (void)
{
    _HYTextBox * ob = (_HYTextBox*)hyphyConsoleWindow->GetObject (0);
    _String    * sel;
    ob->StoreText (sel,true);
    if (sel->sLength) {
        _ExecutionList exl (*sel);
        exl.Execute();
        terminateExecution = false;
    }
    DeleteObject (sel);
}

//_________________________________________________________________________

long    FindWindowByName (_String & s)
{
    for (long i=0; i<windowObjects.lLength; i++) {
        _HYWindow* thisWindow = (_HYWindow*)windowObjectRefs(i);
        if (s.Equal(&thisWindow->_GetTitle())) {
            return i;
        }
    }
    return -1;
}

//_________________________________________________________________________

_HYGuiObject*   FindWindowByNameAndOpen (_String & s)
{
    long f = FindWindowByName (s);
    if (f>=0) {
        _HYPlatformWindow* thisWindow = (_HYPlatformWindow*)windowObjects(f);
        thisWindow->_Activate();
        return (_HYGuiObject*)windowObjectRefs(f);
    }
    return nil;
}

//_________________________________________________________________________

_HYGuiObject*   FindWindowByID (long ID)
{
    for (long i=0; i<windowObjects.lLength; i++) {
        _HYWindow* thisWindow = (_HYWindow*)(windowObjectRefs.lData[i]);
        if (thisWindow->MatchID(ID)) {
            return thisWindow;
        }
    }
    return nil;
}
//____________________________________________________________________________________________

#ifdef __MAC__
#include "Timer.h"

void    yieldCPUTime(void)
{
    handleGUI(true);
}
#endif

#ifdef __WINDOZE__
#include        "Windows.h"
#include        "preferences.h"
#include        "HYSharedMain.h"
#include        "HYPlatformWindow.h"

extern  bool    hyphyExiting;

void            yieldCPUTime     (void)
{
    MessageLoop();
    if (hyphyExiting) {
        WritePreferences    ();
        ExitProcess(0);
    }

    while (isSuspended) {
        MessageLoop(false,false);
        if (hyphyExiting) {
            WritePreferences    ();
            ExitProcess         (0);
        }
    }
}
#endif

#ifdef  __HYPHY_GTK__

#include <gtk/gtk.h>
void    yieldCPUTime (void)
{
    while (gtk_events_pending ()) {
        gtk_main_iteration();
    }
}

#endif



//EOF