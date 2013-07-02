/*
    Mac OS Portions of the console window class

    Sergei L. Kosakovsky Pond, December 2003.
*/

#include "HYConsoleWindow.h"
#include "HYTextBox.h"
#include <ControlDefinitions.h>
#include <MacTextEditor.h>
#include "HYPlatformGraphicPane.h"
#include "parser.h"
#include <Timer.h>

extern   PixPatHandle   whiteFill,
         blueFill,
         statusBarFill;

extern   RGBColor       _BLACK_;

extern   CIconHandle    redButtonIcon,
         orangeButtonIcon,
         yellowButtonIcon;

extern  _String         VerbosityLevelString;

UnsignedWide            lastMeasure = {0,0};

//__________________________________________________________________

void _HYConsoleWindow::_SetMenuBar(void)
{
    //BufferToConsole ("_HYConsoleWindow::_SetMenuBar\n");
    _HYWindow::_SetMenuBar();

    MenuHandle  editMenuH  = GetMenuHandle (130);

    if (CountMenuItems (editMenuH) < 13) {
        InsertMenuItem (editMenuH,"\p(-", 6);
        InsertMenuItem (editMenuH,"\pFind.../F", 7);
        InsertMenuItem (editMenuH,"\p(Find Again/G", 8);

        InsertMenuItem (editMenuH,"\p(Redo/Y", 1);

        EnableMenuItem (editMenuH,12);
        EnableMenuItem (editMenuH,13);

        _UpdateEditMenu();
    }
    InvalMenuBar();
}

//__________________________________________________________________

void _HYConsoleWindow::_UnsetMenuBar(void)
{
    //BufferToConsole ("_HYConsoleWindow::_UnsetMenuBar\n");
    MenuHandle      editMenuH  = GetMenuHandle (130);
    if (CountMenuItems (editMenuH) >= 13) {
        DeleteMenuItem (editMenuH,8);
        DeleteMenuItem (editMenuH,8);
        DeleteMenuItem (editMenuH,8);
        DeleteMenuItem (editMenuH,2);
        _HYWindow::_UnsetMenuBar();
    }
}

//__________________________________________________________________

void _HYConsoleWindow::_UpdateEditMenu (void)
{

    _HYTextBox* txb = (_HYTextBox*)GetObject (0);

    bool        haveSelection = !TXNIsSelectionEmpty(txb->txn),
                canPaste      = TXNIsScrapPastable();

    _String     undoAction (mlteCommandString(txb->txn,true)),
                redoAction (mlteCommandString(txb->txn,false));

    if ((((bool)editOptions&HY_CONSOLE_CAN_COPY)!=haveSelection)||
            (((bool)editOptions&HY_CONSOLE_CAN_PASTE)!=canPaste)||
            (((bool)editOptions&HY_CONSOLE_CAN_UNDO)!=(bool)undoAction.sLength)||
            (((bool)editOptions&HY_CONSOLE_CAN_REDO)!=(bool)redoAction.sLength)) {
        MenuHandle  t   = GetMenuHandle (130);

        editOptions = 0;
        if (haveSelection) {
            editOptions |= HY_CONSOLE_CAN_COPY;
            EnableMenuItem (t,4);
            EnableMenuItem (t,5);
            EnableMenuItem (t,7);
        } else {
            DisableMenuItem (t,4);
            DisableMenuItem (t,5);
            DisableMenuItem (t,7);
        }

        if (canPaste) {
            EnableMenuItem (t,6);
            editOptions |= HY_CONSOLE_CAN_PASTE;
        } else {
            DisableMenuItem (t,6);
        }

        if (undoAction.sLength) {
            editOptions |= HY_CONSOLE_CAN_UNDO;
            EnableMenuItem (t,1);
        } else {
            DisableMenuItem (t,1);
        }

        if (redoAction.sLength) {
            editOptions |= HY_CONSOLE_CAN_REDO;
            EnableMenuItem (t,2);
        } else {
            DisableMenuItem (t,2);
        }
    }
}

//__________________________________________________________________

bool _HYConsoleWindow::_Close (Ptr data)
{
    return _HYPlatformWindow::_Close(data);
}

//__________________________________________________________________

bool        _HYConsoleWindow::_ProcessMenuSelection (long msel)
{
    long        menuChoice = msel&0x0000ffff;

    HiliteMenu(0);
    InvalMenuBar();

    switch (msel/0xffff) {
    case 129: { // undo
        switch (menuChoice) {
        case 4:
            SaveConsole();
            return true;

        case 7:
            _DoPageSetup();
            return true;

        case 8:
            _DoPrint();
            return true;
        }
        break;
    }
    case 130: { //
        if (menuChoice<14) { //
            _HYTextBox * ob = (_HYTextBox*) components (0);
            switch (menuChoice) {
            case 1:
                ob->_DoUndo(true);
                break;
            case 2:
                ob->_DoRedo(true);
                break;
            case 4:
                ob->_DoCut(true);
                break;
            case 5:
                ob->_DoCopy(true);
                break;
            case 6:
                ob->_DoPaste(true);
                break;
            case 7:
                ob->_DoClear(false,true);
                break;
            case 9:
                DoFind (true);
                break;
            case 10:
                DoFind (false);
                break;
            case 12:
                ob->_DoSelectAll(true);
                break;
            case 13:
                ob->_DoClear(true,true);
                break;
            }
            return true;
        }
    }
    }

    return _HYTWindow::_ProcessMenuSelection(msel);
}

//__________________________________________________________________

bool _HYConsoleWindow::_ProcessOSEvent (Ptr vEvent)
{
    //EventRecord*  theEvent = (EventRecord*)vEvent;

    return _HYTWindow::_ProcessOSEvent (vEvent);
}

//__________________________________________________________________
void    _HYConsoleWindow::_DoPageSetup          (void)
{
    TXNPageSetup (((_HYTextBox*)components(0))->txn);
}

//__________________________________________________________________
void    _HYConsoleWindow::_DoPrint          (void)
{
    Deactivate();
    TXNPrint (((_HYTextBox*)components(0))->txn);
}


//__________________________________________________________________

void _HYConsoleWindow::_PaintStatusBar(Ptr, bool force)
{
    _Parameter      vL;
    checkParameter (VerbosityLevelString, vL, 0.0);
    if (vL<-0.5 && !force) {
        UnsignedWide curMeasure;
        Microseconds (&curMeasure);
        _Parameter diff = 0.000001*((curMeasure.hi-lastMeasure.hi)*4294967296.0+(curMeasure.lo-lastMeasure.lo));
        if (diff < -vL) {
            return;
        }
         lastMeasure = curMeasure;
    }

    GrafPtr saveport;

    Rect    box,
            box2,
            box3;

    static  PixPatHandle    blueFill  = GetPixPat(132),
                            whiteFill = GetPixPat(128);


    GWorldPtr   offScreenMap = nil;

    GetPort(&saveport);

    Rect wr;
#ifdef OPAQUE_TOOLBOX_STRUCTS
    GetWindowBounds (theWindow, kWindowGlobalPortRgn, &wr);
    OffsetRect (&wr,-wr.left,-wr.top);
#else
    wr = theWindow->portRect;
#endif

    SetRect(&box,0,wr.bottom-15,wr.right - 15,wr.bottom);

    if (NewGWorld(&offScreenMap,0,&box,0,GetMainDevice(),noNewDevice)!=noErr) {
        return;
    }

    LockPixels (GetGWorldPixMap(offScreenMap));
    box3 = box;

#ifdef OPAQUE_TOOLBOX_STRUCTS
    SetPort (offScreenMap);
#else
    SetPort ((GrafPort*)offScreenMap);
#endif

    box.top = 0;
    box.bottom = 15;
    box2 = box;
    InsetRect(&box, -1, -1);
    InsetRect(&box, 1, 1);
    box.right = 30;
    FillCRect (&box, blueFill);
    box.left = 30;
    box.right = box2.right;
    FillCRect  (&box,statusBarFill);

    TextMode (srcBic);
    TextFont(kFontIDGeneva );
    TextSize(9);
    TextFace(bold);

    MoveTo(4, box.bottom - 3);
    DrawText(cState.sData,0,cState.sLength);

    MoveTo(32,  box.bottom - 3);
    TextMode    (srcOr);
    TextFace    (normal);
    DrawText    (fileName.getStr(),0,fileName.sLength);
    box.left  = 150;
    box.right = 190;
    FillCRect   (&box, blueFill);
    TextMode    (srcBic);
    TextFace    (bold);

    MoveTo(152, box.bottom-3);
    DrawText(cTask.sData,0,cTask.sLength);
    box.left = 190;
    TextMode (srcOr);
    TextFace(normal);
    MoveTo(192, box.bottom-3);

    if (inputStatus == 1) {
        DrawText(cInput.getStr(),0,cInput.sLength);
    } else {
        DrawText(action.getStr(),0,action.sLength);
    }

    box.right = box2.right - 20;
    box.left  = box.right  - 50;
    FillCRect (&box, blueFill);
    TextMode (srcBic);
    MoveTo(box.left+2, box.bottom-3);
    DrawText(timer.getStr(),0,timer.sLength);

    if (percentDone >= 0 || percentDone == -HY_SL_DONE) {
        TextMode (srcOr);
        box.right = box.left-5;
        box.left-=75;
        box.top+=1;
        box.bottom-=2;
        FillCRect (&box,whiteFill);

        ThemeTrackDrawInfo  pB;
        pB.kind     = kThemeProgressBar;
        pB.bounds   = box;
        pB.min = 0;
        pB.max = 100;
        pB.value                    =  percentDone>=0?percentDone:100;
        pB.attributes               = kThemeTrackHorizontal;
        pB.enableState              = kThemeTrackActive;
        pB.trackInfo.progress.phase = 0;
        DrawThemeTrack (&pB,nil,nil,0);
        TextFont(kFontIDTimes);
        MoveTo(box.left+28, box.bottom-2);

        _String pLine;
        if (percentDone >= 0) {
            pLine = _String(percentDone) & "%";
        } else {
            pLine = "DONE";
        }
        DrawText(pLine.getStr(),0,pLine.sLength);
    }

    box = box2;
    box.right -= 3;
    box.left  = box.right - 12;
    box.top +=2;
    box.bottom --;

    CIconHandle cIcon;

    if (echoFileRef) {
        if (echoStatus == 1) {
            cIcon = redButtonIcon;
        } else {
            cIcon = orangeButtonIcon;
        }
    } else {
        cIcon = yellowButtonIcon;
    }

    PlotCIcon (&box,cIcon);


    MoveTo(0, 0);
    Line(wr.right - 15, 0);

    SetPortWindowPort(theWindow);
#ifdef OPAQUE_TOOLBOX_STRUCTS
    CopyBits (GetPortBitMapForCopyBits(offScreenMap),GetPortBitMapForCopyBits(GetWindowPort(theWindow)), &box2, &box3, srcCopy,nil);
    RgnHandle   dirtyRgn = NewRgn ();
    RectRgn (dirtyRgn, &box3);
    QDFlushPortBuffer (GetWindowPort (theWindow), nil);
    DisposeRgn (dirtyRgn);
#else
    CopyBits (&(GrafPtr(offScreenMap)->portBits),&(theWindow->portBits), &box2, &box3, srcCopy,nil);
#endif

    DisposeGWorld (offScreenMap);

    SetPort(saveport);
    yieldCPUTime();
}

//EOF