/*
    Mac OS Portions of the data panel class

    Sergei L. Kosakovsky Pond, Spring 2000 - December 2002.
*/

#include "HYTreePanel.h"
#include "HYDataPanel.h"
#include "HYUtils.h"

#include "likefunc.h"

//__________________________________________________________________


bool        _HYDataPanel::_ProcessMenuSelection (long msel)
{
    if (_HYWindow::_ProcessMenuSelection(msel)) {
        return true;
    }

    long        menuChoice = msel&0x0000ffff;
    MenuHandle  treeMenu;
    _HYSequencePane* sp = (_HYSequencePane*)GetObject (0);
    _HYSequencePane* sp2 =(_HYSequencePane*)GetObject (4);
    _String     prompt;
    bool        done = false;

    switch (msel/0xffff) {
    case 129: { // file menu
        if (menuChoice==4) { // save
            SaveDataPanel (savePath.sLength);
            done = true;
        } else if (menuChoice==8) { // print
            _PrintData();
            done = true;
        }
        HiliteMenu(0);
        if (done) {
            return true;
        }
        break;
    }
    case 130: { // edit
        HiliteMenu(0);
        if (menuChoice==4) { // copy
            _CopySelectionToClipboard   ();
        } else if (menuChoice==6) { // delete selection
            //DeleteCurrentSelection();
            return true;
        } else if (menuChoice==3) { // cut selection
            //CutSelectionToClipboard ();
            return true;
        } else if (menuChoice==5) { // paste selection
            //PasteClipboardTree();
            return true;
        } else if (menuChoice==8) { // select all
            sp->SelectAll(true);
            return true;
        } else if (menuChoice==11) { // Find
            FindFunction();
            return true;
        } else if (menuChoice==12) { // Find
            HandleSearchAndReplace();
            return true;
        }
        return false;
    }
    case HY_DATAPANEL_MENU_ID: {
        switch (menuChoice) {
        case 1:
            SelectPartition();
            break;
        case 2:
            if (sp->selection.lLength) {
                CreatePartition (sp->selection,1,true);
            }
            break;
        case 3:
            InvertSelection();
            break;
        case 11:
            PartitionPropsMenu ();
            break;
        case 12:
            InputPartitionString ();
            break;
        case 15:
            SimulateDataSet (0,true);
            break;
        case 17:
            HandleFontChange    ();
            break;
        }
        HiliteMenu(0);
        return true;
    }
    case HY_DATAPANEL_MENU_ID+1: {
        switch (menuChoice) {
        case 1:
            BuildLikelihoodFunction();
            break;
        case 3:
            OptimizeLikelihoodFunction();
            break;
        case 5:
            DisplayParameterTable ();
            break;
        case 7:
            OpenGeneralBSWindow ();
            break;
        }
        HiliteMenu(0);
        return true;
    }
    case HY_DATAPANEL_HMENU_ID: {
        treeMenu = GetMenuHandle (HY_DATAPANEL_HMENU_ID);
        long             newBlockSize;
        if (menuChoice==1) {
            newBlockSize = 9;
        } else {
            newBlockSize = 10;
        }

        if (sp->blockWidth!=newBlockSize) {
            sp->blockWidth = newBlockSize;
            sp2->blockWidth = newBlockSize;
            sp->BuildPane();
            sp->_MarkForUpdate();
            sp2->BuildPane();
            sp2->_MarkForUpdate();
            SetItemMark(treeMenu,menuChoice==1?2:1,noMark);
            SetItemMark(treeMenu,menuChoice,checkMark);
            InvalMenuBar();
        }
        HiliteMenu(0);
        return true;
    }

    case HY_DATAPANEL_HMENU_ID+1: {
        treeMenu = GetMenuHandle (HY_DATAPANEL_HMENU_ID+1);
        bool       newDisplay;
        if (menuChoice==1) {
            newDisplay = false;
        } else {
            newDisplay = true;
        }

        if (sp->showDots!=newDisplay) {
            sp->showDots = newDisplay;
            sp->BuildPane();
            sp->_MarkForUpdate();
            SetItemMark(treeMenu,menuChoice==1?2:1,noMark);
            SetItemMark(treeMenu,menuChoice,checkMark);
            InvalMenuBar();
        }
        HiliteMenu(0);
        return true;
    }

    case HY_DATAPANEL_HMENU_ID+2: {
        treeMenu = GetMenuHandle (HY_DATAPANEL_HMENU_ID+2);
        if (menuChoice<=3) {
            unsigned char newDisplay;
            switch (menuChoice) {
            case 1:
                newDisplay = HY_SEQUENCE_PANE_NAMES_NONE;
                break;
            case 2:
                newDisplay = HY_SEQUENCE_PANE_NAMES_SHORT;
                break;
            case 3:
                newDisplay = HY_SEQUENCE_PANE_NAMES_ALL;

            }

            if ((sp->nameDisplayFlags&HY_SEQUENCE_PANE_NAMES_MASK)!=newDisplay) {
                SetItemMark(treeMenu,(sp->nameDisplayFlags&HY_SEQUENCE_PANE_NAMES_MASK)+1,noMark);
                sp->SetNameDisplayMode(newDisplay,true);
                sp2->SetNameDisplayMode(newDisplay,true);
                BuildThermometer();
                BuildMarksPane();
                SetItemMark(treeMenu,menuChoice,checkMark);
                InvalMenuBar();
            }
        } else {
            if (menuChoice==5) { // alphabetize
                sp->AlphabetizeSpecies();
            } else if (menuChoice==6) {
                sp->RevertFileOrder();
            } else if (menuChoice==7) {
                sp->CleanUpSequenceNames();
            }
        }
        HiliteMenu(0);
        return true;
    }
    case HY_DATAPANEL_HMENU_ID+3: { // omitted sequences
        RestoreOmittedSequence(menuChoice-3);
        HiliteMenu(0);
        return true;
    }

    case HY_DATAPANEL_HMENU_ID+4: { // status lines
        if (AdjustStatusLine (menuChoice-1)) {
            CheckMenuItem(GetMenuHandle (HY_DATAPANEL_HMENU_ID+4),menuChoice,true);
        } else {
            CheckMenuItem(GetMenuHandle (HY_DATAPANEL_HMENU_ID+4),menuChoice,false);
        }
        HiliteMenu(0);
        return true;
    }

    case HY_DATAPANEL_HMENU_ID+5: { // likelihood display
        ComputeLikelihoodFunction (menuChoice-1);
        HiliteMenu(0);
        return true;
    }

    case HY_DATAPANEL_HMENU_ID+6: { // simulate data set
        SimulateDataSet (menuChoice-1);
        HiliteMenu(0);
        return true;
    }

    case HY_DATAPANEL_HMENU_ID+7: { // save/save as
        SaveDataPanel (menuChoice==2);
        HiliteMenu(0);
        return true;
    }

    case HY_DATAPANEL_HMENU_ID+8: { // infer
        InferTopologies (menuChoice==2);
        _VerifyInferMenu ();
        HiliteMenu(0);
        return true;
    }
    case HY_DATAPANEL_HMENU_ID+9: { // dataProcs
        ExecuteProcessor (menuChoice-1);
        return true;
    }

    HiliteMenu(0);
    InvalMenuBar();
    }

    return false;
}

//__________________________________________________________________

void        _HYDataPanel::_PaintThermRect(bool update)
{
    navRect = ComputeNavRect();
    _HYCanvas* theCanvas = (_HYCanvas*)GetObject (1);
    GrafPtr savedPort;
    GetPort(&savedPort);
#ifdef OPAQUE_TOOLBOX_STRUCTS
    SetPort(GetWindowPort(theWindow));
#else
    SetPort(theWindow);
#endif
    Rect r;
    r.left      = navRect.left+theCanvas->rel.left+thermRect.left+1;
    r.right     = navRect.right+theCanvas->rel.left+thermRect.left-1;
    r.top       = navRect.top+theCanvas->rel.top+thermRect.top+1;
    r.bottom    = navRect.bottom+theCanvas->rel.top+thermRect.top-1;

    RGBColor    saveColor,
                newColor = {255*256,151*256,51*256};
    PenState ps;
    GetPenState (&ps);
    GetForeColor (&saveColor);
    RGBForeColor (&newColor);
    PenSize (2,2);
    FrameRect (&r);
    RGBForeColor (&saveColor);

    if (update) {
        RgnHandle   oldClip, newClip;
        oldClip = NewRgn();
        newClip = NewRgn();
        RectRgn (oldClip,&r);
        InsetRect (&r,2,2);
        RectRgn (newClip,&r);
        DiffRgn (oldClip,newClip,newClip);
        GetClip (oldClip);
        DiffRgn (oldClip,newClip,newClip);
        SetClip (newClip);
        _HYRect rect;
        rect.left = componentL.lData[1];
        rect.right = componentR.lData[1];
        rect.top = componentT.lData[1];
        rect.bottom = componentB.lData[1];
        theCanvas->_Paint((char*)&rect);
        SetClip (oldClip);
        DisposeRgn (oldClip);
        DisposeRgn (newClip);
    }
    SetPenState (&ps);
    _PaintLFStatus ();
    SetPort(savedPort);
}

//__________________________________________________________________

void        _HYDataPanel::_PaintLFStatus(void)
{
    GrafPtr savedPort;
    GetPort(&savedPort);
#ifdef OPAQUE_TOOLBOX_STRUCTS
    SetPort(GetWindowPort(theWindow));
#else
    SetPort(theWindow);
#endif
    if (lfID<0) {
        _SimpleList goodP;


        bool    paintOrange = GenerateGoodPartitions (goodP);


        if (goodP.lLength) {
            _PaintTheCircle (paintOrange?orangeButtonIcon:yellowButtonIcon,theWindow);
        } else {
            _PaintTheCircle (redButtonIcon,theWindow);
        }
    } else {
        _PaintTheCircle (greenButtonIcon,theWindow);
    }
    SetPort(savedPort);
}

//__________________________________________________________________

void        _HYDataPanel::_PrintData(void)
{
    _HYSequencePane* sp = (_HYSequencePane*)GetObject (0);

    GrafPtr     savePort;
#ifdef      TARGET_API_MAC_CARBON
    PMRect prRect;
#else
    TPrStatus   prStatus;
    TPPrPort    printPort;
    OSErr       err;
#endif

#ifdef TARGET_API_MAC_CARBON
    OSStatus theStatus;
    Boolean isAccepted;

    PMPrintSession hyPC;
    theStatus = PMCreateSession(&hyPC);
    if (theStatus != noErr) {
        return;
    }
#endif

    if (!InitPrint(hyPC)) {
        _String errMsg ("Could not initialize printing variables.");
        WarnError (errMsg);
        terminateExecution = false;
        return;
    }

    GetPort(&savePort);

#ifdef TARGET_API_MAC_CARBON
    if (gPrintSettings != kPMNoPrintSettings) {
        theStatus = PMSessionValidatePrintSettings(hyPC,gPrintSettings, kPMDontWantBoolean);
    } else {
        theStatus = PMCreatePrintSettings(&gPrintSettings);

        if ((theStatus == noErr) && (gPrintSettings != kPMNoPrintSettings)) {
            theStatus = PMSessionDefaultPrintSettings(hyPC,gPrintSettings);
        }
    }

    if (theStatus == noErr) {
        theStatus = PMSessionPrintDialog(hyPC,gPrintSettings, gPageFormat, &isAccepted);

        if (isAccepted) {
            theStatus = PMGetAdjustedPageRect(gPageFormat, &prRect);
            if (theStatus != noErr) {
                return;
            }

            theStatus = PMSessionBeginDocument(hyPC,gPrintSettings, gPageFormat);
            if (theStatus != noErr) {
                return;
            }

            long     printW    = prRect.right-prRect.left-2,
                     printH    = prRect.bottom-prRect.top-2;

            UInt32   startPage,
                     endPage;

            PMGetFirstPage (gPrintSettings,&startPage);
            PMGetLastPage  (gPrintSettings,&endPage);
#else
    PrOpen();
    if (err=PrError()) {
        _String errMsg ("Could not print the data set. Error Code:");
        errMsg = errMsg & (long)err;
        WarnError (errMsg);
        terminateExecution = false;
        return;
    }

    if (PrJobDialog(prRecHdl)) {
        printPort = PrOpenDoc(prRecHdl, nil, nil);
        SetPort((GrafPtr)printPort);
        long     printW = (*prRecHdl)->prInfo.rPage.right-2,
                 printH = (*prRecHdl)->prInfo.rPage.bottom-2,
                 startPage = (*prRecHdl)->prJob.iFstPage,
                 endPage = (*prRecHdl)->prJob.iLstPage;
#endif
            long
            vOffset = sp->GetSlotHeight()*sp->speciesIndex.lLength+sp->GetSlotHeight()*3/2+1+20*(dataPartitions.lLength>0),
            cOffset = (printW-sp->headerWidth)/sp->charWidth,
            pageShift = printH/vOffset,
            sC = sp->startColumn,
            lC = sp->endColumn,
            sR = sp->startRow,
            lR = sp->endRow,
            sH = sp->settings.bottom,
            pageCount;

            _HYColor c1 = sp->backColor,
                     c2 = sp->headerColor,
                     bc = {0,0,0};


            sp->backColor   = (_HYColor) {
                255,255,255
            };
            sp->headerColor = (_HYColor) {
                255,255,255
            };
            cOffset     -= ((cOffset/sp->blockWidth)*2)/sp->charWidth;
            cOffset     = (cOffset/sp->blockWidth)*sp->blockWidth;

            pageShift   *= cOffset;

            if (sp->columnStrings.lLength%pageShift==0) {
                pageShift = sp->columnStrings.lLength / pageShift;
            } else {
                pageShift = sp->columnStrings.lLength / pageShift + 1;
            }

            if (endPage > pageShift) {
                endPage = pageShift;
            }

            sp->startColumn = 0;
            sp->endColumn = cOffset;
            sp->startRow = 0;
            sp->endRow = sp->speciesIndex.lLength;
            sp->settings.bottom = vOffset+5+HY_SCROLLER_WIDTH;

            for (pageCount = 1; pageCount<startPage; pageCount++) {
#ifdef TARGET_API_MAC_CARBON
                theStatus = PMSessionBeginPage(hyPC,gPageFormat,NULL);
                if (theStatus != noErr) {
                    break;
                }
                theStatus = PMSessionEndPage(hyPC);
                if (theStatus != noErr) {
                    break;
                }
#else
                PrOpenPage      (printPort, nil);
                PrClosePage     (printPort);
#endif
                sp->endColumn       +=  printH/vOffset * cOffset;
                sp->startColumn     +=  printH/vOffset * cOffset;
            }

            Rect         hangover = {0,
                                     sp->headerWidth + cOffset * sp->charWidth + (cOffset * 2)/sp->blockWidth + 1,
                                     vOffset,
                                     printW
                                    },

                         frame = {0,0,vOffset+1,hangover.left};

            if (sp->startColumn< sp->columnStrings.lLength)
                for (pageCount = startPage; pageCount<=endPage; pageCount++) {
                    pageShift       = vOffset;
#ifdef TARGET_API_MAC_CARBON
                    theStatus = PMSessionBeginPage(hyPC,gPageFormat, NULL);
                    if (theStatus != noErr) {
                        break;
                    }
                    GrafPtr ppPort;
                    PMSessionGetGraphicsContext (hyPC, NULL, (void**)&ppPort);
                    SetPort (ppPort);
#else
                    PrOpenPage      (printPort, nil);
#endif
                    SetOrigin       (0,0);
                    sp->_HYGraphicPane::SetFont     (sp->GetFont());
                    sp->SetColor    (sp->GetColor());
                    while (pageShift < printH) {
                        sp->BuildPane   (false);
                        sp->SetColor (bc);
                        PenSize (1,1);
                        //EraseRect (&hangover);
                        if (dataPartitions.lLength) {
                            _HYRect          daFrame;
                            daFrame.top    = frame.bottom;
                            daFrame.bottom = frame.bottom + 20;
                            daFrame.left   = frame.left   + HY_SEQUENCE_PANE_CHAR_SPACING/2 + sp->headerWidth;
                            daFrame.right  = daFrame.left + cOffset*sp->charWidth + 2*(cOffset/sp->blockWidth) - HY_SCROLLER_WIDTH;
                            BuildThermometer (&daFrame);
                        }
                        //FrameRect (&frame);
                        sp->startColumn  =  sp->endColumn;
                        sp->endColumn   +=  cOffset;
                        SetOrigin           (0,-pageShift);
                        pageShift       +=  vOffset;

                        if (sp->startColumn>=sp->columnStrings.lLength) {
                            break;
                        }
                    }
#ifdef TARGET_API_MAC_CARBON
                    theStatus = PMSessionEndPage(hyPC);
                    if (theStatus != noErr) {
                        break;
                    }
#else
                    PrClosePage     (printPort);
#endif
                }

            sp->startColumn = sC;
            sp->endColumn   = lC;
            sp->startRow    = sR;
            sp->endRow      = lR;
            sp->backColor   = c1;
            sp->headerColor = c2;
            sp->settings.bottom = sH;


#ifdef TARGET_API_MAC_CARBON
            theStatus = PMSessionEndDocument(hyPC);
            SetPort(savePort);
            if (theStatus == noErr) {
                if (gFlattenedFormat != NULL) {
                    DisposeHandle(gFlattenedFormat);
                    gFlattenedFormat = NULL;
                }

                theStatus = PMFlattenPageFormat(gPageFormat, &gFlattenedFormat);
            }

            if (theStatus == noErr) {
                if (gFlattenedSettings != NULL) {
                    DisposeHandle(gFlattenedSettings);
                    gFlattenedSettings = NULL;
                }

                theStatus = PMFlattenPrintSettings(gPrintSettings, &gFlattenedSettings);
            }

            if (gPageFormat != kPMNoPageFormat) {
                theStatus = PMRelease(gPageFormat);
                gPageFormat = kPMNoPageFormat;
            }

            if (gPrintSettings != kPMNoPrintSettings) {
                theStatus = PMRelease(gPrintSettings);
                gPrintSettings = kPMNoPrintSettings;
            }

            theStatus = PMRelease(hyPC);

#else
            PrCloseDoc(printPort);
            if (((*prRecHdl)->prJob.bJDocLoop = bSpoolLoop) && (!PrError() ) ) {
                PrPicFile(prRecHdl, nil, nil, nil, &prStatus);
            }
#endif
        }
#ifdef TARGET_API_MAC_CARBON
        else {
            theStatus = PMRelease(hyPC);
        }
#endif

#ifdef TARGET_API_MAC_CARBON
    }
#else
        PrClose();
        SetPort(savePort);
#endif
}

//__________________________________________________________________

void _HYDataPanel::_VerifyInferMenu(void)
{
    MenuHandle  lfMenu          = GetMenuHandle (HY_DATAPANEL_MENU_ID+1),
                inferSubMenu = GetMenuHandle (HY_DATAPANEL_HMENU_ID+8);

    _SimpleList    gp;

    if (GenerateGoodPartitions(gp)) {
        if (!inferSubMenu) {
            inferSubMenu        = NewMenu(HY_DATAPANEL_HMENU_ID+8,"\p");
            InsertMenu (inferSubMenu,hierMenu);
            SetItemCmd (lfMenu,1,hMenuCmd);
            SetItemMark(lfMenu,1,HY_DATAPANEL_HMENU_ID+8);
            InsertMenuItem (inferSubMenu,"\pInfer Topology/L",10000);
            InsertMenuItem (inferSubMenu,"\pInfer Topology with Constraints",10000);
            SetMenuItemText (lfMenu,1,"\pInference...");
        }
    } else {
        if (inferSubMenu) {
            SetItemMark (lfMenu,1,noMark);
            DeleteMenu  (HY_DATAPANEL_HMENU_ID+8);
            DisposeMenu (inferSubMenu);
            SetMenuItemText (lfMenu,1,"\pBuild Function");
            SetItemCmd (lfMenu,1,'L');
        }
    }
}
//__________________________________________________________________

void _HYDataPanel::_SetMenuBar(void)
{
    //BufferToConsole ("_HYDataPanel::_SetMenuBar\n");
    _HYWindow::_SetMenuBar();

    MenuHandle  t  = GetMenuHandle (130),
                dM = GetMenuHandle (HY_DATAPANEL_MENU_ID);
    _HYWindow::_SetMenuBar();
    EnableMenuItem (t,4);
    EnableMenuItem (t,8);

    if (!dM) {
        MenuHandle dataMenu         = NewMenu(HY_DATAPANEL_MENU_ID,"\pData"),
                   lfMenu             = NewMenu (HY_DATAPANEL_MENU_ID+1,"\pLikelihood"),
                   blockMenu       = NewMenu(HY_DATAPANEL_HMENU_ID,"\p"),
                   repeatCharMenu   = NewMenu(HY_DATAPANEL_HMENU_ID+1,"\p"),
                   nameDisplayMenu  = NewMenu(HY_DATAPANEL_HMENU_ID+2,"\p"),
                   omittedSpecies  = NewMenu(HY_DATAPANEL_HMENU_ID+3,"\p"),
                   additionalInfo  = NewMenu(HY_DATAPANEL_HMENU_ID+4,"\p"),
                   lfDisplayMode   = NewMenu(HY_DATAPANEL_HMENU_ID+5,"\p"),
                   simulateData  = NewMenu(HY_DATAPANEL_HMENU_ID+6,"\p"),
                   saveSubMenu     = NewMenu(HY_DATAPANEL_HMENU_ID+7,"\p"),
                   dataProcMenu     = NewMenu(HY_DATAPANEL_HMENU_ID+9,"\p");

        if (!(dataMenu&&blockMenu&&repeatCharMenu&&nameDisplayMenu&&omittedSpecies&&additionalInfo)) {
            warnError (-108);
        }
        InsertMenuItem (dataMenu,"\p(Partition->Selection/1",10000); // 1
        InsertMenuItem (dataMenu,"\p(Selection->Partition/2",10000); // 2
        InsertMenuItem (dataMenu,"\pInvert Selection/3",10000);      // 3
        InsertMenuItem (dataMenu,"\p(-",10000);                      // 4
        InsertMenuItem (dataMenu,"\pBlock Width",10000);             // 5
        InsertMenuItem (dataMenu,"\pRepeating Characters",10000);    // 6
        InsertMenuItem (dataMenu,"\pName Display",10000);            // 7
        InsertMenuItem (dataMenu,"\pOmitted Sequences",10000);       // 8
        InsertMenuItem (dataMenu,"\pAdditional Info",10000);         // 9
        InsertMenuItem (dataMenu,"\p(-",10000);                      // 10
        InsertMenuItem (dataMenu,"\p(Paritition Properties",10000);  // 11
        InsertMenuItem (dataMenu,"\pInput Partition",10000);         // 12
        InsertMenuItem (dataMenu,"\p(-",10000);                      // 13
        InsertMenuItem (dataMenu,"\p(Simulation",10000);             // 14
        InsertMenuItem (dataMenu,"\p(Ancestors",10000);              // 15
        InsertMenuItem (dataMenu,"\p(-",10000);                      // 16
        InsertMenuItem (dataMenu,"\pFont Options",10000);            // 17
        InsertMenuItem (dataMenu,"\p(-",10000);                      // 18
        InsertMenuItem (dataMenu,"\pData Processing",10000);         // 19

        InsertMenuItem (blockMenu,"\p9",10000);
        InsertMenuItem (blockMenu,"\p10",10000);


        InsertMenuItem (lfMenu,"\pBuild Function/L",10000);

        InsertMenuItem (lfMenu,"\p(Display",10000);
        InsertMenuItem (lfMenu,"\p(Optimize/T",10000);
        InsertMenuItem (lfMenu,"\p(-",10000);
        InsertMenuItem (lfMenu,"\p(Show Parameters/H",10000);
        InsertMenuItem (lfMenu,"\p(-",10000);
        InsertMenuItem (lfMenu,"\p(General Bootstrap/B",10000);
        InsertMenuItem (lfMenu,"\p(-",10000);
        InsertMenuItem (lfDisplayMode,"\pLog-Lkhd Only",10000);
        InsertMenuItem (lfDisplayMode,"\pLog-Lkhd & Parameter Values",10000);
        InsertMenuItem (lfDisplayMode,"\pLog-Lkhd & Concise Trees",10000);
        InsertMenuItem (lfDisplayMode,"\pParameter Listing",10000);
        InsertMenuItem (lfDisplayMode,"\pLog-Lkhd & Complete Trees",10000);
        InsertMenuItem (simulateData,"\pSimulate 1",10000);
        InsertMenuItem (simulateData,"\pSimulate 1 To File",10000);
        InsertMenuItem (simulateData,"\pSimulate Many",10000);
        //SetItemMark    (blockMenu,2,checkMark);
        InsertMenuItem (repeatCharMenu,"\pDisplay Actual Character",10000);
        InsertMenuItem (repeatCharMenu,"\pDisplay '.'",10000);
        //SetItemMark    (repeatCharMenu,1,checkMark);
        InsertMenuItem (nameDisplayMenu,"\pNone",10000);
        InsertMenuItem (nameDisplayMenu,"\pFirst 10 characters",10000);
        _HYSequencePane *sp = (_HYSequencePane*)GetCellObject(2,0);
        if (sp->shortHeaderWidth==sp->fullHeaderWidth) {
            DisableMenuItem (nameDisplayMenu,2);
        }
        InsertMenuItem (nameDisplayMenu,"\pFull Names",10000);
        InsertMenuItem (nameDisplayMenu,"\p(-",10000);
        InsertMenuItem (nameDisplayMenu,"\pAlphabetize names",10000);
        InsertMenuItem (nameDisplayMenu,"\pRevert to file order",10000);
        InsertMenuItem (nameDisplayMenu,"\pClean up sequence names",10000);
        //SetItemMark    (nameDisplayMenu,3,checkMark);
        InsertMenuItem (omittedSpecies,"\pRestore All",10000);
        InsertMenuItem (omittedSpecies,"\p(-",10000);
        InsertMenuItem (additionalInfo,"\pConsensus Sequence",10000);
        InsertMenuItem (additionalInfo,"\p(Rate Class",10000);
        InsertMenuItem (additionalInfo,"\p(Aminoacid Translation",10000);
        InsertMenuItem (additionalInfo,"\pReference Sequence",10000);
        InsertMenuItem (saveSubMenu,"\pSave.../S",10000);
        InsertMenuItem (saveSubMenu,"\pSave As...",10000);

        InsertMenu (dataMenu,132);
        InsertMenu (lfMenu,132);
        InsertMenu (blockMenu,hierMenu);
        InsertMenu (repeatCharMenu,hierMenu);
        InsertMenu (nameDisplayMenu,hierMenu);
        InsertMenu (omittedSpecies,hierMenu);
        InsertMenu (additionalInfo,hierMenu);
        InsertMenu (lfDisplayMode,hierMenu);
        InsertMenu (simulateData,hierMenu);
        InsertMenu (saveSubMenu,hierMenu);
        InsertMenu (dataProcMenu,hierMenu);
        SetItemCmd (lfMenu,2,hMenuCmd);
        SetItemMark(lfMenu,2,HY_DATAPANEL_HMENU_ID+5);
        SetItemCmd (dataMenu,5,hMenuCmd);
        SetItemMark(dataMenu,5,HY_DATAPANEL_HMENU_ID);
        SetItemCmd (dataMenu,6,hMenuCmd);
        SetItemMark(dataMenu,6,HY_DATAPANEL_HMENU_ID+1);
        SetItemCmd (dataMenu,7,hMenuCmd);
        SetItemMark(dataMenu,7,HY_DATAPANEL_HMENU_ID+2);
        SetItemCmd (dataMenu,8,hMenuCmd);
        SetItemMark(dataMenu,8,HY_DATAPANEL_HMENU_ID+3);
        SetItemCmd (dataMenu,9,hMenuCmd);
        SetItemMark(dataMenu,9,HY_DATAPANEL_HMENU_ID+4);
        SetItemCmd (dataMenu,14,hMenuCmd);
        SetItemMark(dataMenu,14,HY_DATAPANEL_HMENU_ID+6);
        SetItemCmd (dataMenu,19,hMenuCmd);
        SetItemMark(dataMenu,19,HY_DATAPANEL_HMENU_ID+9);
        if (omittedSeqs.lLength==0) {
            DisableMenuItem (dataMenu,8);
        } else {
            _OmitSelectedSpecies(omittedSeqs);
        }

        if (dataPanelProcessors.lLength == 0) {
            DisableMenuItem (dataMenu,19);
        } else {
            Str255    buffer;
            for (long k=0; k<dataPanelProcessors.lLength; k++) {
                _String *thisItem = (_String*)dataPanelProcessors (k),
                         chopped = thisItem->Cut (thisItem->FindBackwards (':',0,-1)+1,-1);
                StringToStr255  (chopped,buffer);
                InsertMenuItem  (dataProcMenu, buffer,10000);
            }
        }

        CheckMenuItem (additionalInfo,1,addedLines&HY_DATAPANEL_CONSENSUS);
        CheckMenuItem (additionalInfo,2,addedLines&HY_DATAPANEL_RATECLASS);
        CheckMenuItem (additionalInfo,3,addedLines&HY_DATAPANEL_TRANSLATION);
        CheckMenuItem (additionalInfo,4,addedLines&HY_DATAPANEL_REFERENCE);
        CheckMenuItem (blockMenu,(sp->blockWidth==10)?2:1,true);

        if (sp->nameDisplayFlags&HY_SEQUENCE_PANE_NAMES_ALL) {
            CheckMenuItem (nameDisplayMenu,3,true);
        } else if (sp->nameDisplayFlags&HY_SEQUENCE_PANE_NAMES_SHORT) {
            CheckMenuItem (nameDisplayMenu,2,true);
        } else {
            CheckMenuItem (nameDisplayMenu,1,true);
        }
        if (sp->showDots) {
            CheckMenuItem (repeatCharMenu,2,true);
        } else {
            CheckMenuItem (repeatCharMenu,1,true);
        }
        if (dataType&HY_DATAPANEL_NUCDATA) {
            EnableMenuItem (additionalInfo,3);
        }

        _UpdateLFMenu();

        if (aquaInterfaceOn) {
            InsertMenuItem (t,"\p(-", 9);
            InsertMenuItem (t,"\pFind.../F", 10);
            InsertMenuItem (t,"\pSearch and Replace...", 11);
        } else {
            InsertMenuItem (t,"\pFind.../F", 10);
            InsertMenuItem (t,"\pSearch and Replace...", 11);
            InsertMenuItem (t,"\p(-", 12);
        }
    }
    t = GetMenuHandle (129);
    SetItemCmd (t,4,hMenuCmd);
    SetItemMark(t,4,HY_DATAPANEL_HMENU_ID+7);
    _VerifyInferMenu    ();
    InvalMenuBar();
}

//__________________________________________________________________

void _HYDataPanel::_UpdateLFMenu (void)
{
    MenuHandle lfMenu = GetMenuHandle (HY_DATAPANEL_MENU_ID+1),
               dataMenu = GetMenuHandle (HY_DATAPANEL_MENU_ID),
               addMenu= GetMenuHandle (HY_DATAPANEL_HMENU_ID+4);
    if (lfMenu && dataMenu && addMenu) {
        if (lfID>=0) {
            EnableMenuItem (lfMenu,2);
            EnableMenuItem (lfMenu,3);
            EnableMenuItem (lfMenu,5);
            EnableMenuItem (lfMenu,7);
            EnableMenuItem (dataMenu,14);
            EnableMenuItem (dataMenu,15);

            if (((_LikelihoodFunction*)likeFuncList (lfID))->GetCategoryVars().lLength) {
                EnableMenuItem (addMenu,2);
                return;
            }
        } else {
            DisableMenuItem (lfMenu,3);
            DisableMenuItem (lfMenu,2);
            DisableMenuItem (lfMenu,5);
            DisableMenuItem (lfMenu,7);
            DisableMenuItem (dataMenu,14);
            DisableMenuItem (dataMenu,15);
        }
        DisableMenuItem (addMenu,2);
    }
}

//__________________________________________________________________

void _HYDataPanel::_UpdateSelectionChoices (bool toggle)
{
    MenuHandle dataMenu  = GetMenuHandle (HY_DATAPANEL_MENU_ID);

    if (toggle) {
        EnableMenuItem(dataMenu,2);
        //EnableMenuItem(dataMenu,3);
    } else {
        DisableMenuItem(dataMenu,2);
        //DisableMenuItem(dataMenu,3);
    }
    InvalMenuBar();

}

//__________________________________________________________________

void _HYDataPanel::_CopySelectionToClipboard (void)
{
    _HYSequencePane*    sp = (_HYSequencePane*)GetObject(0);
    _String             cbStr (128L,true);

    if (sp->selection.lLength) {
        for (long m=0; m<sp->speciesIndex.lLength; m++) {
            long idx = sp->speciesIndex.lData[m];
            for (long k=0; k<sp->selection.lLength; k++) {
                cbStr << ((_String*)(sp->columnStrings(sp->selection.lData[k])))->sData[idx];
                if (k&&((k+1)%sp->blockWidth==0)) {
                    cbStr << ' ';
                }
            }
            cbStr << '\r';
        }
    } else if (sp->vselection.lLength)
        for (long m=0; m<sp->vselection.lLength; m++) {
            cbStr << (_String*)(sp->rowHeaders(sp->speciesIndex(sp->vselection.lData[m])));
            cbStr << '\r';
        }

    cbStr.Finalize();

    if (cbStr.sLength) {
        PlaceStringInClipboard (cbStr,nil);
    }
}

//__________________________________________________________________

void _HYDataPanel::_OmitSelectedSpecies (_SimpleList& idx)
{
    MenuHandle dataMenu = GetMenuHandle (HY_DATAPANEL_MENU_ID),
               omittedSpecies = GetMenuHandle (HY_DATAPANEL_HMENU_ID+3);
    if (omittedSpecies) {
        _HYSequencePane*    sp = (_HYSequencePane*)GetObject(0);

        for (long k=0; k<idx.lLength; k++) {
            Str255   buffer;
            _String* thisSpec = (_String*)sp->rowHeaders(idx.lData[k]);
            StringToStr255 (*thisSpec, buffer);
            InsertMenuItem (omittedSpecies,buffer,10000);
        }

        EnableMenuItem (dataMenu,8);
    }
}

//__________________________________________________________________

void _HYDataPanel::_RestoreOmittedSequence (long index)
{
    MenuHandle dataMenu = GetMenuHandle (HY_DATAPANEL_MENU_ID),
               omittedSpecies = GetMenuHandle (HY_DATAPANEL_HMENU_ID+3);

    if (index>=0) {
        DeleteMenuItem (omittedSpecies,index+3);
        if (CountMenuItems(omittedSpecies)==2) {
            DisableMenuItem (dataMenu,8);
        }
    } else {
        for (long k=0; k< omittedSeqs.lLength; k++) {
            DeleteMenuItem (omittedSpecies,3);
        }
        DisableMenuItem (dataMenu,8);
    }
}

//__________________________________________________________________

void _HYDataPanel::_UpdatePartitionOperations (_SimpleList* sl)
{
    MenuHandle dataMenu = GetMenuHandle (HY_DATAPANEL_MENU_ID);

    if (sl->lData[0]) {
        EnableMenuItem (dataMenu,1);
        EnableMenuItem (dataMenu,11);
    } else {
        DisableMenuItem(dataMenu,1);
        DisableMenuItem(dataMenu,11);
    }

    InvalMenuBar();
}

//__________________________________________________________________

void _HYDataPanel::_UnsetMenuBar(void)
{
    //BufferToConsole ("_HYDataPanel::_UnsetMenuBar\n");

    MenuHandle treeMenu         = GetMenuHandle (HY_DATAPANEL_MENU_ID),
               lfMenu             = GetMenuHandle (HY_DATAPANEL_MENU_ID+1),
               blockMenu        = GetMenuHandle (HY_DATAPANEL_HMENU_ID),
               repeatCharMenu   = GetMenuHandle (HY_DATAPANEL_HMENU_ID+1),
               nameDisplayMenu  = GetMenuHandle (HY_DATAPANEL_HMENU_ID+2),
               omittedSequences = GetMenuHandle (HY_DATAPANEL_HMENU_ID+3),
               additionalInfo    = GetMenuHandle (HY_DATAPANEL_HMENU_ID+4),
               lfDisplayOptions = GetMenuHandle (HY_DATAPANEL_HMENU_ID+5),
               simulateData  = GetMenuHandle (HY_DATAPANEL_HMENU_ID+6),
               saveSubMenu     = GetMenuHandle (HY_DATAPANEL_HMENU_ID+7),
               inferSubMenu      = GetMenuHandle (HY_DATAPANEL_HMENU_ID+8),
               dataPanelProc    = GetMenuHandle (HY_DATAPANEL_HMENU_ID+9),
               fMenu            = GetMenuHandle (129);

    DeleteMenu (HY_DATAPANEL_MENU_ID);
    DeleteMenu (HY_DATAPANEL_MENU_ID+1);
    DeleteMenu (HY_DATAPANEL_HMENU_ID);
    DeleteMenu (HY_DATAPANEL_HMENU_ID+1);
    DeleteMenu (HY_DATAPANEL_HMENU_ID+2);
    DeleteMenu (HY_DATAPANEL_HMENU_ID+3);
    DeleteMenu (HY_DATAPANEL_HMENU_ID+4);
    DeleteMenu (HY_DATAPANEL_HMENU_ID+5);
    DeleteMenu (HY_DATAPANEL_HMENU_ID+6);
    DeleteMenu (HY_DATAPANEL_HMENU_ID+7);
    DeleteMenu (HY_DATAPANEL_HMENU_ID+9);
    DisposeMenu (treeMenu);
    DisposeMenu (lfMenu);
    DisposeMenu (blockMenu);
    DisposeMenu (repeatCharMenu);
    DisposeMenu (nameDisplayMenu);
    DisposeMenu (omittedSequences);
    DisposeMenu (additionalInfo);
    DisposeMenu (lfDisplayOptions);
    DisposeMenu (simulateData);
    DisposeMenu (saveSubMenu);
    DisposeMenu (dataPanelProc);
    if (inferSubMenu) {
        DeleteMenu (HY_DATAPANEL_HMENU_ID+8);
        DisposeMenu (inferSubMenu);
    }

    SetItemCmd (fMenu,4,'S');
    SetItemMark(fMenu,4,noMark);

    fMenu = GetMenuHandle (130);

    if (!aquaInterfaceOn) {
        DeleteMenuItem (fMenu,13);
    }

    DeleteMenuItem (fMenu,12);
    DeleteMenuItem (fMenu,11);

    if (aquaInterfaceOn) {
        DeleteMenuItem (fMenu,10);
    }


    _HYWindow::_UnsetMenuBar();
}

//__________________________________________________________________
bool _HYDataPanel::_ProcessOSEvent (Ptr vEvent)
{
    EventRecord* theEvent = (EventRecord*)vEvent;
    static  UInt32  lastClick = 0;
    static  int     lastH = 0, lastV = 0;
    if (!_HYTWindow::_ProcessOSEvent (vEvent)) {
        if (theEvent->what==mouseDown) {
            Point localClick = theEvent->where;
            GrafPtr savedPort;
            GetPort(&savedPort);
#ifdef OPAQUE_TOOLBOX_STRUCTS
            SetPort(GetWindowPort(theWindow));
#else
            SetPort(theWindow);
#endif
            GlobalToLocal (&localClick);
            bool    dblClick = (theEvent->when-lastClick<GetDblTime())&&(abs(localClick.h-lastH)<5)&&(abs(localClick.v-lastV)<5);

            lastClick = theEvent->when;
            lastH = localClick.h;
            lastV = localClick.v;
            int   c = FindClickedCell(localClick.h,localClick.v),ch,cv;
            if (c<0) {
                return false;
            }

            if (c==1) { // navBar
                ch = localClick.h-componentL.lData[1]-thermRect.left;
                cv = localClick.v-componentT.lData[1]-thermRect.top;
                if (dblClick) {
                    NavBarDblClick (ch);
                    return true;
                }
                if (navRect.Contains(ch,cv)) {
                    Point   oldPt, newPt,deltaPt;
                    deltaPt.h = (navRect.right+navRect.left)/2-ch;
                    deltaPt.v = (navRect.top+navRect.bottom)/2-cv;
                    oldPt=localClick;
                    if (StillDown()) {
                        while (WaitMouseUp()) {
                            GetMouse( &newPt);
                            //if ( DeltaPoint(oldPt, newPt) )
                            if ( oldPt.h!=newPt.h ) {
                                oldPt=newPt;
                                ch = newPt.h-componentL.lData[1]+deltaPt.h;
                                cv = newPt.v-componentT.lData[1]+deltaPt.v;
                                ch-=thermRect.left;
                                cv-=thermRect.top;
                                forceUpdateForScrolling = true;
                                SetNavRectCenter (ch,cv);
                                forceUpdateForScrolling = false;
                            }
                        }
                    }
                } else {
                    SetNavRectCenter (ch,cv);
                }
                return true;
            } else if ((c==4)&&(theEvent->modifiers&controlKey)) {
                _HYSequencePane* sp2 = (_HYSequencePane*)components (4);
                sp2->ProcessContextualPopUp (localClick.h,localClick.v);
                return true;
            }
        } else if ((theEvent->what==keyDown) || (theEvent->what==autoKey)) {
            unsigned char keyCode = (theEvent->message&keyCodeMask)>>8;
            if ((keyCode==0x7B)||(keyCode==0x7C)) { // left/right arrow
                _HYSequencePane* sp = (_HYSequencePane*) GetObject (0);
                if ((keyCode==0x7B)&&(sp->startColumn)) {
                    sp->HScrollPane (-1);
                } else if ((keyCode==0x7C)&&(sp->endColumn<sp->columnStrings.lLength)) {
                    sp->HScrollPane (1);
                }
                return true;
            }
        }
    } else {
        return true;
    }
    return false;
}



//EOF