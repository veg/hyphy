/*
    Tree  Panel Object  for Win32 API

    Sergei L. Kosakovsky Pond, June 2000-January 2003
*/

#include "HYTreePanel.h"
#include "HYUtils.h"

_String     saveForTreesPrompt ("Save Tree As:");

HPEN        navPen = CreatePen (PS_SOLID,2,RGB(255,151,51));

#define     HY_TREE_WIN32_MENU_BASE  7000

//__________________________________________________________________

bool _HYTreePanel::_ProcessOSEvent (Ptr vEvent)
{
    static  char draggin = 0;
    _HYWindowsUIMessage *theEvent = (_HYWindowsUIMessage*)vEvent;
    if (theEvent->iMsg==WM_RBUTTONDOWN) {
        if (FindClickedCell(LOWORD(theEvent->lParam),HIWORD(theEvent->lParam)) == 0) {
            HandleContextPopup(LOWORD(theEvent->lParam),HIWORD(theEvent->lParam));
            return true;
        }
    }
    if(!_HYTWindow::_ProcessOSEvent (vEvent)) {
        if ((theEvent->iMsg==WM_LBUTTONDOWN)||(theEvent->iMsg==WM_LBUTTONDBLCLK)) {
            int ch = LOWORD(theEvent->lParam),cv = HIWORD(theEvent->lParam),c;
            c = FindClickedCell(ch,cv);
            if (c<0) {
                return false;
            }
            _HYComponent* thisComponent = (_HYComponent*)components(c);
            if (c==1) { // navBar

                ch -= componentL.lData[1];
                cv -= componentT.lData[1];
                if (navRect.Contains(ch,cv)) {
                    draggin = 1;
                } else {
                    SetNavRectCenter (ch,cv);
                }
                return true;
            } else if (c==0) {
                ch -= componentL.lData[0];
                cv -= componentT.lData[0];
                char shiftFlag = 0;

                if (theEvent->wParam&MK_SHIFT) {
                    shiftFlag |= 0x01;
                }
                if (theEvent->wParam&MK_CONTROL) {
                    shiftFlag |= 0x02;
                }

                if (shiftFlag > 2) {
                    draggin = 2;
                    theEvent->iMsg = WM_MOUSEMOVE;
                } else {
                    if (IsVertical()) {
                        c = ch;
                        ch = cv+thisComponent->vOrigin;
                        cv = c+thisComponent->hOrigin;
                    } else {
                        ch+=thisComponent->hOrigin;
                        cv+=thisComponent->vOrigin;
                    }

                    if ((theEvent->iMsg==WM_LBUTTONDBLCLK)&&(currentSelection.lLength)) {
                        InvokeNodeEditor();
                        return true;
                    }

                    if(FindSelection (ch,cv,shiftFlag)) {
                        _UpdateOperationsMenu();
                        RenderTree();
                    }
                }
            }
        }

        if (theEvent->iMsg==WM_LBUTTONUP) {
            draggin = 0;
            return true;
        } else if (theEvent->iMsg==WM_MOUSEMOVE) {
            if (draggin == 1) {
                int ch = LOWORD(theEvent->lParam)-componentL.lData[1],
                    cv = HIWORD(theEvent->lParam)-componentT.lData[1];
                SetNavRectCenter (ch,cv);
            } else {
                if (draggin == 2) {
                    _HYCanvas   *theTree    = (_HYCanvas*)GetObject (0);

                    int ch = LOWORD(theEvent->lParam),
                        cv = HIWORD(theEvent->lParam);

                    if (IsVertical()) {
                        cv = cv-componentL.lData[0]+theTree->hOrigin;
                        ch = ch-componentT.lData[0]+theTree->vOrigin;
                    } else {
                        ch = ch-componentL.lData[0]+theTree->hOrigin;
                        cv = cv-componentT.lData[0]+theTree->vOrigin;
                    }
                    if ((ch<theTree->_HYComponent::GetMaxW())&&(cv<theTree->_HYComponent::GetMaxH())) {
                        FishEyeProjection (ch,theTree->_HYComponent::GetMaxH()-cv,theTree->_HYComponent::GetMaxW(),
                                           theTree->_HYComponent::GetMaxH(),coordTree);
                        treeFlags |= HY_TREEPANEL_PROJECTION;
                        RenderTree(false);
                    }
                }
            }
            return true;
        } else {
            if (theEvent->iMsg==WM_CHAR) {
                TCHAR  theC = (TCHAR)theEvent->wParam;
                //printf ("%x\n", theC);
                if (theC == VK_BACK) {
                    DeleteCurrentSelection();
                    _UpdateOperationsMenu();
                    return true;
                } else if (theC == VK_RETURN) {
                    InvokeNodeEditor ();
                    _UpdateOperationsMenu();
                    return true;
                }
            }
        }
        return false;
    }
}
//__________________________________________________________________


bool        _HYTreePanel::_ProcessMenuSelection (long msel)
{
    if (_HYWindow::_ProcessMenuSelection(msel)) {
        return true;
    }

    switch (msel) {

    case HY_WINDOW_MENU_ID_FILE+1: { // Save tree
        _String           filePath,
                          dtreeName = treeName,
                          ext;

        bool              good = false;

        _List             treeFormats;

        long pidOptions = GetUniversalSaveOptions (treeFormats);

        ext = "Bitmap (.BMP)";
        treeFormats && & ext;
        ext = "Graphical Metafile (.EMF)";
        treeFormats && & ext;

        long              menuChoice = SaveFileWithPopUp (filePath, saveForTreesPrompt,dtreeName,
                                       empty,
                                       treeFormats);
        if (menuChoice>=0) {
            _HYCanvas   *theTree = (_HYCanvas*)GetObject (0);

            if (menuChoice == pidOptions+1) {
                HDC         mfDC = CreateEnhMetaFile (GetDC(NULL),filePath.getStr(),NULL,NULL), saveDC = theTree->thePane;
                if (!mfDC) {
                    _String errMsg ("Could not create meta file device context. Windows Error:");
                    errMsg = errMsg & _String ((long)GetLastError());
                    WarnError (errMsg);
                    return true;
                }
                theTree->thePane = mfDC;
                RenderTree (true, true);
                theTree->thePane = saveDC;
                DeleteEnhMetaFile(CloseEnhMetaFile (mfDC));
            } else if (menuChoice == pidOptions) {
                HBITMAP hBmp = (HBITMAP)GetCurrentObject (theTree->thePane,OBJ_BITMAP);
                BITMAP bmp;
                PBITMAPINFO pbmi;
                long    rowWidth;
                bool    good = true, ruler = (!IsVertical())&&(scaleVariable.sLength);
                if (!::GetObject(hBmp, sizeof(BITMAP), (LPSTR)&bmp)) {
                    good = false;
                } else {

                    pbmi = (PBITMAPINFO) LocalAlloc(LPTR,sizeof(BITMAPINFOHEADER)+sizeof(RGBQUAD) * 2);
                    pbmi->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
                    pbmi->bmiHeader.biWidth = bmp.bmWidth;
                    pbmi->bmiHeader.biHeight = bmp.bmHeight;
                    pbmi->bmiHeader.biPlanes = 1;
                    pbmi->bmiHeader.biBitCount = 1;
                    pbmi->bmiHeader.biClrUsed = 0;
                    pbmi->bmiHeader.biCompression = BI_RGB;
                    pbmi->bmiHeader.biSizeImage = (pbmi->bmiHeader.biWidth + 7) /8;
                    if ( pbmi->bmiHeader.biSizeImage%4) {
                        pbmi->bmiHeader.biSizeImage = (pbmi->bmiHeader.biSizeImage/4+1)*4;
                    }
                    rowWidth = pbmi->bmiHeader.biSizeImage;
                    pbmi->bmiHeader.biSizeImage*=pbmi->bmiHeader.biHeight;
                    pbmi->bmiHeader.biClrImportant = 0;
                    pbmi->bmiHeader.biXPelsPerMeter = 0x0EC4;
                    pbmi->bmiHeader.biYPelsPerMeter = 0x0EC4;
                    RGBQUAD     t = {0,0,0,0};
                    pbmi->bmiColors[0] = t;
                    t.rgbRed = t.rgbGreen = t.rgbBlue = 255;
                    pbmi->bmiColors[1] = t;

                    HANDLE hf;                 // file handle
                    BITMAPFILEHEADER hdr;       // bitmap file-header
                    PBITMAPINFOHEADER pbih;     // bitmap info-header
                    LPBYTE lpBits;
                    // memory pointer
                    DWORD dwTotal;              // total count of bytes
                    DWORD cb;                   // incremental count of bytes
                    BYTE *hp;                   // byte pointer
                    DWORD dwTmp;

                    pbih = (PBITMAPINFOHEADER) pbmi;
                    lpBits = (LPBYTE) GlobalAlloc(GMEM_FIXED, pbih->biSizeImage);

                    if (!lpBits) {
                        good = false;
                    } else if (!GetDIBits(theTree->thePane, hBmp, 0, (WORD) bmp.bmHeight, lpBits, pbmi,DIB_RGB_COLORS)) {
                        good = false;
                        //printf ("\nStage 3 Error:%d\n",GetLastError());
                    } else {
                        hf = CreateFile(filePath.getStr(),
                                        GENERIC_READ | GENERIC_WRITE,
                                        (DWORD) 0,
                                        NULL,
                                        CREATE_ALWAYS,
                                        FILE_ATTRIBUTE_NORMAL,
                                        (HANDLE) NULL);
                        if (hf == INVALID_HANDLE_VALUE) {
                            good = false;
                        } else {
                            hdr.bfType = 0x4d42;
                            if (ruler)
                                hdr.bfSize = (DWORD) (sizeof(BITMAPFILEHEADER) + pbih->biSize + 2 * sizeof(RGBQUAD) +
                                                      (pbih->biSizeImage/pbmi->bmiHeader.biHeight)*(pbmi->bmiHeader.biHeight+HY_TREEPANEL_RULER_EXPANDED));
                            else {
                                hdr.bfSize = (DWORD) (sizeof(BITMAPFILEHEADER) + pbih->biSize + 2 * sizeof(RGBQUAD) + pbih->biSizeImage);
                            }
                            hdr.bfReserved1 = 0;
                            hdr.bfReserved2 = 0;

                            hdr.bfOffBits = (DWORD) (sizeof(BITMAPFILEHEADER) +
                                                     pbih->biSize + 2 * sizeof (RGBQUAD));

                            if (WriteFile(hf, (LPVOID) &hdr, sizeof(BITMAPFILEHEADER),
                                          (LPDWORD) &dwTmp,  NULL)) {
                                DWORD oldImSize = pbih->biSizeImage;
                                if (ruler) {
                                    pbih->biHeight+=HY_TREEPANEL_RULER_EXPANDED;
                                    pbih->biSizeImage=(pbih->biSizeImage/pbmi->bmiHeader.biHeight)*(pbmi->bmiHeader.biHeight+HY_TREEPANEL_RULER_EXPANDED);
                                }
                                if (WriteFile(hf, (LPVOID) pbih, sizeof(BITMAPINFOHEADER)
                                              + 2 * sizeof (RGBQUAD),
                                              (LPDWORD) &dwTmp, ( NULL))) {
                                    dwTotal = cb = oldImSize;
                                    hp = lpBits;
                                    WriteFile(hf, (LPSTR) hp, (int) cb, (LPDWORD) &dwTmp,NULL);

                                }
                            }
                        }
                    }
                    if (ruler) {

                        _HYCanvas   *theRuler = (_HYCanvas*)GetObject (2);
                        hBmp = (HBITMAP)GetCurrentObject (theRuler->thePane,OBJ_BITMAP);
                        if (::GetObject(hBmp, sizeof(BITMAP), (LPSTR)&bmp)) {
                            pbmi->bmiHeader.biHeight = bmp.bmHeight;
                            pbmi->bmiHeader.biSizeImage = (pbmi->bmiHeader.biWidth + 7) /8;
                            if ( pbmi->bmiHeader.biSizeImage%4) {
                                pbmi->bmiHeader.biSizeImage = (pbmi->bmiHeader.biSizeImage/4+1)*4;
                            }
                            pbmi->bmiHeader.biSizeImage *= bmp.bmHeight;
                            if (!GetDIBits(theRuler->thePane, hBmp, 0, (WORD) bmp.bmHeight, lpBits, pbmi,DIB_RGB_COLORS)) {
                                good = false;
                                //printf ("\nStage 3 Error:%d\n",GetLastError());
                            } else {
                                dwTotal = cb = pbmi->bmiHeader.biSizeImage;
                                hp = lpBits;
                                WriteFile(hf, (LPSTR) hp, (int) cb, (LPDWORD) &dwTmp,NULL);
                            }
                        } else {
                            good = false;
                        }
                    }
                    CloseHandle(hf);
                    GlobalFree((HGLOBAL)lpBits);
                    LocalFree (pbmi);
                }
                if (!good) {
                    _String errMsg = "Couldn't create/write .bmp file";
                    WarnError (errMsg);
                }
            } else {
                HandleTreeSave (menuChoice, filePath);

                /*FILE* theFile = fopen (filePath.sData,"w");
                if (theFile)
                {
                    _String res = GetTreeString();
                    fwrite (res.sData,1,res.sLength,theFile);
                    fclose (theFile);
                }
                else
                {
                    filePath = "Error creating tree file";
                    WarnError (filePath);
                }*/
            }
        }
        return true;
    }
    break;
    case HY_WINDOW_MENU_ID_FILE+2: {
        _PrintTree();
        break;
    }

    case HY_WINDOW_MENU_ID_EDIT:   // Undo
        UndoLastOperation();
        break;

    case HY_WINDOW_MENU_ID_EDIT+1: { // COPY
        if (treeFlags&HY_TREEPANEL_CLIPBOARD_READY) {
            CutSelectionToClipboard (false);
        } else {
            _HYCanvas   *theTree = (_HYCanvas*)GetObject (0);
            bool ruler = (!IsVertical())&&(scaleVariable.sLength);
            if (ruler) {
                _HYCanvas   *theRuler = (_HYCanvas*)GetObject (2);

            } else {
                PlaceBitmapInClipboard ((HBITMAP)GetCurrentObject (theTree->thePane,OBJ_BITMAP),theWindow);
            }
        }
        break;
    }

    case HY_WINDOW_MENU_ID_EDIT+2:   // Cut
        CutSelectionToClipboard ();
        break;

    case HY_WINDOW_MENU_ID_EDIT+3:   // Paste
        PasteClipboardTree();
        break;

    case HY_WINDOW_MENU_ID_EDIT+4:   // Delete
        DeleteCurrentSelection();
        break;

    case HY_WINDOW_MENU_ID_EDIT+5:   // Select All
        SelectAllBranches();
        break;

    case HY_WINDOW_MENU_ID_EDIT+6:   // S & R
        HandleSearchAndReplace(false);
        break;

    case HY_WINDOW_MENU_ID_EDIT+7:   // S & R
        HandleSearchAndReplace(true);
        break;

    case HY_TREE_WIN32_MENU_BASE: {
        unsigned short newF;
        if (treeFlags&HY_TREEPANEL_TIP_LABELS) {
            newF = treeFlags - HY_TREEPANEL_TIP_LABELS;
            SetFlags (newF);
        } else {
            newF = treeFlags + HY_TREEPANEL_TIP_LABELS;
            SetFlags (newF);
        }
        CheckMenuItem (GetSubMenu(windowMenu,2),0,MF_BYPOSITION|((treeFlags&HY_TREEPANEL_TIP_LABELS)?MF_CHECKED:MF_UNCHECKED));
        break;
    }

    case HY_TREE_WIN32_MENU_BASE+1: {
        unsigned short newF;
        if (treeFlags&HY_TREEPANEL_INT_LABELS) {
            newF = treeFlags - HY_TREEPANEL_INT_LABELS;
            SetFlags (newF);
        } else {
            newF = treeFlags + HY_TREEPANEL_INT_LABELS;
            SetFlags (newF);
        }
        CheckMenuItem (GetSubMenu(windowMenu,2),1,MF_BYPOSITION|((treeFlags&HY_TREEPANEL_INT_LABELS)?MF_CHECKED:MF_UNCHECKED));
        break;
    }

    case HY_TREE_WIN32_MENU_BASE+2: {
        SwapSelectedSubTrees ();
        break;
    }

    case HY_TREE_WIN32_MENU_BASE+3: {
        CollapseSelectedBranch ();
        break;
    }

    case HY_TREE_WIN32_MENU_BASE+4: {
        JoinSelectedBranches ();
        break;
    }

    case HY_TREE_WIN32_MENU_BASE+5: {
        GraftATip ();
        break;
    }

    case HY_TREE_WIN32_MENU_BASE+6: {
        RerootTree ();
        break;
    }

    case HY_TREE_WIN32_MENU_BASE+7: {
        FlipSelectedBranches ();
        break;
    }

    case HY_TREE_WIN32_MENU_BASE+8:
    case HY_TREE_WIN32_MENU_BASE+9:
    case HY_TREE_WIN32_MENU_BASE+10: {
        HandleSelection (msel-HY_TREE_WIN32_MENU_BASE-8);
        break;
    }

    case HY_TREE_WIN32_MENU_BASE+17:
    case HY_TREE_WIN32_MENU_BASE+18: {
        HandleSelection (msel-HY_TREE_WIN32_MENU_BASE-14);
        break;
    }

    case HY_TREE_WIN32_MENU_BASE+30:
    case HY_TREE_WIN32_MENU_BASE+31:
    case HY_TREE_WIN32_MENU_BASE+32: {
        HandleSelection (msel-HY_TREE_WIN32_MENU_BASE-25);
        break;
    }

    case HY_TREE_WIN32_MENU_BASE+11: {
        InvokeNodeEditor ();
        break;
    }

    case HY_TREE_WIN32_MENU_BASE+12: {
        RecalculateLikelihood ();
        break;
    }

    case HY_TREE_WIN32_MENU_BASE+13: {
        DisplayParameterTable();
        break;
    }

    case HY_TREE_WIN32_MENU_BASE+14: {
        HandleViewOptions ();
        break;
    }

    case HY_TREE_WIN32_MENU_BASE+15:
    case HY_TREE_WIN32_MENU_BASE+16: {
        ShowModelMatrix (msel-HY_TREE_WIN32_MENU_BASE-15);
        break;
    }

    case HY_TREE_WIN32_MENU_BASE+19:
    case HY_TREE_WIN32_MENU_BASE+20: {
        GenerateDistanceTable (msel-HY_TREE_WIN32_MENU_BASE-19);
        break;
    }

    case HY_TREE_WIN32_MENU_BASE+24:
    case HY_TREE_WIN32_MENU_BASE+25: {
        HandleLabels(msel-HY_TREE_WIN32_MENU_BASE-24);
        break;
    }

    case HY_TREE_WIN32_MENU_BASE+35:
    case HY_TREE_WIN32_MENU_BASE+36:
    case HY_TREE_WIN32_MENU_BASE+37:
    case HY_TREE_WIN32_MENU_BASE+38:
    case HY_TREE_WIN32_MENU_BASE+39: {
        HandleComparison (msel-HY_TREE_WIN32_MENU_BASE-35);
        _UpdateOperationsMenu();
        break;
    }

    case HY_TREE_WIN32_MENU_BASE+40:
    case HY_TREE_WIN32_MENU_BASE+41: {
        HandleSelection(msel-HY_TREE_WIN32_MENU_BASE-33);
        break;
    }

    default: { // proc menu
        if (msel>=HY_TREE_WIN32_MENU_BASE+1000) {
            ExecuteProcessor (msel-HY_TREE_WIN32_MENU_BASE-1000);
            _UpdateOperationsMenu();
            return true;
        }
    }
    }

    if (((msel>=HY_TREE_WIN32_MENU_BASE)&&(msel<HY_TREE_WIN32_MENU_BASE+12))||
            ((msel>=HY_WINDOW_MENU_ID_EDIT)&&(msel<HY_WINDOW_MENU_ID_EDIT+6))) {
        _UpdateOperationsMenu();
    }

    DrawMenuBar(theWindow);
    return true;
}



//__________________________________________________________________

void        _HYTreePanel::_PaintNavRect(void)
{
    navRect = ComputeNavRect();
    _HYCanvas* theCanvas = (_HYCanvas*)GetCellObject (0,0);
    HDC      theContext = GetDC(theWindow);
    RECT r;
    HPEN     savePen = (HPEN)SelectObject (theContext,navPen);
    r.left = navRect.left+theCanvas->rel.left+HY_TREEPANEL_NAVSPACING;
    r.right = navRect.right+theCanvas->rel.left;
    r.top = navRect.top+theCanvas->rel.top+HY_TREEPANEL_NAVSPACING;
    r.bottom = navRect.bottom+theCanvas->rel.top;
    MoveToEx (theContext,r.left,r.top,NULL);
    LineTo(theContext,r.right,r.top);
    LineTo(theContext,r.right,r.bottom);
    LineTo(theContext,r.left,r.bottom);
    LineTo(theContext,r.left,r.top);
    SelectObject (theContext,savePen);
    _PaintLFStatus ((Ptr)theContext);
    ReleaseDC (theWindow, theContext);
    //printf ("\nPaintNavRect");
}

//__________________________________________________________________

void        _HYTreePanel::_PrintTree(long hPages, long vPages)
{
    if ((hPages<0)||(vPages<0)) {
        if (!TreePrintSetup (hPages, vPages, (Ptr)this)) {
            return;
        }
    }

    _HYCanvas* theCanvas = (_HYCanvas*)GetObject (0);
    _HYCanvas* theRuler =  (_HYCanvas*)GetObject (2);

    static DOCINFO  di = {sizeof(DOCINFO), "HYPHY.out", NULL };
    static PRINTDLG pd;
    BOOL            SuccessFlag;

    pd.lStructSize         = sizeof(PRINTDLG);
    pd.hwndOwner           = theWindow;
    pd.hDevMode            = NULL;
    pd.hDevNames           = NULL;
    pd.hDC                 = NULL;
    pd.Flags               = PD_ALLPAGES | PD_COLLATE | PD_RETURNDC | PD_NOSELECTION;
    pd.nFromPage           = 1;
    pd.nToPage             = 0xffff;
    pd.nMinPage            = 1;
    pd.nMaxPage            = 0xffff;
    pd.nCopies             = 1;
    pd.hInstance           = NULL;
    pd.lCustData           = 0L;
    pd.lpfnPrintHook       = NULL;
    pd.lpfnSetupHook       = NULL;
    pd.lpPrintTemplateName = NULL;
    pd.lpSetupTemplateName = NULL;
    pd.hPrintTemplate      = NULL;
    pd.hSetupTemplate      = NULL;

    if (!PrintDlg(&pd)) {
        return;
    }

    if (pd.hDC == NULL) {
        pd.hDC = GetPrinterDeviceContext(theWindow);
    }


    EnableWindow(theWindow, FALSE);

    SuccessFlag   = TRUE;
    UserAbortFlag = FALSE;

    PrintDialogHandle = CreateDialog(GetModuleHandle(NULL), (LPCTSTR)"PrintDlgBox", theWindow,
                                     PrintDialogProc);
    SetDlgItemText(PrintDialogHandle, IDD_FNAME, treeName.getStr());

    SetAbortProc(pd.hDC, AbortProc);

    if (StartDoc(pd.hDC, &di) > 0) {

        HDC         saveDC1 = theCanvas->thePane,
                    saveDC2 = theRuler->thePane;


        long    pageW = GetDeviceCaps(pd.hDC, HORZRES),
                pageH = GetDeviceCaps(pd.hDC, VERTRES),
                hRes = GetDeviceCaps(pd.hDC, LOGPIXELSX),
                vRes = GetDeviceCaps(pd.hDC, LOGPIXELSY),
                screenHRes = GetDeviceCaps(saveDC1, LOGPIXELSX),
                screenVRes = GetDeviceCaps(saveDC1, LOGPIXELSY),
                visW = theCanvas->_HYComponent::GetMaxW(),
                visH = theCanvas->_HYComponent::GetMaxH(),
                printW  = pageW * hPages,
                printH  = pageH * vPages;


        bool    hasRuler = (scaleVariable.sLength)&&(!IsVertical());

        _HYFont oldFont  = treeLabelFont;

        _Parameter  hsc = 1.0,
                    vsc = 1.0,
                    xResC = (_Parameter)screenHRes/hRes,
                    yResC = (_Parameter)screenVRes/vRes;

        if (hPages <= 0) {
            hPages = visW/(pageW*xResC)+1;
            printW= pageW * hPages;
        }

        if (vPages <= 0) {
            vPages = visH/(pageH*yResC)+1;
            printH = pageH * vPages;
        }

        hRes = printW*xResC;
        vRes = printH*yResC;

        screenHRes = printW;
        screenVRes = printH;

        printW = hRes;
        printH = vRes;


        if (visW>printW) {
            hsc = (_Parameter)printW/visW;
        }

        if (hasRuler) {
            printH -= HY_TREEPANEL_RULER_EXPANDED;
        }

        if (visH>printH) {
            vsc = (_Parameter)printH/visH;
        }

        _HYFont     bf1 = branchLabel1,
                    bf2 = branchLabel2;

        if ((hsc<1.0)||(vsc<1.0)) {
            treeLabelFont.size = ((_Parameter)treeLabelFont.size)*MIN(hsc,vsc);
            bf1.size = ((_Parameter)bf1.size)*MIN(hsc,vsc);
            bf2.size = ((_Parameter)bf2.size)*MIN(hsc,vsc);
            ShiftScreenCoordinates  (-windowTextMarginH,-windowTextMarginV,coordTree);
            Convert2ScreenCoordinates (hsc,vsc,0,coordTree);
        }

        theCanvas->thePane = pd.hDC;
        theRuler->thePane  = pd.hDC;

        theCanvas->lastPenSize = 2;

        if (visW<printW) {
            visW = (printW-visW)/2;
        } else {
            visW = 0;
        }
        if (visH<printH) {
            visH = (printH-visH)/2;
        } else {
            visH = 0;
        }

        if (IsVertical()) {
            long t = visH;
            visH = visW;
            visW = t;
        }

        for (long hCount = 0; hCount < hPages; hCount ++)
            for (long vCount = 0; vCount < vPages; vCount ++) {
                if (StartPage(pd.hDC) < 0) {
                    SuccessFlag = FALSE;
                } else {
                    theCanvas->SetColor(theCanvas->color);
                    /*if (visH||visW)
                    {
                        ShiftScreenCoordinates  (visW,hasRuler?visH+HY_TREEPANEL_RULER_EXPANDED/2:visH,coordTree);
                    }
                    else
                        if (hasRuler)
                            ShiftScreenCoordinates  (0,HY_TREEPANEL_RULER_EXPANDED/2,coordTree);
                    if (hasRuler)
                    {
                        RenderRuler (hsc,true,visW,visH-HY_TREEPANEL_RULER_EXPANDED/2);
                    }

                    if (visH||visW)
                    {
                        ShiftScreenCoordinates  (-visW,hasRuler?-visH-HY_TREEPANEL_RULER_EXPANDED/2:-visH,coordTree);
                    }
                    else
                        if (hasRuler)
                            ShiftScreenCoordinates  (0,-HY_TREEPANEL_RULER_EXPANDED/2,coordTree);*/


                    SetMapMode  (pd.hDC, MM_ISOTROPIC);
                    SetWindowExtEx (pd.hDC, hRes, vRes,nil);
                    SetViewportExtEx (pd.hDC, screenHRes, screenVRes, nil);

                    //printf ("\n%d %d %d %d %d %d\n", hsc, vsc, hRes, vRes, screenHRes, screenVRes);

                    if (visH||visW) {
                        ShiftScreenCoordinates  (visW,hasRuler?visH+HY_TREEPANEL_RULER_EXPANDED-5:visH,coordTree);
                    } else if (hasRuler&&(vCount==0)) {
                        ShiftScreenCoordinates  (0,HY_TREEPANEL_RULER_EXPANDED-5,coordTree);
                    }

                    if (hCount||vCount) {
                        ShiftScreenCoordinates  (-pageW*xResC*hCount,-pageH*yResC*vCount, coordTree);
                    }


                    if (hasRuler && (vCount == 0)) {
                        RenderRuler (hsc,true,visW,visH);
                    }

                    theCanvas->SetFont (treeLabelFont);

                    if (treeFlags&HY_TREEPANEL_ARCS) {
                        PaintArcs(theCanvas, coordTree);
                    } else if (treeFlags&(HY_TREEPANEL_STRAIGHT|HY_TREEPANEL_CIRCULAR)) {
                        PaintStraight(theCanvas, coordTree);
                    } else {
                        if (IsVertical()) {
                            PaintVSquare(theCanvas, coordTree);
                        } else {
                            PaintSquare(theCanvas, coordTree);
                            if (treeFlags&HY_TREEPANEL_LABEL1) {
                                theCanvas->SetFont(bf1);
                                PaintSquareBranchLabels (theCanvas,coordTree,true);
                            }
                            if (treeFlags&HY_TREEPANEL_LABEL2) {
                                theCanvas->SetFont(bf2);
                                PaintSquareBranchLabels (theCanvas,coordTree,false);
                            }
                        }
                    }

                    if (visH||visW) {
                        ShiftScreenCoordinates  (-visW,hasRuler?-visH-HY_TREEPANEL_RULER_EXPANDED+5:-visH,coordTree);
                    } else if (hasRuler&&(vCount==0)) {
                        ShiftScreenCoordinates  (0,-HY_TREEPANEL_RULER_EXPANDED+5,coordTree);
                    }

                    if (hCount||vCount) {
                        ShiftScreenCoordinates  (pageW*hCount*xResC, pageH*vCount*yResC,coordTree);
                    }
                    if (EndPage (pd.hDC) < 0) {
                        SuccessFlag = FALSE;
                    }
                }
            }
        //if ((hsc<1.0)||(vsc<1.0))
        {
            ShiftScreenCoordinates  (-windowTextMarginH,-windowTextMarginV,coordTree);
            Convert2ScreenCoordinates (1.0/hsc,1.0/vsc,0,coordTree);
        }

        treeLabelFont = oldFont;
        theCanvas->SetFont (oldFont);

        theCanvas->thePane = saveDC1;
        theRuler->thePane = saveDC2;
    } else {
        SuccessFlag = FALSE;
    }


    if (SuccessFlag) {
        EndDoc(pd.hDC);
    }

    if (!UserAbortFlag) {
        EnableWindow(theWindow, TRUE);
        DestroyWindow(PrintDialogHandle);
    }

    DeleteDC (pd.hDC);

    if (!SuccessFlag && !UserAbortFlag) {
        MessageBox(theWindow, "Could not print the tree", AppName, MB_OK | MB_ICONEXCLAMATION);
    }
}

//__________________________________________________________________

void _HYTreePanel::_SetMenuBar(void)
{
    _HYWindow::_SetMenuBar();

    HMENU            windowMenu =  GetMenu (theWindow),
                     editMenu   =  GetSubMenu(windowMenu,1),
                     treeMenu   =  GetSubMenu(windowMenu,2);

    if (!treeMenu) {
        treeMenu = CreateMenu();

        HMENU          labelMenu  =  CreatePopupMenu(),
                       selectMenu =  CreatePopupMenu(),
                       compMenu   =  CreatePopupMenu(),
                       procMenu   =  CreatePopupMenu();

        checkPointer   (treeMenu);
        checkPointer   (labelMenu);
        checkPointer   (selectMenu);
        checkPointer   (compMenu);
        checkPointer   (procMenu);

        InsertMenu      (labelMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_TREE_WIN32_MENU_BASE+24,
                         "&Above Branches\tCtrl-8");

        InsertMenu      (labelMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_TREE_WIN32_MENU_BASE+25,
                         "&Below Branches\tCtrl-9");

        InsertMenu      (compMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_TREE_WIN32_MENU_BASE+35,
                         "Test For &Equality");

        InsertMenu      (compMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING|MF_GRAYED, HY_TREE_WIN32_MENU_BASE+36,
                         "Find &Subtree In Another Tree");

        InsertMenu      (compMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_TREE_WIN32_MENU_BASE+37,
                         "Find &Maximal Common Subtree");

        InsertMenu      (compMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_TREE_WIN32_MENU_BASE+38,
                         "Find Maximal Common &Forest");


        InsertMenu      (compMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_TREE_WIN32_MENU_BASE+39,
                         "&Match To Tree Pattern");


        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING|((treeFlags&HY_TREEPANEL_TIP_LABELS)?
                         MF_CHECKED:0), HY_TREE_WIN32_MENU_BASE, "Tip &Labels");             //0

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING|((treeFlags&HY_TREEPANEL_INT_LABELS)?
                         MF_CHECKED:0), HY_TREE_WIN32_MENU_BASE+1, "I&nternal Labels"); //1

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_SEPARATOR, 0, nil);                            //2

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING|MF_GRAYED, HY_TREE_WIN32_MENU_BASE+2,   //3
                         "S&wap Subtrees\tCtrl-1");

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING|MF_GRAYED, HY_TREE_WIN32_MENU_BASE+3,   //4
                         "&Collapse Branch\tCtrl-2");

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING|MF_GRAYED, HY_TREE_WIN32_MENU_BASE+4,   //5
                         "&Join\tCtrl-3");

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING|MF_GRAYED, HY_TREE_WIN32_MENU_BASE+5,   //6
                         "&Graft A Tip\tCtrl-4");

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING|MF_GRAYED, HY_TREE_WIN32_MENU_BASE+6,   //7
                         "&Reroot\tCtrl-5");

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING|MF_GRAYED, HY_TREE_WIN32_MENU_BASE+7,
                         "&Flip tip ordering\tCtrl-6"); // 8

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_SEPARATOR, 0, nil); // 9


        InsertMenu      (selectMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING|MF_GRAYED, HY_TREE_WIN32_MENU_BASE+8,
                         "Select &Entire Subtree");
        InsertMenu      (selectMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_TREE_WIN32_MENU_BASE+9,
                         "Select Incom&plete Branches\tCtrl-I");
        InsertMenu      (selectMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_TREE_WIN32_MENU_BASE+10,
                         "Select Bran&ches Without Models");
        InsertMenu      (selectMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_TREE_WIN32_MENU_BASE+17,
                         "Select Branches By &Name");

        InsertMenu      (selectMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_TREE_WIN32_MENU_BASE+18,
                         "Select Branches By &Length");

        InsertMenu      (selectMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_TREE_WIN32_MENU_BASE+30,
                         "In&vert Selection");

        InsertMenu      (selectMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING|MF_GRAYED, HY_TREE_WIN32_MENU_BASE+31,
                         "&Grow Selection");

        InsertMenu      (selectMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING|MF_GRAYED, HY_TREE_WIN32_MENU_BASE+32,
                         "&Map Selection to Datapanel");

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_POPUP, (UINT)selectMenu,
                         "Select Branches"); // 10

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING|MF_GRAYED, HY_TREE_WIN32_MENU_BASE+11,
                         "Edit Prope&rties"); // 11

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_SEPARATOR, 0, nil); // 12

        UINT            flags = MF_BYPOSITION|MF_STRING;

        if (likeFuncID < 0) {
            flags |= MF_GRAYED;
        }

        InsertMenu      (treeMenu,  0xFFFFFFFF, flags, HY_TREE_WIN32_MENU_BASE+12,
                         "Optimi&ze Again\tCtrl-T"); // 13

        InsertMenu      (treeMenu,  0xFFFFFFFF, flags, HY_TREE_WIN32_MENU_BASE+13,
                         "Sh&ow Parameters in Table\tCtrl-H"); //14

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_SEPARATOR, 0, nil); // 15

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_TREE_WIN32_MENU_BASE+14,
                         "Tr&ee Display Options..."); // 16

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_POPUP, (UINT)labelMenu,
                         "Branch Labels"); // 17

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_SEPARATOR, 0, nil); // 18

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING|MF_GRAYED, HY_TREE_WIN32_MENU_BASE+15, // 19
                         "Show Rate Matri&x");
        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING|MF_GRAYED, HY_TREE_WIN32_MENU_BASE+16,
                         "S&how Transition Matrix"); // 20


        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_SEPARATOR, 0, nil); // 21

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_TREE_WIN32_MENU_BASE+19, // 22
                         "Pairwise Distan&ces");

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_TREE_WIN32_MENU_BASE+20, // 23
                         "Branch Length Distributi&on");

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_SEPARATOR, 0, nil); // 24

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_POPUP, (UINT)compMenu,
                         "Tree Compar&ison"); // 25

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_TREE_WIN32_MENU_BASE+40, // 26
                         "Match Leaves To Se&quence Data");

        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_TREE_WIN32_MENU_BASE+41, // 27
                         "Find selection in anothe&r Tree");


        InsertMenu      (treeMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_SEPARATOR, 0, nil); // 12

        if (treeProcessors.lLength == 0) {
            InsertMenu      (treeMenu, 0xFFFFFFFF, MF_BYPOSITION|MF_STRING|MF_GRAYED, 0, "&Additional Tools");
            DestroyMenu     (procMenu);
        } else {
            for (long k=0; k<treeProcessors.lLength; k++) {
                _String *thisItem = (_String*)treeProcessors (k),
                         chopped = thisItem->Cut (thisItem->FindBackwards ('\\',0,-1)+1,-1);

                InsertMenu      (procMenu, 0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_TREE_WIN32_MENU_BASE+1000+k, chopped.sData);
            }
            InsertMenu      (treeMenu, 0xFFFFFFFF, MF_BYPOSITION|MF_POPUP, (UINT)procMenu, "&Additional Tools");
        }


        InsertMenu      (editMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_SEPARATOR, 0, nil); // 24

        InsertMenu      (editMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_WINDOW_MENU_ID_EDIT+6,
                         "Search and Repla&ce\tCtrl-F"); // 25

        InsertMenu      (editMenu,  0xFFFFFFFF, MF_BYPOSITION|MF_STRING, HY_WINDOW_MENU_ID_EDIT+7, // 23
                         "Search and Re&place in Selection");


        accels       << (FCONTROL|FVIRTKEY);
        accels       << '1';
        accels       << HY_TREE_WIN32_MENU_BASE+2;

        accels       << (FCONTROL|FVIRTKEY);
        accels       << '2';
        accels       << HY_TREE_WIN32_MENU_BASE+3;

        accels       << (FCONTROL|FVIRTKEY);
        accels       << '3';
        accels       << HY_TREE_WIN32_MENU_BASE+4;

        accels       << (FCONTROL|FVIRTKEY);
        accels       << '4';
        accels       << HY_TREE_WIN32_MENU_BASE+5;

        accels       << (FCONTROL|FVIRTKEY);
        accels       << '5';
        accels       << HY_TREE_WIN32_MENU_BASE+6;

        accels       << (FCONTROL|FVIRTKEY);
        accels       << '6';
        accels       << HY_TREE_WIN32_MENU_BASE+7;

        accels       << (FCONTROL|FVIRTKEY);
        accels       << 'I';
        accels       << HY_TREE_WIN32_MENU_BASE+9;

        accels       << (FCONTROL|FVIRTKEY);
        accels       << 'T';
        accels       << HY_TREE_WIN32_MENU_BASE+12;

        accels       << (FCONTROL|FVIRTKEY);
        accels       << 'H';
        accels       << HY_TREE_WIN32_MENU_BASE+13;

        accels       << (FCONTROL|FVIRTKEY);
        accels       << '8';
        accels       << HY_TREE_WIN32_MENU_BASE+24;

        accels       << (FCONTROL|FVIRTKEY);
        accels       << '9';
        accels       << HY_TREE_WIN32_MENU_BASE+25;

        accels       << (FCONTROL|FVIRTKEY);
        accels       << 'F';
        accels       << HY_WINDOW_MENU_ID_EDIT+6;

        InsertMenu   (windowMenu, 2, MF_BYPOSITION|MF_POPUP, (UINT) treeMenu , "&Tree");

        _AddStandardAccels();
        _BuildAccelTable  (true);
        accels.Clear();


        treeMenu =  GetSubMenu(windowMenu,1);
        EnableMenuItem (treeMenu, 2, MF_ENABLED|MF_BYPOSITION);
        EnableMenuItem (treeMenu, 6, MF_ENABLED|MF_BYPOSITION);
    }
}

//__________________________________________________________________

void _HYTreePanel::_UnsetMenuBar(void)
{

}

//__________________________________________________________________

void _HYTreePanel::_UpdateOperationsMenu (void)
{
    node<nodeCoord>* node1, *node2, *t;
    HMENU       treeMenu    = GetSubMenu(GetMenu (theWindow),2),
                editMenu  = GetSubMenu(GetMenu (theWindow),1),
                selectMenu  = GetSubMenu(treeMenu, 10),
                compMenu   = GetSubMenu(treeMenu, 25);

    EnableMenuItem(treeMenu, 3 ,MF_GRAYED|MF_BYPOSITION);
    EnableMenuItem(treeMenu, 4 ,MF_GRAYED|MF_BYPOSITION);
    EnableMenuItem(treeMenu, 5 ,MF_GRAYED|MF_BYPOSITION);
    EnableMenuItem(treeMenu, 6 ,MF_GRAYED|MF_BYPOSITION);
    EnableMenuItem(treeMenu, 7 ,MF_GRAYED|MF_BYPOSITION);
    EnableMenuItem(treeMenu, 8 ,MF_GRAYED|MF_BYPOSITION);
    EnableMenuItem(selectMenu,0 ,MF_GRAYED|MF_BYPOSITION);
    EnableMenuItem(selectMenu,6 ,MF_GRAYED|MF_BYPOSITION);
    EnableMenuItem(selectMenu,7 ,MF_GRAYED|MF_BYPOSITION);
    EnableMenuItem(selectMenu,8 ,MF_GRAYED|MF_BYPOSITION);
    EnableMenuItem(treeMenu,11 ,MF_GRAYED|MF_BYPOSITION);
    EnableMenuItem(treeMenu,19 ,MF_GRAYED|MF_BYPOSITION);
    EnableMenuItem(treeMenu,20 ,MF_GRAYED|MF_BYPOSITION);

    EnableMenuItem(editMenu, 8 ,MF_GRAYED|MF_BYPOSITION);
    EnableMenuItem(editMenu, 11,MF_GRAYED|MF_BYPOSITION);
    EnableMenuItem(compMenu, 1 ,MF_GRAYED|MF_BYPOSITION);

    bool  good = true;
    long  k,j;
    if (currentSelection.lLength==2) {
        node1 = (node<nodeCoord>*)currentSelection(0);
        node2 = (node<nodeCoord>*)currentSelection(1);
        t = node1->parent;
        while (t) {
            if (t==node2) {
                good = false;
                break;
            }
            t = t->parent;
        }
        if (good) {
            t = node2->parent;
            while (t) {
                if (t==node1) {
                    good = false;
                    break;
                }
                t = t->parent;
            }
        }
        if (good) {
            EnableMenuItem(treeMenu, 3 ,MF_ENABLED|MF_BYPOSITION);
        }
    }
    if (currentSelection.lLength) {
        EnableMenuItem(treeMenu, 6 ,MF_ENABLED|MF_BYPOSITION);
        EnableMenuItem(treeMenu, 11 ,MF_ENABLED|MF_BYPOSITION);
        EnableMenuItem(selectMenu,6 ,MF_ENABLED|MF_BYPOSITION);
        EnableMenuItem(selectMenu,7 ,MF_ENABLED|MF_BYPOSITION);
        EnableMenuItem(selectMenu,8 ,MF_ENABLED|MF_BYPOSITION);
        EnableMenuItem(editMenu, 11,MF_ENABLED|MF_BYPOSITION);
        for (k=0; k<currentSelection.lLength; k++) {
            node1 = (node<nodeCoord>*)currentSelection(k);
            if (node1->get_num_nodes()) {
                EnableMenuItem(treeMenu, 4 ,MF_ENABLED|MF_BYPOSITION);
                break;
            }
        }
        for (k=0; k<currentSelection.lLength; k++) {
            node1 = (node<nodeCoord>*)currentSelection(k);
            if (!node1->get_num_nodes()) {
                EnableMenuItem(treeMenu, 5 ,MF_ENABLED|MF_BYPOSITION);
                break;
            }
        }
        if (currentSelection.lLength>=2) {
            node1 = (node<nodeCoord>*)currentSelection(0);
            t = node1->parent;
            if (t&&(t->get_num_nodes()>currentSelection.lLength)) {
                for (k=1; k<currentSelection.lLength; k++) {
                    node1 = (node<nodeCoord>*)currentSelection(k);
                    if (node1->parent!=t) {
                        break;
                    }
                }
                if (k==currentSelection.lLength) {
                    EnableMenuItem(treeMenu, 5 ,MF_ENABLED|MF_BYPOSITION);
                }
            }
        } else {
            node1 = (node<nodeCoord>*)currentSelection(0);
            if (node1->parent) {
                EnableMenuItem(treeMenu, 7 ,MF_ENABLED|MF_BYPOSITION);
                ModifyMenu    (treeMenu, 7, MF_BYPOSITION, HY_TREE_WIN32_MENU_BASE+6, "Reroot");
            }
            if (node1->get_num_nodes()>0) {
                EnableMenuItem(treeMenu, 8 ,MF_ENABLED|MF_BYPOSITION);
                EnableMenuItem(selectMenu,0 ,MF_ENABLED|MF_BYPOSITION);
            }
            if (node1->in_object.varRef>=0) {
                _CalcNode* thisCNode = (_CalcNode*)LocateVar(node1->in_object.varRef);

                if (thisCNode&&(thisCNode->GetModelIndex()>=0)) {
                    EnableMenuItem(treeMenu,19 ,MF_ENABLED|MF_BYPOSITION);
                    EnableMenuItem(treeMenu,20 ,MF_ENABLED|MF_BYPOSITION);
                }
            }
        }
    } else {
        _TheTree *me = LocateMyTreeVariable();
        if (me) {
            if (me->RootedFlag()==UNROOTED) {
                EnableMenuItem(treeMenu, 7 ,MF_ENABLED|MF_BYPOSITION);
                ModifyMenu    (treeMenu, 7, MF_BYPOSITION,HY_TREE_WIN32_MENU_BASE+6,"Balance");
            } else {
                EnableMenuItem(treeMenu, 7 ,MF_ENABLED|MF_BYPOSITION);
                ModifyMenu    (treeMenu, 7, MF_BYPOSITION, HY_TREE_WIN32_MENU_BASE+6,"Unroot");
            }
        }
    }

    if (likeFuncID!=-1) {
        EnableMenuItem(treeMenu, 13 ,MF_ENABLED|MF_BYPOSITION);
        EnableMenuItem(treeMenu, 14 ,MF_ENABLED|MF_BYPOSITION);
    } else {
        EnableMenuItem(treeMenu, 13 ,MF_GRAYED|MF_BYPOSITION);
        EnableMenuItem(treeMenu, 14 ,MF_GRAYED|MF_BYPOSITION);
    }

    if (treePanelClipboardRoot&&(currentSelection.lLength==1)) {
        EnableMenuItem(editMenu, 4 ,MF_ENABLED|MF_BYPOSITION);
    } else {
        EnableMenuItem(editMenu, 4 ,MF_GRAYED|MF_BYPOSITION);
    }

    EnableMenuItem(editMenu, 3 ,MF_GRAYED|MF_BYPOSITION);
    // check if can cut/paste
    t = nil;

    for (k=0; k<currentSelection.lLength; k++) {
        node1 = (node<nodeCoord>*)currentSelection.lData[k];
        if (node1->parent) {
            if (currentSelection.Find((long)node1->parent)<0) {
                if (t) {
                    break;
                } else {
                    t = node1;
                }
            }
            for (j=0; j<node1->nodes.length; j++) {
                if (currentSelection.Find((long)node1->nodes.data[j])<0) {
                    break;
                }
            }
            if (j<node1->nodes.length) {
                break;
            }
        } else {
            if (t) {
                break;
            } else {
                t = node1;
            }
        }
    }
    selectionTop = nil;
    treeFlags &= 0xFF7F;
    if (t&&(t->parent!=coordTree)&&(t->parent)) {
        if (k==currentSelection.lLength) {
            EnableMenuItem(editMenu, 3 ,MF_ENABLED|MF_BYPOSITION);
            EnableMenuItem(compMenu, 1 ,MF_ENABLED|MF_BYPOSITION);
            selectionTop = t;
            treeFlags |= HY_TREEPANEL_CLIPBOARD_READY;
        }
    }
    _String undoMessage;
    EnableMenuItem(editMenu, 0 ,MF_ENABLED|MF_BYPOSITION);
    UINT   menuFlag = MF_BYPOSITION;
    switch (undoCode) {
    case 1:
        undoMessage = "Undo Swap";
        break;
    case 2:
        undoMessage = "Undo Flip";
        break;
    case 3:
        undoMessage = "Undo Collapse";
        break;
    case 4:
        undoMessage = "Undo Delete";
        break;
    case 5:
        undoMessage = "Undo Join";
        break;
    case 6:
        undoMessage = "Undo Cut";
        break;
    case 7:
        undoMessage = "Undo Graft";
        break;
    case 8:
        undoMessage = "Undo Paste";
        break;
    case 9:
        undoMessage = "Undo Subtree Move";
        break;
    default:
        undoMessage = "Can't Undo";
        menuFlag |= MF_GRAYED;
    }
    undoMessage = undoMessage & "\tCtrl-Z";
    ModifyMenu (editMenu, 0, menuFlag, HY_WINDOW_MENU_ID_EDIT,undoMessage.sData);
    //printf ("Undo message = %s\n", undoMessage.sData);
    DrawMenuBar(theWindow);
}

//__________________________________________________________________

void _HYTreePanel::_HandleIdleEvent (void)
{
    /*#ifdef TARGET_API_MAC_CARBON
    Point    curMouse;
    GetGlobalMouse (&curMouse);

    unsigned long t;
    GetDateTime(&t);


    if ((abs(curMouse.h-saveMouseH)<=3)
      &&(abs(curMouse.v-saveMouseV)<=3)
      &&(t-lastSave>.5))

    {
        if (!HasToolTip())
        {
            GrafPtr curPort;
            GetPort (&curPort);
            SetPort (GetWindowPort (theWindow));
            _DisplayBranchFloater();
            SetPort (curPort);
        }

        lastSave   = t;
    }

    saveMouseH = curMouse.h;
    saveMouseV = curMouse.v;
    #endif*/
}

//__________________________________________________________________

void        _HYTreePanel::_PaintLFStatus(Ptr p)
{
    HDC theDC = (HDC)p;

    if (!theDC) {
        theDC = GetDC (theWindow);
    }

    if (likeFuncID<0) {
        _PaintTheCircle (redButtonIcon,theWindow, theDC);
    } else {
        if (dubiousNodes.lLength) {
            _PaintTheCircle (yellowButtonIcon,theWindow, theDC);
        } else {
            _PaintTheCircle (greenButtonIcon,theWindow, theDC);
        }
    }
    ReleaseDC (theWindow, theDC);
}

//EOF