#include "HYUtils.h"
#include "HYGraphicPane.h"
#include "batchlan.h"
#include "windows.h"
#include "HYWindow.h"
#include "time.h"
#include "shlobj.h"
#include "HYObjectInspector.h"
#include "HYConsoleWindow.h"
#include <Wininet.h>

extern  _String   lastWinPathUsed,
        dialogPrompt,
        objectInspectorTitle ;

extern  _List     pathNames;

extern  bool      updateTimer;

extern  clock_t   timerStart,
        lastTimer;

long    lastFileTypeSelection = 3;


//__________________________________________________________________________________

HMENU   ListToPopUpMenu (_List& menuOptions, long base)
{
    HMENU        listMenu = CreatePopupMenu ();
    checkPointer (listMenu);

    for (long counter=0; counter<menuOptions.lLength; counter++) {
        _String *postItem = (_String*)(menuOptions(counter));

        if (*postItem==_String("SEPARATOR")) {
            AppendMenu (listMenu, MF_SEPARATOR, 0, nil);
        } else {
            AppendMenu (listMenu, MF_STRING, base+counter+1, postItem->sData);
        }
    }

    return listMenu;
}

//__________________________________________________________________________________
void    PlaceStringInClipboard (_String& res,Ptr powner)
{
    HWND    owner = (HWND)powner;
    HGLOBAL TextHandle = GlobalAlloc(GHND, res.sLength+1);
    if (TextHandle) {
        char *TextPtr = (char*)GlobalLock(TextHandle);
        char *SourcePtr = res.sData;
        for (; *SourcePtr; SourcePtr++) {
            *TextPtr++ = *SourcePtr;
        }
        *TextPtr = '\0';
        GlobalUnlock(TextHandle);
        OpenClipboard(owner);
        EmptyClipboard();
        SetClipboardData(CF_TEXT, TextHandle);
        CloseClipboard();
    }
}

//__________________________________________________________________________________
void    PlaceBitmapInClipboard (HBITMAP res,HWND owner)
{
    if (res) {
        OpenClipboard(owner);
        EmptyClipboard();
        SetClipboardData(CF_BITMAP, res);
        CloseClipboard();
    }
}

//__________________________________________________________________________________

_String HandlePullDown (_List& menuOptions, long l, long t,long startPos)
{
    if (menuOptions.lLength) {
        HMENU daMenu = ListToPopUpMenu (menuOptions);

        SetMenuDefaultItem (daMenu, startPos-1, true);

        long tSel = TrackPopupMenu (daMenu, TPM_LEFTALIGN|TPM_TOPALIGN|TPM_RETURNCMD|TPM_LEFTBUTTON|TPM_NONOTIFY,
                                    l,
                                    t,
                                    0, (HWND)hyphyConsoleWindow->GetOSWindowData(), NULL);

        if (tSel) {
            DestroyMenu (daMenu);
            return     *((_String*)menuOptions(tSel-1));
        } else {
            DestroyMenu (daMenu);
        }
    }

    return empty;
}

//__________________________________________________________________________________

long HandlePullDownWithFont (_List& menuOptions, long l, long t,long startPos,_String,long)
{
    _String res = HandlePullDown (menuOptions, l, t, startPos);

    return menuOptions.Find (&res);
}
//________________________________________________________
Ptr     ProcureIconResource (long iconID)
{
    //return (Ptr)LoadImage  (GetModuleHandle(NULL),MAKEINTRESOURCE(iconID),IMAGE_BITMAP,
    //0,0,
    //LR_DEFAULTSIZE);
    return (Ptr)LoadBitmap (GetModuleHandle(NULL),MAKEINTRESOURCE(iconID));
}


//________________________________________________________
long    GetVisibleStringWidth (_String& s, _HYFont& f)
{
    /*static HWND DumbWindow =
                        CreateWindow("STATIC",
                        "",
                        SS_SIMPLE ,
                        CW_USEDEFAULT,
                        CW_USEDEFAULT,
                        10,
                        10,
                        NULL,
                        NULL,
                        GetModuleHandle(NULL),
                        NULL);        */
    static _HYFont lastFont;
    //HDC   defDC = GetDC(DumbWindow);
    if ((lastFont.size!=f.size)||(lastFont.style!=f.style)||(lastFont.face!=f.face)) {
        HFONT   newFont = CreateFont (f.size,0,0,0,(f.style&HY_FONT_BOLD)?FW_BOLD:FW_NORMAL,f.style&HY_FONT_ITALIC,FALSE,FALSE,ANSI_CHARSET,OUT_DEFAULT_PRECIS,
                                      CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,DEFAULT_PITCH|FF_DONTCARE,f.face.sData);
        DeleteObject(SelectObject (otherDC,newFont));
        lastFont.face = f.face;
        lastFont.size = f.size;
        lastFont.style = f.style;
    }
    SIZE  textSize;
    GetTextExtentPoint32(otherDC,s.sData,s.sLength,&textSize);
    //ReleaseDC (DumbWindow, defDC);
    return textSize.cx;
}

//________________________________________________________

long    SaveFileFunction (_String& fileName, _String& prompt, _String &startName, char* fileFilters, Ptr WindowHNDL)
{
    static char buffer[2049], buffer2[256];
    OPENFILENAME ofn;
    ofn.lStructSize = sizeof (OPENFILENAME);
    ofn.hwndOwner = (HWND)WindowHNDL;

    ofn.lpstrFilter = fileFilters;
    ofn.lpstrCustomFilter = nil;

    ofn.nFilterIndex = 1;

    buffer[0] = 0;
    buffer2[0] = 0;
    ofn.lpstrFile = buffer;
    ofn.nMaxFile = 2048;
    ofn.lpstrFileTitle = buffer2;
    memcpy (buffer,startName.sData,startName.sLength+1);
    ofn.nMaxFileTitle = 255;
    ofn.lpstrInitialDir = nil;
    ofn.lpstrDefExt = nil;
    if (lastWinPathUsed.sLength) {
        ofn.lpstrInitialDir = lastWinPathUsed.getStr();
    } else {
        ofn.lpstrInitialDir = ((_String*)pathNames(0))->getStr();
    }
    ofn.lpstrTitle = prompt.getStr();
    ofn.Flags           = OFN_HIDEREADONLY|OFN_PATHMUSTEXIST|OFN_OVERWRITEPROMPT|OFN_EXPLORER;

    if (!GetSaveFileName (&ofn)) {
        return -1;
    }

    fileName = buffer;
    lastWinPathUsed = fileName.Cut (0,ofn.nFileOffset-1);

    return     ofn.nFilterIndex-1;
}

//_________________________________________________________________________
_HYColor    GetDialogBackgroundColor (void)
{
    COLORREF    backColor = GetSysColor (COLOR_MENU);
    _HYColor    res = {GetRValue (backColor), GetGValue (backColor), GetBValue (backColor)};
    return res;
}

//________________________________________________________
_HYRect     GetScreenDimensions (void)
{
    /*RECT          dskRect;
    GetWindowRect   (GetDesktopWindow(), &dskRect);
    _HYRect         res = {0,0,0,0,0};

    res.right       = dskRect.right-dskRect.left;
    res.bottom      = dskRect.bottom-dskRect.top;

    return res;*/

    HMONITOR        defMonitor = MonitorFromPoint ((POINT) {
        0,0
    },MONITOR_DEFAULTTOPRIMARY);
    MONITORINFO     mInfo;
    mInfo.cbSize = sizeof (MONITORINFO);
    GetMonitorInfo  (defMonitor,&mInfo);

    _HYRect         res = {0,0,0,0,0};

    res.right       = mInfo.rcWork.right-mInfo.rcWork.left+1;
    res.bottom      = mInfo.rcWork.bottom-mInfo.rcWork.top-50;

    /*char buffer [1024];

    sprintf (buffer, "%d %d %d %d; %d %d %d %d\n",
                     mInfo.rcWork.left, mInfo.rcWork.top, mInfo.rcWork.right, mInfo.rcWork.bottom,
                     mInfo.rcMonitor.left, mInfo.rcMonitor.top, mInfo.rcMonitor.right, mInfo.rcMonitor.bottom );
    BufferToConsole (buffer);*/

    return res;
}

//________________________________________________________
void        CenterWindow (_HYGuiObject* g)
{
    _HYWindow* w = (_HYWindow*)g;

    _HYRect   screen = GetScreenDimensions();

    long      cleft = 0, ctop = 0;

    if (screen.right>w->right) {
        cleft = (screen.right-w->right)/2;
    }
    if (screen.bottom>w->bottom) {
        ctop = (screen.bottom-w->bottom)/2;
    }

    w->_SetPosition (cleft,ctop);
}

//_________________________________________________________________________
void    DelayNMs (long ms)
{
    Sleep (ms);
}

//________________________________________________________

long    SaveFileWithPopUp (_String& fileName, _String& prompt, _String& defFileName,
                           _String& , _List& menuOptions)

// TBI - add proper dialog window owner

{
    static char buffer[2049],
           buffer2[256];

    _String     formatString (128L, true);

    long        k;

    _List       extensions;

    for (k=0; k<menuOptions.lLength; k++) {
        _String* menuItem = (_String*)menuOptions(k);
        formatString << menuItem;
        formatString << '\0';
        if ((menuItem->sLength>5)&&(menuItem->sData[menuItem->sLength-1]==')')
                &&(menuItem->sData[menuItem->sLength-6]=='(')) {
            formatString << '*';
            formatString << menuItem->Cut (menuItem->sLength-5,menuItem->sLength-2);
            extensions.AppendNewInstance (new _String(menuItem,menuItem->sLength-5,menuItem->sLength-2));
        } else {
            formatString << "*.*";
            extensions && & empty;
        }
        formatString << '\0';
    }
    formatString.Finalize();

    long res = SaveFileFunction (fileName, prompt, defFileName, formatString.sData, (Ptr)GetForegroundWindow());

    if (res>=0) {
        _String * ext = (_String*)extensions (res);
        if (ext->sLength) {
            if (fileName.FindAnyCase (*ext, fileName.sLength-ext->sLength, -1)<0) {
                fileName = fileName & *ext;
            }
        }
    }

    return res;

    /*OPENFILENAME ofn;
    ofn.lStructSize = sizeof (OPENFILENAME);
    ofn.hwndOwner   = nil;

    ofn.lpstrFilter = formatString.sData;
    ofn.lpstrCustomFilter = nil;

    ofn.nFilterIndex = 0;

    buffer[0]           = 0;
    buffer2[0]          = 0;
    ofn.lpstrFile       = buffer;
    ofn.nMaxFile        = 2048;
    ofn.lpstrFileTitle  = buffer2;
    memcpy (buffer,defFileName.sData,defFileName.sLength+1);

    ofn.nMaxFileTitle   = 255;
    ofn.lpstrInitialDir = nil;
    ofn.lpstrDefExt     = nil;
    //ofn.lpfnHook      = sfwpHook;
    ofn.lpfnHook        = nil;

    if (lastWinPathUsed.sLength)
        ofn.lpstrInitialDir = lastWinPathUsed.getStr();
    else
        ofn.lpstrInitialDir = ((_String*)pathNames(0))->getStr();

    ofn.lpstrTitle      = prompt.getStr();
    ofn.Flags           = OFN_HIDEREADONLY|OFN_PATHMUSTEXIST|OFN_OVERWRITEPROMPT|OFN_EXPLORER;

    if (!GetSaveFileName (&ofn))
        return -1;

    fileName = buffer;
    return     ofn.nFilterIndex-1;*/

}

//________________________________________________________

int CALLBACK fontGeneratorProc (ENUMLOGFONT *lf, NEWTEXTMETRIC* , DWORD , LPARAM rec)
{
    _List * recList = (_List*)rec;
    _String fontName (lf->elfLogFont.lfFaceName);
    recList->BinaryInsert (&fontName);
    return 1;
}

//_________________________________________________________________________

void    GenerateFontList (_List& fonts)
{
    fonts.Clear();
    EnumFontFamilies (GetDC(nil), nil, (FONTENUMPROC)fontGeneratorProc, (LPARAM)&fonts);
}

//__________________________________________________________________________________
char    ScanDirectoryForFileNames (_String& source, _List& rec, bool recurse)
{
    _String searchSpec = source & "\\*";

    WIN32_FIND_DATA      findData;

    HANDLE  searchHandle = FindFirstFile (searchSpec.sData, &findData);

    if (searchHandle!=INVALID_HANDLE_VALUE) {
        while (1) {
            if (findData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
                if (recurse) {
                    _String subDirSpec =  findData.cFileName;
                    if (!subDirSpec.beginswith (".")) {
                        subDirSpec = source & '\\' & subDirSpec;
                        ScanDirectoryForFileNames (subDirSpec, rec, true);
                    }
                }
            } else {
                _String fullFileName = findData.cFileName;
                if (!fullFileName.beginswith (".")) {
                    fullFileName = source & '\\' & fullFileName;
                    rec && & fullFileName;
                }
            }
            if (!FindNextFile (searchHandle, &findData)) {
                break;
            }
        }
        FindClose (searchHandle);
    }
    return '\\';
}

//________________________________________________________

void    StartBarTimer(void)
{
    lastTimer   = time(&timerStart);
    updateTimer = true;
}

//________________________________________________________

void    StopBarTimer(void)
{
    updateTimer = false;
}

//________________________________________________________
void    ToggleAnalysisMenu (bool running)
{
    HMENU anMenu = GetSubMenu (GetMenu((HWND)hyphyConsoleWindow->GetOSWindowData()), 2);
    if (running) {
        EnableMenuItem (anMenu, 0, MF_ENABLED|MF_BYPOSITION);
        EnableMenuItem (anMenu, 1, MF_ENABLED|MF_BYPOSITION);
        EnableMenuItem (anMenu, 5, MF_GRAYED|MF_BYPOSITION);
        EnableMenuItem (anMenu, 6, MF_GRAYED|MF_BYPOSITION);
        EnableMenuItem (anMenu, 7, MF_GRAYED|MF_BYPOSITION);
    } else {
        EnableMenuItem (anMenu, 5, MF_ENABLED|MF_BYPOSITION);
        EnableMenuItem (anMenu, 6, MF_ENABLED|MF_BYPOSITION);
        EnableMenuItem (anMenu, 7, MF_ENABLED|MF_BYPOSITION);
        EnableMenuItem (anMenu, 0, MF_GRAYED|MF_BYPOSITION);
        EnableMenuItem (anMenu, 1, MF_GRAYED|MF_BYPOSITION);
        SetStatusBarValue (-1,1,0);
        SetStatusLine ("Idle");
    }
    DrawMenuBar ((HWND)hyphyConsoleWindow->GetOSWindowData());
}

//________________________________________________________
_HYColor        SelectAColor (_HYColor& currentColor, _String& )
{
    static              COLORREF     customColors[16];
    CHOOSECOLOR  cc;
    cc.lStructSize = sizeof (CHOOSECOLOR);
    cc.hwndOwner = GetFocus ();
    cc.hInstance = nil;
    cc.rgbResult = RGB (currentColor.R, currentColor.G, currentColor.B);
    cc.lpCustColors = customColors;
    cc.Flags = CC_ANYCOLOR|CC_RGBINIT|CC_FULLOPEN;
    cc.lCustData = nil;
    cc.lpfnHook  = nil;
    cc.lpTemplateName = nil;

    if (ChooseColor (&cc)) {
        _HYColor res;
        res.R = GetRValue (cc.rgbResult);
        res.G = GetGValue (cc.rgbResult);
        res.B = GetBValue (cc.rgbResult);
        return res;
    }
    return currentColor;
}

//________________________________________________________
long    GetMaxCharWidth (_HYFont& f)
{
    _String testString ('W');
    return  GetVisibleStringWidth (testString, f)+2;
}

//________________________________________________________

char        YesNoCancelPrompt (_String& prompt)
{
    int res = MessageBox (GetForegroundWindow(), prompt.sData, "HyPhy Prompt", MB_YESNOCANCEL|MB_ICONQUESTION|
                          MB_DEFBUTTON1|MB_APPLMODAL);
    if (res == IDYES) {
        return 1;
    }
    if (res == IDCANCEL) {
        return 2;
    }

    return 3;
}

//________________________________________________________

char    * ReturnFileDialogSelectionWin (bool write, _String* initDir)
{
    static char buffer[2049], buffer2[256];
    OPENFILENAME ofn;
    ofn.lStructSize = sizeof (OPENFILENAME);
    ofn.hwndOwner = GetForegroundWindow();

    char  fileFilters[] = "HYPHY Batch Files\0*.bf\0Data File\0*.dat;*.nuc;*.seq;*.phy;*.nex;\0All Files\0*.*\0\0";
    ofn.lpstrFilter = fileFilters;
    ofn.lpstrCustomFilter = nil;

    ofn.nFilterIndex = 3;

    buffer2[0] = 0;
    if (defFileNameValue.sLength) {
        if (defFileNameValue.sLength > 2048) {
            _String truncString (defFileNameValue,0,2047);
            strcpy (buffer, truncString.sData);
        } else {
            strcpy (buffer, defFileNameValue.sData);
        }
    } else {
        buffer[0] = 0;
    }

    ofn.lpstrFile = buffer;
    ofn.nMaxFile = 2048;
    ofn.lpstrFileTitle = buffer2;
    ofn.nMaxFileTitle = 255;

    if (initDir) {
        ofn.lpstrInitialDir = initDir->sData;
    } else {
        if (lastWinPathUsed.sLength) {
            ofn.lpstrInitialDir = lastWinPathUsed.getStr();
        } else if (pathNames.lLength) {
            ofn.lpstrInitialDir = ((_String*)pathNames(0))->getStr();
        } else {
            ofn.lpstrInitialDir = baseDirectory.getStr();
        }
    }

    ofn.lpstrDefExt = nil;
    ofn.lpstrTitle = dialogPrompt.getStr();
    ofn.Flags = 0;
    ofn.Flags|=OFN_FILEMUSTEXIST;
    ofn.Flags|=OFN_HIDEREADONLY;
    ofn.Flags|=OFN_PATHMUSTEXIST;
    /*char defExt[] = ".bf";*/


    bool         res;

    if (write) {
        res = GetSaveFileName (&ofn);
    } else {
        res = GetOpenFileName (&ofn);
    }

    if (res) {
        lastFileTypeSelection = ofn.nFilterIndex;
        lastWinPathUsed = buffer;
        lastWinPathUsed = lastWinPathUsed.Cut (0,ofn.nFileOffset-1);
        return            buffer;
    } else {
        return            empty.sData;
    }

}

//________________________________________________________

bool    PopUpFileDialog(_String ps, _String* defaultLocation)
{
    _String saveDP = dialogPrompt,
            res;

    dialogPrompt = ps;
    res = ReturnFileDialogSelectionWin (false, defaultLocation);
    if (argFileName) {
        *argFileName = res;
    } else {
        argFileName = new _String (res);
        checkPointer (argFileName);
    }

    dialogPrompt = saveDP;
    return argFileName->sLength;
}

//_________________________________________________________________________

void    MoveConsoleWindow (_HYRect& newLoc)
{
    //SetWindowPos(WindowHandle, nil, newLoc.left, newLoc.top,newLoc.right-newLoc.left,newLoc.bottom-newLoc.top,
    //SWP_NOZORDER);

    hyphyConsoleWindow->SetPosition         (newLoc.left, newLoc.top);
    hyphyConsoleWindow->SetWindowRectangle  (newLoc.top,newLoc.left,newLoc.bottom,newLoc.right,true);
}

//_________________________________________________________________________

_String ChooseAFolder       (_String& prompt)
{
    BROWSEINFO BI;
    ITEMIDLIST *IDL;

    char       cDirName [MAX_PATH+1];

    memset(&BI, 0, sizeof(BI));

    BI.hwndOwner = GetFocus();
    BI.pszDisplayName = cDirName;
    BI.lpszTitle = prompt.sData;
    BI.ulFlags = BIF_RETURNONLYFSDIRS;

    if( IDL = SHBrowseForFolder(&BI) ) {
        LPMALLOC pMalloc;
        BOOL bOK;
        SHGetMalloc(&pMalloc);
        bOK = SHGetPathFromIDList(IDL, cDirName);
        pMalloc->Free(pMalloc);
        pMalloc->Release();
        return cDirName;
    }
    return empty;
}

//_________________________________________________________________________

void    ShowObjectInspector (void)
{
    long f = FindWindowByName (objectInspectorTitle);
    if (f>=0) {
        SetForegroundWindow ((HWND)windowPtrs (f));
    } else {
        _HYObjectInspector* newOI = new _HYObjectInspector ();
        newOI->Activate       ( );
    }
}

//_________________________________________________________________________
void    PositionWindow          (_HYGuiObject* twp, _String* args)
{
    _List * argL = args->Tokenize (",");
    _HYWindow*   tw = (_HYWindow*)twp;
    if (argL->lLength>=4) {
        long R[5],
             k;

        for (k=0; k<4; k++) {
            R[k] = ((_String*)(*argL)(k))->toNum();
        }
        if (argL->lLength>4) {
            R[4] = ((_String*)(*argL)(4))->toNum();
        } else {
            R[4] = 0;
        }

        _HYRect   wR = GetScreenDimensions  (),
                  wiR;
        long      W[4] = {wR.left,wR.top, wR.right, wR.bottom};
        for (k=0; k<4; k++)
            if (R[k]<0) {
                R[k] = W[k] + ((k<2)?-1:1)*R[k];
            }

        wiR.left    = R[0];
        wiR.right   = R[2];
        wiR.top     = R[1];
        wiR.bottom  = R[3];

        if (wiR.left>=wiR.right) {
            wiR.right = 1+wiR.left;
        }
        if (wiR.top>=wiR.bottom) {
            wiR.bottom = 1+wiR.top;
        }

        tw->SetPosition        (wiR.left,wiR.top);
        tw->SetWindowRectangle (0,0,wiR.bottom-wiR.top,wiR.right-wiR.left);

        if (R[4]>0) {
            wiR.top         = wiR.bottom+2;
            wiR.bottom      = wR.bottom - 2;
            wiR.left        = 5;
            wiR.right       = wR.right - 2;
            MoveConsoleWindow (wiR);
        }

    }

    DeleteObject (argL);
}

//_________________________________________________________________________

bool    Get_a_URL (_String& urls, _String* fileName)
{

    DWORD dwSize     = 0,
          dataChunk = 65536,
          totalRead = 0;

    HINTERNET   hRootHandle = InternetOpen   (GetVersionString().sData,INTERNET_OPEN_TYPE_DIRECT ,nil,nil,0),
                hUrlDump;


    FILE*       f = nil;
    _String*     storage = nil;

    if (hRootHandle) {
        if (fileName) {
            f = fopen (fileName->sData,"wb");
            if (!f) {
                urls = "Failed to create/open target file";
                return false;
            }
        } else {
            checkPointer (storage = new _String (dataChunk,true));
        }

        hUrlDump = InternetOpenUrl(hRootHandle, urls.sData, NULL, NULL, INTERNET_FLAG_RAW_DATA|INTERNET_FLAG_NO_CACHE_WRITE, 0);

        SetStatusLine (_String("Retrieving ")& urls);
        do {
            Ptr buffer = MemAllocate (dataChunk+1);
            if(!InternetReadFile(hUrlDump,(LPVOID)buffer,dataChunk,&dwSize)) {
                free (buffer);
                break;
            } else {
                buffer[dwSize]='\0';

                if (dwSize == 0) {
                    if (f) {
                        fclose (f);
                    }

                    if (storage) {
                        storage->Finalize();
                        urls = *storage;
                        DeleteObject (storage);
                    }

                    free (buffer);
                    InternetCloseHandle (hRootHandle);
                    SetStatusLine     ("Idle");
                    return true;
                } else {
                    if (storage) {
                        (*storage) << buffer;
                    } else {
                        fwrite (buffer,1,dwSize,f);
                    }

                    totalRead+=dwSize;

                    SetStatusLine (_String("Read ") & (long)(totalRead/1024) & "KB from " & urls);
                }
            }

        } while (true);

        InternetCloseHandle (hRootHandle);
    }


    SetStatusLine     ("Idle");

    LPVOID lpMsgBuf;

    if (!FormatMessage( FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
                        NULL,GetLastError(),MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPTSTR) &lpMsgBuf,0,NULL )) {
        urls = "Unknown system error";
    }

    if (f) {
        fclose (f);
    }

    if (storage) {
        DeleteObject (storage);
    }

    urls = (char*)lpMsgBuf;
    LocalFree( lpMsgBuf );

    return false;
}