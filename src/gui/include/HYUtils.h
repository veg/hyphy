/*

    Handy OS utils



    Sergei L. Kosakovsky Pond, June 2000-November 2004.

*/



#ifndef _HYOSUTILS_

#define _HYOSUTILS_

//#pragma once

#include "hy_strings.h"

#include "HYBaseGUI.h"

#include "calcnode.h"


Ptr     ProcureIconResource     (long);

long    GetVisibleStringWidth   (_String&, _HYFont&);

long    GetVisibleStringWidth   (_String&);

long    GetMaxCharWidth         (_HYFont&);

_HYRect GetScreenDimensions     (void);

void    CenterWindow            (_HYGuiObject*);

void    DelayNMs                (long);

void    PositionWindow          (_HYGuiObject*,_String*);

void    ToggleAnalysisMenu      (bool);

long    SaveFileWithPopUp       (_String&, _String&, _String&, _String&, _List&);

bool    OpenTreeFile            (void);

bool    OpenTextFile            (void);

bool    OpenTable               (void);

bool    OpenDataFile            (_String* defLoc = nil);

bool    OpenModelFile           (_String* defLoc = nil);

bool    ExecuteBatchFile        (void);

bool    PopUpFileDialog         (_String, _String *defLoc = nil);

void    NewChartWindow          (void);

void    NewModel                (_String*);

bool    OpenTable               (void);

bool    OpenModelFile           (_String*);

void    ShowObjectInspector     (void);

void    PlaceStringInClipboard  (_String&,Ptr);





#ifdef __MAC__

#include <Menus.h>

_String MacSimpleFileSave       (void);
_String MacSimpleFileOpen       (void);
_String ChooseAFolder (_String& promptString);


void    StringToStr255      (_String&, Str255&);

void    Str255ToStr         (_String&, Str255&);

void    DrawMenuPlaceHolder (Rect&, _String&, bool enabled = true);

void    DrawEmbossedBox     (Rect&);

void    DrawInfoBox         (Rect&, _String&, _String&);

void    SetWindowFont       (short,short,Style, bool);

void    TreeDependencies    (_SimpleList&, long);

void    DSDependencies      (_SimpleList&, long);

void    ModelDependencies   (_SimpleList&, long);

void    ListToPopUpMenu     (_List&, MenuHandle);

void    UpdateStatusLine    (Ptr theWindow);

long    HandleListSelection (_List&, _SimpleList&, _SimpleList&, Str63, _SimpleList&, long);

_String NewTreeWindow       (long sourceDF = -1);

void    ConvertMovieFile    (bool);



#ifdef  __HYPHYXCODE__

_String DoMacToPOSIX        (const _String&);

#endif

#endif



#ifdef __WINDOZE__



#include    <Windows.h>

long        SaveFileFunction                (_String&, _String&, _String&, char*, Ptr);

HMENU       ListToPopUpMenu                 (_List&,long base = 0);

char        *ReturnFileDialogSelectionWin   (bool write, _String* initDir = nil);

void        PlaceBitmapInClipboard          (HBITMAP,HWND);



#endif



_String     ChooseAFolder                   (_String&);

char        YesNoCancelPrompt               (_String&);



long        FindWindowByName                (_String&);

_HYGuiObject*

FindWindowByNameAndOpen         (_String&);

_HYGuiObject*

FindWindowByID                  (long);



_HYColor    SelectAColor                    (_HYColor&,_String&);



#endif



void        MoveConsoleWindow               (_HYRect&);

_String     HandlePullDown                  (_List&,long,long,long);

long        HandlePullDownWithFont          (_List&,long,long,long,_String,long);



char        ScanDirectoryForFileNames       (_String&, _List&, bool);

// directory to scan, list to receive full path names for each file,

// scan recursively or not.

_HYColor    GetDialogBackgroundColor        (void);



void        StartBarTimer                   (void);

void        StopBarTimer                    (void);



void        GenerateFontList                (_List&);



extern      _String                         menuSeparator,

            *argFileName;



//EOF
