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

#include "likefunc.h"
#include "preferences.h"
#include "list.h"
#include "avllistx.h"
#include "simplelist.h"
#include "preferences.h"

#ifndef  __HEADLESS__
#include "HYUtils.h"
#include "HYConsoleWindow.h"
#include "HYTextBox.h"
#include "HYCanvas.h"
#endif

// preferences -> index mapping


#ifdef __MAC__
short           prefFileRefNumber   = 0;
#include        <Folders.h>
#include        <Script.h>

#ifndef __HYPHYXCODE__
extern  long    _fcreator,
        _ftype;
#endif

extern      MenuHandle recentFilesMenu;
#endif

#ifdef __WINDOZE__
FILE*           prefFileHandle      = nil;
_String         hyWindowsConsoleFont ("CONSOLE_FONT_INFO"),
                prefFileTitle        ("HYPHYSettings.ini");
#ifndef __HEADLESS__
extern          HMENU                 recentFilesMenu;
#endif
#endif

#ifdef __HYPHY_GTK__
#endif

_List               _hyPreferencesKeysAux;
_AVLListX           _hyPreferencesKeys (&_hyPreferencesKeysAux);

_String             _hyPreferencesFontFace      ("Font Face"),
                    _hyPreferencesFontSize       ("Font Size"),
                    _hyPreferencesFontStyle      ("Font Style"),
                    _hyPreferencesPrecision      ("Precision"),
                    _hyPreferencesPersistence    ("Persistence"),
                    _hyPreferencesInitGuess      ("Initial Guess"),
                    _hyPreferencesStartValue ("Starting Value"),
                    _hyPreferencesDeletions      ("Deletions"),
                    _hyPreferencesFormat     ("Output format"),
                    _hyPreferencesLineWidth      ("Line width"),
                    _hyPreferencesGapWidth       ("Gap width"),
                    _hyPreferencesHetSimulation ("Heterogeneity Simulation"),
                    _hyPreferencesRandomSeed ("Random seed"),
                    _hyPreferencesLikelihoodDisp("Likelihood Display"),
                    _hyPreferencesNumberFormat  ("Number Format"),
                    _hyPreferencesTreeDisplay    ("Tree Display"),
                    _hyPreferencesOptProgress    ("Optimization Progress"),
                    _hyPreferencesStartupDialog  ("Startup Dialog"),
                    _hyPreferencesGapFreqs       ("Include Gaps in Frequency Counts"),
                    _hyPreferencesAutomoveC      ("Automove console"),
                    _hyPreferencesMP         ("CPUs to Load");
//____________________________________________________________________________________________

bool                showDialogAtStartup         = true,
                    doAutoConsoleMove             = true;


#ifndef             __HEADLESS__
_HYRect             consolePositionRectangle = {50,50,700,500,0};
#endif

_String             initialDialogPop            ("SHOW_DIALOG_AT_STARTUP"),
                    windowPositionStr           ("CONSOLE_WINDOW_COORDS"),
                    autoMoveConsole             ("AUTOMATICALLY_RESIZE_CONSOLE"),
                    recentFilesList             ("RECENT_FILES_LIST"),
                    treeDisplayOptions          ("TREE_DISPLAY_OPTIONS");

_List               recentFiles,
                    recentPaths,
                    globalPreferencesList;

//____________________________________________________________________________________________

void  ReadPreferences (void)
{
    _String optionList,
            menuKey,
            comma(",");

    // populate the list of preferences

#if !defined __WINDOZE__ && !defined __HEADLESS__
    _List                fonts;
    AddItemToPreferences (1|8,-1,"Console Font Settings","Choose font face, size and style used for displaying text in the console.","",nil,globalPreferencesList,false);
    GenerateFontList     (fonts);

    _hyPreferencesKeys.Insert(_hyPreferencesFontFace.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_POPUP,_hyPreferencesFontFace,"Select the font used to display text in the console window.","Monaco",&fonts,globalPreferencesList,false);
    optionList = "8,9,10,12,14,18";
    _hyPreferencesKeys.Insert(_hyPreferencesFontSize.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_POPUP,_hyPreferencesFontSize,"Select font size used to display text in the console window.","9",optionList.Tokenize (comma), globalPreferencesList,true);
    optionList = "Plain,Bold,Italic";
    _hyPreferencesKeys.Insert(_hyPreferencesFontStyle.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_POPUP,_hyPreferencesFontStyle,"Select font style used to display text in the console window.","Plain",optionList.Tokenize (comma), globalPreferencesList,true);
#endif

    AddItemToPreferences (1|8,-1,"Optimization Settings","Options affecting the optimization algorithm.","",nil,globalPreferencesList,false);
    _hyPreferencesKeys.Insert(_hyPreferencesPrecision.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_TEXTBOX,_hyPreferencesPrecision,"Desired precision(absolute error) in ln-likelihood value. Settings between 0.1 and 0.000000001 are recommended.",
                          "0.001",nil,globalPreferencesList,false);
    optionList = "Low,Normal,High,Very High";
    _hyPreferencesKeys.Insert(_hyPreferencesPersistence.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_POPUP,_hyPreferencesPersistence,"Controls the number iterations the optimization algorithm will perform before it terminates if the desired precision is not met.",
                          "Normal",optionList.Tokenize (comma), globalPreferencesList,true);
    optionList = "Do not use distances,Use distances";
    _hyPreferencesKeys.Insert(_hyPreferencesInitGuess.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_POPUP,_hyPreferencesInitGuess,"Determines whether distance methods are to be used to obtain intial parameter value guesses. Applies only to nuceleotide models.",
                          "Use distances",optionList.Tokenize (comma), globalPreferencesList,true);
    _hyPreferencesKeys.Insert(_hyPreferencesStartValue.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_TEXTBOX,_hyPreferencesStartValue,"Sets starting values for parameters for optimization routines. If starting values are obtained by distance methods, this option is ignored.",
                          "0.1",nil,globalPreferencesList,false);
    AddItemToPreferences (1|8,-1,"Data Read/Write Settings","Options affecting sequence data files reading and writing.","",nil,globalPreferencesList,false);
    optionList = "Skip Deletions,Keep Deletions";
    _hyPreferencesKeys.Insert(_hyPreferencesDeletions.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_POPUP,_hyPreferencesDeletions,"Choose \"Keep Deletions\" to retain deletions (as ambiguities) for analyses. \"Skip Deletions\" filters deletions out as the data is read.","Keep Deletions",
                          optionList.Tokenize (comma), globalPreferencesList,true);
    optionList = "# sequential,# interleaved,PHYLIP Sequential,PHYLIP Interleaved,NEXUS sequential with labels,NEXUS interleaved with labels,NEXUS sequential without labels,NEXUS interleaved without labels,Comma separated characters,FASTA sequential,FASTA interleaved";
    _hyPreferencesKeys.Insert(_hyPreferencesFormat.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_POPUP,_hyPreferencesFormat,"Choose the default file format for data filters output to files via fprintf.","NEXUS sequential without labels",
                          optionList.Tokenize (comma), globalPreferencesList,true);

    _hyPreferencesKeys.Insert(_hyPreferencesLineWidth.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_TEXTBOX,_hyPreferencesLineWidth,"This options sets how many characters will be printed per line for data filters output to files via fprintf. Only affects interleaved formats.",
                          "50",nil,globalPreferencesList,false);
    _hyPreferencesKeys.Insert(_hyPreferencesGapWidth.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_TEXTBOX,_hyPreferencesGapWidth,"This options sets how many characters will be printed per cluster (clusters are separated by spaces) for data filters output to files via fprintf. Only affects interleaved non-NEXUS formats.",
                          "10",nil,globalPreferencesList,false);

    optionList = "Yes,No";
    _hyPreferencesKeys.Insert(_hyPreferencesGapFreqs.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_POPUP,_hyPreferencesGapFreqs,"Include gaps (as fully unresolved characters) in frequency counts","Yes",
                          optionList.Tokenize (comma), globalPreferencesList,true);

    AddItemToPreferences (1|8,-1,"Simulation Options","Options affecting bootstrapping algorithms.","",nil,globalPreferencesList,false);
    optionList = "Discrete Distribution,Continuous Distribution";
    _hyPreferencesKeys.Insert(_hyPreferencesHetSimulation.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_POPUP,_hyPreferencesHetSimulation,"When bootstrapping models with heterogeneous rates, determines whether rate classes are drawn from the continuous (e.g. gamma) distribution or it's discrete approximation.",
                          "Continuous Distribution", optionList.Tokenize (comma), globalPreferencesList,true);

    _hyPreferencesKeys.Insert(_hyPreferencesRandomSeed.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_TEXTBOX,_hyPreferencesRandomSeed,"Set this parameter to -1 to have HYPHY seed random generator anew every time the program is run. A positive value defines the seed to be used instead. Changes will take effect when HYPHY is restarted.",
                          "-1",nil,globalPreferencesList,false);

    AddItemToPreferences (1|8,-1,"Miscellaneous Options","Variuos, primarily formatting, options.","",nil,globalPreferencesList,false);
    optionList = "Function value only,Complete report as list,Tree with branch lengths,Parameters and Constraints,Batch Language Statement,Batch Language Statement with Trees";
    _hyPreferencesKeys.Insert(_hyPreferencesLikelihoodDisp.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_POPUP,_hyPreferencesLikelihoodDisp,"Various ways to display likelihood function and parameters","Tree with branch lengths",
                          optionList.Tokenize (comma), globalPreferencesList,true);

    optionList = "Short,Normal,Long,Maximally Long";
    _hyPreferencesKeys.Insert(_hyPreferencesNumberFormat.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_POPUP,_hyPreferencesNumberFormat,"Determines how many significant digits are displayed when printing numbers via fprintf.","Normal",
                          optionList.Tokenize (comma), globalPreferencesList,true);

    optionList = "No auto display,Auto display single tree,Auto display all trees";

    _hyPreferencesKeys.Insert(_hyPreferencesTreeDisplay.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_POPUP,_hyPreferencesTreeDisplay,"Should HY-PHY automatically open graphical tree windows upon completion of an analysis.","Auto display single tree",
                          optionList.Tokenize (comma), globalPreferencesList,true);

    optionList = "Silent,Verbose";
    _hyPreferencesKeys.Insert(_hyPreferencesOptProgress.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_POPUP,_hyPreferencesOptProgress,"Triggers the optimization functions to print out progress lines while obtaining MLEs","Silent",
                          optionList.Tokenize (comma), globalPreferencesList,true);
    optionList = "Yes,No";
    _hyPreferencesKeys.Insert(_hyPreferencesStartupDialog.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_POPUP,_hyPreferencesStartupDialog,"Display an action dialog when HyPhy starts up","Yes",
                          optionList.Tokenize (comma), globalPreferencesList,true);
    optionList = "Yes,No";
    _hyPreferencesKeys.Insert(_hyPreferencesAutomoveC.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_POPUP,_hyPreferencesAutomoveC,"Automatically move and resize console window when a data panel is opened","Yes",
                          optionList.Tokenize (comma), globalPreferencesList,true);

#ifdef       __MP__
    AddItemToPreferences (1|8,-1,"Multiple Processors","Multiple Processor Settings","",nil,globalPreferencesList,false);
    _List options;
    for (long cpCount = 1; cpCount <= systemCPUCount; cpCount++) {
        _String cc (cpCount);
        options && & cc;
    }
    _hyPreferencesKeys.Insert(_hyPreferencesMP.makeDynamic(),((_List*)globalPreferencesList.lData[4])->lLength);
    AddItemToPreferences (0,PREFITEM_POPUP,_hyPreferencesMP,"How many CPUs should HYPHY try to load.",_String(systemCPUCount),&options,globalPreferencesList,false);
#endif

    _String * fileContents = nil;

#ifdef __HEADLESS__
    optionList       = baseDirectory & ".init";
    FILE * hpp       = doFileOpen (optionList.getStr(),"rb");
    if (hpp) {
        fileContents = new _String (hpp);
        fclose (hpp);
    }
#else
#ifdef __MAC__
    Str255      buffer="\pHYPHYKernel Preferences";
    short       prefVolID,
                osErr;

    long        prefFolderID = 0;

    FindFolder (kOnSystemDisk,kPreferencesFolderType,kCreateFolder,&prefVolID,&prefFolderID);
    FSSpec fileSpec;
    FSMakeFSSpec(prefVolID,prefFolderID,buffer,&fileSpec);
    if ((osErr=FSpOpenDF (&fileSpec,fsRdWrPerm,&prefFileRefNumber))!=noErr) {
#ifndef __HYPHYXCODE__
        if ((osErr=FSpCreate(&fileSpec,_fcreator,'pref',GetScriptManagerVariable(smSysScript)))!=noErr)
#else
        if ((osErr=FSpCreate(&fileSpec,'MuSe','pref',GetScriptManagerVariable(smSysScript)))!=noErr)
#endif
            prefFileRefNumber = -1;
        else if ((osErr=FSpOpenDF (&fileSpec,fsRdWrPerm,&prefFileRefNumber))!=noErr) {
            prefFileRefNumber = -1;
        }
    }
    if (prefFileRefNumber!=-1) {
        GetEOF(prefFileRefNumber,&prefFolderID);
        if (prefFolderID>0) {
            fileContents = new _String ((unsigned long)prefFolderID);
            FSRead (prefFileRefNumber,&prefFolderID,fileContents->sData);
            fileContents->sData[prefFolderID] = 0;
        }
    }
#endif

#ifdef __WINDOZE__
    _String  prefFileName = baseDirectory&prefFileTitle;

    if ((prefFileHandle = fopen (prefFileName.getStr(),"r+"))) {
        fileContents = new _String (prefFileHandle);
        fclose (prefFileHandle);
    }

#endif

#ifdef __HYPHY_GTK__
    optionList       = baseDirectory & ".hyphyprefs";
    FILE * prefsFile = fopen (optionList.sData,"rb");

    if (prefsFile) {
        fileContents = new _String (prefsFile);
        fclose (prefsFile);
    }
#endif
#endif


    _List terms,
          * availNames = (_List*)globalPreferencesList(1),
            * availValues = (_List*)globalPreferencesList(4);

    if (fileContents) {
        _ElementaryCommand::ExtractConditions(*fileContents,0,terms);
        for (long k=0; k<terms.lLength; k++) {
            _List    theTerms;
            _ElementaryCommand::ExtractConditions(*((_String*)terms(k)),0,theTerms,'=');
            if (theTerms.lLength == 2) {
                _String * prefID = (_String*)theTerms.lData[0];
                long    j = availNames->Find (prefID);

                if (j>=0) {
                    *((_String*)availValues->lData[j]) = *(_String*)theTerms.lData[1];
                } else {
                    if (prefID->Equal (&windowPositionStr)) {
#ifndef __HEADLESS__
                        _List*  rectSizes = ((_String*)theTerms.lData[1])->Tokenize (",");
                        if (rectSizes->lLength==4) {
                            consolePositionRectangle.left   = MAX(((_String*)(*rectSizes)(0))->toNum(),0);
                            consolePositionRectangle.right  = MAX(((_String*)(*rectSizes)(2))->toNum(),consolePositionRectangle.left+50);

                            consolePositionRectangle.top    = MAX(((_String*)(*rectSizes)(1))->toNum(),0);
                            consolePositionRectangle.bottom = MAX(((_String*)(*rectSizes)(3))->toNum(),consolePositionRectangle.top+50);
                        }
                        DeleteObject (rectSizes);
#endif
                    } else if (prefID->Equal (&recentFilesList)) {
                        _List* filePaths = ((_String*)theTerms.lData[1])->Tokenize (",");
                        for (long idx = 0; idx<filePaths->lLength; idx+=2) {
                            _String * sn = (_String*)(*filePaths)(idx),
                                      * sp = (_String*)(*filePaths)(idx+1);

                            sn->StripQuotes();
                            sp->StripQuotes();
                            AddStringToRecentMenu(*sn,*sp);
                        }
                        DeleteObject (filePaths);
                    }
#if defined __WINDOZE__ && !defined __HEADLESS__
                    else if (prefID->Equal (&hyWindowsConsoleFont)) {
                        _List fontList;
                        _ElementaryCommand::ExtractConditions(*(_String*)theTerms.lData[1],0,fontList,',');
                        if (fontList.lLength>=2) {
                            _HYFont cF;
                            cF.face     = *(_String*)fontList.lData[0];
                            cF.size     = ((_String*)fontList.lData[1])->toNum();
                            if (fontList.lLength > 2) {
                                cF.style    = ((_String*)fontList.lData[2])->toNum();
                            } else {
                                cF.style    = HY_FONT_PLAIN;
                            }

                            ((_HYTextBox*)hyphyConsoleWindow->GetObject(0))->SetFont(cF);
                            ((_HYTextBox*)hyphyConsoleWindow->GetObject(1))->SetFont(cF);
                        }
                    }
#endif
                }
            }
        }
    }

    long    aS = ((_String*)(availValues->lData[_hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesRandomSeed))]))->toNum();
    if (aS>=0) {
        init_genrand (aS);
    }

    showDialogAtStartup = *(((_String*)(availValues->lData[_hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesStartupDialog))])))==_String("Yes");
    doAutoConsoleMove   = *(((_String*)(availValues->lData[_hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesAutomoveC))])))==_String("Yes");

    DeleteObject        (fileContents);
}

//____________________________________________________________________________________________

void  SetPreferences (void)
{
#ifndef __HEADLESS__
#if defined __MAC__ or defined __HYPHY_GTK__
    if (!hyphyConsoleWindow) {
        return;
    }

    _String*        prefValue,
                    prefKey;

    _HYFont         newFont;
    _List*          prefValues = (_List*)globalPreferencesList(4);

    newFont.face = *(_String*)prefValues->lData[_hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesFontFace))];
    // font face

    newFont.size = ((_String*)prefValues->lData[_hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesFontSize))])->toNum();
    // font size

    if (newFont.size<4 || newFont.size>100) {
        newFont.size = ((_HYTextBox*)hyphyConsoleWindow->GetObject (0))->GetFont().size;    // use existing font size
    }

    prefValue   = (_String*)prefValues->lData[_hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesFontStyle))];
    // font style

    if (*prefValue == _String ("Bold")) {
        newFont.style = HY_FONT_BOLD;
    } else if (*prefValue == _String ("Italic")) {
        newFont.style = HY_FONT_ITALIC;
    } else {
        newFont.style = HY_FONT_PLAIN;
    }

    ((_HYTextBox*)hyphyConsoleWindow->GetObject (0))->SetFont (newFont); // main panel
    ((_HYTextBox*)hyphyConsoleWindow->GetObject (1))->SetFont (newFont); // user input panel
#endif
#endif
}


//____________________________________________________________________________________________

void  WritePreferences (void)
{
#ifdef __MAC__
    if (prefFileRefNumber!=-1)
#endif

#ifdef __WINDOZE__
        _String  prefFileName = baseDirectory&prefFileTitle;
        prefFileHandle = fopen (prefFileName.sData,"w+");
    if (prefFileHandle)
#endif

#ifdef __HYPHY_GTK__
        _String  prefFileName = baseDirectory & ".hyphyprefs";
    FILE * prefFile = fopen (prefFileName.sData, "wb");
    if (prefFile)
#endif

    {
        _String      spoolResult((unsigned long)256, true);
        _SimpleList* pCodes = (_SimpleList*)globalPreferencesList(0);

        _List *pNames = (_List*)globalPreferencesList(1),
               *pValues = (_List*)globalPreferencesList(4);
        for (long i = 0; i<pCodes->lLength; i++)
            if (pCodes->lData[i]<8) { // not a rubrik header
                spoolResult<<(_String*)pNames->lData[i];
                spoolResult<<'=';
                spoolResult<<(_String*)pValues->lData[i];
                spoolResult<<';';
            }

#ifndef __HEADLESS__
        // convert the coordinate info to a comma separated string
        _HYRect    wr = hyphyConsoleWindow->GetWindowRect();

        spoolResult<<windowPositionStr;
        spoolResult<<'=';
        spoolResult<<_String ((long)wr.left);
        spoolResult<<',';
        spoolResult<<_String ((long)wr.top);
        spoolResult<<',';
        spoolResult<<_String ((long)wr.right);
        spoolResult<<',';
        spoolResult<<_String ((long)wr.bottom);
        spoolResult<<';';

        if (recentFiles.lLength) {
            spoolResult << recentFilesList;
            spoolResult << '=';
            for (long k=0; k<recentFiles.lLength; k++) {
                if (k) {
                    spoolResult << ',';
                }
                spoolResult << '"';
                spoolResult << (_String*)recentFiles(k);
                spoolResult << "\",\"";
                spoolResult << (_String*)recentPaths(k);
                spoolResult << '"';
            }
            spoolResult << ';';
        }
#endif

#if defined __WINDOZE__ && ! defined __HEADLESS__
        _HYFont cF  = ((_HYTextBox*)hyphyConsoleWindow->GetObject(0))->GetFont();
        spoolResult<< hyWindowsConsoleFont;
        spoolResult<< '=';
        spoolResult<< &cF.face;
        spoolResult<< ',';
        spoolResult<< _String((long)cF.size);
        spoolResult<< ',';
        spoolResult<< _String((long)cF.style);
        spoolResult<< ';';
#endif

        spoolResult.Finalize();
        StringToConsole (spoolResult);


#ifdef __MAC__
        SetEOF (prefFileRefNumber,0);
        FSWrite(prefFileRefNumber,(long*)&spoolResult.sLength,spoolResult.sData);
        FSClose(prefFileRefNumber);
#endif

#ifdef __WINDOZE__
        fwrite (spoolResult.sData,1,spoolResult.sLength,prefFileHandle);
        fflush (prefFileHandle);
        fclose (prefFileHandle);
#endif

#ifdef __HYPHY_GTK__
        fwrite (spoolResult.sData, spoolResult.sLength, 1, prefFile);
        fclose (prefFile);
#endif

    }
}

//____________________________________________________________________________________________

void  ApplyPreferences (void)
{
    _List           *pfValues = (_List*)globalPreferencesList.lData[4],
                     *t;

    long            keyIndex = _hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesPersistence));

    setParameter    (optimizationPrecision, ((_String*)pfValues->lData[_hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesPrecision))])->toNum());

    t               = (_List*)((_List*)globalPreferencesList.lData[5])->lData[keyIndex];
    long    f       = t->Find (((_String*)pfValues->lData[keyIndex]));

    switch (f) {
    case 0:
        f = 200;
        break;
    case 1:
        f = -1;
        break;
    case 2:
        f = 2000;
        break;
    case 3:
        f = 50000;
    }

    if (f>0) {
        setParameter (maximumIterationsPerVariable, (_Parameter)f);
    }

    setParameter (globalStartingPoint, ((_String*)pfValues->lData[_hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesStartValue))])->toNum());
    setParameter (useInitialDistanceGuess, (_Parameter)(*((_String*)pfValues->lData[_hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesInitGuess))])==_String("Use distances")));

    setParameter (skipOmissions, (_Parameter)(*((_String*)pfValues->lData[_hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesDeletions))])==_String("Skip Deletions")));

    keyIndex = _hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesFormat));
    t = (_List*)((_List*)globalPreferencesList.lData[5])->lData[keyIndex];
    f = t->Find (((_String*)pfValues->lData[keyIndex]));
    if (f<0) {
        f = 6;
    }
    setParameter (dataFilePrintFormat, (_Parameter)f);
    setParameter (dataFileDefaultWidth , (long)(((_String*)pfValues->lData[_hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesLineWidth))])->toNum()));
    setParameter (dataFileGapWidth ,     (long)(((_String*)pfValues->lData[_hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesGapWidth))])->toNum()));

    setParameter (categorySimulationMethod, (_Parameter)((*((_String*)pfValues->lData[_hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesHetSimulation))])
                  !=_String("Discrete Distribution")))+1.0);

    setParameter (hfCountGap, (_Parameter)((*((_String*)pfValues->lData[_hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesGapFreqs))])
                                            !=_String("No"))));

    f = ((_String*)pfValues->lData[_hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesRandomSeed))])->toNum();
    if (f>=0) {
        setParameter (randomSeed, f);
    }

    keyIndex = _hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesLikelihoodDisp));
    t = (_List*)((_List*)globalPreferencesList.lData[5])->lData[keyIndex];
    f = t->Find (((_String*)pfValues->lData[keyIndex]));
    if (f<0) {
        f = 0;
    }
    setParameter (likefuncOutput, (_Parameter)f);

    keyIndex = _hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesNumberFormat));
    t = (_List*)((_List*)globalPreferencesList.lData[5])->lData[keyIndex];
    f = t->Find (((_String*)pfValues->lData[keyIndex]));
    if (f<0) {
        f = 1;
    }
    switch (f) {
    case 0:
        printDigits = 5;
        break;
    case 2:
        printDigits = 12;
        break;
    case 3:
        printDigits = 15;
    }
    setParameter (printDigitsSpec,printDigits);

    keyIndex = _hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesTreeDisplay));
    t = (_List*)((_List*)globalPreferencesList.lData[5])->lData[keyIndex];
    f = t->Find (((_String*)pfValues->lData[keyIndex]));
    setParameter (treeDisplayOptions, f);

    keyIndex = _hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesOptProgress));
    t = (_List*)((_List*)globalPreferencesList.lData[5])->lData[keyIndex];
    verbosityLevel      = t->Find (((_String*)pfValues->lData[keyIndex]))?5:(-1);
    setParameter        (VerbosityLevelString,verbosityLevel);

    doAutoConsoleMove   = *(((_String*)((_List*)globalPreferencesList.lData[4])->lData[_hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesAutomoveC))]))==_String("Yes");
#ifdef       __MP__
    systemCPUCount= ((_String*)pfValues->lData[_hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesMP))])->toNum();
    if (systemCPUCount<1) {
        systemCPUCount = 1;
    }
#endif

}

//__________________________________________________________________________________
void    AddItemToPreferences (long itemCode, long itemKind, _String itemTitle, _String itemDescription,
                              _String itemValue, _List* itemOptions, _List& theTarget, bool deleteItem)
{
    if (theTarget.lLength!=6) {
        theTarget.Clear();
        _SimpleList dummySList;
        _List       dummyList;
        theTarget&& &dummySList;
        theTarget&& &dummyList;
        theTarget&& &dummyList;
        theTarget&& &dummySList;
        theTarget&& &dummyList;
        theTarget&& &dummyList;
    }

    long f = ((_List*)theTarget(1))->Find (&itemTitle);

    if (f<0) {
        (*((_SimpleList*)theTarget(0)))<<itemCode;
        (*((_SimpleList*)theTarget(3)))<<itemKind;
        (*((_List*)theTarget(1)))&& &itemTitle;
        (*((_List*)theTarget(2)))&& &itemDescription;
        (*((_List*)theTarget(4)))&& &itemValue;
        if (itemOptions) {
            (*((_List*)theTarget(5)))&& itemOptions;
        } else {
            _List dummyList;
            (*((_List*)theTarget(5)))&& & dummyList;
        }
    } else {
        (*((_SimpleList*)theTarget(0))).lData[f]=itemCode;
        (*((_SimpleList*)theTarget(3))).lData[f]=itemKind;
        (*((_List*)theTarget(1))).Replace (f,&itemTitle, true);
        (*((_List*)theTarget(2))).Replace (f,&itemDescription, true);
        (*((_List*)theTarget(4))).Replace (f,&itemValue, true);
        if (itemOptions) {
            (*((_List*)theTarget(5))).Replace (f,itemOptions, true);
        } else {
            _List dummyList;
            (*((_List*)theTarget(5))).Replace (f,&dummyList, true);
        }
    }

    if (deleteItem) {
        DeleteObject (itemOptions);
    }
}

//__________________________________________________________________________________
void    AddStringToRecentMenu (_String& fName, _String& pName)
{
    if (recentPaths.Find (&pName)>=0) {
        return;    // already in the menu - nothing to do!
    }

#if defined __HYPHY_GTK__ && ! defined __HEADLESS__
    if (!hyphyConsoleWindow) {
        return;
    }

    static    GtkItemFactoryEntry *itemFactoryHolder = nil;
    static    _List               itemFactoryStrings;
    for (long mi=0; mi<recentFiles.lLength; mi++) {
        GtkWidget * recFile = gtk_item_factory_get_widget_by_action(hyphyConsoleWindow->menu_items,2000+mi);
        if (recFile) {
            gtk_item_factory_delete_item(hyphyConsoleWindow->menu_items,gtk_item_factory_path_from_widget(recFile));
        }
    }
    if (itemFactoryHolder) {
        delete itemFactoryHolder;
    }
    itemFactoryStrings.Clear();

#endif

    if (recentFiles.lLength==RECENT_FILE_ITEMS) {
#ifdef __MAC__
        DeleteMenuItem (recentFilesMenu, RECENT_FILE_ITEMS);
#endif
        recentFiles.Delete(RECENT_FILE_ITEMS-1);
        recentPaths.Delete(RECENT_FILE_ITEMS-1);
    }
#if defined __WINDOZE__ && ! defined __HEADLESS__
    else if (recentFilesMenu) {
        InsertMenu (recentFilesMenu,-1,MF_BYPOSITION,1999+recentFiles.lLength,(LPSTR)fName.getStr());
    }
#endif

    recentFiles.InsertElement(&fName,0,true);
    recentPaths.InsertElement(&pName,0,true);

#if defined __MAC__ && ! defined __HEADLESS__
    if (recentFilesMenu) {
        Str255          buffer;
        StringToStr255 (fName, buffer);
        InsertMenuItem (recentFilesMenu,buffer,0);
    }
#endif

#if defined __WINDOZE__ && ! defined __HEADLESS__
    if (recentFilesMenu) {
        for (long mi=0; mi<recentFiles.lLength; mi++) {
            ModifyMenu (recentFilesMenu,mi+2,MF_BYPOSITION,2000+mi,((_String*)recentFiles(mi))->getStr());
        }

        DrawMenuBar ((HWND)hyphyConsoleWindow->GetOSWindowData());
    }
#endif

#if defined __HYPHY_GTK__ && ! defined __HEADLESS__
    if (recentFiles.lLength) {
        itemFactoryHolder = new GtkItemFactoryEntry[recentFiles.lLength];
    }
    for (long mi=0; mi<recentFiles.lLength; mi++) {
        itemFactoryStrings.AppendNewInstance(new _String(_String("/File/Open/Open Recent/") & *((_String*)recentFiles(mi))));
        itemFactoryHolder[mi] = (GtkItemFactoryEntry) {
            ((_String*)(itemFactoryStrings(mi)))->sData,NULL,(GtkItemFactoryCallback)hyphy_menu_item_callback,2000+mi,"<Item>"
        };
        gtk_item_factory_create_items (hyphyConsoleWindow->menu_items,  1, &(itemFactoryHolder[mi]), hyphyConsoleWindow);
    }
#endif

}

//__________________________________________________________________________________
char        AutoOpenTreeWindow (void)
{
    _Parameter ad;
    checkParameter (treeDisplayOptions, ad, 0.0);
    return     ad;
}

//__________________________________________________________________________________
void        SetShowDialogAtStartup (bool flag)
{
    if (flag) {
        *(((_String*)((_List*)globalPreferencesList.lData[4])->lData[_hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesStartupDialog))])) = "Yes";
    } else {
        *(((_String*)((_List*)globalPreferencesList.lData[4])->lData[_hyPreferencesKeys.GetXtra(_hyPreferencesKeys.Find(&_hyPreferencesStartupDialog))])) = "No";
    }
}
//EOF
