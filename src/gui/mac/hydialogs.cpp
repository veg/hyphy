
#include "stdio.h"
#include "ctype.h"
#include "string.h"
#include "batchlan.h"
#include "HYTreePanel.h"
#include "HYDataPanel.h"
#include "HYObjectInspector.h"
#include "HYUtils.h"
#ifndef  TARGET_API_MAC_CARBON
#include <StandardFile.h>
#endif
#include "Navigation.h"
#include "TextEdit.h"
#include "HYModelWindow.h"
#include "HYChartWindow.h"

#include "Movies.h"
#include <QuickTimeComponents.h>
#include <ImageCompression.h>
#include "Script.h"


extern  _String*  argFileName,
        dialogPrompt;

extern  _SimpleList
windowPtrs,
windowObjects;

extern  bool      echoPaused;
extern  long      lastMatrixDeclared;

extern  _String   objectInspectorTitle ;


_String MacSimpleFileOpen       (void);
_String MacSimpleFileSave       (void);
void    GetFullPathName         (FSSpec&, _String&);


_List               movieFormats;
_SimpleList         qtMovieGrexComponents;


NavPreviewUPP             previewFunctionHook = nil;

pascal void HYOpenEventProc (NavEventCallbackMessage,NavCBRecPtr,void * );
pascal Boolean HYOpenTextPreview (NavCBRecPtr , void *);
void  findMovieExporterComponents(_List& , _SimpleList& );


//_________________________________________________________________________
Str255      hpName = "\pHY-PHY";

pascal void HYOpenEventProc (NavEventCallbackMessage callBackSelector,
                             NavCBRecPtr callBackParms,
                             void * )
{
    if (callBackSelector==kNavCBEvent) {
        EventRecord* theEvent = callBackParms->eventData.eventDataParms.event;
        if ((theEvent->what == activateEvt)||(theEvent->what == updateEvt)) {
            if (theEvent->message != (long)callBackParms->window) {
                long k = windowPtrs.Find((long)theEvent->message);
                if (k>=0) {
                    _HYPlatformWindow* clickedWindow = (_HYPlatformWindow*)windowObjects (k);
                    clickedWindow->_ProcessOSEvent ((Ptr)theEvent);
                }
            }
        }
    }
}

//_________________________________________________________________________

pascal Boolean HYOpenTextPreview (NavCBRecPtr callBackParms, void *)
{
    AEDescList * fdesc = (AEDesc*)callBackParms->eventData.eventDataParms.param,
                 fspec;
    OSErr  navErr;

    _String    preview ("No preview available");
    if ((navErr=AECoerceDesc(fdesc,typeFSS,&fspec))==noErr) {
        FSSpec*     fSR = (FSSpec*)*(fdesc->dataHandle);
        _String feedback;
        GetFullPathName (*fSR, feedback);
        terminateExecution = false;
        //printf ("%s\n",feedback.getStr());
        FILE * F = doFileOpen (feedback.sData,"r");
        if (F) {
            char buffer [1024];
            buffer[fread (buffer,1,1023,F)]=0;
            fclose (F);
            preview = buffer;
            preview = preview.Replace("\n","\r",true);
        }
    }
    AEDisposeDesc (&fspec);
#ifdef OPAQUE_TOOLBOX_STRUCTS
    GrafPtr thisPort = GetWindowPort(callBackParms->window),
            savePort;
#else
    GrafPtr thisPort = (GrafPtr)callBackParms->window,
            savePort;
#endif
    GetPort (&savePort);
    SetPort (thisPort);
#ifdef OPAQUE_TOOLBOX_STRUCTS
    short   savedFace = GetPortTextFont (thisPort),
            savedSize = GetPortTextSize (thisPort);

    Style   savedStyle = GetPortTextFace (thisPort);
#else
    short   savedFace = thisPort->txFont,
            savedSize = thisPort->txSize;

    Style   savedStyle = thisPort->txFace;
#endif
    RGBColor saveColor,
             black = {0,0,0};
    GetForeColor (&saveColor);
    RGBForeColor (&black);
    TextFont (kFontIDHelvetica);
    TextSize (9);
    TextFace (0);
    TETextBox (preview.sData,preview.sLength,&callBackParms->previewRect,teFlushLeft);
    TextFont (savedFace);
    TextSize (savedSize);
    TextFace (savedStyle);
    RGBForeColor (&saveColor);
    SetPort (savePort);
    return   true;
}

//_________________________________________________________________________

_String MacSimpleFileOpen (void)
{
    if (previewFunctionHook == nil) {
        previewFunctionHook = NewNavPreviewUPP(HYOpenTextPreview);
    }

    _String             feedback;
    OSErr               navErr;
    NavReplyRecord      navRR;
    NavDialogOptions    navDO;
    NavEventUPP         navEF = NewNavEventUPP (HYOpenEventProc);
    NavPreviewUPP       navPF = previewFunctionHook;

    terminateExecution       = true;

    NavTypeListHandle   navLH =  (NavTypeListHandle)NewHandle(sizeof(NavTypeList) + (3 * sizeof(OSType)));

    (*navLH)->componentSignature = 'MuSe';
    (*navLH)->osTypeCount = 2;
    (*navLH)->osType[0] = '****';
    (*navLH)->osType[1] = 'TEXT';

    navDO.version       =   kNavDialogOptionsVersion;
    navDO.location      =   (Point) {
        -1,-1
    };
    navDO.dialogOptionFlags
        =   kNavAllFilesInPopup|kNavSelectAllReadableItem|kNavAllowPreviews|kNavDontAutoTranslate;
    StringToStr255      (dialogPrompt,navDO.windowTitle);
    memcpy (navDO.clientName,hpName,hpName[0]+1);
    navDO.actionButtonLabel[0]
        =   0;
    navDO.cancelButtonLabel[0]
        =   0;
    navDO.message[0]    =   0;
    navDO.preferenceKey =   0;
    navDO.popupExtension=   nil;

    navErr = NavGetFile (nil,
                         &navRR,
                         &navDO,
                         navEF,
                         navPF,
                         nil,
                         navLH,
                         nil);

    if (navErr == noErr) {
        if (navRR.validRecord) {
            long countAED;
            if (noErr==AECountItems(&navRR.selection,&countAED)) {
                if (countAED==1) {
                    char    fileRec [2048];
                    Size    actualSize;
                    AEKeyword   keywd;
                    DescType    returnedType;
                    FSSpec*     fSR;
                    navErr= AEGetNthPtr(&navRR.selection, 1, typeFSS, &keywd,
                                        &returnedType, fileRec,
                                        sizeof(FSSpec), &actualSize);
                    if (navErr==noErr) {
                        fSR = (FSSpec*)fileRec;
                        GetFullPathName (*fSR, feedback);
                        terminateExecution = false;
                    }
                }
            }
        }
        NavDisposeReply (&navRR);
    }
#ifdef TARGET_API_MAC_CARBON
    DisposeNavEventUPP (navEF);
#endif
    DisposeHandle ((Handle)navLH);
    return      feedback;
}

//_________________________________________________________________________

_String ChooseAFolder (_String& promptString)
{
    _String             feedback;
    OSErr               navErr;
    NavReplyRecord      navRR;
    NavDialogOptions    navDO;
#ifdef TARGET_API_MAC_CARBON
    NavEventUPP         navEF = NewNavEventUPP (HYOpenEventProc);
#else
    NavEventUPP         navEF = NewNavEventProc (HYOpenEventProc);
#endif

    terminateExecution       = true;

    navDO.version       =   kNavDialogOptionsVersion;
    navDO.location      =   (Point) {
        -1,-1
    };
    navDO.dialogOptionFlags
        =   kNavAllFilesInPopup|kNavSelectAllReadableItem|kNavAllowPreviews;
    StringToStr255      (promptString,navDO.windowTitle);
    memcpy (navDO.clientName,hpName,hpName[0]+1);

    navDO.actionButtonLabel[0]
        =   0;
    navDO.cancelButtonLabel[0]
        =   0;
    navDO.message[0]    =   0;
    navDO.preferenceKey =   0;
    navDO.popupExtension=   nil;

    navErr = NavChooseFolder (nil,
                              &navRR,
                              &navDO,
                              navEF,
                              nil,
                              nil);

    if (navErr == noErr) {
        if (navRR.validRecord) {
            long countAED;
            if (noErr==AECountItems(&navRR.selection,&countAED)) {
                if (countAED==1) {
                    char    fileRec [2048];
                    Size    actualSize;
                    AEKeyword   keywd;
                    DescType    returnedType;
                    FSSpec*     fSR;
                    navErr= AEGetNthPtr(&navRR.selection, 1, typeFSS, &keywd,
                                        &returnedType, fileRec,
                                        sizeof(FSSpec), &actualSize);
                    if (navErr==noErr) {
                        fSR = (FSSpec*)fileRec;
                        GetFullPathName (*fSR, feedback);
                        terminateExecution = false;
                    }
                }
            }
        }
        NavDisposeReply (&navRR);
    }
#ifdef TARGET_API_MAC_CARBON
    DisposeNavEventUPP (navEF);
#endif
    return      feedback;
}

//_________________________________________________________________________

_String MacSimpleFileSave (void)
{
    _String             feedback;
    OSErr               navErr;
    NavReplyRecord      navRR;
    NavDialogOptions    navDO;
#ifdef TARGET_API_MAC_CARBON
    NavEventUPP         navEF = NewNavEventUPP (HYOpenEventProc);
#else
    NavEventUPP         navEF = NewNavEventProc (HYOpenEventProc);
#endif

    terminateExecution       = true;

    navDO.version       =   kNavDialogOptionsVersion;
    navDO.location      =   (Point) {
        -1,-1
    };
    navDO.dialogOptionFlags
        =   kNavAllFilesInPopup|kNavSelectAllReadableItem|kNavNoTypePopup;

    if (defFileNameValue.sLength) {
        StringToStr255 (defFileNameValue, navDO.savedFileName);
    } else {
        navDO.savedFileName[0]= 0;
    }

    StringToStr255      (dialogPrompt,navDO.windowTitle);
    memcpy (navDO.clientName,hpName,hpName[0]+1);
    navDO.actionButtonLabel[0]
        =   0;
    navDO.cancelButtonLabel[0]
        =   0;
    navDO.message[0]    =   0;
    navDO.preferenceKey =   0;
    navDO.popupExtension=   nil;

    navErr = NavPutFile (nil,
                         &navRR,
                         &navDO,
                         navEF,
                         'TEXT',
                         'MuSe',
                         nil);

    if (navErr == noErr) {
        if (navRR.validRecord) {
            long countAED;
            if (noErr==AECountItems(&navRR.selection,&countAED)) {
                if (countAED==1) {
                    char    fileRec [2048];
                    Size    actualSize;
                    AEKeyword   keywd;
                    DescType    returnedType;
                    FSSpec*     fSR;
                    navErr= AEGetNthPtr(&navRR.selection, 1, typeFSS, &keywd,
                                        &returnedType, fileRec,
                                        sizeof(FSSpec), &actualSize);
                    if (navErr==noErr) {
                        fSR = (FSSpec*)fileRec;
                        GetFullPathName (*fSR, feedback);
                        terminateExecution = false;
                    }
                }
            }
        }
        NavDisposeReply (&navRR);
    }
#ifdef TARGET_API_MAC_CARBON
    DisposeNavEventUPP (navEF);
#endif
    return      feedback;
}

//_________________________________________________________________________
bool    PopUpFileDialog(_String ps, _String* defaultLocation)
{
    Str255 promptS;

    StringToStr255 (ps, promptS);

    previewFunctionHook = NewNavPreviewUPP(HYOpenTextPreview);
    OSErr               navErr;
    NavReplyRecord      navRR;
    NavDialogOptions    navDO;
    NavEventUPP         navEF = NewNavEventUPP (HYOpenEventProc);
    NavPreviewUPP       navPF = previewFunctionHook;

    NavTypeListHandle   navLH =  (NavTypeListHandle)NewHandle(sizeof(NavTypeList) + (3 * sizeof(OSType)));

    (*navLH)->componentSignature = kNavGenericSignature;
    (*navLH)->osTypeCount        = 3;
    (*navLH)->osType[0]          = kNavGenericSignature;
    (*navLH)->osType[1]          = 'TEXT';
    (*navLH)->osType[2]          = 'text';

    navDO.version       =   kNavDialogOptionsVersion;
    navDO.location      =   (Point) {
        -1,-1
    };
    navDO.dialogOptionFlags
        =   kNavAllFilesInPopup|kNavAllowPreviews|kNavSelectAllReadableItem|kNavDontAutoTranslate;
    memcpy (navDO.windowTitle,promptS,promptS[0]+1);
    memcpy (navDO.clientName,hpName,hpName[0]+1);
    navDO.actionButtonLabel[0]
        =   0;
    navDO.cancelButtonLabel[0]
        =   0;
    navDO.message[0]    =   0;
    navDO.preferenceKey =   0;
    navDO.popupExtension=   nil;

    FSSpec              defLocFSS;

    if (defaultLocation) {
        Str255 buffer;
        StringToStr255 (*defaultLocation, buffer);
        FSMakeFSSpec (0,0,buffer,&defLocFSS);
        AEDesc defDesc;
        AECreateDesc (typeFSS,&defLocFSS,sizeof (FSSpec),&defDesc);
        navErr = NavGetFile (&defDesc,
                             &navRR,
                             &navDO,
                             navEF,
                             navPF,
                             nil,
                             navLH,
                             nil);
        AEDisposeDesc (&defDesc);
    } else {
        /*  NavDialogRef           daBox;
            navErr = NavCreateGetFileDialog (
                                    &navDO,
                                    &navRR,
                                    navLH,
                                    navEF,
                                    navPF,
                                    NULL,
                                    NULL,
                                    &daBox
                                   );*/

        navErr = NavGetFile (nil,
                             &navRR,
                             &navDO,
                             navEF,
                             navPF,
                             nil,
                             navLH,
                             nil);
    }

    if (navErr == noErr) {
        if (navRR.validRecord) {
            long countAED;
            if (noErr==AECountItems(&navRR.selection,&countAED)) {
                if (countAED==1) {
                    char        fileRec [2048];
                    Size        actualSize;
                    AEKeyword   keywd;
                    DescType    returnedType;
                    FSSpec*     fSR;
                    // get the direct parameter--a descriptor list--and put
                    // it into docList
                    navErr= AEGetNthPtr(&navRR.selection, 1, typeFSS, &keywd,
                                        &returnedType, fileRec,
                                        sizeof(FSSpec), &actualSize);
                    if (navErr==noErr) {
                        fSR = (FSSpec*)fileRec;
                        fSR->name [fSR->name[0]+1]=0;
                        if (!argFileName) {
                            argFileName = new _String ((char*)(fSR->name+1));
                        } else {
                            *argFileName = _String ((char*)(fSR->name+1));
                        }

                        long       parentDirID = fSR->parID;
                        CInfoPBRec infoRec;
                        HFileInfo* accessInfo = (HFileInfo*)&infoRec;
                        Str63      fName;
                        accessInfo->ioVRefNum = fSR->vRefNum;
                        accessInfo->ioNamePtr = fName;
                        accessInfo->ioFDirIndex = -1;


                        /*FSRefParam             fileBlock;
                        FSCatalogInfo        catInfo;
                        fileBlock.ref        = fSR;
                        fileBlock.whichInfo  =  kFSCatInfoFinderInfo |  kFSCatInfoFinderXInfo;
                        fileBlock.spec       = nil;
                        fileBlock.catInfo    = &catInfo;

                        PBGetCatalogInfoSync (&fileBlock);*/

                        while (parentDirID!=fsRtParID) {
                            accessInfo->ioDirID = parentDirID;
                            if (PBGetCatInfo (&infoRec,false)) {
                                NavDisposeReply (&navRR);
#ifdef TARGET_API_MAC_CARBON
                                DisposeNavEventUPP (navEF);
#endif
                                DisposeHandle ((Handle)navLH);
                                return false;
                            }
                            parentDirID = accessInfo->ioFlParID;
                            accessInfo->ioNamePtr [accessInfo->ioNamePtr[0]+1]=0;
                            *argFileName = _String(((char*)(accessInfo->ioNamePtr+1)))&_String(":")&(*argFileName);
                        }
                        volumeName =  _String(((char*)(accessInfo->ioNamePtr+1))) & ':';
                        NavDisposeReply (&navRR);
#ifdef TARGET_API_MAC_CARBON
                        DisposeNavEventUPP (navEF);
#endif
                        return true;
                    }
                }
            }
        }
    }
    NavDisposeReply (&navRR);
#ifdef TARGET_API_MAC_CARBON
    DisposeNavEventUPP (navEF);
#endif
    DisposeHandle ((Handle)navLH);
    return false;
}

//____________________________________________________________________________________

void    GetFullPathName (FSSpec& theReply, _String& feedback)
{
    theReply.name [theReply.name[0]+1]=0;
    feedback= _String(((char*)(theReply.name+1)));
    // get full path name
    long parentDirID = theReply.parID;
    CInfoPBRec infoRec;
    HFileInfo* accessInfo = (HFileInfo*)&infoRec;
    Str63 fName;
    accessInfo->ioVRefNum = theReply.vRefNum;
    accessInfo->ioNamePtr = fName;
    accessInfo->ioFDirIndex = -1;
    while (parentDirID!=fsRtParID) {
        accessInfo->ioDirID = parentDirID;
        if (PBGetCatInfo (&infoRec,false)!=noErr) {
            //_String warnMsg("Error in file prompt dialog.");
            //acknError (warnMsg);
            //return      feedback;
            feedback = empty;
            return;
        }
        parentDirID = accessInfo->ioFlParID;
        accessInfo->ioNamePtr [accessInfo->ioNamePtr[0]+1]=0;
        feedback = _String(((char*)(accessInfo->ioNamePtr+1)))&_String(":")&(feedback);
    }
    volumeName =  _String(((char*)(accessInfo->ioNamePtr+1))) & ':';
}



//_________________________________________________________________________
void findMovieExporterComponents(_List& compList, _SimpleList& compIndex)
{
    ComponentDescription cd, cd2;
    Component c = 0;

    cd.componentType            = MovieExportType;
    cd.componentSubType         = 0;
    cd.componentManufacturer    = 0;
    cd.componentFlags           = 0;
    cd.componentFlagsMask       = 0;

    _String fileFormat;

    while( ( c = FindNextComponent( c, &cd ) ) != 0 ) {
        Handle     cInfo = NewHandle(256);
        GetComponentInfo (c,&cd2,cInfo,nil,nil);
        (*cInfo)[**cInfo+1] = 0;
        fileFormat = (char*)(*cInfo+1);
        if (fileFormat.sLength) {
            compList&& &fileFormat;
            compIndex << (long)c;
        }
        DisposeHandle(cInfo);
    }
}

//_________________________________________________________________________
void    ConvertMovieFile (bool truncate)
{
    EnterMovies();
    OSErr nErr;

    Movie       movie = nil;
    short       nFileRefNum;
    char        fileRec [2048];
    FSSpec*     fSR = nil;
    Str255      promptS = "\pChoose a movie file to convert";

    OSErr               navErr;
    NavReplyRecord      navRR;
    NavDialogOptions    navDO;
#ifdef TARGET_API_MAC_CARBON
    NavEventUPP         navEF = NewNavEventUPP (HYOpenEventProc);
#else
    NavEventUPP         navEF = NewNavEventProc (HYOpenEventProc);
#endif

    NavPreviewUPP       navPF = nil;
    NavTypeList         navTL;
    NavTypeListPtr      navTP = &navTL;
    NavTypeListHandle   navLH = &navTP;

    navTL.componentSignature = 'MuSe';
    navTL.osTypeCount        = 1;
    navTL.osType[0]          = MovieFileType;

    navDO.version       =   kNavDialogOptionsVersion;
    navDO.location      =   (Point) {
        -1,-1
    };
    navDO.dialogOptionFlags
        =   kNavAllFilesInPopup|kNavSelectAllReadableItem|kNavAllowPreviews|kNavDontAutoTranslate;

    memcpy (navDO.windowTitle,promptS,promptS[0]+1);
    memcpy (navDO.clientName,hpName,hpName[0]+1);
    navDO.actionButtonLabel[0]
        =   0;
    navDO.cancelButtonLabel[0]
        =   0;
    navDO.message[0]    =   0;
    navDO.preferenceKey =   0;
    navDO.popupExtension=   nil;

    navErr = NavGetFile (nil,
                         &navRR,
                         &navDO,
                         navEF,
                         navPF,
                         nil,
                         navLH,
                         nil);

    if (navErr == noErr) {
        if (navRR.validRecord) {
            long countAED;
            if (noErr==AECountItems(&navRR.selection,&countAED)) {
                if (countAED==1) {
                    Size    actualSize;
                    AEKeyword   keywd;
                    DescType    returnedType;
                    // get the direct parameter--a descriptor list--and put
                    // it into docList
                    navErr= AEGetNthPtr(&navRR.selection, 1, typeFSS, &keywd,
                                        &returnedType, fileRec,
                                        sizeof(FSSpec), &actualSize);
                    if (navErr==noErr) {
                        fSR = (FSSpec*)fileRec;
                    } else {
                        WarnError (_String("Error ")& (long)navErr & " in AEGetNthPtr.");
                    }
                }
            } else {
                WarnError (_String("Error ")& (long)navErr & " in AECountItems.");
            }
        }
    } else {
        WarnError (_String("Error ")& (long)navErr & " in NavGetFile.");
    }


    if (fSR) {
        nErr = OpenMovieFile(fSR, &nFileRefNum, fsRdPerm);
        if (nErr == noErr) {
            short nResID = 0;
            Str255 strName;
            Boolean bWasChanged;
            nErr = NewMovieFromFile(&movie, nFileRefNum, &nResID, strName,newMovieActive, &bWasChanged);
            SetMovieProgressProc(movie, (MovieProgressUPP)-1L, 0);

            if (qtMovieGrexComponents.lLength==0) {
                findMovieExporterComponents (movieFormats,qtMovieGrexComponents);
            }

            if (qtMovieGrexComponents.lLength==0) {
                WarnError (_String("Failed to find movie export components."));
            } else {
                _String fName,
                        filePr   = "Export Movie To:",
                        formatPr = "Target Format:",
                        movieName(10L,true);

                for (long k = 1; k<=fSR->name[0]; k++) {
                    movieName << fSR->name[k];
                }

                movieName.Finalize();

                long    menuChoice = SaveFileWithPopUp (fName, filePr,movieName,formatPr,movieFormats);

                if (menuChoice>=0) {
                    FSSpec  fs;
                    Str255  buff;
                    StringToStr255 (fName,buff);
                    FSMakeFSSpec(0,0,buff,&fs);


                    ComponentInstance grexc = OpenComponent ((Component)qtMovieGrexComponents(menuChoice));

                    Boolean canceled;

                    if (truncate) {
                        TimeValue maxDuration = GetMovieDuration(movie);
                        _String chopString,
                                chopPrompt = _String ("From-to, max value = ") & (long)maxDuration;

                        canceled = 1;
                        if (EnterStringDialog (chopString, chopPrompt, nil)) {
                            _List * times = chopString.Tokenize("-");
                            if (times->lLength==2) {
                                TimeValue fromMark = ((_String*)(*times)(0))->toNum(),
                                          toMark = ((_String*)(*times)(1))->toNum();

                                if ((fromMark>=0)&&(toMark>fromMark)&&(toMark<=maxDuration)) {
                                    MovieExportDoUserDialog (grexc, movie, nil, fromMark,toMark-fromMark, & canceled);
                                    canceled = 0;
                                }
                                if (canceled) {
                                    WarnError (_String((long)fromMark)&"-"&_String((long)toMark) & " is not a valid segment specification.");
                                }
                            }
                        }
                    } else {
                        MovieExportDoUserDialog (grexc, movie, nil, 0, GetMovieDuration(movie), & canceled);
                    }

                    if (!canceled) {
                        nErr =  ConvertMovieToFile (movie,0,&fs,'    ','MuSe',smSystemScript,nil,0,grexc);
                        if ((nErr != noErr)&&(nErr != progressProcAborted)) {
                            WarnError (_String("Error ")& (long)nErr & " converting movie to file.");
                        }
                    }
                    CloseComponent (grexc);
                }
            }
            CloseMovieFile(nFileRefNum);

        } else {
            WarnError (_String("Error ")& (long)nErr & " opening movie file.");
        }
    }
    NavDisposeReply (&navRR);
    ExitMovies();
}

//_________________________________________________________________________

void    ShowObjectInspector (void)
{
    long f = FindWindowByName (objectInspectorTitle);
    if (f>=0) {
        ShowWindow ((WindowPtr)windowPtrs (f));
        SelectWindow ((WindowPtr)windowPtrs (f));
    } else {
        _HYObjectInspector* newOI = new _HYObjectInspector ();
        //newOI->BuildListOfObjects (0);
        //newOI->_Zoom(true);
        newOI->Activate       ( );
        //newOI->Show();
    }
}

//_________________________________________________________________________

_String     DoMacToPOSIX (const _String& in)
{
    CFStringRef posixPath    = CFStringCreateWithCString (NULL,in.sData,kCFStringEncodingASCII);
    CFURLRef    convertorURL = CFURLCreateWithFileSystemPath (NULL,posixPath,kCFURLHFSPathStyle,false);
    CFRelease (posixPath);
    posixPath = CFURLCopyFileSystemPath (convertorURL, kCFURLPOSIXPathStyle);
    CFRelease (convertorURL);
    _String  newFNAME  (CFStringGetLength(posixPath),false);
    CFStringGetCString (posixPath, newFNAME.sData, newFNAME.sLength+1, kCFStringEncodingASCII);
    CFRelease (posixPath);
    return newFNAME;
}





