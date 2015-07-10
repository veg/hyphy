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

#ifndef __HYPHYXCODE__
#include <unix.h>
#endif

#ifndef  TARGET_API_MAC_CARBON
#include <Types.h>
#else
#include <MacTextEditor.h>
#endif


#include <events.h>
#include <AppleEvents.h>
#include <Aliases.h>
#include <Dialogs.h>
#include <Files.h>
#include <Script.h>
#include <unistd.h>
#include <Menus.h>
#include <stdlib.h>
#include <Fonts.h>
#include <Resources.h>
#include <TextUtils.h>
#include <Sound.h>
#include <Folders.h>
#include <ToolUtils.h>
#include <resources.h>
#include "likefunc.h"
#include <string.h>
#include <ctype.h>
#include <Memory.h>
#include <Gestalt.h>
#include <Appearance.h>
#include <CarbonEvents.h>
#include "Lists.h"
#include <time.h>
#include "HYWindow.h"
#include "HYTreePanel.h"
#include "HYDataPanel.h"
#include "HYUtils.h"
#include "HYDialogs.h"

#include "preferences.h"

#include "HYSharedMain.h"
#include "HYButtonBar.h"
#include "HYButton.h"

#include "HYConsoleWindow.h"
#include "HYEventTypes.h"
#include "HYTextBox.h"
#include "HYObjectInspector.h"

#include "HYDBWindow.h"

#ifdef   HYPHY_MACH
#include <Carbon/Carbon.h>
#include <console.h>
systemCPUCount = (long)MPProcessors();
#else
#include <windows.h>
#endif


#ifdef  __MACPROFILE__
#include "profiler.h"
#endif


extern long     matrixExpCount,
       taylorTermsCount,
       squaringsCount,
       likeFuncEvalCallCount,
       _fcreator,
       _ftype;



bool            aquaInterfaceOn = false;

OSType          fileCreatorType = 'MuSe';

Pattern         penHatchPattern,
                vertPenHatchPattern;


void    DisplayAbout                (bool = false);
void    HandleSuspend               (void);
void    HandleResume                (void);
void    HandleCancel                (void);
void    ShowMessagesLog             (void);
long    SelectATemplate             (void);
void    RunTemplate                 (long);
void    SetupMenusAndSuch           (void);
void    updateTimerF                (_String&, time_t);
void    SpoolFile                   (void);
bool    HandleFirstDialog           (long& choice);
void    RunStandardAnalyses         (void);
void    NewModel                    (_String*);
void    BuildBalancedTree           (long, long , _String& , _List& );
void    SetStatusBarValue           (long,_Parameter, _Parameter);
pascal  OSErr GUIOpenHandler        (AppleEvent *, AppleEvent*, long);
OSErr   HandleQuit                  (AppleEvent *, AppleEvent*, long);
void    InitToolbox                 (void);
void    CleanupMenusAndSuch         (void);
void    GetFullPathName             (FSSpec&,_String&);


pascal  void        TimerAction                     (EventLoopTimerRef,void*);
pascal  void        DialogTimerAnimateAction        (EventLoopTimerRef,void*);
pascal  void        DialogTimerAnimateCornholio     (EventLoopTimerRef,void*);
pascal OSStatus     PrefAction                      (EventHandlerCallRef, EventRef event, void*);




PixPatHandle    whiteFill,
                blueFill,
                statusBarFill;

Handle          menuBar;

MenuHandle      recentFilesMenu   = NewMenu (200, "\pRun Recent");

time_t          timerStart = 0,
                lastTimer  = 0;

bool            trackMouseMovement = false;

CIconHandle     redButtonIcon,
                yellowButtonIcon,
                orangeButtonIcon,
                greenButtonIcon,
                pullDownArrowsIcon,
                tablePDMenuIcon;

Cursor          hSizeCursor,
                editStateCursor;

CCrsrHandle     pickUpCursor,
                dropOffCursor;


//____________________________________________________________________________________________

PMPageFormat            gPageFormat         = kPMNoPageFormat;
PMPrintSettings         gPrintSettings      = kPMNoPrintSettings;
Handle                  gFlattenedFormat    = NULL;
Handle                  gFlattenedSettings  = NULL;

//____________________________________________________________________________________________

Boolean InitPrint(PMPrintSession theSession)
{
    OSStatus theStatus;

    if (gFlattenedFormat != NULL) {
        theStatus = PMUnflattenPageFormat(gFlattenedFormat, &gPageFormat);

        if (theStatus != noErr) {
            return FALSE;
        }
    }

    if (gFlattenedSettings != NULL) {
        theStatus = PMUnflattenPrintSettings(gFlattenedSettings, &gPrintSettings);

        if (theStatus != noErr) {
            return FALSE;
        }
    }

    if (gPageFormat != kPMNoPageFormat) {
        theStatus = PMSessionValidatePageFormat(theSession,gPageFormat, kPMDontWantBoolean);
    } else {
        theStatus = PMCreatePageFormat(&gPageFormat);

        if ((theStatus == noErr) && (gPageFormat != kPMNoPageFormat)) {
            theStatus = PMSessionDefaultPageFormat(theSession, gPageFormat);
        }
    }

    return (theStatus == noErr);
}

//____________________________________________________________________________________________

RgnHandle       mouseRgn = nil;

pascal void TimerAction (EventLoopTimerRef,void*)
{
    if (updateTimer) {
        time_t  tt,
                time_diff = time(&tt)-lastTimer;

        if (time_diff>=1) { // update dispay
            lastTimer   = tt;

            _String      r;
            updateTimerF(r,tt-timerStart);
            SetStatusLine (empty, empty, r, 0, HY_SL_TIMER);
        }
    }
}

//____________________________________________________________________________________________

struct conrholioTimeScrollData {
    PicHandle  corn1,
               corn2;

    Rect       cornRect,
               bigRect;

    long       callCount;

    DialogPtr  AW;
};


//____________________________________________________________________________________________

struct textScrollTimerAction {
    DialogPtr  AW;

    bool       fadeout,
               morePages;

    long       decrement,
               startingString,
               stringDisp,
               textTime,
               color,
               shutter,
               callCount;

    Rect       textBoxInset,
               textBoxSlide;

    _List*     theText;

};

//____________________________________________________________________________________________

pascal void DialogTimerAnimateAction (EventLoopTimerRef,void*userData)
{
    textScrollTimerAction * tSTA = (textScrollTimerAction*)userData;

    GrafPtr   sPort;
    GetPort  (&sPort);
    SetPort  (GetDialogPort (tSTA->AW));
    Str255    buffer;

    RGBColor  myColor;

    if (!tSTA->fadeout) {
        if (tSTA->color >= 0) {
            myColor.red     = tSTA->color;
            myColor.blue    = tSTA->color;
            myColor.green   = tSTA->color;
            RGBForeColor(&myColor);
            tSTA->color -= tSTA->decrement;
            tSTA->decrement+=100;

            long i,
                 //width = tSTA->textBoxInset.right - tSTA->textBoxInset.left,
                 finalTick = tSTA->textBoxInset.top + 11;

            for (i = tSTA->startingString; i<tSTA->theText->lLength; i++, finalTick += 11) {
                if (finalTick>tSTA->textBoxInset.bottom) {
                    break;
                }

                _String* thisString = (_String*)(*tSTA->theText)(i);

                buffer[0] =thisString ->sLength;
                if (buffer[0]&& (thisString->sData[0]=='\\')) { // new page
                    if (i!=tSTA->startingString) {
                        tSTA->morePages = true;
                        tSTA->stringDisp = i;
                        break;
                    }
                    if (buffer[0]&& (thisString->sData[1]=='_')) {
                        TextFace(1);
                        buffer[0]-=2;
                        memcpy (buffer+1,thisString->sData+2,buffer[0]);
                    } else {
                        TextFace(0);
                        buffer[0]--;
                        memcpy (buffer+1,thisString->sData+1,buffer[0]);
                    }
                } else if (buffer[0]&& (thisString->sData[0]=='_')) {
                    TextFace(1);
                    buffer[0]--;
                    memcpy (buffer+1,thisString->sData+1,buffer[0]);
                } else {
                    TextFace(0);
                    memcpy (buffer+1,thisString->sData,buffer[0]);
                }
                //MoveTo ((width-TextWidth ((Ptr)buffer,1,buffer[0]))/2+textBoxInset.left,finalTick);
                MoveTo (tSTA->textBoxInset.left+5,finalTick);
                DrawString (buffer);
            }

            if ((i==tSTA->theText->lLength)&&(tSTA->color<=0)) {
                tSTA->morePages = false;
                tSTA->startingString = -1;
            }

        } else {
            if (tSTA->morePages) {
                tSTA->startingString = tSTA->stringDisp;
            }

            tSTA->color = 0xddff;
            tSTA->fadeout = true;
            tSTA->decrement = 100;
            tSTA->callCount = 1;
            tSTA->textTime = 0;
            if (tSTA->startingString==-1) {
                tSTA->startingString = 0;
            }
        }
    } else {
        if (tSTA->callCount%100 == 0) {
            if (tSTA->shutter == 0) {
                if (tSTA->textTime<tSTA->textBoxSlide.bottom-tSTA->textBoxSlide.top) {
                    tSTA->textTime+=4;

                    Rect    eraser = tSTA->textBoxSlide;

                    eraser.bottom = eraser.top + tSTA->textTime/2;
                    EraseRect (&eraser);
                    eraser.bottom = tSTA->textBoxSlide.bottom;
                    eraser.top = eraser.bottom - tSTA->textTime/2;
                    EraseRect (&eraser);
                } else {
                    tSTA->shutter = (tSTA->shutter+1)%3;
                    tSTA->fadeout = false;
                }
            } else if (tSTA->shutter == 1) {
                if (tSTA->textTime<tSTA->textBoxSlide.right-tSTA->textBoxSlide.left) {
                    tSTA->textTime+=4;

                    Rect    eraser = tSTA->textBoxSlide;

                    eraser.right = eraser.left + tSTA->textTime/2;
                    EraseRect (&eraser);
                    eraser.right = tSTA->textBoxSlide.right;
                    eraser.left = eraser.right - tSTA->textTime/2;
                    EraseRect (&eraser);
                } else {
                    tSTA->shutter = (tSTA->shutter+1)%3;
                    tSTA->fadeout = false;
                }
            } else {
                if (tSTA->textTime<tSTA->textBoxSlide.right-tSTA->textBoxSlide.left) {
                    tSTA->textTime+=4;

                    Rect    eraser = tSTA->textBoxSlide;

                    eraser.right = eraser.left + tSTA->textTime;
                    EraseRect (&eraser);
                    eraser = tSTA->textBoxSlide;
                    eraser.bottom = eraser.top + tSTA->textTime;
                    if (eraser.bottom<=tSTA->textBoxSlide.bottom) {
                        EraseRect (&eraser);
                    }
                } else {
                    tSTA->shutter = (tSTA->shutter+1)%3;
                    tSTA->fadeout = false;
                }
            }

        } else {
            tSTA->callCount ++;
        }

    }
    SetPort(sPort);
}

//____________________________________________________________________________________________

pascal void DialogTimerAnimateCornholio (EventLoopTimerRef,void* userData)
{
    conrholioTimeScrollData * tSTA = (conrholioTimeScrollData*)userData;

    GrafPtr   sPort;
    GetPort  (&sPort);
    SetPort  (GetDialogPort (tSTA->AW));

    if (tSTA->callCount%2) {
        DrawPicture (tSTA->corn2, &tSTA->cornRect);
    } else {
        DrawPicture (tSTA->corn1, &tSTA->cornRect);
    }

    tSTA->cornRect.left += 3;
    tSTA->cornRect.right+= 3;

    if (tSTA->cornRect.left>tSTA->bigRect.right) {
        long w = tSTA->cornRect.right-tSTA->cornRect.left;
        long h = tSTA->cornRect.bottom-tSTA->cornRect.top;
        tSTA->cornRect.left = tSTA->bigRect.left;
        tSTA->cornRect.right = tSTA->cornRect.left+w;
        tSTA->cornRect.top+=h;
        tSTA->cornRect.bottom+=h;
    }

    if (tSTA->cornRect.bottom>tSTA->bigRect.bottom) {
        long h = tSTA->cornRect.bottom-tSTA->cornRect.top;
        tSTA->cornRect.top = tSTA->bigRect.top;
        tSTA->cornRect.bottom = tSTA->cornRect.top + h;
    }

    tSTA->callCount++;

    SetPort(sPort);
}

//____________________________________________________________________________________________

pascal OSStatus PrefAction (EventHandlerCallRef, EventRef event, void*)
{
    HICommand commandStruct;
    GetEventParameter (event, kEventParamDirectObject,typeHICommand, NULL, sizeof(HICommand),NULL, &commandStruct);
    if (commandStruct.commandID ==  kHICommandPreferences) {
        HandlePreferences (globalPreferencesList,"HYPHY Preferences");
        return noErr;
    }
    return eventNotHandledErr;
}

//__________________________________________________________________________________

pascal OSErr GUIOpenHandler (AppleEvent *guiMessage, AppleEvent* , long )
{
    char        aliasRec [8192];
    AEDescList  docList;
    OSErr       myErr;
    long        itemsInList;
    Size        actualSize;
    AEKeyword   keywd;
    DescType    returnedType;
    FSSpec      resolvedFile;
    unsigned char       hasChanged;
    // get the direct parameter--a descriptor list--and put
    // it into docList
    myErr = AEGetParamDesc(guiMessage, keyDirectObject,
                           typeAEList, &docList);
    if (myErr != noErr) {
        return myErr;
    }

    // count the number of descriptor records in the list
    myErr = AECountItems (&docList,&itemsInList);

    // now get each descriptor record from the list, coerce
    // the returned data to an FSSpec record, and open the
    // associated file
    if (itemsInList>=1) {
        for (long daItem = 1; daItem <= itemsInList; daItem++) {
            myErr = AEGetNthPtr(&docList, daItem, typeAlias, &keywd,
                                &returnedType, aliasRec,
                                sizeof(aliasRec), &actualSize);
            if (myErr != noErr) {
                return myErr;
            }

            AliasPtr fileAlias = (AliasPtr)aliasRec;
            myErr = ResolveAlias((FSSpec*)NULL, &fileAlias,&resolvedFile, &hasChanged);
            if (myErr) {
                return myErr;
            }
            resolvedFile.name [resolvedFile.name[0]+1]=0;
            if (!argFileName) {
                argFileName = new _String ((char*)(resolvedFile.name+1));
            } else {
                *argFileName =  _String ((char*)(resolvedFile.name+1));
            }
            // get full path name
            long parentDirID = resolvedFile.parID, sPDID = parentDirID;
            CInfoPBRec infoRec;
            HFileInfo* accessInfo = (HFileInfo*)&infoRec;
            Str63 fName;
            accessInfo->ioVRefNum = resolvedFile.vRefNum;
            accessInfo->ioNamePtr = fName;
            accessInfo->ioFDirIndex = -1;
            while (parentDirID!=fsRtParID) {
                accessInfo->ioDirID = parentDirID;
                if (PBGetCatInfo (&infoRec,false)) {
                    return -43;
                }
                parentDirID = accessInfo->ioFlParID;
                accessInfo->ioNamePtr [accessInfo->ioNamePtr[0]+1]=0;
                *argFileName = _String(((char*)(accessInfo->ioNamePtr+1)))&_String(":")&(*argFileName);
            }
            volumeName =  _String(((char*)(accessInfo->ioNamePtr+1))) & ':';
            FInfo thisFile;
            myErr = HGetFInfo(resolvedFile.vRefNum,sPDID,resolvedFile.name,&thisFile);
            if (thisFile.fdCreator==fileCreatorType || argFileName->endswith (".bf")) {
                OpenBatchFile(false);
                ExecuteBatchFile();
            } else {
                SpoolFile();
            }
        }
        return  noErr;
    }
    return -1;
}

//__________________________________________________________________________________

OSErr HandleQuit (AppleEvent *, AppleEvent*, long)
{
    highLevelQuit = true;
    terminateExecution = true;
    return  noErr;
}

//____________________________________________________________________________________________

void DisplayAbout (bool splash)
{


    DialogPtr       AW;
    EventRecord     theEvent;
    GrafPtr         sPort;
    WindowPtr       whatWindow;

    short           iType,
                    counter;

    Rect            textBox;

    RgnHandle       borderRegion;

    Handle          iHandle;
    _List           theText;
    _String         theString, theMessage;
    Str255          buffer;


    conrholioTimeScrollData cTSD;



    counter = 1;
    GetIndString(buffer,9000+splash,counter);
    do {
        counter ++;
        buffer[buffer[0]+1]=0;
        if (buffer[0]) {
            theString = _String((char*)buffer+1);
            //if (splash&&(theString.sData[0]=='\\')) break;
            theText&&(&theString);
        }
        GetIndString(buffer,9000+splash,counter);
    } while (buffer[0]);

    Handle         sndHandle = nil;
    SndChannelPtr  sndChn = nil;

    char           whichOne = 0;

    if (!splash) {
        theString = "\\_System Tidbits";
        theText&&(&theString);
        theString = _String("Free HYPHY Memory: ")&_String(FreeMem()/1024)&_String (" K");
        theText&&(&theString);
        if (!aquaInterfaceOn) {
            theString = _String("Free Temporary Memory: ")&_String(TempFreeMem()/1024)&_String (" K");
            theText&&(&theString);
        }
        long    gestResp;
        Gestalt (gestaltProcClkSpeed, &gestResp);
        theString = _String("Processor Speed: ")&_String(gestResp/1000000.0)&" MHz";
        theText&&(&theString);
#ifdef __MP__
        theString = _String("Processor Count: ")&_String((long)MPProcessors());
        theText&&(&theString);
#endif
        theString = MatrixExpCounter();
        if (theString.sLength) {
            theText&&(&theString);
            counter+=5;
        }
    } else {
        whichOne = (genrand_real2()>=.5);
        sndHandle = GetResource ('snd ',whichOne?128:129);
        if (sndHandle) {
            HLock (sndHandle);
            SndNewChannel (&sndChn, 0, 0, nil);
            SndPlay (sndChn, (SndListResource**)sndHandle, true);
        }
    }

    theString = "";

    AW = GetNewDialog (132,nil,(WindowPtr)(-1));
    TransitionWindow (GetDialogWindow(AW), kWindowZoomTransitionEffect,kWindowShowTransitionAction,NULL);
    SetThemeWindowBackground (GetDialogWindow(AW),kThemeBrushDialogBackgroundActive,true);
    ShowWindow(GetDialogWindow(AW));
    GetDialogItem(AW,2,&iType,&iHandle,&textBox);

    textScrollTimerAction tsat;

    tsat.textBoxInset = textBox;
    tsat.textBoxSlide = textBox;
    tsat.textBoxSlide.left+=2;
    tsat.textBoxSlide.right-=2;
    tsat.textBoxSlide.top+=2;
    tsat.textBoxSlide.bottom-=2;
    InsetRect (&tsat.textBoxInset,-1,-1);
    tsat.textBoxInset.top+=4;


    tsat.AW         = AW;
    tsat.fadeout    = false;
    tsat.morePages  = false;
    tsat.decrement  = 100;
    tsat.startingString
        = 0;
    tsat.textTime   = 0;
    tsat.color      =  0xffff;
    tsat.shutter    =  0;

    tsat.theText    = & theText;
    tsat.callCount  = 0;

    borderRegion    = NewRgn ();


    GetPort(&sPort);

#ifdef OPAQUE_TOOLBOX_STRUCTS
    SetPort(GetDialogPort(AW));
#else
    SetPort(AW);
#endif
    DrawDialog (AW);
    TextFont (20);
    TextSize (10);

    if (splash) {
        theMessage = "The Great Conrholio Lives!";
    } else {
        theMessage = GetVersionString();
    }

    StringToStr255 (theMessage, buffer);
    counter = TextWidth ((Ptr)buffer,1,buffer[0]);
#ifdef OPAQUE_TOOLBOX_STRUCTS
    Rect    awr;
    GetWindowBounds (GetDialogWindow(AW),kWindowGlobalPortRgn,&awr);
    OffsetRect (&awr,-awr.left,-awr.top);
    MoveTo ((awr.right-awr.left-counter)/2,textBox.top-10);
#else
    MoveTo ((AW->portRect.right-AW->portRect.left-counter)/2,textBox.top-10);
#endif

#ifdef TARGET_API_MAC_CARBON
    SetThemeWindowBackground (GetDialogWindow(AW),kThemeBrushWhite,false);
#else
    SetThemeWindowBackground (AW,kThemeBrushWhite,false);
#endif

    EraseRect       (&textBox);
    DrawEmbossedBox (textBox);
    DrawString      (buffer);

    EventLoopRef        mainLoop;
    EventLoopTimerUPP   timerUPP,
                        cTimerUPP;

    EventLoopTimerRef   theTimer,
                        cTimer;

    mainLoop = GetMainEventLoop();
    timerUPP = NewEventLoopTimerUPP(DialogTimerAnimateAction);
    InstallEventLoopTimer (mainLoop,0,kEventDurationMillisecond*50,timerUPP,(void*)&tsat,&theTimer);

    if (splash) {
        cTSD.AW         = AW;
        cTSD.corn1   = GetPicture (whichOne?160:128);
        cTSD.corn2   = GetPicture (whichOne?161:129);
        cTSD.cornRect = (*cTSD.corn1)->picFrame;
        cTSD.cornRect.left = -cTSD.cornRect.right;
        cTSD.cornRect.right = 0;
        cTSD.callCount = 0;
        GetDialogItem(AW,1,&iType,&iHandle,&cTSD.bigRect);
        cTimerUPP = NewEventLoopTimerUPP(DialogTimerAnimateCornholio);
        InstallEventLoopTimer (mainLoop,0,kEventDurationMillisecond*60,cTimerUPP,(void*)&cTSD,&cTimer);
    }


    while (1) {
#ifdef TARGET_API_MAC_CARBON
        if (WaitNextEvent(everyEvent, &theEvent, 0x7FFFFFFF, (RgnHandle)nil))
#else
        if (WaitNextEvent(everyEvent, &theEvent, 0x7FFFFFFF, (RgnHandle)nil))
#endif
        {
            FindWindow(theEvent.where,&whatWindow);
            if (theEvent.what == mouseDown)
#ifdef OPAQUE_TOOLBOX_STRUCTS
                if (splash||(whatWindow == GetDialogWindow(AW)))
#else
                if (splash||(whatWindow == AW))
#endif
                    break;
            if (theEvent.what == keyDown) {
                break;
            }

        }
    }

    if (splash) {
        if (sndChn) {
            SndDisposeChannel (sndChn, true);
        }

        if (sndHandle) {
            HUnlock (sndHandle);
            ReleaseResource (sndHandle);
        }

    }

    if (splash) {
        ReleaseResource((Handle)cTSD.corn1);
        ReleaseResource((Handle)cTSD.corn2);
        RemoveEventLoopTimer     (cTimer);
        DisposeEventLoopTimerUPP (cTimerUPP);
    }

    RemoveEventLoopTimer     (theTimer);
    DisposeEventLoopTimerUPP (timerUPP);
    TransitionWindow         ((WindowPtr)AW, kWindowZoomTransitionEffect,kWindowHideTransitionAction,NULL);
    DisposeDialog(AW);
    DisposeRgn (borderRegion);

    SetPort (sPort);
}

//____________________________________________________________________________________________

bool    handleGUI (bool useGet)
{
    static              _HYPlatformWindow*  lastWindow = nil;
    short                whatW,
                         f;

    WindowPtr            whatWindow;
    bool                 gotEvent = false,
                         loopMore = false;

    EventRecord         guiEvent;
    long                menuChoice,
                        mSel;
    do {
        if (GlobalGUIEventQueue.lLength) {
            HandleGlobalQueueEvent ();
        }
        mSel = -1;

        if (useGet && ! isSuspended) {
            gotEvent = GetNextEvent (everyEvent, &guiEvent);
            loopMore = gotEvent;
        } else {
            gotEvent = WaitNextEvent (everyEvent, &guiEvent, 0x7FFFFFFF, mouseRgn);
        }

        if (highLevelQuit) {
            return false;
        }

        if (!gotEvent) {
            continue;
        }

        if (guiEvent.what==osEvt) {
            if (((guiEvent.message&osEvtMessageMask)>>24)==suspendResumeMessage) {
                Cursor arrow;
                SetCursor(GetQDGlobalsArrow(&arrow));
                continue;
            } else if (trackMouseMovement) {
                f = windowPtrs.Find ((long)FrontWindow());
                if (f>=0) {
                    _HYPlatformWindow* clickedWindow = (_HYPlatformWindow*)windowObjects (f);
                    clickedWindow->_ProcessOSEvent ((Ptr)&guiEvent);
                }
                continue;
            }
        }

        if(guiEvent.what == kHighLevelEvent) {
            if ((OSType)guiEvent.message==typeAppleEvent)
                if ((*((OSType*)(&guiEvent.where))==kAEOpenDocuments)
                        || (*((OSType*)(&guiEvent.where))==kAEQuitApplication)) {
                    AEProcessAppleEvent(&guiEvent);
                    if (highLevelQuit) {
                        return false;
                    }
                    continue;
                }
        }

        if ((guiEvent.what == keyDown)||(guiEvent.what == autoKey)) {
            whatWindow = FrontWindow();
        } else {
            whatW = FindWindow(guiEvent.where,&whatWindow);
        }

        if (guiEvent.what == updateEvt || guiEvent.what == activateEvt) {
            f = windowPtrs.Find ((long)guiEvent.message);
        } else {
            f = windowPtrs.Find ((long)whatWindow);
        }

        if ((guiEvent.what != updateEvt)&&(guiEvent.what != activateEvt)&&(guiEvent.what != nullEvent)&&windowPtrs.lLength) {
            WindowPtr fw = FrontWindow ();
            mSel = windowPtrs.Find((long)fw);
            if ((mSel!=f)&&(mSel>=0)) {
                _HYWindow* dW = (_HYWindow*) windowObjectRefs (mSel);
                if (dW->flags & HY_WINDOW_DLOG) {
                    continue;
                }
            }
            mSel = -1;
        }

        if (f>=0) {
            _HYPlatformWindow* clickedWindow = (_HYPlatformWindow*)windowObjects (f);

            if (clickedWindow!=lastWindow) {
#ifdef TARGET_API_MAC_CARBON
                Cursor arrow;
                SetCursor(GetQDGlobalsArrow(&arrow));
#else
                SetCursor (&qd.arrow);
#endif
                lastWindow = clickedWindow;
            }
            if(clickedWindow->_ProcessOSEvent ((Ptr)&guiEvent)) {
                continue;
            }
        }


        if (guiEvent.what == keyDown) {
            if (guiEvent.modifiers&cmdKey) {
                if (guiEvent.modifiers&optionKey) {
                    if ((guiEvent.message&keyCodeMask)==0x1F00) {
                        mSel = (long)150*0x10000+2;
                    } else if ((guiEvent.message&keyCodeMask)==0x2D00) {
                        mSel = (long)151*0x10000+2;
                    }
                } else {
                    if (guiEvent.modifiers&shiftKey) {
                        if (((guiEvent.message&charCodeMask)=='O')||((guiEvent.message&charCodeMask)=='o')) {
                            mSel = (long)150*0x10000+3;
                        }
                    } else {
                        mSel = MenuKey (guiEvent.message&charCodeMask);
                    }
                }
            }
        }


        if (guiEvent.what == mouseDown) {
            switch (whatW) { //what window contained the last apple event
            case inDesk:
            case inSysWindow: {
#ifndef TARGET_API_MAC_CARBON
                SystemClick(&guiEvent, whatWindow);
#endif
                continue;
            }
            case inMenuBar: {
                mSel = MenuSelect(guiEvent.where);
            }
            }
        }

        if (mSel>0) {
            HiliteMenu(0);
            f = windowPtrs.Find ((long)FrontWindow());
            if (f>=0) {
                _HYPlatformWindow* clickedWindow = (_HYPlatformWindow*)windowObjects (f);
                if(clickedWindow->_ProcessMenuSelection (mSel)) {
                    continue;
                }
            }

            menuChoice = mSel&0x0000ffff;
            switch (mSel/0xffff) {
            case 128: { // apple menu
                if ((mSel&0x0000ffff) == 1)
                    // About Window
                {
                    DisplayAbout(guiEvent.modifiers&optionKey);
                    continue;
                }

#ifndef TARGET_API_MAC_CARBON
                GetMenuItemText(GetMenu(128),mSel&0x0000ffff,daName);
                OpenDeskAcc(daName);
#endif
                continue;
            }

            case 129: { // file menu
                if (menuChoice==4) { // save
                    SaveConsole ();
                    continue;
                }

                if (menuChoice==10) { // quit
                    terminateExecution = true;
                    return false;
                }
                continue;
            }

            case 130: {
                switch (menuChoice) {
                case 11:
                    HandlePreferences (globalPreferencesList,"HYPHY Preferences");
                    break;
                }

                continue;
            }

            case 131: { // action menu
                if (menuChoice==2) { // suspend
                    if (!isSuspended) {
                        HandleSuspend();
                        isSuspended = true;
                    } else {
                        HandleResume();
                        isSuspended = false;
                    }
                } else if (menuChoice == 1) { // cancel execution
                    HandleCancel();
                } else if (menuChoice == 4) { // show messages
                    ShowMessagesLog();
                } else if (menuChoice == 6) { // display template files
                    RunStandardAnalyses();
                } else if (menuChoice == 8) { // rerun last analysus
                    OpenBatchFile(false);
                    ExecuteBatchFile();
                } else if (menuChoice == 10) { // calculate expression
                    if (calculatorMode) {
                        _HYTextBox         *ib = (_HYTextBox*)hyphyConsoleWindow->GetObject(1);
                        ib->SetText ("exit");
                        hyphyConsoleWindow->ProcessEvent (generateTextEditChangeEvent(ib->GetID(),2));
                        calculatorMode         = false;
                        //ib->SetText (empty);
                    } else {
                        calculatorMode = true;
                        while(calculatorMode&&ExpressionCalculator()) {}
                        calculatorMode = false;
                        if (terminateExecution) {
                            return false;
                        }
                    }
                } else if (menuChoice == 11) { // calculate expression
                    ExecuteSelection();
                }

                continue;
            }

            case 132: { // window menu
                if (menuChoice == 1) { // console window
                    hyphyConsoleWindow->BringToFront();
                } else if (menuChoice == 2) { // object inspector
                    ShowObjectInspector ();
                } else if (menuChoice == 3) { // cycle thru windows
                    if (windowPtrs.lLength>1) {
                        menuChoice = windowPtrs.Find ((long)FrontWindow());
                        if (menuChoice == windowPtrs.lLength-1) {
                            SelectWindow ((WindowPtr)(windowPtrs.lData[0]));
                        } else {
                            SelectWindow ((WindowPtr)(windowPtrs.lData[menuChoice+1]));
                        }
                    }
                } else if (menuChoice == 4) { // hide windows
                    if ((WindowPtr)hyphyConsoleWindow->GetOSWindowData() != FrontWindow()) {
                        HideWindow (FrontWindow());
                    }
                } else {
                    long f = FindWindowByName (objectInspectorTitle);
                    if((f>=0)&&(menuChoice-6>=f)) {
                        menuChoice++;
                    }
                    ShowWindow   ((WindowPtr)(windowPtrs.lData[menuChoice-6]));
                    SelectWindow ((WindowPtr)(windowPtrs.lData[menuChoice-6]));
                }
                continue;
            }

            case 150: { // open submenu
                if (menuChoice==1) { // open batch file
                    if (OpenBatchFile()) {
                        ExecuteBatchFile ();
                    }
                } else if (menuChoice==2) { // open tree file
                    OpenDataFile();
                } else if (menuChoice==3) { // open tree file
                    OpenTreeFile();
                } else if (menuChoice==7) { // open text file
                    OpenTextFile();
                } else if (menuChoice==8) { // open table
                    OpenTable ();
                } else if (menuChoice==9) { // open table
                    OpenDatabaseFile (nil);
                } else if ((menuChoice==11)||(menuChoice==12)) { // convert movie file
                    ConvertMovieFile (menuChoice-10);
                }
                continue;
            }

            case 151: { // new submenu
                if (menuChoice==1) { // open batch file
                    NewTreeWindow();
                } else if (menuChoice==2) { // open batch file
                    NewModel(nil);
                } else if (menuChoice==3) {
                    NewChartWindow();
                } else if (menuChoice==4) {
                    NewGeneticCodeTable(0);
                } else if (menuChoice==5) {
                    NewDatabaseFile(nil);
                }
                continue;
            }

            case 200: { // recent files menu
                if (!argFileName) {
                    argFileName = new _String (*(_String*)recentPaths(menuChoice-1));
                } else {
                    *argFileName = *(_String*)recentPaths(menuChoice-1);
                }

                volumeName          =   argFileName->Cut(0,argFileName->Find(':'));
                OpenBatchFile       (false);
                ExecuteBatchFile    ();
                continue;

            }

            case 201: { // results menu
                ExecuteAPostProcessor (*(_String*)(*(_List*)availablePostProcessors(menuChoice-1))(1));
                continue;
            }
            }

        }

    } while (isSuspended||loopMore);


    return true;
}

//___________________________________________________________________
void InitToolbox()
{
    FlushEvents(everyEvent,0);
    InitCursor();
    long result;
    Gestalt( gestaltMenuMgrAttr, &result);
    if (result & gestaltMenuMgrAquaLayoutMask) {
        aquaInterfaceOn = true;
    }
}



//_________________________________________________________________________
void    HandleSuspend (void)
{
    SetStatusLine (empty,empty,empty,-1,HY_SL_SUSPEND);
    SetMenuItemText (GetMenuHandle(131),2,"\pResume");
    DisableMenuItem(GetMenuHandle(131),1);
}

//_________________________________________________________________________
void    HandleResume (void)
{
    SetStatusLine (empty,empty,empty,-1,HY_SL_RESUME);
    SetMenuItemText (GetMenuHandle(131),2,"\pSuspend");

    time_t      tt;
    timerStart += time(&tt)-lastTimer;

    EnableMenuItem(GetMenuHandle(131),1);
}

//_________________________________________________________________________
void    HandleCancel (void)
{
    SetStatusLine (empty,"Canceling",empty,-1,HY_SL_TASK);
    terminateExecution = true;
}

//_________________________________________________________________________
void    SetupMenusAndSuch (void)
{
    bool good = true;
    if (AEInstallEventHandler ((OSType)typeAppleEvent,(OSType)kAEOpenDocuments,NewAEEventHandlerUPP((AEEventHandlerProcPtr)GUIOpenHandler),0,false)) {
        good = false;
    }
    if (AEInstallEventHandler ((OSType)typeAppleEvent,(OSType)kAEQuitApplication,NewAEEventHandlerUPP((AEEventHandlerProcPtr)HandleQuit),0,false)) {
        good = false;
    }

    if (!good) {
        _String        errMsg ("Apple Event Initialization Error. Quitting....");
        ProblemReport (errMsg);
        exit(1);
    }

    AppendResMenu(GetMenu(128),'DRVR');

    if (aquaInterfaceOn) {
        DeleteMenuItem (GetMenuHandle (129), 9);
        DeleteMenuItem (GetMenuHandle (129), 9);
        DeleteMenuItem (GetMenuHandle (130), 10);
        DeleteMenuItem (GetMenuHandle (130), 10);
        EnableMenuCommand(NULL, kHICommandPreferences);
        EventTypeSpec   eventTypes[1];
        EventHandlerUPP handlerUPP;
        eventTypes[0].eventClass = kEventClassCommand;
        eventTypes[0].eventKind  = kEventCommandProcess;
        handlerUPP = NewEventHandlerUPP(PrefAction);
        InstallApplicationEventHandler (handlerUPP,1, eventTypes,NULL, NULL);
    }


    Str255  buffer2;
    long counter = 1;


    EnableMenuItem  (GetMenuHandle(130),0);
    DisableMenuItem (GetMenuHandle(130),1);
    DisableMenuItem (GetMenuHandle(130),3);
    DisableMenuItem (GetMenuHandle(130),4);
    DisableMenuItem (GetMenuHandle(130),5);
    DisableMenuItem (GetMenuHandle(130),6);
    //DisableMenuItem (GetMenuHandle(130),8);

    DisableMenuItem (GetMenuHandle(131),1);
    DisableMenuItem (GetMenuHandle(131),2);
    if (!hasTemplates) {
        DisableMenuItem (GetMenuHandle(131),6);
    }

    DisableMenuItem (GetMenuHandle(131),7);
    DisableMenuItem (GetMenuHandle(131),8);

    // set up the 2nd menu
    InsertMenu(GetMenu (150), hierMenu);
    InsertMenu(GetMenu (151), hierMenu);
    SetItemCmd (GetMenuHandle(129),1,hMenuCmd);
    SetItemMark(GetMenuHandle(129),1,151);
    SetItemCmd (GetMenuHandle(129),2,hMenuCmd);
    SetItemMark(GetMenuHandle(129),2,150);
    //SetMenuItemText (GetMenuHandle(150),5,"\pRun Recent");
    MenuHandle openMenu = GetMenuHandle(150);
    SetItemCmd (openMenu,5,hMenuCmd);
    SetItemMark(openMenu,5,200);
    InsertMenu(recentFilesMenu,-1);
    if (!recentPaths.lLength) {
        DisableMenuItem(GetMenuHandle(150),5);
    }

    SetMenuItemModifiers (openMenu, 2, kMenuOptionModifier);
    SetMenuItemModifiers (openMenu, 3, kMenuShiftModifier);
    SetMenuItemModifiers (GetMenuHandle(151), 2, kMenuOptionModifier);
    //SetMenuItemModifiers (openMenu, 7, kMenuShiftModifier|kMenuOptionModifier);

    if (availablePostProcessors.countitems()) {
        MenuHandle resultsMenu = NewMenu (201,"\pResults");
        for (counter=0; counter<availablePostProcessors.countitems(); counter++) {
            _String * postItem = (_String*)(*(_List*)availablePostProcessors(counter))(0);
            if (postItem->Equal(&menuSeparator)) {
                AppendMenu (resultsMenu,"\p(-");
            } else {
                StringToStr255 (*postItem,buffer2);
                AppendMenu (resultsMenu,buffer2);
            }
        }
        SetItemCmd (GetMenuHandle(131),7,hMenuCmd);
        SetItemMark(GetMenuHandle(131),7,201);
        SetMenuItemText (GetMenuHandle(131),7,"\pResults");
        InsertMenu(resultsMenu,hierMenu);
        DisableMenuItem(GetMenuHandle(131),7);
    }

    InvalMenuBar();
    GetIndPattern (&penHatchPattern,0,8);
    GetIndPattern (&vertPenHatchPattern,0,27);
}

//_________________________________________________________________________
void    CleanupMenusAndSuch (void)
{
    DisposeHandle(menuBar);
    DisposeMenu  (recentFilesMenu);
}

//________________________________________________________

extern  ModalFilterUPP myFilterProc;

//________________________________________________________
bool        HandleFirstDialog (long& choice)
{
    bool           checkFlag = false;
    short      iType;
    Handle     iHandle;
    Rect       iRect,
               textBox;
    Str255     buffer;
    GrafPtr    sPort;

    DialogPtr  theDlg = GetNewDialog (143,nil,(WindowPtr)-1);
    if (!theDlg) {
        return false;
    }

    GetPort(&sPort);
#ifdef OPAQUE_TOOLBOX_STRUCTS
    SetPort(GetDialogPort(theDlg));
#else
    SetPort(theDlg);
#endif

    SetDialogDefaultItem (theDlg,1);
    GetDialogItem (theDlg,3,&iType,&iHandle,&textBox);
    GetDialogItem (theDlg,2,&iType,&iHandle,&iRect);
    iType = 0;
#ifdef OPAQUE_TOOLBOX_STRUCTS
    SelectWindow (GetDialogWindow(theDlg));
#else
    SelectWindow (theDlg);
#endif
    SetControlValue ((ControlHandle)iHandle,checkFlag);

    DrawDialog (theDlg);
    TextFont (20);
    TextSize (10);
    _String theMessage = GetVersionString();
    StringToStr255 (theMessage, buffer);
    long            counter = TextWidth ((Ptr)buffer,1,buffer[0]);

#ifdef OPAQUE_TOOLBOX_STRUCTS
    Rect    awr;
    GetWindowBounds (GetDialogWindow(theDlg),kWindowGlobalPortRgn,&awr);
    OffsetRect (&awr,-awr.left,-awr.top);
    MoveTo ((awr.right-awr.left-counter)/2,textBox.top-10);
#else
    MoveTo ((theDlg->portRect.right-theDlg->portRect.left-counter)/2,textBox.top-10);
#endif

    DrawString (buffer);

    while (1) {
        ModalDialog (myFilterProc,&iType);
        if (iType == 2) {
            checkFlag = !checkFlag;
            SetControlValue ((ControlHandle)iHandle,checkFlag);
        } else {
            break;
        }
    }
    if (iType > 1) {
        choice = iType-2;
    } else {
        choice = 0;
    }


    SetPort (sPort);
    DisposeDialog (theDlg);
    return !checkFlag;
}


//_________________________________________________________________________

int main (void)
{
#ifndef __HYPHYXCODE__
    _fcreator = 'MuSe';
    _ftype    = 'TEXT';
#endif

    RegisterAppearanceClient();

    InitToolbox();
    whiteFill       = GetPixPat(128);
    blueFill        = GetPixPat(131);
    statusBarFill   = GetPixPat(134);

    redButtonIcon       = GetCIcon (4000);
    yellowButtonIcon    = GetCIcon (4001);
    greenButtonIcon     = GetCIcon (4002);
    orangeButtonIcon    = GetCIcon (4003);
    pullDownArrowsIcon  = GetCIcon (4010);
    tablePDMenuIcon     = GetCIcon (4020);

    CursHandle          hCursor = GetCursor(132);
    if (hCursor) {
        HLock ((Handle)hCursor);
        hSizeCursor = ** hCursor;
        DisposeHandle ((Handle)hCursor);
    }

    hCursor = GetCursor(iBeamCursor);
    if (hCursor) {
        //HLock ((Handle)hCursor);
        editStateCursor = ** hCursor;
        //DisposeHandle ((Handle)hCursor);
    }

    pickUpCursor    = GetCCursor (128);
    dropOffCursor   = GetCCursor (129);


    char buffer[4096];
#ifdef __HYPHYXCODE__
    {
        CFBundleRef myAppsBundle = CFBundleGetMainBundle();
        CFURLRef myBundleURL = CFBundleCopyBundleURL(myAppsBundle);
        FSRef myBundleRef;
        CFURLGetFSRef(myBundleURL, &myBundleRef);
        CFRelease(myBundleURL);
        FSSpec   baseDirSpec;
        FSGetCatalogInfo(&myBundleRef, kFSCatInfoNone,NULL, NULL, &baseDirSpec, NULL);
        GetFullPathName (baseDirSpec, baseDirectory);
        baseDirectory.Trim(0,baseDirectory.FindBackwards (':',0,-1));
    }
#else
    getcwd (buffer,4095);
    baseDirectory = buffer;
#endif
    pathNames && & baseDirectory;

    libDirectory = baseDirectory & APPNAME ".app" & GetPlatformDirectoryChar() & "Contents" & GetPlatformDirectoryChar() & "Resources" & GetPlatformDirectoryChar() & "HBL" & GetPlatformDirectoryChar();
    pathNames && & libDirectory;

    if (!(menuBar = GetNewMBar(128))) {
        _String        errMsg("Error reading menu resources. Quitting....");
        ProblemReport (errMsg);
        exit(1);
    }
    // fprintf(stderr, "%s\n", libDirectory.sData);

    SetMenuBar(menuBar);
    GlobalStartup();
#ifdef __MP__
    if (!(systemCPUCount = MPProcessors())) {
        systemCPUCount = 1;
    }
#endif



    OSErr           status;
    TXNInitOptions options;
    TXNMacOSPreferredFontDescription  defaults;
    ATSUFontID      theFontID;

    defaults.fontID = ATSUFindFontFromName ("Times Roman",
                                            strlen("Times Roman"),
                                            kFontFullName, kFontNoPlatformCode,
                                            kFontNoScriptCode, kFontNoLanguageCode,
                                            &theFontID);

    defaults.pointSize = kTXNDefaultFontSize;
    defaults.fontStyle = kTXNDefaultFontStyle;
    defaults.encoding  = CreateTextEncoding (kTextEncodingMacRoman,
                         kTextEncodingDefaultVariant,
                         kTextEncodingDefaultFormat);

    if (aquaInterfaceOn) {
        options  =  0;
    } else {
        options  =  kTXNAlwaysUseQuickDrawTextMask;
    }

    status = TXNInitTextension (&defaults, 1, options);
    if (status != noErr) {
        FlagError ("Failed to initialize MLTE. Check the version of your CarbonLib (get the latest)");
        return -1;
    }

    hyphyConsoleWindow = new _HYConsoleWindow ("HYPHY Console");
    ReadPreferences();
#ifdef __MP__
    if (systemCPUCount == 1) {
        BufferToConsole ("One processor detected.\n");
    } else {
        snprintf (buffer, sizeof(buffer),"%ld processors detected.\n", systemCPUCount);
        BufferToConsole (buffer);
    }
#endif
    

    MoveConsoleWindow       (consolePositionRectangle);
    SetPreferences          ();
    ReadInTemplateFiles     ();
    SetStatusLine           ("None","Idle","00:00:00");
    hasTemplates = availableTemplateFiles.lLength;
    SetupMenusAndSuch();
    ReadModelTemplates ();
    ReadGeneticCodes ();
    StringToConsole (hyphyCiteString);

    hyphyConsoleWindow->Activate();

    mouseRgn = NewRgn();
    Rect     dummy = {0,0,1,1};
    RectRgn  (mouseRgn,&dummy);
    trackMouseMovement = true;

    long            choice = 0;

    if (showDialogAtStartup) {
        showDialogAtStartup = HandleFirstDialog (choice);

        if (!showDialogAtStartup) {
            _String wStr ("You can later change startup dialog settings in 'Preferences'");
            ProblemReport (wStr);
        }

        if (choice == 1) {
            RunStandardAnalyses();
        } else {
            _String dL (baseDirectory);
            if (choice == 2) {
                dL = dL & "data:";
                OpenDataFile (&dL);
            } else if (choice == 3) {
                dL = dL & "Saves:";
                if (OpenBatchFile (true,&dL)) {
                    ExecuteBatchFile ();
                }
            }
        }

        SetShowDialogAtStartup (showDialogAtStartup);
    }


#ifdef TARGET_API_MAC_CARBON
    EventLoopRef mainLoop;
    EventLoopTimerUPP timerUPP;
    EventLoopTimerRef theTimer;
    mainLoop = GetMainEventLoop();
    timerUPP = NewEventLoopTimerUPP(TimerAction);
    InstallEventLoopTimer (mainLoop,0,kEventDurationSecond,timerUPP,NULL,&theTimer);
#endif


    bool  canQuit;
    do {
        canQuit = true;
        while (handleGUI()) ; // main loop
        highLevelQuit = false;
        for (long i=windowObjects.countitems()-1; i>=1; i--) {
            _HYPlatformWindow* thisWindow = (_HYPlatformWindow*)windowObjects(i);
            _HYWindow*         pWindow = (_HYWindow*)thisWindow;
            postWindowCloseEvent (pWindow->GetID());
            handleGUI (false);
            //if (!pWindow->Close(nil))
            //break;
        }
        handleGUI (true);
        if (windowObjects.lLength>1) {
            canQuit = false;
        }
    } while (!canQuit);
    // close all the windows with prompts

    GlobalShutdown();
    WritePreferences();
    CleanupMenusAndSuch();

    batchLanguageFunctions.Clear();
    PurgeAll(true);

    if (mouseRgn) {
        DisposeRgn  (mouseRgn);
    }

    DisposePixPat(whiteFill);
    DisposePixPat(blueFill);
    DisposePixPat(statusBarFill);

    DisposeCIcon (redButtonIcon);
    DisposeCIcon (yellowButtonIcon);
    DisposeCIcon (greenButtonIcon);
    DisposeCIcon (pullDownArrowsIcon);
    DisposeCIcon (tablePDMenuIcon);


#ifdef TARGET_API_MAC_CARBON
    RemoveEventLoopTimer (theTimer);
    DisposeEventLoopTimerUPP (timerUPP);
#endif
    //ProfilerDump("\pProfile");
    //ProfilerTerm();
    return 0;
}




