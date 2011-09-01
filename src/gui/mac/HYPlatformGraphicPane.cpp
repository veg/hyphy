/*
    A painting canvas with double buffer. MacOS.

    Sergei L. Kosakovsky Pond, May 2000.
*/

#include "HYGraphicPane.h"
#include "QDOffscreen.h"
#include "Quickdraw.h"
//#include "Textedit.h"
#include "Fonts.h"
#include "stdlib.h"
#include "errorfns.h"
#include "scrap.h"
#include "QuickTimeComponents.h"
#include "ImageCompression.h"
#include "HYUtils.h"
#include "HYPulldown.h"
#include "HYPlatformWindow.h"
#include "Appearance.h"
#include "Resources.h"
#include "MacTextEditor.h"

extern      void StringToStr255 (_String&, Str255&);

extern      Pattern penHatchPattern,vertPenHatchPattern;

void        findGraphicsExporterComponents (_List&, _SimpleList&);

_List       graphicsFormats;

_SimpleList qtGrexComponents;

_String     savePicPrompt               ("Save Picture As:"),
            savePicAs                  ("File Format:");

//__________________________________________________________________

Rect    HYRect2Rect (_HYRect& hr)
{
    Rect r;
    r.left = hr.left;
    r.right = hr.right;
    r.top = hr.top;
    r.bottom = hr.bottom;
    return r;
}

//__________________________________________________________________

void findGraphicsExporterComponents(_List& compList, _SimpleList& compIndex)
{
    ComponentDescription cd, cd2;
    Component c = 0;

    cd.componentType            = GraphicsExporterComponentType;
    cd.componentSubType         = 0;
    cd.componentManufacturer    = 0;
    cd.componentFlags           = 0;
    cd.componentFlagsMask       = graphicsExporterIsBaseExporter;

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

//__________________________________________________________________

void InitializeQTExporters(void)
{
    if (!graphicsFormats.lLength) {
        qtGrexComponents.Clear();
        findGraphicsExporterComponents (graphicsFormats, qtGrexComponents);
    }
}

//__________________________________________________________________

_HYPlatformGraphicPane::_HYPlatformGraphicPane(int h, int w, int d)
{
    Rect  bRect;
    bRect.left = bRect.top = 0;
    bRect.right = w;
    bRect.bottom = h;

    short errCode;

    if (d>1) {
        errCode = NewGWorld (&thePane,d,&bRect,0,GetMainDevice(),noNewDevice);

        if (errCode == -108) { // no memory
            errCode = NewGWorld (&thePane,d,&bRect,0,GetMainDevice(),noNewDevice|useTempMem);
        }
    } else {
        errCode = NewGWorld (&thePane,d,&bRect,0,nil,0);

        if (errCode == -108) { // no memory
            errCode = NewGWorld (&thePane,d,&bRect,0,nil,useTempMem);
        }
    }

    fillColor = NewPixPat();
    //backColor = NewPixPat();
    if (errCode||(!fillColor)) {
        _String errMsg ("MacOS Error ");
        errMsg = errMsg & (long)errCode &" while trying to allocate memory for GraphicPane";
        FlagError (errMsg);
    }
    savedPort = nil;
}

//__________________________________________________________________

_HYPlatformGraphicPane::~_HYPlatformGraphicPane(void)
{
    DisposeGWorld (thePane);
    DisposePixPat (fillColor);
    //DisposePixPat (backColor);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_SetPaneSize  (int h,int w, int d)
{
    DisposeGWorld (thePane);
    Rect  bRect;
    bRect.left = bRect.top = 0;
    bRect.right = w;
    bRect.bottom = h;
    short errCode;
    if (d>1) {
        errCode = NewGWorld (&thePane,d,&bRect,0,GetMainDevice(),noNewDevice);

        if (errCode == -108) { // no memory
            errCode = NewGWorld (&thePane,d,&bRect,0,GetMainDevice(),noNewDevice|useTempMem);
        }
    } else {
        errCode = NewGWorld (&thePane,d,&bRect,0,nil,0);

        if (errCode == -108) { // no memory
            errCode = NewGWorld (&thePane,d,&bRect,0,nil,useTempMem);
        }
    }

    if (errCode) {
        _String errMsg ("MacOS Error ");
        errMsg = errMsg & (long)errCode &" while trying to allocate memory for GraphicPane";
        FlagError (errMsg);
    }
}


//__________________________________________________________________

void _HYPlatformGraphicPane::_DrawLine (_HYRect lineDesc)
{
    PenSize (lineDesc.width, lineDesc.width);
    MoveTo  (lineDesc.left, lineDesc.top);
    LineTo  (lineDesc.right,lineDesc.bottom);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_DrawHatchedLine (_HYRect lineDesc)
{
    PenState savePen;
    GetPenState (&savePen);
    if (abs(lineDesc.left-lineDesc.right)>5) {
        PenPat (&penHatchPattern);
    } else {
        PenPat (&vertPenHatchPattern);
    }
    PenSize (lineDesc.width, lineDesc.width);
    MoveTo  (lineDesc.left, lineDesc.top);
    LineTo  (lineDesc.right,lineDesc.bottom);
    SetPenState (&savePen);
}


//__________________________________________________________________
void _HYPlatformGraphicPane::_DisplayText    (_String theText,int t, int l, bool dir)
{
    MoveTo (l,t);
    if (!dir) {
        _HYGraphicPane* theParent = (_HYGraphicPane*)this;
        long            fontSize = theParent->font.size+1;

        t+=fontSize;
        for (long k = 0; k < theText.sLength; k++, t += fontSize) {
            DrawChar (theText.sData[k]);
            MoveTo (l, t);
        }
    } else {
        DrawText (theText.sData,0,theText.sLength);
    }
}

//__________________________________________________________________
void _HYPlatformGraphicPane::_DisplayText    (_String& theText,_HYRect& r, char align)
{
    _EraseRect (r);
    Rect    tb = HYRect2Rect (r);
    CFStringRef cTextRef = CFStringCreateWithCString (nil, theText.sData, kCFStringEncodingMacRoman);
    TXNTextBoxOptionsData opts;
    opts.optionTags = kTXNSetFlushnessMask;
    opts.flushness  = (align==HY_ALIGN_LEFT)?kATSUStartAlignment:(align==HY_ALIGN_RIGHT?kATSUEndAlignment:kATSUCenterAlignment);
    TXNDrawCFStringTextBox (cTextRef, &tb, nil, &opts);
    CFRelease (cTextRef);
    //short saveTextMode = thePane->txMode;
    //TextMode (srcOr);
    //TETextBox(theText.sData, theText.sLength, &tb, (align==HY_ALIGN_LEFT)?0:(align==HY_ALIGN_RIGHT?-1:1));
    //TextMode (saveTextMode);
}

//__________________________________________________________________
void _HYPlatformGraphicPane::_DisplayChar  (char c,int t, int l)
{
    MoveTo (l,t);
    DrawChar (c);
}

//__________________________________________________________________
void _HYPlatformGraphicPane::_SlidePane  (int dv, int dh)
{
    _HYGraphicPane* theParent = (_HYGraphicPane*)this;
    Rect r = {0,0,theParent->h,theParent->w};
    RgnHandle dummy = NewRgn();
    checkPointer (dummy);
    ScrollRect(&r,dh,dv,dummy);
    DisposeRgn (dummy);
}

//__________________________________________________________________
void _HYPlatformGraphicPane::_SlideRect (_HYRect& rct, int dv, int dh)
{
    Rect r;
    r.left = rct.left;
    r.top = rct.top;
    r.bottom = rct.bottom;
    r.right = rct.right;
    RgnHandle dummy = NewRgn();
    checkPointer (dummy);
    ScrollRect(&r,dh,dv,dummy);
    DisposeRgn (dummy);
}

//__________________________________________________________________
void _HYPlatformGraphicPane::_InvertRect (_HYRect& rct)
{
    Rect r = HYRect2Rect (rct);
    InvertRect (&r);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_DrawRect (_HYRect rectDesc)
{
    PenSize (rectDesc.width, rectDesc.width);
    Rect    r= HYRect2Rect (rectDesc);
    FrameRect (&r);

}
//__________________________________________________________________

void _HYPlatformGraphicPane::_FillRect (_HYRect rectDesc)
{
    Rect    r= HYRect2Rect (rectDesc);
    FillCRect (&r,fillColor);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_EraseRect (_HYRect rectDesc)
{
    Rect    r= HYRect2Rect (rectDesc);
    EraseRect (&r);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_DrawOval (_HYRect rectDesc)
{
    PenSize (rectDesc.width, rectDesc.width);
    Rect    r= HYRect2Rect (rectDesc);
    FrameOval (&r);

}
//__________________________________________________________________

void _HYPlatformGraphicPane::_FillOval (_HYRect rectDesc)
{
    Rect    r= HYRect2Rect (rectDesc);
    FillCOval (&r,fillColor);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_EraseOval (_HYRect rectDesc)
{
    Rect    r= HYRect2Rect (rectDesc);
    EraseOval (&r);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_DrawArc (_HYRect rectDesc, int s, int f)
{
    PenSize (rectDesc.width, rectDesc.width);
    Rect    r= HYRect2Rect (rectDesc);
    FrameArc (&r,s,f);

}
//__________________________________________________________________

void _HYPlatformGraphicPane::_FillArc (_HYRect rectDesc, int s, int f)
{
    Rect    r= HYRect2Rect (rectDesc);
    FillCArc (&r,s,f,fillColor);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_EraseArc (_HYRect rectDesc, int s, int f)
{
    Rect    r= HYRect2Rect (rectDesc);
    EraseArc (&r,s,f);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_SetColor  (_HYColor c)
{
    RGBColor      sysColor;
    sysColor.red   = c.R*256;
    sysColor.green = c.G*256;
    sysColor.blue  = c.B*256;
    MakeRGBPat (fillColor,&sysColor);
    RGBForeColor (&sysColor);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_SetBColor  (_HYColor c)
{
    RGBColor         sysColor;
    sysColor.red   = c.R*256;
    sysColor.green = c.G*256;
    sysColor.blue  = c.B*256;

    if (c.R+c.B+(long)c.G==765)
        sysColor = (RGBColor) {
        0xffff,0xffff,0xffff
    };

    RGBBackColor (&sysColor);
    //MakeRGBPat (backColor,&sysColor);
    //BackPat      (backColor);
}


//__________________________________________________________________

void _HYPlatformGraphicPane::_SetFont (_HYFont f)
{
    Str255  fontFace;
    StringToStr255 (f.face,fontFace);
    short fNum;
    GetFNum (fontFace,&fNum);
    TextFont (fNum);
    TextSize (f.size);
    TextFace (f.style);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_DrawPicRes (_HYRect& r, long id)
{
    PicHandle  aPic = GetPicture (id);
    if (aPic) {
        Rect        aRect = HYRect2Rect (r);
        PictInfo    pInfo;
        GetPictInfo (aPic,&pInfo,0,0,0,0);

        if (aRect.right-aRect.left<=0) {
            r.right = aRect.right = aRect.left + pInfo.sourceRect.right - pInfo.sourceRect.left;
        }

        if (aRect.bottom-aRect.top<=0) {
            r.bottom = aRect.bottom = aRect.top + pInfo.sourceRect.bottom - pInfo.sourceRect.top;
        }

        DrawPicture (aPic, &aRect);

        ReleaseResource ((Handle)aPic);
    } else {
        _String errMsg = _String ("No picture resource with ID ") & id;
        ReportWarning (errMsg);
    }
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_SetFontSize (long s)
{
    TextSize (s);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_SetDialogBG (void)
{
    if (aquaInterfaceOn) {
        SetThemeBackground (kThemeBrushDialogBackgroundActive, 32, true);
    } else {
        ((_HYGraphicPane*)this)->SetBColor         (GetDialogBackgroundColor());
    }
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_StartDraw  (void)
{
    _HYGraphicPane * parent = (_HYGraphicPane*)this;
    GetGWorld (&savedPort,&savedDevice);
    ::GetForeColor (&saveFG);
    ::GetBackColor (&saveBG);
    SetGWorld (thePane,nil);
    LockPixels (GetGWorldPixMap(thePane));
    RGBColor c  = {256*parent->bColor.R,256*parent->bColor.G,256*parent->bColor.B};
    if (parent->bColor.R+parent->bColor.B+(long)parent->bColor.G==765)
        c = (RGBColor) {
        0xffff,0xffff,0xffff
    };
    RGBBackColor (&c);
    c = (RGBColor) {
        256*parent->color.R,256*parent->color.G,256*parent->color.B
    };
    RGBForeColor (&c);
    //if (parent->depth>1)
    TextMode     (srcOr);
    //else
    //TextMode   (srcCopy);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_EndDraw    (void)
{
    UnlockPixels (GetGWorldPixMap(thePane));
    SetGWorld    (savedPort,savedDevice);
    RGBForeColor (&saveFG);
    RGBBackColor (&saveBG);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_SetPort    (Ptr nP)
{
    thePane = (GWorldPtr)nP;
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_CopyToClipboard   (void)
{
    _HYGraphicPane* parent = (_HYGraphicPane*)this;
#ifdef TARGET_API_MAC_CARBON
    ClearCurrentScrap();
#else
    ZeroScrap();
#endif
    Rect  bRect;

    bRect.left          = bRect.top = 0;
    bRect.right         = parent->w;
    bRect.bottom        = parent->h;

    PicHandle    pic    = OpenPicture (&bRect);

    GrafPtr      topPort;
    GetPort      (&topPort);

    LockPixels (GetGWorldPixMap(thePane));
#ifdef OPAQUE_TOOLBOX_STRUCTS
    CopyBits (GetPortBitMapForCopyBits(thePane),GetPortBitMapForCopyBits(topPort),
              &bRect,&bRect,srcCopy,(RgnHandle)nil);
#else
    CopyBits ((BitMap*)*GetGWorldPixMap(thePane),
              (BitMap*)&(topPort->portBits),&bRect,&bRect,
              srcCopy,(RgnHandle)nil);
#endif
    UnlockPixels (GetGWorldPixMap(thePane));

    ClosePicture ();
    HLock   ((Handle)pic);

#ifdef TARGET_API_MAC_CARBON
    ScrapRef         theScrapRef;
    GetCurrentScrap(&theScrapRef);
    PutScrapFlavor(theScrapRef, 'PICT', kScrapFlavorMaskNone,GetHandleSize((Handle)pic),*pic);
#else
    PutScrap (GetHandleSize((Handle)pic),'PICT',*pic);
#endif
    KillPicture (pic);
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_SavePicture   (_String prompt)
{
    InitializeQTExporters  ();
    if (graphicsFormats.lLength) {
        _String filePath;
        long    menuChoice = SaveFileWithPopUp (filePath,
                                                savePicPrompt,prompt,savePicAs,graphicsFormats);

        if (menuChoice>=0) {
            ComponentInstance grexc = OpenComponent ((Component)qtGrexComponents(menuChoice));
            GraphicsExportSetInputGWorld (grexc,thePane);
            FSSpec  fs;
            Str255  buff;
            StringToStr255 (filePath,buff);
            FSMakeFSSpec(0,0,buff,&fs);
            GraphicsExportSetOutputFile (grexc,&fs);
            GraphicsExportRequestSettings (grexc,nil,nil);
            unsigned long dummy;
            OSType t,c;
            GraphicsExportGetDefaultFileTypeAndCreator (grexc,&t,&c);
            GraphicsExportSetOutputFileTypeAndCreator (grexc,t,c);
            GraphicsExportDoExport (grexc,&dummy);
            CloseComponent (grexc);
        }
    }
}

//__________________________________________________________________

Ptr _HYPlatformGraphicPane::_DefinePolygon  (_SimpleList& points)
{
    if ((points.lLength>=6)&&(points.lLength%2==0)) {
        PolyHandle   rgn = OpenPoly  ();
        checkPointer (rgn);

        MoveTo       (points.lData[0], points.lData[1]);

        for (long k=2; k<points.lLength; k+=2) {
            LineTo (points.lData[k], points.lData[k+1]);
        }

        LineTo       (points.lData[0], points.lData[1]);
        ClosePoly    ();
        return       (Ptr)rgn;

    }
    return nil;
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_KillPolygon   (Ptr rgn)
{
    if (rgn) {
        KillPoly ((PolyHandle)rgn);
    }
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_DrawPolygon (Ptr rgn, long width)
{
    if (rgn) {
        PenSize (width, width);
        FramePoly ((PolyHandle)rgn);
    }
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_FillPolygon (Ptr rgn)
{
    if (rgn) {
        FillCPoly ((PolyHandle)rgn,fillColor);
    }
}

//__________________________________________________________________

void _HYPlatformGraphicPane::_ErasePolygon (Ptr rgn)
{
    if (rgn) {
        ErasePoly ((PolyHandle)rgn);
    }
}



//EOF