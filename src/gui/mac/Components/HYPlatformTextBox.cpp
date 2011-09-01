/*
    Text input box for Mac OS API

    Sergei L. Kosakovsky Pond, May 2000-December 2002

    Revised in December 2003 to use MLTE
*/

#include "errorfns.h"
#include "HYTextbox.h"
#include "HYUtils.h"
#include "HYEventTypes.h"
#include "HYGraphicPane.h"

#include "ToolUtils.h"
#include "Appearance.h"
#include "string.h"

#include "HYDialogs.h"

//#include "MacTextEditor.h"


//__________________________________________________________________

_HYPlatformTextBox::_HYPlatformTextBox  (void)
{
    backFill  = NewPixPat();
    backTFill = NewPixPat();

    if (!(backFill&&backTFill)) {
        warnError (-108);
    }

#ifdef          _HY_USE_MLTE_
    txn     = nil;
#else
    te      = nil;
#endif

    RGBColor    wht = {0xffff,0xffff,0xffff};
    MakeRGBPat (backFill,&wht);
    MakeRGBPat (backTFill,&wht);

    textBoxRect.left    = textBoxRect.top   = 0;
    textBoxRect.bottom  = textBoxRect.right = 100;

}

//__________________________________________________________________

_HYPlatformTextBox::~_HYPlatformTextBox (void)
{
    if (backFill) {
        DisposePixPat (backFill);
    }

    if (backTFill) {
        DisposePixPat (backTFill);
    }

#ifdef          _HY_USE_MLTE_
    if (txn) {
        TXNDeleteObject (txn);
    }
#else
    if (te) {
        TEDispose(te);
    }
#endif
}

//__________________________________________________________________

void    _HYPlatformTextBox::_SetBackColor (_HYColor& c)
{
    RGBColor newBG;
    newBG.red   = c.R*256;
    newBG.blue  = c.B*256;
    newBG.green = c.G*256;
    MakeRGBPat (backFill,&newBG);
}

//__________________________________________________________________

void    _HYPlatformTextBox::_SetBackTColor (_HYColor& c)
{
    RGBColor newBG;
    newBG.red   = c.R*256;
    newBG.blue  = c.B*256;
    newBG.green = c.G*256;
    MakeRGBPat (backTFill,&newBG);
}

//__________________________________________________________________

long    _HYPlatformTextBox::_GetCharacterCount (void)
{
#ifdef          _HY_USE_MLTE_
    if (txn) {
        return TXNDataSize (txn);
    }
#endif
    return 0;
}


//__________________________________________________________________

void    _HYPlatformTextBox::_SetForeColor (_HYColor& c)
{
#ifdef          _HY_USE_MLTE_
    if (txn) {
        RGBColor newC = {c.R*256,c.B*256,c.G*256};

        TXNTypeAttributes  attr[1];

        attr[0].tag  = kTXNQDFontColorAttribute;
        attr[0].size = kTXNQDFontColorAttributeSize;
        attr[0].data.dataPtr = &newC;
        TXNSetTypeAttributes (txn,1,attr,kTXNStartOffset,kTXNEndOffset);


    }
#else
    if (te) {
        _HYTextBox *  parent = (_HYTextBox*)this;

        TextStyle   newStyle;
        newStyle.tsColor.red   = c.R*256;
        newStyle.tsColor.blue  = c.B*256;
        newStyle.tsColor.green = c.G*256;

        if ((*te)->teLength) {
            short       ss = (*te)->selStart,
                        se = (*te)->selEnd;


            TESetSelect (0,30000,te);
            TESetStyle  (doColor,&newStyle,true,te);
            TESetSelect (ss,se,te);
        } else {
            TESetStyle  (doColor,&newStyle,true,te);
        }

    }
#endif
}

//__________________________________________________________________

void    _HYPlatformTextBox::_SetFont (_HYFont& f)
{
#ifdef          _HY_USE_MLTE_
    if (txn) {
        Str255          fName;
        StringToStr255 (f.face,fName);
        //short             fNum=0;
        //GetFNum      (fName,&fNum);

        TXNTypeAttributes  attr[3];

        attr[0].tag  = kTXNQDFontNameAttribute;
        attr[0].size = kTXNQDFontNameAttributeSize;
        attr[0].data.dataPtr = (Ptr)fName;

        attr[1].tag  = kTXNQDFontStyleAttribute;
        attr[1].size = kTXNQDFontStyleAttributeSize;
        attr[1].data.dataValue = normal;
        if (f.style & HY_FONT_BOLD) {
            attr[1].data.dataValue |= bold;
        }
        if (f.style & HY_FONT_ITALIC) {
            attr[1].data.dataValue |= italic;
        }


        attr[2].tag  = kTXNQDFontSizeAttribute;
        attr[2].size = kTXNFontSizeAttributeSize;
        attr[2].data.dataValue = f.size<<16;

        TXNSetTypeAttributes (txn,3,attr,kTXNStartOffset,kTXNEndOffset);
    }
#else
    if (te) {
        Str255 fName;
        StringToStr255 (f.face,fName);
        short  fNum=0;
        GetFNum (fName,&fNum);

        TextStyle   newStyle;
        newStyle.tsFont = fNum;
        newStyle.tsFace = f.style;
        newStyle.tsSize = f.size;

        if ((*te)->teLength) {
            short       ss = (*te)->selStart,
                        se = (*te)->selEnd;

            TESetSelect (0,30000,te);
            TESetStyle  (doFont|doFace|doSize,&newStyle,false,te);
            TESetSelect (ss,se,te);
        } else {
            TESetStyle  (doFont|doFace|doSize,&newStyle,false,te);
        }
    }
#endif
}

//__________________________________________________________________
void        _HYPlatformTextBox::_Update (Ptr p)
{
    _Paint (p);
}

//__________________________________________________________________
void        _HYPlatformTextBox::_SetDimensions (_HYRect r, _HYRect rel)
{
    _HYTextBox* theParent = (_HYTextBox *) this;
    theParent->_HYPlatformComponent::_SetDimensions (r,rel);
    _SetVisibleSize (rel);
}


//__________________________________________________________________
void        _HYPlatformTextBox::_SetVisibleSize (_HYRect rel)
{
    _HYTextBox *theParent = (_HYTextBox*) this;

    textBoxRect.top     = rel.top   + theParent->margins.top;
    textBoxRect.left    = rel.left  + theParent->margins.left;
    textBoxRect.bottom  = rel.bottom- theParent->margins.bottom;
    textBoxRect.right   = rel.right - theParent->margins.right;


#ifdef          _HY_USE_MLTE_
    if (txn) {
        Rect        textBoxRect2 = textBoxRect;

        if (!(theParent->boxFlags & (HY_TB_BIGBOX))) {
            InsetRect   (&textBoxRect2,1,1);
        }

        TXNLongRect newLongRect = {textBoxRect2.top, textBoxRect2.left, textBoxRect2.bottom, textBoxRect2.right};

        if (theParent->boxFlags & (HY_TB_BIGBOX|HY_TB_WRAP)) {
            newLongRect.bottom += 0x0fffffff;
            newLongRect.right -= HY_SCROLLER_WIDTH;
        }

        TXNSetRectBounds (txn, &textBoxRect2, &newLongRect, true);

        if (theParent->boxFlags & (HY_TB_BIGBOX|HY_TB_WRAP)) {
            TXNForceUpdate (txn);
        }
    }
#else
    if ((textBoxRect.bottom-textBoxRect.top<theParent->editBoxFont.size+4)||
            (textBoxRect.right-textBoxRect.left<theParent->editBoxFont.size+4)) {
        if (te) {
            TEDispose (te);
            te = nil;
        }
        return;
    }

    if (te) {
        CharsHandle th  = TEGetText  (te);
        Rect        textBoxRect2 = textBoxRect;
        InsetRect   (&textBoxRect2,2,2);
        Rect        textBoxRect3 = textBoxRect2;
        if (!theParent->wrapToView) {
            textBoxRect3.right += 10000;
        }
        GrafPtr thisPort;
        GetPort (&thisPort);
#ifdef OPAQUE_TOOLBOX_STRUCTS
        SetPort (GetWindowPort(theParent->parentWindow));
#else
        SetPort (theParent->parentWindow);
#endif

        TEHandle nte = TEStyleNew (&textBoxRect3, &textBoxRect2);
        checkPointer(nte);
        SetPort (thisPort);
        HLock ((Handle)te);
        HLock (th);
        short       ss = (*te)->selStart,
                    se = (*te)->selEnd;
        TESetText   (*th,(*te)->teLength,nte);
        TESetSelect (ss,se,nte);
        HUnlock     (th);
        HUnlock     ((Handle)te);
        TEDispose   (te);
        te = nte;
        _CreateTE();
    }
#endif
}

//__________________________________________________________________
void        _HYPlatformTextBox::_SetAlignFlags (unsigned char f)
{
#ifdef          _HY_USE_MLTE_
    if (txn) {
        TXNControlTag  tags[1];
        TXNControlData td  [1];

        tags[0]  = kTXNJustificationTag;

        if (f&HY_ALIGN_LEFT) {
            td[0].uValue  = kTXNFlushLeft;
        } else if (f&HY_ALIGN_RIGHT) {
            td[0].uValue  = kTXNFlushRight;
        } else {
            td[0].uValue  = kTXNCenter;
        }

        TXNSetTXNObjectControls (txn,false,1,tags,td);

    }
#else
    if (te) {
        if (f&HY_ALIGN_LEFT) {
            TESetAlignment (teFlushLeft,te);
        } else if (f&HY_ALIGN_RIGHT) {
            TESetAlignment (teFlushRight,te);
        } else {
            TESetAlignment (teCenter,te);
        }
    }
#endif
}

//__________________________________________________________________
_String     _HYPlatformTextBox::_GetText (void)
{
#ifdef          _HY_USE_MLTE_
    if (txn) {
        if (_GetCharacterCount()) {
            Handle      th;
            TXNGetDataEncoded
            (txn,kTXNStartOffset,kTXNEndOffset,&th,kTXNTextData);

            HLock       (th);
            _String     result (GetHandleSize(th),false);
            memcpy      (result.sData, *th, result.sLength);
            HUnlock     (th);
            DisposeHandle (th);
            //printf        ("%s\n", result.sData);
            return      result;
        }
    }
#else
    if (te&&(*te)->teLength) {
        _String     result ((unsigned long)(*te)->teLength,false);
        CharsHandle th = TEGetText (te);
        HLock       (th);
        memcpy      (result.sData, *th, result.sLength);
        HUnlock     (th);
        return      result;
    }
#endif
    return empty;
}

//__________________________________________________________________
void        _HYPlatformTextBox::_StoreText (_String*& res, bool selOnly)
{
#ifdef          _HY_USE_MLTE_
    if (txn) {
        TXNOffset       s,e;
        if (selOnly) {
            TXNGetSelection (txn, &s, &e);
            if (e<s) {
                res = new _String;
                return;
            }
        } else if (!_GetCharacterCount()) {
            res = new _String;
            return;
        } else {
            s = kTXNStartOffset;
            e = kTXNEndOffset;
        }

        Handle      th;
        TXNGetDataEncoded  (txn,s,e,&th,kTXNTextData);

        HLock       (th);
        res         = new _String (GetHandleSize(th),false);
        memcpy      (res->sData, *th,res->sLength);
        res->sData[res->sLength] = 0;
        HUnlock     (th);
        DisposeHandle (th);
    }
#endif
}
//__________________________________________________________________
void        _HYPlatformTextBox::_CreateTE (void)
{
#ifdef          _HY_USE_MLTE_
    if (txn) {
        _HYTextBox * theParent = (_HYTextBox*)this;

        _SetFont        (theParent->editBoxFont);
        _SetForeColor   (theParent->foreColor);
        //TEAutoView    (true, te);
        // check this
        _SetAlignFlags  (theParent->alignFlags);
        TXNRecalcTextLayout (txn);

        if ((theParent->boxFlags & HY_TB_ENABLED)&&theParent->activationFlag) {
            TXNActivate (txn,txnFrame,true);
        }

        TXNDraw   (txn, nil);
    }
#else
    if (te) {
        _HYTextBox * theParent = (_HYTextBox*)this;

        _SetFont        (theParent->editBoxFont);
        _SetForeColor   (theParent->foreColor);
        TEAutoView      (true, te);
        TECalText       (te);
        _SetAlignFlags  (theParent->alignFlags);
        if (theParent->isEnabled&&theParent->activationFlag) {
            TEActivate (te);
        }
        Rect            textBoxRect2 = textBoxRect;
        InsetRect       (&textBoxRect2,2,2);
        TEUpdate        (&textBoxRect2,te);
    }
#endif
}

//__________________________________________________________________
void        _HYPlatformTextBox::_SetText (const _String& editBoxText)
{
#ifdef          _HY_USE_MLTE_
    if (!txn) {
        _HYTextBox * parent = ((_HYTextBox*)this);

        Rect        textBoxRect2 = textBoxRect;
        if (parent->boxFlags & (HY_TB_BIGBOX|HY_TB_WRAP)) {
            textBoxRect.right -= HY_SCROLLER_WIDTH;
        } else {
            InsetRect   (&textBoxRect2,1,1);
        }

        Rect        textBoxRect3 = textBoxRect2;


        if (!(parent->boxFlags & (HY_TB_BIGBOX|HY_TB_WRAP))) {
            textBoxRect3.right += 0x0fffffff;
        }

        OSStatus err;

        if (parent->boxFlags & HY_TB_BIGBOX)
            err =  TXNNewObject (nil, ((_HYTextBox*)this)->parentWindow, &textBoxRect3,
                                 kTXNDontDrawCaretWhenInactiveMask|kTXNDontDrawSelectionWhenInactiveMask|kTXNWantVScrollBarMask|kTXNAlwaysWrapAtViewEdgeMask,
                                 kTXNTextEditStyleFrameType,
                                 kTXNTextFile,
                                 kTXNUnicodeEncoding,
                                 &txn,
                                 &txnFrame,
                                 0);
        else
            err =  TXNNewObject (nil, ((_HYTextBox*)this)->parentWindow, &textBoxRect3,
                                 kTXNDontDrawCaretWhenInactiveMask|kTXNDontDrawSelectionWhenInactiveMask|((parent->boxFlags&HY_TB_WRAP)?(kTXNWantVScrollBarMask|kTXNAlwaysWrapAtViewEdgeMask):0),
                                 kTXNTextEditStyleFrameType,
                                 kTXNTextFile,
                                 kTXNUnicodeEncoding,
                                 &txn,
                                 &txnFrame,
                                 0);


        if (err!=noErr) {
            _String errMsg = _String ("System error ") & (long)err & " while trying to create MLTE object";
            WarnError (errMsg);
            txn = nil;
            return;
        }


        TXNCarbonEventInfo      carbonEventInfo;
        TXNControlTag           iControlTags[] = { kTXNUseCarbonEvents };
        TXNControlData          iControlData[1];
        iControlData[0].uValue = (UInt32) &carbonEventInfo;

        carbonEventInfo.useCarbonEvents = false;
        carbonEventInfo.filler = 0;
        carbonEventInfo.flags = 0;
        carbonEventInfo.fDictionary = NULL;
        iControlData[0].uValue = (UInt32) &carbonEventInfo;
        TXNSetTXNObjectControls(txn,
                                false,
                                1,
                                iControlTags,
                                iControlData
                               );

        TXNBackground    tbg;

        tbg.bgType   = kTXNBackgroundTypeRGB;
        tbg.bg.color = (RGBColor) {
            0xffff,0xffff,0xffff
        };

        err = TXNSetBackground (txn,&tbg);

        TXNSetViewRect (txn, & textBoxRect2);
    }

    TXNSetData (txn,kTXNTextData,editBoxText.sData,editBoxText.sLength,kTXNStartOffset,kTXNEndOffset);

#else
    if (!te) {
        _HYTextBox * theParent = (_HYTextBox*)this;
        Rect        textBoxRect2 = textBoxRect;
        InsetRect   (&textBoxRect2,2,2);
        Rect        textBoxRect3 = textBoxRect2;
        if (!theParent->wrapToView) {
            textBoxRect3.right += 10000;
        }
        GrafPtr thisPort;
        GetPort (&thisPort);
#ifdef OPAQUE_TOOLBOX_STRUCTS
        SetPort (GetWindowPort(((_HYTextBox*)this)->parentWindow));
#else
        SetPort (((_HYTextBox*)this)->parentWindow);
#endif
        te = TEStyleNew (&textBoxRect3, &textBoxRect2);
        checkPointer(te);
        SetPort (thisPort);
    }
    TESetText   (editBoxText.sData,editBoxText.sLength,te);
#endif
    _CreateTE();
}

//__________________________________________________________________
void        _HYPlatformTextBox::_InsertText (const _String& editBoxText, bool append)
{
#ifdef          _HY_USE_MLTE_
    if (!txn) {
        _SetText (editBoxText);
    } else {
        //printf ("Insert text\n");
        TXNSetData (txn,kTXNTextData,editBoxText.sData,editBoxText.sLength,append?kTXNEndOffset:kTXNUseCurrentSelection,append?kTXNEndOffset:kTXNUseCurrentSelection);
    }
#else
    if (!te) {
        _SetText (editBoxText);
    } else {
        TEInsert    (editBoxText.sData,editBoxText.sLength,te);
    }
#endif
}

//__________________________________________________________________
void        _HYPlatformTextBox::_Paint (Ptr p)
{
    _HYTextBox * theParent = (_HYTextBox*)this;
    _HYRect    * relRect   = (_HYRect*)p;
    Rect         cRect;

    cRect.left   = relRect->left;
    cRect.right  = relRect->right;
    cRect.top    = relRect->top;
    if (theParent->margins.top) {
        cRect.bottom = relRect->top+theParent->margins.top;
        if (!(theParent->settings.width&HY_COMPONENT_TRANSP_BG)) {
            FillCRect (&cRect,backFill);
        } else {
            EraseRect (&cRect);
        }
    }
    cRect.bottom = relRect->bottom;
    cRect.top = relRect->bottom-theParent->margins.bottom;
    if (theParent->margins.bottom) {
        if (!(theParent->settings.width&HY_COMPONENT_TRANSP_BG)) {
            FillCRect (&cRect,backFill);
        } else {
            EraseRect (&cRect);
        }
    }
    cRect.top    = relRect->top;
    cRect.bottom = relRect->bottom;

    if (theParent->margins.left) {
        cRect.right  = relRect->left+theParent->margins.left;
        if (!(theParent->settings.width&HY_COMPONENT_TRANSP_BG)) {
            FillCRect (&cRect,backFill);
        } else {
            EraseRect (&cRect);
        }
    }

    if (theParent->margins.right) {
        cRect.left  = relRect->right-theParent->margins.right;
        cRect.right  = relRect->right;
        if (!(theParent->settings.width&HY_COMPONENT_TRANSP_BG)) {
            FillCRect (&cRect,backFill);
        } else {
            EraseRect (&cRect);
        }
    }

#ifdef          _HY_USE_MLTE_
    if (txn) {
        if (theParent->boxFlags & HY_TB_BIGBOX) {
            //printf ("Text box update\n");
            TXNDrawObject   (txn, nil, kTXNDrawItemAllMask);
        } else {
            Rect        textBoxRect2 = textBoxRect;

            InsetRect   (&textBoxRect2,1,1);

            DrawThemeEditTextFrame (&textBoxRect2,theParent->activationFlag?kThemeStateActive:kThemeStateInactive);

            TXNDraw   (txn, nil);

            if (theParent->activationFlag&&(theParent->boxFlags & HY_TB_FOCUSED)) {
                //InsetRect          (&textBoxRect2,1,1);
                DrawThemeFocusRect (&textBoxRect2,true);
            }
        }
    }
#else
    if (te) {
        FillCRect(&textBoxRect,backTFill);
        Rect        textBoxRect2 = textBoxRect;
        InsetRect   (&textBoxRect2,1,1);
        DrawThemeEditTextFrame (&textBoxRect2,theParent->activationFlag?kThemeStateActive:kThemeStateInactive);
        InsetRect   (&textBoxRect2,1,1);
        if (theParent->settings.width&HY_COMPONENT_TRANSP_BG) {
            SetThemeWindowBackground (theParent->parentWindow,kThemeBrushWhite,false);
            TEUpdate    (&textBoxRect2,te);
            SetThemeWindowBackground (theParent->parentWindow,kThemeBrushDialogBackgroundActive,false);
        } else {
            TEUpdate    (&textBoxRect2,te);
        }

        if (theParent->activationFlag&&theParent->isFocused) {
            InsetRect   (&textBoxRect2,1,1);
            DrawThemeFocusRect (&textBoxRect2,true);
        }
    }
#endif


    (*theParent)._HYPlatformComponent::_Paint(p);
}


//__________________________________________________________________

void        _HYPlatformTextBox::_EnableTextBox (bool e)
{
#ifdef          _HY_USE_MLTE_
    if (txn) {
        TXNActivate (txn,txnFrame,e);
    }

#else
    if (te)
        if (e) {
            TEActivate (te);
        } else {
            TEDeactivate (te);
        }
#endif
}

//__________________________________________________________________

void        _HYPlatformTextBox::_FocusComponent (void)
{
#ifdef          _HY_USE_MLTE_
    if (txn) {
        //_String ac = _String ("Focused ") & (long)((_HYTextBox*)this)->GetID() & '\n';
        //StringToConsole (ac);
        TXNFocus    (txn,true);
        TXNActivate (txn,txnFrame,true);
    }
#else
    if (te) {
        TEActivate (te);
    }
#endif
}

//__________________________________________________________________

void        _HYPlatformTextBox::_UnfocusComponent (void)
{
#ifdef          _HY_USE_MLTE_
    if (txn) {
        //_String ac = _String ("Unfocused ") & (long)((_HYTextBox*)this)->GetID() & '\n';
        //StringToConsole (ac);
        TXNFocus    (txn,false);
        TXNActivate (txn,txnFrame,false);
    }
#else
    if (te) {
        TEDeactivate (te);
    }
#endif
}

//__________________________________________________________________

void        _HYTextBox::_Activate (void)
{
    if (!activationFlag)
#ifdef          _HY_USE_MLTE_
        if (txn&&(boxFlags & HY_TB_ENABLED)) {
            //_String ac = _String ("Activated ") & (long)GetID() & '\n';
            //StringToConsole (ac);
            TXNActivate (txn,txnFrame,true);
        }
#else
        if (te&&isEnabled) {
            TEActivate (te);
        }
#endif

    _HYPlatformComponent::_Activate();
}

//__________________________________________________________________

void        _HYTextBox::_IdleHandler (void)
{
#ifdef          _HY_USE_MLTE_
    if (txn&&(boxFlags & HY_TB_FOCUSED)) {
        TXNIdle (txn);
    }
#else
    if (te&&isFocused) {
        TEIdle (te);
    }
#endif
}

//__________________________________________________________________

void        _HYTextBox::_Deactivate (void)
{
    if (activationFlag)
#ifdef          _HY_USE_MLTE_
        if (txn&&(boxFlags & HY_TB_ENABLED)) {
            //_String ac = _String ("Deactivated ") & (long)GetID() & '\n';
            //StringToConsole (ac);
            TXNActivate (txn,txnFrame,false);
        }
#else
        if (te&&isEnabled) {
            TEDeactivate (te);
        }
#endif

    _HYPlatformComponent::_Deactivate();
}

//__________________________________________________________________

_String mlteCommandString (TXNObject& txn, bool undoOrRedo)
{
    TXNActionKey tak;

    if (undoOrRedo) {
        if (!TXNCanUndo (txn,&tak)) {
            return empty;
        }
    } else if (!TXNCanRedo (txn,&tak)) {
        return empty;
    }

    switch (tak) {
    case kTXNTypingAction:
        return "Typing";

    case kTXNCutAction:
        return "Cut";

    case kTXNPasteAction:
        return "Paste";

    case kTXNClearAction:
        return "Clear";

    case kTXNDropAction:
        return "Drop";
    }

    return empty;
}


//__________________________________________________________________

void _HYPlatformTextBox::_HandleContextMenu (Point& where)
{
    _HYTextBox* parent = (_HYTextBox*)this;
    _List       menuOptions;
    _String     command;

    command = mlteCommandString (txn, true);
    if (command.sLength) {
        command = _String ("Undo ") & command;
    } else {
        command = "(Can't Undo";
    }

    menuOptions && & command;

    command = mlteCommandString (txn, false);
    if (command.sLength) {
        command = _String ("Redo ") & command;
    } else {
        command = "(Can't Redo";
    }

    menuOptions && & command;
    menuOptions && & menuSeparator;

    if (!TXNIsSelectionEmpty(txn)) {
        command = "Cut";
        menuOptions && & command;
        command = "Copy";
    } else {
        command = "(Cut";
        menuOptions && & command;
        command = "(Copy";
    }

    menuOptions && & command;

    if (TXNIsScrapPastable()) {
        command = "Paste";
    } else {
        command = "(Paste";
    }

    menuOptions && & command;

    if (_GetCharacterCount()) {
        command = "Select All";
    } else {
        command = "(Select All";
    }
    menuOptions && & command;

    command = HandlePullDown (menuOptions, where.h, where.v, -1);

    bool      sendMessage = true;

    long      f = menuOptions.Find (&command);

    if (f>=0) {
        switch (f) {
        case 0:
            parent->_DoUndo();
            break;
        case 1:
            parent->_DoRedo();
            break;
        case 3:
            parent->_DoCut();
            break;
        case 4:
            parent->_DoCopy();
            //sendMessage = false;
            break;
        case 5:
            parent->_DoPaste();
            break;
        case 6:
            parent->_DoSelectAll ();
            //sendMessage = false;
            break;
        default:
            return;

        }

        if (((_HYTextBox*)this)->messageRecipient&&sendMessage) {
            ((_HYTextBox*)this)->messageRecipient->ProcessEvent (generateTextEditChangeEvent (((_HYTextBox*)this)->GetID(),1));
        }
    }
}

//__________________________________________________________________

void _HYTextBox::_DoUndo (bool m)
{
#ifdef          _HY_USE_MLTE_
    TXNUndo (txn);
    if (m&&messageRecipient) {
        messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),1));
    }
#endif
}

//__________________________________________________________________

void _HYTextBox::_DoRedo (bool m)
{
#ifdef          _HY_USE_MLTE_
    TXNRedo (txn);
    if (m&&messageRecipient) {
        messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),1));
    }
#endif
}

//__________________________________________________________________

void _HYTextBox::_DoCut (bool m)
{
#ifdef          _HY_USE_MLTE_
    TXNCut(txn);
    if (m&&messageRecipient) {
        messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),1));
    }
#endif
}

//__________________________________________________________________

void _HYTextBox::_DoCopy (bool )
{
#ifdef          _HY_USE_MLTE_
    TXNCopy(txn);
#endif
}

//__________________________________________________________________

void _HYTextBox::_DoPaste (bool m)
{
#ifdef          _HY_USE_MLTE_
    TXNPaste(txn);
    if (m&&messageRecipient) {
        messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),1));
    }
#endif
}

//__________________________________________________________________

void _HYTextBox::_DoSelectAll (bool m)
{
    ((_HYTextBox*)this)->SetSelection (0,0x7fffffff);
    if (m&&messageRecipient) {
        messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),0));
    }
}

//__________________________________________________________________

void _HYTextBox::_DoClear (bool doAll, bool m)
{
#ifdef          _HY_USE_MLTE_
    if (doAll) {
        _DoSelectAll();
    }
    TXNClear(txn);
    if (doAll) {
        _SetFont (editBoxFont);
    }
    if (m&&messageRecipient) {
        messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),1));
    }

#endif
}

//__________________________________________________________________

void _HYTextBox::_DoFind (_String & st)
{
#ifdef          _HY_USE_MLTE_
    if (txn && (boxFlags & HY_TB_BIGBOX)) {
        TXNMatchTextRecord tmr;
        tmr.iTextPtr            = st.sData;
        tmr.iTextToMatchLength  = st.sLength;
        tmr.iTextEncoding       = CreateTextEncoding (kTextEncodingMacRoman,
                                  kTextEncodingDefaultVariant,
                                  kTextEncodingDefaultFormat);

        TXNOffset start = 0,
                  end;

        TXNGetSelection (txn, &start, &end);
        end   = kTXNEndOffset;

        TXNFind (txn,
                 &tmr,
                 kTXNTextData|kTXNUnicodeTextData,
                 0,
                 start,
                 end,
                 nil,
                 0,
                 &start,
                 &end);

        if ((start == end)&&(start == kTXNUseCurrentSelection)) {
            _String eMsg = _String ("Could not find '") & st & "'.";
            ProblemReport (eMsg, (Ptr)this);
        } else {
            TXNSetSelection (txn, start, end);
            TXNShowSelection (txn, false);
        }

    }
#endif
}


//__________________________________________________________________

bool _HYTextBox::_ProcessOSEvent (Ptr vEvent)
{
    EventRecord*    theEvent = (EventRecord*)vEvent;
    WindowPtr       dummy;
#ifdef          _HY_USE_MLTE_
    //printf ("%d\n",theEvent->what);
    if (txn) {
        switch (theEvent->what) {
        case mouseDown: {
            long evtType = FindWindow (theEvent->where,&dummy);
            if (evtType == inContent) {
                if (boxFlags & HY_TB_FOCUSED) {
                    Point localClick = theEvent->where;
                    if (theEvent->modifiers&controlKey) {
                        _HandleContextMenu (localClick);
                    } else {
                        GlobalToLocal (&localClick);
                        if (PtInRect (localClick,&textBoxRect)) {
                            TXNClick (txn, theEvent);
                            if (messageRecipient) {
                                messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),0));
                            }
                        }
                    }
                } else if (messageRecipient) {
                    messageRecipient->ProcessEvent(generateKeyboardFocusEvent(GetID()));
                }
            }
            return true;
        }
        case keyDown:
        case autoKey: {
            if (boxFlags & HY_TB_FOCUSED) {

                char c = theEvent->message&charCodeMask,
                     k = (theEvent->message&keyCodeMask)>>8;

                bool sendMessage = true;

                if ((((_HYTextBox*)this)->boxFlags&HY_TB_ARROWS) && messageRecipient) {
                    if (k==125) {
                        messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),3));
                        return true;
                    } else if (k==126) {
                        messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),4));
                        return true;
                    }
                }

                if (!(((_HYTextBox*)this)->boxFlags&HY_TB_BIGBOX)) {
                    if (theEvent->modifiers&cmdKey) {
                        bool done = true;
                        if ((c=='c')||(c=='C')) {
                            TXNCopy (txn);
                            sendMessage = false;
                        } else if ((c=='x')||(c=='X')) {
                            TXNCut (txn);
                        } else if ((c=='v')||(c=='V')) {
                            if (TXNIsScrapPastable()) {
                                TXNPaste (txn);
                            }
                        } else if ((c=='a')||(c=='A')) {
                            SetSelection (0,0x7fffffff);
                            sendMessage = false;
                        } else if ((c=='z')||(c=='Z')) {
                            if (mlteCommandString (txn,true).sLength) {
                                TXNUndo (txn);
                            }
                        } else if ((c=='y')||(c=='Y')) {
                            if (mlteCommandString (txn,false).sLength) {
                                TXNRedo (txn);
                            }
                        } else {
                            done = false;
                        }
                    } else {
                        if ((k==0x24)||(k==0x4C)) {
                            if (messageRecipient) {
                                messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),2));
                            }
                            return true;
                        } else {
                            TXNKeyDown (txn,theEvent);
                        }
                    }
                } else {
                    if (theEvent->modifiers&cmdKey) {
                        return false;
                    }

                    TXNKeyDown (txn,theEvent);
                }
                if (messageRecipient&&sendMessage) {
                    messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),1));
                }
                return true;
            }
        }
        }
    }
#else
    if (te) {
        switch (theEvent->what) {
        case mouseDown: {
            long evtType = FindWindow (theEvent->where,&dummy);
            if (evtType == inContent) {
                if (isFocused) {
                    Point localClick = theEvent->where;
                    GlobalToLocal (&localClick);
                    if (PtInRect (localClick,&textBoxRect)) {
                        if (settings.width&HY_COMPONENT_TRANSP_BG) {
                            SetThemeWindowBackground (parentWindow,kThemeBrushWhite,false);
                        }
                        TEClick (localClick, theEvent->modifiers&shiftKey, te);
                        if (settings.width&HY_COMPONENT_TRANSP_BG) {
                            SetThemeWindowBackground (parentWindow,kThemeBrushDialogBackgroundActive,false);
                        }
                        if (messageRecipient) {
                            messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),0));
                        }
                    }
                } else if (messageRecipient) {
                    messageRecipient->ProcessEvent(generateKeyboardFocusEvent(GetID()));
                }
            }
            return true;
        }
        case keyDown:
        case autoKey: {
            if (isFocused) {
                if (settings.width&HY_COMPONENT_TRANSP_BG) {
                    SetThemeWindowBackground (parentWindow,kThemeBrushWhite,false);
                }

                char c = theEvent->message&charCodeMask,
                     k = (theEvent->message&keyCodeMask)>>8;

                bool sendMessage = true;
                if (theEvent->modifiers&cmdKey) {
                    if ((c=='c')||(c=='C')) {
                        TECopy (te);
                        sendMessage = false;
                    } else if ((c=='x')||(c=='X')) {
                        TECut (te);
                    } else if ((c=='v')||(c=='V')) {
                        TEPaste (te);
                    } else if ((c=='a')||(c=='A')) {
                        SetSelection (0,30000);
                        sendMessage = false;
                    }
                } else {
                    if ((k==0x24)||(k==0x4C)) {
                        if (messageRecipient) {
                            messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),2));
                        }
                        return true;
                    } else {
                        TEKey (c,te);
                    }
                }
                if (settings.width&HY_COMPONENT_TRANSP_BG) {
                    SetThemeWindowBackground (parentWindow,kThemeBrushDialogBackgroundActive,false);
                }
                if (messageRecipient&&sendMessage) {
                    messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),1));
                }
                return true;
            }
        }
        }
    }
#endif
    return _HYPlatformComponent::_ProcessOSEvent (vEvent);
}

//__________________________________________________________________

void    _HYTextBox::_SetSelection (long s, long e)
{
#ifdef          _HY_USE_MLTE_
    if (txn&&(s>=0)&&(e>=s)) {
        TXNSetSelection (txn,s,e);
    }
#else
    if (te&&(s>=0)&&(e>=s)) {
        TESetSelect (s,e,te);
    }
#endif
}

//__________________________________________________________________

bool    _HYTextBox::_IsEmpty (void)
{
#ifdef          _HY_USE_MLTE_
    return _GetCharacterCount() == 0;
#else
    return !_GetText().sLength;
#endif
}

//__________________________________________________________________

void    _HYTextBox::_MarkForUpdate(void)
{
    _HYPlatformComponent::_MarkForUpdate();
}





//EOF