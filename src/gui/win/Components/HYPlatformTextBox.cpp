/*
    Text input box for Win32 API

    Sergei L. Kosakovsky Pond, May 2000-December 2002
              Revised in December 2003 to use RichEdit
*/

#include "errorfns.h"
#include "HYTextbox.h"
#include "HYUtils.h"
#include "HYEventTypes.h"
#include "HYGraphicPane.h"
#include "HYPlatformWindow.h"

#include "HYDialogs.h"
//__________________________________________________________________

static LRESULT  CALLBACK editSubclassHandler (HWND WindowHand, UINT iMsg, WPARAM wParam, LPARAM lParam)
{
    _HYTextBox * theParent = (_HYTextBox*)GetWindowLongPtr (WindowHand, GWLP_USERDATA);

    if (iMsg == WM_CHAR) {
        int     keyCode = wParam;
        if (!(theParent->boxFlags & HY_TB_BIGBOX)) {
            if ((keyCode==VK_RETURN)||(keyCode==VK_ESCAPE)) {
                if (theParent->messageRecipient) {
                    theParent->messageRecipient->ProcessEvent (generateTextEditChangeEvent (theParent->GetID(),2));
                }

                SendMessage (theParent->parentWindow, WM_CHAR, wParam, lParam);
                return true;
            } else if (keyCode==VK_TAB) {
                SendMessage (theParent->parentWindow, WM_CHAR, wParam, lParam);
                return true;
            }
        }
    } else if (iMsg == WM_KEYDOWN) {
        if (theParent->messageRecipient && (theParent->boxFlags & HY_TB_ARROWS)) {
            if (wParam==VK_DOWN) {
                theParent->messageRecipient->ProcessEvent (generateTextEditChangeEvent (theParent->GetID(),3));
                //SendMessage (theParent->parentWindow, WM_CHAR, wParam, lParam);
                return true;
            } else if (wParam==VK_UP) {
                theParent->messageRecipient->ProcessEvent (generateTextEditChangeEvent (theParent->GetID(),4));
                //SendMessage (theParent->parentWindow, WM_CHAR, wParam, lParam);
                return true;
            }
        }
    }

    return CallWindowProc ((WNDPROC)theParent->mainHandler, WindowHand, iMsg, wParam, lParam);
}

//__________________________________________________________________

_HYPlatformTextBox::_HYPlatformTextBox  (void)
{
    backFill   = CreateSolidBrush (RGB (255,255,255));
    checkPointer ((Ptr)backFill);
    pLabelFont = nil;


    textBoxRect.left    = textBoxRect.top   = 0;
    textBoxRect.bottom  = textBoxRect.right = 100;
    textColor           = RGB(0,0,0);
    te                  = nil;

    if (((_HYTextBox*)this)->boxFlags & HY_TB_BIGBOX) {
        backTFill           = nil;
        isSingleLine        = false;
        te = CreateWindowEx (0,RICHEDIT_CLASS,"",ES_MULTILINE|WS_CHILD,30000,30000,100,100,((_HYTextBox*)this)->parentWindow,NULL,ProgramInstance, NULL);
        if (!te) {
            long    errCode = GetLastError();
            LPVOID lpMsgBuf;

            if (FormatMessage( FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
                               NULL,errCode,MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPTSTR) &lpMsgBuf,0,NULL )) {
                FlagError ((LPTSTR) lpMsgBuf);
            }
        }
        SendMessage (te,EM_SETUNDOLIMIT,10,0);
        SendMessage (te,EM_SETEVENTMASK,0,ENM_CHANGE|ENM_SELCHANGE);
        SendMessage (te,EM_SETBKGNDCOLOR, 0, (LPARAM)RGB(255,255,255));
    } else {
        backTFill  = CreateSolidBrush (RGB (255,255,255));
        isSingleLine        = true;
        te = CreateWindow ("EDIT","",/*WS_VISIBLE|*/WS_CHILD,30000,30000,100,100,((_HYTextBox*)this)->parentWindow,NULL,ProgramInstance, NULL);
    }

    _HYTextBox * theParent = (_HYTextBox*)this;
    SetWindowLongPtr (te,GWLP_USERDATA, (LONG_PTR)theParent);
    mainHandler = SetWindowLongPtr (te,GWLP_WNDPROC,(LONG_PTR)editSubclassHandler);
}

//__________________________________________________________________

_HYPlatformTextBox::~_HYPlatformTextBox (void)
{
    DeleteObject (backFill);
    DeleteObject (backTFill);
    DeleteObject (pLabelFont);
}

//__________________________________________________________________

void    _HYPlatformTextBox::_SetBackColor (_HYColor& c)
{
    DeleteObject (backFill);
    backFill = CreateSolidBrush (RGB(c.R,c.G,c.B));
    checkPointer (backFill);
}

//__________________________________________________________________

void    _HYPlatformTextBox::_SetBackTColor (_HYColor& c)
{
    if (((_HYTextBox*)this)->boxFlags & HY_TB_BIGBOX) {
        SendMessage (te,EM_SETBKGNDCOLOR,0,(LPARAM)RGB(c.R,c.G,c.B));
    } else {
        DeleteObject (backTFill);
        backTFill = CreateSolidBrush (RGB(c.R,c.G,c.B));
        checkPointer (backTFill);
    }
}

//__________________________________________________________________

void    _HYPlatformTextBox::_SetForeColor (_HYColor& c)
{
    textColor = RGB(c.R,c.G,c.B);
}

//__________________________________________________________________

void    _HYPlatformTextBox::_SetFont (_HYFont& f)
{
    if (pLabelFont) {
        DeleteObject(pLabelFont);
    }

    pLabelFont = CreateFont (f.size,0,0,0,(f.style&HY_FONT_BOLD)?FW_BOLD:FW_NORMAL,f.style&HY_FONT_ITALIC,FALSE,FALSE,ANSI_CHARSET,OUT_DEFAULT_PRECIS,
                             CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,DEFAULT_PITCH|FF_DONTCARE,f.face.sData);

    if (te) {
        _CreateTE();
        if (((_HYTextBox *) this)->boxFlags & HY_TB_ENABLED) {
            ((_HYTextBox *) this)->Activate();
        }
    }
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

    if ((textBoxRect.bottom-textBoxRect.top<theParent->editBoxFont.size+4)||
            (textBoxRect.right-textBoxRect.left<theParent->editBoxFont.size+4)) {
        if (te && (! theParent->boxFlags & HY_TB_BIGBOX)) {
            DestroyWindow (te);
            te = nil;
        }
        return;
    }

    if (te) {
        bool          newSL = !_NeedMultiLines ();
        if (newSL!=isSingleLine) {
            _CreateTE     ();
        } else {
            SetWindowPos(te,nil, textBoxRect.left,textBoxRect.top,textBoxRect.right-textBoxRect.left+1,textBoxRect.bottom-textBoxRect.top+1,SWP_NOZORDER);
        }

        if (theParent->activationFlag) {
            ShowWindow(te,SW_SHOWNA);
        }
    }
    //else
    //_CreateTE ();

}

//__________________________________________________________________
void        _HYPlatformTextBox::_SetAlignFlags (unsigned char f)
{
    if (te)
        if (((_HYTextBox*)this)->boxFlags & HY_TB_BIGBOX) {
            PARAFORMAT pf;
            pf.cbSize = sizeof(PARAFORMAT);
            pf.dwMask = PFM_ALIGNMENT;
            switch (f) {
            case HY_ALIGN_LEFT:
                pf.wAlignment = PFA_LEFT;
                break;

            case HY_ALIGN_RIGHT:
                pf.wAlignment = PFA_RIGHT;
                break;

            default:
                pf.wAlignment = PFA_CENTER;
            }
            SendMessage (te,EM_SETPARAFORMAT,0,(LPARAM)&pf);
        } else {
            _CreateTE     ();
        }
}
//__________________________________________________________________

_String                 retrieveEditControlText (HWND te)
{
    if (te) {
        long stringLength = SendMessage (te,WM_GETTEXTLENGTH,0,0);
        if (stringLength) {
            _String     res (stringLength, false);
            SendMessage (te,WM_GETTEXT, stringLength+1,(LPARAM)res.sData);
            res.sData[stringLength] = 0;
            return      res;
        }
    }
    return empty;

}

//__________________________________________________________________

void                    retrieveEditControlText (HWND te, _String*& res)
{
    if (te) {
        long stringLength = SendMessage (te,WM_GETTEXTLENGTH,0,0);
        if (stringLength) {
            res = new _String(stringLength, false);
            SendMessage (te,WM_GETTEXT, stringLength+1,(LPARAM)res->sData);
        } else {
            res = new _String;
        }
    }
}


//__________________________________________________________________
_String     _HYPlatformTextBox::_GetText (void)
{
    if ((((_HYTextBox*)this)->boxFlags & HY_TB_BIGBOX)) {
        _String * r;
        _StoreText (r,false);
        return _String (r);
    } else {
        return retrieveEditControlText(te);
    }
}

//__________________________________________________________________
void    _HYPlatformTextBox::_StoreText (_String*& rec, bool selOnly)
{
    if ((((_HYTextBox*)this)->boxFlags & HY_TB_BIGBOX)) {
        if (selOnly) {
            long        start,
                        end;
            SendMessage (te,EM_GETSEL,(WPARAM)&start,(LPARAM)&end);
            if (end>=start) {
                rec = new _String (end-start+1,false);
                checkPointer (rec);
                rec->sData[SendMessage (te,EM_GETSELTEXT,0,(LPARAM)rec->sData)] = 0;
            } else {
                rec = new _String;
            }
        } else {
            GETTEXTLENGTHEX gtl = {GTL_NUMCHARS,CP_ACP};
            long l = SendMessage (te,EM_GETTEXTLENGTHEX,(WPARAM)&gtl,0);
            if (l) {
                rec = new _String (l,false);
                checkPointer (rec);
                TEXTRANGE tr = {{0,l},rec->sData};
                rec->sData[SendMessage (te,EM_GETTEXTRANGE,0,(LPARAM)&tr)] = 0;
            } else {
                rec = new _String;
            }
        }
    } else {
        retrieveEditControlText(te, rec);
    }
}

//__________________________________________________________________
bool        _HYPlatformTextBox::_NeedMultiLines (void)
{
    if (((_HYTextBox*)this)->boxFlags & HY_TB_BIGBOX) {
        return 1;
    } else {
        return (textBoxRect.bottom-textBoxRect.top+1>((_HYTextBox*)this)->editBoxFont.size * 2.5);
    }
}

//__________________________________________________________________
void        _HYPlatformTextBox::_CreateTE (void)
{
    _String      currentText;

    if (te) {
        currentText   = _GetText ();
        DestroyWindow (te);
    }

    _HYTextBox * theParent = (_HYTextBox*)this;

    DWORD        windowStyle = 0;

    if (theParent->alignFlags&HY_ALIGN_LEFT) {
        windowStyle = ES_LEFT;
    } else if (theParent->alignFlags&HY_ALIGN_RIGHT) {
        windowStyle = ES_RIGHT;
    } else {
        windowStyle = ES_CENTER;
    }

    if (_NeedMultiLines ()) {
        isSingleLine = false;
        windowStyle |= ES_MULTILINE|ES_WANTRETURN|ES_AUTOVSCROLL;
        if (!(theParent->boxFlags & (HY_TB_WRAP||HY_TB_BIGBOX))) {
            windowStyle |= ES_AUTOHSCROLL;
        }
        if (theParent->boxFlags & HY_TB_BIGBOX) {
            windowStyle |= WS_VSCROLL;
        }
    } else {
        isSingleLine = true;
        windowStyle |= ES_AUTOHSCROLL;
    }

    if ((theParent->boxFlags & HY_TB_ENABLED) == 0) {
        windowStyle |= ES_READONLY;
    }

    if (theParent->boxFlags & HY_TB_BIGBOX) {
        te = CreateWindowEx (0,RICHEDIT_CLASS,currentText.sData,WS_CHILD|windowStyle,
                             textBoxRect.left,textBoxRect.top,textBoxRect.right-textBoxRect.left+1,textBoxRect.bottom-textBoxRect.top+1
                             ,theParent->parentWindow,NULL,ProgramInstance, NULL);

        checkPointer   (te);
        SetWindowLongPtr (te,GWLP_USERDATA, (LONG_PTR)theParent);
        mainHandler = SetWindowLongPtr (te,GWLP_WNDPROC,(LONG_PTR)editSubclassHandler);

        SendMessage (te,EM_SETUNDOLIMIT,10,0);
        SendMessage (te,EM_SETEVENTMASK,0,ENM_CHANGE|ENM_SELCHANGE);
        SendMessage (te,EM_SETBKGNDCOLOR, 0, (LPARAM)RGB(255,255,255));

    } else {

        te = CreateWindow ("EDIT",currentText.sData,/*WS_VISIBLE|*/WS_CHILD|windowStyle,
                           textBoxRect.left,textBoxRect.top,textBoxRect.right-textBoxRect.left+1,textBoxRect.bottom-textBoxRect.top+1
                           ,theParent->parentWindow,NULL,ProgramInstance, NULL);

        checkPointer   (te);
        SetWindowLongPtr (te,GWLP_USERDATA, (LONG_PTR)theParent);
        mainHandler = SetWindowLongPtr (te,GWLP_WNDPROC,(LONG_PTR)editSubclassHandler);


        SendMessage (te, EM_SETMARGINS, EC_LEFTMARGIN, 3);
        SendMessage (te, EM_SETMARGINS, EC_RIGHTMARGIN, 3);
    }


    if (pLabelFont) {
        SendMessage (te, WM_SETFONT, (WPARAM)pLabelFont, 1);
    }


}


//__________________________________________________________________
void        _HYPlatformTextBox::_SetText (const _String& editBoxText)
{
    _HYTextBox * theParent = (_HYTextBox*)this;
    if (!te) {
        _CreateTE();
    }

    _HYGuiObject* stashRec = theParent->messageRecipient;
    theParent->SetMessageRecipient (nil);
    if (theParent->boxFlags & HY_TB_BIGBOX) {
        SETTEXTEX     stext = {ST_DEFAULT,CP_ACP};
        SendMessage   (te,EM_SETTEXTEX,(WPARAM)&stext,(LPARAM)editBoxText.sData);
    } else {
        SendMessage   (te,WM_SETTEXT,0,(LPARAM)editBoxText.sData);
    }
    theParent->SetMessageRecipient (stashRec);
}

//__________________________________________________________________
void        _HYPlatformTextBox::_InsertText (const _String& editBoxText, bool append)
{
    if (!te) {
        _SetText (editBoxText);
    } else {
        if (append) {
            SendMessage (te, EM_SETSEL, 0x7fffffff, 0x7ffffff);
        }

        SendMessage (te, EM_REPLACESEL, 0, (LPARAM)editBoxText.sData);

        if (append) {
            SendMessage (te, WM_VSCROLL, SB_BOTTOM, 0);
        }

    }
}

//__________________________________________________________________
void        _HYPlatformTextBox::_Paint (Ptr p)
{
    _HYTextBox * theParent = (_HYTextBox*)this;
    _HYRect    * relRect   = (_HYRect*)p;

    RECT         cRect;
    HDC          hdc       = (HDC)relRect->width;

    if (!(theParent->settings.width&HY_COMPONENT_TRANSP_BG)) {
        cRect.left   = relRect->left;
        cRect.right  = relRect->right;
        cRect.top    = relRect->top;
        cRect.bottom = relRect->top+theParent->margins.top;

        if (theParent->margins.top) {
            FillRect (hdc,&cRect,backFill);
        }


        cRect.bottom = relRect->bottom;
        cRect.top = relRect->bottom-theParent->margins.bottom;

        if (theParent->margins.bottom) {
            FillRect (hdc,&cRect,backFill);
        }

        cRect.top    = relRect->top;
        cRect.bottom = relRect->bottom;
        cRect.right  = relRect->left+theParent->margins.left;

        if (theParent->margins.left) {
            FillRect (hdc,&cRect,backFill);
        }

        cRect.left  = relRect->right-theParent->margins.right;
        cRect.right  = relRect->right;

        if (theParent->margins.right) {
            FillRect (hdc,&cRect,backFill);
        }
    }

    if (te) {
        if (theParent->boxFlags & HY_TB_ENABLED) {
            UpdateWindow (te);
        }
    }

    if ((theParent->boxFlags & HY_TB_BIGBOX) == 0) {
        cRect = textBoxRect;
        InflateRect (&cRect,2,2);
        DrawEdge (hdc,&cRect,EDGE_BUMP, BF_RECT);
    }

    (*theParent)._HYPlatformComponent::_Paint(p);
}


//__________________________________________________________________

void        _HYPlatformTextBox::_EnableTextBox (bool e)
{
    if (te)
        if (e) {
            SendMessage (te, EM_SETREADONLY, 0, 0);
        } else {
            SendMessage (te, EM_SETREADONLY, 1, 0);
        }
}

//__________________________________________________________________

void        _HYPlatformTextBox::_FocusComponent (void)
{
    if (te) {
        SetFocus (te);
    }
}

//__________________________________________________________________

void        _HYPlatformTextBox::_UnfocusComponent (void)
{
}

//__________________________________________________________________

void        _HYTextBox::_Activate (void)
{
    _HYPlatformComponent::_Activate();
    //ReportWarning (_String ("Activated Component ") & (long)GetID());
    if (te) {
        //ReportWarning (_String ("Component Alive "));
        ShowWindow (te, SW_SHOW);
    }
}

//__________________________________________________________________

void        _HYTextBox::_IdleHandler (void)
{
}

//__________________________________________________________________

void        _HYTextBox::_Deactivate (void)
{
}

//__________________________________________________________________

bool _HYTextBox::_ProcessOSEvent (Ptr vEvent)
{
    _HYWindowsUIMessage * theEvent = (_HYWindowsUIMessage *)vEvent;


    switch (theEvent->iMsg) {
    case WM_COMMAND: {
        if ((theEvent->lParam== 0 )&&(boxFlags & HY_TB_FOCUSED)) {
            switch (LOWORD (theEvent->wParam)) {
            case HY_WINDOW_MENU_ID_EDIT: { // undo
                SendMessage (te, EM_UNDO, 0, 0);
                return true;
            }
            case HY_WINDOW_MENU_ID_EDIT+1: { // copy
                SendMessage (te, WM_COPY, 0, 0);
                return true;
            }
            case HY_WINDOW_MENU_ID_EDIT+2: { // cut
                SendMessage (te, WM_CUT, 0, 0);
                return true;
            }
            case HY_WINDOW_MENU_ID_EDIT+3: { // cut
                SendMessage (te, WM_PASTE, 0, 0);
                return true;
            }
            case HY_WINDOW_MENU_ID_EDIT+5: { // select all
                if (boxFlags & HY_TB_BIGBOX) {
                    CHARRANGE cr = {0,1};
                    SendMessage (te, EM_EXSETSEL, 0, (LPARAM)&cr);
                } else {
                    SetSelection (0,30000);
                }
                return true;
            }
            }
            break;
        } else {
            WORD       wP = HIWORD(theEvent->wParam);
            if ((wP==EN_SETFOCUS)||(wP==EN_CHANGE)) {
                HWND    cbh = (HWND)theEvent->lParam;
                if ((cbh == te)&&(boxFlags & HY_TB_ENABLED)) {
                    if (messageRecipient) {
                        if (wP==EN_SETFOCUS) {
                            messageRecipient->ProcessEvent(generateKeyboardFocusEvent(GetID()));
                        } else {
                            messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),1));
                        }

                    }
                    return true;
                }
            }
        }
        break;
    }

    case WM_NOTIFY: {
        NMHDR * nS = (NMHDR*)theEvent->lParam;
        if (nS->hwndFrom==te && (boxFlags & HY_TB_ENABLED)) {
            if (messageRecipient && nS->code == EN_SELCHANGE) {
                messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),0));
            }
        }
        break;
    }

    case WM_CTLCOLOREDIT:
    case WM_CTLCOLORSTATIC: {
        if (((HWND)theEvent->lParam == te) && (boxFlags & HY_TB_BIGBOX == 0)) {
            HDC         hdc = (HDC)theEvent->wParam;
            SetTextColor (hdc,textColor);
            theEvent->res = (LRESULT)backTFill;
            return true;
        }
        break;
    }

    }
    return _HYPlatformComponent::_ProcessOSEvent (vEvent);
}

//__________________________________________________________________

void    _HYTextBox::_SetSelection (long s, long e)
{
    if (te&&(s>=0)&&(e>=s)) {
        SendMessage (te, EM_SETSEL, s, e);
        SendMessage (te, EM_SCROLLCARET, 0, 0);
    }

}

//__________________________________________________________________

void    _HYTextBox::_MarkForUpdate (void)
{
    if (te) {
        RECT tR;
        GetClientRect (te, &tR);
        InvalidateRect(te, &tR, true);
    }
    _HYPlatformComponent::_MarkForUpdate();
}


//__________________________________________________________________

bool    _HYTextBox::_IsEmpty (void)
{
    if (te) {
        return !SendMessage (te,WM_GETTEXTLENGTH,0,0);
    }
    return true;
}

//__________________________________________________________________

void _HYTextBox::_DoUndo (bool m)
{
    SendMessage (te, EM_UNDO, 0, 0);
    if (m&&messageRecipient) {
        messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),1));
    }
}

//__________________________________________________________________

void _HYTextBox::_DoRedo (bool m)
{
    SendMessage (te, EM_REDO, 0, 0);
    if (m&&messageRecipient) {
        messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),1));
    }
}

//__________________________________________________________________

void _HYTextBox::_DoCut (bool m)
{
    SendMessage (te, WM_CUT, 0, 0);
    if (m&&messageRecipient) {
        messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),1));
    }
}

//__________________________________________________________________

void _HYTextBox::_DoCopy (bool )
{
    SendMessage (te, WM_COPY, 0, 0);
}

//__________________________________________________________________

void _HYTextBox::_DoPaste (bool m)
{
    SendMessage (te, WM_PASTE, 0, 0);
    if (m&&messageRecipient) {
        messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),1));
    }
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
    if (doAll) {
        SendMessage (te,WM_SETTEXT,(WPARAM)0,(LPARAM)empty.sData);
    } else {
        SendMessage (te,EM_REPLACESEL,(WPARAM)true,(LPARAM)empty.sData);
    }
    if (m&&messageRecipient) {
        messageRecipient->ProcessEvent (generateTextEditChangeEvent (GetID(),1));
    }

}

//__________________________________________________________________

void _HYTextBox::_DoFind (_String & st)
{
    //if (te && (boxFlags & HY_TB_BIGBOX))
    //{

    //long match = SendMessage (te,
}


//EOF