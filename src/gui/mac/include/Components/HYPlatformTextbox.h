/*

    A text input box object for MacOS.



    Sergei L. Kosakovsky Pond, May 2000-December 2002.

*/



#ifndef _HYPLTEXTBOX_

#define _HYPLTEXTBOX_



#ifdef   TARGET_API_MAC_CARBON

#include "CarbonEvents.h"

#endif



#include "HYPlatformComponent.h"

#include <MacTextEditor.h>



#define  _HY_USE_MLTE_



//__________________________________________________________________



class _HYPlatformTextBox

{



public:



    _HYPlatformTextBox (void);



    virtual ~_HYPlatformTextBox(void);



    virtual void            _SetBackColor    (_HYColor&);

    virtual void            _SetBackTColor   (_HYColor&);

    virtual void            _SetForeColor    (_HYColor&);

    virtual void            _SetDimensions   (_HYRect,_HYRect);

    virtual void            _SetVisibleSize  (_HYRect);

    virtual void            _SetFont         (_HYFont&);

    void            _SetText         (const _String&);

    void            _InsertText      (const _String&, bool);

    _String         _GetText         (void);

    void            _StoreText       (_String*&, bool);

    void            _CreateTE        (void);

    void            _EnableTextBox   (bool);

    void            _SetMargins      (_HYRect&);

    void            _SetAlignFlags   (unsigned char);

    long            _GetCharacterCount

    (void);



    virtual void            _Paint (Ptr p);

    virtual void            _Update(Ptr p);



    virtual void            _FocusComponent   (void);

    virtual void            _UnfocusComponent (void);



    PixPatHandle    backFill,

                    backTFill;





    Rect            textBoxRect;

#ifdef          _HY_USE_MLTE_



    virtual void            _HandleContextMenu (Point&);



    TXNObject       txn;

    TXNFrameID      txnFrame;

#else

    TEHandle        te;

#endif



};



#ifdef   _HY_USE_MLTE_



_String  mlteCommandString (TXNObject&, bool);



#endif



#endif

