/*

    A general label object.



    Sergei L. Kosakovsky Pond, May 2000.



    December 2001, added _HYCheckbox

*/



#ifndef _HYLABEL_

#define _HYLABEL_

//#pragma once

#define  HY_LABEL_SHADOW 1



#include "HYComponent.h"

#include "HYPlatformLabel.h"

#include "HYPlatformCheckbox.h"



//__________________________________________________________________



class _HYLabel: public _HYComponent, public _HYPlatformLabel

{



public:



    _HYLabel(_HYRect,Ptr);



    virtual ~_HYLabel(void);



    virtual void        SetBackColor (_HYColor);

    virtual _HYColor&   GetBackColor (void);

    virtual void        SetForeColor (_HYColor);

    virtual _HYColor&   GetForeColor (void);

    virtual _String&    GetText (void);

    virtual void        SetFont (_HYFont&);

    virtual _HYFont&    GetFont (void);

    virtual void        SetText (_String);

    virtual void        Paint (Ptr p) {

        _HYPlatformLabel::_Paint(p);

    }

    virtual void        Update(Ptr p) {

        _HYPlatformLabel::_Update(p);

    }

    virtual void        SetDimensions(_HYRect r, _HYRect rel) {

        _HYComponent::SetDimensions(r,rel);

        _HYPlatformLabel::_SetDimensions(r,rel);

    }

    virtual void        SetVisibleSize(_HYRect);

    virtual _HYRect     _SuggestDimensions (void);

    virtual void        SetAlignFlags (unsigned char f) {

        alignFlags = f;

    }

    virtual unsigned char

    GetAlignFlags (void) {

        return alignFlags;

    }

    void        SetShadow(bool);

    bool        HasShadow(void) {

        return visFlags&HY_LABEL_SHADOW;

    }

private:



    _String         labelText;

    _HYColor        backColor, foreColor;

    _HYFont         labelFont;

    unsigned char   alignFlags, visFlags;

};



//__________________________________________________________________



class _HYCheckbox: public _HYLabel, public _HYPlatformCheckbox

{



public:



    _HYCheckbox(_HYRect,Ptr,bool = false);



    virtual ~_HYCheckbox(void);



    virtual bool        GetState   (void)                   {

        return checkState;

    }

    virtual long        GetSpacing (void)                   {

        return checkSpacing;

    }

    virtual void        SetState   (bool, bool = false);

    virtual void        SetSpacing (long);



    virtual void        Enable     (bool);



    virtual void        Paint (Ptr p) {

        _HYPlatformLabel::_Paint(p);

        _HYPlatformCheckbox::_Paint(p);

    }



    virtual void        Update(Ptr p)   {

        _HYPlatformLabel::_Update(p);

        _HYPlatformCheckbox::_Update(p);

    }



    virtual void        SetVisibleSize     (_HYRect);



    virtual _HYRect     _SuggestDimensions (void);

    virtual bool        _ProcessOSEvent     (Ptr);

    virtual void        _Activate           (void);

    virtual void        _Deactivate         (void);





    bool            checkState,

                    isEnabled;

    long            checkSpacing;

};



#endif



//EOF
