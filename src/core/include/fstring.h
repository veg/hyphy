#ifndef     __FSTRING__
#define     __FSTRING__

#include "mathobj.h"

//__________________________________________________________________________________

class _FString : public _MathObject   // strings encountered in formulas
{

public:

    _FString (_String&, bool = true);
    _FString (long);
    _FString (_String*);
    _FString (void);
    virtual  ~_FString ();
//  ~_Constant (void);

    virtual BaseRef   makeDynamic       (void);
    virtual void      Duplicate         (BaseRef);
    virtual _PMathObj Add               (_PMathObj);
    virtual long      AddOn             (_PMathObj);
    virtual _PMathObj AreEqual          (_PMathObj);
    virtual _PMathObj AreEqualCIS       (_PMathObj);
    virtual _PMathObj Less              (_PMathObj);
    virtual _PMathObj LessEq            (_PMathObj);
    virtual _PMathObj Greater           (_PMathObj);
    virtual _PMathObj GreaterEq         (_PMathObj);
    virtual _PMathObj NotEqual          (_PMathObj);
    virtual _PMathObj RerootTree        (void);
    virtual _PMathObj EqualAmb          (_PMathObj);
    virtual _PMathObj EqualRegExp       (_PMathObj,bool = false);
    virtual _PMathObj ReplaceReqExp     (_PMathObj);
    virtual _PMathObj CountGlobalObjects(void);
    virtual _PMathObj FileExists        (void);
    virtual _PMathObj Evaluate          (void);
    virtual _PMathObj Join              (_PMathObj);
    virtual _PMathObj Differentiate     (_PMathObj);
    virtual long      ObjectClass       (void) {
        return STRING;
    }
    virtual _PMathObj Compute           (void) {
        return this;
    }

    virtual _PMathObj MapStringToVector (_PMathObj);
    virtual _PMathObj CharAccess        (_PMathObj,_PMathObj);
    virtual _PMathObj Execute           (long opCode, _MathObject* p = nil , _MathObject* p2 = nil);
    virtual BaseRef   toStr             (void);

    virtual bool      IsVariable        (void) {
        return true;
    }

    virtual bool      HasChanged        (void) {
        return true;
    }

    virtual bool      IsEmpty           (void) {
        return !theString || theString->sLength == 0;
    }
    // SLKP 20100907: a simple utility function to check if the object is an empty string

    _String*          theString;

};

