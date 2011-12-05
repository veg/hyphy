#ifndef     __CONSTANT__
#define     __CONSTANT__

#include "parser.h"

class _Constant : public _MathObject   // a numerical constant
{

public:

    _Constant (_Parameter);
    _Constant (_String&);
    _Constant (void);
    ~_Constant (void) {}

    virtual _PMathObj Add           (_PMathObj);
    virtual _PMathObj Sub           (_PMathObj);
    virtual _PMathObj Minus         (void) ;
    virtual _PMathObj Sum           (void) ;
    virtual _PMathObj Mult          (_PMathObj);
    virtual _PMathObj Div           (_PMathObj);
    virtual _PMathObj lDiv          (_PMathObj);
    virtual _PMathObj longDiv       (_PMathObj);
    virtual _PMathObj Raise         (_PMathObj);
    virtual void      Assign        (_PMathObj);
    virtual bool      Equal         (_PMathObj);
    virtual _PMathObj Abs           (void);
    virtual _PMathObj Sin           (void);
    virtual _PMathObj Cos           (void);
    virtual _PMathObj Tan           (void);
    virtual _PMathObj Exp           (void);
    virtual _PMathObj Log           (void);
    virtual _PMathObj Sqrt          (void);
    virtual _PMathObj Time          (void);
    virtual _PMathObj Arctan        (void);
    virtual _PMathObj Gamma         (void);
    virtual _PMathObj LnGamma       (void);         /* <- added by afyp, February 8, 2007 */
    virtual _PMathObj Beta          (_PMathObj);
    virtual _PMathObj Min           (_PMathObj);
    virtual _PMathObj Max           (_PMathObj);
    virtual _PMathObj GammaDist     (_PMathObj,_PMathObj);
    virtual _PMathObj CGammaDist    (_PMathObj,_PMathObj);
    virtual _PMathObj IBeta         (_PMathObj,_PMathObj);
    virtual _PMathObj IGamma        (_PMathObj);
    virtual _PMathObj CChi2         (_PMathObj);
    virtual _PMathObj InvChi2       (_PMathObj);
    virtual _PMathObj Erf           (void);
    virtual _PMathObj ZCDF          (void);
    virtual _PMathObj Less          (_PMathObj);
    virtual _PMathObj Greater       (_PMathObj);
    virtual _PMathObj LessEq        (_PMathObj);
    virtual _PMathObj GreaterEq     (_PMathObj);
    virtual _PMathObj AreEqual      (_PMathObj);
    virtual _PMathObj NotEqual      (_PMathObj);
    virtual _PMathObj LAnd          (_PMathObj);
    virtual _PMathObj LOr           (_PMathObj);
    virtual _PMathObj LNot          ();
    virtual _PMathObj Random        (_PMathObj);
    virtual _Parameter
    Value       (void);
    virtual _PMathObj FormatNumberString
    (_PMathObj,_PMathObj);
    virtual _PMathObj Compute       (void) {
        return this;
    };

    virtual   void    Initialize            (void);
    virtual   void    Duplicate             (BaseRef);
    virtual   BaseRef makeDynamic           (void);
    virtual   BaseRef toStr                 (void);
    virtual   long    ObjectClass           (void) {
        return NUMBER;
    }
    virtual   void    SetValue              (_Parameter pl) {
        theValue = pl;
    }

public:

    _Parameter theValue;

};

