#ifndef     __MATHOBJ__
#define     __MATHOBJ__

#include "baseobj.h"
#include "defines.h"
#include "errorfns.h"
#include "hy_lists.h"
#include "hy_strings.h"

class   _MathObject : public BaseObj  //abstract math operations class
{

public:

    virtual _MathObject* Add        (_MathObject*)     {
        warnError (-666); ;
        return new _MathObject;
    }
    virtual _MathObject* Sub        (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Minus      (void)             {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Sum        (void)             {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Mult       (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Div        (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* lDiv       (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* longDiv    (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Raise      (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual void         Assign     (_MathObject*)     {
        warnError (-666);
    }
    virtual bool         Equal      (_MathObject*)     {
        warnError (-666);
        return false;
    }
    virtual _MathObject* Abs        (void)             {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Sin        (void)             {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Cos        (void)             {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Tan        (void)             {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Exp        (void)             {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Log        (void)             {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Sqrt       (void)             {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Gamma      (void)             {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Erf        (void)             {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* LnGamma    (void)             {
        warnError (-666);    // <-- added by afyp, February 7, 2007
        return new _MathObject;
    }
    virtual _MathObject* Beta       (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* IGamma     (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* CChi2      (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* IBeta      (_MathObject*,_MathObject*) {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Simplex    (void)             {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Min        (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Max        (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* InvChi2    (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* ZCDF       (void)             {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Time       (void)             {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Arctan     (void)             {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Less       (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Random     (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Greater    (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* LessEq     (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* GreaterEq  (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* AreEqual   (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* NotEqual   (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* LAnd       (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* LOr        (_MathObject*)     {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* GammaDist  (_MathObject*,_MathObject*) {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* CGammaDist (_MathObject*,_MathObject*) {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* LNot       (void)             {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* TipCount   (void)             {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* BranchCount (void)            {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* TipName     (_MathObject*)    {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* BranchName  (_MathObject*)    {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* BranchLength(_MathObject*)    {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* RerootTree  (_MathObject*)    {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* TEXTreeString(_MathObject*) {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* Type                          (void);
    virtual _MathObject* PlainTreeString(_MathObject*) {
        warnError (-666);
        return new _MathObject;
    }
    virtual _MathObject* FormatNumberString (_MathObject*,_MathObject*) {
        warnError (-666);
        return new _MathObject;
    }
    virtual _Parameter   Value (void)              {
        warnError (-666);
        return 0.0;
    }
    virtual _MathObject* Compute (void)            {
        return new _MathObject;
    }
    virtual void         ScanForVariables (_AVLList&,bool = false)
    {}

    virtual      BaseRef makeDynamic               (void);
    virtual bool         IsVariable (void)         {
        return false;
    }
    virtual bool         IsObjectEmpty (void)      {
        return true;
    }
    virtual bool         IsPrintable (void)        {
        return false;
    }

    virtual bool         IsIndependent (void)       {
        return true;
    }
    virtual long         ObjectClass (void)         {
        return HY_UNDEFINED;
    }
    // returns a unique ID for this object
    // 0 - undefined
    // 1 - number
    // 4 - matrix

    virtual _MathObject* Execute (long opCode, _MathObject* p = nil , _MathObject* p2 = nil);
    // execute this operation with the list of Args

    virtual bool         HasChanged (void) {
        return false;
    }

    virtual   bool       IsConstant (void) {
        return true;
    }
};

// pointer to a math object
typedef _MathObject* _PMathObj ;

#endif
