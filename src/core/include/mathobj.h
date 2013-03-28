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

#ifndef     __MATHOBJ__
#define     __MATHOBJ__

#include "baseobj.h"
#include "defines.h"
#include "errorfns.h"
#include "list.h"
#include "avllistx.h"
#include "hy_strings.h"
#include "_hyExecutionContext.h"

class   _MathObject : public BaseObj  //abstract math operations class
{

public:

    virtual _MathObject* Add        (_MathObject*)     {
        warnError (-666); ;
        return this;
    }
    virtual _MathObject* Sub        (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Minus      (void)             {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Sum        (void)             {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Mult       (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Div        (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* lDiv       (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* longDiv    (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Raise      (_MathObject*)     {
        warnError (-666);
        return this;
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
        return this;
    }
    virtual _MathObject* Sin        (void)             {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Cos        (void)             {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Tan        (void)             {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Exp        (void)             {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Log        (void)             {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Sqrt       (void)             {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Gamma      (void)             {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Erf        (void)             {
        warnError (-666);
        return this;
    }
    virtual _MathObject* LnGamma    (void)             {
        warnError (-666);    // <-- added by afyp, February 7, 2007
        return this;
    }
    virtual _MathObject* Beta       (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* IGamma     (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* CChi2      (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* IBeta      (_MathObject*,_MathObject*) {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Simplex    (void)             {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Min        (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Max        (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* InvChi2    (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* ZCDF       (void)             {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Time       (void)             {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Arctan     (void)             {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Less       (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Random     (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Greater    (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* LessEq     (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* GreaterEq  (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* AreEqual   (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* NotEqual   (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* LAnd       (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* LOr        (_MathObject*)     {
        warnError (-666);
        return this;
    }
    virtual _MathObject* GammaDist  (_MathObject*,_MathObject*) {
        warnError (-666);
        return this;
    }
    virtual _MathObject* CGammaDist (_MathObject*,_MathObject*) {
        warnError (-666);
        return this;
    }
    virtual _MathObject* LNot       (void)             {
        warnError (-666);
        return this;
    }
    virtual _MathObject* TipCount   (void)             {
        warnError (-666);
        return this;
    }
    virtual _MathObject* BranchCount (void)            {
        warnError (-666);
        return this;
    }
    virtual _MathObject* TipName     (_MathObject*)    {
        warnError (-666);
        return this;
    }
    virtual _MathObject* BranchName  (_MathObject*)    {
        warnError (-666);
        return this;
    }
    virtual _MathObject* BranchLength(_MathObject*)    {
        warnError (-666);
        return this;
    }
    virtual _MathObject* RerootTree  (_MathObject*)    {
        warnError (-666);
        return this;
    }
    virtual _MathObject* TEXTreeString(_MathObject*) {
        warnError (-666);
        return this;
    }
    virtual _MathObject* Type                          (void);
    virtual _MathObject* PlainTreeString(_MathObject*) {
        warnError (-666);
        return this;
    }
    virtual _MathObject* FormatNumberString (_MathObject*,_MathObject*) {
        warnError (-666);
        return this;
    }
    virtual _Parameter   Value (void)              {
        warnError (-666);
        return 0.0;
    }
    virtual _MathObject* Compute (void)            {
        return this;
    }
    virtual void         ScanForVariables (_AVLList&,bool = false, _AVLListX* = nil, long = 0)
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
    virtual unsigned long  ObjectClass (void)       {
        return HY_UNDEFINED;
    }
    // returns a unique ID for this object
    // 0 - undefined
    // 1 - number
    // 4 - matrix

    virtual _MathObject* Execute (long opCode, _MathObject* p = nil , _MathObject* p2 = nil, _hyExecutionContext* context = _hyDefaultExecutionContext);
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
