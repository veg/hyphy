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

#include "likefunc.h"
#include "customfunction.h"
#include "hy_strings.h"
#include "formula.h"

extern long likeFuncEvalCallCount;

_CustomFunction::_CustomFunction(_String *arg) {
  _String body(*arg), errMsg;

  _FormulaParsingContext fpc(&errMsg, nil);

  if (Parse(&myBody, body, fpc) == HY_FORMULA_EXPRESSION) {
    _SimpleList myVars;
    {
      _AVLList al(&myVars);
      myBody.ScanFForVariables(al, true, false, false);
      al.ReorderList();
    }
    for (unsigned long k = 0; k < myVars.lLength; k++)
      if (LocateVar(myVars.lData[k])->IsIndependent()) {
        GetIndependentVars() << myVars.lData[k];
      }
  } else {
    WarnError(_String(
        "An invalid expression supplied for formula-based custom LF: '") &
              errMsg & '\'');
  }
}

//______________________________________________________________________________
_Parameter _CustomFunction::Compute(void) {

  likeFuncEvalCallCount++;

  _SimpleList *iv = &GetIndependentVars();

  for (unsigned long i = 0; i < iv->lLength; i++) {
    _Parameter result = GetIthIndependent(i);

    if (result < GetIthIndependentBound(i, true) ||
        result > GetIthIndependentBound(i, false)) {
      return -A_LARGE_NUMBER;

    }
  }

  _PMathObj res = myBody.Compute();
  if (res) {
    return res->Value();
  }
  return 0.0;
}
