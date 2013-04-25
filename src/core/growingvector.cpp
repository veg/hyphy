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

#include "growingvector.h"

_GrowingVector::_GrowingVector(bool iscol) : _Matrix(64, 1, false, true) {
  used = 0;
  isColumn = iscol;
}

//______________________________________________________________________________
BaseRef _GrowingVector::makeDynamic(void) {
  _GrowingVector *result = (_GrowingVector *)checkPointer(new _GrowingVector);
  result->_Matrix::Duplicate(this);
  result->used = used;
  result->vDim = 1;
  result->isColumn = isColumn;
  return result;
}

//______________________________________________________________________________
long _GrowingVector::Store(_Parameter toStore) {
  if (used < hDim) {
    theData[used++] = toStore; // increment AFTER argument is sent to function
    return used - 1;
  } else {
    Resize(used + MAX(used / 8, 64)); // allocate another block of 64
    return Store(toStore);
  }
}

//______________________________________________________________________________
void _GrowingVector::operator<<(const _SimpleList &theSource) {
  for (long k = 0; k < theSource.lLength; k++) {
    Store(theSource.lData[k]);
  }
}

//______________________________________________________________________________
void _GrowingVector::Clear(void) {
  _Matrix::Clear();
  ZeroUsed();
  vDim = 1;
}

//______________________________________________________________________________
void _GrowingVector::Duplicate(BaseRef obj) {
  _Matrix::Duplicate(obj);
  used = ((_GrowingVector *)obj)->used;
  isColumn = ((_GrowingVector *)obj)->isColumn;
}
