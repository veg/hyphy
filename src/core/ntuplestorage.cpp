/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2009
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon              (apoon@cfenet.ubc.ca)

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

#include "matrix.h"
#include "ntuplestorage.h"

//______________________________________________________________________________
_NTupleStorage::_NTupleStorage(unsigned long N, unsigned long K) {
  storageN = N;
  storageK = (K > N) ? MIN(N, 1) : K;
  unsigned long matrixDimension = 1;

  // now compute what dimension of the matrix will be required
  // handle the special cases first

  if (storageK) { 
    // not just the empty set
    // populate C_NK_Lookup
    for (long i2 = 0; i2 <= storageN; i2++) {
      // N choose 0 is always 1 for every N
      C_NK_Lookup << 1; 
    }
    for (long i = 1; i <= storageK; i++) {
      for (long filler = 0; filler < i; filler++, C_NK_Lookup << 0); 
      // N choose K where K>N is invalid
      C_NK_Lookup << 1;
      // K choose K is 1
      for (long j = i + 1; j <= storageN; j = j + 1) {
        // N choose K = N/(N-K) times (N-1) choose K
        C_NK_Lookup << C_NK_Lookup.lData[C_NK_Lookup.lLength - 1] * j / (j - i);
      }
    }
  }

  //for (long i=0; i<=storageK; i++)
  //  {
  //      long offset = i*(storageN+1);
  //      for (long j=0; j<=storageN; j++)
  //          printf ("(%d,%d) = %d\n", j,i, C_NK_Lookup.lData[offset+j]);
  //  }


  matrixDimension = C_NK_Lookup.lData[C_NK_Lookup.lLength - 1];
  CreateMatrix(this, 1, matrixDimension, false, true, false);
}

//______________________________________________________________________________
BaseRef _NTupleStorage::makeDynamic(void) {
  _NTupleStorage *copy = new _NTupleStorage;
  checkPointer(copy);
  copy->_Matrix::Duplicate(this);
  copy->storageN = storageN;
  copy->storageK = storageK;
  copy->C_NK_Lookup.Duplicate(&C_NK_Lookup);
  return copy;
}

//______________________________________________________________________________
bool _NTupleStorage::CheckKTuple(_SimpleList &kTuple) {
  if (kTuple.lLength == storageK) {
    if (storageK) {
      kTuple.Sort();
      for (long k = 0; k < kTuple.lLength; k++)
        if (kTuple.lData[k] < 0 || kTuple.lData[k] >= storageN ||
            (k && kTuple.lData[k] == kTuple.lData[k - 1])) {
          return false;
        }
    }
    return true;
  }
  return false;
}

//______________________________________________________________________________
unsigned long _NTupleStorage::Index(_SimpleList &kTuple) {
  unsigned long myIndex = 0;
  if (storageK)
    for (long k = kTuple.lLength - 1; k >= 0; k--) {
      myIndex += C_NK_Lookup.lData[(k + 1) * (1 + storageN) + kTuple.lData[k]];
    }
  return myIndex;
}

//______________________________________________________________________________
_Parameter _NTupleStorage::DirectIndex(unsigned long directIndex) {
  return theData[directIndex];
}

//______________________________________________________________________________
unsigned long _NTupleStorage::Store(_Parameter value, _SimpleList &kTuple) {
  unsigned long myIndex = Index(kTuple);
  theData[myIndex] = value;
  return myIndex;
}

//______________________________________________________________________________
_Parameter _NTupleStorage::Retrieve(_SimpleList &kTuple) {
  unsigned long myIndex = Index(kTuple);
  return theData[myIndex];
}

//______________________________________________________________________________
void _NTupleStorage::IndexToTuple(unsigned long directIndex,
                                  _SimpleList &kTuple) {

  kTuple.Clear();
  if (storageK && directIndex < C_NK_Lookup.lData[C_NK_Lookup.lLength - 1]) {
    long currentN = storageN - 1;
    for (long k = storageK; k > 0; k--) {
      long i = currentN, lookup_offset = k * (storageN + 1);

      while (C_NK_Lookup.lData[lookup_offset + i] > directIndex) {
        //printf ("(%d, %d) -> %d\n", k, i, C_NK_Lookup.lData[lookup_offset+i]);
        i--;
      }

      //printf ("(%d, %d) -> %d\n", k, i, C_NK_Lookup.lData[lookup_offset+i]);
      kTuple << i;
      //printf ("Stored %d\n", i);

      currentN = i - 1;
      directIndex -= C_NK_Lookup.lData[lookup_offset + i];
    }
  }
  kTuple.Flip();
}
