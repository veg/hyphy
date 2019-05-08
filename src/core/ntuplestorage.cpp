/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon42@uwo.ca)
 Steven Weaver (sweaver@temple.edu)
 
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

#include "ntuplestorage.h"

/*--------------------------------------------------------------------------------------------------------------------------------*/
// _NTupleStorage
/*--------------------------------------------------------------------------------------------------------------------------------*/

_NTupleStorage::_NTupleStorage (unsigned long N, unsigned long K) {
    storageN        = N;
    storageK        = (K>N)?MIN(N,1):K;
    unsigned long   matrixDimension = 1UL;

    // now compute what dimension of the matrix will be required
    // handle the special cases first

    if (storageK) { // not just the kEmptyString set
        // populate C_NK_Lookup
        C_NK_Lookup.Populate(storageN+1L, 1L, 0);
        // N choose 0 is always 1 for every N
        
        
        for (long i=1L; i<=storageK; i++) {
            for (long filler = 0; filler < i; filler++, C_NK_Lookup<<0); // N choose K where K>N is invalid
            C_NK_Lookup << 1L;
            // K choose K is 1
            for (long j=i+1L; j<=storageN; j++) {
                // N choose K = N/(N-K) times (N-1) choose K
                C_NK_Lookup << C_NK_Lookup.get (C_NK_Lookup.countitems()-1L) * j/(j-i);
            }
        }
    }


    matrixDimension = C_NK_Lookup.get (C_NK_Lookup.countitems()-1);
    CreateMatrix (this, 1, matrixDimension, false, true);
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

BaseRef _NTupleStorage::makeDynamic (void) const {
    _NTupleStorage* copy = new _NTupleStorage;
    copy->_Matrix::Duplicate (this);
    copy->storageN = storageN;
    copy->storageK = storageK;
    copy->C_NK_Lookup.Duplicate (&C_NK_Lookup);
    return copy;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

bool    _NTupleStorage::CheckKTuple (_SimpleList& kTuple) {
    if (kTuple.countitems() == storageK) {
        if (storageK) {
            kTuple.Sort();
            for (long k=0L; k< storageK; k++) {
                long const item = kTuple.get(k);
                if (item < 0L || item >= storageN || (k && item == kTuple.get(k-1L))) {
                    return false;
                }
            }
        }
        return true;
    }
    return false;
}


/*--------------------------------------------------------------------------------------------------------------------------------*/

unsigned long   _NTupleStorage::Index (_SimpleList const& kTuple) {
    unsigned long myIndex = 0UL;
    if (storageK)
        for (long k=kTuple.lLength-1L; k >=0; k--) {
            myIndex += C_NK_Lookup.get((k+1L)*(1L+storageN) + kTuple.get(k));
        }
    return myIndex;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

hyFloat  _NTupleStorage::DirectIndex (unsigned long directIndex) {
    return theData[directIndex];
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

unsigned long   _NTupleStorage::Store (hyFloat value, _SimpleList& kTuple) {
    unsigned long myIndex = Index (kTuple);
    theData[myIndex] = value;
    return myIndex;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

hyFloat  _NTupleStorage::Retrieve (_SimpleList const& kTuple) {
    return theData[Index (kTuple)];
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void    _NTupleStorage::IndexToTuple (unsigned long directIndex, _SimpleList& kTuple) {
    kTuple.Clear();
    if (storageK && directIndex < C_NK_Lookup.get(C_NK_Lookup.countitems()-1)) {
        long currentN = storageN-1L;
        for (long k=storageK; k > 0L; k--) {
            long  i             = currentN,
                  lookup_offset = k*(storageN+1);

            while (C_NK_Lookup.list_data[lookup_offset+i] > directIndex) {
                //printf ("(%d, %d) -> %d\n", k, i, C_NK_Lookup.list_data[lookup_offset+i]);
                i--;
            }
            //printf ("(%d, %d) -> %d\n", k, i, C_NK_Lookup.list_data[lookup_offset+i]);
            kTuple << i;
            //printf ("Stored %d\n", i);

            currentN     = i-1L;
            directIndex -= C_NK_Lookup.list_data[lookup_offset+i];
        }
    }
    kTuple.Flip();
}

