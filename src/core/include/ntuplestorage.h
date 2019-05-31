/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon    (apoon42@uwo.ca)
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

#ifndef     __NTUPLESTORAGE__
#define     __NTUPLESTORAGE__

#include   "matrix.h"
/*__________________________________________________________________________________________________________________________________________ */

class   _NTupleStorage: public _Matrix
// a way to associate floating points
// with ordered K-tuples indices of the form 0<=i_1 < i_2 < ... < i_K < N, where N>=K
// see http://en.wikipedia.org/wiki/Combinadic
{

public:
    _NTupleStorage  (void) {
        storageN = 0;
        storageK = 0;
    }
    // null constructor; does next to nothing
    _NTupleStorage  (unsigned long,unsigned long);
    // [long] - what is the maximum value for an index (N>=0; N=0 - empty set)
    // [long] - how many indices in an ordered K-tuple (K), must be <=N; if K>N, then it will be reset to Min(1,N)
    // will allocate an array of the appropriate dimension and zero it
    // _NTupleStorage will set N=K=0 if the memory required to store the array exceeds MAX_LONG bytes (a good proxy for
    // the maximum addressable space in the system)

    virtual     ~_NTupleStorage (void) {};

    virtual     BaseRef     makeDynamic (void) const;
    // create a copy of the object dynamic copy

    // getters for N and K
    unsigned long get_N      (void) const {
        return storageN;
    }
    unsigned long get_K      (void) const {
        return storageK;
    }

    bool            CheckKTuple     (_SimpleList&);
    // modify (in-place) the argument to make it an ordered K-tuple; will also check dimensions and such
    // returns true if the k-tuple if valid; false otherwise

    unsigned long   Index           (_SimpleList const&);
    // return an index into the linear array pointed to by the K-tuple in the argument (ASSUMED to be valid here!)

    hyFloat      DirectIndex     (unsigned long);
    // retrieve a value using a direct index for the K-tuple (computed by Index)

    unsigned long   Store           (hyFloat, _SimpleList&);
    // associate a value with the K-tuple; returns the direct index for the value

    hyFloat      Retrieve        (_SimpleList const&);
    // return the value associated with the K-tuple

    void            IndexToTuple    (unsigned long, _SimpleList&);
    // given a direct index in 0..(NcK-1); return a K-tuple corresponding to that index


private:

    unsigned long   storageN,
             storageK;

    _SimpleList     C_NK_Lookup; // an (K+1)x(N+1) linear array which stores I choose J in element (J,I)

};


#endif
