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

#ifndef     __VECTOR__
#define     __VECTOR__

#include "hy_strings.h"
#include "matrix.h"

/*__________________________________________________________________________________________________________________________________________ */

class   _Vector: public _Matrix
// automatically growing matrix class
{

public:
    _Vector  (bool is_col = true);
    virtual     ~_Vector (void) {};

    virtual     BaseRef     makeDynamic (void) const; // duplicate this object into a dynamic copy
    virtual     void        Duplicate   (BaseRefConst); // duplicate an object from reference

    virtual     void        Clear (bool complete = true);

    virtual     unsigned long        GetHDim                     (void) const {
        if (is_column) {
            return get_used();
        }
        return 1UL;
    }
    virtual     long        GetVDim                     (void) const {
        if (!is_column) {
            return get_used();
        }
        return 1UL;
    }
  
    void   Trim             (void);
  
    long   Store            (hyFloat);
    long   get_used          (void) const {
        return used;
    }
    void     ZeroUsed       (void) {
        used = 0UL;
    }
    
    void        Delete (unsigned long index);
    void        DeleteList (_SimpleList const& );

    void    operator <<     (const _SimpleList&);
    _Vector&    operator <<     (hyFloat p) {
        Store (p);
        return *this;
    }
    _Vector&   AppendRange (unsigned long how_many, hyFloat start, hyFloat add) {
        for (unsigned long i = 0; i < how_many; i++) {
            Store (start);
            start += add;
        }
        return *this;
    }

private:
    unsigned long   used;
    bool   is_column;
};

#endif
