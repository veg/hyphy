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

#include "global_things.h"
#include "function_templates.h"
#include "vector.h"

#define  VECTOR_DEFAULT_INITIAL_ALLOCATION  64L


using namespace hy_global;


_Vector::_Vector (bool iscol) : _Matrix (VECTOR_DEFAULT_INITIAL_ALLOCATION,1,false,true) {
    used = 0UL;
    is_column = iscol;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void _Vector::Trim (void) {
  Resize (used);
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

BaseRef     _Vector::makeDynamic (void) const {
    _Vector * result = new _Vector;
    result->_Matrix::Duplicate (this);
    result->used = used;
    result->is_column = is_column;
    return result;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

long        _Vector::Store (hyFloat value) {
    if (used < hDim) {
        theData[used++] = value;  // increment AFTER argument is sent to function
        return used-1UL;
    } else {
        Resize (used + MAX (used>>3,VECTOR_DEFAULT_INITIAL_ALLOCATION));    // allocate another block of 64
        return Store (value);
    }
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void        _Vector::Delete (unsigned long index) {
  if (index + 1UL < used) {
    for (unsigned long i = index+1UL; i < used; i++) {
      theData[i-1] = theData[i];
    }
    --used;
  }
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void        _Vector::DeleteList (_SimpleList const& toDelete) {
    
    unsigned long delete_count = toDelete.countitems();
    
    if (delete_count) {
        unsigned long k = 0;
        for (unsigned long i = 0UL; i<used; i++) {
            if (k<delete_count && i==toDelete.get(k)) {
                k++;
            } else {
                theData[i-k] = theData[i];
            }
        }
        used -= delete_count;
    }
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void        _Vector::operator << (const _SimpleList& list) {
    list.Each([this] (long value, long) -> void {
        this->Store(value);
    });
    
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void        _Vector::Clear (bool complete_clear) {
    _Matrix::Clear(complete_clear);
    ZeroUsed();
    vDim = 1UL;
}

//_____________________________________________________________________________________________
void _Vector::Duplicate (BaseRefConst obj) {
    _Matrix::Duplicate (obj);
    used = ((_Vector*)obj)->used;
    is_column = ((_Vector*)obj)->is_column;
}


