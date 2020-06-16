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

#ifndef _HBASEOBJ_
#define _HBASEOBJ_
//#pragma once



#ifdef      __HEADLESS__
#include "THyPhy.h"
#endif

#include <stdio.h>

#include "defines.h"
#include "hy_types.h"


class BaseObj {
    
    /** 
        This is a legacy implementation of reference - counted objects.
        The idea is that when an object reaches the reference count of 0,
        it can be deleted, because it is no longer being used. 
     
        Storing objects by reference or shallow cloning will increase reference counts,
        inverse operations will decrease said counts.
     
        Also defines mandatory overload functions
     
     */
private:
    long reference_counter;
    
public:
    
    BaseObj();
    
    virtual ~BaseObj(void) {}
    
    virtual BaseObj *toStr (unsigned long padding = 0UL);
        /** string representation for the object 
            @param padding is used to allow 'pretty' rendering of nested objects, like dictss
         */
    
    virtual BaseObj *toErrStr () ;
        /** error string representation for the object 
            (default = same as toStr) 
         */

    virtual void toFileStr (FILE *, unsigned long padding = 0UL) ;
        /** file representation for the object 
            @param padding is used to allow 'pretty' rendering of nested objects, like dictss
         */
    
    virtual BaseObj *makeDynamic(void) const = 0;
    
    virtual void Initialize(bool = true) { reference_counter = 1L; }
    
    virtual void Duplicate(BaseObj const * ref) = 0;
    
    inline void AddAReference(void) { reference_counter++; }
    
    inline void RemoveAReference(void) { reference_counter--; }
    
    inline bool CanFreeMe (void)  const { return reference_counter <= 1L; }
    
    inline bool SingleReference (void)  const { return reference_counter == 1L; }
    // comparison functions
    
    
};

typedef BaseObj*        BaseRef;
typedef BaseObj const * BaseRefConst;


bool    DeleteObject (BaseRef object);
/** 
    Remove one reference count from an object, and if it is no
    longer pointed to by anything (counter = 0), then delete it.
 
    @param object the object to "delete"
    
    @retrun true if object was deleted, otherwise false
 
 */



#endif

//EOF
