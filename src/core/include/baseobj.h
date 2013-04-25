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

#ifndef _HBASEOBJ_
#define _HBASEOBJ_
//#pragma once

typedef char *Ptr;
typedef double _Parameter;

#include "defines.h"

#ifdef __HEADLESS__
#include "THyPhy.h"
#endif

#ifdef __GNUC__
#define _hprestrict_ __restrict
#else
#define _hprestrict_
#endif

#include "stdio.h"

class BaseObj {

  //base object class
public:

  BaseObj();

  virtual ~BaseObj(void) {}

  virtual BaseObj *toStr(void);

  virtual BaseObj *toErrStr(void);

  virtual void toFileStr(FILE *);

  virtual BaseObj *makeDynamic(void);

  virtual long FreeUpMemory(long) { return 0L; }

  virtual void Initialize(void) { nInstances = 1; }

  virtual void Duplicate(BaseObj *ref) { nInstances = ++ref->nInstances; }

  virtual void AddAReference(void) { nInstances++; }

  virtual void RemoveAReference(void) { nInstances--; }

  long nInstances;

};

typedef BaseObj *BaseRef;

#ifdef _SLKP_USE_SSE_INTRINSICS
#include <pmmintrin.h>
#endif

#ifdef _SLKP_USE_AVX_INTRINSICS
#include <immintrin.h>
#endif

#endif

//EOF
