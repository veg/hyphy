/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2006
Primary Development:
  Sergei L Kosakovsky Pond (sergeilkp@mac.com)
Significant contributions from:
  Spencer V Muse (muse@stat.ncsu.edu)
  Simon DW Frost (sdfrost@ucsd.edu)
  Art FY Poon    (apoon@biomail.ucsd.edu)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#ifndef _HBASEOBJ_
#define _HBASEOBJ_
//#pragma once


typedef char* Ptr;

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE  1
#endif

#ifndef nil
#define nil   NULL
#endif

#ifdef      __HEADLESS__
#include "THyPhy.h"
#endif

#ifdef __GNUC__
#define _hprestrict_ __restrict
#else
#define _hprestrict_
#endif

#include "stdio.h"

class BaseObj
{

    //base object class
public:

    BaseObj();

    virtual ~BaseObj(void) {}

    virtual BaseObj* toStr (void);

    virtual BaseObj* toErrStr (void);

    virtual void     toFileStr (FILE*);

    /*virtual operator const char*(void);*/

    virtual BaseObj* makeDynamic (void);

    virtual long     FreeUpMemory (long) {
        return 0;
    }

    virtual void     Initialize (void) {
        nInstances=1;
    }

    virtual void     Duplicate (BaseObj* ref) {
        nInstances=++ref->nInstances;
    }

    virtual void     AddAReference (void)     {
        nInstances ++;
    }

    long             nInstances;


};

typedef BaseObj*  BaseRef;

extern  void      DeleteObject (BaseRef); // delete a dynamic object


#ifdef  __HYPHYDMALLOC__
#define MemAllocate(X)      malloc(X)
#define MemReallocate(X,Y)  realloc(X,Y)
#else
char*   MemAllocate (long);
char*   MemReallocate (Ptr, long);
#endif

bool    GlobalStartup();
bool    GlobalShutdown();


extern  FILE*   globalErrorFile;
extern  FILE*   globalMessageFile;
extern  bool    terminateExecution,
        skipWarningMessages;

extern long     systemCPUCount;

void          PurgeAll                  (bool   all = true);
void          init_genrand              (unsigned long);
unsigned long genrand_int32             (void);
double        genrand_real2             (void);
FILE*         doFileOpen                (const char *, const char *, bool = false);
// 20110324: SLKP added the bool flag to allow automatic "Can't open file" error reports
double        TimerDifferenceFunction   (bool);

#define       RAND_MAX_32               4294967295.0
#define       USE_AVL_NAMES
#define       HY_WIDTH_OF_LONG          ((long)(sizeof(long)*8))

#ifdef __HYALTIVEC__
#include    "altivec.h"
char*   VecMemAllocate (long);
typedef     float  _Parameter; // standard number type - used everywhere in matrices.
// Adjust if extra precision or memory savings are in order
#else
#ifdef      __USE_LONG_DOUBLE__
#define     PRINTF_FORMAT_STRING    "%Lg"
typedef     long double  _Parameter; // standard number type - used everywhere in matrices.
#else
#define     PRINTF_FORMAT_STRING    "%g"
typedef     double       _Parameter; // standard number type - used everywhere in matrices.
#endif
#endif

#if !defined __UNIX__ || defined __HEADLESS__
void    yieldCPUTime        (void);
bool    handleGUI           (bool = false);
#endif

#endif

//EOF
