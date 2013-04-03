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

#ifndef     __HELPERS__
#define     __HELPERS__

#include "baseobj.h"

void            DeleteObject (BaseRef); // delete a dynamic object

char*           MemAllocate (const long);
char*           MemReallocate (Ptr, const long);

bool            GlobalStartup();
bool            GlobalShutdown();

unsigned long   bitStringToLong (const long *, const unsigned long );
void            longToBitString      (long *, const unsigned long, const unsigned long);


void            PurgeAll                  (bool   all = true);
void            init_genrand              (unsigned long);
unsigned long   genrand_int32             (void);
double          genrand_real2             (void);
FILE*           doFileOpen                (const char *, const char *, bool = false);
// 20110324: SLKP added the bool flag to allow automatic "Can't open file" error reports
double          TimerDifferenceFunction   (bool);



#if !defined __UNIX__ || defined __HEADLESS__ || defined __HYPHY_GTK__ || defined __HYPHYQT__
    void    yieldCPUTime        (void);
    bool    handleGUI           (bool = false);
#endif

template <class multipliable>
multipliable compute_power (multipliable base, unsigned long power) {
    multipliable result = 1.;
    unsigned long bit = 0x7FFFFFFF;
    
    while ((bit & power) == 0UL && bit > 0) {
        bit = bit >> 1;
    }
    
    while (power > 0) {
        result *= result;
        if (bit & power) {
            result *= base;
        }
        power = power >> 1;
        bit = bit >> 1;
    }
    
    return result;
}

#endif
