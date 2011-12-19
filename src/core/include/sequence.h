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

#ifndef _HXSTRINGS_
#define _HXSTRINGS_
//#pragma once
#include "hy_strings.h"

#define    NOCOMPRESSION   0
#define    LZWCOMPRESSION  1
#define    FREQCOMPRESSION 2

#define    NUCLEOTIDEALPHABET 128
#define    CODONALPHABET      64
#define    FULLALPHABET       32
#define    FULLNUCLALPHABET   16



//_________________________________________________________

class _CString:public _String   // compressible string
{

public:
    _CString (void);
    //does nothing
    _CString (_String&);
    // string contructor
    _CString (char*);
    // data constructor
    _CString (unsigned long l, bool flag);
    // data constructor
    _CString (char);
    // data constructor
    virtual     ~_CString(void);
    //destructor

    virtual void operator << (_String*);
    // append into operator

    virtual void operator << (char);
    // append into operator

    virtual void Finalize (void);

    virtual     BaseRef makeDynamic (void);
    // create a dynamic copy of this object

    virtual     long    FreeUpMemory (long);

    void        SetFlag (unsigned char flag) {
        compressionType|=flag;
    }
    bool        IsFlag  (unsigned char flag) {
        return (bool)(compressionType&flag);
    }
    static      _String*    SelectAlpha (unsigned char);

    _Parameter      LZWCompress (unsigned char theAlpha); // returns compression ratio
    _Parameter      FrequencyCompress(unsigned char theAlpha, bool  doit = true);
    _Parameter      BestCompress(unsigned char theAlpha, long triggerSize = 25);

    _String*        Decompress (void);

    bool        IsCompressed (void) {
        return (IsFlag (LZWCOMPRESSION)||IsFlag(FREQCOMPRESSION));
    }

    void        SetDecompressed (void) {
        compressionType&=0xf0;
    }

    virtual     void    Duplicate (BaseRef ref);


    _String*        DecompressFrequency (void);
    _String*        DecompressLZW (void);

    // compression flag

    unsigned long            allocatedSpace;
private:

    unsigned char    compressionType;
};

//_________________________________________________________


void        SetAlphabet (_String&);

#endif
