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
