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

#ifndef     __HY_TYPES__
#define     __HY_TYPES__

#ifdef __ZLIB__
    #include <zlib.h>
#endif

#include <stdio.h>

/**
    Put application-wide typedefs and enums here
*/

typedef     void*        hyPointer;
    // generic pointer
typedef     double       hyFloat;
    // standard floating type



enum hyBLFunctionType  {
  kBLFunctionAlwaysUpdate,
  kBLFunctionSkipUpdate,
  kBLFunctionLocal
};

enum hyBLFunctionArgumentType {
  kBLFunctionArgumentNormal,
  kBLFunctionArgumentReference
};

enum hyTreeDefinitionPhase {
  kTreeNotBeingDefined,
  kTreeIsBeingParsed,
  kTreeNodeBeingCreated
};

enum hyComparisonType {
    kCompareLess    = -1,
    kCompareEqual   = 0,
    kCompareGreater = 1,
    kCompareUndefined      = 0xff
};

enum hyFileOpenMode {
    kFileRead,
    kFileReadBinary,
    kFileReadWrite,
    kFileWrite,
    kFileWriteBinary,
    kFileAppend
};

class hyFile {
    public:
        hyFile (void) {
            _fileReference = NULL;
            #ifdef __ZLIB__
                _fileReferenceDirect = NULL;
            #endif
        }
        static hyFile* openFile (const char * file_path, hyFileOpenMode mode , bool error = false, bool compress = false, long buffer = 1024*128);
        void lock (void);
        void unlock (void);
        void rewind (void);
        void seek (long, int);
        void flush (void);
        int  close ();
        bool feof (void);
        unsigned long read (void* buffer, unsigned long size, unsigned long items);
        size_t fwrite( const void* buffer, size_t size, size_t count);
        int    puts(const char *str);
        int    fputc(int chr);
        size_t tell ();
        int getc ();
    #ifdef __ZLIB__
        inline  bool valid (void) const {return _fileReference != NULL || _fileReferenceDirect != NULL;}
        gzFile _fileReference;
        FILE * _fileReferenceDirect;
        bool    is_compressed (void) const { return _fileReference != NULL;}
    #else
        inline  bool valid (void) const {return _fileReference != NULL;}
        bool    is_compressed (void) const { return false;}
        FILE* _fileReference;
    #endif
};

#endif
