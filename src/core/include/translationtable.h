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

#ifndef _HYTRANSLATIONTABLE_
#define _HYTRANSLATIONTABLE_
//#pragma once

#include "list.h"

#define  HY_TRANSLATION_TABLE_NONSTANDARD         0x000
#define  HY_TRANSLATION_TABLE_STANDARD_BINARY     0x001
#define  HY_TRANSLATION_TABLE_STANDARD_NUCLEOTIDE 0x002
#define  HY_TRANSLATION_TABLE_STANDARD_PROTEIN    0x004
#define  HY_TRANSLATION_TABLE_ANY_STANDARD        (HY_TRANSLATION_TABLE_STANDARD_BINARY|HY_TRANSLATION_TABLE_STANDARD_NUCLEOTIDE|HY_TRANSLATION_TABLE_STANDARD_PROTEIN)

//_________________________________________________________
class _TranslationTable:public BaseObj {
    


public:

    _TranslationTable                       (void);
    _TranslationTable                       (const _String&);
    /* 20100618: SLKP

            - new constructor (needed to handle ExecuteCase52 / Simulate properly)
              which takes an alphabet string and checks to see if it's a standard one
              DNA/RNA/Protein or Binary

     */
    _TranslationTable                       (_TranslationTable&);
    virtual ~_TranslationTable              (void) {
        if (checkTable) {
            free (checkTable);
        }
    }
    virtual BaseRef  makeDynamic            (void);
    virtual void     Duplicate              (BaseRef);

    const unsigned long    TokenCode        (const char) const;
    char    CodeToLetter                    (long*) const;

    void    AddBaseSet                      (const _String&);
    bool    TokenCode                       (const char, long*, const bool = true) const;
    void    SplitTokenCode                  (long, long*) const;

    void    AddTokenCode                    (const char, _String&) ;
    void    PrepareForChecks                (void);
    const bool    IsCharLegal               (const char) ;
    const char    GetSkipChar               (void) const;
    const char    GetGapChar                (void) const;
    _String ConvertCodeToLetters            (long, const char);
    const unsigned long    LengthOfAlphabet (void) const;
    inline const unsigned long Dimension    (void) const {return baseLength;}
    const   _String *           RetrieveCharacters
                                            (void) const;
    const   _String&            RetrieveAddedTokens (void) const {return tokensAdded;}
    
    
    void    Clear                           (void);

    void    SetStandardType                 (const unsigned char);
    bool    CheckType                       (const unsigned char) const;
    
    const unsigned char     DetectType      (void) const;
    
    _TranslationTable*
    MergeTables                     (_TranslationTable*);
    
    static const _String * GetDefaultAlphabet (const long);

private:
    
    static bool   CheckValidAlphabet        (const _String &);
    
    unsigned long                           baseLength;
    // number of "fundamental" tokens
    //(4 for nucl, ACGT; 20 for amino acids)


    _String                                 tokensAdded,
                                            baseSet;

    _SimpleList                             translationsAdded;
    char*                                   checkTable;
    // if null - then assume default translation table;
};

extern          _TranslationTable       defaultTranslationTable;
extern          _String                 aminoAcidOneCharCodes,
                                        dnaOneCharCodes,
                                        rnaOneCharCodes,
                                        binaryOneCharCodes;

#endif
