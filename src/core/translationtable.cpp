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

#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "hy_strings.h"
#include "list.h"
#include "translationtable.h"
#include "site.h"


_TranslationTable      defaultTranslationTable;

_String           aminoAcidOneCharCodes    ("ACDEFGHIKLMNPQRSTVWY"),
                  dnaOneCharCodes          ("ACGT"),
                  rnaOneCharCodes          ("ACGU"),
                  binaryOneCharCodes       ("01");



//_________________________________________________________
_TranslationTable::_TranslationTable (void) {
    baseLength = 4;
    checkTable = NULL;
}

//_________________________________________________________
_TranslationTable::_TranslationTable (const char baseL) {
    baseLength = (baseL==20)?20:4;
    checkTable = NULL;
}

//_________________________________________________________
_TranslationTable::_TranslationTable (_TranslationTable& t) {
    Duplicate(&t);
}

//_________________________________________________________
_TranslationTable::_TranslationTable (_String& alphabet) {
    baseLength = alphabet.sLength;
    checkTable = NULL;
    if (!(alphabet.Equal (&dnaOneCharCodes) || alphabet.Equal (&rnaOneCharCodes) ||
            alphabet.Equal (&binaryOneCharCodes) || alphabet.Equal (&aminoAcidOneCharCodes))) {
        AddBaseSet (alphabet);
    }
}

//_________________________________________________________
void     _TranslationTable::Duplicate (BaseRef obj) {
    _TranslationTable * tt = (_TranslationTable *) obj;
    baseLength             = tt->baseLength;
    checkTable             = nil;
    tokensAdded            = tt->tokensAdded;
    baseSet                = tt->baseSet;
    translationsAdded.Duplicate (&tt->translationsAdded);
}

//_________________________________________________________
BaseRef     _TranslationTable::makeDynamic (void) {
    _TranslationTable * r = new _TranslationTable;
    r->Duplicate(this);
    return r;
}

//_________________________________________________________
const unsigned long    _TranslationTable::TokenCode (const char token) const {
    // standard translations
    long receptacle  [HY_WIDTH_OF_LONG];
    TokenCode               (token,receptacle);
    return bitStringToLong(receptacle, baseLength);
}

//_________________________________________________________
char    _TranslationTable::CodeToLetter (long* split) const {
// assumes a non-unique translation of split
// for unique - use ConvertCodeToLetters
    const unsigned long trsl = bitStringToLong (split, LengthOfAlphabet ());


    if (baseSet.sLength == 0) {
        // one of the standard alphabers
        if (baseLength==4) {
            // nucleotides
            switch (trsl) {
            case 3:
                return 'M';
            case 5:
                return 'S';
            case 6:
                return 'R';
            case 7:
                return 'V';
            case 9:
                return 'W';
            case 10:
                return 'Y';
            case 11:
                return 'H';
            case 12:
                return 'K';
            case 14:
                return 'B';
            }
        } else if (baseLength==20) {
            // amino acids
            switch (trsl) {
            case 2052:
                return 'B';
            case 8200:
                return 'Z';
            }
        }
    } else if (tokensAdded.sLength) {
        long f = translationsAdded.Find(trsl);
        // linear search for (binary) translations
        if (f>=0) {
            return tokensAdded.sData[f];
        }
    }
    return '?';
}

//_________________________________________________________
void    _TranslationTable::SplitTokenCode (const long code, long* receptacle) const {
    longToBitString (receptacle, code, baseLength);
}

//_________________________________________________________
const unsigned long    _TranslationTable::LengthOfAlphabet (void) const {
    return baseSet.sLength?baseSet.sLength:baseLength;
}

//_________________________________________________________

bool    _TranslationTable::TokenCode (const char token, long* receptacle, const bool gapToOnes) const {

    long f = tokensAdded.sLength?tokensAdded.Find (token):HY_NOT_FOUND;
    // check for custom translations
    // OPTIMIZE FLAG linear search:
    // SLKP 20071002 should really be a 256 char lookup table

    if (f != HY_NOT_FOUND) {
        SplitTokenCode(translationsAdded.lData[f], receptacle);
        return true;
    }

    if (baseSet.sLength) {
        // custom base alphabet
 
        memset(receptacle, 0, baseLength*sizeof(long));
        f = baseSet.Find(token);
        // OPTIMIZE FLAG linear search:
        // SLKP 20071002 should really be a 256 char lookup table

        if (f!= HY_NOT_FOUND) {
            receptacle[f] = 1;
        }

        return true;
    }

    if (baseLength==4) {
        // standard nucleotide
        memset(receptacle, 0, 4L*sizeof(long));

        switch (token) {
        case 'A':
            receptacle[0]=1;
            break;

        case 'C':
            receptacle[1]=1;
            break;

        case 'G':
            receptacle[2]=1;
            break;

        case 'T':
        case 'U':
            receptacle[3]=1;
            break;

        case 'Y':
            receptacle[3]=1;
            receptacle[1]=1;
            break;

        case 'R':
            receptacle[0]=1;
            receptacle[2]=1;
            break;

        case 'W':
            receptacle[3]=1;
            receptacle[0]=1;
            break;

        case 'S':
            receptacle[1]=1;
            receptacle[2]=1;
            break;

        case 'K':
            receptacle[3]=1;
            receptacle[2]=1;
            break;

        case 'M':
            receptacle[1]=1;
            receptacle[0]=1;
            break;

        case 'B':
            receptacle[1]=1;
            receptacle[2]=1;
            receptacle[3]=1;
            break;

        case 'D':
            receptacle[0]=1;
            receptacle[2]=1;
            receptacle[3]=1;
            break;

        case 'H':
            receptacle[1]=1;
            receptacle[0]=1;
            receptacle[3]=1;
            break;

        case 'V':
            receptacle[1]=1;
            receptacle[2]=1;
            receptacle[0]=1;
            break;

        case 'X':
        case 'N':
        case '?':
        case '.':
        case '*':
            receptacle[1]=1;
            receptacle[2]=1;
            receptacle[3]=1;
            receptacle[0]=1;
            break;

        case '-':
            if (gapToOnes) {
                receptacle[1]=1;
                receptacle[2]=1;
                receptacle[3]=1;
                receptacle[0]=1;
                break;
            }
        }
    } else {
        if (baseLength==20) {
            memset(receptacle, 0, 20L*sizeof(long));

            switch (token) {
            case 'A':
                receptacle[0]=1;
                break;

            case 'B':
                receptacle[2]=1;
                receptacle[11]=1;
                break;

            case 'C':
                receptacle[1]=1;
                break;

            case 'D':
                receptacle[2]=1;
                break;

            case 'E':
                receptacle[3]=1;
                break;

            case 'F':
                receptacle[4]=1;
                break;

            case 'G':
                receptacle[5]=1;
                break;

            case 'H':
                receptacle[6]=1;
                break;

            case 'I':
                receptacle[7]=1;
                break;

            case 'K':
                receptacle[8]=1;
                break;

            case 'L':
                receptacle[9]=1;
                break;

            case 'M':
                receptacle[10]=1;
                break;

            case 'N':
                receptacle[11]=1;
                break;

            case 'P':
                receptacle[12]=1;
                break;

            case 'Q':
                receptacle[13]=1;
                break;

            case 'R':
                receptacle[14]=1;
                break;

            case 'S':
                receptacle[15]=1;
                break;

            case 'T':
                receptacle[16]=1;
                break;

            case 'V':
                receptacle[17]=1;
                break;

            case 'W':
                receptacle[18]=1;
                break;

            case 'Y':
                receptacle[19]=1;
                break;

            case 'Z':
                receptacle[3]=1;
                receptacle[13]=1;
                break;

            case 'X':
            case '?':
            case '.':
            case '*': {
                for (int j = 0; j<20; j++) {
                    receptacle[j] = 1;
                }
            }
            break;
            case '-': {
                if (gapToOnes)
                    for (int j = 0; j<20; j++) {
                        receptacle[j] = 1;
                    }
            }
            break;
            }
        } else {
            // binary
            receptacle[0] = 0;
            receptacle[1] = 0;
            switch (token) {
            case '0':
                receptacle[0]=1;
                break;

            case '1':
                receptacle[1]=1;
                break;

            case 'X':
            case '?':
            case '.':
            case '*': {
                receptacle[0] = 1;
                receptacle[1] = 1;
            }
            break;
            case '-': {
                if (gapToOnes) {
                    receptacle[0] = 1;
                    receptacle[1] = 1;
                }
            }
            break;
            }

        }
    }
    return false;

}
//_________________________________________________________
void    _TranslationTable::PrepareForChecks (void) {
    if (checkTable == NULL) {
        checkTable = MemAllocate (256);
    }
    
    memset (checkTable, 0, 256);

    _String checkSymbols;

    if (baseSet.sLength) {
        checkSymbols = baseSet&tokensAdded;
    } else if (baseLength == 2) {
        checkSymbols = _String("01*?-.")&tokensAdded;
    } else {
        checkSymbols = _String("ABCDEFGHIJKLMNOPQRSTUVWXYZ*?-.")&tokensAdded;
    }

    for (unsigned long i=0; i<checkSymbols.sLength; i++) {
        checkTable[checkSymbols.getChar(i)] = 1;
    }
}

//_________________________________________________________
const bool    _TranslationTable::IsCharLegal (const char c)
{
    if (!checkTable) {
        PrepareForChecks();
    }
    return checkTable[c];
}
//___________________________________________

void    _TranslationTable::AddTokenCode (const char token, _String& code)
{
    long    f,
            newCode = 0L;

    bool    reset_baseset = false;

    if (baseSet.sLength==0) {
        // fill in baseSet for standard alphabets
        if (baseLength == 4) {
            baseSet = dnaOneCharCodes;
        } else if (baseLength == 20) {
            baseSet = aminoAcidOneCharCodes;
        } else {
            baseSet = binaryOneCharCodes;
        }
        reset_baseset = true;
    }


    if (baseSet.sLength) {
        long shifter = 1;
        for (unsigned long j = 0; j<baseSet.sLength; j++, shifter <<= 1)
            if (code.Find (baseSet.sData[j]) != HY_NOT_FOUND) {
                newCode += shifter;
            }
    }

    f = baseSet.Find (token);
    if (reset_baseset) {
        baseSet = empty;
    }
    
    if (f != HY_NOT_FOUND) {
        return;
    }
    // see if the character being added is a base
    // character; those cannot be redefined

    f = tokensAdded.Find (token,0,-1);
    // new definition or redefinition?

    if (f == HY_NOT_FOUND) { // new
        tokensAdded             = tokensAdded&token;
        translationsAdded       << 0L;
        f                       = tokensAdded.sLength-1L;
    }

    translationsAdded.lData[f] = newCode;
}

//_________________________________________________________

void    _TranslationTable::AddBaseSet (_String& code) {
    baseSet         = code;
    baseSet.StripQuotes();
    baseLength      = baseSet.sLength;
    if (baseLength > HY_WIDTH_OF_LONG) {
        // longer than the bit size of 'long'
        // can't handle those
        WarnError (_String ("Alphabets with more than ")
                   & HY_WIDTH_OF_LONG &
                   " characters are not supported");
    }

}

//_________________________________________________________

const char    _TranslationTable::GetSkipChar (void) const {
    if ( DetectType() != HY_TRANSLATION_TABLE_ANY_STANDARD ) {
        return '?';    // this is the default
    }

    // see if there is a symbol
    // which maps to all '1'

    long    all     = 0,
            shifter = 1;
    
    unsigned long   ul       = LengthOfAlphabet();

    for  (unsigned long f=0; f<ul; f++, shifter <<= 1) {
        all |= shifter;
    }

    if  ((all = translationsAdded.Find(all))==HY_NOT_FOUND) {
        return '?';
    } else {
        return tokensAdded.getChar(all);
    }
}

//_________________________________________________________

const char    _TranslationTable::GetGapChar (void) const {
    if ( DetectType() != HY_TRANSLATION_TABLE_ANY_STANDARD ) {
        return '-';    // default gap character
    }

    long f = translationsAdded.Find(0L);

    if  (f== HY_NOT_FOUND) {
        return 0;
    } else {
        return tokensAdded.getChar(f);
    }
}

//_________________________________________________________
_String _TranslationTable::ConvertCodeToLetters (long code, const char base) {

    _String res (base,false);
    
    const unsigned long ul = LengthOfAlphabet();
    
    if (code >= 0) {
        // OPTIMIZE FLAG; repeated memory allocation/deallocation
        if (baseSet.sLength) {
            for (long k=1; k<=base; k++, code/=ul) {
                res.sData[base-k]=baseSet.sData[code%ul];
            }
        } else {
            const _String * std_alphabet = _TranslationTable::GetDefaultAlphabet (baseLength);
            if (std_alphabet) {
                for (long k=1; k<=base; k++, code/=ul) {
                    res.sData[base-k] = std_alphabet->getChar (code%ul);
                }
            } else {
                WarnError ("Internal error in _TranslationTable::ConvertCodeToLetters; unsupported standard alphabet");
            }
        }
            
    } else {
        char c = GetGapChar();
        for (long k=0; k<base; k++) {
            res.sData[k] = c;
        }
    }
    return res;
}

//_________________________________________________________
void _TranslationTable::Clear (void) {
    baseLength = 4;
    baseSet = empty;
    tokensAdded = empty;
    translationsAdded.Clear();
    if (checkTable) {
        free (checkTable);
        checkTable = nil;
    }
}


//_________________________________________________________
void _TranslationTable::SetStandardType (unsigned const char type) {
    Clear();
    switch (type) {
        case HY_TRANSLATION_TABLE_STANDARD_BINARY:
            baseLength = 2;
            return;
        case HY_TRANSLATION_TABLE_STANDARD_PROTEIN:
            baseLength = 20;
            return;
    }
    baseLength = 4;
}

//_________________________________________________________
bool _TranslationTable::CheckType (unsigned char pattern) const {
    return DetectType () & pattern;
}
//_________________________________________________________
const unsigned char _TranslationTable::DetectType (void) const {
    if (baseSet.sLength == 0UL && translationsAdded.lLength == 0UL && tokensAdded.sLength == 0UL) {
        switch (baseLength) {
            case 2:
                return HY_TRANSLATION_TABLE_STANDARD_BINARY;
            case 4:
                return HY_TRANSLATION_TABLE_STANDARD_NUCLEOTIDE;
            case 20:
                return HY_TRANSLATION_TABLE_STANDARD_PROTEIN;
        }
    }
    return HY_TRANSLATION_TABLE_NONSTANDARD;
}

//_________________________________________________________

_TranslationTable*  _TranslationTable::MergeTables (_TranslationTable* table2)
// merge the translation tables if they are compatible, return the result,
// otherwise return nil
{
    const unsigned char my_type = DetectType();
    
    if (my_type != table2->DetectType()) {
        return nil;
    }
    
    if ( my_type == HY_TRANSLATION_TABLE_NONSTANDARD && !baseSet.Equal (&table2->baseSet)) {
        return nil;
    }
    
    _TranslationTable* result = new _TranslationTable (*this);

    if (table2->tokensAdded.sLength) {
        for (unsigned long i=0; i<table2->tokensAdded.sLength; i++) {
            long f = tokensAdded.Find (table2->tokensAdded[i]);
            if (f==HY_NOT_FOUND) {
                result->tokensAdded       && table2->tokensAdded[i];
                result->translationsAdded << table2->translationsAdded[i];
            } else if (translationsAdded.lData[f] != table2->translationsAdded.lData[i]) {
                DeleteObject (result);
                return nil;
            }
        }
        return result;
    } else {
        return result;
    }

    return nil;
}

//_________________________________________________________
const _String* _TranslationTable::RetrieveCharacters (void) const {
    if (baseSet.sLength) {
        return &baseSet;
    }
    
    const _String* res = GetDefaultAlphabet (Dimension());
    return res ? res : &empty;
}

//_________________________________________________________
const _String* _TranslationTable::GetDefaultAlphabet (const long size) {
    switch (size) {
        case 2:
            return &binaryOneCharCodes;
        case 4:
            return &dnaOneCharCodes;
        case 20:
            return &aminoAcidOneCharCodes;
    }
    return nil;
}


