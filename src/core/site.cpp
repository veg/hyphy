/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon@cfenet.ubc.ca)
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

#define  HYPHY_SITE_DEFAULT_BUFFER_SIZE 256
#define   DATA_SET_SWITCH_THRESHOLD     100000

#include <string.h>
#include <ctype.h>
#include <stdlib.h>


#include "site.h"
#include "likefunc.h"
#include "function_templates.h"
#include "hbl_env.h"
#include "avllistxl_iterator.h"


#include "math.h"
#include "string_file_wrapper.h"

#include "global_object_lists.h"

using namespace hyphy_global_objects;


#ifdef __MAC__
extern bool handleGUI(bool);
#endif

#if !defined  __UNIX__ && !defined __HEADLESS__
#include "HYDataPanel.h"
#endif

#if !defined  __UNIX__ || defined __HEADLESS__
#include "preferences.h"
#endif

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

#ifdef    __HYPHYMPI__
extern int _hy_mpi_node_rank;
#endif


_TranslationTable      defaultTranslationTable;

//_________________________________________________________

extern _Parameter dFPrintFormat,
       dFDefaultWidth;

//_________________________________________________________


//_________________________________________________________

/* function declarations */

void    checkTTStatus (FileState* fs);
void    processCommand (_String*s, FileState*fs);
void    FilterRawString (_String& s, FileState* fs, _DataSet & ds);
long    ProcessLine (_String&s , FileState *fs, _DataSet& ds);
void    PadLine (FileState& fState, _DataSet& result);
void    ISelector (FileState& fState, _String& CurrentLine, _DataSet& result);
bool    SkipLine (_String& theLine, FileState* fS);
void    TrimPhylipLine (_String& CurrentLine, _DataSet& ds);
void    ProcessTree    (FileState*, FILE*, _String&);
void    ReadNexusFile (FileState& fState, FILE*f, _DataSet& result);

//_________________________________________________________


_List _TranslationTable::_list_of_default_tables ( _List() < "ACGT"
                                           < "ACGU"
                                           < "ACDEFGHIKLMNPQRSTVWY"
                                           < "01");




const _String& _TranslationTable::GetDefaultTable(long tableType) {
  
  switch (tableType) {
    case HY_TRANSLATION_TABLE_BINARY:
      return *(_String*)_list_of_default_tables (3);
    case HY_TRANSLATION_TABLE_RNA:
      return *(_String*)_list_of_default_tables (1);
    case HY_TRANSLATION_TABLE_PROTEIN:
      return *(_String*)_list_of_default_tables (2);
   case HY_TRANSLATION_TABLE_DNA:
      return *(_String*)_list_of_default_tables (0);
  }
  
  return emptyString;
  
}

_TranslationTable::_TranslationTable (void)
{
    baseLength = 4;
    checkTable = NULL;
}

//_________________________________________________________
_TranslationTable::_TranslationTable (unsigned char baseL) {
    baseLength = (baseL==20)?20:4;
    checkTable = NULL;
}

//_________________________________________________________
_TranslationTable::_TranslationTable (_TranslationTable const& t) {
    tokensAdded = t.tokensAdded;
    baseLength = t.baseLength;
    baseSet = t.baseSet;
    translationsAdded << t.translationsAdded;
    checkTable = NULL;
}

//_________________________________________________________
_TranslationTable::_TranslationTable (_String& alphabet) {
    baseLength = alphabet.sLength;
    checkTable = NULL;
    if (_list_of_default_tables.FindObject (&alphabet) < 0L) {
        AddBaseSet (alphabet);
    }
}

//_________________________________________________________
BaseRef     _TranslationTable::makeDynamic (void) {
    _TranslationTable * r = new _TranslationTable;
    checkPointer(r);

    memcpy ((char*)r, (char*)this, sizeof (_TranslationTable));
    r->nInstances = 1;
    r->tokensAdded.Duplicate (&tokensAdded);
    r->baseSet.Duplicate (&baseSet);
    r->translationsAdded.Duplicate (&translationsAdded);
    r->checkTable = NULL;
    return r;
}

//_________________________________________________________
long    _TranslationTable::TokenCode (char token) const
{
    // standard translations
    long   receptacle[256],
           resolution_count = TokenResolutions(token,receptacle);

    long                    theCode         = 0L,
                            shifter         = 1L;

    for (unsigned long i = 0; i < resolution_count; i++) {
      theCode |= (1L << receptacle[i]); // set the right bit
    }
  
    return theCode;
}

//_________________________________________________________
char    _TranslationTable::AmbigToLetter (long* split, unsigned long resolutions) const
// assumes a non-unique translation of split
// for unique - use ConvertCodeToLetters
{
    long     encoding    = 0L;

    for (unsigned long k=0UL; k< resolutions ; k++) {
        encoding |= (1L << split[k]);
    }

    if (baseSet.sLength == 0)
        // one of the standard alphabers
    {
        if (baseLength==4)
            // nucleotides
        {
            switch (encoding) {
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
        } else if (baseLength==20)
            // amino acids
        {
            switch (encoding) {
              case 2052:
                  return 'B';
              case 8200:
                  return 'Z';
            }
        }
    } else if (tokensAdded.sLength) {
        long lookup = translationsAdded.Find(encoding);
        // linear search for (binary) translations
        if (lookup>=0L) {
            return tokensAdded.sData[lookup];
        }
    }
    return '?';
}

//_________________________________________________________
void    _TranslationTable::SplitTokenCode (long code, long* receptacle) const {
    unsigned long shifter = 1L;
    for (unsigned int i=0; i<baseLength; i++) {
        receptacle[i] = ((code&shifter)!=0) ? 1L : 0L;
        shifter >>= 1;
    }
}

//_________________________________________________________
long    _TranslationTable::LengthOfAlphabet (void) const {
    return baseSet.sLength?baseSet.sLength:baseLength;
}

//_________________________________________________________

const _String&   _TranslationTable::ExpandToken            (char token) const {
  long    buf [256];
  
  long resolution_count = TokenResolutions (token, buf);
  _String const * base_set = &GetAlphabetString();
  _String expansion (base_set->sLength, true);
  
  for (unsigned long tc = 0; tc < resolution_count; tc++) {
    expansion << base_set->sData[buf[tc]];
  }
  
  expansion.Finalize();
  return expansion;
}

//_________________________________________________________

long    _TranslationTable::MultiTokenResolutions (_String const& tokens, long* receptacle, bool gapToOnes) const {

  if (tokens.sLength == 1UL) {
    return TokenResolutions (tokens.getChar(0UL), receptacle, gapToOnes);
  } else {
  
    long * large_store,
           large_store_static [HYPHY_SITE_DEFAULT_BUFFER_SIZE];
    
    
    if (baseLength * tokens.sLength + tokens.sLength >= HYPHY_SITE_DEFAULT_BUFFER_SIZE) {
      large_store = new long [baseLength * tokens.sLength + tokens.sLength];
    } else {
      large_store = large_store_static;
    }
    
    /*
     large_store is a linear array which stores the following data
     
     [0,unitLength) -- the number of resolutions for the i-th character
     
     [unitLength,unitLength + baseLength] -- the actual resolutions for the 1st char
     [unitLength + baseLength, unitLength + 2*baseLength] -- the actual resolutions for the 2nd char
     ...
     */
    
    long resolution_count = 1L;
    
    for (unsigned long char_index = 0; char_index < tokens.sLength ; char_index++) {
      large_store [char_index] = TokenResolutions (tokens.sData[char_index], large_store + tokens.sLength + baseLength * char_index, gapToOnes);
      if (gapToOnes && large_store [char_index] == 0) {
        large_store [char_index] = baseLength;
        InitializeArray(large_store + tokens.sLength + baseLength * char_index, baseLength, 1L);
      }
      resolution_count *= large_store [char_index] > 0 ? large_store [char_index] : 0;
    }
    
    if (resolution_count == 1L) {
      for (unsigned long char_index = 0; char_index < tokens.sLength ; char_index++) {
        large_store[char_index] = large_store[tokens.sLength + baseLength * char_index];
      }
      
      if (receptacle) {
        receptacle [0]   = CombineDigits(large_store, tokens.sLength, baseLength);
      } else {
        resolution_count = CombineDigits(large_store, tokens.sLength, baseLength);
      }
    } else {
      if (receptacle) {
        // handle cases of 2 and 3 characters separately since they are the most common
        
        if (resolution_count > HYPHY_SITE_DEFAULT_BUFFER_SIZE) {
          FlagError(_String ("Too many ambiguous states in call to ") & _String (__PRETTY_FUNCTION__).Enquote());
          return -1L;
        }
        
        if (tokens.sLength == 3) {
          long digits[3],
              *resolution_arrays [3] = {large_store + tokens.sLength, large_store + tokens.sLength + baseLength,large_store + tokens.sLength + 2*baseLength},
              resolutions_index = 0L;
          
          for (unsigned long digit1 = 0; digit1 < large_store[0]; digit1 ++) {
            for (unsigned long digit2 = 0; digit2 < large_store[1]; digit2 ++) {
              for (unsigned long digit3 = 0; digit3 < large_store[2]; digit3 ++) {
                receptacle[resolutions_index++] = resolution_arrays[0][digit1] * baseLength * baseLength + resolution_arrays[1][digit2] * baseLength + resolution_arrays[2][digit3];
              }
            }
          }
          
        } else {
          if (tokens.sLength == 2) {
            long digits[2],
                 *resolution_arrays [2] = {large_store + tokens.sLength,large_store + tokens.sLength + baseLength},
                 resolutions_index = 0L;
            
            for (unsigned long digit1 = 0; digit1 < large_store[0]; digit1 ++) {
              for (unsigned long digit2 = 0; digit2 < large_store[1]; digit2 ++) {
                receptacle[resolutions_index++] = resolution_arrays[0][digit1] * baseLength + resolution_arrays[1][digit2];
              }
            }
          } else { // more than 3 tokens [rare!]
            
            if (tokens.sLength >= 32) {
              FlagError(_String ("The token string is too long in call to ") & _String (__PRETTY_FUNCTION__).Enquote());
              return -1L;
            }
            
            long digits[32] {},
                 resolutions_index = 0L;
            
            do {
              // assemble the current token, backwards
              long this_resolution = 0L,
                   weight = 1L;
              for (long digit = tokens.sLength - 1; digit >= 0; digit --) {
                this_resolution += weight * *(large_store + tokens.sLength + baseLength * digit + digits[digit]);
                weight *= tokens.sLength;
              }
              
              receptacle[resolutions_index++] = this_resolution;
              
              for (long digit = tokens.sLength - 1; digit >= 0; digit --) {
                if (++digits[digit] < large_store[digit]) {
                  break;
                }
                if (digit > 0) {
                  digits[digit] = 0L;
                }
              }
              
            } while (digits[0] < large_store[0]);
            
          }
        }
      } else {
        resolution_count = -1L;
      }
    }
    
    
    if (large_store != large_store_static) {
      delete [] large_store;
    }
    
    return resolution_count;
  }

}

//_________________________________________________________

long    _TranslationTable::TokenResolutions (char token, long* receptacle, bool gapToOnes) const {
  
  long custom_code = tokensAdded.sLength ? tokensAdded.Find (token):-1;
  long resolution_counter = -1L;
  
  if (custom_code != -1) {
    resolution_counter = 0L;
    unsigned long shifter = 1L;
    for (unsigned long i=0UL; i<baseLength; i++) {
      if  ((custom_code&shifter)!=0) {
        receptacle[resolution_counter++] = i;
      }
      shifter >>= 1;
    }
  } else {
    
    if (baseSet.sLength) {
      
      long base_char = baseSet.Find(token);
      // OPTIMIZE FLAG linear search:
      // SLKP 20071002 should really be a 256 char lookup table
      
      if ( base_char !=-1 ) {
        resolution_counter = 1;
        receptacle[0] = base_char;
      }
    } else {
      
      if (baseLength==4) {
        
        switch (token) {
          case 'A':
            resolution_counter = 1L;
            receptacle[0]=0;
            break;
            
          case 'C':
            resolution_counter = 1L;
            receptacle[0]=1;
            break;
            
          case 'G':
            resolution_counter = 1L;
            receptacle[0]=2;
            break;
            
          case 'T':
          case 'U':
            resolution_counter = 1L;
            receptacle[0]=3;
            break;
            
          case 'Y':
            resolution_counter = 2L;
            receptacle[0]=1;
            receptacle[1]=3;
            break;
            
          case 'R':
            resolution_counter = 2L;
            receptacle[0]=0;
            receptacle[1]=2;
            break;
            
          case 'W':
            resolution_counter = 2L;
            receptacle[0]=0;
            receptacle[1]=3;
            break;
            
          case 'S':
            resolution_counter = 2L;
            receptacle[0]=1;
            receptacle[1]=2;
            break;
            
          case 'K':
            resolution_counter = 2L;
            receptacle[0]=2;
            receptacle[1]=3;
            break;
            
          case 'M':
            resolution_counter = 2L;
            receptacle[0]=0;
            receptacle[1]=1;
            break;
            
          case 'B':
            resolution_counter = 3L;
            receptacle[0]=1;
            receptacle[1]=2;
            receptacle[2]=3;
            break;
            
          case 'D':
            resolution_counter = 3L;
            receptacle[0]=0;
            receptacle[1]=2;
            receptacle[2]=3;
            break;
            
          case 'H':
            resolution_counter = 3L;
            receptacle[0]=0;
            receptacle[1]=1;
            receptacle[2]=3;
            break;
            
          case 'V':
            resolution_counter = 3L;
            receptacle[0]=0;
            receptacle[1]=1;
            receptacle[2]=2;
           break;
            
          case 'X':
          case 'N':
          case '?':
          case '.':
          case '*':
            resolution_counter = 4L;
            receptacle[0]=0;
            receptacle[1]=1;
            receptacle[2]=2;
            receptacle[3]=3;
            break;
            
          case '-':
            resolution_counter = 0L;
            break;
        }
      } else {
        if (baseLength==20) {
          
          
          switch (token) {
            case 'A':
              resolution_counter = 1L;
              receptacle[0]=0;
              break;
              
            case 'B':
              resolution_counter = 2L;
              receptacle[0]=2;
              receptacle[1]=11;
              break;
              
            case 'C':
              resolution_counter = 1L;
              receptacle[0]=1;
              break;
              
            case 'D':
              resolution_counter = 1L;
              receptacle[0]=2;
              break;
              
            case 'E':
              resolution_counter = 1L;
              receptacle[0]=3;
              break;
              
            case 'F':
              resolution_counter = 1L;
              receptacle[0]=4;
              break;
              
            case 'G':
              resolution_counter = 1L;
              receptacle[0]=5;
              break;
              
            case 'H':
              resolution_counter = 1L;
              receptacle[0]=6;
              break;
              
            case 'I':
              resolution_counter = 1L;
              receptacle[0]=7;
              break;
              
            case 'K':
              resolution_counter = 1L;
              receptacle[0]=8;
              break;
              
            case 'L':
              resolution_counter = 1L;
              receptacle[0]=9;
              break;
              
            case 'M':
              resolution_counter = 1L;
              receptacle[0]=10;
              break;
              
            case 'N':
              resolution_counter = 1L;
              receptacle[0]=11;
              break;
              
            case 'P':
              resolution_counter = 1L;
              receptacle[0]=12;
              break;
              
            case 'Q':
              resolution_counter = 1L;
              receptacle[0]=13;
              break;
              
            case 'R':
              resolution_counter = 1L;
              receptacle[0]=14;
              break;
              
            case 'S':
              resolution_counter = 1L;
              receptacle[0]=15;
              break;
              
            case 'T':
              resolution_counter = 1L;
              receptacle[0]=16;
              break;
              
            case 'V':
              resolution_counter = 1L;
              receptacle[0]=17;
              break;
              
            case 'W':
              resolution_counter = 1L;
              receptacle[0]=18;
              break;
              
            case 'Y':
              resolution_counter = 1L;
              receptacle[0]=19;
              break;
              
            case 'Z':
              resolution_counter = 2L;
              receptacle[0]=3;
              receptacle[1]=13;
              break;
              
            case 'X':
            case '?':
            case '.':
            case '*': {
              resolution_counter = 20L;
              for (unsigned long j = 0UL; j<20UL; j++) {
                receptacle[j] = j;
              }
            }
              break;
            case '-': {
              resolution_counter = 0L;
            }
              break;
          }
        } else
          // binary
        {
          
          switch (token) {
            case '0':
              resolution_counter = 1L;
              receptacle[0]=0;
              break;
              
            case '1':
              resolution_counter = 1L;
              receptacle[0]=1;
              break;
              
            case 'X':
            case '?':
            case '.':
            case '*': {
              resolution_counter = 2L;
              receptacle[0] = 0;
              receptacle[1] = 1;
            }
              break;
            case '-': {
              resolution_counter = 0L;
            }
              break;
          }
          
        }
      }
    }
  }
  
  if (resolution_counter == 0L && gapToOnes) {
    for (unsigned long i = 0UL; i < baseLength; i++) {
      receptacle[i] = i;
    }
    return baseLength;
  }
  
  return resolution_counter;
  
}

//_________________________________________________________
void    _TranslationTable::PrepareForChecks (void) {
    if (checkTable == NULL) {
        checkTable = MemAllocate (256);
    }

    for (long i2=0; i2<256; i2++) {
        checkTable[i2]=0;
    }

    _String checkSymbols;
//  if (baseLength == 4)
//      checkSymbols = _String("ACGTUYRWSKMBDHVXN?O-.")&tokensAdded;
    if (baseSet.sLength) {
        checkSymbols = baseSet&tokensAdded;
    } else if (baseLength == 2) {
        checkSymbols = _String("01*?-.")&tokensAdded;
    } else {
        checkSymbols = _String("ABCDEFGHIJKLMNOPQRSTUVWXYZ*?-.")&tokensAdded;
    }

    for (long i=0; i<checkSymbols.sLength; i++) {
        checkTable[(unsigned char)checkSymbols(i)] = 1;
    }
}

//_________________________________________________________
bool    _TranslationTable::IsCharLegal (char c)
{
    if (!checkTable) {
        PrepareForChecks();
    }
    return checkTable[(unsigned char)c];
}

const _String& _TranslationTable::GetAlphabetString (void) const {
  if (baseSet.sLength) {
    return baseSet;
  }
  
  if (baseLength == 4) {
    return _TranslationTable::GetDefaultTable(HY_TRANSLATION_TABLE_DNA);
  } else if (baseLength == 20) {
    return _TranslationTable::GetDefaultTable(HY_TRANSLATION_TABLE_PROTEIN);
  } else {
    return _TranslationTable::GetDefaultTable(HY_TRANSLATION_TABLE_BINARY);
  }
  
  return emptyString;

}

//___________________________________________

void    _TranslationTable::AddTokenCode (char token, _String& code) {
    long    f,
            newCode = 0;

    bool    killBS = false;

    if (baseSet.sLength==0)
        // fill in baseSet for standard alphabets
    {
        if (baseLength == 4) {
          baseSet = _TranslationTable::GetDefaultTable(HY_TRANSLATION_TABLE_DNA);
        } else if (baseLength == 20) {
          baseSet = _TranslationTable::GetDefaultTable(HY_TRANSLATION_TABLE_PROTEIN);
        } else {
         baseSet = _TranslationTable::GetDefaultTable(HY_TRANSLATION_TABLE_BINARY);
        }
        killBS = true;
    }


    if (baseSet.sLength) {
        long shifter = 1;
        for (int j = 0; j<baseSet.sLength; j++, shifter*=2)
            if (code.Find (baseSet.sData[j])>=0) {
                newCode += shifter;
            }
    }

    f = baseSet.Find (token);
    if (killBS) {
        baseSet = emptyString;
    }
    if (f>=0) {
        return;
    }
    // see if the character being added is a base
    // character; those cannot be redefined

    f = tokensAdded.Find (token,0,-1);
    // new definition or redefinition?

    if (f==-1) { // new
        tokensAdded             = tokensAdded&token;
        translationsAdded       << 0;
        f                       = tokensAdded.sLength-1;
    }

    translationsAdded.lData[f] = newCode;
}

//_________________________________________________________

void    _TranslationTable::AddBaseSet (_String& code)
{
    baseSet         = code;
    baseSet.StripQuotes();
    baseLength      = baseSet.sLength;
    if (baseLength > HY_WIDTH_OF_LONG)
        // longer than the bit size of 'long'
        // can't handle those
    {
        _String err = _String ("Alphabets with more than ")
                      & HY_WIDTH_OF_LONG &
                      " characters are not supported";
        WarnError (err);
    }

}

//_________________________________________________________

char    _TranslationTable::GetSkipChar (void)
{
    if ( baseSet.sLength==0 && translationsAdded.lLength==0 ) {
        return '?';    // this is the default
    }

    // see if there is a symbol
    // which maps to all '1'

    long    all     = 0,
            ul       = baseSet.sLength?baseSet.sLength:baseLength,
            shifter = 1;

    for  (long f=0; f<ul; f++, shifter <<= 1) {
        all |= shifter;
    }

    if  ((all = translationsAdded.Find(all))==-1) {
        return '?';
    } else {
        return tokensAdded[all];
    }
}

//_________________________________________________________

char    _TranslationTable::GetGapChar (void) const {
    if ( baseSet.sLength==0 && translationsAdded.lLength==0 ) {
        return '-';    // default gap character
    }

    long f = translationsAdded.Find(0L);

    return f >= 0 ? tokensAdded[f] : '\0';
}

//_________________________________________________________
const _String _TranslationTable::ConvertCodeToLetters (long code, unsigned char base) const {

    _String res (base,false);
    if (code >= 0) {
        // OPTIMIZE FLAG; repeated memory allocation/deallocation
        if (baseSet.sLength)
            for (long k=1; k<=base; k++, code/=baseLength) {
                res.sData[base-k]=baseSet.sData[code%baseLength];
            }
        else if (baseLength==4) {
            for (long k=1; k<=base; k++, code/=baseLength) {
                switch (code%baseLength) {
                case 0:
                    res[base-k]='A';
                    break;
                case 1:
                    res[base-k]='C';
                    break;
                case 2:
                    res[base-k]='G';
                    break;
                case 3:
                    res[base-k]='T';
                    break;
                }
            }
        } else if (baseLength == 20) {
            for (long k=1; k<=base; k++, code/=baseLength) {
                char out = code%baseLength;
                if (out==0) {
                    res[base-k] = 'A';
                } else if (out<=7) {
                    res[base-k] = 'B'+out;
                } else if (out<=11) {
                    res[base-k] = 'C'+out;
                } else if (out<=16) {
                    res[base-k] = 'D'+out;
                } else if (out<=18) {
                    res[base-k] = 'E'+out;
                } else {
                    res[base-k]='Y';
                }
            }
        } else if (baseLength == 2)
            for (long k=1; k<=base; k++, code/=baseLength) {
                switch (code%baseLength) {
                case 0:
                    res[base-k]='0';
                    break;
                case 1:
                    res[base-k]='1';
                    break;
                }
            }
    } else {
        char c = GetGapChar();
        for (long k=0L; k<base; k++) {
            res.sData[k] = c;
        }
    }
    return res;


}

//_________________________________________________________

_TranslationTable*  _TranslationTable::MergeTables (_TranslationTable const* table2) const
// merge the translation tables if they are compatible, return the result,
// otherwise return nil
{
    if (baseSet.sLength==table2->baseSet.sLength) {
        if (baseSet.sLength==0) { // standard alphabet
            if (baseLength!=table2->baseLength) {
                return nil;
            }
        } else if(!(baseSet.Equal (&table2->baseSet))) {
            return nil;
        }

        _TranslationTable* result = new _TranslationTable (*this);
        if (table2->tokensAdded.sLength) {
            for (long i=0; i<table2->tokensAdded.sLength; i++) {
                long f = tokensAdded.Find (table2->tokensAdded[i]);
                if (f==-1) {
                    result->tokensAdded       && table2->tokensAdded[i];
                    // SLKP 20071002 added the next line;
                    // was not adding the translation for the new token
                    result->translationsAdded << table2->translationsAdded(i);
                } else if (translationsAdded.lData[f] != table2->translationsAdded.lData[i]) {
                    DeleteObject (result);
                    return nil;
                }
            }
            return result;
        } else {
            return result;
        }
    }
    return nil;
}

//_________________________________________________________

_Site::_Site (void):_CString(16,true)
{
    refNo = -1;
}

//_________________________________________________________
_Site::_Site (_String& s):_CString (s.sLength, true)
{
    refNo = -1;
    (*this)<<&s;
}

//_________________________________________________________
_Site::_Site (char s):_CString (16, true)
{
    refNo = -1;
    (*this)<<s;
}

//_________________________________________________________
_Site::_Site (long s)
//:_CString (1, true)
{
    SetRefNo (s);
}

//_________________________________________________________
_Site::~_Site (void)
{}

//_________________________________________________________
void    _Site::Complete (void)
{
    if (refNo==-1) {
        _String::Finalize();
    }

    refNo = refNo<0?-refNo:refNo;
}
//_________________________________________________________
BaseRef _Site::makeDynamic(void)
{
    _Site * r = new _Site;
    checkPointer(r);

    memcpy ((char*)r, (char*)this, sizeof (_Site));
    r->nInstances = 1;
    nInstances++;
    return r;
}

//_______________________________________________________________________
void    _Site::Duplicate (BaseRef ref)
{
    _Site * s = (_Site*)ref;
    sLength = s->sLength;
    if (sData) {
        free(sData);
    }
    sData      = s->sData;
    allocatedSpace = s->allocatedSpace;
    //nInstances = ref->nInstances;
    if (sData) {
        /*long theLength = sLength/storageIncrement;
        if (!sLength||sLength%storageIncrement) theLength++;
        theLength*=storageIncrement;
        checkPointer (sData = (char*)MemAllocate (theLength));
        memcpy (sData, s->sData, sLength);*/
        if (allocatedSpace) {
            checkPointer (sData = (char*)MemAllocate (allocatedSpace*sizeof(char)));
        } else {
            checkPointer (sData = (char*)MemAllocate (sLength*sizeof(char)));
        }
        memcpy (sData, s->sData, sLength);
    }
    refNo = -1;
}

//_______________________________________________________________________
void    _Site::Clear (void)
{
    if (sData) {
        free(sData);
        sData = NULL;

        //nInstances = 0;
    }
    allocatedSpace = 0;
    sLength = 0;
}
//_______________________________________________________________________
void    _Site::PrepareToUse (void)
{
    if (IsCompressed()) {
        _String * s = Decompress();
        DuplicateErasing(s);;
        DeleteObject (s);
        SetDecompressed();
    }
}
//_______________________________________________________________________
void    _Site::Archive (void)
{
    if ((!IsCompressed())&&(GetRefNo()>=0)) {
        BestCompress (NUCLEOTIDEALPHABET);
    }
}


//_______________________________________________________________________
//  _DataSet Functions
//_______________________________________________________________________

_DataSet::_DataSet (void)
{
    theTT                = &defaultTranslationTable;
    streamThrough        = nil;
    dsh                  = nil;
    useHorizontalRep     = false;
}

//_______________________________________________________________________

_DataSet::_DataSet (long l):_List((unsigned long)l),theFrequencies((unsigned long)l) // with estimated number of sites per file
{
    dsh                 = nil;
    streamThrough       = nil;
    theTT               = &defaultTranslationTable;
    useHorizontalRep     = false;
}

//_______________________________________________________________________

_DataSet::_DataSet (FILE *f)
{
    dsh                  = nil;
    useHorizontalRep     = false;
    theTT                = &defaultTranslationTable;
    streamThrough        = f;
    theMap << 0; // current sequence
    theMap << 0; // current site
    theMap << 0; // total sites
}

//_______________________________________________________________________

_DataSet::~_DataSet (void)
{
    if (theTT != &defaultTranslationTable) {
        DeleteObject (theTT);
    }
}

//_______________________________________________________________________

void _DataSet::Clear (bool )
{
    _List::Clear();
    theMap.Clear();
    theFrequencies.Clear();
    theNames.Clear();
    if (theTT != &defaultTranslationTable) {
        DeleteObject (theTT);
        theTT = &defaultTranslationTable;
    }
    noOfSpecies = 0;
    if (dsh) {
        dsh->incompletePatterns->Clear(false);
        delete (dsh);
        dsh = nil;
    }
    useHorizontalRep     = false;
}

//_______________________________________________________________________

BaseRef _DataSet::makeDynamic (void)
{
    _DataSet * r = new _DataSet;
    checkPointer(r);
    memcpy ((char*)r, (char*)this, sizeof (_DataSet));
    r->nInstances = 1;
    r->theMap.Duplicate (&theMap);
    r->theFrequencies.Duplicate (&theFrequencies);
    if (theTT!=&defaultTranslationTable) {
        r->theTT->nInstances++;
    }
    r->theNames.Duplicate (&theNames);
    r->streamThrough = streamThrough;
    nInstances++;
    r->dsh = nil;
    r->useHorizontalRep      = false;
    return r;
}

//_______________________________________________________________________

void     _DataSet::ResetIHelper (void)
{
    if (dsh && dsh->characterPositions.lLength == 256)
        for (long k=0; k<256; k++) {
            dsh->characterPositions.lData[k] = -1;
        }
}

//_______________________________________________________________________

void     _DataSet::ConvertRepresentations (void)
{
    if (useHorizontalRep == false) {
        _List horStrings;

        if (lLength == 0) {
            AppendNewInstance (new _Site);
        } else {
            _Site * aSite = (_Site*)lData[0];

            for (long str = 0; str < aSite->sLength; str++) {
                _String * aString = new _String (DATA_SET_SWITCH_THRESHOLD,true);
                horStrings << aString;
                aString->nInstances --;
            }

            for  (long s = 0; s < lLength; s++) {
                _Site * aSite = (_Site*)lData[s];
                if (aSite->sLength>horStrings.lLength || aSite->GetRefNo() != -1) {
                    FlagError ("Irrecoverable internal error in _DataSet::ConvertRepresentations. Sorry about that.");
                    return;
                }
                aSite->Finalize();
                for (long s2 = 0; s2 < aSite->sLength; s2++) {
                    (*(_String*)horStrings.lData[s2]) << aSite->sData[s2];
                }
            }

            _List::Clear();
            theFrequencies.Clear();
            {
                for  (long s = 0; s < horStrings.lLength; s++) {
                    (*this) << horStrings(s);
                }
            }
        }
        useHorizontalRep = true;
    }
}

//_______________________________________________________________________

void     _DataSet::AddSite (char c)
{
    if (streamThrough) {
        if (theMap.lData[0] == 0) {
            if (theMap.lData[1] == 0) {
                if (theNames.lLength) {
                    fprintf (streamThrough,">%s\n",((_String*)theNames(0))->getStr());
                } else {
                    fprintf (streamThrough,">Sequence 1\n");
                }
                (*this) && & emptyString;
            }

            theMap.lData[1]++;
            theMap.lData[2]++;
            fputc (c,streamThrough);
        } else {
            WarnError ("Can't add more sites to a file based data set, when more that one sequence has been written!");
        }
    } else {
        if (useHorizontalRep == false) {
            if (lLength < DATA_SET_SWITCH_THRESHOLD) {
                _Site* nC = new _Site(c);
                theFrequencies<<1L;
                AppendNewInstance(nC);
                return;
            } else {
                ConvertRepresentations ();
            }
        }

        (*((_String*)lData[0])) << c;

        /*long  f;

        if (!dsh)
        {
            checkPointer (dsh = new _DSHelper);
            for (f=0; f<256; f++)
                dsh->characterPositions << -1;
        }

        f = dsh->characterPositions.lData[c];

        if (f!=-1)
        {
            _Site* nC = new _Site(f);
            checkPointer(nC);
            theFrequencies[f]++;
            theFrequencies<<0;
            (*this)<<nC;
            nC->nInstances --;
        }
        else
        {
            dsh->characterPositions.lData[c] = lLength;*/
        //}
    }
}
//_______________________________________________________________________

void     _DataSet::Write2Site (long index, char c)
{
    if (streamThrough) {
        if (index == 0) {
            if (theMap.lData[2] == theMap.lData[1]) {
                theMap.lData[0] ++;

                if (theNames.lLength > theMap.lData[0]) {
                    fprintf (streamThrough,"\n>%s\n",((_String*)theNames(theMap.lData[0]))->getStr());
                } else {
                    fprintf (streamThrough,"\n>Sequence %ld\n",theMap.lData[0]+1);
                }

                theMap.lData[1] = 0;
            } else {
                WarnError ("Can't write sequences of unequal lengths to a file based data set.");
                return;
            }
        } else if (index != theMap.lData[1]) {
            WarnError ("Can't write sites which are not consecutive to a file based data set.");
            return;
        }

        theMap.lData[1] ++;
        fputc (c,streamThrough);
    } else {
        /*if (!dsh)
        {
            WarnError ("Internal Error in 'Write2Site' - called Write2Site before any AddSite calls");
            return;
        }*/

        if (useHorizontalRep) {
            long     currentWritten = ((_String*)lData[0])->sLength;

            if (index>=currentWritten) {
                WarnError ("Internal Error in 'Write2Site' - index is too high (using compact representation)");
                return;
            } else {
                if (index == 0) {
                    _String * newString = new _String (currentWritten,true);
                    (*newString) << c;
                    (*this) << newString;
                    newString->nInstances --;
                } else {
                    long s = 1;
                    for (; s<lLength; s++) {
                        _String *aString = (_String*)lData[s];
                        if (aString->sLength == index) {
                            (*aString) << c;
                            break;
                        }
                    }
                    if (s == lLength) {
                        WarnError ("Internal Error in 'Write2Site' - no appropriate  string to write too (compact representation)");
                        return;
                    }
                }
            }
        } else {
            if (index>=lLength) {
                WarnError ("Internal Error in 'Write2Site' - index is too high");
                return;
            }
            _Site* s = (_Site*)lData[index];
            long rN = s->GetRefNo();
            if (rN==-1) { // independent site
                //dsh->incompletePatterns->Delete (s,false);
                (*s)<<c;
                //dsh->incompletePatterns->Insert (s,index);
            } else {
                _Site *ss = (_Site*)lData[rN];
                long sL = ss->sLength-1;
                if (ss->sData[sL] != c) { // appending distinct char
                    s->Duplicate (ss);
                    s->sData[sL]=c;
                    theFrequencies.lData[rN]--;

                    rN = dsh->incompletePatterns->Find (s);
                    if (rN>=0) {
                        rN =  dsh->incompletePatterns->GetXtra (rN);
                        /*_Site* s2 = (_Site*)lData[rN];
                        if (s2->GetRefNo() != -1 || !s->Equal(s2))
                        {
                            WarnError ("Mapping Error");
                        }*/
                        theFrequencies[rN]++;
                        s->Clear();
                        s->SetRefNo(rN);
                    } else {
                        theFrequencies[index]++;
                        s->SetRefNo(-1);
                        dsh->incompletePatterns->Insert (s,index);
                    }
                }
            }
        }
    }
}



//_______________________________________________________________________

void     _DataSet::CheckMapping (long index)
{
    if (index>=lLength) {
        FlagError ("Internal Error in 'CheckMapping' - index is too high");
    }

    _Site* s = (_Site*)lData[index];

    for (long k = 0; k < index; k ++) {
        _Site* ss = (_Site*)lData[k];
        if (ss->GetRefNo() == -1) {
            if (s->Equal(ss)) {
                theFrequencies[index]--;
                theFrequencies[k]++;
                s->Clear();
                s->SetRefNo(k);
            }
        }
    }
}

//_______________________________________________________________________

unsigned long     _DataSet::GetCharDimension         (void) const {  // return the size of the alphabet space
    return theTT->baseLength;
}

//_______________________________________________________________________

long     _DataSet::GetNoTypes (void) const // return the number of unique columns
{
    return theMap.countitems ();
}
//_______________________________________________________________________

unsigned long     _DataSet::GetFreqType (long index) const  { // return the frequency of a site
    return theFrequencies(theMap(index));
}
//_______________________________________________________________________

void    _DataSet:: SetTranslationTable (_DataSet *  newTT )
{
    if (theTT&&(theTT!= &defaultTranslationTable)) {
        DeleteObject (theTT);
    }
    theTT = (_TranslationTable*)newTT->theTT->makeDynamic();

}

//_______________________________________________________________________

void    _DataSet:: SetTranslationTable (_TranslationTable *  newTT )
{
    if (theTT&&(theTT!= &defaultTranslationTable)) {
        DeleteObject (theTT);
    }
    theTT = (_TranslationTable*)newTT->makeDynamic();

}
//_______________________________________________________________________
void    _DataSet::Finalize (void)
{
    if (streamThrough) {
        fclose (streamThrough);
        streamThrough = nil;
        theMap.Clear();
    } else {
        if (useHorizontalRep) {
            bool  good = true;
            for (long s = 0; s<lLength; s++) {
                ((_String*)lData[s])->Finalize();
                good = good && ((_String*)lData[0])->sLength == ((_String*)lData[s])->sLength;
            }

            if (!good) {
                Clear();
                WarnError ("Internal Error in _DataSet::Finalize. Unequal sequence lengths in compact representation");
                return;
            }

            _List               dups;
            _List               uniquePats;
            _AVLListX           dupsAVL (&dups);

            long  siteCounter = ((_String*)lData[0])->sLength;

            for (long i1 = 0; i1<siteCounter; i1++) {
                _Site * tC = new _Site ();
                checkPointer (tC);

                for (long i2 = 0; i2 < lLength; i2++) {
                    (*tC) << ((_String*)lData[i2])->sData[i1];
                }

                tC->Finalize();

                long ff = dupsAVL.Find(tC);
                if (ff<0) {
                    uniquePats << tC;
                    dupsAVL.Insert(tC,theFrequencies.lLength);
                    theMap         << theFrequencies.lLength;
                    theFrequencies << 1;
                } else {
                    ff =  dupsAVL.GetXtra(ff);
                    theMap << ff;
                    theFrequencies.lData[ff] ++;
                }

                DeleteObject (tC);
            }
            dupsAVL.Clear(false);
            _List::Clear();
            _List::Duplicate(&uniquePats);
        } else {
            long j,
                 k;

            _Site         *tC;
            {
                _List          dups;
                _AVLListX      dupsAVL (&dups);

                for (long i1 = 0; i1<lLength; i1++) {
                    tC = (_Site*)lData[i1];
                    long ff = dupsAVL.Find(tC);
                    if (ff<0) {
                        dupsAVL.Insert(tC,i1);
                    } else {
                        ff =  dupsAVL.GetXtra(ff);
                        tC->Clear();
                        tC->SetRefNo(ff);
                        theFrequencies.lData[ff] ++;
                    }
                }
                dupsAVL.Clear(false);
            }

            _SimpleList  refs(lLength),
                         toDelete(lLength);
            j = 0;

            for (long i1 = 0; i1<lLength; i1++) {
                tC = (_Site*)(*(_List*)this)(i1);
                k = tC->GetRefNo();
                if (k==-1) {
                    refs << j++;
                } else {
                    toDelete << i1;
                    refs << -1;
                }
            }

            for (long i2=0; i2<lLength; i2++) {
                tC = (_Site*)(*(_List*)this)(i2);
                k = tC->GetRefNo();
                if (k>=0) {
                    j = refs.lData[k];
                    if (j<0) {
                        warnError (-171);
                    } else {
                        refs.lData[i2]=j;
                    }
                }
            }

            theMap.Clear();
            theMap.Duplicate (&refs);
            DeleteList (toDelete);
            theFrequencies.DeleteList (toDelete);

            for (long i3 = 0; i3<lLength; i3++) {
                tC = (_Site*)(*(_List*)this)(i3);
                tC->SetRefNo (0);
                tC->Finalize();
            }
            if (dsh) {
                dsh->incompletePatterns->Clear(false);
                delete (dsh);
                dsh = nil;
            }
        }
    }
}
//_______________________________________________________________________
void    _DataSet::Compact (long index)
{
    if (useHorizontalRep) {
        WarnError ("Internal Error: _DataSet::Compact called with compact represntation");
        return;
    }

    _Site* tC = (_Site*)(*(_List*)this)(index);
    if (tC->GetRefNo()!=-1)
        // take care of double referencing
    {
        _Site*tCC=tC;
        long lastRef, count = 0;
        do {
            lastRef = tCC->GetRefNo();
            count++;
            tCC = (_Site*)(*(_List*)this)(tCC->GetRefNo());
        } while (tCC->GetRefNo()!=-1);
        if (count>1) {
            theFrequencies[lastRef]++;
        }

        tC->SetRefNo(lastRef);
    }
    /*if (tC->GetRefNo()==-1)
    {
     long f = dsh->incompletePatterns->Find(tC);
     if (f >= 0)
     {
            f = dsh->incompletePatterns->GetXtra (f);
            if (f<index)
            {
            theFrequencies[f]++;
            tC->Clear();
            tC->SetRefNo(f);
        }
        else
            tC->Finalize();
     }
    }*/
}


//_______________________________________________________________________
inline char  _DataSet::operator () (unsigned long site, unsigned long pos, unsigned int) const {
    return (((_String**)lData)[theMap.lData[site]])->sData[pos];
}


//_________________________________________________________
long _DataSet::ComputeSize(void)
{
    long res = sizeof (_DataSet);

    res+=(theMap.lLength+lLength+theFrequencies.lLength)*sizeof(long);
    res+=lLength*sizeof(_Site);

    for (long i=0; i<lLength; i++) {
        res+= ((_Site*)(*(_List*)this)(i))->sLength;
    }

    return res;
}

//_________________________________________________________
_Parameter _DataSet::CheckAlphabetConsistency(void)
{
    long        charsIn = 0,
                gaps    = 0,
                total   = 0;

    bool        checks    [256];
  
    char        gapChar = theTT->GetGapChar();

    _String     baseSymbols;

    if (theTT->baseSet.sLength) {
        baseSymbols = theTT->baseSet;
    } else if (theTT->baseLength == 4) {
        baseSymbols = "ACGUT";
    } else if (theTT->baseLength == 20) {
      baseSymbols = _TranslationTable::GetDefaultTable(HY_TRANSLATION_TABLE_PROTEIN);
    } else {
      baseSymbols = _TranslationTable::GetDefaultTable(HY_TRANSLATION_TABLE_BINARY);
    }



    for (; charsIn<256; charsIn++) {
        checks[charsIn] = false;
    }

    for (charsIn=0; charsIn<baseSymbols.sLength; charsIn++) {
        checks[(unsigned char)baseSymbols.sData[charsIn]] = true;
    }

    charsIn = 0;

    for (long i = 0; i<lLength; i++) {
        _String* thisColumn = (_String*)lData[i];
        long     w = theFrequencies.lData[i];
        for (long j = 0; j<thisColumn->sLength; j++)
            if (checks[(unsigned char)thisColumn->sData[j]]) {
                charsIn+=w;
            } else if (gapChar == thisColumn->sData[j]) {
                gaps += w;
            }

        total += w*thisColumn->sLength;
    }

    return (_Parameter)charsIn/(total-gaps+1.);

}

//___________________________________________________

BaseRef _DataSet::toStr (unsigned long)
{
    _String * s = new _String(NoOfSpecies()*30, true),
    *str;

    checkPointer (s);
    (*s) << _String ((long)NoOfSpecies());
    (*s) << " species:";

    str = (_String*)theNames.toStr();
    (*s) << *str;
    DeleteObject(str);

    (*s) << ";\nTotal Sites:";
    (*s) << _String((long)GetNoTypes());
    (*s) << ";\nDistinct Sites:";
    (*s) << _String((long)theFrequencies.lLength);

    s->Finalize();

    return s;
}

//___________________________________________________

void    _DataSet::toFileStr (FILE* dest, unsigned long padding)
{
    fprintf (dest, "%ld species: ",NoOfSpecies());
    theNames.toFileStr(dest, padding);

    fprintf (dest, ";\nTotal Sites: %ld",GetNoTypes()) ;
    fprintf (dest, ";\nDistinct Sites: %ld",theFrequencies.lLength);

    /*  fprintf (dest,"\n");
        for (long j=0; j<noOfSpecies;j++)
        {
            fprintf (dest,"\n");
            for (long i=0; i<theMap.lLength; i++)
            {
                fprintf (dest,"%c",(*this)(i,j,1));
            }
        }*/

}


//_________________________________________________________

void    _DataSet::AddName (_String const& s) {
    theNames.AppendNewInstance (new _String (s,0,s.FirstNonSpaceIndex (0,-1,-1)));
}

  //_________________________________________________________

void    _DataSet::InsertName (_String const& name, long where ) {
  theNames.InsertElement (new _String (name), where, false);
}




//_________________________________________________________

void    _DataSet::MatchIndices (_Formula&f, _SimpleList& receptacle, bool isVert, long limit, _String const* scope) const {
    _String     varName  = isVert ? "siteIndex" : "speciesIndex";
    varName = AppendContainerName(varName, scope);
    _Variable   *v       = CheckReceptacle (&varName, emptyString, false);
  
    //fprintf (stderr, "\n_DataSet::MatchIndices %d %s [%s] %s\n", isVert, scope ? scope->sData : "none", varName.sData, ((_String*)f.toStr())->sData);

    for (long i=0L; i<limit; i++) {
        v->SetValue (new _Constant((_Parameter)i), nil);
        _PMathObj res = f.Compute();
        //fprintf (stderr, "%ld %g\n", i, res->Compute()->Value());
        if (res && !CheckEqual(res->Value(),0.0)) {
            receptacle<<i;
        }
    }
    v->SetValue (new _Constant(0.0), nil);
}

//_________________________________________________________

void    _DataSet::FindAllSitesLikeThisOne (long index, _SimpleList& receptacle)
{
    if (index>=0 && index<theMap.lLength) {
        index = theMap.lData[index];
        for (long k=0; k<theMap.lLength; k++)
            if (theMap.lData[k]==index) {
                receptacle << k;
            }
    }
}

//_________________________________________________________

_TranslationTable*      _DataSet::CheckCompatibility (_SimpleList const & ref, char concatOrCombine) {
    _DataSet* currentSet = (_DataSet*)dataSetList(ref.Element (0));

    _TranslationTable* theEnd = new _TranslationTable (*(currentSet->theTT));
    checkPointer(theEnd);
    long    refNo     =  concatOrCombine?currentSet->NoOfSpecies():currentSet->NoOfColumns();
    char    emptyStringChar = theEnd->GetSkipChar();

    for (long k=1; k<ref.lLength; k++) {
        currentSet = (_DataSet*)dataSetList(ref.Element(k));

        _TranslationTable* tryMe = theEnd->MergeTables (currentSet->theTT);

        if (tryMe) {
            if (emptyStringChar) {
                DeleteObject (theEnd);
                theEnd = tryMe;
                continue;
            } else {
                if ((concatOrCombine&&(currentSet->NoOfSpecies()==refNo))||(!concatOrCombine&&(currentSet->NoOfColumns()==refNo))) {
                    DeleteObject (theEnd);
                    theEnd = tryMe;
                    continue;
                }
            }
        }
        _String warningMessage ("The data set:");
        warningMessage = warningMessage & *((_String*)dataSetNamesList(ref.Element (k))) & _String (" was found incompatible with one of the following data sets:");
        for (long i=0; i<k; i++) {
            warningMessage = warningMessage & *((_String*)dataSetNamesList(ref.Element (k))) & _String (",");
        }
        WarnError (warningMessage);
        DeleteObject (theEnd);
        return nil;
    }

    return theEnd;
}

//_________________________________________________________

_DataSet*   _DataSet::Concatenate (_SimpleList const & ref)

// concatenates (adds columns together) several datasets
// in case the number of species in the datasets are different the deficiencies will be padded
// by omission symbols
// in case translation tables are different, they will be merged, provided it can be done,
// otherwise the incompatible datasets will be ignored during this operation.

{
    _TranslationTable  * jointTable =  CheckCompatibility (ref,1);
    if (!jointTable) {
      return new _DataSet;
    }

  
  
    _DataSet           * bigDataSet = new _DataSet;
    checkPointer(bigDataSet);

    bigDataSet->theTT = jointTable;

    // pass one - determine the max max number of species present and what dataset are they coming from

    long      maxSpecies=0,
              maxDataSet=0,
              siteIndex;

    _DataSet *currentSet;

    char     emptyStringSlot = jointTable->GetSkipChar();

    for (long i=0; i<ref.lLength; i++) {
        currentSet = (_DataSet*)dataSetList(ref(i));

        long       specCount = currentSet->NoOfSpecies(),
                   siteCount = currentSet->NoOfColumns();


        if (specCount>maxSpecies) {
            maxSpecies = specCount;
            maxDataSet = i;
        }
        for (long j=0; j<siteCount; j++) {
            bigDataSet->AddSite((*currentSet)(j,0,1));
        }
    }

    for (long k=1; k<maxSpecies; k++) {
        siteIndex = 0;
        for (long i=0; i<ref.lLength; i++) {
            currentSet = (_DataSet*)dataSetList(ref.lData[i]);

            long       cns = currentSet->NoOfSpecies(),
                       cnc = currentSet->NoOfColumns();

            if (cns<=k)
                for (long j=0; j< cnc; j++, siteIndex++) {
                    bigDataSet->Write2Site(siteIndex,emptyStringSlot);
                }
            else
                for (long j=0; j< cnc; j++, siteIndex++) {
                    bigDataSet->Write2Site(siteIndex,(*currentSet)(j,k,1));
                }
        }
    }

    currentSet = (_DataSet*)dataSetList(ref(maxDataSet));
    for (long i=0L; i<maxSpecies; i++) {
        bigDataSet->AddName (*currentSet->GetSequenceName(i));
    }

    bigDataSet->Finalize();
    bigDataSet->SetNoSpecies (maxSpecies);
    return bigDataSet;
}

//_________________________________________________________

_DataSet*   _DataSet::Combine (_SimpleList const& ref) {

// combines (adds rows together) several datasets
// in case the number of species in the datasets are different the deficiencies will be padded
// by omission symbols
// in case translation tables are different, they will be merged, provided it can be done,
// otherwise the incompatible datasets will be ignored during this operation.

    _TranslationTable  * joint_table = CheckCompatibility (ref,0);

    if (!joint_table) {
      return new _DataSet;
    }
  
    _DataSet           * combined_data = new _DataSet;
    combined_data->theTT = joint_table;

    // pass one - determine the max max number of sites present and what dataset are they coming from

    unsigned long     max_sites = 0UL,
                      total_species_count = 0UL;

  

    char     emptyStringSlot = joint_table->GetSkipChar();

    for (unsigned long set_index =0UL; set_index <ref.lLength;  set_index++) {
        _DataSet const *current_data_set = (_DataSet const*)dataSetList(ref.Element (set_index));
        StoreIfGreater (max_sites, current_data_set->NoOfColumns());
        total_species_count += current_data_set->NoOfSpecies();
    }


    for (unsigned long set_index = 0UL; set_index < ref.lLength; set_index++) {
        _DataSet const *current_data_set = (_DataSet const*)dataSetList(ref.Element (set_index));
        unsigned long sites_in_this_set = current_data_set->NoOfColumns(),
                      sequences_in_this_set = current_data_set->NoOfSpecies();
      
        for (unsigned long seq_index = 0UL; seq_index < sequences_in_this_set; seq_index++) {
            combined_data->AddName (*current_data_set->GetSequenceName(seq_index));
            if (seq_index == 0UL && set_index == 0UL) {
                /* use AddSite write out the first sequence */
                unsigned long site_index = 0UL;
                for (site_index = 0UL; site_index < sites_in_this_set; site_index ++) {
                    combined_data->AddSite ((*current_data_set)(site_index,0UL,1));
                }
                for (; site_index < max_sites ; site_index++) {
                    combined_data->AddSite (emptyStringSlot);
                }
            } else {
                /* use Write2Site to create subsequence sequences */
                unsigned long site_index = 0UL;
                for (site_index = 0UL; site_index < sites_in_this_set; site_index ++) {
                  combined_data->Write2Site (site_index, (*current_data_set)(site_index,seq_index,1));
                }
                for (; site_index < max_sites ; site_index++) {
                  combined_data->Write2Site (site_index, emptyStringSlot);
                }
            }
        }
    }

    combined_data->Finalize();
    combined_data->SetNoSpecies(total_species_count);
    return combined_data;
}

//_________________________________________________________
// Data Set Filter/Numeric
//_________________________________________________________

_DataSetFilter::_DataSetFilter (void) {
    unitLength = 0;
    theData = NULL;
    accessCache = nil;
}
//_________________________________________________________
_DataSetFilter::_DataSetFilter (_DataSet* ds, char, _String&) {
    theData     = ds;
    accessCache = nil;
}
//_________________________________________________________
_DataSetFilter::~_DataSetFilter (void) {
    DeleteObject (accessCache);
}


//_________________________________________________________
_DataSetFilterNumeric::~_DataSetFilterNumeric (void) {
  const _DataSet* dummy_dataset = GetData();
    // SLKP TODO 20161028; this is very ugly but needs to be done to remove
    // dangling data set references
    long ds_id = dataSetList._SimpleList::Find ((long)dummy_dataset);
  if (ds_id >= 0) {
    KillDataSetRecord(ds_id);
  } else {
    WarnError ("Internal error in ~_DataSetFilterNumeric");
  }
}

//_______________________________________________________________________

_DataSetFilterNumeric::_DataSetFilterNumeric (_Matrix* freqs, _List& values, _DataSet *ds, long cc)
{
    unitLength      = 1;
    categoryCount   = cc;

    SetData         (ds);

    _SimpleList     baseFreqs;

    freqs->ConvertToSimpleList (baseFreqs);
    dimension =                ((_Matrix*)values(0))->GetVDim();

    theNodeMap.Populate         (ds->GetNames().lLength,0,1);
    theOriginalOrder.Populate   (((_Matrix*)values(0))->GetHDim()/categoryCount,0,1);

    //theMap.Populate           (theFrequencies.lLength,0,1);
    //theOriginalOrder.Populate     (theFrequencies.lLength,0,1);
    //duplicateMap.Populate     (theFrequencies.lLength,0,1);


    /*CreateMatrix (&probabilityVectors, theNodeMap.lLength, shifter,false,true, false);

    _Parameter   *storeHere = probabilityVectors.theData;
    for (long spec = 0; spec < theNodeMap.lLength; spec++)
    {
        _Matrix * specMatrix = (_Matrix*)values(spec);
        for (long site = 0; site < theFrequencies.lLength; site++)
            for (long state = 0; state < dimension; state++,storeHere++)
                //probabilityVectors.theData [shifter*spec + site*dimension+state]
                *storeHere= specMatrix->theData[site*dimension+state];
    }*/

    _List        siteScores;
    _AVLListXL   siteIndices(&siteScores);

    duplicateMap.RequestSpace (baseFreqs.lLength+1);

    //bool       startD = false;

    char buffer[255];

    for (long site =0; site <baseFreqs.lLength; site++) {
        _Parameter      testV = 0.0;

        for (long k=0; k<theNodeMap.lLength; k++) // sweep down the columns
            for (long state = 0; state < dimension; state++) {
                testV += ((_Matrix*)(((_Matrix**)values.lData)[k]))->theData[site*dimension+state];
            }

        snprintf (buffer, sizeof(buffer), "%20.18g", testV);
        _String     testS (buffer);
        long        f = siteIndices.Find (&testS);

        _SimpleList * sameScore = nil;

        if (f>=0) {
            sameScore = (_SimpleList*)siteIndices.GetXtra (f);
            for (long k = 0; k<sameScore->lLength; k++) {
                bool fit = true;
                f        = sameScore->lData[k];


                for (long spec=0; spec<theNodeMap.lLength && fit; spec++) { // sweep down the columns
                    _Matrix * specMatrix =(_Matrix*)(((_Matrix**)values.lData)[spec]);
                    for (long state = 0; state < dimension; state++)
                        if (specMatrix->theData[site*dimension+state]!=specMatrix->theData[theMap.lData[f]*dimension+state]) {
                            fit = false;
                            break;
                        }
                }

                if (fit) {
                    theFrequencies[f]+=baseFreqs[site];
                    duplicateMap<<f;
                    f = 0;
                    break;
                } else {
                    f = -1;
                }
            }
        }
        if (f==-1) {
            if (!sameScore) {
                checkPointer        (sameScore = new _SimpleList);
                if (siteIndices.Insert  (testS.makeDynamic(),(long)sameScore,false) < 0) {
                    _String eh ("WTF?");
                    StringToConsole(eh);
                }
            }

            (*sameScore) << theFrequencies.lLength;

            duplicateMap<<theFrequencies.lLength;
            theFrequencies<<baseFreqs[site];
            theMap<<site;
        }
    }

    siteIndices.Clear(true);
    shifter         = theFrequencies.lLength*dimension;
    categoryShifter = shifter*theNodeMap.lLength;

    CreateMatrix (&probabilityVectors, theNodeMap.lLength, shifter*categoryCount,false,true, false);
    _Parameter   *storeHere    = probabilityVectors.theData;

    long      refShifter = 0;
    for (long cc = 0; cc < categoryCount; cc++, refShifter += theOriginalOrder.lLength * dimension) {
        for (long spec = 0; spec < theNodeMap.lLength; spec++) {
            _Matrix * specMatrix = (_Matrix*)values(spec);
            for (long site = 0; site < theFrequencies.lLength; site++)
                for (long state = 0; state < dimension; state++,storeHere++) {
                    *storeHere = specMatrix->theData[refShifter + theMap.lData[site]*dimension+state];
                }
        }
    }
}


//_______________________________________________________________________

void _DataSetFilter::CopyFilter (_DataSetFilter *copyFrom) {
    memcpy ((char*)this, (char*)copyFrom, sizeof (_DataSetFilter));

    theFrequencies.Duplicate        (&copyFrom->theFrequencies);
    theNodeMap.Duplicate            (&copyFrom->theNodeMap);
    theMap.Duplicate                (&copyFrom->theMap);
    theOriginalOrder.Duplicate      (&copyFrom->theOriginalOrder);
    conversionCache.Duplicate       (&copyFrom->conversionCache);
    duplicateMap.Duplicate          (&copyFrom->duplicateMap);

    nInstances              = 1;
    dimension               = copyFrom->dimension;
    undimension             = copyFrom->undimension;
    unitLength              = copyFrom->unitLength;
    accessCache             = nil;

}


//_______________________________________________________________________

BaseRef _DataSetFilter::makeDynamic (void) {
    _DataSetFilter * r = new _DataSetFilter;
    r->CopyFilter   (this);

    return r;
}


//_______________________________________________________________________

BaseRef _DataSetFilterNumeric::makeDynamic (void) {
    _DataSetFilterNumeric * r = new _DataSetFilterNumeric();
    r->CopyFilter           (this);
    r->probabilityVectors.Duplicate(&probabilityVectors);
    return r;
}

//_______________________________________________________________________

_Parameter * _DataSetFilterNumeric::getProbabilityVector (long spec, long site, long categoryID) {
    return probabilityVectors.theData + categoryID * categoryShifter + spec * shifter + site * dimension;
}

//_______________________________________________________________________
long    _DataSetFilter::FreeUpMemory (long requestedBytes)
{
    long res = 0;
    for (long i=0; (i<theMap.lLength)&&(res<requestedBytes); i++) {
        res+=(theData->GetSite(theMap[i]))->FreeUpMemory(requestedBytes-res);
    }
    return res;
}


//_______________________________________________________________________
void    _DataSetFilter::SetDimensions (void)
{
    dimension   = GetDimension(true);
    undimension = GetDimension(false);
}

//_______________________________________________________________________
unsigned long    _DataSetFilter::FindUniqueSequences  (_SimpleList& indices, _SimpleList& map, _SimpleList& counts, short mode) const {
    indices.Clear(); map.Clear(); counts.Clear();
    
    unsigned long             sites  = theMap.lLength,
    seqs   = theNodeMap.lLength,
    unit   = GetUnitLength();
    
    if (mode == 0) {
        _SimpleList hashSupport;        
        _AVLListXL  sequenceHashes     (&hashSupport);
        
        for (unsigned long sequenceIndex = 0; sequenceIndex < seqs; sequenceIndex ++){
            _String * thisSequence = GetSequenceCharacters (sequenceIndex);
            
            long     sequenceHash   = thisSequence->Adler32(),
            f              = sequenceHashes.Find ((BaseRef)sequenceHash),
            rawSequenceIdx = theNodeMap.lData[sequenceIndex];
            
            DeleteObject (thisSequence);
            
            _SimpleList * sameScore = nil;
            if (f>=0) {
                sameScore = (_SimpleList*)sequenceHashes.GetXtra (f);
                for (long k = 0; k<sameScore->lLength; k++) {
                    bool fit = true;
                    f = sameScore->lData[k];
                    
                    long fRaw = theNodeMap.lData[indices.lData[f]];
                    
                    for (unsigned long site = 0; site < sites && fit; site++){
                        for (unsigned long unitIndex = 0; unitIndex < unit; unitIndex ++){
                            _Site * thisSite = theData->GetSite(theMap.lData[unit*site+unitIndex]);
                            if (thisSite->sData[fRaw] != thisSite->sData[rawSequenceIdx]){
                                fit = false;
                                break;
                            }
                        }
                    }
                    
                    if (fit) {
                        map << f;
                        counts.lData[f] ++;
                        
                    } else {
                        f = -1;
                    }
                }
            }
            if (f==-1) { // fit failed or unique site
                if (!sameScore) {
                    sameScore = (_SimpleList*)checkPointer(new _SimpleList);
                    sequenceHashes.Insert ((BaseRef)sequenceHash,(long)sameScore,false);
                }
                
                (*sameScore) << indices.lLength;
                map     << indices.lLength;
                indices << sequenceIndex;
                counts  << 1;
            }
        }
        
        }
    else{
        long             vd  = GetDimension(true);
        
        _Parameter      *translatedVector = (_Parameter*)checkPointer(new _Parameter [vd]),
        *translatedVector2= (_Parameter*)checkPointer(new _Parameter [vd]);
        
        _String         state1 (unit,false),
        state2 (unit,false);
        
        for (long sequenceIndex = 0; sequenceIndex < seqs; sequenceIndex++) {
            bool checkState = false;
            for (long idx=0; idx<indices.countitems(); idx++) {
                for (long m=0; m<sites; m++) {
                    RetrieveState (m,sequenceIndex, state1,false);
                    RetrieveState (m,indices.lData[idx], state2,false);
                    
                    checkState = true;
                    long idx1 = Translate2Frequencies (state1, translatedVector,  true),
                    idx2 = Translate2Frequencies (state2, translatedVector2, true);
                    
                    //printf ("(%ld, %ld) %ld = %ld %ld\n", sequenceIndex, indices.lData[idx], m, idx1, idx2); 
                    
                    if (idx2 >=0 && idx1 >=0) {
                        if (idx1==idx2) {
                            continue;
                        } else {
                            checkState = false;
                            break;
                        }
                    } else {
                        
                        // check for equal ambigs 
                        long k = 0;
                        for (; k < vd; k++){
                            if (translatedVector[k] != translatedVector2[k]){
                                break;
                            }
                        }
                            
                        if (k == vd)
                            continue;
                            
                        if (mode == 1){
                            
                                long count1 = 0,
                                     count2 = 0;
                                
                                for (long t = 0; t<vd; t++) {
                                    count1 += translatedVector[t]>0.0;
                                    count2 += translatedVector2[t]>0.0;
                                }
                                
                                if (count1 < vd && count2 < vd) {
                                    checkState = false;
                                    break;
                                }
                                
                            } else {
                                bool first  = mode==2,
                                second = mode==2;
                                if (mode == 2){
                                    for (long t = 0; (first||second)&&(t<vd); t++) {
                                        if (translatedVector[t]>0.0) {
                                            second &= (translatedVector2[t]>0.0);
                                        }
                                        if (translatedVector2[t]>0.0) {
                                            first  &= (translatedVector[t]>0.0);
                                        }
                                    }
                                    if (!(first||second)) {
                                        checkState = false;
                                        break;
                                    }
                                } else {
                                    for (long t = 0; t<vd; t++) {
                                        if (translatedVector[t]>0.0) {
                                            second |= (translatedVector2[t]>0.0);
                                        }
                                        if (translatedVector2[t]>0.0) {
                                            first  |= (translatedVector[t]>0.0);
                                        }
                                    }                       
                                    if (!(first&&second)) {
                                        checkState = false;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                
                if (checkState) {
                    map << idx;
                    counts.lData[idx] ++;
                    break;
                }
            }
            
            if (!checkState){
                map     << indices.lLength;
                indices << sequenceIndex;
                counts  << 1;
            }
            
        }
        
        delete [] translatedVector;
        delete [] translatedVector2;   
    }
    
    
    return indices.lLength;
}

//_______________________________________________________________________

void    _DataSetFilter::SetFilter (_DataSet const * ds, unsigned char unit, _SimpleList& horizontalList, _SimpleList& verticalList, bool isFilteredAlready) {
    // we must step thru the underlying dataset and recompute the frequenices
    // we will store the vertical map in theMap
    // and the horizontal map in theNodeMap
    // theFrequencies will be store the new frequencies
    // theOriginalOrder is the receptacle for the original site order in the data filter

    bool            copiedSelf = false; // tag if refiltering self

    _DataSetFilter* firstOne = nil;
    if (isFilteredAlready) {
        if ((Ptr)this == (Ptr)ds) {
            firstOne = (_DataSetFilter*)makeDynamic();
            copiedSelf = true;
        } else {
            firstOne = (_DataSetFilter*)ds;
        }
        ds       = firstOne->theData;
    }

    theMap.Clear();
    theNodeMap.Clear();
    theOriginalOrder.Clear();
    theFrequencies.Clear();
    theExclusions.Clear();
    conversionCache.Clear();
    duplicateMap.Clear();

    theData     = (_DataSet*)ds;
    unitLength  = unit;

    long        i,
                j;

     // security checks
    if (horizontalList.empty() || verticalList.lLength<unit) {
        ReportWarning (_String("Row and/or column partition is emptyString. All the data will be used by default."));
        if (horizontalList.empty()) {
            horizontalList.Populate (isFilteredAlready ? firstOne->theNodeMap.lLength : ds->NoOfSpecies(),0,1);
        }
        if (verticalList.lLength<unit) {
            verticalList.Clear();
            verticalList.Populate (isFilteredAlready ? firstOne->GetSiteCount() : ds->GetNoTypes(),0,1);
        }
    }

    if (!isFilteredAlready) {
        theNodeMap.Clear();
        theNodeMap.Duplicate (&horizontalList);
    } else {
        for (unsigned long k = 0UL; k<horizontalList.lLength; k++) {
            theNodeMap<<firstOne->theNodeMap.lData[horizontalList.lData[k]];
        }

        horizontalList.Clear();
        horizontalList.Duplicate(&verticalList);
        verticalList.Clear();
        verticalList.RequestSpace(firstOne->theOriginalOrder.lLength);

        for (i = 0; i<horizontalList.lLength; i++) {
            j = horizontalList.lData[i];
            if (j>=0 && j<firstOne->theOriginalOrder.lLength) {
                verticalList<<firstOne->theOriginalOrder.lData[j];
            } else {
                _String tooBig (j);
                if (j<0) {
                    ReportWarning  (tooBig &" is a negative site index and is ignored");
                } else {
                    ReportWarning  (tooBig &" exceeds the number of sites in the underlying data filter and is ignored");
                }
            }
        }
    }

    j = ds->NoOfSpecies();

    for (i=0; i<theNodeMap.lLength; i++) {
        if (theNodeMap.lData[i]>=j) {
            _String invalid(theNodeMap.lData[i]);
            ReportWarning ((invalid&" exceeds the number of species in the underlying dataset and is ignored"));
            theNodeMap.Delete(i);
            i--;
        }
    }

    j = ds->GetNoTypes();
    for (i=0; i<verticalList.lLength; i++)
        if (verticalList.lData[i]>=j) {
            _String invalid(verticalList.lData[i]);
            ReportWarning ((invalid&" exceeds the number of sites in the underlying dataset and is ignored"));
            verticalList.Delete(i);
            i--;
        }

    if (verticalList.lLength%unit!=0) {
        ReportWarning (_String("Number of sites in datasetfilter is not divisible by the unit - will truncate to the nearest integer"));
        while(verticalList.lLength%unit) {
            verticalList.Delete(verticalList.lLength-1);
        }
    }


    theOriginalOrder.Duplicate (&verticalList);


    // done with security checks

    _SimpleList indices;        // numeric indices intended to facilitate the reindexing
    _AVLListXL  siteIndices     (&indices);


    // sweep through the columns left to right

    duplicateMap.RequestSpace (verticalList.lLength/unit+1);

    _String      siteHolder   (unit*theNodeMap.lLength,false);

    //bool       startD = false;

    for (i=0; i<verticalList.lLength; i+=unit) {
        long colIndex = 0;

        for (j=0; j<unit; j++) // sweep within one block
            for (long k=0; k<theNodeMap.lLength; k++) // sweep down the columns
                //colIndex+=
                //(((_String**)ds->lData)[ds->theMap.lData[verticalList.lData[i+j]]])->sData[theNodeMap.lData[k]];
            {
                siteHolder[colIndex++] = (((_String**)ds->lData)[ds->theMap.lData[verticalList.lData[i+j]]])->sData[theNodeMap.lData[k]];
            }

        colIndex = siteHolder.Adler32();

        long        f = siteIndices.Find ((BaseRef)colIndex);
        _SimpleList * sameScore = nil;
        
        if (f>=0) {
            sameScore = (_SimpleList*)siteIndices.GetXtra (f);
            for (long k = 0; k<sameScore->lLength; k++) {
                bool fit = true;
                f = sameScore->lData[k];
                for (long j=0; fit&&(j<unit); j++) { // sweep within one block
                    _Site * site1 = ds->GetSite(verticalList.lData[i+j]),
                          * site2 = ds->GetSite(theMap.lData[unit*f+j]); 
                            
                    for (long k=0; k<theNodeMap.lLength; k++) // sweep down the columns
                        if (site1->sData[theNodeMap.lData[k]]!=site2->sData[theNodeMap.lData[k]]) {
                            fit = false;
                            break;
                        }
                }

                if (fit) {
                    theFrequencies[f]++;
                    duplicateMap<<f;
                    f = 0;
                    break;
                } else {
                    f = -1;
                }
            }
        }
        if (f==-1) { // fit failed or unique site
            if (!sameScore) {
                sameScore = (_SimpleList*)checkPointer(new _SimpleList);
                siteIndices.Insert ((BaseRef)colIndex,(long)sameScore,false);
             }

            (*sameScore) << theFrequencies.lLength;
            duplicateMap<<theFrequencies.lLength;
            theFrequencies<<1;
            for (j=0; j<unit; j++) {
                theMap<<verticalList.lData[i+j];
            }
        }
    }

    siteIndices.Clear();

    duplicateMap.TrimMemory();
    theOriginalOrder.TrimMemory();

    if (copiedSelf) {
        DeleteObject (firstOne);
    }

    SetDimensions();
    FilterDeletions();

}
//_______________________________________________________________________
long    _DataSetFilter::FindSpeciesName (_List& s, _SimpleList& r) const {
    // MOD 12/16/03
    r.Clear();

    _List           newNames;
    _AVLListX       matched (&newNames);

    for (unsigned long k=0UL; k<theNodeMap.lLength; k++) {
        long i = theNodeMap.lData[k];
        _String * uC = new _String (*(_String*)theData->theNames (i));
        uC->UpCase();
        matched.Insert (uC,i);
    }

    for (unsigned long m = 0UL; m < s.lLength; m++) {
        _String ts (*(_String*)s.GetItem (m));
        ts.UpCase();
        long f = matched.Find (&ts);
        if (f>=0L) {
            r << matched.GetXtra (f);
        } else {
            break;
        }
    }
    
    return r.lLength;
}

//_______________________________________________________________________
extern _String skipOmissions;

void    _DataSetFilter::FilterDeletions(_SimpleList *theExc)
{
    _Parameter      skipo;
    checkParameter (skipOmissions,skipo,0.0);

    if (skipo>.5 || theExc ) { // delete omissions
        //build up the list of "bad" sites
        _SimpleList sitesWithDeletions;
        if (!theExc) {
            for (long i=0; i<theFrequencies.lLength; i++)
                if (HasDeletions(i)) {
                    sitesWithDeletions<<i;
                }
        } else {
            _Parameter   *store_vec = (_Parameter*)checkPointer(new _Parameter [GetDimension(false)]);

            for (long i=0; i<theFrequencies.lLength; i++) {
                long pos = HasExclusions(i,theExc,store_vec);
                if  (pos != -1) {
                    sitesWithDeletions<<i;
                    _String warnMsg ((*this)(i,pos));
                    warnMsg = warnMsg & " was encountered in sequence "& *GetSequenceName (pos) & " at site pattern " & i
                              & ". All corresponding alignment columns will be removed from subsequent analyses.";
                    ReportWarning (warnMsg);
                }
            }

            delete [] store_vec;
        }

        if (sitesWithDeletions.lLength==theFrequencies.lLength) {
            _String errMsg ("All the sites in the datafilter have deletions and removing them creates an emptyString filter");
            ReportWarning(errMsg);
        }

        _SimpleList allDeleted,
                    dupDeletes;

        for (long k=0; k < duplicateMap.lLength; k++)
            if (sitesWithDeletions.BinaryFind (duplicateMap.lData[k]) >= 0) {
                dupDeletes << k;
                for (long j = 0; j < unitLength; j++ ) {
                    allDeleted << k*unitLength + j;
                }
            }

        duplicateMap.DeleteList (dupDeletes);
        dupDeletes.Clear();
        theOriginalOrder.DeleteList (allDeleted);
        theFrequencies.DeleteList (sitesWithDeletions);


        for (long i=0; i<sitesWithDeletions.lLength; i++) {
            long sitePos = sitesWithDeletions.lData[i];

            for (long j=0; j<unitLength; j++) {
                theMap.lData[sitePos*unitLength+j]=-1;
                dupDeletes << sitePos*unitLength+j;
            }
        }


        if (allDeleted.lLength) {
            /*allDeleted.Sort();*/

            _String     warnMsg ("The following sites are being omitted:"),
                        *s = (_String*)allDeleted.toStr();

            if (!theExc) {
                warnMsg = warnMsg & "(b/c of deletions/omissions)";
            }

            warnMsg = warnMsg&*s;
            DeleteObject(s);
            ReportWarning(warnMsg);

            _SimpleList shiftIdxBy (sitesWithDeletions.lLength+theFrequencies.lLength);

            long        shiftBy = sitesWithDeletions.lLength,
                        marker  = sitesWithDeletions.lData[sitesWithDeletions.lLength-1],
                        markerI = sitesWithDeletions.lLength-2;

            shiftIdxBy.lLength = sitesWithDeletions.lLength+theFrequencies.lLength;

            for (long i=shiftIdxBy.lLength-1; i>=0; i--) {
                if (i==marker) {
                    shiftBy--;
                    if (markerI>=0) {
                        marker = sitesWithDeletions.lData[markerI];
                        markerI --;
                    } else {
                        marker = -1;
                    }
                }
                shiftIdxBy.lData[i] = shiftBy;
            }
            {
                for (long i=0; i<duplicateMap.lLength; i++) {
                    duplicateMap.lData[i] -= shiftIdxBy.lData[duplicateMap.lData[i]];
                }
            }
        }

        // one final pass on theMap to clear it out
        /*for (long i=theMap.lLength-1;i>=0;i--)
            if (theMap(i)<0)
                theMap.Delete(i);*/
        _SimpleList saveMap (theMap);
        theMap.DeleteList (dupDeletes);
        {
            for (long k=0; k<theMap.lLength; k++)
                if (theMap.lData[k] < 0) {
                    saveMap.DeleteList (dupDeletes);
                    WarnError ("Internal Error in _DataSetFilter::FilterDeletions");
                }
        }
    }

}
//_______________________________________________________________________
_DataSetFilter*  _DataSetFilter::PairFilter (long index1, long index2, _DataSetFilter* result)
{
    _SimpleList species;
    species<<theNodeMap(index1);
    species<<theNodeMap(index2);
    result->SetFilter (theData,unitLength,species,theMap);
    if (theExclusions.lLength) {
        _String* s = (_String*)theExclusions.toStr();
        *s = s->Cut (1,s->Length()-2);
        result->SetExclusions   (s);
        DeleteObject(s);
    }
    return result;
}

//_________________________________________________________

void    _DataSetFilter::MatchStartNEnd (_SimpleList& order, _SimpleList& positions, _SimpleList* parent) const {
    if (order.lLength == 0) {
        return;
    }

    long p0 = order.lData[0];

    _Parameter uth;
    checkParameter (useTraversalHeuristic,uth,1.0);

    if (uth>.5) {
        if (parent)
            for (long i = 1; i < order.lLength; i++) {
                unsigned long
                j       = 0,
                n       = theNodeMap.lLength-1,
                p0  = parent->lData[i],
                p1  = order.lData[i];

                while (CompareTwoSites(p0,p1,j)) {
                    j++;
                }
                while (CompareTwoSites(p0,p1,n)) {
                    n--;
                }
                n = (n<<16) + j;
                positions << n;
            }
        else
            for (long i = 1; i < order.lLength; i++) {
                unsigned long j = 0,
                              n = theNodeMap.lLength-1,
                              p1 = order.lData[i];

                while (CompareTwoSites(p0,p1,j)) {
                    j++;
                }
                while (CompareTwoSites(p0,p1,n)) {
                    n--;
                }
                n = (n<<16) + j;
                positions << n;
                p0 = p1;
            }
    } else
        for (long i = 1; i < order.lLength; i++) {
            unsigned long j = 0,
                          n = theNodeMap.lLength-1;

            n = (n<<16) + j;
            positions << n;
        }

}

//_______________________________________________________________________

void    _DataSetFilter::SetExclusions (_String* theList, bool filter)
{

    theExclusions.Clear();
    theList->StripQuotes();

    if (theList->sLength == 0) {
        return;
    }

    _List        tokens (theList->Tokenize(','));
    _SimpleList  holder;
    _AVLList     exclusions (&holder);

    for (long k = 0; k < tokens.lLength; k++) {
      
        _String* kth_token = (_String*)tokens.GetItem(k);
      
        long posMarker = MapStringToCharIndex(*kth_token);

        if (posMarker < 0) {
            ReportWarning (_String("Exclusion request for '") & *kth_token &"' does not represent a unique state and will therefore be ignored.");
        } else {
            if (exclusions.Insert((BaseRef)posMarker) < 0) {
                ReportWarning (_String("Exclusion symbol for '") & *kth_token &"' is included more than once.");
            }
        }
    }

    exclusions.ReorderList();

    if (filter) {
        FilterDeletions (&holder);
    }

    theExclusions<<holder;
}

//_______________________________________________________________________

_String*    _DataSetFilter::GetExclusions (void) const {
    _String * res = new _String (16L, true);
  
    if (theExclusions.lLength) {
        for (long k=0; k<theExclusions.lLength-1; k++) {
            (*res) << ConvertCodeToLetters (theExclusions.lData[k], unitLength);
            (*res) << ',';
        }

        (*res) << ConvertCodeToLetters (theExclusions.lData[theExclusions.lLength-1], unitLength);
    }

    res->Finalize();

    return res;
}

//_______________________________________________________________________

unsigned long    _DataSetFilter::GetDimension (bool correct) const {
    unsigned long result = ComputePower (theData->theTT->baseLength, unitLength);
    return correct ? result - theExclusions.lLength : result;
}

//____________________________________________________________________________________
//  20110610: SLKP, some cleanup and refactoring

void    _DataSet::ProcessPartition (_String const & input2 , _SimpleList & target , bool isVertical, _SimpleList const* additionalFilter, _SimpleList const* otherDimension, _String const* scope) const {
    if (!input2.sLength) {
        return;
    }
    // decide if the input is an enumeration or a formula
    long totalLength;

    if (additionalFilter) {
        totalLength = additionalFilter->lLength;
    } else {
        totalLength = isVertical?theMap.lLength:noOfSpecies;
    }

    _String input (input2);

    if (!input.IsALiteralArgument(true)) { // not a literal argument
        _Formula fmla, lhs;
        _FormulaParsingContext fpc;
        fpc.setScope (scope);
      
        long     outcome = Parse (&fmla, input, fpc,&lhs);

        if (outcome!=HY_FORMULA_EXPRESSION) {
            WarnError (input & _String(" is an invalid partition specification"));
            return;
        }
        _PMathObj   fV = fmla.Compute();
        if (fV && fV->ObjectClass()==STRING) {
            _String newSpec (128L, true);
            newSpec << '"'
                    << ((_FString*)fV)->theString
                    << '"';
            newSpec.Finalize();
            ProcessPartition (newSpec, target, isVertical, additionalFilter, nil, scope);
        } else {
            _DataSet::MatchIndices (fmla, target, isVertical, totalLength, scope);
        }
    } else { // an explicit enumeration or a regular expression
        if (input.getChar(0)=='/' && input.getChar(input.sLength-1)=='/')
            // a regular expression
        {
            input.Trim(1,input.sLength-2);
            int   errCode;
            Ptr   regex = PrepRegExp (&input, errCode, true);
            if (errCode) {
                WarnError(GetRegExpError(errCode));
                return;
            }
            // now set do the matching
            // using only the sites that are specced in the additionalFilter

            if (!isVertical) {

                _SimpleList matches;
                for (long specCount = 0L; specCount < (additionalFilter ? additionalFilter->lLength : totalLength); specCount++) {
                    _String pattern (theMap.lLength, false);
                    long    seqPos = additionalFilter ? additionalFilter->Element (specCount) : specCount;

                    if (otherDimension)
                        for (long seqSlider = 0L; seqSlider < otherDimension->lLength; seqSlider ++) {
                            pattern.sData[seqSlider] =  GetSite(otherDimension->Element(seqSlider))->getChar (seqPos);
                        }
                    else
                        for (long seqSlider = 0L; seqSlider < theMap.lLength; seqSlider ++) {
                            pattern.sData[seqSlider] =  GetSite(seqSlider)->getChar(seqPos);
                        }

                    matches.Clear();
                    pattern.RegExpMatch (regex, matches);
                    if (matches.lLength) {
                        target << specCount;
                    }
                }
            } else {
                bool         *eligibleMarks = new bool[lLength];
                checkPointer (eligibleMarks);

                for (long fillerID = 0; fillerID < lLength; fillerID++) {
                    eligibleMarks [fillerID] = false;
                }

                if (additionalFilter)
                    for (long siteIndex = 0; siteIndex < additionalFilter->lLength; siteIndex ++) {
                        eligibleMarks[theMap.lData[additionalFilter->lData[siteIndex]]] = true;
                    }
                else
                    for (long siteIndex = 0; siteIndex < lLength; siteIndex ++) {
                        eligibleMarks[siteIndex] = true;
                    }

                _SimpleList matches;
                _String     *tempString = nil;
                if (otherDimension) {
                    tempString = new _String (otherDimension->lLength,false);
                }

                for (long siteCounter = 0; siteCounter < lLength; siteCounter ++)
                    if (eligibleMarks[siteCounter]) {
                        matches.Clear();
                        if (otherDimension) {
                            _Site * aSite = ((_Site**)lData)[siteCounter];
                            for (long tc = 0; tc < otherDimension->lLength; tc++) {
                                tempString->sData[tc] = aSite->sData[otherDimension->lData[tc]];
                            }
                            tempString->RegExpMatch (regex, matches);
                        } else {
                            ((_Site**)lData)[siteCounter]->RegExpMatch (regex, matches);
                        }
                        if (matches.lLength == 0) {
                            eligibleMarks[siteCounter] = false;
                        }
                    }

                DeleteObject (tempString);

                if (additionalFilter) {
                    for (long afi = 0; afi < additionalFilter->lLength; afi++)
                        if (eligibleMarks[theMap.lData[additionalFilter->lData[afi]]]) {
                            target << afi;
                        }
                } else {
                    for (long afi = 0; afi < theMap.lLength; afi++)
                        if (eligibleMarks[theMap.lData[afi]]) {
                            target << afi;
                        }
                }
                delete eligibleMarks;
            }
            FlushRegExp (regex);
        } else {
            input.KillSpaces (input);
            // now process the string
            long count = 0,anchor,k;

            _SimpleList numbers,
                        links;

            numbers.RequestSpace (1024);
            links.RequestSpace (1024);

            // first check if it is has a comb filter
            if ((input.sData[0]=='<')&&(input.sData[input.sLength-1]=='>')) {
                for (count=1; count<input.sLength-1; count++) {
                    if (input.sData[count]!='0') {
                        numbers<<count-1;
                    }
                }
                if (numbers.lLength) {
                    k = input.sLength-2; // step size
                    anchor = 0;
                    if (totalLength == -1) {
                        totalLength = theMap.lLength;
                    }
                    while (anchor<totalLength-k) {
                        for (count = 0; count< numbers.lLength; count++) {
                            target<<anchor+numbers.lData[count];
                        }
                        anchor+=k;
                    }
                    if ( (k=totalLength-1-anchor) ) {
                        for (count = 0; count< numbers.lLength; count++) {
                            if (numbers.lData[count]>k) {
                                break;
                            }
                            target<<anchor+numbers.lData[count];
                        }
                    }
                    return;

                }
            }

            while (count<input.sLength) {
                anchor = count;
                for (; (count<input.sLength)&&(isdigit(input[count])); count++) ;
                long    aNumber = (input.Cut (anchor,count-1)).toNum();
                if (aNumber < 0) {
                    _String warnMsg ("A negative number was found in partition specification: ");
                    ReportWarning (warnMsg & input.Cut (0,anchor-1) & '?' & input.Cut (anchor,-1));
                    target.Clear();
                    return;
                }
                numbers<< aNumber;

                if ((input[count]=='<')||(input[count]=='>')) {
                    _String warnMsg ("A comb partition cannot be combined with other types. The entire partition is reset to first..last");
                    ReportWarning (warnMsg & input.Cut (0,anchor-1) & '?' & input.Cut (anchor,-1));
                    target.Clear();
                    return;
                }
                if (input[count]=='&') {
                    links << numbers.lLength;
                }
                if ((input[count]==',')||(count==input.sLength)) { // wrap it up dude
                    if (numbers.lLength==1) {
                        target<<numbers(0);
                    } else {
                        if (links.lLength==0) {
                            if (numbers[0]>numbers[1]) { // backward order
                                for (k = numbers[0]; k>=numbers[1]; k--) {
                                    target<<k;
                                }
                            } else {
                                for (k = numbers[0]; k<=numbers[1]; k++) {
                                    target<<k;
                                }
                            }
                        } else {
                            // linked locations
                            if (links.lLength!=(numbers.lLength-2)/2) {
                                _String errMsg ("A part of the partition specification has not been understood and has been skipped.");
                                ReportWarning (errMsg);
                                target.Clear();
                                return;
                            } else {
                                _SimpleList signs;
                                signs<<(numbers(0)<numbers(1)?1:-1);
                                for (k = 0; k<links.lLength; k++) {
                                    signs<<(numbers(links(k))<numbers(links(k)+1)?1:-1);
                                }

                                long l,m;
                                for (k=numbers(0), l=0 ; signs(0)*k<=signs(0)*numbers(1); k+=signs(0), l++) {
                                    target<<numbers(0)+l*signs(0);
                                    for (m=0; m<links.lLength; m++) {
                                        target<<numbers(links(m))+l*signs(m+1);
                                    }
                                }

                            }
                        }
                    }
                    numbers.Clear();
                    links.Clear();
                }
                count++;
            }
        }
    }
}
//_______________________________________________________________________

void     _DataSetFilter::SetMap  (_String&s)
{
    theNodeMap.Clear();
    if (s.Length()) {
        long f,g=0;
        //_String sc(",");
        f = s.Find(',');
        while (f!=-1) {
            theNodeMap<<(long)s.Cut(g,f-1).toNum();
            g = f+1;
            f = s.Find (',',f+1,-1);
        }
        theNodeMap<<(long)s.Cut(g,-1).toNum();
    }
}

//_________________________________________________________

void    _DataSetFilter::FindAllSitesLikeThisOne (long index, _SimpleList& receptacle)
{
    long   oindex = theOriginalOrder.Find(index),m;

    if (oindex<0) {
        return;
    }

    if (theData->NoOfSpecies()==theNodeMap.lLength) {
        long *matchMap = new long[unitLength];

        for (m=0; m<unitLength; m++)
            //matchMap[m] = theData->theMap.lData[theOriginalOrder.lData[oindex+1]];
        {
            matchMap[m] = theData->theMap.lData[theOriginalOrder.lData[oindex+m]];
        }


        for (long k=0; k<theOriginalOrder.lLength; k+=unitLength) {
            for (m=0; m<unitLength; m++) {
                if (theData->theMap.lData[theOriginalOrder.lData[k+m]]!=matchMap[m]) {
                    break;
                }
            }
            if (m==unitLength)
                for (m=0; m<unitLength; m++) {
                    receptacle<<theOriginalOrder.lData[k+m];
                }
        }

        delete [] matchMap;
    } else {
        char ** matchMap = (char**)MemAllocate (sizeof (char*) * unitLength);
        checkPointer (matchMap);

        for (m=0; m<unitLength; m++) {
            matchMap[m] = ((_Site*)(((BaseRef*)theData->lData)[theData->theMap.lData[oindex+m]]))->sData;
        }
        for (long k=0; k<theOriginalOrder.lLength; k+=unitLength) {
            for (m=0; m<unitLength; m++) {
                char* checkStr = ((_Site*)(((BaseRef*)theData->lData)[theData->theMap.lData[k+m]]))->sData;
                long t;
                for (t = 0; t<theNodeMap.lLength; t++) {
                    if (checkStr[t]!=matchMap[m][t]) {
                        break;
                    }
                }
                if (t<theNodeMap.lLength) {
                    break;
                }
            }
            if (m==unitLength)
                for (m=0; m<unitLength; m++) {
                    receptacle<<theOriginalOrder.lData[k+m];
                }
        }
        delete matchMap;
    }
}

//_______________________________________________________________________
_String* _DataSetFilter::MakeSiteBuffer (void) const {
  return new _String ((unsigned long)unitLength, false);
}

//_______________________________________________________________________

_String&     _DataSetFilter::operator () (unsigned long site, unsigned long pos) {
    if (!accessCache || accessCache->sLength != unitLength) {
        if (accessCache) {
            DeleteObject (accessCache);
        }
        accessCache = new _String ((unsigned long)unitLength, false);
    }

    long vIndex = theNodeMap.lData[pos];
    if (unitLength==1) {
        accessCache->sData[0]=(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site]]])->sData[vIndex];
    } else {
        site*=unitLength;
        for (int k = 0; k<unitLength; k++) {
            accessCache->sData[k] = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site++]]])->sData[vIndex];
        }
    }
    return *accessCache;
}

//_______________________________________________________________________

const _String     _DataSetFilter::RetrieveState (unsigned long site, unsigned long pos) const {
  _String state ((unsigned long)unitLength, false);
  RetrieveState (site, pos, state, false);
  return state;
  
}

//_______________________________________________________________________

void     _DataSetFilter::RetrieveState (unsigned long site, unsigned long pos, _String& reply, bool map) const
{
    long vIndex = theNodeMap.lData[pos];
    if (map) {
        if (unitLength==1) {
            reply.sData[0]=(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[duplicateMap.lData[site]]]])->sData[vIndex];
        } else {
            site = unitLength*duplicateMap.lData[site];
            for (int k = 0; k<unitLength; k++) {
                reply.sData[k] = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site++]]])->sData[vIndex];
            }
        }
    } else {
        if (unitLength==1) {
            reply.sData[0]=(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site]]])->sData[vIndex];
        } else
            site*=unitLength;
            for (int k = 0; k<unitLength; k++) {
                reply.sData[k] = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site++]]])->sData[vIndex];
            }
    }
}

//_______________________________________________________________________

void _DataSetFilter::GrabSite (unsigned long site, unsigned long pos, _String& storage)
{

    long vIndex = theNodeMap.lData[pos];
    if (unitLength==1) {
        storage.sData[0]=(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site]]])->sData[vIndex];
    } else {
        site*=unitLength;
        for (int k = 0; k<unitLength; k++) {
            storage.sData[k] = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site++]]])->sData[vIndex];
        }
    }
}

//_______________________________________________________________________

void _DataSetFilter::GrabSite (unsigned long site, unsigned long pos, char * s)
{
    long vIndex = theNodeMap.lData[pos];
    if (unitLength==1) {
        s[0]=(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site]]])->sData[vIndex];
    } else {
        site*=unitLength;
        for (int k = 0; k<unitLength; k++) {
            s[k] = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site++]]])->sData[vIndex];
        }
    }
}

//_______________________________________________________________________

_SimpleList* _DataSetFilter::CountAndResolve (long pattern, _Parameter * storage, bool randomly)
// last cell in the list contains the count of distinct characters in the column
{
    _SimpleList* resList = new _SimpleList (theNodeMap.lLength+1,0,0),
    counts (dimension,0,0);

    checkPointer (resList);

    _List        ambStates;
    _String      aState  (unitLength, false);

    _Parameter*  freqStorage = storage;

    if (!freqStorage) {
        freqStorage = new _Parameter [undimension];
    }

    long    normalizingSum = 0,
            charCount      = 0;

    for (long k=0; k<theNodeMap.lLength; k++) {
        GrabSite (pattern, k, aState);
        long      characterRes = Translate2Frequencies (aState, freqStorage, true);
        if (characterRes>=0) {
            resList->lData[k] = characterRes;

            if (characterRes >= dimension) {
                WarnError (_String("Internal error in _DataSetFilter::CountAndResolve\n"));
            }

            if ((counts.lData[characterRes]++) == 0) {
                normalizingSum ++;
            }

            charCount ++;
        } else {
            _SimpleList * possibleResolutions = new _SimpleList;
            if (!possibleResolutions) {
                checkPointer (possibleResolutions);
            }

            (*possibleResolutions) << k;

            for (long m=0; m<dimension; m++)
                if (freqStorage[m]>0.) {
                    (*possibleResolutions) << m;
                }

            ambStates.AppendNewInstance (possibleResolutions);
        }
    }

    if (normalizingSum > 0) {
        if (ambStates.lLength) {
            _SimpleList  ambResolutions (dimension,0,0);
            for (long t=0; t<ambStates.lLength; t++) {
                _SimpleList * stateResolutions = (_SimpleList*)ambStates(t);

                if (!randomly) {
                    long          totalSum = 0,
                                  idx = 0;

                    for (long l=1; l<stateResolutions->lLength; l++) {
                        long tmp = counts.lData[stateResolutions->lData[l]];
                        if (tmp>totalSum) {
                            idx = l;
                            totalSum = tmp;
                        }
                    }
                    if (idx > 0)
                        // if no resolutions, resolve randomly
                    {
                        idx = stateResolutions->lData[idx];
                        resList->lData[stateResolutions->lData[0]] = idx;
                        ambResolutions.lData [idx] ++;
                        continue;
                    }

                }

                long          totalSum = 0;
                for (long l=1; l<stateResolutions->lLength; l++) {
                    totalSum += counts.lData[stateResolutions->lData[l]];
                }

                if (totalSum > 0) {
                    long          randomN = genrand_real2() * totalSum - counts.lData[stateResolutions->lData[1]],
                                  ind = 1;

                    while (randomN > 0) {
                        randomN -= counts.lData[stateResolutions->lData[++ind]];
                    }

                    totalSum = stateResolutions->lData[ind];
                } else {
                    long          randomN = genrand_real2() * charCount - counts.lData[0],
                                  ind = 0;

                    while (randomN > 0) {
                        randomN -= counts.lData[++ind];
                    }
                }
                resList->lData[stateResolutions->lData[0]] = totalSum;
                ambResolutions.lData [totalSum] ++;
            }

            for (long l=0; l<dimension; l++)
                if (ambResolutions.lData[l] && !counts.lData[l]) {
                    normalizingSum ++;
                }
        }
    }

    resList->lData[theNodeMap.lLength] = normalizingSum;

    if (freqStorage != storage) {
        delete freqStorage;
    }

    return       resList;
}

//_______________________________________________________________________

_Matrix* _DataSetFilter::PairwiseCompare (_SimpleList* s1, _SimpleList *s2, _List* labels)
// s1 and s2 are the lists produced by CountAndResolve
// if labels is not nil, then it will receive row and column labels in the contigency table
// the result matrix has rows labeled by states in s1, and columns - by states in s2
{
    long    * sort1 = new long[dimension],
    * sort2 = new long[dimension],
    c = s2->lData[s2->lLength-1];

    _Matrix * res   = new _Matrix (s1->lData[s1->lLength-1],c,false,true);

    if (sort1 && sort2 && res) {
        for (long k = 0; k<dimension; k++) {
            sort1[k] = -1;
            sort2[k] = -1;
        }

        long idx1 = 0,
             idx2 = 0;

        _SimpleList  *lbl1 = nil,
                      *lbl2 = nil;

        if (labels) {
            lbl1 = new _SimpleList;
            lbl2 = new _SimpleList;

            checkPointer (lbl1);
            checkPointer (lbl2);

            (*labels) << lbl1;
            (*labels) << lbl2;

            DeleteObject (lbl1);
            DeleteObject (lbl2);
        }

        for (long k2 = 0; k2 < s1->lLength-1; k2++) {
            long c1 = s1->lData[k2],
                 c2 = s2->lData[k2];

            if (sort1[c1] < 0) {
                sort1[c1] = idx1;
                if (lbl1) {
                    (*lbl1) << c1;
                }
                c1 = idx1++;
            } else {
                c1 = sort1[c1];
            }

            if (sort2[c2] < 0) {
                sort2[c2] = idx2;
                if (lbl2) {
                    (*lbl2) << c2;
                }
                c2 = idx2++;
            } else {
                c2 = sort2[c2];
            }

            /*if ((c1>=res->GetHDim())||(c2>=res->GetVDim()))
            {
                printf ("\nInternal Error\n");
            }*/

            res->theData[c1*c+c2] += 1.;
        }

        delete [] sort1;
        delete [] sort2;
    } else {
        checkPointer (nil);
    }

    return res;
}

//_______________________________________________________________________

_List *  _DataSetFilter::ComputePatternToSiteMap (void) const {
    _List * result = new _List ();
  
    for (unsigned long k = 0UL; k < theFrequencies.lLength; k++) {
        (*result) < new _SimpleList;
    }
    for (unsigned long s = 0UL; s < duplicateMap.lLength; s++) {
        *((_SimpleList**)result->lData)[duplicateMap.lData[s]] << s;
    }
    return result;
}

//_______________________________________________________________________

char     _DataSetFilter::GetChar (unsigned long site, unsigned long pos)
{
    //long vIndex = theNodeMap.lLength?theNodeMap.lData[pos]:pos;
    return (*theData)(theMap.lData[site],theNodeMap.lData[pos],1);
}
//_______________________________________________________________________

bool     _DataSetFilter::CompareTwoSites (unsigned long site1, unsigned long site2, unsigned long pos1) const {
    pos1 = theNodeMap.lData[pos1];

    if (unitLength == 3) { // codon
        site1*=3;
        site2*=3;
        return
            ((((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1]]])->sData[pos1]==
             ( ((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site2]]])->sData[pos1])
            &&((((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1+1]]])->sData[pos1]==
               (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site2+1]]])->sData[pos1])
            &&((((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1+2]]])->sData[pos1]==
               (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site2+2]]])->sData[pos1]);
    } else {
        site1*=unitLength;
        site2*=unitLength;
        unsigned long k = 0UL;

        /*if ((((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1]]])->sLength<=pos1)
        {
            printf ("(%d)%s\n(%d)%s\n",site1,(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1]]])->sData,
                    site2,(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site2]]])->sData);
            FlagError ("Internal DataSetFilter bug\n");
        }*/

        for (; k<unitLength; k++) {
            if ((((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1]]])->sData[pos1]!=
                    (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site2]]])->sData[pos1]) {
                break;
            }
            site1++;
            site2++;
        }
        if (k==unitLength) {
            return true;
        }
    }
    return false;
}

//_______________________________________________________________________

bool     _DataSetFilterNumeric::CompareTwoSites (unsigned long, unsigned long, unsigned long) const {
    return false;
}


//_______________________________________________________________________

bool     _DataSetFilter::CompareTwoSitesChar (unsigned long site1, unsigned long site2, unsigned long pos1) const {

//  long *fI = theMap.quickArrayAccess();

//  if (theNodeMap.lLength)
//  {
    pos1 = theNodeMap(pos1);
//  }
//  else
//  {
//      vIndex1=pos1;
//  }
//  return ((*theData)(fI[site1],pos1, 1)==(*theData)(fI[site2],pos1, 1));
    return ((*theData)(theMap.lData[site1],pos1, 1)==(*theData)(theMap.lData[site2],pos1, 1));
}
//_______________________________________________________________________
long    _DataSetFilter::SiteFrequency (unsigned long site)
{
    return theFrequencies.lData[site];
}

//_______________________________________________________________________
bool    _DataSetFilter::HasDeletions (unsigned long site, _AVLList* storage)
{
    long        loopDim  = GetDimension();
    _Parameter* store    = new _Parameter [loopDim];

    long j,
         upTo = theNodeMap.lLength?theNodeMap.lLength:theData->NoOfSpecies();

    bool outcome = false;

    for (unsigned int k = 0; k<upTo; k++) {
        Translate2Frequencies ((*this)(site,k), store, false);

        bool oneF = false,
             zeroF = false;

        for (j=0; j<loopDim; j++) {
            if (store[j]==0.0) {
                zeroF = true;
            } else if (store[j]==1.0) {
                oneF = true;
            }
        }
        if (!(oneF&&zeroF)) {
            if (storage) {
                outcome = true;
                storage->Insert ((BaseRef)theNodeMap.lData[k]);
            } else {
                delete [] store;
                return true;
            }
        }
    }

    delete [] store;
    return outcome;
}


//_______________________________________________________________________
bool    _DataSetFilter::IsConstant (unsigned long site,bool relaxedDeletions)
{
    _Parameter *store = new _Parameter [GetDimension()],
               *store2 = new _Parameter [GetDimension()];

    unsigned long j,
         upTo = theNodeMap.lLength?theNodeMap.lLength:theData->NoOfSpecies(),
         loopDim = GetDimension();

    Translate2Frequencies ((*this)(site,0), store, false);

    if (relaxedDeletions) {
        for (unsigned long k = 1UL; k<upTo; k++) {
            Translate2Frequencies ((*this)(site,k), store2, false);
            for (j=0UL; j<loopDim; j++) {
                if (store2[j]==0.0) {
                    store[j]=0.0;
                }
            }
        }
        for (j=0UL; j<loopDim; j++)
            if (store[j]!=0.0) {
              break;
            }
      
        delete [] store;
        delete [] store2;
        return j!=loopDim;
      
    } else {
        for (unsigned long k = 1; k<upTo; k++) {
            Translate2Frequencies ((*this)(site,k), store2, false);
            for (j=0UL; j<loopDim; j++)
              
                if (store[j]!=store2[j]) {
                    delete [] store;
                    delete [] store2;
                    return false;
                }
        }
    }
  
    delete [] store;
    delete [] store2;
    return true;
}

//_______________________________________________________________________

_Matrix*        _DataSetFilter::GetFilterCharacters (bool flip) const {
    long        unit_length = GetUnitLength (),
                seq_length  = flip?theFrequencies.lLength:(GetSiteCount () / unitLength),
                seq_count   = NumberSpecies();

    _List       result;

    _String      char_buffer (unit_length ,false);

    if (flip) {
        for (long k=0; k< seq_length; k++) {
            _String *alignment_column = new _String (seq_count+1,true);
            for (long k2=0; k2< seq_count ; k2++) {
                RetrieveState(k,k2,char_buffer,false);
                (*alignment_column) << char_buffer;
            }
            alignment_column->Finalize();
            result  < alignment_column;
        }
    } else
        for (long k=0; k < seq_count; k++) {
            result      < GetSequenceCharacters(k);
        }

    return new _Matrix (result);
}

//_______________________________________________________________________

_String*        _DataSetFilter::GetSequenceCharacters (long seqID)  const{
    long            unitSizeL   = GetUnitLength();
  
    _String * aSequence = new _String (GetSiteCount(),true);

    if (seqID >= 0 && seqID < theNodeMap.lLength) {
        _String      aState (unitSizeL,false);
        unsigned long        upTo = GetSiteCountInUnits();
        for (unsigned long k2=0UL; k2<upTo; k2++) {
            RetrieveState(k2,seqID,aState);
            (*aSequence) << aState;
        }
    }
    aSequence->Finalize();
    return aSequence;
}

  //_______________________________________________________________________

_String*        _DataSet::GetSequenceCharacters (long seqID)  const{
  
  unsigned long        upTo = NoOfColumns();
  _String * aSequence = new _String (upTo,true);
  
  if (seqID >= 0 && seqID < noOfSpecies) {
    for (unsigned long k2=0UL; k2<upTo; k2++) {
      (*aSequence) << GetSite (k2)->getChar (seqID);
    }
  }
  aSequence->Finalize();
  return aSequence;
}


//_______________________________________________________________________
long    _DataSetFilter::HasExclusions (unsigned long site, _SimpleList* theExc, _Parameter*store )
{
    long   filterDim = GetDimension(false);

    if (theNodeMap.lLength)
        for (unsigned long k = 0; k<theNodeMap.lLength; k++) {
            Translate2Frequencies   ((*this)(site,k), store, false);

            long                    j                       = 0,
                                    s                      = 0;

            for (j=0; j<filterDim; j++)
                if (store[j] > 0.0) {
                    s++;
                    if (theExc->Find(j) < 0) {
                        break;
                    }
                }

            if (j == filterDim && s) {
                return k;
            }
        }

    return -1;
}
//_______________________________________________________________________
void    _DataSetFilter::Freeze (long site)
{
    for (int k = 0; k<unitLength; k++) {
        _Site* tC = (_Site*)((*(_List*)theData)(theData->theMap(this->theMap(site*unitLength+k))));
        tC->SetRefNo(-1);
        tC->PrepareToUse();
    }
}

//_______________________________________________________________________
void    _DataSetFilter::UnFreeze (long site)
{
    for (int k = 0; k<unitLength; k++) {
        _Site* tC = (_Site*)((*(_List*)theData)(theData->theMap(this->theMap(site*unitLength+k))));
        tC->SetRefNo(0);
//      tC->Archive();
    }
}

//_________________________________________________________
_Matrix* _DataSetFilter::ComputePairwiseDifferences (long i, long j, _hy_dataset_filter_ambiguity_resolution resolution_option) const {
  
  try {
    
    if (unitLength > 3) {
      throw _String("ComputePairwiseDifferences is not implemented for data filters with unit size > 3");
    }
    
    long    mxDim      = GetDimension (true);
    
    _Matrix     *res   = new _Matrix  (mxDim,mxDim,false,true);
    
    _Parameter  *sm1   = new _Parameter[mxDim],
    *sm2   = new _Parameter[mxDim];
    
    
    
    _String      state1 (unitLength,false),
    state2 (unitLength,false);
    
    
    if (conversionCache.lLength == 0) {
      throw _String ("ComputePairwiseDifferences called on a filter with emptyString conversionCache");
    }
    
    long        *tcodes  = conversionCache.lData+89,
    *ccodes  = conversionCache.lData+1,
    ccount   = conversionCache.lData[0];
    
    for (unsigned long site_pattern = 0UL; site_pattern < theFrequencies.lLength; site_pattern++) {
      long s1 = -1, s2 = -1;
      
      int c1, c2;
      
      c1 = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[unitLength*site_pattern]]])->sData[theNodeMap.lData[i]],
      c2 = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[unitLength*site_pattern]]])->sData[theNodeMap.lData[j]];
      
      if (unitLength == 1) {
        s1 = conversionCache.lData[(c1-40)*(undimension+1)+undimension],
        s2 = conversionCache.lData[(c2-40)*(undimension+1)+undimension];
      } else {
        int         c12 = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[unitLength*site_pattern+1]]])->sData[theNodeMap.lData[i]],
        c22 = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[unitLength*site_pattern+1]]])->sData[theNodeMap.lData[j]];
        
        
        state1.sData[0] = c1;
        state1.sData[1] = c12;
        
        state2.sData[0] = c2;
        state2.sData[1] = c22;
        
        c1  = ccodes[c1-40];
        c12 = ccodes[c12-40];
        
        c2  = ccodes[c2-40];
        c22 = ccodes[c22-40];
        
        if (unitLength == 2) {
          if ((c1>=0)&&(c12>=0)) {
            s1 = tcodes[c1*ccount+c12];
          }
          
          if ((c2>=0)&&(c22>=0)) {
            s2 = tcodes[c2*ccount+c22];
          }
        } else {
          int         c13 = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[unitLength*site_pattern+2]]])->sData[theNodeMap.lData[i]],
          c23 = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[unitLength*site_pattern+2]]])->sData[theNodeMap.lData[j]];
          
          //printf ("\n%c %c", c13, c23);
          
          state1.sData[2] = c13;
          state2.sData[2] = c23;
          
          c13 = ccodes[c13-40];
          c23 = ccodes[c23-40];
          
          //printf (" %d %d %s %s\n", c13, c23, state1.sData, state2.sData);
          
          if ((c1>=0)&&(c12>=0)&&(c13>=0)) {
            s1 = tcodes[ccount*(c1*ccount+c12)+c13];
          }
          
          if ((c2>=0)&&(c22>=0)&&(c23>=0)) {
            s2 = tcodes[ccount*(c2*ccount+c22)+c23];
          }
        }
      }
      /*
      if (pcAmbiguitiesAverage.Equal (resFlag)) {
        res = dsf->ComputePairwiseDifferences (seq,site,1);
      } else if (pcAmbiguitiesResolve.Equal (resFlag)) {
        res = dsf->ComputePairwiseDifferences (seq,site,2);
      } else if (pcAmbiguitiesSkip.Equal (resFlag)) {
        res = dsf->ComputePairwiseDifferences (seq,site,3);
      } else {
        res = dsf->ComputePairwiseDifferences (seq,site,0);
      }*/

      
      if (s1>=0 && s2>=0) { // one to one
        res->theData[s1*mxDim+s2] += theFrequencies.lData[site_pattern];
      } else {
        if (resolution_option != kAmbiguityHandlingSkip) {
          _Matrix * freqsAtSite = nil;
          
          if (resolution_option != kAmbiguityHandlingResolve) {
            _SimpleList   //seqList,
            siteList;
            
            
            for (long si = 0; si < unitLength; si++) {
              siteList << theMap.lData[unitLength*site_pattern+si];
            }
            
            _SimpleList copy_node_oder (theNodeMap);
            freqsAtSite     = theData->HarvestFrequencies (unitLength, unitLength, 0, copy_node_oder, siteList);
            if (theExclusions.lLength) {
              long k = 0,
                   u = GetDimension (false);
              
              for (long i = 0; i<u; i++) {
                if (i==theExclusions.lData[k] && k<theExclusions.lLength) {
                  k++;
                  continue;
                }
                freqsAtSite->theData[i-k] = freqsAtSite->theData[i];
              }
            }
            //XferwCorrection (freqsAtSite->theData, freqsAtSite->theData, mxDim);
          }
          
          if (s1>=0) {
            // one to many
            if (unitLength>1) {
              Translate2Frequencies (state2,sm1,false);
            } else {
              Translate2Frequencies (c2,sm1,false);
            }
            
            if (freqsAtSite) {
              if (resolution_option == kAmbiguityHandlingAverageFrequencyAware) {
                _Parameter totalW = 0.0;
                
                for  (long m=0; m<mxDim; m++)
                  if (sm1[m]>0.0) {
                    totalW += freqsAtSite->theData[m];
                  }
                
                if (totalW>0.0) {
                  s1 = s1*mxDim;
                  
                  for  (long m=0; m<mxDim; m++,s1++)
                    if (sm1[m]>0.0) {
                      res->theData[s1] += theFrequencies.lData[site_pattern]*freqsAtSite->theData[m]/totalW;
                    }
                }
                
              } else {
                _Parameter maxW   = 0.0;
                long       maxIdx = -1;
                
                for  (long m=0; m<mxDim; m++) {
                  if (sm1[m]>0.0) {
                    _Parameter myWeight = freqsAtSite->theData[m];
                    if (myWeight > maxW) {
                      maxW = myWeight;
                      maxIdx = m;
                    }
                  }
                }
                
                if (maxIdx>=0) {
                  res->theData[s1*mxDim+maxIdx] += theFrequencies.lData[site_pattern];
                }
              }
            } else {
              /* adopt the following convention here:
               - if ambig resolves to one s1 - count as a match
               - otherwise - count all contributions equally
               */
              
              if (sm1[s1] > 0.0) {
                res->theData[s1*mxDim+s1] += theFrequencies.lData[site_pattern];
              } else {
                long ambCount = 0;
                for  (long m=0; m<mxDim; m++) {
                  if (sm1[m]>0.0) {
                    ambCount ++;
                  }
                }
                
                s1 *= mxDim;
                
                _Parameter addFac = theFrequencies.lData[site_pattern]/(_Parameter)ambCount;
                
                for  (long m=0; m<mxDim; m++,s1++)
                  if (sm1[m]>0.0) {
                    res->theData[s1] += addFac;
                  }
              }
            }
          } else {
            if (s2>=0)
              // many to one
            {
              if (unitLength>1) {
                Translate2Frequencies (state1,sm1,false);
              } else {
                Translate2Frequencies (c1,sm1,false);
              }
              
              if (freqsAtSite) {
                if (resolution_option == kAmbiguityHandlingAverageFrequencyAware) {
                  _Parameter totalW = 0.0;
                  
                  for  (long m=0; m<mxDim; m++)
                    if (sm1[m]>0.0) {
                      totalW += freqsAtSite->theData[m];
                    }
                  
                  if (totalW>0.0) {
                    for  (long m=0; m<mxDim; m++,s2+=mxDim)
                      if (sm1[m]>0.0) {
                        res->theData[s2] += theFrequencies.lData[site_pattern]*freqsAtSite->theData[m]/totalW;
                      }
                  }
                  
                } else {
                  _Parameter maxW   = 0.0;
                  long       maxIdx = -1;
                  
                  for  (long m=0; m<mxDim; m++) {
                    if (sm1[m]>0.0) {
                      _Parameter myWeight = freqsAtSite->theData[m];
                      if (myWeight > maxW) {
                        maxW = myWeight;
                        maxIdx = m;
                      }
                    }
                  }
                  
                  if (maxIdx>=0) {
                    res->theData[maxIdx*mxDim+s2] += theFrequencies.lData[site_pattern];
                  }
                }
              } else {
                if (sm1[s2] > 0.0) {
                  res->theData[s2*mxDim+s2] += theFrequencies.lData[site_pattern];
                } else {
                  long ambCount = 0;
                  for  (long m=0; m<mxDim; m++)
                    if (sm1[m]>0.0) {
                      ambCount ++;
                    }
                  
                  _Parameter addFac = theFrequencies.lData[site_pattern]/(_Parameter)ambCount;
                  {
                    for  (long m=0; m<mxDim; m++,s2+=mxDim)
                      if (sm1[m]>0.0) {
                        res->theData[s2] += addFac;
                      }
                  }
                }
              }
            } else
              // many to many
            {
              if (unitLength>1) {
                Translate2Frequencies (state1,sm1,false);
                Translate2Frequencies (state2,sm2,false);
              } else {
                Translate2Frequencies (c1,sm1,false);
                Translate2Frequencies (c2,sm2,false);
              }
              
              if (freqsAtSite) {
                if (resolution_option == kAmbiguityHandlingAverageFrequencyAware) {
                  _Parameter totalW = 0.0;
                  
                  for  (long m=0; m<mxDim; m++)
                    if (sm1[m]>0)
                      for  (long m2=0; m2<mxDim; m2++)
                        if (sm2[m2]>0) {
                          totalW += freqsAtSite->theData[m]*freqsAtSite->theData[m2];
                        }
                  
                  if (totalW>0.0) {
                    for  (long m=0; m<mxDim; m++)
                      if (sm1[m]>0)
                        for  (long m2=0; m2<mxDim; m2++)
                          if (sm2[m2]>0) {
                            res->theData[m*mxDim+m2] += theFrequencies.lData[site_pattern]*freqsAtSite->theData[m]*freqsAtSite->theData[m2]/totalW;
                          }
                  }
                  
                } else {
                  _Parameter maxW   = 0.0;
                  long       maxIdx  = -1,
                  maxIdx2 = -1;
                  
                  for  (long m=0; m<mxDim; m++)
                    if (sm1[m]>0)
                      for  (long m2=0; m2<mxDim; m2++)
                        if (sm2[m2]>0) {
                          _Parameter myWeight = freqsAtSite->theData[m]*freqsAtSite->theData[m2];
                          if (myWeight > maxW) {
                            maxW = myWeight;
                            maxIdx  = m;
                            maxIdx2 = m2;
                          }
                        }
                  
                  if (maxIdx>=0) {
                    res->theData[maxIdx*mxDim+maxIdx2] += theFrequencies.lData[site_pattern];
                  }
                }
              } else {
                long ambCount  = 0,
                ambCount2 = 0,
                m         = 0;
                
                for  (; m<mxDim; m++) {
                  if (sm1[m]>0.0) {
                    if (sm2[m]>0.0) {
                      break;
                    } else {
                      ambCount ++;
                    }
                  } else if (sm2[m]>0.0) {
                    ambCount2 ++;
                  }
                }
                
                if (m==mxDim) {
                  _Parameter addFac = theFrequencies.lData[site_pattern]/(_Parameter)(ambCount*ambCount2);
                  
                  for  (long m=0; m<mxDim; m++)
                    if (sm1[m]>0)
                      for  (long m2=0; m2<mxDim; m2++)
                        if (sm2[m2]>0) {
                          res->theData[m*mxDim+m2] += addFac;
                        }
                }
              }
            }
          }
          DeleteObject (freqsAtSite);
        }
      }
    }
    
    delete[] sm1;
    delete[] sm2;
    
    return res;
  }
  catch (const _String error) {
    WarnError (error);
    return    new _Matrix (1,1,false,true);
    
  }
}

//_________________________________________________________

void _DataSetFilter::ComputePairwiseDifferences (_Matrix& target, long i, long j) const
// matrix of dimension nx4n containing pairwise distances as follows (n=number of species)
// first lower diag - count the same (AA,CC,GG,TT)
// first upper diag - count AC,CA
// 2nd   lower diag - count AG,GA
// 2nd   upper diag - count AT,TA
// 3rd   lower diag - count CG,GC
// 3rd   upper diag - count CT,TC
// 4th   lower diag - count GT,TG
{
    if ((target.GetHDim()!=1)||(target.GetVDim()!=7)) {
        CreateMatrix (&target,1,7,false,true,false);
    }

    if (!theData->theTT->IsStandardNucleotide()) {
        return;
    }
    long k,l;

    for (k=0; k<7; k++) {
        target[k] = 0;
    }
  
    k = theNodeMap.lData[i];
    l = theNodeMap.lData[j];

    if (l>k) {
        EXCHANGE (k,l);
    }

    for (unsigned long m=0; m < theMap.lLength; m++) {
        char const * thisSite = GetColumn (m);
        char a = thisSite[k],
             b = thisSite[l];

        long fc = theFrequencies.lData[m/unitLength];

        if (a>b) {
          EXCHANGE (a,b);
        }
      
        if (a==b) {
            target[0]+=fc;
        } else {
            if (a=='A') {
                switch (b) {
                case 'C': {
                    target[1]+=fc;
                    break;
                }
                case 'G': {
                    target[2]+=fc;
                    break;
                }
                case 'T': {
                    target[3]+=fc;
                    break;
                }
                }
            } else if (a=='C') {
                switch (b) {
                case 'G': {
                    target[4]+=fc;
                    break;
                }
                case 'T': {
                    target[5]+=fc;
                    break;
                }
                }
            } else if (a=='G') {
                if (b=='T') {
                    target[6]+=fc;
                }
            }

        }
    }
}

//_________________________________________________________

_Matrix * _DataSetFilter::HarvestFrequencies (char unit, char atom, bool posSpec, bool countGaps) const {
    _SimpleList copy_seqs (theNodeMap), copy_sites (theOriginalOrder);
    return theData->HarvestFrequencies (unit,atom, posSpec, copy_seqs, copy_sites, countGaps);
}


/*
//_______________________________________________________________________
long    _DataSetFilter::GetVectorCode(long site,long seq)
{
    if (!symbolVector) return -1;
    long* fi = symbolVector->quickArrayAccess();
    return fi[*fi*site+seq+1];
}
//_______________________________________________________________________

void    _DataSetFilter::ProduceSymbolVector(bool smear)
{
    // compute the size of the vector cells
    _Parameter cellSize=log((_Parameter)theData->theTT->LengthOfAlphabet())*_Parameter(unitLength)/log(128.0);
    if (cellSize>2.0)
    {
        _String errMsg ("DataSetFilter has more than 32767 states, which is currently unsupported");
        FlagError(errMsg);
    }
    long intCellSize = cellSize>1.0?2:1;
    // now produce the conversion vector
    long sites = theMap.lLength, species= theNodeMap.lLength?theNodeMap.lLength:theData->NoOfSpecies();
    symbolVector = new _SimpleList ();
    checkPointer(symbolVector);
//  (*symbolVector)<<intCellSize;
    (*symbolVector)<<species;
    // we will now speciate into byte and word size cases
    // the data will be stored column by column
    // if there is a unique code translation, we then store that code in the symbol vector for faster
    // processing during tree pruning business.
    // use a standard convert to frequencies function to check whether a character has a unique convertion
    if (intCellSize==1) // char based storage
    {
        union
        {
            long composite;
            char bytes[sizeof(long)];
        } converterb;
        char byteposition = 0, bytesPerLong = sizeof(long);
        for (long i=0;i<sites;i++)
        {
            for (long j=0; j<species; j++)
            {
//              if (byteposition==bytesPerLong)
//              {
//                  byteposition = 0;
//                  (*symbolVector)<<converterb.composite;
//              }
//              converterb.bytes[byteposition]=(char)Translate2Frequencies((*this)(i,j),nil,smear,false);
//              byteposition++;
                (*symbolVector)<<Translate2Frequencies((*this)(i,j),nil,smear,false);
            }
        }
        if (byteposition)
        {
            for(long i=sizeof(long)-1;i>=byteposition;i--)
            {
                converterb.bytes[i]=0;
            }
            (*symbolVector)<<converterb.composite;
        }

    }
    else
    {
        union
        {
            long composite;
            short int words[sizeof(long)/2];
        } converterw;
        char wordposition = 0, wordsPerLong = sizeof(long);
        for (long i=0;i<sites;i++)
        {
            for (long j=0; j<species; j++)
            {
                if (wordposition==wordsPerLong)
                {
                    wordposition = 0;
                    (*symbolVector)<<converterw.composite;
                }
                converterw.words[wordposition]=(char)Translate2Frequencies((*this)(i,j),nil,smear,false);
            }
        }
        if (wordposition)
        {
            for(long i=sizeof(long)/2-1;i>=wordposition;i--)
            {
                converterw.words[i]=0;
            }
            (*symbolVector)<<converterw.composite;
        }
    }
}*/

//_______________________________________________________________________

long    _DataSetFilter::CorrectCode (long code) const {
    if (theExclusions.lLength == 0) {
      return code;
    }
    return theExclusions.SkipCorrect (code);
}



//_______________________________________________________________________
long    _DataSetFilter::Translate2Frequencies (_String const& str, _Parameter* parvect, bool smear) const {
    long  store      [HYPHY_SITE_DEFAULT_BUFFER_SIZE],
          resolution_count  = -1L;

  
    InitializeArray(parvect, dimension, 0.);
  
    if (unitLength == 1) {
        resolution_count = theData->theTT->TokenResolutions (str.sData[0],store,smear);
    } else {
        resolution_count = theData->theTT->MultiTokenResolutions(str,store, smear);
    }
  
    long mapped_resolution_count = theExclusions.lLength ? theExclusions.CorrectForExclusions(store, resolution_count) : resolution_count;
  
    /* handle the cases when no unambiguous resolutions were available */
    for (long i = 0L; i < mapped_resolution_count; i++) {
        parvect[store[i]] = 1.;
    }
  
    if (mapped_resolution_count == 1L) {
      return store[0];
    }
    if (mapped_resolution_count == 0L && resolution_count == 0L && smear) {
       InitializeArray(parvect, dimension, 1.);
    }
  
    return -1L;
}


//_______________________________________________________________________

long    _DataSetFilter::MapStringToCharIndex (_String& str) const {
  
  long  store      [HYPHY_SITE_DEFAULT_BUFFER_SIZE],
        resolution_count  = -1L;
  
  
  if (unitLength == 1) {
    resolution_count = theData->theTT->TokenResolutions (str.sData[0],store);
  } else {
    resolution_count = theData->theTT->MultiTokenResolutions(str,store);
  }
  
  long mapped_resolution_count = theExclusions.lLength ? theExclusions.CorrectForExclusions(store, resolution_count) : resolution_count;
  
  if (mapped_resolution_count == 1L) {
    return store[0];
  }
  
  return -1L;
}


//_______________________________________________________________________
long    _DataSetFilter::Translate2Frequencies (char s, _Parameter* parvect, bool smear) const {
  long  store      [HYPHY_SITE_DEFAULT_BUFFER_SIZE],
  resolution_count  = theData->theTT->TokenResolutions (s,store,smear);
  
  long mapped_resolution_count = theExclusions.lLength ? theExclusions.CorrectForExclusions(store, resolution_count) : resolution_count;
  
  if (mapped_resolution_count == 0L) {
    if (smear) {
      InitializeArray(parvect, dimension, 1.);
      return -1;
    }
  }
  
  InitializeArray (parvect, dimension, 0.);
  
  for (long i = 0L; i < mapped_resolution_count; i++) {
    parvect[store[i]] = 1.;
  }
  
  return resolution_count==1L?1L:-1L;
}

//_______________________________________________________________________
long    _DataSetFilter::LookupConversion (char s, _Parameter* parvect) const
{
    if (undimension==4) {
        long* cCache = conversionCache.lData+(s-40)*5;
        parvect[0] = cCache[0];
        parvect[1] = cCache[1];
        parvect[2] = cCache[2];
        parvect[3] = cCache[3];
        return cCache[4];

    } else {
        int idx = (s-40)*(undimension+1);
        for (long i=0; i<undimension; parvect[i++] = conversionCache.lData[idx++]) ;
        return conversionCache.lData[idx];
    }
}
//_______________________________________________________________________
bool   _DataSetFilter::ConfirmConversionCache() const {
  return conversionCache.lLength || unitLength > 3;
}

//_______________________________________________________________________
void    _DataSetFilter::SetupConversion (void)
{
    if (conversionCache.lLength) {
        return;
    }

    if ( unitLength==1 ) { // do stuff
        char c = 40;
        _Parameter *temp    = new _Parameter [undimension+1UL];

        while(c<127) {
            //InitializeArray(temp, undimension + 1UL, 0.0);
          
            Translate2Frequencies(c, temp, true);
          
            long resolution_count = -1;
          
          
            for (unsigned long i=0UL; i<undimension; i++) {
                long character_code_resolution =  (long)temp[i];
                conversionCache << character_code_resolution;
                if (character_code_resolution) {
                    if (resolution_count == -1) {
                        resolution_count = i;
                    } else {
                        resolution_count = -2;
                    }
                }
            }
          
            conversionCache<<resolution_count;
            c++;
        }
        delete[] temp;
    } else {
      
        if (unitLength==2 || unitLength==3) {
          
            _String alphabet = theData->theTT->GetAlphabetString();
            unsigned long alphabet_dim = alphabet.sLength;
          

            long  ccache [88],
                  uncorrected_dimension = GetDimension(false) ;

            conversionCache.RequestSpace (89+uncorrected_dimension);
            conversionCache << alphabet.sLength;

            for (unsigned long i=0UL; i<88; i++) {
                ccache[i] = -1;
            }
            for (unsigned long i=0UL; i<alphabet.sLength; i++) {
                ccache [alphabet.sData[i]-40] = i;
            }
            for (unsigned long i=0UL; i<88; i++) {
                conversionCache << ccache[i];
            }

             _String s (unitLength,false);
            for (unsigned long char_index = 0; char_index < uncorrected_dimension; char_index++ ) {
              _SimpleList components = SplitIntoDigits (char_index, unitLength, alphabet_dim);
              for (unsigned long position = 0; position < unitLength; position ++) {
                s.sData[position] = alphabet.sData[components (position)];
              }
              conversionCache << MapStringToCharIndex(s);
            }
        }
    }
}



//_________________________________________________________
//_________________________________________________________
// reading the data set file in here



//_________________________________________________________
void    checkTTStatus (FileState* fs) // check whether the translation table needs to be refreshed
{
    if (fs->translationTable == &defaultTranslationTable) {
        fs->translationTable =  (_TranslationTable*)defaultTranslationTable.makeDynamic();
    }
}
//_________________________________________________________
void    processCommand (_String*s, FileState*fs)
{
    // loop thru understood values of commands
    static _List CommandList;
    if (CommandList.lLength == 0)
        // first time in, should init commands
    {
        _String command ("BASESET");
        CommandList&& & command;
        command="FORMAT";
        CommandList&& & command;
        command="RAWLINE";
        CommandList&& & command;
        command="REPEAT";
        CommandList&& & command;
        command="TOKEN";
        CommandList&& & command;
    }

    long f = -1;
    long i,k = 0,l = 0,m;
    for (i=0; (i<CommandList.lLength); i++) {
        f = s->Find (*(_String*)CommandList(i));
        if (f!=-1) {
            break;
        }
    }

    if (f==-1) { // unrecognized command
        return;
    } else {
        // trim the string
        s->Trim (f+((_String*)CommandList(i))->Length(),-1);
        f = s->Find (":");
        if (f==-1) { // poorly formed command
            return;
        } else {
            s->Trim (f+1,-1);
        }

        if ((i>=1)&&(i<=3)) {
            k = s->Find ('\"');
            if (k==-1) {
                return;
            }
            l = s->Find ('\"', k+1,-1);
            if (l==-1) {
                return;
            }
            if (l<=k) {
                return;
            }
            s->Trim (k+1,l-1);
        }

        switch (i) {
            char c;
        case 4: // new token
            checkTTStatus (fs);
            // attempt to extract a token. Looking for (e.g):   "c" = "AC"
            k = s->Find ('"');
            if (k==-1) {
                return;
            }
            if ((*s)[k+2]!='"') {
                return;
            }
            l = s->Find ('"',k+3,-1);
            m = s->Find ('"',l+1,-1);
            if ((l==-1)||(m==-1)) {
                return;
            }

            c = (*s)[k+1];
            s->Trim (l+1,m-1);
            fs->translationTable->AddTokenCode (c,*s);
            break;


        case 0:// new code set, e.g  "ACGU"
            checkTTStatus(fs);
            // erase previous char definitions
            fs->translationTable->translationsAdded.Clear();
            fs->translationTable->tokensAdded = "";
            if (*s!=_String("BASE20")) {
                fs->translationTable->AddBaseSet (*s);
            } else {
                fs->translationTable->AddBaseSet (emptyString);
                fs->translationTable->baseLength = 20;
            }
            break;

        case 1: //FORMAT
            if (*s==_String("PHYLIPI")) { // PHYLIP Interleaved
                fs->fileType = 1;
                fs->interleaved = TRUE;
            } else if (*s==_String("PHYLIPS")) { // PHYLIP sequential
                fs->fileType = 1;
                fs->interleaved = FALSE;
            }
            if (*s==_String("RAW")) { // RAW Sequential Data (as in NEXUS)
                fs->fileType = 2;
                fs->interleaved = FALSE;
            }
            fs->autoDetect = false;
            break;

        case 3: // REPEAT CHAR
            fs->repeat = s->getChar(0);
            break;

        case 2: // RAWLINE template e.g 1,-1 skips one word at the beginning and one word at the end
            _List chips (s,',');
            for (int i = 0; i<chips.lLength; i++) {
                fs->rawLinesFormat<<(long)(((_String*)chips(i))->toNum());
            }

        }
    }
}
//_________________________________________________________

void    FilterRawString (_String& s, FileState* fs, _DataSet & ds)
{
    int i;
    for (i = 0; i<fs->rawLinesFormat.lLength; i++) {
        long f = fs->rawLinesFormat (i),p=0,l=0;
        if (f>0) {
            for (int j = 0; (j<f)&&(p>=0)&&(l>=0); j++) {
                p = s.FirstNonSpaceIndex(l,-1,1);
                l = s.FirstSpaceIndex(p,-1,1);
            }
            if (l<0) {
                break;
            }
            p = s.FirstNonSpaceIndex(l,-1,1);
            s.Trim(p,-1);
        } else {
            if (f!=0) {
                p = 0;
                l = 0;
                for (int j = 0; (j>f)&&(p>=0)&&(l>=0); j--) {
                    p = s.FirstNonSpaceIndex(p,-1,-1);
                    l = s.FirstSpaceIndex(0,p,-1);
                }
                if (l<0) {
                    break;
                }
                p = s.FirstNonSpaceIndex(0,l,-1);
                s.Trim(0,p);
            } else {
                // Name
                p = s.FirstNonSpaceIndex();
                l = s.FirstSpaceIndex(p+1,-1,1);
                if ((p<0)||(l<0)) {
                    break;
                }
                _String Name = s.Cut (p,l-1);
                ds.AddName (Name);
                s.Trim (s.FirstNonSpaceIndex(l,-1,1),-1);
            }
        }

    }
    if (i!=fs->rawLinesFormat.lLength) {
        s = "";
    }
}
//_________________________________________________________
void    ProcessTree (FileState *fState, FILE* f, _String& CurrentLine)
{
    long j = 0, i=0; // parenthesis balance
    char c;
    _String treeString ((unsigned long)10, true);
    do {
        for (i=0; i<CurrentLine.sLength; i++) {
            c = CurrentLine.sData[i];
            if (!isspace (c)) {
                treeString<<c;
                if (c==')') {
                    j--;
                    if (!j) {
                        break;
                    }
                } else if (c=='(') {
                    j++;
                }
            }

        }
        ReadNextLine (f,&CurrentLine,fState, false);
    } while (j&&CurrentLine.sLength);

    if (j) {
        _String errMsg ("Tree string found in data file had unbalanced parentheses: ");
        if (j>0) {
            errMsg = errMsg & j & " too few closing parentheses.";
        } else {
            errMsg = errMsg & (-j) & " too many closing parentheses.";
        }
        ReportWarning (errMsg);
    } else {
        treeString.Finalize();
        setParameter (dataFileTree,1.0,fState->theNamespace);
        setParameter (dataFileTreeString, new _FString (treeString), nil, false);
    }

}

//_________________________________________________________

long    ProcessLine (_String&s , FileState *fs, _DataSet& ds)
{
    long sitesAttached = 0,sL=s.Length();

    for (long l = 0; l<sL; l++) {
        // see if it is a legal char
        char c = toupper (s.sData[l]);
        if (fs->translationTable->IsCharLegal(c)) { // go on
            if (fs->curSpecies==0) { // add new column
                ds.AddSite (c);
                sitesAttached++;
            } else { //append to exisiting column
                //if (c == fs->skip) continue;
                // check to see if this species needs to be padded

                if (c == fs->repeat) {
                    if (fs->curSite+sitesAttached >= ds.lLength) { // a dot not matched by a previously read character; ignore
                        return sitesAttached;
                    }

                    c = ((_Site*)(ds._List::operator () (fs->curSite+sitesAttached)))->getChar(0);
                    if (c==0)
                        c = ((_Site*)(ds._List::operator ()
                                      (((_Site*)(ds._List::operator () (fs->curSite+sitesAttached)))->GetRefNo())))->getChar(0);
                }

                if (fs->curSite+sitesAttached+1>fs->totalSitesRead) {
                    // pad previous species to full length
                    _Site * newS = new _Site (fs->skip);
                    checkPointer (newS);

                    for (long j = 1; j<fs->curSpecies; j++) {
                        (*newS) << fs->skip;
                    }

                    (*newS) << c;

                    /*long rN = ds.dsh->incompletePatterns->Find (newS);

                    if (rN>=0)
                    {
                        rN =  ds.dsh->incompletePatterns->GetXtra (rN);
                        ds.theFrequencies[rN]++;
                        newS->Clear();
                        newS->SetRefNo(rN);
                        ds.theFrequencies << 0;
                    }
                    else
                    {*/
                    ds.theFrequencies << 1;
                    newS->SetRefNo(-1);
                    //}

                    ds << newS;
                    newS->nInstances --;

                    fs->totalSitesRead++;
                } else {
                    ds.Write2Site (fs->curSite+sitesAttached, c);
                }

                sitesAttached++;
            }
        }
    }
    // make sure that this species has enough data in it, and if not - pad it with '?'

    if ((fs->curSite+sitesAttached<fs->totalSitesRead)&&(fs->interleaved)) {
        // pad this species to full length
        for (long j = fs->curSite+sitesAttached; j<fs->totalSitesRead; j++) {
            ds.Write2Site (j, fs->skip);
        }
    }
    if (!fs->curSpecies) {
        fs->totalSitesRead+=sitesAttached;
    }
    return sitesAttached;
}

//_________________________________________________________
void    PadLine (FileState& fState, _DataSet& result) // make sure that there is enough data in this line
// and if not - "pad" it with '?''s
{
    if (fState.curSite<fState.totalSitesRead) // pad line if needed
        for (long j = fState.curSite; j<fState.totalSitesRead; j++) {
            result.Write2Site (j, fState.skip);
        }
}

//_________________________________________________________
void    ISelector (FileState& fState, _String& CurrentLine, _DataSet& result)
{
    if (fState.interleaved) { // interleaved file
        if (fState.curSpecies&&(!((fState.curSpecies)%fState.totalSpeciesExpected))) { // read a chunk of all species
            if (fState.totalSitesRead && !result.InternalStorageMode()) {
                for (long i = fState.curSite; i<fState.totalSitesRead; i++) {
                    result.Compact(i);
                }

                result.ResetIHelper();

            }
            fState.curSite = fState.totalSitesRead;
            fState.curSpecies = 0;
            ProcessLine (CurrentLine, &fState, result);
            fState.curSpecies = 1;
            if (!fState.curSite) {
                fState.totalSpeciesRead++;
            }
        } else {
            ProcessLine (CurrentLine, &fState, result);
            if (!fState.curSite) {
                fState.totalSpeciesRead++;
            }
            fState.curSpecies++;
        }
    } else {
        if (fState.curSpecies+1<fState.totalSpeciesExpected) {
            fState.curSpecies++;
        }
        if (fState.curSpecies == fState.totalSpeciesRead) {
            PadLine (fState, result);
            fState.curSite = 0;
        }
        if (fState.totalSpeciesRead<fState.totalSpeciesExpected) {
            fState.totalSpeciesRead++;
        }

        fState.curSite+=ProcessLine (CurrentLine, &fState, result);

    }
}

//_________________________________________________________
bool SkipLine (_String& theLine, FileState* fS)
{
    if ( theLine.sData[0]=='/' && theLine.sData[1]=='/' ) {
        return true;
    }

    char c = theLine.FirstNonSpace();

    if (c&& (!((c=='$')&&(!fS->acceptingCommands))) ) {
        return false;
    }

    return true;
}

//_________________________________________________________

#define     READ_NEXT_LINE_BUFFER_SIZE      1024*1024


//_________________________________________________________
void ReadNextLine (FILE* fp, _String *s, FileState* fs, bool, bool upCase)
{
    _String  tempBuffer (1024L, true);

    char lastc;

    if (fp) {
        lastc = fgetc(fp);
    } else {
        lastc = fs->pInSrc<fs->theSource->sLength?fs->theSource->sData[fs->pInSrc++]:0;
    }

    if (fs->fileType != 3) { // not NEXUS - do not skip [..]
        if (fp)
            while ( !feof(fp) && lastc!=10 && lastc!=13 ) {
                if (lastc) {
                    tempBuffer << lastc;
                }

                lastc = fgetc(fp);
            }
        else
            while (lastc && lastc!=10 && lastc!=13 ) {
                tempBuffer << lastc;
                lastc = fs->theSource->sData[fs->pInSrc++];
            }

    } else {
        if (upCase) {
            lastc = toupper(lastc);
        }

        while (((fp&&(!feof(fp)))||(fs->theSource&&(fs->pInSrc<=fs->theSource->sLength))) && lastc!=10 && lastc!=13) {
            if (lastc=='[') {
                if (fs->isSkippingInNEXUS) {
                    ReportWarning ("Nested comments in NEXUS really shouldn't be used.");
                } else {
                    fs->isSkippingInNEXUS = true;
                }
            }
            if (fs->isSkippingInNEXUS) {
                if (lastc==']') {
                    fs->isSkippingInNEXUS = false;
                    tempBuffer << ' ';
                }
            } else {
                tempBuffer << lastc;
            }

            if (fp) {
                if (upCase) {
                    lastc = toupper(fgetc(fp));
                } else {
                    lastc = fgetc(fp);
                }
            } else {
                if (upCase) {
                    lastc = toupper(fs->theSource->sData[fs->pInSrc++]);
                } else {
                    lastc = fs->theSource->sData[fs->pInSrc++];
                }
            }

        }

        if ( lastc==10 || lastc==13 ) {
            tempBuffer << ' ';
        }
    }

    tempBuffer.Finalize();

    if ( (fp && feof(fp)) || (fs->theSource && fs->pInSrc >= fs->theSource->sLength) )
        if (tempBuffer.sLength == 0) {
            *s = "";
            return;
        }

    if (s->nInstances > 1) {
        *s = tempBuffer;
    } else {
        Ptr         saveData = s->sData;
        s->sData    = tempBuffer.sData;
        tempBuffer.sData = saveData;

        s->sLength  = tempBuffer.sLength;
    }

    if (SkipLine (*s, fs)) {
        ReadNextLine(fp,s,fs,false,upCase);
    }

    if (s->sLength && s->sData[s->sLength-1]== '\n') {
        s->Trim (0,s->sLength-2);
    }
}
//_________________________________________________________
void    TrimPhylipLine (_String& CurrentLine, _DataSet& ds)
{
    int  fNS      = CurrentLine.FirstNonSpaceIndex(),
         space2   = CurrentLine.FirstSpaceIndex (fNS + 1);
    
    // hack for PAML support
    if (space2 > fNS && isspace(CurrentLine.getChar (space2+1))) {
        _String     sequence_name (CurrentLine,fNS, space2);
        CurrentLine.Trim(space2+2,-1); // chop out the name
        ds.AddName(sequence_name);        
    } else {
        _String     sequence_name (CurrentLine,fNS, fNS+9);
        CurrentLine.Trim(fNS+10,-1); // chop out the name
        ds.AddName(sequence_name);
    }
}


//_________________________________________________________
_DataSet* ReadDataSetFile (FILE*f, char execBF, _String* theS, _String* bfName, _String* namespaceID, _TranslationTable* dT)
{

    bool     doAlphaConsistencyCheck = true;
    _String::storageIncrement = 16;
    _DataSet* result = new _DataSet;


    _String         CurrentLine = dataFilePartitionMatrix & "={{}};",
                    savedLine;

    if (1) {
        _ExecutionList reset (CurrentLine);
        reset.Execute();
#ifdef __HYPHYMPI__
        if (_hy_mpi_node_rank == 0)
#endif
            terminateExecution = false;
    }

    // initialize the instance of a file state variable
    setParameter(dataFileTree, 0.0);
    FileState   fState;
    fState.translationTable =  dT;
    fState.curSpecies =
        fState.totalSpeciesRead =
            fState.totalSitesRead =
                fState.totalSpeciesExpected =
                    fState.totalSitesExpected =
                        fState.curSite =
                            fState.maxStringLength   = 0;
    fState.acceptingCommands = true;
    fState.allSpeciesDefined = false;
    fState.interleaved       = false;
    fState.isSkippingInNEXUS = false;
    fState.autoDetect        = true;
    fState.fileType          = -1;
    fState.baseLength        = 4;
    fState.repeat            = '.',
           fState.skip            = 0;
    fState.theSource         = theS;
    fState.pInSrc            = 0;
    fState.theNamespace      = namespaceID;

    if (!(f||theS)) {
        CurrentLine = "ReadDataSetFile received null file AND string references. At least one must be specified";
        warnError     (CurrentLine);
    }
    // done initializing

    long     fileLength = 0,
             lastDone = 10,
             cDone;

#ifdef __HYPHYMPI__
    if (_hy_mpi_node_rank == 0) {
#endif
        if       (f) {
            fseek    (f,0,SEEK_END);
            fileLength = ftell(f);
            rewind  (f);
        } else {
            fileLength = theS->sLength;
        }

#ifdef __HYPHYMPI__
    }
#endif



    //if (f==NULL) return (_DataSet*)result.makeDynamic();
    // nothing to do

    CurrentLine = emptyString;

    ReadNextLine (f,&CurrentLine,&fState);
    if (!CurrentLine.sLength) {
        CurrentLine = "emptyString File Encountered By ReadDataSet.";
        WarnError (CurrentLine);
        return result;
    } else {
        if (CurrentLine.beginswith ("#NEXUS",false)) {
            ReadNexusFile (fState,f,(*result));
            doAlphaConsistencyCheck = false;
        } else {
            long i,j,k, filePosition = -1, saveSpecExpected = 0x7FFFFFFF;
            char c;
            while (CurrentLine.sLength) { // stuff to do
                // check if the line has a command in it
#ifdef __HYPHYMPI__
                if (_hy_mpi_node_rank == 0) {
#endif
                    if (f) {
                        cDone = ftell (f)*100./fileLength;
                    } else {
                        cDone = fState.pInSrc*100./fileLength;
                    }

                    if (cDone>lastDone) {
                        SetStatusBarValue (lastDone,1,0);
#ifdef __MAC__
                        handleGUI(true);
#endif
                        lastDone+=10;
                    }
#ifdef __HYPHYMPI__
                }
#endif

                c = CurrentLine.FirstNonSpace();
                while (1) {
                    if (fState.acceptingCommands) {
                        if (c == '$') { // command line
                            processCommand(&CurrentLine, &fState);
                            break;
                        }
                    }

                    if (!fState.skip) {
                        fState.skip = fState.translationTable->GetSkipChar();
                    }
                    fState.acceptingCommands = FALSE;

                    if (fState.fileType==-1) { // undecided file type - assume it is PHYLIP sequential
                        if ((c == '#')||(c=='>')) { // hash-mark format
                            fState.fileType = 0;
                        } else { // assume this is a sequential PHYLIP file
                            fState.fileType = 1;
                            fState.interleaved = false;
                        }

                    }
                    // decide what to do next
                    // if format is PHYLIP and we do not know the expected dimensions,
                    //   we must read those in first
                    if (fState.fileType==1) { // PHYLIP
                        if ((filePosition<0)&&(fState.autoDetect)) {
                            filePosition = (f?
                                            ftell (f)
#ifdef __WINDOZE__
                                            -1
#endif
                                            :fState.pInSrc);
                            savedLine = CurrentLine;
                        }

                        if ((fState.totalSitesExpected==0)||(fState.totalSpeciesExpected==0)) { // must read dimensions first
                            i = CurrentLine.FirstNonSpaceIndex();
                            j = CurrentLine.FirstSpaceIndex(i,-1,1);
                            if (j>=0) {
                                k = CurrentLine.FirstNonSpaceIndex(j,-1,1);
                                if (k>=0) { // could have dimensions
                                    saveSpecExpected = fState.totalSpeciesExpected=CurrentLine.Cut(i,j-1).toNum();
                                    fState.totalSitesExpected=CurrentLine.Cut(k,-1).toNum();
                                }
                                if (CurrentLine.Find ('I', k, -1)>=0) { // interleaved
                                    fState.interleaved = true;
                                }
                            }
                        } else
                            // now for the data crunching part
                            // detect a line, diagnose it and dispatch accordingly
                        {
                            if (fState.interleaved) {
                                if (fState.totalSpeciesRead<fState.totalSpeciesExpected) {
                                    TrimPhylipLine (CurrentLine, (*result));
                                }
                                if (fState.curSite && fState.curSpecies >= saveSpecExpected &&
                                        fState.totalSitesRead >= fState.totalSitesExpected) {
                                    // reached the end of the data - see maybe there is a tree
                                    ReadNextLine (f,&CurrentLine,&fState);
                                    if (CurrentLine.sLength) {
                                        if (CurrentLine.FirstNonSpace()=='(') { // could be a tree string
                                            ProcessTree (&fState,f, CurrentLine);
                                        }
                                    }
                                    break;
                                }

                            } else {
                                if (fState.totalSitesRead > fState.totalSitesExpected)
                                    // oops - autodetect incorrectly assumed that the file was sequential
                                {
                                    fState.curSpecies =
                                        fState.totalSpeciesRead =
                                            fState.totalSitesRead =
                                                fState.curSite =
                                                    fState.totalSpeciesExpected =
                                                        fState.totalSitesExpected =
                                                            fState.maxStringLength = 0;
                                    fState.allSpeciesDefined = false;
                                    fState.interleaved = true;
                                    fState.autoDetect = true;

                                    if(f) {
                                        fseek (f, filePosition, SEEK_SET);
                                    } else {
                                        fState.pInSrc = filePosition;
                                    }

                                    CurrentLine = savedLine;
                                    for (long idx = 0; idx < (*result).lLength; idx++) {
                                        ((_Site*)(*result).lData[idx])->Finalize();
                                    }
                                    (*result).theNames.Clear();
                                    (*result).theMap.Clear();
                                    (*result).Clear();
                                    (*result).theFrequencies.Clear();
                                    if ((*result).dsh) {
                                        (*result).dsh->incompletePatterns->Clear(false);
                                        delete ((*result).dsh);
                                        (*result).dsh = nil;
                                    }
                                    continue;
                                }
                                if (fState.totalSpeciesRead==0) {
                                    fState.totalSpeciesExpected = 1;
                                    if (!fState.curSite) {
                                        TrimPhylipLine (CurrentLine, (*result));
                                    }
                                }

                                else if (fState.curSite>=fState.totalSitesExpected) {
                                    fState.totalSpeciesExpected++;
                                    if (fState.totalSpeciesExpected>saveSpecExpected) {
                                        // reached the end of the data - see maybe there is a tree
                                        ReadNextLine (f,&CurrentLine,&fState);
                                        if (CurrentLine.sLength) {
                                            if (CurrentLine.FirstNonSpace()=='(') { // could be a tree string
                                                ProcessTree (&fState,f, CurrentLine);
                                            }
                                        }
                                        break;
                                    }
                                    TrimPhylipLine (CurrentLine, (*result));
                                }
                            }

                            ISelector (fState, CurrentLine, (*result));
                        }
                        break;
                    }
                    // that's all for PHYLIP

                    // now handle raw data case
                    if (fState.fileType == 2) { // raw data
                        FilterRawString(CurrentLine, &fState, (*result));
                        if (CurrentLine.sLength) {
                            break;
                        }
                        if (ProcessLine (CurrentLine, &fState, (*result))) {
                            fState.curSpecies++;
                            fState.totalSpeciesRead++;
                        }
                        break;
                    }

                    // lastly, handle the auto-detect standard case

                    // check to see if the string defines a name
                    if ((c=='#')||(c=='>')) { // a name it is
                        if (fState.allSpeciesDefined) { // can't define the species after data
                            break;
                        } else {
                            if ((!fState.totalSpeciesRead)&&(fState.totalSpeciesExpected>=1)) {
                                fState.interleaved = TRUE;
                            } else {
                                fState.interleaved = FALSE;
                            }
                            fState.totalSpeciesExpected++;
                            CurrentLine.Trim(CurrentLine.FirstNonSpaceIndex(1,-1,1),-1);
                            if ((CurrentLine.sData[0]=='#')||(CurrentLine.sData[0]=='>')) {
                                CurrentLine = _String("Species")&_String(fState.totalSpeciesExpected);
                            }
                            (*result).AddName (CurrentLine);
                        }
                        break;
                    }
                    // check to see if the string defines a tree
                    if (c=='(') {
                        ProcessTree (&fState,f, CurrentLine);
                        ReadNextLine (f,&CurrentLine,&fState);
                    }

                    // check to see where to stick the incoming line

                    if (!fState.totalSpeciesExpected)
                        // raw data fed before names defined - skip
                    {
                        break;
                    }
                    if((fState.totalSpeciesExpected>1)&&(!fState.totalSpeciesRead)) {
                        fState.allSpeciesDefined = TRUE;
                    }

                    // repeat the structure of PHYLIP reader

                    ISelector (fState, CurrentLine, (*result));

                    break;
                }

                ReadNextLine (f,&CurrentLine,&fState);

            }
        }
    }



    if (fState.totalSitesRead && fState.interleaved && !result->InternalStorageMode()) {
        for (int i = fState.curSite; i<fState.totalSitesRead; i++) {
            (*result).Compact(i);
        }
        (*result).ResetIHelper();
    }

    if ((!fState.interleaved)&&(fState.fileType!=2)) {
        PadLine (fState, (*result));
    }


#ifdef __HYPHYMPI__
    if (_hy_mpi_node_rank == 0) {
#endif
        SetStatusBarValue (-1,1,0);
#ifdef __MAC__
        handleGUI(true);
#endif
#ifdef __HYPHYMPI__
    }
#endif

    // make sure interleaved duplications are handled correctly

    (*result).Finalize();
    _String::storageIncrement   = 64;
    (*result).noOfSpecies       = fState.totalSpeciesRead;
    (*result).theTT             = fState.translationTable;

    // check to see if result may be an amino-acid data
    if (doAlphaConsistencyCheck && result->theTT == &defaultTranslationTable) {
        if (result->GetNoTypes() == 0)
            // emptyString data set
            // try binary data
        {
            _TranslationTable *trialTable = new _TranslationTable (defaultTranslationTable);
            trialTable->baseLength = 2;
            _DataSet * res2 = ReadDataSetFile (f, execBF, theS, bfName, namespaceID, trialTable);
            if (res2->GetNoTypes()) {
                DeleteObject (result);
                return res2;
            }
            DeleteObject (trialTable);
        } else
            // check it out
            if (result->CheckAlphabetConsistency()<0.5)
                // less than 50% of the data in the alphabet is not in the basic alphabet
            {
                _TranslationTable trialTable (defaultTranslationTable);
                trialTable.baseLength = 20;
                (*result).theTT = &trialTable;
                if ((*result).CheckAlphabetConsistency()<0.5) {
                    CurrentLine = "More than 50% of characters in the data are not in the alphabet.";
                    (*result).theTT =  &defaultTranslationTable;
                    ReportWarning (CurrentLine);
                } else {
                    (*result).theTT = (_TranslationTable*)trialTable.makeDynamic();
                }

            }

    }
    if (nexusBFBody.sLength) {
        if (execBF == 1) {
            lastNexusDataMatrix = result;

            long            bfl = GetBFFunctionCount ();

            _ExecutionList nexusBF (nexusBFBody,namespaceID);
            if (bfName) {
                nexusBF.sourceFile = *bfName;
            }

#ifndef __UNIX__
#ifdef __HYPHYMPI__
            if (_hy_mpi_node_rank == 0)
#endif
                ApplyPreferences();
#endif

            nexusBF.ExecuteAndClean(bfl);
            nexusBFBody         = emptyString;
        } else if (execBF == 0) {
            nexusBFBody         = emptyString;
        }
    }

    return result;
}

//_________________________________________________________

BaseRef _DataSetFilter::toStr (unsigned long)
{
    //return new _String("DataSetFilters only print to files");
    _String * res = new _String (4096L, true);
    internalToStr (nil,*res);
    res->Finalize();
    return res;
}

//_________________________________________________________

void    _DataSetFilter::PatternToSiteMapper (void* source, void* target, char mode, long padup) const {
  
  unsigned long site_count = GetSiteCountInUnits();
  
  switch (mode) {
    case 0: {
      _Parameter * target_array = (_Parameter*) target,
                 * source_array = (_Parameter*) source;
      
      for (unsigned site = 0UL; site < site_count; site++ ) {
        target_array [site] = source_array [duplicateMap.lData[site]];
      }
      for (long site = duplicateMap.lLength; site < padup; site++) {
        target_array [site] = 1.;
      }

      break;
    }
    case 1: {
      long * target_array = (long*) target,
           * source_array = (long*) source;
      
      for (unsigned site = 0UL; site < site_count; site++ ) {
        target_array [site] = source_array [duplicateMap.lData[site]];
      }
      for (long site = duplicateMap.lLength; site < padup; site++) {
        target_array [site] = 0;
      }
      
      break;
    }
    case 2: {
      long * target_array = (long*) target;
      _Parameter       * source_array = (_Parameter*) source;
      
      for (unsigned site = 0UL; site < site_count; site++ ) {
        target_array [site] = source_array [duplicateMap.lData[site]];
      }
      break;
    }
  }
}


//_________________________________________________________

long    _DataSetFilter::GetOriginalToShortMap(long index)
{
    long pos1=theData->theMap.lData[theOriginalOrder.lData[index]],pos2;
    pos2 = theMap.Find(pos1);
    if (pos2==-1) {
        for (long i=theData->theMap.lLength-1; i>=0; i--) {
            if (theData->theMap.lData[i]==pos1) {
                pos2 = theMap.Find(i);
                if (pos2!=-1) {
                    break;
                }
            }
        }
    }
    return pos2;
}

//_________________________________________________________

_String const &_DataSetFilter::GenerateConsensusString (_SimpleList* majority) const {
  
    if (unitLength > 3) {
        return emptyString;
    }

    _String     result ((unsigned long)theOriginalOrder.lLength),
                pattern_consensus  ((unsigned long)(unitLength*theFrequencies.lLength));

    long        char_states         = GetDimension(false),
                *translation_buffer = new long [char_states];
 
    _Parameter* count_buffer = new _Parameter [char_states];
  
     for (unsigned long site_pattern = 0UL; site_pattern<theFrequencies.lLength; site_pattern ++) {
        long    index_in_dataset = theMap.lData[site_pattern];

        InitializeArray (count_buffer, char_states, 0.);

        for (unsigned long sequence_index =0UL; sequence_index < theNodeMap.lLength; sequence_index ++) {
            long resolution_count = theData->theTT->TokenResolutions ((*theData)(index_in_dataset, theNodeMap.lData[sequence_index],1),translation_buffer, false);
          
          
            if (resolution_count>1L) {
                _Parameter equal_weight = 1./resolution_count;
                for (long resolution_index = 0L; resolution_index < resolution_count; resolution_index++) {
                  count_buffer [translation_buffer[resolution_index]] += equal_weight;
                }
            } else {
                if (resolution_count == 1) {
                    count_buffer [translation_buffer[0]] += 1.;
                }
            }
        }

        // find the residue with the highest frequency
       
        _Parameter       max_weight      = -1.;
        InitializeArray (translation_buffer, char_states, 0L);
       long             max_char_count  = 0L;
       
        for (unsigned long char_index = 0UL; char_index < char_states; char_index++) {
          if (StoreIfGreater(max_weight, count_buffer[char_index])) {
            max_char_count = 1;
            translation_buffer [0] = char_index;
          } else {
            if (count_buffer[char_index] == max_weight) {
              translation_buffer [max_char_count ++] = char_index;
            }
          }
        }
         
         if (max_char_count > 1L) {
            pattern_consensus.sData[site_pattern]=theData->theTT->AmbigToLetter(translation_buffer, max_char_count);
        } else {
            pattern_consensus.sData[site_pattern] = theData->theTT->ConvertCodeToLetters(translation_buffer[0],1) [0];
        }
         if (majority) {
           (*majority) << max_weight;
         }
   }

    delete [] count_buffer;
    delete [] translation_buffer;

    for (unsigned long m=0UL; m<theOriginalOrder.lLength; m++) {
        result.sData[m] = pattern_consensus.sData[duplicateMap.lData[m]];
    }

    return result;
}


//_________________________________________________________
void    _DataSetFilter::toFileStr (FILE*dest, unsigned long) {
// write out the file with this dataset filter
    if (dest) {
      _String       dummy;
      internalToStr (dest,dummy);
   }
}

//_________________________________________________________
void    _DataSetFilter::ConvertCodeToLettersBuffered (long code, char unit, char* storage, _AVLListXL* lookup) const {
    // write out the file with this dataset filter
    long            lookupC     = lookup->Find ((BaseRef)code);
    const char      *lookupV;
    if (lookupC>=0) {
        lookupV = ((_String*)lookup->GetXtra(lookupC))->sData;
    } else {
        _String * newT = new _String (ConvertCodeToLetters (code,unit));
        lookup->Insert ((BaseRef)code, (long)newT, false);
        lookupV = newT->sData;
    }
  
    if (unit == 1) {
      storage[0] = lookupV[0];
    } else {
      for (unsigned long k = 0UL; k < unit; k++) {
          storage[k] = lookupV[k];
      }
    }
}






//_________________________________________________________

void    _DataSetFilter::internalToStr (FILE * file ,_String& string_buffer) {
  
  
  auto trim_to_10 = [] (const _String& seq_name) -> _String const& {
    if (seq_name.Length() >= 10) {
      return seq_name.Cut (0,9) & ' ';
    }
    return seq_name & _String (_String (" "), 11-seq_name.Length());
  };
  
  // write out the file with this dataset filter
  checkParameter (dataFilePrintFormat,dFPrintFormat,6.0);
  checkParameter (dataFileDefaultWidth,dFDefaultWidth,50.0);
  _Parameter  gW;
  
  long outputFormat = dFPrintFormat,
  printWidth   = dFDefaultWidth,
  gapWidth;
  
  unsigned long sequence_count = NumberSpecies(),
  site_count     = GetSiteCount();
  
  checkParameter (dataFileGapWidth,gW,10.0);
  if(!printWidth) {
    printWidth = 50;
  }
  
  gapWidth = gW;
  if (gapWidth<=0) {
    gapWidth = printWidth;
  }
  
  StringFileWrapper write_here (file ? nil : & string_buffer, file);
  
  if (outputFormat < 4 || outputFormat > 8) {
    // not NEXUS or serial
    if (!(theData->theTT->IsStandardNucleotide() || theData->theTT->IsStandardAA())) {
      _String * bSet = &theData->theTT->baseSet;
      
      write_here << "$BASESET:\""
      << *bSet
      << "\"\n";
      
      if (theData->theTT->tokensAdded.sLength) {
        for (long at = 0; at < theData->theTT->tokensAdded.sLength; at++) {
          write_here << "$TOKEN:\""
          << theData->theTT->tokensAdded.sData[at]
          << "\" = \""
          << theData->theTT->ExpandToken (theData->theTT->tokensAdded.sData[at])
          << "\"\n";
        }
      }
    }
  }
  
  switch (outputFormat) {
    case 1: // hash-mark interleaved
    case 10: { // FASTA interleaved
      
      long sitesDone    = 0,
      upTo;
      
      char seqDelimiter = (outputFormat==1)?'#':'>';
      
      for (unsigned long i = 0UL; i<theNodeMap.lLength; i++) {
        write_here << seqDelimiter
        << GetSequenceName(i)
        << kStringFileWrapperNewLine;
      }
      
      while (sitesDone<theOriginalOrder.lLength) {
        
        write_here << kStringFileWrapperNewLine
        << kStringFileWrapperNewLine;
        
        
        upTo = sitesDone+printWidth;
        if (upTo>theOriginalOrder.lLength) {
          upTo = theOriginalOrder.lLength;
        }
        
        for (unsigned long i = 0UL; i<theNodeMap.lLength; i++) {
          for (unsigned long j = sitesDone; j<upTo; j++) {
            if ((j-sitesDone)%gapWidth==0) {
              write_here << ' ';
            }
            write_here << (*theData)(theOriginalOrder.lData[j],theNodeMap.lData[i],1);
          }
          
          write_here << kStringFileWrapperNewLine;
        }
        
        sitesDone = upTo;
      }
      break;
    }
      
    case 2:     // PHYLIP sequential
    case 11:    // PAML
    {
      
      write_here << _String((long)theNodeMap.lLength)
      << kStringFileWrapperTab
      << _String(theOriginalOrder.lLength)
      << kStringFileWrapperNewLine;
      
      // proceed to spool out the data
      for (unsigned long i = 0UL; i<theNodeMap.lLength; i++) {
        _String const * sequence_name = GetSequenceName(i);
        _String sequence_name_10;
        
        if (outputFormat == 2) { // PHYLIP
          sequence_name_10 = trim_to_10 (*sequence_name);
        } else {
          sequence_name_10 = *sequence_name & "  ";
        }
        
        write_here << sequence_name_10;
        
        for (unsigned long site_index = 0; site_index<theOriginalOrder.lLength; site_index++) {
          if ((site_index%printWidth==0)&&site_index) {
            write_here << "\n           ";
          }
          
          write_here << (*theData)(theOriginalOrder(site_index),theNodeMap(i),1);
          if (site_index%gapWidth==gapWidth-1) {
            write_here << ' ';
          }
        }
        write_here << kStringFileWrapperNewLine;
        
      }
      break;
    }
      
    case 3: { // phylip interleaved
              // print PHYLIP format header
              //fprintf (dest,"$FORMAT:\"PHYLIPI\"\n");
              // print number of species and sites
      
      write_here << _String((long)theNodeMap.lLength)
      << kStringFileWrapperTab
      << _String(theOriginalOrder.lLength)
      << kStringFileWrapperNewLine;
      
      // proceed to spool out the data
      for (unsigned long i = 0UL; i<theNodeMap.lLength; i++) {
        write_here << trim_to_10 (*GetSequenceName(i));
        
        for (unsigned long j = 0UL; j<theOriginalOrder.lLength; j++) {
          if (j==printWidth) {
            write_here << kStringFileWrapperNewLine;
          } else {
            if (j%gapWidth==0) {
              write_here << ' ';
            }
            write_here << (*theData)(theOriginalOrder.lData[j],theNodeMap.lData[i],1);
          }
        }
      }
      
      
      
      unsigned long completed = printWidth;
      
      while (completed<theOriginalOrder.lLength-1) {
        long upTo = completed+printWidth<theOriginalOrder.lLength?completed+printWidth:theOriginalOrder.lLength;
        for (unsigned long i = 0UL; i<theNodeMap.lLength; i++) {
          write_here << "\n           ";
          for (unsigned long j = completed; j<upTo; j++) {
            if ((j-completed)%gapWidth==0) {
              write_here <<  ' ';
            }
            write_here << (*theData)(theOriginalOrder.lData[j],theNodeMap.lData[i],1);
          }
        }
        completed+=printWidth;
        write_here << kStringFileWrapperNewLine;
      }
      
      
      break;
    }
      
      // various flavors of NEXUS
      
    case 4: // labels, sequential
    case 5: // labels, interleaved
    case 6: // no labels, sequential
    case 7: { // no labels, interleaved
              // write out the header
      
      write_here << "#NEXUS\n\nBEGIN TAXA;\n\tDIMENSIONS NTAX = "
      << _String ((long)sequence_count)
      << ";\n\tTAXLABELS\n\t\t";
      
      for (unsigned long i=0UL; i< sequence_count; i++) {
         write_here << GetSequenceName(i)->Enquote('\'') << ' ';
      }
      
      write_here << ";\nEND;\n\nBEGIN CHARACTERS;\n\tDIMENSIONS NCHAR = "
      << _String((long)theOriginalOrder.lLength)
      << ";\n\tFORMAT\n\t\t";
      
      if (theData->theTT->IsStandardNucleotide()) {
        write_here << "DATATYPE = DNA\n";
      } else {
        if (theData->theTT->IsStandardAA()) {
          write_here << "DATATYPE = PROTEIN\n";
        } else if (theData->theTT->IsStandardBinary()) {
          write_here << "DATATYPE = BINARY\n";
        } else {
          long alphabet_length = theData->theTT->baseSet.sLength;
          
          write_here << "\t\tSYMBOLS = \"";
          for (unsigned long bc = 0UL; bc < alphabet_length-1; bc++) {
            write_here << theData->theTT->baseSet.getChar (bc)
            << ' ';
          }
          write_here << theData->theTT->baseSet.getChar (alphabet_length-1)
          << "\"\n";
          
          if (theData->theTT->tokensAdded.sLength)
            for (long at = 0; at < theData->theTT->tokensAdded.sLength; at++) {
              write_here << "\nEQUATE =\""
              << theData->theTT->tokensAdded.sData[at]
              << " = "
              << theData->theTT->ExpandToken(theData->theTT->tokensAdded.sData[at])
              << "\"";
            }
        }
      }
      if (theData->theTT->GetGapChar()) {
        write_here << "\t\tGAP=" << theData->theTT->GetGapChar();
      }
      if (theData->theTT->GetSkipChar()) {
        write_here << "\n\t\tMISSING=" << theData->theTT->GetSkipChar();
      }
      if (outputFormat>5) {
        write_here << "\n\t\tNOLABELS";
      }
      if (outputFormat%2) {
        write_here << "\n\t\tINTERLEAVE";
      }
      
      write_here << "\n\t;\n\nMATRIX";
      
      
      
      //compute space alignment for different taxa names
      // two passes - one to locate the max length and 2nd to compute padding lengths
      
      unsigned long max_length = 0UL;
      
      for (unsigned long i=0UL; i<sequence_count; i++) {
        StoreIfGreater (max_length, GetSequenceName(i)->sLength);
      }
      
      _SimpleList taxaNamesPadding;
      
      for (unsigned long i=0UL; i<sequence_count; i++) {
        taxaNamesPadding <<  max_length - GetSequenceName(i)->sLength;
      }
      
      
      if (outputFormat%2==0) { // sequential
        for (unsigned long i=0UL; i< sequence_count; i++) {
          if (outputFormat == 4) { // labels
            write_here << "\n\t'"
            << GetSequenceName(i)
            << '\''
            << _String (" ", taxaNamesPadding (i));
            
            
          } else {
            write_here << kStringFileWrapperNewLine;
          }
          write_here << ' ';
          for (long site_index = 0UL; site_index < site_count; site_index++) {
            write_here << (*theData)(theOriginalOrder.lData[site_index],theNodeMap.lData[i],1);
          }
        }
      } else {
        long  sitesDone = 0, upTo;
        
        while (sitesDone< site_count) {
          upTo = sitesDone+printWidth;
          
          if (upTo>site_count) {
            upTo = site_count;
          }
          
          
          for (unsigned long i=0UL; i< sequence_count; i++) {
            if (outputFormat == 5) { // labels
              write_here << "\n\t'"
              << GetSequenceName(i)
              << '"'
              << _String (" ", taxaNamesPadding (i));
            } else {
              write_here << kStringFileWrapperNewLine;
            }
            
            write_here << ' ';
            for (long site_index = sitesDone; site_index < upTo; site_index++) {
              write_here << (*theData)(theOriginalOrder.lData[site_index],theNodeMap.lData[i],1);
            }
            
          }
          write_here << kStringFileWrapperNewLine << kStringFileWrapperNewLine;
          sitesDone = upTo;
        }
        
      }
      write_here << ";\nEND;";
      break;
    }
      
    case 8: {
      for (unsigned long i = 0UL; i< sequence_count; i++) {
        write_here << (*theData)(theOriginalOrder(0),theNodeMap(i),1);
        for (unsigned long j = 1UL; j<site_count; j++) {
          write_here << ',' << (*theData)(theOriginalOrder(j),theNodeMap(i),1);
        }
        write_here << kStringFileWrapperNewLine;
      }
      break;
    }
      
    default: { // hash-mark sequential
      char seqDelimiter = (outputFormat==9)?'>':'#';
      
      for (unsigned long i = 0UL; i< sequence_count; i++) {
        write_here << seqDelimiter << GetSequenceName(i);
        for (unsigned long j = 0UL; j<site_count; j++) {
          if (j % printWidth == 0) {
            write_here << kStringFileWrapperNewLine;
          }
          write_here << (*theData)(theOriginalOrder(j),theNodeMap(i),1);
        }
        write_here << kStringFileWrapperNewLine;
      }
      
      // finally see if we need to write out a tree
      
    }
  }
  
  if (outputFormat != 8) {
    _Parameter  treeDefined;
    checkParameter (dataFileTree, treeDefined,0.0);
    if (treeDefined) {
      _Variable *treeVar = FetchVar(LocateVarByName (dataFileTreeString));
      if (treeVar) {
        _String* treeString = (_String*)(treeVar->Compute())->toStr();
        switch (outputFormat) {
          case 0:
          case 1:
          case 9:
          case 10: {
            write_here << kStringFileWrapperNewLine
            << kStringFileWrapperNewLine
            << *treeString;
            break;
          }
          case 2:
          case 3: {
            write_here << "\n1\n" << *treeString;
            break;
          }
          default: {
            write_here << "\n\nBEGIN TREES;\n\tTREE tree = "
            << *treeString
            << ";\nEND;";
          }
        }
        DeleteObject (treeString);
      }
    }
  }
}

//_________________________________________________________

bool    StoreADataSet (_DataSet* ds, _String* setName) {
    if (!setName->IsValidIdentifier (true)) {
        WarnError (*setName & " is not a valid identifier while constructing a DataSet");
        return false;
    }

    long pos = FindDataSetName (*setName);

    if (pos==-1) {
        dataSetNamesList << setName;
        dataSetList < ds;
    } else {
#if !defined __UNIX__ && ! defined __HEADLESS__
        if (!RequestDataSetReplace (pos)) {
            terminateExecution = true;
            DeleteObject (ds);
            return false;
        }
#endif

        _DataSet* existingDS = (_DataSet*)dataSetList (pos);

        bool isDifferent = existingDS->NoOfSpecies () != ds->NoOfSpecies() ||
                           existingDS->NoOfColumns () != ds->NoOfColumns() ||
                           existingDS->NoOfUniqueColumns () != ds->NoOfUniqueColumns() ||
                           existingDS->GetTT () != ds->GetTT();

        
      
        for (AVLListXLIteratorKeyValue filter_key_value : ObjectIndexer (HY_BL_DATASET_FILTER)) {
          _DataSetFilter * filter = (_DataSetFilter*) filter_key_value.get_object();
          if (filter->GetData() == existingDS) {
            if (isDifferent) {
              ReportWarning (_String("Overwriting dataset '") & *setName & "' caused DataSetFilter " & GetFilterName(filter_key_value.get_index())->Enquote('\'') & " to be deleted");
              DeleteDataFilter(filter_key_value.get_index());
            } else {
              filter->SetData(ds);
            }
          }
        }
      
        dataSetList.Replace(pos,ds,false);
    }
  
    _Parameter normalizeSeqNames = 1.;
    checkParameter (normalizeSequenceNames, normalizeSeqNames, 1.0);
  
    CheckReceptacleAndStore (*setName&".mapping",emptyString,false, new _MathObject, false);
    if (normalizeSeqNames > 0.1) {
      _List _id_mapping;
      _AVLListXL id_mapping (&_id_mapping);
      bool       did_something = false;
      
      for (unsigned long i = 0UL; i < ds->NoOfSpecies(); i ++) {
        _String * old_name = new _String (*ds->GetSequenceName (i));
        if (! old_name->IsValidIdentifier(false) ) {
          ds->GetSequenceName (i)->ConvertToAnIdent(false);
          did_something = true;
        }
        if (id_mapping.Find (ds->GetSequenceName (i)) >= 0) {
          _String new_name (*ds->GetSequenceName (i));
          long suffix = 1L;
          do {
            new_name = *ds->GetSequenceName (i) & "_" & suffix++;
          } while (id_mapping.Find (&new_name) >= 0);
          *ds->GetSequenceName (i) = new_name;
          did_something = true;
        }
        
        ds->GetSequenceName (i)->AddAReference();
        id_mapping.Insert (ds->GetSequenceName (i), (long)old_name, false, false);
      }
      
      if (did_something) {
        _AssociativeList * mapping = new _AssociativeList();
        
        _SimpleList history;
        long t,
             current_index = id_mapping.Traverser(history, t, id_mapping.GetRoot());

        while (current_index >= 0L) {
          mapping->MStore(*(_String*)_id_mapping.GetItem (current_index), *(_String*)id_mapping.GetXtra(current_index));
          current_index = id_mapping.Traverser(history, t);
        }
        
        CheckReceptacleAndStore (*setName&".mapping",emptyString,false, mapping, false);
     }
    }

    CheckReceptacleAndStore (*setName&".species",emptyString,false, new _Constant (ds->NoOfSpecies()), false);
    CheckReceptacleAndStore (*setName&".sites",emptyString,false, new _Constant (ds->NoOfColumns()), false);
    CheckReceptacleAndStore (*setName&".unique_sites",emptyString,false, new _Constant (ds->NoOfUniqueColumns()), false);

    return true;
}

//_________________________________________________________

_Matrix * _DataSet::HarvestFrequencies (unsigned char unit, unsigned char atom, bool posSpec, _SimpleList& hSegmentation, _SimpleList& vSegmentation, bool countGaps) const
{

            
    if (hSegmentation.lLength == 0L || vSegmentation.lLength<unit) { // revert to default (all data)
        if (hSegmentation.lLength==0) {
            hSegmentation.Populate (NoOfSpecies(),0,1);
        }
        if (vSegmentation.lLength<unit) {
            vSegmentation.Clear();
            vSegmentation.Populate (GetNoTypes(),0,1);
        }
    }

    if (unit%atom > 0) { // 20120814 SLKP: changed this behavior to throw errors
        WarnError (_String("Atom should divide unit in ") & _String (__PRETTY_FUNCTION__).Enquote() &" call");
        return new _Matrix (1,1);
    }


    _Matrix   *  out = new _Matrix (ComputePower (theTT->baseLength, atom),
                                    posSpec?unit/atom:1,
                                    false,
                                    true);

    long     positions  =   unit/atom,
             static_store [HYPHY_SITE_DEFAULT_BUFFER_SIZE];
  

    _String unit_for_counting (atom, false);
  
    for (unsigned long site_pattern = 0UL; site_pattern <vSegmentation.lLength;  site_pattern +=unit) { // loop over the set of segments
        // make sure the partition is kosher

        if (site_pattern + unit > vSegmentation.lLength) {
            break;
        }

        for (unsigned long primary_site = site_pattern; primary_site < site_pattern+unit; primary_site += atom) {
          
            long   index_in_pattern = (primary_site-site_pattern)/atom;
            
            for (unsigned long sequence_index = 0; sequence_index <hSegmentation.lLength; sequence_index ++) {
                // loop down each column
              
                unsigned long mapped_sequence_index = hSegmentation.lData[sequence_index];
                // build atomic probabilities
              
                for (unsigned long m = 0; m<atom; m++ ) {
                  unit_for_counting.setChar(m, (*this)(vSegmentation.lData[primary_site+m],mapped_sequence_index,atom));
                }
              
                long resolution_count = theTT->MultiTokenResolutions(unit_for_counting, static_store, countGaps);
              
                if (resolution_count > 0UL) {
          
                  _Parameter    normalized = 1./resolution_count;

                  for (long resolution_index = 0UL; resolution_index < resolution_count; resolution_index ++) {
                    out->theData[posSpec? static_store[resolution_index]*positions+index_in_pattern: static_store[resolution_index]] += normalized;
                  }
                }
            }
        }
    }

     //scale the matrix now

    unsigned long row_count    = out->GetHDim(),
                  column_count = out->GetVDim();
  
    for (unsigned long column =0UL; column < column_count; column++) { // normalize each _column_ to sum to 1.
        _Parameter sum = 0.0;

        for (unsigned long row = 0UL; row < row_count; row++) {
          sum += out->theData [row*column_count + column];
        }

        for (unsigned long row = 0UL; row < row_count; row++) {
          out->theData [row*column_count + column] /= sum;
        }
    }


    return out;
}



