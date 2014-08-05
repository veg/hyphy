/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2009
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon              (apoon@cfenet.ubc.ca)

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

#include "hy_strings.h"
#include "hy_globals.h"
#include "hy_list_numeric.h"
#include "hy_list_reference.h"



#ifndef __HYPHYXCODE__
#include "gnuregex.h"
#else
#include "regex.h"
#include <unistd.h>
#endif


#ifdef __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>



char _hyStringDefaultReturn = 0;
unsigned long _String::storageIncrement = 32;

/*extern int _hy_mpi_node_rank;
long  loopCount      = 3;
char* addrBreak      = 0x583088;*/

struct _hyValidIDCharsType {
  bool valid_chars[256];
  _hyValidIDCharsType(void) {
    for (int c = 0; c < 256; c++) {
      valid_chars[c] = false;
    }
    {
      for (unsigned char c = 'a'; c <= 'z'; c++) {
        valid_chars[c] = true;
      }
    }
    {
      for (unsigned char c = 'A'; c <= 'Z'; c++) {
        valid_chars[c] = true;
      }
    }
    {
      for (unsigned char c = '0'; c <= '9'; c++) {
        valid_chars[c] = true;
      }
    }
    valid_chars[(unsigned char) '_'] = true;
  }
} _hyValidIDChars;

/*
==============================================================
Constructors/Initializers/Copiers/Destructors
==============================================================
*/

//Does nothing
_String::_String(void) {
  Initialize();
}

void _String::Initialize (bool) {
  sLength = 0UL;
  sData = nil;
}

//Length constructor
_String::_String(const unsigned long sL) {

    sLength = sL;
    sData = (char*)checkPointer(MemAllocate(sL + 1UL));
    memset (sData, 0, sL + 1UL);
}

//Length constructor
_String::_String(long sL) {
  char s[32];
  snprintf(s, sizeof(s), "%ld", sL);
  for (sLength = 0; s[sLength]; sLength++)
    ;
  checkPointer(sData = (char *)MemAllocate(sLength + 1));
  memcpy(sData, s, sLength + 1);
}

_String::_String(const _String &source, long from, long to) {
  if (source.sLength) {
    
    long requested_range = source.NormalizeRange(from, to);
    
    if (requested_range > 0L) {
      sLength = requested_range;
      sData = (char *)MemAllocate(sLength + 1UL);

      if (sLength > 32UL) {
        memcpy(sData, source.sData + from, sLength);
      } else
        for (unsigned long k = 0UL; k < sLength; k++) {
          sData[k] = source.sData[k + from];
        }

      sData[sLength] = 0;
      return;
      
    }
  } 
    
  sLength = 0UL;
  sData = (char *)MemAllocate(1UL);
  sData[0] = 0;

}

//Stack copy contructor
_String::_String(const _String &s) { Duplicate((BaseRef) & s); }

_String::_String(_String *s) { CopyDynamicString(s, false); }

//Data constructor
_String::_String(const char *s) {
  // room for the null terminator
  sLength = strlen(s);
  sData = (char *)MemAllocate(sLength + 1UL);
  memcpy(sData, s, sLength + 1UL);
}


//Data constructor
_String::_String(const wchar_t *wc) {
    unsigned long allocated = wcslen (wc);
    sLength = 0;
    checkPointer(sData = (char *)MemAllocate(allocated + 1));
    for (unsigned long cid = 0; cid < allocated; cid ++) {
        int this_char = wctob (wc[cid]);
        if (this_char != WEOF) {
            sData[sLength++] = (char) this_char;
        }
    }
    if (sLength != allocated) {
        sData = (char *)MemReallocate((char *)sData, (sLength+1) * sizeof(char));
    }
    sData[sLength] = 0;
}


//Data constructor
_String::_String(const char s) {
  sLength = 1UL;
  sData = (char *)MemAllocate(2UL);
  sData[0] = s;
  sData[1] = 0;
}

//Data constructor
_String::_String(_Parameter val, const char *format) {
  char s_val[128];
  sLength = snprintf(s_val, 128, format ? format : PRINTF_FORMAT_STRING, val);
  sData = (char *)MemAllocate(sLength + 1UL);
  for (unsigned long k = 0; k <= sLength; k++) {
    sData[k] = s_val[k];
  }
}

//READ FROM FILE
_String::_String(FILE *F) {
  sLength = 0;
  sData = nil;
  if (F) {
    fseek(F, 0, SEEK_END);
    sLength = (unsigned long) ftell(F);
    sData = (char *)MemAllocate(sLength + 1);
    rewind(F);
    fread(sData, 1, sLength, F);
    sData[sLength] = 0;
  }
}

//Destructor
_String::~_String(void) {
  if (sData) {
    free(sData);
  }
}

void _String::CopyDynamicString(_String *s, bool flushMe) {
  if (flushMe && sData) {
    free(sData);
  }
  sLength = s->sLength;
  if (s->SingleReference ()) {
    sData = s->sData;
    s->sData = nil;
    DeleteObject(s);
  } else {
    checkPointer(sData = (char *)MemAllocate(sLength + 1UL));
    if (s->sData) {
      memcpy(sData, s->sData, sLength + 1);
    } else {
      sData[0] = 0;
    }
    s->RemoveAReference();
  }
}

void _String::Duplicate(BaseRefConst ref) {
  _String const *s = (_String const*)ref;
  sLength = s->sLength;

  if (s->sData) {
    sData = (char *)MemAllocate(sLength + 1UL);
    memcpy(sData, s->sData, sLength + 1UL);
  } else {
    sData = nil;
  }

}

void _String::DuplicateErasing(BaseRefConst ref) {
  if (sData) {
    free(sData);
  }
  Duplicate(ref);

}

//Make dynamic copy
BaseRef _String::makeDynamic(void) const{
  return new _String (*this);
}


/*
 ==============================================================
 Accessors
 ==============================================================
*/

// String length
unsigned long _String::Length(void) const { return sLength; }

void _String::setChar(unsigned long index, char c) {
  if (index < sLength) {
    sData[index] = c;
  }
}

/*
 ==============================================================
 Represent as STR
 ==============================================================
*/

BaseRef _String::toStr(void) const{
  return new _String (*this);
}

/*
==============================================================
Operator Overloads
==============================================================
*/

// Element location functions
char &_String::operator[](long index) {
  if (((unsigned long) index) < sLength) {
    return sData[index];
  }
  return _hyStringDefaultReturn;
}

//Element location functions
char _String::operator()(long index) const {
  return getChar(index);
}

  //Element location functions
const char _String::getChar(long index) const {
  if (index >= 0L && index < sLength) {
    return sData[index];
  }
  return _hyStringDefaultReturn;
}


// Assignment operator
void _String::operator=(const _String& s) {
  if (sData) {
    free(sData);
  }
  Duplicate(&s);
}

_Parameter _String::toNum(void) const{
    if (sLength == 0UL) {
        return 0.;
    }
    char *endP;
    return strtod(sData, &endP);
}



//Append operator
const _String _String::operator&(const _String& s)  const {
  unsigned long combined_length = sLength + s.sLength;

  if (combined_length == 0UL) {
    return empty;
  }

  _String res(combined_length);

  if (sLength) {
    memcpy(res.sData, sData, sLength);
  }

  if (s.sLength) {
    memcpy(res.sData + sLength, s.sData, s.sLength);
  }

  res.sData[res.sLength] = 0;
  return res;
}


//Return good ole char*
_String::operator const char *(void) const { return sData; }


// lexicographic comparisons
bool _String::operator==(const _String& s) const { return Compare (&s) == 0; }
bool _String::operator>(const _String & s) const { return Compare (&s) > 0; }
bool _String::operator<=(const _String & s) const { return Compare (&s) <= 0; }
bool _String::operator>=(const _String & s) const { return Compare (&s) >= 0; }
bool _String::operator!=(const _String & s) const { return Compare (&s) != 0; }
bool _String::operator<(const _String & s) const { return Compare (&s) < 0; }

/*
 ==============================================================
 Lexicographic Comparison Methods
 ==============================================================
 */

bool _String::contains (const _String& s) const { return Find(s) != -1; }
bool _String::contains (char c) const    { return Find(c) != HY_NOT_FOUND; }

char _String::Compare(_String const * s) const {
  unsigned long upTo = MIN (sLength, s->sLength);
  
  for (unsigned long i = 0UL; i < upTo; i++) {
    if (sData[i] < s->sData[i]) {
      return -1;
    }
    if (sData[i] > s->sData[i]) {
      return 1;
    }
  }
  
  if (sLength == s->sLength) {
    return 0;
  }
  
  return sLength < s->sLength ? -1 : 1;
}

bool _String::Equal(const _String *s) const {
  return Compare (s) == 0;
}

bool _String::Equal(const char c) const {
  return sLength == 1 && sData[0] == c;
}



bool _String::EqualWithWildChar(const _String &s, const char wildchar) const {
  // wildcards only matter in the second string
  
  if (s.sLength > 0UL && wildchar != '\0') {
    unsigned long   match_this_char = 0UL;
    bool            is_wildcard = s.sData[match_this_char] == wildchar,
                    scanning_s = is_wildcard;
    
    unsigned long i = 0UL;
    while (1) {
      
      if (scanning_s) { // skip consecutive wildcards in "s"
        scanning_s = s.sData[++match_this_char] == wildchar;
      } else {
        if (sData[i] == s.sData[match_this_char]) {
              // try character match
              // note that the terminal '0' characters will always match, so
              // this is where we terminate
          i++;
          match_this_char++;
          if (i > sLength || match_this_char > s.sLength) {
            break;
          }
          is_wildcard =  s.sData[match_this_char] == wildchar;
          scanning_s = is_wildcard;
        } else { // match wildcard
          if (!is_wildcard) {
            return false;
          }
          scanning_s = false;
          i++;
        }
      }
    }
    
    return match_this_char > s.sLength;
  } else {
    return sLength == 0UL;
  }
  
  return false;
}

  //Begins with string
bool _String::startswith(const _String& s, bool caseSensitive) const{
  return (caseSensitive ? Find (s, 0L, s.sLength - 1L)
          : FindAnyCase (s, 0L, s.sLength -1L)) == 0L;
}

bool _String:: endswith(const _String& s, bool caseSensitive) const{
  if (sLength >= s.sLength) {
    return (caseSensitive ? Find (s, sLength - s.sLength)
            : FindAnyCase (s, sLength - s.sLength)) == sLength - s.sLength;
  }
  return false;
}

/*
==============================================================
Manipulators
==============================================================
*/

// a convenicence function not to write (const char*)(*this)[]

const char *_String::getStr(void) const { return sData; }


long  _String::NormalizeRange(long & from, long & to) const {

  if (sLength == 0UL) {
    return 0L;
  }

  if (from < 0L) {
    from = 0L;
  }
  
  if (to < 0L || to >= sLength ) {
    to = sLength - 1UL;
  }
  
  return to - from + 1L;
  
}

//s[0]...s[sLength-1] => s[sLength-1]...s[0]
void _String::Flip(void) {
  for (unsigned long i = 0UL; i < (sLength >> 1); i++) {
    char c;
    SWAP  (sData[i], sData[sLength - 1 - i], c);
  }
}

const _String _String::Chop(long from, long to) const {
  
  long resulting_length = NormalizeRange(from,to);
  
  if (resulting_length > 0L) {
    _String res((unsigned long)(sLength - resulting_length));
    if (from > 0L) {
      memcpy(res.sData, sData, from);
    }
    if (to + 1L < sLength) {
      memcpy(res.sData + from, sData + to + 1L, sLength - to - 1L);
    }
    
    return res;
  }
  
  return *this;
 
}

//Cut string from, to (-1 for any means from beginning/to end)
const _String _String::Cut(long from, long to) const {
  return _String (*this, from, to);
}

//Delete range char operator
void _String::Delete(long from, long to) {
  long resulting_length = NormalizeRange(from,to);

  if (resulting_length > 0L) {
    if (to < (long)sLength - 1UL) {
      memmove(sData + from, sData + to + 1, sLength - to - 1);
    }
    sLength -= resulting_length;
    sData = MemReallocate(sData, sizeof(char) * (sLength + 1UL));
    sData[sLength] = 0;
  }
}

//Insert a single char operator
void _String::Insert(char c, long pos) {
  if (pos < 0L || pos >= sLength) {
    pos = sLength;
  }

  sData = MemReallocate(sData, sizeof(char) * (sLength + 2));

  if (pos < sLength) {
    memmove(sData + pos + 1, sData + pos, sLength - pos);
  }

  sData[pos] = c;
  sData[++sLength] = 0;
}

//Cut string from, to (-1 for any means from beginning/to end)
void _String::Trim(long from, long to) {
  
  long resulting_length = NormalizeRange(from, to);
  
  if (resulting_length > 0L) {
      if (from > 0L) {
        memmove(sData, sData + from, resulting_length);
      }
      sLength = resulting_length;
      sData = MemReallocate(sData, resulting_length + 1UL);
      sData[resulting_length] = 0;
  } else {
      sLength = 0UL;
      sData = MemReallocate(sData, 1UL);
      sData[0] = 0;
  }
}


void _String::StripQuotes(char open_char, char close_char) {
  if (sLength >= 2UL) {
    if (getChar(0) == open_char && getChar (sLength - 1UL) == close_char) {
      Trim (1, sLength - 2UL);
    }
  }
}

/*
==============================================================
Case Methods
==============================================================
*/

const _String _String::UpCase(void) const {
  _String res (sLength);
  for (unsigned long i = 0UL; i < sLength; i++) {
    res.sData[i] = toupper(sData[i]);
  }
  return res;
}

const _String _String::LoCase(void) const {
  _String res (sLength);
  for (unsigned long i = 0UL; i < sLength; i++) {
    res.sData[i] = tolower(sData[i]);
  }
  return res;
}

/*
==============================================================
Search and replace
==============================================================
*/

long _String::Find(const _String& s, long from, long to) const {

  if (s.sLength) {
    long span = NormalizeRange(from, to);
    if (span >= (long) s.sLength) {
      const unsigned long upper_bound = to - s.sLength + 2L;
      for (unsigned long i = from; i < upper_bound ; i++) {
        unsigned long j = 0UL;
        for (;j < s.sLength; j++) {
          if (sData[i + j] != s.sData[j]) {
            break;
          } 
        }
        if (j == s.sLength) {
          return i;
        }
      }
    }
  }
  return HY_NOT_FOUND;
}

  //Find first occurence of the string between from and to
long _String::Find(const char s, long from, long to) const {
  if (sLength) {
    long span = NormalizeRange(from, to);
    if (span > 0L) {
      char sentinel = sData[to+1];
      sData[to+1] = s;
      long index = from;
      while (sData[index] != s) {
        index++;
      }
      sData[to+1] = sentinel;
      return index <= to ? index : HY_NOT_FOUND;
    }
  }
  
  return HY_NOT_FOUND;
}


//Find first occurence of the string between from and to
long _String::FindBackwards(const _String & s, long from, long to) const {

  if (s.sLength) {
    long span = NormalizeRange(from, to);
    if (span >= (long) s.sLength) {
      const long upper_bound = to - s.sLength + 1L;
      for (long i = upper_bound; i >= from; i--) {
        unsigned long j = 0UL;
        for (;j < s.sLength; j++) {
          if (sData[i + j] != s.sData[j]) {
            break;
          } 
        }
        if (j == s.sLength) {
          return i;
        }
      }
    }
  }
  return HY_NOT_FOUND;
}


// find first occurence of the string between from and to
// case insensitive
long _String::FindAnyCase(const _String& s, long from, long to) const {
    return UpCase().Find (s.UpCase(), from, to);
}

// find first occurence of the string between from and to
bool _String::ContainsSubstring(const _String &s) const {
  return Find (s) != HY_NOT_FOUND;
}


//Replace string 1 with string 2, all occurences true/false
const _String _String::Replace(const _String& s, const _String& d, bool replace_all) const {

  if (sLength < s.sLength || s.sLength == 0UL) {
    return *this;
  }

  _StringBuffer replacementBuffer;
  unsigned long anchor_index = 0UL;
  for (; anchor_index <= sLength - s.sLength; anchor_index ++) {
    unsigned long search_index = 0UL;
    for (; search_index < s.sLength; search_index++) {
      if (sData[anchor_index + search_index] != s.sData[search_index]) {
        break;
      }
    }

    if (search_index == s.sLength) {
      replacementBuffer << d;
      anchor_index += s.sLength - 1UL;
      if (replace_all == false) {
        anchor_index ++;
        break;
      }
    } else {
      replacementBuffer << sData[anchor_index];
    }
  }
  
  replacementBuffer.AppendSubstring(*this, anchor_index, HY_NOT_FOUND);
  return replacementBuffer;
}

const _List _String::Tokenize(const _String& splitter) const {
  _List tokenized;
  
  long cp = 0L, cpp;
  
  while ((cpp = Find(splitter, cp, HY_NOT_FOUND)) != HY_NOT_FOUND) {
    if (cpp > cp) {
      tokenized.append (new _String(*this, cp, cpp - 1L));
    } else {
      tokenized.append (new _String);
    }
    
    cp = cpp + splitter.sLength;
  }
  
  tokenized.append(new _String(*this, cp, HY_NOT_FOUND));
  return tokenized;
}

/*
==============================================================
Formatters
==============================================================
*/



// Format second difference as HHH..H:MM:SS 
const _String _String::FormatTimeString(long time_diff){

  long fields [3] = {time_diff / 3600L, time_diff / 60L % 60L, time_diff % 60L};
  
  _StringBuffer time_string;
  
  for (unsigned long l = 0; l < 3UL; l++) {
    if (l) {
      time_string << ':';
    }
    if (fields[l] < 10L) {
      time_string << '0';
    }
    time_string << _String (fields[l]);
  }
  
  return time_string;
}


/*
 ==============================================================
 Utility functions with no clear category
 ==============================================================
 */



const _String _String::Random(const unsigned long length, const _String *alphabet) {
  _StringBuffer random(length + 1UL);
  
  unsigned int alphabet_length = alphabet ? alphabet->sLength : 127;
  
  if (length > 0UL && alphabet_length > 0UL) {
    for (unsigned long c = 0UL; c < length; c++) {
      unsigned long idx = genrand_int32() % alphabet_length;
      if (alphabet) {
        random << alphabet->sData[idx];
      } else {
        random << (char)(1UL + idx);
      }
    }
  }
  
  return random;
}

long _String::Adler32(void) const {
  
  unsigned long len = sLength,
                a = 1UL,
                b = 0UL,
                i = 0UL;
  
  while (len) {
    unsigned long tlen = len > 5550UL ? 5550UL : len;
    len -= tlen;
    do {
      a += sData[i++];
      b += a;
    } while (--tlen);
    a = (a & 0xffff) + (a >> 16) * (65536UL - HY_STRING_MOD_ADLER);
    b = (b & 0xffff) + (b >> 16) * (65536UL - HY_STRING_MOD_ADLER);
  }
  
  if (a >= HY_STRING_MOD_ADLER) {
    a -= HY_STRING_MOD_ADLER;
  }
  
  b = (b & 0xffff) + (b >> 16) * (65536UL - HY_STRING_MOD_ADLER);
  
  if (b >= HY_STRING_MOD_ADLER) {
    b -= HY_STRING_MOD_ADLER;
  }
  
  return b << 16 | a;
}

_String const _String::Sort(_SimpleList *index) const {
    
    if (index) {
        index->Clear();
    }
    
    if (sLength > 0UL) {
        _hyListOrderable<char> sorted;
        
        if (index) {
            for (unsigned long i = 0UL; i < sLength; i++) {
                sorted << sData[i];
                (*index) << i;
            }
            sorted.Sort (true, index);
        } else {
            for (unsigned long i = 0; i < sLength; i++) {
                sorted << sData[i];
            }
            sorted.Sort();
        }
        _String result (sLength);
        
        for (unsigned long i = 0UL; i < sLength; i++) {
            result.setChar (i, sorted.AtIndex(i));
        }
        
        return result;
    }
    
    return empty;
}

/*
 ==============================================================
 Identifier Methods
 ==============================================================
 */

bool _String::IsValidIdentifier(bool allow_compounds) const {
    // 201407
    
    return sLength > 0UL && _IsValidIdentifierAux (allow_compounds) == Length() - 1UL;

#ifndef HY_2014_REWRITE_MASK
    // TO DO -- MOVE THIS CHECK ELSEWHERE?
    return hyReservedWords.Find(this) == HY_NOT_FOUND;
#endif
    
    
}

long _String::_IsValidIdentifierAux(bool allow_compounds, char wildcard) const {
    
    unsigned long current_index = 0UL;
    
    bool          first     = true;
    
    for (; current_index < sLength; current_index ++) {
        char current_char = getChar (current_index);
        if (first) {
            if ( ! (isalpha (current_char) || current_char == '_')) {
                break;
            }
            first = false;
        } else {
            if ( ! (isalnum (current_char) || current_char == '_' || current_char == wildcard)) {
                if (allow_compounds && current_char == '.') {
                    first = true;
                } else {
                    break;
                }
            }
        }
    }
    
    if (current_index) {
        return current_index - 1UL;
    }
    
    return HY_NOT_FOUND;
}

const _String _String::ShortenVarID(_String const &containerID) const {
    if (startswith(containerID)) {
        unsigned long prefix_length = containerID.Length();
        
        if (getChar(prefix_length) == '.' && sLength > prefix_length + 1L) {
            return Cut(prefix_length + 1L, HY_NOT_FOUND);
        }
    }
    return *this;
}

/*
 ==============================================================
 Space Methods
 ==============================================================
 */

//Replace all space runs with a single space
const _String _String::CompressSpaces(void) const {
    _StringBuffer temp(sLength + 1UL);
    bool skipping = false;
    
    for (unsigned long k = 0UL; k < sLength; k++) {
        if (!isspace(sData[k])) {
            temp << sData[k];
            skipping = false;
        } else {
            if (!skipping) {
                skipping = true;
                temp << ' ';
            }
        }
    }
    return temp;
}

//Remove all spaces
const _String _String::KillSpaces(void) const {
    _StringBuffer temp(sLength + 1UL);
    for (unsigned long k = 0UL; k < sLength; k++) {
        if (!isspace(sData[k])) {
            temp << sData[k];
        }
    }
    return temp;
}



  //Locate the first non-space charachter of the string
long _String::FirstNonSpaceIndex(long from, long to, unsigned char direction) const {
  
  long requested_range = NormalizeRange(from, to);
  
  if (requested_range > 0L) {
    if (direction == HY_STRING_DIRECTION_FORWARD) {
      for (; from <= to; from++) {
        if (!isspace(sData[from])) {
          return from;
        }
      }
    } else {
      for (; to>=from; to--) {
        if (!isspace(sData[to])) {
          return to;
        }
      }
    }
  }
  
  return HY_NOT_FOUND;
}


  //Locate the first non-space charachter of the string
long _String::FirstSpaceIndex(long from, long to, unsigned char direction) const {
  long requested_range = NormalizeRange(from, to);
  
  if (requested_range > 0L) {
    if (direction == HY_STRING_DIRECTION_FORWARD) {
      for (; from <= to; from++) {
        if (isspace(sData[from])) {
          return from;
        }
      }
    } else {
      for (; to>=from; to--) {
        if (isspace(sData[to])) {
          return to;
        }
      }
    }
  }
  
  return HY_NOT_FOUND;
}

  //Locate the first non-space charachter of the string
char _String::FirstNonSpace(long start, long end, unsigned char direction) const {
  long r = FirstNonSpaceIndex(start, end, direction);
  return r == HY_NOT_FOUND ? _hyStringDefaultReturn : sData[r];
}

/*!!!!!!!!!!!!!!!!!!
 
 DONE UP TO HERE
 
 !!!!!!!!!!!!!!!!!!!!*/


bool _String::IsValidRefIdentifier(void) const {
    if (sLength < 2) {
        return false;
    }
    if (sData[sLength - 1] == '&') {
        return Cut(0, sLength - 2).IsValidIdentifier();
    }
    return false;
}

//Convert a string to a valid ident
void _String::ConvertToAnIdent(bool strict) {
    _StringBuffer *result = new _StringBuffer((unsigned long) sLength + 1UL);
    
    if (sLength) {
        if (strict) {
            if (((sData[0] >= 'a') && (sData[0] <= 'z')) ||
                ((sData[0] >= 'A') && (sData[0] <= 'Z')) || (sData[0] == '_')) {
                (*result) << sData[0];
            } else {
                (*result) << '_';
            }
        } else {
            if (((sData[0] >= 'a') && (sData[0] <= 'z')) ||
                ((sData[0] >= 'A') && (sData[0] <= 'Z')) || (sData[0] == '_') ||
                ((sData[0] >= '0') && (sData[0] <= '9'))) {
                (*result) << sData[0];
            } else {
                (*result) << '_';
            }
        }
        
        long l = 0;
        for (unsigned long k = 1UL; k < sLength; k++) {
            unsigned char c = sData[k];
            if (_hyValidIDChars.valid_chars[c]) {
                (*result) << c;
                l++;
            } else if (result->sData[l] != '_') {
                (*result) << '_';
                l++;
            }
        }
    }
    
    CopyDynamicString(result, true);
}





long _String::FindEndOfIdent(long start, long end, char wild) {
  if (sLength == 0) {
    return HY_NOT_FOUND;
  }

  if (start == -1) {
    start = ((long) sLength) - 1;
  }
  if (end == -1) {
    end = ((long) sLength) - 1;
  }

  long i = start;

  for (; i <= end; i++)
    if (!(isalnum(sData[i]) || sData[i] == '.' || sData[i] == wild ||
          sData[i] == '_')) {
      break;
    }

  if (i > start + 2 && sData[i - 1] == '_' && sData[i - 2] == '_') {
    return i - 3;
  }

  return i - 1;
}
 

long _String::ExtractEnclosedExpression(long &from, char open, char close,
                                        bool respectQuote, bool respectEscape) {
  long currentPosition = from, currentLevel = 0;

  bool isQuote = false, doEscape = false;

  while (currentPosition < sLength) {
    char thisChar = sData[currentPosition];

    if (!doEscape) {
      if (thisChar == '"' && respectQuote && !doEscape) {
        isQuote = !isQuote;
      } else if (thisChar == open && !isQuote) {
        // handle the case when close and open are the same
        if (currentLevel == 1 && open == close && from < currentPosition) {
          return currentPosition;
        }
        currentLevel++;
        if (currentLevel == 1) {
          from = currentPosition;
        }
      } else if (thisChar == close && !isQuote) {
        currentLevel--;
        if (currentLevel == 0 && from < currentPosition) {
          return currentPosition;
        }
        if (currentLevel < 0) {
          return HY_NOT_FOUND;
        }
      } else if (thisChar == '\\' && respectEscape && isQuote && !doEscape) {
        doEscape = true;
      }
    } else {
      doEscape = false;
    }

    currentPosition++;
  }

  return HY_NOT_FOUND;
}




long _String::FindTerminator(long from, _String &terminators) {
  long currentPosition = from, currentCurly = 0, currentSquare = 0,
       currentParen = 0;

  bool isQuote = false, doEscape = false;

  while (currentPosition < sLength) {
    char thisChar = sData[currentPosition];
    if (!doEscape) {
      if (thisChar == '"' && !doEscape) {
        isQuote = !isQuote;
      } else {
        if (!isQuote) {
          if (thisChar == '{') {
            currentCurly++;
          } else if (thisChar == '[') {
            currentSquare++;
          } else if (thisChar == '(') {
            currentParen++;
          }
          if (currentCurly > 0 && thisChar == '}') {
            currentCurly--;
          } else if (currentSquare > 0 && thisChar == ']') {
            currentSquare--;
          } else if (currentParen > 0 && thisChar == ')') {
            currentParen--;
          } else if (currentParen == 0 && currentSquare == 0 &&
                     currentCurly == 0)
            for (long s = 0; s < terminators.sLength; s++)
              if (thisChar == terminators.sData[s]) {
                return currentPosition;
              }
        } else {
          if (thisChar == '\\' && isQuote && !doEscape) {
            doEscape = true;
          }
        }
      }
    } else {
      doEscape = false;
    }

    currentPosition++;
  }

  return HY_NOT_FOUND;
}








long _String::LempelZivProductionHistory(_SimpleList *rec) {
  if (rec) {
    rec->Clear();
  }

  if (sLength == 0) {
    return 0;
  }

  if (rec) {
    (*rec) << 0;
  }

  long cp = 1, pH = 1;

  while (cp < sLength) {
    long maxExtension = 0;

    for (long ip = 0; ip < cp; ip++) {
      long sp = ip, mp = cp;

      while ((mp < sLength) && (sData[mp] == sData[sp])) {
        mp++;
        sp++;
      }

      if (mp == sLength) {
        maxExtension = sLength - cp;
        break;
      } else {
        if ((mp = mp - cp + 1) > maxExtension) {
          maxExtension = mp;
        }
      }
    }

    cp = cp + maxExtension;
    if (rec) {
      (*rec) << cp - 1;
    } else {
      pH++;
    }
  }

  if (rec) {
    return rec->countitems();
  }

  return pH;
}











_Parameter _String::ProcessTreeBranchLength(void) {
  _Parameter res = -1.;

  if (sLength) {
    if (sData[0] == ':') {
      res = Cut(1, -1).toNum();
    } else {
      res = toNum();
    }

    if (res < 1e-10) {
      res = 1e-10;
    }
  }

  return res;
}

bool _String::IsALiteralArgument(bool stripQuotes) {
  if (sLength >= 2) {
    long from = 0, to = ExtractEnclosedExpression(from, '"', '"', false, true);

    if (from == 0 && to == sLength - 1) {
      if (stripQuotes) {
        Trim(1, sLength - 2);
      }
      return true;
    }
  }
  return false;
}

//TODO: This is a global function.
bool hyIDValidator(_String *s) { return s->IsValidIdentifier(false); }





void _String::ProcessParameter(void) {
#ifndef HY_2014_REWRITE_MASK
  if (Equal(&getDString)) {
    *this = ReturnDialogInput();
  }
#endif
}

//==============================================================
//Filename and Platform Methods
//==============================================================

bool _String::ProcessFileName(bool isWrite, bool acceptStringVars, Ptr theP,
                              bool assume_platform_specific,
                              _ExecutionList *caller) {
#ifndef HY_2014_REWRITE_MASK
  _String errMsg;

  try {
    if (Equal(&getFString) || Equal(&tempFString)) { // prompt user for file
      if (Equal(&tempFString)) {
#if not defined __MINGW32__ &&not defined __WINDOZE__
#ifdef __MAC__
        char tmpFileName[] = "HYPHY-XXXXXX";
#else
        char tmpFileName[] = "/tmp/HYPHY-XXXXXX";
#endif

        int fileDescriptor = mkstemp(tmpFileName);
        if (fileDescriptor == -1) {
          throw("Failed to create a temporary file name");
        }
        *this = tmpFileName;
        CheckReceptacleAndStore(&useLastFString, empty, false,
                                new _FString(*this, false), false);
        close(fileDescriptor);
        return true;
#else
        throw(tempFString & " is not implemented for this platform");
#endif
      } else {
        if (!isWrite) {
          *this = ReturnFileDialogInput();
        } else {
          *this = WriteFileDialogInput();
        }
      }
      ProcessFileName(false, false, theP,
#if defined __MAC__ || defined __WINDOZE__
                      true
#else
                      false
#endif
                      ,
                      caller);

      CheckReceptacleAndStore(&useLastFString, empty, false,
                              new _FString(*this, false), false);
      return true;
    }

    if (acceptStringVars) {
      *this = ProcessLiteralArgument(this, (_VariableContainer *)theP, caller);
      if (caller && caller->IsErrorState()) {
        return false;
      }
    } else {
      StripQuotes();
    }

    if (!sLength) {
      return true;
    }
  }

  catch (_String errmsg) {
    if (caller) {
      caller->ReportAnExecutionError(errMsg);
    } else {
      WarnError(errMsg);
    }
    return false;
  }

#if (defined __UNIX__ || defined __HYPHY_GTK__) && !defined __MINGW32__
  //UNIX LINES HERE
  if (Find('\\') != -1) { // DOS (ASSUME RELATIVE) PATH
    *this = Replace("\\", "/", true);
  } else if (Find(':') != -1) { // Mac (Assume Relative) PATH
    *this = Replace("::", ":../", true);
    if (getChar(0) == ':') {
      Trim(1, -1);
    }
    *this = Replace(':', '/', true);
  }

  if (getChar(0) != '/') { // relative path
    if (pathNames.lLength) {
      _String *lastPath = (_String *)pathNames(pathNames.lLength - 1);
      long f = lastPath->sLength - 2, k = 0;

      // check the last stored absolute path and reprocess this relative path
      // into an absolute.
      while (beginswith("../")) {
        if ((f = lastPath->FindBackwards('/', 0, f) - 1) == -1) {
          return true;
        }
        Trim(3, -1);
        k++;
      }
      if (k == 0) {
        *this = *lastPath & (*this);
      } else {
        *this = lastPath->Cut(0, f + 1) & (*this);
      }
    }
  }
#endif

#if defined __WINDOZE__ || defined __MINGW32__ // WIN/DOS code
  if (Find('/') != -1) {                       // UNIX PATH
    if (getChar(0) == '/') {
      Trim(1, -1);
    }
    *this = Replace("/", "\\", true);
  } else {
    if (Find('\\') == -1) {
      // check to see if this is a relative path
      *this = Replace("::", ":..\\", true);
      if ((sData[0] == ':')) {
        Trim(1, -1);
      }
      *this = Replace(':', '\\', true);
    }
  }

  if (Find(':') == -1 && Find("\\\\", 0, 1) == -1) { // relative path

    if (pathNames.lLength) {
      _String *lastPath = (_String *)pathNames(pathNames.lLength - 1);
      long f = lastPath->sLength - 2, k = 0;
      // check the last stored absolute path and reprocess this relative path
      // into an absolute.
      while (beginswith("..\\")) {
        f = lastPath->FindBackwards('\\', 0, f) - 1;
        if (f == -1) {
          return false;
        }
        Trim(3, -1);
        k++;
      }
      if (k == 0) {
        if (lastPath->sData[lastPath->sLength - 1] != '\\') {
          *this = *lastPath & '\\' & (*this);
        } else {
          *this = *lastPath & (*this);
        }
      } else {
        *this = lastPath->Cut(0, f + 1) & (*this);
      }
    }

  }

  _String escapedString(sLength, true);
  for (long stringIndex = 0; stringIndex < sLength; stringIndex++) {
    char currentChar = getChar(stringIndex);
    //char b[256];
    //snprintf (b, sizeof(b),"%c %d\n", currentChar, currentChar);
    //BufferToConsole (b);
    switch (currentChar) {
    case '\t':
      escapedString << '\\';
      escapedString << 't';
      break;
    case '\n':
      escapedString << '\\';
      escapedString << 'n';
      break;
    default:
      escapedString << currentChar;
    }
  }
  escapedString.Finalize();
  (*this) = escapedString;

#endif

#ifdef __MAC__
  if (!assume_platform_specific && Find('/') != -1) { // UNIX PATH
    bool rootPath = false;
    if (sData[0] == '/') {
      rootPath = true;
      *this = volumeName & Cut(1, -1);
    }

    if (beginswith("..")) {
      *this = _String('/') & Cut(2, -1);
    }

    *this = Replace("/", ":", true);
    *this = Replace("..", "", true);

    if (sData[0] != ':' && !rootPath) {
      *this = _String(':') & *this;
    }
  } else {
    if (!assume_platform_specific &&
        Find('\\') != -1) { // DOS PATH (ASSUME PARTIAL)
      if (beginswith("..")) {
        *this = _String('\\') & Cut(2, -1);
      }
      *this = Replace("\\", ":", true);
      *this = Replace("..", "", true);
      if (Find(':') != -1) {
        *this = _String(':') & *this;
      }
    } else { // MAC PATH
      if (Find(':') != -1) {
        if (sData[0] != ':') {
          if (!beginswith(volumeName)) {
            if (pathNames.lLength) {
              _String *lastPath = (_String *)pathNames(pathNames.lLength - 1);
              if (!beginswith(lastPath->Cut(0, lastPath->Find(':')))) {
                *this = _String(':') & *this;
              }
            } else {
              *this = _String(':') & *this;
            }
          }
        }
      } else {
        *this = _String(':') & *this;
      }
    }
  }

  if (sData[0] == ':') { // relative path
    long f = -1, k = 0;
    if (pathNames.lLength) {
      _String *lastPath = (_String *)pathNames(pathNames.lLength - 1);
      // check the last stored absolute path and reprocess this relative path
      // into an absolute.
      while (sData[k] == ':') {
        f = lastPath->FindBackwards(':', 0, f) - 1;
        if (f == -1) {
          return;
        }
        k++;
      }
      *this = lastPath->Cut(0, f + 1) & Cut(k, -1);
    } else {
      *this = empty;
    }
  }
#endif
#endif
  return true;
}

//Compose two UNIX paths (abs+rel)
_String _String::PathComposition(_String relPath) {
  if (relPath.sData[0] != '/') { // relative path
    long f = -1, k = 0;
    f = sLength - 2;
    _String result = *this;

    while (relPath.startswith("../")) {

      //Cut Trim relPath
      f = FindBackwards('/', 0, f) - 1;

      relPath = relPath.Chop(0, 2);
      result.Trim(0, f + 1);

      if (f == -1) {
        return empty;
      }
      k++;

    }

    return result & relPath;
  } else {
    return relPath;
  }
  return empty;
}

//Mac only so far
_String _String::PathSubtraction(_String &p2, char) {
  _String result;
  char separator = GetPlatformDirectoryChar();

  //if (pStyle == 0)
  //    separator = ':';
  long k;
  for (k = 0;(k < sLength) && (k < p2.sLength) && (sData[k] == p2.sData[k]);
       k++)
    ;
  if (k > 0) {
    while (sData[k] != separator) {
      k--;
    }
    if (k > 0) {
      long m = k + 1, levels = 0;
      for (; m < sLength; m++)
        if (sData[m] == separator) {
          levels++;
        }
      if (levels) {
        result = separator;
        while (levels) {
          result.Insert(separator, -1);
          levels--;
        }
      }
      result = result & p2.Cut(k + 1, -1);
      return result;
    }
  }
  return empty;
}



/*
==============================================================
Regular Expression Methods
==============================================================
*/

_String GetRegExpError(int error) {
  char buffer[512];
  buffer[regerror(error, nil, buffer, 511)] = 0;
  return _String("Regular Expression error:") & buffer;
}

void FlushRegExp(Ptr regExpP) {
  regex_t *regEx = (regex_t *)regExpP;
  regfree(regEx);
  delete regEx;
}

Ptr PrepRegExp(_String *source, int &errCode, bool caseSensitive) {
  regex_t *res = new regex_t;
  checkPointer(res);

  errCode = regcomp(res, source->sData,
                    REG_EXTENDED | (caseSensitive ? 0 : REG_ICASE));

  if (errCode) {
    FlushRegExp((Ptr) res);
    return nil;
  }
  return (Ptr) res;
}

void _String::RegExpMatch(Ptr pattern, _SimpleList &matchedPairs) {
  if (sLength) {
    regex_t *regEx = (regex_t *)pattern;

    regmatch_t *matches = new regmatch_t[regEx->re_nsub + 1];
    int errNo = regexec(regEx, sData, regEx->re_nsub + 1, matches, 0);
    if (errNo == 0) {
      for (long k = 0; k <= regEx->re_nsub; k++) {
        matchedPairs << matches[k].rm_so;
        matchedPairs << matches[k].rm_eo - 1;
      }
    }
    delete[] matches;
  }
}

void _String::RegExpMatchAll(Ptr pattern, _SimpleList &matchedPairs) {
  if (sLength) {
    regex_t *regEx = (regex_t *)pattern;

    regmatch_t *matches = new regmatch_t[regEx->re_nsub + 1];
    int errNo = regexec(regEx, sData, regEx->re_nsub + 1, matches, 0);
    while (errNo == 0) {
      long offset = matchedPairs.countitems()
                        ? matchedPairs.Element(-1) + 1
                        : 0;

      matchedPairs << matches[0].rm_so + offset;
      matchedPairs << matches[0].rm_eo - 1 + offset;

      offset += matches[0].rm_eo;
      if (offset < sLength) {
        errNo = regexec(regEx, sData + offset, regEx->re_nsub + 1, matches, 0);
      } else {
        break;
      }
    }
    delete[] matches;
  }
}

void _String::RegExpMatchOnce(_String *pattern, _SimpleList &matchedPairs,
                              bool caseSensitive, bool handleErrors) {
  if (sLength) {
    int errNo = 0;
    Ptr regex = PrepRegExp(pattern, errNo, caseSensitive);
    if (regex) {
      RegExpMatch(regex, matchedPairs);
      FlushRegExp(regex);
    } else if (handleErrors) {
      WarnError(GetRegExpError(errNo));
    }
  }
}



unsigned char _String::ProcessVariableReferenceCases(_String &referenced_object,
                                                     _String *context) {
#ifndef HY_2014_REWRITE_MASK
  char first_char = getChar(0);
  bool is_func_ref = getChar(sLength - 1) == '&';

  if (first_char == '*' || first_char == '^') {
    if (is_func_ref) {
      referenced_object = empty;
      return HY_STRING_INVALID_REFERENCE;
    }
    bool is_global_ref = first_char == '^';
    _String choppedVarID(*this, 1, -1);
    if (context) {
      choppedVarID = *context & '.' & choppedVarID;
    }
    _FString *dereferenced_value =
        (_FString *)FetchObjectFromVariableByType(&choppedVarID, STRING);
    if (dereferenced_value &&
        dereferenced_value->theString->ProcessVariableReferenceCases(
            referenced_object) == HY_STRING_DIRECT_REFERENCE) {
      if (!is_global_ref && context) {
        referenced_object = *context & '.' & referenced_object;
      }
      return is_global_ref ? HY_STRING_GLOBAL_DEREFERENCE
                           : HY_STRING_LOCAL_DEREFERENCE;
    }
  }

  if (is_func_ref) {
    referenced_object = Cut(0, sLength - 2);
    if (referenced_object.IsValidIdentifier()) {
      referenced_object = (context ? (*context & '.' & referenced_object)
                                   : (referenced_object)) & '&';
      return HY_STRING_DIRECT_REFERENCE;
    }
  } else {
    if (IsValidIdentifier()) {
      referenced_object = context ? (*context & '.' & *this) : (*this);
      return HY_STRING_DIRECT_REFERENCE;
    }
  }

  referenced_object = empty;
#endif
  return HY_STRING_INVALID_REFERENCE;
}


