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
#include "errorfns.h"






#ifdef __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <unistd.h>



char _hyStringDefaultReturn = 0;
unsigned long _String::storageIncrement = 32;

/*extern int _hy_mpi_node_rank;
long  loopCount      = 3;
char* addrBreak      = 0x583088;*/

struct _hyValidIDCharsType {
  unsigned char valid_chars[256];
  
  bool canBeFirst (unsigned char c) {
    return valid_chars[c] == 2;
  }

  bool isValidChar (unsigned char c) {
    return valid_chars[c] > 0;
  }
  
  _hyValidIDCharsType(void) {
    for (int c = 0; c < 256; c++) {
      valid_chars[c] = 0;
    }
    for (unsigned char c = 'a'; c <= 'z'; c++) {
      valid_chars[c] = 2;
    }
    for (unsigned char c = 'A'; c <= 'Z'; c++) {
      valid_chars[c] = 2;
    }
    for (unsigned char c = '0'; c <= '9'; c++) {
      valid_chars[c] = 1;
    }
    valid_chars[(unsigned char) '_'] = 2;
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

_Parameter _String::ProcessTreeBranchLength(_Parameter min_value) const {
  _Parameter res = -1.;
  
  if (sLength) {
    if (sData[0] == ':') {
      res = Cut(1, -1).toNum();
    } else {
      res = toNum();
    }
    
    return res < min_value ? min_value : res;
  }
  
  return res;
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

  //Convert a string to a valid ident
const _String  _String::ConvertToAnIdent(bool strict) const {
  _StringBuffer converted;
  const char default_placeholder = '_';
  
  if (sLength) {
    unsigned long index = 0UL;
    if (strict) {
      converted << (_hyValidIDChars.canBeFirst (sData[0]) ? sData[0] : default_placeholder);
      index ++;
    }
    for (;index < sLength; index++) {
      if (_hyValidIDChars.isValidChar(sData[index])) {
        converted << sData[index];
      } else {
        if (index && converted.getChar(converted.Length()-1UL) != default_placeholder)  {
          converted << default_placeholder;
        }
      }
    }
  } else {
    converted << default_placeholder;
  }
  
  return converted;
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


/*
 ==============================================================
 Regular Expression Methods
 ==============================================================
 */

const _String GetRegExpError(int error) {
  char buffer[512];
  buffer[regerror(error, nil, buffer, 511)] = 0;
  return _String("Regular Expression error:") & buffer;
}

void FlushRegExp(regex_t* regExpP) {
  regfree(regExpP);
  delete regExpP;
}

regex_t* PrepRegExp(const _String& source, int &errCode, bool caseSensitive) {
  regex_t *res = new regex_t;
  checkPointer(res);
  
  errCode = regcomp(res, source.getStr(),
                    REG_EXTENDED | (caseSensitive ? 0 : REG_ICASE));
  
  if (errCode) {
    FlushRegExp(res);
    return nil;
  }
  return res;
}

const _SimpleList _String::RegExpMatch(regex_t const* regEx) const {
  _SimpleList matchedPairs;
  
  if (sLength) {
    regmatch_t *matches = new regmatch_t[regEx->re_nsub + 1];
    int errNo = regexec(regEx, sData, regEx->re_nsub + 1, matches, 0);
    if (errNo == 0) {
      for (long k = 0L; k <= regEx->re_nsub; k++) {
        matchedPairs << matches[k].rm_so;
        matchedPairs << matches[k].rm_eo - 1;
      }
    }
    delete[] matches;
  }
  
  return matchedPairs;
}

const _SimpleList _String::RegExpMatchAll(regex_t const* regEx) const {
  _SimpleList matchedPairs;
  
  if (sLength) {
    
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
  return matchedPairs;
}

const _SimpleList _String::RegExpMatchOnce(const _String & pattern,
                              bool caseSensitive, bool handleErrors) const {
  if (sLength) {
    int errNo = 0;
    regex_t* regex = PrepRegExp(pattern, errNo, caseSensitive);
    if (regex) {
      _SimpleList hits = RegExpMatch(regex);
      FlushRegExp(regex);
      return hits;
    } else if (handleErrors) {
      warnError(GetRegExpError(errNo));
    }
  }
  return _SimpleList();
}


/*!!!!!!!!!!!!!!!!!!
 
 DONE UP TO HERE
 
 !!!!!!!!!!!!!!!!!!!!*/



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









