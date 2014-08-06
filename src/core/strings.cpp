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
  s_length = 0UL;
  s_data = nil;
}

//Length constructor
_String::_String(const unsigned long sL) {

    s_length = sL;
    s_data = (char*)checkPointer(MemAllocate(sL + 1UL));
    memset (s_data, 0, sL + 1UL);
}

//Length constructor
_String::_String(long sL) {
  char s[32];
  snprintf(s, sizeof(s), "%ld", sL);
  for (s_length = 0; s[s_length]; s_length++)
    ;
  checkPointer(s_data = (char *)MemAllocate(s_length + 1));
  memcpy(s_data, s, s_length + 1);
}

_String::_String(const _String &source, long from, long to) {
  if (source.s_length) {
    
    long requested_range = source.NormalizeRange(from, to);
    
    if (requested_range > 0L) {
      s_length = requested_range;
      s_data = (char *)MemAllocate(s_length + 1UL);

      if (s_length > 32UL) {
        memcpy(s_data, source.s_data + from, s_length);
      } else
        for (unsigned long k = 0UL; k < s_length; k++) {
          s_data[k] = source.s_data[k + from];
        }

      s_data[s_length] = 0;
      return;
      
    }
  } 
    
  s_length = 0UL;
  s_data = (char *)MemAllocate(1UL);
  s_data[0] = 0;

}

//Stack copy contructor
_String::_String(const _String &s) { Duplicate((BaseRef) & s); }

_String::_String(_String *s) { CopyDynamicString(s, false); }

//Data constructor
_String::_String(const char *s) {
  // room for the null terminator
  s_length = strlen(s);
  s_data = (char *)MemAllocate(s_length + 1UL);
  memcpy(s_data, s, s_length + 1UL);
}


//Data constructor
_String::_String(const wchar_t *wc) {
    unsigned long allocated = wcslen (wc);
    s_length = 0;
    checkPointer(s_data = (char *)MemAllocate(allocated + 1));
    for (unsigned long cid = 0; cid < allocated; cid ++) {
        int this_char = wctob (wc[cid]);
        if (this_char != WEOF) {
            s_data[s_length++] = (char) this_char;
        }
    }
    if (s_length != allocated) {
        s_data = (char *)MemReallocate((char *)s_data, (s_length+1) * sizeof(char));
    }
    s_data[s_length] = 0;
}


//Data constructor
_String::_String(const char s) {
  s_length = 1UL;
  s_data = (char *)MemAllocate(2UL);
  s_data[0] = s;
  s_data[1] = 0;
}

//Data constructor
_String::_String(_Parameter val, const char *format) {
  char s_val[128];
  s_length = snprintf(s_val, 128, format ? format : PRINTF_FORMAT_STRING, val);
  s_data = (char *)MemAllocate(s_length + 1UL);
  for (unsigned long k = 0; k <= s_length; k++) {
    s_data[k] = s_val[k];
  }
}

//READ FROM FILE
_String::_String(FILE *F) {
  s_length = 0;
  s_data = nil;
  if (F) {
    fseek(F, 0, SEEK_END);
    s_length = (unsigned long) ftell(F);
    s_data = (char *)MemAllocate(s_length + 1);
    rewind(F);
    fread(s_data, 1, s_length, F);
    s_data[s_length] = 0;
  }
}

//Destructor
_String::~_String(void) {
  if (s_data) {
    free(s_data);
  }
}

void _String::CopyDynamicString(_String *s, bool flushMe) {
  if (flushMe && s_data) {
    free(s_data);
  }
  s_length = s->s_length;
  if (s->SingleReference ()) {
    s_data = s->s_data;
    s->s_data = nil;
    DeleteObject(s);
  } else {
    checkPointer(s_data = (char *)MemAllocate(s_length + 1UL));
    if (s->s_data) {
      memcpy(s_data, s->s_data, s_length + 1);
    } else {
      s_data[0] = 0;
    }
    s->RemoveAReference();
  }
}

void _String::Duplicate(BaseRefConst ref) {
  _String const *s = (_String const*)ref;
  s_length = s->s_length;

  if (s->s_data) {
    s_data = (char *)MemAllocate(s_length + 1UL);
    memcpy(s_data, s->s_data, s_length + 1UL);
  } else {
    s_data = nil;
  }

}

void _String::DuplicateErasing(BaseRefConst ref) {
  if (s_data) {
    free(s_data);
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
unsigned long _String::Length(void) const { return s_length; }

void _String::setChar(unsigned long index, char c) {
  if (index < s_length) {
    s_data[index] = c;
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
  if (((unsigned long) index) < s_length) {
    return s_data[index];
  }
  return _hyStringDefaultReturn;
}

//Element location functions
char _String::operator()(long index) const {
  return getChar(index);
}

  //Element location functions
const char _String::getChar(long index) const {
  if (index >= 0L && index < s_length) {
    return s_data[index];
  }
  return _hyStringDefaultReturn;
}


// Assignment operator
void _String::operator=(const _String& s) {
  if (s_data) {
    free(s_data);
  }
  Duplicate(&s);
}

_Parameter _String::toNum(void) const{
    if (s_length == 0UL) {
        return 0.;
    }
    char *endP;
    return strtod(s_data, &endP);
}



//Append operator
const _String _String::operator&(const _String& s)  const {
  unsigned long combined_length = s_length + s.s_length;

  if (combined_length == 0UL) {
    return empty;
  }

  _String res(combined_length);

  if (s_length) {
    memcpy(res.s_data, s_data, s_length);
  }

  if (s.s_length) {
    memcpy(res.s_data + s_length, s.s_data, s.s_length);
  }

  res.s_data[res.s_length] = 0;
  return res;
}


//Return good ole char*
_String::operator const char *(void) const { return s_data; }


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
  unsigned long upTo = MIN (s_length, s->s_length);
  
  for (unsigned long i = 0UL; i < upTo; i++) {
    if (s_data[i] < s->s_data[i]) {
      return -1;
    }
    if (s_data[i] > s->s_data[i]) {
      return 1;
    }
  }
  
  if (s_length == s->s_length) {
    return 0;
  }
  
  return s_length < s->s_length ? -1 : 1;
}

bool _String::Equal(const _String *s) const {
  return Compare (s) == 0;
}

bool _String::Equal(const char c) const {
  return s_length == 1 && s_data[0] == c;
}



bool _String::EqualWithWildChar(const _String &s, const char wildchar) const {
  // wildcards only matter in the second string
  
  if (s.s_length > 0UL && wildchar != '\0') {
    unsigned long   match_this_char = 0UL;
    bool            is_wildcard = s.s_data[match_this_char] == wildchar,
                    scanning_s = is_wildcard;
    
    unsigned long i = 0UL;
    while (1) {
      
      if (scanning_s) { // skip consecutive wildcards in "s"
        scanning_s = s.s_data[++match_this_char] == wildchar;
      } else {
        if (s_data[i] == s.s_data[match_this_char]) {
              // try character match
              // note that the terminal '0' characters will always match, so
              // this is where we terminate
          i++;
          match_this_char++;
          if (i > s_length || match_this_char > s.s_length) {
            break;
          }
          is_wildcard =  s.s_data[match_this_char] == wildchar;
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
    
    return match_this_char > s.s_length;
  } else {
    return s_length == 0UL;
  }
  
  return false;
}

  //Begins with string
bool _String::startswith(const _String& s, bool caseSensitive) const{
  return (caseSensitive ? Find (s, 0L, s.s_length - 1L)
          : FindAnyCase (s, 0L, s.s_length -1L)) == 0L;
}

bool _String:: endswith(const _String& s, bool caseSensitive) const{
  if (s_length >= s.s_length) {
    return (caseSensitive ? Find (s, s_length - s.s_length)
            : FindAnyCase (s, s_length - s.s_length)) == s_length - s.s_length;
  }
  return false;
}

/*
==============================================================
Manipulators
==============================================================
*/

// a convenicence function not to write (const char*)(*this)[]

const char *_String::getStr(void) const { return s_data; }


long  _String::NormalizeRange(long & from, long & to) const {

  if (s_length == 0UL) {
    return 0L;
  }

  if (from < 0L) {
    from = 0L;
  }
  
  if (to < 0L || to >= s_length ) {
    to = s_length - 1UL;
  }
  
  return to - from + 1L;
  
}

//s[0]...s[s_length-1] => s[s_length-1]...s[0]
void _String::Flip(void) {
  for (unsigned long i = 0UL; i < (s_length >> 1); i++) {
    char c;
    SWAP  (s_data[i], s_data[s_length - 1 - i], c);
  }
}

const _String _String::Chop(long from, long to) const {
  
  long resulting_length = NormalizeRange(from,to);
  
  if (resulting_length > 0L) {
    _String res((unsigned long)(s_length - resulting_length));
    if (from > 0L) {
      memcpy(res.s_data, s_data, from);
    }
    if (to + 1L < s_length) {
      memcpy(res.s_data + from, s_data + to + 1L, s_length - to - 1L);
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
    if (to < (long)s_length - 1UL) {
      memmove(s_data + from, s_data + to + 1, s_length - to - 1);
    }
    s_length -= resulting_length;
    s_data = MemReallocate(s_data, sizeof(char) * (s_length + 1UL));
    s_data[s_length] = 0;
  }
}

//Insert a single char operator
void _String::Insert(char c, long pos) {
  if (pos < 0L || pos >= s_length) {
    pos = s_length;
  }

  s_data = MemReallocate(s_data, sizeof(char) * (s_length + 2));

  if (pos < s_length) {
    memmove(s_data + pos + 1, s_data + pos, s_length - pos);
  }

  s_data[pos] = c;
  s_data[++s_length] = 0;
}

//Cut string from, to (-1 for any means from beginning/to end)
void _String::Trim(long from, long to) {
  
  long resulting_length = NormalizeRange(from, to);
  
  if (resulting_length > 0L) {
      if (from > 0L) {
        memmove(s_data, s_data + from, resulting_length);
      }
      s_length = resulting_length;
      s_data = MemReallocate(s_data, resulting_length + 1UL);
      s_data[resulting_length] = 0;
  } else {
      s_length = 0UL;
      s_data = MemReallocate(s_data, 1UL);
      s_data[0] = 0;
  }
}


void _String::StripQuotes(char open_char, char close_char) {
  if (s_length >= 2UL) {
    if (getChar(0) == open_char && getChar (s_length - 1UL) == close_char) {
      Trim (1, s_length - 2UL);
    }
  }
}

/*
==============================================================
Case Methods
==============================================================
*/

const _String _String::UpCase(void) const {
  _String res (s_length);
  for (unsigned long i = 0UL; i < s_length; i++) {
    res.s_data[i] = toupper(s_data[i]);
  }
  return res;
}

const _String _String::LoCase(void) const {
  _String res (s_length);
  for (unsigned long i = 0UL; i < s_length; i++) {
    res.s_data[i] = tolower(s_data[i]);
  }
  return res;
}

/*
==============================================================
Search and replace
==============================================================
*/

long _String::Find(const _String& s, long from, long to) const {

  if (s.s_length) {
    long span = NormalizeRange(from, to);
    if (span >= (long) s.s_length) {
      const unsigned long upper_bound = to - s.s_length + 2L;
      for (unsigned long i = from; i < upper_bound ; i++) {
        unsigned long j = 0UL;
        for (;j < s.s_length; j++) {
          if (s_data[i + j] != s.s_data[j]) {
            break;
          } 
        }
        if (j == s.s_length) {
          return i;
        }
      }
    }
  }
  return HY_NOT_FOUND;
}

  //Find first occurence of the string between from and to
long _String::Find(const char s, long from, long to) const {
  if (s_length) {
    long span = NormalizeRange(from, to);
    if (span > 0L) {
      char sentinel = s_data[to+1];
      s_data[to+1] = s;
      long index = from;
      while (s_data[index] != s) {
        index++;
      }
      s_data[to+1] = sentinel;
      return index <= to ? index : HY_NOT_FOUND;
    }
  }
  
  return HY_NOT_FOUND;
}


//Find first occurence of the string between from and to
long _String::FindBackwards(const _String & s, long from, long to) const {

  if (s.s_length) {
    long span = NormalizeRange(from, to);
    if (span >= (long) s.s_length) {
      const long upper_bound = to - s.s_length + 1L;
      for (long i = upper_bound; i >= from; i--) {
        unsigned long j = 0UL;
        for (;j < s.s_length; j++) {
          if (s_data[i + j] != s.s_data[j]) {
            break;
          } 
        }
        if (j == s.s_length) {
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

  if (s_length < s.s_length || s.s_length == 0UL) {
    return *this;
  }

  _StringBuffer replacementBuffer;
  unsigned long anchor_index = 0UL;
  for (; anchor_index <= s_length - s.s_length; anchor_index ++) {
    unsigned long search_index = 0UL;
    for (; search_index < s.s_length; search_index++) {
      if (s_data[anchor_index + search_index] != s.s_data[search_index]) {
        break;
      }
    }

    if (search_index == s.s_length) {
      replacementBuffer << d;
      anchor_index += s.s_length - 1UL;
      if (replace_all == false) {
        anchor_index ++;
        break;
      }
    } else {
      replacementBuffer << s_data[anchor_index];
    }
  }
  
  replacementBuffer.appendSubstring(*this, anchor_index, HY_NOT_FOUND);
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
    
    cp = cpp + splitter.s_length;
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
  
  unsigned int alphabet_length = alphabet ? alphabet->s_length : 127;
  
  if (length > 0UL && alphabet_length > 0UL) {
    for (unsigned long c = 0UL; c < length; c++) {
      unsigned long idx = genrand_int32() % alphabet_length;
      if (alphabet) {
        random << alphabet->s_data[idx];
      } else {
        random << (char)(1UL + idx);
      }
    }
  }
  
  return random;
}

long _String::Adler32(void) const {
  
  unsigned long len = s_length,
                a = 1UL,
                b = 0UL,
                i = 0UL;
  
  while (len) {
    unsigned long tlen = len > 5550UL ? 5550UL : len;
    len -= tlen;
    do {
      a += s_data[i++];
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
    
    if (s_length > 0UL) {
        _hyListOrderable<char> sorted;
        
        if (index) {
            for (unsigned long i = 0UL; i < s_length; i++) {
                sorted << s_data[i];
                (*index) << i;
            }
            sorted.Sort (true, index);
        } else {
            for (unsigned long i = 0; i < s_length; i++) {
                sorted << s_data[i];
            }
            sorted.Sort();
        }
        _String result (s_length);
        
        for (unsigned long i = 0UL; i < s_length; i++) {
            result.setChar (i, sorted.AtIndex(i));
        }
        
        return result;
    }
    
    return empty;
}

_Parameter _String::ProcessTreeBranchLength(_Parameter min_value) const {
  _Parameter res = -1.;
  
  if (s_length) {
    if (s_data[0] == ':') {
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
    
    return s_length > 0UL && _IsValidIdentifierAux (allow_compounds) == Length() - 1UL;

#ifndef HY_2014_REWRITE_MASK
    // TO DO -- MOVE THIS CHECK ELSEWHERE?
    return hyReservedWords.Find(this) == HY_NOT_FOUND;
#endif
    
    
}

long _String::_IsValidIdentifierAux(bool allow_compounds, char wildcard) const {
    
    unsigned long current_index = 0UL;
    
    bool          first     = true;
    
    for (; current_index < s_length; current_index ++) {
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
        
        if (getChar(prefix_length) == '.' && s_length > prefix_length + 1L) {
            return Cut(prefix_length + 1L, HY_NOT_FOUND);
        }
    }
    return *this;
}

  //Convert a string to a valid ident
const _String  _String::ConvertToAnIdent(bool strict) const {
  _StringBuffer converted;
  const char default_placeholder = '_';
  
  if (s_length) {
    unsigned long index = 0UL;
    if (strict) {
      converted << (_hyValidIDChars.canBeFirst (s_data[0]) ? s_data[0] : default_placeholder);
      index ++;
    }
    for (;index < s_length; index++) {
      if (_hyValidIDChars.isValidChar(s_data[index])) {
        converted << s_data[index];
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
    _StringBuffer temp(s_length + 1UL);
    bool skipping = false;
    
    for (unsigned long k = 0UL; k < s_length; k++) {
        if (!isspace(s_data[k])) {
            temp << s_data[k];
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
    _StringBuffer temp(s_length + 1UL);
    for (unsigned long k = 0UL; k < s_length; k++) {
        if (!isspace(s_data[k])) {
            temp << s_data[k];
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
        if (!isspace(s_data[from])) {
          return from;
        }
      }
    } else {
      for (; to>=from; to--) {
        if (!isspace(s_data[to])) {
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
        if (isspace(s_data[from])) {
          return from;
        }
      }
    } else {
      for (; to>=from; to--) {
        if (isspace(s_data[to])) {
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
  return r == HY_NOT_FOUND ? _hyStringDefaultReturn : s_data[r];
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
  
  if (s_length) {
    regmatch_t *matches = new regmatch_t[regEx->re_nsub + 1];
    int errNo = regexec(regEx, s_data, regEx->re_nsub + 1, matches, 0);
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
  
  if (s_length) {
    
    regmatch_t *matches = new regmatch_t[regEx->re_nsub + 1];
    int errNo = regexec(regEx, s_data, regEx->re_nsub + 1, matches, 0);
    while (errNo == 0) {
      long offset = matchedPairs.countitems()
      ? matchedPairs.Element(-1) + 1
      : 0;
      
      matchedPairs << matches[0].rm_so + offset;
      matchedPairs << matches[0].rm_eo - 1 + offset;
      
      offset += matches[0].rm_eo;
      if (offset < s_length) {
        errNo = regexec(regEx, s_data + offset, regEx->re_nsub + 1, matches, 0);
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
  if (s_length) {
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

  if (s_length == 0) {
    return 0;
  }

  if (rec) {
    (*rec) << 0;
  }

  long cp = 1, pH = 1;

  while (cp < s_length) {
    long maxExtension = 0;

    for (long ip = 0; ip < cp; ip++) {
      long sp = ip, mp = cp;

      while ((mp < s_length) && (s_data[mp] == s_data[sp])) {
        mp++;
        sp++;
      }

      if (mp == s_length) {
        maxExtension = s_length - cp;
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









