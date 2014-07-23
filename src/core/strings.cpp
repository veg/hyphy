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

#ifndef __HYPHYXCODE__
#include "gnuregex.h"
#else
#include "regex.h"
#include <unistd.h>
#endif

#ifdef __UNIX__
#if !defined __MINGW32__
#include <sys/utsname.h>
#endif
#include <unistd.h>
#endif

#ifdef __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

#include "hy_list_numeric.h"
#include "hy_list_reference.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>


_String compileDate = __DATE__,
        __KERNEL__VERSION__ =
            _String("2.13") & compileDate.Cut(7, 10) &
            compileDate.Cut(0, 2).Replace("Jan", "01", true)
                .Replace("Feb", "02", true).Replace("Mar", "03", true)
                .Replace("Apr", "04", true).Replace("May", "05", true)
                .Replace("Jun", "06", true).Replace("Jul", "07", true)
                .Replace("Aug", "08", true).Replace("Sep", "09", true)
                .Replace("Oct", "10", true).Replace("Nov", "11", true)
                .Replace("Dec", "12", true) &
            compileDate.Cut(4, 5).Replace(" ", "0", true) & "beta";

_String empty(""), emptyAssociativeList("{}"),
    hyphyCiteString(
        "\nPlease cite S.L. Kosakovsky Pond, S. D. W. Frost and S.V. Muse. "
        "(2005) HyPhy: hypothesis testing using phylogenies. Bioinformatics "
        "21: 676-679 if you use HyPhy in a publication\nIf you are a new HyPhy "
        "user, the tutorial located at http://www.hyphy.org/docs/HyphyDocs.pdf "
        "may be a good starting point.\n");

char defaultReturn = 0;
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
  return defaultReturn;
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
  return defaultReturn;
}

// Assignment operator
void _String::operator=(const _String& s) {
  if (s_data) {
    free(s_data);
  }
  Duplicate(&s);
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


/*!!!!!!!!!!!!!!!!!!
 
 DONE UP TO HERE
 
 !!!!!!!!!!!!!!!!!!!!*/



long _String::FindEndOfIdent(long start, long end, char wild) {
  if (s_length == 0) {
    return HY_NOT_FOUND;
  }

  if (start == -1) {
    start = ((long) s_length) - 1;
  }
  if (end == -1) {
    end = ((long) s_length) - 1;
  }

  long i = start;

  for (; i <= end; i++)
    if (!(isalnum(s_data[i]) || s_data[i] == '.' || s_data[i] == wild ||
          s_data[i] == '_')) {
      break;
    }

  if (i > start + 2 && s_data[i - 1] == '_' && s_data[i - 2] == '_') {
    return i - 3;
  }

  return i - 1;
}
 

long _String::ExtractEnclosedExpression(long &from, char open, char close,
                                        bool respectQuote, bool respectEscape) {
  long currentPosition = from, currentLevel = 0;

  bool isQuote = false, doEscape = false;

  while (currentPosition < s_length) {
    char thisChar = s_data[currentPosition];

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

  while (currentPosition < s_length) {
    char thisChar = s_data[currentPosition];
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
            for (long s = 0; s < terminators.s_length; s++)
              if (thisChar == terminators.s_data[s]) {
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






_String *_String::Sort(_SimpleList *index) {
  if (index) {
    index->Clear();
  }

  if (s_length) {
    _SimpleList charList(s_length);
    if (index) {
      for (unsigned long i = 0UL; i < s_length; i++) {
        charList << s_data[i];
        (*index) << i;
      }
      charList.Sort (true, index);
    } else {
      for (unsigned long i = 0; i < s_length; i++) {
        charList << s_data[i];
      }

      charList.Sort();
    }
    _String *sorted = new _String(s_length);
    for (unsigned long i = 0; i < s_length; i++) {
      sorted->s_data[i] = charList.AtIndex(i);
    }

    return sorted;
  }

  return new _String;
}





_Parameter _String::toNum(void) const{
  if (s_length == 0UL) {
    return 0.;
  }
  char *endP;
  return strtod(s_data, &endP);
}



_Parameter _String::ProcessTreeBranchLength(void) {
  _Parameter res = -1.;

  if (s_length) {
    if (s_data[0] == ':') {
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
  if (s_length >= 2) {
    long from = 0, to = ExtractEnclosedExpression(from, '"', '"', false, true);

    if (from == 0 && to == s_length - 1) {
      if (stripQuotes) {
        Trim(1, s_length - 2);
      }
      return true;
    }
  }
  return false;
}

//TODO: This is a global function.
bool hyIDValidator(_String *s) { return s->IsValidIdentifier(false); }

/*
==============================================================
Space Methods
==============================================================
*/

//Replace all space runs with a single space
void _String::CompressSpaces(void) {
  _StringBuffer temp(s_length + 1);
  bool skipping = false;

  for (unsigned long k = 0UL; k < s_length; k++)
    if (isspace(s_data[k])) {
      if (!skipping) {
        skipping = true;
        temp << ' ';
      }
    } else {
      temp << s_data[k];
      skipping = false;
  }
  *this = temp;
}

//Locate the first non-space charachter of the string
char _String::FirstNonSpace(long start, long end, char direction) {
  long r = FirstNonSpaceIndex(start, end, direction);
  return r == -1 ? 0 : s_data[r];
}

//Locate the first non-space charachter of the string
long _String::FirstNonSpaceIndex(long start, long end, char direction) {
  if (start == -1) {
    start = ((long) s_length) - 1;
  }
  if (end == -1) {
    end = ((long) s_length) - 1;
  }
  if (direction < 0) {
    //long t = start;
    start = end;
    end = start;
  }
  if (s_length && (start < s_length) && (!isspace(s_data[start]))) {
    return start; // first char is non-space
  }
  char *str = s_data + start;
  for (int i = start; i <= end; i += direction, str += direction)
    if (!(((*str >= 9) && (*str <= 13)) || (*str == ' '))) {
      return i;
    }

  return HY_NOT_FOUND;
}

//Locate the first non-space charachter of the string
long _String::FirstSpaceIndex(long start, long end, char direction) {
  if (start == -1) {
    start = ((long) s_length) - 1;
  }
  if (end == -1) {
    end = ((long) s_length) - 1;
  }
  if (direction < 0) {
    //long t = start;
    start = end;
    end = start;
  }
  if (s_length && (isspace(s_data[start]))) {
    return start; // first char is non-space
  }
  char *str = s_data + start;
  for (int i = start; i <= end; i += direction, str += direction)
    if ((((*str >= 9) && (*str <= 13)) || (*str == ' '))) {
      return i;
    }

  return HY_NOT_FOUND;
}

//Remove all spaces
void _String::KillSpaces(_String &result) {
  _StringBuffer temp(s_length + 1);
  for (unsigned long k = 0UL; k < s_length; k++)
    if (!isspace(s_data[k])) {
      temp << s_data[k];
    }
  result = temp;
}





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

    if (!s_length) {
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
      long f = lastPath->s_length - 2, k = 0;

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
      if ((s_data[0] == ':')) {
        Trim(1, -1);
      }
      *this = Replace(':', '\\', true);
    }
  }

  if (Find(':') == -1 && Find("\\\\", 0, 1) == -1) { // relative path

    if (pathNames.lLength) {
      _String *lastPath = (_String *)pathNames(pathNames.lLength - 1);
      long f = lastPath->s_length - 2, k = 0;
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
        if (lastPath->s_data[lastPath->s_length - 1] != '\\') {
          *this = *lastPath & '\\' & (*this);
        } else {
          *this = *lastPath & (*this);
        }
      } else {
        *this = lastPath->Cut(0, f + 1) & (*this);
      }
    }

  }

  _String escapedString(s_length, true);
  for (long stringIndex = 0; stringIndex < s_length; stringIndex++) {
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
    if (s_data[0] == '/') {
      rootPath = true;
      *this = volumeName & Cut(1, -1);
    }

    if (beginswith("..")) {
      *this = _String('/') & Cut(2, -1);
    }

    *this = Replace("/", ":", true);
    *this = Replace("..", "", true);

    if (s_data[0] != ':' && !rootPath) {
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
        if (s_data[0] != ':') {
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

  if (s_data[0] == ':') { // relative path
    long f = -1, k = 0;
    if (pathNames.lLength) {
      _String *lastPath = (_String *)pathNames(pathNames.lLength - 1);
      // check the last stored absolute path and reprocess this relative path
      // into an absolute.
      while (s_data[k] == ':') {
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
  if (relPath.s_data[0] != '/') { // relative path
    long f = -1, k = 0;
    f = s_length - 2;
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
  for (k = 0;(k < s_length) && (k < p2.s_length) && (s_data[k] == p2.s_data[k]);
       k++)
    ;
  if (k > 0) {
    while (s_data[k] != separator) {
      k--;
    }
    if (k > 0) {
      long m = k + 1, levels = 0;
      for (; m < s_length; m++)
        if (s_data[m] == separator) {
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

//TODO: These are global methods. Should they even be here?
char GetPlatformDirectoryChar(void) {
  char c = '/';
#ifdef __MAC__
  c = ':';
#endif
#if defined __WINDOZE__ || defined __MINGW32__
  c = '\\';
#endif

  return c;
}

_String GetVersionString(void) {
  _String theMessage = _String("HYPHY ") & __KERNEL__VERSION__;
#ifdef __MP__
  theMessage = theMessage & "(MP)";
#endif
#ifdef __HYPHYMPI__
  theMessage = theMessage & "(MPI)";
#endif
  theMessage = theMessage & " for ";
#ifdef __MAC__
  theMessage = theMessage & "MacOS";
#ifdef __HYPHYXCODE__
  theMessage = theMessage & "(Universal Binary)";
#else
#ifdef TARGET_API_MAC_CARBON
  theMessage = theMessage & "(Carbon)";
#endif
#endif
#endif
#ifdef __WINDOZE__
  theMessage = theMessage & "Windows (Win32)";
#endif
#ifdef __UNIX__
#if !defined __HEADLESS_WIN32__ && !defined __MINGW32__
  struct utsname name;
  uname(&name);
  theMessage = theMessage & name.sysname & " on " & name.machine;
#endif
#if defined __MINGW32__
  theMessage = theMessage & "MinGW "; // " & __MINGW32_VERSION;
#endif
#endif
  return theMessage;
}

_String GetTimeStamp(bool doGMT) {
  time_t cTime;
  time(&cTime);

  if (doGMT) {
    tm *gmt = gmtime(&cTime);
    return _String((long) 1900 + gmt->tm_year) & '/' &
           _String(1 + (long) gmt->tm_mon) & '/' &
           _String((long) gmt->tm_mday) & ' ' & _String((long) gmt->tm_hour) &
           ':' & _String((long) gmt->tm_min);
  }

  tm *localTime = localtime(&cTime);

  return asctime(localTime);

}

/*
==============================================================
Identifier Methods
==============================================================
*/

bool _String::IsValidIdentifier(bool strict) const {
#ifndef HY_2014_REWRITE_MASK
  if (s_length == 0) {
    return false;
  }

  if (strict) {
    if (!(isalpha(s_data[0]) || s_data[0] == '_')) {
      return false;
    }
  } else if (!(isalnum(s_data[0]) || s_data[0] == '_')) {
    return false;
  }

  for (unsigned long p = 1; p < s_length; p++) {
    char c = s_data[p];
    if (!(isalnum(c) || c == '_' || (strict && c == '.'))) {
      return false;
    }
  }

  // check to see if it's not a keyword / function name etc
  
  return hyReservedWords.Find(this) == HY_NOT_FOUND;
#else
  return false;
#endif

  
}

bool _String::IsValidRefIdentifier(void) const {
  if (s_length < 2) {
    return false;
  }
  if (s_data[s_length - 1] == '&') {
    return Cut(0, s_length - 2).IsValidIdentifier();
  }
  return false;
}

//Convert a string to a valid ident
void _String::ConvertToAnIdent(bool strict) {
  _StringBuffer *result = new _StringBuffer((unsigned long) s_length + 1UL);

  if (s_length) {
    if (strict) {
      if (((s_data[0] >= 'a') && (s_data[0] <= 'z')) ||
          ((s_data[0] >= 'A') && (s_data[0] <= 'Z')) || (s_data[0] == '_')) {
        (*result) << s_data[0];
      } else {
        (*result) << '_';
      }
    } else {
      if (((s_data[0] >= 'a') && (s_data[0] <= 'z')) ||
          ((s_data[0] >= 'A') && (s_data[0] <= 'Z')) || (s_data[0] == '_') ||
          ((s_data[0] >= '0') && (s_data[0] <= '9'))) {
        (*result) << s_data[0];
      } else {
        (*result) << '_';
      }
    }

    long l = 0;
    for (unsigned long k = 1UL; k < s_length; k++) {
      unsigned char c = s_data[k];
      if (_hyValidIDChars.valid_chars[c]) {
        (*result) << c;
        l++;
      } else if (result->s_data[l] != '_') {
        (*result) << '_';
        l++;
      }
    }
  }

  CopyDynamicString(result, true);
}

_String _String::ShortenVarID(_String &containerID) {
  long matched = -1,
       upTo = s_length < containerID.s_length ? s_length : containerID.s_length, k;

  for (k = 0; k < upTo; k++) {
    if (s_data[k] != containerID.s_data[k]) {
      break;
    } else if (s_data[k] == '.') {
      matched = k;
    }
  }

  if ((upTo == containerID.s_length) && (upTo < s_length) && (k == upTo) &&
      (s_data[upTo] == '.')) {
    matched = upTo;
  }

  return Cut(matched + 1, -1);
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

  errCode = regcomp(res, source->s_data,
                    REG_EXTENDED | (caseSensitive ? 0 : REG_ICASE));

  if (errCode) {
    FlushRegExp((Ptr) res);
    return nil;
  }
  return (Ptr) res;
}

void _String::RegExpMatch(Ptr pattern, _SimpleList &matchedPairs) {
  if (s_length) {
    regex_t *regEx = (regex_t *)pattern;

    regmatch_t *matches = new regmatch_t[regEx->re_nsub + 1];
    int errNo = regexec(regEx, s_data, regEx->re_nsub + 1, matches, 0);
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
  if (s_length) {
    regex_t *regEx = (regex_t *)pattern;

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
}

void _String::RegExpMatchOnce(_String *pattern, _SimpleList &matchedPairs,
                              bool caseSensitive, bool handleErrors) {
  if (s_length) {
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
  bool is_func_ref = getChar(s_length - 1) == '&';

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
    referenced_object = Cut(0, s_length - 2);
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


