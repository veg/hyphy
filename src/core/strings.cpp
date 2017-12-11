/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon42@uwo.ca)
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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>

#include "global_things.h"
#include "hy_strings.h"
#include "batchlan.h"
#include "mersenne_twister.h"
#include "function_templates.h"
#include "hy_string_buffer.h"


_String   compileDate = __DATE__,
          __HYPHY__VERSION__ = _String ("2.3.4") & compileDate.Cut (7,10) & compileDate.Cut (0,2).Replace("Jan", "01", true).
                                                                                                  Replace("Feb", "02", true).
                                                                                                  Replace("Mar", "03", true).
                                                                                                  Replace("Apr", "04", true).
                                                                                                  Replace("May", "05", true).
                                                                                                  Replace("Jun", "06", true).
                                                                                                  Replace("Jul", "07", true).
                                                                                                  Replace("Aug", "08", true).
                                                                                                  Replace("Sep", "09", true).
                                                                                                  Replace("Oct", "10", true).
                                                                                                  Replace("Nov", "11", true).
                                                                                                  Replace("Dec", "12", true)
                                                                                                  & compileDate.Cut (4,5).Replace (" ", "0", true) & "beta";

using namespace hy_global;


 

struct _hy_Valid_ID_Chars_Type {
  unsigned char valid_chars[256];
  
  inline bool is_valid_first (unsigned char c) const {
    return valid_chars[c] == 2;
  }
  
  inline bool is_valid (unsigned char c) const {
    return valid_chars[c] > 0;
  }
  
  _hy_Valid_ID_Chars_Type (void) {
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
} hy_Valid_ID_Chars;


/*
==============================================================
Constructors/Destructors/Copiers
==============================================================
*/

_String::_String (void) {
  Initialize();
}

//=============================================================

void _String::Initialize (bool) {
    s_length = 0UL;
    s_data = nil;
}

//=============================================================

void _String::Clear (void) {
  s_length = 0UL;
  if (s_data) {
    free (s_data);
    s_data = nil;
  }
}

//=============================================================

_String::_String(long const number) {
    char s[64];
    s_length = snprintf(s, sizeof(s), "%ld", number);
    AllocateAndCopyString (s, s_length);
}

//=============================================================

_String::_String(const unsigned long sL) {
    s_length = sL;
    s_data = (char*) MemAllocate(sL + 1L, true);
}

//=============================================================

_String::_String(const hyFloat val, const char *format) {
    char s_val[128];
    s_length = snprintf(s_val, 128, format ? format : PRINTF_FORMAT_STRING, val);
    AllocateAndCopyString (s_val, s_length);
}

//=============================================================

_String::_String(const _String &s) {
    Initialize ();
    Duplicate(& s);
}

//=============================================================

_String::_String(_String *s) {
    if (s->CanFreeMe ()) {
        s_data       = s->s_data;
        s_length     = s->s_length;
        s->s_data    = nil;
        DeleteObject (s);
    } else {
        AllocateAndCopyString (s_data, s_length);
        s->RemoveAReference();
    }
}


//=============================================================
_String::_String(const _String &source, long start, long end) {
    if (source.s_length) {
        
        long requested_range = source.NormalizeRange(start, end);
        
        if (requested_range > 0L) {
            AllocateAndCopyString (source.s_data + start, requested_range);
            return;
            
        }
    }
    
    s_length = 0UL;
    s_data = (char *)MemAllocate(1UL);
    s_data[0] = '\0';
    
}

//=============================================================

_String::_String(const char *c_string) {
    AllocateAndCopyString (c_string, strlen(c_string));
}

//=============================================================
_String::_String(const wchar_t *wc_string) {
    unsigned long allocated = wcslen (wc_string);
    s_length = 0UL;
    s_data = (char *)MemAllocate(allocated + 1UL);
    for (unsigned long cid = 0UL; cid < allocated; cid ++) {
        int this_char = wctob (wc_string[cid]);
        if (this_char != WEOF) {
            s_data[s_length++] = (char) this_char;
        }
    }
    if (s_length != allocated) {
        s_data = (char *)MemReallocate((char *)s_data, (s_length+1) * sizeof(char));
    }
    s_data[s_length] = '\0';
}

//=============================================================
_String::_String(const char c) {
    s_length = 1UL;
    s_data = (char *)MemAllocate(2UL);
    s_data[0] = c;
    s_data[1] = '\0';
}

//=============================================================
_String::_String (const _String& str, unsigned long copies) {
    s_length = copies * str.s_length;
    s_data = (char*)MemAllocate (s_length+1UL);
    if (s_length > 0UL) {
        for (unsigned long i = 0UL; i < copies; i++) {
            memcpy (s_data + i * str.s_length, s_data, str.s_length);
        }
    }
    s_data[s_length]='\0';
}

//=============================================================
_String::_String(FILE * file, long read_this_many) {
    Initialize ();
    if (file) {
        if (read_this_many < 0) {
          fseek(file, 0, SEEK_END);
          s_length = (unsigned long) ftell(file);
          rewind(file);
        } else {
          s_length = read_this_many;
        }
        s_data = (char *)MemAllocate(s_length + 1UL);
        unsigned long read_items = fread(s_data, 1, s_length, file);
        if (read_items < s_length) {
          s_data = (char*)MemReallocate(s_data,read_items+1);
          s_length = read_items;
        }
        s_data[s_length] = '\0';
    }
}

//=============================================================
_String::~_String(void)
{
    if (CanFreeMe()) {
        if (s_data) {
            free (s_data);
            s_data = nil;
        }
        s_length = 0UL;
    } else {
        RemoveAReference();
    }
}

//=============================================================
BaseRef _String::makeDynamic (void) const {
    _String * r = new _String;
    r->Duplicate(this);
    return r;
}

//=============================================================
void    _String::Duplicate (BaseRefConst ref) {
    if (s_data) {
        free (s_data);
    }
    
    _String const * s = (_String const*)ref;
    
    s_length = s->s_length;
    s_data   = s->s_data;
    
    if (s_data) {
        AllocateAndCopyString (s->s_data, s_length);
    }
}

//=============================================================
void _String::operator = (_String const& s) {
    Duplicate (&s);
}



/*
 ==============================================================
 Private helpers
 ==============================================================
 */

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

//=============================================================

void _String::AllocateAndCopyString (const char * source_string, unsigned long length) {
    s_length = length;
    s_data = (char*) MemAllocate (length+1UL);
    if (s_length) {
        memcpy (s_data, source_string, length);
    }
    s_data [length] = '\0';
}


/*
==============================================================
Getters and setters
==============================================================
*/

char& _String::operator [] (long index) {
    if (index < s_length && index >= 0L) {
        return s_data[index];
    }
    HandleApplicationError (_String ("Internal error at ") & __PRETTY_FUNCTION__ & ": an invalid index requested");
}

//=============================================================


char _String::operator () (long index) const {
  if (index >= 0L && index < s_length) {
    return s_data[index];
  }
  if (index < 0L && -index <= s_length) {
    return s_data[s_length + index];
  }
  return default_return;
}

//=============================================================

void _String::set_char (unsigned long index, char const data) {
    if (index < s_length) {
        s_data[index] = data;
    }
}

//=============================================================

const char *_String::get_str(void) const { return s_data; }

/*
 ==============================================================
 Type conversions
 ==============================================================
 */


_String::operator const char *(void) const { return s_data; }

//=============================================================

hyFloat _String::to_float (void) const{
  if (s_length == 0UL) {
    return 0.;
  }
  char *endP;
  return strtod(s_data, &endP);
}

//=============================================================

long _String::to_long (void) const {
  if (s_length == 0UL) {
    return 0L;
  }
  char * endP;
  return strtol(s_data,&endP,10);
}

//=============================================================

BaseRef _String::toStr (unsigned long) {
  AddAReference();
  return this;
}

//=============================================================

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
 Comparisons
 ==============================================================
 */

hyComparisonType _String::Compare(_String const& rhs) const {
    unsigned long up_to = MIN (s_length, rhs.s_length);
    
    for (unsigned long i = 0UL; i < up_to; i++) {
        if (s_data[i] < rhs.s_data[i]) {
            return kCompareLess;
        }
        if (s_data[i] > rhs.s_data[i]) {
            return kCompareGreater;
        }
    }

    if (s_length == rhs.s_length) {
        return kCompareEqual;
    }
    
    return s_length < rhs.s_length ? kCompareLess : kCompareGreater;
}

//=============================================================

hyComparisonType _String::CompareIgnoringCase(_String const& rhs) const {
    unsigned long up_to = MIN (s_length, rhs.s_length);
    
    for (unsigned long i = 0UL; i < up_to; i++) {
       
        char llhs = tolower (s_data[i]), lrhs = tolower (rhs.s_data[i]);
        
        
        if (llhs < lrhs) {
            return kCompareLess;
        }
        if (llhs > lrhs) {
            return kCompareGreater;
        }
    }
    
    if (s_length == rhs.s_length) {
        return kCompareEqual;
    }
    
    return s_length < rhs.s_length ? kCompareLess : kCompareGreater;
}

//=============================================================

bool _String::operator==(const _String& s) const { return Compare  (s) == kCompareEqual; }
bool _String::operator>(const _String & s) const { return Compare  (s) == kCompareGreater; }
bool _String::operator<=(const _String & s) const { return Compare (s) != kCompareGreater; }
bool _String::operator>=(const _String & s) const { return Compare (s) != kCompareLess; }
bool _String::operator!=(const _String & s) const { return Compare (s) != kCompareEqual; }
bool _String::operator<(const _String & s) const { return Compare  (s) == kCompareLess; }

bool _String::Equal(const _String& s) const { return Compare  (s) == kCompareEqual; };
bool _String::EqualIgnoringCase(const _String& s) const { return CompareIgnoringCase  (s) == kCompareEqual; };
bool _String::Equal(const char c) const {
    return s_length == 1UL && s_data[0] == c;
}

//=============================================================


bool _String::EqualWithWildChar(const _String& pattern, const char wildchar, unsigned long start_this, unsigned long start_pattern) const {
    // wildcards only matter in the second string
    
    if (pattern.s_length > start_pattern && wildchar != '\0') {
        unsigned long   match_this_char = start_pattern;
        // the position we are currently trying to match in the pattern
        
        bool            is_wildcard = pattern.s_data[match_this_char] == wildchar,
        scanning_pattern = is_wildcard;
        
        unsigned long i = start_this;
        // the position we are currently trying to match in *this
        
        while (i <= s_length) {
            
            if (scanning_pattern) { // skip consecutive wildcards in "pattern"
                scanning_pattern = pattern.s_data[++match_this_char] == wildchar;
            } else {
                if (s_data[i] == pattern.s_data[match_this_char]) {
                    
                    if (is_wildcard) {
                        // could either match the next char or consume it into the wildcard
                        if (EqualWithWildChar (pattern, wildchar, i, match_this_char)) {
                            // matching worked
                            return true;
                        } else {
                            i++;
                            continue;
                        }
                    } else {
                        // try character match
                        // note that the terminal '0' characters will always match, so
                        // this is where we terminate
                        i++;
                        match_this_char++;
                        if (i > s_length || match_this_char > pattern.s_length) {
                            break;
                        }
                        is_wildcard =  pattern.s_data[match_this_char] == wildchar;
                        scanning_pattern = is_wildcard;
                    }
                } else { // match wildcard
                    if (!is_wildcard) {
                        return false;
                    }
                    scanning_pattern = false;
                    i++;
                }
            }
        }
        
        return match_this_char > pattern.s_length;
    } else {
        return s_length == start_this;
    }
    
    return false;
}


/*
 ==============================================================
 Content-modification and extraction methods
 ==============================================================
 */


//=============================================================


//Append operator
const _String _String::operator & (const _String& rhs) const {
    unsigned long combined_length = s_length + rhs.s_length;
    
    if (combined_length == 0UL) {
        return kEmptyString;
    }
    
    _String res(combined_length);
    
    if (s_length) {
        memcpy(res.s_data, s_data, s_length);
    }
    
    if (rhs.s_length) {
        memcpy(res.s_data + s_length, rhs.s_data, rhs.s_length);
    }
    
    res.s_data[res.s_length] = '\0';
    return res;
}

//=============================================================

const _String _String::Chop(long start, long end) const{
    
    long resulting_length = NormalizeRange(start,end);
    
    if (resulting_length > 0L) {
        _String res((unsigned long)(s_length - resulting_length));
        if (start > 0L) {
            memcpy(res.s_data, s_data, start);
        }
        if (end + 1L < s_length) {
            memcpy(res.s_data + start, s_data + end + 1L, s_length - end - 1L);
        }
        
        return res;
    }
    
    return *this;
    
}

//=============================================================

const _String _String::Cut(long start, long end) const {
    return _String (*this, start, end);
}

//=============================================================

void _String::Delete(long start, long end) {
    long resulting_length = NormalizeRange(start,end);
    
    if (resulting_length > 0L) {
        if (end < (long)s_length - 1UL) {
            memmove(s_data + start, s_data + end + 1L, s_length - end - 1L);
        }
        s_length -= resulting_length;
        s_data = (char*)MemReallocate(s_data, sizeof(char) * (s_length + 1UL));
        s_data[s_length] = '\0';
    }
}

//=============================================================

void _String::Flip(void) {
    for (unsigned long i = 0UL; i < (s_length >> 1); i++) {
        char c;
        SWAP  (s_data[i], s_data[s_length - 1 - i], c);
    }
}

//=============================================================

const _String _String::Reverse(void) const {
    
    _String result (*this);
    for (unsigned long s = 0UL, e = s_length - 1L;  s < s_length; s++, e--) {
        result.s_data[s] = s_data[e];
    }
    return result;
}

//=============================================================

void _String::Insert(char c, long where) {
    if (where < 0L || where >= s_length) {
        where = s_length;
    }
    
    s_data = (char*)MemReallocate(s_data, sizeof(char) * (s_length + 2UL));
    
    if (where < s_length) {
        memmove(s_data + where + 1UL, s_data + where, s_length - where);
    }
    
    s_data[where] = c;
    s_data[++s_length] = '\0';
}

//=============================================================

void _String::Trim(long start, long end) {
    
    long resulting_length = NormalizeRange(start, end);
    
    if (resulting_length > 0L) {
        if (start > 0L) {
            memmove(s_data, s_data + start, resulting_length);
        }
        if (s_length != resulting_length) {
            s_length = resulting_length;
            s_data = (char*)MemReallocate(s_data, resulting_length + 1UL);
            s_data[resulting_length] = '\0';
        }
    } else {
        s_length = 0UL;
        s_data = (char*)MemReallocate(s_data, 1UL);
        s_data[0] = '\0';
    }
}

//=============================================================

const _String    _String::ChangeCase (hy_string_case conversion_type) const {
  _String result (s_length);
  
  auto conversion_function = conversion_type == kStringUpperCase ? toupper : tolower;
  
  for (unsigned long i = 0UL; i<s_length; i++) {
    result.s_data [i] = conversion_function (s_data[i]);
  }
}

//=============================================================

const _List _String::Tokenize(const _String& splitter) const {
  _List tokenized;
  
  long cp = 0L, cpp;
  while ((cpp = Find(splitter, cp)) != kNotFound) {
    if (cpp > cp) {
      tokenized < new _String(*this, cp, cpp - 1L);
    } else {
      tokenized < new _String;
    }
    
    cp = cpp + splitter.s_length;
  }
  
  tokenized < new _String(*this, cp);
  return tokenized;
}

  //=============================================================

const _List _String::Tokenize(bool const splitter[256]) const {
  _List tokenized;
  
  long cp = 0L, cpp;
  while ((cpp = Find(splitter, cp)) != kNotFound) {
    if (cpp > cp) {
      tokenized < new _String(*this, cp, cpp - 1L);
    } else {
      tokenized < new _String;
    }
    
    cp = cpp + 1;
  }
  
  tokenized < new _String(*this, cp);
  return tokenized;
}

//=============================================================

const _String _String::Enquote (char quote_char) const {
  return _StringBuffer (2UL + s_length) << quote_char  << *this << quote_char;
}

//=============================================================

const _String _String::Enquote (char open_char, char close_char) const {
    return _StringBuffer (2UL + s_length) << open_char  << *this << close_char;
}

//=============================================================

const _String _String::KillSpaces(void) const {
  _StringBuffer temp(s_length + 1UL);
  for (unsigned long k = 0UL; k < s_length; k++) {
    if (!isspace(s_data[k])) {
      temp << s_data[k];
    }
  }
  return temp;
}

//=============================================================

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


/*
 ==============================================================
 Search Functions
 ==============================================================
*/

long _String::Find(const _String& pattern, long start, long end) const {
  
  if (pattern.s_length) {
    long span = NormalizeRange(start, end);
    if (span >= (long) pattern.s_length) {
      const unsigned long upper_bound = end - pattern.s_length + 2L;
      for (unsigned long i = start; i < upper_bound ; i++) {
        unsigned long j = 0UL;
        for (;j < pattern.s_length; j++) {
          if (s_data[i + j] != pattern.s_data[j]) {
            break;
          }
        }
        if (j == pattern.s_length) {
          return i;
        }
      }
    }
  }
  return kNotFound;
}

//=============================================================

long _String::FindBackwards(const _String& pattern, long start, long end) const {
  
  if (pattern.s_length) {
    long span = NormalizeRange(start, end);
    if (span >= (long) pattern.s_length) {
      const long upper_bound = end - pattern.s_length + 1L;
      for (long i = upper_bound; i >= start; i--) {
        unsigned long j = 0UL;
        for (;j < pattern.s_length; j++) {
          if (s_data[i + j] != pattern.s_data[j]) {
            break;
          }
        }
        if (j == pattern.s_length) {
          return i;
        }
      }
    }
  }
  return kNotFound;
}

//=============================================================

long _String::Find(const char p, long start, long end) const {
  if (s_length) {
    long span = NormalizeRange(start, end);
    if (span > 0L) {
      char sentinel = s_data[end+1];
      s_data[end+1L] = p;
      long index = start;
      while (s_data[index] != p) {
        index++;
      }
      s_data[end+1L] = sentinel;
      return index <= end ? index : kNotFound;
    }
  }
  
  return kNotFound;
}

  //=============================================================

long _String::Find(const bool lookup[256], long start, long end) const {
  if (s_length) {
    long span = NormalizeRange(start, end);
    if (span > 0L) {
      
      for (unsigned long index = start; index <= end; index ++ ) {
        if (lookup [s_data[index]]) {
          return index;
        }
      }
    }
  }
  
  return kNotFound;
}


//=============================================================
long _String::FindAnyCase (const _String& pattern, long start, long end) const {
  
  if (pattern.s_length) {
    long span = NormalizeRange(start, end);
    if (span >= (long) pattern.s_length) {
      const unsigned long upper_bound = end - pattern.s_length + 2L;
      for (unsigned long i = start; i < upper_bound ; i++) {
        unsigned long j = 0UL;
        for (;j < pattern.s_length; j++) {
          if (toupper (s_data[i + j]) != toupper (pattern.s_data[j])) {
            break;
          }
        }
        if (j == pattern.s_length) {
          return i;
        }
      }
    }
  }
  return kNotFound;
}
//=============================================================

const _String _String::Replace(const _String& pattern, const _String replace, bool replace_all) const {
  
  if (s_length < pattern.s_length || pattern.s_length == 0UL) {
    return *this;
  }
  
  _StringBuffer replacement_buffer;
  unsigned long anchor_index = 0UL;
  for (; anchor_index <= s_length - pattern.s_length; anchor_index ++) {
    unsigned long search_index = 0UL;
    for (; search_index < pattern.s_length; search_index++) {
      if (s_data[anchor_index + search_index] != pattern.s_data[search_index]) {
        break;
      }
    }
    
    if (search_index == pattern.s_length) {
      replacement_buffer << replace;
      anchor_index += pattern.s_length - 1UL;
      if (replace_all == false) {
        anchor_index ++;
        break;
      }
    } else {
      replacement_buffer << s_data[anchor_index];
    }
  }
  
  return replacement_buffer.AppendSubstring(*this, anchor_index, kNotFound);
}
//=============================================================

long _String::FirstNonSpaceIndex(long start, long end, hy_string_search_direction direction) const {
  return _FindFirstIndexCondition(start, end, direction, [] (char c) -> bool {return !isspace (c);});
}

char _String::FirstNonSpace(long start, long end, hy_string_search_direction direction) const {
  long r = FirstNonSpaceIndex(start, end, direction);
  return r == kNotFound ? default_return : s_data[r];
}

//=============================================================

long _String::FirstSpaceIndex(long start, long end, hy_string_search_direction direction) const {
  return _FindFirstIndexCondition(start, end, direction, isspace);
}

//=============================================================

long _String::FirstNonSpaceFollowingSpace(long start, long end, hy_string_search_direction direction) const {
  long first_space = FirstSpaceIndex (start, end, direction);
  if (first_space != kNotFound) {
    if (direction == kStringDirectionForward) {
      first_space = FirstNonSpaceIndex (first_space, end, direction);
    } else {
      first_space = FirstNonSpaceIndex (start, first_space, direction);
    }
  }
  return first_space;
}

//=============================================================
bool _String::BeginsWith (_String const& pattern, bool case_sensitive, unsigned long from) const{
  if (pattern.s_length + from < s_length) {
    return (case_sensitive ? Find (pattern, from)
            : FindAnyCase (pattern, from)) == from;
  }
  
 
  return false;
}



//Ends with string
bool _String::EndsWith (_String const& pattern, bool case_sensitive) const{
  if (s_length >= pattern.s_length) {
    unsigned long length_difference = s_length - pattern.s_length;
    return (case_sensitive ? Find (pattern, length_difference)
            : FindAnyCase (pattern, length_difference)) == length_difference;
  }
  return false;
}



//Begins with string
bool _String::BeginsWithAndIsNotAnIdent (_String const& pattern) const {
  
  if (BeginsWith (pattern)) {
    if (s_length > pattern.s_length) {
      char next_char = get_char (pattern.s_length);
      if (isalnum(next_char) || next_char == '.' || next_char == '_' || next_char == '&') {
        // TODO SLKP 20170616: what is the use case for next_char == '&'?
        return false;
      }
    }
    return true;
  }
  return false;
}

/*
 ==============================================================
 Parser-related functions
 TODO: possible deprecate when the move to the grammar is effected
 ==============================================================
*/


//=============================================================


long _String::ExtractEnclosedExpression (long& from, char open, char close, int options) const {
  long   current_position = from,
         current_level    = 0L;
  
  bool       respect_quote = options & fExtractRespectQuote,
             respect_escape = options & fExtractRespectEscape,
             do_escape = false;
  
  char       quote_state = '\0';
  
  while (current_position < s_length) {
    char this_char = s_data[current_position];
    
    if (do_escape) {
      do_escape = false;
    } else {
      if ((this_char == '"' || this_char == '\'') && respect_quote && !do_escape) {
        if (quote_state == '\0') {
          quote_state = this_char;
        } else {
          if (this_char == quote_state) {
            quote_state = '\0';
          }
        }
      } else if (this_char == open && quote_state == '\0') {
        // handle the case when close and open are the same
        if (current_level == 1L && open == close && from < current_position) {
          return current_position;
        }
        current_level++;
        if (current_level == 1L) {
          from = current_position;
        }
      } else if (this_char == close && quote_state == '\0') {
        current_level--;
        if (current_level == 0L && from < current_position) {
          return current_position;
        }
        if (current_level < 0L) {
          return kNotFound;
        }
      } else if (this_char == '\\' && respect_escape && quote_state != '\0' && !do_escape) {
        do_escape = true;
      }
    }
    
    current_position++;
  }
  
  return kNotFound;
}

//=============================================================


long _String::FindTerminator (long start, _String const& terminator) const {
  
  long    current_position  = start;
  
  
  long   curly_depth = 0L,
         square_depth = 0L,
         paren_depth = 0L;
  
  bool   do_escape = false;
  char   quote_state = '\0';
  
  while (current_position < s_length) {
    char this_char = s_data[current_position];
    if (do_escape) {
      do_escape = false;
    } else {
      if ((this_char == '"' || this_char == '\'') && !do_escape) {
        if (quote_state == '\0') {
          quote_state = this_char;
        } else {
          if (this_char == quote_state) {
            quote_state = '\0';
          }
        }
      } else {
        if (quote_state == '\0') {
          
          switch (this_char) {
            case '(':
              paren_depth ++;
              break;
            case ')':
              if (paren_depth > 0L) {
                paren_depth --;
              }
              break;
            case '[':
              square_depth++;
              break;
            case ']':
              if (square_depth > 0L) {
                square_depth --;
              }
              break;
            case '{':
              curly_depth++;
              break;
            case '}':
              if (curly_depth > 0L) {
                curly_depth --;
              }
              break;
            default:
              if (curly_depth == 0L && square_depth == 0L && paren_depth == 0L) {
                if (BeginsWith (terminator, true, current_position)) {
                  return current_position;
                }
              }
          }
        } else {
          if (this_char == '\\' && quote_state != '\0' && !do_escape) {
            do_escape = true;
          }
        }
      }
    }
    current_position++;
  }
  
  return kNotFound;
}

//=============================================================

bool _String::StripQuotes(char open_char, char close_char) {
  if (s_length >= 2UL) {
    if (s_data [0] == open_char && s_data [s_length - 1UL] == close_char) {
      Trim (1, s_length - 2UL);
        return true;
    }
  }
    return false;
}


//=============================================================

bool _String::IsValidIdentifier(int options) const {
  return s_length > 0UL && _IsValidIdentifierAux (options & fIDAllowCompound, options & fIDAllowFirstNumeric) == s_length - 1UL && hyReservedWords.FindObject (this) == kNotFound;
}

//=============================================================

const _String  _String::ConvertToAnIdent(int options) const {
  _StringBuffer converted;
  
  const char default_placeholder = '_';

  char       last_char = '\0';
  bool       allow_compounds = options | fIDAllowCompound,
             allow_first_numeric = options | fIDAllowFirstNumeric;
  
  
  unsigned long current_index = 0UL;
  
  bool          first     = true;
  
  for (; current_index < s_length; current_index ++) {
    char current_char =  s_data[current_index];
    if (first) {
      if ( hy_Valid_ID_Chars.is_valid_first (current_char) || allow_first_numeric && hy_Valid_ID_Chars.is_valid (current_char)) {
        converted << current_char;
      } else {
        if (last_char != default_placeholder) {
          converted << default_placeholder;
        }
      }
      first = false;
    } else {
      if ( ! hy_Valid_ID_Chars.is_valid(current_char)) {
        if (allow_compounds && current_char == '.') {
          first = true;
          converted << current_char;
        } else {
          if (last_char != default_placeholder) {
            converted << default_placeholder;
          }
        }
      }
    }
    last_char = converted.char_at (converted.length () - 1UL);
  }
  
  return converted;
}


//=============================================================

long _String::_IsValidIdentifierAux(bool allow_compounds, bool allow_first_numeric, char wildcard) const {
  
  unsigned long current_index = 0UL;
  
  bool          first     = true;
  
  for (; current_index < s_length; current_index ++) {
    char current_char =  s_data[current_index];
    if (first) {
      if ( ! (hy_Valid_ID_Chars.is_valid_first (current_char) || allow_first_numeric && hy_Valid_ID_Chars.is_valid (current_char))) {
        break;
      }
      first = false;
    } else {
      if ( ! hy_Valid_ID_Chars.is_valid(current_char)) {
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
  
  return kNotFound;
}

//=============================================================

bool    _String::IsALiteralArgument (bool strip_quotes) {
  if (s_length >= 2UL) {
    char quotes [2] = {'"','\''};
    for (char quote : quotes) {
      long from = 0L,
           to = ExtractEnclosedExpression (from,quote,quote, fExtractRespectEscape);
      
      if (from == 0L && to == s_length - 1L) {
        if (strip_quotes){
          Trim (1L, s_length-2L);
        }
        return true;
      }
    }
  }
  return false;
}

//=============================================================


hy_reference_type _String::ProcessVariableReferenceCases (_String& referenced_object, _String const * context) const {
  char first_char    = get_char(0);
  bool is_func_ref   = get_char(s_length-1) == '&';
  
  if (first_char == '*' || first_char == '^') {
    if (is_func_ref) {
      referenced_object = kEmptyString;
      return kStringInvalidReference;
    }
    bool is_global_ref = first_char == '^';
    _String   plain_name (*this, 1, -1);
    
    if (plain_name.IsValidIdentifier(fIDAllowCompound)) {
      if (context) {
        plain_name = *context & '.' & plain_name;
      }
      _FString * dereferenced_value = (_FString*)FetchObjectFromVariableByType(&plain_name, STRING);
      if (dereferenced_value && dereferenced_value->get_str().ProcessVariableReferenceCases (referenced_object) == kStringDirectReference) {
        if (!is_global_ref && context) {
          referenced_object = *context & '.' & referenced_object;
        }
        return is_global_ref?kStringGlobalDeference:kStringLocalDeference;
      }
    } else {
      
      _String try_as_expression;
      if (context) {
        _VariableContainer ctxt (*context);
        try_as_expression = ProcessLiteralArgument (&plain_name, &ctxt);
      } else {
        try_as_expression = ProcessLiteralArgument (&plain_name, nil);
      }
      if (try_as_expression.ProcessVariableReferenceCases (referenced_object) == kStringDirectReference) {
        if (!is_global_ref && context) {
          referenced_object = *context & '.' & try_as_expression;
        }
        
        return is_global_ref?kStringGlobalDeference:kStringLocalDeference;
      }
    }
  }
  
  if (is_func_ref) {
    referenced_object = Cut (0, s_length-2UL);
    if (referenced_object.IsValidIdentifier(fIDAllowCompound)) {
      referenced_object = (context? (*context & '.' & referenced_object): (referenced_object)) & '&';
      return kStringDirectReference;
    }
  }
  else {
    if (IsValidIdentifier(fIDAllowCompound)) {
      if (context) {
        _String cdot = *context & '.';
        referenced_object = BeginsWith(cdot) ? *this : (cdot & *this);
      } else {
        referenced_object = *this;
      }
      return kStringDirectReference;
    }
  }
  
  referenced_object = kEmptyString;
  return kStringInvalidReference;
}




/*
 ==============================================================
 Regular Expression Methods
 ==============================================================
 */

const _String _String::GetRegExpError(int error) {
  char buffer[512];
  buffer[regerror(error, nil, buffer, 511)] = 0;
  return _String("Regular Expression error:") & _String (buffer).Enquote();
}

//=============================================================

void _String::FlushRegExp(regex_t* re) {
  regfree(re);
  delete re;
}

//=============================================================

regex_t* _String::PrepRegExp(const _String& pattern, int &error_code, bool case_sensitive) {
  regex_t *res = new regex_t;
  
  error_code = regcomp(res, pattern.get_str(),
                    REG_EXTENDED | (case_sensitive ? 0 : REG_ICASE));
  
  if (error_code) {
    FlushRegExp(res);
    return nil;
  }
  return res;
}

//=============================================================

const _SimpleList _String::RegExpMatch(regex_t const* re, unsigned long start ) const {
  _SimpleList matched_pairs;
  
  if (s_length && start < s_length) {
    regmatch_t *matches = new regmatch_t[re->re_nsub + 1];
    int error_code = regexec(re, s_data + start, re->re_nsub + 1, matches, 0);
    if (error_code == 0) {
      for (long k = 0L; k <= re->re_nsub; k++) {
        matched_pairs << matches[k].rm_so + start
                      << matches[k].rm_eo - 1 + start;
      }
    }
    delete[] matches;
  }
  
  return matched_pairs;
}

//=============================================================


const _SimpleList _String::RegExpAllMatches(regex_t const* re) const {
  _SimpleList matched_pairs;
  
  if (s_length) {
    
    regmatch_t *matches = new regmatch_t[re->re_nsub + 1];
    int error_code = regexec(re, s_data, re->re_nsub + 1, matches, 0);
    while (error_code == 0) {
      long offset = matched_pairs.countitems()
      ? matched_pairs.Element(-1) + 1
      : 0;
      
      matched_pairs << matches[0].rm_so + offset
                    << matches[0].rm_eo - 1 + offset;
      
      offset += matches[0].rm_eo;
      if (offset < s_length) {
        error_code = regexec(re, s_data + offset, re->re_nsub + 1, matches, 0);
      } else {
        break;
      }
    }
    delete[] matches;
  }
  return matched_pairs;
}

//=============================================================

const _SimpleList _String::_IntRegExpMatch (const _String & pattern,
                                           bool case_sensitive, bool handle_errors, bool match_all) const {
  if (s_length) {
    int err_code = 0;
    regex_t* regex = PrepRegExp(pattern, err_code, case_sensitive);
    if (regex) {
      _SimpleList hits = match_all ? RegExpAllMatches(regex) : RegExpMatch(regex);
      FlushRegExp(regex);
      return hits;
    } else if (handle_errors) {
      HandleApplicationError(GetRegExpError(err_code));
    }
  }
  return _SimpleList();
}

//=============================================================

const _SimpleList _String::RegExpMatch (const _String & pattern,
                                        bool case_sensitive, bool handle_errors) const {
  return _IntRegExpMatch (pattern, case_sensitive, handle_errors, false);
}

//=============================================================

const _SimpleList _String::RegExpAllMatches (const _String & pattern,
                                        bool case_sensitive, bool handle_errors) const {
  return _IntRegExpMatch (pattern, case_sensitive, handle_errors, true);
}

/*
==============================================================
Methods
==============================================================
*/


// Compute Adler-32 CRC for a string
// Implementation shamelessly lifted from http://en.wikipedia.org/wiki/Adler-32
long _String::Adler32(void) const {
  
  const static unsigned long  MOD_ADLER = 65521UL;

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
    a = (a & 0xffff) + (a >> 16) * (65536UL - MOD_ADLER);
    b = (b & 0xffff) + (b >> 16) * (65536UL - MOD_ADLER);
  }
  
  if (a >= MOD_ADLER) {
    a -= MOD_ADLER;
  }
  
  b = (b & 0xffff) + (b >> 16) * (65536UL - MOD_ADLER);
  
  if (b >= MOD_ADLER) {
    b -= MOD_ADLER;
  }
  
  return b << 16 | a;
}

//=============================================================


_String const _String::Random(const unsigned long length, const _String * alphabet) {
  _String random (length);
  
  unsigned long alphabet_length = alphabet?alphabet->s_length:127UL;
  
  if (length > 0UL && alphabet_length > 0UL) {
    for (unsigned long c = 0UL; c < length; c++) {
      unsigned long idx = genrand_int32 () % alphabet_length;
      if (alphabet) {
        random.set_char (c, alphabet->get_char(idx));
      } else {
        random.set_char (c,(char)(1UL+idx));
      }
    }
  }
  
  return random;
}

//=============================================================

unsigned long _String::LempelZivProductionHistory (_SimpleList* rec) const {
    if (rec) {
        rec->Clear();
    }

    if (empty()) {
        return 0UL;
    }

    if (rec) {
        (*rec) << 0;
    }

    unsigned long   current_position = 1UL,
                    production_history  = 1UL;

    while (current_position<s_length) {
      
        unsigned long max_extension = 0UL;

        for (unsigned long ip = 0; ip < current_position; ip++) {
            long sp = ip,
                 mp = current_position;

            while (mp < s_length && s_data [mp] == s_data[sp]) {
                mp++;
                sp++;
            }

            if (mp==s_length) {
                max_extension = s_length-current_position;
                break;
            } else {
                if ((mp = mp - current_position + 1UL) > max_extension) {
                    max_extension = mp;
                }
            }
        }

        current_position += max_extension;
      
        if (rec) {
            (*rec) << current_position - 1UL;
        } else {
            production_history ++;
        }
    }

    if (rec) {
        return rec->lLength;
    }

    return production_history;
}














