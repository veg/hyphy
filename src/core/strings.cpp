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

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "batchlan.h"
#include "function_templates.h"
#include "global_things.h"
#include "hy_string_buffer.h"

#include "hy_strings.h"
#include "mersenne_twister.h"
#include <memory>

_String compileDate = __DATE__;

using namespace hy_global;

struct _hy_Valid_ID_Chars_Type {
  unsigned char valid_chars[256];

  inline bool is_valid_first(unsigned char c) const {
    return valid_chars[c] == 2;
  }

  inline bool is_valid(unsigned char c) const { return valid_chars[c] > 0; }

  _hy_Valid_ID_Chars_Type(void) {
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
    valid_chars[(unsigned char)'_'] = 2;
  }
} hy_Valid_ID_Chars;

/*
==============================================================
Constructors/Destructors/Copiers
==============================================================
*/

_String::_String(void) { _String::Initialize(); }

//=============================================================

void _String::Initialize(bool) {
  s_length = 0UL;
  s_data = nil;
}

//=============================================================

void _String::Clear(void) {
  s_length = 0UL;
  if (s_data) {
    free(s_data);
    s_data = nil;
  }
}

//=============================================================

_String::_String(long const number) {
  char s[64];
  s_length = snprintf(s, sizeof(s), "%ld", number);
  AllocateAndCopyString(s, s_length);
}

//=============================================================

_String::_String(const unsigned long sL, char *buffer) {
  s_length = sL;
  if (buffer) {
    s_data = buffer;
    AddAReference();
  } else {
    s_data = (char *)MemAllocate(sL + 1L, true);
  }
  s_data[sL] = (char)0;
}

//=============================================================

_String::_String(const hyFloat val, const char *format) {
  char s_val[128];
  s_length = snprintf(s_val, 128, format ? format : PRINTF_FORMAT_STRING, val);
  AllocateAndCopyString(s_val, s_length);
}

//=============================================================

_String::_String(const hyFloat val, unsigned char digits) {
  char s_val[128];
  if (digits > 0) {
    s_length = snprintf(s_val, 128, "%.*g", MIN(digits, 20), val);
  } else {
    s_length = snprintf(s_val, 128, "%g", val);
  }
  AllocateAndCopyString(s_val, s_length);
}

//=============================================================

_String::_String(const _String &s) {
  _String::Initialize();
  _String::Duplicate(&s);
}

//=============================================================

_String::_String(_String &&s) {
  s_length = s.s_length;
  s_data = s.s_data;
  s.Initialize();
}

//=============================================================

_String::_String(_StringBuffer &&s) {
  s_length = s.s_length;
  s.TrimSpace();
  s_data = s.s_data;
  s._String::Initialize();
  s.Initialize();
}

//=============================================================

_String::_String(_String *s, bool dynamic) {
  if (s->CanFreeMe()) {
    s_data = s->s_data;
    s_length = s->s_length;
    s->s_data = nil;
    if (dynamic) {
      DeleteObject(s);
    }
  } else {
    AllocateAndCopyString(s->s_data, s->s_length);
    if (dynamic) {
      s->RemoveAReference();
    }
  }
}

//=============================================================
_String::_String(const _String &source, long start, long end) {
  if (source.s_length) {

    long requested_range = source.NormalizeRange(start, end);

    if (requested_range > 0L) {
      AllocateAndCopyString(source.s_data + start, requested_range);
      return;
    }
  }

  s_length = 0UL;
  s_data = (char *)MemAllocate(1UL);
  s_data[0] = '\0';
}

//=============================================================

_String::_String(const char *c_string) {
  AllocateAndCopyString(c_string, strlen(c_string));
}

//=============================================================
_String::_String(const wchar_t *wc_string) {
  unsigned long allocated = wcslen(wc_string);
  s_length = 0UL;
  s_data = (char *)MemAllocate(allocated + 1UL);
  for (unsigned long cid = 0UL; cid < allocated; cid++) {
    int this_char = wctob(wc_string[cid]);
    if (this_char != WEOF) {
      s_data[s_length++] = (char)this_char;
    }
  }
  if (s_length != allocated) {
    s_data =
        (char *)MemReallocate((char *)s_data, (s_length + 1) * sizeof(char));
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
_String::_String(const _String &str, unsigned long copies) {
  s_length = copies * str.s_length;
  s_data = (char *)MemAllocate(s_length + 1UL);
  if (s_length > 0UL) {
    for (unsigned long i = 0UL; i < copies; i++) {
      memcpy(s_data + i * str.s_length, str.s_data, str.s_length);
    }
  }
  s_data[s_length] = '\0';
}

//=============================================================
_String::_String(FILE *file, long read_this_many) {
  _String::Initialize();
  if (file) {
    if (read_this_many < 0) {
      fseek(file, 0, SEEK_END);
      s_length = (unsigned long)ftell(file);
      rewind(file);
    } else {
      s_length = read_this_many;
    }
    s_data = (char *)MemAllocate(s_length + 1UL);
    unsigned long read_items = fread(s_data, 1, s_length, file);
    if (read_items < s_length) {
      s_data = (char *)MemReallocate(s_data, read_items + 1);
      s_length = read_items;
    }
    s_data[s_length] = '\0';
  }
}

//=============================================================
_String::_String(hyFile *file, long read_this_many) {
  _String::Initialize();
  if (file) {
    if (read_this_many < 0) {
      file->seek(0, SEEK_END);
      s_length = (unsigned long)file->tell();
      file->rewind();
    } else {
      s_length = read_this_many;
    }
    s_data = (char *)MemAllocate(s_length + 1UL);
    unsigned long read_items = file->read(s_data, 1, s_length);
    if (read_items < s_length) {
      s_data = (char *)MemReallocate(s_data, read_items + 1);
      s_length = read_items;
    }
    s_data[s_length] = '\0';
  }
}
//=============================================================
_String::~_String(void) {
  if (CanFreeMe()) {
    if (s_data) {
      free(s_data);
      s_data = nil;
    }
    s_length = 0UL;
  } else {
    RemoveAReference();
  }
}

//=============================================================
BaseRef _String::makeDynamic(void) const {
  _String *r = new _String;
  r->Duplicate(this);
  return r;
}

//=============================================================
void _String::Duplicate(BaseRefConst ref) {
  if (s_data) {
    free(s_data);
  }

  _String const *s = (_String const *)ref;

  s_length = s->s_length;
  s_data = s->s_data;

  if (s_data) {
    AllocateAndCopyString(s->s_data, s_length);
  }
}

//=============================================================
_String &_String::operator=(_String const &s) {
  if (&s != this)
    Duplicate(&s);
  return *this;
}

//=============================================================
_String &_String::operator=(_String &&rhs) {
  if (this != &rhs) {
    if (s_data) {
      free(s_data);
    }
    s_data = rhs.s_data;
    s_length = rhs.s_length;
    rhs.s_data = nil;
  }
  return *this;
}

/*
 ==============================================================
 Private helpers
 ==============================================================
 */

/**
 * @brief Normalizes a given range [from, to] to be within the bounds of the
 * string.
 * @param from The starting index of the range (passed by reference and modified
 * in place).
 * @param to The ending index of the range (passed by reference and modified in
 * place).
 * @return The length of the normalized range. Returns 0 if the string is empty.
 * @note This function modifies the 'from' and 'to' parameters.
 * @example
 * @code
 *   _String my_string("hello world");
 *   long from = -2, to = 20;
 *   long length = my_string.NormalizeRange(from, to);
 *   // from is now 0, to is now 10, and length is 11.
 * @endcode
 */
long _String::NormalizeRange(long &from, long &to) const {

  if (s_length == 0UL) {
    return 0L;
  }

  if (from < 0L) {
    from = 0L;
  }

  if (to < 0L || to >= (long)s_length) {
    to = s_length - 1UL;
  }

  return to - from + 1L;
}

//=============================================================

/**
 * @brief Allocates memory for the string and copies the source string into it.
 * @param source_string The C-style string to copy.
 * @param length The length of the source string.
 */
void _String::AllocateAndCopyString(const char *source_string,
                                    unsigned long length) {
  s_length = length;
  s_data = (char *)MemAllocate(length + 1UL);
  // if (s_length) {
  memcpy(s_data, source_string, length);
  //}
  s_data[length] = '\0';
}

/*
==============================================================
Getters and setters
==============================================================
*/

/**
 * @brief Access a character in the string by its index.
 * @param index The index of the character to access.
 * @return A reference to the character at the specified index.
 * @throws HandleApplicationError if the index is out of bounds.
 * @example
 * @code
 *  _String s ("test");
 *  s[0] = 'T'; // s is now "Test"
 * @endcode
 */
char &_String::operator[](long index) {
  if (index < (long)s_length && index >= 0L) {
    return s_data[index];
  }
  HandleApplicationError(_String("Internal error at ") & __PRETTY_FUNCTION__ &
                         ": an invalid index requested");
  return s_data[0];
}

//=============================================================

/**
 * @brief Access a character in the string by its index (const version).
 * @param index The index of the character to access. Can be negative to index
 * from the end of the string.
 * @return The character at the specified index. Returns `default_return` if the
 * index is out of bounds.
 * @example
 * @code
 *  _String s ("test");
 *  char c = s(0);  // c is 't'
 *  char d = s(-1); // d is 't'
 * @endcode
 */
char _String::operator()(long index) const {
  if (index >= 0L && index < (long)s_length) {
    return s_data[index];
  }
  if (index < 0L && -index <= (long)s_length) {
    return s_data[s_length + index];
  }
  return default_return;
}

//=============================================================

/**
 * @brief Sets a character at a specific index.
 * @param index The index of the character to set.
 * @param data The new character.
 * @note This function does nothing if the index is out of bounds.
 */
void _String::set_char(unsigned long index, char const data) {
  if (index < s_length) {
    s_data[index] = data;
  }
}

//=============================================================

/**
 * @brief Sets a character at a specific index without checking bounds.
 * @param index The index of the character to set.
 * @param data The new character.
 * @warning This function does not perform bounds checking. Use with caution.
 */
void _String::set_char_no_check(unsigned long index, char const data) {
  s_data[index] = data;
}
//=============================================================

/**
 * @brief Get the underlying C-style string.
 * @return A const pointer to the C-style string.
 */
const char *_String::get_str(void) const { return s_data; }

/*
 ==============================================================
 Type conversions
 ==============================================================
 */

/**
 * @brief Implicitly convert the _String to a const char*.
 * @return A const pointer to the C-style string.
 */
_String::operator const char *(void) const { return s_data; }

//=============================================================

/**
 * @brief Converts the string to a floating-point number.
 * @return The float value of the string. Returns 0.0 if the string is empty.
 */
hyFloat _String::to_float(void) const {
  if (s_length == 0UL) {
    return 0.;
  }
  char *endP;
  return strtod(s_data, &endP);
}

//=============================================================

/**
 * @brief Converts the string to a long integer.
 * @return The long integer value of the string. Returns 0 if the string is
 * empty.
 */
long _String::to_long(void) const {
  if (s_length == 0UL) {
    return 0L;
  }
  char *endP;
  return strtol(s_data, &endP, 10);
}

//=============================================================

/**
 * @brief Converts the object to a string representation.
 * @param arg Unused.
 * @return A reference to the string object itself.
 */
BaseRef _String::toStr(unsigned long) {
  AddAReference();
  return this;
}

//=============================================================

/**
 * @brief Formats a time difference in seconds into a HH:MM:SS string.
 * @param time_diff The time difference in seconds.
 * @return A _String object representing the formatted time.
 * @example
 * @code
 *  _String formatted_time = _String::FormatTimeString(3661);
 *  // formatted_time is "01:01:01"
 * @endcode
 */
const _String _String::FormatTimeString(long time_diff) {

  long fields[3] = {time_diff / 3600L, time_diff / 60L % 60L, time_diff % 60L};

  _StringBuffer time_string;

  for (unsigned long l = 0; l < 3UL; l++) {
    if (l) {
      time_string << ':';
    }
    if (fields[l] < 10L) {
      time_string << '0';
    }
    time_string << _String(fields[l]);
  }

  return time_string;
}

/*
 ==============================================================
 Comparisons
 ==============================================================
 */

/**
 * @brief Compare this string with another string.
 * @param rhs The string to compare against.
 * @return hyComparisonType indicating the relationship between the strings
 * (kCompareLess, kCompareEqual, kCompareGreater).
 */
hyComparisonType _String::Compare(_String const &rhs) const {

  unsigned long min_len = MIN(s_length, rhs.s_length);
  int result = memcmp(s_data, rhs.s_data, min_len);

  if (result < 0) {
    return kCompareLess;
  }
  if (result > 0) {
    return kCompareGreater;
  }

  if (s_length < rhs.s_length) {
    return kCompareLess;
  }
  if (s_length > rhs.s_length) {
    return kCompareGreater;
  }

  return kCompareEqual;
}

//=============================================================

/**
 * @brief Compare this string with another string, ignoring case.
 * @param rhs The string to compare against.
 * @return hyComparisonType indicating the relationship between the strings
 * (kCompareLess, kCompareEqual, kCompareGreater).
 */
hyComparisonType _String::CompareIgnoringCase(_String const &rhs) const {
  unsigned long up_to = MIN(s_length, rhs.s_length);

  for (unsigned long i = 0UL; i < up_to; i++) {

    char llhs = (char)tolower(s_data[i]), lrhs = (char)tolower(rhs.s_data[i]);

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

/** @brief Equality operator. */
bool _String::operator==(const _String &s) const {
  return Compare(s) == kCompareEqual;
}
/** @brief Greater than operator. */
bool _String::operator>(const _String &s) const {
  return Compare(s) == kCompareGreater;
}
/** @brief Less than or equal to operator. */
bool _String::operator<=(const _String &s) const {
  return Compare(s) != kCompareGreater;
}
/** @brief Greater than or equal to operator. */
bool _String::operator>=(const _String &s) const {
  return Compare(s) != kCompareLess;
}
/** @brief Inequality operator. */
bool _String::operator!=(const _String &s) const {
  return Compare(s) != kCompareEqual;
}
/** @brief Less than operator. */
bool _String::operator<(const _String &s) const {
  return Compare(s) == kCompareLess;
}

/** @brief Check for equality with another string. */
bool _String::Equal(const _String &s) const {
  return Compare(s) == kCompareEqual;
};
/** @brief Check for equality with another string, ignoring case. */
bool _String::EqualIgnoringCase(const _String &s) const {
  return CompareIgnoringCase(s) == kCompareEqual;
};
/** @brief Check for equality with a single character. */
bool _String::Equal(const char c) const {
  return s_length == 1UL && s_data[0] == c;
}

//=============================================================

/**
 * @brief Check for equality with a pattern string containing a wildcard
 * character.
 * @param pattern The pattern string to match against.
 * @param wildchar The wildcard character.
 * @param start_this The starting index in this string.
 * @param start_pattern The starting index in the pattern string.
 * @param wildchar_matches A list to store the indices of wildcard matches.
 * @return True if the string matches the pattern, false otherwise.
 */
bool _String::EqualWithWildChar(const _String &pattern, const char wildchar,
                                unsigned long start_this,
                                unsigned long start_pattern,
                                _SimpleList *wildchar_matches) const {
  // wildcards only matter in the second string

  if (pattern.s_length > start_pattern && wildchar != '\0') {
    unsigned long match_this_char = start_pattern;
    // the position we are currently trying to match in the pattern

    bool is_wildcard = pattern.s_data[match_this_char] == wildchar,
         scanning_pattern = is_wildcard;

    unsigned long i = start_this;
    // the position we are currently trying to match in *this
    long last_matched_char = (long)start_this - 1L;
    // the index of the last character in *this that was matched to something
    // other than the wildcard

    while (i <= s_length) {
      if (scanning_pattern) { // skip consecutive wildcards in "pattern"
        scanning_pattern = pattern.s_data[++match_this_char] == wildchar;
      } else {
        if (s_data[i] == pattern.s_data[match_this_char]) {
          if (is_wildcard) {
            // could either match the next character or consume it into the
            // wildcard
            if (wildchar_matches) {
              // record the current wildcard match
              // if this is the last (0) character, return true
              long rollback_checkpoint = wildchar_matches->countitems();
              if (last_matched_char + 1 <
                  (long)i) { // something get matched to the wildchard
                *wildchar_matches << (last_matched_char + 1) << (i - 1);
              }
              if (i == s_length) {
                return true;
              }
              if (EqualWithWildChar(pattern, wildchar, i, match_this_char,
                                    wildchar_matches)) {
                // matching worked
                return true;
              } else { // consume the character into the wildcard
                i++;
                for (long k =
                         wildchar_matches->countitems() - rollback_checkpoint;
                     k >= 0; k--) {
                  wildchar_matches->Pop();
                }
                continue;
              }
            } else {
              if (EqualWithWildChar(pattern, wildchar, i, match_this_char)) {
                // matching worked
                return true;
              } else { // consume the character into the wildcard
                i++;
                continue;
              }
            }
          } else {
            // try character match
            // note that the terminal '0' characters will always match, so
            // this is where we terminate
            if (wildchar_matches) {
              if (last_matched_char + 1 <
                  (long)i) { // something get matched to the wildchard
                *wildchar_matches << (last_matched_char + 1) << (i - 1);
              }
              last_matched_char = i;
            }
            i++;
            match_this_char++;
            if (i > s_length || match_this_char > pattern.s_length) {
              break;
            }
            // TODO check to see if this will return true strings that match the
            // pattern and have some left-over stuff, like "tree.node.a.b" might
            // incorrectly match "tree.?.a"
            is_wildcard = pattern.s_data[match_this_char] == wildchar;
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

    if (wildchar_matches) {
      if (last_matched_char + 1 <
          (long)i) { // something get matched to the wildchard
        *wildchar_matches << (last_matched_char + 1) << (i - 1);
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

/**
 * @brief Append a string to this string.
 * @param rhs The string to append.
 * @return A new _String object containing the concatenated strings.
 * @example
 * @code
 *  _String s1 ("Hello, ");
 *  _String s2 ("World!");
 *  _String s3 = s1 & s2; // s3 is "Hello, World!"
 * @endcode
 */
_String _String::operator&(const _String &rhs) const {
  unsigned long combined_length = s_length + rhs.s_length;

  if (combined_length == 0UL) {
    return kEmptyString;
  }

  _String res(combined_length);

  if (s_length && s_data) {
    memcpy(res.s_data, s_data, s_length);
  }

  if (rhs.s_length && rhs.s_data) {
    memcpy(res.s_data + s_length, rhs.s_data, rhs.s_length);
  }

  res.s_data[res.s_length] = '\0';
  return res;
}

//=============================================================

/**
 * @brief Removes a substring from the string.
 * @param start The starting index of the substring to remove.
 * @param end The ending index of the substring to remove.
 * @return A new _String object with the specified substring removed.
 * @example
 * @code
 *  _String s ("Hello, World!");
 *  _String s2 = s.Chop(5, 11); // s2 is "Hello!"
 * @endcode
 */
_String _String::Chop(long start, long end) const {

  long resulting_length = NormalizeRange(start, end);

  if (resulting_length > 0L) {
    _String res((unsigned long)(s_length - resulting_length));
    if (start > 0L) {
      memcpy(res.s_data, s_data, start);
    }
    if (end + 1L < (long)s_length) {
      memcpy(res.s_data + start, s_data + end + 1L, s_length - end - 1L);
    }

    return res;
  }

  return *this;
}

//=============================================================

/**
 * @brief Extracts a substring from the string.
 * @param start The starting index of the substring.
 * @param end The ending index of the substring.
 * @return A new _String object containing the extracted substring.
 * @example
 * @code
 *  _String s ("Hello, World!");
 *  _String s2 = s.Cut(7, 11); // s2 is "World"
 * @endcode
 */
_String _String::Cut(long start, long end) const {
  return _String(*this, start, end);
}

//=============================================================

/**
 * @brief Deletes a substring from the string in place.
 * @param start The starting index of the substring to delete.
 * @param end The ending index of the substring to delete.
 * @example
 * @code
 *  _String s ("Hello, World!");
 *  s.Delete(5, 11); // s is now "Hello!"
 * @endcode
 */
void _String::Delete(long start, long end) {
  long resulting_length = NormalizeRange(start, end);

  if (resulting_length > 0L) {
    if (end < (long)s_length - 1) {
      memmove(s_data + start, s_data + end + 1L, s_length - end - 1L);
    }
    s_length -= resulting_length;
    s_data = (char *)MemReallocate(s_data, sizeof(char) * (s_length + 1UL));
    s_data[s_length] = '\0';
  }
}

//=============================================================

/**
 * @brief Reverses the string in place.
 * @example
 * @code
 *  _String s ("abcde");
 *  s.Flip(); // s is now "edcba"
 * @endcode
 */
void _String::Flip(void) {
  for (unsigned long i = 0UL; i < (s_length >> 1); i++) {
    char c;
    SWAP(s_data[i], s_data[s_length - 1 - i], c);
  }
}

//=============================================================

/**
 * @brief Returns a reversed copy of the string.
 * @return A new _String object that is the reverse of this string.
 */
_String _String::Reverse(void) const {

  _String result(*this);
  for (unsigned long s = 0UL, e = s_length - 1L; s < s_length; s++, e--) {
    result.s_data[s] = s_data[e];
  }
  return result;
}

//=============================================================

/**
 * @brief Inserts a character into the string at a specified position.
 * @param c The character to insert.
 * @param where The index at which to insert the character.
 * @example
 * @code
 *  _String s ("Hllo");
 *  s.Insert('e', 1); // s is now "Hello"
 * @endcode
 */
void _String::Insert(char c, long where) {
  if (where < 0L || where >= (long)s_length) {
    where = s_length;
  }

  s_data = (char *)MemReallocate(s_data, sizeof(char) * (s_length + 2UL));

  if (where < (long)s_length) {
    memmove(s_data + where + 1UL, s_data + where, s_length - where);
  }

  s_data[where] = c;
  s_data[++s_length] = '\0';
}

//=============================================================

/**
 * @brief Trims the string to a specified range in place.
 * @param start The starting index of the range to keep.
 * @param end The ending index of the range to keep.
 * @example
 * @code
 *  _String s ("  Hello  ");
 *  s.Trim(2, 6); // s is now "Hello"
 * @endcode
 */
void _String::Trim(long start, long end) {
  long resulting_length = NormalizeRange(start, end);
  /*if (s_length >= 5000 && start > 0) {
      printf ("\nLong trim %d %d %d\n", s_length, start, end);
  }*/

  if (resulting_length > 0L) {
    if (start > 0L) {
      memmove(s_data, s_data + start, resulting_length);
    }
    if ((long)s_length != resulting_length) {
      s_length = resulting_length;
      s_data = (char *)MemReallocate(s_data, resulting_length + 1UL);
      s_data[resulting_length] = '\0';
    }
  } else {
    s_length = 0UL;
    s_data = (char *)MemReallocate(s_data, 1UL);
    s_data[0] = '\0';
  }
}

//=============================================================

/**
 * @brief Changes the case of the string.
 * @param conversion_type The type of case conversion (kStringUpperCase or
 * kStringLowerCase).
 * @return A new _String object with the case changed.
 */
const _String _String::ChangeCase(hy_string_case conversion_type) const {
  _String result(s_length);

  auto conversion_function =
      conversion_type == kStringUpperCase ? toupper : tolower;

  for (unsigned long i = 0UL; i < s_length; i++) {
    result.s_data[i] = (char)conversion_function(s_data[i]);
  }

  return result;
}

//=============================================================

/**
 * @brief Changes the case of the string in place.
 * @param conversion_type The type of case conversion (kStringUpperCase or
 * kStringLowerCase).
 */
void _String::ChangeCaseInPlace(hy_string_case conversion_type) {

  auto conversion_function =
      conversion_type == kStringUpperCase ? toupper : tolower;

  for (unsigned long i = 0UL; i < s_length; i++) {
    s_data[i] = (char)conversion_function(s_data[i]);
  }
}

//=============================================================

/**
 * @brief Tokenizes the string based on a splitter string.
 * @param splitter The string to use as a delimiter.
 * @return A _List of _String objects.
 * @example
 * @code
 *  _String s ("a,b,c");
 *  _List tokens = s.Tokenize(",");
 *  // tokens contains ["a", "b", "c"]
 * @endcode
 */
const _List _String::Tokenize(const _String &splitter) const {
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

  tokenized < new _String(*this, cp, kStringEnd);
  return tokenized;
}

//=============================================================

/**
 * @brief Tokenizes the string based on a set of splitter characters.
 * @param splitter An array of booleans indicating which characters are
 * splitters.
 * @return A _List of _String objects.
 */
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

/**
 * @brief Encloses the string in a pair of identical quotes.
 * @param quote_char The character to use for quoting.
 * @return A new _String object with the quotes added.
 */
const _String _String::Enquote(char quote_char) const {
  return _StringBuffer(2UL + s_length) << quote_char << *this << quote_char;
}

//=============================================================

/**
 * @brief Encloses the string in a pair of opening and closing quotes.
 * @param open_char The opening quote character.
 * @param close_char The closing quote character.
 * @return A new _String object with the quotes added.
 */
const _String _String::Enquote(char open_char, char close_char) const {
  return _StringBuffer(2UL + s_length) << open_char << *this << close_char;
}

//=============================================================

/**
 * @brief Removes all whitespace characters from the string.
 * @return A new _String object with whitespace removed.
 */
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

/**
 * @brief Compresses consecutive whitespace characters into a single space.
 * @return A new _String object with whitespace compressed.
 */
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

/**
 * @brief Finds the first occurrence of a pattern in the string.
 * @param pattern The pattern to search for.
 * @param start The starting index for the search.
 * @param end The ending index for the search.
 * @return The index of the first character of the found pattern, or kNotFound
 * if not found.
 * @example
 * @code
 *  _String s ("hello world");
 *  long index = s.Find("world"); // index is 6
 * @endcode
 */
long _String::Find(const _String &pattern, long start, long end) const {

  if (pattern.s_length) {
    long span = NormalizeRange(start, end);
    if (span >= (long)pattern.s_length) {
      const unsigned long upper_bound = end - pattern.s_length + 2L;
      for (unsigned long i = start; i < upper_bound; i++) {
        unsigned long j = 0UL;
        for (; j < pattern.s_length; j++) {
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

/**
 * @brief Finds the last occurrence of a pattern in the string.
 * @param pattern The pattern to search for.
 * @param start The starting index for the search.
 * @param end The ending index for the search.
 * @return The index of the first character of the found pattern, or kNotFound
 * if not found.
 * @example
 * @code
 *  _String s ("hello world world");
 *  long index = s.FindBackwards("world"); // index is 12
 * @endcode
 */
long _String::FindBackwards(const _String &pattern, long start,
                            long end) const {

  if (pattern.s_length) {
    long span = NormalizeRange(start, end);
    if (span >= (long)pattern.s_length) {
      const long upper_bound = end - pattern.s_length + 1L;
      for (long i = upper_bound; i >= start; i--) {
        unsigned long j = 0UL;
        for (; j < pattern.s_length; j++) {
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

/**
 * @brief Finds the first occurrence of a character in the string.
 * @param p The character to search for.
 * @param start The starting index for the search.
 * @param end The ending index for the search.
 * @return The index of the found character, or kNotFound if not found.
 */
long _String::Find(const char p, long start, long end) const {
  if (s_length) {
    long span = NormalizeRange(start, end);
    if (span > 0L) {
      char sentinel = s_data[end + 1];
      s_data[end + 1L] = p;
      long index = start;
      while (s_data[index] != p) {
        index++;
      }
      s_data[end + 1L] = sentinel;
      return index <= end ? index : kNotFound;
    }
  }

  return kNotFound;
}

//=============================================================

/**
 * @brief Finds the first occurrence of any character from a given set.
 * @param lookup A boolean array representing the set of characters to search
 * for.
 * @param start The starting index for the search.
 * @param end The ending index for the search.
 * @return The index of the first character found, or kNotFound if not found.
 */
long _String::Find(const bool lookup[256], long start, long end) const {
  if (s_length) {
    long span = NormalizeRange(start, end);
    if (span > 0L) {

      for (long index = start; index <= end; index++) {
        if (lookup[(unsigned char)s_data[index]]) {
          return index;
        }
      }
    }
  }

  return kNotFound;
}

//=============================================================

/**
 * @brief Finds the first occurrence of any character from a given set, ignoring
 * case.
 * @param lookup A boolean array representing the set of characters to search
 * for (should be lowercase).
 * @param start The starting index for the search.
 * @param end The ending index for the search.
 * @return The index of the first character found, or kNotFound if not found.
 */
long _String::FindAnyCase(const bool lookup[256], long start, long end) const {
  if (s_length) {
    long span = NormalizeRange(start, end);
    if (span > 0L) {

      for (long index = start; index <= end; index++) {
        if (lookup[tolower(s_data[index])] || lookup[toupper(s_data[index])]) {
          return index;
        }
      }
    }
  }

  return kNotFound;
}

//=============================================================
/**
 * @brief Finds the first occurrence of a pattern in the string, ignoring case.
 * @param pattern The pattern to search for.
 * @param start The starting index for the search.
 * @param end The ending index for the search.
 * @return The index of the first character of the found pattern, or kNotFound
 * if not found.
 */
long _String::FindAnyCase(const _String &pattern, long start, long end) const {

  if (pattern.s_length) {
    long span = NormalizeRange(start, end);
    if (span >= (long)pattern.s_length) {
      const unsigned long upper_bound = end - pattern.s_length + 2L;
      for (unsigned long i = start; i < upper_bound; i++) {
        unsigned long j = 0UL;
        for (; j < pattern.s_length; j++) {
          if (toupper(s_data[i + j]) != toupper(pattern.s_data[j])) {
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

/**
 * @brief Replaces occurrences of a pattern with a replacement string.
 * @param pattern The pattern to search for.
 * @param replace The string to replace the pattern with.
 * @param replace_all If true, all occurrences are replaced; otherwise, only the
 * first.
 * @return A new _String object with the replacements made.
 * @example
 * @code
 *  _String s ("hello world world");
 *  _String s2 = s.Replace("world", "earth", true); // s2 is "hello earth earth"
 * @endcode
 */
const _String _String::Replace(const _String &pattern, const _String &replace,
                               bool replace_all) const {

  if (s_length < pattern.s_length || pattern.s_length == 0UL) {
    return *this;
  }

  _StringBuffer replacement_buffer;
  unsigned long anchor_index = 0UL;
  for (; anchor_index <= s_length - pattern.s_length; anchor_index++) {
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
        anchor_index++;
        break;
      }
    } else {
      replacement_buffer << s_data[anchor_index];
    }
  }

  return replacement_buffer.AppendSubstring(*this, anchor_index, kNotFound);
}
//=============================================================

/**
 * @brief Finds the index of the first non-whitespace character.
 * @param start The starting index for the search.
 * @param end The ending index for the search.
 * @param direction The direction to search in (forward or backward).
 * @return The index of the first non-whitespace character, or kNotFound if not
 * found.
 */
long _String::FirstNonSpaceIndex(long start, long end,
                                 hy_string_search_direction direction) const {
  return _FindFirstIndexCondtion(start, end, direction,
                                 [](char c) -> bool { return !isspace(c); });
}

/**
 * @brief Finds the first non-whitespace character.
 * @param start The starting index for the search.
 * @param end The ending index for the search.
 * @param direction The direction to search in (forward or backward).
 * @return The first non-whitespace character, or `default_return` if not found.
 */
char _String::FirstNonSpace(long start, long end,
                            hy_string_search_direction direction) const {
  long r = FirstNonSpaceIndex(start, end, direction);
  return r == kNotFound ? default_return : s_data[r];
}

//=============================================================

/**
 * @brief Finds the index of the first whitespace character.
 * @param start The starting index for the search.
 * @param end The ending index for the search.
 * @param direction The direction to search in (forward or backward).
 * @return The index of the first whitespace character, or kNotFound if not
 * found.
 */
long _String::FirstSpaceIndex(long start, long end,
                              hy_string_search_direction direction) const {
  return _FindFirstIndexCondtion(start, end, direction, isspace);
}

//=============================================================

/**
 * @brief Finds the first non-whitespace character that follows a whitespace
 * character.
 * @param start The starting index for the search.
 * @param end The ending index for the search.
 * @param direction The direction to search in (forward or backward).
 * @return The index of the found character, or kNotFound if not found.
 */
long _String::FirstNonSpaceFollowingSpace(
    long start, long end, hy_string_search_direction direction) const {
  long first_space = FirstSpaceIndex(start, end, direction);
  if (first_space != kNotFound) {
    if (direction == kStringDirectionForward) {
      first_space = FirstNonSpaceIndex(first_space, end, direction);
    } else {
      first_space = FirstNonSpaceIndex(start, first_space, direction);
    }
  }
  return first_space;
}

/**
 * @brief Checks if the string begins with a given pattern.
 * @param pattern The pattern to check for.
 * @param case_sensitive If true, the comparison is case-sensitive.
 * @param startfrom The index to start the comparison from.
 * @return True if the string begins with the pattern, false otherwise.
 */
bool _String::BeginsWith(_String const &pattern, bool case_sensitive,
                         unsigned long startfrom) const {
  if (s_length >= pattern.s_length + startfrom) {
    if (case_sensitive) {
      for (unsigned long idx = 0; idx < pattern.s_length; idx++) {
        if (s_data[idx + startfrom] != pattern.s_data[idx]) {
          return false;
        }
      }
    } else {
      for (unsigned long idx = 0; idx < pattern.s_length; idx++) {
        if (tolower(s_data[idx + startfrom]) != tolower(pattern.s_data[idx])) {
          return false;
        }
      }
    }
    return true;
  }
  return false;
}

/**
 * @brief Checks if the string begins with any character from a given set.
 * @param pattern A boolean array representing the set of characters to check
 * for.
 * @param case_sensitive If true, the comparison is case-sensitive.
 * @param startfrom The index to start the comparison from.
 * @return True if the string begins with a character from the set, false
 * otherwise.
 */
bool _String::BeginsWith(const bool pattern[256], bool case_sensitive,
                         unsigned long startfrom) const {
  if (s_length >= 1UL + startfrom) {
    if (case_sensitive) {
      for (unsigned long idx = 0; idx < 256; idx++) {
        if (pattern[(unsigned char)s_data[startfrom]]) {
          return true;
        }
      }
    } else {
      for (unsigned long idx = 0; idx < 256; idx++) {
        if (pattern[tolower(s_data[startfrom])] ||
            pattern[toupper(s_data[startfrom])]) {
          return true;
        }
      }
    }
  }
  return false;
}

/**
 * @brief Checks if the string ends with a given pattern.
 * @param pattern The pattern to check for.
 * @param case_sensitive If true, the comparison is case-sensitive.
 * @return True if the string ends with the pattern, false otherwise.
 */
bool _String::EndsWith(_String const &pattern, bool case_sensitive) const {
  if (s_length >= pattern.s_length) {
    long length_difference = s_length - pattern.s_length;
    return (case_sensitive
                ? Find(pattern, length_difference)
                : FindAnyCase(pattern, length_difference)) == length_difference;
  }
  return false;
}

/**
 * @brief Checks if the string begins with a given pattern and is not followed
 * by an identifier character.
 * @param pattern The pattern to check for.
 * @return True if the condition is met, false otherwise.
 */
bool _String::BeginsWithAndIsNotAnIdent(_String const &pattern) const {

  if (BeginsWith(pattern)) {
    if (s_length > pattern.s_length) {
      char next_char = char_at(pattern.s_length);
      if (isalnum(next_char) || next_char == '.' || next_char == '_' ||
          next_char == '&') {
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

/**
 * @brief Removes matching opening and closing characters from the beginning and
 * end of the string.
 * @param open_char The opening character to strip.
 * @param close_char The closing character to strip.
 * @return True if the quotes were stripped, false otherwise.
 */
bool _String::StripQuotes(char open_char, char close_char) {
  if (s_length >= 2UL) {
    if (s_data[0] == open_char && s_data[s_length - 1UL] == close_char) {
      Trim(1, s_length - 2UL);
      return true;
    }
  }
  return false;
}

//=============================================================

/**
 * @brief Removes matching opening and closing characters from a set of possible
 * pairs.
 * @param open_chars A C-string of possible opening characters.
 * @param close_chars A C-string of corresponding closing characters.
 * @return True if any pair of quotes was stripped, false otherwise.
 */
bool _String::StripQuotes(char const *open_chars, char const *close_chars) {
  if (s_length >= 2UL) {
    size_t count = strlen(open_chars);
    for (size_t i = 0; i < count; i++) {
      if (s_data[0] == open_chars[i] &&
          s_data[s_length - 1UL] == close_chars[i]) {
        Trim(1, s_length - 2UL);
        return true;
      }
    }
  }
  return false;
}

//=============================================================

/**
 * @brief Checks if the string is a valid identifier.
 * @param options A bitmask of options for validation.
 * @return True if the string is a valid identifier, false otherwise.
 */
bool _String::IsValidIdentifier(int options) const {
  return s_length > 0UL &&
         _IsValidIdentifierAux(options & fIDAllowCompound,
                               options & fIDAllowFirstNumeric) ==
             (long)s_length - 1 &&
         hyReservedWords.FindObject(this) == kNotFound;
}

//=============================================================

/**
 * @brief Converts the string to a valid identifier.
 * @param options A bitmask of options for conversion.
 * @return A new _String object that is a valid identifier.
 */
const _String _String::ConvertToAnIdent(int options) const {
  _StringBuffer converted;

  const char default_placeholder = '_';

  char last_char = '\0';
  bool allow_compounds = options & fIDAllowCompound,
       allow_first_numeric = options & fIDAllowFirstNumeric;

  unsigned long current_index = 0UL;

  bool first = true;

  for (; current_index < s_length; current_index++) {
    char current_char = s_data[current_index];
    if (first) {
      if (hy_Valid_ID_Chars.is_valid_first(current_char) ||
          (allow_first_numeric && hy_Valid_ID_Chars.is_valid(current_char))) {
        converted << current_char;
      } else {
        if (last_char != default_placeholder) {
          converted << default_placeholder;
        }
      }
      first = false;
    } else {
      if (hy_Valid_ID_Chars.is_valid(current_char)) {
        converted << current_char;
      } else {
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
    last_char = converted.char_at(converted.length() - 1UL);
  }

  return converted;
}

//=============================================================

/**
 * @brief Auxiliary function to check for a valid identifier.
 * @param allow_compounds If true, allows compound identifiers (e.g., with
 * dots).
 * @param allow_first_numeric If true, allows the first character to be a
 * number.
 * @param unused Unused parameter.
 * @return The length of the valid identifier prefix, or kNotFound.
 */
long _String::_IsValidIdentifierAux(bool allow_compounds,
                                    bool allow_first_numeric, char) const {

  unsigned long current_index = 0UL;

  bool first = true;

  for (; current_index < s_length; current_index++) {
    char current_char = s_data[current_index];
    if (first) {
      if (!(hy_Valid_ID_Chars.is_valid_first(current_char) ||
            (allow_first_numeric &&
             hy_Valid_ID_Chars.is_valid(current_char)))) {
        break;
      }
      first = false;
    } else {
      if (!hy_Valid_ID_Chars.is_valid(current_char)) {
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

/**
 * @brief Checks if the string is a literal argument (i.e., enclosed in quotes).
 * @param strip_quotes If true, the quotes are removed from the string.
 * @return True if the string is a literal argument, false otherwise.
 */
bool _String::IsALiteralArgument(bool strip_quotes) {
  if (s_length >= 2UL) {
    char quotes[2] = {'"', '\''};
    for (char quote : quotes) {
      long from = 0L, to = ExtractEnclosedExpression(from, quote, quote,
                                                     fExtractRespectEscape);

      if (from == 0L && to == (long)s_length - 1L) {
        if (strip_quotes) {
          Trim(1L, s_length - 2L);
        }
        return true;
      }
    }
  }
  return false;
}

//=============================================================

/**
 * @brief Processes variable reference cases, resolving dereferences and
 * context.
 * @param referenced_object The string to store the resolved reference.
 * @param context The context to resolve the reference in.
 * @return The type of reference found.
 */
hy_reference_type
_String::ProcessVariableReferenceCases(_String &referenced_object,
                                       _String const *context) const {
  const static _String kDot(".");

  if (nonempty()) {

    char first_char = char_at(0);
    bool is_func_ref = char_at(s_length - 1) == '&';

    if (first_char == '*' || first_char == '^') {
      if (is_func_ref) {
        referenced_object = kEmptyString;
        return kStringInvalidReference;
      }
      bool is_global_ref = first_char == '^';
      _StringBuffer plain_name(length());
      plain_name.AppendSubstring(*this, 1, -1);

      if (plain_name.IsValidIdentifier(fIDAllowCompound |
                                       fIDAllowFirstNumeric)) {
        if (context) {
          plain_name.Clear();
          (plain_name << *context << '.').AppendSubstring(*this, 1, -1);
        }
        _FString *dereferenced_value =
            (_FString *)FetchObjectFromVariableByType(&plain_name, STRING);
        if (dereferenced_value &&
            dereferenced_value->get_str().ProcessVariableReferenceCases(
                referenced_object) == kStringDirectReference) {
          if (!is_global_ref && context) {
            referenced_object = (_StringBuffer(context->length() + 1UL +
                                               referenced_object.length())
                                 << *context << '.' << referenced_object);
          }
          return is_global_ref ? kStringGlobalDeference : kStringLocalDeference;
        }
      } else {

        _String try_as_expression;
        if (context) {
          _VariableContainer ctxt(*context);
          try_as_expression = ProcessLiteralArgument(&plain_name, &ctxt);
        } else {
          try_as_expression = ProcessLiteralArgument(&plain_name, nil);
        }
        if (try_as_expression.ProcessVariableReferenceCases(
                referenced_object) == kStringDirectReference) {
          if (!is_global_ref && context) {
            // referenced_object = *context & '.' & try_as_expression;
            referenced_object = (_StringBuffer(context->length() + 1UL +
                                               try_as_expression.length())
                                 << *context << '.' << try_as_expression);
          }

          return is_global_ref ? kStringGlobalDeference : kStringLocalDeference;
        }
      }
    }

    if (is_func_ref) {
      referenced_object = Cut(0, s_length - 2UL);
      if (referenced_object.IsValidIdentifier(fIDAllowCompound |
                                              fIDAllowFirstNumeric)) {
        referenced_object = (context ? (*context & '.' & referenced_object)
                                     : (referenced_object)) &
                            '&';
        return kStringDirectReference;
      }
    } else {
      if (IsValidIdentifier(fIDAllowCompound | fIDAllowFirstNumeric)) {
        if (context) {
          if (BeginsWith(*context) &&
              BeginsWith(kDot, true, context->length())) {
            referenced_object = *this;
          } else {
            referenced_object = (_StringBuffer(context->length() + length() + 1)
                                 << *context << '.' << *this);
          }

          //_String cdot = *context & '.';
          // referenced_object = BeginsWith(cdot) ? *this : (cdot & *this);
        } else {
          referenced_object = *this;
        }
        return kStringDirectReference;
      }
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

/**
 * @brief Gets a formatted error string for a regular expression error code.
 * @param error The error code.
 * @return A _String object containing the formatted error message.
 */
const _String _String::GetRegExpError(int error) {
  return _String("Regular Expression error: ") &
         _String(std::to_string(error).c_str()).Enquote();
}

//=============================================================

/**
 * @brief Deletes a compiled regular expression object.
 * @param re A pointer to the std::regex object to delete.
 */

//=============================================================

std::regex *_String::PrepRegExp(const _String &pattern, int &error_code,
                                bool case_sensitive, bool throw_errors) {

  std::regex_constants::syntax_option_type flags =
      std::regex_constants::ECMAScript;
  if (!case_sensitive) {
    flags |= std::regex_constants::icase;
  }

  try {
    std::regex *res = new std::regex(pattern.get_str(), flags);
    error_code = 0;
    return res;
  } catch (const std::regex_error &e) {
    error_code = e.code();
    if (throw_errors) {
      throw(_String(e.what()));
    }
    return nil;
  }
}

//=============================================================

/**
 * @brief Finds the first match of a compiled regular expression in the string.
 * @param re A pointer to the compiled std::regex object.
 * @param start The starting index for the search.
 * @return A _SimpleList of pairs, where each pair is the start and end index of
 * a match group.
 */
const _SimpleList _String::RegExpMatch(std::regex const *re,
                                       unsigned long start) const {
  _SimpleList matched_pairs;

  if (s_length && start < s_length) {
    std::cmatch static_matches;
    if (std::regex_search(s_data + start, static_matches, *re)) {
      for (auto const &match : static_matches) {
        matched_pairs << match.first - s_data << (match.second - s_data) - 1;
      }
    }
  }

  return matched_pairs;
}

//=============================================================

/**
 * @brief Finds all matches of a compiled regular expression in the string.
 * @param re A pointer to the compiled std::regex object.
 * @return A _SimpleList of pairs, where each pair is the start and end index of
 * a match.
 */
const _SimpleList _String::RegExpAllMatches(std::regex const *re) const {
  _SimpleList matched_pairs;

  if (s_length) {

    std::cregex_iterator it(s_data, s_data + s_length, *re);
    std::cregex_iterator end;

    while (it != end) {
      matched_pairs << it->position() << (it->position() + it->length() - 1);
      ++it;
    }
  }
  return matched_pairs;
}

//=============================================================

/**
 * @brief Internal helper function for regular expression matching.
 * @param pattern The regular expression pattern.
 * @param case_sensitive If true, the match is case-sensitive.
 * @param handle_errors If true, application errors are handled.
 * @param match_all If true, all matches are found; otherwise, only the first.
 * @return A _SimpleList of matching indices.
 */
const _SimpleList _String::_IntRegExpMatch(const _String &pattern,
                                           bool case_sensitive,
                                           bool handle_errors,
                                           bool match_all) const {
  if (s_length) {
    int err_code = 0;
    std::regex *regex = PrepRegExp(pattern, err_code, case_sensitive);
    if (regex) {
      _SimpleList hits =
          match_all ? RegExpAllMatches(regex) : RegExpMatch(regex);
      FlushRegExp(regex);
      return hits;
    } else if (handle_errors) {
      HandleApplicationError(
          _String("Regular Expression error: ") &
          _String(std::to_string(err_code).c_str()).Enquote());
    }
  }
  return _SimpleList();
}

//=============================================================

/**
 * @brief Finds the first match of a regular expression pattern in the string.
 * @param pattern The regular expression pattern.
 * @param case_sensitive If true, the match is case-sensitive.
 * @param handle_errors If true, application errors are handled.
 * @return A _SimpleList of pairs for the first match and its capture groups.
 */
const _SimpleList _String::RegExpMatch(const _String &pattern,
                                       bool case_sensitive,
                                       bool handle_errors) const {
  return _IntRegExpMatch(pattern, case_sensitive, handle_errors, false);
}

//=============================================================

/**
 * @brief Finds all matches of a regular expression pattern in the string.
 * @param pattern The regular expression pattern.
 * @param case_sensitive If true, the match is case-sensitive.
 * @param handle_errors If true, application errors are handled.
 * @return A _SimpleList of pairs for all matches.
 */
const _SimpleList _String::RegExpAllMatches(const _String &pattern,
                                            bool case_sensitive,
                                            bool handle_errors) const {
  return _IntRegExpMatch(pattern, case_sensitive, handle_errors, true);
}

//=============================================================

void _String::FlushRegExp(std::regex *re) { delete re; }

/*
==============================================================
Methods
==============================================================
*/

/**
 * @brief Computes the Adler-32 checksum of the string.
 * @return The Adler-32 checksum as a long integer.
 */
long _String::Adler32(void) const {

  const static unsigned long MOD_ADLER = 65521UL;

  unsigned long len = s_length, a = 1UL, b = 0UL, i = 0UL;

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

/**
 * @brief Generates a random string of a given length.
 * @param length The desired length of the random string.
 * @param alphabet An optional string specifying the set of characters to use.
 * If nil, characters from 1 to 127 are used.
 * @return A new _String object containing the random string.
 */
_String const _String::Random(const unsigned long length,
                              const _String *alphabet) {
  _String random(length);

  unsigned long alphabet_length = alphabet ? alphabet->s_length : 127UL;

  if (length > 0UL && alphabet_length > 0UL) {
    for (unsigned long c = 0UL; c < length; c++) {
      unsigned long idx = genrand_int32() % alphabet_length;
      if (alphabet) {
        random.set_char(c, alphabet->char_at(idx));
      } else {
        random.set_char(c, (char)(1UL + idx));
      }
    }
  }

  return random;
}

//=============================================================

/**
 * @brief Computes the Lempel-Ziv production history of the string.
 * @param rec An optional pointer to a _SimpleList to store the production
 * history.
 * @return The number of production steps.
 */
unsigned long _String::LempelZivProductionHistory(_SimpleList *rec) const {
  if (rec) {
    rec->Clear();
  }

  if (empty()) {
    return 0UL;
  }

  if (rec) {
    (*rec) << 0;
  }

  unsigned long current_position = 1UL, production_history = 1UL;

  while (current_position < s_length) {

    unsigned long max_extension = 0UL;

    for (unsigned long ip = 0; ip < current_position; ip++) {
      long sp = ip, mp = current_position;

      while (mp < (long)s_length && s_data[mp] == s_data[sp]) {
        mp++;
        sp++;
      }

      if (mp == (long)s_length) {
        max_extension = s_length - current_position;
        break;
      } else {
        if ((mp = mp - current_position + 1) > (long)max_extension) {
          max_extension = mp;
        }
      }
    }

    current_position += max_extension;

    if (rec) {
      (*rec) << current_position - 1UL;
    } else {
      production_history++;
    }
  }

  if (rec) {
    return rec->lLength;
  }

  return production_history;
}
