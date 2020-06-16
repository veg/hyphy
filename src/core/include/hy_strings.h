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

#ifndef _HY_STRINGS_
#define _HY_STRINGS_

#include "baseobj.h"
#include "hy_types.h"
#include "regex.h"
#include "wchar.h"

#define fExtractRespectQuote  0x01
#define fExtractRespectEscape 0x02
#define fExtractOneLevelOnly  0x04

#define fIDAllowFirstNumeric 0x01
#define fIDAllowCompound 0x02

#define kStringEnd (-1L)

enum hy_reference_type {
  kStringInvalidReference = 0x00,
  kStringDirectReference = 0x01,
  kStringLocalDeference = 0x02,
  kStringGlobalDeference = 0x03
};

enum hy_string_case { kStringUpperCase, kStringLowerCase };

enum hy_string_search_direction {
  kStringDirectionForward,
  kStringDirectionBackward
};

class _SimpleList;
class _List;
class _ExecutionList;
class _StringBuffer;

class _String : public BaseObj {
    
protected:
    char          *s_data;
    unsigned long s_length;
    
    /** this value is returned for "failed"
     access operations that don't throw errors, e.g. getChar */
    const static char default_return = '\0';

public:
    
    
  /*
   ==============================================================
   Constructors/Destructors/Copiers
   ==============================================================
   */

  /**
   * The default constuctor
   * which creates an empty string

   * Revision history
   - SLKP 20170517 porting from v3 branch
   */
  _String(void); // v3;

  /**
   * Standard initalization to 0 length and empty data
   * which creates an empty string

   * Revision history
   - SLKP 20170517 porting from v3 branch
   */
  virtual void Initialize(bool = true);

  /**
   * Clear the string (delete allocated memory)
   * which creates an empty string

   * Revision history
   - SLKP 20170612 iniital implementation
   */
  virtual void Clear(void);

  /**
   * Construct a string representation of a long interger
   * @param number: the number to convert to a string

   * Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   */
  _String(long const number);

  /**
   * Construct a string long enough to hold the specified # of chars
   * Contents will be initialized to 0
   * @param lengths: the number of chars to store

   * Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   */
  _String(const unsigned long sL);

  /**
   * Construct a string representation of a hyFloat(double) to string,
   * using a format string (default is to use PRINTF_FORMAT_STRING formatting)
   * @param number : The floating number to convert to string
   * @param format : The C-style format string to use for the conversion

   * Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   */
  _String(const hyFloat number, const char *format = nil);

    /**
     * Construct a string representation of a hyFloat(double) to string,
     * using with the required digits of precision ("%.[N]g") specified
     * @param number : The floating number to convert to string
     * @param unsigned char : The number of significant digits
     
     * Revision history
     - SLKP 20181009 initial implementation
     */
    _String(const hyFloat number, unsigned char digits_of_precision);

  /**
   * A RHS copy constructor
   * @param str : the string to copy from

   * Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   */
  _String(const _String &str);

  /**
   * A RHS move constructor
   * @param str : the string to copy from
   
   * Revision history
   - SLKP 20180920 initial implementation
   */
   _String(_String && str);

  /**
   * A RHS move constructor for string buffer
   * @param str : the string to copy from
   
   * Revision history
   - SLKP 20180920 initial implementation
   */
   _String(_StringBuffer && str);

  /**
   * The purpose of this constructor is a "move" contents from a dynamically
   * allocated string to a new string variable; it does so without allocating
   * memory (this is a hack for C++ move semantics)
   * After a call to this dynamic_string will be DELETED, so it CANNOT be used
   * again
   * @param dynamic_string: the source string to move data from

   * Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   */
  _String(_String *dynamic_string, bool dynamic = true);

  /**
   * Copy a part of another string into this string
   *

   * @param str   : The source string
   * @param start : Start of the range to copy
   * @param end   : End of the range to copy
   * @sa NormalizeRange for a discussion on ranges

   * Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   */
  _String(const _String &str, long start, long end);

  /**
   * Create a string with the contents of a C-style (0-terminated)
   * char array (they are copied)

   * @param c_string   : The source C char array
   * Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   */
  _String(const char *c_string); // v3

  /**
   * Create a string with the contents of a C-style (0-terminated)
   * wide-char array (they are copied); only single byte characters
   * are copied

   * @param wc_string   : The source C wchar_t char array
   * Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   */

  _String(const wchar_t *wc_string);
  /**
   * Create a string with the from a single charcater
   * @param c   : The source character
   * Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   */
  _String(const char c);

  /**
   * Create a string with several consecutive copies of the source string
   * @param str    : the source string
   * @param copies : the number of copies
   * Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   */
  _String(const _String &str, unsigned long copies);

  /**
   * Create a string with the contents of an open file
   * the file will be rewound and is assumed to be open for reading

   * @param file    : the source file handle
   * @param read_this_many: if -1, then rewind the file and read all of its
   contents, otherwise read 'read_this_many' characters from current position
   * Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   - SLKP 20170623 added the option to read a specified number of chars
                   from the current position of an open file (to handle fscanf
                   specifically); also added a check that the # of chars read
                   was the same as the one requested.
   */
  _String(FILE *file, long read_this_many = -1L);

  /**
   *  A desctructor which respects reference counts
   *  Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   */
  virtual ~_String(void);

  /**
   * Create a dynamically allocated (shallow) copy of this object
   * @return a shallow copy of this object (for strings, shallow == deep copy)

   *  Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   */
  virtual BaseRef makeDynamic(void) const;

  /** Create a shallow copy of the argument (assumed castable to _String*)
   in this object; this will be cleared out prior to this operation


   @param source: the string to duplicate

   *  Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   [CHANGE-NOTE SLKP, this behavior may not be consistently enforced in old
   code]

   */
  virtual void Duplicate(BaseRefConst source);

  /** Create a shallow copy of the argument

   @param rhs : the right hand side of the assignment

   *  Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   [CHANGE-NOTE SLKP, changed parameter type from _String to _String const&]

   */
  void operator=(_String const &rhs);

  void operator=(_String &&rhs);

  /*
   ==============================================================
   Getters and setters
   ==============================================================
   */

  /**
   * Retrieve a writable element at index x.
   * Internal error results if [] is called on an invalid index

   * @param index : the index (0-based) of a character to retrieve
   * @return      : reference to the character at the specified index
   *  Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   [CHANGE-NOTE SLKP 20170517, used to ignore errored indices]
   */
  virtual char &operator[](long index);

  /**
   * Retrieve a read-only element at index x. If the index is invalid,
   return default_return (\0)

   * @param index : the index (0-based) of a character to retrieve
                    if index < 0, return a character this far from the end;
                    e.g. -1 returns  the last character (for non-empty strings)
                    -2 : the second to the last character (for strings with 2 or
   more chars), etc
   * @return      : the character at the specified index or default_return
   * @sa get_char
   *  Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   [CHANGE-NOTE SLKP 20170517, used to have unsigned long argument]
   - SLKP 20170623 handling negative indices; SEMANTICS CHANGE

   */
  char operator()(long index) const;

  /**
   * Retrieve a read-only element at index x.
   * same as s(i), but with this function you don't have to write (*s)(i) for
   pointers

   * @param index : the index (0-based) of a character to retrieve
   * @return      : the character at the specified index or default_return
   * @sa operator ()
   *  Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   */
  virtual char get_char(long index) const {
    if (index >= 0L && index < s_length) {
      return s_data[index];
    }
    return _String::default_return;
  }

  /**
   * Retrieve a read-only element at index x.
   * WITHOUT ANY RANGE CHECKING
   * @param index : the index (0-based) of a character to retrieve
   * @return      : the character at the specified index or default_return
   * @sa operator ()
   * @sa get_char
   *  Revision history
   - SLKP 20170616 initial implementation
   */
  inline char char_at(unsigned long idx) const { return s_data[idx]; }

  /** The sole purpose of this function is to allow warning-free compilation of
   calls like array [string.getUChar (i)], otherwise you'd get warnings about
   atypical indexing types

   *  Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   */
  inline unsigned char get_uchar(long i) const {
    return (unsigned char)s_data[i];
  }

  /** Get the length of this string
   @return the length of the string
   *  Revision history
      - SLKP 20170517 reviewed while porting from v3 branch
   */

  inline unsigned long length(void) const { return s_length; }

  /** Check if the string is emtpy
   *  Revision history
      - SLKP 20170615 initial implementation
   */

  inline bool empty(void) const { return s_length == 0UL || s_data == nil; }

  /** Check if the string is non-emtpy
   *  Revision history
   - SLKP 20170621 initial implementation
   */

  inline bool nonempty(void) const { return !empty(); }

  /** Store the supplied character in a given index; functionally almost the
   same as str[index] = date, but neater to write than (*str)[index] = data, and
   this also ignores invalid indices

   *  Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   [CHANGE-NOTE SLKP 20170517, used to have 'long' argument]
   */
  void set_char(unsigned long index, char const data);

  /** Retrieve the read-only char * for the string contents
   A convenience function to avoid writing (const char*) (*this)

   @return string data (could be null!, no checks performed)
   @sa operator char *

   *  Revision history
   - SLKP 20170608 reviewed while porting from v3 branch
   */
  const char *get_str(void) const;

  /*
   ==============================================================
   Type conversions
   ==============================================================
   */

  /** Retrieve the read-only char * for the string contents

   @return string data (could be null!, no checks performed)
   @sa get_str

   *  Revision history
   - SLKP 20170608 reviewed while porting from v3 branch
   */
  operator const char *(void)const;

  /**
   * Converts a string of form "[\d\.]\+" into a floating point number
   * via a call to strtod
   * \n\n \b Example: "3.14" becomes 3.14

   *  Revision history
   - SLKP 20170608 reviewed while porting from v3 branch
   */

  hyFloat to_float(void) const;

  /**
   * Converts a string into an integer number
   * via a call to strtol
   * \n\n \b Example: "3.14" becomes 3

   *  Revision history
   - SLKP 20170608 reviewed; was not in v3 branch
   */

  long to_long(void) const;

  /**
   * Obtain a string representation of this string
   * Add a reference counter and return 'this'
     @return this string with an extra reference counter
   *  Revision history
   - SLKP 20170608 reviewed while porting from v3 branch
  */
  virtual BaseRef toStr(unsigned long = 0UL);

  /**
   * Turns seconds into a time string in the form "hh:mm:ss"
   * \n\n \b Example:
   * \code
   * long time_diff = 459132;
   * _String("").FormatTimeString(time_diff);
   * \endcode
   * @param time_diff Seconds of time
   * @return duration string to "127:32:12" in the example.
   *  Revision history
   - SLKP 20170616; reviewed while porting from the v3 branch
   */

  static const _String FormatTimeString(long const);

  /*
   ==============================================================
   Comparisons
   ==============================================================
  */

  /** Perform a lexicographic comparison of two strings
   @param rhs right hand side of the comparison
   @returns less, equal, greater
   *  Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   [CHANGE-NOTE SLKP 20170517,
   return type from char to hyComparisonType
   argument from _String const* to _String const & ]

   */
  hyComparisonType Compare(_String const &rhs) const;

  /** Perform a lexicographic comparison of two strings ignoring case.
   Same as casting both strings to lower case and running Compare

   @param rhs right hand side of the comparison

   @returns less, equal, greater
   *  Revision history
   - SLKP 20170517 initial implementation

   */
  hyComparisonType CompareIgnoringCase(_String const &rhs) const;

  /** Obvious lexicographic comparisons, mostly making calls to Compare
   *  Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   */
  bool operator==(const _String &) const;
  bool operator>(const _String &) const;
  bool operator<(const _String &) const;
  bool operator>=(const _String &) const;
  bool operator<=(const _String &) const;
  bool operator!=(const _String &) const;
  bool Equal(const _String &) const;
  bool EqualIgnoringCase(const _String &) const;
  bool Equal(const char) const;

  /** match this string to a shell style pattern where the wildchar specifies
   "match zero or more of anything"

   @param pattern : the pattern to match
   @param wildchar : the charcter to treat as a wild char
   @param start_this : start matching at this position in "this"
   @param start_pattern : start matching at this position in *pattern*
   @param wildchar_matches: if given, push indices of ranges that matched wildcards

   @return did the string match the pattern

   *  Revision history
   - SLKP 20170517 reviewed while porting from v3 branch
   [CHANGE-NOTE SLKP 20170517 change pattern type to _String const& from _String
   const *]
   - SLKP 20181024 add the optional _SimpleList argument to store the index ranges
                   which matched the wildcards

   */
  bool EqualWithWildChar(_String const &pattern, char const wildchar = '*',
                         unsigned long start_this = 0UL,
                         unsigned long start_pattern = 0UL,
                         _SimpleList * wildchar_matches = nil) const;

  /*
   ==============================================================
   Content-modification and extraction methods
   ==============================================================
   */

  /**
  * String concatenation operator, returns "thisrhs"
  * \n\n \b Example: \code _String new_string = _String("A") & _String("B")
  \endcode
  * @param  rhs : the suffix to concatenate to this
  * @return "AB"
* @sa EscapeAndAppend()

   *  Revision history
      - SLKP 20170519 reviewed while porting from v3 branch
  */
  _String operator&(const _String &rhs) const;

  /**
   * Removes part of string that is between the two specified indices
   * \n\n \b Example: \code _String new_string = _String("AAABBBCCC").Chop(3,5)
   \endcode
   * @param start The starting index to chop from
   * @param end The ending index to chop from
   * @return "AAACCC"
   * @sa Cut()
   * @sa Trim()
   *  Revision history
      - SLKP 20170519 reviewed while porting from v3 branch
   */
   _String Chop(long start, long end) const;

  /**
   * Cuts part of string that is between the two specified indices (0-bases,
   inclusive)
   * \n\n \b Example: \code _String new_string = _String("AAABBBCCC").Cut(3,5)
   \endcode
   * @param start The starting index to cut from
   * @param end The ending index to cut from
   * @return "BBB"
   * @sa Chop()
   * @sa Trim()
   *  Revision history
      - SLKP 20170519 reviewed while porting from v3 branch
   */
   _String Cut(long, long) const;

  /**
   * Delete a range of chars from the string (0-based, inclusive indices)
   * \n\n \b Example: \code _String("AAABBBCCC").Delete(3,5) \endcode
   * @param start The starting index to delete from
   * @param end   The ending index to delete to
   * @return Transforms String to "AAACCC"
   * @sa Chop()
   *  Revision history
      - SLKP 20170519 reviewed while porting from v3 branch
   */
  void Delete(long, long);

  /**
   *
   * In-place reversed string
   * \n s[0]...s[sLength-1] => s[sLength-1]...s[0]
   * \n\n \b Example: \code _String("ABC").Flip() \endcode
   * @return nothing
   * @sa Reverse
   *  Revision history
      - SLKP 20170519 reviewed while porting from v3 branch
   */
  void Flip(void);

  /**
   *
   * Return a reversed string, leaving the original unchanged
   * \n s[0]...s[sLength-1] => s[sLength-1]...s[0]
   * \n\n \b Example: \code _String("ABC").Reverse() \endcode
   * @return "CBA"
   *  Revision history
    - SLKP 20170519 reviewed ; (was missing in v3)
   */
   _String Reverse(void) const;

  /**
   * Insert a char at a given position
   * \n\n \b Example: \code _String("AA").insert('C',0) \endcode
   * @param c Character to insert
   * @param where The position (0-based) to insert the character into,
   values less than 0 append to the string
   * @return "CAA"
   *  Revision history
     - SLKP 20170519 reviewed while porting from v3 branch
  */

  void Insert(char, long);

  /**
   * Trim the string in place to retain characters beween the two indices
   (0-bases, inclusive)
   * \n\n \b Example: \code _String("AAABBBCCC").Trim(3,5) \endcode
   * @param start The starting index to cut from
   * @param end  The ending index to cut from
   * @return Transforms string to "BBB"
   * @sa Cut()
   * @sa Chop()
   *  Revision history
      - SLKP 20170519 reviewed while porting from v3 branch
      [CHANGE-NOTE SLKP 20170519 remove the bool argument for memory handling]
   */

  virtual void Trim(long, long);

  /**
   * Converts string to a particular case
   @param conversion_type: which case ? kStringUpperCase or kStringLowerCase

   *  Revision history
      -SLKP 20170614 reviewed while porting from v3 branch

      [CHANGE-NOTE SLKP 20170614 consolidated LoCase and UpCase;
       changed behavior from in-place to returning a modified string
      ]
   */
  const _String ChangeCase(hy_string_case conversion_type) const;
  void   ChangeCaseInPlace(hy_string_case conversion_type);

  /**
   * Returns a list from a string split by a substr
   * \n\n \b Example: _String("hyphy, gattaca, protease").Tokenize(",") will
   create a list {"hyphy","gattaca","protease"}
   * @param splitter The substring to split the string by
   * @return A point to a *_List that holds a list of the resultant strings.
   Retrieve one by list->lData[i]
   *  Revision history
    -SLKP 20170615 reviewed while porting from v3 branch; previous
   impelementation would not handle empty string splitter;
     ]
   */
  const _List Tokenize(_String const &splitter) const;

    /**
     * Returns a list from a string split by a any of the valid chars
     * @param splitter a look table of characters
     * @return A point to a *_List that holds a list of the resultant strings. Retrieve one by list->lData[i]
     *  Revision history
     -SLKP 20170912 initial impementation
     
     */
    const _List  Tokenize (const bool splitter[256]) const;
    /**
   * Decorates the string with quotes

   * @param quote_char which character to use as a "quote"
   * @return quote_char + *this + quote_char
   *  Revision history
      -SLKP 20170616 reviewed while porting from v2.3 branch
      -
   */
    
    
  const _String Enquote(char quote_char = '\'') const;

  /**
   * Decorates the string with open/close chars

   * @param quote_char which character to use as a "quote"
   * @return open_char + *this + close_char
   *  Revision history
   -SLKP 20170626 initial implementation
   -
   */
  const _String Enquote(char open_char, char close_char) const;

  /**
   * Returns a copy of the string with all spaces removed
   * \n\n \b Example: \code _String("   h  y p    h  y").KillSpaces \endcode
   * @param result The string that will have stripped spaces.
   * @sa CompressSpaces()
   * @return The example would return "hyphy"
   *  Revision history
      -SLKP 20170616 reviewed while porting from v3 branch; changed from in
   place to return by value
   */
  const _String KillSpaces(void) const;

  /**
   * Replaces all runs of white spaces with a single ' ' character
   * \n\n \b Example: \code _String("   h  y p    h  y").CompressSpaces()
   \endcode
   * @return Example would return the string to " h y p h y"
   * @sa KillSpaces()
   *  Revision history
      -SLKP 20170616 reviewed while porting from v3 branch; changed from in
   place to return by value
   */
  const _String CompressSpaces(void) const;

  /*
  ==============================================================
  Search functions
  ==============================================================
  */

  /**
   * Find first occurence of the string between "start" and "end" (inclusive)
   * \n\n \b Example: \code _String ("AABBCC").Find("B")\endcode
   * @param pattern The substring to find
   * @param start The 0-based index to start searching from
   * @param end   The 0-based index to search to (inclusive); -1 : end of string
   * @return Returns the index of the first instance of the pattern, kNotFound
   (<0) if not found. 2 in the example
     @sa FindBackwards
   *  Revision history
   - SLKP 20170608 reviewed while porting from v3 branch
   */
  long Find(const _String &pattern, long start = 0L,
            long end = kStringEnd) const;

  /**
   * Find first occurence of the string between "start" and "end" (inclusive)
   * looking backwards (i.e. last occurrence reported)
   * \n\n \b Example: \code _String ("AABBCC").Find("B")\endcode
   * @param pattern The substring to find
   * @param start The 0-based index to start searching from
   * @param end   The 0-based index to search to (inclusive); -1 : end of string
   * @return Returns the index of the first instance of the pattern, kNotFound
   (<0) if not found. 3 in the example
   @sa Find
   *  Revision history
   - SLKP 20170608 reviewed while porting from v3 branch
   */

  long FindBackwards(const _String &pattern, long start = 0L,
                     long end = kStringEnd) const;
  /**
  * Find first occurence of the character between "start" and "end" (inclusive)
  * Uses a sentinel linear search
  * \n\n \b Example: \code _String ("AABBCC").Find('B')\endcode
  * @param p The character to find
  * @param start The 0-based index to start searching from
  * @param end   The 0-based index to search to (inclusive); -1 : end of string
  * @return Returns the index of the first instance of the pattern, kNotFound
  (<0) if not found. 2 in the example

  *  Revision history
  - SLKP 20170608 reviewed while porting from v3 branch
  */
  long Find(const char p, long start = 0L, long to = kStringEnd) const;

/**
 * Find first occurence of the any of the characters marked in the lookup buffer (0/1) between "start" and "end" (inclusive)
 * Uses a sentinel linear search
 * \n\n \b Example: \code _String ("AABBCC").Find('B')\endcode
 * @param lookup The lookup table whioch marks which characters are value
 * @param start The 0-based index to start searching from
 * @param end   The 0-based index to search to (inclusive); -1 : end of string
 * @return Returns the index of the first instance of the pattern, kNotFound (<0) if not found. 2 in the example
 
 *  Revision history
 - SLKP 20170912 introduced
 */

   long    Find (const bool lookup[256] , long start = 0L, long to = kStringEnd) const ;
   long    FindAnyCase (const bool lookup[256] , long start = 0L, long to = kStringEnd) const ;
/**
  * Find first occurence of the string between "start" and "end" (inclusive)
  * @see Find() for parameter explanation
  *  Revision history
   - SLKP 20170612; reviewed and modifed to be the same as Find with case
  normalization while porting from the v3 branch
  */
    
    

  long FindAnyCase(_String const &pattern, long start = 0L,
                   long to = kStringEnd) const;

  /**
   * Replace string `pattern` with string `replace`, all occurences true/false
   * \n\n \b Example: \code _String("AAABBBCCCBBB").Replace("BBB","ZZ",true)
   \endcode
   * @param pattern The substring to replace
   * @param replace The substring to replace the value with
   * @param flag If true, replace all.
   * @return "AAAZZCCCZZ"

   *  Revision history
     - SLKP 20170614; reviewed while porting from the v3 branch
   */

  const _String Replace(const _String &pattern, const _String& replace,
                        bool replace_all) const;

  /**
   * Locate the first non-space character of the string
   * \n\n \b Example: \code _String ("    hyphy").FirstNonSpaceIndex()\endcode
   * @param start Beginning of string search
   * @param end End of string search
   * @param direction Choose between kStringDirectionForward and
   kStringDirectionBackwards
   * @return The char of the first non-space, in the example, 'h'.
   * @see FirstNonSpaceIndex()

   *  Revision history
   - SLKP 20170614; reviewed while porting from the v3 branch
     [CHANGE-NOTE SLKP 20170614 changed to a call to _FindFirstIndexCondtion]

   */

  char FirstNonSpace(
      long start = 0, long end = kStringEnd,
      hy_string_search_direction direction = kStringDirectionForward) const;

  /**
   * Locate the first non-space character of the string
   * \n\n \b Example: \code _String ("    hyphy").FirstNonSpaceIndex()\endcode
   * @param start Beginning of string search
   * @param end End of string search
   * @param direction Choose between kStringDirectionForward and
   kStringDirectionBackwards
   * @return The index of the first non-space, in the example, 4.
   * @see FirstNonSpaceIndex()

   *  Revision history
   - SLKP 20170614; reviewed while porting from the v3 branch

   */
  long FirstNonSpaceIndex(
      long start = 0, long end = kStringEnd,
      hy_string_search_direction direction = kStringDirectionForward) const;

  /**
   * Locate the first space character of the string
   * \n Returns index of first space character
   * \n\n \b Example: \code _String ("h yphy").FirstSpaceIndex()\endcode
   * @param start starting index
   * @param end ending index to search
   * @param direction Choose between kStringDirectionForward and
   kStringDirectionBackwards
   * @return Returns the index of the first non-space. 1 in the example.
   * @sa FirstSpaceIndex()

   *  Revision history
    - SLKP 20170614; reviewed while porting from the v3 branch
      [CHANGE-NOTE SLKP 20170614 changed to a call to _FindFirstIndexCondtion]
   */
  long FirstSpaceIndex(
      long start = 0, long end = kStringEnd,
      hy_string_search_direction direction = kStringDirectionForward) const;

  /**
   * Locate the first non-space character of the string following one or more
   spaces
   * \n Returns index of first space character
   * \n\n \b Example: \code _String ("h yphy").FirstSpaceIndex()\endcode
   * @param start starting index
   * @param end ending index to search
   * @param direction Choose between kStringDirectionForward and
   kStringDirectionBackwards
   * @return Returns the index of the first non-space. 1 in the example.
   * @sa FirstSpaceIndex()
   *  Revision history
   - SLKP 20170614; reviewed while porting from the v3 branch
   [CHANGE-NOTE SLKP 20170614 seems that the search in reverse direction was not
   implemented correctly]
  */

  long FirstNonSpaceFollowingSpace(
      long start = 0, long end = kStringEnd,
      hy_string_search_direction direction = kStringDirectionForward) const;

  /**
   * Checks to see if String begins with substring
   * \n\n \b Example: \code _String("hyphy").BeginsWith("h")\endcode
   * @param pattern Substring
   * @param case_sensitive If true, it will be case sensitive. Default is case
   sensitive.
   * @param from: start matching *this at this position
   * @return true if string begins with substring. Example returns true
   * @sa EndsWith()
   *  Revision history
   - SLKP 20170615; reviewed while porting from the v3 branch, renamed to camel
   case (not cheap) added the third argument to check for match from a given
   position in this
  */
  
  bool BeginsWith (_String const& pattern, bool case_sensitive = true, unsigned long from = 0UL) const;
  bool BeginsWith (bool const lookup[256], bool case_sensitive = true, unsigned long from = 0UL) const;

  /**
   * Checks to see if String ends with substring
   * \n\n \b Example: \code _String("hyphy").EndsWith("hy")\endcode
   * @param pattern Substring
   * @param case_sensitive If true, it will be case sensitive. Default is case
   sensitive.
   * @return true if string ends with substring. Example returns true
   * @sa BeginsWith()
   *  Revision history
      - SLKP 20170616; reviewed while porting from the v3 branch, renamed to
   camel case (not cheap)
   */
  bool EndsWith(_String const &pattern, bool case_sensitive = true) const;

  /**
   * Checks to see if String starts with substring and it can't be extended to
   make a valid ident
   * by checking the next character only
   * \n\n \b Example: \code
   _String("return;").StarsWithAndIsNotAnIdent("return");
   _String("return_me").StarsWithAndIsNotAnIdent("return")\endcode
   * @param pattern the prefix pattern
   * @return true if string starts with substring and can't be extended to a
   identifier. Example 1 would return true, and example 2 would return false
   *  Revision history
      - SLKP 20170616; reviewed while porting from the v2.3 branch, renamed to
   camel case (not cheap)
   * @sa BeginsWith()
   */
  bool BeginsWithAndIsNotAnIdent(_String const &) const;
  /*
   ==============================================================
   Parser-related functions
   TODO: possible deprecate when the move to the grammar is effected
   ==============================================================
   */

  /**
   * Starting at index [argument 1],
   * find a span that encloses an expression (nested) delimited by char[argument
   2]
   * and char[argument 3] (e.g. {}, ()) respecting quotes (argument 4), and
   allowing
   * escaped characters (argument 5)
   * \n SLKP 20090803
   *
   * @param &from The starting position of the segment will be stored here
   * @param open The first character to look for. For example, and open bracket
   '[' or open paranthesis '('
     Can also be any object that supports char == object checks
   * @param close The first character to look for. For example, and open bracket
   ']' or open paranthesis ')'
     Can also be any object that supports char == object checks
   * @param options: a bitmask of options, if fExtractRespectQuote is mixed in
   then do not look within enquoted parts of the string if set if
   fExtractRespectEscape is mixed in do not consider \char as matches to char
   when searching
   *
   * @return Ending position is returned
   *   kNotFound is returned if the starting character could not be found or the
   expression did not terminate before the end of the string
   *
   *  Revision history
     - SLKP 20170614; reviewed while porting from the v2.3 branch; convered the
   two bool flags to a bit-mask so that the calls can be more explict
     - SLKP 20170615; included support for singly quoted literals
     - SLKP 20171211: added support for generic callbacks to check whether or not the final character has been found
  */

    //=============================================================
  
  
  template <class DELIM> long ExtractEnclosedExpression (long& from, DELIM open, DELIM close, int options) const {
    long   current_position = from,
    current_level    = 0L;
    
    bool       respect_quote = options & fExtractRespectQuote,
               respect_escape = options & fExtractRespectEscape,
               one_level_only = options & fExtractOneLevelOnly,
               do_escape = false;
    
    char       quote_state = '\0',
               this_char = get_char (current_position);
      
    while (this_char) {
      bool       check_quote = false;
        
      if (do_escape) {
        do_escape = false;
      } else {
        // also need to handle cases when quotes are in the open / close set
        
        if ((this_char == '"' || this_char == '\'') && respect_quote && !do_escape) {
          if (quote_state == '\0') {
            check_quote = true;
            quote_state = this_char;
          } else {
            if (this_char == quote_state) {
              check_quote = true;
              quote_state = '\0';
            }
          }
        }
        if (open == this_char && (check_quote || quote_state == '\0')) {
            // handle the case when close and open are the same
          if (current_level == 1L && close == this_char && from < current_position) {
            return current_position;
          }
          if (current_level == 0L) {
            from = current_position;
            current_level++;
          } else {
            if (!one_level_only) {
              current_level++;
            }
          }
          
        } else if (close == this_char && (check_quote || quote_state == '\0')) {
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
      
      this_char = get_char (++current_position);
        
    }
      
    // check if \0 is a valid terminator
      
   if (close == this_char) {
       if (current_level == 1L && from < current_position) {
           return current_position;
       }
   }
    
    return kNotFound;
  }
  
  /**
   * Starting at a 0-based index [argument 1],
   * find a span that terminates in one of the characters in [argument 2], while
   * respecting (), [], {}, "" and escapes
   * \n SLKP 20090805
   * @param start the index to start the search from
   * @param terminator The terminator to find
   * @return kNotFound is returned if the starting character could not be found
   or the expression did not terminate before the end of the string
   * @sa IsALiteralArgument()
   *  Revision history
      - SLKP 20170615   reviewed while porting from the v2.3 branch;
                        for the string; included support for singly quoted
                        literals; cleaned up the logic, and fixed broken logic for terminator > 1
                        char long
   
      - SLKP 20180921   converted into a template to make it possible to search
                        for multiple terminators
   */

  template <typename TERMINATOR> long FindTerminator(long start, TERMINATOR const &terminator) const{
    
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
                current_position++;
                continue;
              case ')':
                if (paren_depth > 0L) {
                  paren_depth --;
                  current_position++;
                  continue;
                }
                break;
              case '[':
                square_depth++;
                current_position++;
                continue;
              case ']':
                if (square_depth > 0L) {
                  square_depth --;
                  current_position++;
                  continue;
                }
                break;
              case '{':
                curly_depth++;
                current_position++;
                continue;
              case '}':
                if (curly_depth > 0L) {
                  curly_depth --;
                  current_position++;
                  continue;
                }
                break;
            }
            
            if (curly_depth == 0L && square_depth == 0L && paren_depth == 0L) {
              if (BeginsWith (terminator, true, current_position)) {
                return current_position;
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

  /**
   * Strips quotes from around the string if present (in place)
   * \n\n \b Example: \code _String("\"hyphy\"").StripQuotes("")\endcode
   * @param open_char : the opening quote char
   * @param close_char : the closing quote char
   * @return : true if the string was enquoted and the quotes had been stripped

   *  Revision history
      - SLKP 20170616   reviewed while porting from the v3 branch
      - SLKP 20170702   return TRUE if successfully stripped quotes

   */
  bool StripQuotes(char open_char = '"', char close_char = '"');

    /**
     * Strips quotes from around the string if present (in place) for multiple delimiters at once
     * \n\n \b Example: \code _String("\"hyphy\"").StripQuotes("\"'","\"'")\endcode
     * @param open_char : the opening quote chars (paired with close_char)
     * @param close_char : the closing quote char (paired with open char)
     * @return : true if the string was enquoted and the quotes had been stripped

     *  Revision history
        - SLKP 20200508  initial

     */
  bool StripQuotes(char const *, char const *);

  /**
   * Checks if String is valid ident
   * \n A valid ident is any alphanumeric or '_'
   * \n\n \b Example: '$hyphy' is not legal.  'hy_phy' is legal.
   * @param options if fIDAllowCompound is set, treat 'x.y.z' as a valid
   identifier, if fIDAllowFirstNumeric is set, consider '2x' a valid identifier
   * @sa ConvertToAnIdent();
   *  Revision history
      - SLKP 20170616   reviewed while porting from the v3 branch
                        changed the argument to bitmask, added
   fIDAllowFirstNumeric

   */
  bool IsValidIdentifier(int options = fIDAllowCompound) const;

  /**
   * Converts a string to a valid ident
   * \n A valid ident is any alphanumeric or '_'
   * \n\n \b Example: \code _String("$hyphy") \endcode
   * @param strict If strict, only alphabetic, no numerals.
   * @param options if fIDAllowCompound is set, treat 'x.y.z' as a valid
   identifier, if fIDAllowFirstNumeric is set, consider '2x' a valid identifier
   * @sa IsValidIdentifier();
   * @return the example would return "_hyphy"

   *  Revision history
   - SLKP 20170616   new implementation based on _IsValidIdentifierAux
                     changed the argument to bitmask, added fIDAllowFirstNumeric
                     changed from in-place modification to returning a modified
   string this function actually respects fIDAllowCompound now
   */
  const _String ConvertToAnIdent(int options = fIDAllowCompound) const;

  /**
   * If it is enclosed in quotes, then it is a literal argument
   * \n \n \b Example: "\"hyphy \"quote\"\"" is a literal argument;
   * @param strip_quotes if set to TRUE and the expression is a literal, trim
   the quotes
   *  Revision history
   - SLKP 20170616   reviewed while porting from the v3 branch
                     added support for single quotes in addition to double
   quotes
   */

  bool IsALiteralArgument(bool strip_quotes = false);

  /**
   * Examine the string argument contained in this object, decide what it is,
   and process accordingly
   * \n\n \bExample: \code 'hyphy'.ProcessVariableReferenceCases (object)
   \endcode is a direct reference to object hyphy
   * \n\n \bExample: \code '\"hy\"+\"phy\"'.ProcessVariableReferenceCases
   (object) \endcode is a direct reference to object hyphy
   * \n\n \bExample: \code '*hyphy'.ProcessVariableReferenceCases (object)
   \endcode is a reference to the object whose name is stored in the string
   variable hyphy
   * \n\n \bExample: \code '**hyphy'.ProcessVariableReferenceCases (object)
   \endcode is a reference to the object whose name is stored in the string
   variable hyphy in the global context
   * @param referenced_object will store the handled variable ID
   * @param context is the namespace of the referenced object; could be nil
   * @return one of HY_STRING_INVALID_REFERENCE    HY_STRING_DIRECT_REFERENCE
   HY_STRING_LOCAL_DEREFERENCE    HY_STRING_GLOBAL_DEREFERENCE
   * @see IsValidIdentifier()
   - SLKP 20170616   reviewed while porting from the v2.3 branch
   */

  hy_reference_type
  ProcessVariableReferenceCases(_String &referenced_object,
                                _String const *context = nil) const;

  /*
  ==============================================================
  METHODS
  ==============================================================
  */
    
  /** a by-character iterator
   
   
   @param  cb : a void (char c, unsigned long index) callback argument
   @param  start_at : start the iteration at this position in the string

       - SLKP 20171008   introduced this function

   */

    template <typename CALLBACK> void Each (CALLBACK cb, unsigned long start_at = 0) const {
        for (unsigned long i = start_at; i<s_length; i++) {
            cb ( s_data[i], i );
        }
    }

/** a by-character matching iterator
     
     
     @param  cb : a void (char c, unsigned long index) callback argument
     @param  start_at : start the iteration at this position in the string

         - SLKP 20171008   introduced this function

     */

      template <typename CALLBACK> long Any (CALLBACK cb, unsigned long start_at = 0) const {
          for (unsigned long i = start_at; i<s_length; i++) {
              if (cb ( s_data[i], i )) return i;
          }
          return kNotFound;
      }

  /**
  * Compute Adler-32 CRC for a string
  * \n\n \b Example: \code _String result = new _String ("Wikipedia"); \endcode
  * \n Implementation shamelessly lifted from
  http://en.wikipedia.org/wiki/Adler-32
  * @return the Adler32 checksum. 300286872 returns in the Example

   *  Revision history
   - SLKP 20170614; reviewed while porting from the v3 branch
  */
  long Adler32(void) const;

  /**
   * Generate a random string on
   * @param len (>0) The desired length of the string
   * @param alphabet Which alphabet do the random charcters come from; in nil,
   then this will be generated from 1-128 ASCII codes
   * @return the random string
   *  Revision history
    - SLKP 20170616; reviewed while porting from the v2.3 branch
   */
  static _String const Random(const unsigned long len,
                              const _String *alphabet = nil);

  /**
   * Computes Lempel-Ziv complexity of the string, i.e. roughly the size of the
   substring table
   * that would have been computed using the LZW algorithm
   * @param rec if provided, will store the indices of substrings mapped to
   unique codes
   * @return string complexity (less compressible == higher complexity)
   * \n Example: 1001111011000010 = 6 because subset the input could be reduced
   to ~6 codes
   * The contents of 'rec' would be 0,1,3,7,11,15, implying that the encoded
   substrings would be [0:0] = 1 [1:1] = 0 [2:3] = 01 [4:7] = 1110 [8:11] = 1100
     [12:15] = 0010
    *  Revision history
    - SLKP 20170616; reviewed while porting from the v2.3 branch, not sure
  */
  unsigned long LempelZivProductionHistory(_SimpleList *rec = nil) const;

  /*
   ==============================================================
   Regular Expression Methods
   ==============================================================
   */
  /**
   * Compile a regular expression represented by a _String object.
   * @param pattern the regular expression to compile
   * @param error_code will receive compilation error codes if any
   * @param case_sensitive controls whether or not the RE is case sensitive
   * @param throw_errors if set, errors will result in thrown excptions (_String const type)
   * @return the resulting (opaque) RE datastructure, or NULL if
             compilation failed

   * @sa FlushRegExp
   * @sa GetRegExpError
   *  Revision history
   - SLKP 20170616; reviewed while porting from the v3 branch
                    maded static member of the class, changed argument 1 to
                    const &
   - SLKP 20180803; added the option for automatic error decoding
   */
  static regex_t *PrepRegExp(_String const &pattern, int &error_code,
                             bool case_sensitive, bool throw_errors = false);

  /**
   * Free a reg_exp datastructure previously returned by PrepRegExp
   * @param re the (opaque) data structure for the regular expression
   *  Revision history
   * @sa PrepRegExp
   * @sa GetRegExpError
   - SLKP 20170616; reviewed while porting from the v3 branch
                    maded static member of the class
   */
  static void FlushRegExp(regex_t *re);

  /**
   * Convert internal regexp code into a string message
   * @param code error code
   * @return the string with the decoded error message
   * @sa PrepRegExp
   * @sa FlushRegExp
   *  Revision history
    - SLKP 20170616; reviewed while porting from the v3 branch
                     maded static member of the class
    */
  static const _String GetRegExpError(int code);

  /**
   * Search this string for the first match to regular expression and
   subexpressions return a list of hits (possibly empty) as pairs of ranges; for
   example "hyphy".RegExpMatch("([^y]+).") -> 0,1,0,0, meaning that the entire
   expression matches to [0:1] and the first subexpression matches to [0:0]
   * @param re the regular expression previously compiled by PrepRegExp
   * @param start start matching the string at this position
   * @return the coordinates of matches for the entire expression (first pair),
   and all subexpressions (left to right); empty if no match

   *  Revision history
          - SLKP 20170616; reviewed while porting from the v3 branch
                           return by value vs writing to argument
          - SLKP 20170623; added the option to search from a given start
   position

   * @sa RegExpAllMatches()
   */

  _SimpleList const RegExpMatch(regex_t const *re,
                                unsigned long start = 0) const;

  /**
   * Search this string for the ALL matches to a regular expression (ignoring
   subexpressions) return a list of hits (possibly empty) as pairs of ranges;
   for example "hyphy".RegExpMatch("([^y]+).") -> 0,1,2,4, meaning that [0:1]
   (hy) and [2:4] (phy) match the pattern
   * @param re the regular expression previously compiled by PrepRegExp
   * @return the coordinates of all matches for the entire expression left to
   right; empty if no match

   *  Revision history
   - SLKP 20170616; reviewed while porting from the v3 branch
        return by value vs writing to argument

   * @sa RegExpMatch
    */

  _SimpleList const RegExpAllMatches(regex_t const *re) const;

  /**
     Convenience wrappers for RegExpMatch and RegExpAllMatches taking in regex_t
    arguments where the regular expression is compiled and disposed of
    internally
     @param pattern the regular expression to match
     @param case_sensitive whether to compile the RE as case sensitive or not
     @param handle_errors if set, call application wide error handlers on
    errors, otherwise ignore errors and treat them as a missing match

    * @sa RegExpMatch
    * @sa RegExpAllMatches
    * Revision history
       - SLKP 20170616;  initial implementation
    *
  */
  _SimpleList const RegExpMatch(_String const &pattern, bool case_sensitive,
                                bool handle_errors) const;
  _SimpleList const RegExpAllMatches(_String const &pattern,
                                     bool case_sensitive,
                                     bool handle_errors) const;
    /** given coordinates start and end, converts then to valid string indices
     if called on an empty string, returns 0 and does not change start and end
     if start < 0 it is reset to 0
     if end < 0 or >= string length it is reset to (string length) - 1

     @param start: start of the range (0-based)
     @param end  : end of the range
     @return     : the length of the range

     * Revision history
     - SLKP 20170517 porting from v3 branch
     */
    long NormalizeRange(long &start, long &end) const;


private:
  /** Find the length of the maximum prefix that forms a valid ID

   @param allow_compounds : treat '.' as a valid identifier character (e.g.
   x.y.z)
   @param allow_first_numeric : allow idents that start with a digit (e.g. 2x)
   @param wildcard : treat this character as a valid identifier character (e.g.
   this1.?.x)

   @return the 0-based index of the end of the valid ID prefix (-1 if the prefix
   is empty)

   * Revision history
   - SLKP 20170616 reviewed while porting from the v3 branch
   */

  long _IsValidIdentifierAux(bool allow_compounds, bool allow_first_numeric,
                             char wildcard = '\0') const;

  /** Find the first character in a range that meets a particular condition

     @param start : start of the range to search (0-based)
     @param end : end of the range to search (0-based)
     @direction : forwards or backwards search
     @comparison_function: a function that takes a single argument (char) and
    returns true if it "passes"

    * Revision history
    - SLKP 20170614 factored out common functions for conditional index finding
  */

  template <class CF>
  long _FindFirstIndexCondtion(long start, long end,
                               hy_string_search_direction direction,
                               CF comparison_function) const {
    long requested_range = NormalizeRange(start, end);

    if (requested_range > 0L) {
      if (direction == kStringDirectionForward) {
        for (; start <= end; start++) {
          if (comparison_function(s_data[start])) {
            return start;
          }
        }
      } else {
        for (; end >= start; end--) {
          if (comparison_function(s_data[end])) {
            return end;
          }
        }
      }
    }

    return kNotFound;
  }


  /** this is a utility function which allocates length+1 chars for s_data,
  copies the data from source_string, and sets the terminating 0

  * Revision history
  - SLKP 20170517 factoring repeated functionality

  */
  inline void AllocateAndCopyString(const char *source_string,
                                    unsigned long length);


  /** Factored out core of RegExpMatch and RegExpAllMatches
   * Revision history
    - SLKP 20170616; initial implementation
   */
  const _SimpleList _IntRegExpMatch(const _String &pattern, bool case_sensitive,
                                    bool handle_errors, bool match_all) const;
};

// _______________________________________________________________________

void SetStatusBarValue(long, hyFloat, hyFloat);
void SetStatusLine(_String);
void SetStatusLine(_String, _String, _String, long l);
void SetStatusLine(_String, _String, _String);
void SetStatusLine(_String, _String, _String, long, char);

void SetStatusLineUser(_String const);

void StringToConsole(_String const &, void *extra = nil);
void BufferToConsole(const char *, void *extra = nil);
void NLToConsole(void *extra = nil);
void ObjectToConsole(BaseRef, void *extra = nil);

_String *StringFromConsole(void);

#endif
