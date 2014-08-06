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

#ifndef _HY_STRINGS_
#define _HY_STRINGS_
//#pragma once

#define HY_2014_REWRITE_MASK

#include "baseobj.h"
#include "wchar.h"

#include "regex.h"

class _ExecutionList;     // forward declaration
template<typename PAYLOAD>
class _hyListNumeric; 

template<typename PAYLOAD>
class _hyListReference;

typedef _hyListNumeric <long> _SimpleList;

class _String;

typedef _hyListReference <_String> _List;



#define HY_STRING_INVALID_REFERENCE 0x00
#define HY_STRING_DIRECT_REFERENCE 0x01
#define HY_STRING_LOCAL_DEREFERENCE 0x02
#define HY_STRING_GLOBAL_DEREFERENCE 0x03

#define HY_STRING_MOD_ADLER 65521
#define HY_STRING_DIRECTION_FORWARD   0x01
#define HY_STRING_DIRECTION_BACKWARD  0x02



class _String : public BaseObj {

public:

  /**
  * A constructor.
  * A simple constructor that does nothing
  */
  _String(void);

  /**
  * A constructor that converts a long to string.
  * @param sL The number to convert to string
  */
  _String(long);

  /**
   * A constructor that creates a string ready to hold the specified number of chars.
   * @param sL The number of chars to hold
   */
  _String(const unsigned long sL);
  

  /**
  * A constructor that converts a _Parameter(double) to string.
  * @param sL         The floating number to convert to string
  * @param format     The C-style format string to use for the conversion
  */
  _String(_Parameter, const char * = nil);

  /**
  * A constructor that duplicates from another string.
  */
  _String(const _String &);

  /**
  * A constructor that duplicates from another string.
  */
  _String(_String *);

  /**
  * A constructor that duplicates from another string range.
  * \n Cut a range of the original string
  * @param from The starting position of the string to copy
  * @param to The ending position of the string to copy
  */
  _String(const _String &, long, long);

  /**
  * A constructor that copies from a standard string.
  */
  _String(const char *);

 /**
  * A constructor that copies from a standard wide character string (only single byte characters are actually stored).
  */
  _String(const wchar_t *);

 /**
  * A constructor that copies from a single char.
  */
  _String(const char);

  /**
  * A constructor that reads from an open file.
  * @param F open file stream buffer to copy from
  */
  _String(FILE *);

  /**
  * A destructor.
  */

  virtual ~_String(void);

  /**
  * Returns a dynamic string of the current instance.
  * \n Usage: stringInstance.makeDynamic();
  * @return BaseRef
  */
  virtual BaseRef makeDynamic(void) const;

  /**
  * Initializes _String object to 0 length and 0 s_data
  */
  virtual void Initialize(bool = false);

  /**
  * Duplicates a string
  * \n Usage: string.Duplicate(&existing_string)
  * @param ref A pointer to the string to be duplicated
  * @sa DuplicateErasing()
  * @sa CopyDynamicString()
  */
  virtual void Duplicate(BaseRefConst);

  /**
  * Erases old data and duplicates a string
  * \n Usage: string.DuplicateErasing(&existing_string)
  * @param ref A pointer to the string to be duplicated
  * @sa Duplicate()
  * @sa CopyDynamicString()
  */
  virtual void DuplicateErasing(BaseRefConst);

  /**
  * Returns the character at a specified index
  * \n Usage: char c = string.getChar(5);
  * @param index The int location of the char in the string
  * @return char located at specified location
  */
  const char getChar(long) const;

  /**
  * Sets the character of the string instance at a specified index
  * Avoids (*string)[5] = 'e', and ignores assignments with invalid indices
  * \n Usage: string->setChar(5, 'e');
  * @param index The int location of the char to be replaced
  * @param c The character to set the location with
  * @return Nothing. Changes string instance
  */
  void setChar(unsigned long, char);

  /**
  * Copies a string dynamically and deletes the original string
  * @param s A pointer to the string to be duplicated
  * @param flushMe If true, free existing data before copying
  * @sa Duplicate()
  * @sa DuplicateErasing()
  */
  void CopyDynamicString(_String *, bool = true);


  
  /**
  * Element location function
  * @see getChar()
  */
  virtual char &operator[](long);

  /**
  * Element location function
  * @see getChar()
  */
  char operator()(long) const;

  /**
  * Sets string
  * \n\n \b Example: \code _String str = _String("hyphy"); \endcode
  * @see Duplicate()
  */
  void operator=(const _String &);

  /**
  * Returns the length of the string
  * \n\n \b Example: \code long l = string.Length(); \endcode
  * @return Length of string
  */
  unsigned long Length(void) const;

  /**
  * Append operator
  * \n\n \b Example: \code _String new_string = _String("A") & _String("B")
  * \endcode
  * @return "AB"
  * @sa EscapeAndAppend()
  */
  const _String operator&(const _String& ) const;



  /**
  * Checks if the string contains the substring
  * \n\n \b Example: \code bool contains =
  * _String("AABBCC").ContainsSubstring("BB")); \endcode
  * @param s The substring to check
  * @return Returns true if string contains substring
  */
  bool ContainsSubstring(const _String &) const;

  /**
  * Converts a good ole char*
  * \n\n \b Example: \code string.toStr(); \endcode
  */
  virtual BaseRef toStr(void) const;

  /**
  * Return good ole char*
  */
  virtual operator const char *(void) const;

  /**
  * Returns a good ole char*
  * \n\n \b Example: \code char * new_str = string.getStr(); \endcode
  * @return Returns a good ole char*
  */
  const char *getStr(void) const;

  /**
  * Removes part of string that is between the two specified indices
  * \n\n \b Example: \code _String new_string = _String("AAABBBCCC").Chop(3,5)
  * \endcode
  * @param from The starting index to chop from
  * @param to The ending index to chop from
  * @return "AAACCC"
  * @sa Cut()
  * @sa Trim()
  */
  const _String Chop(long, long) const;

  /**
  * Cuts part of string that is between the two specified indices
  * \n\n \b Example: \code _String new_string = _String("AAABBBCCC").Cut(3,5)
  * \endcode
  * @param from The starting index to cut from
  * @param to The ending index to cut from
  * @return "BBB"
  * @sa Chop()
  * @sa Trim()
  */
  const _String Cut(long, long) const;

  /**
  *
  * In-place flip string
  * \n s[0]...s[s_length-1] => s[s_length-1]...s[0]
  * \n\n \b Example: \code _String("ABC").Flip() \endcode
  * @return Transforms string to "CBA"
  */
  void Flip(void);

  /**
  * Trim the string between from and to
  * \n Cut string from, to (-1 for any means from beginning/to end)
  * \n\n \b Example: \code _String("AAABBBCCC").Trim(3,5) \endcode
  * @param from The starting index to cut from
  * @param to The ending index to cut from
  * @return Transforms string to "BBB"
  * @sa Cut()
  * @sa Chop()
  */

  void Trim(long, long);

  /**
  * Insert a char at a given position (-1 - append)
  * \n\n \b Example: \code _String("AA").insert('C',0) \endcode
  * @param c Character to insert
  * @param pos The position to insert the character into
  * @return "CAA"
  */

  void Insert(char, long);

  /**
  * Delete a range of chars from the string
  * \n\n \b Example: \code _String("AAABBBCCC").Delete(3,5) \endcode
  * @param from The starting index to delete from
  * @param to The ending index to delete to
  * @return Transforms String to "AAACCC"
  * @sa Chop()
  */
  void Delete(long, long);

  /**
  * Replace string 1 with string 2, all occurences true/false
  * \n\n \b Example: \code _String("AAABBBCCC").Replace("BBB","ZZZ") \endcode
  * @param s The substring to replace
  * @param d The substring to replace the value with
  * @param flag If true, replace all.
  * @return "AAAZZZCCC"
  */

  const _String Replace(const _String&, const _String&, bool replace_all) const;

  /**
  * Locate the first non-space character of the string
  * \n\n \b Example: \code _String ("    hyphy").FirstNonSpaceIndex()\endcode
  * @param start Beginning of string search
  * @param end End of string search
  * @param direction Choose between backwards and forwards
  * @return The char of the first non-space, in the example, 'h'.
  * @see FirstNonSpaceIndex()
  */
  char FirstNonSpace(long start = 0, long end = HY_NOT_FOUND, unsigned char direction = HY_STRING_DIRECTION_FORWARD) const;

  /**
  * Locate the first non-space character of the string
  * \n\n \b Example: \code _String ("    hyphy").FirstNonSpaceIndex()\endcode
  * @param start Beginning of string search
  * @param end End of string search
  * @param direction Choose between backwards and forwards
  * @return The index of the first non-space, in the example, 4.
  * @see FirstNonSpaceIndex()
  */
  long FirstNonSpaceIndex(long start = 0, long end = HY_NOT_FOUND, unsigned char direction = HY_STRING_DIRECTION_FORWARD) const;

  /**
  * Locate the first space character of the string
  * \n Returns index of first space character
  * \n\n \b Example: \code _String ("h yphy").FirstSpaceIndex()\endcode
  * @param start starting index
  * @param end ending index to search
  * @param direction If the direction is less than 0, then go backwards
  * @return Returns the index of the first non-space. 1 in the example.
  * @sa FirstSpaceIndex()
  */
  long FirstSpaceIndex(long start = 0, long end = HY_NOT_FOUND, unsigned char direction = HY_STRING_DIRECTION_FORWARD) const;


  /**
  * Find first occurence of the string between from and to
  * \n\n \b Example: \code _String ("AABBCC").Find("B")\endcode
  * @param s The substring to find
  * @param from The index to start searching from
  * @param to The index to search to
  * @return Returns the index of the first instance of the substr, -1 if not
  * found. 2 in the example
  */
  long Find(const _String& s, long from = 0, long to = HY_NOT_FOUND) const;

  /**
  *  @see Find()
  */
  long Find(const char s, long from = 0, long to = HY_NOT_FOUND) const;

  /**
  * Case insensitive Find
  * @see Find()
  */

  long FindAnyCase(const _String&, long from = 0, long to = HY_NOT_FOUND) const;

  /**
  * Find a string
  * @see Find()
  */
  long FindBackwards(const _String&, long from = 0, long to = HY_NOT_FOUND) const;

  /**
  * Compute Adler-32 CRC for a string
  * \n\n \b Example: \code _String result = new _String ("Wikipedia"); \endcode
  * \n Implementation shamelessly lifted from
  * http://en.wikipedia.org/wiki/Adler-32
  * @return the Adler32 checksum. 300286872 returns in the Example
  */
  long Adler32(void) const;

  /**
  * Turns seconds into a time string in the form "hh:mm:ss"
  * \n\n \b Example:
  * \code
  * long time_diff = 459132;
  * _String("").FormatTimeString(time_diff);
  * \endcode
  * @param time_diff Seconds of time
  * @return returns "127:32:12" for the example.
  */

  static const _String FormatTimeString(long time_difference);


  /**
  * Lexicographic comparison
  * \n Checks if Strings are equal lexicographic
  * @param s Second string to compare
  * @return true if strings are equal
  * @sa Compare()
  */
  bool Equal(const _String *) const;

  /**
  * Lexicographic comparison
  * \n Checks if a string is equal to one character
  * @param c the character to compare to
  * @return true if the string is equal to the character
  * @sa Compare()
  */
  bool Equal(const char c) const;

  /**
  * TODO: Lexicographic comparison
  * \n Checks if Strings are equal lexicographic
  * @param s Second string to compare
  * @return 1 if strings are equal, -1 if strings are not
  * @sa Equal()
  */
  char Compare(_String const *) const;

  /**
  * Lexicographic comparison with a wild character
  * \n Checks if Strings are equal lexicographically
  * @param s Second string to compare
  * @param wildchar The wildcharacter
  * @return true if strings are equal
  * @sa Equal()
  */
  bool EqualWithWildChar(const _String&, const char wildChar = '*') const ;

  /**
  * Checks if String is lexicographically greater
  * @see Greater()
  */
  bool operator>(const _String&) const;

  /**
  * Checks if String is lexicographically less
  * @see Less()
  */
  bool operator<(const _String&) const;
  
  /**
  * Checks if String is lexicographically greater or equal
  * @see Greater()
  */
  bool operator>=(const _String&) const;

  /**
  * Checks if String is lexicographically less or equal
  * @see Less()
  */
  bool operator<=(const _String&) const ;

  /**
  * Checks if strings are not lexicographically equal
  * @see Equal()
  */
  bool operator!=(const _String&) const;

  /**
  * Checks if string are lexicographically equal
  * @see Equal()
  */
  bool operator==(const _String&) const;

  /**
  * Checks to see if string contains substring
  * \n\n \b Example: \code _String("hyphy").contains("h")\endcode
  * @return Returns true if string contains substring. Example returns true
  * @see Find()
  */
  bool contains(const _String&) const;

  /**
  * Checks to see if string contains character
  * @return Returns true if string contains character
  * @see contains()
  */
  bool contains(char) const;

  /**
  * Checks to see if String begins with substring
  * \n\n \b Example: \code _String("hyphy").beginswith("h")\endcode
  * @param s Substring
  * @param caseSensitive If true, it will be case sensitive. Default is case
  * sensitive.
  * @return true if string begins with substring. Example returns true
  * @sa contains()
  * @sa startswith()
  * @sa endswith()
  */
  bool startswith(const _String& , bool = true) const;


  /**
  * Checks to see if String ends with substring
  * \n\n \b Example: \code _String("hyphy").endswith("y")\endcode
  * @param s Substring
  * @param caseSensitive If true, it will be case sensitive. Default is case
  * sensitive.
  * @return true if string ends with substring. Example would return true.
  * @sa contains()
  * @sa beginswith()
  * @sa startswith()
  */
  bool endswith(const _String& , bool = true) const;

  /**
  * Converts string to upper case
  * @sa LoCase()
  */
  const _String UpCase(void) const;

  /**
  * Converts string to lower case
  * @sa UpCase()
  */
  const _String LoCase(void) const;



  /**
  * Returns a list from a string split by a substr
  * \n\n \b Example: _String("hyphy, gattaca, protease").Tokenize(",") will
  * create a list {"hyphy","gattaca","protease"}
  * @param s The substring to split the string by
  * @return The list that holds the resultant strings.
  */
  const _List Tokenize(const _String&) const;


  /**
  * Returns a list from a string split by a substr.
  * \n\n \b Example: \code _String("\"hyphy\"").StripQuotes("")\endcode
  * @return string with no quotes. "hyphy" in this example.
  */
  void StripQuotes(char open_char = '"', char close_char = '"');

  /**
  * Checks if String is valid ident
  * \n A valid ident is any alphanumeric or '_'
  * \n\n \b Example: '$hyphy' is not legal.  'hy_phy' is legal.
  * @param strict If strict, only alphabetic, no numerals.
  * @sa ConvertToAnIdent();
  */
  bool IsValidIdentifier(bool = true) const;


  /**
  * Converts a string to a valid ident
  * \n A valid ident is any alphanumeric or '_'
  * \n\n \b Example: \code _String("$hyphy") \endcode
  * @param strict If strict, the first character cannot be a number.
  * @sa IsValidIdentifier();
  * @return the example would return "_hyphy"
  */
  const _String ConvertToAnIdent(bool = true) const;

  /**
  * Removes all spaces in a string
  * \n\n \b Example: \code _String("   h  y p    h  y").KillSpaces \endcode
  * @sa CompressSpaces()
  * @return The string without spaces (e.g. "hyphy" in the example)
  */
   const _String KillSpaces(void) const;

  /**
  * Compresses all runs of spaces in a string into a single space
  * \n\n \b Example: \code _String("   h  y p    h  y").CompressSpaces()
  * \endcode
  * @return "h y p h y"
  * @sa KillSpaces()
  */
  const _String CompressSpaces(void) const;

  /**
  * Shorten the var id by removing the matching beginning portion of the passed
  * in String separated by a "."
  * \n\n \b Example: _String("house.room")._String("house")
  * @param containerID
  * @return _String("house.room")._String("house.") will return "room". However,
  * _String("houseroom") will return "houseroom"
  */
  const _String ShortenVarID(_String const &) const;


  static unsigned long storageIncrement;

  /**
  * A regular expression match
  * @param pattern the compiled regex pattern
  * @return A list that holds the start and end indexes of matches
  * @sa RegExpMatchAll()
  * @sa RegExpMatchOnce()
  */
  const _SimpleList RegExpMatch(regex_t const*) const;

  /**
  * A regular expression match
  * @param pattern the compiled regex pattern
  * @return A list that holds the start and end indexes of matches
  * @sa RegExpMatch()
  * @sa RegExpMatchOnce()
  */
  const _SimpleList RegExpMatchAll(regex_t const *) const;

  /**
  * A regular expression match that only matches once
  * @param pattern A string that holds the regex pattern
  * @preturn A list that holds the start and end indexes of matches
  * @sa RegExpMatch()
  * @sa RegExpMatchAll()
  */
  const _SimpleList RegExpMatchOnce(const _String&, bool, bool) const;

  /**
  * Lexicographically sorts characters string
  * @param index if supplied, then the the i-th element of the list will store the \n
  * original index (in the unsorted string) of the i-th character in the sorted string
  * @return sorted string
  */
  const _String Sort(_SimpleList * = nil) const;

  /**
   * Generate a random string on
   * @param len (>0) The desired length of the string
   * @param alphabet Which alphabet do the random charcters come from; in nil,
   * then this will be generated from 1-128 ASCII codes
   * @return the random string
   */
  static const _String Random(const unsigned long len, const _String *alphabet = nil);

  /**
  * Computes Lempel-Ziv complexity of the string.
  * \n The Lempel-Ziv complexity computes the number of separate substrings in a
  * given string
  * \n Example: 1001111011000010 = 6 because subset of all strings = {1, 0, 01,
  * 1110, 1100, 0010 }
  */
  long LempelZivProductionHistory(_SimpleList * = nil);

  /**
  * Converts a string of form ":[\d\.]\+" into a double
  * \n SLKP 20100831: a utility function to handle the
  * conversion of branch length strings to parameters
  * \n\n \b Example: string = ":3.14" returns 3.14
  * \n If it is not of correct form, it will return 1e-10
  * @sa toNum()
  */
  _Parameter ProcessTreeBranchLength(_Parameter min_value = 1e-10) const;

  /**
  * Converts a string of form "[\d\.]\+" into a double
  * \n\n \b Example: "3.14" becomes 3.14
  * @sa ProcessTreeBranchLength()
  */

  _Parameter toNum(void) const;

  //TODO: Make these protected
  unsigned long s_length;
  char * s_data;


  // Data Fields
  protected:

  
  private:
  
    long NormalizeRange (long & from, long & to) const;
    
    // returns the length of the prefix of the string which constitutes a valid ID
    long _IsValidIdentifierAux (bool = true, char = '\0') const;

};

// _______________________________________________________________________



regex_t* PrepRegExp(const _String& pattern, int & error, bool respect_case);
void     FlushRegExp(regex_t*);
const    _String GetRegExpError(int error);


#endif
