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
#include "wchar.h"
#include "regex.h"
#include "hy_types.h"

#define kStringInvalidReference                 0x00
#define kStringDirectReference                  0x01
#define kStringLocalDeference                   0x02
#define kStringGlobalDeference                  0x03


#define kAppendAnAssignmentToBufferFree         0x01
#define kAppendAnAssignmentToBufferQuote        0x02
#define kAppendAnAssignmentToBufferAssignment   0x04
#define kAppendAnAssignmentToBufferGlobal       0x08

class _SimpleList;
class _List;
class _ExecutionList;


class _String:public BaseObj {

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
    _String (void); // v3;

    /**
     * Standard initalization to 0 length and empty data
     * which creates an empty string
     
     * Revision history
     - SLKP 20170517 porting from v3 branch
     */
    virtual void Initialize(bool = false); // v3; Last reviewed SLKP 20170517

 
    /**
     * Construct a string representation of a long interger
     * @param number: the number to convert to a string
     
     * Revision history
     - SLKP 20170517 reviewed while porting from v3 branch
     */
    _String (long const number);

    /**
     * Construct a string long enough to hold the specified # of chars
     * Contents will be initialized to 0
     * @param lengths: the number of chars to store
     
     * Revision history
     - SLKP 20170517 reviewed while porting from v3 branch
     */
    _String(const unsigned long sL);

    /**
     * Construct a string representation of a hy_float(double) to string,
     * using a format string (default is to use PRINTF_FORMAT_STRING formatting)
     * @param number : The floating number to convert to string
     * @param format : The C-style format string to use for the conversion
     
     * Revision history
     - SLKP 20170517 reviewed while porting from v3 branch
     */
    _String (const hy_float number , const char * format = nil);

    /**
     * A RHS copy constructor
     * @param str : the string to copy from
     
     * Revision history
     - SLKP 20170517 reviewed while porting from v3 branch
     */
     _String (const _String& str);

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
    _String (_String * dynamic_string);

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
    _String (const _String& str, long start, long end);

    /**
     * Create a string with the contents of a C-style (0-terminated)
     * char array (they are copied)
     
     * @param c_string   : The source C char array
     * Revision history
     - SLKP 20170517 reviewed while porting from v3 branch
     */
     _String (const char* c_string); // v3

    /**
     * Create a string with the contents of a C-style (0-terminated)
     * wide-char array (they are copied); only single byte characters
     * are copied
     
     * @param wc_string   : The source C wchar_t char array
     * Revision history
     - SLKP 20170517 reviewed while porting from v3 branch
     */
    
    _String (const wchar_t * wc_string);
    /**
     * Create a string with the from a single charcater
     * @param c   : The source character
     * Revision history
     - SLKP 20170517 reviewed while porting from v3 branch
     */
    _String (const char c);

    /**
     * Create a string with several consecutive copies of the source string
     * @param str    : the source string
     * @param copies : the number of copies
     * Revision history
     - SLKP 20170517 reviewed while porting from v3 branch
     */
    _String (const _String& str, unsigned long copies);

    /**
     * Create a string with the contents of an open file
     * the file will be rewound and is assumed to be open for reading
     
     * @param file    : the source file handle
     * Revision history
     - SLKP 20170517 reviewed while porting from v3 branch
     */
    _String (FILE* file);


    /**
     *  A desctructor which respects reference counts
     *  Revision history
     - SLKP 20170517 reviewed while porting from v3 branch
     */
    virtual     ~_String(void);
    
    /**
     * Create a dynamically allocated (shallow) copy of this object
     * @return a shallow copy of this object (for strings, shallow == deep copy)
     
     *  Revision history
     - SLKP 20170517 reviewed while porting from v3 branch
     */
    virtual     BaseRef makeDynamic (void) const;


    /** Create a shallow copy of the argument (assumed castable to _String*)
     in this object; this will be cleared out prior to this operation
     
     
     @param source: the string to duplicate
     
     *  Revision history
     - SLKP 20170517 reviewed while porting from v3 branch
     [CHANGE-NOTE SLKP, this behavior may not be consistently enforced in old code]
     
     */
    virtual     void    Duplicate (BaseRefConst source);
    
    /** Create a shallow copy of the argument
     
     @param rhs : the right hand side of the assignment
     
     *  Revision history
     - SLKP 20170517 reviewed while porting from v3 branch
     [CHANGE-NOTE SLKP, changed parameter type from _String to _String const&]
     
     */
    void operator = (_String const & rhs);


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
    virtual     char& operator [] (long index);

    
    /**
     * Retrieve a read-only element at index x. If the index is invalid,
     return default_return (\0)
     
     * @param index : the index (0-based) of a character to retrieve
     * @return      : the character at the specified index or default_return
     * @sa get_char
     *  Revision history
     - SLKP 20170517 reviewed while porting from v3 branch
     [CHANGE-NOTE SLKP 20170517, used to have unsigned long argument]
     
     
     */
    char operator  () (long index);
    
    
     const       char    get_char          (long) const;
        /**
         * Retrieve a read-only element at index x. 
         * same as s(i), but with this function you don't have to write (*s)(i) for pointers
         
         * @param index : the index (0-based) of a character to retrieve
         * @return      : the character at the specified index or default_return
         * @sa operator ()
         *  Revision history
           - SLKP 20170517 reviewed while porting from v3 branch
         */

    /** The sole purpose of this function is to allow warning-free compilation of
     calls like array [string.getUChar (i)], otherwise you'd get warnings about
     atypical indexing types
     
     *  Revision history
     - SLKP 20170517 reviewed while porting from v3 branch
     */
     inline const       unsigned char    get_uchar          (long i) const {
      return (unsigned char)s_data[i];
    }
    
    /** Get the length of this string 
     
     @return the length of the string

     *  Revision history
        - SLKP 20170517 reviewed while porting from v3 branch
     */
    
    inline unsigned long length (void) const {
        return s_length;
    }


    /** Store the supplied character in a given index; functionally almost the same as
     str[index] = date, but neater to write than (*str)[index] = data, and this also
     ignores invalid indices
     
     *  Revision history
     - SLKP 20170517 reviewed while porting from v3 branch
     [CHANGE-NOTE SLKP 20170517, used to have 'long' argument]
     */
    void    set_char          (unsigned long index , char const data );

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
     return type from char to _HY_COMPARISON_TYPE
     argument from _String const* to _String const & ]
     
     */
    _HY_COMPARISON_TYPE Compare (_String const& rhs) const;

    /** Perform a lexicographic comparison of two strings ignoring case.
     Same as casting both strings to lower case and running Compare
     
     @param rhs right hand side of the comparison
     
     @returns less, equal, greater
     *  Revision history
     - SLKP 20170517 initial implementation
     
     */
    _HY_COMPARISON_TYPE CompareIgnoringCase (_String const& rhs) const;

    /** Obvious lexicographic comparisons, mostly making calls to Compare
     *  Revision history
     - SLKP 20170517 reviewed while porting from v3 branch
     */
    bool operator==(const _String&) const;
    bool operator> (const _String&) const;
    bool operator< (const _String&) const;
    bool operator>=(const _String&) const;
    bool operator<=(const _String&) const ;
    bool operator!=(const _String&) const;
    bool Equal     (const _String&) const;
    bool EqualIgnoringCase     (const _String&) const;
    bool Equal     (const char) const;

    bool EqualWithWildChar (_String const& pattern, char const wildchar = '*', unsigned long start_this = 0UL, unsigned long start_pattern = 0UL) const;
    
    /** match this string to a shell style pattern where the wildchar specifies "match zero or more of anything"
     
        @param pattern : the pattern to match
        @param wildchar : the charcter to treat as a wild char
        @param start_this : start matching at this position in "this"
        @param start_pattern : start matching at this position in *pattern*
    
        @return did the string match the pattern

         *  Revision history
            - SLKP 20170517 reviewed while porting from v3 branch
            [CHANGE-NOTE SLKP 20170517 change pattern type to _String const& from _String const *]

     */
    
    /*
     ==============================================================
     Content-modification and extraction methods
     ==============================================================
     */
 
    /**
    * String concatenation operator, returns "thisrhs"
    * \n\n \b Example: \code _String new_string = _String("A") & _String("B") \endcode
    * @param  rhs : the suffix to concatenate to this
    * @return "AB"
    * @sa EscapeAndAppend()

     *  Revision history
        - SLKP 20170519 reviewed while porting from v3 branch
    */
    const _String operator & (const _String & rhs) const;

    /**
     * Removes part of string that is between the two specified indices
     * \n\n \b Example: \code _String new_string = _String("AAABBBCCC").Chop(3,5) \endcode
     * @param start The starting index to chop from
     * @param end The ending index to chop from
     * @return "AAACCC"
     * @sa Cut()
     * @sa Trim()
     *  Revision history
        - SLKP 20170519 reviewed while porting from v3 branch
     */
    const _String Chop(long start, long end) const;

    
    /**
     * Cuts part of string that is between the two specified indices (0-bases, inclusive)
     * \n\n \b Example: \code _String new_string = _String("AAABBBCCC").Cut(3,5) \endcode
     * @param start The starting index to cut from
     * @param end The ending index to cut from
     * @return "BBB"
     * @sa Chop()
     * @sa Trim()
     *  Revision history
        - SLKP 20170519 reviewed while porting from v3 branch
     */
    const _String Cut (long, long) const;

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
    void    Delete (long, long);

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
    void    Flip(void);

    /**
     *
     * Return a reversed string, leaving the original unchanged
     * \n s[0]...s[sLength-1] => s[sLength-1]...s[0]
     * \n\n \b Example: \code _String("ABC").Reverse() \endcode
     * @return "CBA"
     *  Revision history
      - SLKP 20170519 reviewed ; (was missing in v3)
     */
    const _String    Reverse(void) const;
    
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
    
    void    Insert (char, long);

    /**
     * Trim the string in place to retain characters beween the two indices (0-bases, inclusive)
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
    
    void    Trim(long, long);

    /*
     ==============================================================
     Search functions
     ==============================================================
     */

    /**
    * Checks if the string contains the substring
    * \n\n \b Example: \code bool contains = _String("AABBCC").ContainsSubstring("BB")); \endcode
    * @param s The substring to check
    * @return Returns true if string contains substring
    */
    bool ContainsSubstring (_String&);

    /**
    * Converts a good ole char*
    * \n\n \b Example: \code string.toStr(); \endcode
    */
    virtual     BaseRef toStr (unsigned long = 0UL);

    /**
    * Return good ole char*
    */
    virtual     operator const char* (void) const;

    /**
    * Returns a good ole char*
    * \n\n \b Example: \code char * new_str = string.getStr(); \endcode
    * @return Returns a good ole char*
    */
    const char*    getStr(void) const;

 







    /**
    * Replace string 1 with string 2, all occurences true/false
    * \n\n \b Example: \code _String("AAABBBCCC").Replace("BBB","ZZZ") \endcode
    * @param s The substring to replace
    * @param d The substring to replace the value with
    * @param flag If true, replace all.
    * @return "AAAZZZCCC"
    */

    const _String Replace(const _String, const _String, bool) const;


    /**
    * Locate the first non-space character of the string
    * \n\n \b Example: \code _String ("    hyphy").FirstNonSpaceIndex()\endcode
    * @param start Beginning of string search
    * @param end End of string search
    * @param direction Choose between backwards and forwards
    * @return The char of the first non-space, in the example, 'h'.
    * @see FirstNonSpaceIndex()
    */
    char    FirstNonSpace(long start = 0, long end = -1, char direction = 1);

    /**
    * Locate the first non-space character of the string
    * \n\n \b Example: \code _String ("    hyphy").FirstNonSpaceIndex()\endcode
    * @param start Beginning of string search
    * @param end End of string search
    * @param direction Choose between backwards and forwards
    * @return The index of the first non-space, in the example, 4.
    * @see FirstNonSpaceIndex()
    */
    long    FirstNonSpaceIndex(long start = 0, long end = -1, char direction = 1) const;


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
    long    FirstSpaceIndex(long start = 0, long end = -1, char direction = 1) const;

    /**
     * Locate the first non-space character of the string following one or more spaces
     * \n Returns index of first space character
     * \n\n \b Example: \code _String ("h yphy").FirstSpaceIndex()\endcode
     * @param start starting index
     * @param end ending index to search
     * @param direction If the direction is less than 0, then go backwards
     * @return Returns the index of the first non-space. 1 in the example.
     * @sa FirstSpaceIndex()
     */
    long    FirstNonSpaceFollowingSpace(long start = 0, long end = -1, char direction = 1) const;

    /**
    * Finds end of an ID. An ID is made up of alphanumerics, periods, or '_'
    * \n\n \b Example: \code _String ("AA$AAA").FindEndOfIdent()\endcode
    * @param start Where to start looking
    * @param end Where to end looking, -1 is the end of the string
    * @param wild Wild character to skip as well
    * @return Position after the end of the identifier. 3 in the example
    */
    long    FindEndOfIdent(long start = 0, long end = -1, char wild = '*');

    /**
    * Find first occurence of the string between from and to
    * \n\n \b Example: \code _String ("AABBCC").Find("B")\endcode
    * @param s The substring to find
    * @param from The index to start searching from
    * @param to The index to search to
    * @return Returns the index of the first instance of the substr, -1 if not found. 2 in the example
    * @sa FindKMP()
    */
    long    Find(const _String s, long from = 0, long to = -1) const;

    /**
    *  @see Find()
    */
    long    Find(char s, long from = 0, long to = -1) const ;

    /**
    * Find first occurence of the string between from and to
    * @param s The substring to find
    * @param from The index to start searching from
    * @param to The index to search to
    * @return Returns the index of the first instance of the substr, -1 if not found
    * @sa Find()
    * @sa buildKmpTable()
    */
    long    FindKMP(_String s, long from = 0, long to = -1);

    /**
    * Builds a KMP table for use with FindKMP
    */
    void    buildKmpTable(_String s);

    /**
    * Case insensitive Find
    * @see Find()
    */

    long    FindAnyCase (_String, long from = 0, long to = -1);

    /**
    * Backwards Find
    * @see Find()
    */
    long    FindBackwards(_String const&, long, long) const;

    /**
    * Binary searches for a char inside of a string
    * \n\n \b Example: \code _String ("AABBCC").FindBinary('B')\endcode
    * @param s The char to look for inside of the string
    * @return The location of the char, -1 if doesn't exist. 3 in the example.
    */
    long    FindBinary(char);

    /**
    * Compute Adler-32 CRC for a string
    * \n\n \b Example: \code _String result = new _String ("Wikipedia"); \endcode
    * \n Implementation shamelessly lifted from http://en.wikipedia.org/wiki/Adler-32
    * @return the Adler32 checksum. 300286872 returns in the Example
    */
    long    Adler32 (void);

    /**
    * Turns seconds into a time string in the form "hh:mm:ss"
    * \n\n \b Example:
    * \code
    * long time_diff = 459132;
    * _String("").FormatTimeString(time_diff);
    * \endcode
    * @param time_diff Seconds of time
    * @return Transforms string to "127:32:12" in the example.
    */

    void    FormatTimeString (long);


    /**
    * Case Insensitive Lexicographic comparison
    * \n Checks if Strings are equal lexicographic
    * @param s Second string to compare
    * @return true if strings are equal
    * @sa Compare()
    */
    bool iEqual   (_String*);



 

     /**
    * Checks to see if string contains substring
    * \n\n \b Example: \code _String("hyphy").contains("h")\endcode
    * @return Returns true if string contains substring. Example returns true
    * @see Find()
    */
    bool contains (_String);

    /**
    * Checks to see if string contains character
    * @return Returns true if string contains character
    * @see contains()
    */
    bool contains (char);

    /**
    * Checks to see if String begins with substring
    * \n\n \b Example: \code _String("hyphy").beginswith("h")\endcode
    * @param s Substring
    * @param caseSensitive If true, it will be case sensitive. Default is case sensitive.
    * @return true if string begins with substring. Example returns true
    * @sa contains()
    * @sa startswith()
    * @sa endswith()
    */
    bool beginswith (_String const, bool = true) const;

    /**
    * Checks to see if String starts with substring
    * \n\n \b Example: \code _String("hyphy").startswith("h")\endcode
    * @param s Substring
    * @return true if string starts with substring. Example would return true
    * @sa contains()
    * @sa beginswith()
    * @sa endswith()
    */
    bool startswith (_String const&) const;

    /**
     * Checks to see if String starts with substring and it can't be extended to make a valid ident
     * \n\n \b Example: \code _String("return;").startswith("return"); _String("return_me").startswith("return")\endcode
     * @param s Substring
     * @return true if string starts with substring. Example 1 would return true, and example 2 would return false
     * @sa contains()
     * @sa beginswith()
     * @sa endswith()
     */
    bool startswith_noident (_String const&) const;

    /**
    * Checks to see if String ends with substring
    * \n\n \b Example: \code _String("hyphy").endswith("y")\endcode
    * @param s Substring
    * @param caseSensitive If true, it will be case sensitive. Default is case sensitive.
    * @return true if string ends with substring. Example would return true.
    * @sa contains()
    * @sa beginswith()
    * @sa startswith()
    */
    bool endswith (_String, bool = true);

    /**
    * Converts string to upper case
    * @sa LoCase()
    */
    void    UpCase (void);

    /**
    * Converts string to lower case
    * @sa UpCase()
    */
    void    LoCase (void);

 

    /**
    * Returns a list from a string split by a substr
    * \n\n \b Example: _String("hyphy, gattaca, protease").Tokenize(",") will create a list {"hyphy","gattaca","protease"}
    * @param s The substring to split the string by
    * @return A point to a *_List that holds a list of the resultant strings. Retrieve one by list->lData[i]
    */
    const _List  Tokenize (_String const&) const;

    /**
    * TODO: With batchlan
    */
    bool    ProcessFileName (bool isWrite = false, bool acceptStringVars = false, hy_pointer = nil, bool assume_platform_specific = false, _ExecutionList * caller = nil);

    /**
    * TODO: With batchlan
    */
    void    ProcessParameter (void);

    /**
    * Compose two UNIX paths (abs+rel)
    * \n\n \b Example: \code _String("/home/sergei/hyphy").PathComposition("../datamonkey")\endcode
    * @param relPath The relative path to change to
    * @return New File Path, Example would return "/home/sergei/datamonkey"
    */
    _String const PathComposition (_String const) const;

    /**
    * Subtracts the string from the string passed.
    * \n\n \b Example: \code _String("/home/sergei/hyphy").PathSubtraction("/home/sergei/")\endcode
    * @param s String that will be subtracted
    * @return Example would return "hyphy"
    */
    _String const PathSubtraction (_String const, char) const;

    /**
    * Strips quotes from around the string if present (in place)
    * \n\n \b Example: \code _String("\"hyphy\"").StripQuotes("")\endcode
    */
    void    StripQuotes(void);

    /**
     * Decorates the string with quotes

     * @param quote_char which character to use as a "quote"
     * @return quote_char + *this + quote_char
     */
    const    _String Enquote (char quote_char = '\'') const;

    /**
    * Checks if String is valid ident
    * \n A valid ident is any alphanumeric or '_'
    * \n\n \b Example: '$hyphy' is not legal.  'hy_phy' is legal.
    * @param strict If strict, only alphabetic, no numerals.
    * @sa ConvertToAnIdent();
    */
    bool    IsValidIdentifier(bool = true) const;

    /**
    * Same as IsValidIdentifier, but must end with a '&'
    * \n\n \bExample: 'hyphy&' is a valid ref identifier
    * @see IsValidIdentifier()
    */
    bool    IsValidRefIdentifier(void) const;

    /**
    * If it is enclosed in quotes, then it is a literal argument
    * \n \n \b Example: "\"hyphy\"" is a literal argument
    */
    bool    IsALiteralArgument  (bool stripQuotes = FALSE);

    /**
    * Converts a string to a valid ident
    * \n A valid ident is any alphanumeric or '_'
    * \n\n \b Example: \code _String("$hyphy") \endcode
    * @param strict If strict, only alphabetic, no numerals.
    * @sa IsValidIdentifier();
    * @return the example would return "_hyphy"
    */
    void    ConvertToAnIdent (bool = true);

    /**
    * Removes all spaces in a string
    * \n\n \b Example: \code _String("   h  y p    h  y").KillSpaces \endcode
    * @param result A reference to the string that will have stripped spaces. The original string will not be changed.
    * @sa CompressSpaces()
    * @return The example would return "hyphy"
    */
    void    KillSpaces       (_String&);

    /**
    * Removes all spaces in a string
    * \n\n \b Example: \code _String("   h  y p    h  y").CompressSpaces() \endcode
    * @return Example would transform the string to "h y p h y"
    * @sa KillSpaces()
    */
    void    CompressSpaces   (void);

    /**
    * Shorten the var id by removing the matching beginning portion of the passed in String separated by a "."
    * \n\n \b Example: _String("house.room")._String("house")
    * @param containerID
    * @return _String("house.room")._String("house.") will return "room". However, _String("houseroom") will return "houseroom"
    */
    _String ShortenVarID     (_String&);

    /**
    * Examine the string argument contained in this object, decide what it is, and process accordingly
    * \n\n \bExample: \code 'hyphy'.ProcessVariableReferenceCases (object) \endcode is a direct reference to object hyphy
    * \n\n \bExample: \code '\"hy\"+\"phy\"'.ProcessVariableReferenceCases (object) \endcode is a direct reference to object hyphy
    * \n\n \bExample: \code '*hyphy'.ProcessVariableReferenceCases (object) \endcode is a reference to the object whose name is stored in the string variable hyphy 
    * \n\n \bExample: \code '**hyphy'.ProcessVariableReferenceCases (object) \endcode is a reference to the object whose name is stored in the string variable hyphy in the global context
    * @param referenced_object will store the handled variable ID
    * @param context is the namespace of the referenced object; could be nil
    * @return one of HY_STRING_INVALID_REFERENCE    HY_STRING_DIRECT_REFERENCE   HY_STRING_LOCAL_DEREFERENCE    HY_STRING_GLOBAL_DEREFERENCE 
    * @see IsValidIdentifier()
    */

    unsigned char  ProcessVariableReferenceCases (_String& referenced_object, _String const * context = nil) const;

    

    /**
    * A regular expression match
    * @param pattern A string(A hy_pointer is a char*) that holds the regex pattern
    * @param matchedPairs A list that holds the start and end indexes of matches
    * @sa RegExpMatchAll()
    * @sa RegExpMatchOnce()
    */
    void    RegExpMatch      (hy_pointer, _SimpleList&);

    /**
    * A regular expression match
    * @param pattern A string that holds the regex pattern
    * @param matchedPairs A list that holds the start and end indexes of matches
    * @sa RegExpMatch()
    * @sa RegExpMatchOnce()
    */
    void    RegExpMatchAll   (hy_pointer, _SimpleList&);

    /**
    * A regular expression match that only matches once
    * @param pattern A string that holds the regex pattern
    * @param matchedPairs A list that holds the start and end indexes of matches
    * @sa RegExpMatch()
    * @sa RegExpMatchAll()
    */
    void    RegExpMatchOnce  (_String*, _SimpleList&, bool, bool);

    /**
    * Lexicographically sorts the string
    * @param index Needs a list to act as an index
    * @return sorted string
    */
    _String*Sort             (_SimpleList* = nil);

    /**
     * Generate a random string on 
     * @param len (>0) The desired length of the string
     * @param alphabet Which alphabet do the random charcters come from; in nil, then this will be generated from 1-128 ASCII codes 
     * @return the random string
     */
    static _String Random             (const unsigned long len, const _String * alphabet = nil);

    /**
    * Computes Lempel-Ziv complexity of the string.
    * \n The Lempel-Ziv complexity computes the number of separate substrings in a given string
    * \n Example: 1001111011000010 = 6 because subset of all strings = {1, 0, 01, 1110, 1100, 0010 }
    */
    long    LempelZivProductionHistory
    (_SimpleList* = nil);

    /**
    * Converts a string of form ":[\d\.]\+" into a double
    * \n SLKP 20100831: a utility function to handle the
    * conversion of branch length strings to parameters
    * \n\n \b Example: string = ":3.14" returns 3.14
    * \n If it is not of correct form, it will return 1e-10
    * @sa toNum()
    */
    hy_float
    ProcessTreeBranchLength ();



    /**
    * Converts a string of form "[\d\.]\+" into a double
    * \n\n \b Example: "3.14" becomes 3.14
    * @sa ProcessTreeBranchLength()
    */

    hy_float      toNum (void) const;

    /**
     * Converts a into a long
     * \n\n \b Example: "3.14" becomes 3
     */
    
    long      toLong (void) const;

  /**
    * Sets Length
    */
    void    SetLength (unsigned long nl) {
        s_length=nl;
    }

    /**
    * Starting at index [argument 1],
    * find a span that encloses an expression (nested) delimited by char[argument 2]
    * and char[argument 3] (e.g. {}, ()) respecting quotes (argument 4), and allowing
    * escaped characters (argument 5)
    * \n SLKP 20090803
    *
    * @param &from The starting position of the segment will be stored here
    * @param open The first character to look for. For example, and open bracket '[' or open paranthesis '('
    * @param close The first character to look for. For example, and open bracket ']' or open paranthesis ')'
    * @param respectQuote
    * @param respectEscape
    *
    * @return Ending position is returned
    *-1 is returned if the starting character could not be found or the expression did not terminate before the end of the string
    *
    */

    long    ExtractEnclosedExpression (long&, char, char, bool, bool);

    /**
    * Starting at index [argument 1],
    * find a span that terminates in one of the characters in [argument 2], while
    * respecting (), [], {}, "" and escapes
    * \n SLKP 20090805
    * @param s The terminator to find
    * @return -1 is returned if the starting character could not be found or the expression did not terminate before the end of the string
    * @sa IsALiteralArgument()
    *
    */

    long    FindTerminator          (long, _String const&) const;

    
protected:
    unsigned long s_length;
    char*         s_data;
    
private:
    
    const static  char   default_return = '\0';
        /** this value is returned for "failed" 
            access operations that don't throw errors, e.g. getChar */
    
    inline void AllocateAndCopyString (const char * source_string, unsigned long length);
        /** this is a utility function which allocates length+1 chars for s_data, copies
            the data from source_string, and sets the terminating 0 
         
         * Revision history
            - SLKP 20170517 factoring repeated functionality
         
         */
         
    long NormalizeRange (long & start, long & end) const;
    
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

};

/** DEPRECATED
 - virtual     void    DuplicateErasing (BaseRef);
 SLKP 20170517 ::Duplicate now clears *this always
*/


// _______________________________________________________________________



extern _String 
       emptyAssociativeList,
       hyphyCiteString;


void    SetStatusBarValue           (long,hy_float,hy_float);
void    SetStatusLine               (_String);
void    SetStatusLine               (_String, _String, _String, long l);
void    SetStatusLine               (_String, _String, _String);
void    SetStatusLine               (_String, _String, _String, long, char);

void    SetStatusLineUser           (_String const);


hy_pointer     PrepRegExp                  (_String*, int&, bool);
void    FlushRegExp                 (hy_pointer);
_String GetRegExpError              (int);
_String GetVersionString            (void);
_String GetTimeStamp                (bool = false);

void    StringToConsole             (_String const, void * extra = nil);
void    BufferToConsole             (const char*, void * extra = nil);
void    NLToConsole                 (void * extra = nil);

_String*StringFromConsole           (bool=true);

char    GetPlatformDirectoryChar    (void);


extern  _String                     __HYPHY__VERSION__;

#ifdef __UNIX__
	extern bool	needExtraNL;
#endif

typedef bool (*_hyStringValidatorType) (_String*);
bool    hyIDValidator (_String*);


/** REMOVED
 - virtual     void    DuplicateErasing (BaseRef);
        SLKP 20170517 ::Duplicate now clears *this always
 - void    CopyDynamicString (_String* s, bool = true);
    SLKP 20170517 the same can be accomplished by 'x = s' and constructor elision
 - bool    iEqual 
    SLKP 20170517 replace with a more general CompareIgnoringCase
 

 */

/**
 * Append operator
 * \n\n \b Example: \code _String new_string = _String("A") & _String("B") \endcode
 * @return "AB"
 * @sa EscapteAndAppend()
 */
// _String & operator << (const _String*); // MOVE TO STRING BUFFER

/**
 * Append operator
 */
//_String & operator << (const _String&); // MOVE TO STRING BUFFER

/**
 * Append operator
 * \n\n \b Example: \code _String new_string = _String("A") & _String("B") \endcode
 * @return "AB"
 * @sa EscapteAndAppend()
 */
//void    AppendNewInstance (_String*); // MOVE TO STRING BUFFER

/**
 * Append multiple copies of the same string to the buffer
 * @param value the string to copy
 * @param copies how many copies to make
 */
//void    AppendNCopies   (_String const& value, unsigned long copies); // MOVE TO STRING BUFFER

/**
 * Append operator
 * \n\n \b Example: \code _String new_string = _String("A") & _String("B") \endcode
 * @return "AB"
 * @sa AppendNewInstance()
 */
//_String& operator << (const char); // MOVE TO STRING BUFFER

/**
 * Escape all characters in a string and append to this string
 * \n\n \b Example: \code _String("AB").EscapeAndAppend('<',4); \endcode
 * \n Above code will transform string to "AB&lt;"
 * @param c The character to escape and append
 * @param mode What sort of escaping
 * \n mode = 0 : normal "text" escaping
 * \n mode = 1: PostScript escaping
 * \n mode = 2: SQLite escaping
 * \n mode = 3: SQLite escaping
 * \n mode = 4: HTML escaping
 * \n mode = 5: Regexp escaping
 */
//virtual void EscapeAndAppend (const char, char); // MOVE TO STRING BUFFER

/**
 * Escape all characters in a string and append to this string
 * \n\n \b Example: \code _String("AB").EscapeAndAppend('<',4); \endcode
 * \n Above code will transform string to "AB&lt;"
 * @param s The string to escape and append
 * @param mode What sort of escaping
 * @see EscapeAndAppend(const char, char)
 */
//virtual void EscapeAndAppend (const _String &, char mode = 0); // MOVE TO STRING BUFFER

/**
 * Append into operator
 */
//_String& operator << (const char*); // MOVE TO STRING BUFFER

/**
 * Finalizes a string by putting a 0 at the end of the string.
 */
//virtual void Finalize (void); // MOVE TO STRING BUFFER

/**
 * SLKP 20090817: A utility function to append a statement of the form
 * \n\n \b Example: _String("hyphy").AppendAnAssignmentToBuffer("12","12",false,false,false) makes "hyphy12=12"
 * @param id = value; to the current string assumed to be in the buffer form
 * @param flags: a bitwise combination of flags; set kAppendAnAssignmentToBufferFree to free 'value'; \\
 set kAppendAnAssignmentToBufferQuote to put quotes around the value \\
 set kAppendAnAssignmentToBufferAssignment to use ':=' instead of '=' \\
 default is to use kAppendAnAssignmentToBufferFree
 * @sa AppendNewInstance()
 * @sa AppendVariableValueAVL()
 */

// void    AppendAnAssignmentToBuffer (_String*, _String*, unsigned long = kAppendAnAssignmentToBufferFree); // MOVE TO STRING BUFFER

/**
 * SLKP 20090817:
 * A utility function to append a statement of the form
 * id["varname"] = varvalue; for each variable in the SimpleList arguments
 * for String valued variables, their values are properly quoted
 * @param id = value; to the current string assumed to be in the buffer form
 * @param doFree free the 2nd string argument when done
 * @param doQuotes put quotes around the value
 * @param doBind use := instead of =
 * @sa AppendNewInstance()
 * @sa AppendAnAssignmentToBuffer()
 */

//void    AppendVariableValueAVL (_String*, _SimpleList&);// MOVE TO STRING BUFFER
#endif
