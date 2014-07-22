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

#ifndef _HY_STRING_BUFFER_
#define _HY_STRING_BUFFER_
//#pragma once

#define   HY_STRING_BUFFER_ALLOCATION_CHUNK 16UL

#include "hy_strings.h"

enum _hyStringBufferEscapeMode   {HY_ESCAPE_NORMAL, // used to be 0
                                  HY_ESCAPE_POSTSCRIPT, // 1
                                  HY_ESCAPE_SQLITE, // 2
  HY_ESCAPE_UNUSED, //3
  HY_ESCAPE_HTML, //4
  HY_ESCAPE_REGEXP//5
};


class _StringBuffer : public _String {

private:

  unsigned long saLength;

  void AllocateBufferSpace (const unsigned long character_count);
  void ResizeString        (void);
  void PushChar            (const char);
  void PushCharBuffer      (const char*, const unsigned long);

public:

  /**
   * A constructor that creates a string buffer of a default size.
   */
  _StringBuffer(void);

  /**
   * A constructor that creates a string buffer of a given size.
   * @param character_count The number to convert to string
   */
  _StringBuffer(const unsigned long character_count);

  /**
  * Stack copy.
  */
  _StringBuffer(const _StringBuffer &);


  /**
  * A constructor that copies from a char string.
   * @param buffer create
  */
  _StringBuffer(const char * buffer);

  /**
   * A constructor that copies from a standard string.
   * @param buffer create
   */
  _StringBuffer(const _String& buffer);

  /**
   * Initializes _String object to 0 allocated length
   */
  virtual void Initialize(void);

  virtual ~_StringBuffer(void) {};

  /**
  * Returns a dynamic copy of the current instance.
  * \n Usage: stringInstance.makeDynamic();
  * @return BaseRef
  */
  virtual BaseRef makeDynamic(void);

  /**
  * Duplicates a string
  * \n Usage: string.Duplicate(&existing_string)
  * @param ref A pointer to the string to be duplicated
  */
  virtual void Duplicate(BaseRefConst);


  /**
  * Append all characters in the argument string to the buffer
  * @param buffer append characters from here
  */
  virtual void operator<<(const _String * buffer);

  /**
   * Append all characters in the argument string to the buffer
   * @param buffer append characters from here
  */
  virtual void operator<<(const _String &);


  /**
   * Append all characters in the argument string to the buffer
   * delete the buffer object afterwards
   * @param buffer append characters from here
   */
  void AppendNewInstance(_String *);

  /**
   * Append a single char to the buffer
   * @param buffer append characters from here
   */
  virtual void operator<<(const char);

  void sanitizeForSQLAndAppend(const char);
  void sanitizeForSQLAndAppend(const _String&);
  void sanitizeForHTMLAndAppend(const char);
  void sanitizeForHTMLAndAppend(const _String&);
  void sanitizeAndAppend(const char);
  void sanitizeAndAppend(const _String&);
  void sanitizeForPostScriptAndAppend(const char);
  void sanitizeForPostScriptAndAppend(const _String&);
  void sanitizeForRegExAndAppend(const char);
  void sanitizeForRegExAndAppend(const _String&);

  /**
  * Escape all characters in a string and append to this string
  * \n\n \b Example: \code _StringBuffer("AB").EscapeAndAppend('<',4); \endcode
  * \n Above code will transform string to "AB&lt;"
  * @param c The character to escape and append
  * @param mode What sort of escaping (see _hyStringBufferEscapeMode)
  */
  virtual void EscapeAndAppend(const char, const _hyStringBufferEscapeMode);

  /**
  * Escape all characters in a string and append to this string
  * \n\n \b Example: \code _String("AB").EscapeAndAppend('<',4); \endcode
  * \n Above code will transform string to "AB&lt;"
  * @param s The string to escape and append
  * @param mode What sort of escaping
  * @see EscapeAndAppend(const char, char)
  */
  virtual void EscapeAndAppend(const _String &, const _hyStringBufferEscapeMode = HY_ESCAPE_NORMAL);

  /**
   * Append all chars in the string buffer to this string
   * @param buffer append characters from here
   */
  virtual void operator<<(const char *);

  /**
   * SLKP 20090817: A utility function to append a statement of the form
   * \n\n \b Example:
   * _String("hyphy").AppendAnAssignmentToBuffer("12","12",false,false,false)
   * makes "hyphy12=12"
   * @param id = value; to the current string assumed to be in the buffer form
   * @param doFree free the 2nd string argument when done
   * @param doQuotes put quotes around the value
   * @param doBind use := instead of =
   * @sa AppendNewInstance()
   * @sa AppendVariableValueAVL()
   */

  void AppendAnAssignmentToBuffer(_String *, _String *, bool = true,
                                  bool = false, bool = false);

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

  //void AppendVariableValueAVL(_String *, _SimpleList &);
  void AppendVariableValueAVL(_String *, _List &);

};

#endif

