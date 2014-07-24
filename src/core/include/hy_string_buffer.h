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

class _StringBuffer : public _String {

private:

  unsigned long sa_length;

  void allocateBufferSpace(const unsigned long character_count);
  void resizeString(void);
  void pushChar(const char);
  void pushCharBuffer(const char*, const unsigned long);

public:

  /**
   * A constructor that creates a string buffer of a default size.
   */
  _StringBuffer(void);

  /**
   * A constructor that creates a string buffer of a given size.
   * @param character_count The number of characters worth of space to
   * allocate
   */
  _StringBuffer(const unsigned long character_count);

  /**
  * Stack copy.
  */
  _StringBuffer(const _StringBuffer &);

  /**
  * A constructor that copies from a char string.
   * @param buffer Create buffer from provided char array
  */
  _StringBuffer(const char * buffer);

  /**
   * A constructor that copies from a standard string.
   * @param buffer Create buffer from provided HyPhy _String
   */
  _StringBuffer(const _String& buffer);

  /**
   * Initializes _String object to 0 allocated length
   * @param Allocate memory or not
   */
  virtual void initialize(bool=false);

  virtual ~_StringBuffer(void) {};

  /**
  * Returns a dynamic copy of the current instance.
  * \n Usage: string_buffer_instance.makeDynamic();
  * @return BaseRef
  */
  virtual BaseRef makeDynamic(void) const;

  /**
  * Duplicates a _StringBuffer
  * \n Usage: string_buffer.duplicate(&existing_string_buffer)
  * @param ref A pointer to the _StringBuffer to be duplicated
  */
  virtual void duplicate(BaseRefConst);

  /**
  * Append all characters in the argument string to the buffer
  * @param buffer Append characters from standard hyphy _String
  */
  virtual void operator<<(const _String *buffer);

  /**
   * Append all characters in the argument string to the buffer
   * @param buffer append characters from standard HyPhy _String
  */
  virtual void operator<<(const _String &);

  /**
   * Append a single char to the buffer
   * @param buffer append characters from character
   */
  virtual void operator<<(const char);

  /**
   * Append all chars in the string buffer to this string
   * @param buffer append characters from char array
   */
  virtual void operator<<(const char *);

  /**
   * Append all characters in the argument string to the buffer
   * delete the buffer object afterwards
   * @param buffer append characters from stanard HyPhy _String
   */
  void appendNewInstance(_String *);

  /**
   * Append a substring of the source string to this buffer
   */
  void appendSubstring(const _String &, long from, long to);

  /**
  * MDS 20140722: Sanitize (escape) all characters in a string and append to
  * this string \n\n \b Example: \code
  * _StringBuffer("AB").sanitizeForHTMLAndAppend('<'); \endcode \n Above
  * code will transform string to "AB&lt;"
  * @param c The character to escape and append
  */
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
   * SLKP 20090817: A utility function to append a statement of the form
   * \n\n \b Example:
   * _String("hyphy").AppendAnAssignmentToBuffer("12","12",false,false,false)
   * makes "hyphy12=12"
   * @param id = value; to the current string assumed to be in the buffer form
   * @param do_free free the 2nd string argument when done
   * @param do_quotes put quotes around the value
   * @param do_bind use := instead of =
   * @sa appendNewInstance()
   * @sa appendVariableValueAVL()
   */
  void appendAnAssignmentToBuffer(_String *, _String *, bool = true,
                                  bool = false, bool = false);
};

#endif

