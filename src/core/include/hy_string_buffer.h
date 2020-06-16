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

#include "hy_strings.h"


#define   HY_STRING_BUFFER_ALLOCATION_CHUNK 16UL
#define kAppendAnAssignmentToBufferPlain        0x00
#define kAppendAnAssignmentToBufferFree         0x01
#define kAppendAnAssignmentToBufferQuote        0x02
#define kAppendAnAssignmentToBufferAssignment   0x04
#define kAppendAnAssignmentToBufferGlobal       0x08

#define _HY_STRING_BUFFER_PREALLOCATE_SLOTS     1024

class _StringBuffer : public _String {

private:

  /**
      Allocated buffer size (>= s_length)
   
   */
  unsigned long sa_length;

   /** Allocate buffer to hold a specified number of characters
   
      @param character_count: the number of characters in the buffer
   *  Revision history
      - SLKP 20170613 reviewed while porting from v3 branch
   */
  void AllocateBufferSpace(const unsigned long character_count);

  /** Trim unused characters in the buffer
   
   *  Revision history
   - SLKP 20170613 reviewed while porting from v3 branch
   */
  void ResizeString(void);

  /** Add a character to buffer; resize if needed
   
      @param c: the character to add
   *  Revision history
   - SLKP 20170613 reviewed while porting from v3 branch
   */
  void PushChar(const char c);

  /** Add a character string to buffer; resize if needed
   
   @param str: the char* to add to buffer
   @param size: the length of the buffer (no checks performed)
   
   *  Revision history
   - SLKP 20170613 reviewed while porting from v3 branch
   */
  void PushCharBuffer(const char* str, const unsigned long size);

public:

  /**
   * A constructor that creates a string buffer of a default size
     (HY_STRING_BUFFER_ALLOCATION_CHUNK)
 
   *  Revision history
      - SLKP 20170613 reviewed while porting from v3 branch
   */
  _StringBuffer(void);

  /**
   * A constructor that creates a string buffer of a given size.
   * @param character_count The number of characters worth of space to
   * allocate
 
   *  Revision history
   - SLKP 20170613 reviewed while porting from v3 branch
   */
  _StringBuffer(const unsigned long character_count);

  /**
  * Stack copy.
   
   *  Revision history
      - SLKP 20170613 reviewed while porting from v3 branch
  */
  _StringBuffer(const _StringBuffer &);

  /**
  * A constructor that copies from a char string.
   * @param buffer Create buffer from provided char array
 
   *  Revision history
   - SLKP 20170613 reviewed while porting from v3 branch
   */
  _StringBuffer(const char * buffer);

  /**
   * A constructor that copies from a standard string.
   * @param buffer Create buffer from provided HyPhy _String

   *  Revision history
   - SLKP 20170613 reviewed while porting from v3 branch
   */
  _StringBuffer(const _String& buffer);

    /**
     * A constructor that moves from a standard string.
     * @param buffer Create buffer from provided HyPhy _String
     
     *  Revision history
     - SLKP 20190507 initial implementation
     */
   _StringBuffer(_String&& buffer);

  /**
   * A constructor that moves data from a standard string.
   * @param buffer This string will be deleted upon return
   
   *  Revision history
   - SLKP 20170920 initial implementation
   */
  _StringBuffer(_String* buffer);

  /**
   *  Default constructor initializer
   *  Revision history
   - SLKP 20170613 reviewed while porting from v3 branch
   */
  virtual void Initialize(bool = true);

  /**
   * Clears this object
   *  Revision history
   - SLKP 20170613 reviewed while porting from v3 branch
  */
  virtual void Clear();

  /**
   * Empties this object without clearing the memory
   *  Revision history
   - SLKP 20200224 Initial Implementation
  */
    
  virtual void Reset();

  /**
   * Standard destructor
   *  Revision history
   - SLKP 20170614 reviewed while porting from v3 branch
   */

  virtual ~_StringBuffer(void);

  /**
  * Returns a dynamic copy of the current instance.
  * \n Usage: string_buffer_instance.makeDynamic();
  * @return BaseRef
  *  Revision history
    - SLKP 20170613 reviewed while porting from v3 branch
  */
  virtual BaseRef makeDynamic (void) const;

  /**
  * Duplicates a _StringBuffer;
  * \n Usage: string_buffer.duplicate(&existing_string_buffer)
  * @param ref A pointer to the _StringBuffer to be duplicated
   *  Revision history
   - SLKP 20170614 reviewed while porting from v3 branch
  */
  void Duplicate (BaseRefConst);
    
    
  virtual void Trim(long, long);

    
 /**
  * Move semantics for buffer assignment
  * @param rhs A pointer to the _StringBuffer to be moved from
  *  Revision history
   - SLKP 20190507 initial implementations
  */
  _StringBuffer& operator = (_StringBuffer&&rhs);

  /**
  * Append all characters in the argument string to the buffer
  * @param buffer Append characters from standard hyphy _String
           if buffer == null, does nothing
   *  Revision history
   - SLKP 20170614 reviewed while porting from v3 branch
     [CHANGE-NOTE SLKP 20170614 all << operators return *this for chaining]
  */
  _StringBuffer& operator<<(const _String *buffer);

  /**
   * Append all characters in the argument string to the buffer
   * @param buffer append characters from standard HyPhy _String
   *  Revision history
    - SLKP 20170614 reviewed while porting from v3 branch
      [CHANGE-NOTE SLKP 20170614 all << operators return *this for chaining]
  */
  _StringBuffer& operator<<(const _String & buffer);

  /**
   * Append a single char to the buffer
   * @param buffer append this character
   *  Revision history
   - SLKP 20170614 reviewed while porting from v3 branch
     [CHANGE-NOTE SLKP 20170614 all << operators return *this for chaining]
   */
  _StringBuffer& operator<<(const char);

  /**
   * Append all chars in the string buffer to this string
   * @param buffer append characters from char array (assumed non-null and \0 terminated
   *  Revision history
   - SLKP 20170614 reviewed while porting from v3 branch
    [CHANGE-NOTE SLKP 20170614 all << operators return *this for chaining]
  */
  _StringBuffer& operator<<(const char *);

  /**
   * Append all characters in the argument string to the buffer
   * delete the buffer object afterwards (respecting reference counts)
   * @param buffer append characters from stanard HyPhy _String
   *  Revision history
   - SLKP 20170614 reviewed while porting from v3 branch
   [CHANGE-NOTE SLKP 20170614 return *this for chaining]
   */
  _StringBuffer& AppendNewInstance(_String *);

  /**
   * Append multiple copies of the same string to the buffer
   * @param value the string to copy
   * @param copies how many copies to make
   *  Revision history
   - SLKP 20170614 reviewed while porting from v2.3 branch
   [CHANGE-NOTE SLKP 20170614 return *this for chaining]
   */
  _StringBuffer& AppendNCopies   (_String const& value, unsigned long copies);


  /**
   * Append a substring of the source string to this buffer
   * @param source the string to copy from
   * @param start 0-based index to start the substring at
   * @param end 0-based index to end the substring at
   
   *  Revision history
   - SLKP 20170614 reviewed while porting from v2.3 branch
   [CHANGE-NOTE SLKP 20170614 return *this for chaining]
   */
  _StringBuffer& AppendSubstring(_String const& source, long start, long end);

  
  
  /**
   A suite of functions useful for escaping specific characters and pushing them
   onto a buffer, either a single character at a time, or an entire string at a time
   
   *  Revision history
   - SLKP 20170614 reviewed while porting from v2.3 branch
   [CHANGE-NOTE SLKP 20170614 return *this for chaining]
  */
  _StringBuffer& SanitizeForSQLAndAppend(const char);
  _StringBuffer& SanitizeForSQLAndAppend(const _String&);
  _StringBuffer& SanitizeForHTMLAndAppend(const char);
  _StringBuffer& SanitizeForHTMLAndAppend(const _String&);
  _StringBuffer& SanitizeAndAppend(const char);
  _StringBuffer& SanitizeAndAppend(const _String&);
  _StringBuffer& SanitizeForPostScriptAndAppend(const char);
  _StringBuffer& SanitizeForPostScriptAndAppend(const _String&);
  _StringBuffer& SanitizeForRegExAndAppend(const char);
  _StringBuffer& SanitizeForRegExAndAppend(const _String&);

  /**
   * A utility function to append a statement of the form
   * \n\n \b Example: _String("hyphy").AppendAnAssignmentToBuffer("12","12",false,false,false) makes "hyphy12=12"
   * @param id = value; to the current string assumed to be in the buffer form
   * @param flags: a bitwise combination of flags; set kAppendAnAssignmentToBufferFree to free 'value'; \\
                   set kAppendAnAssignmentToBufferQuote to put quotes around the value \\
                   set kAppendAnAssignmentToBufferAssignment to use ':=' instead of '=' \\
                   default is to use kAppendAnAssignmentToBufferFree
   * @sa AppendNewInstance()
   * @sa AppendVariableValueAVL()
   *  Revision history
   - SLKP 20170614 reviewed while porting from v3 branch
   [CHANGE-NOTE SLKP 20170614 added the kAppendAnAssignmentToBufferPlain option]
   */
   void    AppendAnAssignmentToBuffer (_String const* id , _String * value, unsigned long = kAppendAnAssignmentToBufferFree);

  /**
   * A utility function to append a statement of the form
   * id["varname"] = varvalue; for each variable (index) in the SimpleList arguments
   * for String valued variables, their values are properly quoted
   * @sa AppendNewInstance()
   * @sa AppendAnAssignmentToBuffer()
   *  Revision history
   - SLKP 20170614 reviewed while porting from v3 branch
   [CHANGE-NOTE SLKP 20170614 changed first argument type to const]
   */
  
   void    AppendVariableValueAVL (_String const*, _SimpleList const&);

    /**
     *  Trim the buffer allocation to space used. 
     *  Revision history
     -  SLKP 20170923 initial implementation
     */
    virtual void TrimSpace (void);
    
    
    /// memory buffering
    void * operator new       (size_t size);
    void   operator delete    (void * p);

    static  _SimpleList                    free_slots;
    static  unsigned char                  preallocated_buffer[];

    
};

#endif

