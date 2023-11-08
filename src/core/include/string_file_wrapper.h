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

#ifndef _StringFileWrapper_
#define _StringFileWrapper_


#include "hy_strings.h"
#include "hy_string_buffer.h"

enum StringFileWrapperConstants {
    kStringFileWrapperNewLine,
    kStringFileWrapperLinefeed,
    kStringFileWrapperTab
};


class StringFileWrapper {
  /** This is a simple convenience flag that unifies << operations
    *  when the LHS is either a _String, or a FILE*
   */
public:
  
  StringFileWrapper (_StringBuffer * string, hyFile * file);
  /** Create a wrapper around around a string / file pair 
      If both arguments are null, the wrapper will simply "eat" the 
      bufferring operations (/dev/null equivalent). If both arguments
      are NON-null, string takes precedence.
      
      @param string if not NULL, the wrapper will paste all the arguments
      into the string
      @param file if not NULL, the wrapper will write all the arguments to file
   
   */
  
  ~StringFileWrapper () {}
  
  StringFileWrapper & operator << (const char* buffer);
  /** Write a literal string to the underlying buffer 
        
      @param buffer the buffer to write
      @return this for chaining
   */
  
  StringFileWrapper & operator << (const char letter);
  /** Write a character to the underlying buffer
   
   @param letter the character to write
   @return this for chaining
   */
  
  StringFileWrapper & operator << (const _String& buffer);
  /** Write _String to the underlying buffer
   
   @param buffer the string to write
   @return this for chaining
   */
  
  StringFileWrapper & operator << (const _String* buffer);
  /** Write _String* to the underlying buffer
   
   @param buffer the string to write
   @return this for chaining
   */

  StringFileWrapper & operator << (const StringFileWrapperConstants& constant);
  /** Write a special character to the underlying buffer
   
   @param the designation of the character to write
   @return this for chaining
   */
  
private:
  _StringBuffer * string_buffer;
  hyFile*     file_buffer;
  
};

#endif
