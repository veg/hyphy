/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
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

#include "hy_string_buffer.h"
#include "global_things.h"
#include "parser.h"


using namespace hy_global;

#include <string.h> // for strlen


/*
==============================================================
Utility Functions
==============================================================
*/

void _StringBuffer::AllocateBufferSpace(const unsigned long character_count) {
  sa_length = character_count;
  if (s_data) {
    s_data = (char *)MemReallocate(s_data, sa_length + 1UL);
  } else {
    s_data = (char *)MemAllocate(sa_length + 1UL);
  }

}

//=============================================================

void _StringBuffer::ResizeString(void) {
  if (sa_length < s_length) {
    unsigned long inc_by = s_length - sa_length;

    if (inc_by < HY_STRING_BUFFER_ALLOCATION_CHUNK) {
      inc_by = HY_STRING_BUFFER_ALLOCATION_CHUNK;
    }

    if (inc_by < s_length >> 4) {
      inc_by = s_length >> 4;
    }

    this->AllocateBufferSpace(sa_length + inc_by);
  }
}

//=============================================================

void _StringBuffer::Initialize(bool) {
  sa_length = 0UL;
}

//=============================================================

void _StringBuffer::TrimSpace (void) {
    if (sa_length > s_length) {
        AllocateBufferSpace (s_length);
    }
}

//=============================================================

void _StringBuffer::Clear() {
  _String::Clear();
  sa_length = 0UL;
}


/*
==============================================================
Constructors
==============================================================
*/

_StringBuffer::_StringBuffer(void) : _String () {
  this->AllocateBufferSpace(HY_STRING_BUFFER_ALLOCATION_CHUNK);
  s_data[0] = '\0';
}

//=============================================================

_StringBuffer::_StringBuffer(const unsigned long character_count)
                            : _String () {
  this->AllocateBufferSpace(character_count);
  s_data[0] = '\0';
}

//=============================================================

_StringBuffer::_StringBuffer(const char * buffer) : _String () {
  this->Initialize();
  (*this) << buffer;
}


//=============================================================

_StringBuffer::_StringBuffer(const _String& buffer) : _String () {
  this->Initialize();
  (*this) << buffer;
}

//=============================================================
_StringBuffer::_StringBuffer(const _StringBuffer &s): _String () {
  this->Initialize();
  this->Duplicate (&s);
}

  //=============================================================
_StringBuffer::_StringBuffer(_String* buffer): _String (buffer) {
  sa_length = s_length;
}


/*
==============================================================
Cloners and Copiers
==============================================================
*/

void _StringBuffer::Duplicate (BaseRefConst ref) {
  _String::Duplicate (ref);
  this->sa_length = this->s_length;
}

//=============================================================


BaseRef _StringBuffer::makeDynamic (void) const {
  _StringBuffer * r = new _StringBuffer;
  r->Duplicate(this);
  return r;
}

/*
==============================================================
Operator Overloads
==============================================================
*/


// Append operator
_StringBuffer& _StringBuffer::operator<<(const _String *s) {
  if (s && s->length()) {
    this->PushCharBuffer(s->get_str(), s->length());
  }
  return *this;
}

//=============================================================
_StringBuffer& _StringBuffer::operator<<(const _String &s) {
  (*this) << &s;
  return *this;
}

//=============================================================
_StringBuffer& _StringBuffer::operator<<(const char *str) {
  this->PushCharBuffer(str, strlen(str));
  return *this;
}

//=============================================================
_StringBuffer& _StringBuffer::operator<<(const char c) {
  this->PushChar(c);
  return *this;
}

/*
==============================================================
Methods
==============================================================
*/

void _StringBuffer::PushChar(const char c){
  s_length++;
  this->ResizeString();
  s_data[s_length-1UL] = c;
  s_data[s_length] = '\0';
}

//=============================================================

void _StringBuffer::PushCharBuffer( const char* buffer,
                                    const unsigned long buffer_l) {
  if (buffer_l) {
    unsigned long offset = s_length;
    s_length += buffer_l;
    this->ResizeString();

    for (unsigned long k = 0UL; k < buffer_l; k++) {
      s_data[offset + k] = buffer[k];
    }
    s_data[s_length] = '\0';
  }
}

//=============================================================

_StringBuffer& _StringBuffer::AppendNewInstance(_String *s) {
  (*this) << s;
  DeleteObject(s);
  return *this;
}

//=============================================================

_StringBuffer&    _StringBuffer::AppendNCopies   (_String const& value, unsigned long copies) {
  for (unsigned long i = 0UL; i < copies; i++) {
    (*this) << value;
  }
  return *this;
}

//=============================================================

_StringBuffer& _StringBuffer::AppendSubstring(const _String& source, long start, long end) {
  (*this) << _String (source, start, end);
}

 

//=============================================================

_StringBuffer& _StringBuffer::SanitizeForSQLAndAppend(const char c) {
  this->PushChar(c);
  if (c == '\'') {
    this->PushChar(c);
  }
  return *this;
}

_StringBuffer& _StringBuffer::SanitizeForSQLAndAppend(const _String& s) {
  unsigned long sl = s.length ();
  for (unsigned long i = 0UL; i < sl; i++) {
    this->SanitizeForSQLAndAppend(s.get_char(i));
  }
  return *this;
}

// MDS 20140722: modified this quite a bit, previously the switch
// would fall through to the generic sanitize. Not sure if that was intended.
// If so, call sanitize and append instead of pushChar as default
_StringBuffer& _StringBuffer::SanitizeForPostScriptAndAppend(const char c) {
  switch (c) {
    case '(':
    case ')':
    case '%':
      this->PushChar('\\');
      this->PushChar(c);
      break;
    default:
      this->PushChar(c);
  }
  return *this;
}

_StringBuffer& _StringBuffer::SanitizeForPostScriptAndAppend(const _String& s) {
  unsigned long sl = s.length ();
  for (unsigned long i = 0UL; i < sl; i++) {
    this->SanitizeForPostScriptAndAppend(s.get_char(i));
  }
}

_StringBuffer& _StringBuffer::SanitizeForHTMLAndAppend(const char c) {
  switch (c) {
    case '"':
      (*this) << "&quot;";
      break;
    case '\'':
      (*this) << "&apos;";
      break;
    case '<':
      (*this) << "&lt;";
      break;
    case '>':
      (*this) << "&gt;";
      break;
    case '&':
      (*this) << "&amp;";
      break;
    default:
      this->PushChar(c);
  }
  return *this;
}

_StringBuffer& _StringBuffer::SanitizeForHTMLAndAppend(const _String& s) {
  unsigned long sl = s.length ();
  for (unsigned long i = 0UL; i < sl; i++) {
    this->SanitizeForHTMLAndAppend(s.get_char(i));
  }
  return *this;
}

_StringBuffer& _StringBuffer::SanitizeForRegExAndAppend(const char c) {
  switch (c) {
    case '[':
    case '^':
    case '$':
    case '.':
    case '|':
    case '?':
    case '*':
    case '+':
    case '(':
    case ')':
      this->PushChar('\\');
      this->PushChar(c);
      break;
    case '\\':
      (*this) << "\\\\";
      break;
    default:
      this->PushChar(c);
  }
  return *this;
}

_StringBuffer& _StringBuffer::SanitizeForRegExAndAppend(const _String& s) {
  unsigned long sl = s.length ();
  for (unsigned long i = 0UL; i < sl; i++) {
    this->SanitizeForHTMLAndAppend(s.get_char(i));
  }
  return *this;
}

_StringBuffer& _StringBuffer::SanitizeAndAppend(const char c) {
  switch (c) {
    case '\n':
    case '\t':
    case '"':
    case '\\':
      this->PushChar('\\');
      this->PushChar(c);
      break;
    default:
      this->PushChar(c);
  }
  return *this;
}

_StringBuffer& _StringBuffer::SanitizeAndAppend(const _String& s) {
  unsigned long sl = s.length ();
  for (unsigned long i = 0UL; i < sl; i++) {
    this->SanitizeAndAppend(s.get_char(i));
  }
  return *this;
}



void _StringBuffer::AppendAnAssignmentToBuffer(_String* id, _String *value, unsigned long flags) {
  
  if (flags & kAppendAnAssignmentToBufferGlobal) {
    (*this) << "global ";
  }
  
  (*this) << id;
  
  if (flags & kAppendAnAssignmentToBufferAssignment) {
    (*this) << ':';
  }
  
  (*this) << '=';
  
  if (flags & kAppendAnAssignmentToBufferQuote) {
    (*this) << '"';
  }
  
  (*this) << value;
  
  if (flags & kAppendAnAssignmentToBufferQuote) {
    (*this) << '"';
  }
  
  (*this) << ";\n";
  
  if (flags & kAppendAnAssignmentToBufferFree) {
    DeleteObject (value);
  }
}

void _StringBuffer::AppendVariableValueAVL (_String const* id, _SimpleList const& var_numbers) {
  
  for (unsigned long k=0UL; k<var_numbers.countitems(); k++) {
    _Variable *tiv = LocateVar(var_numbers.lData[k]);
    if (tiv) {
      (*this) << id
              << "[\""
              << tiv->GetName()
              << "\"]=";
      
      _PMathObj varValue = tiv->Compute();
      switch (varValue->ObjectClass()) {
        case NUMBER:
          (*this) << _String (varValue->Value());
          break;
        case STRING:
          (*this) << '"';
          SanitizeAndAppend (((_FString*)varValue)->get_str());
          (*this) << '"';
          break;
        default:
          AppendNewInstance ((_String*)(varValue->toStr()));
          break;
          
      }
      (*this) << ";\n";
    }
  }
}



