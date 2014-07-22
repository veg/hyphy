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

#include "hy_string_buffer.h"
#include "helperfunctions.h"
#include "errorfns.h"

#include <string.h> // for strlen


/*
==============================================================
Utility Functions
==============================================================
*/

void _StringBuffer::allocateBufferSpace(const unsigned long character_count) {
  sa_length = character_count;
  if (s_data) {
    checkPointer(s_data = (char *)MemReallocate(s_data, sa_length + 1UL));
  } else {
    checkPointer(s_data = (char *)MemAllocate(sa_length + 1UL));
  }

}

void _StringBuffer::resizeString(void) {
  if (sa_length < s_length) {
    unsigned long inc_by = s_length - sa_length;

    if (inc_by < HY_STRING_BUFFER_ALLOCATION_CHUNK) {
      inc_by = HY_STRING_BUFFER_ALLOCATION_CHUNK;
    }

    if (inc_by < s_length >> 4) {
      inc_by = s_length >> 4;
    }

    this->allocateBufferSpace(sa_length + inc_by);
  }
}

void _StringBuffer::initialize(bool p) {
  _String::Initialize(p);
  sa_length = 0UL;
}


/*
==============================================================
Constructors
==============================================================
*/

_StringBuffer::_StringBuffer(void) : _String () {
  this->allocateBufferSpace(HY_STRING_BUFFER_ALLOCATION_CHUNK);
  s_data[0] = 0;
}

_StringBuffer::_StringBuffer(const unsigned long character_count)
                            : _String () {
  this->allocateBufferSpace(character_count);
  s_data[0] = 0;
}

_StringBuffer::_StringBuffer(const char * buffer) {
  this->initialize();
  (*this) << buffer;
}

_StringBuffer::_StringBuffer(const _String& buffer) {
  this->initialize();
  (*this) << buffer;
}

// Stack copy contructor
_StringBuffer::_StringBuffer(const _StringBuffer &s) {
  this->duplicate(&s);
}

/*
==============================================================
Cloners and Copiers
==============================================================
*/

void _StringBuffer::duplicate(BaseRefConst src_obj) {
  _String::Duplicate(src_obj);
  sa_length = ((_StringBuffer*)src_obj)->sa_length;
}

BaseRef _StringBuffer::makeDynamic(void) const {
  return new _StringBuffer(*this);
}

/*
==============================================================
Operator Overloads
==============================================================
*/


// Append operator
void _StringBuffer::operator<<(const _String *s) {
  if (s && s->s_length) {
    this->pushCharBuffer(s->s_data, s->s_length);
  }
}

// Append operator
void _StringBuffer::operator<<(const _String &s) {
  (*this) << &s;
}

// Append operator
void _StringBuffer::operator<<(const char *str) {
  this->pushCharBuffer(str, strlen(str));
}

// Append operator
void _StringBuffer::operator<<(const char c) {
  this->pushChar(c);
}

/*
==============================================================
Methods
==============================================================
*/

void _StringBuffer::pushChar(const char c){
  s_length++;
  this->resizeString();
  s_data[s_length-1UL] = c;
  s_data[s_length] = 0;
}

void _StringBuffer::pushCharBuffer( const char* buffer,
                                    const unsigned long buffer_l) {
  if (buffer_l) {
    unsigned long offset = s_length;
    s_length += buffer_l;
    this->resizeString();

    for (unsigned long k = 0UL; k < buffer_l; k++) {
      s_data[offset + k] = buffer[k];
    }
    s_data[s_length] = 0;
  }
}

// Append and delete
void _StringBuffer::appendNewInstance(_String *s) {
  (*this) << s;
  DeleteObject(s);
}

void _StringBuffer::sanitizeForSQLAndAppend(const char c) {
  this->pushChar(c);
  if (c == '\'') {
    this->pushChar(c);
  }
}

void _StringBuffer::sanitizeForSQLAndAppend(const _String &s) {
  for (unsigned long i = 0UL; i < s.s_length; i++) {
    this->sanitizeForSQLAndAppend(s.s_data[i]);
  }
}

// MDS 20140722: modified this quite a bit, previously the switch
// would fall through to the generic sanitize. Not sure if that was intended.
// If so, call sanitize and append instead of pushChar as default
void _StringBuffer::sanitizeForPostScriptAndAppend(const char c) {
  switch (c) {
    case '(':
    case ')':
    case '%':
      this->pushChar('\\');
      this->pushChar(c);
      break;
    default:
      this->pushChar(c);
  }
}

void _StringBuffer::sanitizeForPostScriptAndAppend(const _String &s) {
  for (unsigned long i = 0UL; i < s.s_length; i++) {
    this->sanitizeForPostScriptAndAppend(s.s_data[i]);
  }
}

void _StringBuffer::sanitizeForHTMLAndAppend(const char c) {
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
      this->pushChar(c);
  }
}

void _StringBuffer::sanitizeForHTMLAndAppend(const _String &s) {
  for (unsigned long i = 0UL; i < s.s_length; i++) {
    this->sanitizeForHTMLAndAppend(s.s_data[i]);
  }
}

void _StringBuffer::sanitizeForRegExAndAppend(const char c) {
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
      this->pushChar('\\');
      this->pushChar(c);
      break;
    case '\\':
      (*this) << "\\\\";
      break;
    default:
      this->pushChar(c);
  }
}

void _StringBuffer::sanitizeForRegExAndAppend(const _String &s) {
  for (unsigned long i = 0UL; i < s.s_length; i++) {
    this->sanitizeForRegExAndAppend(s.s_data[i]);
  }
}

void _StringBuffer::sanitizeAndAppend(const char c) {
  switch (c) {
    case '\n':
      this->pushChar('\\');
      this->pushChar('n');
      break;
    case '\t':
      this->pushChar('\\');
      this->pushChar('t');
      break;
    case '"':
      this->pushChar('\\');
      this->pushChar('"');
      break;
    case '\\':
      this->pushChar('\\');
      this->pushChar('\\');
      break;
    default:
      this->pushChar(c);
  }
}

void _StringBuffer::sanitizeAndAppend(const _String &s) {
  for (unsigned long i = 0UL; i < s.s_length; i++) {
    this->sanitizeAndAppend(s.s_data[i]);
  }
}

// Append and delete
void _StringBuffer::appendSubstring(const _String& s, long from, long to) {
  (*this) << _String (s, from, to);
}

// Special purpose append
void _StringBuffer::appendAnAssignmentToBuffer( _String *id,
                                                _String *value,
                                                bool do_free,
                                                bool do_quotes,
                                                bool do_bind) {
  (*this) << id;
  if (do_bind) {
    (*this) << ':';
  }
  (*this) << '=';
  if (do_quotes) {
    (*this) << '"';
  }
  (*this) << value;
  if (do_quotes) {
    (*this) << '"';
  }
  (*this) << ";\n";
  if (do_free) {
    DeleteObject(value);
  }
}

// Special purpose append
void _StringBuffer::appendVariableValueAVL( _String *id,
                                            _List &var_numbers) {
#ifndef HY_2014_REWRITE_MASK
  for (long k = 0; k < var_numbers.l_length; k++) {
    _Variable *tiv = LocateVar(var_numbers.l_data[k]);
    if (tiv) {
      (*this) << id;
      (*this) << "[\"";
      (*this) << tiv->GetName();
      (*this) << "\"]=";
      _PMathObj var_value = tiv->Compute();
      switch (var_value->ObjectClass()) {
        case NUMBER:
          (*this) << _String(var_value->Value());
          break;
        case STRING:
          (*this) << '"';
          this->EscapeAndAppend(*((_FString *)var_value)->theString);
          (*this) << '"';
          break;
        default:
          this->appendNewInstance((_String *)(var_value->toStr()));
          break;

      }
      (*this) << ";\n";
    }
  }
#endif
}


