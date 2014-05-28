/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2009
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon              (apoon@cfenet.ubc.ca)

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

  // utility functions
void _StringBuffer::AllocateBufferSpace (const unsigned long character_count) {
  saLength = character_count;
  if (sData) {
    checkPointer(sData = (char *)MemReallocate(sData, saLength + 1UL));
  } else {
    checkPointer(sData = (char *)MemAllocate(saLength + 1UL));
  }
  
}

void _StringBuffer::ResizeString (void) {
  if (saLength < sLength) {
    unsigned long incBy = sLength - saLength;
    
    if (incBy < HY_STRING_BUFFER_ALLOCATION_CHUNK) {
      incBy = HY_STRING_BUFFER_ALLOCATION_CHUNK;
    }
    
    if (incBy < sLength >> 4) {
      incBy = sLength >> 4;
    }
    
    AllocateBufferSpace (saLength + incBy);
  }
}

void _StringBuffer::Initialize (bool p) {
  _String::Initialize(p);
  saLength = 0UL;
}


  // constructors

_StringBuffer::_StringBuffer (void) : _String () {
  AllocateBufferSpace (HY_STRING_BUFFER_ALLOCATION_CHUNK);
  sData[0] = 0;
}

_StringBuffer::_StringBuffer (const unsigned long character_count) : _String () {
  AllocateBufferSpace (character_count);
  sData[0] = 0;
}

_StringBuffer::_StringBuffer (const char * buffer) {
  Initialize ();
  (*this) << buffer;
}

_StringBuffer::_StringBuffer (const _String& buffer) {
  Initialize ();
  (*this) << buffer;
}


//Stack copy contructor
_StringBuffer::_StringBuffer(const _StringBuffer &s) { Duplicate(& s); }

/*
 ==============================================================
cloners/copiers
 ==============================================================
 */

void _StringBuffer::Duplicate (BaseRefConst src_obj) {
  _String::Duplicate (src_obj);
  saLength = ((_StringBuffer*)src_obj)->saLength;
}

BaseRef _StringBuffer::makeDynamic (void) const {
  return new _StringBuffer (*this);
}

/*
==============================================================
Operator Overloads
==============================================================
*/


void _StringBuffer::operator<<(const _String *s) {
  if (s && s->sLength) {
    PushCharBuffer (s->sData, s->sLength);
  }
}

// append operator
void _StringBuffer::operator<<(const _String &s) { (*this) << &s; }

//Append operator
void _StringBuffer::operator<<(const char *str) {
  PushCharBuffer(str, strlen (str));
}

//Append operator
void _StringBuffer::operator<<(const char c) {
  PushChar (c);
}
/*
==============================================================
Methods
==============================================================
*/

void _StringBuffer::PushChar(const char c){
  sLength ++;
  ResizeString ();
  sData[sLength-1UL] = c;
  sData[sLength] = 0;
}

void _StringBuffer::PushCharBuffer(const char* buffer, const unsigned long buffer_l){
  unsigned long offset = sLength;
  sLength += buffer_l;
  ResizeString ();
  
  for (unsigned long k = 0UL; k < buffer_l; k++) {
    sData[offset + k] = buffer[k];
  }
  sData[sLength] = 0;
}

//Append and delete operator
void _StringBuffer::AppendNewInstance(_String *s) {
  (*this) << s;
  DeleteObject(s);
}

  //Append operator
void _StringBuffer::EscapeAndAppend(const char c, const _hyStringBufferEscapeMode mode) {
  if (mode == HY_ESCAPE_SQLITE) {
    PushChar(c);
    switch (c) {
      case '\'':
        PushChar(c);
    }
    return;
  } else {
    if (mode == HY_ESCAPE_POSTSCRIPT) {
      switch (c) {
        case '(':
        case ')':
        case '%':
          PushChar ('\\');
          PushChar(c);
          return;
      }
    } else {
      if (mode == HY_ESCAPE_HTML) {
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
            PushChar(c);
        }
        return;
      } else {
        if (mode == HY_ESCAPE_REGEXP) { // regexp
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
              PushChar('\\');
              PushChar(c);
              break;
            case '\\':
              (*this) << "\\\\";
              break;
            default:
              PushChar(c);
          }
          return;
          
        }
      }
    }
  }
  switch (c) {
    case '\n':
      PushChar('\\');
      PushChar('n');
      break;
    case '\t':
      PushChar('\\');
      PushChar('t');
      break;
    case '"':
      PushChar('\\');
      PushChar('"');
      break;
    case '\\':
      PushChar('\\');
      PushChar('\\');
      break;
    default:
      PushChar(c);
  }
}

  //Append operator
void _StringBuffer::EscapeAndAppend(const _String &s, const _hyStringBufferEscapeMode mode) {
  for (unsigned long i = 0UL; i < s.sLength; i++) {
    EscapeAndAppend(s.sData[i], mode);
  }
}

void _StringBuffer::AppendAnAssignmentToBuffer(_String *id, _String *value,
                                         bool doFree, bool doQuotes,
                                         bool doBind) {
  (*this) << id;
  if (doBind) {
    (*this) << ':';
  }
  (*this) << '=';
  if (doQuotes) {
    (*this) << '"';
  }
  (*this) << value;
  if (doQuotes) {
    (*this) << '"';
  }
  (*this) << ";\n";
  if (doFree) {
    DeleteObject(value);
  }
}


void _StringBuffer::AppendVariableValueAVL(_String *id, _SimpleList &varNumbers) {
#ifndef HY_2014_REWRITE_MASK
  for (long k = 0; k < varNumbers.lLength; k++) {
    _Variable *tiv = LocateVar(varNumbers.lData[k]);
    if (tiv) {
      (*this) << id;
      (*this) << "[\"";
      (*this) << tiv->GetName();
      (*this) << "\"]=";
      _PMathObj varValue = tiv->Compute();
      switch (varValue->ObjectClass()) {
        case NUMBER:
          (*this) << _String(varValue->Value());
          break;
        case STRING:
          (*this) << '"';
          EscapeAndAppend(*((_FString *)varValue)->theString);
          (*this) << '"';
          break;
        default:
          AppendNewInstance((_String *)(varValue->toStr()));
          break;
          
      }
      (*this) << ";\n";
    }
  }
#endif
}


