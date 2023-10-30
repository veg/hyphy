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
#include <utility>  // for std::move


#ifndef _USE_EMSCRIPTEN_
    unsigned char               _StringBuffer::preallocated_buffer [_HY_STRING_BUFFER_PREALLOCATE_SLOTS*sizeof (_StringBuffer)];
    _SimpleList                 _StringBuffer::free_slots;
#endif

_List      _hyTerminalColors_aux;

_AVLListXL _hyTerminalColors (&_hyTerminalColors_aux, _List (
    new _String("RED"), new _String ("\033[31m"),
    new _String("GREEN"), new _String ("\033[32m"),
    new _String("NONE"), new _String ("\033[0m")
));

/*
==============================================================
Utility Functions
==============================================================
*/

void _StringBuffer::AllocateBufferSpace(const unsigned long character_count) {
  sa_length = character_count;
  if (s_data) {
    if (s_data  != allocated_ptr) {
        char * t =  (char *)MemAllocate(sa_length + 1UL);
        memcpy (t, s_data, s_length );
        free (allocated_ptr);
        allocated_ptr = s_data = t;
    } else {
        allocated_ptr = s_data = (char *)MemReallocate(allocated_ptr, sa_length + 1UL);
    }
  } else {
    allocated_ptr = s_data = (char *)MemAllocate(sa_length + 1UL);
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
  allocated_ptr = s_data;
}

//=============================================================

void _StringBuffer::TrimSpace (void) {
    if (sa_length > s_length) {
        AllocateBufferSpace (s_length);
    }
}

//=============================================================

void _StringBuffer::Clear() {
  s_data = allocated_ptr;
  _String::Clear();
  sa_length = 0UL;
  allocated_ptr = s_data;
}

//=============================================================

void _StringBuffer::Reset() {
    s_length = 0L;
    if (s_data) {
        s_data[0] = '\0';
    }
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
_StringBuffer::_StringBuffer(const _StringBuffer &source, long start, long end) {
    if (source.s_length) {
        
        long requested_range = source.NormalizeRange(start, end);
        
        if (requested_range > 0L) {
            AllocateBufferSpace (requested_range);
            s_length = requested_range;
            memcpy (s_data, source.s_data + start, sizeof(char)*requested_range);
            s_data[requested_range] = '\0';
            return;
            
        }
    }
    
    s_length = 0UL;
    sa_length = 0UL;
    allocated_ptr = s_data = (char *)MemAllocate(1UL);
    s_data[0] = '\0';
    
}

//=============================================================

_StringBuffer::_StringBuffer(_String&& buffer) : _String (std::move (buffer)) {
    sa_length = s_length;
    allocated_ptr = s_data;
}

//=============================================================
_StringBuffer::_StringBuffer(const _StringBuffer &s): _String () {
  this->Initialize();
  this->Duplicate (&s);
}

//=============================================================

_StringBuffer::_StringBuffer(_String* buffer): _String (buffer) {
  sa_length = s_length;
  allocated_ptr = s_data;
}

//=============================================================

_StringBuffer::~_StringBuffer (void ){
    if (s_data && allocated_ptr)
        s_data = allocated_ptr;
    sa_length = 0L;
}

//=============================================================
_StringBuffer::_StringBuffer(hyFile* file): _String () {
    const unsigned long buffer_size = 65535;
    this->Initialize();
    char buffer [buffer_size+1L];
    unsigned long items_read;
    _String buffer_str (buffer_size, buffer);
    do {
        items_read = file->read (buffer, 1, buffer_size);
        if (items_read < buffer_size) break;
        (*this) << buffer_str;
    } while (items_read == buffer_size);
    if (items_read) {
        buffer[items_read] = 0;
        (*this) << buffer;
    }
}

/*
==============================================================
Cloners and Copiers
==============================================================
*/

void _StringBuffer::Duplicate (BaseRefConst ref) {
  if (s_data && s_data != allocated_ptr) {
    s_data = allocated_ptr;
  }
  _String::Duplicate (ref);
  this->sa_length = this->s_length;
  this->allocated_ptr = this->s_data;
}

//=============================================================


BaseRef _StringBuffer::makeDynamic (void) const {
  _StringBuffer * r = new _StringBuffer;
  r->Duplicate(this);
  return r;
}

//=============================================================


_StringBuffer& _StringBuffer::operator = (_StringBuffer && rhs) {
    if (&rhs != this) {
        Clear();
        s_data = rhs.s_data;
        s_length = rhs.s_length;
        sa_length = rhs.sa_length;
        allocated_ptr = rhs.allocated_ptr;
        rhs.allocated_ptr = nil;
        rhs._String::Initialize();
    }
    return *this;
}



//=============================================================


_StringBuffer& _StringBuffer::operator = (_StringBuffer const & rhs) {
    if (&rhs != this) {
        Clear();
        Initialize();
        AllocateBufferSpace(rhs.s_length);
        (*this) << rhs;
    }
    return *this;
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
 long requested_range = source.NormalizeRange(start, end);
 
  if (requested_range > 0L) {
      for (long i = start; i <= end; i++) {
          (*this) << source.char_at(i);
      }
  }
  return (*this);
}

//=============================================================

void _StringBuffer::Trim(long start, long end) {
    
    
    long resulting_length = NormalizeRange(start, end);

    if (resulting_length > 0L) {
        if (start > 0L) {
             s_data    += start;
             sa_length -= start;
            //memmove(s_data, s_data + start, resulting_length);
        }
        if (s_length != resulting_length) {
            s_length = resulting_length;
            s_data[resulting_length] = '\0';
        }
    } else {
        s_length = 0UL;
        s_data[0] = '\0';
    }
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
    this->SanitizeForSQLAndAppend(s.char_at(i));
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
    this->SanitizeForPostScriptAndAppend(s.char_at(i));
  }
  return *this;
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
    this->SanitizeForHTMLAndAppend(s.char_at(i));
  }
  return *this;
}

_StringBuffer& _StringBuffer::ConvertToTerminalColor(const _String& s) {
#ifdef _USE_EMSCRIPTEN_
    static const _String kNONE ("NONE");
#endif
    BaseRef color_lookup = _hyTerminalColors.GetDataByKey (&s);
    if (color_lookup) {
#ifdef _USE_EMSCRIPTEN_
        if (s == kNONE) {
            (*this) << "</span>";
        } else {
            (*this) << "<span style='color:" << s << "'>";
        }
#else
        (*this) << *(const _String*)color_lookup;
#endif
    } else {
        (*this) << s;
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
    this->SanitizeForRegExAndAppend(s.char_at(i));
  }
  return *this;
}

_StringBuffer& _StringBuffer::SanitizeAndAppend(const char c) {
  switch (c) {
    case '\n':
          this->PushChar ('\\');
          this->PushChar ('n');
          break;
    case '\t':
        this->PushChar ('\\');
        this->PushChar ('t');
        break;

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
    this->SanitizeAndAppend(s.char_at(i));
  }
  return *this;
}



void _StringBuffer::AppendAnAssignmentToBuffer(_String const* id, _String *value, unsigned long flags) {
  
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
    _Variable *tiv = LocateVar(var_numbers.list_data[k]);
    if (tiv) {
      (*this) << id
              << "[\""
              << tiv->GetName()
              << "\"]=";
      
      HBLObjectRef varValue = tiv->Compute();
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

//__________________________________________________________________________________
void * _StringBuffer::operator new (size_t size) {
#ifndef _USE_EMSCRIPTEN_
    if (_StringBuffer::free_slots.nonempty()) {
        _StringBuffer * result = (_StringBuffer *)(_StringBuffer::preallocated_buffer) + _StringBuffer::free_slots.Pop();
        //printf ("Allocate slot %ld\n", _StringBuffer::free_slots.get (_StringBuffer::free_slots.countitems()));
        //result->Initialize();
        return result;
    }
#endif
    //printf ("Allocating string buffer %ld\n", size);
    return ::operator new(size);
}


//__________________________________________________________________________________
void  _StringBuffer::operator delete (void * p) {
#ifndef _USE_EMSCRIPTEN_
    _StringBuffer * sb = (_StringBuffer*) p,
                  * sbb = (_StringBuffer *)_StringBuffer::preallocated_buffer;
    
    if (sb >= sbb &&  sb < (sbb+ _HY_STRING_BUFFER_PREALLOCATE_SLOTS)) {
        //printf ("Free slot %ld\n", ((long)(sb - _StringBuffer::preallocated_buffer)));
        //sb->_StringBuffer::Clear();
        free_slots << ((long)(sb - sbb));
    } else {
#endif
        //printf ("%p\n", p);
        ::operator delete (p);
#ifndef _USE_EMSCRIPTEN_
    }
#endif
}


