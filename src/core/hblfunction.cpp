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

#include "hblfunction.h"
#include "trieiterator.h"

_HBLFunction :: _HBLFunction (_Trie& names, _List& values) : _ExecutionList () {
  if (names.countitems() != values.countitems()) {
    WarnError ("Internal Error in _HBLFunction :: _HBLFunction; the lengths of argument arrays must be the same");
  }
  parameters.Duplicate    (&names);
  default_values.Duplicate(&default_values);
}


_HBLFunction::~_HBLFunction (void) {
}

_PMathObj _HBLFunction::Call (const _Trie* names_arguments, const _List* values, _hyExecutionContext * parent_context) {
  return NULL;
}

_String* _hyMapArgumentValues (const _Trie* named_args,
                               const _Trie* func_args,
                               const _List* passed_values,
                               const _List* default_values,
                               _List & mapped_call
                           ) {
  
  
    /* call logic:
        1). traverse all the named_args, and assign supplied values
            to corresponding function parameters;
              throws an error if 
                (a) there is no function argument by such name
        2). take all remaining (positional) values and assign them to 
              functional arguments
        3). if any unassigned arguments remain, check to see if they 
              have default values; if not throw an error
     
     return NULL if all arguments got mapped, otherwise the error message
     
    */
  
  mapped_call.Clear(false);
  
  _SimpleList argument_disposition (func_args->countitems(),HY_NOT_FOUND, 0L);
  
  /* this list will have values >=0L mapped to passed arguments and 
     < -1L for default arguments. If any of the values are left are HY_NOT_FOUND, 
     this means the mapping has failed */
  
  if (named_args && named_args->countitems()) {
    _TrieIterator ti (*named_args);
    _String       *a_name = ti.First();
    while (a_name) {
      long arg_id = func_args->Find (*a_name);
      if (arg_id == HY_TRIE_NOTFOUND) {
        *a_name = _String ("There is no function argument named '") & *a_name & "'";
        return a_name;
      }
      argument_disposition.lData[func_args->GetValue(arg_id)] = named_args->GetValue (ti.CurrentIndex());
      a_name = ti.Next();
    }
    DeleteObject (a_name);
  }
  
  for (long k = 0L; k < argument_disposition.lLength; k+=1) {
    if (argument_disposition.lData[k] == HY_NOT_FOUND) {
      if (passed_values->lLength <= k) {
        if ((BaseRef)default_values->GetElement(k) != NULL) {
          argument_disposition.lData[k] = -k-2;
        } else {
          return new _String(_String ("No value supplied for a positional argument with no default value (index ") & k & ')');
        }
      }
    }
  }
  
  return NULL;
}


