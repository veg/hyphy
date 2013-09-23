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

#ifndef __HY_HBL_FUNCTION__
#define __HY_HBL_FUNCTION__

#include "executionlist.h"
#include "trie.h"
#include "_hyExecutionContext.h"


/* the HBL function object
  the function object is simply an execution list (base object), 
  together with a list of parameters, and their associated default values.
 
  The _Variable base class stores the function ID (in the name variable, i.e.
  in the global namespace)
 
  The 'parameters' field stores parameter names along with their positions 
  in the function declaration
  
  The 'default_values' field stores default parameter values (or None objects)
 
 
*/
 
//____________________________________________________________________________________

class _HBLFunction : public virtual _MathObject, _ExecutionList  {
  
  _HBLFunction          (_Trie& names, _List& values);
  virtual               ~_HBLFunction (void);
  
  _PMathObj             Call (const _Trie*, const _List *, _hyExecutionContext* = _hyDefaultExecutionContext);
  virtual               unsigned long        ObjectClass (void) {return HBL_FUNCTION; }
  virtual               void Duplicate (BaseRef);
  virtual               BaseRef toStr();
  virtual               void    toFileStr(FILE*);
  virtual               BaseRef makeDynamic ();
  
  
private:
  
  _Trie parameters;
  _List default_values;
  
};

_String* _hyMapArgumentValues (
                           const _Trie* named_args,
                           const _Trie* func_args,
                           const _List* passed_values,
                           const _List* default_values,
                           _List & mapped_call
                           );

#endif
