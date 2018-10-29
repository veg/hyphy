/*

 HyPhy - Hypothesis Testing Using Phylogenies.

 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon42@uwo.ca)
 Steven Weaver (sweaver@temple.edu)

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

#include "formula_parsing_context.h"
#include "variablecontainer.h"

//__________________________________________________________________________________
_FormulaParsingContext::_FormulaParsingContext (_String* err, _VariableContainer const* scope) {
    assignment_ref_id   = -1;
    assignment_ref_type = kStringDirectReference;
    is_volatile = false;
    in_assignment = false;
    build_complex_objects = true;
    allow_template = '\0';
    err_msg = err;
    formula_scope = scope;
}

//__________________________________________________________________________________
_String const _FormulaParsingContext::contextualizeRef (_String& ref) {
    if (formula_scope) {
        return *formula_scope->GetName () & '.' & ref;
    }
    return ref;
}

//__________________________________________________________________________________
void _FormulaParsingContext::setScope(const _String *scope) {
  if (scope && scope->nonempty()) {
    _VariableContainer vc (*scope);
    formula_scope = (_VariableContainer*)FetchVar(vc.get_index());
  } else {
    formula_scope = nil;
  }
}
