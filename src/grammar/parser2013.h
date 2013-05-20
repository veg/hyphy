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

#ifndef __PARSER2013_
#define __PARSER2013__

#include "legacy_parser.h"
#include "batchlan.h"
#include "executionlist.h"
#include "formula.h"
#include "wchar.h"

// parser support functions

void _parser2013_pushNumber         (_Formula& f, _FormulaParsingContext& fpc, const wchar_t* value);
    // convert a floating point number stored in value and push it onto 'f' 

void _parser2013_pushString         (_Formula& f, _FormulaParsingContext& fpc, const wchar_t* value);
    // process a string stored in value and push it onto 'f'

void _parser2013_pushNone        (_Formula& f, _FormulaParsingContext& fpc);
// process None onto 'f'

void _parser2013_pushIdentifier     (_Formula& f, _FormulaParsingContext& fpc, const wchar_t* value, bool globalKey, bool takeVarReference);
// process a variable identifier (take __references & */^ references if needed) and push onto 'f'

void _parser2013_pushOp           (_Formula& f, _FormulaParsingContext& fpc, long op_code, long num_terms);
// push an operation given it's internal op_code and the number of terms

void _parser2013_pushFunctionCall           (_Formula& f, _FormulaParsingContext& fpc, _String& fundID, const _List&);
// push a function call given the function id and the list of 'named' arguments
// a blank entry in the list implies that the argument is positional
// generally speaking, all positional arguments need to preceed all named arguments

// grammar conflict resolvers

bool    IsFollowedByAnOpenParenthesis (void * p);
bool    IsSimpleStatement             (void * p);

#endif
