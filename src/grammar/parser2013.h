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

void _parser2013_pushNumber         (void * vp, _Formula& f, _FormulaParsingContext& fpc, const wchar_t* value);
    // convert a floating point number stored in value and push it onto 'f' 

void _parser2013_pushString         (void * vp, _Formula& f, _FormulaParsingContext& fpc, const wchar_t* value);
    // process a string stored in value and push it onto 'f'

void _parser2013_pushObject        (void * vp, _Formula& f, _FormulaParsingContext& fpc, _PMathObj);
    // push an object onto the formula

void _parser2013_pushNone        (void * vp, _Formula& f, _FormulaParsingContext& fpc);
// process None onto 'f'

void _parser2013_pushIdentifier     (void * vp, _Formula& f, _FormulaParsingContext& fpc, const wchar_t* value, bool globalKey, bool takeVarReference);
// process a variable identifier (take __references & */^ references if needed) and push onto 'f'

void _parser2013_pushOp           (void * vp, _Formula& f, _FormulaParsingContext& fpc, long op_code, long num_terms);
// push an operation given it's internal op_code and the number of terms

void _parser2013_pushFunctionCall           (void * vp, _Formula& f, _FormulaParsingContext& fpc, _String& fundID, const _List&);
// push a function call given the function id and the list of 'named' arguments
// a blank entry in the list implies that the argument is positional
// generally speaking, all positional arguments need to preceed all named arguments

void _parser2013_matrix_checkRowLengths (void *vp, _FormulaParsingContext& fpc, unsigned long & global_count, unsigned long& local_count);

long _parser2013_checkLvalue (void *vp, _Formula &f, _FormulaParsingContext& fpc);

void _parser2013_pushSparseElementEntry (void *vp, _FormulaParsingContext& fpc, _SimpleList&, _Formula*, _Formula*, _Formula*, bool & );

void _parser2013_createSparseMatrix (void* vp, _Formula&, _FormulaParsingContext&, 
          _Formula*, _Formula*, _SimpleList*, bool); 
          



_Matrix*  _parser2013_createDenseMatrix (void* vp, _FormulaParsingContext& fpc, _SimpleList* entries, 
      const unsigned long n_rows, const unsigned long n_cols, const bool is_const); 
/*

  Create a _Matrix object from the list of entries, stored row by row
  elements row by row
    
  is_const is set to true is each entry of the matrix is a number
  
  
*/

void  _parser2013_handleAssignment (void* vp, _Formula& lhs, _Formula &rhs, 
                                        _FormulaParsingContext& fpc, long assignment_type,
                                        long op_code, long lvalue_index);


void _parser2013_addADictionaryElement (void* vp, _SimpleList& dictionary_entries, _Formula* key, _Formula *value, 
                                        _FormulaParsingContext& fpc, bool & is_const);
                                        
void _parser2013_createDictionary (void* vp, _Formula &f, _FormulaParsingContext& fpc, 
                                  _SimpleList& dictionary_entries, bool is_const);
                                        

void _parser2013_add_matrix_entry (void *vp, _SimpleList& matrix_entries, _Formula* f, _FormulaParsingContext& fpc, bool & is_const);

void _parser2013_pushStatementOntoList (void *vp, _ExecutionList& current_command_stream, _Formula* f);
void _parser2013_pushJumpOntoList (void *vp, _ExecutionList& current_command_stream, _Formula* f);
void _parser2013_pushSetJumpCommmandIndices (void *vp, _ExecutionList&, long, long);
void  _parser2013_addLoopContext (void *vp);
void  _parser2013_popLoopContext (void *vp, _ExecutionList&, long, long);
void  _parser2013_handleContinueBreak (void *vp, _ExecutionList&, bool);

// grammar conflict resolvers

bool    _parser2013_IdentFollowedByAnOpenParenthesis (void * p);
bool    _parser2013_TwoOpenBraces                   (void * p);
bool    _parser2013_isSimpleStatement             (void * p);
bool    _parser2013_isFollowedByAnCommaOrClosingBrace (void *p);
bool    _parser2013_StringAndColon (void *p);


// Utility functions

void    _parser2013_reportError                   (void * vp, const _String, _FormulaParsingContext&);
bool    _parser2013_errorFree                     (void * vp);
long    _parser2013_expected_arguments            (const long);
#endif
