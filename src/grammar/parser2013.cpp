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

#include "Parser.h"
#include "hy_globals.h"


  // parser support functions

void _parser2013_pushNumber (void * vp, _Formula& f, _FormulaParsingContext& fpc, const wchar_t* value) {
  if (_parser2013_errorFree(vp) == false) return;
  wchar_t * endptr = NULL;
  double numeric_value = wcstod (value, &endptr);
  f.Push (new _Operation (new _Constant (numeric_value)));
}

void _parser2013_pushString (void * vp, _Formula& f, _FormulaParsingContext& fpc, const wchar_t* value) {
  if (_parser2013_errorFree(vp) == false) return;
  _String literal (value);
  literal.Trim(1, literal.sLength - 2);
  f.Push (new _Operation (new _FString (literal)));
}

void _parser2013_pushNone (void * vp, _Formula& f, _FormulaParsingContext& fpc) {
  if (_parser2013_errorFree(vp) == false) return;
  f.Push (new _Operation ( noneToken,false));
}

void _parser2013_pushOp (void * vp, _Formula& f, _FormulaParsingContext& fpc, long op_code, long num_terms) {
  if (_parser2013_errorFree(vp) == false) return;
  f.Push (new _Operation ( op_code,num_terms));
}

void _parser2013_pushObject (void * vp, _Formula& f, _FormulaParsingContext& fpc, _PMathObj obj) {
  if (_parser2013_errorFree(vp) == false) return;
  f.Push (new _Operation (obj));
}

void _parser2013_pushFunctionCall (void * vp, _Formula& f, _FormulaParsingContext& fpc, _String& funcId, const _List& argumentNames) {
  if (_parser2013_errorFree(vp) == false) return;
  long arg_count   = argumentNames.countitems(),
  built_in_id = BuiltInFunctions.GetValueFromString(funcId);
  
  if (built_in_id >= 0) {
    long expect_arg = _parser2013_expected_arguments (built_in_id);
    if (arg_count != expect_arg) {
      _parser2013_reportError(vp, _String(arg_count) &" arguments passed to '" & funcId & "', expected " & _String(expect_arg) & '.', fpc);
    }
    f.Push (new _Operation ( built_in_id,arg_count));
  }
  else {
    f.Push (new _Operation ( true, funcId,arg_count));
  }
}

void _parser2013_pushIdentifier (void* vp, _Formula& f, _FormulaParsingContext& fpc, const wchar_t* value, bool globalKey, bool takeVarReference) {
  if (_parser2013_errorFree(vp) == false) return;
  _String ident (value);
  long curOpl = ident.sLength;
  if (curOpl > 2 && ident[curOpl - 1] == '_' &&
      ident[curOpl - 2] == '_') { // instant variable refrence
    _String realVarName(ident, 0, curOpl - 3);
    
    realVarName = fpc.contextualizeRef(realVarName);
    
    long realVarLoc = LocateVarByName(realVarName);
    if (realVarLoc < 0) { // bad instant variable reference
      _parser2013_reportError(vp, "Attempted to take value of undeclared variable ", fpc);
    }
    _Operation * theVar = new _Operation( realVarName, 
                                         true, globalKey, fpc.formulaScope(), false, true);
    f.Push(theVar);
  } else {
    if (fpc.formulaScope() &&
        _hyApplicationGlobals.Find(&ident) >= 0) {
      f.Push(new _Operation( ident, true, globalKey, nil, takeVarReference));
    } else {
      _Operation * id_op = new _Operation( ident, true, globalKey, fpc.formulaScope(), takeVarReference);
      if (fpc.isReference()) {
        id_op->SetAttribute (_HY_OPERATION_TOGGLE);
      }
      f.Push (id_op);
    }
  }
}

_Matrix* _parser2013_createDenseMatrix (void* vp, _FormulaParsingContext& fpc, _SimpleList* entries, const unsigned long n_rows, 
                                        const unsigned long n_cols, const bool is_const) {
  
  if (_parser2013_errorFree(vp) == false) return new _Matrix;
  _Matrix * m = nil;
  
  if (n_cols * n_rows != entries->countitems()) {
    _parser2013_reportError(vp, 
                            "Internal error in _parser2013_createMatrixObject (incompatible entries and column_count arguments)",
                            fpc);
    return nil;
  }
  
  if (is_const) {
    m = new _Matrix(n_rows, n_cols, false, true, false);
  } else {
    m = new _Matrix(n_rows, n_cols, false, false, true);     
  }
  
  unsigned long overall_index = 0L;
  
  if (is_const) {
    for (unsigned long r = 0UL; r < n_rows; r ++) {
      for (unsigned long c = 0UL; c < n_cols; c++, overall_index++) {
        m->Store(r, c, ((_Formula*)entries->GetElement (overall_index))->Compute()->Value());
      }
    }
    entries->ClearFormulasInList();
  } else {
    overall_index = entries->countitems();
    _Formula ** fp = (_Formula**) m->theData;
    for (unsigned long c = 0UL; c < overall_index; c++) {
      fp[c] = (_Formula*)entries->GetElement (c);
    }
  }
  
  return m;
}

void _parser2013_createSparseMatrix (void* vp, _Formula& f, _FormulaParsingContext& fpc, 
                                     _Formula* hd, _Formula* vd, _SimpleList* entries, bool is_const) {
  
  if (_parser2013_errorFree(vp) == false) return;
  if (vd->IsEmpty()) {
    vd->Duplicate((BaseRef)hd);
  }
  
  if (is_const) {
    is_const = hd->IsConstant() && vd->IsConstant();
  }
  
  
  if (is_const) {
    _List* constants = new _List (2UL + entries->lLength);
    constants->AppendNewInstance (new _Constant (hd->Compute()->Value()));
    constants->AppendNewInstance (new _Constant (vd->Compute()->Value()));
    delete (hd); delete (vd); 
    for (unsigned long k = 0; k < entries->lLength; k++) {
      constants->AppendNewInstance (new _Constant (((_Formula*)entries->GetElement(k))->Compute()->Value()));
    }
    entries->ClearFormulasInList();
    DeleteObject(entries);
    f.Push (new _Operation (new _Matrix ((_PMathObj)constants, true)));
    DeleteObject (constants);
    
  } else {
    entries->InsertElement ((BaseRef)hd, 0L, false, false);
    entries->InsertElement ((BaseRef)vd, 1L, false, false);
    f.Push (new _Operation (_HY_OPERATION_SPARSE_MATRIX, _HY_OPERATION_INVALID_REFERENCE,
                            _HY_OPERATION_INVALID_REFERENCE, (_PMathObj) entries));
  }
  
}


void _parser2013_matrix_checkRowLengths (void *vp, _FormulaParsingContext& fpc, 
                                         unsigned long & global_count, unsigned long& local_count) {
  if (_parser2013_errorFree(vp) == false) return;
  if (global_count == 0L) {
    global_count = local_count;
  } else {
    if (global_count != local_count) {
      _parser2013_reportError(vp, _String ("Rows of unequal dimensions in a matrix constructor: ") & (long) global_count 
                              & " vs " & (long) local_count, fpc);
      
    }
  }
}

long _parser2013_checkLvalue (void *vp, _Formula &f, _FormulaParsingContext& fpc) {
  if (_parser2013_errorFree(vp) == false) return HY_NOT_FOUND;
  return f.LValueIndex (0L, true);
}

void _parser2013_add_matrix_entry (void* vp, _SimpleList& matrix_entries, _Formula* f, _FormulaParsingContext& fpc, bool & is_const) {
  if (_parser2013_errorFree(vp) == false) return;
  f->SimplifyConstants();
  if (is_const) {
    is_const = !f->IsEmpty () && f->IsConstant();
  }
  matrix_entries << (long)f;
}

void  _parser2013_handleAssignment (void* vp, _Formula& lhs, _Formula &rhs, 
                                    _FormulaParsingContext& fpc, long assignment_type,
                                    long op_code, long lvalue_index) {

  if (_parser2013_errorFree(vp) == false) return;
  if (lvalue_index == HY_NOT_FOUND) {
    _parser2013_reportError(vp, _String ("Invalid LHS in assignment"), fpc);
  } else {

    if (assignment_type == _HY_OPERATION_ASSIGNMENT_VALUE 
        || assignment_type == _HY_OPERATION_ASSIGNMENT_EXPRESSION 
        || assignment_type == _HY_OPERATION_ASSIGNMENT_BOUND 
        ) {

      long reference = lhs.PrepareLHS (lvalue_index);
      if (assignment_type != _HY_OPERATION_ASSIGNMENT_EXPRESSION) {
        lhs.GetList() << rhs.GetList(); 
        lhs.Push (new _Operation (assignment_type, reference, op_code, NULL));
      } else {
        lhs.Push (new _Operation (assignment_type, reference, HY_OP_CODE_NONE, (_PMathObj)&rhs));
          //printf ("\n%s\n", _String ((_String*)lhs.GetList().toStr()).sData);
        return;
      }
        //printf ("\n%s\n", _String ((_String*)lhs.GetList().toStr()).sData);
    } else {
      
    }
  }

  delete (&rhs);
  
}


void _parser2013_addADictionaryElement (void* vp, _SimpleList& dictionary_entries, _Formula* key, _Formula *value, _FormulaParsingContext& fpc, bool & is_const) {
  
  if (_parser2013_errorFree(vp) == false) return;
  key->SimplifyConstants();
  value->SimplifyConstants();
  
  if (is_const) {
    is_const = key->IsConstant() && value->IsConstant();
  }
  
  dictionary_entries << (long)key;
  dictionary_entries << (long)value;
}


void _parser2013_createDictionary (void* vp, _Formula &f, _FormulaParsingContext& fpc, 
                                   _SimpleList& dictionary_entries, bool is_const) {
  
  if (_parser2013_errorFree(vp) == false) return;
  
  if (is_const) {
    f.Push (new _Operation (new _AssociativeList ((_PMathObj) &dictionary_entries)));
    dictionary_entries.ClearFormulasInList();
    DeleteObject (&dictionary_entries);
  } else {
    f.Push (new _Operation (_HY_OPERATION_DICTIONARY, _HY_OPERATION_INVALID_REFERENCE,
                            _HY_OPERATION_INVALID_REFERENCE, (_PMathObj) &dictionary_entries));    
  }
  
  
}

void _parser2013_pushSparseElementEntry (void* vp, _FormulaParsingContext& fpc,
                                         _SimpleList& matrix_entries, _Formula* r, _Formula* c, _Formula* d, bool & is_const ) {
  if (_parser2013_errorFree(vp) == false) return;
  _Formula* f[3] = {r,c,d};
  for (long k = 0; k < 3; k ++) {
    f[k] ->SimplifyConstants();
    if (is_const) {
      is_const = f[k]->IsConstant();
    }
    matrix_entries << (long)(f[k]);
  }
}

void  _parser2013_addLoopContext (void *vp) {
  Parser* p = (Parser*)vp;
  p->loop_contexts.AppendNewInstance(new _SimpleList);
}


void  _parser2013_popLoopContext (void *vp, 
                                  _ExecutionList&current_command_stream, 
                                  long loop_start, long loop_end) {
  if (_parser2013_errorFree(vp) == false) return;
  Parser * p = (Parser*) vp;
  if (p->loop_contexts.countitems() > 0L) {
    _SimpleList * last_loop_context = (_SimpleList *)p->loop_contexts.Pop();
    for (long command = 0; command < last_loop_context->countitems(); command += 2) {
      _parser2013_pushSetJumpCommmandIndices (vp, current_command_stream,
                                              last_loop_context->GetElement(command),
                                              last_loop_context->GetElement(command+1)?
                                              loop_start: loop_end);
    }
  } else {
    _FormulaParsingContext fpc;
    _parser2013_reportError(vp, "_parser2013_popLoopContext called without a loop context", fpc);
  }
}



  // LL(1) resolvers

bool    _parser2013_IdentFollowedByAnOpenParenthesis (void * vp) {
  Parser* p = (Parser*)vp;
  p->scanner->ResetPeek();
  if (p->la->kind == p->_IDENTIFIER && p->scanner->Peek()->kind == p->_OPEN_PARENTHESIS) {
    return true;
  }
  return false;
}

bool    _parser2013_TwoOpenBraces (void * vp) {
  Parser* p = (Parser*)vp;
  p->scanner->ResetPeek();
  if (p->la->kind == p->_OPEN_BRACE && p->scanner->Peek()->kind == p->_OPEN_BRACE) {
    return true;
  }
  return false;
}

bool    _parser2013_isSimpleStatement (void * vp) {
  Parser* p = (Parser*)vp;
  p->scanner->ResetPeek();
  
  if (p->la->kind == p->_IDENTIFIER) {
    int next_token_kind = p->scanner->Peek()->kind;
    if (next_token_kind == p->_EQUAL || next_token_kind == p->_ASSIGN) {
      return false;
    }
  }
  return true;
}

bool    _parser2013_StringAndColon (void * vp) {
  Parser* p = (Parser*)vp;
  p->scanner->ResetPeek();
  if (p->la->kind == p->_CLOSE_BRACE) return true;
  if (p->la->kind == p->_SINGLE_QUOTE_STRING || p->la->kind == p->_DOUBLE_QUOTE_STRING) {
    return p->scanner->Peek()->kind == p->_COLON;
  }
  return false;
}

void    _parser2013_reportError    (void * vp, const _String err, _FormulaParsingContext& fpc){
  Parser* p = (Parser*)vp;
  
  _String * buffer_errors = fpc.errMsg();
  
  if (buffer_errors) {
    if (buffer_errors->sLength) {
      *buffer_errors = *buffer_errors & '\n';
    } 
    *buffer_errors = *buffer_errors & "Error in line " & (long)p->la->line
    & " column " & (long)p->la->col & ':' & 
    err;
    p->errors->count++;
  } else {
    wchar_t * buffer = new wchar_t [1+err.sLength];
    swprintf(buffer, err.sLength+1, L"%hs", err.sData);
    p->errors->Error(p->la->line, p->la->col, buffer);
    delete[] buffer;
  }
}

bool _parser2013_isFollowedByAnCommaOrClosingBrace (void *vp) {
  Parser* p = (Parser*)vp;
  p->scanner->ResetPeek();
  int la2 = p->scanner->Peek()->kind;
  if (p->la->kind == p->_MULTIPLY && (la2 == p->_COMMA || la2 == p->_CLOSE_BRACE)) {
    return true;
  }
  return false;
  
}

bool    _parser2013_errorFree                     (void * vp){
  return ((Parser*)vp)->errors->count == 0;
}


long    _parser2013_expected_arguments (const long funcId) {
  long function_list_index = _builtInArgumentCounts.FindLong (funcId);
  if (function_list_index >= 0) {
    return _builtInArgumentCounts.GetXtra(function_list_index);
  }
  return 1L;
}

void _parser2013_pushStatementOntoList (void *vp, _ExecutionList& current_command_stream, _Formula* f) {
  if (_parser2013_errorFree(vp) == false) return;
  _ElementaryCommand * a_statement = new _ElementaryCommand (_HY_HBL_COMMAND_SIMPLE_STATEMENT);  
  a_statement->AppendToSimpleParameters((long)f);
  current_command_stream.Place (a_statement);
}

void  _parser2013_handleContinueBreak (void *vp, _ExecutionList&current_command_stream, bool is_continue) {
  if (_parser2013_errorFree(vp) == false) return;
  Parser * p = (Parser*) vp;
  
  if (p->loop_contexts.countitems() == 0L) {
    _FormulaParsingContext fpc;
    _parser2013_reportError(vp, "break or continue statement outside a loop", fpc);
    return;
  } else {
    _SimpleList * most_recent_loop = (_SimpleList *)p->loop_contexts.GetElement(-1L);
    (*most_recent_loop) << current_command_stream.countitems();
    (*most_recent_loop) << is_continue;
  }
  
  _ElementaryCommand * a_statement = new _ElementaryCommand (_HY_HBL_COMMAND_JUMP_STATEMENT);
  a_statement->AppendToSimpleParameters(current_command_stream.countitems());  
  current_command_stream.Place (a_statement);
}


void _parser2013_pushJumpOntoList (void *vp, _ExecutionList& current_command_stream, _Formula* f) {
  if (_parser2013_errorFree(vp) == false) return;
  _ElementaryCommand * a_statement = new _ElementaryCommand (_HY_HBL_COMMAND_JUMP_STATEMENT); 
  a_statement->AppendToSimpleParameters(current_command_stream.countitems());
  if (f) {
    a_statement->AppendToSimpleParameters((long)f);
  }
  current_command_stream.Place (a_statement);
}

void	_parser2013_pushSetJumpCommmandIndices (void *vp, _ExecutionList& current_command_stream,
                                              long index_if, long index_then) {
  if (_parser2013_errorFree(vp) == false) return;
  _ElementaryCommand* jump_command = (_ElementaryCommand*) current_command_stream.GetElement(index_if);
  jump_command->SetSimpleParameter (0L, index_then);
                                                 
}



