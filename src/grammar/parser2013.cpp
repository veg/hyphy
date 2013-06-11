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


// parser support functions

void _parser2013_pushNumber (void *, _Formula& f, _FormulaParsingContext& fpc, const wchar_t* value) {
    wchar_t * endptr = NULL;
    double numeric_value = wcstod (value, &endptr);
    f.Push (new _Operation (new _Constant (numeric_value)));
}

void _parser2013_pushString (void *, _Formula& f, _FormulaParsingContext& fpc, const wchar_t* value) {
    _String literal (value);
    literal.Trim(1, literal.sLength - 2);
    f.Push (new _Operation (new _FString (literal)));
}

void _parser2013_pushNone (void *, _Formula& f, _FormulaParsingContext& fpc) {
    f.Push (new _Operation ( noneToken,false));
}

void _parser2013_pushOp (void *, _Formula& f, _FormulaParsingContext& fpc, long op_code, long num_terms) {
    f.Push (new _Operation ( op_code,num_terms));
}

void _parser2013_pushObject (void *, _Formula& f, _FormulaParsingContext& fpc, _PMathObj obj) {
    f.Push (new _Operation (obj));
}

void _parser2013_pushFunctionCall (void * vp, _Formula& f, _FormulaParsingContext& fpc, _String& funcId, const _List& argumentNames) {
    long arg_count   = argumentNames.countitems(),
         built_in_id = BuiltInFunctions.BinaryFind(&funcId);
         
    if (built_in_id >= 0) {
        long expect_arg = ExepectedBuiltInArguments (funcId);
        if (arg_count != expect_arg) {
            _parser2013_reportError(vp, _String(arg_count) &" arguments passed to '" & funcId & "', expected " & _String(expect_arg) & '.');
        }
        f.Push (new _Operation ( built_in_id,arg_count));
    }
    else {
        f.Push (new _Operation ( true, funcId,arg_count));
    }
}

void _parser2013_pushIdentifier (void* vp, _Formula& f, _FormulaParsingContext& fpc, const wchar_t* value, bool globalKey, bool takeVarReference) {
    _String ident (value);
    long curOpl = ident.sLength;
    if (curOpl > 2 && ident[curOpl - 1] == '_' &&
        ident[curOpl - 2] == '_') { // instant variable refrence
        _String realVarName(ident, 0, curOpl - 3);
        
        realVarName = fpc.contextualizeRef(realVarName);
        
        long realVarLoc = LocateVarByName(realVarName);
        if (realVarLoc < 0) { // bad instant variable reference
            _parser2013_reportError(vp, "Attempted to take value of undeclared variable ");
        }
        _Operation * theVar = new _Operation( realVarName, 
          true, globalKey, fpc.formulaScope(), false, true);
        f.Push(theVar);
    } else {
        if (fpc.formulaScope() &&
                _hyApplicationGlobals.Find(&ident) >= 0) {
            f.Push(new _Operation( ident, true, globalKey, nil, takeVarReference));
        } else {
            f.Push(new _Operation( ident, true, globalKey, fpc.formulaScope(), takeVarReference));
        }
    }
}

_Matrix* _parser2013_createDenseMatrix (void* vp, _SimpleList* entries, const unsigned long n_rows, 
        const unsigned long n_cols, const bool is_const) {
   
   _Matrix * m = nil;
   
    if (n_cols * n_rows != entries->countitems()) {
        _parser2013_reportError(vp, 
          "Internal error in _parser2013_createMatrixObject (incompatible entries and column_count arguments)");
        return nil;
    }
 
    if (is_const) {
      m = new _Matrix(n_rows, n_cols, false, true, false);
    } else {
      m = new _Matrix(n_rows, n_cols, false, false, true);     
      m->Convert2Formulas (); 
    }
    
    long overall_index = 0L;
    
    if (is_const) {
      for (unsigned long r = 0UL; r < n_rows; r ++) {
        for (unsigned long c = 0UL; c < n_cols; c++, overall_index++) {
          m->Store(r, c, ((_Formula*)entries->GetElement (overall_index))->Compute()->Value());
        }
      }
      entries->ClearFormulasInList();
    } else {
      for (unsigned long r = 0UL; r < n_rows; r ++) {
        for (unsigned long c = 0UL; c < n_cols; c++, overall_index++) {
          m->StoreFormula(r,c, *(_Formula*)entries->GetElement (overall_index),false);
        }
      }      
    }
      
   return m;
}

void _parser2013_matrix_checkRowLengths (void *vp, unsigned long & global_count, unsigned long& local_count) {
    if (global_count == 0L) {
      global_count = local_count;
    } else {
      if (global_count != local_count) {
          _parser2013_reportError(vp, _String ("Rows of unequal dimensions in a matrix constructor: ") & (long) global_count 
              & " vs " & (long) local_count);

      }
    }
}


void _parser2013_add_matrix_entry (_SimpleList& matrix_entries, _Formula* f, _FormulaParsingContext& fpc, bool & is_const) {
  f->SimplifyConstants();
  if (is_const) {
    is_const = f->IsConstant();
  }
  matrix_entries << (long)f;
}

// LL(1) resolvers

bool    _parser2013_isFollowedByAnOpenParenthesis (void * vp) {
    Parser* p = (Parser*)vp;
    p->scanner->ResetPeek();
    if (p->la->kind == p->_IDENTIFIER && p->scanner->Peek()->kind == p->_OPEN_PARENTHESIS) {
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

void    _parser2013_reportError    (void * vp, const _String err){
    Parser* p = (Parser*)vp;
    wchar_t * buffer = new wchar_t [1+err.sLength];
    swprintf(buffer, err.sLength+1, L"%hs", err.sData);
    p->errors->Error(p->la->line, p->la->col, buffer);
    delete[] buffer;
}

// utility functions 


