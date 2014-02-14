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

#ifndef __EXECUTIONLIST__
#define __EXECUTIONLIST__

#include <stdio.h>
#include "legacy_parser.h"
#include "site.h"
#include "trie.h"

extern _ExecutionList *currentExecutionList;

//____________________________________________________________________________________
struct _CELInternals {
  _SimpleFormulaDatum *values, *stack;

  _SimpleList varList, storeResults;

};

#define _HY2EXECUTIONLIST(X) (dynamic_cast<_ExecutionList*>(X))

class _ExecutionList : public _List // a sequence of commands to be executed
                       {
public:
  _ExecutionList(); // doesn't do much
  _ExecutionList(_String &, _String * = nil, bool = false, bool * = nil);

  virtual ~_ExecutionList(void);

  virtual BaseRef makeDynamic(void);

  virtual BaseRef toStr(void);

  virtual void Duplicate(BaseRef);
  bool BuildList(_String &, _SimpleList * = nil, bool = false, bool = false);

  _PMathObj Execute(void); // run this execution list
  _PMathObj GetResult(void) { return result; }
  
  void ExecuteSimple(void); // run a simple compiled list
  bool
  TryToMakeSimple(void);    // see if a list can be made into a compiled version

  long ExecuteAndClean(long, _String * = nil);

  void ResetFormulae(void); // decompile formulas (for reference functions)
  void ResetNameSpace(void);
  void SetNameSpace(_String);
  
  _hyExecutionContext * GetExecutionContext (void) {return &execution_context;}
  
    
  _String GetFileName(void);
  _String *GetNameSpace(void);
  _String AddNameSpaceToID(_String &, _String * = nil);
  _String TrimNameSpaceFromID(_String &);
  _String *FetchFromStdinRedirect(void);
  _ElementaryCommand *FetchLastCommand(void) {
    if (currentCommand - 1 < lLength && currentCommand > 0) {
      return (_ElementaryCommand *)(*this)(currentCommand - 1);
    }
    return nil;
  }

  void GoToLastInstruction(void) {
    currentCommand = MAX(currentCommand, lLength - 1);
  }

  bool IsErrorState(void) { return errorState; }

  void ReportAnExecutionError(_String errMsg, bool doCommand = true,
                              bool appendToExisting = false);
  /**
   * Handle an error message according to the reporting policy of this execution
   * list (defined by errorHandlingMode)
   * @param errMsg -- the current command text stream
   * @param doCommand -- add standard text about the current command
   * @param appendToExisting -- append text to existing error
   
   */

  // data fields
  // _____________________________________________________________

  long currentCommand;
  char doProfile;
  int errorHandlingMode; // how does this execution list handle errors
  bool errorState;

  _PMathObj result;
  
  _hyExecutionContext execution_context;

  _VariableContainer *nameSpacePrefix;

  _AVLListXL *stdinRedirect;

  _List *stdinRedirectAux;

  _String sourceFile, sourceText;

  _SimpleList callPoints, lastif;

  _Matrix *profileCounter;

  _CELInternals *cli;

};

#endif
