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

#ifndef __ELEMENTARYCOMMAND__
#define __ELEMENTARYCOMMAND__

#include "legacy_parser.h"
#include "trie.h"
#include <stdio.h>

struct _HBLCommandExtras {
  long cut_string;
  char extract_condition_separator;
  _SimpleList extract_conditions;
  _List command_invocation;

  bool do_trim, is_assignment, needs_verb;
};

class _ElementaryCommand :
    public _String // string contains the literal for this command
    {
public:

  _ElementaryCommand(void); //dummy default constructor
  _ElementaryCommand(long); // with operation code
  _ElementaryCommand(
      _String &command);    // process this string (and maybe an entire scope)
                            // starting at a given position
  virtual ~_ElementaryCommand(void);

  virtual BaseRef makeDynamic(void);
  virtual void Duplicate(BaseRef);
  virtual BaseRef toStr(void);

  bool Execute(_ExecutionList &); // perform this command in a given list
  void ExecuteCase0(_ExecutionList &);
  void ExecuteCase4(_ExecutionList &);
  void ExecuteCase5(_ExecutionList &);
  void ExecuteDataFilterCases(_ExecutionList &);
  void ExecuteCase11(_ExecutionList &);
  void ExecuteCase12(_ExecutionList &);
  void ExecuteCase21(_ExecutionList &);
  void ExecuteCase25(_ExecutionList &, bool = false); // fscanf
  void ExecuteCase26(_ExecutionList &);               // ReplicateConstraint
  void ExecuteCase31(_ExecutionList &);               // model construction
  void ExecuteCase32(_ExecutionList &);               // list selection handler
  void ExecuteCase34(_ExecutionList &);               // CovarianceMatrix
  void ExecuteCase36(_ExecutionList &);               // OpenDataPanel
  void ExecuteCase37(_ExecutionList &);               // GetInformation
  void ExecuteCase38(_ExecutionList &, bool);         // Reconstruct Ancestors
  void ExecuteCase39(_ExecutionList &);               // Execute Commands
  void ExecuteCase40(_ExecutionList &);               // Open Window
  void ExecuteCase41(_ExecutionList &);               // Spawn LF
  void ExecuteCase43(_ExecutionList &);               // FindRoot
  void ExecuteCase44(_ExecutionList &);               // MPISend
  void ExecuteCase45(_ExecutionList &);               // MPIReceive
  void ExecuteCase46(_ExecutionList &);               // GetDataInfo
  void ExecuteCase47(_ExecutionList &);               // ConstructStateCounter
  void ExecuteCase52(_ExecutionList &);               // Simulate
  void ExecuteCase53(_ExecutionList &);               // DoSQL
  void ExecuteCase54(_ExecutionList &);               // Topology
  void ExecuteCase55(_ExecutionList &);               // AlignSequences
  void ExecuteCase57(_ExecutionList &);               // GetNeutralNull
  void ExecuteCase58(_ExecutionList &);               // Profile Code
  void ExecuteCase61(_ExecutionList &);               // SCFG
  void ExecuteCase63(_ExecutionList &); // NN; currently not functional
  void ExecuteCase64(_ExecutionList &); // BGM

  bool HandleFprintf(_ExecutionList &);
  bool HandleHarvestFrequencies(_ExecutionList &);
  bool HandleOptimizeCovarianceMatrix(_ExecutionList &, bool);
  bool HandleComputeLFFunction(_ExecutionList &);
  bool HandleSelectTemplateModel(_ExecutionList &);
  bool HandleUseModel(_ExecutionList &);
  bool HandleSetParameter(_ExecutionList &);
  bool HandleAssert(_ExecutionList &);
  bool HandleRequireVersion(_ExecutionList &);
  bool HandleDeleteObject(_ExecutionList &);
  bool HandleClearConstraints(_ExecutionList &);
  bool HandleMolecularClock(_ExecutionList &);
  bool HandleGetURL(_ExecutionList &);
  bool HandleGetString(_ExecutionList &);
  bool HandleExport(_ExecutionList &);
  bool HandleDifferentiate(_ExecutionList &);
  long GetCode(void) { return code; }
  ;

  static _String FindNextCommand(_String &, bool = false);
  // finds & returns the next command block in input
  // chops the input to remove the newly found line

  static long ExtractConditions(_String &, long, _List &, char delimeter = ';',
                                bool includeEmptyConditions = true);
  // used to extract the loop, if-then conditions

  static bool ExtractValidateAddHBLCommand(_String &current_stream,
                                           const long command_code,
                                           _List *pieces,
                                           _HBLCommandExtras *command_spec,
                                           _ExecutionList &command_list);
  /**
   * Take a command from the current command stream, extract it, make an
   * _ElementaryCommand and add it to the execution list
   * @param current_stream -- the current command text stream
   * @param command_code   -- the numerical code (from HY_HBL_COMMAND_*)
   * @param pieces         -- the list of parameters extracted from the () part
   * of the command
   * @param command_spec   -- command specification structure
   * @param command_list   -- the command list object to append the command to
   * @return success/failure.
   */

  static bool BuildFor(_String &, _ExecutionList &, _List &);
  // builds the for loop starting from
  // the beginning of input
  // this will process the loop header
  // and the entire scope afterwards

  static bool BuildIfThenElse(_String &, _ExecutionList &, _SimpleList *);
  // builds the if-then-else construct starting from
  // the beginning of input
  // this will process the loop header
  // and the entire scope afterwards

  static bool BuildWhile(_String &, _ExecutionList &, _List &);
  // builds the while(..) construct starting from
  // the beginning of input
  // this will process the loop header
  // and the entire scope afterwards

  static bool BuildDoWhile(_String &, _ExecutionList &);
  // builds the do {} while(..); construct starting from
  // the beginning of input
  // this will process the loop header
  // and the entire scope afterwards

  static bool ProcessInclude(_String &, _ExecutionList &);
  // processes the include command

  static bool ConstructDataSet(_String &, _ExecutionList &);
  // construct a dataset from the string

  static bool ConstructExport(_String &, _ExecutionList &);
  // construct a matrix export command

  static bool ConstructGetString(_String &, _ExecutionList &);
  // construct a matrix import command

  static bool ConstructDataSetFilter(_String &, _ExecutionList &);
  // construct a dataset filter from the string

  static bool ConstructTree(_String &, _ExecutionList &);
  // construct a tree

  static bool ConstructFscanf(_String &, _ExecutionList &);
  // construct a fscanf command

  static bool ConstructExecuteCommands(_String &, _ExecutionList &);
  // construct a fscanf command

  static bool ConstructReplicateConstraint(_String &, _ExecutionList &);
  // construct a replicate constraint command

  static bool ConstructLF(_String &, _ExecutionList &);
  // construct a likelihood function

  static bool ConstructFunction(_String &, _ExecutionList &);
  // construct a fprintf command

  static bool ConstructReturn(_String &, _ExecutionList &);
  // construct a fprintf command

  static bool ConstructSetParameter(_String &, _ExecutionList &);
  // construct a set parameter clause

  static bool ConstructCategory(_String &, _ExecutionList &);
  // construct a category variable

  static bool ConstructChoiceList(_String &, _ExecutionList &);
  // construct a category variable

  static bool ConstructCategoryMatrix(_String &, _ExecutionList &);
  // construct a category matrix for the optimized like func

  static bool ConstructOpenDataPanel(_String &, _ExecutionList &);
  // open data panel with given settings

  static bool ConstructOpenWindow(_String &, _ExecutionList &);

  static bool ConstructSpawnLF(_String &, _ExecutionList &);

  static bool ConstructFindRoot(_String &, _ExecutionList &);

  static bool ConstructGetInformation(_String &, _ExecutionList &);

  static bool ConstructModel(_String &, _ExecutionList &);

  static bool ConstructMPISend(_String &, _ExecutionList &);

  static bool ConstructMPIReceive(_String &, _ExecutionList &);

  static bool ConstructGetDataInfo(_String &, _ExecutionList &);

  static bool ConstructStateCounter(_String &, _ExecutionList &);

  static bool ConstructDoSQL(_String &, _ExecutionList &);

  static bool ConstructAlignSequences(_String &, _ExecutionList &);

  static bool ConstructGetNeutralNull(_String &, _ExecutionList &);

  static bool ConstructProfileStatement(_String &, _ExecutionList &);

  static bool ConstructDeleteObject(_String &, _ExecutionList &);

  static bool ConstructSCFG(_String &, _ExecutionList &);

  static bool ConstructNN(_String &, _ExecutionList &);

  static bool ConstructBGM(_String &, _ExecutionList &);

  static bool ConstructAssert(_String &, _ExecutionList &);

  static bool SelectTemplateModel(_String &, _ExecutionList &);

  static bool MakeGeneralizedLoop(_String *, _String *, _String *, bool,
                                  _String &, _ExecutionList &);

protected:

  bool MakeJumpCommand(_String *, long, long, _ExecutionList &);
  // internal command used
  // to build a jump command
  // with two branches
  // and a condition

  void addAndClean(_ExecutionList &, _List * = nil, long = 0);

  friend class _ExecutionList;
  friend void DeleteVariable(long, bool);
  friend void UpdateChangingFlas(long);
  friend void UpdateChangingFlas(_SimpleList &);

protected: // data members

  _List parameters;             // a list of parameters
  _SimpleList simpleParameters; // a list of numeric parameters
  int code;                     // code describing this command

};

#endif
