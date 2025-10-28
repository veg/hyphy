/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon    (apoon42@uwo.ca)
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

#ifndef __BATCHLANGUAGE__
#define __BATCHLANGUAGE__

#include "associative_list.h"
#include "dataset.h"
#include "global_things.h"
#include "hy_string_buffer.h"
#include "parser.h"
#include "site.h"
#include "trie.h"

#include <stdio.h>

#define HY_BL_ERROR_HANDLING_DEFAULT 0
#define HY_BL_ERROR_HANDLING_SOFT 1

//____________________________________________________________________________________
/**
 * @brief Internals for a compiled execution list
 *
 */
struct _CELInternals {
  _SimpleFormulaDatum *values, *stack;

  bool *is_compiled;

  _SimpleList varList, storeResults;
};

//____________________________________________________________________________________
/**
 * @brief Extra information for a command
 *
 */
struct _HBLCommandExtras {
  long cut_string;
  char extract_condition_separator;
  _SimpleList extract_conditions;
  _List command_invocation;

  bool do_trim, is_assignment, needs_verb;
};

class _ElementaryCommand;

//____________________________________________________________________________________
class _ExecutionList : public _List // a sequence of commands to be executed
{

public:
  /**
   * @brief Construct a new _ExecutionList object.
   */
  _ExecutionList();
  /**
   * @brief Construct a new _ExecutionList object from a string buffer.
   *
   * @param sb The string buffer to construct from.
   * @param s An optional string.
   * @param b1 An optional boolean.
   * @param b2 An optional boolean.
   */
  _ExecutionList(_StringBuffer &, _String * = nil, bool = false, bool * = nil);
  /**
   * @brief Initialize the execution list.
   *
   * @param s An optional string.
   */
  void Init(_String * = nil);

  /**
   * @brief Destroy the _ExecutionList object
   *
   */
  virtual ~_ExecutionList(void);

  /**
   * @brief Make a dynamic copy of the execution list
   *
   * @return BaseRef
   */
  virtual BaseRef makeDynamic(void) const;
  /**
   * @brief Convert the execution list to a string
   *
   * @param ul
   * @return BaseRef
   */
  virtual BaseRef toStr(unsigned long = 0UL);
  /**
   * @brief Duplicate the execution list from a reference
   *
   * @param brc
   */
  virtual void Duplicate(BaseRefConst);
  /**
   * @brief Build the execution list from a string buffer
   *
   * @param sb
   * @param sl
   * @param b1
   * @param b2
   * @return true
   * @return false
   */
  bool BuildList(_StringBuffer &, _SimpleList * = nil, bool = false,
                 bool = false);

  /**
   * @brief Set the KWArgs object
   *
   * @param al
   */
  void SetKWArgs(_AssociativeList *);

  /**
   * @brief Execute the execution list
   *
   * @param parent
   * @param ignore_parent_kwargs
   * @return HBLObjectRef
   */
  HBLObjectRef Execute(_ExecutionList *parent = nil,
                       bool ignore_parent_kwargs = false);
  // if parent is specified, copy stdin redirects from it
  // run this execution list
  /**
   * @brief Get the Result object
   *
   * @return HBLObjectRef
   */
  HBLObjectRef GetResult(void) { return result; }
  /**
   * @brief Execute a simple compiled list
   *
   * @param parent
   */
  void
  ExecuteSimple(_ExecutionList *parent = nil); // run a simple compiled list
  /**
   * @brief Try to make a simple compiled list
   *
   * @param partial_ok
   * @return true
   * @return false
   */
  bool TryToMakeSimple(bool partial_ok = false); // see if a list can be made
                                                 // into a compiled version

  /**
   * @brief Execute and clean the list
   *
   * @param l
   */
  void ExecuteAndClean(long);

  /**
   * @brief Decompile formulas (for reference functions)
   *
   */
  void ResetFormulae(void); // decompile formulas (for reference functions)
  /**
   * @brief Reset the name space
   *
   */
  void ResetNameSpace(void);
  /**
   * @brief Set the Name Space object
   *
   * @param s
   */
  void SetNameSpace(_String const &);
  /**
   * @brief Get the File Name object
   *
   * @return _String
   */
  _String const GetFileName(void) const;
  /**
   * @brief Get the Name Space object
   *
   * @return _String*
   */
  _String *GetNameSpace(void);
  /**
   * @brief Add a name space to an ID
   *
   * @param s1
   * @param s2
   * @return _String
   */
  _String const AddNameSpaceToID(_String const &, _String const * = nil);
  /**
   * @brief Trim a name space from an ID
   *
   * @param s
   * @return _String
   */
  _String TrimNameSpaceFromID(_String &);

  /**
   * @brief Check if the list has a stdin redirect
   *
   * @return true
   * @return false
   */
  bool has_stdin_redirect(void) const { return stdinRedirect != nil; }
  /**
   * @brief Check if the list has keyword arguments
   *
   * @return true
   * @return false
   */
  bool has_keyword_arguments(void) const {
    return (kwargs && kwargs->countitems()) ||
           (kwarg_tags && kwarg_tags->countitems());
  }

  /**
   * @brief Fetch from the stdin redirect
   *
   * @param dialog_tag
   * @param handle_multi_choice
   * @param do_echo
   * @return _String*
   */
  _String *FetchFromStdinRedirect(_String const *dialog_tag = nil,
                                  bool handle_multi_choice = false,
                                  bool do_echo = false);

  /**
   * @brief Fetch the last command
   *
   * @return _ElementaryCommand*
   */
  _ElementaryCommand *FetchLastCommand(void) {
    if (currentCommand - 1L < (long)lLength && currentCommand > 0L) {
      return (_ElementaryCommand *)(*this)(currentCommand - 1L);
    }
    return nil;
  }
  /**
   * @brief Get the Ith Command object
   *
   * @param i
   * @return _ElementaryCommand*
   */
  _ElementaryCommand *GetIthCommand(long i) const {
    return (_ElementaryCommand *)GetItem(i);
  }

  /**
   * @brief Go to the last instruction
   *
   */
  void GoToLastInstruction(void) {
    currentCommand = MAX(currentCommand, (long)lLength - 1L);
  }

  /**
   * @brief Generate a help message
   *
   * @param options
   * @param inputs
   * @param scanned_functions
   * @return _StringBuffer
   */
  _StringBuffer const GenerateHelpMessage(_List *options = nil,
                                          _List *inputs = nil,
                                          _AVLListXL *scanned_functions = nil);

  /**
   * @brief Check if the list is in an error state
   *
   * @return true
   * @return false
   */
  bool IsErrorState(void) { return errorState; }

  /**
   * @brief Check if this is a dry run
   *
   * @return true
   * @return false
   */
  bool IsDryRun(void) { return dry_run; }

  /**
   * @brief Set the Dry Run object
   *
   * @param dr
   */
  void SetDryRun(bool dr) { dry_run = dr; }

  /**
   * @brief Clear the execution list
   *
   */
  void ClearExecutionList(void);

  /**
   * @brief Report an execution error
   *
   * @param errMsg
   * @param doCommand
   * @param appendToExisting
   * @param isAssertion
   */
  void ReportAnExecutionError(_String errMsg, bool doCommand = true,
                              bool appendToExisting = false,
                              bool isAssertion = false);
  /**
   * Handle an error message according to the reporting policy of this execution
   list (defined by errorHandlingMode)
   * @param errMsg -- the current command text stream
   * @param doCommand -- add standard text about the current command
   * @param appendToExisting -- append text to existing error
   * @param isAssertion -- is this call triggered by an assertion (this will
   result in different error reporting)

   */

  /**
   * @brief Scan the body of this function/code for dependancies on various
   * objects (currently only supports HBL functions), and store them in
   * `collection`.
   *
   * @param collection
   * @param recursive if true, then new HBL functions will be scanned for
   * dependancies as well
   * @param help_mode
   */
  void BuildListOfDependancies(_AVLListX &collection, bool recursive = true,
                               bool help_mode = false);

  /**
   * @brief Advance program counter
   *
   */
  void advance(void) { currentCommand++; }

  /**
   * @brief Check if the list is compiled
   *
   * @param idx
   * @return true
   * @return false
   */
  bool is_compiled(long idx = -1) const {
    if (cli) {
      if (idx < 0L)
        return true;
      else
        return cli->is_compiled[idx];
    }
    return false;
  }

  /**
   * @brief Copy CLI to variables
   *
   */
  void CopyCLIToVariables(void);
  /**
   * @brief Start profiling
   *
   */
  void StartProfile(void);
  /**
   * @brief Collect profiling data
   *
   * @return _AssociativeList*
   */
  _AssociativeList *CollectProfile(void);

  // data fields
  // _____________________________________________________________

  long currentCommand, currentKwarg;

  char doProfile;
  int errorHandlingMode; // how does this execution list handle errors
  bool errorState;
  bool dry_run;

  HBLObjectRef result;

  _VariableContainer *nameSpacePrefix;
  _AssociativeList *kwargs;

  _AVLListXL *stdinRedirect;
  _List *stdinRedirectAux, *kwarg_tags;

  /** SLKP 20190223
      kwarg_tags, if set, stores the ordered list of tagged inputs, which are
     created by invoking the KeywordArgument ("keyword", "description",
     "default"); procedure. They are inherited down the call chain, like
     stdinRedirect and stdinRedirectAux
   */

  _String sourceFile, sourceText, enclosingNamespace;

  _SimpleList callPoints, lastif;

  _Matrix *profileCounter;

  _CELInternals *cli;

protected:
  void BuildExecuteCommandInstruction(_List *pieces, long code);
  void BuildFscanf(_List *pieces, long code);
  void BuildChoiceList(_List *pieces, long code);
};

//____________________________________________________________________________________
// an elementary command

/**
 * @brief An elementary command
 *
 */
class _ElementaryCommand
    : public _String // string contains the literal for this command
{
public:
  /**
   * @brief Construct a new _ElementaryCommand object.
   */
  _ElementaryCommand(void);
  /**
   * @brief Construct a new _ElementaryCommand object with an operation code.
   *
   * @param l The operation code.
   */
  _ElementaryCommand(long);
  /**
   * @brief Construct a new _ElementaryCommand object from a string.
   * This will process the string and create the command.
   *
   * @param command The string to construct from.
   */
  _ElementaryCommand(_String &command);
  // starting at a given position
  /**
   * @brief Destroy the _ElementaryCommand object
   *
   */
  virtual ~_ElementaryCommand(void);

  /**
   * @brief Make a dynamic copy of the command
   *
   * @return BaseRef
   */
  virtual BaseRef makeDynamic(void) const;
  /**
   * @brief Duplicate the command from a reference
   *
   * @param brc
   */
  virtual void Duplicate(BaseRefConst);
  /**
   * @brief Convert the command to a string
   *
   * @param ul
   * @return BaseRef
   */
  virtual BaseRef toStr(unsigned long = 0UL);

  /**
   * @brief Execute the command
   *
   * @param el
   * @return true
   * @return false
   */
  bool Execute(_ExecutionList &); // perform this command in a given list
  /**
   * @brief @sergeilkp
   *
   * @param el
   */
  void ExecuteCase5(_ExecutionList &);
  /**
   * @brief @sergeilkp
   *
   * @param el
   */
  void ExecuteDataFilterCases(_ExecutionList &);
  /**
   * @brief @sergeilkp
   *
   * @param el
   */
  void ExecuteCase11(_ExecutionList &);
  /**
   * @brief @sergeilkp
   *
   * @param el
   */
  void ExecuteCase12(_ExecutionList &);
  /**
   * @brief @sergeilkp
   *
   * @param el
   */
  void ExecuteCase31(_ExecutionList &); // model construction
  /**
   * @brief @sergeilkp
   *
   * @param el
   * @param b
   */
  void ExecuteCase38(_ExecutionList &, bool); // Reconstruct Ancestors
  /**
   * @brief @sergeilkp
   *
   * @param el
   */
  void ExecuteCase52(_ExecutionList &); // Simulate
  /**
   * @brief @sergeilkp
   *
   * @param el
   */
  void ExecuteCase54(_ExecutionList &); // Topology
  /**
   * @brief @sergeilkp
   *
   * @param el
   */
  void ExecuteCase58(_ExecutionList &); // Profile Code
  /**
   * @brief @sergeilkp
   *
   * @param el
   */
  void ExecuteCase61(_ExecutionList &); // SCFG
  /**
   * @brief @sergeilkp
   *
   * @param el
   */
  void ExecuteCase63(_ExecutionList &); // NN; currently not functional
  /**
   * @brief @sergeilkp
   *
   * @param el
   */
  void ExecuteCase64(_ExecutionList &); // BGM

  /**
   * @brief Handle a conditional branch
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleConditionalBranch(_ExecutionList &);
  /**
   * @brief Handle a generic expression
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleGenericExpression(_ExecutionList &);
  /**
   * @brief Handle a replicate constraint
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleReplicateConstraint(_ExecutionList &);
  /**
   * @brief Handle aligning sequences
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleAlignSequences(_ExecutionList &);
  /**
   * @brief Handle constructing a category matrix
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleConstructCategoryMatrix(_ExecutionList &);
  /**
   * @brief Handle getting data info
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleGetDataInfo(_ExecutionList &);
  /**
   * @brief Handle getting information
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleGetInformation(_ExecutionList &);
  /**
   * @brief Handle fprintf
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleFprintf(_ExecutionList &);
  /**
   * @brief Handle harvesting frequencies
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleHarvestFrequencies(_ExecutionList &);
  /**
   * @brief Handle optimizing a covariance matrix
   *
   * @param el
   * @param b
   * @return true
   * @return false
   */
  bool HandleOptimizeCovarianceMatrix(_ExecutionList &, bool);
  /**
   * @brief Handle computing a LF function
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleComputeLFFunction(_ExecutionList &);
  /**
   * @brief Handle selecting a template model
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleSelectTemplateModel(_ExecutionList &);
  /**
   * @brief Handle using a model
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleUseModel(_ExecutionList &);
  /**
   * @brief Handle setting a parameter
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleSetParameter(_ExecutionList &);
  /**
   * @brief Handle an assert
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleAssert(_ExecutionList &);
  /**
   * @brief Handle requiring a version
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleRequireVersion(_ExecutionList &);
  /**
   * @brief Handle deleting an object
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleDeleteObject(_ExecutionList &);
  /**
   * @brief Handle clearing constraints
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleClearConstraints(_ExecutionList &);
  /**
   * @brief Handle a molecular clock
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleMolecularClock(_ExecutionList &);
  /**
   * @brief Handle getting a URL
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleGetURL(_ExecutionList &);
  /**
   * @brief Handle getting a string
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleGetString(_ExecutionList &);
  /**
   * @brief Handle a keyword argument
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleKeywordArgument(_ExecutionList &);
  /**
   * @brief Handle an export
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleExport(_ExecutionList &);
  /**
   * @brief Handle a differentiate
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleDifferentiate(_ExecutionList &);
  /**
   * @brief Handle finding a root or integrating
   *
   * @param el
   * @param do_integrate
   * @return true
   * @return false
   */
  bool HandleFindRootOrIntegrate(_ExecutionList &, bool do_integrate = false);
  /**
   * @brief Handle converting a branch length
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleConvertBranchLength(_ExecutionList &);
  /**
   * @brief Handle an MPI send
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleMPISend(_ExecutionList &);
  /**
   * @brief Handle an MPI receive
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleMPIReceive(_ExecutionList &);
  /**
   * @brief Handle executing commands cases
   *
   * @param el
   * @param do_load_from_file
   * @param do_load_library
   * @return true
   * @return false
   */
  bool HandleExecuteCommandsCases(_ExecutionList &,
                                  bool do_load_from_file = false,
                                  bool do_load_library = false);
  /**
   * @brief Handle fscanf
   *
   * @param el
   * @param is_sscanf
   * @return true
   * @return false
   */
  bool HandleFscanf(_ExecutionList &, bool is_sscanf = false);
  /**
   * @brief Handle a choice list
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleChoiceList(_ExecutionList &);
  /**
   * @brief Handle initializing an iterator
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleInitializeIterator(_ExecutionList &);
  /**
   * @brief Handle advancing an iterator
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleAdvanceIterator(_ExecutionList &);
  /**
   * @brief Handle a nested list
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleNestedList(_ExecutionList &);
  /**
   * @brief Handle tree construction
   *
   * @param el
   * @return true
   * @return false
   */
  bool HandleTreeConstruction(_ExecutionList &);

  /**
   * @brief Get the code of the command
   *
   * @return long
   */
  long get_code(void) const { return code; };
  /**
   * @brief Get the number of parameters
   *
   * @return unsigned
   */
  unsigned long parameter_count(void) const { return parameters.countitems(); }

  /**
   * @brief Find the next command
   *
   * @param sb1
   * @param sb2
   */
  static void FindNextCommand(_StringBuffer &, _StringBuffer &);
  // finds & stores the next command from _String into _StringBuffer
  // chops the input to remove the newly found line

  /**
   * @brief Extract conditions from a string
   *
   * @param sb
   * @param l
   * @param list
   * @param delimeter
   * @param includeEmptyConditions
   * @return long
   */
  static long ExtractConditions(_StringBuffer const &, long, _List &,
                                char delimeter = ';',
                                bool includeEmptyConditions = true);
  // used to extract the loop, if-then conditions

  /**
   * @brief Take a command from the current command stream, extract it, make an
   * _ElementaryCommand and add it to the execution list
   *
   * @param current_stream the current command text stream
   * @param command_code the numerical code (from HY_HBL_COMMAND_*)
   * @param pieces the list of parameters extracted from the () part of the
   * command
   * @param command_spec command specification structure
   * @param command_list the command list object to append the command to
   * @return true
   * @return false
   */
  static bool ExtractValidateAddHBLCommand(_StringBuffer &current_stream,
                                           const long command_code,
                                           _List *pieces,
                                           _HBLCommandExtras *command_spec,
                                           _ExecutionList &command_list);

  /**
   * @brief Build a for loop
   *
   * @param sb
   * @param el
   * @param l
   * @return true
   * @return false
   */
  static bool BuildFor(_StringBuffer &, _ExecutionList &, _List *);
  // builds the for loop starting from
  // the beginning of input
  // this will process the loop header
  // and the entire scope afterwards

  /**
   * @brief Build an if-then-else construct
   *
   * @param sb
   * @param el
   * @param sl
   * @return true
   * @return false
   */
  static bool BuildIfThenElse(_StringBuffer &, _ExecutionList &, _SimpleList *);
  // builds the if-then-else construct starting from
  // the beginning of input
  // this will process the loop header
  // and the entire scope afterwards

  /**
   * @brief Build a while loop
   *
   * @param sb
   * @param el
   * @param l
   * @return true
   * @return false
   */
  static bool BuildWhile(_StringBuffer &, _ExecutionList &, _List *);
  // builds the while(..) construct starting from
  // the beginning of input
  // this will process the loop header
  // and the entire scope afterwards

  /**
   * @brief Build a do-while loop
   *
   * @param sb
   * @param el
   * @return true
   * @return false
   */
  static bool BuildDoWhile(_StringBuffer &, _ExecutionList &);
  // builds the do {} while(..); construct starting from
  // the beginning of input
  // this will process the loop header
  // and the entire scope afterwards

  /**
   * @brief Process an include command
   *
   * @param sb
   * @param el
   * @return true
   * @return false
   */
  static bool ProcessInclude(_StringBuffer &, _ExecutionList &);
  // processes the include command

  /**
   * @brief Construct a dataset from a string
   *
   * @param sb
   * @param el
   * @return true
   * @return false
   */
  static bool ConstructDataSet(_StringBuffer &, _ExecutionList &);
  // construct a dataset from the string

  /**
   * @brief Construct a dataset filter from a string
   *
   * @param sb
   * @param el
   * @return true
   * @return false
   */
  static bool ConstructDataSetFilter(_StringBuffer &, _ExecutionList &);
  // construct a dataset filter from the string

  /**
   * @brief Construct a tree from a string
   *
   * @param sb
   * @param el
   * @return true
   * @return false
   */
  static bool ConstructTree(_StringBuffer &, _ExecutionList &);
  // construct a tree

  /**
   * @brief Construct a likelihood function from a string
   *
   * @param sb
   * @param el
   * @return true
   * @return false
   */
  static bool ConstructLF(_StringBuffer &, _ExecutionList &);
  // construct a likelihood function

  /**
   * @brief Construct a function from a string
   *
   * @param sb
   * @param el
   * @return true
   * @return false
   */
  static bool ConstructFunction(_StringBuffer &, _ExecutionList &);
  // construct a fprintf command

  /**
   * @brief Construct a return command
   *
   * @param sb
   * @param el
   * @return true
   * @return false
   */
  static bool ConstructReturn(_StringBuffer &, _ExecutionList &);
  // construct a fprintf command

  /**
   * @brief Construct a category variable from a string
   *
   * @param sb
   * @param el
   * @return true
   * @return false
   */
  static bool ConstructCategory(_StringBuffer &, _ExecutionList &);
  // construct a category variable

  /**
   * @brief Construct a choice list from a string
   *
   * @param sb
   * @param el
   * @return true
   * @return false
   */
  static bool ConstructChoiceList(_StringBuffer &, _ExecutionList &);
  // construct a category variable

  /**
   * @brief Construct a model from a string
   *
   * @param sb
   * @param el
   * @return true
   * @return false
   */
  static bool ConstructModel(_StringBuffer &, _ExecutionList &);

  /**
   * @brief Construct a profile statement from a string
   *
   * @param sb
   * @param el
   * @return true
   * @return false
   */
  static bool ConstructProfileStatement(_StringBuffer &, _ExecutionList &);

  /**
   * @brief Construct a SCFG from a string
   *
   * @param sb
   * @param el
   * @return true
   * @return false
   */
  static bool ConstructSCFG(_StringBuffer &, _ExecutionList &);

  /**
   * @brief Construct a BGM from a string
   *
   * @param sb
   * @param el
   * @return true
   * @return false
   */
  static bool ConstructBGM(_StringBuffer &, _ExecutionList &);

  /**
   * @brief Make a generalized loop
   *
   * @param s1
   * @param s2
   * @param s3
   * @param b
   * @param sb
   * @param el
   * @return true
   * @return false
   */
  static bool MakeGeneralizedLoop(_String *, _String *, _String *, bool,
                                  _StringBuffer &, _ExecutionList &);

  /**
   * @brief Decompile formulae
   *
   * @return true
   * @return false
   */
  bool DecompileFormulae(void);

  /**
   * @brief Check this command for (currently only supports HBL functions), and
   * store them in `collection`.
   *
   * @param collection
   * @param recursive if true, then new HBL functions will be scanned for
   * dependancies as well
   * @param chain
   * @param help_mode
   */
  void BuildListOfDependancies(_AVLListX &collection, bool recursive,
                               _ExecutionList const &chain,
                               bool help_mode = false);

  static const _List fscanf_allowed_formats;

protected:
  /**
   * @brief Scan a string expression for HBL functions
   *
   * @param s
   * @param elc
   * @param b
   * @param avl
   * @param help_mode
   */
  static void ScanStringExpressionForHBLFunctions(_String *,
                                                  _ExecutionList const &, bool,
                                                  _AVLListX &,
                                                  bool help_mode = false);

  /**
   * @brief Get the Ith Parameter object
   *
   * @param i
   * @param range_check
   * @return _String*
   */
  _String *GetIthParameter(unsigned long i, bool range_check = true) const {
    BaseRef p = parameters.GetItemRangeCheck(i);
    if (!p && range_check) {
      hy_global::HandleApplicationErrorAndExit(
          "Internal error in ElemenaryCommand::GetIthParameter");
      return nil;
    }
    return (_String *)p;
  }

  /**
   * @brief Make a jump command
   *
   * @param s
   * @param l1
   * @param l2
   * @param el
   * @return true
   * @return false
   */
  bool MakeJumpCommand(_String *, long, long, _ExecutionList &);
  // internal command used
  // to build a jump command
  // with two branches
  // and a condition

  /**
   * @brief Add and clean
   *
   * @param el
   * @param l
   * @param c
   */
  void addAndClean(_ExecutionList &, _List * = nil, long = 0);
  /**
   * @brief Append compiled formulae
   *
   * @param f1
   * @param f2
   */
  void appendCompiledFormulae(_Formula *, _Formula * = nil);

  friend class _ExecutionList;
  friend void DeleteVariable(long, bool);
  friend void UpdateChangingFlas(long);
  friend void UpdateChangingFlas(_SimpleList const &);

private:
  /**
   * @brief Validate a storage variable
   *
   * @param program
   * @param argument_index
   * @return _Variable*
   */
  _Variable *_ValidateStorageVariable(_ExecutionList &program,
                                      unsigned long argument_index = 0UL) const;

  /**
      Extract the identifier from an expression like
      Type <id> = (...);

      @param ([in] _String) the source text to scan
      @param ([out] long) a positive integer to the first character following
     the '='; kNotFound if missing
      @param ([in] bool) if true, validate the identifier
      @param ([in] bool) if true, throw _String exceptions on errors
      @param ([in] long) start searching at this position in the string

      @return the identifier, or empty string if failed
      @version 0.1 SLKP 20190211
  */
  static const _String ExtractStatementAssignment(_String const &source,
                                                  long &end_at,
                                                  const bool validate = true,
                                                  const bool exceptions = true,
                                                  const long offset = 0L);

  /**
   Process declaration of the form
   Type <id> = procedure (...);

   @param ([in] _String) the source text to scan
   @param ([out] procedure) the string for the procedure name (e.g.
   CreateFilter)
   @param ([out] pieces) comma separated arguments from the parentheses

   @return the identifier, or empty string if failed

   Throws _String exceptions

   @version 0.1 SLKP 20190212
   */
  static const _String ProcessProcedureCall(_String const &source,
                                            _String &procedure, _List &pieces);

protected:                      // data members
  _List parameters;             // a list of parameters
  _SimpleList simpleParameters; // a list of numeric parameters
  long code;                    // code describing this command
};

//____________________________________________________________________________________

_ElementaryCommand *makeNewCommand(long);

//____________________________________________________________________________________

#ifdef __HYPHYMPI__
#include <mpi.h>

#define HYPHY_MPI_SIZE_TAG 111
#define HYPHY_MPI_STRING_TAG 112
#define HYPHY_MPI_DONE_TAG 113
#define HYPHY_MPI_VARS_TAG 114
#define HYPHY_MPI_DATA_TAG 115

#define HYPHY_MPI_DIE_TAG 666

void ReportMPIError(int, bool);
void MPISendString(_String const &, long, bool = false);
_String *MPIRecvString(long, long &);

#endif
//____________________________________________________________________________________

//  EXTERN GLOBALS TODO SLKP 20180920 : these need to be reviewed and removed

extern _List

    dataSetList,
    dataSetNamesList, likeFuncList, templateModelList, scfgNamesList, scfgList,
    batchLanguageFunctions,

    bgmNamesList, // modified by afyp
    bgmList,

    likeFuncNamesList, modelNames, executionStack, compiledFormulaeParameters,
    standardLibraryPaths, standardLibraryExtensions;

extern _SimpleList modelMatrixIndices, modelTypeList,
    // SLKP: 20100313 this list stores 0 for  normal (rate-matrix based models),
    //       vs expression based matrices, for which the dimension is stored.
    modelFrequenciesIndices, listOfCompiledFormulae;

extern _String

    useLastFString,
    mpiMLELFValue, lf2SendBack, defFileNameValue, defFileString, blConstructCM,
    blFprintfRedirect, blFprintfDevNull, globalPolynomialCap,
    enforceGlobalPolynomialCap, dropPolynomialTerms, maxPolyTermsPerVariable,
    maxPolyExpIterates, polyExpPrecision, systemVariableDump, selfDump,
    printDigitsSpec, explicitFormMExp, multByFrequencies, getDString,
    useLastFString, getFString, defFileString, VerbosityLevelString, clearFile,
    keepFileOpen, closeFile, useLastDefinedMatrix, MessageLogging,
    selectionStrings, stdoutDestination, messageLogDestination,
    dataPanelSourcePath, windowTypeTree, windowTypeClose, windowTypeTable,
    windowTypeDistribTable, windowTypeDatabase, screenWidthVar, screenHeightVar,
    useNexusFileData, mpiMLELFValue, lf2SendBack, lfStartCompute, lfDoneCompute,
    getURLFileFlag, versionString, timeStamp, listLoadedLibraries,
    simulationFilter, prefixDS, prefixDF, prefixLF, replaceTreeStructure,
    hyphyBaseDirectory, platformDirectorySeparator, covarianceParameterList,
    matrixEvalCount, scfgCorpus, pathToCurrentBF, errorReportFormatExpression,
    errorReportFormatExpressionStr, errorReportFormatExpressionStack,
    errorReportFormatExpressionStdin, deferConstrainAssignment, kBGMData,
    bgmConstraintMx, bgmParameters, bgmWeights, assertionBehavior, dialogPrompt,
    _hyLastExecutionError, _hyExecutionErrorMode, blReturn, blDataSet,
    blDataSetFilter, blLF, blLF3, blTree, blTopology, blSCFG;

extern _ExecutionList *currentExecutionList;

extern _AVLList loadedLibraryPaths;
extern _AVLListX _HY_HBLCommandHelper, _HY_GetStringGlobalTypes;

extern _Trie _HY_ValidHBLExpressions, _HY_HBL_Namespaces,
    _HY_HBL_KeywordsPreserveSpaces;

const _String &GetBFFunctionNameByIndex(long);
long GetBFFunctionArgumentCount(long);
_List &GetBFFunctionArgumentList(long);
_SimpleList &GetBFFunctionArgumentTypes(long);
hyBLFunctionType GetBFFunctionType(long);
_ExecutionList &GetBFFunctionBody(long);

_String const ExportBFFunction(long, bool = true, _AVLList * = nil);

void ClearBFFunctionLists(long = -1L);
bool IsBFFunctionIndexValid(long);
long GetBFFunctionCount(void);

void ScanModelForVariables(long modelID, _AVLList &theReceptacle, bool inclG,
                           long modelID2, bool inclCat);
/* 20100316 SLKP:
    factored out a function call to scan a particular model
    for variables to permit the use of explicit (formula-based) model
   definitions
 */

void ReadBatchFile(_String &, _ExecutionList &);
_String const ReturnDialogInput(bool dispPath = false,
                                _String const *rel_path = nil);
_String const ReturnFileDialogInput(_String const *rel_path = nil);
_String const WriteFileDialogInput(_String const *rel_path = nil);

hyFloat _ProcessNumericArgumentWithExceptions(_String &,
                                              _VariableContainer const *);

hyFloat ProcessNumericArgument(_String *, _VariableContainer const *,
                               _ExecutionList * = nil);
const _String ProcessLiteralArgument(_String const *,
                                     _VariableContainer const *,
                                     _ExecutionList * = nil);
_AssociativeList *ProcessDictionaryArgument(_String *data,
                                            _VariableContainer *theP,
                                            _ExecutionList * = nil);

const _String GetStringFromFormula(_String const *, _VariableContainer *);

void SerializeModel(_StringBuffer &, long, _AVLList * = nil, bool = false,
                    _AssociativeList *options = nil);
bool Get_a_URL(_String &, _String * = nil);

long AddDataSetToList(_String &, _DataSet *);
void KillLFRecord(long, bool = true);
void KillDataSetRecord(long);
void KillModelRecord(long);
void KillExplicitModelFormulae(void);
bool PushFilePath(_String &, bool = true, bool process = true);
_String const PopFilePath(void);
_String const *PeekFilePath(void);
_String const GetPathStack(const _String &spacer = ",");

void RetrieveModelComponents(long, _Matrix *&, _Matrix *&, bool &);
void RetrieveModelComponents(long, _Variable *&, _Variable *&, bool &);
long RetrieveModelFreq(long);

void ReadModelList(void);
_String ProcessStringArgument(_String *data);

const _String _hblCommandAccessor(_ExecutionList *, long);
_String const _HYGenerateANameSpace(void);
void _HYClearANameSpace(_String const &);

/**
  Check if a valid (defined) HBL function with the specified number of
  arguments. If not, exit with an error, or throw an exception otherwise (with
  the error message) return a formula for the callback. Not found will not
  return (exception or exit)
 */

_Formula *ValidateCallbackFunctionArgument(_String const &function_id,
                                           unsigned long argument_count,
                                           bool trap_errors = true,
                                           _String const *error_msg = nullptr);
HBLObjectRef ExecuteCallbackFunction(_Formula *callback, _List const &arguments,
                                     unsigned long valid_type = HY_ANY_OBJECT);

HBLObjectRef ProcessAnArgumentByType(_String const *,
                                     _VariableContainer const *, long,
                                     _ExecutionList * = nil);

void _HBL_Init_Const_Arrays(void);

/**
 An accessor function which attempts to retrieve a reference to a HyPhy Batch
 Language Object by name. A list of acceptable object classes can be specified
 in the type parameter. Note that types will be searched in the following order:

 HY_BL_DATASET,HY_BL_DATASET_FILTER,HY_BL_LIKELIHOOD_FUNCTION,HY_BL_SCFG,HY_BL_BGM,HY_BL_MODEL,HY_BL_HBL_FUNCTION

 i.e. if there is a dataset named 'foo' and a likelihood function named 'foo',
 then the dataset will be returned.



 @param   name provides a string with the name of the object to be retrieved.
 @param   type [in] which types of objects will be searched.
 [out] which type of object was retrieved (HY_BL_NOT_DEFINED if not found)
 @param   index (if not nil) will receive the index of the found object in the
 corresponding array
 @param   errMsg if set to True, will cause the function to report an error if
 no object of corresponding type could be found
 @param   tryLiteralLookup if set to True, will cause the function to, upon a
 failed lookup, to also try interpreting name as a string variable ID
 @return  pointer to the retrieved object or nil if not found
 @author  SLKP
 @version 20120324
 */

_String const _HYHBLTypeToText(long type);

_HBLCommandExtras *
_hyInitCommandExtras(const long = 0, const long = 0,
                     const _String & = hy_global::kEmptyString,
                     const char = ';', const bool = true, const bool = false,
                     const bool = false, _SimpleList * = nil);

extern bool numericalParameterSuccessFlag;
extern hyFloat messageLogFlag;

extern enum _hy_nested_check {
  _HY_NO_FUNCTION,
  _HY_FUNCTION,
  _HY_NAMESPACE
} isInFunction;

#endif
