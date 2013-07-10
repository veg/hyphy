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

#include "alignment.h"
#include "batchlan.h"
#include "batchlan_globals.h"
#include "baseobj.h"
#include "bayesgraph.h"
#include "elementarycommand.h"
#include "executionlist.h"
#include "hy_globals.h"
#include "HYNetInterface.h"
#include "thetree.h"
#include "datasetfilter.h"
#include "likefunc.h"
#include "ctype.h"
#include "scfg.h"
#include "customfunction.h"
#include "thetree.h"

#ifndef __HYPHY_NO_SQLITE__
#include "sqlite3.h"
#endif

#if !defined __HEADLESS__ && !defined __UNIX__
#include      "HYUtils.h"
#endif

#ifdef __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

int _HYSQLCallBack(void *data, int callCount);

int _HYSQLCallBack(void *exL, int cc, char **rd, char **cn) {
  _ExecutionList *exList = (_ExecutionList *)exL;

  if (!terminateExecution)
    if (exList && cc && exList->lLength) {
      _List rowData, columnNames;

      for (long cnt = 0; cnt < cc; cnt++) {
        if (rd[cnt]) {
          rowData.AppendNewInstance(new _String(rd[cnt]));
        } else {
          rowData.AppendNewInstance(new _String);
        }

        if (cn[cnt]) {
          columnNames.AppendNewInstance(new _String(cn[cnt]));
        } else {
          columnNames.AppendNewInstance(new _String);
        }
      }

      _Matrix *rowDataM = new _Matrix(rowData),
              *columnNamesM = new _Matrix(columnNames);

      if (!(rowDataM && columnNamesM)) {
        checkPointer(nil);
      }

      _Variable *rdv = CheckReceptacle(&sqlRowData, blDoSQL, false),
                *cnv = CheckReceptacle(&sqlColNames, blDoSQL, false);

      rdv->SetValue(rowDataM, false);
      cnv->SetValue(columnNamesM, false);

      exList->Execute();

    }

  return 0;

}

//______________________________________________________________________________
_ElementaryCommand::_ElementaryCommand(void) { code = -1; }

//______________________________________________________________________________
_ElementaryCommand::_ElementaryCommand(long ccode) { code = ccode; }

//______________________________________________________________________________
_ElementaryCommand::_ElementaryCommand(_String &s) {
  code = -1;
  _String::Duplicate(&s);
}



//______________________________________________________________________________
_ElementaryCommand::~_ElementaryCommand(void) {
  if (nInstances == 1) {
      switch (code) {
        case _HY_HBL_COMMAND_SIMPLE_STATEMENT:
          delete (_Formula*) simpleParameters.GetElement(0L);
          break;
        case _HY_HBL_COMMAND_JUMP_STATEMENT:
          if (simpleParameters.lLength > 1) {
            delete (_Formula*) simpleParameters.GetElement(1L);
          }

      }
  }

}

//______________________________________________________________________________
void    _ElementaryCommand::AppendToSimpleParameters (const long sp) {
  simpleParameters << sp;
}

//______________________________________________________________________________
void    _ElementaryCommand::SetSimpleParameter (const long index, const long sp) {
  simpleParameters[index] = sp;
}


//______________________________________________________________________________
void    _ElementaryCommand::AppendToParameters (const BaseRef p) {
  parameters.Place (p);
}



//______________________________________________________________________________
BaseRef _ElementaryCommand::makeDynamic(void) {
  _ElementaryCommand *nec = new _ElementaryCommand;
  if (!nec) {
    return nil;
  }
  nec->code = code;

  //memcpy ((char*)Res, (char*)this, sizeof (_ElementaryCommand));
  //if (code == 4 || code==6 || code==9 || code==0)
  //nInstances++;

  nec->nInstances = 1;
  nec->Duplicate(this);
  return nec;
}

//______________________________________________________________________________
void _ElementaryCommand::Duplicate(BaseRef source) {

  _ElementaryCommand *sec = (_ElementaryCommand *)source;

  //if (code == 0)
  //_String::DuplicateErasing (source);
  //else
  _String::Duplicate(source);

  parameters.Duplicate(&sec->parameters);
  if (code != 0) {
    simpleParameters.Duplicate(&sec->simpleParameters);
  }
}

//______________________________________________________________________________
BaseRef _ElementaryCommand::toStr(void) {
  _String result, *converted = nil;
  long k;
  switch (code) {

    case _HY_HBL_COMMAND_SIMPLE_STATEMENT: {
      _Formula *f = (_Formula*) simpleParameters.GetElement(0L);
      result = _String("A simple HBL statement")
        & _String ((_String*)f->GetList().toStr());
    }
    break;

  case _HY_HBL_COMMAND_JUMP_STATEMENT:

    result = _String ("Jump to ") & simpleParameters.GetElement (0);
    if (simpleParameters.lLength > 1) {
      result = result & " subject to a condition";
    } else {
      result = result & " unconditionally";
    }
    break;

  case 5: // data set contruction
    converted = (_String *)parameters(0)->toStr();
    result = _String("Read Data Set ") & (*converted);
    DeleteObject(converted);
    converted = (_String *)parameters(1)->toStr();
    result = result & _String(" from file ") & (*converted);

    break;

  case 6:                                  // data set filter construction
  case HY_HBL_COMMAND_HARVEST_FREQUENCIES: // or HarvestFrequencies
    if (code == HY_HBL_COMMAND_HARVEST_FREQUENCIES) {
      result = "Harvest Frequencies into matrix ";
      converted = (_String *)(parameters(0)->toStr());
      result = result & *converted;
      DeleteObject(converted);
      converted = (_String *)(parameters(1)->toStr());
      result = result & " from DataSet " & *converted;
    } else {
      result = "Build Data Set Filter ";
      converted = (_String *)(parameters(0)->toStr());
      result = result & *converted;
      DeleteObject(converted);
      converted = (_String *)(parameters(1)->toStr());
      result = result & (_String)(" from DataSet ") & *converted;
    }

    DeleteObject(converted);

    if (parameters.lLength > 2) {
      result = result & " with unit size = " & *(_String *)parameters(2);
    }

    if (code == 9) {

      result = result & " with atom size = " & *(_String *)parameters(3);
      result = result & " with position specific flag = " &
               *(_String *)parameters(4);
    }

    result = result & " Partition:";
    for (k = code == 6 ? 3 : 5; k < parameters.countitems(); k++) {
      result = result & "," & *(_String *)parameters(k);
    }

    converted = nil;

    break;

  case 7:  // build a tree
  case 54: // build a tree

    converted = (_String *)parameters(0)->toStr();
    if (code == 7) {
      result = _String("Build Tree ") & *converted;
    } else {
      result = _String("Build Topology ") & *converted;
    }

    DeleteObject(converted);

    converted = (_String *)parameters(1)->toStr();
    result = result & " from string " & *converted;
    break;

  case HY_HBL_COMMAND_FPRINTF: // print stuff to file (or stdout)

    converted = (_String *)parameters(0)->toStr();
    result = _String("fprintf(") & (*converted);
    DeleteObject(converted);

    converted = nil;
    for (k = 1; k < parameters.countitems(); k++) {
      result = result & *"," & *((_String *)parameters(k));
    }
    result = result & ")";

    break;

  case HY_HBL_COMMAND_OPTIMIZE:          // optimize the likelihood function
  case HY_HBL_COMMAND_COVARIANCE_MATRIX: // compute the covariance matrix
    converted = (_String *)parameters(0)->toStr();
    if (code == HY_HBL_COMMAND_OPTIMIZE) {
      result = _String("Optimize storing into, ");
    } else {
      result = _String("Calculate the Covariance Matrix storing into, ");
    }

    result = result & (*converted) & ", the following likelihood function:";
    DeleteObject(converted);

    converted = nil;
    for (k = 1; k < parameters.countitems(); k++) {
      result = result & *((_String *)parameters(k)) & " ; ";
    }

    break;

  case 11: // build the likelihood function

    converted = (_String *)parameters(0)->toStr();
    result = _String("Construct the following likelihood function:");
    DeleteObject(converted);

    converted = nil;
    for (k = 1; k < parameters.countitems(); k++) {
      result = result & *((_String *)parameters(k)) & " ; ";
    }

    break;

  case 12: // data set simulation
    converted = (_String *)parameters(0)->toStr();
    result = _String("Simulate Data Set ") & (*converted);
    DeleteObject(converted);
    converted = (_String *)parameters(1)->toStr();
    result = result & _String(" from the likelihood function ") & (*converted);

    break;

  case 13: // a function
    converted = (_String *)parameters(0)->toStr();
    result = _String("Function ") & (*converted);
    DeleteObject(converted);
    converted = new _String(simpleParameters(0));
    checkPointer(converted);
    result = result & _String(" of ") & (*converted) & _String(" parameters:{");
    DeleteObject(converted);
    converted = (_String *)parameters(simpleParameters(0) - 1)->toStr();
    result = result & _String("}");

    break;

  case 14: // return statement
    result = _String("A return statement with:");
    converted = new _String(simpleParameters(0));
    checkPointer(converted);
    result = result & *converted;

    break;

  case 16: // data set merger
    converted = (_String *)parameters(0)->toStr();
    result = _String("Build dataset") & (*converted) & _String(" by ");
    if (abs(simpleParameters(0) == 1)) {
      result = result & _String(" concatenating ");
    } else {
      result = result & _String(" combining ");
    }
    if (simpleParameters(0) < 0) {
      result = result & _String("(deleting arguments upon completion)");
    }
    for (k = 1; k < parameters.countitems(); k++) {
      result = result & *((_String *)parameters(k));
      if (k < parameters.countitems() - 1) {
        result = result & _String(",");
      }
    }
    break;

  case HY_HBL_COMMAND_EXPORT:
    converted = (_String *)parameters(1)->toStr();
    result = _String("Export ") & (*converted);
    DeleteObject(converted);
    converted = (_String *)parameters(0)->toStr();
    checkPointer(converted);
    result = result & _String(" to ") & *converted;
    break;

  case HY_HBL_COMMAND_MOLECULAR_CLOCK: // a call to MolecularClock
    converted = (_String *)parameters(0)->toStr();
    result = ",";
    result =
        _String("Molecular clock imposed starting at ") & (*converted) &
        " on variables " & _String((_String *)parameters.Join(&result, 1, -1));
    break;

  case 20: // category variable construction
    converted = (_String *)parameters.toStr();
    result = _String("Category variable: ") & (*converted);
    break;

  case 21: // optimize the likelihood function
           {
    converted = (_String *)parameters(0)->toStr();
    result = _String("Construct the category matrix, ") & (*converted) &
             ", for the likelihood function:";
    DeleteObject(converted);
    result = result & *((_String *)parameters(1)) & " in ";
    converted = nil;
    _String testAgainst("COMPLETE");
    if (parameters.countitems() > 2) {
      if (((_String *)parameters(2))->Equal(&testAgainst)) {
        result = result & "complete format";
        break;
      }
    }
    result = result & "short format";

    break;
  }
  case HY_HBL_COMMAND_CLEAR_CONSTRAINTS: { // clear constraints
    converted = (_String *)parameters.toStr();
    result = _String("Clear contstraints on: ") & (*converted);
    break;
  }
  case HY_HBL_COMMAND_SET_DIALOG_PROMPT: { // set dialog prompt
    converted = (_String *)parameters.toStr();
    result = _String("Set dialog prompt to: ") & (*converted);
    break;
  }
  case 24: { // select standard model
    converted = (_String *)parameters.toStr();
    result = _String("Select Template Model for ") & (*converted);
    break;
  }
  case 25:   // fscanf
  case 56: { // sscanf
    converted = (_String *)parameters(0);
    result = _String("Read the following data from ") & *converted & " : ";
    long shift = 1;
    for (long p = 0; p < simpleParameters.lLength; p++) {
      long theFormat = simpleParameters(p);
      if (theFormat < 0) {
        shift--;
        result = result & " REWIND THE FILE";
      } else {
        result = result & *(_String *)allowedFormats(theFormat) & " into " &
                 *(_String *)parameters(p + shift) & ";";
      }
    }
    converted = nil;
    break;
  }
  case 26: { // replicate constraint
    converted = (_String *)parameters(0);
    result = _String("Replicate the following constraint ") & *converted &
             " using these variables:";
    for (long p = 1; p < parameters.lLength; p++) {
      result = result & *(_String *)parameters(p) & ";";
    }
    converted = nil;
    break;
  }

  case HY_HBL_COMMAND_USE_MODEL: { // use matrix
    converted = (_String *)parameters(0);
    result = _String("Use the matrix ") & *converted;
    converted = nil;
    break;
  }

  case 31: { // define a model
    converted = (_String *)parameters(0);
    result = _String("Define the model ") & *converted;
    converted = (_String *)parameters(1);
    result = result & _String(" using transition matrix '") & *converted & "'";
    if (parameters.lLength > 2) {
      converted = (_String *)parameters(2);
      result =
          result & _String(" and equilibrium frequencies '") & *converted & "'";
    }
    converted = nil;
    break;
  }

  case 32: { // choice list
    converted = (_String *)parameters(1)->toStr();
    result = _String("Choice List for ") & *converted;
    DeleteObject(converted);
    converted = (_String *)parameters(3)->toStr();
    result = result & _String(" with choice list:") & *converted;
    DeleteObject(converted);
    converted = (_String *)parameters(0);
    result = result & _String(". Store result in ") & *converted;
    converted = nil;
    break;
  }

  case HY_HBL_COMMAND_GET_STRING: { // get string from object
    converted = (_String *)parameters(2)->toStr();
    result = _String("Get string ") & *converted;
    DeleteObject(converted);
    converted = (_String *)parameters(1)->toStr();
    result = result & _String(" from object:") & *converted;
    DeleteObject(converted);
    converted = (_String *)parameters(0);
    result = result & _String(". Store result in ") & *converted;
    converted = nil;
    break;
  }
  case HY_HBL_COMMAND_SET_PARAMETER: { // set parameter value
    converted = (_String *)parameters(1)->toStr();
    result = _String("Set parameter ") & *converted;
    DeleteObject(converted);
    converted = (_String *)parameters(0)->toStr();
    result = result & _String(" of ") & *converted;
    DeleteObject(converted);
    converted = (_String *)parameters(2);
    result = result & _String(" to ") & *converted;
    converted = nil;
    break;
  }
  case 36: { // open data panel
    converted = (_String *)parameters(0)->toStr();
    result = _String("Open data window for the data set ") & *converted;
    DeleteObject(converted);
    converted = (_String *)parameters(1)->toStr();
    result = result & _String(" list species ") & *converted;
    DeleteObject(converted);
    converted = (_String *)parameters(2);
    result = result & _String(" with window settings ") & *converted;
    converted = (_String *)parameters(3);
    result = result & _String(" with partition settings ") & *converted;
    converted = nil;
    break;
  }
  case 37: { // open data panel
    converted = (_String *)parameters(0)->toStr();
    result = _String("Get Information for ") & *converted;
    DeleteObject(converted);
    converted = (_String *)parameters(1)->toStr();
    result = result & _String(" storing in ") & *converted;
    break;
  }

  case 38: { // reconsruct ancestors
    converted = (_String *)parameters(0)->toStr();
    result = _String("Reconstruct Ancestors into ") & (*converted);
    DeleteObject(converted);
    converted = (_String *)parameters(1)->toStr();
    result = result & _String(" from the likelihood function ") & (*converted);
    break;
  }

  case 39:   // execute commands
  case 62:   // execute a file
  case 66: { // load function library
    converted = (_String *)parameters(0)->toStr();
    if (code == 39) {
      result = _String("ExecuteCommands in string ");
    } else if (code == 62) {
      result = _String("ExecuteAFile from file ");
    } else {
      result = _String("LoadFunctionLibrary from file");
    }

    result = result & *converted & _String(" using basepath ") &
             ((_String *)parameters(1)->toStr()) & '.';
    if (simpleParameters.lLength) {
      result = result & "\nCompiled.";
    } else if (parameters.lLength > 2) {
      _String inputName((_String *)parameters(2)->toStr());
      result = result & " reading input from " & inputName;
      _AssociativeList *inputValues =
          (_AssociativeList *)FetchObjectFromVariableByType(&inputName,
                                                            ASSOCIATIVE_LIST);
      if (inputValues) {
        result = result & '\n' & _String((_String *)inputValues->toStr());
      }

      if (parameters.lLength > 3) {
        result = result & " using name space prefix " &
                 _String((_String *)parameters(3)->toStr());
      }
    }
    break;
  }

  case 40: { // open window
    converted = (_String *)parameters(0)->toStr();
    result = _String("Open window of type ") & (*converted);
    DeleteObject(converted);
    converted = (_String *)parameters(1)->toStr();
    result = result & _String(" with parameters from ") & (*converted);
    break;
  }

  case 41: { // spawn LF
    converted = (_String *)parameters(0)->toStr();
    result = _String("Spawn Likelihood Function ") & (*converted);
    DeleteObject(converted);
    converted = (_String *)parameters(1)->toStr();
    result = result & _String(" with tree ") & (*converted);
    DeleteObject(converted);
    converted = (_String *)parameters(2)->toStr();
    result = result & _String(" from ") & (*converted);
    DeleteObject(converted);
    converted = (_String *)parameters(3)->toStr();
    result = result & _String(" using sequences ids in ") & (*converted);
    break;
  }

  case HY_HBL_COMMAND_DIFFERENTIATE: { // Differentiate
    converted = (_String *)parameters(1)->toStr();
    result = _String("Differentiate '") & (*converted);
    DeleteObject(converted);
    converted = (_String *)parameters(2)->toStr();
    result = result & _String("' on ") & (*converted);
    DeleteObject(converted);
    if (parameters.lLength == 4) {
      converted = (_String *)parameters(3)->toStr();
      result = result & _String(" ") & (*converted) & " times ";
      DeleteObject(converted);
    }
    converted = (_String *)parameters(2)->toStr();
    result = result & _String(" storing result in ") & (*converted);
    break;
  }
  case 43:   // Find Root
  case 48: { // Integrate
    converted = (_String *)parameters(1)->toStr();
    if (code == 43) {
      result = _String("Find the root of ") & (*converted);
    } else {
      result = _String("Integrate ") & (*converted);
    }
    DeleteObject(converted);
    converted = (_String *)parameters(2)->toStr();
    result = result & _String(" in ") & (*converted);
    DeleteObject(converted);
    converted = (_String *)parameters(3)->toStr();
    result = result & _String(" between ") & (*converted);
    DeleteObject(converted);
    converted = (_String *)parameters(4)->toStr();
    result = result & _String(" and ") & (*converted);
    DeleteObject(converted);
    converted = (_String *)parameters(0)->toStr();
    result = result & _String(" storing the result in ") & (*converted);
    break;
  }
  case 44: { // MPISend
    converted = (_String *)parameters(1)->toStr();
    result = _String("MPI Send ") & (*converted);
    DeleteObject(converted);
    if (parameters.lLength > 2) {
      converted = (_String *)parameters(2)->toStr();
      result = _String(" with input from ") & (*converted);
      DeleteObject(converted);
    }
    converted = (_String *)parameters(0)->toStr();
    result = result & _String(" to MPI Node ") & (*converted);
    break;
  }
  case 45: { // MPIReceive
    converted = (_String *)parameters(0)->toStr();
    result = _String("MPI Receive from ") & (*converted);
    DeleteObject(converted);
    converted = (_String *)parameters(1)->toStr();
    result =
        result & _String(" storing actual sender node into ") & (*converted);
    DeleteObject(converted);
    converted = (_String *)parameters(2)->toStr();
    result =
        result & _String(" and storing the string result into ") & (*converted);
    break;
  }
  case 46: { //GetDataInfo
    converted = (_String *)parameters(1)->toStr();
    if (parameters.lLength > 4) {
      result = _String("Compute pairwise differences from ") & (*converted);
    } else {
      result = _String("GetDataInfo from ") & (*converted);
    }
    DeleteObject(converted);
    if (parameters.lLength > 2) {
      converted = (_String *)parameters(2)->toStr();
      result = result & _String(" for sequence ") & (*converted);
      DeleteObject(converted);
      if (parameters.lLength > 3) {
        converted = (_String *)parameters(3)->toStr();
        if (parameters.lLength > 4) {
          result = result & _String(" and site ") & (*converted);
        } else {
          result = result & _String(" and sequence ") & (*converted);
        }
      }
      DeleteObject(converted);
    }
    converted = (_String *)parameters(0)->toStr();
    result = result & _String(" and storing the result in ") & (*converted);
    break;
  }
  case 47: { //GetDataInfo
    converted = (_String *)parameters(0)->toStr();
    result = _String("StateCounter on ") & (*converted);
    DeleteObject(converted);
    converted = (_String *)parameters(1)->toStr();
    result = result & _String(" using the callback: ") & (*converted);
    break;
  }
  case HY_HBL_COMMAND_LFCOMPUTE: { //Compute LF
    converted = (_String *)parameters(0)->toStr();
    result = blLFCompute & (*converted);
    DeleteObject(converted);
    converted = (_String *)parameters(1)->toStr();
    result = result & ',' & (*converted) & ')';
    break;
  }
  case HY_HBL_COMMAND_GET_URL: { //GetURL
    converted = (_String *)parameters(0)->toStr();
    result = blGetURL & (*converted);
    DeleteObject(converted);
    converted = (_String *)parameters(1)->toStr();
    result = result & ',' & (*converted);
    if (parameters.lLength > 2) {
      DeleteObject(converted);
      converted = (_String *)parameters(2)->toStr();
      result = result & ',' & (*converted);
    }
    result = result & ')';
    break;
  }
  case 52: { //Simulate
    converted = (_String *)parameters(0)->toStr();
    result = blDataSet & (*converted) & "=Simulate(";
    DeleteObject(converted);
    converted = (_String *)parameters(1)->toStr();
    result = result & (*converted);
    for (long i = 2; i < parameters.lLength; i++) {
      DeleteObject(converted);
      converted = (_String *)parameters(i)->toStr();
      result = result & ',' & (*converted);
    }
    result = result & ')';
    break;
  }
  case 53:   //DoSQL
  case 55:   //AlignSequences
  case 57: { //GetNeutralNull
    converted = (_String *)parameters(0)->toStr();
    result = (code == 53 ? blDoSQL
                         : (code == 55 ? blAlignSequences : blGetNeutralNull)) &
             (*converted) & '(';
    for (long i = 1; i < parameters.lLength; i++) {
      DeleteObject(converted);
      converted = (_String *)parameters(i)->toStr();
      result = result & ',' & (*converted);
    }
    result = result & ')';
    break;
  }
  case 58: {
    converted = (_String *)parameters(0)->toStr();
    result = blHBLProfile & " " & *converted;
    break;
  }
  case HY_HBL_COMMAND_DELETE_OBJECT: {
    converted = (_String *)parameters.toStr();
    result = blDeleteObject & '(' & *converted & ')';
    break;
  }
  case HY_HBL_COMMAND_REQUIRE_VERSION: {
    converted = (_String *)parameters(0)->toStr();
    result = blRequireVersion & '(' & *converted & ')';
    break;
  }
  case 61:
  case 63: {
    converted = (_String *)parameters(0)->toStr();
    result = (code == 61 ? blSCFG : blNN) & *converted & "=(";
    for (long i = 1; i < parameters.lLength; i++) {
      DeleteObject(converted);
      converted = (_String *)parameters(i)->toStr();
      if (i > 1) {
        result = result & ',' & (*converted);
      } else {
        result = result & (*converted);
      }
    }
    result = result & ')';
    break;
  }
  case 64:
    converted = (_String *)parameters(0)->toStr();
    result = blBGM & *converted & "=(";
    for (long p = 1; p < parameters.lLength; p++) {
      DeleteObject(converted);
      converted = (_String *)parameters(p)->toStr();
      if (p > 1) {
        result = result & ',' & (*converted);
      } else {
        result = result & (*converted); // first argument
      }
    }
    result = result & ')';
    break;
  case HY_HBL_COMMAND_ASSERT: {
    converted = (_String *)parameters(0)->toStr();
    result = _String("Assert ") & "'" & *converted & "'";
    break;
  }
  }

  DeleteObject(converted);
  return result.makeDynamic();
}

//______________________________________________________________________________
/*
void _ElementaryCommand::ExecuteCase0(_ExecutionList &chain) {
  chain.currentCommand++;

  if (chain.cli) {
    _Parameter result = ((_Formula *)simpleParameters.lData[1])
        ->ComputeSimple(chain.cli->stack, chain.cli->values);
    long sti = chain.cli->storeResults.lData[chain.currentCommand - 1];
    if (sti >= 0) {
      chain.cli->values[sti].value = result;
    }
    return;
  }


  ExecuteFormula((_Formula *)simpleParameters.lData[1],
                 (_Formula *)simpleParameters.lData[2],
                 simpleParameters.lData[0], simpleParameters.lData[3],
                 chain.nameSpacePrefix, simpleParameters.lData[4]);

  if (terminateExecution) {
    WarnError(_String("Problem occurred in line: ") & *this);
    return;
  }
}
*/

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase4(_ExecutionList &chain) {
  chain.currentCommand++;
  if (simpleParameters.lLength == 2) {

  }

  if (simpleParameters.lLength == 3 || parameters.lLength) {
    if (parameters.lLength && simpleParameters.lLength < 3) {
      _Formula f;
      //printf ("Namespace: %x\nCode: %s\n", chain.nameSpacePrefix,
      //((_String*)parameters(0))->sData);

      _FormulaParsingContext fpc(nil, chain.nameSpacePrefix);
      long status = Parse(&f, *(_String *)parameters(0), fpc);

      //printf ("Print formula: %s\n", _String((_String*)f.toStr()).sData);

      if (status == HY_FORMULA_EXPRESSION) {
        simpleParameters << long(f.makeDynamic());
      } else {
        return;
      }
    }

    if (chain.cli) {
      if (((_Formula *)simpleParameters(2))
              ->ComputeSimple(chain.cli->stack, chain.cli->values) == 0.0) {
        chain.currentCommand = simpleParameters.lData[1];
        return;
      }
    } else {
      _PMathObj result = ((_Formula *)simpleParameters(2))->Compute();
      if (!result) {
        WarnError("Condition Evaluation Failed");
        return;
      }

      if (terminateExecution) {
        subNumericValues = 2;
        _String *s = (_String *)((_Formula *)simpleParameters(2))->toStr();
        subNumericValues = 0;
        _String err =
            _String("Failed while evaluating: ") &
            _String((_String *)((_Formula *)simpleParameters(2))->toStr()) &
            " - " & *s;
        DeleteObject(s);
        WarnError(err);
        return;
      }

      bool conditionFalse = false;

      switch (result->ObjectClass()) {
      case NUMBER:
        conditionFalse = result->Value() == 0.0;
        break;
      case STRING:
        conditionFalse = ((_FString *)result)->IsEmpty();
        break;
      default:
        WarnError("Condition evaluation result be be a number or a string");
        return;

      }

      if (conditionFalse) {
        chain.currentCommand = simpleParameters.lData[1];
        return;
      }
    }
  }
  chain.currentCommand = simpleParameters.lData[0];

  if (chain.currentCommand == -1) {
    terminateExecution = true;
    chain.currentCommand = chain.lLength;
  }
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase5(_ExecutionList &chain) {
  chain.currentCommand++;
  FILE *df;
  _String fName = *(_String *)parameters(1);
  _DataSet *ds;

  if (simpleParameters.lLength == 1) {
    fName =
        GetStringFromFormula((_String *)parameters(1), chain.nameSpacePrefix);
    ds = ReadDataSetFile(
        nil, 0, &fName, nil,
        chain.nameSpacePrefix ? chain.nameSpacePrefix->GetName() : nil);
  } else {
    if (fName.Equal(&useNexusFileData)) {
      if (!lastNexusDataMatrix) {
        _String errMsg = useNexusFileData & " was used in ReadDataFile, and no "
                                            "NEXUS data matrix was available.";
        acknError(errMsg);
        return;
      }
      ds = lastNexusDataMatrix;
    } else {
#if defined __WINDOZE__ && !defined __HEADLESS__
      lastFileTypeSelection = 1;
#endif

      fName.ProcessFileName(false, false, (Ptr) chain.nameSpacePrefix);
      if (terminateExecution) {
        return;
      }
      SetStatusLine("Loading Data");

      df = doFileOpen(fName.getStr(), "rb");
      if (df == nil) {
        // try reading this file as a string formula
        fName = GetStringFromFormula((_String *)parameters(1),
                                     chain.nameSpacePrefix);
        fName.ProcessFileName(false, false, (Ptr) chain.nameSpacePrefix);

        if (terminateExecution) {
          return;
        }

        df = doFileOpen(fName.getStr(), "rb");
        if (df == nil) {
          _String errMsg("Could not find source dataset file:");
          errMsg = errMsg & *(_String *)parameters(1) & " Path stack: " &
                   _String((_String *)pathNames.toStr());
          WarnError(errMsg);
          return;
        }
      }
      ds = ReadDataSetFile(
          df, 0, nil, nil,
          chain.nameSpacePrefix ? chain.nameSpacePrefix->GetName() : nil);
      fclose(df);
    }
  }

  // 20110802: need to check that this data set is not empty

  if (ds->NoOfSpecies() && ds->NoOfColumns()) {
    _String *dsID =
        new _String(chain.AddNameSpaceToID(*(_String *)parameters(0)));
    StoreADataSet(ds, dsID);
    DeleteObject(dsID);
  } else {
    DeleteObject(ds);
    WarnError("The format of the sequence file has not been recognized and may "
              "be invalid");
  }

  //StoreADataSet (ds, (_String*)parameters(0));
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase11(
    _ExecutionList &chain) /*
                            code cleanup SLKP 20090316
                           */
    {
  chain.currentCommand++;

  _String parm, errMsg;

  bool
  explicitFreqs = simpleParameters.lLength,
  // if false then the likelihood function will be of the form
  // (filter1,tree1,filter2,tree2,...,filterN,treeN)
      // if true then the likelihood function will be of the form
      // (filter1,tree1,freq1,filter2,tree2,freq2,...,filterN,treeN,freqN)
      assumeList = parameters.lLength > 2;
  // if there is only one parameter to the function constructor, it is assumed
  // to be a string matrix
  // otherwise it is expected to be a collection of literals

  _List *likelihoodFunctionSpec = nil, passThisToLFConstructor;

  if (assumeList) {
    likelihoodFunctionSpec = new _List(parameters, 1, -1);
  } else {
    _Matrix *matrixOfStrings = (_Matrix *)ProcessAnArgumentByType(
        (_String *)parameters(1), chain.nameSpacePrefix, MATRIX);
    if (matrixOfStrings && matrixOfStrings->IsAStringMatrix()) {
      likelihoodFunctionSpec = new _List;
      matrixOfStrings->FillInList(*likelihoodFunctionSpec);
      if (likelihoodFunctionSpec->lLength == 0) {
        DeleteObject(likelihoodFunctionSpec);
        likelihoodFunctionSpec = nil;
      }
    }
    if (likelihoodFunctionSpec == nil) {
      WarnError(_String("Not a valid string matrix object passed to a "
                        "_LikelihoodFunction constructor: ") &
                *(_String *)parameters(1));
      return;
    }
  }

  long i = 0, stepper = explicitFreqs ? 3 : 2;

  for (; i <= (long) likelihoodFunctionSpec->lLength - stepper; i += stepper) {
    _String *dataset = (_String *)(*likelihoodFunctionSpec)(i),
            *tree = (_String *)(*likelihoodFunctionSpec)(i + 1),
            *freq = explicitFreqs ? (_String *)(*likelihoodFunctionSpec)(i + 2)
                                  : nil;

    if (FindDataSetFilterName(
            AppendContainerName(*dataset, chain.nameSpacePrefix)) != -1) {
      _TheTree *thisTree = (_TheTree *)FetchObjectFromVariableByType(
          &AppendContainerName(*tree, chain.nameSpacePrefix), TREE);
      if (thisTree) {
        _CalcNode *thisNode = thisTree->DepthWiseTraversal(true);
        if (!freq) { // no explicit frequency parameter; grab one from the tree
          long theFreqID = -1, theModelID = -1, finalFreqID = -1;
          bool done = false;

          while (1) {
            if ((theModelID = thisNode->GetModelIndex()) ==
                HY_NO_MODEL) { // this node has no model
              done = false;
              break;
            }
            theFreqID = modelFrequenciesIndices.lData[theModelID];
            thisNode = thisTree->DepthWiseTraversal();
            while (thisNode) {
              theModelID = thisNode->GetModelIndex();
              if (theModelID == HY_NO_MODEL) { // no model
                done = false;
                break;
              }
              if (modelFrequenciesIndices.lData[theModelID] != theFreqID) {
                done = true;
                break;
              }
              if (thisTree->IsCurrentNodeTheRoot()) {
                break;
              }
              thisNode = thisTree->DepthWiseTraversal();
            }
            if (theFreqID < 0) {
              finalFreqID = -theFreqID - 1;
            } else {
              finalFreqID = theFreqID;
            }
            break;
          }

          if (finalFreqID >= 0) {
            _String freqID =
                chain.TrimNameSpaceFromID(*LocateVar(finalFreqID)->GetName());
            passThisToLFConstructor &&dataset;
            passThisToLFConstructor &&tree;
            passThisToLFConstructor &&freqID;
            continue;
          } else {
            if (!done) {
              errMsg = (((_String)(
                  "LF: Not a well-defined tree/model combination: ") & *tree));
            } else {
              errMsg = (((_String)("LF: All models in the tree: ") & *tree &
                         _String(" must share the same frequencies vector")));
            }
          }
        } else {
          if (FetchObjectFromVariableByType(
                  &AppendContainerName(*freq, chain.nameSpacePrefix), MATRIX)) {
            passThisToLFConstructor &&dataset;
            passThisToLFConstructor &&tree;
            passThisToLFConstructor &&freq;
            continue;
          }
          errMsg =
              (((_String)("LF: Not a valid frequency matrix ID: ") & *freq));
        }
      } else {
        errMsg = (((_String)("LF: Not a valid tree ID: ") & *tree));
      }

    } else {
      errMsg = (((_String)("LF: Not a valid dataset filter: ") & *dataset));
    }

    if (errMsg.sLength) {
      DeleteObject(likelihoodFunctionSpec);
      WarnError(errMsg);
      return;
    }
  }

  if (i == likelihoodFunctionSpec->lLength - 1) { // computing template
    passThisToLFConstructor &&*((_String *)(*likelihoodFunctionSpec)(i));
  }

  DeleteObject(likelihoodFunctionSpec);

  _String lfID = chain.AddNameSpaceToID(
      *(_String *)parameters(0)); // the ID of the likelihood function
  long likeFuncObjectID = FindLikeFuncName(lfID);
  if (likeFuncObjectID == HY_NOT_FOUND)
      // not an existing LF ID
      {
    _LikelihoodFunction *lkf = new _LikelihoodFunction();
    if (!lkf->Construct(passThisToLFConstructor, chain.nameSpacePrefix))
        // constructor failed
        {
      DeleteObject(lkf);
    } else {
      likeFuncObjectID = likeFuncNamesList.Find(&empty);
      // see if there are any vacated spots in the list

      if (likeFuncObjectID < 0) {
        likeFuncList << lkf;
        likeFuncNamesList &&(&lfID);
        DeleteObject(lkf);
      } else {
        likeFuncNamesList.Replace(likeFuncObjectID, &lfID, true);
        likeFuncList.lData[likeFuncObjectID] = (long) lkf;
      }
    }
  } else {
    _LikelihoodFunction *lkf =
        (_LikelihoodFunction *)likeFuncList(likeFuncObjectID);
    if (!lkf->Construct(passThisToLFConstructor, chain.nameSpacePrefix)) {
      KillLFRecord(likeFuncObjectID, false);
    }
  }

}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase12(_ExecutionList &chain) {
  chain.currentCommand++;
  SetStatusLine("Simulating Data");

  _String likefID = chain.AddNameSpaceToID(*(_String *)parameters(1)),
          tempString = ProcessStringArgument(&likefID), errMsg;

  if (tempString.sLength) {
    likefID = tempString;
  }

  long f = FindLikeFuncName(likefID), s2 = FindSCFGName(likefID);

  if (f == -1 && s2 == -1) {
    WarnError(_String("Likelihood Function (or SCFG)") & likefID &
              " has not been initialized");
    return;
  }

  if (f >= 0) {
    _DataSet *ds = new _DataSet;
    checkPointer(ds);

    _List theExclusions;

    if (parameters.lLength > 2) // there is a list of exclusions there
        // ';'-sep for different partititons
        // ','-sep for states in a given partition
        {
      // SLKP mod 20070622 to allow string expressions as well
      _String theExc(ProcessLiteralArgument((_String *)parameters(2),
                                            chain.nameSpacePrefix));
      if (theExc.sLength) {
        long f = theExc.Find(';'), g = 0;

        while (1) {
          _String subExc(theExc, g, (f == -1) ? (-1) : (f - 1));
          long h = subExc.Find(','), l = 0;
          _List myExc;

          while (1) {
            _String excludeMe(subExc, l, (h == -1) ? (-1) : (h - 1));
            myExc &&&excludeMe;
            if (h == -1) {
              break;
            }
            l = h + 1;
            h = subExc.Find(',', h + 1, -1);
          }
          theExclusions &&&myExc;
          if (f == -1) {
            break;
          }
          g = f + 1;
          f = theExc.Find(';', f + 1, -1);
        }
      }

    }

    _Matrix *catValues = nil, *catNames = nil;

    _Variable *catValVar = nil, *catNameVar = nil;

    if (parameters.lLength > 3)
        // a matrix to store simulated category values
        {
      _String matrixName(chain.AddNameSpaceToID(*(_String *)parameters(3)));
      if (!(catValVar =
                CheckReceptacle(&matrixName, blSimulateDataSet, true))) {
        return;
      } else {
        checkPointer(catValues = new _Matrix(1, 1, false, true));
      }
    }

    if (parameters.lLength > 4)
        // a matrix to store simulated category values
        {
      _String matrixName(chain.AddNameSpaceToID(*(_String *)parameters(4)));
      if (!(catNameVar =
                CheckReceptacle(&matrixName, blSimulateDataSet, true))) {
        return;
      } else {
        checkPointer(catNames = new _Matrix(1, 1, false, true));
      }
    }

    _String *resultingDSName =
        new _String(chain.AddNameSpaceToID(*(_String *)parameters(0)));

    if (!resultingDSName->IsValidIdentifier(true)) {
      errMsg = *resultingDSName &
               " is not a valid receptacle identifier in call to " &
               blSimulateDataSet;
      DeleteObject(resultingDSName);
      WarnError(errMsg);
      return;
    }

    ((_LikelihoodFunction *)likeFuncList(f))
        ->Simulate(*ds, theExclusions, catValues, catNames);

    if (catValues) {
      catValVar->SetValue(catValues, false);
    }
    if (catNames) {
      catNameVar->SetValue(catNames, false);
    }

    StoreADataSet(ds, resultingDSName);
    DeleteObject(resultingDSName);
  } else {
    _String newCorpus = chain.AddNameSpaceToID(*(_String *)parameters(0));
    CheckReceptacleAndStore(
        &newCorpus, " SimulateDataSet (SCFG)", true,
        new _FString(((Scfg *)scfgList(s2))->SpawnRandomString()), false);
  }
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase38(_ExecutionList &chain, bool sample) {
  chain.currentCommand++;
  SetStatusLine("Reconstructing Ancestors");

  _String *likef = (_String *)parameters(1),
          tempString = ProcessStringArgument(likef), errMsg;

  if (tempString.sLength) {
    likef = &tempString;
  }

  _String name2lookup = AppendContainerName(*likef, chain.nameSpacePrefix);
  long objectID = FindLikeFuncName(name2lookup);
  if (objectID >= 0) {
    _DataSet *ds = (_DataSet *)checkPointer(new _DataSet);
    _String *dsName = new _String(
        AppendContainerName(*(_String *)parameters(0), chain.nameSpacePrefix));
    _LikelihoodFunction *lf = ((_LikelihoodFunction *)likeFuncList(objectID));

    _Matrix *partitionList = nil;
    if (parameters.lLength > 2) {
      _String secondArg = *(_String *)parameters(2);
      partitionList = (_Matrix *)ProcessAnArgumentByType(
          &secondArg, chain.nameSpacePrefix, MATRIX);
    }
    _SimpleList partsToDo;
    if (lf->ProcessPartitionList(partsToDo, partitionList,
                                 " ancestral reconstruction")) {
      lf->ReconstructAncestors(*ds, partsToDo, *dsName, sample,
                               simpleParameters.Find(-1) >= 0,
                               simpleParameters.Find(-2) >= 0);
    }
    StoreADataSet(ds, dsName);
    DeleteObject(dsName);
  } else {
    objectID = FindSCFGName(name2lookup);
    if (objectID >= 0) /* reconstruct best parse tree for corpus using SCFG */
        {
      CheckReceptacleAndStore(
          &AppendContainerName(*(_String *)parameters(0),
                               chain.nameSpacePrefix),
          " ReconstructAncestors (SCFG)", true,
          new _FString(((Scfg *)scfgList(objectID))->BestParseTree()), false);
    } else {
      errMsg = (((_String)("Likelihood Function/SCFG") & *likef &
                 _String(" has not been initialized")));
      WarnError(errMsg);
      return;
    }
  }
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase39(_ExecutionList &chain) {
  chain.currentCommand++;

  _String *commands, theCommand, *namespc = nil;

  if (code == 39) {
    commands = ProcessCommandArgument((_String *)parameters(0));
  } else {
    _String filePath = GetStringFromFormula((_String *)parameters(0),
                                            chain.nameSpacePrefix),
            originalPath = filePath;

    FILE *commandSource = nil;
    if (code == 66) {
      bool hasExtension = filePath.FindBackwards(".", 0, -1) > 0;

      for (unsigned long p = 0;
           !commandSource && p < standardLibraryPaths.lLength; p++) {
        for (unsigned long e = 0;
             !commandSource && e < standardLibraryExtensions.lLength; e++) {
          _String tryPath = *((_String *)standardLibraryPaths(p)) & filePath &
                            *((_String *)standardLibraryExtensions(e));

          // printf ("%s\n", tryPath.sData);
          _Parameter reload = 0.;
          checkParameter(alwaysReloadLibraries, reload, 0.);

          if (loadedLibraryPaths.Find(&tryPath) >= 0 &&
              parameters.lLength == 2 && reload < 0.5) {
            ReportWarning(_String("Already loaded '") & originalPath &
                          "' from " & tryPath);
            return;
          }
          if ((commandSource = doFileOpen(tryPath.getStr(), "rb"))) {
            filePath = tryPath;
            break;
          }
          if (hasExtension) {
            break;
          }
        }
      }

    }

    if (commandSource == nil) {
      filePath.ProcessFileName(false, false, (Ptr) chain.nameSpacePrefix);
      if ((commandSource = doFileOpen(filePath.getStr(), "rb")) == nil) {
        WarnError(_String(
            "Could not read command file in ExecuteAFile.\nOriginal path: '") &
                  originalPath & "'.\nExpanded path: '" & filePath & "'");
        return;
      }
    }

    if (code == 66 && commandSource) {
      ReportWarning(_String("Loaded '") & originalPath & "' from " & filePath);
      loadedLibraryPaths.Insert(filePath.makeDynamic(), 0, false, true);
    }

    commands = new _String(commandSource);
    if (fclose(commandSource)) { // failed to fclose
      DeleteObject(commands);
      WarnError(_String("Internal error: failed in a call to fclose ") &
                filePath);
    }
    PushFilePath(filePath);
  }

  if (!commands) {
    return;
  }

  if (code == 39) {
    theCommand = ProcessLiteralArgument(commands, chain.nameSpacePrefix);
  } else {
    theCommand = commands;
  }

  if (theCommand.sLength == 0) {
    WarnError(_String("Invalid string argument '") & *commands &
              "' in call to ExecuteCommands/ExecuteAFile.");
    return;
  }

  if (code == 39 && ((_String *)parameters(1))->sLength) {
    pathNames << (_String *)parameters(1);
  }

  _AVLListXL *inArg = nil;
  _List *inArgAux = nil;

  if (parameters.lLength >= 3)
      // stdin redirect (and/or name space prefix)
      {
    _PMathObj inAVL = ProcessDictionaryArgument((_String *)parameters(2),
                                                chain.nameSpacePrefix);

    if (!inAVL) {
      if (parameters.lLength == 3) {
        WarnError(
            _String("Not a valid associative array index passed as input "
                    "redirect argument to ExecuteCommands/ExecuteAFile: )") &
            *(_String *)parameters(2));
        return;
      }
    } else {
      _AssociativeList *stdinRedirect = (_AssociativeList *)inAVL;

      checkPointer(inArgAux = new _List);
      checkPointer(inArg = new _AVLListXL(inArgAux));

      _List *stdKeys = stdinRedirect->GetKeys();

      for (long kid = 0; kid < stdKeys->lLength; kid++) {
        _String *aKey = (_String *)(*stdKeys)(kid);
        if (aKey) {
          _FString *aString =
              (_FString *)stdinRedirect->GetByKey(*aKey, STRING);
          if (!aString) {
            WarnError(_String(
                "All entries in the associative array used as input redirect "
                "argument to ExecuteCommands/ExecuteAFile must be strings. The "
                "following key was not: ") & *aKey);
            DeleteObject(inAVL);
            return;
          }
          inArg->Insert(aKey->makeDynamic(),
                        (long) new _String(*aString->theString), false);
        }
      }
    }

    DeleteObject(inAVL);

    if (parameters.lLength > 3) {
      _String nameSpaceID = ProcessLiteralArgument((_String *)parameters(3),
                                                   chain.nameSpacePrefix);
      if (!nameSpaceID.IsValidIdentifier(true)) {
        WarnError(_String(
            "Invalid namespace ID in call to ExecuteCommands/ExecuteAFile: ") &
                  *(_String *)parameters(3));
        return;
      }
      namespc = new _String(nameSpaceID);
    }
  }

  if (parameters.lLength < 4 && chain.nameSpacePrefix) {
    namespc = new _String(*chain.nameSpacePrefix->GetName());
  }

  if (theCommand.beginswith("#NEXUS")) {
    ReadDataSetFile(nil, 1, &theCommand, nil, namespc);
  } else {
    bool result = false;
    _ExecutionList exc(theCommand, namespc, false, &result);

    if (!result) {
      chain.ReportAnExecutionError("Encountered an error while parsing HBL",
                                   false, true);
    } else {

      exc.stdinRedirectAux = inArgAux ? inArgAux : chain.stdinRedirectAux;
      exc.stdinRedirect = inArg ? inArg : chain.stdinRedirect;

      if (simpleParameters.lLength && exc.TryToMakeSimple()) {
        ReportWarning("Successfully compiled an execution list.");
        exc.ExecuteSimple();
      } else {
        exc.Execute();
      }

      exc.stdinRedirectAux = nil;
      exc.stdinRedirect = nil;
      if (exc.result) {
        DeleteObject(chain.result);
        chain.result = exc.result;
        exc.result = nil;
      }
    }
  }

  if (inArg) {
    DeleteObject(inArg);
    DeleteObject(inArgAux);
  }

  DeleteObject(namespc);

  if (code == 62) {
    PopFilePath();
  } else if (((_String *)parameters(1))->sLength) {
    pathNames.Delete(pathNames.lLength - 1);
  }
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase40(_ExecutionList &chain) {
  chain.currentCommand++;
  _String errMsg;

#if !defined __UNIX__ && !defined __HEADLESS__

  bool done = false;

  _String *windowType = ((_String *)parameters(0)), *positionString = nil;

  _HYWindow *resultWindow = nil;

  long positionSpec[] = { -1, -1, -1, -1 };

  if (parameters.lLength == 3) {
    positionString = ((_String *)parameters(2));
    _List *positionList = positionString->Tokenize(";");

    if ((positionList->lLength >= 2) && (positionList->lLength <= 4)) {
      _HYRect screenDim = GetScreenDimensions();
      setParameter(screenWidthVar, screenDim.right);
      setParameter(screenHeightVar, screenDim.bottom);

      for (long i = 0; i < positionList->lLength; i++) {
        _String *thisSpec = (_String *)(*positionList)(i);
        if (thisSpec->sLength) {
          _Formula fla(*thisSpec, nil, false);
          _PMathObj flav = fla.Compute();
          if (flav && (flav->ObjectClass() == NUMBER)) {
            positionSpec[i] = flav->Value();
          }
        }
      }
    }
    DeleteObject(positionList);
  }

  if (windowType->Equal(&windowTypeClose)) {
    _String wName =
        GetStringFromFormula((_String *)parameters(1), chain.nameSpacePrefix);
    long wID = FindWindowByName(wName);
    if (wID >= 0) {
      postWindowCloseEvent(((_HYWindow *)windowObjectRefs(wID))->GetID());
    }
    done = true;
  } else if (windowType->Equal(&windowTypeTree) ||
             windowType->Equal(&windowTypeTable) ||
             windowType->Equal(&windowTypeDistribTable) ||
             windowType->Equal(&windowTypeDatabase)) {
    _String *windowOptions = ((_String *)parameters(1));
    _Matrix *theOptions = (_Matrix *)FetchObjectFromVariableByType(
        &AppendContainerName(*windowOptions, chain.nameSpacePrefix), MATRIX);
    if (!theOptions) {
      if ((windowOptions->sLength > 3) && (windowOptions->sData[0] == '{') &&
          (windowOptions->sData[windowOptions->sLength - 1] == '}')) {
        theOptions = new _Matrix(*windowOptions);
      } else {
        WarnError(*windowOptions &
                  " is not a valid inline matrix specification/variable "
                  "reference in call to OpenWindow.");
        return;
      }
    } else {
      theOptions->nInstances++;
    }

    if (theOptions->MatrixType() != 2 || theOptions->GetVDim() > 1 ||
        theOptions->GetHDim() < 1) {
      DeleteObject(theOptions);
      WarnError(*windowOptions &
                " is not a valid options matrix in call to OpenWindow.");
      return;
    }

    if (windowType->Equal(&windowTypeTree)) {
      // row 1: a tree var reference
      _PMathObj arg = theOptions->GetFormula(0, 0)->Compute();
      if (!(arg && (arg->ObjectClass() == STRING))) {
        DeleteObject(theOptions);
        errMsg = "The first entry in the specification matrix must be a tree "
                 "variable ID string in call to OpenWindow.";
        WarnError(errMsg);
        return;
      }

      windowOptions = ((_FString *)arg)->theString;

      long f = LocateVarByName(
          AppendContainerName(*windowOptions, chain.nameSpacePrefix));
      if (f >= 0) {
        _Variable *treeVar = FetchVar(f);
        if (treeVar->ObjectClass() == TREE) {
          _String windowName = _String("Tree ") & *treeVar->GetName();
          long fw = FindWindowByName(windowName);
          _HYTreePanel *myTP;
          if (fw >= 0) {
            myTP = (_HYTreePanel *)windowObjectRefs(fw);
            myTP->SetVariableReference(*treeVar->GetName());
            myTP->BuildTree(true);
          } else {
            myTP = new _HYTreePanel(*treeVar->GetName(), *treeVar->GetName());
            checkPointer(myTP);
          }
          resultWindow = myTP;
          for (long k = 1; k < theOptions->GetHDim(); k++) {
            _PMathObj arg = theOptions->GetFormula(k, 0)->Compute();
            if (!(arg && (arg->ObjectClass() == STRING))) {
              DeleteObject(theOptions);
              errMsg = "All entries in the specification matrix must be "
                       "strings in call to OpenWindow.";
              WarnError(errMsg);
              return;
            }
            windowOptions = ((_FString *)arg)->theString;
            if (windowOptions->sLength) {
              switch (k) {
              case 1: { // view options
                unsigned int treeFlags = windowOptions->toNum();
                if ((treeFlags & HY_TREEPANEL_SCALE_TO_WINDOW) !=
                    (myTP->GetFlags() & HY_TREEPANEL_SCALE_TO_WINDOW)) {
                  myTP->ToggleScaleOption();
                }
                myTP->SetFlags(treeFlags);
                break;
              }
              case 2: { // window position spec
                PositionWindow(myTP, windowOptions);
                break;
              }
              case 3: { // scale variable
                myTP->SetScaleVariable(*windowOptions);
                break;
              }
              case 4: { // select branches
                _List *nodes2Select = windowOptions->Tokenize(",");
                if (nodes2Select->lLength) {
                  myTP->SelectRangeAndScroll(*nodes2Select);
                }
                DeleteObject(nodes2Select);
                break;
              }
              }
            }
          }
        } else {
          f = -1;
        }
      }
      if (f < 0) {
        errMsg = *windowOptions &
                 " is not the ID of an existing tree in call to OpenWindow.";
        WarnError(errMsg);
      }

      DeleteObject(theOptions);
      done = true;
    } else if (windowType->Equal(&windowTypeDatabase)) {
      _PMathObj arg = theOptions->GetFormula(0, 0)->Compute();
      if (!(arg && (arg->ObjectClass() == STRING))) {
        DeleteObject(theOptions);
        errMsg = "The first entry in the specification matrix must be a file "
                 "name string in call to OpenWindow.";
        WarnError(errMsg);
        return;
      }

      windowOptions = ((_FString *)arg)->theString;

      _String windowName = _String("SQLite DB: ") & *windowOptions;
      long fw = FindWindowByName(windowName);
      _HYDBWindow *myDW;
      if (fw >= 0) {
        myDW = (_HYDBWindow *)windowObjectRefs(fw);
        myDW->SetDB(windowOptions);
      } else {
        myDW = new _HYDBWindow(windowName, windowOptions);
        checkPointer(myDW);
      }
      resultWindow = myDW;

      DeleteObject(theOptions);
      done = true;
    } else {
      bool isD = windowType->Equal(&windowTypeDistribTable);

      if (theOptions->GetHDim() < 3) {
        DeleteObject(theOptions);

        errMsg = "The matrix specification for a chart window must have at "
                 "least 3 entries in call to OpenWindow.";
        WarnError(errMsg);
        return;

      }
      _PMathObj arg = theOptions->GetFormula(1, 0)->Compute(),
                arg2 = theOptions->GetFormula(2, 0)->Compute(),
                arg3 = theOptions->GetFormula(0, 0)->Compute();

      if (!(arg && (arg->ObjectClass() == STRING) && arg2 &&
            (arg2->ObjectClass() == STRING) && arg3 &&
            (arg3->ObjectClass() == STRING))) {
        DeleteObject(theOptions);

        errMsg =
            "The first two entries in the specification matrix must be strings "
            "with matrix ID (e.g. \"matrixID\") in call to OpenWindow.";
        WarnError(errMsg);
        return;
      }

      windowOptions = ((_FString *)arg)->theString;

      _String *options2 = ((_FString *)arg2)->theString;
      long f = LocateVarByName(
               AppendContainerName(*windowOptions, chain.nameSpacePrefix)),
           f2 = LocateVarByName(
               AppendContainerName(*options2, chain.nameSpacePrefix));

      if ((f >= 0) && (f2 >= 0)) {
        _Variable *labelMatrix = FetchVar(f), *dataMatrix = FetchVar(f2);

        if ((labelMatrix->ObjectClass() == MATRIX) &&
            (dataMatrix->ObjectClass() == MATRIX)) {
          _Matrix *lMatrix = (_Matrix *)labelMatrix->GetValue(),
                  *dMatrix = (_Matrix *)dataMatrix->GetValue();

          _String windowName = *((_FString *)arg3)->theString;
          long fw = FindWindowByName(windowName);
          _List columnHeaders;

          for (f2 = 0; f2 < lMatrix->GetVDim(); f2++) {
            _Formula *fff = lMatrix->GetFormula(0, f2);
            if (fff) {
              _PMathObj arg = fff->Compute();
              if (arg->ObjectClass() == STRING) {
                columnHeaders &&((_FString *)arg)->theString;
                continue;
              }
            }
            columnHeaders &&&empty;

          }

          if ((dMatrix->GetVDim() != columnHeaders.lLength) &&
              (dMatrix->GetVDim() != columnHeaders.lLength - 1)) {
            errMsg = "The number of columns in the data matrix must match the "
                     "dimension of the header matrix in call to OpenWindow "
                     "(CHART_WINDOW).";
            WarnError(errMsg);
            return;
          }

          if ((dMatrix->GetHDim() == 0) || (lMatrix->GetVDim() == 0)) {
            errMsg = "Empty data matrix (or label matrix) in call to "
                     "OpenWindow (CHART_WINDOW).";
            WarnError(errMsg);
            return;
          }

          _String mL[14];

          for (long k = 3; k < theOptions->GetHDim(); k++) {
            _PMathObj arg = theOptions->GetFormula(k, 0)->Compute();
            if (!(arg && (arg->ObjectClass() == STRING))) {
              DeleteObject(theOptions);
              errMsg = "All entries in the specification matrix must be "
                       "strings in call to OpenWindow.";
              WarnError(errMsg);
              return;
            }
            windowOptions = ((_FString *)arg)->theString;
            mL[k - 3] = *windowOptions;

            if (k == 16) {
              break;
            }
          }

          _HYChartWindow *myTP;

          if (isD) {
            _List *varInfo = mL[13].Tokenize(";"), cInfo, dInfo;

            bool firstDerived = false;
            _String errMsg;

            for (long k = 0; k < varInfo->lLength; k++) {
              _List *item = ((_String *)(*varInfo)(k))->Tokenize(":");
              if (item->lLength == 1) {
                if (cInfo.lLength) {
                  dInfo << (*item)(0);
                } else {
                  errMsg = "Derived variables appear before the atom variables";
                }
                firstDerived = true;
              } else {
                if (firstDerived) {
                  errMsg = "Some atom variables appear after derived variables";
                } else if ((item->lLength >= 3) && (item->lLength % 2 == 1)) {
                  _String *vName = (_String *)(*item)(0);
                  if (vName->IsValidIdentifier()) {
                    long catCount = (item->lLength - 1) / 2;

                    _Matrix wts(1, catCount, false, true),
                        prb(1, catCount, false, true);

                    for (long k2 = 0, k3 = catCount + 1; k2 < catCount;
                         k2++, k3++) {
                      wts.theData[k2] = ((_String *)(*item)(k2 + 1))->toNum();
                      prb.theData[k2] = ((_String *)(*item)(k3))->toNum();
                    }

                    _List anItem;
                    anItem << vName;
                    anItem &&&wts;
                    anItem &&&prb;
                    cInfo &&&anItem;
                  } else {
                    errMsg = *vName & " is an invalid atom variable identifier";
                  }
                } else {
                  errMsg = "Invalid item count in atom variable specification";
                }
              }

              if (errMsg.sLength) {
                DeleteObject(item);
                DeleteObject(varInfo);
                errMsg = errMsg & " in call to OpenWindow";
                ProblemReport(errMsg);
                return;
              }
            }

            DeleteObject(varInfo);

            if (fw >= 0) {
              myTP = (_HYChartWindow *)windowObjectRefs(fw);
              myTP->SetTable(columnHeaders, *dMatrix);
              ((_HYDistributionChartWindow *)myTP)->SetAtoms(*dMatrix, cInfo);
            } else {
              myTP = new _HYDistributionChartWindow(windowName, columnHeaders,
                                                    *dMatrix, cInfo, nil);
              checkPointer(myTP);
            }

            for (long dc = 0; dc < dInfo.lLength; dc++) {
              ((_HYDistributionChartWindow *)myTP)
                  ->AddVariable((_String *)dInfo(dc));
            }
          } else {
            if (fw >= 0) {
              myTP = (_HYChartWindow *)windowObjectRefs(fw);
              myTP->SetTable(columnHeaders, *dMatrix);
            } else {
              myTP =
                  new _HYChartWindow(windowName, columnHeaders, *dMatrix, nil);
              checkPointer(myTP);
            }
          }

          resultWindow = myTP;

          myTP->SetChartType(mL[0], mL[1], mL[2], false);
          myTP->ToggleSuspend(true);

          long kk;

          _Parameter uMin = 0.0, uMax = uMin;

          _List *ubounds = mL[12].Tokenize(",");
          if (ubounds->lLength == 3) {
            mL[12] = *(_String *)(*ubounds)(0);
            uMin = ((_String *)(*ubounds)(1))->toNum();
            uMax = ((_String *)(*ubounds)(2))->toNum();
          }
          DeleteObject(ubounds);

          if ((kk = mL[8].Find(';')) > 0) {
            myTP->SetLabels(mL[3], mL[4], mL[5], mL[6].toNum(), mL[7],
                            mL[8].Cut(0, kk - 1).toNum() - 1,
                            mL[8].Cut(kk + 1, -1).toNum() - 1, mL[12].toNum(),
                            uMin, uMax);
          } else {
            myTP->SetLabels(mL[3], mL[4], mL[5], mL[6].toNum(), mL[7], -1, -1,
                            mL[12].toNum(), uMin, uMax);
          }

          _List *optList = mL[9].Tokenize(";");

          if (optList->lLength == 3) {
            myTP->SetProjection(((_String *)(*optList)(0))->toNum(),
                                ((_String *)(*optList)(1))->toNum(),
                                ((_String *)(*optList)(2))->toNum());
          }

          DeleteObject(optList);

          optList = mL[10].Tokenize(";");
          myTP->SetFonts(optList);
          DeleteObject(optList);

          optList = mL[11].Tokenize(";");
          myTP->SetColors(optList);
          DeleteObject(optList);

          myTP->ToggleSuspend(false);
          myTP->DrawChart();
        } else {
          f = -1;
        }
      }

      if (f < 0 || f2 < 0) {
        errMsg = *windowOptions & " and " & *options2 &
                 " must both refer to existing matrices in call to OpenWindow.";
        WarnError(errMsg);
      }

      DeleteObject(theOptions);
      done = true;
    }
  }

  if (resultWindow) {
    if ((positionSpec[0] > -1) && (positionSpec[1] > -1)) {
      resultWindow->SetWindowRectangle(0, 0, positionSpec[1], positionSpec[0]);
    }
    if ((positionSpec[2] > -1) && (positionSpec[3] > -1)) {
      resultWindow->SetPosition(positionSpec[2], positionSpec[3]);
    }

    resultWindow->BringToFront();
    handleGUI(true);
  }

  if (!done) {
    errMsg = *windowType & " is not a valid window type in call to OpenWindow.";
    WarnError(errMsg);
  }
#endif
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase41(_ExecutionList &chain) {
  chain.currentCommand++;

#if !defined __UNIX__ && !defined __HEADLESS__
  _String lfID, dsWindow, errMsg;

  long treeID = -1;

  _SimpleList subsetSpec;

  lfID =
      ProcessLiteralArgument((_String *)parameters(0), chain.nameSpacePrefix);

  if (!lfID.IsValidIdentifier()) {
    errMsg = lfID & " is not a valid likelihood function identifier.";
    WarnError(errMsg);
    return;
  }

  dsWindow =
      ProcessLiteralArgument((_String *)parameters(1), chain.nameSpacePrefix);

  treeID = LocateVarByName(dsWindow);

  if (treeID < 0 || FetchVar(treeID)->ObjectClass() != TREE) {
    errMsg = dsWindow & " is not the name of an existing tree.";
    WarnError(errMsg);
    return;
  } else {
    treeID = variableNames.GetXtra(treeID);
  }

  dsWindow =
      ProcessLiteralArgument((_String *)parameters(2), chain.nameSpacePrefix);
  long k = FindWindowByName(dsWindow);

  if ((k < 0) || (((_HYWindow *)windowObjectRefs(k))->WindowKind() !=
                  HY_WINDOW_KIND_DATAPANEL)) {
    errMsg = dsWindow & " is not the name of an open data panel window.";
    WarnError(errMsg);
    return;
  }

  _HYDataPanel *theDP = (_HYDataPanel *)windowObjectRefs(k);
  _Matrix *theMatrix = (_Matrix *)FetchObjectFromVariableByType(
      &AppendContainerName(*(_String *)parameters(3), chain.nameSpacePrefix),
      MATRIX);

  if (theMatrix) {
    for (long i1 = 0; i1 < theMatrix->GetHDim(); i1++)
      for (long i2 = 0; i2 < theMatrix->GetVDim(); i2++) {
        subsetSpec << (*theMatrix)(i1, i2);
      }

    theDP->BuildLikelihoodFunction(&lfID, &subsetSpec, treeID);
    return;
  }

  errMsg =
      (*(_String *)parameters(3)) & " is not the name of an existing matrix.";
  WarnError(errMsg);

#endif
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase25(_ExecutionList &chain, bool issscanf) {
  chain.currentCommand++;
  // first of all obtain the string to be parsed
  // either read the file into a string or get a string from standard input
  _String currentParameter = *(_String *)parameters(0), *data = nil;

  long p, p2 = 0, r, q, v, t, shifter = simpleParameters.lData[0] < 0;

  bool skipDataDelete = false;
  _Variable *iseof = CheckReceptacle(&hasEndBeenReached, empty, false);

  if (currentParameter == _String("stdin")) { //
    if (chain.stdinRedirect) {
      data = chain.FetchFromStdinRedirect();
    } else {
      if (!CheckEqual(iseof->Compute()->Value(), 0) &&
          currentParameter.Equal(&scanfLastFilePath)) {
        WarnError("Ran out of standard input\n");
        return;
      }
      checkPointer((Ptr)(data = new _String(StringFromConsole())));
    }
  } else {
    if (issscanf) {
      currentParameter = chain.AddNameSpaceToID(currentParameter);
      _FString *sscanfData =
          (_FString *)FetchObjectFromVariableByType(&currentParameter, STRING);
      if (!sscanfData) {
        WarnError(currentParameter &
                  " does not refer to a string variable in call to sscanf");
        return;
      }
      data = sscanfData->theString;
      skipDataDelete = true;

      if (iseof->Compute()->Value() > 0.) {
        scanfLastFilePath = empty;
      }

      if (!currentParameter.Equal(&scanfLastFilePath) || shifter) {
        scanfLastFilePath = currentParameter;
        p = scanfLastReadPosition = 0;
      } else {
        p = p2 = scanfLastReadPosition;
        if (p >= data->sLength) {
          iseof->SetValue(new _Constant(1.0), false);
          return;
        }
      }
    } else {
      FILE *inputBuffer;
      if (currentParameter.Find('"') == -1) {
        currentParameter =
            GetStringFromFormula(&currentParameter, chain.nameSpacePrefix);
      }

      currentParameter.ProcessFileName(false, false,
                                       (Ptr) chain.nameSpacePrefix);
      if (terminateExecution) {
        return;
      }
      inputBuffer = doFileOpen(currentParameter.getStr(), "rb");
      if (!inputBuffer) {
        WarnError(currentParameter &
                  " could not be opened for reading by fscanf. Path stack: " &
                  _String((_String *)pathNames.toStr()));
        return;
      }

      if (iseof->Compute()->Value() > 0) {
        scanfLastFilePath = empty;
      }

      if (!currentParameter.Equal(&scanfLastFilePath) || shifter) {
        scanfLastFilePath = currentParameter;
        scanfLastReadPosition = 0;
      }

      fseek(inputBuffer, 0, SEEK_END);
      p = ftell(inputBuffer);
      p -= scanfLastReadPosition;

      if (p <= 0) {
        iseof->SetValue(new _Constant(1.0), false);
        fclose(inputBuffer);
        return;
      }

      data = (_String *)checkPointer(new _String((unsigned long) p));
      rewind(inputBuffer);
      fseek(inputBuffer, scanfLastReadPosition, SEEK_SET);
      fread(data->sData, 1, p, inputBuffer);
      fclose(inputBuffer);
    }
  }
  // now that the string has been read in, read in all the terms, ignoring all
  // the characters in between
  if (!skipDataDelete) {
    p = 0; // will be used to keep track of the position in the string
  }

  q = 0;
  r = shifter;

  while (r < simpleParameters.lLength && p < data->sLength) {
    _String *currentParameter = ProcessCommandArgument(
        (_String *)parameters(r + 1 - shifter)); // name of the receptacle
    if (!currentParameter) {
      DeleteObject(data);
      return;
    }
    if (!currentParameter->IsValidIdentifier()) {
      WarnError(_String('\\') & *currentParameter &
                "\" is not a valid identifier in call to fscanf.");
      DeleteObject(data);
      return;
    }
    _String namespacedParameter(chain.AddNameSpaceToID(*currentParameter));

    v = LocateVarByName(namespacedParameter);
    if (v < 0) {
      if (simpleParameters.lData[r] != 2) {
        v = CheckReceptacle(&namespacedParameter, empty, false)->GetAVariable();
      }
    } else {
      if (simpleParameters.lData[r] == 2)
        if (FetchVar(v)->ObjectClass() == TREE) {
          DeleteVariable(*FetchVar(v)->GetName());
        }
    }

    _Variable *theReceptacle = FetchVar(v); //this will return nil for TREE

    if (simpleParameters.lData[r] == 0) { // number
      q = p;
      while (!(((('0' <= data->sData[q]) && (data->sData[q] <= '9')) ||
                (data->sData[q] == '-') || (data->sData[q] == '+') ||
                (data->sData[q] == '.') || (data->sData[q] == 'e') ||
                (data->sData[q] == 'E'))) && (q < data->sLength)) {
        q++;
      }
      p = q;
      while (((('0' <= data->sData[q]) && (data->sData[q] <= '9')) ||
              (data->sData[q] == '-') || (data->sData[q] == '+') ||
              (data->sData[q] == '.') || (data->sData[q] == 'e') ||
              (data->sData[q] == 'E')) && (q < data->sLength)) {
        q++;
      }
      theReceptacle->SetValue(new _Constant(data->Cut(p, q - 1).toNum()),
                              false);
      while ((q < data->sLength - 1) && (isspace(data->sData[q + 1]))) {
        q++;
      }
    } else {
      if (simpleParameters.lData[r] == 3) { // string
        q = 0;
        bool startFound = false;
        while (q + p < data->sLength) {
          char c = data->sData[q + p];
          if (!startFound) {
            if (!isspace(c)) {
              p += q;
              startFound = true;
              q = 0;
            }
          } else if (c == '\n' || c == '\r' || c == '\t') {
            break;
          }
          q++;
        }
        if (startFound) {
          theReceptacle->SetValue(
              new _FString(new _String(*data, p, q + p - 1)), false);
        } else {
          theReceptacle->SetValue(new _FString, false);
        }

        p += q;
        r++;
        continue;
      } else if (simpleParameters.lData[r] == 5) { // raw
        theReceptacle->SetValue(new _FString(new _String(*data, p, -1)), false);
        p = data->sLength;
        r++;
        continue;
      } else {
        if (simpleParameters.lData[r] == 6) { // lines

          _String inData(*data, p, -1);

          _List lines;

          long lastP = 0, loopP = 0;

          for (loopP = 0; loopP < inData.sLength; loopP++) {
            if (inData.sData[loopP] == '\r' || inData.sData[loopP] == '\n') {
              if (lastP < loopP) {
                lines.AppendNewInstance(new _String(inData, lastP, loopP - 1));
              } else {
                lines &&&empty;
              }

              lastP = loopP + 1;

              if (lastP < inData.sLength && (inData.sData[lastP] == '\r' ||
                                             inData.sData[lastP] == '\n') &&
                  (inData.sData[lastP] != inData.sData[lastP - 1])) {
                lastP++;
              }
              loopP = lastP - 1;
            }
          }

          if (lastP < inData.sLength && lastP < loopP) {
            lines.AppendNewInstance(new _String(inData, lastP, loopP - 1));
          } else if (lines.lLength == 0) {
            lines &&&empty;
          }

          theReceptacle->SetValue(new _Matrix(lines), false);
          p = data->sLength;
          r++;
          continue;

        } else {
          char delimiter1 = (simpleParameters.lData[r] == 2) ? '(' : '{',
               delimiter2 = delimiter1 == '{' ? '}' : ')';
          q = data->Find(delimiter1, p, -1);
          if (q == -1) {
            break;
          }
          p = q;
          t = 0;
          do {
            if (data->sData[q] == delimiter1) {
              t++;
            } else if (data->sData[q] == delimiter2) {
              t--;
            }
            q++;
          } while (t && (q < data->sLength));
          if (t) {
            break;
          }
          if ((simpleParameters.lData[r] == 1) || (simpleParameters.lData[r] == 4)) {
            // matrix
            _String localData(*data, p, q - 1);
            _Matrix *newMatrixValue =
                new _Matrix(localData, simpleParameters.lData[r] == 4);
            checkPointer(newMatrixValue);
            theReceptacle->SetValue(newMatrixValue, false);
          } else {
            long varID = LocateVarByName(namespacedParameter);
            if (varID >= 0) {
              if (FetchVar(varID)->ObjectClass() == TREE) {
                DeleteVariable(*FetchVar(varID)->GetName());
              }
            }
            _String treeString(*data, p, q - 1);
            _TheTree dummyTree(namespacedParameter, treeString);
          }
        }
      }
    }
    p = q + 1;
    r++;
  }

  if (r < simpleParameters.lLength) {
    ReportWarning("fscanf could not read all the parameters requested.");
    iseof->SetValue(new _Constant(1.0), false);
  } else {
    iseof->SetValue(new _Constant(0.0), false);
  }

  if (skipDataDelete) {
    scanfLastReadPosition += p - p2;
  } else {
    scanfLastReadPosition += p;
    DeleteObject(data);
  }
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase31(_ExecutionList &chain) {

  // 20100312 SLKP: added matrix-expression based model
  // definitions
  chain.currentCommand++;

  // first check to see if matrix parameters here are valid
  bool usingLastDefMatrix = false, doExpressionBased = false;

  _Formula *isExpressionBased = nil;

  _String *parameterName, errMsg,
      arg0 = chain.AddNameSpaceToID(*(_String *)parameters(0));

  long f, f2 = -1, matrixDim, f3, multFreqs = 1;

  if (parameters.lLength > 3) {
    parameterName = (_String *)parameters.lData[3];
    if (parameterName->Equal(&ModelTrainNNFlag)) {
      _String arg1 = chain.AddNameSpaceToID(*(_String *)parameters(1));
      TrainModelNN(&arg0, &arg1);
      return;
    } else if (parameterName->Equal(&explicitFormMExp)) {
      doExpressionBased = true;
      multFreqs = 0;
    } else {
      multFreqs = ProcessNumericArgument(parameterName, chain.nameSpacePrefix);
    }
  }

  _Matrix *checkMatrix = nil;

  parameterName = (_String *)parameters.lData[1];

  if (parameterName->Equal(&useLastDefinedMatrix)) {
    if (lastMatrixDeclared < 0) {
      errMsg = "First Call to Model. USE_LAST_DEFINED_MATRIX is meaningless.";
      acknError(errMsg);
      return;
    }
    f3 = lastMatrixDeclared;
    f = modelMatrixIndices[f3];
    usingLastDefMatrix = true;
  } else {
    if (doExpressionBased) {
      _String matrixExpression(ProcessLiteralArgument(
          (_String *)parameters.lData[1], chain.nameSpacePrefix)),
          defErrMsg =
              _String("The expression for the explicit matrix exponential "
                      "passed to Model must be a valid matrix-valued HyPhy "
                      "formula that is not an assignment.") & ':' &
              matrixExpression;
      // try to parse the expression, confirm that it is a square  matrix,
      // and that it is a valid transition matrix
      isExpressionBased = (_Formula *)checkPointer(new _Formula);
      _FormulaParsingContext fpc(nil, chain.nameSpacePrefix);
      long parseCode = Parse(isExpressionBased, matrixExpression, fpc);
      if (parseCode != HY_FORMULA_EXPRESSION ||
          isExpressionBased->ObjectClass() != MATRIX) {
        WarnError(defErrMsg);
        return;
      }
      //for (unsigned long k = 0; k < isExpressionBased
      checkMatrix = (_Matrix *)isExpressionBased->Compute();

    } else {

      parameterName = (_String *)parameters.lData[1];
      _String augName(chain.AddNameSpaceToID(*parameterName));
      f = LocateVarByName(augName);

      if (f < 0) {
        WarnError(*parameterName &
                  " has not been defined prior to the call to Model = ...");
        return;
      }

      _Variable *checkVar = usingLastDefMatrix ? LocateVar(f) : FetchVar(f);
      if (checkVar->ObjectClass() != MATRIX) {
        WarnError(*parameterName &
                  " must refer to a matrix in the call to Model = ...");
        return;
      }

      checkMatrix = (_Matrix *)checkVar->GetValue();

    }
  }

  matrixDim = checkMatrix->GetHDim();
  if (matrixDim != checkMatrix->GetVDim() || matrixDim < 2) {
    WarnError(
        *parameterName &
        " must be a square matrix of dimension>=2 in the call to Model = ...");
    return;
  }

  // so far so good
  // this is the frequency matrix (if there is one!)
  parameterName = (_String *)parameters.lData[2]; 
  _String freqNameTag(chain.AddNameSpaceToID(*parameterName));

  f2 = LocateVarByName(freqNameTag);

  if (f2 < 0) {
    WarnError(*parameterName &
              " has not been defined prior to the call to Model = ...");
    return;
  }

  _Variable *checkVar = FetchVar(f2);
  if (checkVar->ObjectClass() != MATRIX) {
    WarnError(*parameterName &
              " must refer to a column/row vector in the call to Model = ...");
    return;
  }

  checkMatrix = (_Matrix *)checkVar->GetValue();
  if (checkMatrix->GetVDim() == 1) {
    if (checkMatrix->GetHDim() != matrixDim) {
      WarnError(*parameterName &
                " must be a column vector of the same dimension as the model "
                "matrix in the call to Model = ...");
      return;
    }
  } else if (checkMatrix->GetHDim() == 1) {
    if (checkMatrix->GetVDim() != matrixDim) {
      WarnError(*parameterName &
                " must be a row vector of the same dimension as the model "
                "matrix in the call to Model = ...");
      return;
    }

    errMsg = *parameterName &
             " has been transposed to the default column vector setting ";
    checkMatrix->Transpose();
    ReportWarning(errMsg);

  } else {
    WarnError(*parameterName &
              " must refer to a column/row vector in the call to Model = ...");

    return;
  }

  if (usingLastDefMatrix) {
    if (modelFrequenciesIndices[f3] < 0) {
      f2 = -f2 - 1;
    }
  } else if (multFreqs == 0) { // optional flag present
    f2 = -f2 - 1;
  }

  long existingIndex = modelNames.Find(&arg0);

  if (existingIndex == -1) { // name not found
    lastMatrixDeclared = modelNames.Find(&empty);

    if (lastMatrixDeclared >= 0) {
      modelNames.Replace(lastMatrixDeclared, &arg0, true);
      modelTypeList.lData[lastMatrixDeclared] =
          isExpressionBased ? matrixDim : 0;
      if (isExpressionBased) {
        modelMatrixIndices.lData[lastMatrixDeclared] = (long) isExpressionBased;
      } else {
        modelMatrixIndices.lData[lastMatrixDeclared] =
            (usingLastDefMatrix ? f : variableNames.GetXtra(f));
      }

      if (f2 >= 0) {
        modelFrequenciesIndices.lData[lastMatrixDeclared] =
            variableNames.GetXtra(f2);
      } else {
        modelFrequenciesIndices.lData[lastMatrixDeclared] =
            -variableNames.GetXtra(-f2 - 1) - 1;
      }
    } else {
      modelNames &&&arg0;
      modelTypeList << (isExpressionBased ? matrixDim : 0);
      if (isExpressionBased) {
        modelMatrixIndices << (long) isExpressionBased;
      } else {
        modelMatrixIndices
            << (usingLastDefMatrix ? f : variableNames.GetXtra(f));
      }
      if (f2 >= 0) {
        modelFrequenciesIndices << variableNames.GetXtra(f2);
      } else {
        modelFrequenciesIndices << -variableNames.GetXtra(-f2 - 1) - 1;
      }
      lastMatrixDeclared = modelNames.lLength - 1;
    }
  } else {
    modelNames.Replace(existingIndex, &arg0, true);
    if (modelTypeList.lData[existingIndex]) {
      delete ((_Formula *)modelMatrixIndices[existingIndex]);
    }

    modelTypeList.lData[existingIndex] = isExpressionBased ? matrixDim : 0;
    if (isExpressionBased) {
      modelMatrixIndices[existingIndex] = (long) isExpressionBased;
    } else {
      modelMatrixIndices[existingIndex] =
          usingLastDefMatrix ? f : variableNames.GetXtra(f);
    }

    if (f2 >= 0) {
      modelFrequenciesIndices[existingIndex] = variableNames.GetXtra(f2);
    } else {
      modelFrequenciesIndices[existingIndex] =
          -variableNames.GetXtra(-f2 - 1) - 1;
    }

    lastMatrixDeclared = existingIndex;
  }
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase32(_ExecutionList &chain) {

  chain.currentCommand++;
  // first check to see if matrix parameters here are valid
  long f = LocateVarByName(AppendContainerName(*(_String *)parameters(3),
                                               chain.nameSpacePrefix)),
       fixedLength = ProcessNumericArgument((_String *)parameters(2),
                                            chain.nameSpacePrefix);

  _String saveTheArg;
  _SimpleList sel, exclusions;

  _Variable *holder;

  if (fixedLength < 0) {
    fixedLength = 1;
    saveTheArg = *(_String *)parameters(2) &
                 " should represent a non-negative integer in call to "
                 "ChoiceList. The value was reset to 1";
    ReportWarning(saveTheArg);
  }

  if (f >= 0) {
    holder = FetchVar(f);
    if (holder->ObjectClass() == NUMBER) {
      if ((f = holder->Value()) >= 0) {
        exclusions << f;
      }
    } else if (holder->ObjectClass() == MATRIX) {
      _Matrix *theExcl = (_Matrix *)holder->GetValue()->Compute();
      for (long k = theExcl->GetHDim() * theExcl->GetVDim() - 1; k >= 0; k--) {
        f = (*theExcl)[k];
        if (f >= 0) {
          exclusions << f;
        }
      }
      exclusions.Sort();
    }
  }

  holder = CheckReceptacle(
      &AppendContainerName(*(_String *)parameters(0), chain.nameSpacePrefix),
      "Choice List", true);
  holder->SetBounds(-2.0, holder->GetUpperBound());

  bool validChoices = simpleParameters.lData[0] == 0;

  // some data structure present - process accordingly
  if (simpleParameters.lData[0]) {
    saveTheArg = *(_String *)parameters(4);
    // see if there is a "standard argument"
    _List choices;
    if (saveTheArg == _String("LikelihoodFunction")) {
      parameters.Delete(4);
      for (f = 0; f < likeFuncList.lLength; f++) {
        if (exclusions.BinaryFind(f) >= 0) {
          continue;
        }

        if (likeFuncList.lData[f]) {
          _List thisPair;
          thisPair << likeFuncNamesList(f);
          _String likeFuncDesc("Likelihood Function \"");
          likeFuncDesc =
              likeFuncDesc & *(_String *)likeFuncNamesList(f) & ("\".");

          thisPair &&&likeFuncDesc;
          choices &&&thisPair;
        }
      }
      validChoices = true;
      parameters &&&choices;
    } else {
      _String nmspName = AppendContainerName(saveTheArg, chain.nameSpacePrefix);
      f = FindDataSetFilterName(nmspName);
      if (f >= 0) {
        parameters.Delete(4);
        _DataSetFilter *theFilter = (_DataSetFilter *)dataSetFilterList(f);
        for (f = 0; f < theFilter->NumberSpecies(); f++) {
          if (exclusions.BinaryFind(f) >= 0) {
            continue;
          }

          _List thisPair;
          thisPair << theFilter->GetData()->GetNames()(f);
          _String spNumber("Taxon ");
          spNumber = spNumber & (f + 1) & '(' &
                     *(_String *)theFilter->GetData()->GetNames()(f) & ')';
          thisPair &&&spNumber;
          choices &&&thisPair;
        }
        validChoices = true;
        parameters &&&choices;
      } else {
        f = FindDataSetName(nmspName);
        if (f >= 0) {
          parameters.Delete(4);
          _DataSet *theSet = (_DataSet *)dataSetList(f);
          for (f = 0; f < theSet->NoOfSpecies(); f++) {
            if (exclusions.BinaryFind(f) >= 0) {
              continue;
            }
            _List thisPair;
            thisPair << theSet->GetNames()(f);
            _String spNumber("Taxon ");
            spNumber = spNumber & (f + 1) & '(' &
                       *(_String *)theSet->GetNames()(f) & ')';
            thisPair &&&spNumber;
            choices &&&thisPair;
          }
          validChoices = true;
          parameters &&&choices;
        } else {
          if (saveTheArg == lastModelParameterList) {
            f = lastMatrixDeclared;
          } else {
            f = modelNames.Find(&nmspName);
          }

          if (f >= 0) {
            parameters.Delete(4);
            _Variable *theSet = LocateVar(modelMatrixIndices.lData[f]);
            _SimpleList modelParms;
            _String ts("All Parameters");
            _List tl;
            tl &&&ts;
            ts = "All local model parameters are constrained";
            tl &&&ts;
            choices &&&tl;
            _AVLList modelParmsA(&modelParms);
            theSet->ScanForVariables(modelParmsA, false);
            modelParmsA.ReorderList();
            for (f = 0; f < modelParms.lLength; f++) {
              if (exclusions.BinaryFind(f) >= 0) {
                continue;
              }

              _List thisPair;
              thisPair << LocateVar(modelParms.lData[f])->GetName();
              _String spNumber("Constrain parameter ");
              spNumber = spNumber & *LocateVar(modelParms.lData[f])->GetName();
              thisPair &&&spNumber;
              choices &&&thisPair;
            }
            validChoices = true;
            parameters &&&choices;
          } else {
            f = LocateVarByName(nmspName);
            if (f >= 0) {
              _Variable *theV = FetchVar(f);
              if (theV->ObjectClass() == MATRIX) {
                _Matrix *vM = (_Matrix *)theV->GetValue();
                if (vM->IsAStringMatrix() && (vM->GetVDim() == 2)) {
                  parameters.Delete(4);
                  for (f = 0; f < vM->GetHDim(); f++) {
                    if (exclusions.BinaryFind(f) < 0) {
                      _Formula *f1 = vM->GetFormula(f, 0),
                               *f2 = vM->GetFormula(f, 1);

                      if (f1 && f2) {
                        _PMathObj p1 = f1->Compute(), p2 = f2->Compute();

                        if (p1 && p2 && (p1->ObjectClass() == STRING) &&
                            (p2->ObjectClass() == STRING)) {
                          _List thisPair;
                          thisPair << ((_FString *)p1)->theString;
                          thisPair << ((_FString *)p2)->theString;
                          choices &&&thisPair;
                        }
                      }
                    }
                  }
                  validChoices = true;
                  parameters &&&choices;
                }
              }
            }
          }
        }
      }
    }
  }

  if (validChoices) {
    long choice = -1;
    _List *theChoices = (_List *)parameters(4);
    if (fixedLength > theChoices->lLength) {
      _String e = "List of selections is too short in ChoiceList";
      acknError(e);
    } else {
      if (chain.stdinRedirect) {
        if (fixedLength == 1) {
          _String buffer(chain.FetchFromStdinRedirect());
          for (choice = 0; choice < theChoices->lLength; choice++)
            if (buffer.Equal((_String *)(*(_List *)(*theChoices)(choice))(0))) {
              break;
            }
          if (choice == theChoices->lLength) {
            choice = -1;
            WarnError(_String("Not a valid option: '") & buffer &
                      "' passed to Choice List '" &
                      ((_String *)parameters(1))->sData &
                      "' using redirected stdin input");
            return;
          }
        } else {
          if (fixedLength > 0) {
            while (sel.lLength < fixedLength) {
              _String buffer(chain.FetchFromStdinRedirect());
              for (choice = 0; choice < theChoices->lLength; choice++)
                if (buffer.Equal(
                        (_String *)(*(_List *)(*theChoices)(choice))(0))) {
                  break;
                }
              if (choice < theChoices->lLength && sel.Find(choice) == -1) {
                sel << choice - 1;
              } else {
                break;
              }
            }
            if (sel.lLength < fixedLength) {
              WarnError("Failed to make the required number of choices in "
                        "ChoiceList using redirected stdin input.");
              return;
            }
          } else
            while (1) {
              _String buffer(chain.FetchFromStdinRedirect());
              if (buffer.sLength) {
                for (choice = 0; choice < theChoices->lLength; choice++)
                  if (buffer.Equal(
                          (_String *)(*(_List *)(*theChoices)(choice))(0))) {
                    break;
                  }

                if (choice < theChoices->lLength && sel.Find(choice) == -1) {
                  sel << choice;
                } else {
                  WarnError(
                      _String("Not a valid (or duplicate) option: '") & buffer &
                      "' passed to ChoiceList (with multiple selections) '" &
                      ((_String *)parameters(1))->sData &
                      "' using redirected stdin input");
                  return;
                }
              } else {
                break;
              }
            }
        }
      } else {
#ifdef __HEADLESS__
        WarnError("Unhandled request for data from standard input in "
                  "ChoiceList in headless HyPhy");
        return;
#else
#if defined __HYPHYQT__
        SetStatusLine("Waiting for user selection.");
        _String *param = (_String *)parameters(1);

        _SimpleList std(2, 0, 1), all(theChoices->lLength, 0, 1);

        choice =
            HandleListSelection(*theChoices, std, all, *param, sel, fixedLength,
                                (Ptr) _hyPrimaryConsoleWindow);
#else
        _String *param = (_String *)parameters(1);
        printf("\n\n\t\t\t+");

        for (f = 1; f < param->sLength + 1; f++) {
          printf("-");
        }

        printf("+\n\t\t\t|%s|\n\t\t\t+", param->getStr());

        for (f = 1; f < param->sLength + 1; f++) {
          printf("-");
        }

        printf("+\n\n");

        long loopits = 1;

        if (fixedLength == 1) {
          while (choice == -1) {
            for (choice = 0; choice < theChoices->lLength; choice++) {
              printf(
                  "\n\t(%ld):[%s] %s", choice + 1,
                  ((_String *)(*(_List *)(*theChoices)(choice))(0))->getStr(),
                  ((_String *)(*(_List *)(*theChoices)(choice))(1))->getStr());
            }

            printf("\n\n Please choose an option (or press q to cancel "
                   "selection):");
            _String buffer(StringFromConsole());
            if (buffer.sData[0] == 'q' || buffer.sData[0] == 'Q') {
              choice = -1;
              break;
            }
            choice = buffer.toNum();
            if (choice < 1 || choice > theChoices->lLength) {
              choice = -1;
              if (loopits++ > 10) {
                FlagError("Failed to make a valid selection in ChoiceList "
                          "after 10 tries");
                return;
              }
            } else {
              choice--;
            }
          }
        } else {
          if (fixedLength > 0)
            while (sel.lLength < fixedLength) {
              for (choice = 0; choice < theChoices->lLength; choice++) {
                if (sel.Find(choice) == -1) {
                  printf("\n\t(%ld):%s", choice + 1,
                         ((_String *)(*(_List *)(*theChoices)(choice))(1))
                             ->getStr());
                }
              }
              printf("\n\n Please choose option %ld of %ld (or press q to "
                     "cancel selection):",
                     sel.lLength + 1, fixedLength);
              _String buffer(StringFromConsole());
              if (buffer.sData[0] == 'q' || buffer.sData[0] == 'Q') {
                choice = -1;
                break;
              }
              choice = buffer.toNum();
              if ((choice >= 1) && (choice <= theChoices->lLength)) {
                if (sel.Find(choice) == -1) {
                  sel << choice - 1;
                }
              } else {
                if (loopits++ > 10) {
                  FlagError("Failed to make a valid selection in ChoiceList "
                            "after 10 tries");
                  return;
                }
              }
            }
          else
            while (1) {
              for (choice = 0; choice < theChoices->lLength; choice++) {
                if (sel.Find(choice) == -1) {
                  printf("\n\t(%ld):[%s] %s", choice + 1,
                         ((_String *)(*(_List *)(*theChoices)(choice))(0))
                             ->getStr(),
                         ((_String *)(*(_List *)(*theChoices)(choice))(1))
                             ->getStr());
                }
              }
              printf("\n\n Please choose option %ld, enter d to complete "
                     "selection, enter q to cancel:",
                     sel.lLength + 1);
              _String buffer(StringFromConsole());
              if (buffer.sData[0] == 'q' || buffer.sData[0] == 'Q') {
                choice = -1;
                break;
              }
              if (buffer.sData[0] == 'd' || buffer.sData[0] == 'D') {
                break;
              }

              choice = buffer.toNum();
              if ((choice >= 1) && (choice <= theChoices->lLength)) {
                if (sel.Find(choice) == -1) {
                  sel << choice - 1;
                }
              } else {
                if (loopits++ > 10) {
                  FlagError("Failed to make a valid selection in ChoiceList "
                            "after 10 tries");
                  return;
                }
              }
            }
        }
#endif
#endif
      }

      _Variable *sStrV = CheckReceptacle(&selectionStrings, empty, false);

      if (fixedLength == 1) {
        if (choice >= 0) {
          _FString choiceString(*(_String *)((_List *)(*theChoices)(choice))
                                    ->lData[0]);
          sStrV->SetValue(&choiceString);
          for (long k = 0; k < exclusions.lLength; k++) {
            if (choice >= exclusions.lData[k]) {
              choice++;
            } else {
              break;
            }
          }
        }
        holder->SetValue(new _Constant(choice), false);
      } else {
        if (fixedLength == 0) {
          fixedLength = sel.lLength;
          if (fixedLength == 0) {
            fixedLength = 1;
          }
        }
        sel.Sort();
        _Matrix selVector(1, fixedLength, false, true);
        _Matrix selMatrix(1, fixedLength, false, true);
        if (choice == -1) {
          selVector[0] = -1;
        } else {
          for (f = 0; f < fixedLength; f++) {
            choice = sel.lData[f];

            _FString *choiceString = new _FString(
                (*(_String *)((_List *)(*theChoices)(choice))->lData[0]));
            _Formula sf(choiceString);
            selMatrix.MStore(0, f, sf);
            for (long k = 0; k < exclusions.lLength; k++) {
              if (choice >= exclusions.lData[k]) {
                choice++;
              } else {
                break;
              }
            }
            selVector[f] = choice;
            //DeleteObject (choiceString);
          }
          sStrV->SetValue(&selMatrix);
        }
        holder->SetValue(&selVector);
      }

      if (choice < 0) {
        terminateExecution = true;
      }
    }
  } else {
    WarnError("List of selections is invalid in ChoiceList");
  }

  if (simpleParameters.lData[0]) {
    parameters.Delete(4);
    parameters &&&saveTheArg;
  }
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase36(_ExecutionList &chain) {
  chain.currentCommand++;
  _String *currentArgument = (_String *)parameters(0), errMsg, result;

  long f = dataSetNamesList.Find(
           &AppendContainerName(*currentArgument, chain.nameSpacePrefix)),
       s, k, m;

  if (f < 0) {
    ReportWarning(*currentArgument &
                  " is not a valid data set in call to OpenDataPanel");
    return;
  }
  _DataSet *theDS = (_DataSet *)dataSetList(f);

  // process species list
  result =
      ProcessLiteralArgument((_String *)parameters(1), chain.nameSpacePrefix);
  if (result.sLength) {
    result.Insert('"', 0);
    result.Insert('"', -1);
  }
  _SimpleList speciesList;
  theDS->ProcessPartition(result, speciesList, true);

  // check the validity of list entries
  s = theDS->NoOfSpecies();

  for (m = speciesList.lLength - 1; m >= 0; m--) {
    k = speciesList.lData[m];
    if (k < 0 || k >= s) {
      speciesList.Delete(m);
      m--;
    }
  }

  if (speciesList.lLength == s) {
    speciesList.Clear();
  }

#if !defined __UNIX__ && !defined __HEADLESS__

  _HYDataPanel *newDP = new _HYDataPanel(empty, empty);
  if (speciesList.lLength) {
    newDP->SetDataSetReference(*(_String *)dataSetNamesList(f), &speciesList);
  } else {
    newDP->SetDataSetReference(*(_String *)dataSetNamesList(f), nil);
  }
  result = dataPanelSourcePath;
  result = ProcessLiteralArgument(&result, chain.nameSpacePrefix);
  if (result.sLength) {
    newDP->SetFilePath(result);
  }
  currentArgument = (_String *)parameters(3);
  newDP->RestorePartInfo(currentArgument);
  currentArgument = (_String *)parameters(2);
  newDP->RestorePanelSettings(currentArgument);
  currentArgument = (_String *)parameters(1);
  newDP->SetSavePath(chain.sourceFile);
  newDP->BringToFront();
  if (parameters.lLength > 4) {
    newDP->BuildLikelihoodFunction((_String *)parameters(4));
    newDP->RestoreSavedLFs();
  }
#endif
}

//______________________________________________________________________________
// GetInformation()
void _ElementaryCommand::ExecuteCase37(_ExecutionList &chain) {
  chain.currentCommand++;

  _String matrixName = chain.AddNameSpaceToID(*(_String *)parameters(0)),
          *objectName = (_String *)parameters(1);

  long sID;
  if (parameters.lLength > 2) {
    sID =
        ProcessNumericArgument((_String *)parameters(2), chain.nameSpacePrefix);
  }

  _Matrix *result = nil;

  // object is a non-empty string
  if (objectName->sLength > 2 && objectName->sData[0] == '"' && objectName->sData[objectName->sLength - 1] == '"') {
   // regular expression
    _String regExp = GetStringFromFormula(objectName, chain.nameSpacePrefix);
    int errNo = 0;
    Ptr regex = PrepRegExp(&regExp, errNo, true);
    if (regex) {

      _List matches;
      _SimpleList tcache;
      long iv, k = variableNames.Traverser(tcache, iv, variableNames.GetRoot());

      for (; k >= 0; k = variableNames.Traverser(tcache, iv)) {
        _String *vName = (_String *)variableNames.Retrieve(k);
        _SimpleList mtch;
        vName->RegExpMatch(regex, mtch);
        if (mtch.lLength) {
          matches << vName;
        }

      }

      if (matches.lLength) {
        result = new _Matrix(matches);
      }

      FlushRegExp(regex);
    } else {
      WarnError(GetRegExpError(errNo));
    }
  } else { // object is not a string, is some kind of variable
    _String objectNameID = chain.AddNameSpaceToID(*objectName);
    long f = LocateVarByName(objectNameID);
    if (f >= 0) { // it's a numeric variable
      _Variable *theObject = FetchVar(f);
      if (theObject->ObjectClass() == STRING) {
        objectNameID = _String((_String *)theObject->Compute()->toStr());
        theObject = FetchVar(LocateVarByName(objectNameID));
      }
      if (theObject) {
        if (theObject->IsCategory()) {
          _CategoryVariable *thisCV = (_CategoryVariable *)theObject;
          thisCV->Refresh();

          _Matrix *values = thisCV->GetValues(),
                  *weights = thisCV->GetWeights(!thisCV->IsUncorrelated());

          f = values->GetHDim() * values->GetVDim();
          checkPointer(result = new _Matrix(2, f, false, true));

          for (long k = 0; k < f; k++) {
            result->theData[k] = values->theData[k];
            result->theData[f + k] = weights->theData[k];
          }
        } else {
          if (theObject->ObjectClass() == TREE_NODE) {
            _CalcNode *theNode = (_CalcNode *)theObject;
            if (theNode->GetModelIndex() != HY_NO_MODEL) {
              checkPointer(result = new _Matrix);
              theNode->RecomputeMatrix(0, 1, result);
            }
          }

          if ((!result) && theObject->ObjectClass() == NUMBER) {
            checkPointer(result = new _Matrix(1, 3, false, true));
            result->theData[0] = theObject->Compute()->Value();
            result->theData[1] = theObject->GetLowerBound();
            result->theData[2] = theObject->GetUpperBound();
          }
        }
      }
    } else {
      f = likeFuncNamesList.Find(&objectNameID);
      if (f >= 0) { // it's a likelihood function
        _LikelihoodFunction *lf = (_LikelihoodFunction *)likeFuncList(f);
        f = lf->GetCategoryVars().lLength;
        if (f == 0) {
          f++;
        }

        _List catVars;

        for (long k = 0; k < lf->GetCategoryVars().lLength; k++) {
          _String varName =
              *LocateVar(lf->GetCategoryVars().lData[k])->GetName();
          catVars &&&varName;
        }
        result = (_Matrix *)checkPointer(new _Matrix(catVars));
      } else {
        if ((f = dataSetFilterNamesList.Find(&objectNameID)) >= 0) {
          // return a vector of strings - each with actual characters of the
          // corresponding sequence
          _DataSetFilter *daFilter = (_DataSetFilter *)dataSetFilterList(f);
          result = daFilter->GetFilterCharacters();
        } else {
          // it's a tree node with a rate matrix assigned
          f = FindModelName(objectNameID);

          // for models, return the list of variables in the model
          if (f >= 0) {
            _SimpleList modelParms;
            _AVLList modelParmsA(&modelParms);
            LocateVar(modelMatrixIndices.lData[f])->ScanForVariables(modelParmsA, false);
            _List modelPNames;

            for (unsigned long vi = 0; vi < modelParms.lLength; vi++) {
              modelPNames << LocateVar(modelParms.lData[vi])->GetName();
            }
            result = new _Matrix(modelPNames);
          }
        }
      }
    }
  }

  if (!result) {
    result = new _Matrix(0, 0, false, false);
  }

  CheckReceptacleAndStore(&matrixName, empty, true, result, true);
  DeleteObject(result);

}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase43(_ExecutionList &chain) {
  chain.currentCommand++;

  _String *currentArgument = (_String *)parameters(0), result;

  _Variable *theReceptacle = CheckReceptacle(
      &AppendContainerName(*currentArgument, chain.nameSpacePrefix),
      code == 43 ? blFindRoot : blIntegrate, true);

  if (theReceptacle) {
    _String exprString = *(_String *)parameters(1);
    _Formula theExpression(exprString);

    currentArgument = (_String *)parameters(2);
    long f = LocateVarByName(
        AppendContainerName(*currentArgument, chain.nameSpacePrefix));

    if (f < 0) {
      ReportWarning(*currentArgument & " is not an existing variable to solve "
                                       "for in call to FindRoot/Integrate.");
      return;
    }

    if (terminateExecution) {
      return;
    }

    _Formula *dF = (code == 43) ? theExpression.Differentiate(
                                      *(_String *)parameters(2), false)
                                : nil;

    _Parameter lb = ProcessNumericArgument((_String *)parameters(3),
                                           chain.nameSpacePrefix),
               ub = ProcessNumericArgument((_String *)parameters(4),
                                           chain.nameSpacePrefix);

    if (ub <= lb && code == 48) {
      ReportWarning(
          _String('[') & lb & ',' & ub &
          "] is not a valid search interval in call to FindRoot/Integrate");
      return;
    }

    if (code == 43) {
      if (dF) {
        theReceptacle->SetValue(
            new _Constant(theExpression.Newton(*dF, FetchVar(f), 0.0, lb, ub)),
            false);
      } else {
        theReceptacle->SetValue(
            new _Constant(theExpression.Brent(FetchVar(f), lb, ub)), false);
      }
    } else {
      theReceptacle->SetValue(new _Constant(theExpression.Integral(
                                  FetchVar(f), lb, ub, ub - lb > 1e10)),
                              false);
    }

    if (dF) {
      delete (dF);
    }
  }
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase44(_ExecutionList &chain) {

  chain.currentCommand++;

#ifdef __HYPHYMPI__
  _String *arg1 = (_String *)parameters(0), *arg2 = (_String *)parameters(1),
          *arg3 = parameters.lLength > 2 ? (_String *)parameters(2) : nil,
          *theMessage = nil;

  _Parameter nodeCount;
  checkParameter(mpiNodeCount, nodeCount, 1);

  long destID = ProcessNumericArgument(arg1, chain.nameSpacePrefix), g;

  if (!numericalParameterSuccessFlag || destID < 0 || destID >= nodeCount) {
    WarnError(*arg1 & " is not a valid MPI node ID in call to MPISend.");
    return;
  }

  if (arg3) {
    _AssociativeList *ar = (_AssociativeList *)FetchObjectFromVariableByType(
        &AppendContainerName(*arg3, chain.nameSpacePrefix), ASSOCIATIVE_LIST);
    if (!ar) {
      WarnError(*arg3 & " is not a valid associative array for input options "
                        "in call to MPISend.");
      return;
    }
    theMessage = new _String(256L, true);
    checkPointer(theMessage);
    _String arrayID("_HYPHY_MPI_INPUT_ARRAY_");
    (*theMessage) << arrayID;
    (*theMessage) << '=';
    arg3 = ar->Serialize(arrayID);
    (*theMessage) << arg3;
    DeleteObject(arg3);
    (*theMessage) << ';';
    arrayID = *arg2;
    arrayID.ProcessFileName(false, true, (Ptr) chain.nameSpacePrefix);
    (*theMessage) << "\nExecuteAFile (\"";
    (*theMessage) << arrayID;
    (*theMessage) << "\",_HYPHY_MPI_INPUT_ARRAY_);";
    theMessage->Finalize();
  } else if ((g = FindLikeFuncName(
                 AppendContainerName(*arg2, chain.nameSpacePrefix))) >= 0) {
    checkPointer(theMessage = new _String(1024L, true));
    ((_LikelihoodFunction *)likeFuncList(g))
        ->SerializeLF(*theMessage, _hyphyLFSerializeModeOptimize);
    theMessage->Finalize();
  } else {
    theMessage =
        new _String(ProcessLiteralArgument(arg2, chain.nameSpacePrefix));
  }

  if (theMessage == nil || theMessage->sLength == 0) {
    WarnError(
        *arg2 &
        " is not a valid (or is an empty) string (LF ID) in call to MPISend.");
  } else {
    MPISendString(*theMessage, destID);
  }

  DeleteObject(theMessage);

#else
  WarnError("MPISend can't be used by non-MPI versions of HyPhy.");
#endif

}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase45(_ExecutionList &chain) {
  chain.currentCommand++;

#ifdef __HYPHYMPI__
  _String *arg1 = (_String *)parameters(0), *arg2 = (_String *)parameters(1),
          *arg3 = (_String *)parameters(2);

  _Parameter nodeCount;
  checkParameter(mpiNodeCount, nodeCount, 1);

  long srcT = ProcessNumericArgument(arg1, chain.nameSpacePrefix), srcID, g;

  if ((!numericalParameterSuccessFlag) || (srcT < -1) || (srcT >= nodeCount)) {
    WarnError(*arg1 & " is not a valid MPI node ID in call to MPIReceive.");
    return;
  }

  _Variable *idVar = CheckReceptacle(
                &AppendContainerName(*arg2, chain.nameSpacePrefix),
                "MPIReceive"),
            *mVar = CheckReceptacle(
                &AppendContainerName(*arg3, chain.nameSpacePrefix),
                "MPIReceive");

  if (!(idVar && mVar)) {
    return;
  }

  _FString *theMV = new _FString(MPIRecvString(srcT, srcID));
  checkPointer(theMV);
  idVar->SetValue(new _Constant(srcID), false);
  mVar->SetValue(theMV, false);
#else
  WarnError("MPIReceive can't be used by non-MPI versions of HyPhy.");
#endif

}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase46(_ExecutionList &chain) {

  chain.currentCommand++;

  _String *arg1 = (_String *)parameters(1), *arg2 = (_String *)parameters(0),
          errMsg;

  long k = dataSetFilterNamesList.Find(
      &AppendContainerName(*arg1, chain.nameSpacePrefix));

  if (k < 0) {
    errMsg = *arg1 & " is not a defined data set filter ID ";
  } else {
    _DataSetFilter *dsf = (_DataSetFilter *)dataSetFilterList(k);
    _Variable *stVar = CheckReceptacle(
        &AppendContainerName(*arg2, chain.nameSpacePrefix), "GetDataInfo");

    if (stVar) {
      if (parameters.lLength == 2) {
        _Matrix *res = new _Matrix(1, dsf->duplicateMap.lLength, false, true);
        checkPointer(res);
        for (k = 0; k < dsf->duplicateMap.lLength; k++) {
          res->theData[k] = dsf->duplicateMap.lData[k];
        }
        stVar->SetValue(res, false);
      } else {
        if (parameters.lLength == 3) {
          _String checker = ProcessLiteralArgument((_String *)parameters(2),
                                                   chain.nameSpacePrefix);
          if (checker == _String("CHARACTERS")) {
            _List characters;
            k = dsf->GetDimension(true);
            long fd = dsf->GetUnitLength();
            for (long idx = 0; idx < k; idx++) {
              characters.AppendNewInstance(new _String(
                  dsf->ConvertCodeToLetters(dsf->CorrectCode(idx), fd)));
            }

            stVar->SetValue(new _Matrix(characters), false);
          } else if (checker == _String("PARAMETERS")) {
            _AssociativeList *parameterInfo = new _AssociativeList;
            parameterInfo->MStore("ATOM_SIZE",
                                  new _Constant(dsf->GetUnitLength()), false);
            parameterInfo->MStore("EXCLUSIONS",
                                  new _FString(dsf->GetExclusions()), false);
            parameterInfo->MStore(
                "SITES_STRING",
                new _FString(
                    (_String *)dsf->theOriginalOrder.ListToPartitionString()),
                false);
            parameterInfo->MStore(
                "SEQUENCES_STRING",
                new _FString(
                    (_String *)dsf->theNodeMap.ListToPartitionString()),
                false);
            stVar->SetValue(parameterInfo, false);

          } else if (checker == _String("CONSENSUS")) {
            stVar->SetValue(
                new _FString(new _String(dsf->GenerateConsensusString())),
                false);
          } else {
            long seqID = ProcessNumericArgument((_String *)parameters(2),
                                                chain.nameSpacePrefix);
            if (seqID >= 0 && seqID < dsf->NumberSpecies()) {
              stVar->SetValue(new _FString(dsf->GetSequenceCharacters(seqID)),
                              false);
            } else {
              // 20110916 SLKP : the option for filtering duplicate sequences
              if (seqID >= -4 && seqID <= -1) {
                _SimpleList indices, map, counts;
                long uniqueSequences =
                    dsf->FindUniqueSequences(indices, map, counts, -seqID - 1);
                _AssociativeList *parameterInfo = new _AssociativeList;
                parameterInfo->MStore("UNIQUE_SEQUENCES",
                                      new _Constant(uniqueSequences), false);
                parameterInfo->MStore("UNIQUE_INDICES", new _Matrix(indices),
                                      false);
                parameterInfo->MStore("SEQUENCE_MAP", new _Matrix(map), false);
                parameterInfo->MStore("UNIQUE_COUNTS", new _Matrix(counts),
                                      false);
                stVar->SetValue(parameterInfo, false);
              }
            }
          }
        } else {
          long seq = ProcessNumericArgument((_String *)parameters(2),
                                            chain.nameSpacePrefix),
               site = ProcessNumericArgument((_String *)parameters(3),
                                             chain.nameSpacePrefix);

          if (parameters.lLength == 4) {
            if ((seq >= 0) && (site >= 0) && (seq < dsf->NumberSpecies()) &&
                (site < dsf->NumberDistinctSites())) {
              _Matrix *res = (_Matrix *)checkPointer(
                  new _Matrix(dsf->GetDimension(true), 1, false, true));

              _Parameter onlyTheIndex = 0.0;
              checkParameter(getDataInfoReturnsOnlyTheIndex, onlyTheIndex, 0.0);

              long theValue = dsf->Translate2Frequencies((*dsf)(site, seq),
                                                         res->theData, true);

              if (onlyTheIndex > 0.5) {
                stVar->SetValue(new _Constant(theValue), false);
                DeleteObject(res);
              } else {
                stVar->SetValue(res, false);
              }
            } else {
              errMsg = _String(seq) & "," & _String(site) &
                       " is an invalid site index ";
            }
          } else {
            if ((seq >= 0) && (site >= 0) && (seq < dsf->NumberSpecies()) &&
                (site < dsf->NumberSpecies())) {
              _String *resFlag = (_String *)parameters(4);
              _Matrix *res;

              if (pcAmbiguitiesAverage.Equal(resFlag)) {
                res = dsf->ComputePairwiseDifferences(seq, site, 1);
              } else if (pcAmbiguitiesResolve.Equal(resFlag)) {
                res = dsf->ComputePairwiseDifferences(seq, site, 2);
              } else if (pcAmbiguitiesSkip.Equal(resFlag)) {
                res = dsf->ComputePairwiseDifferences(seq, site, 3);
              } else {
                res = dsf->ComputePairwiseDifferences(seq, site, 0);
              }

              stVar->SetValue(res, false);
            } else {
              errMsg = _String(seq) & "," & _String(site) &
                       " is an invalid sequence pair specification.";
            }

          }
        }
      }
    }
  }

  if (errMsg.sLength) {
    errMsg = errMsg & " in call to GetDataInfo ";
    WarnError(errMsg);
  }

}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase47(_ExecutionList &chain) {

  chain.currentCommand++;

  _String *arg1 = (_String *)parameters(0), *arg2 = (_String *)parameters(1),
          errMsg;

  long k = FindLikeFuncName(AppendContainerName(*arg1, chain.nameSpacePrefix));

  if (k < 0) {
    _String litArg = ProcessLiteralArgument(arg1, chain.nameSpacePrefix);
    k = FindLikeFuncName(litArg);
    if (k < 0) {
      errMsg = *arg1 & " is not a defined likelihood function ID ";
    }
  }

  if (errMsg.sLength == 0) {
    _LikelihoodFunction *lf = (_LikelihoodFunction *)likeFuncList(k);
    _String callBack = ProcessLiteralArgument(arg2, chain.nameSpacePrefix);

    k = batchLanguageFunctionNames.Find(&callBack);

    if (k < 0) {
      errMsg = *arg2 & " is not a defined user batch language function ";
    } else {
      if (batchLanguageFunctionParameters.lData[k] != 2) {
        errMsg = *arg2 & " callback function must depend on 2 parameters ";
      } else {
        lf->StateCounter(k);
      }
    }
  }

  if (errMsg.sLength) {
    errMsg = errMsg & " in call to StateCounter.";
    WarnError(errMsg);
  }
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase52(_ExecutionList &chain) {

  chain.currentCommand++;
  _String errMsg;

  // check validity of an alphabet
  long siteCount =
      ProcessNumericArgument((_String *)parameters(4), chain.nameSpacePrefix);
  _String givenState;

  if (siteCount < 1) {
    givenState =
        ProcessLiteralArgument((_String *)parameters(4), chain.nameSpacePrefix);
    siteCount = givenState.sLength;
  }

  if (siteCount < 1) {
    errMsg = *(_String *)parameters(4) &
             " must either evaluate to a positive integer or be a non-empty "
             "string of root states";
    WarnError(errMsg);
    return;
  }

  _Variable *alphabet = FetchVar(LocateVarByName(AppendContainerName(
                *(_String *)parameters(3), chain.nameSpacePrefix))),
            *treeVar = FetchVar(LocateVarByName(AppendContainerName(
                *(_String *)parameters(1), chain.nameSpacePrefix))),
            *freqVar = FetchVar(LocateVarByName(AppendContainerName(
                *(_String *)parameters(2), chain.nameSpacePrefix)));

  if (alphabet && treeVar && freqVar) {
    if (alphabet->ObjectClass() == MATRIX) {
      _Matrix *alphabetMatrix = (_Matrix *)alphabet->GetValue();

      if (alphabetMatrix->IsAStringMatrix() && alphabetMatrix->GetHDim() == 2 &&
          alphabetMatrix->GetVDim() > 1) {
        _String baseSet;

        for (long k = 0; k < alphabetMatrix->GetVDim(); k++) {
          _FString *aState =
              (_FString *)alphabetMatrix->GetFormula(0, k)->Compute();
          if (aState) {
            if (aState->theString->sLength == 1) {
              char c = aState->theString->sData[0];
              if (baseSet.Find(c) == -1) {
                baseSet = baseSet & c;
              } else {
                break;
              }
            } else {
              break;
            }
          } else {
            break;
          }
        }

        if (baseSet.sLength == alphabetMatrix->GetVDim()) {
          long unitSize = ((_FString *)alphabetMatrix->GetFormula(1, 0)
                               ->Compute())->theString->toNum();

          if (unitSize >= 1) {
            _Formula *exclusionFormula = alphabetMatrix->GetFormula(1, 1);
            _String *theExclusions = &empty;

            if (exclusionFormula)
              theExclusions =
                  ((_FString *)exclusionFormula->Compute())->theString;

            if (treeVar->ObjectClass() == TREE) {
              if (freqVar->ObjectClass() == MATRIX) {
                _TheTree *spawningTree = (_TheTree *)treeVar;

                if (!(parameters.lLength > 6 &&
                      (spawningTree->CountTreeCategories() > 1))) {

                  if (givenState.sLength > 1)
                      // root state
                      {
                    if ((givenState.sLength >= unitSize) &&
                        (givenState.sLength % unitSize == 0)) {
                      siteCount = givenState.sLength / unitSize;
                    } else {
                      errMsg = "Root state string is either too short or has "
                               "length which is not divisible by the unit size";
                    }
                  }
                  if (errMsg.sLength == 0) {
                    _TranslationTable newTT(baseSet);
                    _DataSet *ds = (_DataSet *)checkPointer(new _DataSet);

                    if (!newTT.CheckType(
                            HY_TRANSLATION_TABLE_STANDARD_NUCLEOTIDE)) {
                      ds->SetTranslationTable(
                          &newTT); // mod 20060113 to properly deal with
                                   // non-standard alphabets
                    }
                    // make a dummy
                    spawningTree->AddNodeNamesToDS(ds, true, false, 1);

                    char c = baseSet.sData[0];
                    long s = ds->GetNames().lLength;

                    if (s < 2) {
                      _String rt("Root");
                      ds->GetNames().InsertElement(&rt, 0, true);
                      s++;
                    }

                    unsigned long ssi = _String::storageIncrement;

                    if (s > ssi) {
                      _String::storageIncrement = s;
                    }

                    ds->AddSite(c);
                    for (long u = 1; u < s; u++) {
                      ds->Write2Site(0, c);
                    }
                    ds->Finalize();

                    _String::storageIncrement = ssi;
                    ds->SetNoSpecies(s);

                    _SimpleList *theMap = &ds->GetTheMap();
                    theMap->RequestSpace(siteCount * unitSize);
                    for (long filler = 0; filler < siteCount * unitSize;
                         filler++) {
                      theMap->lData[filler] = 0;
                    }

                    theMap->lLength = siteCount * unitSize;

                    _DataSetFilter *newFilter = new _DataSetFilter();
                    checkPointer(newFilter);
                    _SimpleList h, v;

                    newFilter->SetFilter(ds, unitSize, h, v, false);
                    newFilter->SetExclusions(theExclusions, true);
                    newFilter->SetupConversion();

                    /*char buffer[255];
                    snprintf (buffer, sizeof(buffer),"%d %d\n",siteCount,
                    newFilter->GetFullLengthSpecies(),unitSize);
                    BufferToConsole (buffer);
                    */
                    _Matrix *rootStates = nil;
                    if (givenState.sLength >= unitSize) {
                      rootStates = new _Matrix(1, siteCount, false, true);
                      checkPointer(rootStates);
                      _Parameter *holder =
                          new _Parameter[newFilter->GetDimension(false)];
                      checkPointer(holder);

                      for (long cc = 0; cc < siteCount; cc++) {
                        _String aState(givenState.Cut(cc * unitSize,
                                                      (cc + 1) * unitSize - 1));
                        long stateV = newFilter->Translate2Frequencies(
                            aState, holder, false);
                        if (stateV < 0) {
                          errMsg =
                              aState &
                              " found in the root state string at position " &
                              cc * unitSize & " is an invalid state";
                          break;
                        } else {
                          rootStates->theData[cc] = stateV;
                        }
                      }
                      delete[] holder;
                    }

                    if (errMsg.sLength == 0) {

                      long filterID =
                          AddFilterToList(simulationFilter, newFilter);

                      spawningTree->SetUp();
                      spawningTree->InitializeTreeFrequencies(
                          (_Matrix *)freqVar->Compute(), true);
                      errMsg =
                          *(_String *)dataSetFilterNamesList(filterID) & ',' &
                          *spawningTree->GetName() & ',' & *freqVar->GetName();

                      _LikelihoodFunction lf(errMsg, nil);

                      if (terminateExecution) {
                        return;
                      }

                      bool doInternals = false;

                      if (parameters.lLength > 5) {
                        doInternals = (ProcessNumericArgument(
                            (_String *)parameters(5), chain.nameSpacePrefix) >
                                       0.5);
                      }

                      _String spoolFile;

                      FILE *mainFile = nil;

                      errMsg = empty;

                      if (parameters.lLength > 6) {
                        spoolFile = ProcessLiteralArgument(
                            (_String *)parameters(6), chain.nameSpacePrefix);
                        spoolFile.ProcessFileName();
                        mainFile = doFileOpen(spoolFile.sData, "w");
                        if (!mainFile) {
                          errMsg = _String("Failed to open ") & spoolFile &
                                   " for writing";
                        }

                        if (doInternals) {
                          spoolFile = spoolFile & ".anc";
                        }

                      }

                      if (errMsg.sLength == 0) {
                        _DataSet *simDataSet;

                        if (mainFile) {
                          simDataSet = new _DataSet(mainFile);
                        } else {
                          simDataSet = new _DataSet(siteCount);
                        }

                        checkPointer(simDataSet);

                        _List exclusions;

                        _String *simName = new _String(AppendContainerName(
                            *(_String *)parameters(0), chain.nameSpacePrefix));
                        _String mxName = *simName & ".rates";
                        setParameter(mxName, 0.0);
                        _Variable *catValVar =
                            FetchVar(LocateVarByName(mxName));
                        _Matrix *catValues = new _Matrix(1, 1, false, true);
                        checkPointer(catValues);

                        mxName = *simName & ".rateVars";
                        setParameter(mxName, 0.0);
                        _Variable *catNameVar =
                            FetchVar(LocateVarByName(mxName));
                        _Matrix *catNames = new _Matrix(1, 1, false, true);

                        SetStatusLine("Simulating Data");
                        lf.Simulate(
                            *simDataSet, exclusions, catValues, catNames,
                            rootStates,
                            doInternals ? (mainFile ? &spoolFile : &empty)
                                        : nil);

                        SetStatusLine("Idle");

                        catValVar->SetValue(catValues, false);
                        catNameVar->SetValue(catNames, false);

                        StoreADataSet(simDataSet, simName);
                        DeleteObject(simName);
                        KillDataFilterRecord(filterID);
                        errMsg = empty;
                      }
                    }
                    DeleteObject(ds);
                    if (rootStates) {
                      DeleteObject(rootStates);
                    }

                    if (errMsg.sLength == 0) {
                      return;
                    }
                  }
                } else {
                  errMsg = "Can't use spool to file option in Simulate when "
                           "the tree depends on category variables.";
                }
              } else {
                errMsg =
                    *(_String *)parameters(2) & " must be an existing matrix";
              }
            } else {
              errMsg = *(_String *)parameters(1) & " must be an existing tree";
            }
          } else {
            errMsg = "Invalid unit length specification (must be >=1)";
          }
        } else {
          errMsg = "Invalid alphabet character specification";
        }
      }
    }
    if (errMsg.sLength == 0) {
      errMsg = _String("Alphabet specification variable ") &
               *(_String *)parameters(3) &
               " must be a string matrix with 2 rows and at least 2 columns";
    }
  } else {
    long i = 0;
    if (!alphabet) {
      i = 3;
    } else {
      if (!treeVar) {
        i = 1;
      } else if (!freqVar) {
        i = 2;
      }
    }

    errMsg = _String("Variable ") & *(_String *)parameters(i) &
             " has not been defined";

  }

  if (errMsg.sLength) {
    errMsg = errMsg & " in Simulate.";
    WarnError(errMsg);
  }

}

//______________________________________________________________________________
bool _ElementaryCommand::ExecuteSimpleStatement (_ExecutionList& chain) {
  _Formula *f = (_Formula*) simpleParameters.GetElement(0L);  
  _PMathObj statement_result = f->Compute(0L, chain.GetExecutionContext());
  if (statement_result) {
    printf ("_ElementaryCommand::ExecuteSimpleStatement %s\n", _String((_String*)statement_result->toStr()).sData);
  }
  return statement_result != NULL;
}

//______________________________________________________________________________
bool _ElementaryCommand::ExecuteJumpStatement (_ExecutionList& chain) {
  bool take_the_jump = true;
  
  if (simpleParameters.countitems() > 1) {
    _Formula *f = (_Formula*) simpleParameters.GetElement(1L);  
    _PMathObj statement_result = f->Compute(0L, chain.GetExecutionContext(), nil, NUMBER);
    if (statement_result != NULL) {
      take_the_jump = CheckEqual(statement_result->Value(), 0.0);
    } else {
      return false;
    }
  }
  if (take_the_jump) {
    chain.currentCommand = simpleParameters.GetElement(0L);
  } else {
      chain.currentCommand++;
  }
  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::Execute(
    _ExecutionList &chain) // perform this command in a given list
    {
  _String errMsg;

        //printf ("Executing command code %ld, simpleParameters[0] = %ld\n", code, simpleParameters.GetElement(0L));

  switch (code) {

  case _HY_HBL_COMMAND_SIMPLE_STATEMENT: // a simple statement
    chain.currentCommand++;
    return ExecuteSimpleStatement(chain);
  case _HY_HBL_COMMAND_JUMP_STATEMENT: // a jump command
    return ExecuteJumpStatement(chain);
  default:
    chain.currentCommand++;

  }
  
  return false;

}

//______________________________________________________________________________
_String _ElementaryCommand::FindNextCommand(_String &input, bool useSoftTrim) {

  bool isString = false, skipping = false;
  char isComment = 0;
  long scopeIn = 0, matrixScope = 0, parenIn = 0, bracketIn = 0, index,
       saveSI = _String::storageIncrement;

  _SimpleList isDoWhileLoop;

  if (input.sLength / 4 > saveSI) {
    _String::storageIncrement = input.sLength / 4;
  }

  _String result(128L, true);

  char lastChar = 0;

  index = input.Length();

  if (!index) {
    result.Finalize();
    return empty;
  }

  // non printable characters at the end ?
  while (index >= 0 && !isprint(input[--index]))
    ;
  input.Trim(0, index, useSoftTrim);

  for (index = 0; index < input.Length(); index++) {
    char c = input.sData[index];

    if (!isString && c == '\t') {
      c = ' ';
    }

    // check for comments
    if (isComment) {
      if (isComment == 1) {
        if (c == '/' && input.sData[index - 1] == '*') {
          isComment = 0;
        }
      } else if (c == '\r' || c == '\n') {
        isComment = 0;
      }

      lastChar = 0;
      continue;
    } else {
      if (!isString && c == '/') {
        switch (input.getChar(index + 1)) {
        case '*':
          isComment = 1;
          break;
        case '/':
          isComment = 2;
        }

        if (isComment) {
          lastChar = 0;
          index++;
          continue;
        }
      }
    }

    // skip spaces
    if (!isString && isspace(c)) {
      if (index >= 6 && input.getChar(index - 1) == 'n' &&
          input.getChar(index - 2) == 'r' && input.getChar(index - 3) == 'u' &&
          input.getChar(index - 4) == 't' && input.getChar(index - 5) == 'e' &&
          input.getChar(index - 6) == 'r') {
        if (index == 6 || (index > 6 && !(isalnum(input.getChar(index - 7)) ||
                                          input.getChar(index - 7) == '_'))) {
          result << ' ';
        }
      }

      skipping = true;
      continue;
    }

    if (skipping && (isalpha(c) || c == '_') &&
        (isalnum(lastChar) || lastChar == '_')) {
      result << ' ';
    }

    skipping = false;

    result << c;

    if (isString && c == '\\') {
      result << input.getChar(++index);
      continue;
    }

    // are we inside a string literal?
    if (c == '"') {
      isString = !isString;
      lastChar = 0;
      continue;
    }

    if (isString) {
      continue;
    }

    // maybe we are done?
    if (c == ';' && scopeIn == 0 && matrixScope == 0 && parenIn <= 0 &&
        bracketIn <= 0) {
      break;
    }

    // check to see whether we are defining a matrix
    if (c == '(') {
      parenIn++;
      lastChar = 0;
      continue;
    }

    if (c == ')') {
      parenIn--;
      if (parenIn < 0) {
        WarnError(_String("Too many closing ')' near '") &
                  input.Cut(MAX(0, index - 32), index) & "'.");
        input = empty;
        return empty;
      }
      lastChar = 0;
      continue;
    }

    if (c == '[') {
      bracketIn++;
      lastChar = 0;
      continue;
    }

    if (c == ']') {
      bracketIn--;
      lastChar = 0;
      continue;
    }

    if (c == '{') {
      if (matrixScope) {
        matrixScope++;
      } else if (lastChar == '=') { // a matrix def
        matrixScope++;
      } else {
        scopeIn++;
        if (index >= 2) {
          long t = input.FirstNonSpaceIndex(0, index - 1, -1);
          if (t >= 1) {
            if (input.getChar(t) == 'o' && input.getChar(t - 1) == 'd') {
              isDoWhileLoop << scopeIn - 1;
              //printf ("%d\n%s\n\n", isDoWhileLoop, input.Cut (t,-1).sData);
            }
          }
        }
      }
      lastChar = 0;
      continue;
    }

    if (c == '}') {
      if (matrixScope) {
        matrixScope--;
      } else {
        scopeIn--;
        if (!parenIn && !bracketIn) {
          if (scopeIn >= 0 && isDoWhileLoop.lLength &&
              isDoWhileLoop.lData[isDoWhileLoop.lLength - 1] == scopeIn) {
            isDoWhileLoop.Delete(isDoWhileLoop.lLength - 1);
          } else if (scopeIn == 0) {
            break;
          }
        }

      }
      lastChar = 0;
      continue;
    }

    lastChar = c;
  }

  result.Finalize();
  _String::storageIncrement = saveSI;

  if (scopeIn || isString || isComment == 1 || parenIn || matrixScope) {
    if (result != '}') {
      WarnError(
          _String("Expression appears to be incomplete/syntax error. Scope: ") &
          scopeIn & ", paretheses depth: " & parenIn & ", matrix scope: " &
          matrixScope & '.' & matrixScope & '.' &
          (isString ? "In a literal. " : empty) &
          (isComment == 1 ? "In a comment " : empty) & '\n' & input);
      input = empty;
      return empty;
    } else {
      result = empty;
    }
  }

  lastChar = 0;
  while (result.getChar(lastChar) == '{') {
    lastChar++;
  }

  if (lastChar) {
    long index2 = result.sLength - 1;

    while (result[index2] == '}') {
      index2--;
    }

    if (result.sLength - index2 - 1 < lastChar) {
      ReportWarning((_String)("Expression appears to be incomplete/syntax "
                              "error and will be ignored:") & input);
      result.DuplicateErasing(&empty);
    } else {
      result.Trim(lastChar, result.sLength - 1 - lastChar);
    }
  }

  if (index < input.Length() - 1) {
    input.Trim(index + 1, -1, useSoftTrim);
  } else if (useSoftTrim) {
    input.sLength = 0;
  } else {
    input.DuplicateErasing(&empty);
  }

  return result;
}

//______________________________________________________________________________
long _ElementaryCommand::ExtractConditions(_String &source, long startwith,
                                           _List &receptacle, char delimeter,
                                           bool includeEmptyConditions) {
  long parenLevel = 1, lastsemi = startwith, index, quote = 0, curlyLevel = 0;

  for (index = startwith; index < source.sLength; index++) {
    char c = source.sData[index];
    if (quote == 0) {
      if (c == '(') {
        parenLevel++;
        continue;
      }
      if (c == '{') {
        curlyLevel++;
        continue;
      }
      if (c == '}') {
        curlyLevel--;
        continue;
      }
      if (c == ')') {
        parenLevel--;
        if (!parenLevel) {
          break;
        }
        continue;
      }
    }
    if (c == '"') {
      if (index == startwith || source.sData[index - 1] != '\\') {
        quote += quote ? -1 : 1;
      }
      continue;
    }
    if (c == delimeter) {
      if (parenLevel > 1 || quote || curlyLevel) {
        continue;
      }

      _String *term =
          (_String *)checkPointer(new _String(source, lastsemi, index - 1));
      receptacle.AppendNewInstance(term);
      lastsemi = index + 1;
      continue;
    }
  }

  if (includeEmptyConditions || lastsemi <= index - 1) {
    receptacle.AppendNewInstance(new _String(source, lastsemi, index - 1));
  }
  return index + 1;
}

//______________________________________________________________________________
bool _ElementaryCommand::MakeGeneralizedLoop(_String *p1, _String *p2,
                                             _String *p3, bool fb,
                                             _String &source,
                                             _ExecutionList &target) {

  // extract the for enclosure
  long beginning = target.lLength, forreturn = target.lLength, index;

  int success = 1;
  bool hasIncrement = false;

  _SimpleList bc;

  while (1) {

    if (p1 && p1->Length()) { // initialization stage
      forreturn++;
      success *= target.BuildList(*p1, nil, true); // add init step
    }

    // append condition now

    if (!success) {
      break;
    }

    if (fb)
      if (p2 && p2->Length()) { // condition stage
        _ElementaryCommand condition(*p2);
        target &&(&condition);
      }

    if (source.getChar(0) == '{') {
      source.Trim(1, -1);
    }

    if ((success *= target.BuildList(source, &bc)) ==
        0) { // construct the main body
      break;
    }

    if (p3 && p3->Length()) {                      // increment stage
      success *= target.BuildList(*p3, nil, true); // add increment step
      hasIncrement = true;
    }

    if (!success) {
      break;
    }

    if (fb) {
      _ElementaryCommand loopback;
      success *= loopback.MakeJumpCommand(nil, forreturn, 0, target);
      target &&(&loopback);
      if (p2 && p2->Length()) {
        success *= ((_ElementaryCommand *)(target(forreturn)))
            ->MakeJumpCommand(p2, forreturn + 1, target.lLength, target);
      }
    } else {
      if (p2) {
        _ElementaryCommand *loopback = new _ElementaryCommand;
        checkPointer(loopback);
        success *= loopback->MakeJumpCommand(p2, forreturn, target.lLength + 1,
                                             target);
        target.AppendNewInstance(loopback);
      }
    }
    break;

  }

  if (!success) { // clean up
    for (index = beginning; index < target.lLength; index++) {
      target.Delete(beginning);
    }
    return false;
  } else {
    // write out the breaks and continues
    for (index = 0; index < bc.lLength; index++) {
      long loc = bc(index);
      if (loc > 0) { // break
        ((_ElementaryCommand *)(target(loc)))
            ->MakeJumpCommand(nil, target.lLength, 0, target);
      } else { // continue
        ((_ElementaryCommand *)(target(-loc)))->MakeJumpCommand(
            nil, target.lLength - (hasIncrement ? 2 : 1), 0, target);
      }
    }
  }

  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::BuildFor(_String &source, _ExecutionList &target,
                                  _List &pieces) {
  // the for loop becomes this:
  // initialize
  // if (condition) then
  // else go to end+2
  // end
  // increment
  // goto if(condition)
  return MakeGeneralizedLoop((_String *)pieces(0), (_String *)pieces(1),
                             (_String *)pieces(2), true, source, target);
}

//______________________________________________________________________________
bool _ElementaryCommand::BuildWhile(_String &source, _ExecutionList &target,
                                    _List &pieces) {
  return MakeGeneralizedLoop(nil, (_String *)pieces(0), nil, true, source,
                             target);
}

//______________________________________________________________________________
bool _ElementaryCommand::BuildIfThenElse(_String &source,
                                         _ExecutionList &target,
                                         _SimpleList *bc) {
  _List pieces;
  long upto = ExtractConditions(source, 3, pieces), beginning = target.lLength;
  target.lastif << target.lLength;
  int success = 1, intIfs = target.lastif.lLength;

  {
    if (pieces.lLength != 1) {
      WarnError("'if' header makes no sense");
    }

    source.Trim(upto, -1);
    target.AppendNewInstance(new _ElementaryCommand);

    _String nextCommand(FindNextCommand(source));
    success *= target.BuildList(nextCommand, bc, true);

  }

  if (!success) { // clean up
    for (unsigned long index = beginning; index < target.lLength; index++) {
      target.Delete(beginning);
    }
    return false;
  } else {
    _ElementaryCommand *ec = (_ElementaryCommand *)(target(beginning));
    ((_ElementaryCommand *)(target(beginning)))->MakeJumpCommand(
        ((_String *)pieces(0)), beginning + 1,
        (ec->simpleParameters.lLength < 2) ? target.lLength
                                           : ec->simpleParameters(1),
        target);
  }

  while (target.lastif.lLength > intIfs) {
    target.lastif.Delete(target.lastif.lLength - 1);
  }

  return target.BuildList(source, bc, true);

}

//______________________________________________________________________________
bool _ElementaryCommand::BuildDoWhile(_String &source, _ExecutionList &target) {
  long upto = source.FindBackwards(_String('}'), 0, -1);
  if (upto >= 0) {
    _String clipped(source, upto + 1, -1);
    if (clipped.beginswith(blWhile)) {
      source.Trim(blDo.sLength, upto);
      _List pieces;
      ExtractConditions(clipped, blWhile.sLength, pieces);
      if (pieces.lLength != 1) {
        WarnError("Malformed while clause in a do-while loop");
        return false;
      }

      if (!MakeGeneralizedLoop(nil, (_String *)pieces(0), nil, false, source,
                               target)) {
        return false;
      }

      return true;
    }
  }
  WarnError(
      "Could not find a matching 'while' in the definition of a do-while loop");

  return false;
}

//______________________________________________________________________________
bool _ElementaryCommand::ProcessInclude(_String &source,
                                        _ExecutionList &target) {

  _String fileName(source, blInclude.sLength, source.sLength - 2);
  fileName = ProcessLiteralArgument(&fileName, target.nameSpacePrefix);
  if (fileName.sLength == 0) {
    WarnError(_String("#include missing a meaningful filename. Check that "
                      "there is a ';' at the end of the statement. Had ") &
              source.Cut(8, source.sLength - 2));
    return false;
  }

  fileName.ProcessFileName(false, false, (Ptr) target.nameSpacePrefix);
  if (terminateExecution) {
    return false;
  }

  PushFilePath(fileName);
  ReadBatchFile(fileName, target);
  PopFilePath();

  return true;
}

void _ElementaryCommand::addAndClean(_ExecutionList &target, _List *parList,
                                     long parFrom) {
  if (parList)
    for (long k = parFrom; k < parList->lLength; k++) {
      parameters &&(*parList)(k);
    }
  target << this;
  DeleteObject(this);
}

//______________________________________________________________________________
// DataSet    dataSetid = ReadDataFile ("..");
// or
// DataSet    dataSetid = SimulateDataSet (likeFunc);
// or
// DataSet    dataSetid = Concatenate (<purge>,list of DataSets);
// or
// DataSet    dataSetid = Combine (<purge>,list of DataSets);
// or
// DataSet    dataSetid = ReconstructAncestors (lf)
// or
// DataSet    dataSetid = SampleAncestors (lf)
// or
// DataSet    dataSetid = Simulate (tree, freqs, alphabet, <store internal
// nodes, root vector>)
// or
// DataSet    dataSetid = ReadFromString (string);
bool _ElementaryCommand::ConstructDataSet(_String &source,
                                          _ExecutionList &target) {

  // first we must segment out the data set name
  // then the ReadDataFile command
  // then the data set file name
  // look for the data set name first

  long mark1 = source.FirstSpaceIndex(0, -1, 1),
       mark2 = source.Find('=', mark1, -1);

  _String dsID(source, mark1 + 1, mark2 - 1);

  if (mark1 == -1 || mark2 == -1 || dsID.Length() == 0) {
    WarnErrorWhileParsing("DataSet declaration missing a valid identifier",
                          source);
    return false;
  }

  // now look for the opening paren
  mark1 = source.Find('(', mark2, -1);

  _ElementaryCommand dsc;
  _String oper(source, mark2 + 1, mark1 - 1);

  if (oper == _String("ReadDataFile") ||
      oper == _String("ReadFromString")) { // a switch statement if more than 1
    _List pieces;
    mark2 = ExtractConditions(source, mark1 + 1, pieces, ',');
    if (pieces.lLength != 1) {
      WarnErrorWhileParsing("DataSet declaration missing a valid filename",
                            source);
      return false;
    }

    _ElementaryCommand *dsc = makeNewCommand(5);

    dsc->parameters &&(&dsID);
    dsc->parameters &&(pieces(0));

    if (oper == _String("ReadFromString")) {
      dsc->simpleParameters << 1;
    }

    dsc->addAndClean(target);
    return true;

  } else if (oper.Equal(&blSimulateDataSet)) {
    _List pieces;
    mark2 = ExtractConditions(source, mark1 + 1, pieces, ',');
    if (pieces.lLength > 4 || pieces.lLength == 0) {
      WarnErrorWhileParsing(
          blSimulateDataSet &
              "expects 1-4 parameters: likelihood function ident (needed), a "
              "list of excluded states, a matrix to store random rates in, and "
              "a matrix to store the order of random rates in (last 3 - "
              "optional).",
          source);
      return false;
    }

    dsc.code = 12;
    dsc.parameters &&(&dsID);
    dsc.parameters &&(pieces(0));
    for (mark2 = 1; mark2 < pieces.lLength; mark2++) {
      dsc.parameters &&(pieces(mark2));
    }

    target &&(&dsc);
    return true;

  } else if (oper == _String("Concatenate") || oper == _String("Combine")) {
    _List pieces;
    mark2 = ExtractConditions(source, mark1 + 1, pieces, ',');
    if (pieces.lLength == 0) {
      WarnErrorWhileParsing(
          "DataSet merging operation missing a valid list of arguments.",
          source);
      return false;
    }

    dsc.code = 16;
    dsc.parameters &&(&dsID);

    long i = 0;

    dsc.simpleParameters << ((oper == _String("Concatenate")) ? 1 : 2);

    _String purge("purge");
    if (purge.Equal((_String *)pieces(0))) {
      dsc.simpleParameters[0] *= -1;
      i++;
    }

    for (; i < pieces.lLength; i++) {
      dsc.parameters << pieces(i);
    }

    if (dsc.parameters.lLength <= 1) {
      WarnErrorWhileParsing(
          "DataSet merging operation missing a valid list of arguments.",
          source);
      return false;
    }

    target &&(&dsc);
    return true;

  } else {
    if (oper == _String("ReconstructAncestors") ||
        oper == _String("SampleAncestors")) {
      _List pieces;
      mark2 = ExtractConditions(source, mark1 + 1, pieces, ',');
      if (pieces.lLength > 3 || pieces.lLength == 0) {
        WarnErrorWhileParsing(
            "ReconstructAncestors and SampleAncestors expects 1-4 parameters: "
            "likelihood function ident (mandatory), an matrix expression to "
            "specify the list of partition(s) to reconstruct/sample from "
            "(optional), and, for ReconstructAncestors, an optional MARGINAL "
            "flag, plus an optional DOLEAVES flag.",
            source);
        return false;
      }

      dsc.code = (oper == _String("ReconstructAncestors")) ? 38 : 50;
      dsc.parameters &&(&dsID);
      dsc.parameters << pieces(0);
      for (long optP = 1; optP < pieces.lLength; optP++)
        if (((_String *)pieces(optP))->Equal(&marginalAncestors)) {
          dsc.simpleParameters << -1;
        } else if (((_String *)pieces(optP))->Equal(&doLeavesAncestors)) {
          dsc.simpleParameters << -2;
        } else {
          dsc.parameters << pieces(optP);
        }

      target &&(&dsc);
      return true;
    } else if (oper == _String("Simulate")) {
      _List pieces;
      mark2 = ExtractConditions(source, mark1 + 1, pieces, ',');
      if ((pieces.lLength > 7) || (pieces.lLength < 4)) {
        WarnErrorWhileParsing(
            "Simulate expects 4-6 parameters: tree with attached models, "
            "equilibrium frequencies, character map, number of sites|root "
            "sequence, <save internal node sequences>, <file name for direct "
            "storage>",
            source);
        return false;
      }

      dsc.code = 52;
      dsc.parameters &&(&dsID);

      for (mark2 = 0; mark2 < pieces.lLength; mark2++) {
        dsc.parameters &&(pieces(mark2));
      }

      target &&(&dsc);
      return true;
    } else {
      WarnErrorWhileParsing(
          "Expected DataSet ident = ReadDataFile(filename); or DataSet ident = "
          "SimulateDataSet (LikelihoodFunction); or DataSet ident = Combine "
          "(list of DataSets); or DataSet ident = Concatenate (list of "
          "DataSets); or DataSet ident = ReconstructAnscetors (likelihood "
          "function); or DataSet ident = SampleAnscetors (likelihood function) "
          "or DataSet	  dataSetid = ReadFromString (string);",
          source);
    }
  }

  return false;
}

//______________________________________________________________________________
// category <id> = (number of int, weights, method for representation,
// density, cumulative, left bound, right bound);
bool _ElementaryCommand::ConstructCategory(_String &source,
                                           _ExecutionList &target) {

  long mark1 = source.FirstSpaceIndex(0, -1, 1),
       mark2 = source.Find('=', mark1, -1);

  _String catID(source, mark1 + 1, mark2 - 1);

  if (mark1 == -1 || mark2 == -1 || catID.Length() == 0) {
    _String errMsg("Category variable declaration missing a valid identifier");
    WarnError(errMsg);
    return false;
  }

  // now look for the opening paren
  mark1 = source.Find('(', mark2, -1);

  if (mark1 != -1) {
    mark2 = source.FindBackwards(')', mark1 + 1, -1);
    if (mark2 != -1) {
      source = source.Cut(mark1 + 1, mark2 - 1);
      _List args;
      mark2 = ExtractConditions(source, 0, args, ',');
      if (args.lLength >= 7) {
        _ElementaryCommand *cv = new _ElementaryCommand(20);
        checkPointer(cv);
        cv->parameters &&(&catID);
        cv->addAndClean(target, &args, 0);
        return true;
      }
    }
  }
  _String errMsg(
      "Expected: category <id> = (number of intervals, weights, method for "
      "representation, density, cumulative, left bound, right bound,<optional "
      "mean cumulative function>,<optional hidden markov matrix>);");
  WarnError(errMsg);
  return false;
}

//______________________________________________________________________________
// UseMatrix (matrixIdent)
bool _ElementaryCommand::ConstructStateCounter(_String &source,
                                               _ExecutionList &target) {
  _List args;
  ExtractConditions(source, blStateCounter.sLength, args, ',');
  if (args.lLength != 2) {
    WarnError("Expected: StateCounter(likefuncID, callback function ID)");
    return false;
  }
  _ElementaryCommand *sc = new _ElementaryCommand(47);
  sc->addAndClean(target, &args, 0);
  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::ConstructChoiceList(_String &source,
                                             _ExecutionList &target) {
  _List args;
  ExtractConditions(source, blChoiceList.sLength, args, ',');
  if (args.lLength < 5) {
    WarnError("ChoiceList needs at least 5 arguments");
    return false;
  }
  _ElementaryCommand *cv = new _ElementaryCommand(32);

  cv->parameters << args(0);
  ((_String *)args.lData[1])->StripQuotes();
  cv->parameters << args(1);
  cv->parameters << args(2);
  cv->parameters << args(3);

  if (args.lLength > 5) {
    _List choices;
    for (long k = 4; k < args.lLength - 1; k += 2) {
      ((_String *)args.lData[k])->StripQuotes();
      ((_String *)args.lData[k + 1])->StripQuotes();
      _List thisChoice;
      thisChoice << args(k);
      thisChoice << args(k + 1);
      choices &&&thisChoice;
    }
    cv->parameters &&&choices;
    cv->simpleParameters << 0;
  } else {
    cv->parameters << args(4);
    cv->simpleParameters << 1;
  }

  cv->addAndClean(target, nil, 0);
  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::ConstructReplicateConstraint(_String &source,
                                                      _ExecutionList &target) {

  // ReplicateConstraint ("constraint to be replicated in terms of
  // this1,...,thisn and"
  // list of n variables to put in place of this1,
  // this2, ... thisn);
  // this1 .. etc are all expected to be either trees of nodes of trees with
  // wildcards.

  _List args;
  ExtractConditions(source, 20, args, ',');
  if (args.lLength < 2) {
    _String errMsg(
        "Expected: ReplicateConstraint (\"constraint to be replicated in terms "
        "of this1,...,thisn and wildcard *\", list of n variables to put in "
        "place of this1, this2, ... thisn);");
    acknError(errMsg);
    return false;
  }

  /*_String *theConstraint = (_String*)args(0), thisString;
    long k = 0;
    theConstraint->StripQuotes();
    do
    {
        k++;
        thisString  = _String("this")&_String(k);
    }
    while (theConstraint->Find(thisString)!=-1);

    if (args.lLength!=k)
    {
        _String errMsg ("Replicate constraint could not match the number of
  'this' arguments with actual variables");
        acknError (errMsg);
        return false;
    }*/

  _ElementaryCommand cv;
  cv.code = 26;

  for (long k = 0; k < args.lLength; k++) {
    cv.parameters << args(k);
  }

  target &&&cv;
  return true;
}

//______________________________________________________________________________

bool _ElementaryCommand::ConstructTree(_String &source, _ExecutionList &target) {
  // Tree   treeid = (...) or Topology = (...);
  long mark1 = source.FirstSpaceIndex(0, -1, 1), mark2, mark3;
  mark2 = source.Find('=', mark1, -1);
  mark3 = mark2;

  if ((mark1 == -1) || (mark2 == -1) || (mark1 + 1 > mark2 - 1)) {
    _String errMsg("Tree declaration missing a valid identifier");
    acknError(errMsg);
    return false;
  }

  _String dsID = source.Cut(mark1 + 1, mark2 - 1);
  // now look for the opening paren

  mark1 = source.Find('(', mark2, -1);
  mark2 = source.FindBackwards(')', mark1, -1);
  if ((mark1 == -1) || (mark2 == -1) || (mark2 < mark1)) {
    if (source.Find(getDString) == -1) {
      mark1 = mark3 + 1;
      mark2 = source.Find(';', mark3, -1) - 1;
    } else {
      source = getDString;
      mark1 = 0;
      mark2 = -1;
    }
  }

  _ElementaryCommand *dsc =
      new _ElementaryCommand(source.startswith(blTree) ? 7 : 54);
  checkPointer(dsc);
  dsc->parameters &&(&dsID);
  dsc->parameters.AppendNewInstance(new _String(source, mark1, mark2));
  dsc->addAndClean(target, nil, 0);
  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::ConstructDataSetFilter(_String &source,
                                                _ExecutionList &target) {

  // DataSetFilter      dataSetFilterid = CreateFilter
  // (datasetid;unit;vertical partition; horizontal partition; alphabet
  // exclusions);
  // first we must segment out the data set name

  long mark1 = source.FirstSpaceIndex(0, -1, 1),
       mark2 = source.Find('=', mark1, -1);

  _String dsID(source, mark1 + 1, mark2 - 1), command;

  if (mark1 == -1 || mark2 == -1 || dsID.Length() == 0) {
    _String errMsg("DataSetFilter declaration missing a valid identifier");
    acknError(errMsg);
    return false;
  }

  // now look for the opening paren
  mark1 = source.Find('(', mark2, -1);
  command = source.Cut(mark2 + 1, mark1 - 1);
  _ElementaryCommand *dsf;
  _List pieces;

  if (command == _String("CreateFilter")) {
    dsf = new _ElementaryCommand(6);
  } else if (command == _String("Permute")) {
    dsf = new _ElementaryCommand(27);
  } else if (command == _String("Bootstrap")) {
    dsf = new _ElementaryCommand(28);
  } else {
    _String errMsg("Expected: DataSetFilter	  dataSetFilterid = CreateFilter "
                   "(datasetid,unit,vertical partition,horizontal "
                   "partition,alphabet exclusions); or Permute/Bootstrap "
                   "(dataset/filter,<atom>,<column partition>)");
    acknError(errMsg);
    return false;
  }

  mark2 = ExtractConditions(source, mark1 + 1, pieces, ',');
  if (!(pieces.lLength >= 2 || (pieces.lLength == 1 && dsf->code == 6))) {
    _String errMsg("Parameter(s) missing in DataSetFilter definition.");
    acknError(errMsg);
    return false;
  }

  dsf->parameters &&(&dsID);
  dsf->addAndClean(target, &pieces);
  return true;
}

//______________________________________________________________________________
// Model ID = (inst transition matrix ident, equilibrium frequencies ident,
// <multiply by frequencies>);
// if the third parameter is explicitFormMExp, then inst transition matrix
// ident is expected to be an explicit matrix exponential
// EXPRESSION
bool _ElementaryCommand::ConstructModel(_String &source, _ExecutionList &target) {

  // first we must segment out the data set name
  long mark1 = source.FirstSpaceIndex(0, -1, 1),
       mark2 = source.Find('=', mark1, -1);

  _String modelID(source, mark1 + 1, mark2 - 1);

  if (mark1 == -1 || mark2 == -1 || !modelID.IsValidIdentifier()) {
    _String errMsg("Model declaration missing a valid identifier.");
    acknError(errMsg);
    return false;
  }

  // now look for the opening paren
  mark1 = source.Find('(', mark2, -1);
  _List pieces;
  mark2 = ExtractConditions(source, mark1 + 1, pieces, ',');

  if (pieces.lLength < 2) {
    _String errMsg("Parameter(s) missing in Model definition. Must have a "
                   "matrix and a compatible eqiulibrium frequencies vector.");

    acknError(errMsg);
    return false;

  } else {
    if (pieces.lLength > 3) {
      _String errMsg("Too many parameters (3 max) in Model definition");
      acknError(errMsg);
      return false;
    }
  }

  _ElementaryCommand *model = new _ElementaryCommand(31);
  checkPointer(model);
  model->parameters &&(&modelID);
  model->addAndClean(target, &pieces, 0);
  return true;

}

//______________________________________________________________________________
/*bool    _ElementaryCommand::ConstructFprintf (_String&source,
_ExecutionList&target)
{

    _ElementaryCommand  *fpr = (_ElementaryCommand*)checkPointer(new
_ElementaryCommand (8));

    long     lastStart = 8;
    bool     done      = false;

    _String  comma (",");

    while (!done) {
        long lastEnd = source.FindTerminator(lastStart,comma);
        if (lastEnd < 0) {
            lastEnd = source.sLength-2;
            done = true;
        }
        _String *thisArgument = new _String (source, lastStart, lastEnd-1);

        if (fpr->parameters.lLength && thisArgument->IsALiteralArgument(true)) {
            fpr->simpleParameters << fpr->parameters.lLength;
            _FString converted (*thisArgument, true);
            fpr->parameters << converted.theString;
            DeleteObject (thisArgument);
        } else {
            fpr->parameters.AppendNewInstance (thisArgument);
        }
        lastStart = lastEnd + 1;
    }

    fpr->addAndClean(target, nil, 0);
    return true;
}*/

//______________________________________________________________________________
// syntax:
// fscanf (stdin or "file name" or PROMPT_FOR_FILE, "argument descriptor",
// list of arguments to be read);
// argument descriptor is a comma separated list of one of the three
// constants
// "number", "matrix", "tree"
// list of arguments to be read specifies which variables will receive the
// values
bool _ElementaryCommand::ConstructFscanf(_String &source,
                                         _ExecutionList &target) {

  if (!allowedFormats.lLength) {
    allowedFormats.AppendNewInstance(new _String("Number"));
    allowedFormats.AppendNewInstance(new _String("Matrix"));
    allowedFormats.AppendNewInstance(new _String("Tree"));
    allowedFormats.AppendNewInstance(new _String("String"));
    allowedFormats.AppendNewInstance(new _String("NMatrix"));
    allowedFormats.AppendNewInstance(new _String("Raw"));
    allowedFormats.AppendNewInstance(new _String("Lines"));
  }

  _ElementaryCommand *fscan =
      new _ElementaryCommand(source.startswith(blsscanf) ? 56 : 25);
  _List arguments, argDesc;
  long f, p, shifter = 0;
  ExtractConditions(source, 7, arguments, ',');

  if (arguments.lLength < 3) {
    WarnError(_String("Too few arguments in call to fscanf or sscanf"));
    DeleteObject(fscan);
    return false;
  }

  fscan->parameters << arguments(0);

  if (((_String *)arguments(1))->Equal(&blScanfRewind)) {
    fscan->simpleParameters << -1;
    shifter = 1;
  }

  ((_String *)arguments(1 + shifter))->StripQuotes();
  ExtractConditions(*((_String *)arguments(1 + shifter)), 0, argDesc, ',');

  for (f = 0; f < argDesc.lLength; f++) {
    p = allowedFormats.Find(argDesc(f));
    if (p == -1) {
      WarnError(
          *((_String *)argDesc(f)) &
          " is not a valid type descriptor for fscanf. Allowed ones are:" &
          _String((_String *)allowedFormats.toStr()));
      DeleteObject(fscan);
      return false;
    } else {
      fscan->simpleParameters << p;
    }
  }

  if (arguments.lLength != fscan->simpleParameters.lLength + 2) {
    WarnError(_String("fscanf passed ") &
              _String((long)(fscan->simpleParameters.lLength - shifter)) &
              " parameter type descriptors and " &
              _String((long)(arguments.lLength - 2 - shifter)) &
              " actual arguments");
    DeleteObject(fscan);
    return false;
  }

  for (f = 2 + shifter; f < arguments.lLength; f++) {
    _String *thisArg = (_String *)arguments(f);
    if (thisArg->IsValidIdentifier()) {
      fscan->parameters << thisArg;
    } else {
      WarnError(_String("fscanf passed an invalid variable identifier: ") &
                *thisArg);
      DeleteObject(fscan);
      return false;
    }
  }

  fscan->addAndClean(target, nil, 0);
  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::MakeJumpCommand(_String *source, long branch1,
                                         long branch2,
                                         _ExecutionList &parentList) {
  long oldFla = 0;
  code = 4;

  if (simpleParameters.lLength == 3) {
    if (source) {
      _Formula *f = (_Formula *)simpleParameters(2);
      delete (f);
    } else {
      oldFla = simpleParameters(2);
    }
  }

  if (branch1 == -1) {
    if (simpleParameters.lLength == 0) {
      WarnError("An if-then-else scoping error. Check opening and closing "
                "brackets and double else's.");
      return false;
    }
    branch1 = simpleParameters[0];
  }

  simpleParameters.Clear();
  simpleParameters << branch1;
  simpleParameters << branch2;
  if (source) {
    /*_Formula f;
    long status = Parse (&f, *source, parentList.nameSpacePrefix,nil);
    if (status==-1)
        simpleParameters<<long(f.makeDynamic());*/
    parameters &&source;
  } else if (oldFla) {
    simpleParameters << oldFla;
  }
  return true;
}

//______________________________________________________________________________
// syntax: FindRoot (receptacle, expression, variable, left bound, right
// bound)
// or: Integrate (receptacle, expression, variable, left bound, right bound
bool _ElementaryCommand::ConstructFindRoot(_String &source,
                                           _ExecutionList &target) {
  _List pieces;
  long mark1 = source.Find('(');
  _String oper(source, 0, mark1);
  source.Trim(ExtractConditions(source, mark1 + 1, pieces, ','), -1);
  if (pieces.lLength != 5) {
    WarnError("Expected: FindRoot|Integrate (receptacle, expression, variable, "
              "left bound, right bound).");
    return false;
  }

  _ElementaryCommand *fri =
      new _ElementaryCommand(oper.Equal(&blFindRoot) ? 43 : 48);
  fri->addAndClean(target, &pieces, 0);
  return true;

}

//______________________________________________________________________________
// syntax: MPISend (numeric node ID, string with HBL code <or> a LF ID,
// <input redirect target>);
bool _ElementaryCommand::ConstructMPISend(_String &source,
                                          _ExecutionList &target) {

  _List pieces;
  ExtractConditions(source, blMPISend.sLength, pieces, ',');
  if (pieces.lLength != 2 && pieces.lLength != 3) {
    WarnError("Expected: MPISend (numeric node ID, string with HBL code <or> a "
              "LF ID).");
    return false;
  }
  _ElementaryCommand *mpiSend = makeNewCommand(44);
  mpiSend->addAndClean(target, &pieces, 0);
  return true;
}

//______________________________________________________________________________
// syntax: MPIReceive (can receive from node, received from node, receptacle
// for the string result);
bool _ElementaryCommand::ConstructMPIReceive(_String &source,
                                             _ExecutionList &target) {
  _List pieces;
  ExtractConditions(source, blMPIReceive.sLength, pieces, ',');
  if (pieces.lLength != 3) {
    WarnError("Expected: MPIReceive (can receive from node, received from "
              "node, receptacle for the string result).");
    return false;
  }

  _ElementaryCommand *mpiRecv = makeNewCommand(45);
  mpiRecv->addAndClean(target, &pieces, 0);
  return true;
}

//______________________________________________________________________________
// syntax: ConstructCategoryMatrix (receptacle, likelihood function,
// [COMPLETE/SHORT, which partitions to include -- a matrix agrument] )
bool _ElementaryCommand::ConstructCategoryMatrix(_String &source,
                                                 _ExecutionList &target) {
  _List pieces;
  ExtractConditions(source, blConstructCM.sLength, pieces, ',');
  if (pieces.lLength < 2) {
    WarnError("Expected: ConstructCategoryMatrix (receptacle, likelihood "
              "function,COMPLETE/SHORT/WEIGHTS [optional; default is "
              "COMPLETE], [optional matrix argument with partitions to "
              "include; default is to include all]");
    return false;
  }

  _ElementaryCommand *constuctCatMatrix = makeNewCommand(21);
  constuctCatMatrix->addAndClean(target, &pieces, 0);
  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::ConstructExecuteCommands(_String &source,
                                                  _ExecutionList &target) {

  _List pieces;
  char execAFile = source.startswith(blExecuteAFile)
                       ? 1
                       : (source.startswith(blLoadFunctionLibrary) ? 2 : 0);

  long code = 39;

  switch (execAFile) {
  case 0:
    ExtractConditions(source, blExecuteCommands.sLength, pieces, ',');
    break;

  case 1:
    ExtractConditions(source, blExecuteAFile.sLength, pieces, ',');
    code = 62;
    break;

  case 2:
    ExtractConditions(source, blLoadFunctionLibrary.sLength, pieces, ',');
    code = 66;
    break;
  }

  if (pieces.lLength < 1 || pieces.lLength > 3) {
    WarnError(
        "Expected: ExecuteCommands (identifier, <compiled|(input "
        "redirect<,string prefix>)>) or ExecuteAFile (path name, "
        "<compiled|(input redirect<,string prefix>)> or LoadFunctionLibrary "
        "(path name, <compiled|(input redirect<,string prefix>)>)");
    return false;
  }

  _ElementaryCommand *exc =
      (_ElementaryCommand *)checkPointer(new _ElementaryCommand(code));

  exc->parameters << pieces(0);

  if (pathNames.lLength) {
    exc->parameters &&pathNames(pathNames.lLength - 1);
  } else {
    exc->parameters &&&empty;
  }

  if (pieces.lLength > 1) {
    if (*(_String *)pieces(1) == _String("compiled")) {
      exc->simpleParameters << 1;
    } else {
      exc->parameters << pieces(1);
      if (pieces.lLength > 2) {
        exc->parameters << pieces(2);
      }
    }
  }

  exc->addAndClean(target);
  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::ConstructOpenWindow(_String &source,
                                             _ExecutionList &target) {

  _List pieces;
  ExtractConditions(source, blOpenWindow.sLength, pieces, ',');
  if (pieces.lLength < 2 || pieces.lLength > 3) {
    WarnError("Expected: OpenWindow (window type,parameter matrix,<position "
              "string>)");
    return false;
  }

  if (pieces.lLength == 3) {
    ((_String *)pieces(2))->StripQuotes();
  }

  _ElementaryCommand *exc = new _ElementaryCommand(40);
  exc->addAndClean(target, &pieces, 0);
  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::ConstructSpawnLF(_String &source,
                                          _ExecutionList &target) {
  _List pieces;
  ExtractConditions(source, blSpawnLF.sLength, pieces, ',');
  if (pieces.lLength != 4) {
    _String errMsg("Expected: SpawnLikelihoodFunction (likeFuncID, treeID, "
                   "window ID, subset matrix)");
    acknError(errMsg);
    return false;
  }

  _ElementaryCommand *exc = new _ElementaryCommand(41);
  exc->addAndClean(target, &pieces, 0);
  return true;
}

//______________________________________________________________________________
// syntax: GetDataInfo(matrixID, dataFilterID,<sequence ref, site ref |
// sequence 1 , sequence 2, DISTANCES>)
bool _ElementaryCommand::ConstructGetDataInfo(_String &source,
                                              _ExecutionList &target) {
  _List pieces;
  ExtractConditions(source, blGetDataInfo.sLength, pieces, ',');
  if (pieces.lLength < 2 || pieces.lLength > 5) {
    WarnError("Expected: syntax: GetDataInfo(matrix ID, dataFilterID,<sequence "
              "ref, site ref | sequence 1 , sequence 2, DISTANCES>)");
    return false;
  }

  _ElementaryCommand *gdi = new _ElementaryCommand(46);
  gdi->addAndClean(target, &pieces, 0);
  return true;
}

//______________________________________________________________________________
// syntax: OpenDataPanel(dataSetID,"species order","display
// settings","partition settings")

bool _ElementaryCommand::ConstructOpenDataPanel(_String &source,
                                                _ExecutionList &target) {
  _List pieces;
  ExtractConditions(source, blOpenDataPanel.sLength, pieces, ',');
  if (pieces.lLength != 4 && pieces.lLength != 5) {
    ReportWarning(
        "Expected: syntax: OpenDataPanel(dataSetID,\"species order\",\"display "
        "settings\",\"partition settings\"),[likefunc ID]");
    return false;
  }

  _ElementaryCommand *sp = new _ElementaryCommand(36);
  sp->addAndClean(target, &pieces, 0);
  return true;
}

//______________________________________________________________________________
// syntax: GetInformation(object,receptacle)
bool _ElementaryCommand::ConstructGetInformation(_String &source,
                                                 _ExecutionList &target) {

  _List pieces;
  ExtractConditions(source, blGetInformation.sLength, pieces, ',');

  if (pieces.lLength < 2) {
    _String errMsg("Expected at least 2 arguments: "
                   "GetInformation(object,receptacle,...);");
    WarnError(errMsg);
    return false;
  } else {

    _String *s0 = (_String *)pieces(0), *s1 = (_String *)pieces(1);

    if (!(s0->IsValidIdentifier() &&
          ((s1->sLength > 2 &&
            s1->getChar(0) == '"' & s1->getChar(s1->sLength - 1) == '"') ||
           s1->IsValidIdentifier()))) {
      WarnError(_String("Both ") & *s0 & " and " & *s1 &
                " must be valid identifiers in call to GetInformation.");
      return false;
    }
  }

  _ElementaryCommand *sp = makeNewCommand(37);
  sp->addAndClean(target, &pieces, 0);
  return true;
}

//______________________________________________________________________________
// syntax: LikelihoodFunction id = (filter1, tree1, ..., filterN, treeN,
// optional compute template)
// or LikelihoodFunction3 id = (filter1, tree1, freq1, ... filterN, treeN,
// freqN, optional compute template)
bool _ElementaryCommand::ConstructLF(_String &source, _ExecutionList &target) {
  long mark1 = source.FirstSpaceIndex(0, -1, 1),
       mark2 = source.Find('=', mark1, -1);

  if (mark1 == -1 || mark2 == -1 || mark1 + 1 > mark2 - 1) {
    _String errMsg(
        "Likelihood function declaration missing a valid identifier");
    acknError(errMsg);
    return false;
  }

  _String lfID(source, mark1 + 1, mark2 - 1);

  // now look for the opening paren
  _List pieces;
  mark1 = source.Find('(', mark2, -1);
  mark2 = source.FindBackwards(')', mark1, -1);
  ExtractConditions(source, mark1 + 1, pieces, ',');

  if (mark1 == -1 || mark2 == -1 || mark2 < mark1) {
    WarnError(
        "Expected: Likelihood Function ident = (tree1, datasetfilter1,...)");
    return false;
  }

  _ElementaryCommand *dsc =
      (_ElementaryCommand *)checkPointer(new _ElementaryCommand(11));
  dsc->parameters &&(&lfID);

  if (source.startswith(blLF3)) {
    dsc->simpleParameters << 1;
  }

  dsc->addAndClean(target, &pieces, 0);
  return true;
}

//______________________________________________________________________________
// syntax: function <ident> (comma separated list of parameters) {body}
bool _ElementaryCommand::ConstructFunction(_String &source,
                                           _ExecutionList &chain) {
  if (isInFunction) {
    WarnError("Nested function declarations are not allowed");
    return false;
  }

  isInFunction = true;

  bool isFFunction = source.beginswith(blFFunction),
       isLFunction = source.beginswith(blLFunction);

  long mark1 = source.FirstNonSpaceIndex(
                         (isFFunction || isLFunction) ? blFFunction.sLength
                                                      : blFunction.sLength,
                         -1, 1),
       mark2 = source.Find('(', mark1, -1);

  if (mark1 == -1 || mark2 == -1 || mark1 + 1 > mark2 - 1) {

    WarnError("Function declaration missing a valid function identifier or "
              "parameter list.");

    isInFunction = false;
    return false;
  }

  _String *funcID =
      (_String *)checkPointer(new _String(source.Cut(mark1, mark2 - 1)));
  *funcID = chain.AddNameSpaceToID(*funcID);

  // now look for the opening paren
  if ((mark1 = batchLanguageFunctionNames.Find(funcID)) != -1) {
    ReportWarning(_String("Overwritten previously defined function:'") &
                  *funcID & '\'');
  }

  _List pieces;
  long upto = ExtractConditions(source, mark2 + 1, pieces, ',', false);
  if (upto == source.sLength || source[upto] != '{' || source[source.sLength - 1] != '}') {
    WarnError(
        _String("Function declaration is missing a valid function body."));
    isInFunction = false;
    return false;
  }

  _String extraNamespace;
  if (isLFunction) {
    extraNamespace = _HYGenerateANameSpace();
  }

  for (long k = 0; k < pieces.lLength; k++) {
    pieces.Replace(k, new _String(chain.AddNameSpaceToID(*(_String *)pieces(k),
                                                         &extraNamespace)),
                   false);
  }

  _String sfunctionBody(source, upto + 1, source.Length() - 2);
  _ExecutionList *functionBody;
  if (isLFunction) {
    _String *existing_namespace = chain.GetNameSpace();
    if (existing_namespace) {
      extraNamespace = *existing_namespace & '.' & extraNamespace;
    }
    functionBody = new _ExecutionList(sfunctionBody, &extraNamespace, true);
  } else {
    functionBody =
        new _ExecutionList(sfunctionBody, chain.GetNameSpace(), true);
  }

  //  take care of all the return statements
  while (returnlist.lLength) {
    ((_ElementaryCommand *)(*functionBody)(returnlist(0)))->simpleParameters
        << functionBody->lLength;
    returnlist.Delete(0);
  }

  if (mark1 >= 0) {
    batchLanguageFunctions.Replace(mark1, functionBody, false);
    batchLanguageFunctionNames.Replace(mark1, funcID, false);
    batchLanguageFunctionParameterLists.Replace(mark1, &pieces, true);
    batchLanguageFunctionParameters.lData[mark1] = pieces.lLength;
    batchLanguageFunctionClassification.lData[mark1] =
        isFFunction ? BL_FUNCTION_NORMAL_UPDATE : BL_FUNCTION_ALWAYS_UPDATE;
  } else {
    batchLanguageFunctions.AppendNewInstance(functionBody);
    batchLanguageFunctionNames.AppendNewInstance(funcID);
    batchLanguageFunctionParameterLists &&(&pieces);
    batchLanguageFunctionParameters << pieces.lLength;
    batchLanguageFunctionClassification
        << (isFFunction ? BL_FUNCTION_NORMAL_UPDATE
                        : BL_FUNCTION_ALWAYS_UPDATE);
  }

  isInFunction = false;
  return true;

}

//______________________________________________________________________________
// syntax: return <statement>
bool _ElementaryCommand::ConstructReturn(_String &source,
                                         _ExecutionList &target) {

  /*if (!isInFunction)
  {
      _ElementaryCommand exit;
      exit.code = 4;
      exit.simpleParameters<<-1;
      target&&(&exit);
      return true;
  }*/

  long mark1 = source.FirstNonSpaceIndex(blReturn.sLength, -1, 1);
  _ElementaryCommand ret;
  ret.code = 14;

  if (mark1 != -1) {
    _String cut_s;
    if (source.sData[source.sLength - 1] == ';') {
      cut_s = source.Cut(mark1, source.sLength - 2);
    } else {
      cut_s = source.Cut(mark1, -1);
    }
    ret.parameters &&(&cut_s);
  }

  if (isInFunction) {
    returnlist << target.lLength;
  } else {
    ret.simpleParameters << -1;
  }

  target &&(&ret);
  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::ExtractValidateAddHBLCommand(
    _String &current_stream, const long command_code, _List *pieces,
    _HBLCommandExtras *command_spec, _ExecutionList &command_list) {
  if (command_spec->is_assignment) {
    // TBA
  } else {
    // by default push all of the 'pieces' arguments to the "argument" list
    _ElementaryCommand *cv = new _ElementaryCommand(command_code);
    cv->addAndClean(command_list, pieces, 0);
  }

  return true;
}

//______________________________________________________________________________
// syntax: DoSQL (dbID,action string|file name,<callback ID>)
bool _ElementaryCommand::ConstructDoSQL(_String &source, _ExecutionList &target)
{
  _List pieces;
  ExtractConditions(source, blDoSQL.sLength, pieces, ',');
  if (pieces.lLength != 3) {
    WarnError(_String("Expected syntax:") & blDoSQL & "(dbID|" & sqlOpen & '|' &
              sqlClose & ",transaction string|file name,callback ID for an SQL "
                         "transaction|where to store DB numeric ID)");
    return false;
  }

  _ElementaryCommand *dsql = new _ElementaryCommand(53);
  dsql->addAndClean(target, &pieces, 0);
  return true;
}

//______________________________________________________________________________
// syntax: #profile START|PAUSE|RESUME|indetifier to dump in
bool _ElementaryCommand::ConstructProfileStatement(_String &source,
                                                   _ExecutionList &target) {

  _List pieces;
  ExtractConditions(source, blHBLProfile.sLength + 1, pieces, ';');
  if (pieces.lLength != 2) {
    WarnError(_String("Expected syntax:") & blHBLProfile &
              " START|PAUSE|RESUME|where to store)");
    return false;
  }

  _ElementaryCommand *sp = new _ElementaryCommand(58);
  sp->addAndClean(target, &pieces, 0);
  return true;
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteDataFilterCases(_ExecutionList &chain) {
  chain.currentCommand++;

  _String dataObjectID = chain.AddNameSpaceToID(*(_String *)parameters(1));

  long dsID = (parameters.lLength > 2) ? FindDataSetName(dataObjectID) : -1;
  bool isFilter = false;

  if (dsID == -1) {
    dsID = (parameters.lLength > 2) ? FindDataSetFilterName(dataObjectID) : -1;
    if (dsID == -1) {
      _AssociativeList *numericFilter =
          (_AssociativeList *)FetchObjectFromVariableByType(&dataObjectID,
                                                            ASSOCIATIVE_LIST);
      if (numericFilter) {
        _String errCode;

        long categoryCount = 1;

        // have multiple categories
        if (parameters.lLength > 2) {
          categoryCount =
              (long) ProcessNumericArgument((_String *)parameters(2), nil);
        }

        _String namesKey("FILTER_NAMES"), dataKey("FILTER_ARRAYS"),
            freqKey("FILTER_FREQS");

        _Matrix *sequenceNames =
            (_Matrix *)numericFilter->GetByKey(namesKey, MATRIX);
        _List seqNames;
        if (sequenceNames) {
          sequenceNames->FillInList(seqNames);
        }
        if (!sequenceNames || seqNames.lLength == 0) {
          errCode = _String("Expected a non-empty string matrix as the ") &
                    namesKey & " argument in call to CreateFilter";
        } else {
          _AssociativeList *dataList =
              (_AssociativeList *)numericFilter->GetByKey(dataKey,
                                                          ASSOCIATIVE_LIST);
          _Matrix *freqList =
              (_Matrix *)numericFilter->GetByKey(freqKey, MATRIX);

          if (dataList && freqList) {
            _List goodSeqs;
            long sitePatterns = freqList->GetVDim(), categDim = -1;

            if (freqList->GetHDim() != 1 || sitePatterns < 1 ||
                freqList->MatrixType() != 1) {
              errCode =
                  _String("Expected a non-empty numeric ROW matrix as the ") &
                  freqKey & " argument in call to CreateFilter";
            } else {
              for (long k = 0; k < seqNames.lLength; k = k + 1) {
                _Matrix *dataMx = (_Matrix *)dataList->GetByKey(k, MATRIX);
                if (dataMx && dataMx->MatrixType() == 1) {
                  if (categDim < 0) {
                    categDim = dataMx->GetVDim();
                    if (categDim < 1) {
                      break;
                    }
                  } else if (dataMx->GetVDim() != categDim) {
                    break;
                  }
                  if (dataMx->GetHDim() != sitePatterns * categoryCount) {
                    break;
                  }
                  goodSeqs << dataMx;
                  continue;
                }
                break;
              }

              if (goodSeqs.lLength == seqNames.lLength) {
                _DataSet *dummyDS = new _DataSet;
                dummyDS->SetNoSpecies(seqNames.lLength);
                dummyDS->GetNames().Duplicate(&seqNames);
                dummyDS->GetTheMap().Populate(sitePatterns, 0, 1);
                errCode = (*(_String *)parameters(0)) & "_internal_ds";
                dsID = FindDataSetName(errCode);
                if (dsID < 0) {
                  dataSetList << dummyDS;
                  DeleteObject(dummyDS);
                  dataSetNamesList &&&errCode;
                } else {
                  dataSetList.Replace(dsID, dummyDS, false);
                }

                errCode = (*(_String *)parameters(0));
                _DataSetFilterNumeric *dsn = new _DataSetFilterNumeric(
                    freqList, goodSeqs, dummyDS, categoryCount);
                checkPointer(dsn);
                dsID = FindDataSetFilterName(errCode);

                if (dsID < 0) {
                  dataSetFilterList << dsn;
                  DeleteObject(dsn);
                  dataSetFilterNamesList &&&errCode;
                } else {
                  dataSetFilterList.Replace(dsID, dsn, false);
                }
                return;

              } else {
                errCode =
                    _String("Site frequency patterns/numeric vectors did not "
                            "pass dimension checks in call to CreateFilter");
              }
            }
          }

        }
        if (errCode) {
          WarnError(errCode);
          return;
        }

      }
      _String errMsg =
          (((_String)("DataSet(Filter)/Associative Array ") & dataObjectID &
            _String(" has not been properly initialized")));
      WarnError(errMsg);
      return;
    }
    isFilter = true;
  }

  // build the formula from the 2nd parameter (unit size)
  char unit =
      ProcessNumericArgument((_String *)parameters(2), chain.nameSpacePrefix);

  // here's our unit
  _String dataFilterID(chain.AddNameSpaceToID(*(_String *)parameters(0))),
      hSpecs, vSpecs;

  long status = FindDataSetFilterName(dataFilterID);

  _DataSetFilter *thedf;

  if (status != -1) {
    thedf = (_DataSetFilter *)dataSetFilterList(status);
  } else {
    thedf = new _DataSetFilter();
    checkPointer(thedf);
    AddFilterToList(dataFilterID, thedf, false);
  }

  if (parameters.lLength > 3) {
    vSpecs = *(_String *)parameters(3);
  }
  if (parameters.lLength > 4) {
    hSpecs = *(_String *)parameters(4);
  } else {
    hSpecs = empty;
  }

  _DataSet *dataset;

  _SimpleList hL, vL;

  hL.RequestSpace(1024);
  vL.RequestSpace(1024);

  if (!isFilter) {
    dataset = (_DataSet *)dataSetList(dsID);
    dataset->ProcessPartition(hSpecs, hL, false);
    if (code != 6 && vSpecs.sLength == 0) {
      vSpecs = _String("0-") & _String(dataset->NoOfColumns() - 1);
    }
    dataset->ProcessPartition(vSpecs, vL, true);
  } else {
    _DataSetFilter *dataset1 = (_DataSetFilter *)dataSetFilterList(dsID);
    dataset1->GetData()->ProcessPartition(
        hSpecs, hL, false, &dataset1->theNodeMap, &dataset1->theOriginalOrder);

    if (code != 6 && vSpecs.sLength == 0) {
      vSpecs = _String("0-") & _String(dataset1->GetFullLengthSpecies() - 1);
    }

    dataset1->GetData()->ProcessPartition(
        vSpecs, vL, true, &dataset1->theOriginalOrder, &dataset1->theNodeMap);
    dataset = (_DataSet *)dataset1;
  }

  if (code != 6) {
    if (vL.lLength % unit) {
      vSpecs = (_String) "Unit size of " & unit &
               " doesn't divide the length of specified partition in call to ";
      if (code == 27) { // Permute
        vSpecs = vSpecs & "Permute";
      } else {
        vSpecs = vSpecs & "Bootstrap";
      }

      vSpecs = vSpecs & ". The partition has been trimmed at the end.";
      ReportWarning(vSpecs);
      for (status = vL.lLength % unit; status > 0; status--) {
        vL.Delete(vL.lLength - 1);
      }
    }
    if (code == 27) {
      vL.Permute(unit);
    } else {
      vL.PermuteWithReplacement(unit);
    }

  }

  thedf->SetFilter(dataset, unit, hL, vL, isFilter);

  if (parameters.lLength > 5) {
    hSpecs =
        GetStringFromFormula((_String *)parameters(5), chain.nameSpacePrefix);
    thedf->SetExclusions(&hSpecs);
  } else if ((code != 6) && isFilter) {
    _DataSetFilter *df1 = (_DataSetFilter *)dataSetFilterList(dsID);
    if (df1->theExclusions.lLength) {
      thedf->theExclusions.Duplicate(&df1->theExclusions);
      thedf->SetDimensions();
    }
  }

  thedf->SetDimensions();
  thedf->SetupConversion();

  SetDataFilterParameters(dataFilterID, thedf, true);
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase21(_ExecutionList &chain) {
  chain.currentCommand++;

  SetStatusLine(_hyStatusConditionProbsMatrix);
  _String errMsg,
      objectName = chain.AddNameSpaceToID(*(_String *)parameters(1)),
      resultID = chain.AddNameSpaceToID(*(_String *)parameters(0));

  long objectID = FindLikeFuncName(objectName, true);
  _PMathObj ob = nil;

  if (objectID >= 0) { // likelihood function
    _Matrix *partitionList = nil;
    if (parameters.lLength > 3) {
      _String secondArg = *(_String *)parameters(3);
      partitionList = (_Matrix *)ProcessAnArgumentByType(
          &secondArg, chain.nameSpacePrefix, MATRIX);
    }
    _SimpleList partsToDo;
    _LikelihoodFunction *lf = (_LikelihoodFunction *)likeFuncList(objectID);
    if (lf->ProcessPartitionList(partsToDo, partitionList,
                                 _hyStatusConditionProbsMatrix)) {
      char runMode = _hyphyLFConstructCategoryMatrixConditionals;
      if (parameters.lLength > 2) {
        if (((_String *)parameters(2))->Equal(&completeFlag)) {
          runMode = _hyphyLFConstructCategoryMatrixConditionals;
        } else if (((_String *)parameters(2))->Equal(&conditionalWeights)) {
          runMode = _hyphyLFConstructCategoryMatrixWeights;
        } else if (((_String *)parameters(2))->Equal(&siteProbabilities)) {
          runMode = _hyphyLFConstructCategoryMatrixSiteProbabilities;
        } else {
          runMode = _hyphyLFConstructCategoryMatrixClasses;
        }
      }
      ob = lf->ConstructCategoryMatrix(partsToDo, runMode, true, &resultID);
    }
  } else {
    _TheTree *testTree =
        (_TheTree *)FetchObjectFromVariableByType(&objectName, TREE);
    if (testTree) {
      long pid = 0;
      objectID = testTree->IsLinkedToALF(pid);
      if (objectID >= 0) {
        _LikelihoodFunction *anLF =
            (_LikelihoodFunction *)likeFuncList(objectID);
        _DataSetFilter *dsf =
            (_DataSetFilter *)dataSetFilterList(anLF->GetTheFilters()(pid));
        anLF->PrepareToCompute();
        anLF->Compute();
        objectID = dsf->NumberDistinctSites();

        _Matrix *condMx =
            new _Matrix(2 * objectID * (testTree->GetLeafCount() +
                                        testTree->GetINodeCount()) *
                            testTree->categoryCount,
                        testTree->GetCodeBase(), false, true);

        _List leafNames, inodeNames;
        testTree->DepthWiseT(true);

        while (testTree->currentNode) {
          _String *bs = new _String;
          testTree->GetNodeName(testTree->currentNode, *bs);
          if (testTree->IsCurrentNodeATip()) {
            leafNames << bs;
          } else {
            inodeNames << bs;
          }
          DeleteObject(bs);
          testTree->DepthWiseT(false);
        }

        leafNames << inodeNames;

        _Matrix *nodeNames = new _Matrix(leafNames);

        for (long siteC = 0; siteC < objectID; siteC++) {
          testTree->RecoverNodeSupportStates(dsf, siteC, siteC - 1, *condMx);
        }

        anLF->DoneComputing();
        _AssociativeList *retMe = new _AssociativeList;
        retMe->MStore("Nodes", nodeNames, false);
        retMe->MStore("Values", condMx, false);
        ob = retMe;
      }
    }
  }

  if (ob) {
    CheckReceptacleAndStore(&resultID, blConstructCM, true, ob, false);
  } else {
    WarnError(objectName & " must be either a likelihood function or a tree "
                           "variable tied to a likelihood function.");
  }

}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase53(_ExecutionList &chain) {
  chain.currentCommand++;

#ifdef __HYPHY_NO_SQLITE__
  _String errStr("SQLite commands can not be used in a HyPhy version built "
                 "with the __HYPHY_NO_SQLITE__ flag");
  WarnError(errStr);
#else

  _String arg1(*(_String *)parameters(0));

  char *errMsg = nil;
  _String errStr;

  if (arg1.Equal(&sqlOpen)) {
    _Variable *dbVar = CheckReceptacle((_String *)parameters(2), blDoSQL);

    if (dbVar) {
      _String arg2(*(_String *)parameters(1));
      arg2.ProcessFileName(true, true, (Ptr) chain.nameSpacePrefix);
      int errCode = SQLITE_OK;
      sqlite3 *aDB = nil;
#ifdef __HYPHYXCODE__
      errCode = sqlite3_open(DoMacToPOSIX(arg2).getStr(), &aDB);
#else
      errCode = sqlite3_open(arg2.sData, &aDB);
#endif
      if (errCode == SQLITE_OK) {
        errCode = sqlite3_exec(aDB, "SELECT COUNT(*) FROM sqlite_master",
                               _HYSQLCallBack, nil, nil);
      }
      if (errCode != SQLITE_OK) {
        WarnError(sqlite3_errmsg(aDB));
        sqlite3_close(aDB);
        return;
      } else {
        long f = sqlDatabases.Find(0);
        if (f < 0) {
          f = sqlDatabases.lLength;
          sqlDatabases << (long) aDB;
        } else {
          sqlDatabases.lData[f] = (long) aDB;
        }

        sqlite3_busy_timeout(aDB, 5000);

        dbVar->SetValue(new _Constant(f), false);
      }
    }
  } else {
    bool doClose = arg1.Equal(&sqlClose);

    long dbIdx = ProcessNumericArgument(
        doClose ? (_String *)parameters(2) : &arg1, chain.nameSpacePrefix);

    if (dbIdx < 0 || dbIdx >= sqlDatabases.lLength ||
        sqlDatabases.lData[dbIdx] == 0) {
      errStr = _String(dbIdx) & " is an invalid database index";
    } else {
      if (doClose) {
        sqlite3_close((sqlite3 *)sqlDatabases.lData[dbIdx]);
        sqlDatabases.lData[dbIdx] = 0;
      } else {
        _String arg3(ProcessLiteralArgument((_String *)parameters(2),
                                            chain.nameSpacePrefix));

        _ExecutionList sqlProcessor(
            arg3,
            chain.nameSpacePrefix ? (chain.nameSpacePrefix->GetName()) : nil);
        if (!terminateExecution) {
          _String arg2(ProcessLiteralArgument((_String *)parameters(1),
                                              chain.nameSpacePrefix));

          if (sqlite3_exec((sqlite3 *)sqlDatabases.lData[dbIdx], arg2.sData,
                           _HYSQLCallBack, (Ptr) & sqlProcessor, &errMsg) !=
              SQLITE_OK) {
            WarnError(sqlite3_errmsg((sqlite3 *)sqlDatabases.lData[dbIdx]));
            return;
          }
        }
      }
    }

  }

  if (errStr.sLength) {
    errStr = errStr & " in call to DoSQL";
    WarnError(errStr);
  }

#endif
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase54(_ExecutionList &chain) {
  chain.currentCommand++;

  SetStatusLine(_String("Constructing Topology ") & *(_String *)parameters(0));

  _String *treeSpec = ((_String *)parameters(1));
  treeSpec->ProcessParameter();
  _TreeTopology *tr = nil;

  if (treeSpec->sLength) {
    if (treeSpec->sData[0] != '(') {
      _Variable *testTree = FetchVar(LocateVarByName(
          AppendContainerName(*treeSpec, chain.nameSpacePrefix)));
      if (testTree && testTree->ObjectClass() == TREE) {
        tr = new _TreeTopology((_TheTree *)testTree);
      } else {
        _String flaData(*treeSpec);
        _Formula nameForm(flaData, chain.nameSpacePrefix);
        _PMathObj formRes = nameForm.Compute();
        if (formRes && formRes->ObjectClass() == STRING)
          tr = new _TreeTopology(AppendContainerName(*(_String *)parameters(0),
                                                     chain.nameSpacePrefix),
                                 *((_FString *)formRes)->theString, false);
      }
    } else
      tr = new _TreeTopology(
          AppendContainerName(*(_String *)parameters(0), chain.nameSpacePrefix),
          *(_String *)parameters(1), false);
  }

  if (!tr) {
    WarnError("Illegal right hand side in call to Topology id = ...; it must "
              "be a string, a Newick tree spec or a topology");
  }
}

//______________________________________________________________________________
bool _ElementaryCommand::ConstructAlignSequences(_String &source,
                                                 _ExecutionList &target)
    // syntax: AlignSequences (result, input string matrix,  options matrix)
    {
  _List pieces;
  ExtractConditions(source, blAlignSequences.sLength, pieces, ',');
  if (pieces.lLength != 3) {
    WarnError("Expected syntax: AlignSequences(result, input string matrix, "
              "options list);");
    return false;
  }

  _ElementaryCommand *as = new _ElementaryCommand(55);
  as->addAndClean(target, &pieces, 0);
  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::ConstructGetNeutralNull(_String &source,
                                                 _ExecutionList &target)
    // syntax: GetNeutralNull (result, likelihood function, syn sub count
    // matrix, non-syn sub count matrix, iterations per root state)
    {
  _List pieces;
  ExtractConditions(source, blGetNeutralNull.sLength, pieces, ',');
  if (pieces.lLength != 5) {
    WarnError(
        "Expected syntax: GetNeutralNull (result, likelihood function, syn sub "
        "count matrix, non-syn sub count matrix, iterations per root state);");
    return false;
  }

  _ElementaryCommand *gnn = new _ElementaryCommand(57);
  gnn->addAndClean(target, &pieces, 0);
  return true;
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase55(_ExecutionList &chain) {
  chain.currentCommand++;

  _String errStr;

  _Variable *storeResultIn = CheckReceptacle(
      &AppendContainerName(*(_String *)parameters(0), chain.nameSpacePrefix),
      blAlignSequences, true);

  if (storeResultIn) {
    _Matrix *inStrings = CheckMatrixArg(
        &AppendContainerName(*(_String *)parameters(1), chain.nameSpacePrefix),
        true);
    if (inStrings && (inStrings->GetHDim() == 1 || inStrings->GetVDim() == 1)) {
      _AssociativeList *mappingTable =
          CheckAssociativeListArg(&AppendContainerName(
              *(_String *)parameters(2), chain.nameSpacePrefix));
      if (mappingTable) {

        // check for required parameters
        _FString *charVector =
            (_FString *)mappingTable->GetByKey(seqAlignMap, STRING);

        long charCount = 0;

        _SimpleList ccount(256, -1, 0);

        if (charVector) {
          for (long cc = 0; cc < charVector->theString->sLength; cc++)
            if (ccount.lData[(unsigned char) charVector->theString->sData[cc]] >= 0) {
              // this is an error condition for
              // duplicate characters in the string
              charCount = 0;               
              break;
            } else {
              ccount.lData[(unsigned char) charVector->theString->sData[cc]] = cc;
              charCount++;
            }
        }

        if (charVector && charCount) {
          // now check that all characters
          bool doLocal = false, doAffine = false, doLinear = true,
               doCodon = false, doFullLocal = false;

          long codonCount = charCount * charCount * charCount;

          _PMathObj c = mappingTable->GetByKey(seqAlignCodonAlign, NUMBER);
          if (c) {
            doCodon = c->Compute()->Value() > 0.5;
          }

          _Matrix *scoreMatrix =
              (_Matrix *)mappingTable->GetByKey(seqAlignScore, MATRIX);
          if (scoreMatrix && scoreMatrix->GetHDim() ==
                                 (doCodon ? codonCount + 1 : charCount) &&
              scoreMatrix->GetVDim() == scoreMatrix->GetHDim()) {
            scoreMatrix = (_Matrix *)scoreMatrix->ComputeNumeric();
            scoreMatrix->CheckIfSparseEnough(true);

            char gapCharacter = '-';
            _FString *gapC =
                (_FString *)mappingTable->GetByKey(seqAlignGapChar, STRING);

            _Matrix *codon3x5 = nil, *codon3x4 = nil, *codon3x2 = nil,
                    *codon3x1 = nil;

            if (doCodon) {
              codon3x5 = (_Matrix *)mappingTable->GetByKey(seqAlignGapCodon3x5,
                                                           MATRIX);
              codon3x4 = (_Matrix *)mappingTable->GetByKey(seqAlignGapCodon3x4,
                                                           MATRIX);
              codon3x2 = (_Matrix *)mappingTable->GetByKey(seqAlignGapCodon3x2,
                                                           MATRIX);
              codon3x1 = (_Matrix *)mappingTable->GetByKey(seqAlignGapCodon3x1,
                                                           MATRIX);
              if (codon3x5 && codon3x4 && codon3x2 && codon3x1 &&
                  codon3x5->GetHDim() == codonCount + 1 &&
                  codon3x4->GetHDim() == codonCount + 1 &&
                  codon3x2->GetHDim() == codonCount + 1 &&
                  codon3x1->GetHDim() == codonCount + 1 &&
                  codon3x5->GetVDim() ==
                      charCount * charCount * charCount * 10 &&
                  codon3x4->GetVDim() ==
                      charCount * charCount * charCount * 4 &&
                  codon3x2->GetVDim() == charCount * charCount * 3 &&
                  codon3x1->GetVDim() == charCount * 3) {
                codon3x5 = (_Matrix *)codon3x5->ComputeNumeric();
                codon3x5->CheckIfSparseEnough(true);
                codon3x4 = (_Matrix *)codon3x4->ComputeNumeric();
                codon3x4->CheckIfSparseEnough(true);
                codon3x2 = (_Matrix *)codon3x2->ComputeNumeric();
                codon3x2->CheckIfSparseEnough(true);
                codon3x1 = (_Matrix *)codon3x1->ComputeNumeric();
                codon3x1->CheckIfSparseEnough(true);
              } else {
                errStr =
                    (seqAlignGapCodon3x5 & ", " & seqAlignGapCodon3x4 & ", " &
                     seqAlignGapCodon3x2 & ", or " & seqAlignGapCodon3x1 &
                     " matrices are missing or have incorrect dimensions");
              }

            }

            if (errStr.sLength == 0) {
              _String settingReport(128L, true);

              settingReport
                  << "Running sequence alignment with the following options:";

              if (gapC && gapC->theString->sLength == 1) {
                gapCharacter = gapC->theString->sData[0];
              }

              settingReport << "\n\tGap character:";
              settingReport << gapCharacter;

              _Parameter gapOpen = 15., gapOpen2 = 15., gapExtend = 1.,
                         gapExtend2 = 1., gapFrameshift = 50.;

              c = mappingTable->GetByKey(seqAlignGapOpen, NUMBER);
              if (c) {
                gapOpen = c->Compute()->Value();
              }

              settingReport << "\n\tGap open cost:";
              settingReport << _String(gapOpen);

              gapOpen2 = gapOpen;
              c = mappingTable->GetByKey(seqAlignGapOpen2, NUMBER);
              if (c) {
                gapOpen2 = c->Compute()->Value();
              }

              settingReport << "\n\tGap open cost 2:";
              settingReport << _String(gapOpen2);

              c = mappingTable->GetByKey(seqAlignGapExtend, NUMBER);
              if (c) {
                gapExtend = c->Compute()->Value();
              }

              settingReport << "\n\tGap extend cost:";
              settingReport << _String(gapExtend);

              gapExtend2 = gapExtend;
              c = mappingTable->GetByKey(seqAlignGapExtend2, NUMBER);
              if (c) {
                gapExtend2 = c->Compute()->Value();
              }

              settingReport << "\n\tGap extend cost 2:";
              settingReport << _String(gapExtend2);

              c = mappingTable->GetByKey(seqAlignFrameShift, NUMBER);
              if (c) {
                gapFrameshift = c->Compute()->Value();
              }

              settingReport << "\n\tCodon frameshift cost:";
              settingReport << _String(gapFrameshift);

              c = mappingTable->GetByKey(seqAlignGapLocal, NUMBER);
              if (c) {
                doLocal = c->Compute()->Value() > 0.5;
              }

              settingReport << "\n\tIgnore terminal gaps: ";
              settingReport << (doLocal ? "Yes" : "No");

              settingReport
                  << "\n\tUse codon alignment with frameshift routines: ";
              if (doCodon) {
                for (long i = 0; i < 256; i++)
                  if (ccount.lData[i] < 0) {
                    ccount.lData[i] = -codonCount - 1;
                  }
                settingReport << "Yes";
                doLinear = false;

                settingReport
                    << "\n\t Linear space routines  are not implemented";
              } else {
                settingReport << "No";
              }

              c = mappingTable->GetByKey(seqAlignDoLocal, NUMBER);
              if (c) {
                doFullLocal = c->Compute()->Value() > 0.5;
              }
              settingReport << "\n\tLocal alignment: ";
              settingReport << (doFullLocal ? "Yes" : "No");
              if (!doCodon && doFullLocal) {
                settingReport << "\n\t Local alignment is currently available "
                                 "for the codon aligner only.";
              }

              c = mappingTable->GetByKey(seqAlignGapAffine, NUMBER);
              if (c) {
                doAffine = c->Compute()->Value() > 0.5;
              }
              settingReport << "\n\tAffine gap costs: ";
              settingReport << (doAffine ? "Yes" : "No");

              c = mappingTable->GetByKey(seqAlignGapLinearSpace, NUMBER);
              if (c) {
                doLinear = c->Compute()->Value() > 0.5;
              }

              settingReport << "\n\tUse linear space routines: ";
              settingReport << (doLinear ? "Yes" : "No");

              settingReport.Finalize();
              ReportWarning(settingReport);

              long stringCount = inStrings->GetHDim() * inStrings->GetVDim();

              _AssociativeList *alignedStrings = new _AssociativeList;
              checkPointer(alignedStrings);

              for (long s1 = 0; s1 < stringCount; s1++) {
                _String *str1 = ((_FString *)inStrings->GetFormula(0, s1)
                                     ->Compute())->theString;
                if (!str1) {
                  errStr = _String("The ") & (s1 + 1) &
                           "-th argument is not a string";
                  break;
                }

                for (long s2 = s1 + 1; s2 < stringCount; s2++) {
                  _String *string2 = ((_FString *)inStrings->GetFormula(0, s2)
                                          ->Compute())->theString;
                  if (!string2) {
                    errStr = _String("The ") & (s2 + 1) &
                             "-th argument is not a string";
                    break;
                  }

                  _AssociativeList *pairwiseComp = new _AssociativeList;
                  checkPointer(pairwiseComp);

                  _Parameter score = 0.0;

                  if (doLinear == false) {
                    char *str1r = NULL, *str2r = NULL;
                    _List store;
                    score = AlignStrings(
                        str1->sData, string2->sData, str1r, str2r, ccount.lData,
                        scoreMatrix->theData, scoreMatrix->GetVDim(),
                        gapCharacter, gapOpen, gapExtend, gapOpen2, gapExtend2,
                        gapFrameshift, doLocal, doAffine, doCodon, charCount,
                        codon3x5->theData, codon3x4->theData, codon3x2->theData,
                        codon3x1->theData, doFullLocal);

                    if (str1r && str2r) {
                      _String *r_res =
                                  (_String *)checkPointer(new _String(str1r)),
                              *q_res =
                                  (_String *)checkPointer(new _String(str2r));
                      delete[] str1r;
                      delete[] str2r;
                      r_res->Finalize();
                      q_res->Finalize();
                      store.AppendNewInstance(r_res);
                      store.AppendNewInstance(q_res);
                    } else {
                      WarnError("Internal Error in AlignStrings");
                    }

                    store.bumpNInst();

                    if (store.lLength == 0) {
                      errStr = "Unspecified error in AlignStrings";
                      DeleteObject(pairwiseComp);
                      s1 = stringCount;
                      break;
                    } else {
                      pairwiseComp->MStore(
                          "1", new _FString((_String *)store(0)), false);
                      pairwiseComp->MStore(
                          "2", new _FString((_String *)store(1)), false);
                    }
                  } else {
                    _Matrix scoreM(string2->sLength + 1, 1, false, true),
                        scoreM2(string2->sLength + 1, 1, false, true),
                        gap1Matrix(string2->sLength + 1, 1, false, true),
                        gap2Matrix(string2->sLength + 1, 1, false, true),
                        gap1Matrix2(string2->sLength + 1, 1, false, true),
                        gap2Matrix2(string2->sLength + 1, 1, false, true),
                        *buffers[6];

                    char *alignmentRoute = new char[2 * (string2->sLength + 1)];

                    alignmentRoute[0] = alignmentRoute[string2->sLength + 1] =
                        0;
                    buffers[0] = &scoreM;
                    buffers[1] = &gap1Matrix;
                    buffers[2] = &gap2Matrix;
                    buffers[3] = &scoreM2;
                    buffers[4] = &gap1Matrix2;
                    buffers[5] = &gap2Matrix2;
                    _SimpleList ops(str1->sLength + 2, -2, 0);
                    ops.lData[str1->sLength + 1] = string2->sLength;
                    ops.lData[0] = -1;

                    score = LinearSpaceAlign(
                        str1, string2, ccount, scoreMatrix, gapOpen, gapExtend,
                        gapOpen2, gapExtend2, doLocal, doAffine, ops, score, 0,
                        str1->sLength, 0, string2->sLength, buffers, 0,
                        alignmentRoute);

                    delete[] alignmentRoute;

                    _String *result1 = new _String(str1->sLength + 1, true),
                            *result2 = new _String(string2->sLength + 1, true);

                    long last_column = ops.lData[ops.lLength - 1];

                    for (long position = str1->sLength - 1; position >= 0;
                         position--) {
                      long current_column = ops.lData[position + 1];

                      if (current_column < 0) {
                        //|| (current_column == -3 &&last_column == string2->sLength)
                        if (current_column == -2 ) {
                          current_column = last_column;
                        } else if (current_column == -3) {
                          // find the next matched char or a -1
                          long p = position, s2p;
                          while ((ops.lData[p + 1]) < -1) {
                            p--;
                          }

                          s2p = ops.lData[p + 1];
                          //if (last_column == string2->sLength)
                          //  last_column = string2->sLength-1;

                          //if (s2p < 0)
                          //  s2p = 0;

                          for (long j = last_column - 1; j > s2p;) {
                            (*result1) << gapCharacter;
                            (*result2) << string2->sData[j--];
                          }

                          last_column = s2p + 1;

                          for (; position > p; position--) {
                            (*result2) << gapCharacter;
                            (*result1) << str1->sData[position];
                          }
                          position++;
                          continue;
                        } else {
                          for (last_column--; last_column >= 0; last_column--) {
                            (*result1) << gapCharacter;
                            (*result2) << string2->sData[last_column];
                          }
                          while (position >= 0) {
                            (*result1) << str1->sData[position--];
                            (*result2) << gapCharacter;
                          }
                          break;
                        }
                      }

                      if (current_column == last_column) { 
                        // insert in sequence 2
                        (*result1) << str1->sData[position];
                        (*result2) << gapCharacter;
                      } else {

                        last_column--;
                        for (; last_column > current_column;
                             last_column--) { // insert in column 1
                          (*result2) << string2->sData[last_column];
                          (*result1) << gapCharacter;
                        }
                        (*result1) << str1->sData[position];
                        (*result2) << string2->sData[current_column];
                      }
                      //printf ("%s\n%s\n", result1->sData, result2->sData);
                    }

                    for (last_column--; last_column >= 0; last_column--) {
                      (*result1) << gapCharacter;
                      (*result2) << string2->sData[last_column];
                    }

                    result1->Finalize();
                    result1->Flip();
                    result2->Finalize();
                    result2->Flip();
                    pairwiseComp->MStore("1", new _FString(result1), false);
                    pairwiseComp->MStore("2", new _FString(result2), false);
                  }
                  /*
                                    long gap1c = 0,
                                         gap2c = 0;

                                     _Parameter scoreCheck =
                  0.;

                                     for (long sp = 0;
                  sp<result1->sLength; sp++)
                                     {
                                         char cs1 =
                  result1->sData[sp],
                                              cs2 =
                  result2->sData[sp];

                                         if (cs1 ==
                  gapCharacter)
                                         {
                                             if (gap1c &&
                  doAffine)
                                                 scoreCheck -=
                  gapExtend;
                                             else
                                                 scoreCheck -=
                  gapOpen;
                                             gap2c = 0;
                                             gap1c++;
                                         }
                                         else
                                         if (cs2 ==
                  gapCharacter)
                                         {
                                             if (gap2c &&
                  doAffine)
                                                 scoreCheck -=
                  gapExtend2;
                                             else
                                                 scoreCheck -=
                  gapOpen2;
                                             gap1c = 0;
                                             gap2c++;
                                         }
                                         else
                                         {
                                             gap1c = 0;
                                             gap2c = 0;
                                             long code1 =
                  ccount.lData[cs1],
                                                  code2 =
                  ccount.lData[cs2];

                                             if (code1 >=0 &&
                  code2 >=0 )
                                                 scoreCheck +=
                  (*scoreMatrix)(code1,code2);
                                         }
                                     }
                                     if (doLocal)
                                     {
                                        for (long k = 0;
                  result1->sData[k] == gapCharacter; k++)
                                            if (doAffine)
                                                scoreCheck +=
                  k?gapExtend:gapOpen;
                                            else
                                                scoreCheck +=
                  gapOpen;
                                         for (long k = 0;
                  result2->sData[k] == gapCharacter; k++)
                                             if (doAffine)
                                                 scoreCheck +=
                  k?gapExtend2:gapOpen2;
                                             else
                                                 scoreCheck +=
                  gapOpen2;
                                         for (long k =
                  result1->sLength-1; result1->sData[k] == gapCharacter; k--)
                                             if (doAffine)
                                                 scoreCheck +=
                  k==result1->sLength-1?gapOpen:gapExtend;
                                             else
                                                 scoreCheck +=
                  gapOpen;
                                         for (long k =
                  result2->sLength-1; result2->sData[k] == gapCharacter; k--)
                                             if (doAffine)
                                                 scoreCheck +=
                  k==result2->sLength-1?gapOpen2:gapExtend2;
                                             else
                                                 scoreCheck +=
                  gapOpen2;
                                     }*/

                  pairwiseComp->MStore("0", new _Constant(score), false);
                  /*pairwiseComp->MStore ("3", new _Constant (score2), false);
                  pairwiseComp->MStore ("4", new _FString(result1), false);
                  pairwiseComp->MStore ("5", new _FString(result2), false);
                  pairwiseComp->MStore ("6", new
                  _FString((_String*)ops.toStr()), false);
                  pairwiseComp->MStore ("7", new _Constant (scoreCheck),
                  false);*/
                  alignedStrings->MStore(_String(s1), pairwiseComp, false);
                }
              }

              storeResultIn->SetValue(alignedStrings, false);
            }

          } else {
            errStr = seqAlignScore &
                     " is a required option, which must be a square matrix "
                     "with dimension matching the size of " & seqAlignMap;
          }
        } else {
          errStr =
              seqAlignMap & " is a required option, which must be a non-empty "
                            "string without repeating characters ";
        }
      } else {
        errStr =
            *(_String *)parameters(2) &
            " was expected to be an associative array of alignment options";
      }
    } else {
      errStr =
          *(_String *)parameters(1) & " was expected to be a vector of strings";
    }
  }

  if (errStr.sLength) {
    errStr = errStr & " in call to " & blAlignSequences;
    WarnError(errStr);
  }
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase26(_ExecutionList &chain) {

  chain.currentCommand++;
  // we have to build a list of _CalcNodes to deal with
  // all of the trees/nodes in ReplicateConstraint must be of the same topology
  // the constraint will be processed by trying all of the subnodes of the given
  // node
  // and within each - trying all of the variables to see if the constraint is
  // matched
  // exactly the same operation will be repeated on each of the other parameters

  _String *replicateSource, thisS,
      prStr =
          GetStringFromFormula((_String *)parameters(0), chain.nameSpacePrefix);

  replicateSource = &prStr;

  _List parts, theConstraints;

  _SimpleList thisIndex, thisArgs;

  long ind1 = replicateSource->Find("this"), ind2, ind3, ind4;

  if (ind1 < 0) {
    WarnError(*(_String *)parameters(0) &
              " has no 'this' references in call to ReplicateConstraint!");
    return;
  }

  _SimpleList thisHits(parameters.lLength - 1, 0, 0);

  while (ind1 >= 0) { // references to 'this' still exist
    ind2 = ind1 + 4;  // look forward to the number of 'this'
    while ('0' <= replicateSource->sData[ind2] &&
           replicateSource->sData[ind2] <= '9') {
      ind2++;
    }

    ind3 = replicateSource->Cut(ind1 + 4, ind2 - 1).toNum();
    ind2 = replicateSource->FindEndOfIdent(ind1, -1, '?');
    // now ind1-ind2 contains a reference with this...
    _String newS(*replicateSource, ind1, ind2);
    thisS = _String("this") & _String(ind3);
    if ((ind4 = ((_String *)parameters(ind3))->Find('.')) >= 0) { 
      // branch argument
      newS = newS.Replace(thisS, ((_String *)parameters(ind3))->Cut(0, ind4 - 1), true);
    } else { 
      // non-branch argument
      newS = newS.Replace(thisS, *((_String *)parameters(ind3)), true);
    }
    parts &&&newS;
    ind3--;
    thisIndex << ind3; // sequence of references to this

    if (ind3 < 0 || ind3 >= thisHits.lLength) {
      WarnError(_String("Invalid reference to ") & thisS &
                " in the constraint specification");
      return;
    }

    thisHits.lData[ind3] = 1;
    if (ind2 >= replicateSource->sLength - 1) {
      break;
    }
    ind1 = replicateSource->Find("this", ind2 + 1, -1);
    if (ind1 == -1) {
      newS = replicateSource->Cut(ind2 + 1, -1);
    } else {
      newS = replicateSource->Cut(ind2 + 1, ind1 - 1);
    }
    parts &&&newS;
    thisIndex << -1;
  }

  // now that the string is conveniently partritioned into blocks
  // we will check the arguments and store references
  for (ind1 = 1; ind1 < parameters.lLength; ind1++) {
    if (thisHits.lData[ind1 - 1] == 0) {
      WarnError(_String("Unused ") & ind1 & "-th reference variable: " &
                *(_String *)parameters(ind1));
      return;
    }

    ind2 = LocateVarByName(*(_String *)parameters(ind1));
    if (ind2 < 0) {
      _String newS = *(_String *)parameters(ind1) &
                     " is undefined in call to ReplicateConstraint.";
      acknError(newS);
      return;
    }

    _Variable *thisNode = FetchVar(ind2);
    if (thisNode->ObjectClass() == TREE_NODE) {
      thisArgs << (long)((_CalcNode *)thisNode)->LocateMeInTree();
    } else if (thisNode->ObjectClass() == TREE) {
      thisArgs << (long) & ((_TheTree *)thisNode)->GetRoot();
    } else {
      WarnError(
          *(_String *)parameters(ind1) &
          " is neither a tree nor a tree node in call to ReplicateConstraint.");
      return;
    }
  }

  // now with this list ready we can recurse down the tree and produce the
  // contsraints
  if (RecurseDownTheTree(thisArgs, parameters, theConstraints, parts,
                         thisIndex)) {
    if (theConstraints.lLength) {
      ReportWarning(_String(
          "\nReplicateConstraint generated the following contsraints:"));
      _Parameter doDeferSet;
      checkParameter(deferConstrainAssignment, doDeferSet, 0.0);
      bool applyNow = CheckEqual(doDeferSet, 0.0);
      _String *constraintAccumulator =
          (_String *)checkPointer(new _String(128L, true));

      if (applyNow) {
        deferSetFormula = new _SimpleList;
        checkPointer(deferSetFormula);
      }

      for (ind1 = 0; ind1 < theConstraints.lLength; ind1++) {
        replicateSource = (_String *)(theConstraints(ind1)->toStr());
        if (applyNow) {
          _Formula rhs, lhs;
          _FormulaParsingContext fpc(nil, chain.nameSpacePrefix);
          ind2 = Parse(&rhs, *replicateSource, fpc, &lhs);
          //ExecuteFormula(&rhs, &lhs, ind2, fpc.assignmentRefID(),
          //               chain.nameSpacePrefix, fpc.assignmentRefType());
        }

        (*constraintAccumulator) << replicateSource;
        (*constraintAccumulator) << ';';
        (*constraintAccumulator) << '\n';
        //ReportWarning (*replicateSource);
        DeleteObject(replicateSource);
      }
      constraintAccumulator->Finalize();
      ReportWarning(*constraintAccumulator);
      CheckReceptacleAndStore(&lastSetOfConstraints, "ReplicateConstraint",
                              false, new _FString(constraintAccumulator),
                              false);
      if (applyNow) {
        FinishDeferredSF();
      }
    }
  }
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase57(_ExecutionList &chain) {

  chain.currentCommand++;
  _String errStr;

  _Variable *storeResultIn =
                CheckReceptacle(&AppendContainerName(*(_String *)parameters(0),
                                                     chain.nameSpacePrefix),
                                blGetNeutralNull, true),
            *sv = FetchVar(LocateVarByName(AppendContainerName(
                *(_String *)parameters(2), chain.nameSpacePrefix))),
            *nsv = FetchVar(LocateVarByName(AppendContainerName(
                *(_String *)parameters(3), chain.nameSpacePrefix)));

  _Parameter itCountV =
      ProcessNumericArgument((_String *)parameters(4), chain.nameSpacePrefix);

  _String *lfName = (_String *)parameters(1);

  long f =
      FindLikeFuncName(AppendContainerName(*lfName, chain.nameSpacePrefix));

  if (f >= 0) {
    if (sv && sv->ObjectClass() == MATRIX) {
      if (nsv && nsv->ObjectClass() == MATRIX) {
        _Matrix *sMatrix =
            (_Matrix *)((_Matrix *)sv->Compute())->ComputeNumeric();
        _Matrix *nsMatrix =
            (_Matrix *)((_Matrix *)nsv->Compute())->ComputeNumeric();

        sMatrix->CheckIfSparseEnough(true);
        nsMatrix->CheckIfSparseEnough(true);

        if (sMatrix->GetHDim() == sMatrix->GetVDim() &&
            nsMatrix->GetHDim() == nsMatrix->GetVDim() &&
            sMatrix->GetHDim() == nsMatrix->GetVDim()) {
          _LikelihoodFunction *theLF = (_LikelihoodFunction *)likeFuncList(f);

          if (((_DataSetFilter *)dataSetFilterList(theLF->GetTheFilters()(0)))
                  ->GetDimension(true) == sMatrix->GetHDim()) {
            long itCount = itCountV;
            if (itCount > 0) {
              _AssociativeList *res = theLF->SimulateCodonNeutral(
                  (_Matrix *)sMatrix, (_Matrix *)nsMatrix, itCount);
              storeResultIn->SetValue(res, false);
            } else {
              errStr = "Invalid iterations per character state";
            }
          } else {
            errStr = "Incompatible data and cost matrices";
          }
        } else {
          errStr = "Incompatible syn and non-syn cost matrix dimensions";
        }
      } else {
        errStr = "Invalid non-syn cost matrix argument";
      }
    } else {
      errStr = "Invalid syn cost matrix argument";
    }

  } else {
    errStr =
        _String("Likelihood function ") & *lfName & " has not been defined";
  }

  if (errStr.sLength) {
    errStr = errStr & " in call to " & blGetNeutralNull;
    WarnError(errStr);
  }
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase58(_ExecutionList &chain) {

  chain.currentCommand++;

  _String errStr;
  _String *profileCode = (_String *)parameters(0);

  if (*profileCode == _String("START")) {
    if (chain.profileCounter) {
      DeleteObject(chain.profileCounter);
    }
    checkPointer(chain.profileCounter =
                     new _Matrix(chain.lLength, 2, false, true));
    chain.doProfile = 1;
  } else if (*profileCode == _String("PAUSE")) {
    chain.doProfile = 2;
  } else if (*profileCode == _String("RESUME")) {
    chain.doProfile = 1;
  } else {
    _Variable *outVar = CheckReceptacle(
        &AppendContainerName(*profileCode, chain.nameSpacePrefix), blHBLProfile,
        true);
    if (outVar) {
      if (chain.profileCounter) {
        _AssociativeList *profileDump = new _AssociativeList;
        checkPointer(profileDump);

        _SimpleList instructions;
        _List descriptions;

        for (long k = 1; k < 2 *chain.lLength; k += 2) {
          if (chain.profileCounter->theData[k] > 0.0) {
            instructions << k / 2;
            _String *desc =
                (_String *)((_ElementaryCommand *)chain(k / 2))->toStr();
            descriptions << desc;
            DeleteObject(desc);
          }
        }

        _Matrix *execProfile =
                    new _Matrix(instructions.lLength, 2, false, true),
                *instCounter = new _Matrix(instructions),
                *descList = new _Matrix(descriptions);

        checkPointer(execProfile);
        checkPointer(instCounter);
        checkPointer(descList);

        long k2 = 0;
        for (long m = 1; m < 2 *chain.lLength; m += 2) {
          if (chain.profileCounter->theData[m] > 0.0) {
            execProfile->theData[k2++] = chain.profileCounter->theData[m];
            execProfile->theData[k2++] = chain.profileCounter->theData[m - 1];
          }
        }

        _FString aKey;
        *aKey.theString = "INSTRUCTION INDEX";
        profileDump->MStore(&aKey, instCounter, false);
        *aKey.theString = "INSTRUCTION";
        profileDump->MStore(&aKey, descList, false);
        *aKey.theString = "STATS";
        profileDump->MStore(&aKey, execProfile, false);
        outVar->SetValue(profileDump, false);
        chain.doProfile = 0;
        DeleteObject(chain.profileCounter);
        chain.profileCounter = nil;
      } else {
        errStr = "Profiler dump invoked before #profile START; ";
      }
    }
  }

}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase61(_ExecutionList &chain) {

  chain.currentCommand++;
  _PMathObj avl1 = FetchObjectFromVariableByType(
                &AppendContainerName(*(_String *)parameters(1),
                                     chain.nameSpacePrefix),
                ASSOCIATIVE_LIST),
            avl2 = FetchObjectFromVariableByType(
                &AppendContainerName(*(_String *)parameters(2),
                                     chain.nameSpacePrefix),
                ASSOCIATIVE_LIST),
            start = parameters.lLength > 3
                        ? FetchObjectFromVariableByType(
                              &AppendContainerName(*(_String *)parameters(3),
                                                   chain.nameSpacePrefix),
                              NUMBER)
                        : nil;

  if (!(avl1 && avl2)) {
    WarnError(
        _String("Both arguments (") & *(_String *)parameters(1) & " and " &
        *(_String *)parameters(2) &
        " in a call to SCFG = ... must be evaluate to associative arrays");
  } else {
    Scfg *scfg = new Scfg((_AssociativeList *)avl1, (_AssociativeList *)avl2,
                          start ? start->Value() : 0);
    _String scfgName =
        AppendContainerName(*(_String *)parameters(0), chain.nameSpacePrefix);
    long f = FindSCFGName(scfgName);

    if (f == -1) {
      for (f = 0; f < scfgNamesList.lLength; f++)
        if (((_String *)scfgNamesList(f))->sLength == 0) {
          break;
        }

      if (f == scfgNamesList.lLength) {
        scfgList << scfg;
        scfgNamesList &&(&scfgName);
        DeleteObject(scfg);
      } else {
        scfgNamesList.Replace(f, &scfgName, true);
        scfgList.lData[f] = (long) scfg;
      }
    } else {
      scfgNamesList.Replace(f, &scfgName, true);
      scfgList.Replace(f, scfg, false);
    }
  }
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase63(_ExecutionList &chain) {

  chain.currentCommand++;

  /*_PMathObj         avl1    = FetchObjectFromVariableByType
  ((_String*)parameters(1),ASSOCIATIVE_LIST),
                        avl2    = FetchObjectFromVariableByType
  ((_String*)parameters(2),ASSOCIATIVE_LIST),
                        start   =
  parameters.lLength>3?FetchObjectFromVariableByType
  ((_String*)parameters(3),NUMBER):nil;

    if (! (avl1 && avl2))
    {
        _String errMsg = _String ("Both arguments (") &
  *(_String*)parameters(1) & " and " & *(_String*)parameters(2) & " in a call to
  SCFG = ... must be evaluate to associative arrays";
        WarnError (errMsg);
    }
    else
    {
        Scfg    * scfg   = new Scfg
  ((_AssociativeList*)avl1,(_AssociativeList*)avl2,start?start->Value():0);
        _String * str    = (_String*)parameters(0);
        long    f        = FindSCFGName (*str);

        if (f==-1)
        {
            for (f=0; f<scfgNamesList.lLength; f++)
                if (((_String*)scfgNamesList(f))->sLength==0)
                    break;

            if (f==scfgNamesList.lLength)
            {
                scfgList << scfg;
                scfgNamesList&&(str);
                DeleteObject (scfg);
            }
            else
            {
                scfgNamesList.Replace(f,str,true);
                scfgList.lData[f] = (long)scfg;
            }
        }
        else
        {
            scfgNamesList.Replace(f,str,true);
            scfgList.Replace(f,scfg,false);
        }
    }   */
}

//______________________________________________________________________________
void _ElementaryCommand::ExecuteCase64(_ExecutionList &chain) {
  ReportWarning(_String("ExecuteCase64()"));
  chain.currentCommand++;

  _PMathObj avl1 = FetchObjectFromVariableByType(
      &AppendContainerName(*(_String *)parameters(1), chain.nameSpacePrefix),
      ASSOCIATIVE_LIST);

  if (!(avl1)) {
    WarnError(_String("Argument (") & *(_String *)parameters(1) &
              " in call to BGM = ... must evaluate to associative array");
  } else {
    _BayesianGraphicalModel *bgm =
        new _BayesianGraphicalModel((_AssociativeList *)avl1);

    _String bgmName =
        AppendContainerName(*(_String *)parameters(0), chain.nameSpacePrefix);
    long bgmIndex = FindBgmName(bgmName);

    if (bgmIndex == -1) { // not found
      for (bgmIndex = 0; bgmIndex < bgmNamesList.lLength; bgmIndex++) {
        // locate empty strings in name list
        if (((_String *)bgmNamesList(bgmIndex))->sLength == 0) {
          break;
        }
      }

      if (bgmIndex == bgmNamesList.lLength) {
        // reached end of list without finding empty string, append new string
        bgmList.AppendNewInstance(bgm);
        bgmNamesList &&(&bgmName);
      } else {
        // replace empty string in list
        bgmNamesList.Replace(bgmIndex, &bgmName, true);
        bgmList.Replace(bgmIndex, bgm, false);
      }
    } else { 
      // 20070626: SLKP edit to deal with already existing BGMs
      bgmNamesList.Replace(bgmIndex, &bgmName, true);
      bgmList.Replace(bgmIndex, bgm, false);
    }

    ReportWarning(_String("Created BGM ") & bgmName & " at index " & bgmIndex);
  }
}

//______________________________________________________________________________
// syntax: SCFG ident = (Rules1, Rules2 <,start>)
bool _ElementaryCommand::ConstructSCFG(_String &source, _ExecutionList &target) {

  long mark1 = source.FirstSpaceIndex(0, -1, 1),
       mark2 = source.Find('=', mark1, -1);

  _String scfgID(source, mark1 + 1, mark2 - 1);

  if (mark1 == -1 || mark2 == -1 || mark1 + 1 > mark2 - 1 ||
      !scfgID.IsValidIdentifier(true)) {
    WarnError("SCFG declaration missing a valid identifier");
    return false;
  }

  _List pieces;

  mark1 = source.Find('(', mark2, -1);
  if (mark1 >= 0) {
    ExtractConditions(source, mark1 + 1, pieces, ',');
  }

  if (pieces.lLength != 2 && pieces.lLength != 3) {
    WarnError("Expected: SCFG ident = (Rules1, Rules2 <,start>)");
    return false;
  }

  _ElementaryCommand *scfg = new _ElementaryCommand(61);

  scfg->parameters &&(&scfgID);
  scfg->addAndClean(target, &pieces, 0);
  return true;
}

//______________________________________________________________________________
// syntax: NeuralNet ident = (InMatrix,OutMatrix,HiddenNodes)
bool _ElementaryCommand::ConstructNN(_String &source, _ExecutionList &target) {

  long mark1 = source.FirstSpaceIndex(0, -1, 1),
       mark2 = source.Find('=', mark1, -1);

  _String nnID(source, mark1 + 1, mark2 - 1);

  if (mark1 == -1 || mark2 == -1 || mark1 + 1 > mark2 - 1 ||
      !nnID.IsValidIdentifier(true)) {
    WarnError("NeutalNet declaration missing a valid identifier");
    return false;
  }

  _List pieces;

  mark1 = source.Find('(', mark2, -1);
  if (mark1 >= 0) {
    ExtractConditions(source, mark1 + 1, pieces, ',');
  }

  if (pieces.lLength != 3) {
    WarnError("NeuralNet ident = (InMatrix,OutMatrix,HiddenNodes)");
    return false;
  }

  _ElementaryCommand *nn = new _ElementaryCommand(63);
  nn->parameters &&(&nnID);
  nn->addAndClean(target, &pieces, 0);
  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::ConstructBGM(_String &source, _ExecutionList &target) {

  // syntax: BGM ident = (<nodes>)
  ReportWarning(_String("ConstructBGM()"));
  // locate ident in HBL string
  long mark1 = source.FirstSpaceIndex(0, -1, 1),
       mark2 = source.Find('=', mark1, -1);

  // assign ident to _String variable
  _String bgmID(source, mark1 + 1, mark2 - 1);

  if (mark1 == -1 || mark2 == -1 || mark1 + 1 > mark2 - 1 ||
      !bgmID.IsValidIdentifier(true)) {
    WarnError("BGM declaration missing a valid identifier");
    return false;
  }

  // extract arguments from remainder of HBL string
  _List pieces;

  mark1 = source.Find('(', mark2, -1);
  if (mark1 >= 0) {
    ExtractConditions(source, mark1 + 1, pieces, ',');
  }

  if (pieces.lLength != 1) {
    WarnError("Expected: BGM ident = (<nodes>)");
    return false;
  }

  _ElementaryCommand *bgm = new _ElementaryCommand(64);
  bgm->parameters &&(&bgmID);
  bgm->addAndClean(target, &pieces, 0);

  return true;
}

//______________________________________________________________________________
bool
_ElementaryCommand::HandleHarvestFrequencies(_ExecutionList &currentProgram) {
  currentProgram.currentCommand++;

  _String freqStorageID = *(_String *)parameters(0),
          dataID = currentProgram.AddNameSpaceToID(*(_String *)parameters(1)),
          errMsg;

  _Variable *theReceptacle = CheckReceptacleCommandID(
      &AppendContainerName(freqStorageID, currentProgram.nameSpacePrefix),
      HY_HBL_COMMAND_HARVEST_FREQUENCIES, true, false, &currentProgram);

  if (!theReceptacle) {
    return false;
  }

  SetStatusLine("Gathering Frequencies");

  long objectType = HY_BL_DATASET | HY_BL_DATASET_FILTER;
  BaseRef sourceObject =
      _HYRetrieveBLObjectByName(dataID, objectType, nil, false);
      
  long unit = ProcessNumericArgument((_String *)parameters(2),currentProgram.nameSpacePrefix),
       posspec = ProcessNumericArgument((_String *)parameters(4), currentProgram.nameSpacePrefix),
       atom = ProcessNumericArgument((_String *)parameters(3), currentProgram.nameSpacePrefix);

  _Matrix *receptacle = nil;

  _Parameter cghf = 1.0;
  checkParameter(hfCountGap, cghf, 1.0, currentProgram.nameSpacePrefix);

  if (objectType == HY_BL_DATASET) { // harvest directly from a DataSet

    _String vSpecs, hSpecs;

    if (parameters.lLength > 5) {
      vSpecs = *(_String *)parameters(5);
    }

    if (parameters.lLength > 6) {
      hSpecs = *(_String *)parameters(6);
    }

    _DataSet *dataset = (_DataSet *)sourceObject;
    _SimpleList hL, vL;
    dataset->ProcessPartition(hSpecs, hL, false);
    dataset->ProcessPartition(vSpecs, vL, true);

    receptacle = dataset->HarvestFrequencies(unit, atom, posspec, hL, vL, cghf > 0.5);

  } else { // harvest from a DataSetFilter
    if (objectType == HY_BL_DATASET_FILTER) {
      receptacle = ((_DataSetFilter *)sourceObject)->HarvestFrequencies(unit, atom, posspec, cghf > 0.5);
    } else {
      errMsg = _String("'") & dataID & "' is neither a DataSet nor a DataSetFilter";
    }
  }

  SetStatusLine(empty);

  if (errMsg.sLength || receptacle == nil) {
    DeleteObject(receptacle);
    currentProgram.ReportAnExecutionError(errMsg);
    theReceptacle->SetValue(new _MathObject, false);
    return false;
  }

  theReceptacle->SetValue(receptacle, false);
  return true;

  //CheckReceptacleCommandIDAndStore
  //(&freqStorageID,HY_HBL_COMMAND_HARVEST_FREQUENCIES,true, receptacle, false);
}

//______________________________________________________________________________
bool _ElementaryCommand::HandleOptimizeCovarianceMatrix(
    _ExecutionList &currentProgram, bool doOptimize) {

  currentProgram.currentCommand++;

  // construct the receptacle matrix
  _String lfResName(currentProgram.AddNameSpaceToID(*(_String *)parameters(0))),
      lfNameID(currentProgram.AddNameSpaceToID(*(_String *)parameters(1)));

  _Variable *result = CheckReceptacleCommandID(
      &lfResName,
      doOptimize ? HY_HBL_COMMAND_OPTIMIZE : HY_HBL_COMMAND_COVARIANCE_MATRIX,
      true);

  // Handle string variables passed as likefunc IDs?
  _String temp = ProcessLiteralArgument(&lfNameID, currentProgram.nameSpacePrefix);
  if (temp.sLength) {
    lfNameID = temp;
  }

  long objectType = HY_BL_LIKELIHOOD_FUNCTION | HY_BL_SCFG | HY_BL_BGM;

  _LikelihoodFunction *lkf = (_LikelihoodFunction *)_HYRetrieveBLObjectByName(lfNameID, objectType, nil, doOptimize == false);

  if (lkf == nil) { // will only happen if the object is a custom function
    lkf = (_LikelihoodFunction *)checkPointer(new _CustomFunction(&lfNameID));
  }

  if (!doOptimize) {

    // COVARIANCE_MATRIX
    SetStatusLine(_String("Finding the cov. matrix/profile CI for ") & lfNameID);
    _String cpl = currentProgram.AddNameSpaceToID(covarianceParameterList);
    _Variable *restrictVariable = FetchVar(LocateVarByName(cpl));
    _SimpleList *restrictor = nil;

    if (objectType == HY_BL_LIKELIHOOD_FUNCTION || objectType == HY_BL_SCFG) {
      // not a BGM
      if (restrictVariable) { // only consider some variables
        _SimpleList variableIDs;

        // a list of variables stored as keys in an associative array
        if (restrictVariable->ObjectClass() == ASSOCIATIVE_LIST) {
          checkPointer(restrictor = new _SimpleList);
          _List *restrictedVariables = ((_AssociativeList *)restrictVariable->GetValue())->GetKeys();
          for (unsigned long iid = 0; iid < restrictedVariables->lLength; iid++) {
            _String varID = currentProgram.AddNameSpaceToID(*(_String *)(*restrictedVariables)(iid));
            variableIDs << LocateVarByName(varID);
          }

        } else if (restrictVariable->ObjectClass() == STRING) {
          // a single variable stored in a string
          _String varID = currentProgram.AddNameSpaceToID(
              *((_FString *)restrictVariable->Compute())->theString);
          variableIDs << LocateVarByName(varID);
        }
        if (variableIDs.lLength > 0) {
          checkPointer(restrictor = new _SimpleList());
          for (unsigned long var_index = 0; var_index < variableIDs.lLength; var_index++) {
            long vID = lkf->GetIndependentVars().Find(variableIDs.lData[var_index]);
            if (vID >= 0) {
              (*restrictor) << vID;
            }
          }
          if (restrictor->lLength == 0) {
            DeleteObject(restrictor);
            restrictor = nil;
          }
        }
      }
      result->SetValue((_Matrix *)lkf->CovarianceMatrix(restrictor), false);
      DeleteObject(restrictor);
    } else {
      // BGM
      _Matrix *optRes = (_Matrix *)lkf->CovarianceMatrix(nil);
      if (optRes) {
        result->SetValue(optRes, false);
      }
    }
  } else {
    // OPTIMIZE
    if (objectType != HY_BL_NOT_DEFINED) {
      SetStatusLine(_String("Optimizing ") & _HYHBLTypeToText(objectType) & ' ' & lfNameID);
    } else {
      SetStatusLine(_String("Optimizing user function ") & lfNameID);
    }
    result->SetValue(lkf->Optimize(), false);
  }

  if (objectType == HY_BL_NOT_DEFINED) {
    DeleteObject(lkf); // delete the custom function object
  }

  SetStatusLine("Finished with the optimization");

  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::HandleComputeLFFunction(_ExecutionList &currentProgram) {

  currentProgram.currentCommand++;

  _String *arg1 = (_String *)parameters(0), *arg2 = (_String *)parameters(1),
          name2Find = AppendContainerName(*arg1, currentProgram.nameSpacePrefix);

  // bool isSCFG  = false;

  long objectType = HY_BL_LIKELIHOOD_FUNCTION | HY_BL_SCFG | HY_BL_BGM;
  _LikelihoodFunction *lf = (_LikelihoodFunction *)_HYRetrieveBLObjectByName(name2Find, objectType, nil, true, true);

  if (*arg2 == lfStartCompute) {
    lf->PrepareToCompute(true);
  } else if (*arg2 == lfDoneCompute) {
    lf->DoneComputing(true);
  } else {
    if (!lf->HasBeenSetup()) {
      WarnError(_String("Please call LFCompute (lf_id, ") & lfStartCompute & ") before evaluating the likelihood function");
      return false;
    } else {

      _Variable *rec = CheckReceptacleCommandID(
          &AppendContainerName(*arg2, currentProgram.nameSpacePrefix),
          HY_HBL_COMMAND_LFCOMPUTE, true);

      if (!rec) {
        return false;
      }

      rec->SetValue(new _Constant(lf->Compute()), false);

    }
  }

  return true;

}

//______________________________________________________________________________
bool
_ElementaryCommand::HandleSelectTemplateModel(_ExecutionList &currentProgram) {

  currentProgram.currentCommand++;

  SetStatusLine("Waiting for model selection");

  _String modelFile, errMsg;

  ReadModelList();

  if (((_String *)parameters(0))->Equal(&useLastModel)) {
    if (lastModelUsed.sLength) {
      PushFilePath(lastModelUsed);
    } else {
      WarnError(_String("First call to SelectTemplateModel. ") & useLastModel & " is meaningless.");
      return false;
    }
  } else {
    _String filterName(
        currentProgram.AddNameSpaceToID(*(_String *)parameters(0)));
    long objectType = HY_BL_DATASET_FILTER;
    _DataSetFilter *thisDF = (_DataSetFilter *)_HYRetrieveBLObjectByName(
        filterName, objectType, nil, true);
    // decide what this DF is comprised of

    _String dataType;
    long dataDimension = thisDF->GetDimension(),
         unitLength = thisDF->GetUnitLength();

    _TranslationTable *thisTT = thisDF->GetData()->GetTT();

    if (unitLength == 1) {
      if (thisTT->CheckType(HY_TRANSLATION_TABLE_STANDARD_NUCLEOTIDE)) {
        dataType = "nucleotide";
      } else if (thisTT->CheckType(HY_TRANSLATION_TABLE_STANDARD_PROTEIN)) {
        dataType = "aminoacid";
      }
    } else {
      if (thisTT->CheckType(HY_TRANSLATION_TABLE_STANDARD_NUCLEOTIDE)) {
        if (unitLength == 3) {
          dataType = "codon";
        } else {
          if (unitLength == 2) {
            dataType = "dinucleotide";
          }
        }
      }
    }

    if (!dataType.sLength) {
      WarnError(_String("DataSetFilter '") & filterName &
                "' contains non-standard data and SelectTemplateModel is not "
                "applicable.");
      return false;
    }

    _SimpleList matchingModels;

    for (unsigned long model_index = 0; model_index < templateModelList.lLength;
         model_index++) {
      _List *model_components = (_List *)templateModelList(model_index);

      if (dataType.Equal((_String *)model_components->GetItem(3))) {
        _String *dim = (_String *)model_components->GetItem(2);
        if (*dim == _String("*") || dataDimension == dim->toNum()) {
          matchingModels << model_index;
        }
      }
    }

    if (!matchingModels.lLength) {
      WarnError((_String)("DataSetFilter '") & filterName &
                "' could not be matched with any template models.");
      return false;
    }
    unsigned long model_id = HY_MAX_LONG_VALUE;

    if (currentProgram.stdinRedirect) {
      errMsg = currentProgram.FetchFromStdinRedirect();
      for (model_id = 0; model_id < matchingModels.lLength; model_id++)
        if (errMsg.Equal((_String *)(*(_List *)templateModelList(matchingModels(model_id)))(0))) {
          break;
        }

      if (model_id >= matchingModels.lLength) {
        WarnError(errMsg & " is not a valid model (with input redirect) in "
                           "call to SelectTemplateModel");
        return false;
      }
    } else {
#ifdef __HEADLESS__
      WarnError("Unhandled standard input interaction in SelectTemplateModel "
                "for headless HyPhy");
      return false;
#else
#if defined __UNIX__ && !defined __HYPHYQT__
      while (model_id == HY_MAX_LONG_VALUE) {
        printf("\n\n               +--------------------------+\n");
            printf("               | Select a standard model. |\n");
            printf("               +--------------------------+\n\n\n");

        for (model_id = 0; model_id < matchingModels.lLength; model_id++) {
          printf("\n\t(%s):%s", ((_String *)(*(_List *)templateModelList(matchingModels(model_id)))(0))->getStr(),
                 ((_String *)(*(_List *)templateModelList(matchingModels(model_id)))(1))->getStr());
        }
        printf("\n\n Please type in the abbreviation for the model you want to use:");
        dataType.CopyDynamicString(StringFromConsole());
        dataType.UpCase();
        for (model_id = 0; model_id < matchingModels.lLength; model_id++) {
          if (dataType.Equal((_String *)(
                  *(_List *)templateModelList(matchingModels(model_id)))(0))) {
            break;
          }
        }
        if (model_id == matchingModels.lLength) {
          model_id = HY_MAX_LONG_VALUE;
        }
      }
#endif
#if !defined __UNIX__ || defined __HYPHYQT__
      _SimpleList choiceDummy(2, 0, 1), selDummy;
      model_id = HandleListSelection(
          templateModelList, choiceDummy, matchingModels,
          "Choose one of the standard substitution models", selDummy, 1, nil);
      if (model_id == -1) {
        terminateExecution = true;
        return false;
      }
#endif
#endif
    }
    modelFile = _HYStandardDirectory(HY_HBL_DIRECTORY_TEMPLATE_MODELS) & *((_String *)(*(_List *)templateModelList(matchingModels(model_id)))(4));
    PushFilePath(modelFile, false);
  }

  _ExecutionList stdModel;
  if (currentProgram.nameSpacePrefix) {
    stdModel.SetNameSpace(*currentProgram.nameSpacePrefix->GetName());
  }

  ReadBatchFile(modelFile, stdModel);
  PopFilePath();
  lastModelUsed = modelFile;

  stdModel.stdinRedirectAux = currentProgram.stdinRedirectAux;
  stdModel.stdinRedirect = currentProgram.stdinRedirect;
  stdModel.Execute();
  stdModel.stdinRedirectAux = nil;
  stdModel.stdinRedirect = nil;

  return true;

}

//______________________________________________________________________________
bool _ElementaryCommand::HandleUseModel(_ExecutionList &currentProgram) {

  currentProgram.currentCommand++;
  _String namedspacedMM(currentProgram.AddNameSpaceToID(*(_String *)parameters(0)));
  long mID = FindModelName(namedspacedMM);

  if (mID < 0 && !useNoModel.Equal((_String *)parameters(0))) {
    WarnError(
        *(_String *)parameters(0) &
        " does not refer to a valid defined substitution model in call to " &
        _HY_ValidHBLExpressions.RetrieveKeyByPayload(HY_HBL_COMMAND_USE_MODEL));
    return false;
  } else {
    lastMatrixDeclared = mID;
  }

  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::HandleSetParameter(_ExecutionList &currentProgram) {

  currentProgram.currentCommand++;
  /*
      first check to see if matrix parameters here are valid
  */

  _String *currentArgument = (_String *)parameters(0),
          nmspc = AppendContainerName(*currentArgument,
                                      currentProgram.nameSpacePrefix),
          errMsg, result;

  if (currentArgument->Equal(&randomSeed)) {
    globalRandSeed = ProcessNumericArgument((_String *)parameters(1), currentProgram.nameSpacePrefix);
    init_genrand(globalRandSeed);
    setParameter(randomSeed, ((long) globalRandSeed));
    return true;
  }

  if (currentArgument->Equal(&randomSeed)) {
    globalRandSeed = ProcessNumericArgument((_String *)parameters(1), currentProgram.nameSpacePrefix);
    init_genrand(globalRandSeed);
    setParameter(randomSeed, ((long) globalRandSeed));
    return true;
  }

  if (currentArgument->Equal(&deferConstrainAssignment)) {
    bool on = ProcessNumericArgument((_String *)parameters(1),currentProgram.nameSpacePrefix);
    if (on) {
      deferSetFormula = (_SimpleList *)checkPointer(new _SimpleList);
    } else if (deferSetFormula) {
      FinishDeferredSF();
    }
    return true;
  }

  if (currentArgument->Equal(&_hyExecutionErrorMode)) {
    currentProgram.errorHandlingMode = ProcessNumericArgument((_String *)parameters(1), currentProgram.nameSpacePrefix);
    return true;
  }

  if (currentArgument->Equal(&statusBarProgressValue)) {
#if !defined __UNIX__
    SetStatusLine(empty, empty, empty,
                  ProcessNumericArgument((_String *)parameters(1),
                                         currentProgram.nameSpacePrefix),
                  HY_SL_PERCENT);
#endif
    return true;
  }

  if (currentArgument->Equal(&statusBarUpdateString)) {
    _String sbar_value = ProcessLiteralArgument((_String *)parameters(1),
                                                currentProgram.nameSpacePrefix);

#if defined __UNIX__
#if not defined __HYPHY_GTK__ &&not defined __HEADLESS__
    SetStatusLineUser(sbar_value);
#else
    SetStatusLine(sbar_value);
#endif
#else
    SetStatusLine(empty, sbar_value, empty, 0, HY_SL_TASK);
#endif
    return true;
  }

  long objectIndex, typeFlag = HY_BL_ANY;

  BaseRef theObject = _HYRetrieveBLObjectByName(nmspc, typeFlag, &objectIndex);

  switch (typeFlag) {
  case HY_BL_BGM: { // BGM Branch
    currentArgument = (_String *)parameters(1);

    _BayesianGraphicalModel *lkf = (_BayesianGraphicalModel *)theObject;
    // set data matrix
    if (currentArgument->Equal(&bgmData)) {
      _Matrix *dataMx = (_Matrix *)FetchObjectFromVariableByType(
          &AppendContainerName(*(_String *)parameters(2),currentProgram.nameSpacePrefix),
          MATRIX, HY_HBL_COMMAND_SET_PARAMETER);

      if (dataMx) {
        long num_nodes = ((_BayesianGraphicalModel *)lkf)->GetNumNodes();

        if (dataMx->GetVDim() == num_nodes) {
          ((_BayesianGraphicalModel *)lkf)->SetDataMatrix((_Matrix *)dataMx);
        } else {
          currentProgram.ReportAnExecutionError(
              _String("Data matrix columns (") & dataMx->GetVDim() & " ) does not match number of nodes in graph (" & num_nodes & ")");
          return false;
        }
      } else {
        return false;
      }

    } else if (currentArgument->Equal(&bgmScores)) {
        // restore node score cache
      _AssociativeList *cacheAVL =
          (_AssociativeList *)FetchObjectFromVariableByType(
              &AppendContainerName(*(_String *)parameters(2),currentProgram.nameSpacePrefix),ASSOCIATIVE_LIST, HY_HBL_COMMAND_SET_PARAMETER);
      if (cacheAVL) {
        ((_BayesianGraphicalModel *)lkf)->ImportCache(cacheAVL);
      } else {
        return false;
      }
    } else if (currentArgument->Equal(&bgmGraph)) {
        // set structure to user-specified adjacency matrix
        _Matrix *graphMx = (_Matrix *)FetchObjectFromVariableByType(
            &AppendContainerName(*(_String *)parameters(2), currentProgram.nameSpacePrefix),
            MATRIX, HY_HBL_COMMAND_SET_PARAMETER);

      if (graphMx) {
        long num_nodes = ((_BayesianGraphicalModel *)lkf)->GetNumNodes();

        if (graphMx->GetHDim() == num_nodes && graphMx->GetVDim() == num_nodes) {
            ((_BayesianGraphicalModel *)lkf)->SetStructure((_Matrix *)graphMx->makeDynamic());

        } else {
          currentProgram.ReportAnExecutionError(
              "Dimension of graph does not match current graph");
          return false;
        }
      } else {
        return false;
      }

    } else if (currentArgument->Equal(&bgmConstraintMx)) {
      // set constraint matrix
      _Matrix *constraintMx = (_Matrix *)FetchObjectFromVariableByType(
          &AppendContainerName(*(_String *)parameters(2), currentProgram.nameSpacePrefix),
          MATRIX, HY_HBL_COMMAND_SET_PARAMETER);

      if (constraintMx) {
        long num_nodes = ((_BayesianGraphicalModel *)lkf)->GetNumNodes();

        if (constraintMx->GetHDim() == num_nodes && constraintMx->GetVDim() == num_nodes) {
          ((_BayesianGraphicalModel *)lkf)->SetConstraints((_Matrix *)constraintMx->makeDynamic());
        } else {
          currentProgram.ReportAnExecutionError(
              "Dimensions of constraint matrix do not match current graph");
          return false;
        }
      } else {
        return false;
      }
    } else if (currentArgument->Equal(&bgmNodeOrder)) {
       // set node order
      _Matrix *orderMx = (_Matrix *)FetchObjectFromVariableByType(
          &AppendContainerName(*(_String *)parameters(2), currentProgram.nameSpacePrefix),
          MATRIX, HY_HBL_COMMAND_SET_PARAMETER);

      if (orderMx) {
        // UNDER DEVELOPMENT  April 17, 2008 afyp
        long num_nodes = ((_BayesianGraphicalModel *)lkf)->GetNumNodes();

        _SimpleList *orderList = new _SimpleList();

        orderList->Populate(num_nodes, 0, 0);

        if (orderMx->GetVDim() == num_nodes) {
          for (long i = 0; i < num_nodes; i++) {
            orderList->lData[i] = (long)((*orderMx)(0, i));
          }

          ((_BayesianGraphicalModel *)lkf)->SetNodeOrder((_SimpleList *)orderList->makeDynamic());

        } else {
          currentProgram.ReportAnExecutionError("Length of order vector doesn't match number of nodes in graph");
          return false;
        }
      } else {
        return false;
      }
    } else if (currentArgument->Equal(&bgmParameters)) {
        // set network parameters
      _AssociativeList *inAVL =
          (_AssociativeList *)FetchObjectFromVariableByType(
              &AppendContainerName(*(_String *)parameters(2),
                                   currentProgram.nameSpacePrefix),
              ASSOCIATIVE_LIST, HY_HBL_COMMAND_SET_PARAMETER);
      if (inAVL) {
        ((_BayesianGraphicalModel *)lkf)->ImportCache(inAVL);
      } else {
        return false;
      }
    } else {
      // anything else
      currentProgram.ReportAnExecutionError(*currentArgument & " is not a valid BGM parameter");
      return false;
    }
  } // end BGM
      break;

  case HY_BL_SCFG:
  case HY_BL_LIKELIHOOD_FUNCTION: {
    currentArgument = (_String *)parameters(1);
    if (typeFlag == HY_BL_SCFG && currentArgument->Equal(&scfgCorpus)) {
      ((Scfg *)theObject)->SetStringCorpus((_String *)parameters(2));
    } else {
      _LikelihoodFunction *lkf = (_LikelihoodFunction *)theObject;
      currentArgument = (_String *)parameters(1);
      long g = ProcessNumericArgument(currentArgument, currentProgram.nameSpacePrefix);

      if (g < 0 || g >= lkf->GetIndependentVars().lLength) {
        currentProgram.ReportAnExecutionError(*currentArgument & " (=" & g & ") is not a valid parameter index");
        return false;
      }
      currentArgument = (_String *)parameters(2);
      lkf->SetIthIndependent(g, ProcessNumericArgument(currentArgument, currentProgram.nameSpacePrefix));
    }
  } break;
  // end SCFG and LF

  case HY_BL_DATASET:
  case HY_BL_DATASET_FILTER: {
    _DataSet *ds = nil;
    long f = ProcessNumericArgument((_String *)parameters(1), currentProgram.nameSpacePrefix);
    if (typeFlag == HY_BL_DATASET) {
      ds = (_DataSet *)theObject;
    } else {
      _DataSetFilter *dsf = (_DataSetFilter *)theObject;
      ds = dsf->GetData();
      if (f >= 0 && f < dsf->theNodeMap.lLength) {
        f = dsf->theNodeMap.lData[f];
      } else
        f = -1;
    }

    _List *dsNames = &ds->GetNames();

    if (f < 0 || f >= dsNames->lLength) {
      currentProgram.ReportAnExecutionError(*((_String *)parameters(1)) &
                                            " (=" & f &
                                            ") is not a valid sequence index");
      return false;
    }

    dsNames->Replace(f, new _String(ProcessLiteralArgument((_String *)parameters(2), currentProgram.nameSpacePrefix)), false);
  } // end data set and data set filter
      break;

  // Dataset and Datasetfilter
  default:
    // check to see if this is a calcnode
    _CalcNode *treeNode = (_CalcNode *)FetchObjectFromVariableByType(&nmspc, TREE_NODE);
    if (treeNode) {
      if (*((_String *)parameters(1)) == _String("MODEL")) {

        _String modelName = AppendContainerName(*((_String *)parameters(2)), currentProgram.nameSpacePrefix);
        long modelType = HY_BL_MODEL, modelIndex;
        BaseRef modelObject = _HYRetrieveBLObjectByName(modelName, modelType, &modelIndex, true);

        if (modelObject) {
          _VariableContainer *parentTree = treeNode->ParentTree();
          if (!parentTree) {
            currentProgram.ReportAnExecutionError(
                *((_String *)parameters(0)) &
                " is an orphaned tree node (the parent tree has been deleted)");
            return false;
          }
          long pID,
              lfID = ((_TheTree *)parentTree->Compute())->IsLinkedToALF(pID);
          if (lfID >= 0) {
            currentProgram.ReportAnExecutionError(
                (*parentTree->GetName()) &
                " is linked to a likelihood function (" &
                *_HBLObjectNameByType(HY_BL_LIKELIHOOD_FUNCTION, lfID) &
                ") and cannot be modified ");
            return false;
          }

          treeNode->ReplaceModel(modelName, parentTree);
          break;

        } else {
          currentProgram.ReportAnExecutionError(*((_String *)parameters(2)) & " does not appear to be a valid model name");
          return false;
        }
      } else {
        currentProgram.ReportAnExecutionError(
            *((_String *)parameters(1)) &
            " is not a supported parameter type for a tree node argument");
        return false;
      }
    }

    currentProgram.ReportAnExecutionError(
        *currentArgument & " is not a valid likelihood function/data set "
                           "filter/tree topology/tree node");
    return false;

  } // end cases
  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::HandleAssert(_ExecutionList &currentProgram) {
  currentProgram.currentCommand++;

  _String assertion(*(_String *)parameters(0));

  _Formula rhs, lhs;
  _FormulaParsingContext fpc(nil, currentProgram.nameSpacePrefix);
  if (Parse(&rhs, assertion, fpc, &lhs) == HY_FORMULA_EXPRESSION) {
    _PMathObj assertionResult = rhs.Compute();
    if (assertionResult && assertionResult->ObjectClass() == NUMBER) {
      if (CheckEqual(assertionResult->Value(), 0.0)) {
        _Parameter whatToDo;
        checkParameter(assertionBehavior, whatToDo, 0.0);

        _String errMsg;

        if (parameters.lLength == 1) {
          errMsg =
              _String("Assertion '") & *(_String *)parameters(0) & "' failed.";
        } else {
          errMsg = ProcessLiteralArgument((_String *)parameters(1),
                                          currentProgram.nameSpacePrefix);
        }

        if (CheckEqual(whatToDo, 1.)) {
          StringToConsole(errMsg);
          NLToConsole();
          currentProgram.GoToLastInstruction();
        } else {
          currentProgram.ReportAnExecutionError(errMsg);
          return false;
        }
      }
      return true;
    }
  }
  currentProgram.ReportAnExecutionError(
      _String("Assertion statement '") & *(_String *)parameters(0) &
      "' could not be computed or was not numeric.");

  return false;
}

//______________________________________________________________________________
bool _ElementaryCommand::HandleRequireVersion(_ExecutionList &currentProgram) {
  currentProgram.currentCommand++;
  _String theVersion = ProcessLiteralArgument((_String *)parameters(0), currentProgram.nameSpacePrefix);

  if (__KERNEL__VERSION__.toNum() < theVersion.toNum()) {
    currentProgram.ReportAnExecutionError(
        _String("Current batch file requires at least version :") & theVersion &
        " of HyPhy. Please download an updated version from "
        "http://www.hyphy.org and try again.");

    return false;
  }

  return true;

}

//______________________________________________________________________________
bool _ElementaryCommand::HandleDeleteObject(_ExecutionList &currentProgram) {

  currentProgram.currentCommand++;
  for (unsigned long objCount = 0; objCount < parameters.lLength; objCount++) {

    long objectType = HY_BL_LIKELIHOOD_FUNCTION, f = -1;

    BaseRef sourceObject = _HYRetrieveBLObjectByName(
        AppendContainerName(*(_String *)parameters(objCount), currentProgram.nameSpacePrefix),
        objectType, &f, false);

    if (sourceObject) {
      KillLFRecord(f, true);
    }

  }

  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::HandleClearConstraints(_ExecutionList &currentProgram) {
  currentProgram.currentCommand++;
  for (unsigned long i = 0; i < parameters.lLength; i++) {
    _String cName(currentProgram.AddNameSpaceToID(*(_String *)parameters(i)));
    long cID = LocateVarByName(cName);
    if (cID >= 0) { // variable exists
      FetchVar(cID)->ClearConstraints();
    }
  }
  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::HandleMolecularClock(_ExecutionList &currentProgram) {
  currentProgram.currentCommand++;

  _String theBaseNode(
      currentProgram.AddNameSpaceToID(*(_String *)parameters(0))),
      treeName;

  _Variable *theObject = FetchVar(LocateVarByName(theBaseNode));

  if (!theObject || (theObject->ObjectClass() != TREE &&
                     theObject->ObjectClass() != TREE_NODE)) {

    WarnError(_String("Not a defined tree/tree node object '") & theBaseNode &
              "' in call to " & _HY_ValidHBLExpressions.RetrieveKeyByPayload(
                                    HY_HBL_COMMAND_MOLECULAR_CLOCK));
    return false;
  }

  _TheTree *theTree = nil;
  if (theObject->ObjectClass() == TREE_NODE) {
    theTree = (_TheTree *)((_VariableContainer *)theObject)->GetTheParent();
    if (!theTree) {
      WarnError(_String("Internal error - orphaned tree node '") & theBaseNode &
                "' in call to " & _HY_ValidHBLExpressions.RetrieveKeyByPayload(
                                      HY_HBL_COMMAND_MOLECULAR_CLOCK));
      return false;

    }
    treeName = *theTree->GetName();
    theBaseNode = theObject->GetName()->Cut(treeName.sLength + 1, -1);
  } else {
    treeName = *theObject->GetName();
    theTree = (_TheTree *)theObject;
    theBaseNode = empty;
  }

  theTree->MolecularClock(theBaseNode, parameters);
  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::HandleGetURL(_ExecutionList &currentProgram) {
  currentProgram.currentCommand++;

  _String url(ProcessLiteralArgument((_String *)parameters(1),
                                     currentProgram.nameSpacePrefix)),
      *arg1 = (_String *)parameters(0),
      *act = parameters.lLength > 2 ? (_String *)parameters(2) : nil, errMsg;

  if (act == nil) {
    _Variable *rec = CheckReceptacleCommandID(
        &AppendContainerName(*arg1, currentProgram.nameSpacePrefix),
        HY_HBL_COMMAND_GET_URL, true, false, &currentProgram);

    if (!rec) {
      return false;
    }

    if (Get_a_URL(url)) {
      rec->SetValue(new _FString(url, false), false);
    } else {
      errMsg = _String("Could not fetch '") & url & "'";
    }

  } else {
    if (act->Equal(&getURLFileFlag)) {
      _String fileName(
          ProcessLiteralArgument(arg1, currentProgram.nameSpacePrefix));
      fileName.ProcessFileName(true, false,
                               (Ptr) currentProgram.nameSpacePrefix);
      if (!Get_a_URL(url, &fileName)) {
        errMsg = _String("Could not fetch '") & url & "'";
      }
    } else {
      errMsg = "Unknown action flag";
    }
  }

  if (errMsg.sLength) {
    currentProgram.ReportAnExecutionError(errMsg);
    return false;
  }

  return true;

}

//______________________________________________________________________________
bool _ElementaryCommand::HandleGetString(_ExecutionList &currentProgram) {
  currentProgram.currentCommand++;

  _String errMsg, *result = nil;

  long f, sID, sID2 = -1;

  _Variable *theReceptacle = CheckReceptacleCommandID(
      &AppendContainerName(*((_String *)parameters(0)),
                           currentProgram.nameSpacePrefix),
      HY_HBL_COMMAND_GET_STRING, true, false, &currentProgram);

  if (!theReceptacle) {
    return false;
  }

  sID = ProcessNumericArgument((_String *)parameters(2), currentProgram.nameSpacePrefix);

  if (parameters.lLength > 3) {
    sID2 = ProcessNumericArgument((_String *)parameters(3), currentProgram.nameSpacePrefix);
  }

  f = _HY_GetStringGlobalTypes.Find((_String *)parameters(1));

  if (f >= 0) {
    f = _HY_GetStringGlobalTypes.GetXtra(f);
  }

  switch (f) {

  case HY_BL_LIKELIHOOD_FUNCTION: // LikelihoodFunction
  case HY_BL_DATASET:
  case HY_BL_DATASET_FILTER:
  case HY_BL_SCFG:
  case HY_BL_BGM: {
    result = (_String *)_HBLObjectNameByType(f, sID);
    if (result) {
      result = (_String *)result->makeDynamic();
      //ReportWarning(_String((const char*)"In HandleGetString(): ") & result);
    }
    break;
  }

  case HY_BL_HBL_FUNCTION: // UserFunction
    result = (_String *)_HBLObjectNameByType(HY_BL_HBL_FUNCTION, sID);
    if (result) {
      _AssociativeList *resAVL =
          (_AssociativeList *)checkPointer(new _AssociativeList);
      resAVL->MStore("ID", new _FString(*result), false);
      resAVL->MStore(
          "Arguments",
          new _Matrix(*(_List *)batchLanguageFunctionParameterLists(sID)),
          false);
      theReceptacle->SetValue(resAVL, false);
      return true;
    }
    break;

  case HY_BL_TREE: { // Tree
                     // 20110608 SLKP: REFACTOR into a separate function
    // I am sure this is used somewhere else (perhaps for other types)
    result = FetchMathObjectNameOfTypeByIndex(TREE, sID);
    if (result) {
      result = (_String *)result->makeDynamic();
    }
    break;
  }

  default: { // everything else...
             // decide what kind of object current argument represents

    _String *currentArgument = (_String *)parameters(1),
            nmspaced = AppendContainerName(*currentArgument,
                                           currentProgram.nameSpacePrefix);
    long typeFlag = HY_BL_ANY, index = -1;

    BaseRef theObject = _HYRetrieveBLObjectByName(nmspaced, typeFlag, &index);

    if (theObject) {
      switch (typeFlag) {
      case HY_BL_DATASET: {
        _DataSet *dataSetObject = (_DataSet *)theObject;
        if (sID >= 0 && sID < dataSetObject->NoOfSpecies()) {
          result = (_String *)(dataSetObject->GetNames())(sID)->makeDynamic();
        } else {
          if (sID < 0) {
            theReceptacle->SetValue(new _Matrix(dataSetObject->GetNames()), false);
            return true;
          }
        }
        break;
      }
      case HY_BL_DATASET_FILTER: {
        _DataSetFilter *dataSetFilterObject = (_DataSetFilter *)theObject;

        if (sID >= 0 && sID < dataSetFilterObject->NumberSpecies()) {
          result = (_String *)(dataSetFilterObject->GetData()->GetNames())(
              dataSetFilterObject->theNodeMap(sID))->makeDynamic();
        } else {
          if (sID < 0) {
            _List filterSeqNames,
                *originalNames = &dataSetFilterObject->GetData()->GetNames();
            for (long seqID = 0; seqID < dataSetFilterObject->NumberSpecies(); seqID++) {
              filterSeqNames << (*originalNames)(dataSetFilterObject->theNodeMap(seqID));
            }
            theReceptacle->SetValue(new _Matrix(filterSeqNames), false);
            return true;
          }
        }
        break;
      }
      case HY_BL_BGM: {
        //ReportWarning(_String("In HandleGetString() for case HY_BL_BGM"));
        _BayesianGraphicalModel *this_bgm = (_BayesianGraphicalModel *)theObject;

        switch (sID) {
        case HY_HBL_GET_STRING_BGM_SCORE:
            { // return associative list containing node score cache
          _AssociativeList *export_alist = new _AssociativeList;

          if (this_bgm->ExportCache(export_alist)) {
            theReceptacle->SetValue(export_alist, false);
            return true;
          } else {
            DeleteObject(export_alist);
            errMsg = _String("Failed to export node score cache for BGM '") &
                     nmspaced & "'";
          }
          break;
        }
        case HY_HBL_GET_STRING_BGM_SERIALIZE:
            { // return associative list with network structure and parameters
          result = new _String(1024L, true);
          this_bgm->SerializeBGM(*result);
          result->Finalize();
          theReceptacle->SetValue(new _FString(result), false);
          return true;
        }
        default: {
          errMsg = _String("Unrecognized index ") & sID & " for a BGM object";
        }
        }
      }
      case HY_BL_LIKELIHOOD_FUNCTION:
      case HY_BL_SCFG: {

        _LikelihoodFunction *lf = (_LikelihoodFunction *)theObject;
        if (sID >= 0) {
          if (sID < lf->GetIndependentVars().lLength) {
            result = (_String *)(LocateVar(lf->GetIndependentVars().lData[sID])->GetName())->makeDynamic();
          } else {
            if (sID < lf->GetIndependentVars().lLength + lf->GetDependentVars().lLength) {
              result = (_String *)(LocateVar(lf->GetDependentVars().lData[sID - lf->GetIndependentVars().lLength])->GetName())->makeDynamic();
            }
          }
        } else {
          _AssociativeList *resList = lf->CollectLFAttributes();
          if (typeFlag == HY_BL_SCFG) {
            ((Scfg *)lf)->AddSCFGInfo(resList);
          }
          theReceptacle->SetValue(resList, false);
          return true;
        }
        break;
      }

      case HY_BL_MODEL: {
        if (sID >= 0) {
          // check to make see if the
          if (sID2 < 0) { // get the sID's parameter name
            _SimpleList modelP;
            _AVLList modelPA(&modelP);
            ScanModelForVariables(index, modelPA, false, -1, false);
            modelPA.ReorderList();
            if (sID < modelP.lLength) {
              result = (_String *)LocateVar(modelP.lData[sID])->GetName()->makeDynamic();
            }

          } else { // get the formula for cell (sID, sID2)
            if (!IsModelOfExplicitForm(index)) {
              _Variable *theMx = (_Variable *)theObject;
              _Formula *cellFla = ((_Matrix *)theMx->GetValue())->GetFormula(sID, sID2);
              if (cellFla) {
                result = new _String((_String *)cellFla->toStr());
              }
            }
          }

        } else {
          _Variable *tV, *tV2;
          bool mByF;
          RetrieveModelComponents(index, tV, tV2, mByF);

          if (tV) {
            if (sID == -1) { // branch length expression
              result = ((_Matrix *)tV->GetValue())->BranchLengthExpression((_Matrix *)tV2->GetValue(), mByF);
            } else /*
                       returns an AVL with keys
                       "RATE_MATRIX" - the ID of the rate matrix
                       "EQ_FREQS"    - the ID of eq. freq. vector
                       "MULT_BY_FREQ" - a 0/1 flag to determine which format the
                       matrix is in.
                   */
                {
              _AssociativeList *resList = new _AssociativeList;
              resList->MStore("RATE_MATRIX", new _FString(*tV->GetName()),false);
              resList->MStore("EQ_FREQS", new _FString(*tV2->GetName()), false);
              resList->MStore("MULT_BY_FREQ", new _Constant(mByF), false);
              theReceptacle->SetValue(resList, false);
              return true;
            }
          }
        }
        break;
      }
      case HY_BL_HBL_FUNCTION: {
        _AssociativeList *resAVL = (_AssociativeList *)checkPointer(new _AssociativeList);
        resAVL->MStore("ID", new _FString(*_HBLObjectNameByType(HY_BL_HBL_FUNCTION, index, false)), false);
        resAVL->MStore("Arguments",new _Matrix(*(_List *)batchLanguageFunctionParameterLists(index)),false);
        resAVL->MStore("Body", new _FString(((_ExecutionList *)batchLanguageFunctions(index))->sourceText, false), false);
        theReceptacle->SetValue(resAVL, false);
        return true;
      }
      } // end of "switch"
    } else {
      if (currentArgument->Equal(&versionString)) {
        if (sID > 1.5) {
#ifdef __HEADLESS__
          result = new _String(_String("Library version ") & __KERNEL__VERSION__);
#else
#ifdef __MAC__
        result = new _String(_String("Macintosh ") & __KERNEL__VERSION__);
#else
#ifdef __WINDOZE__
        result = new _String(_String("Windows ") & __KERNEL__VERSION__);
#else
        result = new _String(_String("Source ") & __KERNEL__VERSION__);
#endif
#endif
#endif
        } else if (sID > 0.5) {
          result = new _String(GetVersionString());
        }
        else {
          result = new _String(__KERNEL__VERSION__);
        }
      } else if (currentArgument->Equal(&timeStamp)) {
        result = new _String(GetTimeStamp(sID < 0.5));
      } else {
        _Variable *theVar = FetchVar(LocateVarByName(*currentArgument));
        if (theVar) {
          if (theVar->IsIndependent()) {
            result = (_String *)theVar->toStr();
          } else {
            if (sID == -1)
                // list of variables
                {
              _SimpleList vL;
              _AVLList vAVL(&vL);
              theVar->ScanForVariables(vAVL, true);
              vAVL.ReorderList();
              _AssociativeList *resL =
                  (_AssociativeList *)checkPointer(new _AssociativeList);
              _List splitVars;
              SplitVariableIDsIntoLocalAndGlobal(vL, splitVars);
              InsertVarIDsInList(resL, "Global", *(_SimpleList *)splitVars(0));
              InsertVarIDsInList(resL, "Local", *(_SimpleList *)splitVars(1));

              theReceptacle->SetValue(resL, false);
              return true;
            } else { // formula string
              result = (_String *)theVar->GetFormulaString();
            }
          }
        } else {
          errMsg = _String("'") & *currentArgument &
                   "' is not an allowed argument type ";
        }
      }
    }
  }
  }

  if (errMsg.sLength) {
    currentProgram.ReportAnExecutionError(errMsg);
    DeleteObject(result);
    result = nil;
  }

  if (result) {
    theReceptacle->SetValue(new _FString(result), false);
    return true;
  }

  theReceptacle->SetValue(new _MathObject(), false);

  return false;
}

//______________________________________________________________________________
bool _ElementaryCommand::HandleExport(_ExecutionList &currentProgram) {

  currentProgram.currentCommand++;

  _String objectID(currentProgram.AddNameSpaceToID(*(_String *)parameters(1))),
      arg1(currentProgram.AddNameSpaceToID(*(_String *)parameters(0))), errMsg;

  _Variable *theReceptacle = CheckReceptacleCommandID(
      &AppendContainerName(arg1, currentProgram.nameSpacePrefix),
      HY_HBL_COMMAND_EXPORT, true, false, &currentProgram);
  if (!theReceptacle) {
    return false;
  }

  _FString *outLF = new _FString(new _String(8192L, 1));
  checkPointer(outLF);
  long typeFlag =
           HY_BL_MODEL | HY_BL_LIKELIHOOD_FUNCTION | HY_BL_DATASET_FILTER,
       index;

  BaseRef objectToExport = _HYRetrieveBLObjectByName(objectID, typeFlag, &index);
  if (!objectToExport) {
    errMsg = _String("'") & objectID & "' is not a supported type";
  } else {
    switch (typeFlag) {
    case HY_BL_LIKELIHOOD_FUNCTION: {
      ((_LikelihoodFunction *)objectToExport)->SerializeLF(*outLF->theString);
      outLF->theString->Finalize();
      break;
    }
    case HY_BL_DATASET_FILTER: {
      outLF->theString->Finalize();
      DeleteObject(outLF->theString);
      checkPointer(outLF->theString = new _String(
          (_String *)((_DataSetFilter *)objectToExport)->toStr()));
      break;
    }
    case HY_BL_MODEL: {
      SerializeModel(*outLF->theString, index, nil, true);
      outLF->theString->Finalize();
      break;
    }
    }
  }

  if (errMsg.sLength) {
    outLF->theString->Finalize();
    DeleteObject(outLF);
    currentProgram.ReportAnExecutionError(errMsg);
    theReceptacle->SetValue(new _MathObject, false);
    return false;
  }

  theReceptacle->SetValue(outLF, false);
  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::HandleDifferentiate(_ExecutionList &currentProgram) {

  currentProgram.currentCommand++;

  _String arg1(currentProgram.AddNameSpaceToID(*(_String *)parameters(0))),
      errMsg, expressionToParse = *(_String *)parameters(1);

  _Formula *theResult = nil;

  _Variable *theReceptacle = CheckReceptacleCommandID(
      &AppendContainerName(arg1, currentProgram.nameSpacePrefix),
      HY_HBL_COMMAND_DIFFERENTIATE, true, false, &currentProgram);
  if (!theReceptacle) {
    return false;
  }

  _Formula theExpression(expressionToParse, currentProgram.nameSpacePrefix, &errMsg);

  if (!theExpression.IsEmpty() && errMsg.sLength == 0) {
    long times = 1;
    if (parameters.lLength == 4) {
      times = ProcessNumericArgument((_String *)parameters(3), currentProgram.nameSpacePrefix, &currentProgram);
      if (!numericalParameterSuccessFlag) {
        return false;
      }
    }
    if (times <= 0) {
      errMsg = "The number of times to differentiate must be a non-negative integer";
    }

    theResult = theExpression.Differentiate(*(_String *)parameters(2), false);
    for (; times > 1 && theResult; times--) {
      _Formula *temp = theResult->Differentiate(*(_String *)parameters(2));
      delete (theResult);
      theResult = temp;
    }
  }

  if (errMsg.sLength || theResult == nil) {
    if (theResult) {
      delete (theResult);
    } else {
      errMsg = _String("Differentiation of '") & *(_String *)parameters(1) & "' failed";
    }
    currentProgram.ReportAnExecutionError(errMsg);
    theReceptacle->SetValue(new _MathObject, false);
    return false;
  }

  theReceptacle->SetFormula(*theResult);
  if (theResult)
    delete (theResult);

  return true;
}

//______________________________________________________________________________
bool _ElementaryCommand::HandleFprintf(_ExecutionList &currentProgram) {
  currentProgram.currentCommand++;
  _String *targetName = (_String *)parameters(0), fnm;

  bool doClose = true, print_to_stdout = false;

  FILE *dest = nil;

  try {
    bool skipFilePathEval = false;

    if (targetName->Equal(&stdoutDestination)) {
      _FString *redirect = (_FString *)FetchObjectFromVariableByType(&blFprintfRedirect, STRING);
      if (redirect && redirect->theString->sLength) {
        if (redirect->theString->Equal(&blFprintfDevNull)) {
          return true; // "print" to /dev/null
        } else {
          skipFilePathEval = true;
          targetName = redirect->theString;
        }
      } else {
        print_to_stdout = true;
      }
    }

    checkParameter(printDigitsSpec, printDigits, 0);

    if (!print_to_stdout) {
      fnm = *targetName;
      if (fnm.Equal(&messageLogDestination)) {
        if ((dest = globalMessageFile) == nil) {
          return true; // requested print to MESSAGE_LOG, but it does not exist
          // (e.g. remote MPI nodes, or running from a read only location
        }
      } else {
        if (skipFilePathEval == false && !fnm.IsALiteralArgument()) {
          fnm = GetStringFromFormula(&fnm, currentProgram.nameSpacePrefix);
        }

        if (!fnm.ProcessFileName(true, false, (Ptr) currentProgram.nameSpacePrefix, false, &currentProgram)) {
          return false;
        }

        long k = openFileHandles.Find(&fnm);
        doClose = k < 0;

        if (!doClose) {
          dest = (FILE *)openFileHandles.GetXtra(k);
        } else {
          if ((dest = doFileOpen(fnm.getStr(), "a")) == nil)
            throw(_String("Could not create/open output file at path '") & fnm & "'.");
        }
      }
    }

    for (long i = 1; i < parameters.lLength; i++) {
      _String *varname = (_String *)parameters(i);

      BaseRef thePrintObject = nil;
      _Formula f;

      if (varname->Equal(&clearFile)) {
        if (!print_to_stdout && dest) {
          fclose(dest);
          dest = doFileOpen(fnm.getStr(), "w");
          long k = openFileHandles.Find(&fnm);
          if (k >= 0) {
            openFileHandles.SetXtra(k, (long) dest);
          }
        }
      } else if (varname->Equal(&keepFileOpen) && !print_to_stdout) {
        if (openFileHandles.Find(&fnm) < 0) {
          openFileHandles.Insert(fnm.makeDynamic(), (long) dest);
        }
        doClose = false;
      } else if (varname->Equal(&closeFile) && !print_to_stdout) {
        openFileHandles.Delete(&fnm, true);
        doClose = true;
      } else if (varname->Equal(&systemVariableDump)) {
        thePrintObject = &variableNames;
      } else if (varname->Equal(&selfDump)) {
        thePrintObject = &currentProgram;
      } else {
        // check for possible string reference
        _String temp = ProcessStringArgument(varname), nmspace;
        if (temp.sLength > 0) {
          nmspace = AppendContainerName(temp, currentProgram.nameSpacePrefix);
          if (nmspace.IsValidIdentifier()) {
            thePrintObject = FetchObjectFromVariableByType(&nmspace, HY_ANY_OBJECT);
          }
        } else {
          nmspace = AppendContainerName(*varname, currentProgram.nameSpacePrefix);
        }

        if (thePrintObject == nil) {
          long typeFlag = HY_BL_ANY;

          thePrintObject = _HYRetrieveBLObjectByName(nmspace, typeFlag);

          if (!thePrintObject) {
            _String argCopy = *varname, errMsg;

            _FormulaParsingContext fpc(&errMsg, currentProgram.nameSpacePrefix);
            if (Parse(&f, argCopy, fpc) == HY_FORMULA_EXPRESSION) {
              _hyExecutionContext localContext (currentProgram.nameSpacePrefix);
              thePrintObject = f.Compute(0, &localContext);
            } else {
              if (errMsg.sLength)
                throw(errMsg);
              else
                throw(_String("Argument ") & i & " is not a simple expression");
            }
          }
        }
      }

      if (thePrintObject) {
        if (!print_to_stdout) {
          thePrintObject->toFileStr(dest);
        } else {
          _String outS((_String *)thePrintObject->toStr());
          StringToConsole(outS);
        }
      }
    }
  }
  catch (_String errMsg) {
    currentProgram.ReportAnExecutionError(errMsg);
  }

#if !defined __UNIX__ || defined __HEADLESS__ || defined __HYPHYQT__
  if (print_to_stdout) {
    yieldCPUTime();
  }
#endif
  if (dest && dest != globalMessageFile && doClose) {
    fclose(dest);
  }

  return !currentProgram.IsErrorState();
}
