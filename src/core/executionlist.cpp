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

#include "executionlist.h"
#include "likefunc.h"
#include "batchlan.h"
#include "string.h"
#include "ctype.h"
#include "polynoml.h"
#include "time.h"
#include "scfg.h"
#include "HYNetInterface.h"
#include "hy_globals.h"
#include "batchlan_globals.h"

#include "bayesgraph.h"

#if not defined __HEADLESS__ && defined __HYPHYQT__
    #include "HYSharedMain.h"
    #include "hyphy_qt_helpers.h"
#endif


_ExecutionList* currentExecutionList = nil;

_ExecutionList::_ExecutionList ()
{
    result              = nil;
    currentCommand      = 0;
    cli                 = nil;
    profileCounter      = nil;
    stdinRedirect       = nil;
    stdinRedirectAux    = nil;
    doProfile           = 0;
    nameSpacePrefix     = nil;
    
    if (currentExecutionList) {
        errorHandlingMode  = currentExecutionList->errorHandlingMode;
        errorState         = currentExecutionList->errorState;
    } else {
        errorHandlingMode = HY_BL_ERROR_HANDLING_DEFAULT;
        errorState = false;
    }

} // doesn't do much

//____________________________________________________________________________________
_ExecutionList::_ExecutionList (_String& source, _String* namespaceID, bool copySource, bool* successFlag)
{
    currentCommand = 0;
    result         = nil;
    cli            = nil;
    profileCounter = nil;
    doProfile      = 0;
    stdinRedirect  = nil;
    stdinRedirectAux = nil;
    nameSpacePrefix = nil;
    
    if (namespaceID) {
        SetNameSpace (*namespaceID);
    }
    if (copySource) {
        sourceText.Duplicate (&source);
    }
    if (currentExecutionList) {
        errorHandlingMode  = currentExecutionList->errorHandlingMode;
        errorState         = currentExecutionList->errorState;
    } else {
        errorHandlingMode = HY_BL_ERROR_HANDLING_DEFAULT;
        errorState = false;
    }

    bool result = BuildList (source, nil, false, true);
    if (successFlag) {
        *successFlag = result;
    }
}

//____________________________________________________________________________________

_ExecutionList::~_ExecutionList (void)
{
    if (cli) {
        delete [] cli->values;
        delete [] cli->stack;
        delete cli;
        cli = nil;
    }

    if (profileCounter) {
        DeleteObject (profileCounter);
        profileCounter = nil;
    }

    DeleteObject (stdinRedirect);
    DeleteObject (stdinRedirectAux);
    DeleteObject (nameSpacePrefix);

    ResetFormulae();
    DeleteObject (result);
}

//____________________________________________________________________________________

BaseRef     _ExecutionList::makeDynamic (void)
{
    _ExecutionList * Res = (_ExecutionList*)checkPointer(new _ExecutionList);

    memcpy ((char*)Res, (char*)this, sizeof (_ExecutionList));

    Res->nInstances         = 1;
    Res->Duplicate          (this);
    Res->cli                = nil;
    Res->profileCounter     = nil;
    Res->doProfile          = doProfile;
    Res->errorHandlingMode  = errorHandlingMode;
    Res->errorState         = errorState;

    if(result) {
        Res->result = (_PMathObj)result->makeDynamic();
    }

    return Res;
}

//____________________________________________________________________________________

void        _ExecutionList::Duplicate   (BaseRef source)
{
    _List::Duplicate    (source);

    _ExecutionList* s = (_ExecutionList*)source;

    if (s->result) {
        s->result=(_PMathObj)result->makeDynamic();
    }

    errorHandlingMode  = s->errorHandlingMode;
    errorState         = s->errorState;
}


//____________________________________________________________________________________
void    _ExecutionList::ReportAnExecutionError (_String errMsg, bool doCurrentCommand, bool appendToExisting) {
    if (doCurrentCommand) {
        _ElementaryCommand *theCommand = FetchLastCommand();
        if (theCommand) {
            errMsg = errMsg & " in call to " & _HY_ValidHBLExpressions.RetrieveKeyByPayload(theCommand->GetCode());
        }
    }
    errorState = true;
    switch (errorHandlingMode) {
        case HY_BL_ERROR_HANDLING_SOFT:
            if (appendToExisting) {
              _FString * existing = (_FString*) FetchObjectFromVariableByType(&_hyLastExecutionError, STRING);
              if (existing) {
                errMsg = *existing->theString & '\n' & errMsg;
              }
            }
            setParameter(_hyLastExecutionError, new _FString (errMsg, false), false);
            
            break;
        default: 
            WarnError (errMsg);
    }
}

//____________________________________________________________________________________
_String*    _ExecutionList::FetchFromStdinRedirect (void)
// grab a string from the front of the input queue
// complain if nothing is left
{
    if (!stdinRedirect) {
        WarnError ("No input buffer was given for a redirected standard input read.");
        return new _String;
    }
    long d = stdinRedirect->First();
    if (d<0) {
        WarnError ("Ran out of input in buffer during a redirected standard input read.");
        return new _String;
    }
    _String *sendBack = (_String*)stdinRedirect->GetXtra (d);
    sendBack->nInstances++;
    stdinRedirect->Delete ((*(_List*)stdinRedirect->dataList)(d),true);
    return sendBack;
}

//____________________________________________________________________________________

_String       _ExecutionList::GetFileName     (void)  {
    if (sourceFile.sLength) {
        return sourceFile;
    } else {
        if (pathNames.lLength)
            return *(_String*)pathNames.GetElement (-1);
    }
    return empty;
}
// doesn't do much
//____________________________________________________________________________________

_PMathObj       _ExecutionList::Execute     (void)      // run this execution list
{

    setParameter(_hyLastExecutionError, new _MathObject, false);
    
    _ExecutionList*      stashCEL = currentExecutionList;
    callPoints << currentCommand;
    executionStack       << this;

    _String             dd (GetPlatformDirectoryChar());

    _FString            bp  (baseDirectory, false),
                        lp  (libDirectory, false),
                        ds  (dd),
                        cfp (pathNames.lLength?*(_String*)pathNames(pathNames.lLength-1):empty),
                        * stashed = (_FString*)FetchObjectFromVariableByType (&pathToCurrentBF, STRING);

    setParameter        (platformDirectorySeparator, &ds);
    setParameter        (hyphyBaseDirectory, &bp);
    setParameter        (hyphyLibDirectory, &lp);

    if (stashed) {
        stashed = (_FString*)stashed->makeDynamic();
    }
    setParameter        (pathToCurrentBF,&cfp);

    DeleteObject        (result);
    result               = nil;
    currentExecutionList = this;
    currentCommand       = 0;

    terminateExecution  = false;
    skipWarningMessages = false;

    while (currentCommand<lLength) {
        if (doProfile == 1 && profileCounter) {
            long        instCounter = currentCommand;
            _Parameter  timeDiff    = 0.0;

            TimerDifferenceFunction (false);
            (((_ElementaryCommand**)lData)[currentCommand])->Execute(*this);
            timeDiff   = TimerDifferenceFunction(true);

            if (profileCounter) {
                profileCounter->theData[instCounter*2]   += timeDiff;
                profileCounter->theData[instCounter*2+1] += 1.0;
            }
        } else {
            (((_ElementaryCommand**)lData)[currentCommand])->Execute(*this);
        }

        if (terminateExecution) {
            break;
        }
    }
    currentCommand = callPoints.lData[callPoints.lLength-1];
    callPoints.Delete (callPoints.lLength-1);
    currentExecutionList = stashCEL;

    if (stashed) {
        setParameter        (pathToCurrentBF,stashed,false);
    }

    executionStack.Delete (executionStack.lLength-1);
    if (result == nil) {
        result = new _MathObject();
    }

    return result;
}

//____________________________________________________________________________________

long        _ExecutionList::ExecuteAndClean     (long g, _String* fName)        // run this execution list
{
    long    f = -1;
    Execute ();

    if (fName && !terminateExecution) {
        f = batchLanguageFunctionNames.Find (fName);
    }

    terminateExecution      = false;
    skipWarningMessages     = false;

    while (g<batchLanguageFunctionNames.lLength) {
        batchLanguageFunctionNames.Delete           (g);
        batchLanguageFunctionParameters.Delete      (g);
        batchLanguageFunctions.Delete               (g);
        batchLanguageFunctionClassification.Delete  (g);
        batchLanguageFunctionParameterLists.Delete  (g);
    }
    return f;
}

//____________________________________________________________________________________

bool        _ExecutionList::TryToMakeSimple     (void)
{
    _SimpleList     varList,
                    formulaeToConvert,
                    parseCodes;

    long            stackDepth  = 0;

    bool            status      = true;

    for (unsigned long k = 0; k<lLength && status; k++) {
        _ElementaryCommand * aStatement = (_ElementaryCommand*)(*this)(k);
        switch (aStatement->code) {
        case 0: {
            _String * formulaString = (_String*)aStatement->parameters(0);

            if (formulaString->sData[formulaString->sLength-1]!='}') {
                _Formula *f  = new _Formula,
                *f2 = new _Formula;

                checkPointer ((BaseRef)(f&&f2));

                _FormulaParsingContext fpc (nil, nameSpacePrefix);

                long          parseCode = Parse(f,*formulaString,fpc,f2);

                if (parseCode == HY_FORMULA_EXPRESSION || parseCode == HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT) {
                    if (f->AmISimple(stackDepth,varList)) {
                        aStatement->simpleParameters<<parseCode;
                        aStatement->simpleParameters<<(long)f;
                        aStatement->simpleParameters<<(long)f2;
                        aStatement->simpleParameters<<fpc.assignmentRefID();

                        formulaeToConvert << (long)f;

                        if (HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT) {
                            parseCodes        << fpc.assignmentRefID();
                        } else {
                            parseCodes        << -1;
                        }
                        break;
                    }
                }

                delete f;
                delete f2;
            }
            status = false;
            break;
        }

        case 4:
            parseCodes        << -1;
            if (aStatement->simpleParameters.lLength == 3 || aStatement->parameters.lLength) {
                if (aStatement->parameters.lLength) {
                    _Formula f;
                    _FormulaParsingContext fpc (nil, nameSpacePrefix);

                    long status = Parse (&f, *(_String*)aStatement->parameters(0), fpc);

                    if (status== HY_FORMULA_EXPRESSION) {
                        aStatement->simpleParameters<<long(f.makeDynamic());
                    }
                }


                _Formula *cf  = ((_Formula*)aStatement->simpleParameters(2));
                if (cf->AmISimple(stackDepth,varList)) {
                    formulaeToConvert << (long)cf;
                } else {
                    status = false;
                }
            }
            break;

        default:
            status = false;
        }
        if (status == false) {
            ReportWarning (_String ("Failed to compile an execution list: offending command was\n") & _String (((_String*)aStatement->toStr())));
        }
    }

    if (status) {
        cli = new _CELInternals;
        checkPointer (cli);
        checkPointer (cli->values = new _SimpleFormulaDatum[varList.lLength+1]);
        checkPointer (cli->stack  = new _SimpleFormulaDatum[stackDepth+1]);

        _SimpleList  avlData;
        _AVLListX    avlList (&avlData);

        for (unsigned long fi = 0; fi < formulaeToConvert.lLength; fi++) {
            ((_Formula*)formulaeToConvert(fi))->ConvertToSimple (varList);
        }

        for (unsigned long vi = 0; vi < varList.lLength; vi++) {
            avlList.Insert ((BaseRef)varList.lData[vi], vi);
        }

        for (unsigned long ri = 0; ri<parseCodes.lLength; ri++) {
            if (parseCodes.lData[ri] < 0) {
                cli->storeResults << -1;
            } else {
                cli->storeResults << avlList.GetXtra (avlList.Find ((BaseRef) parseCodes.lData[ri]));
            }
        }
        cli->varList.Duplicate(&varList);
    }

    return status;
}

//____________________________________________________________________________________

void        _ExecutionList::ExecuteSimple       (void)
{
    PopulateArraysForASimpleFormula (cli->varList, cli->values);
    Execute();

    for (long vi2 = 0; vi2 < cli->varList.lLength; vi2++) {
        _Variable * mv = LocateVar(cli->varList.lData[vi2]);
        if (mv->ObjectClass() == NUMBER) {
            mv->SetValue (new _Constant (cli->values[vi2].value),false);
        }
    }
}

//____________________________________________________________________________________

void        _ExecutionList::ResetFormulae       (void)      // run this execution list
{
    currentCommand = 0;
    while (currentCommand<lLength) {
        _ElementaryCommand* thisCommand = ((_ElementaryCommand**)lData)[currentCommand];
        if (thisCommand->code==0) {
            if (thisCommand->simpleParameters.lLength) {
                //printf ("[ResetFormulae] %s\n", thisCommand->sData);
                _Formula* f = (_Formula*)
                              thisCommand->simpleParameters.lData[1],
                              *f2 = (_Formula*)
                                    thisCommand->simpleParameters.lData[2] ;
                if (f) {
                    delete f;
                }
                if (f2) {
                    delete f2;
                }
                thisCommand->simpleParameters.Clear();
                long k = listOfCompiledFormulae.Find((long)thisCommand);
                if (k >= 0) {
                    listOfCompiledFormulae.Delete(k);
                    //printf ("[ResetFormulae:listOfCompiledFormulae %d]\n",k);
                    compiledFormulaeParameters.Delete(k);
                    //printf ("[ResetFormulae:compiledFormulaeParameters %d]\n",k);
                }
            }
        } else {
            if (thisCommand->code==4) {
                if (thisCommand->parameters.lLength && thisCommand->simpleParameters.lLength == 3) {
                    _Formula* f = (_Formula*)thisCommand->simpleParameters.lData[2];
                    if (f) {
                        delete f;
                    }
                    thisCommand->simpleParameters.Delete (2);
                }
            }
        }
        currentCommand++;
    }
}
//____________________________________________________________________________________

BaseRef  _ExecutionList::toStr (void)
{
    _String *result = new _String (1,true),
    step ("\n\nStep"),
    dot (".");

    for (unsigned long i=0; i<countitems(); i++) {
        (*result) << &step;
        _String lineNumber (i);
        (*result)<< &lineNumber;
        (*result)<< '.';
        result->AppendNewInstance ((_String*)(*this)(i)->toStr());
    }
    result->Finalize();
    return result;
}

//____________________________________________________________________________________

void     _ExecutionList::ResetNameSpace (void)
{
    DeleteObject (nameSpacePrefix);
    nameSpacePrefix = nil;
}

//____________________________________________________________________________________

void     _ExecutionList::SetNameSpace (_String nID)
{
    ResetNameSpace ();
    nameSpacePrefix = new _VariableContainer(nID);
    checkPointer(nameSpacePrefix);
}

//____________________________________________________________________________________

_String*     _ExecutionList::GetNameSpace ()
{
    if (nameSpacePrefix) {
        return nameSpacePrefix->GetName();
    }
    return nil;
}

//____________________________________________________________________________________

_String  _ExecutionList::AddNameSpaceToID (_String& theID, _String * extra)
{
    _String check_dereferences,
            name_space;
            
    if (extra && extra->sLength) {
        if (nameSpacePrefix) {
            name_space = (*nameSpacePrefix->GetName())&'.'& *extra;
        } else {
            name_space = *extra;
        }
    } else {
        if (nameSpacePrefix) {
            name_space = (*nameSpacePrefix->GetName());        
        }
    }
            
    return AppendContainerName (theID, &name_space);
}

//____________________________________________________________________________________

_String  _ExecutionList::TrimNameSpaceFromID (_String& theID)
{
    if (nameSpacePrefix) {
        if (theID.startswith(*nameSpacePrefix->GetName())) {
            return theID.Cut(nameSpacePrefix->GetName()->sLength+1,-1);
        }
    }
    return theID;
}

//____________________________________________________________________________________

bool        _ExecutionList::BuildList   (_String& s, _SimpleList* bc, bool processed, bool empty_is_success)
{
    if (terminateExecution) {
        return false;
    }

    char * savePointer = s.sData;
    
    _SimpleList          triePath;

    while (s.Length()) { // repeat while there is stuff left in the buffer
        _String currentLine (_ElementaryCommand::FindNextCommand (s,true));

        if (currentLine.getChar(0)=='}') {
            currentLine.Trim(1,-1);
        }

        if (!currentLine.Length()) {
            continue;
        }
        
        triePath.Clear(false);
        long prefixTreeCode = _HY_ValidHBLExpressions.Find (currentLine, &triePath, true);
        
        _List *pieces = nil;
        _HBLCommandExtras *commandExtraInfo = nil;
        
        if (prefixTreeCode != HY_TRIE_NOTFOUND) {
            prefixTreeCode = _HY_ValidHBLExpressions.GetValue(prefixTreeCode);
            long commandExtra = _HY_HBLCommandHelper.FindLong (prefixTreeCode);
            if (commandExtra >= 0) { // pre-trim all strings as needed
                commandExtraInfo = (_HBLCommandExtras*)_HY_HBLCommandHelper.GetXtra (commandExtra);
                if (commandExtraInfo->extract_conditions.lLength > 0) {
                    pieces = new _List;
                    long upto = _ElementaryCommand::ExtractConditions (currentLine, commandExtraInfo->cut_string,*pieces,commandExtraInfo->extract_condition_separator),
                         condition_index_match = commandExtraInfo->extract_conditions.Find(pieces->lLength);
                    if (condition_index_match < 0) {
                        // try to see if the command accepts a variable number of arguments (at least X)
                       _String parseFail;
                       if (commandExtraInfo->extract_conditions.lLength == 1 && commandExtraInfo->extract_conditions.lData[0] < 0) {
                            if (pieces->lLength < -commandExtraInfo->extract_conditions.lData[0]) {
                                 parseFail = _String("Incorrect number of arguments (") & (long) pieces->lLength & ") supplied: expected at least " & _String (-commandExtraInfo->extract_conditions.lData[0]) & ", while processing '"& currentLine.Cut (0, upto) & "'. ";
                             }
                        } else {
                            parseFail = _String("Incorrect number of arguments (") & (long) pieces->lLength & ") supplied: expected one of " & _String ((_String*)commandExtraInfo->extract_conditions.toStr()) & ", while processing '"& currentLine.Cut (0, upto) & "'. ";
                        }
                        if (parseFail.sLength) {
                            if (currentExecutionList) {
                                currentExecutionList->ReportAnExecutionError(parseFail, false, true);
                            } else {
                                acknError (parseFail);
                            }
                            DeleteObject (pieces);
                            return false;  
                        }
                    }
                    if (commandExtraInfo->do_trim) {
                        currentLine.Trim (upto, -1);
                    }
                }
            }
        }
        
        bool handled = false;
               
        switch (prefixTreeCode) {
            case HY_HBL_COMMAND_FOR:
                _ElementaryCommand::BuildFor (currentLine, *this, *pieces);
                handled = true;
                break;
            case HY_HBL_COMMAND_WHILE:
                _ElementaryCommand::BuildWhile (currentLine, *this, *pieces);
                handled = true;
                break;
            case HY_HBL_COMMAND_BREAK:
            case HY_HBL_COMMAND_CONTINUE:
                if (bc) {
                    AppendNewInstance(new _ElementaryCommand);
                    (*bc) << ((prefixTreeCode == HY_HBL_COMMAND_BREAK) ? (countitems()-1) : (-(long)countitems()+1));
                } else {
                    WarnError (currentLine & " only makes sense in the context of a loop.");
                    return false;
                }
                handled = true;
                break;
            case HY_HBL_COMMAND_SET_DIALOG_PROMPT:
            case HY_HBL_COMMAND_HARVEST_FREQUENCIES:
            case HY_HBL_COMMAND_OPTIMIZE:
            case HY_HBL_COMMAND_COVARIANCE_MATRIX:
            case HY_HBL_COMMAND_LFCOMPUTE:
            case HY_HBL_COMMAND_SELECT_TEMPLATE_MODEL:
            case HY_HBL_COMMAND_USE_MODEL:
            case HY_HBL_COMMAND_SET_PARAMETER:
            case HY_HBL_COMMAND_ASSERT:
            case HY_HBL_COMMAND_REQUIRE_VERSION:
            case HY_HBL_COMMAND_DELETE_OBJECT:
            case HY_HBL_COMMAND_CLEAR_CONSTRAINTS:
            case HY_HBL_COMMAND_MOLECULAR_CLOCK:
            case HY_HBL_COMMAND_GET_URL:
            case HY_HBL_COMMAND_GET_STRING:
            case HY_HBL_COMMAND_EXPORT:
            case HY_HBL_COMMAND_DIFFERENTIATE:
            case HY_HBL_COMMAND_FPRINTF:
                _ElementaryCommand::ExtractValidateAddHBLCommand (currentLine, prefixTreeCode, pieces, commandExtraInfo, *this);
                handled = true;
                break;
                
        }
        
        if (handled)
            DeleteObject (pieces);
        
        // 20111212: this horrendous switch statement should be replaced with a 
        // prefix tree lookup 

        if (!handled) {
            if (currentLine.startswith (blFunction)||currentLine.startswith (blFFunction)||currentLine.startswith (blLFunction)) { // function declaration
                _ElementaryCommand::ConstructFunction (currentLine, *this);
            } else if (currentLine.startswith (blReturn) || currentLine.startswith (blReturn2)) { // function return statement
                _ElementaryCommand::ConstructReturn (currentLine, *this);
            } else if (currentLine.startswith (blIf)) { // if-then-else statement
                _ElementaryCommand::BuildIfThenElse (currentLine, *this, bc);
            } else if (currentLine.startswith (blElse)) { // else clause of an if-then-else statement
                if (lastif.countitems()) {
                    long    temp = countitems(),
                            lc   = lastif.countitems(),
                            lif  = lastif.lData[lc-1];

                    _ElementaryCommand      * stuff = new _ElementaryCommand ();
                    stuff->MakeJumpCommand  (nil,0,0,*this);
                    AppendNewInstance       (stuff);
                    currentLine.Trim        (4,-1);

                    long  index         = currentLine.Length()-1,
                          scopeIn     = 0;

                    while (currentLine.sData[scopeIn]=='{' && currentLine.sData[index]=='}') {
                        scopeIn++;
                        index--;
                    }

                    if (scopeIn) {
                        currentLine.Trim (scopeIn,index);
                    }

                    BuildList (currentLine,bc,true);

                    if (lif<0 || lif>=lLength) {
                        WarnError ("'else' w/o an if to latch on to...");
                        return false;
                    }

                    ((_ElementaryCommand*)((*this)(lif)))->MakeJumpCommand(nil,-1,temp+1,*this);
                    ((_ElementaryCommand*)(*this)(temp))->simpleParameters[0]=countitems();

                    while (lastif.countitems()>=lc) {
                        lastif.Delete(lastif.countitems()-1);
                    }
                } else {
                    WarnError ("'else' w/o an if to latch on to...");
                    return false;
                }

            } else if (currentLine.startswith (blDo)) { // do {} while statement
                _ElementaryCommand::BuildDoWhile (currentLine, *this);
            }  else if (currentLine.startswith (blInclude)) { // #include
                _ElementaryCommand::ProcessInclude (currentLine, *this);
            } else if (currentLine.startswith (blDataSet)) { // data set definition
                _ElementaryCommand::ConstructDataSet (currentLine, *this);
            } else if (currentLine.startswith (blDataSetFilter)) { // data set filter definition
                _ElementaryCommand::ConstructDataSetFilter (currentLine, *this);
            } else if (currentLine.startswith (blConstructCM)) { // construct category assignments matrix
                _ElementaryCommand::ConstructCategoryMatrix (currentLine, *this);
            } else if (currentLine.startswith (blTree) || currentLine.startswith (blTopology)) { // tree definition
                _ElementaryCommand::ConstructTree (currentLine, *this);
            } else if (currentLine.startswith (blLF) || currentLine.startswith (blLF3)) { // LF definition
                _ElementaryCommand::ConstructLF (currentLine, *this);
            } else if (currentLine.startswith (blfscanf) || currentLine.startswith (blsscanf)) { // fscanf call
                _ElementaryCommand::ConstructFscanf (currentLine, *this);
            } else if (currentLine.startswith (blReplicate)) { // replicate constraint statement
                _ElementaryCommand::ConstructReplicateConstraint (currentLine, *this);
            } else if (currentLine.startswith (blCategory)) { // category variable declaration
                _ElementaryCommand::ConstructCategory (currentLine, *this);
            } else if (currentLine.startswith (blGetNeutralNull)) { // select a template model
                _ElementaryCommand::ConstructGetNeutralNull (currentLine, *this);
            } else if (currentLine.startswith (blModel)) { // Model declaration
                _ElementaryCommand::ConstructModel (currentLine, *this);
            } else if (currentLine.startswith (blChoiceList)) { // choice list
                _ElementaryCommand::ConstructChoiceList (currentLine, *this);
            } else if (currentLine.startswith (blOpenDataPanel)) { // open data panel window
                _ElementaryCommand::ConstructOpenDataPanel (currentLine, *this);
            } else if (currentLine.startswith (blGetInformation)) { // get information
                _ElementaryCommand::ConstructGetInformation (currentLine, *this);
            } else if (currentLine.startswith (blExecuteCommands) || currentLine.startswith (blExecuteAFile) || currentLine.startswith (blLoadFunctionLibrary))
                // execute commands
            {
                _ElementaryCommand::ConstructExecuteCommands (currentLine, *this);
            } else if (currentLine.startswith (blOpenWindow)) { // execute commands
                _ElementaryCommand::ConstructOpenWindow (currentLine, *this);
            } else if (currentLine.startswith (blSpawnLF)) { // execute commands
                _ElementaryCommand::ConstructSpawnLF (currentLine, *this);
            } else if (currentLine.startswith (blFindRoot)||currentLine.startswith (blIntegrate))
                // find a root of an expression in an interval
                // or an integral
            {
                _ElementaryCommand::ConstructFindRoot (currentLine, *this);
            } else if (currentLine.startswith (blMPISend)) { // MPI Send
                _ElementaryCommand::ConstructMPISend (currentLine, *this);
            } else if (currentLine.startswith (blMPIReceive)) { // MPI Receive
                _ElementaryCommand::ConstructMPIReceive (currentLine, *this);
            } else if (currentLine.startswith (blGetDataInfo)) { // Get Data Info
                _ElementaryCommand::ConstructGetDataInfo (currentLine, *this);
            } else if (currentLine.startswith (blStateCounter)) { // Get Data Info
                _ElementaryCommand::ConstructStateCounter (currentLine, *this);
            } else if (currentLine.startswith (blDoSQL)) { // Do SQL
                _ElementaryCommand::ConstructDoSQL (currentLine, *this);
            } else if (currentLine.startswith (blAlignSequences)) { // Do AlignSequences
                _ElementaryCommand::ConstructAlignSequences (currentLine, *this);
            } else if (currentLine.startswith (blHBLProfile)) { // #profile
                _ElementaryCommand::ConstructProfileStatement (currentLine, *this);
            } else if (currentLine.startswith (blSCFG)) { // SCFG definition
                _ElementaryCommand::ConstructSCFG (currentLine, *this);
            } else if (currentLine.startswith (blNN)) { // Neural Net definition
                _ElementaryCommand::ConstructNN (currentLine, *this);
            } else if (currentLine.startswith (blBGM)) {    // Bayesian Graphical Model definition
                _ElementaryCommand::ConstructBGM (currentLine, *this);
            } 
            // plain ol' formula - parse it as such!
            else {
                _String checker (currentLine);
                if (_ElementaryCommand::FindNextCommand (checker).Length()==currentLine.Length()) {
                    if (currentLine.Length()>1)
                        while (currentLine[currentLine.Length()-1]==';') {
                            currentLine.Trim (0,currentLine.Length()-2);
                        }
                    else {
                        continue;
                    }
                    _ElementaryCommand* oddCommand = new _ElementaryCommand(currentLine);
                    oddCommand->code = 0;
                    oddCommand->parameters&&(&currentLine);
                    AppendNewInstance (oddCommand);
                } else {
                    while (currentLine.Length()) {
                        _String part (_ElementaryCommand::FindNextCommand (currentLine));
                        BuildList (part,bc,processed);
                    }
                }
            }
            
            /*if (currentLine.sLength > 1 || currentLine.sLength == 1 && currentLine.getChar(0) != ';'){
                WarnError (_String ("Missing semicolon before ") & currentLine);
                return false;
            }*/
        }
    }
    s.sData = savePointer;
    s.DuplicateErasing (&empty);
    return empty_is_success || countitems();
}


