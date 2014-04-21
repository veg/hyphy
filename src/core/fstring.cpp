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

#include "defines.h"
#include "fstring.h"
#include "constant.h"
#include "matrix.h"
#include "calcnode.h"

extern long lastMatrixDeclared;
extern _AVLListX _HY_GetStringGlobalTypes;

extern _List likeFuncList,
             batchLanguageFunctionNames,
             dataSetFilterList,
             dataSetList,
             scfgList;

extern _SimpleList modelMatrixIndices;
extern _String lastModelParameterList;

_String internalRerootTreeID("_INTERNAL_REROOT_TREE_");
//__________________________________________________________________________________
_FString::_FString (void)
{
    theString = new _String;
}

//__________________________________________________________________________________
_FString::~_FString ()
{
    if (nInstances<=1) {
        DeleteObject (theString);
    } else {
        nInstances--;
    }

}

//__________________________________________________________________________________
_FString::_FString (_String* data)
{
    theString = data;
}

//__________________________________________________________________________________
_FString::_FString (long inData)
{
    checkPointer (theString = new _String (inData));
}

//__________________________________________________________________________________
_FString::_FString (_String& data, bool meta)
{
    if (meta) {
        unsigned long ssi = _String::storageIncrement;

        if (data.sLength>ssi) {
            _String::storageIncrement = data.sLength;
        }

        theString = new _String (data.sLength,true);
        for (unsigned long k=0; k<data.sLength; k++) {
            char c = data.sData[k];
            if (c=='\\') {
                if (k<data.sLength-1) {
                    c = data.sData[k+1];
                    switch (c) {
                    case 'n':
                        (*theString)<<'\n';
                        break;
                    case 't':
                        (*theString)<<'\t';
                        break;
                    case '\\':
                        (*theString)<<'\\';
                        break;
                    default:
                        (*theString)<<'\\';
                        (*theString)<<c;
                    }
                    k++;
                    continue;
                }
            }
            (*theString)<<c;
        }
        _String::storageIncrement = ssi;
        theString->Finalize();
    } else {
        theString = new _String (data);
    }
}

//__________________________________________________________________________________
BaseRef _FString::makeDynamic (void)
{
    //_FString* res = new _FString;
    //res->theString->Duplicate(theString);
    return new _FString(*theString, false);
}

//__________________________________________________________________________________
void _FString::Duplicate (BaseRef o)
{
    DeleteObject (theString);
    _FString* m = (_FString*)o;
    theString = (_String*)m->theString->makeDynamic();
}

//__________________________________________________________________________________
_PMathObj _FString::Add (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        _String * res    = new _String (theString->sLength+theStr->theString->sLength,true);
        if (res) {
            (*res) << theString;
            (*res) << theStr->theString;
            res->Finalize();
            return new _FString (res);
        } else {
            checkPointer (nil);
        }
    }
    _String* convStr = (_String*)p->toStr();
    _String res (*theString& *((_String*)convStr));
    DeleteObject (convStr);
    return new _FString (res, false);
}

//__________________________________________________________________________________
long _FString::AddOn (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        (*theString) << theStr->theString;
        return theStr->theString->sLength;
    } else if (p->ObjectClass()==NUMBER) {
        long s = p->Value();
        if (s) {
            delete theString;
            checkPointer (theString = new _String (s, true));
            return s;
        } else {
            theString->Finalize();
        }
    } else {
        WarnError ("Invalid 2nd argument in call to string*number");
    }
    return 0;
}

//__________________________________________________________________________________
_PMathObj _FString::AreEqual (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        bool     equal = theString->Equal(theStr->theString);
        return new _Constant ((_Parameter)equal);
    } else {
        /*_String* convStr = (_String*)p->toStr();
         bool     equal = theString->Equal(convStr);
         DeleteObject (convStr);
         return new _Constant ((_Parameter)equal);*/
         return new HY_CONSTANT_FALSE;
    }
}

//__________________________________________________________________________________
_PMathObj _FString::AreEqualCIS (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        _String   t1 (*theString),
                  t2 (*theStr->theString);
        t1.UpCase();
        t2.UpCase();
        bool     equal = t1.Equal(&t2);
        return new _Constant ((_Parameter)equal);
    } else {
        return AreEqual (p);
    }
}

//__________________________________________________________________________________
_PMathObj _FString::Join (_PMathObj p)
{
    _List theStrings;

    if (p->ObjectClass()==MATRIX) {
        ((_Matrix*)(p->Compute()))->FillInList (theStrings,true);
    } else if (p->ObjectClass()==ASSOCIATIVE_LIST) {
        ((_AssociativeList*)(p->Compute()))->FillInList (theStrings);
    }

    return new _FString((_String*)theStrings.Join(theString));
}

//__________________________________________________________________________________
_PMathObj _FString::EqualAmb (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        bool     equal = theString->EqualWithWildChar(theStr->theString);
        return new _Constant ((_Parameter)equal);
    } else {
        _String  convStr      ((_String*)p->toStr());
        return   new _Constant(theString->EqualWithWildChar(&convStr));
    }
}

//__________________________________________________________________________________

_PMathObj _FString::EqualRegExp (_PMathObj p, bool matchAll)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        _SimpleList matches;

        if (matchAll) {
            int errNo = 0;
            Ptr regex = PrepRegExp (theStr->theString, errNo, true);
            if (regex) {
                theString->RegExpMatchAll(regex, matches);
                FlushRegExp (regex);
            } else {
                WarnError (GetRegExpError (errNo));
            }
        } else {
            theString->RegExpMatchOnce(theStr->theString, matches, true, true);
        }

        if (matches.lLength == 0) {
            matches << -1;
            matches << -1;
        }
        _Matrix * res = new _Matrix (matches);
        res->Transpose();
        return res;
    } else {
        WarnError ("Invalid 2nd argument in call to string$reg.exp.");
        return new _Matrix (2,1,false,true);
    }
}

//__________________________________________________________________________________
_PMathObj _FString::ReplaceReqExp (_PMathObj p)
{
    if (theString && theString->sLength) {
        if (p->ObjectClass()==MATRIX) {
            _Matrix * m = (_Matrix*)p;

            if (m->IsAStringMatrix() && m->GetHDim() * m->GetVDim() >= 2) {
                _FString* theStr  = (_FString*)m->GetFormula(0,0)->Compute(),
                          * repWith = (_FString*)m->GetFormula(1,-1)->Compute();

                _SimpleList matches;

                int errNo = 0;
                Ptr regex = PrepRegExp (theStr->theString, errNo, true);

                if (!regex) {
                    WarnError (GetRegExpError (errNo));
                    return new _FString (empty);
                }

                theString->RegExpMatchAll(regex, matches);
                _FString * res;
                if (matches.lLength) {
                    _String * newString = new _String (theString->sLength+1,true);
                    checkPointer (newString);
                    long      idx  = matches.lData[0],
                              midx = 0;

                    for (long k=0; k<theString->sLength;) {
                        if (k==idx) {
                            (*newString) << repWith->theString;
                            k = matches.lData[midx+1]+1;
                            midx += 2;
                            if (midx == matches.lLength) {
                                idx = -1;
                            } else {
                                idx = matches.lData[midx];
                            }
                        } else {
                            (*newString) << theString->sData[k++];
                        }
                    }
                    newString->Finalize();
                    res = new _FString (newString);
                } else {
                    res = new _FString (*theString,false);
                }

                FlushRegExp (regex);
                return res;
            }
        }

        WarnError ("Invalid 2nd argument in call to string^{{pattern,replacement}}");
    }
    return new _FString (empty,false);
}

//__________________________________________________________________________________
_PMathObj _FString::NotEqual (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        bool     equal = theString->Equal(theStr->theString);
        return new _Constant ((_Parameter)!equal);
    } else {
        //_String* convStr = (_String*)p->toStr();
        //bool     equal = theString->Equal(convStr);
        //DeleteObject (convStr);
        //return new _Constant ((_Parameter)!equal);
        return new HY_CONSTANT_TRUE;
    }
}

//__________________________________________________________________________________
_PMathObj _FString::Less (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        bool     equal = theString->Less(theStr->theString);
        return new _Constant ((_Parameter)equal);
    } else {
        _String* convStr = (_String*)p->toStr();
        bool     equal = theString->Less(convStr);
        DeleteObject (convStr);
        return new _Constant ((_Parameter)equal);
    }
}

//__________________________________________________________________________________
_PMathObj _FString::LessEq (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        bool     equal = theString->Less(theStr->theString)||theString->Equal(theStr->theString);
        return new _Constant ((_Parameter)equal);
    } else {
        _String* convStr = (_String*)p->toStr();
        bool     equal = theString->Less(convStr)||theString->Equal(convStr);
        DeleteObject (convStr);
        return new _Constant ((_Parameter)equal);
    }
}

//__________________________________________________________________________________
_PMathObj _FString::Differentiate (_PMathObj p)
{
    _Formula F;

    _String  *X,
             *DFDX = nil;

    bool     deleteX = false;

    if (p->ObjectClass()==STRING) {
        X = ((_FString*)p)->theString;
    } else {
        deleteX = true;
        X = (_String*)p->toStr();
    }

    _String copyMe (*theString);
    _FormulaParsingContext fpc;
    if (Parse (&F,copyMe, fpc, nil) == HY_FORMULA_EXPRESSION) {
        _Formula *DF = F.Differentiate (*X,true);
        if (DF) {
            DFDX = (_String*)DF->toStr();
        }

    }

    if (deleteX) {
        DeleteObject (X);
    }

    return new _FString (DFDX?DFDX:new _String());

}

//__________________________________________________________________________________
_PMathObj _FString::GreaterEq (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        bool     equal = theString->Greater(theStr->theString)||theString->Equal(theStr->theString);
        return new _Constant ((_Parameter)equal);
    } else {
        _String* convStr = (_String*)p->toStr();
        bool     equal = theString->Greater(convStr)||theString->Equal(convStr);
        DeleteObject (convStr);
        return new _Constant ((_Parameter)equal);
    }
}

//__________________________________________________________________________________
_PMathObj _FString::Greater (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        bool     equal = theString->Greater(theStr->theString);
        return new _Constant ((_Parameter)equal);
    } else {
        _String* convStr = (_String*)p->toStr();
        bool     equal = theString->Greater(convStr);
        DeleteObject (convStr);
        return new _Constant ((_Parameter)equal);
    }
}

//__________________________________________________________________________________
BaseRef  _FString::toStr()
{
    theString->nInstances++;
    return theString;
}

//__________________________________________________________________________________
_PMathObj _FString::RerootTree (void)
{
    long     stashedModelID = lastMatrixDeclared,
             totalNodeCount = 0;

    lastMatrixDeclared      = HY_NO_MODEL;
    /* unset current model; do not want the internal tree to have an attached model */

    _TheTree    rTree (internalRerootTreeID,*theString);

    if (rTree.IsDegenerate()) { // no need to reroot two-sequence trees
        lastMatrixDeclared = stashedModelID;
        DeleteVariable  (internalRerootTreeID);
        return new _FString (*theString, false);
    }

    if (terminateExecution) {
        lastMatrixDeclared = stashedModelID;
        DeleteVariable  (internalRerootTreeID);
        return new _FString;
    }

    _CalcNode   *iterator = rTree.DepthWiseTraversal (true),
                 *rerootAt;

    node<long>  *cNode;

    _GrowingVector  valueCache;

    while       (iterator)
        // count the number of descendants of a given node, store as numeric value of the CalcNode
    {
        cNode    = &rTree.GetCurrentNode();
        valueCache.Store(iterator->Value());
        if (long myNodeCount = cNode->get_num_nodes()) {
            _Parameter tNodeCount = 0.0;

            for (long k = 1; k <= myNodeCount; k++) {
                tNodeCount += ((_CalcNode*)LocateVar(cNode->go_down(k)->in_object))->Value();
            }

            iterator->SetNumericValue(tNodeCount+1.0);
        } else {
            iterator->SetNumericValue(1.0);
        }

        iterator = rTree.DepthWiseTraversal (false);
        totalNodeCount ++;
    }

    iterator = rTree.DepthWiseTraversal (true);

    long        maxMin = 0;
    _Parameter  bRatio  = 0.0;

    while       (iterator) {
        _Parameter      nodeMin   = totalNodeCount-iterator->Value(),
                        thisRatio = nodeMin/(_Parameter)iterator->Value();

        if (thisRatio>1.0) {
            thisRatio = 1./thisRatio;
        }

        cNode    = &rTree.GetCurrentNode();
        if (cNode->get_num_nodes()) {
            for (long k = cNode->get_num_nodes(); k; k--) {
                long tt = ((_CalcNode*)LocateVar(cNode->go_down(k)->in_object))->Value();
                if (tt<nodeMin) {
                    nodeMin = tt;
                }
            }
        } else {
            nodeMin = 1;
        }

        if ((nodeMin>maxMin)||((nodeMin==maxMin)&&(thisRatio>bRatio))) {
            bRatio = thisRatio;
            maxMin = nodeMin;
            rerootAt = iterator;
            if (!cNode->get_parent()) {
                rerootAt = nil;
            }
        }
        iterator = rTree.DepthWiseTraversal (false);
    }

    iterator        = rTree.DepthWiseTraversal (true);
    totalNodeCount  = 0;
    while       (iterator)
        // restore branch lengths
    {
        iterator->SetNumericValue(valueCache.theData[totalNodeCount]);
        iterator = rTree.DepthWiseTraversal (false);
        totalNodeCount ++;
    }

    _FString* res;
    if (rerootAt) {
        _String stringCopy = *rerootAt->GetName();
        stringCopy.Trim (stringCopy.FindBackwards ('.',0,-1)+1,-1);
        _FString    rAt  (stringCopy);
        res = (_FString*)rTree.RerootTree (&rAt);
    } else {
        res = new _FString (*theString, false);
    }

    DeleteVariable  (internalRerootTreeID);
    lastMatrixDeclared = stashedModelID;

    return      res;
}

//__________________________________________________________________________________

_PMathObj _FString::Evaluate (_hyExecutionContext* context)
{
    if (theString && theString->sLength) {
        _String     s (*theString);
        _Formula    evaluator (s, (_VariableContainer*)context->GetContext());
        _PMathObj   evalTo = evaluator.Compute(0,(_VariableContainer*)context->GetContext());

        if (evalTo && !terminateExecution) {
            evalTo->AddAReference();
            return evalTo;
        }
    }
    return new _Constant (.0);
}

//__________________________________________________________________________________

_PMathObj _FString::Dereference(bool ignore_context, _hyExecutionContext* context, bool return_variable_ref) {
    _String referencedVariable = *theString;
    if (!ignore_context && context) {
        referencedVariable = AppendContainerName(referencedVariable, (_VariableContainer *)context->GetContext());
    }
    if (return_variable_ref) {
        return FetchVar (LocateVarByName(referencedVariable));
    }
    _PMathObj result = FetchObjectFromVariableByType(&referencedVariable, HY_ANY_OBJECT); 
    //printf ("\n\nDereferencing %s in this context %x\n\n", referencedVariable.sData, context);
    if (!result) {
        _String errM = _String("Failed to dereference '") & referencedVariable & "'";
        if (context) {
            context->ReportError(errM);
        } else {
            WarnError (errM);
        }
        result = new _FString;
    } else {
        result->AddAReference();
    }
    return result;
}

//__________________________________________________________________________________


_PMathObj _FString::Execute (long opCode, _PMathObj p, _PMathObj p2, _hyExecutionContext* context)   // execute this operation with the second arg if necessary
{
    switch (opCode) {
    case HY_OP_CODE_NOT: // !
        return FileExists();
    case HY_OP_CODE_NEQ: // !=
        return NotEqual(p);
        break;
    case HY_OP_CODE_IDIV: // $ match regexp
        return EqualRegExp(p);
        break;
    case HY_OP_CODE_MOD: // % equal case insenstive
        return AreEqualCIS(p);
        break;
    case HY_OP_CODE_AND: { // && upcase or lowercase
        _Parameter pVal = 0.0;
        if (p->ObjectClass () == NUMBER) {
            pVal = p->Value();
        }

        if (pVal < 0.0) {
            return (_PMathObj)makeDynamic();
        } else {
            _String * t = nil;


            if (CheckEqual(pVal,2.0) || CheckEqual(pVal,3.0) || CheckEqual(pVal,4.0) || CheckEqual(pVal,5.0) || CheckEqual(pVal,6.0)) {
                checkPointer (t = new _String (theString->sLength+1,true));
                t->EscapeAndAppend (*theString, CheckEqual(pVal,3.0) + 2*CheckEqual(pVal,4.0) + 4*CheckEqual(pVal,5.0) + 5*CheckEqual(pVal,6.0));
                t->Finalize();
            } else {
                t = new _String (*theString);
                checkPointer (t);
                if (CheckEqual(pVal,1.0)) {
                    t->UpCase();
                } else {
                    t->LoCase();
                }
            }

            return new _FString (t);
        }
    }
    break;
    case HY_OP_CODE_MUL: // *
        if (p) {
         // NOT a dereference
            if (p->ObjectClass() == MATRIX) {
                return      MapStringToVector (p);
            } else {
                return new _Constant(AddOn(p));
            }
        } else {
            return Dereference(false, context);
        }
        break;
    case HY_OP_CODE_ADD: // +
        if (p) {
            return Add(p);
        } else {
            return Sum();
        }
        break;
    case HY_OP_CODE_DIV: // /
        return EqualAmb(p);
        break;
    case HY_OP_CODE_LESS: // <
        return Less(p);
        break;
    case HY_OP_CODE_LEQ: // <=
        return LessEq(p);
        break;
    case HY_OP_CODE_EQ: // ==
        return AreEqual(p);
        break;
    case HY_OP_CODE_GREATER: // >
        return Greater(p);
        break;
    case HY_OP_CODE_GEQ: // >=
        return GreaterEq(p);
        break;
    case HY_OP_CODE_ABS: // Abs
        return new _Constant (theString->sLength);
        break;
    case HY_OP_CODE_DIFF: // Differentiate
        return Differentiate(p);
        break;
    case HY_OP_CODE_EVAL: // Eval
        return Evaluate(context);
        break;
    case HY_OP_CODE_EXP: // Exp
        return new _Constant (theString->LempelZivProductionHistory(nil));
        break;
    case HY_OP_CODE_FORMAT: { // Format
        _String cpyString (*theString);
        _Formula f (cpyString);
        _PMathObj fv = f.Compute();
        if (fv && fv->ObjectClass () == NUMBER) {
            return ((_Constant*)fv)->FormatNumberString (p,p2);
        } else {
            ReportWarning (_String("Failed to evaluate ")& *theString & " to a number in call to Format (string...)");
            return new _FString();
        }
    }
    break;
    case HY_OP_CODE_INVERSE: { // Inverse
        _FString * res = new _FString (*theString, false);
        checkPointer (res);
        for (long i1 = 0, i2 = theString->sLength-1; i1<theString->sLength; i1++, i2--) {
            res->theString->sData[i1] = theString->sData[i2];
        }

        return res;
    }
    break;
    case HY_OP_CODE_JOIN: // Inverse
        return Join (p);

    case HY_OP_CODE_LOG: // Log - check sum
        return new _Constant (theString->Adler32());
    case HY_OP_CODE_MACCESS: // MAccess
        return CharAccess(p,p2);
        break;
    case HY_OP_CODE_REROOTTREE: // RerootTree
        return RerootTree ();
        break;
    case HY_OP_CODE_ROWS: // Count Objects of given type
        return CountGlobalObjects();
        break;
    case HY_OP_CODE_TYPE: // Type
        return Type();
        break;
    case HY_OP_CODE_POWER: // Replace (^)
        if (p)
            return ReplaceReqExp (p);
        return Dereference(true, context);
        break;
    case HY_OP_CODE_OR: // Match all instances of the reg.ex (||)
        return EqualRegExp (p, true);
        break;
    }

    WarnNotDefined (this, opCode,context);
    return new _FString;

}

//__________________________________________________________________________________
_PMathObj   _FString::MapStringToVector (_PMathObj p)
{
    if (theString->sLength && p->ObjectClass () == MATRIX) {
        _Matrix         * factoringMatrix = (_Matrix *)p;

        if (factoringMatrix->IsAVector () && factoringMatrix->IsAStringMatrix()) {
            long            mapper [255],
                            keys    = factoringMatrix->GetHDim() * factoringMatrix->GetVDim(),
                            byRows   = factoringMatrix->IsAVector (HY_MATRIX_COLUMN_VECTOR);

            for (long c = 0; c < 255; c++) {
                mapper[c] = -1;
            }

            for (long r = 0; r < keys; r++) {
                _FString* aKey = (_FString*)factoringMatrix->GetFormula(byRows?r:0,byRows?0:r)->Compute();
                if (aKey->theString->sLength == 1) {
                    unsigned char thisChar = aKey->theString->sData[0];
                    if (mapper[thisChar] < 0) {
                        mapper[thisChar] = r;
                    }
                }
            }

            _SimpleList mapped;
            for (long s = 0; s < theString->sLength; s++) {
                mapped << mapper[(unsigned char)theString->sData[s]];
            }

            return new _Matrix (mapped);
        }
    }

    return new _Matrix;
}

//__________________________________________________________________________________
_PMathObj   _FString::CharAccess (_PMathObj p,_PMathObj p2)
{
    unsigned long index = p->Value();
    _String res;

    if (p2) {
        unsigned long index2 = p2->Value();
        res = theString->Cut (index,index2);
    } else if (index<theString->sLength) {
        res = theString->sData[index];
    }

    return new _FString (res);
}
//__________________________________________________________________________________
_PMathObj   _FString::FileExists (void)
{
    _Constant  * retValue = new _Constant (0.0);
    if (theString) {
        _String cpy (*theString);
        cpy.ProcessFileName();
        FILE * test = doFileOpen (cpy, "rb");
        if (test) {
            retValue->SetValue (1.0);
            fclose (test);
        }
    }
    return retValue;
}

//__________________________________________________________________________________
_PMathObj   _FString::CountGlobalObjects (void)
{
    _Parameter res = 0.0;

    long      standardType = _HY_GetStringGlobalTypes.Find(theString);
    if (standardType >=0 ) {
        standardType = _HY_GetStringGlobalTypes.GetXtra (standardType);
    }

    switch (standardType) {
    case HY_BL_LIKELIHOOD_FUNCTION:
        return new _Constant (likeFuncList.lLength);
    case HY_BL_DATASET:
        return new _Constant (dataSetList.lLength);
    case HY_BL_DATASET_FILTER:
        return new _Constant (dataSetFilterList.lLength);
    case HY_BL_HBL_FUNCTION:
        return new _Constant (batchLanguageFunctionNames.lLength);
    case HY_BL_TREE: {
        _SimpleList tc;
        long        si,
                    vi = variableNames.Traverser (tc,si,variableNames.GetRoot());

        for (; vi >= 0; vi = variableNames.Traverser (tc,si))
            if (((_Variable*)FetchVar(vi))->ObjectClass () == TREE) {
                res += 1.;
            }

        break;
    }

    case HY_BL_SCFG:
        return new _Constant (scfgList.lLength);
    case HY_BL_VARIABLE:
        return new _Constant (variableNames.countitems());

    }

    if (standardType < 0) {
        if ((*theString)==lastModelParameterList) {
            if (lastMatrixDeclared>=0) {
                _SimpleList p;
                _Variable *theM = LocateVar (modelMatrixIndices.lData[lastMatrixDeclared]);
                {
                    _AVLList pA (&p);
                    theM->ScanForVariables (pA,false);
                    pA.ReorderList();
                }
                res = p.lLength;
            }
        }
    }
    return new _Constant (res);
}
