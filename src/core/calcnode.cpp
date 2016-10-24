/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon@cfenet.ubc.ca)
 Steven Weaver (sweaver@temple.edu)
 
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

#include "math.h"
#include "ctype.h"
#include "string.h"
#include "calcnode.h"
#include "scfg.h"
#include "parser.h"
#include "function_templates.h"

#include "category.h"
#include "batchlan.h"
#include "likefunc.h"
#include "float.h"

//#ifndef     __ALTIVEC__
//#define     ALMOST_ZERO  1e-305

//#define _UBER_VERBOSE_MX_UPDATE_DUMP
#define _UBER_VERBOSE_MX_UPDATE_DUMP_LF_EVAL 1


#define     ANCESTRAL_SCALING_MAX 16
#define     ALMOST_ZERO           0.0
//#else
//#define     ALMOST_ZERO  1e-35
//#endif


#define     TREE_V_SHIFT            8.0
#define     TREE_H_SHIFT            10.0

#define     LIKELIHOOD_SCALER           1.0
#define     LIKELIHOOD_SCALER_INT       1.0

#define     DEGREES_PER_RADIAN          57.29577951308232286465


extern      _Parameter   explicitFormMatrixExponential;
extern      _String      VerbosityLevelString,
            BRANCH_LENGTH_STENCIL;

long*       nonZeroNodes = nil,
            nonZeroNodesDim = 0;
            
extern      long    likeFuncEvalCallCount;

#ifdef      __MP__
    #include <pthread.h>
    struct   ThreadMatrixTask {
        long   cID,
               tcat,
               startAt,
               endAt;
        _SimpleList* updateCN;

    };
    pthread_t*       matrixThreads = nil;
    ThreadMatrixTask*matrixTasks = nil;
    pthread_mutex_t  matrixMutex = PTHREAD_MUTEX_INITIALIZER;
#endif

char        isDefiningATree         = 0,
            takeBranchLengths       = 0,
            autoSolveBranchLengths  = 0;

_Parameter  ignoringInternalNames   = 0.0;

_Parameter  treeLayoutVert,
            scalingLogConstant    = log (1.e300);

_SimpleList convertedMatrixExpressionsL;
_AVLListX   convertedMatrixExpressions (&convertedMatrixExpressionsL);

_String     expectedNumberOfSubs  = "EXPECTED_NUMBER_OF_SUBSTITUTIONS",
            stringSuppliedLengths = "STRING_SUPPLIED_LENGTHS",
            noInternalLabels      = "NO_INTERNAL_LABELS",
            includeModelSpecs      = "INCLUDE_MODEL_SPECS",
            acceptRootedTrees    = "ACCEPT_ROOTED_TREES",
            acceptBranchLengths    = "ACCEPT_BRANCH_LENGTHS",
            autoConvertBL         = "AUTOMATICALLY_CONVERT_BRANCH_LENGTHS",
            splitNodeNames         = "SPLIT_NODE_NAMES",
            internalNodePrefix     = "INTERNAL_NODE_PREFIX",
            ignoreUserINames       = "IGNORE_INTERNAL_NODE_LABELS",
            cotNode               = "COT_NODE",
            cotSplit             = "COT_SPLIT",
            cotBranchLength         = "COT_BRANCH_LENGTH",
            cotDistance              = "COT_DISTANCE",
            cotCDF                 = "COT_CDF",
            cotSamples           = "COT_SAMPLES",
            cotSampler            = "COT_SAMPLER",
            cotToNode              = "COT_TO_NODE",
            treeOutputAVL          = "TREE_OUTPUT_OPTIONS",
            treeOutputBackground  = "TREE_OUTPUT_BACKGROUND",
            treeOutputRightMargin = "TREE_OUTPUT_RIGHT_MARGIN",
            treeOutputEmbed       = "TREE_OUTPUT_EMBED",
            treeOutputXtraMargin  = "TREE_OUTPUT_XTRA_MARGIN",
            treeOutputSplit        = "TREE_OUTPUT_BRANCH_SPLIT",
            treeOutputNotchesColor= "TREE_OUTPUT_BRANCH_NOTCHES_COLOR",
            treeOutputNotches   = "TREE_OUTPUT_BRANCH_NOTCHES",
            treeOutputLabel       = "TREE_OUTPUT_BRANCH_LABEL",
            treeOutputTLabel   = "TREE_OUTPUT_BRANCH_TLABEL",
            treeOutputColor         = "TREE_OUTPUT_BRANCH_COLOR",
            treeOutputThickness      = "TREE_OUTPUT_BRANCH_THICKNESS",
            treeOutputLinecap    = "TREE_OUTPUT_BRANCH_LINECAP",
            treeOutputDash         = "TREE_OUTPUT_BRANCH_DASH",
            treeOutputOLabel   = "TREE_OUTPUT_OVER_BRANCH",
            treeOutputSymbols   = "TREE_OUTPUT_SYMBOLS",
            treeOutputSymbolSize  = "TREE_OUTPUT_SYMBOL_SIZE",
            treeOutputExtraPS     = "TREE_OUTPUT_EXTRA_POSTSCRIPT",
            treeOutputPrefixPS      = "TREE_OUTPUT_PREFIX_POSTSCRIPT",
            treeOutputLayout   = "TREE_OUTPUT_LAYOUT",
            treeOutputNNPlaceH      = "__NODE_NAME__",
            treeOutputFSPlaceH     = "__FONT_SIZE__",
            largeMatrixBranchLengthDimension
            = "LARGE_MATRIX_BRANCH_LENGTH_MODIFIER_DIMENSION",
            largeMatrixBranchLength
            = "LARGE_MATRIX_BRANCH_LENGTH_MODIFIER",
            newNodeGraftName      = "NAME",
            newNodeGraftWhere      = "WHERE",
            newNodeGraftParent   = "PARENT",
            newNodeGraftLength   = "LENGTH",
            newNodeGraftParentLength = "PARENT_LENGTH",
            eqWithReroot        = "Equal with reroot at ",
            eqWithoutReroot       = "Equal without rerooting",
            iNodePrefix;

_Parameter  _timesCharWidths[256]= { // Hardcoded relative widths of all 255 characters in the Times font, for the use of PSTreeString
    0,0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0,0.25098,0.721569,0.721569,0.721569,0,0.721569,0.721569,
    0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0.721569,0,0.721569,0.721569,
    0.25098,0.333333,0.407843,0.501961,0.501961,0.831373,0.776471,0.180392,0.333333,0.333333,0.501961,0.564706,0.25098,0.333333,0.25098,0.278431,
    0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,0.278431,0.278431,0.564706,0.564706,0.564706,0.443137,
    0.921569,0.721569,0.666667,0.666667,0.721569,0.611765,0.556863,0.721569,0.721569,0.333333,0.388235,0.721569,0.611765,0.890196,0.721569,0.721569,
    0.556863,0.721569,0.666667,0.556863,0.611765,0.721569,0.721569,0.945098,0.721569,0.721569,0.611765,0.333333,0.278431,0.333333,0.470588,0.501961,
    0.333333,0.443137,0.501961,0.443137,0.501961,0.443137,0.333333,0.501961,0.501961,0.278431,0.278431,0.501961,0.278431,0.776471,0.501961,0.501961,
    0.501961,0.501961,0.333333,0.388235,0.278431,0.501961,0.501961,0.721569,0.501961,0.501961,0.443137,0.478431,0.2,0.478431,0.541176,0.721569,
    0.721569,0.721569,0.666667,0.611765,0.721569,0.721569,0.721569,0.443137,0.443137,0.443137,0.443137,0.443137,0.443137,0.443137,0.443137,0.443137,
    0.443137,0.443137,0.278431,0.278431,0.278431,0.278431,0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,0.501961,
    0.501961,0.4,0.501961,0.501961,0.501961,0.34902,0.454902,0.501961,0.760784,0.760784,0.980392,0.333333,0.333333,0.54902,0.890196,0.721569,
    0.713725,0.54902,0.54902,0.54902,0.501961,0.576471,0.494118,0.713725,0.823529,0.54902,0.27451,0.27451,0.309804,0.768627,0.666667,0.501961,
    0.443137,0.333333,0.564706,0.54902,0.501961,0.54902,0.611765,0.501961,0.501961,1,0.25098,0.721569,0.721569,0.721569,0.890196,0.721569,
    0.501961,1,0.443137,0.443137,0.333333,0.333333,0.54902,0.494118,0.501961,0.721569,0.168627,0.745098,0.333333,0.333333,0.556863,0.556863,
    0.501961,0.25098,0.333333,0.443137,1,0.721569,0.611765,0.721569,0.611765,0.611765,0.333333,0.333333,0.333333,0.333333,0.721569,0.721569,
    0.788235,0.721569,0.721569,0.721569,0.721569,0.278431,0.333333,0.333333,0.333333,0.333333,0.333333,0.333333,0.333333,0.333333,0.333333,0.333333
},
_maxTimesCharWidth = 0.980392;


#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif


//__________________________________________________________________________________

#define MIN_TEX_HEIGHT   50
#define MIN_TEX_WIDTH    50
#define MAX_TEX_WIDTH    160
#define MAX_TEX_HEIGHT   150
#define WIDTH_PER_BRANCH 10


//__________________________________________________________________________________

_Parameter  computeChordLength (_Parameter l, _Parameter angle, _Parameter* maxCoord = nil)
{
  _Parameter sinV        = sin(angle),
  cosV         = cos(angle);
  
  if (maxCoord) {
    maxCoord[0] = MAX (maxCoord[0], cosV*l);
    maxCoord[1] = MIN (maxCoord[1], cosV*l);
    maxCoord[2] = MAX (maxCoord[2], sinV*l);
    maxCoord[3] = MIN (maxCoord[3], sinV*l);
  }
  
  return l/MAX(fabs(sinV),fabs(cosV));
}

_String const&  getINodePrefix (void) {
  iNodePrefix = "Node";
  _PMathObj iv = FetchObjectFromVariableByType(&internalNodePrefix,STRING);
  if (iv) {
    iNodePrefix = *((_FString*)iv)->theString;
  }
  checkParameter(ignoreUserINames, ignoringInternalNames, 0.0);
  return iNodePrefix;
}


//_______________________________________________________________________________________________

_CalcNode::_CalcNode    ()
{
    theProbs = nil;    // default constructor, doesn't do much
    compExp = nil;
    matrixCache = nil;
    referenceNode=-1;
}

//_______________________________________________________________________________________________
_CalcNode::_CalcNode    (_String name, _String parms, int codeBase, _VariableContainer* theP, _AVLListXL* aCache):_VariableContainer (name, "", theP)
    // construct a node from a string of the form
    // matrix name, <optional comma separated variable declarations, inititalizations>
    // also should be passed the pointer to a container tree
{
    InitializeCN (parms, codeBase, theP, aCache);
}

//_______________________________________________________________________________________________
_CalcNode::_CalcNode    (_CalcNode* sourceNode, _VariableContainer* theP):_VariableContainer (sourceNode->ContextFreeName(), "", theP) {
    _String model = *sourceNode->GetModelName();
    InitializeCN (model, 0, theP);
    if (model.sLength) { // copy model parameter values 
        CopyMatrixParameters(sourceNode, true);
    }
}

//_______________________________________________________________________________________________
void    _CalcNode::InitializeCN     ( _String& parms, int, _VariableContainer* theP, _AVLListXL* aCache)
{
    if (theIndex < 0) return;
    
    cBase         = 0;
    theProbs      = nil;
    compExp       = nil;
    matrixCache   = nil;
    referenceNode = -1;
    slaveNodes    = 0;

    long f = parms.Find(','),
         g = -1;

    _String            matrixName (parms,0, f>=0?f-1:-1);

    InitializeVarCont (emptyString, matrixName, theP, aCache);

    if (GetModelIndex() == HY_NO_MODEL && parms.Length()) {
        f = 0;
    }

    while (f!=-1) {
        g = parms.Find (',',f+1,-1);
        if (f==0) {
            f=-1;
        }

        if (g!=-1) {
            _String  paramName (parms,f+1,g-1);
            _Formula fg (paramName, this);
        } else {
            _String  paramName (parms,f+1,g);
            _Formula fg (paramName, this);
        }
        f = g;
    }

    // attach the variables to the lists inside the node
    ScanAndAttachVariables();
    
    // check for category variables
    if (iVariables) {
        for (f = iVariables->lLength-2; f>=0 && iVariables->lData[f+1] >= 0; f-=2) {
            _Variable *theV = LocateVar(iVariables->lData[f+1]);
            if (theV->IsCategory()) { 
                        /* this has to do with local category variables; 
                           NOT TESTED and could be BROKEN
                        */
                _CategoryVariable* theCV = (_CategoryVariable*)theV;

                _Formula           newDensity,
                                   newCumulative;

                _SimpleList        iv,
                                   iv2,
                                   dv,
                                   dv2;

                for (unsigned long k = 0; k<iVariables->lLength; k+=2) {
                    iv  << iVariables->lData[k];
                    iv2 << iVariables->lData[k+1];
                }

                if (dVariables)
                    for (unsigned long k = 0; k<dVariables->lLength; k+=2) {
                        dv  << dVariables->lData[k];
                        dv2 << dVariables->lData[k+1];
                    }


                newDensity.LocalizeFormula    (theCV->GetDensity(),   *GetName(), iv, iv2, dv,dv2);
                newCumulative.LocalizeFormula (theCV->GetCumulative(),*GetName(), iv, iv2, dv,dv2);

                _CategoryVariable newCV;
                newCV.Duplicate (theCV);
                newCV.GetDensity().Duplicate((BaseRef)&newDensity);
                newCV.GetCumulative().Duplicate((BaseRef)&newCumulative);

                theV = LocateVar(iVariables->lData[f]);
                newCV.GetName()->Duplicate (theV->GetName());
                ReplaceVar(&newCV);

                categoryVariables<<iVariables->lData[f];
                categoryIndexVars<<iVariables->lData[f+1];
                iVariables->Delete(f);
                iVariables->Delete(f);
            }
        }

        if (iVariables->lLength) {
            iVariables->TrimMemory();
        } else {
            delete (iVariables);
            iVariables = nil;
        }
    }
    if (gVariables) {
        for (f = gVariables->lLength-1; f>=0; f--) {
            _Variable *theV = LocateVar(gVariables->lData[f]);
            if (theV->IsCategory()) {
                categoryVariables<<gVariables->lData[f];
                categoryIndexVars<<-1;
                gVariables->Delete(f);
            }
        }
        if (gVariables->lLength) {
            gVariables->TrimMemory();
        } else {
            delete (gVariables);
            gVariables = nil;
        }
    }

    BaseRef temp =  (variablePtrs(theIndex));
    variablePtrs[theIndex]=this->makeDynamic();
    DeleteObject(temp);

}

//__________________________________________________________________________________
void    _CalcNode::SetModel (long modelID, _AVLListXL* varCache)
{
   _VariableContainer::SetModel (modelID, varCache);
}

//_______________________________________________________________________________________________

long      _CalcNode::SetDependance (long varIndex)
{
    varIndex = _VariableContainer::SetDependance (varIndex);
    if (varIndex >= 0) {
        _SimpleList checkVars;
        _AVLList    myVars (&checkVars);
        LocateVar (varIndex)->ScanForVariables (myVars,true);

        for (long k=0; k<checkVars.lLength; k++)
            if (LocateVar(checkVars.lData[k])->IsCategory() &&(categoryVariables >> checkVars.lData[k])) {
                categoryIndexVars<<-1;
            }

    }
    return varIndex;
}

//_______________________________________________________________________________________________

void    _CalcNode::SetCodeBase (int codeBase)
{
    if (codeBase>0) {
        if ((codeBase != cBase)||!theProbs) {
#ifndef __HYALTIVEC__
            if (theProbs) {
                delete theProbs;
            }
            theProbs = new _Parameter [codeBase];
#else
            if (theProbs) {
                vec_free(theProbs);
            }
            theProbs = (_Parameter*)VecMemAllocate(codeBase*sizeof(_Parameter));
#endif
            cBase = codeBase;
            theProbs[0]=1.0;
        } else {
            theProbs[0]=1.0;
        }
    }
}

//_______________________________________________________________________________________________
void    _CalcNode::SetCompMatrix (long categID)
{
    compExp = GetCompExp (categID);
}

//_______________________________________________________________________________________________

_CalcNode::~_CalcNode (void)
{

#ifndef __HYALTIVEC__
    if (theProbs) {
        delete [] theProbs;
    }
#else
    if (theProbs) {
        vec_free(theProbs);
    }
#endif
    if (compExp && referenceNode < 0) {
        DeleteObject (compExp);
    }
}

//_______________________________________________________________________________________________

long    _CalcNode::FreeUpMemory (long)
{
    long res = 0;
    if (compExp && referenceNode < 0) {
        res = compExp->GetMySize();
        DeleteObject (compExp);
        compExp = nil;
    }
    return res;
}

//__________________________________________________________________________________

void _CalcNode::RemoveModel (void)
{
  
    if (compExp && referenceNode < 0) {
        DeleteObject (compExp);
        compExp = nil;
    }
  
    if (matrixCache) {
      
    }

    categoryVariables.Clear();
    categoryIndexVars.Clear();
    remapMyCategories.Clear();

    Clear();

}

//__________________________________________________________________________________

void _CalcNode::ReplaceModel (_String& modelName, _VariableContainer* theConext)
{
    RemoveModel    ();
    DeleteVariable (theIndex, false); // this will clean up all the node.xx variables
    InitializeCN   (modelName, 0 , theConext);
}

//_______________________________________________________________________________________________
bool    _CalcNode::MatchSubtree (_CalcNode* mNode)
{
    node <long>* myNode    = LocateMeInTree (),
                 * matchNode = mNode->LocateMeInTree ();
    if (myNode&&matchNode) {
        return myNode->compare_subtree(matchNode);
    }
    return      false;
}

//_______________________________________________________________________________________________

_Parameter  _CalcNode::ComputeBranchLength (void)
{

    if (theModel < 0) {
        return Value();
    }
  
    {
      _FString   *stencil = (_FString*)FetchObjectFromVariableByType (&BRANCH_LENGTH_STENCIL,STRING);

      if (stencil && stencil->theString->Equal (&stringSuppliedLengths)) {
          return Value();
      }
    }
    {
      _AssociativeList *lookup = (_AssociativeList*)FetchObjectFromVariableByType (&BRANCH_LENGTH_STENCIL,ASSOCIATIVE_LIST);
      if (lookup) {
        _String lookup_name = ContextFreeName();
        _Constant * value = (_Constant*)lookup->GetByKey (lookup_name, NUMBER);
        if (value) {
          return value->Value();
        }
      }
    }

  

    _Matrix     *freqMx,
                *theMx;

    bool        mbf;

    RetrieveModelComponents (theModel, theMx, freqMx, mbf);

    if (!freqMx || !theMx) {
        return Value();
    }

    _Parameter              weight = 1.0,
                            result = 0.0;

    long                    categoryCounter,
                            totalCategs = 1;

    _CategoryVariable* cVar = nil;

    if (categoryVariables.lLength) {
        for (categoryCounter = 0; categoryCounter<categoryVariables.lLength; categoryCounter++) {
            cVar = (_CategoryVariable*)LocateVar (categoryVariables.lData[categoryCounter]);
            cVar->Refresh();
            totalCategs *= cVar->GetNumberOfIntervals();
        }
    }

    freqMx = (_Matrix*)freqMx->ComputeNumeric();
    categoryCounter = 0;

    do {
        if (categoryVariables.lLength) {
            long c = categoryCounter;
            weight = 1.0;
            for (long k=categoryVariables.lLength-1; k>=0; k--) {
                cVar = (_CategoryVariable*)LocateVar (categoryVariables.lData[k]);
                long t = cVar->GetNumberOfIntervals();
                cVar->SetIntervalValue(c%t);
                weight*=cVar->GetIntervalWeight(c%t);
                c/=t;
            }
        }

        _Matrix*    theMx   = ComputeModelMatrix();
        _Parameter  expSubs = theMx->ExpNumberOfSubs (freqMx, mbf);

        _Parameter divisor;
        checkParameter (largeMatrixBranchLengthDimension, divisor, 20.);

        if (theMx->GetHDim()>divisor) {
            checkParameter (largeMatrixBranchLength, divisor, 3.);
            expSubs /= divisor;
        }

        categoryCounter++;
        result += fabs(expSubs)*weight;
    } while (categoryCounter<totalCategs);

    return result;
}


//_______________________________________________________________________________________________

_Parameter& _CalcNode::operator[] (unsigned long i)
{
    return theProbs [i];
}

//_______________________________________________________________________________________________

BaseRef _CalcNode::toStr (unsigned long)
{
    _String * res = new _String (16L, true);
    checkPointer (res);
    (*res) << theName;
    (*res) << '(';

    if (iVariables) {
        _String tempS = (long)(iVariables->lLength/2);
        (*res) << &tempS;
    } else {
        (*res) << '0';
    }
    (*res) << ',';
    if (dVariables) {
        _String tempS = (long)(dVariables->lLength/2);
        (*res) << &tempS;
    } else {
        (*res) << '0';
    }

    (*res) << ')';
    res->Finalize();
    return res;
}

//__________________________________________________________________________________

void    _CalcNode::Duplicate (BaseRef theO)
{
    _VariableContainer::Duplicate (theO);
    cBase         = 0;
    compExp       = nil;
    matrixCache   = nil;
    theProbs      = nil;
    lastState     = -1;
    referenceNode = -1;
    slaveNodes    = 0;
}

//_______________________________________________________________________________________________

bool        _CalcNode::HasChanged(bool)
{
    if (_VariableContainer::HasChanged()) {
        return true;
    }
    for (long i = 0; i<categoryVariables.lLength; i++)
        if (LocateVar (categoryVariables.lData[i])->HasChanged()) {
            return true;
        }
    return false;
}

//_______________________________________________________________________________________________

long        _CalcNode::CheckForReferenceNode(void)
{
    long rN = -1,
         idx = 0;

    // check if all independents are global first
    long modIdx = GetModelIndex();

    if (modIdx != HY_NO_MODEL) {

        if (iVariables && iVariables->lLength) {
            return -1;
        }

        if (dVariables)
            for (idx = 0; idx < dVariables->lLength; idx+=2) {
                if (dVariables->lData[idx+1]>=0) {
                    bool       good = false;
                    _Variable* thisDep = LocateVar (dVariables->lData[idx]);
                    //while (thisDep->NumberOperations() == 1)
                    while (thisDep->varFormula && thisDep->varFormula->NumberOperations () == 1) {
                        _Operation* op = (_Operation*)thisDep->varFormula->GetList() (0);
                        long        isVar = op->GetAVariable();
                        if (isVar >= 0) {
                            thisDep = LocateVar (isVar);
                            if (thisDep->IsIndependent()) {
                                good = true;
                                break;
                            }
                        } else {
                            break;
                        }
                    }

                    if (good) {
                        if (thisDep->IsGlobal()) {
                            continue;
                        } else {
                            _String varName = *thisDep->GetName();
                            long    dot = varName.FindBackwards ('.',0,-1);
                            if (dot > 0) {
                                varName.Trim (0,dot-1);
                                dot = LocateVarByName (varName);
                                if (dot < 0) {
                                    break;
                                }

                                if (rN == -1) {
                                    thisDep = FetchVar (dot);

                                    if (thisDep->ObjectClass () != TREE_NODE) {
                                        break;
                                    }

                                    if (((_CalcNode*)thisDep)->GetModelIndex() != modIdx) {
                                        break;
                                    }

                                    rN = thisDep->GetAVariable();
                                } else {
                                    if (rN != variableNames.GetXtra(dot)) {
                                        break;
                                    }
                                }
                            } else {
                                break;
                            }
                        }
                    }
                }
            }
    }
    return rN;
}

//_______________________________________________________________________________________________

bool        _CalcNode::NeedNewCategoryExponential(long catID) const
{
    if (isInOptimize && referenceNode>=0) {
        return ((_CalcNode*)LocateVar(referenceNode))->NeedNewCategoryExponential(catID);
    }

    if (_VariableContainer::NeedToExponentiate(catID>=0)) {
        return true;
    }

    if (catID==-1) {
        if (!compExp) {
            return true;
        }

        for (unsigned long i = 0; i<categoryVariables.lLength; i++)
            if (LocateVar (categoryVariables.lData[i])->HasChanged()) {
                return true;
            }
    } else {
        if (!GetCompExp(catID)) {
            return true;
        }

        for (unsigned long i = 0; i<categoryVariables.lLength; i++)
            if (((_CategoryVariable*)LocateVar (categoryVariables.lData[i]))->HaveParametersChanged(remapMyCategories.lData[catID*(categoryVariables.lLength+1)+i+1])) {
                return true;
            }
    }
    return false;
}

//_______________________________________________________________________________________________
bool        _CalcNode::RecomputeMatrix  (long categID, long totalCategs, _Matrix* storeRateMatrix, _List* queue, _SimpleList* tags, _List* bufferedOps)
{
    // assumed that NeedToExponentiate was called prior to this function

    if (isInOptimize) {
      if (referenceNode >= 0) {
        //fprintf (stderr, "\n\n********** REFERENCE NODE ******************\n\n");
        _CalcNode* rN = (_CalcNode*)LocateVar(referenceNode);
        rN->RecomputeMatrix (categID, totalCategs, storeRateMatrix);
        
        if (totalCategs>1) {
          matrixCache[categID] = rN->matrixCache[categID];
          compExp = matrixCache[categID];
        } else {
          compExp = rN->compExp;
        }
        
        return false;
      } else {
        if (referenceNode<-1) {
          slaveNodes++;
          if (slaveNodes>1) {
            if (slaveNodes == -referenceNode) {
              slaveNodes = 0;
            }
            return false;
          }
        }
      }
    }


    _Variable* curVar, *locVar;


    if (iVariables)
        for (unsigned long i=0; i<iVariables->lLength; i+=2)
            if (iVariables->lData[i+1]>=0) {
                curVar = LocateVar (iVariables->lData[i+1]);
                if (curVar->IsIndependent()) {
                    locVar = LocateVar (iVariables->lData[i]);
                    curVar->SetValue(locVar->Compute());
                    #ifdef _UBER_VERBOSE_MX_UPDATE_DUMP
                      if (1 || likeFuncEvalCallCount == _UBER_VERBOSE_MX_UPDATE_DUMP_LF_EVAL) {
                        fprintf (stderr, "[_CalcNode::RecomputeMatrix] Node %s, var %s, value = %15.12g\n", GetName()->sData, curVar->GetName()->sData, curVar->Compute()->Value());
                      }
                    #endif
                }
            }

    if (dVariables)
        for (unsigned long i=0; i<dVariables->lLength; i+=2)
            if (dVariables->lData[i+1]>=0) {
                curVar = LocateVar (dVariables->lData[i+1]);
                if (curVar->IsIndependent()) {
                    locVar = LocateVar (dVariables->lData[i]);
                    curVar->SetValue(locVar->Compute());
                    #ifdef _UBER_VERBOSE_MX_UPDATE_DUMP
                      if (1 || likeFuncEvalCallCount == _UBER_VERBOSE_MX_UPDATE_DUMP_LF_EVAL) {
                        fprintf (stderr, "[_CalcNode::RecomputeMatrix] Node %s, var %s, value = %15.12g\n", GetName()->sData, curVar->GetName()->sData, curVar->Compute()->Value());
                      }
                    #endif
                }
            }
    
    #ifdef _UBER_VERBOSE_MX_UPDATE_DUMP
      if (1|| likeFuncEvalCallCount == _UBER_VERBOSE_MX_UPDATE_DUMP_LF_EVAL && gVariables) {
        for (unsigned long i=0; i<gVariables->lLength; i++) {
          _Variable* curVar = LocateVar(gVariables->GetElement(i));
          fprintf (stderr, "[_CalcNode::RecomputeMatrix] Node %s, var %s, value = %15.12g\n", GetName()->sData, curVar->GetName()->sData, curVar->Compute()->Value());
        }
      }
    #endif


    for (unsigned long i=0; i<categoryVariables.lLength; i++) {
        if (categoryIndexVars.lData[i]<0) {
            continue;
        }
        curVar = LocateVar (categoryIndexVars.lData[i]);
        locVar = LocateVar (categoryVariables.lData[i]);
        curVar->SetValue(locVar->Compute());
    }

    if (!storeRateMatrix) {
      if (totalCategs>1) {
  #ifdef _UBER_VERBOSE_MX_UPDATE_DUMP
        fprintf (stderr, "[_CalcNode::RecomputeMatrix] Deleting category %ld for node %s at %p\n", categID, GetName()->sData, GetCompExp(categID));
  #endif
        DeleteObject(GetCompExp(categID, true));
        
      } else {
        if (compExp) {
          DeleteObject (compExp);
          compExp = nil;
        }
      }
    }

    bool    isExplicitForm  = HasExplicitFormModel ();
  
    if (isExplicitForm && bufferedOps) {
        _Matrix * bufferedExp = (_Matrix*)GetExplicitFormModel()->Compute (0,nil, bufferedOps);
        #ifdef _UBER_VERBOSE_MX_UPDATE_DUMP
            fprintf (stderr, "[_CalcNode::RecomputeMatrix] Setting (buffered) category %ld/%ld for node %s\n", categID, totalCategs, GetName()->sData);
         #endif
        SetCompExp ((_Matrix*)bufferedExp->makeDynamic(), totalCategs>1?categID:-1);
        return false;
    }

    unsigned long      previous_length = queue && tags ? queue->lLength: 0;
    _Matrix * myModelMatrix = GetModelMatrix(queue,tags); 
    
    if (isExplicitForm && !myModelMatrix) { // matrix exponentiations got cached
        if (queue && queue->lLength > previous_length) {
            return true;
        } else {
            WarnError ("Internal error");
            return false;
        }
    } else {

        if (myModelMatrix->MatrixType()!=_POLYNOMIAL_TYPE) {
            _Matrix *temp = nil;
            if (isExplicitForm) {
                temp = (_Matrix*)myModelMatrix->makeDynamic();
            } else {
                temp = (_Matrix*)myModelMatrix->MultByFreqs(theModel);
            }

            if (dVariables)
                for (unsigned long i=0; i<dVariables->lLength; i+=2)
                    if (dVariables->lData[i+1]>=0) {
                        curVar = LocateVar (dVariables->lData[i+1]);
                        if (!curVar->IsIndependent()) {
                            locVar = LocateVar (dVariables->lData[i]);
                            if (locVar->IsIndependent()) {
                                locVar->SetValue (curVar->Compute());
                            }
                        }
                    }

            if (storeRateMatrix) {
                storeRateMatrix->Duplicate(temp);
                return isExplicitForm;
            }


            if (queue) {
                (*queue) << temp;
                if (tags) {
                    (*tags) << (isExplicitForm);
                }
                return isExplicitForm;
            }

            #ifdef _UBER_VERBOSE_MX_UPDATE_DUMP
                fprintf (stderr, "[_CalcNode::RecomputeMatrix] Setting category %ld/%ld for node %s\n", categID, totalCategs, GetName()->sData);
            #endif
            SetCompExp ((_Matrix*)(isExplicitForm?temp:temp->Exponentiate()), totalCategs>1?categID:-1);

        } else {
            compExp = (_Matrix*)myModelMatrix->Evaluate(false);
        }
    }
    return false;
}

//_______________________________________________________________________________________________
void        _CalcNode::SetCompExp  (_Matrix* m, long catID)
{
    compExp = m;
    if (catID >= 0 && matrixCache) {
        if (remapMyCategories.lLength) {
            catID = remapMyCategories.lData[catID*(categoryVariables.lLength+1)];
        }
        matrixCache[catID] = compExp;
    }
}
//_______________________________________________________________________________________________
_Matrix*        _CalcNode::ComputeModelMatrix  (bool)
{
    // assumed that NeedToExponentiate was called prior to this function
    _Variable   * curVar,
                *locVar;

    if (iVariables)
        for (unsigned long i=0; i<iVariables->lLength; i+=2)
            if (iVariables->lData[i+1]>=0) {
                curVar = LocateVar (iVariables->lData[i+1]);
                if (curVar->IsIndependent()) {
                    locVar = LocateVar (iVariables->lData[i]);
                    curVar->SetValue(locVar->Compute());
                }
            }

    if (dVariables)
        for (unsigned long i=0; i<dVariables->lLength; i+=2)
            if (dVariables->lData[i+1]>=0) {
                curVar = LocateVar (dVariables->lData[i+1]);
                if (curVar->IsIndependent()) {
                    locVar = LocateVar (dVariables->lData[i]);
                    curVar->SetValue(locVar->Compute());
                }
            }

    _Matrix * modelMx = GetModelMatrix();
    if (modelMx && modelMx->ObjectClass()==MATRIX && modelMx->MatrixType()!=_POLYNOMIAL_TYPE) {
        return (_Matrix*)modelMx->ComputeNumeric();
    }

    return nil;
}

//_______________________________________________________________________________________________

_Matrix*    _CalcNode::GetCompExp       (long catID, bool doClear) const
{
    if (catID==-1) {
        return compExp;
    } else {
        if (remapMyCategories.lLength) {
            catID = remapMyCategories.lData[catID * (categoryVariables.lLength+1)];
        }
      //if (matrixCache)
      //    printf ("%d %d %d\n", remapMyCategories[0*catID * (categoryVariables.lLength+1)], remapMyCategories[catID * (categoryVariables.lLength+1)], remapMyCategories[2*catID * (categoryVariables.lLength+1)]);
        _Matrix* ret = matrixCache?matrixCache[catID]:compExp;
      if (doClear && matrixCache) {
        matrixCache[catID] = nil;
      }
      return ret;
    }
}

//_______________________________________________________________________________________________

BaseRef     _CalcNode::makeDynamic(void)
{
    _CalcNode* res = new (_CalcNode);
    checkPointer(res);
    res->_VariableContainer::Duplicate (this);
    res->categoryVariables.Duplicate ((BaseRef)&categoryVariables);
    //res->randomVariables.Duplicate ((BaseRef)&randomVariables);
    res->categoryIndexVars.Duplicate ((BaseRef)&categoryIndexVars);
    //res->randomIndexVars.Duplicate ((BaseRef)&randomIndexVars);
    res->theValue = theValue;
    res->cBase = cBase;
    if (cBase) {
        res->theProbs = new _Parameter [cBase];
        checkPointer(res->theProbs);
        memcpy (res->theProbs, theProbs, sizeof(_Parameter)*cBase);
    } else {
        res->theProbs = nil;
    }
    res->compExp = compExp;
    if (compExp) {
        compExp->nInstances++;
    }
    res->referenceNode = referenceNode;
    res->slaveNodes    = slaveNodes;
    return res;
}

//_______________________________________________________________________________________________

_TreeTopology::_TreeTopology ()
{
    rooted = UNROOTED;
}

//_______________________________________________________________________________________________

_TheTree::_TheTree ()
{
    categoryCount           = 1;
    rootIChildrenCache      = nil;
    marginalLikelihoodCache = nil;
    nodeMarkers             = nil;
    nodeStates              = nil;
    aCache                  = nil;
}       // default constructor - doesn't do much

//_______________________________________________________________________________________________

_TreeTopology::~_TreeTopology (void)
{
    if (theRoot) {
        theRoot->delete_tree();
        delete (theRoot);
        theRoot = nil;
    }
    if (compExp) {
        DeleteObject (compExp);
        compExp = nil;
    }
}

//_______________________________________________________________________________________________

_TheTree::~_TheTree (void)
{
    if (rootIChildrenCache) {
        free (rootIChildrenCache);
        rootIChildrenCache = nil;
    }
    if (marginalLikelihoodCache) {
        free (marginalLikelihoodCache);
        marginalLikelihoodCache = nil;
    }
    if (nodeMarkers) {
        free (nodeMarkers);
        nodeMarkers = nil;
    }
    if (nodeStates) {
        free (nodeStates);
        nodeMarkers = nil;
    }
    DeleteObject (aCache);

}

//_______________________________________________________________________________________________

void    _TheTree::PurgeTree (void) {
    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER);
  
    while (_CalcNode * iterator = ti.Next()) {
        DeleteVariable (*iterator->GetName());
    }
  
    theRoot->delete_tree (true);
}



//_______________________________________________________________________________________________

void    _TheTree::PreTreeConstructor (bool)
{
    rooted                  = UNROOTED;
    rootIChildrenCache      = nil;
    marginalLikelihoodCache = nil;
    nodeMarkers             = nil;
    nodeStates              = nil;
    categoryCount           = 1;

    aCache                  = new _AVLListXL (new _SimpleList);

    convertedMatrixExpressionsL.ClearFormulasInList();
    convertedMatrixExpressions.Clear();

    getINodePrefix();
  
}

//_______________________________________________________________________________________________

void    _TreeTopology::PreTreeConstructor (bool) {
    rooted                  = UNROOTED;
    compExp                 = (_Matrix*)checkPointer(new _GrowingVector);

    getINodePrefix ();
}

//_______________________________________________________________________________________________

void    _TheTree::PostTreeConstructor (bool dupMe)
{
    _Parameter acceptRTs = 0.0;
    checkParameter (acceptRootedTrees,acceptRTs, 0.0);

    DeleteObject (aCache->dataList);
    DeleteObject (aCache);
    aCache = nil;

    convertedMatrixExpressionsL.ClearFormulasInList();
    convertedMatrixExpressions.Clear();

    while (theRoot->get_num_nodes() == 1) { // dumb tree w/ an extra top level node
        // remove the root
        node<long> *node_temp = theRoot->go_down(1);
        if (!node_temp) {
            WarnError (_String("Vacuos Tree Supplied"));
            isDefiningATree = 0;
            return;
        }
        if (node_temp->get_num_nodes()) {
            _String pp = *LocateVar(theRoot->get_data())->theName;
            DeleteVariable(pp);
            delete node_temp->get_parent();
            node_temp->detach_parent();
            theRoot = node_temp;
        } else {
            break;
        }
    }

    if (theRoot->get_num_nodes() == 2) { // rooted tree - check
        if (acceptRTs<0.1) {
            long  i;
            for (i = 1; i<=2; i++) {
                node<long> *node_temp = theRoot->go_down(i);
                if (node_temp->get_num_nodes()) { // an internal node - make it a root
                    node_temp->detach_parent();
                    //node_temp = theRoot->go_down((i==1)?2:1);
                    _CalcNode * thecn = (_CalcNode*)LocateVar(theRoot->get_data());
                    _String pp = *thecn->theName;
                    DeleteVariable(pp);
                    if (i==1) {
                        node_temp->add_node(*theRoot->go_down(2));
                        delete theRoot;
                        theRoot = node_temp;
                        rooted = ROOTED_LEFT;
                    } else {
                        node_temp->prepend_node(*theRoot->go_down(1));
                        delete theRoot;
                        theRoot = node_temp;
                        rooted = ROOTED_RIGHT;
                    }
                    thecn = (_CalcNode*)LocateVar(theRoot->get_data());
                    pp=*thecn->theName;
                    DeleteVariable(pp, false);
                    if (i==1) {
                        ReportWarning (_String("Rooted tree. Removing one branch - the left root child has been promoted to be the new root"));
                    } else {
                        ReportWarning (_String("Rooted tree. Removing one branch - the right root child has been promoted to be the new root"));
                    }
                    break;
                }
            }
            if (i==3) {
                ReportWarning ((_String("One branch tree supplied - hopefully this IS what you meant to do.")));
                node<long> *node_temp = theRoot->go_down(1);
                node_temp->detach_parent();
                _CalcNode * thecn = (_CalcNode*)LocateVar(theRoot->get_data());
                _String pp = *thecn->theName;
                DeleteVariable(pp);
                node_temp->add_node(*theRoot->go_down(2));
                delete theRoot;
                theRoot = node_temp;
                rooted = ROOTED_LEFT;
                thecn = (_CalcNode*)LocateVar(theRoot->get_data());
                pp=*thecn->theName;
                DeleteVariable(pp, false);
                ReportWarning (_String("Rooted tree. Removing one branch - the left root child has been promoted to be the new root"));
                //PurgeTree();
            }
        }
    }

    if (!theRoot) {
        WarnError ("Invalid tree/topology string specification.");
    } else {
        /*_CalcNode* theN = StepWiseTraversal(TRUE);
        while (theN)
        {
            theN->SetCodeBase(4);
            theN = StepWiseTraversal();
        }*/
        BaseRef temp =  variablePtrs(theIndex);

        if (dupMe) {
            variablePtrs[theIndex]=this->makeDynamic();
        } else {
            variablePtrs[theIndex]=this;
        }

        DeleteObject(temp);
    }
}
//_______________________________________________________________________________________________

void    _TreeTopology::PostTreeConstructor (bool dupMe)
{
    BaseRef temp =  variablePtrs(theIndex);
    variablePtrs[theIndex]=dupMe ? this->makeDynamic() : this;
    DeleteObject(temp);
}

//_______________________________________________________________________________________________
_TheTree::_TheTree              (_String name, _String& parms, bool dupMe):_TreeTopology (&name)
{
    PreTreeConstructor   (dupMe);
    if (MainTreeConstructor  (parms)) {
        PostTreeConstructor  (dupMe);
    }
}

//_______________________________________________________________________________________________
_TheTree::_TheTree              (_String name, _TreeTopology* top):_TreeTopology (&name)
{
    PreTreeConstructor   (false);
    if (top->theRoot) {
        isDefiningATree         = 1;
        theRoot                 = top->theRoot->duplicate_tree ();
      
        node_iterator<long>  tree_iterator (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
      
        while (node<long>*topTraverser = tree_iterator.Next()) {
            _Parameter   nodeVal = top->compExp->theData[topTraverser->in_object];
            _String      nodeVS,
                         nodeName       (*(_String*)top->flatTree(topTraverser->in_object)),
                         nodeParams     (*(_String*)top->flatCLeaves(topTraverser->in_object));

            if (!nodeName.IsValidIdentifier(false)) {
                nodeName.ConvertToAnIdent (false);
            }

            if (nodeVal != 0.0) {
                nodeVS = nodeVal;
            }
            FinalizeNode (topTraverser, 0, nodeName, nodeParams, nodeVS);
        }
      
        isDefiningATree         = 0;
        PostTreeConstructor      (false);
    } else {
        WarnError ("Can't create an emptyString tree");
        return;
    }
}

//_______________________________________________________________________________________________
_TheTree::_TheTree              (_String name, _TheTree* otherTree):_TreeTopology (&name)
{
    PreTreeConstructor   (false);
    if (otherTree->theRoot) {
        isDefiningATree         = 1;
        theRoot                 = otherTree->theRoot->duplicate_tree ();
        
        node_iterator<long>  tree_iterator (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
        while (node<long>*topTraverser = tree_iterator.Next()) {
            _CalcNode   *sourceNode = (_CalcNode*)LocateVar(topTraverser->in_object),
                          copiedNode (sourceNode, this);
            topTraverser->init (copiedNode.theIndex);
        }
        
        isDefiningATree         = 0;
        PostTreeConstructor      (false);
    } else {
        WarnError ("Can't create an emptyString tree");
        return;
    }
}


  //_______________________________________________________________________________________________
_TreeTopology::_TreeTopology (_TheTree *top):_CalcNode (*top->GetName(), emptyString) {
  PreTreeConstructor   (false);
  if (top->theRoot) {
    isDefiningATree         = 1;
    theRoot                 = top->theRoot->duplicate_tree ();
    
    node_iterator<long>  tree_iterator (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
    while (node<long>* iterator = tree_iterator.Next()) {
      if (!iterator->is_root()) {
        _String              nodeVS,
                             *nodeSpec,
                             nodeName;
        
        top->GetBranchValue (iterator,nodeVS);
        top->GetNodeName    (iterator,nodeName);
        nodeSpec = map_node_to_calcnode(iterator)->GetBranchSpec();
        
        FinalizeNode (iterator, 0, top->GetNodeName    (iterator), *nodeSpec, nodeVS);
        DeleteObject (nodeSpec);
      }
    }
    isDefiningATree         = 0;
  } else {
    WarnError ("Can't create an emptyString tree");
    return;
  }
}

//_______________________________________________________________________________________________
_TreeTopology::_TreeTopology    (_String name, _String& parms, bool dupMe, _AssociativeList* mapping):_CalcNode (name,emptyString)
// builds a tree from a string
{
    PreTreeConstructor   (dupMe);
    if (MainTreeConstructor  (parms, false, mapping)) {
        PostTreeConstructor  (dupMe);
    } else {
        DeleteObject     (compExp);
        compExp = nil;
    }
}

//_______________________________________________________________________________________________
_TreeTopology::_TreeTopology    (_String* name):_CalcNode (*name,emptyString) {

}

//_______________________________________________________________________________________________

bool _MainTreeConstructor_error (const _String& error, const _String& tree_string, unsigned long index) {
  isDefiningATree = 0;
  WarnError (   error & ", in the following string context " &
                tree_string.Cut(index>31L?index-32L:0L,index)&
                "<ERROR HERE>"&
                tree_string.Cut(index+1L,tree_string.sLength-index>32L?index+32L:-1L)
             );
  
  return false;
}

//_______________________________________________________________________________________________
bool    _TreeTopology::MainTreeConstructor  (_String& parms, bool checkNames, _AssociativeList* mapping) {
    long i,
         nodeCount=0,
         lastNode;
    
    _Parameter      checkABL;

    checkParameter  (acceptBranchLengths, checkABL, 1.0);
    takeBranchLengths = !CheckEqual (checkABL, 0.0);

    checkParameter  (autoConvertBL, checkABL, 0.0);
    autoSolveBranchLengths = CheckEqual (checkABL, 1.0);


    _SimpleList nodeStack,
                nodeNumbers;

    _String     nodeName,
                nodeParameters,
                nodeValue,
                nodeComment;

    char        lastChar    = 0;
    bool        isInLiteral = false;

    node<long>* currentNode = theRoot = nil,
              * newNode     = nil,
              * parentNode  = nil;

    isDefiningATree         = 1;

    for (i=0; i<parms.sLength; i++) {
        switch (parms[i]) {
        case '(': { // creating a new internal node one level down
            // a new node
            newNode = new node<long>;

           if (lastChar == '(' || lastChar == ',') {
                 if (!currentNode) {
                   return _MainTreeConstructor_error (_String ("Unexpected '") & lastChar & "'", parms, i);
                 }
                currentNode->add_node (*newNode);
            } else {
                if (theRoot) {
                    parentNode = currentNode->get_parent();
                    if (!parentNode) {
                      return _MainTreeConstructor_error (_String ("Unexpected '") & lastChar & "'", parms, i);
                    } else {
                        parentNode->add_node (*newNode);
                    }
                } else {
                    theRoot = newNode;
                }
                currentNode = newNode;
                nodeStack<<(long)currentNode;
                nodeNumbers<<nodeCount;
                newNode = new node<long>;
                currentNode->add_node (*newNode);
                nodeCount++;
            }
            currentNode = newNode;
            nodeStack<<(long)currentNode;
            nodeNumbers<<nodeCount;
            nodeCount++;
            break;
        }

        case ',':
        case ')': { // creating a new node on the same level and finishes updating the list of parameters
            lastNode = nodeStack.lLength-1;
            if (lastNode<0) {
                return _MainTreeConstructor_error (_String ("Unexpected '") & parms[i] & "'", parms, i);
            }
            parentNode = (node<long>*)nodeStack(lastNode);
            if (mapping) {
              _FString * mapped_name = (_FString*)mapping->GetByKey (nodeName, STRING);
              if (mapped_name) {
                nodeName = *mapped_name->theString;
              }
            }
            FinalizeNode (parentNode, nodeNumbers(lastNode), nodeName, nodeParameters, nodeValue, &nodeComment);
            nodeName = emptyString;
            nodeStack.Delete(lastNode, false);
            nodeNumbers.Delete(lastNode, false);

            if (parms[i]==',') { // also create a new node on the same level
                checkPointer (newNode = new node<long>);
                if (!(parentNode = parentNode->get_parent())) {
                    return _MainTreeConstructor_error (_String ("Unexpected '") & parms[i] & "'", parms, i);
                }
                currentNode = newNode;
                parentNode->add_node(*currentNode);
                nodeStack<<(long)currentNode;
                nodeNumbers<<nodeCount;
                nodeCount++;
            }
            break;
        }

        case '{' : { // parameter list definition
            lastNode = parms.Find ("}",i+1,-1);
            if (lastNode<0) {
                return _MainTreeConstructor_error (_String ("'{' has no matching '}'"), parms, i);
            }
            nodeParameters = parms.Cut (i+1,lastNode-1);
            i = lastNode;
            break;
        }

         case '[' : { // hackish Newick annotation
            lastNode = parms.Find ("]",i+1,-1);
            if (lastNode<0) {
                return _MainTreeConstructor_error (_String ("'[' has no matching ']'"), parms, i);
            }
            nodeComment = parms.Cut (i+1,lastNode-1);
            i = lastNode;
            break;
        }
        
        case ':' : { // tree branch definition
            lastNode = i+1;
            char c = parms[lastNode];

            while ( isspace (c) )
                if (lastNode<parms.sLength) {
                    c = parms[++lastNode];
                } else {
                    break;
                }

            if ( lastNode<parms.sLength )
                while ( (c<='9' && c>='0') || c=='.' ||c=='-' ||c=='+' || c=='e' || c=='E') {
                    if (lastNode<parms.sLength) {
                        c = parms[++lastNode];
                    } else {
                        break;
                    }
                }

            nodeValue = parms.Cut(i,lastNode-1);
            i = lastNode-1;
            break;
        }

        default: { // node name
            lastNode = i;

            char c = parms[lastNode];

            if (c==';') {
                break;
            }

            if (isspace(c)) {
                continue;
            }

            if (c == '\'') {
                c = parms[lastNode++];
                isInLiteral = true;
                i++;
            }

            if (!isInLiteral && !(isalnum(c)|| c=='_')) {
                return _MainTreeConstructor_error (_String("Node names should begin with a letter, a number, or an underscore"), parms, i);
            }

            if (checkNames)
                while (isalnum(c)||c=='_')
                    if (lastNode<parms.sLength) {
                        lastNode++;
                        c = parms[lastNode];
                    } else {
                        break;
                    }
            else
                while (isInLiteral || !(c==',' || c==':' || c==')' || c=='(' || c=='{' ||c== '}' || isspace(c)))
                    if (lastNode<parms.sLength) {
                        lastNode++;
                        c = parms[lastNode];
                        if (c == '\'') {
                            if (isInLiteral) {
                                break;
                            } else {
                              return _MainTreeConstructor_error (_String ("Unexpected '") & c & "'", parms, i);
                            }
                        }
                    } else {
                        break;
                    }

            if (isInLiteral) {
                if (c!='\'') {
                  return _MainTreeConstructor_error (_String ("Unterminated literal"), parms, i);
                }
                nodeName        = parms.Cut(i,lastNode-1);
                i               = lastNode;
                isInLiteral     = false;
            } else {
                nodeName        = parms.Cut(i,lastNode-1);
                i               = lastNode-1;
            }
            break;
        }
        }

        if ((lastChar = parms[i])==';') {
            break;
        }
    }

    if (!theRoot) {
      isDefiningATree = 0;
      WarnError ("Can't create emptyString trees.");
      return false;
    }

    lastNode = nodeStack.lLength-1;

    while (lastNode>=0) {
        parentNode = (node<long>*)nodeStack(lastNode);
        if (mapping) {
          _FString * mapped_name = (_FString*)mapping->GetByKey (nodeName, STRING);
          if (mapped_name) {
            nodeName = *mapped_name->theString;
          }
        }
        FinalizeNode (parentNode, nodeNumbers(lastNode), nodeName, nodeParameters, nodeValue, &nodeComment);
        lastNode--;
    }


    isDefiningATree = 0;
    return true;

}
//_______________________________________________________________________________________________


bool    _TheTree::FinalizeNode (node<long>* nodie, long number , _String nodeName, _String& nodeParameters, _String& nodeValue, _String* nodeComment)
{
    bool isAutoGenerated = (nodeName.sLength == 0 || (!CheckEqual(ignoringInternalNames,0.0) && nodie->get_num_nodes()>0));
    if (isAutoGenerated) {
        nodeName = iNodePrefix & number;
    } else {
        
        if (!nodeName.IsValidIdentifier(false)) {
            _String oldName (nodeName);
            nodeName.ConvertToAnIdent();
            ReportWarning (_String ("Automatically renamed ") & oldName & " to " & nodeName & " in order to create a valid HyPhy identifier");
        }
    }
    if (nodie == theRoot) {
        nodeParameters = emptyString;
        nodeValue      = emptyString;
    } else {
        if (!nodeParameters.sLength && lastMatrixDeclared!=-1) {
            nodeParameters=*(((_String**)modelNames.lData)[lastMatrixDeclared]);
        }

        if (nodeParameters.sLength) {
            ReportWarning ((_String("Model ")&nodeParameters&_String(" assigned to ")& nodeName));
        } else {
            ReportWarning (_String("No nodel was assigned to ")& nodeName);
        }

    }

    char saveIDT = isDefiningATree;
    isDefiningATree = 2;
    _CalcNode cNt (nodeName,nodeParameters, 4, this, aCache);
    isDefiningATree = saveIDT;
    nodie->init (cNt.theIndex);


    _Constant val (nodeValue.ProcessTreeBranchLength());
    
    if (nodeValue.Length() && takeBranchLengths) {
        if (cNt.iVariables && cNt.iVariables->lLength == 2) { // can assign default values
            bool setDef = true;
            if (autoSolveBranchLengths) {
                long nodeModelID = cNt.GetModelIndex();
                if (nodeModelID != HY_NO_MODEL) {
                    _Formula * expressionToSolveFor = nil;
                    long alreadyConverted = convertedMatrixExpressions.Find ((BaseRef)nodeModelID);
                    if (alreadyConverted < 0) {
                        _Variable   * tV, * tV2;
                        bool         mByF;
                        RetrieveModelComponents (nodeModelID, tV, tV2, mByF);
                        _String * result = ((_Matrix*)tV->GetValue())->BranchLengthExpression((_Matrix*)tV2->GetValue(),mByF);
                        if (result->sLength) {
                            expressionToSolveFor = new _Formula (*result);
                            for (unsigned long cc = 0; cc < cNt.categoryVariables.lLength; cc++) {
                                _CategoryVariable * thisCC = (_CategoryVariable *)LocateVar(cNt.categoryVariables.lData[cc]);
                                thisCC -> SetValue (new _Constant(thisCC->Mean()), false);

                            }
                        }
                        DeleteObject (result);
                    } else {
                        expressionToSolveFor = (_Formula*)convertedMatrixExpressions.GetXtra(alreadyConverted);
                    }

                    if (expressionToSolveFor != nil) {
                        _Variable * solveForMe = LocateVar (cNt.iVariables->lData[1]);
                        _Parameter modelP = expressionToSolveFor->Brent (solveForMe,solveForMe->GetLowerBound(), solveForMe->GetUpperBound(), 1e-6, nil, val.Value());
                        ReportWarning (_String("Branch parameter of ") & nodeName&" set to " & modelP);
                        LocateVar (cNt.iVariables->lData[0])->SetValue(new _Constant (modelP), false);
                        setDef = false;
                    }
                }
            }

            if (setDef) {
                LocateVar (cNt.iVariables->lData[0])->SetValue (&val);
                ReportWarning (_String("Branch parameter of ") & nodeName&" set to " & nodeValue);
            }
        } else {
            ReportWarning (nodeName&" has "& _String((long)(cNt.iVariables?cNt.iVariables->lLength/2:0)) & " parameters - branch length not assigned");
        }
    }

    _CalcNode *nodeVar = (_CalcNode*)LocateVar(cNt.theIndex);
    
    if (nodeVar == NULL) return false;

    nodeVar->SetValue (&val);

    nodeName       = emptyString;
    nodeParameters = emptyString;
    nodeValue      = emptyString;
    if (nodeComment && nodeComment->sLength)
    {
        _String commentName = *nodeVar->GetName() & "._comment";
        CheckReceptacleAndStore(&commentName, emptyString, false, new _FString (*nodeComment));
        *nodeComment    = emptyString;
    }

    nodeVar->categoryVariables.TrimMemory();
    nodeVar->categoryIndexVars.TrimMemory();
    nodeVar->_VariableContainer::TrimMemory();


    return true;
}

//_______________________________________________________________________________________________

bool    _TreeTopology::FinalizeNode (node<long>* nodie, long number , _String nodeName, _String& nodeParameters, _String& nodeValue, _String* nodeComment)
{
    if (!nodeName.sLength || (!CheckEqual(ignoringInternalNames,0.0) && nodie->get_num_nodes()>0)) {
        nodeName = iNodePrefix & number;
    }

    if (nodie==theRoot) {
        nodeParameters = "";
        nodeValue = "";
    }

    nodie->in_object = flatTree.lLength;
    flatTree          && & nodeName;
    flatCLeaves       && & nodeParameters;

    ((_GrowingVector*)compExp)->Store (nodeValue.ProcessTreeBranchLength());

    nodeName       = emptyString;
    nodeParameters = emptyString;
    nodeValue      = emptyString;
    if (nodeComment)
        *nodeComment    = emptyString;

    return true;
}

//_______________________________________________________________________________________________

node<long>* _TreeTopology::FindNodeByName (_String const* match) const {
  
    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
    while (node<long>* iterator = ni.Next()) {
        if (*match == GetNodeName  (iterator)) {
            return iterator;
        }
    }

    return nil;

}

//_______________________________________________________________________________________________

void    _TreeTopology::RemoveANode (_PMathObj nodeName) {
    if (nodeName->ObjectClass () == STRING) {
        _FString         * removeMe     = (_FString*)nodeName;
        
        node<long>* removeThisNode = FindNodeByName (removeMe->theString),
                  * parentOfRemoved;
        
        if (!removeThisNode || ( parentOfRemoved = removeThisNode->get_parent()) == nil) {
            WarnError ("Node not found in the tree or is the root node call to _TreeTopology::RemoveANode");
            return;
        }
        
        _SimpleList cleanIndices;
        
        while (parentOfRemoved) {
            cleanIndices << removeThisNode->in_object;
            parentOfRemoved->detach_child(removeThisNode->get_child_num());
            
          
            for (int orphans = 1; orphans <= removeThisNode->get_num_nodes(); orphans++) {
                parentOfRemoved->add_node(*removeThisNode->go_down(orphans));
            }
            delete removeThisNode;
            removeThisNode = parentOfRemoved;
            parentOfRemoved = parentOfRemoved->get_parent();
            
            if (parentOfRemoved == nil && removeThisNode->get_num_nodes() == 1) {
                removeThisNode = removeThisNode->go_down (1);
                parentOfRemoved = removeThisNode->get_parent();
               
            }
            
        }
        
        cleanIndices.Sort();
        flatTree.DeleteList(cleanIndices);
        flatCLeaves.DeleteList(cleanIndices);
        
        
        
        cleanIndices << flatTree.lLength + 16L; // this is so that we don't get a index error thrown at GetElement
        
        //printf ("%s\n", ((_String*)cleanIndices.toStr())->getStr());
        
        _GrowingVector* blengths = ((_GrowingVector*)compExp);
        long offset = 0L;
        _SimpleList     oldToNew;
        
        for (long k = 0L; k < blengths->used; k++) {
            if (k == cleanIndices.GetElement(offset)){
                oldToNew << -1;
                offset++;
            } else {
                blengths->theData[k-offset] = blengths->theData[k];
                oldToNew << k - offset;
            }
        }
        
        blengths->used -= offset - 1;
        node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
      
        //printf ("%s\n", ((_String*)oldToNew.toStr())->getStr());

        while (node<long>* iterator = ni.Next()) {
            //printf ("[%ld->%ld]\n", currentNode->in_object, oldToNew.GetElement(currentNode->in_object));
            iterator->in_object = oldToNew.GetElement(iterator->in_object);
        }
       
        
        
    } else {
        WarnError ("An invalid argument (not a string) supplied to _TreeTopology::RemoveANode");
    }

}


//_______________________________________________________________________________________________

void    _TreeTopology::AddANode (_PMathObj newNode)
{
    if (newNode->ObjectClass () == ASSOCIATIVE_LIST) {
        _AssociativeList * newNodeSpec = (_AssociativeList*)newNode;
        _FString         * newName     = (_FString*)newNodeSpec->GetByKey (newNodeGraftName, STRING),
                         * newLocation = (_FString*)newNodeSpec->GetByKey (newNodeGraftWhere, STRING),
                         * newParent   = (_FString*)newNodeSpec->GetByKey (newNodeGraftParent, STRING);
                         
        _Constant        * branchLengthSelf = (_Constant*)newNodeSpec->GetByKey (newNodeGraftLength, NUMBER),
                         * branchLengthParent = (_Constant*)newNodeSpec->GetByKey (newNodeGraftParentLength, NUMBER);


        if (!newLocation) {
            WarnError (_String("Missing/invalid mandatory argument (\"")&newNodeGraftWhere&"\") in call to _TreeTopology::AddANode");
            return;
        }
        
        
        if (! (newName || newParent)) {
            WarnError (_String("At least one of '") & newNodeGraftName&"', '"& newNodeGraftParent &"') must be specified in call to _TreeTopology::AddANode");
            return;
        }

        node<long>* graftAt = FindNodeByName (newLocation->theString);
        if (!graftAt || graftAt->get_parent() == nil) {
            WarnError ("Attachment node must be an exiting non-root node in call to _TreeTopology::AddANode");
            return;
        }

        node<long>* newp = newParent ? new node<long> : nil,
                  * curp = graftAt->get_parent();

        if (newp) {
          newp->set_parent  (*curp);
          newp->add_node    (*graftAt);
          curp->replace_node(graftAt,newp);
        } 

        if (newName && !newName->IsEmpty()) {
            node<long>* newt = (node<long>*) checkPointer(new node<long>);
            if (newp) {
              newp->add_node(*newt);
            } else {
              graftAt->add_node (*newt);
            }
            if (branchLengthSelf) {
              _String bl (branchLengthSelf->Value());
              FinalizeNode (newt, 0, *newName->theString,   emptyString, bl);
            } else {
              FinalizeNode (newt, 0, *newName->theString,   emptyString, emptyString);
            }
        }

        if (newp && ! newParent->IsEmpty()) {
            if (branchLengthParent) {
              _String bl (branchLengthParent->Value());
               FinalizeNode (newp, 0, *newParent->theString, emptyString, bl);
            } else {
              FinalizeNode (newp, 0, *newParent->theString, emptyString, emptyString);
            }
          
        }

    } else {
        WarnError ("An invalid argument (not an associative array) supplied to _TreeTopology::AddANode");
    }

}
//_______________________________________________________________________________________________

BaseRef  _TreeTopology::makeDynamic (void)
{
    _TreeTopology* res = new _TreeTopology;
    checkPointer(res);
    res->_CalcNode::Duplicate (this);

    res->flatTree.Duplicate (&flatTree);
    res->flatCLeaves.Duplicate (&flatCLeaves);
    if (compExp) {
        res->compExp = (_Matrix*)compExp->makeDynamic();
    } else {
        res->compExp = nil;
    }

    res->rooted = rooted;
    res->theRoot = CopyTreeStructure (theRoot,true);
    return res;
}

//_______________________________________________________________________________________________

BaseRef  _TheTree::makeDynamic (void)
{
    _TheTree* res = new _TheTree;
    checkPointer(res);
    res->_CalcNode::Duplicate (this);

    res->rooted = rooted;
    res->categoryCount = 1;
    res->theRoot = CopyTreeStructure (theRoot,true);
    return res;
}


//_______________________________________________________________________________________________

BaseRef  _TheTree::makeDynamicCopy (_String* replacementName)
{
    _TheTree* res = new _TheTree;
    checkPointer(res);

    res->rooted = rooted;
    if (theRoot) {
        _String rn = *replacementName&'.';
        res->theRoot = DuplicateTreeStructure (theRoot, &rn, true);
    } else {
        res->theRoot = nil;
    }

    res->SetIndex(variableNames.GetXtra(LocateVarByName (*replacementName)));
    res->theName = replacementName;
    res->theName->nInstances++;
    return res;
}



//_______________________________________________________________________________________________

_PMathObj       _CalcNode::Compute  (void)
{
    return this;
}


//_______________________________________________________________________________________________

node<long>*  _TheTree::DuplicateTreeStructure (node <long>* theNode, _String* replacementName,  bool)
{
    long i,j;
    node<long>* locNode = new node<long>;
    for (i=0; i<theNode->get_num_nodes(); i++) {
        locNode->add_node(*DuplicateTreeStructure (theNode->go_down(i+1), replacementName , false));
    }
    if (1)
        // process this node now
    {
        _String    replacedName = *GetName()&'.',*temp;
        _CalcNode* sourceNode = (_CalcNode*)LocateVar(theNode->get_data());
        sourceNode = (_CalcNode*)sourceNode->makeDynamic();
        _String newNodeName = (LocateVar(sourceNode->GetAVariable())->GetName())->Replace(replacedName,*replacementName,true);
        _Variable dummyVar (newNodeName);
        j = dummyVar.GetAVariable();
        temp = sourceNode->GetName();
        DeleteObject (temp);
        sourceNode->theName = dummyVar.GetName();
        sourceNode->theName->nInstances++;
        ReplaceVar (sourceNode);
        DeleteObject(sourceNode);
        sourceNode = (_CalcNode*)LocateVar (j);
        locNode->init (j);

#ifndef USE_POINTER_VC
        for (i=0; i<sourceNode->independentVars.lLength; i++) {
            newNodeName = (LocateVar(sourceNode->independentVars.lData[i])->GetName())->Replace(replacedName,*replacementName,true);
            _Variable dummyVar (newNodeName);
#ifndef USE_AVL_NAMES
            sourceNode->independentVars.lData[i] = variableReindex.lData[LocateVarByName (newNodeName)];
#else
            sourceNode->independentVars.lData[i] = variableNames.GetXtra(LocateVarByName (newNodeName));
#endif
        }

        // done with independents - now set the dependancies
        for (i=0; i<sourceNode->dependentVars.lLength; i++) {
            newNodeName = (LocateVar(sourceNode->dependentVars.lData[i])->GetName())->Replace(replacedName,*replacementName,true);
            _Variable dummyVar (newNodeName);
#ifndef USE_AVL_NAMES
            sourceNode->dependentVars.lData[i] = variableReindex.lData[LocateVarByName (newNodeName)];
#else
            sourceNode->dependentVars.lData[i] = variableNames.GetXtra(LocateVarByName (newNodeName));
#endif
            _String* newFormula = LocateVar(sourceNode->dependentVars.lData[i])->GetFormulaString();
            *newFormula = newFormula->Replace(replacedName,*replacementName,true);
            _Formula dummyF (*newFormula);
            LocateVar(sourceNode->dependentVars.lData[i])->SetFormula(dummyF);
            DeleteObject (newFormula);
        }
#else
        if (sourceNode->iVariables)
            for (i=0; i<sourceNode->iVariables->lLength; i+=2) {
                newNodeName = (LocateVar(sourceNode->iVariables->lData[i])->GetName())->Replace(replacedName,*replacementName,true);
                _Variable dummyVar (newNodeName);
#ifndef USE_AVL_NAMES
                sourceNode->iVariables->lData[i] = variableReindex.lData[LocateVarByName (newNodeName)];
#else
                sourceNode->iVariables->lData[i] = variableNames.GetXtra(LocateVarByName (newNodeName));
#endif
            }

        // done with independents - now set the dependancies

        if (sourceNode->dVariables)
            for (i=0; i<sourceNode->dVariables->lLength; i+=2) {
                newNodeName = (LocateVar(sourceNode->dVariables->lData[i])->GetName())->Replace(replacedName,*replacementName,true);
                _Variable dummyVar (newNodeName);
#ifndef USE_AVL_NAMES
                sourceNode->dVariables->lData[i] = variableReindex.lData[LocateVarByName (newNodeName)];
#else
                sourceNode->dVariables->lData[i] = variableNames.GetXtra(LocateVarByName (newNodeName));
#endif
                _String* newFormula = LocateVar(sourceNode->dVariables->lData[i])->GetFormulaString();
                *newFormula = newFormula->Replace(replacedName,*replacementName,true);
                _Formula dummyF (*newFormula);
                LocateVar(sourceNode->dVariables->lData[i])->SetFormula(dummyF);
                DeleteObject (newFormula);
            }

#endif
    }
    return locNode;
}

//_______________________________________________________________________________________________

node<long>*  _TreeTopology::CopyTreeStructure (node <long>* theNode,  bool)
{
    long i;
    node<long>* locNode = new node<long>;
    for (i=0; i<theNode->get_num_nodes(); i++) {
        locNode->add_node(*CopyTreeStructure (theNode->go_down(i+1), false));
    }

    locNode->init (theNode->in_object);
    return locNode;
}




//_______________________________________________________________________________________________

const _String    _TreeTopology::GetNodeName (node<long>* n, bool fullName) const {
    if (fullName) {
        return *GetName()&'.'&*(_String*)(flatTree.GetItem(n->in_object));
    }
    return *(_String*)(flatTree.GetItem (n->in_object));
 }

//_______________________________________________________________________________________________

_String const*    _TreeTopology::GetNodeModel (node<long>* n) const {
  return (_String*)flatCLeaves.GetItem(n->in_object);
}

//_______________________________________________________________________________________________

_String const*    _TheTree::GetNodeModel (node<long>* n) const {
    return ((_VariableContainer*)LocateVar (n->in_object))->GetModelName ();
}

//_______________________________________________________________________________________________

void    _TreeTopology::SetLeafName(long res, _String* newName)
{
    node_iterator<long> ni (theRoot);
    unsigned long leaf_count      = 0UL;
  
    while (node<long>* iterator = ni.Next()) {
      if (iterator->is_leaf()) {
        if (res == leaf_count) {
          flatTree.Replace (iterator->in_object, newName, true);
          return;
        }
        leaf_count ++;
      }
    }
}
//_______________________________________________________________________________________________

const _String    _TheTree::GetNodeName      (node<long>* n, bool fullName) const
{
    if (fullName) {
        return *map_node_to_calcnode(n)->GetName();
    }
    return map_node_to_calcnode(n)->GetName()->Cut (GetName()->sLength+1,-1);
}


//_______________________________________________________________________________________________

BaseRef     _TheTree::toStr (unsigned long) {
  return _TreeTopology::toStr();
}

  //_______________________________________________________________________________________________

_List*     _TreeTopology::MapNodesToModels (void) {
  _List* map = new _List;
  
  node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
  
  while (node<long>* iterator = ni.Next()) {
    if (iterator->is_root())
      break;
    _List * node_record = new _List;
    (*node_record) < new _String (GetNodeName(iterator)) < new _String(*GetNodeModel(iterator));
    (*map) < node_record;
  }
  return map;
}

  //_______________________________________________________________________________________________

_String const  _TreeTopology::GetNodeStringForTree                (node<long> * n , int flags) const {

  _String node_desc;
  
  if (flags & kGetNodeStringForTreeName) {
    node_desc = GetNodeName (n);
  }
  
  if (flags & kGetNodeStringForTreeModel) {
    _String const *mSpec = GetNodeModel (n);
    if (mSpec->sLength) {
      node_desc = node_desc & '{' & *mSpec & '}';
    }
  }
  return node_desc;
}

  //_______________________________________________________________________________________________

BaseRef     _TreeTopology::toStr (unsigned long) {
  _String     * res = new _String((unsigned long)128,true),
  num;
  
  _Parameter    skipILabels,
  includeMSP;
  
  checkParameter (noInternalLabels, skipILabels, 0.0);
  checkParameter (includeModelSpecs, includeMSP , 0.0);
  
  int            leaf_flag      = kGetNodeStringForTreeName,
                 inode_flag     = skipILabels < 0.1 ? kGetNodeStringForTreeName : 0;
  
  if (includeMSP > 0.5) {
    leaf_flag   = leaf_flag | kGetNodeStringForTreeModel;
    inode_flag  = inode_flag | kGetNodeStringForTreeModel;
  }
  
  if (IsDegenerate ()) {
    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
    (*res)<<'(';
    
    while (node<long>* iterator = ni.Next()) {
      (*res) << GetNodeStringForTree (iterator,  leaf_flag );
      (*res)<< (iterator->is_root() ? ')' : ',');
    }
  } else {
    
    long        level       =   0L,
    myLevel     =   0L,
    lastLevel   =   0L;
    
    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
    
    node<long>*     curNode   =   ni.Next(),
              *     nextNode;
    
    level                     = ni.Level();
    
    nextNode = ni.Next();
    
    while (nextNode) {
      if (level>lastLevel) {
        if (lastLevel) {
          (*res)<<',';
        }
        res->AppendNCopies('(', level-lastLevel);
      } else if (level<lastLevel) {
        res->AppendNCopies(')', lastLevel-level);
      } else {
        (*res)<<',';
      }
      
      (*res) << GetNodeStringForTree (curNode,  curNode->is_leaf() ? leaf_flag : inode_flag);
      
      
      lastLevel = level;
      curNode   = nextNode;
      level     = ni.Level();
      nextNode  = ni.Next();
    }
    
    res->AppendNCopies(')', lastLevel-level);
  }
  (*res)<<';';
  res->Finalize();
  return res;
}

//__________________________________________________________________________________
void _TreeTopology::toFileStr(FILE* f, unsigned long padding) {
    _String * s = (_String*)toStr(padding);
    fprintf (f, "%s", s->sData);
    DeleteObject(s);
}

//__________________________________________________________________________________

void _TheTree::CompileListOfModels (_SimpleList& l) {
  
    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER);
    while (_CalcNode* iterator = ti.Next()) {
        long    modelID = iterator->GetModelIndex();
        if (modelID != HY_NO_MODEL && l.Find(modelID) == -1 ) {
            l << modelID;
        }
    }
}

//_______________________________________________________________________________________________
void    _TheTree::SetCompMatrices (long catID) const {
    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER | _HY_TREE_TRAVERSAL_SKIP_ROOT);
    while   (_CalcNode* iterator = ti.Next()) {
      iterator->SetCompMatrix (catID);
    }
}

//__________________________________________________________________________________

void _TheTree::SetUp (void) {

    if (marginalLikelihoodCache) {
        free (marginalLikelihoodCache);
        marginalLikelihoodCache = nil;
    }

    if (nodeMarkers) {
        free (nodeMarkers);
        nodeMarkers = nil;
    }

    if (nodeStates) {
        free (nodeStates);
        nodeMarkers = nil;
    }

    flatTree.Clear();
    flatNodes.Clear();
    flatLeaves.Clear();
    flatCLeaves.Clear();

    flatParents.Clear();
    _SimpleList flatINodeParents;

    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER);

    while   (_CalcNode* iterator = ti.Next()) {
        if (ti.IsAtLeaf()) {
          flatCLeaves << iterator;
          flatLeaves << (long)ti.GetNode();
          flatParents << (long)ti.GetNode()->get_parent();
        } else {
            flatTree<<iterator;
            flatNodes<< (long)ti.GetNode();
            iterator->lastState = -1;
            flatINodeParents << (long)ti.GetNode()->get_parent();
        }
    }

    flatParents << flatINodeParents;
  
    _SimpleList parentlist (flatNodes),
                indexer (flatNodes.lLength,0,1);
  
    SortLists   (&parentlist,&indexer);
    for (long k=0L; k<flatParents.lLength; k++)
        if (flatParents.lData[k]) {
            flatParents.lData[k] = indexer.lData[parentlist.BinaryFind(flatParents.lData[k])];
        } else {
            flatParents.lData[k] = -1;
        }


    if (cBase>0) {
        marginalLikelihoodCache = (_Parameter*)MemAllocate ((flatNodes.lLength+flatLeaves.lLength)*sizeof (_Parameter)*cBase*systemCPUCount);
    }
    nodeStates                  = (long*)MemAllocate ((flatNodes.lLength+flatLeaves.lLength)*sizeof (long)*systemCPUCount);
    nodeMarkers                 = (char*)MemAllocate (flatNodes.lLength*sizeof (char)*systemCPUCount);

    unsigned long    iNodeCounter = 0UL,
                     leafCounter  = 0UL;
  
    ti.Reset();

    while   (_CalcNode* iterator = ti.Next()) {
        if (ti.IsAtLeaf()) {
            iterator->nodeIndex = leafCounter++;
        } else {
            nodeMarkers[iNodeCounter] = -1;
            for (long k=1; k<systemCPUCount; k++) {
                nodeMarkers[iNodeCounter+k*flatNodes.lLength] = -1;
            }
            iterator->nodeIndex = flatLeaves.lLength+iNodeCounter++;
            nodeStates[iterator->nodeIndex]=-1;
            for (long m=1; m<systemCPUCount; m++) {
                nodeStates[iterator->nodeIndex+m*(flatNodes.lLength+flatLeaves.lLength)] = -1;
            }
        }
    }

    BuildINodeDependancies();
}

//__________________________________________________________________________________

bool _TheTree::AllBranchesHaveModels (long matchSize) const {
  
  _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER | _HY_TREE_TRAVERSAL_SKIP_ROOT);

  while (_CalcNode* iterator = ti.Next()) {
    _Matrix * model = iterator->GetModelMatrix();
    if (model && (matchSize < 0L || model->GetHDim()==matchSize)) {
      continue;
    }
    return false;
  }
 
  return true;
}

//__________________________________________________________________________________

_String*    _TheTree::TreeUserParams (void) const {
    _String * result = new _String (16L, true);

    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER);
    while (_CalcNode* iterator = ti.Next()) {
        result->AppendNewInstance(iterator->GetSaveableListOfUserParameters());
    }

    result->Finalize();
    return result;
}

//__________________________________________________________________________________

_PMathObj _TreeTopology::ExecuteSingleOp (long opCode, _List* arguments, _hyExecutionContext* context) {
  
  
  switch (opCode) { // first check operations without arguments
    case HY_OP_CODE_ABS: // Abs
      return FlatRepresentation();
    case HY_OP_CODE_BRANCHCOUNT: //BranchCount
      return BranchCount();
    case HY_OP_CODE_TIPCOUNT: // TipCount
      return TipCount();
    case HY_OP_CODE_TYPE: // Type
      return Type();
  }
  
  _MathObject * arg0 = _extract_argument (arguments, 0UL, false);
  
  switch (opCode) { // next check operations without arguments or with one argument
    case HY_OP_CODE_ADD: // +
      if (!arg0) {
        return Sum();
      }
      AddANode (arg0);
      return new _Constant (0.0);
    case HY_OP_CODE_SUB:
      if (!arg0) {
        return new _MathObject;
      }
      RemoveANode (arg0);
      return new _Constant (0.0);
  }
  
  if (arg0) {
    switch (opCode) { // operations that require exactly one argument
      case HY_OP_CODE_MUL: // compute the strict consensus between T1 and T2
        return SplitsIdentity (arg0);
        
      case HY_OP_CODE_LEQ: { // MatchPattern (<=)
        
        if (arg0->ObjectClass()!=TREE && arg0->ObjectClass()!=TOPOLOGY) {
          context->ReportError ("Invalid (not a tree/topology) 2nd argument is call to <= for trees/topologies.");
          return new _MathObject;
        }
        _String  res (((_TreeTopology*)arg0)->MatchTreePattern (this));
        return new _Constant (!res.beginswith ("Unequal"));
      }
      case HY_OP_CODE_EQ: // ==
        return  Compare(arg0);
      case HY_OP_CODE_BRANCHLENGTH: //BranchLength
        return BranchLength(arg0);
      case HY_OP_CODE_BRANCHNAME: //BranchName
        return TreeBranchName(arg0);
      case HY_OP_CODE_TIPNAME: //BranchName
        return TipName(arg0);
      case HY_OP_CODE_MIN: // COT (Min)
        return FindCOT (arg0);
      case HY_OP_CODE_REROOTTREE: // RerootTree
        return RerootTree(arg0);
      case HY_OP_CODE_POWER: //^
        return AVLRepresentation (arg0);
      case HY_OP_CODE_IDIV: { // Split ($) - 2nd argument
        if (arg0->ObjectClass()!=NUMBER) {
          context->ReportError ("Invalid (not a number) 2nd argument is call to $ for trees.");
          return new _MathObject;
        }
        _Constant*  cc     = (_Constant*)TipCount();
        long        size   = cc->Value()/arg0->Value();
        
        if  ((size<=4)||(size>cc->Value()/2)) {
          context->ReportError ("Poor choice of the 2nd numeric agrument in to $ for tree. Either the resulting cluster size is too big(>half of the tree), or too small (<4)!");
          return new _MathObject;
        }
        
        long        checkSize = 1,
        tol       = 0;
        
        while (tol<size-2) {
          _List*      resL   = SplitTreeIntoClusters (size,tol);
          
          checkSize = cc->Value();
          
          if (resL->lLength) {
            _Matrix*    mRes   = new _Matrix (resL->lLength, 2, false, true);
            checkPointer (mRes);
            
            for (long k = 0; k < resL->lLength; k++) {
              _List* thisList = (_List*)(*resL)(k);
              long   nL       = ((_Constant*)(*thisList)(1))->Value();
              mRes->Store (k,0, nL);
              mRes->Store (k,1, thisList->lLength-2);
              checkSize -= nL;
            }
            
            
            if (checkSize == 0) {
              DeleteObject (cc);
              _Matrix     selMatrix (1,resL->lLength,false,true);
              _List       sortedList;
              for (long k = 0; k < resL->lLength; k++) {
                _List* thisList = (_List*)(*resL)(k);
                sortedList << (_String*)(*thisList)(0);
                _FString  *choiceString = new _FString (*(_String*)(*thisList)(0));
                _Formula  sf (choiceString);
                selMatrix.MStore(0,k,sf);
              }
              sortedList.Sort();
              for (long m = 0; m < sortedList.lLength; m++) {
                _FString  *choiceString = new _FString (*(_String*)sortedList(m));
                _Formula  sf (choiceString);
                selMatrix.MStore(0,m,sf);
              }
              
              CheckReceptacle (&splitNodeNames, emptyString, false)->SetValue (&selMatrix);
              DeleteObject (resL);
              return mRes;
            }
            
            DeleteObject (mRes);
          }
          
          DeleteObject (resL);
          tol ++;
        }
        
        DeleteObject (cc);
        return new _Matrix (1,1,false, true);
      }
    }
    
    _MathObject * arg1 = _extract_argument (arguments, 1UL, false);
    
    switch (opCode) {
      case HY_OP_CODE_MACCESS: // MAccess
        return TreeBranchName (arg0,true, arg1);
    }
    
    
    if (arg1) {
      switch (opCode) {
        case HY_OP_CODE_FORMAT: { // Format
          _String  *tStr = new _String  ((unsigned long)1024,true);
          SubTreeString (theRoot, *tStr, arg0->Compute()->Value() > 0.1 , arg1->Compute()->Value() > 0.1 ? -3:-1, nil);
          tStr->Finalize();
          return new _FString (tStr);
        }
          
      }
    }
    
  }
  
  switch (opCode) {
    case HY_OP_CODE_MUL: // compute the strict consensus between T1 and T2
    case HY_OP_CODE_LEQ: // MatchPattern (<=)
    case HY_OP_CODE_EQ: // ==
    case HY_OP_CODE_BRANCHLENGTH: //BranchLength
    case HY_OP_CODE_BRANCHNAME: //BranchName
    case HY_OP_CODE_MIN: // COT (Min)
    case HY_OP_CODE_REROOTTREE: // RerootTree
    case HY_OP_CODE_POWER: //^
    case HY_OP_CODE_IDIV:  // Split ($) - 2nd argument
    case HY_OP_CODE_MACCESS: // MAccess
    case HY_OP_CODE_FORMAT:  // Format
      WarnWrongNumberOfArguments (this, opCode,context, arguments);
      break;
    default:
      WarnNotDefined (this, opCode,context);
  }
  
  return new _MathObject;
}


//__________________________________________________________________________________

_PMathObj _TheTree::ExecuteSingleOp (long opCode, _List* arguments, _hyExecutionContext* context) {

    switch (opCode) {
       case HY_OP_CODE_TYPE: // Type
        return Type();
        break;
    }

    _MathObject * arg0 = _extract_argument (arguments, 0UL, false);
  
    if (arg0) {
      switch (opCode) {
        case HY_OP_CODE_TEXTREESTRING: // TEXTreeString
          return TEXTreeString(arg0);
       }
      
      _MathObject * arg1 = _extract_argument (arguments, 1UL, false);
 
      if (arg1) {
        switch (opCode) {
          case HY_OP_CODE_PSTREESTRING: //PlainTreeString
            return PlainTreeString(arg0,arg1);
        }
       
      }
      
    }
  
    switch (opCode) {
      case HY_OP_CODE_TEXTREESTRING:
      case HY_OP_CODE_PSTREESTRING: 
        WarnWrongNumberOfArguments (this, opCode,context, arguments);
        return new _MathObject;
    }

    return  _TreeTopology::ExecuteSingleOp (opCode,arguments,context);

}

  //_______________________________________________________________________________________________

const _List     _TreeTopology::RetrieveNodeNames (bool doTips, bool doInternals, int traversalType) const {
  _List result;
  
  node_iterator<long> ni (theRoot, traversalType);
  
  while (node<long> * iterator = ni.Next()) {
    if (iterator->is_leaf() && doTips || doInternals && !iterator->is_leaf()) {
      result < new _String (GetNodeName(iterator));
    }
  }
  
  return result;
}

//__________________________________________________________________________________

void _TreeTopology::FindCOTHelper (node<long>* aNode, long parentIndex, _Matrix& distances, _Matrix& rootDistances, _Matrix& branchLengths, _List& childLists, _AVLListX& addressToIndexMap2, _Parameter d)
{
    long          myIndex     = addressToIndexMap2.GetXtra(addressToIndexMap2.Find((BaseRef)aNode)),
                  leafCount   = distances.GetVDim();
    //bool        isRoot      = parentIndex < 0;

    _SimpleList * childLeaves = (_SimpleList*)childLists(myIndex);

    _Matrix     * lookup = parentIndex>=0?&distances:&rootDistances;

    if (parentIndex < 0) {
        parentIndex = 0;
    }

    long ci2 = 0;

    _Parameter myLength = branchLengths.theData [myIndex];

    for (long ci = 0; ci < leafCount; ci++) {
        if (ci == childLeaves->lData[ci2]) {
            if (ci2 < childLeaves->lLength - 1) {
                ci2++;
            }
        } else {
            distances.Store (myIndex, ci, (*lookup)(parentIndex,ci) + myLength);
        }
    }


    for (long ci3 = aNode->get_num_nodes(); ci3; ci3--) {
        FindCOTHelper (aNode->go_down (ci3), myIndex, distances, rootDistances, branchLengths, childLists, addressToIndexMap2, myLength);
    }
}

//__________________________________________________________________________________

void _TreeTopology::FindCOTHelper2 (node<long>* aNode, _Matrix& branchSpans, _Matrix& branchLengths, _AVLListX& addressToIndexMap2, node<long>* referrer, _Parameter d)
{
    long          myIndex     = aNode->parent?addressToIndexMap2.GetXtra(addressToIndexMap2.Find((BaseRef)aNode)):-1;
    _Parameter    myLength    = myIndex>=0?branchLengths.theData [myIndex]:0.0;

    for (long ci = aNode->get_num_nodes(); ci; ci--) {
        node <long>* daChild = aNode->go_down (ci);
        if (daChild != referrer) {
            FindCOTHelper2 (daChild,  branchSpans, branchLengths, addressToIndexMap2, nil, d+myLength);
        }
    }
    if (aNode->parent) { // not the root
        if (d>=0.0) {
            branchSpans.Store (myIndex,0,d);
        } else {
            branchSpans.Store (myIndex,0,0.);
        }

        branchSpans.Store (myIndex,1,d+myLength);
        if (referrer) {
            FindCOTHelper2 (aNode->parent, branchSpans, branchLengths, addressToIndexMap2, aNode, d+myLength);
        }
    }
}

//__________________________________________________________________________________

_AssociativeList* _TreeTopology::FindCOT (_PMathObj p) {
    // Find the Center of the Tree (COT) location
    // using an L_p metric (L_2 works well)
  
    _Parameter         power           = p->Compute()->Value(),
                       totalTreeLength = 0.0;

    _AssociativeList * resList = new _AssociativeList;

    if (power<=0.) {
        WarnError (_String("Invalid power argument in call to COT finder (Min on trees). Must be positive, had :") & power);
        return resList;
    }

    _SimpleList        avlSL,
                       avlSL2,
                       listOfNodes;

    _List              avlSL3,
                       childLists;

    _AVLListX          addressToIndexMap  (&avlSL),         // leaves only
                       addressToIndexMap2 (&avlSL2),        // complete index
                       lengthToIndexMap   (&avlSL3);

    long               leafCount   = 0,
                       branchCount = 0,
                       tIndex      = 0;

    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);

    while (node<long>* iterator = ni.Next()) {
        if (iterator->is_leaf()) {
            addressToIndexMap.Insert ((BaseRef)iterator, leafCount++);
        } else {
            branchCount++;
        }

        addressToIndexMap2.Insert ((BaseRef)iterator, tIndex++);
    }

    // allocate the matrix of path lengths with hardwired (traversal order) indices
    // also allocate a list of sorted lists to store children nodes
    // and a map of (longed) node addresses to post order traversal indices

    _Matrix      distances          (branchCount+leafCount, leafCount, false, true),
                 rootDistances        (1, leafCount,  false, true),
                 branchLengths        (1, branchCount+leafCount,  false, true),
                 branchSpans      (branchCount+leafCount+1,2,false,true);

     // pass 1: fill up the nodes up to the root (i.e. below any internal node)

    ni.Reset (theRoot);
    tIndex            = 0;

    while (node<long>* iterator = ni.Next()) {
      if (iterator->is_root()) {
        break;
      }
      _Parameter          myLength = GetBranchLength     (iterator);
      lengthToIndexMap.Insert (new _String(totalTreeLength), tIndex, false, true);
      totalTreeLength      += myLength;

      branchLengths.Store (0, tIndex++, myLength);
      listOfNodes << (long)iterator;
      _SimpleList         *childIndices = new _SimpleList;

      if (iterator->is_leaf()) {
          (*childIndices) << addressToIndexMap.GetXtra(addressToIndexMap.Find((BaseRef)iterator));
      } else {
          long           myIndex = addressToIndexMap2.GetXtra(addressToIndexMap2.Find((BaseRef)iterator));
          _SimpleList    mappedLeaves (leafCount,0,0);
        
          for (long ci = iterator->get_num_nodes(); ci; ci--) {
              long          childIndex = addressToIndexMap2.GetXtra(addressToIndexMap2.Find((BaseRef)iterator->go_down (ci)));
              _SimpleList * childLeaves = (_SimpleList*)childLists(childIndex);

              myLength = branchLengths.theData[childIndex];

              for (long ci2 = 0; ci2 < childLeaves->lLength; ci2++) {
                  long ttIndex = childLeaves->lData[ci2];
                  mappedLeaves.lData[ttIndex] = 1;
                  distances.Store (myIndex, ttIndex, distances (childIndex, ttIndex) + myLength);
              }
          }

          for (long ci2 = 0; ci2 < leafCount; ci2++) {
                  if (mappedLeaves.lData[ci2]) {
                      (*childIndices) << ci2;
                  }
          }
      }

        childLists.AppendNewInstance(childIndices);
    }

    // pass 2: fill the root vector

    //nodeName = "COT_DM1";
    //setParameter (nodeName, &distances);

    for (long ci = theRoot->get_num_nodes(); ci; ci--) {
        long          childIndex = addressToIndexMap2.GetXtra(addressToIndexMap2.Find((BaseRef)theRoot->go_down (ci)));
        _SimpleList * childLeaves = (_SimpleList*)childLists(childIndex);
        _Parameter       myLength = branchLengths.theData[childIndex];
        for (long ci2 = 0; ci2 < childLeaves->lLength; ci2++) {
            tIndex = childLeaves->lData[ci2];
            rootDistances.Store (0, tIndex, distances (childIndex, tIndex) + myLength);
            //printf ("root->%s = %g\n", ((_String*)leafNames(tIndex))->sData, distances (childIndex, tIndex) + myLength);
        }
    }


    // pass 3: fill in the "other site" branch lengths

    for (long ci3 = theRoot->get_num_nodes(); ci3; ci3--) {
        FindCOTHelper (theRoot->go_down (ci3), -1, distances, rootDistances, branchLengths, childLists, addressToIndexMap2, 0);
    }

    //nodeName = "COT_DM2";
    //setParameter (nodeName, &distances);



    // now traverse the tree and look for the maximum

    tIndex                         = 0;         // stores the index of current min

    _Parameter  currentMin         = 1e100,
                currentBranchSplit = 0;


    for (long ci = distances.GetHDim()-1; ci>=0; ci--) {
        _Parameter    T           = branchLengths.theData[ci];
        _SimpleList * childLeaves = (_SimpleList*)childLists(ci);
        long          ci2         = 0;

        if (CheckEqual (power,2.0)) {
            _Parameter    sumbT  = 0.,
                          sumbT2 = 0.,
                          suma   = 0.,
                          suma2  = 0.;



            for (long ci3 = 0; ci3 < leafCount; ci3++) {
                _Parameter tt = distances(ci,ci3);

                /*
                printf ("%s->%s = %g\n", ttt2.sData, ((_String*)leafNames(ci3))->sData, tt);
                */

                if (ci3 == childLeaves->lData[ci2]) {
                    if (ci2 < childLeaves->lLength-1) {
                        ci2++;
                    }

                    suma   += tt;
                    suma2  += tt*tt;
                } else {
                    sumbT  += tt;
                    sumbT2 += tt*tt;
                }
            }


            _Parameter tt = (sumbT-suma)/leafCount;/*(sumbT-suma)/leafCount*/;
            if (tt < 0.0) {
                tt = 0.;
            } else if (tt > T) {
                tt = T;
            }

            sumbT = tt*tt*leafCount + 2*tt*(suma-sumbT) + suma2 + sumbT2;

            if (sumbT < currentMin) {
                tIndex             = ci;
                currentBranchSplit = tt;
                currentMin         = sumbT;
            }
        } else {
            _Parameter  step        = T>0.0?T*0.0001:0.1,
                        currentT    = 0.;

            while (currentT<T) {
                _Parameter dTT = 0.0;

                ci2 = 0;

                for (long ci3 = 0; ci3 < leafCount; ci3++) {
                    _Parameter tt = distances(ci,ci3);
                    if (ci3 == childLeaves->lData[ci2]) {
                        if (ci2 < childLeaves->lLength-1) {
                            ci2++;
                        }

                        dTT   += pow(tt+currentT,power);
                    } else {
                        dTT   += pow(tt-currentT,power);
                    }
                }

                if (dTT < currentMin) {
                    tIndex             = ci;
                    currentBranchSplit = currentT;
                    currentMin         = dTT;
                }
                currentT += step;
            }
        }
    }

    node <long>*    cotBranch = (node<long>*)listOfNodes.lData[tIndex];
 
    resList->MStore (cotNode,  new _FString( GetNodeName     (cotBranch),false),false);
    resList->MStore (cotSplit, new _Constant (currentBranchSplit), false);
    resList->MStore (cotDistance, new _Constant (currentMin), false);
    resList->MStore (cotBranchLength, new _Constant (branchLengths.theData[tIndex]), false);

    //  compute the distribution of lengths away from COT

    FindCOTHelper2  (cotBranch, branchSpans, branchLengths, addressToIndexMap2, nil, currentBranchSplit-branchLengths.theData [tIndex]);
    if (cotBranch->parent) {
        _Parameter adjuster = branchLengths.theData [tIndex]-currentBranchSplit;
        branchSpans.Store (branchCount+leafCount,1,adjuster);
        node <long>* cotParent = cotBranch->parent;
        if (cotParent->parent) {
            adjuster -= branchLengths.theData[addressToIndexMap2.GetXtra(addressToIndexMap2.Find((BaseRef)cotParent))];
        }
        FindCOTHelper2  (cotParent, branchSpans, branchLengths, addressToIndexMap2, cotBranch, adjuster);
    }

    _List        timeSplits;
    _AVLListX    timeSplitsAVL (&timeSplits);

    for (long sc=0; sc <= branchCount+leafCount; sc++) {
        char       buffer[256];
        snprintf  (buffer, 256, "%.15f", branchSpans(sc,0));
        _String nodeName = buffer;
        branchSpans.Store(sc,0,nodeName.toNum());
        timeSplitsAVL.Insert (nodeName.makeDynamic(),0,false,true);
        snprintf  (buffer, 256, "%.15f", branchSpans(sc,1));
        nodeName = buffer;
        branchSpans.Store(sc,1,nodeName.toNum());
        timeSplitsAVL.Insert (nodeName.makeDynamic(),0,false,true);
    }

    _Matrix       cotCDFPoints (timeSplitsAVL.countitems(),3,false,true);

    _AssociativeList  * ctl = (_AssociativeList  *)checkPointer(new _AssociativeList ());

    _SimpleList tcache;

    long        iv,
                k = timeSplitsAVL.Traverser (tcache, iv, timeSplitsAVL.GetRoot());

    for (long pc = 0; k>=0; k = timeSplitsAVL.Traverser (tcache, iv), pc++) {
        timeSplitsAVL.SetXtra (k, pc);
        cotCDFPoints.Store (pc,0,((_String*)(*((_List*)timeSplitsAVL.dataList))(k))->toNum());
    }


    for (long mxc=0; mxc <= branchCount+leafCount; mxc++) {
        _Parameter T0 =  branchSpans(mxc,0),
                   T1 =  branchSpans(mxc,1);

        if (mxc<branchCount+leafCount) {
            _String nodeName = GetNodeName     ((node<long>*)listOfNodes(mxc));
            ctl->MStore (nodeName, new _Constant (T1), false);
        }

        char       buffer[256];
        snprintf  (buffer, 256, "%.15f", T0);
        _String nn = buffer;
        tcache.Clear();
        long       startingPos =  timeSplitsAVL.Find (&nn,tcache),
                   k           = timeSplitsAVL.Next (startingPos,tcache);

        for (long pc = timeSplitsAVL.GetXtra (k); k>=0; k = timeSplitsAVL.Next (k,tcache), pc++) {
            _Parameter ub = ((_String*)(*((_List*)timeSplitsAVL.dataList))(k))->toNum();
            if (ub < T1 || CheckEqual (T1, ub)) {
                cotCDFPoints.Store (pc,1,cotCDFPoints(pc,1)+1);
            } else {
                break;
            }
        }
    }

    for (long mxc=1; mxc<cotCDFPoints.GetHDim() ; mxc++) {
        cotCDFPoints.Store (mxc,1,cotCDFPoints(mxc,1)*(cotCDFPoints(mxc,0)-cotCDFPoints(mxc-1,0))/totalTreeLength);
        cotCDFPoints.Store (mxc,2,cotCDFPoints(mxc,1)+cotCDFPoints(mxc-1,2));
    }
    timeSplitsAVL.Clear(true);


    resList->MStore (cotCDF, &cotCDFPoints, true);
    resList->MStore (cotToNode, ctl, false);

    //  sample  random branch placement
    if (totalTreeLength > 0.0) {
        _Parameter     sampler;
        checkParameter (cotSamples, sampler, 0.0);

        tIndex = sampler;
        if (tIndex >= 1) {
            _Matrix sampledDs  (tIndex, 1, false, true);

            for (long its = 0; its < tIndex; its ++) {
                _Parameter tSample = (genrand_real2 () * totalTreeLength);
                _String nn = tSample;
                long branchIndex = 0;
                if (lengthToIndexMap.FindBest (&nn,branchIndex)<=0) {
                    branchIndex --;
                }

                _Parameter    T         = branchLengths.theData[branchIndex];
                _SimpleList * childLeaves = (_SimpleList*)childLists(branchIndex);

                _Parameter dTT = 0.0;
                tSample -= ((_String*)(*(_List*)lengthToIndexMap.dataList)(branchIndex))->toNum();

                long          ci2 = 0;
                for (long ci3 = 0; ci3 < leafCount; ci3++) {
                    _Parameter tt = distances(branchIndex,ci3);
                    if (ci3 == childLeaves->lData[ci2]) {
                        if (ci2 < childLeaves->lLength-1) {
                            ci2++;
                        }

                        dTT   += pow(tt+tSample,power);
                    } else {
                        dTT   += pow(tt+T-tSample,power);
                    }
                }

                sampledDs.Store (its,0,dTT);
            }
            setParameter (cotSampler, &sampledDs);
        }

    }
    return resList;
}

//__________________________________________________________________________________

_FString*    _TreeTopology::Compare (_PMathObj p)
// compare tree topologies
{
    _FString * res = new _FString;

    long objClass = p->ObjectClass();

    if (objClass==TREE || objClass==TOPOLOGY) {
        _String cmp = CompareTrees ((_TreeTopology*)p);
        if (cmp.startswith (eqWithReroot)) {
            (*res->theString) = cmp.Cut(eqWithReroot.sLength + ((_TreeTopology*)p)->GetName()->sLength + 1, cmp.sLength-2);
        } else if (cmp.startswith(eqWithoutReroot)) {
            (*res->theString) = _String (' ');
        }
    }
    return res;
}


//__________________________________________________________________________________

char     _TreeTopology::internalNodeCompare (node<long>* n1, node<long>* n2, _SimpleList& subTreeMap, _SimpleList* reindexer, bool cangoup, long totalSize, node<long>* n22, _TreeTopology const* tree2, bool isPattern) const
// compares whether the nodes create the same subpartition
{
    // count the number of children
    long  nc1 = n1->get_num_nodes(),
          nc2 = n2->get_num_nodes() + (cangoup&&n2->parent) - (n22&&(!n2->parent));


    if ((nc1 == nc2)||(isPattern&&(nc2<nc1))) {
        // now see if the descendants in all directions can be matched
        // first prepare the list of subtrees in the 2nd tree

        _List           nodeMap;

        _SimpleList*    complement = nil;
        _SimpleList     stSizes;

        if ((cangoup||n22)&&n2->parent) {
            complement = new _SimpleList ((unsigned long)totalSize);
            checkPointer (complement);
            complement->lLength = totalSize;

            for (long k3 = 0; k3 < totalSize; k3++) {
                complement->lData[k3] = 1;
            }
        }

        nc2 = n2->get_num_nodes();

        long  skippedChild = -1,
              complementCorrection = 0;

      for (long k=1; k <= nc2; k++) {
        node<long>* meChild = n2->go_down(k);
        _SimpleList * childLeaves   = (_SimpleList*)meChild->in_object;
        if (meChild == n22) {
          if (complement) {
            if (reindexer)
              for (long k2 = 0; k2 < childLeaves->lLength; k2++) {
                long lidx = reindexer->lData[childLeaves->lData[k2]];
                if (lidx>=0) {
                  complement->lData[lidx]    = 0;
                } else {
                  complementCorrection ++;
                }
              }
            else
              for (long k2 = 0; k2 < childLeaves->lLength; k2++) {
                complement->lData[childLeaves->lData[k2]] = 0;
              }
          }
          
          skippedChild = k-1;
        } else {
                _SimpleList * subTreeLeaves = new _SimpleList ((unsigned long)totalSize);
                checkPointer (subTreeLeaves);
                subTreeLeaves->lLength = totalSize;

                //subTreeMap << ((_SimpleList*)n1->go_down(k)->in_object)->lLength;

                if (complement) {
                    if (reindexer) {
                        for (long k2 = 0; k2 < childLeaves->lLength; k2++) {
                            long lidx = reindexer->lData[childLeaves->lData[k2]];
                            if (lidx < 0) {
                                delete complement;
                                delete subTreeLeaves;
                                return 0;
                            }
                            subTreeLeaves->lData[lidx] = 1;
                            complement->lData[lidx]    = 0;
                        }
                    } else {
                        for (long k2 = 0; k2 < childLeaves->lLength; k2++) {
                            subTreeLeaves->lData[childLeaves->lData[k2]] = 1;
                            complement->lData[childLeaves->lData[k2]] = 0;
                        }
                    }
                } else if (reindexer) {
                    for (long k2 = 0; k2 < childLeaves->lLength; k2++) {
                        long lidx = reindexer->lData[childLeaves->lData[k2]];
                        if (lidx < 0) {
                            delete subTreeLeaves;
                            return 0;
                        }
                        subTreeLeaves->lData[lidx] = 1;
                    }
                } else {
                    for (long k2 = 0; k2 < childLeaves->lLength; k2++) {
                        subTreeLeaves->lData[childLeaves->lData[k2]] = 1;
                    }
                }

                stSizes << childLeaves->lLength;

                nodeMap << subTreeLeaves;

                DeleteObject (subTreeLeaves);
            }
        }

        if (complement) {
            nodeMap << complement;
            DeleteObject (complement);
            stSizes << totalSize;

            for (long k4 = 0; k4 < stSizes.lLength-1; k4++) {
                stSizes.lData[stSizes.lLength-1] -= stSizes.lData[k4];
            }

            if (n22) {
                stSizes.lData[stSizes.lLength-1] -= ((_SimpleList*)n22->in_object)->lLength-complementCorrection;
            }
        }

        // now go through the subtrees of the first tree and match up (if we can)
        // with the 2nd 1

        _SimpleList     unmatchedPatterns;

        for (long k5=1; k5 <= nc1; k5++) {
            _SimpleList * childNodes = (_SimpleList*) n1->go_down(k5)->in_object;
            long k6;
            for (k6=0; k6 < stSizes.lLength; k6++)
                if (stSizes.lData[k6] == childNodes->lLength)
                    // potential subtree match
                {
                    _SimpleList* potMap = (_SimpleList*)nodeMap(k6);
                    long k7;
                    for (k7=0; k7<childNodes->lLength; k7++)
                        if (potMap->lData[childNodes->lData[k7]] == 0) {
                            break;
                        }

                    if (k7 == childNodes->lLength) {
                        stSizes.lData[k6] = -1;

                        if (complement && (k6 == stSizes.lLength-1)) {
                            subTreeMap << -1;
                        } else {
                            if (n22&&(k6>=skippedChild)) {
                                subTreeMap << k6+1;
                            } else {
                                subTreeMap << k6;
                            }
                        }
                        break;
                    }
                }

                if (k6 == stSizes.lLength) {
                  if (isPattern) {
                    unmatchedPatterns << k5;
                    subTreeMap << -2;
                  } else {
                    return 0;
                  }
                }
        }

        if (isPattern&&unmatchedPatterns.lLength) {
            _SimpleList rematchedPatterns;
            for (long k7 = 0; k7 < unmatchedPatterns.lLength; k7++) {
                rematchedPatterns << -1;
            }

            for (long k8 = 0; k8 < stSizes.lLength; k8++) {
                if (stSizes.lData[k8]>0) {
                    _SimpleList* potMap = (_SimpleList*)nodeMap(k8);
                    for (long k9 = 0; k9 < unmatchedPatterns.lLength; k9++) {
                        if (rematchedPatterns.lData[k9] < 0) {
                            _SimpleList * childNodes = (_SimpleList*) n1->go_down(unmatchedPatterns.lData[k9])->in_object;
                            long k10 = 0;
                            for (; k10 < childNodes->lLength; k10++)
                                if (potMap->lData[childNodes->lData[k10]] == 0) {
                                    break;
                                }
                            if (k10 == childNodes->lLength) {
                                rematchedPatterns.lData[k9] = k8;
                                if (complement && (k8 == stSizes.lLength-1)) {
                                    subTreeMap.lData [unmatchedPatterns.lData[k9]-1] = -(0xfffffff);
                                } else {
                                    if (n22&&(k8>=skippedChild)) {
                                        subTreeMap.lData [unmatchedPatterns.lData[k9]-1] = -k8-3;
                                    } else {
                                        subTreeMap.lData [unmatchedPatterns.lData[k9]-1] = -k8-2;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    stSizes.lData[k8] = 0;
                }
            }

            if (rematchedPatterns.Find(-1)>=0) {
                return 0;
            }

            for (long k11 = 0; k11 < rematchedPatterns.lLength; k11++) {
                stSizes.lData[rematchedPatterns.lData[k11]] -= ((_SimpleList*) n1->go_down(unmatchedPatterns.lData[k11])->in_object)->lLength;
            }

            for (long k12 = 0; k12 < stSizes.lLength; k12++)
                if (stSizes.lData[k12]) {
                    return 0;
                }
        }
        return 1;
    }

    return 0;
}


//long   itcCount = 0;
//__________________________________________________________________________________

char     _TreeTopology::internalTreeCompare (node<long>* n1, node<long>* n2, _SimpleList* reindexer, char compMode, long totalSize, node<long>* n22, _TreeTopology const* tree2, bool isPattern) const
// compare tree topologies
{
    if (n1->get_num_nodes() == 0) {
        return 1;
    } else {
        _SimpleList mapper;
        if (internalNodeCompare (n1,n2, mapper, reindexer, compMode, totalSize, n22, tree2, isPattern)) {
            long nc1 = n1->get_num_nodes();

            _SimpleList     furtherMatchedPatterns;
            _List           patternList;

            for (long k1=0; k1 < nc1; k1++) {
                if (mapper.lData[k1]>=0) {
                    if (internalTreeCompare (n1->go_down(k1+1),n2->go_down(mapper.lData[k1]+1), reindexer, 0, totalSize, nil, tree2, isPattern)<1) {
                        return -1;
                    }
                } else if (mapper.lData[k1] == -1) {
                    if (internalTreeCompare (n1->go_down(k1+1),n2->parent, reindexer, 0, totalSize, n2, tree2, isPattern)<1) {
                        return -1;
                    }
                } else {
                    long idx = -mapper.lData[k1]-2,
                         k;

                    _SimpleList * patched;
                    if ((k=furtherMatchedPatterns.Find(idx))<0) {
                        k = furtherMatchedPatterns.lLength;
                        furtherMatchedPatterns << idx;
                        patched = new _SimpleList;
                        checkPointer (patched);
                        patternList << patched;
                        DeleteObject (patched);
                    }

                    patched = (_SimpleList*)patternList(k);
                    (*patched) << k1;
                }
            }
            for (long k2=0; k2 < furtherMatchedPatterns.lLength; k2++) {
                node <long>* dummy = new node<long>;
                checkPointer (dummy);
                dummy->parent = n1->parent;
                _SimpleList * children = (_SimpleList*)patternList (k2),
                              * newLeaves = new _SimpleList;

                dummy->in_object = (long)newLeaves;


                for (long k3 = 0; k3 < children->lLength; k3++) {
                    node<long>* aChild = n1->go_down(children->lData[k3]+1);
                    dummy->nodes.add (aChild);
                    _SimpleList  t;
                    t.Union (*newLeaves, *(_SimpleList*)aChild->in_object);
                    newLeaves->Clear();
                    newLeaves->Duplicate(&t);
                }


                long k4 = furtherMatchedPatterns.lData[k2];
                char res = 1;

                if (newLeaves->lLength > 1) {
                    if (k4<n2->get_num_nodes()) {
                        res = internalTreeCompare(dummy, n2->go_down(k4+1),  reindexer, 0, totalSize, nil, tree2, isPattern);
                    } else {
                        res = internalTreeCompare (dummy,n2->parent, reindexer, 0, totalSize, n2, tree2, isPattern);
                    }
                }

                DeleteObject (newLeaves);
                delete (dummy);

                if (res<1) {
                    return -1;
                }

            }
        } else {
            return 0;
        }
    }
    return 1;
}

//__________________________________________________________________________________

const _String&  _TheTree::FindMaxCommonSubTree (_TheTree const*  compareTo, long& sizeVar, _List* forest) const{
    _List           myLeaves,
                    otherLeaves,
                    sharedLeaves;

    _SimpleList     indexer,
                    otherIndexer,
                    sharedLeavesIDin1,
                    sharedLeavesIDin2;

    node<long>      *myCT,
         *otherCT;

    _String         rerootAt;

    myCT            = prepTree4Comparison(myLeaves, indexer);
    otherCT         = compareTo->prepTree4Comparison(otherLeaves, otherIndexer);

    sharedLeaves.Intersect (otherLeaves,myLeaves,&sharedLeavesIDin1,&sharedLeavesIDin2);

    if (sharedLeaves.lLength>1) { // more than one common leaf
        // now we need to map shared leaves to a common indexing space

        _SimpleList    reindexer ((unsigned long)otherLeaves.lLength),
                       lidx1     ((unsigned long)myLeaves.lLength),
                       lidx2     ((unsigned long)otherLeaves.lLength),
                       ldx1,
                       ldx2;

        reindexer.lLength       = (unsigned long)otherLeaves.lLength;
        lidx1.lLength = myLeaves.lLength;
        lidx2.lLength = otherLeaves.lLength;

        for (long k=0; k<otherLeaves.lLength; k++) {
            lidx2.lData[otherIndexer.lData[k]] = k;
        }

        for (long k0=0; k0<myLeaves.lLength; k0++) {
            lidx1.lData[indexer.lData[k0]] = k0;
        }

        for (long k1=0; k1<reindexer.lLength; k1++) {
            reindexer.lData[k1] = -1;
        }

        for (long k2=0; k2<sharedLeaves.lLength; k2++) {
            reindexer.lData[lidx2.lData[sharedLeavesIDin2.lData[k2]]] = lidx1.lData[sharedLeavesIDin1.lData[k2]];
        }

        // now we map actual leaf structures to their respective leaf indices

        node_iterator<long> ni (myCT, _HY_TREE_TRAVERSAL_POSTORDER);
      

        while (node<long>* iterator = ni.Next()) {
            if (iterator->is_leaf()) {
                ldx1 << (long)iterator;
            }
        }

        ni.Reset (otherCT);
        while (node<long>* iterator = ni.Next()) {
            if (iterator->is_leaf()) {
                ldx2 << (long)iterator;
            }
        }

        // now we loop through the list of leaves and try to match them all up

        _SimpleList     matchedTops,
                        matchedSize;

        for (long k3=0; k3<sharedLeaves.lLength-1; k3++) {
            if (reindexer.lData[lidx2.lData[sharedLeavesIDin2.lData[k3]]]>=0)
                // leaf still available
            {
                node<long>*      ln1 = (node<long>*)ldx1.lData[lidx1.lData[sharedLeavesIDin1.lData[k3]]],
                                 *       ln2 = (node<long>*)ldx2.lData[lidx2.lData[sharedLeavesIDin2.lData[k3]]],
                                         *      p1  = ln1->parent,
                                                *         p2  = ln2->parent;

                char             cRes = 0;

                while ((internalTreeCompare (p1,p2,&reindexer,0,myLeaves.lLength,p2->parent?nil:ln2,compareTo) == 1)&&p1&&p2) {
                    ln1 = p1;
                    ln2 = p2;
                    p1=p1->parent;
                    p2=p2->parent;

                    cRes = 1;
                }
                if (cRes) {
                    _SimpleList* matchedLeaves = (_SimpleList*)ln2->in_object;

                    matchedTops << (long)ln1;
                    matchedSize << matchedLeaves->lLength;

                    for (long k4=0; k4<matchedLeaves->lLength; k4++) {
                        reindexer.lData[matchedLeaves->lData[k4]] = -1;
                    }
                }
            }
        }

        if (matchedSize.lLength) {
            if (forest) {
                sizeVar = 0;
                for (long k6=0; k6<matchedSize.lLength; k6++) {
                    long maxSz = 0;

                  
                    ni.Reset(myCT);
                    node<long>*   mNode  = (node<long>*)matchedTops.lData[k6];

                    while (node<long>* iterator = ni.Next()) {
                        if (!iterator->is_leaf()) {
                            if (iterator == mNode) {
                                break;
                            }
                            maxSz ++;
                        }
                    }

                    ni.Reset(theRoot);

                    while (node<long>* iterator = ni.Next()) {
                      if (!iterator->is_leaf()) {
                            if (maxSz == 0) {
                                (*forest) << LocateVar(iterator->in_object)->GetName();
                                break;
                            }
                            maxSz--;
                        }
                    }
                    sizeVar += matchedSize.lData[k6];
                }
            } else {
                long maxSz = -1,
                     maxIdx  = 0;

                for (long k5=0; k5<matchedSize.lLength; k5++)
                    if (matchedSize.lData[k5]>maxSz) {
                        maxSz = matchedSize.lData[k5];
                        maxIdx  = k5;
                    }

                sizeVar = maxSz;

                maxSz = 0;

                node<long>*   mNode  = (node<long>*)matchedTops.lData[maxIdx];

                ni.Reset (myCT);
                while (node<long>* iterator = ni.Next()) {
                    if (iterator->get_num_nodes()) {
                        if (iterator == mNode) {
                            break;
                        }
                        maxSz ++;
                    }
                }

                ni.Reset (theRoot);

              while (node<long>* iterator = ni.Next()) {
                  if (iterator->get_num_nodes()) {
                       if (maxSz == 0) {
                            return *LocateVar(iterator->in_object)->GetName();
                        }
                        maxSz--;
                    }
                }
            }
        }
    }
    return emptyString;
}


//__________________________________________________________________________________

const _String&  _TheTree::CompareSubTrees (_TheTree* compareTo, node<long>* topNode) {
    // compare tree topologies

    _List           myLeaves,
                    otherLeaves,
                    sharedLeaves;

    _SimpleList     indexer,
                    otherIndexer,
                    sharedLeavesID;

    node<long>      *myCT,
         *otherCT;

    _String         rerootAt;

    myCT            = prepTree4Comparison(myLeaves, indexer);
    otherCT         = compareTo->prepTree4Comparison(otherLeaves, otherIndexer, topNode);

    sharedLeaves.Intersect (myLeaves, otherLeaves,&sharedLeavesID);

    // first compare the inclusion for the set of leaf labels

    if (sharedLeavesID.lLength == otherLeaves.lLength) {
        _SimpleList ilist ((unsigned long)myLeaves.lLength);
        ilist.lLength = myLeaves.lLength;

        //for BCC
        for (long k2 = 0; k2 < ilist.lLength; k2++) {
            ilist.lData[k2] = -1;
        }
            //for BCC
        for (long k2 = 0; k2 < otherIndexer.lLength; k2++) {
            ilist.lData[sharedLeavesID.lData[otherIndexer.lData[k2]]] = k2;
        }

        for (long k2 = 0; k2<indexer.lLength; k2++) {
            long         lidx      = ilist.lData[indexer.lData[k2]];
            indexer.lData[k2] = lidx >=0L ? lidx : -1L;
        }

        _SimpleList *reindexer = &indexer;

        // now compare explore possible subtree matchings
        // for all internal nodes of this tree except the root

        char   compRes = 0;

        long   tCount = 1L,
               nc2 = topNode->get_num_nodes();
      
        node_iterator<long> ni (myCT, _HY_TREE_TRAVERSAL_POSTORDER);
        ni.Next();
        node<long>* iterator = ni.Next();

        while (iterator!=myCT) {
          long nc = iterator->get_num_nodes();
          if (nc == nc2) {
              long kk;
              for (kk = 1; kk <= nc; kk++) {
                  compRes = internalTreeCompare (otherCT,iterator, reindexer, 0, otherLeaves.lLength, iterator->go_down(kk),compareTo);
                  if (compRes) {
                      if (compRes == -1) {
                          iterator = myCT;
                      }
                      break;
                  }
              }
              if (kk>nc) {
                  compRes = internalTreeCompare (otherCT,iterator, reindexer, 0, otherLeaves.lLength, nil, compareTo);
                  if (compRes) {
                      if (compRes == -1) {
                          iterator = myCT;
                      }
                      break;
                  }
              } else {
                  break;
              }
          }
          tCount ++;
          iterator = ni.Next();
        }

        if (iterator != myCT) {
            ni.Reset (theRoot);
            iterator = ni.Next ();
            while (iterator != theRoot) {
                if (tCount==0) {
                    rerootAt = _String("Matched at the ") & *LocateVar (iterator->in_object)->GetName() & '.';
                    break;
                } else {
                    tCount --;
                }

                iterator = ni.Next();
            }
        } else {
            long nc = myCT->get_num_nodes();

            if (nc == nc2+1)
                for (long kk = 1; kk <= nc; kk++) {
                    compRes = internalTreeCompare (otherCT,myCT, reindexer, 0, otherLeaves.lLength, myCT->go_down(kk),compareTo);
                    if (compRes==1) {
                        break;
                    }
                }

            if (compRes == 1) {
                rerootAt = "Matched at the root.";
            }
        }

        if (!rerootAt.sLength) {
            rerootAt = "No match: Different topologies (matching label sets).";
        }
    } else {
        rerootAt = "No match: Unequal label sets.";
    }

    destroyCompTree (myCT);
    destroyCompTree (otherCT);

    return          rerootAt;
}

//__________________________________________________________________________________
void _TreeTopology::EdgeCount (long& leaves, long& internals) {
    leaves    = 0L;
    internals = 0L;
  
    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
    while (node<long>* iterator = ni.Next()) {
        if (iterator->is_leaf()) {
            leaves ++;
        } else {
            internals ++;
        }
    }

}


//__________________________________________________________________________________
_PMathObj _TreeTopology::TipCount (void)
{
    long leaves, ints;
    EdgeCount (leaves, ints);
    return new _Constant (leaves);
}

//__________________________________________________________________________________
_PMathObj _TreeTopology::BranchCount (void)
{
    long leaves, ints;
    EdgeCount (leaves, ints);
    return new _Constant (ints-1);
}

//__________________________________________________________________________________
_PMathObj _TreeTopology::FlatRepresentation (void)
{
    _SimpleList     flatTree;

    unsigned long      count = 0UL;
    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);

    while     (node<long>* iterator = ni.Next()) {
        flatTree << iterator->in_object;
        iterator->in_object = count++;
    }


    _Matrix * res = new _Matrix (1,count, false, true);
 
    ni.Reset (theRoot);
    count = 0UL;

    while     (node<long>* iterator = ni.Next()) {
        if (iterator->is_root()) {
           res->theData[count] = -1;
        } else {
           res->theData[count] = iterator->parent->in_object;
        }

        iterator->in_object = flatTree.lData[count++];
    }
    return res;
}

//__________________________________________________________________________________
_PMathObj _TreeTopology::AVLRepresentation (_PMathObj layoutOption) {

    if (layoutOption->ObjectClass () == NUMBER) {
        bool               preOrder = layoutOption->Compute()->Value()>0.5;

        _AssociativeList * masterList = (_AssociativeList * ) checkPointer(new _AssociativeList ());
        //             arrayKey;


        long         rootIndex = 0;
      
        _SimpleList  nodeList;
        _AVLListX    nodeIndexList (&nodeList);

        node_iterator<long> ni (theRoot, preOrder ? _HY_TREE_TRAVERSAL_PREORDER : _HY_TREE_TRAVERSAL_POSTORDER);
      
        while     (node<long>* iterator = ni.Next()) {
            nodeIndexList.Insert ((BaseObj*)iterator, nodeIndexList.countitems()+1L);

            if (iterator->is_root()) {
                rootIndex = nodeIndexList.countitems();
            }
         }

        ni.Reset (theRoot);
      
        while     (node<long>* iterator = ni.Next()) {
            _AssociativeList * nodeList = new _AssociativeList ();
            nodeList->MStore ("Name", new _FString (GetNodeName (iterator)), false);
            nodeList->MStore ("Length", new _Constant (GetBranchLength (iterator)), false);
            nodeList->MStore ("Depth", new _Constant (ni.Level()), false);
            if (! iterator->is_root()) {
                nodeList->MStore ("Parent", new _Constant(nodeIndexList.GetXtra (nodeIndexList.Find((BaseObj*)iterator->parent))), false);
            }

            long nCount = iterator->get_num_nodes();
            if (nCount) {
                _AssociativeList * childList = new _AssociativeList ();
                for (long k = 1; k<=nCount; k++ ) {
                    childList->MStore (_String((long)(k-1)),new _Constant(nodeIndexList.GetXtra (nodeIndexList.Find((BaseObj*)iterator->go_down(k)))) , false);
                }
                nodeList->MStore ("Children", childList, false);
            }
            masterList->MStore (_String((long)nodeIndexList.GetXtra (nodeIndexList.Find((BaseObj*)iterator))), nodeList, false);
        }

        _AssociativeList * headerList = new _AssociativeList ();

        headerList->MStore ("Name", new _FString (*GetName()), false);
        headerList->MStore ("Root", new _Constant(rootIndex), false);
        masterList->MStore ("0", headerList, false);

        return masterList;
    }
    return new _MathObject;
}

//__________________________________________________________________________________
_PMathObj _TreeTopology::TipName (_PMathObj p) {
  
    if (p&& p->ObjectClass()==NUMBER) {
        long tip_index        = p->Value(),
             count            = -1L;

        if (tip_index < 0L) {
          return new _Matrix (RetrieveNodeNames(true, false, _HY_TREE_TRAVERSAL_POSTORDER));
        }
      
        node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
        while (node<long>* iterator = ni.Next()) {
          if (iterator->is_leaf()) {
            count++;
            if (count == tip_index) {
              return new _FString (GetNodeName (iterator));
            }
          }
        }
    }
    return new _FString (emptyString);
}

//__________________________________________________________________________________
_PMathObj _TreeTopology::BranchLength (_PMathObj p) {
  _Parameter resValue = HY_INVALID_RETURN_VALUE;
  
  if (p) {
    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
    if (p->ObjectClass()==NUMBER) {
      long res        = p->Value();
      
      if (res < 0L) {
          // get ALL branch lengths
        _GrowingVector * branch_lengths = new _GrowingVector;
        
        while (node<long>* iterator = ni.Next()) {
          if (!iterator->is_root()){
            branch_lengths->Store(GetBranchLength (iterator));
          }
        }
        branch_lengths->Store (0.0); // backward compatibility with storing 0 at the root
        branch_lengths->Trim();
        branch_lengths->Transpose();
        return branch_lengths;
      } else {
          // get a specific branch length
        long count = -1L;
        while (node<long>* iterator = ni.Next()) {
          if (!iterator->is_root()) {
            if (++count == res) {
              resValue = GetBranchLength(iterator);
              break;
            }
          }
        }
      }
    } else {
      if (p->ObjectClass()==STRING) {
        _List twoIDs = ((_FString*)p->Compute())->theString->Tokenize(";");
        
        if (twoIDs.lLength == 2 || twoIDs.lLength == 1) {
          
          _String * nodes[2] = {(_String*)twoIDs.GetItem (0),
                                 twoIDs.lLength>1?(_String*)twoIDs.GetItem (1):nil};
          
          
          node<long>* node_objects[2] = {nil,nil};
          long levels[2] = {0L,0L};
          
          
          while (node<long>* iterator = ni.Next()) {
            _String nn = GetNodeName(iterator);
            for (long i = 0L; i < 2; i++) {
              if (nodes[i] && nn.Equal(nodes[i])) {
                node_objects[i] = iterator;
                levels[i] = ni.Level();
                if (node_objects[1-i]) {
                  break;
                }
              }
            }
          }
          
          if (node_objects[0] && node_objects[1]) {
            resValue = 0.;
            
            for (long i = 0L; i < 2L; i++) { // walk up to the same depth
              while (levels[1-i] < levels[i]) {
                resValue      += GetBranchLength(node_objects[i]);
                node_objects[i] = node_objects[i]->get_parent();
                levels[i]--;
              }
            }
            
            while (node_objects[0] != node_objects[1]) {
              resValue += GetBranchLength (node_objects[0]) + GetBranchLength (node_objects[1]);
              node_objects[0] = node_objects[0]->parent;
              node_objects[1] = node_objects[1]->parent;
            }
          } else if (node_objects[0]) {
            if (nodes[1]) {
              if (nodes[0]->Equal(nodes[1])) {
                resValue = 0.0;
              } else if (nodes[1]->Equal (&expectedNumberOfSubs)) {
                _String bl;
                GetBranchLength (node_objects[0], bl, true);
                if (bl.sLength) {
                  return new _FString (bl);
                }
              }
            } else {
              resValue = GetBranchLength(node_objects[0]);
            }
          }
        }
      }
    }
  }
  
  if (isnan (resValue)) {
    return new _MathObject ();
  }
  
  return new _Constant (resValue);
  
}
//__________________________________________________________________________________
_PMathObj _TreeTopology::TreeBranchName (_PMathObj p, bool subtree, _PMathObj p2) {
  _String resString;
  
  if (p) {
    if (p->ObjectClass()==NUMBER) {
      node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
     
      long argument      = p->Value(),
      count         = -1L;
      
      if (argument>=0L) { // get a specific internal node name/subtree
        while (node<long>* iterator = ni.Next()) {
          
          if (iterator->is_root()) break;
          
          if (!iterator->is_leaf()) {
            count++;
          }
          
          if (argument == count) {
            if (subtree) {
              char            mapMode  = -1;
              if (p2) {
                _String * t = (_String*)p2->Compute()->toStr();
                DetermineBranchLengthMappingMode (t,mapMode);
                DeleteObject (t);
                switch (mapMode) {
                  case 3:
                    mapMode = -1;
                    break;
                  case 1:
                    mapMode = -3;
                    break;
                  case 2:
                    mapMode = -2;
                    break;
                }
              }
              _String         st (128L, true);
              SubTreeString   (iterator, st,true,mapMode);
              st.Finalize();
              resString = st;
            } else {
              resString = GetNodeName (iterator);
            }
            break;
          }
        }
      } else {
        
        _List branch_lengths;
        
         while (node<long>* iterator = ni.Next()) {
          branch_lengths.AppendNewInstance(new _String (GetNodeName (iterator)));
        }
        return new _Matrix (branch_lengths);
      }
    } else {
      if (p->ObjectClass()==STRING) {
        _List twoIDs = ((_FString*)p->Compute())->theString->Tokenize(";");
        

        if (twoIDs.lLength == 2UL || twoIDs.lLength == 1UL) {
          
          _String * nodes[2] = {(_String*)twoIDs.GetItem(0),
                                (_String*)(twoIDs.lLength >= 1L?twoIDs.GetItem(1):nil)};
          
          
          
          if (twoIDs.lLength == 1UL) {
            _AssociativeList * resList = new _AssociativeList;
            long            masterLevel = 0L;
            

            node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_PREORDER);
            
            while (node<long>* iterator = ni.Next()) {
                //
              _String node_name = GetNodeName   (iterator);
              if (node_name == *nodes[0]) {
                masterLevel = ni.Level();
                  //resList->MStore(node_name,new _Constant (iterator->get_num_nodes()));
                do {
                  // 20151203: SLKP this will store the root name twice; commented the line above
                  resList->MStore(GetNodeName   (iterator),new _Constant (1L+(iterator->get_num_nodes()>0L)));
                  iterator = ni.Next();
                } while (iterator && ni.Level() > masterLevel);
                
                break;
              }
            }
            if (resList->avl.countitems() == 0) {
              // did not find the target node
              DeleteObject (resList);
              return new _MathObject;
            }
            return resList;
          } else {
              // this returns the sequence of nodes between node 1 and node 2
            node<long>* node_objects[2]= {nil, nil};
            long levels[2] = {0L, 0L};
            
            node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
 
            while (node<long>* iterator = ni.Next()) {
              _String nn = GetNodeName(iterator);
              for (long i = 0L; i < 2; i++) {
                if (nodes[i] && nn == nodes[i]) {
                  node_objects[i] = iterator;
                  levels[i] = ni.Level();
                  if (node_objects[1-i]) {
                    break;
                  }
                }
              }
            }
            
            if (node_objects [0] && node_objects [1]) {
              _List partial_paths [2];
              
              
              for (long i = 0L; i < 2L; i++) { // walk up to the same depth
                while (levels[1-i] < levels[i]) {
                  partial_paths[i].AppendNewInstance(new _String (GetNodeName(node_objects[i])));
                  node_objects[i] = node_objects[i]->get_parent();
                  levels[i]--;
                }
              }
              
              
              
              while (node_objects[0] != node_objects[1]) {
                for (long i = 0L; i < 2; i++) {
                  node_objects[i] = node_objects[i]->parent;
                  partial_paths[i].AppendNewInstance(new _String (GetNodeName(node_objects[i])));
               }
             }

              partial_paths[1].Flip();
              partial_paths[0] << partial_paths[1];
              return new _Matrix(partial_paths[0]);
            } else if (node_objects[0]) {
              return new _Matrix();
            }
            return new _MathObject();
          }
        }
      }
    }
  }
  return new _FString (resString);
  
}

  //__________________________________________________________________________________
void _TreeTopology::SubTreeString (node<long>* root, _String&res, bool allNames,long branchLengths,_AVLListXL* subs) const {
  long    last_level        = 0L;
  
  
  node_iterator<long> ni (root, _HY_TREE_TRAVERSAL_POSTORDER);
  
  while (node <long>* iterator = ni.Next()) {
    //printf ("[%s] %s\n", GetNodeName(root).sData, GetNodeName(iterator).sData);
    if (ni.Level() > last_level) {
      if (last_level) {
        res<<',';
      }
      res.AppendNCopies('(', ni.Level()-last_level);
    } else if (ni.Level() < last_level) {
      res.AppendNCopies(')', last_level-ni.Level());
    } else if (last_level) {
      res<<',';
    }
    
    last_level = ni.Level();

    _String node_name = GetNodeName (iterator);
    if (subs) {
      long mapIdx = subs->Find (&node_name);
      if (mapIdx >= 0) {
        node_name = *(_String*)subs->GetXtra (mapIdx);
      }
    }
    
    if (!iterator->is_root()) {
      if (allNames || (!node_name.startswith (iNodePrefix))) {
        res << node_name;
      }
      PasteBranchLength (iterator,res,branchLengths);
    }
  }
 }


//__________________________________________________________________________________
void _TreeTopology::RerootTreeInternalTraverser (node<long>* iterator, long originator, bool passedRoot, _String&res, long blOption, bool firstTime) const {
  
  //printf ("[RerootTreeInternalTraverser]%s %ld %ld\n", GetNodeName(iterator).sData, originator, passedRoot);
  
    if (passedRoot) {
        SubTreeString (iterator, res, false, blOption);
    } else {
        // move to parent now
        node<long>*     iterator_parent = iterator->get_parent();
      
        if (iterator != theRoot) { // not root yet
            res<<'(';
            long the_index_of_this_child = iterator->get_child_num();
            RerootTreeInternalTraverser (iterator_parent, the_index_of_this_child ,false,res,blOption);
          
            if (iterator_parent->get_parent()) {
            
              for (long i = 1; i<=iterator_parent->get_num_nodes(); i++) {
                  if (i!=the_index_of_this_child) {
                    res<<',';
                    SubTreeString (iterator_parent->go_down(i),res, false,blOption);
                  }
              }
             }
            res<<')';
            if (!firstTime) {
              _String node_name = GetNodeName (iterator);
              if (!node_name.startswith(iNodePrefix)) {
                res<<node_name;
              }
            }
            PasteBranchLength (iterator,res,blOption);
        } else {
            /* passing old root
               create a new root with >=2 children nodes - this node,
               and one more containing all other children (>=2 of them)
            */
            long count               = 0L,
                 root_children_count = theRoot->get_num_nodes();

            if (root_children_count > 2) {
                res << '(';
            }

            node<long>* stash_originator = nil;

            for (long k = 1; k<=theRoot->get_num_nodes(); k++) {
                if (k==originator) {
                    stash_originator = theRoot->go_down(k);
                    continue;
                }
                if (count) {
                    res<<',';
                }
                count++;
                SubTreeString (theRoot->go_down(k), res,false,blOption);
            }
          
            if (!stash_originator) {
              WarnError ("Internal error in RerootTreeInternalTraverser");
              return;
            }
            
            if (root_children_count > 2) {
                res<<')';
            }

            PasteBranchLength (stash_originator,res,blOption);
        }
    }
}


//__________________________________________________________________________________
void            _TreeTopology::PasteBranchLength (node<long>* iterator, _String& res, long branchLengths, _Parameter factor) const {
    if (branchLengths!=-1) {
        _String t;
        if (branchLengths==-2) {
            GetBranchValue (iterator,t);
        } else if (branchLengths==-3) {
            t = GetBranchLength (iterator);
        } else {
            GetBranchVarValue (iterator, t, branchLengths);
        }

        if (t.sLength) {
            res<<':';
            res<< _String (t.toNum()*factor);
        }
    }

}

//__________________________________________________________________________________
void            _TreeTopology::GetBranchLength (node<long> * n, _String& r, bool getBL) const
{
    if (getBL) {
        r = emptyString;
    } else {
        r = compExp->theData[n->in_object];
    }
}

//__________________________________________________________________________________
_Parameter            _TreeTopology::GetBranchLength (node<long> * n) const {
    return compExp->theData[n->in_object];
}

//__________________________________________________________________________________
void            _TheTree::GetBranchLength (node<long> * n, _String& r, bool getBL) const
{
    if (getBL) {
        bool    mbf;

        _Matrix *mm,
                *fv;

        RetrieveModelComponents(((_CalcNode*)(((BaseRef*)variablePtrs.lData)[n->in_object]))->GetModelIndex(), mm, fv, mbf);

        if (mm && fv) {
            r.CopyDynamicString(mm->BranchLengthExpression(fv,mbf), true);
        } else {
            r = emptyString;
        }

    } else {
        r = ((_CalcNode*)(((BaseRef*)variablePtrs.lData)[n->in_object]))->ComputeBranchLength();
    }
}

//__________________________________________________________________________________
_Parameter            _TheTree::GetBranchLength (node<long> * n) const {
    return ((_CalcNode*)(((BaseRef*)variablePtrs.lData)[n->in_object]))->ComputeBranchLength();
}

//__________________________________________________________________________________
void            _TreeTopology::GetBranchValue (node<long> *, _String& r) const{
    r = emptyString;
}

//__________________________________________________________________________________
void            _TheTree::GetBranchValue (node<long> * n, _String& r) const
{
    _Parameter t = ((_CalcNode*)(((BaseRef*)variablePtrs.lData)[n->in_object]))->Value();
    if (t != -1.) {
        r = t;
    } else {
        r = emptyString;
    }
}

//__________________________________________________________________________________
_String*            _CalcNode::GetBranchSpec (void) {
  
    _String * res = new _String (32L, true);
    *res << GetModelName();

    if (iVariables && iVariables->lLength) {
        (*res) << (res->sLength ? ',' : '{');
      

        for (unsigned long k=0UL; k < iVariables->lLength; k+=2UL) {
            if (k) {
                (*res) << ',';
            }
          
            _Variable * av = LocateVar (iVariables->lData[k]);
            if (iVariables->lData[k+1UL] >= 0L) {
                res->AppendAnAssignmentToBuffer(LocateVar (iVariables->lData[k+1UL])->GetName(),
                                                new _String (av->Value()));
            } else {
                res->AppendAnAssignmentToBuffer(av->GetName(),
                                                 new _String (av->Value()));
            }
        }
    }

    if (dVariables && dVariables->lLength) {
        for (unsigned long k=0UL; k < dVariables->lLength; k+=2UL) {
            if (dVariables->lData[k+1UL] <= 0L) {
                (*res) << (res->sLength ? ',' : '{');
              
                _Variable * av = LocateVar (dVariables->lData[k]);
                res->AppendAnAssignmentToBuffer(av->GetName(),
                                                av->GetFormulaString(),
                                                kAppendAnAssignmentToBufferFree | kAppendAnAssignmentToBufferAssignment);
                                                //true, false, true);

            }
        }
    }

    if (res->sLength) {
        (*res) << '}';
    }

    res->Finalize();
    return res;
}

//__________________________________________________________________________________
void            _TreeTopology::GetBranchVarValue (node<long> *, _String& r, long) const {
    r = emptyString;
}

//__________________________________________________________________________________
void            _TheTree::GetBranchVarValue (node<long> * n, _String& r, long idx) const {
    _CalcNode * travNode = (_CalcNode*)(((BaseRef*)variablePtrs.lData)[n->in_object]);
    long iVarValue = travNode->iVariables->FindStepping(idx,2,1);
    if (iVarValue>0)
        // user variable, but not in the model
    {
        r = _String(LocateVar (travNode->iVariables->lData[iVarValue-1])->Value());
    } else {
        _String query = _String('.') & *LocateVar(idx)->GetName();
        for (long k = 0; k<travNode->iVariables->lLength; k+=2) {
            _Variable *localVar = LocateVar(travNode->iVariables->lData[k]);
            if (localVar->GetName()->endswith (query)) {
                r = _String(localVar->Value());
                return;
            }
        }
    }
}

//__________________________________________________________________________________
_PMathObj _TreeTopology::RerootTree (_PMathObj p)
{
    _String * res = new _String ((unsigned long)256, true);

    getINodePrefix ();
  
    if (p&& p->ObjectClass()==STRING) {
        if (rooted == UNROOTED) {
            ReportWarning ("Reroot was called with an unrooted tree. Rerooting was still performed.");
        }

        _String tNodeN = (_String*)p->toStr();

        node<long>* reroot_at = FindNodeByName (&tNodeN);
      
 
        if (reroot_at) { // good node name, can reroot
            if (reroot_at->is_root()) {
                SubTreeString (theRoot, *res,0,-2);
            } else {
              (*res)<<'('; // opening tree (
              RerootTreeInternalTraverser (reroot_at, reroot_at->get_child_num(),false,*res,-2,true);
              (*res)<<',';
              SubTreeString (reroot_at, *res,0,-2);
              (*res)<<')';
            }
        }
    } else {
        _String errMsg ("Reroot Tree was passed an invalid branch argument.");
        WarnError (errMsg);
    }

    res->Finalize();
  
    //printf ("%s\n", res->sData);
  
    return new _FString (res);
}
//__________________________________________________________________________________

void    _TheTree::AlignNodes (node<nodeCoord>* theNode) const {
    long k = theNode->get_num_nodes();
    if (k) {
        theNode->in_object.v = (theNode->go_down(1)->in_object.v+theNode->go_down(k)->in_object.v)/2.0;
        theNode->in_object.h = 0;
        for (; k; k--) {
            _Parameter t = theNode->go_down(k)->in_object.h;
            if (t<theNode->in_object.h) {
                theNode->in_object.h = t;
            }
        }
        theNode->in_object.h -= TREE_H_SHIFT;
    } else {
        theNode->in_object.v = 0;
        theNode->in_object.h = 0;
    }
}
//__________________________________________________________________________________

node<nodeCoord>* _TheTree::AlignedTipsMapping (node<long>* iterator, bool first, bool respectRoot) const{
    if (first) {
        treeLayoutVert = 0.0;
        long descendants = theRoot->get_num_nodes();
        node<nodeCoord>* aRoot = new node<nodeCoord>; // reslove rootedness here
        aRoot->in_object.varRef = -1;
        aRoot->go_next(); // dumbass initialization for later
        if (rooted == UNROOTED || !respectRoot) {
            for (long k=1L; k<=descendants; k++) {
                aRoot->add_node(*AlignedTipsMapping(theRoot->go_down(k)));
            }
            AlignNodes (aRoot);
            return aRoot;
        } else {
            node<nodeCoord>* aChild = new node<nodeCoord>;
            aChild->in_object.varRef = -1;
            if (rooted == ROOTED_LEFT) {
                aRoot->add_node (*aChild);
                for (long k=1L; k<descendants; k++) {
                   aChild->add_node(*AlignedTipsMapping(theRoot->go_down(k)));
                }
                aRoot->add_node(*AlignedTipsMapping(theRoot->go_down(descendants)));
            } else {
                aRoot->add_node(*AlignedTipsMapping(theRoot->go_down(1)));
                for (long k=2L; k<=descendants; k++) {
                    aChild->add_node(*AlignedTipsMapping(theRoot->go_down(k)));
                }
                aRoot->add_node (*aChild);
            }
            AlignNodes (aChild);
            AlignNodes (aRoot);
            return      aRoot;
        }
    } else {
        node<nodeCoord>*aNode = new node<nodeCoord>;
        long descendants = iterator->get_num_nodes();
        for (long k=1L; k<=descendants; k++) {
            aNode->add_node(*AlignedTipsMapping(iterator->go_down(k)));
        }
        if (!descendants) { // terminal node
            aNode->in_object.v = treeLayoutVert;
            aNode->in_object.h = 0;
            treeLayoutVert+=TREE_V_SHIFT;
        } else {
            AlignNodes (aNode);
        }
        aNode->in_object.varRef = iterator->in_object;
        return aNode;
    }
}


//__________________________________________________________________________________

void _TheTree::ScaledBranchReMapping (node<nodeCoord>* theNode, _Parameter tw) const
{
    theNode->in_object.h -= tw;
    for (long k=1; k<=theNode->get_num_nodes(); k++) {
        ScaledBranchReMapping (theNode->go_down(k), tw);
    }
}
//__________________________________________________________________________________

_String      _TreeTopology::DetermineBranchLengthMappingMode (_String* param, char& mapMode)
{
    mapMode = 3;
    if (param) {
        if (param->Equal(&expectedNumberOfSubs)) {
            mapMode = 1;
        } else if (param->Equal(&stringSuppliedLengths)) {
            mapMode = 2;
        } else {
            mapMode = 0;
            return _String('.') & *param;
        }

    }
    return emptyString;
}
//__________________________________________________________________________________

_Parameter       _TheTree::DetermineBranchLengthGivenScalingParameter (long varRef, _String& matchString, char mapMode) const
{
    if (mapMode == 3) {
        return 1.;
    }

    _CalcNode * travNode = (_CalcNode*)LocateVar(varRef);

    _Parameter branchLength = HY_REPLACE_BAD_BRANCH_LENGTH_WITH_THIS;

    if (mapMode==1) {
        return travNode->ComputeBranchLength();
    } else if (mapMode==2) {
        branchLength = travNode->Value();
        if (branchLength<=0.0) {
            branchLength = HY_REPLACE_BAD_BRANCH_LENGTH_WITH_THIS;
        }
    } else {
        long j;
        if (travNode->iVariables)
            for (j=0; j<travNode->iVariables->lLength; j+=2) {
                _Variable* curVar  = LocateVar (travNode->iVariables->lData[j]);
                if (curVar->GetName()->endswith (matchString)) {
                    branchLength = curVar->Compute()->Value();
                    if (branchLength<=0.0) {
                        branchLength = HY_REPLACE_BAD_BRANCH_LENGTH_WITH_THIS;
                    } else {
                        break;
                    }
                }
            }

        if (((!travNode->iVariables) || j == travNode->iVariables->lLength) && travNode->dVariables)
            for (j=0; j<travNode->dVariables->lLength; j+=2) {
                _Variable* curVar = LocateVar (travNode->dVariables->lData[j]);
                if (curVar->GetName()->endswith (matchString)) {
                    branchLength = curVar->Compute()->Value();
                    if (branchLength<=0.0) {
                        branchLength = HY_REPLACE_BAD_BRANCH_LENGTH_WITH_THIS;
                    } else {
                        break;
                    }
                }
            }
    }

    return branchLength;
}


//__________________________________________________________________________________

node<nodeCoord>* _TheTree::ScaledBranchMapping (node<nodeCoord>* theParent, _String* scalingParameter, long locDepth, long& depth, char mapMode) const{
    // run a pass of aligned tip mapping then perform one more pass from the root to the children
    // pre-order to remap the length of branches.

    static  _Parameter treeWidth;
    bool     wasRoot = !theParent;

    if (!theParent) {
        theParent = AlignedTipsMapping (theRoot, true,true);
        theParent->in_object.h = 0.0;
        treeWidth = 0;
    }

    node<nodeCoord>* currentN;
    long descendants = theParent->get_num_nodes(),
         k           = 1,
         j,
         b           = -1;

    _Parameter  branchLength = HY_REPLACE_BAD_BRANCH_LENGTH_WITH_THIS;

    for  (; k<=descendants; k++) {
        currentN = theParent->go_down(k);
        j        = currentN->in_object.varRef;

        if (j>=0) {

            branchLength  = currentN->in_object.bL = DetermineBranchLengthGivenScalingParameter(j,*scalingParameter,mapMode);
            branchLength += theParent->in_object.h;

            if (branchLength>treeWidth) {
                treeWidth = branchLength;
            }

            theParent->go_down (k)->in_object.h = branchLength;
            ScaledBranchMapping (theParent->go_down(k), scalingParameter, locDepth+1, depth, mapMode);

        } else {
            theParent->go_down (k)->in_object.h = 0;
            b = k;
        }

    }

    if (k==descendants+1) {
        if (locDepth>=depth) {
            depth = locDepth+1;
        }
    }

    if (wasRoot) {
        if (b>0 && descendants==2) {
            j = (b==1)? 2 : 1;

            //branchLength = theParent->go_down(j)->in_object.h/2;
            //treeWidth -= branchLength;
            //branchLength = theParent->go_down(j)->in_object.h;
            ScaledBranchReMapping (theParent->go_down(j),0);
            theParent->go_down(b)->in_object.h = 0;
            ScaledBranchMapping (theParent->go_down(b), scalingParameter, locDepth, depth, mapMode);
        }

        ScaledBranchReMapping (theParent, treeWidth);
        return theParent;
    }
    return nil;
}

//__________________________________________________________________________________

node<nodeCoord>* _TheTree::RadialBranchMapping (node<long>* referenceNode, node<nodeCoord>* parentNode, _String* scalingParameter, _Parameter anglePerTip, long& currentTipID, _Parameter& maxRadius, char mapMode)
{
    // label 1 stores current radial distance from the root
    // label 2 stores the angle of the line to this node
    // h and v store the Cartesian coordinates

    node <nodeCoord>* current_node = new node <nodeCoord>;

    _Parameter          branchL     = 0.,
                        referenceL   = 0.;

    if  (parentNode == nil) {
        current_node->in_object.label1 = 0.0;
        current_node->in_object.label2 = 0.0;
    } else {
        referenceL = parentNode->in_object.label1;
        branchL = DetermineBranchLengthGivenScalingParameter(referenceNode->in_object,*scalingParameter,mapMode);
    }

    long                children    = referenceNode->get_num_nodes();

    current_node->in_object.label1 = referenceL + branchL;
    if (children == 0) {
        current_node->in_object.label2 = anglePerTip * currentTipID++;
        //printf ("%d %g\n",currentTipID, current_node->in_object.label2);
    } else {
        _Parameter angleSum = 0.;
        for (long n = 1; n <= children; n++) {
            node<nodeCoord>* newChild = RadialBranchMapping (referenceNode->go_down(n), current_node, scalingParameter, anglePerTip, currentTipID, maxRadius, mapMode);
            current_node->add_node(*newChild);
            angleSum += newChild->in_object.label2;
        }
        current_node->in_object.label2 = angleSum / children;
    }

    current_node->in_object.h        = current_node->in_object.label1*cos(current_node->in_object.label2);
    current_node->in_object.v        = current_node->in_object.label1*sin(current_node->in_object.label2);
    maxRadius = MAX(maxRadius, current_node->in_object.label1);
    current_node->in_object.varRef  = referenceNode->in_object;
    current_node->in_object.bL      = branchL;

    return current_node;
}



//__________________________________________________________________________________

void _TheTree::AssignLabelsToBranches (node<nodeCoord>* theParent, _String* scalingParameter, bool below)
{

    bool    wasRoot = !(theParent->parent);

    node<nodeCoord>* currentN;
    long descendants = theParent->get_num_nodes(),k=1,j,b=-1;

    _Parameter  branchLength = HY_REPLACE_BAD_BRANCH_LENGTH_WITH_THIS;

    char        mapMode;
    _String     matchString = DetermineBranchLengthMappingMode(scalingParameter, mapMode);

    for  (; k<=descendants; k++) {
        currentN = theParent->go_down(k);
        j = currentN->in_object.varRef;
        if (j>=0) {

            branchLength = DetermineBranchLengthGivenScalingParameter(j, matchString, mapMode);
            if (below) {
                currentN->in_object.label2 = branchLength;
            } else {
                currentN->in_object.label1 = branchLength;
            }

            AssignLabelsToBranches (theParent->go_down(k), scalingParameter, below);

        } else {
            if (below) {
                currentN->in_object.label2 = 0.;
            } else {
                currentN->in_object.label1 = 0.;
            }
            b = k;
            AssignLabelsToBranches (theParent->go_down(k), scalingParameter, below);
        }

    }

    if (wasRoot) {
        if ((b>0)&&(descendants==2)) {
            if (b==1) {
                j=2;
            } else {
                j=1;
            }
            if (below) {
                theParent->in_object.label2 = theParent->go_down(j)->in_object.label2/2.;
                theParent->go_down(j)->in_object.label2/=2.;
            } else {
                theParent->in_object.label1 = theParent->go_down(j)->in_object.label1/2.;
                theParent->go_down(j)->in_object.label1/=2.;
            }
        }
    }
}

//__________________________________________________________________________________

_Parameter _TheTree::PSStringWidth (_String& s)
{
    _Parameter nnWidth = 0.;
    for (long cc = 0; cc < s.sLength; cc++) {
        nnWidth += _timesCharWidths[(int)s.getChar(cc)];
    }
    return nnWidth;
}


//__________________________________________________________________________________
_PMathObj _TheTree::PlainTreeString (_PMathObj p, _PMathObj p2)
{
    _String *     res        = new _String ((unsigned long)512, true);
    long          treeLayout = 0;

    if (p&&(p->ObjectClass()==STRING)) {
        if (p2&&(p2->ObjectClass()==MATRIX)) {
            node<nodeCoord>* newRoot,
                 *currentNd;

            bool     doEmbed = false;
            bool     doSymbol = false;
            _FString *extraPS   = nil,
                      *prefixPS   = nil;
            long    symbolsize = 3;


            _AssociativeList * toptions  = (_AssociativeList*)FetchObjectFromVariableByType (&treeOutputAVL,ASSOCIATIVE_LIST);

            if (toptions) {
                _PMathObj lc = toptions->GetByKey (treeOutputLayout, NUMBER);
                if (lc) {
                    treeLayout = lc->Value();
                }
                lc = toptions->GetByKey (treeOutputEmbed, NUMBER);
                if (lc) {
                    doEmbed = lc->Value();
                }
                lc = toptions->GetByKey(treeOutputSymbols, NUMBER);
                if ( lc ) {
                    doSymbol = lc->Value();
                }

                lc = toptions->GetByKey(treeOutputExtraPS, STRING);
                if ( lc ) {
                    extraPS = (_FString*)lc->Compute();
                }

                lc = toptions->GetByKey(treeOutputPrefixPS, STRING);
                if ( lc ) {
                    prefixPS = (_FString*)lc->Compute();
                }

                lc = toptions->GetByKey(treeOutputSymbolSize, NUMBER);
                if ( lc ) {
                    symbolsize = lc->Value();
                }
            }

            _String*        theParam = (_String*) p->toStr(),
                            t;

            bool            scaling         = theParam->sLength,
                            doLabelWidth  = true;

            long            tipCount        = 0,
                            fontSize       = -1;

            _Matrix*        dimMatrix           = ((_Matrix*)p2->Compute());


            _Parameter      hScale              = 1.0,
                            vScale                = 1.0,
                            labelWidth          = 0.,

                            treeHeight           = (*dimMatrix)(0,1),
                            treeWidth         = (*dimMatrix)(0,0),

                            treeRotation      = (dimMatrix->GetVDim()>2)?(*dimMatrix)(0,2):0.0,
                            treeRadius          ,
                            treeArcAngle        ,
                            totalTreeL          = 0.,

                            mappedTreeHeight = 0.0,
                            mappedTreeHeight2 = 1e+100;

            // letter size in points
            if (treeLayout == 1) {
                treeRadius      = treeWidth;
                treeArcAngle    = treeHeight;

                if (treeRadius <= 0.0) {
                    treeRadius = 300.;
                }
                if (treeArcAngle <= 0.0 || treeArcAngle > 360.) {
                    treeArcAngle = 360.0;
                }

                treeWidth = treeHeight = treeRadius * 2.;

                treeRotation = (treeRotation - floor(treeRotation/360.)*360.)/DEGREES_PER_RADIAN;
            } else {
                if (treeHeight <= 0.0) {
                    treeHeight = 792.;
                }
                if (treeWidth  <= 0.0) {
                    treeWidth  = 576.;
                }
            }



            t = _String("%!\n%% PS file for the tree '")&*GetName()&"'.\n%% Generated by "& GetVersionString() & " on " & GetTimeStamp ();
            (*res)<<&t;
            if (treeLayout == 1) {
                (*res) << "% Radial layout\n";
                (*res) << "/righttext  {dup newpath 0 0 moveto false charpath closepath pathbbox pop exch pop exch sub       4 -1 roll exch sub 3 -1 roll newpath moveto show} def\n";
            }
            if (!doEmbed) {
                (*res)<<"<< /PageSize [";

                (*res)<<_String(treeWidth+15); /*wayne added 15 to make trees sit inside the page */
                (*res)<<' ';
                (*res)<<_String(treeHeight+15);

                (*res)<<"] >> setpagedevice\n";
            }

 
            long        xtraChars = 0;

            if (toptions) {
                _Matrix* rgbColor = (_Matrix*)(toptions)->GetByKey (treeOutputBackground, MATRIX);
                if (rgbColor) {
                    t = _String((*rgbColor)(0,0)) & " " & _String((*rgbColor)(0,1)) & " " & _String((*rgbColor)(0,2)) & " setrgbcolor\nnewpath\n";
                    (*res) << t;
                    (*res) << "0 0 moveto\n";
                    (*res) << "0 ";
                    (*res) << _String(treeHeight);
                    (*res) << " lineto\n";
                    (*res) << _String(treeWidth);
                    (*res) << ' ';
                    (*res) << _String(treeHeight);
                    (*res) << " lineto\n";
                    (*res) << _String(treeWidth);
                    (*res) << " 0 lineto\n";
                    (*res) << "closepath\nfill\n0 0 0 setrgbcolor\n";
                }

                if ( doSymbol ) {
                    /*add some symbol drawing postscript functions */
                    (*res) << "/size ";
                    (*res) << _String(symbolsize);
                    (*res) << " def\n";
                    (*res) << "/box { 0 -0.5 size mul rmoveto\n";
                    (*res) << "1 size mul 0 size mul rlineto\n";
                    (*res) << "0 size mul 1 size mul rlineto\n";
                    (*res) << "-1 size mul 0 size mul rlineto\n";
                    (*res) << "closepath\n";
                    (*res) << "} def\n";

                    (*res) << "/triangle { size size 0.5 mul rlineto 0 size -1 mul rlineto closepath } def\n";

                    (*res) << "/circle {currentpoint exch 0.5 size mul add exch 0.5 size mul 180 540 arc\n";
                    (*res) << "closepath\n";
                    (*res) << "} def\n";

                    (*res) << "/diamond { 0 -0.5 size mul rmoveto 0.5 size mul 0 rmoveto 45 rotate 0.707107 size mul 0 rlineto 0 size 0.707107 mul rlineto -0.707107 size mul 0 rlineto -45 rotate  closepath} def\n";
                }


                _Constant* fontSizeIn = (_Constant*)(toptions)->GetByKey (treeOutputFSPlaceH, NUMBER);
                if (fontSizeIn) {
                    fontSize = fontSizeIn->Value();
                }

                fontSizeIn = (_Constant*)(toptions)->GetByKey (treeOutputRightMargin, NUMBER);
                if (fontSizeIn) {
                    treeWidth = MAX(treeWidth/4+1,treeWidth-fontSizeIn->Value());
                    doLabelWidth = false;
                }

                if ((fontSizeIn = (_Constant*)(toptions)->GetByKey (treeOutputXtraMargin, NUMBER))) {
                    xtraChars = fontSizeIn->Value();
                }
            }

            char    mapMode            = 0;
            _String scalingStringMatch;

            if (scaling) {
                scalingStringMatch = DetermineBranchLengthMappingMode (theParam, mapMode);
            }

            if (treeLayout == 1) {
                _PMathObj tempValue = TipCount();
                tipCount            = tempValue->Value();
                DeleteObject        (tempValue);
                long                tipID = 0;
                hScale              = 0.;
                newRoot             = RadialBranchMapping (theRoot,nil,&scalingStringMatch,treeArcAngle*pi_const/(180.*tipCount),tipID,hScale,mapMode);
                totalTreeL          = hScale;
            } else {
                if (scaling) {
                    newRoot     = ScaledBranchMapping (nil, &scalingStringMatch, 0, tipCount,mapMode);
                } else {
                    newRoot     = AlignedTipsMapping  (theRoot, true);
                }

                hScale    = -treeWidth/newRoot->in_object.h;
            }

            currentNd = NodeTraverser (newRoot);

            while (currentNd) {
                if (currentNd->in_object.v > mappedTreeHeight) {
                    mappedTreeHeight = currentNd->in_object.v;
                }
                if (currentNd->in_object.v < mappedTreeHeight2) {
                    mappedTreeHeight2 = currentNd->in_object.v;
                }

                if (!currentNd->get_num_nodes()) {
                    tipCount++;
                }

                currentNd = NodeTraverser ((node<nodeCoord>*)nil);
            }

            // compute the size of the font, based on the spacing between branches.
            // 36 is max and 2 is min;

            if (fontSize<0) {
                fontSize = treeHeight/(tipCount+1+(scaling>0)*1.5);

                if (fontSize>36) {
                    fontSize = 36;
                } else if (fontSize<2) {
                    fontSize = 2;
                }
            }

            // now recompute the width of the of the tree, including the labels
            // do not allow names to take up more that 1/2 of the width
            // try to keep full length but if the names are too wide,
            // shrink the size of the font

            currentNd = NodeTraverser (newRoot);

            _Parameter  plotBounds[4];

            if (treeLayout == 1) {
                plotBounds [0] = plotBounds[1] = plotBounds [2] = plotBounds[3] = 0.;
            }

            while (currentNd) {
                if (currentNd->in_object.varRef>=0) {
                    _String nodeName (*LocateVar (currentNd->in_object.varRef)->GetName());
                    nodeName.Trim    (nodeName.FindBackwards ('.',0,-1)+1,-1);

                    _AssociativeList * nodeOptions = nil;
                    if (toptions) {
                        nodeOptions = (_AssociativeList *) toptions->GetByKey (nodeName, ASSOCIATIVE_LIST);
                    }


                    _PMathObj nodeLabel     = nodeOptions?nodeOptions->GetByKey (treeOutputLabel,STRING):nil,
                              nodeTLabel = nodeOptions?nodeOptions->GetByKey (treeOutputTLabel,STRING):nil;

                    _Parameter      nnWidth = 0.0;

                    if (nodeLabel) {
                        nnWidth = 0.0;
                    } else if (nodeTLabel) {
                        nnWidth = 1.+PSStringWidth (*((_FString*)nodeTLabel->Compute())->theString);
                        //printf ("%g\n", nnWidth);
                    } else if (currentNd->get_num_nodes() == 0) {
                        nnWidth = 1.+PSStringWidth (nodeName);
                    }

                    nnWidth += _maxTimesCharWidth * xtraChars;
                    nnWidth *= fontSize;

                    if (treeLayout == 1) {
                        currentNd->in_object.label2 += treeRotation;
                        _Parameter chordLength =  computeChordLength (treeRadius, currentNd->in_object.label2,plotBounds),
                                   overflow    =  MAX(0., treeRadius +
                                                      (nnWidth - chordLength) *
                                                      hScale / MAX(HY_REPLACE_BAD_BRANCH_LENGTH_WITH_THIS,currentNd->in_object.label1));

                        //nnWidth + currentNd->in_object.label1 * vScale-chordLength);



                        if (overflow > treeRadius*.5) { // shift too big
                            chordLength -= currentNd->in_object.label1/(2.*hScale) * treeRadius;

                            chordLength = chordLength / nnWidth;
                            fontSize = fontSize * chordLength;

                            if (fontSize<2) {
                                fontSize = 2;
                            }

                            nnWidth  = treeRadius*.5;
                        } else {
                            nnWidth = overflow;
                        }
                    } else {
                        if (nnWidth>treeWidth*.5) {
                            fontSize = fontSize / (2.*nnWidth/treeWidth);
                            if (fontSize<2) {
                                fontSize = 2;
                            }
                            nnWidth  = treeWidth*.5;
                        }
                    }
                    if (nnWidth > labelWidth) {
                        labelWidth = nnWidth;
                    }
                }
                currentNd = NodeTraverser ((node<nodeCoord>*)nil);
            }

            if (!doLabelWidth) {
                labelWidth = 0;
            }

            if (scaling) {
                treeHeight -= 3.0*fontSize;
            }

            if (treeLayout == 1) {
                hScale      = (treeRadius-labelWidth)/hScale;
                vScale      = 0.;
            } else {
                hScale      = -(treeWidth-labelWidth)/newRoot->in_object.h;
                vScale      = (treeHeight-2*fontSize)/(mappedTreeHeight - mappedTreeHeight2);
            }
            t = _String ("/Times-Roman findfont\n") & fontSize & " scalefont\nsetfont\n";
            (*res)<<&t;

            nodeCoord dummy;

            long   lw = fontSize/6+1;

            (*res) << _String(lw);
            (*res) << " setlinewidth\n1 setlinecap\n";

             if (prefixPS) {
                (*res) << prefixPS->theString->Replace (treeOutputFSPlaceH, _String(fontSize), true);
            }

            if (treeLayout == 1) {
                //newRoot->in_object.h = -plotBounds[1];
                //newRoot->in_object.v = -plotBounds[3];
                //hScale                 *= 2.*treeRadius / MAX(plotBounds[0]-plotBounds[1],plotBounds[2]-plotBounds[3]);
                newRoot->in_object.h = treeRadius;
                newRoot->in_object.v = treeRadius;
                vScale               = 0.;

                TreePSRecurse (newRoot, (*res), hScale, vScale, ceil(treeWidth), ceil(treeHeight), fontSize/2, labelWidth-fontSize/2, toptions, 1, &treeRadius);
            } else {
                TreePSRecurse (newRoot, (*res), hScale, vScale, ceil(treeWidth), ceil(treeHeight), fontSize/2, labelWidth-fontSize/2, toptions);
            }

            if (scaling) { /* ruler */
                if (fontSize < 8) {
                    fontSize = 8;
                    (*res) << (_String ("/Times-Roman findfont\n") & fontSize & " scalefont\nsetfont\n");
                }

                if (treeWidth < 4*fontSize) { // enforce minimal ruler width
                    treeWidth = 4*fontSize;
                }

                vScale = exp (log(10.)*floor (log10(treeLayout==1?totalTreeL:treeWidth/hScale))); // compute the scaling factor
                if (vScale*hScale < 10.0) {
                    vScale *= 10.;
                }

                _String rulerLabel (vScale, "%5.2g");

                while (vScale*hScale > (treeLayout==1?treeRadius/3:treeWidth-newRoot->in_object.h)) {
                    vScale     *= 0.5;
                    rulerLabel = vScale;
                }

                while (PSStringWidth(rulerLabel)*fontSize>vScale*hScale-3) {
                    vScale      *= 2.0;
                    rulerLabel = vScale;
                }

                long   lm = newRoot->in_object.h,
                       rm = lw+vScale*hScale;

                if (treeLayout == 1) {
                    treeHeight = 2*treeRadius - 2*fontSize;
                    lm = treeWidth - fontSize - rm;
                    rm = treeWidth - fontSize - lw;
                }

                (*res) << "newpath\n";
                t  = _String (lm) & ' ' & _String ((long)treeHeight+2*lw) & " moveto\n";
                (*res)<<&t;
                t =  _String (rm) & ' ' & (long)(treeHeight+2*lw) & " lineto\nstroke\nnewpath\n"; // main horizontal bar
                (*res)<<&t;
                t = _String (lm) & ' ' & _String ((long)treeHeight) & " moveto\n";
                (*res)<<&t;
                t =  _String (lm) & ' ' & (long)(treeHeight+4*lw) & " lineto\nstroke\nnewpath\n"; // left notch
                (*res)<<&t;
                t =  _String (rm) & ' ' & (long)(treeHeight) & " moveto\n";
                (*res)<<&t;
                t =  _String (rm) & ' ' & (long)(treeHeight+4*lw) & " lineto\nstroke\nnewpath\n"; // right notch
                (*res)<<&t;
                t =  _String (lm+lw) & ' ' & _String ((long)treeHeight+3*lw) & " moveto\n"; // main text
                (*res)<<&t;
                t =  _String ('(') & _String (rulerLabel) & ") show\n";
                (*res)<<&t;
            }

            newRoot->delete_tree ();
            delete  newRoot;

            if (extraPS) {
                (*res) << extraPS->theString->Replace (treeOutputFSPlaceH, _String(fontSize), true);
            }

            if (!doEmbed) {
                t = "showpage";
                (*res)<<&t;
            }
            DeleteObject (theParam);
        } else {
            _String errMsg ("An invalid 3rd parameter was passed to PSTreeString.");
            ReportWarning (errMsg);
        }
    } else {
        _String errMsg ("An invalid 2nd parameter was passed to PSTreeString.");
        ReportWarning (errMsg);
    }

    res->Finalize();
    return new _FString (res);
}

//_______________________________________________________________________________________________


void    _TheTree::BuildINodeDependancies (void) { // Possible deprecation?
  
  _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER);
  unsigned long        iNodeCounter  = 0UL;
  unsigned long         leafCounter = 0UL;
  
  leftiNodes.Clear ();
  topLevelNodes.Clear();
  
  while (_CalcNode * iterator = ti.Next()) {
    if (ti.IsAtLeaf()) {
      leftiNodes << iNodeCounter;
      leafCounter ++;
    } else {
      iNodeCounter++;
    }
  }
}

//_______________________________________________________________________________________________

/*

 // 20151203: SLKP seems unused (and hackish), slated for deprecation
 
void    _TheTree::BuildTopLevelCache (void) {
  
    unsigned long        iNodeCounter    = 0UL;
     unsigned long       leafCounter     = 0UL;

 
    topLevelNodes.Clear();
    topLevelLeftL.Clear();
    topLevelRightL.Clear();

    // use cBase to count the number of leaves at or below a given node
    // use categoryIndexVars to store left and right leaves

    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER);
  
    while (_CalcNode * iterator = ti.Next()) {
        if (ti.IsAtLeaf()) {
            iterator->categoryIndexVars<<leafCounter;
            iterator->categoryIndexVars<<leafCounter++;
            iterator->cBase = 1L;
        } else {
            iterator->cBase = 0L;
            node<long>* node_object = ti.GetNode();
          
            unsigned long children_count = node_object->get_num_nodes();
            for (unsigned long k = 1UL; k <= children_count; k++) {
                iterator->cBase += map_node_to_calcnode(node_object->go_down (k))->cBase;
            }
            iterator->categoryIndexVars << map_node_to_calcnode(node_object->go_down (1))->categoryIndexVars.Element(-2);
            iterator->categoryIndexVars << map_node_to_calcnode(node_object->go_down (children_count))->categoryIndexVars.Element(-1);
            iterator->lastState = iNodeCounter++;
        }
    }
  
    long threshold = (4*leafCounter)/5;

    for (unsigned long level = 1UL; level <= theRoot->nodes.length; level++) {
        node<long>* child_node_object = theRoot->go_down (level);
        _CalcNode*  child_node        = map_node_to_calcnode(child_node_object);
      
        if (child_node->cBase>1L) { // an internal node
            topLevelNodes << child_node->lastState;
            topLevelLeftL << child_node->categoryIndexVars.Element(-2L);
            topLevelRightL<< child_node->categoryIndexVars.Element  (-1L);
            if (child_node->cBase>threshold) {
                // one i-node hogging all the descendants
                _SimpleList sndLevel;
              
                for (long k = 0; k < np->nodes.length; k++) {
                    np2 = np->nodes.data[k];
                    if (np2->nodes.length) {
                        sndLevel << (long)np2;
                    }
                }
                if (sndLevel.lLength>1) {
                    topLevelLeftL.Delete (topLevelNodes.lLength-1);
                    topLevelRightL.Delete (topLevelNodes.lLength-1);
                    topLevelNodes.Delete (topLevelNodes.lLength-1);
                    for (long k=0; k<sndLevel.lLength; k++) {
                        travNode =  (_CalcNode*)LocateVar(((node<long>*)sndLevel.lData[k])->in_object);
                        topLevelNodes << travNode->lastState;
                        topLevelLeftL << travNode->categoryIndexVars[travNode->categoryIndexVars.lLength-2];
                        topLevelRightL<< travNode->categoryIndexVars[travNode->categoryIndexVars.lLength-1];
                    }
                    break;
                }
            }
        }
    }

    // restore the settings of CalcNodes.
    travNode = DepthWiseTraversal (true);
    while (travNode) {
        if (!IsCurrentNodeATip()) {
            travNode->cBase = cBase;
            travNode->lastState = -1;
        }
        long k=travNode->categoryIndexVars.lLength-2;
        travNode->categoryIndexVars.Delete (k);
        travNode->categoryIndexVars.Delete (k,true);
        travNode = DepthWiseTraversal ();
    }

    if (topLevelNodes.lLength) {
        topLevelNodes   << 0;
        topLevelLeftL   << leafCounter;
        topLevelRightL  << leafCounter-1;
    }
}
 
 */

//_______________________________________________________________________________________________


void    _TheTree::KillTopLevelCache (void)
{
    topLevelNodes.Clear();
    if (rootIChildrenCache) {
        free (rootIChildrenCache);
    }
    rootIChildrenCache = nil;
}
//_______________________________________________________________________________________________

long    _TheTree::CountTreeCategories (void) {
    categoryVariables.Clear();
    {
        _AVLList           cVA (&categoryVariables);
        ScanForCVariables (cVA);
        cVA.ReorderList   ();
    }
    categoryCount = 1L;
    for (unsigned long k=0UL; k<categoryVariables.lLength; k++) {
        categoryCount *= ((_CategoryVariable*)LocateVar(categoryVariables.lData[k]))->GetNumberOfIntervals();
    }
    return categoryCount;
}

//_______________________________________________________________________________________________

void    _TheTree::AllocateResultsCache (long size)
{
    if (rootIChildrenCache) {
        free (rootIChildrenCache);
    }
    rootIChildrenCache = nil;
    //topLevelNodes.Clear();
    size *= categoryCount;
    if (topLevelNodes.lLength) {
        rootIChildrenCache = (_Parameter*) MemAllocate (size*cBase*(topLevelNodes.lLength-1)*sizeof (_Parameter));
    }
}


//__________________________________________________________________________________

nodeCoord   _TheTree::TreeTEXRecurse (node<nodeCoord>* iterator, _String&res, _Parameter hScale, _Parameter vScale, long hSize, long vSize) const {

    _Parameter h, v;
  
    if (iterator->is_leaf()) { // terminal node
        v = vSize - iterator->in_object.v*vScale;
        h = hSize + iterator->in_object.h*hScale;
        res<< (_String ("\n\\put(")&h&','&v&"){\\circle*{2}}");
        res<< (_String ("\n\\put(")&(h+2.)&','&(v-1.)&"){\\makebox{\\tiny{");
        res<< map_node_to_calcnode(iterator)->ContextFreeName();
        res.AppendNCopies('}',3);
    } else {
        v = vSize-iterator->in_object.v*vScale;
        h = hSize + iterator->in_object.h*hScale;
      
        _Parameter leftmost_leaf_v, rightmost_leaf_v;
      
        for (long k = 1UL; k<=iterator->get_num_nodes(); k++) {
            node<nodeCoord>* child_node = iterator->go_down(k);
            TreeTEXRecurse (child_node, res, hScale, vScale, hSize, vSize);

            long        v_child = vSize - child_node->in_object.v*vScale,
                        h_child = hSize + child_node->in_object.h*hScale;

            res< (_String ("\n\\put(")&h&','&v_child&"){\\line(1,0){"&(h_child-h)&"}}");
            if (k==1UL) {
                leftmost_leaf_v  = v_child;
            } else if (k==iterator->get_num_nodes()) {
                rightmost_leaf_v = v_child;
            }
        }
        res<< (_String ("\n\\put(")&h&','&rightmost_leaf_v&"){\\line(0,1){"&(leftmost_leaf_v-rightmost_leaf_v)&"}}");
        res<< (_String ("\n\\put(")&h&','&v&"){\\circle{2}}");
      
        if ( ! iterator->is_root()) {
            res<< (_String ("\n\\put(")&(h+2.)&','&(v-1.)&"){\\makebox{\\tiny{");
          
            _String node_name = map_node_to_calcnode(iterator)->ContextFreeName();
          
            if (node_name.beginswith(iNodePrefix)) {
                node_name.Trim (iNodePrefix.Length(), -1);
            }
            res<< node_name;
            res.AppendNCopies('}',3);
        }
    }

    nodeCoord resC;
    resC.h = h;
    resC.v = v;
    return resC;
}

//__________________________________________________________________________________

void    _TheTree::TreePSRecurse (node<nodeCoord>* iterator, _String&res, _Parameter hScale, _Parameter vScale,
                                 long hSize, long vSize, long halfFontSize, long shift, _AssociativeList* outOptions,
                                 char layout, _Parameter * xtra) const
{
    unsigned long               descendants = iterator->get_num_nodes();
    bool                        is_leaf     = descendants == 0UL;
    //lineW    = halfFontSize/3+1;

    _Parameter         vc,
                       hc,
                       vcl,
                       hcl,
                       hc1,
                       hc2;

    _String            t,
                       varName,
                       colorString ("0 0 0 setrgbcolor\n");

    if (!iterator->is_root()) {
        if (iterator->in_object.varRef >=0) {
          varName = *map_node_to_calcnode(iterator)->GetName();
        }
    } else {
        hc = GetRoot().in_object;
        if (hc >= 0.) {
            varName = (*LocateVar(hc)->GetName());
        }
        if (layout == 1) {
            res <<  (_String (iterator->in_object.h) & ' ' & _String (iterator->in_object.v) & " translate\n");
        }
    }

    varName.Trim(varName.FindBackwards ('.',0L,-1L)+1L,-1L);

    _AssociativeList * nodeOptions = outOptions ? (_AssociativeList *) outOptions->GetByKey (varName, ASSOCIATIVE_LIST) : nil;

    _PMathObj nodeLabel  = nodeOptions?nodeOptions->GetByKey (treeOutputLabel,STRING):nil,
              nodeTLabel = nodeOptions?nodeOptions->GetByKey (treeOutputTLabel,STRING):nil;

    if (layout == 1) {
        hcl = iterator->in_object.h*hScale;
        vcl = iterator->in_object.v*hScale;
    } else {
        vcl = vSize-iterator->in_object.v*vScale,
        hcl = hSize+iterator->in_object.h*hScale-shift;
    }

    if (is_leaf || nodeLabel)
        // terminal node or default label
    {
        t = emptyString;
        _Parameter myAngle = layout==1?iterator->in_object.label2*DEGREES_PER_RADIAN:0.0;
        if (layout == 1) {
            res << (_String(myAngle) & " rotate\n");
            vc = 0;
            hcl = (vScale + iterator->in_object.bL)*hScale;
            hc  = hcl + halfFontSize;
        } else {
            vc = vcl-halfFontSize;
            hc = hcl+halfFontSize;
        }

        if (!nodeTLabel) {
            res<<"newpath\n";

            if (nodeLabel) {
                t = _String(hcl) & ' ' & _String (vc) & " moveto\n";

                res << &t;
                t = *((_FString*)nodeLabel->Compute())->theString;
                t = t.Replace (treeOutputNNPlaceH, varName, true);
                t = t.Replace (treeOutputFSPlaceH, _String(halfFontSize*2), true) & '\n';
            }

            if (is_leaf && t.sLength == 0)
                // generate the default label
            {
                if (layout == 1 && myAngle > 90. && myAngle < 270.) {
                    _Parameter xt = hc-halfFontSize/2,
                               yt = vc-2*halfFontSize/3;
                    t = _String (xt) & _String (" 0 translate 180 rotate 0 ") & _String (yt) & _String ('(') & varName & ") righttext 180 rotate -" & xt & " 0 translate\n";
                } else {
                    t = _String(hc-halfFontSize/2) & ' ' & _String (vc-2*halfFontSize/3) & " moveto\n";
                    res<<&t;
                    t = _String ('(') & varName & ") show\n";
                }
            }

            res<<&t;
        }

        if (is_leaf) {
            iterator->in_object.h = hc-halfFontSize;
        }

        if (layout == 1) {
            res << (_String(-iterator->in_object.label2*DEGREES_PER_RADIAN) & " rotate\n");
        }
    }

    long  minChildHC = 0x0fffffff,
          newV       = 0;

    if (!is_leaf) {
        vc = vSize - iterator->in_object.v*vScale;
        hc = hSize + iterator->in_object.h*hScale-shift;

        nodeCoord childCoord;
        for (long k = 1; k<=descendants; k++) {
            node<nodeCoord>* child = iterator->go_down(k);
            TreePSRecurse (child, res, hScale, (layout==1)?vScale+iterator->in_object.bL:vScale, hSize, vSize,halfFontSize,shift,outOptions,layout,xtra);
            if (k==1) {
                hc1 = layout==1?child->in_object.label2:child->in_object.v;
            }
            if (k==descendants) {
                hc2 = layout==1?child->in_object.label2:child->in_object.v;
            }
        }

        char doVLines = 3;

        for (long k = 1; k<=descendants; k++) {
            node<nodeCoord>* child = iterator->go_down(k);

            if (child->in_object.varRef>=0) {
                t = map_node_to_calcnode(child)->ContextFreeName();
            } else {
                t = emptyString;
            }

            newV += child->in_object.v;

            _AssociativeList * childOptions = nil;

            _Parameter         splitBranch = -1.,
                               lineWP = 0.0;

            _String            childColor,
                               notchColor,
                               *blabelString = nil,
                                childDash,
                                linewidth1,
                                linewidth2,
                                linecap1,
                                linecap2;

            _Matrix         *  notches    = nil,
                               *  multiColor = nil;

            if (layout == 1) {
                res << (_String(child->in_object.label2*DEGREES_PER_RADIAN) & " rotate\n");
            }

            if (outOptions) {
                childOptions = (_AssociativeList *) outOptions->GetByKey (t, ASSOCIATIVE_LIST);
                if (childOptions) {
                    _PMathObj keyVal = childOptions->GetByKey (treeOutputThickness,NUMBER);
                    if (keyVal) {
                        lineWP = keyVal->Compute()->Value();
                        linewidth1 = _String("currentlinewidth ") & lineWP & " setlinewidth\n";
                        linewidth2 = "setlinewidth\n";
                    }
                    keyVal = childOptions->GetByKey (treeOutputLinecap,NUMBER);
                    if (keyVal) {
                        linecap1 = _String("currentlinecap ") & (long)keyVal->Compute()->Value() & " setlinecap\n";
                        linecap2 = "setlinecap\n";
                    }
                    keyVal = childOptions->GetByKey (treeOutputSplit,NUMBER);
                    if (keyVal) {
                        splitBranch = keyVal->Compute()->Value();
                    }

                    keyVal = childOptions->GetByKey (treeOutputNotches,MATRIX);
                    if (keyVal) {
                        notches = (_Matrix*)(((_Matrix*)keyVal->Compute())->ComputeNumeric())->makeDynamic();
                    }

                    keyVal = childOptions->GetByKey (treeOutputColor,MATRIX);
                    if (keyVal) {
                        _Matrix* rgbColor = (_Matrix*)keyVal->Compute();
                        if (rgbColor->GetHDim() > 1 && rgbColor->GetVDim () == 4 && layout != 1) {
                            multiColor = (_Matrix*)rgbColor->makeDynamic();
                        } else {
                            childColor = _String((*rgbColor)(0,0)) & " " & _String((*rgbColor)(0,1)) & " " & _String((*rgbColor)(0,2)) & " setrgbcolor\n";
                        }
                    }

                    keyVal = childOptions->GetByKey (treeOutputNotchesColor,MATRIX);
                    if (keyVal) {
                        _Matrix* rgbColor = (_Matrix*)keyVal->Compute();
                        notchColor = _String((*rgbColor)(0,0)) & " " & _String((*rgbColor)(0,1)) & " " & _String((*rgbColor)(0,2)) & " setrgbcolor\n";
                    }

                    keyVal = childOptions->GetByKey (treeOutputDash,MATRIX);
                    if (keyVal) {
                        _Matrix* dash = (_Matrix*)keyVal->Compute();
                        childDash = _String('[') & _String((*dash)(0,0)) & " " & _String((*dash)(0,1)) & "] " & _String((*dash)(0,2)) & " setdash\n";
                    }
                    keyVal = childOptions->GetByKey (treeOutputOLabel,STRING);
                    if (keyVal) {
                        blabelString = ((_FString*)keyVal->Compute())->theString;
                    }
                }
            }


            if (blabelString) {
                if (layout == 1) {
                    t = _String(iterator->in_object.label1*hScale) & " 0 moveto\n";
                } else {
                    t = _String(hc) & ' ' & _String (child->in_object.v) & " moveto\n";
                }
                res<<&t;
                *blabelString = blabelString->Replace (treeOutputNNPlaceH, varName, true).Replace (treeOutputFSPlaceH, _String(halfFontSize*2), true) & '\n';

                res<<blabelString;
                res<<'\n';
            }

            res << &childDash;
            res << &childColor;
            res << &linewidth1;
            res << &linecap1;

            if (layout == 1) {
                if (child->get_num_nodes()) { // internal node
                    res << "newpath\n";
                    res << (_String("0 0 ") & child->in_object.label1*hScale & ' ' & child->in_object.v & ' ' & (child->in_object.auxD) & " arc \n");
                    res << "stroke\n";
                }
            } else {
                if (childOptions && descendants == 2) {
                    res << "newpath\n";
                    res << (_String(hc) & ' ' & _String (0.5*(hc1+hc2)) & " moveto\n");
                    res << (_String(hc) & ' ' & _String (k==1?hc1:hc2) & " lineto\n");
                    res << "stroke\n";
                    doVLines -= k;
                }

            }

            res << "newpath\n";
            if (layout == 1) {
                res << (_String(child->in_object.label1*hScale) & " 0 moveto\n");
                t = _String(-child->in_object.bL*hScale) & " 0 rlineto\n";
            } else {

                _Parameter lineWidthInset = 0.0;

                //if (lineWP > 0.0)
                //  lineWidthInset = (lineWP-lineW)*.5;

                if (multiColor) {
                    _Parameter span        = child->in_object.h - hc + 2*lineWidthInset,
                               currentX    = hc - lineWidthInset;

                    res <<  (_String(currentX) & ' ' & _String (child->in_object.v) & " moveto\n");
                    for (long seg = 0; seg < multiColor->GetHDim(); seg++) {
                        res << (_String((*multiColor)(seg,0)) & " " & _String((*multiColor)(seg,1)) & " " & _String((*multiColor)(seg,2)) & " setrgbcolor\n");
                        _Parameter mySpan = span*(*multiColor)(seg,3);
                        res << (_String(mySpan) & " 0 rlineto\n");
                        if (seg < multiColor->GetHDim()-1) {
                            currentX += mySpan;
                            res << "stroke\nnewpath\n";
                            res <<  (_String(currentX) & ' ' & _String (child->in_object.v) & " moveto\n");
                        }
                    }
                    DeleteObject (multiColor);
                    t = emptyString;
                } else {
                    res<< (_String(hc - lineWidthInset) & ' ' & _String (child->in_object.v) & " moveto\n");
                    t = _String(child->in_object.h) & ' ' & _String (child->in_object.v + lineWidthInset) & " lineto\n";
                }
                if (layout == 1) {
                    minChildHC = MIN(child->in_object.label1,minChildHC);
                } else {
                    minChildHC = MIN(child->in_object.h,minChildHC);
                }
            }
            res<<&t;
            res << "stroke\n";
            res << &linecap2;
            res << &linewidth2;

            if (childDash.sLength) {
                res << "[] 0 setdash\n";
            }

            if (splitBranch >= 0.0 && splitBranch <= 1.0) {
                res << "newpath\n";
                _Parameter x,
                           y;

                if (layout == 1) {
                    x = (child->in_object.label1 - child->in_object.bL*splitBranch)*hScale;
                    y = 0.;
                } else {
                    x = hc+(child->in_object.h-hc)*(1.-splitBranch);
                    y = child->in_object.v;
                }

                res<< (_String(x) & ' ' & _String (y) & " " & halfFontSize  & " 0 360 arc\n");
                res << "fill\n";
            }

            if (notches) {
                notches->CheckIfSparseEnough(true);
                res << notchColor;
                for (long l = 0; l < notches->GetSize(); l++) {
                    _Parameter aNotch = (*notches)[l];
                    if (aNotch >= 0. && aNotch <= 1.) {
                        res << "newpath\n";
                        _Parameter x,
                                   y;

                        if (layout == 1) {
                            x = (child->in_object.label1 - child->in_object.bL*aNotch)*hScale;
                            y = 0.;
                        } else {
                            x = hc+(child->in_object.h-hc)*(1.-aNotch);
                            y = child->in_object.v;
                        }

                        res << (_String(x-0.5*halfFontSize) & ' ' & _String (y-0.5*halfFontSize) & " moveto ");
                        res << (_String(x+0.5*halfFontSize) & ' ' & _String (y+0.5*halfFontSize) & " lineto\n");
                        res << (_String(x-0.5*halfFontSize) & ' ' & _String (y+0.5*halfFontSize) & " moveto ");
                        res << (_String(x+0.5*halfFontSize) & ' ' & _String (y-0.5*halfFontSize) & " lineto\n");
                        res << "stroke\n";
                    }
                }
                DeleteObject (notches);
            }

            res << &colorString;
            if (layout == 1) {
                res << (_String(-child->in_object.label2*DEGREES_PER_RADIAN) & " rotate\n");
            }
        }

        if (layout == 0 && doVLines) {
            _String linewidth1,
                    linewidth2,
                    linecap1,
                    linecap2;

            if (nodeOptions) {
                _PMathObj keyVal = nodeOptions->GetByKey (treeOutputThickness,NUMBER);
                if (keyVal) {
                    _Parameter lineWP = keyVal->Compute()->Value();
                    linewidth1 = _String("currentlinewidth ") & lineWP & " setlinewidth\n";
                    linewidth2 = "setlinewidth\n";
                }
                keyVal = nodeOptions->GetByKey (treeOutputLinecap,NUMBER);
                if (keyVal) {
                    linecap1 = _String("currentlinecap ") & (long)keyVal->Compute()->Value() & " setlinecap\n";
                    linecap2 = "setlinecap\n";
                }
            }

            res << &linecap1;
            res << &linewidth1;
            res << "newpath\n";

            if (doVLines == 3) {
                t = _String(hc) & ' ' & _String (hc2) & " moveto\n";
            } else {
                t = _String(hc) & ' ' & _String (0.5*(hc1+hc2)) & " moveto\n";
            }
            res<<&t;
            if (doVLines == 3) {
                t = _String(hc) & ' ' & _String (hc1) & " lineto\n";
            } else {
                t = _String(hc) & ' ' & _String (doVLines==1?hc1:hc2) & " lineto\n";
            }
            res<<&t;
            res << "stroke\n";
            res << &linewidth2;
            res << &linecap2;
        }

        if (layout == 0) {
            iterator->in_object.h = hc;
            iterator->in_object.v = newV/descendants;
        } else {
            iterator->in_object.auxD = (hc2-iterator->in_object.label2)*DEGREES_PER_RADIAN;
            if (!(iterator->is_root())) {
                iterator->in_object.v    = (hc1-iterator->in_object.label2)*DEGREES_PER_RADIAN;
            }
        }
    } else {
        iterator->in_object.v = vc;
    }

    if (nodeTLabel) {
        t = *((_FString*)nodeTLabel->Compute())->theString;
        if (t.sLength) {
            _Parameter    scF     = 2.*halfFontSize;

            if (layout == 1) {
                res << (_String(iterator->in_object.label2*DEGREES_PER_RADIAN) & " rotate\n");
            } else {
                if (!is_leaf) {
                    vcl     = iterator->in_object.v;
                }
            }


            if (scF != 2.*halfFontSize) {
                res << (_String("/Times-Roman findfont\n")& scF & " scalefont\nsetfont\n");
            }

            res << "newpath\n";
            if (layout == 1) {
                res << (_String(iterator->in_object.label1 * hScale + halfFontSize) & ' ' & _String (-2*halfFontSize/3) & " moveto\n");
            } else {
                if (!is_leaf) {
                    hc = hcl + scF*0.5;
                    vc = vcl - scF*0.33;

                } else {
                    hc -= scF*0.25;
                    vc -= scF*0.33;
                }
                res << (_String(hc) & ' ' & _String (vc) & " moveto\n");
            }

            res << '(';
            res << t;
            res << ") show \n";

            if (scF != 2.*halfFontSize) {
                res << (_String("/Times-Roman findfont\n")& halfFontSize*2 & " scalefont\nsetfont\n");
            }

            if (layout == 1) {
                res << (_String(-iterator->in_object.label2*DEGREES_PER_RADIAN) & " rotate\n");
            }
        }
    }

    if (colorString.sLength) {
        res << "0 0 0 setrgbcolor\n";
    }

    if (iterator->is_root() == nil && layout == 1) {
        res <<  (_String (-iterator->in_object.h) & ' ' & _String (-iterator->in_object.v) & " translate\n");
    }
}

//__________________________________________________________________________________
_PMathObj _TheTree::TEXTreeString (_PMathObj p) const {
    _String * res = new _String ((unsigned long)10, true);
    if (p&&(p->ObjectClass()==STRING)) {
        node<nodeCoord>*    newRoot;
        _String             *theParam = (_String*) p->toStr(),
                             t;

        bool                scaling = theParam->sLength;

        long                tipCount = 0;

        node<nodeCoord>*    currentNd;

        _Parameter          hScale = 1.0,
                            vScale = 1.0,
                            treeHeight = 0.0,
                            treeWidth;
        if (scaling) {
            newRoot = ScaledBranchMapping (nil, theParam, 0, tipCount, 0);

            treeWidth = tipCount*WIDTH_PER_BRANCH;

            if (treeWidth<MIN_TEX_WIDTH) {
                treeWidth = MIN_TEX_WIDTH;
            } else if (treeWidth>MAX_TEX_WIDTH) {
                treeWidth = MAX_TEX_WIDTH;
            }

            hScale = -treeWidth/newRoot->in_object.h;
        } else {
            newRoot = AlignedTipsMapping(theRoot, true);
            treeWidth = - newRoot->in_object.h;
            if (treeWidth<MIN_TEX_WIDTH) {
                hScale = MIN_TEX_WIDTH/treeWidth;
                treeWidth = MIN_TEX_WIDTH;
            } else if (treeWidth>MAX_TEX_WIDTH) {
                hScale = treeWidth/MAX_TEX_WIDTH;
                treeWidth = MAX_TEX_WIDTH;
            }

        }
        currentNd = newRoot;

        tipCount = newRoot->get_num_nodes();

        while (tipCount) {
            currentNd = currentNd->go_down(1);
            tipCount = currentNd->get_num_nodes();
        }

        treeHeight = currentNd->in_object.v;
        tipCount = newRoot->get_num_nodes();
        currentNd = newRoot;

        while (tipCount) {
            currentNd = currentNd->go_down(tipCount);
            tipCount = currentNd->get_num_nodes();
        }

        treeHeight = currentNd->in_object.v - treeHeight;


        if (treeHeight<MIN_TEX_HEIGHT) {
            vScale = (MIN_TEX_HEIGHT)/treeHeight;
            treeHeight = MIN_TEX_HEIGHT;
        } else if (treeHeight>MAX_TEX_HEIGHT) {
            vScale = treeHeight/(MAX_TEX_HEIGHT);
            treeHeight = MAX_TEX_HEIGHT;
        }

        t = "\n\\setlength{\\unitlength}{1mm}\n\\begin{picture}(";
        (*res)<<&t;
        t = (long)(treeWidth+5);
        (*res)<<&t;
        (*res)<<',';
        t = (long)(treeHeight+5);
        (*res)<<&t;
        (*res)<<')';

        getINodePrefix();
        TreeTEXRecurse (newRoot, (*res), hScale, vScale, ceil(treeWidth), ceil(treeHeight));
        newRoot->delete_tree ();
        delete  newRoot;

        t = "\n\\end{picture}";
        (*res)<<&t;

        DeleteObject (theParam);
    } else {
        _String errMsg ("An invalid 2nd parameter was passed to TEXTreeString");
    }

    res->Finalize();
    return new _FString (res);
}

//__________________________________________________________________________________

void _TheTree::SetUpMatrices (long categCount) {
  //fprintf (stderr, "[_TheTree::SetUpMatrices] %ld\n", categCount);
  
    categoryCount = Maximum (categCount,1L);

    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER);
 
    while   (_CalcNode* iterator = ti.Next()) {
        if (iterator->IsConstant()) {
            iterator->varFlags |= HY_VC_NO_CHECK;
        }
        iterator->ConvertToSimpleMatrix();
      
        if (categoryCount==1L) {
            iterator->matrixCache = nil;
        } else {
            iterator->matrixCache = (_Matrix**)MemAllocate (categoryCount*sizeof(_Matrix*));
            InitializeArray (iterator->matrixCache, categCount, (_Matrix*)nil);
        }
    }
}


//__________________________________________________________________________________

void _TheTree::CleanUpMatrices (void) {
  
    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER);
  
    if (categoryCount == 1) {
        while   (_CalcNode* iterator = ti.Next()) {

            // mod 05/03/2003 - uncomment next 5 lines
            // this breaks after ReplicateConstraint or MolecularClock is called
            // WTF?

            iterator->ConvertFromSimpleMatrix();

            if (iterator->referenceNode>=0) {
                iterator->SetRefNode (-1);
                iterator->compExp = nil;
            } else {
                if (iterator->referenceNode < -1) {
                    iterator->SetRefNode (-1);
                }
            }
            if (iterator->compExp) {
                DeleteObject (iterator->compExp);
                iterator->compExp = nil;
            }

            iterator->varFlags &= HY_VC_CLR_NO_CHECK;
        }
    } else {
        while   (_CalcNode* iterator = ti.Next()) {
            iterator->ConvertFromSimpleMatrix();
            if (iterator->referenceNode>=0) {
                iterator->SetRefNode (-1);
            } else
                for (long i=0; i<categoryCount; i++) {
                    DeleteObject(iterator->matrixCache[i]);
                }

            free (iterator->matrixCache);
            iterator->matrixCache = nil;
            iterator->compExp = nil;
            iterator->varFlags &= HY_VC_CLR_NO_CHECK;
        }
        categoryCount = 1;
    }

}

//__________________________________________________________________________________

void _TheTree::RemoveModel (void) {
    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER);
    while   (_CalcNode* iterator = ti.Next()) {
        iterator->RemoveModel();
    }
    categoryCount = 1;
}

  //__________________________________________________________________________________

bool _TheTree::FindScalingVariables (_SimpleList& variable_list) const {
  variable_list.Clear();
  
  _TreeIterator ti (this, _HY_TREE_TRAVERSAL_PREORDER | _HY_TREE_TRAVERSAL_SKIP_ROOT);
  _CalcNode     *iterator = ti.Next();
  
  if (iterator) {
    _SimpleList * node_variables[2] = {iterator->iVariables, iterator->dVariables};
    
    for (long k = 0L; k < 2L; k++) {
      if (node_variables[k]) {
        for (unsigned long i=1UL; i< node_variables[k]->lLength; i+=2)
          if (node_variables[k]->lData[i]>=0) {
            variable_list<< node_variables[k]->lData[i];
          }
      }
    }
    
    while (variable_list.lLength && (iterator = ti.Next())) {
      for (long i=0UL; i<variable_list.countitems(); i++) {
        if ( iterator->iVariables &&  iterator->iVariables->FindStepping(variable_list.lData[i],2,1) >= 0L  ||
             iterator->dVariables &&  iterator->dVariables->FindStepping(variable_list.lData[i],2,1) >= 0L ) {
          continue;
        }
        variable_list.Delete (i--);
      }
    }
  }
  
  return variable_list.lLength;
}

//__________________________________________________________________________________

bool _TheTree::HaveStringBranchLengths (void) const{
    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_SKIP_ROOT | _HY_TREE_TRAVERSAL_POSTORDER);
  
    while   (_CalcNode* iterator = ti.Next()) {
        if (iterator->Value() < -0.9) {
            return false;
        }
    }
    return true;
}

//__________________________________________________________________________________

void _TheTree::ScanContainerForVariables (_AVLList& l,_AVLList& l2, _AVLListX * tagger, long weight) const {
    unsigned long traversal_order = 0UL;
    _TreeIterator ti (this,  _HY_TREE_TRAVERSAL_POSTORDER);
    while   (_CalcNode* iterator = ti.Next()) {
        iterator->ScanContainerForVariables(l,l2, tagger, weight +  flatNodes.lLength + flatLeaves.lLength - traversal_order);
        traversal_order ++;
    }
}

//__________________________________________________________________________________

void _TheTree::ScanAndAttachVariables (void) const {
    _TreeIterator ti (this,  _HY_TREE_TRAVERSAL_POSTORDER);
    while   (_CalcNode* iterator = ti.Next()) {
        iterator->ScanAndAttachVariables();
    }
}

//__________________________________________________________________________________

void _TheTree::ScanForDVariables (_AVLList& l,_AVLList& l2) const {
    _TreeIterator ti (this,  _HY_TREE_TRAVERSAL_POSTORDER);
    while   (_CalcNode* iterator = ti.Next()) {
      iterator->ScanForDVariables(l,l2);
    }
}

//__________________________________________________________________________________

void _TheTree::ScanForGVariables (_AVLList& li, _AVLList& ld, _AVLListX * tagger, long weight) const {
    _SimpleList cL;
    _AVLList    cLL (&cL);

    _TreeIterator ti (this,  _HY_TREE_TRAVERSAL_POSTORDER);
    while   (_CalcNode* iterator = ti.Next()) {
    
        _Formula *explicitFormMExp = iterator->GetExplicitFormModel ();
        _Matrix  *modelM = explicitFormMExp?nil:iterator->GetModelMatrix();
        
        if ((explicitFormMExp && cLL.Find ((BaseRef)explicitFormMExp) < 0) || (modelM && cLL.Find(modelM) < 0)) {
            _SimpleList temp;
            {
                _AVLList tempA (&temp);
                if (modelM) {
                    modelM->ScanForVariables(tempA, true);
                } else {
                    explicitFormMExp->ScanFForVariables(tempA, true, false, true, true);
                }
                tempA.ReorderList();
            }
            for (unsigned long i=0; i<temp.lLength; i++) {
                long p = temp.lData[i];
                _Variable* v = LocateVar (p);
                if (v&&v->IsGlobal()) {
                    if(v->IsIndependent()) {
                        li.Insert ((BaseRef)p);
                        if (tagger) {
                            tagger->UpdateValue((BaseRef)p, weight, 0);
                        }
                    } else {
                        ld.Insert ((BaseRef)p);
                    }
                }
            }
            cLL.Insert (modelM?(BaseRef)modelM:(BaseRef)explicitFormMExp);
        }
        iterator -> ScanForGVariables(li,ld);
     }
}

//__________________________________________________________________________________

void _TheTree::ScanForCVariables (_AVLList& lcat) const {
    _TreeIterator ti (this,  _HY_TREE_TRAVERSAL_POSTORDER);
    while   (_CalcNode* iterator = ti.Next()) {
        for (unsigned long i = 0UL; i < iterator->categoryVariables.lLength; i++) {
            lcat.Insert ((BaseRef)iterator->categoryVariables.Element(i));
        }
    }
}

//__________________________________________________________________________________

bool _TheTree::HasChanged (bool) {
    _TreeIterator ti (this,  _HY_TREE_TRAVERSAL_POSTORDER);
    while   (_CalcNode* iterator = ti.Next()) {
        if (iterator->HasChanged()) {
            return true;
        }
    }
    return false;
}

  //__________________________________________________________________________________

bool _TheTree::HasChanged2 (void) {
  
  for (unsigned long k = 0; k < categoryVariables.lLength;  k++) {
    if (((_CategoryVariable*)LocateVar(categoryVariables.Element(k)))->HaveParametersChanged()) {
      return true;
    }
  }
  
  _TreeIterator ti (this,  _HY_TREE_TRAVERSAL_POSTORDER);
  while   (_CalcNode* iterator = ti.Next()) {
    if (iterator->_VariableContainer::HasChanged()) {
      return true;
    }
  }
  return false;
}


//_______________________________________________________________________________________________
void     _TheTree::InitializeTreeFrequencies (_Matrix *mx, bool setDim) {
// this will take the  matrix of frequencies and
// 1) use its dimensions to initialize tree freq holders
// 2) place global frequencies into the tree holder for later use by the pruning algo
// must be called before any tree pruning computations are started

    unsigned long vecDim = mx->GetHDim()*mx->GetVDim();
    // theModel = mx;

    if  (setDim) {
        SetTreeCodeBase (vecDim);
    } else {
        for (unsigned long i=0UL; i<vecDim; i++) {
            theProbs[i] = mx->theData[i];
        }
    }
}
  //_______________________________________________________________________________________________
void     _TheTree::SetTreeCodeBase (long b) {
  // this will take the  matrix of frequencies and
  // 1) use its dimensions to initialize tree freq holders
  // 2) place global frequencies into the tree holder for later use by the pruning algo
  // must be called before any tree pruning computations are started
  SetCodeBase (b);
  if (marginalLikelihoodCache) {
    free (marginalLikelihoodCache);
    marginalLikelihoodCache = nil;
  }
  if (cBase>0)
    marginalLikelihoodCache =
    (_Parameter*)MemAllocate ((flatNodes.lLength+flatLeaves.lLength)*sizeof (_Parameter)*cBase*systemCPUCount);
  
  _TreeIterator ti (this,  _HY_TREE_TRAVERSAL_POSTORDER);
  while   (_CalcNode* iterator = ti.Next()) {
    iterator -> SetCodeBase (b);
  }
  
}

  //_______________________________________________________________________________________________
long     _TheTree::IsLinkedToALF (long& pid) const {
  for (long lfID = 0; lfID < likeFuncList.lLength; lfID ++)
    if (likeFuncList.lData[lfID] && (pid = ((_LikelihoodFunction*)likeFuncList(lfID))->DependOnTree (*GetName())) >= 0) {
      return lfID;
    }
  return -1;
}



//_______________________________________________________________________________________________

bool     _TheTree::IntPopulateLeaves (_DataSetFilter const* dsf, long site_index, long) {
// assign proper values to leaf conditional probability vectors
    bool site_has_all_gaps = true;
  
    _String* buffer = dsf->MakeSiteBuffer();

    for (long leaf_index = 0; leaf_index<flatLeaves.lLength; leaf_index++) {
 
        _CalcNode * iterator = (_CalcNode*)flatCLeaves.GetItem(leaf_index);
        dsf->RetrieveState(site_index, leaf_index, *buffer, false);
      
      site_has_all_gaps &= ((iterator->lastState = dsf->Translate2Frequencies (*buffer, iterator->theProbs, true))<0); // ambig
      site_has_all_gaps &= (!ArrayAny (iterator->theProbs, cBase, [](_Parameter x, unsigned long) {return x == 0.0; })); //completely unresolved
      map_node_to_calcnode (((node <long>*)flatLeaves.GetElement(leaf_index))->parent)->cBase = -1;
    }

    DeleteObject (buffer);
    return site_has_all_gaps;
}


//_______________________________________________________________________________________________

void     _TheTree::RecoverNodeSupportStates (_DataSetFilter const* dsf, long index, long lastIndex, _Matrix& resultMatrix)
// assume current values of all parameters
// return 2 sets of vectors for each branch
//   - top-down  conditionals
//   - bottom-up conditionals
//   resultMatrix is assumed to contain
//      uniqueSites X (flatLeaves.lLength+flatTree.lLength)*cBase*2 X categoryCount
{

    long      globalShifter        = (flatLeaves.lLength+flatTree.lLength)*cBase,
              catShifer             = dsf->GetPatternCount() * 2 * globalShifter;

    IntPopulateLeaves (dsf, index,lastIndex);

    /* pass 1; populate top-down vectors */
    /* ugly top-bottom algorithm for debuggability and compactness */

    for (long catCount   = 0; catCount < categoryCount; catCount ++) {
        _Parameter* currentStateVector = resultMatrix.theData + 2*globalShifter*index + catShifer*catCount,
                    *   vecPointer         = currentStateVector;

        for (long nodeCount = 0; nodeCount<flatCLeaves.lLength; nodeCount++) {
            _Parameter *leafVec     = ((_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]))->theProbs;

            for (long cc = 0; cc < cBase; cc++,vecPointer++) {
                *vecPointer = leafVec[cc];
            }
        }

        for (long iNodeCount = 0; iNodeCount < flatTree.lLength; iNodeCount++) {
            node<long>* thisINode       = (node<long>*)flatNodes.lData[iNodeCount];

            for (long cc = 0; cc < cBase; cc++,vecPointer++) {
                _Parameter      tmp = 1.0;
                for (long nc = 0; nc < thisINode->nodes.length; nc++) {
                    _Parameter  tmp2 = 0.0;
                    _CalcNode   * child         = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[thisINode->nodes.data[nc]->in_object]);
                    _Parameter  * childSupport  = currentStateVector + child->nodeIndex*cBase,
                                  * transMatrix   = child->GetCompExp(categoryCount>1?catCount:(-1))->theData + cc*cBase;

                    for (long cc2 = 0; cc2 < cBase; cc2++) {
                        tmp2 += transMatrix[cc2] * childSupport[cc2];
                    }

                    tmp *= tmp2;
                }
                *vecPointer = tmp;
            }
        }
        RecoverNodeSupportStates2 (&GetRoot(),currentStateVector+globalShifter,currentStateVector,categoryCount>1?catCount:(-1));
    }
    /* pass 2; populate bottom-up vectors */
    /* for this we need to traverse the tree pre-order */
    /* because speed is not much of a concern, use a recursive call for compactness */

}

//_______________________________________________________________________________________________

void     _TheTree::RecoverNodeSupportStates2 (node<long>* thisNode, _Parameter* resultVector, _Parameter* forwardVector, long catID)
{
    _CalcNode   * thisNodeC     = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[thisNode->in_object]);
    _Parameter  * vecPointer    = resultVector + thisNodeC->nodeIndex * cBase;

    if (thisNode->parent) {
        if (thisNode->parent->parent) {
            for (long cc = 0; cc < cBase; cc++,vecPointer++) {
                _Parameter tmp = 1.0;
                for (long nc = 0; nc < thisNode->parent->nodes.length; nc++) {
                    _Parameter  tmp2            = 0.0;
                    _CalcNode   * child         = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[thisNode->parent->nodes.data[nc]->in_object]);
                    bool          invert        = (child == thisNodeC);;
                    if (invert) {
                        child = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[thisNode->parent->in_object]);
                    }

                    _Parameter  * childSupport  = invert?resultVector + cBase*child->nodeIndex
                                                  :forwardVector + child->nodeIndex*cBase,
                                                  * transMatrix   = child->GetCompExp(catID)->theData + cc*cBase;

                    for (long cc2 = 0; cc2 < cBase; cc2++) {
                        tmp2 += transMatrix[cc2] * childSupport[cc2];
                    }

                    tmp *= tmp2;
                }
                *vecPointer = tmp;
            }
        } else {
            for (long cc = 0; cc < cBase; cc++,vecPointer++) {
                _Parameter tmp = 1.0;
                for (long nc = 0; nc < thisNode->parent->nodes.length; nc++) {
                    _Parameter  tmp2            = 0.0;
                    _CalcNode   * child         = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[thisNode->parent->nodes.data[nc]->in_object]);
                    if (child != thisNodeC) {
                        _Parameter  * childSupport  = forwardVector + child->nodeIndex*cBase,
                                      * transMatrix   = child->GetCompExp(catID)->theData + cc*cBase;

                        for (long cc2 = 0; cc2 < cBase; cc2++) {
                            tmp2 += transMatrix[cc2] * childSupport[cc2];
                        }

                        tmp *= tmp2;
                    }
                }
                *vecPointer = tmp;
            }
        }
    } else
        for (long cc=0; cc<cBase; cc++) {
            vecPointer [cc] = 1.0;
        }

    for (long nc = 0; nc < thisNode->nodes.length; nc++) {
        RecoverNodeSupportStates2 (thisNode->nodes.data[nc],resultVector,forwardVector,catID);
    }
}



//_______________________________________________________________________________________________
_AVLListX*  _TheTree::ConstructNodeToIndexMap (bool doINodes) const {
    _SimpleList * nodes  = new _SimpleList;
    const _SimpleList * whichL = doINodes?&flatNodes:&flatLeaves;
    _AVLListX   * result = new _AVLListX (nodes);

    for (unsigned long   pistolero = 0; pistolero < whichL->lLength; pistolero++) {
        result->Insert ((BaseRef)whichL->lData[pistolero], pistolero, false);
    }

    return        result;

}

  //_______________________________________________________________________________________________
void _TheTree::MapPostOrderToInOderTraversal (_SimpleList& storeHere, bool doINodes) const {
  _AVLListX*          nodeMapper    = ConstructNodeToIndexMap (doINodes);
  
  _TreeIterator       ti (this, doINodes ? _HY_TREE_TRAVERSAL_PREORDER : _HY_TREE_TRAVERSAL_POSTORDER);
  
  unsigned long                allNodeCount = 0UL;
  
  storeHere.Populate (doINodes?flatTree.lLength:flatLeaves.lLength, 0, 0);
  
  while (_CalcNode* iterator = ti.Next()) {
    bool isTip = ti.IsAtLeaf();
    if ( isTip && !doINodes  || !isTip && doINodes) {
      storeHere.lData[nodeMapper->GetXtra (nodeMapper->Find((BaseRef)(ti.GetNode())))] = allNodeCount++;
    }
  }
  
  nodeMapper->DeleteAll(false);
  DeleteObject (nodeMapper);
}

  //_______________________________________________________________________________________________
void    _TheTree::MarkDone (void) {
  _TreeIterator       ti (this, _HY_TREE_TRAVERSAL_POSTORDER);
  while (_CalcNode* iterator = ti.Next()) {
    iterator -> MarkDone();
  }
}

//_______________________________________________________________________________________________

long    _TheTree::ComputeReleafingCostChar (_DataSetFilter const* dsf, long firstIndex, long secondIndex) const {

    const char *pastState = dsf->GetColumn(firstIndex),
               *thisState = dsf->GetColumn(secondIndex);

  
    _SimpleList markedNodes (flatTree.lLength, 0, 0);

    for (long nodeID = 0; nodeID<flatLeaves.lLength; nodeID++) {
        long f = dsf->theNodeMap.lData[nodeID];
        if (thisState[f] != pastState[f]) {
            markedNodes.lData[flatParents.lData[nodeID]] = 1;
        }
    }

    long theCost = 0;

    for (long i=0; i<flatTree.lLength; i++) {
        if (markedNodes.lData[i]) {
            long myParent = flatParents.lData[i + flatLeaves.lLength];
            if (myParent >= 0) {
                markedNodes.lData[myParent] = 1;
            }
            theCost += ((node <long>*)(flatNodes.lData[i]))->nodes.length;
        }
    }

    return theCost;

}

//_______________________________________________________________________________________________

void    _TheTree::ClearConstraints (void) {
    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER);
    while (_CalcNode* iterator = ti.Next()) {
        iterator->ClearConstraints();
    }
}


//_______________________________________________________________________________________________

long    _TheTree::ComputeReleafingCost (_DataSetFilter const* dsf, long firstIndex, long secondIndex, _SimpleList* traversalTags, long orderIndex) const {

    long        filterL = dsf->GetPatternCount();

    _SimpleList     markedNodes (flatTree.lLength,0,0);

    for (long leafID = 0; leafID<flatLeaves.lLength; leafID++)
        if (!dsf->CompareTwoSites(firstIndex,secondIndex,leafID)) {
            markedNodes.lData [flatParents.lData[leafID]] = 1;
        }

    // now compute the cost

    long theCost = 0;


    for (long i=0; i<flatTree.lLength; i++) {
        if (markedNodes.lData[i]) {
            long myParent    =  flatParents.lData[flatLeaves.lLength + i];
            if (myParent >= 0) {
                markedNodes.lData[myParent] = 1;
            }
            theCost         += ((node <long>*)(flatNodes.lData[i]))->nodes.length;
        } else if (traversalTags && orderIndex) {
            long theIndex = filterL * i + orderIndex;
            traversalTags->lData[theIndex/_HY_BITMASK_WIDTH_] |= bitMaskArray.masks[theIndex%_HY_BITMASK_WIDTH_];
            /*printf ("%d %d %d %d %d %d\n", firstIndex, secondIndex, orderIndex, i, theIndex/_HY_BITMASK_WIDTH_, theIndex%_HY_BITMASK_WIDTH_);
            char b[4]; b[3] = 0;
            for (long k = 0; k < flatLeaves.lLength; k++)
            {
                dsf->GrabSite (firstIndex, dsf->theNodeMap.lData[k], b);
                printf ("%s ", b);
                dsf->GrabSite (secondIndex, dsf->theNodeMap.lData[k], b);
                printf ("%s\n", b);
            }*/
        }
    }

    return theCost;

}

//_______________________________________________________________________________________________

void    _TheTree::MarkMatches (_DataSetFilter* dsf, long firstIndex, long secondIndex) const{
 
    //_CalcNode* iterator ;

    for (unsigned long n = 0UL; n<flatLeaves.lLength; n++) {
        if (!dsf->CompareTwoSites(firstIndex,secondIndex,n)) {
            map_node_to_calcnode(((node <long>*)flatLeaves.Element(n))->parent)->cBase = -1L;
        }
    }
  
    for (unsigned long n=0UL; n<flatTree.lLength; n++) {
        _CalcNode * iterator = (_CalcNode*)flatTree.GetItem(n);
        if (iterator->cBase == -1) {
            iterator = map_node_to_calcnode(((node <long>*)(flatNodes.Element (n)))->parent);
            if (iterator) {
                iterator->cBase = -1L;
            }
        }
    }
  
    for (unsigned long n=0UL; n<flatTree.lLength; n++) {
        _CalcNode * iterator = (_CalcNode*)flatTree.GetItem (n);
        if (iterator->cBase != -1L) {
            iterator->lastState = -2L;
        } else {
            iterator->cBase = cBase;
        }
    }
}

  //_______________________________________________________________________________________________

long    _TheTree::GetLowerBoundOnCost(_DataSetFilter* dsf, _SimpleList* sl) const {
  unsigned long    theCost = 0UL;
  
  _CalcNode* travNode ;
  
  for (long siteIndex = 0; siteIndex<dsf->theFrequencies.lLength; siteIndex++) {
    
    for (unsigned long n = 0UL; n<flatTree.lLength; n++) {
      ((_CalcNode*)flatTree.GetItem(n))->lastState = -1L;
    }
    
    for (unsigned long matchIndex = 0UL; matchIndex< (sl ? siteIndex : dsf->theFrequencies.lLength); matchIndex++) {
      if (matchIndex!=siteIndex) {
        MarkMatches (dsf,sl ? sl->Element (siteIndex) : siteIndex, sl ? sl->Element(matchIndex) : matchIndex);
      }
    }
    
    for (unsigned long n = 0UL; n<flatTree.lLength; n++) {
      _CalcNode * iterator = (_CalcNode*)flatTree.GetItem(n);
      if (iterator->lastState != -2) {
        theCost += ((node <long>*)(flatNodes.lData[n]))->get_num_nodes();
      }
      iterator->lastState = -1;
    }
  }
  return theCost;
}

//_______________________________________________________________________________________________

void    _TheTree::MolecularClock (_String const& baseNode, _List& varsToConstrain) const {
    node<long>* topNode = nil;
  
    if (baseNode.sLength == 0UL) { // called Molecular Clock on the entire tree
        topNode = &GetRoot();
        _String*  childNameP;
        if (rooted == ROOTED_LEFT) { // run separate constraint on the right child of the root
           MolecularClock (map_node_to_calcnode(theRoot->go_down (theRoot->get_num_nodes()))->ContextFreeName(), varsToConstrain);
          
        } else if (rooted == ROOTED_RIGHT) {
           MolecularClock (map_node_to_calcnode(theRoot->go_down (1))->ContextFreeName(), varsToConstrain);
        }
    } else {
        topNode = FindNodeByName (&baseNode);
    }

    if (!topNode) {
        WarnError (_String ("Molecular clock constraint has failed, since node '")
                   &baseNode
                   &"' is not a part of tree '"
                   &*GetName() & "'");
    } else
        for (unsigned long k=1UL; k<varsToConstrain.lLength; k++) {
            long varIndex = LocateVarByName (*(_String*)varsToConstrain (k));
          
            if (varIndex<0) {
                WarnError (_String ("Molecular clock constraint has failed, since variable' ") &*(_String*)varsToConstrain (k) &"' is undefined.");
                return ;
            }
            map_node_to_calcnode(topNode)->RecurseMC (variableNames.GetXtra(varIndex), topNode, true, rooted);
        }
}



//_______________________________________________________________________________________________

node<long>* _CalcNode::LocateMeInTree (void) const {

    _String parentName = ParentObjectName (),
            myName     = ContextFreeName();
  
    return  ((_TreeTopology*)FetchVar(LocateVarByName(parentName)))->FindNodeByName(&myName);
  
}

//_______________________________________________________________________________________________

void _CalcNode::ConvertToSimpleMatrix (void) const {
    _Formula * mf = GetExplicitFormModel();
    if (mf) {
        mf->ConvertMatrixArgumentsToSimpleOrComplexForm (false);
    } else {
        _Matrix * mm = GetModelMatrix();
        if (mm) {
            mm->MakeMeSimple();
        }

        mm = GetFreqMatrix();
        if (mm) {
            mm->MakeMeSimple();
        }
    }
}

//_______________________________________________________________________________________________

void _CalcNode::ConvertFromSimpleMatrix (void) {
    _Formula * mf = GetExplicitFormModel();
    if (mf) {
        mf->ConvertMatrixArgumentsToSimpleOrComplexForm (true);
    } else {
        _Matrix * mm = GetModelMatrix();
        if (mm) {
            mm->MakeMeGeneral();
        }

        mm = GetFreqMatrix();

        if (mm) {
            mm->MakeMeGeneral();
        }
    }
}

//_______________________________________________________________________________________________

_Formula*   _CalcNode::RecurseMC (long varToConstrain, node<long>* whereAmI, bool first, char rooted) {
    long descendants = whereAmI->get_num_nodes(),
         f = iVariables?iVariables->FindStepping(varToConstrain,2,1):-1,
         k,
         l,
         start = 0;

    if (f<0 && !first) {
        WarnError (_String ("Molecular clock constraint has failed, since variable '")
                    &*LocateVar(varToConstrain)->GetName()
                    &"' is not an independent member of the node '"
                    & *GetName()
                    & '\''
                   );
        return nil;
    }


    if (descendants == 0) {
      if (first) {
        return nil;
      } else {
        return new _Formula (LocateVar(iVariables->lData[f-1]),true);
      }
    }

    if (first && (!whereAmI->get_parent()) && (rooted == ROOTED_LEFT)) {
        descendants --;
    }
    if (first && (!whereAmI->get_parent()) && (rooted == ROOTED_RIGHT)) {
        start++;
    }

    // internal node - must do some work

    _Formula**  nodeConditions = new _Formula * [descendants-start];

    for (k=start+1; k<=descendants; k++) {
        node<long>* downWeGo = whereAmI->go_down(k);
        if (!(nodeConditions[k-1-start] = map_node_to_calcnode(downWeGo)->RecurseMC (varToConstrain, downWeGo))) {
            for (long f2 = 0; f2 < k-start-1; f2++) {
                delete nodeConditions[f2];
            }

            delete[] nodeConditions;
            return nil;
        }
    }

    // all the conditions have been written. now check how we should resolve them

    for (k=0; k<descendants-start; k++)
        if ((nodeConditions[k])->GetList().lLength>1) {
            break;
        }

    if (k==descendants-start) { // all underlying branches are "simple"
        for (k=1; k<descendants-start; k++) {
            LocateVar (((_Operation*)((*(nodeConditions[k])).GetList()(0)))->GetAVariable())->SetFormula (*nodeConditions[0]);
            delete (nodeConditions[k]);
            nodeConditions[k] = nil;
        }
        k = 0;
    } else {

        for (l=k+1; l<descendants-start; l++)
            if (nodeConditions[l]->GetList().lLength>1) {
                break;
            }

        if (l==descendants-start) // all but one underlying branches are "simple"
            for (l=0; l<descendants-start; l++) {
                if (l==k) {
                    continue;
                }
                LocateVar (nodeConditions[l]->GetIthTerm(0)->GetAVariable())->SetFormula (*nodeConditions[k]);
                delete (nodeConditions[l]);
                nodeConditions[l] = nil;
            }
        // really bad bongos! must solve for non-additive constraint
        else
            for (l=0; l<descendants-start; l++) {
                if (l==k) {
                    continue;
                }
                if (nodeConditions[l]->GetList().lLength==1) {
                    LocateVar (nodeConditions[l]->GetIthTerm(0)->GetAVariable())->SetFormula (*nodeConditions[k]);
                } else { // solve for a non-additive constraint
                    _Variable* nonAdd = LocateVar (nodeConditions[l]->GetIthTerm(0)->GetAVariable());
                    nodeConditions[l]->GetList().Delete(0);
                    _Formula  newConstraint;
                    newConstraint.Duplicate((BaseRef)nodeConditions[k]);
                    for (long m=0; m<nodeConditions[l]->GetList().lLength; m++) {
                        _Operation* curOp = (_Operation*)(*nodeConditions[l]).GetList()(m);
                        if (curOp->GetNoTerms()) {
                            newConstraint.GetList().AppendNewInstance(new _Operation ('-', 2L));
                        } else {
                            newConstraint.GetList()<<curOp;
                        }
                    }
                    delete (nodeConditions[l]);
                    nodeConditions[l] = nil;
                    nonAdd->SetFormula(newConstraint);
                }
            }
    }


    if(!first) {
        _Formula     *result = nodeConditions[k];

        _Operation   *newVar = new _Operation;
 
        newVar->SetAVariable (iVariables->lData[f-1]);

        result->GetList().AppendNewInstance( newVar);
        result->GetList().AppendNewInstance(new _Operation ('+', 2L));


        delete [] nodeConditions;
        return result;
    }

    for (k=0; k<descendants-start; k++)
        if (nodeConditions[k]) {
            delete nodeConditions[k];
        }

    delete [] nodeConditions;
    return nil;
}

//_______________________________________________________________________________________________

void        _TheTree::ScanSubtreeVars  (_List& rec, char flags, _CalcNode* startAt) const {
    // flags = 1 - do ind
    // flags = 2 - do dep
  
  _SimpleList scanVars;
  _VariableContainer*  thisV;
  _String     chop;
  
  _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER | _HY_TREE_TRAVERSAL_SKIP_ROOT);
  
  if (startAt) {
    thisV = startAt;
  } else {
    thisV = ti.Next();
  }
  
  _AVLList scanVarsA (&scanVars);
  if (flags&0x01) {
    thisV->ScanContainerForVariables (scanVarsA,scanVarsA);
  }
  if (flags&0x02) {
    thisV->ScanForDVariables(scanVarsA,scanVarsA);
  }
  
  scanVarsA.ReorderList();
  
  
  for (unsigned long k=0UL; k<scanVars.lLength; k++) {
    rec.AppendNewInstance(new _String(((_VariableContainer*)LocateVar (scanVars.lData[k]))->ContextFreeName()));
  }
  
  if (startAt) {
    
    _TreeIterator poi (this, _HY_TREE_TRAVERSAL_PREORDER);
    
    while (_CalcNode* iterator = poi.Next()) {
      if (iterator == startAt) {
        break;
      }
    }
    
    long level = poi.Depth();
    
    if (level >= 0) {
      
      while (_CalcNode *iterator = poi.Next()) {
        if (poi.Depth() <= level) {
          break;
        }
        iterator->MatchParametersToList(rec,true,(flags&0x02)!=0);
      }
      
    }
    rec.Clear();
  } else {
    while (_CalcNode * iterator = ti.Next()) {
      iterator->MatchParametersToList(rec,true,(flags&0x02)!=0);
    }
  }
}

  //_______________________________________________________________________________________________

void     _TheTree::AddNodeNamesToDS (_DataSet* ds, bool doTips, bool doInternals, char dOrS) const {
  
  if (dOrS == 2 && doTips && doInternals) {
    AddNodeNamesToDS (ds, false, true, 0);
    AddNodeNamesToDS (ds, true, false, 0);
    return;
  }
  
  const _List nodenames = RetrieveNodeNames (doTips, doInternals,
                                             dOrS ? _HY_TREE_TRAVERSAL_POSTORDER : _HY_TREE_TRAVERSAL_PREORDER);
  
  for (unsigned long i = 0; i < nodenames.countitems(); i++) {
    ds->AddName(* (_String const*)nodenames.GetItem(i));
  }
  
  
}

//_______________________________________________________________________________________________

_String             _TreeTopology::CompareTrees      (_TreeTopology* compareTo) const {
    _List           myLeaves,
                    otherLeaves;

    _SimpleList     indexer,
                    otherIndexer;

    node<long>      *myCT,
         *otherCT;

    _String         rerootAt;

    myCT            = prepTree4Comparison(myLeaves, indexer);
    otherCT         = compareTo->prepTree4Comparison(otherLeaves, otherIndexer);

    // first compare the set of leaf labels

    if (myLeaves.Equal (otherLeaves)) {
        _SimpleList * reindexer = nil;

        if (!indexer.Equal (otherIndexer)) {
            _SimpleList ilist ((unsigned long)myLeaves.lLength);
            ilist.lLength = myLeaves.lLength;

            for (long k2 = 0; k2 < indexer.lLength; k2++) {
                ilist.lData[indexer.lData[k2]] = k2;
            }

            for (long k3 = 0; k3<otherIndexer.lLength; k3++) {
                otherIndexer.lData[k3] = ilist.lData[otherIndexer.lData[k3]];
            }

            reindexer = &otherIndexer;
        }

        char compRes;

        if ((compRes=internalTreeCompare (myCT, otherCT, reindexer, 1, myLeaves.lLength, nil, compareTo))>0) {
            rerootAt = eqWithoutReroot;
        } else {
            long   tCount = 0;

            node_iterator<long> ni (otherCT, _HY_TREE_TRAVERSAL_POSTORDER);
            node<long> * meNode;
          
            while (meNode = ni.Next()) {
              if (meNode == otherCT) {
                break;
              }
              if (meNode->get_num_nodes()) {
                  compRes = internalTreeCompare (myCT, meNode, reindexer, 1, myLeaves.lLength, nil, compareTo);
                  if (compRes>0) {
                      break;
                  } else if (compRes) {
                      meNode = otherCT;
                      break;
                  }
              }

              tCount ++;
            }

            if (meNode!=otherCT) {
                node_iterator<long> ni (compareTo->theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
                while (node<long> * meNode = ni.Next()) {
                    if (meNode == theRoot) {
                      break;
                    }
                    if (tCount==1) {
                        rerootAt = eqWithReroot & *LocateVar (meNode->in_object)->GetName() & '.';
                        break;
                    } else {
                        tCount --;
                    }
                }
            }
        }
        if (!rerootAt.sLength) {
            rerootAt = "Unequal topologies (matching label sets).";
        }
    } else {
        rerootAt = "Unequal label sets.";
    }

    destroyCompTree (myCT);
    destroyCompTree (otherCT);

    return          rerootAt;
}

//_______________________________________________________________________________________________

node<long>* _TreeTopology::prepTree4Comparison (_List& leafNames, _SimpleList& mapping, node<long>* topNode) const {
  
    node<long>* res     = topNode?topNode->duplicate_tree():theRoot->duplicate_tree();
  
    node_iterator<long> ni (res, _HY_TREE_TRAVERSAL_POSTORDER);

    _SimpleList     indexer;

    while (node<long> * iterator = ni.Next()) {
      
        _SimpleList * descendants = new _SimpleList;

        if (!iterator->is_leaf()) {
            for (unsigned long k = 1UL; k <= iterator->get_num_nodes(); k++) {
                (*descendants) << *(_SimpleList*)iterator->go_down(k)->in_object;
            }
        } else {
            (*descendants) << leafNames.lLength;
            indexer        << leafNames.lLength;
            leafNames.AppendNewInstance (new _String (GetNodeName (iterator)));
        }

        iterator->in_object = (long)descendants;
    }

    // now sort leaf names and build the indexer

    mapping.Clear();
    mapping.Duplicate (&indexer);

    SortLists (&leafNames, &indexer);
    SortLists (&indexer, &mapping);

    return res;
}

//_______________________________________________________________________________________________

void    _TreeTopology::destroyCompTree (node<long>* compRoot) const {
    for (unsigned long i=1UL; i<=compRoot->get_num_nodes(); i++) {
        destroyCompTree (compRoot->go_down(i));
    }

    DeleteObject ((BaseRef)compRoot->in_object);
    delete (compRoot);
}

//_______________________________________________________________________________________________


bool    _TheTree::MatchLeavesToDF (_SimpleList& tipMatches, _DataSetFilter* df, bool doNumeric) const {
    // SLKP 20151206 mark for deprecation; only used in HYDataPanel.cpp
    tipMatches.Clear();

    _List           tips = RetrieveNodeNames (true, false, _HY_TREE_TRAVERSAL_PREORDER);

    long number_matched = df->FindSpeciesName (tips, tipMatches);

    if (doNumeric) {
        if (number_matched!=tips.lLength) {
          
            for (unsigned long j=0UL; j<tips.lLength; j++) {
                _String *thisName = (_String*)tips(j);
                long k = atoi (thisName->sData);
                _String tryAgain (k);
                if ( tryAgain.Equal(thisName) && k<=tips.lLength ) {
                    tipMatches<<k;
                } else {
                    return false;
                }
            }
          
            if (tipMatches.Find (0L) < 0) {
              tipMatches.Offset(-1);
            }
        }
    }

    return (number_matched == tips.lLength);
}
//_______________________________________________________________________________________________

_List*  _TreeTopology::SplitTreeIntoClustersInt (node<long>* , _List* , _AVLListX& , long , long ) const
// returns the list of nodes included in the splits BELOW trav1
{
    /*_List * dependancies = nil;

    for (long k = trav1->get_num_nodes(); k; k--)
    {
        _List * me = SplitTreeIntoClustersInt(trav1->go_down(k), res, counts, size, tol);
        if (me)
        {
            if (!dependancies)
                checkPointer (dependancies = new _List);

            (*dependancies) << (*me);
            DeleteObject (me);
        }
    }

    long distinctElements = counts.GetXtra (counts.Find((BaseRef)trav1->in_object)),
         dLength = dependancies?dependancies->lLength:0;

    if (((trav1->parent==nil)||(distinctElements+dLength>=size-tol))&&(distinctElements+dLength<=size+tol))
    {
        _List   entry;
        _String nName;
        GetNodeName (trav1, nName);
        entry && & nName;
        entry && & _Constant (distinctElements);
        if (dependancies)
            entry << (*dependancies);
        (*res) && & entry;

        trav1 = trav1->parent;
        while (trav1)
        {
            long aidx = counts.Find((BaseRef)trav1->in_object);
            counts.SetXtra(aidx,counts.GetXtra (aidx) - distinctElements);
            trav1 = trav1->parent;
        }

        if (dependancies)
            dependancies->Clear();
        else
            checkPointer (dependancies = new _List);

        (*dependancies) && & nName;
    }

    return dependancies;*/

    // decide if child or i-node

    /*if (trav1->get_num_nodes())
    {

    }
    else
    {
        node <long>* pn = trav1->parent;

        while (pn)
        {
            long pidx = counts.Find((BaseRef)pn->in_object);

        }
    }*/
    // TBI
    return new _List;
}

//_______________________________________________________________________________________________

_List*  _TreeTopology::SplitTreeIntoClusters (unsigned long size, unsigned long tol) const {
// assume fixed rooting for now
// returns the list of list speccing pointers to calcnodes rooting the splits
// the member list contains 2 or more entries for each node:
// itself, the number of (independent)leaves encompassed by that node,
// and an optional references to the node whose
// results it depends on.

// assumes that size-tol>=2

    _SimpleList     counts;
    _AVLListX       cavl (&counts);

    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);
    while (node<long> * iterator = ni.Next()) {
        long nC = iterator->get_num_nodes();
        if (nC) {
            long c = 0;
            for (long k=1; k<=nC; k++) {
                c += counts.lData[iterator->go_down(k)->in_object];
            }
            cavl.Insert((BaseRef)iterator->in_object,c);
        } else {
            cavl.Insert((BaseRef)iterator->in_object,1);
        }
    }

    _List*                    result = new _List;

    DeleteObject(SplitTreeIntoClustersInt (theRoot, result, cavl, size, tol));

    return           result;
}

  //_______________________________________________________________________________________________

const _String _TreeTopology::MatchTreePattern (_TreeTopology const* compareTo) const {
    // the pattern is this
    // compare to is the potential match
  _List           myLeaves,
  otherLeaves,
  overlappedLeaves;
  
  _SimpleList     indexer,
  otherIndexer;
  
  node<long>      *myCT,
  *otherCT;
  
  _String         rerootAt;
  
  myCT            = prepTree4Comparison(myLeaves, indexer);
  otherCT         = compareTo->prepTree4Comparison(otherLeaves, otherIndexer);
  
  
  _SimpleList      matchedLeaves;
  
  overlappedLeaves.Intersect (myLeaves, otherLeaves,&matchedLeaves);
  
  if ((myLeaves.lLength>=otherLeaves.lLength)&&(overlappedLeaves.lLength == otherLeaves.lLength)) {
    if (myLeaves.lLength>otherLeaves.lLength) {
        // map leaves back to ordering space
      
      
      
        // trim the pattern tree
      _SimpleList     allowedLeaves  ((unsigned long)myLeaves.lLength),
      recordTransfer ((unsigned long)myLeaves.lLength),
      invertedMap    ((unsigned long)myLeaves.lLength);
      
      allowedLeaves.lLength  = myLeaves.lLength;
      recordTransfer.lLength = myLeaves.lLength;
      invertedMap.lLength    = myLeaves.lLength;
      
      for (long giantsSuck = 0; giantsSuck < myLeaves.lLength; giantsSuck++) {
        invertedMap.lData[indexer.lData[giantsSuck]] = giantsSuck;
      }
      
      for (long padresSuck = 0; padresSuck < matchedLeaves.lLength; padresSuck ++) {
        allowedLeaves.lData[invertedMap.lData[matchedLeaves.lData[padresSuck]]] = 1;
      }
      
      long slider = 0;
      
      for (long dodgersSuck = 0; dodgersSuck < recordTransfer.lLength; dodgersSuck ++)
        if (allowedLeaves.lData[dodgersSuck]) {
          recordTransfer [dodgersSuck] = slider++;
        }
      
        // pass 1 - delete all the superfluous leaves
      
      node_iterator<long> ni (myCT, _HY_TREE_TRAVERSAL_POSTORDER);
      while (node<long>* iterator = ni.Next()) {
        if (iterator != myCT && !iterator->is_leaf()) {
          _SimpleList*  descendants = ((_SimpleList*)iterator->in_object);
          
          if ((descendants&&(allowedLeaves.lData[descendants->lData[0]] == 0))||(!descendants)) {
              // mark this node for deletion
            if (descendants) {
              DeleteObject(descendants);
            }
            
            node<long>* sacLamb = iterator;
            iterator = ni.Next();
            
            if (sacLamb->parent->get_num_nodes()==1) {
              DeleteObject((BaseRef)sacLamb->parent->in_object);
              sacLamb->parent->in_object = 0L;
            }
            
            sacLamb->parent->detach_child (sacLamb->get_child_num());
            delete (sacLamb);
            continue;
          }
        }
      }
      
        // pass 2 - prune internal nodes with exactly one child
      
        // >O
      ni.Reset (myCT);
      
      while (node<long>* iterator = ni.Next()) {
        long nn = iterator->get_num_nodes();
        if (nn == 1) {
          node<long>* loneChild = iterator->go_down(1);
          DeleteObject((_SimpleList*)iterator->in_object);
          iterator->in_object = loneChild->in_object;
          iterator->detach_child (1);
          delete (loneChild);
        } else {
          if (nn > 1) {
            _SimpleList * myDescs = (_SimpleList*)iterator->in_object;
            myDescs->Clear();
            
            for (long cc = 1; cc <= nn; cc++) {
              _SimpleList temp;
              
              temp.Union (*myDescs, *(_SimpleList*)iterator->go_down(cc)->in_object);
              
              myDescs->Clear();
              myDescs->Duplicate (&temp);
            }
            
            /*for (long cc2 = 0; cc2 < myDescs->lLength; cc2++)
             myDescs->lData[cc2] = recordTransfer.lData[myDescs->lData[cc2]];*/
          } else {
            ((_SimpleList*)iterator->in_object)->lData[0] = recordTransfer.lData[((_SimpleList*)iterator->in_object)->lData[0]];
          }
          
        }
      }
      
      _List   newLeaves;
      recordTransfer.Clear();
      invertedMap.Clear();
      
      for (long lc = 0; lc < allowedLeaves.lLength; lc++)
        if (allowedLeaves.lData[lc]) {
          recordTransfer << newLeaves.lLength;
          invertedMap << newLeaves.lLength;
          newLeaves << myLeaves (indexer.lData[lc]);
        }
      
      SortLists (&newLeaves, &recordTransfer);
      SortLists (&recordTransfer,&invertedMap);
      indexer.Clear();
      indexer.Duplicate (&invertedMap);
      myLeaves.Clear();
      myLeaves.Duplicate (&newLeaves);
      
        // finally check whether the root still has 3 or more children
      
      if (myCT->get_num_nodes()==2) {
        node <long>* promoteMe = nil;
        
        if (myCT->go_down(1)->get_num_nodes ()) {
          promoteMe = myCT->go_down(1);
        } else {
          promoteMe = myCT->go_down(2);
        }
        
        long nn = promoteMe->get_num_nodes();
        if (nn) {
          for (long cc = 1; cc <= nn; cc++) {
            myCT->add_node (*promoteMe->go_down(cc));
          }
          
          myCT->detach_child (promoteMe->get_child_num());
          DeleteObject ((BaseRef)promoteMe->in_object);
          delete promoteMe;
        } else {
          WarnError ("Internal tree pattern error in MatchTreePattern");
          return   "Unequal: Error";
        }
      }
    }
    
    _SimpleList * reindexer = nil;
    
    if (!indexer.Equal (otherIndexer)) {
      _SimpleList ilist ((unsigned long)myLeaves.lLength);
      ilist.lLength = myLeaves.lLength;
      
      for (long k2 = 0; k2 < indexer.lLength; k2++) {
        ilist.lData[indexer.lData[k2]] = k2;
      }
      
      for (long k3 = 0; k3<otherIndexer.lLength; k3++) {
        otherIndexer.lData[k3] = ilist.lData[otherIndexer.lData[k3]];
      }
      
      reindexer = &otherIndexer;
    }
    
    char compRes;
    
    if ((compRes=internalTreeCompare (myCT, otherCT, reindexer, 1, myLeaves.lLength, nil, compareTo, true))>0) {
      rerootAt = "Equal w/o rerooting.";
    } else {
      long   tCount = 0;
      
      
      node_iterator <long> ni (otherCT, _HY_TREE_TRAVERSAL_POSTORDER);
      node<long> * iterator = ni.Next();
      
      while (iterator!=otherCT) {
        if (!iterator->is_leaf()) {
          compRes = internalTreeCompare (myCT, iterator, reindexer, 1, myLeaves.lLength, nil, compareTo, true);
          if (compRes>0) {
            break;
          } else if (compRes) {
            iterator = otherCT;
            break;
          }
        }
        iterator = ni.Next();
        
        tCount ++;
      }
      
      if (iterator!=otherCT) {
        ni.Reset (compareTo->theRoot);
        iterator = ni.Next();
        while (iterator!=theRoot) {
          if (tCount==1) {
            rerootAt = _String("Equal with reroot at ") & *LocateVar (iterator->in_object)->GetName() & '.';
            break;
          } else {
            tCount --;
          }
          
          iterator = ni.Next();
        }
      }
    }
    if (!rerootAt.sLength) {
      rerootAt = "Unequal topologies (matching label sets).";
    }
  } else {
    rerootAt = "Unequal label sets.";
  }
  
  destroyCompTree (myCT);
  destroyCompTree (otherCT);
  return          rerootAt;
}





