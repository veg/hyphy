/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2009
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon              (apoon@cfenet.ubc.ca)

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

#include "category.h"
#include "batchlan.h"
#include "likefunc.h"
#include "float.h"

//#ifndef     __ALTIVEC__
//#define     ALMOST_ZERO  1e-305

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

char        isDefiningATree         = false,
            takeBranchLengths       = false,
            autoSolveBranchLengths  = false;

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

//__________________________________________________________________________________

#define MIN_TEX_HEIGHT   50
#define MIN_TEX_WIDTH    50
#define MAX_TEX_WIDTH    160
#define MAX_TEX_HEIGHT   150
#define WIDTH_PER_BRANCH 10


//_______________________________________________________________________________________________

_CalcNode::_CalcNode    ()
{
    theProbs = nil;    // default constructor, doesn't do much
    compExp = nil;
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
void    _CalcNode::InitializeCN     ( _String& parms, int, _VariableContainer* theP, _AVLListXL* aCache)
{
    cBase         = 0;
    theProbs      = nil;
    compExp       = nil;
    referenceNode = -1;
    slaveNodes    = 0;

    long f = parms.Find(','),
         g = -1;

    _String            matrixName (parms,0, f>=0?f-1:-1);

    InitializeVarCont (empty, matrixName, theP, aCache);

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

                for (long k = 0; k<iVariables->lLength; k+=2) {
                    iv  << iVariables->lData[k];
                    iv2 << iVariables->lData[k+1];
                }

                if (dVariables)
                    for (long k = 0; k<dVariables->lLength; k+=2) {
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
    Clear();
    if ((compExp) && (referenceNode < 0)) {
        DeleteObject (compExp);
        compExp = nil;
    }
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

_Parameter  _CalcNode::BranchLength (void)
{

    if (theModel < 0) {
        return Value();
    }

    _FString   *stencil = (_FString*)FetchObjectFromVariableByType (&BRANCH_LENGTH_STENCIL,STRING);

    if (stencil && stencil->theString->Equal (&stringSuppliedLengths)) {
        return Value();
    }


    _Matrix     *freqMx,
                *theMx;

    bool        mbf;

    RetrieveModelComponents (theModel, theMx, freqMx, mbf);

    if (!freqMx && !theModel) {
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

BaseRef _CalcNode::toStr (void)
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

bool        _CalcNode::HasChanged(void)
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

bool        _CalcNode::NeedToExponentiate(long catID)
{
    if (isInOptimize&&(referenceNode>=0)) {
        return ((_CalcNode*)LocateVar(referenceNode))->NeedToExponentiate(catID);
    }

    if (_VariableContainer::NeedToExponentiate(catID>=0)) {
        return true;
    }

    if (catID==-1) {
        if (!compExp) {
            return true;
        }

        for (long i = 0; i<categoryVariables.lLength; i++)
            if (LocateVar (categoryVariables.lData[i])->HasChanged()) {
                return true;
            }
    } else {
        if (!GetCompExp(catID)) {
            return true;
        }

        for (long i = 0; i<categoryVariables.lLength; i++)
            if (((_CategoryVariable*)LocateVar (categoryVariables.lData[i]))->HaveParametersChanged(remapMyCategories.lData[catID*(categoryVariables.lLength+1)+i+1])) {
                return true;
            }
    }
    return false;
}

//_______________________________________________________________________________________________
void        _CalcNode::RecomputeMatrix  (long categID, long totalCategs, _Matrix* storeRateMatrix, _List* queue, _SimpleList* tags)
{
    // assumed that NeedToExponentiate was called prior to this function

    if (isInOptimize)
        if (referenceNode >= 0) {
            _CalcNode* rN = (_CalcNode*)LocateVar(referenceNode);
            rN->RecomputeMatrix (categID, totalCategs, storeRateMatrix);

            if (totalCategs>1) {
                matrixCache[categID] = rN->matrixCache[categID];
                compExp = matrixCache[categID];
            } else {
                compExp = rN->compExp;
            }

            return;
        } else {
            if (referenceNode<-1) {
                slaveNodes++;
                if (slaveNodes>1) {
                    if (slaveNodes == -referenceNode) {
                        slaveNodes = 0;
                    }
                    return;
                }
            }
        }

#ifdef __MP__
    if (matrixTasks) {
        pthread_mutex_lock(&matrixMutex);
    }
#endif
    _Variable* curVar, *locVar;

    long i;

    if (iVariables)
        for (i=0; i<iVariables->lLength; i+=2)
            if (iVariables->lData[i+1]>=0) {
                curVar = LocateVar (iVariables->lData[i+1]);
                if (curVar->IsIndependent()) {
                    locVar = LocateVar (iVariables->lData[i]);
                    curVar->SetValue(locVar->Compute());
                }
            }

    if (dVariables)
        for (i=0; i<dVariables->lLength; i+=2)
            if (dVariables->lData[i+1]>=0) {
                curVar = LocateVar (dVariables->lData[i+1]);
                if (curVar->IsIndependent()) {
                    locVar = LocateVar (dVariables->lData[i]);
                    curVar->SetValue(locVar->Compute());
                }
            }

    for (i=0; i<categoryVariables.lLength; i++) {
        if (categoryIndexVars.lData[i]<0) {
            continue;
        }
        curVar = LocateVar (categoryIndexVars.lData[i]);
        locVar = LocateVar (categoryVariables.lData[i]);
        curVar->SetValue(locVar->Compute());
    }

    if (!storeRateMatrix)
        if (totalCategs>1) {
            DeleteObject(GetCompExp(categID));
        } else if (compExp) {
            DeleteObject (compExp);
        }

    bool    isExplicitForm  = HasExplicitFormModel ();

    _Matrix * myModelMatrix = GetModelMatrix();

    if (myModelMatrix->MatrixType()!=_POLYNOMIAL_TYPE) {
        _Matrix *temp = (_Matrix*)(isExplicitForm?myModelMatrix->makeDynamic():myModelMatrix->MultByFreqs(theModel));

        if (dVariables)
            for (i=0; i<dVariables->lLength; i+=2)
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
            return;
        }

#ifdef __MP__
        if (matrixTasks) {
            temp = (_Matrix*) temp->makeDynamic();
            pthread_mutex_unlock(&matrixMutex);
        }
#endif

#ifdef  _SLKP_LFENGINE_REWRITE_
        if (queue) {
            (*queue) << temp;
            if (tags) {
                (*tags) << (isExplicitForm == 0);
            }
            return;
        }
#endif

        SetCompExp ((_Matrix*)(isExplicitForm?temp:temp->Exponentiate()), totalCategs>1?categID:-1);
#ifdef __MP__
        if (matrixTasks) {
            DeleteObject (temp);
        }
#endif

    } else {
#ifdef __MP__
        if (matrixTasks) {
            pthread_mutex_unlock(&matrixMutex);
        }
#endif
        compExp = (_Matrix*)myModelMatrix->Evaluate(false);
    }
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
        for (long i=0; i<iVariables->lLength; i+=2)
            if (iVariables->lData[i+1]>=0) {
                curVar = LocateVar (iVariables->lData[i+1]);
                if (curVar->IsIndependent()) {
                    locVar = LocateVar (iVariables->lData[i]);
                    curVar->SetValue(locVar->Compute());
                }
            }

    if (dVariables)
        for (long i=0; i<dVariables->lLength; i+=2)
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

_Matrix*    _CalcNode::GetCompExp       (long catID)
{
    if (catID==-1) {
        return compExp;
    } else {
        if (remapMyCategories.lLength) {
            catID = remapMyCategories.lData[catID * (categoryVariables.lLength+1)];
        }

        return matrixCache?matrixCache[catID]:compExp;
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
#if USE_SCALING_TO_FIX_UNDERFLOW
    scalingForUnderflow = nil;
#endif
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

void    _TheTree::PurgeTree (void)
{
    _CalcNode* curNode = DepthWiseTraversal (TRUE), *nextNode;
    nextNode = DepthWiseTraversal ();
    // this loop deletes the data & the structure
    while (nextNode) {
        DeleteVariable (*curNode->GetName());
        curNode = nextNode;
        nextNode = DepthWiseTraversal ();
        delete currentNode;
    }

    DeleteObject (curNode); // error checking is implicit
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

    iNodePrefix = "Node";
    _PMathObj iv = FetchObjectFromVariableByType(&internalNodePrefix,STRING);
    if (iv) {
        iNodePrefix = *((_FString*)iv)->theString;
    }
}

//_______________________________________________________________________________________________

void    _TreeTopology::PreTreeConstructor (bool)
{
    rooted                  = UNROOTED;
    compExp                 = (_Matrix*)checkPointer(new _GrowingVector);

    iNodePrefix = "Node";
    _PMathObj iv = FetchObjectFromVariableByType(&internalNodePrefix,STRING);
    if (iv) {
        iNodePrefix = *((_FString*)iv)->theString;
    }
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
            isDefiningATree = false;
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

    if (dupMe) {
        variablePtrs[theIndex]=this->makeDynamic();
    } else {
        variablePtrs[theIndex]=this;
    }

    DeleteObject(temp);
}

//_______________________________________________________________________________________________
_TheTree::_TheTree              (_String name, _String& parms, bool dupMe):_TreeTopology (&name)
{
#if USE_SCALING_TO_FIX_UNDERFLOW
    scalingForUnderflow = nil;
#endif
    PreTreeConstructor   (dupMe);
    if (MainTreeConstructor  (parms)) {
        PostTreeConstructor  (dupMe);
    }
}

//_______________________________________________________________________________________________
_TheTree::_TheTree              (_String name, _TreeTopology* top):_TreeTopology (&name)
{
#if USE_SCALING_TO_FIX_UNDERFLOW
    scalingForUnderflow = nil;
#endif
    PreTreeConstructor   (false);
    if (top->theRoot) {
        isDefiningATree         = true;
        theRoot                 = top->theRoot->duplicate_tree ();
        node<long>*topTraverser = DepthWiseStepTraverser (theRoot);
        while (topTraverser) {
            _Parameter   nodeVal = top->compExp->theData[topTraverser->in_object];
            _String      nodeVS,
                         nodeName       (*(_String*)top->flatTree(topTraverser->in_object)),
                         nodeParams     (*(_String*)top->flatCLeaves(topTraverser->in_object));

            if (!nodeName.IsValidIdentifier(true)) {
                nodeName.ConvertToAnIdent (true);
            }

            if (nodeVal != 0.0) {
                nodeVS = nodeVal;
            }



            FinalizeNode (topTraverser, 0, nodeName, nodeParams, nodeVS);
            topTraverser = DepthWiseStepTraverser ((node<long>*)nil);
        }
        isDefiningATree         = false;
        PostTreeConstructor      (false);
    } else {
        WarnError ("Can't create an empty tree");
        return;
    }
}

//_______________________________________________________________________________________________
_TreeTopology::_TreeTopology (_TheTree *top):_CalcNode (*top->GetName(), empty)
{
    PreTreeConstructor   (false);
    if (top->theRoot) {
        isDefiningATree         = true;
        theRoot                 = top->theRoot->duplicate_tree ();
        node<long>*topTraverser = DepthWiseStepTraverser (theRoot);
        while (topTraverser && topTraverser->parent) {
            _String              nodeVS,
                                 *nodeSpec,
                                 nodeName;

            top->GetBranchValue (topTraverser,nodeVS);
            top->GetNodeName    (topTraverser,nodeName);
            nodeSpec = top->GetBranchSpec   (topTraverser);

            FinalizeNode (topTraverser, 0, nodeName, *nodeSpec, nodeVS);
            DeleteObject (nodeSpec);
            topTraverser = DepthWiseStepTraverser ((node<long>*)nil);
        }
        isDefiningATree         = false;
    } else {
        WarnError ("Can't create an empty tree");
        return;
    }
}

//_______________________________________________________________________________________________
_TreeTopology::_TreeTopology    (_String name, _String& parms, bool dupMe):_CalcNode (name,empty)
// builds a tree from a string
{
    PreTreeConstructor   (dupMe);
    if (MainTreeConstructor  (parms, false)) {
        PostTreeConstructor  (dupMe);
    } else {
        DeleteObject     (compExp);
        compExp = nil;
    }
}

//_______________________________________________________________________________________________
_TreeTopology::_TreeTopology    (_String* name):_CalcNode (*name,empty)
{

}

//_______________________________________________________________________________________________
bool    _TreeTopology::MainTreeConstructor  (_String& parms, bool checkNames)
{
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
                * newNode   = nil,
                  * parentNode  = nil;

    isDefiningATree         = true;

    for (i=0; i<parms.sLength; i++) {
        switch (parms[i]) {
        case '(': { // creating a new internal node one level down
            // a new node
            newNode = new node<long>;
            checkPointer(newNode);
            if (lastChar == '(' || lastChar == ',') {
                currentNode->add_node (*newNode);
            } else {
                if (theRoot) {
                    parentNode = currentNode->get_parent();
                    if (!parentNode) {
                        WarnError ((_String("'(' is out of context: ...")&parms.Cut(i>31?i-32:0,i)&"?"&parms.Cut(i+1,parms.sLength-i>32?i+32:-1)));
                        isDefiningATree = false;
                        return false;
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
                checkPointer(newNode);
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
                WarnError ((_String(parms[i])&_String(" is out of context:")&parms.Cut(i>31?i-32:0,i)&"?"&parms.Cut(i+1,parms.sLength-i>32?i+32:-1)));
                //PurgeTree();
                isDefiningATree = false;
                return false;
            }
            parentNode = (node<long>*)nodeStack(lastNode);
            FinalizeNode (parentNode, nodeNumbers(lastNode), nodeName, nodeParameters, nodeValue, &nodeComment);
            nodeStack.Delete(lastNode, false);
            nodeNumbers.Delete(lastNode, false);

            if (parms[i]==',') { // also create a new node on the same level
                checkPointer (newNode = new node<long>);
                if (!(parentNode = parentNode->get_parent())) {
                    WarnError ((_String("',' is out of context:")&parms.Cut(i>31?i-32:0,i)&"?"&parms.Cut(i+1,parms.sLength-i>32?i+32:-1)));
                    isDefiningATree = false;
                    return false;
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
                WarnError ((_String("'{' has no matching '}':")&parms.Cut(i>31?i-32:0,i)&"?"&parms.Cut(i+1,parms.sLength-i>32?i+32:-1)));
                isDefiningATree = false;
                return false;
            }
            nodeParameters = parms.Cut (i+1,lastNode-1);
            i = lastNode;
            break;
        }

         case '[' : { // hackish Newick annotation
            lastNode = parms.Find ("]",i+1,-1);
            if (lastNode<0) {
                WarnError ((_String("'[' has no matching ']':")&parms.Cut(i>31?i-32:0,i)&"?"&parms.Cut(i+1,parms.sLength-i>32?i+32:-1)));
                isDefiningATree = false;
                return false;
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
                WarnError ((_String("Node names should begin with a letter, a number, or an underscore. Had:")&parms.Cut(i>31?i-32:0,i)&"?"&parms.Cut(i+1,parms.sLength-i>32?i+32:-1)));
                isDefiningATree = false;
                return false;
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
                        if (c == '\'')
                            if (isInLiteral) {
                                break;
                            } else {
                                WarnError ((_String("Unxpected \'. Had:")&parms.Cut(lastNode>31?lastNode-32:0,lastNode)&"?"&parms.Cut(lastNode+1,parms.sLength-lastNode>32?lastNode+32:-1)));
                                isDefiningATree = false;
                                return false;
                            }
                    } else {
                        break;
                    }

            if (isInLiteral) {
                if (c!='\'') {
                    WarnError ((_String("Unterminated \'. Had:")&parms.Cut(lastNode>31?lastNode-32:0,lastNode)&"?"&parms.Cut(lastNode+1,parms.sLength-lastNode>32?lastNode+32:-1)));
                    isDefiningATree = false;
                    return false;
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

    lastNode = nodeStack.lLength-1;

    while (lastNode>=0) {
        parentNode = (node<long>*)nodeStack(lastNode);
        FinalizeNode (parentNode, nodeNumbers(lastNode), nodeName, nodeParameters, nodeValue, &nodeComment);
        lastNode--;
    }

    if (!theRoot) {
        isDefiningATree = false;
        WarnError ("Can't create empty trees.");
        return false;
    }

    isDefiningATree = false;
    return true;

}
//_______________________________________________________________________________________________


bool    _TheTree::FinalizeNode (node<long>* nodie, long number , _String& nodeName, _String& nodeParameters, _String& nodeValue, _String* nodeComment)
{
    bool isAutoGenerated = (nodeName.sLength == 0);
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
        nodeParameters = empty;
        nodeValue      = empty;
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

    isDefiningATree = isAutoGenerated?2:1;
    _CalcNode cNt (nodeName,nodeParameters, 4, this, aCache);
    isDefiningATree = 1;
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

    nodeVar->SetValue (&val);

    nodeName       = empty;
    nodeParameters = empty;
    nodeValue      = empty;
    if (nodeComment && nodeComment->sLength)
    {
        _String commentName = *nodeVar->GetName() & "._comment";
        CheckReceptacleAndStore(&commentName, empty, false, new _FString (*nodeComment));
        *nodeComment    = empty;
    }

    nodeVar->categoryVariables.TrimMemory();
    nodeVar->categoryIndexVars.TrimMemory();
    nodeVar->_VariableContainer::TrimMemory();


    return true;
}

//_______________________________________________________________________________________________

bool    _TreeTopology::FinalizeNode (node<long>* nodie, long number , _String& nodeName, _String& nodeParameters, _String& nodeValue, _String* nodeComment)
{
    if (!nodeName.sLength) {
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

    nodeName       = empty;
    nodeParameters = empty;
    nodeValue      = empty;
    if (nodeComment)
        *nodeComment    = empty;

    return true;
}

//_______________________________________________________________________________________________

node<long>* _TreeTopology::FindNodeByName (_String* match)
{
    DepthWiseT(true);

    _String       nn;
    while (currentNode) {
        GetNodeName  (currentNode, nn);
        if (match->Equal(&nn)) {
            return currentNode;
        }
        DepthWiseT();
    }

    return nil;

}

//_______________________________________________________________________________________________

void    _TreeTopology::AddANode (_PMathObj newNode)
{
    if (newNode->ObjectClass () == ASSOCIATIVE_LIST) {
        _AssociativeList * newNodeSpec = (_AssociativeList*)newNode;
        _FString         * newName     = (_FString*)newNodeSpec->GetByKey (newNodeGraftName, STRING),
                           * newLocation = (_FString*)newNodeSpec->GetByKey (newNodeGraftWhere, STRING),
                             * newParent   = (_FString*)newNodeSpec->GetByKey (newNodeGraftParent, STRING);


        if (!newName) {
            WarnError (_String("Missing/invalid mandatory argument (\"")&newNodeGraftName&"\") in call to _TreeTopology::AddANode");
            return;
        }
        if (!newLocation) {
            WarnError (_String("Missing/invalid mandatory argument (\"")&newNodeGraftWhere&"\") in call to _TreeTopology::AddANode");
            return;
        }
        if (!newParent) {
            WarnError (_String("Missing/invalid mandatory argument (\"")&newNodeGraftParent&"\") in call to _TreeTopology::AddANode");
            return;
        }

        node<long>* graftAt = FindNodeByName (newLocation->theString);
        if (!graftAt || graftAt->get_parent() == nil) {
            WarnError ("Attachment node must be an exiting non-root node in call to _TreeTopology::AddANode");
            return;
        }

        node<long>* newp = (node<long>*) checkPointer(new node<long>),
                    * curp = graftAt->get_parent();

        newp->set_parent  (*curp);
        newp->add_node    (*graftAt);
        curp->replace_node(graftAt,newp);

        if (!newName->IsEmpty()) {
            node<long>* newt = (node<long>*) checkPointer(new node<long>);
            newp->add_node(*newt);
            FinalizeNode (newt, 0, *newName->theString,   empty, empty);
        }

        FinalizeNode (newp, 0, *newParent->theString, empty, empty);

    } else {
        WarnError ("An invalid argument (not an associative array) supplied to _TreeTopology::AddANode");
    }

}


//_______________________________________________________________________________________________

_CalcNode* _TheTree::DepthWiseTraversal (bool init)
{
    DepthWiseT (init);

    if (currentNode) {
        return (_CalcNode*)(((BaseRef*)variablePtrs.lData)[currentNode->in_object]);
    }

    return nil;
}

//_______________________________________________________________________________________________

void _TreeTopology::DepthWiseT (bool init, _HYTopologyTraversalFunction* handler, Ptr extra)
{
    if (init) {
        currentNode =  DepthWiseStepTraverser (theRoot);
    } else {
        currentNode =(DepthWiseStepTraverser((node<long>*)nil));
    }

    if (handler)
        if (!(*handler)(currentNode, extra)) {
            currentNode = nil;
        }

}

//_______________________________________________________________________________________________

_CalcNode* _TheTree::DepthWiseTraversalRight (bool init)
{
    DepthWiseTRight (init);

    if (currentNode) {
        return (_CalcNode*)(((BaseRef*)variablePtrs.lData)[currentNode->in_object]);
    }

    return nil;
}

//_______________________________________________________________________________________________

void _TreeTopology::DepthWiseTRight (bool init)
{
    if (init) {
        currentNode =  DepthWiseStepTraverserRight (theRoot);
    } else {
        currentNode =(DepthWiseStepTraverserRight((node<long>*)nil));
    }
}

//_______________________________________________________________________________________________

_CalcNode* _TheTree::LeafWiseTraversal (bool init)
{
    LeafWiseT (init);
    if (currentNode) {
        return (_CalcNode*)(((BaseRef*)variablePtrs.lData)[currentNode->in_object]);
    }
    return nil;
}

//_______________________________________________________________________________________________

void _TreeTopology::LeafWiseT (bool init)
{
    if (init) {
        currentNode =  DepthWiseStepTraverser (theRoot);
    } else {
        currentNode = DepthWiseStepTraverser((node<long>*)nil);
    }

    while (currentNode && currentNode->get_num_nodes()) {
        currentNode = DepthWiseStepTraverser((node<long>*)nil);
    }
}

//_______________________________________________________________________________________________

_CalcNode* _TheTree::StepWiseTraversal (bool init)
{
    StepWiseT (init);

    if (currentNode) {
        return (_CalcNode*)(((BaseRef*)variablePtrs.lData)[currentNode->in_object]);
    }

    return nil;
}

//_______________________________________________________________________________________________

void _TreeTopology::StepWiseT (bool init,_HYTopologyTraversalFunction* handler, Ptr extra)
{
    if (init) {
        currentNode =  StepWiseTraverser (theRoot);
    } else {
        currentNode =(StepWiseTraverser((node<long>*)nil));
    }

    if (handler)
        if (!(*handler)(currentNode, extra)) {
            currentNode = nil;
        }
}

//_______________________________________________________________________________________________

void _TreeTopology::StepWiseTLevel (long& level, bool init)
{
    if (init) {
        currentNode =  StepWiseTraverserLevel (level, theRoot);
    } else {
        currentNode =(StepWiseTraverserLevel (level, (node<long>*)nil));
    }

}

//_______________________________________________________________________________________________

_CalcNode* _TheTree::StepWiseTraversalLevel (long& level, bool init)
{
    StepWiseTLevel (level, init);

    if (currentNode) {
        return (_CalcNode*)(((BaseRef*)variablePtrs.lData)[currentNode->in_object]);
    }
    return nil;
}

//_______________________________________________________________________________________________

void _TreeTopology::DepthWiseTLevel (long& level, bool init)
{
    if (init) {
        currentNode =  DepthWiseStepTraverserLevel (level, theRoot);
    } else {
        currentNode =(DepthWiseStepTraverserLevel (level, (node<long>*)nil));
    }
}
//_______________________________________________________________________________________________

_CalcNode* _TheTree::DepthWiseTraversalLevel (long& level, bool init)
{
    DepthWiseTLevel (level, init);

    if (currentNode) {
        return (_CalcNode*)(((BaseRef*)variablePtrs.lData)[currentNode->in_object]);
    }
    return nil;
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

    res->currentNode = currentNode;
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

    res->currentNode = currentNode;
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

_PMathObj       _TreeTopology::Compute  (void)
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

bool     _TreeTopology::IsCurrentNodeATip (void)
{
    return currentNode?currentNode->get_num_nodes()==0:false;
}

//_______________________________________________________________________________________________

bool     _TreeTopology::IsCurrentNodeTheRoot (void)
{
    return currentNode->parent==nil;
}

//_______________________________________________________________________________________________

bool     _TreeTopology::IsDegenerate (void)
{
    return theRoot->get_num_nodes()==1;
}

//_______________________________________________________________________________________________

void    _TreeTopology::GetNodeName (node<long>* n, _String& r, bool fullName)
{
    if (fullName) {
        r = *GetName()&'.'&*(_String*)(flatTree (n->in_object));
    } else {
        r = *(_String*)(flatTree (n->in_object));
    }
}

//_______________________________________________________________________________________________

void    _TreeTopology::SetLeafName(long res, _String* newName)
{
    long count      = 0;

    LeafWiseT(true);
    while (currentNode) {
        if (res==count) {
            flatTree.Replace (currentNode->in_object, newName, true);
            break;
        }
        count++;
        LeafWiseT(false);
    }
}
//_______________________________________________________________________________________________

void    _TheTree::GetNodeName      (node<long>* n, _String& r, bool fullName)
{
    if (fullName) {
        r = *((_CalcNode*)(((BaseRef*)variablePtrs.lData)[n->in_object]))->GetName();
    } else {
        r = ((_CalcNode*)(((BaseRef*)variablePtrs.lData)[n->in_object]))->GetName()->Cut (GetName()->sLength+1,-1);
    }
}


//_______________________________________________________________________________________________

BaseRef     _TheTree::toStr (void)
{
    _String     * res = new _String((unsigned long)128,true),
    num;

    _Parameter    skipILabels,
                  includeMSP;

    checkParameter (noInternalLabels, skipILabels, 0.0);
    checkParameter (includeModelSpecs, includeMSP , 0.0);

    if (IsDegenerate ()) {
        _CalcNode*  curNode  = DepthWiseTraversal(true),
                    *  nextNode = DepthWiseTraversal(false);

        long     l1             = GetName()->Length();

        (*res)<<'(';
        num = nextNode->GetName()->Cut(l1+1,-1);
        (*res)<<&num;
        if (includeMSP>0.5) {
            long midx = curNode->GetModelIndex();
            if (midx != HY_NO_MODEL) {
                (*res) << '{';
                (*res) << (_String*)modelNames (midx);
                (*res) << '}';
            }
        }
        (*res)<<',';
        num = curNode->GetName()->Cut(l1+1,-1);
        (*res)<<&num;
        if (includeMSP>0.5) {
            long midx = nextNode->GetModelIndex();
            if (midx != HY_NO_MODEL) {
                (*res) << '{';
                (*res) << (_String*)modelNames (midx);
                (*res) << '}';
            }
        }
        (*res)<<')';
    } else {

        long        level       =   0,
                    myLevel     =   0,
                    lastLevel    = 0,
                    l1          = GetName()->Length(),
                    j;


        _CalcNode*  curNode=DepthWiseTraversalLevel(myLevel,true),
                    *  nextNode;

        level       = myLevel;

        bool        isCTip = IsCurrentNodeATip(),
                    isCTip2;

        nextNode    = DepthWiseTraversalLevel(myLevel) ;

        isCTip2     = IsCurrentNodeATip();

        while (nextNode) {
            if (level>lastLevel) {
                if (lastLevel) {
                    (*res)<<',';
                }
                for (j=0; j<level-lastLevel; j++) {
                    (*res)<<'(';
                }
            } else if (level<lastLevel) {
                for (j=0; j<lastLevel-level; j++) {
                    (*res)<<')';
                }
            } else {
                (*res)<<',';
            }

            if ((skipILabels<0.1)||isCTip) {
                num = curNode->GetName()->Cut(l1+1,-1);
                (*res)<<&num;
            }
            if (includeMSP>0.5) {
                long midx = curNode->GetModelIndex();
                if (midx != HY_NO_MODEL) {
                    (*res) << '{';
                    (*res) << (_String*)modelNames (midx);
                    (*res) << '}';
                }
            }

            lastLevel = level;
            level     = myLevel;
            curNode   = nextNode;
            isCTip    = isCTip2;
            nextNode  = DepthWiseTraversalLevel(myLevel) ;
            isCTip2   = IsCurrentNodeATip();

        }
        for (j=0; j<lastLevel-level; j++) {
            (*res)<<')';
        }
    }
    (*res)<<';';
    (*res).Finalize();
    return res;
}


//_______________________________________________________________________________________________

BaseRef     _TreeTopology::toStr (void)
{
    _String     * res = new _String((unsigned long)128,true),
    num;

    _Parameter    skipILabels,
                  includeMSP;

    checkParameter (noInternalLabels, skipILabels, 0.0);
    checkParameter (includeModelSpecs, includeMSP , 0.0);

    if (IsDegenerate ()) {

        DepthWiseT (true);

        (*res)<<'(';

        GetNodeName (theRoot, num);

        (*res)<<&num;
        if (includeMSP>0.5) {
            _String *mSpec = (_String*)flatCLeaves(theRoot->in_object);
            if (mSpec->sLength) {
                (*res) << '{';
                (*res) << mSpec;
                (*res) << '}';
            }
        }
        (*res)<<',';
        GetNodeName (currentNode, num);
        (*res)<<&num;
        if (includeMSP>0.5) {
            _String *mSpec = (_String*)flatCLeaves(currentNode->in_object);
            if (mSpec->sLength) {
                (*res) << '{';
                (*res) << mSpec;
                (*res) << '}';
            }
        }
        (*res)<<')';
    } else {

        long        level       =   0,
                    myLevel     =   0,
                    lastLevel    = 0,
                    j;

        DepthWiseTLevel(myLevel,true);

        node<long>*     curNode= currentNode,
                        *     nextNode;

        level       = myLevel;

        bool        isCTip = IsCurrentNodeATip(),
                    isCTip2;

        DepthWiseTLevel(myLevel) ;
        nextNode = currentNode;

        isCTip2     = IsCurrentNodeATip();

        while (nextNode) {
            if (level>lastLevel) {
                if (lastLevel) {
                    (*res)<<',';
                }
                for (j=0; j<level-lastLevel; j++) {
                    (*res)<<'(';
                }
            } else if (level<lastLevel) {
                for (j=0; j<lastLevel-level; j++) {
                    (*res)<<')';
                }
            } else {
                (*res)<<',';
            }

            if ((skipILabels<0.1)||isCTip) {
                GetNodeName (curNode, num);
                (*res)<<&num;
            }

            if (includeMSP>0.5) {
                _String *mSpec = (_String*)flatCLeaves(curNode->in_object);
                if (mSpec->sLength) {
                    (*res) << '{';
                    (*res) << mSpec;
                    (*res) << '}';
                }
            }

            lastLevel = level;
            level     = myLevel;
            curNode   = nextNode;
            isCTip    = isCTip2;
            DepthWiseTLevel(myLevel) ;
            nextNode = currentNode;
            isCTip2   = IsCurrentNodeATip();

        }
        for (j=0; j<lastLevel-level; j++) {
            (*res)<<')';
        }
    }
    (*res)<<';';
    (*res).Finalize();
    return res;
}

//__________________________________________________________________________________
void _TreeTopology::toFileStr(FILE* f)
{
    _String * s = (_String*)toStr();
    fprintf (f, "%s", s->sData);
    DeleteObject(s);
}

//__________________________________________________________________________________

void _TheTree::CompileListOfModels (_SimpleList& l)
{
    _CalcNode* curNode = DepthWiseTraversal (true);
    while (curNode) {
        long    modelID = curNode->GetModelIndex();
        if (modelID != HY_NO_MODEL && l.Find(modelID) == -1 ) {
            l << modelID;
        }
        curNode = DepthWiseTraversal (false);
    }
}
//_______________________________________________________________________________________________
void    _TheTree::SetCompMatrices (long catID)
{
    _CalcNode* travNode = DepthWiseTraversal (TRUE);

    while   (!IsCurrentNodeTheRoot()) {
        travNode->SetCompMatrix (catID);
        travNode = DepthWiseTraversal ();
    }
}

//__________________________________________________________________________________

void _TheTree::SetUp (void)
{
    _CalcNode* travNode = DepthWiseTraversal (TRUE);

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


    while   (travNode) {
        if (!IsCurrentNodeATip()) {
            flatTree<<travNode;
            flatNodes<<(long)(currentNode);
            travNode->lastState = -1;
            flatINodeParents << (long)(currentNode->parent);
        } else {
            flatLeaves << (long)(currentNode);
            flatCLeaves << travNode;
            flatParents << (long)(currentNode->parent);
        }
        travNode = DepthWiseTraversal ();
    }

    flatParents << flatINodeParents;
    _SimpleList parentlist (flatNodes), indexer (flatNodes.lLength,0,1);
    SortLists   (&parentlist,&indexer);
    for (long k=0; k<flatParents.lLength; k++)
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

    long    iNodeCounter = 0,
            leafCounter = 0;

    travNode = DepthWiseTraversal (TRUE);

    while   (travNode) {
        if (IsCurrentNodeATip()) {
            travNode->nodeIndex = leafCounter++;
        } else {
            nodeMarkers[iNodeCounter] = -1;
            for (long k=1; k<systemCPUCount; k++) {
                nodeMarkers[iNodeCounter+k*flatNodes.lLength] = -1;
            }
            travNode->nodeIndex = flatLeaves.lLength+iNodeCounter++;
            nodeStates[travNode->nodeIndex]=-1;
            for (long m=1; m<systemCPUCount; m++) {
                nodeStates[travNode->nodeIndex+m*(flatNodes.lLength+flatLeaves.lLength)] = -1;
            }
        }
        travNode = DepthWiseTraversal ();
    }

    BuildINodeDependancies();
}

//__________________________________________________________________________________

bool _TheTree::AllBranchesHaveModels (long matchSize)
{
    _CalcNode* travNode;
    travNode = DepthWiseTraversal (TRUE);
    if (matchSize>0) {
        while   (!IsCurrentNodeTheRoot()) {
            long  mID = travNode->GetTheModelID();
            if (mID<0) {
                return false;
            }
            travNode = DepthWiseTraversal ();
        }
    } else {
        while   (!IsCurrentNodeTheRoot()) {
            long  mID = travNode->GetTheModelID();
            if (mID<0) {
                return false;
            } else {
                if (travNode->GetModelMatrix()->GetHDim()!=matchSize) {
                    return false;
                }
            }
            travNode = DepthWiseTraversal ();
        }

    }
    return true;
}

//__________________________________________________________________________________

_String*    _TheTree::TreeUserParams (void)
{
    _String * result = new _String (16L, true);
    checkPointer (result);

    _CalcNode* travNode;
    travNode = DepthWiseTraversal (TRUE);
    while   (travNode) {
        _String * nodeString = travNode->GetSaveableListOfUserParameters();
        if (nodeString->sLength) {
            *result << nodeString;
        }
        DeleteObject (nodeString);
        travNode = DepthWiseTraversal ();
    }

    result->Finalize();
    return result;
}

//__________________________________________________________________________________

_PMathObj _TreeTopology::Execute (long opCode, _PMathObj p, _PMathObj p2)   // execute this operation with the second arg if necessary
{

    switch (opCode) {
    case HY_OP_CODE_IDIV: { // Split ($) - 2nd argument
        if (p->ObjectClass()!=NUMBER) {
            _String errMsg ("Invalid (not a number) 2nd argument is call to $ for trees.");
            WarnError (errMsg);
        }
        _Constant*  cc     = (_Constant*)TipCount();
        long        size   = cc->Value()/p->Value();

        if  ((size<=4)||(size>cc->Value()/2)) {
            _String errMsg ("Poor choice of the 2nd numeric agrument in to $ for tree. Either the resulting cluster size is too big(>half of the tree), or too small (<4)!");
            WarnError (errMsg);
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

                    CheckReceptacle (&splitNodeNames, empty, false)->SetValue (&selMatrix);
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
    break;

    case HY_OP_CODE_MUL: // compute the strict consensus between T1 and T2
        return SplitsIdentity (p);
        break;

    case HY_OP_CODE_ADD: // +
        if (!p) {
            return Sum();
        }
        AddANode (p);
        return new _Constant (0.0);
        break;

    case HY_OP_CODE_LEQ: { // MatchPattern (<=)
        if ((p->ObjectClass()!=TREE)&&(p->ObjectClass()!=TOPOLOGY)) {
            _String errMsg ("Invalid (not a tree/topology) 2nd argument is call to <= for trees/topologies.");
            WarnError (errMsg);
        }
        _String  res (((_TreeTopology*)p)->MatchTreePattern (this));
        return new _Constant (!res.beginswith ("Unequal"));
        break;
    }
    case HY_OP_CODE_EQ: // ==
        return  Compare(p);
        break;
    case HY_OP_CODE_ABS: // Abs
        return FlatRepresentation();
        break;
    case HY_OP_CODE_BRANCHCOUNT: //BranchCount
        return BranchCount();
        break;
    case HY_OP_CODE_BRANCHLENGTH: //BranchLength
        return BranchLength(p);
        break;
    case HY_OP_CODE_BRANCHNAME: //BranchName
        return BranchName(p);
        break;
    case HY_OP_CODE_FORMAT: { // Format
        currentNode = theRoot;
        _String  *tStr = new _String  ((unsigned long)1024,true);
        SubTreeString (*tStr, p->Compute()->Value() > 0.1 , p2->Compute()->Value() > 0.1 ? -3:-1, nil);
        tStr->Finalize();
        return new _FString (tStr);
    }
    case HY_OP_CODE_MACCESS: // MAccess
        return BranchName (p,true, p2);
        break;
    case HY_OP_CODE_MIN: // COT (Min)
        return FindCOT (p);
        break;
    case HY_OP_CODE_REROOTTREE: // RerootTree
        return RerootTree(p);
        break;
    case HY_OP_CODE_TEXTREESTRING: // TEXTreeString
        //return TEXTreeString(p);
        break;
    case HY_OP_CODE_TIPCOUNT: // TipCount
        return TipCount();
        break;
    case HY_OP_CODE_TIPNAME: // TipName
        return TipName(p);
        break;
    case HY_OP_CODE_TYPE: // Type
        return Type();
        break;
    case HY_OP_CODE_POWER: //^
        return AVLRepresentation (p);
        break;
    }

    WarnNotDefined (this, opCode);
    return nil;

}


//__________________________________________________________________________________

_PMathObj _TheTree::Execute (long opCode, _PMathObj p, _PMathObj p2)   // execute this operation with the second arg if necessary
{

    switch (opCode) {
    case HY_OP_CODE_PSTREESTRING: //PlainTreeString
        return PlainTreeString(p,p2);
        break;
    case HY_OP_CODE_TEXTREESTRING: // TEXTreeString
        return TEXTreeString(p);
        break;
    case HY_OP_CODE_TYPE: // Type
        return Type();
        break;
    }

    return  _TreeTopology::Execute (opCode,p,p2);

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

_AssociativeList* _TreeTopology::FindCOT (_PMathObj p)
// Find the Center of the Tree (COT) location
// using an L_p metric (L_2 works well)
{
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

    DepthWiseT        (true);

    while (currentNode->parent) {
        if (IsCurrentNodeATip ()) {
            addressToIndexMap.Insert ((BaseRef)currentNode, leafCount++);
        } else {
            branchCount++;
        }

        addressToIndexMap2.Insert ((BaseRef)currentNode, tIndex++);
        DepthWiseT (false);
    }

    // allocate the matrix of path lengths with hardwired (traversal order) indices
    // also allocate a list of sorted lists to store children nodes
    // and a map of (longed) node addresses to post order traversal indices

    _Matrix      distances          (branchCount+leafCount, leafCount, false, true),
                 rootDistances        (1, leafCount,  false, true),
                 branchLengths        (1, branchCount+leafCount,  false, true),
                 branchSpans      (branchCount+leafCount+1,2,false,true);

    _String      nodeName;

    // pass 1: fill up the nodes up to the root (i.e. below any internal node)

    DepthWiseT        (true);
    tIndex            = 0;

    while (currentNode->parent) {

        _Parameter          myLength;
        GetBranchLength     (currentNode, myLength);
        lengthToIndexMap.Insert (new _String(totalTreeLength), tIndex, false, true);
        totalTreeLength      += myLength;

        branchLengths.Store (0, tIndex++, myLength);
        listOfNodes << (long)currentNode;
        _SimpleList         *childIndices = (_SimpleList*)checkPointer(new _SimpleList);

        if (IsCurrentNodeATip ()) {
            (*childIndices) << addressToIndexMap.GetXtra(addressToIndexMap.Find((BaseRef)currentNode));
            //leafNames<<LocateVar(currentNode->in_object)->GetName();
        } else {
            long           myIndex = addressToIndexMap2.GetXtra(addressToIndexMap2.Find((BaseRef)currentNode));

            _SimpleList    mappedLeaves (leafCount,0,0);
            for (long ci = currentNode->get_num_nodes(); ci; ci--) {
                long          childIndex = addressToIndexMap2.GetXtra(addressToIndexMap2.Find((BaseRef)currentNode->go_down (ci)));
                _SimpleList * childLeaves = (_SimpleList*)childLists(childIndex);

                myLength = branchLengths.theData[childIndex];

                for (long ci2 = 0; ci2 < childLeaves->lLength; ci2++) {
                    long ttIndex = childLeaves->lData[ci2];
                    mappedLeaves.lData[ttIndex] = 1;
                    distances.Store (myIndex, ttIndex, distances (childIndex, ttIndex) + myLength);
                }
            }

            for (long ci2 = 0; ci2 < leafCount; ci2++)
                if (mappedLeaves.lData[ci2]) {
                    (*childIndices) << ci2;
                }
        }

        childLists.AppendNewInstance(childIndices);
        DepthWiseT (false);
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
    GetNodeName     (cotBranch,nodeName);


    resList->MStore (cotNode, new _FString(nodeName,false),false);
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
        sprintf  (buffer, "%.15f", branchSpans(sc,0));
        nodeName = buffer;
        branchSpans.Store(sc,0,nodeName.toNum());
        timeSplitsAVL.Insert (nodeName.makeDynamic(),0,false,true);
        sprintf  (buffer, "%.15f", branchSpans(sc,1));
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
            GetNodeName     ((node<long>*)listOfNodes(mxc), nodeName);
            ctl->MStore (nodeName, new _Constant (T1), false);
        }

        char       buffer[256];
        sprintf  (buffer, "%.15f", T0);
        nodeName = buffer;
        tcache.Clear();
        long       startingPos =  timeSplitsAVL.Find (&nodeName,tcache),
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
                nodeName = tSample;
                long branchIndex = 0;
                if (lengthToIndexMap.FindBest (&nodeName,branchIndex)<=0) {
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

char     _TreeTopology::internalNodeCompare (node<long>* n1, node<long>* n2, _SimpleList& subTreeMap, _SimpleList* reindexer, bool cangoup, long totalSize, node<long>* n22, _TreeTopology* tree2, bool isPattern)
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
                if (complement)
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

            if (k6 == stSizes.lLength)
                if (isPattern) {
                    unmatchedPatterns << k5;
                    subTreeMap << -2;
                } else {
                    return 0;
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

        /*_String nm1, nm2;
        tree2->GetNodeName (n1,nm1);
        GetNodeName (n2,nm2);
        printf ("Node match %s, %s\n",  nm1.sData, nm2.sData, "\n");*/

        return 1;
    }

    return 0;
}


//long   itcCount = 0;
//__________________________________________________________________________________

char     _TreeTopology::internalTreeCompare (node<long>* n1, node<long>* n2, _SimpleList* reindexer, char compMode, long totalSize, node<long>* n22, _TreeTopology* tree2, bool isPattern)
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

                checkPointer (newLeaves);
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

_String  _TheTree::FindMaxCommonSubTree (_TheTree* compareTo, long& sizeVar, _List* forest)
{
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

        node<long>* meNode = DepthWiseStepTraverser (myCT);

        while (meNode) {
            if (meNode->get_num_nodes() == 0) {
                ldx1 << (long)meNode;
            }
            meNode = DepthWiseStepTraverser ((node<long>*)nil);
        }

        meNode = DepthWiseStepTraverser (otherCT);

        while (meNode) {
            if (meNode->get_num_nodes() == 0) {
                ldx2 << (long)meNode;
            }
            meNode = DepthWiseStepTraverser ((node<long>*)nil);
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

                    node<long>* meNode = DepthWiseStepTraverser (myCT),
                                *   mNode  = (node<long>*)matchedTops.lData[k6];

                    while (meNode) {
                        if (meNode->get_num_nodes()) {
                            if (meNode == mNode) {
                                break;
                            }
                            maxSz ++;
                        }
                        meNode = DepthWiseStepTraverser ((node<long>*)nil);
                    }

                    meNode = DepthWiseStepTraverser (theRoot);

                    while (meNode) {
                        if (meNode->get_num_nodes()) {
                            if (maxSz == 0) {
                                (*forest) << LocateVar(meNode->in_object)->GetName();
                                break;
                            }
                            maxSz--;
                        }
                        meNode = DepthWiseStepTraverser ((node<long>*)nil);
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

                node<long>* meNode = DepthWiseStepTraverser (myCT),
                            *   mNode  = (node<long>*)matchedTops.lData[maxIdx];

                while (meNode) {
                    if (meNode->get_num_nodes()) {
                        if (meNode == mNode) {
                            break;
                        }
                        maxSz ++;
                    }
                    meNode = DepthWiseStepTraverser ((node<long>*)nil);
                }

                meNode = DepthWiseStepTraverser (theRoot);

                while (meNode) {
                    if (meNode->get_num_nodes()) {
                        if (maxSz == 0) {
                            return *LocateVar(meNode->in_object)->GetName();
                        }
                        maxSz--;
                    }
                    meNode = DepthWiseStepTraverser ((node<long>*)nil);
                }
            }
        }
    }
    return empty;
}


//__________________________________________________________________________________

_String  _TheTree::CompareSubTrees (_TheTree* compareTo, node<long>* topNode)
// compare tree topologies
{
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

        {
            //for BCC
            for (long k2 = 0; k2 < ilist.lLength; k2++) {
                ilist.lData[k2] = -1;
            }
        }
        {
            //for BCC
            for (long k2 = 0; k2 < otherIndexer.lLength; k2++) {
                ilist.lData[sharedLeavesID.lData[otherIndexer.lData[k2]]] = k2;
            }
        }
        for (long k2 = 0; k2<indexer.lLength; k2++) {
            long         lidx      = ilist.lData[indexer.lData[k2]];

            if (lidx>=0) {
                indexer.lData[k2] = lidx;
            } else {
                indexer.lData[k2] = -1;
            }
        }

        _SimpleList *reindexer = &indexer;

        // now compare explore possible subtree matchings
        // for all internal nodes of this tree except the root

        char   compRes = 0;

        long   tCount = 1,
               nc2 = topNode->get_num_nodes();

        node<long>* meNode = DepthWiseStepTraverser (myCT);
        meNode = DepthWiseStepTraverser ((node<long>*)nil);

        while (meNode!=myCT) {
            long nc = meNode->get_num_nodes();
            if (nc == nc2) {
                long kk;
                for (kk = 1; kk <= nc; kk++) {
                    compRes = internalTreeCompare (otherCT,meNode, reindexer, 0, otherLeaves.lLength, meNode->go_down(kk),compareTo);
                    if (compRes) {
                        if (compRes == -1) {
                            meNode = myCT;
                        }
                        break;
                    }
                }
                if (kk>nc) {
                    compRes = internalTreeCompare (otherCT,meNode, reindexer, 0, otherLeaves.lLength, nil, compareTo);
                    if (compRes) {
                        if (compRes == -1) {
                            meNode = myCT;
                        }
                        break;
                    }
                } else {
                    break;
                }
            }
            tCount ++;
            meNode = DepthWiseStepTraverser ((node<long>*)nil);
        }

        if (meNode!=myCT) {
            meNode = DepthWiseStepTraverser (theRoot);
            while (meNode!=theRoot) {
                if (tCount==0) {
                    rerootAt = _String("Matched at the ") & *LocateVar (meNode->in_object)->GetName() & '.';
                    break;
                } else {
                    tCount --;
                }

                meNode = DepthWiseStepTraverser ((node<long>*)nil);
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
void _TreeTopology::EdgeCount (long& leaves, long& internals)
{
    leaves    = 0;
    internals = 0;
    DepthWiseT(true);
    while (currentNode) {
        if (IsCurrentNodeATip()) {
            leaves ++;
        } else {
            internals ++;
        }

        DepthWiseT(false);
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

    node      <long>* tNode = DepthWiseStepTraverser (theRoot);

    long      count = 0;

    while     (tNode) {
        flatTree << tNode->in_object;
        tNode->in_object = count++;
        tNode = DepthWiseStepTraverser ((node<long>*)nil);
    }


    _Matrix * res = new _Matrix (1,count, false, true);
    checkPointer (res);

    tNode = DepthWiseStepTraverser (theRoot);
    count = 0;

    while     (tNode) {
        if (tNode->parent) {
            res->theData[count] = tNode->parent->in_object;
        } else {
            res->theData[count] = -1;
        }

        tNode->in_object = flatTree.lData[count++];
        tNode = DepthWiseStepTraverser ((node<long>*)nil);
    }

    return res;
}

//__________________________________________________________________________________
_PMathObj _TreeTopology::AVLRepresentation (_PMathObj layoutOption)
{

    if (layoutOption->ObjectClass () == NUMBER) {
        bool               preOrder = layoutOption->Compute()->Value()>0.5;

        _AssociativeList * masterList = (_AssociativeList * ) checkPointer(new _AssociativeList ());
        _FString           nameHolder;
        //             arrayKey;

        _Constant          lengthHolder;

        long         rootIndex = 0,
                     nodeLevel = 0;

        _SimpleList  nodeList;
        _AVLListX    nodeIndexList (&nodeList);

        node      <long>* tNode = preOrder?StepWiseTraverser(theRoot):DepthWiseStepTraverser (theRoot);

        while     (tNode) {
            nodeIndexList.Insert ((BaseObj*)tNode, nodeIndexList.countitems()+1);

            if (tNode->parent == nil) {
                rootIndex = nodeIndexList.countitems();
            }

            tNode = preOrder?StepWiseTraverser((node<long>*)nil):DepthWiseStepTraverser ((node<long>*)nil);
        }

        tNode = preOrder?StepWiseTraverserLevel(nodeLevel,theRoot):DepthWiseStepTraverserLevel (nodeLevel,theRoot);

        while     (tNode) {
            _AssociativeList * nodeList = (_AssociativeList * ) checkPointer(new _AssociativeList ());
            GetNodeName (tNode, *nameHolder.theString);
            nodeList->MStore ("Name", &nameHolder, true);
            GetBranchLength (tNode, lengthHolder.theValue);
            nodeList->MStore ("Length", &lengthHolder, true);
            lengthHolder.theValue = nodeLevel;
            nodeList->MStore ("Depth", new _Constant (nodeLevel), false);
            if (tNode->parent) {
                nodeList->MStore ("Parent", new _Constant(nodeIndexList.GetXtra (nodeIndexList.Find((BaseObj*)tNode->parent))), false);
            }

            long nCount = tNode->get_num_nodes();
            if (nCount) {
                _AssociativeList * childList = new _AssociativeList ();
                checkPointer (childList);
                for (long k = 1; k<=nCount; k=k+1) {
                    childList->MStore (_String((long)(k-1)),new _Constant(nodeIndexList.GetXtra (nodeIndexList.Find((BaseObj*)tNode->go_down(k)))) , false);
                }
                nodeList->MStore ("Children", childList, false);
            }
            masterList->MStore (_String((long)nodeIndexList.GetXtra (nodeIndexList.Find((BaseObj*)tNode))), nodeList, false);
            tNode = preOrder?StepWiseTraverserLevel(nodeLevel,(node<long>*)nil):DepthWiseStepTraverserLevel (nodeLevel,(node<long>*)nil);
        }

        _AssociativeList * headerList = new _AssociativeList ();
        checkPointer (headerList);

        headerList->MStore ("Name", new _FString (*GetName()), false);
        headerList->MStore ("Root", new _Constant(rootIndex), false);
        masterList->MStore ("0", headerList, false);

        return masterList;
    }
    return new _Constant (0.0);
}

//__________________________________________________________________________________
_PMathObj _TreeTopology::TipName (_PMathObj p)
{
    _String resString;

    if (p&& p->ObjectClass()==NUMBER) {
        long res        = p->Value(),
             count         = 0;

        _List           * allLeaves = nil;

        if (res < 0) {
            allLeaves = (_List*)checkPointer(new _List);
        }

        LeafWiseT(true);

        while (currentNode) {
            if (res < 0) {
                GetNodeName (currentNode,resString);
                (*allLeaves) && & resString;
            } else if (res==count) {
                //resString = travNode->GetName()->Cut(travNode->GetName()->Find ('.')+1,-1);
                GetNodeName (currentNode,resString);
                break;
            }
            count++;
            LeafWiseT(false);
        }

        if (res < 0) {
            _Matrix * res = new _Matrix (*allLeaves);
            DeleteObject (allLeaves);
            return res;
        }
    }

    return new _FString (resString,false);
}

//__________________________________________________________________________________
_PMathObj _TreeTopology::BranchLength (_PMathObj p)
{
    _Parameter resValue = HY_INVALID_RETURN_VALUE;


    if (p) {
        if (p->ObjectClass()==NUMBER) {
            long res        = p->Value(),
                 count         = 0;

            if (res < 0)
                // get ALL branch lengths
            {
                EdgeCount (count, res);
                _Matrix*   branchLengths = (_Matrix*) checkPointer(new _Matrix (1,count+res,false,true));

                count       = 0;
                DepthWiseT  (true);
                while (!IsCurrentNodeTheRoot()) {
                    GetBranchLength (currentNode, branchLengths->theData[count++]);
                    DepthWiseT      (false);
                }
                return branchLengths;
            } else
                // get a branch length
            {
                DepthWiseT(true);
                while (currentNode && res != count) {
                    count++;
                    DepthWiseT(false);
                }
                if (currentNode && !IsCurrentNodeTheRoot()) {
                    GetBranchLength (currentNode,resValue);
                }
            }
        } else {
            if (p->ObjectClass()==STRING) {
                _List * twoIDs = ((_FString*)p->Compute())->theString->Tokenize(";");
                if (twoIDs->lLength == 2 || twoIDs->lLength == 1) {
                    _String * node1 = (_String*)(*twoIDs) (0),
                              * node2 = twoIDs->lLength>1?(_String*)(*twoIDs) (1):nil;

                    node<long>* n1 = nil,
                                * n2 = nil;

                    long        l1 = 0,
                                l2 = 0,
                                l  = 0;

                    DepthWiseTLevel(l,true);

                    _String      cBranchName;

                    while (currentNode && (!n1 || !n2)) {
                        GetNodeName (currentNode, cBranchName);
                        if (cBranchName.Equal (node1)) {
                            n1 = currentNode;
                            l1 = l;
                        } else if (node2 && cBranchName.Equal (node2)) {
                            n2 = currentNode;
                            l2 = l;
                        }

                        DepthWiseTLevel  (l,false);
                    }

                    if (n1 && n2) {
                        _Parameter p1 = 0,
                                   p2 = 0,
                                   p;

                        while (l1<l2) {
                            GetBranchLength (n2,p);
                            p2 += p;
                            n2 = n2->parent;
                            l2--;
                        }

                        while (l2<l1) {
                            GetBranchLength (n1,p);
                            p1 += p;
                            n1 = n1->parent;
                            l1--;
                        }

                        while (n1!=n2) {
                            GetBranchLength (n1,p);
                            p1 += p;
                            GetBranchLength (n2,p);
                            p2 += p;
                            n2 = n2->parent;
                            n1 = n1->parent;
                        }

                        resValue = p1+p2;
                    } else if (n1)
                        if (node2) {
                            if (node1->Equal(node2)) {
                                resValue = 0.0;
                            } else if (node2->Equal (&expectedNumberOfSubs)) {
                                _String bl;
                                GetBranchLength (n1, bl, true);
                                if (bl.sLength) {
                                    DeleteObject(twoIDs);
                                    return new _FString (bl);
                                }
                            }
                        } else {
                            GetBranchLength(n1,resValue);
                        }
                }
                DeleteObject (twoIDs);
            }
        }
    }

    if (isnan (resValue)) {
        return new _MathObject ();
    }

    return new _Constant (resValue);

}
//__________________________________________________________________________________
_PMathObj _TreeTopology::BranchName (_PMathObj p, bool subtree, _PMathObj p2)
{
    _String resString;

    if (p) {
        if (p->ObjectClass()==NUMBER) {
            long res    = p->Value(),
                 count    = -1;

            if (res>=0) {
                DepthWiseT(true);
                while (currentNode) {
                    if (!IsCurrentNodeATip()) {
                        count++;
                    }

                    if (res==count) {
                        if (subtree) {
                            _String         st (128L, true);
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
                            SubTreeString   (st,true,mapMode);
                            st.Finalize();
                            resString = st;
                        } else
                            //resString = travNode->GetName()->Cut(travNode->GetName()->Find ('.')+1,-1);
                        {
                            GetNodeName (currentNode, resString);
                        }
                        break;
                    }
                    DepthWiseT(false);
                    if (IsCurrentNodeTheRoot()) {
                        return new _MathObject();
                    }
                }
            } else {
                count = 0;
                DepthWiseT(true);
                while (currentNode) {
                    count++;
                    DepthWiseT(false);
                }

                _Matrix* branchLengths = new _Matrix (1,count,false,true);
                checkPointer (branchLengths);
                branchLengths->Convert2Formulas ();

                count = 0;
                DepthWiseT(true);

                //long  cutAt = GetName()->sLength+1;

                while (currentNode) {
                    _String       bs;   //(travNode->GetName()->Cut(cutAt,-1));
                    GetNodeName  (currentNode, bs);
                    _FString*     bName = new _FString (bs);
                    checkPointer (bName);
                    _Formula      bNamef (bName,false);

                    branchLengths->StoreFormula (0,count++,bNamef);
                    DepthWiseT(false);
                }
                return branchLengths;
            }
        } else {
            if (p->ObjectClass()==STRING) {
                _List * twoIDs = ((_FString*)p->Compute())->theString->Tokenize(";");


                if (twoIDs->lLength == 2 || twoIDs->lLength == 1) {

                    _String * node1 = (_String*)(*twoIDs) (0),
                              * node2 = twoIDs->lLength>1?(_String*)(*twoIDs) (1):nil;



                    if (twoIDs->lLength == 1) {
                        _AssociativeList * resList = (_AssociativeList*) checkPointer(new _AssociativeList);
                        long            level       = 0,
                                        masterLevel = 0;

                        StepWiseTLevel   (level,true);

                        _String       givenNodeName;

                        while (currentNode) {
                            GetNodeName   (currentNode,givenNodeName);
                            if (givenNodeName.Equal (node1)) {
                                masterLevel = level;
                                resList->MStore(givenNodeName,new _Constant (currentNode->get_num_nodes()));
                                do {
                                    GetNodeName   (currentNode,givenNodeName);
                                    resList->MStore(givenNodeName,new _Constant (1+(currentNode->get_num_nodes()>0)));
                                    StepWiseTLevel(level,false);
                                } while (currentNode && level > masterLevel);

                                break;
                            }
                            StepWiseTLevel(level,false);
                        }
                        if (resList->avl.countitems() == 0)
                            // fail
                        {
                            DeleteObject (resList);
                            return new _MathObject;
                        }
                        return resList;
                    } else {
                        node<long>* n1 = nil,
                                    * n2 = nil;

                        long        l1 = 0,
                                    l2 = 0,
                                    l  = 0;

                        DepthWiseTLevel(l,true);

                        _String      cBranchName;

                        while (currentNode && (!n1 || !n2)) {
                            GetNodeName (currentNode, cBranchName);
                            if (cBranchName.Equal (node1)) {
                                n1 = currentNode;
                                l1 = l;
                            } else if (node2 && cBranchName.Equal (node2)) {
                                n2 = currentNode;
                                l2 = l;
                            }

                            DepthWiseTLevel  (l,false);
                        }

                        if (n1 && n2) {
                            _List prefix,
                                  suffix;

                            while (l1<l2) {
                                GetNodeName              (n2,cBranchName);
                                suffix.AppendNewInstance (cBranchName.makeDynamic());
                                n2 = n2->parent;
                                l2--;
                            }

                            while (l2<l1) {
                                GetNodeName             (n1,cBranchName);
                                prefix.AppendNewInstance (cBranchName.makeDynamic());
                                n1 = n1->parent;
                                l1--;
                            }

                            while (n1!=n2) {
                                GetNodeName     (n2,cBranchName);
                                suffix.AppendNewInstance (cBranchName.makeDynamic());
                                GetNodeName     (n1,cBranchName);
                                prefix.AppendNewInstance (cBranchName.makeDynamic());
                                n2 = n2->parent;
                                n1 = n1->parent;
                            }

                            suffix.Flip();
                            prefix << suffix;
                            return new _Matrix(prefix);
                        } else if (n1) {
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
void _TreeTopology::RerootTreeInternalTraverser (long originator, bool passedRoot, _String&res, long blOption, bool firstTime)
{
    if (passedRoot) {
        SubTreeString (res);
    } else {
        // move to parent now
        node<long>*     myParent = currentNode->get_parent(), *saveCurrent;
        _String t;
        if (myParent->get_parent()) { // not root yet
            res<<'(';
            saveCurrent = currentNode;
            currentNode = myParent;
            RerootTreeInternalTraverser (currentNode->get_child_num(),false,res,blOption);
            for (long i = 1; i<=myParent->get_num_nodes(); i++) {
                if (i==originator) {
                    continue;
                }
                currentNode = myParent->go_down(i);
                res<<',';
                SubTreeString (res,false,blOption);
            }
            res<<')';
            currentNode = saveCurrent;
            if (!firstTime) {
                GetNodeName (currentNode/*myParent*/, t);
                if (!t.startswith(iNodePrefix)) {
                    res<<&t;
                }
            }
            PasteBranchLength (currentNode,res,blOption);
        } else { // passing old root
            // create a new root with >=2 children nodes - this, and one more containing all other children (>=2 of them)
            long count        = 0,
                 rootChildren = theRoot->get_num_nodes();

            if (rootChildren > 2) {
                res << '(';
            }

            node<long>* stashOriginator;

            for (long k = 1; k<=theRoot->get_num_nodes(); k++) {
                currentNode = theRoot->go_down(k);
                if (k==originator) {
                    stashOriginator = currentNode;
                    continue;
                }
                if (count) {
                    res<<',';
                }
                count++;
                SubTreeString (res,false,blOption);
            }
            if (rootChildren > 2) {
                res<<')';
            }

            /*if (stashOriginator->get_num_nodes())
            {
                GetNodeName (stashOriginator, t);
                if (!t.startswith(iNodePrefix))
                    res<<&t;
            }*/
            PasteBranchLength (stashOriginator,res,blOption);
        }
    }
}

//__________________________________________________________________________________
void _TreeTopology::SubTreeString (_String&res, bool allNames,long branchLengths,_AVLListXL* subs)
{
    long    level           = 0,
            lastLevel        = 0,
            j;

    _String t;

    node<long>* saveCurrent = currentNode;

    currentNode = DepthWiseStepTraverserLevel(level,currentNode);

    while (currentNode) {
        if (level>lastLevel) {
            if (lastLevel) {
                res<<',';
            }
            for (j=0; j<level-lastLevel; j++) {
                res<<'(';
            }
        } else if (level<lastLevel) {
            for (j=0; j<lastLevel-level; j++) {
                res<<')';
            }
        } else if (lastLevel) {
            res<<',';
        }

        GetNodeName (currentNode, t);
        if (subs) {
            long mapIdx = subs->Find (&t);
            if (mapIdx >= 0) {
                t = *(_String*)subs->GetXtra (mapIdx);
            }
        }

        lastLevel = level;
        if (!IsCurrentNodeTheRoot()) {
            if (allNames || (!t.startswith (iNodePrefix))) {
                res<< &t;
            }
            PasteBranchLength (currentNode,res,branchLengths);

        }
        currentNode = DepthWiseStepTraverserLevel(level,(node<long>*)nil);
    }

    currentNode = saveCurrent;
}

//__________________________________________________________________________________
void            _TreeTopology::PasteBranchLength (node<long>* currentNode, _String& res, long branchLengths, _Parameter factor)
{
    if (branchLengths!=-1) {
        _String t;
        if (branchLengths==-2) {
            GetBranchValue (currentNode, t);
        } else if (branchLengths==-3) {
            GetBranchLength (currentNode,t);
        } else {
            GetBranchVarValue (currentNode, t, branchLengths);
        }

        if (t.sLength) {
            t = t.toNum()*factor;
            res<<':';
            res<<&t;
        }
    }

}

//__________________________________________________________________________________
void            _TreeTopology::GetBranchLength (node<long> * n, _String& r, bool getBL)
{
    if (getBL) {
        r = empty;
    } else {
        r = compExp->theData[n->in_object];
    }
}

//__________________________________________________________________________________
void            _TreeTopology::GetBranchLength (node<long> * n, _Parameter& r)
{
    r = compExp->theData[n->in_object];
}

//__________________________________________________________________________________
void            _TheTree::GetBranchLength (node<long> * n, _String& r, bool getBL)
{
    if (getBL) {
        bool    mbf;

        _Matrix *mm,
                *fv;

        RetrieveModelComponents(((_CalcNode*)(((BaseRef*)variablePtrs.lData)[n->in_object]))->GetModelIndex(), mm, fv, mbf);

        if (mm && fv) {
            r.CopyDynamicString(mm->BranchLengthExpression(fv,mbf), true);
        } else {
            r = empty;
        }

    } else {
        r = ((_CalcNode*)(((BaseRef*)variablePtrs.lData)[n->in_object]))->BranchLength();
    }
}

//__________________________________________________________________________________
void            _TheTree::GetBranchLength (node<long> * n, _Parameter& r)
{
    r = ((_CalcNode*)(((BaseRef*)variablePtrs.lData)[n->in_object]))->BranchLength();
}

//__________________________________________________________________________________
void            _TreeTopology::GetBranchValue (node<long> *, _String& r)
{
    r = empty;
}

//__________________________________________________________________________________
void            _TheTree::GetBranchValue (node<long> * n, _String& r)
{
    _Parameter t = ((_CalcNode*)(((BaseRef*)variablePtrs.lData)[n->in_object]))->Value();
    if (t != -1.) {
        r = t;
    } else {
        r = empty;
    }
}

//__________________________________________________________________________________
_String*            _TheTree::GetBranchSpec (node<long> * n)
{
    _CalcNode* cn = (_CalcNode*)(((BaseRef*)variablePtrs.lData)[n->in_object]);

    _String * res = new _String (32L, true);
    long      nModel = cn->GetTheModelID();
    if (nModel >= 0) {
        (*res) << '{';
        (*res) << LocateVar(modelMatrixIndices.lData[nModel])->GetName();
    }

    if (iVariables && iVariables->lLength) {
        if (res->sLength) {
            (*res) << ',';
        } else {
            (*res) << '{';
        }

        for (long k=0; k < iVariables->lLength; k+=2) {
            if (k) {
                (*res) << ',';
            }
            _Variable * av = LocateVar (iVariables->lData[k]);
            if (iVariables->lData[k+1] >= 0) {
                (*res) << LocateVar (iVariables->lData[k+1])->GetName();
            } else {
                (*res) << av->GetName();
            }

            (*res) << '=';
            _String   vval (av->Value());
            (*res) << &vval;
        }
    }

    if (dVariables && dVariables->lLength) {
        long m = 0;
        for (long k=0; k < dVariables->lLength; k+=2) {
            if (dVariables->lData[k+1] < 0) {
                if (m) {
                    (*res) << ',';
                } else if (res->sLength) {
                    (*res) << ',';
                } else {
                    (*res) << '{';
                }

                m++;
                _Variable * av = LocateVar (dVariables->lData[k]);
                (*res) << av->GetName();
                (*res) << ":=";
                (*res) << '=';
                _String   * vval = av->GetFormulaString();
                (*res) << vval;
                DeleteObject (vval);
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
void            _TreeTopology::GetBranchVarValue (node<long> *, _String& r, long)
{
    r = empty;
}

//__________________________________________________________________________________
void            _TheTree::GetBranchVarValue (node<long> * n, _String& r, long idx)
{
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

    iNodePrefix = "Node";
    _PMathObj iv = FetchObjectFromVariableByType(&internalNodePrefix,STRING);
    if (iv) {
        iNodePrefix = *((_FString*)iv)->theString;
    }

    if (p&& p->ObjectClass()==STRING) {
        if (rooted == UNROOTED) {
            ReportWarning ("Reroot was called with an unrooted tree. Rerooting was still performed.");
        }

        _String *tNodeN    = (_String*)p->toStr();

        node<long>* rerootNode = FindNodeByName (tNodeN);

        if (rerootNode) { // good node name, can reroot
            if (rerootNode->parent) {
                node <long> * saveCN = rerootNode;
                (*res)<<'('; // opening tree (
                RerootTreeInternalTraverser (rerootNode->get_child_num(),false,*res,-2,true);
                (*res)<<',';
                currentNode = saveCN;
                SubTreeString (*res,0,-2);
                (*res)<<')';
            } else {
                SubTreeString (*res,0,-2);
            }
        }

        DeleteObject (tNodeN);
    } else {
        _String errMsg ("Reroot Tree was passed an invalid branch argument.");
        WarnError (errMsg);
    }

    res->Finalize();
    return new _FString (res);
}
//__________________________________________________________________________________

void    _TheTree::AlignNodes (node<nodeCoord>* theNode)
{
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

node<nodeCoord>* _TheTree::AlignedTipsMapping (bool first, bool respectRoot)
{
    long  k;
    long descendants;
    if (first) {
        treeLayoutVert = 0.0;
        descendants = theRoot->get_num_nodes();
        node<nodeCoord>* aRoot = new node<nodeCoord>; // reslove rootedness here
        aRoot->in_object.varRef = -1;
        aRoot->go_next(); // dumbass initialization for later
        if ((rooted == UNROOTED)||(!respectRoot)) {
            for (k=1; k<=descendants; k++) {
                currentNode = theRoot->go_down(k);
                aRoot->add_node(*AlignedTipsMapping());
            }
            AlignNodes (aRoot);
            return aRoot;
        } else {
            node<nodeCoord>* aChild = new node<nodeCoord>;
            aChild->in_object.varRef = -1;
            if (rooted == ROOTED_LEFT) {
                aRoot->add_node (*aChild);
                for (k=1; k<descendants; k++) {
                    currentNode = theRoot->go_down(k);
                    aChild->add_node(*AlignedTipsMapping());
                }
                currentNode = theRoot->go_down(k);
                aRoot->add_node(*AlignedTipsMapping());
            } else {
                currentNode = theRoot->go_down(1);
                aRoot->add_node(*AlignedTipsMapping());
                for (k=2; k<=descendants; k++) {
                    currentNode = theRoot->go_down(k);
                    aChild->add_node(*AlignedTipsMapping());
                }
                aRoot->add_node (*aChild);
            }
            AlignNodes (aChild);
            AlignNodes (aRoot);
            return      aRoot;
        }
    } else {
        node<nodeCoord>*aNode = new node<nodeCoord>;
        node<long>* saveCurrent = currentNode;
        descendants = saveCurrent->get_num_nodes();
        for (k=1; k<=descendants; k++) {
            currentNode = saveCurrent->go_down(k);
            aNode->add_node(*AlignedTipsMapping());
        }
        if (!descendants) { // terminal node
            aNode->in_object.v = treeLayoutVert;
            aNode->in_object.h = 0;
            treeLayoutVert+=TREE_V_SHIFT;
        } else {
            AlignNodes (aNode);
        }
        aNode->in_object.varRef = saveCurrent->in_object;
        currentNode = saveCurrent;
        return aNode;
    }
}

//__________________________________________________________________________________

#define     BAD_BRANCH_LENGTH  0.000000001

//__________________________________________________________________________________

void _TheTree::ScaledBranchReMapping (node<nodeCoord>* theNode, _Parameter tw)
{
    theNode->in_object.h -= tw;

    long descendants=theNode->get_num_nodes();

    for (long k=1; k<=descendants; k++) {
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
    return empty;
}
//__________________________________________________________________________________

_Parameter       _TheTree::DetermineBranchLengthGivenScalingParameter (long varRef, _String& matchString, char mapMode)
{
    if (mapMode == 3) {
        return 1.;
    }

    _CalcNode * travNode = (_CalcNode*)LocateVar(varRef);

    _Parameter branchLength = BAD_BRANCH_LENGTH;

    if (mapMode==1) {
        return travNode->BranchLength();
    } else if (mapMode==2) {
        branchLength = travNode->Value();
        if (branchLength<=0.0) {
            branchLength = BAD_BRANCH_LENGTH;
        }
    } else {
        long j;
        if (travNode->iVariables)
            for (j=0; j<travNode->iVariables->lLength; j+=2) {
                _Variable* curVar  = LocateVar (travNode->iVariables->lData[j]);
                if (curVar->GetName()->endswith (matchString)) {
                    branchLength = curVar->Compute()->Value();
                    if (branchLength<=0.0) {
                        branchLength = BAD_BRANCH_LENGTH;
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
                        branchLength = BAD_BRANCH_LENGTH;
                    } else {
                        break;
                    }
                }
            }
    }

    return branchLength;
}


//__________________________________________________________________________________

node<nodeCoord>* _TheTree::ScaledBranchMapping (node<nodeCoord>* theParent, _String* scalingParameter, long locDepth, long& depth, char mapMode)
{
    // run a pass of aligned tip mapping then perform one more pass from the root to the children
    // pre-order to remap the length of branches.

    static  _Parameter treeWidth;
    bool    wasRoot = !theParent;

    if (!theParent) {
        theParent = AlignedTipsMapping (true,true);
        theParent->in_object.h = 0.0;
        treeWidth = 0;
    }

    node<nodeCoord>* currentN;
    long descendants = theParent->get_num_nodes(),
         k           = 1,
         j,
         b           = -1;

    _Parameter  branchLength = BAD_BRANCH_LENGTH;

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

    _Parameter  branchLength = BAD_BRANCH_LENGTH;

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
                    newRoot     = AlignedTipsMapping  (true);
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
                vScale      = treeRadius/hScale;
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
                                                      hScale / MAX(BAD_BRANCH_LENGTH,currentNd->in_object.label1));

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

                _String rulerLabel (vScale);

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


void    _TheTree::BuildINodeDependancies (void)
{
    _CalcNode * travNode = DepthWiseTraversal (true);
    leftiNodes.Clear ();
    long        iNodeCounter  = 0,
                leafCounter = 0;

    topLevelNodes.Clear();

    while (travNode) {
        if (IsCurrentNodeATip()) {
            leftiNodes << iNodeCounter;
            leafCounter ++;
        } else {
            iNodeCounter++;
        }
        travNode = DepthWiseTraversal ();
    }
}

//_______________________________________________________________________________________________


void    _TheTree::BuildTopLevelCache (void)
{
    _CalcNode * travNode         = DepthWiseTraversal (true);
    long        iNodeCounter     = 0,
                leafCounter     = 0;

    node <long> *np,
         *np2;

    topLevelNodes.Clear();
    topLevelLeftL.Clear();
    topLevelRightL.Clear();

    // use cBase to count the number of leaves at or below a given node
    // use categoryIndexVars to store left and right leaves

    while (travNode) {
        if (IsCurrentNodeATip()) {
            travNode->categoryIndexVars<<leafCounter;
            travNode->categoryIndexVars<<leafCounter;
            leafCounter ++;
            travNode->cBase = 1;
        } else {
            travNode->cBase = 0;
            for (long k = 0; k < currentNode->nodes.length; k++) {
                _CalcNode *cNode = (_CalcNode*)LocateVar(currentNode->nodes.data[k]->in_object);
                if (k==0) {
                    travNode->categoryIndexVars << cNode->categoryIndexVars[cNode->categoryIndexVars.lLength-2];
                }
                if (k==currentNode->nodes.length-1) {
                    travNode->categoryIndexVars << cNode->categoryIndexVars[cNode->categoryIndexVars.lLength-1];
                }

                travNode->cBase += cNode->cBase;
            }
            travNode->lastState = iNodeCounter;
            iNodeCounter++;
        }
        travNode = DepthWiseTraversal ();
    }

    for (long level = 0; level < theRoot->nodes.length; level++) {
        np = theRoot->nodes.data[level];
        travNode =  (_CalcNode*)LocateVar(np->in_object);
        if (travNode->cBase>1) {
            topLevelNodes << travNode->lastState;
            topLevelLeftL << travNode->categoryIndexVars[travNode->categoryIndexVars.lLength-2];
            topLevelRightL<< travNode->categoryIndexVars[travNode->categoryIndexVars.lLength-1];
            if (travNode->cBase>4*leafCounter/5) {
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

long    _TheTree::CountTreeCategories (void)
{
    categoryVariables.Clear();
    {
        _AVLList           cVA (&categoryVariables);
        ScanForCVariables (cVA);
        cVA.ReorderList   ();
    }
    categoryCount = 1;
    for (unsigned long k=0; k<categoryVariables.lLength; k++) {
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

nodeCoord   _TheTree::TreeTEXRecurse (node<nodeCoord>* currNode, _String&res, _Parameter hScale, _Parameter vScale, long hSize, long vSize)
{
    long  descendants = currNode->get_num_nodes(), vc, hc, hc1, hc2;

    _String t;

    if (descendants==0) { // terminal node
        vc = vSize-currNode->in_object.v*vScale;
        hc = hSize + currNode->in_object.h*hScale;
        t = _String ("\n\\put(")&hc&','&vc&"){\\circle*{2}}";
        res<<&t;
        t = _String ("\n\\put(")&hc+2&','&vc-1&"){\\makebox{\\tiny{";
        res<<&t;
        t = *(LocateVar(currNode->in_object.varRef)->GetName());
        t = t.Cut (t.Find('.')+1,-1);
        res<<&t;
        res<<'}';
        res<<'}';
        res<<'}';
    } else {
        vc = vSize-currNode->in_object.v*vScale;
        hc = hSize + currNode->in_object.h*hScale;
        for (long k = 1; k<=descendants; k++) {
            node<nodeCoord>* childNode = currNode->go_down(k);
            TreeTEXRecurse (childNode, res, hScale, vScale, hSize, vSize);

            long        vcc = vSize-childNode->in_object.v*vScale,
                        hcc = hSize + childNode->in_object.h*hScale;

            t = _String ("\n\\put(")&hc&','&vcc&"){\\line(1,0){"&(hcc-hc)&"}}";
            res<<&t;
            if (k==1) {
                hc1 = vcc;
            } else if (k==descendants) {
                hc2 = vcc;
            }
        }
        t = _String ("\n\\put(")&hc&','&hc2&"){\\line(0,1){"&(hc1-hc2)&"}}";
        res<<&t;
        t = _String ("\n\\put(")&hc&','&vc&"){\\circle{2}}";
        res<<&t;
        if (currNode->get_parent()) {
            t = _String ("\n\\put(")&hc+2&','&vc-1&"){\\makebox{\\tiny{";
            res<<&t;
            t = *(LocateVar(currNode->in_object.varRef)->GetName());
            t = t.Cut (t.Find('.')+1,-1);
            if (t.beginswith ("Node")) {
                t = t.Cut(4,-1);
            }
            res<<&t;
            res<<'}';
            res<<'}';
            res<<'}';
        }
    }

    nodeCoord resC;
    resC.h = hc;
    resC.v = vc;
    return resC;
}

//__________________________________________________________________________________

void    _TheTree::TreePSRecurse (node<nodeCoord>* currNode, _String&res, _Parameter hScale, _Parameter vScale,
                                 long hSize, long vSize, long halfFontSize, long shift, _AssociativeList* outOptions,
                                 char layout, _Parameter * xtra)
{
    long               descendants = currNode->get_num_nodes();
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

    if (currNode->parent) {
        if (currNode->in_object.varRef >=0) {
            varName = (*LocateVar(currNode->in_object.varRef)->GetName());
        }
    } else {
        hc = GetRoot().in_object;
        if (hc >=0) {
            varName = (*LocateVar(hc)->GetName());
        }
        if (layout == 1) {
            res <<  (_String (currNode->in_object.h) & ' ' & _String (currNode->in_object.v) & " translate\n");
        }
    }

    varName.Trim(varName.Find('.')+1,-1);

    _AssociativeList * nodeOptions = nil;
    if (outOptions) {
        nodeOptions = (_AssociativeList *) outOptions->GetByKey (varName, ASSOCIATIVE_LIST);
    }


    _PMathObj nodeLabel     = nodeOptions?nodeOptions->GetByKey (treeOutputLabel,STRING):nil,
              nodeTLabel = nodeOptions?nodeOptions->GetByKey (treeOutputTLabel,STRING):nil;

    if (layout == 1) {
        hcl = currNode->in_object.h*hScale;
        vcl = currNode->in_object.v*hScale;
    } else {
        vcl = vSize-currNode->in_object.v*vScale,
        hcl = hSize+currNode->in_object.h*hScale-shift;
    }

    if (descendants==0 || nodeLabel)
        // terminal node or default label
    {
        t = empty;
        _Parameter myAngle = layout==1?currNode->in_object.label2*DEGREES_PER_RADIAN:0.0;
        if (layout == 1) {
            res << (_String(myAngle) & " rotate\n");
            vc = 0;
            hcl = (vScale + currNode->in_object.bL)*hScale;
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

            if (descendants==0 && t.sLength == 0)
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

        if (descendants==0) {
            currNode->in_object.h = hc-halfFontSize;
        }

        if (layout == 1) {
            res << (_String(-currNode->in_object.label2*DEGREES_PER_RADIAN) & " rotate\n");
        }
    }

    long  minChildHC = 0x0fffffff,
          newV       = 0;

    if (descendants) {
        vc = vSize - currNode->in_object.v*vScale;
        hc = hSize + currNode->in_object.h*hScale-shift;

        nodeCoord childCoord;
        for (long k = 1; k<=descendants; k++) {
            node<nodeCoord>* child = currNode->go_down(k);
            TreePSRecurse (child, res, hScale, (layout==1)?vScale+currNode->in_object.bL:vScale, hSize, vSize,halfFontSize,shift,outOptions,layout,xtra);
            if (k==1) {
                hc1 = layout==1?child->in_object.label2:child->in_object.v;
            }
            if (k==descendants) {
                hc2 = layout==1?child->in_object.label2:child->in_object.v;
            }
        }

        char doVLines = 3;

        for (long k = 1; k<=descendants; k++) {
            node<nodeCoord>* child = currNode->go_down(k);

            if (child->in_object.varRef>=0) {
                t = *LocateVar(child->in_object.varRef)->GetName();
                t.Trim (t.Find('.')+1,-1);
            } else {
                t = empty;
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
                    t = _String(currNode->in_object.label1*hScale) & " 0 moveto\n";
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
                    t = _String(hc) & ' ' & _String (0.5*(hc1+hc2)) & " moveto\n";
                    res<<&t;
                    t = _String(hc) & ' ' & _String (k==1?hc1:hc2) & " lineto\n";
                    res<<&t;
                    res << "stroke\n";
                    doVLines -= k;
                }

            }

            res << "newpath\n";
            if (layout == 1) {
                t = _String(child->in_object.label1*hScale) & " 0 moveto\n";
                res << &t;
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
                    t = empty;
                } else {
                    t = _String(hc - lineWidthInset) & ' ' & _String (child->in_object.v) & " moveto\n";
                    res<<&t;
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

                t = _String(x) & ' ' & _String (y) & " " & halfFontSize  & " 0 360 arc\n";
                res<<&t;
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

                        t = _String(x-0.5*halfFontSize) & ' ' & _String (y-0.5*halfFontSize) & " moveto ";
                        res << &t;
                        t = _String(x+0.5*halfFontSize) & ' ' & _String (y+0.5*halfFontSize) & " lineto\n";
                        res << &t;
                        t = _String(x-0.5*halfFontSize) & ' ' & _String (y+0.5*halfFontSize) & " moveto ";
                        res << &t;
                        t = _String(x+0.5*halfFontSize) & ' ' & _String (y-0.5*halfFontSize) & " lineto\n";
                        res << &t;
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
            currNode->in_object.h = hc;
            currNode->in_object.v = newV/descendants;
        } else {
            currNode->in_object.auxD = (hc2-currNode->in_object.label2)*DEGREES_PER_RADIAN;
            if (currNode->parent) {
                currNode->in_object.v    = (hc1-currNode->in_object.label2)*DEGREES_PER_RADIAN;
            }
        }
    } else {
        currNode->in_object.v = vc;
    }

    if (nodeTLabel) {
        t = *((_FString*)nodeTLabel->Compute())->theString;
        if (t.sLength) {
            _Parameter    nnWidth = 1.+PSStringWidth (t),
                          scF     = 2.*halfFontSize;

            if (layout == 1) {
                res << (_String(currNode->in_object.label2*DEGREES_PER_RADIAN) & " rotate\n");
            } else {
                if (descendants) {
                    nnWidth = (minChildHC - hc)/nnWidth;
                    vcl     = currNode->in_object.v;
                    //if (nnWidth < 2.*halfFontSize)
                    //scF = nnWidth;
                    //else
                    //nnWidth = 1.;
                }
            }


            if (scF != 2.*halfFontSize) {
                res << (_String("/Times-Roman findfont\n")& scF & " scalefont\nsetfont\n");
            }

            res << "newpath\n";
            if (layout == 1) {
                res << (_String(currNode->in_object.label1 * hScale + halfFontSize) & ' ' & _String (-2*halfFontSize/3) & " moveto\n");
            } else {
                if (descendants) {
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
                res << (_String(-currNode->in_object.label2*DEGREES_PER_RADIAN) & " rotate\n");
            }
        }
    }

    if (colorString.sLength) {
        res << "0 0 0 setrgbcolor\n";
    }

    if (currNode->parent == nil && layout == 1) {
        res <<  (_String (-currNode->in_object.h) & ' ' & _String (-currNode->in_object.v) & " translate\n");
    }
}

//__________________________________________________________________________________
_PMathObj _TheTree::TEXTreeString (_PMathObj p)
{
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
            newRoot = AlignedTipsMapping(true);
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

        tipCount = 0;

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

void _TheTree::SetUpMatrices (long categCount)
{
    _CalcNode* travNode;
    categoryCount = categCount>=1?categCount:1;
    travNode = DepthWiseTraversal (TRUE);
    while   (travNode) {
        if (travNode->IsConstant()) {
            travNode->varFlags |= HY_VC_NO_CHECK;
        }
        travNode->ConvertToSimpleMatrix();
        if (categoryCount==1) {
            travNode->matrixCache = nil;
        } else {
            checkPointer(travNode->matrixCache = (_Matrix**)MemAllocate (categoryCount*sizeof(_Matrix*)));
            for (long i=0; i<categoryCount; i++) {
                travNode->matrixCache[i] = nil;
            }
        }
        travNode = DepthWiseTraversal ();
    }
}

//__________________________________________________________________________________

#if USE_SCALING_TO_FIX_UNDERFLOW
void _TheTree::AllocateUnderflowScalers (long sites)
{
    DeleteObject (scalingForUnderflow);
    scalingForUnderflow = new _Matrix (sites, 1, false, true);
    for (long k=0; k<sites; k++) {
        scalingForUnderflow->theData[k] = 1.0;    //genrand_real2()*2.0;
    }
}

//__________________________________________________________________________________

void _TheTree::DeallocateUnderflowScalers (void)
{
    DeleteObject (scalingForUnderflow);
    scalingForUnderflow = nil;
}
#endif


//__________________________________________________________________________________

void _TheTree::CleanUpMatrices (void)
{
    _CalcNode* travNode;
    travNode = DepthWiseTraversal (TRUE);
    if (categoryCount == 1) {
        while   (travNode) {

            // mod 05/03/2003 - uncomment next 5 lines
            // this breaks after ReplicateConstraint or MolecularClock is called
            // WTF?

            travNode->ConvertFromSimpleMatrix();

            if (travNode->referenceNode>=0) {
                travNode->SetRefNode (-1);
                travNode->compExp = nil;
            } else {
                if (travNode->referenceNode < -1) {
                    travNode->SetRefNode (-1);
                }
            }
            if (travNode->compExp) {
                DeleteObject (travNode->compExp);
                travNode->compExp = nil;
            }

            travNode->varFlags &= HY_VC_CLR_NO_CHECK;
            travNode = DepthWiseTraversal ();
        }
    } else {
        while   (travNode) {
            travNode->ConvertFromSimpleMatrix();
            if (travNode->referenceNode>=0) {
                travNode->SetRefNode (-1);
            } else
                for (long i=0; i<categoryCount; i++) {
                    DeleteObject(travNode->matrixCache[i]);
                }

            free (travNode->matrixCache);
            travNode->matrixCache = nil;
            travNode->compExp = nil;
            travNode->varFlags &= HY_VC_CLR_NO_CHECK;
            travNode = DepthWiseTraversal ();
        }
        categoryCount = 1;
    }

}

//__________________________________________________________________________________

void _TheTree::RemoveModel (void)
{
    _CalcNode* travNode;
    travNode = DepthWiseTraversal (TRUE);
    while   (travNode) {
        travNode->RemoveModel();
        travNode = DepthWiseTraversal ();
    }
    categoryCount = 1;
}

//__________________________________________________________________________________

bool _TheTree::FindScalingVariables (_SimpleList& rec)
{
    long i;
    rec.Clear();
    _CalcNode* travNode;
    travNode = StepWiseTraversal (TRUE);
    travNode = StepWiseTraversal ();
    if (travNode) {
        if (travNode->iVariables)
            for (i=1; i<travNode->iVariables->lLength; i+=2)
                if (travNode->iVariables->lData[i]>=0) {
                    rec<<travNode->iVariables->lData[i];
                }

        if (travNode->dVariables)
            for (i=1; i<travNode->dVariables->lLength; i+=2)
                if (travNode->dVariables->lData[i]>=0) {
                    rec<<travNode->dVariables->lData[i];
                }

    }
    if (rec.lLength==0) {
        return false;
    }
    while   (travNode) {
        for (i=0; i<rec.countitems(); i++) {
            if ( ((!travNode->iVariables) || travNode->iVariables->FindStepping(rec.lData[i],2,1)<0)  &&
                    ((!travNode->dVariables) || travNode->dVariables->FindStepping(rec.lData[i],2,1)<0) ) {
                rec.Delete(i);
                if (rec.lLength==0) {
                    break;
                }
                i--;
            }
        }

        if (!( (travNode->iVariables && travNode->iVariables->lLength) || (travNode->dVariables && travNode->dVariables->lLength)
                || (travNode->gVariables && travNode->gVariables->lLength))) {
            rec.Clear();
            return false;
        }

        travNode = StepWiseTraversal ();
    }
    return true;
}

//__________________________________________________________________________________

bool _TheTree::HaveStringBranchLengths (void)
{
    _CalcNode* travNode = DepthWiseTraversal (TRUE);
    while   (travNode && !IsCurrentNodeTheRoot()) {
        if (travNode->Value()<-0.9) {
            return false;
        }
        travNode = DepthWiseTraversal ();
    }
    return true;
}

//__________________________________________________________________________________

void _TheTree::ScanForVariables (_AVLList& l,_AVLList& l2)
{
    _CalcNode* curNode = DepthWiseTraversal (true);
    while (curNode) {
        curNode->ScanForVariables(l,l2);
        curNode = DepthWiseTraversal();
    }
}

//__________________________________________________________________________________

void _TheTree::ScanAndAttachVariables (void)
{
    _CalcNode* curNode = DepthWiseTraversal (true);
    while (curNode) {
        curNode->ScanAndAttachVariables();
        curNode = DepthWiseTraversal();
    }
}

//__________________________________________________________________________________

void _TheTree::ScanForDVariables (_AVLList& l,_AVLList& l2)
{
    _CalcNode* curNode = DepthWiseTraversal (true);
    while (curNode) {
        curNode->ScanForDVariables(l,l2);
        curNode = DepthWiseTraversal();
    }
}



//__________________________________________________________________________________

void _TheTree::ScanForGVariables (_AVLList& li, _AVLList& ld)
{
    _CalcNode*  curNode = DepthWiseTraversal (true);
    _SimpleList cL;
    _AVLList    cLL (&cL);

    while (curNode) {
        _Matrix *modelM = curNode->GetModelMatrix();
        if (modelM && cLL.Find(modelM) < 0) {
            _SimpleList temp;
            {
                _AVLList tempA (&temp);
                (modelM)->ScanForVariables(tempA, true);
                tempA.ReorderList();
            }
            for (long i=0; i<temp.lLength; i++) {
                long p = temp.lData[i];
                _Variable* v = LocateVar (p);
                if (v&&v->IsGlobal()) {
                    if(v->IsIndependent()) {
                        li.Insert ((BaseRef)p);
                    } else {
                        ld.Insert ((BaseRef)p);
                    }
                }
            }
            cLL.Insert (modelM);
        }
        curNode -> ScanForGVariables(li,ld);
        curNode = DepthWiseTraversal();
    }
}

//__________________________________________________________________________________

void _TheTree::ScanForCVariables (_AVLList& lcat)
{
    _CalcNode* curNode = DepthWiseTraversal (true);
    while (curNode) {
        for (long i = curNode->categoryVariables.lLength-1; i>=0; i--) {
            lcat.Insert ((BaseRef)curNode->categoryVariables(i));
        }

        curNode = DepthWiseTraversal();

    }
}

//__________________________________________________________________________________

bool _TheTree::HasChanged (void)
{
    _CalcNode* curNode = StepWiseTraversal (true);
    while (curNode) {
        if (curNode->HasChanged()) {
            return true;
        }
        curNode = StepWiseTraversal();
    }
    return false;
}

//__________________________________________________________________________________

bool _TheTree::HasChanged2 (void)
{
    for (long k = 0; k < categoryVariables.lLength;  k++)
        if (((_CategoryVariable*)LocateVar(categoryVariables.lData[k]))->HaveParametersChanged()) {
            return true;
        }
    _CalcNode* curNode = StepWiseTraversal (true);
    while (curNode) {
        if (curNode->_VariableContainer::HasChanged()) {
            return true;
        }
        curNode = StepWiseTraversal();
    }
    return false;
}

//_______________________________________________________________________________________________

_Parameter  _TheTree::Probij  (long i, long j, _CalcNode * childNode)
{
    // assumed that ExpMatrix has already been called for that node
    // range checking is implicit in the call of matrix (.,.) function
    if (childNode) {
        if(!childNode->GetCompExp()) {
            childNode->RecomputeMatrix();
        }
//      if(childNode->GetCompExp())
        return  (*childNode->GetCompExp())(i,j);
    }
    return 0;
}

//_______________________________________________________________________________________________
void     _TheTree::InitializeTreeFrequencies (_Matrix *mx, bool setDim)
// this will take the  matrix of frequencies and
// 1) use its dimensions to initialize tree freq holders
// 2) place global frequencies into the tree holder for later use by the pruning algo
// must be called before any tree pruning computations are started
{
    long vecDim = mx->GetHDim()*mx->GetVDim();
    // theModel = mx;

    if  (setDim) {
        SetTreeCodeBase (vecDim);
    } else
        for (long i=0; i<vecDim; i++) {
            theProbs[i] = mx->theData[i];
        }
}
//_______________________________________________________________________________________________
void     _TheTree::SetTreeCodeBase (long b)
// this will take the  matrix of frequencies and
// 1) use its dimensions to initialize tree freq holders
// 2) place global frequencies into the tree holder for later use by the pruning algo
// must be called before any tree pruning computations are started
{
    SetCodeBase (b);
    if (marginalLikelihoodCache) {
        free (marginalLikelihoodCache);
        marginalLikelihoodCache = nil;
    }
    if (cBase>0)
        marginalLikelihoodCache =
            (_Parameter*)MemAllocate ((flatNodes.lLength+flatLeaves.lLength)*sizeof (_Parameter)*cBase*systemCPUCount);

    _CalcNode*  travNode = StepWiseTraversal (TRUE);
    while (travNode) {
        travNode -> SetCodeBase (b);
        travNode = StepWiseTraversal();
    }

}

//_______________________________________________________________________________________________
long     _TheTree::IsLinkedToALF (long& pid)
{
    for (long lfID = 0; lfID < likeFuncList.lLength; lfID ++)
        if (likeFuncList.lData[lfID] && (pid = ((_LikelihoodFunction*)likeFuncList(lfID))->DependOnTree (*GetName())) >= 0) {
            return lfID;
        }
    return -1;
}

#ifdef __MP__
//_______________________________________________________________________________________________

void*   MatrixUpdateFunction (void* arg)
{
    ThreadMatrixTask* theTask = (ThreadMatrixTask*)arg;
    for (long k=theTask->startAt; k<theTask->endAt; k++) {
        ((_CalcNode*)(theTask->updateCN->lData[k]))->RecomputeMatrix(theTask->cID, theTask->tcat);
    }
    return nil;
}
#endif

//_______________________________________________________________________________________________

_Parameter   _TheTree::ReleafTreeAndCheck (_DataSetFilter* dsf, long index, bool cache, long categID)
// set leaf value and reexp if needed
// for first entry into a datafilter
{

#if defined     __MP__
    if (systemCPUCount>1) {
        ThreadMatrixUpdate (categID, cache);
    } else
#endif
        SerialMatrixUpdate (categID, cache);

    if (cache) {
        MatrixCacheUpdate ();
    }

    if (flatLeaves.lLength == 1) {
        return ReleafTreeDegenerate (dsf,index);
    }

    if (cache) {
        return ThreadReleafTreeCache (dsf,index,-1,0,flatLeaves.lLength-1,categID>=0?categID:0);
    } else {
        return ReleafTree (dsf,index,-1,0,flatLeaves.lLength-1);
    }
}

//_______________________________________________________________________________________________

void        _TheTree::ThreadMatrixUpdate (long categID, bool cache)
{
#ifdef __MP__
    _CalcNode       *travNode;
    node <long>     *nodeChild;
    _SimpleList     *taintedNodes = new _SimpleList;

    for (long nodeCount = 0; nodeCount<flatLeaves.lLength; nodeCount++) {
        travNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
        bool reexpnt = travNode->NeedToExponentiate(categID);
        if (reexpnt&&travNode->GetModelMatrix()) {
            (*taintedNodes) << (long)travNode;
            if (cache) {
                travNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                            [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);
                travNode->cBase = -1;
            }
        } else if (categID>=0) {
            travNode->SetCompMatrix (categID);
        }

    }

    for (long nodeCount=0; nodeCount<flatTree.lLength; nodeCount++) {
        travNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[nodeCount]);
        bool reexpnt = travNode->NeedToExponentiate(categID);
        if (reexpnt&&travNode->GetModelMatrix()) {
            (*taintedNodes) << (long)travNode;
            if (cache) {
                travNode->cBase = -1;
            }
        } else if (categID>=0) {
            travNode->SetCompMatrix (categID);
        }

        if (cache&(travNode->cBase==-1)) {
            nodeChild = ((node <long>*)(flatNodes.lData[nodeCount]))->parent;
            if (nodeChild) {
                travNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->in_object]);
                travNode->cBase = -1;
            }
        }
    }

    if ((*taintedNodes).lLength>1) {
        long tStep = (*taintedNodes).lLength/systemCPUCount,
             threadCount = 0,
             k,
             errCode;
        if (tStep>0) {
            threadCount = systemCPUCount-1;
        } else {
            tStep = 1;
            threadCount = (*taintedNodes).lLength-1;
        }
        matrixTasks   = new ThreadMatrixTask[threadCount];
        matrixThreads = new pthread_t [threadCount];

        for (k=0; k<threadCount; k++) {
            matrixTasks[k].cID = categID;
            matrixTasks[k].tcat = categoryCount;
            matrixTasks[k].startAt = tStep*(k+1);
            matrixTasks[k].endAt = tStep*(k+2);
            if (k==threadCount-1) {
                matrixTasks[k].endAt = (*taintedNodes).lLength;
            }
            matrixTasks[k].updateCN = taintedNodes;

            if ( pthread_create( matrixThreads+k, NULL, MatrixUpdateFunction ,(void*)(matrixTasks+k)))
            {
                FlagError("Failed to initialize a POSIX thread in ReleafTreeAndCheck()");
                exit(1);
            }
        }

        for (k=0; k<tStep; k++) {
            ((_CalcNode*)((*taintedNodes).lData[k]))->RecomputeMatrix(categID, categoryCount);
        }

        for (k=0; k<threadCount; k++)
            if ( (errCode = pthread_join ( matrixThreads[k], NULL )) ) {
                FlagError(_String("Failed to join POSIX threads in ReleafTreeAndCheck(). Error Code=")&errCode);
                exit(1);
            }

        delete matrixTasks;
        delete matrixThreads;
        matrixTasks = nil;
    } else if ((*taintedNodes).lLength==1) {
        ((_CalcNode*)((*taintedNodes).lData[0]))->RecomputeMatrix(categID, categoryCount);
    }

    delete taintedNodes;
#endif
}
//_______________________________________________________________________________________________

void        _TheTree::SerialMatrixUpdate (long categID, bool cache)
{
    _CalcNode  *travNode;
    node <long>*nodeChild;

    for (long nodeCount = 0; nodeCount<flatLeaves.lLength; nodeCount++) {
        travNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
        bool reexpnt = travNode->NeedToExponentiate(categID);
        // flag all the nodes up the chain - to the root

        if (reexpnt&&travNode->GetModelMatrix()) {
            travNode->RecomputeMatrix(categID, categoryCount);
            if (cache) {
                travNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                            [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);
                travNode->cBase = -1;
            }
        } else if (categID>=0) {
            travNode->SetCompMatrix (categID);
        }

    }

    {
        //for BCC
        for (long nodeCount=0; nodeCount<flatTree.lLength; nodeCount++) {
            travNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[nodeCount]);
            bool reexpnt = travNode->NeedToExponentiate(categID);
            // flag all the nodes up the chain - to the root

            if (reexpnt&&travNode->GetModelMatrix()) {
                travNode->RecomputeMatrix(categID, categoryCount);
                if (cache) {
                    travNode->cBase = -1;
                }
            } else if (categID>=0) {
                travNode->SetCompMatrix (categID);
            }

            if (cache&(travNode->cBase==-1)) {
                nodeChild = ((node <long>*)(flatNodes.lData[nodeCount]))->parent;
                if (nodeChild) {
                    travNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->in_object]);
                    travNode->cBase = -1;
                }
            }
        }
    }
}

//_______________________________________________________________________________________________

void        _TheTree::MatrixCacheUpdate (void)
{
    long c = 0, off = 1;
    _CalcNode * travNode;
    for (long nodeCount=0; nodeCount<topLevelNodes.lLength-1; nodeCount++, off<<=1) {
        travNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[topLevelNodes.lData[nodeCount]]);
        if (travNode->cBase<=-1) {
            c |= off;
        }
    }
    topLevelNodes.lData[topLevelNodes.lLength-1] = c;
    for (long nodeCount2=0; nodeCount2<flatTree.lLength; nodeCount2++) {
        travNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[nodeCount2]);
        travNode->cBase = cBase;
    }
}


//_______________________________________________________________________________________________

_Parameter   _TheTree::ReleafTreeAndCheckChar4 (_DataSetFilter* dsf, long index, bool cache, long categID)
{

    long    nodeCount  = 0,
            f;

    _Parameter *mCache = marginalLikelihoodCache;

    if (dsf->IsNormalFilter()) {
        char       *thisState = dsf->GetColumn(index);

        for (; nodeCount<flatLeaves.lLength; nodeCount++,mCache+=4) {
            _CalcNode * travNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
            f                    = dsf->theNodeMap.lData[nodeCount];
            long* cCache         = dsf->conversionCache.lData+(thisState[f]-40)*5;
#if USE_SCALING_TO_FIX_UNDERFLOW
            _Parameter scalingFactor = scalingForUnderflow->theData[index];
            mCache[0] = travNode->theProbs[0] = *(cCache++)*scalingFactor;
            mCache[1] = travNode->theProbs[1] = *(cCache++)*scalingFactor;
            mCache[2] = travNode->theProbs[2] = *(cCache++)*scalingFactor;
            mCache[3] = travNode->theProbs[3] = *(cCache++)*scalingFactor;
#else
            mCache[0] = travNode->theProbs[0] = *(cCache++);
            mCache[1] = travNode->theProbs[1] = *(cCache++);
            mCache[2] = travNode->theProbs[2] = *(cCache++);
            mCache[3] = travNode->theProbs[3] = *(cCache++);
#endif
            nodeStates[nodeCount] = travNode->lastState = *cCache;
        }
    } else {
        _DataSetFilterNumeric * dsfN = (_DataSetFilterNumeric*)dsf;
        for (; nodeCount<flatLeaves.lLength; nodeCount++,mCache+=4) {
            _CalcNode * travNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
            //_Parameter*        pv = dsfN->getProbabilityVector (dsf->theNodeMap.lData[nodeCount],index);
            _Parameter*          pv = (categID<0)
                                      ?(dsfN->probabilityVectors.theData + dsf->theNodeMap.lData[nodeCount]*dsfN->shifter + index*4)
                                      :(dsfN->categoryShifter*categID + dsfN->probabilityVectors.theData + dsf->theNodeMap.lData[nodeCount]*dsfN->shifter + index*4)
                                      ;
            mCache[0] = travNode->theProbs[0] = pv[0];
            mCache[1] = travNode->theProbs[1] = pv[1];
            mCache[2] = travNode->theProbs[2] = pv[2];
            mCache[3] = travNode->theProbs[3] = pv[3];
            nodeStates[nodeCount] = travNode->lastState = -1;
        }
    }

    if (flatLeaves.lLength==1) {
        _CalcNode* travNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[theRoot->nodes.data[0]->in_object]);
        bool        reexpnt = travNode->NeedToExponentiate(categID);

        if (reexpnt&&travNode->GetModelMatrix()) {
            travNode->RecomputeMatrix(categID, categoryCount);
        } else if (categID>=0) {
            travNode->SetCompMatrix (categID);
        }
        return ReleafTreeChar4Degenerate (dsf,index);
    }
    if (cache) {
        PruneTreeChar4Cache(categID);
        return ThreadReleafTreeChar4 (dsf,index,-1,0,flatLeaves.lLength-1,categID<0?0:categID);
        //return res;
    }

    return PruneTreeChar4(categID);
}
//_______________________________________________________________________________________________

_Parameter   _TheTree::ReleafTreeDegenerate (_DataSetFilter* dsf, long index)
// set leaf value and re-prune w/o reexping
// for subsequent entries into a datafilter
{
    _CalcNode* rt  = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[theRoot->in_object]),
               * tip = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[theRoot->nodes.data[0]->in_object]);

    _Parameter reslt = 0.;
    // sum over one branch in the direction from the root to the single leaf

    long    state1 = dsf->Translate2Frequencies ((*dsf)(index,0), rt->theProbs, true),
            state2 = dsf->Translate2Frequencies ((*dsf)(index,1), tip->theProbs, true);

    // now perform one loop
    _Parameter* fastIdx =  tip->GetCompExp()->theData;
    //*nodeProbs;
    // 4 cases are possible:
    // 1-1
    // 1-many
    // many-1
    // many-many
    if (state1>=0 && state2 >= 0) { // 1-1
        reslt = theProbs[state1]*fastIdx[state1*cBase+state2];
    } else if (state1 >= 0) { // 1-many
        _Parameter tmp = 0.;

        fastIdx += state1*cBase;

        for (long i=0; i<cBase; i++) {
            tmp += fastIdx[i] * tip->theProbs[i];
        }

        reslt = theProbs[state1] * tmp;
    } else if (state2 >= 0) { // many to 1
        fastIdx  += state2;

        for (long i=0; i<cBase; i++, fastIdx+=cBase) {
            reslt += rt->theProbs[i] * *fastIdx * theProbs[i];
        }

    } else // many to many
        for (long i=0; i<cBase; i++) {
            _Parameter tmp = 0.0;
            for (long j=0; j<cBase; j++,fastIdx++) {
                tmp += *fastIdx * tip->theProbs[j];
            }

            reslt += tmp*rt->theProbs[i] * theProbs[i];
        }

    return reslt <= 0.0 ? ALMOST_ZERO : reslt;

    /*_CalcNode* rt  = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[theRoot->in_object]),
             * tip = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[theRoot->nodes.data[0]->in_object]);

    _Parameter reslt, tmp;
    // sum over one branch in the direction from the root to the single leaf
    reslt = dsf->Translate2Frequencies ((*dsf)(index,0), rt->theProbs, true);
    tmp   = dsf->Translate2Frequencies ((*dsf)(index,1), tip->theProbs, true);
    // now perform one loop
    _Matrix* theTP = tip->GetCompExp();
    _Parameter* fastIdx = theTP->theData, *nodeProbs;
    // 4 cases are possible:
    // 1-1
    // 1-many
    // many-1
    // many-many
    if ((reslt>=0.0)&&(tmp>=0.0)) // 1-1
    {
        tmp = fastIdx[(long)reslt*cBase+(long)tmp];
        reslt = theProbs[(long)reslt]*tmp;
    }
    else
        if ((reslt>=0.0)&&(tmp<0.0)) // 1-many
        {
            tmp = 0;
            fastIdx += (long)reslt*cBase;
            nodeProbs = tip->GetProbs();
            for (long i=0; i<cBase; i++, fastIdx++,nodeProbs++)
            {
                tmp+=*fastIdx* *nodeProbs;
            }
            reslt = theProbs[(long)reslt]*tmp;
        }
        else
            if ((reslt<0.0)&&(tmp>=0.0)) // many to 1
            {
                fastIdx += (long)tmp;
                tmp = 0;
                nodeProbs  = rt->GetProbs();
                long i;
                for (i=0; i<cBase; i++, fastIdx+=cBase,nodeProbs++)
                {
                    *nodeProbs *=*fastIdx;
                }

                nodeProbs-=cBase;
                reslt = 0;
                fastIdx = GetProbs();
                for (i=0; i<cBase; i++, nodeProbs++, fastIdx++)
                {
                    reslt+=*fastIdx * *nodeProbs;
                }
            }
            else // many to many
            {
                nodeProbs  = rt->GetProbs();
                long i,j;
                for (i=0; i<cBase; i++, nodeProbs++)
                {
                    tmp = 0;
                    for (j=0; j<cBase; j++,fastIdx++)
                        tmp += *fastIdx* tip->GetProbs(j);
                    *nodeProbs*=tmp;
                }

                reslt = 0;
                fastIdx = GetProbs();
                nodeProbs-=cBase;
                for (i=0; i<cBase; i++, nodeProbs++, fastIdx++)
                {
                    reslt+=*fastIdx * *nodeProbs;
                }
            }
    if ((reslt<0.0)||(reslt==.0)) reslt = ALMOST_ZERO;
    return reslt;*/
}

//_______________________________________________________________________________________________

_Parameter   _TheTree::ReleafTreeCharDegenerate (_DataSetFilter* dsf, long index)
// set leaf value and re-prune w/o reexping
// for subsequent entries into a datafilter
{

    _CalcNode* rt  = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[theRoot->in_object]),
               * tip = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[theRoot->nodes.data[0]->in_object]);

    _Parameter reslt = 0.;
    // sum over one branch in the direction from the root to the single leaf

    char       *thisState = dsf->GetColumn(index);

    long       state1 = dsf->LookupConversion (thisState[dsf->theNodeMap.lData[0]], rt->theProbs),
               state2 = dsf->LookupConversion (thisState[dsf->theNodeMap.lData[1]], tip->theProbs);

    // now perform one loop
    _Parameter* fastIdx =  tip->GetCompExp()->theData;
    //*nodeProbs;
    // 4 cases are possible:
    // 1-1
    // 1-many
    // many-1
    // many-many
    if (state1>=0 && state2 >= 0) { // 1-1
        reslt = theProbs[state1]*fastIdx[state1*cBase+state2];
    } else if (state1 >= 0) { // 1-many
        _Parameter tmp = 0.;

        fastIdx += state1*cBase;

        for (long i=0; i<cBase; i++) {
            tmp += fastIdx[i] * tip->theProbs[i];
        }

        reslt = theProbs[state1] * tmp;
    } else if (state2 >= 0) { // many to 1
        fastIdx  += state2;

        for (long i=0; i<cBase; i++, fastIdx+=cBase) {
            reslt += rt->theProbs[i] * *fastIdx * theProbs[i];
        }

    } else // many to many
        for (long i=0; i<cBase; i++) {
            _Parameter tmp = 0.0;
            for (long j=0; j<cBase; j++,fastIdx++) {
                tmp += *fastIdx * tip->theProbs[j];
            }

            reslt += tmp*rt->theProbs[i] * theProbs[i];
        }

    return reslt <= 0.0 ? ALMOST_ZERO : reslt;

    /*_CalcNode* rt  = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[theRoot->in_object]),
             * tip = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[theRoot->nodes.data[0]->in_object]);

    _Parameter reslt, tmp;
    // sum over one branch in the direction from the root to the single leaf

    char  *thisState = dsf->GetColumn(index);

    reslt = dsf->LookupConversion (thisState[dsf->theNodeMap.lData[0]], rt->theProbs);
    tmp   = dsf->LookupConversion (thisState[dsf->theNodeMap.lData[1]], tip->theProbs);
    //reslt = dsf->Translate2Frequencies ((*dsf)(index,0), rt->theProbs, true);
    //tmp   = dsf->Translate2Frequencies ((*dsf)(index,1), tip->theProbs, true);

    // now perform one loop
    _Matrix   * theTP = tip->compExp;
    _Parameter* fastIdx = theTP->theData,
              * nodeProbs;
    // 4 cases are possible:
    // 1-1
    // 1-many
    // many-1
    // many-many
    if ((reslt>=0.0)&&(tmp>=0.0)) // 1-1
    {
        tmp = fastIdx[(long)reslt*cBase+(long)tmp];
        reslt = theProbs[(long)reslt]*tmp;
    }
    else
        if ((reslt>=0.0)&&(tmp<0.0)) // 1-many
        {
            tmp = 0;
            fastIdx += (long)reslt*cBase;
            nodeProbs = tip->theProbs;
            for (long i=0; i<cBase; i++, fastIdx++,nodeProbs++)
            {
                tmp+=*fastIdx* *nodeProbs;
            }
            reslt = theProbs[(long)reslt]*tmp;
        }
        else
            if ((reslt<0.0)&&(tmp>=0.0)) // many to 1
            {
                fastIdx += (long)tmp;
                tmp = 0;
                nodeProbs  = rt->theProbs;
                long i;
                for (i=0; i<cBase; i++, fastIdx+=cBase,nodeProbs++)
                {
                    *nodeProbs *=*fastIdx;
                }

                nodeProbs-=cBase;
                reslt = 0;
                fastIdx = theProbs;
                for (i=0; i<cBase; i++, nodeProbs++, fastIdx++)
                {
                    reslt+=*fastIdx * *nodeProbs;
                }
            }
            else // many to many
            {
                nodeProbs  = rt->theProbs;
                long i,j;
                for (i=0; i<cBase; i++, nodeProbs++)
                {
                    tmp = 0;
                    for (j=0; j<cBase; j++,fastIdx++)
                        tmp += *fastIdx* tip->theProbs[j];
                    *nodeProbs*=tmp;
                }

                reslt = 0;
                fastIdx = theProbs;
                nodeProbs-=cBase;
                for (i=0; i<cBase; i++, nodeProbs++, fastIdx++)
                {
                    reslt+=*fastIdx * *nodeProbs;
                }
            }

    if ((reslt<0.0)||(reslt==.0)) reslt = ALMOST_ZERO;
    return reslt;*/
}

//_______________________________________________________________________________________________

_Parameter   _TheTree::ReleafTree (_DataSetFilter* dsf, long index, long lastIndex, long startLeaf, long endLeaf)
// set leaf value and re-prune w/o reexping
// for subsequent entries into a datafilter
{

    long nodeCount;
    _CalcNode* travNode,* theChildNode ;
    _Parameter* fastIndex,*theProbbs, *stopper;
    node<long>* nodeChild;

    if (dsf->GetUnitLength() == 3) {
        char *c1, *c2, *c3,
             *o1, *o2, *o3;

        c1 = dsf->GetColumn (3*index);
        c2 = dsf->GetColumn (3*index+1);
        c3 = dsf->GetColumn (3*index+2);

        if (lastIndex>=0) {
            o1 = dsf->GetColumn (3*lastIndex);
            o2 = dsf->GetColumn (3*lastIndex+1);
            o3 = dsf->GetColumn (3*lastIndex+2);
        }

        long   ccount  = dsf->conversionCache.lData[0],
               ccount2 = ccount*ccount,
               *ccodes  = dsf->conversionCache.lData+1,
                *tcodes  = dsf->conversionCache.lData+89;

        for (nodeCount = startLeaf; nodeCount<=endLeaf; nodeCount++) {
            long       nMap = dsf->theNodeMap.lData[nodeCount];
            travNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
            char       A = c1[nMap],
                       B = c2[nMap],
                       C = c3[nMap];

            if ((lastIndex<0)||(A!=o1[nMap])||(B!=o2[nMap])||(C!=o3[nMap])) {
                A = ccodes[A-40];
                B = ccodes[B-40];
                C = ccodes[C-40];

                if ((A==-1)||(B==-1)||(C==-1)) {
                    travNode->lastState = dsf->Translate2Frequencies ((*dsf)(index,nodeCount), travNode->theProbs, true);
                } else {
                    travNode->lastState = tcodes[A*ccount2+B*ccount+C];
                    if (travNode->lastState < 0) {
                        travNode->lastState = dsf->Translate2Frequencies ((*dsf)(index,nodeCount), travNode->theProbs, true);
                    } else {
                        for (nMap=0; nMap<travNode->lastState; nMap++) {
                            travNode->theProbs [nMap] = 0.0;
                        }
                        travNode->theProbs [nMap++] = 1.0;
                        for (; nMap<cBase; nMap++) {
                            travNode->theProbs [nMap] = 0.0;
                        }
                    }
                }
                theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                                [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);
                theChildNode->cBase = -1;
            }
        }
    } else {
        if (dsf->GetUnitLength() == 2) {
            char *c1, *c2,
                 *o1, *o2;

            c1 = dsf->GetColumn (2*index);
            c2 = dsf->GetColumn (2*index+1);

            if (lastIndex>=0) {
                o1 = dsf->GetColumn (2*lastIndex);
                o2 = dsf->GetColumn (2*lastIndex+1);
            }

            long   ccount  = dsf->conversionCache.lData[0],
                   *ccodes  = dsf->conversionCache.lData+1,
                    *tcodes  = dsf->conversionCache.lData+89;

            for (nodeCount = startLeaf; nodeCount<=endLeaf; nodeCount++) {
                long       nMap = dsf->theNodeMap.lData[nodeCount];
                travNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
                char       A = c1[nMap],
                           B = c2[nMap];

                if ((lastIndex<0)||(A!=o1[nMap])||(B!=o2[nMap])) {
                    A = ccodes[A-40];
                    B = ccodes[B-40];

                    if ((A==-1)||(B==-1)) {
                        travNode->lastState = dsf->Translate2Frequencies ((*dsf)(index,nodeCount), travNode->theProbs, true);
                    } else {
                        travNode->lastState = tcodes[A*ccount+B];
                        if (travNode->lastState < 0) {
                            travNode->lastState = dsf->Translate2Frequencies ((*dsf)(index,nodeCount), travNode->theProbs, true);
                        } else {
                            for (nMap=0; nMap<travNode->lastState; nMap++) {
                                travNode->theProbs [nMap] = 0.0;
                            }
                            travNode->theProbs [nMap++] = 1;
                            for (; nMap<cBase; nMap++) {
                                travNode->theProbs [nMap] = 0.0;
                            }
                        }
                    }
                    theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                                    [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);
                    theChildNode->cBase = -1;
                }
            }
        } else {
            for (nodeCount = startLeaf; nodeCount<=endLeaf; nodeCount++) {
                travNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
                if ((lastIndex<0)||(!dsf->CompareTwoSites(lastIndex,index,nodeCount))) {
                    travNode->lastState = dsf->Translate2Frequencies ((*dsf)(index,nodeCount), travNode->theProbs, true);
                    theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                                    [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);
                    theChildNode->cBase = -1;
                }
            }
        }
    }

    for (nodeCount = leftiNodes.lData[startLeaf]; nodeCount<flatTree.lLength; nodeCount++) {
        theChildNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[nodeCount]);
        // we can immediately place the changes in the parent node
        if (theChildNode->cBase<0) {
            nodeChild = (node<long>*)flatNodes.lData[nodeCount];
            for (long i=0; i<cBase; i++) {
                theChildNode->theProbs[i] = LIKELIHOOD_SCALER;
            }

            for (long k=0; k<nodeChild->nodes.length; k++) {
                travNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->nodes.data[k]->in_object]);
                fastIndex = travNode->compExp->theData;
                theProbbs = travNode->theProbs;
                long      nZ = travNode->lastState;

                if (nZ>=0) {
                    _Parameter a = theProbbs[nZ];
                    theProbbs = theChildNode->theProbs;
                    fastIndex+=nZ;
                    stopper = theProbbs+cBase-cBase%4;
                    for (; theProbbs!=stopper;) {
                        *(theProbbs++)*=a**fastIndex;
                        fastIndex += cBase;
                        *(theProbbs++)*=a**fastIndex;
                        fastIndex += cBase;
                        *(theProbbs++)*=a**fastIndex;
                        fastIndex += cBase;
                        *(theProbbs++)*=a**fastIndex;
                        fastIndex += cBase;
                    }
                    switch (cBase%4) {
                    case 1: {
                        *theProbbs*=a**fastIndex;
                        break;
                    }
                    case 2: {
                        *(theProbbs++)*=a**fastIndex;
                        fastIndex += cBase;
                        *theProbbs*=a**fastIndex;
                        break;
                    }
                    case 3: {
                        *(theProbbs++)*=a**fastIndex;
                        fastIndex += cBase;
                        *(theProbbs++)*=a**fastIndex;
                        fastIndex += cBase;
                        *theProbbs*=a**fastIndex;
                        break;
                    }
                    }
                } else {
                    for (long i=0; i<cBase; i++) {
                        _Parameter tmp = 0.0;
                        stopper = theProbbs+cBase-cBase%4;
                        // loop unrolled to depth 4
                        for (; theProbbs!=stopper;) {
                            tmp += *(theProbbs++)**(fastIndex++);
                            tmp += *(theProbbs++)**(fastIndex++);
                            tmp += *(theProbbs++)**(fastIndex++);
                            tmp += *(theProbbs++)**(fastIndex++);
                        }
                        switch (cBase%4) {
                        case 1: {
                            tmp += *(theProbbs)**(fastIndex++);
                            break;
                        }
                        case 2: {
                            tmp += *(theProbbs++)**(fastIndex++);
                            tmp += *(theProbbs)**(fastIndex++);
                            break;
                        }
                        case 3: {
                            tmp += *(theProbbs++)**(fastIndex++);
                            tmp += *(theProbbs++)**(fastIndex++);
                            tmp += *(theProbbs)**(fastIndex++);
                            break;
                        }
                        }
                        theChildNode->theProbs[i]*=tmp;
                        theProbbs = travNode->theProbs;
                    }
                }
            }
            theChildNode->cBase = cBase;
            nodeChild = nodeChild->parent;
            if (nodeChild) {
                ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->in_object])->cBase = -1;
            }
        }
    }
    _Parameter result = 0;

    for (long i=0; i<cBase; i++) {
        result+= theProbs[i]*theChildNode->theProbs[i];
    }

    if (result<=0.0) {
        /*_Matrix stashZeroValues (cBase, flatTree.lLength, false, true);

        _List   labels;

        for (nodeCount = 0; nodeCount<flatTree.lLength; nodeCount++)
        {
            theChildNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[nodeCount]);
            labels <<  theChildNode->GetName();
            for (long ii=0; ii<cBase; ii++)
                stashZeroValues.Store (ii,nodeCount, theChildNode->theProbs[ii]);
        }

        _String cwt ("Report for site ");
        cwt = cwt & index;
        _HYChartWindow *chart = new _HYChartWindow(cwt, labels, stashZeroValues);
        chart->BringToFront();*/

        return ALMOST_ZERO;
    }

    return result;
}


//_______________________________________________________________________________________________

_Parameter   _TheTree::ReleafTreeCache (_DataSetFilter* dsf, long index, long lastIndex, long startLeaf, long endLeaf, long position)
// set leaf value and re-prune w/o reexping
// for subsequent entries into a datafilter
{

    long nodeCount;
    _CalcNode* travNode,* theChildNode ;
    _Parameter* fastIndex,*theProbbs, *stopper;
    node<long>* nodeChild;

    if (dsf->GetUnitLength() == 3) {
        char *c1, *c2, *c3,
             *o1, *o2, *o3;

        c1 = dsf->GetColumn (3*index);
        c2 = dsf->GetColumn (3*index+1);
        c3 = dsf->GetColumn (3*index+2);

        if (lastIndex>=0) {
            o1 = dsf->GetColumn (3*lastIndex);
            o2 = dsf->GetColumn (3*lastIndex+1);
            o3 = dsf->GetColumn (3*lastIndex+2);
        }

        long   ccount  = dsf->conversionCache.lData[0],
               ccount2 = ccount*ccount,
               *ccodes  = dsf->conversionCache.lData+1,
                *tcodes  = dsf->conversionCache.lData+89;

        for (nodeCount = startLeaf; nodeCount<=endLeaf; nodeCount++) {
            long       nMap = dsf->theNodeMap.lData[nodeCount];
            travNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
            char       A = c1[nMap],
                       B = c2[nMap],
                       C = c3[nMap];

            if ((lastIndex<0)||(A!=o1[nMap])||(B!=o2[nMap])||(C!=o3[nMap])) {
                A = ccodes[A-40];
                B = ccodes[B-40];
                C = ccodes[C-40];

                if ((A==-1)||(B==-1)||(C==-1)) {
                    travNode->lastState = dsf->Translate2Frequencies ((*dsf)(index,nodeCount), travNode->theProbs, true);
                } else {
                    travNode->lastState = tcodes[A*ccount2+B*ccount+C];
                    if (travNode->lastState < 0) {
                        travNode->lastState = dsf->Translate2Frequencies ((*dsf)(index,nodeCount), travNode->theProbs, true);
                    } else {
                        for (nMap=0; nMap<travNode->lastState; nMap++) {
                            travNode->theProbs [nMap] = 0.0;
                        }
                        travNode->theProbs [nMap++] = 1.0;
                        for (; nMap<cBase; nMap++) {
                            travNode->theProbs [nMap] = 0.0;
                        }
                    }
                }
                theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                                [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);
                theChildNode->cBase = -1;
            }
        }
    } else {
        if (dsf->GetUnitLength() == 2) {
            char *c1, *c2,
                 *o1, *o2;

            c1 = dsf->GetColumn (2*index);
            c2 = dsf->GetColumn (2*index+1);

            if (lastIndex>=0) {
                o1 = dsf->GetColumn (2*lastIndex);
                o2 = dsf->GetColumn (2*lastIndex+1);
            }

            long   ccount  = dsf->conversionCache.lData[0],
                   *ccodes  = dsf->conversionCache.lData+1,
                    *tcodes  = dsf->conversionCache.lData+89;

            for (nodeCount = startLeaf; nodeCount<=endLeaf; nodeCount++) {
                long       nMap = dsf->theNodeMap.lData[nodeCount];
                travNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
                char       A = c1[nMap],
                           B = c2[nMap];

                if ((lastIndex<0)||(A!=o1[nMap])||(B!=o2[nMap])) {
                    A = ccodes[A-40];
                    B = ccodes[B-40];

                    if ((A==-1)||(B==-1)) {
                        travNode->lastState = dsf->Translate2Frequencies ((*dsf)(index,nodeCount), travNode->theProbs, true);
                    } else {
                        travNode->lastState = tcodes[A*ccount+B];
                        if (travNode->lastState < 0) {
                            travNode->lastState = dsf->Translate2Frequencies ((*dsf)(index,nodeCount), travNode->theProbs, true);
                        } else {
                            for (nMap=0; nMap<travNode->lastState; nMap++) {
                                travNode->theProbs [nMap] = 0.0;
                            }
                            travNode->theProbs [nMap++] = 1;
                            for (; nMap<cBase; nMap++) {
                                travNode->theProbs [nMap] = 0.0;
                            }
                        }
                    }
                    theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                                    [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);
                    theChildNode->cBase = -1;
                }
            }
        } else {
            for (nodeCount = startLeaf; nodeCount<=endLeaf; nodeCount++) {
                travNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
                if ((lastIndex<0)||(!dsf->CompareTwoSites(lastIndex,index,nodeCount))) {
                    travNode->lastState = dsf->Translate2Frequencies ((*dsf)(index,nodeCount), travNode->theProbs, true);
                    theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                                    [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);
                    theChildNode->cBase = -1;
                }
            }
        }
    }

    long f = topLevelNodes.lLength-1,
         mm = topLevelNodes.lData[f],
         startINode = leftiNodes.lData[startLeaf];

    for (long nI = 0; nI < f; nI++,mm>>=1) {
        if (mm&1)
            //if (1)
            // marked for recalculation
        {
            long startAt = 0;

            if (nI) {
                startAt = topLevelNodes.lData[nI-1]+1;
            }

            if (startAt<startINode) {
                startAt = startINode;
            }

            for (nodeCount=startAt; nodeCount<=topLevelNodes.lData[nI]; nodeCount++) {
                theChildNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[nodeCount]);
                // we can immediately place the changes in the parent node
                if (theChildNode->cBase<0) {
                    nodeChild = (node<long>*)flatNodes.lData[nodeCount];
                    for (long i=0; i<cBase; i++) {
                        theChildNode->theProbs[i] = LIKELIHOOD_SCALER;
                    }

                    for (long k=0; k<nodeChild->nodes.length; k++) {
                        travNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->nodes.data[k]->in_object]);
                        fastIndex = travNode->compExp->theData;
                        theProbbs = travNode->theProbs;
                        long      nZ = travNode->lastState;
                        if (nZ>=0) {
                            _Parameter a = theProbbs[nZ];
                            theProbbs = theChildNode->theProbs;
                            fastIndex+=nZ;
                            stopper = theProbbs+cBase;
                            for (; theProbbs!= stopper; fastIndex+=cBase, theProbbs++) {
                                *theProbbs*=a**fastIndex;
                            }

                        } else {
                            /*for (long i=0; i<cBase; i++)
                            {
                                _Parameter tmp = 0.0;
                                for (long j=0; j< cBase; j++)
                                    tmp += theProbbs[j] * fastIndex [j];

                                theChildNode->theProbs[i]*=tmp;
                                fastIndex += cBase;
                            }*/

                            long rem = cBase%4;
                            stopper = theProbbs+(cBase-rem);
                            if (rem==1) {
                                for (long i=0; i<cBase; i++) {
                                    _Parameter tmp = 0.0;
                                    for (; theProbbs!=stopper;) {
                                        tmp += *theProbbs ** fastIndex;
                                        theProbbs++;
                                        fastIndex++;
                                        tmp += *theProbbs ** fastIndex;
                                        theProbbs++;
                                        fastIndex++;
                                        tmp += *theProbbs ** fastIndex;
                                        theProbbs++;
                                        fastIndex++;
                                        tmp += *theProbbs ** fastIndex;
                                        theProbbs++;
                                        fastIndex++;
                                    }
                                    tmp += *(theProbbs)**(fastIndex++);
                                    theChildNode->theProbs[i]*=tmp;
                                    theProbbs = travNode->theProbs;
                                }
                            } else if (rem==2) {
                                for (long i=0; i<cBase; i++) {
                                    _Parameter tmp = 0.0;
                                    for (; theProbbs!=stopper;) {
                                        tmp += *theProbbs ** fastIndex;
                                        theProbbs++;
                                        fastIndex++;
                                        tmp += *theProbbs ** fastIndex;
                                        theProbbs++;
                                        fastIndex++;
                                        tmp += *theProbbs ** fastIndex;
                                        theProbbs++;
                                        fastIndex++;
                                        tmp += *theProbbs ** fastIndex;
                                        theProbbs++;
                                        fastIndex++;
                                    }
                                    tmp+= *(theProbbs++)**(fastIndex++);
                                    tmp+= *(theProbbs)**(fastIndex++);
                                    theChildNode->theProbs[i]*=tmp;
                                    theProbbs = travNode->theProbs;
                                }

                            } else if (rem==3) {
                                for (long i=0; i<cBase; i++) {
                                    _Parameter tmp = 0.0;
                                    for (; theProbbs!=stopper;) {
                                        tmp += *theProbbs ** fastIndex;
                                        theProbbs++;
                                        fastIndex++;
                                        tmp += *theProbbs ** fastIndex;
                                        theProbbs++;
                                        fastIndex++;
                                        tmp += *theProbbs ** fastIndex;
                                        theProbbs++;
                                        fastIndex++;
                                        tmp += *theProbbs ** fastIndex;
                                        theProbbs++;
                                        fastIndex++;
                                    }
                                    tmp+= *(theProbbs++)**(fastIndex++);
                                    tmp+= *(theProbbs++)**(fastIndex++);
                                    tmp+=*(theProbbs)**(fastIndex++);
                                    theChildNode->theProbs[i]*=tmp;
                                    theProbbs = travNode->theProbs;
                                }

                            } else {
                                for (long i=0; i<cBase; i++) {
                                    _Parameter tmp = 0.0;
                                    for (; theProbbs!=stopper;) {
                                        tmp += *theProbbs ** fastIndex;
                                        theProbbs++;
                                        fastIndex++;
                                        tmp += *theProbbs ** fastIndex;
                                        theProbbs++;
                                        fastIndex++;
                                        tmp += *theProbbs ** fastIndex;
                                        theProbbs++;
                                        fastIndex++;
                                        tmp += *theProbbs ** fastIndex;
                                        theProbbs++;
                                        fastIndex++;
                                    }
                                    theChildNode->theProbs[i]*=tmp;
                                    theProbbs = travNode->theProbs;
                                }
                            }
                        }
                    }
                    theChildNode->cBase = cBase;
                    nodeChild = nodeChild->parent;
                    if (nodeChild) {
                        ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->in_object])->cBase = -1;
                    }
                }
            }

            theChildNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[topLevelNodes.lData[nI]]);
            _Parameter* nodeProbs = theChildNode->theProbs;
            fastIndex = rootIChildrenCache+position*cBase*(topLevelNodes.lLength-1)+cBase*nI;
            for (long i=0; i<cBase; i++) {
                *(fastIndex++) = *(nodeProbs++);
            }
        } else {
            _Parameter* nodeProbs = ((_CalcNode*)(((BaseRef*)flatTree.lData)[topLevelNodes.lData[nI]]))->theProbs;
            fastIndex = rootIChildrenCache+position*cBase*(topLevelNodes.lLength-1)+cBase*nI;
            for (long i=0; i<cBase; i++) {
                *(nodeProbs++) = *(fastIndex++);
            }

            nodeChild = ((node <long>*)(flatNodes.lData[topLevelNodes.lData[nI]]))->parent;
            ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->in_object])->cBase = -1;
        }
    }

    for (nodeCount=topLevelNodes.lData[f-1]+1; nodeCount<flatTree.lLength; nodeCount++) {
        theChildNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[nodeCount]);
        if (theChildNode->cBase<0) {
            nodeChild = (node<long>*)flatNodes.lData[nodeCount];
            for (long i=0; i<cBase; i++) {
                theChildNode->theProbs[i] = LIKELIHOOD_SCALER;
            }

            for (long k=0; k<nodeChild->nodes.length; k++) {
                travNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->nodes.data[k]->in_object]);
                fastIndex = travNode->compExp->theData;
                theProbbs = travNode->theProbs;
                long      nZ = travNode->lastState;

                if (nZ>=0) {
                    _Parameter a = theProbbs[nZ];
                    theProbbs = theChildNode->theProbs;
                    fastIndex+=nZ;
                    stopper = theProbbs+cBase;
                    for (; theProbbs!=stopper; theProbbs++, fastIndex+=cBase) {
                        *theProbbs*=a**fastIndex;
                    }
                } else {
                    /*for (long i=0; i<cBase; i++)
                    {
                        _Parameter tmp = 0.0;
                        for (long j=0;j<cBase;j++)
                            tmp += theProbbs[j] *fastIndex[j];

                        theChildNode->theProbs[i]*=tmp;
                        fastIndex += cBase;
                    }*/

                    long rem = cBase%4;
                    stopper = theProbbs+(cBase-rem);
                    if (rem==1) {
                        for (long i=0; i<cBase; i++) {
                            _Parameter tmp = 0.0;
                            for (; theProbbs!=stopper;) {
                                tmp += *theProbbs ** fastIndex;
                                theProbbs++;
                                fastIndex++;
                                tmp += *theProbbs ** fastIndex;
                                theProbbs++;
                                fastIndex++;
                                tmp += *theProbbs ** fastIndex;
                                theProbbs++;
                                fastIndex++;
                                tmp += *theProbbs ** fastIndex;
                                theProbbs++;
                                fastIndex++;
                            }
                            tmp += *(theProbbs)**(fastIndex++);
                            theChildNode->theProbs[i]*=tmp;
                            theProbbs = travNode->theProbs;
                        }
                    } else if (rem==2) {
                        for (long i=0; i<cBase; i++) {
                            _Parameter tmp = 0.0;
                            for (; theProbbs!=stopper;) {
                                tmp += *theProbbs ** fastIndex;
                                theProbbs++;
                                fastIndex++;
                                tmp += *theProbbs ** fastIndex;
                                theProbbs++;
                                fastIndex++;
                                tmp += *theProbbs ** fastIndex;
                                theProbbs++;
                                fastIndex++;
                                tmp += *theProbbs ** fastIndex;
                                theProbbs++;
                                fastIndex++;
                            }
                            tmp+= *(theProbbs++)**(fastIndex++);
                            tmp+= *(theProbbs)**(fastIndex++);
                            theChildNode->theProbs[i]*=tmp;
                            theProbbs = travNode->theProbs;
                        }

                    } else if (rem==3) {
                        for (long i=0; i<cBase; i++) {
                            _Parameter tmp = 0.0;
                            for (; theProbbs!=stopper;) {
                                tmp += *theProbbs ** fastIndex;
                                theProbbs++;
                                fastIndex++;
                                tmp += *theProbbs ** fastIndex;
                                theProbbs++;
                                fastIndex++;
                                tmp += *theProbbs ** fastIndex;
                                theProbbs++;
                                fastIndex++;
                                tmp += *theProbbs ** fastIndex;
                                theProbbs++;
                                fastIndex++;
                            }
                            tmp+= *(theProbbs++)**(fastIndex++);
                            tmp+= *(theProbbs++)**(fastIndex++);
                            tmp+=*(theProbbs)**(fastIndex++);
                            theChildNode->theProbs[i]*=tmp;
                            theProbbs = travNode->theProbs;
                        }

                    } else {
                        for (long i=0; i<cBase; i++) {
                            _Parameter tmp = 0.0;
                            for (; theProbbs!=stopper;) {
                                tmp += *theProbbs ** fastIndex;
                                theProbbs++;
                                fastIndex++;
                                tmp += *theProbbs ** fastIndex;
                                theProbbs++;
                                fastIndex++;
                                tmp += *theProbbs ** fastIndex;
                                theProbbs++;
                                fastIndex++;
                                tmp += *theProbbs ** fastIndex;
                                theProbbs++;
                                fastIndex++;
                            }
                            theChildNode->theProbs[i]*=tmp;
                            theProbbs = travNode->theProbs;
                        }
                    }
                }
            }
            theChildNode->cBase = cBase;
            nodeChild = nodeChild->parent;
            if (nodeChild) {
                ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->in_object])->cBase = -1;
            }
        }
    }

    _Parameter result = 0.;

    for (long i=0; i<cBase; i++) {
        result+= theProbs[i]*theChildNode->theProbs[i];
    }

    if (result<=0.0) {
        return ALMOST_ZERO;
    }

    return result;
}

//_______________________________________________________________________________________________

#if USE_SCALING_TO_FIX_UNDERFLOW
_Parameter   _TheTree::ThreadReleafTreeCache (_DataSetFilter* dsf, long index, long lastIndex, long startLeaf, long endLeaf, long position, long offset, long fixAttempt, _Parameter localScalingFactor)
#else
_Parameter   _TheTree::ThreadReleafTreeCache (_DataSetFilter* dsf, long index, long lastIndex, long startLeaf, long endLeaf, long position, long offset)

#endif
// set leaf value and re-prune w/o reexping
// for subsequent entries into a datafilter
{

    long        nodeCount,
                * ns  = nodeStates+(flatLeaves.lLength+flatNodes.lLength)*offset;

    _CalcNode * travNode,
              * theChildNode;

    _Parameter* fastIndex,
                * theProbbs,
                * stopper,
                * mlc = marginalLikelihoodCache+cBase*(flatLeaves.lLength+flatNodes.lLength)*offset;

    node<long>* nodeChild;

    char      * nm  = nodeMarkers+(flatNodes.lLength)*offset,
                unitSize = dsf->GetUnitLength ();

#if USE_SCALING_TO_FIX_UNDERFLOW
    _Parameter scalingFactor = scalingForUnderflow->theData[index];
    if (lastIndex>=0 && scalingFactor != scalingForUnderflow->theData[lastIndex]) {
        lastIndex = -1;
        startLeaf = 0;
        endLeaf = flatLeaves.lLength-1;
    }
    bool overflowFlag = fixAttempt<0;
    if (overflowFlag) {
        fixAttempt = -fixAttempt;
    }

#endif

    if (unitSize == 3) {
        char *c1, *c2, *c3,
             *o1, *o2, *o3;

        c1 = dsf->GetColumn (3*index);
        c2 = dsf->GetColumn (3*index+1);
        c3 = dsf->GetColumn (3*index+2);

        if (lastIndex>=0) {
            o1 = dsf->GetColumn (3*lastIndex);
            o2 = dsf->GetColumn (3*lastIndex+1);
            o3 = dsf->GetColumn (3*lastIndex+2);
        }

        long   ccount  = dsf->conversionCache.lData[0],
               ccount2 = ccount*ccount,
               *ccodes  = dsf->conversionCache.lData+1,
                *tcodes  = dsf->conversionCache.lData+89;

        for (nodeCount = startLeaf; nodeCount<=endLeaf; nodeCount++) {
            long       nMap = dsf->theNodeMap.lData[nodeCount];
            char       A = c1[nMap],
                       B = c2[nMap],
                       C = c3[nMap];

            if ( lastIndex<0 || A!=o1[nMap] || B!=o2[nMap] || C!=o3[nMap]) {
                A = ccodes[A-40];
                B = ccodes[B-40];
                C = ccodes[C-40];
                fastIndex = mlc+cBase*nodeCount;

                if (A==-1 || B==-1 || C==-1) {
                    _String codon (3,false);
                    dsf->GrabSite (index,nodeCount,codon);
                    ns[nodeCount] = dsf->Translate2Frequencies (codon, fastIndex, true);
#if USE_SCALING_TO_FIX_UNDERFLOW
                    if (scalingFactor != 1.0)
                        for (nMap = 0; nMap < cBase; nMap ++) {
                            fastIndex[nMap] *= scalingFactor;
                        }
#endif
                } else {
                    long nState = ns[nodeCount] = tcodes[A*ccount2+B*ccount+C];
                    if (nState < 0) {
                        _String codon (3,false);
                        dsf->GrabSite (index,nodeCount,codon);
                        ns[nodeCount] = dsf->Translate2Frequencies (codon,fastIndex, true);
#if USE_SCALING_TO_FIX_UNDERFLOW
                        if (scalingFactor != 1.0)
                            for (nMap = 0; nMap < cBase; nMap ++) {
                                fastIndex[nMap] *= scalingFactor;
                            }
#endif
                    } else {
                        for (long z = 0; z< cBase; z++) {
                            fastIndex[z] = 0.0;
                        }
#if USE_SCALING_TO_FIX_UNDERFLOW
                        fastIndex [nState] = scalingFactor;
#else
                        fastIndex [nState] = 1.0;
#endif
                    }
                }
                theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                                [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);
                nm[theChildNode->nodeIndex-flatLeaves.lLength] = -1;
            }
        }
    } else {
        if (unitSize == 2) {
            char *c1, *c2,
                 *o1, *o2;

            c1 = dsf->GetColumn (2*index);
            c2 = dsf->GetColumn (2*index+1);

            if (lastIndex>=0) {
                o1 = dsf->GetColumn (2*lastIndex);
                o2 = dsf->GetColumn (2*lastIndex+1);
            }

            long   ccount  = dsf->conversionCache.lData[0],
                   *ccodes  = dsf->conversionCache.lData+1,
                    *tcodes  = dsf->conversionCache.lData+89;

            for (nodeCount = startLeaf; nodeCount<=endLeaf; nodeCount++) {
                long       nMap = dsf->theNodeMap.lData[nodeCount];
                //travNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
                char       A = c1[nMap],
                           B = c2[nMap];

                if ((lastIndex<0)||(A!=o1[nMap])||(B!=o2[nMap])) {
                    A = ccodes[A-40];
                    B = ccodes[B-40];

                    fastIndex = mlc+cBase*nodeCount;

                    if ((A==-1)||(B==-1)) {
                        _String di (2,false);
                        dsf->GrabSite (index,nodeCount,di);
                        ns[nodeCount] = dsf->Translate2Frequencies (di, fastIndex, true);
                    } else {
                        long nState = ns[nodeCount] = tcodes[A*ccount+B];

                        if (nState < 0) {
                            _String di (2,false);
                            dsf->GrabSite (index,nodeCount,di);
                            ns[nodeCount] = dsf->Translate2Frequencies (di, fastIndex, true);
                        } else {
                            for (nMap=0; nMap<nState; nMap++) {
                                fastIndex [nMap] = 0.0;
                            }
                            fastIndex [nMap++] = 1.;
                            for (; nMap<cBase; nMap++) {
                                fastIndex [nMap] = 0.0;
                            }
                        }
                    }
                    theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                                    [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);
                    nm[theChildNode->nodeIndex-flatLeaves.lLength] = -1;
                }
            }
        } else {
            for (nodeCount = startLeaf; nodeCount<=endLeaf; nodeCount++) {
                //travNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
                if ((lastIndex<0)||(!dsf->CompareTwoSites(lastIndex,index,nodeCount))) {
                    _String         dPoint (unitSize,false);
                    dsf->GrabSite (index,nodeCount,dPoint);
                    ns[nodeCount] = dsf->Translate2Frequencies (dPoint, mlc+cBase*nodeCount, true);
                    theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                                    [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);
                    nm[theChildNode->nodeIndex-flatLeaves.lLength] = -1;
                }
            }
        }
    }

    long f          = topLevelNodes.lLength-1,
         mm       = topLevelNodes.lData[f],
         startINode = leftiNodes.lData[startLeaf],
         rem       = cBase%4,
         upTo        = cBase-rem;


    for (long nI = 0; nI < f; nI++,mm>>=1) {
        if (mm&1)
            //if (1)
            // marked for recalculation
        {
            long startAt = 0;

            if (nI) {
                startAt = topLevelNodes.lData[nI-1]+1;
            }

            if (startAt<startINode) {
                startAt = startINode;
            }

            for (nodeCount=startAt; nodeCount<=topLevelNodes.lData[nI]; nodeCount++) {
                if (nm[nodeCount]<0) {
                    nodeChild = (node<long>*)flatNodes.lData[nodeCount];
                    _Parameter* theseProbs = mlc + cBase*(nodeCount+flatLeaves.lLength);
                    for (long i=0; i<cBase; i++) {
                        theseProbs[i] = LIKELIHOOD_SCALER;
                    }

                    for (long k=0; k<nodeChild->nodes.length; k++) {
                        travNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->nodes.data[k]->in_object]);
                        fastIndex = travNode->compExp->theData;
                        _Parameter * cProbs =  mlc + cBase*(travNode->nodeIndex);
                        theProbbs = cProbs;
                        long      nZ = ns[travNode->nodeIndex];
                        if (nZ>=0) {
                            _Parameter a = theProbbs[nZ];
                            theProbbs       = theseProbs;
                            fastIndex       +=nZ;
                            stopper          = theProbbs+cBase;
                            for (; theProbbs!= stopper; fastIndex+=cBase, theProbbs++)
#if USE_SCALING_TO_FIX_UNDERFLOW
                                *theProbbs*=a**fastIndex;// * scalingFactor;
#else
                                *theProbbs*=a**fastIndex;
#endif
                        } else {
                            stopper = theProbbs+(cBase-rem);
                            if (rem==1) {
                                for (long i=0; i<cBase; i++) {
                                    _Parameter tmp1  = 0.0,
                                               tmp2  = 0.0,
                                               tmp3  = 0.0,
                                               tmp4  = 0.0;

                                    for (long k=0; k<upTo; k+=4) {
                                        tmp1 += cProbs[k]  *fastIndex[k];
                                        tmp2 += cProbs[k+1]*fastIndex[k+1];
                                        tmp3 += cProbs[k+2]*fastIndex[k+2];
                                        tmp4 += cProbs[k+3]*fastIndex[k+3];
                                    }

                                    theseProbs[i] *= tmp1 + cProbs[upTo]*fastIndex[upTo] + tmp2+tmp3+tmp4;
                                    fastIndex += cBase;
                                }
                            } else if (rem==2) {
                                for (long i=0; i<cBase; i++) {
                                    _Parameter tmp1  = 0.0,
                                               tmp2  = 0.0,
                                               tmp3  = 0.0,
                                               tmp4  = 0.0;

                                    for (long k=0; k<upTo; k+=4) {
                                        tmp1 += cProbs[k]  *fastIndex[k];
                                        tmp2 += cProbs[k+1]*fastIndex[k+1];
                                        tmp3 += cProbs[k+2]*fastIndex[k+2];
                                        tmp4 += cProbs[k+3]*fastIndex[k+3];
                                    }

                                    theseProbs[i] *= tmp1+tmp2+tmp3+tmp4+
                                                     + cProbs[upTo]*fastIndex[upTo] +
                                                     + cProbs[upTo+1]*fastIndex[upTo+1];
                                    fastIndex += cBase;
                                }
                            } else if (rem==3)
                                for (long i=0; i<cBase; i++) {
                                    _Parameter tmp1  = 0.0,
                                               tmp2  = 0.0,
                                               tmp3  = 0.0,
                                               tmp4  = 0.0;

                                    for (long k=0; k<upTo; k+=4) {
                                        tmp1 += cProbs[k]  *fastIndex[k];
                                        tmp2 += cProbs[k+1]*fastIndex[k+1];
                                        tmp3 += cProbs[k+2]*fastIndex[k+2];
                                        tmp4 += cProbs[k+3]*fastIndex[k+3];
                                    }

                                    theseProbs[i] *= tmp1+tmp2+tmp3+tmp4+
                                                     + cProbs[upTo]*fastIndex[upTo] +
                                                     + cProbs[upTo+1]*fastIndex[upTo+1]+
                                                     + cProbs[upTo+2]*fastIndex[upTo+2];

                                    fastIndex += cBase;
                                }
                            else
                                for (long i=0; i<cBase; i++) {
                                    _Parameter tmp1  = 0.0,
                                               tmp2  = 0.0,
                                               tmp3  = 0.0,
                                               tmp4  = 0.0;

                                    for (long k=0; k<upTo; k+=4) {
                                        tmp1 += cProbs[k]  *fastIndex[k];
                                        tmp2 += cProbs[k+1]*fastIndex[k+1];
                                        tmp3 += cProbs[k+2]*fastIndex[k+2];
                                        tmp4 += cProbs[k+3]*fastIndex[k+3];
                                    }

                                    theseProbs[i] *= tmp1+tmp2+tmp3+tmp4;
                                    fastIndex     += cBase;
                                }
                        }
#if USE_SCALING_TO_FIX_UNDERFLOW
                        long nodeCheck = 0;
                        for (nodeCheck = 0; nodeCheck < cBase; nodeCheck++)
                            if (theseProbs[nodeCheck] > 0.0 || theseProbs[nodeCheck] == HUGE_VAL) {
                                break;
                            }

                        if (nodeCheck == cBase || theseProbs[nodeCheck] == HUGE_VAL) {
                            return doScaling (dsf, index, position, offset,fixAttempt,localScalingFactor,overflowFlag,nodeCheck < cBase);
                        }
#endif
                    }
                    nm[nodeCount] = 0;
                    nodeChild = nodeChild->parent;
                    if (nodeChild) {
                        nm[((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->in_object])->nodeIndex-flatLeaves.lLength] = -1;
                    }
                }
            }

            theChildNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[topLevelNodes.lData[nI]]);
            _Parameter* nodeProbs = mlc+cBase*(theChildNode->nodeIndex);
            fastIndex = rootIChildrenCache+position*cBase*(topLevelNodes.lLength-1)+cBase*nI;
            for (long i=0; i<cBase; i++) {
                *(fastIndex++) = *(nodeProbs++);
            }
        } else {
            _Parameter* nodeProbs = mlc+cBase*
                                    (((_CalcNode*)(((BaseRef*)flatTree.lData)[topLevelNodes.lData[nI]]))->nodeIndex);

            fastIndex = rootIChildrenCache+position*cBase*(topLevelNodes.lLength-1)+cBase*nI;
            for (long i=0; i<cBase; i++) {
                *(nodeProbs++) = *(fastIndex++);
            }

            nodeChild = ((node <long>*)(flatNodes.lData[topLevelNodes.lData[nI]]))->parent;
            nm[((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->in_object])->nodeIndex-flatLeaves.lLength] = -1;
        }
    }

    for (nodeCount=topLevelNodes.lData[f-1]+1; nodeCount<flatTree.lLength; nodeCount++) {
        if (nm[nodeCount]<0) {
            nodeChild = (node<long>*)flatNodes.lData[nodeCount];
            _Parameter* theseProbs = mlc + cBase*(nodeCount+flatLeaves.lLength);
            for (long i=0; i<cBase; i++) {
                theseProbs[i] = LIKELIHOOD_SCALER;
            }

            for (long k=0; k<nodeChild->nodes.length; k++) {
                travNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->nodes.data[k]->in_object]);
                fastIndex = travNode->compExp->theData;
                _Parameter * cProbs =  mlc + cBase*(travNode->nodeIndex);
                theProbbs = cProbs;
                long      nZ = ns[travNode->nodeIndex];
                if (nZ>=0) {
                    _Parameter a = theProbbs[nZ];
                    theProbbs = theseProbs;
                    fastIndex+=nZ;
                    stopper = theProbbs+cBase;
                    for (; theProbbs!= stopper; fastIndex+=cBase, theProbbs++)
#if USE_SCALING_TO_FIX_UNDERFLOW
                        *theProbbs*=a**fastIndex;// * scalingFactor;
#else
                        *theProbbs*=a**fastIndex;
#endif
                } else {
                    stopper = theProbbs+(cBase-rem);
                    if (rem==1) {
                        for (long i=0; i<cBase; i++) {
                            _Parameter tmp1  = 0.0,
                                       tmp2  = 0.0,
                                       tmp3  = 0.0,
                                       tmp4  = 0.0;

                            for (long k=0; k<upTo; k+=4) {
                                tmp1 += cProbs[k]  *fastIndex[k];
                                tmp2 += cProbs[k+1]*fastIndex[k+1];
                                tmp3 += cProbs[k+2]*fastIndex[k+2];
                                tmp4 += cProbs[k+3]*fastIndex[k+3];
                            }

                            theseProbs[i] *= tmp1 + cProbs[upTo]*fastIndex[upTo] + tmp2+tmp3+tmp4;
                            fastIndex += cBase;
                        }
                    } else if (rem==2) {
                        for (long i=0; i<cBase; i++) {
                            _Parameter tmp1  = 0.0,
                                       tmp2  = 0.0,
                                       tmp3  = 0.0,
                                       tmp4  = 0.0;

                            for (long k=0; k<upTo; k+=4) {
                                tmp1 += cProbs[k]  *fastIndex[k];
                                tmp2 += cProbs[k+1]*fastIndex[k+1];
                                tmp3 += cProbs[k+2]*fastIndex[k+2];
                                tmp4 += cProbs[k+3]*fastIndex[k+3];
                            }

                            theseProbs[i] *= tmp1 + cProbs [upTo]   * fastIndex [upTo]
                                             + cProbs [upTo+1] * fastIndex [upTo+1]
                                             + tmp2+tmp3+tmp4;
                            fastIndex += cBase;
                        }

                    } else if (rem==3)
                        for (long i=0; i<cBase; i++) {
                            _Parameter tmp1  = 0.0,
                                       tmp2  = 0.0,
                                       tmp3  = 0.0,
                                       tmp4  = 0.0;

                            for (long k=0; k<upTo; k+=4) {
                                tmp1 += cProbs[k]  *fastIndex[k];
                                tmp2 += cProbs[k+1]*fastIndex[k+1];
                                tmp3 += cProbs[k+2]*fastIndex[k+2];
                                tmp4 += cProbs[k+3]*fastIndex[k+3];
                            }

                            theseProbs[i] *= tmp1 + cProbs [upTo]   * fastIndex [upTo]
                                             + cProbs [upTo+1] * fastIndex [upTo+1]
                                             + cProbs [upTo+2] * fastIndex [upTo+2]
                                             + tmp2+tmp3+tmp4;
                            fastIndex += cBase;
                        }
                    else {
                        for (long i=0; i<cBase; i++) {
                            _Parameter tmp1  = 0.0,
                                       tmp2  = 0.0,
                                       tmp3  = 0.0,
                                       tmp4  = 0.0;

                            for (long k=0; k<upTo; k+=4) {
                                tmp1 += cProbs[k]  *fastIndex[k];
                                tmp2 += cProbs[k+1]*fastIndex[k+1];
                                tmp3 += cProbs[k+2]*fastIndex[k+2];
                                tmp4 += cProbs[k+3]*fastIndex[k+3];
                            }

                            theseProbs[i] *= tmp1 + tmp2+tmp3+tmp4;
                            fastIndex += cBase;
                        }
                    }
                }
            }
            nm[nodeCount] = 0;
            nodeChild = nodeChild->parent;
            if (nodeChild) {
                nm[((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->in_object])->nodeIndex-flatLeaves.lLength] = -1;
            }

#if USE_SCALING_TO_FIX_UNDERFLOW
            long nodeCheck = 0;
            for (nodeCheck = 0; nodeCheck < cBase; nodeCheck++)
                if (theseProbs[nodeCheck] > 0.0 || theseProbs[nodeCheck] == HUGE_VAL) {
                    break;
                }

            if (nodeCheck == cBase || theseProbs[nodeCheck] == HUGE_VAL) {
                return doScaling (dsf, index, position, offset,fixAttempt,localScalingFactor,overflowFlag,nodeCheck < cBase);
            }
#endif

        }

    }

    _Parameter result = 0;

    fastIndex = mlc+cBase*(flatNodes.lLength+flatLeaves.lLength-1);

    for (long i=0; i<cBase; i++) {
        result+= theProbs[i]*fastIndex[i];
    }

#if USE_SCALING_TO_FIX_UNDERFLOW
    if (result <= 0.0) { // underflow
        return doScaling (dsf, index, position, offset,fixAttempt,localScalingFactor,overflowFlag,false);
    } else if (result == HUGE_VAL) {
        return doScaling (dsf, index, position, offset,fixAttempt,localScalingFactor,overflowFlag,true);
    }
#else
    if (result<=0.0) {
        return ALMOST_ZERO;
    }
#endif

    return result;
}

//_______________________________________________________________________________________________

_Parameter   _TheTree::ReleafTreeChar4 (_DataSetFilter* dsf, long index, long lastIndex, long startLeaf, long endLeaf, long position)

// set leaf value and re-prune w/o reexping
// for subsequent entries into a datafilter
{

    _CalcNode               * travNode,
                            * theChildNode;

    long                    nodeCount,
                            f;

    _Parameter              * fastIndex;

    node<long>              * nodeChild;

    char                    *pastState = dsf->GetColumn(lastIndex),
                             *thisState = dsf->GetColumn(index);

    for (nodeCount = startLeaf; nodeCount<=endLeaf; nodeCount++) {
        f = dsf->theNodeMap.lData[nodeCount];
        if (thisState[f]!=pastState[f]) {
            travNode     = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
            long* cCache = dsf->conversionCache.lData+(thisState[f]-40)*5;
            travNode->theProbs[0] = *(cCache++);
            travNode->theProbs[1] = *(cCache++);
            travNode->theProbs[2] = *(cCache++);
            travNode->theProbs[3] = *(cCache++);
            travNode->lastState = *cCache;
            theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                            [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);
            if (theChildNode->cBase>0) {
                theChildNode->cBase = -1;
            }
        }
    }

    f = topLevelNodes.lLength-1;
    long mm = topLevelNodes.lData[f],
         startINode = leftiNodes.lData[startLeaf];

    for (long nI = 0; nI < f; nI++,mm>>=1) {
        if (mm&1)
            //if (1)
            // marked for recalculation
        {
            long startAt = 0;

            if (nI) {
                startAt = topLevelNodes.lData[nI-1]+1;
            }

            if (startAt<startINode) {
                startAt = startINode;
            }

            for (nodeCount=startAt; nodeCount<=topLevelNodes.lData[nI]; nodeCount++) {
                theChildNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[nodeCount]);
                if (theChildNode->cBase == -1) { // marked for update
                    nodeChild = ((node <long>*)(flatNodes.lData[nodeCount]));
                    theChildNode->theProbs[0]=LIKELIHOOD_SCALER;
                    theChildNode->theProbs[1]=LIKELIHOOD_SCALER;
                    theChildNode->theProbs[2]=LIKELIHOOD_SCALER;
                    theChildNode->theProbs[3]=LIKELIHOOD_SCALER;
                    theChildNode->cBase = 4;

                    for (long k=0; k <nodeChild->nodes.length; k++) {
                        travNode =
                            ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->nodes.data[k]->in_object]);
                        fastIndex = travNode->compExp->theData;
                        if (travNode->lastState<0) {
                            _Parameter a = *(travNode->theProbs),
                                       b = travNode->theProbs[1],
                                       c = travNode->theProbs[2],
                                       d = travNode->theProbs[3];

                            _Parameter tmp = a * *(fastIndex++);
                            tmp +=b* *(fastIndex++);
                            tmp +=c* *(fastIndex++);
                            theChildNode->theProbs[0]*=tmp+d* *(fastIndex++);

                            tmp = a * *(fastIndex++);
                            tmp +=b* *(fastIndex++);
                            tmp +=c* *(fastIndex++);
                            theChildNode->theProbs[1]*=tmp+d* *(fastIndex++);

                            tmp = a * *(fastIndex++);
                            tmp +=b* *(fastIndex++);
                            tmp +=c* *(fastIndex++);
                            theChildNode->theProbs[2]*=tmp+d* *(fastIndex++);

                            tmp = a * *(fastIndex++);
                            tmp +=b* *(fastIndex++);
                            tmp +=c* *(fastIndex++);
                            theChildNode->theProbs[3]*=tmp+d* *(fastIndex++);
                        } else {
                            fastIndex+=travNode->lastState;
                            theChildNode->theProbs[0]*=*fastIndex;
                            theChildNode->theProbs[1]*=fastIndex[4];
                            theChildNode->theProbs[2]*=fastIndex[8];
                            theChildNode->theProbs[3]*=fastIndex[12];
                        }
                    }
                    ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->parent->in_object])->cBase = -1;

                }
            }
            theChildNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[topLevelNodes.lData[nI]]);
            _Parameter* nodeProbs = theChildNode->theProbs;
            fastIndex = rootIChildrenCache+position*4*(topLevelNodes.lLength-1)+4*nI;
            fastIndex[0] = nodeProbs[0];
            fastIndex[1] = nodeProbs[1];
            fastIndex[2] = nodeProbs[2];
            fastIndex[3] = nodeProbs[3];
        } else {
            _Parameter* nodeProbs = ((_CalcNode*)(((BaseRef*)flatTree.lData)[topLevelNodes.lData[nI]]))->theProbs;
            fastIndex = rootIChildrenCache+position*4*(topLevelNodes.lLength-1)+4*nI;
            nodeProbs[0] = fastIndex[0];
            nodeProbs[1] = fastIndex[1];
            nodeProbs[2] = fastIndex[2];
            nodeProbs[3] = fastIndex[3];
            nodeChild = ((node <long>*)(flatNodes.lData[topLevelNodes.lData[nI]]))->parent;
            ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->in_object])->cBase = -1;
        }
    }

    for (nodeCount=topLevelNodes.lData[f-1]+1; nodeCount<flatTree.lLength; nodeCount++)
        //for (nodeCount=startINode; nodeCount<flatTree.lLength; nodeCount++)
    {
        theChildNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[nodeCount]);
        if (theChildNode->cBase == -1) { // marked for update
            nodeChild = ((node <long>*)(flatNodes.lData[nodeCount]));
            theChildNode->theProbs[0]=LIKELIHOOD_SCALER;
            theChildNode->theProbs[1]=LIKELIHOOD_SCALER;
            theChildNode->theProbs[2]=LIKELIHOOD_SCALER;
            theChildNode->theProbs[3]=LIKELIHOOD_SCALER;
            theChildNode->cBase = 4;

            for (long k=0; k <nodeChild->nodes.length; k++) {
                travNode =
                    ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->nodes.data[k]->in_object]);
                fastIndex = travNode->compExp->theData;
                if (travNode->lastState<0) {
                    _Parameter a = *(travNode->theProbs),
                               b = travNode->theProbs[1],
                               c = travNode->theProbs[2],
                               d = travNode->theProbs[3];

                    _Parameter tmp = a * *(fastIndex++);
                    tmp +=b* *(fastIndex++);
                    tmp +=c* *(fastIndex++);
                    theChildNode->theProbs[0]*=tmp+d* *(fastIndex++);

                    tmp = a * *(fastIndex++);
                    tmp +=b* *(fastIndex++);
                    tmp +=c* *(fastIndex++);
                    theChildNode->theProbs[1]*=tmp+d* *(fastIndex++);

                    tmp = a * *(fastIndex++);
                    tmp +=b* *(fastIndex++);
                    tmp +=c* *(fastIndex++);
                    theChildNode->theProbs[2]*=tmp+d* *(fastIndex++);

                    tmp = a * *(fastIndex++);
                    tmp +=b* *(fastIndex++);
                    tmp +=c* *(fastIndex++);
                    theChildNode->theProbs[3]*=tmp+d* *(fastIndex++);
                } else {
                    fastIndex+=travNode->lastState;
                    theChildNode->theProbs[0]*=*fastIndex;
                    theChildNode->theProbs[1]*=fastIndex[4];
                    theChildNode->theProbs[2]*=fastIndex[8];
                    theChildNode->theProbs[3]*=fastIndex[12];
                }
            }
            nodeChild = nodeChild->parent;
            if (nodeChild) {
                ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->in_object])->cBase = -1;
            }
        }
    }

    theChildNode->cBase = 4;

    _Parameter  result = theProbs[0]*theChildNode->theProbs[0]+
                         theProbs[1]*theChildNode->theProbs[1]+
                         theProbs[2]*theChildNode->theProbs[2]+
                         theProbs[3]*theChildNode->theProbs[3];

    if (result<=0.0) {
        return ALMOST_ZERO;
    }

    return result;
}

//_______________________________________________________________________________________________
//_______________________________________________________________________________________________

#if USE_SCALING_TO_FIX_UNDERFLOW
_Parameter   _TheTree::ThreadReleafTreeChar4
(_DataSetFilter* dsf, long index, long lastIndex, long startLeaf, long endLeaf, long position, long offset, long fixAttempt, _Parameter localScalingFactor)
#else
_Parameter   _TheTree::ThreadReleafTreeChar4
(_DataSetFilter* dsf, long index, long lastIndex, long startLeaf, long endLeaf, long position, long offset)
#endif
{

    _CalcNode               * travNode,
                            * theChildNode;

    long                    nodeCount,
                            f = topLevelNodes.lLength-1,
                            * ns  = nodeStates+(flatLeaves.lLength+flatNodes.lLength)*offset;

    _Parameter              * fastIndex,
                            * mlc = marginalLikelihoodCache+cBase*(flatLeaves.lLength+flatNodes.lLength)*offset;

    node<long>              * nodeChild;

    char                    *pastState ,
                            *thisState      =   dsf->GetColumn(index),
                             *nm             =   nodeMarkers+(flatNodes.lLength)*offset;

    if (lastIndex>=0) {
        pastState = dsf->GetColumn(lastIndex);
    }

#if USE_SCALING_TO_FIX_UNDERFLOW
    if (lastIndex>=0 && scalingForUnderflow->theData[index] != scalingForUnderflow->theData[lastIndex]) {
        lastIndex = -1;
        startLeaf = 0;
        endLeaf   = flatLeaves.lLength-1;
    }
    _Parameter scalingFactor = scalingForUnderflow->theData[index];
    bool            overflowFlag = fixAttempt<0;
    if (overflowFlag) {
        fixAttempt = -fixAttempt;
    }
#endif

    long mm         = topLevelNodes.lData[f],
         startINode = leftiNodes.lData[startLeaf],
         fromL        = startLeaf,
         toL        = 0;

    // mod SLKP 20070209 to skip over constant parts of the tree

    for (long tlc = 0; tlc < topLevelNodes.lLength; tlc++,mm>>=1) {
        if (mm&1) {
            toL = topLevelRightL.lData[tlc];
        } else {
            toL = topLevelLeftL.lData[tlc]-1;
        }

        if (fromL<=toL) {
            long l1 = MAX(startLeaf,fromL),
                 l2 = MIN(endLeaf,  toL);

            for (nodeCount = l1; nodeCount<=l2; nodeCount++) {
                long f2 = dsf->theNodeMap.lData[nodeCount];
                if ( lastIndex==-1 || thisState[f2]!=pastState[f2]
#if USE_SCALING_TO_FIX_UNDERFLOW
                        || scalingFactor != scalingForUnderflow->theData[lastIndex]
#endif
                   ) {
                    travNode     = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
                    long* cCache = dsf->conversionCache.lData+(thisState[f2]-40)*5;

                    _Parameter* fi = mlc+nodeCount*4;

                    fi[0] = cCache[0];
                    fi[1] = cCache[1];
                    fi[2] = cCache[2];
                    fi[3] = cCache[3];

#if USE_SCALING_TO_FIX_UNDERFLOW
                    if (scalingFactor != 1.0) {
                        fi[0] *= scalingFactor;
                        fi[1] *= scalingFactor;
                        fi[2] *= scalingFactor;
                        fi[3] *= scalingFactor;
                    }
#endif

                    ns[nodeCount] = cCache[4];

                    theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                                    [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);

                    f2  = theChildNode->nodeIndex-flatLeaves.lLength;

                    if (nm[f2]>=0) {
                        nm[f2] = -1;
                    }
                }
            }
        }
        fromL = topLevelRightL.lData[tlc]+1;
        if (fromL > endLeaf) {
            break;
        }
    }

    // end 20070209 mod

    mm  = topLevelNodes.lData[f];

    for (long nI = 0; nI < f; nI++,mm>>=1) {
        if (mm&1)
            // marked for recalculation
        {
            long startAt = 0;

            if (nI) {
                startAt = topLevelNodes.lData[nI-1]+1;
            }

            if (startAt<startINode) {
                startAt = startINode;
            }

            for (nodeCount=startAt; nodeCount<=topLevelNodes.lData[nI]; nodeCount++) {
                if (nm[nodeCount] == -1) { // marked for update
                    nodeChild = ((node <long>*)(flatNodes.lData[nodeCount]));

                    _Parameter * theseProbs = mlc + 4*(nodeCount+flatLeaves.lLength);
                    theseProbs[0]=LIKELIHOOD_SCALER_INT;
                    theseProbs[1]=LIKELIHOOD_SCALER_INT;
                    theseProbs[2]=LIKELIHOOD_SCALER_INT;
                    theseProbs[3]=LIKELIHOOD_SCALER_INT;
                    nm[nodeCount] = 0;

                    for (long k=0; k <nodeChild->nodes.length; k++) {
                        travNode =
                            ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->nodes.data[k]->in_object]);

                        _Parameter * childProbs = mlc + 4*(travNode->nodeIndex);
                        long         nState = ns[travNode->nodeIndex];

                        fastIndex = travNode->compExp->theData;
                        if (nState<0) {
                            //marginalLFEvalsAmb ++;

                            _Parameter a = *(childProbs),
                                       b = childProbs[1],
                                       c = childProbs[2],
                                       d = childProbs[3];

                            _Parameter tmp  =    a*  fastIndex[0]
                                                 +b*  fastIndex[1]
                                                 +c*  fastIndex[2]
                                                 +d*  fastIndex[3],
                                                 tmp2 =  a*  fastIndex[4]
                                                         +b*  fastIndex[5]
                                                         +c*  fastIndex[6]
                                                         +d*  fastIndex[7],
                                                         tmp3 =     a*  fastIndex[8]
                                                                    +b*  fastIndex[9]
                                                                    +c*  fastIndex[10]
                                                                    +d*  fastIndex[11],
                                                                    tmp4 =    a*  fastIndex[12]
                                                                            +b*  fastIndex[13]
                                                                            +c*  fastIndex[14]
                                                                            +d*  fastIndex[15];

                            theseProbs[0] *= tmp;
                            theseProbs[1] *= tmp2;
                            theseProbs[2] *= tmp3;
                            theseProbs[3] *= tmp4;

                        } else {
                            //marginalLFEvals ++;
                            fastIndex     += nState;
#if USE_SCALING_TO_FIX_UNDERFLOW
                            theseProbs[0] *= fastIndex[0] * scalingFactor;
                            theseProbs[1] *= fastIndex[4] * scalingFactor;
                            theseProbs[2] *= fastIndex[8] * scalingFactor;
                            theseProbs[3] *= fastIndex[12]* scalingFactor;

#else
                            theseProbs[0] *= fastIndex[0];
                            theseProbs[1] *= fastIndex[4];
                            theseProbs[2] *= fastIndex[8];
                            theseProbs[3] *= fastIndex[12];
#endif
                        }
                    }
                    nm[((_CalcNode*)((BaseRef*)variablePtrs.lData)
                        [nodeChild->parent->in_object])->nodeIndex-flatLeaves.lLength] = -1;

#if USE_SCALING_TO_FIX_UNDERFLOW
                    if (theseProbs[0] == 0.0 && theseProbs[1] == 0.0 && theseProbs[2] == 0.0 && theseProbs[3] == 0.0) {
                        return doChar4Scaling (dsf, index,  position, offset,fixAttempt,localScalingFactor,overflowFlag,false);
                    } else if (theseProbs[0] >= 1e500 || theseProbs[1] >= 1e500 || theseProbs[2] >= 1e500 || theseProbs[3] >= 1e500 ) {
                        return doChar4Scaling (dsf, index,  position, offset,fixAttempt,localScalingFactor,overflowFlag,true);
                    }
#endif

                }
            }

            theChildNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[topLevelNodes.lData[nI]]);

            _Parameter* nodeProbs =  mlc+4*theChildNode->nodeIndex;

            fastIndex = rootIChildrenCache+position*4*(topLevelNodes.lLength-1)+4*nI;

            fastIndex[0] = nodeProbs[0];
            fastIndex[1] = nodeProbs[1];
            fastIndex[2] = nodeProbs[2];
            fastIndex[3] = nodeProbs[3];
        } else {
            _Parameter* nodeProbs = mlc+
                                    4*((_CalcNode*)(((BaseRef*)flatTree.lData)[topLevelNodes.lData[nI]]))->nodeIndex;
            fastIndex = rootIChildrenCache+position*4*(topLevelNodes.lLength-1)+4*nI;
            nodeProbs[0] = fastIndex[0];
            nodeProbs[1] = fastIndex[1];
            nodeProbs[2] = fastIndex[2];
            nodeProbs[3] = fastIndex[3];
            nodeChild = ((node <long>*)(flatNodes.lData[topLevelNodes.lData[nI]]))->parent;
            nm[((_CalcNode*)((BaseRef*)variablePtrs.lData)
                [nodeChild->in_object])->nodeIndex-flatLeaves.lLength] = -1;
        }
    }

    for (nodeCount=topLevelNodes.lData[f-1]+1; nodeCount<flatTree.lLength; nodeCount++) {
        if (nm[nodeCount] == -1) { // marked for update
            nodeChild = ((node <long>*)(flatNodes.lData[nodeCount]));

            _Parameter * theseProbs = mlc + 4*(nodeCount+flatLeaves.lLength);
            theseProbs[0]=LIKELIHOOD_SCALER;
            theseProbs[1]=LIKELIHOOD_SCALER;
            theseProbs[2]=LIKELIHOOD_SCALER;
            theseProbs[3]=LIKELIHOOD_SCALER;
            nm[nodeCount] = 0;

            for (long k=0; k <nodeChild->nodes.length; k++) {
                travNode =
                    ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->nodes.data[k]->in_object]);

                _Parameter * childProbs = mlc + 4*(travNode->nodeIndex);
                long         nState = ns[travNode->nodeIndex];

                fastIndex = travNode->compExp->theData;
                if (nState<0) {
                    //marginalLFEvalsAmb++;

                    _Parameter a = *(childProbs),
                               b = childProbs[1],
                               c = childProbs[2],
                               d = childProbs[3];

                    _Parameter tmp  =    a*  fastIndex[0]
                                         +b*  fastIndex[1]
                                         +c*  fastIndex[2]
                                         +d*  fastIndex[3],
                                         tmp2 =  a*  fastIndex[4]
                                                 +b*  fastIndex[5]
                                                 +c*  fastIndex[6]
                                                 +d*  fastIndex[7],
                                                 tmp3 =     a*  fastIndex[8]
                                                            +b*  fastIndex[9]
                                                            +c*  fastIndex[10]
                                                            +d*  fastIndex[11],
                                                            tmp4 =    a*  fastIndex[12]
                                                                    +b*  fastIndex[13]
                                                                    +c*  fastIndex[14]
                                                                    +d*  fastIndex[15];

                    theseProbs[0] *= tmp;
                    theseProbs[1] *= tmp2;
                    theseProbs[2] *= tmp3;
                    theseProbs[3] *= tmp4;
                } else {
                    //marginalLFEvals++;

                    fastIndex+=nState;
#if USE_SCALING_TO_FIX_UNDERFLOW
                    theseProbs[0] *= fastIndex[0] * scalingFactor;
                    theseProbs[1] *= fastIndex[4] * scalingFactor;
                    theseProbs[2] *= fastIndex[8] * scalingFactor;
                    theseProbs[3] *= fastIndex[12]* scalingFactor;

#else
                    theseProbs[0] *= fastIndex[0];
                    theseProbs[1] *= fastIndex[4];
                    theseProbs[2] *= fastIndex[8];
                    theseProbs[3] *= fastIndex[12];
#endif
                }
            }

#if USE_SCALING_TO_FIX_UNDERFLOW
            if (theseProbs[0] == 0.0 && theseProbs[1] == 0.0 && theseProbs[2] == 0.0 && theseProbs[3] == 0.0) {
                return doChar4Scaling (dsf, index,  position, offset,fixAttempt,localScalingFactor,overflowFlag,false);
            } else if (theseProbs[0] >= 1e500 || theseProbs[1] >= 1e500 || theseProbs[2] >= 1e500 || theseProbs[3] >= 1e500 ) {
                return doChar4Scaling (dsf, index , position, offset,fixAttempt,localScalingFactor,overflowFlag,true);
            }
#endif

            nodeChild = nodeChild->parent;
            if (nodeChild)
                nm[((_CalcNode*)((BaseRef*)variablePtrs.lData)
                    [nodeChild->in_object])->nodeIndex-flatLeaves.lLength] = -1;
        }
    }

    //theChildNode->cBase = 4;

    fastIndex = mlc+4*(flatLeaves.lLength+flatNodes.lLength-1);

    _Parameter  result = theProbs[0]*fastIndex[0]+
                         theProbs[1]*fastIndex[1]+
                         theProbs[2]*fastIndex[2]+
                         theProbs[3]*fastIndex[3];

#if USE_SCALING_TO_FIX_UNDERFLOW
    if (result <= 0.0) { // underflow
        return doChar4Scaling (dsf, index, position, offset,fixAttempt,localScalingFactor,overflowFlag,false);
    } else if (result >= 1e500) { // overflow
        return doChar4Scaling (dsf, index, position, offset,fixAttempt,localScalingFactor,overflowFlag,true);
    }

#else
    if (result<=0.0) {
        return ALMOST_ZERO;
    }
#endif

    return result;
}
//_______________________________________________________________________________________________

#if USE_SCALING_TO_FIX_UNDERFLOW
_Parameter   _TheTree::ReleafTreeChar4 (_DataSetFilter* dsf, long index, long lastIndex, long startLeaf, long endLeaf, long fixAttempt, _Parameter localScalingFactor)
#else
_Parameter   _TheTree::ReleafTreeChar4 (_DataSetFilter* dsf, long index, long lastIndex, long startLeaf, long endLeaf)
#endif

// set leaf value and re-prune w/o reexping
// for subsequent entries into a datafilter
{
    _CalcNode* travNode,
               * theChildNode ;

    long nodeCount,f;

    _Parameter* fastIndex;
    node<long>* nodeChild;

    char *pastState = nil,
          *thisState = dsf->GetColumn(index);

    if (lastIndex>=0) {
        pastState = dsf->GetColumn(lastIndex);
    }

#if USE_SCALING_TO_FIX_UNDERFLOW
    if (lastIndex>=0 && scalingForUnderflow->theData[index] != scalingForUnderflow->theData[lastIndex]) {
        lastIndex = -1;
        startLeaf = 0;
        endLeaf   = flatLeaves.lLength-1;
    }
    _Parameter scalingFactor = scalingForUnderflow->theData[index];
    bool            overflowFlag = fixAttempt<0;
    if (overflowFlag) {
        fixAttempt = -fixAttempt;
    }
#endif

    for (nodeCount = startLeaf; nodeCount<=endLeaf; nodeCount++) {
        f = dsf->theNodeMap.lData[nodeCount];
        if ( lastIndex < 0 || thisState[f]!=pastState[f]
#if USE_SCALING_TO_FIX_UNDERFLOW
                || scalingFactor != scalingForUnderflow->theData[lastIndex]
#endif
           ) {
            travNode     = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
            long* cCache = dsf->conversionCache.lData+(thisState[f]-40)*5;

#if USE_SCALING_TO_FIX_UNDERFLOW
            travNode->theProbs[0] = *(cCache++)*scalingFactor;
            travNode->theProbs[1] = *(cCache++)*scalingFactor;
            travNode->theProbs[2] = *(cCache++)*scalingFactor;
            travNode->theProbs[3] = *(cCache++)*scalingFactor;

#else
            travNode->theProbs[0] = *(cCache++);
            travNode->theProbs[1] = *(cCache++);
            travNode->theProbs[2] = *(cCache++);
            travNode->theProbs[3] = *(cCache++);
#endif

            travNode->lastState = *cCache;

            theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                            [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);
            if (theChildNode->cBase>0) {
                theChildNode->cBase = -1;
            }
        }
    }

    for (nodeCount= leftiNodes.lData[startLeaf]; nodeCount<flatTree.lLength; nodeCount++) {
        theChildNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[nodeCount]);
        if (theChildNode->cBase == -1) { // marked for update
            nodeChild = ((node <long>*)(flatNodes.lData[nodeCount]));
            theChildNode->theProbs[0]=LIKELIHOOD_SCALER;
            theChildNode->theProbs[1]=LIKELIHOOD_SCALER;
            theChildNode->theProbs[2]=LIKELIHOOD_SCALER;
            theChildNode->theProbs[3]=LIKELIHOOD_SCALER;
            theChildNode->cBase = 4;

            for (long k=0; k <nodeChild->nodes.length; k++) {
                travNode =
                    ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->nodes.data[k]->in_object]);
                fastIndex = travNode->compExp->theData;
                if (travNode->lastState<0) {
                    _Parameter a = *(travNode->theProbs),
                               b = travNode->theProbs[1],
                               c = travNode->theProbs[2],
                               d = travNode->theProbs[3];

                    _Parameter tmp = a * *(fastIndex++);
                    tmp +=b* *(fastIndex++);
                    tmp +=c* *(fastIndex++);
                    theChildNode->theProbs[0]*=tmp+d* *(fastIndex++);

                    tmp = a * *(fastIndex++);
                    tmp +=b* *(fastIndex++);
                    tmp +=c* *(fastIndex++);
                    theChildNode->theProbs[1]*=tmp+d* *(fastIndex++);

                    tmp = a * *(fastIndex++);
                    tmp +=b* *(fastIndex++);
                    tmp +=c* *(fastIndex++);
                    theChildNode->theProbs[2]*=tmp+d* *(fastIndex++);

                    tmp = a * *(fastIndex++);
                    tmp +=b* *(fastIndex++);
                    tmp +=c* *(fastIndex++);
                    theChildNode->theProbs[3]*=tmp+d* *(fastIndex++);
                } else {
                    fastIndex+=travNode->lastState;
#if USE_SCALING_TO_FIX_UNDERFLOW
                    theChildNode->theProbs[0]*=*fastIndex   * scalingFactor;
                    theChildNode->theProbs[1]*=fastIndex[4] * scalingFactor;
                    theChildNode->theProbs[2]*=fastIndex[8] * scalingFactor;
                    theChildNode->theProbs[3]*=fastIndex[12]* scalingFactor;
#else
                    theChildNode->theProbs[0]*=*fastIndex;
                    theChildNode->theProbs[1]*=fastIndex[4];
                    theChildNode->theProbs[2]*=fastIndex[8];
                    theChildNode->theProbs[3]*=fastIndex[12];
#endif
                }
            }
            travNode->cBase = 4;
            nodeChild = nodeChild->parent;
            if (nodeChild) {
                ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->in_object])->cBase = -1;
            }

#if USE_SCALING_TO_FIX_UNDERFLOW
            if (theChildNode->theProbs[0] == 0.0 && theChildNode->theProbs[1] == 0.0 && theChildNode->theProbs[2] == 0.0 && theChildNode->theProbs[3] == 0.0) {
                return doChar4Scaling_nc (dsf, index, fixAttempt,localScalingFactor,overflowFlag,false);
            } else if (theChildNode->theProbs[0] == HUGE_VAL || theChildNode->theProbs[1] == HUGE_VAL || theChildNode->theProbs[2] == HUGE_VAL || theChildNode->theProbs[3] == HUGE_VAL ) {
                return doChar4Scaling_nc (dsf, index,fixAttempt,localScalingFactor,overflowFlag,true);
            }
#endif


        }
    }

    theChildNode->cBase = 4;
    _Parameter  result = theProbs[0]*theChildNode->theProbs[0]+
                         theProbs[1]*theChildNode->theProbs[1]+
                         theProbs[2]*theChildNode->theProbs[2]+
                         theProbs[3]*theChildNode->theProbs[3];


#if USE_SCALING_TO_FIX_UNDERFLOW
    if (result == 0.0) {
        return doChar4Scaling_nc (dsf, index, fixAttempt,localScalingFactor,overflowFlag,false);
    } else if (result == HUGE_VAL) {
        return doChar4Scaling_nc (dsf, index,fixAttempt,localScalingFactor,overflowFlag,true);
    }
#else
    if (result<=0.0) {
        return ALMOST_ZERO;
    }
#endif
    return result;
}



//_______________________________________________________________________________________________

_Parameter   _TheTree::ReleafTreeChar4Degenerate (_DataSetFilter* dsf, long index)
// set leaf value and re-prune w/o reexping
// for subsequent entries into a datafilter
{

    _CalcNode* travNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[theRoot->in_object]),
               * theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[theRoot->nodes.data[0]->in_object]);

    char     *thisState = dsf->GetColumn(index);
    long     * rootState = dsf->conversionCache.lData+(thisState[dsf->theNodeMap.lData[0]]-40)*5,
               * tipState = dsf->conversionCache.lData+(thisState[dsf->theNodeMap.lData[1]]-40)*5,
                 nodeCount = rootState[4],
                 f = tipState[4];

    _Matrix   * theTP = theChildNode->GetCompExp();

    _Parameter reslt,
               tmp,
               * fastIdx = theTP->fastIndex(),
                 * nodeProbs;

    if ((nodeCount>=0)&&(f>=0)) { // 1-1
        reslt =  fastIdx[nodeCount*4+f]*theProbs[nodeCount];
    } else if (nodeCount>=0) { // 1-many
        fastIdx += nodeCount*cBase;

        tmp=*fastIdx* *tipState;
        tmp+=fastIdx[1]*tipState[1];
        tmp+=fastIdx[2]*tipState[2];
        tmp+=fastIdx[3]*tipState[3];

        reslt = theProbs[nodeCount]*tmp;
    } else if (f>=0) { // many to 1
        fastIdx += f;
        nodeProbs  = travNode->theProbs;
        *nodeProbs = *fastIdx** rootState;
        nodeProbs[1] = fastIdx[4]* rootState[1];
        nodeProbs[2] = fastIdx[8]* rootState[2];
        nodeProbs[3] = fastIdx[12]* rootState[3];

        reslt = *theProbs ** nodeProbs+
                theProbs[1]*nodeProbs[1]+
                theProbs[2]*nodeProbs[2]+
                theProbs[3]*nodeProbs[3];

    } else { // many to many
        nodeProbs  = travNode->theProbs;
        nodeProbs[0] = (fastIdx[0]*tipState[0]+
                        fastIdx[1]*tipState[1]+
                        fastIdx[2]*tipState[2]+
                        fastIdx[3]*tipState[3])*rootState[0];
        nodeProbs[1] = (fastIdx[4]*tipState[0]+
                        fastIdx[5]*tipState[1]+
                        fastIdx[6]*tipState[2]+
                        fastIdx[7]*tipState[3])*rootState[1];
        nodeProbs[2] = (fastIdx[8]*tipState[0]+
                        fastIdx[9]*tipState[1]+
                        fastIdx[10]*tipState[2]+
                        fastIdx[11]*tipState[3])*rootState[2];
        nodeProbs[3] = (fastIdx[12]*tipState[0]+
                        fastIdx[13]*tipState[1]+
                        fastIdx[14]*tipState[2]+
                        fastIdx[15]*tipState[3])*rootState[3];

        reslt = *theProbs ** nodeProbs+
                theProbs[1]*nodeProbs[1]+
                theProbs[2]*nodeProbs[2]+
                theProbs[3]*nodeProbs[3];
    }
    if (reslt<=0.0) {
        reslt = ALMOST_ZERO;
    }
    //printf ("%d\t%g\n", index, reslt);
    return reslt;
}

//_______________________________________________________________________________________________

_Parameter   _TheTree::ThreadReleafTreeCharCache (_DataSetFilter* dsf, long index, long lastIndex, long startLeaf, long endLeaf,
        long position,long offset)
// set leaf value and re-prune w/o reexping
// for subsequent entries into a datafilter
{


    _CalcNode *travNode,
              *theChildNode;

    long       nodeCount,
               f,
               *ns  = nodeStates+(flatLeaves.lLength+flatNodes.lLength)*offset;

    _Parameter*fastIndex,
               *stopper,
               *mlc = marginalLikelihoodCache+cBase*(flatLeaves.lLength+flatNodes.lLength)*offset;

    node<long>* nodeChild;


    char *pastState = nil,
          *thisState = dsf->GetColumn(index),
           *nm  = nodeMarkers+flatNodes.lLength*offset;

    if (lastIndex>=0) {
        pastState = dsf->GetColumn(lastIndex);
    }

    for (nodeCount = startLeaf; nodeCount<=endLeaf; nodeCount++) {
        f = dsf->theNodeMap.lData[nodeCount];
        if ((lastIndex==-1)||(thisState[f]!=pastState[f])) {
            ns[nodeCount] = dsf->LookupConversion (thisState[f], mlc+cBase*nodeCount);
            theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                            [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);
            nm[theChildNode->nodeIndex-flatLeaves.lLength] = -1;
        }
    }


    long lm1 = topLevelNodes.lLength-1,
         mm = topLevelNodes.lData[lm1],
         startINode = leftiNodes.lData[startLeaf];

    for (long nI = 0; nI < lm1; nI++,mm>>=1) {
        if (mm&1)
            // marked for recalculation
        {
            long startAt = 0;

            if (nI) {
                startAt = topLevelNodes.lData[nI-1]+1;
            }

            if (startAt<startINode) {
                startAt = startINode;
            }

            for (nodeCount=startAt; nodeCount<=topLevelNodes.lData[nI]; nodeCount++) {
                if (nm[nodeCount]==-1) { // marked for update
                    _Parameter * theseProbs = fastIndex = mlc+cBase*(nodeCount+flatLeaves.lLength);
                    stopper =   fastIndex+cBase;

                    for (; fastIndex!=stopper; *(fastIndex++)=1.) ;

                    nodeChild = ((node <long>*)(flatNodes.lData[nodeCount]));

                    nm[nodeCount] = 0;

                    for (long k=0; k <nodeChild->nodes.length; k++) {
                        travNode =
                            ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->nodes.data[k]->in_object]);

                        _Parameter* cProbs = mlc+cBase*travNode->nodeIndex,
                                    * stopper2 = cProbs+cBase;
                        fastIndex = travNode->compExp->theData;

                        f = ns[travNode->nodeIndex];
                        if (f<0) {
                            for (long i=0; i<cBase; i++) {
                                _Parameter tmp = *(cProbs) * *fastIndex;
                                fastIndex++;
                                stopper = cProbs+1;
                                for (; stopper!=stopper2; fastIndex++, stopper++) {
                                    tmp += *stopper * (*fastIndex);
                                }
                                theseProbs[i]*=tmp;
                            }
                        } else {
                            fastIndex+=f;
                            _Parameter tmp = cProbs[f];
                            for (long i=0; i<cBase; i++,fastIndex+=cBase) {
                                theseProbs[i]*=tmp* (*fastIndex);
                            }
                        }
                    }
                    nodeChild = nodeChild->parent;
                    if (nodeChild) {
                        nm[((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->in_object])->nodeIndex-flatLeaves.lLength] = -1;
                    }
                }
            }

            theChildNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[topLevelNodes.lData[nI]]);
            _Parameter* nodeProbs = mlc+cBase*theChildNode->nodeIndex;
            fastIndex = rootIChildrenCache+position*cBase*(topLevelNodes.lLength-1)+cBase*nI;
            for (long i=0; i<cBase; i++) {
                fastIndex[i] = nodeProbs[i];
            }
        } else {
            _Parameter* nodeProbs = mlc+cBase*((_CalcNode*)(((BaseRef*)flatTree.lData)[topLevelNodes.lData[nI]]))->nodeIndex;
            fastIndex = rootIChildrenCache+position*cBase*(topLevelNodes.lLength-1)+cBase*nI;
            for (long i=0; i<cBase; i++) {
                nodeProbs[i] =fastIndex[i];
            }
            nodeChild = ((node <long>*)(flatNodes.lData[topLevelNodes.lData[nI]]))->parent;
            nm[((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->in_object])->nodeIndex-flatLeaves.lLength] = -1;
        }
    }

    for (nodeCount=topLevelNodes.lData[lm1-1]+1; nodeCount<flatTree.lLength; nodeCount++) {
        if (nm[nodeCount]==-1) { // marked for update
            _Parameter * theseProbs = fastIndex = mlc+cBase*(nodeCount+flatLeaves.lLength);

            for (long zeroPop = 0; zeroPop < cBase; zeroPop++) {
                fastIndex[zeroPop] = 1.;
            }

            nodeChild = ((node <long>*)(flatNodes.lData[nodeCount]));

            nm[nodeCount] = 0;

            for (long k=0; k <nodeChild->nodes.length; k++) {
                travNode =
                    ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->nodes.data[k]->in_object]);

                _Parameter* cProbs = mlc+cBase*travNode->nodeIndex,
                            * stopper2 = cProbs+cBase;

                fastIndex = travNode->compExp->theData;

                f = ns[travNode->nodeIndex];
                if (f<0) {
                    for (long i=0; i<cBase; i++) {
                        _Parameter tmp = *(cProbs) * *fastIndex;
                        fastIndex++;
                        stopper = cProbs+1;
                        for (; stopper!=stopper2; fastIndex++, stopper++) {
                            tmp += *stopper * (*fastIndex);
                        }
                        theseProbs[i]*=tmp;
                    }
                    /*if (rem==1)
                    {
                        for (long i=0; i<cBase; i++)
                        {
                            _Parameter tmp1  = 0.0,
                                       tmp2  = 0.0,
                                       tmp3  = 0.0,
                                       tmp4  = 0.0;

                            for (long k=0; k<upTo; k+=4)
                            {
                                tmp1 += cProbs[k]  *fastIndex[k];
                                tmp2 += cProbs[k+1]*fastIndex[k+1];
                                tmp3 += cProbs[k+2]*fastIndex[k+2];
                                tmp4 += cProbs[k+3]*fastIndex[k+3];
                            }

                            theseProbs[i] *= tmp1 + cProbs[upTo]*fastIndex[upTo] + tmp2+tmp3+tmp4;
                            fastIndex += cBase;
                        }
                    }
                    else
                        if (rem==2)
                        {
                            for (long i=0; i<cBase; i++)
                            {
                                _Parameter tmp1  = 0.0,
                                           tmp2  = 0.0,
                                           tmp3  = 0.0,
                                           tmp4  = 0.0;

                                for (long k=0; k<upTo; k+=4)
                                {
                                    tmp1 += cProbs[k]  *fastIndex[k];
                                    tmp2 += cProbs[k+1]*fastIndex[k+1];
                                    tmp3 += cProbs[k+2]*fastIndex[k+2];
                                    tmp4 += cProbs[k+3]*fastIndex[k+3];
                                }

                                theseProbs[i] *= tmp1+tmp2+tmp3+tmp4+
                                                 + cProbs[upTo]*fastIndex[upTo] +
                                                 + cProbs[upTo+1]*fastIndex[upTo+1];
                                fastIndex += cBase;
                            }
                        }
                        else
                            if (rem==3)
                                for (long i=0; i<cBase; i++)
                                {
                                    _Parameter tmp1  = 0.0,
                                               tmp2  = 0.0,
                                               tmp3  = 0.0,
                                               tmp4  = 0.0;

                                    for (long k=0; k<upTo; k+=4)
                                    {
                                        tmp1 += cProbs[k]  *fastIndex[k];
                                        tmp2 += cProbs[k+1]*fastIndex[k+1];
                                        tmp3 += cProbs[k+2]*fastIndex[k+2];
                                        tmp4 += cProbs[k+3]*fastIndex[k+3];
                                    }

                                    theseProbs[i] *= tmp1+tmp2+tmp3+tmp4+
                                                     + cProbs[upTo]*fastIndex[upTo] +
                                                     + cProbs[upTo+1]*fastIndex[upTo+1];
                                                     + cProbs[upTo+2]*fastIndex[upTo+2];

                                    fastIndex += cBase;
                                }
                            else
                                for (long i=0; i<cBase; i++)
                                {
                                    _Parameter tmp1  = 0.0,
                                               tmp2  = 0.0,
                                               tmp3  = 0.0,
                                               tmp4  = 0.0;

                                    for (long k=0; k<upTo; k+=4)
                                    {
                                        tmp1 += cProbs[k]  *fastIndex[k];
                                        tmp2 += cProbs[k+1]*fastIndex[k+1];
                                        tmp3 += cProbs[k+2]*fastIndex[k+2];
                                        tmp4 += cProbs[k+3]*fastIndex[k+3];
                                    }

                                    theseProbs[i] *= tmp1+tmp2+tmp3+tmp4;
                                    fastIndex     += cBase;
                                }       */

                } else {
                    fastIndex+=f;
                    _Parameter tmp = cProbs[f];
                    for (long i=0; i<cBase; i++,fastIndex+=cBase) {
                        theseProbs[i]*=tmp* (*fastIndex);
                    }
                }
            }
            nodeChild = nodeChild->parent;
            if (nodeChild) {
                nm[((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->in_object])->nodeIndex-flatLeaves.lLength] = -1;
            }
        }
    }

    _Parameter result = 0.0;
    fastIndex = mlc+cBase*(flatLeaves.lLength+flatTree.lLength-1);

    for (long i=0; i<cBase; i++) {
        result += theProbs[i] * fastIndex[i];
    }

    if (result<=0.0) {
        return ALMOST_ZERO;
    }
    return result;
}

//_______________________________________________________________________________________________

bool     _TheTree::IntPopulateLeaves (_DataSetFilter* dsf, long index, long)
// assign proper values to leaf conditional probability vectors
{
    bool allGaps = true;

    for (long nodeCount = 0; nodeCount<flatLeaves.lLength; nodeCount++)
        //if (lastIndex<0 || !dsf->CompareTwoSites(lastIndex,index,nodeCount))
    {
        _CalcNode * travNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
        allGaps &= ((travNode->lastState = dsf->Translate2Frequencies ((*dsf)(index,nodeCount), travNode->theProbs, true))<0);
        if (allGaps) // check to see if all
            for (long b = 0; b < cBase; b++)
                if (travNode->theProbs[b] == 0.0) {
                    allGaps = false;
                    break;
                }

        _CalcNode * theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                                    [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);
        theChildNode->cBase = -1;
    }

    return allGaps;
}


//_______________________________________________________________________________________________

void     _TheTree::RecoverNodeSupportStates (_DataSetFilter* dsf, long index, long lastIndex, _Matrix& resultMatrix)
// assume current values of all parameters
// return 2 sets of vectors for each branch
//   - top-down  conditionals
//   - bottom-up conditionals
//   resultMatrix is assumed to contain
//      uniqueSites X (flatLeaves.lLength+flatTree.lLength)*cBase*2 X categoryCount
{

    long      globalShifter        = (flatLeaves.lLength+flatTree.lLength)*cBase,
              catShifer             = dsf->NumberDistinctSites() * 2 * globalShifter;

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

_List* _TheTree::MapCBaseToCharacters (_DataSetFilter* dsf, bool normalanc)
// 20090203/SLKP: clean up
{
    _List               * reply = new _List;
    checkPointer         (reply);

    long unitLength = dsf->GetUnitLength();

    for (long pad = 0; pad<unitLength ; pad ++) {
        reply->AppendNewInstance (new _String (5,true));
    }

    _CalcNode*  theChildNode = StepWiseTraversal (true);
    _String     rootValue = dsf->ConvertCodeToLetters (dsf->CorrectCode(theChildNode->cBase), unitLength);

    for (long f=0; f<unitLength; f++) {
        *((_String*)(*reply)(f)) << rootValue[f];
    }

    theChildNode = StepWiseTraversal (false);

    while (theChildNode) {
        if (!IsCurrentNodeATip()) {
            if (normalanc) {
                _CalcNode* travNode              =  (_CalcNode*)((BaseRef*)variablePtrs.lData)[currentNode->parent->in_object];
                theChildNode->cBase  = theChildNode->theProbs[travNode->cBase];
                theChildNode->categoryVariables.Delete(theChildNode->categoryVariables.lLength-1);
            }

            _String  letterValue = dsf->ConvertCodeToLetters (dsf->CorrectCode(theChildNode->cBase), unitLength);
            for (long i=0; i<unitLength; i++) {
                *((_String*)(*reply)(i)) << letterValue[i];
            }
        }
        theChildNode = StepWiseTraversal (false);
    }

    for (long j = 0; j<unitLength ; j ++) {
        ((_String*)(*reply)(j))->Finalize();
    }

    return reply;
}

//_______________________________________________________________________________________________

_Parameter  _TheTree::ConditionalNodeLikelihood   (node<long>* parentNode, node<long>* thisNode, _Parameter* scoreBelow, _Parameter* myScore, long myState, long offset)
{
    if (parentNode) {
        _Parameter       conditionalLikelihood = LIKELIHOOD_SCALER;

        for (long child = 0; child<thisNode->nodes.length; child++) {
            _CalcNode *                childCNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[thisNode->nodes.data[child]->in_object]);
            conditionalLikelihood      *= *(childCNode->compExp->theData + myState*cBase + childCNode->cBase) * childCNode->theValue;
        }

        myScore[myState] = conditionalLikelihood;

        return ConditionalBranchLikelihood (parentNode, thisNode, myScore, scoreBelow, -1, offset);
    } else {
        _Parameter      tLik = theProbs[myState];

        for (long child = 0; child<thisNode->nodes.length; child++) {
            _CalcNode *      childCNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[thisNode->nodes.data[child]->in_object]);
            tLik      *= *(childCNode->compExp->theData + myState*cBase + childCNode->cBase) * childCNode->theValue;
        }

        return tLik;
    }
}

//_______________________________________________________________________________________________

_Parameter  _TheTree::ConditionalBranchLikelihood   (node<long>* parentNode, node<long>* thisNode, _Parameter* scoreBelow, _Parameter* myScore, long myState, long offset)
{
    for (long pstate = (myState>=0?myState:0); pstate < (myState>=0?myState+1:cBase); pstate ++) {
        _Parameter       conditionalLikelihood = LIKELIHOOD_SCALER;
        for (long child = 0; child<parentNode->nodes.length; child++) {
            _Parameter       childContribution = 0.0;

            node<long>*      childNode  = parentNode->nodes.data[child];
            _CalcNode *      childCNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[childNode->in_object]);

            _Parameter*      conds,
                             *       transitions = childCNode->compExp->theData + pstate*cBase,
                                     *      stopper;

            if (childNode==thisNode) { // use scoreBelow
                conds = scoreBelow;
            } else { // use existing conditionals
                if (offset>=0) {
                    conds = marginalLikelihoodCache+cBase*(flatLeaves.lLength+flatNodes.lLength)*offset + cBase*((long)childCNode->theProbs[0]);
                } else {
                    conds = childCNode->theProbs;
                }
            }

            long      remaind = cBase%4;

            if (remaind) {
                stopper = conds + (cBase-remaind);
                for (; conds!= stopper; ) {
                    childContribution += *transitions * *conds;
                    transitions++;
                    conds++;
                    childContribution += *transitions * *conds;
                    transitions++;
                    conds++;
                    childContribution += *transitions * *conds;
                    transitions++;
                    conds++;
                    childContribution += *transitions * *conds;
                    transitions++;
                    conds++;
                }

                switch (remaind) {
                case 1: {
                    childContribution += *(transitions) * *(conds);
                    break;
                }
                case 2: {
                    childContribution += *transitions * *conds;
                    transitions++;
                    conds++;
                    childContribution += *transitions * *conds;
                    break;
                }
                case 3: {
                    childContribution += *transitions * *conds;
                    transitions++;
                    conds++;
                    childContribution += *transitions * *conds;
                    transitions++;
                    conds++;
                    childContribution += *transitions * *conds;
                    break;
                }
                }
            } else {
                stopper = conds + cBase;
                for (; conds!= stopper; ) {
                    childContribution += *transitions * *conds;
                    transitions++;
                    conds++;
                    childContribution += *transitions * *conds;
                    transitions++;
                    conds++;
                    childContribution += *transitions * *conds;
                    transitions++;
                    conds++;
                    childContribution += *transitions * *conds;
                    transitions++;
                    conds++;
                }
            }

            conditionalLikelihood *= childContribution;

            if (conditionalLikelihood == 0.0) {
                if (myState >= 0) {
                    return 0.0;
                }
                break;
            }
        }
        myScore[pstate] = conditionalLikelihood;
    }

    if (parentNode->parent) {
        return ConditionalBranchLikelihood (parentNode->parent, parentNode, myScore, scoreBelow, -1, offset);
    }

    if (myState < 0) {
        _Parameter       treeLikelihood = 0.0;

        for (long pstate = 0; pstate < cBase; pstate ++) {
            treeLikelihood += theProbs[pstate] * myScore [pstate];
        }

        return          treeLikelihood;
    }
    return theProbs[myState]*myScore[myState];
}

//_______________________________________________________________________________________________

void     _TheTree::WeightedCharacterDifferences (_Parameter siteLikelihood, _Matrix* res1, _Matrix* res2, long offset)
// presumes that compExp is populated and so are all the state vectors in internal nodes
{
    if (cBase>128) {
        WarnError ("State spaces with more than 128 states are not supported in WeightedCharacterDifferences");
        return;
    }

    _Parameter vec1[128],
               vec2[128];

    //checkPointer (vec1);
    //checkPointer (vec2);

    //checkParameter (VerbosityLevelString,verbLevel, 0.0);
    //long          branchCount = 0;

    for (long       nodeCount=0; nodeCount<flatTree.lLength; nodeCount++) {
        node<long>* disNode = ((node <long>*)(flatNodes.lData[nodeCount]));


        for (long   childCount = disNode->nodes.length-1; childCount >= 0; childCount --) {
            node<long>* childNode = disNode->nodes.data[childCount];

            _CalcNode*  cNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[childNode->in_object]);
            _Parameter  weightFactor = cNode->Value(),
                        consistencyCheck = 0.0;


            for (long pState = 0; pState < cBase; pState ++) {
                for (long cState = 0; cState < cBase; cState ++) {
                    for (long eraser = 0; eraser < cBase; eraser++) {
                        vec1[eraser] = 0.;
                        vec2[eraser] = 0.;
                    }

                    if (offset>=0) {
                        vec1[cState] = *(marginalLikelihoodCache+cBase*(flatLeaves.lLength+flatNodes.lLength)*offset + cBase*((long)cNode->theProbs[0])+cState);
                    } else {
                        vec1[cState] = cNode->theProbs[cState];
                    }


                    _Parameter  cLik = (ConditionalBranchLikelihood (disNode, childNode, vec1, vec2, pState, offset)/siteLikelihood);
                    res1->theData[pState*cBase+cState] += cLik;
                    res2->theData[pState*cBase+cState] += cLik*weightFactor;
                    consistencyCheck += cLik;
                }
            }

#if !defined __UNIX__ || defined __HEADLESS__
            if ((cBase>=20)&&(offset<1)) {
                yieldCPUTime();
                if (terminateExecution) {
                    return;
                }
            }
#endif

            if (!CheckEqual(consistencyCheck, 1.)) {
                _String consistencyErra ("Failed Internal Consistency Check In WeightedCharacterDifferences at ");
                consistencyErra = consistencyErra & * LocateVar(disNode->in_object)->GetName() & " and " & *cNode->GetName() &
                                  ". Summed RLS to " & consistencyCheck;
                WarnError (consistencyErra);
            }
        }
    }

    //delete vec1;
    //delete vec2;
}

//_______________________________________________________________________________________________
_AVLListX*  _TheTree::ConstructNodeToIndexMap (bool doINodes)
{
    _SimpleList * nodes  = new _SimpleList,
    * whichL = doINodes?&flatNodes:&flatLeaves;
    _AVLListX   * result = new _AVLListX (nodes);

    for (long   pistolero = 0; pistolero < whichL->lLength; pistolero++) {
        result->Insert ((BaseRef)whichL->lData[pistolero], pistolero, false);
    }

    return        result;

}

//_______________________________________________________________________________________________
void _TheTree::MapPostOrderToInOderTraversal (_SimpleList& storeHere, bool doINodes)
{
    _AVLListX*          nodeMapper    = ConstructNodeToIndexMap (doINodes);
    _CalcNode*          traversalNode = doINodes?StepWiseTraversal(true):DepthWiseTraversal(true);

    long                allNodeCount = 0;

    storeHere.Populate (doINodes?flatTree.lLength:flatLeaves.lLength, 0, 0);

    while (traversalNode) {
        bool isTip = IsCurrentNodeATip();
        if ((isTip && !doINodes) || (!isTip && doINodes)) {
            storeHere.lData[nodeMapper->GetXtra (nodeMapper->Find((BaseRef)(&GetCurrentNode())))] = allNodeCount++;
        }

        traversalNode = doINodes?StepWiseTraversal(false):DepthWiseTraversal(false);
    }

    nodeMapper->DeleteAll(false);
    DeleteObject (nodeMapper);


}

//_______________________________________________________________________________________________
void    _TheTree::MarkDone (void)
{
    _CalcNode*  travNode = StepWiseTraversal (TRUE);

    while (travNode) {
        travNode -> MarkDone();
        travNode = StepWiseTraversal();
    }
}

//_______________________________________________________________________________________________

_Parameter  _TheTree::PruneTree (long categID)
// notice that the limiting atomic probabilites are assumed to be stored in the
// model matrix of the tree itself

{
    // do the depth first traversal to start at the leaves

    _CalcNode* curNode = DepthWiseTraversal (true);
    while (curNode) {
        bool reexpnt = curNode->NeedToExponentiate(categID);
        // flag all the nodes up the chain - to the root

        if (reexpnt&&curNode->GetModelMatrix()) {
            curNode->RecomputeMatrix(categID, categoryCount);
        } else if (categID>=0) {
            curNode->SetCompMatrix (categID);
        }
        // the probabilities below in the tree have already been computed
        long    nNodes = currentNode->get_num_nodes();
        if (nNodes) {
            long j;
            for (j = 0; j<nNodes; j++) {
                _CalcNode* tC = ((_CalcNode*)(LocateVar(currentNode->get_node(j+1)->get_data())));
                if(!tC->GetCompExp(categID)) {
                    tC->RecomputeMatrix(categID, categoryCount);
                } else if (categID>=0) {
                    tC->SetCompMatrix(categID);
                }
            }
        }
        // the probabilities have now been set and the flag erased even if it had been raised before

        curNode = DepthWiseTraversal();
    }

    return 0.0;
}

//_______________________________________________________________________________________________

_Parameter  _TheTree::PruneTreeChar4 (long categID)
// notice that the limiting atomic probabilites are assumed to be stored in the
// model matrix of the tree itself

{

    long            nodeCount;
    _Parameter*     fastIndex;
    node<long>*     nodeChild;
    bool            reexpnt;
    _CalcNode*      theChildNode,
                    *travNode;


    for (nodeCount = 0; nodeCount<flatLeaves.lLength; nodeCount++) {
        theChildNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
        reexpnt = theChildNode->NeedToExponentiate(categID);
        if (reexpnt&&theChildNode->GetModelMatrix()) {
            theChildNode->RecomputeMatrix(categID, categoryCount);
        } else if (categID>=0) {
            theChildNode->SetCompMatrix (categID);
        }
    }

    for (nodeCount=0; nodeCount<flatTree.lLength; nodeCount++) {
        theChildNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[nodeCount]);
        theChildNode->theProbs[0]=LIKELIHOOD_SCALER;
        theChildNode->theProbs[1]=LIKELIHOOD_SCALER;
        theChildNode->theProbs[2]=LIKELIHOOD_SCALER;
        theChildNode->theProbs[3]=LIKELIHOOD_SCALER;

        reexpnt = theChildNode->NeedToExponentiate(categID);
        if (reexpnt&&theChildNode->GetModelMatrix()) {
            theChildNode->RecomputeMatrix(categID, categoryCount);
            theChildNode->cBase = -1;
        } else {
            if (categID>=0) {
                theChildNode->SetCompMatrix (categID);
            }
        }
    }

    for (nodeCount = 0; nodeCount<flatLeaves.lLength; nodeCount++) {
        travNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
        theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                        [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);

        fastIndex = travNode->compExp->theData;
        if (travNode->lastState>=0) {
            fastIndex+=travNode->lastState;
            theChildNode->theProbs[0]*=*fastIndex;
            theChildNode->theProbs[1]*=fastIndex[4];
            theChildNode->theProbs[2]*=fastIndex[8];
            theChildNode->theProbs[3]*=fastIndex[12];
        } else {
            _Parameter a = *(travNode->theProbs),
                       b = travNode->theProbs[1],
                       c = travNode->theProbs[2],
                       d = travNode->theProbs[3];

            _Parameter tmp = a * *(fastIndex++);
            tmp +=b* *(fastIndex++);
            tmp +=c* *(fastIndex++);
            theChildNode->theProbs[0]*=tmp+d* *(fastIndex++);

            tmp = a * *(fastIndex++);
            tmp +=b* *(fastIndex++);
            tmp +=c* *(fastIndex++);
            theChildNode->theProbs[1]*=tmp+d* *(fastIndex++);

            tmp = a * *(fastIndex++);
            tmp +=b* *(fastIndex++);
            tmp +=c* *(fastIndex++);
            theChildNode->theProbs[2]*=tmp+d* *(fastIndex++);

            tmp = a * *(fastIndex++);
            tmp +=b* *(fastIndex++);
            tmp +=c* *(fastIndex++);
            theChildNode->theProbs[3]*=tmp+d* *(fastIndex++);
        }
    }

    for (nodeCount=0; nodeCount<flatTree.lLength-1; nodeCount++) {
        travNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[nodeCount]);
        nodeChild = ((node <long>*)(flatNodes.lData[nodeCount]))->parent;
        theChildNode =
            ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->in_object]);

        fastIndex = travNode->compExp->theData;

        _Parameter a = *(travNode->theProbs),
                   b = travNode->theProbs[1],
                   c = travNode->theProbs[2],
                   d = travNode->theProbs[3];

        _Parameter tmp = a * *(fastIndex++);
        tmp +=b* *(fastIndex++);
        tmp +=c* *(fastIndex++);
        theChildNode->theProbs[0]*=tmp+d* *(fastIndex++);

        tmp = a * *(fastIndex++);
        tmp +=b* *(fastIndex++);
        tmp +=c* *(fastIndex++);
        theChildNode->theProbs[1]*=tmp+d* *(fastIndex++);

        tmp = a * *(fastIndex++);
        tmp +=b* *(fastIndex++);
        tmp +=c* *(fastIndex++);
        theChildNode->theProbs[2]*=tmp+d* *(fastIndex++);

        tmp = a * *(fastIndex++);
        tmp +=b* *(fastIndex++);
        tmp +=c* *(fastIndex++);
        theChildNode->theProbs[3]*=tmp+d* *(fastIndex++);
    }


    theChildNode = (_CalcNode*)(((BaseRef*)variablePtrs.lData)[theRoot->in_object]);

    _Parameter  result = theProbs[0]*theChildNode->theProbs[0]+
                         theProbs[1]*theChildNode->theProbs[1]+
                         theProbs[2]*theChildNode->theProbs[2]+
                         theProbs[3]*theChildNode->theProbs[3];

    if (result<=0.0) {
        return ALMOST_ZERO;
    }
    return result;
}

//_______________________________________________________________________________________________

_Parameter  _TheTree::PruneTreeChar4Cache (long categID)
// notice that the limiting atomic probabilites are assumed to be stored in the
// model matrix of the tree itself

{

    long            nodeCount;
    _Parameter*     fastIndex;
    node<long>*     nodeChild;
    bool            reexpnt;
    _CalcNode*      theChildNode,
                    *travNode;


    for (nodeCount = 0; nodeCount<flatLeaves.lLength; nodeCount++) {
        theChildNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
        reexpnt = theChildNode->NeedToExponentiate(categID);
        if (reexpnt&&theChildNode->GetModelMatrix()) {
            theChildNode->RecomputeMatrix(categID, categoryCount);
            theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                            [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);
            theChildNode->cBase = -1;
        } else if (categID>=0) {
            theChildNode->SetCompMatrix (categID);
        }
    }

    for (nodeCount=0; nodeCount<flatTree.lLength; nodeCount++) {
        theChildNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[nodeCount]);
        theChildNode->theProbs[0]=LIKELIHOOD_SCALER;
        theChildNode->theProbs[1]=LIKELIHOOD_SCALER;
        theChildNode->theProbs[2]=LIKELIHOOD_SCALER;
        theChildNode->theProbs[3]=LIKELIHOOD_SCALER;
        reexpnt = theChildNode->NeedToExponentiate(categID);
        if (reexpnt&&theChildNode->GetModelMatrix()) {
            theChildNode->RecomputeMatrix(categID, categoryCount);
            theChildNode->cBase = -1;
        } else {
            if (categID>=0) {
                theChildNode->SetCompMatrix (categID);
            }
        }
        if (theChildNode->cBase==-1) {
            nodeChild = ((node <long>*)(flatNodes.lData[nodeCount]))->parent;
            if (nodeChild) {
                theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->in_object]);
                theChildNode->cBase = -1;
            }
        }
    }

    for (nodeCount = 0; nodeCount<flatLeaves.lLength; nodeCount++) {
        travNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
        theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                        [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);

        fastIndex = travNode->GetCompExp()->fastIndex();
        if (travNode->lastState>=0) {
            fastIndex+=travNode->lastState;
            theChildNode->theProbs[0]*=*fastIndex;
            theChildNode->theProbs[1]*=fastIndex[4];
            theChildNode->theProbs[2]*=fastIndex[8];
            theChildNode->theProbs[3]*=fastIndex[12];
        } else {
            _Parameter a = *(travNode->theProbs),
                       b = travNode->theProbs[1],
                       c = travNode->theProbs[2],
                       d = travNode->theProbs[3];

            _Parameter tmp = a * *(fastIndex++);
            tmp +=b* *(fastIndex++);
            tmp +=c* *(fastIndex++);
            theChildNode->theProbs[0]*=tmp+d* *(fastIndex++);

            tmp = a * *(fastIndex++);
            tmp +=b* *(fastIndex++);
            tmp +=c* *(fastIndex++);
            theChildNode->theProbs[1]*=tmp+d* *(fastIndex++);

            tmp = a * *(fastIndex++);
            tmp +=b* *(fastIndex++);
            tmp +=c* *(fastIndex++);
            theChildNode->theProbs[2]*=tmp+d* *(fastIndex++);

            tmp = a * *(fastIndex++);
            tmp +=b* *(fastIndex++);
            tmp +=c* *(fastIndex++);
            theChildNode->theProbs[3]*=tmp+d* *(fastIndex++);
        }
    }

    for (nodeCount=0; nodeCount<flatTree.lLength-1; nodeCount++) {
        travNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[nodeCount]);
        nodeChild = ((node <long>*)(flatNodes.lData[nodeCount]))->parent;
        theChildNode =
            ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->in_object]);

        fastIndex = travNode->GetCompExp()->fastIndex();

        _Parameter a = *(travNode->theProbs),
                   b = travNode->theProbs[1],
                   c = travNode->theProbs[2],
                   d = travNode->theProbs[3];

        _Parameter tmp = a * *(fastIndex++);
        tmp +=b* *(fastIndex++);
        tmp +=c* *(fastIndex++);
        theChildNode->theProbs[0]*=tmp+d* *(fastIndex++);

        tmp = a * *(fastIndex++);
        tmp +=b* *(fastIndex++);
        tmp +=c* *(fastIndex++);
        theChildNode->theProbs[1]*=tmp+d* *(fastIndex++);

        tmp = a * *(fastIndex++);
        tmp +=b* *(fastIndex++);
        tmp +=c* *(fastIndex++);
        theChildNode->theProbs[2]*=tmp+d* *(fastIndex++);

        tmp = a * *(fastIndex++);
        tmp +=b* *(fastIndex++);
        tmp +=c* *(fastIndex++);
        theChildNode->theProbs[3]*=tmp+d* *(fastIndex++);
    }

    long c = 0, off = 1;
    for (nodeCount=0; nodeCount<topLevelNodes.lLength-1; nodeCount++, off<<=1) {
        travNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[topLevelNodes.lData[nodeCount]]);
        if (travNode->cBase<=-1) {
            c |= off;
        }
    }
    topLevelNodes.lData[topLevelNodes.lLength-1] = c;
    for (nodeCount=0; nodeCount<flatTree.lLength; nodeCount++) {
        theChildNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[nodeCount]);
        theChildNode->cBase = cBase;
    }

    _Parameter result = 0.0;

    theChildNode = (_CalcNode*)(((BaseRef*)variablePtrs.lData)[theRoot->in_object]);
    //((_CalcNode*)(LocateVar(theRoot->get_data())));

    for (nodeCount=0; nodeCount<cBase; nodeCount++) {
        result+= theProbs[nodeCount]*theChildNode->theProbs[nodeCount];
    }

    if (result<=0.0) {
        return ALMOST_ZERO;
    }

    return result;
}

//_______________________________________________________________________________________________

_Parameter  _TheTree::PruneTreeChar (long categID)
// notice that the limiting atomic probabilites are assumed to be stored in the
// model matrix of the tree itself

{

    long            nodeCount,j;
    _Parameter*     fastIndex;
    node<long>*     nodeChild;
    bool            reexpnt;
    _CalcNode*      theChildNode,
                    *travNode;


    for (nodeCount=0; nodeCount<flatTree.lLength; nodeCount++) {
        theChildNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[nodeCount]);
        for (j=0; j<cBase; j++) {
            theChildNode->theProbs[j]=LIKELIHOOD_SCALER;
        }
        reexpnt = theChildNode->NeedToExponentiate(categID);
        if (reexpnt&&theChildNode->GetModelMatrix()) {
            theChildNode->RecomputeMatrix(categID, categoryCount);
        } else if (categID>=0) {
            theChildNode->SetCompMatrix (categID);
        }
    }

    for (nodeCount = 0; nodeCount<flatLeaves.lLength; nodeCount++) {
        theChildNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
        reexpnt = theChildNode->NeedToExponentiate(categID);
        if (reexpnt&&theChildNode->GetModelMatrix()) {
            theChildNode->RecomputeMatrix(categID, categoryCount);
        } else if (categID>=0) {
            theChildNode->SetCompMatrix (categID);
        }
    }

    for (nodeCount = 0; nodeCount<flatLeaves.lLength; nodeCount++) {
        travNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[nodeCount]);
        theChildNode = ((_CalcNode*)((BaseRef*)variablePtrs.lData)
                        [((node <long>*)(flatLeaves.lData[nodeCount]))->parent->in_object]);

        fastIndex = travNode->GetCompExp()->fastIndex();
        if (travNode->lastState>=0) {
            fastIndex+=travNode->lastState;
            _Parameter tmp = travNode->theProbs[travNode->lastState];
            for (long i=0; i<cBase; i++,fastIndex+=cBase) {
                theChildNode->theProbs[i]*=tmp* (*fastIndex);
            }
        } else {
            for (j=0; j<cBase; j++) {
                _Parameter tmp = *(travNode->theProbs) * *fastIndex;
                fastIndex++;
                for (long k = 1; k<cBase; k++,fastIndex++) {
                    tmp +=travNode->theProbs[k]* (*fastIndex);
                }
                theChildNode->theProbs[j]*=tmp;
            }
        }
    }

    for (nodeCount=0; nodeCount<flatTree.lLength; nodeCount++) {
        travNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[nodeCount]);
        nodeChild = ((node <long>*)(flatNodes.lData[nodeCount]))->parent;
        if (nodeChild) { // not a root yet
            theChildNode =
                ((_CalcNode*)((BaseRef*)variablePtrs.lData)[nodeChild->in_object]);

            fastIndex = travNode->GetCompExp()->fastIndex();

            for (j=0; j<cBase; j++) {
                _Parameter tmp = *(travNode->theProbs) * *fastIndex;
                fastIndex++;
                for (long k = 1; k<cBase; k++,fastIndex++) {
                    tmp +=travNode->theProbs[k]* (*fastIndex);
                }
                theChildNode->theProbs[j]*=tmp;
            }

        }
    }


    _Parameter result = 0;

    theChildNode = (_CalcNode*)(((BaseRef*)variablePtrs.lData)[theRoot->in_object]);
    //((_CalcNode*)(LocateVar(theRoot->get_data())));

    for (nodeCount=0; nodeCount<cBase; nodeCount++) {
        result+= theProbs[nodeCount]*theChildNode->theProbs[nodeCount];
    }

    if (result<=0.0) {
        return ALMOST_ZERO;
    }
    return result;
}

//_______________________________________________________________________________________________

long    _TheTree::ComputeReleafingCostChar (_DataSetFilter* dsf, long firstIndex, long secondIndex)
{

    char *pastState = dsf->GetColumn(firstIndex),
          *thisState = dsf->GetColumn(secondIndex);

    _SimpleList markedNodes (flatTree.lLength, 0, 0);

    for (long nodeID = 0; nodeID<flatLeaves.lLength; nodeID++) {
        long f = dsf->theNodeMap.lData[nodeID];
        if (thisState[f]!=pastState[f]) {
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

void    _TheTree::ClearConstraints (void)
{
    _CalcNode* travNode = StepWiseTraversal(true);
    while(travNode) {
        travNode->ClearConstraints();
        travNode = StepWiseTraversal();
    }
}


//_______________________________________________________________________________________________

long    _TheTree::ComputeReleafingCost (_DataSetFilter* dsf, long firstIndex, long secondIndex, _SimpleList* traversalTags, long orderIndex)
{

    //static _SimpleList flatLeaves, nodeCount;
    //static _List     flatTree;
    long        filterL = dsf->NumberDistinctSites();

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

void    _TheTree::MarkMatches (_DataSetFilter* dsf, long firstIndex, long secondIndex)
{

    long n = 0,f;

    _CalcNode* travNode ;

    for (n = 0; n<flatLeaves.lLength; n++) {
        travNode = (_CalcNode*)(((BaseRef*)flatCLeaves.lData)[n]);
        if (!dsf->CompareTwoSites(firstIndex,secondIndex,n)) {
            node <long>* theTreeNode = ((node <long>*)(flatLeaves.lData[n]))->parent;
            _CalcNode* cN = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[theTreeNode->in_object]);
            cN->cBase = -1;
        }
    }
    n = 0;
    for (f=0; f<flatTree.lLength; f++) {
        travNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[f]);
        if (travNode->cBase == -1) {
            node <long>* theTreeNode = ((node <long>*)(flatNodes.lData[f]))->parent;
            if (theTreeNode) {
                _CalcNode* cN = ((_CalcNode*)((BaseRef*)variablePtrs.lData)[theTreeNode->in_object]);
                cN->cBase = -1;
            }
        }
    }
    for (f=0; f<flatTree.lLength; f++) {
        travNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[f]);
        if (travNode->cBase != -1) {
            travNode->lastState = -2;
        } else {
            travNode->cBase = cBase;
        }
    }
}

//_______________________________________________________________________________________________

long    _TheTree::GetLowerBoundOnCost(_DataSetFilter* dsf)
{
    long n = 0,
         theCost = 0;

    _CalcNode* travNode ;

    for (long siteIndex = 0; siteIndex<dsf->theFrequencies.lLength; siteIndex++) {
        for (n = 0; n<flatTree.lLength; n++) {
            travNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[n]);
            travNode->lastState = -1;
        }
        for (long matchIndex = 0; matchIndex<dsf->theFrequencies.lLength; matchIndex++)
            if (matchIndex!=siteIndex) {
                MarkMatches (dsf,siteIndex,matchIndex);
            }
        for (n = 0; n<flatTree.lLength; n++) {
            travNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[n]);
            if (travNode->lastState != -2) {
                theCost += ((node <long>*)(flatNodes.lData[n]))->nodes.length;
            }
            travNode->lastState = -1;
        }
    }
    return theCost;
}

//_______________________________________________________________________________________________

long    _TheTree::GetLowerBoundOnCostWithOrder (_DataSetFilter* dsf, _SimpleList* sl)
{
    long n = 0,
         theCost = 0;

    _CalcNode* travNode ;

    for (long siteIndex = 0; siteIndex<dsf->theFrequencies.lLength; siteIndex++) {
        for (n = 0; n<flatTree.lLength; n++) {
            travNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[n]);
            travNode->lastState = -1;
        }
        for (long matchIndex = 0; matchIndex<siteIndex; matchIndex++)
            if (matchIndex!=siteIndex) {
                MarkMatches (dsf,sl->lData[siteIndex],sl->lData[matchIndex]);
            }

        for (n = 0; n<flatTree.lLength; n++) {
            travNode = (_CalcNode*)(((BaseRef*)flatTree.lData)[n]);
            if (travNode->lastState != -2) {
                theCost += ((node <long>*)(flatNodes.lData[n]))->nodes.length;
            }
            travNode->lastState = -1;
        }
    }
    return theCost;
}

//_______________________________________________________________________________________________

void    _TheTree::DumpingOrder (_DataSetFilter* dsf, _SimpleList& receptacle)
{

    _SimpleList flatLeaves, nodeCount;
    _List       flatTree;
    long        theCost,firstIndex,secondIndex;

    _CalcNode* travNode ;

    travNode = StepWiseTraversal (TRUE);

    while   (travNode) {
        travNode->GetProbs()[1]=1;
        flatTree<<travNode;
        nodeCount<<currentNode->get_num_nodes();
        travNode = StepWiseTraversal ();
        receptacle<<receptacle.lLength;
    }

    flatLeaves.Clear();
    travNode = LeafWiseTraversal (TRUE);

    while   (travNode) {
        flatLeaves<<long(currentNode);
        travNode = LeafWiseTraversal ();
    }

    for (firstIndex=0,secondIndex=1; secondIndex<dsf->NumberDistinctSites(); firstIndex=secondIndex++) {
        for (theCost = 0; theCost<flatLeaves.lLength; theCost++) {
            if ((*dsf)(firstIndex,theCost)!=(*dsf)(secondIndex,theCost)) {
                // mark all the nodes up the ladder as "tainted"
                node <long>* theTreeNode = (node <long>*)(flatLeaves(theCost));
                while (theTreeNode) {
                    ((_CalcNode*)(LocateVar(theTreeNode->get_data())))->SetSummedFlag();
                    theTreeNode=theTreeNode->get_parent();
                }
            }
        }

        // now compute the cost

        theCost = 0;

        for (long i=0; i<flatTree.lLength; i++) {
            travNode = (_CalcNode*)flatTree (i);
            if (travNode->IsSummedFlagged()) {
                travNode->RemoveSummedFlag();
                travNode->GetProbs()[1]++;
            }
        }
    }

    _SimpleList ref;
    for (theCost=0; theCost<flatTree.lLength; theCost++) {
        ref<<((_CalcNode*)flatTree(theCost))->GetProbs()[1];
    }
    SortLists (&ref, &receptacle);

}


//_______________________________________________________________________________________________

void    _TheTree::MolecularClock (_String& baseNode, _List& varsToConstrain)
{
    node<long>* topNode = nil;
    _CalcNode * curNode = StepWiseTraversal (true);
    if (!baseNode.Length()) { // called Molecular Clock on the entire tree
        topNode = &GetRoot();
        _String*  childNameP;
        if (rooted == ROOTED_LEFT) { // run separate constraint on the right child of the root
            childNameP = ((_CalcNode*)(((BaseRef*)variablePtrs.lData)[theRoot->go_down(theRoot->get_num_nodes())->in_object]))->GetName();
            _String childName = childNameP->Cut(childNameP->Find('.')+1,-1);
            MolecularClock (childName, varsToConstrain);
        } else if (rooted == ROOTED_RIGHT) {
            childNameP = ((_CalcNode*)(((BaseRef*)variablePtrs.lData)[theRoot->go_down(1)->in_object]))->GetName();
            _String childName = childNameP->Cut(childNameP->Find('.')+1,-1);
            MolecularClock (childName, varsToConstrain);
        }
    } else {
        baseNode = _String(".")&baseNode;
        while (curNode) {
            if (curNode->GetName()->endswith(baseNode)) {
                topNode = currentNode;
                break;
            }
            curNode = StepWiseTraversal();
        }
    }

    if (!topNode) {
        WarnError (_String ("Molecular clock constraint has failed, since node '")
                   &baseNode
                   &"' is not a part of tree '"
                   &*GetName() & "'");
    } else
        for (unsigned long k=1; k<varsToConstrain.lLength; k++) {
            long varIndex = LocateVarByName (*(_String*)varsToConstrain (k));
            if (varIndex<0) {
                WarnError (_String ("Molecular clock constraint has failed, since variable' ") &*(_String*)varsToConstrain (k) &"' is undefined.");
                return ;
            }
            curNode->RecurseMC (variableNames.GetXtra(varIndex), topNode, true, rooted);
        }
}



//_______________________________________________________________________________________________

node<long>* _CalcNode::LocateMeInTree (void)
{

    _String     baseNode = GetName()->Cut(0,GetName()->Find('.')-1);
    _TheTree    *parentTree = (_TheTree*)FetchVar(LocateVarByName(baseNode));
    _CalcNode   *curNode = parentTree->StepWiseTraversal (true);

    baseNode = GetName()->Cut(GetName()->FindBackwards('.',0,-1),-1);
    while (curNode) {
        if (curNode->GetName()->endswith(baseNode)) {
            return &parentTree->GetCurrentNode();
        }
        curNode = parentTree->StepWiseTraversal();
    }
    return nil;
}

//_______________________________________________________________________________________________

long _CalcNode::ConvertToSimpleMatrix (void)
{
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

    return 0;
}

//_______________________________________________________________________________________________

void _CalcNode::ConvertFromSimpleMatrix (void)
{
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

_Formula*   _CalcNode::RecurseMC (long varToConstrain, node<long>* whereAmI, bool first, char rooted)
{
    long descendants = whereAmI->get_num_nodes(),
         f = iVariables?iVariables->FindStepping(varToConstrain,2,1):-1,
         k,
         l,
         m,
         start = 0;

    if ((f<0)&&(!first)) { // bad bongos!
        _String errMsg ("Molecular clock constraint has failed, since variable ");
        errMsg = errMsg&*LocateVar(varToConstrain)->GetName();
        errMsg = errMsg&" is not an independent member of the node ";
        errMsg = errMsg&*GetName();
        WarnError (errMsg);
        return nil;
    }


    if (descendants == 0) // leaf node
        if (first) {
            return nil;
        } else {
            return new _Formula (LocateVar(iVariables->lData[f-1]),true);
        }

    if (first && (!whereAmI->get_parent()) && (rooted == ROOTED_LEFT)) {
        descendants --;
    }
    if (first && (!whereAmI->get_parent()) && (rooted == ROOTED_RIGHT)) {
        start++;
    }

    // internal node - must do some work

    _Formula**  nodeConditions = (_Formula**)MemAllocate((descendants-start)*sizeof(_Formula*));

    for (k=start+1; k<=descendants; k++) {
        node<long>* downWeGo = whereAmI->go_down(k);
        if (!(nodeConditions[k-1-start] = ((_CalcNode*)LocateVar(downWeGo->get_data()))->RecurseMC (varToConstrain, downWeGo))) {
            for (long f2 = 0; f2 < k-start-1; f2++) {
                delete nodeConditions[f2];
            }

            free (nodeConditions);
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
                LocateVar (((_Operation*)((*(nodeConditions[l])).GetList()(0)))->GetAVariable())->SetFormula (*nodeConditions[k]);
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
                    LocateVar (((_Operation*)((*(nodeConditions[l])).GetList()(0)))->GetAVariable())->SetFormula (*nodeConditions[k]);
                } else { // solve for a non-additive constraint
                    _Variable* nonAdd = LocateVar (((_Operation*)((*(nodeConditions[l])).GetList()(0)))->GetAVariable());
                    nodeConditions[l]->GetList().Delete(0);
                    _Formula  newConstraint;
                    newConstraint.Duplicate((BaseRef)nodeConditions[k]);
                    _String mns ("-");
                    _Operation mins (mns,2L);
                    for (m=0; m<nodeConditions[l]->GetList().lLength; m++) {
                        _Operation* curOp = (_Operation*)(*nodeConditions[l]).GetList()(m);
                        if (curOp->GetNoTerms()) {
                            newConstraint.GetList()&& &mins;
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

        _String      pls ('+');

        _Operation   *newVar = new _Operation,
        *plus   = new _Operation (pls, 2L);

        if (!(newVar && plus)) {
            checkPointer (nil);
        }

        newVar->SetAVariable (iVariables->lData[f-1]);

        result->GetList() << newVar;
        result->GetList() << plus;

        DeleteObject (newVar);
        DeleteObject (plus);

        free (nodeConditions);
        return result;
    }

    for (k=0; k<descendants-start; k++)
        if (nodeConditions[k]) {
            delete nodeConditions[k];
        }

    free (nodeConditions);
    return nil;
}

//_______________________________________________________________________________________________

void        _TheTree::ScanSubtreeVars  (_List& rec, char flags, _CalcNode* startAt)
// flags = 1 - do ind
// flags = 2 - do dep
{
    _SimpleList scanVars;
    _VariableContainer*  thisV;
    _String     chop;
    long        k,
                f;

    if (startAt) {
        thisV = startAt;
    } else {
        thisV = DepthWiseTraversal (true);
    }

    {
        _AVLList scanVarsA (&scanVars);
        if (flags&0x01) {
            thisV->ScanForVariables (scanVarsA,scanVarsA);
        }
        if (flags&0x02) {
            thisV->ScanForDVariables(scanVarsA,scanVarsA);
        }

        scanVarsA.ReorderList();
    }

    for (k=0; k<scanVars.lLength; k++) {
        thisV = (_VariableContainer*)LocateVar (scanVars.lData[k]);
        f = thisV->GetName()->FindBackwards('.',0,-1);
        if (f>=0) {
            chop = thisV->GetName()->Cut (f+1,-1);
            rec&& & chop;
        }
    }

    if (startAt) {
        thisV = StepWiseTraversalLevel(k,true);
        while (thisV&&(thisV!=startAt)) {
            thisV = StepWiseTraversalLevel(k);
        }
        if (thisV) {
            f = k;
            while ((k==f)&&thisV) {
                thisV = StepWiseTraversalLevel(k);
            }
            if (thisV) {
                while ((k>f)&&rec.lLength) {
                    thisV->MatchParametersToList(rec,true,(flags&0x02)!=0);
                    thisV = StepWiseTraversalLevel(k);
                }
                return;
            }
        }
        rec.Clear();
    } else {
        thisV = DepthWiseTraversal();
        while (thisV&&rec.lLength&&(!IsCurrentNodeTheRoot())) {
            thisV->MatchParametersToList(rec,true,(flags&0x02)!=0);
            thisV = DepthWiseTraversal();
        }
    }
}

//_______________________________________________________________________________________________

_String             _TreeTopology::CompareTrees      (_TreeTopology* compareTo)
{
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

            node<long>* meNode = DepthWiseStepTraverser (otherCT);

            while (meNode!=otherCT) {
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
                meNode = DepthWiseStepTraverser ((node<long>*)nil);
            }

            if (meNode!=otherCT) {
                meNode = DepthWiseStepTraverser (compareTo->theRoot);
                while (meNode!=theRoot) {
                    if (tCount==1) {
                        rerootAt = eqWithReroot & *LocateVar (meNode->in_object)->GetName() & '.';
                        break;
                    } else {
                        tCount --;
                    }

                    meNode = DepthWiseStepTraverser ((node<long>*)nil);
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

node<long>* _TreeTopology::prepTree4Comparison (_List& leafNames, _SimpleList& mapping, node<long>* topNode)
{
    node<long>* res     = topNode?topNode->duplicate_tree():theRoot->duplicate_tree(),
                * meNode;

    checkPointer (res);

    meNode =        DepthWiseStepTraverser (res);

    _SimpleList     indexer;

    while (meNode) {
        long        numChildren = meNode->get_num_nodes();

        _SimpleList * descendants = new _SimpleList;
        checkPointer (descendants);
        if (numChildren) {
            for (long k = 1; k <= numChildren; k++) {
                node <long> * aChild = meNode->go_down(k);
                (*descendants) << *(_SimpleList*)aChild->in_object;
            }
        } else {
            (*descendants) << leafNames.lLength;
            indexer        << leafNames.lLength;
            _String     * dd = (_String*)checkPointer(new _String);// (*LocateVar (meNode->in_object)->GetName(),treeNameLength,-1);
            GetNodeName (meNode, *dd);
            leafNames.AppendNewInstance (dd);
        }

        meNode->in_object = (long)descendants;
        meNode = DepthWiseStepTraverser ((node<long>*)nil);
    }

    // now sort leaf names and build the indexer

    mapping.Clear();
    mapping.Duplicate (&indexer);

    SortLists (&leafNames, &indexer);
    SortLists (&indexer, &mapping);

    return res;
}

//_______________________________________________________________________________________________

void    _TreeTopology::destroyCompTree (node<long>* compRoot)
{
    long    nc = compRoot->get_num_nodes();
    for (int i=1; i<=nc; i++) {
        destroyCompTree (compRoot->go_down(i));
    }

    DeleteObject ((BaseRef)compRoot->in_object);
    delete (compRoot);
}

//_______________________________________________________________________________________________


bool    _TheTree::MatchLeavesToDF (_SimpleList& tipMatches, _DataSetFilter* df, bool doNumeric)
{
    tipMatches.Clear();
    _CalcNode*      travNode = StepWiseTraversal(true);
    _List           tips;
    long            j,k;

    while (travNode) {
        if (IsCurrentNodeATip()) {
            _String tipName (travNode->GetName()->Cut(travNode->GetName()->FindBackwards('.',0,-1)+1,-1));
            tips&& &tipName;
        }
        travNode = StepWiseTraversal(false);
    }

    /*for (j=0;j<tips.lLength;j++)
    {
        k = df->FindSpeciesName((_String*)tips(j));
        if (k==-1) break;
        tipMatches<<k;
    }*/

    j = df->FindSpeciesName (tips, tipMatches);

    if (doNumeric) {
        if (j!=tips.lLength) {
            long sj = j;
            for (j=0; j<tips.lLength; j++) {
                _String *thisName = (_String*)tips(j);
                k = atoi (thisName->sData);
                _String tryAgain (k);
                if ((tryAgain.Equal(thisName))&&(k<=tips.lLength)) {
                    tipMatches<<k;
                } else {
                    break;
                }
            }
            if (j==tips.lLength) {
                if (tipMatches.Find(0)==-1) // map to indexing from 0
                    for (j=0; j<tips.lLength; j++) {
                        tipMatches.lData[j]--;
                    }
            } else {
                j=sj;
            }
        }
    }

    return (j == tips.lLength);
}
//_______________________________________________________________________________________________

_List*  _TreeTopology::SplitTreeIntoClustersInt (node<long>* , _List* , _AVLListX& , long , long )
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

_List*  _TreeTopology::SplitTreeIntoClusters (unsigned long size, unsigned long tol)
// assume fixed rooting for now
// returns the list of list speccing pointers to calcnodes rooting the splits
// the member list contains 2 or more entries for each node:
// itself, the number of (independent)leaves encompassed by that node,
// and an optional references to the node whose
// results it depends on.

// assumes that size-tol>=2

{
    _SimpleList     counts;
    _AVLListX       cavl (&counts);

    DepthWiseT      (true);

    while (currentNode) {
        long nC = currentNode->get_num_nodes();
        if (nC) {
            long c = 0;
            for (long k=1; k<=nC; k++) {
                c += counts.lData[currentNode->go_down(k)->in_object];
            }
            cavl.Insert((BaseRef)currentNode->in_object,c);
        } else {
            cavl.Insert((BaseRef)currentNode->in_object,1);
        }

        DepthWiseT (false);
    }

    _List*                    result = new _List;
    checkPointer             (result);

    DeleteObject(SplitTreeIntoClustersInt (theRoot, result, cavl, size, tol));

    return           result;
}

//_______________________________________________________________________________________________

_String _TreeTopology::MatchTreePattern (_TreeTopology* compareTo)
// the pattern is this
// compare to is the potential match
{
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
            node<long> * padresStillSuck = DepthWiseStepTraverser (myCT);
            while (padresStillSuck != myCT) {
                if (padresStillSuck->get_num_nodes() == 0) {
                    _SimpleList*  descendants = ((_SimpleList*)padresStillSuck->in_object);

                    if ((descendants&&(allowedLeaves.lData[descendants->lData[0]] == 0))||(!descendants)) {
                        // mark this node for deletion
                        if (descendants) {
                            DeleteObject(descendants);
                        }

                        node<long>* sacLamb = padresStillSuck;
                        padresStillSuck = DepthWiseStepTraverser ((node<long>*)nil);

                        if (sacLamb->parent->get_num_nodes()==1) {
                            DeleteObject((BaseRef)sacLamb->parent->in_object);
                            sacLamb->parent->in_object = nil;
                        }

                        sacLamb->parent->detach_child (sacLamb->get_child_num());
                        delete (sacLamb);
                        continue;
                    }
                }
                padresStillSuck = DepthWiseStepTraverser ((node<long>*)nil);
            }

            // pass 2 - prune internal nodes with exactly one child

            // >O
            padresStillSuck = DepthWiseStepTraverser (myCT);

            while (padresStillSuck) {
                long nn = padresStillSuck->get_num_nodes();
                if (nn == 1) {
                    node<long>* loneChild = padresStillSuck->go_down(1);
                    DeleteObject((_SimpleList*)padresStillSuck->in_object);
                    padresStillSuck->in_object = loneChild->in_object;
                    padresStillSuck->detach_child (1);
                    delete (loneChild);
                } else {
                    if (nn > 1) {
                        _SimpleList * myDescs = (_SimpleList*)padresStillSuck->in_object;
                        myDescs->Clear();

                        for (long cc = 1; cc <= nn; cc++) {
                            _SimpleList temp;

                            temp.Union (*myDescs, *(_SimpleList*)padresStillSuck->go_down(cc)->in_object);

                            myDescs->Clear();
                            myDescs->Duplicate (&temp);
                        }

                        /*for (long cc2 = 0; cc2 < myDescs->lLength; cc2++)
                            myDescs->lData[cc2] = recordTransfer.lData[myDescs->lData[cc2]];*/
                    } else {
                        ((_SimpleList*)padresStillSuck->in_object)->lData[0] = recordTransfer.lData[((_SimpleList*)padresStillSuck->in_object)->lData[0]];
                    }

                }
                padresStillSuck = DepthWiseStepTraverser ((node<long>*)nil);
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

            node<long>* meNode = DepthWiseStepTraverser (otherCT);

            while (meNode!=otherCT) {
                if (meNode->get_num_nodes()) {
                    compRes = internalTreeCompare (myCT, meNode, reindexer, 1, myLeaves.lLength, nil, compareTo, true);
                    if (compRes>0) {
                        break;
                    } else if (compRes) {
                        meNode = otherCT;
                        break;
                    }
                }

                tCount ++;
                meNode = DepthWiseStepTraverser ((node<long>*)nil);
            }

            if (meNode!=otherCT) {
                meNode = DepthWiseStepTraverser (compareTo->theRoot);
                while (meNode!=theRoot) {
                    if (tCount==1) {
                        rerootAt = _String("Equal with reroot at ") & *LocateVar (meNode->in_object)->GetName() & '.';
                        break;
                    } else {
                        tCount --;
                    }

                    meNode = DepthWiseStepTraverser ((node<long>*)nil);
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

void     _TheTree::AddNodeNamesToDS (_DataSet* ds, bool doTips, bool doInternals, char dOrS)
{
    if (dOrS == 2 && doTips && doInternals) {
        AddNodeNamesToDS (ds, false, true, 0);
        AddNodeNamesToDS (ds, true, false, 0);
        return;
    }

    _CalcNode*iNodeTraverser = dOrS?DepthWiseTraversal (true) :StepWiseTraversal (true);

    long j = GetName()->sLength+1;


    while (iNodeTraverser) {
        if (IsCurrentNodeATip()) {
            if (doTips) {
                ds->GetNames().AppendNewInstance(new _String (*iNodeTraverser->GetName(), j, -1));
            }
        } else if (doInternals) {
            ds->GetNames().AppendNewInstance(new _String (*iNodeTraverser->GetName(), j, -1));
        }

        iNodeTraverser = dOrS?DepthWiseTraversal (false) :StepWiseTraversal (false);
    }
}
