/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon42@uwo.ca)
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

#include  <math.h>

#include  "category.h"
#include  "function_templates.h"
#include  "global_things.h"
#include  "global_object_lists.h"

using namespace hy_global;
using namespace hyphy_global_objects;


//___________________________________________________________________________________________

const _String    defaultEqual         ("EQUAL"),
           medianRep          ("MEDIAN"),
           scaledMedianRep        ("SCALED_MEDIAN"),
           maxCatIvals            ("MAX_CATEGORY_INTERVALS"),
           constantOnPartition    ("CONSTANT_ON_PARTITION");

unsigned long maxCategoryIntervals = 100UL;

#ifdef     _SLKP_LFENGINE_REWRITE_
#define    SLIGHT_SHIFT  0.
#else
#define    SLIGHT_SHIFT  1e-150
#endif

extern     _List        modelNames;
extern     _SimpleList  modelMatrixIndices,
           modelFrequenciesIndices;


//___________________________________________________________________________________________

_CategoryVariable::_CategoryVariable (_String& name, _List* parms, _VariableContainer* theP):_Variable (name)
{
    values = intervalEnds = weights = conditionalWeights = nil;
    hiddenMarkovModel = HY_NO_MODEL;
    flags             = 0;
    covariant        = -1;
    intervalSplitter = -1;
    Construct (*parms, theP);
}

//___________________________________________________________________________________________
bool _CategoryVariable::checkWeightMatrix(_Matrix& w, long row) {
    bool    check = true;
    _Constant iterate;
    hyFloat sumCheck = 0;
    if (row>=0) {
        long shift = w.GetVDim()*row;
        for (long i=0; i<intervals; i++) {
            sumCheck+=w[shift+i];
        }
        if (fabs(sumCheck-1.0)>=1e-8) {
            if (sumCheck<1.0) {
                sumCheck = 1.0/sumCheck;
            }
            for (long k=0; k<intervals; k++) {
                w[shift+k]/=sumCheck;
            }
        }
    } else {
        for (long i=0L; i<intervals; i++) {
            sumCheck+=w.theData[i];
        }
        if (fabs(sumCheck-1.0)>=1e-8) {
            if (sumCheck<1.0) {
                sumCheck = 1.0/sumCheck;
            }
            for (long k=0; k<intervals; k++) {
                w.theData[k]/=sumCheck;
            }
        }
    }
    return check;
}

//___________________________________________________________________________________________

void _CategoryVariable::Construct (_List& parameters, _VariableContainer *theP) {
// main workhorse constructor
// expects a list of strings containing the following
// number of intervals
// types of representation (MEAN or MEDIAN)
// probability weight formula or "EQUAL"
// probability density function (must contain at least hy_x_variable)
// (optional) cumulative distribution f-n (must contain at least hy_x_variable)
// range for hy_x_variable


    _String     errorMsg = _String ("While attempting to construct category variable ") & *GetName() & ": ";

    _SimpleList scannedVarsList,
                variableDependanceAllocationsAux;

    _AVLListXL  variableDependanceAllocations (&variableDependanceAllocationsAux);

    bool    check,
            covariantVar = false;

    long    f;

    Clear(); // clear this variable if needed


    checkParameter (maxCatIvals, maxCategoryIntervals, 100UL);
    // set up the number of intervals and the matrices
    _String*            param = (_String*)parameters(0);
    intervals                 = ProcessNumericArgument(param,theP);
    if (intervals<=0) {
        HandleApplicationError (errorMsg & _String("Category variable must have a positive number of classes - had ")
                   & *param);
        return;
    }

    if (intervals>maxCategoryIntervals) {
        intervals = maxCategoryIntervals;
        ReportWarning (errorMsg & _String("Category variable cannot have more than ")&maxCatIvals&" classes - had "
                       & *param&". Reset to "& _String(intervals));
    }

    values             = new _Matrix (intervals, 1, false, true);
    intervalEnds       = new _Matrix (intervals, 1, false, true);
    weights            = new _Matrix (intervals, 1, false, true);

    // construct the formula for interval weights
    param = (_String*)parameters (1);
    check = false;
    if (*param != defaultEqual) {
        // do something here, otherwise they are just equal
        _String             splitterName (AppendContainerName(*param,theP));
        f = LocateVarByName (splitterName);

        if (f>=0 && FetchVar(f)->IsCategory()) {
            _CategoryVariable * iSplitter = (_CategoryVariable*)FetchVar(f);
            if (!CheckEqual (iSplitter->GetMinX(),SLIGHT_SHIFT) ||
                    !CheckEqual (iSplitter->GetMaxX(),1.0) ||
                    *theName == splitterName ||
                    (intervals = iSplitter->GetNumberOfIntervals()+1) < 2) {
                HandleApplicationError (errorMsg & _String("Category variables which specify interval splitting options must be supported on [0,1], and not result in circular dependance"));
                return;
            }

            intervalSplitter = iSplitter->get_index();

            _AVLList      ivl (&scannedVarsList);
            iSplitter->ScanForVariables (ivl, true);
            ivl.ReorderList();
            
            BatchDelete(values, intervalEnds, weights);

            values             = new _Matrix (intervals, 1, false, true);
            intervalEnds       = new _Matrix (intervals, 1, false, true);
            weights            = new _Matrix (intervals, 1, false, true);
            check = true;
        } else {
            _Formula      probabilities(*param,theP);
            {
                _AVLList      sv (&scannedVarsList);
                probabilities.ScanFForVariables (sv, true);
                sv.ReorderList();
            }

            _Matrix *tryMatrix = (_Matrix*)probabilities.GetTheMatrix();
            // check to see if it is a matrix spec
            if (tryMatrix) {
                _Matrix* weightMatrix = (_Matrix*)tryMatrix;

                if (!( ((weightMatrix->GetHDim()==1)&&(weightMatrix->GetVDim()==intervals))||
                        ((weightMatrix->GetHDim()==intervals)&&(weightMatrix->GetVDim()==1)))) {
                    if (weightMatrix->GetVDim() == intervals)
                        // covariance structure here
                    {
                        if (weightMatrix->GetHDim() > 1) {
                            check        = true;
                            covariantVar = true;
                            if (weightMatrix->IsIndependent())
                                for (long k=0; check && (k<weightMatrix->GetHDim()); k++) {
                                    check = check && checkWeightMatrix (*weightMatrix, k);
                                }

                            if (check) {
                                weights->Duplicate(weightMatrix);
                            }
                        } else {
                            check = false;
                        }
                    } else {
                        check = false;
                    }
                } else
                    // indepenent category variable
                {
                    if (weightMatrix->IsIndependent()) {
                        check=checkWeightMatrix(*weightMatrix);
                    } else {
                        for (long k = 0; k < weightMatrix->GetHDim() * weightMatrix->GetVDim (); k++) {
                            _Formula* thisCell = weightMatrix->GetFormula (k,-1);
                            if (thisCell) {
                                _SimpleList   probVars;
                                _AVLList      sv (&probVars);
                                thisCell->ScanFForVariables (sv, true);
                                sv.ReorderList();
                                for (long v = 0; v < probVars.lLength; v++) {
                                    long f = variableDependanceAllocations.Find ((BaseRef)probVars.list_data[v]);
                                    if (f < 0) {
                                        f = variableDependanceAllocations.Insert ((BaseRef)probVars.list_data[v], (long)(new _SimpleList (intervals,0,0)),false);
                                    }
                                    ((_SimpleList*) variableDependanceAllocations.GetXtra (f))->list_data[k] = 1;
                                }
                            }
                        }
                        check = true;
                    }
                    if (check) {
                        weights->Duplicate(weightMatrix);
                    }
                }
            } else {
                if (scannedVarsList.lLength) {
                    if(scannedVarsList.lLength==1) {
                        if (scannedVarsList[0]==hy_n_variable->get_index()) {
                              for (unsigned long i=0UL; i<intervals; i++) {
                                hy_n_variable->SetValue(new _Constant ((hyFloat)i), false, true, NULL);
                                (*weights)[i]= probabilities.Compute()->Value();
                            }
                            check = checkWeightMatrix (*weights);
                        }
                    }
                    if (!check) {
                        ReportWarning (errorMsg & _String("Interval weights must be specified in terms of ") & hy_n_variable->GetName()->Enquote());
                    }
                    scannedVarsList.Clear();
                }
            }
        }
    } else {
        for (long i=0; i<intervals; i++) {
            (*weights)[i] = 1.0/intervals;
        }
        check = true;
    }

    if (!check) {
        ReportWarning (errorMsg & _String("Invalid weight matrix for a category variable. Reset to EQUAL"));
        covariantVar = false;
        for (long i=0; i<intervals; i++) {
            (*weights)[i] = 1.0/intervals;
        }
    }

    // set the representation mode

    param = (_String*)parameters(2);
    if (medianRep == *param) {
        representation = MEDIAN;
    } else if (scaledMedianRep == *param) {
        representation = SCALED_MEDIAN;
    } else {
        representation = MEAN;
    }

    // get the probability density

    param = (_String*)parameters(3);

    if (!covariantVar) {
        _FormulaParsingContext fpc (nil, theP);
        Parse (&density, *param, fpc, nil); // check if the formula is good
    }

    if (!density.IsEmpty()) {
        if (covariantVar) {
            HandleApplicationError (errorMsg & "Continuous distributions are not supported by non-independent category variables - specify a discrete range.");
            return;
        }

        {
            _SimpleList   densityVars,
                          existingVars (scannedVarsList);

            _AVLList      sv (&densityVars);
            density.ScanFForVariables (sv, true);
            sv.ReorderList();
            scannedVarsList.Union (densityVars,existingVars);
        }

        f = scannedVarsList.Find(hy_x_variable->get_index());
        if (f!=-1) { // no dummy variable
            check = true;
            scannedVarsList.Delete(f);
        } else {
            //scannedVars.Clear();
            check = false;
        }

        // get the cumulative probability

        param = (_String*)parameters(4);
        if (param->empty ()) { // no cumul. dist'n specified - integration is yet to be implemented
            ReportWarning (errorMsg & _String("Runtime integration of probability density can be _very_ slow - please provide the analytic form for cumulative distribution if known."));
        } else {
            if(check) {
                _FormulaParsingContext fpc (nil, theP);
                Parse(&cumulative,*param,fpc, nil);
                {
                    _SimpleList   densityVars,
                                  existingVars (scannedVarsList);

                    _AVLList      sv (&densityVars);
                    cumulative.ScanFForVariables (sv, true);
                    sv.ReorderList();
                    scannedVarsList.Union (densityVars,existingVars);
                }
                f = scannedVarsList.BinaryFind(hy_x_variable->get_index());
                if (f<0) { // no dummy variable
                    HandleApplicationError (errorMsg & _String("Cumulative distribution must be specified in terms of ") & hy_x_variable->GetName()->Enquote() & ". Had: " &param->Enquote());
                    return;
                } else {
                    scannedVarsList.Delete(f);
                }
            }

        }

        // get the bounds here
        param = (_String*)parameters(5);
        x_min = param->to_float()+SLIGHT_SHIFT;

        param = (_String*)parameters(6);
        x_max = param->to_float()-SLIGHT_SHIFT;

        if (x_max<=x_min) {
            HandleApplicationError (errorMsg & _String("Bad variable bounds. Had:")&*(_String*)parameters(5)&" and "&*param);
            return;
        }

        if (x_max >= INFINITE_BOUND) {
            x_max = INFINITE_BOUND;
        }

        if (!check) { // uniform distribution
            if (x_max==INFINITE_BOUND) {
                HandleApplicationError ( errorMsg & _String("Can't have uniform distributions over infinite intervals. "));
                return;
            } else {
                errorMsg = errorMsg & _String("Since density ")&*(_String*)parameters(3)& " contains no " & kXVariableName.Enquote() &", the distribution was set to uniform over ["&_String(x_min)&","&_String(x_max)&"]";
                ReportWarning (errorMsg);
                density.Clear();
                hyFloat dns = 1.0/(x_max-x_min);
                errorMsg = _String(dns);
                _FormulaParsingContext fpc;
                Parse(&density, errorMsg,fpc, nil);
                errorMsg = _String(dns)&"*(" & kXVariableName &"-"&_String(x_min)&")";
                Parse(&cumulative, errorMsg,fpc, nil);
            }
        }
    } else { // enumerated interval parameters
        x_min = 0.;
        x_max = 1.;
        // expect a matrix for the cumulative distribution
        if (covariantVar) {
            param               = (_String*)parameters(3);
            _String             splitterName (AppendContainerName(*param,theP));
            f = LocateVarByName (splitterName);
            if (f>=0) {
                _Variable * cbase = FetchVar (f);
                if (cbase->IsCategory()) {
                    check = ((_CategoryVariable*)cbase)->IsUncorrelated();
                } else {
                    check = false;
                }
            } else {
                check = false;
            }

            if (!check) {
                HandleApplicationError ( errorMsg & *param & " must be the identifier of an existing, independent category variable.");
                return;
            } else {
                covariant = variableNames.GetXtra (f);
                if (((_CategoryVariable*)FetchVar (f))->GetNumberOfIntervals() != weights->GetHDim()) {
                    HandleApplicationError (errorMsg & *param & " is incompatible with the conditional probability matrix for " & *GetName() &". The number of possible values of " & *param &" must match the row count of the matrix.");
                    return;
                }
            }
        }
        param = (_String*)parameters(4);
        _Formula cumulative(*param,theP);
        {
            _SimpleList   densityVars,
                          existingVars (scannedVarsList);

            _AVLList      sv (&densityVars);
            cumulative.ScanFForVariables (sv, true);
            sv.ReorderList();
            scannedVarsList.Union (densityVars,existingVars);
        }
        // check to see if it is a matrix spec
        HBLObjectRef  tryMatrix = cumulative.GetTheMatrix();
        if (tryMatrix) {
            _Matrix* catMatrix = (_Matrix*)tryMatrix;
            if (!( ((catMatrix->GetHDim()==1)&&(catMatrix->GetVDim()==intervals))||
                    ((catMatrix->GetHDim()==intervals)&&(catMatrix->GetVDim()==1)))) {
                HandleApplicationError (errorMsg & ("Dimension of category representatives matrix is not the same as the number of categories"));
                return;
            } else {
                values->Duplicate(catMatrix);
                if (!catMatrix->IsIndependent()) { // not a constant matrix
                    for (long k = 0; k < catMatrix->GetHDim() * catMatrix->GetVDim (); k++) {
                        _Formula* thisCell = catMatrix->GetFormula (k,-1);
                        if (thisCell) {
                            _SimpleList   densityVars,
                                          existingVars (scannedVarsList);

                            _AVLList      sv (&densityVars);
                            thisCell->ScanFForVariables (sv, true);
                            sv.ReorderList();
                            for (long v = 0; v < densityVars.lLength; v++) {
                                long f = variableDependanceAllocations.Find ((BaseRef)densityVars.list_data[v]);
                                if (f < 0) {
                                    f = variableDependanceAllocations.Insert ((BaseRef)densityVars.list_data[v], (long)(new _SimpleList (intervals,0,0)),false);
                                }
                                ((_SimpleList*) variableDependanceAllocations.GetXtra (f))->list_data[k] = 1;
                            }

                            scannedVarsList.Union (densityVars,existingVars);
                        }
                    }
                }
            }
        } else {
            HandleApplicationError (errorMsg & ("Expected an explicit enumeration of category representatives in place of cumulative distribution: ") & ((_String*)cumulative.toStr(kFormulaStringConversionNormal))->Enquote() );
            return;
        }
    }

    // disallow category -> category dependance
    for (long i=0; i<scannedVarsList.lLength; i++) {
        _Variable * curVar = (_Variable*)variablePtrs (scannedVarsList.list_data[i]);
        if (curVar->IsCategory()) {
            HandleApplicationError (errorMsg & _String("Can't create a category variable that depends on another category variable, ") & curVar->GetName()->Enquote());
            return;
        }
    }



    hiddenMarkovModel = HY_NO_MODEL;

    parameterList.Duplicate  (&scannedVarsList);
    // finally go thru all the variables and put them where they belong in dependance containers

    _SimpleList     exclude;

    if (parameters.countitems()>7) { // aux mean formula
        param = (_String*)parameters(7);
        _FormulaParsingContext fpc (nil, theP);
        Parse    (&meanC,*param,fpc, nil);

        if (parameters.lLength>8) {
            if (((_String*)parameters(8))->nonempty()) {
                _String hmmModelName = AppendContainerName(*(_String*)parameters(8),theP);
                f = FindModelName(hmmModelName);
                if (f==-1) {
                    if (constantOnPartition == *(_String*)parameters (8)) {
                        flags = CONSTANT_ON_PARTITION;
                    } else {
                        HandleApplicationError (errorMsg & (*(_String*)parameters(8))& " is not an existing model identifier in call to 'category'");
                        return;
                    }
                } else {
                    if (covariantVar) {
                        HandleApplicationError (errorMsg & "Non-independent random variables can't also be hidden Markov.");
                        return;
                    }
                    long mindex = f;
                    _Matrix * hmm,
                            * freq;

                    bool    mbf;

                    RetrieveModelComponents (mindex, hmm, freq, mbf);
                    mbf = false;

                    if (hmm) {

                        f = weights->GetHDim()*weights->GetVDim();

                        if (hmm->GetHDim()==f && hmm->GetVDim()==f) {
                            _SimpleList   hmmVars,
                                          existingVars (parameterList);

                            _AVLList      sv (&hmmVars);

                            hmm->ScanForVariables (sv,true);
                            freq->ScanForVariables (sv,true);


                            sv.ReorderList();
                            parameterList.Union (hmmVars,existingVars);
                            exclude.Subtract (hmmVars,existingVars);

                            hiddenMarkovModel = mindex;
                            mbf = true;
                        }
                    }

                    if (!mbf) {
                        HandleApplicationError (errorMsg & (*(_String*)parameters(8))& " is not a valid HMM-component model (square matrix of dimension "&_String (f) & ") identifier in call to 'category'");
                    }
                }
            }
        }
    }

    for (long vid = 0; vid < parameterList.lLength; vid ++) {
        long vf = variableDependanceAllocations.Find ((BaseRef)parameterList.list_data[vid]);
        if (vf >= 0) {
            affectedClasses << (_SimpleList*)(variableDependanceAllocations.GetXtra (vf));
        } else if (exclude.Find (parameterList.list_data[vid]) >= 0) {
            affectedClasses.AppendNewInstance (new _SimpleList (intervals,0,0));
        } else {
            affectedClasses.AppendNewInstance (new _SimpleList (intervals,1,0));
        }

        _String vlog = _String ("Variable ") & *LocateVar(parameterList.list_data[vid])->GetName() & " mapped to class " & _String(((_String*)affectedClasses(vid)->toStr()));
        ReportWarning (vlog);
    }


    if (covariant >= 0) {
        conditionalWeights = new _Matrix (intervals, 1, false, true);
    }

    /*{
        _SimpleList tl;
        _AVLList    test (&tl);
        ScanForVariables (test, true);

        _String str (128L, true);

        _SimpleList  hist;
        long         ls, cn;

        cn = test.Traverser (hist,ls,test.GetRoot());

        while (cn>=0)
        {
            long keyVal = (long)test.Retrieve (cn);
            str << *LocateVar(keyVal)->GetName();
            if (LocateVar(keyVal)->IsGlobal())
                str << "  Global";
            str << '\n';
            cn = test.Traverser (hist,ls);
        }

        StringToConsole (str);
    }*/
}

//___________________________________________________________________________________________

bool _CategoryVariable::IsUncorrelated (void)
{
    return covariant == -1;
}

//___________________________________________________________________________________________

bool _CategoryVariable::IsLayered (void)
{
    return intervalSplitter >= 0;
}

//___________________________________________________________________________________________

void _CategoryVariable::ChangeNumberOfIntervals (long newi) {
    if (newi==intervals) {
        return;
    }

    DeleteObject (values);
    DeleteObject (intervalEnds);
    DeleteObject (weights);
    intervals = newi;
    values = new _Matrix (intervals, 1, false, true);
    intervalEnds = new _Matrix (intervals, 1, false, true);
    weights = new _Matrix (intervals, 1, false, true);
    covariant = -1;
    intervalSplitter = -1;

    for (long i=0; i<intervals; i++) {
        (*weights)[i] = 1.0/intervals;
    }

    UpdateIntervalsAndValues();
}

//___________________________________________________________________________________________

BaseRef     _CategoryVariable::makeDynamic(void) const {
    _CategoryVariable* result = new _CategoryVariable();
    result->Duplicate(this);
    return result;
}
//___________________________________________________________________________________________
void        _CategoryVariable::Duplicate(BaseRefConst s) {
    _CategoryVariable const* cv = (_CategoryVariable const*)s;
    Clear();
    intervals = cv->intervals;
    density.Duplicate (&cv->density);
    cumulative.Duplicate (&cv->cumulative);
    meanC.Duplicate (&cv->meanC);
    representation = cv->representation;
    x_min = cv->x_min;
    x_max = cv->x_max;
    if (cv->values) {
        values = (_Matrix*)cv->values->makeDynamic();
    } else {
        values = nil;
    }
    if (cv->intervalEnds) {
        intervalEnds = (_Matrix*)cv->intervalEnds->makeDynamic();
    } else {
        intervalEnds = nil;
    }
    if (cv->weights) {
        weights = (_Matrix*)cv->weights->makeDynamic();
    } else {
        weights = nil;
    }

    if (cv->conditionalWeights) {
        conditionalWeights = (_Matrix*)cv->conditionalWeights->makeDynamic();
    } else {
        conditionalWeights = nil;
    }

    covariant = cv->covariant;
    intervalSplitter = cv->intervalSplitter;
    hiddenMarkovModel = cv->hiddenMarkovModel;
    flags = cv->flags;
    parameterList.Duplicate (&cv->parameterList);
    affectedClasses.Duplicate (&cv->affectedClasses);
    this->_Variable::Duplicate (s);
}

//___________________________________________________________________________________________
void    _CategoryVariable::Clear (void)
{
    density.Clear       ();
    cumulative.Clear    ();
    DeleteObject        (values);
    DeleteObject        (intervalEnds);
    DeleteObject        (weights);
    DeleteObject        (conditionalWeights);
    covariant         = -1;
    intervalSplitter  = -1;
    hiddenMarkovModel = HY_NO_MODEL;
    flags             = 0;
    parameterList.Clear();
    affectedClasses.Clear();
}

//___________________________________________________________________________________________
BaseRef _CategoryVariable::toStr (unsigned long)
{
    UpdateIntervalsAndValues(true);
    _StringBuffer * result  = new _StringBuffer (256UL);
    if (weights) {
        *result<< "\nClass weights are:";
        _Matrix* cw =(_Matrix*)weights->ComputeNumeric();
        checkWeightMatrix(*cw);
        result->AppendNewInstance((_String*)cw->toStr());
        *result<<'\n';
    }
    if (values) {
        *result<<"Classes represented by:";
        result->AppendNewInstance((_String*)values->toStr());
    }
    if (intervalEnds) {
        *result<< "Interval ends:";
        result->AppendNewInstance((_String*)intervalEnds->toStr());
    }
    if (!density.IsEmpty()) {
        *result << "\nSupported on [";
        *result << _String(x_min);
        *result << ',';
        *result << _String(x_max);
        *result << "]\n";
    }
    result->TrimSpace ();
    return result;
}

//___________________________________________________________________________________________
hyFloat  _CategoryVariable::SetIntervalValue (long ival, bool recalc)
{
    hyFloat newIntervalValue;
    if (recalc) {
        newIntervalValue  = GetValues()->theData[ival];
    } else {
        newIntervalValue = ((_Matrix*)values->RetrieveNumeric())->theData[ival];
    }
    SetValue (new _Constant(newIntervalValue),false,true,NULL);
    /*if (ival == 0) {
        printf ("\n\n");
    }
    printf ("%ld => %g\n", ival, newIntervalValue);*/
    return newIntervalValue;
}

//___________________________________________________________________________________________
hyFloat  _CategoryVariable::GetIntervalValue (long ival)
{
    if (values) {
        return GetValues()->theData[ival];
    } else {
        return 0.0;
    }
}

//___________________________________________________________________________________________
hyFloat  _CategoryVariable::GetIntervalWeight (long ival)
{
    if (weights) {
        if (covariant >= 0 || intervalSplitter >= 0) {
            return GetWeights()->theData[ival];
        }

        if (weights->IsIndependent()) {
            return ((_Matrix*)weights->ComputeNumeric())->theData[ival];
        } else {
            _Matrix* cw = ((_Matrix*)weights->ComputeNumeric());
            checkWeightMatrix(*cw);
            return cw->theData[ival];
        }
    } else {
        return 0.0;
    }
}

//___________________________________________________________________________________________
_Matrix*    _CategoryVariable::GetValues (void)
{
    return (_Matrix*)values->ComputeNumeric();
}

//___________________________________________________________________________________________
long        _CategoryVariable::GetCurrentState (void)
{
    _Matrix         *v  = GetValues();
    hyFloat      cv  = Compute()->Value();

    for (long res = 0; res < intervals; res ++)
        if (CheckEqual (cv, v->theData[res])) {
            return res;
        }

    return 0;
}


//___________________________________________________________________________________________
_Matrix*    _CategoryVariable::GetWeights (bool covAll)
{
    _Matrix * cw;

    if (intervalSplitter>=0) {
        _CategoryVariable * iSplitter = (_CategoryVariable*)LocateVar (intervalSplitter);
        cw = iSplitter->GetValues();
        hyFloat      minusMe = 0.0;
        for (long k=0; k<intervals-1; k++) {
            weights->theData[k] = cw->theData[k] - minusMe;
            minusMe = cw->theData[k];
        }
        weights->theData[intervals-1] = 1.-minusMe;
        return weights;
    }


    if (weights->IsIndependent()) {
        cw = (_Matrix*)weights->ComputeNumeric();
    } else {
        cw = ((_Matrix*)weights->ComputeNumeric());
        if (covariant < 0) {
            checkWeightMatrix(*cw);
        }
    }

    if (covariant >= 0) {
        _CategoryVariable * cv = (_CategoryVariable*)LocateVar (covariant);
        if (covAll) {
            long iv2 = cv->GetNumberOfIntervals();

            for (long m=0; m<iv2; m++) {
                checkWeightMatrix (*cw, m);
            }

            _Matrix * cw2 = cv->GetWeights ();

            for (long k=0; k<intervals; k++) {
                hyFloat sum = 0.0;
                for (long j=0; j<iv2; j++) {
                    sum += cw2->theData[j]* (*cw)(j,k);
                }

                conditionalWeights->theData[k] = sum;
            }
            cw = conditionalWeights;
        } else {
            long        rowIdx = cv->GetCurrentState()*cw->GetVDim();
            for (long k=0; k<intervals; k++) {
                conditionalWeights->theData[k] = cw->theData[rowIdx+k];
            }
            cw = conditionalWeights;
            checkWeightMatrix (*cw);
        }
    }

    return cw;
}

//___________________________________________________________________________________________
_Matrix*    _CategoryVariable::GetIntervalEnds (void)
{
    return (_Matrix*)intervalEnds->ComputeNumeric();
}

//___________________________________________________________________________________________
_Matrix*    _CategoryVariable::ComputeHiddenMarkov (void) {
    _Matrix *hmmr = (_Matrix*)((_Matrix*)LocateVar (modelMatrixIndices.list_data[hiddenMarkovModel])->GetValue())->ComputeNumeric();
    if (!hmmr->IsValidTransitionMatrix()) {
        HandleApplicationError(_String ("Hidden Markov Model transition matrix variable did not compute to a valid transition matrix"));
    }
    return hmmr;
}

//___________________________________________________________________________________________
_Matrix*    _CategoryVariable::ComputeHiddenMarkovFreqs (void)
{
    long       fIndex = modelFrequenciesIndices.list_data[hiddenMarkovModel];
    if (fIndex<0) {
        fIndex = -fIndex-1;
    }
    _Matrix * hmmf = (_Matrix*)((_Matrix*)LocateVar (fIndex)->GetValue())->ComputeNumeric();
    checkWeightMatrix(*hmmf);
    return hmmf;
}

//___________________________________________________________________________________________
_Matrix*    _CategoryVariable::GetHiddenMarkov (void) const {
    _Variable* theMX = LocateVar (modelMatrixIndices.list_data[hiddenMarkovModel]);
    return (_Matrix*)theMX->GetValue();
}

//___________________________________________________________________________________________
_Matrix*    _CategoryVariable::GetHiddenMarkovFreqs (void) const {
    long       fIndex = modelFrequenciesIndices.list_data[hiddenMarkovModel];
    if (fIndex<0) {
        fIndex = -fIndex-1;
    }
    _Variable* theMX = LocateVar (fIndex);
    return (_Matrix*)theMX->GetValue();
}

//___________________________________________________________________________________________
bool    _CategoryVariable::HaveParametersChanged (long catID) {
    for (unsigned long i=0UL; i < parameterList.lLength; i++) {
        _Variable * p = LocateVar(parameterList.list_data[i]);
        if (p->HasChanged())
            if (catID == -1 || ((_SimpleList**)affectedClasses.list_data)[i]->list_data[catID]) {
                return true;
            }
    }

    if (catID == -1 && is_hidden_markov()) {
        return GetHiddenMarkov()->HasChanged() || GetHiddenMarkovFreqs() -> HasChanged();
    }
    
    return false;
}

//___________________________________________________________________________________________
bool    _CategoryVariable::IsConstant (void)
{
    for (unsigned long i=0; i<parameterList.lLength; i++)
        if (LocateVar(parameterList.list_data[i])->IsConstant() == false) {
            return false;
        }

    return true;
}

//___________________________________________________________________________________________
bool    _CategoryVariable::IsGlobal (void)
{
    return true;
}

//___________________________________________________________________________________________
void      _CategoryVariable::ScanForVariables (_AVLList& l, bool globals, _AVLListX * tagger, long weight) const {
    density.ScanFForVariables(l,true, false, true, false, tagger, weight);
    weights->ScanForVariables(l,true,tagger, weight);
    values->ScanForVariables(l,true,tagger, weight);

    if (hiddenMarkovModel != HY_NO_MODEL) {
        GetHiddenMarkov()->ScanForVariables (l,true, tagger, weight);
        GetHiddenMarkovFreqs()->ScanForVariables (l,true, tagger, weight);
    }
    if (intervalSplitter != HY_NO_MODEL) {
        LocateVar(intervalSplitter)->ScanForVariables (l, globals, tagger, weight);
    }

    if (globals) {
        l.Delete ((BaseRef)(hy_x_variable->get_index()));
    }

}

//___________________________________________________________________________________________
void      _CategoryVariable::ScanForGVariables (_AVLList& l)
{
    _SimpleList temp;
    {
        _AVLList    tempA (&temp);
        density.ScanFForVariables(tempA,true);
        weights->ScanForVariables(tempA,true);
        values->ScanForVariables(tempA,true);
        if (hiddenMarkovModel != HY_NO_MODEL) {
            GetHiddenMarkov()->ScanForVariables (tempA,true);
            GetHiddenMarkovFreqs()->ScanForVariables (tempA,true);
        }

        tempA.ReorderList();
    }

    long xi = hy_x_variable->get_index();

    for (long i=0; i<temp.lLength; i++) {
        if (temp.list_data[i]!=xi) {
            _Variable* theV = LocateVar(temp.list_data[i]);

            if (theV->IsGlobal()&& theV->IsIndependent()) {
                l.Insert ((BaseRef)temp.list_data[i]);
            }
        }
    }
}

//___________________________________________________________________________________________
hyFloat      _CategoryVariable::Mean (void)
{
    hyFloat mean = 0.;
    _Matrix * wts = GetWeights(),
              * val = GetValues();

    for (long ii = 0; ii < intervals; ii++) {
        mean += wts->theData[ii] * val->theData[ii];
    }

    return mean;
}

//___________________________________________________________________________________________
bool        _CategoryVariable::UpdateIntervalsAndValues (bool force) {
    if (density.IsEmpty()) {
        return false;
    }

    long i=0;

    if (intervalSplitter >= 0) {
        force = ((_CategoryVariable*)LocateVar(intervalSplitter))->UpdateIntervalsAndValues(force);
    }

    if (!force) {
        force = density.HasChanged();
        if (!force) {
            force = (*values)[0]==0.0;
        }
    }

    if (force) { // stuff to do
        _Matrix* ew;

        if (intervalSplitter >= 0) {
            ew = GetWeights();
        } else {
            ew = (_Matrix*)weights->ComputeNumeric();
        }

        if (!weights->IsIndependent() && !checkWeightMatrix(*ew)) {
            HandleApplicationError (_String("Matrix of category weights invalid at runtime: ") & _String((_String*)ew->toStr()));
        }

        hyFloat currentBase = 0.0,
                   currentLeft = x_min;

        for (i = 0; i<intervals-1; i++) {
            hyFloat    currentProb  = (*ew)[i];
            currentBase+=currentProb;

            if (!cumulative.IsEmpty()) {
                (*intervalEnds)[i] = cumulative.Newton (density, currentBase, currentLeft,x_max, hy_x_variable);    // get the next interval point
            } else {
                (*intervalEnds)[i] = density.Newton ( hy_x_variable, currentBase,x_min, currentLeft);    // get the next interval point
            }

            if (currentProb) {
                if (representation == MEAN) { // compute the MEAN
                    // need to perform integration here of p(x) dx
                    if (meanC.IsEmpty()) {
                        values->theData[i] = density.MeanIntegral (hy_x_variable,currentLeft,(*intervalEnds)[i])/(*ew)[i];
                    } else {
                        _Constant    currentRight ((*intervalEnds)[i]);
                        hy_x_variable->SetValue(&currentRight, true, true, NULL);
                        values->theData[i] = meanC.Compute()->Value();
                        currentRight.SetValue(currentLeft);
                        hy_x_variable->SetValue(&currentRight, true, true, NULL);
                        values->theData[i] = x_min+((*values)[i]- meanC.Compute()->Value())/(*ew)[i];
                        if (values->theData[i]>x_max) {
                            values->theData[i] = x_max;
                        }
                        if (i) {
                            if (values->theData[i]<intervalEnds->theData[i-1]) {
                                values->theData[i]=intervalEnds->theData[i-1];
                            }
                        } else if (values->theData[i]<x_min) {
                            values->theData[i]=x_min;
                        }
                    }

                } else { // compute the MEDIAN
                    if (!cumulative.IsEmpty()) {
                        (*values)[i] = cumulative.Newton (density,currentBase-(*ew)[i]/2, currentLeft, x_max ,hy_x_variable);
                    } else {
                        (*values)[i] = density.Newton (hy_x_variable, currentBase-(*ew)[i]/2,x_min, currentLeft);
                    }
                }
            } else {
                (*values)[i] = (*intervalEnds)[i];
            }

            currentLeft = (*intervalEnds)[i];
        }
        // finally do something special for the last interval, since it may be over (a,infinity)

        hyFloat    lastProb  = (*ew)[i];
        if (lastProb) {
            if (representation == MEAN) { // compute the MEAN
                // need to perform integration here of p(x) dx
                if (meanC.IsEmpty()) {
                    (*values)[i] = density.MeanIntegral (hy_x_variable,currentLeft,(*intervalEnds)[i],true)/(*ew)[i];
                } else {
                    _Constant    currentRight (x_max);
                    hy_x_variable->SetValue(&currentRight, true, false, NULL);
                    values->theData[i] = meanC.Compute()->Value();
                    currentRight.SetValue(currentLeft);
                    hy_x_variable->SetValue(&currentRight, true, false, NULL);
                    values->theData[i] = x_min+((*values)[i]- meanC.Compute()->Value())/(*ew)[i];
                    if (values->theData[i]>x_max) {
                        values->theData[i] = x_max;
                    }
                    if (i) {
                        if (values->theData[i]<intervalEnds->theData[i-1]) {
                            values->theData[i]=intervalEnds->theData[i-1];
                        }
                    } else {
                        if (values->theData[i]<x_min) {
                            values->theData[i]=x_min;
                        }
                    }
                }
            } else { // compute the MEDIAN
                if (!cumulative.IsEmpty()) {
                    (*values)[i] = cumulative.Newton (density,currentBase+(*ew)[i]/2, currentLeft, x_max, hy_x_variable);
                } else {
                    (*values)[i] = density.Newton (hy_x_variable, currentBase+(*ew)[i]/2,x_min, currentLeft);
                }
            }
        } else {
            (*values)[i] = currentLeft;
        }

        if (representation == SCALED_MEDIAN) {

            hyFloat distMean,discMean = 0;

            if (meanC.IsEmpty()) {
                distMean = density.MeanIntegral (hy_x_variable,x_min,x_max, true);
            } else {
                hy_x_variable->SetValue(new _Constant (x_max), false, false, NULL);
                distMean = meanC.Compute()->Value();
            }

            for (i=0; i<intervals; i++) {
                discMean +=(*values)[i]*(*ew)[i];
            }
            discMean = distMean/discMean;
            for (i=0; i<intervals; i++) {
                (*values)[i]*=discMean;
            }
        }
    }
    hy_x_variable->MarkDone();

    return force;
}


//___________________________________________________________________________________________
void _CategoryVariable::SerializeCategory (_StringBuffer & rec)
{
    _String     weightNames = *GetName()&'.'&"weights",
                catNames    = *GetName()&'.'&"points";

    if (intervalSplitter>=0) {
        ((_CategoryVariable*)LocateVar(intervalSplitter))->SerializeCategory(rec);
    }

    bool        hasDensity = (!density.IsEmpty());

    rec << '\n';
    if (intervalSplitter==-1) {
        weights->Serialize (rec, weightNames);
    }

    rec << '\n';
    if (!hasDensity) {
        values->Serialize (rec, catNames);
    }

    rec << '\n';
    if (hiddenMarkovModel != HY_NO_MODEL) {
        SerializeModel (rec,hiddenMarkovModel);
    }

    rec << "\ncategory ";
    rec << *GetName();
    rec << "=(";
    rec << _String ((long)intervals);
    rec << ',';

    if (intervalSplitter==-1) {
        rec << weightNames;
    } else {
        rec << *LocateVar(intervalSplitter)->GetName();
    }

    rec << ',';
    switch (representation) {
    case MEDIAN:
        rec << medianRep;
        break;
    case SCALED_MEDIAN:
        rec << scaledMedianRep;
        break;
    default:
        rec << "MEAN";
        break;
    }
    rec << ',';
    if (hasDensity) {
        rec.AppendNewInstance((_String*)density.toStr(kFormulaStringConversionNormal));
        rec << ',';
        rec.AppendNewInstance((_String*)cumulative.toStr(kFormulaStringConversionNormal));
    } else {
        if (IsUncorrelated()) {
            rec << ',';
        } else {
            rec << LocateVar (covariant)->GetName();
            rec << ',';
        }
        rec << catNames;
    }
    rec << ',';
    rec << _String(x_min-SLIGHT_SHIFT);
    rec << ',';
    rec << _String(x_max+SLIGHT_SHIFT);
    rec << ',';
    rec.AppendNewInstance ((_String*)meanC.toStr(kFormulaStringConversionNormal));

    if ((hiddenMarkovModel != HY_NO_MODEL)||(flags&CONSTANT_ON_PARTITION)) {
        rec << ',';
        if (hiddenMarkovModel != HY_NO_MODEL) {
            rec << *(_String*)modelNames (hiddenMarkovModel);
        }

        if (flags&CONSTANT_ON_PARTITION) {
            rec << constantOnPartition;
        }
    }
    rec << ");\n";
}
