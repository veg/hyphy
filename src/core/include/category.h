/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2006
Primary Development:
  Sergei L Kosakovsky Pond (sergeilkp@mac.com)
Significant contributions from:
  Spencer V Muse (muse@stat.ncsu.edu)
  Simon DW Frost (sdfrost@ucsd.edu)
  Art FY Poon    (apoon@biomail.ucsd.edu)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#ifndef     __CATEGORY__
#define     __CATEGORY__

#include  "parser.h"
#include  "batchlan.h"
#include  "matrix.h"

#define    MEAN                     0
#define    MEDIAN                   1
#define    SCALED_MEDIAN            2
#define    CONSTANT_ON_PARTITION    1

#ifndef    __HYALTIVEC__
#define    INFINITE_BOUND 1e50
#else
#define    INFINITE_BOUND 1e10
#endif

//__________________________________________________________________________________

class     _CategoryVariable: public _Variable
{

public:
    // c&d

    _CategoryVariable () {
        values = nil;
        intervalEnds = nil;
        weights = nil;
        conditionalWeights = nil;
    };
    _CategoryVariable (_CategoryVariable& cv) {
        Duplicate (&cv);
    }
    _CategoryVariable (_String& name, _List* parms, _VariableContainer*);

    // std functions
    virtual
    ~_CategoryVariable () {
        DeleteObject (values);
        DeleteObject (intervalEnds);
        DeleteObject (weights);
    };
    virtual
    BaseRef     makeDynamic             (void);
    virtual
    void        Duplicate               (BaseRef);
    virtual
    BaseRef     toStr                   (void);

    virtual
    bool        IsGlobal                (void);

    virtual
    bool        IsConstant              (void);

    virtual
    bool        IsCategory              (void) {
        return true;
    }

    virtual
    void       ScanForVariables         (_AVLList&, bool = false);

    virtual
    void       ScanForGVariables        (_AVLList&);

    bool       HaveParametersChanged    (long = -1);

    /*virtual
        bool       IsIndependent (void) { return false;} */

    // access functions

    long        GetNumberOfIntervals () {
        return intervals;
    }

    char        GetRepresentationType () {
        return representation;
    }

    _Parameter  SetIntervalValue (long, bool recacl = true);
    // set interval value is returned

    _Parameter  Mean (void);

    _Parameter  GetIntervalValue (long);

    _Parameter  GetIntervalWeight(long);

    _Parameter* GetIntervalWeights(void) {
        return weights->fastIndex();
    }

    _Matrix*    GetValues (void);

    _Matrix*    GetWeights(bool = false);

    _Matrix*    GetIntervalEnds (void);

    _Matrix*    ComputeHiddenMarkov (void);
    _Matrix*    ComputeHiddenMarkovFreqs (void);
    _Matrix*    GetHiddenMarkov (void);
    _Matrix*    GetHiddenMarkovFreqs (void);

    _Formula&   GetDensity(void) {
        return density;
    }

    _Formula&   GetCumulative(void) {
        return cumulative;
    }


    bool        Refresh(bool force=false) {
        return UpdateIntervalsAndValues(force);
    }

    _Parameter  GetMinX (void)  {
        return x_min;
    }
    _Parameter  GetMaxX (void)  {
        return x_max;
    }
    bool        IsHiddenMarkov
    (void)  {
        return (hiddenMarkovModel!=-1);
    }

    bool        IsConstantOnPartition
    (void)  {
        return (flags&CONSTANT_ON_PARTITION);
    }

    void        ChangeNumberOfIntervals (long);
    // assumes a 'standard' category variable - i.e.
    // EQUAL freqs, and density/cumulative

    void        SerializeCategory       (_String&);

    long        GetCurrentState         (void);
    bool        IsUncorrelated          (void);
    bool        IsLayered               (void);

private:

    bool        UpdateIntervalsAndValues (bool force = false);
    void        Construct   (_List&, _VariableContainer*);
    void        Clear       (void);
    bool        checkWeightMatrix (_Matrix&, long = -1);

    // data members
private:

    long        intervals,  // number of intervals
                hiddenMarkovModel,
                covariant,
                intervalSplitter;


    char        flags;
    _Formula    density,
                cumulative,
                meanC;
    // weights of intervals
    // probability density function, in terms of parameters and _x_
    // cumulative prob function, in terms of parameters and _x_
    char        representation;
    // how to represent intervals, by means or medians

    _Matrix     *values,
                *intervalEnds,
                *weights,
                *conditionalWeights;

    _Parameter  x_min,
                x_max;      // distribution range

    _SimpleList parameterList;
    _List       affectedClasses;



};

//__________________________________________________________________________________

extern  _Parameter  maxCategoryIntervals;
extern  _Variable*  _x_, *_n_;

#endif

