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
    
    _CategoryVariable (_CategoryVariable const & cv) {
        Duplicate (&cv);
    }
    
    _CategoryVariable (_String& name, _List* parms, _VariableContainer*);
    
    _CategoryVariable const & operator = (_CategoryVariable const& cv) {
        if (this != &cv) {
            Duplicate (&cv);
        }
        return *this;
    }

    // std functions
    virtual
    ~_CategoryVariable () {
        BatchDelete(values, intervalEnds, weights, conditionalWeights);
    };
    virtual
    BaseRef     makeDynamic             (void) const;
    virtual
    void        Duplicate               (BaseRefConst);
    virtual
    BaseRef     toStr                   (unsigned long = 0UL);

    virtual
    bool        IsGlobal                (void);

    virtual
    bool        IsConstant              (void);

    virtual
    bool        IsCategory              (void) {
        return true;
    }

    virtual
    void       ScanForVariables         (_AVLList&, bool = false, _AVLListX* tagger = nil, long weight = 0) const;

    virtual
    void       ScanForGVariables        (_AVLList&);

    bool       HaveParametersChanged    (long = -1);

    /*virtual
        bool       IsIndependent (void) { return false;} */

    // access functions

    long        GetNumberOfIntervals () const {
        return intervals;
    }

    char        GetRepresentationType () const {
        return representation;
    }

    hyFloat  SetIntervalValue (long, bool recacl = true);
    // set interval value is returned

    hyFloat  Mean (void);

    hyFloat  GetIntervalValue (long);

    hyFloat  GetIntervalWeight(long);

    hyFloat* GetIntervalWeights(void) {
        return weights->fastIndex();
    }

    _Matrix*    GetValues (void);

    _Matrix*    GetWeights(bool = false);

    _Matrix*    GetIntervalEnds (void);

    _Matrix*    ComputeHiddenMarkov (void);
    _Matrix*    ComputeHiddenMarkovFreqs (void);
    _Matrix*    GetHiddenMarkov (void) const;
    _Matrix*    GetHiddenMarkovFreqs (void) const;

    _Formula&   GetDensity(void) {
        return density;
    }

    _Formula&   GetCumulative(void) {
        return cumulative;
    }


    bool        Refresh(bool force=false) {
        return UpdateIntervalsAndValues(force);
    }

    hyFloat  GetMinX (void)  const {
        return x_min;
    }
    hyFloat  GetMaxX (void)  const {
        return x_max;
    }
    bool        is_hidden_markov (void)  const {
        return (hiddenMarkovModel!=-1);
    }

    bool        is_constant_on_partition (void) const {
        return (flags&CONSTANT_ON_PARTITION);
    }

    void        ChangeNumberOfIntervals (long);
    // assumes a 'standard' category variable - i.e.
    // EQUAL freqs, and density/cumulative

    void        SerializeCategory       (_StringBuffer &);

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

    hyFloat  x_min,
                x_max;      // distribution range

    _SimpleList parameterList;
    _List       affectedClasses;



};

template <typename ACTION> void IntergrateOverAssignments (_SimpleList const& indices, bool refresh, ACTION&& do_this) {
    
    long                    current_category = 0L,
                            total_cat_count  = 1L,
                            cat_var_count    = indices.countitems();
    
    hyFloat                 weight           = 1.;
    
    
    indices.Each ( [&] (long cat_var_idx, unsigned long) -> void {
        _CategoryVariable*      cat_var = (_CategoryVariable*)LocateVar (cat_var_idx);
        if (refresh) {
            cat_var->Refresh();
        }
        total_cat_count *= cat_var->GetNumberOfIntervals();
    });

    do {
        if (indices.nonempty()) {
            long c = current_category;
            weight = 1.0;
            for (long k=cat_var_count-1; k>=0; k--) {
                _CategoryVariable*      cat_var = (_CategoryVariable*)LocateVar (indices.get (k));
                
                long t          = cat_var->GetNumberOfIntervals(),
                     this_index = c % t;
                
                cat_var->SetIntervalValue(this_index);
                weight *= cat_var->GetIntervalWeight(this_index);
                c/=t;
            }
        }

        do_this (current_category, weight);
        current_category++;
    } while (current_category<total_cat_count);
}

//__________________________________________________________________________________

extern  unsigned long  maxCategoryIntervals;

#endif
