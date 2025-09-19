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

/**
 * @brief A category variable
 *
 */
class     _CategoryVariable: public _Variable
{

public:
    // c&d

    /**
     * @brief Construct a new _CategoryVariable object
     *
     */
    _CategoryVariable () {
        values = nil;
        intervalEnds = nil;
        weights = nil;
        conditionalWeights = nil;
    };
    
    /**
     * @brief Construct a new _CategoryVariable object
     *
     * @param cv
     */
    _CategoryVariable (_CategoryVariable const & cv) {
        Duplicate (&cv);
    }
    
    /**
     * @brief Construct a new _CategoryVariable object
     *
     * @param name
     * @param parms
     * @param vc
     */
    _CategoryVariable (_String& name, _List* parms, _VariableContainer*);
    
    /**
     * @brief Assignment operator
     *
     * @param cv
     * @return _CategoryVariable const&
     */
    _CategoryVariable const & operator = (_CategoryVariable const& cv) {
        if (this != &cv) {
            Duplicate (&cv);
        }
        return *this;
    }

    // std functions
    /**
     * @brief Destroy the _CategoryVariable object
     *
     */
    virtual
    ~_CategoryVariable () {
        BatchDelete(values, intervalEnds, weights, conditionalWeights);
    };
    /**
     * @brief Make a dynamic copy of the object
     *
     * @return BaseRef
     */
    virtual
    BaseRef     makeDynamic             (void) const;
    /**
     * @brief Duplicate the object from a reference
     *
     * @param brc
     */
    virtual
    void        Duplicate               (BaseRefConst);
    /**
     * @brief Convert the object to a string
     *
     * @param ul
     * @return BaseRef
     */
    virtual
    BaseRef     toStr                   (unsigned long = 0UL);

    /**
     * @brief Check if the variable is global
     *
     * @return true
     * @return false
     */
    virtual
    bool        IsGlobal                (void);

    /**
     * @brief Check if the variable is constant
     *
     * @return true
     * @return false
     */
    virtual
    bool        IsConstant              (void);

    /**
     * @brief Check if the variable is a category
     *
     * @return true
     * @return false
     */
    virtual
    bool        IsCategory              (void) {
        return true;
    }

    /**
     * @brief Scan for variables
     *
     * @param avl
     * @param b
     * @param tagger
     * @param weight
     */
    virtual
    void       ScanForVariables         (_AVLList&, bool = false, _AVLListX* tagger = nil, long weight = 0) const;

    /**
     * @brief Scan for global variables
     *
     * @param avl
     */
    virtual
    void       ScanForGVariables        (_AVLList&);

    /**
     * @brief Check if the parameters have changed
     *
     * @param l
     * @return true
     * @return false
     */
    bool       HaveParametersChanged    (long = -1);

    /*virtual
        bool       IsIndependent (void) { return false;} */

    // access functions

    /**
     * @brief Get the Number Of Intervals
     *
     * @return long
     */
    long        GetNumberOfIntervals () const {
        return intervals;
    }

    /**
     * @brief Get the Representation Type
     *
     * @return char
     */
    char        GetRepresentationType () const {
        return representation;
    }

    /**
     * @brief Set the Interval Value
     *
     * @param l
     * @param recacl
     * @return hyFloat
     */
    hyFloat  SetIntervalValue (long, bool recacl = true);
    // set interval value is returned

    /**
     * @brief Get the Mean
     *
     * @return hyFloat
     */
    hyFloat  Mean (void);

    /**
     * @brief Get the Interval Value
     *
     * @param l
     * @return hyFloat
     */
    hyFloat  GetIntervalValue (long);

    /**
     * @brief Get the Interval Weight
     *
     * @param l
     * @return hyFloat
     */
    hyFloat  GetIntervalWeight(long);

    /**
     * @brief Get the Interval Weights
     *
     * @return hyFloat*
     */
    hyFloat* GetIntervalWeights(void) {
        return weights->fastIndex();
    }

    /**
     * @brief Get the Values
     *
     * @return _Matrix*
     */
    _Matrix*    GetValues (void);

    /**
     * @brief Get the Weights
     *
     * @param b
     * @return _Matrix*
     */
    _Matrix*    GetWeights(bool = false);

    /**
     * @brief Get the Interval Ends
     *
     * @return _Matrix*
     */
    _Matrix*    GetIntervalEnds (void);

    /**
     * @brief Compute the Hidden Markov model
     *
     * @return _Matrix*
     */
    _Matrix*    ComputeHiddenMarkov (void);
    /**
     * @brief Compute the Hidden Markov frequencies
     *
     * @return _Matrix*
     */
    _Matrix*    ComputeHiddenMarkovFreqs (void);
    /**
     * @brief Get the Hidden Markov model
     *
     * @return _Matrix*
     */
    _Matrix*    GetHiddenMarkov (void) const;
    /**
     * @brief Get the Hidden Markov frequencies
     *
     * @return _Matrix*
     */
    _Matrix*    GetHiddenMarkovFreqs (void) const;

    /**
     * @brief Get the Density
     *
     * @return _Formula&
     */
    _Formula&   GetDensity(void) {
        return density;
    }

    /**
     * @brief Get the Cumulative
     *
     * @return _Formula&
     */
    _Formula&   GetCumulative(void) {
        return cumulative;
    }


    /**
     * @brief Refresh the variable
     *
     * @param force
     * @return true
     * @return false
     */
    bool        Refresh(bool force=false) {
        return UpdateIntervalsAndValues(force);
    }

    /**
     * @brief Get the Min X
     *
     * @return hyFloat
     */
    hyFloat  GetMinX (void)  const {
        return x_min;
    }
    /**
     * @brief Get the Max X
     *
     * @return hyFloat
     */
    hyFloat  GetMaxX (void)  const {
        return x_max;
    }
    /**
     * @brief Check if the model is a hidden markov model
     *
     * @return true
     * @return false
     */
    bool        is_hidden_markov (void)  const {
        return (hiddenMarkovModel!=-1);
    }

    /**
     * @brief Check if the model is constant on partition
     *
     * @return true
     * @return false
     */
    bool        is_constant_on_partition (void) const {
        return (flags&CONSTANT_ON_PARTITION);
    }

    /**
     * @brief Change the number of intervals
     *
     * @param l
     */
    void        ChangeNumberOfIntervals (long);
    // assumes a 'standard' category variable - i.e.
    // EQUAL freqs, and density/cumulative

    /**
     * @brief Serialize the category
     *
     * @param sb
     */
    void        SerializeCategory       (_StringBuffer &);

    /**
     * @brief Get the Current State
     *
     * @return long
     */
    long        GetCurrentState         (void);
    /**
     * @brief Check if the model is uncorrelated
     *
     * @return true
     * @return false
     */
    bool        IsUncorrelated          (void);
    /**
     * @brief Check if the model is layered
     *
     * @return true
     * @return false
     */
    bool        IsLayered               (void);

private:

    /**
     * @brief Update the intervals and values
     *
     * @param force
     * @return true
     * @return false
     */
    bool        UpdateIntervalsAndValues (bool force = false);
    /**
     * @brief Construct the object
     *
     * @param l
     * @param vc
     */
    void        Construct   (_List&, _VariableContainer*);
    /**
     * @brief Clear the object
     *
     */
    void        Clear       (void);
    /**
     * @brief Check the weight matrix
     *
     * @param m
     * @param l
     * @return true
     * @return false
     */
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

/**
 * @brief Integrate over all possible assignments of a list of category variables.
 *
 * @tparam ACTION A function object that takes the current category and the weight as arguments.
 * @param indices A list of indices of the category variables.
 * @param refresh Whether to refresh the category variables.
 * @param do_this The function to execute for each assignment.
 */
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
