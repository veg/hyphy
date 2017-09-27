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

#if not defined __AFYP_REWRITE_BGM__

#include "simplelist.h"
#include "list.h"
#include "classes.h"
#include "likefunc.h"
#include "parser.h"
#include <math.h>
#include "matrix.h"
#include "baseobj.h"
#include "batchlan.h"
// #include "HYUtils.h"


/*SLKP 20070926; include progress report updates */
#if !defined __UNIX__ && !defined __HEADLESS__
#include "HYConsoleWindow.h"
#endif
/*SLKP*/


#if defined __AFYP_DEVELOPMENT__ && defined __HYPHYMPI__
#include "mpi.h"
#endif


#define     LOG_SCALING_FACTOR          64.
#define     LARGE_NEGATIVE_NUMBER       -99999.0
#define     MAX_LSE_SCALING_ATTEMPTS    4096

#define     MAX_FAIL_RANDOMIZE          1000
#define     RANDOMIZE_PROB_SWAP         0.1

#define     DIRICHLET_FLATTENING_CONST  0.5




/* ____________________________________________________________________________________________________________________________________ */
class Bgm : public _LikelihoodFunction
{
public:
    Bgm ()  {
        /* default constructor does diddly squat */
    }

    Bgm (_AssociativeList *, _AssociativeList *);       // constructor

    virtual ~Bgm (void);        // destructor

    virtual _Parameter      Compute (void);         // function polymorphism, no argument returns likelihood of network
    // specified by _Matrix object dag





    virtual _Matrix *       Optimize ();    // estimate maximum posterior network via greedy hill-climbing heuristic

    //virtual _Matrix *     Optimize (_SimpleList *);

    virtual _PMathObj       CovarianceMatrix (_SimpleList *);   // wrapper function for MCMC, to explore posterior distribution



    virtual void            SetDataMatrix (_Matrix *);

    void            SetWeightMatrix (_Matrix *),
                    SetGraphMatrix (_Matrix *),
                    SetBanMatrix (_Matrix *),
                    SetEnforceMatrix (_Matrix *),
                    SetBestOrder (_SimpleList *);

    _Matrix *       ExportNodeScores (void);
    _Matrix *       ExportGraph (void);

    void            SerializeBgm (_String &);

    void            ImportNodeScores (_Matrix *);


    _Parameter      ComputeDiscreteScore (long node_id),    // compute K2 or BDeu scoring metric for a discrete node in network
                    // with discrete-valued parents, given data and prior
                    // use BDeu if prior_sample_size > 0.
                    ComputeDiscreteScore (long, _Matrix *),
                    ComputeDiscreteScore (long, _SimpleList &),

                    ComputeContinuousScore (long node_id),  // compute scoring metric for continuous node with both discrete
                    // and continuous parents.
                    ComputeContinuousScore (long, _Matrix *),
                    ComputeContinuousScore (long, _SimpleList &, _SimpleList &);

    virtual _Parameter      ImputeDiscreteScore (long, _SimpleList &);

    _Parameter      GibbsApproximateDiscreteScore (long, _SimpleList &),
                    K2Score (long, _Matrix &, _Matrix &),
                    BDeScore (long, _Matrix &, _Matrix &);

    long            GetNumNodes (void)      {
        return num_nodes;
    }
    long            GetNumCases (void)      {
        return (obsData ? obsData->GetHDim() : 0);
    }
    _SimpleList     GetLastOrder (void)     {
        return last_node_order;
    }

protected:
    _Parameter      Compute (_SimpleList*, _List*),         // return likelihood of node ordering, model averaging over all
                    // possible networks that are consistent with the order.
                    Compute (_Matrix *),
                    ComputeDynamic (_SimpleList *, _List *),    // compute for DBNs
                    ComputeDynamic2 (_SimpleList *, _List *);


    virtual void            CacheNodeScores (void);
    void            ReleaseNodeScores (void);

    void            MPICacheNodeScores (long);
    void            MPIReceiveScores (_Matrix *, bool, long);
    void            CacheNetworkParameters (void);


    void            InitComputeLists (_List *);
    void            DumpComputeLists (_List *);



    _Parameter      LnGamma (_Parameter),
                    LogSumExpo (_GrowingVector *);

    _Matrix *       RunColdChain (_SimpleList *, long, long);
    void            RunHotChain (_SimpleList *, long, long, _Parameter);

    _Matrix *       GraphMCMC (bool);

    _Parameter      TryEdge (long, long, long, _Parameter);     // DEPRECATE

    void            RandomizeDag (long),
                    RandomizeGraph (_Matrix *, _SimpleList *, long, bool),
                    PrintGraph (_Matrix *),
                    ResetGraph (_Matrix *);

    long            MarginalSum (bool, long);   // DEPRECATE

    bool            IsCyclic (void);    // DEPRECATE



    long            num_nodes,          // number of variables in network

                    max_max_parents;    // the global maximum number of parents per node


    _Constant   *   calc_bit;           // provides access to mathematical functionality of _Constant


    _Matrix         dag,                // i.e. directed acyclic graph as a square adjacency matrix,
                    // where 1 indicates presence of edge (row --> col)


                    banned_edges,       // complements DAG, 1 entry indicates banned edge
                    enforced_edges,

                    prior_sample_size,  // prior distributions for hyperparameters
                    prior_mean,
                    prior_precision;


    _SimpleList     is_discrete,        // vector containing 1 for discrete node, 0 for continuous

                    num_levels,         // vector containing number of levels for discrete nodes

                    max_parents,        // maximum number of parents allowed per node

                    has_missing,        // data matrix for i-th node contains missing values

                    last_node_order,    // store node orderings from order-MCMC
                    best_node_order;

    _List           node_scores,
                    network_parameters;


    _Matrix     *   obsData,            // data matrix read from file

                * obsWeights;         // vector of weights corresponding to each observation in data matrix


    bool            scores_cached;
};



/* ____________________________________________________________________________________________________________________________________ */

class _DynamicBgm : public Bgm
{

public:
    _DynamicBgm (_AssociativeList *, _AssociativeList *);

    virtual ~_DynamicBgm (void)         { } // destructor

    virtual void        CacheNodeScores     (void);
    void        CollapseNodeSpace   (void);

    void        CollapseDynamicGraph (void);
private:

};



#endif


