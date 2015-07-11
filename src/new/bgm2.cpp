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

#include "bgm.h"



// this old method is WAY too slow, constantly allocating and freeing _Constant objects

/*
inline _Parameter Bgm::LnGamma(_Constant * calculator, _Parameter x)
{
    // wrapper function for _Constant member function
    calculator->SetValue (x);
    calculator = (_Constant *) calculator->LnGamma();
    _Parameter rv = calculator->Value();
    DeleteObject(calculator);
    return rv;
}
*/

//___________________________________________________________________________________________

/* maxParentString ("Bgm_MAXIMUM_PARENTS"), */
extern _String      mcemMaxSteps,   // modify string constant name to avoid conflicts with bgm.cpp
       mcemBurnin,
       mcemSampleSize;

_String     _HYBgm_NThreads ("BGM_NTHREADS");


#ifdef __MP__
#include <pthread.h>

struct  ThreadCacheTask {
    Bgm     * b;

    long    startNode,
            nextNode,   // one past the last node
            numNodes;

    _List   *       nodeScores;

    _SimpleList     maxParents,     // pass member variables to global scope
                    isDiscrete;
};

pthread_t *         BgmThreads  = nil;
ThreadCacheTask *   BgmTasks    = nil;
void *              CacheNodeScoreThread (void *);
#endif


//___________________________________________________________________________________________
void    Bgm::ImportNodeScores (_Matrix * import_scores)
{
    /* UNDER DEVELOPMENT!!! */
}



//___________________________________________________________________________________________
//  This doesn't work yet -- afyp :-/
_Matrix *   Bgm::ExportNodeScores (void)
{
    // determine the required number of rows
    long        nrows = 0;

    for (long child = 0; child < num_nodes; child++) {
        long    maxp = max_parents.lData[child];
        nrows++;
        if (maxp > 0) {
            nrows += num_nodes - 1;
            if (maxp > 1) {
                _List   * childs_list   = (_List *) node_scores.lData[child];
                for (long np = 2; np < maxp; np++) {
                    _NTupleStorage  * nts   = (_NTupleStorage *) childs_list->lData[np];
                    nrows += nts->GetSize();
                }
            }
        }
    }


    _Matrix *   export_mx = new _Matrix (nrows, 4, false, true);        // child ID, #parents, parent combo #, score


    for (long index = 0, child = 0; child < num_nodes; child++) {
        _List   * childs_list   = (_List *) node_scores.lData[child];

        for (long np = 0; np < max_parents.lData[child]; np++) {
            if (np == 0) {
                _Constant   *   orphan_score = (_Constant *) childs_list->lData[0];

                export_mx->Store (index, 0, child);
                export_mx->Store (index, 1, np);
                export_mx->Store (index, 2, 0);
                export_mx->Store (index, 3, orphan_score->Value());

                index++;
            } else if (np == 1) {
                _Matrix *   single_parent_scores = (_Matrix *) childs_list->lData[1];

                for (long parent = 0; parent < num_nodes; parent++) {
                    if (parent == child) {
                        continue;
                    }
                    export_mx->Store (index, 0, child);
                    export_mx->Store (index, 1, np);
                    export_mx->Store (index, 2, parent);
                    export_mx->Store (index, 3, (*single_parent_scores) (parent, 0));
                    index++;
                }
            } else {
                _NTupleStorage *    family_scores = (_NTupleStorage *) childs_list->lData[np];

                for (long family = 0; family < family_scores->GetSize(); family++) {
                    export_mx->Store (index, 0, child);
                    export_mx->Store (index, 1, np);
                    export_mx->Store (index, 2, family);    // direct index
                    export_mx->Store (index, 3, family_scores->DirectIndex(family));
                    index++;
                }
            }
        }
    }


    return (_Matrix *)(export_mx->makeDynamic());
}



//___________________________________________________________________________________________
//  THIS DOES NOT WORK YET --- need to provide global scope for member variable node_scores,
//  through accessor function maybe?
#ifdef __NOT_DEFINED_MP__
void * CacheNodeScoreThread (void * arg)
{
    ThreadCacheTask *   theTask = (ThreadCacheTask *) arg;

    long                num_nodes   = theTask->numNodes;
    _SimpleList         is_discrete = theTask->isDiscrete;

    for (long node_id = theTask->startNode; node_id < theTask->nextNode; node_id++) {
        long        maxp        = theTask->maxParents.lData[node_id];
        _List   *   this_list   = (_List *) (theTask->nodeScores)->lData[node_id];


        this_list->Clear();

        /*
         for (long par = 0; par < num_nodes; par++) // reset local graph at this node
         {
         dag.Store (par, node_id, 0);
         }
         */

        // handle case of no parents
        _SimpleList parents;
        _Parameter  score = is_discrete.lData[node_id] ? theTask->b->ComputeDiscreteScore (node_id, parents) : theTask->b->ComputeContinuousScore (node_id);
        _Constant   orphan_score (score);

        (*this_list) && (&orphan_score);        // _List Append() specific to objects in BaseObj hierarchy


        char buf [256];
        snprintf (buf, sizeof(buf), "orphan score = %f\n", score);
        BufferToConsole (buf);


        // handle case of one parent
        if (maxp > 0) {
            _Matrix     single_parent_scores (num_nodes, 1, false, true);

            for (long par = 0; par < num_nodes; par++) {
                if (par == node_id) {   // unused
                    continue;
                }

                //dag.Store (par, node_id, 1);
                parents << par;
                single_parent_scores.Store (par, 0, is_discrete.lData[node_id] ? theTask->b->ComputeDiscreteScore (node_id, parents) :
                                            theTask->b->ComputeContinuousScore (node_id));
                parents.Clear();
                //dag.Store (par, node_id, 0);  // reset
            }
            (*this_list) && (&single_parent_scores);
        }


        // handle cases of more than one parent using (n,k)-tuple indexing
        if (maxp > 1) {
            _SimpleList     all_but_one (num_nodes-1, 0, 1),    // 0, 1, 2, ... , n-1
                            aux_list,
                            nk_tuple;

            for (long np = 2; np <= maxp; np++) {
                _NTupleStorage  family_scores (num_nodes-1, np);

                if (all_but_one.NChooseKInit (aux_list, nk_tuple, np, false)) {
                    bool    remaining;
                    long    res;

                    do {
                        remaining = all_but_one.NChooseK (aux_list, nk_tuple);


                        for (long par_idx = 0; par_idx < np; par_idx++) {
                            long par = nk_tuple.lData[par_idx];

                            if (par >= node_id) {
                                par++;    // map from (n-1) parent space to nodespace
                            }
                            parents << par;
                            //  dag.Store (par, node_id, 1);
                        }


                        score = is_discrete.lData[node_id] ? theTask->b->ComputeDiscreteScore (node_id, parents) :
                                theTask->b->ComputeContinuousScore (node_id);

                        res = family_scores.Store (score, nk_tuple);

                        parents.Clear();
                        /*
                         for (long par = 0; par < num_nodes; par++)
                         {
                         dag.Store (par, node_id, 0);   // reset graph
                         }
                         */

                    } while (remaining);

                } else {
                    _String oops ("Failed to initialize _NTupleStorage object in Bgm::CacheNodeScores().\n");
                    WarnError(oops);
                }

                (*this_list) && (&family_scores);   // append duplicate to storage
            }
        }
    }
    pthread_exit (NULL);
}
#endif



//___________________________________________________________________________________________
_Parameter  Bgm::ImputeDiscreteScore (long node_id, _SimpleList & parents)
{
    char            buf [256];

    //  Use Monte Carlo over assignments to m_ijk to calculate expectation, which can be used for EM (MCEM)
    //      but m_ijk are constrained by partial observations, a hassle for bookkeeping  :-/
    //      Alternatively, we can perform MCMC over imputations of missing values; which is faster?

    _SimpleList     multipliers ((long)1);

    _Matrix     n_ijk, n_ij;

    long        num_parent_combos   = 1,
                r_i                   = num_levels.lData[node_id],
                mcem_interval;

    _Parameter  log_score   = 0,
                mcem_max_steps,
                mcem_burnin,
                mcem_sample_size;


    // set MCMC parameters from batch language definitions
    checkParameter (mcemMaxSteps, mcem_max_steps, 0);
    checkParameter (mcemBurnin, mcem_burnin, 0);
    checkParameter (mcemSampleSize, mcem_sample_size, 0);

    if (mcem_max_steps == 0 || mcem_sample_size == 0) {
        _String oops ("Did you forget to specify a value for BGM_MCEM_MAXSTEPS or BGM_MCEM_SAMPLES?\n");
        WarnError (oops);
    }

    mcem_interval = (long) (mcem_max_steps / mcem_sample_size);


    // count number of parent state combinations
    for (long par = 0; par < parents.lLength; par++) {
        num_parent_combos *= num_levels.lData[parents.lData[par]];
        multipliers << num_parent_combos;
    }

    CreateMatrix (&n_ijk, num_parent_combos, r_i, false, true, false);
    CreateMatrix (&n_ij, num_parent_combos, 1, false, true, false);



    _GrowingVector      impute,     // to store a copy of incomplete cases
                        is_missing; // boolean values mapping to missing values in [impute]

    _GrowingVector  *   post_chain = new _GrowingVector(),
    *   prior_denom = new _GrowingVector();


    // prepare containers for keeping track of state-frequencies per node
    // for later use in imputation

    long                max_num_levels  = num_levels.lData [node_id],
                        family_size        = parents.lLength + 1;

    _SimpleList         family_nlevels (max_num_levels);

    _Matrix             m_ijk,      // let m_ijk denote missing values such that n_ijk + m_ijk sums to N_i for all j, k.
                        m_ij,
                        observed_freqs;


    CreateMatrix (&m_ijk, num_parent_combos, r_i, false, true, false);  // note these are populated with 0.0's
    CreateMatrix (&m_ij, num_parent_combos, 1, false, true, false);


    // evaluate the number of levels for each node in family
    for (long par = 0; par < parents.lLength; par++) {
        long    par_nlevels = num_levels.lData [ parents.lData[par] ];

        family_nlevels << par_nlevels;
        if (par_nlevels > max_num_levels) {
            max_num_levels = par_nlevels;
        }
    }

    CreateMatrix (&observed_freqs, family_size, max_num_levels, false, true, false);



#ifdef __DEBUG_IDS__
    snprintf (buf, sizeof(buf), "family_nlevels: ");
    BufferToConsole (buf);
    for (long bug = 0; bug < family_nlevels.lLength; bug++) {
        snprintf (buf, sizeof(buf), "%d ", family_nlevels.lData[bug]);
        BufferToConsole (buf);
    }
    NLToConsole ();
#endif



    // tally complete cases
    for (long index, obs = 0; obs < obsData->GetHDim(); obs++) {
        long    child_state = (*obsData) (obs, node_id);

        index = 0;

        if (child_state > -1) {
            observed_freqs.Store (0, (long)child_state, observed_freqs (0, (long)child_state) + 1);

            for (long par = 0; par < parents.lLength; par++) {
                long    this_parent         = parents.lData[par],
                        this_parent_state = (long) ((*obsData) (obs, this_parent));

                if (this_parent_state < 0) {
                    index = -1; // flag observation to skip
                    break;
                } else {
                    index += this_parent_state * multipliers.lData[par];
                    observed_freqs.Store (par+1, this_parent_state, observed_freqs (par+1, this_parent_state) + 1);
                }
            }
        } else {
            index = -1;
        }

        if (index > -1) {   // update n_ijk's with complete case
            n_ijk.Store ((long) index, (long)child_state, n_ijk(index, (long)child_state) + 1);
            n_ij.Store ((long) index, 0, n_ij(index, 0) + 1);
        } else {
            // store family in impute, always in the order (child, parent0, parent1, ...)
            impute.Store ( (_Parameter) child_state);
            is_missing.Store ( (_Parameter) ((child_state < 0) ? 1 : 0) );

            for (long par = 0; par < parents.lLength; par++) {
                long    parent_state    = (long) ((*obsData) (obs, parents.lData[par]));

                impute.Store ( (_Parameter) parent_state);
                is_missing.Store ( (_Parameter) ((parent_state < 0) ? 1 : 0) );
            }
        }
    }


#ifdef __DEBUG_IDS__
    snprintf (buf, sizeof(buf), "complete cases N_ijk: \n");
    BufferToConsole (buf);
    for (long j = 0; j < num_parent_combos; j++) {
        for (long k = 0; k < family_nlevels.lData[0]; k++) {
            snprintf (buf, sizeof(buf), "%d ", (long) n_ijk (j,k));
            BufferToConsole (buf);
        }
        NLToConsole ();
    }
    NLToConsole ();

    snprintf (buf, sizeof(buf), "impute: \n");
    BufferToConsole (buf);
    for (long bug = 0; bug < impute.GetUsed(); bug++) {
        snprintf (buf, sizeof(buf), "%d ", (long) impute (0, bug));
        BufferToConsole (buf);

        if (bug % family_size == family_size - 1) {
            NLToConsole ();
        }
    }
    NLToConsole ();
#endif



    // convert observed states into frequencies
    for (long obs_total, i = 0; i < family_size; i++) {
        obs_total = 0;
        for (long obs_state = 0; obs_state < family_nlevels.lData[i]; obs_state++) {
            obs_total += observed_freqs (i, obs_state);
        }

        for (long obs_state = 0; obs_state < family_nlevels.lData[i]; obs_state++) {
            observed_freqs.Store (i, obs_state, (_Parameter) (observed_freqs (i, obs_state) / obs_total));
        }
    }


#ifdef __DEBUG_IDS__
    snprintf (buf, sizeof(buf), "observed_freqs: \n");
    BufferToConsole (buf);

    for (long fam = 0; fam < family_size; fam++) {
        for (long lev = 0; lev < family_nlevels.lData[fam]; lev++) {
            snprintf (buf, sizeof(buf), "%f ", observed_freqs (fam, lev));
            BufferToConsole (buf);
        }
        NLToConsole ();
    }
    NLToConsole ();
#endif



    // Initial imputation with random assignments to missing values, based on observed frequencies
    //  We can always burn-in the chain to converge to more realistic assignments
    //  before taking expectation.
    //  Also use this loop to calculate m_ijk's, i.e. cell counts for cases with missing values
    for (long i = 0; i < impute.GetUsed(); i++) {
        if (is_missing(0,i)) {
            double  urn         = genrand_real2();
            long    member      = i % family_size;

            for (long lev = 0; lev < family_nlevels.lData[member]; lev++) {
                if (urn < observed_freqs (member, lev)) {
                    impute._Matrix::Store (0, i, (_Parameter)lev);
                    break;
                } else {
                    urn -= observed_freqs (member, lev);    // rescale random float
                }
            }
        }
    }



    // compute m_ijk's and m_ij's
    long        num_incomplete_cases    = impute.GetUsed() / family_size;

    for (long ic = 0; ic < num_incomplete_cases; ic++) {
        long    index           = 0,
                child_state     = impute (0, ic * family_size);

        for (long par = 0; par < parents.lLength; par++) {
            long    parent_state    = impute (0, ic*family_size + par+1);
            index += parent_state * multipliers.lData[par];
        }

        m_ijk.Store ((long) index, child_state, m_ijk(index, child_state) + 1);
        m_ij.Store ((long) index, 0, m_ij(index, 0) + 1);
    }



#ifdef __DEBUG_IDS__
    snprintf (buf, sizeof(buf), "impute: \n");
    BufferToConsole (buf);
    for (long bug = 0; bug < impute.GetUsed(); bug++) {
        snprintf (buf, sizeof(buf), "%d ", (long) impute (0, bug));
        BufferToConsole (buf);

        if (bug % family_size == family_size - 1) {
            NLToConsole ();
        }
    }
    NLToConsole ();

    snprintf (buf, sizeof(buf), "m_ijk: \n");
    BufferToConsole (buf);
    for (long j = 0; j < num_parent_combos; j++) {
        for (long k = 0; k < family_nlevels.lData[0]; k++) {
            snprintf (buf, sizeof(buf), "%d ", (long) m_ijk (j,k));
            BufferToConsole (buf);
        }

        snprintf (buf, sizeof(buf), "| %d", (long) m_ij (j,0));
        BufferToConsole (buf);

        NLToConsole ();
    }
    NLToConsole ();
#endif


    // Calculate probability of imputation {m_ijk}
    //  For the time being, use multinomial with Dirichlet prior based on non-missing cases (n_ijk)
    //  for each node, i.e., network parameters for empty graph (independence prior).
    //  Later, we will want to allow the user to submit their own prior distribution (e.g., after an iteration
    //  of structural EM) - afyp (March 28, 2008)

    // add Jeffrey's invariance constant (0.5) to every cell n_ijk to avoid zero counts?

    _Parameter  log_prior_impute        = 0.0;

    for (long j = 0; j < num_parent_combos; j++) {

        log_prior_impute += LnGamma(n_ij(j,0) + DIRICHLET_FLATTENING_CONST * r_i)
                            - LnGamma(n_ij(j,0) + DIRICHLET_FLATTENING_CONST * r_i + m_ij(j,0));
        /*
         log_prior_impute += LnGamma(m_ij(j,0)+1) + LnGamma(n_ij(j,0)+DIRICHLET_FLATTENING_CONST*r_i)
         - LnGamma(m_ij(j,0)+n_ij(j,0)+DIRICHLET_FLATTENING_CONST*r_i);
         */
        for (long k = 0; k < r_i; k++) {

            log_prior_impute += LnGamma(n_ijk(j,k) + DIRICHLET_FLATTENING_CONST + m_ijk(j,k))
                                - LnGamma(n_ijk(j,k) + DIRICHLET_FLATTENING_CONST);
            /*
             log_prior_impute += LnGamma(m_ijk(j,k) + n_ijk(j,k)+DIRICHLET_FLATTENING_CONST)
             - LnGamma(m_ijk(j,k)+1) - LnGamma(n_ijk(j,k)+DIRICHLET_FLATTENING_CONST);
             */
        }
    }



    // Calculate posterior probability for initial imputation
    _Parameter  last_score          = log_prior_impute,
                last_prior_impute    = log_prior_impute;


    for (long j = 0; j < num_parent_combos; j++) {
        last_score += LnGamma(num_levels.lData[node_id]);   // (r-1)!
        last_score -= LnGamma(n_ij(j, 0) + m_ij(j,0) + num_levels.lData[node_id]);  // (N+r-1)!

        for (long k = 0; k < num_levels.lData[node_id]; k++) {
            last_score += LnGamma(n_ijk (j,k) + m_ijk (j,k) + 1);   // (N_ijk)!
        }
    }



#ifdef __DEBUG_IDS__
    snprintf (buf, sizeof(buf), "log_prior_impute = %f\nlast_score = %f\n", log_prior_impute, last_score);
    BufferToConsole (buf);
#endif

    /********
     /  MCMC  /
     ********/

    // Run MCMC over imputations, use as proposal function a single case reassignment (guaranteed to modify m_ijk's)
    //      Take into account case being reassigned so that we don't have to re-sum m_ijk every time
    //      use burn-in and thinning of the sample when calculating the expectation for family posterior probability
    _Parameter  lk_ratio,
                expect_log_score = 0.;

    _Matrix     last_impute ((_Matrix &) impute),       // duplicate matrix constructors
                last_m_ijk (m_ijk),
                last_m_ij (m_ij);


    for (long step = 0; step < mcem_max_steps + mcem_burnin; step++) {
        /*  PROPOSAL  */

        // choose a random case from the incomplete set and pick a random variable that is missing
        long    imp,
                imp_case,
                imp_var,
                last_state,
                child_state,
                pick;



        while (1) { // this could be optimized by pre-computing a list of missing value indices
            imp = genrand_int32 () % impute.GetUsed();

            if (is_missing(0, imp)) {
                break;
            }
        }

        imp_case = (long) (imp / family_size);      // which incomplete case?
        imp_var  = imp % family_size;       // which variable?
        last_state = impute (0, imp);



#ifdef __DEBUG_IDS__
        snprintf (buf, sizeof(buf), "imp = %d\nimp_case = %d\nimp_var = %d\nnum_incomplete_cases=%d\n", imp, imp_case, imp_var, num_incomplete_cases);
        BufferToConsole (buf);
#endif


        long    old_index = 0,
                new_index;

        for (long par = 0; par < parents.lLength; par++) {
            long    parent_state    = impute (0, imp_case*family_size + par+1);
            old_index += parent_state * multipliers.lData[par];
        }



        // randomly assign a new state to the variable

        if (family_nlevels.lData[imp_var] < 2) {
            _String oops ("WARNING: Attempting to impute missing values for variable having fewer than two known levels..\n");
            WarnError(oops);

            continue;   // drop this step from the chain
        } else if (family_nlevels.lData[imp_var] == 2) {
            // by convention, levels are zero-indexed; if not, this won't work
            impute._Matrix::Store (0, imp, (_Parameter) (last_state ? 0 : 1));
        } else {    // randomly assign a state (other than the current one) with uniform probability
            pick = genrand_int32() % (family_nlevels.lData[imp_var]-1);

            if (pick >= last_state) {
                pick++;     // skip current state
            }

            impute._Matrix::Store (0, imp_var, (_Parameter)pick);
        }


#ifdef __DEBUG_IDS__
        snprintf (buf, sizeof(buf), "proposed impute: \n");
        BufferToConsole (buf);
        for (long bug = 0; bug < impute.GetUsed(); bug++) {
            snprintf (buf, sizeof(buf), "%d ", (long) impute (0, bug));
            BufferToConsole (buf);

            if (bug % family_size == family_size - 1) {
                NLToConsole ();
            }
        }
        NLToConsole ();
#endif

        // given the case just modified, update m_ijk's.

        child_state = impute(0, imp_case*family_size);

        if (imp_var > 0) {
            // changed a parental variable, update index
            new_index = 0;

            for (long par = 0; par < parents.lLength; par++) {
                long    parent_state    = impute (0, imp_case*family_size + par+1);
                new_index += parent_state * multipliers.lData[par];
            }

#ifdef __DEBUG_IDS__
            snprintf (buf, sizeof(buf), "old_index = %d\nnew_index = %d\n", old_index, new_index);
            BufferToConsole (buf);
#endif
            //  It shouldn't be necessary to re-calculate log_prior_impute from scratch every time
            //  but this doesn't work properly yet :-P afyp

            /*
             log_prior_impute -= LnGamma(n_ijk(old_index,child_state) + DIRICHLET_FLATTENING_CONST
             + m_ijk(old_index, child_state));
             log_prior_impute -= LnGamma(n_ijk(new_index,child_state) + DIRICHLET_FLATTENING_CONST
             + m_ijk(new_index, child_state));
             */
            m_ijk.Store (old_index, child_state, m_ijk (old_index, child_state) - 1);
            m_ijk.Store (new_index, child_state, m_ijk (new_index, child_state) + 1);
            /*
             log_prior_impute += LnGamma(n_ijk(old_index,child_state) + DIRICHLET_FLATTENING_CONST
             + m_ijk(old_index, child_state));
             log_prior_impute += LnGamma(n_ijk(new_index,child_state) + DIRICHLET_FLATTENING_CONST
             + m_ijk(new_index, child_state));
             */

            // non-integer flattening constant prevents us from operating directly on factorial terms
            /*
             log_prior_impute -= log ( n_ijk (old_index, child_state) + m_ijk (old_index, child_state) + DIRICHLET_FLATTENING_CONST - 1 );
             log_prior_impute += log ( n_ijk (new_index, child_state) + m_ijk (new_index, child_state) + DIRICHLET_FLATTENING_CONST );
             */

            /*
             log_prior_impute -= - LnGamma(n_ij(old_index,0) + DIRICHLET_FLATTENING_CONST * r_i + m_ij(old_index,0));
             log_prior_impute -= - LnGamma(n_ij(new_index,0) + DIRICHLET_FLATTENING_CONST * r_i + m_ij(new_index,0));
             */
            m_ij.Store (old_index, 0, m_ij (old_index, 0) - 1);
            m_ij.Store (new_index, 0, m_ij (new_index, 0) + 1);
            /*
             log_prior_impute += - LnGamma(n_ij(old_index,0) + DIRICHLET_FLATTENING_CONST * r_i + m_ij(old_index,0));
             log_prior_impute += - LnGamma(n_ij(new_index,0) + DIRICHLET_FLATTENING_CONST * r_i + m_ij(new_index,0));
             */
        } else {
            /*
             log_prior_impute -= LnGamma(n_ijk(old_index,child_state) + DIRICHLET_FLATTENING_CONST
             + m_ijk(old_index, last_state));
             log_prior_impute -= LnGamma(n_ijk(old_index,child_state) + DIRICHLET_FLATTENING_CONST
             + m_ijk(old_index, last_state));
             */
            m_ijk.Store (old_index, last_state, m_ijk (old_index, last_state) - 1);
            m_ijk.Store (old_index, child_state, m_ijk (old_index, child_state) + 1);
            /*
             log_prior_impute += LnGamma(n_ijk(old_index,child_state) + DIRICHLET_FLATTENING_CONST
             + m_ijk(old_index, child_state));
             log_prior_impute += LnGamma(n_ijk(old_index,child_state) + DIRICHLET_FLATTENING_CONST
             + m_ijk(old_index, child_state));
             */

            /*
             log_prior_impute -= log ( n_ijk (old_index, last_state) + m_ijk (old_index, last_state) - 1);
             log_prior_impute += log ( n_ijk (old_index, child_state) + m_ijk (old_index, child_state) );
             */
            // m_ij's unchanged, parental combination not modified
        }


        // re-compute prior
        last_prior_impute       = log_prior_impute;
        log_prior_impute        = 0.0;
        /*
         for (long j = 0; j < num_parent_combos; j++)
         {
         log_prior_impute += LnGamma(m_ij(j,0)+1) + LnGamma(n_ij(j,0)+DIRICHLET_FLATTENING_CONST*r_i)
         - LnGamma(m_ij(j,0)+n_ij(j,0)+DIRICHLET_FLATTENING_CONST*r_i);
         for (long k = 0; k < r_i; k++)
         {
         log_prior_impute += LnGamma(m_ijk(j,k) + n_ijk(j,k)+DIRICHLET_FLATTENING_CONST)
         - LnGamma(m_ijk(j,k)+1) - LnGamma(n_ijk(j,k)+DIRICHLET_FLATTENING_CONST);
         }
         }
         */
        for (long j = 0; j < num_parent_combos; j++) {
            log_prior_impute += LnGamma(n_ij(j,0) + DIRICHLET_FLATTENING_CONST * r_i)
                                - LnGamma(n_ij(j,0) + DIRICHLET_FLATTENING_CONST * r_i + m_ij(j,0));

            for (long k = 0; k < r_i; k++) {
                log_prior_impute += LnGamma(n_ijk(j,k) + DIRICHLET_FLATTENING_CONST + m_ijk(j,k))
                                    - LnGamma(n_ijk(j,k) + DIRICHLET_FLATTENING_CONST);
            }
        }


#ifdef __DEBUG_IDS__
        snprintf (buf, sizeof(buf), "proposed m_ijk: \n");
        BufferToConsole (buf);
        for (long j = 0; j < num_parent_combos; j++) {
            for (long k = 0; k < family_nlevels.lData[0]; k++) {
                snprintf (buf, sizeof(buf), "%d ", (long) m_ijk (j,k));
                BufferToConsole (buf);
            }
            NLToConsole ();
        }
        NLToConsole ();
#endif


        // compute posterior probability for this step, given N_ijk = n_ijk + m_ijk
        log_score = log_prior_impute;

        for (long j = 0; j < num_parent_combos; j++) {
            log_score += LnGamma(num_levels.lData[node_id]);    // (r-1)!
            log_score -= LnGamma(n_ij(j, 0) + m_ij(j,0) + num_levels.lData[node_id]);   // (N+r-1)!

            for (long k = 0; k < num_levels.lData[node_id]; k++) {
                log_score += LnGamma( n_ijk (j,k) + m_ijk (j,k) + 1);   // (N_ijk)!
            }
        }

        /* I have a hunch this can all be stream-lined... - afyp April 4, 2008 */


        // accept step?     (Metropolis-Hastings)
        lk_ratio    = exp(log_score - last_score);

#ifdef __DEBUG_IDS__
        snprintf (buf, sizeof(buf), "prior = %f, log_score = %f, last_score = %f, lk_ratio = %f\n", log_prior_impute, log_score, last_score, lk_ratio);
        BufferToConsole (buf);
#endif
        if (lk_ratio > 1. || genrand_real2() < lk_ratio) {  // accept proposed imputation
            last_impute = (_Matrix &)impute;
            last_score = log_score;
            last_prior_impute   = log_prior_impute;
            last_m_ijk = m_ijk;
            last_m_ij = m_ij;
#ifdef __DEBUG_IDS__
            snprintf (buf, sizeof(buf), "accept\n\n");
            BufferToConsole (buf);
#endif
        } else {    // revert to original imputation
            for (long i = 0; i < impute.GetUsed(); i++) {
                impute._Matrix::Store (0, i, last_impute (0, i));
            }
            log_score = last_score;
            log_prior_impute    = last_prior_impute;
            m_ijk = last_m_ijk;
            m_ij = last_m_ij;
#ifdef __DEBUG_IDS__
            snprintf (buf, sizeof(buf), "reject\n\n");
            BufferToConsole (buf);
#endif
        }



        // handle sampling of chain
        if (step >= mcem_burnin) {
            if ( (step - (long)mcem_burnin) % (long)mcem_interval == 0) {
                prior_denom->Store (log_prior_impute);  // append to _GrowingVector objects
                post_chain->Store (log_score);
#ifdef __DEBUG_IDS__
                snprintf (buf, sizeof(buf), "%f, %f\n", log_prior_impute, log_score);
                BufferToConsole (buf);
#endif
            }
        }
    }
    // exit MCEM for-loop

    expect_log_score = LogSumExpo (post_chain) - LogSumExpo (prior_denom);



    snprintf (buf, sizeof(buf), "Node %ld, parents(%ld): ", node_id, parents.lLength);
    BufferToConsole (buf);
    for (long par = 0; par < parents.lLength; par++) {
        snprintf (buf, sizeof(buf), " %ld", parents.lData[par]);
        BufferToConsole (buf);
    }


    snprintf (buf, sizeof(buf), "  sum(prior) = %f  E[log score] = %f\n", LogSumExpo (prior_denom), expect_log_score);
    BufferToConsole (buf);

    DeleteObject (prior_denom);
    DeleteObject (post_chain);

    return expect_log_score;
}





//___________________________________________________________________________________________
_Parameter Bgm::ComputeContinuousScore (long node_id)
{
    _SimpleList parents,
                continuous_parents;

    // generate lists of discrete or continuous parents from graph
    for (long par = 0; par < num_nodes; par++) {
        if (dag(par, node_id) == 1) {
            if (is_discrete.lData[par]) {
                parents << par;
            } else {
                continuous_parents << par;
            }
        }
    }

    return ComputeContinuousScore (node_id, parents, continuous_parents);
}


_Parameter Bgm::ComputeContinuousScore (long node_id, _Matrix *g)
{
    _SimpleList     dparents,
                    cparents;

    for (long par = 0; par < num_nodes; par++) {
        // support for discrete parents only
        if ((*g)(par, node_id) == 1) {
            if (is_discrete.lData[par]) {
                dparents << par;
            } else {
                cparents << par;
            }
        }
    }
    return ComputeContinuousScore (node_id, dparents, cparents);
}


//___________________________________________________________________________________________

_Parameter Bgm::ComputeContinuousScore (long node_id, _SimpleList & parents, _SimpleList & continuous_parents)
{
    /* WARNING, untested function! */
    /* --------------------------- */

    // mm = prior estimate of unconditional mean at continuous node (i.e. intercept)
    // phi = scale parameter of inverse gamma prior for variance, uninformative at low values
    //                      (Kohn, Smith, and Chan, 2001 Stat Comput 11: 313-322)

    _Parameter  log_score = 0.;

    long        num_parent_combos = 1,  // i.e. 'q'
                k = continuous_parents.lLength;     // current number of continuous parents

    _Matrix     n_ij,
                pa_indexing,        // track discrete parent combinations per observation
                mu,
                tau;

    _SimpleList     multipliers ((long)1);



    // set location hyperparameter for Gaussian prior
    CreateMatrix (&mu, k+1, 1, false, true, false);
    mu.Store (0, 0, prior_mean (node_id, 0));               // prior intercept
    for (long i = 1; i < mu.GetHDim(); i++) {
        mu.Store (i, 0, 0);     // set prior expectation of regression coefficients to zero
    }



    // how many combinations of parental states are there?
    for (long par = 0; par < parents.lLength; par++) {
        num_parent_combos *= num_levels.lData[parents.lData[par]];
        multipliers << num_parent_combos;
    }

    // set prior degrees of freedom (for inverse gamma / scaled inverse chi-square)
    _Parameter      rho = prior_sample_size (node_id, 0) > 0 ? (prior_sample_size (node_id, 0) / num_parent_combos) : 1.0;


    // set precision hyperparameter for Gaussian prior
    CreateMatrix (&tau, k+1, k+1, false, true, false);
    for (long row = 0; row < tau.GetHDim(); row++) {
        for (long col = 0; col < tau.GetVDim(); col++) {
            if (row == col) {
                tau.Store (row, col, rho);
            } else {
                tau.Store (row, col, 0.);    // zero off-diagonal entries
            }
        }
    }


    // count up number of data points per parent combination
    CreateMatrix (&n_ij, num_parent_combos, 1, false, true, false);
    CreateMatrix (&pa_indexing, obsData->GetHDim(), 1, false, true, false);

    for (long obs = 0; obs < obsData->GetHDim(); obs++) {
        long    index       = 0,
                multiplier = 1;
        //child_state = (*obsData)(obs, node_id);

        for (long par = 0; par < parents.lLength; par++) {
            long    this_parent     = parents.lData[par];
            // max index = sum (parent * levels(parent))
            index += (*obsData)(obs, this_parent) * multiplier;
            multiplier *= num_levels.lData[this_parent];
        }

        pa_indexing.Store (obs, 0, index);
        n_ij.Store (index, 0, n_ij(index, 0) + 1);
    }



    // for every parent combination, calculate contribution to score
    for (long pa = 0; pa < num_parent_combos; pa++) {
        _Parameter  pa_log_score = 0.;

        _Matrix     zbpa (n_ij(pa, 0), continuous_parents.lLength + 1, false, true),
                    yb (n_ij(pa, 0), 1, false, true),
                    scale;

        long        count_n     = 0;
        // number of data points with this parent combination.


        // populate zbpa matrix with (1, y_1, y_2, ..., y_m) entries for this parental combo
        for (long obs = 0; obs < obsData->GetHDim(); obs++) {
            if (pa_indexing(obs, 0) == pa) {    // this observation has the current parent combo
                // load observed states at continuous parents into matrix
                //   I'm sure there's a faster way to do this! - afyp

                zbpa.Store (count_n, 0, 1);     // intercept

                for (long parent = 0; parent < continuous_parents.lLength; parent++) {
                    zbpa.Store (count_n, parent+1, (*obsData)(obs, continuous_parents.lData[parent]));
                }

                // state vector for this continuous node
                yb.Store (count_n, 0, (*obsData)(obs, node_id));

                count_n++;
            }
        }

        // calculate scale parameter for non-central t distribution
        // from S. Bottcher (2001) p.25
        scale = zbpa;
        scale *=  * (_Matrix *) tau.Inverse();
        zbpa.Transpose();
        scale *= zbpa;
        zbpa.Transpose();
        for (long row = 0; row < scale.GetHDim(); row++) {  // add identity matrix
            scale.Store (row, row, scale(row, row)+(_Parameter)1.);
        }
        scale *= (_Parameter) (prior_precision (node_id, 0) / rho);


        // calculate the determinant of scale parameter matrix
        _Matrix         temp_mat (scale);
        _Parameter      pi_const = 3.141592653589793;

        temp_mat *= (_Parameter) (pi_const * rho);

        _AssociativeList *  eigen       = (_AssociativeList *) temp_mat.Eigensystem();
        _Matrix *           eigenvalues = (_Matrix *)eigen->GetByKey(0, MATRIX);
        _Parameter          det         = 1.;

        // determinant is product of eigenvalues (should be > 0 for positive definite matrices)
        for (long i = 0; i < eigenvalues->GetHDim(); i++) {
            det *= (_Parameter)(*eigenvalues) (i,0);
        }


        // calculate first term of score
        pa_log_score += LnGamma((rho + n_ij(pa, 0))/2.);
        pa_log_score -= LnGamma(rho/2.) + 0.5 * log(det);


        // calculate second term of score
        _Matrix     next_mat;
        zbpa *= mu;
        yb -= zbpa;
        temp_mat = yb;
        next_mat = temp_mat;
        next_mat *= * (_Matrix *) scale.Inverse();
        temp_mat.Transpose();
        next_mat *= temp_mat;
        next_mat *= (_Parameter) (1./prior_precision (node_id, 0)); // should be a 1-element matrix

        pa_log_score += -(rho + n_ij(pa,0))/2. * (next_mat(0,0) + 1.);
        log_score += pa_log_score;
    }


    return log_score;
}



//__________________________________________________________________________________________//
//                                                                                          //
//                              DYNAMIC BAYESIAN NETWORKS                                   //
//__________________________________________________________________________________________//
//                                                                                          //
//      A(t)    B(t)        "2 parents" means  A(t+1) <-- A(t) + B(t) + B(t+1)              //
//       |     / |                                                                          //
//       |   /   |          i.e., there are never actually two parents..                    //
//       V L     V                                                                          //
//   A(t+1) <-- B(t+1)                                                                      //
//__________________________________________________________________________________________//

_DynamicBgm::_DynamicBgm (_AssociativeList * dnodes, _AssociativeList * cnodes) : Bgm (dnodes, cnodes)
{
    // call base class constructor
    if (num_nodes % 2 > 0) {
        _String errorMsg ("ERROR: Expected even number of nodes in network.");
        WarnError (errorMsg);
    }

    // check that variables retain type across time slices, i.e. a discrete node doesn't become continuous
    for (long node = 0; node < num_nodes; node += 2) {
        if (is_discrete.lData[node] != is_discrete.lData[node+1]) {
            _String errorMsg ("ERROR: Inconsistent node type across time slices.");
            WarnError (errorMsg);
        }
    }
}

//___________________________________________________________________________________________
#ifdef __NEVER_DEFINED__
_Parameter  _DynamicBgm::Compute (_SimpleList * node_order, _List * results)
{
    _Parameter          log_likel   = 0.;
    _GrowingVector      *gv1, *gv2;

    for (long i = 0; i < num_nodes * num_nodes; i++) {
        gv1 = (_GrowingVector *) results->lData[i];
        gv1 -> ZeroUsed();
    }

    for (long nodeIndex = 0; nodeIndex < node_order->lLength; nodeIndex++) {
        long                child_node      = node_order->lData[nodeIndex],
                            maxp            = max_parents.lData[child_node];

        _List           *   score_lists     = (_List *) node_scores.lData[child_node];

        gv1 = (_GrowingVector *) results->lData[child_node * num_nodes + child_node];
        gv1->ZeroUsed();

        if (child_node % 2 == 1) {
            _Constant       *   orphan_score    = (_Constant *) (score_lists->lData[0]);
            gv1 -> Store (orphan_score->Value());
        } else {
            _Matrix *   single_parent_scores    = (_Matrix *) (score_lists->lData[1]);
        }

        long        self_parent             = node_order->lData[nodeIndex + 1];


        gv1 -> Store ((*single_parent_scores) (self_parent, 0));

        gv2 = (_GrowingVector *) results->lData[child_node * num_nodes + self_parent];
        gv2 -> Store ((*single_parent_scores) (self_parent, 0));

        if (maxp > 1) {
            _SimpleList         parents,
                                eligible_parents;
            _Parameter          tuple_score;
            _NTupleStorage *    family_scores;

            for (long parIndex = nodeIndex + 3; parIndex < node_order->lLength; parIndex = parIndex + 2) {
                long    par = node_order->lData[parIndex];

                if (banned_edges (par, child_node) == 0) {
                    eligible_parents << par;
                }
            }

            family_scores = (_NTupleStorage *) (score_lists->lData[2]);
            parents.Populate (2, 0, 0);
            for (long ep = 0; ep < eligible_parents.lLength; ep++) {
                long    this_parent = eligible_parents.lData[ep];
                if (this_parent > child_node) {
                    parents.lData[0] = self_parent-1;
                    parents.lData[1] = this_parent-1;
                } else {
                    parents.lData[0] = this_parent;
                    parents.lData[1] = self_parent-1;
                }
                tuple_score = family_scores -> Retrieve (parents);

                gv1 -> Store (tuple_score);

                gv2 = (_GrowingVector *) results->lData[child_node * num_nodes + self_parent];
                gv2 -> Store (tuple_score);
                gv2 = (_GrowingVector *) results->lData[child_node * num_nodes + this_parent];
                gv2 -> Store (tuple_score);
            }


            if (maxp > 2) {
                _SimpleList         indices (eligible_parents.lLength, 0, 1);

                for (long nparents = 3; nparents <= maxp; nparents++) {
                    _SimpleList     subset,
                                    auxil;
                    bool            not_finished;

                    if (nparents-1 > eligible_parents.lLength) {
                        break;
                    }

                    if (indices.NChooseKInit (auxil, subset, nparents-1, false)) {
                        parents.Populate(nparents, 0, 0);
                        parents.lData[0] = self_parent;

                        family_scores = (_NTupleStorage *) (score_lists->lData[nparents]);

                        do {
                            not_finished = indices.NChooseK (auxil, subset);

                            for (long i = 0; i < nparents-1; i++) {
                                long    realized = eligible_parents.lData[subset.lData[i]];
                                if (realized >= child_node) {
                                    realized--;
                                }
                                parents.lData[i] = realized;
                            }
                            parents.Sort(TRUE);

                            tuple_score = family_scores -> Retrieve (parents);
                            gv1 -> Store (tuple_score);
                            gv2 = (_GrowingVector *) results -> lData [child_node*num_nodes + self_parent];
                            gv2 -> Store (tuple_score);

                            for (long i = 0; i < nparents-1; i++) {
                                gv2 = (_GrowingVector *) results->lData[child_node * num_nodes + eligible_parents.lData[subset.lData[i]]];
                                gv2 -> Store (tuple_score);
                            }
                        } while (not_finished);
                    }
                }
            }
        }
    }
    gv1 -> _Matrix::Store (0, 0, LogSumExpo(gv1));
    log_likel += (*gv1)(0, 0);
}
return log_likel;
}
#endif

#ifdef __NEVER_DEFINED__
//___________________________________________________________________________________________
void    _DynamicBgm::SetDataMatrix (_Matrix * data)
{
    scores_cached   = FALSE;

    if (obsData)    {
        DeleteObject (obsData);    // dispense with previous data matrix
    }
    obsData         = data;
    obsData->CheckIfSparseEnough(TRUE);


    // reset data-dependent member variables
    for (long node = 0; node < num_nodes; node++) {
        has_missing.lData[node] = 0;
        num_levels.lData[node] = 0;
    }


    // check for missing data and compute levels
    if (obsData->GetVDim() == 2*num_nodes) {
        long    nrows = obsData->GetHDim(),
                obs,
                node;

        for (long col = 0; col < obsData->GetVDim(); col++) {
            node = col / 2;

            if (is_discrete.lData[node]) {
                if (col % 2 == 0) {
                    num_levels.lData[node] = 1;    // reset for new variable
                }

                for (long row = 0; row < nrows; row++) {
                    obs = (*obsData)(row,col);

                    if (has_missing[node] == 0 && obs < 0) {
                        has_missing.lData[node] = 1;
                        continue;   // skip next step to check levels
                    }

                    obs++;  // adjust for zero-indexing
                    if (obs > num_levels.lData[node]) {
                        num_levels.lData[node]++;
                    }
                }
            } else {
                // continuous node defined to have no levels
                num_levels.lData[node] = 0;
            }
        }
    } else {
        _String errorMsg ("Number of variables in data matrix do not match number of nodes in graph.");
        WarnError (errorMsg);       // note, this produces a crash because the batch file proceeds to execute
        // BGM routines without data. - AFYP
    }


    last_node_order.Clear();        // forget the last step taken in MCMC chain
    best_node_order.Clear();

    CacheNodeScores();
}
#endif


//___________________________________________________________________________________________
void    _DynamicBgm::CacheNodeScores (void)
{
    //  derived CacheNodeScores() imposes dBGM structural constraints when computing
    //  node scores, which are then re-indexed to group like nodes across time slices,
    //  halving [num_nodes].
    //  will have to reset member variables.
    if (scores_cached) {
        return;
    }


    _SimpleList     parents, empty_list;
    _Parameter      score;
    long            dbn_nodes   = num_nodes / 2;

    for (long node_id = 0; node_id < num_nodes; node_id++) {
        long        maxp        = max_parents.lData[node_id];

        _List   *   this_list   = (_List *) node_scores.lData[node_id]; // retrieve pointer to list of scores for this node
        this_list->Clear();

        parents.Clear();

        if (node_id % 2 == 0) { // cache only even-numbered nodes in current time slice
            // will always have at least one parent (node+1)
            parents << node_id + 1;
            score = is_discrete.lData[node_id]  ? ComputeDiscreteScore (node_id, parents)
                    : ComputeContinuousScore (node_id);
            _Constant   orphan_score (score);
            (*this_list) && (&orphan_score);


            // 1 parent --> 3 DBN-parents
            _Matrix     single_parent_scores (dbn_nodes, 1, false, true);

            parents.Populate (3, 0, 0);
            for (long par = 0; par < num_nodes; par = par + 2) {    // traverse even-nodes
                if (par == node_id) {
                    continue;
                }

                if (par < node_id) {    // parents need to be in numerical order
                    parents.lData[0] = par;
                    parents.lData[1] = par+1;
                    parents.lData[2] = node_id+1;
                } else {
                    parents.lData[0] = node_id+1;
                    parents.lData[1] = par;
                    parents.lData[2] = par+1;
                }

                single_parent_scores.Store (par/2, 0, is_discrete.lData[node_id]    ? ComputeDiscreteScore (node_id, parents)
                                            : ComputeContinuousScore (node_id));
            }

            (*this_list) && (&single_parent_scores);


            // more than one other parent
            if (maxp > 1) {
                _SimpleList     all_but_one (dbn_nodes-1, 0, 1),    // group nodes into time-slice pairs
                                aux_list,
                                nk_tuple;

                for (long unused, np = 2; np <= maxp; np++) {
                    _NTupleStorage  family_scores (dbn_nodes-1, np);

                    parents.Populate (2*np + 1, 0, 0);
                    parents.lData[0] = node_id+1;   // always depends on itself!

                    if (all_but_one.NChooseKInit (aux_list, nk_tuple, np, false)) {
                        bool    remaining = TRUE;
                        do {
                            remaining = all_but_one.NChooseK (aux_list, nk_tuple);
                            for (long par, par_idx = 0; par_idx < np; par_idx++) {
                                par = 2 * nk_tuple.lData[par_idx];  // convert back to num_nodes index
                                if (par >= node_id) {
                                    par += 2;
                                }
                                parents.lData[2*par_idx+1] = par;
                                parents.lData[2*par_idx+2] = par + 1;
                            }
                            parents.Sort();
                            score = is_discrete.lData[node_id]  ? ComputeDiscreteScore (node_id, parents)
                                    : ComputeContinuousScore (node_id);
                            unused = family_scores.Store (score, nk_tuple);
                            parents.lData[0] = node_id+1;   // required because of prior sort
                        } while (remaining);
                    } else {
                        _String oops ("Failed to initialize _NTupleStorage object in Bgm::CacheNodeScores().\n");
                        WarnError(oops);
                    }
                    (*this_list) && (&family_scores);   // append duplicate to storage
                }
            }
        }
    }

    scores_cached = TRUE;
    CollapseDynamicGraph(); // rescale graph so that A(t)->A(t+1) is treated as a single node
}



//___________________________________________________________________________________________
void    _DynamicBgm::CollapseDynamicGraph (void)
{
    if (scores_cached) {
        _String oops ("Attempted to call _DynamicBgm::CollapseDynamicGraph() before CacheNodeScores().\n");
        WarnError(oops);
    } else if (!obsData) {
        _String oops ("No data matrix assigned to _DynamicBgm object in CollapseDynamicGraph().\n");
        WarnError(oops);
    } else if (num_nodes < obsData->GetVDim()) {
        _String oops ("In CollapseDynamicGraph(): Node space already collapsed!\n");
        WarnError(oops);
    } else {
        long            dbn_nodes   = num_nodes / 2;

        // overwrite entries in class member _Matrix and _SimpleList objects
        for (long row = 0; row < num_nodes; row += 2) {
            for (long col = 0; col < num_nodes; col += 2) {
                /*
                if (banned_edges(row,col) == 1 || banned_edges);    // (P-/->C)
                {

                }
                 */
            }
        }
        num_nodes  = dbn_nodes;
    }
}

//___________________________________________________________________________________________

/*  WARNING:  DO NOT activate these debugging flags unless the number of iterations is very LOW */

//#define   __DEBUG_GADS__
//#define __DEBUG_GADS2__
_Parameter  Bgm::GibbsApproximateDiscreteScore (long node_id, _SimpleList & parents)
{
    //char          buf [255];
    ReportWarning (_String("Imputing missing values for node ") & node_id & " with parents " & (_String *) parents.toStr());

    long            num_parent_combos   = 1,
                    r_i                   = num_levels.lData[node_id],
                    max_num_levels        = r_i,
                    family_size           = parents.lLength + 1;

    _SimpleList     multipliers ((long) 1), // length constructor, populator
                    family_nlevels (family_size, 0, 0),
                    is_missing;             // store linear indices to missing entries

    is_missing.RequestSpace((long) family_size * obsData->GetHDim());

    _Matrix         n_ijk, n_ij,
                    m_ijk, m_ij,
                    data_deep_copy,         // make a deep copy of the relevant columns of the data matrix
                    reassign_probs;

    _GrowingVector  * vector_of_scores  = new _GrowingVector();     // for storing scores sampled during re-assignment of missing values

    _Parameter      log_score           = 0,
                    max_iterations,
                    parent_state, child_state,
                    denom,
                    this_prob;

    double          urn;            // uniform random number


    // get MCEM settings
    checkParameter (mcemMaxSteps, max_iterations, 0.);
    if (max_iterations < 1) {
        _String oops ("Problem in BGM::GibbsApproximateDiscreteScore(): Did you forget to specify a value for BGM_MCEM_MAXSTEPS?\n");
        WarnError (oops);
    }



    // to index into N_ij/k matrices by parent combination
    for (long par = 0; par < parents.lLength; par++) {
        num_parent_combos *= num_levels.lData[parents.lData[par]];
        multipliers << num_parent_combos;
    }


    // get number of levels per family member
    family_nlevels.lData[0] = r_i;
    for (long par = 0; par < parents.lLength; par++) {
        family_nlevels.lData[par+1] = num_levels.lData [ parents.lData[par] ];
        if (family_nlevels.lData[par+1] > max_num_levels) {
            max_num_levels = family_nlevels.lData[par+1];
        }
    }



    // allocate space to matrices
    CreateMatrix (&n_ijk, num_parent_combos, r_i, false, true, false);
    CreateMatrix (&n_ij, num_parent_combos, 1, false, true, false);
    CreateMatrix (&data_deep_copy, obsData->GetHDim(), family_size, false, true, false);
    CreateMatrix (&reassign_probs, max_num_levels, 1, false, true, false);



    // make deep copy and annotate which entries are missing - for use during Gibbs sampling
    // simultaneously tally N_ijk's for complete cases
    _SimpleList is_complete ((unsigned long) data_deep_copy.GetHDim());
    long        n_complete_cases = 0;

    for (long pa_index, row = 0; row < obsData->GetHDim(); row++) {
        pa_index                = 0;
        child_state             = (*obsData) (row, node_id);
        data_deep_copy.Store (row, 0, child_state);
        is_complete.lData[row]  = 1;

        if (child_state < 0) {
            is_missing << row * family_size;
            is_complete.lData[row] = 0;
        }

        for (long par = 0; par < parents.lLength; par++) {
            parent_state = (*obsData) (row, parents.lData[par]);
            data_deep_copy.Store (row, par+1, parent_state);

            if (parent_state < 0) {
                is_missing << row * family_size + par+1;
                is_complete.lData[row] = 0;
            } else {
                pa_index += ((long)parent_state) * multipliers.lData[par];
            }
        }

        if (is_complete.lData[row]) {
            n_complete_cases += 1;
            n_ijk.Store (pa_index, (long) child_state, n_ijk(pa_index, (long) child_state) + 1);
            n_ij.Store  (pa_index, 0, n_ij(pa_index, 0) + 1);
        }
    }


    // convert N_ijk's into probability matrix - equation 16 from Cooper and Herskovits
    // reset N_ijk and N_ij while we're looping
    _Matrix prob_n_ijk, prob_n_ij;

    CreateMatrix (&prob_n_ijk, num_parent_combos, r_i, false, true, false);
    CreateMatrix (&prob_n_ij, num_parent_combos, 1, false, true, false);

    for (long j = 0; j < num_parent_combos; j++) {
        prob_n_ij.Store ( j, 0, (n_ij(j,0) + 1) / (n_complete_cases + num_parent_combos) ); // is this the correct estimator?

        for (long k = 0; k < r_i; k++) {
            prob_n_ijk.Store( j, k, prob_n_ij(j,0) * (n_ijk(j,k) + 1) / (n_ij(j,0) + r_i) );
            n_ijk.Store( j, k, 0.);
        }
        n_ij.Store ( j, 0, 0.);
    }


    // cache marginalizations over this probability matrix as determined by which variables are missing
    // this would be nice, but let's brute-force it for now :-P


    // initialize missing entries to random assignments based on complete cases (N_ijk)
    if (family_size > 1) {
        for (long row = 0; row < data_deep_copy.GetHDim(); row++) {
            if (!is_complete.lData[row]) {  // this row contains a missing value
                _Parameter  denominator = 0.;
                _SimpleList keep_pa_indices;

                // determine which parental state combinations (indexed by j) are possible
                // store those indices in a list and accumulate the probabilities in a denominator
                for (long j = 0; j < num_parent_combos; j++) {
                    bool keep = TRUE;
                    for (long par = 1; par <= parents.lLength; par++) {
                        if ( data_deep_copy(row,par) >= 0 ) {
                            if ( (long)(j/multipliers.lData[par-1]) % num_levels.lData[parents.lData[par-1]] != (long)data_deep_copy(row,par) ) {
                                // observed parent state is incompatible with this pa combo
                                keep = FALSE;
                                break;
                            }
                        }
                    }

                    if (keep) {
                        keep_pa_indices << j;

                        long    k = (long) data_deep_copy(row,0);
                        if ( k >= 0 ) {
                            denominator += prob_n_ijk (j, k);   // child state is observed, append from this column only
                        } else {
                            for (k = 0; k < r_i; k++) {
                                denominator += prob_n_ijk (j, k);    // sum across columns (unknown child state)
                            }
                        }
                    }
                }


                // select random values (j, k) subject to constraints
                long    k       = (long) data_deep_copy(row,0),
                        pick_k,
                        pick_j;

                urn = genrand_real2();

                for (long pa = 0; pa < keep_pa_indices.lLength; pa++) {
                    pick_j = keep_pa_indices.lData[pa];

                    if (k >= 0 && k < r_i) {    /* k is a fixed non-missing value */
                        if ( urn < (prob_n_ijk (pick_j,k) / denominator) ) {    /*** LEAK ***/
                            break;
                        }
                        urn -= prob_n_ijk (pick_j,k) / denominator; /*** LEAK ***/
                    } else if (k == -1) {   // child is missing
                        for (pick_k = 0; pick_k < r_i; pick_k++) {
                            if ( urn < (prob_n_ijk (pick_j,pick_k) / denominator) ) {
                                break;
                            }
                            urn -= prob_n_ijk (pick_j,pick_k) / denominator;
                        }
                        if (pick_k < r_i) {
                            // break out of outer loop carrying current values [j] and [pick_k]
                            break;
                        }
                        // otherwise go to top of outer loop to try a different [j]
                    } else {
                        WarnError(_String("This shouldn't happen! : row = ") & row & "  k = " & k & "  data_deep_copy(row,0) = " & data_deep_copy(row,0) & "  (long)data_deep_copy(row,0) = " & (long)data_deep_copy(row,0));
                    }
                }

                //ReportWarning (_String ("j,k=") & j & "," & k);

                // convert j, k to data
                if ( k < 0 ) {
                    data_deep_copy.Store(row, 0, pick_k);
                }

                for (long par = 1; par <= parents.lLength; par++) {
                    if ( data_deep_copy(row,par) < 0 ) {
                        data_deep_copy.Store(row, par, (long)(pick_j/multipliers.lData[par-1]) % num_levels.lData[parents.lData[par-1]]);
                    }
                }

            }
        }
    } else {    // family_size == 1
        for (long row = 0; row < data_deep_copy.GetHDim(); row++) {
            if (data_deep_copy(row,0) < 0) {
                urn = genrand_real2();
                for (child_state = 0; child_state < r_i; child_state++) {
                    this_prob = prob_n_ijk(0,(long)child_state);
                    if (urn < this_prob) {
                        break;
                    }
                    urn -= this_prob;
                }

                data_deep_copy.Store (row, 0, child_state);
            }
        }
    }


    //ReportWarning (_String("Initialized data : ") & (_String *) data_deep_copy.toStr());


    // recount N_ij and N_ijk from entire data set (with imputed cases)
    n_complete_cases = data_deep_copy.GetHDim();    // variable re-use!
    for (long pa_index, row = 0; row < n_complete_cases; row++) {
        pa_index        = 0;
        child_state     = data_deep_copy (row, 0);

        for (long par = 0; par < parents.lLength; par++) {
            parent_state = data_deep_copy (row, par+1);
            pa_index += ((long)parent_state) * multipliers.lData[par];
        }

        n_ijk.Store (pa_index, (long) child_state, n_ijk(pa_index, (long) child_state) + 1);
        n_ij.Store  (pa_index, 0, n_ij(pa_index, 0) + 1);
    }


    //ReportWarning (_String("N_ijk:") & (_String *) n_ijk.toStr());


    //  Cycle through random permutation of missing entries and reassign by marginal probability.
    //      Based on family configuration of the corresponding case, adjust N_ijk to compute reassignment probs.
    //      If accept reassignment, update N_ijk values.

    // log_score = K2Score (r_i, n_ij, n_ijk);  // compute score from complete cases

    for (long iter = 0; iter < max_iterations; iter++) {
        // shuffle list of missing entries
        // is_missing.Permute(1);
        // ReportWarning (_String("is_missing Permute:") & (_String *) is_missing.toStr());

        for (long row, col, pa_index, missing_idx = 0; missing_idx < is_missing.lLength; missing_idx++) {
            row         = is_missing.lData[missing_idx] / family_size;
            col         = is_missing.lData[missing_idx] % family_size;
            pa_index    = 0;

            // determine parent combination for this row
            for (long par = 0; par < parents.lLength; par++) {
                pa_index += ((long) data_deep_copy (row, par+1)) * multipliers.lData[par];
            }


            if (col == 0) { // child node, reassign N_ijk's, keep N_ij's constant
                child_state = data_deep_copy (row, col);
                n_ijk.Store (pa_index, (long) child_state, n_ijk(pa_index, (long) child_state) - 1);

                // compute probabilities for all instantiations of child node
                denom = 0.;
                for (long k = 0; k < r_i; k++) {
                    denom += this_prob = (n_ijk(pa_index,k) + 1) / (n_ij(pa_index,0) + r_i);
                    reassign_probs.Store (k, 0, this_prob);
                }


                // random reassignment
                urn = genrand_real2 ();
                for (long lev = 0; lev < r_i; lev++) {
                    this_prob = reassign_probs(lev,0) / denom;
                    //ReportWarning(_String("child: urn, this_prob = ") & urn & "," & this_prob);
                    if ( urn < this_prob ) {
                        n_ijk.Store (pa_index, lev, n_ijk(pa_index,lev) + 1);
                        data_deep_copy.Store (row, col, lev);
                        break;
                    }
                    urn -= this_prob;
                }
            } else {    // col > 0, parent node, reassign both N_ij's AND N_ijk's
                parent_state    = data_deep_copy (row, col);
                child_state     = data_deep_copy (row, 0);

                n_ij.Store  (pa_index, 0, n_ij(pa_index,0) - 1);
                n_ijk.Store (pa_index, (long) child_state, n_ijk(pa_index, (long) child_state) - 1);

                pa_index    -= ((long)parent_state) * multipliers.lData[col-1]; // marginalize index
                denom       = 0.;
                for (long pa_temp, lev = 0; lev < family_nlevels.lData[col]; lev++) {
                    // parental index for this instantiation of parent node given other parents (fixed)
                    pa_temp = pa_index + lev * multipliers.lData[col-1];

                    // by Bayes' Theorem, Pr(P|C,P') is proportional to Pr(C|P,P') x Pr(P,P')
                    denom += this_prob = (n_ijk(pa_temp,(long)child_state) + 1) / (n_ij(pa_temp,0) + r_i)
                                         * (n_ij(pa_temp,0)+1) / (n_complete_cases + num_parent_combos);
                    reassign_probs.Store (lev, 0, this_prob);
                }


                // re-assign parent state
                urn = genrand_real2();  // a uniform random number within interval [0,1)

                for (long lev = 0; lev < family_nlevels.lData[col]; lev++) {
                    this_prob = reassign_probs (lev,0) / denom;
                    //ReportWarning(_String("parent: urn, this_prob = ") & urn & "," & this_prob);
                    if (urn < this_prob) {
                        pa_index += lev * multipliers.lData[col-1]; // note that pa_index was decremented above

                        n_ij.Store (pa_index, 0, n_ij(pa_index,0) + 1);
                        n_ijk.Store (pa_index, (long) child_state, n_ijk(pa_index, (long) child_state) + 1);
                        data_deep_copy.Store (row, col, lev);
                        break;
                    }

                    urn -= this_prob;
                }

            }
            // end if-else
            //ReportWarning (_String ("N_ijk at step ") & missing_idx & " : " & (_String *) n_ijk.toStr() );
        }
        // end loop over missing values
        vector_of_scores->Store ( K2Score(r_i, n_ij, n_ijk) );
    }
    // end loop over iterations

    /*
    ReportWarning (_String("Score sample: "));
    for (long s = 0; s < vector_of_scores->GetUsed(); s++)
    {
        ReportWarning (_String(" ") & (* (_Matrix *)vector_of_scores) (s,0));
    }
    */

    // compute the average of sampled log-likelihoods
    log_score = LogSumExpo (vector_of_scores) - log(max_iterations);

    DeleteObject (vector_of_scores);

    ReportWarning (_String("Returning expected log score: ") & log_score);

    return (log_score);
}



//___________________________________________________________________________________________
_Parameter Bgm::K2Score (long r_i, _Matrix & n_ij, _Matrix & n_ijk)
{
    _Parameter log_score = 0.;

    for (long j = 0; j < n_ij.GetHDim(); j++) {
        log_score += LnGamma(r_i);  // (r-1)!
        log_score -= LnGamma(n_ij(j, 0) + r_i); // (N+r-1)!

        for (long k = 0; k < r_i; k++) {
            log_score += LnGamma(n_ijk(j,k) + 1);   // (N_ijk)!
        }
    }

    return log_score;
}


//___________________________________________________________________________________________
//  Later we'll need to expand this function to take arbitrary alpha_ijk values.
/*
_Parameter Bgm::BDeScore (long r_i, _Matrix & n_ij, _Matrix & n_ijk)
{
    _Parameter  n_prior_ij      = prior_sample_size (node_id, 0) / num_parent_combos,
                n_prior_ijk     = n_prior_ij / num_levels.lData[node_id],
                log_score       = 0.;

    for (long j = 0; j < num_parent_combos; j++)
    {
        log_score += LnGamma(n_prior_ij) - LnGamma(n_prior_ij + n_ij(j,0));

        for (long k = 0; k < num_levels.lData[node_id]; k++)
        {
            log_score += LnGamma(n_prior_ijk + n_ijk(j,k)) - LnGamma(n_prior_ijk);
        }
    }

    return log_score;
}
*/

//___________________________________________________________________________________________
//___________________________________________________________________________________________
//___________________________________________________________________________________________
//___________________________________________________________________________________________
//___________________________________________________________________________________________
#endif

