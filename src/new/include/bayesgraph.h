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


#ifndef __BAYES_GRAPH_H__

#define __BAYES_GRAPH_H__
#include "list.h"
#include "simplelist.h"
#include "classes.h"
#include "likefunc.h"
#include "parser.h"
#include <math.h>
#include "matrix.h"
#include "baseobj.h"
#include "batchlan.h"
// #include "HYUtils.h"



#if defined __AFYP_DEVELOPMENT__ && defined __HYPHYMPI__
#include "mpi.h"
#endif



/* 
	Quantity to add to all entries in contingency table 
	to avoid zeroes that would cause an error when taking
	natural log.
	see Fienberg and Holland (1972) J Multivar Anal 2: 127-134.
*/
#define   kBGMDirichletFlattening  0.5
#define		kBGMMinSize       				5


class _BayesianGraphicalModel : public _LikelihoodFunction {
public:
    /* constructors */
    _BayesianGraphicalModel () { }
    _BayesianGraphicalModel (_AssociativeList *);

    /* destructor */
    virtual ~_BayesianGraphicalModel (void);


    /* network initialization */
    bool            SetDataMatrix   (_Matrix const *),    // via SetParameter HBL
                    SetWeightMatrix (_Matrix const *),
                    SetConstraints    (_Matrix const *),    //  "       "
                    SetStructure  (_Matrix const *),
                    SetParameters (_AssociativeList const *),
                    SetNodeOrder  (_SimpleList const *);


    /* computation */
    virtual hyFloat      Compute (void);	// compute likelihood of current network structure
    hyFloat              Compute (_Matrix const &),	// compute likelihood of given network structure
                         Compute (_SimpleList const &, _List *);	// compute likelihood of given node order
                    									//	return edge marginal probabilities at pointer
                    									
    virtual _Matrix *       Optimize (_AssociativeList const * options = nil);	// generic wrapper from HBL to different optimization methods
    										// e.g., K2, structural MCMC, order MCMC (see next functions)
	
    _Matrix*        GraphMetropolis (bool, long, long, long, hyFloat);
    _Matrix*        OrderMetropolis (bool, long, long, hyFloat);
    _Matrix*        K2Search (bool, long, long);


    void            CacheNodeScores (void);	// MPI enabled
    void            MPIReceiveScores (_Matrix *, bool, long);
    void            ReleaseCache (void);

    hyFloat         ComputeDiscreteScore (long node_id),
                    ComputeDiscreteScore (long, _Matrix const &),
                    ComputeDiscreteScore (long, _SimpleList const &),

                    ComputeContinuousScore (long node_id),
                    ComputeContinuousScore (long, _Matrix const &),
                    ComputeContinuousScore (long, _SimpleList const &);


    hyFloat         ImputeDiscreteNodeScore (long, _SimpleList const &) const,	// use Gibbs sampling to compute expectation over missing data
                    ImputeCGNodeScore (long, _SimpleList const &) const;		// arguments: node ID, parent ID's

    void            ComputeParameters (void),	// UNDER DEVELOPMENT - and I think I ended up using HBL instead
                    ComputeParameters (_Matrix *);




    /* input/output */
    void                SerializeBGM (_StringBuffer &);	// export network structure and parameters to HBL
    bool                ImportModel (_AssociativeList *),	// THIS HAS NOT BEEN WRITTEN
                        ExportCache (_AssociativeList *) const,	// send node score cache to HBL as associative list
                        ImportCache (_AssociativeList *);	// set node score cache to HBL associative list


    /* utility */
    void            InitMarginalVectors (_List *) const;
    void            DumpMarginalVectors (_List *) const;

    void            SerializeBGMtoMPI (_StringBuffer &);	// pass network object to compute node as HBL

    void            RandomizeGraph (_Matrix *, _SimpleList *, hyFloat, long, long, bool);
    _SimpleList const   GetOrderFromGraph (_Matrix const &) const;
    bool            GraphObeysOrder (_Matrix &, _SimpleList const &) const;

    void            UpdateDirichletHyperparameters (long , _SimpleList const &, _Matrix * , _Matrix * );

    hyFloat          K2Score (long, _Matrix const &, _Matrix const &) const,
                    BDeScore (long, _Matrix const&, _Matrix const&) const,
                    BottcherScore (_Matrix const &, _Matrix const &, _Matrix const &, _Matrix const &, hyFloat, hyFloat, long) const;

    long            GetNumNodes (void)  {
        return num_nodes;
    }
    long            GetNumCases (void)  {
        return theData.GetHDim();
    }

    void            GetNodeOrder (_Matrix * order) const;
    void            GetStructure (_Matrix * graph) const;
    _Matrix*        GetConstraints (void) const {
        return (_Matrix*)constraint_graph.makeDynamic();
      
    }

protected:

    bool            is_node_continuous (long node) const {return node_type.get (node) == 1L;}
    bool            is_node_discrete   (long node) const {return node_type.get (node) == 0L;}
    const _SimpleList
                    ComputeParentMultipliers (_SimpleList const&) const;
    long            ReindexParentObservations(_SimpleList const& parents, _SimpleList& n_ij, _SimpleList& pa_indexing) const;
  
    long            num_nodes;

    /* ------------------------------------------- */

    _Matrix         theData,
                    theWeights;

    _List           node_names;     // list of strings

    _SimpleList     node_type,      // boolean, 0 = discrete, 1 = continuous
                    num_levels,     // integer, if discrete, number of levels
                    max_parents,    // integer, maximum number of parents
                    has_missing;    // boolean, 0 = complete data, 1 = missing, (2 = latent, i.e., all missing)

    _Matrix         prior_sample_size,
                    prior_mean,         // for continuous (Gaussian) nodes
                    prior_precision,
                    prior_scale;

    hyFloat      continuous_missing_value;       // some arbitrary value set in HBL to indicate just that

    /* ------------------------------------------- */

    _Matrix         theStructure;

    _Matrix         constraint_graph;   // integer, 0 = no constraint, -1 = banned edge, 1 = enforced edge

    _List           node_score_cache;
    bool            scores_cached;

    _SimpleList     node_order_arg;     // provides access to node ordering functionality as HBL argument

    /* ------------------------------------------- */
  
    static          hyFloat      LogSumExpo (_Vector * log_values);

};

#endif

