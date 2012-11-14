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


/*SLKP 20070926; include progress report updates */
#if !defined __UNIX__ && !defined __HEADLESS__
#include "HYConsoleWindow.h"
#endif
/*SLKP*/


#if defined __AFYP_DEVELOPMENT__ && defined __HYPHYMPI__
#include "mpi.h"
#endif



/* 
	Quantity to add to all entries in contingency table 
	to avoid zeroes that would cause an error when taking
	natural log.
	see Fienberg and Holland (1972) J Multivar Anal 2: 127-134.
*/
#define     DIRICHLET_FLATTENING_CONST  0.5
#define		MIN_SAMPLE_SIZE				5


class _BayesianGraphicalModel : public _LikelihoodFunction
{
public:
    /* constructors */
    _BayesianGraphicalModel () { }
    _BayesianGraphicalModel (_AssociativeList *);

    /* destructor */
    virtual ~_BayesianGraphicalModel (void);


    /* network initialization */
    bool            SetDataMatrix   (_Matrix *),    // via SetParameter HBL
                    SetWeightMatrix (_Matrix *),
                    SetConstraints    (_Matrix *),    //  "       "
                    SetStructure  (_Matrix *),
                    SetParameters (_AssociativeList *),
                    SetNodeOrder  (_SimpleList *);


    /* computation */
    virtual _Parameter      Compute (void);	// compute likelihood of current network structure
    _Parameter      Compute (_Matrix &),	// compute likelihood of given network structure
                    Compute (_SimpleList &, _List *);	// compute likelihood of given node order
                    									//	return edge marginal probabilities at pointer
                    									
    virtual _Matrix *       Optimize ();	// generic wrapper from HBL to different optimization methods
    										// e.g., K2, structural MCMC, order MCMC (see next functions)
	
    void            GraphMetropolis (bool, long, long, long, _Parameter, _Matrix *),
                    OrderMetropolis (bool, long, long, _Parameter, _Matrix *),
                    K2Search (bool, long, long, _Matrix *);


    void            CacheNodeScores (void);	// MPI enabled
    void            MPIReceiveScores (_Matrix *, bool, long);
    void            ReleaseCache (void);

    _Parameter      ComputeDiscreteScore (long node_id),
                    ComputeDiscreteScore (long, _Matrix &),
                    ComputeDiscreteScore (long, _SimpleList &),

                    ComputeContinuousScore (long node_id),
                    ComputeContinuousScore (long, _Matrix &),
                    ComputeContinuousScore (long, _SimpleList &);


    _Parameter      ImputeDiscreteNodeScore (long, _SimpleList &),	// use Gibbs sampling to compute expectation over missing data
    				ImputeCGNodeScore (long, _SimpleList &);		// arguments: node ID, parent ID's

    void            ComputeParameters (void),	// UNDER DEVELOPMENT - and I think I ended up using HBL instead
                    ComputeParameters (_Matrix *);




    /* input/output */
    void                SerializeBGM (_String &);	// export network structure and parameters to HBL
    bool                ImportModel (_AssociativeList *),	// THIS HAS NOT BEEN WRITTEN
                        ExportCache (_AssociativeList *),	// send node score cache to HBL as associative list
                        ImportCache (_AssociativeList *);	// set node score cache to HBL associative list


    /* utility */
    void            InitMarginalVectors (_List *);
    void            DumpMarginalVectors (_List *);

    void            SerializeBGMtoMPI (_String &);	// pass network object to compute node as HBL

    void            RandomizeGraph (_Matrix *, _SimpleList *, _Parameter, long, long, bool);
    _SimpleList *   GetOrderFromGraph (_Matrix &);
    bool            GraphObeysOrder (_Matrix &, _SimpleList &);

    void            UpdateDirichletHyperparameters (long , _SimpleList &, _Matrix * , _Matrix * );

    _Parameter      K2Score (long, _Matrix &, _Matrix &),
                    BDeScore (long, _Matrix &, _Matrix &),
                    BottcherScore (_Matrix &, _Matrix &, _Matrix &, _Matrix &, _Parameter, _Parameter, long);

    long            GetNumNodes (void)  {
        return num_nodes;
    }
    long            GetNumCases (void)  {
        return theData.GetHDim();
    }

    void            GetNodeOrder (_Matrix * order);
    void            GetStructure (_Matrix * graph);
    void            GetConstraints (_Matrix * graph) {
        graph = (_Matrix *) constraint_graph.makeDynamic();
    }

protected:

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

    _Parameter      continuous_missing_value;       // some arbitrary value set in HBL to indicate just that

    /* ------------------------------------------- */

    _Matrix         theStructure;

    _Matrix         constraint_graph;   // integer, 0 = no constraint, -1 = banned edge, 1 = enforced edge

    _List           node_score_cache;
    bool            scores_cached;

    _SimpleList     node_order_arg;     // provides access to node ordering functionality as HBL argument

    /* ------------------------------------------- */

};



//______________________________________________________________________________________________
#ifdef __NEVER_DEFINED__
class _DynamicBayesGraph : public _BayesianGraphicalModel
{
public:

protected:

};
#endif


