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
 
 _AVLList	structure inspired by the excellent documentation of 
 GNU libavl 2.0.1 by Ben Pfaff (http://www.msu.edu/~pfaffben/avl/index.html)
								
								*/

/*
 *  Bgm.cpp
 *  HYPHY_XCode
 *
 *  Created by Art Poon on 2/5/07.
 *	Based on source code written in C by Fraser Iain Lewis.
 *
 */

#define		__OMIT_MISSING_DATA__		1

#define		LOG_SCALING_FACTOR			64.
#define		LARGE_NEGATIVE_NUMBER		-99999.0
#define		MAX_LSE_SCALING_ATTEMPTS	4096

#include "bgm.h"

/*SLKP 20070926; include progress report updates */
	#if !defined __UNIX__ && !defined __HEADLESS__
			#include "HYConsoleWindow.h"
	#endif
/*SLKP*/

_String		_HYBgm_BAN_PARENT_KEY	("BanParent"),
			_HYBgm_BAN_CHILD_KEY	("BanChild"),
			
			_HYBgm_NODE_INDEX_KEY	("NodeID"),
			_HYBgm_PRIOR_SIZE_KEY	("PriorSize"),
			_HYBgm_MAX_PARENT_KEY	("MaxParents"),

			_HYBgm_PRIOR_MEAN_KEY	("PriorMean"),		/* for continuous (Gaussian) nodes */
			_HYBgm_PRIOR_PREC_KEY	("PriorVar"),
			
/*SLKP 20070926; add string constants for progress report updates */
			_HYBgm_STATUS_LINE_MCMC			("Running Bgm MCMC"),
			_HYBgm_STATUS_LINE_MCMC_DONE	("Finished Bgm MCMC"),
/*SLKP*/
			/* maxParentString ("Bgm_MAXIMUM_PARENTS"), */
			maxNumRestart	("BGM_NUM_RESTARTS"),
			numRandomize	("BGM_NUM_RANDOMIZE"),
			
			mcmcNumChains	("BGM_MCMC_NCHAINS"),		// re-use parameter
			mcmcTemperature ("BGM_MCMC_TEMPERATURE"),
			mcmcSteps		("BGM_MCMC_DURATION"),
			mcmcBurnin		("BGM_MCMC_BURNIN"),
			mcmcSamples		("BGM_MCMC_SAMPLES"),

			isDynamicGraph	("BGM_DYNAMIC"),

			_HYBgm_NThreads	("BGM_NTHREADS");


#ifdef __MP__
	#ifndef __MACHACKMP__
		#include <pthread.h>
	#else
		#include "mypthread.h"
	#endif
	struct	ThreadCacheTask {
		Bgm		* b;
		
		long	startNode,
				nextNode,	// one past the last node
				numNodes;
		
		_List	*		nodeScores;
		
		_SimpleList		maxParents,		// pass member variables to global scope
						isDiscrete;
	};
	
	pthread_t *			BgmThreads	= nil;
	ThreadCacheTask *	BgmTasks	= nil;
	void *				CacheNodeScoreThread (void *);
#endif


#ifdef 		__UNIX__

void		ConsoleBGMStatus (_String, _Parameter, _String * fileName = nil);

//__________________________________________________________________________________

void		ConsoleBGMStatus (_String statusLine, _Parameter percentDone, _String * fileName = nil)
{
	FILE		   *outFile = fileName?doFileOpen (fileName->sData,"w"):nil;
	
	_String		   reportLine (statusLine);
	
	

	
	if (percentDone >= 0.0)
		reportLine = reportLine & ". " & percentDone & "% done.";

	if (outFile)
		fprintf (outFile,"%s", reportLine.sData);			
	else
		if (verbosityLevel == 1)
			printf ("\033\015 %s", reportLine.sData);

	if (percentDone < -1.5)
	{
		printf ("\033\015 ");
		setvbuf (stdout,nil,_IOLBF,1024);
	}
	else
		if (percentDone < -0.5)
			setvbuf (stdout,nil, _IONBF,1);
	if (outFile)
		fclose (outFile);

}

#endif

//___________________________________________________________________________________________

long integerPower (long base, long exponent)
{
	long	result = 1,
			mask   = 1<<(sizeof(long)*8-2);	// left shift to left-most position of binary sequence for long integer
											// e.g. 100...0 (30 zeroes for signed long)
	
	while ((exponent & mask) == 0) mask >>= 1;	// bitwise AND, right-shift mask until overlaps with first '1'
	
	while (mask)
	{
		result *= result;
		if (exponent & mask)
			result = result * base;
		mask >>= 1;
	}
	return result;
}



//___________________________________________________________________________________________

Bgm::Bgm(_AssociativeList * dnodes, _AssociativeList * cnodes, _AssociativeList * banlist)
{	
	char		buf [255];
	_String		errorMessage;
	
	long		num_discrete	= dnodes->avl.countitems(),
				num_continuous	= cnodes->avl.countitems(),
				node_index;
	
	_Constant	* node_id,						/* pointers into HBL associative array arguments */
				* mp,
				* size, * mean, * precision,
				* banParent, * banChild;
	
	_AssociativeList	* discreteNode, * continuousNode, * bannedEdge;
	
	
	
	// set member variables to default values
	calc_bit			= new _Constant ();			// for calculating LnGamma in member functions
	max_max_parents		= 0;						// record the most parents a node can have
	obsData				= NULL;
	obsWeights			= NULL;
	
	
	// is this a dynamic Bayesian network?
	_Parameter		dynamicArg;
	checkParameter (isDynamicGraph, dynamicArg, 0.);
	if (dynamicArg > 0)
		is_dynamic_graph = TRUE;
	else
		is_dynamic_graph = FALSE;
	
	
	// extract number of variables in data
	num_nodes = num_discrete + num_continuous;
	if (num_nodes < 2)
		WarnError (_String("ERROR: Attempting to construct Bayesian network on fewer than two variables."));
	
	//sprintf (buf, "Initializing Bgm with %d nodes...\n", num_nodes);
	//BufferToConsole (buf);	
	
	
	CreateMatrix (&dag, num_nodes, num_nodes, false, true, false);		// allocate space for matrices
	CreateMatrix (&banned_edges, num_nodes, num_nodes, false, true, false);
	
	CreateMatrix (&prior_sample_size, num_nodes, 1, false, true, false);
	CreateMatrix (&prior_mean, num_nodes, 1, false, true, false);
	CreateMatrix (&prior_precision, num_nodes, 1, false, true, false);
	
	
	is_discrete.Populate (num_nodes, 0, 0);		// allocate space for _SimpleList objects
	num_levels.Populate (num_nodes, 0, 0);
	max_parents.Populate (num_nodes, 0, 0);
	
	// node_scores.Populate (num_nodes, 0, 0);
	
	
	// convert banlist to matrix form
	for (long banitem = 0; banitem < banlist->avl.countitems(); banitem++)
	{
		bannedEdge	= (_AssociativeList *) (banlist->GetByKey (banitem, ASSOCIATIVE_LIST));
		
		if (bannedEdge)
		{
			banParent	= (_Constant *) (bannedEdge->GetByKey (_HYBgm_BAN_PARENT_KEY, NUMBER));
			banChild	= (_Constant *) (bannedEdge->GetByKey (_HYBgm_BAN_CHILD_KEY, NUMBER));
			
			if (banParent && banChild)
			{
				banned_edges.Store ((long)banParent->Value(), (long)banChild->Value(), 1.);
			}
			else
			{
				errorMessage = _String("Banned edge ") & banitem & " must specify both a parent and a child node.";
				break;
			}
		}
	}
	
	
	// initialize discrete nodes
	for (long dn = 0; dn < num_discrete; dn++)
	{
		discreteNode	= (_AssociativeList *) (dnodes->GetByKey (dn, ASSOCIATIVE_LIST));
		
		if (discreteNode)
		{
			node_id		= (_Constant *) (discreteNode->GetByKey (_HYBgm_NODE_INDEX_KEY, NUMBER));
			mp			= (_Constant *) (discreteNode->GetByKey (_HYBgm_MAX_PARENT_KEY, NUMBER));
			size		= (_Constant *) (discreteNode->GetByKey (_HYBgm_PRIOR_SIZE_KEY, NUMBER));
			
			if (node_id && mp && size)
			{
				node_index	= (long) (node_id->Value());
				
				is_discrete.lData[node_index] = 1.;
				max_parents.lData[node_index] = (long) mp->Value();
				
				if ((long) mp->Value() > max_max_parents)
					max_max_parents = (long) mp->Value();
				
				prior_sample_size.Store(node_index, 0, (long) size->Value());
			}
			else
			{
				errorMessage = _String ("Missing key (NodeID, MaxParents, PriorSize) in associative array for discrete node ") 
								& dn;
				break;
			}
		}
		else
		{
			errorMessage = _String("Failed to retrieve discrete node specification at index ") & dn;
			break;
		}
	}
	
	if (errorMessage.sLength)
	{
		WarnError (errorMessage);
		errorMessage = _String();	// reset to empty string
	}
	
	
	
	// initialize continuous nodes
	for (long cn = 0; cn < num_continuous; cn++)
	{
		continuousNode	= (_AssociativeList *) (cnodes->GetByKey (cn, ASSOCIATIVE_LIST));
		
		if (continuousNode)
		{
			node_id		= (_Constant *) (continuousNode->GetByKey (_HYBgm_NODE_INDEX_KEY, NUMBER));
			mp			= (_Constant *) (continuousNode->GetByKey (_HYBgm_MAX_PARENT_KEY, NUMBER));
			size		= (_Constant *) (continuousNode->GetByKey (_HYBgm_PRIOR_SIZE_KEY, NUMBER));
			mean		= (_Constant *) (continuousNode->GetByKey (_HYBgm_PRIOR_MEAN_KEY, NUMBER));
			precision	= (_Constant *) (continuousNode->GetByKey (_HYBgm_PRIOR_PREC_KEY, NUMBER));
			
			if (node_id && mp && size && mean && precision)
			{
				node_index = (long) (node_id->Value());
				
				is_discrete.lData[node_index] = 0.;
				max_parents.lData[node_index] = (long) mp->Value();
				
				if (mp->Value() > max_max_parents)
					max_max_parents = (long) mp->Value();
				
				prior_sample_size.Store (node_index, 0, (long) size->Value());
				prior_mean.Store (node_index, 0, (long) mean->Value());
				prior_precision.Store (node_index, 0, (long) precision->Value());
				
			}
			else
			{
				errorMessage = _String ("Missing key (NodeID, MaxParents, PriorSize, PriorMean, PriorVar)")
										& "in associative array for continuous node " & cn;
				break;
			}
		}
		else
		{
			errorMessage = _String ("Failed to retrieve continuous node specification at index ") & cn;
			break;
		}
	}
	
	if (errorMessage.sLength)
	{
		WarnError (errorMessage);
		errorMessage = _String();	// reset to empty string
	}
	
	
	
	// append pointer to _List object that will store floats or _NTupleStorage objects
	_List		emptyList (max_max_parents+1);
	
	for (long node = 0; node < num_nodes; node++)
	{
		node_scores && (&emptyList);
	}
	
	
	// require Bgm to recalculate node scores
	scores_cached = FALSE;
	
}



//___________________________________________________________________________________________
Bgm::~Bgm (void)
{
	DeleteObject (calc_bit);
}



//___________________________________________________________________________________________
void Bgm::SetDataMatrix (_Matrix * data)
{
	// reset cached node scores and edge posteriors
	// DumpComputeLists ();
	scores_cached	= FALSE;
	
	if (obsData)	DeleteObject (obsData);		// dispense with previous data matrix
	
	obsData			= data;
	
	obsData->CheckIfSparseEnough(TRUE);
	
#ifdef DEBUG_SDM
	char buf [256];
	for (long i = 0; i < obsData->GetHDim(); i++)
	{
		for (long j = 0; j < obsData->GetVDim(); j++)
		{
			sprintf (buf, "%d ", (long) (*obsData)(i,j));
			BufferToConsole (buf);
		}
		NLToConsole ();
	}
#endif
	
	if (obsData->GetVDim() == num_nodes)
	{
		for (long node = 0; node < num_nodes; node++)
		{
			if (is_discrete.lData[node])
			{
				// calculate number of levels for discrete node
				num_levels.lData[node] = 1;
				
				for (long obs = 0; obs < obsData->GetHDim(); obs++)
				{
					// adjust for zero-indexing
					if ( ((*obsData)(obs,node) + 1) > num_levels.lData[node])
						num_levels.lData[node] = num_levels.lData[node] + 1;
				}
			}
			else
			{
				// continuous node defined to have no levels
				num_levels.lData[node] = 0;
			}
		}
	}
	else
	{
		_String errorMsg ("Number of variables in data matrix do not match number of nodes in graph.");
		WarnError (errorMsg);
	}
	
#ifdef DEBUG_SDM
	sprintf (buf, "Levels: ");
	BufferToConsole (buf);
	for (long i = 0; i < num_nodes; i++)
	{
		sprintf (buf, "%d ", num_levels.lData[i]);
		BufferToConsole (buf);
	}
	NLToConsole ();
#endif
	
	CacheNodeScores();
	
	// re-allocate memory to lists
	// InitComputeLists ();
}



//___________________________________________________________________________________________
void	Bgm::SetWeightMatrix (_Matrix * weights)
{
	if (obsData && obsData->GetHDim() == weights->GetHDim())
	{
		obsWeights = weights;
	}
	else
	{
		_String errorMsg ("Number of weights does not match number of observations in current data set.");
		WarnError (errorMsg);
	}
	
	scores_cached = FALSE;
}



//___________________________________________________________________________________________
void Bgm::PrintDag (void)
{
	char	buf [255];
	/* for DEBUGGING */
	
	for (long row = 0; row < dag.GetHDim(); row++)
	{
		for (long col = 0; col < dag.GetVDim(); col++)
		{
			sprintf (buf, "%d ", (long)dag(row,col));
			BufferToConsole (buf);
		}
		sprintf (buf, "\n");
		BufferToConsole (buf);
	}
	sprintf (buf, "\n");
	BufferToConsole(buf);
}


//___________________________________________________________________________________________
inline void	Bgm::ResetDag (void)
{
	for (long row = 0; row < num_nodes; row++)
		for (long col = 0; col < num_nodes; col++)
			dag.Store (row,col,0.);
}



//___________________________________________________________________________________________
long	Bgm::MarginalSum (bool over_rows, long offset)
{
	long	res		= 0;
	
	if (over_rows)
		for (long i = 0; i < num_nodes; i++) res += dag(offset,i);	// sum of n-th row
	else
		for (long i = 0; i < num_nodes; i++) res += dag(i,offset);	// sum of n-th column
	
	return res;
}



//___________________________________________________________________________________________
bool	Bgm::IsCyclic (void)
{
	_Matrix		nodes_left (num_nodes, 1, false, true);
	
	for (long i = 0; i < num_nodes; i++)		// populate vector with 1's
		nodes_left.Store (i,0,1.);
	
	long		one_count = num_nodes,
				last_count,
				num_children;
	
	do	// until no more leaves are removed
	{
		last_count = one_count;
		one_count = 0;
		
		// find and remove all nodes with no children (i.e. leaf nodes)
		for (long parent = 0; parent < num_nodes; parent++)
		{
			if (nodes_left(parent,0) == 1.)
			{
				one_count++;
				
				num_children = 0;
				// count all existing nodes that have this parent
				for (long child = 0; child < num_nodes; child++)
					if (nodes_left(child,0) == 1. && child != parent)
						if (dag(parent,child) == 1.)
							num_children++;
				
				if (num_children == 0)		// count number of nodes with this node as parent
				{
					nodes_left.Store (parent, 0, 0.);
					one_count--;
				}
			}
		}
	} while (one_count < last_count);
	
	
	if (one_count > 0) 
	{
		return TRUE;	// dag is cyclic
	}
	else
	{
		return FALSE;	// removed all nodes as leaves, dag is acyclic
	}
}



//___________________________________________________________________________________________
void	Bgm::RandomizeDag (long num_steps)
{	
	
	for (long step = 0; step < num_steps; step++)
	{
		// store current graph in case altered graph becomes cyclic
		_Matrix reset_dag (num_nodes, num_nodes, false, true);
		for (long h = 0; h < num_nodes; h++)
			for (long v = 0; v < num_nodes; v++)
				reset_dag.Store (h,v,dag(h,v));
		
		do
		{
			{
			for (long row = 0; row < num_nodes; row++)
				for (long col = 0; col < num_nodes; col++)
					dag.Store (row,col,reset_dag(row,col));
			}
			// PrintDag();
			
			
			// select a random arc, discarding cycles and banned edges
			long	row		= 0,
					col		= row;
			
			while (col == row && banned_edges(row,col) == 1)
			{
				row = genrand_int32() % num_nodes;
				col = genrand_int32() % num_nodes;
			}
			
			// printf ("RandomizeDag() picked %d,%d\n", row, col);
			
			if (dag(row,col) == 0.)
			{
				if (dag(col,row) == 0.)	// no arc
				{
					// add an arc -- sum over n-th column gives number of parents for n-th child
					if (MarginalSum(FALSE,col) < max_parents.lData[col])
						dag.Store (row,col,1.);
				}
				else
				{
					// reverse or remove arc
					dag.Store (col,row,0.);
					if (genrand_int32() % 2 == 0)
						dag.Store (row,col,1.);
				}
			}
			else
			{
				if (dag(col,row) == 0.)
				{
					// reverse or remove
					dag.Store (row,col,0.);
					if (genrand_int32() % 2 == 0)
						dag.Store (col,row,1.);
				}
				else
				{
					// double arc, set to one of three legal arcs
					long	option = genrand_int32() % 3;
					if (option == 0)
						dag.Store (row,col,0.);
					else if (option == 1)
						dag.Store (col,row,0.);
					else
					{
						dag.Store (row,col,0.);
						dag.Store (col,row,0.);
					}
				}
			}
			
			// PrintDag();
			
		} while ( IsCyclic() );
	}
	// end loop over steps
}





//___________________________________________________________________________________________

inline _Parameter Bgm::LnGamma(_Constant * calculator, long x)
{
	// wrapper function for _Constant member function
	calculator->SetValue (x);
	calculator = (_Constant *) calculator->LnGamma();
	_Parameter rv = calculator->Value();
	DeleteObject(calculator);
	return rv;
}



//___________________________________________________________________________________________
//	Wrapper to retain original functionality.
_Parameter	Bgm::ComputeDiscreteScore (long node_id)
{
	_SimpleList		parents;
	
	for (long par = 0; par < num_nodes; par++)
	{
		if (dag(par, node_id) == 1 && is_discrete.lData[par])
			parents << par;
	}
	
	return ComputeDiscreteScore (node_id, parents);
}



//___________________________________________________________________________________________
_Parameter	Bgm::ComputeDiscreteScore (long node_id, _SimpleList & parents)
{
	_SimpleList		multipliers ((long)1);
	
	
	// generate list of discrete parents from graph
	/*
	for (long par = 0; par < num_nodes; par++)
	{
		if (dag(par, node_id) == 1 && is_discrete.lData[par])
			parents << par;
	}
	 */
	
	
	// use cached node scores if possible
	
	if (scores_cached)
	{
		_List *		scores	= (_List *) node_scores.lData[node_id];
		
		if (parents.lLength == 0)
		{
			_Constant *		orphan_score = (_Constant *) scores->lData[0];
			return (_Parameter) orphan_score->Value();
		}
		else if (parents.lLength == 1)
		{
			_Matrix *	single_parent_scores = (_Matrix *) scores->lData[1];
			return (_Parameter) (*single_parent_scores) (parents.lData[0], 0);
		}
		else
		{
			_NTupleStorage *	family_scores	= (_NTupleStorage *) scores->lData[parents.lLength];
			_SimpleList			nktuple;
			
			for (long i = 0; i < parents.lLength; i++)		// map parents back to nk-tuple
			{
				long	par = parents.lData[i];
				if (par > node_id) par--;
				nktuple << par;
			}
			return (_Parameter) family_scores->Retrieve (nktuple);	// using nk-tuple
		}
	}
	
	
	
	// else need to compute de novo
	_Matrix		n_ijk,
				n_ij;		// [i] indexes child nodes, [j] indexes combinations of values for parents of i-th node, [k] indexes values of i-th node
	
	long		num_parent_combos = 1;	// i.e. 'q'
	
	_Parameter	n_prior_ijk	= 0,
				n_prior_ij	= 0,
				log_score	= 0;
	
	
	// how many combinations of parental states are there?
	for (long par = 0; par < parents.lLength; par++)
	{
		num_parent_combos *= num_levels.lData[parents.lData[par]];
		multipliers << num_parent_combos;		// unused!
	}
	
	
	/* count observations by parent combination, using direct indexing */
	CreateMatrix (&n_ijk, num_parent_combos, num_levels.lData[node_id], false, true, false);
	CreateMatrix (&n_ij, num_parent_combos, 1, false, true, false);
	
	for (long obs = 0; obs < obsData->GetHDim(); obs++)
	{
		long	index		= 0,
				multiplier	= 1,
				child_state = (*obsData)(obs, node_id);
#ifdef __OMIT_MISSING_DATA__
		bool	skip_obs	= FALSE;
		
		if (child_state < 0)	// encode missing data using negative integer values
			continue;
#endif
		
		for (long par = 0; par < parents.lLength; par++)
		{
			long	this_parent			= parents.lData[par],
					this_parent_state	= (*obsData)(obs, this_parent);
			
#ifdef __OMIT_MISSING_DATA__
			if (this_parent_state < 0)
			{
				skip_obs = TRUE;
				break;
			}
#endif
			
			index += this_parent_state * multiplier;
			multiplier *= num_levels.lData[this_parent];
		}
		
#ifdef __OMIT_MISSING_DATA__
		if (!skip_obs)
#endif
		{
			n_ijk.Store ((long) index, child_state, n_ijk(index, child_state) + 1);
			n_ij.Store ((long) index, 0, n_ij(index, 0) + 1);
		}
	}
	
	
	
	/* compute scoring metric */
	if (prior_sample_size (node_id, 0) == 0)
	{
		/* assume no prior information, use K2 metric */
		for (long j = 0; j < num_parent_combos; j++)
		{
			log_score += LnGamma(calc_bit, num_levels.lData[node_id]);	// (r-1)!
			log_score -= LnGamma(calc_bit, (long) n_ij(j, 0) + num_levels.lData[node_id]);	// (N+r-1)!
			
			for (long k = 0; k < num_levels.lData[node_id]; k++)
				log_score += LnGamma (calc_bit, ((long)n_ijk(j,k) + 1));	// (N_ijk)!
			
#ifdef BGM_DEBUG_CDS
			char buf [256];
			
			sprintf (buf, "Node %d, parent(s) ", node_id);
			BufferToConsole (buf);
			
			for (long k = 0; k < parents.lLength; k++)
			{
				sprintf (buf, "%d ", parents.lData[k]);
				BufferToConsole (buf);
			}
			NLToConsole();
			
			sprintf (buf, "\tlog(r-1)! = log %d! = %lf\n", num_levels.lData[node_id] - 1, LnGamma(calc_bit, num_levels.lData[node_id]));
			BufferToConsole (buf);
			
			sprintf (buf, "\tlog(N+r-1)! = log %d! = %lf\n", (long) n_ij(j, 0) + num_levels.lData[node_id] - 1, 
					LnGamma(calc_bit, (long) n_ij(j, 0) + num_levels.lData[node_id]));
			BufferToConsole (buf);
			
			for (long k = 0; k < num_levels.lData[node_id]; k++)
			{
				sprintf (buf, "\tlog (N_ijk)! = log %d! = %lf\n", ((long)n_ijk(j,k)), LnGamma (calc_bit, ((long)n_ijk(j,k) + 1)));
				BufferToConsole (buf);
			}
			
			sprintf (buf, "j = %d\tcumulative log score = %lf\n", j, log_score);
			BufferToConsole (buf);
#endif
		}
	} 
	else 
	{
		/* calculate Bayesian Dirichlet metric (BDeu) for this node */
		/* see p212 in Heckerman, Geiger, and Chickering (1995) Machine Learning 20, 197-243 */
		n_prior_ij = prior_sample_size (node_id, 0) / num_parent_combos;
		n_prior_ijk = n_prior_ij / num_levels.lData[node_id];
		
		for (long j = 0; j < num_parent_combos; j++)
		{
			log_score += LnGamma (calc_bit, n_prior_ij) - LnGamma (calc_bit, n_prior_ij + n_ij(j,0));
			
			for (long k = 0; k < num_levels.lData[node_id]; k++)
				log_score += LnGamma (calc_bit, n_prior_ijk + n_ijk(j,k)) - LnGamma (calc_bit, n_prior_ijk);
		}
	}
#ifdef BGM_DEBUG_CDS
	char	buf[255];
	sprintf (buf, "Log score = %f\n\n", log_score);
	BufferToConsole (buf);
#endif
	
	return (log_score);
}



//___________________________________________________________________________________________

_Parameter Bgm::ComputeContinuousScore (long node_id)
{
	/* WARNING, untested function! */
	/* --------------------------- */
	
	// mm = prior estimate of unconditional mean at continuous node (i.e. intercept)
	// phi = scale parameter of inverse gamma prior for variance, uninformative at low values 
	//						(Kohn, Smith, and Chan, 2001 Stat Comput 11: 313-322)
	
	_Parameter	log_score = 0.;
	
	long		num_parent_combos = 1;	// i.e. 'q'
	
	_Matrix		n_ij,
				pa_indexing,		// track discrete parent combinations per observation
				mu,
				tau;
	
	_SimpleList		parents,
					continuous_parents,
					multipliers ((long)1);
	
	
	
	// generate lists of discrete or continuous parents from graph
	for (long par = 0; par < num_nodes; par++)
	{
		if (dag(par, node_id) == 1)
		{
			if (is_discrete.lData[par])
				parents << par;
			else
				continuous_parents << par;
		}
	}
	
	// current number of continuous parents
	long	k = continuous_parents.lLength;
	
	
	// set location hyperparameter for Gaussian prior
	CreateMatrix (&mu, k+1, 1, false, true, false);
	mu.Store (0, 0, prior_mean (node_id, 0));				// prior intercept
	for (long i = 1; i < mu.GetHDim(); i++)
	{
		mu.Store (i, 0, 0);		// set prior expectation of regression coefficients to zero
	}
	
	
	
	// how many combinations of parental states are there?
	{
	for (long par = 0; par < parents.lLength; par++)
	{
		num_parent_combos *= num_levels.lData[parents.lData[par]];
		multipliers << num_parent_combos;
	}
	}
	
	// set prior degrees of freedom (for inverse gamma / scaled inverse chi-square)
	_Parameter		rho	= prior_sample_size (node_id, 0) > 0 ? (prior_sample_size (node_id, 0) / num_parent_combos) : 1.0;
	
	
	// set precision hyperparameter for Gaussian prior
	CreateMatrix (&tau, k+1, k+1, false, true, false);
	for (long row = 0; row < tau.GetHDim(); row++)
	{
		for (long col = 0; col < tau.GetVDim(); col++)
		{
			if (row == col) tau.Store (row, col, rho);
			else tau.Store (row, col, 0.);	// zero off-diagonal entries
		}
	}
	
	
	// count up number of data points per parent combination
	CreateMatrix (&n_ij, num_parent_combos, 1, false, true, false);
	CreateMatrix (&pa_indexing, obsData->GetHDim(), 1, false, true, false);
	for (long obs = 0; obs < obsData->GetHDim(); obs++)
	{
		long	index		= 0,
				multiplier	= 1,
				child_state = (*obsData)(obs, node_id);
		
		for (long par = 0; par < parents.lLength; par++)
		{
			long	this_parent		= parents.lData[par];
			// max index = sum (parent * levels(parent))
			index += (*obsData)(obs, this_parent) * multiplier;
			multiplier *= num_levels.lData[this_parent];
		}
		
		pa_indexing.Store (obs, 0, index);
		n_ij.Store (index, 0, n_ij(index, 0) + 1);
	}
	
	
	
	// for every parent combination, calculate contribution to score
	for (long pa = 0; pa < num_parent_combos; pa++)
	{
		_Parameter	pa_log_score = 0.;
		
		_Matrix		zbpa (n_ij(pa, 0), continuous_parents.lLength + 1, false, true),
					yb (n_ij(pa, 0), 1, false, true),
					scale;
		
		long		count_n		= 0;	
								// number of data points with this parent combination.
		
		
		// populate zbpa matrix with (1, y_1, y_2, ..., y_m) entries for this parental combo
		for (long obs = 0; obs < obsData->GetHDim(); obs++)
		{
			if (pa_indexing(obs, 0) == pa)		// this observation has the current parent combo
			{									// load observed states at continuous parents into matrix
												//   I'm sure there's a faster way to do this! - afyp
				
				zbpa.Store (count_n, 0, 1);		// intercept
				
				for (long parent = 0; parent < continuous_parents.lLength; parent++)
					zbpa.Store (count_n, parent+1, (*obsData)(obs, continuous_parents.lData[parent]));
				
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
		for (long row = 0; row < scale.GetHDim(); row++)	// add identity matrix
		{
			scale.Store (row, row, scale(row, row)+(_Parameter)1.);
		}
		scale *= (_Parameter) (prior_precision (node_id, 0) / rho);
		
		
		// calculate the determinant of scale parameter matrix
		_Matrix			temp_mat (scale);
		_Parameter		pi_const = 3.141592653589793;
		
		temp_mat *= (_Parameter) (pi_const * rho);
		
		_AssociativeList *	eigen		= (_AssociativeList *) temp_mat.Eigensystem();
		_Matrix *			eigenvalues = (_Matrix *)eigen->GetByKey(0, MATRIX);
		_Parameter			det			= 1.;
		
		// determinant is product of eigenvalues (should be > 0 for positive definite matrices)
		for (long i = 0; i < eigenvalues->GetHDim(); i++)
			det *= (_Parameter)(*eigenvalues) (i,0);
		
		
		// calculate first term of score
		pa_log_score += LnGamma (calc_bit, (rho + n_ij(pa, 0))/2.);
		pa_log_score -= LnGamma (calc_bit, rho/2.) + 0.5 * log(det);
		
		
		// calculate second term of score
		_Matrix		next_mat;
		zbpa *= mu;
		yb -= zbpa;
		temp_mat = yb;
		next_mat = temp_mat;
		next_mat *= * (_Matrix *) scale.Inverse();
		temp_mat.Transpose();
		next_mat *= temp_mat;
		next_mat *= (_Parameter) (1./prior_precision (node_id, 0));	// should be a 1-element matrix
		
		pa_log_score += -(rho + n_ij(pa,0))/2. * (next_mat(0,0) + 1.);
		log_score += pa_log_score;
	}
	
	
	return log_score;
}



//___________________________________________________________________________________________
//  THIS DOES NOT WORK YET --- need to provide global scope for member variable node_scores, 
//	through accessor function maybe?
#ifdef __NOT_DEFINED_MP__
void * CacheNodeScoreThread (void * arg)
{
	ThreadCacheTask	*	theTask	= (ThreadCacheTask *) arg;
	
	long				num_nodes	= theTask->numNodes;
	_SimpleList			is_discrete	= theTask->isDiscrete;
	
	for (long node_id = theTask->startNode; node_id < theTask->nextNode; node_id++)
	{
		long		maxp		= theTask->maxParents.lData[node_id];
		_List	*	this_list	= (_List *) (theTask->nodeScores)->lData[node_id];
		
		
		this_list->Clear();
		
		/*
		for (long par = 0; par < num_nodes; par++)	// reset local graph at this node
		{
			dag.Store (par, node_id, 0);
		}
		*/
		
		// handle case of no parents
		_SimpleList	parents;
		_Parameter	score = is_discrete.lData[node_id] ? theTask->b->ComputeDiscreteScore (node_id, parents) : theTask->b->ComputeContinuousScore (node_id);
		_Constant	orphan_score (score);
		
		(*this_list) && (&orphan_score);		// _List Append() specific to objects in BaseObj hierarchy
		
		
		char buf [256];
		sprintf (buf, "orphan score = %f\n", score);
		BufferToConsole (buf);
		
		
		// handle case of one parent
		if (maxp > 0)
		{
			_Matrix		single_parent_scores (num_nodes, 1, false, true);
			
			for (long par = 0; par < num_nodes; par++)
			{
				if (par == node_id)		// unused
					continue;
				
				//dag.Store (par, node_id, 1);
				parents << par;
				single_parent_scores.Store (par, 0, is_discrete.lData[node_id] ? theTask->b->ComputeDiscreteScore (node_id, parents) : 
											theTask->b->ComputeContinuousScore (node_id));
				parents.Clear();
				//dag.Store (par, node_id, 0);	// reset
			}
			(*this_list) && (&single_parent_scores);
		}
		
		
		// handle cases of more than one parent using (n,k)-tuple indexing
		if (maxp > 1)
		{
			_SimpleList		all_but_one (num_nodes-1, 0, 1),	// 0, 1, 2, ... , n-1
							aux_list,
							nk_tuple;
			
			for (long np = 2; np <= maxp; np++)
			{
				_NTupleStorage	family_scores (num_nodes-1, np);
				
				if (all_but_one.NChooseKInit (aux_list, nk_tuple, np, false))
				{
					bool	remaining;
					long	res;
					
					do 
					{
						remaining = all_but_one.NChooseK (aux_list, nk_tuple);
						
						
						for (long par_idx = 0; par_idx < np; par_idx++)
						{
							long par = nk_tuple.lData[par_idx];
							
							if (par >= node_id) par++;	// map from (n-1) parent space to nodespace
							parents << par;
						//	dag.Store (par, node_id, 1);
						}
						 
						
						score = is_discrete.lData[node_id] ? theTask->b->ComputeDiscreteScore (node_id, parents) : 
															theTask->b->ComputeContinuousScore (node_id);
						
						res = family_scores.Store (score, nk_tuple);
						
						parents.Clear();
						/*
						for (long par = 0; par < num_nodes; par++)
						{
							dag.Store (par, node_id, 0);	// reset graph
						}
						 */
						
					} while (remaining);
					
				} else {
					_String	oops ("Failed to initialize _NTupleStorage object in Bgm::CacheNodeScores().\n");
					WarnError(oops);
				}
				
				(*this_list) && (&family_scores);	// append duplicate to storage
			}
		}
	}	
	pthread_exit (NULL);
}
#endif


//___________________________________________________________________________________________
void Bgm::CacheNodeScores (void)
{
	if (scores_cached)
		return;
	
	
	//char buf [255];
	//sprintf (buf, "\nCaching node scores...\n");
	//BufferToConsole (buf);
	
	
	//_SimpleList		* args	= new _SimpleList((unsigned long) 2);

	
#ifdef __NO_DEFINE_MP__
	_Parameter		nthreads;
	checkParameter (_HYBgm_NThreads, nthreads, 1);
	
	sprintf (buf, "Creating %d threads...\n", (long) nthreads);
	BufferToConsole (buf);
	
	BgmThreads	=	new pthread_t [(long) nthreads];
	
	if (nthreads < 1)
	{
		_String oops ("Dude, I can't run CacheNodeScores() with less than one pthread.\n");
		WarnError (oops);
		
		scores_cached = FALSE;
		return;
	}
	else if (nthreads > 1)
	{
		long	step_size	= (long) (num_nodes / nthreads),
				step_remain	= num_nodes % (long) nthreads;
		
		for (long first_node = 0, tc = 0; tc < (long) nthreads; tc++)
		{
			args->lData[0] = first_node;
			args->lData[1] = first_node + step_size;	// next node
			
			
			sprintf (buf, "Thread %d assigned nodes %d to %d-1", tc, args->lData[0], args->lData[1]);
			BufferToConsole (buf);
			NLToConsole();
			
			
			if (tc == (long)nthreads - 1)	// last thread gets the extra nodes
				args->lData[1] += step_remain;
#ifdef __MACHACKMP__
			if (pthread_create (BgmThreads+tc, NULL, CacheNodeScoreThreadHook, (void *) args))
#else
			if (pthread_create (BgmThreads+tc, NULL, CacheNodeScoreThread, (void *) args))
#endif
			{
				FlagError ("Failed to initialize a POSIX thread in Bgm::CacheNodeScores()");
				exit(1);
			}
			
			first_node += step_size;
		}
		
		for (long tc = 0; tc < (long) nthreads; tc++)
		{
			if ( pthread_join (BgmThreads[tc], NULL) )
			{
				FlagError ("Failed to join POSIX threads in CacheNodeScores()");
				exit (1);
			}
		}
		delete BgmThreads;
		
		
	}
	else
#else
		// carry out as member function
		for (long node_id = 0; node_id < num_nodes; node_id++)
		{
			long		maxp		= max_parents.lData[node_id];
			_List	*	this_list	= (_List *) node_scores.lData[node_id];
			
			this_list->Clear();
			
			_SimpleList	parents;
			_Parameter	score = is_discrete.lData[node_id] ? ComputeDiscreteScore (node_id, parents) : ComputeContinuousScore (node_id);
			_Constant	orphan_score (score);
			
			(*this_list) && (&orphan_score);
			
			if (maxp > 0)
			{
				_Matrix		single_parent_scores (num_nodes, 1, false, true);
				for (long par = 0; par < num_nodes; par++)
				{
					if (par == node_id)
						continue;
					
					parents << par;
					single_parent_scores.Store (par, 0, is_discrete.lData[node_id] ? ComputeDiscreteScore (node_id, parents) : 
												ComputeContinuousScore (node_id));
					parents.Clear();
				}
				(*this_list) && (&single_parent_scores);
			}
			
			if (maxp > 1)
			{
				_SimpleList		all_but_one (num_nodes-1, 0, 1),	// 0, 1, 2, ... , n-1
				aux_list,
				nk_tuple;
				
				for (long np = 2; np <= maxp; np++)
				{
					_NTupleStorage	family_scores (num_nodes-1, np);
					
					if (all_but_one.NChooseKInit (aux_list, nk_tuple, np, false))
					{
						bool	remaining;
						long	res;
						do 
						{
							remaining = all_but_one.NChooseK (aux_list, nk_tuple);
							for (long par_idx = 0; par_idx < np; par_idx++)
							{
								long par = nk_tuple.lData[par_idx];
								if (par >= node_id) par++;
								parents << par;
							}
							score = is_discrete.lData[node_id]	?	ComputeDiscreteScore (node_id, parents) : 
																	ComputeContinuousScore (node_id);
							res = family_scores.Store (score, nk_tuple);
							parents.Clear();						
						} 
						while (remaining);
					} else {
						_String	oops ("Failed to initialize _NTupleStorage object in Bgm::CacheNodeScores().\n");
						WarnError(oops);
					}
					(*this_list) && (&family_scores);	// append duplicate to storage
				}
			}
		}
#endif
	
	
	scores_cached = TRUE;
}



//___________________________________________________________________________________________
_Matrix *	Bgm::ExportNodeScores (void)
{
	// determine the required number of rows
	long		nrows = 0;
	
	for (long child = 0; child < num_nodes; child++)
	{
		long	maxp = max_parents.lData[child];
		nrows++;
		if (maxp > 0)
		{
			nrows += num_nodes - 1;
			if (maxp > 1)
			{
				_List	* childs_list	= (_List *) node_scores.lData[child];
				for (long np = 2; np < maxp; np++)
				{
					_NTupleStorage	* nts	= (_NTupleStorage *) childs_list->lData[np];
					nrows += nts->GetSize();
				}
			}
		}
	}
	
	
	_Matrix *	export_mx = new _Matrix (nrows, 4, false, true);		// child ID, #parents, parent combo #, score
	
	
	for (long index = 0, child = 0; child < num_nodes; child++)
	{
		_List	* childs_list	= (_List *) node_scores.lData[child];
		
		for (long np = 0; np < max_parents.lData[child]; np++)
		{
			if (np == 0)
			{
				_Constant	*	orphan_score = (_Constant *) childs_list->lData[0];
				
				export_mx->Store (index, 0, child);
				export_mx->Store (index, 1, np);
				export_mx->Store (index, 2, 0);
				export_mx->Store (index, 3, orphan_score->Value());
				
				index++;
			}
			else if (np == 1)
			{
				_Matrix	*	single_parent_scores = (_Matrix *) childs_list->lData[1];
				
				for (long parent = 0; parent < num_nodes; parent++)
				{
					if (parent == child) continue;
					export_mx->Store (index, 0, child);
					export_mx->Store (index, 1, np);
					export_mx->Store (index, 2, parent);
					export_mx->Store (index, 3, (*single_parent_scores) (parent, 0));
					index++;
				}
			}
			else
			{
				_NTupleStorage *	family_scores = (_NTupleStorage *) childs_list->lData[np];
				
				for (long family = 0; family < family_scores->GetSize(); family++)
				{
					export_mx->Store (index, 0, child);
					export_mx->Store (index, 1, np);
					export_mx->Store (index, 2, family);	// direct index
					export_mx->Store (index, 3, family_scores->DirectIndex(family));
					index++;
				}
			}
		}
	}
	
	
	return export_mx;
}



//___________________________________________________________________________________________
void	Bgm::ImportNodeScores (_Matrix * import_scores)
{
	
}



#ifdef __NEVER_DEFINED__
//___________________________________________________________________________________________
unsigned long Bgm::IndexIntoCache (long node_id, _SimpleList parents)
{
	/* m0							= index of child node				*/
	/* m1, ..., m[max_parents]		= indices of parent nodes			*/
	/* N							= total number of nodes				*/
	/* e.g. Index = (m0 * N^3) + (m1 * N^2) + (m2 * N) + m3				*/
	/*	if m0 == m1 then child node has only two parents;				*/
	/*	if m0 == m1 == m2 then child node has only one parent (m3);		*/
	/*	if m0 == m1 == m2 == m3, child node has no parents.				*/
	

	unsigned long	index = node_id * integerPower(num_nodes, max_parents);
	long			placer = max_parents - 1;
	
	if (parents.lLength < max_parents)
	{
		// absent parents are replaced by node_id
		for (long par = 0; par < max_parents - parents.lLength; par++)
		{
			index += node_id * integerPower(num_nodes, placer);
			placer--;
		}
	}
	
	for (long par = 0; par < parents.lLength; par++)
	{
		index += parents.lData[par] * integerPower (num_nodes, placer);
		placer--;
	}
	
	if (index == HUGE_VAL)
	{
		_String	oops ("Bgm::IndexIntoCache() overflowed.  Try reducing maximum number of parents.\n");
		WarnError(oops);
	}
	return index;
}


//___________________________________________________________________________________________
//	Return child and parent node id's given index; complements IndexIntoCache()

void	Bgm::IndicesFromCache (long index)
{
	long	nn_pow = 1;
	
	for (long node = max_parents; node >= 0; node--)	// starting from last parent (child = node 0)
	{
		node_ids -> lData[node] = (long int) ((index % (nn_pow * num_nodes)) / nn_pow);
		nn_pow *= num_nodes;
	}
	
	node_ids->Flip();		// reverse order such that child node is first
}

#endif




//___________________________________________________________________________________________
//	Allocate memory for storing node scores and edge posteriors

void	Bgm::InitComputeLists (_List * compute_list)
{
	// adaptive storage for floating point numbers
	_GrowingVector *	newstore;
	checkPointer (newstore = new _GrowingVector);
	
	for (long i = 0; i < num_nodes * num_nodes; i++)
		(*compute_list) && newstore;
	
	DeleteObject (newstore);
}



//___________________________________________________________________________________________
//	Free up memory allocated to cached node scores and edge posteriors

void		Bgm::DumpComputeLists (_List * compute_list)
{	
	for (long i = 0; i < compute_list->lLength; i++)
		((_GrowingVector *) compute_list->lData[i]) -> Clear();
	
	compute_list->Clear();
}



//___________________________________________________________________________________________
_Parameter	Bgm::Compute (_SimpleList * node_order, _List * results)
{
	/*	Calculate equation (8) from Friedman and Koller (2003), i.e. joint probability 
		of data by summing all families at i-th node that are consistent with node order, 
		then taking product across all nodes in the network.  Traverse using AVL indexing.
	
		Node order is a vector of length [num_nodes] containing rank of each node, i.e. 
		a position in an ordered sequence.  In addition, compute marginal posteriors of
		each potential edge in the graph, by Proposition 3.2.								*/
	
	
	_Parameter			log_likel	= 0.;
	
	_GrowingVector		*gv1, *gv2;
	
#ifdef DEBUG_COMPUTE
	char buf [256];
	
	sprintf (buf, "\nCompute(): node order {%d", node_order->lData[0]);
	BufferToConsole (buf);
	
	for (long i = 1; i < num_nodes; i++)
	{
		sprintf (buf, ",%d", node_order->lData[i]);
		BufferToConsole (buf);
	}
	sprintf (buf, "}\n");
	BufferToConsole (buf);
#endif
	
	
	// reset _GrowingVector objects stored in _List object
	for (long i = 0; i < num_nodes * num_nodes; i++)
	{
		gv1 = (_GrowingVector *) results->lData[i];
		gv1 -> ZeroUsed();
	}

	
	if (is_dynamic_graph)	// DBN
	{
		// loop through node order
		for (long nodeIndex = 0; nodeIndex < node_order->lLength; nodeIndex++)
		{
			long				child_node		= node_order->lData[nodeIndex],
								maxp			= max_parents.lData[child_node];
			
			_List			*	score_lists		= (_List *) node_scores.lData[child_node];
			
			
#ifdef DEBUG_COMPUTE
			sprintf (buf, "\tNode %d:\n", child_node);
			BufferToConsole (buf);
#endif
			
			
			gv1 = (_GrowingVector *) results->lData[child_node * num_nodes + child_node];	// store denominator (all scores) in diagonal
			gv1->ZeroUsed();
			
			
			if (child_node % 2 == 1)	// previous time slice, no parents
			{
				_Constant		*	orphan_score	= (_Constant *) (score_lists->lData[0]);
				gv1 -> Store (orphan_score->Value());
				
#ifdef DEBUG_COMPUTE
				sprintf (buf, "\t\torphan score: %f\n", orphan_score->Value());
				BufferToConsole (buf);
#endif
			}
			else
			{
				// handle trivial case of one parent
				_Matrix *	single_parent_scores	= (_Matrix *) (score_lists->lData[1]);
				long		self_parent				= node_order->lData[nodeIndex + 1];
				
				
				gv1 -> Store ((*single_parent_scores) (self_parent, 0));
				
				gv2 = (_GrowingVector *) results->lData[child_node * num_nodes + self_parent];
				gv2 -> Store ((*single_parent_scores) (self_parent, 0));
				
#ifdef DEBUG_COMPUTE
				sprintf (buf, "\t\tsingle parent score: %f\n", (*single_parent_scores) (self_parent, 0));
				BufferToConsole (buf);
#endif
				
				if (maxp > 1)
				{
					// always include self from previous time point, plus...
					_SimpleList			parents,
										eligible_parents;
					
					_Parameter			tuple_score;
					
					_NTupleStorage *	family_scores;
					
					
					// skip every other node in ordering to stay in parental time slice
#ifdef DEBUG_COMPUTE
					sprintf (buf, "\t\teligible non-self parents:");
					BufferToConsole (buf);
#endif
					for (long parIndex = nodeIndex + 3; parIndex < node_order->lLength; parIndex = parIndex + 2)
					{
						long	par = node_order->lData[parIndex];
						
						if (banned_edges (par, child_node) == 0)
						{
							eligible_parents << par;
#ifdef DEBUG_COMPUTE
							sprintf (buf, " %d", par);
							BufferToConsole (buf);
#endif
						}
					}
#ifdef DEBUG_COMPUTE
					NLToConsole ();
#endif
					
					// two parents
					family_scores = (_NTupleStorage *) (score_lists->lData[2]);
					
					for (long ep = 0; ep < eligible_parents.lLength; ep++)
					{
						long	this_parent	= eligible_parents.lData[ep];
						
						// prepare 2-tuple, adjusting for k-tuple indexing (skips child node)
						parents.Clear();
						parents << self_parent-1;
						
						if (this_parent > child_node)
						{
							parents << this_parent - 1;
						}
						else
						{
							parents << this_parent;
						}
						
						parents.Sort(TRUE);
						tuple_score	= family_scores -> Retrieve (parents);
						
						gv1 -> Store (tuple_score);
						
						gv2 = (_GrowingVector *) results->lData[child_node * num_nodes + self_parent];
						gv2 -> Store (tuple_score);
						gv2 = (_GrowingVector *) results->lData[child_node * num_nodes + this_parent];
						gv2 -> Store (tuple_score);
#ifdef DEBUG_COMPUTE
						sprintf (buf, "\t\t2-parent score (%d,%d) = %f\n", parents.lData[0], parents.lData[1], tuple_score);
						BufferToConsole (buf);
#endif
					}
					
					
					if (maxp > 2)	// more than two parents, need to use k-tuple objects
					{
						_SimpleList			indices (eligible_parents.lLength, 0, 1);
						
						for (long nparents = 3; nparents <= maxp; nparents++)
						{
							_SimpleList		subset,
											auxil;
							
							bool			not_finished;
							
							
							if (nparents-1 > eligible_parents.lLength)
								break;
							
							if (indices.NChooseKInit (auxil, subset, nparents-1, false))
							{
								family_scores = (_NTupleStorage *) (score_lists->lData[nparents]);
								
								do
								{
									parents.Clear();
									parents << self_parent;
									
									not_finished = indices.NChooseK (auxil, subset);	// cycle through index combinations
									
									
									for (long i = 0; i < nparents-1; i++)	// convert indices to parent IDs (skipping child)
									{
										long	realized = eligible_parents.lData[subset.lData[i]];
										
										if (realized >= child_node) realized--;
										parents << realized;
									}
									parents.Sort(TRUE);
									
									
									tuple_score	= family_scores -> Retrieve (parents);
									
									gv1 -> Store (tuple_score);
									
									gv2 = (_GrowingVector *) results -> lData [child_node*num_nodes + self_parent];
									gv2 -> Store (tuple_score);
									
									for (long i = 0; i < nparents-1; i++)
									{
										gv2 = (_GrowingVector *) results->lData[child_node * num_nodes + eligible_parents.lData[subset.lData[i]]];
										gv2 -> Store (tuple_score);
									}
#ifdef DEBUG_COMPUTE
									sprintf (buf, "\t\t%d-parent score = %f\n", nparents, tuple_score);
									BufferToConsole (buf);
#endif
								}
								while (not_finished);
							}
						}
					}
				}
			}
			
			gv1 -> _Matrix::Store (0, 0, LogSumExpo(gv1));	// replace first entry with sum, i.e. marginal log-likelihood of child node
			log_likel += (*gv1)(0, 0);
		} // end loop over node order
	}
	else	// static Bayesian network
	{
		for (long nodeIndex = 0; nodeIndex < node_order->lLength; nodeIndex++)
		{
			long				child_node		= node_order->lData[nodeIndex],
								maxp			= max_parents.lData[child_node];
			
			_List			*	score_lists		= (_List *) node_scores.lData[child_node];
			_Constant		*	orphan_score	= (_Constant *) (score_lists->lData[0]);
							
			
			gv1 = (_GrowingVector *) results->lData[child_node * num_nodes + child_node];	// store denominator in diagonal
			gv1->ZeroUsed();
			gv1 -> Store (orphan_score->Value());	// handle case of no parents
			
			
			
			if (maxp > 0)
			{
				// all nodes to the right are potential parents, except banned parents!
				_SimpleList		precedes;
				for (long parIndex = nodeIndex + 1; parIndex < node_order->lLength; parIndex++)
				{
					long	par = node_order->lData[parIndex];
					
					if (banned_edges(par, child_node) == 0)
						precedes << par;
				}
				
				
				// handle trivial case of one parent
				_Matrix *	single_parent_scores	= (_Matrix *) (score_lists->lData[1]);
				
				for (long i = 0; i < precedes.lLength; i++)
				{
					long	par = precedes.lData[i];
					
					gv1 -> Store ((*single_parent_scores) (par, 0));
					gv2 = (_GrowingVector *) results->lData[child_node * num_nodes + par];
					gv2 -> Store ((*single_parent_scores) (par, 0));
				}
				
				
				// more than one parent requires k-tuples
				if (maxp > 1)
				{
					_SimpleList			indices (precedes.lLength, 0, 1);	// populates list with 0, 1, 2, ..., M-1
																			// where M is the number of eligible parents in [precedes]
					_NTupleStorage *	family_scores;
					
					for (long nparents = 2; nparents <= maxp; nparents++)
					{
						_SimpleList		subset,
										auxil;
						
						bool			not_finished;
						
						
						if (nparents > precedes.lLength)	// not enough eligible parents to form tuples!
							break;
						
						
						if (indices.NChooseKInit (auxil, subset, nparents, false))
						{
							_Parameter		tuple_score;
							_SimpleList		parents;
							
							family_scores = (_NTupleStorage *) (score_lists->lData[nparents]);
							
							
							do
							{
								parents.Clear();
								not_finished = indices.NChooseK (auxil, subset);	// cycle through index combinations
								
								
								for (long i = 0; i < nparents; i++)	// convert indices to parent IDs (skipping child)
								{
									long	realized = precedes.lData[subset.lData[i]];
									if (realized >= child_node) realized--;
									parents << realized;
								}
								parents.Sort(TRUE);
								
	#ifdef DEBUG_COMPUTE
								sprintf (buf, "\tnk-tuple: ");
								BufferToConsole (buf);
								for (long k = 0; k < subset.lLength; k++)
								{
									sprintf (buf, "%d ", subset.lData[k]);
									BufferToConsole (buf);
								}
								NLToConsole();
								
								sprintf (buf, "\tparents: ");
								BufferToConsole (buf);
								for (long k = 0; k < parents.lLength; k++)
								{
									sprintf (buf, "%d ", parents.lData[k]);
									BufferToConsole (buf);
								}
								NLToConsole();
	#endif		
								
								tuple_score	= family_scores -> Retrieve (parents);
								
								gv1 -> Store (tuple_score);
								
								for (long i = 0; i < nparents; i++)
								{
									gv2 = (_GrowingVector *) results->lData[child_node * num_nodes + precedes.lData[subset.lData[i]]];
									gv2 -> Store (tuple_score);
								}
							} 
							while (not_finished);
						}
					}
				}
				
	#ifdef DEBUG_COMPUTE
				sprintf (buf, "\tnode %d: ", child_node);
				BufferToConsole (buf);
				
				for (long i = 0; i < gv1->GetUsed(); i++)
				{
					sprintf (buf, "%f ", (_Parameter)(*gv1)(i,0));
					BufferToConsole (buf);
				}
				NLToConsole ();
	#endif
				
			}
			
			gv1 -> _Matrix::Store (0, 0, LogSumExpo(gv1));	// replace first entry with sum, i.e. marginal log-likelihood of child node
			log_likel += (*gv1)(0, 0);
			
		}
		// end loop over child nodes
	}
#ifdef DEBUG_COMPUTE
	sprintf (buf, "log L = %f\n", log_likel);
	BufferToConsole (buf);
#endif
	
	return log_likel;
}



//___________________________________________________________________________________________
_Parameter Bgm::Compute (void)
{
	//CacheNodeScores();
	
	// return posterior probability of a given network defined by 'dag' matrix
	_Parameter	log_score = 0.;
	
	for (long node_id = 0; node_id < num_nodes; node_id++)
		log_score += is_discrete.lData[node_id] ? ComputeDiscreteScore (node_id) : ComputeContinuousScore (node_id);
	
	return log_score;
}



//___________________________________________________________________________________________
_Matrix *	Bgm::Optimize (void)
{
	#ifdef DEBUG_OPTIMIZE
		char		buf [255];
	#endif
	
	CacheNodeScores();
	ResetDag();
	
	_Parameter	num_restarts,			// HBL settings
				num_randomize;
	
	checkParameter (maxNumRestart, num_restarts, 1.);
	checkParameter (numRandomize, num_randomize, num_nodes);
	
	
	_Matrix		log_scores (num_nodes, 1, false, true),
				best_dag (num_nodes, num_nodes, false, true),
				initial_dag (num_nodes, num_nodes, false, true);
	
	_Parameter	log_score,
				old_score,
				best_score = LARGE_NEGATIVE_NUMBER;
	
	
	// load DAG into initial DAG
	for (long row = 0; row < num_nodes; row++)
		for (long col = 0; col < num_nodes; col++)
			initial_dag.Store (row,col,dag(row,col));
	
	
	// greedy hill-climbing algorithm with random restarts
	for (long restart = 0; restart < num_restarts; restart++)
	{
		for (long row = 0; row < num_nodes; row++)
			for (long col = 0; col < num_nodes; col++)
				dag.Store (row,col,initial_dag(row,col));
		//dag = initial_dag;
		
		RandomizeDag(num_nodes);	// randomly modify as many edges as there are nodes in the network (ARBITRARY)
		log_score = Compute();
		
#ifdef DEBUG_OPTIMIZE
		
		sprintf (buf, "Randomized graph...\n");
		BufferToConsole (buf);
		PrintDag();
		
		sprintf (buf, "with initial score = %f\n", log_score);
		BufferToConsole (buf);
#endif
		
		do {
			old_score = log_score;
			
			for (long one_node = 0; one_node < num_nodes; one_node++)
			{
				for (long two_node = one_node + 1; two_node < num_nodes; two_node++)
				{
					long	option = 0;
					
					if (one_node == two_node) continue;
					
					if (dag(one_node, two_node) == 0)
					{
						if (dag(two_node, one_node) == 0)		// no arc
						{
							log_score = TryEdge (one_node, two_node, 0, log_score);	// add1
							log_score = TryEdge (two_node, one_node, 0, log_score);	// add2, will negate add1
						}
						else	// 2-->1
						{
							log_score = TryEdge (one_node, two_node, 1, log_score);	// reverse
							log_score = TryEdge (one_node, two_node, 2, log_score);	// delete, will negate reverse
						}
					}
					else	// 1-->2 exists
					{
						if (dag(one_node, two_node) == 0)
						{
							log_score = TryEdge (one_node, two_node, 1, log_score);	// reverse to 2-->1
							log_score = TryEdge (two_node, one_node, 2, log_score);	// delete 1-->2
						}
						else
						{
							// bidirectional arc should not occur
							log_score = TryEdge (one_node, two_node, 0, log_score);	// add 2-->1 only
							log_score = TryEdge (two_node, one_node, 0, log_score);	// add 1-->2 only
							log_score = TryEdge (one_node, two_node, 2, log_score);	// delete both
						}
					}
					// exit with best improvement
				}
			}
#ifdef DEBUG_OPTIMIZE
			sprintf (buf, "\nnext log score = %f\n", log_score);
			BufferToConsole (buf);
			PrintDag();
#endif
		}
		while (log_score > old_score);	// until no further gain is possible
		
		if (log_score > best_score)
		{
			best_score = log_score;
			for (long row = 0; row < num_nodes; row++)	// best_dag = dag;
				for (long col = 0; col < num_nodes; col++)
					best_dag.Store (row,col,dag(row,col));
		}
	}
	// end random restarts
	{
	for (long row = 0; row < num_nodes; row++)	// dag = best_dag;
		for (long col = 0; col < num_nodes; col++)
			dag.Store (row,col,best_dag(row,col));
	}
	
	// return NULL;
	_Matrix	* result = new _Matrix(num_nodes * num_nodes, 2, false, true);
	result->Store (0,0,best_score);
	{
	for (long row = 0; row < num_nodes; row++)
		for (long col = 0; col < num_nodes; col++)
			result->Store (row*num_nodes+col,1,dag(row,col));
	}

	return (_Matrix*)(result->makeDynamic());
	
}



//___________________________________________________________________________________________
_Parameter	Bgm::TryEdge (long child, long parent, long operation, _Parameter old_score)
{
	_Parameter	log_score,
				last_state1	= dag (child, parent),
				last_state2	= dag (parent, child);
	
#ifdef _NEVER_DEFINED_

	/*
	// PrintDag();
	if (is_marginal)	// modification affects only child node, don't recalculate other scores
		log_marginal = is_discrete(child,0) ? ComputeDiscreteScore (child) : ComputeContinuousScore (child);
	 */
	
	switch (operation)
	{
		case 0:
			if (MarginalSum (FALSE, child) < max_parents && banned_edges(parent,child) == 0)
			{
				dag.Store (parent, child, 1.);	// add an edge [parent] --> [child]
				dag.Store (child, parent, 0.);
			} 
			else
				return old_score;
			break;
		
		case 1:	// reverse an edge
			if (MarginalSum (FALSE, parent) < max_parents && banned_edges(child,parent) == 0)
			{
				dag.Store (child, parent, last_state2);	// switch states
				dag.Store (parent, child, last_state1);
			} 
			else
				return old_score;
			break;
			
		case 2:	// delete an edge
			dag.Store (child, parent, 0.);
			dag.Store (parent, child, 0.);
			break;
		
		default:
			return old_score;
			break;
	}
	
	// PrintDag();
	
	if (IsCyclic())	// new graph is cyclic, reject
	{
		dag.Store (child, parent, last_state1);
		dag.Store (parent, child, last_state2);
		return old_score;
	}
	

	/*
	if (is_marginal)		// update one term only
		log_score = old_score - log_marginal + (is_discrete(child,0) ? ComputeDiscreteScore (child) : ComputeContinuousScore (child));
	else 
		log_score = Compute();
	 */
	
	
	log_score = Compute();
	
	if (log_score > old_score)	// keep new graph
		return log_score;
	else {
		// restore previous settings
		dag.Store (child, parent, last_state1);
		dag.Store (parent, child, last_state2);
	}
	
	return old_score;
#endif
	return 0;
}




//___________________________________________________________________________________________
//	Calculating equation (8) of Friedman and Koller (2003) requires the term:
//		log ( Sum ( score(x_i, U | D) ) ) 
//	where x_i are very small numbers stored as log(x_i)'s but taking the exponential of 
//	the log values can result in numerical underflow.

_Parameter		Bgm::LogSumExpo (_GrowingVector * log_values)
{
	long		size			= log_values->GetUsed();
	
	_Parameter	sum_exponents			= 0.,
				scaling_factor			= 0.,
				largest_negative_log	= (*log_values) (0, 0) + 1.;
	
	char		buf [255];
	
	
	// handle some trivial cases
	if (size == 0)
		return 0.;	// log(exp(0)) = log(1) = 0
	else if (size == 1)
		return (*log_values)(0,0);	// log(exp(log(x)))
	
	
	// find the largest (least negative) log-value
	_Parameter		max_log	 = (*log_values) (0, 0),
					this_log;
	
	for (long val = 1; val < size; val++)
	{
		this_log = (*log_values) (val, 0);
		
		if (this_log > max_log)
		{
			max_log = this_log;
		}
	}
	
	
	// go back through log values and increment by max log value
	//		This will cause some underflow for the smallest values, but we can
	//	use this approximation to handle very large ranges of values.
	//	NOTE: subtracting a negative value.
	for (long val = 0; val < size; val++)
	{
		sum_exponents += exp( (*log_values) (val, 0) - max_log );
	}
	
	return (log(sum_exponents) + max_log);
	
	
#ifdef __NEVER_DEFINED__
	/*
	 // this code is slightly faster, but crashes the HyPhy GUI for some reason... :-/
	 
	_SimpleList		copy_of_logs (size);
	
	for (long val = 0; val < size; val++)
		copy_of_logs.lData[val] = (*log_values) (val,0);
	
	copy_of_logs.Sort(TRUE);	// in ascending order
	scaling_factor = copy_of_logs.lData[size/2];	// take the middle value to adjust others
	
	for (long val = 0; val < size; val++)
		sum_exponents += exp(copy_of_logs.lData[val] - scaling_factor);
	 
	return (log(sum_exponents) + scaling_factor);
	*/
	
	// determine required extent of scaling
	for (long val = 0; val < size; val++)
	{
		_Parameter	this_log		= (*log_values) (val, 0);
		long		num_attempts	= 0;
		
		if (this_log < largest_negative_log)
		{
			largest_negative_log = this_log;
			while (exp(this_log + scaling_factor) == 0. && num_attempts < MAX_LSE_SCALING_ATTEMPTS)		// underflowed
			{
				scaling_factor	+= LOG_SCALING_FACTOR;	// to apply to sum of exponents
				num_attempts++;
			}
		}
		
		if (exp(this_log + scaling_factor) == 0. && num_attempts >= MAX_LSE_SCALING_ATTEMPTS)
			ReportWarning (_String("Bailing out of LogSumExpo() after ") & num_attempts & " attempts to adjust " & this_log & ", last scaling factor = " & scaling_factor);
	}
	
	
	// apply scaling factor to computing the sum of exponents
	for (long val = 0; val < log_values->GetUsed(); val++)
	{
		sum_exponents += exp((*log_values) (val, 0) + scaling_factor);
		/*
		if (sum_exponents > 2147483647)
		{
			_String oops ("Uh-oh. Sum of exponentiated log values in LogSumExpo() is about to overflow.  Bad scaling factor :-/\n");
			WarnError (oops);
			printf ("val = %d, sum_exponents = %lf, scaling_factor = %lf\n", val, sum_exponents, scaling_factor);
		}
		 */
	}
	
	if (log(sum_exponents) - scaling_factor > 99999)
	{
		_String oops ("Nonsense value in LogSumExpo (");
		
		WarnError (oops & (log(sum_exponents) - scaling_factor) & "), bailing out...\n");
		
		for (long foo = 0; foo < size; foo++)
		{
			printf ("%d: %lf + %lf\n", foo, (*log_values) (foo,0), scaling_factor);
		}
	}
	
	
	
	return (log(sum_exponents) - scaling_factor);
#endif
}


//___________________________________________________________________________________________
//	Markov chain Monte Carlo for Bgm
_PMathObj Bgm::CovarianceMatrix (_SimpleList * unused)
{
	_Parameter		mcmc_steps,
					mcmc_burnin,
					mcmc_samples,
					mcmc_nchains,
					mcmc_temp;
	
	_SimpleList *	node_order	= new _SimpleList (num_nodes, 0, 1);
	
	_Matrix *		mcmc_output;
	
	
	// acquisition of HBL arguments with a few sanity checks
	checkParameter (mcmcSteps, mcmc_steps, 0);	
	checkParameter (mcmcBurnin, mcmc_burnin, 0);
	checkParameter (mcmcSamples, mcmc_samples, 0);	
	checkParameter (mcmcNumChains, mcmc_nchains, 1);
	checkParameter (mcmcTemperature, mcmc_temp, 1);
	
	if (mcmc_steps == 0)
	{
		_String oops ("You asked HyPhy to run MCMC with zero steps in the chain! Did you forget to set Bgm_MCMC_STEPS?\n");
		WarnError (oops);
	}
	
	if (mcmc_samples == 0)
	{
		_String oops ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
		WarnError (oops);
	}
	
	
	
	// allocate space for output
	
	
	
	if (mcmc_nchains == 1)
	{
		if (is_dynamic_graph)
		{
			node_order->Permute(2);
		}
		else
		{
			node_order->Permute(1);		// initialize with randomized node ordering
		}
		
		mcmc_output = RunColdChain (node_order, (long) mcmc_burnin, -1);	// negative argument indicates no sampling
		delete (mcmc_output);		// discard burn-in
		mcmc_output = RunColdChain (node_order, (long) mcmc_steps, (long int) (mcmc_steps / mcmc_samples) );
	}
	else		// run coupled MCMC
	{
		node_order->Permute(1);
		// RunHotChain (node_order, (long) 
		
		/* STILL UNDER DEVELOPMENT! */
	}
	
	DeleteObject (node_order);
	
	return mcmc_output;
}



//___________________________________________________________________________________________
_Matrix *	Bgm::RunColdChain (_SimpleList * current_order, long nsteps, long sample_lag)
{
	/*	Execute Metropolis sampler using a swap of two nodes in an ordered sequence				
		as a proposal function.  The posterior probabilities of edges in the network are 
		stored as a member matrix [edge_posteriors].  Note that the total number of 
		possible orderings (i.e. permutations of a sequence of length N) is factorial(N), 
		and can possibly be computed exactly for N < 8.										*/
	
#ifdef DEBUG_RCC
	char		buf [255];
#endif
	
	long		row, 
				first_node, second_node,
				mcmc_samples = (long int) (nsteps / sample_lag);
	
	_SimpleList	proposed_order;
	
	_Matrix	*	mcmc_output		= new _Matrix (mcmc_samples > (num_nodes*num_nodes) ? mcmc_samples : (num_nodes*num_nodes), 2, false, true);
	_List	*	clist			= new _List ();		// compute list
	
	
	InitComputeLists (clist);		// storage for marginal node- and edge-scores
	
	
	
	_Parameter	lk_ratio,
				prob_current_order	= Compute (current_order, clist),
				prob_proposed_order = 0.;
	
	
	/* SLKP 20070926 
		Add user feedback via the console window status bar
	*/
	VerbosityLevel();
	long howManyTimesUpdated = 0; // how many times has the line been updated; is the same as the # of seconds
	TimerDifferenceFunction(false); // save initial timer; will only update every 1 second
	#ifdef __HEADLESS__
		SetStatusLine (_HYBgm_STATUS_LINE_MCMC & (sample_lag < 0 ? _String(" burnin"):empty));
	#else
		#ifndef __UNIX__
			SetStatusLine 	  (empty,_HYBgm_STATUS_LINE_MCMC & (sample_lag < 0 ? _String(" burnin"):empty),empty,0,HY_SL_TASK|HY_SL_PERCENT);
		#else
			_String			* progressReportFile = NULL;
			_Variable 	    * progressFile = CheckReceptacle (&optimizationStatusFile, empty, false);
			
			if (progressFile->ObjectClass () == STRING)
				progressReportFile = ((_FString*)progressFile->Compute())->theString;
			ConsoleBGMStatus (_HYBgm_STATUS_LINE_MCMC & (sample_lag < 0 ?  _String(" burnin"):empty), -1., progressReportFile);
		#endif
	#endif
	
	/* SLKP */
	
	
	
	for (long step = 0; step < nsteps; step++)
	{
		// copy contents of current order to proposed ordering for permutation
		proposed_order.Clear();
		proposed_order.Duplicate (current_order);
		
		
		// swap random nodes in ordered sequence
		if (is_dynamic_graph)
		{
			first_node	= genrand_int32() % (num_nodes / 2);
			second_node = genrand_int32() % (num_nodes / 2);
			while (first_node == second_node)
				second_node = genrand_int32() % (num_nodes / 2);
			
			// swap tuples --- note that this enforces causal edge X(t-1) --> X(t)
			proposed_order.Swap (first_node*2, second_node*2);
			proposed_order.Swap (first_node*2 + 1, second_node*2 + 1);
		}
		else
		{
			first_node	= genrand_int32() % num_nodes;
			second_node = genrand_int32() % num_nodes;
			
			while (first_node == second_node)	// keep sampling until pick different node
				second_node = genrand_int32() % num_nodes;
			
			proposed_order.Swap (first_node, second_node);
		}
		
		
		
		// calculate probability of proposed order --- edge posterior probs loaded into member variable
		prob_proposed_order = Compute (&proposed_order, clist);
		
		
		
		lk_ratio	= exp(prob_proposed_order - prob_current_order);
		
		if (lk_ratio > 1. || genrand_real2() < lk_ratio)	// then set current to proposed order
		{
			(*current_order).Clear();
			(*current_order) << proposed_order;
			prob_current_order = prob_proposed_order;
		}
		
		
		if (sample_lag >= 0)
		{
			// write Markov chain state to output matrix
			
			if (step % sample_lag == 0)
			{
				_GrowingVector *	gv;
				_Parameter			denom;
				
				row = (long int) (step / sample_lag);
				
				// report log likelhood
				mcmc_output->Store (row, 0, prob_current_order);
				
				
#ifdef DEBUG_RCC
				sprintf (buf, "Contents of clist:\n");
				BufferToConsole (buf);
				
				for (long i = 0; i < num_nodes; i++)
				{
					for (long j = 0; j < num_nodes; j++)
					{
						gv = (_GrowingVector *) clist->lData[i * num_nodes + j];
						
						sprintf (buf, "i=%d, j=%d, n=%d: ", i, j, gv->GetUsed());
						BufferToConsole (buf);
						
						for (long k = 0; k < gv->GetUsed(); k++)
						{
							sprintf (buf, "%f ", (*gv)(k,0));
							BufferToConsole (buf);
						}
						
						NLToConsole ();
					}
				}
#endif
				
				
				
				// compute marginal edge posteriors
				for (long child = 0; child < num_nodes; child++)
				{
					gv		= (_GrowingVector *) clist->lData[child * num_nodes + child];
					denom	= (*gv)(0, 0);	// first entry holds node marginal post pr
					
					
					for (long edge, parent = 0; parent < num_nodes; parent++)
					{
						if (parent == child)
							continue;
						
						edge	= child * num_nodes + parent;
						gv		= (_GrowingVector *) clist->lData[edge];
						
						
						if (gv->GetUsed() > 0)
						{
#ifdef DEBUG_RCC
							sprintf (buf, "edge %d: pp = %f, gv = ", edge, (*mcmc_output)(edge, 1));
							BufferToConsole (buf);
							
							for (long i = 0; i < gv->GetUsed(); i++)
							{
								sprintf (buf, "%f ", (*gv) (i,0));
								BufferToConsole (buf);
							}
							
							sprintf (buf, "sum to %f", LogSumExpo (gv));
							BufferToConsole (buf);
#endif
							
							mcmc_output->Store (edge, 1, (*mcmc_output)(edge, 1) + exp (LogSumExpo(gv) - denom));
	
#ifdef DEBUG_RCC
							sprintf (buf, "= %f\n", (*mcmc_output)(edge, 1));
							BufferToConsole (buf);
#endif
						}
						
						
					}
				}
			}
		}
		
		
		/*SLKP 20070926; include progress report updates */
		if (TimerDifferenceFunction(true)>1.0) // time to update
		{
			howManyTimesUpdated ++;
			_String statusLine = _HYBgm_STATUS_LINE_MCMC & (sample_lag < 0 ? _String(" burnin"):empty) & " " & (step+1) & "/" & nsteps 
									& " steps (" & (1.0+step)/howManyTimesUpdated & "/second)";
			#if  defined __HEADLESS__
				SetStatusLine (statusLine);
			#else
				#if !defined __UNIX__
					SetStatusLine 	  (empty,statusLine,empty,100*step/(nsteps),HY_SL_TASK|HY_SL_PERCENT);
					yieldCPUTime(); // let the GUI handle user actions
					if (terminateExecution) // user wants to cancel the analysis
						break;
				#endif
				#if defined __UNIX__ && ! defined __HEADLESS__
					ConsoleBGMStatus (statusLine, 100.*step/(nsteps), progressReportFile);
				#endif
			#endif

			TimerDifferenceFunction(false); // reset timer for the next second
		}
		/* SLKP */
		
	}
	// end loop over steps
	
	
	if (sample_lag >= 0)
		// convert sums of edge posterior probs. in container to means
		for (long edge = 0; edge < num_nodes * num_nodes; edge++)
			mcmc_output->Store (edge, 1, (*mcmc_output)(edge,1) / mcmc_samples);
	
	// release cached scores (assuming that this is the last analysis to be done!)
	DumpComputeLists (clist);
	
	/*SLKP 20070926; include progress report updates */
	#if !defined __UNIX__ || defined __HEADLESS__
		SetStatusLine 	  (_HYBgm_STATUS_LINE_MCMC_DONE & (sample_lag < 0 ? _String(" burnin"):empty));
	#endif
	#if defined __UNIX__ && ! defined __HEADLESS__
		ConsoleBGMStatus (_HYBgm_STATUS_LINE_MCMC_DONE & (sample_lag < 0 ? _String(" burnin"):empty), -1.0, progressReportFile);
	#endif
	/* SLKP */
	
	return mcmc_output;
}



//___________________________________________________________________________________________
#ifdef __NEVER_DEFINED__
void	Bgm::RunHotChain (_SimpleList * current_order, long nsteps, long sample_lag, _Parameter temper)
{
	/*	Execute Metropolis sampler using a swap of two nodes in an ordered sequence				
	 as a proposal function.  The posterior probabilities of edges in the network are 
	 stored as a member matrix [edge_posteriors].  Note that the total number of 
	 possible orderings (i.e. permutations of a sequence of length N) is factorial(N), 
	 and can possibly be computed exactly for N < 8.										*/
	
	_List	*	clist;
	
	_Parameter	prob_current_order = Compute (current_order, clist),
				prob_proposed_order = 0.;
	
	static long	step_counter = 0;
	
	long		sample_interval = (long int) ((nsteps - mcmc_burnin) / mcmc_samples);
	
	_SimpleList	proposed_order;
	
	bool		accept_step = FALSE;
	
	/* SLKP 20070926 
	 Add user feedback via the console window status bar
	 */
#if !defined __UNIX__ || defined __HEADLESS__
	long howManyTimesUpdated = 0; // how many times has the line been updated; is the same as the # of seconds
	TimerDifferenceFunction(false); // save initial timer; will only update every 1 second
	SetStatusLine 	  (empty,_HYBgm_STATUS_LINE_MCMC,empty,0,HY_SL_TASK|HY_SL_PERCENT);
#endif
	/* SLKP */
	
	
	for (long step = 0; step < mcmc_steps; step++)
	{
		// copy contents of current order to proposed ordering for permutation
		proposed_order.Clear();
		proposed_order.Duplicate (current_order);
		
		
		// swap random nodes in ordered sequence
		long	first_node	= genrand_int32 () % proposed_order.lLength,
		second_node = genrand_int32 () % proposed_order.lLength;
		
		while (first_node == second_node)	// keep sampling until pick different node
			second_node = genrand_int32 () % proposed_order.lLength;
		
		proposed_order.Swap(first_node, second_node);
		
		
		// calculate probability of proposed order --- edge posterior probs loaded into member variable
		prob_proposed_order = Compute (&proposed_order, FALSE);
		
		
		_Parameter	lk_ratio	= exp(prob_proposed_order - prob_current_order);
		
		accept_step = FALSE;
		
		if (lk_ratio > 1.)	accept_step = TRUE;	// accept greater likelihood
		else if (genrand_real2() < lk_ratio) accept_step = TRUE;
		
		
		if (accept_step)	// then set current to proposed order
		{
			(*current_order).Clear();
			(*current_order) << proposed_order;
			prob_current_order = prob_proposed_order;
		}
		
		// write Markov chain state to output matrix
		if (step % sample_interval == 0 && step >= mcmc_burnin)	
		{
			long	row = (long int) ((step-mcmc_burnin) / sample_interval);
			
			mcmc_chain->Store (row, 0, prob_current_order);
			
			for (long edge = 0; edge < num_nodes*num_nodes; edge++)
			{
				long	child	= edge % num_nodes,
				parent	= edge / num_nodes;
				
				mcmc_chain->Store (edge, 1, (*mcmc_chain)(edge, 1) + edge_posteriors(child, parent));
			}
		}
		
		/*SLKP 20070926; include progress report updates */
#if !defined __UNIX__ || defined __HEADLESS__
		if (TimerDifferenceFunction(true)>1.0) // time to update
		{
			howManyTimesUpdated ++;
			_String statusLine = _HYBgm_STATUS_LINE_MCMC & " " & (step+1) & "/" & mcmc_steps 
			& " steps (" & (1.0+step)/howManyTimesUpdated & "/second)";
			SetStatusLine 	  (empty,statusLine,empty,100*step/(mcmc_steps),HY_SL_TASK|HY_SL_PERCENT);
			TimerDifferenceFunction(false); // reset timer for the next second
			yieldCPUTime(); // let the GUI handle user actions
			if (terminateExecution) // user wants to cancel the analysis
				break;
		}
#endif
		/* SLKP */
	}
	// end loop over steps
	
	
	for (long edge = 0; edge < num_nodes * num_nodes; edge++)
	{
		// convert sums of edge posterior probs. in container to means
		mcmc_chain->Store (edge, 1, (*mcmc_chain)(edge,1) / mcmc_samples);
	}
	
	
	// release cached scores (assuming that this is the last analysis to be done!)
	DumpComputeLists ();
	
	/*SLKP 20070926; include progress report updates */
#if !defined __UNIX__ || defined __HEADLESS__
	SetStatusLine 	  (_HYBgm_STATUS_LINE_MCMC_DONE);
#endif
	/* SLKP */
	
}

#endif



//__________________________________________________________________________________
//	note: original version of function found in baseobj.cpp (made polymorphic by afyp)
bool 	ReadDataFromFile (char delimiter, _Matrix& data)
{
	_List		readStrings;
	long		columns		= 0,
				lastRead;
	
	_String	prompt_filename = ReturnFileDialogInput();
	if (!prompt_filename)
		return false;
	
	FILE *f = doFileOpen (prompt_filename.getStr(), "r");
	
	
	if (!f)
	{
		_String errMsg ("Failed to open file ");
		errMsg = errMsg & prompt_filename;
		WarnError (errMsg);
		return false;
	}
	
	
	int  c = fgetc (f);
	while (!feof(f))
	{
		lastRead = readStrings.lLength;
		
		_String  *currentTerm = new _String (16L, true);
		checkPointer (currentTerm);
		
		while ((c!='\n')&&(c!='\r')&&(c!=-1))
		{
			if (c==delimiter)
			{
				currentTerm->Finalize();
				readStrings << currentTerm;
				currentTerm = new _String (16L, true);
				checkPointer (currentTerm);
			}
			else	
				(*currentTerm) << (char)c;
			c = fgetc (f);
		}
		
		currentTerm->Finalize();
		if ((readStrings.lLength>lastRead)||(currentTerm->FirstNonSpaceIndex(0,-1,1)>=0))
			readStrings << currentTerm;
		else
			DeleteObject (currentTerm);
		
		c = fgetc (f);
		while ((c=='\n')||(c=='\r'))
			c = fgetc(f);
		
		if (lastRead==0)
			columns = readStrings.lLength;
	}
	
	fclose (f);
	
	if (columns&&(readStrings.lLength>=2*columns)&&(readStrings.lLength%columns==0))
	{
		data.Clear();
		CreateMatrix(&data,readStrings.lLength/columns,columns,false,true,false);
		
		_Constant h, v;
		
		// treat all lines as containing data
		for (lastRead = 0; lastRead < readStrings.lLength/columns; lastRead++)
		{
			h.SetValue (lastRead);
			for (long k=0; k<columns; k++)
			{
				_String * thisString = (_String*)readStrings (lastRead*columns+k);
				if ((thisString->sLength)&&(thisString->FirstNonSpaceIndex (0,-1)>=0))
				{
					_Formula f (*thisString,nil,false);
					v.SetValue (k);
					if (!f.IsEmpty())
						data.MStore (&h,&v,f);
				}
			}
		}
		
		/*
		 for (lastRead = 0; lastRead < columns; lastRead++)
		 names << readStrings (lastRead);
		 */
	}
	else
	{
		_String errMsg ("HyPhy didn't find a well-defined '");
		errMsg = errMsg & delimiter & "'-separated data table in the input file.";
		if (columns == 0)
			errMsg = errMsg & " Input file was empty...";
		else
			if (readStrings.lLength<2*columns)
				errMsg = errMsg & " Input file contained less than a single valid entry line";
		else
			errMsg = errMsg & " The number of fields read (" &(long)readStrings.lLength& ") was not divisible by the number of columns ("& columns &").";
		
		WarnError (errMsg);
		return false;
	}
	
	return true;
}


//___________________________________________________________________________________________
//___________________________________________________________________________________________
//___________________________________________________________________________________________
//___________________________________________________________________________________________
//___________________________________________________________________________________________



