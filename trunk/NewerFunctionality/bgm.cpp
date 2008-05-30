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

#define		__AFYP_DEVELOPMENT__
#define		__MISSING_DATA__

#define		LOG_SCALING_FACTOR			64.
#define		LARGE_NEGATIVE_NUMBER		-99999.0
#define		MAX_LSE_SCALING_ATTEMPTS	4096

#define		DIRICHLET_FLATTENING_CONST	0.5
#define		MAX_FAIL_RANDOMIZE			1000
#define		RANDOMIZE_PROB_SWAP			0.1

#include "bgm.h"

#ifdef __AFYP_DEVELOPMENT__
	#include "mpi.h"
#endif

/*SLKP 20070926; include progress report updates */
	#if !defined __UNIX__ && !defined __HEADLESS__
			#include "HYConsoleWindow.h"
	#endif
/*SLKP*/

_String		_HYBgm_BAN_PARENT_KEY	("BanParent"),
			_HYBgm_BAN_CHILD_KEY	("BanChild"),
			_HYBgm_ENFORCE_PARENT_KEY	("EnforceParent"),
			_HYBgm_ENFORCE_CHILD_KEY	("EnforceChild"),
			
			_HYBgm_NODE_INDEX_KEY	("NodeID"),
			_HYBgm_PRIOR_SIZE_KEY	("PriorSize"),
			_HYBgm_MAX_PARENT_KEY	("MaxParents"),

			_HYBgm_PRIOR_MEAN_KEY	("PriorMean"),		/* for continuous (Gaussian) nodes */
			_HYBgm_PRIOR_PREC_KEY	("PriorVar"),
			
/*SLKP 20070926; add string constants for progress report updates */
			_HYBgm_STATUS_LINE_MCMC			("Running Bgm MCMC"),
			_HYBgm_STATUS_LINE_MCMC_DONE	("Finished Bgm MCMC"),
			_HYBgm_STATUS_LINE_CACHE		("Caching Bgm scores"),
			_HYBgm_STATUS_LINE_CACHE_DONE	("Done caching Bgm scores"),
/*SLKP*/
			/* maxParentString ("Bgm_MAXIMUM_PARENTS"), */
			maxNumRestart	("BGM_NUM_RESTARTS"),
			numRandomize	("BGM_NUM_RANDOMIZE"),
			useNodeOrder	("BGM_USE_NODEORDER"),
			bgmOptimizationMethod ("BGM_OPTIMIZATION_METHOD"),
			
			mcmcNumChains	("BGM_MCMC_NCHAINS"),		// re-use parameter
			mcmcTemperature ("BGM_MCMC_TEMPERATURE"),
			mcmcSteps		("BGM_MCMC_DURATION"),
			mcmcBurnin		("BGM_MCMC_BURNIN"),
			mcmcSamples		("BGM_MCMC_SAMPLES"),
			
			mcemMaxSteps	("BGM_MCEM_MAXSTEPS"),
			mcemBurnin		("BGM_MCEM_BURNIN"),
			mcemSampleSize	("BGM_MCEM_SAMPLES"),
			
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

Bgm::Bgm(_AssociativeList * dnodes, _AssociativeList * cnodes)
{	
	char		buf [255];
	_String		errorMessage;
	
	long		num_discrete	= dnodes->avl.countitems(),
				num_continuous	= cnodes->avl.countitems(),
				node_index;
	
	_Constant	* node_id,						/* pointers into HBL associative array arguments */
				* mp,
				* size, * mean, * precision;
	
	_AssociativeList	* discreteNode, * continuousNode;
	
	
	
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
	CreateMatrix (&enforced_edges, num_nodes, num_nodes, false, true, false);
	CreateMatrix (&prior_sample_size, num_nodes, 1, false, true, false);
	CreateMatrix (&prior_mean, num_nodes, 1, false, true, false);
	CreateMatrix (&prior_precision, num_nodes, 1, false, true, false);
	
	
	is_discrete.Populate (num_nodes, 0, 0);		// allocate space for _SimpleList objects
	num_levels.Populate (num_nodes, 0, 0);
	max_parents.Populate (num_nodes, 0, 0);
	has_missing.Populate (num_nodes, 0, 0);
	
	
	
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
	
	obsData->CheckIfSparseEnough(TRUE);		// check if we can use a more compact representation of the
											// matrix; before including this command, we suffered a
											// noticeable slow-down in this routine. - AFYP
	
	
	// reset data-dependent member variables
	for (long node = 0; node < num_nodes; node++)
	{
		has_missing.lData[node] = 0;
		num_levels.lData[node] = 0;
	}
	
	
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
	
	
	
	if (obsData->GetVDim() == num_nodes)	// make sure the data matrix is compatible with graph
	{
#ifdef __MISSING_DATA__
		long	nrows = obsData->GetHDim();
		
		for (long node = 0; node < num_nodes; node++)
		{
			if (is_discrete.lData[node])
			{
				num_levels.lData[node] = 1;
				
				for (long obs, row = 0; row < nrows; row++)
				{
					obs = (*obsData)(row, node);
					
					if (has_missing[node] == 0 && obs < 0)	// use negative integer values to annotate missing data
					{
						has_missing.lData[node] = 1;
						continue;	// skip next step to check levels
					}
					
					if (obs+1 > num_levels.lData[node])
					{
						num_levels.lData[node] = num_levels.lData[node] + 1;
					}
				}
			}
			else
			{
				num_levels.lData[node] = 0;
				
				// not implementing missing data for continuous nodes yet
				//  not until I've decided on annotation anyhow :-/  afyp
			}
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
		
		sprintf (buf, "Missing (0=FALSE, 1=TRUE): ");
		BufferToConsole (buf);
		for (long i = 0; i < num_nodes; i++)
		{
			sprintf (buf, "%d ", has_missing.lData[i]);
			BufferToConsole (buf);
		}
		NLToConsole ();
	#endif
	
#else
		for (long node = 0; node < num_nodes; node++)	// for every column in matrix
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
#endif
	}
	else
	{
		_String errorMsg ("Number of variables in data matrix do not match number of nodes in graph.");
		WarnError (errorMsg);		// note, this produces a crash because the batch file proceeds to execute
									// BGM routines without data. - AFYP
	}
	
	
	last_node_order.Clear();		// forget the last step taken in MCMC chain
	best_node_order.Clear();
	
	CacheNodeScores();
	
	// re-allocate memory to lists
	// InitComputeLists ();
}



//___________________________________________________________________________________________
void	Bgm::SetWeightMatrix (_Matrix * weights)
{
	if (obsData && (long) (obsData->GetHDim()) == (long) (weights->GetHDim()))
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
void	Bgm::SetGraphMatrix (_Matrix *graph)
{
	/*
	char	bug [255];
	
	sprintf (bug, "Entered Bgm::SetGraphMatrix()\n");
	BufferToConsole (bug);
	*/
	dag = (_Matrix &) (*graph);	// matrix assignment
}


void	Bgm::SetBanMatrix (_Matrix *banMx)
{
	banned_edges = (_Matrix &) (*banMx);
}

void	Bgm::SetEnforceMatrix (_Matrix *enforceMx)
{
	enforced_edges = (_Matrix &) (*enforceMx);
}

void	Bgm::SetBestOrder (_SimpleList * orderList)
{
	best_node_order.Populate(num_nodes, 0, 0);
	
	for (long i = 0; i < num_nodes; i++)
	{
		best_node_order.lData[i] = orderList->lData[i];
	}
}

//___________________________________________________________________________________________
//	For debugging..
void Bgm::PrintGraph (_Matrix * g)
{
	char	buf [255];
	
	if (g)
	{
		for (long row = 0; row < g->GetHDim(); row++)
		{
			for (long col = 0; col < g->GetVDim(); col++)
			{
				sprintf (buf, "%d ", (long) (*g)(row,col));
				BufferToConsole (buf);
			}
			sprintf (buf, "\n");
			BufferToConsole (buf);
		}
	}
	else
	{
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
	}
	NLToConsole();
}


//___________________________________________________________________________________________
void	Bgm::ResetGraph (_Matrix * g)
{
	if (g)
	{
		for (long row = 0; row < num_nodes; row++)
		{
			for (long col = 0; col < num_nodes; col++)
			{
				g->Store (row, col, enforced_edges (row, col));
			}
		}
	}
	else
	{
		for (long row = 0; row < num_nodes; row++)
		{
			for (long col = 0; col < num_nodes; col++)
			{
				// enforced edges must always appear in DAG
				// if none specified, all entries are zeroed
				dag.Store (row, col, enforced_edges (row, col));
			}
		}
	}
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
void	Bgm::RandomizeGraph (_Matrix * graphMx, _SimpleList * order, long num_steps, bool fixed_order)
{
#ifdef __DEBUG_RG__
	char	bug [255];
#endif
	long	step = 0, fail = 0, 
			child, parent, 
			child_idx, parent_idx;
	
	
	// convert order into matrix format of edge permissions
	_Matrix orderMx (num_nodes, num_nodes, false, true);
	
	for (long p_index = 0; p_index < num_nodes; p_index++)
	{
		for (long par = order->lData[p_index], c_index = 0; c_index < num_nodes; c_index++)
		{
			orderMx.Store (par, order->lData[c_index], (p_index > c_index) ? 1 : 0);
		}
	}
	
	
	
	// before we do anything, make sure that graph complies with order /* DEBUGGING */
	for (long parent = 0; parent < num_nodes; parent++)
	{
		for (long child = 0; child < num_nodes; child++)
		{
			if ( (*graphMx)(parent, child)==1 && orderMx(parent, child)==0 )
			{
				PrintGraph (&orderMx);
				PrintGraph (graphMx);
				
				_String oops ("Order-network discrepancy at start of RandomizeGraph(), exiting..");
				WarnError (oops);
				break;
			}
		}
	}
	
	
	// calculate number of parents for each child for rapid access later
	_SimpleList		num_parents;
	// long			graph_density = 0;
	
	num_parents.Populate (num_nodes, 0, 0);
	
	for (child = 0; child < num_nodes; child++)
	{
		for (parent = 0; parent < num_nodes; parent++)
		{
			if ( (*graphMx)(parent, child) > 0)
			{
				num_parents.lData[child]++;
				// graph_density++;
			}
		}
		
		if (num_parents.lData[child] > max_parents.lData[child])
		{
			_String oops ("Number of parents exceeds maximum BEFORE randomization of graph.");
			WarnError (oops);
			break;
		}
	}
	
	
	
	// randomize graph
	do
	{
		// a fail-safe to avoid infinite loops
		if (fail > MAX_FAIL_RANDOMIZE)
		{
			PrintGraph (graphMx);
			PrintGraph (&orderMx);
			
			_String oops ("Failed to modify the graph in GraphMCMC() after MAX_FAIL_RANDOMIZE attempts.");
			WarnError (oops);
			break;
		}
		
		// pick a random edge
		parent_idx	= (genrand_int32() % (num_nodes-1)) + 1;	// shift right to not include lowest node in order
		child_idx	= genrand_int32() % parent_idx;
		
		child = order->lData[child_idx];
		parent = order->lData[parent_idx];
		
#ifdef __DEBUG_RG__
		sprintf (bug, "Try edge %d-->%d\n", parent, child);
		BufferToConsole (bug);
#endif
		
		if (fixed_order || genrand_real2() > RANDOMIZE_PROB_SWAP)
		{
			if ( (*graphMx)(parent,child) == 0 && !banned_edges(parent,child))		// add an edge
			{
				if (num_parents.lData[child] == max_parents.lData[child])
				{
					// move an edge from an existing parent to target parent
					long	deadbeat_dad	= genrand_int32() % max_parents.lData[child];
					long	social_worker	= 0;
					for (; social_worker < num_nodes; social_worker++)
					{
						if ( (*graphMx)(social_worker,child) == 1)
						{
							if (deadbeat_dad)
							{
								deadbeat_dad--;
							}
							else
							{
								break;
							}
						}
					}
					
					if ( !enforced_edges(social_worker, child) )
					{
						graphMx->Store (social_worker, child, 0.);
						num_parents.lData[child]--;
#ifdef __DEBUG_RG__
						sprintf (bug, "Remove and ");
						BufferToConsole (bug);
#endif
					}
					else
					{
						fail++;
					}
				}
				else
				{
					graphMx->Store (parent, child, 1.);
					num_parents.lData[child]++;
					step++;
#ifdef __DEBUG_RG__
					sprintf (bug, "Add\n");
					BufferToConsole(bug);
#endif
				}
			}
			else if ( (*graphMx)(parent,child) == 1 && !enforced_edges(parent,child))		// delete an edge
			{
				graphMx->Store (parent, child, 0.);
				num_parents.lData[child]--;
				step++;
#ifdef __DEBUG_RG__
				sprintf (bug, "Remove\n");
				BufferToConsole (bug);
#endif
			}
			else
			{
				fail++;
			}
		}
		else	// swap nodes in ordering and flip edge if present
		{
			long	buzz = 0;
			
			for (long bystander, i=child_idx+1; i < parent_idx; i++)
			{
				bystander = order->lData[i];
				if (enforced_edges (parent, bystander) || enforced_edges (bystander, child))
				{
					fail++;
					buzz = 1;
				}
			}
			
			if ( buzz == 0 && 
				 ( (*graphMx)(parent,child) == 0
				|| ( (*graphMx)(parent,child) == 1 && !enforced_edges(parent,child) && num_parents.lData[parent] < max_parents.lData[parent] )
				))
			{
				// flip the target edge
				if ( (*graphMx)(parent,child) == 1)
				{
					graphMx->Store (parent,child, 0);
					graphMx->Store (child,parent,1);
					num_parents.lData[child]--;
					num_parents.lData[parent]++;
				}
				
				// flip the other edges affected by node swap
				for (long bystander, i = child_idx+1; i < parent_idx; i++)
				{
					bystander = order->lData[i];
					if ( (*graphMx)(bystander, child) == 1)
					{
						graphMx->Store (bystander, child, 0);
						num_parents.lData[child]--;
						
						if (num_parents.lData[bystander] < max_parents.lData[bystander])
						{
							graphMx->Store (child, bystander, 1);
							num_parents.lData[bystander]++;
						}
					}
					
					if ( (*graphMx)(parent,bystander) == 1)
					{
						graphMx->Store (parent, bystander, 0);
						num_parents.lData[bystander]--;
						
						if (num_parents.lData[parent] < max_parents.lData[parent])
						{
							graphMx->Store (bystander, parent, 1);
							num_parents.lData[parent]++;
						}
					}
				}
				
				// swap nodes in order
				order->lData[parent_idx] = child;
				order->lData[child_idx] = parent;
				
				// refresh order matrix
				for (long p_index = 0; p_index < num_nodes; p_index++)
				{
					for (long par = order->lData[p_index], c_index = 0; c_index < num_nodes; c_index++)
					{
						orderMx.Store (par, order->lData[c_index], (p_index > c_index) ? 1 : 0);
					}
				}
				
				step++;
				
			}
			else
			{
				fail++;
			}	
		}
	}
	while (step < num_steps);
	
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
			for (long row = 0; row < num_nodes; row++)
			{
				for (long col = 0; col < num_nodes; col++)
				{
					dag.Store (row,col,reset_dag(row,col));
				}
			}
			
			
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
			
			// PrintGraph(nil);
			
		} while ( IsCyclic() );
	}
	// end loop over steps
}





//___________________________________________________________________________________________

inline _Parameter Bgm::LnGamma(_Constant * calculator, _Parameter x)
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

_Parameter	Bgm::ComputeDiscreteScore (long node_id, _Matrix * g)
{
	_SimpleList		parents;
	
	for (long par = 0; par < num_nodes; par++)
	{
		if ((*g)(par, node_id) == 1 && is_discrete.lData[par])
			parents << par;
	}
	
	return ComputeDiscreteScore (node_id, parents);
}


//___________________________________________________________________________________________
_Parameter	Bgm::ComputeDiscreteScore (long node_id, _SimpleList & parents)
{
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
	
	
	
#ifdef __MISSING_DATA__
	//	Is node with missing data in Markov blanket of focal node?
	if (has_missing.lData[node_id])
	{
		return (ImputeDiscreteScore (node_id, parents));
	}
	else
	{
		for (long par = 0; par < parents.lLength; par++)
		{
			if (has_missing.lData[parents.lData[par]])
			{
				return (ImputeDiscreteScore (node_id, parents));
			}
		}
	}
#endif
	
	
	_SimpleList		multipliers ((long)1);
	
	// else need to compute de novo
	_Matrix		n_ijk,
				n_ij;		// [i] indexes child nodes, 
							// [j] indexes combinations of values for parents of i-th node, 
							// [k] indexes values of i-th node
	
	long		num_parent_combos	= 1,					// i.e. 'q'
				r_i					= num_levels.lData[node_id];
	
	_Parameter	n_prior_ijk	= 0,
				n_prior_ij	= 0,
				log_score	= 0;
	
	
	
	// how many combinations of parental states are there?
	for (long par = 0; par < parents.lLength; par++)
	{
		num_parent_combos *= num_levels.lData[parents.lData[par]];
		multipliers << num_parent_combos;
	}
	
	
#ifdef __DEBUG_IDS__
	char			buf [255];
	
	sprintf (buf, "Multipliers: ");
	BufferToConsole (buf);
	for (long i = 0; i < multipliers.lLength; i++)
	{
		sprintf (buf, "%d ", multipliers.lData[i]);
		BufferToConsole (buf);
	}
	NLToConsole();
#endif
	
	
	/* count observations by parent combination, using direct indexing */
	CreateMatrix (&n_ijk, num_parent_combos, r_i, false, true, false);
	CreateMatrix (&n_ij, num_parent_combos, 1, false, true, false);
	
	
#ifdef __MISSING_DATA__
	/*  METHODS FOR COMPLETE DATA  */
	for (long obs = 0; obs < obsData->GetHDim(); obs++)
	{
		long	index			= 0,
				//multiplier	= 1,
				child_state		= (*obsData)(obs, node_id);
		
		for (long par = 0; par < parents.lLength; par++)
		{
			long	this_parent			= parents.lData[par],
			this_parent_state	= (*obsData)(obs, this_parent);
			
			index += this_parent_state * multipliers.lData[par];
		}
		
		n_ijk.Store ((long) index, child_state, n_ijk(index, child_state) + 1);
		n_ij.Store ((long) index, 0, n_ij(index, 0) + 1);
	}
	
	if (prior_sample_size (node_id, 0) == 0)	// K2
	{
		for (long j = 0; j < num_parent_combos; j++)
		{
			log_score += LnGamma(calc_bit, num_levels.lData[node_id]);	// (r-1)!
			log_score -= LnGamma(calc_bit, n_ij(j, 0) + num_levels.lData[node_id]);	// (N+r-1)!
			
			for (long k = 0; k < num_levels.lData[node_id]; k++)
				log_score += LnGamma (calc_bit, n_ijk(j,k) + 1);	// (N_ijk)!
		}
	} 
	else	// BDe
	{
		n_prior_ij = prior_sample_size (node_id, 0) / num_parent_combos;
		n_prior_ijk = n_prior_ij / num_levels.lData[node_id];
		
		for (long j = 0; j < num_parent_combos; j++)
		{
			log_score += LnGamma (calc_bit, n_prior_ij) - LnGamma (calc_bit, n_prior_ij + n_ij(j,0));
			
			for (long k = 0; k < num_levels.lData[node_id]; k++)
				log_score += LnGamma (calc_bit, n_prior_ijk + n_ijk(j,k)) - LnGamma (calc_bit, n_prior_ijk);
		}
	}
	
#else
	for (long obs = 0; obs < obsData->GetHDim(); obs++)
	{
		long	index		= 0,
				multiplier	= 1,
				child_state = (*obsData)(obs, node_id);
		
		for (long par = 0; par < parents.lLength; par++)
		{
			long	this_parent			= parents.lData[par],
					this_parent_state	= (*obsData)(obs, this_parent);
				
					
			/*		// this system doesn't work with unequal levels!  :-P  afyp March 31, 2008
			index += this_parent_state * multiplier;		
			multiplier *= num_levels.lData[this_parent];
			 */
			index += this_parent_state * multipliers.lData[par];
		}
		
		n_ijk.Store ((long) index, child_state, n_ijk(index, child_state) + 1);
		n_ij.Store ((long) index, 0, n_ij(index, 0) + 1);
		
	}
	
	
	
	/* compute scoring metric */
	if (prior_sample_size (node_id, 0) == 0)
	{
		/* assume no prior information, use K2 metric */
		for (long j = 0; j < num_parent_combos; j++)
		{
			log_score += LnGamma(calc_bit, r_i);	// (r-1)!
			log_score -= LnGamma(calc_bit, n_ij(j, 0) + r_i);	// (N+r-1)!
			
			for (long k = 0; k < num_levels.lData[node_id]; k++)
				log_score += LnGamma (calc_bit, n_ijk(j,k) + 1);	// (N_ijk)!
			
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
			
			sprintf (buf, "\tlog(N+r-1)! = log %d! = %lf\n", n_ij(j, 0) + num_levels.lData[node_id] - 1, 
					LnGamma(calc_bit, n_ij(j, 0) + num_levels.lData[node_id]));
			BufferToConsole (buf);
			
			for (long k = 0; k < num_levels.lData[node_id]; k++)
			{
				sprintf (buf, "\tlog (N_ijk)! = log %d! = %lf\n", ((long)n_ijk(j,k)), LnGamma (calc_bit, n_ijk(j,k) + 1));
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
#endif
	
	
#ifdef BGM_DEBUG_CDS
	char	buf[255];
	sprintf (buf, "Log score = %f\n\n", log_score);
	BufferToConsole (buf);
#endif
	
	return (log_score);
}



//___________________________________________________________________________________________
//	Use Monte Carlo over assignments to m_ijk to calculate expectation, which can be used for EM (MCEM)
//		but m_ijk are constrained by partial observations, a hassle for bookkeeping  :-/
//		Alternatively, we can perform MCMC over imputations of missing values; which is faster?
_Parameter	Bgm::ImputeDiscreteScore (long node_id, _SimpleList & parents)
{
	_SimpleList		multipliers ((long)1);
	
	_Matrix		n_ijk,
				n_ij;
	
	long		num_parent_combos	= 1,
				r_i					= num_levels.lData[node_id],
				mcem_interval;
	
	_Parameter	log_score	= 0,
				mcem_max_steps,
				mcem_burnin,
				mcem_sample_size;
	
	
	// set MCMC parameters from batch language definitions
	checkParameter (mcemMaxSteps, mcem_max_steps, 0);
	checkParameter (mcemBurnin, mcem_burnin, 0);
	checkParameter (mcemSampleSize, mcem_sample_size, 0);
	
	if (mcem_max_steps == 0 || mcem_sample_size == 0)
	{
		_String oops ("Did you forget to specify a value for BGM_MCEM_MAXSTEPS or BGM_MCEM_SAMPLES?\n");
		WarnError (oops);
	}
	
	mcem_interval = (long) (mcem_max_steps / mcem_sample_size);
	
	
	// count number of parent state combinations
	for (long par = 0; par < parents.lLength; par++)
	{
		num_parent_combos *= num_levels.lData[parents.lData[par]];
		multipliers << num_parent_combos;
	}
	
	CreateMatrix (&n_ijk, num_parent_combos, r_i, false, true, false);
	CreateMatrix (&n_ij, num_parent_combos, 1, false, true, false);
	
	
	
	_GrowingVector		impute,		// to store a copy of incomplete cases
						is_missing; // boolean values mapping to missing values in [impute]
	
	_GrowingVector	*	post_chain = new _GrowingVector(),
					*	prior_denom = new _GrowingVector();
	
	
	// prepare containers for keeping track of state-frequencies per node
	// for later use in imputation
	
	long				max_num_levels	= num_levels.lData [node_id],
						family_size		= parents.lLength + 1;
	
	_SimpleList			family_nlevels (max_num_levels);
	
	_Matrix				m_ijk,		// let m_ijk denote missing values such that n_ijk + m_ijk sums to N_i for all j, k.
						m_ij,
						observed_freqs;
	
	
	CreateMatrix (&m_ijk, num_parent_combos, r_i, false, true, false);	// note these are populated with 0.0's
	CreateMatrix (&m_ij, num_parent_combos, 1, false, true, false);
	
	
	// evaluate the number of levels for each node in family
	for (long par = 0; par < parents.lLength; par++)
	{
		long	par_nlevels = num_levels.lData [ parents.lData[par] ];
		
		family_nlevels << par_nlevels;
		if (par_nlevels > max_num_levels)
		{
			max_num_levels = par_nlevels;
		}
	}
	
	CreateMatrix (&observed_freqs, family_size, max_num_levels, false, true, false);
	
	
	
#ifdef __DEBUG_IDS__
	char	buf[255];
	sprintf (buf, "family_nlevels: ");
	BufferToConsole (buf);
	for (long bug = 0; bug < family_nlevels.lLength; bug++)
	{
		sprintf (buf, "%d ", family_nlevels.lData[bug]);
		BufferToConsole (buf);
	}
	NLToConsole ();
#endif
	
	
	
	// tally complete cases
	for (long index, obs = 0; obs < obsData->GetHDim(); obs++)
	{
		long	child_state = (*obsData) (obs, node_id);
		
		index = 0;
		
		if (child_state > -1)
		{
			observed_freqs.Store (0, child_state, observed_freqs (0, child_state) + 1);
			
			for (long par = 0; par < parents.lLength; par++)
			{
				long	this_parent			= parents.lData[par],
				this_parent_state	= (*obsData) (obs, this_parent);
				
				if (this_parent_state < 0) 
				{
					index = -1;	// flag observation to skip
					break;
				}
				else
				{
					index += this_parent_state * multipliers.lData[par];
					observed_freqs.Store (par+1, this_parent_state, observed_freqs (par+1, this_parent_state) + 1);
				}
			}
		}
		else
		{
			index = -1;
		}
		
		if (index > -1)	// update n_ijk's with complete case
		{
			n_ijk.Store ((long) index, child_state, n_ijk(index, child_state) + 1);
			n_ij.Store ((long) index, 0, n_ij(index, 0) + 1);
		}
		else
		{
			// store family in impute, always in the order (child, parent0, parent1, ...)
			impute.Store ( (_Parameter) child_state);				
			is_missing.Store ( (_Parameter) ((child_state < 0) ? 1 : 0) );
			
			for (long par = 0; par < parents.lLength; par++)
			{
				long	parent_state	= (*obsData) (obs, parents.lData[par]);
				
				impute.Store ( (_Parameter) parent_state);
				is_missing.Store ( (_Parameter) ((parent_state < 0) ? 1 : 0) );
			}
		}
	}
	
	
#ifdef __DEBUG_IDS__
	sprintf (buf, "complete cases N_ijk: \n");
	BufferToConsole (buf);
	for (long j = 0; j < num_parent_combos; j++)
	{
		for (long k = 0; k < family_nlevels.lData[0]; k++)
		{
			sprintf (buf, "%d ", (long) n_ijk (j,k));
			BufferToConsole (buf);
		}
		NLToConsole ();
	}
	NLToConsole ();
	
	sprintf (buf, "impute: \n");
	BufferToConsole (buf);
	for (long bug = 0; bug < impute.GetUsed(); bug++)
	{
		sprintf (buf, "%d ", (long) impute (0, bug));
		BufferToConsole (buf);
		
		if (bug % family_size == family_size - 1)
		{
			NLToConsole ();
		}
	}
	NLToConsole ();
#endif
	
	
	
	// convert observed states into frequencies
	for (long obs_total, i = 0; i < family_size; i++)
	{
		obs_total = 0;
		for (long obs_state = 0; obs_state < family_nlevels.lData[i]; obs_state++)
		{
			obs_total += observed_freqs (i, obs_state);
		}
		
		for (long obs_state = 0; obs_state < family_nlevels.lData[i]; obs_state++)
		{
			observed_freqs.Store (i, obs_state, (_Parameter) (observed_freqs (i, obs_state) / obs_total));
		}
	}
	
	
#ifdef __DEBUG_IDS__
	sprintf (buf, "observed_freqs: \n");
	BufferToConsole (buf);
	
	for (long fam = 0; fam < family_size; fam++)
	{
		for (long lev = 0; lev < family_nlevels.lData[fam]; lev++)
		{
			sprintf (buf, "%f ", observed_freqs (fam, lev));
			BufferToConsole (buf);
		}
		NLToConsole ();
	}
	NLToConsole ();
#endif
	
	
	
	// Initial imputation with random assignments to missing values, based on observed frequencies
	//	We can always burn-in the chain to converge to more realistic assignments
	//	before taking expectation.
	//	Also use this loop to calculate m_ijk's, i.e. cell counts for cases with missing values
	for (long i = 0; i < impute.GetUsed(); i++)
	{
		if (is_missing(0,i))
		{
			double	urn			= genrand_real2();
			long	member		= i % family_size;
			
			for (long lev = 0; lev < family_nlevels.lData[member]; lev++)
			{
				if (urn < observed_freqs (member, lev))
				{
					impute._Matrix::Store (0, i, (_Parameter)lev);
					break;
				}
				else
				{
					urn -= observed_freqs (member, lev);	// rescale random float
				}
			}
		}
	}
	
	
	
	// compute m_ijk's and m_ij's
	long		num_incomplete_cases	= impute.GetUsed() / family_size;
	
	for (long ic = 0; ic < num_incomplete_cases; ic++)
	{
		long	index			= 0,
		child_state		= impute (0, ic * family_size);
		
		for (long par = 0; par < parents.lLength; par++)
		{
			long	parent_state	= impute (0, ic*family_size + par+1);		
			index += parent_state * multipliers.lData[par];
		}
		
		m_ijk.Store ((long) index, child_state, m_ijk(index, child_state) + 1);
		m_ij.Store ((long) index, 0, m_ij(index, 0) + 1);
	}
	
	
	
#ifdef __DEBUG_IDS__
	sprintf (buf, "impute: \n");
	BufferToConsole (buf);
	for (long bug = 0; bug < impute.GetUsed(); bug++)
	{
		sprintf (buf, "%d ", (long) impute (0, bug));
		BufferToConsole (buf);
		
		if (bug % family_size == family_size - 1)
		{
			NLToConsole ();
		}
	}
	NLToConsole ();
	
	sprintf (buf, "m_ijk: \n");
	BufferToConsole (buf);
	for (long j = 0; j < num_parent_combos; j++)
	{
		for (long k = 0; k < family_nlevels.lData[0]; k++)
		{
			sprintf (buf, "%d ", (long) m_ijk (j,k));
			BufferToConsole (buf);
		}
		
		sprintf (buf, "| %d", (long) m_ij (j,0));
		BufferToConsole (buf);
		
		NLToConsole ();
	}
	NLToConsole ();
#endif
	
	
	// Calculate probability of imputation {m_ijk}
	//	For the time being, use multinomial with Dirichlet prior based on non-missing cases (n_ijk)
	//	for each node, i.e., network parameters for empty graph (independence prior).
	//	Later, we will want to allow the user to submit their own prior distribution (e.g., after an iteration
	//	of structural EM) - afyp (March 28, 2008)
	
	// add Jeffrey's invariance constant (0.5) to every cell n_ijk to avoid zero counts?
	
	_Parameter	log_prior_impute		= 0.0;
	
	for (long j = 0; j < num_parent_combos; j++)
	{
		
		log_prior_impute += LnGamma (calc_bit, n_ij(j,0) + DIRICHLET_FLATTENING_CONST * r_i)
		- LnGamma (calc_bit, n_ij(j,0) + DIRICHLET_FLATTENING_CONST * r_i + m_ij(j,0));
		/*
		log_prior_impute += LnGamma (calc_bit, m_ij(j,0)+1) + LnGamma (calc_bit, n_ij(j,0)+DIRICHLET_FLATTENING_CONST*r_i)
		- LnGamma (calc_bit, m_ij(j,0)+n_ij(j,0)+DIRICHLET_FLATTENING_CONST*r_i);
		 */
		for (long k = 0; k < r_i; k++)
		{
			
			log_prior_impute += LnGamma (calc_bit, n_ijk(j,k) + DIRICHLET_FLATTENING_CONST + m_ijk(j,k))
			- LnGamma (calc_bit, n_ijk(j,k) + DIRICHLET_FLATTENING_CONST);
			/*
			log_prior_impute += LnGamma (calc_bit, m_ijk(j,k) + n_ijk(j,k)+DIRICHLET_FLATTENING_CONST) 
			- LnGamma (calc_bit, m_ijk(j,k)+1) - LnGamma (calc_bit, n_ijk(j,k)+DIRICHLET_FLATTENING_CONST);
			 */
		}
	}
	
	
	
	// Calculate posterior probability for initial imputation		
	_Parameter	last_score	= log_prior_impute,
				last_prior_impute = log_prior_impute;
	
	for (long j = 0; j < num_parent_combos; j++)
	{
		last_score += LnGamma(calc_bit, num_levels.lData[node_id]);	// (r-1)!
		last_score -= LnGamma(calc_bit, n_ij(j, 0) + m_ij(j,0) + num_levels.lData[node_id]);	// (N+r-1)!
		
		for (long k = 0; k < num_levels.lData[node_id]; k++)
		{
			last_score += LnGamma (calc_bit, n_ijk (j,k) + m_ijk (j,k) + 1);	// (N_ijk)!
		}
	}
	
	
	
#ifdef __DEBUG_IDS__
	sprintf (buf, "log_prior_impute = %f\nlast_score = %f\n", log_prior_impute, last_score);
	BufferToConsole (buf);
#endif
	
	  /********
	 /  MCMC  /
	 ********/
	
	// Run MCMC over imputations, use as proposal function a single case reassignment (guaranteed to modify m_ijk's)
	//		Take into account case being reassigned so that we don't have to re-sum m_ijk every time
	//		use burn-in and thinning of the sample when calculating the expectation for family posterior probability
	_Parameter	lk_ratio,
				expect_log_score = 0.;
	
	_Matrix		last_impute ((_Matrix &) impute),		// duplicate matrix constructors
				last_m_ijk (m_ijk),
				last_m_ij (m_ij);
	
	
	for (long step = 0; step < mcem_max_steps + mcem_burnin; step++)
	{
		/*	PROPOSAL  */
		
		// choose a random case from the incomplete set and pick a random variable that is missing
		long	imp,
				imp_case,
				imp_var,
				last_state,
				child_state,
				pick;
		
		
		
		while (1)	// this could be optimized by pre-computing a list of missing value indices
		{
			imp = genrand_int32 () % impute.GetUsed();
			
			if (is_missing(0, imp))
			{
				break;
			}
		}
		
		imp_case = (long) (imp / family_size);		// which incomplete case?
		imp_var	 = imp % family_size;		// which variable?
		last_state = impute (0, imp);
		
		
		
#ifdef __DEBUG_IDS__
		sprintf (buf, "imp = %d\nimp_case = %d\nimp_var = %d\nnum_incomplete_cases=%d\n", imp, imp_case, imp_var, num_incomplete_cases);
		BufferToConsole (buf);
#endif
		
		
		long	old_index = 0,
		new_index;
		
		for (long par = 0; par < parents.lLength; par++)
		{
			long	parent_state	= impute (0, imp_case*family_size + par+1);		
			old_index += parent_state * multipliers.lData[par];
		}
		
		
		
		// randomly assign a new state to the variable
		
		if (family_nlevels.lData[imp_var] < 2)
		{
			_String	oops ("WARNING: Attempting to impute missing values for variable having fewer than two known levels..\n");
			WarnError(oops);
			
			continue;	// drop this step from the chain
		}
		else if (family_nlevels.lData[imp_var] == 2)
		{
			// by convention, levels are zero-indexed; if not, this won't work
			impute._Matrix::Store (0, imp, (_Parameter) (last_state ? 0 : 1));
		}
		else		// randomly assign a state (other than the current one) with uniform probability
		{
			pick = genrand_int32() % (family_nlevels.lData[imp_var]-1);
			
			if (pick >= last_state)
			{
				pick++;		// skip current state
			}
			
			impute._Matrix::Store (0, imp_var, (_Parameter)pick);
		}
		
		
#ifdef __DEBUG_IDS__
		sprintf (buf, "proposed impute: \n");
		BufferToConsole (buf);
		for (long bug = 0; bug < impute.GetUsed(); bug++)
		{
			sprintf (buf, "%d ", (long) impute (0, bug));
			BufferToConsole (buf);
			
			if (bug % family_size == family_size - 1)
			{
				NLToConsole ();
			}
		}
		NLToConsole ();
#endif
		
		// given the case just modified, update m_ijk's.
		
		child_state = impute(0, imp_case*family_size);
		
		if (imp_var > 0)
		{
			// changed a parental variable, update index
			new_index = 0;
			
			for (long par = 0; par < parents.lLength; par++)
			{
				long	parent_state	= impute (0, imp_case*family_size + par+1);		
				new_index += parent_state * multipliers.lData[par];
			}
			
#ifdef __DEBUG_IDS__
			sprintf (buf, "old_index = %d\nnew_index = %d\n", old_index, new_index);
			BufferToConsole (buf);
#endif
			//	It shouldn't be necessary to re-calculate log_prior_impute from scratch every time
			//	but this doesn't work properly yet :-P afyp
			
			/*
			 log_prior_impute -= LnGamma (calc_bit, n_ijk(old_index,child_state) + DIRICHLET_FLATTENING_CONST 
			 + m_ijk(old_index, child_state));
			 log_prior_impute -= LnGamma (calc_bit, n_ijk(new_index,child_state) + DIRICHLET_FLATTENING_CONST 
			 + m_ijk(new_index, child_state));
			 */
			m_ijk.Store (old_index, child_state, m_ijk (old_index, child_state) - 1);
			m_ijk.Store (new_index, child_state, m_ijk (new_index, child_state) + 1);
			/*
			 log_prior_impute += LnGamma (calc_bit, n_ijk(old_index,child_state) + DIRICHLET_FLATTENING_CONST 
			 + m_ijk(old_index, child_state));
			 log_prior_impute += LnGamma (calc_bit, n_ijk(new_index,child_state) + DIRICHLET_FLATTENING_CONST 
			 + m_ijk(new_index, child_state));
			 */
			
			// non-integer flattening constant prevents us from operating directly on factorial terms
			/*
			 log_prior_impute -= log ( n_ijk (old_index, child_state) + m_ijk (old_index, child_state) + DIRICHLET_FLATTENING_CONST - 1 );
			 log_prior_impute += log ( n_ijk (new_index, child_state) + m_ijk (new_index, child_state) + DIRICHLET_FLATTENING_CONST );				
			 */
			
			/*
			 log_prior_impute -= - LnGamma (calc_bit, n_ij(old_index,0) + DIRICHLET_FLATTENING_CONST * r_i + m_ij(old_index,0));
			 log_prior_impute -= - LnGamma (calc_bit, n_ij(new_index,0) + DIRICHLET_FLATTENING_CONST * r_i + m_ij(new_index,0));
			 */
			m_ij.Store (old_index, 0, m_ij (old_index, 0) - 1);
			m_ij.Store (new_index, 0, m_ij (new_index, 0) + 1);
			/*
			 log_prior_impute += - LnGamma (calc_bit, n_ij(old_index,0) + DIRICHLET_FLATTENING_CONST * r_i + m_ij(old_index,0));
			 log_prior_impute += - LnGamma (calc_bit, n_ij(new_index,0) + DIRICHLET_FLATTENING_CONST * r_i + m_ij(new_index,0));
			 */
		}
		else
		{
			/*
			 log_prior_impute -= LnGamma (calc_bit, n_ijk(old_index,child_state) + DIRICHLET_FLATTENING_CONST 
			 + m_ijk(old_index, last_state));
			 log_prior_impute -= LnGamma (calc_bit, n_ijk(old_index,child_state) + DIRICHLET_FLATTENING_CONST 
			 + m_ijk(old_index, last_state));
			 */
			m_ijk.Store (old_index, last_state, m_ijk (old_index, last_state) - 1);
			m_ijk.Store (old_index, child_state, m_ijk (old_index, child_state) + 1);
			/*
			 log_prior_impute += LnGamma (calc_bit, n_ijk(old_index,child_state) + DIRICHLET_FLATTENING_CONST 
			 + m_ijk(old_index, child_state));
			 log_prior_impute += LnGamma (calc_bit, n_ijk(old_index,child_state) + DIRICHLET_FLATTENING_CONST 
			 + m_ijk(old_index, child_state));
			 */
			
			/*
			 log_prior_impute -= log ( n_ijk (old_index, last_state) + m_ijk (old_index, last_state) - 1);
			 log_prior_impute += log ( n_ijk (old_index, child_state) + m_ijk (old_index, child_state) );
			 */
			// m_ij's unchanged, parental combination not modified
		}
		
		
		// re-compute prior
		last_prior_impute		= log_prior_impute;
		log_prior_impute		= 0.0;
		/*
		for (long j = 0; j < num_parent_combos; j++)
		{
			log_prior_impute += LnGamma (calc_bit, m_ij(j,0)+1) + LnGamma (calc_bit, n_ij(j,0)+DIRICHLET_FLATTENING_CONST*r_i)
			- LnGamma (calc_bit, m_ij(j,0)+n_ij(j,0)+DIRICHLET_FLATTENING_CONST*r_i);
			for (long k = 0; k < r_i; k++)
			{
				log_prior_impute += LnGamma (calc_bit, m_ijk(j,k) + n_ijk(j,k)+DIRICHLET_FLATTENING_CONST) 
				- LnGamma (calc_bit, m_ijk(j,k)+1) - LnGamma (calc_bit, n_ijk(j,k)+DIRICHLET_FLATTENING_CONST);
			}
		}
		*/
		for (long j = 0; j < num_parent_combos; j++)
		{
			log_prior_impute += LnGamma (calc_bit, n_ij(j,0) + DIRICHLET_FLATTENING_CONST * r_i)
			- LnGamma (calc_bit, n_ij(j,0) + DIRICHLET_FLATTENING_CONST * r_i + m_ij(j,0));
			
			for (long k = 0; k < r_i; k++)
			{
				log_prior_impute += LnGamma (calc_bit, n_ijk(j,k) + DIRICHLET_FLATTENING_CONST + m_ijk(j,k))
				- LnGamma (calc_bit, n_ijk(j,k) + DIRICHLET_FLATTENING_CONST);
			}
		}
		
		
#ifdef __DEBUG_IDS__
		sprintf (buf, "proposed m_ijk: \n");
		BufferToConsole (buf);
		for (long j = 0; j < num_parent_combos; j++)
		{
			for (long k = 0; k < family_nlevels.lData[0]; k++)
			{
				sprintf (buf, "%d ", (long) m_ijk (j,k));
				BufferToConsole (buf);
			}
			NLToConsole ();
		}
		NLToConsole ();
#endif
		
		
		// compute posterior probability for this step, given N_ijk = n_ijk + m_ijk
		log_score = log_prior_impute;
		
		for (long j = 0; j < num_parent_combos; j++)
		{
			log_score += LnGamma(calc_bit, num_levels.lData[node_id]);	// (r-1)!
			log_score -= LnGamma(calc_bit, n_ij(j, 0) + m_ij(j,0) + num_levels.lData[node_id]);	// (N+r-1)!
			
			for (long k = 0; k < num_levels.lData[node_id]; k++)
			{
				log_score += LnGamma (calc_bit,  n_ijk (j,k) + m_ijk (j,k) + 1);	// (N_ijk)!
			}
		}
		
		/* I have a hunch this can all be stream-lined... - afyp April 4, 2008 */
		
		
		// accept step?		(Metropolis-Hastings)
		lk_ratio	= exp(log_score - last_score);
		
#ifdef __DEBUG_IDS__
		sprintf (buf, "prior = %f, log_score = %f, last_score = %f, lk_ratio = %f\n", log_prior_impute, log_score, last_score, lk_ratio);
		BufferToConsole (buf);
#endif
		if (lk_ratio > 1. || genrand_real2() < lk_ratio)	// accept proposed imputation
		{
			last_impute = (_Matrix &)impute;
			last_score = log_score;
			last_prior_impute	= log_prior_impute;
			last_m_ijk = m_ijk;
			last_m_ij = m_ij;
#ifdef __DEBUG_IDS__			
			sprintf (buf, "accept\n\n");
			BufferToConsole (buf);
#endif
		}
		else		// revert to original imputation
		{
			for (long i = 0; i < impute.GetUsed(); i++)
			{
				impute._Matrix::Store (0, i, last_impute (0, i));
			}
			log_score = last_score;
			log_prior_impute	= last_prior_impute;
			m_ijk = last_m_ijk;
			m_ij = last_m_ij;
#ifdef __DEBUG_IDS__				
			sprintf (buf, "reject\n\n");
			BufferToConsole (buf);
#endif
		}
		
		
		
		// handle sampling of chain
		if (step >= mcem_burnin)
		{
			if ( (step - (long)mcem_burnin) % (long)mcem_interval == 0)
			{
				prior_denom->Store (log_prior_impute);	// append to _GrowingVector objects
				post_chain->Store (log_score);
#ifdef __DEBUG_IDS__
				sprintf (buf, "%f, %f\n", log_prior_impute, log_score);
				BufferToConsole (buf);
#endif
			}
		}
	}
	// exit MCEM for-loop
	
	expect_log_score = LogSumExpo (post_chain) - LogSumExpo (prior_denom);
#ifdef __DEBUG_IDS__		
	sprintf (buf, "expected log score = %f\n\n\n", expect_log_score);
	BufferToConsole (buf);
#endif		
	
	DeleteObject (prior_denom);
	DeleteObject (post_chain);
	
	return expect_log_score;
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
//	Master node (0) farms jobs for each graph node to compute nodes.  
//	Compute nodes receive jobs and 
void Bgm::CacheNodeScores (void)
{
	if (scores_cached)
		return;
	
	/*
	char buf [255];
	sprintf (buf, "\nCaching node scores...\n");
	BufferToConsole (buf);
	*/
	
	
#ifdef __AFYP_DEVELOPMENT__ && __HYPHYMPI__
	
	 char buf [255];
	 sprintf (buf, "\nMPI node score caching\n");
	 BufferToConsole (buf);
	 
	
	// MPI_Init() is called in main()
	int		mpi_size,
			my_rank;
	
	_Matrix		single_parent_scores (num_nodes, 1, false, true);
	
	_SimpleList	parents,
				all_but_one (num_nodes-1, 0, 1),
				aux_list,
				nk_tuple;
	
	_Parameter	score;
	
	char		mpi_message [256];
	
	MPI_Status	status;	// contains source, tag, and error code
	
	MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);	// number of processes
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
	
	
	
	if (my_rank == 0)	// master node handles job farming, collating results into its _List member variable
	{
		_Matrix		mpi_node_status (mpi_size, 2, false, true);		// column 1 stores busy signal, 2 stores node id
		long		senderID;
		
		for (long child = 0; child < num_nodes; child++)
		{
			long		maxp			= max_parents.lData[child],
						ntuple_receipt,
						mpi_node		= 0;
			
			bool		remaining;
			
			_List	*	this_list	= (_List *) node_scores.lData[child];
			
			_Parameter	score		= is_discrete.lData[child] ? ComputeDiscreteScore (child, parents) : ComputeContinuousScore (child);
			_Constant	orphan_score (score);
			
			
			this_list->Clear();
			(*this_list) && (&orphan_score);	// handle orphan score locally
			
			
			if (maxp > 0)	// don't bother to farm out trivial cases
			{
				// look for idle nodes
				do
				{
					if (mpi_node_status(mpi_node+1, 0) == 0)
					{
						ReportMPIError(MPI_Send(&child, 1, MPI_LONG, mpi_node+1, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD), true);
						
						sprintf (buf, "Node 0 sent child %d to node %d\n", child, mpi_node+1);
						BufferToConsole (buf);
						
						mpi_node_status.Store (mpi_node+1, 0, 1);	// set busy signal
						mpi_node_status.Store (mpi_node+1, 1, child);
						
						break;
					}
					mpi_node++;
				}
				while (mpi_node < mpi_size-1);
			}
			
			
			if (mpi_node == mpi_size-1)	// all nodes are busy, wait for one to send results
			{
				ReportMPIError (MPI_Recv ((double *)(single_parent_scores.theData), num_nodes, MPI_DOUBLE, MPI_ANY_SOURCE, HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD, &status), false);
				
				(*this_list) && (&single_parent_scores);
				
				senderID = (long) status.MPI_SOURCE;
				mpi_node_status.Store (senderID, 0, 0);	// reset busy signal
				maxp = max_parents.lData[(long)mpi_node_status(senderID, 1)];
				
				sprintf (buf, "Node 0 received scores for child %d from node %d\n", (long)mpi_node_status(senderID, 1), senderID);
				BufferToConsole (buf);
				
				for (long np = 2; np <= maxp; np++)
				{
					_NTupleStorage	family_scores (num_nodes-1, np);
					
					if (all_but_one.NChooseKInit (aux_list, nk_tuple, np, false))
					{
						_Matrix		scores_to_store (nk_tuple.lLength, 1, false, true);
						long		score_index = 0;
						
						// receive nk-tuple indexed node scores from same compute node
						ReportMPIError (MPI_Recv (scores_to_store.theData, nk_tuple.lLength, MPI_DOUBLE, senderID, HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD, &status), false);
						
						do
						{
							remaining = all_but_one.NChooseK (aux_list, nk_tuple);	// update nk-tuple in aux_list
							ntuple_receipt = family_scores.Store (scores_to_store(score_index, 0), nk_tuple);
							score_index++;
						}
						while (remaining);
					}
					else 
					{
						_String	oops ("Failed to initialize _NTupleStorage object in Bgm::CacheNodeScores().\n");
						WarnError(oops);
					}
					
					(*this_list) && (&family_scores);
				}
			}
		}
		// end loop over child nodes
		
		// shut down compute nodes
		for (long shutdown = -1, mpi_node = 0; mpi_node < mpi_size-1; mpi_node++)
		{
			ReportMPIError(MPI_Send(&shutdown, 1, MPI_LONG, mpi_node+1, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD), true);
			
			sprintf (buf, "Node 0 sending shutdown signal to node %d\n", mpi_node+1);
			BufferToConsole (buf);
		}
	}
	
	else			// compute node
	{
		long		node_id,
					maxp;
		_List		list_of_matrices;
		
		while (1)
		{
			ReportMPIError (MPI_Recv (&node_id, 1, MPI_LONG, 0, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD, &status), false);
			
			sprintf (buf, "Node %d received child %d from node %d", my_rank, node_id, status.MPI_SOURCE);
			BufferToConsole (buf);
			
			if (node_id < 0)
			{
				break;	// received shutdown message (-1)
			}
			
			maxp = max_parents.lData[node_id];
			parents.Populate (1,0,0);
			
			for (long par = 0; par < num_nodes; par++)
			{
				if (par == node_id)		// child cannot be its own parent, except in Kansas
					continue;
				
				parents.lData[0] = par;
				single_parent_scores.Store (par, 0, is_discrete.lData[node_id] ? ComputeDiscreteScore (node_id, parents) : 
											ComputeContinuousScore (node_id));
			}
			
			parents.Clear();
			
			for (long np = 2; np <= maxp; np++)
			{
				parents.Populate (np, 0, 0);
				
				if (all_but_one.NChooseKInit (aux_list, nk_tuple, np, false))
				{
					bool		remaining;
					long		tuple_index = 0;
					_Matrix		tuple_scores (nk_tuple.lLength, 1, false, true);
					
					do
					{
						remaining = all_but_one.NChooseK (aux_list, nk_tuple);
						for (long par_idx = 0; par_idx < np; par_idx++)
						{
							long par = nk_tuple.lData[par_idx];
							if (par >= node_id) par++;
							parents.lData[par_idx] = par;
						}
						score = is_discrete.lData[node_id] ? ComputeDiscreteScore (node_id, parents) : 
									ComputeContinuousScore (node_id);
						tuple_scores.Store (tuple_index, 0, (double)score);
					} 
					while (remaining);
					
					list_of_matrices && (&tuple_scores);		// append duplicate
				}
				else 
				{
					_String	oops ("Failed to initialize _NTupleStorage object in Bgm::CacheNodeScores().\n");
					WarnError(oops);
				}
			}
			
			// send results to master node
			ReportMPIError (MPI_Send (&single_parent_scores.theData, num_nodes, MPI_DOUBLE, 0, HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD), true);
			for (long np = 2; np <= maxp; np++)
			{
				_Matrix	* storedMx = (_Matrix *) list_of_matrices.lData[np-2];
				ReportMPIError (MPI_Send (storedMx->theData, storedMx->GetHDim(), MPI_DOUBLE, 0, HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD), true);
			}
			
			sprintf (buf, "Node %d sent scores for child %d to node 0\n", my_rank, node_id);
			BufferToConsole (buf);
		}
	}
	
	
// end HYPHYMPI
#else
	
	
	
#ifdef __NEVER_DEFINE_MP__		// this totally crashes :-/  - AFYP
	_SimpleList		* args	= new _SimpleList((unsigned long) 2);
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
		
		return;
	}
#else		// carry out as member function
	
#if !defined __UNIX__ || defined __HEADLESS__
	TimerDifferenceFunction(false); // save initial timer; will only update every 1 second
	SetStatusLine 	  (empty,_HYBgm_STATUS_LINE_CACHE, empty, 0, HY_SL_TASK|HY_SL_PERCENT);
	
	_Parameter	seconds_accumulator = .0,
				temp;
#endif
	
	for (long node_id = 0; node_id < num_nodes; node_id++)
	{
		long		maxp		= max_parents.lData[node_id];
		_List	*	this_list	= (_List *) node_scores.lData[node_id];	// retrieve pointer to list of scores for this child node
		
		this_list->Clear();		// reset the list
		
		// prepare some containers
		_SimpleList	parents;
		_Parameter	score = is_discrete.lData[node_id] ? ComputeDiscreteScore (node_id, parents) : ComputeContinuousScore (node_id);
		_Constant	orphan_score (score);
		
		(*this_list) && (&orphan_score);
		
#if !defined __UNIX__ || defined __HEADLESS__
		temp = .0;
#endif
		
		if (maxp > 0)
		{
			_Matrix		single_parent_scores (num_nodes, 1, false, true);
			for (long par = 0; par < num_nodes; par++)
			{
				if (par == node_id)		// child cannot be its own parent, except in Kansas
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
		
#if !defined __UNIX__ || defined __HEADLESS__
		if ((temp=TimerDifferenceFunction(true))>1.0) // time to update
		{
			seconds_accumulator += temp;
			
			_String statusLine = _HYBgm_STATUS_LINE_CACHE & " " & (node_id+1) & "/" & num_nodes 
									& " nodes (" & (1.0+node_id)/seconds_accumulator & "/second)";
			SetStatusLine (empty,statusLine,empty,100*(float)node_id/(num_nodes),HY_SL_TASK|HY_SL_PERCENT);
			TimerDifferenceFunction (false); // reset timer for the next second
			yieldCPUTime (); // let the GUI handle user actions
			
			if (terminateExecution) // user wants to cancel the analysis
				break;
		}
#endif
	}
#endif

#endif
	
#if !defined __UNIX__ || defined __HEADLESS__
	SetStatusLine 	  (_HYBgm_STATUS_LINE_CACHE_DONE);
#endif
	
	scores_cached = TRUE;
}



//___________________________________________________________________________________________
_Matrix *	Bgm::ExportGraph (void)
{
	_Matrix	* export_graph = new _Matrix ((_Matrix &)dag);	// duplicator
	
	return (_Matrix *) (export_graph->makeDynamic());
}



//___________________________________________________________________________________________
//	This doesn't work yet -- afyp :-/
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
	
	
	return (_Matrix *)(export_mx->makeDynamic());
}



//___________________________________________________________________________________________
void	Bgm::ImportNodeScores (_Matrix * import_scores)
{
	
}





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
							
							parents.Populate (nparents, 0, 0);	// allocate memory
							
							family_scores = (_NTupleStorage *) (score_lists->lData[nparents]);
							
							
							do
							{
								//parents.Clear();
								not_finished = indices.NChooseK (auxil, subset);	// cycle through index combinations
								
								
								for (long i = 0; i < nparents; i++)	// convert indices to parent IDs (skipping child)
								{
									long	realized = precedes.lData[subset.lData[i]];
									if (realized >= child_node) realized--;
									parents.lData[i] = realized;
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
_Parameter Bgm::Compute (_Matrix * g)
{
	//CacheNodeScores();
	
	// return posterior probability of a given network defined by 'dag' matrix
	_Parameter	log_score = 0.;
	
	for (long node_id = 0; node_id < num_nodes; node_id++)
		log_score += is_discrete.lData[node_id] ? ComputeDiscreteScore (node_id, g) : 0 /*: ComputeContinuousScore (node_id, g)*/;
	
	return log_score;
}




//___________________________________________________________________________________________
_Matrix *	Bgm::Optimize (void)
{
	#ifdef DEBUG_OPTIMIZE
		char		buf [255];
	#endif
	
	if (!scores_cached) 
	{
		CacheNodeScores();
	}
	
	
	_Parameter	num_restarts,			// HBL settings
				num_randomize;
	
	checkParameter (maxNumRestart, num_restarts, 1.);
	checkParameter (numRandomize, num_randomize, num_nodes);
	//checkParameter (useNodeOrder, use_node_order, 0);
	
	
	
	// char		bug [255];
	_Parameter		optMethod;	/* 0 = K2 fixed order; 1 = K2 shuffle order with restarts */
								/* 2 = MCMC fixed order; 3 = MCMC over networks and orders */
	checkParameter (bgmOptimizationMethod, optMethod, 0.);
	
	if (optMethod > 1)
	{
		return GraphMCMC( (optMethod < 3) ? TRUE : FALSE);		// use MCMC to evaluate posterior probability of network space
	}
	
	
	
	_Parameter	this_score, best_score, next_score;
	
	_Matrix		orderMx (num_nodes, num_nodes, false, true),
				best_dag (num_nodes, num_nodes, false, true);
	
	bool		reshuffle_order = FALSE;
	
	
	// if best node order hasn't been estimated, 
	if (optMethod == 1 && best_node_order.lLength == 0)
	{
		reshuffle_order = TRUE;
		best_node_order.Populate (num_nodes, 0, 1);
		best_node_order.Permute (1);
	}
	
	
	//	Convert node order to binary matrix where edge A->B is permitted if
	//	orderMx[B][A] = 1, i.e. B is to the right of A in node order.
	for (long i = 0; i < num_nodes; i++)
	{
		long	child = best_node_order.lData[i];
		
		for (long j = 0; j < num_nodes; j++)
		{
			long	parent = best_node_order.lData[j];
			
			orderMx.Store (parent, child, (j > i) ? 1 : 0);
		}
	}
	
	
	// greedy hill-climbing algorithm (K2)
	ResetGraph (nil);
	best_score = Compute();		// best over all node orders (if we're reshuffling)
	
	
	for (long iter = 0; iter < num_restarts; iter++)
	{
		next_score = Compute();		// reset to empty graph score
		
		
		for (long child = 0; child < num_nodes; child++)
		{
			long		num_parents = 0,
						improvement_flag = 0,
						next_parent_to_add;
			
			do
			{
				for (long parent = 0; parent < num_nodes; parent++)
				{
					if ( (parent != child)						// must meet all conditions!
						&& (dag(parent, child) == 0) 
						&& (orderMx (parent, child) == 1) 
						&& banned_edges(parent, child) == 0)
					{
						dag.Store (parent, child, 1);
						this_score = Compute();
						
						if (this_score > next_score)
						{
							improvement_flag	= 1;
							next_score			= this_score;
							next_parent_to_add	= parent;
						}
						dag.Store (parent, child, 0);	// revert
					}
				}
				
				if (improvement_flag)	// adding another parent improves network score
				{
					dag.Store (next_parent_to_add, child, 1);
					num_parents++;
					improvement_flag = 0;	// reset for next parent
				}
				else
				{
					break;	// unable to improve further
				}
			}
			while (num_parents < max_parents.lData[child]);
		}
		
		
		if (reshuffle_order)
		{
			this_score = Compute();
			
			if (this_score > best_score)
			{
				best_score = this_score;
				best_dag = (_Matrix &) dag;		// store graph optimized from the current ordering
			}
			
			ResetGraph (nil);
			
			best_node_order.Permute (1);
			
			for (long i = 0; i < num_nodes; i++)
			{
				long	child = best_node_order.lData[i];
				for (long j = 0; j < num_nodes; j++)
				{
					long	parent = best_node_order.lData[j];
					orderMx.Store (parent, child, (j > i) ? 1 : 0);
				}
			}
		}
		else
		{
			break;	// only one iteration when node ordering is known
		}
	}
	
	if (reshuffle_order)
	{
		dag = (_Matrix &) best_dag;
		best_node_order.Clear();	// reset to initial state
	}
	
	
	_Matrix	* result = new _Matrix(num_nodes * num_nodes, 2, false, true);
	result->Store (0, 0, Compute());
	
	for (long row = 0; row < num_nodes; row++)
	{
		for (long col = 0; col < num_nodes; col++)
		{
			result->Store (row*num_nodes+col, 1, dag(row, col));
		}
	}
	
	return (_Matrix*)(result->makeDynamic());
}



//___________________________________________________________________________________________
_Matrix *	Bgm::GraphMCMC (bool fixed_order)
{
#ifdef __DEBUG_GMCMC__	
	char		bug [255];
	sprintf (bug, "Entered GraphMCMC()\n");
	BufferToConsole (bug);
#endif
	
	
	_Matrix		*	proposed_graph	= new _Matrix (num_nodes, num_nodes, false, true),
				*	orderMx			= new _Matrix (num_nodes, num_nodes, false, true);
	
	_Matrix			current_graph (num_nodes, num_nodes, false, true),
					best_graph (num_nodes, num_nodes, false, true);
	
	_Parameter		mcmc_steps, mcmc_burnin, mcmc_samples,
					current_score, proposed_score, best_score,
					num_randomize,
					lk_ratio;
	
	long			sampling_interval;
	
	_SimpleList	*	proposed_order	= new _SimpleList();
	_SimpleList		current_order;
	
	
	
	checkParameter (mcmcSteps, mcmc_steps, 0);			// parse HBL settings
	if (mcmc_steps == 0)
	{
		_String oops ("You asked HyPhy to run MCMC with zero steps in the chain! Did you forget to set Bgm_MCMC_STEPS?\n");
		WarnError (oops);
	}
	
	checkParameter (mcmcBurnin, mcmc_burnin, 0);
	checkParameter (mcmcSamples, mcmc_samples, 0);	
	if (mcmc_samples == 0)
	{
		_String oops ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
		WarnError (oops);
	}
	
	sampling_interval = (long) mcmc_steps / (long) mcmc_samples;
	
	checkParameter (numRandomize, num_randomize, num_nodes*num_nodes);
	
	
	
	_Matrix		* result;
	
	checkPointer (result = new _Matrix ( (mcmc_samples > num_nodes*num_nodes) ? mcmc_samples : num_nodes*num_nodes, 3, false, true));
								//	Contents:	(1) chain trace; 
								//				(2) model-averaged edge probabilities; 
								//				(3) best graph;
	
	
	
	if (best_node_order.lLength == 0)
	{
		if (fixed_order)
		{
			_String	oops ("Cannot run fixed-order graph MCMC without defined order (via order-MCMC). Run CovarianceMatrix(receptacle, BGM) first.");
			WarnError (oops);
			
			result = new _Matrix();		// return an empty matrix
		}
		proposed_order->Populate (num_nodes, 0, 1);
		proposed_order->Permute (1);
	}
	else
	{
		proposed_order->Populate(num_nodes, 0, 0);	// allocate memory
		
		for (long i = 0; i < num_nodes; i++)
		{
			proposed_order->lData[i] = best_node_order.lData[i];
		}
	}
	
	
	// make sure that the initial order complies with enforced edges
	for (long pindex = 0; pindex < num_nodes-1; pindex++)
	{
		for (long parent = proposed_order->lData[pindex], cindex = pindex+1; cindex < num_nodes; cindex++)
		{
			long	child	= proposed_order->lData[cindex];
			
			if (enforced_edges(parent,child)==1 && pindex < cindex)
			{
				if (fixed_order)
				{
					_String oops ("Forced to modify fixed order to comply with enforced edge matrix.\n");
					WarnError (oops);
				}
				
				proposed_order->lData[pindex] = child;		// swap ordering
				proposed_order->lData[cindex] = parent;
			}
		}
	}
	
	
	
	//	Populate order matrix -- and while we're at it, set proposed order to current order
	for (long i = 0; i < num_nodes; i++)
	{
		long	child = proposed_order->lData[i];
		
		for (long j = 0; j < num_nodes; j++)
		{
			//	if A precedes B, then set entry (A,B) == 1, permitting edge A->B
			orderMx->Store (proposed_order->lData[j], child, (j > i) ? 1 : 0);
		}
	}
	
	
	
	// status line
#if !defined __UNIX__ || defined __HEADLESS__
	long	updates = 0;
	TimerDifferenceFunction (false);
	SetStatusLine (empty, _HYBgm_STATUS_LINE_MCMC, empty, 0, HY_SL_TASK|HY_SL_PERCENT);
#endif
	
	
	
	//  Reset the graph to either an empty graph, or enforced edges if defined
	ResetGraph (proposed_graph);
	
#ifdef __DEBUG_GMCMC__
	sprintf (bug, "Reset graph to:\n");
	BufferToConsole (bug);
	PrintGraph (proposed_graph);
#endif
	
	RandomizeGraph (proposed_graph, proposed_order, (long)num_randomize, fixed_order);
	
#ifdef __DEBUG_GMCMC__
	sprintf (bug, "Randomized graph:\n");
	BufferToConsole (bug);
	PrintGraph (proposed_graph);
	
	sprintf (bug, "Initial node order:");
	BufferToConsole (bug);
	for (long rex = 0; rex < num_nodes; rex++)
	{
		sprintf (bug, "%d ", proposed_order->lData[rex]);
		BufferToConsole (bug);
	}
	NLToConsole ();
#endif
	
	current_order = (*proposed_order);
	current_graph = (_Matrix &) (*proposed_graph);
	best_graph = (_Matrix &) (*proposed_graph);
	best_score = proposed_score = current_score = Compute(proposed_graph);
	
	
	
	
	// MAIN LOOP
	for (long step = 0; step < mcmc_steps + mcmc_burnin; step++)
	{			
		RandomizeGraph (proposed_graph, proposed_order, 1, fixed_order);
		
		proposed_score = Compute(proposed_graph);
		
		
#ifdef __DEBUG_GMCMC__	
		sprintf (bug, "Current score = %f\n", current_score);
		BufferToConsole (bug);
		sprintf (bug, "Propose graph:\n");
		BufferToConsole (bug);
		PrintGraph (proposed_graph);
		sprintf (bug, "Proposed score = %f\n", proposed_score);
		BufferToConsole (bug);
#endif
		
		
		lk_ratio = exp(proposed_score - current_score);
		
		if (lk_ratio > 1. || genrand_real2() < lk_ratio)		// Metropolis-Hastings
		{
#ifdef __DEBUG_GMCMC__
			sprintf (bug, "Accept move\n");
			BufferToConsole (bug);
#endif
			current_graph		= (_Matrix &) (*proposed_graph);		// accept proposed graph
			current_score	= proposed_score;
			
			for (long foo = 0; foo < num_nodes; foo++)
			{
				current_order.lData[foo] = proposed_order->lData[foo];
			}
			
			if (current_score > best_score)
			{
				best_score = current_score;
				best_graph = (_Matrix &) current_graph;			// keep track of best graph ever visited by chain
				
#ifdef __DEBUG_GMCMC__
				PrintGraph (&best_graph);
				
				sprintf (bug, "Update best node ordering: ");
				BufferToConsole (bug);
				for (long deb = 0; deb < num_nodes; deb++)
				{
					sprintf (bug, "%d ", current_order.lData[deb]);
					BufferToConsole (bug);
				}
				NLToConsole ();
#endif
			}
		}
		else
		{
#ifdef __DEBUG_GMCMC__
			sprintf (bug, "Reject move\n");
			BufferToConsole (bug);
#endif
			// revert proposal
			for (long row = 0; row < num_nodes; row++)
			{
				proposed_order->lData[row] = current_order.lData[row];
				
				for (long col = 0; col < num_nodes; col++)
				{
					proposed_graph->Store(row, col, current_graph(row, col));
				}
			}
		}
		
		
		
		if (step >= mcmc_burnin)		// handle output
		{
			if ( (step-(long)mcmc_burnin) % sampling_interval == 0)
			{
				long	entry	= (long int) ((step-mcmc_burnin) / sampling_interval);
				
				result->Store (entry, 0, current_score);		// update chain trace
				
				for (long row = 0; row < num_nodes; row++)		// update edge tallies
				{
					for (long offset=row*num_nodes, col = 0; col < num_nodes; col++)
					{
						result->Store (offset+col, 1, (*result)(offset+col,1) + current_graph(row, col));
					}
				}
			}
		}
		
#if !defined __UNIX__ || defined __HEADLESS__
		if (TimerDifferenceFunction(true)>1.0) // time to update
		{
			updates ++;
			_String statusLine = _HYBgm_STATUS_LINE_MCMC & " " & (step+1) & "/" & (mcmc_steps + mcmc_burnin) 
									& " steps (" & (1.0+step)/updates & "/second)";
			SetStatusLine 	  (empty,statusLine,empty,100*step/(mcmc_steps + mcmc_burnin),HY_SL_TASK|HY_SL_PERCENT);
			TimerDifferenceFunction(false); // reset timer for the next second
			yieldCPUTime(); // let the GUI handle user actions
			if (terminateExecution) // user wants to cancel the analysis
				break;
		}
#endif
	}
	
	// convert edge tallies to frequencies, and report best graph
	for (long row = 0; row < num_nodes; row++)
	{
		for (long offset=row*num_nodes, col = 0; col < num_nodes; col++)
		{
			result->Store (offset+col, 1, (*result)(offset+col,1)/mcmc_samples);
			result->Store (offset+col, 2, best_graph(row, col));
		}
	}
	
	
	
	
	delete (proposed_graph);
	delete (orderMx);
	delete (proposed_order);
	
	return (_Matrix *) (result->makeDynamic());
}



//___________________________________________________________________________________________
_Parameter	Bgm::TryEdge (long child, long parent, long operation, _Parameter old_score)
{
#ifdef _NEVER_DEFINED_
	_Parameter	log_score,
				last_state1	= dag (child, parent),
				last_state2	= dag (parent, child);
	
	/*
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
	
	// PrintGraph(nil);
	
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
	_Parameter	sum_exponents	= 0.;
	
	
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
}


//___________________________________________________________________________________________
//	Markov chain Monte Carlo for Bgm
_PMathObj Bgm::CovarianceMatrix (_SimpleList * unused)
{
	//char	bug [255];
	
	
	_Parameter		mcmc_steps,
					mcmc_burnin,
					mcmc_samples,
					mcmc_nchains,
					mcmc_temp;
	
	_SimpleList *	node_order;
	
	if (last_node_order.lLength == 0)
	{
		node_order = new _SimpleList (num_nodes, 0, 1);		// populate with series (0, 1, ..., num_nodes)
		
		if (is_dynamic_graph)
		{
			node_order->Permute(2);
		}
		else
		{
			node_order->Permute(1);		// initialize with randomized node ordering
		}
	}
	else
	{
		node_order = new _SimpleList (last_node_order);
	}
	
	
	
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
	
	_Matrix	*	mcmc_output		= new _Matrix (mcmc_samples > (num_nodes*num_nodes) ? mcmc_samples : (num_nodes*num_nodes), 4, false, true);
											//	Contents:	(1) chain trace; (2) edge marginal posterior probabilities;
											//				(3) best node ordering; (4) last node ordering (for restarting chains).
	
	_List	*	clist			= new _List ();		// compute list
	
	
	InitComputeLists (clist);		// storage for marginal node- and edge-scores
	
	
	
	_Parameter	lk_ratio,
				prob_current_order	= Compute (current_order, clist),
				prob_proposed_order = 0.,
				best_prob = prob_current_order;
	
	
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
	
	
	best_node_order = (_SimpleList &) (*current_order);		// initialize class member variable
	
#ifdef DEBUG_RCC
	sprintf (buf, "Initial node order (log L = %f): ", prob_current_order);
	BufferToConsole (buf);
	for (long i = 0; i < num_nodes; i++)
	{
		sprintf (buf, "%d ", best_node_order.lData[i]);
		BufferToConsole (buf);
	}
	NLToConsole();
#endif
	
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
		
#ifdef DEBUG_RCC
		sprintf (buf, "Proposed node order (log L = %f): ", prob_proposed_order);
		BufferToConsole (buf);
		for (long i = 0; i < num_nodes; i++)
		{
			sprintf (buf, "%d ", proposed_order.lData[i]);
			BufferToConsole (buf);
		}
		NLToConsole();
#endif
		
		
		lk_ratio	= exp(prob_proposed_order - prob_current_order);
		
		if (lk_ratio > 1. || genrand_real2() < lk_ratio)	// then set current to proposed order
		{
			(*current_order).Clear();
			(*current_order) << proposed_order;
			prob_current_order = prob_proposed_order;
			
#ifdef DEBUG_RCC
			sprintf (buf, "Accept\n");
			BufferToConsole (buf);
#endif
			
			if (prob_proposed_order > best_prob)
			{
				best_prob = prob_proposed_order;		// update best node ordering
				best_node_order = proposed_order;
			}
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
				
				
#ifdef _DEBUG_RCC
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
	{
		// convert sums of edge posterior probs. in container to means
		for (long edge = 0; edge < num_nodes * num_nodes; edge++)
		{
			mcmc_output->Store (edge, 1, (*mcmc_output)(edge,1) / mcmc_samples);
		}
	}
	
	
	// export node ordering info
	for (long node = 0; node < num_nodes; node++)
	{
		mcmc_output->Store (node, 2, (_Parameter) (best_node_order.lData[node]));
		mcmc_output->Store (node, 3, (_Parameter) (current_order->lData[node]));
	}
	
	last_node_order = (_SimpleList &) (*current_order);
	
	
	
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



//___________________________________________________________________________________________
#ifdef __AFYP_DEVELOPMENT__
void		Bgm::CacheNetworkParameters (void)
{
	_SimpleList		multipliers,	// for indexing parent combinations (j)
					parents;
	
	long			num_parent_combos,
					r_i;
	
	_Matrix			n_ijk,
					theta_ijk;
	
	network_parameters.Clear();		// reset cache
	
	
	for (long i = 0; i < num_nodes; i++)
	{
		// reset variables to defaults
		multipliers.Clear();
		multipliers << 1;
		
		n_ijk.Clear();
		theta_ijk.Clear();
		
		num_parent_combos	= 1;
		r_i					= num_levels.lData[i];
		
		
		// assemble parents based on current graph
		for (long par = 0; par < num_nodes; par++)
		{
			if (dag(par, i) == 1 && is_discrete.lData[par])
			{
				parents << par;	// list of parents will be in ascending numerical order
				num_parent_combos *= num_levels.lData[par];
				multipliers << num_parent_combos;
			}
		}
		
		
		// re-allocate space in matrix
		CreateMatrix (&n_ijk, num_parent_combos, r_i+1, false, true, false);
		CreateMatrix (&theta_ijk, num_parent_combos, r_i, false, true, false);
		
		
		// tally observations
		for (long obs = 0; obs < obsData->GetHDim(); obs++)
		{
			long	index			= 0,
					child_state		= (*obsData)(obs, i);
			
			// is there a faster algorithm for this? - afyp
			for (long par = 0; par < parents.lLength; par++)
			{
				long	this_parent			= parents.lData[par],
						this_parent_state	= (*obsData)(obs, this_parent);
				
				index += this_parent_state * multipliers.lData[par];
			}
			
			n_ijk.Store ((long) index, child_state, n_ijk(index, child_state) + 1);
			n_ijk.Store ((long) index, r_i, n_ijk(index, r_i) + 1);	// store n_ij sum in last entry
		}
		
		
		// compute expected network conditional probabilities (Cooper and Herskovits 1992 eq.16)
		for (long j = 0; j < num_parent_combos; j++)
		{
			_Parameter	denom = n_ijk(j,r_i+1) + r_i;
			
			for (long k = 0; k < r_i; k++)
			{
				theta_ijk.Store (j, k, (n_ijk(j,k)+1.) / denom);
			}
		}
		
		
		// append parameters to member list
		network_parameters && (&theta_ijk);
	}
}
#endif


//___________________________________________________________________________________________
#ifdef __AFYP_DEVELOPMENT2__
_Matrix *	Bgm::SimulateDiscreteCases (long ncases)
{
	/*	Identify nodes without parents in given graph.
		Draw random state assignments based on learned parameters (orphan score).
		Use graph matrix transpose to index child nodes.
		Visit children and randomly assign states dependent on parent assignment.
		Rinse and repeat.
	*/
	CacheNetworkParameters();
	
	_Matrix		gad		= dag.Transpose();
	_Matrix *	cases	= new _Matrix (num_nodes, ncases, false, true);
	
	// proceed using postorder traversal (check for parents before resolving node)
	
}
#endif



//___________________________________________________________________________________________
#ifdef __AFYP_DEVELOPMENT2__
void		PostOrderInstantiate (long node_id, _SimpleList & results)
{
	_SimpleList		parents;
	
	// assign parents first
	for (long par = 0; par < num_nodes; par++)
	{
		if (par != node_id && dag (par, node_id) == 1)
		{
			PostOrderInstantiate (par, results);
			parents << par;
		}
	}
	
	// assign random state given parents
	for (long pindex = 0; pindex < parents.lLength; pindex++)
	{
		long	par = parents.lData[pindex];
		
	}
}


long		Bgm::PredictNode (long node_id, _Matrix * assignments)
{
	
}
#endif





//___________________________________________________________________________________________
//___________________________________________________________________________________________
//___________________________________________________________________________________________
//___________________________________________________________________________________________
//___________________________________________________________________________________________



