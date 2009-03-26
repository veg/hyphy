/*
 *  bayesgraph.cpp
 *  HyPhyXCode
 *
 *  Created by Art Poon on 2/19/09.
 *  Copyright 2009 UCSD AVRC. All rights reserved.
 *
 */

#if defined __AFYP_REWRITE_BGM__

#include "bayesgraph.h"



_String		_HYBgm_NODE_INDEX	("NodeID"),
			_HYBgm_NODETYPE		("NodeType"),
			_HYBgm_MAX_PARENT	("MaxParents"),
			_HYBgm_PRIOR_SIZE	("PriorSize"),
			_HYBgm_PRIOR_MEAN	("PriorMean"),		/* for continuous (Gaussian) nodes */
			_HYBgm_PRIOR_PRECISION	("PriorPrecision"),

			/*SLKP 20070926; add string constants for progress report updates */
			_HYBgm_STATUS_LINE_MCMC			("Running Bgm MCMC"),
			_HYBgm_STATUS_LINE_MCMC_DONE	("Finished Bgm MCMC"),
			_HYBgm_STATUS_LINE_CACHE		("Caching Bgm scores"),
			_HYBgm_STATUS_LINE_CACHE_DONE	("Done caching Bgm scores"),
			/*SLKP*/

			_HYBgm_METHOD_KEY ("BGM_OPTIMIZATION_METHOD"),
			
			_HYBgm_K2_RESTARTS ("BGM_K2_RESTARTS"),
			_HYBgm_K2_RANDOMIZE ("BGM_K2_RANDOMIZE"),
			
			_HYBgm_MCMC_NCHAINS		("BGM_MCMC_NCHAINS"),
			_HYBgm_MCMC_TEMP		("BGM_MCMC_TEMPERATURE"),
			_HYBgm_MCMC_MAXSTEPS	("BGM_MCMC_MAXSTEPS"),
			_HYBgm_MCMC_BURNIN		("BGM_MCMC_BURNIN"),
			_HYBgm_MCMC_SAMPLES		("BGM_MCMC_SAMPLES"),

			_HYBgm_MCMC_PROBSWAP	("BGM_MCMC_PROBSWAP"),
			_HYBgm_MCMC_MAXFAILS	("BGM_MCMC_MAXFAILS"),

			_HYBgm_IMPUTE_MAXSTEPS	("BGM_IMPUTE_MAXSTEPS"),
			_HYBgm_IMPUTE_BURNIN	("BGM_IMPUTE_BURNIN"),
			_HYBgm_IMPUTE_SAMPLES	("BGM_IMPUTE_SAMPLES");


//__________________________________________________________________________________________________________
#ifdef 		__UNIX__

void		ConsoleBGMStatus (_String, _Parameter, _String * fileName = nil);


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



//__________________________________________________________________________________________________________
long		integerPower (long base, long exponent)
{
	//	Rapid computation of an integer power.
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


//__________________________________________________________________________________________________________
_Parameter	LnGamma(_Parameter theValue)
{
	//	Returns Log(gamma(x))
	if (theValue <= 0)
	{
		_String oops ("ERROR: Requested Bgm::LnGamma(x) for x <= 0.");
		WarnError (oops);
		
		return 0;
	}
	
	static _Parameter lngammaCoeff [6] = {	 76.18009172947146,
											-86.50532032941677,
											 24.01409824083091,
											 -1.231739572450155,
											  0.1208650973866179e-2,
											 -0.5395239384953e-5	};
	
	static _Parameter lookUpTable [20] = {	 0.       ,  0.       ,  0.6931472,  1.7917595,  3.1780538,
											 4.7874917,  6.5792512,  8.5251614, 10.6046029, 12.8018275, 
											15.1044126, 17.5023078, 19.9872145, 22.5521639, 25.1912212,
											27.8992714, 30.6718601, 33.5050735, 36.3954452, 39.3398842};
	
	// use look-up table for small integer values
	if (theValue <= 20 && (theValue - (long)theValue) > 0.)	
		return (lookUpTable [(long) theValue - 1]);
	
	
	// else do it the hard way
	_Parameter	x, y, tmp, ser;
	
	y = x = theValue;
	tmp = x + 5.5;
	tmp -= (x+0.5) * log(tmp);
	ser = 1.000000000190015;
	
	for (long j = 0; j <= 5; j++) ser += lngammaCoeff[j] / ++y;
	
	return (-tmp + log(2.506628274631005*ser/x));
}


//__________________________________________________________________________________________________________
_Parameter		LogSumExpo (_GrowingVector * log_values)
{
	//	Computes the sum of a vector whose values are stored as log-transforms,
	//	such that exponentiating the vector entries would result in numerical underflow.
	
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




//__________________________________________________________________________________________________________
//__________________________________________________________________________________________________________
//__________________________________________________________________________________________________________
_BayesianGraphicalModel::_BayesianGraphicalModel (_AssociativeList * nodes)
{
	_String				errorMessage;
	_AssociativeList *	this_avl;
	long				global_max_parents	= 0;
	_Constant *			avl_val;
	
	
	
	theData = NULL;
	
	
	// the fundamental defining characteristic of _BayesianGraphicalModel objects
	num_nodes			= nodes->avl.countitems(),
	
	
	// allocate space to class matrices
	CreateMatrix (&theStructure, num_nodes, num_nodes, false, true, false);
	
	CreateMatrix (&prior_sample_size, num_nodes, 1, false, true, false);
	CreateMatrix (&prior_mean, num_nodes, 1, false, true, false);
	CreateMatrix (&prior_precision, num_nodes, 1, false, true, false);
	
	CreateMatrix (&constraint_graph, num_nodes, num_nodes, false, true, false);
	
	// allocate space for _SimpleList objects
	data_type.Populate (num_nodes, 0, 0);
	max_parents.Populate (num_nodes, 0, 0);
	num_levels.Populate(num_nodes, 0, 0);
	has_missing.Populate(num_nodes, 0, 0);
	
	
	for (long node = 0; node < num_nodes; node++)
	{
		this_avl = (_AssociativeList *) (nodes->GetByKey (node, ASSOCIATIVE_LIST));
		
		// AVL should include following items:	"NodeID"		- variable name
		//										"NodeType"		- 0 = discrete, 1 = continuous (Gaussian)
		//										"MaxParents"	- maximum number of parents
		//										"PriorSize"		- hyperparameter for discrete node (BDe)
		//										"PriorMean"		- hyperparameter for Gaussian node
		//										"PriorPrecision" - hyperparameter for Gaussian node
		
		// node type (0 = discrete, 1 = continuous)
		if (avl_val = (_Constant *) (this_avl->GetByKey (_HYBgm_NODETYPE, NUMBER)))
		{
			data_type.lData[node] = (long)(avl_val->Value());
		}
		else
		{
			errorMessage = _String ("Missing NodeType in associative array for node ") & node;
			break;
		}
		
		// max parents
		if (avl_val = (_Constant *) (this_avl->GetByKey (_HYBgm_MAX_PARENT, NUMBER)))
		{
			max_parents.lData[node] = (long)(avl_val->Value());
			if (max_parents.lData[node] > global_max_parents)
			{
				global_max_parents = max_parents.lData[node];
			}
		}
		else
		{
			errorMessage = _String ("Missing MaxParents in associative array for node ") & node;
			break;
		}
		
		// prior sample size
		if (avl_val = (_Constant *) (this_avl->GetByKey (_HYBgm_PRIOR_SIZE, NUMBER)))
		{
			prior_sample_size.Store (node, 0, (_Parameter) (avl_val->Value()));
		}
		else
		{
			errorMessage = _String ("Missing PriorSize in associative array for node ") & node;
			break;
		}
		
		// prior mean (Gaussian)
		if (avl_val = (_Constant *) (this_avl->GetByKey (_HYBgm_PRIOR_MEAN, NUMBER)))
		{
			prior_mean.Store (node, 0, (_Parameter) (avl_val->Value()));
		}
		else if (!avl_val && data_type.lData[node] == 1)
		{
			errorMessage = _String ("Missing PriorMean in associative array for node ") & node;
			break;
		}
		
		// prior precision (Gaussian)
		if (avl_val = (_Constant *) (this_avl->GetByKey (_HYBgm_PRIOR_PRECISION, NUMBER)))
		{
			prior_precision.Store (node, 0, (_Parameter) (avl_val->Value()));
		}
		else if (!avl_val && data_type.lData[node] == 1)
		{
			errorMessage = _String ("Missing PriorPrecision in associative array for node ") & node;
			break;
		}
	}
	
	
	if (errorMessage.sLength)
	{
		WarnError (errorMessage);
	}
	
	
	// allocate node score cache
	_List	emptyList (global_max_parents + 1);
	
	for (long node = 0; node < num_nodes; node++)
	{
		node_score_cache && (&emptyList);
	}
	
	scores_cached = FALSE;
	
	ReportWarning (_String ("Constructed BGM with ") & num_nodes & " nodes.");
}



//__________________________________________________________________________________________________________
_BayesianGraphicalModel::~_BayesianGraphicalModel (void)
{
	
}




//__________________________________________________________________________________________________________
bool _BayesianGraphicalModel::SetDataMatrix (_Matrix * data)
{
	ReportWarning (_String ("Entered SetDataMatrix()."));
	
	if (data->GetVDim() == num_nodes)
	{
		if (theData) DeleteObject (theData);	// purge duplicate from memory (makeDynamic)
		
		theData = data;	// reassign pointer
		theData->CheckIfSparseEnough (TRUE);
		
		scores_cached = FALSE;
		
		
		for (long nrows = theData->GetHDim(), node = 0; node < num_nodes; node++)
		{
			// if discrete node, compute number of levels and missingness
			if (data_type.lData[node] == 0)
			{
				num_levels.lData[node] = 1;
				
				for (long val, row = 0; row < nrows; row++)
				{
					val = (*theData) (row, node);
					
					if (val < 0 && has_missing.lData[node] == 0)
					{
						has_missing.lData[node] = 1;
						continue;
					}
					
					if (val + 1 > num_levels.lData[node])
					{
						num_levels.lData[node]++;
					}
				}
			}
			
			// continuous
			else if (data_type.lData[node] == 1)
			{
				for (long val, row = 0; row < nrows; row++)
				{
					if (val == MISSING_CONTINUOUS_VALUE && has_missing.lData[node] == 0)
					{
						has_missing.lData[node] = 1;
						break;
					}
				}
			}
		}
	}
	else
	{
		_String errorMsg ("ERROR: Number of variables in data does not match number of nodes in graph.");
		WarnError (errorMsg);
		return (FALSE);
	}
	
	CacheNodeScores();
	
	ReportWarning (_String ("Set data matrix."));
	return (TRUE);
}



//__________________________________________________________________________________________________________
bool _BayesianGraphicalModel::SetConstraints (_Matrix * constraints)
{
	if (constraints->GetHDim() == num_nodes)
	{
		constraint_graph = (_Matrix &) (*constraints);
		return (TRUE);
	}
	
	_String errorMsg ("ERROR: Constraint matrix incompatible dimensions to graph.");
	WarnError (errorMsg);
	return (FALSE);
}



//__________________________________________________________________________________________________________
bool _BayesianGraphicalModel::SetStructure (_Matrix * structure)
{
	//	Set [theStructure] to matrix argument from HBL SetParameter().
	//	-------------------------------------------------------------
	if (structure->GetHDim() == num_nodes)
	{
		// check graph against constraint matrix
		for (long row = 0; row < num_nodes; row++)
		{
			for (long col = 0; col < num_nodes; col++)
			{
				if (constraint_graph(row,col) < 0 && (*structure)(row,col) == 1)
				{
					_String errorMsg ("ERROR: Structure contains banned edge: ");
					errorMsg = errorMsg & _String (row) & _String ("->") & _String (col);
					WarnError (errorMsg);
					return (FALSE);
				}
				
				if (constraint_graph(row,col) > 0 && (*structure)(row,col) == 0)
				{
					_String errorMsg ("ERROR: Structure lacks enforced edge:");
					errorMsg = errorMsg & _String (row) & _String ("->") & _String (col);
					WarnError (errorMsg);
					return (FALSE);
				}
			}
		}
		
		
		// is new structure compatible with node_order?
		if (!GraphObeysOrder (theStructure, node_order))
		{
			// need to reset node_order
			_SimpleList * templist;
			
			templist = GetOrderFromGraph (theStructure);
			node_order = (_SimpleList &) (*templist);
			DeleteObject (templist);
		}
		
		theStructure = (_Matrix &) (*structure);
		return (TRUE);
	}
	
	_String errorMsg ("ERROR: Structure incompatible dimensions to graph.");
	WarnError (errorMsg);
	return (FALSE);
}


//__________________________________________________________________________________________________________
bool _BayesianGraphicalModel::SetNodeOrder (_SimpleList * order)
{
	//	Set node ordering to vector argument from HBL SetParameter().
	//	-------------------------------------------------------------
	if (order->lLength == num_nodes)
	{
		if (GraphObeysOrder (theStructure, (_SimpleList &) *order))
		{
			node_order.Populate(num_nodes, 0, 0);	// reset member variable
			
			for (long i = 0; i < num_nodes; i++)
			{
				node_order.lData[i] = order->lData[i];
			}
			
			return (TRUE);
		}
		else
		{
			_String errorMsg ("ERROR: Node order incompatible with current graph.");
			WarnError (errorMsg);
			return (FALSE);
		}
	}
	else
	{
		_String errorMsg ("ERROR: Node order argument incorrect length.");
		WarnError (errorMsg);
		return (FALSE);
	}
}



//__________________________________________________________________________________________________________
_Parameter _BayesianGraphicalModel::Compute (void)
{
	//	Return posterior probability of current network structure.
	//	-------------------------------------------------------------
	
	_Parameter	log_score = 0.;
	
	for (long node_id = 0; node_id < num_nodes; node_id++)
	{
		log_score += data_type.lData[node_id] ? ComputeDiscreteScore (node_id) : ComputeContinuousScore (node_id);
	}
	
	return log_score;
}



//__________________________________________________________________________________________________________
_Parameter _BayesianGraphicalModel::Compute (_Matrix & g)
{
	//	Return posterior probability of given network structure.
	//	-------------------------------------------------------------
	
	_Parameter	log_score = 0.;
	
	for (long node_id = 0; node_id < num_nodes; node_id++)
		log_score += data_type.lData[node_id] ? ComputeDiscreteScore (node_id, g) : ComputeContinuousScore (node_id, g);
	
	return log_score;
}



//__________________________________________________________________________________________________________
_Parameter	_BayesianGraphicalModel::Compute (_SimpleList & node_order, _List * marginals)
{
	//	Return posterior probability of given node order by integrating over all compatible
	//		structures; also return marginal posterior probabilities for edges.
	//	-------------------------------------------------------------
	
	_Parameter			log_likel	= 0.;
	_GrowingVector		*gv1, *gv2;
	
	
	// reset _GrowingVector objects stored in _List object
	for (long i = 0; i < num_nodes * num_nodes; i++)
	{
		gv1 = (_GrowingVector *) marginals->lData[i];
		gv1 -> ZeroUsed();
	}
	
	
	for (long nodeIndex = 0; nodeIndex < node_order.lLength; nodeIndex++)
	{
		long				child_node		= node_order.lData[nodeIndex],
							maxp			= max_parents.lData[child_node];
		
		_List			*	score_lists		= (_List *) node_score_cache.lData[child_node];
		_Constant		*	orphan_score	= (_Constant *) (score_lists->lData[0]);
		
		
		gv1 = (_GrowingVector *) marginals->lData[child_node * num_nodes + child_node];	// store denominator in diagonal
		gv1->ZeroUsed();
		gv1 -> Store (orphan_score->Value());	// handle case of no parents
		
		
		
		if (maxp > 0)
		{
			// all nodes to the right are potential parents, except banned parents!
			_SimpleList		precedes;
			for (long parIndex = nodeIndex + 1; parIndex < node_order.lLength; parIndex++)
			{
				long	par = node_order.lData[parIndex];
				
				if (constraint_graph(par, child_node) >= 0)	// not banned
					precedes << par;
			}
			
			
			// handle trivial case of one parent
			_Matrix *	single_parent_scores	= (_Matrix *) (score_lists->lData[1]);
			
			for (long i = 0; i < precedes.lLength; i++)
			{
				long	par = precedes.lData[i];
				
				gv1 -> Store ((*single_parent_scores) (par, 0));
				gv2 = (_GrowingVector *) marginals->lData[child_node * num_nodes + par];
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
							
							tuple_score	= family_scores -> Retrieve (parents);
							
							gv1 -> Store (tuple_score);
							
							for (long i = 0; i < nparents; i++)
							{
								gv2 = (_GrowingVector *) marginals->lData[child_node * num_nodes + precedes.lData[subset.lData[i]]];
								gv2 -> Store (tuple_score);
							}
						} 
						while (not_finished);
					}
				}
			}
		}
		
		gv1 -> _Matrix::Store (0, 0, LogSumExpo(gv1));	// replace first entry with sum, i.e. marginal log-likelihood of child node
		log_likel += (*gv1)(0, 0);
		
	}
	// end loop over child nodes
	
	return log_likel;	
}



//__________________________________________________________________________________________________________
_Parameter	_BayesianGraphicalModel::ComputeDiscreteScore (long node_id)
{
	_SimpleList		parents;
	
	for (long par = 0; par < num_nodes; par++)
	{
		if (theStructure(par, node_id) == 1 && data_type.lData[par] == 1)
			parents << par;
	}
	
	return ComputeDiscreteScore (node_id, parents);
}


//__________________________________________________________________________________________________________
_Parameter	_BayesianGraphicalModel::ComputeDiscreteScore (long node_id, _Matrix & g)
{
	_SimpleList		parents;
	
	for (long par = 0; par < num_nodes; par++)
	{
		if ( g(par, node_id) == 1 && data_type.lData[par] == 1)
			parents << par;
	}
	
	return ComputeDiscreteScore (node_id, parents);
}


//__________________________________________________________________________________________________________
_Parameter	_BayesianGraphicalModel::ComputeDiscreteScore (long node_id, _SimpleList & parents)
{
	//	Currently, discrete nodes may only have discrete parents, this should be enforced 
	//		via wrapper functions (above) and/or constraint matrix.
	//	---------------------------------------------------------------------------------
	
	// use cached node scores if available
	
	ReportWarning (_String ("Called ComputeDiscreteScore with ") & node_id);
	
	
	if (scores_cached)
	{
		_List *		scores	= (_List *) node_score_cache.lData[node_id];
		
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
	
	
	
	// impute score if missing data
	if (has_missing.lData[node_id])
	{
		return (GibbsImputeScore (node_id, parents));
	}
	else
	{
		for (long par = 0; par < parents.lLength; par++)
		{
			if (has_missing.lData[parents.lData[par]])
			{
				return (GibbsImputeScore (node_id, parents));
			}
		}
	}
	
	
	
	// compute score for complete data
	_SimpleList		multipliers ((long)1);
	
	_Matrix			n_ijk,
					n_ij;		// [i] indexes child nodes, 
								// [j] indexes combinations of values for parents of i-th node, 
								// [k] indexes values of i-th node
	
	long			num_parent_combos	= 1,					// i.e. 'q'
					r_i					= num_levels.lData[node_id];
	
	
	
	
	// how many combinations of parental states are there?
	for (long par = 0; par < parents.lLength; par++)
	{
		num_parent_combos *= num_levels.lData[parents.lData[par]];
		multipliers << num_parent_combos;
	}
	
	
	
	// count observations by parent combination, using direct indexing
	CreateMatrix (&n_ijk, num_parent_combos, r_i, false, true, false);
	CreateMatrix (&n_ij, num_parent_combos, 1, false, true, false);
	
	
	for (long obs = 0; obs < theData->GetHDim(); obs++)
	{
		long	index			= 0,
				child_state		= (*theData)(obs, node_id);
		
		for (long par = 0; par < parents.lLength; par++)
		{
			long	this_parent			= parents.lData[par],
					this_parent_state	= (*theData)(obs, this_parent);
			
			index += this_parent_state * multipliers.lData[par];
		}
		
		n_ijk.Store ((long) index, child_state, n_ijk(index, child_state) + 1);
		n_ij.Store ((long) index, 0, n_ij(index, 0) + 1);
	}
	
	
	
	return ( (prior_sample_size (node_id, 0) == 0) ? 
				K2Score (node_id, n_ij, n_ijk) : 
				BDeScore (node_id, n_ij, n_ijk) );
}




//___________________________________________________________________________________________
_Parameter _BayesianGraphicalModel::K2Score (long node_id, _Matrix & n_ij, _Matrix & n_ijk)
{
	_Parameter	log_score	= 0.;
	long		r_i			= num_levels.lData[node_id];
	
	for (long j = 0; j < n_ij.GetHDim(); j++)
	{
		log_score += LnGamma(r_i);	// (r-1)!
		log_score -= LnGamma(n_ij(j, 0) + r_i);	// (N+r-1)!
		
		for (long k = 0; k < r_i; k++)
		{
			log_score += LnGamma(n_ijk(j,k) + 1);	// (N_ijk)!
		}
	}
	
	return log_score;
}



//___________________________________________________________________________________________
_Parameter _BayesianGraphicalModel::BDeScore (long node_id, _Matrix & n_ij, _Matrix & n_ijk)
{
	_Parameter	n_prior_ij		= prior_sample_size (node_id, 0) / n_ij.GetHDim(),
				n_prior_ijk		= n_prior_ij / num_levels.lData[node_id],
				log_score		= 0.;
	
	for (long j = 0; j < n_ij.GetHDim(); j++)
	{
		log_score += LnGamma(n_prior_ij) - LnGamma(n_prior_ij + n_ij(j,0));
		
		for (long k = 0; k < num_levels.lData[node_id]; k++)
		{
			log_score += LnGamma(n_prior_ijk + n_ijk(j,k)) - LnGamma(n_prior_ijk);
		}
	}
	
	return log_score;
}
 


//___________________________________________________________________________________________
void	_BayesianGraphicalModel::InitMarginalVectors (_List * compute_list)
{
	_GrowingVector * newstore;
	checkPointer (newstore = new _GrowingVector);
	
	for (long i = 0; i < num_nodes * num_nodes; i++)
		(*compute_list) && newstore;
	
	DeleteObject (newstore);
}



//___________________________________________________________________________________________
void	_BayesianGraphicalModel::DumpMarginalVectors (_List * compute_list)
{	
	for (long i = 0; i < compute_list->lLength; i++)
		((_GrowingVector *) compute_list->lData[i]) -> Clear();
	
	compute_list->Clear();
}




//___________________________________________________________________________________________
void	_BayesianGraphicalModel::CacheNodeScores (void)
{
	ReportWarning (_String ("Entered CacheNodeScores()"));
	
	if (scores_cached)
		return;
	
	
#if defined __AFYP_DEVELOPMENT__ && defined __HYPHYMPI__
	// MPI_Init() is called in main()
	int			size,
				rank;
	
	long		mpi_node;
	
	_Matrix		single_parent_scores (num_nodes, 1, false, true);
	
	_SimpleList	parents,
				all_but_one (num_nodes-1, 0, 1),
				aux_list,
				nk_tuple;
	
	_Parameter	score;
	
	char		mpi_message [256];
	
	MPI_Status	status;	// contains source, tag, and error code
	
	
	
	MPI_Comm_size (MPI_COMM_WORLD, &size);	// number of processes
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	
	if (rank == 0)
	{
		_String		bgmSwitch ("_BGM_SWITCH_"),
					bgmStr;
		
		_List	*	this_list;
		
		SerializeBgm (bgmStr);
		
		_Matrix		* mpi_node_status = new _Matrix ((long)size, 2, false, true);
		long		senderID;
		
		
		
		// switch compute nodes from mpiNormal to bgm loop context
		for (long ni = 1; ni < size; ni++)
		{
			MPISendString (bgmSwitch, ni);
		}
		
		
		
		// receive confirmation of successful switch
		for (long ni = 1; ni < size; ni++)
		{
			long fromNode = ni;
			
			_String t (MPIRecvString (ni,fromNode));
			if (!t.Equal (&bgmSwitch))
			{
				WarnError (_String("Failed to confirm MPI mode switch at node ") & ni);
				return;
			}
			else
			{
				ReportWarning (_String("Successful mode switch to Bgm MPI confirmed from node ") & ni);
				MPISendString (bgmStr, ni);
			}
		}
		
		
		
		// farm out jobs to idle nodes until none are left
		for (long node_id = 0; node_id < num_nodes; node_id++)
		{
			long		maxp			= max_parents.lData[node_id],
			ntuple_receipt,
			this_node;
			
			bool		remaining;
			
			_String		mxString,
			mxName;
			
			
			this_list	= (_List *) node_scores.lData[node_id];
			
			// [_SimpleList parents] should always be empty here
			_Parameter	score;
			
			if (data_type.lData[node_id] == 0)
			{
				ComputeDiscreteScore (node_id, parents);
			}
			else
			{
				ComputeContinuousScore (node_id, parents);
			}
			
			_Constant	orphan_score (score);
			
			
			this_list->Clear();
			(*this_list) && (&orphan_score);	// handle orphan score locally
			
			
			if (maxp > 0)	// don't bother to farm out trivial cases
			{
				// look for idle nodes
				mpi_node = 1;
				do
				{
					if ((*mpi_node_status)(mpi_node, 0) == 0)
					{
						ReportMPIError(MPI_Send(&node_id, 1, MPI_LONG, mpi_node, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD), true);
						ReportWarning (_String ("Sent child ") & node_id & " to node " & mpi_node);
						
						mpi_node_status->Store (mpi_node, 0, 1);	// set busy signal
						mpi_node_status->Store (mpi_node, 1, node_id);
						
						break;
					}
					mpi_node++;
				}
				while (mpi_node < size);
			}
			
			
			if (mpi_node == size)	// all nodes are busy, wait for one to send results
			{
				MPIReceiveScores (mpi_node_status, true, node_id);
			}
		}
		// end loop over child nodes
		
		
		// collect remaining jobs
		while (1)
		{
			// look for a busy node
			mpi_node = 1;
			do
			{
				if ( (*mpi_node_status)(mpi_node,0) == 1)
				{
					break;
				}
				mpi_node++;
			}
			while (mpi_node < size);
			
			
			if (mpi_node < size)	// at least one node is busy
			{
				MPIReceiveScores (mpi_node_status, false, node_id);
			}
			else
			{
				break;
			}
		}
		
		
		// shut down compute nodes
		for (long shutdown = -1, mpi_node = 1; mpi_node < size; mpi_node++)
		{
			ReportMPIError(MPI_Send(&shutdown, 1, MPI_LONG, mpi_node, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD), true);
			ReportWarning (_String ("Node 0 sending shutdown signal to node ") & mpi_node);
		}
		
	}
	
	else	// compute node
	{
		long		node_id,
		maxp;
		
		_List		list_of_matrices;
		
		while (1)
		{
			_String		mxString,
						mxName;
			
			list_of_matrices.Clear();
			
			
			// wait for master node to issue node ID to compute
			ReportMPIError (MPI_Recv (&node_id, 1, MPI_LONG, 0, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD, &status), false);
			ReportWarning (_String("Node") & (long)rank & " received child " & (long)node_id & " from node " & (long)status.MPI_SOURCE);
			
			if (node_id < 0)
			{
				break;	// received shutdown message (-1)
			}
			
			
			maxp = max_parents.lData[node_id];
			
			parents.Clear();
			parents.Populate (1,0,0);
			
			
			// compute single parent scores
			for (long par = 0; par < num_nodes; par++)
			{
				if (par == node_id)		// child cannot be its own parent, except in Kansas
				{
					single_parent_scores.Store (par, 0, 0.);
				}
				else
				{
					parents.lData[0] = par;
					single_parent_scores.Store (par, 0, data_type.lData[node_id] ? ComputeDiscreteScore (node_id, parents) : 
												ComputeContinuousScore (node_id, parents));
				}
			}
			
			
			// compute multiple parents cores
			for (long np = 2; np <= maxp; np++)
			{
				parents.Clear();
				parents.Populate (np, 0, 0);
				
				if (all_but_one.NChooseKInit (aux_list, nk_tuple, np, false))
				{
					bool		remaining;
					long		tuple_index		= 0,
					num_ktuples		= exp(LnGamma(num_nodes) - LnGamma(num_nodes - np) - LnGamma(np+1));
					
					
					_Matrix		tuple_scores (num_ktuples, 1, false, true);
					
					for (long tuple_index = 0; tuple_index < num_ktuples; tuple_index++)
					{
						remaining = all_but_one.NChooseK (aux_list, nk_tuple);
						
						if (!remaining && tuple_index < num_ktuples-1)
						{
							ReportWarning (_String ("ERROR: Ran out of (n,k)tuples in CacheNodeScores()."));
						}
						
						for (long par_idx = 0; par_idx < np; par_idx++)
						{
							long par = nk_tuple.lData[par_idx];
							if (par >= node_id) par++;
							parents.lData[par_idx] = par;
						}
						
						score = data_type.lData[node_id] ?	ComputeDiscreteScore (node_id, parents) : 
															ComputeContinuousScore (node_id, parents);
						
						tuple_scores.Store (tuple_index, 0, (double)score);
					}
					
					if (remaining)
					{
						ReportWarning (_String ("ERROR: Did not compute all nk-tuples in CacheNodeScores()"));
					}
					
					list_of_matrices && (&tuple_scores);		// append duplicate
				}
				else 
				{
					ReportWarning (_String ("Failed to initialize _NTupleStorage object in Bgm::CacheNodeScores()"));
				}
				
			}
			
			// send results to master node
			
			ReportMPIError (MPI_Send (single_parent_scores.theData, num_nodes, MPI_DOUBLE, 0, HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD), true);
			
			
			for (long np = 2; np <= maxp; np++)
			{
				_Matrix	* storedMx		= (_Matrix *) list_of_matrices.lData[np-2];
				
				ReportMPIError (MPI_Send (storedMx->theData, storedMx->GetHDim(), MPI_DOUBLE, 0, HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD), true);
			}
		}
	}
	
#else
	
		
#if !defined __UNIX__ || defined __HEADLESS__
	TimerDifferenceFunction(false); // save initial timer; will only update every 1 second
#if !defined __HEADLESS__
	SetStatusLine 	  (empty,_HYBgm_STATUS_LINE_CACHE, empty, 0, HY_SL_TASK|HY_SL_PERCENT);
#else
	SetStatusLine 	  (_HYBgm_STATUS_LINE_CACHE);
#endif	
	_Parameter	seconds_accumulator = .0,
	temp;
#endif
	
	for (long node_id = 0; node_id < num_nodes; node_id++)
	{
		long		maxp		= max_parents.lData[node_id];
		_List	*	this_list	= (_List *) node_score_cache.lData[node_id];	// retrieve pointer to list of scores for this child node
		
		this_list->Clear();		// reset the list
		
		// prepare some containers
		_SimpleList	parents;
		_Parameter	score;
		
		if (data_type.lData[node_id] == 0)
		{
			score = ComputeDiscreteScore (node_id, parents);
		}
		else 
		{
			score = ComputeContinuousScore (node_id, parents);
		}
		
		_Constant	orphan_score (score);	// score computed with no parents (initialized empty list)
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
				
				if (data_type.lData[node_id] == 0)
				{
					score = ComputeDiscreteScore (node_id, parents);
				}
				else 
				{
					score = ComputeContinuousScore (node_id, parents);
				}
				
				single_parent_scores.Store (par, 0, score);
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
						
						if (data_type.lData[node_id] == 0)
						{
							score = ComputeDiscreteScore (node_id, parents);
						}
						else 
						{
							score = ComputeContinuousScore (node_id, parents);
						}
						
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
			
#if defined __HEADLESS__
			SetStatusLine (statusLine);
#else
			SetStatusLine (empty,statusLine,empty,100*(float)node_id/(num_nodes),HY_SL_TASK|HY_SL_PERCENT);
#endif
			TimerDifferenceFunction (false); // reset timer for the next second
			yieldCPUTime (); // let the GUI handle user actions
			
			if (terminateExecution) // user wants to cancel the analysis
				break;
		}
#endif
	}
#endif

#if !defined __UNIX__ || defined __HEADLESS__
	SetStatusLine 	  (_HYBgm_STATUS_LINE_CACHE_DONE);
#endif
	
	scores_cached = TRUE;
	
	ReportWarning (_String ("Cached node scores."));
}



//________________________________________________________________________________________________________
#if defined __AFYP_DEVELOPMENT__ && defined __HYPHYMPI__
void _BayesianGraphicalModel::MPIReceiveScores (_Matrix * mpi_node_status, bool sendNextJob, long node_id)
{
	_Matrix		single_parent_scores (num_nodes, 1, false, true);
	MPI_Status	status;
	
	ReportMPIError (MPI_Recv (single_parent_scores.theData, num_nodes, MPI_DOUBLE, MPI_ANY_SOURCE, 
							  HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD, &status), false);
	
	
	long		senderID	= (long) status.MPI_SOURCE,
				this_node	= (long) (*mpi_node_status) (senderID, 1),
				maxp		= max_parents.lData[this_node];
	
	_List	*	this_list	= (_List *) node_scores.lData[this_node];
	
	
	_String		mxString,
				mxName;
	
	
	mpi_node_status->Store (senderID, 0, 0);	// set node status to idle
	ReportWarning (_String("Received scores for child ") & this_node & " from node " & senderID);
	
	(*this_list) && (&single_parent_scores);
	
	
	_SimpleList		parents,
					all_but_one (num_nodes-1, 0, 1),
					aux_list,
					nk_tuple;
	
	_Parameter		score;
	
	// parse nk-tuple results
	for (long np = 2; np <= maxp; np++)
	{
		_NTupleStorage	family_scores (num_nodes-1, np);
		bool			remaining;
		
		
		if (all_but_one.NChooseKInit (aux_list, nk_tuple, np, false))
		{
			long		score_index		= 0,
			num_ktuples		= exp(LnGamma(num_nodes) - LnGamma(num_nodes - np) - LnGamma(np+1)),
			ntuple_receipt;
			
			_Matrix		scores_to_store (num_ktuples, 1, false, true);
			
			// receive nk-tuple indexed node scores from same compute node
			ReportMPIError (MPI_Recv (scores_to_store.theData, num_ktuples, MPI_DOUBLE, senderID, 
									  HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD, &status), false);
			
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
	
	
	// send next job
	if (sendNextJob)
	{
		ReportMPIError(MPI_Send(&node_id, 1, MPI_LONG, senderID, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD), true);
		mpi_node_status->Store (senderID, 0, 1);	// reset busy signal
		mpi_node_status->Store (senderID, 1, node_id);
		ReportWarning (_String ("Sent child ") & node_id & " to node " & senderID);
	}
}


void _BayesianGraphicalModel::SerializeBgm (_String & rec)
{
	_String	*	bgmName = (_String *) bgmNamesList (bgmList._SimpleList::Find((long)this));
	_String		dataStr,
				dataName ("bgmData");
	
	_Parameter	mcem_max_steps, mcem_burnin, mcem_sample_size;	// permit user to reset impute settings
	
	checkParameter (mcemMaxSteps, mcem_max_steps, 10000);
	checkParameter (mcemBurnin, mcem_burnin, 1000);
	checkParameter (mcemSampleSize, mcem_sample_size, 100);
	
	// write utility functions
	rec << "function make_dnode (id,n,maxp)\n";
	rec << "{\ndnode={};\n";
	rec << "dnode[\"NodeID\"]=id;\n";
	rec << "dnode[\"PriorSize\"]=n;\n";
	rec << "dnode[\"MaxParents\"]=maxp;\n";
	rec << "dnode[\"PriorMean\"]=0;\n";
	rec << "dnode[\"PriorVar\"]=0;\n";
	rec << "return dnode;\n}\n";
	
	rec << "function make_cnode (id,n,maxp,mean,prec)\n";
	rec << "{\ncnode={};\n";
	rec << "cnode[\"NodeID\"]=id;\n";
	rec << "cnode[\"PriorSize\"]=n;\n";
	rec << "cnode[\"MaxParents\"]=maxp;\n";
	rec << "cnode[\"PriorMean\"]=mean;\n";
	rec << "cnode[\"PriorVar\"]=prec;\n";
	rec << "return cnode;\n}\n";
	
	
	// prepare assortative lists
	rec << "dnodes={};\ncnodes={};\n";
	
	for (long node_id = 0; node_id < num_nodes; node_id++)
	{
		if (is_discrete.lData[node_id])
		{
			rec << "nodes[Abs(nodes)]=make_dnode(";
			sprintf (buf, "%d,%d,%d", node_id, (long)prior_sample_size(node_id,0), (long)max_parents.lData[node_id]);
			rec << buf;
			rec << ");\n";
		}
		else
		{
			rec << "nodes[Abs(nodes)]=make_cnode(";
			sprintf (buf, "%d,%d,%d,%f,%f", node_id, (long)prior_sample_size(node_id,0), (long)max_parents.lData[node_id],
					 prior_mean(node_id,0), prior_precision(node_id,0));
			rec << buf;
			rec << ");\n";
		}
	}
	
	// write BGM constructor
	rec << "BGM ";
	rec << bgmName;
	rec << "=(nodes);\n";
	
	// missing data imputation settings
	rec << "BGM_MCEM_MAXSTEPS = ";
	rec << (long) mcem_max_steps;
	rec << ";\nBGM_MCEM_BURNIN = ";
	rec << (long) mcem_burnin;
	rec << ";\nBGM_MCEM_SAMPLES = ";
	rec << ";\n";
	
	// serialize data matrix and assign to BGM
	obsData->Serialize (dataStr, dataName);
	rec << dataStr;
	rec << "\n";
	rec << "SetParameter(";
	rec << bgmName;
	rec << ",BGM_DATA_MATRIX,bgmData);\n";
}
#endif



//_________________________________________________________________________________________________________
_Matrix *	_BayesianGraphicalModel::Optimize (void)
{
	//	Wrapper function to access various optimization methods (K2, structural MCMC, order MCMC)
	//		and implement MPI functionality.
	
	if (!scores_cached) 
	{
		CacheNodeScores();
	}
	
	
	_Parameter		optMethod;		// 0 = K2 fixed order
									// 1 = K2 shuffle order with restarts 
									// 2 = structural MCMC, fixed order
									// 3 = structural MCMC
									// 4 = order MCMC
	
	_Matrix *	output_matrix;
	
	
	checkParameter (_HYBgm_METHOD_KEY, optMethod, 0.);
	
	if (optMethod < 2)
	{
		_Parameter	num_restarts,			// HBL settings
					num_randomize;
		
		checkParameter (_HYBgm_K2_RESTARTS, num_restarts, 1.);
		checkParameter (_HYBgm_K2_RANDOMIZE, num_randomize, num_nodes);
		
		checkPointer (output_matrix =  new _Matrix (num_nodes * num_nodes, 2, false, true));
		K2Search (optMethod, num_restarts, num_randomize, output_matrix);
	}
	else 
	{
		_String			oops;
		_Parameter		mcmc_steps, mcmc_burnin, mcmc_samples,
						mcmc_nchains, mcmc_dtemp;
		
		
		// acquisition of HBL arguments with a few sanity checks
		checkParameter (_HYBgm_MCMC_MAXSTEPS, mcmc_steps, 0);
		if (mcmc_steps <= 0)
		{
			oops = _String ("You asked HyPhy to run MCMC with zero steps in the chain! Did you forget to set Bgm_MCMC_STEPS?\n");
			WarnError (oops);
		}
		
		checkParameter (_HYBgm_MCMC_BURNIN, mcmc_burnin, 0);
		if (mcmc_burnin < 0)
		{
			oops = _String("You can't have a negative burn-in (_HYBgm_MCMC_BURNIN)!\n");
			WarnError (oops);
		}
		
		checkParameter (_HYBgm_MCMC_SAMPLES, mcmc_samples, 0);
		if (mcmc_samples < 0)
		{
			oops = _String ("You can't have a negative sample size!");
			WarnError (oops);
		}
		
#if defined __AFYP_DEVELOPMENT__ && defined __HYPHYMPI__
		int			size,
					rank;
		long		mpi_node;
		char		mpi_message [256];
		
		MPI_Status	status;	// contains source, tag, and error code
		
		MPI_Comm_size (MPI_COMM_WORLD, &size);	// number of processes
		MPI_Comm_rank (MPI_COMM_WORLD, &rank);
		
		
		checkParameter (_HYBgm_MCMC_NCHAINS, mcmc_nchains, 1);
		if (mcmc_nchains <= 0)
		{
			oops = _String ("You must run at least one chain in MCMC.");
			WarnError (oops);
			
			return nil;
		}
		
		
		
		// parallel coupled chains (MPI)
		if (mcmc_nchains > 1)
		{
			checkParameter (_HYBgm_MCMC_TEMP, mcmc_dtemp, 1);
			if (mcmc_nchains > 0 && mcmc_dtemp <= 1)
			{
				oops = _String ("Temperature increment must be greater than 1.");
				WarnError (oops);
			}
			
			_Parameter	chain_T = 1. / (1. + mcmc_dtemp*rank);
			
			
			/* note that GraphMetropolis already returns structure and last posterior prob as matrix entries */
			/* use sampling interval to swap chains */
				/* each chain needs a matrix pointer for output */
			if (optMethod < 4)
			{
				checkPointer (output_matrix = new _Matrix ( (mcmc_samples > num_nodes*num_nodes) ? mcmc_samples : num_nodes*num_nodes, 4, false, true));
				GraphMetropolis (optMethod==2, mcmc_burnin, mcmc_steps, mcmc_samples, chain_T, output_matrix);
				
				for (long samp = 0; samp < mcmc_samples; samp++)
				{
					/*** UNDER DEVELOPMENT ***/
					if (rank > 0)
					{
						
					}
					else
					{
						
					}
				}
			}
		}
		
		// single chain
		else
		{
#endif
			checkPointer (output_matrix = new _Matrix ( (mcmc_samples > num_nodes*num_nodes) ? mcmc_samples : num_nodes*num_nodes, 4, false, true));
			
			if (optMethod < 4)
			{
				GraphMetropolis ( (optMethod == 2), mcmc_burnin, mcmc_steps, mcmc_samples, 1., output_matrix);
			}
			else
			{
				OrderMetropolis (mcmc_burnin, mcmc_steps, mcmc_samples, 1., output_matrix);
			}
			
#if defined __AFYP_DEVELOPMENT__ && defined __HYPHYMPI__
		}
#endif
		
	}
	
	return (_Matrix *) (output_matrix->makeDynamic());
}



//_________________________________________________________________________________________________________
void _BayesianGraphicalModel::K2Search (bool do_permute_order, long n_restart, long n_randomize, _Matrix * result)
{
	//  THIS NEEDS UPDATING
#ifdef __NEVER_DEFINED__
	_Parameter	this_score, best_score, next_score;
	
	_Matrix		order_matrix (num_nodes, num_nodes, false, true),
				best_dag (num_nodes, num_nodes, false, true);
	
	_SimpleList	best_node_order;
	
	
	
	best_node_order.Populate (num_nodes, 0, 1);
	best_node_order.Permute (1);
	
	
	
	//	Convert node order to binary matrix where edge A->B is permitted if
	//	order_matrix[B][A] = 1, i.e. B is to the right of A in node order.
	for (long i = 0; i < num_nodes; i++)
	{
		long	child = best_node_order.lData[i];
		
		for (long j = 0; j < num_nodes; j++)
		{
			long	parent = best_node_order.lData[j];
			
			order_matrix.Store (parent, child, (j > i) ? 1 : 0);
		}
	}
	
	
	// greedy hill-climbing algorithm (K2)
	// ResetGraph (nil);
	best_score = Compute();		// best over all node orders (if we're reshuffling)
	
	
	for (long iter = 0; iter < n_restart; iter++)
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
						&& (theStructure(parent, child) == 0) 
						&& (order_matrix (parent, child) == 1) 
						&& constraint_graph(parent, child) >= 0)
					{
						theStructure.Store (parent, child, 1);
						this_score = Compute();
						
						if (this_score > next_score)
						{
							improvement_flag	= 1;
							next_score			= this_score;
							next_parent_to_add	= parent;
						}
						theStructure.Store (parent, child, 0);	// revert
					}
				}
				
				if (improvement_flag)	// adding another parent improves network score
				{
					theStructure.Store (next_parent_to_add, child, 1);
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
		
		
		if (do_permute_order)
		{
			this_score = Compute();
			
			if (this_score > best_score)
			{
				best_score = this_score;
				best_dag = (_Matrix &) theStructure;		// store graph optimized from the current ordering
			}
			
			// ResetGraph (nil);
			
			best_node_order.Permute (1);
			
			for (long i = 0; i < num_nodes; i++)
			{
				long	child = best_node_order.lData[i];
				for (long j = 0; j < num_nodes; j++)
				{
					long	parent = best_node_order.lData[j];
					order_matrix.Store (parent, child, (j > i) ? 1 : 0);
				}
			}
		}
		else
		{
			break;	// only one iteration when node ordering is known
		}
	}
	
	if (do_permute_order)
	{
		theStructure = (_Matrix &) best_dag;
		best_node_order.Clear();	// reset to initial state
	}
	
	
	result->Store (0, 0, Compute());
	
	for (long row = 0; row < num_nodes; row++)
	{
		for (long col = 0; col < num_nodes; col++)
		{
			result->Store (row*num_nodes+col, 1, theStructure(row, col));
		}
	}
#endif
}



//_______________________________________________________________________________________________________________________
bool	_BayesianGraphicalModel::GraphObeysOrder (_Matrix & graph, _SimpleList & order)
{
	_Matrix order_matrix (num_nodes, num_nodes, false, true);
	
	// convert order vector to matrix form
	for (long p_index = 0; p_index < num_nodes; p_index++)
	{
		for (long par = order.lData[p_index], c_index = 0; c_index < num_nodes; c_index++)
		{
			order_matrix.Store (par, order.lData[c_index], (p_index > c_index) ? 1 : 0);
		}
	}
	
	// loop thru graph, check that edges are compatible with node order
	for (long parent = 0; parent < num_nodes; parent++)
	{
		for (long child = 0; child < num_nodes; child++)
		{
			if (graph(parent, child)==1 && order_matrix(parent, child)==0)
			{
				return 0;
			}
		}
	}
	
	return 1;
}



//_______________________________________________________________________________________________________________________
_SimpleList *	_BayesianGraphicalModel::GetOrderFromGraph (_Matrix & graph)
{
	//  To quickly generate a node order based on graph argument.
	//	Loop through nodes in graph and rank according to parentage.
	
	_SimpleList	*	new_order  = new _SimpleList (1, 0, 0);	// initialize with single entry, [0]
	
	for (long left_of, node = 1; node < num_nodes; node++)
	{
		// loop through nodes in list looking for parents
		for (left_of = 0; left_of < new_order->lLength; left_of++)
		{
			if (graph (left_of, node))
			{
				new_order->InsertElement ((BaseRef) node, left_of);
				break;
			}
		}
		
		// if we reach end of list, append
		if (left_of == new_order->lLength)  (*new_order) << node;
		
	}
	
	return (_SimpleList *) new_order->makeDynamic();
}



//_______________________________________________________________________________________________________________________
void	_BayesianGraphicalModel::GraphMetropolis (bool fixed_order, long mcmc_burnin, long mcmc_steps, long mcmc_samples, 
												  _Parameter chain_t, _Matrix * result)
{
	//	Performs MCMC over structures.  Initialize chain using class member [theStructure] (accessible by HBL).
	//		If [fixed_order]=TRUE, only explores the subset of graphs compatible with [node_order].
	//	Returns pointer to matrix containing:	(1) vector of posterior probabiities
	//											(2) linearized matrix of edge marginal posteriors
	//											(3)		"		"	for best graph
	//											(4)		"		"	for last graph visited in chain
	
	
	_Matrix		*	proposed_graph	= new _Matrix (num_nodes, num_nodes, false, true);
	
	_Matrix			current_graph (num_nodes, num_nodes, false, true),
					best_graph (num_nodes, num_nodes, false, true);
	
	_Parameter		current_score, proposed_score, best_score,
					lk_ratio,
					prob_swap;
	
	long			sampling_interval = mcmc_steps / mcmc_samples,
					max_fails;
	
	_SimpleList	*	proposed_order	= new _SimpleList();
	_SimpleList		current_order;
	
	
	// parse HBL settings
	checkParameter (_HYBgm_MCMC_PROBSWAP, prob_swap, 0);
	if (prob_swap < 0 || prob_swap > 1.)
	{
		_String oops ("Bgm_PROB_SWAP must be assigned a value between 0 and 1.  Exiting.\n");
		WarnError (oops);
	}
	
	
	
	
	if (fixed_order && node_order.lLength == 0)
	{
		_String	oops ("ERROR: Missing node order when attempting structural-MCMC with fixed order argument.");
		WarnError (oops);
	}
	
	
	// initialize state
	if (GraphObeysOrder (theStructure, node_order))
	{
		(*proposed_graph) = (_Matrix &) theStructure;
		(*proposed_order) = (_SimpleList &) node_order;
	}
	else
	{
		_String	oops ("ERROR: Structure does not match order, aborting GraphMetropolis().");
		WarnError (oops);
	}
	
	
	// status line
#if !defined __UNIX__ || defined __HEADLESS__
	long	updates = 0;
	TimerDifferenceFunction (false);
	#if defined __HEADLESS__
	SetStatusLine (_HYBgm_STATUS_LINE_MCMC);
	#else	
	SetStatusLine (empty, _HYBgm_STATUS_LINE_MCMC, empty, 0, HY_SL_TASK|HY_SL_PERCENT);
	#endif
#endif
	
	
	current_order = (*proposed_order);
	
	current_graph = (_Matrix &) (*proposed_graph);
	best_graph = (_Matrix &) (*proposed_graph);
	
	best_score = proposed_score = current_score = Compute((_Matrix &) *proposed_graph);
	
	
	// main loop
	for (long step = 0; step < mcmc_steps + mcmc_burnin; step++)
	{
		RandomizeGraph (proposed_graph, proposed_order, prob_swap, 1, max_fails, fixed_order);
		
		proposed_score = Compute((_Matrix &) *proposed_graph);
		
		lk_ratio = exp(proposed_score - current_score);
		
		if (lk_ratio > 1. || genrand_real2() < lk_ratio)		// Metropolis-Hastings
		{
			current_graph	= (_Matrix &) (*proposed_graph);		// accept proposed graph
			current_score	= proposed_score;
			
			for (long foo = 0; foo < num_nodes; foo++)
			{
				current_order.lData[foo] = proposed_order->lData[foo];
			}
			
			if (current_score > best_score)
			{
				best_score = current_score;
				best_graph = (_Matrix &) current_graph;			// keep track of best graph ever visited by chain
			}
		}
		else
		{
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
		
		// chain sampling
		if (step >= mcmc_burnin)
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
	#if defined __HEADLESS__
			SetStatusLine 	  (statusLine);
	#else	
			SetStatusLine 	  (empty,statusLine,empty,100*step/(mcmc_steps + mcmc_burnin),HY_SL_TASK|HY_SL_PERCENT);
	#endif
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
	
	
	// release memory
	delete (proposed_graph);
	delete (proposed_order);
}



//_____________________________________________________________________________________________________________
void	_BayesianGraphicalModel::RandomizeGraph (_Matrix * graph, _SimpleList * order, _Parameter prob_swap, long num_steps, long max_fails, bool fixed_order)
{
	//	Modify a graph by adding, removing, or reversing an edge, so long as that modification
	//	complies with the constraint matrix and (when fixed_order=TRUE) node order.
	
	long	step = 0, fail = 0, 
			child, parent, 
			child_idx, parent_idx;
	
	
	// convert order into matrix format of edge permissions
	_Matrix order_matrix (num_nodes, num_nodes, false, true);
	
	for (long p_index = 0; p_index < num_nodes; p_index++)
	{
		for (long par = order->lData[p_index], c_index = 0; c_index < num_nodes; c_index++)
		{
			order_matrix.Store (par, order->lData[c_index], (p_index > c_index) ? 1 : 0);
		}
	}
	
	
	
	// calculate number of parents for each child for rapid access later
	_SimpleList		num_parents;
	
	num_parents.Populate (num_nodes, 0, 0);
	
	for (child = 0; child < num_nodes; child++)
	{
		for (parent = 0; parent < num_nodes; parent++)
		{
			if ( (*graph)(parent, child) > 0)
			{
				num_parents.lData[child]++;
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
		if (fail > max_fails)
		{
			_String oops ("Failed to modify the graph in GraphMCMC() after max_fails attempts.");
			WarnError (oops);
			break;
		}
		
		// pick a random edge
		parent_idx	= (genrand_int32() % (num_nodes-1)) + 1;	// shift right to not include lowest node in order
		child_idx	= genrand_int32() % parent_idx;
		
		child = order->lData[child_idx];
		parent = order->lData[parent_idx];
		
		
		if (fixed_order || genrand_real2() > prob_swap)
		{
			if ( (*graph)(parent,child) == 0 && constraint_graph(parent,child) >= 0)		// add an edge
			{
				if (num_parents.lData[child] == max_parents.lData[child])
				{
					// move an edge from an existing parent to target parent
					long	deadbeat_dad	= genrand_int32() % max_parents.lData[child];
					long	social_worker	= 0;
					
					for (; social_worker < num_nodes; social_worker++)
					{
						if ( (*graph)(social_worker,child) == 1)
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
					
					if (constraint_graph(social_worker, child) <= 0)	// not an enforced edge
					{
						graph->Store (social_worker, child, 0.);
						num_parents.lData[child]--;
					}
					else
					{
						fail++;
					}
				}
				else
				{
					graph->Store (parent, child, 1.);
					num_parents.lData[child]++;
					step++;
				}
			}
			else if ( (*graph)(parent,child) == 1 && constraint_graph(parent,child) <= 0)		// delete an edge
			{
				graph->Store (parent, child, 0.);
				num_parents.lData[child]--;
				step++;
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
				if (constraint_graph(parent,bystander) > 0 || constraint_graph(bystander,child) > 0)
				{
					fail++;
					buzz = 1;
				}
			}
			
			if ( buzz == 0 && 
				( (*graph)(parent,child) == 0
				 || ( (*graph)(parent,child) == 1 && constraint_graph(parent,child) <= 0 
					 && num_parents.lData[parent] < max_parents.lData[parent])
				 ))
			{
				// flip the target edge
				if ( (*graph)(parent,child) == 1)
				{
					graph->Store (parent,child, 0);
					graph->Store (child,parent,1);
					num_parents.lData[child]--;
					num_parents.lData[parent]++;
				}
				
				// flip the other edges affected by node swap
				for (long bystander, i = child_idx+1; i < parent_idx; i++)
				{
					bystander = order->lData[i];
					if ( (*graph)(bystander, child) == 1)
					{
						graph->Store (bystander, child, 0);
						num_parents.lData[child]--;
						
						if (num_parents.lData[bystander] < max_parents.lData[bystander])
						{
							graph->Store (child, bystander, 1);
							num_parents.lData[bystander]++;
						}
					}
					
					if ( (*graph)(parent,bystander) == 1)
					{
						graph->Store (parent, bystander, 0);
						num_parents.lData[bystander]--;
						
						if (num_parents.lData[parent] < max_parents.lData[parent])
						{
							graph->Store (bystander, parent, 1);
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
						order_matrix.Store (par, order->lData[c_index], (p_index > c_index) ? 1 : 0);
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



//_______________________________________________________________________________________________________________________
void	_BayesianGraphicalModel::OrderMetropolis (long burn_in, long n_steps, long sample_size, _Parameter chain_t, 
												  _Matrix * result)
{
	
}



#endif



