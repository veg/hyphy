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
#if defined __AFYP_REWRITE_BGM__

#include "bayesgraph.h"

extern _String		_HYBgm_IMPUTE_MAXSTEPS,
					_HYBgm_IMPUTE_BURNIN,
					_HYBgm_IMPUTE_SAMPLES;

extern _Parameter	LnGamma (_Parameter),
					LogSumExpo (_GrowingVector *);

//___________________________________________________________________________________________
_Parameter	gaussianDeviate (void)
{
	/* 
		Use Box-Muller transform to generate random deviates from Gaussian distribution with
		zero mean and unit variance (Numerical Recipes).
	 */
	
	static int		iset = 0;
	static double	gset;
	double			fac, rsq, v1, v2;
	
	if (iset == 0)
	{
		do
		{
			v1 = 2.0 * genrand_real2() - 1.0;	// uniform random number on (0,1), i.e. no endpoints
			v2 = 2.0 * genrand_real2() - 1.0;
			rsq = v1*v1 + v2*v2;
		} 
		while (rsq >= 1.0 || rsq == 0.0);
		
		fac = sqrt(-2.0 * log(rsq)/rsq);
		gset = v1 * fac;
		iset = 1;			// set flag to indicate second deviate available
		
		return (_Parameter) (v2 * fac);
	} 
	else
	{
		iset = 0;
		return (_Parameter) gset;		// use second deviate
	}
}


//___________________________________________________________________________________________
_AssociativeList *	_BayesianGraphicalModel::ExportModel (void)
{
	return nil;
}


//___________________________________________________________________________________________
bool _BayesianGraphicalModel::ImportCache (_AssociativeList * cache_import)
{
	return 0;	// failed
}


//___________________________________________________________________________________________
bool _BayesianGraphicalModel::ExportCache (_AssociativeList * cache_export)
{
	
	
	// Export an associative list containing node_score_cache contents to HBL
	//	entries are indexed by key string referring to child node and number of parents, 
	//		e.g. "N0P2" refers to scores associated with child node 0 and two parents
	_String			keyString;
	
	_List			* this_list;
	_Constant		* orphan_score;
	_Matrix			* single_parent_scores;
	_NTupleStorage	* family_scores;
	
	if (scores_cached)
	{
		ReportWarning (_String("Exporting cache with ") & num_nodes & " nodes");
		
		for (long node = 0; node < num_nodes; node++)
		{
			this_list = (_List *) node_score_cache.lData[node];
			
			for (long npar = 0; npar <= max_parents.lData[node]; npar++)
			{
				keyString = _String ("N") & node & "P" & npar;
				_FString	aKey (keyString, false);
				
				ReportWarning (_String("Inserting with key ") & keyString);
				
				if (npar == 0)
				{
					orphan_score = (_Constant *) this_list->lData[npar];
					cache_export->MStore (&aKey, orphan_score, true);	// make a dynamic copy
				}
				else if (npar == 1)
				{
					single_parent_scores = (_Matrix *) this_list->lData[npar];
					cache_export->MStore (&aKey, single_parent_scores, true);	// make a dynamic copy
				}
				else
				{
					family_scores = (_NTupleStorage *) this_list->lData[npar];
					cache_export->MStore (&aKey, family_scores, true);	// make a dynamic copy					
				}
			}
		}
		
		return TRUE;
	}
	else
	{
		WarnError (_String ("Unable to export node score cache, no cache exists!"));
		return FALSE;
	}
}




//___________________________________________________________________________________________
// #define _DEBUG_CCS_
_Parameter _BayesianGraphicalModel::ComputeContinuousScore (long node_id)
{
	_SimpleList parents;
	
	// generate lists of discrete or continuous parents from graph
	for (long par = 0; par < num_nodes; par++)
	{
		if (theStructure(par, node_id) == 1)
		{
			parents << par;
		}
	}
	
	return ComputeContinuousScore (node_id, parents);
}


_Parameter _BayesianGraphicalModel::ComputeContinuousScore (long node_id, _Matrix &dag)
{
	_SimpleList		parents;
	
	for (long par = 0; par < num_nodes; par++)
	{
		// support for discrete parents only
		if (dag(par, node_id) == 1)
		{
			parents << par;
		}
	}
	return ComputeContinuousScore (node_id, parents);
}


//___________________________________________________________________________________________

_Parameter _BayesianGraphicalModel::ComputeContinuousScore (long node_id, _SimpleList & parents)
{
	/* ---------------------------------------------------------------------------
		ComputeContinuousScore()
	 
			Compute the conditional Gaussian network score of node [c] according to 
			S. Bottcher's formulation.
			
			Let [y] denote the set of other continuous nodes.
			Let [i] denote the set of discrete nodes.
			Let [pa] define a subset of continuous or discrete nodes on [c]
			
			Assume local probability distribution is Gaussian local regression on 
			continuous parents, with parameters conditional on discrete parent 
			configuration [pa_i].  This has parameters:
			
				beta (c | pa_i(c)) = vector of regression coefficients over pa_y (continuous parents)
				m (c | pa_i(c)) = regression intercept
				sigma^2 (c | pa_i(c)) = variance conditional on discrete parent configuration
		 
			Hence, for a case {y_1, y_2, ..., y_c-1, y_c+1, ..., y_m, i_1, i_2, ..., i_n},
				(y_c | pa_i(c) ~ Gaussian(m + beta_1 y_1 + ... + beta_m y_m, sigma^2)
			
			We need conjugate priors for these parameters, which are given by:
				(m, beta | sigma^2) ~ Gaussian (mu, sigma^2 tau^-1)
				(sigma^2) ~ Inverse-Gamma ( rho_{y|pa_i(y)} / 2, phi_{y|pa_i(y)} / 2)
			
			In other words, the intercept and regression coefficients of continuous 
			node (c) on continuous parents pa_y(c) given discrete parent configuration pa_i(y) 
		 
			We integrate over these priors given the hyperparameters:
				mu		- prior mean vector (intercept, regression coefficients...)
				tau		- prior precision
				rho		- degrees of freedom for inverse gamma conjugate prior
				phi		- scale parameter for inverse gamma conjugate prior
		 
			Typically, [rho] and [phi] are given small values (0.001) that are meant 
			to render an uninformative prior (flatter).  However, Gelman asserts that 
			this practice is not optimal.
		 
			[tau] is a precision matrix whose off-diagonal entries are usually zeroed 
			out.  Hence, the values set by an AVL upon constructing the BGM fill the 
			diagonal.  (I should probably implement some access functions so that
			these can be changed.)
	   ----------------------------------------------------------------------------- */
	
	ReportWarning (_String ("Called ComputeContinuousScore with ") & node_id & " <- " & (_String *) parents.toStr());
	
#ifdef _DEBUG_CCS_
	char debug [255];
	sprintf (debug, "\nComputeContinuousScore (%d,[", node_id);
	BufferToConsole (debug);
	
	for (long i = 0; i < parents.lLength; i++)
	{
		sprintf (debug, "%d,", parents.lData[i]);
		BufferToConsole (debug);
	}
	sprintf (debug, "])\n");
	BufferToConsole (debug);
	
	sprintf (debug, "prior precision: ");
	BufferToConsole (debug);
	for (long i = 0; i < num_nodes; i++)
	{
		sprintf (debug, "%1.3f ", prior_precision (i,0));
		BufferToConsole (debug);
	}
	NLToConsole ();
#endif
	
	/* WARNING, untested function! */
	/* --------------------------- */
	
	// mm = prior estimate of unconditional mean at continuous node (i.e. intercept)
	// phi = scale parameter of inverse gamma prior for variance, uninformative at low values 
	//						(Kohn, Smith, and Chan, 2001 Stat Comput 11: 313-322)
	
	_Parameter		log_score = 0.;
	
	long			num_parent_combos = 1,	// i.e. 'q'
					k;						// current number of continuous parents
				
	_Matrix			mu,
					tau;
	
	_SimpleList		multipliers ((long)1),
					c_parents,
					d_parents;
	
	
	// partition parent list into discrete and continuous
	for (long p = 0; p < parents.lLength; p++)
	{
		if (data_type.lData[parents.lData[p]] == 0)	// discrete
		{
			d_parents << parents.lData[p];
		}
		else
		{
			c_parents << parents.lData[p];
		}
	}
	
	k = c_parents.lLength;
	
	
	// how many combinations of parental states are there?
	for (long par = 0; par < d_parents.lLength; par++)
	{
		num_parent_combos *= num_levels.lData[d_parents.lData[par]];
		multipliers << num_parent_combos;
	}
	
	
	// set location hyperparameter for Gaussian prior
#ifdef _DEBUG_CCS_
	sprintf (debug, "\tlocation hyperparameter (mu): ");
	BufferToConsole (debug);
#endif
	CreateMatrix (&mu, k+1, 1, false, true, false);
	mu.Store (0, 0, prior_mean (node_id, 0));				// prior intercept
	for (long i = 1; i < mu.GetHDim(); i++)
	{
		mu.Store (i, 0, 0);		// set prior expectation of regression coefficients to zero
	}
#ifdef _DEBUG_CCS_
	for (long i = 0; i < mu.GetHDim(); i++)
	{
		sprintf (debug, "%1.3f ", mu(i,0));
		BufferToConsole (debug);
	}
	NLToConsole ();
#endif
	
	// set prior degrees of freedom (for inverse gamma / scaled inverse chi-square)
	//	and scale of variance prior
	_Parameter		rho	= prior_sample_size (node_id, 0) > 0 ? (prior_sample_size (node_id, 0) / num_parent_combos) : 1.0;
	_Parameter		phi = prior_scale (node_id, 0);
	
#ifdef _DEBUG_CCS_
	sprintf (debug, "\trho = %1.3f\n", rho);
	BufferToConsole (debug);
#endif
	
	// set precision hyperparameter for Gaussian prior
#ifdef _DEBUG_CCS_
	sprintf (debug, "\ttau (precision):\n\t\t");
	BufferToConsole (debug);
#endif
	CreateMatrix (&tau, k+1, k+1, false, true, false);		// for continuous nodes only, including child
	for (long row = 0; row < k+1; row++)
	{
		for (long col = 0; col < k+1; col++)
		{
			if (row == col)		// set diagonal entries
			{
				if (row == 0)
				{
					tau.Store (0, 0, prior_precision (node_id, 0));
				}
				else
				{
					tau.Store (row, col, prior_precision (c_parents.lData[row],0));
				}
			}
			
			else tau.Store (row, col, 0.);	// zero off-diagonal entries
#ifdef _DEBUG_CCS_		
			sprintf (debug, "%1.3f ", tau(row,col));
			BufferToConsole (debug);
		}
		sprintf (debug, "\n\t\t");
		BufferToConsole (debug);
	}
	NLToConsole ();
#else
		}
	}
#endif
	


	// count up number of data points per parent combination
	_SimpleList		n_ij,				// CAVEAT: other functions use _Matrix for this object
					pa_indexing;		// track discrete parent combinations per observation
	
	n_ij.Populate (num_parent_combos, 0, 0);
	pa_indexing.Populate (theData->GetHDim(), 0, 0);
	/*
	CreateMatrix (&n_ij, num_parent_combos, 1, false, true, false);		// why use matrices?  _SimpleList would be fine.. -afyp
	CreateMatrix (&pa_indexing, theData->GetHDim(), 1, false, true, false);
	*/

	for (long obs = 0; obs < theData->GetHDim(); obs++)
	{
		long	index		= 0,
				multiplier	= 1;
		
		for (long par = 0; par < d_parents.lLength; par++)
		{
			long	this_parent		= parents.lData[par];
			// max index = sum (parent * levels(parent))
			index += (*theData)(obs, this_parent) * multiplier;
			multiplier *= num_levels.lData[this_parent];
		}
		
		pa_indexing.lData[obs] = index;
		n_ij.lData[index] += 1;
	}
	
#ifdef _DEBUG_CCS_
	sprintf (debug, "\tn_ij: ");
	BufferToConsole (debug);
	for (long i = 0; i < n_ij.lLength; i++)
	{
		sprintf (debug, "%d ", n_ij.lData[i]);
		BufferToConsole (debug);
	}
	NLToConsole ();
	
	sprintf (debug, "\tpa_indexing: ");
	BufferToConsole (debug);
	for (long i = 0; i < pa_indexing.lLength; i++)
	{
		sprintf (debug, "%d ", pa_indexing.lData[i]);
		BufferToConsole (debug);
	}
	NLToConsole ();
#endif
	
	// for every parent combination, calculate contribution to score
	for (long count_n, pa = 0; pa < num_parent_combos; pa++)
	{
#ifdef _DEBUG_CCS_
		sprintf (debug, "parent combo %d of %d\n", pa, num_parent_combos);
		BufferToConsole (debug);
#endif
		
		_Matrix		zbpa (n_ij.lData[pa], k+1, false, true),
					yb (n_ij.lData[pa], 1, false, true);
		
		count_n		= 0;	
		// number of data points with this parent combination.
		
		
		// populate zbpa matrix with (1, y_1, y_2, ..., y_m) entries for this parental combo
		for (long obs = 0; obs < theData->GetHDim(); obs++)
		{
			if (pa_indexing.lData[obs] == pa)		// this observation has the current parent combo
			{									// load observed states at continuous parents into matrix
				//   I'm sure there's a faster way to do this! - afyp
				
				zbpa.Store (count_n, 0, 1);		// intercept
				
				for (long parent = 0; parent < k; parent++)
				{
					zbpa.Store (count_n, parent+1, (*theData)(obs, c_parents.lData[parent]));
				}
				
				// state vector for this continuous node
				yb.Store (count_n, 0, (*theData)(obs, node_id));
				
				count_n++;
			}
		}
		
#ifdef _DEBUG_CCS_	
		sprintf (debug, "\tzbpa:\n\t\t");
		BufferToConsole (debug);
		for (long i = 0; i < zbpa.GetHDim(); i++)
		{
			for (long j = 0; j < zbpa.GetVDim(); j++)
			{
				sprintf (debug, "%1.3f ", zbpa(i,j));
				BufferToConsole (debug);
			}
			sprintf (debug, "\n\t\t");
			BufferToConsole(debug);
		}
		NLToConsole();
		
		sprintf (debug, "\tyb: ");
		BufferToConsole (debug);
		for (long i = 0; i < yb.GetHDim(); i++)
		{
			sprintf (debug, "%1.3f ", yb(i,0));
			BufferToConsole (debug);
		}
		NLToConsole();
#endif
		
			
		log_score += BottcherScore (yb, zbpa, mu, tau, rho, phi, n_ij.lData[pa]);

	}
	
	return log_score;
}



//___________________________________________________________________________________________________

_Parameter	_BayesianGraphicalModel::BottcherScore (_Matrix & yb, _Matrix & zbpa, _Matrix & mu, _Matrix & tau, _Parameter rho, _Parameter phi, long batch_size)
{
	// calculate scale parameter for non-central t distribution
	// from S. Bottcher (2001) p.25
	_Parameter	pa_log_score = 0.;
	_Matrix		scale;
	
	scale = zbpa;
	scale *=  * (_Matrix *) tau.Inverse();	// n_ij x k+1
	zbpa.Transpose();	// k+1 x n_ij
	scale *= zbpa;		// n_ij x n_ij
	zbpa.Transpose();
	for (long row = 0; row < scale.GetHDim(); row++)	// add identity matrix
	{
		scale.Store (row, row, scale(row, row)+(_Parameter)1.);
	}
	scale *= (_Parameter) (phi / rho);
	
	
#ifdef _DEBUG_CCS_
	char debug [255];
	
	sprintf (debug, "\tscale:\n\t\t");
	BufferToConsole (debug);
	for (long i = 0; i < scale.GetHDim(); i++)
	{
		for (long j = 0; j < scale.GetVDim(); j++)
		{
			sprintf (debug, "%1.3f ", scale(i,j));
			BufferToConsole (debug);
		}
		sprintf (debug, "\n\t\t");
		BufferToConsole(debug);
	}
	NLToConsole();
#endif
	
	
	// calculate the determinant of scale parameter matrix
	_Matrix			temp_mat (scale);
	_Parameter		pi_const = 3.141592653589793;
	
	temp_mat *= (_Parameter) (pi_const * rho);
	
	_AssociativeList *	eigen		= (_AssociativeList *) temp_mat.Eigensystem();
	_Matrix *			eigenvalues = (_Matrix *)eigen->GetByKey(0, MATRIX);
	_Parameter			det			= 1.;
	
#ifdef _DEBUG_CCS_
	sprintf (debug, "\teigenvalues: ");
	BufferToConsole (debug);
	for (long i = 0; i < eigenvalues->GetHDim(); i++)
	{
		sprintf (debug, "%1.3f ", (*eigenvalues)(i,0));
		BufferToConsole (debug);
	}
	NLToConsole();
#endif
	
	// determinant is product of eigenvalues (should be > 0 for positive definite matrices)
	for (long i = 0; i < eigenvalues->GetHDim(); i++)
	{
		det *= (_Parameter)(*eigenvalues) (i,0);
	}
	
	
	
	
	// calculate first term of score
	
	pa_log_score += LnGamma((rho + batch_size)/2.);
	pa_log_score -= LnGamma(rho/2.) + 0.5 * log(det);
	
#ifdef _DEBUG_CCS_
	sprintf (debug, "\tdet: %1.3f\n", det);
	BufferToConsole (debug);
	sprintf (debug, "pa_log_score: %1.3f -> ", pa_log_score);
	BufferToConsole (debug);
	sprintf (debug, "%1.3f %1.3f %1.3f ", pa_log_score, LnGamma(rho/2.), log(det));
	BufferToConsole (debug);
	sprintf (debug, "%1.3f\n", pa_log_score);
	BufferToConsole (debug);
#endif
	
	// calculate second term of score
	_Matrix		next_mat;
	
	zbpa *= mu;		// n_ij x (k+1)  *  (k+1) x 1  --> n_ij x 1
	
#ifdef _DEBUG_CCS_
	sprintf (debug, "\tzbpa*mu:\n\t\t");
	BufferToConsole (debug);
	for (long i = 0; i < zbpa.GetHDim(); i++)
	{
		for (long j = 0; j < zbpa.GetVDim(); j++)
		{
			sprintf (debug, "%1.3f ", zbpa(i,j));
			BufferToConsole (debug);
		}
		sprintf (debug, "\n\t\t");
		BufferToConsole(debug);
	}
	NLToConsole();
#endif
	
	yb -= zbpa;		// n_ij x 1
	next_mat = yb;
	
	
#ifdef _DEBUG_CCS_
	sprintf (debug, "\t(yb-zbpa):\n\t\t");
	BufferToConsole (debug);
	for (long i = 0; i < next_mat.GetHDim(); i++)
	{
		for (long j = 0; j < next_mat.GetVDim(); j++)
		{
			sprintf (debug, "%1.3f ", next_mat(i,j));
			BufferToConsole (debug);
		}
		sprintf (debug, "\n\t\t");
		BufferToConsole(debug);
	}
	NLToConsole();
#endif
	
	
	temp_mat = * (_Matrix *) scale.Inverse();
	temp_mat *= next_mat;						// left multiply by (yb-zbpa)
	next_mat = temp_mat;
	
#ifdef _DEBUG_CCS_
	sprintf (debug, "\tnext_mat = scale^(-1) (yb-zbpa)^T :\n\t\t");
	BufferToConsole (debug);
	for (long i = 0; i < next_mat.GetHDim(); i++)
	{
		for (long j = 0; j < next_mat.GetVDim(); j++)
		{
			sprintf (debug, "%1.3f ", next_mat(i,j));
			BufferToConsole (debug);
		}
		sprintf (debug, "\n\t\t");
		BufferToConsole(debug);
	}
	NLToConsole();
#endif
	
	yb.Transpose();
	yb *= next_mat;
	next_mat = yb;
	
#ifdef _DEBUG_CCS_	
	sprintf (debug, "\tnext_mat = (yb-zbpa) scale^(-1) (yb_zbpa)^T:\n\t\t");
	BufferToConsole (debug);
	for (long i = 0; i < next_mat.GetHDim(); i++)
	{
		for (long j = 0; j < next_mat.GetVDim(); j++)
		{
			sprintf (debug, "%1.3f ", next_mat(i,j));
			BufferToConsole (debug);
		}
		sprintf (debug, "\n\t\t");
		BufferToConsole(debug);
	}
	NLToConsole();
#endif
	
	next_mat *= (_Parameter) (1./rho);	// should be a 1-element matrix
	
#ifdef _DEBUG_CCS_
	sprintf (debug, "\tnext_mat = (yb-zbpa) scale^(-1) (yb_zbpa)^T  / rho:\n\t\t");
	BufferToConsole (debug);
	for (long i = 0; i < next_mat.GetHDim(); i++)
	{
		for (long j = 0; j < next_mat.GetVDim(); j++)
		{
			sprintf (debug, "%1.3f ", next_mat(i,j));
			BufferToConsole (debug);
		}
		sprintf (debug, "\n\t\t");
		BufferToConsole(debug);
	}
	NLToConsole();
#endif
	
	pa_log_score += -(rho + batch_size)/2. * (next_mat(0,0) + 1.);
	
#ifdef _DEBUG_CSS_
	sprintf (debug, "pa_log_score = %f\n", pa_log_score);
	BufferToConsole (debug);
#endif
	
	return pa_log_score;
}



//___________________________________________________________________________________________________
_Parameter _BayesianGraphicalModel::ImputeNodeScore (long node_id, _SimpleList & parents)
{
	/* ---------------------------------------------------------------------------------------
		ImputeNodeScore ()
			Impute (i.e., "fill in") missing values in the data.  More exactly, take the
			expectation of node score over imputations of the data sampled using Gibbs
			sampling.  Missing values are flagged for discrete data by negative values.
			For continuous data, I haven't settled on an unambiguous flag and I'm leaving this
			up to the user for now (BGM_CONTINUOUS_MISSING_VALUE).
	   --------------------------------------------------------------------------------------- */
	
	long			num_parent_combos	= 1,
					r_i					= num_levels.lData[node_id],
					max_num_levels		= r_i,
					family_size			= parents.lLength + 1,
					sampling_interval;
	
	
	_SimpleList		multipliers (1, 1, 0),	// length constructor, populator
	
					dparents, cparents,
					parents_by_nodetype,	// indexes parents separately in mixed family (discete and CG)
					
					family_nlevels (family_size, 0, 0),		// for fast look-up
					family_datatype (family_size, 0, 0),
					
					is_missing,				// store linear indices to missing entries in deep copy
	
					pa_indices;
	

	_Matrix			n_ijk, n_ij,		// summary statistics for discrete nodes
	
					yb, zbpa,			// summary statistics for Gaussian nodes
					mu, tau,			// prior hyperparameters for CG nodes
	
					data_deep_copy,		// make a deep copy of the relevant columns of the data matrix
					observed_values;
	
	
	_GrowingVector	* vector_of_scores	= new _GrowingVector(),		// store scores sampled during imputation
					* reassign_probs	= new _GrowingVector();
	
	_Parameter		log_score			= 0,
	
					impute_maxsteps, impute_burnin, impute_samples,	// HBL settings
	
					parent_state, child_state,
					denom,
					this_prob,
														// prior hyperparameters for CG nodes
					rho	= prior_sample_size (node_id, 0) > 0 ? (prior_sample_size (node_id, 0) / num_parent_combos) : 1.0,
					phi = prior_scale (node_id, 0);
	
	
	double			urn;			// uniform random number
	
	
	
	// set Gibbs sampler parameters from batch language definitions
	checkParameter (_HYBgm_IMPUTE_MAXSTEPS, impute_maxsteps, 0);
	checkParameter (_HYBgm_IMPUTE_BURNIN, impute_burnin, 0);
	checkParameter (_HYBgm_IMPUTE_SAMPLES, impute_samples, 0);
	
	if (impute_maxsteps <= 0 || impute_burnin <= 0 || impute_samples <= 0 || impute_samples > impute_maxsteps)
	{
		WarnError (_String("ERROR: Invalid IMPUTE setting(s) in ImputeNodeScore()"));
		return 0.;
	}
	
	
	// number of steps between samples
	sampling_interval = (long) (impute_maxsteps / impute_samples);
	
	
	// partition parent list into discrete and continuous; record stuff while we're at it
	num_parent_combos = 1;
	family_nlevels.lData[0] = r_i;
	
	for (long nlev, p = 0; p < parents.lLength; p++)
	{
		if (data_type.lData[parents.lData[p]] == 0)	// discrete
		{
			dparents << parents.lData[p];
			parents_by_nodetype << dparents.lLength;
			
			nlev = family_nlevels.lData[p+1] = num_levels.lData[parents.lData[p]];
			
			num_parent_combos *= nlev;
			multipliers << num_parent_combos;
			
			if (nlev > max_num_levels)	max_num_levels = nlev;
		}
		else
		{
			cparents << parents.lData[p];
			parents_by_nodetype << cparents.lLength;
			family_nlevels.lData[p+1] = 0;
		}
	}
	
	
	
	
	// allocate space to matrices
	CreateMatrix (&n_ijk, num_parent_combos, r_i, false, true, false);
	CreateMatrix (&n_ij, num_parent_combos, 1, false, true, false);
	CreateMatrix (&data_deep_copy, theData->GetHDim(), family_size, false, true, false);
	
	pa_indices.Populate(theData->GetHDim(), 0, 0);
	
	
	// allocate space to _GrowingVector object -- used for imputing discrete nodes only
	for (long pa = 0; pa < num_parent_combos; pa++)	reassign_probs->Store (0.);
	
	
	
	// for discrete nodes, the j-th column tallies the number of observed instances of j-th level
	// for continuous nodes, store summary statistics (mean, variance, sample size)
	CreateMatrix (&observed_values, family_size, (max_num_levels > 3) ? max_num_levels : 3, false, true, false);
	
	
	
	// make deep copy, annotate which entries are missing, and tally observed states
	for (long pa_index, row = 0; row < theData->GetHDim(); row++)
	{
		pa_index	= 0;
		child_state	= (*theData) (row, node_id);
		
		data_deep_copy.Store (row, 0, child_state);	// first entry in row always stores child state
		
		// findex = family-wise index
		for (long fnode, fstate, findex = 0; findex < family_size; findex++)
		{
			fnode = (findex == 0) ? node_id : parents.lData[findex-1];	// parent node, get network-wide index
			fstate = (*theData) (row, fnode);
			data_deep_copy.Store (row, findex, fstate);
			
			if (data_type.lData[fnode] == 0)	// discrete node
			{
				if (fstate < 0)
				{
					is_missing << row * family_size + findex;
				}
				else
				{
					observed_values.Store (findex, fstate, observed_values(findex,fstate) + 1);
				}
			}
			else								// continuous node
			{
				if (fstate == continuous_missing_value)	
				{
					is_missing << row * family_size + findex;
				}
				else	
				{
					// column 0 stores number of observations
					observed_values.Store (findex, 0, observed_values(findex,0) + 1);
					
					// column 1 stores sum for computing mean
					observed_values.Store (findex, 1, observed_values(findex,1) + fstate);
					
					// column 2 stores sum of squares for computing variance
					observed_values.Store (findex, 2, observed_values(findex,2) + fstate*fstate);
				}
			}
		}
	}
	
	
	
	// convert observed state vectors into proportions
	for (long obs_total, fnode, findex = 0; findex < family_size; findex++)
	{
		obs_total = 0;
		fnode = (findex == 0) ? node_id : parents.lData[findex-1];
		
		if (data_type.lData[fnode] == 0)	// normalize discrete levels
		{
			for (long lev = 0; lev < family_nlevels.lData[findex]; lev++)
				obs_total += observed_values (findex, lev);	// sum tallies across levels
			
			if (obs_total == 0)	// uh-oh, no observations!  This is a latent variable, assume uniform "prior"
			{
				for (long lev = 0; lev < family_nlevels.lData[findex]; lev++)
					observed_values.Store (findex, lev, 1./family_nlevels.lData[findex]);
			}
			else
			{
				for (long lev = 0; lev < family_nlevels.lData[findex]; lev++)
					observed_values.Store (findex, lev, observed_values(findex, lev) / (_Parameter) obs_total);
			}
		}
		else	// compute mean and standard deviation -- if no observations, use prior parameter settings
		{
			if (observed_values(findex,1) > 0)
			{
				observed_values.Store (findex, 1, observed_values(findex,1) / observed_values(findex,0));
				observed_values.Store (findex, 1, sqrt(observed_values(findex,2) / observed_values(findex,0) 
												- observed_values(findex,1) * observed_values(findex,1)) );
			}
			else
			{
				observed_values.Store (findex, 1, prior_mean(fnode,0));
				observed_values.Store (findex, 2, sqrt(1./prior_precision(fnode,0)) );
			}
		}
	}
	
	
	
	// initialize missing entries to random assignments based on observed cases or prior info
	for (long fnode, row, col, missing_idx = 0; missing_idx < is_missing.lLength; missing_idx++)
	{
		row		= is_missing.lData[missing_idx] / family_size;
		col		= is_missing.lData[missing_idx] % family_size;
		fnode	= (col == 0) ? node_id : parents.lData[col];
		
		if (data_type.lData[fnode] == 0)	// discrete node
		{
			urn	= genrand_real2 ();
			
			for (long level = 0; level < family_nlevels.lData[col]; level++)
			{
				if (urn < observed_values (col, level))
				{
					data_deep_copy.Store (row, col, (_Parameter) level);
					break;
				}
				else
				{
					urn -= observed_values (col, level);
				}
			}
		}
		else								// continuous node
		{
			data_deep_copy.Store (row, col, gaussianDeviate() + observed_values(col,1) * observed_values(col,2));
		}
	}
	
	
	
	// the rest depends on whether child node is discrete or CG
	if (data_type.lData[node_id] == 0)
	{
		// discrete child node --- use Gibbs sampling
		
		// tally N_ijk summary statistics
		for (long pa_index, row = 0; row < data_deep_copy.GetHDim(); row++)
		{
			pa_index	= 0;
			child_state	= data_deep_copy (row, 0);
			
			for (long dpar = 0; dpar < dparents.lLength; dpar++)
			{
				pa_index += (long) data_deep_copy (row, dpar+1) * multipliers[dpar];
			}
			
			n_ijk.Store (pa_index, (long) child_state, n_ijk(pa_index, (long) child_state) + 1);
			n_ij.Store (pa_index, 0, n_ij(pa_index, 0) + 1);
		}
		
		
		// Gibbs sampler
		for (long iter = 0; iter < impute_maxsteps + impute_burnin; iter++)
		{
			is_missing.Permute(1);	// shuffle list of missing entries
			
			log_score = K2Score (r_i, n_ij, n_ijk);	// compute initial score
			
			for (long row, col, pa_index, missing_idx = 0; missing_idx < is_missing.lLength; missing_idx++)
			{
				row			= is_missing.lData[missing_idx] / family_size;
				col			= is_missing.lData[missing_idx] % family_size;
				pa_index	= 0;
				denom		= 0.;
				
				
				// determine parent combination for this row -- use [pa_indices] object instead?  AFYP
				for (long dpar = 0; dpar < dparents.lLength; dpar++)
					pa_index += ((long) data_deep_copy (row, dpar+1)) * multipliers.lData[dpar];
				
				
				reassign_probs->ZeroUsed();	// reset _GrowingVector
				
				
				if (col == 0)	// child node, reassign N_ijk's, keep N_ij's constant
				{				
					child_state = data_deep_copy (row, col);
					log_score -= log(n_ijk(pa_index, child_state));
					
					n_ijk.Store (pa_index, (long) child_state, n_ijk(pa_index, (long) child_state) - 1);
					
					// compute probabilities for all instantiations of child node
					for (long k = 0; k < r_i; k++)
					{
						reassign_probs->Store (log_score + log (n_ijk(pa_index, k) + 1));
					}
					
					denom = LogSumExpo (reassign_probs);
					
					// random reassignment
					urn = genrand_real2 ();
					
					for (long lev = 0; lev < r_i; lev++)
					{
						if ( urn  <  ( this_prob = exp ((* (_Matrix *) reassign_probs) (lev,0) - denom) ))
						{
							n_ijk.Store (pa_index, lev, n_ijk(pa_index,lev) + 1);
							data_deep_copy.Store (row, col, lev);
							
							log_score = (* (_Matrix *) reassign_probs) (lev, 0);
							
							if (missing_idx == is_missing.lLength - 1) vector_of_scores->Store ((* (_Matrix *) reassign_probs) (lev,0));
							
							break;
						}
						
						urn -= this_prob;
					}
					
				}
				else	// parent node, reassign both N_ij's AND N_ijk's
				{
					parent_state	= data_deep_copy (row, col);
					child_state		= data_deep_copy (row, 0);
					
					log_score	-= log (n_ijk (pa_index, (long) child_state));
					log_score	+= log (n_ij  (pa_index, 0) + 1);
					
					n_ij.Store  (pa_index, 0, n_ij(pa_index,0) - 1);
					n_ijk.Store (pa_index, (long) child_state, n_ijk(pa_index, (long) child_state) - 1);
					
					pa_index	-= parent_state * multipliers.lData[col-1];
					
					
					for (long pa_temp, lev = 0; lev < family_nlevels.lData[col]; lev++)
					{
						pa_temp = pa_index + lev * multipliers.lData[col-1];
						reassign_probs->Store (log_score + log (n_ijk (pa_temp, (long) child_state) + 1)
											   - log (n_ij (pa_temp, 0) + 2) );
					}
					
					
					// compute denominator and convert log-values to probability vector
					denom = LogSumExpo (reassign_probs);
					
					// re-assign state
					urn = genrand_real2();	// a uniform random number within interval [0,1)
					
					for (long lev = 0; lev < family_nlevels.lData[col]; lev++)
					{
						this_prob = exp((* (_Matrix *) reassign_probs) (lev, 0) - denom);
						
						if (urn < this_prob)
						{
							pa_index += lev * multipliers.lData[col-1];
							
							n_ij.Store (pa_index, 0, n_ij(pa_index,0) + 1);
							n_ijk.Store (pa_index, (long) child_state, n_ijk(pa_index, (long) child_state) + 1);
							
							log_score = (* (_Matrix *) reassign_probs) (lev, 0);
							
							data_deep_copy.Store (row, col, lev);
							if (missing_idx == is_missing.lLength - 1) vector_of_scores->Store ((* (_Matrix *) reassign_probs) (lev,0));
							
							break;
						}
						
						urn -= this_prob;
					}
					
				}
				
			}
			// end loop over missing values	
			
			// record current log score
			if (iter > impute_burnin && ((iter-impute_burnin) % sampling_interval == 0))
			{
				vector_of_scores->Store (log_score);
			}
		}
	}
	
	else		// child node is CG
	{
		long			k_i			= cparents.lLength,
						pa_index;
		
		_Parameter		lk_ratio, 
						next_log_score;		// log-likelihood of the proposed step
		
		_Matrix			log_scores_by_pa (num_parent_combos, 1, false, true),	// track separately to make updating easier
						next_log_score_by_pa, num_parent_combos, 1, false, true);
		
		
		
		// set mean intercept/regression coeff. hyperparameter vector for CG prior
		CreateMatrix (&mu, k_i+1, 1, false, true, false);
		
		mu.Store (0, 0, prior_mean (node_id, 0));
		for (long i = 1; i < k_i+1; i++) 
			mu.Store (i,0,0);	// set prior expectation of regression coeffs. to 0
		
		
		// set precision hyperparameter for CG prior
		CreateMatrix (&tau, k_i+1, k_i+1, false, true, false);
		
		tau.Store (0, 0, prior_precision (node_id, 0));
		for (long i = 1; i < k_i+1; i++)
			tau.Store (i, i, prior_precision(cparents.lData[i],0));
		
		
		// partition data set by discrete parent state combination
		for (long multiplier, obs = 0; obs < theData->GetHDim(); obs++)
		{
			pa_index	= 0;
			multiplier	= 1;	// takes place of num_parent_combos
			
			for (long this_parent, dpar = 0; dpar < dparents.lLength; dpar++)
			{
				this_parent = dparents.lData[dpar];
				pa_index += data_deep_copy(obs, dpar) * multiplier;
				multiplier *= num_levels.lData[this_parent];
			}
			
			pa_indices << pa_index;
			n_ij.Store (pa_index, 0, n_ij(pa_index,0) + 1);
		}
		
		
		
		// calculate CG summary statistics and compute initial scores
		log_score = 0.;
		
		for (long batch_count, pa = 0; pa < num_parent_combos; pa++)
		{
			// these are single-use matrices, reset with every loop
			CreateMatrix (&zbpa, n_ij(pa,0), k_i+1, false, true, false);	// (1, CG parent node values) in batch
			CreateMatrix (&yb, n_ij(pa,0), 1, false, true, false);			// CG child node values in batch
			
			batch_count = 0;	// number of cases with this parent combo
			
			// populate zbpa with continuous parent states
			for (long obs = 0; obs < pa_indices.lLength; obs++)
			{
				if (pa_indices.lData[obs] == pa)
				{
					// (1, x_0, x_1, ..., x_N) 
					zbpa.Store (batch_count, 0, 1);	// intercept
					
					for (long findex = 1; findex < family_size; findex++)
					{
						if (data_type.lData[parents.lData[findex-1]] == 1)	// CG parent
						{
							zbpa.Store (batch_count, parents_by_nodetype.lData[findex-1], data_deep_copy(obs, findex));
						}
					}
					
					yb.Store (batch_count, 0, data_deep_copy(obs, 0));	// child node state
					
					batch_count++;
				}
			}
			
			log_scores_by_pa.Store (pa, 0, BottcherScore (yb, zbpa, mu, tau, rho, phi, n_ij(pa,0)));
			log_score += log_scores_by_pa (pa, 0);
		}
		
		
		
		// Metropolis-Hastings sampling, but loop through all missing entries at each step
		for (long iter = 0; iter < impute_burnin + impute_maxsteps; iter++)
		{
			is_missing.Permute(1);
			
			for (long batch_count, row, col, missing_idx = 0; missing_idx < is_missing.lLength; missing_idx++)
			{
				// indices into data_deep_copy for this missing value
				row			= is_missing.lData[missing_idx] / family_size;
				col			= is_missing.lData[missing_idx] % family_size;
				denom		= 0.;
				pa_index	= pa_indices.lData[row];
				
				
				// if child node, then row remains in same 'pa' batch, Z^b_pa and N_ij unchanged
				if (col == 0)
				{
					// random draw from Gaussian distribution centered at previous value
					child_state = data_deep_copy (row, col);
					data_deep_copy.Store (row, col, observed_values(0,2) * (gaussianDeviate() + data_deep_copy(row,col)) );
					
					
					// re-calculate summary stats
					CreateMatrix (&zbpa, n_ij(pa_index,0), k_i+1, false, true, false);
					CreateMatrix (&yb, n_ij(pa_index,0), 1, false, true, false);			
					
					batch_count = 0;
					
					for (long obs = 0; obs < pa_indices.lLength; obs++)
					{
						if (pa_indices.lData[obs] == pa_index)
						{
							zbpa.Store (batch_count, 0, 1);
							
							for (long findex = 1; findex < family_size; findex++)
							{
								if (data_type.lData[parents.lData[findex-1]] == 1)
								{
									zbpa.Store (batch_count, parents_by_nodetype.lData[findex-1], data_deep_copy(obs, findex));
								}
							}
							
							yb.Store (batch_count, 0, data_deep_copy(obs, 0));
							
							batch_count++;
						}
					}
					
					next_log_score_by_pa.Store (pa_index, 0, BottcherScore (yb, zbpa, mu, tau, rho, phi, (long)n_ij(pa_index,0)));
					next_log_score	= log_score - log_scores_by_pa (pa_index, 0) + next_log_score_by_pa(pa_index, 0);
					
					lk_ratio		= exp(log_score - next_log_score);
					
					
					// accept this step?
					if (lk_ratio >= 1. || genrand_real2() < lk_ratio)
					{
						log_score = next_log_score;
						
						// update log score contribution for this parent combo (pa)
						log_scores_by_pa.Store (pa_index, 0, next_log_score_by_pa(pa_index, 0));
					}
					else
					{
						// revert to previous state
						data_deep_copy.Store (row, col, child_state);
					}
				}
				
				/* ---------- MISSING VALUE AT PARENT NODE (col > 0) ---------- */
				else
				{
					// if parent is discrete, then this case is moved to another batch (pa) -- MORE THAN ONE BATCH IS AFFECTED
					if (data_type.lData[parents.lData[col-1]] == 0)
					{
						parent_state	= data_deep_copy (row, col);
						
						
						// first compute the likelihood of the original batch minus this case
						n_ij.Store  (pa_index, 0, n_ij(pa_index,0) - 1);
						CreateMatrix (&zbpa, n_ij(pa_index,0), k_i+1, false, true, false);
						CreateMatrix (&yb, n_ij(pa_index,0), 1, false, true, false);
						
						for (long obs = 0; obs < data_deep_copy.GetHDim(); obs++)
						{
							if (obs == row) continue;	// skip this case
							
							if (pa_indices.lData[obs] == pa_index)
							{
								zbpa.Store (batch_count, 0, 1);
								
								for (long findex = 1; findex < family_size; findex++)
								{
									if (data_type.lData[parents.lData[findex-1]] == 1)
									{
										zbpa.Store (batch_count, parents_by_nodetype.lData[findex-1], data_deep_copy(obs, findex));
									}
								}
								
								yb.Store (batch_count, 0, data_deep_copy(obs, 0));
								
								batch_count++;
							}
						}
						
						next_log_score_by_pa.Store (pa_index, 0, BottcherScore (yb, zbpa, mu, tau, rho, phi, (long) n_ij (pa_index, 0)));
						next_log_score = log_score - log_scores_by_pa (pa_index,0) + next_log_score_by_pa (pa_index,0);
						
						
						// random re-assignment of discrete parent node, compute likelihood of second batch
						urn = genrand_real2() - observed_values(col, parent_state);
						
						for (long lev = 0; lev < family_nlevels.lData[col]; lev++)
						{
							if (lev == parent_state) continue;
							
							if (urn < observed_values(col, lev))
							{
								deep_data_copy.Store (row, col, lev);
								pa_index += (lev - parent_state) * multipliers.lData[parents_by_nodetype.lData[col-1]];
								
								
								// compute summary stats
								n_ij.Store  (pa_index, 0, n_ij(pa_index,0) + 1);
								batch_count = 0;
								CreateMatrix (&zbpa, n_ij(pa_index,0), k_i+1, false, true, false);
								CreateMatrix (&yb, n_ij(pa_index,0), 1, false, true, false);
								
								for (long obs = 0; obs < data_deep_copy.GetHDim(); obs++)
								{
									// we won't update pa_indices[] for this case, already have pa_index to work with
									if (obs == row || pa_indices.lData[obs] == pa_index)
									{
										zbpa.Store (batch_count, 0, 1);
										
										for (long findex = 1; findex < family_size; findex++)
										{
											if (data_type.lData[parents.lData[findex-1]] == 1)
											{
												zbpa.Store (batch_count, parents_by_nodetype.lData[findex-1], data_deep_copy(obs, findex));
											}
										}
										yb.Store (batch_count, 0, data_deep_copy(obs, 0));
										batch_count++;
									}
								}
								
								next_log_score_by_pa.Store (pa_index, 0, BottcherScore (yb, zbpa, mu, tau, rho, phi, (long) n_ij (pa_index, 0)));
								break;
							}
							else
							{
								urn -= observed_values(col, lev);
							}
						}
						
						next_log_score -= log_scores_by_pa (pa_index, 0);	// remove old score for second batch
						next_log_score += next_log_score_by_pa (pa_index, 0);				// add new scores for both batches
						
						lk_ratio		= exp(log_score - next_log_score);
						if (lk_ratio >= 1. || genrand_real2() < lk_ratio)
						{
							log_score = next_log_score;
							log_scores_by_pa.Store (pa_index, 0, next_log_score_by_pa (pa_index,0) );
							log_scores_by_pa.Store (pa_indices.lData[row], 0, next_log_score_by_pa (pa_indices.lData[row],0) );	// contains the old index
							
							pa_indices.lData[row] = pa_index;	// replace old pa index with new
						}
						else
						{
							// revert to previous state
							data_deep_copy.Store (row, col, parent_state);
							
							n_ij.Store (pa_index, 0, n_ij(pa_index,0) - 1);
							n_ij.Store (pa_indices.lData[row], 0, n_ij(pa_indices.lData[row],0) + 1);
							pa_index = pa_indices.lData[row];
						}
					}
					
					// otherwise, if parent is continuous then just update Z_bpa; case stays in the same batch
					else
					{
						// random draw from Gaussian distribution centered at previous value
						data_deep_copy.Store (row, col, (gaussianDeviate() + data_deep_copy(row,col)) * observed_values(0,2));
						
						CreateMatrix (&zbpa, n_ij(pa_index,0), k_i+1, false, true, false);
						CreateMatrix (&yb, n_ij(pa_index,0), 1, false, true, false);
						
						for (long obs = 0; obs < data_deep_copy.GetHDim(); obs++)
						{
							if (pa_indices.lData[obs] == pa_index)
							{
								zbpa.Store (batch_count, 0, 1);
								
								for (long findex = 1; findex < family_size; findex++)
								{
									if (data_type.lData[parents.lData[findex-1]] == 1)
									{
										zbpa.Store (batch_count, parents_by_nodetype.lData[findex-1], data_deep_copy(obs, findex));
									}
								}
								
								yb.Store (batch_count, 0, data_deep_copy(obs, 0));
								
								batch_count++;
							}
						}
						
						next_log_score_by_pa.Store (pa_index, 0, BottcherScore (yb,zbpa,mu,tau,rho,phi, (long)n_ij(pa_index,0)));
						next_log_score = log_score - log_scores_by_pa(pa_index,0) + next_log_scores_by_pa(pa_index,0);
						
						lk_ratio = exp(log_score - next_log_score);
						
						if (lk_ratio >= 1. || genrand_real2() < lk_ratio)
						{
							log_score = next_log_score;
							log_scores_by_pa.Store (pa_index, 0, next_log_score_by_pa(pa_index,0));
							pa_indices.lData[row] = pa_index;
						}
						else
						{
							data_deep_copy.Store (row, col, parent_state);
						}
					}
				}
				
				
			}
			// end loop over missing entries
			
			if (iter > impute_burnin && ((iter - impute_burnin) % sampling_interval))
			{
				vector_of_scores.Store (log_score);
			}
		}
		// end sampler
	}
	
	
	// compute the average of sampled log-likelihoods
	log_score = LogSumExpo (vector_of_scores) - log((double)vector_of_scores->GetUsed());
	
	DeleteObject (vector_of_scores);
	DeleteObject (reassign_probs);
	
	return (log_score);
}

#endif






















