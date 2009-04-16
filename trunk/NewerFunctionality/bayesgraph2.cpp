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

extern _Parameter LnGamma (_Parameter);


_AssociativeList *	_BayesianGraphicalModel::ExportModel (void)
{
	return nil;
}



bool _BayesianGraphicalModel::ImportCache (_AssociativeList *)
{
	return 0;
}


//___________________________________________________________________________________________
void _BayesianGraphicalModel::ExportCache (_AssociativeList * cache_export)
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
	}
	else
	{
		WarnError (_String ("Unable to export node score cache, no cache exists!"));
	}
}




//___________________________________________________________________________________________
#define _DEBUG_CCS_
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
	/*
		Compute the conditional Gaussian network score of node [c] according to S. Bottcher's formulation.
		
		Let [y] denote the set of other continuous nodes.
		Let [i] denote the set of discrete nodes.
		Let [pa] define a subset of continuous or discrete nodes on [c]
		
		Assume local probability distribution is Gaussian local regression on continuous parents, with 
		parameters conditional on discrete parent configuration [pa_i].  This has parameters:
			beta (c | pa_i(c)) = vector of regression coefficients over pa_y (continuous parents)
			m (c | pa_i(c)) = regression intercept
			sigma^2 (c | pa_i(c)) = variance conditional on discrete parent configuration
	 
		Hence, for a case {y_1, y_2, ..., y_c-1, y_c+1, ..., y_m, i_1, i_2, ..., i_n},
			(y_c | pa_i(c) ~ Gaussian(m + beta_1 y_1 + ... + beta_m y_m, sigma^2)
		
		We need conjugate priors for these parameters, which are given by:
			(m, beta | sigma^2) ~ Gaussian (mu, sigma^2 tau^-1)
			(sigma^2) ~ Inverse-Gamma ( rho_{y|pa_i(y)} / 2, phi_{y|pa_i(y)} / 2)
		
		In other words, the intercept and regression coefficients of continuous node (c) on continuous
		parents pa_y(c) given discrete parent configuration pa_i(y) 
	 
		We integrate over these priors given the hyperparameters:
			mu		- prior mean vector (intercept, regression coefficients...)
			tau		- prior precision
			rho		- degrees of freedom for inverse gamma conjugate prior
			phi		- scale parameter for inverse gamma conjugate prior
	 
		Typically, [rho] and [phi] are given small values (0.001) that are meant to render an
		uninformative prior (flatter).  However, Gelman asserts that this practice is not optimal.
	 
		[tau] is a precision matrix whose off-diagonal entries are usually zeroed out.  Hence, the
		values set by an AVL upon constructing the BGM fill the diagonal.  (I should probably implement
		some access functions so that these can be changed.)
	*/
	
	ReportWarning (_String ("ComputeContinuousScore ") & node_id);
	
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


	_SimpleList		n_ij,
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
	for (long pa = 0; pa < num_parent_combos; pa++)
	{
#ifdef _DEBUG_CCS_
		sprintf (debug, "parent combo %d of %d\n", pa, num_parent_combos);
		BufferToConsole (debug);
#endif
		
		_Parameter	pa_log_score = 0.;
		
		_Matrix		zbpa (n_ij.lData[pa], k+1, false, true),
					yb (n_ij.lData[pa], 1, false, true),
					scale;
		
		long		count_n		= 0;	
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
		
		
		// calculate scale parameter for non-central t distribution
		// from S. Bottcher (2001) p.25
		
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
		
		pa_log_score += LnGamma((rho + n_ij.lData[pa])/2.);
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
		
		pa_log_score += -(rho + n_ij.lData[pa])/2. * (next_mat(0,0) + 1.);
		
#ifdef _DEBUG_CCS_
		sprintf (debug, "pa_log_score: %1.3f\n", pa_log_score);
		BufferToConsole (debug);
		
		log_score += pa_log_score;
#endif
	}
	
	return log_score;
}



_Parameter _BayesianGraphicalModel::GibbsImputeScore (long node_id, _SimpleList & parents)
{
	return 0;
}

#endif
