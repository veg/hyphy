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
	ReportWarning (_String ("ComputeContinuousScore ") & node_id);
	
	/* WARNING, untested function! */
	/* --------------------------- */
	
	// mm = prior estimate of unconditional mean at continuous node (i.e. intercept)
	// phi = scale parameter of inverse gamma prior for variance, uninformative at low values 
	//						(Kohn, Smith, and Chan, 2001 Stat Comput 11: 313-322)
	
	_Parameter		log_score = 0.;
	
	long			num_parent_combos = 1,	// i.e. 'q'
					k;						// current number of continuous parents
	_Matrix			n_ij,
					pa_indexing,		// track discrete parent combinations per observation
					mu,
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
	CreateMatrix (&mu, k+1, 1, false, true, false);
	mu.Store (0, 0, prior_mean (node_id, 0));				// prior intercept
	for (long i = 1; i < mu.GetHDim(); i++)
	{
		mu.Store (i, 0, 0);		// set prior expectation of regression coefficients to zero
	}
	
	
	// set prior degrees of freedom (for inverse gamma / scaled inverse chi-square)
	_Parameter		rho	= prior_sample_size (node_id, 0) > 0 ? (prior_sample_size (node_id, 0) / num_parent_combos) : 1.0;
	
	
	
	// set precision hyperparameter for Gaussian prior
	CreateMatrix (&tau, k+1, k+1, false, true, false);
	for (long row = 0; row < k+1; row++)
	{
		for (long col = 0; col < k+1; col++)
		{
			if (row == col) tau.Store (row, col, rho);	// set diagonal entries to rho
			else tau.Store (row, col, 0.);	// zero off-diagonal entries
		}
	}
	
	
	
	// count up number of data points per parent combination
	CreateMatrix (&n_ij, num_parent_combos, 1, false, true, false);
	CreateMatrix (&pa_indexing, theData->GetHDim(), 1, false, true, false);
	
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
		
		pa_indexing.Store (obs, 0, index);
		n_ij.Store (index, 0, n_ij(index, 0) + 1);
	}
	
	
	
	// for every parent combination, calculate contribution to score
	for (long pa = 0; pa < num_parent_combos; pa++)
	{
		_Parameter	pa_log_score = 0.;
		
		_Matrix		zbpa (n_ij(pa, 0), k+1, false, true),
					yb (n_ij(pa, 0), 1, false, true),
					scale;
		
		long		count_n		= 0;	
		// number of data points with this parent combination.
		
		
		// populate zbpa matrix with (1, y_1, y_2, ..., y_m) entries for this parental combo
		for (long obs = 0; obs < theData->GetHDim(); obs++)
		{
			if (pa_indexing(obs, 0) == pa)		// this observation has the current parent combo
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
		
		
		ReportWarning (_String ("before matrix inverse"));
		
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
		scale *= (_Parameter) (prior_precision (node_id, 0) / rho);
		
		
		ReportWarning (_String ("before matrix determinant"));
		
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
		pa_log_score += LnGamma((rho + n_ij(pa, 0))/2.);
		pa_log_score -= LnGamma(rho/2.) + 0.5 * log(det);
		
		
		// calculate second term of score
		_Matrix		next_mat;
		ReportWarning (_String ("before zbpa * mu"));
		zbpa *= mu;		// n_ij x (k+1)  *  (k+1) x 1  --> n_ij x 1
		ReportWarning (_String ("before yb - zbpa"));
		yb -= zbpa;		// n_ij x 1
		temp_mat = yb;
		next_mat = temp_mat;
		ReportWarning (_String ("before yb * scale-1"));
		next_mat *= * (_Matrix *) scale.Inverse();	// incompatible matrix dimensions  n_ij x 1  * n_ij x n_ij
		ReportWarning (_String ("scale dimensions: ") & scale.GetHDim() & " x " & scale.GetVDim());
		ReportWarning (_String ("next_mat dimensions: ") & next_mat.GetHDim() & " x " & next_mat.GetVDim());
		temp_mat.Transpose();
		ReportWarning (_String ("before yb * scale-1 * yb^T"));
		next_mat *= temp_mat;
		next_mat *= (_Parameter) (1./prior_precision (node_id, 0));	// should be a 1-element matrix
		
		pa_log_score += -(rho + n_ij(pa,0))/2. * (next_mat(0,0) + 1.);
		log_score += pa_log_score;
	}
	
	ReportWarning (_String ("returning from ComputeContinuousScore()"));
	
	return log_score;
}



_Parameter _BayesianGraphicalModel::GibbsImputeScore (long node_id, _SimpleList & parents)
{
	return 0;
}

#endif
