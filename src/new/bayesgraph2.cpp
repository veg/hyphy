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


#include "bayesgraph.h"

extern _String      _HYBgm_IMPUTE_MAXSTEPS,
       _HYBgm_IMPUTE_BURNIN,
       _HYBgm_IMPUTE_SAMPLES;

extern _Parameter   lnGamma (_Parameter),
       gaussDeviate (void),
       LogSumExpo (_GrowingVector *);




//___________________________________________________________________________________________
void _BayesianGraphicalModel::SerializeBGM (_String & rec)
{
    /* --------------------------------------------------------------------
        ExportModel()
            Export the network structure and parameters as an associative
            array.
       -------------------------------------------------------------------- */

    ReportWarning (_String("Entered _BayesianGraphicalModel::ExportModel()"));

    if (theData.GetHDim() > 0) {
        _String *   bgmName = (_String *) bgmNamesList (bgmList._SimpleList::Find((long)this));
        _String     bgmNameLocal;
        _String     nodeKey;

        bgmNameLocal.Duplicate(bgmName);
        bgmNameLocal << "_export";

        rec << bgmNameLocal;
        rec << "={};\n";

        for (long node = 0; node < num_nodes; node++) {
            nodeKey = bgmNameLocal & "[\"" & node & "\"]";

            rec << nodeKey;
            rec << " = {};\n";

            // record node type
            rec << '(';
            rec << nodeKey;
            rec << ")[\"NodeType\"] = ";
            rec << _String(node_type.lData[node]);
            rec << ";\n";

            // record parents as _Matrix object
            _SimpleList parents,
                        d_parents, c_parents;

            for (long pnode = 0; pnode < num_nodes; pnode++) {
                if (pnode == node) {
                    continue;
                }
                if (theStructure(pnode,node) == 1) {
                    parents << pnode;

                    if (node_type.lData[pnode] == 0) {
                        d_parents << pnode;
                    } else {
                        c_parents << pnode;
                    }
                }
            }

            rec << '(';
            rec << nodeKey;
            rec << ")[\"Parents\"] = {";
            rec << _String((_String *)parents.toStr());
            rec << "};\n";

            // record network parameters as conditional probability density functions
            rec << '(';
            rec << nodeKey;
            rec << ")[\"CPDFs\"] = {{";


            // posterior probability density functions (PDFs) conditional on child node state (indexed by i)
            if (node_type.lData[node] == 0) {
                // discrete child node, compute Dirichlet hyperparameters
                _Matrix     n_ij, n_ijk;


                // tally cases by discrete parent combination and child state
                UpdateDirichletHyperparameters (node, d_parents, &n_ij, &n_ijk);


                // report posterior Dirichlet PDFs
                for (long j = 0; j < n_ij.GetHDim(); j++) {
                    rec << "\"Random({{";
                    for (long k = 0; k < num_levels.lData[node]; k++) { // k indexes child node states
                        rec << _String(n_ijk(j,k));
                        if (k < num_levels.lData[node]-1) {
                            rec << ',';
                        }
                    }
                    rec << "}},{\\\"PDF\\\":\\\"Dirichlet\\\"})\"";
                    if (j < n_ij.GetHDim()-1) {
                        rec << ',';
                    }
                }
            }

            else if (node_type.lData[node] == 1) {  // CG child node
                long            k                   = c_parents.lLength,        // a useful short-hand
                                num_parent_combos = 1;

                _SimpleList     multipliers ((long)1),
                                n_ij, pa_indexing;

                _Parameter      rho_prior   = prior_sample_size (node, 0) > 0 ? (prior_sample_size (node, 0) / num_parent_combos) : 1.0,
                                phi_prior = prior_scale (node, 0);

                _Matrix         mu_prior (k+1, 1, false, true),
                                tau_prior (k+1, k+1, false, true);


                // index cases by discrete parent configurations
                for (long par = 0; par < d_parents.lLength; par++) {
                    num_parent_combos *= num_levels.lData[d_parents.lData[par]];
                    multipliers << num_parent_combos;
                }

                n_ij.Populate (num_parent_combos, 0, 0);
                pa_indexing.Populate (theData.GetHDim(), 0, 0);

                if (d_parents.lLength > 0) {
                    for (long obs = 0; obs < theData.GetHDim(); obs++) {
                        long    index       = 0,
                                multiplier = 1;

                        for (long par = 0; par < d_parents.lLength; par++) {
                            long    this_parent     = parents.lData[par];

                            index += theData(obs, this_parent) * multiplier;    // max index = sum (parent * levels(parent))
                            multiplier *= num_levels.lData[this_parent];
                        }

                        pa_indexing.lData[obs] = index;
                        n_ij.lData[index] += 1;
                    }
                } else {
                    n_ij.lData[0] = theData.GetHDim();
                }


                // set CG hyperparameter priors (independent of discrete parent config)
                for (long row = 0; row < k+1; row++) {
                    for (long col = 0; col < k+1; col++) {
                        if (row == col) {   // set diagonal entries
                            if (row == 0) {
                                tau_prior.Store (0, 0, prior_precision (node, 0));
                            } else {
                                tau_prior.Store (row, col, prior_precision (c_parents.lData[row],0));
                            }
                        } else {
                            tau_prior.Store (row, col, 0.);    // zero off-diagonal entries
                        }
                    }
                }


                mu_prior.Store (0, 0, prior_mean (node, 0));                // prior intercept
                for (long i = 1; i < k+1; i++) {
                    mu_prior.Store (i, 0, 0);    // set prior expectation of regression coefficients to zero
                }


                // for every discrete parent config (on which CG marginal distribution is conditioned)
                for (long pa = 0; pa < num_parent_combos; pa++) {
                    // compute summary statistics
                    _Matrix     zbpa (n_ij.lData[pa], k+1, false, true),
                                yb (n_ij.lData[pa], 1, false, true);

                    for (long count_n = 0, obs = 0; obs < theData.GetHDim(); obs++) {
                        if (pa_indexing.lData[obs] == pa) {
                            zbpa.Store (count_n, 0, 1);
                            for (long cpar = 0; cpar < k; cpar++) {
                                zbpa.Store (count_n, cpar+1, theData(obs, c_parents.lData[cpar]));
                            }

                            yb.Store (count_n, 0, theData(obs, node));

                            count_n++;
                        }
                    }


                    // update CG hyperparameters with summary statistics (y_b, Z_bpa, n_ij)
                    _Matrix     t_zbpa = zbpa;
                    t_zbpa.Transpose();
                    t_zbpa *= zbpa;     // Z^T x Z

                    _Matrix     tau (tau_prior);    // duplicator
                    tau += t_zbpa;

                    _Matrix     temp (tau_prior);
                    temp *= mu_prior;
                    _Matrix     mu (temp);
                    temp = zbpa;
                    temp.Transpose();
                    temp *= yb;
                    mu += temp;

                    _Matrix     * tauinv = (_Matrix *) tau.Inverse();
                    temp = *tauinv;
                    temp *= mu;
                    mu = temp;

					DeleteObject (tauinv);
					
                    _Parameter  rho = rho_prior + n_ij.lData[pa],
                                phi = phi_prior;

                    temp = zbpa;
                    temp *= mu;
                    temp *= -1;
                    temp += yb;
                    temp.Transpose();

                    temp *= yb;
                    phi += temp(0,0);

                    temp = mu_prior;
                    temp -= mu;
                    temp.Transpose();
                    temp *= tau_prior;
                    temp *= mu_prior;

                    phi += temp(0,0);

                    // report CG posterior distributions to output string
                    rec << "\"Random(";
                    rec << _String((_String *)mu.toStr());
                    rec << ",{\\\"PDF\\\":\\\"Gaussian\\\",\\\"ARG0\\\":Random({{";
                    rec << _String(phi);
                    rec << "}},{\\\"PDF\\\":\\\"InverseWishart\\\",\\\"ARG0\\\":{{";
                    rec << _String(rho);
                    rec << "}}})*";
                    rec << _String((_String*)tau.Inverse()->toStr());
                    rec << "})\"";

                    if (pa < num_parent_combos-1) {
                        rec << ',';
                    }
                }

            }

            rec << "}};\n"; // end matrix entry of conditional PDFs.
        }
    } else {
        WarnError (_String("Cannot export network parameters, this _BayesianGraphicalModel object has no data!"));
    }
}




//___________________________________________________________________________________________
bool _BayesianGraphicalModel::ImportCache (_AssociativeList * cache_import)
{
    //  Unpack contents of associative list sent from HBL.
    //   use to restore node score cache as an alternative to computing scores
    //   de novo from data -- particularly useful for time-consuming imputation
    //   in cases of missing data
    ReportWarning (_String("Entered ImportCache() with avl: ") & (_String *) cache_import->toStr());

    _String         keyString;
    _PMathObj       valuePtr;

    if (scores_cached) {
        ReportWarning (_String("WARNING: Overwriting pre-existing node score cache in bayesgraph2.cpp:ImportCache()"));
    }

    for (long node = 0; node < num_nodes; node++) {
        _String     errMsg;
        _List   *   node_scores = (_List *) node_score_cache.lData[node];

        node_scores->Clear();

        for (long num_parents = 0; num_parents <= max_parents.lData[node]; num_parents++) {
            keyString = _String("Node") & node & "NumParents" & num_parents;

            if (num_parents == 0) {
                if ((valuePtr = cache_import->GetByKey(keyString, NUMBER))) {
                    (*node_scores) && (_Constant *) valuePtr;    // append duplicate
                } else {
                    errMsg = _String ("Expecting numerical value in associative list for key ") & keyString;
                }
            } else if (num_parents == 1) {
                if ((valuePtr = cache_import->GetByKey(keyString, MATRIX))) {
                    (*node_scores) && (_Matrix *) valuePtr;
                } else {
                    errMsg = _String ("Expecting matrix in associative list for key ") & keyString;
                }
            } else {
                if ((valuePtr = cache_import->GetByKey(keyString, MATRIX))) {
                    (*node_scores) && (_NTupleStorage *) valuePtr;
                } else {
                    errMsg = _String("Expecting matrix (_NTupleStorage) object in associative list for key ") & keyString;
                }
            }
        }

        if (errMsg) {
            WarnError (errMsg);
            return false;
        }
    }

    scores_cached = TRUE;
    return true;
}


//___________________________________________________________________________________________
bool _BayesianGraphicalModel::ExportCache (_AssociativeList * cache_export)
{
    // Export an associative list containing node_score_cache contents to HBL
    //  entries are indexed by key string referring to child node and number of parents,
    //      e.g. "N0P2" refers to scores associated with child node 0 and two parents
    _String         keyString;

    _List           * this_list;
    _Constant       * orphan_score;
    _Matrix         * single_parent_scores;
    _NTupleStorage  * family_scores;

    if (scores_cached) {
        ReportWarning (_String("Exporting cache with ") & num_nodes & " nodes");

        for (long node = 0; node < num_nodes; node++) {
            this_list = (_List *) node_score_cache.lData[node];

            for (long npar = 0; npar <= max_parents.lData[node]; npar++) {
                keyString = _String ("Node") & node & "NumParents" & npar;
                _FString    aKey (keyString, false);

                ReportWarning (_String("Inserting with key ") & keyString);

                if (npar == 0) {
                    orphan_score = (_Constant *) this_list->lData[npar];
                    cache_export->MStore (&aKey, orphan_score, true);   // make a dynamic copy
                } else if (npar == 1) {
                    single_parent_scores = (_Matrix *) this_list->lData[npar];
                    cache_export->MStore (&aKey, single_parent_scores, true);   // make a dynamic copy
                } else {
                    family_scores = (_NTupleStorage *) this_list->lData[npar];
                    cache_export->MStore (&aKey, family_scores, true);  // make a dynamic copy
                }
            }
        }

        return TRUE;
    } else {
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
    for (long par = 0; par < num_nodes; par++) {
        if (theStructure(par, node_id) == 1) {
            parents << par;
        }
    }

    return ComputeContinuousScore (node_id, parents);
}


_Parameter _BayesianGraphicalModel::ComputeContinuousScore (long node_id, _Matrix &dag)
{
    _SimpleList     parents;

    for (long par = 0; par < num_nodes; par++) {
        // support for discrete parents only
        if (dag(par, node_id) == 1) {
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
            Returns the network score based on Bottcher's formulation.
            Note, does not update posterior distributions of hyperparameters,
            this requires a separate albeit highly similar function.
       ----------------------------------------------------------------------------- */

    //ReportWarning (_String ("Called ComputeContinuousScore with ") & node_id & " <- " & (_String *) parents.toStr());


    long            k;
    _Parameter      log_score = 0.;
    _SimpleList     c_parents, d_parents;


    // use cached node scores if available
    if (scores_cached) {
        _List *     scores  = (_List *) node_score_cache.lData[node_id];

        if (parents.lLength == 0) {
            _Constant *     orphan_score = (_Constant *) scores->lData[0];
            return (_Parameter) orphan_score->Value();
        } else if (parents.lLength == 1) {
            _Matrix *   single_parent_scores = (_Matrix *) scores->lData[1];
            return (_Parameter) (*single_parent_scores) (parents.lData[0], 0);
        } else {
            _NTupleStorage *    family_scores   = (_NTupleStorage *) scores->lData[parents.lLength];
            _SimpleList         nktuple;

            for (long i = 0; i < parents.lLength; i++) {    // map parents back to nk-tuple
                long    par = parents.lData[i];
                if (par > node_id) {
                    par--;
                }
                nktuple << par;
            }
            return (_Parameter) family_scores->Retrieve (nktuple);  // using nk-tuple
        }
    }


    if (theData.GetHDim() == 0) {
        WarnError (_String ("Uh-oh, there's no node score cache nor is there any data matrix to compute scores from!"));
        return 0.;
    }

	
	// impute score if missing data
	if (has_missing.lData[node_id]) {
		return ImputeCGNodeScore (node_id, parents);
	} else {
		for (long par = 0; par < parents.lLength; par++) {
			if (has_missing.lData[parents.lData[par]]) {
				return ImputeCGNodeScore (node_id, parents);
			}
		}
	}
	

    // partition parent list into discrete and continuous
    for (long p = 0; p < parents.lLength; p++) {
        if (node_type.lData[parents.lData[p]] == 0) {   // discrete
            d_parents << parents.lData[p];	// note these are network-level node indices
        } else {
            c_parents << parents.lData[p];
        }
    }
    k = c_parents.lLength;


    // if there are discrete parents then integrate over marginal posteriors conditioned on discrete parent states

	long            num_parent_combos = 1;
	_SimpleList     multipliers ((long)1);


	// how many combinations of parental states are there?
	for (long par = 0; par < d_parents.lLength; par++) {
		num_parent_combos *= num_levels.lData[d_parents.lData[par]];
		multipliers << num_parent_combos;
	}



	// count up number of data points per discrete parent combination
	_SimpleList     n_ij,               // CAVEAT: other functions use _Matrix for this object
					pa_indexing;        // track discrete parent combinations per observation

	n_ij.Populate (num_parent_combos, 0, 0);
	pa_indexing.Populate (theData.GetHDim(), 0, 0);

	for (long obs = 0; obs < theData.GetHDim(); obs++) {
		long    index       = 0,
				multiplier = 1;

		for (long par = 0; par < d_parents.lLength; par++) {
			long    this_parent     = d_parents.lData[par];

			index += theData(obs, this_parent) * multiplier;    // max index = sum (parent * levels(parent))
			multiplier *= num_levels.lData[this_parent];
		}

		pa_indexing.lData[obs] = index;
		n_ij.lData[index] += 1;
	}


	// set hyperparameter priors -- they will not be modified on pass to BottcherScore() so keep outside loop
	_Matrix     tau (k+1, k+1, false, true),
				mu (k+1, 1, false, true);

	_Parameter  rho     = prior_sample_size (node_id, 0) > 0 ? (prior_sample_size (node_id, 0) / num_parent_combos) : 1.0,
				phi      = prior_scale (node_id, 0);

	for (long row = 0; row < k+1; row++) {
		for (long col = 0; col < k+1; col++) {
			if (row == col) {   // set diagonal entries
				if (row == 0) {
					tau.Store (0, 0, prior_precision (node_id, 0));
				} else {
					tau.Store (row, col, prior_precision (c_parents.lData[row-1],0));
				}
			} else {
				tau.Store (row, col, 0.);    // zero off-diagonal entries
			}
		}
	}

	mu.Store (0, 0, prior_mean (node_id, 0));               // prior intercept
	for (long i = 1; i < mu.GetHDim(); i++) {
		mu.Store (i, 0, 0);    // set prior expectation of regression coefficients to zero
	}


	// for every discrete parent combination (on which Gaussian node state is conditioned)
	for (long pa = 0; pa < num_parent_combos; pa++) {
		// collect cases into a batch defined by discrete parent states
		_Matrix     zbpa (n_ij.lData[pa], k+1, false, true),
					yb (n_ij.lData[pa], 1, false, true);

		for (long count_n = 0, obs = 0; obs < theData.GetHDim(); obs++) {
			if (pa_indexing.lData[obs] == pa) {
				// populate zbpa matrix with (1, y_1, y_2, ..., y_m) entries for this parental combo
				zbpa.Store (count_n, 0, 1);
				for (long cpar = 0; cpar < k; cpar++) {
					zbpa.Store (count_n, cpar+1, theData(obs, c_parents.lData[cpar]));
				}

				yb.Store (count_n, 0, theData(obs, node_id));

				count_n++;
			}
		}
		/*
		ReportWarning (_String("Calling BottcherScore() for node ") & node_id & ", parents " & (_String *)parents.toStr() & ", pa=" & pa);
		ReportWarning (_String("yb=") & (_String *) yb.toStr());
		ReportWarning (_String("zpba=") & (_String *) zbpa.toStr());
		ReportWarning (_String("tau=") & (_String *) tau.toStr());
		ReportWarning (_String("mu=") & (_String *) mu.toStr());
		ReportWarning (_String("rho=") & rho);
			*/
		log_score += BottcherScore (yb, zbpa, tau, mu, rho, phi, n_ij.lData[pa]);
	}


    return log_score;
}



//___________________________________________________________________________________________________

_Parameter  _BayesianGraphicalModel::BottcherScore (_Matrix & yb, _Matrix & zbpa, _Matrix & tau, _Matrix & mu, _Parameter rho, _Parameter phi, long batch_size)
{
    /* -------------------------------------------------------------------------------
     Compute the conditional Gaussian network score of node [c] according to
     Susanne G. Bottcher's formulation.

     Let [y] denote the set of other continuous nodes.
     Let [i] denote the set of discrete nodes.
     Let [pa] define a subset of continuous or discrete nodes on [c]

     Assume local probability distribution is Gaussian local regression on
     continuous parents, with parameters conditional on discrete parent
     configuration [pa_i].  This has parameters:

     beta (c | pa_i(c))     = vector of regression coefficients over
                                pa_y (continuous parents)
     m (c | pa_i(c))        = regression intercept
     sigma^2 (c | pa_i(c))  = variance conditional on discrete parent
                                configuration

     Hence, for a case {y_1, y_2, ..., y_c-1, y_c+1, ..., y_m, i_1, i_2, ..., i_n},
     (y_c | pa_i(c) ~ Gaussian(m + beta_1 y_1 + ... + beta_m y_m, sigma^2)

     but conditional variance (sigma^2) is unknown.

     Posterior distributions are:
        (m, beta | sigma^2) ~ Gaussian (mu, sigma^2 tau^-1)
        (sigma^2) ~ Inverse-Gamma ( rho_{y|pa_i(y)} / 2, phi_{y|pa_i(y)} / 2)

     given the hyperparameters:
     mu     - mean vector (intercept, regression coefficients...) [(k+1) x 1 matrix]
     tau    - precision [(k+1) x (k+1) matrix]
     rho    - degrees of freedom for inverse Gamma conjugate prior [scalar]
     phi    - scale parameter for inverse Gamma conjugate prior [scalar]
    --------------------------------------------------------------------------------- */

    // ReportWarning (_String("Entered BottcherScore()"));

    // calculate scale parameter for non-central t distribution
    // from S. Bottcher (2001) p.25


    // zbpa is a {N x (k+1)} matrix where 'k' is the number of continuous parents
    //  so if there are no parents then it is a column vector of length N (number of cases in data)
    _Matrix temp_mat (zbpa);

    //ReportWarning (_String("zbpa = ") & (_String *) temp_mat.toStr());

	if (tau.GetHDim() == 1 && tau.GetVDim() == 1) {
		temp_mat *= 1/(tau(0,0)); // just use scalar multiplication of the reciprocal
	} else {
		_Matrix * tauinv = (_Matrix *) tau.Inverse();
		temp_mat *= *tauinv;
		DeleteObject (tauinv);
	}

    //ReportWarning (_String("... %*% tau^(-1) = ") & (_String *) temp_mat.toStr());

    _Matrix t_zbpa (zbpa);
    t_zbpa.Transpose();

    temp_mat *= t_zbpa;

    //ReportWarning (_String("... %*% zbpa^T = ") & (_String *) temp_mat.toStr());

    for (long row = 0; row < temp_mat.GetHDim(); row++) {
        temp_mat.Store (row, row, temp_mat(row,row)+(_Parameter)1.);    // add identity matrix
    }

    //ReportWarning (_String("... + I = ") & (_String *) temp_mat.toStr());

    _Matrix scale (temp_mat);
    scale *= phi / rho;

    //ReportWarning (_String("scale = ") & (_String *) scale.toStr());

    // calculate the determinant of scale parameter matrix
    _Parameter      pi_const = 3.141592653589793;

    temp_mat = scale;
    temp_mat *= (_Parameter) (pi_const * rho);
	
	//ReportWarning (_String("BottcherScore() calling Eigensystem on matrix ") & (_String *) temp_mat.toStr() );
	
    _AssociativeList *  eigen       = (_AssociativeList *) temp_mat.Eigensystem();
	
	// sometimes the eigendecomposition fails
	if ( (eigen->GetKeys())->lLength == 0 ) {
		WarnError (_String("Eigendecomposition failed in bayesgraph2.cpp BottcherScore()."));
		return -A_LARGE_NUMBER;
	}
	
    _Matrix *           eigenvalues = (_Matrix *)eigen->GetByKey(0, MATRIX);
    _Parameter          log_det         = 0.;   // compute determinant on log scale to avoid overflow

    // determinant is product of eigenvalues (should be > 0 for positive definite matrices)
    for (long i = 0; i < eigenvalues->GetHDim(); i++) {
        log_det += log((*eigenvalues) (i,0));
    }
	
	DeleteObject (eigen);
	
    //ReportWarning (_String("log(det) = ") & log_det);

    // calculate first term of score
    _Parameter  pa_log_score = 0.;

    pa_log_score += lnGamma((rho + batch_size)/2.);
    pa_log_score -= lnGamma(rho/2.) + 0.5 * log_det;


    // calculate second term of score
    _Matrix     next_mat (yb);  // N x 1

    temp_mat = zbpa;        // N x (k+1)
    temp_mat *= mu;         // (k+1) x 1
    next_mat -= temp_mat;   // (yb - zbpa %*% mu)
    next_mat.Transpose();

    temp_mat = next_mat;
	
	_Matrix * scaleinv = (_Matrix *) scale.Inverse();
	
    temp_mat *= *scaleinv;  // N x N
	
	DeleteObject (scaleinv);
	
    next_mat.Transpose();
    temp_mat *= next_mat;
    temp_mat *= (_Parameter) 1./rho;

    //ReportWarning (_String ("2nd term = ") & (_Parameter) temp_mat(0,0));

    pa_log_score += -(rho + batch_size)/2. * log(1. + temp_mat(0,0));


    //ReportWarning (_String("BottcherScore() returning with log L = ") & pa_log_score);

    return pa_log_score;
}



//___________________________________________________________________________________________________
_Parameter _BayesianGraphicalModel::ImputeDiscreteNodeScore (long node_id, _SimpleList & parents)
{
    /* ---------------------------------------------------------------------------------------
	 ImputeDiscreteNodeScore
		Impute (i.e., "fill in") missing values in the data.  More exactly, take the
	 expectation of node score over imputations of the data sampled using Gibbs
	 sampling.  Missing values are flagged for discrete data by negative values.
	 --------------------------------------------------------------------------------------- */
	
	long            num_parent_combos   = 1,
					r_i                 = num_levels.lData[node_id],	// number of levels for child node
					max_num_levels		= r_i,
					family_size			= parents.lLength + 1,
					sampling_interval,
					parent_state,
					child_state;
	
	_SimpleList     multipliers (1, 1, 0),  // length constructor, populator
					family_nlevels (family_size, 0, 0),
					is_missing,		// store linear indices to missing entries in deep copy
					pa_indices;
	
	
	_Matrix         n_ijk, n_ij,	// summary statistics for discrete node
					data_deep_copy,
					observed_values;
	
	
    _GrowingVector  * vector_of_scores  = new _GrowingVector(),     // store scores sampled during imputation
					* reassign_probs    = new _GrowingVector();
    
	
	_Parameter      log_score           = 0,
					denom, this_prob,
					impute_maxsteps, impute_burnin, impute_samples; // HBL settings
    
	double			urn;	// uniform random number
	
	
	
	// set Gibbs sampler parameters from batch language definitions
    checkParameter (_HYBgm_IMPUTE_MAXSTEPS, impute_maxsteps, 0);
    checkParameter (_HYBgm_IMPUTE_BURNIN, impute_burnin, 0);
    checkParameter (_HYBgm_IMPUTE_SAMPLES, impute_samples, 0);
	
    if (impute_maxsteps <= 0 || impute_burnin < 0 || impute_samples <= 0 || impute_samples > impute_maxsteps) {
        WarnError (_String("ERROR: Invalid IMPUTE setting(s) in ImputeNodeScore()"));
        return 0.;
    }

	// number of steps between samples
    sampling_interval = (long) (impute_maxsteps / impute_samples);
	
	
	family_nlevels.lData[0] = r_i; // for more convenient look-up
	
    for (long nlev, p = 0; p < parents.lLength; p++) {
		nlev = family_nlevels.lData[p+1] = num_levels.lData[parents.lData[p]];
		num_parent_combos *= nlev;
		multipliers << num_parent_combos;
		if (nlev > max_num_levels)  {
			max_num_levels = nlev;
		}
	}
	
	
	
    // allocate space to matrices
    CreateMatrix (&n_ijk, num_parent_combos, r_i, false, true, false);
    CreateMatrix (&n_ij, num_parent_combos, 1, false, true, false);
    CreateMatrix (&data_deep_copy, theData.GetHDim(), family_size, false, true, false);
	
    pa_indices.Populate(theData.GetHDim(), 0, 0);
	
	
    // allocate space to _GrowingVector object -- used for imputing discrete nodes only
    for (long pa = 0; pa < num_parent_combos; pa++) {
        reassign_probs->Store (0.);
    }
	
	
	
    // for discrete nodes, the j-th column tallies the number of observed instances of j-th level
    CreateMatrix (&observed_values, family_size, (max_num_levels > 3) ? max_num_levels : 3, false, true, false);
	
	
    // make deep copy, annotate which entries are missing, and record empirical distributions
    for (long row = 0; row < theData.GetHDim(); row++) {
        // findex = family-wise index
        for (long fnode, fstate, findex = 0; findex < family_size; findex++) {
            fnode = (findex == 0) ? node_id : parents.lData[findex-1];  // parent node, get network-wide index
            fstate = theData (row, fnode);
            data_deep_copy.Store (row, findex, fstate);
			
			if (fstate < 0) {
				// store location of missing entry by linear index
				is_missing << row * family_size + findex;
			} else {
				// use state to index into vector and increment counter
				observed_values.Store (findex, fstate, observed_values(findex,fstate) + 1);
			}
        }
    }
	
	
	
    // make a decision whether to use empirical or prior distribution to do initial imputations
    for (long obs_total, fnode, findex = 0; findex < family_size; findex++) {
        obs_total = 0;
        fnode = (findex == 0) ? node_id : parents.lData[findex-1];
		
        // count non-missing entries across states for this node
		for (long lev = 0; lev < family_nlevels.lData[findex]; lev++) {
			obs_total += observed_values (findex, lev);    
		}
		
		if (obs_total < MIN_SAMPLE_SIZE) {   
			// not enough observations - sample from prior given hyperparameters
			_Matrix imaginary_sample_counts;
			
			CreateMatrix (&imaginary_sample_counts, family_nlevels.lData[findex], 1, false, true, false);
			
			for (long lev = 0; lev < family_nlevels.lData[findex]; lev++) {
				imaginary_sample_counts.Store (lev, 0, prior_sample_size (fnode, 0) / family_nlevels.lData[findex] + DIRICHLET_FLATTENING_CONST);
			}
			
			// random deviate from Dirichlet distribution parameterized by imaginary sample size hyperparameters
			_Matrix * result = (_Matrix *) imaginary_sample_counts.DirichletDeviate();
			
			for (long lev = 0; lev < family_nlevels.lData[findex]; lev++) {
				observed_values.Store (findex, lev, (*result) (lev,0));
			}
			
			DeleteObject (result);
		} else {
			// normalize counts to get empirical distribution
			for (long lev = 0; lev < family_nlevels.lData[findex]; lev++) {
				observed_values.Store (findex, lev, observed_values(findex, lev) / (_Parameter) obs_total);
			}
		}
    }
	
	
	
    // initialize missing entries to random assignments based on observed cases or prior info
    for (long fnode, row, col, missing_idx = 0; missing_idx < is_missing.lLength; missing_idx++) {
        row     = is_missing.lData[missing_idx] / family_size;
        col     = is_missing.lData[missing_idx] % family_size;
        fnode   = (col == 0) ? node_id : parents.lData[col-1];
		
		urn = genrand_real2 ();
		
		for (long level = 0; level < family_nlevels.lData[col]; level++) {
			if (urn < observed_values (col, level)) {
				data_deep_copy.Store (row, col, (_Parameter) level);
				break;
			} else {
				urn -= observed_values (col, level);
			}
		}
    }
	
	
	// tally N_ijk summary statistics for imputed data
	for (long pa_index, row = 0; row < data_deep_copy.GetHDim(); row++) {
		pa_index    = 0;
		child_state = data_deep_copy (row, 0);
		for (long par = 0; par < parents.lLength; par++) {
			pa_index += (long) data_deep_copy (row, par+1) * multipliers[par];
		}
		n_ijk.Store (pa_index, (long) child_state, n_ijk(pa_index, (long) child_state) + 1);
		n_ij.Store (pa_index, 0, n_ij(pa_index, 0) + 1);
	}
	
	
	// Gibbs sampler
	for (long iter = 0; iter < impute_maxsteps + impute_burnin; iter++) {
		// this line was commented out in bgm2.cpp, but I can't remember why - afyp
		//is_missing.Permute(1);  // shuffle list of missing entries
		
		log_score = K2Score (r_i, n_ij, n_ijk); // compute initial score
		
		for (long row, col, pa_index, missing_idx = 0; missing_idx < is_missing.lLength; missing_idx++) {
			row         = is_missing.lData[missing_idx] / family_size;
			col         = is_missing.lData[missing_idx] % family_size;
			pa_index    = 0;
			denom       = 0.;
			
			
			// determine parent combination for this row -- use [pa_indices] object instead?  AFYP
			for (long par = 0; par < parents.lLength; par++) {
				pa_index += ((long) data_deep_copy (row, par+1)) * multipliers.lData[par];
			}
			
			reassign_probs->ZeroUsed(); // reset _GrowingVector
			
			if (col == 0) { // child node, reassign N_ijk's, keep N_ij's constant
				child_state = data_deep_copy (row, col);
				log_score -= log(n_ijk(pa_index, child_state));
				
				n_ijk.Store (pa_index, (long) child_state, n_ijk(pa_index, (long) child_state) - 1);
				
				// compute probabilities for all instantiations of child node
				for (long k = 0; k < r_i; k++) {
					reassign_probs->Store (log_score + log (n_ijk(pa_index, k) + 1));
				}
				
				denom = LogSumExpo (reassign_probs);
				
				// random reassignment
				urn = genrand_real2 ();
				
				for (long lev = 0; lev < r_i; lev++) {
					if ( urn  <  ( this_prob = exp ((* (_Matrix *) reassign_probs) (lev,0) - denom) )) {
						n_ijk.Store (pa_index, lev, n_ijk(pa_index,lev) + 1);
						data_deep_copy.Store (row, col, lev);
						
						log_score = (* (_Matrix *) reassign_probs) (lev, 0);
						
						/* // I don't see why this is necessary - this would short-circuit the burnin!
						if (missing_idx == is_missing.lLength - 1) {
							vector_of_scores->Store ((* (_Matrix *) reassign_probs) (lev,0));
						}
						 */
						
						break;
					}
					
					urn -= this_prob;
				}
				
			} else {    // parent node, reassign both N_ij's AND N_ijk's
				parent_state    = data_deep_copy (row, col);
				child_state     = data_deep_copy (row, 0);
				
				log_score   -= log (n_ijk (pa_index, (long) child_state));
				log_score   += log (n_ij  (pa_index, 0) + 1);
				
				n_ij.Store  (pa_index, 0, n_ij(pa_index,0) - 1);
				n_ijk.Store (pa_index, (long) child_state, n_ijk(pa_index, (long) child_state) - 1);
				
				pa_index    -= parent_state * multipliers.lData[col-1];
				
				
				for (long pa_temp, lev = 0; lev < family_nlevels.lData[col]; lev++) {
					pa_temp = pa_index + lev * multipliers.lData[col-1];
					reassign_probs->Store (log_score + log (n_ijk (pa_temp, (long) child_state) + 1)
										   - log (n_ij (pa_temp, 0) + 1) ); // why was this adding 2?
				}
				
				
				// compute denominator and convert log-values to probability vector
				denom = LogSumExpo (reassign_probs);
				
				// re-assign state
				urn = genrand_real2();  // a uniform random number within interval [0,1)
				
				for (long lev = 0; lev < family_nlevels.lData[col]; lev++) {
					this_prob = exp((* (_Matrix *) reassign_probs) (lev, 0) - denom);
					
					if (urn < this_prob) {
						pa_index += lev * multipliers.lData[col-1];
						
						n_ij.Store (pa_index, 0, n_ij(pa_index,0) + 1);
						n_ijk.Store (pa_index, (long) child_state, n_ijk(pa_index, (long) child_state) + 1);
						
						log_score = (* (_Matrix *) reassign_probs) (lev, 0);
						
						data_deep_copy.Store (row, col, lev);
						
						/*
						if (missing_idx == is_missing.lLength - 1) {
							vector_of_scores->Store ((* (_Matrix *) reassign_probs) (lev,0));
						}
						 */
						break;
					}
					urn -= this_prob;
				}
			}
		}
		// end loop over missing values
		
		// record current log score
		if (iter > impute_burnin && ((long)(iter-impute_burnin) % sampling_interval == 0)) {
			vector_of_scores->Store (log_score);
		}
	}
	
	
	ReportWarning (_String("ImputeDiscreteNodeScore(") & node_id & "<-" & (_String *) parents.toStr() & ":\n");
	long gv_used = vector_of_scores->GetUsed();
	for (long i = 0; i < gv_used; i++) {
		ReportWarning ( _String (",") & (*vector_of_scores)(i,0) );
		log_score += (*vector_of_scores)(i,0);
	}
	
	log_score /= gv_used;
	
	
    // compute the average of sampled log-likelihoods
    //log_score = LogSumExpo (vector_of_scores) - log((double)vector_of_scores->GetUsed());
	
	
    DeleteObject (vector_of_scores);
    DeleteObject (reassign_probs);
	
    return (log_score);
}



//___________________________________________________________________________________________________
_Parameter _BayesianGraphicalModel::ImputeCGNodeScore (long node_id, _SimpleList & parents)
{
    /* ---------------------------------------------------------------------------------------
        ImputeNodeScore ()
            Impute (i.e., "fill in") missing values in the data.  More exactly, take the
            expectation of node score over imputations of the data sampled using Gibbs
            sampling.  Missing values are flagged for discrete data by negative values.
            For continuous data, I haven't settled on an unambiguous flag and I'm leaving this
            up to the user for now (BGM_CONTINUOUS_MISSING_VALUE).
       --------------------------------------------------------------------------------------- */

	ReportWarning (_String("ImputeCGNodeScore(") & node_id & "<-" & (_String *) parents.toStr() & ":\n");
	
    long            num_parent_combos   = 1,
                    r_i                   = num_levels.lData[node_id],	// number of levels for child node
                    max_num_levels        = r_i,
                    family_size           = parents.lLength + 1,
                    sampling_interval,
					batch_count;


    _SimpleList     multipliers (1, 1, 0),  // length constructor, populator

                    dparents, cparents,
                    parents_by_nodetype,    // indexes parents separately in mixed family (discete and CG)

                    family_nlevels (family_size, 0, 0),     // for fast look-up
                    family_datatype (family_size, 0, 0),

                    is_missing,             // store linear indices to missing entries in deep copy

                    pa_indices;				// store discrete parent combination per row


    _Matrix         n_ijk, n_ij,        // summary statistics for discrete nodes
                    mu, tau,            // prior hyperparameters for CG nodes

                    data_deep_copy,     // make a deep copy of the relevant columns of the data matrix
                    observed_values;


    _GrowingVector  * vector_of_scores  = new _GrowingVector(),     // store scores sampled during imputation
					* reassign_probs    = new _GrowingVector();

    _Parameter      log_score           = 0,

                    impute_maxsteps, impute_burnin, impute_samples, // HBL settings

                    parent_state, child_state,
                    denom,
                    // prior hyperparameters for CG nodes
                    rho = prior_sample_size (node_id, 0) > 0 ? (prior_sample_size (node_id, 0) / num_parent_combos) : 1.0,
                    phi = prior_scale (node_id, 0);


    double          urn;            // uniform random number

	
	// parameters for CG
	_Matrix rho_mx (1, 1, false, true);
	rho_mx.Store(0, 0, rho);
	ReportWarning (_String("rho_mx: ") & (_String *) rho_mx.toStr());
	
	_Matrix phi_mx (1, 1, false, true);
	phi_mx.Store(0, 0, phi);
	
	_Matrix * iw_ptr;
	_Parameter iw_deviate;
	

    // set Gibbs sampler parameters from batch language definitions
    checkParameter (_HYBgm_IMPUTE_MAXSTEPS, impute_maxsteps, 0);
    checkParameter (_HYBgm_IMPUTE_BURNIN, impute_burnin, 0);
    checkParameter (_HYBgm_IMPUTE_SAMPLES, impute_samples, 0);

    if (impute_maxsteps <= 0 || impute_burnin < 0 || impute_samples <= 0 || impute_samples > impute_maxsteps) {
        WarnError (_String("ERROR: Invalid IMPUTE setting(s) in ImputeNodeScore()"));
        return 0.;
    }


    // number of steps between samples
    sampling_interval = (long) (impute_maxsteps / impute_samples);
	//ReportWarning(_String("ImputeCGNodeScore sampling interval set to ") & sampling_interval);


    // partition parent list into discrete and continuous; record stuff while we're at it
    num_parent_combos = 1;
    family_nlevels.lData[0] = r_i;

    for (long nlev, p = 0; p < parents.lLength; p++) {
        if (node_type.lData[parents.lData[p]] == 0) {
			// parent is discrete-valued
            dparents << parents.lData[p];
            parents_by_nodetype << (dparents.lLength-1);	// stores index into dparents

            nlev = family_nlevels.lData[p+1] = num_levels.lData[parents.lData[p]];

            num_parent_combos *= nlev;
            multipliers << num_parent_combos;

            if (nlev > max_num_levels)  {
                max_num_levels = nlev;
            }
        } else { 
			// parent is continuous-valued
            cparents << parents.lData[p];
            parents_by_nodetype << (cparents.lLength-1);
            family_nlevels.lData[p+1] = 0;
        }
    }



    // allocate space to matrices
    CreateMatrix (&n_ijk, num_parent_combos, r_i, false, true, false);
    CreateMatrix (&n_ij, num_parent_combos, 1, false, true, false);
    CreateMatrix (&data_deep_copy, theData.GetHDim(), family_size, false, true, false);

    pa_indices.Populate(theData.GetHDim(), 0, 0);


    // allocate space to _GrowingVector object -- used for imputing discrete nodes only
    for (long pa = 0; pa < num_parent_combos; pa++) {
        reassign_probs->Store (0.);
    }



    // for discrete nodes, the j-th column tallies the number of observed instances of j-th level
    // for continuous nodes, store summary statistics (sample size, mean, variance)
	// therefore, the minimum number of columns is 3
    CreateMatrix (&observed_values, family_size, (max_num_levels > 3) ? max_num_levels : 3, false, true, false);

	
    // make deep copy, annotate which entries are missing, and record empirical distributions
    for (long row = 0; row < theData.GetHDim(); row++) {
        // findex = family-wise index
        for (long fnode, findex = 0; findex < family_size; findex++) {
            fnode = (findex == 0) ? node_id : parents.lData[findex-1];  // parent node, get network-wide index
            _Parameter fstate = theData (row, fnode);
            data_deep_copy.Store (row, findex, fstate);

            if (node_type.lData[fnode] == 0) {  // discrete node
                if (fstate < 0) {
                    is_missing << row * family_size + findex; // note this is family level indexing
                } else {
					// use state into index into vector and increment counter
                    observed_values.Store (findex, fstate, observed_values(findex,fstate) + 1);
                }
            } else {                            // continuous node
                if (fstate == continuous_missing_value) {
                    is_missing << row * family_size + findex;
                } else {
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

	//ReportWarning (_String ("is_missing = ") & (_String *) is_missing.toStr() );
	ReportWarning (_String ("observed values matrix: ") & (_String *) observed_values.toStr() );


    // make a decision whether to use empirical or prior distribution to do initial imputations
    for (long obs_total, fnode, findex = 0; findex < family_size; findex++) {
        obs_total = 0;
        fnode = (findex == 0) ? node_id : parents.lData[findex-1];

        if (node_type.lData[fnode] == 0) {
			// DISCRETE NODE
            for (long lev = 0; lev < family_nlevels.lData[findex]; lev++) {
                obs_total += observed_values (findex, lev);    // sum tallies across levels
            }

            if (obs_total < MIN_SAMPLE_SIZE) {   
				// not enough observations - sample from prior given hyperparameters
				_Matrix imaginary_sample_counts;
				
				CreateMatrix (&imaginary_sample_counts, family_nlevels.lData[findex], 1, false, true, false);
				
				for (long lev = 0; lev < family_nlevels.lData[findex]; lev++) {
					imaginary_sample_counts.Store (lev, 0, prior_sample_size (fnode, 0) / family_nlevels.lData[findex] + DIRICHLET_FLATTENING_CONST);
				}
				
				// random deviate from Dirichlet distribution parameterized by imaginary sample size hyperparameters
				_Matrix * result = (_Matrix *) imaginary_sample_counts.DirichletDeviate();
				
                for (long lev = 0; lev < family_nlevels.lData[findex]; lev++) {
                    observed_values.Store (findex, lev, (*result) (lev,0));
                }
				
				DeleteObject (result);
            } else {
				// normalize counts to get empirical distribution
                for (long lev = 0; lev < family_nlevels.lData[findex]; lev++) {
                    observed_values.Store (findex, lev, observed_values(findex, lev) / (_Parameter) obs_total);
                }
            }
        } else {
			// CONTINUOUS NODE
            if (observed_values(findex, 0) < MIN_SAMPLE_SIZE) {	
				
				// unknown variance drawn from scaled inverse chi-square distribution
				_Matrix		temp1, temp2;
				CreateMatrix (&temp1, 1, 1, false, true, false);
				temp1.Store (0, 0, prior_scale(fnode,0));
				CreateMatrix (&temp2, 1, 1, false, true, false);
				temp2.Store (0, 0, prior_sample_size(fnode,0));
				
				_Matrix *	result = (_Matrix *) temp1.InverseWishartDeviate (temp2);
				
                observed_values.Store (findex, 2, (*result)(0,0));
				observed_values.Store (findex, 1, gaussDeviate() * sqrt( (*result)(0,0) ) + prior_mean(fnode,0));
				
				DeleteObject (result);
				
			} else { 
				// compute empirical mean and s.d.
				
					// E[x] = sum(x) / len(x)
                observed_values.Store (findex, 1, observed_values(findex,1) / observed_values(findex,0));
					// Var[x] = E[x^2] - (E[x])^2
                observed_values.Store (findex, 2, sqrt(observed_values(findex,2) / observed_values(findex,0)
                                                       - observed_values(findex,1) * observed_values(findex,1)) );
            }
        }
    }

	ReportWarning (_String ("observed values matrix 2: ") & (_String *) observed_values.toStr() );
	
	ReportWarning (_String ("data before impute:\n") & (_String *) data_deep_copy.toStr());

    // initialize missing entries to random assignments based on observed cases or prior info
    for (long fnode, row, col, missing_idx = 0; missing_idx < is_missing.lLength; missing_idx++) {
        row     = is_missing.lData[missing_idx] / family_size;
        col     = is_missing.lData[missing_idx] % family_size;
        fnode   = (col == 0) ? node_id : parents.lData[col-1];

        if (node_type.lData[fnode] == 0) {  // discrete node
            urn = genrand_real2 ();

            for (long level = 0; level < family_nlevels.lData[col]; level++) {
                if (urn < observed_values (col, level)) {
                    data_deep_copy.Store (row, col, (_Parameter) level);
                    break;
                } else {
                    urn -= observed_values (col, level);
                }
            }
        } else {                            // continuous node
            data_deep_copy.Store (row, col, gaussDeviate() * observed_values(col,2) + observed_values(col,1) );
        }
    }

	ReportWarning (_String ("data after impute:\n") & (_String *) data_deep_copy.toStr());


	long            k           = cparents.lLength,
					pa_index;

	_Parameter      lk_ratio,
					next_log_score;     // log-likelihood of the proposed step

	_Matrix         log_scores_by_pa (num_parent_combos, 1, false, true),   // track separately to make updating easier
					next_log_score_by_pa (num_parent_combos, 1, false, true);



	// set mean intercept/regression coeff. hyperparameter vector for CG prior
	CreateMatrix (&mu, k+1, 1, false, true, false);

	mu.Store (0, 0, prior_mean (node_id, 0));   // intercept
	for (long i = 1; i < k+1; i++) {
		mu.Store (i,0,0);    // set prior expectation of regression coeffs. to 0
	}


	// set precision hyperparameter for CG prior
	CreateMatrix (&tau, k+1, k+1, false, true, false);

	tau.Store (0, 0, prior_precision (node_id, 0));
	for (long i = 1; i < k+1; i++) {
		tau.Store (i, i, prior_precision(cparents.lData[i-1],0));
	}




	// partition data set by discrete parent state combination
	for (long obs = 0; obs < theData.GetHDim(); obs++) {
		pa_index = 0;
		
		for (long this_parent, findex = 1; findex < family_size; findex++)
		{
			this_parent = parents.lData[findex-1];
			
			if (node_type.lData[this_parent] == 0) {
				// discrete parent
				pa_index += data_deep_copy (obs, findex) * multipliers.lData[parents_by_nodetype.lData[findex-1]];
			}
		}

		pa_indices.lData[obs] = pa_index;
		n_ij.Store (pa_index, 0, n_ij(pa_index,0) + 1);
	}



	// calculate CG summary statistics and compute initial scores
	log_score = 0.;

	for (long pa = 0; pa < num_parent_combos; pa++) {
		if (n_ij(pa,0) == 0) {
			continue; // this batch is empty
		}
		
		// these are single-use matrices, reset with every loop
		_Matrix zbpa (n_ij(pa,0), k+1, false, true);	// (1, CG parent node values) in batch
		_Matrix yb (n_ij(pa,0), 1, false, true);		// CG child node values in batch
		

		// populate zbpa with continuous parent states
		batch_count = 0;
		for (long obs = 0; obs < pa_indices.lLength; obs++) {
			if (pa_indices.lData[obs] == pa) {
				// (1, x_0, x_1, ..., x_N)
				zbpa.Store (batch_count, 0, 1); // intercept

				for (long findex = 1; findex < family_size; findex++) {
					if (node_type.lData[parents.lData[findex-1]] == 1) {    // CG parent
						zbpa.Store (batch_count, parents_by_nodetype.lData[findex-1], data_deep_copy(obs, findex));
					}
				}

				yb.Store (batch_count, 0, data_deep_copy(obs, 0));  // child node state

				batch_count++;
			}
		}

		log_scores_by_pa.Store (pa, 0, BottcherScore (yb, zbpa, tau, mu, rho, phi, n_ij(pa,0)));
		log_score += log_scores_by_pa (pa, 0);
	}



	// Gibbs sampling over missing entries
	for (long iter = 0; iter < impute_burnin + impute_maxsteps; iter++) {
		//is_missing.Permute(1);

		for (long row, col, missing_idx = 0; missing_idx < is_missing.lLength; missing_idx++) {
			// indices into data_deep_copy for this missing value
			row         = is_missing.lData[missing_idx] / family_size;
			col         = is_missing.lData[missing_idx] % family_size;
			denom       = 0.;
			pa_index    = pa_indices.lData[row];


			// if child node, then row remains in same 'pa' batch, Z^b_pa and N_ij unchanged
			if (col == 0) {
				
				// draw new state from posterior of CG node
								
				iw_ptr = (_Matrix *) phi_mx.InverseWishartDeviate (rho_mx);
				iw_deviate = (*iw_ptr)(0,0);
				ReportWarning(_String("sigma=") & iw_deviate);
				

				_Matrix zbpa (n_ij(pa_index,0), k+1, false, true);
				_Matrix yb (n_ij(pa_index,0), 1, false, true);
				
				// collect all continuous parent states under given discrete parent combination (pa_index)
				batch_count = 0;
				for (long obs = 0; obs < pa_indices.lLength; obs++) {
					if (obs == row) continue; // skip the current case
					if (pa_indices.lData[obs] == pa_index) {
						zbpa.Store (batch_count, 0, 1);
						for (long findex = 1; findex < family_size; findex++) {
							if (node_type.lData[parents.lData[findex-1]] == 1) {
								zbpa.Store (batch_count, parents_by_nodetype.lData[findex-1], data_deep_copy(obs, findex));
							}
						}
						yb.Store (batch_count, 0, data_deep_copy(obs, 0));
						batch_count++;
					}
				}
				
				ReportWarning(_String("zbpa=") & (_String *) zbpa.toStr());
				zbpa.Transpose();
				ReportWarning(_String("zbpa=") & (_String *) zbpa.toStr());
				zbpa *= yb;
				ReportWarning(_String("zbpa=") & (_String *) zbpa.toStr());
				

				next_log_score_by_pa.Store (pa_index, 0, BottcherScore (yb, zbpa, tau, mu, rho, phi, (long)n_ij(pa_index,0)));
				next_log_score  = log_score - log_scores_by_pa (pa_index, 0) + next_log_score_by_pa(pa_index, 0);

				lk_ratio        = exp(log_score - next_log_score);


				// accept this step?
				if (lk_ratio >= 1. || genrand_real2() < lk_ratio) {
					log_score = next_log_score;

					// update log score contribution for this parent combo (pa)
					log_scores_by_pa.Store (pa_index, 0, next_log_score_by_pa(pa_index, 0));
				} else {
					// revert to previous state
					data_deep_copy.Store (row, col, child_state);
				}
			}

			/* ---------- MISSING VALUE AT PARENT NODE (col > 0) ---------- */
			else {
				// if parent is discrete, then this case is moved to another batch (pa) -- MORE THAN ONE BATCH IS AFFECTED
				if (node_type.lData[parents.lData[col-1]] == 0) {
					
					parent_state    = data_deep_copy (row, col);

					// first compute the likelihood of the original batch minus this case
					n_ij.Store  (pa_index, 0, n_ij(pa_index,0) - 1);
					
					_Matrix zbpa (n_ij(pa_index,0), k+1, false, true);
					_Matrix yb (n_ij(pa_index,0), 1, false, true);

					batch_count = 0;
					for (long obs = 0; obs < data_deep_copy.GetHDim(); obs++) {
						if (obs == row) {
							continue;    // skip this case
						}

						if (pa_indices.lData[obs] == pa_index) {
							zbpa.Store (batch_count, 0, 1);

							for (long findex = 1; findex < family_size; findex++) {
								if (node_type.lData[parents.lData[findex-1]] == 1) {
									zbpa.Store (batch_count, parents_by_nodetype.lData[findex-1], data_deep_copy(obs, findex));
								}
							}

							yb.Store (batch_count, 0, data_deep_copy(obs, 0));

							batch_count++;
						}
					}

					next_log_score_by_pa.Store (pa_index, 0, BottcherScore (yb, zbpa, tau, mu, rho, phi, (long) n_ij (pa_index, 0)));
					next_log_score = log_score - log_scores_by_pa (pa_index,0) + next_log_score_by_pa (pa_index,0);


					// random re-assignment of discrete parent node, compute likelihood of second batch
					
					urn = genrand_real2() * (1.0 - observed_values(col, parent_state)); // rescale to omit original state

					for (long lev = 0; lev < family_nlevels.lData[col]; lev++) {
						if (lev == parent_state) {
							continue;
						}

						if (urn < observed_values(col, lev)) {
							data_deep_copy.Store (row, col, lev);
							
							/* BREAK */
							
							pa_index += (lev - parent_state) * multipliers.lData[parents_by_nodetype.lData[col-1]];


							// compute summary stats
							n_ij.Store  (pa_index, 0, n_ij(pa_index,0) + 1);
							
							_Matrix zbpa (n_ij(pa_index,0), k+1, false, true);
							_Matrix yb (n_ij(pa_index,0), 1, false, true);

							batch_count = 0;
							
							for (long obs = 0; obs < data_deep_copy.GetHDim(); obs++) {
								// we won't update pa_indices[] for this case, already have pa_index to work with
								if (obs == row || pa_indices.lData[obs] == pa_index) {
									zbpa.Store (batch_count, 0, 1);

									for (long findex = 1; findex < family_size; findex++) {
										if (node_type.lData[parents.lData[findex-1]] == 1) {
											zbpa.Store (batch_count, parents_by_nodetype.lData[findex-1], data_deep_copy(obs, findex));
										}
									}
									yb.Store (batch_count, 0, data_deep_copy(obs, 0));
									batch_count++;
								}
							}

							next_log_score_by_pa.Store (pa_index, 0, BottcherScore (yb, zbpa, tau, mu, rho, phi, (long) n_ij (pa_index, 0)));
							break;
						} else {
							urn -= observed_values(col, lev);
						}
					}

					next_log_score -= log_scores_by_pa (pa_index, 0);   // remove old score for second batch
					next_log_score += next_log_score_by_pa (pa_index, 0);               // add new scores for both batches

					lk_ratio        = exp(log_score - next_log_score);
					if (lk_ratio >= 1. || genrand_real2() < lk_ratio) {
						log_score = next_log_score;
						log_scores_by_pa.Store (pa_index, 0, next_log_score_by_pa (pa_index,0) );
						log_scores_by_pa.Store (pa_indices.lData[row], 0, next_log_score_by_pa (pa_indices.lData[row],0) ); // contains the old index

						pa_indices.lData[row] = pa_index;   // replace old pa index with new
					} else {
						// revert to previous state
						data_deep_copy.Store (row, col, parent_state);

						n_ij.Store (pa_index, 0, n_ij(pa_index,0) - 1);
						n_ij.Store (pa_indices.lData[row], 0, n_ij(pa_indices.lData[row],0) + 1);
						pa_index = pa_indices.lData[row];
					}
				}

				// otherwise, if parent is continuous then just update Z_bpa; case stays in the same batch
				else {
					// random draw from Gaussian distribution centered at previous value
					data_deep_copy.Store (row, col, (gaussDeviate() + data_deep_copy(row,col)) * observed_values(0,2));

					_Matrix zbpa (n_ij(pa_index,0), k+1, false, true);
					_Matrix yb (n_ij(pa_index,0), 1, false, true);

					batch_count = 0;
					
					for (long obs = 0; obs < data_deep_copy.GetHDim(); obs++) {
						if (pa_indices.lData[obs] == pa_index) {
							zbpa.Store (batch_count, 0, 1);

							for (long findex = 1; findex < family_size; findex++) {
								if (node_type.lData[parents.lData[findex-1]] == 1) {
									zbpa.Store (batch_count, parents_by_nodetype.lData[findex-1], data_deep_copy(obs, findex));
								}
							}

							yb.Store (batch_count, 0, data_deep_copy(obs, 0));

							batch_count++;
						}
					}

					next_log_score_by_pa.Store (pa_index, 0, BottcherScore (yb,zbpa,tau,mu,rho,phi, (long)n_ij(pa_index,0)));
					next_log_score = log_score - log_scores_by_pa(pa_index,0) + next_log_score_by_pa(pa_index,0);

					lk_ratio = exp(log_score - next_log_score);

					if (lk_ratio >= 1. || genrand_real2() < lk_ratio) {
						log_score = next_log_score;
						log_scores_by_pa.Store (pa_index, 0, next_log_score_by_pa(pa_index,0));
						pa_indices.lData[row] = pa_index;
					} else {
						data_deep_copy.Store (row, col, parent_state);
					}
				}
			}
		}
		// end loop over missing entries

		
		if (iter > impute_burnin && ((long)(iter - impute_burnin) % sampling_interval == 0) ) {
			vector_of_scores->Store (log_score);
			//ReportWarning (_String("ImputeCGNodeScore: ") & log_score);
		}
	}
	// end sampler

	
	long gv_used = vector_of_scores->GetUsed();
	for (long i = 0; i < gv_used; i++) {
		//ReportWarning ( _String (",") & (*vector_of_scores)(i,0) );
		log_score += (*vector_of_scores)(i,0);
	}
	ReportWarning(_String("chain = ") & (_String *) vector_of_scores->toStr());
	
	log_score /= gv_used;
    // compute the average of sampled log-likelihoods
    //log_score = LogSumExpo (vector_of_scores) - log((double)vector_of_scores->GetUsed());

    DeleteObject (vector_of_scores);
    DeleteObject (reassign_probs);

    return (log_score);
}




#ifdef __NEVER_DEFINED__
// __________________________________________________________________________________________
void _BayesianGraphicalModel::ComputeParameters(_Matrix * structure)
{
    /* -----------------------------------------------------------------------
        ComputeParameters()

            Store _Formula objects into _AssociativeArray, each represents
            a posterior probability distribution over network parameters for
            each node in the network, given network structure.

            These distributions are defined by hyperparameters that in turn
            are defined by prior assumptions and data.

            Discrete node has conditional probability table (π):
                π_ijk ~ Dirichlet(n_ijk + a_ijk)

            where n_ijk is from data and a_ijk is prior.

            Continuous node has intercept (m) and k-vector of regression
            coefficients (ß), where k is number of continuous parents:
                (m,ß) ~ k-Gauss(mu', sigma/tau')
                sigma ~ Inv-Gamma(rho'/2, phi'/2)


       ----------------------------------------------------------------------- */

    _Matrix     n_ijk, n_ij,
                params;

    _Formula    * aFormula;



    // loop through nodes in the network
    for (long node = 0; node < num_nodes; node++) {
        long            r_i         = num_levels.lData[node];

        _SimpleList     parents;


        // find parents in graph
        for (long other_node = 0; other_node < num_nodes; node++) {
            if (other_node == node) {
                continue;
            }

            if (theStructure(other_node, node) == 1) {
                parents << other_node;
            }
        }


        // this is a root node without parents
        if (parents.lLength == 0) {
            if (node_type.lData[node] == 0) {   /* discrete node */
                long            a_ijk   = (prior_sample_size(node,0) == 0) ? 1 : prior_sample_size(node,0);
                _Parameter      n_i     = 0;


                // tally Dirichlet hyperparameters
                CreateMatrix (&n_ijk, 1, r_i, false, true, false);

                for (long k, obs = 0; obs < theData.GetHDim(); obs++) {
                    k = theData(obs, node); // child state
                    n_ijk.Store (0, k, n_ijk(0,k) + a_ijk);
                    n_i += n_ijk(0,k);
                }


                // store posterior distribution
                CreateMatrix (&params, 1, r_i, false, true, false);
                params.Convert2Formulas();

                for (long k = 0; k < r_i; k++) {
                    checkPointer (aFormula = new _Formula (new _Constant (n_ijk(0,k)/n_i), false));
                    params.StoreFormula(0, k, *aFormula);
                }
            } else if (node_type.lData[node] == 1) {

            } else {
                WarnError (_String ("Unsupported data type detected for node ") & node & " in ComputeParameters()");
            }
        }

        else {
            long    num_parent_combos   = 1,
                    pa_index;

            _SimpleList multipliers ((long)1);

            // tabulate parent combinations
            for (long par = 0; par < parents.lLength; par++) {
                num_parent_combos *= num_levels.lData[parents.lData[par]];
                multipliers << num_parent_combos;
            }

            CreateMatrix (&n_ijk, num_parent_combos, r_i, false, true, false);
            CreateMatrix (&n_ij, num_parent_combos, 1, false, true, false);


            // tally N_ijk and N_ij's
            for (long child_state, obs = 0; obs < theData.GetHDim(); obs++) {
                child_state =  theData(obs, node);
                pa_index    = 0;

                for (long par = 0; par < parents.lLength; par++) {
                    pa_index += theData(obs, parents.lData[par]) * multipliers.lData[par];
                }

                n_ijk.Store (pa_index, child_state, n_ijk(pa_index, child_state) + 1);
                n_ij.Store (pa_index, 0, n_ij(pa_index,0) + 1);
            }


            if (node_type.lData[node] == 0) {
                // that's enough information to compute parameters on discrete nodes

            } else if (node_type.lData[node] == 1) {

            } else {
                WarnError (_String ("Unsupported data type detected for node ") & node & " in ComputeParameters()");
            }
        }
    }
}
#endif







