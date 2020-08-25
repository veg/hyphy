/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon42@uwo.ca)
 Steven Weaver (sweaver@temple.edu)
 
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


#include "function_templates.h"
#include "global_things.h"
#include "bayesgraph.h"

#include "ntuplestorage.h"

using namespace hy_global;


//___________________________________________________________________________________________
void _BayesianGraphicalModel::SerializeBGM (_StringBuffer & rec) {
    /* --------------------------------------------------------------------
        ExportModel()
            Export the network structure and parameters as an associative
            array.
       -------------------------------------------------------------------- */
  
    static const _String kNodeType ("NodeType") , kParents ("Parents"), kCPDFs ("CPDFs");

    ReportWarning (_String("Entered _BayesianGraphicalModel::ExportModel()"));

    if (theData.GetHDim() > 0) {
        _String     bgm_name = *(_String *) bgmNamesList (bgmList._SimpleList::Find((long)this)) & "_export";

        rec << bgm_name << "={};\n";
 
        for (long node = 0; node < num_nodes; node++) {
          
           // start the dictionary entry for this node
           rec << bgm_name << "['" << node << "'] = {\n"
               << kNodeType.Enquote() << ':' << _String (node_type.get (node)) << ",\n"
               << kParents.Enquote() << ':' ;
          
            _SimpleList parents, d_parents, c_parents;
          
            for (long pnode = 0; pnode < num_nodes; pnode++) {
              if (pnode == node) {
                continue;
              }
              if (theStructure(pnode,node) == 1) {
                parents << pnode;
                
                if (is_node_discrete(pnode)) {
                  d_parents << pnode;
                } else {
                  c_parents << pnode;
                }
              }
            }
          
            rec << '{';
            rec.AppendNewInstance((_String*)parents.toStr()) << "},\n";
            rec << kCPDFs.Enquote() << " : {{";
 
 
            // posterior probability density functions (PDFs) conditional on child node state (indexed by i)
            if (is_node_discrete (node)) {
                // discrete child node, compute Dirichlet hyperparameters
                _Matrix     n_ij, n_ijk;


                // tally cases by discrete parent combination and child state
                UpdateDirichletHyperparameters (node, d_parents, &n_ij, &n_ijk);


                // report posterior Dirichlet PDFs
                for (long j = 0; j < n_ij.GetHDim(); j++) {
                    rec << "\"Random({{";
                    for (long k = 0; k < num_levels.list_data[node]; k++) { // k indexes child node states
                        rec << _String(n_ijk(j,k));
                        if (k < num_levels.list_data[node]-1) {
                            rec << ',';
                        }
                    }
                    rec << "}},{'PDF':'Dirichlet'})\"";
                    if (j < n_ij.GetHDim()-1) {
                        rec << ',';
                    }
                }
            } else if (is_node_continuous(node)) {  // CG child node
 
 
                 long            k                 = c_parents.countitems(),        // a useful short-hand
                                 num_parent_combos = 1L;
              
                hyFloat         rho_prior   = prior_sample_size (node, 0) > 0. ? (prior_sample_size (node, 0) / num_parent_combos) : 1.0,
                                phi_prior = prior_scale (node, 0);

                _Matrix         mu_prior  (k+1L, 1L, false, true),
                                tau_prior (k+1L, k+1L, false, true);

 
              _SimpleList     multipliers (ComputeParentMultipliers(d_parents));
              if (multipliers.nonempty()) {
                num_parent_combos = multipliers.GetElement (-1L);
              }
              
              _SimpleList     n_ij,
                              pa_indexing;
              
                num_parent_combos = ReindexParentObservations (d_parents, n_ij, pa_indexing);

                tau_prior.Store (0L, 0L, prior_precision (node, 0));
                for (long row = 1L; row <= k; row++) {
                  tau_prior.Store (row, row, prior_precision (c_parents.get(row-1L),0L));
                }
        


                mu_prior.Store (0, 0, prior_mean (node, 0));                // prior intercept
 

                // for every discrete parent config (on which CG marginal distribution is conditioned)
                for (long pa = 0; pa < num_parent_combos; pa++) {
                    // compute summary statistics
                    _Matrix     zbpa (n_ij.get(pa), k+1, false, true),
                                yb (n_ij.get(pa), 1, false, true);

                    for (long count_n = 0L, obs = 0L; obs < theData.GetHDim(); obs++) {
                        if (pa_indexing.get(obs) == pa) {
                            zbpa.Store (count_n, 0, 1);
                            for (long cpar = 0; cpar < k; cpar++) {
                                zbpa.Store (count_n, cpar+1, theData(obs, c_parents.get(cpar)));
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

                    _Matrix     * tauinv = (_Matrix *) tau.Inverse(nil);
                    temp = *tauinv;
                    temp *= mu;
                    mu = temp;

                    DeleteObject (tauinv);
					
                    hyFloat  rho = rho_prior + n_ij.list_data[pa],
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
                    rec << ",{'PDF':'Gaussian','ARG0':Random({{";
                    rec << _String(phi);
                    rec << "}},{'PDF':'InverseWishart','ARG0':{{";
                    rec << _String(rho);
                    rec << "}}})*";
                    rec << _String((_String*)tau.Inverse(nil)->toStr());
                    rec << "})\"";

                    if (pa < num_parent_combos-1) {
                        rec << ',';
                    }
                }

            }

            rec << "}}\n" << "};\n";
        }
    } else {
        HandleApplicationError (_String("Cannot export network parameters, this _BayesianGraphicalModel object has no data!"));
    }
}




//___________________________________________________________________________________________
bool _BayesianGraphicalModel::ImportCache (_AssociativeList * cache_import) {
    //  Unpack contents of associative list sent from HBL.
    //   use to restore node score cache as an alternative to computing scores
    //   de novo from data -- particularly useful for time-consuming imputation
    //   in cases of missing data
    //ReportWarning (_String("Entered ImportCache() with avl: ") & (_String *) cache_import->toStr());


    if (scores_cached) {
        ReportWarning (_String("WARNING: Overwriting pre-existing node score cache in ") & __PRETTY_FUNCTION__);
    }

    try {
      for (long node = 0; node < num_nodes; node++) {
          _List   *   node_scores = (_List *) node_score_cache.GetItem (node);
          node_scores->Clear();

          for (long num_parents = 0; num_parents <= max_parents.get(node); num_parents++) {
            
              const _String key = _String("Node") & node & "NumParents" & num_parents;
            
              HBLObjectRef cache_value = (HBLObjectRef) cache_import->GetByKey (key);
              if (!cache_value) {
                throw key.Enquote() & " not present in the cache dictionary";
              }

              if (num_parents == 0) {
                  if (cache_value->ObjectClass () == NUMBER) {
                      (*node_scores) && cache_value;
                  } else {
                    throw key.Enquote() & " was expected to be a number";
                 }
              } else {
                if (cache_value->ObjectClass () == MATRIX) {
                  (*node_scores) && cache_value;
                } else {
                  throw key.Enquote() & " was expected to be a matrix";
                }
              }
          }
      }
    } catch (const _String &err) {
      HandleApplicationError(err);
      return false;
    }
  

    scores_cached = true;
    return true;
}


//___________________________________________________________________________________________
bool _BayesianGraphicalModel::ExportCache (_AssociativeList * cache_export) const {
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
            this_list = (_List *) node_score_cache.GetItem (node);

            for (long npar = 0; npar <= max_parents.list_data[node]; npar++) {
                keyString = _String ("Node") & node & "NumParents" & npar;
                _FString    aKey (keyString, false);
              
              

                //ReportWarning (_String("Inserting with key ") & keyString);

                if (npar == 0) {
                    orphan_score = (_Constant *) this_list->list_data[npar];
                    cache_export->MStore (&aKey, orphan_score, true);   // make a dynamic copy
                } else if (npar == 1) {
                    single_parent_scores = (_Matrix *) this_list->list_data[npar];
                    cache_export->MStore (&aKey, single_parent_scores, true);   // make a dynamic copy
                } else {
                    family_scores = (_NTupleStorage *) this_list->list_data[npar];
                    cache_export->MStore (&aKey, family_scores, true);  // make a dynamic copy
                }
            }
        }

        return TRUE;
    } else {
        HandleApplicationError (_String ("Unable to export node score cache, no cache exists!"));
        return FALSE;
    }
}




//___________________________________________________________________________________________
// #define _DEBUG_CCS_
hyFloat _BayesianGraphicalModel::ComputeContinuousScore (long node_id) {
    return ComputeContinuousScore (node_id, theStructure);
}


hyFloat _BayesianGraphicalModel::ComputeContinuousScore (long node_id, _Matrix const &dag) {
    _SimpleList     parents;

    for (long par = 0L; par < num_nodes; par++) {
        // support for discrete parents only
        if (dag(par, node_id) == 1L) {
            parents << par;
        }
    }
    return ComputeContinuousScore (node_id, parents);
}

//___________________________________________________________________________________________

const _SimpleList _BayesianGraphicalModel::ComputeParentMultipliers (_SimpleList const& parents) const {
  long            num_parent_combos = 1L;
  _SimpleList     multipliers ((long)1L);
  multipliers.RequestSpace(parents.countitems());
  
  parents.Each([this, &multipliers, &num_parent_combos] (long parent_index, unsigned long) -> void {
    multipliers << (num_parent_combos *= num_levels.get (parent_index));
  });
  
  return multipliers;
}

//___________________________________________________________________________________________

long            _BayesianGraphicalModel::ReindexParentObservations(_SimpleList const & parents, _SimpleList& n_ij, _SimpleList& pa_indexing) const {
  long obs_count = theData.GetHDim();
  _SimpleList multipliers = ComputeParentMultipliers (parents);
  
  long parent_combos = parents.nonempty() ? multipliers.GetElement(-1L) : 1;
  
  n_ij.Populate (parent_combos, 0, 0);
  pa_indexing.Populate (obs_count, 0, 0);
  
  for (long obs = 0L; obs < obs_count; obs++) {
    
    long     index       = 0L;
    
    parents.Each ([this, obs, &index, &multipliers] (long parent, unsigned long parent_index) -> void {
      index      += (long)this->theData(obs, parent) * multipliers.get (parent_index);
    });
    
    pa_indexing[obs] = index;
    n_ij [index] += 1;
  }
  return parent_combos;
}


//___________________________________________________________________________________________

hyFloat _BayesianGraphicalModel::ComputeContinuousScore (long node_id, _SimpleList const & parents) {
    /* ---------------------------------------------------------------------------
        ComputeContinuousScore()
            Returns the network score based on Bottcher's formulation.
            Note, does not update posterior distributions of hyperparameters,
            this requires a separate albeit highly similar function.
     ----------------------------------------------------------------------------- */

    //ReportWarning (_String ("Called ComputeContinuousScore with ") & node_id & " <- " & (_String *) parents.toStr());
  hyFloat      log_score = 0.;

  try {
         //long            k;
        //_SimpleList     c_parents, d_parents;


        // use cached node scores if available
        if (scores_cached) {
            _List *     scores  = (_List *) node_score_cache.GetItem (node_id);

            if (parents.empty()) {
                 return ((_Constant *) scores->GetItem(0))->Value();
            } else if (parents.lLength == 1) {
                return (*(_Matrix *) scores->GetItem(1)) (parents.get (0), 0);
            } else {
              return ((_NTupleStorage *) scores->GetItem(parents.countitems()))->Retrieve (parents.MapList([node_id] (long par, unsigned long) -> long {
                return par > node_id ? par - 1L : par;
              }));  // using nk-tuple
            }
        }

        if (theData.is_empty() == 0) {
            throw _String("Uh-oh, there's no node score cache nor is there any data matrix to compute scores from!");
        }

    
      // impute score if missing data
      if (has_missing.get(node_id)) {
        return ImputeCGNodeScore (node_id, parents);
      } else {
        if (parents.Any ([this] (long par, unsigned long) -> bool {
          return has_missing.get (par) > 0;
        })) {
          return ImputeCGNodeScore (node_id, parents);
        }
      }
    
        _SimpleList d_parents (parents.Filter([this] (long node, unsigned long)->bool {return is_node_discrete(node);}));
        _SimpleList c_parents (parents.Filter([this] (long node, unsigned long)->bool {return is_node_continuous(node);})),
                    n_ij, pa_indexing;


    
        long continuous_parent_count = c_parents.countitems(),
             num_parent_combos       = ReindexParentObservations (d_parents, n_ij, pa_indexing);

        // if there are discrete parents then integrate over marginal posteriors conditioned on discrete parent states
  
       // set hyperparameter priors -- they will not be modified on pass to BottcherScore() so keep outside loop
      _Matrix     tau (continuous_parent_count+1, continuous_parent_count+1, false, true),
                  mu (continuous_parent_count+1, 1, false, true);

      hyFloat  rho     = prior_sample_size (node_id, 0) > 0 ? (prior_sample_size (node_id, 0) / num_parent_combos) : 1.0,
               phi     = prior_scale (node_id, 0);

      // the matrices are already zero-ed so only set the diagonal
    
      tau.Store (0L, 0L, prior_precision (node_id, 0));
      for (long row = 1L; row <= continuous_parent_count; row++) {
        tau.Store (row, row, prior_precision (c_parents.get(row-1L),0L));
      }
    

      mu.Store (0, 0, prior_mean (node_id, 0));               // prior intercept


      // for every discrete parent combination (on which Gaussian node state is conditioned)
      for (long pa = 0; pa < num_parent_combos; pa++) {
        // collect cases into a batch defined by discrete parent states
        _Matrix     zbpa (n_ij.get(pa), continuous_parent_count+1, false, true),
                    yb (n_ij.get(pa), 1, false, true);

        for (long count_n = 0, obs = 0; obs < theData.GetHDim(); obs++) {
          if (pa_indexing.get(obs) == pa) {
            // populate zbpa matrix with (1, y_1, y_2, ..., y_m) entries for this parental combo
            zbpa.Store (count_n, 0, 1);
            for (long cpar = 0; cpar < continuous_parent_count; cpar++) {
              zbpa.Store (count_n, cpar+1, theData(obs, c_parents.list_data[cpar]));
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
        log_score += BottcherScore (yb, zbpa, tau, mu, rho, phi, n_ij.list_data[pa]);
      }
    } catch (const _String & err) {
      HandleApplicationError(err);
    }

    return log_score;
}



//___________________________________________________________________________________________________

hyFloat  _BayesianGraphicalModel::BottcherScore (_Matrix const& yb, _Matrix const& zbpa, _Matrix const& tau, _Matrix const& mu, hyFloat rho, hyFloat phi, long batch_size) const
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
    try {
      _Matrix temp_mat (zbpa);

      //ReportWarning (_String("zbpa = ") & (_String *) temp_mat.toStr());

      if (tau.check_dimension(1L,1L)) {
        temp_mat *= 1./(tau(0,0)); // just use scalar multiplication of the reciprocal
      } else {
        _Matrix * tauinv = (_Matrix *) tau.Inverse(nil);
        temp_mat *= *tauinv;
        DeleteObject (tauinv);
      }

      //ReportWarning (_String("... %*% tau^(-1) = ") & (_String *) temp_mat.toStr());

      _Matrix t_zbpa (zbpa);
      t_zbpa.Transpose();
      temp_mat *= t_zbpa;

      //ReportWarning (_String("... %*% zbpa^T = ") & (_String *) temp_mat.toStr());

      for (long row = 0; row < temp_mat.GetHDim(); row++) {
          temp_mat.set (row, row) += 1.;
      }

      //ReportWarning (_String("... + I = ") & (_String *) temp_mat.toStr());

      _Matrix scale (temp_mat);
      scale *= phi / rho;

      //ReportWarning (_String("scale = ") & (_String *) scale.toStr());

      // calculate the determinant of scale parameter matrix

      temp_mat = scale;
      temp_mat *= (hyFloat) (pi_const * rho);
    
    //ReportWarning (_String("BottcherScore() calling Eigensystem on matrix ") & (_String *) temp_mat.toStr() );
    
      _AssociativeList *  eigen       = (_AssociativeList *) temp_mat.Eigensystem(nil);
    
    // sometimes the eigendecomposition fails
      if ( eigen->countitems() == 0 ) {
        throw _String("Eigendecomposition failed in bayesgraph2.cpp BottcherScore().");
      }
    
      _Matrix *        eigenvalues = (_Matrix *)eigen->GetByKey(0, MATRIX);
      hyFloat          log_det         = 0.;   // compute determinant on log scale to avoid overflow
      eigenvalues->ForEachCellNumeric([&log_det] (hyFloat cell, long, long, long) -> void {
        log_det += log (cell);
      });
      DeleteObject (eigen);
    
      //ReportWarning (_String("log(det) = ") & log_det);

      // calculate first term of score
      hyFloat  pa_log_score = 0.;

      pa_log_score += _ln_gamma((rho + batch_size) * 0.5);
      pa_log_score -= _ln_gamma(rho * 0.5) + 0.5 * log_det;


      // calculate second term of score
      _Matrix     next_mat (yb),  // N x 1
                  copy_mu (mu);

      
      temp_mat = zbpa;        // N x (k+1)
      temp_mat *= copy_mu;         // (k+1) x 1
      next_mat -= temp_mat;   // (yb - zbpa %*% mu)
      next_mat.Transpose();

      temp_mat = next_mat;
    
      _Matrix * scaleinv = (_Matrix *) scale.Inverse(nil);
      temp_mat *= *scaleinv;  // N x N
      DeleteObject (scaleinv);
    
      next_mat.Transpose();
      temp_mat *= next_mat;
      temp_mat *= (hyFloat) 1./rho;

      //ReportWarning (_String ("2nd term = ") & (hyFloat) temp_mat(0,0));

      pa_log_score += -(rho + batch_size)/2. * log(1. + temp_mat(0,0));


      //ReportWarning (_String("BottcherScore() returning with log L = ") & pa_log_score);

      return pa_log_score;
    } catch (_String const & err) {
      HandleApplicationError(err);
    }
  return -INFINITY;
}

//___________________________________________________________________________________________________

void _check_impute_settings (long& impute_maxsteps, long& impute_burnin, long& impute_samples) {
  static const _String          kHYBgm_IMPUTE_MAXSTEPS    ("BGM_IMPUTE_MAXSTEPS"),
                                kHYBgm_IMPUTE_BURNIN      ("BGM_IMPUTE_BURNIN"),
                                kHYBgm_IMPUTE_SAMPLES     ("BGM_IMPUTE_SAMPLES");

  impute_maxsteps = hy_env :: EnvVariableGetNumber(kHYBgm_IMPUTE_MAXSTEPS, 0.);
  impute_burnin   = hy_env :: EnvVariableGetNumber(kHYBgm_IMPUTE_BURNIN, 0.);
  impute_samples  = hy_env :: EnvVariableGetNumber(kHYBgm_IMPUTE_SAMPLES, 0.);

  if (impute_maxsteps <= 0 || impute_burnin < 0 || impute_samples <= 0 || impute_samples > impute_maxsteps) {
    throw _String("Invalid IMPUTE setting(s)");
  }
}

void _sample_discrete_node (long levels, hyFloat prior_sample_size, _Matrix& observed_values, long row_index) {
// not enough observations - sample from prior given hyperparameters
  
  long obs_total = 0L;
  for (long lev = 0; lev < levels; lev++) {
    obs_total += observed_values (row_index, lev);    // sum tallies across levels
  }
  
  if (obs_total < kBGMMinSize) {
    _Matrix imaginary_sample_counts (levels, 1, false, true);

    prior_sample_size /= levels;
    
    for (long lev = 0; lev < levels; lev++) {
      imaginary_sample_counts.set (lev, 0) = prior_sample_size  + kBGMDirichletFlattening;
    }

    // random deviate from Dirichlet distribution parameterized by imaginary sample size hyperparameters
    _Matrix * result = (_Matrix *) imaginary_sample_counts.DirichletDeviate();

    for (long lev = 0; lev < levels; lev++) {
      observed_values.set (row_index, lev) = result->get (lev, 0);
    }

    DeleteObject (result);
  } else {
    prior_sample_size = 1./obs_total;
    for (long lev = 0; lev < levels; lev++) {
      observed_values.set (levels, lev) *= prior_sample_size;
    }
  }
}

//___________________________________________________________________________________________________
hyFloat _BayesianGraphicalModel::ImputeDiscreteNodeScore (long node_id, _SimpleList const &  parents) const {
  /* ---------------------------------------------------------------------------------------
   ImputeDiscreteNodeScore
   Impute (i.e., "fill in") missing values in the data.  More exactly, take the
   expectation of node score over imputations of the data sampled using Gibbs
   sampling.  Missing values are flagged for discrete data by negative values.
   --------------------------------------------------------------------------------------- */
  
  
  
  
  hyFloat log_score = -INFINITY;
  
  try {
    
    
    long impute_maxsteps = 0L, impute_burnin = 0L, impute_samples = 0L;
    _check_impute_settings (impute_maxsteps, impute_burnin,impute_samples);
    // set Gibbs sampler parameters from batch language definitions
    
    
    
    // number of steps between samples
    
    long sampling_interval    = Maximum(impute_maxsteps / impute_samples, 1L),
    r_i                  = num_levels.get(node_id), // number of levels for child node
    max_num_levels       = r_i,
    family_size          = parents.countitems() + 1L,
    num_parent_combos    = 1;
    
    
    _SimpleList                family_nlevels (family_size, 0, 0),
                               family_index   (family_size, 0, 0),
                               multipliers,
                               is_missing;
    
    _List                      heap_tracker;
    
    family_nlevels[0]         = r_i;
    family_index[0]           = node_id;
    
    
    for (long nlev, p = 0; p < parents.lLength; p++) {
      long parent_idx = parents.get (p);
      family_index [p+1] = parent_idx;
      nlev = family_nlevels[p+1] = num_levels.get(parent_idx);
      num_parent_combos *= nlev;
      multipliers << num_parent_combos;
      if (nlev > max_num_levels)  {
        max_num_levels = nlev;
      }
    }
    
    _Matrix n_ijk (num_parent_combos, r_i, false, true),
    n_ij  (num_parent_combos, 1, false, true),
    data_deep_copy (theData.GetHDim(), family_size, false, true),
    observed_values (family_size, Maximum (max_num_levels, 3L), false, true);
    
    
    // allocate space to matrices
    
    //pa_indices.Populate(theData.GetHDim(), 0, 0);
    
    _Vector * vector_of_scores = new _Vector (),
            * reassign_probs   = new _Vector ();
    
    heap_tracker < vector_of_scores < reassign_probs;
    
    // allocate space to _GrowingVector object -- used for imputing discrete nodes only
    for (long pa = 0L; pa < num_parent_combos; pa++) {
      (*reassign_probs) << 0.;
    }
    
    
    // make deep copy, annotate which entries are missing, and record empirical distributions
    for (long row = 0; row < theData.GetHDim(); row++) {
      // findex = family-wise index
      for (long findex = 0; findex < family_size; findex++) {
        long fnode = family_index.get (findex);  // parent node, get network-wide index
        hyFloat fstate = theData (row, fnode);
        data_deep_copy.set (row, findex) = fstate;
        
        if (fstate < 0L) {
          // store location of missing entry by linear index
          is_missing << row * family_size + findex;
        } else {
          // use state to index into vector and increment counter
          observed_values.set (findex, fstate) += 1.;
        }
      }
    }
    
    
    
    // make a decision whether to use empirical or prior distribution to do initial imputations
    for (long findex = 0L; findex < family_size; findex++) {
       long fnode = family_index.get (findex);
      
      // count non-missing entries across states for this node
      _sample_discrete_node (family_nlevels.get(findex), prior_sample_size (fnode, 0), observed_values, findex);
    }
    
    
    
    // initialize missing entries to random assignments based on observed cases or prior info
    for (long missing_idx = 0; missing_idx < is_missing.countitems(); missing_idx++) {
      long row     = is_missing.get(missing_idx) / family_size;
      long col     = is_missing.get(missing_idx) % family_size;
      //fnode   = (col == 0) ? node_id : parents.list_data[col-1];
      
      
      hyFloat urn = genrand_real2 ();
      
      for (long level = 0; level < family_nlevels.get(col); level++) {
        if (urn <= observed_values (col, level)) {
          data_deep_copy.Store (row, col, (hyFloat) level);
          break;
        } else {
          urn -= observed_values (col, level);
        }
      }
    }
    
    
    // tally N_ijk summary statistics for imputed data
    for (long row = 0; row < data_deep_copy.GetHDim(); row++) {
      long pa_index    = 0;
      long child_state = data_deep_copy (row, 0);
      for (long par = 0; par < parents.countitems(); par++) {
        pa_index += (long) data_deep_copy (row, par+1) * multipliers[par];
      }
      n_ijk.Store (pa_index, child_state, n_ijk(pa_index, child_state) + 1);
      n_ij.Store  (pa_index, 0, n_ij(pa_index, 0) + 1);
    }
    
    
    // Gibbs sampler
    for (long iter = 0; iter < impute_maxsteps + impute_burnin; iter++) {
      // this line was commented out in bgm2.cpp, but I can't remember why - afyp
      //is_missing.Permute(1);  // shuffle list of missing entries
      
      log_score = K2Score (r_i, n_ij, n_ijk); // compute initial score
      
      for (long row, col, pa_index, missing_idx = 0; missing_idx < is_missing.countitems(); missing_idx++) {
        row         = is_missing.list_data[missing_idx] / family_size;
        col         = is_missing.list_data[missing_idx] % family_size;
        pa_index    = 0;
        
        
        // determine parent combination for this row -- use [pa_indices] object instead?  AFYP
        for (long par = 0; par < parents.countitems(); par++) {
          pa_index += data_deep_copy (row, par+1) * multipliers.get(par);
        }
        
        reassign_probs->ZeroUsed(); // reset _GrowingVector
        
        if (col == 0) { // child node, reassign N_ijk's, keep N_ij's constant
          long child_state = data_deep_copy (row, col);
          log_score -= log(n_ijk(pa_index, child_state));
          
          n_ijk.Store (pa_index, (long) child_state, n_ijk(pa_index, (long) child_state) - 1);
          
          // compute probabilities for all instantiations of child node
          for (long k = 0; k < r_i; k++) {
            reassign_probs->Store (log_score + log (n_ijk(pa_index, k) + 1));
          }
          
          hyFloat denom = LogSumExpo (reassign_probs);
          
          // random reassignment
          hyFloat urn = genrand_real2 ();
          
          for (long lev = 0; lev < r_i; lev++) {
            hyFloat this_prob = exp ((* (_Matrix *) reassign_probs) (lev,0) - denom);
            
            if ( urn  <  this_prob ) {
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
          long parent_state    = data_deep_copy (row, col);
          long child_state     = data_deep_copy (row, 0);
          
          log_score   -= log (n_ijk (pa_index, (long) child_state));
          log_score   += log (n_ij  (pa_index, 0) + 1);
          
          n_ij.Store  (pa_index, 0, n_ij(pa_index,0) - 1);
          n_ijk.Store (pa_index, (long) child_state, n_ijk(pa_index, (long) child_state) - 1);
          
          pa_index    -= parent_state * multipliers.list_data[col-1];
          
          
          for (long pa_temp, lev = 0; lev < family_nlevels.get(col); lev++) {
            pa_temp = pa_index + lev * multipliers.get(col-1);
            reassign_probs->Store (log_score + log (n_ijk (pa_temp, (long) child_state) + 1)
                                   - log (n_ij (pa_temp, 0) + 1) ); // why was this adding 2?
          }
          
          
          // compute denominator and convert log-values to probability vector
          hyFloat denom = LogSumExpo (reassign_probs);
          
          // re-assign state
          hyFloat urn = genrand_real2();  // a uniform random number within interval [0,1)
          
          for (long lev = 0; lev < family_nlevels.list_data[col]; lev++) {
            hyFloat this_prob = exp((* (_Matrix *) reassign_probs) (lev, 0) - denom);
            
            if (urn < this_prob) {
              pa_index += lev * multipliers.get(col-1);
              
              n_ij.Store (pa_index, 0, n_ij(pa_index,0) + 1);
              n_ijk.Store (pa_index, (long) child_state, n_ijk(pa_index, (long) child_state) + 1);
              
              log_score = (* (_Matrix *) reassign_probs) (lev, 0);
              
              data_deep_copy.Store (row, col, lev);
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
    
    
    log_score = vector_of_scores->MaxElement(1) / vector_of_scores->get_used(); // mean of sampled scores
    
    
    // compute the average of sampled log-likelihoods
    //log_score = LogSumExpo (vector_of_scores) - log((double)vector_of_scores->GetUsed());
    
  } catch (const _String & err) {
    HandleApplicationError(err);
  }
  
  return (log_score);
}



//___________________________________________________________________________________________________
hyFloat _BayesianGraphicalModel::ImputeCGNodeScore (long node_id, _SimpleList const & parents) const {
  /* ---------------------------------------------------------------------------------------
   ImputeNodeScore ()
   Impute (i.e., "fill in") missing values in the data.  More exactly, take the
   expectation of node score over imputations of the data sampled using Gibbs
   sampling.  Missing values are flagged for discrete data by negative values.
   For continuous data, I haven't settled on an unambiguous flag and I'm leaving this
   up to the user for now (BGM_CONTINUOUS_MISSING_VALUE).
   --------------------------------------------------------------------------------------- */
  
  
  hyFloat      log_score           = 0.;
  try {
    
    long impute_maxsteps = 0L, impute_burnin = 0L, impute_samples = 0L;
    _check_impute_settings (impute_maxsteps, impute_burnin,impute_samples);

    
    long sampling_interval    = Maximum(impute_maxsteps / impute_samples, 1L),
         r_i                  = num_levels.get(node_id), // number of levels for child node
         max_num_levels       = r_i,
         family_size          = parents.countitems() + 1L,
          num_parent_combos    = 1;
    
    _SimpleList                family_nlevels (family_size, 0, 0),
                               family_index   (family_size, 0, 0),
                               multipliers,
                               is_missing,
                               dparents,
                               cparents,
                               parents_by_nodetype;
    
    family_nlevels[0]         = r_i;
    family_index[0]           = node_id;
    _List                      heap_tracker;
    // set Gibbs sampler parameters from batch language definitions

    
    
    
    // partition parent list into discrete and continuous; record stuff while we're at it
    
    for (long p = 0; p < parents.countitems(); p++) {
      long parent_index = parents.get (p);
      if (is_node_discrete(p)) {
        // parent is discrete-valued
        dparents << parent_index;
        parents_by_nodetype << (dparents.countitems()-1);	// stores index into dparents
        
        long nlev = family_nlevels[p+1] = num_levels.get (parent_index);
        
        num_parent_combos *= nlev;
        multipliers << num_parent_combos;
        
        if (nlev > max_num_levels)  {
          max_num_levels = nlev;
        }
      } else {
        // parent is continuous-valued
        cparents << parent_index;
        parents_by_nodetype << (cparents.countitems()-1);
      }
    }
    
    hyFloat rho = prior_sample_size (node_id, 0) > 0. ? (prior_sample_size (node_id, 0) / num_parent_combos) : 1.0,
            phi = prior_scale (node_id, 0);

    // parameters for CG
    _Matrix rho_mx (rho, 1, 1),
            phi_mx (phi, 1, 1);
    
    _Matrix * iw_ptr;
    hyFloat iw_deviate;

    
    // allocate space to matrices
    
    _Matrix n_ijk (num_parent_combos, r_i, false, true),
    n_ij  (num_parent_combos, 1, false, true),
    data_deep_copy (theData.GetHDim(), family_size, false, true),
    observed_values (family_size, Maximum (max_num_levels, 3L), false, true);

    
    //pa_indices.Populate(theData.GetHDim(), 0, 0);
    
    
    // allocate space to _GrowingVector object -- used for imputing discrete nodes only
    _Vector * vector_of_scores = new _Vector (),
            * reassign_probs   = new _Vector ();
    
    heap_tracker < vector_of_scores < reassign_probs;
    
    // allocate space to _GrowingVector object -- used for imputing discrete nodes only
    for (long pa = 0L; pa < num_parent_combos; pa++) {
      (*reassign_probs) << 0.;
    }

    
    
    
    // make deep copy, annotate which entries are missing, and record empirical distributions
    for (long row = 0; row < theData.GetHDim(); row++) {
      // findex = family-wise index
      for (long findex = 0; findex < family_size; findex++) {
        long fnode = family_index.get (findex);  // parent node, get network-wide index
        hyFloat fstate = theData (row, fnode);
        data_deep_copy.set (row, findex) = fstate;
        
        if (is_node_discrete(fnode)) {  // discrete node
          if (fstate < 0) {
            is_missing << row * family_size + findex; // note this is family level indexing
          } else {
            // use state into index into vector and increment counter
            observed_values.set (findex, fstate) += 1.;
          }
        } else {                            // continuous node
          if (fstate == continuous_missing_value) {
            is_missing << row * family_size + findex;
          } else {
            // column 0 stores number of observations
            observed_values.set (findex, 0) += 1.;
            
            // column 1 stores sum for computing mean
            observed_values.set (findex, 1) +=  fstate;
            
            // column 2 stores sum of squares for computing variance
            observed_values.set (findex, 2) += fstate*fstate;
          }
        }
      }
    }
    
    //ReportWarning (_String ("is_missing = ") & (_String *) is_missing.toStr() );

    
    // make a decision whether to use empirical or prior distribution to do initial imputations
    for (long findex = 0L; findex < family_size; findex++) {
      long fnode = family_index.get (findex);
      
      if (is_node_discrete(fnode)) {
        // DISCRETE NODE
        _sample_discrete_node (family_nlevels.get (findex), prior_sample_size (fnode, 0), observed_values, findex);

      } else {
        // CONTINUOUS NODE
        if (observed_values(findex, 0) < kBGMMinSize) {
          
          // unknown variance drawn from scaled inverse chi-square distribution
          _Matrix		temp1 (prior_scale(fnode,0),1,1),
                    temp2 (prior_sample_size(fnode,0),1,1);
          
          _Matrix *	result = (_Matrix *) temp1.InverseWishartDeviate (temp2);
          
          observed_values.set (findex, 2) = (*result)(0,0);
          observed_values.set (findex, 1) = gaussDeviate() * sqrt( (*result)(0,0) ) + prior_mean(fnode,0);
          
          DeleteObject (result);
          
        } else {
          // compute empirical mean and s.d.
          
          // E[x] = sum(x) / len(x)
          observed_values.set (findex, 1) /= observed_values(findex,0);
          // Var[x] = E[x^2] - (E[x])^2
          observed_values.set (findex, 2) = sqrt(observed_values(findex,2) / observed_values(findex,0)
                                               - observed_values(findex,1) * observed_values(findex,1)) ;
        }
      }
    }
    
    // initialize missing entries to random assignments based on observed cases or prior info
    for (long missing_idx = 0; missing_idx < is_missing.countitems(); missing_idx++) {
      long row     = is_missing.get(missing_idx) / family_size;
      long col     = is_missing.get(missing_idx) % family_size;
      long fnode =  family_index.get (col);
      
      if (is_node_discrete(fnode)) {  // discrete node
        hyFloat urn = genrand_real2 ();
        
        for (long level = 0; level < family_nlevels.get(col); level++) {
          if (urn < observed_values (col, level)) {
            data_deep_copy.Store (row, col, (hyFloat) level);
            break;
          } else {
            urn -= observed_values (col, level);
          }
        }
      } else {                            // continuous node
        data_deep_copy.set (row, col) = gaussDeviate() * observed_values(col,2) + observed_values(col,1);
      }
    }
    
    
    long            continuous_parents   = cparents.countitems();
    
    hyFloat      lk_ratio,
    next_log_score;     // log-likelihood of the proposed step
    
    _Matrix         log_scores_by_pa (num_parent_combos, 1, false, true),   // track separately to make updating easier
                    next_log_score_by_pa (num_parent_combos, 1, false, true),
                    mu (continuous_parents + 1, 1, false, true),
                    tau (continuous_parents+1, continuous_parents+1, false, true);
    
    
    
    // set mean intercept/regression coeff. hyperparameter vector for CG prior
    mu.set (0, 0) = prior_mean (node_id, 0);   // intercept
    
    
    // set precision hyperparameter for CG prior
    
    
    tau.set (0, 0) = prior_precision (node_id, 0);
    for (long i = 1; i < continuous_parents+1; i++) {
      tau.set (i, i) = prior_precision(cparents.get(i-1),0);
    }
    
    _SimpleList pa_indices (theData.GetHDim(), 0, 0);
     // partition data set by discrete parent state combination
    for (long obs = 0; obs < theData.GetHDim(); obs++) {
      long pa_index = 0;
      
      for (long findex = 1; findex < family_size; findex++) {
        long this_parent = family_index.get (findex);
        
        if (is_node_discrete(this_parent)) {
          // discrete parent
          pa_index += data_deep_copy (obs, findex) * multipliers.get (parents_by_nodetype.get (findex-1));
        }
      }
      
      pa_indices[obs] = pa_index;
      n_ij.set (pa_index, 0) += 1.;
    }
    
    
    
    // calculate CG summary statistics and compute initial scores
    log_score = 0.;
    
    for (long pa = 0; pa < num_parent_combos; pa++) {
      if (n_ij(pa,0) == 0) {
        continue; // this batch is empty
      }
      
      // these are single-use matrices, reset with every loop
      _Matrix zbpa (n_ij(pa,0), continuous_parents+1, false, true);	// (1, CG parent node values) in batch
      _Matrix yb (n_ij(pa,0), 1, false, true);		// CG child node values in batch
      
      
      // populate zbpa with continuous parent states
      long batch_count = 0;
      for (long obs = 0; obs < pa_indices.countitems(); obs++) {
        if (pa_indices.get(obs) == pa) {
          // (1, x_0, x_1, ..., x_N)
          zbpa.Store (batch_count, 0, 1); // intercept
          
          for (long findex = 1; findex < family_size; findex++) {
            if (is_node_continuous(family_index.get (findex))) {
              zbpa.Store (batch_count, parents_by_nodetype.get(findex-1), data_deep_copy(obs, findex));
            }
          }
          
          yb.set (batch_count, 0) = data_deep_copy(obs, 0);  // child node state
          
          batch_count++;
        }
      }
      
      log_scores_by_pa.Store (pa, 0, BottcherScore (yb, zbpa, tau, mu, rho, phi, n_ij(pa,0)));
      log_score += log_scores_by_pa (pa, 0);
    }
    
    
    
    // Gibbs sampling over missing entries
    for (long iter = 0; iter < impute_burnin + impute_maxsteps; iter++) {
      //is_missing.Permute(1);
      
      for (long missing_idx = 0; missing_idx < is_missing.lLength; missing_idx++) {
        // indices into data_deep_copy for this missing value
        long row         = is_missing.get(missing_idx) / family_size;
        long col         = is_missing.get(missing_idx) % family_size;
        long pa_index    = pa_indices.get(row);
        
        
        // if child node, then row remains in same 'pa' batch, Z^b_pa and N_ij unchanged
        if (col == 0) {
          
          // draw new state from posterior of CG node
          
          iw_ptr = (_Matrix *) phi_mx.InverseWishartDeviate (rho_mx);
          iw_deviate = (*iw_ptr)(0,0);
          //ReportWarning(_String("sigma=") & iw_deviate);
          
          
          _Matrix zbpa (n_ij(pa_index,0), continuous_parents+1, false, true);
          _Matrix yb (n_ij(pa_index,0), 1, false, true);
          
          // collect all continuous parent states under given discrete parent combination (pa_index)
          long batch_count = 0;
          for (long obs = 0; obs < pa_indices.countitems(); obs++) {
            if (obs == row) continue; // skip the current case
            if (pa_indices.get(obs) == pa_index) {
              zbpa.Store (batch_count, 0, 1);
              for (long findex = 1; findex < family_size; findex++) {
                if (is_node_continuous(family_index.get (findex))) {
                  zbpa.Store (batch_count, parents_by_nodetype.get(findex-1), data_deep_copy(obs, findex));
                }
              }
              yb.Store (batch_count, 0, data_deep_copy(obs, 0));
              batch_count++;
            }
          }
          
          //ReportWarning(_String("zbpa=") & (_String *) zbpa.toStr());
          zbpa.Transpose();
          //ReportWarning(_String("zbpa=") & (_String *) zbpa.toStr());
          zbpa *= yb;
          //ReportWarning(_String("zbpa=") & (_String *) zbpa.toStr());
          
          
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
            //data_deep_copy.Store (row, col, child_state);
            // TODO SLKP 20180919 :: THERE WAS NO CHILD_STATE SET
          }
        }
        
        /* ---------- MISSING VALUE AT PARENT NODE (col > 0) ---------- */
        else {
          // if parent is discrete, then this case is moved to another batch (pa) -- MORE THAN ONE BATCH IS AFFECTED
          if (is_node_discrete(family_index.get (col))) {
            
            long parent_state    = data_deep_copy (row, col);
            
            // first compute the likelihood of the original batch minus this case
            n_ij.set  (pa_index, 0) -= 1.;
            
            _Matrix zbpa (n_ij(pa_index,0), continuous_parents+1L, false, true);
            _Matrix yb (n_ij(pa_index,0), 1, false, true);
            
            long batch_count = 0;
            for (long obs = 0; obs < data_deep_copy.GetHDim(); obs++) {
              if (obs == row) {
                continue;    // skip this case
              }
              
              if (pa_indices.list_data[obs] == pa_index) {
                zbpa.Store (batch_count, 0, 1);
                
                for (long findex = 1; findex < family_size; findex++) {
                   if (is_node_continuous(family_index.get (findex))) {
                    zbpa.Store (batch_count, parents_by_nodetype.list_data[findex-1], data_deep_copy(obs, findex));
                  }
                }
                
                yb.Store (batch_count, 0, data_deep_copy(obs, 0));
                
                batch_count++;
              }
            }
            
            next_log_score_by_pa.Store (pa_index, 0, BottcherScore (yb, zbpa, tau, mu, rho, phi, (long) n_ij (pa_index, 0)));
            next_log_score = log_score - log_scores_by_pa (pa_index,0) + next_log_score_by_pa (pa_index,0);
            
            
            // random re-assignment of discrete parent node, compute likelihood of second batch
            
            hyFloat urn = genrand_real2() * (1.0 - observed_values(col, parent_state)); // rescale to omit original state
            
            for (long lev = 0; lev < family_nlevels.list_data[col]; lev++) {
              if (lev == parent_state) {
                continue;
              }
              
              if (urn < observed_values(col, lev)) {
                data_deep_copy.Store (row, col, lev);
                
                /* BREAK */
                
                pa_index += (lev - parent_state) * multipliers.list_data[parents_by_nodetype.list_data[col-1]];
                
                
                // compute summary stats
                n_ij.Store  (pa_index, 0, n_ij(pa_index,0) + 1);
                
                _Matrix zbpa (n_ij(pa_index,0), continuous_parents+1, false, true);
                _Matrix yb (n_ij(pa_index,0), 1, false, true);
                
                batch_count = 0;
                
                for (long obs = 0; obs < data_deep_copy.GetHDim(); obs++) {
                  // we won't update pa_indices[] for this case, already have pa_index to work with
                  if (obs == row || pa_indices.get(obs) == pa_index) {
                    zbpa.Store (batch_count, 0, 1);
                    
                    for (long findex = 1; findex < family_size; findex++) {
                      if (is_node_continuous(family_index.get (findex))) {
                        zbpa.Store (batch_count, parents_by_nodetype.list_data[findex-1], data_deep_copy(obs, findex));
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
              log_scores_by_pa.Store (pa_indices.list_data[row], 0, next_log_score_by_pa (pa_indices.list_data[row],0) ); // contains the old index
              
              pa_indices.list_data[row] = pa_index;   // replace old pa index with new
            } else {
              // revert to previous state
              data_deep_copy.Store (row, col, parent_state);
              
              n_ij.Store (pa_index, 0, n_ij(pa_index,0) - 1);
              n_ij.Store (pa_indices.list_data[row], 0, n_ij(pa_indices.list_data[row],0) + 1);
              // pa_index = pa_indices.list_data[row];
            }
          }
          
          // otherwise, if parent is continuous then just update Z_bpa; case stays in the same batch
          else {
            hyFloat parent_state = data_deep_copy(row,col);
            // random draw from Gaussian distribution centered at previous value
            data_deep_copy.Store (row, col, (gaussDeviate() + data_deep_copy(row,col)) * observed_values(0,2));
            
            _Matrix zbpa (n_ij(pa_index,0), continuous_parents+1, false, true);
            _Matrix yb (n_ij(pa_index,0), 1, false, true);
            
            long batch_count = 0;
            
            for (long obs = 0; obs < data_deep_copy.GetHDim(); obs++) {
              if (pa_indices.get(obs) == pa_index) {
                zbpa.Store (batch_count, 0, 1);
                
                for (long findex = 1; findex < family_size; findex++) {
                  if (is_node_continuous(family_index.get (findex))) {
                    zbpa.Store (batch_count, parents_by_nodetype.list_data[findex-1], data_deep_copy(obs, findex));
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
              pa_indices.list_data[row] = pa_index;
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
    
    
    
    log_score = vector_of_scores->MaxElement(1) / vector_of_scores->get_used();
    // compute the average of sampled log-likelihoods
    //log_score = LogSumExpo (vector_of_scores) - log((double)vector_of_scores->GetUsed());
  
  } catch (const _String & err) {
    HandleApplicationError(err);
    return -INFINITY;
  }
  
  return (log_score);
}








