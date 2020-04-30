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


#include "bayesgraph.h"
#include "function_templates.h"
#include "time_difference.h"
#include "hbl_env.h"
#include "ntuplestorage.h"

#include "global_things.h"

using namespace hy_global;

const   _String  kBGMContraintMatrix ("BGM_CONSTRAINT_MATRIX");

/*
_String

            // SLKP 20070926; add string constants for progress report updates
            _HYBgm_STATUS_LINE_MCMC           ("Running Bgm MCMC"),
            _HYBgm_STATUS_LINE_MCMC_DONE  ("Finished Bgm MCMC"),
            _HYBgm_STATUS_LINE_CACHE      ("Caching Bgm scores"),
            _HYBgm_STATUS_LINE_CACHE_DONE ("Done caching Bgm scores"),
 

            ,
 

            _HYBgm_K2_RESTARTS ("BGM_K2_RESTARTS"),
            _HYBgm_K2_RANDOMIZE ("BGM_K2_RANDOMIZE"),

 
 
  
            _HYBgm_CONTINUOUS_MISSING_VALUE ("BGM_CONTINUOUS_MISSING_VALUE");
*/


//__________________________________________________________________________________________________________
#ifdef      __UNIX__

void        ConsoleBGMStatus (_String const statusLine, hyFloat percentDone, _String const * fileName) {
    FILE           *outFile = fileName?doFileOpen (fileName->get_str(),"w"):nil;
    _String        reportLine (statusLine);

    if (percentDone >= 0.0) {
        reportLine = reportLine & ". " & percentDone & "% done.";
    }

    if (outFile) {
        fprintf (outFile,"%s", reportLine.get_str());
    } else if (verbosity_level == 1) {
        printf ("\033\015 %s", reportLine.get_str());
    }

    if (percentDone < -1.5) {
        printf ("\033\015 ");
        setvbuf (stdout,nil,_IOLBF,1024);
    } else if (percentDone < -0.5) {
        setvbuf (stdout,nil, _IONBF,1);
    }
  
    if (outFile) {
        fclose (outFile);
    }
}

#endif


//__________________________________________________________________________________________________________
hyFloat      _BayesianGraphicalModel::LogSumExpo (_Vector * log_values) {
  
    //  Computes the log of a sum of a vector whose values are stored as log-transforms,
    //  such that exponentiating the vector entries would result in numerical underflow.

    long        size            = log_values->get_used();
    hyFloat     sum_exponents   = 0.;


    // handle some trivial cases
    if (size == 0L) {
        return 0.;    // log(exp(0)) = log(1) = 0
    } else if (size == 1) {
        return log_values->directIndex(0L);    // log(exp(log(x)))
    }

    // find the largest (least negative) log-value
  
    hyFloat      max_log  = (*log_values) (0, 0);
    for (long idx = 1L; idx < size; idx++) {
        StoreIfGreater (max_log, log_values->directIndex(idx));
    }


    // go back through log values and increment by max log value
    //      This will cause some underflow for the smallest values, but we can
    //  use this approximation to handle very large ranges of values.
    //  NOTE: subtracting a negative value.
    for (long idx = 0L; idx < size; idx++) {
        sum_exponents += exp( log_values->directIndex(idx) - max_log );
    }

    return (log(sum_exponents) + max_log);
}




//__________________________________________________________________________________________________________
//__________________________________________________________________________________________________________
//__________________________________________________________________________________________________________
_BayesianGraphicalModel::_BayesianGraphicalModel (_AssociativeList * nodes) {
    
    static const  _String _HYBgm_NODE_INDEX   ("NodeID"),
                          _HYBgm_NODETYPE       ("NodeType"),
                          _HYBgm_NUM_LEVELS ("NumLevels"),
                          _HYBgm_MAX_PARENT ("MaxParents"),
                          _HYBgm_PRIOR_SIZE ("PriorSize"),
                          _HYBgm_PRIOR_MEAN ("PriorMean"),      // for continuous (Gaussian) nodes
                          _HYBgm_PRIOR_PRECISION    ("PriorPrecision"),
                          _HYBgm_PRIOR_SCALE    ("PriorScale");
    
    long                global_max_parents  = 0L;

    node_names.Clear();

    
    try {

    // the fundamental defining characteristic of _BayesianGraphicalModel objects
        num_nodes           = nodes->countitems();


        // allocate space to class matrices
        _Matrix::CreateMatrix (&theStructure, num_nodes, num_nodes, false, true, false);
        _Matrix::CreateMatrix (&constraint_graph, num_nodes, num_nodes, false, true, false);

        _Matrix::CreateMatrix (&prior_sample_size, num_nodes, 1, false, true, false);
        _Matrix::CreateMatrix (&prior_mean, num_nodes, 1, false, true, false);
        _Matrix::CreateMatrix (&prior_precision, num_nodes, 1, false, true, false);
        _Matrix::CreateMatrix (&prior_scale, num_nodes, 1, false, true, false);

 

        // allocate space for _SimpleList objects
        node_type.Populate (num_nodes, 0, 0);
        max_parents.Populate (num_nodes, 0, 0);
        num_levels.Populate(num_nodes, 0, 0);
        has_missing.Populate(num_nodes, 0, 0);


        for (long node = 0L; node < num_nodes; node++) {
            
            const _AssociativeList * this_avl = (_AssociativeList *) (nodes->GetByKey (node, ASSOCIATIVE_LIST));

            // AVL should include following items:  "NodeID"        - variable name
            //                                      "NodeType"      - 0 = discrete, 1 = continuous (Gaussian)
            //                                      "NLevels"       - (discrete only)
            //                                      "MaxParents"    - maximum number of parents
            //                                      "PriorSize"     - hyperparameter for multinomial-Dirichlet
            //                                                      - also used for degrees of freedom hyperparameter for Gaussian node
            //                                      "PriorMean"     - hyperparameter for Gaussian node
            //                                      "PriorPrecision" - hyperparameter for Gaussian node
            //                                      "PriorScale"    - fourth hyperparameter for Gaussian node

            _FString *name = (_FString*)this_avl->GetByKey (_HYBgm_NODE_INDEX, STRING);
            if (!name){
                throw _String("Invalid node name (expected a string) passed to a BGM constructor");
            }
            
            node_names < new _StringBuffer (name->get_str());  // append a pointer to _String duplicate

            // DEBUGGING
            /* wuz: ReportWarning (_String("node_name[") & node & "]=" & (_String *)node_names.list_data[node]);
             20111210 SLKP : _String (_String*) contstructor will actually assume that the argument is
                           : 'unencumbered' i.e. not a member of _Lists etc
                           : this function call will create a stack copy of node_names.list_data[node]
                           : print it to messages.log and then kill the dynamic portion of the object (sData)
                           : this will create all kinds of havoc downstream
             */
            // bug fix:
            ReportWarning (_String("node_name[") & node & "]=" & *((_String *)node_names.GetItem(node)));
            
            // node type (0 = discrete, 1 = continuous)
            node_type[node] = (long)this_avl->GetNumberByKey(_HYBgm_NODETYPE);
            if (node_type.get(node) < 0L || node_type.get(node) > 2L) {
                throw _String("Unsupported NodeType ") & node_type.get(node) & " for node " & node;
            }

            // number of levels (discrete nodes)
            num_levels[node] = (long)this_avl->GetNumberByKey(_HYBgm_NUM_LEVELS);
            if (num_levels.get(node) <= 1) {
                throw _String("NLevels must be greater than 1, received ") & num_levels.get(node) & " for node " & node;
            }
            
            // max parents
            max_parents[node] = (long)this_avl->GetNumberByKey(_HYBgm_MAX_PARENT);
            global_max_parents = Maximum(global_max_parents, max_parents.get (node));
            if (max_parents.get (node) <= 0L || max_parents.get(node) >= num_nodes) {
                throw _String ("MaxParents must be in [1,") & num_nodes & "], received " & max_parents.get (node) & " for node " & node;
            }
            
            prior_sample_size.Store (node, 0, this_avl->GetNumberByKey (_HYBgm_PRIOR_SIZE));
            if (prior_sample_size(node,0) < 0) {
                throw _String ("PriorSampleSize must be at least zero, received ") & prior_sample_size(node,0) & " for node " & node ;
            }
            
            
            
            if (is_node_continuous(node)) { // the next set of options only makes sense for continuous
                // prior mean (Gaussian)
                prior_mean.Store (node, 0,this_avl->GetNumberByKey (_HYBgm_PRIOR_MEAN));

                prior_precision.Store (node, 0, this_avl->GetNumberByKey (_HYBgm_PRIOR_PRECISION));
                if (prior_precision (node, 0) <= 0.) {
                    throw _String ("PriorPrecision must be greater than zero, received ") & prior_precision(node,0) & " for node " & node;
                }

                prior_scale.Store (node, 0, this_avl->GetNumberByKey (_HYBgm_PRIOR_SCALE));
                
                if (prior_scale (node, 0) <= 0.) {
                    throw _String ("PriorScale must be greater than zero, received ") & prior_scale(node,0) & " for node " & node;
                }

                // ReportWarning (_String("prior_precision[") & node & "] set to " & prior_precision(node,0) );
            }
        }


        // allocate node score cache
        
        for (long node = 0L; node < num_nodes; node++) {
            node_score_cache < new _List ((unsigned long)global_max_parents + 1UL);	// appends a dynamic copy of _List object to start
        }
        scores_cached = FALSE;
        ReportWarning (_String ("Constructed BayesianGraphicalModel with ") & num_nodes & " nodes.");
    } catch (const _String & err) {
        HandleApplicationError(err);
    }
}


//__________________________________________________________________________________________________________
_BayesianGraphicalModel::~_BayesianGraphicalModel (void) {
    /* destructor */
}




//__________________________________________________________________________________________________________
bool _BayesianGraphicalModel::SetDataMatrix (_Matrix const * data) {
    
    const static _String _HYBgm_CONTINUOUS_MISSING_VALUE ("BGM_CONTINUOUS_MISSING_VALUE");
    
    /* ------------------------------------------------------------------------------------------------
        SetDataMatrix()
            Takes a matrix pointer passed from HBL and assign it to class member variable.
            Checks that number of levels for discrete nodes agrees with value set at construction.
            Checks for missing values in each column and annotate the [has_missing] list accordingly.
            Missing values for discrete nodes are flagged by any negative integer.
            Missing values for continuous nodes are flagged by a value set by the user; otherwise,
            it defaults to a fun value I picked arbitrarily.
       ------------------------------------------------------------------------------------------------ */


    _SimpleList data_nlevels;

	/* reset missing value indicators to 0, in case we 
		are replacing a data matrix with missing values */
    
	for (long node = 0L; node < num_nodes; node++) {
		has_missing[node] = 0L;
	}
	
    /* check for user assignment of continuous missing value indicator */
    
    continuous_missing_value = hy_env::EnvVariableGetNumber(_HYBgm_CONTINUOUS_MISSING_VALUE, -666.0);
    ReportWarning (_String ("Entered SetDataMatrix() with missing CG flag: ") & continuous_missing_value & " and node types" & (_String *) node_type.toStr());
	
    try {
        if (data->GetVDim() != num_nodes) {
            throw _String("Number of variables in data (") & data->GetVDim() & ") does not match number of nodes in graph (" & num_nodes & ")";
        }

        data_nlevels.Populate (num_nodes, 1, 0);

            // if (theData) DeleteObject (theData);
            // purge duplicate from memory (makeDynamic)

        theData.Duplicate(data);        // duplicate matrix to member variable
        theData.CheckIfSparseEnough (TRUE); // make dense

        scores_cached = FALSE;


        for (long nrows = theData.GetHDim(), node = 0; node < num_nodes; node++) {
            // if discrete node, compute number of levels and missingness
            // missingness is encoded by negative values
            // otherwise the range is assumed to be 0 : max_value (I think (?) SLKP)
            if (is_node_discrete(node)) {
                data_nlevels[node] = 1;

                for (long row = 0L; row < nrows; row++) {
                    long val = theData (row, node);
                    
                    if (val < 0) {
                        has_missing[node] = 1L;
                    } else {
                        if (val + 1L > data_nlevels.get (node)) {
                            ++data_nlevels[node];
                        }
                    }
                }

                if (data_nlevels.get (node) != num_levels.get (node)) {
                    throw (_String (" Number of levels in data (") & data_nlevels.get(node) & ") for discrete node "
                               & node & " is not compatible with node setting (" & num_levels.get(node)
                               & ").  Check your data or reset the BayesianGraphicalModel.");

                }
            } else if (is_node_continuous(node)) {
                for (long row = 0L; row < nrows; row++) {
                    hyFloat val = theData (row, node);
                    if (val == continuous_missing_value) {
                        has_missing[node] = 1;
                        ReportWarning (_String("Detected missing continuous value at row ") & row);
                        break; // TODO SLKP : discard partial row?
                    }
                }
            }
        }
        ReportWarning (_String ("Set data matrix to:\n") & _String((_String *)theData.toStr()) & "\n and missing values at " & _String((_String *) has_missing.toStr()));
        
    } catch (const _String & err) {
        HandleApplicationError (err);
        return false;
    }

    // compute node scores and store in cache
    CacheNodeScores();
    return true;
}

//__________________________________________________________________________________________________________

bool _BayesianGraphicalModel::SetWeightMatrix (_Matrix const * weights) {
    if (weights->GetHDim() == theData.GetHDim() && weights->GetHDim() == num_nodes) {
        theWeights.Duplicate (weights);
        ReportWarning(_String("Assigned weight matrix:\n") & *(_String *) theWeights.toStr());
        return true;
    }

    HandleApplicationError ("Incompatible matrix dimensions in SetWeightMatrix().");
    return false;
}


//__________________________________________________________________________________________________________
bool _BayesianGraphicalModel::SetConstraints (_Matrix const * constraints)
{
    /* -----------------------------------------------------------------
        SetConstraints()
            Assign pointer to _Matrix object passed from batchlan.cpp
            A constraint matrix is N x N where N is the number of nodes
            in the network.
            1 entry indicates an enforced edge (always in network)
            -1 entry indicates a banned edge (never in network)
            0 entry indicates no constraint
       ----------------------------------------------------------------- */
    if (constraints->check_dimension(num_nodes, num_nodes)) {
        constraint_graph.Duplicate(constraints);
        ReportWarning (_String("Assigned constraint matrix:\n ") & *(_String*) constraint_graph.toStr() );
        return true;
    }

    HandleApplicationError ("ERROR: Constraint matrix incompatible dimensions to graph.");
    return (FALSE);
}



//__________________________________________________________________________________________________________

bool _BayesianGraphicalModel::SetStructure (_Matrix const * structure) {
    /* -------------------------------------------------------------------
        SetStructure()
            Assign pointer to _Matrix object from batchlan:SetParameter()
            that specificies a network structure.
            Check the structure against constraint matrix and node
            ordering if the latter is set.
            If node order is incompatible, reset the node order.
       ------------------------------------------------------------------- */

    try {
    
        if (!structure->check_dimension(num_nodes, num_nodes)) {
            throw _String("Structure incompatible dimensions to graph.");
        }
        
        structure->ForEachCellNumeric ([this] (hyFloat value, long index, long row, long col) -> void {
            if (value == 1.) {
                if (this->constraint_graph (row, col) < 0.) {
                     throw (_String ("Structure contains banned edge: ") & _String (row) & _String ("->") & _String (col));
                }
            } else {
                if (value == 0.) {
                    if (this->constraint_graph (row, col) > 0.) {
                        throw (_String ("Structure lacks enforced edge: ") & _String (row) & _String ("->") & _String (col));
                    }
                }
            }
        });

        // is new structure compatible with node order set in HBL?
        if (node_order_arg.lLength == num_nodes && !GraphObeysOrder (theStructure, node_order_arg)) {
            // need to reset node_order
            node_order_arg = GetOrderFromGraph (theStructure);
            ReportWarning ("Structure is incompatible with existing node order, resetting order.");
        }

        theStructure.Duplicate(structure);
        return true;

    } catch (const _String & err) {
        HandleApplicationError (err);
        return false;
    }
    return true;
}

//__________________________________________________________________________________________________________

void _BayesianGraphicalModel::GetStructure (_Matrix * graph)  const{
    graph->Duplicate (&theStructure);
    ReportWarning (_String("GetStructure() copied graph ") & *(_String *) graph->toStr());
}

//__________________________________________________________________________________________________________
bool _BayesianGraphicalModel::SetNodeOrder (_SimpleList const * order) {
    /* -----------------------------------------------------------------
        SetNodeOrder()
            Set node ordering to vector argument from HBL SetParameter().
       ----------------------------------------------------------------- */

    try {
        if (order->countitems() != num_nodes) {
            throw _String("Node order argument has incorrect length.");
        }
        if (!GraphObeysOrder (theStructure, *order)) {
            throw _String("Node order incompatible with current graph.");
        }
        
        node_order_arg.Duplicate(order);
        ReportWarning (_String("BayesianGraphicalModel node order arg set to ") & *(_String *) node_order_arg.toStr());

    } catch (_String const & err) {
        HandleApplicationError(err);
        return false;
    }
    return true;
}

//__________________________________________________________________________________________________________

void _BayesianGraphicalModel::GetNodeOrder (_Matrix * receptacle) const {
    if (node_order_arg.countitems() == num_nodes) {
        node_order_arg.Each ([receptacle] (long value, unsigned long index) -> void {
            receptacle->Store (0, index, value);
        });
    }
}


//__________________________________________________________________________________________________________
hyFloat _BayesianGraphicalModel::Compute (void) {
    /* --------------------------------------------------------------
        Compute()   [POLYMORPHIC]
            Return posterior probability of [CURRENT] network
            structure.
       -------------------------------------------------------------- */

    hyFloat  log_score = 0.;

    for (long node_id = 0L; node_id < num_nodes; node_id++) {
        log_score += is_node_continuous(node_id) ? ComputeContinuousScore (node_id) : ComputeDiscreteScore (node_id);
    }

    return log_score;
}



//__________________________________________________________________________________________________________
hyFloat _BayesianGraphicalModel::Compute (_Matrix const & g) {
    /* --------------------------------------------------------------
        Compute()   [POLYMORPHIC]
            Return posterior probability of [GIVEN] network
            structure argument.
       -------------------------------------------------------------- */

    hyFloat  log_score = 0.;

    for (long node_id = 0; node_id < num_nodes; node_id++) {    // discrete nodes are 0 (FALSE)
        log_score += is_node_continuous(node_id)? ComputeContinuousScore (node_id, g) : ComputeDiscreteScore (node_id, g);
    }

    return log_score;
}



//__________________________________________________________________________________________________________
hyFloat  _BayesianGraphicalModel::Compute (_SimpleList const & node_order, _List * marginals) {
    /* --------------------------------------------------------------
        Compute()   [POLYMORPHIC]
            Return posterior probability of given node order by
            integrating over all compatible structures.
            Also return marginal posterior probabilities for edges.
       -------------------------------------------------------------- */

    hyFloat          log_likel   = 0.;


    // reset _GrowingVector objects stored in _List object
    for (long i = 0; i < num_nodes * num_nodes; i++) {
        ((_Vector *) marginals->GetItem(i))->ZeroUsed();
    }


    for (long nodeIndex = 0L; nodeIndex < node_order.countitems(); nodeIndex++) {
        
        long                child_node      = node_order.get(nodeIndex),
                            maxp            = max_parents.get(child_node);

        _List           *   score_lists     = (_List *) node_score_cache.GetItem(child_node);
        _Constant       *   orphan_score    = (_Constant *) (score_lists->GetItem(0));


        _Vector * gv1 = (_Vector *) marginals->GetItem (child_node * num_nodes + child_node); // store denominator in diagonal
        gv1 -> ZeroUsed();
        gv1 -> Store (orphan_score->Value());   // append score - note gv1 does not change within



        if (maxp > 0L) {
            // all nodes to the right are potential parents, except banned parents!
            // SLKP. To the right? Not the left?
            
            _SimpleList precedes = node_order.Filter([child_node,this] (long value, unsigned long) -> bool {
                return this->constraint_graph (value, child_node) >= 0.;
            }, nodeIndex+1UL);
            

            // handle trivial case of one parent
            _Matrix *   single_parent_scores    = (_Matrix *) (score_lists->GetItem(1));
            precedes.Each ([this, gv1, single_parent_scores,marginals,child_node] (long parent_index, unsigned long) -> void {
                hyFloat score = (*single_parent_scores) (parent_index, 0);
                gv1->Store (score);
                ((_Vector *) marginals->GetItem(child_node * this->num_nodes + parent_index))->Store (score);
            });


            // more than one parent requires k-tuples
            if (maxp > 1) {
                _SimpleList         indices (precedes.countitems(), 0, 1);   // populates list with 0, 1, 2, ..., M-1
                // where M is the number of eligible parents in [precedes]

                for (long nparents = 2L; nparents <= maxp && nparents <= precedes.countitems(); nparents++) {
                    _SimpleList     subset,
                                    auxil;

                    bool            not_finished;

            
                    if (indices.NChooseKInit (auxil, subset, nparents, false)) {
                        hyFloat      tuple_score;
                        _SimpleList     parents;

                        parents.Populate (nparents, 0, 0);  // allocate memory
                        _NTupleStorage * family_scores = (_NTupleStorage *) (score_lists->GetItem(nparents));

                        do {
                            //parents.Clear();
                            not_finished = indices.NChooseK (auxil, subset);    // cycle through index combinations

                            for (long i = 0L; i < nparents; i++) {   // convert indices to parent IDs (skipping child)
                                long    realized = precedes.get(subset.get(i));
                                if (realized >= child_node) {
                                    realized--;
                                }
                                parents[i] = realized;
                            }
                            parents.Sort();
                            tuple_score = family_scores -> Retrieve (parents);

                            gv1 -> Store (tuple_score); // append to denominator

                            for (long i = 0; i < nparents; i++) {
								// update parent combination-specific numerator
                                ((_Vector *) marginals->GetItem(child_node * num_nodes + precedes.get(subset.get(i))))-> Store (tuple_score);
                            }
                        } while (not_finished);
                    }
                }
            }
        }

        hyFloat log_sum = LogSumExpo(gv1);
        gv1 -> _Matrix::Store (0, 0, log_sum);  // replace first entry with sum, i.e. marginal log-likelihood of child node
        log_likel += log_sum;

    }
    // end loop over child nodes

    return log_likel;
}



//__________________________________________________________________________________________________________
hyFloat  _BayesianGraphicalModel::ComputeDiscreteScore (long node_id) {
    return ComputeDiscreteScore (node_id, theStructure);
}

//__________________________________________________________________________________________________________
hyFloat  _BayesianGraphicalModel::ComputeDiscreteScore (long node_id, _Matrix const & g) {
    
    _SimpleList     parents;

    // SLKP TODO 20180902 : do we want to potentially include the node itself?
    for (long par = 0L; par < num_nodes; par++) {
        if ( g(par, node_id) == 1. && is_node_discrete(par)) {
            parents << par;
        }
    }

    return ComputeDiscreteScore (node_id, parents);
}


//__________________________________________________________________________________________________________
hyFloat  _BayesianGraphicalModel::ComputeDiscreteScore (long node_id, _SimpleList const & parents) {
    /* --------------------------------------------------------------------
        ComputeDiscreteScore()
            Returns posterior probability of local network structure
            centered at discrete child node.
            Only discrete nodes may be parents of a discrete child node.
       -------------------------------------------------------------------- */


    // use cached node scores if available
    if (scores_cached) {
        _List *     scores  = (_List *) node_score_cache.GetItem (node_id);

        if (parents.empty()) {
            return ((_Constant *)scores->GetItem(0))->Value();
        } else if (parents.countitems() == 1UL) {
            return (* ((_Matrix *) scores->GetItem(1))) (parents.get(0), 0);
        } else {
            _NTupleStorage *    family_scores   = (_NTupleStorage *) scores->GetItem (parents.countitems());
        
            return (hyFloat) family_scores->Retrieve (parents.MapList( [node_id] (long value, unsigned long) -> long {
                return value > node_id ? value - 1L : value;
            }));
            // using nk-tuple
        }
    }


    //ReportWarning (_String ("Non-cached call of ComputeDiscreteScore with ") & node_id & " <- " & (_String *) parents.toStr());

    // impute score if missing data
    if (has_missing.get(node_id)) {
        return ImputeDiscreteNodeScore (node_id, parents);
    } else {
        if (parents.Any ([this] (long value, unsigned long) -> bool {
            return this->has_missing.get (value);
        })) {
            return ImputeDiscreteNodeScore (node_id, parents);
        }
    }

    _Matrix         n_ijk,
                    n_ij;
    // [i] indexes child nodes,
    // [j] indexes combinations of values for parents of i-th node,
    // [k] indexes values of i-th node

    UpdateDirichletHyperparameters (node_id, parents, &n_ij, &n_ijk);

    return ( (prior_sample_size (node_id, 0) == 0.) ?
             K2Score (node_id, n_ij, n_ijk) :
             BDeScore (node_id, n_ij, n_ijk) );
}


//_______________________________________________________________
void    _BayesianGraphicalModel::UpdateDirichletHyperparameters (long dnode, _SimpleList const &dparents, _Matrix * n_ij, _Matrix * n_ijk) {
    if (! is_node_discrete(dnode)) {
        ReportWarning ("ERROR: UpdateDirichletHyperparameters() called on non-discrete node!  That sucks!");
        return;
    }

    if (dparents.nonempty ()) {
        _SimpleList     multipliers (ComputeParentMultipliers(dparents));
        long            num_parent_combos = multipliers.GetElement(-1L);
 
        _Matrix::CreateMatrix (n_ij, num_parent_combos, 1, false, true, false);
        _Matrix::CreateMatrix (n_ijk, num_parent_combos, num_levels.list_data[dnode], false, true, false);

        hyFloat norm  = prior_sample_size(dnode,0)/num_parent_combos,
                norm2 = 1./num_levels.get(dnode);

        // initialize matrices with prior hyperparameters (a_ijk) -- if using K2, entries will remain zero
        for (long j = 0; j < num_parent_combos; j++) {
            n_ij->Store (j, 0, norm);

            for (long k = 0; k < num_levels.list_data[dnode]; k++) {
                n_ijk->Store (j, k, (*n_ij)(j,0) * norm2);
            }
        }


        // update hyperparameters with data
        /*
            What should we do with incomplete cases?  Currently using Gibbs sampling to impute
                an average node score - should we export one of these samples?
            For now, just export complete cases - BUT THIS WILL CAUSE PROBLEMS IF NONE OF
                THE CASES IS COMPLETE! - afyp, 2010/02/25
         */
        for (long obs = 0L; obs < theData.GetHDim(); obs++) {
            long    index           = 0,
                    child_state     = theData(obs, dnode);

            if (child_state < 0L) {
                continue;    /* missing observation */
            }

            for (long par = 0L; par < dparents.countitems(); par++) {
                long    this_parent         = dparents.get(par),
                        this_parent_state = theData(obs, this_parent);

                if (this_parent_state < 0L) {
                    index = -1; /* missing observation */
                    break;
                }
                index += this_parent_state * multipliers.get(par);
            }

            if (index >= 0) {
				// this is where we would modify the increment value (1) if we are weighting cases...
                n_ijk->Store ((long) index, child_state, (*n_ijk)(index, child_state) + 1);
                n_ij->Store ((long) index, 0, (*n_ij)(index, 0) + 1);
            }
        }
    } else {
        // conditionally independent node
        _Matrix::CreateMatrix (n_ij, 1, 1, false, true, false);
        _Matrix::CreateMatrix (n_ijk, 1, num_levels.get(dnode), false, true, false);

        // prior hyperparameters (uniform distribution)
        hyFloat flat_value = prior_sample_size(dnode,0) / num_levels.get(dnode);
        for (long k = 0; k < num_levels.get(dnode); k++) {
            n_ijk->Store(0, k, flat_value);
        }

        // update with data
        for (long obs = 0; obs < theData.GetHDim(); obs++) {
            long child_state = theData(obs, dnode);

            if (child_state < 0L) {
                continue;
            }

			// this is where we would modify the increment value (1) if we are weighting cases...
            n_ijk->Store (0, child_state, (*n_ijk)(0, child_state) + 1);
            n_ij->Store  (0, 0, (*n_ij)(0,0) + 1);
        }
    }
}



//___________________________________________________________________________________________
hyFloat _BayesianGraphicalModel::K2Score (long node_id, _Matrix const & n_ij, _Matrix const & n_ijk) const {
    hyFloat  log_score   = 0.;
    long        r_i         = num_levels.get(node_id);

    for (long j = 0L; j < n_ij.GetHDim(); j++) {
        log_score += _ln_gamma(r_i);  // (r-1)!
        log_score -= _ln_gamma(n_ij(j, 0) + r_i); // (N+r-1)!

        for (long k = 0L; k < r_i; k++) {
            log_score += _ln_gamma(n_ijk(j,k) + 1);   // (N_ijk)!
        }
    }

    return log_score;
}



//___________________________________________________________________________________________
hyFloat _BayesianGraphicalModel::BDeScore (long node_id, _Matrix const & n_ij, _Matrix const & n_ijk) const {
    // note that n_ij and n_ijk already contain prior counts, updated by data
    hyFloat  n_prior_ij      = prior_sample_size (node_id, 0) / n_ij.GetHDim(),
                n_prior_ijk     = n_prior_ij / num_levels.list_data[node_id],
                log_score        = 0.;

    for (long j = 0L; j < n_ij.GetHDim(); j++) {
        log_score += _ln_gamma(n_prior_ij) - _ln_gamma(n_ij(j,0));

        for (long k = 0L; k < num_levels.get(node_id); k++) {
            log_score += _ln_gamma(n_ijk(j,k)) - _ln_gamma(n_prior_ijk);
        }
    }

    return log_score;
}



//___________________________________________________________________________________________
void    _BayesianGraphicalModel::InitMarginalVectors (_List * compute_list) const {
    /* ------------------------------------------------------------------------------------
        InitMarginalVectors()
            Allocate storage of node and edge marginal posterior probabilities accumulated
            during order-MCMC.  Off-diagonals correspond to entries of an adjacency matrix
            (edges), whereas diagonal entries are used to store node marginals.
       ------------------------------------------------------------------------------------ */

    for (long i = 0; i < num_nodes * num_nodes; i++) {
        (*compute_list) < new _Vector;
    }
}



//___________________________________________________________________________________________
void    _BayesianGraphicalModel::DumpMarginalVectors (_List * compute_list)  const {
    DeleteObject (compute_list);
}

//___________________________________________________________________________________________
void    _BayesianGraphicalModel::CacheNodeScores (void) {
    /*  -----------------------------------------------------------------------------------
        CacheNodeScores() loops through nodes in the network and calls compute functions
        to calculate local network scores given the number and identity of parents for that node.

        The node scores are cached in a custom data structure that comprises a _List object
        that stores pointers to the following objects:

            0 parents   _Constant       a single float for node score without parents
            1           _Matrix         a vector where i-th entry holds score with parent i
            2 or more   _NTupleStorage  a vector indexed by combinadics, i.e., a (n,k) tuple
                                        mapped to integer space

        Each of these _List object are in turn stored in a _List object that is a class member
        variable [node_score_cache].
        ----------------------------------------------------------------------------------- */

    //ReportWarning (_String ("Entered CacheNodeScores()"));
  
    static const _String kHYBgm_MPI_CACHING ("USE_MPI_CACHING");

    if (scores_cached) {
        return;
    }
  
    try {

      _SimpleList     all_but_one (num_nodes-1L, 0, 1),
                      nk_tuple;
      
      _Matrix         single_parent_scores (num_nodes, 1, false, true);
      
      auto            handle_node_score = [this, &single_parent_scores,&all_but_one,&nk_tuple] (long node_id) -> _List * {
                          _SimpleList parents (0L),
                                      aux_list;
        
                          _List*      list_of_matrices = new _List;;
                          long        maxp = max_parents.get(node_id);
                          hyFloat     score;
        
                            // compute single parent scores
                          for (long par = 0L; par < num_nodes; par++) {
                            if (par == node_id) {   // child cannot be its own parent
                              continue;
                            }
                            
                            parents[0UL] = par;
                            
                            if (is_node_discrete(node_id) ) {
                                // discrete-valued child node cannot have continuous parent
                              if (is_node_discrete (par)) {
                                score = ComputeDiscreteScore (node_id, parents);
                              }
                            } else {
                                // continuous child can have discrete or continuous parents
                              score = ComputeContinuousScore (node_id, parents);
                            }
                            
                            single_parent_scores.Store (par, 0UL, score);
                          }
                            // compute multiple parent scores
                          for (long np = 2; np <= maxp; np++) {
                            _NTupleStorage * family_scores = new _NTupleStorage (num_nodes-1, np);
                            
                            parents << 0L;   // allocate space for one additional parent
                            
                            if (all_but_one.NChooseKInit (aux_list, nk_tuple, np, false)) {
                              bool    remaining, family_is_discrete;
                              
                              do {
                                remaining = all_but_one.NChooseK (aux_list, nk_tuple);
                                
                                if (is_node_discrete(node_id)) {
                                  family_is_discrete = TRUE;
                                  for (long par, par_idx = 0; par_idx < np; par_idx++) {
                                    par = nk_tuple.get(par_idx);
                                    if (par >= node_id) {
                                      par++;
                                    }
                                    
                                    if (is_node_discrete (par) == false) {    // all parents must be discrete
                                      family_is_discrete = false;
                                      break;
                                    }
                                    
                                    parents[par_idx] = par;
                                  }
                                  
                                  if (family_is_discrete) {
                                    score = ComputeDiscreteScore (node_id, parents);
                                  }
                                } else {
                                  for (long par_idx = 0; par_idx < np; par_idx++) {
                                    long par = nk_tuple.list_data[par_idx];
                                    if (par >= node_id) {
                                      par++;
                                    }
                                    parents.list_data[par_idx] = par;
                                  }
                                  
                                  score = ComputeContinuousScore (node_id, parents);
                                }
                                
                                family_scores->Store (score, nk_tuple);
                              } while (remaining);
                            } else {
                              throw _String("Failed to initialize _NTupleStorage object in Bgm::CacheNodeScores().\n");
                            }
                            
                            list_of_matrices->AppendNewInstance(family_scores);   // append duplicate to storage
                          }
                          return list_of_matrices;
        
                      };

      auto compute_singleton_scores = [this] (long node_id) -> hyFloat {
        _SimpleList  parents;
        if (is_node_discrete(node_id)) {
          return ComputeDiscreteScore (node_id, parents);    // [_SimpleList parents] should always be empty for master node
        } else {
          return ComputeContinuousScore (node_id, parents);
        }
      };

#if defined __HYPHYMPI__
    bool  use_mpi_caching = hy_env::EnvVariableTrue(kHYBgm_MPI_CACHING);
 
    if (use_mpi_caching) {
        ReportWarning (_String ("Using MPI to cache node scores."));

        // MPI_Init() is called in main()
        int             size, rank;
        long            mpi_node;
      

        MPI_Status      status; // contains source, tag, and error code


        MPI_Comm_size (MPI_COMM_WORLD, &size);  // number of processes
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);


        if (rank == 0) {
            _String  const   bgmSwitch ("_BGM_SWITCH_");
            _StringBuffer    bgmStr (2048UL);

            _List       * this_list;

            _Matrix     * mpi_node_status = new _Matrix ((long)size, 2, false, true);   // semaphore


            // prepare message exporting _BayesianGraphicalModel object to compute nodes
            SerializeBGMtoMPI (bgmStr);
            ReportWarning (_String("Serialized _BayesianGraphicalModel object as:\n") & bgmStr);


            // switch compute nodes from mpiNormal to bgm loop context
            for (long ni = 1L; ni < size; ni++) {
                MPISendString (bgmSwitch, ni);
            }


            // receive confirmation of successful switch
            for (long ni = 1L; ni < size; ni++) {
                long fromNode = ni;

                _String t (MPIRecvString (ni,fromNode));
              
                if (t != bgmSwitch) {
                    throw _String ("Failed to confirm MPI mode switch at node ") & ni;
                    return;
                } else {
                    ReportWarning (_String("Successful mode switch to Bgm MPI confirmed from node ") & ni);
                    MPISendString (bgmStr, ni);
                }
            }


            // farm out jobs to idle compute nodes until none are left
            for (long node_id = 0; node_id < num_nodes; node_id++) { // node_id iterates over network nodes
                long        maxp            = max_parents.get(node_id);

                _String     mxString,
                            mxName;
              
                hyFloat     score = compute_singleton_scores (node_id);



                //_Constant   orphan_score (score);


                this_list   = (_List *) node_score_cache.GetItem(node_id);
                this_list->Clear();
                 // handle orphan score locally
              
              

                if (maxp > 0L) { // don't bother to farm out trivial cases
                    // look for idle nodes
                    mpi_node = 1;
                    do {
                        if ((*mpi_node_status)(mpi_node, 0) == 0) {
                            ReportMPIError(MPI_Send(&node_id, 1, MPI_LONG, mpi_node, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD), true);
                            ReportWarning (_String ("Sent child ") & node_id & " to node " & mpi_node);

                            mpi_node_status->Store (mpi_node, 0, 1);    // set busy signal
                            mpi_node_status->Store (mpi_node, 1, node_id);

                            break;
                        }
                     } while (++mpi_node < size);
                }

                if (mpi_node == size) { // all nodes are busy, wait for one to send results
                    MPIReceiveScores (mpi_node_status, true, node_id);
                }
            }
            // end loop over child nodes


            // collect remaining jobs
            mpi_node = 0L;
            mpi_node_status->ForEachCellNumeric([&mpi_node] (hyFloat value, long, long, long column)-> void {
              if (column == 0L) {
                mpi_node += value > 0.0 ? 1 : 0;
              }
            });
          
            while (mpi_node > 0L) {
              MPIReceiveScores (mpi_node_status, false, 0);
              mpi_node --;
            }
          

            // shut down compute nodes
            for (long shutdown = -1, mpi_node = 1L; mpi_node < size; mpi_node++) {
                ReportMPIError(MPI_Send(&shutdown, 1, MPI_LONG, mpi_node, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD), true);
                ReportWarning (_String ("Node 0 sending shutdown signal to node ") & mpi_node);
            }
        } else {  // compute node
 
            while (1) { // wait for the next message from the master
                _String     mxString,
                            mxName;

                long        node_id;  // update in while loop
              
                // wait for master node to issue node ID to compute
                ReportMPIError (MPI_Recv (&node_id, 1, MPI_LONG, 0, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD, &status), false);
                ReportWarning (_String("Node") & (long)rank & " received child " & (long)node_id & " from node " & (long)status.MPI_SOURCE);

                if (node_id < 0) {
                    ReportWarning (_String("Node") & (long)rank & " recognizes BGM loop shutdown signal.\n");
                    break;  // received shutdown message (-1)
                }
              
              
                _List * list_of_matrices = handle_node_score (node_id);
                _List   memory_cleanup;
                memory_cleanup.AppendNewInstance (list_of_matrices);
                ReportMPIError (MPI_Send (single_parent_scores.theData, num_nodes, MPI_DOUBLE, 0, HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD), true);
                list_of_matrices->ForEach ([&] (BaseRef mx, unsigned long idx) -> void {
                    _Matrix * storedMx = (_Matrix *) mx;
                    ReportMPIError (MPI_Send (storedMx->theData, storedMx->GetHDim(), MPI_DOUBLE, 0, HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD), true);
                });
              
            }
            // end while, received shutdown from master node
        }
        // end else if compute node
    } else {  // perform single-threaded
#else

    for (long node_id = 0L; node_id < num_nodes; node_id++) {
         _List   *   this_list   = (_List *) node_score_cache.GetItem (node_id);    // retrieve pointer to list of scores for this child node

        this_list->Clear();     // reset the list

        // prepare some containers
        hyFloat  score = compute_singleton_scores (node_id);
        this_list->AppendNewInstance(new _Constant (score));
      
      
      

        if (max_parents.get(node_id) > 0) {
            _List *     scores = handle_node_score (node_id); // this will also populate single_parent_scores
            (*this_list) && & single_parent_scores;
            (*this_list) << *scores;
            delete scores;
        }


    }
#endif

#if defined __HYPHYMPI__
    }   // completes else statement
#endif

    } catch (const _String & e) {
      HandleApplicationError(e);
      return;
    }
      
      
    scores_cached = TRUE;
    ReportWarning (_String ("Cached node scores."));
      
}



//________________________________________________________________________________________________________
#if defined __HYPHYMPI__
/*
	Head node processes results from compute node.
 */
void _BayesianGraphicalModel::MPIReceiveScores (_Matrix * mpi_node_status, bool sendNextJob, long node_id)
{
    _Matrix     single_parent_scores (num_nodes, 1, false, true);
    MPI_Status  status;

    ReportMPIError (MPI_Recv (single_parent_scores.theData, num_nodes, MPI_DOUBLE, MPI_ANY_SOURCE,
                              HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD, &status), false);


    long        senderID    = (long) status.MPI_SOURCE,
                this_node    = (long) (*mpi_node_status) (senderID, 1),
                maxp       = max_parents.list_data[this_node];

    _List   *   this_list   = (_List *) node_score_cache.list_data[this_node];


    _String     mxString,
                mxName;


    mpi_node_status->Store (senderID, 0, 0);    // set node status to idle
    ReportWarning (_String("Received scores for child ") & this_node & " from node " & senderID);

    (*this_list) && (&single_parent_scores);


    _SimpleList     parents,
                    all_but_one (num_nodes-1, 0, 1),
                    aux_list,
                    nk_tuple;

    hyFloat      score;

    // parse nk-tuple results
    for (long np = 2; np <= maxp; np++) {
        _NTupleStorage  family_scores (num_nodes-1, np);
        bool            remaining;


        if (all_but_one.NChooseKInit (aux_list, nk_tuple, np, false)) {
            long        score_index     = 0,
                        num_ktuples      = exp(_ln_gamma(num_nodes) - _ln_gamma(num_nodes - np) - _ln_gamma(np+1)),
                        ntuple_receipt;

            _Matrix     scores_to_store (num_ktuples, 1, false, true);

            // receive nk-tuple indexed node scores from same compute node
            ReportMPIError (MPI_Recv (scores_to_store.theData, num_ktuples, MPI_DOUBLE, senderID,
                                      HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD, &status), false);

            do {
                remaining = all_but_one.NChooseK (aux_list, nk_tuple);  // update nk-tuple in aux_list
                ntuple_receipt = family_scores.Store (scores_to_store(score_index, 0), nk_tuple);
                score_index++;
            } while (remaining);
        } else {
            HandleApplicationError("Failed to initialize _NTupleStorage object in Bgm::CacheNodeScores().\n");
        }

        (*this_list) && (&family_scores);
    }


    // send next job
    if (sendNextJob) {
        ReportMPIError(MPI_Send(&node_id, 1, MPI_LONG, senderID, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD), true);
        mpi_node_status->Store (senderID, 0, 1);    // reset busy signal
        mpi_node_status->Store (senderID, 1, node_id);
        ReportWarning (_String ("Sent child ") & node_id & " to node " & senderID);
    }
}


//___________________________________________________________________________________________

/* ------------------------------------------------------------------------------------
	Pass current network structure, data, constraint graph, and other analysis settings
	to compute node as HBL.
   ------------------------------------------------------------------------------------ */
void _BayesianGraphicalModel::SerializeBGMtoMPI (_StringBuffer & rec) {
  
    const static _String           kHYBgm_IMPUTE_MAXSTEPS    ("BGM_IMPUTE_MAXSTEPS"),
                                   kHYBgm_IMPUTE_BURNIN  ("BGM_IMPUTE_BURNIN"),
                                   kHYBgm_IMPUTE_SAMPLES ("BGM_IMPUTE_SAMPLES");

  
    rec << "\
USE_MPI_CACHING=1;\
PRINT_DIGITS=-1;\
function add_discrete_node (id,maxp,size,nlev) {\
  return {\
    'NodeID' : id,\
    'NodeType' : 0,\
    'MaxParents' : maxp.\
    'PriorSize'  : size,\
    'NumLevels' : nlev\
  };\
}\
function add_gaussian_node (id,maxp,size,mean,prec,scale) {\
  return {\
  'NodeID' : id,\
  'NodeType' : 1,\
  'MaxParents' : maxp.\
  'PriorVar'  : prec,\
  'NumLevels' : nlev,\
  'PriorScale' : scale\
  };\
}\
nodes = {};";

    for (long node_id = 0L; node_id < num_nodes; node_id++) {
        if (is_node_discrete (node_id)) {
            rec << "nodes + add_discrete_node("
                << *(_String*)node_names.GetItem (node_id) << max_parents.get (node_id) << ',' << prior_sample_size (node_id, 0) << ',' << num_levels.get (node_id) << ")'\n";
            // I can't figure out how to write the _String object from [node_names] to file :-P
         } else {
            rec << "nodes + add_gaussian_node("
              << *(_String*)node_names.GetItem (node_id) << max_parents.get (node_id) << ',' << prior_sample_size (node_id, 0) << ',' << prior_mean (node_id, 0) << ',' << prior_precision (node_id, 0) << ',' << prior_scale (node_id, 0) << ")'\n";
         }
    }

    // write BGM constructor
  
    _String const * bgm_name =(_String *) bgmNamesList (bgmList._SimpleList::Find((long)this));
  
    rec << "BayesianGraphicalModel " << *bgm_name << "(nodes);\n"
        << kHYBgm_IMPUTE_MAXSTEPS << '=' << hy_env :: EnvVariableGetNumber (kHYBgm_IMPUTE_MAXSTEPS,10000.) << ";\n"
        << kHYBgm_IMPUTE_BURNIN << '=' << hy_env :: EnvVariableGetNumber (kHYBgm_IMPUTE_BURNIN,1000.) << ";\n"
        << kHYBgm_IMPUTE_SAMPLES << '=' << hy_env :: EnvVariableGetNumber (kHYBgm_IMPUTE_SAMPLES,1000.) << ";\n";


  
    // serialize constraint matrix
    rec << "bgmConstraints = " << (_String *)constraint_graph.toStr() << ";\n";
    rec << "SetParameter(" << *bgm_name << ',' << kBGMContraintMatrix << ",bgmConstraints);\n";

    // serialize data matrix and assign to BGM
    rec << kBGMData << '=' << (_String *)theData.toStr() << ";\n";
    rec << "SetParameter(" << *bgm_name << ',' << kBGMData << ',' << kBGMData << ");\n";

}
#endif



//_________________________________________________________________________________________________________
_Matrix *   _BayesianGraphicalModel::Optimize (_AssociativeList const * options ) {
    /* ---------------------------------------------------------------------------------------------
        OPTIMIZE()
        Wrapper function to access various optimization methods (K2, structural MCMC, order MCMC)
            and implement MPI functionality.
       --------------------------------------------------------------------------------------------- */
    static const _String kHYBgm_METHOD_KEY      ("BGM_OPTIMIZATION_METHOD"),
                         kHYBgm_K2_RESTARTS     ("BGM_K2_RESTARTS"),
                         kHYBgm_K2_RANDOMIZE    ("BGM_K2_RANDOMIZE"),
                         kHYBgm_MCMC_NCHAINS       ("BGM_MCMC_NCHAINS"),
                         kHYBgm_MCMC_TEMP      ("BGM_MCMC_TEMPERATURE"),
                         kHYBgm_MCMC_MAXSTEPS  ("BGM_MCMC_MAXSTEPS"),
                         kHYBgm_MCMC_BURNIN        ("BGM_MCMC_BURNIN"),
                         kHYBgm_MCMC_SAMPLES       ("BGM_MCMC_SAMPLES");

    ReportWarning (_String("Entered _BayesianGraphicalModel::Optimize()"));

    if (!scores_cached) {
        CacheNodeScores();
    }
    


    hyFloat      optMethod =   hy_env::EnvVariableGetNumber(kHYBgm_METHOD_KEY, 0.0);    // 0 = K2 fixed order
    // 1 = K2 shuffle order with restarts
    // 2 = structural MCMC, fixed order
    // 3 = structural MCMC
    // 4 = order MCMC

    _Matrix *   output_matrix = nil;      // for order-MCMC and graph-MCMC has 4 columns: (1) chain trace;
    //  (2) edge marginal posteriors; (3) best state (node ordering or graph);
    //  (4) last state visited in chain

    ReportWarning (_String("... optimization method set to ") & optMethod);
  
    try {

      if (optMethod < 2.) {
          ReportWarning (_String("... starting K2 algorithm"));
        
        
          //output_matrix =  new _Matrix (num_nodes * num_nodes, 2, false, true);

        
           output_matrix = K2Search (optMethod,
                     hy_env::EnvVariableGetNumber(kHYBgm_K2_RESTARTS, 1.0),
                     hy_env::EnvVariableGetNumber(kHYBgm_K2_RANDOMIZE, num_nodes));
      } else {
        long      mcmc_steps = hy_env :: EnvVariableGetNumber(kHYBgm_MCMC_MAXSTEPS, 0.),
                  mcmc_burnin = hy_env :: EnvVariableGetNumber(kHYBgm_MCMC_BURNIN, 0.),
                  mcmc_samples = hy_env :: EnvVariableGetNumber(kHYBgm_MCMC_SAMPLES, 0.);

          if (mcmc_steps <= 0)    {
              throw (kHYBgm_MCMC_MAXSTEPS & " must be positive\n");
          }
          if (mcmc_burnin < 0 || mcmc_burnin >= mcmc_steps)    {
              throw(kHYBgm_MCMC_BURNIN & " must be positive and less than " & kHYBgm_MCMC_MAXSTEPS);
          }

          if (mcmc_samples < 0 || mcmc_samples > (mcmc_steps - mcmc_burnin))   {
              throw(kHYBgm_MCMC_SAMPLES & " must be positive and less than " & kHYBgm_MCMC_MAXSTEPS & " - " & kHYBgm_MCMC_BURNIN);
          }


  #if defined __HYPHYMPI__ && defined __AFYP_DEVELOPMENT__
          int         size,
                      rank;
          long        mpi_node;
          char        mpi_message [256];

          MPI_Status  status; // contains source, tag, and error code

          MPI_Comm_size (MPI_COMM_WORLD, &size);  // number of processes
          MPI_Comm_rank (MPI_COMM_WORLD, &rank);


           long mcmc_nchains = hy_env :: EnvVariableGetNumber(kHYBgm_MCMC_NCHAINS, 1.);
           if (mcmc_nchains <= 0) {
              throw _String("You must run at least one chain in MCMC.");
          }



          // parallel coupled chains (MPI)
          if (mcmc_nchains > 1) {
              hyFloat mcmc_dtemp = hy_env :: EnvVariableGetNumber (kHYBgm_MCMC_TEMP, 1.);
              if (mcmc_dtemp <= 1.) {
                  throw _String("Temperature increment must be greater than 1.");
              }

              hyFloat  chain_T = 1. / (1. + mcmc_dtemp*rank);


              /* note that GraphMetropolis already returns structure and last posterior prob as matrix entries */
              /* use sampling interval to swap chains */
              /* each chain needs a matrix pointer for output */
              if (optMethod < 4) {
                  output_matrix = ;
                  output_matrix = GraphMetropolis (optMethod==2, mcmc_burnin, mcmc_steps, mcmc_samples, chain_T);

                  for (long samp = 0; samp < mcmc_samples; samp++) {
                      /*** UNDER DEVELOPMENT ***/
                      if (rank > 0) {

                      } else {

                      }
                  }
              }
          }

          // single chain
          else {
  #endif
 
              if (optMethod < 4) {
                  ReportWarning (_String("... starting graph-mcmc"));
                  output_matrix = GraphMetropolis ( (optMethod == 2), mcmc_burnin, mcmc_steps, mcmc_samples, 1.);
              } else {
                  ReportWarning (_String("... starting order-mcmc"));
                  if (mcmc_burnin > 0) {
                      ReportWarning (_String("Executing order-MCMC for burn-in period of ") & mcmc_burnin & " steps");
                      output_matrix = OrderMetropolis (FALSE, mcmc_burnin, mcmc_samples, 1.);

                      ReportWarning (_String("Automatically reset node_order_arg to best order visited in order-MCMC burn-in:\n "));
                      if (node_order_arg.lLength == 0) {
                          node_order_arg.Populate (num_nodes, 0, 0);
                      }
                      for (long i = 0; i < num_nodes; i++) {
                          node_order_arg.list_data[i] = (*output_matrix) (i,3);
                      }
                      ReportWarning    (_String((_String*)node_order_arg.toStr()));
                      delete output_matrix;
                  }
                  ReportWarning (_String("Executing order-MCMC for ") & mcmc_steps & " steps, sample size " & mcmc_samples);
                  output_matrix = OrderMetropolis (TRUE, mcmc_steps, mcmc_samples, 1.);
              }

  #if defined __AFYP_DEVELOPMENT__ && defined __HYPHYMPI__
          }
  #endif

      }
    } catch (const _String & err) {
      HandleApplicationError(err);
      return new _Matrix;
    }


    return (_Matrix *) output_matrix;
}



//_________________________________________________________________________________________________________
_Matrix* _BayesianGraphicalModel::K2Search (bool do_permute_order, long n_restart, long n_randomize) {
    //  THIS NEEDS UPDATING
    return new _Matrix;
#ifdef __NEVER_DEFINED__
    hyFloat  this_score, best_score, next_score;

    _Matrix     order_matrix (num_nodes, num_nodes, false, true),
                best_dag (num_nodes, num_nodes, false, true);

    _SimpleList best_node_order;



    best_node_order.Populate (num_nodes, 0, 1);
    best_node_order.Permute (1);



    //  Convert node order to binary matrix where edge A->B is permitted if
    //  order_matrix[B][A] = 1, i.e. B is to the right of A in node order.
    for (long i = 0; i < num_nodes; i++) {
        long    child = best_node_order.list_data[i];

        for (long j = 0; j < num_nodes; j++) {
            long    parent = best_node_order.list_data[j];

            order_matrix.Store (parent, child, (j > i) ? 1 : 0);
        }
    }


    // greedy hill-climbing algorithm (K2)
    // ResetGraph (nil);
    best_score = Compute();     // best over all node orders (if we're reshuffling)


    for (long iter = 0; iter < n_restart; iter++) {
        next_score = Compute();     // reset to empty graph score


        for (long child = 0; child < num_nodes; child++) {
            long        num_parents = 0,
                        improvement_flag = 0,
                        next_parent_to_add;

            do {
                for (long parent = 0; parent < num_nodes; parent++) {
                    if ( (parent != child)                      // must meet all conditions!
                            && (theStructure(parent, child) == 0)
                            && (order_matrix (parent, child) == 1)
                            && constraint_graph(parent, child) >= 0) {
                        theStructure.Store (parent, child, 1);
                        this_score = Compute();

                        if (this_score > next_score) {
                            improvement_flag    = 1;
                            next_score          = this_score;
                            next_parent_to_add  = parent;
                        }
                        theStructure.Store (parent, child, 0);  // revert
                    }
                }

                if (improvement_flag) { // adding another parent improves network score
                    theStructure.Store (next_parent_to_add, child, 1);
                    num_parents++;
                    improvement_flag = 0;   // reset for next parent
                } else {
                    break;  // unable to improve further
                }
            } while (num_parents < max_parents.list_data[child]);
        }


        if (do_permute_order) {
            this_score = Compute();

            if (this_score > best_score) {
                best_score = this_score;
                best_dag = (_Matrix &) theStructure;        // store graph optimized from the current ordering
            }

            // ResetGraph (nil);

            best_node_order.Permute (1);

            for (long i = 0; i < num_nodes; i++) {
                long    child = best_node_order.list_data[i];
                for (long j = 0; j < num_nodes; j++) {
                    long    parent = best_node_order.list_data[j];
                    order_matrix.Store (parent, child, (j > i) ? 1 : 0);
                }
            }
        } else {
            break;  // only one iteration when node ordering is known
        }
    }

    if (do_permute_order) {
        theStructure = (_Matrix &) best_dag;
        best_node_order.Clear();    // reset to initial state
    }


    result->Store (0, 0, Compute());

    for (long row = 0; row < num_nodes; row++) {
        for (long col = 0; col < num_nodes; col++) {
            result->Store (row*num_nodes+col, 1, theStructure(row, col));
        }
    }
#endif
}



//_______________________________________________________________________________________________________________________
bool    _BayesianGraphicalModel::GraphObeysOrder (_Matrix & graph, _SimpleList const & order) const {
    _Matrix order_matrix (num_nodes, num_nodes, false, true);

    // convert order vector to matrix form
    for (long p_index = 0; p_index < num_nodes; p_index++) {
        for (long par = order.list_data[p_index], c_index = 0; c_index < num_nodes; c_index++) {
            order_matrix.Store (par, order.list_data[c_index], (p_index > c_index) ? 1 : 0);
        }
    }

    // loop thru graph, check that edges are compatible with node order
    for (long parent = 0; parent < num_nodes; parent++) {
        for (long child = 0; child < num_nodes; child++) {
            if (graph(parent, child)==1 && order_matrix(parent, child)==0) {
                return 0;
            }
        }
    }

    return 1;
}



//_______________________________________________________________________________________________________________________
_SimpleList const   _BayesianGraphicalModel::GetOrderFromGraph (_Matrix const & graph) const
{
    /* ---------------------------------------------------------------------------------
        GetOrderFromGraph()
            To quickly generate a node order based on graph argument.
            Loop through nodes in graph and insert into list according to parentage.
			*** Nodes can only be parents of nodes that appear to their left ***
            For an empty graph, this should return (0,1,2,...,num_nodes-1)
     ----------------------------------------------------------------------------------- */

    _SimpleList  new_order  (1, 0, 0); // initialize with single entry, [0]

    for (long left_of, node = 1; node < num_nodes; node++) {
        // loop through nodes in list looking for parents
        for (left_of = 0L; left_of < new_order.countitems(); left_of++) {
            // if the node is the child,
            if (graph (left_of, node)) {
                new_order.InsertElement ((BaseRef) node, left_of, false, false);
                break;
            }
        }

        // if we reach end of list, append
        if (left_of == new_order.countitems()) {
            new_order << node;
        }

    }
	ReportWarning(_String("Constructed node order from graph:\n") & _String ((_String *)new_order.toStr()) & "\n");
    return new_order;
}



//_______________________________________________________________________________________________________________________
_Matrix*    _BayesianGraphicalModel::GraphMetropolis (bool fixed_order, long mcmc_burnin, long mcmc_steps, long mcmc_samples,
        hyFloat chain_t)
{
    /* --------------------------------------------------------------------------------------------------------
        GraphMetropolis()

        Performs MCMC over structures.  Initialize chain using class member [theStructure] (accessible by HBL).
            If theStructure is empty graph, initialize order with random permutation.
            If [fixed_order]=TRUE, only explores the subset of graphs compatible with [node_order_arg] if specified.

        returns                 :   (1) vector of posterior probabiities
                                                (2) linearized matrix of edge marginal posteriors
                                                (3)     "       "   for best graph
                                                (4)     "       "   for last graph visited in chain
     
        throws const _String exceptions
       --------------------------------------------------------------------------------------------------------- */
  
   const static _String           kHYBgm_MCMC_PROBSWAP  ("BGM_MCMC_PROBSWAP"),
                                  kHYBgm_MCMC_MAXFAILS  ("BGM_MCMC_MAXFAILS");


    _Matrix     *   proposed_graph  = new _Matrix (num_nodes, num_nodes, false, true);
    _Matrix     *   result = new _Matrix ( (mcmc_samples > num_nodes*num_nodes) ? mcmc_samples : num_nodes*num_nodes, 4, false, true);

    _Matrix         current_graph (num_nodes, num_nodes, false, true),
                    best_graph    (num_nodes, num_nodes, false, true);

    hyFloat         current_score, proposed_score, best_score,
                    lk_ratio,
                    prob_swap = hy_env::EnvVariableGetNumber(kHYBgm_MCMC_PROBSWAP, 0.1);
  

    long            sampling_interval = mcmc_steps / mcmc_samples,
                    max_fails = hy_env::EnvVariableGetNumber(kHYBgm_MCMC_MAXFAILS, 100);

    _SimpleList *   proposed_order  = new _SimpleList();
    _SimpleList     current_order;
  
    _List           track_local_allocations; // use this to dispose of local allocations, including thrown exceptions
    track_local_allocations < proposed_graph < result < proposed_order;

    // parse HBL settings
    if (prob_swap < 0. || prob_swap > 1.) {
        throw (kHYBgm_MCMC_PROBSWAP & " must be assigned a value between 0 and 1");
    }

    if (max_fails <= 0.) {
        throw  (kHYBgm_MCMC_MAXFAILS & " must be assigned a value greater than 0");

    }



    // initialize chain state
    if (fixed_order) {
      if (node_order_arg.nonempty() && GraphObeysOrder (theStructure, node_order_arg) ) {
          (*proposed_graph) = theStructure;
          (*proposed_order) = node_order_arg;

          ReportWarning (_String ("Starting GraphMetropolis() using node_order_arg:\n ") & _String((_String *) proposed_order->toStr()));
        } else {
              throw _String("ERROR: Structure does not match order, aborting GraphMetropolis().");
        }
    } else {
      *proposed_order = GetOrderFromGraph (*proposed_graph);
    }

    
    // randomize the initial graph
    RandomizeGraph (proposed_graph, proposed_order, prob_swap, num_nodes*num_nodes, max_fails, fixed_order);
    ReportWarning (_String ("seeding with randomized graph:\n") & _String ((_String *) proposed_graph->toStr()));



    current_order = (*proposed_order);

    current_graph = (*proposed_graph);
    best_graph    = (*proposed_graph);

    best_score = proposed_score = current_score = Compute((_Matrix &) *proposed_graph);


    // main loop
    long entry = 0L;
    for (long step = 0; step < mcmc_steps + mcmc_burnin; step++) {
        //ReportWarning (_String ("current graph (") & current_score & "): \n" & (_String *) current_graph.toStr());
		
		// modify the graph by one edge
        RandomizeGraph (proposed_graph, proposed_order, prob_swap, 1, max_fails, fixed_order);


        proposed_score = Compute((_Matrix &) *proposed_graph);
        //ReportWarning (_String ("propose graph: (") & proposed_score & "):\n" & (_String *) proposed_graph->toStr());

        lk_ratio = exp(proposed_score - current_score);

        if (lk_ratio > 1. || genrand_real2() < lk_ratio) {  // Metropolis-Hastings
            current_graph   = *proposed_graph;        // accept proposed graph
            current_score   = proposed_score;
            current_order   = *proposed_order;
          

            if (StoreIfGreater(best_score, current_score)) {
                best_graph = current_graph;         // keep track of best graph ever visited by chain
            }
        } else {
            // revert proposal
            *proposed_order = current_order;
            *proposed_graph = current_graph;
        }

        // chain sampling
        if (step >= mcmc_burnin) {
            if ( (step- mcmc_burnin) % sampling_interval == 0L) {
              
                result->Store (entry++, 0, current_score);        // update chain trace
              
                for (long row = 0L; row < num_nodes; row++) {    // update edge tallies
                    for (long offset=row*num_nodes, col = 0L; col < num_nodes; col++) {
                        result->get_dense_numeric_cell (offset+col, 1) += current_graph(row, col);
                    }
                }
            }
        }
        // TODO SLKP 20180917 : implement progress updates in the new API, maybe use VERBOSITY_LEVEL
        /*
        if (step % ((mcmc_steps + mcmc_burnin)/100) == 0) {
             (_String ("graphMCMC at step ") & step & " of " & (mcmc_steps+mcmc_burnin));
        }
        */
      

    }
  
    hyFloat normalizer = 1./mcmc_samples;

    // convert edge tallies to frequencies, and report best and last graph
    for (long row = 0L; row < num_nodes; row++) {
        for (long offset=row*num_nodes, col = 0L; col < num_nodes; col++) {
            result->get_dense_numeric_cell (offset+col, 1) *= normalizer;
            result->get_dense_numeric_cell (offset+col, 2) = best_graph(row, col);
            result->get_dense_numeric_cell (offset+col, 3) = current_graph(row,col);
        }
    }


    // update theStructure with last graph visited by chain
    theStructure = (_Matrix &) current_graph;
    ReportWarning (_String("On exiting GraphMetropolic() assigned last state to theStructure: ") & (_String *) theStructure.toStr());

    result->AddAReference(); // do not clear track_local_allocations is deleted
    return result;
}



//_____________________________________________________________________________________________________________
void    _BayesianGraphicalModel::RandomizeGraph (_Matrix * graph, _SimpleList * order, hyFloat prob_swap, long num_steps, long max_fails, bool fixed_order) {
    //  Modify a graph by adding, removing, or reversing an edge, so long as that modification
    //  complies with the constraint matrix and (when fixed_order=TRUE) node order.
    //  throws const _String exceptions
  
  
      auto remove_random_parent_edge = [this, graph] (long node_index) -> bool {
      _SimpleList     removeable_edges;
      
        // build a list of current parents
      for (long par = 0L; par < num_nodes; par++)
          // par is parent of child AND the edge is NOT enforced
        if ((*graph)(par,node_index) && constraint_graph(par,node_index) <=0.) {
          removeable_edges << par;
        }
      
      if (removeable_edges.nonempty()) {
          // shuffle the list and remove the first parent
        graph->Store (removeable_edges.get (removeable_edges.Choice()), node_index, 0.);
        return true;
      } else {
          // none of the edges can be removed
        return false;
      }

    };

    if (num_nodes > 1) { // nothing to do
     
      long    step = 0L, fail = 0L;


      // calculate number of parents for each child for rapid access later
      _SimpleList     num_parents;

      num_parents.Populate (num_nodes, 0, 0);
      

      for (long child = 0; child < num_nodes; child++) {
          for (long parent = 0; parent < num_nodes; parent++) {
              if ( (*graph)(parent, child) > 0.) {
                  num_parents[child]++;
              }
          }

          if (num_parents.get(child) > max_parents (child)) {
              throw (_String ("Number of parents exceeds maximum BEFORE randomization of graph at node ") & child & " (" & num_parents.get(child) & " > " & max_parents.get(child) & ")\n" );
          }
      }


      // randomize graph
       do {
          // a fail-safe to avoid infinite loops
          if (fail > max_fails) {
              throw (_String ("Failed to modify the graph in RandomizeGraph() after ") & max_fails & " attempts.");
          }

          // pick a random edge
          long p_rank  = (genrand_int32() % (num_nodes-1)) + 1L,    // shift right to not include lowest node in order
               c_rank  =  genrand_int32() % p_rank,
               child  = order->get(c_rank),
               parent = order->get(p_rank);


          if (fixed_order || genrand_real2() > prob_swap) {
              if ( (*graph)(parent,child) == 0. && constraint_graph(parent,child) >= 0.) {
                  // add an edge if it is absent and not banned
                  if (num_parents.get(child) == max_parents.get (child)) {
                    if (remove_random_parent_edge (child)) {
                      graph->Store (parent, child, 1.);
                      step ++;
                    } else {
                      fail ++;
                    }
                  } else {
            // child can accept one more edge
                      graph->Store (parent, child, 1.);
                      num_parents[child]++;
                      step++;
                  }
              } else if ( (*graph)(parent,child) == 1. && constraint_graph(parent,child) <= 0.) {
          // delete an edge if it is present and not enforced
                  graph->Store (parent, child, 0.);
                  num_parents[child]--;
                  step++;
              } else {
                  fail++;
              }
          } else {    // swap nodes in ordering and flip edge if present
  
              //
              bool    ok_to_go = !((*graph)(parent,child) == 1.  &&  ( constraint_graph(child,parent)<0. || constraint_graph(parent,child)>0. )) &&
                        //P->C cannot be flipped if C->P is banned or P->C is enforced
                      !order->Any ([parent, child, graph, this] (long intermediate, unsigned long) -> bool {
                              return ( (*graph)(parent,intermediate)==1. && constraint_graph (parent, intermediate)>0. )  ||
                                     ( (*graph)(intermediate,child) ==1. && constraint_graph (intermediate, child)>0. )  ||
                                      constraint_graph (intermediate,parent) <0.  ||
                                      constraint_graph (child,intermediate)  <0.  ;

                      }, c_rank+1L)
            
                // by flipping the parent->child edge, we would screw up ordering for
                // an enforced edge:  C < B < P   becomes  P < B < C  where B->C or P->B is enforced.
                // also, a banned edge implies a node ordering constraint

             ;
          
            

              // if everything checks out OK
              if ( ok_to_go ) {
                  if ( (*graph)(parent,child) == 1 && constraint_graph(child,parent)>=0 ) {
            // flip the target edge
  
                      if (num_parents.get(parent) == max_parents.get(parent)) {
                        if (!remove_random_parent_edge (parent)) {
                          fail ++;
                          continue;
                        } else {
                          num_parents[parent]--;
                        }
                      }
                      graph->Store (parent, child, 0.);
                      graph->Store (child, parent, 1.);
                      num_parents[child]--;
                      num_parents[parent]++;

                      step++;
                  }
          // if the number of parents for parent node will exceed maximum, 
          // then the edge is deleted instead of flipped (delete without add)

                  // swap nodes in order
                  (*order)[p_rank] = child;
                  (*order)[c_rank] = parent;	// remember to update the order matrix also!

                  // update edges affected by node swap
                  //  child <-  N   ...   N <- parent                   _________________,
                  //                                    becomes        |                 v
                  //                                                 parent    N         N    child
                  //                                                           ^                |
                  //                                                           `----------------+
                  for (long bystander, i = c_rank+1; i < p_rank; i++) {
                      bystander = order->get(i);

                      if ( (*graph)(bystander, child) == 1. ) {
                          graph->Store (bystander, child, 0.);
                          num_parents[child]--;
                          graph->Store (child, bystander, 1.);
                          num_parents[bystander]++;
                        
                  
                          if (num_parents.get(bystander) > max_parents.get(bystander)) {
                              // number of parents for bystander node now exceeds maximum,
                              //  select a random edge to remove - including the edge C->B we just made!
                              remove_random_parent_edge (bystander);
                              // guaranteed to work because at least C->B is removable
                              num_parents[bystander]--;
                          }
                      }

                      if ( (*graph)(parent,bystander) == 1.) {
                          // flip edge from parent to bystander (P->B)
                          graph->Store (parent, bystander, 0.);
                          num_parents[bystander]--;

                          graph->Store (bystander, parent, 1.);
                          num_parents[parent]++;

                          if (num_parents.get(parent) > max_parents.get(parent)) {
                              // remove excess edge from X to parent other than bystander (or bystander as well)
                              remove_random_parent_edge (parent);
                              num_parents[parent]--;
                          }
                      }
                  }

                  step++;
              } else {
                  fail++;
              }
          }
      } while (step < num_steps);
    }
}



//_______________________________________________________________________________________________________________________
_Matrix *    _BayesianGraphicalModel::OrderMetropolis (bool do_sampling, long n_steps, long sample_size, hyFloat chain_t)
{
    /* ----------------------------------------------------------------------------------------
        OrderMetropolis()

        Execute Metropolis sampler using a swap of two nodes in an ordered sequence
        as a proposal function.  The posterior probabilities of edges in the network are
        stored as a member matrix [edge_posteriors].  Note that the total number of
        possible orderings (i.e. permutations of a sequence of length N) is factorial(N),
        and can possibly be computed exactly for N < 8.
       ---------------------------------------------------------------------------------------- */
   _Matrix     *   result = new _Matrix ( Maximum(n_steps , num_nodes*num_nodes), 4, false, true);
  
    long            sample_lag = n_steps / sample_size;

    hyFloat         lk_ratio,
                    prob_current_order, prob_proposed_order, best_prob,
                    denom;

    _SimpleList     current_order, proposed_order, best_node_order;

    _List           * marginals         = new _List ();

    _Vector  * gv;


    InitMarginalVectors (marginals);    // allocate storage




    // initialize node ordering
    if (node_order_arg.nonempty()) {
        current_order = node_order_arg;
    } else {
        current_order = GetOrderFromGraph (theStructure);
    }

    best_prob = prob_current_order = Compute (current_order, marginals);
    best_node_order = current_order;

    proposed_order.Populate (num_nodes, 0, 1);

    // chain
    for (long step = 0L; step < n_steps; step++) {
        // copy over current order to proposed order
        proposed_order = current_order;
        _SimpleList      pair_to_swap = proposed_order.Sample (2);

        proposed_order.Swap (pair_to_swap.get(0), pair_to_swap.get (1));


        // compute likelihood ratio of proposed : current orders
        prob_proposed_order = Compute (proposed_order, marginals);
        lk_ratio            = exp(prob_proposed_order - prob_current_order);

        if (lk_ratio > 1. || genrand_real2() < lk_ratio) {  // then set current to proposed order
            current_order = proposed_order;
            prob_current_order = prob_proposed_order;
          
            if (StoreIfGreater(best_prob, prob_proposed_order)) {
                best_node_order = proposed_order;
            }
        }


        // if past burn-in period and at sampling step, record trace and marginals
        if (do_sampling) {
            if (step % sample_lag == 0L) {
                ReportWarning(_String("At step ") & step & " order: " & _String((_String *) current_order.toStr()));
                result->Store (step / sample_lag, 0, prob_current_order);

                for (long child = 0; child < num_nodes; child++) {
                    // retrieve information from Compute()
                    gv      = (_Vector *) marginals->GetItem (child * num_nodes + child);
                    denom   = (*gv)(0, 0);  // first entry holds node marginal post pr

                    for (long edge, parent = 0; parent < num_nodes; parent++) {
                        if (parent == child) {
                            continue;
                        }

                        edge    = child * num_nodes + parent;
                        gv      = (_Vector *) marginals->GetItem (edge);

                        // not all _GrowingVector entries in marginals are being used,
                        //      i.e., edges incompatible with order
                        if (gv->get_used() > 0) {
                            // store as transpose so that row = parent and column = child, same as GraphMCMC
                            result->get_dense_numeric_cell(parent*num_nodes+child, 1) += exp (LogSumExpo(gv) - denom);
                        }
                    }
                }
            }
        }
    }


    // chain terminates


    // convert sums of edge posterior probs. in container to means
    hyFloat norm = 1./ sample_size;
    for (long edge = 0; edge < num_nodes * num_nodes; edge++) {
        result->get_dense_numeric_cell(edge, 1)*= norm;
    }


    // export node ordering info
    for (long node = 0; node < num_nodes; node++) {
        result->Store (node, 2, best_node_order.get(node));
        result->Store (node, 3, current_order.get(node));
    }


    DumpMarginalVectors (marginals);

    // update node order
    node_order_arg = current_order;
    ReportWarning (_String("Set node_order_arg to last order visited in orderMCMC:\n") & _String((_String*)node_order_arg.toStr()));
  
    return result;
}






