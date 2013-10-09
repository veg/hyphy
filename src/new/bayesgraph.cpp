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

#ifdef __HYPHYQT__
    #include "hyphymain.h"
#endif

#if defined __MAC__ || defined __WINDOZE__ || defined __HYPHY_GTK__
    #include "HYConsoleWindow.h"
    #include "HYDialogs.h"

#endif

extern  _Parameter  lnGamma (_Parameter);


_String     _HYBgm_NODE_INDEX   ("NodeID"),
            _HYBgm_NODETYPE       ("NodeType"),
            _HYBgm_NUM_LEVELS ("NumLevels"),
            _HYBgm_MAX_PARENT ("MaxParents"),
            _HYBgm_PRIOR_SIZE ("PriorSize"),
            _HYBgm_PRIOR_MEAN ("PriorMean"),      /* for continuous (Gaussian) nodes */
            _HYBgm_PRIOR_PRECISION    ("PriorPrecision"),
            _HYBgm_PRIOR_SCALE    ("PriorScale"),

            /*SLKP 20070926; add string constants for progress report updates */
            _HYBgm_STATUS_LINE_MCMC           ("Running Bgm MCMC"),
            _HYBgm_STATUS_LINE_MCMC_DONE  ("Finished Bgm MCMC"),
            _HYBgm_STATUS_LINE_CACHE      ("Caching Bgm scores"),
            _HYBgm_STATUS_LINE_CACHE_DONE ("Done caching Bgm scores"),
            /*SLKP*/

            _HYBgm_METHOD_KEY ("BGM_OPTIMIZATION_METHOD"),
            _HYBgm_MPI_CACHING ("USE_MPI_CACHING"),

            _HYBgm_K2_RESTARTS ("BGM_K2_RESTARTS"),
            _HYBgm_K2_RANDOMIZE ("BGM_K2_RANDOMIZE"),

            _HYBgm_MCMC_NCHAINS       ("BGM_MCMC_NCHAINS"),
            _HYBgm_MCMC_TEMP      ("BGM_MCMC_TEMPERATURE"),
            _HYBgm_MCMC_MAXSTEPS  ("BGM_MCMC_MAXSTEPS"),
            _HYBgm_MCMC_BURNIN        ("BGM_MCMC_BURNIN"),
            _HYBgm_MCMC_SAMPLES       ("BGM_MCMC_SAMPLES"),

            _HYBgm_MCMC_PROBSWAP  ("BGM_MCMC_PROBSWAP"),
            _HYBgm_MCMC_MAXFAILS  ("BGM_MCMC_MAXFAILS"),

            _HYBgm_IMPUTE_MAXSTEPS    ("BGM_IMPUTE_MAXSTEPS"),
            _HYBgm_IMPUTE_BURNIN  ("BGM_IMPUTE_BURNIN"),
            _HYBgm_IMPUTE_SAMPLES ("BGM_IMPUTE_SAMPLES"),

            _HYBgm_CONTINUOUS_MISSING_VALUE ("BGM_CONTINUOUS_MISSING_VALUE");


//__________________________________________________________________________________________________________
#ifdef      __UNIX__

void        ConsoleBGMStatus (_String, _Parameter, _String * fileName = nil);


void        ConsoleBGMStatus (_String statusLine, _Parameter percentDone, _String * fileName)
{
    FILE           *outFile = fileName?doFileOpen (fileName->sData,"w"):nil;
    _String        reportLine (statusLine);


    if (percentDone >= 0.0) {
        reportLine = reportLine & ". " & percentDone & "% done.";
    }

    if (outFile) {
        fprintf (outFile,"%s", reportLine.sData);
    } else if (verbosityLevel == 1) {
        printf ("\033\015 %s", reportLine.sData);
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
long        integerPower (long base, long exponent)
{
    //  Rapid computation of an integer power.
    long    result = 1,
            mask   = 1L<<(sizeof(long)*8-2); // left shift to left-most position of binary sequence for long integer
    // e.g. 100...0 (30 zeroes for signed long)

    while ((exponent & mask) == 0) {
        mask >>= 1;    // bitwise AND, right-shift mask until overlaps with first '1'
    }

    while (mask) {
        result *= result;
        if (exponent & mask) {
            result = result * base;
        }
        mask >>= 1;
    }
    return result;
}





//__________________________________________________________________________________________________________
_Parameter      LogSumExpo (_GrowingVector * log_values)
{
    //  Computes the sum of a vector whose values are stored as log-transforms,
    //  such that exponentiating the vector entries would result in numerical underflow.

    long        size            = log_values->GetUsed();
    _Parameter  sum_exponents   = 0.;


    // handle some trivial cases
    if (size == 0) {
        return 0.;    // log(exp(0)) = log(1) = 0
    } else if (size == 1) {
        return (*log_values)(0,0);    // log(exp(log(x)))
    }


    // find the largest (least negative) log-value
    _Parameter      max_log  = (*log_values) (0, 0),
                    this_log;

    for (long val = 1; val < size; val++) {
        this_log = (*log_values) (val, 0);

        if (this_log > max_log) {
            max_log = this_log;
        }
    }


    // go back through log values and increment by max log value
    //      This will cause some underflow for the smallest values, but we can
    //  use this approximation to handle very large ranges of values.
    //  NOTE: subtracting a negative value.
    for (long val = 0; val < size; val++) {
        sum_exponents += exp( (*log_values) (val, 0) - max_log );
    }

    return (log(sum_exponents) + max_log);
}




//__________________________________________________________________________________________________________
//__________________________________________________________________________________________________________
//__________________________________________________________________________________________________________
_BayesianGraphicalModel::_BayesianGraphicalModel (_AssociativeList * nodes)
{
    _String             errorMessage;
    _AssociativeList *  this_avl;
    long                global_max_parents  = 0;
    _Constant *         avl_val;

    node_names.Clear();


    // the fundamental defining characteristic of _BayesianGraphicalModel objects
    num_nodes           = nodes->avl.countitems();


    // allocate space to class matrices
    CreateMatrix (&theStructure, num_nodes, num_nodes, false, true, false);

    CreateMatrix (&prior_sample_size, num_nodes, 1, false, true, false);
    CreateMatrix (&prior_mean, num_nodes, 1, false, true, false);
    CreateMatrix (&prior_precision, num_nodes, 1, false, true, false);
    CreateMatrix (&prior_scale, num_nodes, 1, false, true, false);

    CreateMatrix (&constraint_graph, num_nodes, num_nodes, false, true, false);


    // allocate space for _SimpleList objects
    node_type.Populate (num_nodes, 0, 0);
    max_parents.Populate (num_nodes, 0, 0);
    num_levels.Populate(num_nodes, 0, 0);
    has_missing.Populate(num_nodes, 0, 0);


    for (long node = 0; node < num_nodes; node++) {
        this_avl = (_AssociativeList *) (nodes->GetByKey (node, ASSOCIATIVE_LIST));

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
            WarnError("Invalid node name (expected a string) passed to a BGM constructor");
            return;
        }
        
        node_names.AppendNewInstance(new _String (*name->theString));  // append a pointer to _String duplicate

        // DEBUGGING
        /* wuz: ReportWarning (_String("node_name[") & node & "]=" & (_String *)node_names.lData[node]);
         20111210 SLKP : _String (_String*) contstructor will actually assume that the argument is 
                       : 'unencumbered' i.e. not a member of _Lists etc
                       : this function call will create a stack copy of node_names.lData[node]
                       : print it to messages.log and then kill the dynamic portion of the object (sData)
                       : this will create all kinds of havoc downstream
         */
        // bug fix:
        ReportWarning (_String("node_name[") & node & "]=" & *((_String *)node_names.lData[node]));
        

        // node type (0 = discrete, 1 = continuous)
        if ((avl_val = (_Constant *) (this_avl->GetByKey (_HYBgm_NODETYPE, NUMBER)))) {
            node_type.lData[node] = (long)(avl_val->Value());
            if (node_type.lData[node] < 0 || node_type.lData[node] > 2) {
                errorMessage = _String("Unsupported NodeType ") & node_type.lData[node] & " for node " & node;
                break;
            }
        } else {
            errorMessage = _String ("Missing NodeType in associative array for node ") & node;
            break;
        }


        // number of levels (discrete nodes)
        if ((avl_val = (_Constant *) (this_avl->GetByKey (_HYBgm_NUM_LEVELS, NUMBER)))) {
            num_levels.lData[node] = (long)(avl_val->Value());

            if (num_levels.lData[node] <= 1) {
                errorMessage = _String("NLevels must be greater than 1, received ") & num_levels.lData[node] & " for node " & node;
                break;
            }
        } else {
            if (!avl_val && node_type.lData[node] == 0) {
                errorMessage = _String ("Missing NumLevels in associative array for node ") & node;
                break;
            }
        }


        // max parents
        if ((avl_val = (_Constant *) (this_avl->GetByKey (_HYBgm_MAX_PARENT, NUMBER)))) {
            max_parents.lData[node] = (long)(avl_val->Value());

            if (max_parents.lData[node] > global_max_parents) {
                global_max_parents = max_parents.lData[node];
            }

            if (max_parents.lData[node] <= 0) {
                errorMessage = _String ("MaxParents must be greater than zero, received ") & max_parents.lData[node] & " for node " & node;
                break;
            }

            if (max_parents.lData[node] >= num_nodes) {
                errorMessage = _String ("MaxParents cannot equal or exceed the number of nodes in the network, see node ") & node;
                break;
            }
        } else {
            errorMessage = _String ("Missing MaxParents in associative array for node ") & node;
            break;
        }


        // prior sample size
        if ((avl_val = (_Constant *) (this_avl->GetByKey (_HYBgm_PRIOR_SIZE, NUMBER)))) {
            prior_sample_size.Store (node, 0, (_Parameter) (avl_val->Value()));

            if (prior_sample_size(node,0) < 0) {
                errorMessage = _String ("PriorSampleSize must be at least zero, received ") & prior_sample_size(node,0) & " for node " & node;
                break;
            }
        } else {
            errorMessage = _String ("Missing PriorSize in associative array for node ") & node;
            break;
        }


        // prior mean (Gaussian)
        if ((avl_val = (_Constant *) (this_avl->GetByKey (_HYBgm_PRIOR_MEAN, NUMBER)))) {
            prior_mean.Store (node, 0, (_Parameter) (avl_val->Value()));
        } else if (!avl_val && node_type.lData[node] == 1) {
            errorMessage = _String ("Missing PriorMean in associative array for node ") & node;
            break;
        }


        // prior precision (Gaussian)
        if ((avl_val = (_Constant *) (this_avl->GetByKey (_HYBgm_PRIOR_PRECISION, NUMBER)))) {
            prior_precision.Store (node, 0, (_Parameter) (avl_val->Value()));

            if (avl_val <= 0) {
                errorMessage = _String ("PriorPrecision must be greater than zero, received ") & prior_precision(node,0) & " for node " & node;
                break;
            }
            // ReportWarning (_String("prior_precision[") & node & "] set to " & prior_precision(node,0) );
        } else if (!avl_val && node_type.lData[node] == 1) {
            errorMessage = _String ("Missing PriorPrecision in associative array for node ") & node;
            break;
        }


        // prior scale (Gaussian)
        if ((avl_val = (_Constant *) (this_avl->GetByKey (_HYBgm_PRIOR_SCALE, NUMBER)))) {
            prior_scale.Store (node, 0, (_Parameter) (avl_val->Value()));

            if (avl_val <= 0) {
                errorMessage = _String ("PriorScale must be greater than zero, received ") & prior_scale(node,0) & " for node " & node;
                break;
            }
        } else if (!avl_val && node_type.lData[node] == 1) {
            errorMessage = _String ("Missing PriorScale in associative array for node ") & node;
            break;
        }
    }


    if (errorMessage.sLength) {
        WarnError (errorMessage);
    }


    // allocate node score cache
    _List   emptyList (global_max_parents + 1);

    for (long node = 0; node < num_nodes; node++) {
        node_score_cache && (&emptyList);	// appends a dynamic copy of _List object to start
    }

    scores_cached = FALSE;

    ReportWarning (_String ("Constructed BayesianGraphicalModel with ") & num_nodes & " nodes.");
}



//__________________________________________________________________________________________________________
_BayesianGraphicalModel::~_BayesianGraphicalModel (void)
{
    /* destructor */
}




//__________________________________________________________________________________________________________
bool _BayesianGraphicalModel::SetDataMatrix (_Matrix * data)
{
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
	for (long node = 0; node < num_nodes; node++) {
		has_missing.lData[node] = 0;
	}
	
    /* check for user assignment of continuous missing value indicator */
    checkParameter (_HYBgm_CONTINUOUS_MISSING_VALUE, continuous_missing_value, -666.0);
	ReportWarning (_String ("Entered SetDataMatrix() with missing CG flag: ") & continuous_missing_value & " and node types" & (_String *) node_type.toStr());
	
    data_nlevels.Populate (num_nodes, 1, 0);

    if (data->GetVDim() == num_nodes) {
        // if (theData) DeleteObject (theData);
        // purge duplicate from memory (makeDynamic)

        theData = (_Matrix &) *data;        // duplicate matrix to member variable
        theData.CheckIfSparseEnough (TRUE);

        scores_cached = FALSE;


        for (long nrows = theData.GetHDim(), node = 0; node < num_nodes; node++) {
            // if discrete node, compute number of levels and missingness
            if (node_type.lData[node] == 0) {
                data_nlevels.lData[node] = 1;

                for (long val, row = 0; row < nrows; row++) {
                    val = theData (row, node);

                    if (val < 0 && has_missing.lData[node] == 0) {
                        has_missing.lData[node] = 1;
                        continue;
                    }

                    if (val + 1 > data_nlevels.lData[node]) {
                        data_nlevels.lData[node] = data_nlevels.lData[node] + 1;
                    }
                }

                if (data_nlevels.lData[node] != num_levels.lData[node]) {
                    WarnError (_String ("ERROR: Number of levels in data (") & data_nlevels.lData[node] & ") for discrete node "
                               & node & " is not compatible with node setting (" & num_levels.lData[node]
                               & ").  Check your data or reset the BayesianGraphicalModel.");

                    return (FALSE);
                }
            }

            // continuous
            else if (node_type.lData[node] == 1) {
                for (long val, row = 0; row < nrows; row++) {
					val = theData (row, node);
                    if (val == continuous_missing_value && has_missing.lData[node] == 0) {
                        has_missing.lData[node] = 1;
						ReportWarning (_String("Detected missing continuous value at row ") & row);
                        break;
                    }
                }
            }
        }

        ReportWarning (_String ("Set data matrix to:\n") & (_String *)theData.toStr() & "\n" & " and missing values at " & (_String *) has_missing.toStr());
    } else {
        WarnError (_String("ERROR: Number of variables in data (") & data->GetVDim() & ") does not match number of nodes in graph (" & num_nodes & ")");
        return (FALSE);
    }


    // compute node scores and store in cache
    CacheNodeScores();


    return (TRUE);
}



bool _BayesianGraphicalModel::SetWeightMatrix (_Matrix * weights)
{
    if (weights->GetHDim() == theData.GetHDim() && weights->GetHDim() == num_nodes) {
        theWeights = (_Matrix &) *weights;
        ReportWarning(_String("Assigned weight matrix:\n") & (_String *) theWeights.toStr());
        return TRUE;
    }

    WarnError (_String("Incompatible matrix dimensions in SetWeightMatrix()."));
    return FALSE;
}


//__________________________________________________________________________________________________________
bool _BayesianGraphicalModel::SetConstraints (_Matrix * constraints)
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
    if (constraints->GetHDim() == num_nodes) {
        constraint_graph = (_Matrix &) (*constraints);
        ReportWarning (_String("Assigned constraint matrix:\n ") & (_String*) constraint_graph.toStr() );
        return (TRUE);
    }

    _String errorMsg ("ERROR: Constraint matrix incompatible dimensions to graph.");
    WarnError (errorMsg);
    return (FALSE);
}



//__________________________________________________________________________________________________________
bool _BayesianGraphicalModel::SetStructure (_Matrix * structure)
{
    /* -------------------------------------------------------------------
        SetStructure()
            Assign pointer to _Matrix object from batchlan:SetParameter()
            that specificies a network structure.
            Check the structure against constraint matrix and node
            ordering if the latter is set.
            If node order is incompatible, reset the node order.
       ------------------------------------------------------------------- */

    if (structure->GetHDim() == num_nodes) {
        // check graph against constraint matrix
        for (long row = 0; row < num_nodes; row++) {
            for (long col = 0; col < num_nodes; col++) {
                if (constraint_graph(row,col) < 0 && (*structure)(row,col) == 1) {
                    _String errorMsg ("ERROR: Structure contains banned edge: ");
                    errorMsg = errorMsg & _String (row) & _String ("->") & _String (col);
                    WarnError (errorMsg);
                    return (FALSE);
                }

                if (constraint_graph(row,col) > 0 && (*structure)(row,col) == 0) {
                    _String errorMsg ("ERROR: Structure lacks enforced edge:");
                    errorMsg = errorMsg & _String (row) & _String ("->") & _String (col);
                    WarnError (errorMsg);
                    return (FALSE);
                }
            }
        }


        // is new structure compatible with node order set in HBL?
        if (node_order_arg.lLength == num_nodes && !GraphObeysOrder (theStructure, node_order_arg)) {
            // need to reset node_order
            _SimpleList * templist;

            templist = GetOrderFromGraph (theStructure);
            node_order_arg = (_SimpleList &) (*templist);
            DeleteObject (templist);

            ReportWarning (_String ("Structure is incompatible with existing node order, resetting order."));
        }

        theStructure = (_Matrix &) (*structure);
        return (TRUE);
    }

    _String errorMsg ("ERROR: Structure incompatible dimensions to graph.");
    WarnError (errorMsg);
    return (FALSE);
}


void _BayesianGraphicalModel::GetStructure (_Matrix * graph)
{
    for (long row = 0; row < num_nodes; row++)
        for (long col = 0; col < num_nodes; col++) {
            graph->Store (row, col, theStructure(row, col));
        }

    ReportWarning (_String("GetStructure() copied graph ") & (_String *) graph->toStr());
}

//__________________________________________________________________________________________________________
bool _BayesianGraphicalModel::SetNodeOrder (_SimpleList * order)
{
    /* -----------------------------------------------------------------
        SetNodeOrder()
            Set node ordering to vector argument from HBL SetParameter().
       ----------------------------------------------------------------- */

    if (order->lLength == num_nodes) {
        if (GraphObeysOrder (theStructure, (_SimpleList &) *order)) {
            node_order_arg.Populate(num_nodes, 0, 0);   // reset member variable

            for (long i = 0; i < num_nodes; i++) {
                node_order_arg.lData[i] = order->lData[i];
            }

            ReportWarning (_String("BayesianGraphicalModel node order arg set to ") & (_String *) node_order_arg.toStr());

            return (TRUE);
        } else {
            _String errorMsg ("ERROR: Node order incompatible with current graph.");
            WarnError (errorMsg);
            return (FALSE);
        }
    } else {
        _String errorMsg ("ERROR: Node order argument incorrect length.");
        WarnError (errorMsg);
        return (FALSE);
    }
}


void _BayesianGraphicalModel::GetNodeOrder (_Matrix * receptacle)
{
    if (node_order_arg.lLength == num_nodes) {
        for (long node = 0; node < num_nodes; node++) {
            receptacle->Store (0, node, node_order_arg.lData[node]);
        }
    }
}


//__________________________________________________________________________________________________________
_Parameter _BayesianGraphicalModel::Compute (void)
{
    /* --------------------------------------------------------------
        Compute()   [POLYMORPHIC]
            Return posterior probability of [CURRENT] network
            structure.
       -------------------------------------------------------------- */

    _Parameter  log_score = 0.;

    for (long node_id = 0; node_id < num_nodes; node_id++) {
        log_score += node_type.lData[node_id] ? ComputeContinuousScore (node_id) : ComputeDiscreteScore (node_id);
    }

    return log_score;
}



//__________________________________________________________________________________________________________
_Parameter _BayesianGraphicalModel::Compute (_Matrix & g)
{
    /* --------------------------------------------------------------
        Compute()   [POLYMORPHIC]
            Return posterior probability of [GIVEN] network
            structure argument.
       -------------------------------------------------------------- */

    _Parameter  log_score = 0.;

    for (long node_id = 0; node_id < num_nodes; node_id++) {    // discrete nodes are 0 (FALSE)
        log_score += node_type.lData[node_id] ? ComputeContinuousScore (node_id, g) : ComputeDiscreteScore (node_id, g);
    }

    return log_score;
}



//__________________________________________________________________________________________________________
_Parameter  _BayesianGraphicalModel::Compute (_SimpleList & node_order, _List * marginals)
{
    /* --------------------------------------------------------------
        Compute()   [POLYMORPHIC]
            Return posterior probability of given node order by
            integrating over all compatible structures.
            Also return marginal posterior probabilities for edges.
       -------------------------------------------------------------- */

    _Parameter          log_likel   = 0.;
    _GrowingVector      *gv1, *gv2;


    // reset _GrowingVector objects stored in _List object
    for (long i = 0; i < num_nodes * num_nodes; i++) {
        gv1 = (_GrowingVector *) marginals->lData[i];
        gv1 -> ZeroUsed();
    }


    for (long nodeIndex = 0; nodeIndex < node_order.lLength; nodeIndex++) {
        long                child_node      = node_order.lData[nodeIndex],
                            maxp            = max_parents.lData[child_node];

        _List           *   score_lists     = (_List *) node_score_cache.lData[child_node];
        _Constant       *   orphan_score    = (_Constant *) (score_lists->lData[0]);


        gv1 = (_GrowingVector *) marginals->lData[child_node * num_nodes + child_node]; // store denominator in diagonal
        gv1->ZeroUsed();
        gv1 -> Store (orphan_score->Value());   // append score - note gv1 does not change within



        if (maxp > 0) {
            // all nodes to the right are potential parents, except banned parents!
            _SimpleList     precedes;
            for (long parIndex = nodeIndex + 1; parIndex < node_order.lLength; parIndex++) {
                long    par = node_order.lData[parIndex];

                if (constraint_graph(par, child_node) >= 0) {   // not banned
                    precedes << par;
                }
            }


            // handle trivial case of one parent
            _Matrix *   single_parent_scores    = (_Matrix *) (score_lists->lData[1]);

            for (long i = 0; i < precedes.lLength; i++) {
                long    par = precedes.lData[i];

                gv1 -> Store ((*single_parent_scores) (par, 0)); // append to denominator
                gv2 = (_GrowingVector *) marginals->lData[child_node * num_nodes + par]; // update parent-specific numerator
                gv2 -> Store ((*single_parent_scores) (par, 0));
            }


            // more than one parent requires k-tuples
            if (maxp > 1) {
                _SimpleList         indices (precedes.lLength, 0, 1);   // populates list with 0, 1, 2, ..., M-1
                // where M is the number of eligible parents in [precedes]
                _NTupleStorage *    family_scores;

                for (long nparents = 2; nparents <= maxp; nparents++) {
                    _SimpleList     subset,
                                    auxil;

                    bool            not_finished;


                    if (nparents > precedes.lLength) {  // not enough eligible parents to form tuples!
                        break;
                    }


                    if (indices.NChooseKInit (auxil, subset, nparents, false)) {
                        _Parameter      tuple_score;
                        _SimpleList     parents;

                        parents.Populate (nparents, 0, 0);  // allocate memory

                        family_scores = (_NTupleStorage *) (score_lists->lData[nparents]);


                        do {
                            //parents.Clear();
                            not_finished = indices.NChooseK (auxil, subset);    // cycle through index combinations


                            for (long i = 0; i < nparents; i++) {   // convert indices to parent IDs (skipping child)
                                long    realized = precedes.lData[subset.lData[i]];
                                if (realized >= child_node) {
                                    realized--;
                                }
                                parents.lData[i] = realized;
                            }
                            parents.Sort(TRUE);

                            tuple_score = family_scores -> Retrieve (parents);

                            gv1 -> Store (tuple_score); // append to denominator

                            for (long i = 0; i < nparents; i++) {
								// update parent combination-specific numerator
                                gv2 = (_GrowingVector *) marginals->lData[child_node * num_nodes + precedes.lData[subset.lData[i]]];
                                gv2 -> Store (tuple_score);
                            }
                        } while (not_finished);
                    }
                }
            }
        }

        gv1 -> _Matrix::Store (0, 0, LogSumExpo(gv1));  // replace first entry with sum, i.e. marginal log-likelihood of child node
        log_likel += (*gv1)(0, 0);

    }
    // end loop over child nodes

    return log_likel;
}



//__________________________________________________________________________________________________________
_Parameter  _BayesianGraphicalModel::ComputeDiscreteScore (long node_id)
{
    _SimpleList     parents;

    for (long par = 0; par < num_nodes; par++) {
        if (theStructure(par, node_id) == 1 && node_type.lData[par] == 0) {
            parents << par;
        }
    }

    return ComputeDiscreteScore (node_id, parents);
}


//__________________________________________________________________________________________________________
_Parameter  _BayesianGraphicalModel::ComputeDiscreteScore (long node_id, _Matrix & g)
{
    _SimpleList     parents;

    for (long par = 0; par < num_nodes; par++) {
        if ( g(par, node_id) == 1 && node_type.lData[par] == 0) {
            parents << par;
        }
    }

    return ComputeDiscreteScore (node_id, parents);
}


//__________________________________________________________________________________________________________
_Parameter  _BayesianGraphicalModel::ComputeDiscreteScore (long node_id, _SimpleList & parents)
{
    /* --------------------------------------------------------------------
        ComputeDiscreteScore()
            Returns posterior probability of local network structure
            centered at discrete child node.
            Only discrete nodes may be parents of a discrete child node.
       -------------------------------------------------------------------- */


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


    //ReportWarning (_String ("Non-cached call of ComputeDiscreteScore with ") & node_id & " <- " & (_String *) parents.toStr());

    // impute score if missing data
    if (has_missing.lData[node_id]) {
        return ImputeDiscreteNodeScore (node_id, parents);
    } else {
        for (long par = 0; par < parents.lLength; par++) {
            if (has_missing.lData[parents.lData[par]]) {
                return ImputeDiscreteNodeScore (node_id, parents);
            }
        }
    }

    _Matrix         n_ijk,
                    n_ij;       // [i] indexes child nodes,
    // [j] indexes combinations of values for parents of i-th node,
    // [k] indexes values of i-th node

    UpdateDirichletHyperparameters (node_id, parents, &n_ij, &n_ijk);

    return ( (prior_sample_size (node_id, 0) == 0) ?
             K2Score (node_id, n_ij, n_ijk) :
             BDeScore (node_id, n_ij, n_ijk) );
}


//_______________________________________________________________
void    _BayesianGraphicalModel::UpdateDirichletHyperparameters (long dnode, _SimpleList &dparents, _Matrix * n_ij, _Matrix * n_ijk)
{
    if (node_type.lData[dnode] > 0) {
        ReportWarning ("ERROR: UpdateDirichletHyperparameters() called on non-discrete node!  That sucks!");
    }


    if (dparents.lLength > 0) {
        _SimpleList     multipliers ((long)1);
        long            num_parent_combos   = 1;

        for (long par = 0; par < dparents.lLength; par++) {
            num_parent_combos *= num_levels.lData[dparents.lData[par]];
            multipliers << num_parent_combos;
        }

        CreateMatrix (n_ij, num_parent_combos, 1, false, true, false);
        CreateMatrix (n_ijk, num_parent_combos, num_levels.lData[dnode], false, true, false);


        // initialize matrices with prior hyperparameters (a_ijk) -- if using K2, entries will remain zero
        for (long j = 0; j < num_parent_combos; j++) {
            n_ij->Store (j, 0, prior_sample_size(dnode,0) / num_parent_combos);

            for (long k = 0; k < num_levels.lData[dnode]; k++) {
                n_ijk->Store (j, k, (*n_ij)(j,0) / num_levels.lData[dnode]);
            }
        }


        // update hyperparameters with data
        /*
            What should we do with incomplete cases?  Currently using Gibbs sampling to impute
                an average node score - should we export one of these samples?
            For now, just export complete cases - BUT THIS WILL CAUSE PROBLEMS IF NONE OF
                THE CASES IS COMPLETE! - afyp, 2010/02/25
         */
        for (long obs = 0; obs < theData.GetHDim(); obs++) {
            long    index           = 0,
                    child_state     = theData(obs, dnode);

            if (child_state < 0) {
                continue;    /* missing observation */
            }

            for (long par = 0; par < dparents.lLength; par++) {
                long    this_parent         = dparents.lData[par],
                        this_parent_state = theData(obs, this_parent);

                if (this_parent_state < 0) {
                    index = -1; /* missing observation */
                    break;
                }
                index += this_parent_state * multipliers.lData[par];
            }

            if (index >= 0) {
				// this is where we would modify the increment value (1) if we are weighting cases...
                n_ijk->Store ((long) index, child_state, (*n_ijk)(index, child_state) + 1);
                n_ij->Store ((long) index, 0, (*n_ij)(index, 0) + 1);
            }
        }
    }

    else {
        // conditionally independent node
        CreateMatrix (n_ij, 1, 1, false, true, false);
        CreateMatrix (n_ijk, 1, num_levels.lData[dnode], false, true, false);

        // prior hyperparameters (uniform distribution)
        for (long k = 0; k < num_levels.lData[dnode]; k++) {
            n_ijk->Store(0, k, prior_sample_size(dnode,0) / num_levels.lData[dnode]);
        }

        // update with data
        for (long obs = 0; obs < theData.GetHDim(); obs++) {
            long child_state = theData(obs, dnode);

            if (child_state < 0) {
                continue;
            }

			// this is where we would modify the increment value (1) if we are weighting cases...
            n_ijk->Store (0, child_state, (*n_ijk)(0, child_state) + 1);
            n_ij->Store (0, 0, (*n_ij)(0,0) + 1);
        }
    }
}



//___________________________________________________________________________________________
_Parameter _BayesianGraphicalModel::K2Score (long node_id, _Matrix & n_ij, _Matrix & n_ijk)
{
    _Parameter  log_score   = 0.;
    long        r_i         = num_levels.lData[node_id];

    for (long j = 0; j < n_ij.GetHDim(); j++) {
        log_score += lnGamma(r_i);  // (r-1)!
        log_score -= lnGamma(n_ij(j, 0) + r_i); // (N+r-1)!

        for (long k = 0; k < r_i; k++) {
            log_score += lnGamma(n_ijk(j,k) + 1);   // (N_ijk)!
        }
    }

    return log_score;
}



//___________________________________________________________________________________________
_Parameter _BayesianGraphicalModel::BDeScore (long node_id, _Matrix & n_ij, _Matrix & n_ijk)
{
    // note that n_ij and n_ijk already contain prior counts, updated by data
    _Parameter  n_prior_ij      = prior_sample_size (node_id, 0) / n_ij.GetHDim(),
                n_prior_ijk     = n_prior_ij / num_levels.lData[node_id],
                log_score        = 0.;

    for (long j = 0; j < n_ij.GetHDim(); j++) {
        log_score += lnGamma(n_prior_ij) - lnGamma(n_ij(j,0));

        for (long k = 0; k < num_levels.lData[node_id]; k++) {
            log_score += lnGamma(n_ijk(j,k)) - lnGamma(n_prior_ijk);
        }
    }

    return log_score;
}



//___________________________________________________________________________________________
void    _BayesianGraphicalModel::InitMarginalVectors (_List * compute_list)
{
    /* ------------------------------------------------------------------------------------
        InitMarginalVectors()
            Allocate storage of node and edge marginal posterior probabilities accumulated
            during order-MCMC.  Off-diagonals correspond to entries of an adjacency matrix
            (edges), whereas diagonal entries are used to store node marginals.
       ------------------------------------------------------------------------------------ */

    _GrowingVector * newstore;
    checkPointer (newstore = new _GrowingVector);

    for (long i = 0; i < num_nodes * num_nodes; i++) {
        (*compute_list) && newstore;
    }

    DeleteObject (newstore);
}



//___________________________________________________________________________________________
void    _BayesianGraphicalModel::DumpMarginalVectors (_List * compute_list)
{
    for (long i = 0; i < compute_list->lLength; i++) {
        ((_GrowingVector *) compute_list->lData[i]) -> Clear();
    }

    DeleteObject (compute_list);
}




//___________________________________________________________________________________________
void    _BayesianGraphicalModel::CacheNodeScores (void)
{
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

    ReportWarning (_String ("Entered CacheNodeScores()"));

    if (scores_cached) {
        return;
    }


#if defined __HYPHYMPI__
    _Parameter  use_mpi_caching;
    checkParameter (_HYBgm_MPI_CACHING, use_mpi_caching, 0);

    if (use_mpi_caching) {
        ReportWarning (_String ("Using MPI to cache node scores."));

        // MPI_Init() is called in main()
        int             size, rank;
        long            mpi_node;
        _SimpleList     parents,
                        all_but_one (num_nodes-1, 0, 1),
                        aux_list,
                        nk_tuple;

        _Parameter      score;
        _Matrix         single_parent_scores (num_nodes, 1, false, true);

        MPI_Status      status; // contains source, tag, and error code


        MPI_Comm_size (MPI_COMM_WORLD, &size);  // number of processes
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);


        if (rank == 0) {
            _String     bgmSwitch ("_BGM_SWITCH_"),
                        bgmStr;

            _List       * this_list;

            _Matrix     * mpi_node_status = new _Matrix ((long)size, 2, false, true);   // semaphore


            // prepare message exporting _BayesianGraphicalModel object to compute nodes
            SerializeBGMtoMPI (bgmStr);
            ReportWarning (_String("Serialized _BayesianGraphicalModel object as:\n") & bgmStr);


            // switch compute nodes from mpiNormal to bgm loop context
            for (long ni = 1; ni < size; ni++) {
                MPISendString (bgmSwitch, ni);
            }


            // receive confirmation of successful switch
            for (long ni = 1; ni < size; ni++) {
                long fromNode = ni;

                _String t (MPIRecvString (ni,fromNode));
                if (!t.Equal (&bgmSwitch)) {
                    WarnError (_String("Failed to confirm MPI mode switch at node ") & ni);
                    return;
                } else {
                    ReportWarning (_String("Successful mode switch to Bgm MPI confirmed from node ") & ni);
                    MPISendString (bgmStr, ni);
                }
            }


            // farm out jobs to idle nodes until none are left
            for (long node_id = 0; node_id < num_nodes; node_id++) {
                long        maxp            = max_parents.lData[node_id];
                _String     mxString,
                            mxName;
                _Parameter  score;


                this_list   = (_List *) node_score_cache.lData[node_id];

                if (node_type.lData[node_id] == 0) {
                    score = ComputeDiscreteScore (node_id, parents);    // [_SimpleList parents] should always be empty for master node
                } else {
                    score = ComputeContinuousScore (node_id, parents);
                }

                _Constant   orphan_score (score);


                this_list->Clear();
                (*this_list) && (&orphan_score);    // handle orphan score locally


                if (maxp > 0) { // don't bother to farm out trivial cases
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
                        mpi_node++;
                    } while (mpi_node < size);
                }


                if (mpi_node == size) { // all nodes are busy, wait for one to send results
                    MPIReceiveScores (mpi_node_status, true, node_id);
                }
            }
            // end loop over child nodes


            // collect remaining jobs
            while (1) {
                // look for a busy node
                mpi_node = 1;
                do {
                    if ( (*mpi_node_status)(mpi_node,0) == 1) {
                        break;
                    }
                    mpi_node++;
                } while (mpi_node < size);

                if (mpi_node < size) {
                    MPIReceiveScores (mpi_node_status, false, 0);    // at least one node is busy
                } else {
                    break;
                }
            }


            // shut down compute nodes
            for (long shutdown = -1, mpi_node = 1; mpi_node < size; mpi_node++) {
                ReportMPIError(MPI_Send(&shutdown, 1, MPI_LONG, mpi_node, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD), true);
                ReportWarning (_String ("Node 0 sending shutdown signal to node ") & mpi_node);
            }

        }

        else {  // compute node
            long        node_id, maxp;  // update in while loop
            _List       list_of_matrices;

            while (1) {
                _String     mxString,
                            mxName;

                list_of_matrices.Clear();


                // wait for master node to issue node ID to compute
                ReportMPIError (MPI_Recv (&node_id, 1, MPI_LONG, 0, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD, &status), false);
                ReportWarning (_String("Node") & (long)rank & " received child " & (long)node_id & " from node " & (long)status.MPI_SOURCE);

                if (node_id < 0) {
                    ReportWarning (_String("Node") & (long)rank & " recognizes shutdown signal.\n");
                    break;  // received shutdown message (-1)
                }

                maxp = max_parents.lData[node_id];
                parents << 0;   // allocate space for one parent

                if (parents.lLength != 1) { /* DEBUGGING */
                    WarnError (_String("ERROR: _SimpleList parents not expected length (1), found ") & parents.lLength);
                }


                // compute single parent scores
                for (long par = 0; par < num_nodes; par++) {
                    if (par == node_id) {   // child cannot be its own parent
                        continue;
                    }

                    parents.lData[0] = par;

                    if (node_type.lData[node_id] == 0) {
                        // discrete-valued child node cannot have continuous parent
                        if (node_type.lData[par] == 0) {
                            score = ComputeDiscreteScore (node_id, parents);
                        }
                    } else {
                        // continuous child can have discrete or continuous parents
                        score = ComputeContinuousScore (node_id, parents);
                    }

                    single_parent_scores.Store (par, 0, score);
                }


                // compute multiple parent scores
                for (long np = 2; np <= maxp; np++) {
                    _NTupleStorage  family_scores (num_nodes-1, np);

                    parents << 0;   // allocate space for one additional parent

                    if (all_but_one.NChooseKInit (aux_list, nk_tuple, np, false)) {
                        bool    remaining, family_is_discrete;
                        long    res;

                        do {
                            remaining = all_but_one.NChooseK (aux_list, nk_tuple);

                            if (node_type.lData[node_id] == 0) {
                                family_is_discrete = TRUE;
                                for (long par, par_idx = 0; par_idx < np; par_idx++) {
                                    par = nk_tuple.lData[par_idx];
                                    if (par >= node_id) {
                                        par++;
                                    }

                                    if (node_type.lData[par] != 0) {    // all parents must be discrete
                                        family_is_discrete = FALSE;
                                        break;
                                    }

                                    parents.lData[par_idx] = par;
                                }

                                if (family_is_discrete) {
                                    score = ComputeDiscreteScore (node_id, parents);
                                }
                            } else {
                                for (long par_idx = 0; par_idx < np; par_idx++) {
                                    long par = nk_tuple.lData[par_idx];
                                    if (par >= node_id) {
                                        par++;
                                    }
                                    parents.lData[par_idx] = par;
                                }

                                score = ComputeContinuousScore (node_id, parents);
                            }

                            res = family_scores.Store (score, nk_tuple);
                        } while (remaining);
                    } else {
                        _String oops ("Failed to initialize _NTupleStorage object in Bgm::CacheNodeScores().\n");
                        WarnError(oops);
                    }

                    list_of_matrices << (&family_scores);   // append duplicate to storage
                }

                // send results to master node

                ReportMPIError (MPI_Send (single_parent_scores.theData, num_nodes, MPI_DOUBLE, 0, HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD), true);

                for (long np = 2; np <= maxp; np++) {
                    _Matrix * storedMx = (_Matrix *) list_of_matrices.lData[np-2];
                    ReportMPIError (MPI_Send (storedMx->theData, storedMx->GetHDim(), MPI_DOUBLE, 0, HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD), true);
                }
            }
            // end while, received shutdown from master node
        }
        // end else if compute node
    }

    else {  // perform single-threaded
#else

    _SimpleList     parents,
                    all_but_one (num_nodes-1, 0, 1),
                    aux_list,
                    nk_tuple;

    _Matrix         single_parent_scores (num_nodes, 1, false, true);


#if !defined __UNIX__ || defined __HEADLESS__ || defined __HYPHYQT__ || defined __HYPHY_GTK__
    TimerDifferenceFunction(false); // save initial timer; will only update every 1 second
#if !defined __HEADLESS__
    SetStatusLine     (empty,_HYBgm_STATUS_LINE_CACHE, empty, 0, HY_SL_TASK|HY_SL_PERCENT);
#else
    SetStatusLine     (_HYBgm_STATUS_LINE_CACHE);
#endif
    _Parameter  seconds_accumulator = .0,
                temp;
#endif

    for (long node_id = 0; node_id < num_nodes; node_id++) {
        long        maxp        = max_parents.lData[node_id];
        _List   *   this_list   = (_List *) node_score_cache.lData[node_id];    // retrieve pointer to list of scores for this child node

        this_list->Clear();     // reset the list

        // prepare some containers
        _SimpleList parents;
        _Parameter  score;

        if (node_type.lData[node_id] == 0) {    // child node is discrete (multinomial)
            score = ComputeDiscreteScore (node_id, parents);
        } else {                            // child node is continuous (conditional Gaussian)
            score = ComputeContinuousScore (node_id, parents);
        }


        _Constant   orphan_score (score);   // score computed with no parents (initialized empty list)
        (*this_list) && (&orphan_score);

#if !defined __UNIX__ || defined __HEADLESS__
        temp = .0;
#endif

        if (maxp > 0) {
            _Matrix     single_parent_scores (num_nodes, 1, false, true);

            parents << 0;   // allocate space for one parent

            for (long par = 0; par < num_nodes; par++) {
                if (par == node_id) {   // child cannot be its own parent, except morally-challenged time travellers
                    continue;
                }

                parents.lData[0] = par;

                // discrete child cannot have continuous parent
                if (node_type.lData[node_id] == 0) {
                    score = (node_type.lData[par] == 0) ? ComputeDiscreteScore (node_id, parents) : -A_LARGE_NUMBER;
                } else {
                    score = ComputeContinuousScore (node_id, parents);
                }

                single_parent_scores.Store (par, 0, score);
            }
            (*this_list) && (&single_parent_scores);
        }

        for (long np = 2; np <= maxp; np++) {
            _NTupleStorage  family_scores (num_nodes-1, np);

            parents << 0;   // allocate space for one additional parent

            if (all_but_one.NChooseKInit (aux_list, nk_tuple, np, false)) {
                bool    remaining, family_is_discrete;
                long    res;

                do {
                    remaining = all_but_one.NChooseK (aux_list, nk_tuple);

                    if (node_type.lData[node_id] == 0) {
                        family_is_discrete = TRUE;
                        for (long par, par_idx = 0; par_idx < np; par_idx++) {
                            par = nk_tuple.lData[par_idx];
                            if (par >= node_id) {
                                par++;
                            }

                            if (node_type.lData[par] != 0) {    // all parents must be discrete
                                family_is_discrete = FALSE;
                                break;
                            }

                            parents.lData[par_idx] = par;
                        }

                        if (family_is_discrete) {
                            score = ComputeDiscreteScore (node_id, parents);
                        }
                    } else {
                        for (long par_idx = 0; par_idx < np; par_idx++) {
                            long par = nk_tuple.lData[par_idx];
                            if (par >= node_id) {
                                par++;
                            }
                            parents.lData[par_idx] = par;
                        }

                        score = ComputeContinuousScore (node_id, parents);
                    }

                    res = family_scores.Store (score, nk_tuple);
                } while (remaining);
            } else {
                _String oops ("Failed to initialize _NTupleStorage object in Bgm::CacheNodeScores().\n");
                WarnError(oops);
            }
            (*this_list) && (&family_scores);   // append duplicate to storage
        }


#if !defined __UNIX__ || defined __HEADLESS__ || defined __HYPHYQT__ || defined __HYPHY_GTK__
        if ((temp=TimerDifferenceFunction(true))>1.0) { // time to update
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

            if (terminateExecution) { // user wants to cancel the analysis
                break;
            }
        }
#endif
    }
#endif

#if defined __HYPHYMPI__
    }   // completes else statement
#endif

#if !defined __UNIX__ || defined __HEADLESS__
    SetStatusLine     (_HYBgm_STATUS_LINE_CACHE_DONE);
#endif

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
                maxp       = max_parents.lData[this_node];

    _List   *   this_list   = (_List *) node_score_cache.lData[this_node];


    _String     mxString,
                mxName;


    mpi_node_status->Store (senderID, 0, 0);    // set node status to idle
    ReportWarning (_String("Received scores for child ") & this_node & " from node " & senderID);

    (*this_list) && (&single_parent_scores);


    _SimpleList     parents,
                    all_but_one (num_nodes-1, 0, 1),
                    aux_list,
                    nk_tuple;

    _Parameter      score;

    // parse nk-tuple results
    for (long np = 2; np <= maxp; np++) {
        _NTupleStorage  family_scores (num_nodes-1, np);
        bool            remaining;


        if (all_but_one.NChooseKInit (aux_list, nk_tuple, np, false)) {
            long        score_index     = 0,
                        num_ktuples      = exp(lnGamma(num_nodes) - lnGamma(num_nodes - np) - lnGamma(np+1)),
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
            _String oops ("Failed to initialize _NTupleStorage object in Bgm::CacheNodeScores().\n");
            WarnError(oops);
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
void _BayesianGraphicalModel::SerializeBGMtoMPI (_String & rec)
{
    char        buf [255];

    _String *   bgmName = (_String *) bgmNamesList (bgmList._SimpleList::Find((long)this));

    _Parameter  impute_max_steps, impute_burnin, impute_sample_size;    // permit user to reset impute settings

    checkParameter (_HYBgm_IMPUTE_MAXSTEPS, impute_max_steps, 10000);
    checkParameter (_HYBgm_IMPUTE_BURNIN, impute_burnin, 1000);
    checkParameter (_HYBgm_IMPUTE_SAMPLES, impute_sample_size, 100);

    rec << "USE_MPI_CACHING=1;\n";
    rec << "PRINT_DIGITS=-1;\n";

    // write utility functions
    rec << "function add_discrete_node (id,maxp,size,nlev)\n";
    rec << "{\ndnode={};\n";
    rec << "dnode[\"NodeID\"]=id;\n";
    rec << "dnode[\"NodeType\"]=0;\n";
    rec << "dnode[\"MaxParents\"]=maxp;\n";
    rec << "dnode[\"PriorSize\"]=size;\n";
    rec << "dnode[\"NumLevels\"]=nlev;\n";
    rec << "return dnode;\n}\n";

    rec << "function add_gaussian_node (id,maxp,size,mean,prec,scale)\n";
    rec << "{\ngnode={};\n";
    rec << "gnode[\"NodeID\"]=id;\n";
    rec << "gnode[\"NodeType\"]=1;\n";
    rec << "gnode[\"PriorSize\"]=size;\n";
    rec << "gnode[\"MaxParents\"]=maxp;\n";
    rec << "gnode[\"PriorMean\"]=mean;\n";
    rec << "gnode[\"PriorVar\"]=prec;\n";
    rec << "gnode[\"PriorScale\"]=scale;\n";
    rec << "return gnode;\n}\n";


    // prepare assortative lists
    rec << "nodes={};\n";

    for (long node_id = 0; node_id < num_nodes; node_id++) {
        if (node_type.lData[node_id] == 0) {
            rec << "nodes[Abs(nodes)]=add_discrete_node(";
            // I can't figure out how to write the _String object from [node_names] to file :-P
            snprintf (buf, sizeof(buf), "\"Node%d\",%d,%d,%d", node_id, (long)max_parents.lData[node_id], (long)prior_sample_size(node_id,0), num_levels.lData[node_id]);
            rec << buf;
            rec << ");\n";
        } else {
            rec << "nodes[Abs(nodes)]=add_gaussian_node(";
            snprintf (buf, sizeof(buf), "\"Node%d\",%d,%d,%f,%f,%f", node_id, (long)max_parents.lData[node_id], (long)prior_sample_size(node_id,0),
                     prior_mean(node_id,0), prior_precision(node_id,0), prior_scale(node_id,0) );
            rec << buf;
            rec << ");\n";
        }
    }

    // write BGM constructor
    rec << "BayesianGraphicalModel ";
    rec << bgmName;
    rec << "=(nodes);\n";

    // missing data imputation settings
    snprintf (buf, sizeof(buf), "BGM_IMPUTE_MAXSTEPS = %d;\n", (long)impute_max_steps);
    rec << buf;
    snprintf (buf, sizeof(buf), "BGM_IMPUTE_BURNIN = %d;\n", (long)impute_burnin);
    rec << buf;
    snprintf (buf, sizeof(buf), "BGM_IMPUTE_SAMPLES = %d;\n", (long)impute_sample_size);
    rec << buf;

    // serialize constraint matrix
    rec << "bgmConstraints=";
    rec << (_String *)constraint_graph.toStr();
    rec << ";\n";
    rec << "SetParameter(";
    rec << bgmName;
    rec << ",BGM_CONSTRAINT_MATRIX,bgmConstraints);\n";

    // serialize data matrix and assign to BGM
    rec << "bgmData=";
    rec << (_String *)theData.toStr();
    rec << ";\n";
    rec << "SetParameter(";
    rec << bgmName;
    rec << ",BGM_DATA_MATRIX,bgmData);\n";

    /*
    rec << "bgmWeights=";
    rec << (_String *)theWeights.toStr();
    rec << ";\n";
    rec << "SetParameter(";
    rec << bgmName;
    rec << ",BGM_WEIGHT_MATRIX,bgmWeights);\n";
    */
}
#endif



//_________________________________________________________________________________________________________
_Matrix *   _BayesianGraphicalModel::Optimize (void)
{
    /* ---------------------------------------------------------------------------------------------
        OPTIMIZE()
        Wrapper function to access various optimization methods (K2, structural MCMC, order MCMC)
            and implement MPI functionality.
       --------------------------------------------------------------------------------------------- */
    ReportWarning (_String("Entered _BayesianGraphicalModel::Optimize()"));

    if (!scores_cached) {
        CacheNodeScores();
    }


    _Parameter      optMethod;      // 0 = K2 fixed order
    // 1 = K2 shuffle order with restarts
    // 2 = structural MCMC, fixed order
    // 3 = structural MCMC
    // 4 = order MCMC

    _Matrix *   output_matrix;      // for order-MCMC and graph-MCMC has 4 columns: (1) chain trace;
    //  (2) edge marginal posteriors; (3) best state (node ordering or graph);
    //  (4) last state visited in chain


    checkParameter (_HYBgm_METHOD_KEY, optMethod, 0.);

    ReportWarning (_String("... optimization method set to ") & optMethod);

    if (optMethod < 2) {
        ReportWarning (_String("... starting K2 algorithm"));

        _Parameter  num_restarts,           // HBL settings
                    num_randomize;

        checkParameter (_HYBgm_K2_RESTARTS, num_restarts, 1.);
        checkParameter (_HYBgm_K2_RANDOMIZE, num_randomize, num_nodes);

        checkPointer (output_matrix =  new _Matrix (num_nodes * num_nodes, 2, false, true));
        K2Search (optMethod, num_restarts, num_randomize, output_matrix);
    } else {
        _String         oops;
        _Parameter      mcmc_steps, mcmc_burnin, mcmc_samples,
                        mcmc_nchains, mcmc_dtemp;


        // acquisition of HBL arguments with a few sanity checks
        checkParameter (_HYBgm_MCMC_MAXSTEPS, mcmc_steps, 0);
        if (mcmc_steps <= 0)    {
            oops = _String ("You asked HyPhy to run MCMC with zero steps in the chain! Did you forget to set Bgm_MCMC_STEPS?\n");
        }

        checkParameter (_HYBgm_MCMC_BURNIN, mcmc_burnin, 0);
        if (mcmc_burnin < 0)    {
            oops = _String("You can't have a negative burn-in (_HYBgm_MCMC_BURNIN)!\n");
        }

        checkParameter (_HYBgm_MCMC_SAMPLES, mcmc_samples, 0);
        if (mcmc_samples < 0)   {
            oops = _String ("You can't have a negative sample size!");
        }

        if (oops.sLength > 0) {
            WarnError (oops);
            return nil;
        }

#if defined __HYPHYMPI__ && defined __AFYP_DEVELOPMENT__
        int         size,
                    rank;
        long        mpi_node;
        char        mpi_message [256];

        MPI_Status  status; // contains source, tag, and error code

        MPI_Comm_size (MPI_COMM_WORLD, &size);  // number of processes
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);


        checkParameter (_HYBgm_MCMC_NCHAINS, mcmc_nchains, 1);
        if (mcmc_nchains <= 0) {
            oops = _String ("You must run at least one chain in MCMC.");
            WarnError (oops);

            return nil;
        }



        // parallel coupled chains (MPI)
        if (mcmc_nchains > 1) {
            checkParameter (_HYBgm_MCMC_TEMP, mcmc_dtemp, 1);
            if (mcmc_nchains > 0 && mcmc_dtemp <= 1) {
                oops = _String ("Temperature increment must be greater than 1.");
                WarnError (oops);
            }

            _Parameter  chain_T = 1. / (1. + mcmc_dtemp*rank);


            /* note that GraphMetropolis already returns structure and last posterior prob as matrix entries */
            /* use sampling interval to swap chains */
            /* each chain needs a matrix pointer for output */
            if (optMethod < 4) {
                checkPointer (output_matrix = new _Matrix ( (mcmc_samples > num_nodes*num_nodes) ? mcmc_samples : num_nodes*num_nodes, 4, false, true));
                GraphMetropolis (optMethod==2, mcmc_burnin, mcmc_steps, mcmc_samples, chain_T, output_matrix);

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
            checkPointer (output_matrix = new _Matrix ( (mcmc_samples > num_nodes*num_nodes) ? mcmc_samples : num_nodes*num_nodes, 4, false, true));

            if (optMethod < 4) {
                ReportWarning (_String("... starting graph-mcmc"));
                GraphMetropolis ( (optMethod == 2), mcmc_burnin, mcmc_steps, mcmc_samples, 1., output_matrix);
            } else {
                ReportWarning (_String("... starting order-mcmc"));
                if (mcmc_burnin > 0) {
                    ReportWarning (_String("Executing order-MCMC for burn-in period of ") & mcmc_burnin & " steps");
                    OrderMetropolis (FALSE, mcmc_burnin, mcmc_samples, 1., output_matrix);

                    ReportWarning (_String("Automatically reset node_order_arg to best order visited in order-MCMC burn-in:\n "));
                    if (node_order_arg.lLength == 0) {
                        node_order_arg.Populate (num_nodes, 0, 0);
                    }
                    for (long i = 0; i < num_nodes; i++) {
                        node_order_arg.lData[i] = (*output_matrix) (i,3);
                    }
                    ReportWarning    (_String((_String*)node_order_arg.toStr()));
                }
                ReportWarning (_String("Executing order-MCMC for ") & mcmc_steps & " steps, sample size " & mcmc_samples);
                OrderMetropolis (TRUE, mcmc_steps, mcmc_samples, 1., output_matrix);
            }

#if defined __AFYP_DEVELOPMENT__ && defined __HYPHYMPI__
        }
#endif

    }


    return (_Matrix *) output_matrix;
}



//_________________________________________________________________________________________________________
void _BayesianGraphicalModel::K2Search (bool do_permute_order, long n_restart, long n_randomize, _Matrix * result)
{
    //  THIS NEEDS UPDATING
#ifdef __NEVER_DEFINED__
    _Parameter  this_score, best_score, next_score;

    _Matrix     order_matrix (num_nodes, num_nodes, false, true),
                best_dag (num_nodes, num_nodes, false, true);

    _SimpleList best_node_order;



    best_node_order.Populate (num_nodes, 0, 1);
    best_node_order.Permute (1);



    //  Convert node order to binary matrix where edge A->B is permitted if
    //  order_matrix[B][A] = 1, i.e. B is to the right of A in node order.
    for (long i = 0; i < num_nodes; i++) {
        long    child = best_node_order.lData[i];

        for (long j = 0; j < num_nodes; j++) {
            long    parent = best_node_order.lData[j];

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
            } while (num_parents < max_parents.lData[child]);
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
                long    child = best_node_order.lData[i];
                for (long j = 0; j < num_nodes; j++) {
                    long    parent = best_node_order.lData[j];
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
bool    _BayesianGraphicalModel::GraphObeysOrder (_Matrix & graph, _SimpleList & order)
{
    _Matrix order_matrix (num_nodes, num_nodes, false, true);

    // convert order vector to matrix form
    for (long p_index = 0; p_index < num_nodes; p_index++) {
        for (long par = order.lData[p_index], c_index = 0; c_index < num_nodes; c_index++) {
            order_matrix.Store (par, order.lData[c_index], (p_index > c_index) ? 1 : 0);
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
_SimpleList *   _BayesianGraphicalModel::GetOrderFromGraph (_Matrix & graph)
{
    /* ---------------------------------------------------------------------------------
        GetOrderFromGraph()
            To quickly generate a node order based on graph argument.
            Loop through nodes in graph and insert into list according to parentage.
			*** Nodes can only be parents of nodes that appear to their left ***
            For an empty graph, this should return (0,1,2,...,num_nodes-1)
     ----------------------------------------------------------------------------------- */

    _SimpleList *   new_order  = new _SimpleList (1, 0, 0); // initialize with single entry, [0]

    for (long left_of, node = 1; node < num_nodes; node++) {
        // loop through nodes in list looking for parents
        for (left_of = 0; left_of < new_order->lLength; left_of++) {
            // if the node is the child,
            if (graph (left_of, node)) {
                new_order->InsertElement ((BaseRef) node, left_of, false, false);
                break;
            }
        }

        // if we reach end of list, append
        if (left_of == new_order->lLength) {
            (*new_order) << node;
        }

    }
	ReportWarning(_String("Constructed node order from graph:\n") & (_String *)new_order->toStr() & "\n");
    return new_order;
}



//_______________________________________________________________________________________________________________________
void    _BayesianGraphicalModel::GraphMetropolis (bool fixed_order, long mcmc_burnin, long mcmc_steps, long mcmc_samples,
        _Parameter chain_t, _Matrix * result)
{
    /* --------------------------------------------------------------------------------------------------------
        GraphMetropolis()

        Performs MCMC over structures.  Initialize chain using class member [theStructure] (accessible by HBL).
            If theStructure is empty graph, initialize order with random permutation.
            If [fixed_order]=TRUE, only explores the subset of graphs compatible with [node_order_arg] if specified.

        Returns pointer to matrix containing:   (1) vector of posterior probabiities
                                                (2) linearized matrix of edge marginal posteriors
                                                (3)     "       "   for best graph
                                                (4)     "       "   for last graph visited in chain
       --------------------------------------------------------------------------------------------------------- */

    _Matrix     *   proposed_graph  = new _Matrix (num_nodes, num_nodes, false, true);

    _Matrix         current_graph (num_nodes, num_nodes, false, true),
                    best_graph (num_nodes, num_nodes, false, true);

    _Parameter      current_score, proposed_score, best_score,
                    lk_ratio,
                    prob_swap, param_max_fails;

    long            sampling_interval = mcmc_steps / mcmc_samples,
                    max_fails;

    _SimpleList *   proposed_order  = new _SimpleList();
    _SimpleList     current_order;


    // parse HBL settings
    checkParameter (_HYBgm_MCMC_PROBSWAP, prob_swap, 0.1);
    if (prob_swap < 0 || prob_swap > 1.) {
        _String oops ("BGM_MCMC_PROBSWAP must be assigned a value between 0 and 1.  Exiting.\n");
        WarnError (oops);
        return;
    }

    checkParameter (_HYBgm_MCMC_MAXFAILS, param_max_fails, 100);
    if (param_max_fails <= 0.) {
        WarnError ("BGM_MCMC_MAXFAILS must be assigned a value greater than 0");
        return;
    } else {
        max_fails = (long) param_max_fails;
    }




    // initialize chain state
	if (fixed_order)
	{
		if (node_order_arg.lLength > 0 && GraphObeysOrder (theStructure, node_order_arg) )
		{
			(*proposed_graph) = (_Matrix &) theStructure;
            (*proposed_order) = (_SimpleList &) node_order_arg;

            ReportWarning (_String ("Starting GraphMetropolis() using node_order_arg:\n ") & (_String *) proposed_order->toStr());
        } else {
            _String oops ("ERROR: Structure does not match order, aborting GraphMetropolis().");
            WarnError (oops);

            return;
		}
	} else {
		/*
		(*proposed_graph)   = (_Matrix &) theStructure;
		proposed_order      = GetOrderFromGraph (theStructure);
		 */
		proposed_order = GetOrderFromGraph (*proposed_graph);
	}

	
	// randomize the initial graph
	RandomizeGraph (proposed_graph, proposed_order, prob_swap, num_nodes*num_nodes, max_fails, fixed_order);
	ReportWarning (_String ("seeding with randomized graph:\n") & (_String *) proposed_graph->toStr());

    // status line
#if !defined __UNIX__ || defined __HEADLESS__ || defined __HYPHYQT__ || defined __HYPHY_GTK__
    long    updates = 0;
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
    for (long step = 0; step < mcmc_steps + mcmc_burnin; step++) {
        //ReportWarning (_String ("current graph (") & current_score & "): \n" & (_String *) current_graph.toStr());
		
		// modify the graph by one edge
        RandomizeGraph (proposed_graph, proposed_order, prob_swap, 1, max_fails, fixed_order);


        proposed_score = Compute((_Matrix &) *proposed_graph);
        //ReportWarning (_String ("propose graph: (") & proposed_score & "):\n" & (_String *) proposed_graph->toStr());

        lk_ratio = exp(proposed_score - current_score);

        if (lk_ratio > 1. || genrand_real2() < lk_ratio) {  // Metropolis-Hastings
            current_graph   = (_Matrix &) (*proposed_graph);        // accept proposed graph
            current_score   = proposed_score;

            for (long foo = 0; foo < num_nodes; foo++) {
                current_order.lData[foo] = proposed_order->lData[foo];
            }

            if (current_score > best_score) {
                best_score = current_score;
                best_graph = (_Matrix &) current_graph;         // keep track of best graph ever visited by chain
            }
        } else {
            // revert proposal
            for (long row = 0; row < num_nodes; row++) {
                proposed_order->lData[row] = current_order.lData[row];

                for (long col = 0; col < num_nodes; col++) {
                    proposed_graph->Store(row, col, current_graph(row, col));
                }
            }
        }

        // chain sampling
        if (step >= mcmc_burnin) {
            if ( (step-(long)mcmc_burnin) % sampling_interval == 0) {
                long    entry   = (long int) ((step-mcmc_burnin) / sampling_interval);

                result->Store (entry, 0, current_score);        // update chain trace

                for (long row = 0; row < num_nodes; row++) {    // update edge tallies
                    for (long offset=row*num_nodes, col = 0; col < num_nodes; col++) {
                        result->Store (offset+col, 1, (*result)(offset+col,1) + current_graph(row, col));
                    }
                }
            }
        }

#if !defined __UNIX__ || defined __HEADLESS__
        if (TimerDifferenceFunction(true)>1.0) { // time to update
            updates ++;
            _String statusLine = _HYBgm_STATUS_LINE_MCMC & " " & (step+1) & "/" & (mcmc_steps + mcmc_burnin)
                                 & " steps (" & (1.0+step)/updates & "/second)";
#if defined __HEADLESS__
            SetStatusLine     (statusLine);
#else
            SetStatusLine     (empty,statusLine,empty,100*step/(mcmc_steps + mcmc_burnin),HY_SL_TASK|HY_SL_PERCENT);
#endif
            TimerDifferenceFunction(false); // reset timer for the next second
            yieldCPUTime(); // let the GUI handle user actions
            if (terminateExecution) { // user wants to cancel the analysis
                break;
            }
        }
#else
        if (step % ((mcmc_steps + mcmc_burnin)/100) == 0) {
            ReportWarning (_String ("graphMCMC at step ") & step & " of " & (mcmc_steps+mcmc_burnin));
        }
#endif

    }

    // convert edge tallies to frequencies, and report best and last graph
    for (long row = 0; row < num_nodes; row++) {
        for (long offset=row*num_nodes, col = 0; col < num_nodes; col++) {
            result->Store (offset+col, 1, (*result)(offset+col,1)/mcmc_samples);
            result->Store (offset+col, 2, (_Parameter) best_graph(row, col));
            result->Store (offset+col, 3, (_Parameter) current_graph(row,col));
        }
    }


    // update theStructure with last graph visited by chain
    theStructure = (_Matrix &) current_graph;
    ReportWarning (_String("On exiting GraphMetropolic() assigned last state to theStructure: ") & (_String *) theStructure.toStr());

    // release memory
    delete (proposed_graph);
    delete (proposed_order);
}



//_____________________________________________________________________________________________________________
void    _BayesianGraphicalModel::RandomizeGraph (_Matrix * graph, _SimpleList * order, _Parameter prob_swap, long num_steps, long max_fails, bool fixed_order)
{
    //  Modify a graph by adding, removing, or reversing an edge, so long as that modification
    //  complies with the constraint matrix and (when fixed_order=TRUE) node order.

    long    step = 0, fail = 0;


    // calculate number of parents for each child for rapid access later
    _SimpleList     num_parents;

    num_parents.Populate (num_nodes, 0, 0);

    for (long child = 0; child < num_nodes; child++) {
        for (long parent = 0; parent < num_nodes; parent++) {
            if ( (*graph)(parent, child) > 0) {
                num_parents.lData[child]++;
            }
        }

        if (num_parents.lData[child] > max_parents.lData[child]) {
            WarnError (_String ("Number of parents exceeds maximum BEFORE randomization of graph at node ") & child & " (" & num_parents.lData[child] & " > " & max_parents.lData[child] & ")\n" );
            break;
        }
    }


    // randomize graph
    long    p_rank, c_rank, parent, child;

    do {
        // a fail-safe to avoid infinite loops
        if (fail > max_fails) {
            WarnError (_String ("Failed to modify the graph in RandomizeGraph() after ") & max_fails & " attempts.");
            break;
        }

        // pick a random edge
        p_rank  = (genrand_int32() % (num_nodes-1)) + 1;    // shift right to not include lowest node in order
        c_rank  = genrand_int32() % p_rank;

        child = order->lData[c_rank];
        parent = order->lData[p_rank];


        if (fixed_order || genrand_real2() > prob_swap) {
            if ( (*graph)(parent,child) == 0 && constraint_graph(parent,child) >= 0) {
				// add an edge if it is absent and not banned
                if (num_parents.lData[child] == max_parents.lData[child]) {
					// child is already at maximum number of parents
                    // move an edge from an existing parent to target parent
                    _SimpleList     removeable_edges;

                    // build a list of current parents
                    for (long par = 0; par < num_nodes; par++)
						// par is parent of child AND the edge is NOT enforced
                        if ((*graph)(par,child) && constraint_graph(par,child)<=0) {
                            removeable_edges << par;
                        }

                    if (removeable_edges.lLength > 0) {
                        // shuffle the list and remove the first parent
						removeable_edges.Permute(1);
						graph->Store (removeable_edges.lData[0], child, 0.);
						graph->Store (parent, child, 1.);
						
						/*
                        graph->Store (removeable_edges.lData[ (long) (genrand_real2()*removeable_edges.lLength) ], child, 0.);
                        graph->Store (parent, child, 1.);
						 */
                        step++;
                    } else {
                        // none of the edges can be removed
                        fail++;
                    }
                } else {
					// child can accept one more edge
                    graph->Store (parent, child, 1.);
                    num_parents.lData[child]++;
                    step++;
                }
            }

            else if ( (*graph)(parent,child) == 1 && constraint_graph(parent,child) <= 0) {
				// delete an edge if it is present and not enforced
                graph->Store (parent, child, 0.);
                num_parents.lData[child]--;
                step++;
            } else {
                fail++;
            }
        } else {    // swap nodes in ordering and flip edge if present
            bool    ok_to_go = TRUE;

            // P->C cannot be flipped if C->P is banned or P->C is enforced
            if ( (*graph)(parent,child) == 1  &&  ( constraint_graph(child,parent)<0 || constraint_graph(parent,child)>0 )  ) {
                ok_to_go = FALSE;
            }


            // check all other nodes affected by the swap
            if (ok_to_go) {
                for (long bystander, i=c_rank+1; i < p_rank; i++) {
                    bystander = order->lData[i];    // retrieve node id
					
					// by flipping the parent->child edge, we would screw up ordering for
					// an enforced edge:  C < B < P   becomes  P < B < C  where B->C or P->B is enforced.
					// also, a banned edge implies a node ordering constraint
                    if (
                        ( (*graph)(parent,bystander)==1 && constraint_graph (parent, bystander)>0 )  ||
                        ( (*graph)(bystander,child)==1 && constraint_graph (bystander, child)>0 )  ||
                        constraint_graph (bystander,parent)<0  ||  
						constraint_graph (child,bystander)<0
                    ) {
                        ok_to_go = FALSE;
                        break;
                    }
                }
            }


            // if everything checks out OK
            if ( ok_to_go ) {
                if ( (*graph)(parent,child) == 1 && constraint_graph(child,parent)>=0 ) {
					// flip the target edge
                    graph->Store (parent, child, 0);
                    graph->Store (child, parent, 1);
                    num_parents.lData[child]--;
                    num_parents.lData[parent]++;

                    if (num_parents.lData[parent] > max_parents.lData[parent]) {
                        // parent cannot accept any more edges, delete one of the edges at random (including the edge to flip)
                        _SimpleList     removeable_edges;

                        for (long par = 0; par < num_nodes; par++)
                            if (  (*graph)(par, parent)  &&  constraint_graph(par, parent)<=0  ) {
                                removeable_edges << par;
                            }

						removeable_edges.Permute(1);
						graph->Store (removeable_edges.lData[0], parent, 0.);
						/*
                        graph->Store (removeable_edges.lData[ (long) (genrand_real2()*removeable_edges.lLength) ], parent, 0.);
						 */
                        num_parents.lData[parent]--;
                    }

                    step++;
                }
				// if the number of parents for parent node will exceed maximum, 
				// then the edge is deleted instead of flipped (delete without add)

                // swap nodes in order
                order->lData[p_rank] = child;
                order->lData[c_rank] = parent;	// remember to update the order matrix also!

                // update edges affected by node swap
				//  child <-  N   ...   N <- parent                   _________________,
                //                                    becomes        |                 v
                //                                                 parent    N         N    child
                //                                                           ^                |
                //                                                           `----------------+
                for (long bystander, i = c_rank+1; i < p_rank; i++) {
                    bystander = order->lData[i];

                    if ( (*graph)(bystander, child) == 1 ) {
                        graph->Store (bystander, child, 0);
                        num_parents.lData[child]--;

                        graph->Store (child, bystander, 1);
                        num_parents.lData[bystander]++;

                        if (num_parents.lData[bystander] > max_parents.lData[bystander]) {
                            // number of parents for bystander node now exceeds maximum,
                            //  select a random edge to remove - including the edge C->B we just made!
                            _SimpleList     removeable_edges;

                            for (long par = 0; par < num_nodes; par++) {
                                if (  (*graph)(par, bystander)  &&  constraint_graph(par, bystander)<=0  ) {
                                    removeable_edges << par;
                                }
							}
							
							
                            removeable_edges.Permute(1);
							graph->Store (removeable_edges.lData[0], bystander, 0.);
							/*
                            graph->Store (removeable_edges.lData[ (long) (genrand_real2()*removeable_edges.lLength) ], bystander, 0.);
							 */
                            num_parents.lData[bystander]--;
                        }
                    }

                    if ( (*graph)(parent,bystander) == 1) {
						// flip edge from parent to bystander (P->B)
                        graph->Store (parent, bystander, 0);
                        num_parents.lData[bystander]--;

                        graph->Store (bystander, parent, 1);
                        num_parents.lData[parent]++;

                        if (num_parents.lData[parent] > max_parents.lData[parent]) {
							// remove excess edge from X to parent other than bystander
                            _SimpleList     removeable_edges;

                            for (long par = 0; par < num_nodes; par++) {
                                if (  (*graph)(par, parent)  &&  constraint_graph(par, parent)<=0  ) {
                                    removeable_edges << par;
                                }
							}
							
                            removeable_edges.Permute(1);
							graph->Store (removeable_edges.lData[0], parent, 0.);
							/*
                            graph->Store (removeable_edges.lData[ (long) (genrand_real2()*removeable_edges.lLength) ], parent, 0.);
							 */
                            num_parents.lData[parent]--;
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



//_______________________________________________________________________________________________________________________
void    _BayesianGraphicalModel::OrderMetropolis (bool do_sampling, long n_steps, long sample_size, _Parameter chain_t,
        _Matrix * result)
{
    /* ----------------------------------------------------------------------------------------
        OrderMetropolis()

        Execute Metropolis sampler using a swap of two nodes in an ordered sequence
        as a proposal function.  The posterior probabilities of edges in the network are
        stored as a member matrix [edge_posteriors].  Note that the total number of
        possible orderings (i.e. permutations of a sequence of length N) is factorial(N),
        and can possibly be computed exactly for N < 8.
       ---------------------------------------------------------------------------------------- */
    long            first_node, second_node,
                    sample_lag = n_steps / sample_size;

    _Parameter      lk_ratio,
                    prob_current_order, prob_proposed_order, best_prob,
                    denom;

    _SimpleList     current_order, proposed_order, best_node_order,
                    * ptr_to_order;

    _List           * marginals         = new _List ();

    _GrowingVector  * gv;


    InitMarginalVectors (marginals);    // allocate storage



    /* SLKP 20070926
     Add user feedback via the console window status bar
     */
    VerbosityLevel();
    long howManyTimesUpdated = 0; // how many times has the line been updated; is the same as the # of seconds
    TimerDifferenceFunction(false); // save initial timer; will only update every 1 second
#ifdef __HEADLESS__
    SetStatusLine (_HYBgm_STATUS_LINE_MCMC & (do_sampling ? empty : _String(" burnin")));
#else
#ifndef __UNIX__
    SetStatusLine     (empty,_HYBgm_STATUS_LINE_MCMC & (do_sampling ? empty : _String(" burnin")), empty ,0,HY_SL_TASK|HY_SL_PERCENT);
#else
    _String         * progressReportFile = NULL;
    _Variable       * progressFile = CheckReceptacle (&optimizationStatusFile, empty, false);

    if (progressFile->ObjectClass () == STRING) {
        progressReportFile = ((_FString*)progressFile->Compute())->theString;
    }
    ConsoleBGMStatus (_HYBgm_STATUS_LINE_MCMC & (do_sampling ? empty : _String(" burnin")), -1., progressReportFile);
#endif
#endif

    /* SLKP */



    // initialize node ordering
    if (node_order_arg.lLength > 0) {
        current_order = node_order_arg;
    } else {
        ptr_to_order = GetOrderFromGraph (theStructure);
        current_order.Duplicate(ptr_to_order);
        DeleteObject (ptr_to_order);
    }

    best_prob = prob_current_order = Compute (current_order, marginals);
    best_node_order = current_order;

    proposed_order.Populate (num_nodes, 0, 1);

    // chain
    for (long step = 0; step < n_steps; step++) {
        // copy over current order to proposed order
        for (long i = 0; i < proposed_order.lLength; i++) {
            proposed_order.lData[i] = current_order.lData[i];
        }


        // swap random nodes in ordered sequence
        first_node  = genrand_int32() % num_nodes;
        second_node = genrand_int32() % num_nodes;

        while (first_node == second_node) {
            second_node = genrand_int32() % num_nodes;
        }

        proposed_order.Swap (first_node, second_node);


        // compute likelihood ratio of proposed : current orders
        prob_proposed_order = Compute (proposed_order, marginals);
        lk_ratio            = exp(prob_proposed_order - prob_current_order);

        if (lk_ratio > 1. || genrand_real2() < lk_ratio) {  // then set current to proposed order
            current_order = proposed_order;
            prob_current_order = prob_proposed_order;

            if (prob_proposed_order > best_prob) {
                best_prob = prob_proposed_order;        // update best node ordering
                best_node_order = proposed_order;
            }
        }


        // if past burn-in period and at sampling step, record trace and marginals
        if (do_sampling) {
            if (step % sample_lag == 0) {
				ReportWarning(_String("At step ") & step & " order: " & (_String *) current_order.toStr());
				
                result->Store (step / sample_lag, 0, prob_current_order);

                for (long child = 0; child < num_nodes; child++) {
                    // retrieve information from Compute()
                    gv      = (_GrowingVector *) marginals->lData[child * num_nodes + child];
                    denom   = (*gv)(0, 0);  // first entry holds node marginal post pr

                    for (long edge, parent = 0; parent < num_nodes; parent++) {
                        if (parent == child) {
                            continue;
                        }

                        edge    = child * num_nodes + parent;
                        gv      = (_GrowingVector *) marginals->lData[edge];

                        // not all _GrowingVector entries in marginals are being used,
                        //      i.e., edges incompatible with order
                        if (gv->GetUsed() > 0) {
							// store as transpose so that row = parent and column = child, same as GraphMCMC
                            result->Store (parent*num_nodes+child, 1, (*result)(parent*num_nodes+child, 1) + exp (LogSumExpo(gv) - denom));
                        }
                    }
                }
            }
        }


        /*SLKP 20070926; include progress report updates */
        if (TimerDifferenceFunction(true)>1.0) { // time to update
            howManyTimesUpdated ++;
            _String statusLine = _HYBgm_STATUS_LINE_MCMC & (do_sampling ? empty : _String(" burnin")) & " " & (step+1) & "/" & n_steps
                                 & " steps (" & (1.0+step)/howManyTimesUpdated & "/second)";
#if  defined __HEADLESS__
            SetStatusLine (statusLine);
#else
#if !defined __UNIX__
            SetStatusLine     (empty,statusLine,empty,100*step/(n_steps),HY_SL_TASK|HY_SL_PERCENT);
            yieldCPUTime(); // let the GUI handle user actions
            if (terminateExecution) { // user wants to cancel the analysis
                break;
            }
#endif
#if defined __UNIX__ && ! defined __HEADLESS__ 
            ConsoleBGMStatus (statusLine, 100.*step/(n_steps), progressReportFile);
#endif
#endif

            TimerDifferenceFunction(false); // reset timer for the next second
        }
        /* SLKP */

    }
    // chain terminates


    // convert sums of edge posterior probs. in container to means
    for (long edge = 0; edge < num_nodes * num_nodes; edge++) {
        result->Store (edge, 1, (*result)(edge,1) / sample_size);
    }


    // export node ordering info
    for (long node = 0; node < num_nodes; node++) {
        result->Store (node, 2, (_Parameter) (best_node_order.lData[node]));
        result->Store (node, 3, (_Parameter) (current_order.lData[node]));
    }


    DumpMarginalVectors (marginals);


    /*SLKP 20070926; include progress report updates */
#if !defined __UNIX__ || defined __HEADLESS__ || defined __HYPHYQT__ || defined __HYPHY_GTK__
    SetStatusLine     (_HYBgm_STATUS_LINE_MCMC_DONE);
#endif
#if defined __UNIX__ && ! defined __HEADLESS__ && !defined __HYPHYQT__ && !defined __HYPHY_GTK__
    ConsoleBGMStatus (_HYBgm_STATUS_LINE_MCMC_DONE, -1.0, progressReportFile);
#endif
    /* SLKP */


    // update node order
    node_order_arg = current_order;
    ReportWarning (_String("Set node_order_arg to last order visited in orderMCMC:\n") & (_String*)node_order_arg.toStr());
}






