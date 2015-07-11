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


#if not defined __AFYP_REWRITE_BGM__

#ifdef __HYPHYQT__
    #include "hyphymain.h"
#endif

#if defined __MAC__ || defined __WINDOZE__ || defined __HYPHY_GTK__
    #include "HYConsoleWindow.h"
    #include "HYDialogs.h"

#endif


//#define       __AFYP_DEVELOPMENT__
#define     __MISSING_DATA__


#include "bgm.h"

_String     _HYBgm_BAN_PARENT_KEY   ("BanParent"),
            _HYBgm_BAN_CHILD_KEY  ("BanChild"),
            _HYBgm_ENFORCE_PARENT_KEY ("EnforceParent"),
            _HYBgm_ENFORCE_CHILD_KEY  ("EnforceChild"),

            _HYBgm_NODE_INDEX_KEY ("NodeID"),
            _HYBgm_PRIOR_SIZE_KEY ("PriorSize"),
            _HYBgm_MAX_PARENT_KEY ("MaxParents"),

            _HYBgm_PRIOR_MEAN_KEY ("PriorMean"),      /* for continuous (Gaussian) nodes */
            _HYBgm_PRIOR_PREC_KEY ("PriorVar"),

            /*SLKP 20070926; add string constants for progress report updates */
            _HYBgm_STATUS_LINE_MCMC           ("Running Bgm MCMC"),
            _HYBgm_STATUS_LINE_MCMC_DONE  ("Finished Bgm MCMC"),
            _HYBgm_STATUS_LINE_CACHE      ("Caching Bgm scores"),
            _HYBgm_STATUS_LINE_CACHE_DONE ("Done caching Bgm scores"),
            /*SLKP*/

            _HYBgm_MPI_CACHING ("USE_MPI_CACHING"),

            /* maxParentString ("Bgm_MAXIMUM_PARENTS"), */
            maxNumRestart ("BGM_NUM_RESTARTS"),
            numRandomize  ("BGM_NUM_RANDOMIZE"),
            useNodeOrder  ("BGM_USE_NODEORDER"),
            bgmOptimizationMethod ("BGM_OPTIMIZATION_METHOD"),

            useMPIcaching ("USE_MPI_CACHING"),

            mcmcNumChains ("BGM_MCMC_NCHAINS"),       // re-use parameter
            mcmcTemperature ("BGM_MCMC_TEMPERATURE"),
            mcmcSteps     ("BGM_MCMC_DURATION"),
            mcmcBurnin        ("BGM_MCMC_BURNIN"),
            mcmcSamples       ("BGM_MCMC_SAMPLES"),

            mcemMaxSteps  ("BGM_MCEM_MAXSTEPS"),
            mcemBurnin        ("BGM_MCEM_BURNIN"),
            mcemSampleSize    ("BGM_MCEM_SAMPLES");


#ifdef      __UNIX__

void        ConsoleBGMStatus (_String, _Parameter, _String * fileName);

//__________________________________________________________________________________

void        ConsoleBGMStatus (_String statusLine, _Parameter percentDone, _String * fileName = nil)
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


//___________________________________________________________________________________________

long integerPower (long base, long exponent)
{
    long    result = 1L,
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



//___________________________________________________________________________________________
//  afyp Oct. 15, 2008, replaced wrapper function for _Constant member function
//  with this local function.
_Parameter Bgm::LnGamma(_Parameter theValue)
{
    if (theValue <= 0) {


        _String oops ("ERROR: Requested Bgm::LnGamma(x) for x <= 0.");
        WarnError (oops);

        return 0;
    }

    static _Parameter lngammaCoeff [6] = {   76.18009172947146,
                                         -86.50532032941677,
                                         24.01409824083091,
                                         - 1.231739572450155,
                                         0.1208650973866179e-2,
                                         - 0.5395239384953e-5
                                         };

    static _Parameter lookUpTable [20] = {  0., 0., 0.6931472, 1.7917595, 3.1780538,
                                            4.7874917, 6.5792512, 8.5251614, 10.6046029, 12.8018275,
                                            15.1044126, 17.5023078, 19.9872145, 22.5521639, 25.1912212,
                                            27.8992714, 30.6718601, 33.5050735, 36.3954452, 39.3398842
                                         };

    // use look-up table for small integer values
    if (theValue <= 20 && (theValue - (long)theValue) > 0.) {
        return (lookUpTable [(long) theValue - 1]);
    }


    // else do it the hard way
    _Parameter  x, y, tmp, ser;

    y = x = theValue;
    tmp = x + 5.5;
    tmp -= (x+0.5) * log(tmp);
    ser = 1.000000000190015;

    for (long j = 0; j <= 5; j++) {
        ser += lngammaCoeff[j] / ++y;
    }

    return (-tmp + log(2.506628274631005*ser/x));
}



//___________________________________________________________________________________________
Bgm::Bgm(_AssociativeList * dnodes, _AssociativeList * cnodes)
{
    _String     errorMessage;

    long        num_discrete    = dnodes->avl.countitems(),
                num_continuous   = cnodes->avl.countitems(),
                node_index;

    num_nodes   = num_discrete + num_continuous;    // set member variable



    _Constant   * node_id,                      /* pointers into HBL associative array arguments */
                * mp,
                * size, * mean, * precision;

    _AssociativeList    * discreteNode, * continuousNode;



    // set member variables to default values
    calc_bit            = new _Constant ();         // for calculating LnGamma in member functions
    max_max_parents     = 0;                        // record the most parents a node can have
    obsData             = NULL;
    obsWeights          = NULL;



    // extract number of variables in data
    CreateMatrix (&dag, num_nodes, num_nodes, false, true, false);      // allocate space for matrices
    CreateMatrix (&banned_edges, num_nodes, num_nodes, false, true, false);
    CreateMatrix (&enforced_edges, num_nodes, num_nodes, false, true, false);
    CreateMatrix (&prior_sample_size, num_nodes, 1, false, true, false);    // for storing floats as _Parameter objects
    CreateMatrix (&prior_mean, num_nodes, 1, false, true, false);
    CreateMatrix (&prior_precision, num_nodes, 1, false, true, false);


    // allocate space for _SimpleList objects
    is_discrete.Populate (num_nodes, 0, 0);
    num_levels.Populate (num_nodes, 0, 0);
    max_parents.Populate (num_nodes, 0, 0);
    has_missing.Populate (num_nodes, 0, 0);



    // initialize discrete nodes
    for (long dn = 0; dn < num_discrete; dn++) {
        discreteNode    = (_AssociativeList *) (dnodes->GetByKey (dn, ASSOCIATIVE_LIST));

        if (discreteNode) {
            node_id     = (_Constant *) (discreteNode->GetByKey (_HYBgm_NODE_INDEX_KEY, NUMBER));
            mp          = (_Constant *) (discreteNode->GetByKey (_HYBgm_MAX_PARENT_KEY, NUMBER));
            size        = (_Constant *) (discreteNode->GetByKey (_HYBgm_PRIOR_SIZE_KEY, NUMBER));

            if (node_id && mp && size) {
                node_index  = (long) (node_id->Value());

                is_discrete.lData[node_index] = 1.;
                max_parents.lData[node_index] = (long) mp->Value();

                if ((long) mp->Value() > max_max_parents) {
                    max_max_parents = (long) mp->Value();
                }

                prior_sample_size.Store(node_index, 0, (long) size->Value());
            } else {
                errorMessage = _String ("Missing key (NodeID, MaxParents, PriorSize) in associative array for discrete node ")
                               & dn;
                break;
            }
        } else {
            errorMessage = _String("Failed to retrieve discrete node specification at index ") & dn;
            break;
        }
    }
    if (errorMessage.sLength) {
        WarnError (errorMessage);
        errorMessage = _String();   // reset to empty string
    }



    // initialize continuous nodes
    for (long cn = 0; cn < num_continuous; cn++) {
        continuousNode  = (_AssociativeList *) (cnodes->GetByKey (cn, ASSOCIATIVE_LIST));

        if (continuousNode) {
            node_id     = (_Constant *) (continuousNode->GetByKey (_HYBgm_NODE_INDEX_KEY, NUMBER));
            mp          = (_Constant *) (continuousNode->GetByKey (_HYBgm_MAX_PARENT_KEY, NUMBER));
            size        = (_Constant *) (continuousNode->GetByKey (_HYBgm_PRIOR_SIZE_KEY, NUMBER));
            mean        = (_Constant *) (continuousNode->GetByKey (_HYBgm_PRIOR_MEAN_KEY, NUMBER));
            precision   = (_Constant *) (continuousNode->GetByKey (_HYBgm_PRIOR_PREC_KEY, NUMBER));

            if (node_id && mp && size && mean && precision) {
                node_index = (long) (node_id->Value());

                is_discrete.lData[node_index] = 0.;
                max_parents.lData[node_index] = (long) mp->Value();

                if (mp->Value() > max_max_parents) {
                    max_max_parents = (long) mp->Value();
                }

                prior_sample_size.Store (node_index, 0, (_Parameter) size->Value());
                prior_mean.Store (node_index, 0, (_Parameter) mean->Value());
                prior_precision.Store (node_index, 0, (_Parameter) precision->Value());

            } else {
                errorMessage = _String ("Missing key (NodeID, MaxParents, PriorSize, PriorMean, PriorVar)")
                               & "in associative array for continuous node " & cn;
                break;
            }
        } else {
            errorMessage = _String ("Failed to retrieve continuous node specification at index ") & cn;
            break;
        }
    }

    if (errorMessage.sLength) {
        WarnError (errorMessage);
        errorMessage = _String();   // reset to empty string
    }



    // append duplicates of _List object to store pointers to _Constant, _Matrix, _NTupleStorage objects
    _List       emptyList (max_max_parents+1);
    for (long node = 0; node < num_nodes; node++) {
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
    scores_cached   = FALSE;

    if (obsData)    {
        DeleteObject (obsData);    // dispense with previous data matrix
    }

    obsData         = data;

    obsData->CheckIfSparseEnough(TRUE);     // check if we can use a more compact representation of the
    // matrix; before including this command, we suffered a
    // noticeable slow-down in this routine. - AFYP

    /*
    // allocate space to weight matrix
    if (obsWeights)
    {
        DeleteObject(obsWeights);
    }
    obsWeights = new _Matrix (obsData->GetHDim(), obsData->GetVDim(), false, true);

    // by default all cases are weighted equally - require user to call SetWeightMatrix() if otherwise
    for (long row = 0; row < obsData->GetHDim(); row++)
        for (long col = 0; col < obsData->GetVDim(); col++)
            obsWeights->Store(row, col, 1.);

    */

    // reset data-dependent member variables
    for (long node = 0; node < num_nodes; node++) {
        has_missing.lData[node] = 0;
        num_levels.lData[node] = 0;
    }



#ifdef DEBUG_SDM
    char buf [256];
    for (long i = 0; i < obsData->GetHDim(); i++) {
        for (long j = 0; j < obsData->GetVDim(); j++) {
            snprintf (buf, sizeof(buf), "%d ", (long) (*obsData)(i,j));
            BufferToConsole (buf);
        }
        NLToConsole ();
    }
#endif



    if (obsData->GetVDim() == num_nodes) {  // make sure the data matrix is compatible with graph
#ifdef __MISSING_DATA__
        long    nrows = obsData->GetHDim();

        for (long node = 0; node < num_nodes; node++) {
            if (is_discrete.lData[node]) {
                num_levels.lData[node] = 1;

                for (long obs, row = 0; row < nrows; row++) {
                    obs = (*obsData)(row, node);

                    if (has_missing.lData[node] == 0 && obs < 0) {  // use negative integer values to annotate missing data
                        has_missing.lData[node] = 1;
                        continue;   // skip next step to check levels
                    }

                    if (obs+1 > num_levels.lData[node]) {
                        num_levels.lData[node] = num_levels.lData[node] + 1;
                    }
                }
            } else {
                num_levels.lData[node] = 0;

                // not implementing missing data for continuous nodes yet
                //  not until I've decided on annotation anyhow :-/  afyp
            }
        }

#ifdef DEBUG_SDM
        snprintf (buf, sizeof(buf), "Levels: ");
        BufferToConsole (buf);
        for (long i = 0; i < num_nodes; i++) {
            snprintf (buf, sizeof(buf), "%d ", num_levels.lData[i]);
            BufferToConsole (buf);
        }
        NLToConsole ();

        snprintf (buf, sizeof(buf), "Missing (0=FALSE, 1=TRUE): ");
        BufferToConsole (buf);
        for (long i = 0; i < num_nodes; i++) {
            snprintf (buf, sizeof(buf), "%d ", has_missing.lData[i]);
            BufferToConsole (buf);
        }
        NLToConsole ();
#endif

#else
        for (long node = 0; node < num_nodes; node++) { // for every column in matrix
            if (is_discrete.lData[node]) {
                // calculate number of levels for discrete node
                num_levels.lData[node] = 1;

                for (long obs = 0; obs < obsData->GetHDim(); obs++) {
                    // adjust for zero-indexing
                    if ( ((*obsData)(obs,node) + 1) > num_levels.lData[node]) {
                        num_levels.lData[node] = num_levels.lData[node] + 1;
                    }
                }
            } else {
                // continuous node defined to have no levels
                num_levels.lData[node] = 0;
            }
        }
#endif
    } else {
        _String errorMsg ("Number of variables in data matrix do not match number of nodes in graph.");
        WarnError (errorMsg);       // note, this produces a crash because the batch file proceeds to execute
        // BGM routines without data. - AFYP
    }


    last_node_order.Clear();        // forget the last step taken in MCMC chain
    best_node_order.Clear();

    CacheNodeScores();

    // re-allocate memory to lists
    // InitComputeLists ();
}



//___________________________________________________________________________________________
void    Bgm::SetWeightMatrix (_Matrix * weights)
{
    if ( obsData
            && (long) (obsData->GetHDim()) == (long) (weights->GetHDim())
            /*&& (long) (obsData->GetVDim()) == (long) (weights->GetVDim())*/ ) {
        obsWeights = weights;
        ReportWarning( _String("Set weight matrix to ") & (_String *) obsWeights->toStr() );
    } else {
        _String errorMsg ("Number of weights does not match number of observations in current data set.");
        WarnError (errorMsg);
    }

    scores_cached = FALSE;
}




//___________________________________________________________________________________________
void    Bgm::SetGraphMatrix (_Matrix *graph)
{
    /*
    char    bug [255];

    snprintf (bug, sizeof(bug), "Entered Bgm::SetGraphMatrix()\n");
    BufferToConsole (bug);
    */
    dag = (_Matrix &) (*graph); // matrix assignment
    ReportWarning (_String("set graph matrix to:\n") & (_String *) dag.toStr() & "\n" );
}


void    Bgm::SetBanMatrix (_Matrix *banMx)
{
    banned_edges = (_Matrix &) (*banMx);

    ReportWarning (_String("Set ban matrix to ") & (_String *) banned_edges.toStr() & "\n" );
}

void    Bgm::SetEnforceMatrix (_Matrix *enforceMx)
{
    enforced_edges = (_Matrix &) (*enforceMx);
    ReportWarning (_String("Set enforce matrix to ") & (_String *) enforced_edges.toStr() & "\n" );
}

void    Bgm::SetBestOrder (_SimpleList * orderList)
{
    best_node_order.Populate(num_nodes, 0, 0);

    for (long i = 0; i < num_nodes; i++) {
        best_node_order.lData[i] = orderList->lData[i];
    }
}

//___________________________________________________________________________________________
//  For debugging..
void Bgm::PrintGraph (_Matrix * g)
{
    char    buf [255];

    if (g) {
        for (long row = 0; row < g->GetHDim(); row++) {
            for (long col = 0; col < g->GetVDim(); col++) {
                snprintf (buf, sizeof(buf), "%ld ", (long) (*g)(row,col));
                BufferToConsole (buf);
            }
            snprintf (buf, sizeof(buf), "\n");
            BufferToConsole (buf);
        }
    } else {
        for (long row = 0; row < dag.GetHDim(); row++) {
            for (long col = 0; col < dag.GetVDim(); col++) {
                snprintf (buf, sizeof(buf), "%ld ", (long)dag(row,col));
                BufferToConsole (buf);
            }
            snprintf (buf, sizeof(buf), "\n");
            BufferToConsole (buf);
        }
    }
    NLToConsole();
}


//___________________________________________________________________________________________
void    Bgm::ResetGraph (_Matrix * g)
{
    if (g) {
        for (long row = 0; row < num_nodes; row++) {
            for (long col = 0; col < num_nodes; col++) {
                g->Store (row, col, dag (row, col));
            }
        }
        ReportWarning (_String("Reset graph to dag argument\n"));
    } else {
        for (long row = 0; row < num_nodes; row++) {
            for (long col = 0; col < num_nodes; col++) {
                // enforced edges must always appear in DAG
                // if none specified, all entries are zeroed
                dag.Store (row, col, enforced_edges (row, col));
            }
        }
    }
}



//___________________________________________________________________________________________
long    Bgm::MarginalSum (bool over_rows, long offset)
{
    long    res     = 0;

    if (over_rows)
        for (long i = 0; i < num_nodes; i++) {
            res += dag(offset,i);    // sum of n-th row
        }
    else
        for (long i = 0; i < num_nodes; i++) {
            res += dag(i,offset);    // sum of n-th column
        }

    return res;
}



//___________________________________________________________________________________________
//  DEPRECATED
bool    Bgm::IsCyclic (void)
{
    _Matrix     nodes_left (num_nodes, 1, false, true);

    for (long i = 0; i < num_nodes; i++) {  // populate vector with 1's
        nodes_left.Store (i,0,1.);
    }

    long        one_count = num_nodes,
                last_count,
                num_children;

    do {    // until no more leaves are removed
        last_count = one_count;
        one_count = 0;

        // find and remove all nodes with no children (i.e. leaf nodes)
        for (long parent = 0; parent < num_nodes; parent++) {
            if (nodes_left(parent,0) == 1.) {
                one_count++;

                num_children = 0;
                // count all existing nodes that have this parent
                for (long child = 0; child < num_nodes; child++)
                    if (nodes_left(child,0) == 1. && child != parent)
                        if (dag(parent,child) == 1.) {
                            num_children++;
                        }

                if (num_children == 0) {    // count number of nodes with this node as parent
                    nodes_left.Store (parent, 0, 0.);
                    one_count--;
                }
            }
        }
    } while (one_count < last_count);


    if (one_count > 0) {
        return TRUE;    // dag is cyclic
    } else {
        return FALSE;   // removed all nodes as leaves, dag is acyclic
    }
}




//___________________________________________________________________________________________
void    Bgm::RandomizeGraph (_Matrix * graphMx, _SimpleList * order, long num_steps, bool fixed_order)
{
    long    step = 0, fail = 0;


    // convert order into matrix format of edge permissions for more rapid look-up later  /* DEBUGGING ONLY */
#ifdef __DEBUG_RG__
    long    mode = 0;
    _Matrix orderMx (num_nodes, num_nodes, false, true);

    for (long p_rank = 0; p_rank < num_nodes; p_rank++) {
        for (long c_rank = 0; c_rank < num_nodes; c_rank++) {
            // use actual node id's to index into matrix - 1 indicates a permitted edge
            orderMx.Store (order->lData[p_rank], order->lData[c_rank], p_rank > c_rank ? 1 : 0);
        }
    }
#endif


    // before we do anything, make sure that graph complies with order /* DEBUGGING */
    // calculate number of parents for each child for rapid access later
    _SimpleList     num_parents;

    num_parents.Populate (num_nodes, 0, 0);

    for (long child = 0; child < num_nodes; child++) {
        for (long parent = 0; parent < num_nodes; parent++) {
#ifdef __DEBUG_RG__
            if ((*graphMx)(parent,child) > 0 && banned_edges(parent,child) > 0) {
                WarnError (_String("BEFORE RandomizeGraph() there is a banned edge: ") & parent & "->" & child & " mode " & mode);
                break;
            }
            if ((*graphMx)(parent,child) == 0 && enforced_edges(parent,child) > 0) {
                WarnError (_String("BEFORE RandomizeGraph() there is a missing enforced edge: ") & parent & "->" & child);
                break;
            }

            if ( (*graphMx)(parent, child)==1 && orderMx(parent, child)==0 ) {
                WarnError (_String ("Order-network discrepancy at start of RandomizeGraph() at edge ") & parent & "->" & child);
                break;
            }
#endif

            if ( (*graphMx)(parent, child) > 0) {
                num_parents.lData[child]++;
            }
        }
        // end for

        if (num_parents.lData[child] > max_parents.lData[child]) {
            WarnError (_String ("Number of parents exceeds maximum BEFORE randomization of graph at node ") & child & " (" & num_parents.lData[child] & " > " & max_parents.lData[child] & ")\n" );
            break;
        }
    }




    // randomize graph
    long    p_rank, c_rank, parent, child;

    do {
        // a fail-safe to avoid infinite loops
        if (fail > MAX_FAIL_RANDOMIZE) {
            WarnError(_String("Bgm::RandomizeGraph() failed to modify the graph in GraphMCMC() after MAX_FAIL_RANDOMIZE (") & (long)MAX_FAIL_RANDOMIZE & ") attempts.\n");
            break;
        }

        // pick a random edge permitted under ordering
        p_rank  = (genrand_int32() % (num_nodes-1)) + 1;    // shift right to not include lowest node in order
        c_rank  = genrand_int32() % p_rank;

        child = order->lData[c_rank];
        parent = order->lData[p_rank];


        if (fixed_order || genrand_real2() > RANDOMIZE_PROB_SWAP) { // attempt an add or delete
            if ( (*graphMx)(parent,child) == 0 && !banned_edges(parent,child)) {    // add an edge
#ifdef __DEBUG_RG__
                mode = 0;
#endif
                if (num_parents.lData[child] == max_parents.lData[child]) { // child cannot accept any additional edges
                    // move an edge from an existing parent to target parent
                    _SimpleList     removeable_edges;

                    // build a list of current parents
                    for (long par = 0; par < num_nodes; par++)
                        if ((*graphMx)(par,child) && !enforced_edges(par,child)) {
                            removeable_edges << par;
                        }

                    if (removeable_edges.lLength > 0) {
                        // shuffle the list and remove the first parent
                        removeable_edges.Permute(1);
                        graphMx->Store (removeable_edges.lData[0], child, 0.);
                        graphMx->Store (parent, child, 1.);
                        step++;
                    } else {
                        // none of the edges can be removed
                        fail++;
                    }
                } else {
                    // child can accept another edge
                    graphMx->Store (parent, child, 1.);
                    num_parents.lData[child]++;
                    step++;
                }
            } else if ( (*graphMx)(parent,child) == 1 && !enforced_edges(parent,child)) {   // delete an edge
#ifdef __DEBUG_RG__
                mode = 1;
#endif
                graphMx->Store (parent, child, 0.);
                num_parents.lData[child]--;
                step++;
            } else {
                fail++;
            }
        } else {    // swap nodes in ordering and flip edge if present
            long    ok_to_go = 1;
#ifdef __DEBUG_RG__
            mode = 2;
#endif
            // edge cannot be flipped
            if ( (*graphMx)(parent,child) == 1  &&  ( banned_edges(child,parent) || enforced_edges(parent,child) )  ) {
                ok_to_go = 0;
            }


            // check all other nodes affected by the swap
            if (ok_to_go) {
                for (long bystander, i=c_rank+1; i < p_rank; i++) {
                    bystander = order->lData[i];    // retrieve node id

                    if (
                        ( (*graphMx)(parent,bystander)==1 && enforced_edges (parent, bystander) )  ||
                        ( (*graphMx)(bystander,child)==1 && enforced_edges (bystander, child) )  ||
                        banned_edges (bystander,parent)  ||  banned_edges (child,bystander)
                    ) {
                        // by flipping the parent->child edge, we would screw up ordering for
                        // an enforced edge:  C < B < P   becomes  P < B < C  where B->C or P->B is enforced

                        // also, a banned edge implies a node ordering constraint
                        ok_to_go = 0;
                        break;
                    }
                }
            }


            // if everything checks out OK
            if ( ok_to_go ) {
                // flip the target edge
                if ( (*graphMx)(parent,child) == 1 && !banned_edges(child,parent) ) {
                    graphMx->Store (parent, child, 0);
                    graphMx->Store (child, parent, 1);
                    num_parents.lData[child]--;
                    num_parents.lData[parent]++;

                    if (num_parents.lData[parent] > max_parents.lData[parent]) {
                        // parent cannot accept any more edges, delete one of the edges at random (including the edge to flip)
                        _SimpleList     removeable_edges;

                        for (long par = 0; par < num_nodes; par++)
                            if (  (*graphMx)(par, parent)  &&  !enforced_edges(par, parent)  ) {
                                removeable_edges << par;
                            }

                        removeable_edges.Permute(1);
                        graphMx->Store (removeable_edges.lData[0], parent, 0.);
                        num_parents.lData[parent]--;
                    }

                    step++;
                }
                // if number of parents for parent node will exceed maximum, then the edge is deleted instead of flipped

                // swap nodes in order
                order->lData[p_rank] = child;
                order->lData[c_rank] = parent;  // remember to update order matrix also!


                // flip the other edges affected by node swap
                //  child <-  N   ...   N <- parent                   _________________,
                //                                    becomes        |                 v
                //                                                 parent    N         N    child
                //                                                           ^                |
                //                                                           `----------------+
                for (long bystander, i = c_rank+1; i < p_rank; i++) {
                    bystander = order->lData[i];

                    if ( (*graphMx)(bystander, child) == 1 ) {
                        graphMx->Store (bystander, child, 0);
                        num_parents.lData[child]--;

                        graphMx->Store (child, bystander, 1);
                        num_parents.lData[bystander]++;

                        if (num_parents.lData[bystander] > max_parents.lData[bystander]) {
                            _SimpleList     removeable_edges;

                            for (long par = 0; par < num_nodes; par++)
                                if (  (*graphMx)(par, bystander)  &&  !enforced_edges(par, bystander)  ) {
                                    removeable_edges << par;
                                }

                            removeable_edges.Permute(1);
                            graphMx->Store (removeable_edges.lData[0], bystander, 0.);
                            num_parents.lData[bystander]--;
                        }
                    }

                    if ( (*graphMx)(parent,bystander) == 1) {
                        graphMx->Store (parent, bystander, 0);
                        num_parents.lData[bystander]--;

                        graphMx->Store (bystander, parent, 1);
                        num_parents.lData[parent]++;

                        if (num_parents.lData[parent] > max_parents.lData[parent]) {
                            _SimpleList     removeable_edges;

                            for (long par = 0; par < num_nodes; par++)
                                if (  (*graphMx)(par, parent)  &&  !enforced_edges(par, parent)  ) {
                                    removeable_edges << par;
                                }

                            removeable_edges.Permute(1);
                            graphMx->Store (removeable_edges.lData[0], parent, 0.);
                            num_parents.lData[parent]--;
                        }
                    }
                }
#ifdef __DEBUG_RG__
                // refresh order matrix
                for (long p_rank = 0; p_rank < num_nodes; p_rank++) {
                    for (long c_rank = 0; c_rank < p_rank; c_rank++) {
                        orderMx.Store (order->lData[p_rank], order->lData[c_rank], 1);
                    }
                }
#endif
                step++;
            } else {
                fail++;
            }
        }
    } while (step < num_steps);


    // a final check to make sure we haven't screwed up!
#ifdef __DEBUG_RG__
    num_parents.Populate (num_nodes, 0, 0);

    for (long child = 0; child < num_nodes; child++) {
        for (long parent = 0; parent < num_nodes; parent++) {
            if ((*graphMx)(parent,child) > 0 && banned_edges(parent,child) > 0) {
                WarnError (_String("Aw crap!  RandomizeGraph() introduced a banned edge: ") & parent & "->" & child & " mode " & mode);
                break;
            }
            if ((*graphMx)(parent,child) == 0 && enforced_edges(parent,child) > 0) {
                WarnError (_String("Aw crap!  RandomizeGraph() deleted an enforced edge: ") & parent & "->" & child);
                break;
            }
            if ( (*graphMx)(parent, child)==1 && orderMx(parent, child)==0 ) {
                WarnError (_String ("Order-network discrepancy at end of RandomizeGraph(), mode ") & mode);
                break;
            }

            if ( (*graphMx)(parent, child) > 0) {
                num_parents.lData[child]++;
            }
        }

        if (num_parents.lData[child] > max_parents.lData[child]) {
            WarnError (_String ("Exceeded acceptable number of parents for node ") & child &  ", mode " & mode);
            break;
        }
    }
#endif
}



//___________________________________________________________________________________________
void    Bgm::RandomizeDag (long num_steps)
{

    for (long step = 0; step < num_steps; step++) {
        // store current graph in case altered graph becomes cyclic
        _Matrix reset_dag (num_nodes, num_nodes, false, true);
        for (long h = 0; h < num_nodes; h++)
            for (long v = 0; v < num_nodes; v++) {
                reset_dag.Store (h,v,dag(h,v));
            }

        do {
            for (long row = 0; row < num_nodes; row++) {
                for (long col = 0; col < num_nodes; col++) {
                    dag.Store (row,col,reset_dag(row,col));
                }
            }


            // select a random arc, discarding cycles and banned edges
            long    row     = 0,
                    col      = row;

            while (col == row && banned_edges(row,col) == 1) {
                row = genrand_int32() % num_nodes;
                col = genrand_int32() % num_nodes;
            }

            // printf ("RandomizeDag() picked %d,%d\n", row, col);

            if (dag(row,col) == 0.) {
                if (dag(col,row) == 0.) {   // no arc
                    // add an arc -- sum over n-th column gives number of parents for n-th child
                    if (MarginalSum(FALSE,col) < max_parents.lData[col]) {
                        dag.Store (row,col,1.);
                    }
                } else {
                    // reverse or remove arc
                    dag.Store (col,row,0.);
                    if (genrand_int32() % 2 == 0) {
                        dag.Store (row,col,1.);
                    }
                }
            } else {
                if (dag(col,row) == 0.) {
                    // reverse or remove
                    dag.Store (row,col,0.);
                    if (genrand_int32() % 2 == 0) {
                        dag.Store (col,row,1.);
                    }
                } else {
                    // double arc, set to one of three legal arcs
                    long    option = genrand_int32() % 3;
                    if (option == 0) {
                        dag.Store (row,col,0.);
                    } else if (option == 1) {
                        dag.Store (col,row,0.);
                    } else {
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
//  Wrappers to retain original functionality.
_Parameter  Bgm::ComputeDiscreteScore (long node_id)
{
    _SimpleList     parents;

    for (long par = 0; par < num_nodes; par++) {
        if (dag(par, node_id) == 1 && is_discrete.lData[par]) {
            parents << par;
        }
    }

    return ComputeDiscreteScore (node_id, parents);
}

_Parameter  Bgm::ComputeDiscreteScore (long node_id, _Matrix * g)
{
    _SimpleList     parents;

    for (long par = 0; par < num_nodes; par++) {
        if ((*g)(par, node_id) == 1 && is_discrete.lData[par]) {
            parents << par;
        }
    }

    return ComputeDiscreteScore (node_id, parents);
}


//___________________________________________________________________________________________
//#define BGM_DEBUG_CDS
_Parameter  Bgm::ComputeDiscreteScore (long node_id, _SimpleList & parents)
{
    //char          buf [255];

    // use cached node scores if possible
    if (scores_cached) {
        _List *     scores  = (_List *) node_scores.lData[node_id];

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



#ifdef __MISSING_DATA__
    //  Are any of the edges banned?
    for (long par = 0; par < parents.lLength; par++) {
        if (banned_edges(parents.lData[par], node_id) > 0) {
            // score should never be used
            ReportWarning(_String("Skipping node score for family containing banned edge ") & parents.lData[par] & "->" & node_id & "\n");
            return (-A_LARGE_NUMBER);
        }
    }


    //  Is node with missing data in Markov blanket of focal node?
    if (has_missing.lData[node_id]) {
        //return (ImputeDiscreteScore (node_id, parents));
        return (GibbsApproximateDiscreteScore (node_id, parents));
    } else {
        for (long par = 0; par < parents.lLength; par++) {
            if (has_missing.lData[parents.lData[par]]) {
                // return (ImputeDiscreteScore (node_id, parents));
                return (GibbsApproximateDiscreteScore (node_id, parents));
            }
        }
    }
#endif


    _SimpleList     multipliers ((long)1);

    // else need to compute de novo
    _Matrix     n_ijk,
                n_ij;       // [i] indexes child nodes,
    // [j] indexes combinations of values for parents of i-th node,
    // [k] indexes values of i-th node

    long        num_parent_combos   = 1,                    // i.e. 'q'
                r_i                   = num_levels.lData[node_id];

    _Parameter  n_prior_ijk = 0,
                n_prior_ij  = 0,
                log_score  = 0;



    // how many combinations of parental states are there?
    for (long par = 0; par < parents.lLength; par++) {
        num_parent_combos *= num_levels.lData[parents.lData[par]];
        multipliers << num_parent_combos;
    }


#ifdef BGM_DEBUG_CDS
    snprintf (buf, sizeof(buf), "Multipliers: ");
    BufferToConsole (buf);
    for (long i = 0; i < multipliers.lLength; i++) {
        snprintf (buf, sizeof(buf), "%d ", multipliers.lData[i]);
        BufferToConsole (buf);
    }
    NLToConsole();
#endif


    /* count observations by parent combination, using direct indexing */
    CreateMatrix (&n_ijk, num_parent_combos, r_i, false, true, false);
    CreateMatrix (&n_ij, num_parent_combos, 1, false, true, false);


#ifdef __MISSING_DATA__
    /*  METHODS FOR COMPLETE DATA  */
    for (long obs = 0; obs < obsData->GetHDim(); obs++) {
        long    index           = 0,
                //multiplier    = 1,
                child_state     = (*obsData)(obs, node_id);

        for (long par = 0; par < parents.lLength; par++) {
            long    this_parent         = parents.lData[par],
                    this_parent_state = (*obsData)(obs, this_parent);

            index += this_parent_state * multipliers.lData[par];
        }

        n_ijk.Store ((long) index, child_state, n_ijk(index, child_state) + 1 );
        n_ij.Store ((long) index, 0, n_ij(index, 0) + 1 );
    }


#ifdef BGM_DEBUG_CDS
    snprintf (buf, sizeof(buf), "Node %d, parent(s) ", node_id);
    BufferToConsole (buf);

    for (long k = 0; k < parents.lLength; k++) {
        snprintf (buf, sizeof(buf), "%d ", parents.lData[k]);
        BufferToConsole (buf);
    }
    NLToConsole();

    for (long j = 0; j < num_parent_combos; j++) {
        snprintf (buf, sizeof(buf), "N(%d,%d) = %f\n", node_id, j, n_ij(j,0));
        BufferToConsole (buf);

        for (long k = 0; k < r_i; k++) {
            snprintf (buf, sizeof(buf), "N(%d,%d,%d) = %f\n", node_id, j, k, n_ijk(j,k));
            BufferToConsole (buf);
        }
    }
#endif


    if (prior_sample_size (node_id, 0) == 0) {  // K2
        for (long j = 0; j < num_parent_combos; j++) {
            log_score += LnGamma(num_levels.lData[node_id]);    // (r-1)!
            log_score -= LnGamma(n_ij(j, 0) + num_levels.lData[node_id]);   // (N+r-1)!

            for (long k = 0; k < r_i; k++) {
                log_score += LnGamma (n_ijk(j,k) + 1);    // (N_ijk)!
            }

#ifdef BGM_DEBUG_CDS
            snprintf (buf, sizeof(buf), "\tlog(r-1)! = log %d! = %lf\n", num_levels.lData[node_id] - 1, LnGamma(num_levels.lData[node_id]));
            BufferToConsole (buf);

            snprintf (buf, sizeof(buf), "\tlog(N(%d,%d)+r-1)! = log %d! = %lf\n", node_id, j, ((long)n_ij(j, 0)) + r_i - 1,
                     LnGamma(n_ij(j, 0) + num_levels.lData[node_id]));
            BufferToConsole (buf);

            for (long k = 0; k < r_i; k++) {
                snprintf (buf, sizeof(buf), "\tlog (N(%d,%d,%d)!) = log %d! = %lf\n", node_id, j, k, ((long)n_ijk(j,k)), LnGamma(n_ijk(j,k) + 1));
                BufferToConsole (buf);
            }

            snprintf (buf, sizeof(buf), "j = %d\tcumulative log score = %lf\n", j, log_score);
            BufferToConsole (buf);
#endif
        }
    } else {    // BDe
        n_prior_ij = prior_sample_size (node_id, 0) / num_parent_combos;
        n_prior_ijk = n_prior_ij / num_levels.lData[node_id];

        for (long j = 0; j < num_parent_combos; j++) {
            log_score += LnGamma(n_prior_ij) - LnGamma(n_prior_ij + n_ij(j,0));

            for (long k = 0; k < num_levels.lData[node_id]; k++) {
                log_score += LnGamma(n_prior_ijk + n_ijk(j,k)) - LnGamma(n_prior_ijk);
            }
        }
    }

#else
    for (long obs = 0; obs < obsData->GetHDim(); obs++) {
        long    index       = 0,
                //multiplier   = 1,
                child_state = (*obsData)(obs, node_id);

        for (long par = 0; par < parents.lLength; par++) {
            long    this_parent         = parents.lData[par],
                    this_parent_state = (*obsData)(obs, this_parent);


            /*      // this system doesn't work with unequal levels!  :-P  afyp March 31, 2008
            index += this_parent_state * multiplier;
            multiplier *= num_levels.lData[this_parent];
             */
            index += this_parent_state * multipliers.lData[par];
        }

        n_ijk.Store ((long) index, child_state, n_ijk(index, child_state) + 1);
        n_ij.Store ((long) index, 0, n_ij(index, 0) + 1);

    }



    /* compute scoring metric */
    if (prior_sample_size (node_id, 0) == 0) {
        /* assume no prior information, use K2 metric */
        for (long j = 0; j < num_parent_combos; j++) {
            log_score += LnGamma(r_i);  // (r-1)!
            log_score -= LnGamma(n_ij(j, 0) + r_i); // (N+r-1)!

            for (long k = 0; k < num_levels.lData[node_id]; k++) {
                log_score += LnGamma(n_ijk(j,k) + 1);    // (N_ijk)!
            }



            snprintf (buf, sizeof(buf), "Node %d, parent(s) ", node_id);
            BufferToConsole (buf);

            for (long k = 0; k < parents.lLength; k++) {
                snprintf (buf, sizeof(buf), "%d ", parents.lData[k]);
                BufferToConsole (buf);
            }
            NLToConsole();

            snprintf (buf, sizeof(buf), "j = %d\tcumulative log score = %lf\n", j, log_score);
            BufferToConsole (buf);

#ifdef BGM_DEBUG_CDS
            snprintf (buf, sizeof(buf), "\tlog(r-1)! = log %d! = %lf\n", num_levels.lData[node_id] - 1, LnGamma(num_levels.lData[node_id]));
            BufferToConsole (buf);

            snprintf (buf, sizeof(buf), "\tlog(N+r-1)! = log %d! = %lf\n", n_ij(j, 0) + num_levels.lData[node_id] - 1,
                     LnGamma(n_ij(j, 0) + num_levels.lData[node_id]));
            BufferToConsole (buf);

            for (long k = 0; k < num_levels.lData[node_id]; k++) {
                snprintf (buf, sizeof(buf), "\tlog (N_ijk)! = log %d! = %lf\n", ((long)n_ijk(j,k)), LnGamma(n_ijk(j,k) + 1));
                BufferToConsole (buf);
            }


#endif
        }
    } else {
        /* calculate Bayesian Dirichlet metric (BDeu) for this node */
        /* see p212 in Heckerman, Geiger, and Chickering (1995) Machine Learning 20, 197-243 */
        n_prior_ij = prior_sample_size (node_id, 0) / num_parent_combos;
        n_prior_ijk = n_prior_ij / num_levels.lData[node_id];

        for (long j = 0; j < num_parent_combos; j++) {
            log_score += LnGamma(n_prior_ij) - LnGamma(n_prior_ij + n_ij(j,0));

            for (long k = 0; k < num_levels.lData[node_id]; k++) {
                log_score += LnGamma(n_prior_ijk + n_ijk(j,k)) - LnGamma(n_prior_ijk);
            }
        }
    }
#endif

    /*
    snprintf (buf, sizeof(buf), "Node %d, parents (%d): ", node_id, parents.lLength);
    BufferToConsole (buf);

    for (long par = 0; par < parents.lLength; par++)
    {
        snprintf (buf, sizeof(buf), " %d", parents.lData[par]);
        BufferToConsole (buf);
    }
    snprintf (buf, sizeof(buf), " Log score = %f\n", log_score);
    BufferToConsole (buf);
    */

    return (log_score);
}




//___________________________________________________________________________________________
//#define __DEBUG_CNS__
void Bgm::CacheNodeScores (void)
{
    if (scores_cached) {
        return;
    }

    /*
    char buf [255];
    snprintf (buf, sizeof(buf), "\nCaching node scores...\n");
    BufferToConsole (buf);
    */


#if defined __HYPHYMPI_ND__
    _Parameter  use_mpi_caching;
    checkParameter (useMPIcaching, use_mpi_caching, 0);

    if (use_mpi_caching) {

        // MPI_Init() is called in main()
        int     size,
                rank;

        long    mpi_node;

        _Matrix     single_parent_scores (num_nodes, 1, false, true);

        _SimpleList parents,
                    all_but_one (num_nodes-1, 0, 1),
                    aux_list,
                    nk_tuple;

        _Parameter  score;

        char        mpi_message [256];

        MPI_Status  status; // contains source, tag, and error code

        MPI_Comm_size (MPI_COMM_WORLD, &size);  // number of processes
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);


        if (rank == 0) {
            _String     bgmSwitch ("_BGM_SWITCH_"),
                        bgmStr;

            _List   *   this_list;

            SerializeBgm (bgmStr);
#ifdef __DEBUG_CNS__
            ReportWarning (bgmStr);
#endif
            _Matrix     * mpi_node_status = new _Matrix ((long)size, 2, false, true);
            long        senderID;



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


            long node_id;

            // farm out jobs to idle nodes until none are left
            for ( node_id = 0; node_id < num_nodes; node_id++) {
                long        maxp            = max_parents.lData[node_id],
                            ntuple_receipt,
                            this_node;

                bool        remaining;

                _String     mxString,
                            mxName;


                this_list   = (_List *) node_scores.lData[node_id];

                // [_SimpleList parents] should always be empty here
                _Parameter  score       = is_discrete.lData[node_id] ? ComputeDiscreteScore (node_id, parents) : ComputeContinuousScore (node_id);
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


                if (mpi_node < size) {  // at least one node is busy
                    MPIReceiveScores (mpi_node_status, false, 0);
                } else {
                    break;
                }
            }


            // shut down compute nodes
            for (long shutdown = -1, mpi_node = 1; mpi_node < size; mpi_node++) {
                ReportMPIError(MPI_Send(&shutdown, 1, MPI_LONG, mpi_node, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD), true);
                ReportWarning (_String ("Node 0 sending shutdown signal to node ") & mpi_node);
            }

        } else {    // compute node
            long        node_id,
                        maxp;

            _List       list_of_matrices;

            while (1) {
                _String     mxString,
                            mxName;

                list_of_matrices.Clear();


                // wait for master node to issue node ID to compute
                ReportMPIError (MPI_Recv (&node_id, 1, MPI_LONG, 0, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD, &status), false);
                ReportWarning (_String("Node ") & (long)rank & " received child " & (long)node_id & " from node " & (long)status.MPI_SOURCE & "\n");

                if (node_id < 0) {
                    ReportWarning (_String("Node") & (long)rank & " recognizes shutdown signal.\n");
                    break;  // received shutdown message (-1)
                }


                maxp = max_parents.lData[node_id];

                parents.Clear();
                parents.Populate (1,0,0);


                // compute single parent scores
                for (long par = 0; par < num_nodes; par++) {
                    if (par == node_id) {   // child cannot be its own parent, except in Kansas
                        single_parent_scores.Store (par, 0, 0.);
                    } else {
                        parents.lData[0] = par;
                        single_parent_scores.Store (par, 0, is_discrete.lData[node_id] ? ComputeDiscreteScore (node_id, parents) :
                                                    ComputeContinuousScore (node_id));
                    }
                }


                // compute multiple parents cores
                for (long np = 2; np <= maxp; np++) {
                    parents.Clear();
                    parents.Populate (np, 0, 0);

                    if (all_but_one.NChooseKInit (aux_list, nk_tuple, np, false)) {
                        bool        remaining;
                        long        tuple_index     = 0,
                                    num_ktuples      = exp(LnGamma(num_nodes) - LnGamma(num_nodes - np) - LnGamma(np+1));


                        _Matrix     tuple_scores (num_ktuples, 1, false, true);

                        for (long tuple_index = 0; tuple_index < num_ktuples; tuple_index++) {
                            remaining = all_but_one.NChooseK (aux_list, nk_tuple);

                            if (!remaining && tuple_index < num_ktuples-1) {
                                ReportWarning (_String ("ERROR: Ran out of (n,k)tuples in CacheNodeScores()."));
                            }

                            for (long par_idx = 0; par_idx < np; par_idx++) {
                                long par = nk_tuple.lData[par_idx];
                                if (par >= node_id) {
                                    par++;
                                }
                                parents.lData[par_idx] = par;
                            }

                            score = is_discrete.lData[node_id] ? ComputeDiscreteScore (node_id, parents) :
                                    ComputeContinuousScore (node_id);

                            tuple_scores.Store (tuple_index, 0, (double)score);
                        }

                        if (remaining) {
                            ReportWarning (_String ("ERROR: Did not compute all nk-tuples in CacheNodeScores()"));
                        }

                        list_of_matrices && (&tuple_scores);        // append duplicate
                    } else {
                        ReportWarning (_String ("Failed to initialize _NTupleStorage object in Bgm::CacheNodeScores()"));
                    }

                }

                // send results to master node

                ReportMPIError (MPI_Send (single_parent_scores.theData, num_nodes, MPI_DOUBLE, 0, HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD), true);


                for (long np = 2; np <= maxp; np++) {
                    _Matrix * storedMx      = (_Matrix *) list_of_matrices.lData[np-2];
                    ReportMPIError (MPI_Send (storedMx->theData, storedMx->GetHDim(), MPI_DOUBLE, 0, HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD), true);
                }
            }
        }
    } else {    // perform single-threaded
#else


#if !defined __UNIX__ || defined __HEADLESS__
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
        _List   *   this_list   = (_List *) node_scores.lData[node_id]; // retrieve pointer to list of scores for this child node

        this_list->Clear();     // reset the list

        // prepare some containers
        _SimpleList parents;
        _Parameter  score = is_discrete.lData[node_id] ? ComputeDiscreteScore (node_id, parents) : ComputeContinuousScore (node_id);
        _Constant   orphan_score (score);

        (*this_list) && (&orphan_score);

#if !defined __UNIX__ || defined __HEADLESS__
        temp = .0;
#endif

        if (maxp > 0) {
            _Matrix     single_parent_scores (num_nodes, 1, false, true);
            for (long par = 0; par < num_nodes; par++) {
                if (par == node_id) {   // child cannot be its own parent, except in Kansas
                    continue;
                }

                parents << par;
                single_parent_scores.Store (par, 0, is_discrete.lData[node_id] ? ComputeDiscreteScore (node_id, parents) :
                                            ComputeContinuousScore (node_id));
                parents.Clear();
            }
            (*this_list) && (&single_parent_scores);
        }

        if (maxp > 1) {
            _SimpleList     all_but_one (num_nodes-1, 0, 1),    // 0, 1, 2, ... , n-1
                            aux_list,
                            nk_tuple;

            for (long np = 2; np <= maxp; np++) {
                _NTupleStorage  family_scores (num_nodes-1, np);

                if (all_but_one.NChooseKInit (aux_list, nk_tuple, np, false)) {
                    bool    remaining;
                    long    res;
                    do {
                        remaining = all_but_one.NChooseK (aux_list, nk_tuple);
                        for (long par_idx = 0; par_idx < np; par_idx++) {
                            long par = nk_tuple.lData[par_idx];
                            if (par >= node_id) {
                                par++;
                            }
                            parents << par;
                        }
                        score = is_discrete.lData[node_id]  ?   ComputeDiscreteScore (node_id, parents) :
                                ComputeContinuousScore (node_id);
                        res = family_scores.Store (score, nk_tuple);
                        parents.Clear();
                    } while (remaining);
                } else {
                    _String oops ("Failed to initialize _NTupleStorage object in Bgm::CacheNodeScores().\n");
                    WarnError(oops);
                }
                (*this_list) && (&family_scores);   // append duplicate to storage
            }
        }

#if !defined __UNIX__ || defined __HEADLESS__
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

    } // end for loop over nodes
#endif


#if defined __HYPHYMPI_ND__
    }
#endif

#if !defined __UNIX__ || defined __HEADLESS__
    SetStatusLine     (_HYBgm_STATUS_LINE_CACHE_DONE);
#endif

    scores_cached = TRUE;


}




//___________________________________________________________________________________________
#if defined __HYPHYMPI__
void    Bgm::MPIReceiveScores (_Matrix * mpi_node_status, bool sendNextJob, long node_id)
{
    _Matrix     single_parent_scores (num_nodes, 1, false, true);
    MPI_Status  status;

    ReportMPIError (MPI_Recv (single_parent_scores.theData, num_nodes, MPI_DOUBLE, MPI_ANY_SOURCE, HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD, &status), false);


    long        senderID    = (long) status.MPI_SOURCE,
                this_node    = (long) (*mpi_node_status) (senderID, 1),
                maxp       = max_parents.lData[this_node];

    _List   *   this_list   = (_List *) node_scores.lData[this_node];


    mpi_node_status->Store (senderID, 0, 0);    // set node status to idle


    _String     mxString,
                mxName;

    ReportWarning (_String("Received scores for child ") & this_node & " from node " & senderID);
#ifdef __DEBUG_MPIRS__
    mxName = _String ("single_parent_scores");
    single_parent_scores.Serialize (mxString, mxName);
    ReportWarning (mxString);
#endif
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
                        num_ktuples      = exp(LnGamma(num_nodes) - LnGamma(num_nodes - np) - LnGamma(np+1)),
                        ntuple_receipt;

            _Matrix     scores_to_store (num_ktuples, 1, false, true);

            // receive nk-tuple indexed node scores from same compute node
            ReportMPIError (MPI_Recv (scores_to_store.theData, num_ktuples, MPI_DOUBLE, senderID, HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD, &status), false);
#ifdef __DEBUG_MPIRS__
            mxName = _String("tuple_scores") & np;
            mxString.Initialize();
            scores_to_store.Serialize (mxString, mxName);
            ReportWarning (mxString);
#endif

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
#endif



//___________________________________________________________________________________________
_Matrix *   Bgm::ExportGraph (void)
{
    _Matrix * export_graph = new _Matrix ((_Matrix &)dag);  // duplicator

    return (_Matrix *) (export_graph->makeDynamic());
}




//___________________________________________________________________________________________
//  Allocate memory for storing node scores and edge posteriors

void    Bgm::InitComputeLists (_List * compute_list)
{
    // adaptive storage for floating point numbers
    _GrowingVector *    newstore;
    checkPointer (newstore = new _GrowingVector);

    for (long i = 0; i < num_nodes * num_nodes; i++) {
        (*compute_list) && newstore;
    }

    DeleteObject (newstore);
}



//___________________________________________________________________________________________
//  Free up memory allocated to cached node scores and edge posteriors

void        Bgm::DumpComputeLists (_List * compute_list)
{
    for (long i = 0; i < compute_list->lLength; i++) {
        ((_GrowingVector *) compute_list->lData[i]) -> Clear();
    }

    compute_list->Clear();
}



//___________________________________________________________________________________________
_Parameter  Bgm::Compute (_SimpleList * node_order, _List * results)
{
    /*  Calculate equation (8) from Friedman and Koller (2003), i.e. joint probability
        of data by summing all families at i-th node that are consistent with node order,
        then taking product across all nodes in the network.  Traverse using AVL indexing.

        Node order is a vector of length [num_nodes] containing rank of each node, i.e.
        a position in an ordered sequence.  In addition, compute marginal posteriors of
        each potential edge in the graph, by Proposition 3.2.                               */



    _Parameter          log_likel   = 0.;
    _GrowingVector      *gv1, *gv2;


    // reset _GrowingVector objects stored in _List object
    for (long i = 0; i < num_nodes * num_nodes; i++) {
        gv1 = (_GrowingVector *) results->lData[i];
        gv1 -> ZeroUsed();
    }


    for (long nodeIndex = 0; nodeIndex < node_order->lLength; nodeIndex++) {
        long                child_node      = node_order->lData[nodeIndex],
                            maxp            = max_parents.lData[child_node];

        _List           *   score_lists     = (_List *) node_scores.lData[child_node];
        _Constant       *   orphan_score    = (_Constant *) (score_lists->lData[0]);


        gv1 = (_GrowingVector *) results->lData[child_node * num_nodes + child_node];   // store denominator in diagonal
        gv1->ZeroUsed();
        gv1 -> Store (orphan_score->Value());   // handle case of no parents



        if (maxp > 0) {
            // all nodes to the right are potential parents, except banned parents!
            _SimpleList     precedes;
            for (long parIndex = nodeIndex + 1; parIndex < node_order->lLength; parIndex++) {
                long    par = node_order->lData[parIndex];

                if (banned_edges(par, child_node) == 0) {
                    precedes << par;
                }
            }


            // handle trivial case of one parent
            _Matrix *   single_parent_scores    = (_Matrix *) (score_lists->lData[1]);

            for (long i = 0; i < precedes.lLength; i++) {
                long    par = precedes.lData[i];

                gv1 -> Store ((*single_parent_scores) (par, 0));
                gv2 = (_GrowingVector *) results->lData[child_node * num_nodes + par];
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

                            gv1 -> Store (tuple_score);

                            for (long i = 0; i < nparents; i++) {
                                gv2 = (_GrowingVector *) results->lData[child_node * num_nodes + precedes.lData[subset.lData[i]]];
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



//___________________________________________________________________________________________
_Parameter Bgm::Compute (void)
{
    //CacheNodeScores();

    // return posterior probability of a given network defined by 'dag' matrix
    _Parameter  log_score = 0.;

    for (long node_id = 0; node_id < num_nodes; node_id++) {
        log_score += is_discrete.lData[node_id] ? ComputeDiscreteScore (node_id) : ComputeContinuousScore (node_id);
    }

    return log_score;
}



//___________________________________________________________________________________________
_Parameter Bgm::Compute (_Matrix * g)
{
    //CacheNodeScores();

    // return posterior probability of a given network defined by 'dag' matrix
    _Parameter  log_score = 0.;

    for (long node_id = 0; node_id < num_nodes; node_id++) {
        log_score += is_discrete.lData[node_id] ? ComputeDiscreteScore (node_id, g) : ComputeContinuousScore (node_id, g);
    }

    return log_score;
}




//___________________________________________________________________________________________
_Matrix *   Bgm::Optimize (void)
{
#ifdef DEBUG_OPTIMIZE
    char        buf [255];
#endif

    if (!scores_cached) {
        CacheNodeScores();
    }


    _Parameter  num_restarts,           // HBL settings
                num_randomize;

    checkParameter (maxNumRestart, num_restarts, 1.);
    checkParameter (numRandomize, num_randomize, num_nodes);
    //checkParameter (useNodeOrder, use_node_order, 0);



    // char     bug [255];
    _Parameter      optMethod;  /* 0 = K2 fixed order; 1 = K2 shuffle order with restarts */
    /* 2 = MCMC fixed order; 3 = MCMC over networks and orders */
    checkParameter (bgmOptimizationMethod, optMethod, 0.);

    if (optMethod > 1) {
        return GraphMCMC( (optMethod < 3) ? TRUE : FALSE);      // use MCMC to evaluate posterior probability of network space
    }



    _Parameter  this_score, best_score, next_score;

    _Matrix     orderMx (num_nodes, num_nodes, false, true),
                best_dag (num_nodes, num_nodes, false, true);

    bool        reshuffle_order = FALSE;


    // if best node order hasn't been estimated,
    if (optMethod == 1 && best_node_order.lLength == 0) {
        reshuffle_order = TRUE;
        best_node_order.Populate (num_nodes, 0, 1);
        best_node_order.Permute (1);
    }


    //  Convert node order to binary matrix where edge A->B is permitted if
    //  orderMx[B][A] = 1, i.e. B is to the right of A in node order.
    for (long i = 0; i < num_nodes; i++) {
        long    child = best_node_order.lData[i];

        for (long j = 0; j < num_nodes; j++) {
            long    parent = best_node_order.lData[j];

            orderMx.Store (parent, child, (j > i) ? 1 : 0);
        }
    }


    // greedy hill-climbing algorithm (K2)
    ResetGraph (nil);
    best_score = Compute();     // best over all node orders (if we're reshuffling)


    for (long iter = 0; iter < num_restarts; iter++) {
        next_score = Compute();     // reset to empty graph score


        for (long child = 0; child < num_nodes; child++) {
            long        num_parents = 0,
                        improvement_flag = 0,
                        next_parent_to_add;

            do {
                for (long parent = 0; parent < num_nodes; parent++) {
                    if ( (parent != child)                      // must meet all conditions!
                            && (dag(parent, child) == 0)
                            && (orderMx (parent, child) == 1)
                            && banned_edges(parent, child) == 0) {
                        dag.Store (parent, child, 1);
                        this_score = Compute();

                        if (this_score > next_score) {
                            improvement_flag    = 1;
                            next_score          = this_score;
                            next_parent_to_add  = parent;
                        }
                        dag.Store (parent, child, 0);   // revert
                    }
                }

                if (improvement_flag) { // adding another parent improves network score
                    dag.Store (next_parent_to_add, child, 1);
                    num_parents++;
                    improvement_flag = 0;   // reset for next parent
                } else {
                    break;  // unable to improve further
                }
            } while (num_parents < max_parents.lData[child]);
        }


        if (reshuffle_order) {
            this_score = Compute();

            if (this_score > best_score) {
                best_score = this_score;
                best_dag = (_Matrix &) dag;     // store graph optimized from the current ordering
            }

            ResetGraph (nil);

            best_node_order.Permute (1);

            for (long i = 0; i < num_nodes; i++) {
                long    child = best_node_order.lData[i];
                for (long j = 0; j < num_nodes; j++) {
                    long    parent = best_node_order.lData[j];
                    orderMx.Store (parent, child, (j > i) ? 1 : 0);
                }
            }
        } else {
            break;  // only one iteration when node ordering is known
        }
    }

    if (reshuffle_order) {
        dag = (_Matrix &) best_dag;
        best_node_order.Clear();    // reset to initial state
    }


    _Matrix * result = new _Matrix(num_nodes * num_nodes, 2, false, true);
    result->Store (0, 0, Compute());

    for (long row = 0; row < num_nodes; row++) {
        for (long col = 0; col < num_nodes; col++) {
            result->Store (row*num_nodes+col, 1, dag(row, col));
        }
    }

    return (_Matrix*)(result->makeDynamic());
}



//___________________________________________________________________________________________
_Matrix *   Bgm::GraphMCMC (bool fixed_order)
{
    _Matrix     *   proposed_graph  = new _Matrix (num_nodes, num_nodes, false, true),
    *   orderMx         = new _Matrix (num_nodes, num_nodes, false, true);

    _Matrix         current_graph (num_nodes, num_nodes, false, true),
                    best_graph (num_nodes, num_nodes, false, true);

    _Parameter      mcmc_steps, mcmc_burnin, mcmc_samples,
                    current_score, proposed_score, best_score,
                    num_randomize,
                    lk_ratio;

    long            sampling_interval;

    _SimpleList *   proposed_order  = new _SimpleList();
    _SimpleList     current_order;


    // parse HBL settings
    checkParameter (mcmcSteps, mcmc_steps, 0);
    if (mcmc_steps == 0) {
        _String oops ("You asked HyPhy to run MCMC with zero steps in the chain! Did you forget to set Bgm_MCMC_STEPS?\n");
        WarnError (oops);
    }

    checkParameter (mcmcBurnin, mcmc_burnin, 0);
    checkParameter (mcmcSamples, mcmc_samples, 0);
    if (mcmc_samples == 0) {
        _String oops ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
        WarnError (oops);
    }

    sampling_interval = (long) mcmc_steps / (long) mcmc_samples;

    checkParameter (numRandomize, num_randomize, num_nodes*num_nodes);


    //  Contents:   (1) chain trace;
    //              (2) model-averaged edge probabilities;
    //              (3) best graph;
    _Matrix     * result;

    checkPointer (result = new _Matrix ( (mcmc_samples > num_nodes*num_nodes) ? mcmc_samples : num_nodes*num_nodes, 3, false, true));



    //  Set current graph to member object (dag, an empty graph by default)
    ResetGraph (proposed_graph);

    //  does this graph conform to enforced/banned edges, maximum parentage settings?
    _SimpleList pars;
    for (long chi = 0; chi < num_nodes; chi++) {
        pars.Clear();
        for (long par = 0; par < num_nodes; par++) {
            if ( (*proposed_graph)(par,chi) == 1) {
                if ( banned_edges(par,chi) ) {
                    proposed_graph->Store(par, chi, 0);
                    ReportWarning(_String("Deleted banned edge ") & par & "->" & chi & " from graph");
                }
                pars << par;
            } else if ( (*proposed_graph)(par,chi) == 0 && enforced_edges(par,chi) ) {
                proposed_graph->Store(par, chi, 1);
                ReportWarning(_String("Restored enforced edge ") & par & "->" & chi & " to graph");
            }
        }

        if (pars.lLength > max_parents.lData[chi]) {
            ReportWarning(_String("Number of parents (") & (long)pars.lLength & ")exceed maximum setting for child node " & chi & ": " & max_parents.lData[chi] & "\n");
            pars.Permute(1);
            for (long p = 0; p < (long)pars.lLength - max_parents.lData[chi]; p++) {
                proposed_graph->Store(pars.lData[p], chi, 0);
                ReportWarning(_String("Deleted excess edge ") & pars.lData[p] & "->" & chi & " from graph");
            }
        }
    }



    if (fixed_order) {
        // coerce graph to conform to user-defined node order
        if (best_node_order.lLength == 0) {
            _String oops ("Cannot run fixed-order graph MCMC without defined order (via order-MCMC). Run CovarianceMatrix(receptacle, BGM) first.");
            WarnError (oops);

            result = new _Matrix();     // return an empty matrix
        } else {
            proposed_order->Populate(num_nodes, 0, 0);  // allocate memory
            for (long i = 0; i < num_nodes; i++) {
                proposed_order->lData[i] = best_node_order.lData[i];
            }

            ReportWarning (_String("Transferring best node order to proposal: ") & (_String *)proposed_order->toStr());
        }
    } else {
        if (best_node_order.lLength == 0) {
            // quickly generate a node order from graph
            for (long order_index, onode, node = 0; node < num_nodes; node++) {
                // locate the lowest-ranked parent in the order, if any
                for (order_index = 0; order_index < proposed_order->lLength; order_index++) {
                    onode = proposed_order->lData[order_index];
                    if ((*proposed_graph) (onode, node)) {
                        // insert new node immediately left of this parent
                        proposed_order->InsertElement ((BaseRef) node, order_index, false, false);
                        break;
                    }
                }
                // if reached end of order list without finding a parent, push node onto start of list
                if (order_index == proposed_order->lLength) {
                    proposed_order->InsertElement ((BaseRef) node, 0, false, false);
                }
            }
            ReportWarning(_String("Constructed node order from graph:\n") & (_String *)proposed_order->toStr() & "\n");
        } else {
            // restore order
            proposed_order->Populate(num_nodes,0,1);
            for (long i = 0; i < num_nodes; i++) {
                proposed_order->lData[i] = best_node_order.lData[i];
            }
        }
    }


    // randomize graph
    if (num_randomize > 0) {
        RandomizeGraph (proposed_graph, proposed_order, (long)num_randomize, fixed_order);
    } else {
        ReportWarning(_String("NUM_RANDOMIZE set to 0, skipping initial RandomizeGraph()\n"));
    }



    //  Populate order matrix -- and while we're at it, set proposed order to current order
    for (long i = 0; i < num_nodes; i++) {
        long    child = proposed_order->lData[i];

        for (long j = 0; j < num_nodes; j++) {
            //  if A precedes B, then set entry (A,B) == 1, permitting edge A->B
            orderMx->Store (proposed_order->lData[j], child, (j > i) ? 1 : 0);
        }
    }
    ReportWarning(_String("Populated order matrix\n"));


    // status line
#if !defined __UNIX__ || defined __HEADLESS__
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
    best_score = proposed_score = current_score = Compute(proposed_graph);

    ReportWarning(_String("Initiating MCMC with graph:\n") & (_String *) current_graph.toStr() );


    // MAIN LOOP
    for (long step = 0; step < mcmc_steps + mcmc_burnin; step++) {
        RandomizeGraph (proposed_graph, proposed_order, 1, fixed_order);

        proposed_score = Compute(proposed_graph);


#ifdef __DEBUG_GMCMC__
        snprintf (bug, sizeof(bug), "Current score = %f\n", current_score);
        BufferToConsole (bug);
        snprintf (bug, sizeof(bug), "Propose graph:\n");
        BufferToConsole (bug);
        PrintGraph (proposed_graph);
        snprintf (bug, sizeof(bug), "Proposed score = %f\n", proposed_score);
        BufferToConsole (bug);
#endif


        lk_ratio = exp(proposed_score - current_score);

        if (lk_ratio > 1. || genrand_real2() < lk_ratio) {  // Metropolis-Hastings
#ifdef __DEBUG_GMCMC__
            snprintf (bug, sizeof(bug), "Accept move\n");
            BufferToConsole (bug);
#endif
            current_graph       = (_Matrix &) (*proposed_graph);        // accept proposed graph
            current_score   = proposed_score;

            for (long foo = 0; foo < num_nodes; foo++) {
                current_order.lData[foo] = proposed_order->lData[foo];
            }

            if (current_score > best_score) {
                best_score = current_score;
                best_graph = (_Matrix &) current_graph;         // keep track of best graph ever visited by chain

#ifdef __DEBUG_GMCMC__
                PrintGraph (&best_graph);

                snprintf (bug, sizeof(bug), "Update best node ordering: ");
                BufferToConsole (bug);
                for (long deb = 0; deb < num_nodes; deb++) {
                    snprintf (bug, sizeof(bug), "%d ", current_order.lData[deb]);
                    BufferToConsole (bug);
                }
                NLToConsole ();
#endif
            }
        } else {
#ifdef __DEBUG_GMCMC__
            snprintf (bug, sizeof(bug), "Reject move\n");
            BufferToConsole (bug);
#endif
            // revert proposal
            for (long row = 0; row < num_nodes; row++) {
                proposed_order->lData[row] = current_order.lData[row];

                for (long col = 0; col < num_nodes; col++) {
                    proposed_graph->Store(row, col, current_graph(row, col));
                }
            }
        }



        if (step >= mcmc_burnin) {  // handle output
            if ( (step-(long)mcmc_burnin) % sampling_interval == 0) {
                long    entry   = (long int) ((step-mcmc_burnin) / sampling_interval);

                result->Store (entry, 0, current_score);        // update chain trace

                for (long row = 0; row < num_nodes; row++) {    // update edge tallies
                    for (long offset=row*num_nodes, col = 0; col < num_nodes; col++) {
                        // row = parent, col = child
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
        if (step % (long)((mcmc_steps + mcmc_burnin)/100) == 0) {
            ReportWarning (_String ("GraphMCMC at step ") & step & " of " & (mcmc_steps+mcmc_burnin) & " with posterior " & current_score & " and graph:\n" & (_String *)current_graph.toStr() & " and order:\n" & (_String *)current_order.toStr() & "\n");
        }
#endif
    }

    // convert edge tallies to frequencies, and report best graph
    for (long row = 0; row < num_nodes; row++) {
        for (long offset=row*num_nodes, col = 0; col < num_nodes; col++) {
            result->Store (offset+col, 1, (*result)(offset+col,1)/mcmc_samples);
            result->Store (offset+col, 2, best_graph(row, col));
        }
    }


    // set dag member to last visit graph
    dag = (_Matrix &) current_graph;
    best_node_order = current_order;

    delete (proposed_graph);
    delete (orderMx);
    delete (proposed_order);

    return (_Matrix *) (result->makeDynamic());
}



//___________________________________________________________________________________________
//  DEPRECATED
_Parameter  Bgm::TryEdge (long child, long parent, long operation, _Parameter old_score)
{
#ifdef _NEVER_DEFINED_
    _Parameter  log_score,
                last_state1 = dag (child, parent),
                last_state2 = dag (parent, child);

    /*
    if (is_marginal)    // modification affects only child node, don't recalculate other scores
        log_marginal = is_discrete(child,0) ? ComputeDiscreteScore (child) : ComputeContinuousScore (child);
     */

    switch (operation) {
    case 0:
        if (MarginalSum (FALSE, child) < max_parents && banned_edges(parent,child) == 0) {
            dag.Store (parent, child, 1.);  // add an edge [parent] --> [child]
            dag.Store (child, parent, 0.);
        } else {
            return old_score;
        }
        break;

    case 1: // reverse an edge
        if (MarginalSum (FALSE, parent) < max_parents && banned_edges(child,parent) == 0) {
            dag.Store (child, parent, last_state2); // switch states
            dag.Store (parent, child, last_state1);
        } else {
            return old_score;
        }
        break;

    case 2: // delete an edge
        dag.Store (child, parent, 0.);
        dag.Store (parent, child, 0.);
        break;

    default:
        return old_score;
        break;
    }

    // PrintGraph(nil);

    if (IsCyclic()) {   // new graph is cyclic, reject
        dag.Store (child, parent, last_state1);
        dag.Store (parent, child, last_state2);
        return old_score;
    }


    /*
    if (is_marginal)        // update one term only
        log_score = old_score - log_marginal + (is_discrete(child,0) ? ComputeDiscreteScore (child) : ComputeContinuousScore (child));
    else
        log_score = Compute();
     */


    log_score = Compute();

    if (log_score > old_score) {    // keep new graph
        return log_score;
    } else {
        // restore previous settings
        dag.Store (child, parent, last_state1);
        dag.Store (parent, child, last_state2);
    }

    return old_score;
#endif
    return 0;
}




//___________________________________________________________________________________________
//  Calculating equation (8) of Friedman and Koller (2003) requires the term:
//      log ( Sum ( score(x_i, U | D) ) )
//  where x_i are very small numbers stored as log(x_i)'s but taking the exponential of
//  the log values can result in numerical underflow.

_Parameter      Bgm::LogSumExpo (_GrowingVector * log_values)
{
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


//___________________________________________________________________________________________
//  Markov chain Monte Carlo for Bgm
_PMathObj Bgm::CovarianceMatrix (_SimpleList * unused)
{
    //char  bug [255];


    _Parameter      mcmc_steps,
                    mcmc_burnin,
                    mcmc_samples,
                    mcmc_nchains,
                    mcmc_temp;

    _SimpleList *   node_order;

    if (last_node_order.lLength == 0) {
        node_order = new _SimpleList (num_nodes, 0, 1);     // populate with series (0, 1, ..., num_nodes)


        node_order->Permute(1);     // initialize with randomized node ordering

    } else {
        node_order = new _SimpleList (last_node_order);
    }



    _Matrix *       mcmc_output;


    // acquisition of HBL arguments with a few sanity checks
    checkParameter (mcmcSteps, mcmc_steps, 0);
    checkParameter (mcmcBurnin, mcmc_burnin, 0);
    checkParameter (mcmcSamples, mcmc_samples, 0);
    checkParameter (mcmcNumChains, mcmc_nchains, 1);
    checkParameter (mcmcTemperature, mcmc_temp, 1);

    if (mcmc_steps == 0) {
        _String oops ("You asked HyPhy to run MCMC with zero steps in the chain! Did you forget to set Bgm_MCMC_STEPS?\n");
        WarnError (oops);
    }

    if (mcmc_samples == 0) {
        _String oops ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
        WarnError (oops);
    }



    // allocate space for output
    if (mcmc_nchains == 1) {
        mcmc_output = RunColdChain (node_order, (long) mcmc_burnin, -1);    // negative argument indicates no sampling
        delete (mcmc_output);       // discard burn-in

        mcmc_output = RunColdChain (node_order, (long) mcmc_steps, (long int) (mcmc_steps / mcmc_samples) );
    } else {    // run coupled MCMC
        node_order->Permute(1);
        // RunHotChain (node_order, (long)

        /* STILL UNDER DEVELOPMENT! */
    }

    DeleteObject (node_order);

    return mcmc_output;
}



//___________________________________________________________________________________________
_Matrix *   Bgm::RunColdChain (_SimpleList * current_order, long nsteps, long sample_lag)
{
    /*  Execute Metropolis sampler using a swap of two nodes in an ordered sequence
        as a proposal function.  The posterior probabilities of edges in the network are
        stored as a member matrix [edge_posteriors].  Note that the total number of
        possible orderings (i.e. permutations of a sequence of length N) is factorial(N),
        and can possibly be computed exactly for N < 8.                                     */

#ifdef DEBUG_RCC
    char        buf [255];
#endif

    long        row,
                first_node, second_node,
                mcmc_samples = (long int) (nsteps / sample_lag);

    _SimpleList proposed_order;

    _Matrix *   mcmc_output     = new _Matrix (mcmc_samples > (num_nodes*num_nodes) ? mcmc_samples : (num_nodes*num_nodes), 4, false, true);
    //  Contents:   (1) chain trace; (2) edge marginal posterior probabilities;
    //              (3) best node ordering; (4) last node ordering (for restarting chains).

    _List   *   clist           = new _List ();     // compute list


    InitComputeLists (clist);       // storage for marginal node- and edge-scores



    _Parameter  lk_ratio,
                prob_current_order  = Compute (current_order, clist),
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
    SetStatusLine     (empty,_HYBgm_STATUS_LINE_MCMC & (sample_lag < 0 ? _String(" burnin"):empty),empty,0,HY_SL_TASK|HY_SL_PERCENT);
#else
    _String         * progressReportFile = NULL;
    _Variable       * progressFile = CheckReceptacle (&optimizationStatusFile, empty, false);

    if (progressFile->ObjectClass () == STRING) {
        progressReportFile = ((_FString*)progressFile->Compute())->theString;
    }
    ConsoleBGMStatus (_HYBgm_STATUS_LINE_MCMC & (sample_lag < 0 ?  _String(" burnin"):empty), -1., progressReportFile);
#endif
#endif

    /* SLKP */


    best_node_order = (_SimpleList &) (*current_order);     // initialize class member variable

#ifdef DEBUG_RCC
    snprintf (buf, sizeof(buf), "Initial node order (log L = %f): ", prob_current_order);
    BufferToConsole (buf);
    for (long i = 0; i < num_nodes; i++) {
        snprintf (buf, sizeof(buf), "%d ", best_node_order.lData[i]);
        BufferToConsole (buf);
    }
    NLToConsole();
#endif

    for (long step = 0; step < nsteps; step++) {
        // copy contents of current order to proposed ordering for permutation
        proposed_order.Clear();
        proposed_order.Duplicate (current_order);


        // swap random nodes in ordered sequence
        first_node  = genrand_int32() % num_nodes;
        second_node = genrand_int32() % num_nodes;

        while (first_node == second_node) { // keep sampling until pick different node
            second_node = genrand_int32() % num_nodes;
        }

        proposed_order.Swap (first_node, second_node);


        // calculate probability of proposed order --- edge posterior probs loaded into member variable
        prob_proposed_order = Compute (&proposed_order, clist);

#ifdef DEBUG_RCC
        snprintf (buf, sizeof(buf), "Proposed node order (log L = %f): ", prob_proposed_order);
        BufferToConsole (buf);
        for (long i = 0; i < num_nodes; i++) {
            snprintf (buf, sizeof(buf), "%d ", proposed_order.lData[i]);
            BufferToConsole (buf);
        }
        NLToConsole();
#endif


        lk_ratio    = exp(prob_proposed_order - prob_current_order);

        if (lk_ratio > 1. || genrand_real2() < lk_ratio) {  // then set current to proposed order
            (*current_order).Clear();
            (*current_order) << proposed_order;
            prob_current_order = prob_proposed_order;

#ifdef DEBUG_RCC
            snprintf (buf, sizeof(buf), "Accept\n");
            BufferToConsole (buf);
#endif

            if (prob_proposed_order > best_prob) {
                best_prob = prob_proposed_order;        // update best node ordering
                best_node_order = proposed_order;
            }
        }


        if (sample_lag >= 0) {
            // write Markov chain state to output matrix

            if (step % sample_lag == 0 && step > 0) {   // AFYP 2011/3/2 added this second conditional because the zero-th step was getting recorded
                _GrowingVector *    gv;
                _Parameter          denom;

                row = (long int) (step / sample_lag) - 1;   // AFYP 2011/3/2 need to shift over other info to maintain 0-index

                // report log likelhood
                mcmc_output->Store (row, 0, prob_current_order);


#ifdef _DEBUG_RCC
                snprintf (buf, sizeof(buf), "Contents of clist:\n");
                BufferToConsole (buf);

                for (long i = 0; i < num_nodes; i++) {
                    for (long j = 0; j < num_nodes; j++) {
                        gv = (_GrowingVector *) clist->lData[i * num_nodes + j];

                        snprintf (buf, sizeof(buf), "i=%d, j=%d, n=%d: ", i, j, gv->GetUsed());
                        BufferToConsole (buf);

                        for (long k = 0; k < gv->GetUsed(); k++) {
                            snprintf (buf, sizeof(buf), "%f ", (*gv)(k,0));
                            BufferToConsole (buf);
                        }

                        NLToConsole ();
                    }
                }
#endif



                // compute marginal edge posteriors
                for (long child = 0; child < num_nodes; child++) {
                    gv      = (_GrowingVector *) clist->lData[child * num_nodes + child];
                    denom   = (*gv)(0, 0);  // first entry holds node marginal post pr


                    for (long edge, parent = 0; parent < num_nodes; parent++) {
                        if (parent == child) {
                            continue;
                        }

                        edge    = child * num_nodes + parent;
                        gv      = (_GrowingVector *) clist->lData[edge];


                        if (gv->GetUsed() > 0) {
#ifdef DEBUG_RCC
                            snprintf (buf, sizeof(buf), "edge %d: pp = %f, gv = ", edge, (*mcmc_output)(edge, 1));
                            BufferToConsole (buf);

                            for (long i = 0; i < gv->GetUsed(); i++) {
                                snprintf (buf, sizeof(buf), "%f ", (*gv) (i,0));
                                BufferToConsole (buf);
                            }

                            snprintf (buf, sizeof(buf), "sum to %f", LogSumExpo (gv));
                            BufferToConsole (buf);
#endif

                            mcmc_output->Store (edge, 1, (*mcmc_output)(edge, 1) + exp (LogSumExpo(gv) - denom));

#ifdef DEBUG_RCC
                            snprintf (buf, sizeof(buf), "= %f\n", (*mcmc_output)(edge, 1));
                            BufferToConsole (buf);
#endif
                        }


                    }
                }
            }
        }


        /*SLKP 20070926; include progress report updates */
        if (TimerDifferenceFunction(true)>1.0) { // time to update
            howManyTimesUpdated ++;
            _String statusLine = _HYBgm_STATUS_LINE_MCMC & (sample_lag < 0 ? _String(" burnin"):empty) & " " & (step+1) & "/" & nsteps
                                 & " steps (" & (1.0+step)/howManyTimesUpdated & "/second)";
#if  defined __HEADLESS__
            SetStatusLine (statusLine);
#else
#if !defined __UNIX__ || defined __HYPHYQT__ || defined __HYPHY_GTK__
            SetStatusLine     (empty,statusLine,empty,100*step/(nsteps),HY_SL_TASK|HY_SL_PERCENT);
            yieldCPUTime(); // let the GUI handle user actions
            if (terminateExecution) { // user wants to cancel the analysis
                break;
            }
#endif
#if defined __UNIX__ && ! defined __HEADLESS__ && !defined __HYPHYQT__ && !defined __HYPHY_GTK__
            ConsoleBGMStatus (statusLine, 100.*step/(nsteps), progressReportFile);
#endif
#endif

            TimerDifferenceFunction(false); // reset timer for the next second
        }
        /* SLKP */

    }
    // end loop over steps


    if (sample_lag >= 0) {
        // convert sums of edge posterior probs. in container to means
        for (long edge = 0; edge < num_nodes * num_nodes; edge++) {
            mcmc_output->Store (edge, 1, (*mcmc_output)(edge,1) / mcmc_samples);
        }
    }


    // export node ordering info
    for (long node = 0; node < num_nodes; node++) {
        mcmc_output->Store (node, 2, (_Parameter) (best_node_order.lData[node]));
        mcmc_output->Store (node, 3, (_Parameter) (current_order->lData[node]));
    }

    last_node_order = (_SimpleList &) (*current_order);



    // release cached scores (assuming that this is the last analysis to be done!)
    DumpComputeLists (clist);

    /*SLKP 20070926; include progress report updates */
#if !defined __UNIX__ || defined __HEADLESS__ || defined __HYPHYQT__ || defined __HYPHY_GTK__
    SetStatusLine     (_HYBgm_STATUS_LINE_MCMC_DONE & (sample_lag < 0 ? _String(" burnin"):empty));
#endif
#if defined __UNIX__ && ! defined __HEADLESS__ && !defined __HYPHYQT__ && !defined __HYPHY_GTK__
    ConsoleBGMStatus (_HYBgm_STATUS_LINE_MCMC_DONE & (sample_lag < 0 ? _String(" burnin"):empty), -1.0, progressReportFile);
#endif
    /* SLKP */

    return mcmc_output;
}



//___________________________________________________________________________________________
#ifdef __NEVER_DEFINED__
void    Bgm::RunHotChain (_SimpleList * current_order, long nsteps, long sample_lag, _Parameter temper)
{
    /*  Execute Metropolis sampler using a swap of two nodes in an ordered sequence
     as a proposal function.  The posterior probabilities of edges in the network are
     stored as a member matrix [edge_posteriors].  Note that the total number of
     possible orderings (i.e. permutations of a sequence of length N) is factorial(N),
     and can possibly be computed exactly for N < 8.                                        */

    _List   *   clist;

    _Parameter  prob_current_order = Compute (current_order, clist),
                prob_proposed_order = 0.;

    static long step_counter = 0;

    long        sample_interval = (long int) ((nsteps - mcmc_burnin) / mcmc_samples);

    _SimpleList proposed_order;

    bool        accept_step = FALSE;

    /* SLKP 20070926
     Add user feedback via the console window status bar
     */
#if !defined __UNIX__ || defined __HEADLESS__
    long howManyTimesUpdated = 0; // how many times has the line been updated; is the same as the # of seconds
    TimerDifferenceFunction(false); // save initial timer; will only update every 1 second
    SetStatusLine     (empty,_HYBgm_STATUS_LINE_MCMC,empty,0,HY_SL_TASK|HY_SL_PERCENT);
#endif
    /* SLKP */


    for (long step = 0; step < mcmc_steps; step++) {
        // copy contents of current order to proposed ordering for permutation
        proposed_order.Clear();
        proposed_order.Duplicate (current_order);


        // swap random nodes in ordered sequence
        long    first_node  = genrand_int32 () % proposed_order.lLength,
                second_node = genrand_int32 () % proposed_order.lLength;

        while (first_node == second_node) { // keep sampling until pick different node
            second_node = genrand_int32 () % proposed_order.lLength;
        }

        proposed_order.Swap(first_node, second_node);


        // calculate probability of proposed order --- edge posterior probs loaded into member variable
        prob_proposed_order = Compute (&proposed_order, FALSE);


        _Parameter  lk_ratio    = exp(prob_proposed_order - prob_current_order);

        accept_step = FALSE;

        if (lk_ratio > 1.)  {
            accept_step = TRUE;    // accept greater likelihood
        } else if (genrand_real2() < lk_ratio) {
            accept_step = TRUE;
        }


        if (accept_step) {  // then set current to proposed order
            (*current_order).Clear();
            (*current_order) << proposed_order;
            prob_current_order = prob_proposed_order;
        }

        // write Markov chain state to output matrix
        if (step % sample_interval == 0 && step >= mcmc_burnin) {
            long    row = (long int) ((step-mcmc_burnin) / sample_interval);

            mcmc_chain->Store (row, 0, prob_current_order);

            for (long edge = 0; edge < num_nodes*num_nodes; edge++) {
                long    child   = edge % num_nodes,
                        parent    = edge / num_nodes;

                mcmc_chain->Store (edge, 1, (*mcmc_chain)(edge, 1) + edge_posteriors(child, parent));
            }
        }

        /*SLKP 20070926; include progress report updates */
#if !defined __UNIX__ || defined __HEADLESS__
        if (TimerDifferenceFunction(true)>1.0) { // time to update
            howManyTimesUpdated ++;
            _String statusLine = _HYBgm_STATUS_LINE_MCMC & " " & (step+1) & "/" & mcmc_steps
                                 & " steps (" & (1.0+step)/howManyTimesUpdated & "/second)";
            SetStatusLine     (empty,statusLine,empty,100*step/(mcmc_steps),HY_SL_TASK|HY_SL_PERCENT);
            TimerDifferenceFunction(false); // reset timer for the next second
            yieldCPUTime(); // let the GUI handle user actions
            if (terminateExecution) { // user wants to cancel the analysis
                break;
            }
        }
#endif
        /* SLKP */
    }
    // end loop over steps


    for (long edge = 0; edge < num_nodes * num_nodes; edge++) {
        // convert sums of edge posterior probs. in container to means
        mcmc_chain->Store (edge, 1, (*mcmc_chain)(edge,1) / mcmc_samples);
    }


    // release cached scores (assuming that this is the last analysis to be done!)
    DumpComputeLists ();

    /*SLKP 20070926; include progress report updates */
#if !defined __UNIX__ || defined __HEADLESS__
    SetStatusLine     (_HYBgm_STATUS_LINE_MCMC_DONE);
#endif
    /* SLKP */

}

#endif



//___________________________________________________________________________________________
#ifdef __AFYP_DEVELOPMENT2__
void        Bgm::CacheNetworkParameters (void)
{
    _SimpleList     multipliers,    // for indexing parent combinations (j)
                    parents;

    long            num_parent_combos,
                    r_i;

    _Matrix         n_ijk,
                    theta_ijk;

    network_parameters.Clear();     // reset cache


    for (long i = 0; i < num_nodes; i++) {
        // reset variables to defaults
        multipliers.Clear();
        multipliers << 1;

        n_ijk.Clear();
        theta_ijk.Clear();

        num_parent_combos   = 1;
        r_i                 = num_levels.lData[i];


        // assemble parents based on current graph
        for (long par = 0; par < num_nodes; par++) {
            if (dag(par, i) == 1 && is_discrete.lData[par]) {
                parents << par; // list of parents will be in ascending numerical order
                num_parent_combos *= num_levels.lData[par];
                multipliers << num_parent_combos;
            }
        }


        // re-allocate space in matrix
        CreateMatrix (&n_ijk, num_parent_combos, r_i+1, false, true, false);
        CreateMatrix (&theta_ijk, num_parent_combos, r_i, false, true, false);


        // tally observations
        for (long obs = 0; obs < obsData->GetHDim(); obs++) {
            long    index           = 0,
                    child_state     = (*obsData)(obs, i);

            // is there a faster algorithm for this? - afyp
            for (long par = 0; par < parents.lLength; par++) {
                long    this_parent         = parents.lData[par],
                        this_parent_state = (*obsData)(obs, this_parent);

                index += this_parent_state * multipliers.lData[par];
            }

            n_ijk.Store ((long) index, child_state, n_ijk(index, child_state) + 1);
            n_ijk.Store ((long) index, r_i, n_ijk(index, r_i) + 1); // store n_ij sum in last entry
        }


        // compute expected network conditional probabilities (Cooper and Herskovits 1992 eq.16)
        for (long j = 0; j < num_parent_combos; j++) {
            _Parameter  denom = n_ijk(j,r_i+1) + r_i;

            for (long k = 0; k < r_i; k++) {
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
_Matrix *   Bgm::SimulateDiscreteCases (long ncases)
{
    /*  Identify nodes without parents in given graph.
        Draw random state assignments based on learned parameters (orphan score).
        Use graph matrix transpose to index child nodes.
        Visit children and randomly assign states dependent on parent assignment.
        Rinse and repeat.
    */
    CacheNetworkParameters();

    _Matrix     gad     = dag.Transpose();
    _Matrix *   cases   = new _Matrix (num_nodes, ncases, false, true);

    // proceed using postorder traversal (check for parents before resolving node)

}
#endif



//___________________________________________________________________________________________
#ifdef __AFYP_DEVELOPMENT2__
void        Bgm::PostOrderInstantiate (long node_id, _SimpleList & results)
{
    _SimpleList     parents;

    // assign parents first
    for (long par = 0; par < num_nodes; par++) {
        if (par != node_id && dag (par, node_id) == 1) {
            PostOrderInstantiate (par, results);
            parents << par;
        }
    }

    // assign random state given parents
    for (long pindex = 0; pindex < parents.lLength; pindex++) {
        long    par = parents.lData[pindex];

    }
}


long        Bgm::PredictNode (long node_id, _Matrix * assignments)
{

}
#endif


//___________________________________________________________________________________________
#if defined __HYPHYMPI__
void    Bgm::SerializeBgm (_String & rec)
{
    char        buf [255];
    _String *   bgmName = (_String *) bgmNamesList (bgmList._SimpleList::Find((long)this));
    _String     dataStr,
                dataName ("bgmData"),
                banStr,
                banName ("ban_matrix");

    _Parameter  mcem_max_steps, mcem_burnin, mcem_sample_size;

    checkParameter (mcemMaxSteps, mcem_max_steps, 10000);
    checkParameter (mcemBurnin, mcem_burnin, 1000);
    checkParameter (mcemSampleSize, mcem_sample_size, 100);

    rec << "PRINT_DIGITS=-1;\n";
    rec << "USE_MPI_CACHING=1;\n";

    // write utility functions
    rec << "function make_dnode (id,n,maxp)\n";
    rec << "{\ndnode={};\n";
    rec << "dnode[\"NodeID\"]=id;\n";
    rec << "dnode[\"PriorSize\"]=n;\n";
    rec << "dnode[\"MaxParents\"]=maxp;\n";
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

    for (long node_id = 0; node_id < num_nodes; node_id++) {
        if (is_discrete.lData[node_id]) {
            rec << "dnodes[Abs(dnodes)]=make_dnode(";
            snprintf (buf, sizeof(buf), "%ld,%ld,%ld", node_id, (long)prior_sample_size(node_id,0), (long)max_parents.lData[node_id]);
            rec << buf;
            rec << ");\n";
        } else {
            rec << "cnodes[Abs(cnodes)]=make_cnode(";
            snprintf (buf, sizeof(buf), "%ld,%ld,%ld,%f,%f", node_id, (long)prior_sample_size(node_id,0), (long)max_parents.lData[node_id],
                     prior_mean(node_id,0), prior_precision(node_id,0));
            rec << buf;
            rec << ");\n";
        }
    }

    // write BGM constructor
    rec << "BGM ";
    rec << bgmName;
    rec << "=(dnodes,cnodes);\n";

    // missing data imputation settings
    snprintf (buf, sizeof(buf), "BGM_MCEM_MAXSTEPS = %d;\n", (long)mcem_max_steps);
    rec << buf;
    snprintf (buf, sizeof(buf), "BGM_MCEM_BURNIN = %d;\n", (long)mcem_burnin);
    rec << buf;
    snprintf (buf, sizeof(buf), "BGM_MCEM_SAMPLES = %d;\n", (long)mcem_sample_size);
    rec << buf;

    // ban matrix command has to come BEFORE setting data matrix (which calls CacheNodeScores)
    rec << "ban_matrix=";
    rec << (_String *)banned_edges.toStr();
    rec << ";\n";
    rec << "SetParameter(";
    rec << bgmName;
    rec << ",BGM_BAN_MATRIX,ban_matrix);\n";

    // serialize data matrix and assign to BGM
    rec << "bgmData=";
    rec << (_String *)obsData->toStr();
    rec << ";\n";
    rec << "SetParameter(";
    rec << bgmName;
    rec << ",BGM_DATA_MATRIX,bgmData);\n";

}
#endif

//___________________________________________________________________________________________
//___________________________________________________________________________________________
//___________________________________________________________________________________________
//___________________________________________________________________________________________
//___________________________________________________________________________________________


#endif
