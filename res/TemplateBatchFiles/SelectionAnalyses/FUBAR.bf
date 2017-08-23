RequireVersion  ("2.31");

LoadFunctionLibrary("libv3/all-terms.bf");
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary ("libv3/models/codon/MG_REV.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");

LoadFunctionLibrary("modules/io_functions.ibf");
LoadFunctionLibrary("modules/selection_lib.ibf");
LoadFunctionLibrary("libv3/models/codon/BS_REL.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

/*------------------------------------------------------------------------------*/

fubar.json = {
                    terms.json.input: {},
                    terms.json.fits : {},
                    terms.json.timers : {},
                    terms.json.test_results : {}
                };



/*------------------------------------------------------------------------------*/

fubar.analysis_description = {terms.io.info : "Perform a Fast Unbiased AppRoximate Bayesian (FUBAR)
analysis of a coding sequence alignment to determine whether some
sites have been subject to pervasive purifying or diversifying
selection.

Please note that a FUBAR analysis generates a cache and a results JSON file in the same directory as
directory as the original alignment. HyPhy needs to have write
privileges to this directory. For example if the original file is in
/home/sergei/FUBAR/data/pol.nex then at the end of a FUBAR run, there
will also exist files such as /home/sergei/FUBAR/data/pol.nex.json,
/home/sergei/FUBAR/data/pol.nex.cache. They also
provide checkpointing so that a partially completed analysis can be
restarted.",
           terms.io.version : "2.0",
           terms.io.reference : "FUBAR: a fast, unconstrained bayesian approximation for inferring selection (2013), Mol Biol Evol. 30(5):1196-205",
           terms.io.authors : "Sergei L Kosakovsky Pond",
           terms.io.contact : "spond@temple.edu",
           terms.io.requirements : "in-frame codon alignment (possibly partitioned) and a phylogenetic tree (one per partition)"
          };

io.DisplayAnalysisBanner ( fubar.analysis_description );

namespace fubar {
    LoadFunctionLibrary ("modules/shared-load-file.bf");
    load_file ({utility.getGlobalValue("terms.prefix"): "fubar", utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "selection.io.SelectAllBranches"}});
}


fubar.path.base = (fubar.json [terms.json.input])[terms.json.file];
fubar.path.cache = fubar.path.base + ".cache";
fubar.cache = io.LoadCacheFromFile (fubar.path.cache);



namespace terms.fubar {
    namespace cache {
        gtr          = "gtr-fit";
        grid         = "grid";
        conditionals = "conditionals";
        settings     = "settings";
     }
    namespace settings {
        grid_points = "grid points";
    }
}

console.log ( "> FUBAR will write cache and result files to _" + fubar.path.base + "_.extension\n\n");

//----------------------------------------------------------------------------
// PHASE 1: nucleotide fit
//----------------------------------------------------------------------------

if (utility.Has (fubar.cache, terms.fubar.cache.gtr, "AssociativeList")) {
     fubar.gtr_results = fubar.cache[terms.fubar.cache.gtr];
     io.ReportProgressMessageMD("fubar", "nuc-fit", "_Cached GTR fit_ " + selection.io.report_fit (fubar.gtr_results, 0, fubar.codon_data_info[terms.data.sample_size]));
}
else {
    namespace fubar {
        doGTR ("fubar");
    }
    fubar.cache[terms.fubar.cache.gtr] = fubar.gtr_results;
    io.WriteCacheToFile (fubar.path.cache, fubar.cache);
}


utility.ForEachPair (fubar.gtr_results [terms.branch_length], "_key_", "_value_",
                        '
                             _value_lengths_ = utility.Map (_value_, "_info_", "_info_[utility.getGlobalValue (\\"terms.fit.MLE\\")]");
                             io.ReportProgressMessageMD ("fubar", "nuc-fit", "* Tree length for partition " + (1 + _key_) + " : " + (+_value_lengths_));
                        ');



if (utility.Has (fubar.cache, terms.fubar.cache.grid, "Matrix") && utility.Has (fubar.cache, terms.fubar.cache.conditionals, "AssociativeList")) {

} else {
    // define FUBAR model

    fubar.mg_rev = model.generic.DefineModel("models.codon.MG_REV.ModelDescription",
        "fubar.codon_model", {
            "0": "terms.local",
            "1": fubar.codon_data_info [terms.code]
        },
        fubar.filter_names,
        None);

    fubar.mg_rev [utility.getGlobalValue("terms.model.set_branch_length")] = "fubar.scaled.set_branch_length";
    fubar.mg_rev ['approximate ratio'] = 0.5;
    fubar.model_assignment = {
        "default": fubar.mg_rev
    };
    fubar.model_id_to_object = {
        "fubar.codon_model": fubar.mg_rev
    };

    fubar.trees.names = utility.Map (utility.Range (fubar.partition_count, 1, 1), "_index_", "'fubar.codon_tree_' + _index_");
    fubar.lf.components = {fubar.partition_count * 2, 1};
    utility.ForEachPair (fubar.filter_names, "_index_", "_filter_",
        '
            fubar.lf.components [2*(0+_index_)] = _filter_;
            fubar.lf.components [2*(0+_index_) + 1] = fubar.trees.names[_index_];
            model.ApplyModelToTree(fubar.trees.names [_index_], fubar.trees[_index_], fubar.model_assignment, None);
        ');


    LikelihoodFunction fubar.lf.codon = (fubar.lf.components);


    estimators.ApplyExistingEstimates("fubar.lf.codon", fubar.model_id_to_object, fubar.gtr_results, None);
    estimators.TraverseLocalParameters ("fubar.lf.codon", fubar.model_id_to_object, "fubar.set_scale");

    fprintf (stdout, fubar.lf.codon);

    terms.fubar.settings.grid_points = io.PromptUser ("> Number of grid points per dimension (total number is D^2)",20,5,50,TRUE);
}

return 0;

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

lfunction fubar.scaled.set_branch_length (model, value, parameter) {
    io.CheckAssertion ('Type (`&value`) == "Number"', "Internal error in fubar.scaled.set_branch_length");
    local = (model[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.local")];
    alpha = local [utility.getGlobalValue("terms.parameters.synonymous_rate")];
    beta = local [utility.getGlobalValue("terms.parameters.nonsynonymous_rate")];
    parameters.SetConstraint (beta, " " + model['approximate ratio'] + "*" + alpha, "");
    utility.ExecuteInGlobalNamespace ("FindRoot (`&s`,(" + model[utility.getGlobalValue("terms.model.branch_length_string")] + ")-" + 3*value + "," + alpha + ",0,10000)");
    parameters.SetValue ("`parameter`.`alpha`", s);
    parameters.SetValue ("`parameter`.`beta`", ^beta);
}

//----------------------------------------------------------------------------
// PHASE 2: branch length scaling and grid calculation
//----------------------------------------------------------------------------

fprintf (stdout, "\n\n");

_fubarCodonFitLocation = fubar.paths["Base"] + fubar.paths["Codon fit suffix"];
_fubarGridInfoLocation = fubar.paths["Base"] + fubar.paths["Grid information"];

if (fubar.caching && !_fubarGridInfoLocation && !_fubarCodonFitLocation) {
     fprintf (stdout, "[CACHED] FUBAR found the self-contained codon fit file at ", _fubarCodonFitLocation, "\n");
     fprintf (stdout, "[CACHED] FUBAR found the site likelihoods file at ", _fubarGridInfoLocation, "\n");
}
else
{
    fubar.caching = 0;
    _grid_points = prompt_for_a_value ("Number of grid points per dimension (total number is D^2)",20,5,50,1);
    fprintf (stdout, "[DIAGNOSTIC] FUBAR will use a ", _grid_points , "X", _grid_points, " grid\n");

    ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{PATH_TO_CURRENT_BF[0][Abs(PATH_TO_CURRENT_BF)-2],"FUBAR_HBL","FUBAR_PHASE_2.bf"}}), {"0" : fubar.gtr_fit,
                                                                                                                                 "1" : _fubarCodonFitLocation,
                                                                                                                                 "2" : _fubarGridInfoLocation,
                                                                                                                                 "3" : "" + _grid_points});
    fprintf (stdout, "[DIAGNOSTIC] FUBAR wrote the self-contained codon fit file to ", _fubarCodonFitLocation, "\n");
    fprintf (stdout, "[DIAGNOSTIC] FUBAR wrote the the site likelihoods file to ", _fubarGridInfoLocation, "\n");
}

//----------------------------------------------------------------------------
// PHASE 3: MCMC
//----------------------------------------------------------------------------

fprintf (stdout, "\n\n");

_fubarMCMCSamplesLocation = fubar.paths["Base"] + fubar.paths["MCMC samples"];

if (fubar.caching && !_fubarMCMCSamplesLocation) {
     fscanf (_fubarMCMCSamplesLocation, "Number", _fubarChainCount);
     fprintf (stdout, "[CACHED] FUBAR found the MCMC samples based on ", _fubarChainCount, " chains at ", _fubarMCMCSamplesLocation, "\n");
}
else
{
    fubar.caching = 0;
    _fubarChainCount = prompt_for_a_value ("Number of MCMC chains to run",5,2,20,1);
    fprintf (stdout, "[DIAGNOSTIC] FUBAR will use run ", _fubarChainCount, " independent chains\n");
    _fubarChainLength  = prompt_for_a_value ("The length of each chain",2000000,500000,100000000,1);
    fprintf (stdout, "[DIAGNOSTIC] FUBAR will run the chains for ", _fubarChainLength, " steps\n");
    _fubarChainBurnin  = prompt_for_a_value ("Discard this many samples as burn-in",_fubarChainLength$2,_fubarChainLength$20,_fubarChainLength*95$100,1);
    fprintf (stdout, "[DIAGNOSTIC] FUBAR will run discard ", _fubarChainBurnin, " steps as burn-in\n");
    _fubarTotalSamples = prompt_for_a_value ("How many samples should be drawn from each chain",100,10,_fubarChainLength-_fubarChainBurnin,1);
    fprintf (stdout, "[DIAGNOSTIC] FUBAR will run thin each chain down to ", _fubarTotalSamples, " samples\n");
    _fubarPriorShape = prompt_for_a_value ("The concentration parameter of the Dirichlet prior",0.5,0.001,1,0);
    fprintf (stdout, "[DIAGNOSTIC] FUBAR will use the Dirichlet prior concentration parameter of ", _fubarPriorShape, "\n");

    ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{PATH_TO_CURRENT_BF[0][Abs(PATH_TO_CURRENT_BF)-2],"FUBAR_HBL","FUBAR_PHASE_3.bf"}}), {"0" : _fubarMCMCSamplesLocation,
                                                                                                                                 "1" : _fubarGridInfoLocation,
                                                                                                                                 "2" : "" + _fubarChainCount,
                                                                                                                                 "3" : "" + _fubarChainLength,
                                                                                                                                 "4" : "" + _fubarChainBurnin,
                                                                                                                                 "5" : "" + _fubarTotalSamples,
                                                                                                                                 "6" : "" + _fubarPriorShape
                                                                                                                                  });



    fprintf (stdout, "\n[DIAGNOSTIC] FUBAR wrote samples from ", _fubarChainCount, " independent chains to ", _fubarMCMCSamplesLocation, "[0-", _fubarChainCount-1, "]\n");
}

//----------------------------------------------------------------------------
// PHASE 4: PROCESSING & FDR Simulation
//----------------------------------------------------------------------------

fprintf (stdout, "\n\n");

_fubarResultLocation = fubar.paths["Base"] + fubar.paths["Output"];
_fubarSimGrid        = fubar.paths["Base"] + fubar.paths["SimGrid"];
_fubarSimFitFile     = fubar.paths["Base"] + fubar.paths["SimFitFile"];

_fubar_do_simulations = 0;

if (fubar.caching && !_fubarResultLocation && (_fubar_do_simulations == 0 || (!_fubarSimGrid && !_fubarSimFitFile))) {
     fprintf (stdout, "[CACHED] FUBAR found the results file at ",_fubarResultLocation  ,"\n");
}
else
{

    ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{PATH_TO_CURRENT_BF[0][Abs(PATH_TO_CURRENT_BF)-2],"FUBAR_HBL","FUBAR_PHASE_4.bf"}}), {"0" : fubar.gtr_fit,
                                                                                                                                 "1" : _fubarGridInfoLocation,
                                                                                                                                 "2" : _fubarMCMCSamplesLocation,
                                                                                                                                 "3" : "" + _fubarChainCount,
                                                                                                                                 "4" : _fubarResultLocation,
                                                                                                                                 "5" : _fubarSimFitFile,
                                                                                                                                 "6" : _fubarSimGrid,
                                                                                                                                 "7" : _fubarCodonFitLocation
                                                                                                                                  });



    fprintf (stdout, "\n[DIAGNOSTIC] FUBAR wrote the results of its analysis to ", _fubarResultLocation, "\n");
    if (_fubar_do_simulations) {
        fprintf (stdout, "[DIAGNOSTIC] FUBAR wrote FDR simulation data to ", _fubarSimFitFile, "\n");
        fprintf (stdout, "[DIAGNOSTIC] FUBAR wrote FDR grid information to ", _fubarSimFitFile, "\n");
    }
}

fubar_data = (ReadCSVTable (_fubarResultLocation, 1))[1]%4;

ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{PATH_TO_CURRENT_BF[0][Abs(PATH_TO_CURRENT_BF)-2],"FUBAR_HBL","FUBAR_tools.ibf"}}));

fprintf (stdout, "\n[RESULTS] At posterior probability >= 0.9 ");

idx = Rows(fubar_data);
mean_pp = 0;

p_i = {};

while (fubar_data[idx-1][4] >= 0.9 && idx > 0) {
    mean_pp += (1-fubar_data[idx-1][4]);
    p_i + (1-fubar_data[idx-1][4]);
    idx += -1;
}

if (idx == Rows(fubar_data) ) {
    fprintf (stdout, "there were no sites under diversifying positive selection\n");
} else {
    detected = Rows(fubar_data)-idx;
    ci = computeENFP_CI (p_i, 0.05);
    fprintf (stdout, "there were ", detected, " sites under diversifying positive selection, of which ", Format (mean_pp, 5,2), " [", ci[0], " - ", ci[1], "] are expected to be false positives.\n");
    _fubar_did_simulations = Columns(fubar_data) > 9;
    if (_fubar_did_simulations) {
        fprintf (stdout, "\nCodon\tProb[dN/dS>1]\tEBF[dN/dS]>1\tPSRF\tN_eff\tFDR");
        for (idx2 = Rows(fubar_data)-1; idx2 >= idx; idx2 += -1) {
            fprintf (stdout, "\n", fubar_data[idx2][0], "\t",  fubar_data[idx2][4], "\t",  fubar_data[idx2][6], "\t", fubar_data[idx2][7], "\t",  fubar_data[idx2][8], "\t",  fubar_data[idx2][9]);
        }
    } else {
        fprintf (stdout, "\nCodon\tProb[dN/dS>1]\tEBF[dN/dS]>1\tPSRF\tN_eff");
        for (idx2 = Rows(fubar_data)-1; idx2 >= idx; idx2 += -1) {
            fprintf (stdout, "\n", fubar_data[idx2][0], "\t",  fubar_data[idx2][4], "\t",  fubar_data[idx2][6], "\t", fubar_data[idx2][7], "\t",  fubar_data[idx2][8]);
        }
    }
    fprintf (stdout, "\n");
}

//------------------------------------------------------------------------------------------------//

// library functions

//------------------------------------------------------------------------------------------------//

lfunction fubar.DefineAlphaBetaGrid (one_d_points) {
    alphaBetaGrid = {one_d_points^2,2}; // (alpha, beta) pair
    oneDGrid      = {one_d_points,1};

    one_d_points    = Max (one_d_points, 10);
    neg_sel         = 0.7;
    neg_sel_points  = ((one_d_points)*neg_sel+0.5)$1;
    pos_sel_points  = (one_d_points-1)*(1-neg_sel)$1;
    if (neg_sel_points + pos_sel_points != one_d_points) {
        pos_sel_points = one_d_points - neg_sel_points;
    }
    _neg_step = 1/neg_sel_points;
    for (_k = 0; _k < neg_sel_points; _k += 1) {
        oneDGrid [_k][0] =  _neg_step * _k;
    }
    oneDGrid [neg_sel_points-1][0] = 1;
    _pos_step = 49^(1/3)/pos_sel_points;
    for (_k = 1; _k <= pos_sel_points; _k += 1) {
        oneDGrid [neg_sel_points+_k-1][0] = 1+(_pos_step*_k)^3;
    }

    _p = 0;
    for (_r = 0; _r < one_d_points; _r += 1) {
        for (_c = 0; _c < one_d_points; _c += 1) {
           alphaBetaGrid[_p][0] = oneDGrid[_r];
           alphaBetaGrid[_p][1] = oneDGrid[_c];
           _p += 1;
        }
    }

    return alphaBetaGrid;
}
