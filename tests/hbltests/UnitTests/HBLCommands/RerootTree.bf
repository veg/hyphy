ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName ()
{
	return "Reroot Tree";
}

function runTest ()
{

    // Indexed by pre-order traversal
    Tree T1 = (((EELA:0.150276,CONGERA:0.213019):0.230956,(EELB:0.263487,CONGERB:0.202633):0.246917):0.094785,((CAVEFISH:0.451027,(GOLDFISH:0.340495,ZEBRAFISH:0.390163):0.220565):0.067778,((((((NSAM:0.008113,NARG:0.014065):0.052991,SPUN:0.061003,(SMIC:0.027806,SDIA:0.015298,SXAN:0.046873):0.046977):0.009822,(NAUR:0.081298,(SSPI:0.023876,STIE:0.013652):0.058179):0.091775):0.073346,(MVIO:0.012271,MBER:0.039798):0.178835):0.147992,((BFNKILLIFISH:0.317455,(ONIL:0.029217,XCAU:0.084388):0.201166):0.055908,THORNYHEAD:0.252481):0.061905):0.157214,LAMPFISH:0.717196,((SCABBARDA:0.189684,SCABBARDB:0.362015):0.282263,((VIPERFISH:0.318217,BLACKDRAGON:0.109912):0.123642,LOOSEJAW:0.3971):0.287152):0.140663):0.206729):0.222485,(COELACANTH:0.558103,((CLAWEDFROG:0.441842,SALAMANDER:0.299607):0.135307,((CHAMELEON:0.771665,((PIGEON:0.150909,CHICKEN:0.172733):0.082163,ZEBRAFINCH:0.099172):0.272338):0.014055,((BOVINE:0.167569,DOLPHIN:0.15745):0.104783,ELEPHANT:0.166557):0.367205):0.050892):0.114731):0.295021);
    REROOT = RerootTree(T1, "NSAM");
    Tree NEW_T1 = REROOT;
    assert(TipName(NEW_T1, TipCount(NEW_T1) - 1) =="NSAM","Rerooting a tree failed");

    // Read a large tree from a file
    tree_fn = PATH_TO_CURRENT_BF + "res/EU3031.nwk";
    fscanf(tree_fn, "String", treeString); 
    Tree T2 = treeString;

    REROOT = RerootTree(T2, "UNICORN");
    assert(REROOT=="","Rerooting a tree with non-existent tip should result in empty string");


    // TODO: Reroot trees with annotations
    // Read a tree with annotations

    //Tree T3 = (((EELA:0.150276,CONGERA:0.213019):0.230956,(EELB:0.263487,CONGERB:0.202633):0.246917):0.094785,((CAVEFISH{Foreground}:0.451027,(GOLDFISH{Foreground}:0.340495,ZEBRAFISH{Foreground}:0.390163){Foreground}:0.220565){Foreground}:0.067778,((((((NSAM{Foreground}:0.008113,NARG{Foreground}:0.014065){Foreground}:0.052991,SPUN{Foreground}:0.061003,(SMIC{Foreground}:0.027806,SDIA{Foreground}:0.015298,SXAN{Foreground}:0.046873){Foreground}:0.046977){Foreground}:0.009822,(NAUR{Foreground}:0.081298,(SSPI{Foreground}:0.023876,STIE{Foreground}:0.013652){Foreground}:0.058179){Foreground}:0.091775){Foreground}:0.073346,(MVIO{Foreground}:0.012271,MBER{Foreground}:0.039798){Foreground}:0.178835){Foreground}:0.147992,((BFNKILLIFISH{Foreground}:0.317455,(ONIL{Foreground}:0.029217,XCAU{Foreground}:0.084388){Foreground}:0.201166){Foreground}:0.055908,THORNYHEAD{Foreground}:0.252481){Foreground}:0.061905){Foreground}:0.157214,LAMPFISH{Foreground}:0.717196,((SCABBARDA{Foreground}:0.189684,SCABBARDB{Foreground}:0.362015){Foreground}:0.282263,((VIPERFISH{Foreground}:0.318217,BLACKDRAGON{Foreground}:0.109912){Foreground}:0.123642,LOOSEJAW{Foreground}:0.3971){Foreground}:0.287152){Foreground}:0.140663){Foreground}:0.206729):0.222485,(COELACANTH:0.558103,((CLAWEDFROG:0.441842,SALAMANDER:0.299607):0.135307,((CHAMELEON:0.771665,((PIGEON:0.150909,CHICKEN:0.172733):0.082163,ZEBRAFINCH:0.099172):0.272338):0.014055,((BOVINE:0.167569,DOLPHIN:0.15745):0.104783,ELEPHANT:0.166557):0.367205):0.050892):0.114731):0.295021)
    //REROOT = RerootTree(T3, "NSAM");
    //Tree NEW_T3 = REROOT;
    //assert(TipName(NEW_T3, TipCount(NEW_T3) - 1) =="NSAM","Rerooting a tree with annotations failed");

    testResult = 1;
    return testResult;
}
