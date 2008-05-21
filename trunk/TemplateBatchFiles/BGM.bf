RequireVersion  ("0.9920060815");

/* ________________________________________________________________	*/

function make_discrete_node (node_id, sample_size, max_parents)
{
	dnode = {};
	dnode["NodeID"] = node_id;
	dnode["PriorSize"] = sample_size;
	dnode["MaxParents"] = max_parents;
	return dnode;
}

function make_continuous_node (node_id, sample_size, mean, precision)
{
	cnode = {};
	cnode["NodeID"] 	= node_id;
	cnode["PriorSize"] 	= sample_size;
	cnode["MaxParents"] = max_parents;
	cnode["PriorMean"]	= mean;
	cnode["PriorVar"]	= precision;
	return cnode;
}

function make_banned_edge (parent, child)
{
	a_rule = {};
	a_rule["BanParent"] = parent;
	a_rule["BanChild"] = child;
	return a_rule;
}

discreteNodes = {};
continuousNodes = {};
banlist = {};

/* BGM_NTHREADS = 2; */
num_parents = num_parents$1;
num_nodes   = Columns (bgm_data_matrix);


for (k = 0; k < num_nodes; k = k+1)
{
	discreteNodes[Abs(discreteNodes)] = make_discrete_node (k, 0, num_parents);
}

if (BGM_MCMC_DURATION == 0)
{
	BGM_MCMC_DURATION 	= 100000;
	BGM_MCMC_BURNIN		= BGM_MCMC_DURATION $ 10;
	BGM_MCMC_SAMPLES 	= BGM_MCMC_DURATION $ 100;
}


BGM gen_bgm = 		(discreteNodes, continuousNodes, banlist);
SetParameter 		(gen_bgm, BGM_DATA_MATRIX, 	 bgm_data_matrix);
SetParameter 		(gen_bgm, BGM_WEIGHT_MATRIX, num_nodes);
CovarianceMatrix 	(postp, gen_bgm);

nvar 				= Abs		(discreteNodes);
nobs 				= Rows 	(bgm_data_matrix);

/*
edgeSupport = {nvar,nvar} ["postp[_MATRIX_ELEMENT_ROW_*nvar+_MATRIX_ELEMENT_COLUMN_][1]"];
*/
function checkNode (nID)
{
	if (nodeList[nID] == 0)
	{
		nodeList [nID] = 1;
		fprintf (LAST_FILE_PATH, nID, "; ");
	}
	return 0;
}
