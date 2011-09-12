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

/* ________________________________________________________________	*/

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

/* ________________________________________________________________	*/

function make_banned_edge (parent, child)
{
	a_rule = {};
	a_rule["BanParent"] = parent;
	a_rule["BanChild"] = child;
	return a_rule;
}

/* ________________________________________________________________	*/

function checkNode (nID)
{
	if (nodeList[nID] == 0)
	{
		nodeList [nID] = 1;
		fprintf (LAST_FILE_PATH, nID, "; ");
	}
	return 0;
}

/* ________________________________________________________________	*/

function obtainBGMParameters (_lfID)
{
	ChoiceList (ambChoice, "Treatment of Ambiguities",1,SKIP_NONE,
				"Averaged","All possible resolutions are considered and averaged.",
				"Resolved","The most frequent (for that site) resolution is chosen.");

	ExecuteAFile ("Utility/AncestralMapper.bf");
	ExecuteAFile ("Utility/DescriptiveStatistics.bf");
	ExecuteAFile ("Utility/GrabBag.bf");
	
	site_map				= {};
	_SITE_RESULTS			= obtainSubstitutionMatrix ("lf",0,site_map,_OBSERVED_NS_);
	site_map				= {};
	branchCount				= Rows(_SITE_RESULTS);
	nodeCount				= Columns(_SITE_RESULTS);
	
	substitution_counts		= ({1,branchCount}["1"])*_SITE_RESULTS;
	substitution_stats		= GatherDescriptiveStats (substitution_counts);
	
	PrintDescriptiveStats   ("Counts of inferred non-synonymous substitution by site",substitution_stats);
	/* determine the appropriate lower and upper bounds */
	nontrivial_sites		= substitution_counts["_MATRIX_ELEMENT_VALUE_>0"] * ({nodeCount,1})["1"];
	
	if (nontrivial_sites[0] < 2)
	{
		fprintf (stdout, "\nERROR: BGM analysis requires at least 2 sites with non-synonymous susbtitutions\n");
		return site_map;
	}

	for (k = 0+substitution_stats["Max"]; k >= 1; k=k-1)
	{
		nontrivial_sites		= (substitution_counts["_MATRIX_ELEMENT_VALUE_>=k"]) * ({nodeCount,1})["1"];
		if (nontrivial_sites[0] >= 2)
		{
			break;
		}
	}
	
	cutoff					= prompt_for_a_value ("Include only sites with at least this many total substitutions",Max(1,substitution_stats["Median"]), Max(1,substitution_stats["Min"]), k, 1);

	for (h=0; h<nodeCount;h=h+1)
	{
		if (substitution_counts[h]>=cutoff)
		{	
			site_map[Abs(site_map)] = h;
		}
	}
	
	nodeCount = Abs (site_map);
	fprintf (stdout, "\nFound ", nodeCount, " sites with at least one non-synonymous mutation\n");
	ChoiceList  (num_parents,"Maximum parents",1,NO_SKIP,
				/*0*/ "1","Each site can be conditionally dependant on at most ONE other site. This setting permits the processing of large datasets quickly",
				/*1*/ "2","Each site can be conditionally dependant on at most TWO other sites. This setting permits the recovery of more complex dependancies, but is computationally costly. It may be too slow/memory hungry for more than 100 sites.");		
	
	if (num_parents < 0)
	{
		return site_map;
	}
	num_parents = num_parents + 1; 
	BGM_MCMC_DURATION					= prompt_for_a_value ("Run the MCMC chain for this many iterations",100000, 1000, 1e26, 1);
	BGM_MCMC_BURNIN						= prompt_for_a_value ("How many burn-in steps before the main chain is run",BGM_MCMC_DURATION$10, 100, 1e26, 1);
	BGM_MCMC_SAMPLES					= prompt_for_a_value ("Sample from the chain every so many steps",BGM_MCMC_DURATION$100, 10, BGM_MCMC_DURATION, 1);
	
	BGM_MCMC_SAMPLES					= Max(1,(BGM_MCMC_DURATION-BGM_MCMC_BURNIN)$BGM_MCMC_SAMPLES);

	ChoiceList  (resample,"Ancestral Resampling",1,NO_SKIP,
				/*0*/ "No","Base inference on the maximum likelihood ancestal reconstruction only",
				/*1*/ "Yes","In addition to maximum likelihood ancestral states, sample a number (S) of alternative ancestral reconstructions to assess robustness. Runs S additional BGM analyses [MPI Enabled]");		

	if (resample < 0)
	{
		return site_map;
	}
	if (resample > 0)
	{
		resample = prompt_for_a_value ("How many ancestral samples?",100,1,1e26,1);
	}
	
	fprintf (stdout, "\nRunning a BGM on ", nodeCount, " nodes with", 
					 "\n\t", Format(num_parents,20,0), " maximum parents per node",
					 "\n\t", Format(BGM_MCMC_BURNIN,20,0), " burn-in steps",
					 "\n\t", Format(BGM_MCMC_DURATION,20,0), " chain length",
					 "\n\t", Format(BGM_MCMC_SAMPLES,20,0), " samples\n");
			
	if (resample > 0)
	{
		fprintf (stdout, "\nWill generate ", resample, " ancestral samples\n");
	}
	_bgm_data = {};
	_bgm_data ["MAP"]				= site_map;
	_bgm_data ["MATRIX"]			= _SITE_RESULTS;				 		 
	_bgm_data ["BGM_MCMC_DURATION"] = BGM_MCMC_DURATION;				 		 
	_bgm_data ["BGM_MCMC_BURNIN"]   = BGM_MCMC_BURNIN;				 		 
	_bgm_data ["BGM_MCMC_SAMPLES"]  = BGM_MCMC_SAMPLES;				 		 
	_bgm_data ["PARENTS"]			= num_parents;				 		 
	_bgm_data ["RESAMPLE"]			= resample;				 		 
	return _bgm_data;
}

/* ________________________________________________________________	*/

function obtainSubstitutionMatrix (_lfID, sample_flag, site_map, _filterMatrix)
{
	_ancestral_id = _buildAncestralCacheInternal (_lfID, 0, _sample_flag);
	
	_fd					= _filterDimensions (_ancestral_id);
	if (Abs(site_map) == 0)
	{
		for (_k = 0; _k < _fd[0]; _k=_k+1)
		{
			site_map[_k] = _k;
		}
	}
	else
	{
		_fd[0]				= Abs (site_map);
	}
	_theMatrix			= {_fd[1],_fd[0]};
	
	for (_k = 0; _k < Abs(site_map); _k = _k+1)
	{
		_subsitution_matrix = _countSubstitutionsByBranchSite (_ancestral_id,site_map[_k],_filterMatrix);
		for (_j = 0; _j < _fd[1]; _j=_j+1)
		{
			_theMatrix [_j][_k] = _subsitution_matrix[_j];
		}
	}
	
	_destroyAncestralCache (_ancestral_id);
	return _theMatrix;
}

/* ________________________________________________________________	*/

function handleMPIBGM (_bgm_data, jobID)
{
	if (MPI_NODE_COUNT <= 1)
	{
		if (jobID >= 0)
		{				
			_sample_results [jobID] = runBGM(_bgm_data);
		}
	}
	else
	{
		
	}
}

/* ________________________________________________________________	*/

function handleMPIBGM (_bgm_data, jobID)
{
	if (MPI_NODE_COUNT <= 1)
	{
		if (jobID >= 0)
		{				
			_sample_results [jobID] = runBGM(_bgm_data);
			fprintf (stdout, "Ancestral sample ", jobID + 1, "\n");
		}
	}
	else
	{
		mpiNode = 0;
		jobToSend = "";
		if (jobID >= 0)
		{
			bgmFilePath = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "BGM.bf";
			jobToSend * 128;
			jobToSend * ("ExecuteAFile (\""+bgmFilePath+"\");");
			jobToSend * (""+_bgm_data);
			jobToSend * ("; return runBGM(_hyphyAssociativeArray);");
			jobToSend * 0;
			for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode=mpiNode+1)
			{
				if (bgm_MPI[mpiNode] < 0)
				{
					break;
				}
			}
			
		}
		doReceive = (jobID < 0) || (mpiNode == MPI_NODE_COUNT-1);
		if (doReceive)
		{
			MPIReceive		(-1, mpiNode, _jobResult);
			mpiNode			= mpiNode-1;
			receivedID		= bgm_MPI [mpiNode];
			fprintf (stdout, "Ancestral sample ", receivedID + 1, " from node ", mpiNode+1, "\n");
			ExecuteCommands ("_sample_results [" + receivedID + "] = " + _jobResult);
			bgm_MPI[mpiNode] = -1;
		}
		
		if (Abs(jobToSend))
		{
			bgm_MPI[mpiNode] = jobID;
			MPISend (mpiNode+1,jobToSend);
		}
	}
	return 0;
}

/* ________________________________________________________________	*/

function runBGM (_bgm_data)
{
	
	num_nodes			=	Abs (_bgm_data["MAP"]);
	num_parents			=	_bgm_data["PARENTS"];
	branches			=	Rows(_bgm_data["MATRIX"]);
	
	BGM_MCMC_DURATION   = _bgm_data ["BGM_MCMC_DURATION"];				 		 
	BGM_MCMC_BURNIN     = _bgm_data ["BGM_MCMC_BURNIN"];	
	BGM_MCMC_SAMPLES	= _bgm_data ["BGM_MCMC_SAMPLES"];
	
	bgm_data_matrix = {branches,num_nodes};
	
	for (k = 0; k < num_nodes; k=k+1)
	{
		i = (_bgm_data["MAP"])[k];
		for (j = 0; j < branches; j=j+1)
		{
			bgm_data_matrix[j][k] = (_bgm_data["MATRIX"])[j][i];
		}
	}

	discreteNodes   = {};
	continuousNodes = {};

	/* BGM_NTHREADS = 2; */
	num_parents = num_parents$1;

	for (k = 0; k < num_nodes; k = k+1)
	{
		discreteNodes[Abs(discreteNodes)] = make_discrete_node (k, 0, num_parents);
	}

	BGM gen_bgm = 		(discreteNodes, continuousNodes);
	SetParameter 		(gen_bgm, BGM_DATA_MATRIX, 	 bgm_data_matrix);
	SetParameter 		(gen_bgm, BGM_WEIGHT_MATRIX, num_nodes);
	CovarianceMatrix 	(postp, gen_bgm);
	return				 postp;
}
