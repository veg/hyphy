#include "qndhelper1.ibf";

REPLACE_TREE_STRUCTURE = 1;


ChoiceList  (rOptions,"dN/dS bias parameter options",1,NO_SKIP,
			 "Neutral","dN/dS=1",
			 "User","Custom dN/dS value",
			 "Estimate", "Estimate from data with branch corrections(slower).",
			 "Estimate + CI", "Estimate dN/dS and a 95% confidence interval from data (slower).",
			 "Estimate dN/dS only", "Estimate from data without branch corrections."
);
	

if (rOptions < 0)
{
	return;
}

if (rOptions == 1)
{
	RR = -1.0;
	while (RR<=0.0)
	{
		fprintf (stdout, "\nEnter the value for dNdS = dN/dS (>0):");
		fscanf 	(stdin,"Number",RR);
	}
	dNdS = RR;
}

ChoiceList  (cOptions,"Ancestor Counting Options",1,NO_SKIP,
			/*0*/ "Single Ancestor Counting","Use the most likely ancestor state",
			/*1*/ "Weighted Ancestor Counting","Weigh the contributions of every possible ancestral state (slow)",
			/*2*/ "Sample Ancestal States","Generate a sample of ancestral state reconstructions to assess single ancestor reliability",
			/*3*/ "Process Sampled Ancestal States","Process a previously generated sample of ancestral state reconstructions to assess single ancestor reliability",
			/*4*/ "One rate FEL","Fixed effects site-by-site likelihood estimation. dS is held constant across sites.",
			/*5*/ "Two rate FEL","Fixed effects site-by-site likelihood estimation. dS adjusted across sites.",
			/*6*/ "Rate Distribution","Obtain a site-specific distribution of dN and dS, under the assumption that E[dS] = 1. Also obtain an upper bound on rate-variation model likelihood scores.",
			/*7*/ "Full site-by-site LRT","Each site has a separate dN AND branch lengths. Experimental and VERY slow.",
			/*8*/ "Multirate FEL","Fixed effects site-by-site likelihood estimation, where non-synonymous rates are split into several classes",
			/*9*/ "BGM co-evolution","Use a Bayesian graphical model fitted to site substitution patterns to detect co-evolving sites");


if (cOptions < 0)
{
	return;
}

#include "qndhelper2.ibf";

fprintf (stdout, "\n\nPhase 4: Ancestral State Reconstruction and Counting\n\n");


pooledFreqs = {4,1};

for (k=0; k<4; k=k+1)
{
	pooledFreqs[k] = (positionFrequencies[k][0]+positionFrequencies[k][1]+positionFrequencies[k][2])/3;
}

_EFV_MATRIX0_ = {{1,AC__*pooledFreqs[1],pooledFreqs[2],AT__*pooledFreqs[3]}
				{AC__*pooledFreqs[0],1,CG__*pooledFreqs[2],CT__*pooledFreqs[3]}
				{pooledFreqs[0],CG__*pooledFreqs[3],1,GT__*pooledFreqs[3]}
				{AT__*pooledFreqs[0],CT__*pooledFreqs[3],GT__*pooledFreqs[2],1}};

_EFV_MATRIX1_ 			= _EFV_MATRIX0_;
_EFV_MATRIX2_ 			= _EFV_MATRIX0_;
useCustomCountingBias 	= 1;

if (cOptions == 9)
{
		ExecuteAFile ("Distances/CodonToolsMain.def");
		ExecuteAFile ("BGM.bf");
		_bgm_data = obtainBGMParameters ("lf");
		if (Abs(_bgm_data) == 0)
		{
			return 0;
		}
		
		ml_bgm_results			= runBGM (_bgm_data);
		
		_resamples				= _bgm_data ["RESAMPLE"];
		_sample_results			= {};
		_site_map				= _bgm_data["MAP"];
		_null_map				= {};
		for (_it = 0; _it < Abs(_site_map); _it = _it+1)
		{
			_null_map[_it] = _it;
		}
		
		bgm_MPI = {MPI_NODE_COUNT-1,1}["-1"];
		
		if (_resamples>0)
		{
			_bgm_data["MAP"] = _null_map;
			for (_it = 0; _it < _resamples; _it = _it + 1)
			{
				_bgm_data["MATRIX"]   = obtainSubstitutionMatrix("lf", 1, _site_map, _OBSERVED_NS_);
				handleMPIBGM (_bgm_data, _it);
			}
		}
												
		if (_resamples>0 && MPI_NODE_COUNT>1)
		{
			left_to_do = Transpose(bgm_MPI["1"]) * bgm_MPI["_MATRIX_ELEMENT_VALUE_>=0"];
			for (_it = 0; _it < left_to_do[0]; _it = _it + 1)
			{
				handleMPIBGM (0,-1);
			}
		}

		columnHeaders = {{"LogL"}};
		trace = ml_bgm_results[{{0,0}}][{{BGM_MCMC_SAMPLES-1,0}}];	
		traceStats = GatherDescriptiveStats (trace);
		rangeY = traceStats["Max"]-traceStats["Min"];
		if (rangeY == 0)
		{
			traceStats["Min"] = traceStats["Min"]-0.5;
			traceStats["Max"] = traceStats["Max"]+0.5;	
		}
		else
		{
			traceStats["Min"] = traceStats["Min"]-rangeY;
			traceStats["Max"] = traceStats["Max"]+rangeY;	
		}
		
		columnHeaders = {{"LogL"}};
		
		OpenWindow (CHARTWINDOW,{{"MCMC Trace"}
				{"columnHeaders"}
				{"trace"}
				{"Scatterplot"}
				{"Index"}
				{"LogL"}
				{"Sample"}
				{""}
				{"LogL"}
				{"0"}
				{""+traceStats["Mean"]}
				{"-1;-1"}
				{"10;1.309;0.785398"}
				{"Times:12:0;Times:10:0;Times:12:2"}
				{"0;0;16777215;13421772;0;0;16711680;11842740;13158600;14474460;0;3947580;16777215;0;6845928;16771158;2984993;9199669;7018159;1460610;16748822;11184810;14173291"}
				{"16,"+traceStats["Min"]+","+traceStats["Max"]}
				},
				"(SCREEN_WIDTH-30)/3;(SCREEN_HEIGHT-50)/2;10;SCREEN_HEIGHT/2+25");
				
		/* now process the edges */
		
		nodeCount			= Abs(_site_map);
		highProbEdges       = {};
		overHalf		    = {};
		edgeCount			= nodeCount*(nodeCount-1)$2;
		
		allEdges			= {edgeCount,5+3(_resamples > 0)};
		h2 = 0;
		
		SetDialogPrompt		("Spool posterior edge probabilities (as .csv) to");
		DEFAULT_FILE_SAVE_NAME = "Edge_support.csv";
		
		fprintf 			(PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN,"Site1,Site2,Site1->Site2,Site2->Site1,Site1<->Site2");
		
		if (_resamples > 0)
		{
			fprintf (LAST_FILE_PATH, ",Support for Site1->Site2,Support for Site1->Site2,Support for Site1<->Site2");
		}
		
		for (h=0; h<nodeCount; h=h+1)
		{
			for (v=h+1; v<nodeCount; v=v+1)
			{
				allEdges[h2][0] = _site_map[h]+1;
				allEdges[h2][1] = _site_map[v]+1;
				allEdges[h2][2] = ml_bgm_results[h*nodeCount+v][1];
				allEdges[h2][3] = ml_bgm_results[v*nodeCount+h][1];
				allEdges[h2][4] = allEdges[h2][2]+allEdges[h2][3];
				if (_resamples > 0)
				{
					for (_i = 0; _i < _resamples; _i = _i+1)
					{
						_ss   = {{(_sample_results[_i])[h*nodeCount+v][1],
								  (_sample_results[_i])[v*nodeCount+h][1],
								  0}};
						_ss [2] = _ss[0] + _ss[1];
						for (k = 0; k < 3; k=k+1)
						{
							allEdges[h2][5+k] = allEdges[h2][5+k]+_ss[k];
						}
					}
					for (k = 0; k < 3; k=k+1)
					{
						allEdges[h2][5+k] = allEdges[h2][5+k]/_resamples;
					}
				}
				
				if (allEdges[h2][4]>=0.95)
				{
					highProbEdges[Abs(highProbEdges)] = allEdges[h2][-1];
				}
				else
				{
					if (allEdges[h2][4]>=0.5)
					{
						overHalf[Abs(overHalf)] = allEdges[h2][-1];
					}				
				}
				fprintf (LAST_FILE_PATH,"\n",allEdges[h2][0],",",allEdges[h2][1],",",allEdges[h2][2],",",allEdges[h2][3],",",allEdges[h2][4]);
				if (_resamples > 0)
				{
					for (k=5; k<Columns(allEdges);k=k+1)
					{
						fprintf (LAST_FILE_PATH,",",allEdges[h2][k]);
					}
				}
				h2=h2+1;
			}
		}
		fprintf (LAST_FILE_PATH,CLOSE_FILE);
		
		if (Abs(highProbEdges))
		{
			fprintf (stdout, "Found ", Abs(highProbEdges), " with high (>=0.95) posterior edge support\n");
		}
		else
		{
			fprintf (stdout, "Found ", 0, " edges with high (>=0.95) posterior support\n",
							 "Reporting ", Abs(overHalf), " edges with weak (>=0.5) posterior support\n");
			highProbEdges = overHalf;
		}
		
		reportableEdges = Abs(highProbEdges);
		if (reportableEdges)
		{
			if (_resamples)
			{
				fprintf (stdout, "Site 1  Site 2 Edge 1->2 (support) Edge 2->1/ (support) Total; Edge 1<->2 (support)\n");			
			}
			else
			{
				fprintf (stdout, "Site 1  Site 2 Edge 1->2 Edge 2->1 Total; Edge 1<->2\n");
			}
			for (h=0; h<reportableEdges; h=h+1)
			{
				fprintf (stdout, Format ((highProbEdges[h])[0],6,0), " ",
								 Format ((highProbEdges[h])[1],7,0), " ");
								 
				if (Columns(highProbEdges[h])>5)
				{
					for (k=0; k<3; k=k+1)
					{
						fprintf (stdout,Format ((highProbEdges[h])[2+k],9,4), " (", Format ((highProbEdges[h])[5+k],7,3), ")");
					}
					fprintf (stdout, "\n");
				}
				else
				{
					 fprintf (stdout,Format ((highProbEdges[h])[2],9,4), " ",
									 Format ((highProbEdges[h])[3],9,4), " ",
									 Format ((highProbEdges[h])[4],18,4), "\n");
				}
			}
			
			SetDialogPrompt		("Spool a .DOT graph file (display using GraphViz)");
			DEFAULT_FILE_SAVE_NAME = "BGM.dot";
			fprintf (PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN);
			xinches = 8;
			yinches = 11;
			
			isDirected = (num_parents>1);
			
			if (isDirected)
			{
				fprintf (LAST_FILE_PATH, "di");
			}
			
			fprintf (LAST_FILE_PATH, "graph G{\nratio = \"fill\"; size = \"",xinches,",",yinches,"\"; margin = \"0,0\"; page = \"",xinches,",",yinches,"\";\n node [shape=box]; ");
			
			nodeList = {};
			edges    = "";
			edges 	 * 128;
			
			for (k=0; k < reportableEdges; k=k+1)
			{
				node1 = (highProbEdges[k])[0];
				node2 = (highProbEdges[k])[1];
				checkNode (node1);
				checkNode (node2);
				if (isDirected)
				{
					p1 = (highProbEdges[k])[2];
					p2 = (highProbEdges[k])[3];
					
					if (p2>p1)
					{	
						n = node2;
						node2 = node1;
						node1 = n;
					}
					if (_resamples == 0)
					{
						edges * ("" + node1 + " -> " + node2 + " [label = \"" + Format (p1,4,2) + "/" + Format (p2,4,2) + "\"]\n");
					}
					else
					{
						edges * ("" + node1 + " -> " + node2 + " [label = \"" + Format (p1,4,2) + "(" + Format ((highProbEdges[k])[5],4,2) + ")/" + Format (p2,4,2) + "(" + Format ((highProbEdges[k])[6],4,2) + ")\"]\n");					
					}
				}
				else
				{
					if (_resamples == 0)
					{
						edges * ("" + node1 + " -- " + node2 + " [label = \"" + Format ((highProbEdges[k])[4],4,2) + "\"]\n");
					}
					else
					{
						edges * ("" + node1 + " -- " + node2 + " [label = \"" + Format ((highProbEdges[k])[4],4,2) + "(" + Format ((highProbEdges[k])[7],4,2) + ")\"]\n");					
					}
				}	
			} 
			edges 	 * 0;
			fprintf (LAST_FILE_PATH, "\n", edges, "\n}", CLOSE_FILE);
		}
		
		
		colMin = 0;
		colMax = 1;
		numberOfBins = 50;
		matrixOfCounts = {numberOfBins, 5};
		colMax = (colMax-colMin)/numberOfBins;
	
		if (colMax==0)
		{
			colMax = 0;
		}
	
		for (counter=0; counter < edgeCount; counter = counter+1)
		{
			term = Min(((allEdges [counter][4]-colMin)/colMax)$1,numberOfBins-1);
			matrixOfCounts [term][2] = matrixOfCounts [term][2]+1;
		}
		
		term = 0;
	
		for (counter=0; counter < numberOfBins; counter = counter+1)
		{
			matrixOfCounts [counter][0] = colMin;
			term2 = matrixOfCounts [counter][2]/edgeCount;
			matrixOfCounts [counter][3] = term2;
			term = term+term2;
			matrixOfCounts [counter][4] = term;
			colMin = colMin + colMax;
			matrixOfCounts [counter][1] = colMin;
		}
	
		labelMatrix = {{"Bin Left Bound","Bin Right Bound", "Raw Count", "Bin Weight", "Cumulative Weight"}};
		promptString = "Edge Support Posterior Probabilities Histogram";
		OpenWindow (CHARTWINDOW,{{promptString}
								   {"labelMatrix"},
								   {"matrixOfCounts"},
								   {"Bar Chart"},
								   {labelMatrix[0]},
								   {labelMatrix[3]},
								   {"Posterior Edge Probability"},
								   {""},
								   {"Weight"},
									{"0"}
									{""}
									{"-1;-1"}
									{"10;1.309;0.785398"}
									{"Times:12:0;Times:10:0;Times:12:2"}
									{"0;0;16777215;1644825;0;0;6579300;11842740;13158600;14474460;0;3947580;16777215;5000268;6845928;16771158;2984993;9199669;7018159;1460610;16748822;11184810;14173291"}
									{"16,0,0"}
									},
								   "(SCREEN_WIDTH-60)/2;(SCREEN_HEIGHT-50)/2;(SCREEN_WIDTH-60)/2;50");
	return 0;
}

if (cOptions == 0)
{
	synRate = 1;
	dNdS = dNdS+1;
	cBeta = 0;
	for (h=0; h<ModelMatrixDimension; h=h+1)
	{
		cBeta = cBeta - CodonMatrix[h][h]*codonFrequencies[h];
	}
	dNdS 	= dNdS-1;
	cBeta   = cBeta-blCodon;
	cAlpha  = blCodon-dNdS*cBeta;
	
	codonT	= 0;
	v = TipCount (codonTree);
	for (h=0; h<v; h=h+1)
	{
		ExecuteCommands ("codonT=codonT+codonTree."+TipName(codonTree,h)+".synRate;");
	}
	v = BranchCount (codonTree);
	for (h=0; h<v; h=h+1)
	{
		ExecuteCommands ("codonT=codonT+codonTree."+BranchName(codonTree,h)+".synRate;");
	}
	
	
	ChoiceList  (slacOptions,"SLAC Options",1,NO_SKIP,
				 "Full tree","Analyze the entire tree",
				 "Tips vs Internals","Analyze terminal and internal branches separately");
				 
	if (slacOptions < 0)
	{	
		return;
	}

	if (slacOptions)
	{
		ExecuteCommands ("#include \"SGIvL.bf\";");
	}
	else
	{
		ExecuteCommands ("#include \"SGEmulator.bf\";");
		
		ChoiceList  (cOptions,"Rate class estimator",1,NO_SKIP,
				 "Skip","Skip the estimation of number of dS and dN rate classes",
				 "Count","Obtain approximate numbers of dS and dN rate classes supported by the data");
				 
		if (cOptions>0)
		{
			v = Rows(resultMatrix);
			scaledRates = {v,2};
			
			p  = 1/(codonT*cAlpha);
			p2 = 1/(codonT*cBeta);
			
			cd1 = 0;
			cd2 = 0;
			
			for (h=0; h<v; h=h+1)
			{
				cd1 = cd1 + resultMatrix[h][0]*p;
				cd2 = cd2 + resultMatrix[h][1]*p2;
			}
			
			p  = p*v/cd1;
			p2 = p2*dNdS*v/cd2;


			for (h=0; h<v; h=h+1)
			{
				scaledRates [h][0] = resultMatrix[h][0]*p;
				scaledRates [h][1] = resultMatrix[h][1]*p2;
			}
			
			global SMult  = 1;
			global NSMult = dNdS;
			
			dNdS := NSMult/SMult;
			
			ClearConstraints (codonTree);
			ReplicateConstraint ("this1.?.synRate:=SMult*this2.?.synRate__",codonTree,codonTree);

			h = 0.00001;
			LFCompute (lf,LF_START_COMPUTE);
			LFCompute (lf,c1);
			SMult = 1+h;
			LFCompute (lf,c2);
			SMult = 1-h;
			LFCompute (lf,c3);
			ce = -(c2+c3-2*c1)/(4*h*h)/filteredData.sites;
			SMult = 1;
			dNdS = dNdS+h;
			LFCompute (lf,c2);
			dNdS = dNdS-2*h;
			LFCompute (lf,c3);
			dNdS = dNdS+h;
			nce = -(c2+c3-2*c1)/(4*h*h)/filteredData.sites;

			fprintf (stdout, "\nApproximate synonymous Information per site: ", ce , "\n"
							 ,"Approximate non-synonymous Information per site: ", nce , "\n");
							 
			LFCompute (lf,LF_DONE_COMPUTE);
			rateClassCounter (0, ce, "synonymous");
			rateClassCounter (1, nce, "non-synonymous");
		}
	}
}
else
{
	if (cOptions == 1)
	{
		ExecuteAFile ("WANC.bf");
	}
	else
	{
		if (cOptions == 2)
		{
			fprintf (stdout, "How many ancestral samples should be generated?");
			fscanf	(stdin,"Number",sampleCount);
			SetDialogPrompt ("Save data replicates to:");
			fprintf (PROMPT_FOR_FILE,CLEAR_FILE);
			FILE_PATH = LAST_FILE_PATH;
			
			for (sampleCounter = 1; sampleCounter <= sampleCount; sampleCounter = sampleCounter + 1)
			{
				DataSet		  simAnc  = SampleAncestors (lf);
				DataSetFilter ancData = CreateFilter (simAnc,3,"","",GeneticCodeExclusions);
				multFact = ancData.sites*3*ancData.species;
				if (sampleCounter > 1)
				{
					fPath = FILE_PATH + "_" + sampleCounter;
				}
				else
				{
					fPath = FILE_PATH;
				}
				fprintf (fPath, CLEAR_FILE, ancData);
				HarvestFrequencies (nucSimF,  simAnc,1,1,0);
				fprintf (stdout, "\n\nIteration ", sampleCounter, " nucleotide composition\n",
								 "A:",nucSimF[0]*multFact, "\n",
								 "C:",nucSimF[1]*multFact, "\n",
								 "G:",nucSimF[2]*multFact, "\n",
								 "T:",nucSimF[3]*multFact, "\n");
			}
		}
		else
		{
			if (cOptions == 3)
			{
				ExecuteCommands ("#include \"SGASimProcessor.bf\";");
			}
			else
			{
				GetDataInfo    (dupInfo, filteredData);
				alreadyDone = {filteredData.unique_sites,1};
				_lfUBEST	= {filteredData.sites,1};
				
				if (cOptions == 6)
				{
					global			sFactor = 1;
					global			nFactor = 1;
					
					nFactor:<100;
					sFactor:<100;
					
					ExecuteCommands ("#include\"qndhelper3.ibf\";");

					doneSites    = {filteredData.unique_sites,4};
					fullSites    = {filteredData.sites,4};
					Tree		   siteTree = treeString;
					
					ReplicateConstraint ("this1.?.synRate   :=sFactor*this2.?.synRate__",siteTree,codonTree);
					ReplicateConstraint ("this1.?.nonSynRate:=nFactor*this2.?.synRate__",siteTree,codonTree);
					labels = {{"dS","dN","Log[L]"}};
					
					if (MPI_NODE_COUNT<=1)
					{
						for (siteCount = 0; siteCount < filteredData.sites; siteCount = siteCount+1)
						{
							siteMap = dupInfo[siteCount];
							if (alreadyDone[siteMap] == 0)
							{
								filterString = "" + (siteCount*3) + "-" + (siteCount*3+2);
								DataSetFilter siteFilter = CreateFilter (ds,3,filterString,"",GeneticCodeExclusions);
								HarvestFrequencies (f1, siteFilter, 3, 3, 0);
								
								m1 = 0;
								for (mpiNode=0; mpiNode < 64; mpiNode=mpiNode+1)
								{
									if (f1[mpiNode]>0)
									{
										m1=m1+1;
									}
								}	
								
								LikelihoodFunction siteLikelihood = (siteFilter, siteTree);
								
								if (m1>1)
								{															
									sFactor = 1;
									nFactor	= 1;
									Optimize (site_res, siteLikelihood);
									alreadyDone[siteMap]  = 1;
									doneSites[siteMap][0] = sFactor;
									doneSites[siteMap][1] = nFactor;
									doneSites[siteMap][2] = site_res[1][0];
								}
								else
								{
									sFactor = 0;
									nFactor = 0;
									LFCompute (siteLikelihood,LF_START_COMPUTE);
									LFCompute (siteLikelihood,constLF);
									LFCompute (siteLikelihood,LF_DONE_COMPUTE);									
									doneSites[siteMap][0] = sFactor;
									doneSites[siteMap][1] = nFactor;
									doneSites[siteMap][2] = constLF;
								}
							}
							dummy = ReportSite3 (siteCount, siteMap);				 
						}	
					}
					else
					{
						MPINodeState = {MPI_NODE_COUNT-1,2};
						for (siteCount = 0; siteCount < filteredData.sites; siteCount = siteCount+1)
						{
							siteMap = dupInfo[siteCount];
							if (alreadyDone[siteMap] == 0)
							{
								alreadyDone[siteMap] = 1;			
									
								filterString = "" + (siteCount*3) + "-" + (siteCount*3+2);
								DataSetFilter siteFilter = CreateFilter (filteredData,3,filterString,"",GeneticCodeExclusions);
								HarvestFrequencies (f1, siteFilter, 3, 3, 0);
								
								m1 = 0;
								for (mpiNode=0; mpiNode < 64; mpiNode=mpiNode+1)
								{
									if (f1[mpiNode]>0)
									{
										m1=m1+1;
									}
								}	
								
								LikelihoodFunction siteLikelihood = (siteFilter, siteTree);	
								
								if (m1>1)
								{			
									sFactor = 1;
									nFactor	= 1;
									
									for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
									{
										if (MPINodeState[mpiNode][0]==0)
										{
											break;	
										}
									}
									
									if (mpiNode==MPI_NODE_COUNT-1)
									{
										mpiNode = ReceiveJobs3 (1);
									}
									else
									{
										MPISend (mpiNode+1,siteLikelihood);
										MPINodeState[mpiNode][0] = 1;
										MPINodeState[mpiNode][1] = siteCount;
									}
								}
								else
								{
									sFactor = 0;
									nFactor = 0;
									LFCompute (siteLikelihood,LF_START_COMPUTE);
									LFCompute (siteLikelihood,constLF);
									LFCompute (siteLikelihood,LF_DONE_COMPUTE);									
									doneSites[siteMap][0] = sFactor;
									doneSites[siteMap][1] = nFactor;
									doneSites[siteMap][2] = constLF;
									ReportSite3 (siteCount, siteMap);								
								}
							}
						}					
						while (1)
						{
							for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
							{
								if (MPINodeState[nodeCounter][0]==1)
								{
									fromNode = ReceiveJobs3 (0);
									break;	
								}
							}
							if (nodeCounter == MPI_NODE_COUNT-1)
							{
								break;
							}
						}					
						fprintf (stdout, "\n\n\n");
						for (siteCount = 0; siteCount < filteredData.sites; siteCount = siteCount+1)
						{
							siteMap = dupInfo[siteCount];
							dummy = ReportSite3 (siteCount, siteMap);				 
						}
					}
				
					nodeCounter = 0;
					likelihoodBound = 0;

					for (siteCount = 0; siteCount < filteredData.sites; siteCount = siteCount+1)
					{
						nodeCounter 		= nodeCounter + fullSites[siteCount][0];
						likelihoodBound		= likelihoodBound + fullSites[siteCount][2];
					}
					
					nodeCounter = nodeCounter/filteredData.sites; 

					for (siteCount = 0; siteCount < filteredData.sites; siteCount = siteCount+1)
					{
						fullSites[siteCount][0] = fullSites[siteCount][0]/nodeCounter;
						fullSites[siteCount][1] = fullSites[siteCount][1]/nodeCounter;
					}
					
					OpenWindow (CHARTWINDOW,{{"Data Rates"}
											   {"labels"},
											   {"fullSites"},
											   {"Contrast Bars"},
											   {"Index"},
											   {labels[0]+";"+labels[1]},
											   {"Site Index"},
											   {labels[1]},
											   {labels[0]},
											   {"0"}},
											   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");
											   
					SetDialogPrompt ("Save rate results to:");
					
					fprintf (stdout, "\n\nApproximate likelihood upper bound = ", likelihoodBound,"\n");
					siteCount = Columns (fullSites);
					fprintf (PROMPT_FOR_FILE,CLEAR_FILE,labels[0]);
					for (nodeCounter=1; nodeCounter<siteCount; nodeCounter=nodeCounter+1)
					{
						fprintf (LAST_FILE_PATH,",",labels[nodeCounter]);
					}
					
					for (nodeCounter=0; nodeCounter < Rows (fullSites); nodeCounter = nodeCounter+1)
					{
						fprintf (LAST_FILE_PATH,"\n",fullSites[nodeCounter][0]);
						for (mpiNode=1; mpiNode<siteCount; mpiNode=mpiNode+1)
						{
							fprintf (LAST_FILE_PATH,",",fullSites[nodeCounter][mpiNode]);
						}
					}									
				}
				else
				{
					pValue = 0;
					while ((pValue<=0)||(pValue>=1))
					{	
						fprintf (stdout, "\nSignificance level for Likelihood Ratio Tests (between 0 and 1)?");
						fscanf  (stdin,"Number", pValue);
					}
					fprintf (stdout, "\n");
					if (cOptions == 4 || cOptions == 7)
					{
						Tree		   siteTree = treeString;
						doneSites    = {filteredData.unique_sites,4};
						fullSites    = {filteredData.sites,4};
						if (cOptions != 7)
						{
							ReplicateConstraint ("this1.?.synRate:=this2.?.synRate__",siteTree,codonTree);						
						}
						else
						{
							OPTIMIZATION_PRECISION = 0.0001;
						}
						
						labels = {{"dN/dS","LRT","p-value","Log(L)"}};
						
						if (MPI_NODE_COUNT <= 1)
						{
							for (siteCount = 0; siteCount < filteredData.sites; siteCount = siteCount+1)
							{
								siteMap = dupInfo[siteCount];
								if (alreadyDone[siteMap] == 0)
								{
									alreadyDone[siteMap] = 1;					
									filterString = "";
									filterString = filterString + (siteCount*3) + "-" + (siteCount*3+2);
									DataSetFilter siteFilter = CreateFilter (filteredData,3,filterString,"",GeneticCodeExclusions);

									HarvestFrequencies (f1, siteFilter, 3, 3, 0);
									m1 = 0;
									for (mpiNode=0; mpiNode < 64; mpiNode=mpiNode+1)
									{
										if (f1[mpiNode]>0)
										{
											m1=m1+1;
										}
									}	
									
									if (m1>1)
									{
										LikelihoodFunction siteLikelihood = (siteFilter, siteTree);
										if (cOptions < 7)
										{
											dNdS = 1;
											LFCompute (siteLikelihood,LF_START_COMPUTE);
											LFCompute (siteLikelihood,nullLF);
											Optimize (site_res, siteLikelihood);
										}
										else
										{
											ReplicateConstraint ("this1.?.synRate:=this2.?.synRate__",siteTree,codonTree);						
											ClearConstraints (siteTree);
											dNdS := 1;
											Optimize (site_resN, siteLikelihood);
											nullLF = site_resN[1][0];
											dNdS = 1;
											Optimize (site_res, siteLikelihood);
										}
										doneSites[siteMap][0] = dNdS;
										doneSites[siteMap][1] = 2*(site_res[1][0]-nullLF);
										doneSites[siteMap][2] = (1-CChi2(doneSites[siteMap][1],1))/2;	
										doneSites[siteMap][3] = site_res[1][0];	
									}
									else
									{
										doneSites[siteMap][0] = 1;
										doneSites[siteMap][1] = 0;
										doneSites[siteMap][2] = 1.;										
										doneSites[siteMap][3] = 0.;										
									}
								}
								dummy = ReportSite1 (siteCount, siteMap);				 
							}
						}
						else
						{
							MPINodeState = {MPI_NODE_COUNT-1,3};
							for (siteCount = 0; siteCount < filteredData.sites; siteCount = siteCount+1)
							{
								siteMap = dupInfo[siteCount];
								if (alreadyDone[siteMap] == 0)
								{
									alreadyDone[siteMap] = 1;				
									filterString = "";
									filterString = filterString + (siteCount*3) + "-" + (siteCount*3+2);
									DataSetFilter siteFilter = CreateFilter (filteredData,3,filterString,"",GeneticCodeExclusions);
									
									HarvestFrequencies (f1, siteFilter, 3, 3, 0);
									m1 = 0;
									for (mpiNode=0; mpiNode < 64; mpiNode=mpiNode+1)
									{
										if (f1[mpiNode]>0)
										{
											m1=m1+1;
										}
									}	
									
									
									if (m1>1)
									{
										if (cOptions == 7)
										{
											ReplicateConstraint ("this1.?.synRate:=this2.?.synRate__",siteTree,codonTree);						
											ClearConstraints (siteTree);
										}
										LikelihoodFunction siteLikelihood = (siteFilter, siteTree);				
										for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
										{
											if (MPINodeState[mpiNode][0]==0)
											{
												break;	
											}
										}
										
										if (mpiNode==MPI_NODE_COUNT-1)
										/* all nodes busy */
										{
											if (cOptions == 7)
											{
												mpiNode = ReceiveJobs4 (1,1);										
											}
											else
											{
												mpiNode = ReceiveJobs1 (1);
											}
										}
										else
										{
											MPISend (mpiNode+1,siteLikelihood);
											MPINodeState[mpiNode][0] = 1;
											MPINodeState[mpiNode][1] = siteCount;
											MPINodeState[mpiNode][2] = 1;
										}
										
										if (cOptions == 7)
										{
											ReplicateConstraint ("this1.?.synRate:=this2.?.synRate__",siteTree,codonTree);						
											ClearConstraints (siteTree);
											dNdS := 1;
											for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
											{
												if (MPINodeState[mpiNode][0]==0)
												{
													break;	
												}
											}
											
											if (mpiNode==MPI_NODE_COUNT-1)
											/* all nodes busy */
											{
												mpiNode = ReceiveJobs4 (1,0);
											}
											else
											{
												MPISend (mpiNode+1,siteLikelihood);
												MPINodeState[mpiNode][0] = 1;
												MPINodeState[mpiNode][1] = siteCount;
												MPINodeState[mpiNode][2] = 0;
											}
											dNdS = 1;
										}
									}
									else
									{
										doneSites[siteMap][0] = 1;
										doneSites[siteMap][2] = 1;										
									}
								}
							}					
							while (1)
							{
								for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
								{
									if (MPINodeState[nodeCounter][0]==1)
									{
										if (cOptions == 7)
										{
											fromNode = ReceiveJobs4 (0,0);										
										}
										else
										{
											fromNode = ReceiveJobs1 (0);
										}
										break;	
									}
								}
								if (nodeCounter == MPI_NODE_COUNT-1)
								{
									break;
								}
							}					
							fprintf (stdout, "\n\n\n");
							for (siteCount = 0; siteCount < filteredData.sites; siteCount = siteCount+1)
							{
								siteMap = dupInfo[siteCount];
								dummy = ReportSite1 (siteCount, siteMap);				 
							}
						}
					}
					else
					{
						SHORT_MPI_RETURN = 1;
						FEL_RUN_TIMER    = Time(1);
						FEL_MASTER_TIMER = Time(0);
						
						if (cOptions == 5)
						{
							ChoiceList  (brOptions,"Branch Options",1,NO_SKIP,
												 "All","Test for non-neutral evolution on all branches",
												 "Internal Only","Test for non-neutral evolution on internal branches only",
												 "A subtree only","Test for non-neutral evolution on a given subtree only only",
												 "Custom subset","Pick the set of branches which should be tested for presence of selection");
											 
							if (brOptions < 0)
							{
								return 0;
							}
						}
					
						if (cOptions == 5)
						{
							ExecuteCommands ("#include\"qndhelper3.ibf\";");				
							global			sFactor = 1;
							global			nFactor = 1;
							
							if (brOptions > 0)
							{
								doneSites    = {filteredData.unique_sites,8};
								fullSites    = {filteredData.sites,8};						
								labels = {{"dN","dS","dN/dS","dS=dN","LRT","p-value","Full Log(L)","dN_other"}};
							}
							else
							{
								doneSites    = {filteredData.unique_sites,7};
								fullSites    = {filteredData.sites,7};
								labels = {{"dN","dS","dN/dS","dS=dN","LRT","p-value","Full Log(L)"}};
							}
							Tree		   siteTree = treeString;
							
							ReplicateConstraint ("this1.?.synRate   :=sFactor*this2.?.synRate__",siteTree,codonTree);
							if (brOptions == 1)
							{
								global nFactorOther = 1;
								ReplicateConstraint ("this1.Node?.nonSynRate:=nFactor*this2.Node?.synRate__",siteTree,codonTree);
								ReplicateConstraint ("this1.?.nonSynRate:=nFactorOther*this2.?.synRate__",siteTree,codonTree);
							}
							else
							{
								if (brOptions == 2)
								{
									mxTreeSpec = {5,1};
									mxTreeSpec [0] = "givenTree";
									mxTreeSpec [1] = "8240";
									mxTreeSpec [2] = "10,40,-10,-175,1";
									mxTreeSpec [3] = "";
									mxTreeSpec [4] = "";
									OpenWindow (TREEWINDOW, mxTreeSpec);							
									internalNodes = BranchCount(siteTree);
									choiceMatrix = {internalNodes,2};
									for (bc=0; bc<internalNodes; bc=bc+1)
									{
										choiceMatrix[bc][0] = BranchName(siteTree,bc);
										choiceMatrix[bc][1] = "Internal Branch Rooting " + siteTree[bc];
									}
									ChoiceList  (stOption,"Choose root of subtree",1,NO_SKIP,choiceMatrix);

									if (stOption < 0)
									{
										return -1;
									}

									ExecuteCommands("ReplicateConstraint(\"this1.?.nonSynRate:=nFactor*this2.?.synRate__\",siteTree."+BranchName(siteTree,stOption)+",codonTree."+BranchName(givenTree,stOption)+");");																					      
									global nFactorOther = 1;
									ReplicateConstraint ("this1.?.nonSynRate:=nFactorOther*this2.?.synRate__",siteTree,codonTree);
								}
								else
								{
									if (brOptions == 3)
									{
										mxTreeSpec = {5,1};
										mxTreeSpec [0] = "givenTree";
										mxTreeSpec [1] = "8240";
										mxTreeSpec [2] = "10,40,-10,-175,1";
										mxTreeSpec [3] = "";
										mxTreeSpec [4] = "";
										OpenWindow (TREEWINDOW, mxTreeSpec);							
										internalNodes = BranchCount(siteTree);
										leafNodes	  = TipCount (siteTree);
										choiceMatrix = {leafNodes+internalNodes,2};
										for (bc=0; bc<leafNodes; bc=bc+1)
										{
											choiceMatrix[bc][0] = TipName(siteTree,bc);
											choiceMatrix[bc][1] = "Taxon " + choiceMatrix[bc][0];
										}
										for (bc=0; bc<internalNodes; bc=bc+1)
										{
											choiceMatrix[bc+leafNodes][0] = BranchName(siteTree,bc);
											choiceMatrix[bc+leafNodes][1] = "Internal Branch Rooting " + siteTree[bc];
										}
										ChoiceList  (stOption,"Choose branches to test",0,NO_SKIP,choiceMatrix);

										if (stOption[0] < 0)
										{
											return -1;
										}
										
										fprintf (stdout, "\nSite-by-site selection will be tested along the following subset of branches:\n");
										for (bc=0; bc<Columns(stOption) * Rows(stOption); bc=bc+1)
										{
											bc2 = stOption[bc];
											bc2 = choiceMatrix[bc2][0];
											ExecuteCommands("siteTree."+bc2+".nonSynRate:=nFactor*siteTree."+bc2+".synRate;");		
											fprintf (stdout, bc2,"\n");																		      
										}

										global nFactorOther = 1;
										ReplicateConstraint ("this1.?.nonSynRate:=nFactorOther*this2.?.synRate__",siteTree,codonTree);
											
									}
									else
									{
										ReplicateConstraint ("this1.?.nonSynRate:=nFactor*this2.?.synRate__",siteTree,codonTree);											
									}											
								}
							}
						}
						else
						{
							ExecuteCommands ("#include\"qndhelper4.ibf\";");	
							brOptions = 1;									
							global			sFactor      = 1;
							global			nFactor      = 1;
							global 			nFactorOther = 1;
							
							if (brOptions > 0)
							{
								doneSites    = {filteredData.unique_sites,8};
								fullSites    = {filteredData.sites,8};						
								labels = {{"dN","dS","dN/dS","dS=dN","LRT","p-value","Full Log(L)","dN_other"}};
							}

							Tree		   siteTree = treeString;
							
							ReplicateConstraint ("this1.?.synRate   :=sFactor*this2.?.synRate__",siteTree,codonTree);
							ExecuteCommands ("ReplicateConstraint(\"this1.?."+_rateLabelToTest+":=nFactor*this2.?.synRate__\",siteTree,codonTree)");
							for (k=0; k<Abs (aaRateClassIDs)-1; k=k+1)
							{
								ReplicateConstraint ("this1.?.?:=nFactorOther*this2.?.synRate__",siteTree,codonTree);
							}
						}
					
						if (MPI_NODE_COUNT<=1)
						{
							for (siteCount = 0; siteCount < filteredData.sites; siteCount = siteCount+1)
							{
								siteMap = dupInfo[siteCount];
								if (alreadyDone[siteMap] == 0)
								{
									alreadyDone[siteMap] = 1;					
									filterString = "";
									filterString = filterString + (siteCount*3) + "-" + (siteCount*3+2);
									DataSetFilter siteFilter = CreateFilter (ds,3,filterString,"",GeneticCodeExclusions);
									
									HarvestFrequencies (f1, siteFilter, 3, 3, 0);
									m1 = 0;
									for (mpiNode=0; mpiNode < 64; mpiNode=mpiNode+1)
									{
										if (f1[mpiNode]>0)
										{
											m1=m1+1;
										}
									}
									
									if (m1>1)
									{
										LikelihoodFunction siteLikelihood = (siteFilter, siteTree);
										sFactor      = 1;
										nFactor	     = 1;
										nFactorOther = 1;
										Optimize (site_res, siteLikelihood);
										doneSites[siteMap][0] = nFactor;
										doneSites[siteMap][1] = sFactor;
										m1 = sFactor+nFactor;
										sFactor 	= m1/2;
										nFactor := sFactor;
										Optimize (site_resN, siteLikelihood);
										
										doneSites[siteMap][2] = sFactor;
										doneSites[siteMap][3] = 2*(site_res[1][0]-site_resN[1][0]);
										doneSites[siteMap][4] = 1-CChi2(doneSites[siteMap][3],1);															
										doneSites[siteMap][5] = site_res[1][0];															
										doneSites[siteMap][6] = doneSites[siteMap][0]/doneSites[siteMap][1];	
										if (brOptions > 0)
										{
											doneSites[siteMap][7] = nFactorOther;
										}													
									}
									else
									{
										doneSites[siteMap][0] = 0;
										doneSites[siteMap][1] = 0;
										doneSites[siteMap][2] = 0;
										doneSites[siteMap][3] = 0;
										doneSites[siteMap][4] = 1;															
										doneSites[siteMap][5] = 0;
										if (brOptions > 0)
										{
											doneSites[siteMap][7] = 0;
										}
									}	
								}
								dummy = ReportSite2 (siteCount, siteMap);				 
							}	
						}
						else
						{
							MPINodeState = {MPI_NODE_COUNT-1,4};
							
							lfSpawnDone = 0;
							for (siteCount = 0; siteCount < filteredData.sites; siteCount = siteCount+1)
							{
								siteMap = dupInfo[siteCount];
								if (alreadyDone[siteMap] == 0)
								{
									filterString = "";
									filterString = filterString + (siteCount*3) + "-" + (siteCount*3+2);
									DataSetFilter siteFilter = CreateFilter (filteredData,3,filterString,"",GeneticCodeExclusions);

									HarvestFrequencies (f1, siteFilter, 3, 3, 0);
									m1 = 0;
									for (mpiNode=0; mpiNode < 64; mpiNode=mpiNode+1)
									{
										if (f1[mpiNode]>0)
										{
											m1=m1+1;
										}
									}
									
									alreadyDone[siteMap] = 1;				
									if (m1>1)
									{
										if (lfSpawnDone == 0)
										{
											LikelihoodFunction siteLikelihood = (siteFilter, siteTree);	
											lfSpawnDone = 1;
										}
													
										sFactor = 1;
										nFactor	= 1;
										nFactorOther = 1;

										for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
										{
											if (MPINodeState[mpiNode][0]==0)
											{
												break;	
											}
										}
										
										if (mpiNode==MPI_NODE_COUNT-1)
										/* all nodes busy */
										{
											mpiNode = ReceiveJobs2 (1,1);
										}
										else
										{
											MPISend (mpiNode+1,siteLikelihood);
											MPINodeState[mpiNode][0] = 1;
											MPINodeState[mpiNode][1] = siteCount;
											MPINodeState[mpiNode][2] = 1;
											MPINodeState[mpiNode][3] = MPINodeState[mpiNode][3] + 1;
											
										}
										
										m1 = sFactor+nFactor;
										sFactor 	= m1/2;
										nFactor:=sFactor;

										for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
										{
											if (MPINodeState[mpiNode][0]==0)
											{
												break;	
											}
										}
										
										/*DataSetFilter siteFilter = CreateFilter (filteredData,3,filterString,"",GeneticCodeExclusions);
										LikelihoodFunction siteLikelihood = (siteFilter, siteTree);		*/
												
										if (mpiNode==MPI_NODE_COUNT-1)
										/* all nodes busy */
										{
											mpiNode = ReceiveJobs2 (1,0);
										}
										else
										{
											MPISend (mpiNode+1,siteLikelihood);
											MPINodeState[mpiNode][0] = 1;
											MPINodeState[mpiNode][1] = siteCount;
											MPINodeState[mpiNode][2] = 0;
										}
									}
									else
									{
				                        doneSites[siteMap][0] = 0;
				                        doneSites[siteMap][1] = 0;
				                        doneSites[siteMap][2] = 0;
				                        doneSites[siteMap][3] = 0;
				                        doneSites[siteMap][4] = 1;									
									}
								}
							}					
							while (1)
							{
								for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
								{
									if (MPINodeState[nodeCounter][0]==1)
									{
										fromNode = ReceiveJobs2 (0,0);
										break;	
									}
								}
								if (nodeCounter == MPI_NODE_COUNT-1)
								{
									break;
								}
							}					
							fprintf (stdout, "\n\n\n");
							for (siteCount = 0; siteCount < filteredData.sites; siteCount = siteCount+1)
							{
								siteMap = dupInfo[siteCount];
								dummy = ReportSite2 (siteCount, siteMap);				 
							}
						}
					}
					
					OpenWindow (CHARTWINDOW,{{"FEL Results"}
												{"labels"},
												{"fullSites"},
												{"Bar Chart"},
												{"Index"},
												{labels[0]},
												{"Site Index"},
												{""},
												{labels[0]},
												{"0"}
												{""}
												{"-1;-1"}
												{"10;1.309;0.785398"}
												{"Times:12:0;Times:10:0;Times:12:2"}
												{"0;0;16777215;1644825;0;0;6579300;11842740;13158600;14474460;0;3947580;16777215;5000268;6845928;16771158;2984993;9199669;7018159;1460610;16748822;11184810;14173291"}
												{"16,0,0"}
												},
												"(SCREEN_WIDTH-60)/2;(SCREEN_HEIGHT-50)/2;(SCREEN_WIDTH-60)/2;50");
											   
					SHORT_MPI_RETURN = 0;
											   
					fprintf (stdout, "\n\nWall-clock run-time (seconds): ", Time(1) - FEL_RUN_TIMER, 
									 "\nProcess time on the master node: ", Time(0) - FEL_MASTER_TIMER,
									 "\nMPI nodes used: ", MPI_NODE_COUNT, "\n",
									 "\nWall-clock time per site (seconds): ", (Time(1) - FEL_RUN_TIMER)/filteredData.unique_sites, "\n");
									 
							
					/*if (MPI_NODE_COUNT>1)
					{		 
						unitVector = {1,MPI_NODE_COUNT-1};
						unitVector = unitVector["1"]*MPINodeState;				 
						for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
						{
							fprintf (stdout, "Node ", Format(nodeCounter+1,5,0), " processed: ", Format(MPINodeState[nodeCounter][3],5,0), 
												   "(", Format(MPINodeState[nodeCounter][3]*100/unitVector[3],6,3),"%)  jobs\n"); 
						}
					}*/
					
					SetDialogPrompt ("Save site-by-site LRT results to:");
					
					siteCount = Columns (fullSites);
					fprintf (PROMPT_FOR_FILE,CLEAR_FILE,labels[0]);
					
					outString = "";
					outString * 8192;
					
					for (nodeCounter=1; nodeCounter<siteCount; nodeCounter=nodeCounter+1)
					{
						outString * ("," + labels[nodeCounter]);
					}
					
					for (nodeCounter=0; nodeCounter < Rows (fullSites); nodeCounter = nodeCounter+1)
					{
						outString * ("\n"+fullSites[nodeCounter][0]);
						for (mpiNode=1; mpiNode<siteCount; mpiNode=mpiNode+1)
						{
							outString * (","+fullSites[nodeCounter][mpiNode]);
						}
					}
					outString * 0;
					fprintf (LAST_FILE_PATH, outString);
					outString = 0;
				}
			}
		}
	}
}

useCustomCountingBias = 0;

/*****************************************************************************/

function  DoNCategoryFitK (weightMatrix)
{
	breakPointMatrix = (SELECTED_CHART_DATA == weightMatrix);
	return breakPointMatrix[3][0];
}

/*****************************************************************************/

function rateClassCounter (cIdx, varC, descriptor)
{
	rateMx = scaledRates%cIdx;

	dataPoints = Rows(rateMx);
	temp_data_vector = {2, dataPoints};
	
	currentIndex = 0;

	temp_data_vector[0][0] = rateMx[0][cIdx];
	temp_data_vector[1][0] = 1;

	for (nextIndex = 1; nextIndex < dataPoints; nextIndex = nextIndex + 1)
	{
		if (rateMx[nextIndex][cIdx]!=rateMx[nextIndex-1][cIdx])
		{
			currentIndex = currentIndex+1;
			temp_data_vector[0][currentIndex] = rateMx[nextIndex][cIdx];
		}
		temp_data_vector[1][currentIndex] = temp_data_vector[1][currentIndex] + 1;
	}

	SELECTED_CHART_DATA = {2, currentIndex+1};

	for (nextIndex = 0; nextIndex <= currentIndex; nextIndex = nextIndex + 1)
	{
		SELECTED_CHART_DATA[0][nextIndex] = temp_data_vector[0][nextIndex];
		SELECTED_CHART_DATA[1][nextIndex] = temp_data_vector[1][nextIndex];
	}


	temp_data_vector = 0;
	normal_sigma 	 = .25;
	lastMax 		 = -100000000000;

	fprintf (stdout, 	  "\nDoing a very approximate fit with ", descriptor, " rates .\n");

	PROFILE_MEAN_VAR_MULT = 1/varC;

	for (resp = 2; resp <= Columns(SELECTED_CHART_DATA); resp = resp+1)
	{
		freqStrMx    = {resp,1};
		freqStrMx[0] = "APS_1";

		for (mi=1; mi<resp-1; mi=mi+1)
		{
			freqStrMx[mi] = "";
			for (mi2=1;mi2<=mi;mi2=mi2+1)
			{
				freqStrMx[mi] = freqStrMx[mi]+"(1-APS_"+mi2+")";		
			}
			freqStrMx[mi] = freqStrMx[mi]+"APS_"+(mi+1);	
		}	

		freqStrMx[mi] = "";
		for (mi2=1;mi2<mi;mi2=mi2+1)
		{
			freqStrMx[mi] = freqStrMx[mi]+"(1-APS_"+mi2+")";		
		}
		freqStrMx[mi] = freqStrMx[mi]+"(1-APS_"+mi+")";	
		mx = {resp,1};
		
		execString = "";
		for (resp2 = 0; resp2 < resp; resp2 = resp2 + 1)
		{
			if (resp2)
			{
					execString = execString + "global APS_"+resp2+":>0;APS_"+resp2+":<1;";
			}
			execString = execString + "mx[" + resp2 + "]:=" + freqStrMx[resp2] + ";";
		}
		ExecuteCommands (execString);
		
		
		bF = {{0}{-1e100}};
		
		for (k=0; k<10; k=k+1)
		{
			execString = "";
			rM = {resp-1,1};
			for (resp2 = 1; resp2 < resp; resp2 = resp2 + 1)
			{	
				rM [resp2-1] = Random(0,1);
				execString = execString +"APS_"+resp2+"=rM["+(resp2-1)+"];";
			}
			ExecuteCommands (execString);		
			Optimize (bestFit, DoNCategoryFitK(mx));
			if (bestFit[1][0]>bF[1][0])
			{
				bF = bestFit;
				smx = mx;
				sbp = breakPointMatrix;
			}
		}
		
		bestFit = bF;
		mx = smx;
		breakPointMatrix = sbp;
		
		/*fprintf (stdout, Format (resp,3,0), " rate classes: Log(L) = ", bestFit[1][0], " \n");*/
		mi = 0;
		for (resp2 = 0; resp2 < resp; resp2 = resp2+1)
		{
			mi = mi + breakPointMatrix[2][resp2]*mx[resp2];
		}

		cI_R = {2,resp};

		for (resp2 = 0; resp2 < resp; resp2 = resp2+1)
		{
			cI_R[0][resp2] = breakPointMatrix[2][resp2]/mi;
			cI_R[1][resp2] = mx[resp2];
			/*fprintf (stdout, "Class ", Format (resp2+1,3,0), " : ",  Format (cI_R[0][resp2], 8, 5), " weight = ", Format (cI_R[1][resp2], 8, 5), "\n");*/
		}	
		
		if (bestFit[1][0] - lastMax < 2)
		{
			break;
		}
		saveBPM = breakPointMatrix;
		saveMx	= mx;
		lastMax = bestFit[1][0];
	}

	resp = resp-1;

	fprintf (stdout, "\n\nFitted ", resp, " ", descriptor, " categories to the data\n\n");

	if (varC==0)
	{
		mi = 0;
		for (resp2 = 0; resp2 < resp; resp2 = resp2+1)
		{
			mi = mi + saveBPM[2][resp2]*saveMx[resp2];
		}
	}
	else
	{
		mi = mi+1;
	}

	cI_R = {2,resp};

	for (resp2 = 0; resp2 < resp; resp2 = resp2+1)
	{
		cI_R[0][resp2] = saveBPM[2][resp2]/mi;
		cI_R[1][resp2] = saveMx[resp2];
		fprintf (stdout, "Class ", Format (resp2+1,3,0), " : ",  Format (cI_R[0][resp2], 8, 5), " weight = ", Format (cI_R[1][resp2], 8, 5), "\n");
	}
	
	return resp;
}