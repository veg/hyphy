TRY_NUMERIC_SEQUENCE_MATCH = 1;
DISTANCE_PROMPTS = 1;

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function PopulateModelMatrix (ModelMatrixName&, doRV)
{
	modelDefString = "";
	modelDefString*128;
	
	if (doRV)
	{
		incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"defineGamma.mdl";
		ExecuteCommands  ("#include \""+incFileName+"\";");
		blp = "c*t;";
	}	
	else
	{
		blp = "t;";
	}
	
	ModelMatrixName = {4,4};
	
	for (transition = 0; transition < 4; transition = transition + 1)
	{
		for (transition2 = transition + 1; transition2 < 4; transition2 = transition2 + 1)
		{
			modelDefString*("ModelMatrixName["+transition+"]["+transition2+"] := "+_nucBiasTerms[transition][transition2]+blp+
					        "ModelMatrixName["+transition2+"]["+transition+"] := "+_nucBiasTerms[transition2][transition]+blp);
		}
	}
	modelDefString*0;
	ExecuteCommands (modelDefString);
	return 1;
}

/* ________________________________________________________________________________________________*/

function checkEmbedding (_m1, _m2)
{
	for (r=0; r<6; r=r+1)
	{
		if (_m2[r]<_m1[r])
		{
			return 0;
		}
		if (_m2[r]>_m1[r])
		{
			for (r2 = 0; r2 < 6; r2 = r2+1)
			{
				if ((_m2[r2]==_m2[r])&&(_m1[r2]!=_m1[r]))
				{
					return 0;
				}
			}
		}
	}
	return 1;
}

/* ________________________________________________________________________________________________*/

function ReceiveJobs (sendOrNot)
{
	if (MPI_NODE_COUNT>1)
	{
		MPIReceive (-1, fromNode, result_String);
		jobModelNum = MPINodeState[fromNode-1][1];
		vv1 = MPINodeState[fromNode-1][2];
		vv2 = MPINodeState[fromNode-1][3];
		vv3 = MPINodeState[fromNode-1][4];
		vv4 = MPINodeState[fromNode-1][5];
		vv5 = MPINodeState[fromNode-1][6];
		vv6 = MPINodeState[fromNode-1][7];
		if (sendOrNot)
		{
			MPISend (fromNode,lf);
			MPINodeState[fromNode-1][1] = modelNum;		
			MPINodeState[fromNode-1][2] = v1;
			MPINodeState[fromNode-1][3] = v2;
			MPINodeState[fromNode-1][4] = v3;
			MPINodeState[fromNode-1][5] = v4;
			MPINodeState[fromNode-1][6] = v5;
			MPINodeState[fromNode-1][7] = v6;
		}
		else
		{
			MPINodeState[fromNode-1][0] = 0;
			MPINodeState[fromNode-1][1] = -1;		
		}
		
		ExecuteCommands (result_String);
	}
	else
	{
		jobModelNum = modelNum;
	}
	
	if (jobModelNum == 0)
	{
		stdl = lf_MLES[1][0];
		fullnp = lf_MLES[1][1]+totalBranchCount;
		fprintf(stdout,"\n(012345) Full Model ln-lik =  ",stdl,". Parameter Count=",Format(lf_MLES[1][1],0,0)," AIC = ", 2*(fullnp-stdl),"\n\n");


		resultCache [0][0] = 1;
		resultCache [0][1] = 2;
		resultCache [0][2] = 3;
		resultCache [0][3] = 4;
		resultCache [0][4] = 5;
		resultCache [0][5] = lf_MLES[1][0];
		resultCache [0][6] = lf_MLES[1][1]+totalBranchCount;
		resultCache [0][7] = 0;
		resultCache [0][8] = 0;
		
		fprintf (stdout,"\n#   |  Model   | # prm |    lnL    |      LRT       |    AIC     |   P-Value        |");   
		fprintf (stdout,"\n----|----------|-------|-----------|----------------|------------|------------------|");   

		if (MPI_NODE_COUNT>1)
		{
			for (h=1; h<203; h=h+1)
			{
				lnL = resultCache[h][5];
				
				if (lnL<0)
				{
					np = resultCache[h][6];
					LRT = -2*(lnL-stdl);
					if (LRT<0)
					{
						LRT = 0;
					}
					AIC = -2*lnL+2*np;
					PRINT_DIGITS = 3;
					fprintf (stdout,"\n",h);
					PRINT_DIGITS = 1;
					fprintf (stdout," | (",0, resultCache[h][0], resultCache[h][1], resultCache[h][2], resultCache[h][3], resultCache[h][4],") | ");
					fprintf (stdout,Format (np,5,0));
					PRINT_DIGITS = 8;
					fprintf (stdout, " |  ",lnL," | ",Format(LRT,14,3), " |  ", AIC, "  |  ", );
					
					PRINT_DIGITS = 15;
					if (LRT==0)
					{
						pValue = 1;					
					}
					else
					{
						pValue = 1-CChi2(LRT,fullnp-np);
					}
					fprintf (stdout,pValue," |");
					resultCache [jobModelNum][7] = pValue;
					if (pValue<rejectAt)
					{
						rejectCount = rejectCount+1;
						resultCache [jobModelNum][8] = 0;
					}
					else
					{
						resultCache [jobModelNum][8] = 1;					
					}
					
					if (pValue<rejectAt)
					{
						fprintf (stdout,"(*)");
					}				
				}
			}
		}
		
		return fromNode-1;
	}
	else
	{
		if ((MPI_NODE_COUNT>1)&&(resultCache[0][5]>=0))
		{
			resultCache [jobModelNum][0] = vv2;
			resultCache [jobModelNum][1] = vv3;
			resultCache [jobModelNum][2] = vv4;
			resultCache [jobModelNum][3] = vv5;
			resultCache [jobModelNum][4] = vv6;
			resultCache [jobModelNum][5] = lf_MLES[1][0];
			resultCache [jobModelNum][6] = lf_MLES[1][1]+totalBranchCount;
						
			return fromNode - 1;
		}
	}

	np = lf_MLES[1][1]+totalBranchCount;
	lnL = lf_MLES[1][0];
	LRT = -2*(lnL-stdl);
	if (LRT<0)
	{
		LRT = 0;
	}
	AIC = -2*lnL+2*np;
	PRINT_DIGITS = 3;
	fprintf (stdout,"\n",jobModelNum);
	PRINT_DIGITS = 1;
	fprintf (stdout," | (",vv1,vv2,vv3,vv4,vv5,vv6,") | ");
	fprintf (stdout,Format (np,5,0));
	PRINT_DIGITS = 8;
	fprintf (stdout, " |  ",lnL," | ",Format(LRT,14,3), " |  ", AIC, "  |  ", );
		
	PRINT_DIGITS = 15;
	if (LRT==0)
	{
		pValue = 1;					
	}
	else
	{
		pValue = 1-CChi2(LRT,fullnp-np);
	}
	fprintf (stdout,pValue," |");
	
	resultCache [jobModelNum][0] = vv2;
	resultCache [jobModelNum][1] = vv3;
	resultCache [jobModelNum][2] = vv4;
	resultCache [jobModelNum][3] = vv5;
	resultCache [jobModelNum][4] = vv6;
	resultCache [jobModelNum][5] = lf_MLES[1][0];
	resultCache [jobModelNum][6] = lf_MLES[1][1]+totalBranchCount;
	resultCache [jobModelNum][7] = pValue;
	if (pValue<rejectAt)
	{
		rejectCount = rejectCount+1;
		resultCache [jobModelNum][8] = 0;
	}
	else
	{
		resultCache [jobModelNum][8] = 1;					
	}
	
	if (pValue<rejectAt)
	{
		fprintf (stdout,"(*)");
	}
			
	return fromNode-1;
}

/* ________________________________________________________________________________________________*/

function  setElement (h,v,cc)
{
	if (modelType==0)
	{
		if (cc==0)
		{
			m[h][v]:=a;
			m[v][h]:=a;
		}
		else
		{
			if (cc==1)
			{
				m[h][v]:=b;
				m[v][h]:=b;
			}
			else
			{
				if (cc==2)
				{
					m[h][v]:=c;
					m[v][h]:=c;
				}
				else
				{
					if (cc==3)
					{
						m[h][v]:=d;
						m[v][h]:=d;
					}
					else
					{
						if (cc==4)
						{
							m[h][v]:=e;
							m[v][h]:=e;
						}
						else
						{
							m[h][v]:=f;
							m[v][h]:=f;
						}
					}
				}
			}
		}
	}
	return 1;
}

/* ________________________________________________________________________________________________*/

function ModelSelect (modelType, branchLengths, rejectAt)
{
	svl = VERBOSITY_LEVEL;
	VERBOSITY_LEVEL = -1;
	m={4,4};

	m[0][1]:=a;
	m[1][0]:=a;

	if (modelType > 0)
	{
		global AC;
		global AT;
		global CG;
		global CT;
		global GT;
		
		if (modelType == 1)
		{
			m = {{*,AC*t,t,AT*t}
				 {AC*t,*,CG*t,CT*t}
				 {t,CG*t,*,GT*t}
				 {AT*t,CT*t,GT*t,*}};
		}
		else
		{
			if (modelType == 2)
			{
				incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"defineGamma.mdl";
				ExecuteCommands  ("#include \""+incFileName+"\";");
			}
			else
			{
				if (modelType == 3)
				{
					incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"defineHM.mdl";
					ExecuteCommands  ("#include \""+incFileName+"\";");
				}
			}
			m = {{*,AC*c*t,c*t,AT*c*t}
				 {AC*c*t,*,CG*c*t,CT*c*t}
				 {c*t,CG*c*t,*,GT*c*t}
				 {AT*c*t,CT*c*t,GT*c*t,*}};
		}
	}
	else
	{
		r = setElement (0,2,1);
		r = setElement (0,3,2);
		r = setElement (1,2,3);
		r = setElement (1,3,4);
		r = setElement (2,3,5);
	}

	KEEP_OPTIMAL_ORDER = 1;
	MESSAGE_LOGGING    = 0;

	totalBranchCount = 0;
	modelNum	= 0;
	rejectCount = 0;
	resultCache = {203,9};


	Model currentModel = (m,observedFreqs);
	Tree tr = treeString;
	LikelihoodFunction lf = (filteredData, tr);

	if (MPI_NODE_COUNT>1)
	{
		if (branchLengths == 0)
		{
			SHORT_MPI_RETURN = 1;
		}
		MPINodeState = {MPI_NODE_COUNT-1,8};
		OPTIMIZE_SUMMATION_ORDER = 0;
		MPISend (1,lf);
		MPINodeState[0][0] = 1;
		MPINodeState[0][1] = modelNum;
		if (branchLengths)
		{
			dummy = ReceiveJobs (0);
			SHORT_MPI_RETURN = 1;
		}
	}
	else
	{
		Optimize (lf_MLES,lf);
		vv1 = 0;
		vv2 = 0;
		vv3 = 0;
		vv4 = 0;
		vv5 = 0;
		vv6 = 0;
		dummy = ReceiveJobs (0);
	}

	if (branchLengths)
	{
		totalBranchCount     = TipCount(tr) + BranchCount (tr);
		stashedLengths		 = {totalBranchCount,1};
		
		branchNames = BranchName (givenTree,-1);
		
		pia = observedFreqs[0];
		pic = observedFreqs[1];
		pig = observedFreqs[2];
		pit = observedFreqs[3];
		
		global totalFactor := AC*(2*pia__*pic__)+2*pia__*pig__+(2*pia__*pit__)*AT+
								 (2*pic__*pig__)*CG+(2*pic__*pit__)*CT+(2*pig__*pit__)*GT;
								 
		for (v2 = 0; v2 < totalBranchCount; v2 = v2+1)
		{
			stashedLengths [v2] = BranchLength(tr,v2);
		}
	}

	rateBiasTerms = {{"AC","1","AT","CG","CT","GT"}};

	for (v2=0; v2<=1; v2=v2+1)
	{
		for (v3=0; v3<=v2+1; v3=v3+1)
		{
			if (v3>v2)
			{
				ub4 = v3;
			}
			else
			{
				ub4 = v2;
			}
			for (v4=0; v4<=ub4+1; v4=v4+1)
			{
				if (v4>=ub4)
				{
					ub5 = v4;
				}
				else
				{
					ub5 = ub4;
				}
				for (v5=0; v5<=ub5+1; v5=v5+1)
				{
					if (v5>ub5)
					{
						ub6 = v5;
					}
					else
					{
						ub6 = ub5;
					}
					for (v6=0; v6<=ub6+1; v6=v6+1)
					{
						if (v6==5)
						{
							break;
						}
						
						if (modelType > 0)
						{
							paramCount	  = 0;

							modelDesc = "0"+Format(v2,1,0);
							modelDesc = modelDesc+Format(v3,1,0);
							modelDesc = modelDesc+Format(v4,1,0);
							modelDesc = modelDesc+Format(v5,1,0);
							modelDesc = modelDesc+Format(v6,1,0);
							
							modelConstraintString = "";
							
							AC = 1;
							AT = 1;
							CG = 1;
							CT = 1;
							GT = 1;

							for (customLoopCounter2=1; customLoopCounter2<6; customLoopCounter2=customLoopCounter2+1)
							{
								for (customLoopCounter=0; customLoopCounter<customLoopCounter2; customLoopCounter=customLoopCounter+1)
								{
									if (modelDesc[customLoopCounter2]==modelDesc[customLoopCounter])
									{
										if (rateBiasTerms[customLoopCounter2] == "1")
										{
											modelConstraintString = modelConstraintString + rateBiasTerms[customLoopCounter]+":="+rateBiasTerms[customLoopCounter2]+";";
										}
										else
										{
											modelConstraintString = modelConstraintString + rateBiasTerms[customLoopCounter2]+":="+rateBiasTerms[customLoopCounter]+";";			
										}
										break;
									}
								}
							}	

							if (Abs(modelConstraintString))
							{
								ExecuteCommands (modelConstraintString);
							}
						}
						else
						{
							r = setElement 		(0,2,v2);
							r = setElement 		(0,3,v3);
							r = setElement 		(1,2,v4);
							r = setElement 		(1,3,v5);
							r = setElement 		(2,3,v6);
						}
										
						Model currentModel = (m,observedFreqs);
						if (modelType == 0)
						{
							Tree tr = (dummy,dummy2);
						}
						Tree tr = treeString;
						if (branchLengths)
						{
							for (mpiNode = 0; mpiNode < totalBranchCount; mpiNode = mpiNode+1)
							{
								eCommand = "tr."+branchNames[mpiNode]+".t:="+Format(stashedLengths [mpiNode],20,12)+"/totalFactor";
								ExecuteCommands (eCommand);
							}
						}
						
						LikelihoodFunction lf = (filteredData, tr);
						
						modelNum = modelNum+1;
						if (MPI_NODE_COUNT>1)
						{
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
								mpiNode = ReceiveJobs (1);
							}
							else
							{
								MPISend (mpiNode+1,lf);
								MPINodeState[mpiNode][0] = 1;
								MPINodeState[mpiNode][1] = modelNum;
								MPINodeState[mpiNode][2] = v1;
								MPINodeState[mpiNode][3] = v2;
								MPINodeState[mpiNode][4] = v3;
								MPINodeState[mpiNode][5] = v4;
								MPINodeState[mpiNode][6] = v5;
								MPINodeState[mpiNode][7] = v6;
							}
						}
						else
						{
							Optimize (lf_MLES,lf);
							vv1 = v1;
							vv2 = v2;
							vv3 = v3;
							vv4 = v4;
							vv5 = v5;
							vv6 = v6;
							dummy = ReceiveJobs (0);
						}
					}
				}
			}
		}

	}

	if (MPI_NODE_COUNT>1)
	{
		while (1)
		{
			for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
			{
				if (MPINodeState[nodeCounter][0]==1)
				{
					fromNode = ReceiveJobs (0);
					break;	
				}
			}
			if (nodeCounter == MPI_NODE_COUNT-1)
			{
				break;
			}
		}	
		OPTIMIZE_SUMMATION_ORDER = 1;
	}

	PRINT_DIGITS = 0;

	fprintf (stdout, "\n\n--------------------------\n   (*) => p-Value < ", rejectAt, "\nRejected ", rejectCount, " models.\n");


	if (rejectCount<202)
	{
		fprintf (stdout, "\nPerforming nested tests on the remaining models...\n");

		done = 0;
		while (!done)
		{
			done = 1;
			for (v2=1; v2<203; v2=v2+1)
			{
				if (resultCache[v2][8])
				{
					modelString = "0";
					for (v3 = 0; v3<5; v3=v3+1)
					{
						modelString = modelString + resultCache [v2][v3];
					}
					for (v3 = v2+1; v3<203; v3 = v3+1)
					{
						if (resultCache[v3][8])
						{
							modelString2 = "0";
							for (v4 = 0; v4<5; v4=v4+1)
							{
								modelString2 = modelString2 + resultCache [v3][v4];
							}	
							if (checkEmbedding (modelString, modelString2))
							{
								fprintf (stdout,"H: (", modelString,") A: (", modelString2, "). ");
								done = 0;
								LRT = 2*(resultCache[v3][5]-resultCache[v2][5]);
								npd = resultCache[v3][6]-resultCache[v2][6];
								if (LRT<0)
								{
									pValue = 1;
								}
								else
								{
									pValue = 1-CChi2(LRT,npd);
								}
								fprintf (stdout," P-Value=", Format (pValue,10,3));
								if (pValue<rejectAt)
								{
									fprintf (stdout,". Rejected H.\n");
									resultCache[v2][8] = 0;
									break;
								}
								else
								{
									fprintf (stdout,". Failed to reject H. Discarding A.\n");
									resultCache[v3][8] = 0;
								}
							}
						}
					}
				}
			}
		}

		
		fprintf (stdout,"\n\nRemaining models:\n\n#   |  Model   | # prm |    lnL    |      LRT       |    AIC     |   P-Value        |");   
		fprintf (stdout,"\n----|----------|-------|-----------|----------------|------------|------------------|"); 
		
		modelNum = 0;  
		v5 = 1e10;
		v4 = 0;

		for (v2=1; v2<203; v2=v2+1)
		{
			if (resultCache[v2][8])
			{
				modelNum = 0;
				modelString = "0";
				for (v3 = 0; v3<5; v3=v3+1)
				{
					modelString = modelString + resultCache [v2][v3];
				}
				np  = resultCache[v2][6];
				lnL = resultCache[v2][5];
				LRT = -2*(lnL-stdl);
				if (LRT<0)
				{
					LRT = 0;
				}
				AIC = -2*lnL+2*np;
				modelNum = modelNum + 1;
				PRINT_DIGITS = 3;
				fprintf (stdout,"\n",v2);
				PRINT_DIGITS = 1;
				fprintf (stdout," | (",0,resultCache[v2][0],resultCache[v2][1],resultCache[v2][2],resultCache[v2][3],resultCache[v2][4],") | ");
				fprintf (stdout,Format (np,5,0));
				PRINT_DIGITS = 8;
				fprintf (stdout, " |  ",lnL," | ",Format(LRT,14,3), " |  ", AIC, "  |  ", );
				PRINT_DIGITS = 15;
				if (LRT==0)
				{
					pValue = 1;					
				}
				else
				{
					pValue = 1-CChi2(LRT,fullnp-np);
				}
				if (AIC<v5)
				{
					v5 = AIC;
					v4 = v2;
				}
				fprintf (stdout,pValue," |");
				
			}
		}
		
		PRINT_DIGITS = 0;
		modelString = "0";
		for (v3 = 0; v3<5; v3=v3+1)
		{
			modelString = modelString + Format(resultCache [v4][v3],0,0);
		}
		
		fprintf (stdout, "\n\nAIC based winner: (", modelString, ") with AIC = ", v5, "\n\n");	
	}
	else
	{
		fprintf (stdout, "\nGeneral Reversible Model is the winner!\n");
		modelString = "012345";
	}
	VERBOSITY_LEVEL = svl;
	return modelString;
}

/* ________________________________________________________________________________________________*/

function TreeMatrix2TreeString (doLengths)
{
	treeString = "";
	p = 0;
	k = 0;
	m = treeNodes[0][1];
	n = treeNodes[0][0];
	d = treeString*(Rows(treeNodes)*25);

	while (m)
	{	
		if (m>p)
		{
			if (p)
			{
				d = treeString*",";
			}
			for (j=p;j<m;j=j+1)
			{
				d = treeString*"(";
			}
		}
		else
		{
			if (m<p)
			{
				for (j=m;j<p;j=j+1)
				{
					d = treeString*")";
				}
			}	
			else
			{
				d = treeString*",";
			}	
		}
		if (n<ds.species)
		{
			GetString (nodeName, ds, n);
			d = treeString*nodeName;
		}
		if (doLengths>.5)
		{
			nodeName = ":"+treeNodes[k][2];
			d = treeString*nodeName;
		}
		k=k+1;
		p=m;
		n=treeNodes[k][0];
		m=treeNodes[k][1];
	}

	for (j=m;j<p;j=j+1)
	{
		d = treeString*")";
	}
	
	d=treeString*0;
	return treeString;
}

/* ________________________________________________________________________________________________*/

function InferTreeTopology(optionFlag)
{
	distanceMatrix = {ds.species,ds.species};	
	if (optionFlag == 0)
	{
		incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"chooseDistanceFormula.def";
		ExecuteCommands  ("#include \""+incFileName+"\";");
		InitializeDistances (0);
		for (i = 0; i<ds.species; i=i+1)
		{
			for (j = i+1; j<ds.species; j = j+1)
			{
				distanceMatrix[i][j] = ComputeDistanceFormula (i,j);
			}
		}
	}
	else
	{
		svl = VERBOSITY_LEVEL;
		VERBOSITY_LEVEL = -1;
		for (i = 0; i<ds.species; i=i+1)
		{
			for (j = i+1; j<ds.species; j = j+1)
			{
				DataSetFilter	dFilter = CreateFilter(simData,1,"",speciesIndex == i || speciesIndex == j);
				Tree			dTree	= (1,2);
				LikelihoodFunction dLF  = (dFilter, dTree);
				Optimize (dRes, dLF);
				distanceMatrix[i][j] = BranchLength (dTree,0);
			}
		}
		VERBOSITY_LEVEL = svl;
	}

	MESSAGE_LOGGING 		 	= 1;
	cladesMade 					= 1;
	
	if (ds.species == 2)
	{
		d1 = distanceMatrix[0][1]/2;
		treeNodes = {{0,1,d1__},
					 {1,1,d1__},
					 {2,0,0}};
					 
		cladesInfo = {{2,0}};
	}
	else
	{
		if (ds.species == 3)
		{
			d1 = (distanceMatrix[0][1]+distanceMatrix[0][2]-distanceMatrix[1][2])/2;
			d2 = (distanceMatrix[0][1]-distanceMatrix[0][2]+distanceMatrix[1][2])/2;
			d3 = (distanceMatrix[1][2]+distanceMatrix[0][2]-distanceMatrix[0][1])/2;
			treeNodes = {{0,1,d1__},
						 {1,1,d2__},
						 {2,1,d3__}
						 {3,0,0}};
						 
			cladesInfo = {{3,0}};		
		}
		else
		{	
			njm = (distanceMatrix > methodIndex)>=ds.species;
				
			treeNodes 		= {2*(ds.species+1),3};
			cladesInfo	    = {ds.species-1,2};
			
			for (i=Rows(treeNodes)-1; i>=0; i=i-1)
			{
				treeNodes[i][0] = njm[i][0];
				treeNodes[i][1] = njm[i][1];
				treeNodes[i][2] = njm[i][2];
			}

			for (i=Rows(cladesInfo)-1; i>=0; i=i-1)
			{
				cladesInfo[i][0] = njm[i][3];
				cladesInfo[i][1] = njm[i][4];
			}
			
			njm = 0;
		}
	}
	distanceMatrix = 0;
	return TreeMatrix2TreeString (1);
}

/*--------------------------- MAIN SECTION BEGINS HERE -----------------------------*/

SetDialogPrompt 	("Please specify a nucleotide data file:");
DataSet ds 			= ReadDataFile (PROMPT_FOR_FILE);

fprintf (stdout, "\nRead the following data:", ds,"\n\n");

#include "partitionSequences.ibf";

DataSetFilter 	   filteredData = CreateFilter (ds,1);
HarvestFrequencies (observedFreqs,filteredData,1,1,1);

treeString	= 	InferTreeTopology (0);

fprintf (stdout, "\Using this initial NJ topology: ", treeString);

Topology givenTree = treeString;
if (givenTree <= splitTop)
{
	fprintf (stdout, "\n\n\t ** This topology supports the split.\n");
}
else
{
	fprintf (stdout, "\n\n\t ** This topology DOES NOT support the split.\n");
}



ChoiceList  (mOptions,"Model Options",1,NO_SKIP,
			 "Auto-select","Run a model selection procedure to choose the appropriate model.",
			 "User","User chooses a model.");
	

if (mOptions < 0)
{
	return;	
}

global AC = 1;
global AT = 1;
global CG = 1;
global CT = 1;
global GT = 1;

MGCustomRateBiasTerms = {{"AC*","","AT*","CG*","CT*","GT*"}};	

if (mOptions)
{
	done = 0;
	while (!done)
	{
		fprintf (stdout,"\nPlease enter a 6 character model designation (e.g:010010 defines HKY85):");
		fscanf  (stdin,"String", modelDesc);
		if (Abs(modelDesc)==6)
		{	
			done = 1;
		}
	}			
}
else
{
	modelDesc	= 	ModelSelect (1,1,0.0002);
}

paramCount	  = 0;
_nucBiasTerms = {4,4};
_nucBiasTerms[0][0] = "";


if (modelDesc[0]==modelDesc[1])
{
	MGCustomRateBiasTerms[0] = MGCustomRateBiasTerms[1];
}

_nucBiasTerms[1][0] = MGCustomRateBiasTerms[0];
_nucBiasTerms[0][1] = MGCustomRateBiasTerms[0];
_nucBiasTerms[2][0] = MGCustomRateBiasTerms[1];
_nucBiasTerms[0][2] = MGCustomRateBiasTerms[1];

h = 0;
v = 3;

for (customLoopCounter2=2; customLoopCounter2<6; customLoopCounter2=customLoopCounter2+1)
{
	for (customLoopCounter=0; customLoopCounter<customLoopCounter2; customLoopCounter=customLoopCounter+1)
	{
		if (modelDesc[customLoopCounter]==modelDesc[customLoopCounter2])
		{
			_nucBiasTerms[h][v] = MGCustomRateBiasTerms[customLoopCounter];
			_nucBiasTerms[v][h] = MGCustomRateBiasTerms[customLoopCounter];
			break;
		}
	}
	if (customLoopCounter == customLoopCounter2)
	{
		_nucBiasTerms[h][v] = MGCustomRateBiasTerms[customLoopCounter2];
		_nucBiasTerms[v][h] = MGCustomRateBiasTerms[customLoopCounter2];
	}
	
	v = v+1;
	if (v==4)
	{
		h=h+1;
		v=h+1;
	}
}

c44matrix = 0;
MULTIPLY_BY_FREQS 			= PopulateModelMatrix ("c44matrix", 0);
Model c44model		 		= (c44matrix,observedFreqs,MULTIPLY_BY_FREQS);
Tree   NRTree 				= treeString;

LikelihoodFunction lfNR = (filteredData, NRTree);
Optimize (resNR, lfNR);

fprintf (stdout, "\n\nNo rate variation likelihood fit results:", lfNR, "\n");

ChoiceList  (rvOptions,"Rate Variation Options",1,NO_SKIP,
			 "Test","Test for site-to-site rate variation.",
			 "Skip the test","Skip this test.");
	

if (rvOptions < 0)
{
	return;	
}

execString = "Model simModel = (c44matrix,simFreq,1)";		

if (rvOptions == 0)
{
	c44matrixRV = 0;
	MULTIPLY_BY_FREQS 			= PopulateModelMatrix ("c44matrixRV", 1);
	Model c44modelRV		 	= (c44matrixRV,observedFreqs,MULTIPLY_BY_FREQS);
	Tree   RVTree 				= treeString;
	LikelihoodFunction lfRV = (filteredData, RVTree);
	Optimize (resRV, lfRV);
	
	fprintf (stdout, "\n\nSite-to-site rate variation likelihood fit results:", lfRV, "\n");
	
	AICNR = 2*(resNR[1][1]-resNR[1][0]);
	AICRV = 2*(resRV[1][1]-resRV[1][0]);
	
	fprintf (stdout, "\n\n AIC comparison:\n\t No rate variation: ", AICNR, "\n\t Site-to-site rate variation: ", AICRV, "\n");
	
	if (AICRV < AICNR)
	{
		fprintf (stdout, "\n\nChoosing the site-to-site rate variation model.");
		for (siteCount = 0; siteCount < resRV[1][2]; siteCount = siteCount+1)
		{
			GetString (globalVarName,lfRV,siteCount);
			ExecuteCommands (globalVarName+":="+globalVarName+"__;");
			execString = "Model simModel = (c44matrixRV,simFreq,1)";
		}	
	}
	else
	{
		for (siteCount = 0; siteCount < resNR[1][2]; siteCount = siteCount+1)
		{
			GetString (globalVarName,lfNR,siteCount);
			ExecuteCommands (globalVarName+":="+globalVarName+"__;");
		}	
		fprintf (stdout, "\n\nChoosing the homogeneous rates model.");	
	}
}	
else
{
	for (siteCount = 0; siteCount < resNR[1][2]; siteCount = siteCount+1)
	{
		GetString (globalVarName,lfNR,siteCount);
		ExecuteCommands (globalVarName+":="+globalVarName+"__;");
	}	
}

bootstrapSamples = 0;
supportClade	 = 0;

while (bootstrapSamples<=0)
{
	fprintf (stdout,"\nHow many bootstrap replicates (>=1)?");
	fscanf  (stdin,"Number", bootstrapSamples);
}	

SetDialogPrompt ("Spool boostrap trees to:");
fprintf (PROMPT_FOR_FILE,CLEAR_FILE);
bsTrees = LAST_FILE_PATH;
for (itCount = 1; itCount <= bootstrapSamples; itCount = itCount+1)
{
	DataSetFilter simData = Bootstrap(ds,1);
	HarvestFrequencies (simFreq,simData,1,1,0);
	ExecuteCommands (execString);
	SimTreeString	= 	InferTreeTopology (1);
	fprintf (bsTrees,SimTreeString,"\n");
	fprintf (stdout, "Iteration ", Format(itCount,5,0), "/", bootstrapSamples);
	Topology SimTree = SimTreeString;
	if (SimTree<=splitTop)
	{
		fprintf (stdout, " + ");	
		supportClade = supportClade + 1;
	}
	else
	{
		fprintf (stdout, " - ");		
	}
	fprintf (stdout, "\t Running support value:", Format (supportClade/itCount,10,6), "\n");
}
