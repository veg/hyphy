RequireVersion ("0.9920061127");

/* 
	This file takes a nucleotide data set and a tree (either from the data file or from a separate 	file) and computes maximum likelihood estimates for every possible 4x4 reversible model on that data and tree.

		We use the string (v1,v2,v3,v4,v5,v6), where and v1..6 = 0..5
   to encode a 4x4 symmetric transition matrix with entries
   [*  v1  v2  v3]
   [-  *   v4  v5]
   [-  -   *   v6]
   [-  -   -   * ]
   
   		For instance: (010010) encodes HKY85.
   		
   		For each model the following information is reported:
   		- Model string. (e.g. (012345) for the GRM)
   		- Number of model parameters
   		- Max ln-likelihood for the model
   		- Likelihood ratio statistic (as a sub-model of the GRM)
   		- AIC
   		- P-Value for the Likelihood Ratio Test.
   		
   
   Sergei L. Kosakovsky Pond, Summer 2002.
   
*/ 

VERBOSITY_LEVEL = -1;

function checkEmbedding (_m1, _m2)
{
	for (r=0; r<6; r=r+1)
	{
		if (_m2[r]<_m1[r])
		{
			/*fprintf (stdout,_m1," ", _m2, " Reject 1 at position ",r,"\n");*/
			return 0;
		}
		if (_m2[r]>_m1[r])
		{
			for (r2 = 0; r2 < 6; r2 = r2+1)
			{
				if ((_m2[r2]==_m2[r])&&(_m1[r2]!=_m1[r]))
				{
					/*fprintf (stdout,_m1," ", _m2, " Reject 2 at positions ",r,r2,"\n");*/
					return 0;
				}
			}
		}
	}
	return 1;
}

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
		fullnp = lf_MLES[1][1]+totalBranchCount+3;
		fprintf(stdout,"\n(012345) Full Model ln-lik =  ",stdl,". Parameter Count=",Format(fullnp,0,0)," AIC = ", 2*(fullnp-stdl),"\n\n");


		resultCache [0][0] = 1;
		resultCache [0][1] = 2;
		resultCache [0][2] = 3;
		resultCache [0][3] = 4;
		resultCache [0][4] = 5;
		resultCache [0][5] = lf_MLES[1][0];
		resultCache [0][6] = lf_MLES[1][1]+totalBranchCount+3;
		resultCache [0][7] = 0;
		resultCache [0][8] = 0;
		if (modelType)
		{
			resultCache [0][9] = AC;
			resultCache [0][10] = AT;
			resultCache [0][11] = CG;
			resultCache [0][12] = CT;
			resultCache [0][13] = GT;
		}
		
		if (doFitSave)
		{
			slfo = LIKELIHOOD_FUNCTION_OUTPUT;
			LIKELIHOOD_FUNCTION_OUTPUT = 6;
			NEW_FILE_NAME = BASE_PATH+".012345";
			fprintf (NEW_FILE_NAME,CLEAR_FILE,lf);
			LIKELIHOOD_FUNCTION_OUTPUT = slfo;
		}

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
					fprintf (stdout, " |  ",lnL," | ",Format(LRT,14,3), " |  ", AIC, "  |  ");
					
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
			resultCache [jobModelNum][6] = lf_MLES[1][1]+totalBranchCount+3;
			
			if (modelType)
			{
				resultCache [jobModelNum][9] = AC;
				resultCache [jobModelNum][10] = AT;
				resultCache [jobModelNum][11] = CG;
				resultCache [jobModelNum][12] = CT;
				resultCache [jobModelNum][13] = GT;
			}
			/*slfo = LIKELIHOOD_FUNCTION_OUTPUT;
			LIKELIHOOD_FUNCTION_OUTPUT = 6;
			NEW_FILE_NAME = BASE_PATH+".0"+Format(resultCache[h][0],1,0);
			NEW_FILE_NAME = NEW_FILE_NAME+Format(resultCache[h][1],1,0);
			NEW_FILE_NAME = NEW_FILE_NAME+Format(resultCache[h][2],1,0);
			NEW_FILE_NAME = NEW_FILE_NAME+Format(resultCache[h][3],1,0);
			NEW_FILE_NAME = NEW_FILE_NAME+Format(resultCache[h][4],1,0);
			fprintf (NEW_FILE_NAME,CLEAR_FILE,lf);
			LIKELIHOOD_FUNCTION_OUTPUT = slfo;*/
			
			return fromNode - 1;
		}
	}

	np = lf_MLES[1][1]+totalBranchCount+3;
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
	fprintf (stdout, " |  ",lnL," | ",Format(LRT,14,3), " |  ", AIC, "  |  ");
		
	if (doFitSave)
	{
		slfo = LIKELIHOOD_FUNCTION_OUTPUT;
		LIKELIHOOD_FUNCTION_OUTPUT = 6;
		NEW_FILE_NAME = BASE_PATH+".0" + Format(vv2,1,0);
		NEW_FILE_NAME = NEW_FILE_NAME  + Format(vv3,1,0);
		NEW_FILE_NAME = NEW_FILE_NAME  + Format(vv4,1,0);
		NEW_FILE_NAME = NEW_FILE_NAME  + Format(vv5,1,0);
		NEW_FILE_NAME = NEW_FILE_NAME  + Format(vv6,1,0);
		
		fprintf (NEW_FILE_NAME,CLEAR_FILE,lf);
		LIKELIHOOD_FUNCTION_OUTPUT = slfo;
	}

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
	resultCache [jobModelNum][6] = lf_MLES[1][1]+totalBranchCount+3;
	resultCache [jobModelNum][7] = pValue;
	if (modelType)
	{
		if (MPI_NODE_COUNT>1)
		{
			resultCache [jobModelNum][9]  = lf_MLE_VALUES["AC"];
			resultCache [jobModelNum][10] = lf_MLE_VALUES["AT"];
			resultCache [jobModelNum][11] = lf_MLE_VALUES["CG"];
			resultCache [jobModelNum][12] = lf_MLE_VALUES["CT"];
			resultCache [jobModelNum][13] = lf_MLE_VALUES["GT"];
		}
		else
		{
			resultCache [jobModelNum][9] = AC;
			resultCache [jobModelNum][10] = AT;
			resultCache [jobModelNum][11] = CG;
			resultCache [jobModelNum][12] = CT;
			resultCache [jobModelNum][13] = GT;		
		}
	}
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


function printModelMatrix (modelString)
{
	
	if (modelType)
	{
		mstrConv = "1";
		for (v2 = 1; v2 < 6; v2 = v2+1)
		{
			if (modelString[v2]=="0")
			{
				mstrConv = mstrConv+"1";
			}
			else
			{
				if (modelString[v2]=="1")
				{
					mstrConv = mstrConv+"B";
				}
				else
				{
					if (modelString[v2]=="2")
					{
						mstrConv = mstrConv+"C";
					}
					else
					{
						if (modelString[v2]=="3")
						{
							mstrConv = mstrConv+"D";
						}
						else
						{
							if (modelString[v2]=="4")
							{
								mstrConv = mstrConv+"E";
							}
							else
							{
								mstrConv = mstrConv+"F";
							}
						}
					}
				}
			}
		}
		sep = "+---+-----+-----+-----+-----+\n";
		fprintf (stdout, sep,
						 "|   |  A  |  C  |  G  |  T  |\n",
						 sep,
						 "| A |  *  | ", mstrConv[0], "*t | ", mstrConv[1], "*t | ", mstrConv[2], "*t |\n",
						 sep,
						 "| C | ", mstrConv[0], "*t |  *  | ", mstrConv[3], "*t | ", mstrConv[4], "*t |\n",
						 sep,
						 "| G | ", mstrConv[1], "*t | " , mstrConv[3], "*t |  *  | ", mstrConv[5], "*t |\n",
						 sep,
						 "| T | ", mstrConv[2], "*t | " , mstrConv[4], "*t | ", mstrConv[5], "*t |  *  |\n",
						 sep);
						 
	
	}
	else
	{
		mstrConv = "a";
		for (v2 = 1; v2 < 6; v2 = v2+1)
		{
			if (modelString[v2]=="0")
			{
				mstrConv = mstrConv+"a";
			}
			else
			{
				if (modelString[v2]=="1")
				{
					mstrConv = mstrConv+"b";
				}
				else
				{
					if (modelString[v2]=="2")
					{
						mstrConv = mstrConv+"c";
					}
					else
					{
						if (modelString[v2]=="3")
						{
							mstrConv = mstrConv+"d";
						}
						else
						{
							if (modelString[v2]=="4")
							{
								mstrConv = mstrConv+"e";
							}
							else
							{
								mstrConv = mstrConv+"f";
							}
						}
					}
				}
			}
		}
		sep = "+---+-----+-----+-----+-----+\n";
		fprintf (stdout, sep,
						 "|   |  A  |  C  |  G  |  T  |\n",
						 sep,
						 "| A |  *  |  ", mstrConv[0], "  |  ", mstrConv[1], "  |  ", mstrConv[2], "  |\n",
						 sep,
						 "| C |  ", mstrConv[0], "  |  *  |  ", mstrConv[3], "  |  ", mstrConv[4], "  |\n",
						 sep,
						 "| G |  ", mstrConv[1], "  |  " , mstrConv[3], "  |  *  |  ", mstrConv[5], "  |\n",
						 sep,
						 "| T |  ", mstrConv[2], "  |  " , mstrConv[4], "  |  ", mstrConv[5], "  |  *  |\n",
						 sep);
						 
	}
	return 1;
}

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

#include "TemplateModels/modelParameters.mdl";

m={4,4};

m[0][1]:=a;
m[1][0]:=a;

branchLengths = 0;

if (modelType > 0)
{
	ChoiceList (branchLengths,"Estimate Branch Lengths",1,SKIP_NONE,
				"Every Time","Branch lengths are reestimated for every model.",
				"Once","Branch lenghts obtained from the general reversible model are reused for subsequent models."
		       );

	if (branchLengths<0)
	{
		return;
	}

	global AC = 1;
	global AT = 1;
	global CG = 1;
	global CT = 1;
	global GT = 1;
	
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
			/*m[0][1]:=c*t;
			m[1][0]:=c*t;*/
			#include "TemplateModels/defineGamma.mdl";
		}
		else
		{
			if (modelType == 3)
			{
				/*m[0][1]:=c*t;
				m[1][0]:=c*t;*/
				#include "TemplateModels/defineHM.mdl";
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

SetDialogPrompt ("Please specify a nucleotide data file:");

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter nucData = CreateFilter (ds,1);

HarvestFrequencies (pi,nucData,1,1,1);

_DO_TREE_REBALANCE_ = 1;
#include "queryTree.bf";

rejectAt = 0;

while (rejectAt<=0 || rejectAt>=1 )
{
	fprintf (stdout, "\nModel rejection level (e.g. 0.05):");
	fscanf  (stdin,"Number", rejectAt);
}

ChoiceList (doFitSave, "Save each of the 203 fits?", 1, SKIP_NONE, "No", "Do not write out the files", "Yes", "Save each of the 203 files to a separate file");
if (doFitSave<0)
{
	return 0;
}

SetDialogPrompt ("Save summary and results to:");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE);
BASE_PATH = LAST_FILE_PATH;


KEEP_OPTIMAL_ORDER = 1;
MESSAGE_LOGGING    = 0;

totalBranchCount = 0;
modelNum	= 0;
rejectCount = 0;
resultCache = {203,14};


Model currentModel = (m,pi);
Tree tr = treeString;
LikelihoodFunction lf = (nucData, tr, pi);

meanDiv = -1;

if (MPI_NODE_COUNT>1)
{
	MPINodeState = {MPI_NODE_COUNT-1,8};
	SHORT_MPI_RETURN		 = 1;
	OPTIMIZE_SUMMATION_ORDER = 0;
	if (branchLengths > 0)
	{
		AUTO_PARALLELIZE_OPTIMIZE = 1;
		Optimize 	(lf_MLES,lf);
		AUTO_PARALLELIZE_OPTIMIZE = 0;
		smnc = MPI_NODE_COUNT;
		MPI_NODE_COUNT = 1;
		ReceiveJobs (0);
		MPI_NODE_COUNT = smnc;
		
	}
	else
	{
		MPISend (1,lf);
		MPINodeState[0][0] = 1;
		MPINodeState[0][1] = modelNum;
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
	
	pia = pi[0];
	pic = pi[1];
	pig = pi[2];
	pit = pi[3];
	
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
									
					Model currentModel = (m,pi);
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
					
					LikelihoodFunction lf = (nucData, tr);
					
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

	fprintf (BASE_PATH,CLEAR_FILE,"\n\nRemaining models:\n\n#   |  Model   | # prm |    lnL    |      LRT       |    AIC     |   P-Value        |");
	fprintf (BASE_PATH,"\n----|----------|-------|-----------|----------------|------------|------------------|");
	
	fprintf (stdout,"\n\nRemaining models:\n\n#   |  Model   | # prm |    lnL    |      LRT       |    AIC     |   P-Value        |");   
	fprintf (stdout,"\n----|----------|-------|-----------|----------------|------------|------------------|"); 
	
	modelNum = 0;  
	v5 = 1e10;
	v4 = 0;
	
	for (v2=1; v2<203; v2=v2+1)
	{
		if (resultCache[v2][8])
		{
			np  = resultCache[v2][6];
			lnL = resultCache[v2][5];
			AIC = -2*lnL+2*np;
			modelNum = 0;
			modelString = "0";
			for (v3 = 0; v3<5; v3=v3+1)
			{
				modelString = modelString + resultCache [v2][v3];
			}
			LRT = -2*(lnL-stdl);
			if (LRT<0)
			{
				LRT = 0;
			}
			modelNum = modelNum + 1;
			PRINT_DIGITS = 3;
			fprintf (stdout,"\n",v2);
			fprintf (BASE_PATH,"\n",v2);
			PRINT_DIGITS = 1;
			fprintf (stdout," | (",0,resultCache[v2][0],resultCache[v2][1],resultCache[v2][2],resultCache[v2][3],resultCache[v2][4],") | ");
			fprintf (stdout,Format (np,5,0));
			fprintf (BASE_PATH," | (",0,resultCache[v2][0],resultCache[v2][1],resultCache[v2][2],resultCache[v2][3],resultCache[v2][4],") | ");
			fprintf (BASE_PATH,Format (np,5,0));
			PRINT_DIGITS = 8;
			fprintf (stdout, " |  ",lnL," | ",Format(LRT,14,3), " |  ", AIC, "  |  ", );
			fprintf (BASE_PATH, " |  ",lnL," | ",Format(LRT,14,3), " |  ", AIC, "  |  ", );
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
			fprintf (BASE_PATH,pValue," |");	
		}
		
	}
	

	PRINT_DIGITS = 0;
	modelString = "0";
	for (v3 = 0; v3<5; v3=v3+1)
	{
		modelString = modelString + Format(resultCache [v4][v3],0,0);
	}
	
	fprintf (stdout, "\n\nAIC based winner: (", modelString, ") with AIC = ", v5, "\n\n");
	fprintf (BASE_PATH, "\n\nAIC based winner: (", modelString, ") with AIC = ", v5, "\n\n");
	
	bestModelString = modelString;
	
	printModelMatrix (modelString);
	
	modelString2 = "";
	if (modelString == "000000")
	{
		modelString2 = "F81";
	}
	if (modelString == "010010")
	{
		modelString2 = "HKY85";
	}
	if (modelString == "010020")
	{
		modelString2 = "TrN";
	}
	if (Abs(modelString2))
	{
		fprintf (stdout, "\nThis model is better known as:", modelString2, "\n");
		fprintf (BASE_PATH, "\nThis model is better known as:", modelString2, "\n");
	}

}
else
{
	fprintf (stdout, "\nGeneral Reversible Model is the winner!\n");
	fprintf (BASE_PATH, "\nGeneral Reversible Model is the winner!\n");
	bestModelString = "012345";
}

if (modelType)
{
	modelAICs = {203,2};
	modelAICs[0][0] = 2*(-resultCache[0][5]+resultCache[0][6]+Log(nucData.sites));
	modelAICs[0][1] = 0;

	for (v2=1; v2<203; v2=v2+1)
	{
		modelAICs[v2][0] = 2*(resultCache[v2][6]-resultCache[v2][5]+Log(nucData.sites));
		modelAICs[v2][1] = v2;
	}

	modelAICs = modelAICs%0;
	v3 = 1;

	/* compute Akaike weights */

	for (v2=1; v2<203; v2=v2+1)
	{
		v4 = Exp(0.5*(modelAICs[0][0]-modelAICs[v2][0]));
		modelAICs[v2][0] = v4;
		v3 = v3+v4;
	}

	modelAICs[0][0] = 1;

	for (v2=0; v2<203; v2=v2+1)
	{
		modelAICs[v2][0] = modelAICs[v2][0]/v3;
	}




	/* now compute model averaged stuff */

	modelAveragedRates   = {5,1};
	modelAveragedEqual	 = {6,6};

	for (v2=0; v2<203; v2=v2+1)
	{
		v3  = modelAICs[v2][1];
		AIC = modelAICs[v2][0];
		for (v4 = 0; v4 < 5; v4 = v4 + 1)
		{
			aRate = resultCache[v3][9+v4];
			modelAveragedRates[v4] = modelAveragedRates[v4] + aRate*AIC;
			if (aRate == 1)
			{
				modelAveragedEqual[1][v4+(v4>0)] = modelAveragedEqual[1][v4+(v4>0)] + AIC;
				modelAveragedEqual[v4+(v4>0)][1] = modelAveragedEqual[v4+(v4>0)][1] + AIC;
			}
		} 
		for (v4=1; v4<5; v4=v4+1)
		{
			for (v5=0; v5<v4; v5=v5+1)
			{
				if (resultCache[v3][9+v4] == resultCache[v3][9+v5])
				{
					modelAveragedEqual[v4+(v4>0)][v5+(v5>0)] = modelAveragedEqual[v4+(v4>0)][v5+(v5>0)] + AIC;
					modelAveragedEqual[v5+(v5>0)][v4+(v4>0)] = modelAveragedEqual[v5+(v5>0)][v4+(v4>0)] + AIC;
				}
			}	
		}
	}

	fprintf (stdout, "\nModel averaged rates relative to AG (REV estimates):\n\tAC = ", Format(modelAveragedRates[0],8,4),"\t(",Format(resultCache[0][9],8,4),")",
										    "\n\tAT = ", Format(modelAveragedRates[1],8,4),"\t(",Format(resultCache[0][10],8,4),")",
										    "\n\tCG = ", Format(modelAveragedRates[2],8,4),"\t(",Format(resultCache[0][11],8,4),")",
										    "\n\tCT = ", Format(modelAveragedRates[3],8,4),"\t(",Format(resultCache[0][12],8,4),")",
										    "\n\tGT = ", Format(modelAveragedRates[4],8,4),"\t(",Format(resultCache[0][13],8,4),")\n");
										    
	fprintf (BASE_PATH, "\nModel averaged rates relative to AG (REV estimates):\n\tAC = ", Format(modelAveragedRates[0],8,4),"\t(",Format(resultCache[0][9],8,4),")",
										    "\n\tAT = ", Format(modelAveragedRates[1],8,4),"\t(",Format(resultCache[0][10],8,4),")",
										    "\n\tCG = ", Format(modelAveragedRates[2],8,4),"\t(",Format(resultCache[0][11],8,4),")",
										    "\n\tCT = ", Format(modelAveragedRates[3],8,4),"\t(",Format(resultCache[0][12],8,4),")",
										    "\n\tGT = ", Format(modelAveragedRates[4],8,4),"\t(",Format(resultCache[0][13],8,4),")\n");

	rateAssignments = {6,1};
	for (v2 = 0; v2 < 6; v2 = v2+1)
	{
		v4 = rateAssignments[v2];
		
		for (v3 = v2 + 1; v3 < 6; v3 = v3+1)
		{
			v5 = rateAssignments[v3];
			if (v4 == v5)
			{
				if (modelAveragedEqual[v2][v3]<0.0005)
				{
					rateAssignments[v3] = rateAssignments[v3]+1;
				}
			}
		}
	}
	
	modelString = "";
	for (v2=0;v2<6;v2=v2+1)
	{
		modelString = modelString+Format(rateAssignments[v2],1,0);
	}

	fprintf (stdout, "\n\nModel averaged selection: (", modelString, ")\n\n");
	fprintf (BASE_PATH, "\n\nModel averaged selection: (", modelString, ")\n\n");
	printModelMatrix (modelString);
	
}

