/* 

THIS HyPhy BATCH LANGUAGE FILE IS PROVIDED AS AN EXAMPLE/TOOL FOR 
THE 2nd EDITION OF 'The Phylogenetic Handbook'

Written by 
Sergei L Kosakovsky Pond (spond@ucsd.edu)
Art FY Poon				 (apoon@ucsd.edu)
Simon DW Frost			 (sdfrost@ucsd.edu)

*/


RequireVersion ("0.9920060830");

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function BuildCodonFrequencies (obsF)
{
	PIStop = 1.0;
	result = {ModelMatrixDimension,1};
	hshift = 0;

	for (h=0; h<64; h=h+1)
	{
		first = h$16;
		second = h%16$4;
		third = h%4;
		if (_Genetic_Code[h]==10) 
		{
			hshift = hshift+1;
			PIStop = PIStop-obsF[first][0]*obsF[second][1]*obsF[third][2];
			continue; 
		}
		result[h-hshift][0]=obsF[first][0]*obsF[second][1]*obsF[third][2];
	}
	return result*(1.0/PIStop);
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function TreeMatrix2TreeString (levelIndex)
{
	treeString = "";
	p 		   = 0;
	k 		   = 0;
	m 		   = treeNodes[0][levelIndex+1];
	n 		   = treeNodes[0][levelIndex];

	while (m)
	{	
		if (m>p)
		{
			if (p)
			{
				treeString = treeString+",";
			}
			for (j=p;j<m;j=j+1)
			{
				treeString = treeString+"(";
			}
		}
		else
		{
			if (m<p)
			{
				for (j=m;j<p;j=j+1)
				{
					treeString = treeString+")";
				}
			}	
			else
			{
				treeString = treeString+",";
			}	
		}
		if (n<ds.species)
		{
			GetString (nodeName, ds, n);
			treeString = treeString+nodeName+":0.05";
		}
		k=k+1;
		p=m;
		n=treeNodes[k][levelIndex];
		m=treeNodes[k][levelIndex+1];
	}

	for (j=m;j<p;j=j+1)
	{
		treeString = treeString+")";
	}
	
	return treeString;
}


/*---------------------------------------------------------------------------------------------------------------------------------------------*/


fprintf (stdout, "\n\nThis example analysis fits the Muse-Gaut 94 model with the \ntransition/transversion correction to a small codon dataset",
				 "\nusing all (2N-5)!! possible topologies, to investigate the\neffect of incorrect topology on dN/dS rate estimation\n\n");

ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"DescriptiveStatistics.bf");
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def");

ModelMatrixDimension = 64;
for (k=0; k<64; k=k+1)
{
	if (_Genetic_Code[k] == 10)
	{
		ModelMatrixDimension = ModelMatrixDimension -1;
	}
}

ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"2RatesAnalyses"+DIRECTORY_SEPARATOR+"MG94xREV.mdl");

ChoiceList (freqs, "Alignment", 1, SKIP_NONE, "Flu H5N1 HA", 		"Use the example Influenza H5N1 heamagglutinin alignment (5 sequences)",
											  "Drospophila adh", 	"Use the example Drosophila ADH alignment (6 sequences).",
											  "Custom", 			"Load your own 4-7 sequence alignment.");
													 
if (freqs < 0)
{
	return 0;
}
													 
if (freqs == 2)
{
	SetDialogPrompt     ("Choose a nucleotide alignment");
	DataSet ds        = ReadDataFile (PROMPT_FOR_FILE);
}
else
{
	if (freqs == 1)
	{
		GetURL			(dataString, "http://www.hyphy.org/phylohandbook/data/Drosophila_adh.nex");	
	}
	else
	{
		GetURL			(dataString, "http://www.hyphy.org/phylohandbook/data/H5N1_HA_5.nex");
	}	
	DataSet ds		 = ReadFromString (dataString);
}

if (ds.species<4 || ds.species>7)
{
	fprintf (stdout, "\nERROR: the alignment must contain 4-7 sequences\n");
	return 0;
}
			  			  			  			  
HarvestFrequencies  (baseFreqs,ds,3,1,1);
DataSetFilter	  filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
ecf		   		  = BuildCodonFrequencies (baseFreqs);
PopulateModelMatrix ("MGREVMX", baseFreqs, 0);
Model MGREV 	  = (MGREVMX, ecf, 0);

AC					= 1;
AT					:= AC;
CG					:= AC;
CT					:= 1;
GT					:= AC;


treeNodes 		= {2*(ds.species+1),2*(ds.species-2)};
cladesInfo 		= {ds.species,2*(ds.species-2)};
branchIndex		= {ds.species-3,1};
currentLevel 	= 0;
done 			= 0;
treeCounter 	= 0;

i 				= 2*ds.species-5;
totalTreeCount 	= 1;
while (i>1)
{
	totalTreeCount = totalTreeCount*i;
	i = i-2;
}

dNdSArray			=   {totalTreeCount,1};
treeNodes[0][0]		=	0;
treeNodes[0][1]		=	1;
treeNodes[1][0]		=	1;
treeNodes[1][1]		=	1;
treeNodes[2][0]		=	2;
treeNodes[2][1]		=	1;
treeNodes[3][0]		=	ds.species;
treeNodes[3][1]		=	0;
cladesInfo[0][0]	=	0;
cladesInfo[0][1]	=	4;
scores				=   {{0,0}{-1e25,0}};

while (!done)
{
	if (branchIndex[currentLevel]<2*currentLevel+3)
	{
		i = 0;
		shift = 0;
		j = 2*currentLevel;
		k = j+2;
		m = j+1;
		while (treeNodes[i][m])
		{
			/*copy tree from prev level to this level */
			if (i==branchIndex[currentLevel])
			/*insert new branch*/
			{
				shift = 2;
				if (treeNodes[i][j]<ds.species)
				/* simple branch */
				{
					treeNodes[i][k]=treeNodes[i][j];
					treeNodes[i][k+1]=treeNodes[i][m]+1;
					treeNodes[i+1][k]=currentLevel+3;
					treeNodes[i+1][k+1]=treeNodes[i][m]+1;
					treeNodes[i+2][k]=currentLevel+ds.species+1;
					treeNodes[i+2][k+1]=treeNodes[i][m];
					cladesInfo[currentLevel+1][k] = i;
					cladesInfo[currentLevel+1][k+1] = 3;					
				}
				else
				{
					/* update node depths for the entire clade now*/
					l = treeNodes[i][j]-ds.species;
					s = cladesInfo[l][j];
					for (p=s+cladesInfo[l][m]-1; p>=s; p=p-1)
					{
						treeNodes[i][k]=treeNodes[i][j];
						treeNodes[i][k+1]=treeNodes[i][m]+1;						
						i=i-1;
					}
					i=i+cladesInfo[l][m];
					/* new clade record */
					cladesInfo[currentLevel+1][k] = cladesInfo[l][j];
					cladesInfo[currentLevel+1][k+1] = cladesInfo[l][m]+2;
					/* now we need to insert two more nodes */
					treeNodes[i+1][k]=currentLevel+3;
					treeNodes[i+1][k+1]=treeNodes[i][m]+1;
					treeNodes[i+2][k]=currentLevel+ds.species+1;
					treeNodes[i+2][k+1]=treeNodes[i][m];
				}
				for (p=0; p<=currentLevel; p=p+1)
				{
					if (cladesInfo[p][j]>i)
					{
						cladesInfo[p][k] = cladesInfo[p][j]+2;
					}
					else
					{
						cladesInfo[p][k] = cladesInfo[p][j];
					}
					
					if ((cladesInfo[p][j]<=i)&&((cladesInfo[p][j]+cladesInfo[p][m])>i+1))
					{
						cladesInfo[p][k+1] = cladesInfo[p][m]+2;
					}
					else
					{
						cladesInfo[p][k+1] = cladesInfo[p][m];
					}
				}				
			}
			else
			{
				treeNodes[i+shift][k]=treeNodes[i][j];
				treeNodes[i+shift][k+1]=treeNodes[i][m];
			}
			i = i+1;
		}
		treeNodes[i+2][k]=treeNodes[i][j];
		treeNodes[i+2][k+1]=treeNodes[i][j+1];
		if (currentLevel<ds.species-4)
		{
			currentLevel = currentLevel+1;
		}
		else
		{
			R 							= 	1.0;
			AC							= 	1.0;
			thisTree 					= 	TreeMatrix2TreeString (2*(currentLevel+1));
			branchIndex[currentLevel]	=	branchIndex[currentLevel]+1;
			Tree    				  T = thisTree;
			LikelihoodFunction 		 lf = (filteredData, T);
			Optimize 					  (res, lf);
			dNdSArray	[treeCounter]	= R;
			treeCounter 		  		= treeCounter+1;
			fprintf (stdout, "\nTree#",Format(treeCounter,0,0),"/",Format(totalTreeCount,0,0),": ", Format(T,0,1), "\n");
			fprintf (stdout, "Log(L) = ", Format (res[1][0], 10,3), ". Global dN/dS = ", Format (R,10,3), "\n");
			if (res[1][0]<scores[0][0])
			{
				scores[0][0] = res[1][0];
				scores[0][1] = treeCounter-1;
			}
			if (res[1][0]>scores[1][0])
			{
				scores[1][0] = res[1][0];
				scores[1][1] = treeCounter-1;
			}
		}
	}
	else
	{
		branchIndex[currentLevel]=0;
		if (currentLevel==0)
		{
			done = 1;
		}
		else
		{
			currentLevel = currentLevel-1;
			branchIndex[currentLevel]=branchIndex[currentLevel]+1;
		}
	}
}

dsc = GatherDescriptiveStats(dNdSArray);

fprintf (stdout, "\ndN/dS report:\n",
				 "\n\tMean       = ", dsc ["Mean"], 
				 "\n\tVariance   = ", dsc ["Variance"],
				 "\n\tMin        = ", dsc ["Min"],
				 "\n\tMax        = ", dsc ["Max"],
				 "\n\tBest tree  = ", dNdSArray[scores[1][1]],
				 "\n\tWorst tree = ", dNdSArray[scores[0][1]], "\n");