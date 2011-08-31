VERBOSITY_LEVEL = -1;
ALLOW_SEQUENCE_MISMATCH = 1;

fprintf(stdout,"\n ---- RUNNING KH ALTERNATIVE TOPOLOGY TEST ---- \n");

ChoiceList (dataType,"Data type",1,SKIP_NONE,"Nucleotide/Protein","Nucleotide or amino-acid (protein).",
				     "Codon","Codon (several available genetic codes).");

if (dataType<0) 
{
	return;
}

if (dataType)
{
	NICETY_LEVEL = 3;
	SetDialogPrompt ("Please choose a codon data file:");
	#include "TemplateModels/chooseGeneticCode.def";
}
else
{
	SetDialogPrompt ("Please choose a nucleotide or amino-acid data file:");
}

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
fprintf (stdout,"The following data were read:\n",ds,"\n");

ChoiceList (compType,"Comparison Kind",1,SKIP_NONE,"Full Alignment","Compare the trees applied to the full alignment.",
				     "Recombination Test","Compare the trees on 2 parts of the alignment, reversing their roles before and after the breakpoint.");


if (compType<0) 
{
	return;
}

if (dataType)
{
	DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
}
else
{
	DataSetFilter filteredData = CreateFilter (ds,1);
}


SelectTemplateModel(filteredData);

firstTreeString  = "";
secondTreeString = "";

if (Rows(NEXUS_FILE_TREE_MATRIX)>=2)
{
	ChoiceList (tc,"Choose first tree",1,SKIP_NONE,NEXUS_FILE_TREE_MATRIX);
	if (tc>=0)
	{
		firstTreeString = NEXUS_FILE_TREE_MATRIX[tc][1];
		skipMx = {{tc__}};
		ChoiceList (tc,"Choose second tree",1,skipMx,NEXUS_FILE_TREE_MATRIX);
		if (tc>=0)
		{
			secondTreeString = NEXUS_FILE_TREE_MATRIX[tc][1];		
		}
	}
	else
	{
		ChoiceList (tc,"Choose second tree",1,SKIP_NONE,NEXUS_FILE_TREE_MATRIX);	
		if (tc>=0)
		{
			secondTreeString = NEXUS_FILE_TREE_MATRIX[tc][1];		
		}
	}
}
else
{
	SetDialogPrompt ("Load a tree for the first partition");
	fscanf (PROMPT_FOR_FILE, "String", firstTreeString);
	SetDialogPrompt ("Load a tree for the second partition");
	fscanf (PROMPT_FOR_FILE, "String", secondTreeString);
}

if (Abs(firstTreeString) == 0)
{
	SetDialogPrompt ("Locate the file with the first tree:");
	fscanf 			(PROMPT_FOR_FILE,"Tree",firstTree);
}
else
{
	Tree firstTree = firstTreeString;
}

if (Abs(secondTreeString) == 0)
{
	SetDialogPrompt ("Locate the file with the second tree:");
	fscanf 			(PROMPT_FOR_FILE,"Tree",secondTree);
}
else
{
	Tree secondTree = secondTreeString;
}

if (firstTree == secondTree)
{
	fprintf (stdout, "\nERROR: The two trees cannot be equal\n");
	return 0;
}

firstTreeString  = ""+firstTree;
secondTreeString = ""+secondTree;

fprintf (stdout, "\nFirst tree:\n",firstTreeString, "\n\nSecond Tree:\n",secondTreeString,"\n");

itCount = 0;
while (itCount<1)
{
	fprintf (stdout, "\nHow many bootstrap replicates should be done (>0):");
	fscanf (stdin, "Number", itCount);
}

labels = {{"Resampled LRT"}};

if (compType)
{
	breakPoint = 0;
	if (dataType)
	{
		pPrompt = "codons";
	}
	else
	{
		pPrompt = "bp";	
	}
	while (breakPoint<1 || breakPoint > filteredData.sites-2)
	{
		fprintf (stdout, "\nEnter the location of the breakpoint (0-based and measured in ", pPrompt, "). Valid locations are [1-", filteredData.sites-3, "]:");
		fscanf (stdin, "Number", breakPoint);
	}
	fprintf (stdout, "\nUsing the breakpoint at ", breakPoint, "-th ", pPrompt, "\n");
	if (dataType)
	{
		DataSetFilter filteredData1 = CreateFilter (ds,3,siteIndex<3*breakPoint,"",GeneticCodeExclusions);
		DataSetFilter filteredData2 = CreateFilter (ds,3,siteIndex>=3*breakPoint,"",GeneticCodeExclusions);
	}
	else
	{
		DataSetFilter filteredData1 = CreateFilter (ds,1,siteIndex<breakPoint);
		DataSetFilter filteredData2 = CreateFilter (ds,1,siteIndex>=breakPoint);	
	}
	LikelihoodFunction lf1_1 = (filteredData1,firstTree);
	Optimize (res1_1, lf1_1);
	ConstructCategoryMatrix (v1,lf1_1,COMPLETE);
	fprintf (stdout, "\n\n1). FITTING TREE 1 TO PARTITION 1\n", lf1_1);
	LikelihoodFunction lf1_2 = (filteredData1,secondTree);
	Optimize (res1_2, lf1_2);
	ConstructCategoryMatrix (v2,lf1_2,COMPLETE);
	fprintf (stdout, "\n\n2). FITTING TREE 2 TO PARTITION 1\n", lf1_2);
	LikelihoodFunction lf2_1 = (filteredData2,firstTree);
	Optimize (res2_1, lf2_1);
	ConstructCategoryMatrix (v3,lf2_1,COMPLETE);
	fprintf (stdout, "\n\n3). FITTING TREE 1 TO PARTITION 2\n", lf2_1);
	LikelihoodFunction lf2_2 = (filteredData2,secondTree);
	Optimize (res2_2, lf2_2);
	ConstructCategoryMatrix (v4,lf2_2,COMPLETE);
	fprintf (stdout, "\n\n4). FITTING TREE 2 TO PARTITION 2\n", lf2_2);
	
	LRT1 = 2*(res1_1[1][0]-res1_2[1][0]);
	LRT2 = 2*(res2_2[1][0]-res2_1[1][0]);
	
	fprintf (stdout, "\n\nSUMMARY:\n\n\tPARTITION 1 LRT = ",LRT1, "\n\tPARTITION 2 LRT = ", LRT2);
	if (LRT1 < 0 || LRT2 < 0)
	{
		fprintf (stdout, "\n\nERROR: Both LRTs were expeceted to be positive.\nPlease check your trees and partition");
		return 0;
	}
	
	test1 = testLRT (v1,v2)%0;
	
	for (k=0; k<Columns(v1); k=k+1)
	{
		if (test1[k]>0)
		{
			break;
		}
	}	
	
	fprintf (stdout, "\n\nEstimated p-value for the FIRST partition (based on ", itCount," replicates): ", Format(k/itCount,10,4));
	
	
	OpenWindow (CHARTWINDOW,{{"Simulated LRT for the First Partition"}
		{"labels"}
		{"test1"}
		{"Step Plot"}
		{"Index"}
		{labels[0]}
		{"Index"}
		{"LRT"}
		{"LRT"}
		{""}
		{""+LRT1}
		{""}
		{"10;1.309;0.785398"}
		{"Times:12:0;Times:10:0;Times:14:2"}
		{"0;0;16777215;5000268;0;0;6750054;11842740;13158600;14474460;0;3947580;16777215;5000268;11776947;10066329;9199669;7018159;1460610;16748822;11184810;14173291;14173291"}
		{"16"}
		},
		"(SCREEN_WIDTH-30)/2;(SCREEN_HEIGHT-50)/2;10;40");
	
	test2 = testLRT (v4,v3)%0;

	for (k=0; k<Columns(v3); k=k+1)
	{
		if (test2[k]>0)
		{
			break;
		}
	}	
	
	fprintf (stdout, "\n\nEstimated p-value for the SECOND partition (based on ", itCount," replicates): ", Format(k/itCount,10,4),"\n");
	OpenWindow (CHARTWINDOW,{{"Simulated LRT for the Second Partition"}
		{"labels"}
		{"test2"}
		{"Step Plot"}
		{"Index"}
		{labels[0]}
		{"Index"}
		{"LRT"}
		{"LRT"}
		{""}
		{""+LRT2}
		{""}
		{"10;1.309;0.785398"}
		{"Times:12:0;Times:10:0;Times:14:2"}
		{"0;0;16777215;5000268;0;0;6750054;11842740;13158600;14474460;0;3947580;16777215;5000268;11776947;10066329;9199669;7018159;1460610;16748822;11184810;14173291;14173291"}
		{"16"}
		},
		"(SCREEN_WIDTH-30)/2;(SCREEN_HEIGHT-50)/2;20+(SCREEN_WIDTH)/2;40");
}
else
{
	LikelihoodFunction		 lf1 = (filteredData,firstTree);
	Optimize				 (res1, lf1);
	ConstructCategoryMatrix  (v1,lf1,COMPLETE);
	fprintf					 (stdout, "\n\n1). FITTING TREE 1 TO THE DATA\n", lf1);
	LikelihoodFunction lf2 = (filteredData,secondTree);
	Optimize				 (res2, lf2);
	ConstructCategoryMatrix  (v2,lf2,COMPLETE);
	fprintf					 (stdout, "\n\n2). FITTING TREE 2 TO THE DATA\n", lf2);

	LRT1					 = 2*(res1[1][0]-res2[1][0]);
	
	fprintf (stdout, "\n\nSUMMARY:\n\n\tLRT = ",LRT1);
	if (LRT1 < 0)
	{
		fprintf (stdout, "\n\nERROR: The LRTs was expeceted to be positive.\nPlease check your trees and alignment (or swap the trees)");
		return 0;
	}


	test1 = testLRT (v1,v2)%0;
	for (k=0; k<Columns(v1); k=k+1)
	{
		if (test1[k]>0)
		{
			break;
		}
	}	
	
	fprintf (stdout, "\n\nEstimated p-value for  (based on ", itCount," replicates): ", Format(k/itCount,10,4),"\n");
	OpenWindow (CHARTWINDOW,{{"Simulated LRT"}
		{"labels"}
		{"test1"}
		{"Step Plot"}
		{"Index"}
		{labels[0]}
		{"Index"}
		{"LRT"}
		{"LRT"}
		{""}
		{""+LRT1}
		{""}
		{"10;1.309;0.785398"}
		{"Times:12:0;Times:10:0;Times:14:2"}
		{"0;0;16777215;5000268;0;0;6750054;11842740;13158600;14474460;0;3947580;16777215;5000268;11776947;10066329;9199669;7018159;1460610;16748822;11184810;14173291;14173291"}
		{"16"}
		},
		"(SCREEN_WIDTH-30)/2;(SCREEN_HEIGHT-50)/2;20+(SCREEN_WIDTH)/2;40");
}

/*--------------------------------------------------------------------------------------*/

function testLRT (vec1, vec2)
{
	size1 = Columns(vec1);
	
	sumVec1 = {size1,1};
	jvec	= {2,size1};
	resMx1	= {itCount,1};
	resMx2	= {itCount,1};
	
	for (k=0; k<size1; k=k+1)
	{
		sumVec1 [k]	   = 1;
		jvec	[0][k] = Log(vec1[k]);
		jvec	[1][k] = Log(vec2[k]);
	}
	
	
	for (k=0; k<itCount; k=k+1)
	{
		resampled = Random(jvec,1);
		resampled = resampled*sumVec1;
		resMx1[k] = resampled[0];
		resMx2[k] = resampled[1];
	}
	
	return (resMx1-resMx2)*2;
}
