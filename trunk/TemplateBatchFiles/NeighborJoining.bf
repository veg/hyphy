ChoiceList (distanceChoice, "Distance Computation",1,SKIP_NONE,
			"Distance formulae","Use one of the predefined distance measures based on data comparisons. Fast.",
			"Full likelihood","Estimate distances using pairwise MLE. More choices but slow.",
			"Use existing matrix","Load a HyPhy distance matrix for the data file");
			
if (distanceChoice < 0)
{
	return 0;
}

#include "distanceMethodNPBootstrap.bf";

if (distanceChoice == 2)
{
	IS_NPBOOTSTRAP_AVAILABLE = 0;
}

ChoiceList (dataType,"Data type",1,SKIP_NONE,"Nucleotide/Protein","Nucleotide or amino-acid (protein).",
				     "Codon","Codon (several available genetic codes).");

A_DISTANCE_METHOD   = 1;


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

if (dataType)
{
	DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
}
else
{
	DataSetFilter filteredData = CreateFilter (ds,1,"","");
}

if (distanceChoice == 1)
{
	SelectTemplateModel(filteredData);
}

ChoiceList (methodIndex,"Negative Branch Lengths",1,SKIP_NONE,
			"Keep Negative","Negative Branch Lengths are Allowed.",
			"Force Zero","Negative Branch Lengths are Forced to 0.");			
if (methodIndex<0)
{
	return;
}


function InferTreeTopology(verbFlag)
{
	distanceMatrix = {ds.species,ds.species};

	if (distanceChoice)
	{
		if (distanceChoice == 2)
		{
			SetDialogPrompt ("Load the distance matrix");
			fscanf (PROMPT_FOR_FILE,"NMatrix",distanceMatrix);
			if ((Rows(distanceMatrix) != ds.species)||(Columns(distanceMatrix) != ds.species))
			{
				fprintf (stdout, "\nThe dimensions of the distance matrix are incompatible with the data set.\n");
				return  0;
			}
		}
		else
		{
			if (verbFlag)
			{
				fprintf (stdout,"\nHYPHY Kernel is computing pairwise maximum likelihood distance estimates. A total of ", Format(ds.species*(ds.species-1)/2,0,0),
							   " estimations will be performed.\n");
			}

			_pddVF = 1;
			ExecuteAFile ("pairwiseDistanceEstimator.ibf");
		}
	}
	else
	{
		#include "chooseDistanceFormula.def";
		
		dummy = InitializeDistances (0);
		
		for (i = 0; i<ds.species; i=i+1)
		{
			for (j = i+1; j<ds.species; j = j+1)
			{
				distanceMatrix[i][j] = ComputeDistanceFormula (i,j);
			}
		}
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
	
	return 1.0;
}


DISTANCE_PROMPTS  = 1;
p = InferTreeTopology(1.0);
DISTANCE_PROMPTS = 0;

/* now with the treeNodes matrix ready we can convert it into a Newick string */

treeString = TreeMatrix2TreeString (1);

fprintf (stdout,"\n\n --------------------- INFERRED TREE --------------------- \n\n", treeString);

fprintf (stdout, "\n\n***********Save this tree to a file (y/n)?");

fscanf  (stdin, "String", resp);

if ((resp!="n")&&(resp!="N"))
{
	SetDialogPrompt ("Write tree string to:");
	fprintf (PROMPT_FOR_FILE,CLEAR_FILE,treeString,";");
}

UseModel (USE_NO_MODEL);

Tree Inferred_Tree = treeString; 

bestTreeNodes = treeNodes;

