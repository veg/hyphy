RequireVersion ("0.9920061122");

#include "SequentialAddition.ibf";

/*---------------------------------------------------------------*/

ChoiceList (dataType,"Data type",1,SKIP_NONE,"Nucleotide/Protein","Nucleotide or amino-acid (protein).",
				     "Codon","Codon (several available genetic codes).");

if (dataType<0) 
{
	return;
}


MESSAGE_LOGGING = 0;
twHasBeenOpened = 0;

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

ChoiceList (randomOption,"Addition Order",1,SKIP_NONE,
			"Given Order","The sequences will be added in the order they appear in the data file.",
			"Random Order","The order of addition will be random.");


if (randomOption<0)
{
	return;
}

ChoiceList (haveTreeConstraint,"Topology Constraint",1,SKIP_NONE,
			"No Constraint","No restrictions on topology.",
			"Use Constraint","Use a topological constraint during tree searches.");


if (haveTreeConstraint<0)
{
	return;
}

if (haveTreeConstraint)
{
	SetDialogPrompt ("Please select a topology constraint file:");
	fscanf (PROMPT_FOR_FILE, "String", _topologyPatternString);
	Tree _topologyPattern = _topologyPatternString;
}


ChoiceList (doNNIOption,"Branch Swapping",1,SKIP_NONE,
			"No Swapping","No branch swapping is performed.",
			"Complete NNI","Nearest neighbor interchange is performed after EACH sequence is added. Order (sequences)^2 additional trees are examined.",
			"Complete SPR","Subtree pruning and regrafting is performed after EACH sequence is added. Order (sequences)^3 additional trees are examined.",
			"Global NNI","Nearest neighbor interchange is performed after ALL the sequences have been added. Order (sequences)^1 additional trees are examined.",
			"Global SPR","Subtree pruning and regrafting is performed after ALL the sequences have been added. Order (sequences)^2 additional trees are examined.",
			"NNI+SPR","Nearest neighbor interchange is performed after EACH sequence is added. Subtree pruning and regrafting performed on the final tree. Order (sequences)^2 additional trees are examined.");


if (doNNIOption<0)
{
	return;
}

if (doNNIOption == 1)
{
	nniPeriod = 0;
	while (nniPeriod <= 0)
	{
		fprintf (stdout, "\nDo NNI every time this many branches are added (>=1):");
		fscanf  (stdin, "Number", nniPeriod);
	}
}

if (dataType)
{
	DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
}
else
{
	DataSetFilter filteredData = CreateFilter (ds,1,"","");
}

SelectTemplateModel(filteredData);

ExecuteAFile ("globalChecker.ibf");

ChoiceList (methodIndex,"Starting 3 taxa tree",1,SKIP_NONE,
			"First 3","The starting 3 taxa tree will comprise first 3 taxa from the data file.",
			"Choose 3","User selects 3 taxa for the starting 3 taxa tree.",
			"Best 3","The starting 3 taxa tree will be chosen by selecting the ML estimation among all possible 3 taxa trees. Warning: there are O(n^3) 3 taxa trees for n sequences.",
			"Random","Select 3 random starting sequences.");

/* begin by selecting the best 3-taxa tree */
if (methodIndex<0)
{
	return;
}

first3Taxa = {3,1};

if (methodIndex == 0)
{
	first3Taxa = first3Taxa["_MATRIX_ELEMENT_ROW_"];
}
else
{
	if (methodIndex == 1)
	{
		ChoiceList (first3Taxa, "Choose 3 taxa for the starting tree:",3,SKIP_NONE,ds);
		if (first3Taxa[0]<0)
		{
			return ;
		}
	}
}

if (pCount > 0)
{

	ChoiceList (globalParameters,"Global parameters",1,SKIP_NONE,
				"Estimate always", "Re-estimate global model parameters (e.g. rate variation parameters, substitution biases etc) for each tree. This option can be quite slow, and global parameter estimates may be unreliable for small trees, leading to possible biases.",
				"Get from a user tree", "Estimate global model parameters (e.g. rate variation parameters, substitution biases etc) from a user tree (e.g. NJ), and hold them constant for the rest of the search. This option has the advantage of big speed gains, and is based on the assumption that global model parameters are robust to some errors in the tree. This assumption could be wrong in pathological cases, however."
				);
		
			
	if (globalParameters < 0)
	{
		return ;
	}

	if (globalParameters == 1)
	{
		fprintf					(stdout, "\n[WILL USE GLOBAL ESTIMATES FROM A USER-PROVIDED TREE]\n");
		ExecuteAFile 			("queryTree.bf");
		fprintf					(stdout, "\n[OBTAINING GLOBAL PARAMETER ESTIMATES]\n");
		LikelihoodFunction 		apprxLF = (filteredData, givenTree);
		Optimize 			    (arg, apprxLF);
		
		for (k=0; k<pCount; k=k+1)
		{
			ExecuteCommands ("_param_val = " + globalParamList[k] + ";");
			fprintf (stdout, "\t", globalParamList[k], " = ", _param_val, "\n");
			ExecuteCommands (globalParamList[k] + ":=" + _param_val + ";");
		}
	}
}

l = InferTreeTopology (1.0);

if (l)
{
	fprintf (stdout,"\n\n --------------------- RESULTS --------------------- \n\n");
	fprintf (stdout,"BestTree =", bestTree);

	Tree				Inferred_Tree = bestTree;
	LikelihoodFunction  lf 			  = (filteredData, Inferred_Tree);

	Optimize (res,lf);
	fprintf (stdout, "\n",lf,"\n\n***********Save this tree to a file (y/n)?");

	fscanf  (stdin, "String", resp);
	OpenWindow (TREEWINDOW, {{"Inferred_Tree"}});

	if ((resp!="n")&&(resp!="N"))
	{
		SetDialogPrompt ("Write tree string to:");
		fprintf (PROMPT_FOR_FILE,CLEAR_FILE,bestTree,";");
	}

	saveTreeNodes = {2*(ds.species+1),3};
	for (i=2*ds.species+1;i>=0;i=i-1)
	{
		saveTreeNodes[i][0] = bestTreeNodes[i][0];
		saveTreeNodes[i][1] = bestTreeNodes[i][1];
	}
}




