NICETY_LEVEL = 3;

#include "TemplateModels/chooseGeneticCode.def";
#include "simpleBootstrap.bf";

SetDialogPrompt ("Please specify a codon data file:");

COUNT_GAPS_IN_FREQUENCIES = 0;
VERBOSITY_LEVEL = 1;

DataSet 	  ds 		   = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);

fprintf (stdout,"\n______________READ THE FOLLOWING DATA______________\n",ds);

SelectTemplateModel(filteredData);

_DO_TREE_REBALANCE_ = 1;
#include "queryTree.bf";

if (modelType)
{
	ChoiceList (branchLengths, "Branch Lengths", 1, SKIP_NONE,
							   "Estimate", "Estimate branch lengths by ML",
							   "Proportional to input tree", "Branch lengths are proportional to those in input tree");
				 				  
	if (branchLengths < 0)
	{
		return;
	}
	
	if (branchLengths == 1)
	{
		global treeScaler = 1;
		ReplicateConstraint ("this1.?.?:=treeScaler*this2.?.?__", givenTree, givenTree);
	}
}

LikelihoodFunction lf = (filteredData,givenTree);

Optimize	(res,lf);

fprintf (stdout, "\n______________RESULTS______________\n",lf);

/* compute syn and non-syn stencils for current genetic code */

#include "categoryEcho.bf";

GetString 				(sendMeBack,lf,-1);
sendMeBack["LogL"] 		= res[1][0];
sendMeBack["NP"] 		= res[1][1];

return sendMeBack;
