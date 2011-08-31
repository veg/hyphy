#include "TemplateModels/chooseGeneticCode.def";

#include "simpleBootstrap.bf";

SetDialogPrompt ("Please specify a codon data file:");

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);

DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);

fprintf (stdout,"\n______________READ THE FOLLOWING DATA______________\n",ds);

SelectTemplateModel(filteredData);

_DO_TREE_REBALANCE_ = 1;

#include "queryTree.bf";


cs = 0;
lfString = "";
lfString * 128;

partCount = 10;
perProc = filteredData.sites$(partCount);

fprintf (stdout,"\n",perProc, " sites / processor\n");
for (k=0; k<partCount; k=k+1)
{
	if (k < partCount-1)
	{
		upTo = cs+perProc;
	}
	else
	{
		upTo = filteredData.sites;
	}
	fc = "DataSetFilter filteredData_"+ k + "=CreateFilter (ds,3,\""+(cs*3)+"-"+(3*upTo-1)+"\",\"\",GeneticCodeExclusions);";
	ExecuteCommands (fc);
	ExecuteCommands ("Tree tree_"+ k + "=treeString;");
	if (k)
	{
		ExecuteCommands ("ReplicateConstraint(\"this1.?.?:=this2.?.?\",tree_"+k+",tree_0);");
		lfString * ",";
	}	
	cs = cs+perProc;
	lfString * ("filteredData_"+k+",tree_"+k);
}

lfString * 0;

ExecuteCommands ("LikelihoodFunction lf = ("+lfString+");");

mxF81 = {{*,t,t,t}{t,*,t,t}{t,t,*,t}{t,t,t,*}};
HarvestFrequencies (nucFreq,ds,1,1,1);
Model M81 = (mxF81,nucFreq,1);
Tree nucTree = treeString;
DataSetFilter nucData = CreateFilter (ds,1);
LikelihoodFunction nlf = (nucData,nucTree);
Optimize (nuc_ref,nlf);

ReplicateConstraint ("this1.?.synRate:=3*this2.?.t",tree_0,nucTree);
ClearConstraints (tree_0);

USE_LAST_RESULTS = 1;

Optimize (res,lf);

fprintf (stdout, "\n______________RESULTS______________\n",lf);

#include "categoryEcho.bf";
