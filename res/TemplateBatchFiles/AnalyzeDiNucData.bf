SetDialogPrompt ("Please specify a di-nucleotide (e.g. stem RNA) file:");

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter (ds,2);

fprintf (stdout,"\n______________READ THE FOLLOWING DATA______________\n",ds);

SelectTemplateModel(filteredData);

_DO_TREE_REBALANCE_ = 1;

#include "queryTree.bf";

LikelihoodFunction lf = (filteredData,givenTree);

timer = Time(0);

USE_ADAPTIVE_VARIABLE_STEP = 0;

Optimize (res,lf);

timer = Time(0)-timer;


fprintf (stdout, "\n______________RESULTS______________\nTime taken = ", timer, " seconds\nAIC Score = ", 
				  2(res[1][1]-res[1][0]),"\n",lf);

#include "categoryEcho.bf";
