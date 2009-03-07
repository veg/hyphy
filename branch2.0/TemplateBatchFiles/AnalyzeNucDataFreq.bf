SetDialogPrompt ("Please specify a nuleotide data or aminoacid file:");

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter (ds,1);

global estPiA;
global estPiC;
global estPiG;
global estPiT;
global sum;

estPiA:<1;
estPiC:<1;
estPiG:<1;
estPiT:=Abs(1-estPiA-estPiC-estPiG);
sum := estPiA+estPiC+estPiG+estPiT;

EMBED_FREQUENCY_DEPENDENCE = 1;

SelectTemplateModel(filteredData);

estPiA = .25;
estPiC = .25;
estPiG = .25;

vectorOfFrequencies = {{estPiA/sum},{estPiC/sum},{estPiG/sum},{estPiT/sum}};

fprintf (stdout,"\n",ds);

_DO_TREE_REBALANCE_ = 1;

#include "queryTree.bf";

LikelihoodFunction lf = (filteredData,givenTree);

Optimize (res,lf);

fprintf (stdout, "\n",lf);

HarvestFrequencies (obsFreq,filteredData,1,1,1);

fprintf (stdout, "\nEstimated Frequencies:", vectorOfFrequencies, "\nObserved Frequencies", obsFreq);

GetString (mbl, TrNModel,-1);

fprintf (stdout,"\n",mbl,"\n");
