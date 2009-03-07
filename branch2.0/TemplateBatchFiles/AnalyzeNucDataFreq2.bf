SetDialogPrompt ("Please specify a nuleotide data or aminoacid file:");

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter (ds,1);

global estPiA = 0.25;
global estPiC = 0.25;
global estPiG = 0.25;
global estPiT = 0.25;
global sum;

estPiA:<1;
estPiC:<1;
estPiG:<1;
estPiT:=Abs(1-estPiA-estPiC-estPiG);
sum := estPiA+estPiC+estPiG+estPiT;

global estPiAR = 0.25;
global estPiCR = 0.25;
global estPiGR = 0.25;
global estPiTR = 0.25;
global sumR;

estPiAR:<1;
estPiCR:<1;
estPiGR:<1;
estPiTR:=Abs(1-estPiAR-estPiCR-estPiGR);
sumR := estPiAR+estPiCR+estPiGR+estPiTR;

EMBED_FREQUENCY_DEPENDENCE = 1;

SelectTemplateModel(filteredData);

estPiA = .25;
estPiC = .25;
estPiG = .25;

vectorOfFrequencies     = {{estPiA/sum},{estPiC/sum},{estPiG/sum},{estPiT/sum}};
vectorOfFrequenciesRoot = {{estPiAR/sumR},{estPiCR/sumR},{estPiGR/sumR},{estPiTR/sumR}};

fprintf (stdout,"\n",ds);

_DO_TREE_REBALANCE_ = 0;
ACCEPT_ROOTED_TREES = 1;
#include "queryTree.bf";

LikelihoodFunction3 lf = (filteredData,givenTree,vectorOfFrequenciesRoot);

Optimize (res,lf);

fprintf (stdout, "\n",lf);

HarvestFrequencies (obsFreq,filteredData,1,1,1);

fprintf (stdout, "\nEstimated Root Frequencies:\n", vectorOfFrequenciesRoot, 
				 "\nEstimated Matrix Frequencies:\n", vectorOfFrequencies, 
				 "\nObserved Frequencies\n", obsFreq);

