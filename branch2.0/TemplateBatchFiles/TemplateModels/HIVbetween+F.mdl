#include "HIVbetween.ibf";

if (SKIP_FREQ_HARVESTING == 0)
{
	HarvestFrequencies 		 (vectorOfFrequencies,filteredData,1,1,0);
}
HIVBetweenMatrix	    = 0;
MULTIPLY_BY_FREQS 	    = PopulateModelMatrix ("HIVBetweenMatrix",vectorOfFrequencies);
Model HIVBetweenModel 	= (HIVBetweenMatrix, vectorOfFrequencies, MULTIPLY_BY_FREQS);
FREQUENCY_SENSITIVE     = 1;
