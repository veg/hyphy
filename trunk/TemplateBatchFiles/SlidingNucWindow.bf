SetDialogPrompt ("Please choose a nucleotide or amino-acid data file:");
DataSet 	    ds 			 = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter   filteredData = CreateFilter (ds,1);
DataSetFilter   allData		 = CreateFilter (ds,1);

ChoiceList	  (refSequence,"Reference Sequence", 1, SKIP_NONE, filteredData);

if (refSequence < 0)
{
	return 0;
}

windowSize = -1;
while (windowSize < 1 || windowSize > filteredData.sites)
{
	fprintf (stdout, "Sliding Window Size [1-", filteredData.sites, "]:");
	fscanf (stdin,"Number", windowSize);
	windowSize = windowSize $ 1;
}

if (windowSize<filteredData.sites)
{
	windowStride = -1;
	while (windowStride < 1 || windowStride > filteredData.sites-windowSize)
	{
		fprintf (stdout, "Sliding Window Stride [1-", filteredData.sites-windowSize, "]:");
		fscanf (stdin,"Number", windowStride);
		windowStride = windowStride $ 1;
	}
}
else
{
	windowStride = 1;
}

DISTANCE_PROMPTS = 1;
#include "chooseDistanceFormula.def";
dummy = InitializeDistances (0);

totalWindows   = ((filteredData.sites-windowSize)/windowStride+0.5)$1+1;
distanceMatrix = {totalWindows,allData.species};
labelMatrix	   = {1,allData.species};
labelMatrix[0] = "Midpoint";

doneSkip = 0;

for (si = 0; si < allData.species ; si = si+1)
{
	if (si == refSequence)
	{
		doneSkip = 1;
	}
	else
	{
		GetString (seqName, allData, si);
		labelMatrix[si+1-doneSkip] = seqName;
	}
}

currentStart = 0;

for (wc = 0; wc < totalWindows; wc = wc+1)
{
	currentEnd	 = currentStart + windowSize - 1;
	fprintf (stdout, currentStart, "-", currentEnd, "\n");
	distanceMatrix[wc][0] = ((currentEnd+currentStart)*0.5+0.5)$1;
	DataSetFilter filteredData = CreateFilter (ds, 1, siteIndex >= currentStart && siteIndex < currentEnd);
	doneSkip = 0;
	for (si = 0; si < allData.species ; si = si+1)
	{
		if (si == refSequence)
		{
			doneSkip = 1;
		}
		else
		{
			k = ComputeDistanceFormula (si,refSequence);
			distanceMatrix [wc][1+si-doneSkip] = k;
		}
	}
	currentStart = currentStart + windowStride;
}

OpenWindow (CHARTWINDOW,{{"Sliding Window Distance"}
		{"labelMatrix"}
		{"distanceMatrix"}
		{"None"}
		{"Index"}
		{"None"}
		{""}
		{""}
		{""}
		{"0"}
		{""}
		{"0;0"}
		{"10;1.309;0.785398"}
		{"Times:12:0;Times:10:0;Times:12:2"}
		{"0;0;13816530;16777215;0;0;6579300;11842740;13158600;14474460;0;3947580;16777215;15670812;6845928;16771158;2984993;9199669;7018159;1460610;16748822;11184810;14173291"}
		{"16,0,0"}
		},
		"1046;657;70;70");


	
