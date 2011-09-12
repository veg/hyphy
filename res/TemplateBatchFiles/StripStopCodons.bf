#include "TemplateModels/chooseGeneticCode.def";

SetDialogPrompt ("Please choose a codon data file:");

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
fprintf (stdout, "\n\nData Read:\n", ds);

if (IS_TREE_PRESENT_IN_DATA)
{
	fprintf (stdout, "\nTree In Data:", DATAFILE_TREE);
}

DataSetFilter	    all64 = CreateFilter (ds, 3, "", "");

stopCodonTemplate   	= {1,64};
nonStopCodonTemplate	= {1,64};

for (stateCount=0; stateCount<64; stateCount=stateCount+1)
{
	if (_Genetic_Code[stateCount] == 10)
	{
		stopCodonTemplate[stateCount] = 1;
	}
	else
	{
		nonStopCodonTemplate[stateCount] = 1;
	}
}

validSequences 			= "";
dummy					= validSequences*512;

sequenceNames = {all64.species, 1};

for (sequenceIndex = 0; sequenceIndex < all64.species; sequenceIndex = sequenceIndex+1)
{
	GetString (replacementString, all64, sequenceIndex);
	sequenceNames[sequenceIndex] = replacementString;
}

GetInformation (sequenceData, all64);
GetDataInfo    (duplicateMapper, all64);

goodSequences	  = 0;

for (sequenceIndex = 0; sequenceIndex < all64.species; sequenceIndex = sequenceIndex+1)
{
	fprintf (stdout, "\nSequence ", Format(sequenceIndex+1,7,0),"/",all64.species," : "); 
	stopCodonCount     = 0;
	
	for (siteIndex = 0; siteIndex < all64.unique_sites; siteIndex = siteIndex+1)
	{
		GetDataInfo (siteInfo, all64, sequenceIndex, siteIndex);
		siteInfo1 = stopCodonTemplate*siteInfo;
		siteInfo2 = nonStopCodonTemplate*siteInfo;
		if (siteInfo1[0]>0 && siteInfo2[0] == 0)
		{
			stopCodonCount 	= stopCodonCount+1;
			break;
		}
	}
	
	if (stopCodonCount)
	{
		fprintf (stdout, " STOP CODONS.");		
	}
	else
	{
		if (goodSequences == 0)
		{
			dummy	= validSequences*(""+sequenceIndex);
		}
		else
		{
			dummy	= validSequences*(","+sequenceIndex);		
		}
		fprintf (stdout, " CLEAN.");
		goodSequences = goodSequences + 1;	
	}
}

dummy	= validSequences*0;

if (goodSequences<all64.species)
{
	if (goodSequences == 0)
	{
		fprintf (stdout, "\n\nEvery sequence had stop codons!\n\n");
	}
	else
	{
		fprintf (stdout, "\n\nKept ", goodSequences, "/",all64.species," sequences.\n\n");
		DataSetFilter	    all64 = CreateFilter (ds, 1, "", validSequences);
		SetDialogPrompt ("Save cleaned data to:");
		fprintf 		(PROMPT_FOR_FILE, CLEAR_FILE, all64);
	}
}
else
{
	fprintf (stdout, "\n\nNo stop codons found\n\n");
}

sequenceData = 0;
