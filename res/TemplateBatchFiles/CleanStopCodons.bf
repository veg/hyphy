if (onlyFilterSequenceNames != 1)
{
	#include "TemplateModels/chooseGeneticCode.def";

	SetDialogPrompt ("Please choose a codon data file:");

	DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
	fprintf (stdout, "\n\nData Read:\n", ds);

	if (IS_TREE_PRESENT_IN_DATA)
	{
		fprintf (stdout, "\nTree In Data:", DATAFILE_TREE);
	}

	DataSetFilter	    all64 = CreateFilter (ds, 3, "", "");

	stopCodonTemplate   = {1,64};

	stopCodonTemplate   	= {1,64};
	nonStopCodonTemplate	= {1,64};
}
else
{
	SetDialogPrompt ("Please choose a data file:");

	DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
	fprintf (stdout, "\n\nData Read:\n", ds);

	if (IS_TREE_PRESENT_IN_DATA)
	{
		fprintf (stdout, "\nTree In Data:", DATAFILE_TREE);
	}

	DataSetFilter	    all64 = CreateFilter (ds, 1, "", "");
}

ChoiceList (filteringOption,"Filter duplicates/gaps?",1,SKIP_NONE,"No/No",    "Keep all sequences and sites",
                                                                  "No/Yes",   "Keep all sequences, filter sites with nothing but gaps",
															      "Yes/No",   "Filter duplicate sequences but keep all sites",
                                                                  "Yes/Yes",  "Filter duplicate sequences and sites with nothing but gaps",
                                                                  "Disallow stops", "Filter duplicate sequences and all sequences that have stop codons");
															   
if (filteringOption < 0)
{
	return 0;
}
GetDataInfo (filterDimensions,all64,"CHARACTERS");
filterDimensions = Columns(filterDimensions);

if (onlyFilterSequenceNames != 1)
{
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
}
else
{
	nonStopCodonTemplate = {1,filterDimensions}["1"];
	stopCodonTemplate	 = {1,filterDimensions}["0"];
}

sequenceNames = {all64.species, 1};
doSomething 	  = 0;

validID = "[_|a-z|A-Z|0-9]+";

seqNamesList    = {};
seqNamesListAll = {};

for (k=0; k<all64.species;k=k+1)
{
	GetString(seqName, all64, k);
	newName 	= "";
	changed 	= 0;
	k2    		= Abs(seqName);
	if (k2 == 0)
	{
		newName = "Unnamed";
	}
	else
	{
		for (k3 = 0; k3 < k2; k3 = k3+1)
		{
			aChar = seqName[k3];
			if ((aChar$validID)[0] < 0)
			{
				newName = newName+"_";
				changed = 1;
			}
			else
			{
				newName = newName + aChar;
			}
		}
	}
	
	baseName = newName;
	testName = newName && 1;
	k2 = 2;
	while (seqNamesListAll[testName] > 0)
	{
		newName = baseName + "_" + k2;
		testName = newName && 1;
		k2 = k2+1;
		changed = 1;
	}
	
	if (changed)
	{
		doSomething = 1;
		sequenceNames   [k] = newName;
	}
	else
	{
		sequenceNames   [k] = seqName;
	}
	
	seqNamesListAll[testName] = 1;
}

GetInformation (sequenceData, 	 all64);
GetDataInfo    (duplicateMapper, all64);

if (onlyFilterSequenceNames != 1)
{
	replacementString = "---";
}
else
{
	replacementString = "-";
}

notDuplicate 	  = {};
duplicateChecker  = {};
haveInfoAtSites	  = {};

for (sequenceIndex = 0; sequenceIndex < all64.species; sequenceIndex = sequenceIndex+1)
{

	stopCodonCount     = 0;
	sitesWithDeletions = {1,all64.unique_sites};
	
	for (siteIndex = 0; siteIndex < all64.unique_sites; siteIndex = siteIndex+1)
	{
		GetDataInfo (siteInfo, all64, sequenceIndex, siteIndex);
		siteInfo1 = stopCodonTemplate*siteInfo;
		siteInfo2 = nonStopCodonTemplate*siteInfo;
		if (siteInfo1[0]>0 && siteInfo2[0] == 0)
		{
			sitesWithDeletions[siteIndex] = 1;
			stopCodonCount 	= stopCodonCount+1;
		}
		if (filteringOption % 2)
		{
			if (haveInfoAtSites[siteIndex] == 0)
			{
				if (siteInfo1[0]+siteInfo2[0] < stopCodonTemplate)
				{
					haveInfoAtSites[siteIndex] = 1;
				}
			}
		}
	}
	
	if (stopCodonCount > 0)
	{
        if (filterinOption == 4)
        {   
            continue;
        }
		fprintf (stdout, "\nSequence ", sequenceNames[sequenceIndex], ":"); 
		fprintf (stdout, "\n\t", Format(stopCodonCount,8,0), " stop codons found.");
		
		doSomething 		= 1;
		cleanedString 		= "";
		seqString 			= sequenceData[sequenceIndex];
		cleanedString 		* (Abs(seqString)+1);

		for (siteIndex = 0; siteIndex < all64.sites; siteIndex = siteIndex+1)
		{
			stopCodonCount = duplicateMapper[siteIndex];
			if (sitesWithDeletions[stopCodonCount])
			{
				cleanedString * replacementString;
			}
			else
			{
				cleanedString * seqString[3*siteIndex][3*siteIndex+2];
			}
		}
		cleanedString * 0;
		sequenceData[sequenceIndex] = cleanedString;
	}
	
	if (filteringOption >= 2)
	{
		if (duplicateChecker[sequenceData[sequenceIndex]] == 0)
		{	
			duplicateChecker[sequenceData[sequenceIndex]] = 1;
			notDuplicate[sequenceIndex] = 1;
		}
		else
		{
			doSomething 		= 1;
		}
	}
	else
	{
		notDuplicate[sequenceIndex] = 1;
	}
}

if (filteringOption%2)
{
	doSomething = doSomething || (Abs(haveInfoAtSites)<all64.unique_sites);
}

if (doSomething)
{
	SetDialogPrompt ("Save cleaned data to:");
	fprintf (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN);
	for (sequenceIndex = 0; sequenceIndex < all64.species; sequenceIndex = sequenceIndex+1)
	{
		if (notDuplicate[sequenceIndex])
		{
			fprintf (LAST_FILE_PATH, ">", sequenceNames[sequenceIndex], "\n", sequenceData[sequenceIndex], "\n\n");
		}
	}
	fprintf (LAST_FILE_PATH, CLOSE_FILE);
}
else
{
	fprintf (stdout, "\n\nNo stop codons found\n\n");
}

sequenceData = 0;
