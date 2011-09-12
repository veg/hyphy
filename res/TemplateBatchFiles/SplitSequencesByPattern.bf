ExecuteAFile 					(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "ReadDelimitedFiles.bf");

SetDialogPrompt 		("Please choose a data file:");

DataSet 				ds = ReadDataFile (PROMPT_FOR_FILE);
fprintf 				(stdout, "\n\nData Read:\n", ds);
fprintf					(stdout, "\nA regular expression to split with");
fscanf					(stdin,"String",regExp);

DataSetFilter	    	allData = CreateFilter (ds,1);
GetString				(taxonNames,allData,-1);

splitIDs = {};
usedNames = {};

for (k=0; k<allData.species; k=k+1)
{
	newName = normalizeSequenceID (taxonNames[k], "usedNames");
	if (newName != taxonNames[k])
	{
		SetParameter (ds,k,newName);
	}
}

GetString				(taxonNames,allData,-1);

for (k=0; k<allData.species; k=k+1)
{
	match = taxonNames[k]$regExp;
	if (match[0]>=0)
	{
		matchTo = (taxonNames[k])[match[0]][match[1]];
		if (Abs(splitIDs[matchTo]) == 0)
		{
			splitIDs[matchTo] = {};
		}
		(splitIDs[matchTo])[k+1] = 1;
	}
}

basePath = LAST_FILE_PATH;

keys = Rows(splitIDs);

for (k=0; k<Columns(keys); k=k+1)
{
	fileOut = basePath + "." + keys[k];
	DataSetFilter splitDF = CreateFilter (ds,1,"",(splitIDs[keys[k]])[speciesIndex+1]);
	fprintf (fileOut,CLEAR_FILE,splitDF);
	fprintf (stdout, keys[k], " matched ", splitDF.species, " sequences\n");
}
