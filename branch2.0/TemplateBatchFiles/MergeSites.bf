SetDialogPrompt ("Please choose the first data file:");
DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
fprintf (stdout, "\n\nData Set 1:\n", ds);

if (IS_TREE_PRESENT_IN_DATA)
{
	fprintf (stdout, "\nTree In Data:", DATAFILE_TREE);
}

SetDialogPrompt ("Please choose the second data file:");
DataSet ds2 = ReadDataFile (PROMPT_FOR_FILE);

while (ds2.species!=ds.species)
{
	ChoiceList (userChoice,"Unequal number of sequences",1,SKIP_NONE,
				"Choose Again", "Choose another datafile to merge",
				"Merge Anyway", "Pad missing sequences with deletions");
	if (userChoice<0)
	{
		return;
	}
	if (userChoice==1)
	{
		break;
	}
	DataSet ds2 = ReadDataFile (PROMPT_FOR_FILE);
}

fprintf (stdout, "\n\nData Set 2:\n", ds2);

DataSet dsCombined = Concatenate (purge,ds,ds2);
fprintf (stdout, "\n\nMerged Data Set:\n", dsCombined);

SetDialogPrompt ("Save resulting data set to:");
DataSetFilter dsf = CreateFilter (dsCombined,1);

fprintf (PROMPT_FOR_FILE,CLEAR_FILE,dsf);
