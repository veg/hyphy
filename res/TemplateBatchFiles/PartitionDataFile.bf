SetDialogPrompt ("Please choose a data file to partition:");

ChoiceList (SKIP_OMISSIONS,"Deletions",1,SKIP_NONE,"Keep Deletions","Deletions will NOT be filtered from the data.",
			"Skip Deletions","Deletions will be pruned out and NOT included in the converted file.");


if (SKIP_OMISSIONS<0)
{
	fprintf (stdout, "\n\nExecution Canceled...");
}

DataSet bigDataSet = ReadDataFile (PROMPT_FOR_FILE);

if (IS_TREE_PRESENT_IN_DATA)
{
	fprintf (stdout, "\nTree In Data:", DATAFILE_TREE);
}
fprintf (stdout, "\nRead:\n", bigDataSet);

ChoiceList (DATA_FILE_PRINT_FORMAT,"Output Format",1,SKIP_NONE,"# sequential","Sequential format with taxa names preceded by #",
			"# interleaved","Interleaved format with taxa names preceded by #",
			"PHYLIP sequential","PHYLIP Sequential Format. Taxa names are chopped at 10 characters.",
			"PHYLIP interleaved","PHYLIP Interleaved Format. Taxa names are chopped at 10 characters.",
			"NEXUS sequential, labels","NEXUS sequential format with taxa labels present in the data matrix.",
			"NEXUS interleaved, labels","NEXUS interleaved format with taxa labels present in the data matrix.",
			"NEXUS sequential, no labels","NEXUS sequential format without taxa labels in the data matrix.",
			"NEXUS interleaved, no labels","NEXUS interleaved format without taxa labels in the data matrix.");

if (DATA_FILE_PRINT_FORMAT<0)
{
	fprintf (stdout, "\n\nExecution Canceled...");
}
else
{
	fprintf (stdout,"\n\nPartition starts at site (0-based):");
	fscanf  (stdin,"Number",beginPart);
	fprintf (stdout,"\nPartition ends at site (0-based):");
	fscanf  (stdin,"Number",endPart);
	
	DataSetFilter dsf = CreateFilter (bigDataSet,1,(siteIndex>=beginPart)&&(siteIndex<=endPart));

	SetDialogPrompt ("Write to file:");
	
	fprintf (PROMPT_FOR_FILE, CLEAR_FILE, dsf);
}
