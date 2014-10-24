SetDialogPrompt ("Please choose a data file to convert:");

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
DataSetFilter dsf = CreateFilter (bigDataSet,1);

fprintf (stdout, "\nRead:\n", bigDataSet);


ChoiceList (DATA_FILE_PRINT_FORMAT,"Output Format",1,SKIP_NONE,
            /* 0 */ "# sequential","Sequential format with taxa names preceded by #",
			/* 1 */ "# interleaved","Interleaved format with taxa names preceded by #",
			/* 2 */ "PHYLIP sequential","PHYLIP Sequential Format. Taxa names are chopped at 10 characters.",
			/* 3 */ "PHYLIP interleaved","PHYLIP Interleaved Format. Taxa names are chopped at 10 characters.",
			/* 4 */ "NEXUS sequential, labels","NEXUS sequential format with taxa labels present in the data matrix.",
			/* 5 */ "NEXUS interleaved, labels","NEXUS interleaved format with taxa labels present in the data matrix.",
			/* 6 */ "NEXUS sequential, no labels","NEXUS sequential format without taxa labels in the data matrix.",
			/* 7 */ "NEXUS interleaved, no labels","NEXUS interleaved format without taxa labels in the data matrix.",
			/* 8 */ "CSV", "Comma separated characters",
			/* 9 */ "FASTA sequential","FASTA Sequential Format.",
			/* 10 */ "FASTA interleaved","FASTA Interleaved Format.",
			/* 11 */ "PAML compatible", "PAML Compatible PHYLIP-like format");
			
if (DATA_FILE_PRINT_FORMAT<0)
{
	fprintf (stdout, "\n\nExecution Canceled...");
}
else
{
	if (DATA_FILE_PRINT_FORMAT%2 && DATA_FILE_PRINT_FORMAT < 8 || DATA_FILE_PRINT_FORMAT == 10)
	{
		fprintf (stdout,"\n How many sites per line in interleaved file:");
		fscanf  (stdin, "Number", DATA_FILE_DEFAULT_WIDTH);
		if (DATA_FILE_PRINT_FORMAT<4)
		{
			fprintf (stdout,"\n How many sites per cluster:");
			fscanf  (stdin, "Number", DATA_FILE_GAP_WIDTH);
		}
	}
	

	SetDialogPrompt ("Write to file:");
	fprintf (PROMPT_FOR_FILE, CLEAR_FILE, dsf);
}
