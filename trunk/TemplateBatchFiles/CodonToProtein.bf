function 	  setUpCodonToAA (dataSetID)
{
	codonToAAMap = {};
	codeToAA 	 = "FLIMVSPTAYXHQNKDECWRG";
	
	nucChars = "ACGT";
	
	for (p1=0; p1<64; p1=p1+1)
	{
		codon = nucChars[p1$16]+nucChars[p1%16$4]+nucChars[p1%4];
		ccode = _Genetic_Code[p1];
		codonToAAMap[codon] = codeToAA[ccode];
	}
	
	codonToAAMap["---"] = "-";
	
	ExecuteCommands ("DataSetFilter _converterfilteredData 	= CreateFilter	("+dataSetID+",1);");
	ExecuteCommands ("DataSetFilter _converterfilteredDataC 	= CreateFilter	("+dataSetID+",3);");
	GetInformation (theSequences,_converterfilteredData);
	GetDataInfo    (siteToPatternMap,_converterfilteredDataC);
	return 0;
}

/*--------------------------------------------------------------------------------------------*/

function 	  doTheMapping (dummy)
{
	outSequences = "";
	outSequences *  (bigDataSet.sites* bigDataSet.species);
	
	freqCount = {};
	
	for (seqCounter = 0; seqCounter < bigDataSet.species; seqCounter = seqCounter+1)
	{
		
		aSeq = theSequences[seqCounter];
		seqLen = Abs(aSeq)-2;
		GetString (seqName, _converterfilteredData, seqCounter);
		translString = "";
		translString * (seqLen/3+1);
		for (seqPos = 0; seqPos < seqLen; seqPos = seqPos+3)
		{
			codon = aSeq[seqPos][seqPos+2];
			
			gap_count = codon$"[N-]";
			if (gap_count[0] >= 0)
			{	
				/* handle cases where codon contains one or two gap characters - this was yielding 'F' in original script */
				prot = "?";
			}
			else
			{
				prot = codonToAAMap[codon];
				if (Abs(prot) == 0)
				{
					/* 
						see if we can map this presumed ambiguitiy to a single 
						amino-acid
					*/
					GetDataInfo 			(mappedToCodon, _converterfilteredDataC, seqCounter, siteToPatternMap[seqPos$3]);
					resolutionMapping 	  = {21,1};
					for (resID = 0; resID < 64; resID = resID + 1)
					{
						if (mappedToCodon[resID])
						{
							resolutionMapping[_Genetic_Code[resID]] = 1;
						}
					}
					
					if ((+resolutionMapping) == 1)
					{
						prot = codeToAA[((Transpose(resolutionMapping))["_MATRIX_ELEMENT_COLUMN_"])[0]];
					}
					else
					{
						prot = "?";
					}
					
					codonToAAMap[codon] = prot; 
				}
			}
			freqCount[codon] = freqCount[codon]+1;
			translString * prot;
		} 
		translString * 0;
		stopPos = translString||"X";
		if (stopPos[0]>=0)
		{
			stopCount = Rows(stopPos)$2;
			fprintf (stdout, "Warning: sequence ", seqName, " has ", stopCount, " stop codons.\n");  
		}	
		outSequences * (">" + seqName + "\n" + translString + "\n");
	}
	
	outSequences * 0;
	return outSequences;
}

/*--------------------------------------------------------------------------------------------*/

if (_runAsFunctionLibrary != 1)
{
	LoadFunctionLibrary ("chooseGeneticCode.def");
	
	SetDialogPrompt ("Please choose a data file to convert:");
	
	ChoiceList (SKIP_OMISSIONS,"Deletions",1,SKIP_NONE,"Keep Deletions","Deletions will NOT be filtered from the data.",
							   "Skip Deletions","Deletions will be pruned out and NOT included in the converted file.");
	
	if (SKIP_OMISSIONS<0)
	{
		fprintf (stdout, "\n\nExecution Canceled...");
		return 0;
	}
	
	havePartitions = 0;
	
	DataSet        bigDataSet 		= ReadDataFile 	(PROMPT_FOR_FILE);
	_pCount = Rows(NEXUS_FILE_TREE_MATRIX);
	if (_pCount == Columns(DATA_FILE_PARTITION_MATRIX) && Columns(DATA_FILE_PARTITION_MATRIX) && _pCount > 1)
	{
		_runAsFunctionLibrary   = 1;
		ExecuteAFile			("_MFReader_.ibf");
		readNEXUS_partitionInfo (0);
		_runAsFunctionLibrary   = 0;
		havePartitions 			= 1;
	}

	setUpCodonToAA ("bigDataSet");
	
	
	fprintf (stdout, "Read ", bigDataSet.species, " sequences with ", bigDataSet.sites, " sites.");
	fprintf (stdout, "\nRead:\n", bigDataSet);
	
	outSequences = doTheMapping (0);
	 
	sht 		 = IS_TREE_PRESENT_IN_DATA;
	sdt 		 = DATAFILE_TREE;
	
	DataSet        bigDataSet 				= ReadFromString   (outSequences);
	DataSetFilter  _converterfilteredData 	= CreateFilter	   (bigDataSet,1);
	
	DEFAULT_FILE_SAVE_NAME = LAST_FILE_PATH;
	stopPos = DEFAULT_FILE_SAVE_NAME $ (DIRECTORY_SEPARATOR+"([^\\"+DIRECTORY_SEPARATOR+"]+)$");
	if (stopPos [0] >= 0)
	{
		stopPos = stopPos[2];
		DEFAULT_FILE_SAVE_NAME = DEFAULT_FILE_SAVE_NAME[stopPos][Abs(DEFAULT_FILE_SAVE_NAME)] + ".prot";
	}
	
	if (havePartitions)
	{
		DATA_FILE_PRINT_FORMAT  = 6;
		IS_TREE_PRESENT_IN_DATA = 0;
		DATAFILE_TREE			= 0;	
	}
	else
	{
		IS_TREE_PRESENT_IN_DATA = sht;
		DATAFILE_TREE 			= sdt;
	}
	
	
	SetDialogPrompt ("Save translated file to:");
	fprintf (PROMPT_FOR_FILE, CLEAR_FILE, _converterfilteredData);
	
	if (havePartitions && Abs (myTrees) == _pCount && Abs (splitCoordinates) == _pCount)
	{
		fprintf (LAST_FILE_PATH, "\nBEGIN TREES;");
		for (k=0; k<_pCount; k=k+1)
		{
			fprintf (LAST_FILE_PATH, "\n\tTree T_", k+1, " = ", myTrees[k], ";");
		}
		fprintf (LAST_FILE_PATH, "\nEND TREES;\nBEGIN ASSUMPTIONS;\n");
		for (k=0; k<_pCount; k=k+1)
		{
			fprintf (LAST_FILE_PATH, "\n\tCHARSET P_", k+1, " = ", 1+(splitCoordinates[k])[0]$3, "-", 1+(splitCoordinates[k])[1]$3, ";");
		}
		fprintf (LAST_FILE_PATH, "\nEND ASSUMPTIONS;\n");
	}
}
