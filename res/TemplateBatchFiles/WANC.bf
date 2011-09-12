function _siteWeighter_ (treeIndex, siteIndex)
{						
	_SITE_OS_COUNT = matrixTrick*(_OBSERVED_S_$STATE_COUNT_MATRIX)*Transpose(matrixTrick);
	_SITE_ON_COUNT = matrixTrick*(_OBSERVED_NS_$STATE_COUNT_MATRIX)*Transpose(matrixTrick);
	_SITE_ES_COUNT = matrixTrick*(_PAIRWISE_S_$WSTATE_COUNT_MATRIX)*Transpose(matrixTrick);
	_SITE_EN_COUNT = matrixTrick*(_PAIRWISE_NS_$WSTATE_COUNT_MATRIX)*Transpose(matrixTrick);
	
	resultMatrix[siteIndex][0] = _SITE_OS_COUNT[0];
	resultMatrix[siteIndex][1] = _SITE_ON_COUNT[0];
	resultMatrix[siteIndex][2] = _SITE_ES_COUNT[0];
	resultMatrix[siteIndex][3] = _SITE_EN_COUNT[0];
	
	weight = _SITE_EN_COUNT[0]+_SITE_ES_COUNT[0];
	
	p = _SITE_ES_COUNT[0]/weight;
	
	resultMatrix[siteIndex][5] = p;
	
	p2 = resultMatrix[siteIndex][0]+resultMatrix[siteIndex][1];
	
	resultMatrix[siteIndex][7] = resultMatrix[siteIndex][1]/resultMatrix[siteIndex][3];
	
	if (resultMatrix[siteIndex][2])
	{
		resultMatrix[siteIndex][6] = resultMatrix[siteIndex][0]/resultMatrix[siteIndex][2];
		resultMatrix[siteIndex][8] = resultMatrix[siteIndex][7]-resultMatrix[siteIndex][6];	
		resultMatrix[siteIndex][11] = resultMatrix[siteIndex][8]/totalTreeLength;	
	}
	
	if (p2)
	{
		resultMatrix[siteIndex][4]  = resultMatrix[siteIndex][0]/p2;		
		resultMatrix[siteIndex][9]  = extendedBinTail (p2,p,resultMatrix[siteIndex][0]);
		if (resultMatrix[siteIndex][0]>=1)
		{
			resultMatrix[siteIndex][10] = 1-resultMatrix[siteIndex][9]+(resultMatrix[siteIndex][9]-extendedBinTail(p2,p,resultMatrix[siteIndex][0]-1));
		}
		else
		{
			resultMatrix[siteIndex][10] = 1-resultMatrix[siteIndex][9]+(resultMatrix[siteIndex][9]-extendedBinTail (p2,p,0));
		}
	}		


	return 0;
}


/*___________________________________________________________________________________________________________*/

function extendedBinTail (ebn, ebp, ebx)
{
	ebr = ebx$1; /* rounded to nearest integer */
	
	currentBinCoeff = (1-ebp)^ebn; /*compute the first binomial coefficient */
	
	binHead = 0;
	
	for (ebk=0; ebk<=ebr; ebk=ebk+1)
	{
		binHead			= binHead + currentBinCoeff;
		currentBinCoeff = currentBinCoeff * (ebn-ebk) / (ebk+1) * ebp / (1-ebp);
	}
	
	if (ebx <= ebn$1)
	{
		binHead = binHead + currentBinCoeff*(ebx-ebr);
	}
	else
	{
		binHead = binHead + (1-binHead)*(ebx-ebr)/(ebn-ebn$1);	
	}
		
	return binHead;
}

/*___________________________________________________________________________________________________________*/

function PadString (padLength,padChar)
{
	for (padCounter=0;padCounter<padLength;padCounter=padCounter+1)
	{
		fprintf (stdout,padChar);
	}
	return padLength;
}

/*___________________________________________________________________________________________________________*/

function	PrintASCIITable (dataMatrix, titleMatrix)
{
	columnWidths = {1,Columns(titleMatrix)};
	for (counter1=0; counter1<Columns(titleMatrix); counter1 = counter1+1)
	{
		counter2 = Abs (titleMatrix[0][counter1])+2;
		if (counter2<12)
		{
			counter2 = 12;
		}
		columnWidths[0][counter1] = counter2;
	}
	fprintf (stdout, "\n");
	for (counter2 = 0; counter2 < Columns (titleMatrix); counter2 = counter2+1)
	{
		fprintf (stdout,"+-");
		dummy = PadString (columnWidths[0][counter2],"-");
		fprintf (stdout,"-");
	}
	fprintf (stdout,"+\n| ");
	
	for (counter1=0; counter1<Columns(titleMatrix); counter1 = counter1+1)
	{
		fprintf (stdout, titleMatrix[counter1]);
		dummy = PadString (columnWidths[0][counter1]-Abs(titleMatrix[counter1])," ");
		fprintf (stdout, " | ");
	}
	
	fprintf (stdout, "\n");
	
	for (counter1=-1; counter1<Rows(dataMatrix); counter1 = counter1 + 1)
	{
		if (counter1>=0)
		{
			fprintf (stdout,"| ");
			fprintf (stdout,Format(counter1+1,columnWidths[0][0],0));
			for (counter2 = 1; counter2 < Columns (titleMatrix); counter2 = counter2+1)
			{
				fprintf (stdout," | ");
				fprintf (stdout,Format(dataMatrix[counter1][counter2-1],columnWidths[0][counter2],-1));
			}
			fprintf (stdout," ");
			fprintf (stdout, "|\n");
		}
		for (counter2 = 0; counter2 < Columns (titleMatrix); counter2 = counter2+1)
		{
			fprintf (stdout,"+-");
			dummy = PadString (columnWidths[0][counter2],"-");
			fprintf (stdout,"-");
		}
		fprintf (stdout, "+\n");
	}
	
	return 1;
}

/*___________________________________________________________________________________________________________*/

function	PrintTableToFile (dataMatrix, titleMatrix, promptOrNot)
{
	SetDialogPrompt ("Export tab separated data to:");
	
	if (promptOrNot)
	{
		fprintf (PROMPT_FOR_FILE, CLEAR_FILE);
	}
	
	fprintf (LAST_FILE_PATH, titleMatrix[0][0]);

	for (counter1=1; counter1<Columns(titleMatrix); counter1 = counter1+1)
	{
		fprintf (LAST_FILE_PATH, "\t", titleMatrix[0][counter1]);
	}
	
	for (counter1=0; counter1<Rows(dataMatrix); counter1 = counter1 + 1)
	{
		fprintf (LAST_FILE_PATH,"\n",dataMatrix[counter1][0]);
		for (counter2 = 1; counter2 < Columns (dataMatrix); counter2 = counter2+1)
		{
			fprintf (LAST_FILE_PATH,"\t",dataMatrix[counter1][counter2]);
		}
	}
	
	return 1;
}

/*----------------------------------------------------------------------------*/


if (useCustomCountingBias)
{
	incFileName = "Distances/CodonToolsMain.def";
}
else
{
	incFileName = "Distances/CodonTools.def";
}

ExecuteCommands  ("#include \""+incFileName+"\";");

/* first column - the number of synonymous sites
   second column - the number of non-syn sites */

matrixTrick = {1,stateCharCount};

for (h=Columns(matrixTrick)-1; h>=0; h=h-1)
{
	matrixTrick  [h] = 1;
}

resultMatrix = {filteredData.sites,12};

if (pipeThroughFlag == 0)
{
	likelihoodFnChoice = 0;
	if (Rows("LikelihoodFunction")>1)
	{
		ChoiceList  (likelihoodFnChoice,"Choose a Likelihood Function",1,NO_SKIP,LikelihoodFunction);
	}		

	if (likelihoodFnChoice<0)
	{
		return;
	} 
				 
	LIKELIHOOD_FUNCTION_OUTPUT = response;

	if (RESTORE_GLOBALS)
	{
		dumb = RestoreGlobalValues (likelihoodFnChoice);
	}

	GetString(likelihoodFunctionName,LikelihoodFunction,likelihoodFnChoice);
}

StateCounter (likelihoodFunctionName__,"_siteWeighter_");

if (pipeThroughFlag == 0)
{
	labelMatrix =     {1,12};
	labelMatrix[0] = "Observed S Changes";
	labelMatrix[1] = "Observed NS Changes";
	labelMatrix[2] = "E[S Sites]";
	labelMatrix[3] = "E[NS Sites]";
	labelMatrix[4] = "Observed S. Prop.";
	labelMatrix[5] = "P{S}";
	labelMatrix[6] = "dS";
	labelMatrix[7] = "dN";
	labelMatrix[8] = "dN-dS";
	labelMatrix[9] = " P{S leq. observed}";
	labelMatrix[10] = "P{S geq. observed}";
	labelMatrix[11] = "Scaled dN-dS";

	sigLevel = -1;

	while ((sigLevel<=0)||(sigLevel>=1))
	{
		fprintf (stdout, "\nSignificance level for a site to be classified as positively/negatively selected?");
		fscanf  (stdin, "Number", sigLevel);
	}

	posSelected = 0;
	negSelected = 0;

	p = Rows(resultMatrix);

	for (p2=0; p2<p; p2=p2+1)
	{
		if (resultMatrix [p2][0]+resultMatrix [p2][1]>0.5)
		{
			v = resultMatrix [p2][8];
			if (v>0)
			{
				if (resultMatrix [p2][9] < sigLevel)
				{
					posSelected = posSelected+1;
				}
			}
			else
			{
				if (v<0)
				{
					if (resultMatrix [p2][10] < sigLevel)
					{
						negSelected = negSelected+1;
					}
				}
			}
		}
	}

	selLabelMatrix = {{"Index","Site Index","dN-dS","p-value"}};

	if (posSelected)
	{
		psMatrix = {posSelected, 3};
		h = 0;
		for (p2=0; p2<p; p2=p2+1)
		{
			if (resultMatrix [p2][0]+resultMatrix [p2][1]>0.5)
			{
				v = resultMatrix [p2][8];
				if (v>0)
				{
					if (resultMatrix [p2][9] < sigLevel)
					{
						psMatrix[h][0] = p2+1;
						psMatrix[h][1] = v;
						psMatrix[h][2] = resultMatrix [p2][9];
						h = h+1;
					}
				}
			}
		}
		
		fprintf (stdout,"\n******* FOUND ", posSelected, " POSITIVELY SELECTED SITES ********\n\n");
		dummy = PrintASCIITable  (psMatrix, selLabelMatrix);
	}
	else
	{
		fprintf (stdout,"\n******* FOUND NO POSITIVELY SELECTED SITES ********\n\n");
	}

	if (negSelected)
	{
		psMatrix = {negSelected, 3};
		h = 0;
		for (p2=0; p2<p; p2=p2+1)
		{
			if (resultMatrix [p2][0]+resultMatrix [p2][1]>0.5)
			{
				v = resultMatrix [p2][8];
				if (v<0)
				{
					if (resultMatrix [p2][10] < sigLevel)
					{
						psMatrix[h][0] = p2+1;
						psMatrix[h][1] = v;
						psMatrix[h][2] = resultMatrix [p2][10];
						h = h+1;
					}
				}
			}
		}
		
		fprintf (stdout,"\n******* FOUND ", negSelected, " NEGATIVELY SELECTED SITES ********\n\n");
		dummy = PrintASCIITable  (psMatrix, selLabelMatrix);
	}
	else
	{
		fprintf (stdout,"\n******* FOUND NO NEGATIVELY SELECTED SITES ********\n\n");
	}

	ChoiceList (outputChoice, "Output Options",1,SKIP_NONE,
				"ASCII Table",	 "Output is printed to the console as an ASCII table.",
				"Export to File","Output is spooled to a tab separated file.",
				"Chart","A HYPHY chart window is displayed (Mac/Doze Only).");

	if (outputChoice<0)
	{
		return;
	}

	if (outputChoice==0)
	{
		dummy = PrintASCIITable  (resultMatrix, labelMatrix);
	}
	else
	{
		if (outputChoice == 1)
		{
			SetDialogPrompt ("Save result table to");
			dummy = PrintTableToFile  (resultMatrix, labelMatrix,1);
		}
		else
		{
			OpenWindow (CHARTWINDOW,{{"Rates by sites"}
									   {"labelMatrix"},
									   {"resultMatrix"},
									   {"Contrast Bars"},
									   {"Index"},
									   {labelMatrix[6]+";"+labelMatrix[7]},
									   {"Site Index"},
									   {"dN"},
									   {"dS"},
									   {"0"}},
									   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");
		}
	}
}
