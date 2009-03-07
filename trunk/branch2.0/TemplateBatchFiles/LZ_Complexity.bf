SetDialogPrompt ("Please specify a nucleotide or amino-acid data file:");

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter (ds,1);

fprintf (stdout, "\nRead ",  filteredData.species, " sequences.\n");

GetInformation (actualSequences,filteredData);

LZContent = {Rows(actualSequences),3};

labels    = {1,4};
seqNameLabel = "Sequence Name";
labels [0] = "Sequence Length";
labels [1] = "Character entropy";
labels [2] = "Lempel-Ziv Complexity";

maxLength = 14;

for (seqCounter = 0; seqCounter < Rows(actualSequences); seqCounter = seqCounter+1)
{
	GetString (seqName, filteredData, seqCounter);
	seqLength = Abs(seqName);
	if (seqLength > maxLength)
	{
		maxLength = seqLength;
	}
	seqNameLabel = seqNameLabel + ";" + seqName;
	aSeq = actualSequences[seqCounter];
	for (siteCounter = Abs(aSeq)-1; siteCounter>0; siteCounter=siteCounter-1)
	{
		if (aSeq[siteCounter]!="?")
		{
			break;
		}
	}
	DataSetFilter seqFilter = CreateFilter (ds,1,"",speciesIndex==seqCounter);
	HarvestFrequencies (seqFreq,seqFilter,1,1,0);
	seqFreq = Transpose (seqFreq)*(Log(seqFreq)*(1/Log(2)));
	actualSequences[seqCounter] = aSeq[0][siteCounter];
	LZContent [seqCounter][0] = siteCounter;
	LZContent [seqCounter][1] = -seqFreq[0];
	LZContent [seqCounter][2] = Exp (actualSequences[seqCounter]);
}

labels[3] = seqNameLabel;
OpenWindow (CHARTWINDOW,{{"Sequence Complexity"}
						   {"labels"},
						   {"LZContent"},
						   {"Bar Chart"},
						   {"Index"},
						   {labels[2]},
						   {""},
						   {""},
						   {labels[2]},
						   {""}},
						   "SCREEN_WIDTH-300;SCREEN_HEIGHT-200;100;100");


stringOfSpaces = "";
stringOfSpaces * maxLength;
for (seqCounter = 0; seqCounter < maxLength; seqCounter = seqCounter + 1)
{
	stringOfSpaces * " ";
}
stringOfSpaces * 0;

fprintf (stdout, "\nSequence Name", stringOfSpaces[0][maxLength-14], "\tLength\tEntropy\tLZ Complexity\n");

for (seqCounter = 0; seqCounter < Rows(actualSequences); seqCounter = seqCounter+1)
{
	GetString (seqName, filteredData, seqCounter);
	fprintf (stdout, seqName, stringOfSpaces[0][maxLength-Abs(seqName)-1], "\t", 
			Format (LZContent[seqCounter][0],6,0), "\t",
			Format (LZContent[seqCounter][1],7,5), "\t",
			Format (LZContent[seqCounter][2],12,4), "\n"
			);
}
