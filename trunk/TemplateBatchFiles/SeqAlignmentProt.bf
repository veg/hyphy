ExecuteAFile 					("Utility/GrabBag.bf");

/* START ALIGNMENT SETTINGS */

#include "SeqAlignShared.ibf";

/* END ALIGNMENT SETTINGS */


/* REFERENCE SEQUENCES OPTIONS */

SetDialogPrompt 		("Sequence File:");

doLongestSequence = 	(refSeq==1);

/* build codon translation table */




if (refSeq > 1)
{
    LoadFunctionLibrary ("chooseGeneticCode");
    codonToAAMap		= defineCodonToAA();
	DataSet        unal2 			= ReadDataFile 	(PROMPT_FOR_FILE);
	refSeq = ">" + predefSeqNames [refSeq][0] + "\n" + translateCodonToAA (RefSeqs[refSeq-2], codonToAAMap, 0);
	DataSet		   refD = ReadFromString (refSeq);
	DataSet        unal = Combine (refD,unal2);
}
else
{
	ChoiceList (refSeq2,"Insert a coordinate reference sequence?",1,SKIP_NONE,predefSeqNames2);
	if (refSeq2 < 0)
	{
		return 0;
	}
	if (refSeq2)
	{
		skipCodeSelectionStep = 1;
		ApplyGeneticCodeTable (0);
		
		DataSet        unal2 			= ReadDataFile 	(PROMPT_FOR_FILE);
		refSeq = ">" + predefSeqNames [refSeq2][0] + "\n" + translateCodonToAA (RefSeqs[refSeq2-1], codonToAAMap, 0);
		DataSet		   refD = ReadFromString (refSeq);
		DataSet        unal = Combine (unal2,refD);
		if (doLongestSequence)
		{
			doLongestSequence = 2;
		}
	}
	else
	{
		DataSet        unal 			= ReadDataFile 	(PROMPT_FOR_FILE);
	}
}


DataSetFilter  filteredData 	= CreateFilter	(unal,1);

ExecuteAFile ("Utility/GrabBag.bf");

GetInformation (UnalignedSeqs,filteredData);

COUNT_GAPS_IN_FREQUENCIES = 0;
InitializeDistances		   (0);

seqCount = Rows(UnalignedSeqs) * Columns(UnalignedSeqs);
/* preprocess sequences */


refSeqID = 0;
maxL	 = 0;

for (seqCounter = 0; seqCounter < seqCount; seqCounter = seqCounter+1)
{
	aSeq = UnalignedSeqs[seqCounter];
	UnalignedSeqs[seqCounter] = aSeq^{{"[^a-zA-Z]",""}};
	if (doLongestSequence > 0 && seqCounter < seqCount - (refSeq2 > 0))
	{
		if (Abs (UnalignedSeqs[seqCounter]) > maxL)
		{
			maxL	 = Abs (UnalignedSeqs[seqCounter]);
			refSeqID = seqCounter;
		}
	}
}

refLength 		 = Abs(UnalignedSeqs[refSeqID]);
s1 = UnalignedSeqs[refSeqID];

if (doLongestSequence > 0)
{
	GetString (refname, filteredData, refSeqID);
	fprintf (stdout, "Set the reference sequence to ", refname, " (length = ", refLength, ")\n");
}

SeqAlignments 	 = {};
startingPosition = {seqCount,2};
refInsertions	 = {refLength,1};

fprintf (stdout,"\nPerforming pairwise alignment with reference sequences\n");

//pairwiseDistances = {seqCount,1};

for (seqCounter = 0; seqCounter < seqCount; seqCounter = seqCounter+1)
{
	if (seqCounter == refSeqID)
	{
		continue;
	}
	
	s2 			 = UnalignedSeqs[seqCounter];
	inStr 		 = {{s1,s2}};
	AlignSequences(aligned, inStr, alignOptions);
	
	aligned = aligned[0];
	SeqAlignments[seqCounter] = aligned;
	aligned = aligned[1];
	
	myStartingPosition = aligned$"[^-]";
	myEndingPosition  = Abs (aligned)-1;
	while (aligned[myEndingPosition]=="-")
	{
		myEndingPosition = myEndingPosition - 1;
	}
	myStartingPosition = myStartingPosition[0];
	startingPosition[seqCounter][0] = myStartingPosition;
	startingPosition[seqCounter][1] = myEndingPosition;
	aligned = aligned[myStartingPosition][myEndingPosition];
	cleanedAlignedSeq = ((SeqAlignments[seqCounter])[2])[myStartingPosition][myEndingPosition]^{{"N","-"}};
	pairOfSequences = ">1\n" + cleanedAlignedSeq + "\n>2\n" + aligned;
	refInsert = aligned||"-+";
	if (refInsert[0]>0)
	{
		insCount = Rows (refInsert)/2;
		offset = 0;
		for (insN = 0; insN < insCount; insN = insN+1)
		{
			insPos 		= refInsert[insN*2];
			insLength	= refInsert[insN*2+1]-insPos+1;
			insPos 		= insPos-offset;
			if (refInsertions[insPos]<insLength)
			{
				refInsertions[insPos]=insLength;
			}
			offset = offset + insLength;
		}
	}
	SetParameter (STATUS_BAR_STATUS_STRING, "Performing pairwise alignment with reference sequences ("+seqCounter+"/"+seqCount+" done)",0);
}

/* produce a fully gapped reference sequence */

fprintf (stdout,"\nMerging pairwise alignments into a MSA\n");

fullRefSeq = "";
fullRefSeq * refLength;
fullRefSeq * (s1[0]);


for (seqCounter=1;seqCounter<refLength;seqCounter=seqCounter+1)
{
	gapCount = refInsertions[seqCounter];
	for (k=0; k<gapCount;k=k+1)
	{
		fullRefSeq*("-");
	}	
	fullRefSeq  * (s1[seqCounter]);
}

fullRefSeq * 0;

refLength = Abs(fullRefSeq);

SetDialogPrompt ("Save alignment to:");

GetString (seqName,unal,refSeqID);
fprintf (PROMPT_FOR_FILE,CLEAR_FILE,">",seqName,"\n",fullRefSeq);
fName = LAST_FILE_PATH;

for (seqCounter = 0; seqCounter < seqCount; seqCounter = seqCounter+1)
{
	if (seqCounter == refSeqID)
	{
		continue;
	}
	GetString (seqName,unal,seqCounter);
	aligned = SeqAlignments[seqCounter];
	
	aligned1 = aligned[1];
	aligned2 = aligned[2];
	
	s2 = startingPosition[seqCounter][0];
	e2 = startingPosition[seqCounter][1];
	
	gappedSeq = "";
	gappedSeq * Abs(aligned2);

	
	k=0;
	
	while (k<refLength)
	{
		while (fullRefSeq[k]!=aligned1[s2] && k < refLength)
		{
			gappedSeq*("-");
			k=k+1;
		}
		if (k == refLength)
		{
			break;
		}
		gappedSeq*(aligned2[s2]);
		s2=s2+1;
		k=k+1;
	}

	gappedSeq * 0;


	if (refSeq2 && seqCounter == seqCount-1)
	{
		fscanf (fName, "Raw", soFar);
		fprintf (fName, CLEAR_FILE,">",seqName,"\n",gappedSeq,"\n",soFar);
		
	}
	else
	{
		fprintf (fName,"\n>",seqName,"\n",gappedSeq);
	}
}

