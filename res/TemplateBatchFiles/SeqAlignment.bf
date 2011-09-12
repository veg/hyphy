RequireVersion ("2.0020100903");
ExecuteAFile   ("Utility/GrabBag.bf");

#include		"SeqAlignShared.ibf";

ASSERTION_BEHAVIOR = 1;

/* END ALIGNMENT SETTINGS */

SetDialogPrompt 		("Sequence File:");
doLongestSequence = 	(refSeq==1);

if (refSeq > 1)
{
	DataSet        unal2 			= ReadDataFile 	(PROMPT_FOR_FILE);
	refSeq = ">" + predefSeqNames [refSeq][0] + "\n" + RefSeqs[refSeq-2];
	DataSet		   refD = ReadFromString (refSeq);
	DataSet        unal = Combine (refD,unal2);
}
else
{
	ChoiceList (refSeq2,"Insert a coordinate reference sequence?",1,SKIP_NONE,predefSeqNames2);
	assert (refSeq2 >= 0, "Reference sequence selection"); 
	if (refSeq2)
	{
		DataSet        unal2 			= ReadDataFile 	(PROMPT_FOR_FILE);
		refSeq = ">" + predefSeqNames [refSeq2][0] + "\n" + RefSeqs[refSeq2-1];
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




ChoiceList (skipOutliers,"Include outliers?",		1,SKIP_NONE,"No","Skip sequences with unusually poor alignment scores","Yes","Include all alignable sequences");
assert (skipOutliers >= 0, "Outlier handling selection"); 

ChoiceList (doRC,		 "Try reverse complements?",1,SKIP_NONE,"No","Only align in the provided direction","Yes","Try both each sequences and its reverse complement.");
assert	 (doRC >= 0, "Reverse complement selection"); 

// build codon translation table
ExecuteAFile       ("TemplateModels/chooseGeneticCode.def");
codonToAAMap		= defineCodonToAA();

_handleAlignment ("unal",1);
