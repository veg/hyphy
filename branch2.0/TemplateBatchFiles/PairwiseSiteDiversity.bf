distanceChoice = 0;

ChoiceList (dataType,"Data type",1,SKIP_NONE,"Nucleotide/Protein","Nucleotide or amino-acid (protein).",
				     "Codon","Codon (several available genetic codes).");

if (dataType<0) 
{
	return;
}

if (dataType)
{
	NICETY_LEVEL = 3;
	SetDialogPrompt ("Please choose a codon data file:");
	#include "TemplateModels/chooseGeneticCode.def";
}
else
{
	SetDialogPrompt ("Please choose a nucleotide or amino-acid data file:");
}


DataSet ds = ReadDataFile (PROMPT_FOR_FILE);

if (dataType)
{
	DataSetFilter fullData 	   = CreateFilter (ds,3,"","",GeneticCodeExclusions);
	DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
}
else
{
	DataSetFilter fullData = CreateFilter (ds,1);
	DataSetFilter filteredData = CreateFilter (ds,1);
}

sampleCount = 0;

while (sampleCount < 2)
{
	fprintf (stdout,"\nHow many samples are there in the data file(>=2)?");
	fscanf  (stdin, "Number", sampleCount);
}

sampleSize = {sampleCount, 1};

samplesDone  = 0;
totalSamples = 0;

while (samplesDone < sampleCount-1)
{
	thisSample = 0;
	while (thisSample <= 0)
	{
		fprintf (stdout,"\nHow many sequences are in sample ", samplesDone+1, " (between 1 and ", ds.species-totalSamples,")?");
		fscanf  (stdin, "Number", thisSample);
	}
	sampleSize [samplesDone] = thisSample;
	totalSamples = totalSamples + thisSample;
	samplesDone= samplesDone + 1;
}

sampleSize [samplesDone] = ds.species-totalSamples;

fprintf (stdout,"\n\nUsing these sample sizes:\n");

for (samplesDone = 0; samplesDone < sampleCount; samplesDone = samplesDone + 1)
{
	fprintf (stdout, "Sample ", samplesDone + 1, " : ", sampleSize [samplesDone], "\n");
}

DISTANCE_PROMPTS = 1;

withinSampleDiversity = {sampleCount, fullData.sites};

#include "chooseDistanceFormula.def";

fprintf (stdout, "\n\nComputing within sample site-by-site diversities\n\n");

firstSequence = 0;

dummy = InitializeDistances (0);

for (samplesDone = 0; samplesDone < sampleCount; samplesDone = samplesDone + 1)
{
	fprintf (stdout, "\nProcessing sample ", samplesDone+1, "\n");

	seqCount     = sampleSize[samplesDone];
	lastSequence = firstSequence + seqCount;
	
	for (sitesDone = 0; sitesDone < fullData.sites; sitesDone = sitesDone + 1)
	{
		meanDistance = 0;
		seqFilter  = ""+(firstSequence)+"-"+(lastSequence-1);
		
		if (dataType)
		{
			siteFilter = ""+(3*sitesDone)+"-"+(3*sitesDone+2);
			DataSetFilter filteredData = CreateFilter (ds,3,siteFilter,seqFilter,GeneticCodeExclusions);
		}
		else
		{
			siteFilter = ""+sitesDone;
			DataSetFilter filteredData = CreateFilter (ds,1,siteFilter,seqFilter);
		}
		
		for (sequence1 = 0; sequence1 < seqCount-1; sequence1 = sequence1 + 1)
		{
			for (sequence2 = sequence1 + 1; sequence2 < seqCount; sequence2 = sequence2 + 1)
			{
				meanDistance = meanDistance + ComputeDistanceFormula (sequence1, sequence2);
			}
		}
		
		withinSampleDiversity [samplesDone][sitesDone] = 2*meanDistance/seqCount/(seqCount-1);
		/*fprintf (stdout, "\nSite ", Format (sitesDone+1,5,0), ":", Format (withinSampleDiversity [samplesDone][sitesDone],8,4));*/
	}
	
	firstSequence = lastSequence;
}


fprintf (stdout, "\n\nComputing between sample site-by-site diversities\n\n");

firstSequence = 0;

betweenSampleDiversity = {sampleCount*(sampleCount-1)/2, fullData.sites};

pComp = 0;

for (samplesDone = 0; samplesDone < sampleCount-1; samplesDone = samplesDone + 1)
{
	firstSequence2  = firstSequence;
	seqCount 	    = sampleSize[samplesDone];
	lastSequence    = firstSequence + seqCount;
	firstSequence2  = lastSequence;

	for (samplesDone2 = samplesDone+1; samplesDone2 < sampleCount; samplesDone2 = samplesDone2 + 1)
	{
		seqCount2 	   = sampleSize[samplesDone2];
		lastSequence2  = firstSequence2 + seqCount2;
		
		fprintf (stdout, "\nProcessing samples ", samplesDone+1, " and ", samplesDone2 +1 ,"\n");

		lastSequence = firstSequence + seqCount;
		
		for (sitesDone = 0; sitesDone < fullData.sites; sitesDone = sitesDone + 1)
		{
			meanDistance = 0;
			seqFilter  = ""+(firstSequence)+"-"+(lastSequence-1)+","+(firstSequence2)+"-"+(lastSequence2-1);
			
			if (dataType)
			{
				siteFilter = ""+(3*sitesDone)+"-"+(3*sitesDone+2);
				DataSetFilter filteredData = CreateFilter (ds,3,siteFilter,seqFilter,GeneticCodeExclusions);
			}
			else
			{
				siteFilter = ""+sitesDone;
				DataSetFilter filteredData = CreateFilter (ds,1,siteFilter,seqFilter);
			}
			
			for (sequence1 = 0; sequence1 < seqCount; sequence1 = sequence1 + 1)
			{
				for (sequence2 = seqCount; sequence2 < seqCount+seqCount2; sequence2 = sequence2 + 1)
				{
					meanDistance = meanDistance + ComputeDistanceFormula (sequence1, sequence2);
				}
			}
			
			betweenSampleDiversity [pComp][sitesDone] = meanDistance/seqCount/seqCount2;
		}
		firstSequence2 = lastSequence2;
		pComp = pComp+1;
	}
	firstSequence = lastSequence;
}

fprintf (stdout, "\n\nGenerating Results\n\n");

pComp = 0;

outString = "";
p = outString*256000;
p = outString*"Sample1\tSample2\tSize1\tSize2\tCodon\tDiv1\tDiv2\tDivBetween\n";

for (samplesDone = 0; samplesDone < sampleCount-1; samplesDone = samplesDone + 1)
{
	for (samplesDone2 = samplesDone+1; samplesDone2 < sampleCount; samplesDone2 = samplesDone2 + 1)	
	{
		for (codonCount=0; codonCount < fullData.sites; codonCount = codonCount + 1)
		{
			p=outString*(""+(samplesDone+1)+"\t"+(samplesDone2+1)+"\t"+sampleSize[samplesDone]+"\t"+sampleSize[samplesDone2]+"\t"+(codonCount+1)+"\t"+
							 withinSampleDiversity[samplesDone][codonCount]+"\t"+withinSampleDiversity[samplesDone2][codonCount]+"\t"+betweenSampleDiversity[pComp][codonCount]+"\n");
		}
		pComp = pComp+1;
	}
	p=outString*"\n";
}

p = outString*0;
SetDialogPrompt ("Write results to:");
fprintf (PROMPT_FOR_FILE,CLEAR_FILE,outString);
