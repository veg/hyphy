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

SKIP_OMISSIONS = 1;

if (dataType)
{
	DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
}
else
{
	DataSetFilter filteredData = CreateFilter (ds,1,"","");
}

SKIP_OMISSIONS = 0;

if (dataType)
{
	DataSetFilter filteredData2 = CreateFilter (ds,3,"","",GeneticCodeExclusions);
}
else
{
	DataSetFilter filteredData2 = CreateFilter (ds,1,"","");
}

fprintf (stdout, "\nUpper log-likelihood bound on the data set:",LAST_FILE_PATH,
				 "\nSequences : ",filteredData.species,
				 "\nSites     : ",filteredData.sites,
				 "\nSite types: ",filteredData.unique_sites);
				 
if (filteredData2.sites!=filteredData.sites)
{
	fprintf (stdout, "\n\nWARNING: Data set contains deletions. The upper bound is NOT valid for this data file.\n\n");
}
upperLnLikBound = 0;

for (counter = 0; counter < filteredData.unique_sites; counter = counter+1)
{
	upperLnLikBound = upperLnLikBound + filteredData.site_freqs[counter] *
										Log (filteredData.site_freqs[counter]/filteredData.sites);
}

fprintf (stdout,"\nUpper bound on likelihood = ", upperLnLikBound,"\n");
