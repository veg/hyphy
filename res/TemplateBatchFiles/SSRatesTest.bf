/* This batch file will test whether elevated substitution rates are strand specific.
In some ssDNA viruses, the virus is in a ss state sufficiently long for it to be subject to 
oxidative deamination. Oxidative deamination results in elevated rates of 
C-T; G-A and G-T. C-T and G-A are complimentary mutations, but G-T has a C-A complement. A 
test of strand specific substitution rates therefore involves a likelihood ratio test of
a model in which rates are constrained to their complement (ie C-T := G-A) versus a model in
which all rates are independent.
*/

function printModelDescription ( gammaRate, freqsVector, strandGTR )
{
if ( strandGTR ) {
	sep =            "+---+-------------+-------------+-------------+-------------+\n";
	fprintf (stdout, "Rate matrix\n", 
					 sep,
					 "|   |      A      |      C      |      G      |      T      |\n",
					 sep,
					 "| A |      *      | ",Format(AC,11,5)," | ", Format(1.0,11,5), " | ", Format(AT,11,5), " |\n",
					 sep,
					 "| C | ",Format(GT,11,5)," |      *      | ", Format(CG,11,5), " | ", Format(CT,11,5), " |\n",
					 sep,
					 "| G | ",Format(CT,11,5), " | ",  Format(CG,11,5), " |      *      | ", Format(GT,11,5), " |\n",
					 sep,
					 "| T | ",Format(AT,11,5), " | ",  Format(AG,11,5), " | ",  Format(AC,11,5)," |      *      |\n",
					 sep);
						 
	}
	else {
		sep =            "+---+-------------+-------------+-------------+-------------+\n";
	fprintf (stdout, "Rate matrix\n", 
					 sep,
					 "|   |      A      |      C      |      G      |      T      |\n",
					 sep,
					 "| A |      *      | ",Format(AC,11,5)," | ", Format(1.0,11,5), " | ", Format(AT,11,5), " |\n",
					 sep,
					 "| C | ",Format(AC,11,5)," |      *      | ", Format(CG,11,5), " | ", Format(CT,11,5), " |\n",
					 sep,
					 "| G | ",Format(1,11,5), " | ",  Format(CG,11,5), " |      *      | ", Format(GT,11,5), " |\n",
					 sep,
					 "| T | ",Format(AT,11,5), " | ",  Format(CT,11,5), " | ",  Format(GT,11,5)," |      *      |\n",
					 sep);
	}
	fprintf (stdout, "\nBase frequencies",
					 "\nA = ", Format (freqsVector[0],10,4),
					 "\nC = ", Format (freqsVector[1],10,4),
					 "\nG = ", Format (freqsVector[2],10,4),
					 "\nT = ", Format (freqsVector[3],10,4), "\n");
					 
	if (gammaRate>1) /* gamma or HMM*/
	{
		fprintf (stdout, "\nRate variation");
		GetInformation 		(catInfo, c);
		for (idx = 0; idx < Columns (catInfo); idx = idx + 1)
		{
			fprintf (stdout, "\nRate ",Format (idx+1,5,0), " = ", Format (catInfo[0][idx], 10,5), " (weight = ", Format (catInfo[1][idx], 10,5), ")");
		}
		fprintf (stdout, "\n");
	}
	
	return 1;
}

SetDialogPrompt ("Select a nucleotide alignment");
DataSet 		ds = ReadDataFile(PROMPT_FOR_FILE);
fprintf			(stdout, "Read an alignment with ", ds.species, " sequences and ", ds.sites, " sites\n");
ExecuteAFile 	(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "queryTree.bf");

DataSetFilter filteredData = CreateFilter (ds,1);

fprintf ( stdout, "Select parameter settings for the complement constrained GRM model\n" );
incFileName = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR +  "TemplateModels" + DIRECTORY_SEPARATOR + "STGRM.mdl";
ExecuteCommands ( "#include  \"" + incFileName + "\";" );
Export ( Modelstring, USE_LAST_MODEL );
fprintf ( stdout, Modelstring, "\n" );
Tree STGRMTree = treeString;
LikelihoodFunction lf_STGRM = ( filteredData, STGRMTree );
Optimize ( res_STGRM, lf_STGRM );
totalBranchCount = TipCount(STGRMTree) + BranchCount ( STGRMTree );
ll_STGRM = res_STGRM [ 1 ][ 0 ];
aic_STGRM = -2*ll_STGRM + 2*(res_STGRM [ 1 ][ 1 ] + totalBranchCount + 3);

if ( modelType >= 1 ) {
	fprintf ( stdout, "\nstrand specific General reversible model optimised parameters\n" );
	fprintf ( stdout, "---------------------------------------------------------------\n" );
	printModelDescription ( modelType, vectorOfFrequencies, 1 );
	fprintf ( stdout, "\n" );
}

fprintf ( stdout, "Select parameter settings for the NRM model\n" );
incFileName = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR +  "TemplateModels" + DIRECTORY_SEPARATOR + "NRM-Freqs.mdl";
ExecuteCommands ( "#include  \"" + incFileName + "\";" );
Export ( Modelstring, USE_LAST_MODEL );
fprintf ( stdout, Modelstring, "\n" );
Tree NRMTree = treeString;
LikelihoodFunction lf_NRM = ( filteredData, NRMTree );
Optimize ( res_NRM, lf_NRM );
ll_NRM = res_NRM [ 1 ][ 0 ];
aic_NRM = -2*ll_NRM + 2*(res_NRM [ 1 ][ 1 ] + totalBranchCount + 3);

if ( modelType >= 1 ) {
	fprintf ( stdout, "\nstrand specific General reversible model optimised parameters\n" );
	fprintf ( stdout, "---------------------------------------------------------------\n" );
	printModelDescription ( modelType, vectorOfFrequencies, 1 );
	fprintf ( stdout, "\n" );
}

fprintf ( stdout, "Model comparison\n" );
fprintf ( stdout, "--------------------\n" );
fprintf ( stdout, "stGRM: logLk = ", ll_STGRM, "; AIC = ", aic_STGRM, "; Parameters = ", res_STGRM [ 1 ][ 1 ] + totalBranchCount + 3, "\n" );
fprintf ( stdout, "NRM: logLk = ", ll_NRM, "; AIC = ", aic_NRM, "; Parameters = ", res_NRM [ 1 ][ 1 ] + totalBranchCount + 3, "\n" );

fprintf ( stdout, "Likelihood ratio test statistic = ", 2*(ll_NRM-ll_STGRM), "; P = ", 1-CChi2 ( 2*(ll_NRM-ll_STGRM), res_NRM[1][1]-res_STGRM[1][1] ), "\n" ); 




