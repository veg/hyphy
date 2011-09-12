/* will do a test of a standard GRM against a strand complimentary GRM using AIC 
strand complimentary GRM has AC := TG not AC := CA

Wayne Delport
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

fprintf ( stdout, "Select parameter settings for GRM model:\n" );
incFileName = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR +  "TemplateModels" + DIRECTORY_SEPARATOR + "GRM.mdl";
ExecuteCommands ( "#include  \"" + incFileName + "\";" );
Export ( Modelstring, USE_LAST_MODEL );
fprintf ( stdout, Modelstring, "\n" );
Tree GRMTree = treeString;
LikelihoodFunction lf_GRM = ( filteredData, GRMTree );
Optimize ( res_GRM, lf_GRM );
totalBranchCount = TipCount(GRMTree) + BranchCount ( GRMTree );
ll_GRM = res_GRM [ 1 ][ 0 ];
aic_GRM = -2*ll_GRM + 2*(res_GRM [ 1 ][ 1 ] + totalBranchCount + 3);

if ( modelType >= 1 ) {
	fprintf ( stdout, "\nGeneral reversible model optimised parameters\n" );
	fprintf ( stdout,   "---------------------------------------------\n" );
	printModelDescription ( modelType, vectorOfFrequencies, 0 );
	fprintf ( stdout, "\n" );
}

fprintf ( stdout, "Select parameter settings for stGRM model:\n" );
incFileName = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR +  "TemplateModels" + DIRECTORY_SEPARATOR + "STGRM.mdl";
ExecuteCommands ( "#include  \"" + incFileName + "\";" );
Export ( Modelstring, USE_LAST_MODEL );
fprintf ( stdout, Modelstring, "\n" );
Tree STGRMTree = treeString;
LikelihoodFunction lf_STGRM = ( filteredData, STGRMTree );
Optimize ( res_STGRM, lf_STGRM );
ll_STGRM = res_STGRM [ 1 ][ 0 ];
aic_STGRM = -2*ll_STGRM + 2*(res_STGRM [ 1 ][ 1 ] + totalBranchCount + 3);

if ( modelType >= 1 ) {
	fprintf ( stdout, "\nstrand specific General reversible model optimised parameters\n" );
	fprintf ( stdout, "---------------------------------------------------------------\n" );
	printModelDescription ( modelType, vectorOfFrequencies, 1 );
	fprintf ( stdout, "\n" );
}

fprintf ( stdout, "AIC model comparison\n" );
fprintf ( stdout, "--------------------\n" );
fprintf ( stdout, "GRM: logLk = ", ll_GRM, "; AIC = ", aic_GRM, "; Parameters = ", res_GRM [ 1 ][ 1 ] + totalBranchCount + 3, "\n" );
fprintf ( stdout, "stGRM: logLk = ", ll_STGRM, "; AIC = ", aic_STGRM, "; Parameters = ", res_STGRM [ 1 ][ 1 ] + totalBranchCount + 3, "\n" );

