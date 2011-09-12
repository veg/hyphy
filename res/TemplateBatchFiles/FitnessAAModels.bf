aaCodes = {{"A",
"C",
"D",
"E",
"F",
"G",
"H",
"I",
"K",
"L",
"M",
"N",
"P",
"Q",
"R",
"S",
"T",
"V",
"W",
"Y"}};

aaNames = {{"Alanine",
"Cysteine",
"Aspartic Acid",
"Glutamic Acid",
"Phenylalanine",
"Glycine",
"Histidine",
"Isoleucine",
"Lysine",
"Leucine",
"Methionine",
"Asparagine",
"Proline",
"Glutamine",
"Arginine",
"Serine",
"Threonine",
"Valine",
"Tryptophan",
"Tyrosine"}};

function PopulateModelMatrix (ModelMatrixName&, EFV&, classIndex)
{
	ModelMatrixName = {20,20};
	EFV				= {20,1};
	
	commandString = "";
	
	for (rI = 1; rI < 20; rI = rI+1)
	{
		commandString = commandString + "global F_" + aaCodes[rI] + "_" + classIndex + 
						"=0; F_" + aaCodes[rI] + "_" + classIndex + 
						":>-1e10;\n";
	}
	
	ExecuteCommands (commandString);
	
	for (cI = 1; cI < 20; cI = cI+1)
	{
		commandString = "ModelMatrixName[0]["+cI+"]:=t*(1+(F_" + 
						 aaCodes[cI] + "_" + classIndex + "<0)*(Exp(F_"+
						 aaCodes[cI] + "_" + classIndex + ")-1));\n" + 
						 "ModelMatrixName["+cI+"][0]:=t*(1+(0>" + "F_" + 
						 aaCodes[cI] + "_" + classIndex + ")*(Exp(-F_"+
						 aaCodes[cI] + "_" + classIndex + ")-1));\n"; 
		ExecuteCommands	 (commandString);
	}

	for (rI = 1; rI < 19; rI = rI+1)
	{
		for (cI = rI+1; cI < 20; cI = cI+1)
		{
			commandString = "ModelMatrixName["+rI+"]["+cI+"]:=t*(1+(F_" + 
							 aaCodes[cI] + "_" + classIndex + "<" + "F_" + 
							 aaCodes[rI] + "_" + classIndex + ")*(Exp(F_"+
							 aaCodes[cI] + "_" + classIndex + "-F_"+
							 aaCodes[rI] + "_" + classIndex + ")-1));\n" + 
							 "ModelMatrixName["+cI+"]["+rI+"]:=t*(1+(F_" + 
							 aaCodes[rI] + "_" + classIndex + "<" + "F_" + 
							 aaCodes[cI] + "_" + classIndex + ")*(Exp(F_"+
							 aaCodes[rI] + "_" + classIndex + "-F_"+
							 aaCodes[cI] + "_" + classIndex + ")-1));\n"; 
			ExecuteCommands	 (commandString);
		}
	}
	
	commandString = "global EFV_Norm_"+classIndex+":=1";
	for (rI = 1; rI < 20; rI = rI+1)
	{
		commandString = commandString+"+Exp(F_"+aaCodes[rI]+"_"+classIndex+")";
	}
	ExecuteCommands (commandString+";");
	commandString = "EFV[0]:=1/EFV_Norm_"+classIndex+";";
	ExecuteCommands (commandString);
	for (rI = 1; rI < 20; rI = rI+1)
	{
		commandString = "EFV["+rI+"]:=Exp(F_"+aaCodes[rI]+"_"+classIndex+")/EFV_Norm_"+classIndex+";";
		ExecuteCommands (commandString);
	}
	return 0;
}


fprintf (stdout, "\n\n***** RUNNING MULTIPLE FITNESS CLASS AMINOACID MODELS*****\n\nModel Credit: Matthew W. Dimmic, David P. Mindell and Richard A. Goldstein\n\n");

fitnessClassNumber = 0;
while (fitnessClassNumber<1)
{
	fprintf (stdout, "\n\nHow many fitness classes (>=1) do you wish to spawn?");
	fscanf  (stdin, "Number", fitnessClassNumber);
}

/*testMatrix = 0;
testEFV	   = 0;
dummy = PopulateModelMatrix ("testMatrix","testEFV",0);*/

SetDialogPrompt ("Please specify a nucleotide or amino-acid data file:");

DataSet ds 	= ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter (ds,1);
fprintf (stdout, "DATASET:", LAST_FILE_PATH, "\n\n", fitnessClassNumber, " FITNESS CLASSES\n\n");

if (IS_TREE_PRESENT_IN_DATA)
{
	fprintf (stdout, "\n\nA tree was found in the data file:\n",DATAFILE_TREE,"\n\nWould you like to use it:(Y/N)?");
	fscanf (stdin, "String", response);
	if ((response=="n")||(response=="N"))
	{
		IS_TREE_PRESENT_IN_DATA = 0;
	}
	else
	{
		if (_DO_TREE_REBALANCE_)
		{
			treeString = RerootTree (DATAFILE_TREE,0);
		}
		else
		{
			treeString = DATAFILE_TREE;
		}
		
		IS_TREE_PRESENT_IN_DATA = 1;
	}
	fprintf (stdout, "\n\n");

}

if (!IS_TREE_PRESENT_IN_DATA)
{
	SetDialogPrompt ("Please select a tree file for the data:");

	fscanf (PROMPT_FOR_FILE, "String", treeString);
	
	while ((Abs(treeString)==0)||(treeString[0]!="("))
	{
		if (END_OF_FILE)
		{
			fprintf (stdout, "\nThis doesn't seem to be a valid Newick string file.\n").
			return;
		}
		fscanf (LAST_FILE_PATH, "String", treeString);
	}
	
	if (_DO_TREE_REBALANCE_)
	{
		treeString = RerootTree (treeString,0);
	}
}

if (fitnessClassNumber==1)
{
	fitnessRates = 0;
	fitnessEFV_0	 = 0;
	mi = PopulateModelMatrix ("fitnessRates","fitnessEFV_0", 0);
	Model fitnessModel = (fitnessRates, fitnessEFV_0, 0);
	Tree givenTree = treeString;
	LikelihoodFunction lf = (filteredData,givenTree);
}
else
{
	freqStrMx    = {fitnessClassNumber,1};
	freqStrMx[0] = "PS_1";

	for (mi=1; mi<fitnessClassNumber; mi=mi+1)
	{
		ExecuteCommands ("global PS_"+mi+" = 1/"+ ((fitnessClassNumber+1)-mi) + ";\nPS_"+mi+":<1;\n");
	}

	for (mi=1; mi<fitnessClassNumber-1; mi=mi+1)
	{
		freqStrMx[mi] = "";
		for (mi2=1;mi2<=mi;mi2=mi2+1)
		{
			freqStrMx[mi] = freqStrMx[mi]+"(1-PS_"+mi2+")";		
		}
		freqStrMx[mi] = freqStrMx[mi]+"PS_"+(mi+1);	
	}	

	freqStrMx[mi] = "";
	for (mi2=1;mi2<mi;mi2=mi2+1)
	{
		freqStrMx[mi] = freqStrMx[mi]+"(1-PS_"+mi2+")";		
	}
	freqStrMx[mi] = freqStrMx[mi]+"(1-PS_"+mi+")";
	
	lfString 	   = "LikelihoodFunction lf = (";	
	templateString = "Log(";
	
	for (fC = 0; fC < fitnessClassNumber; fC = fC+1)
	{	
		ExecuteCommands ("fitnessRates_"+fC+"=0;fintessEFV_"+fC+"=0;dummy=PopulateModelMatrix(\"fitnessRates_"+
						  fC+"\",\"fitnessEFV_"+fC+"\",fC);Model fitnessModel_"+fC+"=(fitnessRates_"+fC+",fitnessEFV_"+fC+",0);Tree tree_"+
						  fC+"="+treeString+";");
		if (fC)
		{
			lfString = lfString+",";
			templateString = templateString + "+";
		}
		
		lfString = lfString+"filteredData,tree_"+fC;	
		templateString = templateString	+ "SITE_LIKELIHOOD["+fC+"]*"+freqStrMx[fC];
	}
	templateString = templateString+")";
	ExecuteCommands (lfString+",templateString);");
}

ChoiceList  (rI,"Optimize or Restore",1,NO_SKIP,
			 "Optimize","Obtain MLEs of model parameters.",
			 "Restore","Reload MLEs from a previously saved file.");

if (rI<0)
{
	return 0;
}

if (rI==0)
{
	SetDialogPrompt ("Save MLEs to:");
	fprintf (PROMPT_FOR_FILE, CLEAR_FILE);
	Optimize (res,lf);
	mi				   = LIKELIHOOD_FUNCTION_OUTPUT;
	fprintf (stdout, lf);
	LIKELIHOOD_FUNCTION_OUTPUT = 4;
	fprintf (LAST_FILE_PATH, lf);
	LIKELIHOOD_FUNCTION_OUTPUT = mi;
}
else
{
	SetDialogPrompt ("Load MLEs from:");
	fscanf			(PROMPT_FOR_FILE,"String", dummy);
	ExecuteCommands ("#include\""+LAST_FILE_PATH+"\";");
}

estimatedFitnessParameters = {20,fitnessClassNumber};
estimatedEFVs			   = {20,fitnessClassNumber};

for (fC = 0; fC < fitnessClassNumber; fC = fC+1)
{
	for (mi=1; mi<20; mi=mi+1)
	{
		ExecuteCommands ("estimatedFitnessParameters[mi][fC]=F_"+aaCodes[mi]+"_"+fC+";");
	}
	for (mi=0; mi<20; mi=mi+1)
	{
		ExecuteCommands ("estimatedEFVs[mi][fC]=fitnessEFV_"+fC+"[mi];");
	}
}

if (fitnessClassNumber>1) 
/* compute marginals */
{
	fitnessMarginals = {filteredData.sites,fitnessClassNumber+1};
	ConstructCategoryMatrix(flatMarginals, lf, COMPLETE);
	mi = 0;
	for (fC = 0; fC < fitnessClassNumber; fC = fC+1)
	{
		for (rI = 0; rI < filteredData.sites; rI = rI + 1)
		{
			fitnessMarginals [rI][fC] = flatMarginals[mi];
			mi = mi+1;
		}
	}
	classPriors = {fitnessClassNumber,1};
	for (fC = 0; fC < fitnessClassNumber; fC = fC+1)
	{
		ExecuteCommands ("classPriors[fC]="+freqStrMx[fC]+";");	
	}
	for (rI = 0; rI < filteredData.sites; rI = rI + 1)
	{
		maxTerm  = 0;
		maxIndex = 0;
		cI		 = 0;
		
		for (fC = 0; fC < fitnessClassNumber; fC = fC+1)
		{
			flatMarginals = classPriors[fC]*fitnessMarginals[rI][fC];
			cI = cI + flatMarginals;
			if (flatMarginals>maxTerm)
			{
				maxTerm  = flatMarginals;
				maxIndex = fC;
			}
		}
		fitnessMarginals[rI][fitnessClassNumber] = maxIndex;
		for (fC = 0; fC < fitnessClassNumber; fC = fC+1)
		{
			fitnessMarginals[rI][fC] =  classPriors[fC]*fitnessMarginals[rI][fC]/cI;
		}
	}	
}


labelMatrix  = {1,fitnessClassNumber+1};
labelMatrix2 = {1,fitnessClassNumber+1};
labelMatrix3 = {{"Class Prior"}};

seriesString = "";

for (fC = 0; fC < fitnessClassNumber; fC = fC+1)
{
	labelMatrix[fC] = "Class "+fC;
	labelMatrix2[fC] = "Class "+fC;
	if (fC)
	{
		seriesString = seriesString + ";";
	}
	seriesString = seriesString + labelMatrix[fC];
}

aaString = "Aminoacid";

for (fC = 0; fC < 20; fC = fC+1)
{
	aaString = aaString + ";" + aaNames[fC];
}

labelMatrix[fitnessClassNumber] = aaString;

labelMatrix2[fitnessClassNumber] = "Model Assignment";

OpenWindow (CHARTWINDOW,{{"Estimated Fitnesses"}
						   {"labelMatrix"},
						   {"estimatedFitnessParameters"},
						   {"Line Plot"},
						   {"Index"},
						   {seriesString},
						   {"AA Index"},
						   {""},
						   {"Relative Fitness"},
						   {"3"}},
						   "(SCREEN_WIDTH-30)/3;(SCREEN_HEIGHT-50)/2;10;SCREEN_HEIGHT/2+25");

OpenWindow (CHARTWINDOW,{{"Estimated EFV"}
						   {"labelMatrix"},
						   {"estimatedEFVs"},
						   {"Line Plot"},
						   {"Index"},
						   {seriesString},
						   {"AA Index"},
						   {""},
						   {"Estimated Freq"},
						   {"3"}},
						   "(SCREEN_WIDTH-30)/3;(SCREEN_HEIGHT-50)/2;10+(SCREEN_WIDTH-30)/3;SCREEN_HEIGHT/2+25");
						   
if (fitnessClassNumber>1)
{
	OpenWindow (CHARTWINDOW,{{"Class Priors"}
							   {"labelMatrix3"},
							   {"classPriors"},
							   {"Pie Chart"},
							   {"Class Prior"},
							   {"None"},
							   {""},
							   {""},
							   {"Class Prior"},
							   {"3"}},
							   "(SCREEN_WIDTH-30)/3;(SCREEN_HEIGHT-50)/2;12+2*(SCREEN_WIDTH-30)/3;SCREEN_HEIGHT/2+25");


	OpenWindow (CHARTWINDOW,{{"Class Posteriors"}
							   {"labelMatrix2"},
							   {"fitnessMarginals"},
							   {"Stacked Bars"},
							   {"Index"},
							   {seriesString},
							   {"AA Index"},
							   {""},
							   {"Class Posterior"},
							   {"3"}},
							   "SCREEN_WIDTH-30;(SCREEN_HEIGHT-50)/2;15;25");
}


