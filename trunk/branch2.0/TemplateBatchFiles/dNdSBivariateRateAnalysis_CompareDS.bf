RequireVersion ("0.9920060817");

ChoiceList (runType,"Run Type",1,SKIP_NONE,
			"New run","Start a new run",
			"Continued run","Start with an initial point output by a Gateux derivative approximator."
		    );
		    
if (runType < 0)
{
	return 0;
}

randomizeInitValues  = 0;
ModelMatrixDimension = 0;
/*---------------------------------------------------------------------------------------------------------------------------------------------------*/

function BuildCodonFrequencies (obsF)
{
	PIStop = 1.0;
	result = {ModelMatrixDimension,1};
	hshift = 0;

	for (h=0; h<64; h=h+1)
	{
		first = h$16;
		second = h%16$4;
		third = h%4;
		if (_Genetic_Code[h]==10) 
		{
			hshift = hshift+1;
			PIStop = PIStop-obsF[first][0]*obsF[second][1]*obsF[third][2];
			continue; 
		}
		result[h-hshift][0]=obsF[first][0]*obsF[second][1]*obsF[third][2];
	}
	return result*(1.0/PIStop);
}

/*---------------------------------------------------------------------------------------------------------------------------------------------------*/

function ReportDistributionString (rc)
{
	distroString = "";
	distroString * 1024;
	
	reportMx = {rc,4};
	
	for (mi=0; mi<rc; mi=mi+1)
	{
		ExecuteCommands ("reportMx[mi][0]=S_"+mi+"/c_scale;");
		ExecuteCommands ("reportMx[mi][1]=NS_"+mi+"/c_scale;");
		reportMx[mi][2] = reportMx[mi][1]/reportMx[mi][0];
		ExecuteCommands ("reportMx[mi][3]="+freqStrMx[mi]+";");
		distroString * ("Class "+(mi+1)+"\n\tdS    = "+Format(reportMx[mi][0],10,3)
									   +"\n\tdN    = "+Format(reportMx[mi][1],10,3)
									   +"\n\tdN/dS = "+Format(reportMx[mi][2],10,3)
									   +"\n\tProb  = "+Format(reportMx[mi][3],10,3) + "\n");
	}								   
	
	distroString * 0;
	return distroString;
}


	
totalCodonCount1 = 0;
totalCodonCount1  = 0;
totalCodonCount2 = 0;
totalCharCount2  = 0;

/*---------------------------------------------------------------------------------------------------------------------------------------------------*/
if (runType == 0)
{
	ChoiceList (branchLengths,"Branch Lengths",1,SKIP_NONE,
				"Codon Model","Jointly optimize rate parameters and branch lengths (slow and thorough)",
				"Nucleotide Model","Estimate branch lengths once, using an appropriate nucleotide model (quick and dirty)."
			    );

	if (branchLengths<0)
	{
		return;
	}

	fileCount	= 0;
	while (fileCount < 2)
	{
		fprintf (stdout, "How many datafiles are to be analyzed (>=2):?");
		fscanf  (stdin, "Number", fileCount);
		fileCount = fileCount $ 1;
	}
	
	fileCount1 = 0;
	while (fileCount1 < 1 && fileCount1 > fileCount-1)
	{
		fprintf (stdout, "How many parts are there for the 1st data set?[1-",fileCount-1,":");
		fscanf  (stdin, "Number", fileCount1);
		fileCount1 = fileCount1 $ 1;
	}
	
	fileCount2 = fileCount - fileCount1;
	fprintf (stdout, fileCount1, " parts in dataset 1 and ", fileCount2, " - in dataset2\n");
	
	#include "TemplateModels/chooseGeneticCode.def";
	
	for (fileID = 1; fileID <= fileCount; fileID = fileID + 1)
	{
		fprintf (stdout, "\nReading input file ", fileID, "/", fileCount, "\n");
		SetDialogPrompt ("Please specify codon data #" + fileID + ":");

		ExecuteCommands (
		"DataSet				ds_" +fileID + " = ReadDataFile (PROMPT_FOR_FILE);" + 
		"fprintf (stdout,\"\\n______________READ THE FOLLOWING DATA______________\\n\",ds_" + fileID + ");"+
		"DataSetFilter filteredData_" + fileID + " = CreateFilter (ds_"+fileID +",3,\"\",\"\",GeneticCodeExclusions);");
		#include								   "queryTree.bf";
		if (fileID <= fileCount1)
		{
			ExecuteCommands ("treeString_" + fileID + " = treeString;totalCodonCount1=totalCodonCount1+filteredData_"+fileID+".sites;totalCodonCount1=totalCodonCount1+filteredData_"+fileID+".sites*filteredData_"+fileID+".species;");
		}
		else
		{
			ExecuteCommands ("treeString_" + fileID + " = treeString;totalCodonCount2=totalCodonCount2+filteredData_"+fileID+".sites;totalCodonCount2=totalCodonCount2+filteredData_"+fileID+".sites*filteredData_"+fileID+".species;");
		
		}
	}
	

	observedFreq       = {4,3};
	observedFreqSingle = {4,1};
	observedFreq2      = {4,3};
	observedFreqSingle2 = {4,1};
	for (fileID = 1; fileID <= fileCount; fileID = fileID + 1)
	{
		ExecuteCommands 	  ("HarvestFrequencies (tp, filteredData_"+fileID+",3,1,1);HarvestFrequencies (ts, filteredData_"+fileID+",1,1,1);cfs = filteredData_"+fileID+".sites;");
		if (fileID <= fileCount1)
		{
			observedFreq1 		 = observedFreq1 	    + tp*(cfs/totalCodonCount1);
			observedFreqSingle1  = observedFreqSingle1  + ts*(cfs/totalCodonCount1);
		}
		else
		{
			observedFreq2 		 = observedFreq2 	    + tp*(cfs/totalCodonCount2);
			observedFreqSingle2  = observedFreqSingle2  + ts*(cfs/totalCodonCount2);
		}
	}
	
	done = 0;
	while (!done)
	{
		fprintf (stdout,"\nPlease enter a 6 character model designation (e.g:010010 defines HKY85) for data set 1:");
		fscanf  (stdin,"String", modelDesc1);
		if (Abs(modelDesc1)==6)
		{	
			done = 1;
		}
	}			

	ExecuteAFile ("2RatesAnalyses/MG94xREVxBivariate.mdl");

	ModelTitle     = "MG94x"+modelDesc1[0];
	rateBiasTerms  = {{"AC","1","AT","CG","CT","GT"}};
	rateBiasTerms2 = {{"AC2","1","AT2","CG2","CT2","GT2"}};
	paramCount	   = 0;

	modelConstraintString = "";

	for (customLoopCounter2=1; customLoopCounter2<6; customLoopCounter2=customLoopCounter2+1)
	{
		for (customLoopCounter=0; customLoopCounter<customLoopCounter2; customLoopCounter=customLoopCounter+1)
		{
			if (modelDesc[customLoopCounter2]==modelDesc[customLoopCounter])
			{
				ModelTitle  = ModelTitle+modelDesc1[customLoopCounter2];	
				if (rateBiasTerms[customLoopCounter2] == "1")
				{
					modelConstraintString = modelConstraintString + rateBiasTerms[customLoopCounter]+":="+rateBiasTerms[customLoopCounter2]+";";
				}
				else
				{
					modelConstraintString = modelConstraintString + rateBiasTerms[customLoopCounter2]+":="+rateBiasTerms[customLoopCounter]+";";			
				}
				break;
			}
		}
		if (customLoopCounter==customLoopCounter2)
		{
			ModelTitle = ModelTitle+modelDesc1[customLoopCounter2];	
		}
	}	
		
	if (Abs(modelConstraintString))
	{
		ExecuteCommands (modelConstraintString);
	}
	global     AC2 = 1;
	global 	   AT2 = 1;
	global     CG2 = 1;
	global	   CT2 = 1;
	global     GT2 = 1;
	
	done = 0;
	while (!done)
	{
		fprintf (stdout,"\nPlease enter a 6 character model designation (e.g:010010 defines HKY85) for data set 2:");
		fscanf  (stdin,"String", modelDesc2);
		if (Abs(modelDesc2)==6)
		{	
			done = 1;
		}
	}			

	ModelTitle2 = "MG94x"+modelDesc2[0];
	modelConstraintString = "";
	for (customLoopCounter2=1; customLoopCounter2<6; customLoopCounter2=customLoopCounter2+1)
	{
		for (customLoopCounter=0; customLoopCounter<customLoopCounter2; customLoopCounter=customLoopCounter+1)
		{
			if (modelDesc2[customLoopCounter2]==modelDesc2[customLoopCounter])
			{
				ModelTitle2  = ModelTitle2+modelDesc2[customLoopCounter2];	
				if (rateBiasTerms[customLoopCounter2] == "1")
				{
					modelConstraintString = modelConstraintString + rateBiasTerms2[customLoopCounter]+":="+rateBiasTerms2[customLoopCounter2]+";";
				}
				else
				{
					modelConstraintString = modelConstraintString + rateBiasTerms2[customLoopCounter2]+":="+rateBiasTerms2[customLoopCounter]+";";			
				}
				break;
			}
		}
		if (customLoopCounter==customLoopCounter2)
		{
			ModelTitle2 = ModelTitle2+modelDesc2[customLoopCounter2];	
		}
	}	

	mcs2 = modelConstraintString;
	if (Abs(mcs2))
	{
		ExecuteCommands (modelConstraintString);
	}

	ChoiceList (randomizeInitValues, "Initial Value Options",1,SKIP_NONE,
				"Default",	 "Use default inital values for rate distribution parameters.",
				"Randomized",	 "Select initial values for rate distribution parameters at random.",
				"User Values",	 "Set initial values for rate distribution parameters and nucleotide biases");


	if (randomizeInitValues < 0)
	{
		return;
	}

	resp  = 0;

	while (resp<1)
	{
		fprintf (stdout,"Number of rate classes:");
		fscanf  (stdin, "Number", resp);
	}

	if (randomizeInitValues < 2)
	{
		ChoiceList (stratDNDS, "dN/dS Stratification",1,SKIP_NONE,
					"Unconstrained",	 "Do not place constraints on dN/dS classes.",
					"Split",	 "Assign all dN/dS classes to negative/neutral/positive clases.");
		if (stratDNDS<0)
		{
			return 0;
		}
	}
	else
	{
		stratDNDS = 0;
	}

	if (stratDNDS)
	{
		respM = -1;
		respN = 0;
		respP = 0;

		while (respM<0 || respM > resp)
		{
			fprintf (stdout,"Number of negative rate classes [0-",resp,"]:");
			fscanf  (stdin, "Number", respM);
		}

		if (respM < resp)
		{
			respN = -1;
			while (respN<0 || respN > resp-respM)
			{
				fprintf (stdout,"Number of neutral rate classes [0-",resp-respM,"]:");
				fscanf  (stdin, "Number", respN);
			}
		}

		respP = resp-respN-respM;

		fprintf (stdout, "\nUsing\n\t", respM, " negatively selected classes\n\t", respN, " neutrally evolving classes\n\t", 
							respP, " positively selected classes\n\n");
	}

	categDef1 = "";
	categDef1 * 1024;
	
	lfDef	  = {};
	
	for (fileID = 0; fileID < fileCount; fileID = fileID + 1)
	{
		lfDef[fileID] = "";
		lfDef[fileID] * 1024;
		lfDef[fileID] * "Log(";
	}
	

	global		S_0  = 0.5;
	S_0:>0.0000001;
	global		NS_0 = 0.1;

	for (mi=1; mi<resp; mi=mi+1)
	{
		categDef1 * ("global S_"+mi+"=0.5;S_"+mi+":>0.0000001;\nglobal NS_"+mi+";\n");
		if (randomizeInitValues)
		{
			categDef1*("global P_"+mi+" = Random(0.05,0.95);\nP_"+mi+":<1;\n");	
			categDef1*("global S_"+mi+" = Random(0.05,1);");	
		}
		else
		{
			categDef1*("global P_"+mi+" = 1/"+(resp+1-mi)+";\nP_"+mi+":<1;\n");
		}
	}

	freqStrMx    = {resp,1};
	if (resp>1)
	{
		freqStrMx[0] = "P_1";

		for (mi=1; mi<resp-1; mi=mi+1)
		{
			freqStrMx[mi] = "";
			for (mi2=1;mi2<=mi;mi2=mi2+1)
			{
				freqStrMx[mi] = freqStrMx[mi]+"(1-P_"+mi2+")";		
			}
			freqStrMx[mi] = freqStrMx[mi]+"P_"+(mi+1);	
		}	
		freqStrMx[mi] = "";

		for (mi2=1;mi2<mi;mi2=mi2+1)
		{
			freqStrMx[mi] = freqStrMx[mi]+"(1-P_"+mi2+")";		
		}
		freqStrMx[mi] = freqStrMx[mi]+"(1-P_"+mi+")";	
	}
	else
	{
		freqStrMx[0] = "1";
	}

	categDef1*( "\n\nglobal c_scale:=S_0*" + freqStrMx[0]);

	for (mi=1; mi<resp; mi=mi+1)
	{
		categDef1*( "+S_"+mi+"*" + freqStrMx[mi]);
	}

	for (mi=0; mi<resp; mi=mi+1)
	{
		for (fileID = 0; fileID < fileCount; fileID = fileID + 1)
		{
			if (mi)
			{
				lfDef[fileID] * "+";
			}
			lfDef[fileID]*(freqStrMx[mi]+"*SITE_LIKELIHOOD["+(fileID*resp+mi)+"]");
		}
			
	}

	categDef1 * ";";

	for (mi=0; mi<resp; mi=mi+1)
	{
		categDef1 * ("\nglobal R_"+mi+"=1;NS_"+mi+":=R_"+mi+"*S_"+mi+";\n");
	}
	
	if (randomizeInitValues == 2)
	{
		enteredValues = {resp,3};
		
		fprintf (stdout, "AC=");
		fscanf	(stdin,"Number",thisAlpha);
		categDef1 * ("AC:="+thisAlpha+";");
		fprintf (stdout, "AT=");
		fscanf	(stdin,"Number",thisAlpha);
		categDef1 * ("AT:="+thisAlpha+";");
		fprintf (stdout, "CG=");
		fscanf	(stdin,"Number",thisAlpha);
		categDef1 * ("CG:="+thisAlpha+";");
		fprintf (stdout, "CT=");
		fscanf	(stdin,"Number",thisAlpha);
		categDef1 * ("CT:="+thisAlpha+";");
		fprintf (stdout, "GT=");
		fscanf	(stdin,"Number",thisAlpha);
		categDef1 * ("GT:="+thisAlpha+";");
		
		fprintf (stdout, "AC2=");
		fscanf	(stdin,"Number",thisAlpha);
		categDef1 * ("AC2:="+thisAlpha+";");
		fprintf (stdout, "AT2=");
		fscanf	(stdin,"Number",thisAlpha);
		categDef1 * ("AT2:="+thisAlpha+";");
		fprintf (stdout, "CG2=");
		fscanf	(stdin,"Number",thisAlpha);
		categDef1 * ("CG2:="+thisAlpha+";");
		fprintf (stdout, "CT2=");
		fscanf	(stdin,"Number",thisAlpha);
		categDef1 * ("CT2:="+thisAlpha+";");
		fprintf (stdout, "GT2=");
		fscanf	(stdin,"Number",thisAlpha);
		categDef1 * ("GT2:="+thisAlpha+";");

		for (mi = 0; mi < resp; mi = mi+1)
		{
			fprintf (stdout, "Alpha for rate class ", mi+1, ":");
			fscanf	(stdin,"Number",thisAlpha);
			fprintf (stdout, "Beta for rate class ", mi+1, ":");
			fscanf	(stdin,"Number",thisBeta);
			fprintf (stdout, "Weight for rate class ", mi+1, ":");
			fscanf	(stdin,"Number",thisP);
			enteredValues [mi][0] = thisAlpha;
			enteredValues [mi][1] = thisBeta;
			enteredValues [mi][2] = thisP;
		}
		currentMult = 1-enteredValues[0][2];
		for (mi = 1; mi < resp-1; mi = mi+1)
		{
			enteredValues[mi][2] = enteredValues[mi][2]/currentMult;
			currentMult = currentMult*(1-enteredValues[mi][2]);
		}
		
		for (mi = 0; mi < resp-1; mi = mi+1)
		{
			categDef1 * ("S_"+mi+":>0;S_" + mi + ":=" + enteredValues[mi][0] + ";NS_" +mi+ ":=" + enteredValues[mi][1] + ";P_" + (mi+1) + ":=" +enteredValues[mi][2]+";");
		}
		categDef1 * ("S_"+mi+":>0;S_"+mi+":="+enteredValues[mi][0]+";NS_"+mi+":="+enteredValues[mi][1]+";");
	}

	if (stratDNDS)
	{
		for (mi=respM; mi<respN+respM; mi=mi+1)
		{
			categDef1 * ("\nR_"+mi+":=1;");
		}

		if (randomizeInitValues)
		{
			for (mi=0; mi<respM; mi=mi+1)
			{
				categDef1 * ("\nR_"+mi+":<1;R_"+mi+"="+Random(0.05,0.95)+";");
			}


			for (mi=respM+respN; mi<resp; mi=mi+1)
			{
				categDef1 * ("\nR_"+mi+":>1;R_"+mi+"="+Random(1.05,10)+";");
			}
		}
		else
		{
			for (mi=0; mi<respM; mi=mi+1)
			{
				categDef1 * ("\nR_"+mi+":<1;R_"+mi+"=1/"+(1+mi)+";");
			}

			for (mi=respM+respN; mi<resp; mi=mi+1)
			{
				categDef1 * ("\nR_"+mi+":>1;R_"+mi+"=2+"+(mi-respM-respN)+";");
			}
		}
	}
	else
	{
		if (randomizeInitValues == 1)
		{
			for (mi=0; mi<resp; mi=mi+1)
			{
				categDef1 * ("\nR_"+mi+"="+Random(0.05,1.75)+";");
			}
		}
		else
		{
			if (randomizeInitValues == 0)
			{
				for (mi=0; mi<resp; mi=mi+1)
				{
					categDef1 * ("\nR_"+mi+"="+(0.1+0.3*mi)+";");
				}	
			}
		}
	}

	categDef1 * 0;
	lfDef1 = "";
	lfDef1 * 128;
	lfDef1 * "\"";
	
	for (fileID = 0; fileID < fileCount; fileID = fileID + 1)
	{
		lfDef[fileID] * ")";
		lfDef[fileID] * 0;
		lfDef1 * lfDef[fileID];
		if (fileID < fileCount - 1)
		{
			lfDef1 * "+";
		}
	}

	lfDef1 	  * "\"";
	lfDef1	  * 0;

	ExecuteCommands (categDef1);

	ModelMatrixDimension = 64;
	for (h = 0 ;h<64; h=h+1)
	{
		if (_Genetic_Code[h]==10)
		{
			ModelMatrixDimension = ModelMatrixDimension-1;
		}
	}
	
	
	if (randomizeInitValues < 2)
	{
		SetDialogPrompt ("Save resulting fit to:");
		fprintf (PROMPT_FOR_FILE, CLEAR_FILE);
		resToPath = LAST_FILE_PATH;	
	}
	

	vectorOfFrequencies = BuildCodonFrequencies (observedFreq);

	ExecuteCommands (nucModelString+"\nModel nucModel = (nucModelMatrix,observedFreqSingle);");
	
	if (randomizeInitValues == 2)
	{
		fileID = ReportDistributionString(resp);
		fprintf (stdout, "\nUsing the following rate distribution\n", fileID, "\n");
	}

	nlfDef = "";
	nlfDef * 128;
	nlfDef * "LikelihoodFunction nuc_lf = (";

	for (fileID = 1; fileID <= fileCount; fileID = fileID + 1)
	{
		ExecuteCommands ("DataSetFilter nucFilter_"+fileID+" = CreateFilter (filteredData_"+fileID+",1);");
		ExecuteCommands ("Tree  nucTree_"+fileID+" = treeString_"+fileID+";");
		if (fileID > 1)
		{
			nlfDef * ",";
		}
		nlfDef * ("nucFilter_"+fileID+",nucTree_"+fileID);
	}
	nlfDef * 0;
	ExecuteCommands (nlfDef + ");");
	Optimize (nuc_res, nuc_lf);

	computeExpSubWeights (0);
	global codonFactor = 0.33;

	lfParts	= "";
	lfParts * 128;
	lfParts * "LikelihoodFunction lf = (filteredData_1,tree_1_0";

	for (fileID = 1; fileID <= fileCount; fileID = fileID + 1)
	{
		for (part = 0; part < resp; part = part + 1)
		{
			if (fileID == 1)
			{
				ExecuteCommands ("PopulateModelMatrix(\"rate_matrix_"+part+"\",observedFreq,\"S_"+part+"/c_scale\",\"NS_"+part+"/c_scale\");");
				ExecuteCommands ("Model MG94model_"+part+"= (rate_matrix_"+part+",vectorOfFrequencies,0);");
			}
			else
			{
				ExecuteCommands ("UseModel (MG94model_"+part+");");
			}
			
			treeID = "tree_"+fileID+"_"+part;
			ExecuteCommands ("Tree "+treeID+"=treeString_" + fileID +";");
			if (branchLengths)
			{
				ExecuteCommands ("ReplicateConstraint (\"this1.?.synRate:=this2.?.t__/codonFactor\","+treeID+",nucTree_"+fileID+");");
			}
			else
			{
				if (part == 0)
				{	
					ExecuteCommands ("bnames = BranchName(nucTree_" + fileID + ",-1);");
					nlfDef * 128;
					for (lc = 0; lc < Columns (bnames); lc=lc+1)
					{
						nlfDef * (treeID + "." + bnames[lc] + ".synRate = nucTree_" + fileID + "." + bnames[lc] + ".t/codonFactor;");
					}
					nlfDef * 0;
					ExecuteCommands (nlfDef);
				}
				else
				{
					ExecuteCommands ("ReplicateConstraint (\"this1.?.synRate:=this2.?.synRate\","+treeID+",tree_1_0);");
				}
			}
			if (part || fileID > 1)
			{
				lfParts = lfParts + ",filteredData_" + fileID + "," + treeID;
			}
		}
	}
	
	lfParts * 0;
	ExecuteCommands (lfParts + "," + lfDef1 + ");");

	USE_LAST_RESULTS = 1;

}
else
{
	SetDialogPrompt ("Choose a model fit:");
	ExecuteAFile (PROMPT_FOR_FILE);

	GetInformation(vars,"^P_[0-9]+$");
	resp = Columns (vars)+1;
	USE_LAST_RESULTS = 1;
	SKIP_CONJUGATE_GRADIENT = 1;
	
	ChoiceList (runKind, "Optimization type:", 1, SKIP_NONE,
							"Unconstrained", "No constraints on rate classes",
							"Constant dS", "Synonymous rates are constant",
							"No positive selection", "All dS <= dN",
							"HKY85", "Nucleotide substitution rates are forced to conform to HKY85");
					
	GetString (lfInfo, lf, -1);				
	fileCount = Columns (lfInfo["Datafilters"])/resp;	
	
	GetInformation (branchLengths, "^codonFactor$");
	branchLengths = (Columns (branchLengths)>0);
	
	/* do we have codonFactor? */			
							
	for (fileID = 1; fileID <= fileCount; fileID = fileID + 1)
	{
		ExecuteCommands ("totalCodonCount=totalCodonCount+filteredData_"+fileID+".sites;totalCharCount=totalCharCount+filteredData_"+fileID+".sites*filteredData_"+fileID+".species;");
	}


	if (runKind < 0)
	{
		return 0;
	}
	
	
	if (runKind == 0)
	{
		resToPath = LAST_FILE_PATH;
	}
	else
	{
		if (runKind == 1)
		{
			resToPath = LAST_FILE_PATH + ".dN_only";
			
			c_scale := 1;
			for (k = 0; k <resp; k=k+1)
			{
				ExecuteCommands("global R_"+k+"= Max(NS_"+k+",0.00001)/Max(S_"+k+",0.00001);");
				ExecuteCommands("S_"+k+":=1;");
				ExecuteCommands("NS_"+k+"=R_"+k+"*S_"+k+";");
			}
		}
		else
		{
			if (runKind == 2)
			{
				resToPath = LAST_FILE_PATH + ".NN";
				
				for (k = 0; k <resp; k=k+1)
				{
					ExecuteCommands("global R_"+k+";R_"+k+":<1;R_"+k+"=Max(NS_"+k+",0.00001)/Max(S_"+k+",0.00001);");
					ExecuteCommands("NS_"+k+":=R_"+k+"*S_"+k+";");
				}
			}		
			else
			{
				resToPath = LAST_FILE_PATH + ".HKY85";
				
				mi  = (1+CT)/2;
				mi2 = ;
				
				AC = (AC+AT+CG+GT)/(2*(1+CT));
				CT := 1;
				AT := AC;
				CG := AC;
				GT := AC;
			}		
		}
	}
	
	freqStrMx    = {resp,1};
	freqStrMx[0] = "P_1";

	for (mi=1; mi<resp-1; mi=mi+1)
	{
		freqStrMx[mi] = "";
		for (mi2=1;mi2<=mi;mi2=mi2+1)
		{
			freqStrMx[mi] = freqStrMx[mi]+"(1-P_"+mi2+")";		
		}
		freqStrMx[mi] = freqStrMx[mi]+"P_"+(mi+1);	
	}	
	freqStrMx[mi] = "";

	for (mi2=1;mi2<mi;mi2=mi2+1)
	{
		freqStrMx[mi] = freqStrMx[mi]+"(1-P_"+mi2+")";		
	}
	freqStrMx[mi] = freqStrMx[mi]+"(1-P_"+mi+")";
}

Optimize (res,lf);
SKIP_CONJUGATE_GRADIENT= 0;
USE_LAST_RESULTS       = 0;

LIKELIHOOD_FUNCTION_OUTPUT = 2;
sumPath = resToPath + ".summary";

LogL = res[1][0];

/* check for branch length approximator */

GetString (paramList, lf, -1);

degF = Columns(paramList["Global Independent"]) /* this will overcount by one; because the mean alpha := 1 */
	+ 8 /* 9 frequency parameters -1 */
	- branchLengths; /* remove one more for codon scaling factor, if nuc branch lengths are used */

for (fileID = 1; fileID <= fileCount; fileID = fileID + 1)
{
	ExecuteCommands ("degF = degF + BranchCount(tree_" + fileID + "_0) + TipCount(tree_" + fileID + "_0);");
}	

AIC	 = 2*(degF-LogL);
AICc = 2*(degF*totalCharCount/(totalCharCount-degF-1) - LogL);

fprintf (stdout, "\n\nModel fit summary\n", 
				 "\nLog likelihood:", Format (LogL, 15, 5),
				 "\nParameters    :", Format (degF, 15, 0),
				 "\nAIC           :", Format (AIC,  15, 5),
				 "\nc-AIC         :", Format (AICc,  15, 5),"\n"
);

fprintf (sumPath, CLEAR_FILE,"Model fit summary\n", 
				 "\nLog likelihood:", Format (LogL, 15, 5),
				 "\nParameters    :", Format (degF, 15, 0),
				 "\nAIC           :", Format (AIC,  15, 5),
				 "\nc-AIC         :", Format (AICc,  15, 5),"\n");
				 
fileID = ReportDistributionString(resp);

bivariateReturnAVL = {};
bivariateReturnAVL ["LogL"] 		= LogL;
bivariateReturnAVL ["Parameters"] 	= degF;
bivariateReturnAVL ["AIC"] 			= AIC;
bivariateReturnAVL ["cAIC"] 		= AICc;
bivariateReturnAVL ["Rates"] 		= resp;
bivariateReturnAVL ["AC"]			= AC;
bivariateReturnAVL ["AT"]			= AT;
bivariateReturnAVL ["CG"]			= CG;
bivariateReturnAVL ["CT"]			= CT;
bivariateReturnAVL ["GT"]			= GT;
bivariateReturnAVL ["Estimates"]	= reportMx;
				 

if (randomizeInitValues != 2)
{
	fprintf (stdout,  fileID);
	fprintf (sumPath, fileID);

	LIKELIHOOD_FUNCTION_OUTPUT = 7;
	if (runType == 1)
	{
		fprintf (resToPath,CLEAR_FILE,lf);
	}
	else
	{
		fprintf (resToPath,lf);
	}

	/*sumPath = resToPath + ".posterior";
	ConstructCategoryMatrix (catMat, lf, COMPLETE);

	posteriorProbs = {totalCodonCount, resp};

	catMatS = "";
	catMatS * 1024;

	for (mi2 = 0; mi2 < resp; mi2=mi2+1)
	{
		if (mi2)
		{
			catMatS * ",";
		}
		catMatS * ("Rate class " + (mi2+1) + " omega = " + reportMx[mi2][2]);
	}

	siteOffset = 0;
	for (fileID = 1; fileID <= fileCount; fileID = fileID + 1)
	{
		ExecuteCommands ("thisFilterSize=filteredData_" + fileID + ".sites;");
		for (mi=0; mi<thisFilterSize; mi=mi+1)
		{
			rs = 0;
			for (mi2 = 0; mi2 < resp; mi2=mi2+1)
			{
				posteriorProbs[mi+siteOffset][mi2] = catMat[mi+siteOffset+mi2*totalCodonCount] * reportMx[mi2][3];
				rs = rs + posteriorProbs[mi+siteOffset][mi2];
			}
			catMatS * "\n";
			for (mi2 = 0; mi2 < resp; mi2=mi2+1)
			{
				posteriorProbs[mi+siteOffset][mi2] = posteriorProbs[mi+siteOffset][mi2]/rs;
				if (mi2)
				{
					catMatS * ",";
				}
				catMatS * (""+posteriorProbs[mi+siteOffset][mi2]);
			}
		} 
		siteOffset = siteOffset + thisFilterSize;
	}
	catMatS * 0;
	fprintf (sumPath, CLEAR_FILE, catMatS);*/
}


return bivariateReturnAVL;
