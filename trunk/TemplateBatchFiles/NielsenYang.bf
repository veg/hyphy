ModelNames = {{"Neutral",
	  		  "Selection",
			  "Discrete",
			  "Freqs",
			  "Gamma",
			  "2 Gamma",
			  "Beta",
			  "Beta & w",
			  "Beta & Gamma",
			  "Beta & (Gamma+1)",
			  "Beta & (Normal>1)",
			  "0 & 2 (Normal>1)",
			  "3 Normal",
			  "Gamma mod Beta",
			  "Beta & 1"}};
			  			       
MAXIMUM_ITERATIONS_PER_VARIABLE = 2000;
OPTIMIZATION_PRECISION 			= 0.001;

modelLL							= {};
			  
/*-------------------------------------------------------------------------------*/

function PromptForAParameter (promptString, leftBound, rightBound, midx)
{
	if (noMorePrompting)
	{
		return -1;
		
	}
	
	value = leftBound - 1;
	while ((value < leftBound)||(value > rightBound))
	{
		fprintf (stdout, "\n[", ModelNames[midx], "] ", promptString, " in [", leftBound, ",", rightBound, "] (-2 to quit):");
		fscanf	(stdin , "Number", value);
		if (value == (-2))
		{
			noMorePrompting = -1;
			return -1;
		}
	}
	
	return value;
}

/*-------------------------------------------------------------------------------*/

function PromptForWeights (promptString, wc,  midx)
{
	resMx = {wc,1};
	if (noMorePrompting)
	{
		return resMx;
	}
	
	wci	  = 0;
	currentMax = 1;
		
	for (widx = 0; (widx < wc-1)&&(currentMax>0); widx = widx+1)
	{
		resMx[wci] = PromptForAParameter (promptString[wci], 0, currentMax, midx);
		if (resMx[wci] == (-1))
		{
			return resMx;
		}
		currentMax = currentMax - resMx[wci];
		wci = wci+1;
	}
	
	
	leftOver = (1-resMx[0]);
	
	for (widx = 1; widx < wc-1; widx = widx+1)
	{
		if (leftOver == 0)
		{
			return resMx;
		}
		resMx [widx] = resMx [widx]/leftOver;
		leftOver = leftOver*(1-resMx[widx]);
	}
	
	return resMx;
}

/*-------------------------------------------------------------------------------*/

function WriteSnapshot (modelIdx)
{
	likOF = LIKELIHOOD_FUNCTION_OUTPUT;
	LIKELIHOOD_FUNCTION_OUTPUT = 6;
	SPOOL_FILE = SUMMARY_FILE+"_MODEL_"+modelIdx+".nex";
	fprintf (SPOOL_FILE, CLEAR_FILE, lf);
	LIKELIHOOD_FUNCTION_OUTPUT = likOF;
	return 0;
}

/*-------------------------------------------------------------------------------*/

function PromptForModelParameter (modelIdx)
{
	noMorePrompting = 0;
	if (modelIdx == 0)
	{
		global P = PromptForAParameter ("Weight for the '0' rate class", 0, 1, modelIdx);
	}
	else
	{
		if (modelIdx == 1)
		{
			promptMx  = {{"Weight for the '0' rate class", "Weight for the '1' rate class"}};
			wtMx	  = PromptForWeights (promptMx, 3, modelIdx);
			global P1 = wtMx[0];
			global P2 = wtMx[1];
			global W = PromptForAParameter ("Adjustable rate value", 0, 10000 , modelIdx);
		}		
		else
		{
			if (modelIdx == 2)
			{
				promptMx  = {{"Weight for the lowest rate class", "Weight for the middle rate class"}};
				wtMx	  = PromptForWeights (promptMx, 3, modelIdx);
				global P1 = wtMx[0];
				global P2 = wtMx[1];

				global R_1 = PromptForAParameter ("Lowest Rate Value",  0,  10000 , modelIdx);
				global W1  = PromptForAParameter ("Middle Rate Value", R_1, 10000 , modelIdx);
				global R_2 = PromptForAParameter ("Highest Rate Value", W1, 10000 , modelIdx);

				R_1	= R_1/W1;
				R_2 = R_2/W1;
			}
			else
			{
				if (modelIdx == 3)
				{
					promptMx  = {{"Weight for the '0' rate class", "Weight for the '1/3' rate class", "Weight for the '2/3' rate class",  
								  "Weight for the '1' rate class"}};

					wtMx	  = PromptForWeights (promptMx, 5, modelIdx);
					global P1 = wtMx[0];
					global P2 = wtMx[1];
					global P3 = wtMx[2];
					global P4 = wtMx[3];
				}
				else
				{
					if ((modelIdx == 4)||(modelIdx == 5))
					{
						global alpha = PromptForAParameter ("Alpha parameter for the gamma distribution", 0.01, 100 , modelIdx);
						global beta  = PromptForAParameter ("Beta  parameter for the gamma distribution", 0.01, 200 , modelIdx);
						if (modelIdx == 5)
						{
							global alpha2 = PromptForAParameter ("Alpha parameter for the 2nd gamma distribution", 0.01, 100 , modelIdx);
						}
					}
					else
					{
						if ((modelIdx == 6)||(modelIdx == 7)||(modelIdx == 8)||(modelIdx == 9)||(modelIdx == 10)||(modelIdx == 14))
						{
							global betaP = PromptForAParameter ("P parameter for the beta distribution", 0.05, 85 , modelIdx);
							global betaQ = PromptForAParameter ("Q parameter for the beta distribution", 0.05, 85 , modelIdx);
							global P 	 = PromptForAParameter ("Mixture weight for the beta distribution", 0, 1, modelIdx);
							if (modelIdx == 7)
							{
								global W = PromptForAParameter ("Value for omega", 1, 10000 , modelIdx);
							}
							if ((modelIdx == 8)||(modelIdx == 9))
							{
								global alpha = PromptForAParameter ("Alpha parameter for the gamma distribution", 0.01, 100 , modelIdx);
							}
							if (modelIdx == 9)
							{
								global beta = PromptForAParameter ("Beta parameter for the gamma distribution", 0.01, 500 , modelIdx);
							}
							if (modelIdx == 10)
							{
								global mu 	 = PromptForAParameter ("mu parameter for the normal distribution", 0, 10000 , modelIdx);
								global sigma = PromptForAParameter ("sigma parameter for the normal distribution", 0.0001, 10000 , modelIdx);
							}
							if (modelIdx == 14)
							{
								global W := 1;
							}
						}
						else
						{
							if (modelIdx == 11)
							{
								promptMx  = {{"Weight for the '0' rate class", "Weight for the first normal in the mixture"}};
								wtMx	  = PromptForWeights (promptMx, 3, modelIdx);
								P 		  = wtMx[0];
								P1		  = wtMx[1];
								
								global mu 	  = PromptForAParameter ("mu parameter for the first     normal distribution", 0, 10000 , modelIdx);
								global sigma  = PromptForAParameter ("sigma parameter for the first  normal distribution", 0.0001, 10000 , modelIdx);
								global sigma1 = PromptForAParameter ("sigma parameter for the second normal distribution", 0.0001, 10000 , modelIdx);
							}
							else
							{
								if (modelIdx == 12)
								{
									promptMx  = {{"Weight for the first normal in the mixture", "Weight for the second normal in the mixture"}};
									wtMx	  = PromptForWeights (promptMx, 3, modelIdx);
									P 		  = wtMx[0];
									P1		  = wtMx[1];
							
									global mu 	  = PromptForAParameter ("mu parameter for the first     normal distribution", 0, 10000 , modelIdx);
									global sigma  = PromptForAParameter ("sigma parameter for the first  normal distribution", 0.0001, 10000 , modelIdx);
									global sigma1 = PromptForAParameter ("sigma parameter for the second normal distribution", 0.0001, 10000 , modelIdx);
									global sigma2 = PromptForAParameter ("sigma parameter for the third normal distribution", 0.0001, 10000 , modelIdx);
								}
								if (modelIdx == 13)
								{
									global alpha = PromptForAParameter ("Alpha parameter for the gamma distribution", 0.01, 100 , modelIdx);
									global beta  = PromptForAParameter ("Beta  parameter for the gamma distribution", 0.01, 200 , modelIdx);
									global betaP = PromptForAParameter ("P parameter for the beta distribution", 0.05, 85 , modelIdx);
									global betaQ = PromptForAParameter ("Q parameter for the beta distribution", 0.05, 85 , modelIdx);
								}
							}
						}
					}				
				}	
			}	
		}	
	}		
	return noMorePrompting;	
}

#include "TemplateModels/Yang2000Distributions.def";

/* ____________________________________________________________________________________________________________________*/

function FrameText (frameChar,vertChar,parOff,theText)
{
	h = Abs (theText)+4;
	fprintf (stdout,"\n");	
	for (k=0; k<parOff; k=k+1)
	{
		fprintf (stdout," ");
	}
	for (k=0; k<h;k=k+1)
	{
		fprintf (stdout,frameChar);
	}
	fprintf (stdout,"\n");	
	for (k=0; k<parOff; k=k+1)
	{
		fprintf (stdout," ");
	}
	fprintf (stdout,vertChar," ",theText," ",vertChar,"\n");
	for (k=0; k<parOff; k=k+1)
	{
		fprintf (stdout," ");
	}
	for (k=0; k<h;k=k+1)
	{
		fprintf (stdout,frameChar);
	}
	fprintf (stdout,"\n");	
	return 0;
}

/* ____________________________________________________________________________________________________________________*/

function BuildCodonFrequencies4 (obsF)
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
			PIStop = PIStop-obsF[first][0]*obsF[second][0]*obsF[third][0];
			continue; 
		}
		result[h-hshift][0]=obsF[first][0]*obsF[second][0]*obsF[third][0];
	}
	return result*(1.0/PIStop);
}

/* ____________________________________________________________________________________________________________________*/

function BuildCodonFrequencies12 (obsF)
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
/* ____________________________________________________________________________________________________________________*/


function GetDistributionParameters (sigLevel)
{
	GetInformation (distrInfo,c);
	D = Columns(distrInfo);
	E = 0.0;
	T = 0.0;
	sampleVar = 0.0;
	for (k=0; k<D; k=k+1)
	{
		T = distrInfo[0][k]*distrInfo[1][k];
		E = E+T;
		sampleVar = T*distrInfo[0][k]+sampleVar;
	}
	sampleVar = sampleVar-E*E;

	fprintf  (SUMMARY_FILE,"\n\n------------------------------------------------\n\ndN/dS = ",E, " (sample variance = ",sampleVar,")\n");
	for (k=0; k<D; k=k+1)
	{
		fprintf (SUMMARY_FILE,"\nRate[",Format(k+1,0,0),"]=",
				 Format(distrInfo[0][k],12,8), " (weight=", Format(distrInfo[1][k],9,7),")");
	}
	
	for (k=0; k<D; k=k+1)
	{
		if (distrInfo[0][k]>1) break;
	}
	ConstructCategoryMatrix(marginals_1,lf,COMPLETE);
	CC = Columns (marginals_1);
	if (k<D)
	/* have rates > 1 */
	{
		marginals = marginals_1;
		fprintf  (SUMMARY_FILE,"\n\n------------------------------------------------\n\n Sites with dN/dS>1 (Posterior cutoff = ",sigLevel,")\n\n");
		for (v=0; v<CC; v=v+1)
		{
			sampleVar = 0;
			for (h=0; h<D; h=h+1)
			{
				sampleVar = sampleVar+distrInfo[1][h]*marginals[h][v];
			}
			positiveProb = 0;
			for (l=k; l<D; l=l+1)
			{
				positiveProb = positiveProb+distrInfo[1][l]*marginals[l][v];
			}
			positiveProb = positiveProb/sampleVar;
			marginals[0][v] = positiveProb;
			if (positiveProb>=sigLevel)
			{
				fprintf (SUMMARY_FILE,Format (v+1,0,0)," (",positiveProb,")\n");
			}
		}
		fprintf  (SUMMARY_FILE,"\n\n------------------------------------------------\n\n Sites with dN/dS<=1 (Posterior cutoff = ",sigLevel,")\n\n");
		for (v=0; v<CC; v=v+1)
		{
			if (marginals[0][v]<sigLevel)
			{
				fprintf (SUMMARY_FILE,Format (v+1,0,0)," (",marginals[0][v],")\n");
			}
		}
		marginals = 0;
	}
	else
	{
		fprintf  (SUMMARY_FILE,"\n\n------------------------------------------------\n\n No rate classes with dN/dS>1.");
	}
	fprintf  (SUMMARY_FILE,"\n\n------------------------------------------------\n\nTabulated Posteriors for Each Site\nSites in Columns, Rate Classes in Rows\nImport the following part into a data processing program\nfor further analysis\n\n");
	for (v=0; v<CC; v=v+1)
	{
		sampleVar = 0;
		for (h=0; h<D; h=h+1)
		{
			sampleVar = sampleVar+distrInfo[1][h]*marginals_1[h][v];
		}
		for (l=0; l<D; l=l+1)
		{
			marginals_1[l][v] = distrInfo[1][l]*marginals_1[l][v]/sampleVar;
		}
	}
	outString = "Rate/Site";
	for (v=0; v<D; v=v+1)
	{
		outString = outString + "\tRate=" + distrInfo[0][v];
	}	
	
	for (v=0; v<CC; v=v+1)
	{
		outString = outString + "\n" + (v+1);
		for (h=0; h<D; h=h+1)
		{
			outString = outString + "\t" + marginals_1[h][v];
		}
	}

	fprintf  (SUMMARY_FILE,"\n", outString,"\n");
	
	marginals_1 = 0;
	
	return E;
}

/* ____________________________________________________________________________________________________________________*/

function PopulateModelMatrix (ModelMatrixName&, EFV)
{
	ModelMatrixName = {ModelMatrixDimension,ModelMatrixDimension}; 

	hshift = 0;
	
	if (modelType==0)
	{
		for (h=0; h<64; h=h+1)
		{
			if (_Genetic_Code[h]==10) 
			{
				hshift = hshift+1;
				continue; 
			}
			vshift = hshift;
			for (v = h+1; v<64; v=v+1)
			{
				diff = v-h;
				if (_Genetic_Code[v]==10) 
				{
					vshift = vshift+1;
					continue; 
				}
			  	if ((h$4==v$4)||((diff%4==0)&&(h$16==v$16))||(diff%16==0))
			  	{
			  		if (h$4==v$4)
			  		{
			  			transition = v%4;
			  			transition2= h%4;
			  		}
			  		else
			  		{
			  			if(diff%16==0)
			  			{
			  				transition = v$16;
			  				transition2= h$16;
			  			}
			  			else
			  			{
			  				transition = v%16$4;
			  				transition2= h%16$4;
			  			}
			  		}
			  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
			  		{
			  			ModelMatrixName[h-hshift][v-vshift] := t*EFV__[transition__];
			  			ModelMatrixName[v-vshift][h-hshift] := t*EFV__[transition2__];
				  	}
			  		else
			  		{
				  		ModelMatrixName[h-hshift][v-vshift] := c*t*EFV__[transition__];
			  			ModelMatrixName[v-vshift][h-hshift] := c*t*EFV__[transition2__];
		  			}
			  	}
			  }
		}
	}
	else
	{
		if (modelType==1)
		{
			for (h=0; h<64; h=h+1)
			{
				if (_Genetic_Code[h]==10) 
				{
					hshift = hshift+1;
					continue; 
				}
				vshift = hshift;
				for (v = h+1; v<64; v=v+1)
				{
					diff = v-h;
					if (_Genetic_Code[v]==10) 
					{
						vshift = vshift+1;
						continue; 
					}
					nucPosInCodon = 2;
				  	if ((h$4==v$4)||((diff%4==0)&&(h$16==v$16))||(diff%16==0))
				  	{
				  		if (h$4==v$4)
				  		{
				  			transition = v%4;
				  			transition2= h%4;
				  		}
				  		else
				  		{
				  			if(diff%16==0)
				  			{
				  				transition = v$16;
				  				transition2= h$16;
								nucPosInCodon = 0;
				  			}
				  			else
				  			{
				  				transition = v%16$4;
				  				transition2= h%16$4;
								nucPosInCodon = 1;
				  			}
				  		}
				  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
				  		{
				  			ModelMatrixName[h-hshift][v-vshift] := t*EFV__[transition__][nucPosInCodon__];
				  			ModelMatrixName[v-vshift][h-hshift] := t*EFV__[transition2__][nucPosInCodon__];
					  	}
				  		else
				  		{
					  		ModelMatrixName[h-hshift][v-vshift] := c*t*EFV__[transition__][nucPosInCodon__];
				  			ModelMatrixName[v-vshift][h-hshift] := c*t*EFV__[transition2__][nucPosInCodon__];
			  			}
				  	}
				  }
			}
		}
		else
		{
			if (modelType<4)
			{
				for (h=0; h<64; h=h+1)
				{
					if (_Genetic_Code[h]==10) 
					{
						hshift = hshift+1;
						continue; 
					}
					vshift = hshift;
					for (v = h+1; v<64; v=v+1)
					{
						diff = v-h;
						if (_Genetic_Code[v]==10) 
						{
							vshift = vshift+1;
							continue; 
						}
					  	if ((h$4==v$4)||((diff%4==0)&&(h$16==v$16))||(diff%16==0))
					  	{
					  		if (h$4==v$4)
					  		{
					  			transition = v%4;
					  			transition2= h%4;
					  		}
					  		else
					  		{
					  			if(diff%16==0)
					  			{
					  				transition = v$16;
					  				transition2= h$16;
					  			}
					  			else
					  			{
					  				transition = v%16$4;
					  				transition2= h%16$4;
					  			}
					  		}
					  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
					  		{
					  			if (Abs(transition-transition2)%2)
					  			{
					  				ModelMatrixName[h-hshift][v-vshift] := kappa*t;
					  				ModelMatrixName[v-vshift][h-hshift] := kappa*t;
					  			}
					  			else
					  			{
					  				ModelMatrixName[h-hshift][v-vshift] := t;
					  				ModelMatrixName[v-vshift][h-hshift] := t;
					  			}
					  			
						  	}
					  		else
					  		{
					  			if (Abs(transition-transition2)%2)
					  			{
					  				ModelMatrixName[h-hshift][v-vshift] := kappa*c*t;
					  				ModelMatrixName[v-vshift][h-hshift] := kappa*c*t;
					  			}
					  			else
					  			{
					  				ModelMatrixName[h-hshift][v-vshift] := c*t;
					  				ModelMatrixName[v-vshift][h-hshift] := c*t;
					  			}
						  	}
					  	}	
					 }
				}	
			}
			else
			{
				for (h=0; h<64; h=h+1)
				{
					if (_Genetic_Code[h]==10) 
					{
						hshift = hshift+1;
						continue; 
					}
					vshift = hshift;
					for (v = h+1; v<64; v=v+1)
					{
						diff = v-h;
						if (_Genetic_Code[v]==10) 
						{
							vshift = vshift+1;
							continue; 
						}
						nucPosInCodon = 2;
					  	if ((h$4==v$4)||((diff%4==0)&&(h$16==v$16))||(diff%16==0))
					  	{
					  		if (h$4==v$4)
					  		{
					  			transition = v%4;
					  			transition2= h%4;
					  		}
					  		else
					  		{
					  			if(diff%16==0)
					  			{
					  				transition = v$16;
					  				transition2= h$16;
									nucPosInCodon = 0;
					  			}
					  			else
					  			{
					  				transition = v%16$4;
					  				transition2= h%16$4;
									nucPosInCodon = 1;
					  			}
					  		}
					  		
					  		if (modelType == 4)
					  		{
					  			nucPosInCodon = 0;
					  		}
					  		
					  		rateKind = mSpecMatrix[transition][transition2];
					  		
					  		if (rateKind == 1)
					  		{
						  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
						  		{
						  			ModelMatrixName[h-hshift][v-vshift] := AC*synRate*EFV__[transition__][nucPosInCodon__];
						  			ModelMatrixName[v-vshift][h-hshift] := AC*synRate*EFV__[transition2__][nucPosInCodon__];
							  	}
						  		else
						  		{
							  		ModelMatrixName[h-hshift][v-vshift] := AC*c*synRate*EFV__[transition__][nucPosInCodon__];
						  			ModelMatrixName[v-vshift][h-hshift] := AC*c*synRate*EFV__[transition2__][nucPosInCodon__];
					  			}
					  		}
					  		else
					  		{
						  		if (rateKind == 2)
						  		{
							  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
							  		{
							  			ModelMatrixName[h-hshift][v-vshift] := synRate*EFV__[transition__][nucPosInCodon__];
							  			ModelMatrixName[v-vshift][h-hshift] := synRate*EFV__[transition2__][nucPosInCodon__];
								  	}
							  		else
							  		{
								  		ModelMatrixName[h-hshift][v-vshift] := c*synRate*EFV__[transition__][nucPosInCodon__];
							  			ModelMatrixName[v-vshift][h-hshift] := c*synRate*EFV__[transition2__][nucPosInCodon__];
						  			}
						  		}
						  		else
						  		{
							  		if (rateKind == 3)
							  		{
								  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
								  		{
								  			ModelMatrixName[h-hshift][v-vshift] := AT*synRate*EFV__[transition__][nucPosInCodon__];
								  			ModelMatrixName[v-vshift][h-hshift] := AT*synRate*EFV__[transition2__][nucPosInCodon__];
									  	}
								  		else
								  		{
									  		ModelMatrixName[h-hshift][v-vshift] := AT*c*synRate*EFV__[transition__][nucPosInCodon__];
								  			ModelMatrixName[v-vshift][h-hshift] := AT*c*synRate*EFV__[transition2__][nucPosInCodon__];
							  			}
							  		}
							  		else
							  		{
								  		if (rateKind == 4)
								  		{
									  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
									  		{
									  			ModelMatrixName[h-hshift][v-vshift] := CG*synRate*EFV__[transition__][nucPosInCodon__];
									  			ModelMatrixName[v-vshift][h-hshift] := CG*synRate*EFV__[transition2__][nucPosInCodon__];
										  	}
									  		else
									  		{
										  		ModelMatrixName[h-hshift][v-vshift] := CG*c*synRate*EFV__[transition__][nucPosInCodon__];
									  			ModelMatrixName[v-vshift][h-hshift] := CG*c*synRate*EFV__[transition2__][nucPosInCodon__];
								  			}
								  		}
								  		else
								  		{
									  		if (rateKind == 5)
									  		{
										  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
										  		{
										  			ModelMatrixName[h-hshift][v-vshift] := CT*synRate*EFV__[transition__][nucPosInCodon__];
										  			ModelMatrixName[v-vshift][h-hshift] := CT*synRate*EFV__[transition2__][nucPosInCodon__];
											  	}
										  		else
										  		{
											  		ModelMatrixName[h-hshift][v-vshift] := CT*c*synRate*EFV__[transition__][nucPosInCodon__];
										  			ModelMatrixName[v-vshift][h-hshift] := CT*c*synRate*EFV__[transition2__][nucPosInCodon__];
									  			}
									  		}
									  		else
									  		{
										  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
										  		{
										  			ModelMatrixName[h-hshift][v-vshift] := GT*synRate*EFV__[transition__][nucPosInCodon__];
										  			ModelMatrixName[v-vshift][h-hshift] := GT*synRate*EFV__[transition2__][nucPosInCodon__];
											  	}
										  		else
										  		{
											  		ModelMatrixName[h-hshift][v-vshift] := GT*c*synRate*EFV__[transition__][nucPosInCodon__];
										  			ModelMatrixName[v-vshift][h-hshift] := GT*c*synRate*EFV__[transition2__][nucPosInCodon__];
									  			}
									  		}
									  	}
								  	}
							  	}
						  	}
					  	}
					}
				}	
			}
		}
	 }
	 return ((modelType>1)&&(modelType<4));
}

/* ____________________________________________________________________________________________________________________*/

NICETY_LEVEL = 3;

#include "TemplateModels/chooseGeneticCode.def";

ModelMatrixDimension = 64;
for (h = 0 ;h<64; h=h+1)
{
	if (_Genetic_Code[h]==10)
	{
		ModelMatrixDimension = ModelMatrixDimension-1;
	}
}

SetDialogPrompt ("Please specify a codon data file:");

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);

DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);

dummyVar = FrameText ("-","|",2,"READ THE FOLLOWING DATA");

fprintf (stdout,"\n",ds);

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
		treeString = DATAFILE_TREE;
		IS_TREE_PRESENT_IN_DATA = 1;
	}
	fprintf (stdout, "\n\n");

}

if (!IS_TREE_PRESENT_IN_DATA)
{
	SetDialogPrompt ("Please select a tree file for the data:");
	fscanf (PROMPT_FOR_FILE, "String", treeString);
}

treeString = RerootTree (treeString,0);

chosenModelList = {16,1};

ChoiceList (modelType,"Distributions",1,SKIP_NONE,
			"Run All","Run all available dN/dS distributions",
			"Run Custom","Choose from available dN/dS distributions.");
			
if (modelType<0)
{
	return;
}

if (modelType==0)
{
	for (rateType = 0; rateType<16; rateType=rateType+1)
	{
		chosenModelList[rateType][0] = 1;
	}
}
else
{
	ChoiceList (modelTypes,"Distributions",0,SKIP_NONE,
				"Single Rate","Single Rate",
				"Neutral","Neutral",
	  			"Selection","Selection",
			    "Discrete","Discrete",
			    "Freqs","Freqs",
			    "Gamma","Gamma",
			    "2 Gamma","2 Gamma",
			    "Beta","Beta",
			    "Beta & w","Beta & w",
			    "Beta & Gamma","Beta & Gamma",
			    "Beta & (Gamma+1)","Beta & (Gamma+1)",
			    "Beta & (Normal>1)","Beta & (Normal>1)",
			    "0 & 2 (Normal>1)","0 & 2 (Normal>1)",
			    "3 Normal","3 Normal",
			    "Gamma mod Beta","Two parameter gamma partitioned by a two parameter beta",
			    "Beta & 1", "Discrete beta and a point mass at one. The M8a model from PAML.");
			    
	if (modelTypes[0]<0)
	{
		return;
	}
	
	ChoiceList (initValue,"Initial Values",1,SKIP_NONE,
			    "Default", "Use default initial values for distribution parameters",
			    "Custom", "Enter custom initial values for distribution parameters");
			    
	if (initValue<0)
	{
		return;
	}

	for (rateType = 0; rateType < Rows(modelTypes)*Columns(modelTypes); rateType = rateType + 1)
	{
		modelType = modelTypes[rateType];
		chosenModelList[modelType] = 1;
	}
}


if (initValue)
{
	for (rateType = 0; rateType<15; rateType=rateType+1)
	{
		if (chosenModelList[rateType+1]==1)
		{
			if (PromptForModelParameter (rateType))
			{
				return 0;
			}
		}
	}
}


ChoiceList (modelType,"Choose a model",1,SKIP_NONE,
			"MG94 1x4","Muse-Gaut 94 model with 4(-1) nucleotide frequency parameters (intra-codon position independent).",
			"MG94 3x4","Muse-Gaut 94 model with 12(-3) nucleotide frequency parameters (intra-codon position specific).",
			"GY94 1x4","Goldman-Yang 94 model with 4(-1) nucleotide frequency parameters (intra-codon position independent).",
			"GY94 3x4","Goldman-Yang 94 model with 12(-3) nucleotide frequency parameters (intra-codon position specific).",
			"MG94 Custom 1x4", "Muse-Gaut 94 crossed with an arbitrary nucelotide substitution model with 4(-1) nucleotide frequency parameters (intra-codon position independent).",
			"MG94 Custom 3x4", "Muse-Gaut 94 crossed with an arbitrary nucelotide substitution model with 12(-3) nucleotide frequency parameters (intra-codon position specific)."
);


if ((modelType==0)||(modelType==2)||(modelType==4))
{
	freqParamCount = 3;
	HarvestFrequencies (observedFreq,filteredData,1,1,0);
	vectorOfFrequencies = BuildCodonFrequencies4 (observedFreq);
}
else
{
	freqParamCount = 9;
	HarvestFrequencies (observedFreq,filteredData,3,1,1);
	vectorOfFrequencies = BuildCodonFrequencies12 (observedFreq);
}

if (modelType>3)
{
	global AC=1;
	global AT=1;
	global CG=1;
	global CT=1;
	global GT=1;

	done = 0;

	while (!done)
	{
		fprintf (stdout,"\nPlease enter a 6 character model designation (e.g:010010 defines HKY85):");
		fscanf  (stdin,"String", modelDesc);
		if (Abs(modelDesc)==6)
		{	
			done = 1;
		}
	}
	
	done = 0;
	mSpecMatrix = {{*,1,2,3}{1,*,4,5}{2,4,*,6}{3,5,6,*}};

	rateBiasTerms = {{"AC","1","AT","CG","CT","GT"}};
	paramCount	  = 0;

	modelConstraintString = "";

	for (customLoopCounter2=1; customLoopCounter2<6; customLoopCounter2=customLoopCounter2+1)
	{
		for (customLoopCounter=0; customLoopCounter<customLoopCounter2; customLoopCounter=customLoopCounter+1)
		{
			if (modelDesc[customLoopCounter2]==modelDesc[customLoopCounter])
			{
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
	}	

	if (Abs(modelConstraintString))
	{
		ExecuteCommands (modelConstraintString);
	}
}
else
{
	if (modelType>1)
	{	
		global kappa = 2.;
	}
}

fprintf (stdout, "\nChoose the cutoff (0 to 1) for posterior of dN/dS>1 for a site to be considered under selective pressure:");
fscanf  (stdin, "Number",psigLevel);
if ((psigLevel <= 0)||(psigLevel>1))
{
	psigLevel = .95;
}
fprintf (stdout, "\n>Using ", psigLevel , " cutoff\n");

for (rateType = 5; rateType < 16; rateType = rateType + 1)
{
	if (chosenModelList[rateType]==1)
	{
		categCount = -1;
		while (categCount <= 0)
		{
			fprintf (stdout, "\nChoose the number of categories in discretized distributions:");
			fscanf  (stdin, "Number",categCount);
			categCount = categCount$1;
		}
		break;		
	}
}

fprintf (stdout, "\n>Using ", Format (categCount,0,0), " categories.\n");

SetDialogPrompt ("Write detailed results to:");
fprintf 		(PROMPT_FOR_FILE,CLEAR_FILE);
SUMMARY_FILE = LAST_FILE_PATH;

ClearConstraints (c);
global c = 1.;

dummyVar = FrameText ("-","|",2,"SUMMARY TABLE");
tableSeparator =  "+-------------------------+----------------+---------------+---------+\n";
fprintf (stdout, "\n\"p\" is the number of estimated model parameters.\nDetailed results including sites with dN/dS>1 will be written to\n",SUMMARY_FILE,"\n\n");
fprintf (stdout, tableSeparator,
				 "| MODEL (Number & Desc)   | Log likelihood |     dN/dS     |    p    |\n",
				                                             
				 tableSeparator);
				 
cachedBranchLengths = {{-1,-1}};
			
rateType = -1;

if (MPI_NODE_COUNT>1)
{
	MPINodeState = {MPI_NODE_COUNT-1,3};
	OPTIMIZE_SUMMATION_ORDER = 0;
}
				 
if (chosenModelList[0]>0)
{
	timer = Time(1);
	modelMatrix = 0;
	MULTIPLY_BY_FREQS = PopulateModelMatrix ("modelMatrix", observedFreq);
	Model theModel = (modelMatrix,vectorOfFrequencies,MULTIPLY_BY_FREQS);
	Tree  givenTree = treeString;
	DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
	LikelihoodFunction lf 	   = (filteredData,givenTree);
					 
	if (MPI_NODE_COUNT>1)
	{
		for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
		{
			if (MPINodeState[mpiNode][0]==0)
			{
				break;	
			}
		}
		
		if (mpiNode==MPI_NODE_COUNT-1)
		/* all nodes busy */
		{
			mpiNode = ReceiveJobs (1);
		}
		else
		{
			MPISend (mpiNode+1,lf);
			MPINodeState[mpiNode][0] = 1;
			MPINodeState[mpiNode][1] = rateType;
			MPINodeState[mpiNode][2] = Time(1);
		}
	}
	else
	{
		timer = Time(1);
		Optimize (res,lf);
		ReceiveJobs (0);
	}

	/*timer = res[1][1]-res[1][2];
	cachedBranchLengths = {timer,1};
	
	for (rateType = timer; rateType < Columns(cachedBranchLengths); rateType = rateType+1)
	{
		cachedBranchLengths[rateType-timer][0] = res [0][rateType];
	}*/
}

for (rateType = 0; rateType < 15; rateType = rateType + 1)
{
	if (chosenModelList[rateType+1]==0)
	{
		continue;
	}

	SetWDistribution (categCount);
	modelMatrix = 0;
	MULTIPLY_BY_FREQS = PopulateModelMatrix ("modelMatrix", observedFreq);
	Model theModel = (modelMatrix,vectorOfFrequencies,MULTIPLY_BY_FREQS);
	Tree  givenTree = treeString;
	DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);
	LikelihoodFunction lf = (filteredData,givenTree);

	if (MPI_NODE_COUNT>1)
	{
		for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
		{
			if (MPINodeState[mpiNode][0]==0)
			{
				break;	
			}
		}
		
		if (mpiNode==MPI_NODE_COUNT-1)
		/* all nodes busy */
		{
			mpiNode = ReceiveJobs (1);
		}
		else
		{
			MPISend (mpiNode+1,lf);
			MPINodeState[mpiNode][0] = 1;
			MPINodeState[mpiNode][1] = rateType;
			MPINodeState[mpiNode][2] = Time(1);
		}
	}
	else
	{
		timer = Time(1);
		Optimize (res,lf);
		ReceiveJobs (0);
	}

	if (modelType>1)
	{	
		kappa = 2.;
	}
}

if (MPI_NODE_COUNT>1)
{
	while (1)
	{
		for (nodeCounter = 0; nodeCounter < MPI_NODE_COUNT-1; nodeCounter = nodeCounter+1)
		{
			if (MPINodeState[nodeCounter][0]==1)
			{
				fromNode = ReceiveJobs (0);
				break;	
			}
		}
		if (nodeCounter == MPI_NODE_COUNT-1)
		{
			break;
		}
	}	
	OPTIMIZE_SUMMATION_ORDER = 1;
}

/* ____________________________________________________________________________________________________________________*/

function ReceiveJobs (sendOrNot)
{
	if (MPI_NODE_COUNT>1)
	{
		MPIReceive (-1, fromNode, result_String);
		rateType2 = MPINodeState[fromNode-1][1];
		timer     = MPINodeState[fromNode-1][2];
		
		if (sendOrNot)
		{
			MPISend (fromNode,lf);
			MPINodeState[fromNode-1][1] = rateType;
			MPINodeState[fromNode-1][2] = Time(1);
		}
		else
		{
			MPINodeState[fromNode-1][0] = 0;
			MPINodeState[fromNode-1][1] = 0;
			MPINodeState[fromNode-1][2] = 0;
		}
		
		rateType3 = rateType;
		rateType  = rateType2;
		
		if (rateType < 0)
		{
			ClearConstraints (c);
		}
				
		ExecuteCommands (result_String);
				
		res = lf_MLES;
	}
	
	WriteSnapshot (rateType);

	GetString  (lfInfo,lf,-1);

	degFCount = Columns (lfInfo["Global Independent"])+
				Columns (lfInfo["Local Independent"])-1+
				freqParamCount;
	
	if (rateType>=0)
	{
		fprintf (SUMMARY_FILE,"\n*** RUNNING MODEL ", Format(rateType+1,0,0), " (",ModelNames[rateType],") ***\n######################################\n");
		fprintf (SUMMARY_FILE,"\n>Done in ",Time(1)-timer, " seconds \n\n", lf);
		fprintf (stdout, "| ");
		if (rateType<9)
		{
			fprintf (stdout," ");
		}
		fprintf (stdout, Format (rateType+1,0,0), ". ", ModelNames[rateType]);
		for (dummy = Abs(ModelNames[rateType])+5; dummy<25; dummy = dummy+1)
		{
			fprintf (stdout," ");
		}
		fprintf (stdout,"| ",Format (res[1][0],14,6)," | ",Format (GetDistributionParameters(psigLevel),13,8)," |  ",
							 Format(degFCount,5,0),"  |\n",tableSeparator);
	}
	else
	{
		fprintf (SUMMARY_FILE,"\n*** RUNNING SINGLE RATE MODEL ***\n#################################\n");
		fprintf (SUMMARY_FILE,"\n>Done in ", Time(1)-timer, " seconds \n\n");
		fprintf (SUMMARY_FILE,lf,"\n\n-----------------------------------\n\ndN/dS = ",c,"\n\n");

		fprintf (stdout, "|  0. Single Rate Model   | ",Format (res[1][0],14,6)," | ",Format (c,13,8)," |  ",Format(degFCount,5,0),"  |\n",
					 		 tableSeparator);
	}
	
	modelLL [rateType] = res[1][0];
						 
	if (MPI_NODE_COUNT>1)
	{
		rateType  = rateType3;
	}

	return fromNode-1;
}
