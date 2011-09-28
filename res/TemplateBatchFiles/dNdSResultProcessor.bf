LoadFunctionLibrary("PS_Plotters.bf");

distributionM1 		= {{0}};
distributionM2 		= {{0}};
distributionSynM3   = {{0}};
distributionNSM3    = {{0}};
distributionSynM4   = {{0}};
distributionNSM4    = {{0}};
distributionSynM5   = {{0}};
distributionNSM5    = {{0}};
marginalsM1         = {{0}};
marginalsM2         = {{0}};
marginalsM3         = {{0}};
marginalsM4         = {{0}};
marginalsM5         = {{0}};

modelNamesShort = {{"Constant","Proportional","Nonsynonymous","Dual","Lineage Dual"}};


/*___________________________________________________________________________________________________________*/

function  PromptForMarginals (modelID)
{
	promptString = "Marginal Distribution for "+ modelNamesShort [modelID];
	
	SetDialogPrompt (promptString);
	
	if (modelID==1)
	{
		fscanf (PROMPT_FOR_FILE,"Matrix,Matrix",distributionM1,marginalsM1);
	}
	if (modelID==2)
	{
		fscanf (PROMPT_FOR_FILE,"Matrix,Matrix",distributionM2,marginalsM2);
	}
	if (modelID==3)
	{
		fscanf (PROMPT_FOR_FILE,"Matrix,Matrix,Matrix",distributionSynM3,distributionNSM3,marginalsM3);
	}
	if (modelID==4)
	{
		fscanf (PROMPT_FOR_FILE,"Matrix,Matrix,Matrix",distributionSynM4,distributionNSM4,marginalsM4);
	}
	if (modelID==5)
	{
		fscanf (PROMPT_FOR_FILE,"Matrix,Matrix,Matrix",distributionSynM5,distributionNSM5,marginalsM5);
	}
	
	return 1;
}

/*___________________________________________________________________________________________________________*/

function ComputeRateClasses (modelID, binRates)
{
	dummy = PromptForMarginals (modelID);
	
	if (modelID==1)
	{
		marginalMatrix = marginalsM1;
		distribMatrix  = distributionM1;
		divTerm		   = Columns (distributionM1);
	}
	if (modelID==2)
	{
		marginalMatrix = marginalsM2;
		distribMatrix  = distributionM2;
		divTerm		   = Columns (distributionM2);
	}
	if (modelID==3)
	{
		marginalMatrix = marginalsM3;
		distribMatrixS = distributionSynM3;
		distribMatrixN = distributionNSM3;
		divTerm		   = Columns (distributionNSM3);
	}
	if (modelID==4)
	{
		distribMatrixS = distributionSynM4;
		distribMatrixN = distributionNSM4;
		marginalMatrix = marginalsM4;
		divTerm		   = Columns (distributionNSM4);
	}
		
	if ((Rows(marginalMatrix)==0)||(Columns(marginalMatrix)==0))
	{
		fprintf (stdout, "\n***ERROR:Invalid marginal likelihood matrix***\n");
	}
	else
	{
		if (modelID<=2)
		{
			rateAssignmentMatrix = {Columns(marginalMatrix),1};
			for (counter1 = 0; counter1 < Columns (marginalMatrix); counter1=counter1+1)
			{
				columnMax = 0.0;
				maxColumn = 0;
				
				for (counter2 = 0; counter2 < Rows (marginalMatrix); counter2 = counter2+1)
				{
					tempVal = marginalMatrix[counter2][counter1];
					if (tempVal>columnMax)
					{
						columnMax  = tempVal;
						maxColumn  = counter2; 
	 				}
				}
				rateAssignmentMatrix[counter1][0] = maxColumn;
			}
			if (binRates)
			{
				titleMatrix = {divTerm,2};
				for (counter2 = 0; counter2 < Rows (rateAssignmentMatrix); counter2 = counter2 + 1)
				{
				 	counter1 = rateAssignmentMatrix[counter2];
				 	titleMatrix[counter1][1] = titleMatrix[counter1][1] + 1;
				}
				
				for (counter1=0; counter1 < divTerm; counter1 = counter1+1)
				{
					titleMatrix[counter1][0] = distribMatrix[0][counter1];
				}

				rateAssignmentMatrix = titleMatrix;
				titleMatrix	= {{"Rate Class","Rate","Sites In Class"}};
			}
			else
			{
				titleMatrix	= {{"Site Index","Rate Class"}};
			}
		}
		else
		{
			rateAssignmentMatrix = {Columns(marginalMatrix),2};
			for (counter1 = 0; counter1 < Columns (marginalMatrix); counter1=counter1+1)
			{
				columnMax = 0.0;
				maxColumn = 0;
				
				for (counter2 = 0; counter2 < Rows (marginalMatrix); counter2 = counter2+1)
				{
					tempVal = marginalMatrix[counter2][counter1];
					if (tempVal>columnMax)
					{
						columnMax  = tempVal;
						maxColumn  = counter2; 
	 				}
				}
				rateAssignmentMatrix[counter1][0] = maxColumn$divTerm;
				rateAssignmentMatrix[counter1][1] = maxColumn%divTerm;
			}
			if (binRates)
			{
				if (Columns (distribMatrixS) >= Columns (distribMatrixN))
				{
					titleMatrix = {Rows(marginalMatrix)/divTerm,divTerm+2};
					for (counter2 = 0; counter2 < Rows (rateAssignmentMatrix); counter2 = counter2 + 1)
					{
					 	counter1 = rateAssignmentMatrix[counter2][0];
					 	counter3 = rateAssignmentMatrix[counter2][1];
					 	
					 	titleMatrix[counter1][counter3] = titleMatrix[counter1][counter3] + 1;
					}
					rateAssignmentMatrix = titleMatrix;
					titleMatrix	= {1,Columns(rateAssignmentMatrix)+1};
					titleMatrix[0] = "Syn Class";
					for (counter1 = 1; counter1 <= Columns(distribMatrixN); counter1 = counter1+1)
					{
						titleMatrix[counter1] = "NS Rate = "+ distribMatrixN[0][counter1-1];
					}
					titleMatrix[counter1] = "Non-Syn Rates";
					for (counter2 = 0; counter2 < Columns (distribMatrixN); counter2 = counter2+1)
					{
						rateAssignmentMatrix [counter2][counter1-1] = distribMatrixN[0][counter2];
					}
					
					titleMatrix[counter1+1] = "Syn Rates";
					for (counter2 = 0; counter2 < Columns (distribMatrixS); counter2 = counter2+1)
					{
						rateAssignmentMatrix [counter2][counter1] = distribMatrixS[0][counter2];
					}
				}
				else
				{
					divTerm = Columns (distribMatrixS);
					
					titleMatrix = {Rows(marginalMatrix)/divTerm,divTerm+2};
					for (counter2 = 0; counter2 < Rows (rateAssignmentMatrix); counter2 = counter2 + 1)
					{
					 	counter1 = rateAssignmentMatrix[counter2][0];
					 	counter3 = rateAssignmentMatrix[counter2][1];
					 	
					 	titleMatrix[counter3][counter1] = titleMatrix[counter3][counter1] + 1;
					}
					
					rateAssignmentMatrix = titleMatrix;
					titleMatrix	= {1,Columns(rateAssignmentMatrix)+1};
					titleMatrix[0] = "Non-Syn Class";
					for (counter1 = 1; counter1 <= Columns(distribMatrixS); counter1 = counter1+1)
					{
						titleMatrix[counter1] = "Syn Rate = "+ distribMatrixS[0][counter1-1];
					}
					titleMatrix[counter1] = "Syn Rates";
					for (counter2 = 0; counter2 < Columns (distribMatrixS); counter2 = counter2+1)
					{
						rateAssignmentMatrix [counter2][counter1-1] = distribMatrixS[0][counter2];
					}
					
					titleMatrix[counter1+1] = "Non-Syn Rates";
					for (counter2 = 0; counter2 < Columns (distribMatrixN); counter2 = counter2+1)
					{
						rateAssignmentMatrix [counter2][counter1] = distribMatrixN[0][counter2];
					}				
				}
			}
			else
			{
				titleMatrix	= {{"Site Index","Syn Rate Class","Non-syn Rate Class"}};
			}
		}
		
		if (outputChoice == 0)
		{
			dummy = PrintASCIITable (rateAssignmentMatrix, titleMatrix);
		}
		if (outputChoice == 1)
		{
			dummy = PrintTableToFile (rateAssignmentMatrix, titleMatrix,1);
		}
		if (outputChoice == 2)
		{
			labelMatrix = {1,Columns(titleMatrix)-1};
			for (counter1 = 0; counter1 < Columns (labelMatrix); counter1 = counter1+1)
			{
				labelMatrix[0][counter1] = titleMatrix[0][counter1+1];
			}
			
			if (modelID<=2)
			{
				if (binRates)
				{
					promptString = "Class Histogram for " + modelNamesShort[modelID];
					OpenWindow (CHARTWINDOW,{{promptString}
										   {"labelMatrix"}
										   {"rateAssignmentMatrix"}
										   {"Bar Chart"}
										   {labelMatrix[0]}
										   {labelMatrix[1]}
										   {titleMatrix[0]}
										   {""}
										   {titleMatrix[1]}
										   {"0"}},
										   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");
				}				
				else
				{
					promptString = "Class Assignments for " + modelNamesShort[modelID];
					OpenWindow (CHARTWINDOW,{{promptString}
										   {"labelMatrix"}
										   {"rateAssignmentMatrix"}
										   {"Bar Chart"}
										   {"Index"}
										   {labelMatrix[0]}
										   {titleMatrix[0]}
										   {""}
										   {titleMatrix[1]}
										   {"0"}},
										   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");
				}
			}
			else
			{
				if (binRates)
				{
					promptString = "";
					for (counter1 = 1; counter1 < Columns(rateAssignmentMatrix)-2; counter1 = counter1+1)
					{
						promptString = promptString+titleMatrix[counter1]+";";
					}
					promptString = promptString+titleMatrix[Columns(rateAssignmentMatrix)-2];
					scaleString = ""+(Columns (rateAssignmentMatrix)-1)+";"+(Columns (rateAssignmentMatrix));
					OpenWindow (CHARTWINDOW,{{"Class Histogram for " + modelNamesShort[modelID]}
											   {"labelMatrix"},
											   {"rateAssignmentMatrix"},
											   {"3D Histogram"},
											   {"Index"},
											   {promptString},
											   {"A"},
											   {"B"},
											   {"C"},
											   {"3"},
											   {""}
											   {scaleString}},
											   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");
				}
				else
				{
					OpenWindow (CHARTWINDOW,{{"Class Assignments for " + modelNamesShort[modelID]}
											   {"labelMatrix"},
											   {"rateAssignmentMatrix"},
											   {"Contrast Bars"},
											   {"Index"},
											   {"Syn Rate Class;Non-syn Rate Class"},
											   {titleMatrix[0]},
											   {titleMatrix[2]},
											   {titleMatrix[1]},
											   {"0"}},
											   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");
				}
			}
			
		}
	}
	return 1;
}

/*___________________________________________________________________________________________________________*/

function ComputePosteriorRates (modelID)
{
	dummy = PromptForMarginals (modelID);
	
	if (modelID==1)
	{
		marginalMatrix = marginalsM1;
		divTerm		   = Columns (distributionM1);
		distribMX	   = distributionM1;
	}
	if (modelID==2)
	{
		marginalMatrix = marginalsM2;
		divTerm		   = Columns (distributionM2);
		distribMX	   = distributionM2;
	}
	if (modelID==3)
	{
		marginalMatrix = marginalsM3;
		divTerm		   = Columns (distributionNSM3);
		distribMXS	   = distributionSynM3;
		distribMXN	   = distributionNSM3;
	}
	if (modelID==4)
	{
		marginalMatrix = marginalsM4;
		divTerm		   = Columns (distributionNSM4);
		distribMXS	   = distributionSynM4;
		distribMXN	   = distributionNSM4;
	}
		
	rMMX  = Rows(marginalMatrix);
	
	if ((rMMX==0)||(Columns(marginalMatrix)==0))
	{
		fprintf (stdout, "\n***ERROR:Invalid marginal likelihood matrix***\n");
	}
	else
	{
		if (modelID<=2)
		{
			rateAssignmentMatrix = {Columns(marginalMatrix),1};
			for (counter1 = 0; counter1 < Columns (marginalMatrix); counter1=counter1+1)
			{
				columnSum = 0;				
				for (counter2 = 0; counter2 < rMMX; counter2 = counter2+1)
				{
					tempVal    = marginalMatrix [counter2][counter1]*distribMX[1][counter2];
					columnSum += tempVal;
					marginalMatrix [counter2][counter1] = tempVal;
				}
				
				tempVal   = 0;
				
				for (counter2 = 0; counter2 < rMMX; counter2 = counter2+1)
				{
					tempVal   +=  distribMX[0][counter2]*marginalMatrix [counter2][counter1]/columnSum;
				}
				rateAssignmentMatrix[counter1][0] = tempVal;
			}
			if (modelID==1)
			{
				titleMatrix	= {{"Site Index","E[S|i]"}};
			}
			else
			{
				titleMatrix	= {{"Site Index","E[N|i]"}};			
			}
		}
		else
		{
			rateAssignmentMatrix = {Columns(marginalMatrix),3};
			
			for (counter1 = 0; counter1 < Columns (marginalMatrix); counter1+=1)
			{
				columnSum = 0;		
				
				for (counter2 = 0; counter2 < rMMX; counter2 += 1)
				{
					entryWeight    =  marginalMatrix [counter2][counter1] * distribMXS[1][counter2$divTerm] * distribMXN[1][counter2%divTerm];
					columnSum      += entryWeight;
					marginalMatrix [counter2][counter1] = entryWeight;
				}
				
				dnds_post = 0;
				for (counter2 = 0; counter2 < rMMX; counter2 += 1)
				{
					marginalMatrix [counter2][counter1] = marginalMatrix [counter2][counter1]/columnSum;
					dnds_post      += distribMXN[0][counter2%divTerm] / distribMXS[0][counter2$divTerm] * marginalMatrix [counter2][counter1];
				}
				
				rateAssignmentMatrix[counter1][2] = dnds_post;
				
				tempVal = 0;
				
				for (counter2 = 0; counter2 < rMMX; counter2 += divTerm)
				{
					columnSum = 0;
					for (counter3 = counter2; counter3 < counter2+divTerm; counter3 += 1)
					{
						columnSum += marginalMatrix [counter3][counter1];
					}
					tempVal += columnSum * distribMXS[0][counter2$divTerm];
				}
				
				rateAssignmentMatrix[counter1][0] = tempVal;
				
				tempVal = 0;
				
				for (counter2 = 0; counter2 < divTerm; counter2 += 1)
				{
					columnSum = 0;
					for (counter3 = counter2; counter3 < rMMX; counter3 += divTerm)
					{
						columnSum +=  marginalMatrix [counter3][counter1];
					}
					tempVal += columnSum * distribMXN[0][counter2];
				}
				

				rateAssignmentMatrix[counter1][1] = tempVal;
			}
			titleMatrix	= {{"Site Index","E[S|i]","E[N|i]","E[omega|i]"}};			
		}
		
		if (outputChoice == 0)
		{
			PrintASCIITable (rateAssignmentMatrix, titleMatrix);
		}
		if (outputChoice == 1)
		{
			PrintTableToFile (rateAssignmentMatrix, titleMatrix,1);
		}
		if (outputChoice == 2)
		{
			labelMatrix = {1,Columns(titleMatrix)-1};
			for (counter1 = 0; counter1 < Columns (labelMatrix); counter1 = counter1+1)
			{
				labelMatrix[0][counter1] = titleMatrix[0][counter1+1];
			}
			
			promptString = "Expected Posterior Rates for " + modelNamesShort[modelID];
			if (modelID<=2)
			{
				OpenWindow (CHARTWINDOW,{{promptString}
										   {"labelMatrix"},
										   {"rateAssignmentMatrix"},
										   {"Bar Chart"},
										   {"Index"},
										   {labelMatrix[0]},
										   {titleMatrix[0]},
										   {""},
										   {titleMatrix[1]},
										   {"0"}},
										   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");
			}
			else
			{
				promptString2 = labelMatrix[0]+";"+labelMatrix[1];
				OpenWindow (CHARTWINDOW,{{promptString}
										   {"labelMatrix"},
										   {"rateAssignmentMatrix"},
										   {"Contrast Bars"},
										   {"Index"},
										   {promptString2},
										   {titleMatrix[0]},
										   {titleMatrix[2]},
										   {titleMatrix[1]},
										   {"0"}},
										   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");
			}
			
		}
	}
	return 1;
}

/*___________________________________________________________________________________________________________*/

function ComputePositiveSelection (modelID, sitesOnly)
{
	dummy = PromptForMarginals (modelID);
	
	if (modelID==2)
	{
		marginalMatrix = marginalsM2;
		divTerm		   = Columns (distributionM2);
		distribMX	   = distributionM2;
	}
	if (modelID==3)
	{
		marginalMatrix = marginalsM3;
		divTerm		   = Columns (distributionNSM3);
		distribMXS	   = distributionSynM3;
		distribMXN	   = distributionNSM3;
	}
	if (modelID==4)
	{
		marginalMatrix = marginalsM4;
		divTerm		   = Columns (distributionNSM4);
		distribMXS	   = distributionSynM4;
		distribMXN	   = distributionNSM4;
	}
	
	if (sitesOnly > .5)
	{
		if (psChoice)
		{
			thresh = 0;
			while (thresh<=1.)
			{
				fprintf (stdout, "Select the threshold for the Bayes factor on P{N/S>1} (>1) for a site to be marked as positively selected:");
				fscanf  (stdin,"Number",thresh);
			}
			thresh = Log (thresh);
		}
		else
		{
			thresh = -1;
			while ((thresh<0)||(thresh>1))
			{
				fprintf (stdout, "Select the threshold for P{N/S>1} for a site to be marked as positively selected:");
				fscanf  (stdin,"Number",thresh);
			}
		}
	}
		
	rMMX  = Rows(marginalMatrix);
	
	if ((rMMX==0)||(Columns(marginalMatrix)==0))
	{
		fprintf (stdout, "\n***ERROR:Invalid marginal likelihood matrix***\n");
	}
	else
	{
			
		if (sitesOnly > .5)
		{
			if (psChoice==0)
			{
				titleMatrix	= {{"Counter","Site Index","P[NS/S>1|s]"}};
			}			
			else
			{
				titleMatrix	= {{"Counter","Site Index","Log[BF{NS/S>1|s}]"}};
			}
		}		
		else
		{
			if (psChoice==0)
			{
				titleMatrix	= {{"Codon","Log[BF{NS/S>1|s}]","Log[BF{NS/S<=1|s}]"}};
			}			
			else
			{
				titleMatrix	= {{"Codon","P[NS/S>1|s]","P[NS/S<=1|s]"}};
			}
		}
			
		if (modelID==2)
		{
			distribMXR = distribMX;
		}
		else
		{
			D1 = Columns(distribMXN);
			D2 = Columns(distribMXS);
				
			distribMXR = {2,D1*D2};
			
			for (k=0; k<D2; k=k+1)
			{
				E = k*D1;
				for (k2 = 0; k2<D1; k2=k2+1)
				{
					distribMXR [0][E+k2] = distribMXN[0][k2]/distribMXS[0][k];
					distribMXR [1][E+k2] = distribMXN[1][k2]*distribMXS[1][k];
				}
			}
		}
		
		
		counter3 	 = 0;
		ratesOverOne = {Columns(distribMXR),1};
		
		for (counter1=0; counter1<Columns(distribMXR); counter1=counter1+1)
		{
			if (distribMXR [0][counter1]>1)
			{
				ratesOverOne[counter3][0] = counter1;
				counter3 = counter3 + 1;
			}
		}
		
		if (counter3==0)
		{
			fprintf (stdout, "\nThere are no rates with N/S > 1 and thus P{N/S>1|i} = 0 for all sites \n\n");
			return 1;
		}
				
		if (psChoice)
		{
			priorOdds = 0;
			for (k=0; k<Columns(distribMXR); k=k+1)
			{
				if (distribMXR [0][k]>=1)
				{
					priorOdds = priorOdds + distribMXR [1][k];
				}
			}
			priorOdds = priorOdds/(1-priorOdds);
		}
		
		if (sitesOnly>.5)
		{
			E = 0;
			rateAssignmentMatrix = {Columns(marginalMatrix),2};
			for (counter1 = 0; counter1 < Columns (marginalMatrix); counter1=counter1+1)
			{
				columnSum = 0;				
				for (counter2 = 0; counter2 < rMMX; counter2 = counter2+1)
				{
					tempVal   = marginalMatrix [counter2][counter1];
					columnSum = columnSum + tempVal*distribMXR[1][counter2];
					marginalMatrix [counter2][counter1] = tempVal * distribMXR[1][counter2];
				}
				
				tempVal   = 0;
				
				for (counter2 = 0; counter2 < counter3; counter2 = counter2+1)
				{
					indexVal  = ratesOverOne [counter2][0];
					tempVal   = tempVal + marginalMatrix [indexVal][counter1];
				}
				tempVal = tempVal/columnSum;
				
				if (psChoice)
				{
					tempVal = Log(tempVal/(1-tempVal)/priorOdds);
				}
				
				if (tempVal>=thresh)
				{
					rateAssignmentMatrix[E][0] = counter1+1;
					rateAssignmentMatrix[E][1] = tempVal;
					E = E+1;
				}
			}
			
			if (E==0)
			{
				fprintf (stdout, "\nThere are no sites with P{N/S>1|i} >= ",thresh, "\n\n");
				return 1;
			}
			
			ratesOverOne = {E,2};
			for (counter1 = 0; counter1 < E; counter1 = counter1+1)
			{
				ratesOverOne[counter1][0] = rateAssignmentMatrix[counter1][0];
				ratesOverOne[counter1][1] = rateAssignmentMatrix[counter1][1];
			}		
			rateAssignmentMatrix = ratesOverOne;
			ratesOverOne = 0;
		}		
		else
		{
			rateAssignmentMatrix = {Columns(marginalMatrix),2};
			for (counter1 = 0; counter1 < Columns (marginalMatrix); counter1=counter1+1)
			{
				columnSum = 0;				
				for (counter2 = 0; counter2 < rMMX; counter2 = counter2+1)
				{
					tempVal   = marginalMatrix [counter2][counter1];
					columnSum = columnSum + tempVal*distribMXR[1][counter2];
					marginalMatrix [counter2][counter1] = tempVal * distribMXR[1][counter2];
				}
				
				tempVal   = 0;
				
				for (counter2 = 0; counter2 < counter3; counter2 = counter2+1)
				{
					indexVal  = ratesOverOne [counter2][0];
					tempVal   = tempVal + marginalMatrix [indexVal][counter1];
				}
				if (psChoice)
				{
					tempVal = tempVal/columnSum;
                    if (tempVal == 1)
                    {
                        logOdds = 1e26;
                    }
                    else
                    {
                        logOdds = Log(tempVal/(1-tempVal)/priorOdds);
                    }
					rateAssignmentMatrix[counter1][0] = Max(logOdds,0);
					rateAssignmentMatrix[counter1][1] = Max(-logOdds,0);
				}
				else
				{
					rateAssignmentMatrix[counter1][0] = tempVal/columnSum;
					rateAssignmentMatrix[counter1][1] = 1-tempVal/columnSum;
				}
			}
		}
		
		

		if (outputChoice == 0)
		{
			dummy = PrintASCIITable (rateAssignmentMatrix, titleMatrix);
		}
		if (outputChoice == 1)
		{
			dummy = PrintTableToFile (rateAssignmentMatrix, titleMatrix,1);
		}
		if (outputChoice == 2)
		{
			labelMatrix = {1,Columns(titleMatrix)-1};
			for (counter1 = 0; counter1 < Columns (labelMatrix); counter1 = counter1+1)
			{
				labelMatrix[0][counter1] = titleMatrix[0][counter1+1];
			}
			
			if (sitesOnly>.5)
			{
				promptString = "Positively Selected Sites for " + modelNamesShort[modelID];
				OpenWindow (CHARTWINDOW,{{promptString}
										   {"labelMatrix"},
										   {"rateAssignmentMatrix"},
										   {"Bar Chart"},
										   {labelMatrix[0]},
										   {labelMatrix[1]},
										   {titleMatrix[0]},
										   {""},
										   {titleMatrix[1]},
										   {"0"}},
										   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");			
			}			
			else
			{
				if (psChoice)
				{
					promptString = "Log of Bayes Factor for P{N/S>1} for " + modelNamesShort[modelID];				
				}
				else
				{
					promptString = "P{N/S>1|i} for " + modelNamesShort[modelID];
				}
				OpenWindow (CHARTWINDOW,{{promptString}
										   {"labelMatrix"},
										   {"rateAssignmentMatrix"},
										   {"Contrast Bars"},
										   {"Index"},
										   {labelMatrix[0]+";"+labelMatrix[1]},
										   {titleMatrix[0]},
										   {titleMatrix[2]},
										   {titleMatrix[1]},
										   {"0"}},
										   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");			
			}
			
			
		}
	}
	return 1;
}


/*___________________________________________________________________________________________________________*/

function ComputePRatio (modelID)
{
	dummy = PromptForMarginals (modelID);
	
	if (modelID==3)
	{
		marginalMatrix = marginalsM3;
		divTerm		   = Columns (distributionNSM3);
		distribMXS	   = distributionSynM3;
		distribMXN	   = distributionNSM3;
	}
	if (modelID==4)
	{
		marginalMatrix = marginalsM4;
		divTerm		   = Columns (distributionNSM4);
		distribMXS	   = distributionSynM4;
		distribMXN	   = distributionNSM4;
	}
		
	rMMX  = Rows(marginalMatrix);
	
	if ((rMMX==0)||(Columns(marginalMatrix)==0))
	{
		fprintf (stdout, "\n***ERROR:Invalid marginal likelihood matrix***\n");
	}
	else
	{
			
		titleMatrix	= {{"Site Index","P[NS/S>1|i]"}};
		
		D1 = Columns(distribMXN);
		D2 = Columns(distribMXS);
			
		distribMXR = {2,D1*D2};
		distribMXD = {2,D1*D2};
		
		for (k=0; k<D2; k=k+1)
		{
			E = k*D1;
			for (k2 = 0; k2<D1; k2=k2+1)
			{
				distribMXR [0][E+k2] = distribMXN[0][k2]/distribMXS[0][k];
				distribMXR [1][E+k2] = distribMXN[1][k2]*distribMXS[1][k];
				distribMXD [0][E+k2] = distribMXN[0][k2]-distribMXS[0][k];
				distribMXD [1][E+k2] = distribMXN[1][k2]*distribMXS[1][k];
			}
		}
		
		titleMatrix	= {{"Site Index","E[dN/dS|i]","E[dN-dS|i]"}};

		rateAssignmentMatrix = {Columns(marginalMatrix),2};
		
		for (counter1 = 0; counter1 < Columns (marginalMatrix); counter1=counter1+1)
		{
			columnSum  = 0;				

			for (counter2 = 0; counter2 < rMMX; counter2 = counter2+1)
			{
				tempVal   = marginalMatrix [counter2][counter1];
				columnSum = columnSum + tempVal*distribMXR[1][counter2];
				marginalMatrix [counter2][counter1] = tempVal * distribMXR[1][counter2];
			}
			
			tempVal    = 0;
			tempVal2   = 0;
			
			for (counter2 = 0; counter2 < rMMX; counter2 = counter2+1)
			{
				tempVal   = tempVal  + distribMXR[0][counter2]*marginalMatrix [counter2][counter1];
				tempVal2  = tempVal2 + distribMXD[0][counter2]*marginalMatrix [counter2][counter1];
			}
			
			rateAssignmentMatrix[counter1][0] = tempVal/columnSum;
			rateAssignmentMatrix[counter1][1] = tempVal2/columnSum;
		}

		if (outputChoice == 0)
		{
			dummy = PrintASCIITable (rateAssignmentMatrix, titleMatrix);
		}
		if (outputChoice == 1)
		{
			dummy = PrintTableToFile (rateAssignmentMatrix, titleMatrix,1);
		}
		if (outputChoice == 2)
		{
			labelMatrix = {1,Columns(titleMatrix)-1};
			for (counter1 = 0; counter1 < Columns (labelMatrix); counter1 = counter1+1)
			{
				labelMatrix[0][counter1] = titleMatrix[0][counter1+1];
			}
			
			promptString = "Posterior dN dS Ratio and Difference for " + modelNamesShort[modelID];
			OpenWindow (CHARTWINDOW,{{promptString}
									   {"labelMatrix"},
									   {"rateAssignmentMatrix"},
									   {"Bar Chart"},
									   {"Index"},
									   {labelMatrix[0]},
									   {titleMatrix[0]},
									   {""},
									   {titleMatrix[1]},
									   {"0"}},
									   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");			
			
		}
	}
	return 1;
}


/*___________________________________________________________________________________________________________*/

function ComputeCDF (modelID)
{
	dummy = PromptForMarginals (modelID);
		
	if (modelID==1)
	{
		divTerm		   = Columns (distributionM1);
		distribMX	   = distributionM1;
	}
	if (modelID==2)
	{
		divTerm		   = Columns (distributionM2);
		distribMX	   = distributionM2;
	}
	if (modelID==3)
	{
		divTerm		   = Columns (distributionNSM3);
		distribMXS	   = distributionSynM3;
		distribMXN	   = distributionNSM3;
	}
	if (modelID==4)
	{
		divTerm		   = Columns (distributionNSM4);
		distribMXS	   = distributionSynM4;
		distribMXN	   = distributionNSM4;
	}
	
	if (modelID<=2)
	{
		rateAssignmentMatrix = {Columns(distribMX)+1,2};
		for (counter1 = 1; counter1 < Rows (rateAssignmentMatrix); counter1=counter1+1)
		{
			rateAssignmentMatrix[counter1][0] = distribMX[0][counter1-1];
			rateAssignmentMatrix[counter1][1] = rateAssignmentMatrix[counter1-1][1]+distribMX[1][counter1-1];
		}
		titleMatrix	= {{"Class","Rate","CDF"}};
	}
	else
	{
		rateAssignmentMatrixS = {Columns(distribMXS)+1,2};
		for (counter1 = 1; counter1 < Rows (rateAssignmentMatrixS); counter1=counter1+1)
		{
			rateAssignmentMatrixS[counter1][0] = distribMXS[0][counter1-1];
			rateAssignmentMatrixS[counter1][1] = rateAssignmentMatrixS[counter1-1][1]+distribMXS[1][counter1-1];
		}
		rateAssignmentMatrixN = {Columns(distribMXN)+1,2};
		for (counter1 = 1; counter1 < Rows (rateAssignmentMatrixN); counter1=counter1+1)
		{
			rateAssignmentMatrixN[counter1][0] = distribMXN[0][counter1-1];
			rateAssignmentMatrixN[counter1][1] = rateAssignmentMatrixN[counter1-1][1]+distribMXN[1][counter1-1];
		}
		
		D1 = Columns(distribMXN);
		D2 = Columns(distribMXS);
			
		rateMatrixValues = {D1*D2,3};
		epsilon			 = Exp(-4);
		
		distribMXR = {2,D1*D2};
		
		for (k=0; k<D1; k=k+1)
		{
			E = k*D2;
			for (k2 = 0; k2<D2; k2=k2+1)
			{
				distribMXR [0][E+k2] = distribMXN[0][k]/distribMXS[0][k2];
				distribMXR [1][E+k2] = distribMXN[1][k]*distribMXS[1][k2];
				rateMatrixValues [E+k2][0] = Min(4,Log (epsilon + distribMXS[0][k2]));
				rateMatrixValues [E+k2][1] = Min(4,Log (epsilon + distribMXN[0][k]));
				rateMatrixValues [E+k2][2] = distribMXR [1][E+k2];
			}
		}
		
		E = D1*D2;
		
		done = 0;
		
		while (!done)
		{
			done = 1;
			for (k=1; k<E; k=k+1)
			{
				if (distribMXR [0][k] < distribMXR[0][k-1])
				{
					D = distribMXR [0][k];
					distribMXR [0][k] = distribMXR [0][k-1];
					distribMXR [0][k-1] = D;
					D = distribMXR [1][k];
					distribMXR [1][k] = distribMXR [1][k-1];
					distribMXR [1][k-1] = D;
					done = 0;
				}
				else
				{
					if (distribMXR [0][k] == distribMXR[0][k-1])
					{
						distribMXR[1][k-1] = distribMXR[1][k-1]+distribMXR[1][k];
						for (k2=k+1;k2<E;k2=k2+1)
						{
							distribMXR[0][k2-1] = distribMXR[0][k2];
							distribMXR[1][k2-1] = distribMXR[1][k2];
						}
						E = E-1;
					}
				}
			}
		}

		if (E<D1*D2)
		{
			ratioInfoSwap = {2,E};
			for (k=0; k<E; k=k+1)
			{
				ratioInfoSwap[0][k] = distribMXR[0][k];
				ratioInfoSwap[1][k] = distribMXR[1][k];
			}
			distribMXR = ratioInfoSwap;
			ratioInfoSwap = 0;
		}
		
		rateAssignmentMatrixR = {Columns(distribMXR)+1,2};
		RatioCDFString = "";
		RatioCDFString * 128;
		for (counter1 = 1; counter1 < Rows (rateAssignmentMatrixR); counter1=counter1+1)
		{
			rateAssignmentMatrixR[counter1][0] = distribMXR[0][counter1-1];
			rateAssignmentMatrixR[counter1][1] = rateAssignmentMatrixR[counter1-1][1]+distribMXR[1][counter1-1];
			RatioCDFString * ("+"+distribMXR[1][counter1-1]+"*(_x_>="+distribMXR[0][counter1-1]+")");
		}
		RatioCDFString * 0;
		RatioCDFString = RatioCDFString[1][Abs(RatioCDFString)-1];
		
		titleMatrix	= {{"Class","Rate","CDF"}};
	}
	
	if (outputChoice == 0)
	{
		if (modelID<=2)
		{
			dummy = PrintASCIITable  (rateAssignmentMatrix, titleMatrix);
		}
		else
		{
			fprintf (stdout, "\nSynonymous CDF\n\n");
			dummy = PrintASCIITable  (rateAssignmentMatrixS, titleMatrix);
			fprintf (stdout, "\nNonsynonymous CDF\n\n");
			dummy = PrintASCIITable  (rateAssignmentMatrixN, titleMatrix);
			fprintf (stdout, "\nNS/Syn Ratio CDF\n\n");
			dummy = PrintASCIITable  (rateAssignmentMatrixR, titleMatrix);
		}
	}
	if (outputChoice == 1)
	{
		if (modelID<=2)
		{
			dummy = PrintTableToFile (rateAssignmentMatrix, titleMatrix,1);
		}
		else
		{
			SetDialogPrompt ("Export tab separated data to:");
			fprintf (PROMPT_FOR_FILE,CLEAR_FILE, "\nSynonymous CDF\n\n");
			dummy = PrintTableToFile  (rateAssignmentMatrixS, titleMatrix,0);
			fprintf (LAST_FILE_PATH, "\nNonsynonymous CDF\n\n");
			dummy = PrintTableToFile (rateAssignmentMatrixN, titleMatrix,0);
			fprintf (LAST_FILE_PATH, "\nNS/Syn Ratio CDF\n\n");
			dummy = PrintTableToFile  (rateAssignmentMatrixR, titleMatrix,0);
		}
	}
	if (outputChoice == 2)
	{
		labelMatrix = {1,Columns(titleMatrix)-1};
		for (counter1 = 0; counter1 < Columns (labelMatrix); counter1 = counter1+1)
		{
			labelMatrix[0][counter1] = titleMatrix[0][counter1+1];
		}
		
		if (modelID<=2)
		{
			promptString = "Discrete Rate Distribution for " + modelNamesShort[modelID];
			OpenWindow (CHARTWINDOW,{{promptString}
									   {"labelMatrix"},
									   {"rateAssignmentMatrix"},
									   {"Step Plot"},
									   {labelMatrix[0]},
									   {labelMatrix[1]},
									   {labelMatrix[0]},
									   {""},
									   {labelMatrix[1]},
									   {"0"}},
									   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");
		}
		else
		{
			promptString = "Discrete Syn Rate Distribution for " + modelNamesShort[modelID];
			OpenWindow (CHARTWINDOW,{{promptString}
									   {"labelMatrix"},
									   {"rateAssignmentMatrixS"},
									   {"Step Plot"},
									   {labelMatrix[0]},
									   {labelMatrix[1]},
									   {labelMatrix[0]},
									   {""},
									   {labelMatrix[1]},
									   {"0"}},
									   "SCREEN_WIDTH/2-20;SCREEN_HEIGHT/2-25;30;50");
									   
			promptString = "Discrete Non-Syn Rate Distribution for " + modelNamesShort[modelID];
			OpenWindow (CHARTWINDOW,{{promptString}
									   {"labelMatrix"},
									   {"rateAssignmentMatrixN"},
									   {"Step Plot"},
									   {labelMatrix[0]},
									   {labelMatrix[1]},
									   {labelMatrix[0]},
									   {""},
									   {labelMatrix[1]},
									   {"0"}},
									   "SCREEN_WIDTH/2-20;SCREEN_HEIGHT/2-25;SCREEN_WIDTH/2;50");
									   
			promptString = "Discrete Non-Syn/Syn Rate Distribution for " + modelNamesShort[modelID];
			OpenWindow (CHARTWINDOW,{{promptString}
									   {"labelMatrix"},
									   {"rateAssignmentMatrixR"},
									   {"Step Plot"},
									   {labelMatrix[0]},
									   {labelMatrix[1]},
									   {labelMatrix[0]},
									   {""},
									   {labelMatrix[1]},
									   {"0"}},
									   "SCREEN_WIDTH-60;SCREEN_HEIGHT/2-25;30;50+SCREEN_HEIGHT/2");
		}
		
	}
	if (outputChoice == 3)
	{
		psCode = ScaledDensityPlot ("rateMatrixValues", {{-4,4}{-4,4}}, "Courier", {{400,400,14,40}},
										"_dNdSDensityPlot",
										{{"","log (alpha)", "log (beta)"}}, 1, 1);

		SetDialogPrompt ("Write the PostScript plot to:");
		fprintf (PROMPT_FOR_FILE, CLEAR_FILE, psCode);
	}
	return 1;
}

/*___________________________________________________________________________________________________________*/

function PadString (padLength,padChar)
{
	for (padCounter=0;padCounter<padLength;padCounter=padCounter+1)
	{
		fprintf (stdout,padChar);
	}
	return padLength;
}

/*___________________________________________________________________________________________________________*/

function	PrintASCIITable (dataMatrix, titleMatrix)
{
	columnWidths = {1,Columns(titleMatrix)};
	for (counter1=0; counter1<Columns(titleMatrix); counter1 = counter1+1)
	{
		counter2 = Abs (titleMatrix[0][counter1])+2;
		if (counter2<12)
		{
			counter2 = 12;
		}
		columnWidths[0][counter1] = counter2;
	}
	fprintf (stdout, "\n");
	for (counter2 = 0; counter2 < Columns (titleMatrix); counter2 = counter2+1)
	{
		fprintf (stdout,"+-");
		dummy = PadString (columnWidths[0][counter2],"-");
		fprintf (stdout,"-");
	}
	fprintf (stdout,"+\n| ");
	
	for (counter1=0; counter1<Columns(titleMatrix); counter1 = counter1+1)
	{
		fprintf (stdout, titleMatrix[counter1]);
		dummy = PadString (columnWidths[0][counter1]-Abs(titleMatrix[counter1])," ");
		fprintf (stdout, " | ");
	}
	
	fprintf (stdout, "\n");
	
	for (counter1=-1; counter1<Rows(dataMatrix); counter1 = counter1 + 1)
	{
		if (counter1>=0)
		{
			fprintf (stdout,"| ");
			fprintf (stdout,Format(counter1+1,columnWidths[0][0],0));
			for (counter2 = 1; counter2 < Columns (titleMatrix); counter2 = counter2+1)
			{
				fprintf (stdout," | ");
				fprintf (stdout,Format(dataMatrix[counter1][counter2-1],columnWidths[0][counter2],-1));
			}
			fprintf (stdout," ");
			fprintf (stdout, "|\n");
		}
		for (counter2 = 0; counter2 < Columns (titleMatrix); counter2 = counter2+1)
		{
			fprintf (stdout,"+-");
			dummy = PadString (columnWidths[0][counter2],"-");
			fprintf (stdout,"-");
		}
		fprintf (stdout, "+\n");
	}
	
	return 1;
}

/*___________________________________________________________________________________________________________*/

function	PrintTableToFile (dataMatrix, titleMatrix, promptOrNot)
{
	SetDialogPrompt ("Export tab separated data to:");
	
	if (promptOrNot)
	{
		fprintf (PROMPT_FOR_FILE, CLEAR_FILE);
	}
	
	fprintf (LAST_FILE_PATH, Join("\t",titleMatrix), "\n");

	for (counter1=0; counter1<Rows(dataMatrix); counter1 += 1)
	{
		fprintf (LAST_FILE_PATH,counter1+1);
		for (counter2 = 1; counter2 < Columns (titleMatrix); counter2 = counter2+1)
		{
			fprintf (LAST_FILE_PATH,"\t",dataMatrix[counter1][counter2-1]);
		}
		fprintf (LAST_FILE_PATH,"\n");
	}
	
	return 1;
}

/*___________________________________________________________________________________________________________*/

ChoiceList (actionChoice, 					  "Available Actions",1,SKIP_NONE,
			"Rate classes",		 			  "Compute the rate class assignments at each site for the selected model.",
			"Rate histogram",	 			  "Compute the rate class histogram for the selected model.",
			"Posterior rates",	 			  "Compute the posterior distribution of expected rate for the selected model.",
			"Positive Selection",			  "Compute the posterior probability  of {NS Rate > Syn Rate} at each site for the selected model.",			
			"Distributions",	 			  "Display discrete rate distribution functions.",
			"Posterior Ratio",	 			  "Display expected posterior dN/dS ratio for dual rate variation models.",
			"Find Positively Selected Sites", "Find    positively selected sites.");
			
if (actionChoice<0)
{
	return;
}

if (actionChoice==3 || actionChoice==6)
{
	ChoiceList (psChoice, 					  	  "Detection method",1,SKIP_NONE,
				"Absolute threshold",		 	  "Use the absolute value for posterior P{N/S>=1}.",
				"Bayes Factor",	 		  		  "Use the Bayes factor: Posterior Odds P{N/S>=1}/Prior Odds P{N/S>=1}.");
				
	if (psChoice<0)
	{
		return;
	}
}

if (actionChoice!=3 && actionChoice!=5 && actionChoice!=6)
{
	ChoiceList (modelChoice,"Rate Variation Models",1,SKIP_NONE,
				"Proportional","Proportional Variable Rates Model: dS and dN vary along the sequence, but dN = R*dS for every site",
				"Nonsynonymous","Non-synonymous Variable Rates Model: dS = 1 for every site, while dN is drawn from a given distribution.",
				"Dual","Dual Variable Rates Model: dS and dN are drawn from a bivariate distribution (independent or correlated components).",
				"Lineage Dual","Lineage+Dual Variable Rates Model:  dS and dN are drawn from a bivariate distribution (independent or correlated components), plus each lineage has an adjustment factor for the E[dN]/E[dS]."
			    );
	
	if (modelChoice<0)
	{
		return;
	}
}
else
{
	if ((actionChoice == 3)||(actionChoice == 6))
	{
		ChoiceList (modelChoice,"Rate Variation Models",1,SKIP_NONE,
					"Nonsynonymous","Non-synonymous Variable Rates Model: dS = 1 for every site, while dN is drawn from a given distribution.",
					"Dual","Dual Variable Rates Model: dS and dN are drawn from a bivariate distribution (independent or correlated components). Recommened model.",
					"Lineage Dual","Lineage+Dual Variable Rates Model:  dS and dN are drawn from a bivariate distribution (independent or correlated components), plus each lineage has an adjustment factor for the E[dN]/E[dS]."
				    );
	
		if (modelChoice<0)
		{
			return;
		}
		
		modelChoice = modelChoice+1;
	}
	else
	{
		ChoiceList (modelChoice,"Rate Variation Models",1,SKIP_NONE,
					"Dual","Dual Variable Rates Model: dS and dN are drawn from a bivariate distribution (independent or correlated components). Recommened model.",
					"Lineage Dual","Lineage+Dual Variable Rates Model:  dS and dN are drawn from a bivariate distribution (independent or correlated components), plus each lineage has an adjustment factor for the E[dN]/E[dS]."
			    );
	
		if (modelChoice<0)
		{
			return;
		}
		
		modelChoice = modelChoice+2;
	}	
}

if (actionChoice == 4 && modelChoice == 2)
{
	ChoiceList (outputChoice, "Output Options",1,SKIP_NONE,
				"ASCII Table",	 "Output is printed to the console as an ASCII table.",
				"Export to File","Output is spooled to a tab separated file.",
				"Chart","A HYPHY chart window is displayed (GUI versions only).",
				"PostScript","Render a PostScript density plot into a file.");
}
else
{
	ChoiceList (outputChoice, "Output Options",1,SKIP_NONE,
				"ASCII Table",	 "Output is printed to the console as an ASCII table.",
				"Export to File","Output is spooled to a tab separated file.",
				"Chart","A HYPHY chart window is displayed (GUI versions only).");
}

if (outputChoice<0)
{
	return;
}
			
if (actionChoice < 2)
{
	ComputeRateClasses 		(modelChoice+1,actionChoice);
	return;
}

if (actionChoice == 2)
{
	ComputePosteriorRates 	(modelChoice+1);
	return;
}

if (actionChoice == 3)
{
	ComputePositiveSelection(modelChoice+1,0);
	return;
}

if (actionChoice == 4)
{
	ComputeCDF				(modelChoice+1);
	return;
}

if (actionChoice == 5)
{
	ComputePRatio			(modelChoice+1);
	return;
}

if (actionChoice == 6)
{
	ComputePositiveSelection(modelChoice+1,1);
	return;
}
