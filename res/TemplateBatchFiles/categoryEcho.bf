if (categoriesUsed)
{
	ChoiceList  (response,"Rate class Information Dislpay",1,NO_SKIP,
				 "Don't Display","No rate class information is displayed",
				 "Distribution Info","Only the inferred distributions will be shown",
				 "Class Asignments 1","A string of class assignments for sites is displayed for each site from left to right.",
				 "Class Asignments 2","A list of sites, grouped by rate class, is displayed.",
				 "Complete","A table of partial likelihoods is spooled to a file, along with class assignments for sites.");
				 
				 
	if (response>0)
	{
		if (response>3)
		{
			ConstructCategoryMatrix(res,lf,COMPLETE);
			SetDialogPrompt ("Write marginal matrix info to:");
			PRINT_DIGITS = 12;
			fprintf (PROMPT_FOR_FILE,CLEAR_FILE,res);
		}
		else
		{
			ConstructCategoryMatrix(res,lf,SHORT);
			cols = Columns(res);
			GetInformation (distrInfo,c);
			Dvar = Columns(distrInfo);
			Evar = 0.0;
			T = 0.0;
			sampleVar = 0.0;
			for (k=0; k<Dvar; k=k+1)
			{
				T = distrInfo[0][k]*distrInfo[1][k];
				Evar = Evar+T;
				sampleVar = T*distrInfo[0][k]+sampleVar;
			}
			sampleVar = sampleVar-Evar*Evar;
			
			if (categoriesUsed>1)
			{
				fprintf (stdout, "\n", _rateDescriptors[0]);
			}

			fprintf  (stdout,"\n\n------------------------------------------------\n\nSample mean = ",Evar, " (sample variance = ",sampleVar,")\n");
			for (k=0; k<Dvar; k=k+1)
			{
				fprintf (stdout,"\nRate[",Format(k,0,0),"]=",Format(distrInfo[0][k],12,8), " (weight=", 
								  Format(distrInfo[1][k],9,7),")");
			}
			if (categoriesUsed>1)
			{
				GetInformation (distrInfo2,d);
				Dvar = Columns(distrInfo2);
				Evar = 0.0;
				T = 0.0;
				sampleVar = 0.0;
				for (k=0; k<Dvar; k=k+1)
				{
					T = R*distrInfo2[0][k]*distrInfo2[1][k];
					Evar = Evar+T;
					sampleVar = T*R*distrInfo2[0][k]+sampleVar;
				}
				sampleVar = sampleVar-Evar*Evar;

				fprintf (stdout, "\n", _rateDescriptors[1]);
				fprintf  (stdout,"\n\n------------------------------------------------\n\nSample mean = ",Evar, " (sample variance = ",sampleVar,")\n");
				for (k=0; k<Dvar; k=k+1)
				{
					fprintf (stdout,"\nRate[",Format(k,0,0),"]=",Format(R*distrInfo2[0][k],12,8), " (weight=", 
									  Format(distrInfo2[1][k],9,7),")");
				}
			}
			if (response == 2)
			{
				PRINT_DIGITS = 4;
				if (categoriesUsed>1)
				{
					fprintf (stdout,"\n\nClass Assignments (div/mod ",Columns(distrInfo2)," for ",_rateDescriptors[0],"/",_rateDescriptors[1],") :\n\n");
				}
				else
				{
					fprintf (stdout,"\n\nClass Assignments:\n\n");
				}
				fprintf (stdout,"\n");
				PRINT_DIGITS = 3;
				fprintf (stdout,res[0][0]);
				for (n=1;n<cols;n=n+1)
				{
					fprintf (stdout,"\t",res[0][n]);
				}
				_tableHeaders = {1,1};
				_tableHeaders[0] = "Rate Class";
				
				res = Transpose(res);
				
				OpenWindow (CHARTWINDOW,{{"Rate class assignments by site."}
						   {"_tableHeaders"},
						   {"res"},
						   {"Bar Chart"},
						   {"Index"},
						   {_tableHeaders[0]},
						   {"Site Index"},
						   {_tableHeaders[0]},
						   {_tableHeaders[0]},
						   {"0"}},
						   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");
			}
			else
			{
				if (response == 3)
				{
					Dvar = Columns(distrInfo);
					if (categoriesUsed>1)
					{
						Dvar2 = Columns (distrInfo2);
						firstOrNot = {Dvar*Dvar2,1};
						for (k=0; k<Dvar; k=k+1)
						{	
							for (kk = 0; kk<Dvar2; kk=kk+1)
							{
								fprintf (stdout, "\n\nRate Class ",_rateDescriptors[0]," ",Format(k,0,0)," and ",
																   _rateDescriptors[1]," ",Format(kk,0,0),"\n");
																   
								kk2 = k*Dvar2 + kk;
								for (n=0; n<cols; n=n+1)
								{
									if (res[0][n]==kk2)
									{
										if (firstOrNot[kk2]==0)
										{
											firstOrNot[kk2] = 1;
											fprintf (stdout, n);
										}
										else
										{
											fprintf (stdout, ",", n);
										}
									}
								}
							}
						}
					}
					else
					{
						firstOrNot = {Dvar,1};
						for (k=0; k<Dvar; k=k+1)
						{
							fprintf (stdout, "\n\nRate Class ",Format(k,0,0),"\n");
							for (n=0; n<cols; n=n+1)
							{
								if (res[0][n]==k)
								{
									if (firstOrNot[k]==0)
									{
										firstOrNot[k] = 1;
										fprintf (stdout, n);
									}
									else
									{
										fprintf (stdout, ",", n);
									}
								}
							}
						}
					}
					fprintf (stdout, "\n");
				}
			}	
		}
	}
}
