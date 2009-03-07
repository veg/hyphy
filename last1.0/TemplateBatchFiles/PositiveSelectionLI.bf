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
			  "3 Normal"}};
			       
MAXIMUM_ITERATIONS_PER_VARIABLE = 2000;
OPTIMIZATION_PRECISION = 0.001;


/*-------------------------------------------------------------------------------*/

function  setElement (h,v,cc)
{	
	mSpecMatrix[h][v]=cc+1;
	mSpecMatrix[v][h]=cc+1;
	return 1;
}

/*-------------------------------------------------------------------------------*/

function  reportModelType (modelType, prefix, customString)
{	
	modelType = modelType$2;
	if (modelType == 0)
	{
		ms = "MG94";
	}
	else
	{
		if (modelType == 1)
		{
			ms = "GY94";
		}
		else
		{
			ms = "MG94 custom " + customString; 
		}
	}
	fprintf (summaryFile, "Using ", ms, " model for ", prefix,"\n");
	fprintf (stdout, "Using ", ms, " model for ", prefix,"\n");
	return 1;
}

/*-------------------------------------------------------------------------------*/

function  defineCustomModel (dummy)
{
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
	mSpecMatrix = {{*,1,1,1}{1,*,1,1}{1,1,*,1}{1,1,1,*}};

	elementAssignments = {{0,0,0,0,0,0}};
	paramCount = 1;
	r=setElement (0,1,0);

	for (j=1; j<6; j=j+1)
	{
		for (i=0; i<j; i=i+1)
		{
			if (modelDesc[j]==modelDesc[i])
			{
				elementAssignments[j]=elementAssignments[i];
				if (j<3)
				{
					r=setElement (0,j+1,elementAssignments[i]);
				}
				else
				{
					if (j<5)
					{
						r=setElement (1,j-1,elementAssignments[i]);
					}
					else
					{
						r=setElement (2,3,elementAssignments[i]);
					}
				}
				break;
			}
		}
		if (i==j)
		{
			if (j<3)
			{
				r=setElement (0,j+1,paramCount);
			}
			else
			{
				if (j<5)
				{
					r=setElement (1,j-1,paramCount);
				}
				else
				{
					r=setElement (2,3,paramCount);
				}
			}
			elementAssignments[j] = paramCount;
			paramCount = paramCount+1;
		}
	}
	return 0;
}

/*-------------------------------------------------------------------------------*/

function defineNumberOfCategories (rateType)
{
	resp = 0;
	if (rateType > 3)
	{
		while (resp<2)
		{
			fprintf (stdout, "Choose the number of classes in the discretized distribution (>=2):");
			fscanf  (stdin, "Number", resp);
		}
	}
	return resp;
}
			  
/*-------------------------------------------------------------------------------*/

function SetWDistributionLeaves (rateType)
{
	resp = defineNumberOfCategories (rateType);
	
	if (rateType == 0)
	{
		global P_leaves = .5;
		P_leaves:<1;
		categFreqMatrix_leaves = {{P_leaves,1-P_leaves}};
		categRateMatrix_leaves = {{0,1}};
		category category_leaves = (2, categFreqMatrix_leaves , MEAN, ,categRateMatrix_leaves, 0, 1e25);
	}
	else
	{
		if (rateType == 1)
		{
			global P1_leaves = 1/3;
			global P2_leaves = 0;
			
			P1_leaves:<1;
			P2_leaves:<1;
			
			global W_leaves = .5;
			categFreqMatrix_leaves = {{P1_leaves,(1-P1_leaves)*P2_leaves, (1-P1_leaves)*(1-P2_leaves)}} ;
			categRateMatrix_leaves = {{0,1,W_leaves}};
			category category_leaves = (3, categFreqMatrix_leaves , MEAN, ,categRateMatrix_leaves, 0, 1e25);
		}		
		else
		{
			if (rateType == 2)
			{
				global P1_leaves = 1/3;
				global P2_leaves = .5;
				P1_leaves:<1;
				P2_leaves:<1;
				global W1_leaves = .25;
				global R_1_leaves = 4;
				global R_2_leaves = 3;
				R_1_leaves:<1;
				R_2_leaves:>1;
				categFreqMatrix_leaves = {{P1_leaves,(1-P1_leaves)*P2_leaves, (1-P1_leaves)*(1-P2_leaves)}} ;
				categRateMatrix_leaves = {{W1_leaves*R_1_leaves,W1_leaves,W1_leaves*R_2_leaves}};
				category category_leaves = (3, categFreqMatrix_leaves , MEAN, ,categRateMatrix_leaves, 0, 1e25);				
			}
			else
			{
				if (rateType == 3)
				{
					global P1_leaves = 1/5;
					global P2_leaves = 1/4;
					global P3_leaves = 1/3;
					global P4_leaves = 1/2;
					
					P1_leaves:<1;
					P2_leaves:<1;
					P3_leaves:<1;
					P4_leaves:<1;
					
					categFreqMatrix_leaves = {{P1_leaves,
										(1-P1_leaves)P2_leaves,
										(1-P1_leaves)(1-P2_leaves)*P3_leaves,
										(1-P1_leaves)(1-P2_leaves)(1-P3_leaves)P4_leaves,
										(1-P1_leaves)(1-P2_leaves)(1-P3_leaves)(1-P4_leaves)}} ;
					categRateMatrix_leaves = {{0,1/3,2/3,1,3}};
					category category_leaves = (5, categFreqMatrix_leaves , MEAN, ,categRateMatrix_leaves, 0, 1e25);				
				}
				else
				{
					if (rateType == 4)
					{
						global alpha_leaves = .5;
						global beta_leaves = 1;
						alpha_leaves:>0.01;alpha_leaves:<100;
						beta_leaves:>0.01;
						beta_leaves:<200;
						category category_leaves = (resp, EQUAL, MEAN, GammaDist(_x_,alpha_leaves,beta_leaves), CGammaDist(_x_,alpha_leaves,beta_leaves), 0 , 
				 			 		  1e25,CGammaDist(_x_,alpha_leaves+1,beta_leaves)*alpha_leaves/beta_leaves);
					}
					else
					{
						if (rateType == 5)
						{
							global alpha_leaves = .5;
							global beta_leaves  =  1;
							global alpha2_leaves=  .75;
							global P_leaves	 = .5; 
							alpha_leaves:>0.01;alpha_leaves:<100;
							beta_leaves:>0.01;
							beta_leaves:<200;
							P_leaves:<1;
							alpha2_leaves:>0.01;alpha2_leaves:<100;
							category category_leaves = (resp, EQUAL, MEAN, P_leaves*GammaDist(_x_,alpha_leaves,beta_leaves) + (1-P_leaves)*GammaDist(_x_,alpha2_leaves,alpha2_leaves)
														   , P_leaves*CGammaDist(_x_,alpha_leaves,beta_leaves) + (1-P_leaves)*CGammaDist(_x_,alpha2_leaves,alpha2_leaves), 
														   0 , 1e25,
														   P_leaves*CGammaDist(_x_,alpha_leaves+1,beta_leaves)*alpha_leaves/beta_leaves + (1-P_leaves)*CGammaDist(_x_,alpha2_leaves+1,alpha2_leaves));
						}
						else
						{
							if (rateType == 6)
							{
								global betaP_leaves = 1;
								global betaQ_leaves = 1;
								betaP_leaves:>0.05;betaP_leaves:<85;
								betaQ_leaves:>0.05;betaQ_leaves:<85;
								category category_leaves = (resp, EQUAL, MEAN, _x_^(betaP_leaves-1)*(1-_x_)^(betaQ_leaves-1)/Beta(betaP_leaves,betaQ_leaves), IBeta(_x_,betaP_leaves,betaQ_leaves), 0 , 
						 			 		  1,IBeta(_x_,betaP_leaves+1,betaQ_leaves)*betaP_leaves/(betaP_leaves+betaQ_leaves));
							}
							else
							{
								if (rateType == 7)
								{
									global W_leaves = 2;
									W_leaves:>1;
									global P_leaves	 = 1-1/(resp+1);
									global betaP_leaves = 1;
									global betaQ_leaves = 2;
									betaP_leaves:>0.05;
									betaQ_leaves:>0.05;
									betaP_leaves:<85;
									betaQ_leaves:<85;
									P_leaves:>0.0000001;
									P_leaves:<0.9999999;
									categFreqMatrix_leaves = {resp+1,1};
									for (k=0; k<resp; k=k+1)
									{
										categFreqMatrix_leaves[k]:=P_leaves/resp__;
									}
									categFreqMatrix_leaves[resp]:=(1-P_leaves);
									category category_leaves = (resp+1, categFreqMatrix_leaves, MEAN, 
													P_leaves*_x_^(betaP_leaves-1)*(1-Min(_x_,1))^(betaQ_leaves-1)/Beta(betaP_leaves,betaQ_leaves)+W_leaves-W_leaves, 
													P_leaves*IBeta(Min(_x_,1),betaP_leaves,betaQ_leaves)+(1-P_leaves)*(_x_>=W_leaves), 
													0,1e25,
													P_leaves*IBeta(Min(_x_,1),betaP_leaves+1,betaQ_leaves)*betaP_leaves/(betaP_leaves+betaQ_leaves)+(1-P_leaves)*W_leaves*(_x_>=W_leaves));
								}
								else
								{
									if (rateType == 8)
									{
										global P_leaves	 = .5;
										global betaP_leaves = 1;
										global betaQ_leaves = 2;
										betaP_leaves:>0.05;betaP_leaves:<85;
										betaQ_leaves:>0.05;betaQ_leaves:<85;
										global alpha_leaves = .5;
										global beta_leaves  = 1;
										alpha_leaves:>0.01;alpha_leaves:<100;
										beta_leaves:>0.01;										
										beta_leaves:<200;
										P_leaves:<1;
										category category_leaves = (resp, EQUAL, MEAN, 
															P_leaves*_x_^(betaP_leaves-1)*(1-Min(_x_,1))^(betaQ_leaves-1)/Beta(betaP_leaves,betaQ_leaves)+(1-P_leaves)*GammaDist(_x_,alpha_leaves,beta_leaves), 
															P_leaves*IBeta(Min(_x_,1),betaP_leaves,betaQ_leaves)+(1-P_leaves)*CGammaDist(_x_,alpha_leaves,beta_leaves), 
															0,1e25,
															P_leaves*betaP_leaves/(betaP_leaves+betaQ_leaves)*IBeta(Min(_x_,1),betaP_leaves+1,betaQ_leaves)+(1-P_leaves)*alpha_leaves/beta_leaves*CGammaDist(_x_,alpha_leaves+1,beta_leaves));
									}	
									else
									{
										if (rateType == 9)
										{
											global P_leaves	 = .5;
											P_leaves:<1;
											global betaP_leaves = 1;
											betaP_leaves:>0.05;betaP_leaves:<85;
											global betaQ_leaves = 2;
											betaQ_leaves:>0.05;betaQ_leaves:<85;
											global alpha_leaves = .5;
											alpha_leaves:>0.01;alpha_leaves:<100;
											global beta_leaves  = 1;
											beta_leaves:>0.01;beta_leaves:<500;
											category category_leaves = (resp, EQUAL, MEAN, 
																P_leaves*_x_^(betaP_leaves-1)*(1-Min(_x_,1))^(betaQ_leaves-1)/Beta(betaP_leaves,betaQ_leaves)+(1-P_leaves)*(_x_>1)*GammaDist(Max(1e-20,_x_-1),alpha_leaves,beta_leaves), 
																P_leaves*IBeta(Min(_x_,1),betaP_leaves,betaQ_leaves)+(1-P_leaves)*CGammaDist(Max(_x_-1,0),alpha_leaves,beta_leaves), 
																0,1e25,
																P_leaves*betaP_leaves/(betaP_leaves+betaQ_leaves)*IBeta(Min(_x_,1),betaP_leaves+1,betaQ_leaves)+
																		(1-P_leaves)*(alpha_leaves/beta_leaves*CGammaDist(Max(0,_x_-1),alpha_leaves+1,beta_leaves)+CGammaDist(Max(0,_x_-1),alpha_leaves,beta_leaves)));
										}				
										else
										{
											if (rateType == 10)
											{
												global P_leaves	 = .5;
												global betaP_leaves = 1;
												global betaQ_leaves = 2;
												betaP_leaves:>0.05;
												betaQ_leaves:>0.05;
												betaP_leaves:<85;
												betaQ_leaves:<85;
												global mu_leaves = 3;
												global sigma_leaves  = .01;
												sigma_leaves:>0.0001;
												sqrt2pi = Sqrt(8*Arctan(1));
												P_leaves:<1;

												category category_leaves = (resp, EQUAL, MEAN, 
																P_leaves*_x_^(betaP_leaves-1)*(1-Min(_x_,1))^(betaQ_leaves-1)/Beta(betaP_leaves,betaQ_leaves)+
																	(1-P_leaves)*(_x_>=1)*Exp(-(_x_-mu_leaves)(_x_-mu_leaves)/(2*sigma_leaves*sigma_leaves))/(sqrt2pi__*sigma_leaves)/ZCDF((mu_leaves-1)/sigma_leaves), 
																P_leaves*IBeta(Min(_x_,1),betaP_leaves,betaQ_leaves)+(1-P_leaves)*(_x_>=1)*(1-ZCDF((mu_leaves-_x_)/sigma_leaves)/ZCDF((mu_leaves-1)/sigma_leaves)), 
																0,1e25,
																P_leaves*betaP_leaves/(betaP_leaves+betaQ_leaves)*IBeta(Min(_x_,1),betaP_leaves+1,betaQ_leaves)+
																(1-P_leaves)*(_x_>=1)*(mu_leaves*(1-ZCDF((1-mu_leaves)/sigma_leaves)-ZCDF((mu_leaves-_x_)/sigma_leaves))+
																sigma_leaves*(Exp((mu_leaves-1)(1-mu_leaves)/(2*sigma_leaves*sigma_leaves))-Exp((_x_-mu_leaves)(mu_leaves-_x_)/(2*sigma_leaves*sigma_leaves)))/sqrt2pi__)/ZCDF((mu_leaves-1)/sigma_leaves));
											}				
											else
											{
												if (rateType == 11)
												{
													global P_leaves	 = 1/3;
													global P1_leaves    = .5;

													global mu_leaves = 3;
													global sigma_leaves  = .5;
													sigma_leaves:>0.0001;
													global sigma1_leaves  = 1;
													sigma1_leaves:>0.0001;

													sqrt2pi = Sqrt(8*Arctan(1));
													P_leaves:<1;
													P1_leaves:<1;
													
													categFreqMatrix_leaves = {resp+1,1};
													for (k=1; k<=resp; k=k+1)
													{
														categFreqMatrix_leaves[k]:=(1-P_leaves)/resp__;
													}
													categFreqMatrix_leaves[0]:=P_leaves;

													category category_leaves = (resp+1, categFreqMatrix_leaves, MEAN,
																	(1-P_leaves)((1-P1_leaves)*Exp(-(_x_-mu_leaves)(_x_-mu_leaves)/(2*sigma1_leaves*sigma1_leaves))/(sqrt2pi__*sigma1_leaves)/ZCDF(mu_leaves/sigma1_leaves)+
																			  P1_leaves*Exp(-(_x_-1)(_x_-1)/(2*sigma_leaves*sigma_leaves))/(sqrt2pi__*sigma_leaves)/ZCDF(1/sigma_leaves)), 
																	P_leaves+(1-P_leaves)(_x_>1e-20)((1-P1_leaves)(1-ZCDF((mu_leaves-_x_)/sigma1_leaves)/ZCDF(mu_leaves/sigma1_leaves))+
																						P1_leaves*(1-ZCDF((1-_x_)/sigma_leaves)/ZCDF(1/sigma_leaves))), 
																	0,1e25,
																	(1-P_leaves)((1-P1_leaves)(mu_leaves*(1-ZCDF(-mu_leaves/sigma1_leaves)-ZCDF((mu_leaves-_x_)/sigma1_leaves))+
																	sigma1_leaves*(Exp(-mu_leaves*mu_leaves/(2*sigma1_leaves*sigma1_leaves))-Exp((_x_-mu_leaves)(mu_leaves-_x_)/(2*sigma1_leaves*sigma1_leaves)))/sqrt2pi__)/ZCDF(mu_leaves/sigma1_leaves)+
																	P_leaves(1-ZCDF(-1/sigma_leaves)-ZCDF((1-_x_)/sigma_leaves)+
																	sigma_leaves*(Exp(-1/(2*sigma_leaves*sigma_leaves))-Exp((_x_-1)(1-_x_)/(2*sigma_leaves*sigma_leaves)))/sqrt2pi__)/ZCDF(1/sigma_leaves))
																 );
												}
												else		
												{
													if (rateType == 12)
													{
														global P_leaves	 = 1/3;
														global P1_leaves    = .5;

														global mu_leaves = 3;
														global sigma_leaves  = .25;
														global sigma1_leaves = .5;
														global sigma2_leaves = 1;
														sigma_leaves:>0.0001;
														sigma1_leaves:>0.0001;
														sigma2_leaves:>0.0001;

														sqrt2pi = Sqrt(8*Arctan(1));
														P_leaves:<1;
														P1_leaves:<1;

														category category_leaves = (resp, EQUAL , MEAN,
																		2*P_leaves*Exp(-_x_^2/(2*sigma_leaves*sigma_leaves))+
																		(1-P_leaves)((1-P1_leaves)*Exp((_x_-mu_leaves)(mu_leaves-_x_)/(2*sigma2_leaves*sigma2_leaves))/(sqrt2pi__*sigma2_leaves)/ZCDF(mu_leaves/sigma2_leaves)+
																			  P1_leaves*Exp((1-_x_)(_x_-1)/(2*sigma1_leaves*sigma1_leaves))/(sqrt2pi__*sigma1_leaves)/ZCDF(1/sigma1_leaves)), 
																		P_leaves*(1-2*ZCDF(-_x_/sigma_leaves))+
																		(1-P_leaves)((1-P1_leaves)(1-ZCDF((mu_leaves-_x_)/sigma2_leaves)/ZCDF(mu_leaves/sigma2_leaves))+
																			   P1_leaves*(1-ZCDF((1-_x_)/sigma1_leaves)/ZCDF(1/sigma1_leaves))), 
																		0,1e25,
																		2*P_leaves*sigma_leaves*(1-Exp(-_x_*_x_/(2*sigma_leaves*sigma_leaves)))/sqrt2pi__+
																		(1-P_leaves)((1-P1_leaves)(mu_leaves*(1-ZCDF(-mu_leaves/sigma2_leaves)-ZCDF((mu_leaves-_x_)/sigma2_leaves))+
																		sigma2_leaves*(Exp(-mu_leaves*mu_leaves/(2*sigma2_leaves*sigma2_leaves))-Exp((_x_-mu_leaves)(mu_leaves-_x_)/(2*sigma2_leaves*sigma2_leaves)))/sqrt2pi__)/ZCDF(mu_leaves/sigma2_leaves)+
																		P1_leaves(1-ZCDF(-1/sigma1_leaves)-ZCDF((1-_x_)/sigma1_leaves)+
																		sigma1_leaves*(Exp(-1/(2*sigma1_leaves*sigma1_leaves))-Exp((_x_-1)(1-_x_)/(2*sigma1_leaves*sigma1_leaves)))/sqrt2pi__)/ZCDF(mu_leaves/sigma1_leaves))
																		);
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
	}		
	return 0;
}

/*-------------------------------------------------------------------------------*/

function ConstrainRateDistributionLeaves (rateType)
{
	if (rateType == 0)
	{
		P_leaves:=P_inodes;
	}
	else
	{
		if (rateType == 1)
		{
			P1_leaves := P1_inodes;
			P2_leaves := P2_inodes;
		}		
		else
		{
			if (rateType == 2)
			{
				P1_leaves  := P1_inodes;
				P2_leaves  := P2_inodes;
				W1_leaves  := W1_inodes;
				R_1_leaves := R_1_inodes;
				R_2_leaves := R_2_inodes;
			}
			else
			{
				if (rateType == 3)
				{
					P1_leaves := P1_inodes;
					P2_leaves := P2_inodes;
					P3_leaves := P3_inodes;
					P4_leaves := P4_inodes;
				}
				else
				{
					if (rateType == 4)
					{
						alpha_leaves := alpha_inodes;
						beta_leaves  := beta_inodes;
						
					}
					else
					{
						if (rateType == 5)
						{
							alpha_leaves :=  alpha_inodes;
							beta_leaves  :=  beta_inodes;
							alpha2_leaves:=  alpha2_inodes;
							P_leaves	 :=  P_leaves; 
						}
						else
						{
							if (rateType == 6)
							{
								betaP_leaves := betaP_inodes;
								betaQ_leaves := betaQ_inodes;
							}
							else
							{
								if (rateType == 7)
								{
									W_leaves     := W_inodes;
									P_leaves	 := P_inodes;
									betaP_leaves := betaP_inodes;
									betaQ_leaves := betaQ_inodes;
								}
								else
								{
									if (rateType == 8)
									{
										P_leaves	 := P_inodes;
										betaP_leaves := betaP_inodes;
										betaQ_leaves := betaQ_inodes;
										alpha_leaves := alpha_inodes;
										beta_leaves  := beta_inodes;
									}	
									else
									{
										if (rateType == 9)
										{
											P_leaves	 := P_inodes;
											betaP_leaves := betaP_inodes;
											betaQ_leaves := betaQ_inodes;
											alpha_leaves := alpha_inodes;
											beta_leaves  := beta_inodes;
										}				
										else
										{
											if (rateType == 10)
											{
												P_leaves	 := P_inodes;
												betaP_leaves := betaP_inodes;
												betaQ_leaves := betaQ_inodes;
												mu_leaves =  := mu_inodes;
												sigma_leaves := sigma_inodes;
											}				
											else
											{
												if (rateType == 11)
												{
													P_leaves	 := P_inodes;
													P1_leaves    := P1_inodes;
													mu_leaves    := mu_inodes;
													sigma_leaves := sigma_inodes;
													sigma1_leaves:= sigma1_inodes;
												}
												else		
												{
													if (rateType == 12)
													{
														P_leaves	  := P_inodes;
														P1_leaves     := P1_inodes;
														mu_leaves     := mu_inodes;
														sigma_leaves  := sigma_inodes;
														sigma1_leaves := sigma1_inodes;
														sigma2_leaves := sigma2_inodes;
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
	}		
	return 0;
}

/*-------------------------------------------------------------------------------*/

function SetWDistributionINodes (rateType)
{
	resp = defineNumberOfCategories (rateType);
	
	if (rateType == 0)
	{
		global P_inodes = .5;
		P_inodes:<1;
		categFreqMatrix_inodes = {{P_inodes,1-P_inodes}};
		categRateMatrix_inodes = {{0,1}};
		category category_inodes = (2, categFreqMatrix_inodes , MEAN, ,categRateMatrix_inodes, 0, 1e25);
	}
	else
	{
		if (rateType == 1)
		{
			global P1_inodes = 1/3;
			global P2_inodes = 0;
			
			P1_inodes:<1;
			P2_inodes:<1;
			
			global W_inodes = .5;
			categFreqMatrix_inodes = {{P1_inodes,(1-P1_inodes)*P2_inodes, (1-P1_inodes)*(1-P2_inodes)}} ;
			categRateMatrix_inodes = {{0,1,W_inodes}};
			category category_inodes = (3, categFreqMatrix_inodes , MEAN, ,categRateMatrix_inodes, 0, 1e25);
		}		
		else
		{
			if (rateType == 2)
			{
				global P1_inodes = 1/3;
				global P2_inodes = .5;
				P1_inodes:<1;
				P2_inodes:<1;
				global W1_inodes = .25;
				global R_1_inodes = 4;
				global R_2_inodes = 3;
				R_1_inodes:<1;
				R_2_inodes:>1;
				categFreqMatrix_inodes = {{P1_inodes,(1-P1_inodes)*P2_inodes, (1-P1_inodes)*(1-P2_inodes)}} ;
				categRateMatrix_inodes = {{W1_inodes*R_1_inodes,W1_inodes,W1_inodes*R_2_inodes}};
				category category_inodes = (3, categFreqMatrix_inodes , MEAN, ,categRateMatrix_inodes, 0, 1e25);				
			}
			else
			{
				if (rateType == 3)
				{
					global P1_inodes = 1/5;
					global P2_inodes = 1/4;
					global P3_inodes = 1/3;
					global P4_inodes = 1/2;
					
					P1_inodes:<1;
					P2_inodes:<1;
					P3_inodes:<1;
					P4_inodes:<1;
					
					categFreqMatrix_inodes = {{P1_inodes,
										(1-P1_inodes)P2_inodes,
										(1-P1_inodes)(1-P2_inodes)*P3_inodes,
										(1-P1_inodes)(1-P2_inodes)(1-P3_inodes)P4_inodes,
										(1-P1_inodes)(1-P2_inodes)(1-P3_inodes)(1-P4_inodes)}} ;
					categRateMatrix_inodes = {{0,1/3,2/3,1,3}};
					category category_inodes = (5, categFreqMatrix_inodes , MEAN, ,categRateMatrix_inodes, 0, 1e25);				
				}
				else
				{
					if (rateType == 4)
					{
						global alpha_inodes = .5;
						global beta_inodes = 1;
						alpha_inodes:>0.01;alpha_inodes:<100;
						beta_inodes:>0.01;
						beta_inodes:<200;
						category category_inodes = (resp, EQUAL, MEAN, GammaDist(_x_,alpha_inodes,beta_inodes), CGammaDist(_x_,alpha_inodes,beta_inodes), 0 , 
				 			 		  1e25,CGammaDist(_x_,alpha_inodes+1,beta_inodes)*alpha_inodes/beta_inodes);
					}
					else
					{
						if (rateType == 5)
						{
							global alpha_inodes = .5;
							global beta_inodes  =  1;
							global alpha2_inodes=  .75;
							global P_inodes	 = .5; 
							alpha_inodes:>0.01;alpha_inodes:<100;
							beta_inodes:>0.01;
							beta_inodes:<200;
							P_inodes:<1;
							alpha2_inodes:>0.01;alpha2_inodes:<100;
							category category_inodes = (resp, EQUAL, MEAN, P_inodes*GammaDist(_x_,alpha_inodes,beta_inodes) + (1-P_inodes)*GammaDist(_x_,alpha2_inodes,alpha2_inodes)
														   , P_inodes*CGammaDist(_x_,alpha_inodes,beta_inodes) + (1-P_inodes)*CGammaDist(_x_,alpha2_inodes,alpha2_inodes), 
														   0 , 1e25,
														   P_inodes*CGammaDist(_x_,alpha_inodes+1,beta_inodes)*alpha_inodes/beta_inodes + (1-P_inodes)*CGammaDist(_x_,alpha2_inodes+1,alpha2_inodes));
						}
						else
						{
							if (rateType == 6)
							{
								global betaP_inodes = 1;
								global betaQ_inodes = 1;
								betaP_inodes:>0.05;betaP_inodes:<85;
								betaQ_inodes:>0.05;betaQ_inodes:<85;
								category category_inodes = (resp, EQUAL, MEAN, _x_^(betaP_inodes-1)*(1-_x_)^(betaQ_inodes-1)/Beta(betaP_inodes,betaQ_inodes), IBeta(_x_,betaP_inodes,betaQ_inodes), 0 , 
						 			 		  1,IBeta(_x_,betaP_inodes+1,betaQ_inodes)*betaP_inodes/(betaP_inodes+betaQ_inodes));
							}
							else
							{
								if (rateType == 7)
								{
									global W_inodes = 2;
									W_inodes:>1;
									global P_inodes	 = 1-1/(resp+1);
									global betaP_inodes = 1;
									global betaQ_inodes = 2;
									betaP_inodes:>0.05;
									betaQ_inodes:>0.05;
									betaP_inodes:<85;
									betaQ_inodes:<85;
									P_inodes:>0.0000001;
									P_inodes:<0.9999999;
									categFreqMatrix_inodes = {resp+1,1};
									for (k=0; k<resp; k=k+1)
									{
										categFreqMatrix_inodes[k]:=P_inodes/resp__;
									}
									categFreqMatrix_inodes[resp]:=(1-P_inodes);
									category category_inodes = (resp+1, categFreqMatrix_inodes, MEAN, 
													P_inodes*_x_^(betaP_inodes-1)*(1-Min(_x_,1))^(betaQ_inodes-1)/Beta(betaP_inodes,betaQ_inodes)+W_inodes-W_inodes, 
													P_inodes*IBeta(Min(_x_,1),betaP_inodes,betaQ_inodes)+(1-P_inodes)*(_x_>=W_inodes), 
													0,1e25,
													P_inodes*IBeta(Min(_x_,1),betaP_inodes+1,betaQ_inodes)*betaP_inodes/(betaP_inodes+betaQ_inodes)+(1-P_inodes)*W_inodes*(_x_>=W_inodes));
								}
								else
								{
									if (rateType == 8)
									{
										global P_inodes	 = .5;
										global betaP_inodes = 1;
										global betaQ_inodes = 2;
										betaP_inodes:>0.05;betaP_inodes:<85;
										betaQ_inodes:>0.05;betaQ_inodes:<85;
										global alpha_inodes = .5;
										global beta_inodes  = 1;
										alpha_inodes:>0.01;alpha_inodes:<100;
										beta_inodes:>0.01;										
										beta_inodes:<200;
										P_inodes:<1;
										category category_inodes = (resp, EQUAL, MEAN, 
															P_inodes*_x_^(betaP_inodes-1)*(1-Min(_x_,1))^(betaQ_inodes-1)/Beta(betaP_inodes,betaQ_inodes)+(1-P_inodes)*GammaDist(_x_,alpha_inodes,beta_inodes), 
															P_inodes*IBeta(Min(_x_,1),betaP_inodes,betaQ_inodes)+(1-P_inodes)*CGammaDist(_x_,alpha_inodes,beta_inodes), 
															0,1e25,
															P_inodes*betaP_inodes/(betaP_inodes+betaQ_inodes)*IBeta(Min(_x_,1),betaP_inodes+1,betaQ_inodes)+(1-P_inodes)*alpha_inodes/beta_inodes*CGammaDist(_x_,alpha_inodes+1,beta_inodes));
									}	
									else
									{
										if (rateType == 9)
										{
											global P_inodes	 = .5;
											P_inodes:<1;
											global betaP_inodes = 1;
											betaP_inodes:>0.05;betaP_inodes:<85;
											global betaQ_inodes = 2;
											betaQ_inodes:>0.05;betaQ_inodes:<85;
											global alpha_inodes = .5;
											alpha_inodes:>0.01;alpha_inodes:<100;
											global beta_inodes  = 1;
											beta_inodes:>0.01;beta_inodes:<500;
											category category_inodes = (resp, EQUAL, MEAN, 
																P_inodes*_x_^(betaP_inodes-1)*(1-Min(_x_,1))^(betaQ_inodes-1)/Beta(betaP_inodes,betaQ_inodes)+(1-P_inodes)*(_x_>1)*GammaDist(Max(1e-20,_x_-1),alpha_inodes,beta_inodes), 
																P_inodes*IBeta(Min(_x_,1),betaP_inodes,betaQ_inodes)+(1-P_inodes)*CGammaDist(Max(_x_-1,0),alpha_inodes,beta_inodes), 
																0,1e25,
																P_inodes*betaP_inodes/(betaP_inodes+betaQ_inodes)*IBeta(Min(_x_,1),betaP_inodes+1,betaQ_inodes)+
																		(1-P_inodes)*(alpha_inodes/beta_inodes*CGammaDist(Max(0,_x_-1),alpha_inodes+1,beta_inodes)+CGammaDist(Max(0,_x_-1),alpha_inodes,beta_inodes)));
										}				
										else
										{
											if (rateType == 10)
											{
												global P_inodes	 = .5;
												global betaP_inodes = 1;
												global betaQ_inodes = 2;
												betaP_inodes:>0.05;
												betaQ_inodes:>0.05;
												betaP_inodes:<85;
												betaQ_inodes:<85;
												global mu_inodes = 3;
												global sigma_inodes  = .01;
												sigma_inodes:>0.0001;
												sqrt2pi = Sqrt(8*Arctan(1));
												P_inodes:<1;

												category category_inodes = (resp, EQUAL, MEAN, 
																P_inodes*_x_^(betaP_inodes-1)*(1-Min(_x_,1))^(betaQ_inodes-1)/Beta(betaP_inodes,betaQ_inodes)+
																	(1-P_inodes)*(_x_>=1)*Exp(-(_x_-mu_inodes)(_x_-mu_inodes)/(2*sigma_inodes*sigma_inodes))/(sqrt2pi__*sigma_inodes)/ZCDF((mu_inodes-1)/sigma_inodes), 
																P_inodes*IBeta(Min(_x_,1),betaP_inodes,betaQ_inodes)+(1-P_inodes)*(_x_>=1)*(1-ZCDF((mu_inodes-_x_)/sigma_inodes)/ZCDF((mu_inodes-1)/sigma_inodes)), 
																0,1e25,
																P_inodes*betaP_inodes/(betaP_inodes+betaQ_inodes)*IBeta(Min(_x_,1),betaP_inodes+1,betaQ_inodes)+
																(1-P_inodes)*(_x_>=1)*(mu_inodes*(1-ZCDF((1-mu_inodes)/sigma_inodes)-ZCDF((mu_inodes-_x_)/sigma_inodes))+
																sigma_inodes*(Exp((mu_inodes-1)(1-mu_inodes)/(2*sigma_inodes*sigma_inodes))-Exp((_x_-mu_inodes)(mu_inodes-_x_)/(2*sigma_inodes*sigma_inodes)))/sqrt2pi__)/ZCDF((mu_inodes-1)/sigma_inodes));
											}				
											else
											{
												if (rateType == 11)
												{
													global P_inodes	 = 1/3;
													global P1_inodes    = .5;

													global mu_inodes = 3;
													global sigma_inodes  = .5;
													sigma_inodes:>0.0001;
													global sigma1_inodes  = 1;
													sigma1_inodes:>0.0001;

													sqrt2pi = Sqrt(8*Arctan(1));
													P_inodes:<1;
													P1_inodes:<1;
													
													categFreqMatrix_inodes = {resp+1,1};
													for (k=1; k<=resp; k=k+1)
													{
														categFreqMatrix_inodes[k]:=(1-P_inodes)/resp__;
													}
													categFreqMatrix_inodes[0]:=P_inodes;

													category category_inodes = (resp+1, categFreqMatrix_inodes, MEAN,
																	(1-P_inodes)((1-P1_inodes)*Exp(-(_x_-mu_inodes)(_x_-mu_inodes)/(2*sigma1_inodes*sigma1_inodes))/(sqrt2pi__*sigma1_inodes)/ZCDF(mu_inodes/sigma1_inodes)+
																			  P1_inodes*Exp(-(_x_-1)(_x_-1)/(2*sigma_inodes*sigma_inodes))/(sqrt2pi__*sigma_inodes)/ZCDF(1/sigma_inodes)), 
																	P_inodes+(1-P_inodes)(_x_>1e-20)((1-P1_inodes)(1-ZCDF((mu_inodes-_x_)/sigma1_inodes)/ZCDF(mu_inodes/sigma1_inodes))+
																						P1_inodes*(1-ZCDF((1-_x_)/sigma_inodes)/ZCDF(1/sigma_inodes))), 
																	0,1e25,
																	(1-P_inodes)((1-P1_inodes)(mu_inodes*(1-ZCDF(-mu_inodes/sigma1_inodes)-ZCDF((mu_inodes-_x_)/sigma1_inodes))+
																	sigma1_inodes*(Exp(-mu_inodes*mu_inodes/(2*sigma1_inodes*sigma1_inodes))-Exp((_x_-mu_inodes)(mu_inodes-_x_)/(2*sigma1_inodes*sigma1_inodes)))/sqrt2pi__)/ZCDF(mu_inodes/sigma1_inodes)+
																	P_inodes(1-ZCDF(-1/sigma_inodes)-ZCDF((1-_x_)/sigma_inodes)+
																	sigma_inodes*(Exp(-1/(2*sigma_inodes*sigma_inodes))-Exp((_x_-1)(1-_x_)/(2*sigma_inodes*sigma_inodes)))/sqrt2pi__)/ZCDF(1/sigma_inodes))
																 );
												}
												else		
												{
													if (rateType == 12)
													{
														global P_inodes	 = 1/3;
														global P1_inodes    = .5;

														global mu_inodes = 3;
														global sigma_inodes  = .25;
														global sigma1_inodes = .5;
														global sigma2_inodes = 1;
														sigma_inodes:>0.0001;
														sigma1_inodes:>0.0001;
														sigma2_inodes:>0.0001;

														sqrt2pi = Sqrt(8*Arctan(1));
														P_inodes:<1;
														P1_inodes:<1;

														category category_inodes = (resp, EQUAL , MEAN,
																		2*P_inodes*Exp(-_x_^2/(2*sigma_inodes*sigma_inodes))+
																		(1-P_inodes)((1-P1_inodes)*Exp((_x_-mu_inodes)(mu_inodes-_x_)/(2*sigma2_inodes*sigma2_inodes))/(sqrt2pi__*sigma2_inodes)/ZCDF(mu_inodes/sigma2_inodes)+
																			  P1_inodes*Exp((1-_x_)(_x_-1)/(2*sigma1_inodes*sigma1_inodes))/(sqrt2pi__*sigma1_inodes)/ZCDF(1/sigma1_inodes)), 
																		P_inodes*(1-2*ZCDF(-_x_/sigma_inodes))+
																		(1-P_inodes)((1-P1_inodes)(1-ZCDF((mu_inodes-_x_)/sigma2_inodes)/ZCDF(mu_inodes/sigma2_inodes))+
																			   P1_inodes*(1-ZCDF((1-_x_)/sigma1_inodes)/ZCDF(1/sigma1_inodes))), 
																		0,1e25,
																		2*P_inodes*sigma_inodes*(1-Exp(-_x_*_x_/(2*sigma_inodes*sigma_inodes)))/sqrt2pi__+
																		(1-P_inodes)((1-P1_inodes)(mu_inodes*(1-ZCDF(-mu_inodes/sigma2_inodes)-ZCDF((mu_inodes-_x_)/sigma2_inodes))+
																		sigma2_inodes*(Exp(-mu_inodes*mu_inodes/(2*sigma2_inodes*sigma2_inodes))-Exp((_x_-mu_inodes)(mu_inodes-_x_)/(2*sigma2_inodes*sigma2_inodes)))/sqrt2pi__)/ZCDF(mu_inodes/sigma2_inodes)+
																		P1_inodes(1-ZCDF(-1/sigma1_inodes)-ZCDF((1-_x_)/sigma1_inodes)+
																		sigma1_inodes*(Exp(-1/(2*sigma1_inodes*sigma1_inodes))-Exp((_x_-1)(1-_x_)/(2*sigma1_inodes*sigma1_inodes)))/sqrt2pi__)/ZCDF(mu_inodes/sigma1_inodes))
																		);
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
	}		
	return 0;
}

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

function PopulateModelMatrix (ModelMatrixName&, EFV, modelType, catVar&,kappa&,R1&,R2&,R3&,R4&,R5&)
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
				  		ModelMatrixName[h-hshift][v-vshift] := catVar*t*EFV__[transition__];
			  			ModelMatrixName[v-vshift][h-hshift] := catVar*t*EFV__[transition2__];
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
					  		ModelMatrixName[h-hshift][v-vshift] := catVar*t*EFV__[transition__][nucPosInCodon__];
				  			ModelMatrixName[v-vshift][h-hshift] := catVar*t*EFV__[transition2__][nucPosInCodon__];
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
					  				ModelMatrixName[h-hshift][v-vshift] := kappa*catVar*t;
					  				ModelMatrixName[v-vshift][h-hshift] := kappa*catVar*t;
					  			}
					  			else
					  			{
					  				ModelMatrixName[h-hshift][v-vshift] := catVar*t;
					  				ModelMatrixName[v-vshift][h-hshift] := catVar*t;
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
						  			ModelMatrixName[h-hshift][v-vshift] := synRate*EFV__[transition__][nucPosInCodon__];
						  			ModelMatrixName[v-vshift][h-hshift] := synRate*EFV__[transition2__][nucPosInCodon__];
							  	}
						  		else
						  		{
							  		ModelMatrixName[h-hshift][v-vshift] := catVar*synRate*EFV__[transition__][nucPosInCodon__];
						  			ModelMatrixName[v-vshift][h-hshift] := catVar*synRate*EFV__[transition2__][nucPosInCodon__];
					  			}
					  		}
					  		else
					  		{
						  		if (rateKind == 2)
						  		{
							  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
							  		{
							  			ModelMatrixName[h-hshift][v-vshift] := R1*synRate*EFV__[transition__][nucPosInCodon__];
							  			ModelMatrixName[v-vshift][h-hshift] := R1*synRate*EFV__[transition2__][nucPosInCodon__];
								  	}
							  		else
							  		{
								  		ModelMatrixName[h-hshift][v-vshift] := R1*catVar*synRate*EFV__[transition__][nucPosInCodon__];
							  			ModelMatrixName[v-vshift][h-hshift] := R1*catVar*synRate*EFV__[transition2__][nucPosInCodon__];
						  			}
						  		}
						  		else
						  		{
							  		if (rateKind == 3)
							  		{
								  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
								  		{
								  			ModelMatrixName[h-hshift][v-vshift] := R2*synRate*EFV__[transition__][nucPosInCodon__];
								  			ModelMatrixName[v-vshift][h-hshift] := R2*synRate*EFV__[transition2__][nucPosInCodon__];
									  	}
								  		else
								  		{
									  		ModelMatrixName[h-hshift][v-vshift] := R2*catVar*synRate*EFV__[transition__][nucPosInCodon__];
								  			ModelMatrixName[v-vshift][h-hshift] := R2*catVar*synRate*EFV__[transition2__][nucPosInCodon__];
							  			}
							  		}
							  		else
							  		{
								  		if (rateKind == 4)
								  		{
									  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
									  		{
									  			ModelMatrixName[h-hshift][v-vshift] := R3*synRate*EFV__[transition__][nucPosInCodon__];
									  			ModelMatrixName[v-vshift][h-hshift] := R3*synRate*EFV__[transition2__][nucPosInCodon__];
										  	}
									  		else
									  		{
										  		ModelMatrixName[h-hshift][v-vshift] := R3*catVar*synRate*EFV__[transition__][nucPosInCodon__];
									  			ModelMatrixName[v-vshift][h-hshift] := R3*catVar*synRate*EFV__[transition2__][nucPosInCodon__];
								  			}
								  		}
								  		else
								  		{
									  		if (rateKind == 5)
									  		{
										  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
										  		{
										  			ModelMatrixName[h-hshift][v-vshift] := R4*synRate*EFV__[transition__][nucPosInCodon__];
										  			ModelMatrixName[v-vshift][h-hshift] := R4*synRate*EFV__[transition2__][nucPosInCodon__];
											  	}
										  		else
										  		{
											  		ModelMatrixName[h-hshift][v-vshift] := R4*catVar*synRate*EFV__[transition__][nucPosInCodon__];
										  			ModelMatrixName[v-vshift][h-hshift] := R4*catVar*synRate*EFV__[transition2__][nucPosInCodon__];
									  			}
									  		}
									  		else
									  		{
										  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) 
										  		{
										  			ModelMatrixName[h-hshift][v-vshift] := R5*synRate*EFV__[transition__][nucPosInCodon__];
										  			ModelMatrixName[v-vshift][h-hshift] := R5*synRate*EFV__[transition2__][nucPosInCodon__];
											  	}
										  		else
										  		{
											  		ModelMatrixName[h-hshift][v-vshift] := R5*catVar*synRate*EFV__[transition__][nucPosInCodon__];
										  			ModelMatrixName[v-vshift][h-hshift] := R5*catVar_inodes*synRate*EFV__[transition2__][nucPosInCodon__];
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


function GetDistributionParameters (distrInfo, infix)
{
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

	fprintf  (summaryFile,"\n\n---------------",infix,"-------------------\n\nExpected Value = ",E, " (sample variance = ",sampleVar,")\n");
	fprintf  (stdout,"\n\n---------------",infix,"-------------------\n\nExpected Value = ",E, " (sample variance = ",sampleVar,")\n");
	for (k=0; k<D; k=k+1)
	{
		fprintf (stdout,"\nRate[",Format(k+1,0,0),"]=",
				 Format(distrInfo[0][k],12,8), " (weight=", Format(distrInfo[1][k],9,7),")");
		fprintf (summaryFile,"\nRate[",Format(k+1,0,0),"]=",
				 Format(distrInfo[0][k],12,8), " (weight=", Format(distrInfo[1][k],9,7),")");
	}
	
	for (k=0; k<D; k=k+1)
	{
		if (distrInfo[0][k]>1) break;
	}
	
	return (D-k);
}


/* ____________________________________________________________________________________________________________________*/

NICETY_LEVEL = 3;

/* ----- Define genetic code and the appropriate state space dimensions ----- */

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

dataFile   = LAST_FILE_PATH;

DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);

dummyVar = FrameText ("-","|",2,"READ THE FOLLOWING DATA");

fprintf (stdout,"\n",ds);

/* ----- Read the tree and assign different models to leaves and internal nodes ----- */

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

newTreeString = treeString[0];

dummy = 0;

for (sI=1; sI<Abs(treeString); sI = sI+1)
{
	tC = treeString[sI];
	if ((tC == ",")||(tC == ")"))
	{
		if (treeString[sI-1]!=")")
		{
			newTreeString = newTreeString + "{LeafModel}" + tC;
			continue;
		}
		else
		{
			dummy = 1;
		}
	}
	newTreeString = newTreeString+tC;
}

if (dummy == 0)
{
	fprintf (stdout, "\n********* ERROR ***********\nThe tree must have at least one internal branch for this analysis to make sense\n\n");
	return 0;
}

treeString = newTreeString;
newTreeString = "";

fprintf (stdout, "RUNNING BIVARIATE (LEAVES AND INTERNAL NODES) POSITIVE SELECTION ANALYSIS ON\n", dataFile, "\n\n");

/* ----- Output files ----- */

fprintf (stdout, "\nThis analysis will write out two output files, one with MLEs and other summary data,\nand the other with marginal likelihoods for each site suitable for further data processing.\n");

SetDialogPrompt ("Please select a file to write summary data to:");
fprintf (PROMPT_FOR_FILE,CLEAR_FILE,"RUNNING BIVARIATE (LEAVES AND INTERNAL NODES) POSITIVE SELECTION ANALYSIS ON\n", dataFile, "\n\n");
summaryFile = LAST_FILE_PATH;

SetDialogPrompt ("Please select a file to write marginal likelihood data to:");
fprintf (PROMPT_FOR_FILE,CLEAR_FILE);
marginalFile = LAST_FILE_PATH;

/* ----- User selection phase for rate distribution ----- */

ChoiceList (rateSelection,"Rate variation for leaves and internal nodes.",1,SKIP_NONE,
			"Independent","Rate distributions are chosen independently.",
			"Constrained","Rate distributions are the same. Rates are still sampled independently for leaves and internal nodes at each site.");
		    

ChoiceList (rate4Leaves,"Rate distribution for leaves",1,SKIP_NONE,
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
		    "3 Normal","3 Normal");

if (rate4Leaves < 0)
{
	return 0;
}

dummy = SetWDistributionLeaves (rate4Leaves);
fprintf (stdout, "Using \"",ModelNames[rate4Leaves], "\" distribution for leaf nodes\n");
fprintf (summaryFile, "Using \"",ModelNames[rate4Leaves], "\" distribution for leaf nodes\n");

if (resp>0)
{
	fprintf (stdout, \nUsing ", 	   Format (resp,0,0), " discretized categories.\n");
	fprintf (summaryFile, \nUsing ", Format (resp,0,0), " discretized categories.\n");
}

if (rateSelection == 0)
{
	ChoiceList (rate4inodes,"Rate distribution for internal nodes",1,SKIP_NONE,
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
			    "3 Normal","3 Normal");
			    
	if (rate4inodes < 0)
	{
		return 0;
	}

}
else
{
	rate4inodes = rate4Leaves'
}

dummy = SetWDistributionINodes (rate4inodes);
fprintf (stdout, "Using \"",ModelNames[rate4inodes], "\" distribution for internal nodes\n");
fprintf (summaryFile, "Using \"",ModelNames[rate4inodes], "\" distribution for internal nodes\n");

if (resp>0)
{
	fprintf (stdout, "\nUsing ", 	   Format (resp,0,0), " discretized categories.\n");
	fprintf (summaryFile, "\nUsing ", Format (resp,0,0), " discretized categories.\n");
}

if (rateSelection == 1)
{
	dummy = ConstrainRateDistributionLeaves (rate4Leaves);
}

ChoiceList (modelSelection,"Select models for leaves and internal nodes.",1,SKIP_NONE,
			"Independent","Model choices for are made independently. Equilibrium frequencies are the same.",
			"Constrained","Rate matrices AND model parameters are the same.");


if (modelSelection < 0)
{
	return 0;
}

/* ----- User selection phase for the models ----- */

ChoiceList (modelTypeL,"Choose a model for leaves",1,SKIP_NONE,
			"MG94 1x4","Muse-Gaut 94 model with 4(-1) nucleotide frequency parameters (intra-codon position independent).",
			"MG94 3x4","Muse-Gaut 94 model with 12(-3) nucleotide frequency parameters (intra-codon position specific).",
			"GY94 1x4","Goldman-Yang 94 model with 4(-1) nucleotide frequency parameters (intra-codon position independent).",
			"GY94 3x4","Goldman-Yang 94 model with 12(-3) nucleotide frequency parameters (intra-codon position specific).",
			"MG94 Custom 1x4", "Muse-Gaut 94 crossed with an arbitrary nucelotide substitution model with 4(-1) nucleotide frequency parameters (intra-codon position independent).",
			"MG94 Custom 3x4", "Muse-Gaut 94 crossed with an arbitrary nucelotide substitution model with 12(-3) nucleotide frequency parameters (intra-codon position specific)."
);
		    
if (modelTypeL<0)
{
	return 0;
}

if ((modelTypeL==0)||(modelTypeL==2)||(modelTypeL==4))
{
	HarvestFrequencies (observedFreq,filteredData,1,1,0);
	vectorOfFrequencies = BuildCodonFrequencies4 (observedFreq);
}
else
{
	HarvestFrequencies (observedFreq,filteredData,3,1,1);
	vectorOfFrequencies = BuildCodonFrequencies12 (observedFreq);
}

if (modelTypeL>3)
{
	dummy = defineCustomModel (0);
}

LeafModelM = 0;
INodeModelM = 0;

dummy = reportModelType (modelTypeL," leaf nodes ", modelDesc);

if (modelSelection==1)
{
	global R1 = 1;
	global R2 = 1;
	global R3 = 1;
	global R4 = 1;
	global R5 = 1;
	global kappa = 1;
	dummy = PopulateModelMatrix ("LeafModelM",observedFreq,modelTypeL,"category_leaves", "kappa","R1","R2","R3","R4","R5");
	Model 	LeafModel = (LeafModelM,vectorOfFrequencies,dummy);
	dummy = PopulateModelMatrix ("INodeModelM",observedFreq,modelTypeL,"category_inodes", "kappa","R1","R2","R3","R4","R5");
	Model 	INodeModel = (INodeModelM,vectorOfFrequencies,dummy);
	dummy = reportModelType (modelTypeL," internal nodes ", modelDesc);
}
else
{
	global R1_leaves 	= 1;
	global R2_leaves 	= 1;
	global R3_leaves 	= 1;
	global R4_leaves 	= 1;
	global R5_leaves 	= 1;
	global kappa_leaves = 1;

	global R1_inodes 	= 1;
	global R2_inodes 	= 1;
	global R3_inodes 	= 1;
	global R4_inodes 	= 1;
	global R5_inodes 	= 1;
	global kappa_inodes = 1;

	dummy = PopulateModelMatrix ("LeafModelM",observedFreq,modelTypeL,"category_leaves", "kappa_leaves","R1_leaves","R2_leaves","R3_leaves","R4_leaves","R5_leaves");
	Model 	LeafModel = (LeafModelM,vectorOfFrequencies,dummy);

	ChoiceList (modelTypeI,"Choose a model for internal nodes",1,SKIP_NONE,
				"MG94","Muse-Gaut 94 model.",
				"GY94","Goldman-Yang 94 model",
				"MG94 Custom", "Muse-Gaut 94 crossed with an arbitrary nucelotide substitution model.");

	if (modelTypeI<0)
	{
		return 0;
	}
	modelTypeI = modelTypeI * 2;

	if (modelTypeI>3)
	{
		dummy = defineCustomModel (0);
	}

	dummy = PopulateModelMatrix ("INodeModelM",observedFreq,modelTypeI,"category_inodes", "kappa_inodes","R1_inodes","R2_inodes","R3_inodes","R4_inodes","R5_inodes");
	Model 	INodeModel = (INodeModelM,vectorOfFrequencies,dummy);
	dummy = reportModelType (modelTypeI," internal nodes ", modelDesc);
}


fprintf (stdout, "\nChoose the cutoff (0 to 1) for posterior of dN/dS>1 for a site to be considered under selective pressure:");
fscanf  (stdin, "Number",psigLevel);
if ((psigLevel <= 0)||(psigLevel>1))
{
	psigLevel = .95;
}
fprintf (stdout, "\n>Using ", psigLevel , " cutoff\n");
fprintf (summaryFile, "\nUsing ", psigLevel , " cutoff for positive selection.\n");
				 
Tree  givenTree 	  = treeString;
LikelihoodFunction lf = (filteredData,givenTree);
timer = Time(1);
Optimize (res,lf);

dummyVar = FrameText ("-","|",2,"RESULTS");
fprintf (stdout     , lf, "\n");
fprintf (summaryFile,"\n\n*** RESULTS ***\n\n", lf, "\n");

GetInformation (leavesInfo,category_leaves);
GetInformation (inodesInfo,category_inodes);

lCategs    = Columns (leavesInfo);
iCategs    = Columns (inodesInfo);
mCategs	   = lCategs * iCategs;

ConstructCategoryMatrix(marginals,lf,COMPLETE);

GetInformation (categVarIDs,lf);

if (categVarIDs[0]!="category_leaves")
{
	marginalsCorrected = marginals;
	fprintf (MESSAGE_LOG,"\n.Adjusting marginal matrix rows.\n"); 
	for (h=0; h<Columns(marginals); h=h+1)
	{
		transition = 0;
		for (diff=0; diff<iCategs; diff = diff+1)
		{
			for (v=diff; v<mCategs; v=v+iCategs)
			{
				marginalsCorrected[transition][h] = marginals[v][h];
				transition = transition+1;
			}
		}
	}
	marginals = marginalsCorrected;
	marginalsCorrected = 0;
}

fprintf (marginalFile, leavesInfo, inodesInfo, marginals);

leafOver1  = GetDistributionParameters (leavesInfo, " LEAF NODES RATE DISTRIBUTION ");
inodeOver1 = GetDistributionParameters (inodesInfo, " INTERNAL NODES RATE DISTRIBUTION ");

if ((inodeOver1<1)&&(leafOver1<1))
{
	fprintf (stdout, "\nRATES>1 ARE PRESENT NEITHER IN LEAVES NOR IN INTERNAL NODES.\n");
	fprintf (summaryFile, "\nRATES>1 ARE PRESENT NEITHER IN LEAVES NOR IN INTERNAL NODES.\n");
	return 0;
}

colm = Columns(marginals);


resultMatrix = {colm, 3};
weightMatrix = {1, mCategs};

for (c1 = 0; c1 < lCategs; c1 = c1 + 1)
{
	for (c2 = 0; c2 < iCategs; c2 = c2 + 1)
	{
		weightMatrix [c1 * iCategs + c2] = leavesInfo[1][c1] * inodesInfo[1][c2];
	}
}

fprintf (stdout, weightMatrix, "\n\n");

for (c1 = 0; c1 < colm; c1 = c1+1)
{
	E = 0;
	for (c2 = 0; c2 < mCategs; c2 = c2 + 1)
	{
		E = E + weightMatrix[c2] * marginals[c2][c1];
	}
	resultMatrix [c1][2] = E;
}

if (leafOver1>0)
{
	for (c1 = 0; c1 < colm; c1 = c1+1)
	{
		E = 0;
		for (c2 = (lCategs-leafOver1)*iCategs; c2 < mCategs; c2 = c2+1)
		{
			E = E + weightMatrix[c2] * marginals[c2][c1];
		}
		resultMatrix[c1][0] = E/resultMatrix[c1][2];
	}
}
else
{
	fprintf (stdout, "\n\nRATES>1 ARE NOT PRESENT IN LEAVES.\n\n");
	fprintf (summaryFile, "\n\nRATES>1 ARE NOT PRESENT IN LEAVES.\n\n");
}

if (inodeOver1>0)
{
	for (c1 = 0; c1 < colm; c1 = c1+1)
	{
		E = 0;
		for (c2 = (iCategs-inodeOver1); c2 < mCategs; c2 = c2+iCategs)
		{
			E = E + weightMatrix[c2] * marginals[c2][c1];
		}
		resultMatrix[c1][1] = E/resultMatrix[c1][2];
	}
}
else
{
	fprintf (stdout, "\n\nRATES>1 ARE NOT PRESENT IN INTERNAL NODES.\n");
	fprintf (summaryFile, "\n\nRATES>1 ARE NOT PRESENT IN INTERNAL NODES.\n");
}


dummyVar = FrameText ("-","|",2,"POSITIVE SELECTION RESULT TABLE");
fprintf (summaryFile, "\n\nPOSITIVE SELECTION RESULT TABLE\n\n");

fprintf (stdout, "Site  P{Leaf dN/dS >= 1} P{Internal Node dN/dS >= 1} PS (Level ", psigLevel, ")\n",
				 "                                                     LI\n");
fprintf (summaryFile, "Site  P{Leaf dN/dS >= 1} P{Internal Node dN/dS >= 1} PS (Level ", psigLevel, ")\n",
				 "                                                     LI\n");
				 
if (inodeOver1>0)
{
	if (leafOver1>0)
	{
		for (c1 = 0; c1 < colm; c1 = c1+1)
		{
			if (resultMatrix[c1][0]>=psigLevel)
			{
				pss = "+";
			}
			else
			{
				pss = "-";
			}
			if (resultMatrix[c1][1]>=psigLevel)
			{
				pss = pss+"+";
			}
			else
			{
				pss = pss+"-";
			}
			
			fprintf (stdout, Format (c1+1,4,0), " ", Format (resultMatrix[c1][0], 18, 5), "  ", Format (resultMatrix[c1][0], 28, 5), " ", pss, "\n");
			fprintf (summaryFile, Format (c1+1,4,0), " ", Format (resultMatrix[c1][0], 18, 5), "  ", Format (resultMatrix[c1][0], 28, 5), " ", pss, "\n");
		}
	}
	else
	{
		for (c1 = 0; c1 < colm; c1 = c1+1)
		{
			pss = "*";
			
			if (resultMatrix[c1][1]>=psigLevel)
			{
				pss = pss+"+";
			}
			else
			{
				pss = pss+"-";
			}
			
			fprintf (stdout, Format (c1+1,4,0), " ", Format (0, 18, 5), " ", Format (resultMatrix[c1][0], 28, 5), " ", pss, "\n");
			fprintf (summaryFile, Format (c1+1,4,0), " ", Format (0, 18, 5), " ", Format (resultMatrix[c1][0], 28, 5), " ", pss, "\n");
		}
	
	}
}
else
{
	for (c1 = 0; c1 < colm; c1 = c1+1)
	{
		if (resultMatrix[c1][0]>=psigLevel)
		{
			pss = "+";
		}
		else
		{
			pss = "-";
		}
		pss = pss+"*";
		
		fprintf (stdout, Format (c1+1,4,0), " ", Format (resultMatrix[c1][0], 18, 5), " ", Format (0, 28, 5), " ", pss, "\n");
		fprintf (summaryFile, Format (c1+1,4,0), " ", Format (resultMatrix[c1][0], 18, 5), " ", Format (0, 28, 5), " ", pss, "\n");
	}
}


resultMatrix2 = {colm, 2};

for (c2=0; c2<colm; c2=c2+1)
{	
	resultMatrix2[c2][0] = resultMatrix[c2][0];
	resultMatrix2[c2][1] = resultMatrix[c2][1];
}

labelMatrix = {1,2};
labelMatrix[0] = "Leaves";
labelMatrix[1] = "Internal Nodes";

OpenWindow (CHARTWINDOW,{{"Prob dN/dS>1"}
						   {"labelMatrix"},
						   {"resultMatrix2"},
						   {"Contrast Bars"},
						   {"Index"},
						   {labelMatrix[0]+";"+labelMatrix[1]},
						   {"Site"},
						   {""},
						   {"Prob dN/dS>1"},
						   {"0"}},
						   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");

