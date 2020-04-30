fprintf (stdout, "\nRunning an HKY+gamma model fit on a simulated alignment with 8 sequences and 1320 nucleotides with 4 blocks of contiguous rates (HMM test case)\n");

OPTIMIZE_SUMMATION_ORDER = 1;

timer = Time (1);

global SmallCodon_part_Shared_TVTS=1;
global SmallCodon_part_Shape=0.5;
SmallCodon_part_Shape:>0.01;
SmallCodon_part_Shape:<100;
SmallCodon_part_Categ.weights={
{              0.25}
{              0.25}
{              0.25}
{              0.25}
}
;

global						   lambda = 0.25;
lambda :< 1/3;

HMM_transition_matrix 		   = {{1-3*lambda,lambda,lambda,lambda}
								  {lambda,1-3*lambda,lambda,lambda}
								  {lambda,lambda,1-3*lambda,lambda}
								  {lambda,lambda,lambda,1-3*lambda}};
								  
HMM_starting_frequencies	   = {{0.25,0.25,0.25,0.25}};	
Model	HMM_model			   = (HMM_transition_matrix,HMM_starting_frequencies,0);
 

category SmallCodon_part_Categ =(4,SmallCodon_part_Categ.weights,MEAN,
								  GammaDist(_x_,SmallCodon_part_Shape,SmallCodon_part_Shape),
								  CGammaDist(_x_,SmallCodon_part_Shape,SmallCodon_part_Shape),0,1e+25,
								  CGammaDist(_x_,SmallCodon_part_Shape+1,SmallCodon_part_Shape),
								  HMM_model);
								  
SmallCodon_part_HKY85={4,4};
SmallCodon_part_HKY85[0][1]:=t*SmallCodon_part_Shared_TVTS*SmallCodon_part_Categ;
SmallCodon_part_HKY85[0][2]:=t*SmallCodon_part_Categ;
SmallCodon_part_HKY85[0][3]:=t*SmallCodon_part_Shared_TVTS*SmallCodon_part_Categ;
SmallCodon_part_HKY85[1][0]:=t*SmallCodon_part_Shared_TVTS*SmallCodon_part_Categ;
SmallCodon_part_HKY85[1][2]:=t*SmallCodon_part_Shared_TVTS*SmallCodon_part_Categ;
SmallCodon_part_HKY85[1][3]:=t*SmallCodon_part_Categ;
SmallCodon_part_HKY85[2][0]:=t*SmallCodon_part_Categ;
SmallCodon_part_HKY85[2][1]:=t*SmallCodon_part_Shared_TVTS*SmallCodon_part_Categ;
SmallCodon_part_HKY85[2][3]:=t*SmallCodon_part_Shared_TVTS*SmallCodon_part_Categ;
SmallCodon_part_HKY85[3][0]:=t*SmallCodon_part_Shared_TVTS*SmallCodon_part_Categ;
SmallCodon_part_HKY85[3][1]:=t*SmallCodon_part_Categ;
SmallCodon_part_HKY85[3][2]:=t*SmallCodon_part_Shared_TVTS*SmallCodon_part_Categ;


DataSet 	  SmallCodon 	  = ReadDataFile	(PATH_TO_CURRENT_BF + "../data/HMM4_synthetic.fas");
DataSetFilter SmallCodon_part = CreateFilter	(SmallCodon,1,"","4,5,7,6,1,0,2,3");

HarvestFrequencies 								(SmallCodon_part_Freqs,SmallCodon_part,1,1,1);
Model SmallCodon_part_HKY85_model	=			(SmallCodon_part_HKY85,SmallCodon_part_Freqs);
Tree SmallCodon_tree				=			((((D_CD_83_ELI_ACC_K03454,D_CD_83_NDK_ACC_M27323)Node3,D_UG_94_94UG114_ACC_U88824)Node2,D_CD_84_84ZR085_ACC_U88822)Node1,B_US_83_RF_ACC_M17451,((B_FR_83_HXB2_ACC_K03455,B_US_86_JRFL_ACC_U63632)Node10,B_US_90_WEAU160_ACC_U21135)Node9);
LikelihoodFunction SmallCodon_LF 	= 			(SmallCodon_part,SmallCodon_tree);

Optimize(res_SmallCodon_LF,SmallCodon_LF);

timer2 								= Time (1);
expectedLL 							= -3243.9885;
diffLL	   							= Abs(expectedLL - res_SmallCodon_LF[1][0]);
fprintf 							(stdout, SmallCodon_LF, "\nTest optimization took ", timer2-timer, " seconds.\n", diffLL , " difference between obtained and expected likelihood\n\n");

ConstructCategoryMatrix				(mx,  SmallCodon_LF, SHORT);
ConstructCategoryMatrix				(mx2,  SmallCodon_LF);

fprintf (stdout, mx2, "\n");
/* this stores the Viterbi path in mx */

expected = {"0":{{200,3,1}},"1":{{681,1,0}},"2":{{823,0,2}}};

switches = 0;

for (k    = 1; k < Columns (mx); k = k+1)
{
	if (mx[k] != mx[k-1])
	{
		fprintf (stdout, "Switch from ", mx[k-1], " to ", mx[k], " at site ", k, "\n");
		expectedS = expected[switches];
		if (expectedS[0] != k || expectedS[1] != mx[k-1] || expectedS[2] != mx[k])
		{
			fprintf (stdout, "Test FAILED: expected  a switch from ", expectedS[1], " to ", expectedS[2], " at site ", expectedS[0], "\n");
			return 1;
		}
		switches = switches + 1;
	}
}

return 0;
