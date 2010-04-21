fprintf (stdout, "\nRunning an HKY+gamma model fit on a simulated alignment with 8 sequences and 1320 nucleotides with 2 blocks from different trees (equal size)\n");

OPTIMIZE_SUMMATION_ORDER = 1;

timer = Time (1);

global						   lambda = 0.1;
lambda :< 1;
lambda :> 0.000;

HMM_transition_matrix 		   		= {{1-lambda,lambda}
   								  	  {lambda,1-lambda}};
								  
HMM_starting_frequencies	   		= {{0.5,0.5}};	
Model	HMM_model			   		= (HMM_transition_matrix,HMM_starting_frequencies,0);

category HMM_variable			    =(2,{{0.5,0.5}},,{{0,1}},
									  ,0,1e+25,
									  ,
									  HMM_model);
 

DataSet 	  SmallCodon 	  		= ReadDataFile	("../data/HMM2_synthetic.fas");
DataSetFilter SmallCodon_part 		= CreateFilter	(SmallCodon,1);

global		  SmallCodon_part_Shared_TVTS = 1;

SmallCodon_part_HKY85={4,4};
SmallCodon_part_HKY85[0][1]:=t*SmallCodon_part_Shared_TVTS;
SmallCodon_part_HKY85[0][2]:=t;
SmallCodon_part_HKY85[0][3]:=t*SmallCodon_part_Shared_TVTS;
SmallCodon_part_HKY85[1][0]:=t*SmallCodon_part_Shared_TVTS;
SmallCodon_part_HKY85[1][2]:=t*SmallCodon_part_Shared_TVTS;
SmallCodon_part_HKY85[1][3]:=t;
SmallCodon_part_HKY85[2][0]:=t;
SmallCodon_part_HKY85[2][1]:=t*SmallCodon_part_Shared_TVTS;
SmallCodon_part_HKY85[2][3]:=t*SmallCodon_part_Shared_TVTS;
SmallCodon_part_HKY85[3][0]:=t*SmallCodon_part_Shared_TVTS;
SmallCodon_part_HKY85[3][1]:=t;
SmallCodon_part_HKY85[3][2]:=t*SmallCodon_part_Shared_TVTS;


HarvestFrequencies 								(SmallCodon_part_Freqs,SmallCodon_part,1,1,1);
Model SmallCodon_part_HKY85_model	=			(SmallCodon_part_HKY85,SmallCodon_part_Freqs);
Tree tree1							=			((((D_CD_83_ELI_ACC_K03454,D_CD_83_NDK_ACC_M27323)Node3,D_UG_94_94UG114_ACC_U88824)Node2,D_CD_84_84ZR085_ACC_U88822)Node1,B_US_83_RF_ACC_M17451,((B_FR_83_HXB2_ACC_K03455,B_US_86_JRFL_ACC_U63632)Node10,B_US_90_WEAU160_ACC_U21135)Node9);
Tree tree2							=			((((D_CD_83_ELI_ACC_K03454:0.0252441,D_CD_83_NDK_ACC_M27323:0.0135632):0.00881346,D_UG_94_94UG114_ACC_U88824:0.0619783):0.000782367,D_CD_84_84ZR085_ACC_U88822:0.0207434):0.0351827,((B_FR_83_HXB2_ACC_K03455:0.0195062,B_US_86_JRFL_ACC_U63632:0.0157251):0.00159642,B_US_83_RF_ACC_M17451:0.0284274):0,B_US_90_WEAU160_ACC_U21135:0.0245537);
LikelihoodFunction SmallCodon_LF 	= 			(SmallCodon_part,tree1,SmallCodon_part,tree2,"HMM_variable");

Optimize										(res_SmallCodon_LF,SmallCodon_LF);

timer2 								= Time (1);
expectedLL 							= -3580.54867;
diffLL	   							= Abs(expectedLL - res_SmallCodon_LF[1][0]);
fprintf 							(stdout, SmallCodon_LF, "\nTest optimization took ", timer2-timer, " seconds.\n", diffLL , " difference between obtained and expected likelihood\n\n");

ConstructCategoryMatrix				(mx, SmallCodon_LF, SHORT);

expected = {"0":{{731,1,0}}};

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
