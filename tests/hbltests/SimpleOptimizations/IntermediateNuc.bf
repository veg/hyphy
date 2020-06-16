/* test preamble */

	_testDescription 		= " fit the HKY85 model to an Influenza A alignment with 349 sequences and 967 nucleotides";
	_expectedLL 			= -11389.4543728884;
	ExecuteAFile 			("../Shared/TestInstrumentation.bf");
	startTestTimer 			(_testDescription);

/* end test preamble */

global flu_part2_Shared_TVTS=0.25;

flu_part2_HKY85={4,4};
flu_part2_HKY85[0][1]:=t*flu_part2_Shared_TVTS;
flu_part2_HKY85[0][2]:=t;
flu_part2_HKY85[0][3]:=t*flu_part2_Shared_TVTS;
flu_part2_HKY85[1][0]:=t*flu_part2_Shared_TVTS;
flu_part2_HKY85[1][2]:=t*flu_part2_Shared_TVTS;
flu_part2_HKY85[1][3]:=t;
flu_part2_HKY85[2][0]:=t;
flu_part2_HKY85[2][1]:=t*flu_part2_Shared_TVTS;
flu_part2_HKY85[2][3]:=t*flu_part2_Shared_TVTS;
flu_part2_HKY85[3][0]:=t*flu_part2_Shared_TVTS;
flu_part2_HKY85[3][1]:=t;
flu_part2_HKY85[3][2]:=t*flu_part2_Shared_TVTS;

flu_part2_Freqs={
{    0.321488786102}
{    0.221264478507}
{    0.225765445963}
{    0.231481289428}
}
;

Model flu_part2_HKY85_model=(flu_part2_HKY85,flu_part2_Freqs);

UseModel (flu_part2_HKY85_model);
DataSet flu 				= 	ReadDataFile(PATH_TO_CURRENT_BF + "/../data/fluHA.nex");
DataSetFilter flu_part2 	= 	CreateFilter(flu,1,"0-966","101,26,75,103,100,102,49,89,29,9,4,0,19,91,81,57,43,25,54,66,41,30,51,48,52,83,36,73,34,28,45,39,20,6,76,72,69,56,46,62,38,31,78,15,12,105,104,23,97,86,1,77,88,85,60,119,33,107,106,122,121,120,21,93,87,79,27,84,82,63,61,11,67,24,59,32,53,47,42,99,98,96,70,95,117,118,115,114,8,116,113,80,94,92,108,112,111,68,65,110,74,71,37,109,22,90,35,123,124,159,158,157,156,155,153,150,148,151,149,143,58,50,44,142,40,55,154,147,152,146,145,144,141,140,139,138,137,136,135,134,132,133,13,5,131,10,130,127,125,64,129,128,126,160,182,181,16,163,162,2,161,14,180,179,178,177,176,174,175,172,173,170,169,168,166,165,164,171,167,213,212,211,208,195,209,222,217,216,215,207,220,219,218,221,214,210,203,202,201,200,199,198,197,196,206,205,204,191,190,187,185,186,184,183,194,193,192,189,188,241,240,238,239,237,236,235,17,224,223,234,18,233,232,231,230,229,228,227,226,225,250,248,249,247,3,246,244,245,243,242,251,256,255,254,253,252,286,285,284,283,281,280,282,279,275,274,278,276,277,273,272,267,266,264,265,263,262,257,260,259,261,258,271,270,269,268,287,296,295,293,292,294,289,288,291,290,298,297,299,313,312,311,310,309,308,303,302,305,304,301,300,307,306,315,314,327,324,322,321,323,320,319,326,325,328,318,317,316,332,330,331,329,336,335,334,7,333,346,345,344,343,347,348,342,340,341,339,338,337");
Tree Tree_12				=	DATAFILE_TREE;
Tree_12.AG207_796.t        := Tree_12.IN1_1196.t;

LikelihoodFunction flu_LF2  = 	(flu_part2,Tree_12);


VERBOSITY_LEVEL 			  		= 	1;
OPTIMIZATION_METHOD					=   4;
OPTIMIZATION_PROGRESS_QUANTUM 		= 	0.5;
OPTIMIZATION_PROGRESS_STATUS  		= 	"OPTIMIZING THE LIKELIHOOD FUNCTION";
OPTIMIZATION_PROGRESS_TEMPLATE 		= 	"$1 $2 $3% $4 $5 $6";

Optimize								(res_flu_LF2,flu_LF2);

fprintf (stdout, res_flu_LF2[1][0], "\n");

/* test epilogue */
	timeMatrix = endTestTimer 				  (_testDescription);
	if (logTestResult (Abs (res_flu_LF2[1][0] - _expectedLL) < 2*OPTIMIZATION_PRECISION))
	{
		return timeMatrix;
	}
	return 0;
/* end test epilogue */
