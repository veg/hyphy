skipCodeSelectionStep = 1;
ExecuteAFile			("../TemplateModels/chooseGeneticCode.def");
ApplyGeneticCodeTable  (0);

_HXB_env_offset					= 6224;


_HXB_aa_offsets					= 
							    { "pre"		:0,
								  "v1"		:130,
								  "c1"		:149,
								  "v2"		:157,
								  "c2"		:196,
								  "v3"		:295,
								  "c3"		:331,
								  "v4"		:384,
								  "c4"		:418,
								  "v5"		:459,
								  "post"    :471,
								  "gp41ecto":511,
								  "gp41endo":684
								};

_HXB_aa_offset_matrix			= {{0,130,149,157,196,295,331,384,418,459,471,511,684,856}};
								  
_HXB_Annotation = {};
_HXB_Annotation ["pre"] 		= {{_HXB_env_offset__,_HXB_env_offset__+389}};
_HXB_Annotation ["v1"]  		= {{_HXB_env_offset__+390,_HXB_env_offset__+446}};
_HXB_Annotation ["c1"]  		= {{_HXB_env_offset__+447,_HXB_env_offset__+470}};
_HXB_Annotation ["v2"]  		= {{_HXB_env_offset__+471,_HXB_env_offset__+587}};
_HXB_Annotation ["c2"]  		= {{_HXB_env_offset__+588,_HXB_env_offset__+884}};
_HXB_Annotation ["v3"]  		= {{_HXB_env_offset__+885,_HXB_env_offset__+992}};
_HXB_Annotation ["c3"]  		= {{_HXB_env_offset__+993,_HXB_env_offset__+1151}};
_HXB_Annotation ["v4"]  		= {{_HXB_env_offset__+1152,_HXB_env_offset__+1253}};
_HXB_Annotation ["c4"]  		= {{_HXB_env_offset__+1254,_HXB_env_offset__+1376}};
_HXB_Annotation ["v5"]  		= {{_HXB_env_offset__+1377,_HXB_env_offset__+1412}};
_HXB_Annotation ["post"]		= {{_HXB_env_offset__+1413,_HXB_env_offset__+1532}};
_HXB_Annotation ["gp41ecto"]	= {{_HXB_env_offset__+1533,_HXB_env_offset__+2051}};
_HXB_Annotation ["gp41endo"]	= {{_HXB_env_offset__+2052,_HXB_env_offset__+2567}};
_HXB_Annotation ["gp120"]       = {{_HXB_env_offset__,_HXB_env_offset__+1532}};
_HXB_Annotation ["gp160"]       = {{_HXB_env_offset__,_HXB_env_offset__+2567}};
_HXB_Annotation ["gp41"]       = {{_HXB_env_offset__+1533,_HXB_env_offset__+2567}};


_HXB2_Sequence_					= "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGGGCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGCTTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAGAGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCGCTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGGGGAGTGGCGAGCCCTCAGATCCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGGCGCCCGAACAGGGACCTGAAAGCGAAAGGGAAACCAGAGGAGCTCTCTCGACGCAGGACTCGGCTTGCTGAAGCGCGCACGGCAAGAGGCGAGGGGCGGCGACTGGTGAGTACGCCAAAAATTTTGACTAGCGGAGGCTAGAAGGAGAGAGATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATATAGTATGGGCAAGCAGGGAGCTAGAACGATTCGCAGTTAATCCTGGCCTGTTAGAAACATCAGAAGGCTGTAGACAAATACTGGGACAGCTACAACCATCCCTTCAGACAGGATCAGAAGAACTTAGATCATTATATAATACAGTAGCAACCCTCTATTGTGTGCATCAAAGGATAGAGATAAAAGACACCAAGGAAGCTTTAGACAAGATAGAGGAAGAGCAAAACAAAAGTAAGAAAAAAGCACAGCAAGCAGCAGCTGACACAGGACACAGCAATCAGGTCAGCCAAAATTACCCTATAGTGCAGAACATCCAGGGGCAAATGGTACATCAGGCCATATCACCTAGAACTTTAAATGCATGGGTAAAAGTAGTAGAAGAGAAGGCTTTCAGCCCAGAAGTGATACCCATGTTTTCAGCATTATCAGAAGGAGCCACCCCACAAGATTTAAACACCATGCTAAACACAGTGGGGGGACATCAAGCAGCCATGCAAATGTTAAAAGAGACCATCAATGAGGAAGCTGCAGAATGGGATAGAGTGCATCCAGTGCATGCAGGGCCTATTGCACCAGGCCAGATGAGAGAACCAAGGGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACAAATAATCCACCTATCCCAGTAGGAGAAATTTATAAAAGATGGATAATCCTGGGATTAAATAAAATAGTAAGAATGTATAGCCCTACCAGCATTCTGGACATAAGACAAGGACCAAAGGAACCCTTTAGAGACTATGTAGACCGGTTCTATAAAACTCTAAGAGCCGAGCAAGCTTCACAGGAGGTAAAAAATTGGATGACAGAAACCTTGTTGGTCCAAAATGCGAACCCAGATTGTAAGACTATTTTAAAAGCATTGGGACCAGCGGCTACACTAGAAGAAATGATGACAGCATGTCAGGGAGTAGGAGGACCCGGCCATAAGGCAAGAGTTTTGGCTGAAGCAATGAGCCAAGTAACAAATTCAGCTACCATAATGATGCAGAGAGGCAATTTTAGGAACCAAAGAAAGATTGTTAAGTGTTTCAATTGTGGCAAAGAAGGGCACACAGCCAGAAATTGCAGGGCCCCTAGGAAAAAGGGCTGTTGGAAATGTGGAAAGGAAGGACACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAAGATCTGGCCTTCCTACAAGGGAAGGCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTTCCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAACTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAGGCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAACATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCCTGAGTGGGAGTTTGTTAATACCCCTCCCTTAGTGAAATTATGGTACCAGTTAGAGAAAGAACCCATAGTAGGAGCAGAAACCTTCTATGTAGATGGGGCAGCTAACAGGGAGACTAAATTAGGAAAAGCAGGATATGTTACTAATAGAGGAAGACAAAAAGTTGTCACCCTAACTGACACAACAAATCAGAAGACTGAGTTACAAGCAATTTATCTAGCTTTGCAGGATTCGGGATTAGAAGTAAACATAGTAACAGACTCACAATATGCATTAGGAATCATTCAAGCACAACCAGATCAAAGTGAATCAGAGTTAGTCAATCAAATAATAGAGCAGTTAATAAAAAAGGAAAAGGTCTATCTGGCATGGGTACCAGCACACAAAGGAATTGGAGGAAATGAACAAGTAGATAAATTAGTCAGTGCTGGAATCAGGAAAGTACTATTTTTAGATGGAATAGATAAGGCCCAAGATGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAGTGATTTTAACCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGAGAAGCCATGCATGGACAAGTAGACTGTAGTCCAGGAATATGGCAACTAGATTGTACACATTTAGAAGGAAAAGTTATCCTGGTAGCAGTTCATGTAGCCAGTGGATATATAGAAGCAGAAGTTATTCCAGCAGAAACAGGGCAGGAAACAGCATATTTTCTTTTAAAATTAGCAGGAAGATGGCCAGTAAAAACAATACATACTGACAATGGCAGCAATTTCACCGGTGCTACGGTTAGGGCCGCCTGTTGGTGGGCGGGAATCAAGCAGGAATTTGGAATTCCCTACAATCCCCAAAGTCAAGGAGTAGTAGAATCTATGAATAAAGAATTAAAGAAAATTATAGGACAGGTAAGAGATCAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAAATCCACTTTGGAAAGGACCAGCAAAGCTCCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGATAATAGTGACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATTAGGGATTATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGATTAGAACATGGAAAAGTTTAGTAAAACACCATATGTATGTTTCAGGGAAAGCTAGGGGATGGTTTTATAGACATCACTATGAAAGCCCTCATCCAAGAATAAGTTCAGAAGTACACATCCCACTAGGGGATGCTAGATTGGTAATAACAACATATTGGGGTCTGCATACAGGAGAAAGAGACTGGCATTTGGGTCAGGGAGTCTCCATAGAATGGAGGAAAAAGAGATATAGCACACAAGTAGACCCTGAACTAGCAGACCAACTAATTCATCTGTATTACTTTGACTGTTTTTCAGACTCTGCTATAAGAAAGGCCTTATTAGGACACATAGTTAGCCCTAGGTGTGAATATCAAGCAGGACATAACAAGGTAGGATCTCTACAATACTTGGCACTAGCAGCATTAATAACACCAAAAAAGATAAAGCCACCTTTGCCTAGTGTTACGAAACTGACAGAGGATAGATGGAACAAGCCCCAGAAGACCAAGGGCCACAGAGGGAGCCACACAATGAATGGACACTAGAGCTTTTAGAGGAGCTTAAGAATGAAGCTGTTAGACATTTTCCTAGGATTTGGCTCCATGGCTTAGGGCAACATATCTATGAAACTTATGGGGATACTTGGGCAGGAGTGGAAGCCATAATAAGAATTCTGCAACAACTGCTGTTTATCCATTTTCAGAATTGGGTGTCGACATAGCAGAATAGGCGTTACTCGACAGAGGAGAGCAAGAAATGGAGCCAGTAGATCCTAGACTAGAGCCCTGGAAGCATCCAGGAAGTCAGCCTAAAACTGCTTGTACCAATTGCTATTGTAAAAAGTGTTGCTTTCATTGCCAAGTTTGTTTCATAACAAAAGCCTTAGGCATCTCCTATGGCAGGAAGAAGCGGAGACAGCGACGAAGAGCTCATCAGAACAGTCAGACTCATCAAGCTTCTCTATCAAAGCAGTAAGTAGTACATGTAACGCAACCTATACCAATAGTAGCAATAGTAGCATTAGTAGTAGCAATAATAATAGCAATAGTTGTGTGGTCCATAGTAATCATAGAATATAGGAAAATATTAAGACAAAGAAAAATAGACAGGTTAATTGATAGACTAATAGAAAGAGCAGAAGACAGTGGCAATGAGAGTGAAGGAGAAATATCAGCACTTGTGGAGATGGGGGTGGAGATGGGGCACCATGCTCCTTGGGATGTTGATGATCTGTAGTGCTACAGAAAAATTGTGGGTCACAGTCTATTATGGGGTACCTGTGTGGAAGGAAGCAACCACCACTCTATTTTGTGCATCAGATGCTAAAGCATATGATACAGAGGTACATAATGTTTGGGCCACACATGCCTGTGTACCCACAGACCCCAACCCACAAGAAGTAGTATTGGTAAATGTGACAGAAAATTTTAACATGTGGAAAAATGACATGGTAGAACAGATGCATGAGGATATAATCAGTTTATGGGATCAAAGCCTAAAGCCATGTGTAAAATTAACCCCACTCTGTGTTAGTTTAAAGTGCACTGATTTGAAGAATGATACTAATACCAATAGTAGTAGCGGGAGAATGATAATGGAGAAAGGAGAGATAAAAAACTGCTCTTTCAATATCAGCACAAGCATAAGAGGTAAGGTGCAGAAAGAATATGCATTTTTTTATAAACTTGATATAATACCAATAGATAATGATACTACCAGCTATAAGTTGACAAGTTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAGCCAATTCCCATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAATGTAATAATAAGACGTTCAATGGAACAGGACCATGTACAAATGTCAGCACAGTACAATGTACACATGGAATTAGGCCAGTAGTATCAACTCAACTGCTGTTAAATGGCAGTCTAGCAGAAGAAGAGGTAGTAATTAGATCTGTCAATTTCACGGACAATGCTAAAACCATAATAGTACAGCTGAACACATCTGTAGAAATTAATTGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGAGAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGTAACATTAGTAGAGCAAAATGGAATAACACTTTAAAACAGATAGCTAGCAAATTAAGAGAACAATTTGGAAATAATAAAACAATAATCTTTAAGCAATCCTCAGGAGGGGACCCAGAAATTGTAACGCACAGTTTTAATTGTGGAGGGGAATTTTTCTACTGTAATTCAACACAACTGTTTAATAGTACTTGGTTTAATAGTACTTGGAGTACTGAAGGGTCAAATAACACTGAAGGAAGTGACACAATCACCCTCCCATGCAGAATAAAACAAATTATAAACATGTGGCAGAAAGTAGGAAAAGCAATGTATGCCCCTCCCATCAGTGGACAAATTAGATGTTCATCAAATATTACAGGGCTGCTATTAACAAGAGATGGTGGTAATAGCAACAATGAGTCCGAGATCTTCAGACCTGGAGGAGGAGATATGAGGGACAATTGGAGAAGTGAATTATATAAATATAAAGTAGTAAAAATTGAACCATTAGGAGTAGCACCCACCAAGGCAAAGAGAAGAGTGGTGCAGAGAGAAAAAAGAGCAGTGGGAATAGGAGCTTTGTTCCTTGGGTTCTTGGGAGCAGCAGGAAGCACTATGGGCGCAGCCTCAATGACGCTGACGGTACAGGCCAGACAATTATTGTCTGGTATAGTGCAGCAGCAGAACAATTTGCTGAGGGCTATTGAGGCGCAACAGCATCTGTTGCAACTCACAGTCTGGGGCATCAAGCAGCTCCAGGCAAGAATCCTGGCTGTGGAAAGATACCTAAAGGATCAACAGCTCCTGGGGATTTGGGGTTGCTCTGGAAAACTCATTTGCACCACTGCTGTGCCTTGGAATGCTAGTTGGAGTAATAAATCTCTGGAACAGATTTGGAATCACACGACCTGGATGGAGTGGGACAGAGAAATTAACAATTACACAAGCTTAATACACTCCTTAATTGAAGAATCGCAAAACCAGCAAGAAAAGAATGAACAAGAATTATTGGAATTAGATAAATGGGCAAGTTTGTGGAATTGGTTTAACATAACAAATTGGCTGTGGTATATAAAATTATTCATAATGATAGTAGGAGGCTTGGTAGGTTTAAGAATAGTTTTTGCTGTACTTTCTATAGTGAATAGAGTTAGGCAGGGATATTCACCATTATCGTTTCAGACCCACCTCCCAACCCCGAGGGGACCCGACAGGCCCGAAGGAATAGAAGAAGAAGGTGGAGAGAGAGACAGAGACAGATCCATTCGATTAGTGAACGGATCCTTGGCACTTATCTGGGACGATCTGCGGAGCCTGTGCCTCTTCAGCTACCACCGCTTGAGAGACTTACTCTTGATTGTAACGAGGATTGTGGAACTTCTGGGACGCAGGGGGTGGGAAGCCCTCAAATATTGGTGGAATCTCCTACAGTATTGGAGTCAGGAACTAAAGAATAGTGCTGTTAGCTTGCTCAATGCCACAGCCATAGCAGTAGCTGAGGGGACAGATAGGGTTATAGAAGTAGTACAAGGAGCTTGTAGAGCTATTCGCCACATACCTAGAAGAATAAGACAGGGCTTGGAAAGGATTTTGCTATAAGATGGGTGGCAAGTGGTCAAAAAGTAGTGTGATTGGATGGCCTACTGTAAGGGAAAGAATGAGACGAGCTGAGCCAGCAGCAGATAGGGTGGGAGCAGCATCTCGAGACCTGGAAAAACATGGAGCAATCACAAGTAGCAATACAGCAGCTACCAATGCTGCTTGTGCCTGGCTAGAAGCACAAGAGGAGGAGGAGGTGGGTTTTCCAGTCACACCTCAGGTACCTTTAAGACCAATGACTTACAAGGCAGCTGTAGATCTTAGCCACTTTTTAAAAGAAAAGGGGGGACTGGAAGGGCTAATTCACTCCCAAAGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGGGCCAGGGGTCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGATAAGATAGAAGAGGCCAATAAAGGAGAGAACACCAGCTTGTTACACCCTGTGAGCCTGCATGGGATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACGTGGCCCGAGAGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCGCTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGGGGAGTGGCGAGCCCTCAGATCCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA";
_HXB2_Env_Sequence_				= _HXB2_Sequence_ [_HXB_env_offset][_HXB_env_offset+2567];

_HXB2_AA_						= translateCodonToAA (_HXB2_Sequence_,defineCodonToAA(),2);
_HXB2_AA_ENV_					= translateCodonToAA (_HXB2_Env_Sequence_,defineCodonToAA(),0);

_hxb_alignOptions_nuc = {};
_hxb_alignOptions_nuc ["SEQ_ALIGN_CHARACTER_MAP"]="ACGT";
_hxb_alignOptions_nuc ["SEQ_ALIGN_SCORE_MATRIX"] = 	{
{5,-4,-4,-4}
{-4,5,-4,-4}
{-4,-4,5,-4}
{-4,-4,-4,5}
};
_hxb_alignOptions_nuc ["SEQ_ALIGN_GAP_OPEN"]	= 	50;
_hxb_alignOptions_nuc ["SEQ_ALIGN_GAP_OPEN2"]	= 	50;
_hxb_alignOptions_nuc ["SEQ_ALIGN_GAP_EXTEND"]	= 	1;
_hxb_alignOptions_nuc ["SEQ_ALIGN_GAP_EXTEND2"]	= 	1;
_hxb_alignOptions_nuc ["SEQ_ALIGN_AFFINE"]		=   1;
_hxb_alignOptions_nuc ["SEQ_ALIGN_NO_TP"]		=   1;


_hxb_alignOptions_prot = {};
_hxb_alignOptions_prot ["SEQ_ALIGN_CHARACTER_MAP"]		=   "ARNDCQEGHILKMFPSTWYV";
_hxb_alignOptions_prot ["SEQ_ALIGN_GAP_OPEN"]			= 	10;
_hxb_alignOptions_prot ["SEQ_ALIGN_GAP_OPEN2"]			= 	10;
_hxb_alignOptions_prot ["SEQ_ALIGN_GAP_EXTEND"]			= 	5;
_hxb_alignOptions_prot ["SEQ_ALIGN_GAP_EXTEND2"]		= 	5;
_hxb_alignOptions_prot ["SEQ_ALIGN_AFFINE"]				=   1;
_hxb_alignOptions_prot ["SEQ_ALIGN_NO_TP"]				=   1;
_hxb_alignOptions_prot ["SEQ_ALIGN_SCORE_MATRIX"]		=  
	{
	{                 7,                -7,                -7,                -4,               -10,               -11,                -4,                -3,               -10,                -6,                -9,                -9,                -7,               -13,                -3,                -2,                 1,               -16,               -15,                 0,                -5,                -5,                -3,               -17}
	{                -7,                 7,                -5,               -11,                -8,                -2,                -7,                -2,                 0,                -6,                -6,                 2,                -3,               -12,                -4,                -2,                -2,                -5,                -9,               -10,                -7,                -3,                -3,               -17}
	{                -7,                -5,                 8,                 2,                -9,                -6,                -6,                -7,                 0,                -6,               -12,                 0,               -10,               -12,                -9,                 1,                 0,               -17,                -3,               -10,                 6,                -6,                -3,               -17}
	{                -4,               -11,                 2,                 8,               -14,               -10,                 0,                -2,                -3,               -11,               -15,                -7,               -13,               -15,               -13,                -5,                -6,               -16,                -6,                -5,                 7,                 0,                -3,               -17}
	{               -10,                -8,                -9,               -14,                11,               -16,               -15,                -5,                -7,               -11,                -9,               -13,               -14,                 0,               -12,                -1,                -6,                -2,                 0,                -8,               -10,               -16,                -5,               -17}
	{               -11,                -2,                -6,               -10,               -16,                 8,                -2,               -10,                 0,               -12,                -4,                 0,                -8,               -12,                -1,                -9,                -8,               -14,                -9,               -13,                -7,                 6,                -4,               -17}
	{                -4,                -7,                -6,                 0,               -15,                -2,                 7,                -1,                -9,               -12,               -15,                -1,               -10,               -17,               -13,               -11,                -8,               -15,               -12,                -5,                 0,                 6,                -4,               -17}
	{                -3,                -2,                -7,                -2,                -5,               -10,                -1,                 7,               -10,               -11,               -14,                -6,               -12,                -9,               -11,                -1,                -7,                -5,               -14,                -5,                -4,                -3,                -4,               -17}
	{               -10,                 0,                 0,                -3,                -7,                 0,                -9,               -10,                10,               -10,                -4,                -5,               -10,                -6,                -3,                -6,                -6,               -11,                 2,               -14,                -1,                -2,                -3,               -17}
	{                -6,                -6,                -6,               -11,               -11,               -12,               -12,               -11,               -10,                 7,                 0,                -7,                 0,                -2,               -10,                -4,                 0,               -14,                -9,                 2,                -7,               -12,                -2,               -17}
	{                -9,                -6,               -12,               -15,                -9,                -4,               -15,               -14,                -4,                 0,                 6,               -10,                 0,                 0,                -3,                -5,                -8,                -6,                -8,                -4,               -13,                -6,                -4,               -17}
	{                -9,                 2,                 0,                -7,               -13,                 0,                -1,                -6,                -5,                -7,               -10,                 7,                -4,               -14,                -9,                -5,                -1,               -12,               -13,                -9,                -1,                -1,                -2,               -17}
	{                -7,                -3,               -10,               -13,               -14,                -8,               -10,               -12,               -10,                 0,                 0,                -4,                10,                -7,               -11,                -9,                -1,               -11,               -15,                 0,               -11,                -9,                -3,               -17}
	{               -13,               -12,               -12,               -15,                 0,               -12,               -17,                -9,                -6,                -2,                 0,               -14,                -7,                10,               -11,                -5,               -10,                -5,                 1,                -5,               -13,               -14,                -3,               -17}
	{                -3,                -4,                -9,               -13,               -12,                -1,               -13,               -11,                -3,               -10,                -3,                -9,               -11,               -11,                 8,                -1,                -3,               -13,               -11,               -12,               -10,                -3,                -5,               -17}
	{                -2,                -2,                 1,                -5,                -1,                -9,               -11,                -1,                -6,                -4,                -5,                -5,                -9,                -5,                -1,                 8,                 0,               -12,                -6,                -9,                 0,               -10,                -3,               -17}
	{                 1,                -2,                 0,                -6,                -6,                -8,                -8,                -7,                -6,                 0,                -8,                -1,                -1,               -10,                -3,                 0,                 7,               -16,               -10,                -4,                -2,                -8,                -2,               -17}
	{               -16,                -5,               -17,               -16,                -2,               -14,               -15,                -5,               -11,               -14,                -6,               -12,               -11,                -5,               -13,               -12,               -16,                10,                -4,               -16,               -16,               -14,                -8,               -17}
	{               -15,                -9,                -3,                -6,                 0,                -9,               -12,               -14,                 2,                -9,                -8,               -13,               -15,                 1,               -11,                -6,               -10,                -4,                10,               -12,                -4,               -10,                -4,               -17}
	{                 0,               -10,               -10,                -5,                -8,               -13,                -5,                -5,               -14,                 2,                -4,                -9,                 0,                -5,               -12,                -9,                -4,               -16,               -12,                 7,                -7,                -7,                -3,               -17}
	{                -5,                -7,                 6,                 7,               -10,                -7,                 0,                -4,                -1,                -7,               -13,                -1,               -11,               -13,               -10,                 0,                -2,               -16,                -4,                -7,                 7,                -2,                -4,               -17}
	{                -5,                -3,                -6,                 0,               -16,                 6,                 6,                -3,                -2,               -12,                -6,                -1,                -9,               -14,                -3,               -10,                -8,               -14,               -10,                -7,                -2,                 6,                -4,               -17}
	{                -3,                -3,                -3,                -3,                -5,                -4,                -4,                -4,                -3,                -2,                -4,                -2,                -3,                -3,                -5,                -3,                -2,                -8,                -4,                -3,                -4,                -4,                -3,               -17}
	{               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,                 1}
	} [{{0,0}}][{{19,19}}];


scoreMatrix = _hxb_alignOptions_prot ["SEQ_ALIGN_SCORE_MATRIX"];
LoadFunctionLibrary ("SeqAlignmentCodonShared", {"00": "HIV 25%", "01": "No penalty", "02": "First in file" });
_hxb_alignOptions_codon = alignOptions;


/*-------------------------------------------------------------*/
function mapSequenceToHXB2 (seq,option)
/* 
	option 0 - nucleotide alignment
    option 1 - amino-acid alignment
    option 2 - codon alignment
*/
{
 	if (option == 1)
 	{
		return mapSequenceToHXB2Aux (seq, _HXB2_AA_, option);
 	}
	return mapSequenceToHXB2Aux (seq, _HXB2_Sequence_, option);
}

/*-------------------------------------------------------------*/
function mapSequenceToHXB2Aux (seq,ref,option)
/* 
	option 0 - nucleotide alignment
    option 1 - amino-acid alignment
    option 2 - codon alignment
*/
{
	_seqLen	  = Abs(seq);
	_coordMap = {_seqLen,1};
	
	if (option != 1)
	{
		_inStr 		 = {{ref,seq}};
        if (option == 0)
        {
            AlignSequences(aligned, _inStr, _hxb_alignOptions_nuc);
        }
        else
        {
             AlignSequences(aligned, _inStr, _hxb_alignOptions_codon);
        }
    }
	else
	{
		_inStr 		 = {{ref,seq}};
		AlignSequences(aligned, _inStr, _hxb_alignOptions_prot);
	
	}
	
	_alignedHXB  = (aligned[0])[1];
	_alignedQRY  = (aligned[0])[2];
	
	_k				= (_alignedHXB$"^\\-+");
	_referenceSpan	= _k[1]+1;
	
	for (_k = 0; _k < _referenceSpan; _k = _k+1)
	{
		_coordMap[_k] = 0;
	}
	
	_qryCoord = _k;
	_refCoord = 0;

	while (_k < Abs(_alignedQRY))
	{
		if (_alignedQRY[_k] != "-")
		{
			_coordMap[_qryCoord] = _refCoord;
			_qryCoord = _qryCoord + 1;
		}
		if (_alignedHXB[_k] != "-")
		{
			_refCoord = _refCoord + 1;
		}
		_k = _k+1;
	}
	return _coordMap;
}

/*-------------------------------------------------------------*/
function selectHXB2subsequenceAux (seq,theSubset,mode)
{
    _template = {1,Abs(_HXB2_Sequence_)};
	_k2		 = Rows(theSubset)*Columns(theSubset);
	for (_k = 0; _k < _k2; _k = _k+1)
	{
		_span = _HXB_Annotation[theSubset[_k]];
		if (Abs(_span))
		{
			_template += _template["_MATRIX_ELEMENT_COLUMN_>=_span[0]&&_MATRIX_ELEMENT_COLUMN_<=_span[1]"];
		}
	}
	_mappedReference = mapSequenceToHXB2 (seq,0+2*mode);
	_subset			 = ""; _subset * 256;
	_k2 = Rows(_mappedReference);
	_k4 = Columns (_template);
	
    for (_k = 0; _k < _k2; _k = _k+1)
	{
		_k3 = _mappedReference[_k];
		if (_k3 >= _k4)
		{
			break;
		}
		if (_template[_k3])
		{
			_subset * seq[_k];
		}
	}
    
	_subset * 0;
	return _subset;

}

/*-------------------------------------------------------------*/
function selectHXB2subsequence (seq,theSubset)
{
    return selectHXB2subsequenceAux(seq,theSubset,0);
}

/*-------------------------------------------------------------*/
function selectHXB2subsequenceCodon (seq,theSubset)
{
    return selectHXB2subsequenceAux(seq,theSubset,1);
}

//--------------------------------------------------------------------------------

function		isoElectricPoint (seq) {
	COUNT_GAPS_IN_FREQUENCIES = 0;
	
	DataSet 			protSeq = ReadFromString ("$BASESET:BASE20\n>1\n" + seq);
	DataSetFilter		protFil = CreateFilter	 (protSeq,1);
	
	HarvestFrequencies (freqs,protFil,1,1,1);
	
	freqs = freqs*protFil.sites;
	
	expression = ""  + freqs[6 ] + "/(1+10^(pH-6.04))"  + /* H */
				 "+" + freqs[8 ] + "/(1+10^(pH-10.54))" + /* K */
				 "+" + freqs[14] + "/(1+10^(pH-12.48))" + /* R */
				 
				 "-" + freqs[2 ] + "/(1+10^(3.9-pH))"   + /* D */
				 "-" + freqs[3 ] + "/(1+10^(4.07-pH))"   + /* E */
				 "-" + freqs[1 ] + "/(1+10^(8.18-pH))"   + /* C */
				 "-" + freqs[19] + "/(1+10^(10.46-pH))"   ; /* Y */
	
	pH :> 0;
	pH :< 14;
	pH = 6.5;
	

	ExecuteCommands ("function ComputePI (pH){ return -Abs(`expression`); }");
	Optimize 		(res, ComputePI(pH));
		
	return res[0][0];
}


//--------------------------------------------------------------------------------

lfunction		countPNGS		(seq){
    pngs = seq || "N\\-*[^P]\\-*[ST]\\-*[^P]";
	return Rows(pngs)/2 - (pngs[0] < 0) ;
}

/*-------------------------------------------------------------*/
function selectHXB2ENVsubsequence (seq,theSubset, nucOrAA) {
	if (nucOrAA != 1)
	{
		_template = {1,Abs(_HXB2_Env_Sequence_)};
	}
	else
	{
		_template = {1,Abs(_HXB2_AA_ENV_)};	
	}
	
	_k2		  = Rows(theSubset)*Columns(theSubset);
	for (_k = 0; _k < _k2; _k = _k+1)
	{
		_span = _HXB_Annotation[theSubset[_k]] + (- _HXB_env_offset);
		
		if (Abs(_span))
		{
			if (nucOrAA == 0)
			{
				_template += _template["_MATRIX_ELEMENT_COLUMN_>=_span[0]&&_MATRIX_ELEMENT_COLUMN_<=_span[1]"];
			}
			else
			{
				_template += _template["_MATRIX_ELEMENT_COLUMN_>=_span[0]$3&&_MATRIX_ELEMENT_COLUMN_<=_span[1]$3"];			
			}
		}
	}
		
	if (nucOrAA != 1)
	{
		_mappedReference = mapSequenceToHXB2Aux (seq,_HXB2_Env_Sequence_,nucOrAA);
	}
	else
	{
		_mappedReference = mapSequenceToHXB2Aux (seq,_HXB2_AA_ENV_,nucOrAA);
	}
	
	_subset			 = ""; _subset * 256;
	
	_k2 = Rows(_mappedReference);
	_k4 = Columns (_template);
	for (_k = 0; _k < _k2; _k = _k+1)
	{
		_k3 = _mappedReference[_k];
		if (_k3 >= _k4)
		{
			break;
		}
		if (_template[_k3])
		{
			_subset * seq[_k];
		}
	}

	_subset * 0;
	return _subset;
}


/*-------------------------------------------------------------*/
function fractionalHXB2map (seq)
{
	_inStr 		 = {{_HXB2_AA_ENV_,seq}};
	
	AlignSequences(aligned, _inStr, _hxb_alignOptions_prot);
	
		
	_alignedHXB  = (aligned[0])[1];
	_alignedQRY  = (aligned[0])[2];
	
	_seqLen	  = Abs(seq);
	_coordMap = {_seqLen,1};
	
	
	_currentRegionHXBSpan	= 0;
	_currentRegionIndex		= 1;
	
	_currentHXB2Index		= 0;
	_currentQRYIndex		= 0;
	_k2						= Abs (_alignedHXB);
	
	for (_k = 0; _k < _k2; _k+=1)
	{
	
		//fprintf (stdout, _k, ":", _currentHXB2Index,_alignedHXB[_k], "|",  _alignedQRY[_k], _currentQRYIndex,  "\n");
			
		if (_alignedHXB[_k] != "-")
		{
			_currentHXB2Index += 1;
		}
		if (_alignedQRY[_k] != "-")
		{
			_currentQRYIndex += 1;
		}
		
		if (_currentHXB2Index == _HXB_aa_offset_matrix [_currentRegionIndex])
		{
		
			_qrySpan			 = _currentQRYIndex - _lastQRYRegion;
			_hxbSpan			 = _HXB_aa_offset_matrix [_currentRegionIndex] - _HXB_aa_offset_matrix [_currentRegionIndex-1];
			_scaler				 = (_hxbSpan-1)/(_qrySpan-1);
			_hxbSpan			 = _HXB_aa_offset_matrix [_currentRegionIndex-1];
			//fprintf (stdout, "In offset ",_currentRegionIndex,"\n", _scaler, ":", _hxbSpan, ":", _currentQRYIndex, ":", _lastQRYRegion, "\n");
			
			_currentRegionIndex += 1;	
			if (_currentRegionIndex == Columns(_HXB_aa_offset_matrix))
			{
				_k				 = _k2;
				_currentQRYIndex = Abs(seq);
			}
			
			for (_k3 = _lastQRYRegion; _k3 < _currentQRYIndex; _k3 += 1)
			{
				_coordMap [_k3] = _hxbSpan + (_k3-_lastQRYRegion) * _scaler;
				//fprintf (stdout, "\t", _k3, ":", _coordMap [_k3], "\n");
			}
			
			
			_lastQRYRegion       = _currentQRYIndex;
		}

	}
	
	return _coordMap;
}



