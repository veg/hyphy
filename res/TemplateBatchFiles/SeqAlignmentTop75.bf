/* START ALIGNMENT SETTINGS */

alignOptions = {};

alignOptions ["SEQ_ALIGN_CHARACTER_MAP"]="ARNDCQEGHILKMFPSTWYV";

ChoiceList (refSeq,"Scoring Matrix",1,SKIP_NONE,"BLOSUM62","Default BLAST BLOSUM62 matrix",
												"HIV 5%","Empirically derived 5% divergence HIV matrix",
												"HIV 25%","Empirically derived 25% divergence HIV matrix",
												"HIV 50%","Empirically derived 50% divergence HIV matrix");

if (refSeq < 0)
{
	return 0;
	
}
if (refSeq == 0)
{
	scoreMatrix = {
	{7,-3,-3,-3,-1,-2,-2,0,-3,-3,-3,-1,-2,-4,-1,2,0,-5,-4,-1}
	{-3,9,-1,-3,-6,1,-1,-4,0,-5,-4,3,-3,-5,-3,-2,-2,-5,-4,-4}
	{-3,-1,9,2,-5,0,-1,-1,1,-6,-6,0,-4,-6,-4,1,0,-7,-4,-5}
	{-3,-3,2,10,-7,-1,2,-3,-2,-7,-7,-2,-6,-6,-3,-1,-2,-8,-6,-6}
	{-1,-6,-5,-7,13,-5,-7,-6,-7,-2,-3,-6,-3,-4,-6,-2,-2,-5,-5,-2}
	{-2,1,0,-1,-5,9,3,-4,1,-5,-4,2,-1,-5,-3,-1,-1,-4,-3,-4}
	{-2,-1,-1,2,-7,3,8,-4,0,-6,-6,1,-4,-6,-2,-1,-2,-6,-5,-4}
	{0,-4,-1,-3,-6,-4,-4,9,-4,-7,-7,-3,-5,-6,-5,-1,-3,-6,-6,-6}
	{-3,0,1,-2,-7,1,0,-4,12,-6,-5,-1,-4,-2,-4,-2,-3,-4,3,-5}
	{-3,-5,-6,-7,-2,-5,-6,-7,-6,7,2,-5,2,-1,-5,-4,-2,-5,-3,4}
	{-3,-4,-6,-7,-3,-4,-6,-7,-5,2,6,-4,3,0,-5,-4,-3,-4,-2,1}
	{-1,3,0,-2,-6,2,1,-3,-1,-5,-4,8,-3,-5,-2,-1,-1,-6,-4,-4}
	{-2,-3,-4,-6,-3,-1,-4,-5,-4,2,3,-3,9,0,-4,-3,-1,-3,-3,1}
	{-4,-5,-6,-6,-4,-5,-6,-6,-2,-1,0,-5,0,10,-6,-4,-4,0,4,-2}
	{-1,-3,-4,-3,-6,-3,-2,-5,-4,-5,-5,-2,-4,-6,12,-2,-3,-7,-6,-4}
	{2,-2,1,-1,-2,-1,-1,-1,-2,-4,-4,-1,-3,-4,-2,7,2,-6,-3,-3}
	{0,-2,0,-2,-2,-1,-2,-3,-3,-2,-3,-1,-1,-4,-3,2,8,-5,-3,0}
	{-5,-5,-7,-8,-5,-4,-6,-6,-4,-5,-4,-6,-3,0,-7,-6,-5,16,3,-5}
	{-4,-4,-4,-6,-5,-3,-5,-6,3,-3,-2,-4,-3,4,-6,-3,-3,3,11,-3}
	{-1,-4,-5,-6,-2,-4,-4,-6,-5,4,1,-4,1,-2,-4,-3,0,-5,-3,7}
	};
}

if (refSeq == 1)
{
	scoreMatrix = 
{
{                11,                 0,                 2,                 2,                -2,                 3,                -1,                 1,                 0,                -1,                 1,                 0,                 3,                -2,                 0,                 3,                 7,                 5,                -5,                -4}
{                 0,                15,                -2,                -4,                 8,                 2,                 2,                -1,                -2,                 0,                -1,                 0,                -1,                -4,                 0,                 5,                 2,                 0,                 5,                 7}
{                 2,                -2,                12,                 6,                -3,                 4,                 4,                -2,                 1,                -5,                -2,                 7,                -2,                -1,                -1,                 2,                 1,                 2,                -4,                 2}
{                 2,                -4,                 6,                10,                -5,                 4,                 0,                -3,                 4,                -5,                 0,                 1,                -3,                 3,                 0,                -1,                 0,                 1,                -4,                -2}
{                -2,                 8,                -3,                -5,                14,                 0,                 2,                 4,                -3,                 5,                 2,                -1,                -1,                -2,                -2,                 2,                 0,                 2,                 3,                 8}
{                 3,                 2,                 4,                 4,                 0,                10,                -1,                -2,                 1,                -4,                -2,                 0,                -2,                -1,                 3,                 4,                 0,                 1,                 2,                -3}
{                -1,                 2,                 4,                 0,                 2,                -1,                15,                 0,                 2,                 2,                 0,                 6,                 4,                 6,                 6,                 2,                 2,                -3,                 0,                 9}
{                 1,                -1,                -2,                -3,                 4,                -2,                 0,                10,                 0,                 4,                 7,                 1,                -1,                -2,                 0,                 2,                 5,                 7,                -4,                 0}
{                 0,                -2,                 1,                 4,                -3,                 1,                 2,                 0,                11,                -2,                 3,                 5,                 0,                 5,                 7,                 2,                 4,                 0,                -2,                -2}
{                -1,                 0,                -5,                -5,                 5,                -4,                 2,                 4,                -2,                 9,                 5,                -3,                 2,                 2,                 0,                 1,                 0,                 2,                 1,                 0}
{                 1,                -1,                -2,                 0,                 2,                -2,                 0,                 7,                 3,                 5,                15,                 0,                 0,                 0,                 4,                 0,                 5,                 6,                 0,                -2}
{                 0,                 0,                 7,                 1,                -1,                 0,                 6,                 1,                 5,                -3,                 0,                12,                 0,                 1,                 2,                 7,                 5,                -1,                -5,                 4}
{                 3,                -1,                -2,                -3,                -1,                -2,                 4,                -1,                 0,                 2,                 0,                 0,                12,                 5,                 2,                 5,                 3,                -2,                -2,                -1}
{                -2,                -4,                -1,                 3,                -2,                -1,                 6,                -2,                 5,                 2,                 0,                 1,                 5,                11,                 4,                 0,                 0,                -3,                -3,                 0}
{                 0,                 0,                -1,                 0,                -2,                 3,                 6,                 0,                 7,                 0,                 4,                 2,                 2,                 4,                10,                 4,                 3,                -1,                 2,                 0}
{                 3,                 5,                 2,                -1,                 2,                 4,                 2,                 2,                 2,                 1,                 0,                 7,                 5,                 0,                 4,                11,                 6,                 0,                -2,                 2}
{                 7,                 2,                 1,                 0,                 0,                 0,                 2,                 5,                 4,                 0,                 5,                 5,                 3,                 0,                 3,                 6,                11,                 2,                -4,                 0}
{                 5,                 0,                 2,                 1,                 2,                 1,                -3,                 7,                 0,                 2,                 6,                -1,                -2,                -3,                -1,                 0,                 2,                11,                -5,                -2}
{                -5,                 5,                -4,                -4,                 3,                 2,                 0,                -4,                -2,                 1,                 0,                -5,                -2,                -3,                 2,                -2,                -4,                -5,                14,                 3}
{                -4,                 7,                 2,                -2,                 8,                -3,                 9,                 0,                -2,                 0,                -2,                 4,                -1,                 0,                 0,                 2,                 0,                -2,                 3,                14}
};

}

if (refSeq == 2)
{
	scoreMatrix = 
{
{                11,                -2,                 1,                 0,                -5,                 1,                -3,                -1,                -3,                -3,                -2,                -2,                 1,                -4,                -1,                 2,                 5,                 4,                -8,                -7}
{                -2,                16,                -6,                -7,                 6,                 0,                 0,                -4,                -6,                -3,                -5,                -2,                -4,                -7,                -1,                 4,                 0,                 0,                 3,                 6}
{                 1,                -6,                12,                 4,                -6,                 2,                 2,                -5,                -2,                -8,                -5,                 6,                -5,                -4,                -5,                 0,                -1,                 0,                -7,                 0}
{                 0,                -7,                 4,                10,                -8,                 2,                -2,                -6,                 2,                -8,                -2,                -1,                -6,                 1,                -2,                -5,                -2,                 0,                -7,                -4}
{                -5,                 6,                -6,                -8,                14,                -1,                 0,                 2,                -5,                 4,                 0,                -4,                -4,                -5,                -6,                 0,                -4,                 0,                 1,                 7}
{                 1,                 0,                 2,                 2,                -1,                10,                -4,                -6,                 0,                -8,                -5,                -1,                -5,                -4,                 2,                 2,                -1,                 0,                 0,                -7}
{                -3,                 0,                 2,                -2,                 0,                -4,                15,                -3,                 0,                 1,                -2,                 5,                 3,                 4,                 5,                 0,                 0,                -6,                -2,                 7}
{                -1,                -4,                -5,                -6,                 2,                -6,                -3,                10,                -1,                 2,                 5,                 0,                -4,                -6,                 0,                 0,                 4,                 5,                -7,                -2}
{                -3,                -6,                -2,                 2,                -5,                 0,                 0,                -1,                11,                -4,                 1,                 4,                -3,                 4,                 5,                 0,                 3,                -2,                -5,                -6}
{                -3,                -3,                -8,                -8,                 4,                -8,                 1,                 2,                -4,                 9,                 3,                -6,                 1,                 0,                -1,                 0,                -3,                 0,                 0,                -2}
{                -2,                -5,                -5,                -2,                 0,                -5,                -2,                 5,                 1,                 3,                15,                -3,                -4,                -1,                 2,                -3,                 4,                 4,                -2,                -6}
{                -2,                -2,                 6,                -1,                -4,                -1,                 5,                 0,                 4,                -6,                -3,                12,                -3,                 0,                 0,                 5,                 4,                -4,                -8,                 2}
{                 1,                -4,                -5,                -6,                -4,                -5,                 3,                -4,                -3,                 1,                -4,                -3,                12,                 3,                 0,                 3,                 1,                -5,                -4,                -4}
{                -4,                -7,                -4,                 1,                -5,                -4,                 4,                -6,                 4,                 0,                -1,                 0,                 3,                12,                 2,                -2,                -1,                -6,                -5,                -2}
{                -1,                -1,                -5,                -2,                -6,                 2,                 5,                 0,                 5,                -1,                 2,                 0,                 0,                 2,                11,                 2,                 2,                -4,                 0,                -3}
{                 2,                 4,                 0,                -5,                 0,                 2,                 0,                 0,                 0,                 0,                -3,                 5,                 3,                -2,                 2,                12,                 4,                -3,                -5,                 0}
{                 5,                 0,                -1,                -2,                -4,                -1,                 0,                 4,                 3,                -3,                 4,                 4,                 1,                -1,                 2,                 4,                11,                 0,                -7,                -3}
{                 4,                 0,                 0,                 0,                 0,                 0,                -6,                 5,                -2,                 0,                 4,                -4,                -5,                -6,                -4,                -3,                 0,                11,                -8,                -5}
{                -8,                 3,                -7,                -7,                 1,                 0,                -2,                -7,                -5,                 0,                -2,                -8,                -4,                -5,                 0,                -5,                -7,                -8,                14,                 2}
{                -7,                 6,                 0,                -4,                 7,                -7,                 7,                -2,                -6,                -2,                -6,                 2,                -4,                -2,                -3,                 0,                -3,                -5,                 2,                14}
};

}

if (refSeq == 3)
{
	scoreMatrix = 
{
{                10,                 1,                 4,                 3,                 0,                 4,                 0,                 3,                 1,                 0,                 3,                 3,                 4,                 0,                 2,                 5,                 7,                 6,                -2,                -1}
{                 1,                15,                 0,                -2,                 9,                 3,                 4,                 0,                 0,                 1,                 0,                 2,                 0,                -1,                 2,                 6,                 3,                 2,                 6,                 8}
{                 4,                 0,                12,                 7,                -1,                 5,                 5,                 0,                 3,                -3,                 0,                 8,                 0,                 1,                 0,                 4,                 3,                 3,                -2,                 3}
{                 3,                -2,                 7,                10,                -3,                 5,                 1,                 0,                 5,                -3,                 0,                 3,                 0,                 4,                 2,                 0,                 1,                 2,                -2,                 0}
{                 0,                 9,                -1,                -3,                13,                 1,                 5,                 5,                -1,                 6,                 4,                 0,                 1,                 0,                 0,                 3,                 1,                 3,                 4,                 9}
{                 4,                 3,                 5,                 5,                 1,                10,                 1,                 0,                 3,                -2,                 0,                 2,                 0,                 0,                 4,                 5,                 2,                 2,                 3,                -1}
{                 0,                 4,                 5,                 1,                 5,                 1,                14,                 1,                 4,                 3,                 2,                 7,                 5,                 7,                 7,                 4,                 3,                 0,                 1,                10}
{                 3,                 0,                 0,                 0,                 5,                 0,                 1,                 9,                 1,                 5,                 7,                 2,                 0,                 0,                 2,                 3,                 6,                 7,                -1,                 1}
{                 1,                 0,                 3,                 5,                -1,                 3,                 4,                 1,                10,                 0,                 4,                 6,                 1,                 6,                 8,                 3,                 5,                 1,                 0,                 0}
{                 0,                 1,                -3,                -3,                 6,                -2,                 3,                 5,                 0,                 9,                 6,                 0,                 4,                 3,                 1,                 2,                 1,                 3,                 2,                 2}
{                 3,                 0,                 0,                 0,                 4,                 0,                 2,                 7,                 4,                 6,                15,                 2,                 1,                 2,                 5,                 2,                 6,                 7,                 0,                 0}
{                 3,                 2,                 8,                 3,                 0,                 2,                 7,                 2,                 6,                 0,                 2,                11,                 2,                 3,                 4,                 7,                 6,                 0,                -2,                 5}
{                 4,                 0,                 0,                 0,                 1,                 0,                 5,                 0,                 1,                 4,                 1,                 2,                12,                 6,                 3,                 6,                 4,                 0,                -1,                 0}
{                 0,                -1,                 1,                 4,                 0,                 0,                 7,                 0,                 6,                 3,                 2,                 3,                 6,                11,                 5,                 1,                 2,                -1,                -1,                 2}
{                 2,                 2,                 0,                 2,                 0,                 4,                 7,                 2,                 8,                 1,                 5,                 4,                 3,                 5,                10,                 5,                 4,                 0,                 3,                 1}
{                 5,                 6,                 4,                 0,                 3,                 5,                 4,                 3,                 3,                 2,                 2,                 7,                 6,                 1,                 5,                11,                 6,                 1,                 0,                 3}
{                 7,                 3,                 3,                 1,                 1,                 2,                 3,                 6,                 5,                 1,                 6,                 6,                 4,                 2,                 4,                 6,                10,                 4,                -2,                 1}
{                 6,                 2,                 3,                 2,                 3,                 2,                 0,                 7,                 1,                 3,                 7,                 0,                 0,                -1,                 0,                 1,                 4,                10,                -2,                 0}
{                -2,                 6,                -2,                -2,                 4,                 3,                 1,                -1,                 0,                 2,                 0,                -2,                -1,                -1,                 3,                 0,                -2,                -2,                13,                 5}
{                -1,                 8,                 3,                 0,                 9,                -1,                10,                 1,                 0,                 2,                 0,                 5,                 0,                 2,                 1,                 3,                 1,                 0,                 5,                14}
};

}


alignOptions ["SEQ_ALIGN_SCORE_MATRIX"] = 	scoreMatrix;
alignOptions ["SEQ_ALIGN_GAP_OPEN"]		= 	40;
alignOptions ["SEQ_ALIGN_GAP_OPEN2"]	= 	20;
alignOptions ["SEQ_ALIGN_GAP_EXTEND"]	= 	10;
alignOptions ["SEQ_ALIGN_GAP_EXTEND2"]	= 	5;
alignOptions ["SEQ_ALIGN_AFFINE"]		=   1;

ChoiceList (refSeq,"Prefix/Suffix Indels",1,SKIP_NONE,"No penalty","Do not penalize prefix and suffix Indels","Normal penalty","Treat prefix and suffix indels as any other indels");
if (refSeq < 0)
{
	return 0;
}

alignOptions ["SEQ_ALIGN_NO_TP"]		=   1-refSeq;

/* END ALIGNMENT SETTINGS */


/* REFERENCE SEQUENCES OPTIONS */

predefSeqNames = {{"First in file", "Use the first sequence in the data file as a reference"}
				 /*0*/ {"HXB2_env", "Use HIV-1 HXB2 reference strain envelope sequence (K03455)"}
				 /*1*/ {"HXB2_nef", "Use HIV-1 HXB2 reference strain NEF sequence (K03455)"}
				 /*2*/ {"HXB2_gag", "Use HIV-1 HXB2 reference strain gag sequence (K03455)"} 
				 /*3*/ {"HXB2_vpr", "Use HIV-1 HXB2 reference strain vpr sequence (K03455)"} 
				 /*4*/ {"HXB2_vif", "Use HIV-1 HXB2 reference strain vif sequence (K03455)"} 
				 /*5*/ {"HXB2_vpu", "Use HIV-1 HXB2 reference strain vpu sequence (K03455)"} 
				 /*6*/ {"HXB2_pr", "Use HIV-1 HXB2 reference strain protease sequence (K03455)"} 
				 /*7*/ {"HXB2_rt", "Use HIV-1 HXB2 reference strain reverse transcriptase sequence (K03455)"} 
				 /*8*/ {"HXB2_int", "Use HIV-1 HXB2 reference strain integrase sequence (K03455)"} 
				 /*9*/ {"HXB2_rev", "Use HIV-1 HXB2 reference strain rev (exons 1 and 2)sequence (K03455)"} 
				/*10*/ {"HXB2_tat", "Use HIV-1 HXB2 reference strain tat (exons 1 and 2) sequence (K03455)"} 
				/*11*/ {"HXB2_prrt", "Use HIV-1 HXB2 reference strain protease+rt sequence (K03455)"}
				/*12*/ {"NL4_3prrt", "Use HIV-1 NL4-3 reference strain pr+rt sequence"} 
				  				  };
				  
				  
predefSeqNames2 = {{"No", "No reference coordinate sequences"}
				 /*0*/ {"HXB2_env", "Use HIV-1 HXB2 reference strain envelope sequence (K03455)"}
				 /*1*/ {"HXB2_nef", "Use HIV-1 HXB2 reference strain NEF sequence (K03455)"}
				 /*2*/ {"HXB2_gag", "Use HIV-1 HXB2 reference strain gag sequence (K03455)"} 
				 /*3*/ {"HXB2_vpr", "Use HIV-1 HXB2 reference strain vpr sequence (K03455)"} 
				 /*4*/ {"HXB2_vif", "Use HIV-1 HXB2 reference strain vif sequence (K03455)"} 
				 /*5*/ {"HXB2_vpu", "Use HIV-1 HXB2 reference strain vpu sequence (K03455)"} 
				 /*6*/ {"HXB2_pr", "Use HIV-1 HXB2 reference strain protease sequence (K03455)"} 
				 /*7*/ {"HXB2_rt", "Use HIV-1 HXB2 reference strain reverse transcriptase sequence (K03455)"} 
				 /*8*/ {"HXB2_int", "Use HIV-1 HXB2 reference strain integrase sequence (K03455)"} 
				 /*9*/ {"HXB2_rev", "Use HIV-1 HXB2 reference strain rev (exons 1 and 2)sequence (K03455)"} 
				/*10*/ {"HXB2_tat", "Use HIV-1 HXB2 reference strain tat (exons 1 and 2) sequence (K03455)"} 
				/*11*/ {"HXB2_prrt", "Use HIV-1 HXB2 reference strain protease+rt sequence (K03455)"} 
				/*12*/ {"NL4_3prrt", "Use HIV-1 NL4-3 reference strain pr+rt sequence"} 
				  };

RefSeqs = {};
RefSeqs [0] = "ATGAGAGTGAAGGAGAAATATCAGCACTTGTGGAGATGGGGGTGGAGATGGGGCACCATGCTCCTTGGGATGTTGATGATCTGTAGTGCTACAGAAAAATTGTGGGTCACAGTCTATTATGGGGTACCTGTGTGGAAGGAAGCAACCACCACTCTATTTTGTGCATCAGATGCTAAAGCATATGATACAGAGGTACATAATGTTTGGGCCACACATGCCTGTGTACCCACAGACCCCAACCCACAAGAAGTAGTATTGGTAAATGTGACAGAAAATTTTAACATGTGGAAAAATGACATGGTAGAACAGATGCATGAGGATATAATCAGTTTATGGGATCAAAGCCTAAAGCCATGTGTAAAATTAACCCCACTCTGTGTTAGTTTAAAGTGCACTGATTTGAAGAATGATACTAATACCAATAGTAGTAGCGGGAGAATGATAATGGAGAAAGGAGAGATAAAAAACTGCTCTTTCAATATCAGCACAAGCATAAGAGGTAAGGTGCAGAAAGAATATGCATTTTTTTATAAACTTGATATAATACCAATAGATAATGATACTACCAGCTATAAGTTGACAAGTTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAGCCAATTCCCATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAATGTAATAATAAGACGTTCAATGGAACAGGACCATGTACAAATGTCAGCACAGTACAATGTACACATGGAATTAGGCCAGTAGTATCAACTCAACTGCTGTTAAATGGCAGTCTAGCAGAAGAAGAGGTAGTAATTAGATCTGTCAATTTCACGGACAATGCTAAAACCATAATAGTACAGCTGAACACATCTGTAGAAATTAATTGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGAGAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGTAACATTAGTAGAGCAAAATGGAATAACACTTTAAAACAGATAGCTAGCAAATTAAGAGAACAATTTGGAAATAATAAAACAATAATCTTTAAGCAATCCTCAGGAGGGGACCCAGAAATTGTAACGCACAGTTTTAATTGTGGAGGGGAATTTTTCTACTGTAATTCAACACAACTGTTTAATAGTACTTGGTTTAATAGTACTTGGAGTACTGAAGGGTCAAATAACACTGAAGGAAGTGACACAATCACCCTCCCATGCAGAATAAAACAAATTATAAACATGTGGCAGAAAGTAGGAAAAGCAATGTATGCCCCTCCCATCAGTGGACAAATTAGATGTTCATCAAATATTACAGGGCTGCTATTAACAAGAGATGGTGGTAATAGCAACAATGAGTCCGAGATCTTCAGACCTGGAGGAGGAGATATGAGGGACAATTGGAGAAGTGAATTATATAAATATAAAGTAGTAAAAATTGAACCATTAGGAGTAGCACCCACCAAGGCAAAGAGAAGAGTGGTGCAGAGAGAAAAAAGAGCAGTGGGAATAGGAGCTTTGTTCCTTGGGTTCTTGGGAGCAGCAGGAAGCACTATGGGCGCAGCCTCAATGACGCTGACGGTACAGGCCAGACAATTATTGTCTGGTATAGTGCAGCAGCAGAACAATTTGCTGAGGGCTATTGAGGCGCAACAGCATCTGTTGCAACTCACAGTCTGGGGCATCAAGCAGCTCCAGGCAAGAATCCTGGCTGTGGAAAGATACCTAAAGGATCAACAGCTCCTGGGGATTTGGGGTTGCTCTGGAAAACTCATTTGCACCACTGCTGTGCCTTGGAATGCTAGTTGGAGTAATAAATCTCTGGAACAGATTTGGAATCACACGACCTGGATGGAGTGGGACAGAGAAATTAACAATTACACAAGCTTAATACACTCCTTAATTGAAGAATCGCAAAACCAGCAAGAAAAGAATGAACAAGAATTATTGGAATTAGATAAATGGGCAAGTTTGTGGAATTGGTTTAACATAACAAATTGGCTGTGGTATATAAAATTATTCATAATGATAGTAGGAGGCTTGGTAGGTTTAAGAATAGTTTTTGCTGTACTTTCTATAGTGAATAGAGTTAGGCAGGGATATTCACCATTATCGTTTCAGACCCACCTCCCAACCCCGAGGGGACCCGACAGGCCCGAAGGAATAGAAGAAGAAGGTGGAGAGAGAGACAGAGACAGATCCATTCGATTGGTGAACGGATCCTTGGCACTTATCTGGGACGATCTGCGGAGCCTGTGCCTCTTCAGCTACCACCGCTTGAGAGACTTACTCTTGATTGTAACGAGGATTGTGGAACTTCTGGGACGCAGGGGGTGGGAAGCCCTCAAATATTGGTGGAATCTCCTACAGTATTGGAGTCAGGAACTAAAGAATAGTGCTGTTAGCTTGCTCAATGCCACAGCCATAGCAGTAGCTGAGGGGACAGATAGGGTTATAGAAGTAGTACAAGGAGCTTGTAGAGCTATTCGCCACATACCTAGAAGAATAAGACAGGGCTTGGAAAGGATTTTGCTA";
RefSeqs [1] = "ATGGGTGGCAAGTGGTCAAAAAGTAGTGTGATTGGATGGCCTACTGTAAGGGAAAGAATGAGACGAGCTGAGCCAGCAGCAGATAGGGTGGGAGCAGCATCTCGAGACCTGGAAAAACATGGAGCAATCACAAGTAGCAATACAGCAGCTACCAATGCTGCTTGTGCCTGGCTAGAAGCACAAGAGGAGGAGGAGGTGGGTTTTCCAGTCACACCTCAGGTACCTTTAAGACCAATGACTTACAAGGCAGCTGTAGATCTTAGCCACTTTTTAAAAGAAAAGGGGGGACTGGAAGGGCTAATTCACTCCCAAAGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGAT---CAGAACTACACACCAGGGCCAGGGGTCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGATAAGATAGAAGAGGCCAATAAAGGAGAGAACACCAGCTTGTTACACCCTGTGAGCCTGCATGGGATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACGTGGCCCGAGAGCTGCATCCGGAGTACTTCAAGAACTGC";
RefSeqs [2] = "ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATATAGTATGGGCAAGCAGGGAGCTAGAACGATTCGCAGTTAATCCTGGCCTGTTAGAAACATCAGAAGGCTGTAGACAAATACTGGGACAGCTACAACCATCCCTTCAGACAGGATCAGAAGAACTTAGATCATTATATAATACAGTAGCAACCCTCTATTGTGTGCATCAAAGGATAGAGATAAAAGACACCAAGGAAGCTTTAGACAAGATAGAGGAAGAGCAAAACAAAAGTAAGAAAAAAGCACAGCAAGCAGCAGCTGACACAGGACACAGCAATCAGGTCAGCCAAAATTACCCTATAGTGCAGAACATCCAGGGGCAAATGGTACATCAGGCCATATCACCTAGAACTTTAAATGCATGGGTAAAAGTAGTAGAAGAGAAGGCTTTCAGCCCAGAAGTGATACCCATGTTTTCAGCATTATCAGAAGGAGCCACCCCACAAGATTTAAACACCATGCTAAACACAGTGGGGGGACATCAAGCAGCCATGCAAATGTTAAAAGAGACCATCAATGAGGAAGCTGCAGAATGGGATAGAGTGCATCCAGTGCATGCAGGGCCTATTGCACCAGGCCAGATGAGAGAACCAAGGGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACAAATAATCCACCTATCCCAGTAGGAGAAATTTATAAAAGATGGATAATCCTGGGATTAAATAAAATAGTAAGAATGTATAGCCCTACCAGCATTCTGGACATAAGACAAGGACCAAAGGAACCCTTTAGAGACTATGTAGACCGGTTCTATAAAACTCTAAGAGCCGAGCAAGCTTCACAGGAGGTAAAAAATTGGATGACAGAAACCTTGTTGGTCCAAAATGCGAACCCAGATTGTAAGACTATTTTAAAAGCATTGGGACCAGCGGCTACACTAGAAGAAATGATGACAGCATGTCAGGGAGTAGGAGGACCCGGCCATAAGGCAAGAGTTTTGGCTGAAGCAATGAGCCAAGTAACAAATTCAGCTACCATAATGATGCAGAGAGGCAATTTTAGGAACCAAAGAAAGATTGTTAAGTGTTTCAATTGTGGCAAAGAAGGGCACACAGCCAGAAATTGCAGGGCCCCTAGGAAAAAGGGCTGTTGGAAATGTGGAAAGGAAGGACACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAAGATCTGGCCTTCCTACAAGGGAAGGCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAA";				  
RefSeqs [3] = "ATGGAACAAGCCCCAGAAGACCAAGGGCCACAGAGGGAGCCACACAATGAATGGACACTAGAGCTTTTAGAGGAGCTTAAGAATGAAGCTGTTAGACATTTTCCTAGGATTTGGCTCCATGGCTTAGGGCAACATATCTATGAAACTTATGGGGATACTTGGGCAGGAGTGGAAGCCATAATAAGAATTCTGCAACAACTGCTGTTTATCCATTTTCAGAATTGGGTGTCGACATAGCAGAATAGGCGTTACTCGACAGAGGAGAGCAAGAAATGGAGCCAGTAGATCC";				  
RefSeqs [4] = "ATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGATTAGAACATGGAAAAGTTTAGTAAAACACCATATGTATGTTTCAGGGAAAGCTAGGGGATGGTTTTATAGACATCACTATGAAAGCCCTCATCCAAGAATAAGTTCAGAAGTACACATCCCACTAGGGGATGCTAGATTGGTAATAACAACATATTGGGGTCTGCATACAGGAGAAAGAGACTGGCATTTGGGTCAGGGAGTCTCCATAGAATGGAGGAAAAAGAGATATAGCACACAAGTAGACCCTGAACTAGCAGACCAACTAATTCATCTGTATTACTTTGACTGTTTTTCAGACTCTGCTATAAGAAAGGCCTTATTAGGACACATAGTTAGCCCTAGGTGTGAATATCAAGCAGGACATAACAAGGTAGGATCTCTACAATACTTGGCACTAGCAGCATTAATAACACCAAAAAAGATAAAGCCACCTTTGCCTAGTGTTACGAAACTGACAGAGGATAGATGGAACAAGCCCCAGAAGACCAAGGGCCACAGAGGGAGCCACACAATGAATGGACAC";				  
RefSeqs [5] = "ACGCAACCTATACCAATAGTAGCAATAGTAGCATTAGTAGTAGCAATAATAATAGCAATAGTTGTGTGGTCCATAGTAATCATAGAATATAGGAAAATATTAAGACAAAGAAAAATAGACAGGTTAATTGATAGACTAATAGAAAGAGCAGAAGACAGTGGCAATGAGAGTGAAGGAGAAATATCAGCACTTGTGGAGATGGGGGTGGAGATGGGGCACCATGCTCCTTGGGATGTTGATGATC";				  
RefSeqs [6] = "CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT";				  
RefSeqs [7] = "CCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAACTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAGGCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAACATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCCTGAGTGGGAGTTTGTTAATACCCCTCCCTTAGTGAAATTATGGTACCAGTTAGAGAAAGAACCCATAGTAGGAGCAGAAACCTTC";				  
RefSeqs [8] = "TTTTTAGATGGAATAGATAAGGCCCAAGATGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAGTGATTTTAACCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGAGAAGCCATGCATGGACAAGTAGACTGTAGTCCAGGAATATGGCAACTAGATTGTACACATTTAGAAGGAAAAGTTATCCTGGTAGCAGTTCATGTAGCCAGTGGATATATAGAAGCAGAAGTTATTCCAGCAGAAACAGGGCAGGAAACAGCATATTTTCTTTTAAAATTAGCAGGAAGATGGCCAGTAAAAACAATACATACTGACAATGGCAGCAATTTCACCGGTGCTACGGTTAGGGCCGCCTGTTGGTGGGCGGGAATCAAGCAGGAATTTGGAATTCCCTACAATCCCCAAAGTCAAGGAGTAGTAGAATCTATGAATAAAGAATTAAAGAAAATTATAGGACAGGTAAGAGATCAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAAATCCACTTTGGAAAGGACCAGCAAAGCTCCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGATAATAGTGACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATTAGGGATTATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGAT";
RefSeqs [9] = "ATGGCAGGAAGAAGCGGAGACAGCGACGAAGAGCTCATCAGAACAGTCAGACTCATCAAGCTTCTCTAACCCACCTCCCAACCCCGAGGGGACCCGACAGGCCCGAAGGAATAGAAGAAGAAGGTGGAGAGAGAGACAGAGACAGATCCATTCGATTAGTGAACGGATCCTTGGCACTTATCTGGGACGATCTGCGGAGCCTGTGCCTCTTCAGCTACCACCGCTTGAGAGACTTACTCTTGATTGTAACGAGGATTGTGGAACTTCTGGGACGCAGGGGGTGGGAAGCCCTCAAATATTGGTGGAATCTCCTACAGTATTGGAGTCAGGAACTAAAGA";
RefSeqs[10] = "ATGGAGCCAGTAGATCCTAGACTAGAGCCCTGGAAGCATCCAGGAAGTCAGCCTAAAACTGCTTGTACCAATTGCTATTGTAAAAAGTGTTGCTTTCATTGCCAAGTTTGTTTCATAACAAAAGCCTTAGGCATCTCCTATGGCAGGAAGAAGCGGAGACAGCGACGAAGAGCTCATCAGAACAGTCAGACTCATCAAGCTTCTCTATCAAAGCAACCCACCTCCCAACCCCGAGGGGACCCGACAGGCCCGAAGGAATAGAAGAAGAAGGTGGAGAGAGAGACAGAGACAGATCCATTCGAT";
RefSeqs[11] = "CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTTCCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAACTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAGGCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAACATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCCTGAGTGGGAGTTTGTTAATACCCCTCCCTTAGTGAAATTATGGTACCAGTTAGAGAAAGAACCCATAGTAGGAGCAGAAACCTTC";
RefSeqs[12] = "CCTCAGATCACTCTTTGGCAGCGACCCCTCGTCACAATAAAGATAGGGGGGCAATTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAGGACAGTATGATCAGATACTCATAGAAATCTGCGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGCTGCACTTTAAATTTTCCCATTAGTCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAAATGGAAAAGGAAGGAAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGATTTCTGGGAAGTTCAATTAGGAATACCACATCCTGCAGGGTTAAAACAGAAAAAATCAGTAACAGTACTGGATGTGGGCGATGCATATTTTTCAGTTCCCTTAGATAAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAGTGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTCATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGACAACATCTGTTGAGGTGGGGATTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAGGACAGCTGGACTGTCAATGACATACAGAAATTAGTGGGAAAATTGAATTGGGCAAGTCAGATTTATGCAGGGATTAAAGTAAGGCAATTATGTAAACTTCTTAGGGGAACCAAAGCACTAACAGAAGTAGTACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGGGAGATTCTAAAAGAACCGGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAAGGGTGCCCACACTAATGATGTGAAACAATTAACAGAGGCAGTACAAAAAATAGCCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAATTACCCATACAAAAGGAAACATGGGAAGCATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCCTGAGTGGGAGTTTGTCAATACCCCTCCCTTAGTGAAGTTATGGTACCAGTTAGAGAAAGAACCCATAATAGGAGCAGAAACTTTC";

ChoiceList (refSeq,"Choose a reference sequence",1,SKIP_NONE,predefSeqNames);
if (refSeq < 0)
{
	return 0;
}

SetDialogPrompt ("Sequence File:");

if (refSeq > 0)
{
	DataSet        unal2 			= ReadDataFile 	(PROMPT_FOR_FILE);
	refSeq = ">" + predefSeqNames [refSeq][0] + "\n" + RefSeqs[refSeq-1];
	DataSet		   refD = ReadFromString (refSeq);
	DataSet        unal = Combine (refD,unal2);
}
else
{
	
	ChoiceList (refSeq2,"Insert a coordinate reference sequence?",1,SKIP_NONE,predefSeqNames2);
	if (refSeq2 < 0)
	{
		return 0;
	}
	if (refSeq2)
	{
		DataSet        unal2 			= ReadDataFile 	(PROMPT_FOR_FILE);
		refSeq = ">" + predefSeqNames [refSeq2][0] + "\n" + RefSeqs[refSeq2-1];
		DataSet		   refD = ReadFromString (refSeq);
		DataSet        unal = Combine (unal2,refD);
	}
	else
	{
		DataSet        unal 			= ReadDataFile 	(PROMPT_FOR_FILE);
	}
}

ChoiceList (skipOutliers,"Include outliers?",1,SKIP_NONE,"No","Skip sequences with unusually poor alignment scores","Yes","Include all alignable sequences");

if (skipOutliers<0)
{
	return 0;
}

DataSetFilter  filteredData 	= CreateFilter	(unal,1);

GetInformation (UnalignedSeqs,filteredData);
/* preprocess sequences */

unalSequenceCount = Rows(UnalignedSeqs)*Columns(UnalignedSeqs);

for (seqCounter = 0; seqCounter < unalSequenceCount; seqCounter = seqCounter+1)
{
	aSeq = UnalignedSeqs[seqCounter];
	/*for (siteCounter = Abs(aSeq)-1; siteCounter>0; siteCounter=siteCounter-1)
	{
		if (aSeq[siteCounter]!="?")
		{
			break;
		}
	}
	UnalignedSeqs[seqCounter] = aSeq[0][siteCounter];*/
	UnalignedSeqs[seqCounter] = aSeq^{{"[^a-zA-Z]",""}};
	UnalignedSeqs[seqCounter] = UnalignedSeqs[seqCounter]^{{"^N+",""}};
	UnalignedSeqs[seqCounter] = UnalignedSeqs[seqCounter]^{{"N+$",""}};
}

/* build codon translation table */

incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def";
ExecuteCommands  ("#include \""+incFileName+"\";");

codonToAAMap = {};
codeToAA 	 = "FLIMVSPTAYXHQNKDECWRG";

nucChars = "ACGT";

for (p1=0; p1<64; p1=p1+1)
{
	codon = nucChars[p1$16]+nucChars[p1%16$4]+nucChars[p1%4];
	ccode = _Genetic_Code[p1];
	codonToAAMap[codon] = codeToAA[ccode];
}

/* determine reading frames	*/
ProteinSequences = {};
AllTranslations  = {};
ReadingFrames	 = {};
StopCodons		 = {};
StopPositions    = {};

fprintf (stdout, "\nDetecting reading frames for each sequence...\n");
frameCounter  = {3,1};
stillHasStops = {};

for (seqCounter = 0; seqCounter < unalSequenceCount; seqCounter = seqCounter+1)
{
	aSeq = UnalignedSeqs[seqCounter];
	seqLen = Abs(aSeq)-2;
	
	minStops = 1e20;
	tString = "";
	rFrame = 0;
	
	stopPosn = {3,2};
	allTran  = {3,1};
	for (offset = 0; offset < 3; offset = offset+1)
	{
		translString = "";
		translString * (seqLen/3+1);
		for (seqPos = offset; seqPos < seqLen; seqPos = seqPos+3)
		{
			codon = aSeq[seqPos][seqPos+2];
			prot = codonToAAMap[codon];
			if (Abs(prot))
			{
				translString * prot;
			}
			else
			{
				translString * "?";
			}
		} 
		translString * 0;
		translString = translString^{{"X$","?"}};
		stopPos = translString||"X";
		if (stopPos[0]>=0)
		{
			stopCount = Rows(stopPos)$2;
			stopPosn[offset][0] = stopPos[0];
			stopPosn[offset][1] = stopPos[stopCount*2-1];
		}
		else
		{
			stopCount = 0;
		}
		if (stopCount<minStops)
		{
			minStops = stopCount;
			rFrame = offset;
			tString = translString;
		}
		allTran[offset] = translString;
	}
	ReadingFrames[seqCounter] 		= rFrame;
	ProteinSequences[seqCounter]	= tString;
	frameCounter[rFrame] 			= frameCounter[rFrame]+1;
	StopPositions[seqCounter]		= stopPosn;
	AllTranslations [seqCounter]	= allTran;
	
	if (minStops>0)
	{
		stillHasStops[Abs(stillHasStops)] = seqCounter;
		if (seqCounter == 0)
		{
			fprintf (stdout, "Reference sequence must not contain frameshifts\n");
			return 0;
		}
	}
}

s1 = ProteinSequences[0];

fprintf (stdout, "\nFound:\n\t", frameCounter[0], " sequences in reading frame 1\n\t",frameCounter[1], " sequences in reading frame 2\n\t",frameCounter[2], " sequences in reading frame 3\n\nThere were ", Abs(stillHasStops), " sequences with apparent frameshift/sequencing errors\n");

skipSeqs = {};

for (k=0; k<Abs(stillHasStops); k=k+1)
{
	seqCounter = stillHasStops[k];
	GetString (seqName, unal, seqCounter);
	fprintf (stdout,"Sequence ", seqCounter+1, " (", seqName, ") seems to have");
	stopPosn = StopPositions[seqCounter];
	
	fStart = -1;
	fEnd   = -1;
	fMin   = 1e10;
	frame1 = 0;
	frame2 = 0;
	
	checkFramePosition (stopPosn[0][1],stopPosn[1][0],0,1);
	checkFramePosition (stopPosn[1][1],stopPosn[0][0],1,0);
	checkFramePosition (stopPosn[0][1],stopPosn[2][0],0,2);
	checkFramePosition (stopPosn[2][1],stopPosn[0][0],2,0);
	checkFramePosition (stopPosn[2][1],stopPosn[1][0],2,1);
	checkFramePosition (stopPosn[1][1],stopPosn[2][0],1,2);
	
	if (fStart>=0)
	{
		allTran = AllTranslations[seqCounter];
		useq    				   = UnalignedSeqs[seqCounter];
		fprintf (stdout, " a shift from frame ", frame2+1, " to frame ", frame1+1, " between a.a. positions ", fStart, " and ", fEnd, ".");
		fStart2 = Max(fStart-1,0);
		fEnd2   = Min(fEnd+1,Min(Abs(allTran[frame1]),Abs(allTran[frame2]))-1);
		tempString = allTran[frame2];
		fprintf (stdout, "\n\tRegion ", fStart2, "-", fEnd2, " in frame  ", frame2+1, ":\n\t", tempString[fStart2][fEnd2]);
		fprintf (stdout, "\n\t", useq[3*fStart2+frame2][3*fEnd2+frame2-1]);
		tempString = allTran[frame1];
		fprintf (stdout, "\n\tRegion ", fStart2, "-", fEnd2, " in frame  ", frame1+1, ":\n\t", tempString[fStart2][fEnd2]);
		fprintf (stdout, "\n\t", useq[3*fStart2+frame1][3*fEnd2+frame1-1]);
		fprintf (stdout, "\n\t\tAttempting to resolve by alignment to reference. ");
		
		f1s = allTran[frame1];
		f2s = allTran[frame2];
		f1l = Abs(f1s);
		
		bestScore  = -1e10;
		bestSplice = -1;
		
		for (k2=fStart; k2<fEnd; k2=k2+1)
		{
			s2 = f2s[0][k2]+f1s[k2+1][Abs(f1s)];
			inStr = {{s1,s2}};
			AlignSequences(aligned, inStr, alignOptions);
			aligned = aligned[0];
			aligned = aligned[0];
			if (aligned > bestScore)
			{
				bestScore = aligned;
				bestSplice = k2;
				bestString = s2;
			}
		}
		fprintf (stdout, "Best splice site appears to be at a.a. position ", bestSplice, "\n");
		/* update best spliced string */
		
		ProteinSequences[seqCounter] = bestString;
		ReadingFrames[seqCounter]    = 0;
		
		UnalignedSeqs[seqCounter]  = useq[frame2][frame2+3*bestSplice+2] + useq[frame1+3*bestSplice+3][Abs(useq)-1] + "---";
	}
	else
	{
		
		fprintf (stdout, " multiple frameshifts\n");
		skipSeqs[seqCounter] = 1;
	}	
}

SeqAlignments 	 = {};
startingPosition = {unalSequenceCount,2};
refLength = Abs(ProteinSequences[0]);
refInsertions	 = {refLength,1};

fprintf (stdout,"\nPerforming pairwise alignment with reference sequences\n");

alignmentScores = {unalSequenceCount-1,2};

for (seqCounter = 1; seqCounter < unalSequenceCount; seqCounter = seqCounter+1)
{
	if (skipSeqs[seqCounter] == 0)
	{
		s2 			 = ProteinSequences[seqCounter];
		inStr 		 = {{s1,s2}};
		AlignSequences(aligned, inStr, alignOptions);
		aligned = aligned[0];
		SeqAlignments[seqCounter] = aligned;
		alignmentScores[seqCounter-1][0] = aligned[0];
		alignmentScores[seqCounter-1][1] = seqCounter;
		aligned = aligned[1];
		myStartingPosition = aligned$"[^-]";
		myEndingPosition  = Abs (aligned)-1;
		while (aligned[myEndingPosition]=="-")
		{
			myEndingPosition = myEndingPosition - 1;
		}
		myStartingPosition = myStartingPosition[0];
		startingPosition[seqCounter][0] = myStartingPosition;
		startingPosition[seqCounter][1] = myEndingPosition;
		aligned = aligned[myStartingPosition][myEndingPosition];
		refInsert = aligned||"-+";
		if (refInsert[0]>0)
		{
			insCount = Rows (refInsert)/2;
			offset = 0;
			for (insN = 0; insN < insCount; insN = insN+1)
			{
				insPos 		= refInsert[insN*2];
				insLength	= refInsert[insN*2+1]-insPos+1;
				insPos 		= insPos-offset;
				if (refInsertions[insPos]<insLength)
				{
					refInsertions[insPos]=insLength;
				}
				offset = offset + insLength;
			}
		}
	}
}

alignmentScores = alignmentScores%0;

/* produce a fully gapped reference sequence */

fprintf (stdout,"\nMerging pairwise alignments into a MSA\n");

fullRefSeq = "";
fullRefSeq * refLength;
fullRefSeq * (s1[0]);


s1N = UnalignedSeqs[0];

fullRefSeqN = "";
fullRefSeqN * (3*refLength);
fullRefSeqN * (s1N[0][2]);

frameShift = ReadingFrames[0];

for (seqCounter=1;seqCounter<refLength;seqCounter=seqCounter+1)
{
	gapCount = refInsertions[seqCounter];
	for (k=0; k<gapCount;k=k+1)
	{
		fullRefSeq*("-");
		fullRefSeqN*("---");
	}	
	fullRefSeq  * (s1[seqCounter]);
	fullRefSeqN * (s1N[frameShift+seqCounter*3][frameShift+seqCounter*3+2]);
}

fullRefSeq * 0;
fullRefSeqN * 0;

refLength = Abs(fullRefSeq);

SetDialogPrompt ("Save alignment to:");

GetString (seqName,unal,0);
fprintf (PROMPT_FOR_FILE,CLEAR_FILE,">",seqName,"\n",fullRefSeq);
fName = LAST_FILE_PATH;
fNameC = fName+".nuc";
fprintf (fNameC,CLEAR_FILE,">",seqName,"\n",fullRefSeqN);


/*for (seqCounter = 1; seqCounter < unalSequenceCount; seqCounter = seqCounter+1)*/
for (seqID = 0; seqID < 75; seqID = seqID + 1)
{
	seqCounter = alignmentScores[unalSequenceCount-2-seqID][1];
	if (skipSeqs[seqCounter] == 0)
	{
		GetString (seqName,unal,seqCounter);
		aligned = SeqAlignments[seqCounter];
		
		aligned1 = aligned[1];
		aligned2 = aligned[2];
		
		s2 = startingPosition[seqCounter][0];
		e2 = startingPosition[seqCounter][1];
		
		gappedSeq = "";
		gappedSeq * Abs(aligned2);

		
		k=0;
		
		while (k<refLength)
		{
			while (fullRefSeq[k]!=aligned1[s2])
			{
				gappedSeq*("-");
				k=k+1;
			}
			gappedSeq*(aligned2[s2]);
			s2=s2+1;
			k=k+1;
		}

		gappedSeq * 0;

		gappedSeqN = "";
		gappedSeqN * (3*Abs(aligned2));
		
		frameShift = ReadingFrames[seqCounter];

		s1N 	= UnalignedSeqs[seqCounter];
		s2N		= ProteinSequences[seqCounter];
		s2 		= startingPosition[seqCounter][0];
		k 		= 0;
		e2		= Abs(gappedSeq);
		k = 0;
		while  (k<e2)
		{
			while ((s2N[s2]!=gappedSeq[k])&&(k<e2))
			{
				gappedSeqN * ("---");
				k=k+1;
			}
			if (k<e2)
			{
				gappedSeqN * s1N[frameShift+s2*3][frameShift+s2*3+2];
				s2 = s2+1;
				k=k+1;
			}
		}
		gappedSeqN * 0;

		if (refSeq2 && seqCounter == unalSequenceCount-1)
		{
			fscanf (fName, "Raw", soFar);
			fprintf (fName, CLEAR_FILE,">",seqName,"\n",gappedSeq,"\n",soFar);
			fscanf (fNameC, "Raw", soFar);
			fprintf (fNameC,CLEAR_FILE,">",seqName,"\n",gappedSeqN,"\n",soFar);		
			
		}
		else
		{
			fprintf (fName,"\n>",seqName,"\n",gappedSeq);
			fprintf (fNameC,"\n>",seqName,"\n",gappedSeqN);		
		}
	}
}

if (Abs(skipSeqs))
{
	fName = fName+".bad";
	for (seqCounter = 1; seqCounter < unalSequenceCount; seqCounter = seqCounter+1)
	{
		if (skipSeqs[seqCounter])
		{
			GetString (seqName,unal,seqCounter);
			fprintf (fName,">",seqName,"\n",UnalignedSeqs[seqCounter],"\n");
		}
	}
}

function checkFramePosition (pos1, pos2, fr1, fr2)
{
	fSpan  = pos2-pos1;
	
	if (fSpan>1) /* first followed by second*/
	{
		if (fSpan < fMin)
		{
			fMin = fSpan;
			frame1 = fr1;
			frame2 = fr2;
			fStart = pos1+1;
			fEnd   = pos2;
		}
	}	
	return 0;
}

