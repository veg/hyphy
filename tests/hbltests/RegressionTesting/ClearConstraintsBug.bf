#NEXUS

BEGIN TAXA;
	DIMENSIONS NTAX = 11;
	TAXLABELS
		'D1_5' 'D1_10' 'D1_8' 'D1_9' 'D1_7' 'D1_0' 'D1_1' 'D1_2' 'D1_3' 'D1_4' 'D1_6' ;
END;

BEGIN CHARACTERS;
	DIMENSIONS NCHAR = 882;
	FORMAT
		DATATYPE = DNA
		GAP=-
		MISSING=?
	;

MATRIX
	'D1_5'   TGCATTGATGATGTGGTTTGCACTGATAATACTACTAGATCCAATTGTACTTTTTTGACTAATGATAATAGTACCACTAATCTTACTTGGGGAAAAATGGAAAAAGGAGAAATAAAAAATTGCTCTTTCAATGTCACCAGTATAACAAATAAGATGCAGAAAGAATATGCACTTTTTTATAAACTTGATGTAATGCCAATAAATAATGACAGTTATACACTGATAAATTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAACCCATTCCTATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAGTGTAATAATAAGACATTCAATGGATCAGGACCATGTACAAATGTCAGCACAGTACAATGTACACATGGAATTAGGCCAGTAGTGTCAACTCAACTACTGTTAAATGACAGTCTAGCAGAAGAGGAGGTAATGATTAGGTCTGAGAATTTCTCGGACAATGCTAAAATCATAATAGTACAGCTAAATAAAACTGTAGAAATTAACTGTACAAGACCCAATAACAATACAAGAAGAAGTATACATATAGCACCAGGGAGAGCATTCTATGCAACAGATGGTATAATAGGAGATATAAGAAAAGCACATTGTAACATTAGTGCAGCAAAATGGAATAACACTCTAAAACAGATAGTTACAAAATTAAGGGAACAATATGGCAATAAAATAATAGTCTTTAATCAATCCTCAGGAGGGGACCCAGAAATTGTAATGCACAGTTTTAATTGTGGAGGAGAATTTTTTTACTGTAATACATCACAGCTGTTTAATAGTACTTGGCTGTCTAATGGTACTTCTAATCTTGAGGGGCCAGGAAATAGCACAATCACACTTCCATGC
	'D1_10'  TGCATTGATGATGTGGTTTGCACTGATAATACTACTAGATCCAATTGTACTTTTTTGACTAATGATAATAGTACCACTAATCTTACTTGGGGAAAAATGGAAAAAGGAGAAATAAAAAATTGCTCTTTCAATGTCACCAGTATAACAAATAAGATGCAGAAAGAATATGCACTTTTTTATAAACTTGATGTAATGCCAATAAATAATGACAGTTATACACTGATAAATTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAACCCATTCCTATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAGTGTAATAATAAGACATTCAATGGATCAGGACCATGTACAAATGTCAGCACAGTACAATGTACACATGGAATTAGGCCAGTAGTGTCAACTCAACTACTGTTAAATGGCAGTCTAGCAGAAGAGGAGGTAATGATTAGGTCTGAGAATTTCTCGGACAATGCTAAAATCATAATAGTACAGCTAAATAAAACTGTAGAAATTAACTGTACAAGACCCAATAACAATACAAGAAGAAGTATACATATAGCACCAGGGAGAGCATTCTATGCAACAGATGGTATAATAGGAGATATAAGAAAAGCACATTGTAACATTAGTGCAGCAAAATGGAATAACACTCTAAAACAGATAGTTACAAAATTAAGGGAACAATATGGCAATAAAATAATAGTCTTTAATCAATCCTCAGGAGGGGACCCAGAAATTGTAATGCACAGTTTTAATTGTGGAGGAGAATTTTTTTACTGTAATACATCACAGCTGTTTAATAGTACTTGGCTGTCTAATGGTACTTCTAATCTTGAGGGGCCAGGAAATAGCACAATCACACTTCCATGC
	'D1_8'   TGCATTGATGATGTGGTTTGCACTGATAATACTACTAGATCCAATTGTACTTTTTTGACTAATGATAATAGTACCACTAATCTTACTTGGGGAAAAATGGAAAAAGGAGAAATAAAAAATTGCTCTTTCAATGTCACCAGTATAACAAATAAGATGCAGAAAGAATATGCACTTTTTTATAAACTTGATGTAATGCCAATAAATAATGACAGTTATACACTGATAAATTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAACCCATTCCTATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAGTGTAATAATAAGACATTCAATGGATCAGGACCATGTACAAATGTCAGCACAGTACAATGTACACATGGAATTAGGCCAGTAGTGTCAACTCAACTACTGTTAAATGGCAGTCTAGCAGAAGAGGAGGTAATGATTAGGTCTGAGAATTTCTCGGACAATGCTAAAATCATAATAGTACGGCTAAATAAAACTGTAGAAATTAACTGTACAAGACCCAATAACAATACAAGAAGAAGTATACATATAGCACCAGGGAGAGCATTCTATGCAACAGATGGTATAATAGGAGATATAAGAAAAGCACATTGTAACATTAGTGCAGCAAAATGGAATAACACTCTAAAACAGATAGTTACAAAATTAAGGGAACAATATGGCAATAAAATAATAGTCTTTAATCAATCCTCAGGAGGGGACCCAGAAATTGTAATGCACAGTTTTAATTGTGGAGGAGAATTTTTTTACTGTAATACATCACAGCTGTTTAATAGTACTTGGCTGTCTAATGGTACTTCTAATCTTGAGGGGCCAGGAAATAGCACAATCACACTTCCATGC
	'D1_9'   TGCATTGATGATGTGGTTTGCACTGATAATACTACTAGATCCAATTGTACTTTTTTGACTAATGATAATAGTACCACTAATCTTACTTGGGGAAAAATGGAAAAAGGAGAAATAAAAAATTGCTCTTTCAATGTCACCAGTATAACAAATAAGATGCAGAAAGAATATGCACTTTTTTATAAACTTGATGTAATGCCAATAAATAATGACAGTTATACACTGATAAATTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAACCCATTCCTATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAGTGTAATAATAAGACATTCAATGGATCAGGACCATGTACAAATGTCAGCACAGTACAATGTACACATGGAATTAGGCCAGTAGTGTCAACTCAACTACTGTTAAATGGCAGTCTAGCAGAAGAGGAGGTAATGATTAGGTCTGAGAATTTCTCGGACAATGCTAAAATCATAATAGTACAGCTAAATAAAACTGTAGAAATTAACTGTACAAGACCCAATAACAATACAAGAAGAAGTATACATATAGCACCAGGGAGAGCATTCTATGCAACAGATGGTATAATAGGAGATATAAGAAAAGCACATTGTAACATTAGTGCAGCAAAATGGAATAACACTCTAAAACAGATAGTTACAAAATTAAGGGAACAATATGGCAATAAAATAATAGTCTTTAATCAATCCTCAGGAGGGGACCCAGAAATTGTAATGCACAGTTTTAATTGTGGAGGAGAATTTTTTTACTGTAATACATCACAGCTGTTTAATAGTACTTGGCTGTCTAATGGTACTTCTAATCTTGAGGGGCCAGGAAATAGCACAATCACACTTCCATGC
	'D1_7'   TGCATTGATGATGTGGTTTGCACTGATAATACTACTAGATCCAATTGTACTTTTTTGACTAATGATAATAGTACCACTAATCTTACTTGGGGAAAAATGGAAAAAGGAGAAATAAAAAATTGCTCTTTCAATGTCACCAGTATAACAAATAAGATGCAGAAAGAATATGCACTTTTTTATAAACTTGATGTAATGCCAATAAATAATGACAGTTATACACTGATAAATTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAACCCATTCCTATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAGTGTAATAATAAGACATTCAATGGATCAGGACCATGTACAAATGTCAGCACAGTACAATGTACACATGGAATTAGGCCAGTAGTGTCAACTCAACTACTGTTAAATGGCAGTCTAGCAGAAGAGGAGGTAATGATTAGGTCTGAGAATTTCTCGGACAATGCTAAAATCATAATAGTACAGCTAAATAAAACTGTAGAAATTAACTGTACAAGACCCAATAACAATACAAGAAGAAGTATACATATAGCACCAGGGAGAGCATTCTATGCAACAGATGGTATAATAGGAGATATAAGAAAAGCACATTGTAACATTAGTGCAGCAAAATGGAATAACACTCTAAAACAGATAGTTACAAAATTAAGGGAACAATATGGCAATAAAATAATAGTCTTTAATCAATCCTCAGGAGGGGACCCAGAAATTGTAATGCACAGTTTTAATTGTGGAGGAGAATTTTTTTACTGTAATACATCACAGCTGTTTAATAGTACTTGGCTGTCTAATGGTACTTCTAATCTTGAGGGGCCAGGAAATAGCACAATCACACTTCCATGC
	'D1_0'   TGCATTGATGATGTGGTTTGCACTGATAATACTACTAGATCCAATTGTACTTTTTTGACTAATGATAATAGTACCACTAATCTTACTTGGGGAAAAATGGAAAAAGGAGAAATAAAAAATTGCTCTTTCAATGTCACCAGTATAACAAATAAGATGCAGAAAGAATATGCACTTTTTTATAAACTTGATGTAATGCCAATAAATAATGACAGTTATACACTGATAAATTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAACCCATTCCTATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAGTGTAATAATAAGACATTCAATGGATCAGGACCATGTACAAATGTCAGCACAGTACAATGTACACATGGAATTAGGCCAGTAGTGTCAACTCAACTACTGTTAAATGGCAGTCTAGCAGAAGAGGAGGTAATGATTAGGTCTGAGAATTTCTCGGACAATGCTAAAATCATAATAGTACAGCTAAATAAAACTGTAGAAATTAACTGTACAAGACCCAATAACAATACAAGAAGAAGTATACATATAGCACCAGGGAGAGCATTCTATGCAACAGATGGTATAATAGGAGATATAAGAAAAGCACATTGTAACATTAGTGCAGCAAAATGGAATAACACTCTAAAACAGATAGTTACAAAATTAAGGGAACAATATGGCAATAAAATAATAGTCTTTAATCAATCCTCAGGAGGGGACCCAGAAATTGTAATGCACAGTTTTAATTGTGGAGGAGAATTTTTTTACTGTAATACATCACAGCTGTTTAATAGTACTTGGCTGTCTAATGGTACTTCTAATCTTGAGGGGCCAGGAAATAGCACAATCACACTTCCATGC
	'D1_1'   TGCATTGATGATGTGGTTTGCACTGATAATACTACTAGATCCAATTGTACTTTTTTGACTAATGATAATAGTACCACTAATCTTACTTGGGGAAAAATGGAAAAAGGAGAAATAAAAAATTGCTCTTTCAATGTCACCAGTATAACAAATAAGATGCAGAAAGAATATGCACTTTTTTATAAACTTGATGTAATGCCAATAAATAATGACAGTTATACACTGATAAATTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAACCCATTCCCATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAGTGTAATAATAAGACATTCAATGGATCAGGACCATGTACAAATGTCAGTACAGTACAATGTACACATGGAATCAGGCCAGTAGTGTCAACTCAACTACTGTTAAATGGCAGTCTAGCAGAAGAGGAGGTAATGATTAGGTCTGAGAATTTCTCGAGAAATGATAGAATCATAATAGTACAGCTGAATGAAACTGTAGAAATTAATTGTACAAGACCCAATAACAATACAAGAAAAAGTATACATATAGCACCAGGGAGGGCATTCTATGCAACAGATGGTATAATAGGAGATATAAGAAAAGCACATTGTAACATTAGTGAAACAAAATGGAGGAACACTTTAAAACAGATAGCTACAAAATTAAGAGAACAATATGAGAATGAAAAAATAGTCTTTAATCAAACCTCAGGAGGGGACCCAGAAATTGTAATGCACAGTTTTAATTGTGGAGGAGAATTTTTTTACTGTAATACATCACAGCTGTTTAATAGTACTTGGCTGTCTAATGGTACTTCTAATCTTGAGGGGCCAGGAAATAGCACAATCACACTCCCATGC
	'D1_2'   TGCATTGATGATGTGGTTTGCACTGATAATACTACTAGATCCAATTGTACTTTTTTGACTAATGATAATAGTACCACTAATCTTACTTGGGGAAAAATGGAAAAAGGAGAAATAAAAAATTGCTCTTTCAATGTCACCAGTATAACAAATAAGATGCAGAAAGAATATGCACTTTTTTATAAACTTGATGTAATGCCGATAAATAATGACAGTTATACACTGATAAATTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAACCCATTCCCATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAGTGTAATAATAAGACATTCAATGGATCAGGACCATGTACAAATGTCAGTACAGTACAATGTACACATGGAATCAGGCCAGTAGTGTCAACTCAACTACTGTTAAATGGCAGTCTAGCAGAAGAGGAGGTAATGATTAGGTCTGAGAATTTCTCGAGAAATGATAGAATCATAATAGTACAGCTGAATGAAACTGTAGAAATTAATTGTACAAGACCCAATAACAATACAAGAAAAAGTATACATATAGCACCAGGGAGGGCATTCTATGCAACAGATGGTATAATAGGAGATATAAGAAAAGCACATTGTAACATTAGTGAAACAAAATGGAGGAACACTTTAAAACAGATAGCTACAAAATTAAGAGAACAATATGAGAATGAAAAAATAGTCTTTAATCAAACCTCAGGAGGGGACCCAGAAATTGTAATGCACAGTTTTAATTGTGGAGGAGAATTTTTTTACTGTAATACATCACAGCTGTTTAATAGTACTTGGCTGTCTAATGGTACTTCTAATCTTGAGGGGCCAGGAAATAGCACAATCACACTCCCATGC
	'D1_3'   TGCATTGATGATGTGGTTTGCACTGATAATACTACTAGATCCAATTGTACTTTTTTGACTAATGATAATAGTACCACTAATCTTACTTGGGGAAAAATGGAAAAAGGAGAAATAAAAAATTGCTCTTTCAATGTCACCAGTATAACAAATAAGATGCAGAAAGAATCTGCACTTTTTTATAAACTTGATGTAATGCCAATAAATAATGACAGTTATACACTGATAAATTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAACCCATTCCCATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAGTGTAATAATAAGACATTCAATGGATCAGGACCATGTACAAATGTCAGTACAGTACAATGTACACATGGAATCAGGCCAGTAGTGTCAACTCAACTACTGTTAAATGGCAGTCTAGCAGAAGAGGAGGTAATGATTAGGTCTGAGAATTTCTCGAGAAATGATAGAATCATAATAGTACAGCTGAATGAAACTGTAGAAATTAATTGTACAAGACCCAATAACAATACAAGAAAAAGTATACATATAGCACCAGGGAGGGCATTCTATGCAACAGATGGTATAATAGGAGATATAAGAAAAGCACATTGTAACATTAGTGAAACAAAATGGAGGAACACTTTAAAACAGATAGCTACAAAATTAAGAGAACAATATGAGAATGAAAAAATAGTCTTTAATCAAACCTCAGGAGGGGACCCAGAAATTGTAATGCACAGTTTTAATTGTGGAGGAGAATTTTTTTACTGTAATACATCACAGCTGTTTAATAGTACTTGGCTGTCTAATGGTACTTCTAATCTTGAGGGGCCAGGAAATAGCACAATCACACTCCCATGC
	'D1_4'   TGCATTGATGATGTGGTTTGCACTGATAATACTACTAGATCCAATTGTACTTTTTTGACTAATGATAATAGTACCACTAATCTTACTTGGGGAAAAATGGAAAAAGGAGAAATAAAAAATTGCTCTTTCAATGTCACCAGTATAACAAATAAGATGCAGAAAGAATATGCACTTTTTTATAAACTTGATGTAATGCCAATAAATAATGACAGTTATACACTGATAAATTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAACCCATTCCCATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAGTGTAATAATAAGACATTCAATGGATCAGGACCATGTACAAATGTCAGTACAGTACAATGTACACATGGAATCAGGCCAGTAGTGTCAACTCAACTACTGTTAAATGGCAGTCTAGCAGAAGAGGAGGTAATGATTAGGTCTGAGAATTTCTCGAGAAATGATAGAATCATAATAGTACAGCTGAATGAAACTGTAGAAATTAATTGTACAAGACCCAATAACAATACAAGAAAAAGTATACATATAGCACCAGGGAGGGCATTCTATGCAACAGATGGTATAATAGGAGATATAAGAAAAGCACATTGTAACATTAGTGAAACAAAATGGAGGAACACTTTAAAACAGATAGCTACAAAATTAAGAGAACAATATGAGAATGAAAAAATAGTCTTTAATCAAACCTCAGGAGGGGACCCAGAAATTGTAATGCACAGTTTTAATTGTGGAGGAGAATTTTTTTACTGTAATACATCACAGCTGTTTAATAGTACTTGGCTGTCTAATGGTACTTCTAATCTTGAGGGGCCAGGAAATAGCACAATCACACTCCCATGC
	'D1_6'   TGCATTGATGATGTGGTTTGCACTGATAATACTACTAGATCCAATTGTACTTTTTTGACTAATGATAATAGTACCACTAATCTTACTTGGGGAAAAATGGAAAAAGGAGAAATAAAAAATTGCTCTTTCAATGTCACCAGTATAACAAATAAGATGCAGAAAGAATATGCACTTTTTTATAAACTTGATGTAATGCCAATAAATAATGACAGTTATACACTGATAAATTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAACCCATTCCTATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAGTGTAATAATAAGACATTCAATGGATCAGGACCATGTACAAATGTCAGCACAGTACAATGTACACATGGAATTAGGCCAGTAGTGTCAACTCAACTACTGTTAAATGGCAGTCTAGCAGAAGAGGAGGTAATGATTAGGTCTGAGAATTTCTCGGACAATGCTAAAATCATAATAGTACAGCTAAATAAAACTGTAGAAATTAACTGTACAAGACCCAATAACAATACAAGAAGAAGTATACATATAGCACCAGGGAGAGCATTCTATGCAACAGATGGTATAATAGGAGATATAAGAAAAGCACATTGTAACATTAGTGCAGCAAAATGGAATAACACTCTAAAACAGATAGTTACAAAATTAAGGGAACAATATGGCAATAAAATAATAGTCTTTAATCAATCCTCAGGAGGGGACCCAGAAATTGTAATGCACAGTTTTAATTGTGGAGGAGAATTTTTTTACTGTAATACATCACAGCTGTTTAATAGTACTTGGCTGTCTAATGGTACTTCTAATCTTGAGGGGCCAGGAAATAGCACAATCACACTTCCATGC;
END;

BEGIN TREES;
	TREE tree = (D1_5:0.00113652,D1_10:0,D1_8:0.00113575,D1_9:0,D1_7:0,D1_0:0,(D1_1:0,D1_2:0.0011321,D1_3:0.00113262,D1_4:0):0.030296,D1_6:0);
END;

BEGIN HYPHY;


global omega=0.4906081221454927;
global _fasta_part_Shared_GT=0.10978373885463;
global _fasta_part_Shared_AC=0.3275863477719206;
global _fasta_part_Shared_AT=0.09055098248123046;
global _fasta_part_Shared_CG=0.1590071141754616;
global _fasta_part_Shared_CT=0.6330418473697754;
_fasta_part_MG94xREV_3x4={61,61};
_fasta_part_MG94xREV_3x4[0][1]:=_fasta_part_Shared_AC*nonSynRate*0.128015;
_fasta_part_MG94xREV_3x4[0][2]:=synRate*0.12987;
_fasta_part_MG94xREV_3x4[0][3]:=_fasta_part_Shared_AT*nonSynRate*0.38188;
_fasta_part_MG94xREV_3x4[0][4]:=_fasta_part_Shared_AC*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[0][8]:=nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[0][12]:=_fasta_part_Shared_AT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[0][16]:=_fasta_part_Shared_AC*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[0][32]:=nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[1][0]:=_fasta_part_Shared_AC*nonSynRate*0.360235;
_fasta_part_MG94xREV_3x4[1][2]:=_fasta_part_Shared_CG*nonSynRate*0.12987;
_fasta_part_MG94xREV_3x4[1][3]:=_fasta_part_Shared_CT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[1][5]:=_fasta_part_Shared_AC*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[1][9]:=nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[1][13]:=_fasta_part_Shared_AT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[1][17]:=_fasta_part_Shared_AC*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[1][33]:=nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[1][48]:=_fasta_part_Shared_AT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[2][0]:=synRate*0.360235;
_fasta_part_MG94xREV_3x4[2][1]:=_fasta_part_Shared_CG*nonSynRate*0.128015;
_fasta_part_MG94xREV_3x4[2][3]:=_fasta_part_Shared_GT*nonSynRate*0.38188;
_fasta_part_MG94xREV_3x4[2][6]:=_fasta_part_Shared_AC*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[2][10]:=nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[2][14]:=_fasta_part_Shared_AT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[2][18]:=_fasta_part_Shared_AC*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[2][34]:=nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[3][0]:=_fasta_part_Shared_AT*nonSynRate*0.360235;
_fasta_part_MG94xREV_3x4[3][1]:=_fasta_part_Shared_CT*synRate*0.128015;
_fasta_part_MG94xREV_3x4[3][2]:=_fasta_part_Shared_GT*nonSynRate*0.12987;
_fasta_part_MG94xREV_3x4[3][7]:=_fasta_part_Shared_AC*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[3][11]:=nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[3][15]:=_fasta_part_Shared_AT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[3][19]:=_fasta_part_Shared_AC*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[3][35]:=nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[3][49]:=_fasta_part_Shared_AT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[4][0]:=_fasta_part_Shared_AC*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[4][5]:=_fasta_part_Shared_AC*synRate*0.128015;
_fasta_part_MG94xREV_3x4[4][6]:=synRate*0.12987;
_fasta_part_MG94xREV_3x4[4][7]:=_fasta_part_Shared_AT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[4][8]:=_fasta_part_Shared_CG*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[4][12]:=_fasta_part_Shared_CT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[4][20]:=_fasta_part_Shared_AC*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[4][36]:=nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[4][50]:=_fasta_part_Shared_AT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[5][1]:=_fasta_part_Shared_AC*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[5][4]:=_fasta_part_Shared_AC*synRate*0.360235;
_fasta_part_MG94xREV_3x4[5][6]:=_fasta_part_Shared_CG*synRate*0.12987;
_fasta_part_MG94xREV_3x4[5][7]:=_fasta_part_Shared_CT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[5][9]:=_fasta_part_Shared_CG*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[5][13]:=_fasta_part_Shared_CT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[5][21]:=_fasta_part_Shared_AC*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[5][37]:=nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[5][51]:=_fasta_part_Shared_AT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[6][2]:=_fasta_part_Shared_AC*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[6][4]:=synRate*0.360235;
_fasta_part_MG94xREV_3x4[6][5]:=_fasta_part_Shared_CG*synRate*0.128015;
_fasta_part_MG94xREV_3x4[6][7]:=_fasta_part_Shared_GT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[6][10]:=_fasta_part_Shared_CG*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[6][14]:=_fasta_part_Shared_CT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[6][22]:=_fasta_part_Shared_AC*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[6][38]:=nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[6][52]:=_fasta_part_Shared_AT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[7][3]:=_fasta_part_Shared_AC*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[7][4]:=_fasta_part_Shared_AT*synRate*0.360235;
_fasta_part_MG94xREV_3x4[7][5]:=_fasta_part_Shared_CT*synRate*0.128015;
_fasta_part_MG94xREV_3x4[7][6]:=_fasta_part_Shared_GT*synRate*0.12987;
_fasta_part_MG94xREV_3x4[7][11]:=_fasta_part_Shared_CG*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[7][15]:=_fasta_part_Shared_CT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[7][23]:=_fasta_part_Shared_AC*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[7][39]:=nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[7][53]:=_fasta_part_Shared_AT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[8][0]:=nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[8][4]:=_fasta_part_Shared_CG*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[8][9]:=_fasta_part_Shared_AC*nonSynRate*0.128015;
_fasta_part_MG94xREV_3x4[8][10]:=synRate*0.12987;
_fasta_part_MG94xREV_3x4[8][11]:=_fasta_part_Shared_AT*nonSynRate*0.38188;
_fasta_part_MG94xREV_3x4[8][12]:=_fasta_part_Shared_GT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[8][24]:=_fasta_part_Shared_AC*synRate*0.134818;
_fasta_part_MG94xREV_3x4[8][40]:=nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[9][1]:=nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[9][5]:=_fasta_part_Shared_CG*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[9][8]:=_fasta_part_Shared_AC*nonSynRate*0.360235;
_fasta_part_MG94xREV_3x4[9][10]:=_fasta_part_Shared_CG*nonSynRate*0.12987;
_fasta_part_MG94xREV_3x4[9][11]:=_fasta_part_Shared_CT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[9][13]:=_fasta_part_Shared_GT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[9][25]:=_fasta_part_Shared_AC*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[9][41]:=nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[9][54]:=_fasta_part_Shared_AT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[10][2]:=nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[10][6]:=_fasta_part_Shared_CG*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[10][8]:=synRate*0.360235;
_fasta_part_MG94xREV_3x4[10][9]:=_fasta_part_Shared_CG*nonSynRate*0.128015;
_fasta_part_MG94xREV_3x4[10][11]:=_fasta_part_Shared_GT*nonSynRate*0.38188;
_fasta_part_MG94xREV_3x4[10][14]:=_fasta_part_Shared_GT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[10][26]:=_fasta_part_Shared_AC*synRate*0.134818;
_fasta_part_MG94xREV_3x4[10][42]:=nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[10][55]:=_fasta_part_Shared_AT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[11][3]:=nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[11][7]:=_fasta_part_Shared_CG*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[11][8]:=_fasta_part_Shared_AT*nonSynRate*0.360235;
_fasta_part_MG94xREV_3x4[11][9]:=_fasta_part_Shared_CT*synRate*0.128015;
_fasta_part_MG94xREV_3x4[11][10]:=_fasta_part_Shared_GT*nonSynRate*0.12987;
_fasta_part_MG94xREV_3x4[11][15]:=_fasta_part_Shared_GT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[11][27]:=_fasta_part_Shared_AC*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[11][43]:=nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[11][56]:=_fasta_part_Shared_AT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[12][0]:=_fasta_part_Shared_AT*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[12][4]:=_fasta_part_Shared_CT*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[12][8]:=_fasta_part_Shared_GT*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[12][13]:=_fasta_part_Shared_AC*synRate*0.128015;
_fasta_part_MG94xREV_3x4[12][14]:=nonSynRate*0.12987;
_fasta_part_MG94xREV_3x4[12][15]:=_fasta_part_Shared_AT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[12][28]:=_fasta_part_Shared_AC*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[12][44]:=nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[12][57]:=_fasta_part_Shared_AT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[13][1]:=_fasta_part_Shared_AT*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[13][5]:=_fasta_part_Shared_CT*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[13][9]:=_fasta_part_Shared_GT*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[13][12]:=_fasta_part_Shared_AC*synRate*0.360235;
_fasta_part_MG94xREV_3x4[13][14]:=_fasta_part_Shared_CG*nonSynRate*0.12987;
_fasta_part_MG94xREV_3x4[13][15]:=_fasta_part_Shared_CT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[13][29]:=_fasta_part_Shared_AC*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[13][45]:=nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[13][58]:=_fasta_part_Shared_AT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[14][2]:=_fasta_part_Shared_AT*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[14][6]:=_fasta_part_Shared_CT*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[14][10]:=_fasta_part_Shared_GT*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[14][12]:=nonSynRate*0.360235;
_fasta_part_MG94xREV_3x4[14][13]:=_fasta_part_Shared_CG*nonSynRate*0.128015;
_fasta_part_MG94xREV_3x4[14][15]:=_fasta_part_Shared_GT*nonSynRate*0.38188;
_fasta_part_MG94xREV_3x4[14][30]:=_fasta_part_Shared_AC*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[14][46]:=nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[14][59]:=_fasta_part_Shared_AT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[15][3]:=_fasta_part_Shared_AT*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[15][7]:=_fasta_part_Shared_CT*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[15][11]:=_fasta_part_Shared_GT*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[15][12]:=_fasta_part_Shared_AT*synRate*0.360235;
_fasta_part_MG94xREV_3x4[15][13]:=_fasta_part_Shared_CT*synRate*0.128015;
_fasta_part_MG94xREV_3x4[15][14]:=_fasta_part_Shared_GT*nonSynRate*0.12987;
_fasta_part_MG94xREV_3x4[15][31]:=_fasta_part_Shared_AC*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[15][47]:=nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[15][60]:=_fasta_part_Shared_AT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[16][0]:=_fasta_part_Shared_AC*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[16][17]:=_fasta_part_Shared_AC*nonSynRate*0.128015;
_fasta_part_MG94xREV_3x4[16][18]:=synRate*0.12987;
_fasta_part_MG94xREV_3x4[16][19]:=_fasta_part_Shared_AT*nonSynRate*0.38188;
_fasta_part_MG94xREV_3x4[16][20]:=_fasta_part_Shared_AC*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[16][24]:=nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[16][28]:=_fasta_part_Shared_AT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[16][32]:=_fasta_part_Shared_CG*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[17][1]:=_fasta_part_Shared_AC*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[17][16]:=_fasta_part_Shared_AC*nonSynRate*0.360235;
_fasta_part_MG94xREV_3x4[17][18]:=_fasta_part_Shared_CG*nonSynRate*0.12987;
_fasta_part_MG94xREV_3x4[17][19]:=_fasta_part_Shared_CT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[17][21]:=_fasta_part_Shared_AC*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[17][25]:=nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[17][29]:=_fasta_part_Shared_AT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[17][33]:=_fasta_part_Shared_CG*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[17][48]:=_fasta_part_Shared_CT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[18][2]:=_fasta_part_Shared_AC*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[18][16]:=synRate*0.360235;
_fasta_part_MG94xREV_3x4[18][17]:=_fasta_part_Shared_CG*nonSynRate*0.128015;
_fasta_part_MG94xREV_3x4[18][19]:=_fasta_part_Shared_GT*nonSynRate*0.38188;
_fasta_part_MG94xREV_3x4[18][22]:=_fasta_part_Shared_AC*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[18][26]:=nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[18][30]:=_fasta_part_Shared_AT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[18][34]:=_fasta_part_Shared_CG*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[19][3]:=_fasta_part_Shared_AC*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[19][16]:=_fasta_part_Shared_AT*nonSynRate*0.360235;
_fasta_part_MG94xREV_3x4[19][17]:=_fasta_part_Shared_CT*synRate*0.128015;
_fasta_part_MG94xREV_3x4[19][18]:=_fasta_part_Shared_GT*nonSynRate*0.12987;
_fasta_part_MG94xREV_3x4[19][23]:=_fasta_part_Shared_AC*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[19][27]:=nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[19][31]:=_fasta_part_Shared_AT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[19][35]:=_fasta_part_Shared_CG*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[19][49]:=_fasta_part_Shared_CT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[20][4]:=_fasta_part_Shared_AC*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[20][16]:=_fasta_part_Shared_AC*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[20][21]:=_fasta_part_Shared_AC*synRate*0.128015;
_fasta_part_MG94xREV_3x4[20][22]:=synRate*0.12987;
_fasta_part_MG94xREV_3x4[20][23]:=_fasta_part_Shared_AT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[20][24]:=_fasta_part_Shared_CG*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[20][28]:=_fasta_part_Shared_CT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[20][36]:=_fasta_part_Shared_CG*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[20][50]:=_fasta_part_Shared_CT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[21][5]:=_fasta_part_Shared_AC*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[21][17]:=_fasta_part_Shared_AC*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[21][20]:=_fasta_part_Shared_AC*synRate*0.360235;
_fasta_part_MG94xREV_3x4[21][22]:=_fasta_part_Shared_CG*synRate*0.12987;
_fasta_part_MG94xREV_3x4[21][23]:=_fasta_part_Shared_CT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[21][25]:=_fasta_part_Shared_CG*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[21][29]:=_fasta_part_Shared_CT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[21][37]:=_fasta_part_Shared_CG*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[21][51]:=_fasta_part_Shared_CT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[22][6]:=_fasta_part_Shared_AC*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[22][18]:=_fasta_part_Shared_AC*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[22][20]:=synRate*0.360235;
_fasta_part_MG94xREV_3x4[22][21]:=_fasta_part_Shared_CG*synRate*0.128015;
_fasta_part_MG94xREV_3x4[22][23]:=_fasta_part_Shared_GT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[22][26]:=_fasta_part_Shared_CG*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[22][30]:=_fasta_part_Shared_CT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[22][38]:=_fasta_part_Shared_CG*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[22][52]:=_fasta_part_Shared_CT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[23][7]:=_fasta_part_Shared_AC*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[23][19]:=_fasta_part_Shared_AC*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[23][20]:=_fasta_part_Shared_AT*synRate*0.360235;
_fasta_part_MG94xREV_3x4[23][21]:=_fasta_part_Shared_CT*synRate*0.128015;
_fasta_part_MG94xREV_3x4[23][22]:=_fasta_part_Shared_GT*synRate*0.12987;
_fasta_part_MG94xREV_3x4[23][27]:=_fasta_part_Shared_CG*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[23][31]:=_fasta_part_Shared_CT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[23][39]:=_fasta_part_Shared_CG*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[23][53]:=_fasta_part_Shared_CT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[24][8]:=_fasta_part_Shared_AC*synRate*0.443414;
_fasta_part_MG94xREV_3x4[24][16]:=nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[24][20]:=_fasta_part_Shared_CG*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[24][25]:=_fasta_part_Shared_AC*synRate*0.128015;
_fasta_part_MG94xREV_3x4[24][26]:=synRate*0.12987;
_fasta_part_MG94xREV_3x4[24][27]:=_fasta_part_Shared_AT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[24][28]:=_fasta_part_Shared_GT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[24][40]:=_fasta_part_Shared_CG*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[25][9]:=_fasta_part_Shared_AC*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[25][17]:=nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[25][21]:=_fasta_part_Shared_CG*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[25][24]:=_fasta_part_Shared_AC*synRate*0.360235;
_fasta_part_MG94xREV_3x4[25][26]:=_fasta_part_Shared_CG*synRate*0.12987;
_fasta_part_MG94xREV_3x4[25][27]:=_fasta_part_Shared_CT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[25][29]:=_fasta_part_Shared_GT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[25][41]:=_fasta_part_Shared_CG*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[25][54]:=_fasta_part_Shared_CT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[26][10]:=_fasta_part_Shared_AC*synRate*0.443414;
_fasta_part_MG94xREV_3x4[26][18]:=nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[26][22]:=_fasta_part_Shared_CG*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[26][24]:=synRate*0.360235;
_fasta_part_MG94xREV_3x4[26][25]:=_fasta_part_Shared_CG*synRate*0.128015;
_fasta_part_MG94xREV_3x4[26][27]:=_fasta_part_Shared_GT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[26][30]:=_fasta_part_Shared_GT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[26][42]:=_fasta_part_Shared_CG*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[26][55]:=_fasta_part_Shared_CT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[27][11]:=_fasta_part_Shared_AC*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[27][19]:=nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[27][23]:=_fasta_part_Shared_CG*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[27][24]:=_fasta_part_Shared_AT*synRate*0.360235;
_fasta_part_MG94xREV_3x4[27][25]:=_fasta_part_Shared_CT*synRate*0.128015;
_fasta_part_MG94xREV_3x4[27][26]:=_fasta_part_Shared_GT*synRate*0.12987;
_fasta_part_MG94xREV_3x4[27][31]:=_fasta_part_Shared_GT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[27][43]:=_fasta_part_Shared_CG*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[27][56]:=_fasta_part_Shared_CT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[28][12]:=_fasta_part_Shared_AC*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[28][16]:=_fasta_part_Shared_AT*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[28][20]:=_fasta_part_Shared_CT*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[28][24]:=_fasta_part_Shared_GT*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[28][29]:=_fasta_part_Shared_AC*synRate*0.128015;
_fasta_part_MG94xREV_3x4[28][30]:=synRate*0.12987;
_fasta_part_MG94xREV_3x4[28][31]:=_fasta_part_Shared_AT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[28][44]:=_fasta_part_Shared_CG*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[28][57]:=_fasta_part_Shared_CT*synRate*0.183673;
_fasta_part_MG94xREV_3x4[29][13]:=_fasta_part_Shared_AC*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[29][17]:=_fasta_part_Shared_AT*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[29][21]:=_fasta_part_Shared_CT*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[29][25]:=_fasta_part_Shared_GT*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[29][28]:=_fasta_part_Shared_AC*synRate*0.360235;
_fasta_part_MG94xREV_3x4[29][30]:=_fasta_part_Shared_CG*synRate*0.12987;
_fasta_part_MG94xREV_3x4[29][31]:=_fasta_part_Shared_CT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[29][45]:=_fasta_part_Shared_CG*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[29][58]:=_fasta_part_Shared_CT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[30][14]:=_fasta_part_Shared_AC*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[30][18]:=_fasta_part_Shared_AT*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[30][22]:=_fasta_part_Shared_CT*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[30][26]:=_fasta_part_Shared_GT*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[30][28]:=synRate*0.360235;
_fasta_part_MG94xREV_3x4[30][29]:=_fasta_part_Shared_CG*synRate*0.128015;
_fasta_part_MG94xREV_3x4[30][31]:=_fasta_part_Shared_GT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[30][46]:=_fasta_part_Shared_CG*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[30][59]:=_fasta_part_Shared_CT*synRate*0.183673;
_fasta_part_MG94xREV_3x4[31][15]:=_fasta_part_Shared_AC*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[31][19]:=_fasta_part_Shared_AT*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[31][23]:=_fasta_part_Shared_CT*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[31][27]:=_fasta_part_Shared_GT*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[31][28]:=_fasta_part_Shared_AT*synRate*0.360235;
_fasta_part_MG94xREV_3x4[31][29]:=_fasta_part_Shared_CT*synRate*0.128015;
_fasta_part_MG94xREV_3x4[31][30]:=_fasta_part_Shared_GT*synRate*0.12987;
_fasta_part_MG94xREV_3x4[31][47]:=_fasta_part_Shared_CG*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[31][60]:=_fasta_part_Shared_CT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[32][0]:=nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[32][16]:=_fasta_part_Shared_CG*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[32][33]:=_fasta_part_Shared_AC*nonSynRate*0.128015;
_fasta_part_MG94xREV_3x4[32][34]:=synRate*0.12987;
_fasta_part_MG94xREV_3x4[32][35]:=_fasta_part_Shared_AT*nonSynRate*0.38188;
_fasta_part_MG94xREV_3x4[32][36]:=_fasta_part_Shared_AC*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[32][40]:=nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[32][44]:=_fasta_part_Shared_AT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[33][1]:=nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[33][17]:=_fasta_part_Shared_CG*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[33][32]:=_fasta_part_Shared_AC*nonSynRate*0.360235;
_fasta_part_MG94xREV_3x4[33][34]:=_fasta_part_Shared_CG*nonSynRate*0.12987;
_fasta_part_MG94xREV_3x4[33][35]:=_fasta_part_Shared_CT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[33][37]:=_fasta_part_Shared_AC*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[33][41]:=nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[33][45]:=_fasta_part_Shared_AT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[33][48]:=_fasta_part_Shared_GT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[34][2]:=nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[34][18]:=_fasta_part_Shared_CG*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[34][32]:=synRate*0.360235;
_fasta_part_MG94xREV_3x4[34][33]:=_fasta_part_Shared_CG*nonSynRate*0.128015;
_fasta_part_MG94xREV_3x4[34][35]:=_fasta_part_Shared_GT*nonSynRate*0.38188;
_fasta_part_MG94xREV_3x4[34][38]:=_fasta_part_Shared_AC*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[34][42]:=nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[34][46]:=_fasta_part_Shared_AT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[35][3]:=nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[35][19]:=_fasta_part_Shared_CG*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[35][32]:=_fasta_part_Shared_AT*nonSynRate*0.360235;
_fasta_part_MG94xREV_3x4[35][33]:=_fasta_part_Shared_CT*synRate*0.128015;
_fasta_part_MG94xREV_3x4[35][34]:=_fasta_part_Shared_GT*nonSynRate*0.12987;
_fasta_part_MG94xREV_3x4[35][39]:=_fasta_part_Shared_AC*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[35][43]:=nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[35][47]:=_fasta_part_Shared_AT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[35][49]:=_fasta_part_Shared_GT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[36][4]:=nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[36][20]:=_fasta_part_Shared_CG*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[36][32]:=_fasta_part_Shared_AC*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[36][37]:=_fasta_part_Shared_AC*synRate*0.128015;
_fasta_part_MG94xREV_3x4[36][38]:=synRate*0.12987;
_fasta_part_MG94xREV_3x4[36][39]:=_fasta_part_Shared_AT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[36][40]:=_fasta_part_Shared_CG*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[36][44]:=_fasta_part_Shared_CT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[36][50]:=_fasta_part_Shared_GT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[37][5]:=nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[37][21]:=_fasta_part_Shared_CG*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[37][33]:=_fasta_part_Shared_AC*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[37][36]:=_fasta_part_Shared_AC*synRate*0.360235;
_fasta_part_MG94xREV_3x4[37][38]:=_fasta_part_Shared_CG*synRate*0.12987;
_fasta_part_MG94xREV_3x4[37][39]:=_fasta_part_Shared_CT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[37][41]:=_fasta_part_Shared_CG*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[37][45]:=_fasta_part_Shared_CT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[37][51]:=_fasta_part_Shared_GT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[38][6]:=nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[38][22]:=_fasta_part_Shared_CG*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[38][34]:=_fasta_part_Shared_AC*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[38][36]:=synRate*0.360235;
_fasta_part_MG94xREV_3x4[38][37]:=_fasta_part_Shared_CG*synRate*0.128015;
_fasta_part_MG94xREV_3x4[38][39]:=_fasta_part_Shared_GT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[38][42]:=_fasta_part_Shared_CG*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[38][46]:=_fasta_part_Shared_CT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[38][52]:=_fasta_part_Shared_GT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[39][7]:=nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[39][23]:=_fasta_part_Shared_CG*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[39][35]:=_fasta_part_Shared_AC*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[39][36]:=_fasta_part_Shared_AT*synRate*0.360235;
_fasta_part_MG94xREV_3x4[39][37]:=_fasta_part_Shared_CT*synRate*0.128015;
_fasta_part_MG94xREV_3x4[39][38]:=_fasta_part_Shared_GT*synRate*0.12987;
_fasta_part_MG94xREV_3x4[39][43]:=_fasta_part_Shared_CG*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[39][47]:=_fasta_part_Shared_CT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[39][53]:=_fasta_part_Shared_GT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[40][8]:=nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[40][24]:=_fasta_part_Shared_CG*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[40][32]:=nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[40][36]:=_fasta_part_Shared_CG*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[40][41]:=_fasta_part_Shared_AC*synRate*0.128015;
_fasta_part_MG94xREV_3x4[40][42]:=synRate*0.12987;
_fasta_part_MG94xREV_3x4[40][43]:=_fasta_part_Shared_AT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[40][44]:=_fasta_part_Shared_GT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[41][9]:=nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[41][25]:=_fasta_part_Shared_CG*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[41][33]:=nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[41][37]:=_fasta_part_Shared_CG*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[41][40]:=_fasta_part_Shared_AC*synRate*0.360235;
_fasta_part_MG94xREV_3x4[41][42]:=_fasta_part_Shared_CG*synRate*0.12987;
_fasta_part_MG94xREV_3x4[41][43]:=_fasta_part_Shared_CT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[41][45]:=_fasta_part_Shared_GT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[41][54]:=_fasta_part_Shared_GT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[42][10]:=nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[42][26]:=_fasta_part_Shared_CG*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[42][34]:=nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[42][38]:=_fasta_part_Shared_CG*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[42][40]:=synRate*0.360235;
_fasta_part_MG94xREV_3x4[42][41]:=_fasta_part_Shared_CG*synRate*0.128015;
_fasta_part_MG94xREV_3x4[42][43]:=_fasta_part_Shared_GT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[42][46]:=_fasta_part_Shared_GT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[42][55]:=_fasta_part_Shared_GT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[43][11]:=nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[43][27]:=_fasta_part_Shared_CG*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[43][35]:=nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[43][39]:=_fasta_part_Shared_CG*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[43][40]:=_fasta_part_Shared_AT*synRate*0.360235;
_fasta_part_MG94xREV_3x4[43][41]:=_fasta_part_Shared_CT*synRate*0.128015;
_fasta_part_MG94xREV_3x4[43][42]:=_fasta_part_Shared_GT*synRate*0.12987;
_fasta_part_MG94xREV_3x4[43][47]:=_fasta_part_Shared_GT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[43][56]:=_fasta_part_Shared_GT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[44][12]:=nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[44][28]:=_fasta_part_Shared_CG*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[44][32]:=_fasta_part_Shared_AT*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[44][36]:=_fasta_part_Shared_CT*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[44][40]:=_fasta_part_Shared_GT*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[44][45]:=_fasta_part_Shared_AC*synRate*0.128015;
_fasta_part_MG94xREV_3x4[44][46]:=synRate*0.12987;
_fasta_part_MG94xREV_3x4[44][47]:=_fasta_part_Shared_AT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[44][57]:=_fasta_part_Shared_GT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[45][13]:=nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[45][29]:=_fasta_part_Shared_CG*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[45][33]:=_fasta_part_Shared_AT*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[45][37]:=_fasta_part_Shared_CT*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[45][41]:=_fasta_part_Shared_GT*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[45][44]:=_fasta_part_Shared_AC*synRate*0.360235;
_fasta_part_MG94xREV_3x4[45][46]:=_fasta_part_Shared_CG*synRate*0.12987;
_fasta_part_MG94xREV_3x4[45][47]:=_fasta_part_Shared_CT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[45][58]:=_fasta_part_Shared_GT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[46][14]:=nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[46][30]:=_fasta_part_Shared_CG*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[46][34]:=_fasta_part_Shared_AT*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[46][38]:=_fasta_part_Shared_CT*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[46][42]:=_fasta_part_Shared_GT*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[46][44]:=synRate*0.360235;
_fasta_part_MG94xREV_3x4[46][45]:=_fasta_part_Shared_CG*synRate*0.128015;
_fasta_part_MG94xREV_3x4[46][47]:=_fasta_part_Shared_GT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[46][59]:=_fasta_part_Shared_GT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[47][15]:=nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[47][31]:=_fasta_part_Shared_CG*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[47][35]:=_fasta_part_Shared_AT*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[47][39]:=_fasta_part_Shared_CT*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[47][43]:=_fasta_part_Shared_GT*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[47][44]:=_fasta_part_Shared_AT*synRate*0.360235;
_fasta_part_MG94xREV_3x4[47][45]:=_fasta_part_Shared_CT*synRate*0.128015;
_fasta_part_MG94xREV_3x4[47][46]:=_fasta_part_Shared_GT*synRate*0.12987;
_fasta_part_MG94xREV_3x4[47][60]:=_fasta_part_Shared_GT*nonSynRate*0.183673;
_fasta_part_MG94xREV_3x4[48][1]:=_fasta_part_Shared_AT*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[48][17]:=_fasta_part_Shared_CT*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[48][33]:=_fasta_part_Shared_GT*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[48][49]:=_fasta_part_Shared_CT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[48][51]:=_fasta_part_Shared_AC*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[48][54]:=nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[48][58]:=_fasta_part_Shared_AT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[49][3]:=_fasta_part_Shared_AT*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[49][19]:=_fasta_part_Shared_CT*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[49][35]:=_fasta_part_Shared_GT*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[49][48]:=_fasta_part_Shared_CT*synRate*0.128015;
_fasta_part_MG94xREV_3x4[49][53]:=_fasta_part_Shared_AC*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[49][56]:=nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[49][60]:=_fasta_part_Shared_AT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[50][4]:=_fasta_part_Shared_AT*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[50][20]:=_fasta_part_Shared_CT*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[50][36]:=_fasta_part_Shared_GT*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[50][51]:=_fasta_part_Shared_AC*synRate*0.128015;
_fasta_part_MG94xREV_3x4[50][52]:=synRate*0.12987;
_fasta_part_MG94xREV_3x4[50][53]:=_fasta_part_Shared_AT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[50][57]:=_fasta_part_Shared_CT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[51][5]:=_fasta_part_Shared_AT*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[51][21]:=_fasta_part_Shared_CT*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[51][37]:=_fasta_part_Shared_GT*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[51][48]:=_fasta_part_Shared_AC*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[51][50]:=_fasta_part_Shared_AC*synRate*0.360235;
_fasta_part_MG94xREV_3x4[51][52]:=_fasta_part_Shared_CG*synRate*0.12987;
_fasta_part_MG94xREV_3x4[51][53]:=_fasta_part_Shared_CT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[51][54]:=_fasta_part_Shared_CG*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[51][58]:=_fasta_part_Shared_CT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[52][6]:=_fasta_part_Shared_AT*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[52][22]:=_fasta_part_Shared_CT*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[52][38]:=_fasta_part_Shared_GT*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[52][50]:=synRate*0.360235;
_fasta_part_MG94xREV_3x4[52][51]:=_fasta_part_Shared_CG*synRate*0.128015;
_fasta_part_MG94xREV_3x4[52][53]:=_fasta_part_Shared_GT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[52][55]:=_fasta_part_Shared_CG*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[52][59]:=_fasta_part_Shared_CT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[53][7]:=_fasta_part_Shared_AT*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[53][23]:=_fasta_part_Shared_CT*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[53][39]:=_fasta_part_Shared_GT*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[53][49]:=_fasta_part_Shared_AC*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[53][50]:=_fasta_part_Shared_AT*synRate*0.360235;
_fasta_part_MG94xREV_3x4[53][51]:=_fasta_part_Shared_CT*synRate*0.128015;
_fasta_part_MG94xREV_3x4[53][52]:=_fasta_part_Shared_GT*synRate*0.12987;
_fasta_part_MG94xREV_3x4[53][56]:=_fasta_part_Shared_CG*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[53][60]:=_fasta_part_Shared_CT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[54][9]:=_fasta_part_Shared_AT*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[54][25]:=_fasta_part_Shared_CT*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[54][41]:=_fasta_part_Shared_GT*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[54][48]:=nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[54][51]:=_fasta_part_Shared_CG*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[54][55]:=_fasta_part_Shared_CG*nonSynRate*0.12987;
_fasta_part_MG94xREV_3x4[54][56]:=_fasta_part_Shared_CT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[54][58]:=_fasta_part_Shared_GT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[55][10]:=_fasta_part_Shared_AT*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[55][26]:=_fasta_part_Shared_CT*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[55][42]:=_fasta_part_Shared_GT*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[55][52]:=_fasta_part_Shared_CG*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[55][54]:=_fasta_part_Shared_CG*nonSynRate*0.128015;
_fasta_part_MG94xREV_3x4[55][56]:=_fasta_part_Shared_GT*nonSynRate*0.38188;
_fasta_part_MG94xREV_3x4[55][59]:=_fasta_part_Shared_GT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[56][11]:=_fasta_part_Shared_AT*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[56][27]:=_fasta_part_Shared_CT*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[56][43]:=_fasta_part_Shared_GT*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[56][49]:=nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[56][53]:=_fasta_part_Shared_CG*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[56][54]:=_fasta_part_Shared_CT*synRate*0.128015;
_fasta_part_MG94xREV_3x4[56][55]:=_fasta_part_Shared_GT*nonSynRate*0.12987;
_fasta_part_MG94xREV_3x4[56][60]:=_fasta_part_Shared_GT*nonSynRate*0.259431;
_fasta_part_MG94xREV_3x4[57][12]:=_fasta_part_Shared_AT*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[57][28]:=_fasta_part_Shared_CT*synRate*0.134818;
_fasta_part_MG94xREV_3x4[57][44]:=_fasta_part_Shared_GT*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[57][50]:=_fasta_part_Shared_CT*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[57][58]:=_fasta_part_Shared_AC*nonSynRate*0.128015;
_fasta_part_MG94xREV_3x4[57][59]:=synRate*0.12987;
_fasta_part_MG94xREV_3x4[57][60]:=_fasta_part_Shared_AT*nonSynRate*0.38188;
_fasta_part_MG94xREV_3x4[58][13]:=_fasta_part_Shared_AT*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[58][29]:=_fasta_part_Shared_CT*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[58][45]:=_fasta_part_Shared_GT*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[58][48]:=_fasta_part_Shared_AT*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[58][51]:=_fasta_part_Shared_CT*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[58][54]:=_fasta_part_Shared_GT*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[58][57]:=_fasta_part_Shared_AC*nonSynRate*0.360235;
_fasta_part_MG94xREV_3x4[58][59]:=_fasta_part_Shared_CG*nonSynRate*0.12987;
_fasta_part_MG94xREV_3x4[58][60]:=_fasta_part_Shared_CT*synRate*0.38188;
_fasta_part_MG94xREV_3x4[59][14]:=_fasta_part_Shared_AT*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[59][30]:=_fasta_part_Shared_CT*synRate*0.134818;
_fasta_part_MG94xREV_3x4[59][46]:=_fasta_part_Shared_GT*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[59][52]:=_fasta_part_Shared_CT*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[59][55]:=_fasta_part_Shared_GT*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[59][57]:=synRate*0.360235;
_fasta_part_MG94xREV_3x4[59][58]:=_fasta_part_Shared_CG*nonSynRate*0.128015;
_fasta_part_MG94xREV_3x4[59][60]:=_fasta_part_Shared_GT*nonSynRate*0.38188;
_fasta_part_MG94xREV_3x4[60][15]:=_fasta_part_Shared_AT*nonSynRate*0.443414;
_fasta_part_MG94xREV_3x4[60][31]:=_fasta_part_Shared_CT*nonSynRate*0.134818;
_fasta_part_MG94xREV_3x4[60][47]:=_fasta_part_Shared_GT*nonSynRate*0.238095;
_fasta_part_MG94xREV_3x4[60][49]:=_fasta_part_Shared_AT*nonSynRate*0.325294;
_fasta_part_MG94xREV_3x4[60][53]:=_fasta_part_Shared_CT*nonSynRate*0.226964;
_fasta_part_MG94xREV_3x4[60][56]:=_fasta_part_Shared_GT*nonSynRate*0.188312;
_fasta_part_MG94xREV_3x4[60][57]:=_fasta_part_Shared_AT*nonSynRate*0.360235;
_fasta_part_MG94xREV_3x4[60][58]:=_fasta_part_Shared_CT*synRate*0.128015;
_fasta_part_MG94xREV_3x4[60][59]:=_fasta_part_Shared_GT*nonSynRate*0.12987;


_fasta_part_Freqs={
{   0.0542236230763}
{   0.0192691673421}
{    0.019548430637}
{   0.0574816948491}
{   0.0378328320704}
{   0.0134444570619}
{   0.0136393042657}
{    0.040106049448}
{   0.0313899110774}
{   0.0111548696876}
{   0.0113165344657}
{   0.0332760001551}
{   0.0432448857044}
{   0.0153677104563}
{   0.0155904308977}
{    0.045843290854}
{   0.0164864014374}
{  0.00585868686274}
{    0.005943595368}
{   0.0174770006654}
{    0.011502869444}
{  0.00408771497838}
{  0.00414695722444}
{   0.0121940289814}
{  0.00954393391197}
{  0.00339157823138}
{  0.00344073153908}
{   0.0101173891685}
{   0.0131483752909}
{  0.00467246984586}
{  0.00474018680014}
{   0.0139384064242}
{   0.0291158924468}
{   0.0103467634961}
{   0.0104967165903}
{   0.0308653452119}
{   0.0203147006236}
{  0.00721912966365}
{  0.00732375473124}
{   0.0215353264121}
{   0.0168551126427}
{  0.00598971384899}
{  0.00607652129608}
{    0.017867866192}
{   0.0232207545275}
{  0.00825183894796}
{  0.00837143081677}
{   0.0246159929969}
{  0.00798178898273}
{   0.0238104091635}
{    0.015671340481}
{  0.00556904288339}
{  0.00564975364981}
{   0.0166129660893}
{   0.0046206363978}
{  0.00468760214269}
{    0.013783782491}
{   0.0179131534926}
{  0.00636570433128}
{  0.00645796091579}
{   0.0189894803119}
}
;
Model _fasta_part_MG94xREV_3x4_model=(_fasta_part_MG94xREV_3x4,_fasta_part_Freqs,0);

UseModel (_fasta_part_MG94xREV_3x4_model);
Tree _fasta_tree=(D1_5,D1_10,D1_8,D1_9,D1_7,D1_0,(D1_1,D1_2,D1_3,D1_4)Node7,D1_6);

_fasta_tree.D1_3.synRate=0.006957042058010423;
_fasta_tree.D1_2.synRate=0.006951740248521895;
_fasta_tree.D1_1.synRate=0;
_fasta_tree.D1_0.synRate=0;
_fasta_tree.Node7.synRate=0.1871920897518702;
_fasta_tree.D1_4.synRate=0;
_fasta_tree.D1_6.synRate=0;
_fasta_tree.D1_7.synRate=0;
_fasta_tree.D1_5.synRate=0.006995758267155891;
_fasta_tree.D1_10.synRate=0;
_fasta_tree.D1_8.synRate=0.006987551134318812;
_fasta_tree.D1_9.synRate=0;
_fasta_tree.D1_9.nonSynRate:=omega*_fasta_tree.D1_9.synRate;
_fasta_tree.D1_10.nonSynRate:=omega*_fasta_tree.D1_10.synRate;
_fasta_tree.D1_8.nonSynRate:=omega*_fasta_tree.D1_8.synRate;
_fasta_tree.D1_5.nonSynRate:=omega*_fasta_tree.D1_5.synRate;
_fasta_tree.D1_7.nonSynRate:=omega*_fasta_tree.D1_7.synRate;
_fasta_tree.D1_0.nonSynRate:=omega*_fasta_tree.D1_0.synRate;
_fasta_tree.D1_1.nonSynRate:=omega*_fasta_tree.D1_1.synRate;
_fasta_tree.D1_2.nonSynRate:=omega*_fasta_tree.D1_2.synRate;
_fasta_tree.D1_3.nonSynRate:=omega*_fasta_tree.D1_3.synRate;
_fasta_tree.D1_4.nonSynRate:=omega*_fasta_tree.D1_4.synRate;
_fasta_tree.Node7.nonSynRate:=omega*_fasta_tree.Node7.synRate;
_fasta_tree.D1_6.nonSynRate:=omega*_fasta_tree.D1_6.synRate;
DataSet _fasta = ReadDataFile(USE_NEXUS_FILE_DATA);
DataSetFilter _fasta_part = CreateFilter(_fasta,3,"0-881","0-10","TAA,TAG,TGA");
LikelihoodFunction _fasta_LF = (_fasta_part,_fasta_tree);

fprintf (stdout, _fasta_LF);

Tree _fasta_tree2=(D1_5,D1_10,D1_8,D1_9,D1_7,D1_0,(D1_1,D1_2,D1_3,D1_4)Node7,D1_6);
LikelihoodFunction _fasta_LF_2 = (_fasta_part,_fasta_tree2);
_fasta_tree.Node7.synRate=0.1871920897518702;
fprintf (stdout, _fasta_LF);

Tree _fasta_tree3=(D1_5,D1_10,D1_8,D1_9,D1_7,D1_0,D1_6,(D1_1,D1_2,D1_3,D1_4)Node7);
LikelihoodFunction _fasta_LF_3 = (_fasta_part,_fasta_tree3);
/* this should throw an error */
_fasta_tree.Node7.synRate=0.1871920897518702;
fprintf (stdout, _fasta_LF);

END;