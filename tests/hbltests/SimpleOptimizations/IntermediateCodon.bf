/* test preamble */

	_testDescription 		= " fit the MG94xREV model to an Influenza A alignment with 349 sequences and 329 nucleotides";
	_expectedLL 			= -11402.1903626064;
	ExecuteAFile 			("../Shared/TestInstrumentation.bf");
	startTestTimer 			(_testDescription);

/* end test preamble */

global LargeNuc_part_Shared_AC	=	1;
global LargeNuc_part_Shared_AT	=	1;
global LargeNuc_part_Shared_CG	=	1;
global LargeNuc_part_Shared_CT	=	1;
global LargeNuc_part_Shared_GT	=	1;
global LargeNuc_part_Shared_R	=	1;

LargeNuc_part_MG94xREV_3x4={61,61};

LargeNuc_part_MG94xREV_3x4[0][1]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[0][2]:=synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[0][3]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[0][4]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[0][8]:=LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[0][12]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[0][16]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[0][32]:=LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[1][0]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[1][2]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[1][3]:=LargeNuc_part_Shared_CT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[1][5]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[1][9]:=LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[1][13]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[1][17]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[1][33]:=LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[1][48]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[2][0]:=synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[2][1]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[2][3]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[2][6]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[2][10]:=LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[2][14]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[2][18]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[2][34]:=LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[3][0]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[3][1]:=LargeNuc_part_Shared_CT*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[3][2]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[3][7]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[3][11]:=LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[3][15]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[3][19]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[3][35]:=LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[3][49]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[4][0]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[4][5]:=LargeNuc_part_Shared_AC*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[4][6]:=synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[4][7]:=LargeNuc_part_Shared_AT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[4][8]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[4][12]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[4][20]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[4][36]:=LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[4][50]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[5][1]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[5][4]:=LargeNuc_part_Shared_AC*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[5][6]:=LargeNuc_part_Shared_CG*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[5][7]:=LargeNuc_part_Shared_CT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[5][9]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[5][13]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[5][21]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[5][37]:=LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[5][51]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[6][2]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[6][4]:=synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[6][5]:=LargeNuc_part_Shared_CG*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[6][7]:=LargeNuc_part_Shared_GT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[6][10]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[6][14]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[6][22]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[6][38]:=LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[6][52]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[7][3]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[7][4]:=LargeNuc_part_Shared_AT*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[7][5]:=LargeNuc_part_Shared_CT*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[7][6]:=LargeNuc_part_Shared_GT*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[7][11]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[7][15]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[7][23]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[7][39]:=LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[7][53]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[8][0]:=LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[8][4]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[8][9]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[8][10]:=synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[8][11]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[8][12]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[8][24]:=LargeNuc_part_Shared_AC*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[8][40]:=LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[9][1]:=LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[9][5]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[9][8]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[9][10]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[9][11]:=LargeNuc_part_Shared_CT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[9][13]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[9][25]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[9][41]:=LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[9][54]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[10][2]:=LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[10][6]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[10][8]:=synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[10][9]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[10][11]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[10][14]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[10][26]:=LargeNuc_part_Shared_AC*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[10][42]:=LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[10][55]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[11][3]:=LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[11][7]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[11][8]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[11][9]:=LargeNuc_part_Shared_CT*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[11][10]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[11][15]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[11][27]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[11][43]:=LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[11][56]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[12][0]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[12][4]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[12][8]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[12][13]:=LargeNuc_part_Shared_AC*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[12][14]:=LargeNuc_part_Shared_R*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[12][15]:=LargeNuc_part_Shared_AT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[12][28]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[12][44]:=LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[12][57]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[13][1]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[13][5]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[13][9]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[13][12]:=LargeNuc_part_Shared_AC*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[13][14]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[13][15]:=LargeNuc_part_Shared_CT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[13][29]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[13][45]:=LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[13][58]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[14][2]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[14][6]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[14][10]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[14][12]:=LargeNuc_part_Shared_R*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[14][13]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[14][15]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[14][30]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[14][46]:=LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[14][59]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[15][3]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[15][7]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[15][11]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[15][12]:=LargeNuc_part_Shared_AT*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[15][13]:=LargeNuc_part_Shared_CT*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[15][14]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[15][31]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[15][47]:=LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[15][60]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[16][0]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[16][17]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[16][18]:=synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[16][19]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[16][20]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[16][24]:=LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[16][28]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[16][32]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[17][1]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[17][16]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[17][18]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[17][19]:=LargeNuc_part_Shared_CT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[17][21]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[17][25]:=LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[17][29]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[17][33]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[17][48]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[18][2]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[18][16]:=synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[18][17]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[18][19]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[18][22]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[18][26]:=LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[18][30]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[18][34]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[19][3]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[19][16]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[19][17]:=LargeNuc_part_Shared_CT*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[19][18]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[19][23]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[19][27]:=LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[19][31]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[19][35]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[19][49]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[20][4]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[20][16]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[20][21]:=LargeNuc_part_Shared_AC*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[20][22]:=synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[20][23]:=LargeNuc_part_Shared_AT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[20][24]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[20][28]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[20][36]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[20][50]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[21][5]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[21][17]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[21][20]:=LargeNuc_part_Shared_AC*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[21][22]:=LargeNuc_part_Shared_CG*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[21][23]:=LargeNuc_part_Shared_CT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[21][25]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[21][29]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[21][37]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[21][51]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[22][6]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[22][18]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[22][20]:=synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[22][21]:=LargeNuc_part_Shared_CG*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[22][23]:=LargeNuc_part_Shared_GT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[22][26]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[22][30]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[22][38]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[22][52]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[23][7]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[23][19]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[23][20]:=LargeNuc_part_Shared_AT*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[23][21]:=LargeNuc_part_Shared_CT*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[23][22]:=LargeNuc_part_Shared_GT*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[23][27]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[23][31]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[23][39]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[23][53]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[24][8]:=LargeNuc_part_Shared_AC*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[24][16]:=LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[24][20]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[24][25]:=LargeNuc_part_Shared_AC*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[24][26]:=synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[24][27]:=LargeNuc_part_Shared_AT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[24][28]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[24][40]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[25][9]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[25][17]:=LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[25][21]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[25][24]:=LargeNuc_part_Shared_AC*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[25][26]:=LargeNuc_part_Shared_CG*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[25][27]:=LargeNuc_part_Shared_CT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[25][29]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[25][41]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[25][54]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[26][10]:=LargeNuc_part_Shared_AC*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[26][18]:=LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[26][22]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[26][24]:=synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[26][25]:=LargeNuc_part_Shared_CG*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[26][27]:=LargeNuc_part_Shared_GT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[26][30]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[26][42]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[26][55]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[27][11]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[27][19]:=LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[27][23]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[27][24]:=LargeNuc_part_Shared_AT*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[27][25]:=LargeNuc_part_Shared_CT*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[27][26]:=LargeNuc_part_Shared_GT*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[27][31]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[27][43]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[27][56]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[28][12]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[28][16]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[28][20]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[28][24]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[28][29]:=LargeNuc_part_Shared_AC*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[28][30]:=synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[28][31]:=LargeNuc_part_Shared_AT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[28][44]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[28][57]:=LargeNuc_part_Shared_CT*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[29][13]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[29][17]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[29][21]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[29][25]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[29][28]:=LargeNuc_part_Shared_AC*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[29][30]:=LargeNuc_part_Shared_CG*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[29][31]:=LargeNuc_part_Shared_CT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[29][45]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[29][58]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[30][14]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[30][18]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[30][22]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[30][26]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[30][28]:=synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[30][29]:=LargeNuc_part_Shared_CG*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[30][31]:=LargeNuc_part_Shared_GT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[30][46]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[30][59]:=LargeNuc_part_Shared_CT*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[31][15]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[31][19]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[31][23]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[31][27]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[31][28]:=LargeNuc_part_Shared_AT*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[31][29]:=LargeNuc_part_Shared_CT*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[31][30]:=LargeNuc_part_Shared_GT*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[31][47]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[31][60]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[32][0]:=LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[32][16]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[32][33]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[32][34]:=synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[32][35]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[32][36]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[32][40]:=LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[32][44]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[33][1]:=LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[33][17]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[33][32]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[33][34]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[33][35]:=LargeNuc_part_Shared_CT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[33][37]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[33][41]:=LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[33][45]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[33][48]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[34][2]:=LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[34][18]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[34][32]:=synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[34][33]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[34][35]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[34][38]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[34][42]:=LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[34][46]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[35][3]:=LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[35][19]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[35][32]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[35][33]:=LargeNuc_part_Shared_CT*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[35][34]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[35][39]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[35][43]:=LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[35][47]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[35][49]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[36][4]:=LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[36][20]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[36][32]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[36][37]:=LargeNuc_part_Shared_AC*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[36][38]:=synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[36][39]:=LargeNuc_part_Shared_AT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[36][40]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[36][44]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[36][50]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[37][5]:=LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[37][21]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[37][33]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[37][36]:=LargeNuc_part_Shared_AC*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[37][38]:=LargeNuc_part_Shared_CG*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[37][39]:=LargeNuc_part_Shared_CT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[37][41]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[37][45]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[37][51]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[38][6]:=LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[38][22]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[38][34]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[38][36]:=synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[38][37]:=LargeNuc_part_Shared_CG*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[38][39]:=LargeNuc_part_Shared_GT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[38][42]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[38][46]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[38][52]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[39][7]:=LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[39][23]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[39][35]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[39][36]:=LargeNuc_part_Shared_AT*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[39][37]:=LargeNuc_part_Shared_CT*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[39][38]:=LargeNuc_part_Shared_GT*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[39][43]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[39][47]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[39][53]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[40][8]:=LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[40][24]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[40][32]:=LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[40][36]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[40][41]:=LargeNuc_part_Shared_AC*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[40][42]:=synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[40][43]:=LargeNuc_part_Shared_AT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[40][44]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[41][9]:=LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[41][25]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[41][33]:=LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[41][37]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[41][40]:=LargeNuc_part_Shared_AC*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[41][42]:=LargeNuc_part_Shared_CG*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[41][43]:=LargeNuc_part_Shared_CT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[41][45]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[41][54]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[42][10]:=LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[42][26]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[42][34]:=LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[42][38]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[42][40]:=synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[42][41]:=LargeNuc_part_Shared_CG*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[42][43]:=LargeNuc_part_Shared_GT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[42][46]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[42][55]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[43][11]:=LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[43][27]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[43][35]:=LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[43][39]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[43][40]:=LargeNuc_part_Shared_AT*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[43][41]:=LargeNuc_part_Shared_CT*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[43][42]:=LargeNuc_part_Shared_GT*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[43][47]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[43][56]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[44][12]:=LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[44][28]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[44][32]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[44][36]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[44][40]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[44][45]:=LargeNuc_part_Shared_AC*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[44][46]:=synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[44][47]:=LargeNuc_part_Shared_AT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[44][57]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[45][13]:=LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[45][29]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[45][33]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[45][37]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[45][41]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[45][44]:=LargeNuc_part_Shared_AC*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[45][46]:=LargeNuc_part_Shared_CG*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[45][47]:=LargeNuc_part_Shared_CT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[45][58]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[46][14]:=LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[46][30]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[46][34]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[46][38]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[46][42]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[46][44]:=synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[46][45]:=LargeNuc_part_Shared_CG*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[46][47]:=LargeNuc_part_Shared_GT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[46][59]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[47][15]:=LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[47][31]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[47][35]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[47][39]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[47][43]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[47][44]:=LargeNuc_part_Shared_AT*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[47][45]:=LargeNuc_part_Shared_CT*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[47][46]:=LargeNuc_part_Shared_GT*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[47][60]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.176322;
LargeNuc_part_MG94xREV_3x4[48][1]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[48][17]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[48][33]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[48][49]:=LargeNuc_part_Shared_CT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[48][51]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[48][54]:=LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[48][58]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[49][3]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[49][19]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[49][35]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[49][48]:=LargeNuc_part_Shared_CT*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[49][53]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[49][56]:=LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[49][60]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[50][4]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[50][20]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[50][36]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[50][51]:=LargeNuc_part_Shared_AC*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[50][52]:=synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[50][53]:=LargeNuc_part_Shared_AT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[50][57]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[51][5]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[51][21]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[51][37]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[51][48]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[51][50]:=LargeNuc_part_Shared_AC*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[51][52]:=LargeNuc_part_Shared_CG*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[51][53]:=LargeNuc_part_Shared_CT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[51][54]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[51][58]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[52][6]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[52][22]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[52][38]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[52][50]:=synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[52][51]:=LargeNuc_part_Shared_CG*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[52][53]:=LargeNuc_part_Shared_GT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[52][55]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[52][59]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[53][7]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[53][23]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[53][39]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[53][49]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[53][50]:=LargeNuc_part_Shared_AT*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[53][51]:=LargeNuc_part_Shared_CT*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[53][52]:=LargeNuc_part_Shared_GT*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[53][56]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[53][60]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[54][9]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[54][25]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[54][41]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[54][48]:=LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[54][51]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[54][55]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[54][56]:=LargeNuc_part_Shared_CT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[54][58]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[55][10]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[55][26]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[55][42]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[55][52]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[55][54]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[55][56]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[55][59]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[56][11]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[56][27]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[56][43]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[56][49]:=LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[56][53]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[56][54]:=LargeNuc_part_Shared_CT*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[56][55]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[56][60]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.234948;
LargeNuc_part_MG94xREV_3x4[57][12]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[57][28]:=LargeNuc_part_Shared_CT*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[57][44]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[57][50]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[57][58]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[57][59]:=synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[57][60]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[58][13]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[58][29]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[58][45]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[58][48]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[58][51]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[58][54]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[58][57]:=LargeNuc_part_Shared_AC*LargeNuc_part_Shared_R*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[58][59]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.179823;
LargeNuc_part_MG94xREV_3x4[58][60]:=LargeNuc_part_Shared_CT*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[59][14]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[59][30]:=LargeNuc_part_Shared_CT*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[59][46]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[59][52]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[59][55]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[59][57]:=synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[59][58]:=LargeNuc_part_Shared_CG*LargeNuc_part_Shared_R*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[59][60]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.277297;
LargeNuc_part_MG94xREV_3x4[60][15]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.378681;
LargeNuc_part_MG94xREV_3x4[60][31]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.184056;
LargeNuc_part_MG94xREV_3x4[60][47]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.260941;
LargeNuc_part_MG94xREV_3x4[60][49]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.312643;
LargeNuc_part_MG94xREV_3x4[60][53]:=LargeNuc_part_Shared_CT*LargeNuc_part_Shared_R*synRate*0.223626;
LargeNuc_part_MG94xREV_3x4[60][56]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.228782;
LargeNuc_part_MG94xREV_3x4[60][57]:=LargeNuc_part_Shared_AT*LargeNuc_part_Shared_R*synRate*0.292965;
LargeNuc_part_MG94xREV_3x4[60][58]:=LargeNuc_part_Shared_CT*synRate*0.249915;
LargeNuc_part_MG94xREV_3x4[60][59]:=LargeNuc_part_Shared_GT*LargeNuc_part_Shared_R*synRate*0.179823;

LargeNuc_part_Freqs={
{   0.0360502764176}
{   0.0307528786046}
{   0.0221278618943}
{   0.0341222936813}
{   0.0257859197609}
{   0.0219968149738}
{    0.015827542199}
{   0.0244068787914}
{   0.0263804309771}
{   0.0225039659052}
{   0.0161924565185}
{   0.0249695953177}
{   0.0270914342559}
{   0.0231104910055}
{   0.0166288743196}
{    0.025642573866}
{   0.0175220734967}
{    0.014947297294}
{   0.0107551470048}
{    0.016584986224}
{    0.012533129455}
{    0.010691452243}
{  0.00769290516584}
{   0.0118628528407}
{   0.0128220889377}
{   0.0109379506551}
{  0.00787027011728}
{   0.0121363586584}
{   0.0131676688595}
{   0.0112327494317}
{  0.00808238901191}
{   0.0124634568323}
{   0.0248414888717}
{   0.0211911632183}
{   0.0152478452214}
{   0.0235129564318}
{   0.0177685361234}
{   0.0151575435388}
{   0.0109064271477}
{   0.0168182679341}
{   0.0181782013252}
{   0.0155070106017}
{   0.0111578819466}
{   0.0172060240823}
{   0.0186681387624}
{   0.0159249543189}
{    0.011458608294}
{   0.0176697594758}
{   0.0143192328467}
{   0.0158881083871}
{   0.0120065049509}
{   0.0102422124298}
{  0.00736966009109}
{   0.0113643924189}
{   0.0104783533247}
{  0.00753957241628}
{   0.0116264059062}
{   0.0126143819006}
{   0.0107607650706}
{  0.00774277837276}
{    0.011939759874}
}
;
Model LargeNuc_part_MG94xREV_3x4_model=(LargeNuc_part_MG94xREV_3x4,LargeNuc_part_Freqs,0);

UseModel (LargeNuc_part_MG94xREV_3x4_model);

DataSet flu 				= 	ReadDataFile(PATH_TO_CURRENT_BF + "/../data/fluHA.nex");
Tree LargeNuc_tree				=	DATAFILE_TREE;


DataSetFilter LargeNuc_part = CreateFilter(flu,3,"","","TAA,TAG,TGA");
VERBOSITY_LEVEL = 1;
LikelihoodFunction LargeNuc_LF = (LargeNuc_part,LargeNuc_tree);
//AUTO_PARALLELIZE_OPTIMIZE = 1;
OPTIMIZATION_PRECISION    = 0.001;

Optimize(res_LargeNuc_LF,LargeNuc_LF);

/* test epilogue */
	timeMatrix = endTestTimer 				  (_testDescription);
	if (logTestResult (Abs (res_LargeNuc_LF[1][0] - _expectedLL) < 2*OPTIMIZATION_PRECISION))
	{
		return timeMatrix;
	}
	return 0;
/* end test epilogue */
