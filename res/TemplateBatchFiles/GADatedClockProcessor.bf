branchColors     = {};
branchColors [0] = {{0,0,255}}*(1/255);   		/* blue */
branchColors [1] = {{128,0,128}}*(1/255); 		/* purple */
branchColors [2] = {{255,0,0}}*(1/255); 		/* red */


ExecuteAFile(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "GrabBag.bf");
ExecuteAFile(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TreeTools.ibf");
ExecuteAFile(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "WriteDelimitedFiles.bf");
ExecuteAFile(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "DescriptiveStatistics.bf");

SetDialogPrompt ("Locate the results file:");

fscanf (PROMPT_FOR_FILE, "Raw", theSamples);
basePath = LAST_FILE_PATH;
UseModel (USE_NO_MODEL);
sscanf 	(theSamples,"Tree",bestTree);
sscanf 	(theSamples,"Tree",givenTree);
sscanf 	(theSamples,"Tree",singleTree);

TREE_OUTPUT_OPTIONS = {};
tree_lng 			= BranchLength (bestTree,-1);
branchNames			= BranchName   (bestTree,-1);
iNodeCounter		= 0;

fprintf 		(stdout, "\n[PHASE 0.] Making tree plots\n");
TREE_OUTPUT_OPTIONS["__FONT_SIZE__"] = 16;
baseHeight 		= TipCount (givenTree)*52;

drawLetter			= "/drawletter {"+TREE_OUTPUT_OPTIONS["__FONT_SIZE__"]$4+" -"+TREE_OUTPUT_OPTIONS["__FONT_SIZE__"]$2+
					   " rmoveto 1 copy false charpath pathbbox 2 index 3 sub sub exch 3 index 3 sub sub exch  0.85 setgray 4 copy rectfill 0 setgray  3 index 3 index currentlinewidth 0.5 setlinewidth 7 3 roll rectstroke setlinewidth exch 1.5 add exch 1.5 add moveto show} def\n";

psString 		= PSTreeString (givenTree, "STRING_SUPPLIED_LENGTHS",{{612,baseHeight}}) ^ {{"showpage"}{"0 "+(baseHeight-50-TREE_OUTPUT_OPTIONS["__FONT_SIZE__"])+" moveto (Free rates model) drawletter\nshowpage\n"}} ;
psString2 		= PSTreeString (singleTree,"STRING_SUPPLIED_LENGTHS",{{612,baseHeight}}) ^ {{"showpage"}{"0 "+(baseHeight-50-TREE_OUTPUT_OPTIONS["__FONT_SIZE__"])+" moveto (Single model) drawletter\nshowpage\n"}};
psString3 		= PSTreeString (bestTree,  "STRING_SUPPLIED_LENGTHS",{{612,baseHeight}}) ^ {{"showpage"}{"0 "+(baseHeight-50-TREE_OUTPUT_OPTIONS["__FONT_SIZE__"])+" moveto (Multiple rate model) drawletter\nshowpage\n"}};

tName  = basePath + ".nwk";
fName  = basePath + "_treeshapes.ps";
fprintf (fName, CLEAR_FILE, drawLetter, psString, 
							psString2, 
							psString3);
	
fprintf (tName, CLEAR_FILE, Format(givenTree,0,1), "\n", Format(singleTree,0,1), "\n", Format(bestTree,0,1), "\n");


sscanf 	(theSamples,"String",units);
sscanf 	(theSamples,"NMatrix",tipDates);

rateAssignments = {};
internalDates   = {};
AIC				= {};
rateValues 		= {};
modelsRead		= 0;
bestAIC			= 1e100;
bestModelID		= -1;
bestModelRC		= 0;

while (!END_OF_FILE)
{
	modelRates = 0;
	sscanf (theSamples,"NMatrix,Number,NMatrix,NMatrix",modelSpec,anAIC,modelDates,modelRates);
	if (Rows (modelRates))
	{
		AIC[modelsRead] 			= anAIC;
		if (anAIC < bestAIC)
		{
			bestAIC 	= anAIC;
			bestModelID = modelsRead;
			bestModelRC	= Rows(modelRates);
		}
		
		internalDates[modelsRead]	= modelDates;
		rateValues[modelsRead]		= modelRates;
		rateAssignments[modelsRead]	= modelSpec;
		modelsRead 					= modelsRead + 1;
	}
}	

fprintf 		(stdout, "[PHASE 1.] Read ", modelsRead, " models. Best AIC found: ", bestAIC, " using ", bestModelRC," rates.\n");
sortedScores = 	({modelsRead,2}["_MATRIX_ELEMENT_COLUMN_ * _MATRIX_ELEMENT_ROW_ + (_MATRIX_ELEMENT_COLUMN_==0)*Exp(-0.5*(AIC[_MATRIX_ELEMENT_ROW_]-bestAIC))"])%0;

cumulWeight  = 0;
normConst	 = ({1,modelsRead}["1"]*sortedScores[-1][0])[0];
for (upto95p = modelsRead-1; upto95p>=0 && cumulWeight <= 0.95; upto95p = upto95p - 1)
{
	cumulWeight = cumulWeight + sortedScores[upto95p][0]/normConst;
}
ciSize 			= modelsRead-upto95p-1;
fprintf 		(stdout, "[PHASE 2.] ", ciSize , " models in the 95% confidence set.\n");
fprintf 		(stdout, "[PHASE 3.] Tabulating rate and date estimates\n");
branchCount		= BranchCount (givenTree)+1;
allCount		= BranchCount (givenTree) + TipCount (givenTree);

dateEstimates   = {ciSize, branchCount};
rateEstimates	= {ciSize, allCount};
lengthEstimates = {ciSize, branchCount};
dateInfo		= {branchCount,5};
lengthInfo		= {branchCount,5};
rootLengthInfo	= {};
branchSpanInfo	= {};


rateInfo		= {allCount,5};

k4 = 0;
normConst = normConst*cumulWeight;

treeAVL = givenTree ^ 0;
nodeNames		= {};
allNames		= {};
iNodes			= {};
parentNames		= {};
allNamesToI		= {};

iNodeCounter	= 0;
for (k3 = 1; k3 < Abs(treeAVL); k3=k3+1)
{
	if(Abs((treeAVL[k3])["Children"]))
	{
		nodeNames[Abs(nodeNames)] = (treeAVL[k3])["Name"];
		iNodes[(treeAVL[k3])["Name"]] = 1;
		allNamesToI [(treeAVL[k3])["Name"]] = iNodeCounter;
		iNodeCounter				  = iNodeCounter+1;
	}
	if ((treeAVL[k3])["Parent"])
	{
		allNames[Abs(allNames)] = (treeAVL[k3])["Name"];
		parentNames[Abs(parentNames)] = (treeAVL[(treeAVL[k3])["Parent"]])["Name"];
	}
	
}


for (k = modelsRead-1; k > upto95p; k = k - 1)
{
	thisModelID = sortedScores[k][1];
	myAW		= (sortedScores[k][0])/normConst;
	for (k3 = 0; k3 < branchCount; k3=k3+1)
	{
		dateEstimates [k4][k3] = (internalDates[thisModelID])[k3];
		dateInfo	  [k3][0]  = dateInfo[k3][0] + dateEstimates [k4][k3]*myAW;
	}
	rootD		= (internalDates[thisModelID])[branchCount-1];
	iNodeCounter= 0;
	/*tipCounter  = 0;*/
	for (k3 = 0; k3 < allCount; k3=k3+1)
	{
		rateEstimates [k4][k3] 			= (rateValues[thisModelID])[(rateAssignments[thisModelID])[k3]];
		rateInfo 	  [k3][0]  			= rateInfo[k3][0] + rateEstimates [k4][k3]*myAW;
		nodeName						= allNames[k3];
		parentD		  = (internalDates[thisModelID])[allNamesToI[parentNames[k3]]];
		/*if (iNodes[nodeName])
		{
			myD		     							= (internalDates[thisModelID])[iNodeCounter];
			lengthEstimates[k4][iNodeCounter] 		= myD-parentD;
			iNodeCounter 							= iNodeCounter+1;
			lengthInfo [k3][0] 						= lengthInfo [k3][0] + myAW*(myD-parentD);
		}
		else
		{
			myD		   = tipDates[tipCounter];
			tipCounter = tipCounter+1;
		
		}*/
		rootLengthInfo[allNames[k3]]	= (rootLengthInfo[allNames[k3]]) + myAW*(myD-rootD);
		branchSpanInfo[allNames[k3]]	= (branchSpanInfo[allNames[k3]]) + myAW*(myD-parentD);
	}
	k4 = k4 + 1;
}

/* generate date estimates */

fprintf 		(stdout, "[PHASE 4.] Writing CSV files\n");

labels = {{"Branch/Node","Model Averaged Mean", "2.5%", "Mean", "Median", "97.5%"}};

for (k3 = 0; k3 < branchCount; k3=k3+1)
{
	myDS = GatherDescriptiveStats (dateEstimates[-1][k3]);
	dateInfo[k3][1] = myDS["2.5%"];
	dateInfo[k3][2] = myDS["Mean"];
	dateInfo[k3][3] = myDS["Median"];
	dateInfo[k3][4] = myDS["97.5%"];
	myDS = GatherDescriptiveStats (lengthEstimates[-1][k3]);
	lengthEstimates[k3][1] = myDS["2.5%"];
	lengthEstimates[k3][2] = myDS["Mean"];
	lengthEstimates[k3][3] = myDS["Median"];
	lengthEstimates[k3][4] = myDS["97.5%"];
}

WriteSeparatedTable (basePath + "_dates.csv",labels,dateInfo,nodeNames,",");
WriteSeparatedTable (basePath + "_spans.csv",labels,dateInfo,nodeNames,",");

/* generate rate estimates */

for (k3 = 0; k3 < allCount; k3=k3+1)
{
	myDS = GatherDescriptiveStats (rateEstimates[-1][k3]);
	rateInfo[k3][1] = myDS["2.5%"];
	rateInfo[k3][2] = myDS["Mean"];
	rateInfo[k3][3] = myDS["Median"];
	rateInfo[k3][4] = myDS["97.5%"];
}

WriteSeparatedTable (basePath + "_rates.csv",labels,rateInfo,allNames,",");

fprintf (stdout, "[PHASE 5.] Making a PostScript tree with the model averaged rate estimates\n");


md = 0;

for (k3 = 1; k3 < Abs(treeAVL); k3=k3+1)
{
	md = Max (md, (treeAVL[k3])["Depth"]);
}

baseWidth 			= md * 100;

/*
bestModelDates		= internalDates[sortedScores[modelsRead-1][1]];
bestModelRates		= rateAssignments[sortedScores[modelsRead-1][1]];
bestModelRatesV		= rateValues[sortedScores[modelsRead-1][1]];
*/

bestModelDates		= dateInfo[-1][0];
bestModelRates		= rateAssignments[sortedScores[modelsRead-1][1]];
bestModelRatesV		= rateInfo[-1][0];

minR				= 1e100;
maxR				= 0;

iNodeCounter		= 0;

for (k=0; k<Columns(branchNames); k=k+1)
{
	nodeSpec  = {};
	nodeName  = branchNames [k];

	if (iNodes[nodeName])
	{
		nodeSpec ["TREE_OUTPUT_BRANCH_LABEL"] = "(" + Format(bestModelDates[iNodeCounter],6,3) + ") drawletter";
		iNodeCounter = iNodeCounter + 1;
	}
	if (k < allCount)
	{
		minR = Min(minR,bestModelRatesV[k]);
		maxR = Max(maxR,bestModelRatesV[k]);
	}
	TREE_OUTPUT_OPTIONS[nodeName] = nodeSpec;
}

colorByNode	= {};

for (k=0; k<allCount; k=k+1)
{
	myRV	  = bestModelRatesV[k];
	nodeName  = branchNames [k];
	myColor	  = {{(myRV-minR)/(maxR-minR)} /* 1 for maxR, 0 for minR; linear between */
				 {0}
				 {(maxR-myRV)/(maxR-minR)} /* 1 for minR, 0 for maxR; linear between */
				};
				
	colorByNode [nodeName] = myColor;
	(TREE_OUTPUT_OPTIONS[nodeName])["TREE_OUTPUT_BRANCH_COLOR"] = myColor;
	(TREE_OUTPUT_OPTIONS[nodeName])["TREE_OUTPUT_OVER_BRANCH"] = "0 0 0 setrgbcolor\n5 5 rmoveto\n("+Format(myRV/maxR*100,5,2)+"\\%) show";
	
}

psString 			= PSTreeString (givenTree,"",{{baseWidth,baseHeight}});


psLegend = "";
psLegend * 256;

sortedKeys = {{minR,(minR+maxR)/2,maxR}};

currentPosY = baseHeight + 10;
currentPosX = 25;

for (k=2; k>=0; k=k-1)
{
	aKey 		= sortedKeys[k];
	colorMx 	= branchColors[k];
	psLegend * (""+colorMx[0]+" "+colorMx[1]+" "+colorMx[2]+ " setrgbcolor\n"+currentPosX+" "+currentPosY+" moveto\n");
	psLegend * ("("+Format(aKey,10,5)+" subs/site/"+units+") show\n");
	currentPosY = currentPosY + 20;
}

newHeight = baseHeight + (bestModelRC+1)*20;

psLegend * "showpage";
psLegend * 0;
repMx = {{"showpage"}{psLegend}};
psString = psString ^ repMx;
repMx = {{"/PageSize\\ \\[[0-9\\ ]+\\]"}{"/PageSize ["+baseWidth+" "+ newHeight+ "]"}};
psString = psString ^ repMx;
repMx = {{"setfont"}{"setfont\n3 setlinewidth\n1 setlinecap"}};
psString = psString ^ repMx;


fName = basePath + "_ratesDates.ps";
fprintf (fName, CLEAR_FILE,drawLetter,psString);
							
							
fprintf (stdout, "[PHASE 6.] Making a model averaged PostScript tree with dated tips and internal nodes (relative to the root)\n");

distanceAVL			= {};
rootD				= bestModelDates[Rows(bestModelDates)-1];

k2 = 0;
k  = 0;
md = 0;

for (k3 = 1; k3 < Abs(treeAVL); k3=k3+1)
{
	nodeName = (treeAVL[k3])["Name"];
	if(Abs((treeAVL[k3])["Children"]))
	{
		(TREE_OUTPUT_OPTIONS[nodeName])["TREE_OUTPUT_BRANCH_LABEL"] = "(" + nodeName + ") drawletter" ;
	}
	else
	{
		(TREE_OUTPUT_OPTIONS[nodeName]) - "TREE_OUTPUT_BRANCH_LABEL";
	}
	distanceAVL [nodeName] = branchSpanInfo[nodeName];
	(TREE_OUTPUT_OPTIONS[nodeName]) ["TREE_OUTPUT_OVER_BRANCH"] = "0 0 0 setrgbcolor\n5 5 rmoveto\n("+
																	Format(distanceAVL [nodeName],4,1)+") show";	
	md = Max (md, (treeAVL[k3])["Depth"]);
}


ts = PostOrderAVL2StringDistances (treeAVL, distanceAVL);
Tree T = ts;

fprintf (tName, Format(T,0,1),"\n");

psString 		= PSTreeString (T,"STRING_SUPPLIED_LENGTHS",{{baseWidth,baseHeight}});
repMx = {{"showpage"}{psLegend}};
psString = psString ^ repMx;
repMx = {{"/PageSize\\ \\[[0-9\\ ]+\\]"}{"/PageSize [ "+baseWidth+" "+ newHeight+ "]"}};
psString = psString ^ repMx;
repMx = {{"setfont"}{"setfont\n3 setlinewidth\n1 setlinecap"}};
psString = psString ^ repMx;
fName = basePath + "_scaled.ps";
fprintf (fName, CLEAR_FILE,drawLetter,psString);

		
