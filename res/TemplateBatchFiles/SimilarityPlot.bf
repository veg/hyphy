ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" 
								   + DIRECTORY_SEPARATOR + "GrabBag.bf"); 

ExecuteAFile (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" 
								   + DIRECTORY_SEPARATOR + "NJ.bf"); 

SetDialogPrompt ("Please select a data set to process:");
DataSet			 ds  		  = ReadDataFile	(PROMPT_FOR_FILE);
DataSetFilter	 filteredData = CreateFilter	(ds,1);

fprintf (stdout, "[READ ", filteredData.species, " SEQUENCES AND ", filteredData.sites, " SITES]\n");

windowSize  	= prompt_for_a_value ("Window size", 100, Min(50,ds.sites), ds.sites, 1);
windowStride	= prompt_for_a_value ("Window stride", windowSize$5, 1, windowSize, 1);

ChoiceList	 	(whichSequence, "Which Sequence Is The Query?", 1, SKIP_NONE, ds);

if (whichSequence < 0)
{
	return 0;
}

GetString	 			(seqNames, ds, -1);
DISTANCE_PROMPTS 		= 1;
ExecuteAFile 			(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "chooseDistanceFormula.def");

querySequenceName 	= seqNames[whichSequence] && 1;

fprintf		 		(stdout, "\n[SCREENING ", querySequenceName, " USING WINDOW WIDTH ",  windowSize, " AND WINDOW STRIDE ", windowStride, "]\n");
fprintf		 		(stdout, "[DISTANCE METRICS READ FROM ", whichFile, "]\n");


ChoiceList	 		(doNPBS, "Include NJ bootstrap", 1, SKIP_NONE, "No", "Only plot pairwise genetic distances over sliding windows","Yes","For each window, determine if there is sufficient phylogenetic evidence to cluster the query with one of the reference sequences");
if (doNPBS < 0)
{
	return 0;
} 

if (doNPBS == 1)
{
	npbsCutoff  	= prompt_for_a_value ("Bootstrap support", 70, 0, 100, 1);
}

InitializeDistances (0);

steps 	= (ds.sites - windowSize)$windowStride+1;

gdInfo  	  = {steps, ds.species+2};
columnHeaders = {1, ds.species+2};
columnHeaders [0] = "Left";
columnHeaders [1] = "Right";
columnHeaders [2] = "Midpoint";
columnPlots = "";

columnIdx  		 = 3;

for (seqID = 0; seqID < ds.species; seqID = seqID + 1)
{
	if (seqID != whichSequence)
	{
		columnHeaders[columnIdx] = seqNames[seqID];
		columnIdx = columnIdx + 1;
		columnPlots = columnPlots + seqNames[seqID] + ";";
	}
}

DO_NOT_CALL_INITIALIZE_DISTANCES = 1;

if (doNPBS)
{
	fprintf (stdout, "Clustering and bootstrap support on each window interval\n");
}

for 	  (aStep = 0; aStep < steps; aStep = aStep+1)
{
	leftEdge = aStep * windowStride;
	rightEdge = leftEdge + windowSize;
	
	DataSetFilter	 filteredData = CreateFilter	(ds,1,siteIndex >= leftEdge && siteIndex < rightEdge);
	columnIdx  		 = 3;
	gdInfo[aStep][0] = leftEdge;
	gdInfo[aStep][1] = rightEdge;
	gdInfo[aStep][2] = 0.5(leftEdge+rightEdge);
	
	for (seqID = 0; seqID < ds.species; seqID = seqID + 1)
	{
		if (seqID != whichSequence)
		{
			gdInfo[aStep][columnIdx] = ComputeDistanceFormula (seqID, whichSequence);
			columnIdx = columnIdx + 1;
		}
	}
	
	if (doNPBS)
	{
		siblingName 	= getNJNeighbor (InferTreeTopology(1));
		bsSupport		= 0;
		
		for (repID = 0; repID < 100; repID = repID + 1)
		{
			DataSetFilter	 	filteredData = Bootstrap	(ds,1,siteIndex >= leftEdge && siteIndex < rightEdge);
			simSN 				= getNJNeighbor (InferTreeTopology(1));
			bsSupport			= bsSupport + (simSN == siblingName);
		}
		intervalSpec = "[" + leftEdge + "-" + rightEdge + "]";
		if (bsSupport >= npbsCutoff)
		{
			fprintf (stdout, intervalSpec, " : ", siblingName, "\n");
		}
		else
		{
			fprintf (stdout, intervalSpec, " : Ambiguous \n");		
		}
	}
	
}

OpenWindow (CHARTWINDOW,{{"Similarity Plot for " + seqNames[whichSequence]}
						{"columnHeaders"}
						{"gdInfo"}
						{"Scatterplot"}
						{"Midpoint"}
						{columnPlots}
						{"Midpoint, bp"}
						{""}
						{"Genetic distance"}
						{"3"}
						{""}
						{"-1;-1"}
						{"10;1.309;0.785398"}
						{"Times:12:0;Times:10:0;Times:12:2"}
						{"0;0;16777215;1644825;0;0;6579300;11842740;13158600;14474460;0;3947580;16777215;15670812;6845928;16771158;2984993;9199669;7018159;1460610;16748822;11184810;14173291"}
						{"16,0,0"}
						},
						"(SCREEN_WIDTH-60)/2;(SCREEN_HEIGHT-50)/2;100;100");
						
				
/*----------------------------------------------------------------------------*/

function getNJNeighbor (treeString)
{
	Tree njTree = treeString;
	njTreeAVL   = njTree ^ 0;
	
	/* find the query node */
	
	for (k = 1; k < Abs(njTreeAVL); k=k+1)
	{
		nodeName = (njTreeAVL[k])["Name"] && 1;
		if (nodeName == querySequenceName)
		{
			break;
		}
	}
	p 		= (njTreeAVL[k])["Parent"];
	pn 		= njTreeAVL[p];
	if (((pn["Children"])[0]) == k)
	{
		sn = (pn["Children"])[1];
	}
	else
	{
		sn = (pn["Children"])[0];
	}
	pn = njTreeAVL[sn];
	return pn["Name"];
}
