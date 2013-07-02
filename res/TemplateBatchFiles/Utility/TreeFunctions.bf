/*---------------------------------------------------------*/

function computeTotalDivergence (treeID)
{
	bNames = BranchName   (^treeID,-1);
	bLen   = BranchLength (^treeID,-1);
	
	sum  = 0;
	sum2 = +bLen;
	
	for (_k=0; _k<Columns(bNames); _k += 1) {
		sum += bLen[_k]*multFactors[bNames[_k]];
	}	
	return {{sum__,sum2__}};
}

/*---------------------------------------------------------*/

function computeMultFactors (treeID)
{
	treeAVL2 = (^treeID)^ 0;
	leafCount=Max(2,TipCount(^treeID)); 
	
	multFactors = {};
	
	for (_k=1; _k<Abs(treeAVL2); _k += 1)
	{
		aNode			= treeAVL2[_k];
		aNodeName		= aNode["Name"];
		parentIndex		= aNode["Parent"];
		_k2				= Abs(aNode["Children"]);
		if (_k2)
		{
			currentDepth		   = aNode["Below"];
			multFactors[aNodeName] = currentDepth;		
			if (parentIndex > 0)
			{
				(treeAVL2[parentIndex])["Below"] += currentDepth;
			}
		}
		else
		{
			multFactors[aNodeName]		= 1;
			(treeAVL2[parentIndex])["Below"] += 1;
		}
		
	}

	pKeys 			= Rows(multFactors);
	for (_k=0; _k<Columns(pKeys); _k=_k+1)
	{
		aNodeName = pKeys[_k];
		multFactors[aNodeName] = multFactors[aNodeName] * (leafCount-multFactors[aNodeName]);
	}
	
	return 0;
}

