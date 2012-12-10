/* ________________________________________________________________________________________________*/


function InitializeDistances ()
{
	HarvestFrequencies (_dNucFreq,filteredData,1,1,0);
	_d_fR = _dNucFreq[0]+_dNucFreq[2];
	_d_fY = _dNucFreq[1]+_dNucFreq[3];
	
	if (_dNucFreq[0] == 0 || _dNucFreq[1] == 0 || _dNucFreq[2] == 0 || _dNucFreq[3] == 0)
	{
		_useK2P = 1;
	}
	else
	{
		_d_TN_K1 = 2*_dNucFreq[0]*_dNucFreq[2]/_d_fR;
		_d_TN_K2 = 2*_dNucFreq[1]*_dNucFreq[3]/_d_fY;
		_d_TN_K3 = 2*(_d_fR*_d_fY-_dNucFreq[0]*_dNucFreq[2]*_d_fY/_d_fR-_dNucFreq[1]*_dNucFreq[3]*_d_fR/_d_fY);
		_useK2P = 0;
	}
	
	

	return 0;
}

/* ________________________________________________________________________________________________*/

function ComputeDistanceFormulaFromDiffMx (siteDifferenceCount)
{
	totalSitesCompared = +siteDifferenceCount;
    if (_useK2P)
    {
        _dTransitionCounts 	 =    siteDifferenceCount[0][2]+siteDifferenceCount[2][0]  /* A-G and G-A */
                                 +siteDifferenceCount[1][3]+siteDifferenceCount[3][1]; /* C-T and T-C */
                            
        _dTransversionCounts = (siteDifferenceCount[0][0]+siteDifferenceCount[1][1]+siteDifferenceCount[2][2]+siteDifferenceCount[3][3])+_dTransitionCounts;
        
        _dTransitionCounts	 = _dTransitionCounts/totalSitesCompared;
        _dTransversionCounts = 1-_dTransversionCounts/totalSitesCompared;
        
        _d1C = 1-2*_dTransitionCounts-_dTransversionCounts;
        _d2C = 1-2*_dTransversionCounts;
        
        if (_d1C>0 && _d2C>0)
        {
            return -(0.5*Log(_d1C)+.25*Log(_d2C));	
        }
    }
    else
    {
        _dAGCounts 	 =    siteDifferenceCount[0][2]+siteDifferenceCount[2][0]  /* A-G and G-A */;
        _dCTCounts	 = 	  siteDifferenceCount[1][3]+siteDifferenceCount[3][1]; /* C-T and T-C */
                            
        _dTransversionCounts = (siteDifferenceCount[0][0]+siteDifferenceCount[1][1]+siteDifferenceCount[2][2]+
                                siteDifferenceCount[3][3])+_dAGCounts+_dCTCounts;
        
        _dAGCounts	 = _dAGCounts/totalSitesCompared;
        _dCTCounts	 = _dCTCounts/totalSitesCompared;
        
        _dTransversionCounts = 1-_dTransversionCounts/totalSitesCompared;
        
        _d1C = 1-_dAGCounts/_d_TN_K1-0.5*_dTransversionCounts/_d_fR;
        _d2C = 1-_dCTCounts/_d_TN_K2-0.5*_dTransversionCounts/_d_fY;
        _d3C = 1-0.5*_dTransversionCounts/_d_fY/_d_fR;
        
        if ((_d1C>0)&&(_d2C>0)&&(_d3C>0))
        {
            return -_d_TN_K1*Log(_d1C)-_d_TN_K2*Log(_d2C)-_d_TN_K3*Log(_d3C);
        }
    }
	return 1000;
}

/* ________________________________________________________________________________________________*/

function ComputeDistanceFormula (s1,s2)
{
	GetDataInfo (siteDifferenceCount, filteredData, s1, s2, DIST);
	return ComputeDistanceFormulaFromDiffMx(siteDifferenceCount);
}

/* ________________________________________________________________________________________________*/

function TreeMatrix2TreeStringGen (treeNodes,leafCount,leafNames,doLengths)
{
	treeString = "";
	p = 0;
	k = 0;
	m = treeNodes[0][1];
	n = treeNodes[0][0];
	treeString*(Rows(treeNodes)*25);

	while (m)
	{	
		if (m>p)
		{
			if (p)
			{
				treeString*",";
			}
			for (j=p;j<m;j=j+1)
			{
				treeString*"(";
			}
		}
		else
		{
			if (m<p)
			{
				for (j=m;j<p;j=j+1)
				{
					treeString*")";
				}
			}	
			else
			{
				treeString*",";
			}	
		}
		if (n<leafCount)
		{
			if (Abs(leafNames))
			{
				nodeName = ""+leafNames[n];
			}
			else
			{
				GetString (nodeName, filteredData, n);
			}
			treeString*nodeName;
		}
		
		if (doLengths>.5)
		{
			nodeName = ":"+treeNodes[k][2];
			treeString*nodeName;
		}
		
		k=k+1;
		p=m;
		n=treeNodes[k][0];
		m=treeNodes[k][1];
	}

	for (j=m;j<p;j=j+1)
	{
		treeString*")";
	}
	
	treeString*0;
	return treeString;
}

/* ________________________________________________________________________________________________*/

function TreeMatrix2TreeString (doLengths)
{
	return TreeMatrix2TreeStringGen (treeNodes, filteredData.species, 0, doLengths);
}

/* ________________________________________________________________________________________________*/

function InferTreeTopologyFromMatrixGen (distanceMatrix, leafNames, distancesFlag)
{
	MESSAGE_LOGGING 		 	= 1;
	cladesMade 					= 1;
	
	leafCount= Abs(leafNames);
	if (leafCount == 0)
	{
		leafCount = filteredData.species;
	}

	if (leafCount== 2)
	{
		d1 = distanceMatrix[0][1]/2;
		treeNodes = {{0,1,d1__},
					 {1,1,d1__},
					 {2,0,0}};
					 
		cladesInfo = {{2,0}};
	}
	else
	{
		if (leafCount == 3)
		{
			d1 = (distanceMatrix[0][1]+distanceMatrix[0][2]-distanceMatrix[1][2])/2;
			d2 = (distanceMatrix[0][1]-distanceMatrix[0][2]+distanceMatrix[1][2])/2;
			d3 = (distanceMatrix[1][2]+distanceMatrix[0][2]-distanceMatrix[0][1])/2;
			treeNodes = {{0,1,d1__},
						 {1,1,d2__},
						 {2,1,d3__}
						 {3,0,0}};
						 
			cladesInfo = {{3,0}};		
		}
		else
		{	
			njm = (distanceMatrix > methodIndex)>=leafCount;
				
			treeNodes 		= {2*(leafCount+1),3};
			cladesInfo	    = {leafCount-1,2};
			
			for (i=Rows(treeNodes)-1; i>=0; i=i-1)
			{
				treeNodes[i][0] = njm[i][0];
				treeNodes[i][1] = njm[i][1];
				treeNodes[i][2] = njm[i][2];
			}

			for (i=Rows(cladesInfo)-1; i>=0; i=i-1)
			{
				cladesInfo[i][0] = njm[i][3];
				cladesInfo[i][1] = njm[i][4];
			}
			
			njm = 0;
		}
	}
	
	distanceMatrix = 0;
	
	return TreeMatrix2TreeStringGen (treeNodes, leafCount, leafNames, distancesFlag);
}


/* ________________________________________________________________________________________________*/

function InferTreeTopologyFromMatrix (distancesFlag)
{
	return InferTreeTopologyFromMatrixGen (distanceMatrix, 0, distancesFlag);
}



/* ________________________________________________________________________________________________*/

function InferTreeTopology(distancesFlag)
{
	InitializeDistances ();
	distanceMatrix = {filteredData.species,filteredData.species};
		
	for (i = 0; i<filteredData.species; i=i+1)
	{
		for (j = i+1; j<filteredData.species; j = j+1)
		{
			distanceMatrix[i][j] = ComputeDistanceFormula (i,j);
		}
	}

	return InferTreeTopologyFromMatrix (distancesFlag);
}
