function 	_computeSNSSites (_filterName, _genCode, _codonFreqs, callCount)
{
	if (callCount == 0)
	{
		ExecuteCommands ("GetDataInfo (charInfo, `_filterName`, \"CHARACTERS\");_codonCount=`_filterName`.sites;");
		nonStopCount = Columns (charInfo);
		/* make syn and non-syn template matrices */
		ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Distances"+DIRECTORY_SEPARATOR+"CodonTools.def");
	}

	sSites  = 0;
	nsSites = 0;

	synM    = {nonStopCount,nonStopCount};
	nonSynM = {nonStopCount,nonStopCount};

	vertOnes = {nonStopCount,1}["1"];
	horOnes  = {1,nonStopCount}["1"];

	hShift = 0;
	for (h1 = 0; h1 < 64; h1=h1+1)
	{
		gc1 = _genCode[h1];
		if (gc1 == 10)
		{
			hShift = hShift+1;
		}
		else
		{
			sSites  = sSites   + _codonCount * _S_NS_POSITIONS_[0][h1] * _codonFreqs[h1-hShift];
			nsSites = nsSites + _codonCount * _S_NS_POSITIONS_[1][h1] * _codonFreqs[h1-hShift];
			
			vShift = hShift;
			for (v1 = h1+1; v1 < 64; v1=v1+1)
			{
				gc2 = _genCode[v1];
				if (gc2 == 10)
				{
					vShift = vShift + 1;
				}
				else
				{
					if (gc1 == gc2)
					{
						synM [h1-hShift][v1-vShift] = _codonFreqs[h1-hShift];
						synM [v1-vShift][h1-hShift] = _codonFreqs[v1-vShift];
					}
					else
					{
						nonSynM [h1-hShift][v1-vShift] = _codonFreqs[h1-hShift];
						nonSynM [v1-vShift][h1-hShift] = _codonFreqs[v1-vShift];
					}
				}
			}
		}
	}
	
	_returnAVL = {};
	_returnAVL ["Dimension"] = nonStopCount;
	_returnAVL ["Sites"] 	 = _codonCount*3;
	_returnAVL ["NSSites"] 	 = nsSites;
	_returnAVL ["SSites"] 	 = sSites;	
	_returnAVL ["synM"] 	 = synM;	
	_returnAVL ["nonSynM"] 	 = nonSynM;
	
	return 		_returnAVL;
}
