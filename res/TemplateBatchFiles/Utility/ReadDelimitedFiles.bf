/*----------------------------------------------------------------*/

function ReadCSVTable (fileName, haveHeader)
{
	if (Abs(fileName) == 0) {
		fscanf (PROMPT_FOR_FILE, REWIND, "Lines", inData);
	}
	else {
		fscanf (fileName, REWIND, "Lines", inData);		
	}
	if (haveHeader) {
		output = {};
		output[0] = splitOnRegExp (inData[0],"\\,");
	}
	felMXString = "";
	felMXString * 256;
	felMXString * "{";
	if (haveHeader > 1)
	{
		rowHeaders = {Columns(inData) - 1, 1};
		for (lineID = 1; lineID < Columns(inData); lineID = lineID + 1)
		{
			_tempMatrix = inData[lineID]$"\\,";
			if (_tempMatrix[0]>=0)
			{
				rowHeaders [lineID-1] = (inData[lineID])[0][_tempMatrix[0]-1];
				felMXString * ("{" + (inData[lineID])[_tempMatrix[0]+1][Abs(inData[lineID])-1] + "}\n");
			}
		}		
		output [2] = rowHeaders;
	}
	else
	{
		for (lineID = haveHeader; lineID < Columns(inData); lineID += 1) {
			felMXString * ("{" + inData[lineID] + "}\n");
		}
	}
	felMXString * "}";
	felMXString * 0;
	felMXString = Eval (felMXString);
	inData = 0;
	if (haveHeader)
	{
		output[1] = felMXString;
		felMXString = 0;
		return output;
	}
	return felMXString;
}

/*----------------------------------------------------------------*/

lfunction ReadCSVTableText (fileName, haveHeader)
{
	if (Abs(fileName) == 0) {
		fscanf (PROMPT_FOR_FILE, REWIND, "Lines", inData);
	}
	else {
		fscanf (fileName, REWIND, "Lines", inData);		
	}
	if (haveHeader) {
		output = {};
		output[0] = splitOnRegExp (inData[0],"\\,");
	}
	
	_theMatrix = {};
	
	if (haveHeader > 1)
	{
		rowHeaders = {Columns(inData) - 1, 1};
		for (lineID = 1; lineID < Columns(inData); lineID += 1)
		{
			_tempMatrix = inData[lineID]$"\\,";
			if (_tempMatrix[0]>=0)
			{
				rowHeaders [lineID-1] = (inData[lineID])[0][_tempMatrix[0]-1];
				felMXString * ("{" + (inData[lineID])[_tempMatrix[0]+1][Abs(inData[lineID])-1] + "}\n");
			}
		}		
		output [2] = rowHeaders;
	}
	else
	{
		for (lineID = haveHeader; lineID < Columns(inData); lineID += 1) {
			values = splitOnRegExp (inData[lineID], "\\,");
			_theMatrix + {};
			for (columnID = 0; columnID < Abs (values); columnID += 1) {
			    theEntry = values[columnID];
			    isQuoted = Abs(theEntry)>=2 && theEntry[0]=="\"" && theEntry[Abs(theEntry)-1] == "\"";
			    if (isQuoted) {
			        if (Abs(theEntry) == 2) {
			            (_theMatrix[Abs(_theMatrix)-1]) + "";
			        }
			        else {
			            (_theMatrix[Abs(_theMatrix)-1]) + theEntry[1][Abs(theEntry)-2];
			        }
			    } else {
			        (_theMatrix[Abs(_theMatrix)-1]) + theEntry;			    
			    }
			}
		}
	}
	if (haveHeader) {
		output[1] = _theMatrix;
		_theMatrix = 0;
		return output;
	}
	return _theMatrix;
}



/*----------------------------------------------------------------*/

function ReadMatchRegExp (fileName, regExp)
{
	if (Abs(fileName) == 0)
	{
		fscanf (PROMPT_FOR_FILE, REWIND, "Lines", inData);
	}
	else
	{
		fscanf (fileName, REWIND, "Lines", inData);		
	}
	
	linesRead = Columns(inData);
	_tempMatrix = {linesRead,1};
	
	for (lineID = 0; lineID < linesRead; lineID = lineID + 1)
	{
		aLine   = inData[lineID];
		regExpM = aLine $ regExp;
		if (regExpM [0] >= 0)
		{
			sidx = regExpM[0];
			eidx = regExpM[1];
			_tempMatrix[lineID] = 0+aLine[sidx][eidx];
		}
	}
	return _tempMatrix;
}

/*----------------------------------------------------------------*/

function ReadSplitOnRegExp (fileName, regExp)
{
	return ReadSplitOnRegExpCallback (fileName, regExp, "");	
}

/*----------------------------------------------------------------*/

function ReadSplitOnRegExpCallback (fileName, regExp,func)
{
	if (Abs(fileName) == 0)
	{
		fscanf (PROMPT_FOR_FILE, REWIND, "Lines", inData);
	}
	else
	{
		fscanf (fileName, REWIND, "Lines", inData);		
	}
	
	linesRead = Columns(inData);
	if (Abs(func) == 0)
	{
		_tempMatrix = {};
	}
	
	for (lineID = 0; lineID < linesRead; lineID = lineID + 1)
	{
		if (Abs(func))
		{
            thisLine = inData[lineID];
			ExecuteCommands ("_lineStatus = " + func + "(splitOnRegExp (thisLine, regExp),lineID)"); 
			if (_lineStatus < 0)
			{
				return 1;
			}
		}
		else
		{
			_tempMatrix[lineID] = splitOnRegExp (inData[lineID], regExp);
		}
	}
	if (Abs(func) == 0)
	{
		return _tempMatrix;
	}
	return 0;
}

/*----------------------------------------------------------------*/

function ReadTabTable (fileName, haveHeader)
{
	if (Abs(fileName) == 0)
	{
		fscanf (PROMPT_FOR_FILE, REWIND, "Lines", inData);
	}
	else
	{
		fscanf (fileName, REWIND, "Lines", inData);		
	}
	if (haveHeader)
	{
		output = {};
		output[0] = splitStringByTab (inData[0]);
	}
	felMXString = "";
	felMXString * 256;
	felMXString * "_tempMatrix={";
	for (lineID = haveHeader; lineID < Columns(inData); lineID = lineID + 1)
	{
		felMXString * ("{" + inData[lineID]^{{"[^0-9eE\.\+\-]+",","}} + "}\n");
	}
	felMXString * "}";
	felMXString * 0;
	ExecuteCommands (felMXString);
	felMXString = 0;
	inData = 0;
	if (haveHeader)
	{
		output[1] = _tempMatrix;
		_tempMatrix = 0;
		return output;
	}
	return _tempMatrix;
}

/*----------------------------------------------------------------*/

function splitOnRegExp (string, splitter)
{
	matched = string || splitter;
	splitBits = {};
	if (matched [0] < 0)
	{
		splitBits[0] = string;
	}
	else
	{
		mc = 0;
		if (matched[0] == 0)
		{
			fromPos = matched[1]+1;
			mc = 2;
		}
		else
		{
			fromPos = 0;
			toPos	= 0;
		}
		for (; mc < Rows (matched); mc = mc+2)
		{
			toPos = matched[mc]-1;
			splitBits [Abs(splitBits)] = string[fromPos][toPos];
			fromPos    = matched[mc+1]+1;
		}
		splitBits [Abs(splitBits)] = string[fromPos][Abs(string)-1];
	}
	return splitBits;
}


/*----------------------------------------------------------------*/

function extractSubexpressions (string, splitter, merge, spacer)
{
	matched = string $ splitter;
	if (merge)
	{
		splitBits = "";
		splitBits * 128;
	}
	else
	{
		splitBits = {};
	}
	rMatched = Rows(matched);
	if (rMatched>2)
	{
		if (merge)
		{
			for (mc = 2; mc < rMatched; mc = mc+2)
			{
				if (mc>2)
				{
					splitBits * spacer;
				}
				splitBits * string[matched[mc]][matched[mc+1]];
			}
		}		
		else
		{
			for (mc = 2; mc < rMatched; mc = mc+2)
			{
				splitBits [Abs(splitBits)] = string[matched[mc]][matched[mc+1]];
			}
		}
	}
	if (merge)
	{
		splitBits * 0;		
	}
	return splitBits;
}

/*----------------------------------------------------------------*/

function findAllUniqueExpressions (list, regExp, isMatrix)
{
	uniqueExpressions = {};
	if (isMatrix)
	{
		totalCount = Rows(list)*Columns(list);
	}
	else
	{
		totalCount = Abs (list);
	}
	
	for (_k = 0; _k < totalCount; _k += 1)
	{
		matchRes = list[_k] $ regExp;
		if (matchRes[0] >= 0)
		{
			matchRes = (list[_k])[matchRes[0]][matchRes[1]];
			uniqueExpressions [matchRes] += 1;
		}
	}
	
	return uniqueExpressions;
}


/*----------------------------------------------------------------*/

function extractAllExpressions (_string, _splitter, _subexpression)
{
	_matched = _string || _splitter;
	_splitBits = {};
	if (_matched [0] >= 0)
	{
		if (Abs(_subexpression))
		{
			for (_mc = 0; _mc < Rows (_matched); _mc = _mc+2)
			{
				_splitBits [Abs(_splitBits)] = extractSubexpressions(_string[_matched[_mc]][_matched[_mc+1]],_subexpression,1,"");
			}
		}
		else
		{
			for (_mc = 0; _mc < Rows (_matched); _mc = _mc+2)
			{
				_splitBits [Abs(_splitBits)] = _string[_matched[_mc]][_matched[_mc+1]];
			}		
		}
	}
	return _splitBits;
}

/*----------------------------------------------------------------*/

function splitStringByTab (inString)
{
	bits      = {};
	lastIndex = 0;
	strlen 	  = Abs(inString);
	
	for (currentIndex = 0; currentIndex < strlen; currentIndex = currentIndex + 1)
	{
		cchar = inString[currentIndex];
		if (cchar == "\t")
		{
			if (currentIndex == lastIndex)
			{
				bits[Abs(bits)] = "";
			}
			else
			{
				bits[Abs(bits)] = inString[lastIndex][currentIndex-1];
				lastIndex = currentIndex + 1;
			}
		}
	}
	if (currentIndex>lastIndex)
	{
		bits[Abs(bits)] = inString[lastIndex][currentIndex-1];
	}
	
	return bits;
}

/*----------------------------------------------------------------*/

function normalizeSequenceID (seqName, alreadyDefined&)
{
	_testSeqName 	= seqName ^ {{"[^a-zA-Z0-9\\_]","_"}};
	_testSeqName 	= _testSeqName ^ {{"_+","_"}};
	_testSeqName2	= _testSeqName;
	
	for (currentIndex = 2; ;currentIndex=currentIndex+1)
	{	
		_z = alreadyDefined[_testSeqName2];
		if (_z==0)
		{
			alreadyDefined[_testSeqName2] = 1;
			break;
		}
		_testSeqName2 = _testSeqName + "_" + currentIndex;
	}
	return _testSeqName2;
}

/*----------------------------------------------------------------*/

function normalizeSequenceNamesInAFilter (filterName)
{
	return normalizeSequenceNamesInAFilterWithTree (filterName, "");
}

/*----------------------------------------------------------------*/

function normalizeSequenceNamesInAFilterWithTree (filterName, treeName)
{
    if (Abs (treeName) > 0) {
        LoadFunctionLibrary ("TreeTools");
        ExecuteCommands ("_treeAVL = `treeName`^0");
        _node_id_by_name = {};
        for (_nodeID = 1; _nodeID < Abs(_treeAVL); _nodeID += 1) {
            _node_id_by_name [(_treeAVL[_nodeID])["Name"]] = _nodeID;
        }
    }
    
 	ExecuteCommands ("GetString(_seqNames,`filterName`,-1)");
	
	_renamingBuffer = {};
	
	_totalRenamed = 0;
	for (_k = 0; _k < Columns (_seqNames); _k+=1)
	{
		_seqName2 = normalizeSequenceID(_seqNames[_k],"_renamingBuffer");
		if (_seqName2 != _seqNames[_k])
		{
			ExecuteCommands ("SetParameter (`filterName`,_k,_seqName2)");
			if (Abs (treeName) > 0) {
			    (_treeAVL[_node_id_by_name [_seqNames[_k]]])["Name"] = _seqName2;
			}

			_totalRenamed += 1;
		}
	}
	
	if (_totalRenamed && Abs (treeName) > 0) {
	    ExecuteCommands ("Topology `treeName` = " + PostOrderAVL2String (_treeAVL));
	}
	
	return _totalRenamed;
}

/*----------------------------------------------------------------*/

function matchStringToSetOfPatterns (str, patterns)
{
	for (currentIndex = 0; currentIndex < Abs(patterns);currentIndex=currentIndex+1)
	{
		if ((str$patterns[currentIndex])[0] >= 0)
		{
			return currentIndex;
		}
	}
	return -1;
}

/*----------------------------------------------------------------*/

function splitFilterBasedOnRegExp (filterName, patterns)
{
	ExecuteCommands ("GetString(_seqNames,`filterName`,-1)");
	
	sequenceByRegExp = {};
	
	for (_k = 0; _k < Abs(patterns); _k += 1)
	{
		aPattern = patterns["INDEXORDER"][_k];
		subfilterName = patterns[aPattern];
		ExecuteCommands ("DataSetFilter `subfilterName` = CreateFilter (`filterName`,1,\"\",(_seqNames[speciesIndex]$\"`aPattern`\")[0]>=0)");
		sequenceByRegExp[subfilterName] = Eval ("`subfilterName`.species");
	}
	return sequenceByRegExp;
}

/*----------------------------------------------------------------*/

function condenseFilterIntoUniqueSequences (filterName, returnName, matchMode)
{
	ExecuteCommands ("GetDataInfo(_uniquePatterns,`filterName`," + matchMode + ")");
	ExecuteCommands ("GetString(_seqNames,`filterName`,-1)");

	_filterMeIn = {};
	for (_k = 0; _k < Columns (_uniquePatterns["UNIQUE_INDICES"]); _k+=1)
	{
		_filterMeIn[(_uniquePatterns["UNIQUE_INDICES"])[_k]] = (_uniquePatterns["UNIQUE_COUNTS"])[_k];
		ExecuteCommands ("SetParameter(`filterName`,
				(_uniquePatterns[\"UNIQUE_INDICES\"])[_k],
				_seqNames[(_uniquePatterns[\"UNIQUE_INDICES\"])[_k]] + \"_\" + (_uniquePatterns[\"UNIQUE_COUNTS\"])[_k]);");
	}
	ExecuteCommands ("DataSetFilter `returnName` = CreateFilter(`filterName`,1,,Abs(_filterMeIn[speciesIndex])>0)");
	
	return 0;
}
