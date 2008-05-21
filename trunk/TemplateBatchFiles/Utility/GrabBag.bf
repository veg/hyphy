/*---------------------------------------------------------*/
/* Turn the keys of an AVL into a string for labeling 
   chart rows */
   
function avlToLabels (_gb_anAVL,_gb_prefix,_gb_delim)
{
	_gb_resString = "";
	_gb_keys	  = Rows (_gb_anAVL);
	_gb_count	  = Columns (_gb_keys);
	_gb_resString * 128;
	_gb_resString * _gb_prefix;
	if (Abs(_gb_prefix))
	{
		_gb_resString * _gb_delim;
	}
	if (_gb_count)
	{
		_gb_resString * _gb_keys[0];
	}
	for (_gb_idx = 1; _gb_idx < _gb_count; _gb_idx = _gb_idx + 1)
	{
		_gb_resString * (_gb_delim+_gb_keys[_gb_idx]);
	}
	_gb_resString * 0;
	return _gb_resString;
}

/*---------------------------------------------------------*/
/* Turn the keys of an AVL into a numerical column matrix  */
   
function avlKeysToMatrix (_gb_anAVL)
{
	_gb_keys	  = Rows (_gb_anAVL);
	_gb_count	  = Columns (_gb_keys);
	_gb_resMatrix = {_gb_count,1};

	for (_gb_idx = 0; _gb_idx < _gb_count; _gb_idx = _gb_idx + 1)
	{
		_gb_resMatrix[_gb_idx] = 0+_gb_keys[_gb_idx];
	}
	return _gb_resMatrix;
}

/*---------------------------------------------------------*/
/* Assuming that the AVL is 0..N indexed, produce a 
string with AVL entries separated by _gb_delim */

function avlToString (_gb_anAVL,_gb_delim)
{
	_gb_count	  = Abs (_gb_anAVL);
	_gb_resString = "";
	_gb_resString * 128;
	if (_gb_count)
	{
		_gb_resString * (""+_gb_anAVL[0]);
	}
	for (_gb_idx = 1; _gb_idx < _gb_count; _gb_idx = _gb_idx + 1)
	{
		_gb_resString * (_gb_delim+_gb_anAVL[_gb_idx]);
	}
	_gb_resString * 0;
	return _gb_resString;
}


/*---------------------------------------------------------*/
/* Stratify avl by values; store indices for each unique value */

function stratifyAVLByValues (_gb_anAVL)
{
	_gb_count	  = Abs (_gb_anAVL);
	_gb_resAVL	  = {};
	for (_gb_idx = 0; _gb_idx < _gb_count; _gb_idx = _gb_idx + 1)
	{
		_gb_key = _gb_anAVL[_gb_idx];
		if (Abs (_gb_resAVL[_gb_key]) == 0)
		{
			_gb_resAVL[_gb_key] = {};
		}
		(_gb_resAVL[_gb_key])[Abs(_gb_resAVL[_gb_key])] = _gb_idx;
	}
	return _gb_resAVL;
}

/*---------------------------------------------------------*/
/* Assuming that the AVL is 0..N indexed, produce a 
row matrix with AVL entries, using _gb_map to map the values 
and _gb_stride to do the conversion */

function avlToRow (_gb_anAVL, _gb_map, _gb_stride)
{
	_gb_count	  = Abs (_gb_anAVL);
	_gb_matrix	  = {1,_gb_count*_gb_stride};
	
	if (_gb_stride>1)
	{
		for (_gb_idx = 0; _gb_idx < _gb_count; _gb_idx = _gb_idx + 1)
		{
			for (_gb_idx2 = 0; _gb_idx2 < _gb_stride; _gb_idx2 = _gb_idx2 + 1)
			{
				_gb_matrix [_gb_idx*_gb_stride+_gb_idx2] = _gb_map[_gb_stride*_gb_anAVL[_gb_idx]+_gb_idx2];
			}
		}
	}	
	else
	{
		for (_gb_idx = 0; _gb_idx < _gb_count; _gb_idx = _gb_idx + 1)
		{
			_gb_matrix [_gb_idx] = _gb_map[_gb_anAVL[_gb_idx]];
		}
	}
	return _gb_matrix;
}

/*---------------------------------------------------------*/
/* Turn the keys of an AVL into a string for labeling 
   chart rows */
   
function splitFilePath (_filePath)
{
	_splitPath = {};
	_split     = _filePath $ ("[^\\"+DIRECTORY_SEPARATOR+"]+$");
	if (_split[0] == 0 && _split[1] == Abs (_filePath)-1) /* no path, all file name */
	{
		_splitPath ["DIRECTORY"] = "";
	}
	else
	{
		_splitPath ["DIRECTORY"] = _filePath[0][_split[0]-1];
		_filePath = _filePath[_split[0]][Abs(_filePath)];
	}

	_split     = _filePath || "\\.";
	if (_split[0] < 0) /* no extension */
	{
		_splitPath ["EXTENSION"] = "";
		_splitPath ["FILENAME"]  = _filePath;
 	}
	else
	{
		_splitPath ["EXTENSION"] = _filePath[_split[Rows(_split)-1]+1][Abs(_filePath)-1];
		_splitPath ["FILENAME"]  = _filePath[0][_split[Rows(_split)-1]-1];
	}
	return _splitPath;
}

/*---------------------------------------------------------*/
/* fix global variables in a LF at their current values */
   
function fixGlobalParameters (_lfName)
{
	ExecuteCommands ("GetString (_lfInfo," + _lfName + ",-1);");
	_lfInfo = _lfInfo["Global Independent"];
	for (_gb_idx = 0; _gb_idx < Columns (_lfInfo); _gb_idx = _gb_idx + 1)
	{
		ExecuteCommands (_lfInfo[_gb_idx] + ":=" + _lfInfo[_gb_idx] + "__;");
	} 	
	return 0;
}

/*---------------------------------------------------------*/
/* prompt for global variabless in a LF and fix their values */
   
function promptForGlobalParameters (_lfName)
{
	ExecuteCommands ("GetString (_lfInfo," + _lfName + ",-1);");
	_lfInfo = _lfInfo["Global Independent"];
	for (_gb_idx = 0; _gb_idx < Columns (_lfInfo); _gb_idx = _gb_idx + 1)
	{
		fprintf (stdout, "\nEnter a value for ", _lfInfo[_gb_idx], ":");
		fscanf  (stdin, "Number", _tval);
		ExecuteCommands (_lfInfo[_gb_idx] + ":=" + _tval + ";");
	} 	
	return 0;
}

/*---------------------------------------------------------*/
/* prompt for global variabless in a LF and fix their values */
   
function echoGlobalParameters (_lfName)
{
	ExecuteCommands ("GetString (_lfInfo," + _lfName + ",-1);");
	_lfInfo = _lfInfo["Global Independent"];
	for (_gb_idx = 0; _gb_idx < Columns (_lfInfo); _gb_idx = _gb_idx + 1)
	{
		ExecuteCommands ("_tval = "+_lfInfo[_gb_idx]);
		fprintf (stdout, _lfInfo[_gb_idx], " : ", Format (_tval, 12, 4), "\n");
	} 	
	return 0;
}


/*---------------------------------------------------------*/
/* take a snapshot of global parameters */
   
function stashGlobalParameters (_lfName)
{
	ExecuteCommands ("GetString (_lfInfo," + _lfName + ",-1);");
	_lfInfo = _lfInfo["Global Independent"];
	_paramStash = {};
	for (_gb_idx = 0; _gb_idx < Columns (_lfInfo); _gb_idx = _gb_idx + 1)
	{
		ExecuteCommands ("_paramStash[\""+_lfInfo[_gb_idx]+"\"] :=" + _lfInfo[_gb_idx] + ";");
	} 	
	return _paramStash;
}

/*---------------------------------------------------------*/
/* define a global parameter if not already defined */
   
function defineIfNeeded (_parName, _parValue)
{
	ExecuteCommands("GetInformation (_gb_idx, \"^`_parName`$\");");
	if (Rows (_gb_idx) == 0)
	{
		ExecuteCommands ("global `_parName`="+_parValue+";");
		return 0;
	}
	return 1;
}

/*---------------------------------------------------------*/
/* restore values of global parameters */
   

function restoreGlobalParameters (_paramStash)
{
	_stashKeys = Rows(_paramStash);
	for (_gb_idx = 0; _gb_idx < Abs (_paramStash); _gb_idx = _gb_idx + 1)
	{
		ExecuteCommands (_stashKeys[_gb_idx] + "=" + _paramStash[_stashKeys[_gb_idx]] + ";");
	} 	
	return 0;
}

/*---------------------------------------------------------*/
/* take a string row/column matrix and turn it into an AVL of 
   the form avl["matrix entry"] = 1 for each matrix entry */
   
function stringMatrixToAVL (_theList&)
{
	_gb_dim = Rows(_theList)*Columns(_theList);
	_gb_ret = {};
	for (_gb_idx = 0; _gb_idx < _gb_dim; _gb_idx = _gb_idx + 1)
	{
		_gb_ret [_theList[_gb_idx]] = _gb_idx+1;
	} 	
	return _gb_ret;
}

/*---------------------------------------------------------*/
/* take an avl indexed by 0..N-1 and convert to a row matrix */
   
function avlToMatrix (_theList&)
{
	_gb_dim = Abs(_theList);
	_gb_ret = {_gb_dim,1};
	for (_gb_idx = 0; _gb_idx < _gb_dim; _gb_idx = _gb_idx + 1)
	{
		_gb_ret [_gb_idx] = _theList[_gb_idx];
	} 	
	return _gb_ret;
}

/*---------------------------------------------------------*/
/* swap the roles of integer keys (+1) and string values in an avl */
   
function swapKeysAndValues (_theList&)
{
	_gb_ret 	= 	{};
	_gb_dim 	= 	Abs(_theList);
	_gb_keys	= 	Rows (_theList);
	
	for (_gb_idx = 0; _gb_idx < _gb_dim; _gb_idx = _gb_idx + 1)
	{
		_gb_ret [_theList[_gb_keys[_gb_idx]]] = 1+_gb_keys[_gb_idx];
	} 	
	return _gb_ret;
}


/*---------------------------------------------------------*/
/* report a string version of an a/b ratio, handling the cases
   of a and/or b = 0 */
   
function _standardizeRatio (_num, _denom)
{
	if (_denom != 0)
	{
		_ratio = _num/_denom;
		if (_ratio < 10000)
		{
			return Format (_ratio,10,4);
		}
	}
	else
	{
		if (_num == 0)
		{
			return " Undefined";
		}
	}
	return "  Infinite";
}

/*---------------------------------------------------------*/
/* copy branch lengths */
   
function _copyBranchLengths (_treeID1, _treeID2, _op, _suffix)
{
	ExecuteCommands ("_gb_dim=BranchName("+_treeID2+",-1);");
	_gb_ret = "";
	_gb_ret * 128;
	
	for (_gb_idx = 0; _gb_idx < Columns(_gb_dim)-1; _gb_idx = _gb_idx + 1)
	{
		_gb_idx2 = _treeID2 + "." + _gb_dim[_gb_idx] + "." + _suffix;
		ExecuteCommands ("_gb_idx2="+_gb_idx2);
		_gb_ret * (_treeID1 + "." + _gb_dim[_gb_idx] + "." +_suffix + ":=" + _op + _gb_idx2 + ";\n");
	} 	
	_gb_ret * 0;
	return _gb_ret;
}


/*---------------------------------------------*/
/* convert a number into a 3 letter string 
   for initializing stdin lists */
/*---------------------------------------------*/
   
function _mapNumberToString (n)
{
	if (n>=100)
	{
		return "" + n;
	}
	if (n>=10)
	{
		return "0" + n;
	}
	return "00" + n;
}

/*---------------------------------------------
 given two integer vectors of the same length,
cross-tabulate the values of each vector (vec1 
gives the row index, vec2 - the column index)
---------------------------------------------*/

function table (vec1, vec2)
{
	if (Rows(vec1) == Rows(vec2) && Columns (vec1) == 1 && Columns (vec2) == 1)
	{
		maxV = 0;
		rd = Rows(vec1);
		for (_r = 0; _r < rd; _r = _r+1)
		{
			maxV = Max (maxV, Max(vec1[_r],vec2[_r]));
		}	
		tableOut = {maxV+1, maxV+1};
		for (_r = 0; _r < rd; _r = _r+1)
		{
			tableOut[vec1[_r]][vec2[_r]] = tableOut[vec1[_r]][vec2[_r]] + 1;
		}	
		return tableOut;
	}
	else
	{
		return 0;
	} 
}

/*---------------------------------------------
 given two integer vectors of the same length,
cross-tabulate the values of each vector (vec1 
gives the row index, vec2 - the column index)
---------------------------------------------*/

function vector_max (vec1)
{
	rd   = Rows(vec1) * Columns (vec1);
	maxV = -1e100;
	maxI = 0;
	for (_r = 0; _r < rd; _r = _r+1)
	{
		if (vec1[_r] > maxV)
		{
			maxV = vec1[_r];
			maxI = _r;
		}
	}	
	return {{maxV__,maxI__}};
}


/*---------------------------------------------
construct the frequency vector and a LF mixing componenent for 
a general discrete distribution on numberOfRates rates
---------------------------------------------*/

function generate_gdd_freqs (numberOfRates, freqs&, lfMixing&, probPrefix, incrementFlag)
{
	freqs    	 = {numberOfRates,1};
	lfMixing	 = ""; lfMixing * 128; lfMixing * "Log(";
	
	if (numberOfRates == 1)
	{
		freqs[0] = "1";
	}
	else
	{
		if (incrementFlag)
		{
			mi = numberOfRates-1;
			ExecuteCommands ("global "+probPrefix+"_"+mi+":<1;"+probPrefix+"_"+mi+"=1/2;");
			for (mi=1; mi<numberOfRates-1; mi=mi+1)
			{
				ExecuteCommands (""+probPrefix+"_"+mi+" = "+probPrefix+"_"+mi+"*(1-1/numberOfRates);");
			}
			
		}
		else
		{
			for (mi=1; mi<numberOfRates; mi=mi+1)
			{
				ExecuteCommands ("global "+probPrefix+"_"+mi+":<1;"+probPrefix+"_"+mi+" = 1/" + (numberOfRates-mi+1));
			}
		}
		
		freqs[0] 	 = ""+probPrefix+"_1";
		for (mi=1; mi<numberOfRates-1; mi=mi+1)
		{
			freqs[mi] = "";
			for (mi2=1;mi2<=mi;mi2=mi2+1)
			{
				freqs[mi] = freqs[mi]+"(1-"+probPrefix+"_"+mi2+")";		
			}
			freqs[mi] = freqs[mi]+""+probPrefix+"_"+(mi+1);	
		}	
	
		freqs[mi] = "";
		for (mi2=1;mi2<mi;mi2=mi2+1)
		{
			freqs[mi] = freqs[mi]+"(1-"+probPrefix+"_"+mi2+")";		
		}
		freqs[mi] = freqs[mi]+"(1-"+probPrefix+"_"+mi+")";	
	}
	
	lfMixing * ("SITE_LIKELIHOOD[0]*"+freqs[0]);
	for (mi = 1; mi < numberOfRates; mi=mi+1)
	{
		lfMixing * ("+SITE_LIKELIHOOD[" + mi + "]*" + freqs[mi]);
	}
	lfMixing * ")";
	lfMixing * 0;
	return 0;
}