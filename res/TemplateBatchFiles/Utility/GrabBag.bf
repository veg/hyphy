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

function stratifyAVLByValuesAux (key,value)
{
	if (Abs (_gb_resAVL[value]) == 0)
	{
		_gb_resAVL[value] = {};
	}
	
	_gb_resAVL[value] + key;
	
	return 0;
	
}

function stratifyAVLByValues (_gb_anAVL)
{
	_gb_count	  = Abs (_gb_anAVL);
	_gb_resAVL	  = {};
	_gb_anAVL["stratifyAVLByValuesAux"][""];
	return _gb_resAVL;
}

/*---------------------------------------------------------*/
/* Stratify a matrix by values; store indices for each unique value */

function stratifyMatrixByValues (_gb_aMatrix)
{
	_gb_count	  = Rows (_gb_aMatrix)*Columns(_gb_aMatrix);
	_gb_resAVL	  = {};
	for (_gb_idx = 0; _gb_idx < _gb_count; _gb_idx += 1)
	{
		_gb_key = _gb_aMatrix[_gb_idx];
		if (Abs (_gb_resAVL[_gb_key]) == 0)
		{
			_gb_resAVL[_gb_key] = {};
		}
		(_gb_resAVL[_gb_key]) + _gb_idx;
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
/* Split a file path into DIRECTORY, FILENAME and EXTENSION  */
   
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
   
function fixGlobalParameters (_lfName) {
	GetString (_lfInfo,^_lfName,-1);
	_lfInfo = _lfInfo["Global Independent"];
	for (_gb_idx = 0; _gb_idx < Columns (_lfInfo); _gb_idx += 1) {
		ExecuteCommands (_lfInfo[_gb_idx] + ":=" + _lfInfo[_gb_idx] + "__;");
	} 	
	return 0;
}

/*---------------------------------------------------------*/
/* unconstrain global variables in a LF at their current values */
   
function unconstrainGlobalParameters (_lfName) {
	GetString (_lfInfo,^_lfName,-1);
	_lfInfo = _lfInfo["Global Constrained"];
	for (_gb_idx = 0; _gb_idx < Columns (_lfInfo); _gb_idx += 1) {
		ExecuteCommands (_lfInfo[_gb_idx] + "=" + _lfInfo[_gb_idx]);
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
/* echo global parameters */
   
function echoGlobalParameters (_lfName)
{
	ExecuteCommands ("GetString (_lfInfo," + _lfName + ",-1);");
	_lfInfo = _lfInfo["Global Independent"];
	for (_gb_idx = 0; _gb_idx < Columns (_lfInfo); _gb_idx = _gb_idx + 1)
	{
		ExecuteCommands ("_tval = "+_lfInfo[_gb_idx]);
		fprintf (stdout, _lfInfo[_gb_idx], " : ", Format (_tval, 12, 4), "\n");
	} 	
	return Columns (_lfInfo);
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
		ExecuteCommands ("_paramStash[\""+_lfInfo[_gb_idx]+"\"] =" + _lfInfo[_gb_idx] + ";");
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
/* export a list of variables into the string of the form 
   
   varID = varValue; 

*/
   
function exportVarList (_varList)
{
	_exportString = ""; _exportString * 256;
	
	for (_idx = 0; _idx < Columns (_varList) * Rows (_varList); _idx += 1)
	{
		_exportString * (_varList [_idx] + " = " + Eval (_varList[_idx]) + ";\n");
	}
	
	_exportString * 0;
	return _exportString;
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
/* take a string row/columqn matrix and turn it into an AVL of 
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
/* take an avl indexed by 0..N-1 and convert to a choice list matrix */
   
function avlToChoiceMatrix (_theList&)
{
	_gb_dim = Abs(_theList);
	_gb_ret = {_gb_dim,2};
	for (_gb_idx = 0; _gb_idx < _gb_dim; _gb_idx = _gb_idx + 1)
	{
		_gb_ret [_gb_idx][0] = _theList[_gb_idx];
		_gb_ret [_gb_idx][1] = _theList[_gb_idx];
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
 the analog of Python's string.join (list)
---------------------------------------------*/

function string_join (sep, list)
{
	_result = ""; _result * 128;
	
	if (Type (list) == "Matrix")
	{
		_dim = Rows(list)*Columns(list);
	}
	else
	{
		_dim = Abs (list);
	}
	
	if (_dim)
	{
		_result * ("" + list[0]);
		for (_r = 1; _r < _dim; _r = _r + 1)
		{
			_result * (sep + list[_r]);
		
		}
	}
	
	_result * 0;
	return _result;
}

/*---------------------------------------------
 prompt for a value in a given range, given
 a default value. 
---------------------------------------------*/

function prompt_for_a_value (prompt,default,lowerB,upperB,isInteger)
{
	value = lowerB-1;
	while (value < lowerB || value > upperB)
	{
		fprintf (stdout, "<<", prompt, " (permissible range = [", lowerB, ",", upperB, "], default value = ", default);
		if (isInteger)
		{
			fprintf (stdout, ", integer"); 
		}
		fprintf (stdout, ")>>");
		fscanf  (stdin, "String", strVal);
		if (Abs(strVal) == 0)
		{
			value = 0+default;
			break;
		}	
		value = 0+strVal;
	}
	if (isInteger)
		return value$1;
		
	return value;
}

/*---------------------------------------------
take an AVL of the form ["string"] = number
and print it as:

key[_sepChar]+: number (%)

---------------------------------------------*/
	
function _printAnAVL (_theList, _sepChar)
{
	_gb_keys 		= _sortStrings(Rows (_theList));
	_printAnAVLInt (_theList, _gb_keys, _sepChar, 0);
		
	return 0;
}

/*---------------------------------------------
take an AVL of the form ["string"] = number
and print it as:

key[_sepChar]+: number (%)

add a "Total" row

---------------------------------------------*/
	
function _printAnAVLTotal (_theList, _sepChar)
{
	_gb_keys 		= _sortStrings(Rows (_theList));
	_printAnAVLInt (_theList, _gb_keys, _sepChar, 1);
		
	return 0;
}

/*---------------------------------------------
take an AVL of the form [number] = number
and print it as:

key[_sepChar]+: number (%)

---------------------------------------------*/
	
function _printAnAVLNumeric (_theList, _sepChar)
{
	_gb_dim   		= Abs(_theList);
	num_keys		= avlKeysToMatrix (_theList)%0;
	_gb_keys 		= {_gb_dim,1};
	for (_gb_idx = 0; _gb_idx < _gb_dim; _gb_idx = _gb_idx + 1)
	{
		_gb_keys[_gb_idx] = ""+num_keys[_gb_idx];
	}
	
	_printAnAVLInt (_theList, _gb_keys, _sepChar, 0);
		
	return 0;
}

function _printAnAVLNumericTotal (_theList, _sepChar)
{
	_gb_dim   		= Abs(_theList);
	num_keys		= avlKeysToMatrix (_theList)%0;
	_gb_keys 		= {_gb_dim,1};
	for (_gb_idx = 0; _gb_idx < _gb_dim; _gb_idx = _gb_idx + 1)
	{
		_gb_keys[_gb_idx] = ""+num_keys[_gb_idx];
	}
	
	_printAnAVLInt (_theList, _gb_keys, _sepChar, 1);
		
	return 0;
}

/*---------------------------------------------*/

function _printAnAVLInt (_theList, _gb_keys, _sepChar, _doTotal)
{	
	_gb_dim   		= Abs(_theList);
	_gb_total 		= 0;
	_gb_max_key_len = 0;
	if (_doTotal)
	{
		gb_max_key_len = 5;
	}

	for (_gb_idx = 0; _gb_idx < _gb_dim; _gb_idx = _gb_idx + 1)
	{
		_gb_key 		= _gb_keys[_gb_idx];
		_gb_max_key_len = Max (_gb_max_key_len, Abs(_gb_key));
		_gb_total 		= _gb_total + _theList[_gb_key];
	}
	
	fprintf (stdout, "\n");
	for (_gb_idx = 0; _gb_idx < _gb_dim; _gb_idx = _gb_idx + 1)
	{
		_gb_key 		= _gb_keys[_gb_idx];
		fprintf (stdout, _gb_key);
		for (_gb_idx2 = Abs(_gb_key); _gb_idx2 < _gb_max_key_len; _gb_idx2 = _gb_idx2 + 1)
		{
			fprintf (stdout, _sepChar);
		}
		fprintf (stdout, ":", Format (_theList[_gb_key],8,0), " (", Format (100*_theList[_gb_key]/_gb_total,5,2), "%)\n");
	}
	
	if (_doTotal)
	{
		fprintf (stdout, "Total");
		for (_gb_idx2 = 5; _gb_idx2 < _gb_max_key_len; _gb_idx2 = _gb_idx2 + 1)
		{
			fprintf (stdout, _sepChar);
		}
		fprintf (stdout, ":", Format (_gb_total,8,0),"\n");
	}
		
	return 0;
}

/*---------------------------------------------
sort a matrix of strings; return a 
column vector
---------------------------------------------*/
function _sortStringsAux (theKey, theValue)
{
	for (_gb_idx2 = 0; _gb_idx2 < theValue; _gb_idx2=_gb_idx2+1)
	{
		_gb_sortedStrings [_gb_idx] = theKey;
		_gb_idx = _gb_idx + 1;
	}
	return 0;
}

function _sortStrings (_theList)
{
	_gb_dim = Rows (_theList)*Columns (_theList);
	_toSort = {};
	for (_gb_idx = 0; _gb_idx < _gb_dim; _gb_idx = _gb_idx + 1)
	{
		_toSort[_theList[_gb_idx]] = _toSort[_theList[_gb_idx]]+1;
	}
	_gb_sortedStrings = {_gb_dim,1};
	_gb_idx = 0;
	_toSort["_sortStringsAux"][""];
	return _gb_sortedStrings;
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
			ExecuteCommands ("global "+probPrefix+"_"+mi+":>0;");
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
				ExecuteCommands ("global "+probPrefix+"_"+mi+":>0;");
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

/*---------------------------------------------
reverse complement a nucleotide string
---------------------------------------------*/

_nucleotide_rc = {};
_nucleotide_rc["A"] = "T";
_nucleotide_rc["C"] = "G";
_nucleotide_rc["G"] = "C";
_nucleotide_rc["T"] = "A";
_nucleotide_rc["M"] = "K";
_nucleotide_rc["R"] = "Y";
_nucleotide_rc["W"] = "W";
_nucleotide_rc["S"] = "S";
_nucleotide_rc["Y"] = "R";
_nucleotide_rc["K"] = "M";
_nucleotide_rc["B"] = "V";  /* not A */
_nucleotide_rc["D"] = "H";  /* not C */
_nucleotide_rc["H"] = "D";  /* not G */
_nucleotide_rc["V"] = "B";  /* not T */
_nucleotide_rc["N"] = "N";

function nucleotideReverseComplement (seqIn)
{
	_seqOut = "";_seqOut*128;
	_seqL   = Abs(seqIn);
	for (_r = _seqL-1; _r >=0 ; _r = _r-1)
	{
		_seqOut *_nucleotide_rc[seqIn[_r]];
	}
	_seqOut*0;
	return _seqOut;
}

/*---------------------------------------------------------------------*/

function RankMatrix (matrix)
/* take a sorted row matrix and return a rank row matrix averaging ranks for tied values */
{
	lastValue				   			 = matrix[0];
	lastIndex				   			 = 0;
	matrix							 [0] = 0;
	
	sampleCount = Rows (matrix);
	
	for (_i = 1; _i < sampleCount; _i = _i+1)
	{
		if (lastValue != matrix[_i])
		{
			meanIndex = _i - lastIndex;
			lastValue = matrix[_i];
			if (meanIndex > 1)
			{
				meanIndex = (lastIndex + _i - 1) * meanIndex / (2 * meanIndex);
				for (_j = lastIndex; _j < _i; _j = _j + 1)
				{
					matrix[_j] = meanIndex;
				}
			}
			matrix[_i] = _i;
			lastIndex = _i;
		}
	}
	
	meanIndex = _i - lastIndex;
	if (meanIndex > 1)
	{
		meanIndex = (lastIndex + _i - 1) * meanIndex / (2 * meanIndex);
		for (_j = lastIndex; _j < _i; _j = _j + 1)
		{
			matrix[_j] = meanIndex;
		}
	}
	else
	{
		matrix[_i-1] = _i - 1;
	}
	return matrix;
}

/*---------------------------------------------------------------------*/

function mapSets (sourceList,targetList)
// source ID -> target ID (-1 means no correspondence)

{
	targetIndexing = {};
	_d = Rows (targetList) * Columns (targetList);
	
	for (_i = 0; _i < _d; _i += 1)
	{
		targetIndexing [targetList[_i]] = _i + 1;
	}
	_d = Rows (sourceList) * Columns (sourceList);
	mapping 	  = {1,_d};
	for (_i = 0; _i < _d; _i += 1)
	{
		mapping [_i] = targetIndexing[sourceList[_i]] - 1;
	}
	
	return mapping;
}

/*---------------------------------------------------------------------*/

function mapStrings (sourceStr,targetStr)
// source ID -> target ID (-1 means no correspondence)

{
	mapping 	  = {};
	targetIndexing = {};
	_d = Abs(targetStr);
	
	for (_i = 0; _i < _d; _i += 1)
	{
		targetIndexing [targetStr[_i]] = _i + 1;
	}
	_d = Abs (sourceStr);
	for (_i = 0; _i < _d; _i += 1)
	{
		mapping [_i] = targetIndexing[sourceStr[_i]] - 1;
	}
	
	return mapping;
}
/*---------------------------------------------------------------------*/

function remapSequenceCoordinatesToReference (ref,seq)
{
	_seqLen	  = Abs(seq);
	_coordMap = {_seqLen,1}["-1"];
	
		
	_k				= (ref$"^\\-+");
	_referenceSpan	= _k[1]+1;
	
	for (_k = 0; _k < _referenceSpan; _k = _k+1)
	{
		_coordMap[_k] = 0;
	}
	
	_qryCoord = _k;
	_refCoord = 0;

	while (_k < Abs(seq))
	{
		if (seq[_k] != "-")
		{
			_coordMap[_qryCoord] = _refCoord;
			_qryCoord = _qryCoord + 1;
		}
		if (ref[_k] != "-")
		{
			_refCoord = _refCoord + 1;
		}
		_k = _k+1;
	}
	return _coordMap;
}

/*---------------------------------------------------------------------*/

function runModule (module,options,suppressStdout)
{
	if (suppressStdout)
	{
		_gfr = GLOBAL_FPRINTF_REDIRECT;
		GLOBAL_FPRINTF_REDIRECT = "/dev/null";
	}
	LoadFunctionLibrary (module,options);
	if (suppressStdout)
	{
		GLOBAL_FPRINTF_REDIRECT = _gfr;
	}
}

/*---------------------------------------------------------------------*/

function _formatTimeString (secondCount)
{
	_hours = secondCount $3600;
	if (_hours < 10)
	{
		_timeString = "0"+_hours;
	}
	else
	{
		_timeString = ""+_hours;
	}
	_minutes = (secondCount%3600)$60;
	if (_minutes < 10)
	{
		_timeString = _timeString+":0"+_minutes;
	}
	else
	{
		_timeString = _timeString+":"+_minutes;
	}
	_seconds = (secondCount%60);
	if (_seconds < 10)
	{
		_timeString = _timeString+":0"+_seconds;
	}
	else
	{
		_timeString = _timeString+":"+_seconds;
	}
	return _timeString;
}	

/*---------------------------------------------------------------------*/

lfunction _constrainVariablesAndDescendants (variable) {
    GetInformation (allVars, "^" + (variable&&6) + "\\..+$");
    for (k = 0; k < Columns (allVars); k += 1) {
        variableID    = allVars[k];
        current_value = ^variableID;
        ^variableID := current_value__;
    }
    return 0;
}

/*---------------------------------------------------------------------*/

lfunction _unconstrainVariablesAndDescendants (variable) {
    GetInformation (allVars, "^" + (variable&&6) + "\\..+$");
    for (k = 0; k < Columns (allVars); k += 1) {
        variableID    = allVars[k];
        ClearConstraints (^variableID);
    }
    return 0;
}
