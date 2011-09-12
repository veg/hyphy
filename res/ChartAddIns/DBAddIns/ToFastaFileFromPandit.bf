_HY_DBW_OUTPUT_PROCESSOR_NAME_ = "File, FASTA from PANDIT";
ExecuteAFile		(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"ReadDelimitedFiles.bf");
ExecuteAFile		(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"GrabBag.bf");
ExecuteAFile		(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TreeTools.ibf");

VERBOSITY_LEVEL = -1;
/*--------------------------------------------------*/

function _HY_DBW_OUTPUT_PROCESSOR_FUNCTION_ (dummy)
{
	colCount = Columns (SQL_ROW_DATA);
	
	if (colCount == 4)
	{
		Topology		T = SQL_ROW_DATA[1];
		tAVL 			  = T^0;

		usedLeafLabels = {};
		oldToNewMap	   = {};

		for (k=1; k<Abs(tAVL); k=k+1)
		{
			cName 		=	 	(tAVL[k])["Name"];
			nName		= 		normalizeSequenceID (cName, "usedLeafLabels");
			if (nName != cName)
			{
				(tAVL[k])["Name"] = nName;
			}
			oldToNewMap[cName] = nName;
		}

		DATAFILE_TREE 			= PostOrderAVL2String(tAVL);

		DataSet testDS 			= ReadFromString (SQL_ROW_DATA[3]);

		for (k=0; k<testDS.species; k=k+1)
		{
			GetString    (cName, testDS,k);
			SetParameter (testDS,k,oldToNewMap[cName]);
		}			
	
		IS_TREE_PRESENT_IN_DATA = 1;
		
		goodPositions = SQL_ROW_DATA[2]||"x|X";
		includeThese  = stringMatrixToAVL ("goodPositions");
		DataSetFilter theF = CreateFilter (testDS,1,includeThese[siteIndex]);
		
		if (recordCount == 0)
		{
			DEFAULT_FILE_SAVE_NAME		= SQL_ROW_DATA[0] + ".nex";
			SetDialogPrompt 			 ("File to export to:");
			fprintf 				 	 (PROMPT_FOR_FILE,CLEAR_FILE,theF);
			_HY_DBW_TAB_DELIMITED_FILE_ = LAST_FILE_PATH $ ("\\"+DIRECTORY_SEPARATOR+"[^\\"+DIRECTORY_SEPARATOR+"]+$");
			_HY_DBW_TAB_DELIMITED_FILE_	= LAST_FILE_PATH[0][_HY_DBW_TAB_DELIMITED_FILE_[0]];
		}
		else
		{
			_HY_DBW_TAB_DELIMITED_FILE2_ = _HY_DBW_TAB_DELIMITED_FILE_ + SQL_ROW_DATA[0] + ".nex";
			fprintf 				 	 (_HY_DBW_TAB_DELIMITED_FILE2_,CLEAR_FILE,theF);
		}
		recordCount = 1;
	}
	return 0;
}

/*--------------------------------------------------*/

if (_HY_DBW_OUTPUT_RUN_ME_)
{
	recordCount = 0;
	_HY_DBW_TAB_DELIMITED_FILE_ = LAST_FILE_PATH;
	DoSQL (_HY_DBW_OUTPUT_DB_ID_, _HY_DBW_SQL_QUERY, "return _HY_DBW_OUTPUT_PROCESSOR_FUNCTION_ (0);");
}
