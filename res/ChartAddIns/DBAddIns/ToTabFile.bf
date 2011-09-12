_HY_DBW_OUTPUT_PROCESSOR_NAME_ = "File, tab delimited";

/*--------------------------------------------------*/

function _HY_DBW_OUTPUT_PROCESSOR_FUNCTION_ (dummy)
{
	colCount = Columns (SQL_ROW_DATA);
	
	if (colCount)
	{
		cString = "";
		cString * 128;
		if (recordCounter == 0)
		{
			cString  * SQL_COLUMN_NAMES[0];
			for (cc=1; cc<colCount; cc=cc+1)
			{
				cString  * "\t";
				cString  * SQL_COLUMN_NAMES[cc];
			}
			cString * "\n";
		}
		
		cString  * (SQL_ROW_DATA[0]&&2);
		for (cc=1; cc<colCount; cc=cc+1)
		{
			cString  * "\t";
			cString  * (SQL_ROW_DATA[cc]&&2);
		}
		cString * 0;
		
		fprintf (_HY_DBW_TAB_DELIMITED_FILE_,cString,"\n");
		recordCounter = 1;
	}
	return 0;
}

/*--------------------------------------------------*/

if (_HY_DBW_OUTPUT_RUN_ME_)
{
	SetDialogPrompt ("File to export to:");
	fprintf (PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN);
	_HY_DBW_TAB_DELIMITED_FILE_ = LAST_FILE_PATH;
	recordCounter = 0;
	DoSQL (_HY_DBW_OUTPUT_DB_ID_, _HY_DBW_SQL_QUERY, "return _HY_DBW_OUTPUT_PROCESSOR_FUNCTION_ (0);");
	fprintf (_HY_DBW_TAB_DELIMITED_FILE_,CLOSE_FILE);
}
