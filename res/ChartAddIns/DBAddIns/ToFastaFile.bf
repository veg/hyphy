_HY_DBW_OUTPUT_PROCESSOR_NAME_ = "File, FASTA";

/*--------------------------------------------------*/

function _HY_DBW_OUTPUT_PROCESSOR_FUNCTION_ (dummy)
{
	colCount = Columns (SQL_ROW_DATA);
	
	if (colCount >= 2)
	{
		cString = "";
		cString * 128;		
		cString  * (">"+SQL_ROW_DATA[0]+"\n");
		cString  * (SQL_ROW_DATA[1]+"\n");
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
	fprintf (PROMPT_FOR_FILE,CLEAR_FILE);
	_HY_DBW_TAB_DELIMITED_FILE_ = LAST_FILE_PATH;
	recordCounter = 0;
	DoSQL (_HY_DBW_OUTPUT_DB_ID_, _HY_DBW_SQL_QUERY, "return _HY_DBW_OUTPUT_PROCESSOR_FUNCTION_ (0);");
}
