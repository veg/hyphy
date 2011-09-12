_HY_DBW_OUTPUT_PROCESSOR_NAME_ = "File, Frontier Header FASTA";

/*--------------------------------------------------*/

function _HY_DBW_OUTPUT_PROCESSOR_FUNCTION_ (dummy)
{
	colCount = Columns (SQL_ROW_DATA);
	
	if (colCount == 2)
	{
		cString = "";
		cString * 128;		
		headerSplit = splitOnRegExp(SQL_ROW_DATA[0],"\\|");
		
		if (SQL_COLUMN_NAMES[1] == "RAW")
		{
			seqno = 1;
		}
		else
		{
			seqno = 2;
		}
		
		cString  * (">"+headerSplit[0]+" "+headerSplit[1]+" seqno=" + seqno + "\n");
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
	fprintf (PROMPT_FOR_FILE,CLEAR_FILE, KEEP_OPEN);
	_HY_DBW_TAB_DELIMITED_FILE_ = LAST_FILE_PATH;
	recordCounter = 0;
	DoSQL (_HY_DBW_OUTPUT_DB_ID_, _HY_DBW_SQL_QUERY, "return _HY_DBW_OUTPUT_PROCESSOR_FUNCTION_ (0);");
	fprintf (_HY_DBW_TAB_DELIMITED_FILE_,CLOSE_FILE);
}

/*----------------------------------------------------------------*/

function splitOnRegExp (string, splitter)
{
	matched = string || splitter;
	splitBits = {};
	fromPos = 0;
	toPos	= 0;
	for (mc = 0; mc < Rows (matched); mc = mc+2)
	{
		toPos = matched[mc]-1;
		splitBits [Abs(splitBits)] = string[fromPos][toPos];
		fromPos    = matched[mc+1]+1;
	}
	splitBits [Abs(splitBits)] = string[fromPos][Abs(string)-1];
	return splitBits;
}
