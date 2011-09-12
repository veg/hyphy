_HY_DBW_OUTPUT_PROCESSOR_NAME_ = "File, SQL Instructions";

/*--------------------------------------------------*/

function _HY_DBW_OUTPUT_PROCESSOR_FUNCTION_ (dummy)
{
	colCount = Columns (SQL_ROW_DATA);
	
	if (colCount)
	{
		cString = "";
		cString * 128;
		cString * ("\nINSERT INTO "+tableID+ " (");
		cString  * SQL_COLUMN_NAMES[0];
		for (cc=1; cc<colCount; cc=cc+1)
		{
			cString  * ",";
			cString  * SQL_COLUMN_NAMES[cc];
		}
		cString * ") VALUES ('";
		
		cString  * (SQL_ROW_DATA[0]);
		cString  * "'";
		for (cc=1; cc<colCount; cc=cc+1)
		{
			cString  * ", '";
			cString  * (SQL_ROW_DATA[cc]);
			cString  * "'";
		}
		cString * ");\n";
		cString * 0;
		
		fprintf (_HY_DBW_TAB_DELIMITED_FILE_,cString);
	}
	return 0;
}

/*--------------------------------------------------*/

function _HY_DBW_OUTPUT_PROCESSOR_FUNCTION_PRE_ (dummy)
{
	colCount = Columns (SQL_ROW_DATA);
	
	if (colCount)
	{
		cName = SQL_ROW_DATA[0];
		sqliteTables[cName] = SQL_ROW_DATA[1];	
	}
	return 0;
}

/*--------------------------------------------------*/

if (_HY_DBW_OUTPUT_RUN_ME_)
{
	SetDialogPrompt ("SQL text file to export to:");
	fprintf (PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN);
	_HY_DBW_TAB_DELIMITED_FILE_ = LAST_FILE_PATH;
	
	sqliteTables = {};
	DoSQL (_HY_DBW_OUTPUT_DB_ID_, "SELECT name,sql from SQLITE_MASTER WHERE length(sql)>0", "return _HY_DBW_OUTPUT_PROCESSOR_FUNCTION_PRE_ (0);");	
	
	fprintf 		(stdout, "[FOUND ", Abs(sqliteTables), " TABLES]\n");
	
	tableChoice = {Abs(sqliteTables),2};
	tn 			=  Rows(sqliteTables);
	
	for (cc = 0; cc < Abs(sqliteTables); cc = cc+1)
	{
		tableChoice[cc][0] = tn[cc];
		tableChoice[cc][1] = "Table "+ tn[cc];
		
	}	

	ChoiceList		(whichTable,"Which Table?",1, SKIP_NONE,tableChoice);
	if (whichTable < 0)
	{
		return 0;
	}
	
	ChoiceList		(exportTable,"Export Table Definition?",1, SKIP_NONE, "Yes", "Include SQL code to create the table",
																		   "No","Export only record definitions");
	
	if (exportTable<0)
	{
		return 0;
	}
	
	tableID = tableChoice[whichTable][0];
	
	if (exportTable==0)
	{
		fprintf (_HY_DBW_TAB_DELIMITED_FILE_, sqliteTables[tableID], ";\n");
	}
	DoSQL (_HY_DBW_OUTPUT_DB_ID_, "SELECT * from " + tableID, "return _HY_DBW_OUTPUT_PROCESSOR_FUNCTION_ (0);");	
	fprintf (_HY_DBW_TAB_DELIMITED_FILE_,CLOSE_FILE);
}
