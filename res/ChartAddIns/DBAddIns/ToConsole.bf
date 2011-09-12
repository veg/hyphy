_HY_DBW_OUTPUT_PROCESSOR_NAME_ = "Console, by record";

/*--------------------------------------------------*/

function _HY_DBW_OUTPUT_PROCESSOR_FUNCTION_ (dummy)
{
	colCount = Columns (SQL_ROW_DATA);
	
	recordCounter = recordCounter + 1;
	fprintf (stdout, "\nRECORD ", recordCounter, "\n");
	
	for (cc=0; cc<colCount; cc=cc+1)
	{
		cName  = SQL_COLUMN_NAMES[cc];
		cChars = Abs(cName);
		
		underLine = "";
		
		for (ck=0; ck<cChars;ck=ck+1)
		{
			underLine = underLine+"-";
		}	
		
		fprintf (stdout, "\n", SQL_COLUMN_NAMES[cc],"\n",underLine,"\n", SQL_ROW_DATA[cc]&&2,"\n");
	}
	return 0;
}

/*--------------------------------------------------*/

if (_HY_DBW_OUTPUT_RUN_ME_)
{
	recordCounter = 0;
	DoSQL (_HY_DBW_OUTPUT_DB_ID_, _HY_DBW_SQL_QUERY, "return _HY_DBW_OUTPUT_PROCESSOR_FUNCTION_ (0);");
}
