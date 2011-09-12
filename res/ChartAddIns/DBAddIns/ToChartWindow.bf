_HY_DBW_OUTPUT_PROCESSOR_NAME_ = "Chart window";

/*--------------------------------------------------*/

function _HY_DBW_OUTPUT_PROCESSOR_FUNCTION_ (dummy)
{
	colCount = Columns (SQL_ROW_DATA);
	
	if (Abs(_CW_columnHeaders[0]) == 0)
	{
		_CW_columnHeaders = SQL_COLUMN_NAMES;
	}
	
	_CW_dataRows [Abs(_CW_dataRows)] = SQL_ROW_DATA;
	return 0;
}

/*--------------------------------------------------*/

if (_HY_DBW_OUTPUT_RUN_ME_)	
{
	recordCounter = 0;
	
	_CW_columnHeaders = {};
	_CW_dataRows	  = {};
	
	DoSQL (_HY_DBW_OUTPUT_DB_ID_, _HY_DBW_SQL_QUERY, "return _HY_DBW_OUTPUT_PROCESSOR_FUNCTION_ (0);");
	
	rowCount = Abs (_CW_dataRows);
	if (rowCount)
	{
		_CW_datatable = {rowCount,colCount};
		
		for (ri = 0; ri < rowCount; ri = ri+1)
		{
			dRow = _CW_dataRows[ri];
			
			for (ci = 0; ci < colCount; ci = ci+1)
			{
				_CW_datatable [ri][ci] = 0 + dRow [ci];
			}
		}
		
		_CW_dataRows = 0;
		
		OpenWindow (CHARTWINDOW,{{"SQL Table Extract"}
								   {"_CW_columnHeaders"},
								   {"_CW_datatable"},
								   {"None"},
								   {"Index"}},
								   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");
	}
}
