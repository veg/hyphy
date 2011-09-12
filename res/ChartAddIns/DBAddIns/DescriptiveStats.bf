_HY_DBW_OUTPUT_PROCESSOR_NAME_ = "Descriptive Statistics";

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
	colCount = Columns(_CW_columnHeaders);
	
	if (rowCount && colCount == 1)
	{
		_CW_datatable = {rowCount,1}["0 + (_CW_dataRows [_MATRIX_ELEMENT_ROW_])[0]"];
		_CW_dataRows = 0;
		
		LoadFunctionLibrary   ("DescriptiveStatistics");
		LoadFunctionLibrary   ("PS_Plotters.bf");
		
		fprintf (stdout, "Quantity label:");
		fscanf	(stdin, "String", dataLabel);
		histPlot = PSHistogram (_CW_datatable,0,1,"TimesNewRoman",{{400,400,12}},{{0.3,0.3,0.3}},{{"",dataLabel,"Frequency"}},1);
		DEFAULT_FILE_SAVE_NAME = "Histogram.ps";
		SetDialogPrompt ("Write a Postscript histogram file to:");
		fprintf (PROMPT_FOR_FILE, CLEAR_FILE, histPlot);
		
		PrintDescriptiveStats (dataLabel, GatherDescriptiveStats (_CW_datatable));
	}
}
