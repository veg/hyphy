_HY_DBW_OUTPUT_PROCESSOR_NAME_ = "Import Records from SQL file";

if (_HY_DBW_OUTPUT_RUN_ME_)	
{
	SetDialogPrompt ("SQL to import from:");
	fscanf 			(PROMPT_FOR_FILE,"Raw",inSQL);

	DoSQL 			(_HY_DBW_OUTPUT_DB_ID_, inSQL, "return 0;");	
}

/*--------------------------------------------------*/

function _HY_DBW_OUTPUT_PROCESSOR_FUNCTION_ (dummy)
{
	return 0;
}
