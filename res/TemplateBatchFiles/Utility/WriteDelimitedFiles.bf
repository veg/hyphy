function WriteSeparatedTable (fileName, column_headers, numeric_data, first_column, delimiter)
{
	if (Abs(fileName) == 0)
	{
		fprintf (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN);
		_f2r = LAST_FILE_PATH; 
	}
	else
	{
		fprintf (fileName, CLEAR_FILE, KEEP_OPEN);
		_f2r = fileName; 
	}
	
	_rd = Columns (column_headers);
	if (_rd)
	{
		fprintf (_f2r, column_headers[0]);
		for (_gb_idx = 1; _gb_idx < _rd; _gb_idx = _gb_idx + 1)
		{
			fprintf (_f2r, delimiter, column_headers[_gb_idx]);
		}
	}
	_rd2 = Abs (first_column);
	if (_rd2)
	{
		for (_gb_idx = 0; _gb_idx < _rd2; _gb_idx = _gb_idx + 1)
		{
			fprintf (_f2r, "\n", first_column[_gb_idx]);
			for (_gb_idx2 = 1; _gb_idx2 < _rd; _gb_idx2 = _gb_idx2 + 1)
			{
				fprintf (_f2r, delimiter, numeric_data[_gb_idx][_gb_idx2-1]);		
			}
		}
	}
	else
	{
		_rd2 = Rows (numeric_data);
		if (_rd == 0)
		{
			_rd = Columns (numeric_data);
		}
		for (_gb_idx = 0; _gb_idx < _rd2; _gb_idx = _gb_idx + 1)
		{
			fprintf (_f2r, "\n", numeric_data[_gb_idx][0]);
			for (_gb_idx2 = 1; _gb_idx2 < _rd; _gb_idx2 = _gb_idx2 + 1)
			{
				fprintf (_f2r, delimiter, numeric_data[_gb_idx][_gb_idx2]);		
			}
		}
	}
	
	fprintf (_f2r, CLOSE_FILE);
	return 0;
}

