RequireVersion ("0.9920060830");

SetDialogPrompt ("Please select the index file:");
fscanf 			(PROMPT_FOR_FILE,"Lines",inLines);

counter	= Columns (inLines);

if (skipCodeSelectionStep)
{
	stringMatrix 	= {counter$2,1};
	codeTableMatrix = {counter$2,1};
	counter  = 0;
	counter2 = 0;
	for (k=0; k<counter; k=k+2)
	{
		codeTableMatrix [k$2] = 0+inLines[k+1];
		stringMatrix 	[k$2] = inLines[k];
	}
}
else
{
	stringMatrix = Transpose (inLines);
}
