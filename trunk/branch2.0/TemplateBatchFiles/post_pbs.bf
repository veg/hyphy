nOfI = 0;

while (nOfI<=0)
{
	fprintf (stdout,"\nHow many data replicates should be generated?:");
	fscanf  (stdin,"Number",nOfI);
}  

SetDialogPrompt ("Save detailed bootstrap results to:");

fprintf (PROMPT_FOR_FILE,CLEAR_FILE);

MESSAGE_LOGGING = 0;
dumb = BootStrapFunction(nOfI,LAST_FILE_PATH,1);
MESSAGE_LOGGING = 1;
