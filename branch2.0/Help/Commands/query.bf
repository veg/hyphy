classStyle = "\n\\section*{_INSERT_DATA_FOR_Command_HERE_}\\begin{tabular}{lr}\\parbox[t]{2in}{\\textit{_INSERT_DATA_FOR_Class_HERE_\\newline _INSERT_DATA_FOR_Revised_HERE_}}&\\parbox{4in}{\\textsf{_INSERT_DATA_FOR_Description_HERE_}}\\end{tabular}"+
"\n\n\\texttt{_INSERT_DATA_FOR_Syntax_HERE_}\\scriptsize{\\begin{flushleft}_INSERT_DATA_FOR_Notes_HERE_\\end{flushleft}}Example\\linebreak\\linespread{0.3}\\fbox{\\parbox{6in}{\\tiny\\begin{alltt}_INSERT_DATA_FOR_Example_HERE_\\end{alltt}}}\\linespread{1}\n";

latexTop    = "\\documentclass[11pt]{article}\\usepackage{amssymb}\\usepackage{alltt}\\textwidth    = 6.25 in\\textheight  = 9 in\\oddsidemargin = 0.0 in\\evensidemargin = 0.0 in\\topmargin = 0.0 in\\headheight = 0.0 in\\headsep = 0.0 in\\parskip = 0.2in\\parindent = 0.0in\\begin{document}";
latexBottom = "\end{document}";
/*--------------------------------------------------*/

function countHits (dummy)
{
	hitNames[hitCount] = SQL_ROW_DATA[0];
	hitDescs[hitCount] = SQL_ROW_DATA[1];
	hitCount = hitCount + 1;
	return 0;
}

/*--------------------------------------------------*/

function echoQuery (dummy)
{
	colCount = Columns (SQL_ROW_DATA);
	if (doExport)
	{
		cText = classStyle;
		for (cc=0; cc<colCount; cc=cc+1)
		{
			cName  = SQL_COLUMN_NAMES[cc];
			
			inLines = {{""}};
			cr = SQL_ROW_DATA[cc]&&2;
			if (Abs(cr))
			{
				sscanf (cr,"Lines",inLines);
				c = "";
				for (k=0; k<Columns(inLines)-1;k=k+1)
				{
					if (Abs(inLines[k]))
					{
						c = c+inLines[k]+"\\linebreak\n";
					}
					else
					{
						c = c+"\n\n";
					}
				}
				c = c+inLines[k];
				
				cr = "#&$%~_^{}";
				
				for (rcc = 0; rcc < Abs(cr); rcc=rcc+1)
				{
					crp = "\\"+cr[rcc];
					c=c^{{crp,crp}};
				}

				c=c^{{"\\|","$|$"}};
				c=c^{{"\\[","$[$"}};
				c=c^{{"\\]","$]$"}};

				/*cr=cr^{{"\\\\\\\\","$\\backslash$"}};*/
				
				cText  = cText^{{"_INSERT_DATA_FOR_"+cName+"_HERE_",c}};
				sscanf (cName,"Raw",rcc);
			}
			else
			{
				cText  = cText^{{"_INSERT_DATA_FOR_"+cName+"_HERE_","N/A"}};			
			}
		}	
		fprintf (latexFile,"\n", cText);
	}
	else
	{
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
	}
	return 0;
}

/*--------------------------------------------------*/

if (QUERY_FIELD == "Window")
{
	filePath = HYPHY_BASE_DIRECTORY + "Help" + DIRECTORY_SEPARATOR + "Commands"+DIRECTORY_SEPARATOR+"reference.sql";
	OpenWindow (DATABASEWINDOW,{{filePath}},"SCREEN_WIDTH-100;SCREEN_HEIGHT-100;50;50");
	return 0;
}

DoSQL (SQL_OPEN,"reference.sql",dbID);
hitCount        = 0;
hitNames        = {};
hitDescs		= {}; 

if (QUERY_FIELD == "Export")
{
	QUERY_FIELD = "Command";
	doExport = 1;
}
else
{
	doExport = 0;
}

sqlQuery = "select Command, Description from HYPHY_REFERENCE where "+QUERY_FIELD+" glob '*" + QUERY_TERM+"*' order by Command";

DoSQL (dbID,sqlQuery,"return countHits(0);");

if (hitCount == 0)
{
	fprintf (stdout, "\nNo records matched your query\n");
}
else
{
	if (doExport)
	{
		DEFAULT_FILE_SAVE_NAME = "Commands.tex";
		SetDialogPrompt ("Export a LaTex document to:");
		fprintf (PROMPT_FOR_FILE, CLEAR_FILE,latexTop);
		latexFile = LAST_FILE_PATH;
		
	}
	if (hitCount > 1)
	{
		hitRecords = {hitCount,2};
		for (cc=0; cc<hitCount; cc=cc+1)
		{
			hitRecords[cc][0] = hitNames[cc];
			hitRecords[cc][1] = hitDescs[cc];
		}
		
		ChoiceList (hitCount,"Choose a command",0,SKIP_NONE,hitRecords);

		if (hitCount[0]<0)
		{
			DoSQL (SQL_CLOSE,"",dbID);
			return 0;
		}		
		
		for (h=0; h<Columns(hitCount); h=h+1)
		{
			hc = hitCount [h];
			hitNames[0] = hitNames[hc];
			DoSQL (dbID,"select Command, Class, Syntax, Description, Notes, Example, Revised from HYPHY_REFERENCE where Command = '" + hitNames[0]+"' order by Command","return echoQuery(0);");	
		}
		hitDescs    = 0;
		hitRecords  = 0;
	}
	else
	{
		DoSQL (dbID,"select Command, Class, Syntax, Description, Notes, Example, Revised from HYPHY_REFERENCE where Command = '" + hitNames[0]+"' order by Command","return echoQuery(0);");
	}
	hitNames = 0;
	if (doExport)
	{
		fprintf (latexFile,latexBottom);
	}
}

DoSQL (SQL_CLOSE,"",dbID);
