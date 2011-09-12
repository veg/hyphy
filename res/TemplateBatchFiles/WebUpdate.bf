availableFiles = {};
pIndex = "path";
dIndex = "description";

/* -------------------------------------------------------------------------------------- */

function CompareGMTDates (date1, date2)
{
	for (dateComp = 0; dateComp < 5; dateComp = dateComp + 1)
	{
		dateComp2 = date1[dateComp]-date2[dateComp];
		if (dateComp2)
		{
			return dateComp2/Abs(dateComp2);
		}
	}
	return 0;
}

/* -------------------------------------------------------------------------------------- */

function SplitGMTDate (aDate, splitDate&)
{
	chopPoints = aDate || "[^/\\ :]+";
	if (Rows(chopPoints)==10)
	{
		parts = {5,1};
		for (chopPart = 0; chopPart < 10; chopPart = chopPart + 2)
		{
			startOffset = chopPoints[chopPart];
			endOffset	= chopPoints[chopPart+1];
			parts[chopPart/2] = aDate[startOffset][endOffset];
		}
		chopPoints =  parts[0] $ "^[0-9]{4}$";
		if (chopPoints[0]>=0)
		{
			for (chopPart = 1; chopPart < 5; chopPart = chopPart + 1)
			{
				chopPoints =  parts[chopPart] $ "^[0-9][0-9]?$";
				if (chopPoints[0]<0)
				{
					return 0;
				}
			}
			splitDate = {5,1};
			for (chopPart = 0; chopPart < 5; chopPart = chopPart + 1)
			{
				splitDate[chopPart] = 0 + parts[chopPart];
			}
			return 1;
		}
	}
	return 0;
}

/* -------------------------------------------------------------------------------------- */

function SplitInput (fileList)
{
	chopPoints = fileList || "[^\\$]+";
	if  (chopPoints[0]>=0)
	{
		for (entryIndex = 0; entryIndex < Abs(chopPoints); entryIndex = entryIndex + 2)
		{
			startOffset = chopPoints[entryIndex];
			endOffset	= chopPoints[entryIndex+1];
			anEntry 	= fileList[startOffset][endOffset];
			
			splitEntry  = anEntry || "[^\\|]+";
			
			if (Rows(splitEntry) == 8) 
			{
				splitEntryStrings = {4,1};
				for (entryIndex2 = 0; entryIndex2 < 8; entryIndex2 = entryIndex2 + 2)
				{
					startOffset = splitEntry[entryIndex2];
					endOffset	= splitEntry[entryIndex2+1];
					splitEntryStrings[entryIndex2/2] = anEntry[startOffset][endOffset];
				}
				
				entryDate = 0;
				
				if (SplitGMTDate (splitEntryStrings[0],"entryDate"))
				{
					if (CompareGMTDates (lastUpdate, entryDate) < 0)
					{
						newFile = {};
						newFile [pIndex] 		= splitEntryStrings[1];
						newFile [dIndex] = splitEntryStrings[2];
						newFile ["text"]		= ((splitEntryStrings[3]&&1) == "TEXT");
						availableFiles [Abs(availableFiles)] = newFile;
					}
				}
			}
		} 
	}
	return Abs (availableFiles);
}

downloadURLs = {};
downloadURLs [0] = "http://homepage.mac.com/sergeilkp/HyPhy_Distro/HYPHY.dmg";
downloadURLs [1] = "http://homepage.mac.com/sergeilkp/HyPhy_Distro/HYPHY_Win32.zip";
downloadURLs [2] = "http://homepage.mac.com/sergeilkp/HyPhy_Distro/HYPHY_Source.tgz";
downloadURLs [3] = "http://homepage.mac.com/sergeilkp/HyPhy_Distro/HYPHY_UB.dmg";
downloadTargets = {};
downloadTargets [0] = "HYPHY.dmg";
downloadTargets [1] = "HYPHY_Win32.zip";
downloadTargets [2] = "HYPHY_Source.tgz";
downloadTargets [3] = "HYPHY_UB.dmg";

/* -------------------------------------------------------------------------------------- */

fscanf ("last.date","String",lastDate);

if (SplitGMTDate (lastDate, "lastUpdate"))
{
	fprintf (stdout, "\nStarting HyPhy update\n");
	GetString (versionString, HYPHY_VERSION, 2);
	dateComp = versionString$"Macintosh.+";
	if (dateComp[0]>=0)
	{
		versionString = "Macintosh";
		if ((versionString$"Universal")[0]>=0)
		{
			systemKind = 3;		
		}
		else
		{
			systemKind = 0;
		}
	}
	else
	{
		dateComp = versionString$"Windows.+";
		if (dateComp[0]>=0)
		{
			versionString = "Windows";
			systemKind = 1;
		}
		else
		{
			versionString = "Unix";
			systemKind = 2;
		}
	}
	versionString = "http://www.hyphy.org/updates/"+versionString+"/";
	
	GetURL (updateIndex, versionString+"version");
	
	GetString (versionStringS, HYPHY_VERSION, 0);
	
	updateIndex 	= 0+updateIndex;
	versionStringS  = 0+versionStringS;
	
	
	if (versionStringS<updateIndex)
	{
		fprintf (stdout, "\nA newer version of HyPhy(",updateIndex,")  is available for your platform.");
		DEFAULT_FILE_SAVE_NAME = downloadTargets[systemKind];
		SetDialogPrompt ("Download a newer version of HyPhy to:");
		fprintf (PROMPT_FOR_FILE,CLEAR_FILE);
		GetURL (LAST_FILE_PATH,downloadURLs[systemKind],SAVE_TO_FILE);
		fprintf (stdout, "\nA newer version of HyPhy has been downloaded.\nPlease quit HyPhy and install the update.\n");
		return 0;
	}
	else
	{
		fprintf (stdout, "\nYour version of HyPhy is current.\nChecking for individual file updates...\n");	
	}
	
	/* first check the version */
	
	GetURL (updateIndex, versionString+"updates");
	
	
	if (SplitInput (updateIndex))
	{
		for (chopPart = 0; chopPart < Abs(availableFiles); chopPart = chopPart + 1)
		{
			entryData = availableFiles[chopPart];
			filePath = HYPHY_LIB_DIRECTORY+entryData[pIndex]^{{"/"}{DIRECTORY_SEPARATOR}};
			if (entryData["text"])
			{
				fprintf (stdout, "Updating file:", entryData[pIndex], 
							   "\n-------------\n",
							   "\nNotes:", 
							   "\n-----\n", entryData[dIndex], "\n");
				GetURL (theFile,versionString+entryData[pIndex]);
				fprintf (filePath, CLEAR_FILE, theFile);
				fprintf (stdout, "Update Successful\n");
			}
			else
			{
				fprintf (stdout, "Saving file:", entryData[pIndex], 
							   "\n------------\n",
							   "\nNotes:", 
							   "\n-----\n", entryData[dIndex], "\n");
				GetURL (filePath,versionString+availableFiles[pIndex],SAVE_TO_FILE);
				fprintf (stdout, "Download Successful\n");				
			}
		}
	}
	else
	{
		fprintf (stdout, "\nYour distribution seems to be up-to-date\n");
	}
	GetString (gmtString,  TIME_STAMP, 0);
	fprintf ("last.date",CLEAR_FILE,gmtString);
}
else
{
	fprintf (stdout, "\nInternal error : could not read valid date from last.date file\n");
}
