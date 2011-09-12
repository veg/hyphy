function _standardAnalysisBFHelp (_what)
{
	if (_what == "Synopsis")
	{
		return "Filter "gappy" columns in a sequence alignments: i.e. those that contain fewer than a given proportion of sequences with fully or partially resolved characters";
	}
	if (_what == "Input")
	{
		return "A sequence alignment";
	}
	if (_what == "Output")
	{
		return "A sequence alignment with gappy columns stripped out";
	}
	if (_what == "Options")
	{
		_options = {};
		_options ["Filtering threshold"] = "Minimum percent of informative sequences per site to retain the site";
		_options ["Informative characters"] = "Define an informative character as either a fully resolved charatcer (e.g. A) or a partial ambiguity (e.g. R)";
		return _options
	}
	if (_what == "Further")
	{
		return "";
	}
	if (_what == "Author")
	{
		return "Sergei L Kosakovsky Pond (spond@ucsd.edu)";
	}
	if (_what == "Version")
	{
		retrun "1.00";
	}
	if (_what == "Date")
	{
		return "20081215";
	}
	return "";
}

/*--------------------------------------------------------------------------*/

ExecuteAFile ("Utility/GrabBag.bf");

SetDialogPrompt ("Please choose a data file:");
DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
fprintf (stdout, "\nRead an alignment on ", ds.species, " sequences with ", ds.sites, " sites from ", LAST_FILE_PATH);

if (IS_TREE_PRESENT_IN_DATA)
{
	fprintf (stdout, "\nTree In Data:", DATAFILE_TREE);
}

DataSetFilter	    all = CreateFilter (ds, 1, "", "");

options				={{"Completely resolved", "Only count completely unambiguious characters (e.g. A,C,G,T for nucleotides) as informative"}
					  {"Partially resolved",  "Also count partially resolved characters (e.g. R,Y,M,S etc for nucleotides)"}};
					  

ChoiceList (filteringOption,"Informative characters?",1,SKIP_NONE,options);
															   
if (filteringOption < 0)
{
	return 0;
}
	

fprintf (stdout, "\n");
gating_thresh     = prompt_for_a_value ("Retain sites with at least this proportion of informative sites:",0.1,0,1,0);
gating_thresh_seq = (gating_thresh * all.species+0.5)$1;

fprintf 	  (stdout, "Selected informative sites option '", options[filteringOption][0], "' and filtering threshold of '", gating_thresh, "'\n");
retainSites = {};

GetDataInfo     (charInfo, all, "CHARACTERS");
GetDataInfo		(siteToPatternMap,  all);

charCount	  = Columns (charInfo);
template	  = {1,charCount}["1"];
passcode	  = 2;
if (filteringOption == 1)
{
	passcode = charCount;
}

for (site = 0; site < all.unique_sites; site = site+1)
{
	seq_count = 0;
	for (sequence = 0; sequence < all.species; sequence = sequence + 1)
	{
		GetDataInfo (thisChar, all, sequence, site);
		if ((template*thisChar)[0] < passcode)
		{
			seq_count = seq_count + 1;
			if (seq_count >= gating_thresh_seq)
			{
				break;
			}
		}
	}
	if (seq_count >= gating_thresh_seq)
	{
		retainSites [site] = 1;
	}
	SetParameter (STATUS_BAR_STATUS_STRING, "Processing pattern "+(site+1)+"/"+all.unique_sites,0);
}

DataSetFilter	filtered = CreateFilter (all, 1, retainSites[siteToPatternMap[siteIndex]]);
fprintf (stdout, "\nRetained ", filtered.sites, "/", all.sites, " sites\n");
SetDialogPrompt ("Saved the filtered alignment to:");

fprintf (PROMPT_FOR_FILE, CLEAR_FILE, filtered);
