RequireVersion ("2.5.80");
LoadFunctionLibrary("libv3/all-terms.bf"); // must be loaded before CF3x4
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");

utility.SetEnvVariable ("USE_MEMORY_SAVING_DATA_STRUCTURES", 1e8);
utility.SetEnvVariable ("STRICT_ALIGNMENT_VALIDATION_MODE", FALSE);


/*------------------------------------------------------------------------------*/

cln.analysis_description = {
                               terms.io.info : 
"Read a sequence alignment and 'normalize' it, by cleaning  
sequence identifiers, removing duplicates and or gaps, etc",
                               terms.io.version : "0.1",
                               terms.io.reference : "N/A",
                               terms.io.authors : "Sergei L Kosakovsky Pond",
                               terms.io.contact : "spond@temple.edu",
                               terms.io.requirements : "an alignment of nucleotide (coding) sequences",
                               terms.settings: {}
                              };


io.DisplayAnalysisBanner ( cln.analysis_description );


KeywordArgument ("code",        "Genetic code to use", "Universal", "Choose Genetic Code");
KeywordArgument ("alignment",   "Sequence alignment to clean");


cln.code_info = alignments.LoadGeneticCode (None);
cln.alignment_info = alignments.ReadNucleotideDataSet ("cln.sequences", null);
DataSetFilter	    all64 = CreateFilter (cln.sequences, 3, "", "");

_Genetic_Code = cln.code_info[terms.code];

KeywordArgument ("filtering-method",   "How to filter duplicates/gaps?", "No/No");
filteringOption = io.SelectAnOption ({"No/No" :  "Keep all sequences and sites",
                                                                  "No/Yes" :   "Keep all sequences, filter sites with nothing but gaps",
								  	      	    		   	  	  "Yes/No" :  "Filter duplicate sequences but keep all sites",
                                                                  "Yes/Yes" :  "Filter duplicate sequences and sites with nothing but gaps",
                                                                  "Disallow stops" : "Filter duplicate sequences and all sequences that have stop codons"},"Filter duplicates/gaps?");


KeywordArgument ("output", "Write the resulting alignment to");


GetDataInfo (filterDimensions,all64,"CHARACTERS");
filterDimensions = Columns(filterDimensions);

stopCodonTemplate = _Genetic_Code ["_MATRIX_ELEMENT_VALUE_==10"];
nonStopCodonTemplate = _Genetic_Code ["_MATRIX_ELEMENT_VALUE_!=10"];

sequenceNames = {all64.species, 1};
doSomething     = 0;

validID = "[_|a-z|A-Z|0-9]+";

seqNamesList    = {};
seqNamesListAll = {};

for (k=0; k<all64.species; k+= 1) {
	GetString(seqName, all64, k);
	newName 	   = "";
	changed 	   = 0;
	k2    		     = Abs(seqName);
	if (k2 == 0) {
		newName = "Unnamed";
	} else {
	    if (DO_NOT_RENAME_SEQUENCES != TRUE) {
           for (k3 = 0; k3 < k2; k3 += 1) {
                aChar = seqName[k3];
                if ((aChar$validID)[0] < 0) {
                    newName = newName+"_";
                    changed = 1;
                }
                else {
                    newName += aChar;
                }
            }
        } else {
            newName = seqName;
        }
    }

    baseName = newName;
    testName = newName && 1;
    k2 = 2;
    while (seqNamesListAll[testName] > 0) {
        newName = baseName + "_" + k2;
        testName = newName && 1;
        k2 = k2+1;
        changed = 1;
    }

    if (changed) {
       doSomething = 1;
       sequenceNames   [k] = newName;
    }
    else {
        sequenceNames   [k] = seqName;
    }

    seqNamesListAll[testName] = 1;
}


GetInformation (sequenceData,	 all64);
GetDataInfo    (duplicateMapper, all64);

if (onlyFilterSequenceNames != 1) {
	replacementString = "---";
}
else {
	replacementString = "-";
}

notDuplicate	  = {};
duplicateChecker  = {};
haveInfoAtSites	  = {};

all64.unique_sites = Rows (all64.site_freqs) * Columns (all64.site_freqs);

for (sequenceIndex = 0; sequenceIndex < all64.species; sequenceIndex += 1) {

    stopCodonCount     = 0;
    sitesWithDeletions = {1,all64.unique_sites};
    
    COUNT_GAPS_IN_FREQUENCIES = 0;
    
    for (siteIndex = 0; siteIndex < all64.unique_sites; siteIndex += 1) {
    
        GetDataInfo (siteInfo, all64, sequenceIndex, siteIndex);
        
        
        siteInfo1 = stopCodonTemplate*siteInfo;
        siteInfo2 = nonStopCodonTemplate*siteInfo;
    
    
        if (siteInfo1[0]>0 && siteInfo2[0] == 0)    {
            sitesWithDeletions[siteIndex] = 1;
            stopCodonCount  += 1;
        }
        
        if (filteringOption % 2) {
            if (haveInfoAtSites[siteIndex] == 0) {
                if (siteInfo1[0]+siteInfo2[0] > 0) {
                    haveInfoAtSites[siteIndex] = 1;
                }
            }
        }
    }
    
    if (stopCodonCount > 0) {
        if (filterinOption == 4) {
         continue;
        }
        fprintf (stdout, "\nSequence ", sequenceNames[sequenceIndex], ":");
        fprintf (stdout, "\n\t", Format(stopCodonCount,8,0), " stop codons found.");
    
        doSomething		= 1;
        cleanedString		= "";
        seqString			= sequenceData[sequenceIndex];
        cleanedString   * (Abs(seqString)+1);
    
        for (siteIndex = 0; siteIndex < all64.sites; siteIndex += 1) {
            stopCodonCount = duplicateMapper[siteIndex];
            if (sitesWithDeletions[stopCodonCount]) {
                cleanedString * replacementString;
            } else {
                cleanedString * seqString[3*siteIndex][3*siteIndex+2];
            }
        }
        cleanedString * 0;
        sequenceData[sequenceIndex] = cleanedString;
    }
    
    
    
    if (filteringOption >= 2) {
        if (duplicateChecker[sequenceData[sequenceIndex]] == 0) {
            duplicateChecker[sequenceData[sequenceIndex]] = 1;
            notDuplicate[sequenceIndex] = 1;
        } else {
            doSomething   = 1;
        }
    } else {
        notDuplicate[sequenceIndex] = 1;																																														           
    }
}

filterSites = 0;

if (filteringOption%2) {
    filterSites = Abs(haveInfoAtSites)<all64.unique_sites;
	doSomething = doSomething || filterSites;
}

if (!doSomething) {
    fprintf (stdout, "\n\nNo stop codons found\n\n");
}


cln.file_path = io.PromptUserForFilePath ("Save cleaned data to:");

fprintf (cln.file_path, CLEAR_FILE, KEEP_OPEN);
seqLen = Abs (sequenceData[0]);
for (sequenceIndex = 0; sequenceIndex < all64.species; sequenceIndex += 1) {
    if (notDuplicate[sequenceIndex]) {
        if (filterSites) {
            fprintf (cln.file_path, ">", sequenceNames[sequenceIndex], "\n");
            for (s = 0; s < seqLen; s+=3) {
                if (haveInfoAtSites[duplicateMapper[s$3]]) {
                    fprintf (cln.file_path, (sequenceData[sequenceIndex])[s][s+2]);
                }
            }
            fprintf (cln.file_path, "\n\n");
        } else {
            fprintf (cln.file_path, ">", sequenceNames[sequenceIndex], "\n", sequenceData[sequenceIndex], "\n\n");
        }
    }
}
fprintf (cln.file_path, CLOSE_FILE);
