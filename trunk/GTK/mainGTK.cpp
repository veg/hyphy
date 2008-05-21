#include <stdio.h>
#include <gtk/gtk.h>
#include <gdk/gdk.h>
#include <unistd.h>

#include <batchlan.h>

#include "HYDialogs.h"
#include <hy_strings.h>
#include <HYTableWindow.h>
#include <HYCanvas.h>
#include <HYGWindow.h>
#include <HYButton.h>
#include <HYButtonBar.h>
#include <HYTextBox.h>
#include <HYUtils.h>
#include <HYLabel.h>
#include <HYEventTypes.h>
#include <HYSequencePanel.h>
#include <HYTableComponent.h>
#include <HYConsoleWindow.h>
#include <HYTreePanel.h>
#include <HYDataPanel.h>
#include <HYObjectInspector.h>
#include <time.h>

#ifdef   __MP2__
	#include <pthread.h>
#endif

#ifdef	__HYPHYMPI__
	#include "likefunc.h"
	extern _String shortMPIReturn;
	void    mpiOptimizerLoop 			(int, int);
	void 	mpiNormalLoop	 			(int, int, _String*);
	int		_hy_mpi_node_rank;
#endif

//#include <gdk-pixbuf/gdk-pixbuf.h>

						
_List					globalPreferencesList,
						recentFiles,
						recentPaths;
						
extern	_List			availablePostProcessors,
						availableTemplateFiles;

_String					treeDisplayOptions ("TREE_DISPLAY_OPTIONS"),			
						initialDialogPop ("SHOW_DIALOG_AT_STARTUP"),
						windowPositionStr("CONSOLE_WINDOW_COORDS"),
						preserveSlaveNodeState ("PRESERVE_SLAVE_NODE_STATE"),
						recentFilesList	 ("RECENT_FILES_LIST"),
						selectedFileName,
						savedFileName,
						baseArgDir;
						
extern	bool			skipWarningMessages;

bool					showDialogAtStartup = true,
						terminateExecution  = false,
						isSuspended 		= false, 
						hasTemplates 		= false, 
						isRerunAvailable 	= false, 
						updateTimer 		= false, 
						addToRecent			= true,
						echoPaused  		= false;
						
extern  _String			optimizationPrecision, 
						skipOmissions,
						maximumIterationsPerVariable,
						dataFilePrintFormat,
						dataFileDefaultWidth, 
						dataFileGapWidth, 
						useInitialDistanceGuess, 
						globalStartingPoint,
						categorySimulationMethod,
						likefuncOutput,
						printDigitsSpec, 
						randomSeed, 
						VerbosityLevelString,
						objectInspectorTitle,
						baseDirectory,
						*argFileName,
						errorFileName,
						messageFileName;
				

PangoContext*			screenPContext;
GdkPixbuf*				tablePDMenuIcon,
						*redButtonIcon,
						*yellowButtonIcon,
						*orangeButtonIcon,
						*greenButtonIcon;

GdkCursor				*hSizeCursor,
						*pickUpCursor,
						*dropOffCursor;
						
clock_t					timerStart = 0, 
						lastTimer  = 0;

long					leftRectCoord = 0, 
						topRectCoord = 0,
						bottomRectCoord = 0, 
						rightRectCoord = 0;
						
extern	long			systemCPUCount;

double					fontConversionFactor = 1.0;


void					HandleGlobalQueueEvent (void);

_HYConsoleWindow		*hyphyConsoleWindow = nil;

_ExecutionList			ex;

//____________________________________________________________________________________________

void			AddStringToRecentMenu (_String&, _String&);
long  			SelectATemplate 	  (void);
void			RunTemplate 		  (long);
void			ReadInTemplateFiles	  (void);
void			RunStandardAnalyses   (void);
void  			SetPreferences 		  (void);
void			ReadInPostFiles		  (void);

//_________________________________________________________________________
void	updateTimerF (_String& rec, clock_t time_diff)
{
	long secs, mins, hrs;
	secs = time_diff/CLOCKS_PER_SEC;
	mins = secs/60;
	hrs = mins/60;
	mins = mins%60;
	secs = secs%60;
	if (hrs<10)
		rec = _String('0')&hrs;
	else
		rec = _String(hrs);
	rec = rec &':';
	if (mins<10)
		rec = rec&_String('0')&mins;
	else
		rec = rec&_String(mins);
	rec = rec &':';
	if (secs<10)
		rec = rec&_String('0')&secs;
	else
		rec = rec&_String(secs);	
}



//_________________________________________________________________________
void	yieldCPUTime (void)
{
	while (gtk_events_pending ())
		gtk_main_iteration();
}

//____________________________________________________________________________________________

long  SelectATemplate (void)
{	
	_SimpleList std, vc,selection;
	long		i;
	std<<2;
	std<<1;
	for (i=0; i<availableTemplateFiles.lLength; i++)
		vc<<i;
	
	return HandleHierListSelection (availableTemplateFiles, std, vc, "Select a standard analysis to run",selection,1);
}

//________________________________________________________

void	RunTemplate (long idx)
{
	PurgeAll(windowPtrs.lLength==0);
	_String pathName = baseDirectory&"TemplateBatchFiles/";
	pathNames&& &pathName;
	pathName = pathName&*(_String*)(*(_List*)availableTemplateFiles(idx))(2);
	if (!argFileName)
		argFileName = new _String(pathName);
	else
		*argFileName = pathName;
	//ReadBatchFile (pathName,ex);
	selectedFileName = *argFileName;
	ExecuteBatchFile();
}

//________________________________________________________

void		RunStandardAnalyses (void)
{
	long menuChoice=SelectATemplate();
	if (menuChoice >= 0)
		RunTemplate(menuChoice);
}



//_________________________________________________________________________
bool  OpenBatchFile (bool openOrNot, _String* dL)
{
	PurgeAll(windowPtrs.lLength==1);
	
	_String  dialogPrompt = "Open Batch File",
			 fName, 
			 pathName;
			 
	if (openOrNot)
	{
		PopUpFileDialog(dialogPrompt,dL);
		fName = *argFileName;
	}
	else
		fName = selectedFileName;
	
	if (fName.sLength)
	{
		selectedFileName = fName;
		PushFilePath (fName);
		//ExecuteAFile ();
	
		return  true;
	}
	
	return false;
}

//__________________________________________________________________

gboolean GlobalQueueTimer (Ptr* userData)
{
	if (GlobalGUIEventQueue.lLength)
		HandleGlobalQueueEvent ();
	//if (updateTimer)
	//	SendMessage ((HWND)hyphyConsoleWindow->GetOSWindowData(), UPDATE_TIMER, 0, 0);
	
	return true;
}


//__________________________________________________________________

void SetUpStatusBarStuff (GtkWidget* aWindow)
{
	_String			   fName = baseDirectory & "GTKResources/striped.xpm";
	statusBarLayout			 = pango_layout_new (screenPContext);
	statusBarFontDesc		 = pango_font_description_new ();
	stripedFill				 = gdk_pixmap_create_from_xpm (GDK_DRAWABLE(aWindow->window), NULL, NULL, fName.sData);
	stripedFillGC			 = gdk_gc_new (GDK_DRAWABLE(aWindow->window));
	if (stripedFill)
	{
		gdk_gc_set_fill (stripedFillGC,GDK_TILED);
		gdk_gc_set_tile	(stripedFillGC,stripedFill);
	}
	else
	{
		printf ("Failed to load a status bar .xpm from %s\n", fName.sData);
	}
	
	gdk_gc_set_line_attributes		  (stripedFillGC, 1, GDK_LINE_SOLID, GDK_CAP_NOT_LAST, GDK_JOIN_MITER);
	GdkColor saveFG = {0,0,0,0};
	gdk_gc_set_foreground			  (stripedFillGC, &saveFG);

	pango_font_description_set_family (statusBarFontDesc, statusBarFont.face.sData);
	pango_font_description_set_style  (statusBarFontDesc, (statusBarFont.style & HY_FONT_ITALIC) ? PANGO_STYLE_ITALIC : PANGO_STYLE_NORMAL);
	pango_font_description_set_weight (statusBarFontDesc, (statusBarFont.style & HY_FONT_BOLD) ? PANGO_WEIGHT_BOLD : PANGO_WEIGHT_NORMAL);
	pango_font_description_set_size   (statusBarFontDesc, statusBarFont.size*PANGO_SCALE);
	pango_layout_set_font_description (statusBarLayout, statusBarFontDesc ); // ref ?
	pango_layout_set_width			  (statusBarLayout, -1);
	
	redButtonIcon = (GdkPixbuf*)ProcureIconResource(4000);
	yellowButtonIcon = (GdkPixbuf*)ProcureIconResource(4001);
	greenButtonIcon = (GdkPixbuf*)ProcureIconResource(4002);
	orangeButtonIcon = (GdkPixbuf*)ProcureIconResource(4003);
	
}
//____________________________________________________________________________________________

void  ReadPreferences (void)
{
	_String optionList, 
			comma(",");
			
	_List   fonts;
	
	optionList = baseDirectory & ".hyphyprefs";
	FILE * prefsFile = fopen (optionList.sData,"rb");
	
	AddItemToPreferences (1|8,-1,"Console Font Settings","Choose font face, size and style used for displaying text in the console.","",nil);

	GenerateFontList	 (fonts);	
	
	AddItemToPreferences (0,PREFITEM_POPUP,"Font Face","Select the font used to display text in the console window.","Monaco",&fonts);
	fonts.Clear();
	optionList = "8,9,10,12,14,18";
	_List* options = optionList.Tokenize (comma);
	AddItemToPreferences (0,PREFITEM_POPUP,"Font Size","Select font size used to display text in the console window.","9",options);
	DeleteObject(options);
	optionList = "Plain,Bold,Italic";
	options = optionList.Tokenize (comma);
	AddItemToPreferences (0,PREFITEM_POPUP,"Font Style","Select font style used to display text in the console window.","Plain",options);
	DeleteObject(options);
	AddItemToPreferences (1|8,-1,"Optimization Settings","Options affecting the optimization algorithm.","",nil);
	AddItemToPreferences (0,PREFITEM_TEXTBOX,"Precision","Desired precision(absolute error) in ln-likelihood value. Settings between 0.1 and 0.000000001 are recommended.","0.001",nil);	
	optionList = "Low,Normal,High,Very High";
	options = optionList.Tokenize (comma);
	AddItemToPreferences (0,PREFITEM_POPUP,"Persistence","Controls the number iterations the optimization algorithm will perform before it terminates if the desired precision is not met.","Normal",options);	
	DeleteObject(options);
	optionList = "Do not use distances,Use distances";
	options = optionList.Tokenize (comma);
	AddItemToPreferences (0,PREFITEM_POPUP,"Initial Guess","Determines whether distance methods are to be used to obtain intial parameter value guesses. Applies only to nuceleotide models.","Use distances",options);	
	DeleteObject(options);
	AddItemToPreferences (0,PREFITEM_TEXTBOX,"Starting Value","Sets starting values for parameters for optimization routines. If starting values are obtained by distance methods, this option is ignored.","0.1",nil);	
	AddItemToPreferences (1|8,-1,"Data Read/Write Settings","Options affecting sequence data files reading and writing.","",nil);
	optionList = "Skip Deletions,Keep Deletions";
	options = optionList.Tokenize (comma);
	AddItemToPreferences (0,PREFITEM_POPUP,"Deletions","Choose \"Keep Deletions\" to retain deletions (as ambiguities) for analyses. \"Skip Deletions\" filters deletions out as the data is read.","Keep Deletions",options);	
	DeleteObject(options);
	optionList = "# sequential,# interleaved,PHYLIP Sequential,PHYLIP Interleaved,NEXUS sequential with labels,NEXUS interleaved with labels,NEXUS sequential without labels,NEXUS interleaved without labels";
	options = optionList.Tokenize (comma);
	AddItemToPreferences (0,PREFITEM_POPUP,"Output format","Choose the default file format for data filters output to files via fprintf.","NEXUS sequential without labels",options);	
	DeleteObject(options);
	AddItemToPreferences (0,PREFITEM_TEXTBOX,"Line width","This options sets how many characters will be printed per line for data filters output to files via fprintf. Only affects interleaved formats.","50",nil);	
	AddItemToPreferences (0,PREFITEM_TEXTBOX,"Gap width","This options sets how many characters will be printed per cluster (clusters are separated by spaces) for data filters output to files via fprintf. Only affects interleaved non-NEXUS formats.","10",nil);	

	AddItemToPreferences (1|8,-1,"Simulation Options","Options affecting bootstrapping algorithms.","",nil);
	optionList = "Discrete Distribution,Continuous Distribution";
	options = optionList.Tokenize (comma);
	AddItemToPreferences (0,PREFITEM_POPUP,"Heterogeneity Simulation","When bootstrapping models with heterogeneous rates, determines whether rate classes are drawn from the continuous (e.g. gamma) distribution or it's discrete approximation.","Continuous Distribution",options);	
	DeleteObject(options);
	AddItemToPreferences (0,PREFITEM_TEXTBOX,"Random seed","Set this parameter to -1 to have HYPHY seed random generator anew every time the program is run. A positive value defines the seed to be used instead. Changes will take effect when HYPHY is restarted.","-1",nil);	
	AddItemToPreferences (1|8,-1,"Miscellaneous Options","Variuos, primarily formatting, options.","",nil);
	optionList = "Function value only,Complete report as list,Tree with branch lengths,Parameters and Constraints,Batch Language Statement,Batch Language Statement with Trees";
	options = optionList.Tokenize (comma);
	AddItemToPreferences (0,PREFITEM_POPUP,"Likelihood Display","Various ways to display likelihood function and parameters","Tree with branch lengths",options);	
	DeleteObject(options);
	optionList = "Short,Normal,Long,Maximally Long";
	options = optionList.Tokenize (comma);
	AddItemToPreferences (0,PREFITEM_POPUP,"Number Format","Determines how many significant digits are displayed when printing numbers via fprintf.","Normal",options);	
	DeleteObject(options);
	optionList = "No auto display,Auto display single tree,Auto display all trees";
	options = optionList.Tokenize (comma);
	AddItemToPreferences (0,PREFITEM_POPUP,"Tree Display","Should HY-PHY automatically open graphical tree windows upon completion of an analysis.","Auto display single tree",options);	
	DeleteObject(options);
	optionList = "Silent,Verbose";
	options = optionList.Tokenize (comma);
	AddItemToPreferences (0,PREFITEM_POPUP,"Optimization Progress","Triggers the optimization functions to print out progress lines while obtaining MLEs","Silent",options);	
	DeleteObject(options);
	optionList = "Yes,No";
	options = optionList.Tokenize (comma);
	AddItemToPreferences (0,PREFITEM_POPUP,"Startup Dialog","Display an action dialog when HyPhy starts up","Yes",options);	
	DeleteObject(options);
	optionList = "Yes,No";
	options = optionList.Tokenize (comma);
	AddItemToPreferences (0,PREFITEM_POPUP,"Automove console","Automatically move and resize console window when a data panel is opened","Yes",options);	
	DeleteObject(options);
	if (prefsFile)
	{
		long j;
		_String fileContents(prefsFile);
		fclose (prefsFile);
		_List terms, * availNames = (_List*)globalPreferencesList(1), 
					 * availValues = (_List*)globalPreferencesList(4);
		_ElementaryCommand::ExtractConditions(fileContents,0,terms);
		for (long prefFolderID=0;prefFolderID<terms.lLength;prefFolderID++)
		{
			_String* thisTerm = (_String*)terms(prefFolderID);
			_List	 theTerms;
			_ElementaryCommand::ExtractConditions(*thisTerm,0,theTerms,'=');
			if (theTerms.lLength == 2)
			{
				_String * prefID = (_String*)theTerms.lData[0];
				j = availNames->Find (prefID);
				if (j>=0)
					*((_String*)availValues->lData[j]) = *(_String*)theTerms.lData[1];
				else
				{
					if (prefID->Equal (&windowPositionStr))
					{
						_List*	rectSizes = ((_String*)theTerms.lData[1])->Tokenize (",");
						if (rectSizes->lLength==4)
						{
							leftRectCoord   = ((_String*)(*rectSizes)(0))->toNum();
							topRectCoord    = ((_String*)(*rectSizes)(1))->toNum();
							rightRectCoord  = ((_String*)(*rectSizes)(2))->toNum();
							bottomRectCoord = ((_String*)(*rectSizes)(3))->toNum();
							if (leftRectCoord<0)
								leftRectCoord = 0;
							if (rightRectCoord<leftRectCoord)
								rightRectCoord = leftRectCoord+50;
							if (topRectCoord<0)
								topRectCoord = 0;
							if (bottomRectCoord<topRectCoord)
								bottomRectCoord = topRectCoord+50;
						}
						DeleteObject (rectSizes);						
					}
					else
						if (prefID->Equal (&recentFilesList))
						{
							_List* filePaths = ((_String*)theTerms.lData[1])->Tokenize (",");
							for (long idx = 0; idx<filePaths->lLength; idx+=2)
							{
								_String * sf = (_String*)(*filePaths)(idx),
										* sp = (_String*)(*filePaths)(idx+1);
										
								sp->StripQuotes();
								sf->StripQuotes();									
								AddStringToRecentMenu(*sf,*sp);
							}
							DeleteObject (filePaths);
						}
				}
			}
		}
	}
	
	long	aS = ((_String*)((_List*)globalPreferencesList.lData[4])->lData[16])->toNum();
	if (aS>=0)
		init_genrand (aS);
		
	showDialogAtStartup = *(((_String*)((_List*)globalPreferencesList.lData[4])->lData[22]))==_String("Yes");
	doAutoConsoleMove   = *(((_String*)((_List*)globalPreferencesList.lData[4])->lData[23]))==_String("Yes");
}

//____________________________________________________________________________________________

void  ApplyPreferences (void)
{
	_List	*pfValues = (_List*)globalPreferencesList.lData[4], *t;
	setParameter (optimizationPrecision, ((_String*)pfValues->lData[5])->toNum());
	t = (_List*)((_List*)globalPreferencesList.lData[5])->lData[6];
	long	f = t->Find (((_String*)pfValues->lData[6]));
	switch (f)
	{
		case 0:
			f = 200;
			break;
		case 1:
			f = -1;
			break;
		case 2:
			f = 2000;
			break;
		case 3:
			f = 50000;
			break;
	}
	setParameter (globalStartingPoint, ((_String*)pfValues->lData[8])->toNum());
	setParameter (useInitialDistanceGuess, (_Parameter)(*((_String*)pfValues->lData[7])==_String("Use distances")));
	
	setParameter (skipOmissions, (_Parameter)(*((_String*)pfValues->lData[10])==_String("Skip Deletions")));
	if (f>0)
		setParameter (maximumIterationsPerVariable, (_Parameter)f);
	t = (_List*)((_List*)globalPreferencesList.lData[5])->lData[11];
	f = t->Find (((_String*)pfValues->lData[11]));
	if (f<0) f = 0;
	setParameter (dataFilePrintFormat, (_Parameter)f);
	setParameter (dataFileDefaultWidth , (long)(((_String*)pfValues->lData[12])->toNum()));
	setParameter (dataFileGapWidth , (long)(((_String*)pfValues->lData[13])->toNum()));

	setParameter (categorySimulationMethod, (_Parameter)((*((_String*)pfValues->lData[15])!=_String("Discrete Distribution")))+1.0);

	f = ((_String*)pfValues->lData[16])->toNum();
	if (f>=0)
	{
		setParameter (randomSeed, f);
	}

	t = (_List*)((_List*)globalPreferencesList.lData[5])->lData[18];
	f = t->Find (((_String*)pfValues->lData[18]));
	if (f<0) f = 0;
	setParameter (likefuncOutput, (_Parameter)f);
	
	t = (_List*)((_List*)globalPreferencesList.lData[5])->lData[19];
	f = t->Find (((_String*)pfValues->lData[19]));
	if (f<0) f = 1;
	switch (f)
	{
		case 0:
			printDigits = 5;
			break;
		case 2: 
			printDigits = 12;
			break;
		case 3:
			printDigits = 15;
			break;
	}
	setParameter (printDigitsSpec,printDigits);
	t = (_List*)((_List*)globalPreferencesList.lData[5])->lData[20];
	f = t->Find (((_String*)pfValues->lData[20]));
	setParameter (treeDisplayOptions, f);
	t = (_List*)((_List*)globalPreferencesList.lData[5])->lData[21];
	verbosityLevel = 5*t->Find (((_String*)pfValues->lData[21]));
	setParameter (VerbosityLevelString,verbosityLevel);	

	doAutoConsoleMove   = *(((_String*)((_List*)globalPreferencesList.lData[4])->lData[23]))==_String("Yes");
}

//____________________________________________________________________________________________

void  WritePreferences (void)
{
	_String  prefFileName = baseDirectory & ".hyphyprefs";
	FILE * prefFile = fopen (prefFileName.sData, "wb");
	if (prefFile)
	{	
		_String spoolResult((unsigned long)1, true);
		_SimpleList* pCodes = (_SimpleList*)globalPreferencesList(0);
		long i;
		
		_List *pNames = (_List*)globalPreferencesList(1),
			  *pValues = (_List*)globalPreferencesList(4);
		for (i = 0; i<pCodes->lLength; i++)
		{
			if (pCodes->lData[i]>=8) continue;
			spoolResult<<(_String*)pNames->lData[i];
			spoolResult<<'=';
			spoolResult<<(_String*)pValues->lData[i];
			spoolResult<<';';
		}
		
		_HYRect	   wr = hyphyConsoleWindow->GetWindowRect();
		// convert the coordinate info to a comma separated 
		// string and store it as a STR resource
		spoolResult<<windowPositionStr;
		spoolResult<<'=';
		spoolResult<<_String ((long)wr.left);
		spoolResult<<',';
		spoolResult<<_String ((long)wr.top);
		spoolResult<<',';
		spoolResult<<_String ((long)wr.right);
		spoolResult<<',';
		spoolResult<<_String ((long)wr.bottom);
		spoolResult<<';';
		if (recentFiles.lLength)
		{
			spoolResult << recentFilesList;
			spoolResult << '=';
			for (long k=0; k<recentFiles.lLength; k++)
			{
				if (k)
					spoolResult << ',';
				spoolResult << '"';
				spoolResult << (_String*)recentFiles(k);
				spoolResult << '"';
				spoolResult << ',';
				spoolResult << '"';
				spoolResult << (_String*)recentPaths(k);
				spoolResult << '"';
			}
			spoolResult << ';';
		}
		spoolResult.Finalize();
		fwrite (spoolResult.sData, spoolResult.sLength, 1, prefFile);
		fclose (prefFile);
	}
}

//_________________________________________________________________________

void	ShowObjectInspector (void)
{
	long f = FindWindowByName (objectInspectorTitle);
	if (f>=0)
	{
		gtk_window_present(GTK_WINDOW(windowPtrs(f)));	
	}
	else
	{
		_HYObjectInspector* newOI = new _HYObjectInspector ();
		newOI->BringToFront		  ( );
	}
}

//_________________________________________________________________________

void displayAbout (bool splash)
{
	// TBI
}

//_________________________________________________________________________
void   ExecuteAFile (bool spawnEXL)
{
	
	timerStart = clock();
	lastTimer = timerStart;
	updateTimer = true;		

	ToggleAnalysisMenu (true);
	
	_String justTheName = selectedFileName.Cut(selectedFileName.FindBackwards ('/',0,-1)+1,-1);
	AddStringToRecentMenu(justTheName, selectedFileName);
	terminateExecution = false;
	
	if (!spawnEXL)
		ex.Clear();
	SetStatusLine (selectedFileName.Cut (selectedFileName.FindBackwards("/",0,-1)+1,-1), "Loading","00:00:00");
	ApplyPreferences();
	if (!spawnEXL)
	{
		ReadBatchFile (selectedFileName,ex);
		terminateExecution = false;
		skipWarningMessages = false;
		ex.Execute();
	}
	else
	{
		_ExecutionList elist;
		ReadBatchFile (selectedFileName,elist);
		terminateExecution = false;
		skipWarningMessages = false;
		elist.Execute();
	}
	
	setParameter  	  (VerbosityLevelString, 0.0);
	SetStatusLine 	  ("Finished");
	SetStatusBarValue (-1,1,1);
	BufferToConsole ("\n");
	savedFileName = selectedFileName;
	selectedFileName = empty;
	updateTimer = false;

	ToggleAnalysisMenu (false);
	for (long i=0; i<availablePostProcessors.lLength; i++)
	{
		GtkWidget * mItem = gtk_item_factory_get_item_by_action(hyphyConsoleWindow->menu_items, 1000+i);
		if (mItem)
		{
			bool onOff = true;
			_String* condition = (_String*)(*(_List*)availablePostProcessors(i))(2);
			if (condition->sLength)
			{
				_Formula condCheck (*condition,nil);
				_PMathObj condCheckRes = condCheck.Compute();
				if (!condCheckRes || condCheckRes->Value()<.5)
					onOff = false;
			}
			gtk_widget_set_sensitive(mItem, onOff);
		}
	}
	
	if (!spawnEXL)
		ex.ResetFormulae();
}
//_________________________________________________________________________
bool ExecuteBatchFile (void)
{
	ExecuteAFile (false);
	return true;
}

#ifdef __HYPHYMPI__

//__________________________________________________________________________________
void mpiNormalLoop    (int rank, int size, _String & baseDir)
{
	long		 senderID = 0;
	
	ReportWarning ("Entered mpiNormalLoop");

	_String* theMessage = MPIRecvString (-1,senderID),
			* resStr	= nil,
			css("_CONTEXT_SWITCH_MPIPARTITIONS_");
				
	while (theMessage->sLength)
	{
		setParameter    (mpiNodeID, (_Parameter)rank);
		setParameter	(mpiNodeCount, (_Parameter)size);
		//ReportWarning (*theMessage);
		DeleteObject (resStr);
		resStr = nil;
		if (theMessage->Equal (&css) )
		{
			mpiPartitionOptimizer = true;
			ReportWarning ("Switched to mpiOptimizer loop");
			MPISendString(css,senderID);
			mpiOptimizerLoop (rank,size);
			ReportWarning ("Returned from mpiOptimizer loop");
			mpiPartitionOptimizer = false;
			pathNames && & baseDir;
		}
		else
		{
			if (theMessage->beginswith ("#NEXUS"))
			{
				_String		msgCopy (*theMessage);
				ReportWarning ("Received a function to optimize");
				ReadDataSetFile (nil,true,theMessage);
				ReportWarning ("Done with the optimization");
				_Variable*  lfName = nil;
				long		f = LocateVarByName (lf2SendBack);
				if (f>=0)
					lfName = FetchVar(f);
										
				if (!(lfName&&(lfName->ObjectClass()==STRING)))
				{
					_String errMsg ("Malformed MPI likelihood function optimization request - missing LF name to return.\n\n\n");
					errMsg = errMsg & msgCopy;
					FlagError (errMsg);
					break;
				}
				
				f = likeFuncNamesList.Find (((_FString*)lfName->Compute())->theString);
				if (f<0)
				{
					_String errMsg ("Malformed MPI likelihood function optimization request - invalid LF name to return.\n\n\n");
					errMsg = errMsg & msgCopy;
					FlagError (errMsg);
					break;				
				}
				_Parameter pv;
				checkParameter (shortMPIReturn, pv ,0);
				resStr = new _String (1024L,true);
				checkPointer (resStr);
				((_LikelihoodFunction*)likeFuncList (f))->SerializeLF(*resStr,pv>0.5?5:2);
				resStr->Finalize();
			}
			else
			{
				_ExecutionList exL (*theMessage);
				/*printf ("Received:\n %s\n", ((_String*)exL.toStr())->sData);*/
				_PMathObj res = exL.Execute();
				resStr = res?(_String*)res->toStr():new _String ("0");
			}
				
			checkPointer (resStr);
			DeleteObject (theMessage);
			MPISendString(*resStr,senderID);

			_Parameter 	   	keepState = 0.0;
			checkParameter  (preserveSlaveNodeState, keepState, 0.0);
			
			if (keepState < 0.5)
			{
				PurgeAll (true);
				pathNames && & baseDir;		
			}
		}
		theMessage = MPIRecvString (-1,senderID);		
	}
	/*MPISendString(empty,senderID);*/
	DeleteObject (resStr);
	DeleteObject (theMessage);	
}

//__________________________________________________________________________________
void mpiOptimizerLoop (int rank, int size)
{
	long		 senderID = 0;
			
	ReportWarning (_String ("MPI Node:") & (long)rank & " is ready for MPIParallelOptimizer tasks");
				
	if (mpiPartitionOptimizer)
		ReportWarning (_String("MPI Partitions mode"));
	
	//printf ("Node %d waiting for a string\n", rank);
	_String* theMessage = MPIRecvString (-1,senderID);
	while (theMessage->sLength)
	{
		if (theMessage->beginswith ("#NEXUS"))
		{
			ReadDataSetFile (nil,true,theMessage);
			if (likeFuncNamesList.lLength!=1)
			{
				_String errMsg ("Malformed MPI likelihood function paraller optimizer startup command. No valid LF has been defined.n\n\n");
				FlagError (errMsg);
				break;						
			}
			
			// send back the list of independent variables
			
			_LikelihoodFunction * theLF = (_LikelihoodFunction*)likeFuncList (0);
			
			if (mpiParallelOptimizer && theLF->GetCategoryVars().lLength)
			{
				_String errMsg ("Likelihood functions spawned off to slave MPI nodes can't have category variables.n\n\n");
				FlagError (errMsg);
				break;						
			}
			
			_SimpleList* ivl = & theLF->GetIndependentVars();
			
			_String		 variableSpec (128L, true);
			
			(variableSpec) << LocateVar(ivl->lData[0])->GetName();

			for (long kk = 1; kk < ivl->lLength; kk++)
			{
				(variableSpec) << ';';
				(variableSpec) << LocateVar(ivl->lData[kk])->GetName();	
			}
			
			ReportWarning 		  (variableSpec);
			MPISendString		  (variableSpec,senderID);
			theLF->PrepareToCompute();
			theLF->MPI_LF_Compute (senderID, mpiPartitionOptimizer);
			theLF->DoneComputing();
			PurgeAll (true);
		}
		DeleteObject (theMessage);
		theMessage = MPIRecvString (-1,senderID);
	}	
	DeleteObject (theMessage);		
}

#endif


//__________________________________________________________________

int main( int   argc, char *argv[] )
{

	#ifdef	__HYPHYMPI__
		  int 		   rank, 
		  			   size;
		  			   			   			 
		  MPI_Init	   (&argc, &argv);
		  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		  MPI_Comm_size(MPI_COMM_WORLD, &size);
		  
		  _hy_mpi_node_rank = rank;
		  
		  setParameter  (mpiNodeID, (_Parameter)rank);
		  setParameter	(mpiNodeCount, (_Parameter)size);
		  
		  if (rank == 0)
	#endif
   
	gtk_init (&argc, &argv);

	/* set up globals */
	
	char curWd[4096];
	getcwd (curWd,4096);

	_String baseDir (curWd);
	baseDir=baseDir&'/';

	pathNames&& &baseDir;
	baseDirectory = baseDir;
	for (long i=1; i<argc;i++)
	{
		_String thisArg (argv[i]);
		if (thisArg.beginswith ("BASEPATH="))
		{
			baseArgDir = thisArg.Cut(9,-1);
			if (baseArgDir.sLength)
			{
				if (baseArgDir.sData[baseArgDir.sLength-1]!='/')
					baseArgDir = baseArgDir&"/";
					
				baseDirectory = baseArgDir;
			}
		}
		else
			if (thisArg.beginswith ("USEPATH="))
			{
				baseArgDir 			= thisArg.Cut(8,-1);
				errorFileName 		= baseArgDir & errorFileName;
				messageFileName 	= baseArgDir & messageFileName;
				pathNames.Delete 	(0);
				pathNames&& 		&baseDir;
			}
			else
				if (thisArg.beginswith ("CPU="))
				{
					#ifdef __MP__
					_String cpus = thisArg.Cut(4,-1);
					systemCPUCount = cpus.toNum();
					if (systemCPUCount<1)
						systemCPUCount = 1;
					#ifdef __MP2__
						pthread_setconcurrency (systemCPUCount+1);
					#endif
					#endif
				}
				#ifdef	__HYPHYMPI__
					else
						if (thisArg == _String("MPIOPTIMIZER"))
						{
							mpiParallelOptimizer = true;
							setParameter	(mpiNodeCount, 0.0);
						}
						else
							if (thisArg == _String("MPIPARTITIONS"))
							{
								mpiPartitionOptimizer = true;
								setParameter	(mpiNodeCount, 0.0);
							}
				#endif
	}
	
	#ifdef	__HYPHYMPI__
	if (rank == 0)
	#endif
	{
		baseDir = baseDirectory & "GTKResources";
		_List scanRes;
		ScanDirectoryForFileNames(baseDir,scanRes,false);
		if (scanRes.lLength == 0)
		{
			GtkWidget * noRez = gtk_message_dialog_new (NULL, GTK_DIALOG_MODAL, GTK_MESSAGE_ERROR, GTK_BUTTONS_OK, "HYPHY_GTK was unable to find a required GTKResources directory in %s. Please use BASEPATH= command line option to specify where the installation directory of HyPhy can be found.", baseDirectory.sData);
			gtk_dialog_run (GTK_DIALOG (noRez));
			gtk_widget_destroy (noRez);
			return 1;
		}
		_String rcPath = baseDir & "/theme/theme.rc";
		//printf ("Loading res files from %s\n", rcPath.sData);
		gtk_rc_parse (rcPath.sData);
	}
	
 	GlobalStartup();

	#ifdef	__HYPHYMPI__
	if (rank == 0)
	{
	#endif
	GdkDisplay * defDisplay = gdk_screen_get_display (gdk_screen_get_default());
	hSizeCursor = gdk_cursor_new_for_display (defDisplay,GDK_SB_H_DOUBLE_ARROW);
	pickUpCursor = gdk_cursor_new_for_display (defDisplay,GDK_TARGET);
	dropOffCursor = gdk_cursor_new_for_display (defDisplay,GDK_TCROSS);
	
	screenPContext = gdk_pango_context_get_for_screen (gdk_screen_get_default());
	tablePDMenuIcon = (GdkPixbuf*)ProcureIconResource(4020);
		
	/*{
		GdkScreen * defD = gdk_screen_get_default();
		fontConversionFactor = 72.27 / (gdk_screen_get_height (defD) *25.4 / gdk_screen_get_height_mm(defD)); 
		printf ("Pango conversion factor computed at: %g\n", fontConversionFactor);
	}*/

		
	ReadInTemplateFiles ();
		
	hyphyConsoleWindow = new _HYConsoleWindow ("HYPHY Console");
	SetStatusLine ("None","Idle","00:00:00");
	hyphyConsoleWindow->BringToFront();
	while (gtk_events_pending())
		gtk_main_iteration();

	ReadPreferences		();
	ReadGeneticCodes	();	
	ReadModelTemplates	();
	ReadTreeProcessors ();
	StringToConsole (hyphyCiteString);

	#ifdef __HYPHYMPI__ 
	{
		char statBuffer[1024];
		sprintf (statBuffer,"MPI version of HyPhy running on %d nodes (a master and %d compute nodes) in %s mode\n",
							 size, 
							 size-1,
							 mpiPartitionOptimizer?"partition":(mpiParallelOptimizer?"rate heterogeneity":"normal"));
		BufferToConsole (statBuffer);
	}
	#endif
		
	g_timeout_add  (100,GlobalQueueTimer,nil);
	g_timeout_add  (1000,progressTimerFunction,nil);
	gtk_main ();

	WritePreferences();
	#ifdef	__HYPHYMPI__
	}
	else // slave node
	{
		if (mpiParallelOptimizer || mpiPartitionOptimizer)
			mpiOptimizerLoop (rank, size);
		else
			mpiNormalLoop (rank, size, baseDir);
	}
	#endif

	GlobalShutdown();
    return 0;
}

