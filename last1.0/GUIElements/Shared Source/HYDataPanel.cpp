/*
	Data Panel Window
	
	Sergei L. Kosakovsky Pond, Basic version August-December 2000.
	
	Revisions (many not listed)
				  August 2001 - Removed horizontal scroll bar.
				  March  2002 - Replaced a list with a table
				  April	 2006 - Added user plug-ins	
*/


#include "HYDataPanel.h"
#include "HYSharedMain.h"
#include "HYUtils.h"
#include "batchlan.h"
#include "HYDialogs.h"
#include "HYChartWindow.h"
#include "math.h"
#include "HYParameterTable.h"
#include "HYTreePanel.h"
#include "HYButtonBar.h"
#include "HYTableComponent.h"
#include "HYEventTypes.h"

#ifdef 	  __HYPHYDMALLOC__
	#include "dmalloc.h"
#endif

extern	_HYColor navFill,
		 		 navSelect,
		 		 black;
		 
extern	_List	 dataSetList,
				 dataSetNamesList,
				 dataSetFilterNamesList,
				 dataSetFilterList,
				 modelNames,
				 theModelList,
				 likeFuncNamesList,
				 likeFuncList;
				 
extern	_SimpleList	
				 windowObjects,
				 modelMatrixIndices,
				 modelFrequenciesIndices;
				 
				 
_HYColor nucDefaults[5] = {{187,45,63},{40,148,109},{0x80,0x37,0x90},{210,135,0},{180,120,0}},
						   
		 aaDefaults[20] = {{102,51,51},{153,204,153},{102,0,153},{153,0,153},
		 				   {255,51,51},{102,0,153},{255,204,0},{153,255,153},
		 				   {102,153,51},{102,204,102},{0,51,153},{255,102,0},
		 				   {153,153,51},{255,153,0},{153,51,102},{153,153,153},
		 				   {102,102,153},{153,51,51},{204,51,51},{102,153,102}},
		 dataColors[HY_DATAPANEL_DEF_COLORS]={
		 		{220,100,3}, // Red
		 		{0,0,212},   // Blue
		 		{0,100,18},  // Dark Green
		 		{87,43,5},   // Brown
		 		{255,234,23},// Lemon
		 		{107,23,176},// Violet
		 		{255,100,2}, // Orange
		 		{2,170,234}, // Light Blue
		 		{41,25,16},  // Sepia
		 		{242,8,133}},// Pink
		 thermFill   = {245,245,245},
		 labelFColor = {254,242,208},
		 labelBColor = {102,90,88},
		 _hyWhiteColor =  {255,255,255};

_String	 nucDataType 		 ("Nucleotide"),
		 dinucDataType 		 ("Di-Nucl."),
		 codonDataType 		 ("Codon"),
		 proteinDataType	 ("Protein"),
		 binaryDataType		 ("Binary"),
		 disequonDataType	 ("Disequon"),
		 unknownDataType	 ("Unknown"),
		 none 				 ("None"),
		 makeNewTree 		 ("Create New Tree..."),
		 readTreeFile		 ("Read Tree From File..."),
		 inferTreeStr		 ("Infer Tree"),
		 useTreeModels 		 ("Use Tree Models"),
		 overlapWarning		 ("Neither overlaps more than 4 deep, nor any overlaps with more than 255 partitions are allowed."),
		 contextString1 	 ("Nucelotide Frequencies"),
		 contextString1_5 	 ("Binary Frequencies"),
		 contextString2 	 ("Aminoacid Frequencies"),
		 contextString3 	 ("Codon Frequencies"),
		 contextString4 	 ("Select all sites like this one"),
		 contextString5 	 ("Show differences from consensus"),
		 contextString6 	 ("Show proportions of rate classes"),
		 contextString7 	 ("Select all sites in class rate "),
		 contextString8		 ("Show differences from reference"),
		 contextString9		 ("Clean up singletons"),
		 contextString10 	 ("Copy consensus to clipboard"),
		 contextString11 	 ("Spawn a new data panel with selected sequences"),
		 contextString12 	 ("Filter selected sequences by ambiguity content"),
		 contextString13 	 ("Copy selected sites to clipboard as partition string"),
		 contextString14 	 ("Sort on character in this column"),
		 contextString15     ("Change highlight color"),
		 contextString16     ("Select chosen sequences in a tree"),
		 contextString17     ("Copy selected Amino Acid sites to clipboard as CODON partition string"),
		 lfKillWarning  	 ("This operation will invalidate the likelihood function attached to current dataset. You will need to rebuild the likelihood function for further analyses."),
		 lfCantKillWarning   ("Can't proceed with the operation, because the likelihood function is being used in another task."),
		 seqOmitWarning 	 ("You are about to remove all selected species from analyses. To restore them later, choose an appropriate item in the \"Data\" menu."),
		 cantOmitWarning	 ("There must be at least two sequences remaining for meaningful analyses."),
		 statusConsensus	 ("Consensus"),
		 statusGamma		 ("Rate Class"),
		 statusTranslate	 ("Translated"),
		 consensusInfoString ("Consensus"),
		 rateClassInfoString ("Rate Class"),
		 translatedInfoString("Translated"),
		 enterPartitionString("Enter a HYPHY partition specification string"),
		 useOneBasedIndices	 ("Indexing begins at 1 (default = 0)"),
		 invalidPartString 	 ("I didn't understand the partition specification string."),
		 saveDSPrompt 		   ("This data set only exists in memory. I must must write it to disk before any partitioning info can be saved."),
		 modelName 			   ("Model_Name"),
		 modelOptions 		   ("Model_Options"),
		 modelDimension 	   ("Model_Dimension"),
		 modelMatrixDimension  ("ModelMatrixDimension"),
		 modelFunction  	   ("GUIPopulateModelMatrix"),
		 efvFunction		   ("EFVEstimated"),
		 buildCodonFrequencies ("GUIBuildCodonFrequencies"),
		 selectionStatPrefix   (" Current Selection:"),
		 multiplyByFrequencies ("MULTIPLY_BY_FREQS"),
		 savedLFMatrix 		   ("LF_CACHE_MATRIX"),
		 globalPrefix		   ("globalVariable"),
		 categoryPrefix		   ("categoryVariable"),
		 modelEFVVector		   ("Model_EFV_Vector"),
		 dataPartitionIDString ("Data_Partition_ID"),
		 dataSetIDString	   ("Data_Set_ID"),
		 userModelEFV		   ("GUIHandleFrequenciesCollection"),
		 copySeqsToClip		   ("Copy Selected Sequences to Clipboard"),
		 parameterOption[3] =  {"Local",
							    "Global",
							    "Rate Het."},
							    
		 freqOption[5]		=  {"Partition",
		 					    "Dataset",
							    "Equal",
		 					    "Estimate",
		 					    "Model Spec."},
		 					    
		 nullSuffix			(" [null]"),
		 alterSuffix		(" [alternative]"),
		 inferenceDSID		("INFERENCE_DATA_SET"),
		 inferenceDWID		("INFERENCE_DATA_WINDOW"),
		 inferenceNofSeqs	("_NUMBER_OF_SEQUENCES"),
		 inferenceDummyTree ("_Internal_Inferred_Tree"),
		 inferenceDummyLF	("_INTERNAL_HYPHY_DUMMY_LF_"),
		 datapanelProcDF	("_DATAPANEL_DATAFILTER_"),
		 datapanelProcVS	("_DATAPANEL_SELECTED_SEQUENCES_"),
		 datapanelProcHS	("_DATAPANEL_SELECTED_SITES_"),
		 datapanelProcGC	("_DATAPANEL_GENETIC_CODE_"),
		 datapanelProcUS	("_DATAPANEL_UNIT_SIZE_"),
		 datapanelProcEX	("_DATAPANEL_EXCLUSIONS_"),
		 datapanelProcSF	("_DATAPANEL_SELECTED_FILTERS_"),
		 datapanelProcRF	("_DATAPANEL_RETURNED_FILTERS_"),
		 datapanelProcDS	("_DATAPANEL_DATASET_NAME_"),
		 dataPanelSearchTerm,
		 dpsandString1,
		 dpsandString2;
		 
extern	 _String dataFilePrintFormat, 
				 aminoAcidOneCharCodes, 
				 donotWarnAgain, 
				 baseDirectory, 
				 dataPanelSourcePath, 
				 dialogPrompt,
				 likefuncOutput,
				 modelGenCode,
				 blMolClock,
				 blReplicate,
				 useNexusFileData,
				 noInternalLabels,
				 VerbosityLevelString;

				 
extern	_TranslationTable	   defaultTranslationTable;

bool	 warnOmit 			= false,
		 warnKillLF 		= false,
		 autoPopLFTable     = true,
		 regExpSearch	    = false,
		 useOneBased 		= false;

bool	 HandlePreferences (_List&, _String, bool);
_String	 WriteFileDialogInput (void);
		 
_List	 dfFilterFormats,
		 geneticCodes,
		 modelTemplates,
		 dataPanelProcessors;

extern   long  likeFuncEvalCallCount;

#define	 LN_2				0.693147180559945309417232
#define	 DF_COLOR_COLUMN 	0
#define  DF_ID_COLUMN	 	1
#define  DF_TYPE_COLUMN  	2
#define  DF_TREE_COLUMN  	3
#define  DF_MODEL_COLUMN 	4
#define  DF_MDL_OPT_COLUMN	5
#define  DF_MDL_EFV_COLUMN  6
#define  DF_MDL_CLS_COLUMN  7

#define  DF_COLUMN_COUNT    8

long	 findPanelSelection = 0;
//__________________________________________________________

void	allocate_fexact_keys (long,long);
void	free_fexact_keys 	 (void);

//__________________________________________________________
_HYDataPanel::_HYDataPanel(_String& title,_String& argument):_HYTWindow (_String ("DataSet ")&title)
{	
	dataWrapper 		= nil;
	tainted 			= false;
	addedLines 			= 0;
	
	referenceSequence 	= translatedSequence 
						= 0;
						
	genCodeID 			= -1;
	lfID 				= -1;
	cantDeleteLF 		= false;
	
	flags 			   |= HY_WINDOW_STATUS_BAR_LIGHT_LEFT;
	
	_HYRect			canvasSettings = {50,50,50,200,HY_COMPONENT_V_SCROLL|HY_COMPONENT_BORDER_T|HY_COMPONENT_BORDER_B};
	
	_HYSequencePane*sp   = new _HYSequencePane (canvasSettings,GetOSWindowData(),50,200);
	canvasSettings.top   = canvasSettings.bottom = 30;
	canvasSettings.width = HY_COMPONENT_NO_SCROLL;
	_HYStretchCanvas*sc  = new _HYStretchCanvas (canvasSettings,GetOSWindowData(),50,200,32,HY_SCANVAS_HORIZONTAL);
	
	canvasSettings.top   = canvasSettings.bottom = 100;
	canvasSettings.left  = 700;
	canvasSettings.right = 10000;
	canvasSettings.width = HY_COMPONENT_V_SCROLL;
	
	_HYTable* partList    = new	_HYTable (canvasSettings, GetOSWindowData(), 1, DF_COLUMN_COUNT, 70, 100, HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT);
	checkPointer (partList);
	
	partList->selectionType = HY_TABLE_SEL_ROWS|HY_TABLE_DONT_SIZE|HY_TABLE_NODRAG_SELECTION|HY_TABLE_DONT_GROW_VERT;
	canvasSettings.top   = canvasSettings.bottom = 20;
	canvasSettings.width = HY_COMPONENT_NO_SCROLL|HY_COMPONENT_BORDER_T;
	
	_String toolTipText;

	_HYTable* partHead    = new	_HYTable (canvasSettings, GetOSWindowData(), 1, DF_COLUMN_COUNT, 80, 20, HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_BOLD);
	checkPointer (partHead);
	partHead->selectionType = HY_TABLE_SEL_ROWS|HY_TABLE_NODRAG_SELECTION;
	_SimpleList	clrIcon;
	
	clrIcon << (long)ProcureIconResource (HY_DATAPANEL_ICON_ID+12);
	clrIcon << 15;
	clrIcon << 15;

	partHead->SetCellData (&clrIcon,0,DF_COLOR_COLUMN,HY_TABLE_ICON|HY_TABLE_BEVELED|HY_TABLE_CANTSELECT,true);
	toolTipText = "Partition Name";
	partHead->SetCellData (&toolTipText,0,DF_ID_COLUMN,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_BOLD|HY_TABLE_CANTSELECT,true);
	toolTipText = "Partition Type";
	partHead->SetCellData (&toolTipText,0,DF_TYPE_COLUMN,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_BOLD|HY_TABLE_CANTSELECT,true);
	toolTipText = "Tree Topology";
	partHead->SetCellData (&toolTipText,0,DF_TREE_COLUMN,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_BOLD|HY_TABLE_CANTSELECT|HY_TABLE_PULLDOWN,true);
	toolTipText = "Substitution Model";
	partHead->SetCellData (&toolTipText,0,DF_MODEL_COLUMN,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_BOLD|HY_TABLE_CANTSELECT|HY_TABLE_PULLDOWN,true);
	toolTipText = "Parameters";
	partHead->SetCellData (&toolTipText,0,DF_MDL_OPT_COLUMN,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_BOLD|HY_TABLE_CANTSELECT|HY_TABLE_PULLDOWN,true);
	toolTipText = "Equilibrium Freqs.";
	partHead->SetCellData (&toolTipText,0,DF_MDL_EFV_COLUMN,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_BOLD|HY_TABLE_CANTSELECT|HY_TABLE_PULLDOWN,true);
	toolTipText = "Rate Classes";
	partHead->SetCellData (&toolTipText,0,DF_MDL_CLS_COLUMN,HY_TABLE_STATIC_TEXT|HY_TABLE_BEVELED|HY_TABLE_BOLD|HY_TABLE_CANTSELECT,true);


	_HYFont			tableFont;
	
	#ifdef __MAC__
		tableFont.face = "Times";
	#else
		tableFont.face = "MS Sans Serif";
	#endif
	tableFont.size = 12;
	tableFont.style = HY_FONT_PLAIN;
	
	partHead->SetFont (tableFont);
	partHead->AutoFitWidth ();
	partList->SetFont (tableFont);
	partList->AutoFitWidth (*partHead);
	partHead->SetColumnSpacing (DF_MDL_CLS_COLUMN,5000,false);
	partList->SetColumnSpacing (DF_MDL_CLS_COLUMN,5000,false);

	canvasSettings.top   = canvasSettings.bottom = 119;
	canvasSettings.left  = canvasSettings.right = 50;
	canvasSettings.width = HY_COMPONENT_BORDER_R|HY_COMPONENT_BORDER_T;
	
	_HYButtonBar*	b1 = new _HYButtonBar(canvasSettings,GetOSWindowData());
	
	canvasSettings.top	 = canvasSettings.bottom=7;
	canvasSettings.left  = 50;
	canvasSettings.right = 200;
	canvasSettings.width = HY_COMPONENT_NO_SCROLL;
	_HYStretchCanvas*marks = new _HYStretchCanvas (canvasSettings,GetOSWindowData(),5,200,32,HY_SCANVAS_HORIZONTAL);
	
	canvasSettings.left = canvasSettings.top = canvasSettings.bottom = 14;
	canvasSettings.width = HY_COMPONENT_BORDER_B;
	
	
	_HYSequencePane*sp2 = new _HYSequencePane (canvasSettings,GetOSWindowData(),14,200);
	sp2->SetActiveOrPassive (false);
	sp2->SetNumberDisplay   (false);
	
	sc->SetMessageRecipient (this);
	sp->SetMessageRecipient (this);
	b1->SetMessageRecipient (this);
	sp2->SetMessageRecipient(this);
	partList->SetMessageRecipient(this);
	partHead->SetMessageRecipient(this);
	toolTipText = ("Cut partition in 2");
	b1->AddButton (ProcureIconResource(HY_DATAPANEL_ICON_ID),&toolTipText);
	toolTipText = "Join 2 Partitions";
	b1->AddButton (ProcureIconResource(HY_DATAPANEL_ICON_ID+1),&toolTipText);
	toolTipText = "Subtract 2 Overlapping Partitions";
	b1->AddButton (ProcureIconResource(HY_DATAPANEL_ICON_ID+2),&toolTipText);
	toolTipText = "Delete Partition";
	b1->AddButton (ProcureIconResource(HY_DATAPANEL_ICON_ID+3),&toolTipText);
	toolTipText = "Comb Selection/Partition";
	b1->AddButton (ProcureIconResource(HY_DATAPANEL_ICON_ID+4),&toolTipText);
	toolTipText = "Save Partition to Disk";
	b1->AddButton (ProcureIconResource(HY_DATAPANEL_ICON_ID+5),&toolTipText);
	toolTipText = "Interleave Disjoint Nucleotide Partitions";
	b1->AddButton (ProcureIconResource(HY_DATAPANEL_ICON_ID+6),&toolTipText);
	toolTipText = "Data Operations";
	b1->AddButton (ProcureIconResource(HY_DATAPANEL_ICON_ID+11),&toolTipText);
	toolTipText = "Display Table of Parameter Values";
	b1->AddButton (ProcureIconResource(HY_DATAPANEL_ICON_ID+10),&toolTipText);
	toolTipText = "Toggle Coloring Mode";
	b1->AddButton (ProcureIconResource(HY_DATAPANEL_ICON_ID+13),&toolTipText);
	
	b1->EnableButton(0,false);
	b1->EnableButton(1,false);
	b1->EnableButton(2,false);
	b1->EnableButton(3,false);
	b1->EnableButton(4,false);
	b1->EnableButton(5,false);
	b1->EnableButton(6,false);
	b1->EnableButton(7,false);
	b1->EnableButton(8,false);
	b1->MarkAsPullDown (2,true);
	b1->MarkAsPullDown (6,true);
	b1->MarkAsPullDown (7,true);
	b1->SetButtonDim(16);
	b1->SetButtonLayoutW (2);
	b1->SetBackColor(labelBColor);
	

	marks->SetMessageRecipient (this);
	AddObject (sp);      // 0
	AddObject (sc);      // 1
	AddObject (b1);      // 2
	AddObject (marks);   // 3
	AddObject (sp2);     // 4
	AddObject (partList);// 5
	AddObject (partHead);// 6

	SetTableDimensions (5,2);
	SetCell (0,0,sc);
	SetCell (0,1,sc);
	SetCell (1,0,marks);
	SetCell (1,1,marks);
	SetCell (2,0,sp);
	SetCell (2,1,sp);
	SetCell (3,0,b1);
	SetCell (3,1,partHead);
	SetCell (4,0,b1);
	SetCell (4,1,partList);


	DeleteObject (sc);
	DeleteObject (sp);
	DeleteObject (b1);
	DeleteObject (marks);
	DeleteObject (sp2);
	DeleteObject (partList);
	DeleteObject (partHead);
	
	if (dataPanelProcessors.lLength==0)
		ReadDataPanelProcessors ();

	SetDataSetReference(argument);
}

//__________________________________________________________
_HYDataPanel::~_HYDataPanel()
{
	if (lfID>=0)
		PurgeLF(false);
		
	if (dataWrapper)
		DeleteObject (dataWrapper);
}

//__________________________________________________________
void	_HYDataPanel::GenerateStatusLine (void)
{
	_String	 statBar ("Custom Character Data. ");
	if (dataType&HY_DATAPANEL_NUCDATA)
		statBar = "Nucleotide Data. ";
	else
		if (dataType&HY_DATAPANEL_PROTDATA)
			statBar = "Aminoacid Data. ";
		else 
			if (dataType&HY_DATAPANEL_BINARYDATA)
				statBar = "Binary Data. ";
		
	if (!dataWrapper)
	{
		_DataSet *theDS = (_DataSet*)dataSetList(dataSetID);
		statBar = statBar & _String (theDS->NoOfColumns()) &" sites ("&_String (theDS->NoOfUniqueColumns())&" distinct patterns), "&
			_String(theDS->NoOfSpecies())& " species.";
	}
	else
		statBar = statBar & _String (dataWrapper->GetFullLengthSpecies()) &" sites ("&_String (dataWrapper->NumberDistinctSites())&" distinct patterns), "&
			_String(dataWrapper->NumberSpecies())& " species.";
	statBar = statBar & selectionStatPrefix & "empty";
	SetStatusBar(statBar);	
}
//__________________________________________________________

bool	_HYDataPanel::ProcessGEvent (_HYEvent* e)
{	
	_String firstArg,
			secondArg;
	long    k,f;
	
	bool	done = false;
	
	if (e->EventClass()==_hyGlobalLFKillEvent)
	{
		firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
		k = firstArg.toNum();
		if (k!=GetID())
		{
			firstArg 		= e->EventCode().Cut (f+1,-1);
			if (lfID==firstArg.toNum())
			{
				lfID = -1;
				_PaintLFStatus();
			}
		}	
		done = true;
	}
	else
		if (e->EventClass()==_hyGlobalDFKillEvent)
		{
			if (lfID>=0)
			{
				firstArg = e->EventCode().Cut ((f=e->EventCode().Find(','))+1,-1);
				k = firstArg.toNum();
				if (((_LikelihoodFunction*)likeFuncList(lfID))->DependOnDF(k))
				{
					postLFKillEvent (-1,lfID);
					KillLFRecord(lfID);
					lfID = -1;
				}
			}
			done = true;
		}
		else
			if (e->EventClass()==_hyGlobalTreeKillEvent)
			{
				firstArg = e->EventCode().Cut ((f=e->EventCode().Find(','))+1,-1);
				k = firstArg.toNum();
				if ((f = treeVarReferences.Find(k))>=0)
				{
					treeVarReferences.lData[f] = -1;
					_List mOptions;
					GenerateModelList (mOptions,f);
					RefreshPartRow(mOptions,f,true);
				}
				done = true;
			}
			else
				if (e->EventClass()==_hyGlobalDSKillEvent)
				{
					firstArg = e->EventCode().Cut ((f=e->EventCode().Find(','))+1,-1);
					k = firstArg.toNum();
					if (dataSetID==k)
					{
						tainted = false;
						postWindowCloseEvent (GetID());
					}
					done = true;
				}
	if (done)
		return true;
	return _HYWindow::ProcessGEvent (e);
}
//__________________________________________________________

void	_HYDataPanel::RefreshPartRow (_List& modelList, long r, bool doTree, bool doNonTree)
{
	_HYTable*	pl = (_HYTable*)GetObject (5);
	_SimpleList	changedCells;
	long		baseIndex = r*pl->horizontalSpaces.lLength,
				selected = (pl->cellTypes[r*pl->horizontalSpaces.lLength] & HY_TABLE_SELECTED)?HY_TABLE_SELECTED:0,
				color	 = ((_SimpleList*)pl->GetCellData (DF_COLOR_COLUMN,r))->lData[0];

	
	if (doNonTree)
	{
		int			model,
					params,
					freqs,
					rates;
					
		LongToModelData (modelReferences.lData[r],model,params,freqs,rates);
		
		_String *   rightModel  = (_String*)modelList(model),
			    * 	rightParams = &empty,
			    *   rightFreqs  = &empty,
			    	rightRates,
			    
			    *	presModel   = (_String*) pl->GetCellData (DF_MODEL_COLUMN,r),
			    *	presParams  = (_String*) pl->GetCellData (DF_MDL_OPT_COLUMN,r),
			    *	presFreqs   = (_String*) pl->GetCellData (DF_MDL_EFV_COLUMN,r),
			    *	presRates   = (_String*) pl->GetCellData (DF_MDL_CLS_COLUMN,r);
			    
		if (model)
		{
			rightParams = &parameterOption[params];
			rightFreqs  = &freqOption[freqs];
			
			if (rates)
				rightRates = _String((long)rates);
		}
			    	
		if (!presModel->Equal(rightModel))
		{
			changedCells << baseIndex+DF_MODEL_COLUMN;
			pl->SetCellData (rightModel,r,DF_MODEL_COLUMN,selected|HY_TABLE_STATIC_TEXT|HY_TABLE_PULLDOWN,true);
		}
		if (!presParams->Equal(rightParams))
		{
			changedCells << baseIndex+DF_MDL_OPT_COLUMN;
			if (rightParams->sLength)
				pl->SetCellData (rightParams,r,DF_MDL_OPT_COLUMN,selected|HY_TABLE_STATIC_TEXT|HY_TABLE_PULLDOWN,true);
			else
				pl->SetCellData (rightParams,r,DF_MDL_OPT_COLUMN,selected|HY_TABLE_STATIC_TEXT,true);
		}
		if (!presFreqs->Equal(rightFreqs))
		{
			changedCells << baseIndex+DF_MDL_EFV_COLUMN;
			if (rightFreqs->sLength)
				pl->SetCellData (rightFreqs,r,DF_MDL_EFV_COLUMN,selected|HY_TABLE_STATIC_TEXT|HY_TABLE_PULLDOWN,true);
			else
				pl->SetCellData (rightFreqs,r,DF_MDL_EFV_COLUMN,selected|HY_TABLE_STATIC_TEXT,true);
				
		}
		if (!presRates->Equal(&rightRates))
		{
			changedCells << baseIndex+DF_MDL_CLS_COLUMN;
			if (rightRates.sLength)
				pl->SetCellData (&rightRates,r,DF_MDL_CLS_COLUMN,selected|HY_TABLE_EDIT_TEXT,true);		
			else
				pl->SetCellData (&rightRates,r,DF_MDL_CLS_COLUMN,selected|HY_TABLE_STATIC_TEXT,true);
		}
		if (color != partitionColors.lData[r])
		{
			((_SimpleList*)pl->GetCellData (DF_COLOR_COLUMN,r))->lData[0] = partitionColors.lData[r];
			changedCells << baseIndex+DF_COLOR_COLUMN;
		}
	}
	
	if (doTree)
	{
		_String *	rightTree	= &none,
			    *	presTree    = (_String*) pl->GetCellData (DF_TREE_COLUMN,r);
			    
		if (treeVarReferences.lData[r]>=0)
			rightTree = LocateVar(treeVarReferences.lData[r])->GetName();
			
		if (!presTree->Equal(rightTree))
		{
			changedCells << baseIndex+DF_TREE_COLUMN;
			pl->SetCellData (rightTree,r,DF_TREE_COLUMN,selected|HY_TABLE_STATIC_TEXT|HY_TABLE_PULLDOWN,true);
		}
	}
	
	if (changedCells.lLength)
		pl->_MarkRowForUpdate (r);
}

//__________________________________________________________

void  _HYDataPanel::HandleSearchAndReplace (void)
{
	_String p1 ("Search for reg.exp.:"),
			p2 ("Replace with string:");
			
	if (EnterString2Dialog (dpsandString1, dpsandString2, p1, p2, (Ptr)this))
	{
		int errNo = 0;
		Ptr regex = PrepRegExp (&dpsandString1, errNo, true);

		if (!regex)
			WarnError (GetRegExpError (errNo));

		_HYSequencePane* sp = (_HYSequencePane*)components (0);
		
		_List	oldNames (sp->rowHeaders),
				newNames;
				
		bool 	doSomething = false;
				
		for (long k=0; k<oldNames.lLength; k++)
		{
			_SimpleList matches;
			_String * theString = (_String*)oldNames(k);
			theString->RegExpMatchAll(regex, matches);
			if (matches.lLength)
			{
				_String * newString = new _String (theString->sLength+1,true);
				checkPointer (newString);
				
				long	  idx  = matches.lData[0],
						  midx = 0;
				
				for (long k2=0; k2<theString->sLength;)
				{
					if (k2==idx)
					{
						(*newString) << dpsandString2;
						k2 = matches.lData[midx+1]+1;
						midx += 2;
						if (midx == matches.lLength)
							idx = -1;
						else
							idx = matches.lData[midx];
					}
					else
						(*newString) << theString->sData[k2++];
				}
				newString->Finalize();
				newNames << newString;
				DeleteObject (newString);
				doSomething = true;
			}
			else
				newNames && theString;
		}
		
		FlushRegExp (regex);
		if (doSomething)
			sp->BatchRenameSequences (oldNames, newNames);
	}
}
				
//__________________________________________________________

bool	_HYDataPanel::ProcessEvent (_HYEvent* e)
{
	long 		f,i,k,p;
	_String 	firstArg;
	bool		done = false;
	_HYTable*	pl = (_HYTable*)GetObject (5);
	
	if (e->EventClass()==_hyMenuSelChangeEvent)
	{
		firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
		k = firstArg.toNum();
		for (i=0;i<components.lLength;i++)
		{
			if (((_HYGuiObject*)components(i))->MatchID(k))
				break;
		}
		firstArg = e->EventCode().Cut (f+1,-1);
		k = firstArg.toNum();
		if (i==0)
	    // selection change in the sequence pane
		{
			UpdateSelDepPartitionOperations();
			_UpdateSelectionChoices(IsSelectionNonEmpty());
		}		
		tainted = true;
		done = true;
	}
	else
	{
		if (e->EventClass()==_hyRebuildSCanvasEvent)
		{
			k = e->EventCode().toNum();
			for (i=0;i<components.lLength;i++)
			{
				if (((_HYGuiObject*)components(i))->MatchID(k))
					break;
			}
			if (i==1)
				BuildThermometer();
			else
				BuildMarksPane();
			done = true;
		}
		else
		if (e->EventClass()==_hyScrollingEvent)
		{
			_HYSequencePane* sp = (_HYSequencePane*)components (0),
						   * sp2= (_HYSequencePane*)components (4);
				
			if (addedLines)			   
				if (sp2->startColumn!=sp->startColumn)
					sp2->HScrollPane (sp->startColumn-sp2->startColumn);
			BuildMarksPane();
			_PaintThermRect(true);
			done = true;
		}
		else
		if (e->EventClass()==_hyTableChangeSelEvent)
		{
			UpdatePartitionOperations();
			done = true;
		}
		else
		{
			if (e->EventClass()==_hyButtonPushEvent)
			{
				tainted = true;
				firstArg = e->EventCode().Cut (0,(f=e->EventCode().Find(','))-1);
				_HYSequencePane* sp = (_HYSequencePane*)GetObject(0);
				k = firstArg.toNum();
				for (i=0;i<components.lLength;i++)
				{
					if (((_HYGuiObject*)components(i))->MatchID(k))
						break;
				}
				firstArg = e->EventCode().Cut (f+1,-1);
				k = firstArg.toNum();
				if (i==2) // button bar
				{
					_HYButtonBar*   bb  = (_HYButtonBar*)GetObject (2);
					_SimpleList		sel;
					
					pl->GetRowSelection (sel);
					
					int h,v;
					
					bb->GetButtonLoc(k,h,v,true);
					
					switch (k)
					{
						case 0: // split button
						{
							SplitPartition(sel.lData[0],sp->selection.lData[0]);
							break;
						}
						case 1: // join button
						{
							if (sel.lLength>=2)
							{
								h = sel.lData[0];
								for (long pi=1; pi<sel.lLength; pi++)
									JoinPartitions(h,sel.lData[pi]-(pi-1));
							}
							else	
								JoinSpeciesDisplay ();
							break;
						}
						case 2: // 
						{
							if (sel.lLength>=2)
							{
								int i=sel.lData[0],
									j=sel.lData[1];
									
								_String	* p1 = (_String*)pl->GetCellData(DF_ID_COLUMN,i),
									    * p2 = (_String*)pl->GetCellData(DF_ID_COLUMN,j);
								
								firstArg = *p1 & '=' & *p1 & '-' &*p2;
								_List  menuOptions;
								menuOptions && & firstArg;
								firstArg = *p2 & '=' & *p2 & '-' &*p1;
								menuOptions && & firstArg;
								firstArg = HandlePullDown (menuOptions,h,v,0);
								bb->_UnpushButton();
								if (firstArg.sLength)
								{
									if (firstArg.beginswith(*p1&'='))
										SubtractPartitions (i,j);
									else
										SubtractPartitions (j,i);
								}
							}
							break;
						}
						case 3: // kill button
						{
							if (sel.lLength)
								for (long k=sel.lLength-1; k>=0; k--)
									KillPartition(sel.lData[k]);
							else
								OmitSelectedSpecies();
							break;
						}

						case 4: // comb button
						{
							CombPartition((sel.lLength)?sel.lData[0]:-1);
							break;
						}
						
						case 5: // save file
						{
							SavePartition(sel.lData[0]);
							break;
						}
						case 6: // interleave
						{
							if (sel.lLength==2)
							{
								int i=sel.lData[0],
									j=sel.lData[1];
								
								_String	* p1 = (_String*)pl->GetCellData(DF_ID_COLUMN,i),
									    * p2 = (_String*)pl->GetCellData(DF_ID_COLUMN,j);

								firstArg = *p1 & '&' & *p2;
								_List  menuOptions;
								menuOptions && & firstArg;
								firstArg = *p1 & '&' & *p2 & "[Reverse]";
								menuOptions && & firstArg;
								firstArg = *p1 & "[Reverse]&" & *p2;
								menuOptions && & firstArg;
								firstArg = *p1 & "[Reverse]&" & *p2 & "[Reverse]";
								menuOptions && & firstArg;
								bb->_UnpushButton();
								firstArg = HandlePullDown (menuOptions,h,v,0);
								if (firstArg.sLength)
								{
									bool d1 = firstArg.beginswith(*p1&'&'),	
										 d2 = firstArg.endswith ("[Reverse]");
									InterleavePartitions (i,j,!d1,d2);
								}
							}
							break;
						}
						case 7: // deletions/constants
						{
							if (sel.lLength==1)
							{
								_DataSetFilter* df = (_DataSetFilter*)dataSetFilterList(dataPartitions.lData[sel.lData[0]]);
								firstArg = "Constant Sites"; // 0
								_List  menuOptions;
								menuOptions && &firstArg;
								firstArg = "Constant Sites Matching Ambiguities"; // 1
								menuOptions && &firstArg;
								firstArg = menuSeparator;
								menuOptions && &firstArg;
								firstArg = "Sites with Deletions"; // 3
								menuOptions && &firstArg;
								firstArg = "Sites with All Deletions"; // 4
								menuOptions && &firstArg;
								firstArg = "Sequences with Deletions"; //5
								menuOptions && &firstArg;
								firstArg = menuSeparator;
								menuOptions && &firstArg;
								firstArg = "Identical Sequences"; // 7
								menuOptions && &firstArg;
								firstArg = "Identical Sequences Matching Ambiguities"; //8
								menuOptions && &firstArg;
								menuOptions && &menuSeparator;
								firstArg = "Association [Pearson chi2]"; // 10
								menuOptions && &firstArg;
								firstArg = "Association [Fisher Exact]"; // 11
								menuOptions && &firstArg;
								firstArg = "Show Contigency Tables"; // 12
								menuOptions && &firstArg;
								menuOptions && &menuSeparator;
								firstArg = "Character Usage By Site"; // 14
								menuOptions && &firstArg;
								firstArg = "Character Entropy"; // 15
								menuOptions && &firstArg;
								firstArg = "Ambiguous Characters By Sequence"; // 16
								menuOptions && &firstArg;
								if (df->GetUnitLength() == 3) 
								{
									firstArg = "Codon Usage"; // 17
									menuOptions && &firstArg;
								}
								else
									if (df->GetUnitLength() == 1) 
									{
										firstArg = "Count differences from consensus"; // 17
										menuOptions && &firstArg;
										firstArg = "Clean up singletons"; //18
										menuOptions && &firstArg;
									}

								bb->_UnpushButton();
								firstArg = HandlePullDown (menuOptions,h,v,0);
								k = menuOptions.Find (&firstArg);
								switch (k)
								{
									case 0:
										ShowConstantSites (false);
										break;
									case 1:
										ShowConstantSites (false,true);
										break;
									case 3:
										ShowConstantSites (true);
										break;
									case 4:
										ShowConstantSites (true,true);
										break;
									case 5:
										ShowConstantSites (true,true,true);
										break;
									case 7:
										ShowDuplicateSequences (false);
										break;
									case 8:
										ShowDuplicateSequences (true);
										break;
									case 10:
									case 11:
									case 12:
										ShowAssociation (sel.lData[0],k-10);
										break;
									case 14:
										ShowCharacterUsage (sel.lData[0],false);
										break;
									case 15:
										ShowCharacterUsage (sel.lData[0],true);
										break;
									case 17:
										if (df->GetUnitLength() == 3) 
											ShowCodonUsage (sel.lData[0]);
										else
											ShowMutantCount (sel.lData[0]);
										break;
										
									case 18:
										CleanupSingletons (sel.lData[0]);
										break;
								}			
								bb->_UnpushButton();
							}
						}
						break;
														
						case 8: // display parameter table
							DisplayParameterTable ();
							bb->_UnpushButton();
							break;
							
						case 9: // color toggle 
							sp->invertColors = !sp->invertColors;
							bb->ReplaceButton (9,ProcureIconResource(HY_DATAPANEL_ICON_ID+13+(sp->invertColors)),nil);
							sp->BuildPane();
							sp->_MarkForUpdate();
							_HYSequencePane * sp2	   = (_HYSequencePane*)GetObject (4);	
							sp2->invertColors = !sp2->invertColors;
							sp2->BuildPane();
							sp2->_MarkForUpdate();
							break;
					}
				}
				done = true;
			}
			else
			{
				if (e->EventClass()==_hyTableDblClickEvent)
				{
					if (pl->IsRowSelectionSimple())
						EditPartitionProperties(pl->GetFirstRowSelection());
					
					done = true;
				}
				else
					if (e->EventClass()==_hyTablePullDownEvent)
					{
						f = e->EventCode().Find(',',0,-1);
						p = e->EventCode().FindBackwards(',',0,-1);
						k = e->EventCode().Cut (0,f-1).toNum();
						f = e->EventCode().Cut (f+1,p-1).toNum();
						p = e->EventCode().Cut (p+1,-1).toNum();
						
						for (i=0;i<components.lLength;i++)
						{
							if (((_HYGuiObject*)components(i))->MatchID(k))
								break;
						}
						if ((i==5)||(i==6))
						{
							k = p>>16;   // x - coord
							p = p&0xffff;// y - coord
							int 	r = f/pl->horizontalSpaces.lLength,
									c = f%pl->horizontalSpaces.lLength;
							_List	mOptions;
							
							_String * currentType = (_String*)pl->GetCellData (c,r),
									  mChoice;


							switch (c)
							{
								case DF_TYPE_COLUMN: // data type change
								{
									if (i==5)
									{
										ConstructDataTypeOptions (mOptions);
										mChoice = HandlePullDown (mOptions,k,p,mOptions.Find (currentType)+1);
										if (mChoice.sLength&&(!mChoice.Equal(currentType)))
										{
											if (DataTypeChange (r,mOptions.Find (&mChoice)))
											{
												pl->SetCellData (&mChoice,r,c,pl->cellTypes.lData[f],true);
												mOptions.Clear();
												GenerateModelList (mOptions,r);
												RefreshPartRow(mOptions,r,false);
												pl->_MarkCellForUpdate(f);
												tainted = true;
											}
										}
									}
									break;
								}
								case DF_TREE_COLUMN: // tree topology change
								{
									if ((i==6)&&((lfID>=0)||(dataPartitions.lLength==0)))
										break;
									GenerateTreeList (mOptions);
									mChoice = HandlePullDown (mOptions,k,p,mOptions.Find (currentType)+1);
									if (mChoice.sLength)
									{
										_String mRes;
										if (i==5)
										{
											if (mChoice.Equal(currentType))
												break;
												
											mRes = TreeTopologyChange (r,&mChoice);
											pl->SetCellData (&mRes,r,c,pl->cellTypes.lData[f],true);
											pl->_MarkCellForUpdate(f);
											tainted = true;
										}
										else
										{
											for (i=0; i<dataPartitions.lLength; i++)
												if (!mChoice.Equal((_String*)pl->GetCellData (c,i)))
													break;
													
											if (i==dataPartitions.lLength)
												break;
											
											mRes = TreeTopologyChange (0,&mChoice); 
											pl->SetCellData (&mRes,0,c,pl->cellTypes.lData[c],true);
											f = c;
											for (i=1; i<dataPartitions.lLength; i++)
											{
												f += pl->horizontalSpaces.lLength;
												mRes = TreeTopologyChange (i,&mRes); 
												pl->SetCellData (&mRes,i,c,pl->cellTypes.lData[f],true);
											}
											pl->_MarkColumnForUpdate(c);
											tainted = true;
										}
									}
									_VerifyInferMenu();
									break;
								}
								
								case DF_MODEL_COLUMN: // model change
								{
									if ((i==6)&&((lfID>=0)||(dataPartitions.lLength==0)))
										break;
										
									if (i==5)
										GenerateModelList (mOptions,r);
									else
									{
										_String   *firstDataType = (_String*)pl->GetCellData(DF_TYPE_COLUMN,0);
										for (f=1; f<dataPartitions.lLength; f++)
											if (!firstDataType->Equal((_String*)pl->GetCellData(DF_TYPE_COLUMN,f)))
												break;

										if (f==dataPartitions.lLength)
											GenerateModelList (mOptions,0);
										else
											break;
									}
									mChoice = HandlePullDown (mOptions,k,p,mOptions.Find (currentType)+1);
									if (mChoice.sLength)
									{
										if (i==5)
										{
											if (mChoice.Equal(currentType))
												break;
										}
										else
										{
											for (f=0; f<dataPartitions.lLength; f++)
												if (!mChoice.Equal((_String*)pl->GetCellData (c,f)))
													break;
													
											if (f==dataPartitions.lLength)
												break;
										}

										k = mOptions.Find(&mChoice);
										ModelChange (&mChoice,r,k);
										RefreshPartRow(mOptions,r,false);
										if (i==6)
										{
											for (i=1; i<dataPartitions.lLength; i++)
											{
												ModelChange  (&mChoice,i,k);
												RefreshPartRow(mOptions,i,false);
											}
										}
										tainted = true;
									}
									_VerifyInferMenu();
									break;
								}

								case DF_MDL_OPT_COLUMN: 
								case DF_MDL_EFV_COLUMN:
								{
									if ((i==6)&&((lfID>=0)||(dataPartitions.lLength==0)))
										break;
										
									if (c==DF_MDL_OPT_COLUMN)
										GenerateModelPOptionList (mOptions,r);
									else
										GenerateModelFOptionList (mOptions,r);
									
									if (i==6)
									{
										for (f=1; f<dataPartitions.lLength; f++)
											if (!currentType->Equal ((_String*)pl->GetCellData (DF_MODEL_COLUMN,f)))
											{
												_List lOptions;
												if (c==DF_MDL_OPT_COLUMN)
													GenerateModelPOptionList (lOptions,f);
												else
													GenerateModelFOptionList (lOptions,f);
													
												for (long jj = 0; jj<mOptions.countitems();)
													if (lOptions.Find (mOptions (jj))<0)
														mOptions.Delete (jj);
													else
														jj++;
											}
									}
									
									mChoice = HandlePullDown (mOptions,k,p,mOptions.Find (currentType)+1);
									
									if (mChoice.sLength&&(!mChoice.Equal(currentType)))
									{
										if (c==DF_MDL_OPT_COLUMN)
											ModelOptionChange (&mChoice,r);
										else
											ModelFreqChange (&mChoice,r);
											
										mOptions.Clear();
										GenerateModelList (mOptions,r);
										RefreshPartRow(mOptions,r,false);
										
										if (i==6)
											for (f=1; f<dataPartitions.lLength; f++)
											{
												if (c==DF_MDL_OPT_COLUMN)
													ModelOptionChange (&mChoice,f);
												else
													ModelFreqChange (&mChoice,f);
													
												mOptions.Clear();
												GenerateModelList (mOptions,f);
												RefreshPartRow(mOptions,f,false);
											}
										tainted = true;
									}
									break;
								}

								/*case 5: // freq change
								{
									mChoice = HandlePullDown (mOptions,k,p,mOptions.Find (currentType)+1);
									if (mChoice.sLength&&(!mChoice.Equal(currentType)))
									{
										ModelFreqChange (&mChoice,r);
										mOptions.Clear();
										GenerateModelList (mOptions,r);
										RefreshPartRow(mOptions,r,false);
									}
									break;
								}*/
							}
						}
						done = true;
					}
					else
						if (e->EventClass()==_hyContextPopUp)
						{
							f = e->EventCode().Find(',',0,-1);
							p = e->EventCode().FindBackwards(',',0,-1);
							if (p>f+1)
							{
								
								k = e->EventCode().Cut (0,f-1).toNum();
								f = e->EventCode().Cut (f+1,p-1).toNum();
								p = e->EventCode().Cut (p+1,-1).toNum();
									
								for (i=0;i<components.lLength;i++)
								{
									if (((_HYGuiObject*)components(i))->MatchID(k))
										break;
								}
								if (i==0)
									ProcessContextualPopUpMain(f,p);
								else
								if (i==4)
									ProcessContextualPopUpAux(f,p);
									
								
								done = true;
							}
						}
						else
							if (e->EventClass()==_hyTableEditCellEvent)
							{
								f = e->EventCode().Find(',',0,-1);
								k = e->EventCode().Cut (0,f-1).toNum();
								f = e->EventCode().Cut (f+1,-1).toNum();
							
								for (i=0;i<components.lLength;i++)
									if (((_HYGuiObject*)components(i))->MatchID(k))
										break;

								if (i==5)
								{
									int r = f/pl->horizontalSpaces.lLength,
										c = f%pl->horizontalSpaces.lLength;
									
									if (c==DF_MDL_CLS_COLUMN)
									{
										k = ((_String*)pl->GetCellData (c,r))->toNum();
										if  (k>1)
											ModelRateClassChange (r,k);
										else
										{
											_String warnMsg ;
											warnMsg = *(_String*)pl->GetCellData (c,r) & " is not an integer greater than 1.";
											ProblemReport (warnMsg, (Ptr)this);
										}
										_List mOptions;
										GenerateModelList (mOptions,r);
										RefreshPartRow(mOptions,r,false);
									}
								}
								done = true;
							}
							else
								if (e->EventClass()==_hyTableResizeCEvent )
								{
									f = e->EventCode().Find(',',0,-1);
									k = e->EventCode().Find(',',f+1,-1);
									long g = e->EventCode().Cut (0,f-1).toNum();
									f = e->EventCode().Cut (f+1,k-1).toNum();
									k = e->EventCode().Cut (k+1,-1).toNum();
								
									for (i=0;i<components.lLength;i++)
										if (((_HYGuiObject*)components(i))->MatchID(g))
											break;
											
									if (i==6)
									{
										pl->SetColumnSpacing (f,k,true);
									}
									
									done = true;
								}
			}
		}
	}
	if (done) 
	{
		DeleteObject (e);
		return true;
	}
	return _HYTWindow::ProcessEvent(e);	
}

//__________________________________________________________

void	_HYDataPanel::GetMutantCount (long partIndex, _SimpleList& counts, _String* storeConsensus)
{
	_DataSetFilter *df = (_DataSetFilter*)dataSetFilterList(dataPartitions.lData[partIndex]);

	_SimpleList	   				mutantCount;
	
	if (storeConsensus)
		*storeConsensus = df->GenerateConsensusString (&mutantCount);
	else
		df->GenerateConsensusString (&mutantCount);
	
	long			seqCount = df->NumberSpecies();
					//half	 = seqCount/2;
	
	for	(long k=0; k<df->theFrequencies.lLength; k++)
		//if (mutantCount.lData[k]>half)
			mutantCount.lData[k] = seqCount-mutantCount.lData[k];
		//else
			//mutantCount.lData[k] = -1;
			
	for (long kk=0; kk<df->duplicateMap.lLength; kk++)
		counts << mutantCount.lData[df->duplicateMap.lData[kk]];
}

//__________________________________________________________

void	_HYDataPanel::ShowMutantCount (long partIndex)
{
	_SimpleList		  counts;
	GetMutantCount (partIndex, counts, nil);
	
	_List 		 colHeaders;
	_String		 aString ("Mutant Count");
	
	_Matrix		 res (counts.lLength, 1, false, true);
	
	for (long k=0; k<counts.lLength; k++)
		res.theData[k] = counts.lData[k];

	colHeaders && & aString;
		
	aString = _String("Mutant Count For ") & 
				*(_String*)dataSetFilterNamesList(dataPartitions.lData[partIndex]);
		
	long f = FindWindowByName (aString);
	_HYChartWindow * nc;
	if (f>=0)
	{
		nc = (_HYChartWindow*)windowObjectRefs (f);
		nc->SetTable (colHeaders,res);
	}
	else
	{
		nc = new _HYChartWindow (aString, colHeaders, res, nil);
		checkPointer (nc);	
	}
	nc->SetChartType ("Bar Chart","Index","Mutant Count", true);
	nc->BringToFront();
	
}

//__________________________________________________________

void	_HYDataPanel::CleanupSingletons (long partIndex)
{
	_DataSetFilter *df = (_DataSetFilter*)dataSetFilterList(dataPartitions.lData[partIndex]);
	
	if (df->GetUnitLength() == 1)
	{
		_String 		  consensus;
		_SimpleList		  counts;
		
		GetMutantCount (partIndex, counts, &consensus);
		
		_SimpleList	   singletonCount;
		
		for	(long k=0; k<counts.lLength; k++)
			if (counts.lData[k] == 1)
				singletonCount << k;
				
		if (singletonCount.lLength == 0)
		{
			consensus = "No singleton mutations found.";
			ProblemReport (consensus, (Ptr)this);
		}
		else
		{
			_DataSet*		newDS = new _DataSet();
			_String			state(1,false);
			long			k,
							shifter = 0;
							
			singletonCount << -1;
					
			for (k=0; k<consensus.sLength; k++)
			{
				if (k == singletonCount.lData[shifter])
				{
					shifter++;
					newDS->AddSite (consensus.sData[k]);
				}
				else
				{
					df->RetrieveState (k,0,state);
					newDS->AddSite (state.sData[0]);
				}
			}
				
			newDS->AddName(*df->GetSequenceName(0));
						
			for (k=1; k<df->NumberSpecies(); k++)
			{
				newDS->AddName(*df->GetSequenceName(k));
				shifter = 0;
				for (long j=0; j<consensus.sLength; j++)
					if (j == singletonCount.lData[shifter])
					{
						shifter++;
						newDS->Write2Site (j,consensus.sData[j]);
					}
					else
					{
						df->RetrieveState (j,k,state);
						newDS->Write2Site (j,state.sData[0]);
					}
			}
				
			newDS->Finalize();
			newDS->SetTranslationTable (df->GetData()->GetTT());
			newDS->SetNoSpecies (df->NumberSpecies());
			consensus = (*(_String*)dataSetFilterNamesList(partIndex))&("_cleaned");
			AddDataSetToList (consensus, newDS);
			_HYDataPanel* myTWindow = new _HYDataPanel(consensus,consensus);
			myTWindow->BringToFront();
			consensus = _String("Cleaned up ") & (long)singletonCount.lLength & " singletons.";
			ProblemReport (consensus, (Ptr)myTWindow);
		}
	}
}

//__________________________________________________________

void	_HYDataPanel::ShowCharacterUsage (long partIndex, bool useEntropy)
{
	_DataSetFilter *df = (_DataSetFilter*)dataSetFilterList(dataPartitions.lData[partIndex]);
		
	long			stateCount 		= df->GetDimension (true),
					fullCount		= df->GetDimension (false),
					stateSize		= df->GetUnitLength(),
					siteCount  		= df->GetFullLengthSpecies()/stateSize,
					patternCount	= df->NumberDistinctSites();

	_Matrix			res  (useEntropy?1:stateCount, siteCount, false, true),
					res2 (useEntropy?1:stateCount, patternCount, false, true);
					
	_SimpleList		mapper;
	
	if (df->theExclusions.lLength)
	{
		long k = 0;
		for (long i = 0; i<fullCount; i++)
		{	
			if (i==df->theExclusions[k])
			{
				k++;
				continue;
			}
			mapper << i;
		}
	}
	else
		for (long k0=0; k0<stateCount; k0++)
			mapper << k0;
	
	for (long k=0; k<patternCount; k++)
	{						
		_SimpleList			 vSeq;
		
		for (long k2=0; k2<stateSize; k2++)
			vSeq << df->theMap[k*stateSize+k2];
		
		_Matrix* siteFreqs = df->GetData()->HarvestFrequencies (stateSize,stateSize,0,df->theNodeMap, vSeq, false);
		
		if (useEntropy)
		{
			_Parameter entropy = 0.0;
			
			for (long k3=0; k3<stateCount; k3++)
			{
				_Parameter siteFreq = siteFreqs->theData[mapper.lData[k3]];
				if (siteFreq > 0.0)
					entropy -= log(siteFreq)*siteFreq/LN_2;
			}
			
			res2.theData[k] = entropy;
		}
		else
		{
			long k4 = 0;
				
			for (long k3=k; k3<stateCount*patternCount; k3+=patternCount)
				res2.theData[k3] = siteFreqs->theData[mapper.lData[k4++]];
		}
		
		DeleteObject (siteFreqs);
	}
	

	_List colHeaders;
	
	_String	 rowLabels (128L, true);
	
	if (useEntropy)
	{
		for (long k5=0; k5<siteCount; k5++)
			res.theData[k5] = res2.theData[df->duplicateMap[k5]];
			
		rowLabels << "Entropy";
		rowLabels.Finalize();
		
		colHeaders && & rowLabels;

		rowLabels = _String("Character entropy for ") & 
					*(_String*)dataSetFilterNamesList(dataPartitions.lData[partIndex]);
					
	}
	else
	{
		
		for (long k9=1; k9<=siteCount; k9++)
		{
			rowLabels << ";Site ";
			rowLabels << _String(k9);
		}
			
		rowLabels.Finalize();
		for (long k5=0; k5<siteCount; k5++)
		{
			long dIndex = df->duplicateMap[k5];
			
			long k7 = dIndex;
			
			for (long k6=k5; k6<stateCount*siteCount; k6+=siteCount, k7+=patternCount)
				res.theData[k6] = res2.theData[k7];
		}
		
		for (long k8=0; k8<stateCount; k8++)
		{
			_String tts (df->ConvertCodeToLetters (df->CorrectCode(k8), stateSize));
			colHeaders && & tts;
		}
			
		colHeaders && & rowLabels;
		rowLabels = _String("Character Usage Across Sites In ") & 
					*(_String*)dataSetFilterNamesList(dataPartitions.lData[partIndex]);
			
	}
	
	long f = FindWindowByName (rowLabels);
	res.Transpose();

	_HYChartWindow * nc;
	
	if (f>=0)
	{
		nc = (_HYChartWindow*)windowObjectRefs (f);
		nc->SetTable (colHeaders,res);
	}
	else
	{
		nc = new _HYChartWindow (rowLabels, colHeaders, res, nil);
		checkPointer (nc);	
	}
	
	if (useEntropy)
		nc->SetChartType ("Bar Chart","Index","Entropy",true);
		
	nc->BringToFront();
	
}

//__________________________________________________________

void	_HYDataPanel::ShowAssociation (long partIndex, char options)
{
	_String			prompt	("Significance level [0,1]:"),
					cPrompt ("Resolve ambiguities randomly"),
					sig;
					
	bool			randomly = false;
					
	if (options == 2 || EnterStringDialogWithCheckbox (sig, prompt, cPrompt,randomly, (Ptr)this))
	{
		_Parameter		pvalue = sig.toNum();
		
		_DataSetFilter *df = (_DataSetFilter*)dataSetFilterList(dataPartitions.lData[partIndex]);
			
		long			stateSize		= df->GetUnitLength(),
						siteCount  		= df->GetFullLengthSpecies()/stateSize,
						patternCount	= df->NumberDistinctSites();
						
		if (options == 2 && (siteCount > 4 || siteCount < 2))
		{
			prompt = "Pairwise contigency tables only work with filters with 2,3 or 4 sites - otherwise too many windows would need to be open. Chop your data partition into smaller sections."; 
			ProblemReport (prompt, (Ptr)this);
			return;
		}

		_Matrix			patternRes		  (patternCount, patternCount, false, true);
		_List			patterns;		  
		

		_Parameter		* fv = new _Parameter [df->GetDimension (true)];
		checkPointer (fv);
		
		for (long site = 0; site < patternCount; site++)
		{
			_SimpleList* sitePat = df->CountAndResolve (site,fv,randomly);
			patterns << sitePat;
			DeleteObject (sitePat);
		}
		
		delete (fv);
		
		long	totalCount = patternCount*(patternCount-1)/2,
				totalDone  = 0,
				lastDone   = 0;
			
		if (options == 2)
		{	
			for (long site1=0; site1<siteCount-1; site1++)
			{
				for (long site2=site1+1; site2<siteCount;site2++)
				{
					if (df->duplicateMap[site1]!=df->duplicateMap[site2])
					{
						_List		labels,
									listLabels;
						_Matrix		* pc  = df->PairwiseCompare ((_SimpleList*)patterns(df->duplicateMap[site1]),
																 (_SimpleList*)patterns(df->duplicateMap[site2]), 
																 &labels);
						
						_SimpleList *rowLabels = (_SimpleList*)labels(1);
						for (long lc = 0; lc < rowLabels->lLength; lc++)
						{
							_String	   tts (df->ConvertCodeToLetters (df->CorrectCode (rowLabels->lData[lc]),stateSize));
							listLabels && & tts;
						}
						_String temp (128L, true);
						temp<<"C";
						
						rowLabels = (_SimpleList*)labels(0);
						for (long lc = 0; lc < rowLabels->lLength; lc++)
						{
							temp << ';';
							temp << df->ConvertCodeToLetters (df->CorrectCode (rowLabels->lData[lc]),stateSize);
						}
						
						temp.Finalize();
						listLabels && & temp;
						prompt = _String ("Contigency tables for sites ") & (site1+1) & " and " & (site2+1) & " in " &  *(_String*)dataSetFilterNamesList(dataPartitions.lData[partIndex]);
		
						_HYChartWindow * reportChart = (_HYChartWindow*) FindWindowByNameAndOpen (prompt);
						if (reportChart)
							reportChart->SetTable(listLabels, *pc);
						else
							reportChart = new _HYChartWindow (prompt, listLabels, *pc, nil);

						reportChart->BringToFront();	

						DeleteObject (pc);
					}
				}
			}
		}
		else
		{
			allocate_fexact_keys (4096,32);
			SetStatusLine ("Computing Associations");

			for (long site1=0; site1<patternCount-1; site1++)
			{						
				_SimpleList * sp1 = (_SimpleList*) patterns (site1);
				long		  sp1c = sp1->lData[sp1->lLength-1];
				
				for (long site2=site1+1; site2 < patternCount; site2++)
				{
					_SimpleList * sp2 = (_SimpleList*) patterns (site2);
					if ((sp1c > 1) && (sp2->lData[sp2->lLength-1] > 1))
					{
						_Matrix		* pc  = df->PairwiseCompare (sp1,sp2);
						if (options)
							patternRes.theData[site1*patternCount+site2] = pc->FisherExact (5.,80.,1.);
						else
						// use chi^2
						{
							long     rcnt = pc->GetHDim(),
									 ccnt = pc->GetVDim();
									 
									   
							_List* 	   rcc = pc->ComputeRowAndColSums();
							
							_Parameter chi2statistic = 0.0,
									   totals 		 = ((_Constant*)(*rcc)(2))->Value();
									   
							_Matrix*   rowSums 		 = (_Matrix*)(*rcc)(0),
								   *   columnSums 	 = (_Matrix*)(*rcc)(1);

							for (long rows = 0; rows < rcnt; rows++)
							{
								for (long columns = 0; columns < ccnt; columns ++)
								{
									_Parameter eTerm = rowSums->theData[rows]*columnSums->theData[columns]/totals,
										 	   aTerm = (pc->theData[rows*ccnt+columns]-eTerm);
										 	   
									chi2statistic += aTerm * aTerm / eTerm;
								}
							}		
							
							totals = rcnt*ccnt-rcnt-ccnt+1;
							
							_Constant val (chi2statistic),
									  degs (totals),
									  *res = (_Constant*)val.CChi2(&degs);
									  
							patternRes.theData[site1*patternCount+site2] = 1.-res->Value();
							DeleteObject (res);
							DeleteObject (rcc);
							
						}
						patternRes.theData[site2*patternCount+site1] = patternRes.theData[site1*patternCount+site2];
						DeleteObject (pc);
					}
					else
						patternRes.theData[site2*patternCount+site1] = patternRes.theData[site1*patternCount+site2] = 1.;
					totalDone ++;
					long testDone = (100.*totalDone)/totalCount;
					if (testDone > lastDone)
					{
						lastDone = testDone;
						SetStatusBarValue (lastDone,1.,0);
					}
				}
			}
		
			free_fexact_keys ();
			SetStatusBarValue (-1,1.,0);
			SetStatusLine ("Idle");
			
			_SimpleList	   matchedPairs;
			_List		   pValues;
			
			for (long site1=0; site1<siteCount-1; site1++)
			{
				for (long site2=site1+1; site2<siteCount;site2++)
				{
					if (df->duplicateMap[site1]!=df->duplicateMap[site2])
					{
						_Parameter locP = patternRes.theData[df->duplicateMap[site1]*patternCount+df->duplicateMap[site2]];
						if (locP <= pvalue)
						{
							_Constant locPC (locP);
							matchedPairs << site1+1;
							matchedPairs << site2+1;
							pValues && & locPC;
						}
					}
				}
			}
			if (pValues.lLength)
			{
				
				_List colHeaders;
				
				sig = "Site 1";
				colHeaders && & sig;
				
				sig = "Site 2";
				colHeaders && & sig;
				
				sig = "p-value";
				colHeaders && & sig;

				_Matrix res (pValues.lLength,3,false,true);
				
				for (long cntr = 0; cntr < pValues.lLength; cntr++)
				{
					res.Store (cntr,0,matchedPairs[2*cntr]);
					res.Store (cntr,1,matchedPairs[2*cntr+1]);
					res.Store (cntr,2,((_Constant*)(pValues(cntr)))->Value());
				}

				sig = _String("Significant association at level  ") & pvalue & " in " &
							*(_String*)dataSetFilterNamesList(dataPartitions.lData[partIndex]);
					
				long f = FindWindowByName (sig);
				_HYChartWindow * nc;
				if (f>=0)
				{
					nc = (_HYChartWindow*)windowObjectRefs (f);
					nc->SetTable (colHeaders,res);
				}
				else
				{
					nc = new _HYChartWindow (sig, colHeaders, res, nil);
					checkPointer (nc);	
				}
				nc->BringToFront();
				
				_List			mappedAssociations;
				_SimpleList		associatedSite;

				for (long kc = 0; kc < pValues.lLength; kc++)
				{
					long st1 = matchedPairs.lData[kc<<1],
						 st2 = matchedPairs.lData[(kc<<1)+1];
						 
					long f = associatedSite.BinaryFind (st1);
					
					_SimpleList* assSites;
					
					if (f<0)
					{
						f = associatedSite.BinaryInsert (st1);
						_SimpleList d (st1);
						mappedAssociations.InsertElement (&d, f, true);
					}

					assSites = (_SimpleList*)mappedAssociations (f);
					assSites->BinaryInsert(st2);
										
					f = associatedSite.BinaryFind (st2);
					
					if (f<0)
					{
						f = associatedSite.BinaryInsert (st2);
						_SimpleList d (st2);
						mappedAssociations.InsertElement (&d, f, true);
					}

					assSites = (_SimpleList*)mappedAssociations (f);
					assSites->BinaryInsert(st1);
				}
				
				_List	associatedClusters;
				
				for (long ma=0; ma<mappedAssociations.lLength; ma++)
				{
					_SimpleList* mappedAssoc = (_SimpleList*)mappedAssociations (ma);
					if (mappedAssoc->lLength == 2) // special case
					{
						_SimpleList assC;
						assC << *mappedAssoc;
						associatedClusters && & assC;
						
						mappedAssoc = (_SimpleList*)mappedAssociations (associatedSite. BinaryFind (mappedAssoc->lData[1]));
						
						mappedAssoc->Delete (mappedAssoc->BinaryFind (associatedSite.lData[ma]));
					}
					else
					{
						if (mappedAssoc->lLength > 2)
						{
							_SimpleList	   workSpace (*mappedAssoc);
							
							while (mappedAssoc->countitems() > 1)
							{
								long runningCount = 2;
								
								for (long k=0; k<mappedAssoc->lLength; k++)
								{
									if (mappedAssoc->lData[k] != associatedSite.lData[ma])
									{
										_SimpleList temp;
										_SimpleList *otherAssoc = (_SimpleList*) mappedAssociations (associatedSite.BinaryFind(mappedAssoc->lData[k]));
										temp.Intersect (workSpace,*otherAssoc);
										if ((temp.lLength < runningCount)||(k==mappedAssoc->lLength-1))
										{
											// reached end of the line
											associatedClusters && & workSpace;
											for (long kk = 0; kk < workSpace.lLength; kk++)
											{
												long li = associatedSite.BinaryFind(workSpace.lData[kk]);
												temp.Subtract (*(_SimpleList*) mappedAssociations (li), workSpace);
												temp.BinaryInsert (workSpace.lData[kk]);
												((_SimpleList*) mappedAssociations (li))->Duplicate (&temp);
											}
											workSpace.Duplicate (mappedAssoc);
											
											break;
										}
										else
										{
											runningCount ++;
											workSpace.Duplicate(&temp);
										}
											
									}
								}
							}
						}
					}
				}
				
				if (associatedClusters.lLength > 10)
				{
					prompt = _String ((long)associatedClusters.lLength) & " association clusters have been found. Are you sure you want to create a partition for each one of them?";
					if (!ProceedPrompt (prompt, (Ptr)this))
						return;
				}
				
				prompt = *(_String*)dataSetFilterNamesList(dataPartitions.lData[partIndex]) & "_cluster";
				
				bool	 doExclusions = false;
				if ((stateSize == 3)&&(df->theExclusions.lLength))
				{
					char 			a;
					
					LongToPartData (partData.lData[partIndex],a,doExclusions,patternCount);
					if (patternCount>=0)
						doExclusions = true;
				}
				
				
				for (long cc = 0; cc < associatedClusters.lLength; cc++)
				{
					_SimpleList clusterSites;
					_SimpleList * cls = (_SimpleList*)associatedClusters (cc);
					cPrompt = prompt;
					for (long sc=0; sc < cls->lLength; sc++)
						//for (long uc = 0; uc < stateSize; uc++)
						{
							clusterSites << df->theOriginalOrder[(cls->lData[sc]-1)*stateSize];
							cPrompt = cPrompt & "_" & _String(cls->lData[sc]);
						}
							
					
							
					CreatePartition (clusterSites, stateSize, false, & cPrompt);
					if (doExclusions)
					{
						SetCodonExclusions ((_DataSetFilter*)dataSetFilterList(dataPartitions.lData[dataPartitions.lLength-1]),dataPartitions.lLength-1,patternCount);
						partData.lData[dataPartitions.lLength-1] = PartDataToLong (0,false,patternCount);
					}
				}
			}
			else
			{
				sig = _String ("No association with significance at or better than ") & pvalue;
				ProblemReport (sig, (Ptr)this);
			}
		}
	}
	
}

//__________________________________________________________

void	_HYDataPanel::ShowCodonUsage (long partIndex)
{
	// go filter by filter; take each codon filter, apply properties and translate	
	_DataSetFilter *df = (_DataSetFilter*)dataSetFilterList(dataPartitions.lData[partIndex]);
	_Matrix		   *codonUsage = df->HarvestFrequencies (3,3,false);
	
	char 			offset;
	bool 			rev;
	long 			code,
					cCount = df->NumberSpecies ()*df->theOriginalOrder.lLength/3;
	
	LongToPartData  (partData.lData[partIndex], offset, rev, code);
	
	_SimpleList*	gCode = (_SimpleList*)(*((_List*)geneticCodes (code)))(2);
	
	
	code = 0;
	
	char			buffer [1024];
	_String			rec;
	
	sprintf			(buffer,"\nCodon Usage Report\n\nAltogether : %d codons\n\nBy codon:\n\n", cCount);	
	BufferToConsole (buffer);
	
	char			nucs[] = "ACGT";
	
	for	(long k=0; k<gCode->lLength; k++)
	{
		if (gCode->lData[k] == 10)
			code++;
		else
		{
			_Parameter cUsage = codonUsage->theData[k];
			if (cUsage)
			{
				CodeTo3AA (rec, k, gCode);
				sprintf (buffer,"%c%c%c(%s) : %10d (%10.4g %%)\n", nucs[k/16], nucs[(k%16)/4], nucs[k%4], rec.getStr(), (long)(cUsage*cCount), 100.*cUsage);
				BufferToConsole (buffer);
			}
		}
	}
	
	BufferToConsole ("\n\nBy aminoacid:\n\n");	
	
	for	(long kk=0; kk<=20; kk++)
	{
		if (kk != 10)
		{
			bool   t = true;
			code = 0;
			for (long k=0; k<gCode->lLength; k++)
			{
				_Parameter cUsage = codonUsage->theData[k];
				if ((gCode->lData[k] == kk)&&(cUsage))
				{
					if (t)
					{
						CodeTo3AA (rec, k, gCode);
						sprintf (buffer,"%s=>\t", rec.getStr());
						t = false;
						BufferToConsole (buffer);
			
					}
					sprintf (buffer,"%c%c%c:%6d (%6.4g %%)\t", nucs[k/16], nucs[(k%16)/4], nucs[k%4], (long)(cUsage*cCount), 100.*cUsage);
					BufferToConsole (buffer);	
				}
				else
					if (gCode->lData[k] == 10)
						code++;
			}
			if (!t)
				NLToConsole();
		}
	}
	
	NLToConsole();
	DeleteObject    (codonUsage);
	
}

//__________________________________________________________

void	_HYDataPanel::ExecuteProcessor (long procID)
{
	if (procID<dataPanelProcessors.lLength)
	{
		_String * filePathToProc = (_String*)dataPanelProcessors(procID);
		FILE *	  thisFile		 = doFileOpen (filePathToProc->sData,"rb");
		if (thisFile)
		{	
			_HYSequencePane* sp		   = (_HYSequencePane*)	components (0);
		    _HYTable*	 	 pList	   = (_HYTable*)		GetObject  (5);
			_Matrix			 *site_sel = nil,
						     *seqs_sel = nil,
							 *dflt_sel = nil,
							 *dfdm_sel = nil;
							 
			_AssociativeList *genc_sel = new _AssociativeList,
							 *excl_sel = new _AssociativeList;
							 
			long			 createdFilterID = -1;
							 
			_DataSetFilter*  passedDataFilter = nil;
			
			if (sp->vselection.lLength)
				seqs_sel = new _Matrix (sp->vselection);
			else
				if (sp->selection.lLength)
					site_sel = new _Matrix (sp->selection);
				
			if (!site_sel)
				site_sel = new _Matrix;
			if (!seqs_sel)
				seqs_sel = new _Matrix;
				
			_SimpleList	   uti,
						   ld;
						   
			pList->GetRowSelection (uti);
			if (uti.lLength > 0)
			{
				if (uti.lLength == 1)
					passedDataFilter = (_DataSetFilter*)((_DataSetFilter*)dataSetFilterList(dataPartitions.lData[uti.lData[0]]))->makeDynamic();
				else
					if (dataWrapper)
						passedDataFilter = (_DataSetFilter*)dataWrapper->makeDynamic();
					else
					{
						checkPointer (passedDataFilter = new _DataSetFilter);
						_SimpleList  blank, blank2;
						passedDataFilter->SetFilter ((_DataSet*)dataSetList (dataSetID),1,blank,blank2,0);
					}
					
					
				_List	selectedFilters;
				for (long ti = 0; ti < uti.lLength; ti++)
				{
					selectedFilters << dataSetFilterNamesList(dataPartitions.lData[uti.lData[ti]]);
					_Matrix			* gCode = nil;
					_DataSetFilter  * dsf = (_DataSetFilter*)dataSetFilterList(dataPartitions.lData[uti.lData[ti]]);
					_FString		* excl = nil;
					if (dsf->GetUnitLength()==3)
					{
						char 			offset;
						bool 			rev;
						long 			code;
						LongToPartData  (partData.lData[ti], offset, rev, code);
						gCode = new _Matrix(*(_SimpleList*)(*((_List*)geneticCodes (code)))(2));
						excl  = new _FString (GetExclusionsFromCode (code));
					}
					else
					{
						gCode = new _Matrix;
						excl  = new _FString;
					}
						
					ld << dsf->GetUnitLength();
					_String   *keyString = new _String (ti);
					_FString  aKey (keyString);
					genc_sel->MStore (&aKey, gCode, false);
					if (excl)
						excl_sel->MStore (&aKey, excl, false);
				}
				dflt_sel = new _Matrix (selectedFilters);
			}
			else
			{
				if (dataWrapper)
					passedDataFilter = (_DataSetFilter*)dataWrapper->makeDynamic();
				else
				{
					checkPointer (passedDataFilter = new _DataSetFilter);
					_SimpleList  blank;
					passedDataFilter->SetFilter ((_DataSet*)dataSetList (dataSetID),1,uti,blank,0);
				}
				dflt_sel = new _Matrix;
				ld << 1;
			}
			
			if (passedDataFilter)
				createdFilterID = AddFilterToList (datapanelProcDF, passedDataFilter, true);
			
			dfdm_sel = new _Matrix (ld);
			
			setParameter (datapanelProcHS,site_sel,false);
			setParameter (datapanelProcVS,seqs_sel,false);
			setParameter (datapanelProcSF,dflt_sel,false);
			setParameter (datapanelProcGC,genc_sel,false);
			setParameter (datapanelProcUS,dfdm_sel,false);
			setParameter (datapanelProcEX,excl_sel,false);
			setParameter (datapanelProcDS,new _FString (*(_String*)dataSetNamesList (dataSetID), false),false);
				
			_String buffer (thisFile);
			fclose (thisFile);
			
			long   g = batchLanguageFunctionNames.lLength;
			
			_String dpp = *filePathToProc;
			PushFilePath (dpp);				

			_ExecutionList   thisList;
			terminateExecution = false;
			
			thisList.BuildList (buffer);
			thisList.ExecuteAndClean(g);
			terminateExecution = false;
			PopFilePath ();				
		
			if ((excl_sel = CheckAssociativeListArg (&datapanelProcRF))) // returned some filters to make 
			{
				_List * returnedSpecs = excl_sel->GetKeys();
				for (g = 0; g < returnedSpecs->lLength; g++)
				{
					_String * proposedName = (_String*)(*returnedSpecs) (g);
					if ((genc_sel = (_AssociativeList*)excl_sel->GetByKey(*proposedName, ASSOCIATIVE_LIST)))
					{ // should have a _row_ matrices with key "SITES" (for now)
						_String tID (*proposedName);
						tID.ConvertToAnIdent();
						if (tID.IsValidIdentifier (false))
						{
							_String aKey ("SITES");
							site_sel = (_Matrix*)genc_sel->GetByKey (aKey, MATRIX);
							if (site_sel)
							{
								site_sel->ConvertToSimpleList (uti);
								CreatePartition (uti, 1, false, &tID); 
							}							
						}
					}
				}
			}
			
			DeleteVariable (datapanelProcHS);
			DeleteVariable (datapanelProcVS);
			DeleteVariable (datapanelProcSF);
			DeleteVariable (datapanelProcGC);
			DeleteVariable (datapanelProcUS);
			DeleteVariable (datapanelProcEX);
			DeleteVariable (datapanelProcDS);
			DeleteVariable (datapanelProcRF);
			
			if (createdFilterID>=0)
				KillDataFilterRecord (createdFilterID, true);
		}
		else
		{
			_String eMsg = _String("Problem reading file:") & *filePathToProc;
			ProblemReport (eMsg);
		}
	}
}

//__________________________________________________________

bool	_HYDataPanel::EditPartitionProperties (long index, char changeType)
{
	if (dataSetID>=0)
	{
		if ((index<0)||(index>=dataPartitions.lLength)) return false;
		
		long			 dfID = dataPartitions[index], geneticCode = 0;
		_DataSetFilter*  thisDF = (_DataSetFilter*)dataSetFilterList(dfID);
		_String		  	 currentName = *((_String*)dataSetFilterNamesList(dfID)),
					  	 partitionInfo;
		_HYColor	  	 currentColor = LongToHYColor(partitionColors.lData[index]),
					 	 oldColor = currentColor;
		_HYTable*	 	 pList = (_HYTable*)GetObject (5);
		bool		 	 direction = true;
		char		 	 readFrame = 0;
		
		LongToPartData (partData.lData[index],readFrame,direction,geneticCode);

	
		if (changeType)
			partitionInfo = "Changing partition data type";
		else
		{
			partitionInfo = _String (thisDF->NumberSpecies())& " sequences, "&
							_String ((long)(thisDF->theOriginalOrder.lLength/thisDF->GetUnitLength()))&
							" total data sites, with "& _String ((long)thisDF->theFrequencies.lLength)
							&" unique patterns. ";
									
			if (dataType&HY_DATAPANEL_NUCDATA || dataType&HY_DATAPANEL_BINARYDATA )
			{
				switch(thisDF->GetUnitLength())
				{
					case 2:
						partitionInfo = partitionInfo & 
										((dataType&HY_DATAPANEL_NUCDATA)?dinucDataType:disequonDataType);
						break;
					case 3:
						partitionInfo = partitionInfo & codonDataType;
						break;
					default:
						partitionInfo = partitionInfo & 
										((dataType&HY_DATAPANEL_NUCDATA)?nucDataType:binaryDataType);
				}
				
				partitionInfo = partitionInfo & " partition.";
			}
			else
				if (dataType&HY_DATAPANEL_PROTDATA)
					partitionInfo = partitionInfo & proteinDataType &" partition.";
		}
		
		bool		okcancel = false;
		
		_HYPartitionDialog* sd = new _HYPartitionDialog (partitionInfo,&currentName,&currentColor,&direction,&readFrame,&geneticCode,dfID,&okcancel,
							changeType?changeType-1:thisDF->GetUnitLength()-1, (Ptr)this);
		sd->Activate();
		while (windowObjectRefs.Find ((long)sd)>=0)
			handleGUI();
			
		//if (PartitionEditDialog(partitionInfo,currentName,currentColor,direction,readFrame,geneticCode,dfID,
		//	changeType?changeType-1:thisDF->GetUnitLength()-1))
		if (okcancel)
		{
			long cellID = pList->horizontalSpaces.lLength*index;
			if (!(((_String*)dataSetFilterNamesList(dfID))->Equal(&currentName)))
			{
				// rename the filter
				((_String*)dataSetFilterNamesList(dfID))->Duplicate(&currentName);
				pList->SetCellData (&currentName, index, DF_ID_COLUMN, pList->cellTypes[cellID+DF_ID_COLUMN], true);
				pList->_MarkCellForUpdate (cellID+DF_ID_COLUMN);
			}
			if (!(currentColor==oldColor))
			{
				partitionColors.lData[index] = HYColorToLong (currentColor);
				((_SimpleList*)pList->GetCellData (DF_COLOR_COLUMN,index))->lData[0] = partitionColors.lData[index];
				pList->_MarkCellForUpdate (cellID+DF_COLOR_COLUMN);
				BuildThermometer();
				BuildMarksPane();
			}
			long newCode = PartDataToLong (readFrame,direction,geneticCode);
			if (!changeType)
			{
				if (newCode!=partData.lData[index])
				{
					PurgeLFFilter (index);
					if (((partData.lData[index]&HY_DATAPANEL_CODEMASK)>>16)!=geneticCode)
						SetCodonExclusions (thisDF,dfID,geneticCode);
		
					partData.lData[index] = newCode;
					if (addedLines&HY_DATAPANEL_TRANSLATION)
						UpdateTranslationString (partitionInfo,translatedSequence,true);
				}
			}
			else
				partData.lData[index] = newCode;

			tainted = true;
			return true;
		}
	}
	return false;
}

//__________________________________________________________

void	_HYDataPanel::InputPartitionString (void)
{
	if (dataSetID>=0)
	{
		_String  partSpec;
		if (EnterStringDialogWithCheckbox (partSpec,enterPartitionString,useOneBasedIndices,useOneBased, (Ptr)this))
		{
			_SimpleList newPart;
			_DataSet* thisDS = (_DataSet*)dataSetList (dataSetID);
			partSpec.Insert('"',0);
			partSpec.Insert('"',-1);
			thisDS->ProcessPartition (partSpec,newPart,true);
			if (useOneBased)
				for (long k=0; k<newPart.lLength; k++)
					newPart.lData[k]--;
					
			if (newPart.lLength)
				CreatePartition (newPart);
			else
				ProblemReport(invalidPartString,(Ptr)this);
		}
	}
}

//__________________________________________________________

void	_HYDataPanel::SetDataSetReference (_String& varName, _SimpleList* specFilter)
{
	long f = dataSetNamesList.Find(&varName);
	dataType = 0;
	if (f>=0)
	{
		_DataSet		* theDS	   = (_DataSet*) dataSetList (f);
		_HYSequencePane * sp	   = (_HYSequencePane*)GetObject (0);	
		_HYSequencePane * sp2	   = (_HYSequencePane*)GetObject (4);	
		dataSetID = f;
		sp->columnStrings.Clear();
		sp->selection.Clear();
		sp->speciesIndex.Clear();
		
		long upToC = theDS->NoOfColumns();
		siteAssignments.RequestSpace(upToC);
		overlaps.RequestSpace (upToC);
		sp->columnStrings.RequestSpace (upToC);
		for (long k=0;k<upToC;k++)
		{
			sp->InsertColumn (theDS->GetSite(k),-1,false);
			siteAssignments<<-1;
			overlaps<<0;
		}
		sp->SetHeaders 	 (&theDS->GetNames(),false);
		if (specFilter)
		{
			sp->speciesIndex.Clear();
			sp->speciesIndex.Duplicate (specFilter);
		}
		if (theDS->GetTT()->IsStandardNucleotide())
		{
			sp->SetCharColor ('A',nucDefaults[0],false);
			sp->SetCharColor ('C',nucDefaults[1],false);
			sp->SetCharColor ('G',nucDefaults[2],false);
			sp->SetCharColor ('T',nucDefaults[3],false);
			sp->SetCharColor ('U',nucDefaults[4],false);
			sp2->SetCharColor ('A',nucDefaults[0],false);
			sp2->SetCharColor ('C',nucDefaults[1],false);
			sp2->SetCharColor ('G',nucDefaults[2],false);
			sp2->SetCharColor ('T',nucDefaults[3],false);
			sp2->SetCharColor ('U',nucDefaults[4],false);
			dataType = HY_DATAPANEL_NUCDATA;
		}
		else
			if (theDS->GetTT()->IsStandardAA())
			{
				for (long k=0;k<20;k++)
				{
					sp->SetCharColor  (aminoAcidOneCharCodes.sData[k],aaDefaults[k],false);
					sp2->SetCharColor (aminoAcidOneCharCodes.sData[k],aaDefaults[k],false);
				}
				dataType = HY_DATAPANEL_PROTDATA;
			}
			else
				if (theDS->GetTT()->IsStandardBinary())
				{
					sp->SetCharColor (binaryOneCharCodes.sData[0],nucDefaults[0],false);
					sp->SetCharColor (binaryOneCharCodes.sData[1],nucDefaults[1],false);
					dataType = HY_DATAPANEL_BINARYDATA;
				}
				
		dataSetName = (_String*)dataSetNamesList (dataSetID);
		_String     newWindowName ("DataSet ");
		newWindowName = newWindowName & varName;
		SetTitle (newWindowName);
		BuildDataPartitions();
		_HYRect  screenRect = GetScreenDimensions();
		
		SetWindowRectangle (0,0,screenRect.bottom-250,screenRect.right-20);
		
		#ifdef __MAC__
			SetPosition (5,40);
			screenRect.left 	= 	5;
			screenRect.top 		= 	bottom + 52;
			screenRect.bottom 	-= 	10;
			screenRect.right 	= 	right+5;
		#else
			SetPosition (5,5);
			screenRect.left 	= 5;
			screenRect.top 		= bottom + 52;
			screenRect.bottom	-= 10;
			screenRect.right 	= right+5;
		#endif
		if (doAutoConsoleMove)
			MoveConsoleWindow(screenRect);
	}
}


//__________________________________________________________

void	_HYDataPanel::Update (Ptr p)
{
	_HYTWindow::Update(p);
	_HYCanvas* theCanvas = (_HYCanvas*)GetObject (1);
	#ifdef __MAC__
	forceUpdateForScrolling=true;
	theCanvas->_MarkForUpdate();
	forceUpdateForScrolling=false;
	#endif
	_PaintThermRect();
}

//__________________________________________________________

void	_HYDataPanel::Activate (void)
{
	// check the validity of the trees attached
	// to the data partitions
	long 	k, f;
	_String	errMsg;
	_TheTree * thisTree;
	if (lfID>=0)
	{
		for (k=0; k<dataPartitions.lLength; k++)
		{
			f = treeVarReferences.lData[k];
			if (f>=0)
			{
				thisTree = (_TheTree*)LocateVar(f);
				if (!thisTree->AllBranchesHaveModels(((_DataSetFilter*)dataSetFilterList(dataPartitions.lData[k]))->GetDimension(true)))
				{	
					errMsg = _String ("Likelihood function was killed b/c some of the tree models were edited/deleted"); 
					ReportWarning (errMsg);
					terminateExecution = false; 
					PurgeLF();
					break;
				}
			}
		}
	}
	
	for (k=0; k<dataPartitions.lLength; k++)
	{
		f = treeVarReferences.lData[k];
		if (f>=0)
		{
			thisTree = (_TheTree*)LocateVar(f);
			_PMathObj tc = thisTree->TipCount();
			if (tc->Value()!=((_DataSetFilter*)dataSetFilterList(dataPartitions.lData[k]))->NumberSpecies() &&
				(!(tc->Value()==1 && ((_DataSetFilter*)dataSetFilterList(dataPartitions.lData[k]))->NumberSpecies()==2)))
			{
				if (errMsg.sLength)
					errMsg = errMsg&',';
				errMsg = errMsg & *thisTree->GetName() &'(' & _String((long)tc->Value()) &')';
				treeVarReferences.lData[k]=-1;
				_List dummy;
				RefreshPartRow (dummy,k,true,false);
				
			}
			DeleteObject (tc);
		}
	}
	
	if (errMsg.sLength)
	{
		errMsg = _String("HyPhy detected that the following trees had their number of leaves changed since this window was last active:")
				 &errMsg&".\n";
		StringToConsole (errMsg);
	}
	
	_HYTWindow::Activate();
	UpdatePartitionOperations();
	_UpdateSelectionChoices(IsSelectionNonEmpty());
}

//__________________________________________________________

void	_HYDataPanel::Paint (Ptr p)
{
	_HYTWindow::Paint(p);
	_HYCanvas* theCanvas = (_HYCanvas*)GetObject (1);
	#ifdef __MAC__
	forceUpdateForScrolling=true;
	theCanvas->_MarkForUpdate();
	forceUpdateForScrolling=false;
	#endif
	_PaintThermRect();
}


//__________________________________________________________

bool	_HYDataPanel::BuildDataPanel (void)
{
	return true;
}

//__________________________________________________________

void	_HYDataPanel::ComputeLikelihoodFunction (long option)
{
	if (lfID>=0)
	{
		if (lfID != lockedLFID)
		{
			_LikelihoodFunction *lf = (_LikelihoodFunction*)likeFuncList (lfID);
			stashParameter (likefuncOutput,option,true);
			_String* lfStr = (_String*)lf->toStr();
			_String outString (256L,true);
			outString << '\n';
			outString << lfStr;
			outString << '\n';
			outString.Finalize();
			StringToConsole (outString);
			DeleteObject (lfStr);
			stashParameter (likefuncOutput,0.0,false);
		}
		else
			ProblemReport (lfCantKillWarning, (Ptr)this);
	}
}

//__________________________________________________________

void	_HYDataPanel::SimulateDataSet (long option, bool ancestors)
{
	if (lfID>=0)
	{
		if (lfID != lockedLFID)
		{
			_LikelihoodFunction *lf = (_LikelihoodFunction*)likeFuncList (lfID);
			/*if (ancestors)
			{
				if (lf->GetCategoryVars().lLength)
				{
					_String errMsg ("Can't reconstruct ancestors in models with rate variation");
					ProblemReport (errMsg, (Ptr)this);
					return;
				}
			}*/
			_DataSet * target = new _DataSet;
			_List	   emptyList;
			
			if (!ancestors)
				lf->Simulate (*target, emptyList);
			else
				lf->ReconstructAncestors (*target);
				
			if (option==1) // save to file
			{
				dialogPrompt = "Save simulated data set to";
				_String pathName (WriteFileDialogInput());
				if (pathName.sLength)
				{
					FILE * outfile = doFileOpen (pathName.getStr(),"w");
					if (outfile)
					{
						_DataSetFilter dsf;
						_SimpleList e1, e2;
						dsf.SetFilter (target,1,e1,e2,false);
						dsf.toFileStr(outfile);
						fclose (outfile);
					}
					else
						ReportWarning (_String("Could not open file '") & pathName & "' for writing.");
				}
			}
			else
			{
				if (option==2) // save to file
				{
	 				_String promptString ("10"), folderName ("Number of replicates:");
	 				if (EnterStringDialog (promptString,folderName, (Ptr)this))
	 				{
	 					long nIterates = promptString.toNum();
	 					if (nIterates>0)
	 					{
							dialogPrompt = "Save simulated data sets to";
							folderName = ChooseAFolder(dialogPrompt);
							if (folderName.sLength)
								for (long k=1; k<=nIterates; k++)
								{
									if (k>1)
									{
										DeleteObject (target);
										target = new _DataSet;
										lf->Simulate (*target, emptyList);
									}
									_String fileName = folderName & (*(_String*)dataSetNamesList(dataSetID))&("_sim.") & _String (k);
									FILE * outfile = doFileOpen (fileName.getStr(),"w");
									if (outfile)
									{
										_DataSetFilter dsf;
										_SimpleList e1, e2;
										dsf.SetFilter (target,1,e1,e2,false);
										dsf.toFileStr(outfile);
										fclose (outfile);
									}
									else
									{
										ReportWarning (_String("Could not open file '") & fileName & "' for writing.");
										break;
									}
								}
						}
					}
				}
				else
				{
					_String newDataString = (*(_String*)dataSetNamesList(dataSetID));
					if (ancestors)
						newDataString = newDataString & "_ancestors";
					else
						newDataString = newDataString & "_sim";
						
					AddDataSetToList (newDataString, target);
					_HYDataPanel* myTWindow = new _HYDataPanel(newDataString,newDataString);
					myTWindow->BringToFront();
					//myTWindow->Show();	
					return;
				}		
			}
			DeleteObject (target);
		}
		else
			ProblemReport (lfCantKillWarning, (Ptr)this);
	}
}
//__________________________________________________________

void	_HYDataPanel::OptimizeLikelihoodFunction (void)
{
	if (lfID>=0)
	{
		if (lfID!=lockedLFID)
		{
			_LikelihoodFunction *lf = (_LikelihoodFunction*)likeFuncList (lfID);
			
			long				 startingLFEvalsCount = likeFuncEvalCallCount;
			_PMathObj			 startingTime		  = _Constant (0.0).Time(),
								 endingTime;
			
			ToggleAnalysisMenu (true);
			StartBarTimer ();
			_Matrix* res = lf->Optimize();
			StopBarTimer  ();
			ToggleAnalysisMenu (false);
			terminateExecution = false;
			
			startingLFEvalsCount = likeFuncEvalCallCount - startingLFEvalsCount;
			endingTime = _Constant (0.0).Time();
			
			DeleteObject (res);
			ComputeLikelihoodFunction (2);
			char	buffer [1024];
			sprintf (buffer,"\n\tTime taken = %g seconds\n\tLF evaluations/second = % g\n", 
						endingTime->Value()-startingTime->Value(), 
						startingLFEvalsCount/(endingTime->Value()-startingTime->Value()));
			BufferToConsole (buffer);

			postChangeLFEvent (GetID(),lfID);
			if (autoPopLFTable)
				DisplayParameterTable();
				
			ReportAnalysisAsFinished (empty);
				
		}
		else
			ProblemReport (lfCantKillWarning, (Ptr)this);
	}
}

//__________________________________________________________

void	_HYDataPanel::RestoreSavedLFs (void)
{
	if (lfID>=0)
	{
		long k = LocateVarByName (savedLFMatrix);
		if (k>=0)
		{
			_Variable* saveLF = FetchVar (k);
			if (saveLF->ObjectClass() == MATRIX)
			{
				_Matrix* saveLFMatrix = (_Matrix*)saveLF->varValue;
				if (saveLFMatrix->GetVDim()==2)
				{
					for (k=0; k<saveLFMatrix->GetHDim(); k++)
					{
						_PMathObj p1 = ((_Formula**)saveLFMatrix->theData)[2*k]->Compute(),
								  p2 = ((_Formula**)saveLFMatrix->theData)[2*k+1]->Compute();
								  
						if ((p1->ObjectClass() == STRING)&&(p2->ObjectClass() == STRING))
						{
							savedLFNames  << ((_FString*)p1)->theString;
							savedLFStates << ((_FString*)p2)->theString;
						}
					
					}
				}
			}
			DeleteVariable (savedLFMatrix);
		}
	}
}
//__________________________________________________________

bool	_HYDataPanel::GenerateGoodPartitions (_SimpleList& goodPartitions)
{
	bool retval = false;
	long k = 0;
	for (; k<dataPartitions.lLength; k++)
	{		
		if (modelReferences.lData[k] == 0)
			continue;

		if (treeVarReferences.lData[k]<0)
		{
			if (treeVarReferences.lData[k]!=-3)
				continue;
			goodPartitions << k;
			retval = true;
		}
		else
			goodPartitions << k;
	}
	return retval;
}

//__________________________________________________________

void	_HYDataPanel::InferTopologies (bool useConstr)
{
	if (!cantDeleteLF)
	{
		bool	needToDoFixed = false;
		_String pathToFiles = baseDirectory&"TopologyInference";
		
		_List   receptacle,
				fNames,
				pathNameList;
		
		_SimpleList l1,
					l2,
					l3,
					l4;
					
		if (useConstr)
		{
			topConstr 					   = empty;
			_String   		dummyTree 	   = _String ("Tree ")&inferenceDummyTree&" = (0,1);";
			_ExecutionList	spawnDummyTree (dummyTree);
			spawnDummyTree.Execute ();
			
			long 			treeID = LocateVarByName (inferenceDummyTree);
			
			if (treeID < 0)
			{
				pathToFiles = "Internal error in InferTopologies - failed to spawn dummy tree.";
				ProblemReport (pathToFiles);
				return;
			}
			
			#ifndef USE_AVL_NAMES
				treeID = variableReindex.lData[treeID];
			#else
				treeID = variableNames.GetXtra(treeID);
			#endif
			
			_SimpleList sillySubset;
			
			sillySubset << 0;
			sillySubset << 1;
			
			BuildLikelihoodFunction (&inferenceDummyLF, &sillySubset, treeID);
			
			if (lfID < 0)
			{
				pathToFiles = "Internal error in InferTopologies - failed to spawn dummy LF.";
				ProblemReport (pathToFiles);
				return;
			}
			else
			{
				_LikelihoodFunction *me = (_LikelihoodFunction*) likeFuncList (lfID);
				
				_SimpleList	   *myTrees	= &me->GetTheTrees(),
								glVarIDs,
								localVarIDs;
								
				
				_List		   treeNames,
							   globalVars,
							   localVars,
							   localVarID,
							   templateConstraints,
							   templateConstraintStrings;
							  
				dummyTree =  blMolClock &" tree ID or node ID, variable)";
				templateConstraintStrings && & dummyTree;
				dummyTree = "Molecular Clock";
				templateConstraints && & dummyTree;
				
				dummyTree =  blReplicate & "\"this1.?.varID:=this2.?.varID\", tree or node id, tree or node id)";
				templateConstraintStrings && & dummyTree;
				dummyTree = "Replicate Constraints";
				templateConstraints && & dummyTree;

				for (treeID = 0; treeID < myTrees->lLength; treeID ++)
				{
					long	   counter = 0;
					_String    varName;
					
					_TheTree*  aTree = (_TheTree*)LocateVar (myTrees->lData[treeID]);
					treeNames << aTree->GetName();
					glVarIDs.Clear();
					{
						_AVLList gav (&glVarIDs);
						aTree->ScanForGVariables (gav, gav);
						gav.ReorderList ();
					}
					for (counter = 0; counter < glVarIDs.lLength; counter++)
					{
						_Variable* aVar = LocateVar (glVarIDs.lData[counter]);
						if (!aVar->IsCategory())
							globalVars << aVar->GetName();
						
					}
					_SimpleList	 localsForTree;
					aTree->FindScalingVariables (localsForTree);
					
					for (counter = 0; counter < localsForTree.lLength; counter++)
					{
						varName	=  * LocateVar (localsForTree.lData[counter])->GetName();
						localVarID && & varName;
						varName	=  varName & " [" & *aTree->GetName() & ']';
						localVars && & varName;
					}
				}
								
				KillLFRecord (lfID);
				lfID = -1;
				
				bool resB = false;
				
				_HYFont def;
				def.face  = "System Font";
				def.size  = 12;
				def.style = HY_FONT_PLAIN;
				
				_String constrPrompt ("Apply the following tree constraints:");
				
				_HYInferenceConstraints * sd = 
						new _HYInferenceConstraints (constrPrompt,localVars, globalVars, treeNames, templateConstraints,
																  localVarID, globalVars, treeNames, templateConstraintStrings,
																  def,&resB,&topConstr, (Ptr)this);
				sd->Activate();
				while (windowObjectRefs.Find ((long)sd)>=0)
					handleGUI();
					
				if (!resB) return;
				
			}

			//_String constrPrompt ("Apply the following tree constraints:");
			//if (!EnterStringDialog (topConstr, constrPrompt, (Ptr)this))
				//return;
		}

		char del = ScanDirectoryForFileNames (pathToFiles,receptacle,false);
		
		for (long k=0; k<receptacle.lLength;k++)
		{
			FILE * thisFile = doFileOpen (((_String*)receptacle(k))->sData,"r");
			if (thisFile)
			{	
				fclose (thisFile);
				_String   fName = *(_String*)receptacle(k);
				long	  kk = fName.FindBackwards (del,0,-1);
				fName.Trim (kk+1,-1);
				_List	   thisFile;
				
				thisFile && & fName;
				thisFile && & fName;
				
				l1 << fNames.lLength;
				l4 << k;
				fNames && & thisFile;
				_String  tts (*(_String*)receptacle(k),0,kk);
				pathNameList && & tts;
			}
		}
		
		if (fNames.lLength == 0)
		{
			pathToFiles = "Could not infer topologies, because inference modules could not be found. Check you distribution for missing files.";
			ProblemReport (pathToFiles, (Ptr)this);	
			return;
		}
		
		if (omittedSeqs.lLength)
		{
			pathToFiles = "Can't infer topologies when some sequences of the original data set have been omitted. Consider saving the smaller data set as a separate file, and running topology reconstuction on the smaller data set.";
			ProblemReport (pathToFiles, (Ptr)this);	
			return;
		}
		
		l2 << 0;
		l2 << 1;
		
		long		choice = HandleListSelection (fNames, l2, l1, "Inference Method", l3, 1);
		
		if (choice>=0)
		{
			_String infFile (1024L,true);
			infFile << inferenceNofSeqs;
			infFile << "=";
			infFile << _String ((long) ((_DataSet*)dataSetList (dataSetID))->GetNames().lLength);
			infFile << ";\n";
			infFile << inferenceDWID;
			infFile << "=\"";
			infFile << GetTitle();
			infFile << "\";\n";
			
			pathNames << pathNameList (choice);
			
			_SimpleList	 stashedTrees,
						 stashedTreeIDs;
						 
			long		 g;
						 

			FILE * thisFile = doFileOpen (((_String*)receptacle (l4.lData[choice]))->getStr(),"rb");
			if (!thisFile)
			{
				pathToFiles = "Could not infer topologies, because the selected inference module could not be read. It may have been deleted.";
				ProblemReport (pathToFiles, (Ptr)this);	
				infFile.Finalize();
				return;
			}
			infFile << _String (thisFile);
			fclose (thisFile);
			infFile.Finalize();
			infFile = infFile.Replace (inferenceDSID, *(_String*)dataSetNamesList (dataSetID), true);
			
			
			g = batchLanguageFunctionNames.lLength;

			_ExecutionList infL (infFile);
			
			if (terminateExecution)
			{
				terminateExecution = false;
				return;
			}
			
			SetLockState 		(true);
						
			ToggleAnalysisMenu 	(true);
			StartBarTimer 		();
			infL.ExecuteAndClean(g);
			StopBarTimer 		();
			ToggleAnalysisMenu 	(false);
				
			pathNames.Delete 	(pathNames.lLength-1);
			
			if (terminateExecution)
				terminateExecution = false;
			else
			{
				if (inferCache.lLength)
				{
					long iCount = 0;
					_HYTable*		pl = (_HYTable*)GetObject (5);
						
					KillLFRecord (lfID,false);
					
					_String		 alteredTopConstr = topConstr;
					
					for (g=0; g<dataPartitions.lLength; g++)
					{
						if (treeVarReferences[g] == -3)
						{
							#ifndef USE_AVL_NAMES
							treeVarReferences[g] = variableReindex.lData[LocateVarByName(*(_String*)inferCache(iCount))];
							#else
							treeVarReferences[g] = variableNames.GetXtra(LocateVarByName(*(_String*)inferCache(iCount)));
							#endif
							
							pl->SetCellData ((_String*)inferCache(iCount),g,DF_TREE_COLUMN,pl->cellTypes.lData[g*pl->horizontalSpaces.lLength+DF_TREE_COLUMN],true);
							
							pl->_MarkCellForUpdate (g*pl->horizontalSpaces.lLength+DF_TREE_COLUMN);
							
							if (topConstr.sLength)
								alteredTopConstr = alteredTopConstr.Replace (*(_String*)dataSetFilterNamesList (inferCacheDF.lData[iCount]),
														  *(_String*)dataSetFilterNamesList (dataPartitions.lData[g]),
														  true);

							KillDataFilterRecord (inferCacheDF.lData[iCount++]);
						}
					}
					
					lfID = -1;
					
					if (!needToDoFixed)
						for (g = 0; g<stashedTrees.lLength; g++)
							treeVarReferences.lData[stashedTrees.lData[g]] = stashedTreeIDs.lData[g];
						
					stashedTrees.Clear();
						
					BuildLikelihoodFunction ();
					if (topConstr.sLength)
					{
						_String 	   dtc 	  (alteredTopConstr);
						_ExecutionList constr (dtc);
						constr.Execute();
						terminateExecution = false;
					}
					OptimizeLikelihoodFunction ();
					inferCache.Clear();
					inferCacheDF.Clear();
					SetLockState (false);
					return;
				}
				else
				{
					pathToFiles = "Tree inference failed.";
					ProblemReport (pathToFiles, (Ptr)this);	
				}
			}
			
			
			SetLockState (false);
			_HYTable* pl = (_HYTable*)GetObject (5);
			_SimpleList changedCells;
			
			for (long k=DF_ID_COLUMN; k<pl->cellTypes.lLength; k+=pl->horizontalSpaces.lLength)
			{
				if (pl->cellTypes.lData[k] & HY_TABLE_BOLD)
				{
					pl->cellTypes.lData[k] &= (0xffffffff-HY_TABLE_BOLD);
					changedCells << k;
					treeVarReferences.lData[k/pl->horizontalSpaces.lLength] = -3;
				}
			}
			
			pl->_MarkCellsForUpdate (changedCells);
			
			for (g = 0; g<stashedTrees.lLength; g++)
				treeVarReferences.lData[stashedTrees.lData[g]] = stashedTreeIDs.lData[g];
			
			if (lfID>=0)
				KillLFRecord (lfID);
				
			lfID = -1;
		}
	}
}

//__________________________________________________________

void	_HYDataPanel::BuildLikelihoodFunction (_String* lName, _SimpleList* subset, long treeRef)
{
	_String		errMsg;
	if (dataPartitions.lLength==0)
	{
		errMsg = _String("I can't create the likelihood function without at least one partition defined.");
		ProblemReport (errMsg,(Ptr)this);
		return;
	}
	
	long 			k=0,
					g,
					l,
					inferCount  = 0,
					treeRefID   = 0;
					
	
	_SimpleList		goodPartitions;
	
	_List			cachedTreeNames;
	
	bool			infer = GenerateGoodPartitions (goodPartitions);
	
	if (infer && (!subset))
	{
		InferTopologies ();
		return;
	}
	
	if (goodPartitions.lLength==0)
	{
		errMsg = _String ("I can't build the likelihood function because no partition has both a tree and a model attached to it.");
		ProblemReport (errMsg,(Ptr)this);
		return;
	}	
	
	if (lfID>=0)
	{
		if (lfID != lockedLFID)
		{
			errMsg = _String("There already exists a likelihood function for this data set. Would you like to rebuild the likelihood function now?");
			if (treeRef < 0)
			{
				if (ProceedPrompt (errMsg, (Ptr)this))
				{
					bool saveWLF = warnKillLF;
					warnKillLF = true;
					PurgeLF();
					warnKillLF = saveWLF;
				}
				else
					return;
			}
		}
		else
		{
			ProblemReport (lfCantKillWarning, (Ptr)this);
			return;
		}
	}
	// all is good; we can build a likelihood function
	// go through the partitions and define the models one by one
	
	_DataSetFilter* thisPartition;
	
	_List			treesList,
					otherTreeList;
					
	for (l=0; l<likeFuncNamesList.lLength; l++)
		if (((_String*)likeFuncNamesList(l))->sLength)
		{
			_SimpleList * otherTrees = &((_LikelihoodFunction*) likeFuncList (l))->GetTheTrees();
			for (k = 0; k < otherTrees->lLength; k++)
				otherTreeList << LocateVar(otherTrees->lData[k])->GetName();
		}

	for (l=0; l<dataPartitions.lLength; l++)
		if (treeVarReferences.lData[l]>=0)
			cachedTreeNames << LocateVar (treeVarReferences.lData[l])->GetName();
		else
			cachedTreeNames && & empty;
					
	_HYTable*		pl = (_HYTable*)GetObject (5);
	
	for (l=0; l<goodPartitions.lLength; l++)
	{
		k = goodPartitions.lData[l];
		
		treeRefID = treeVarReferences.lData[k];
		
		long		partIDVal;
		
		if (treeRefID >= 0)
		{
			partIDVal = dataPartitions.lData[k];
			thisPartition = (_DataSetFilter*)dataSetFilterList (partIDVal);
		}
		else
		{
			if (inferCount < inferCacheDF.lLength)
				thisPartition = (_DataSetFilter*)dataSetFilterList (inferCacheDF.lData[inferCount]);
			else
			{
				thisPartition = new _DataSetFilter ();
				checkPointer (thisPartition);
				_String 		tempFilterName ("TempFilter");
				inferCacheDF  << AddFilterToList (tempFilterName, thisPartition);
			}
			
			_DataSetFilter	 *  refFilter = (_DataSetFilter*)dataSetFilterList (dataPartitions.lData[k]);
			_SimpleList		 	vList,
								hList (*subset);
			
			thisPartition->SetFilter((_DataSet*)refFilter,refFilter->GetUnitLength(), hList, vList, true);
			
			hList.Duplicate (&refFilter->theExclusions);
			if (hList.lLength)
			{
				_String *string = GetExclusionsFromExcList (&hList);
				thisPartition->SetExclusions (string);
				DeleteObject (string);	
			}	
			
			partIDVal = inferCacheDF.lData[inferCount];
		}
		
		int model, 
			options, 
			freqs, 
			classes;
			
		LongToModelData (modelReferences.lData[k],model,options,freqs,classes);
		_List* 	   theModel = FindModelTemplate (model-1,thisPartition->GetDimension());
		_String*   modelDef = nil;
		if (theModel)
		{
			FILE * thisModel = doFileOpen (((_String*)(*theModel)(2))->getStr(),"rb");
			if (thisModel)
			{
				modelDef = new _String (thisModel);
				fclose (thisModel);
				_String ident, ident2, ident3;
				*modelDef = _String ("modelType=")&_String ((long)options)&';'&*modelDef;
				ident = *(_String*)dataSetFilterNamesList (partIDVal) & "_Shared";
				if (!infer)
					FindUnusedObjectName (errMsg, ident, variableNames,true);
					
				*modelDef = modelDef->Replace (globalPrefix,ident,true);
				if (options)
				{
					if (options>1)
					{
						ident = *(_String*)dataSetFilterNamesList (partIDVal) & "_Categ";
						if (!infer)
							FindUnusedObjectName (errMsg, ident, variableNames,true);
						*modelDef = modelDef->Replace (categoryPrefix,ident,true);	
										
						ident = *(_String*)dataSetFilterNamesList (partIDVal) & "_Shape";
						if (!infer)
							FindUnusedObjectName (errMsg, ident, variableNames,true);
						*modelDef = modelDef->Replace ("shapeParameter",ident,true);
						
						ident = _String ((long)(classes));
						*modelDef = modelDef->Replace ("rateClassCount",ident,true);							
					}
				}
				ident = *(_String*)dataSetFilterNamesList (partIDVal) & "_Freqs";
				
				if (!infer)
					FindUnusedObjectName (errMsg, ident, variableNames,true);
					
				switch (freqs)
				{
					case 0: // partition
						*modelDef = *modelDef & "\nHarvestFrequencies("&ident&','&
									*(_String*)dataSetFilterNamesList (partIDVal)&","& (long)thisPartition->GetUnitLength() &",1,1);";
						break;
					case 1: // dataset
						*modelDef = *modelDef & "\nHarvestFrequencies("&ident&','&
									*(_String*)dataSetNamesList (dataSetID)&","& (long)thisPartition->GetUnitLength() &",1,1);";
						break;		
					case 2: // equal
						{
							long					    colDim =	thisPartition->GetUnitLength(),
														rowDim =	(colDim==1) ? 
																	thisPartition->GetDimension():
																	thisPartition->GetData()->GetCharDimension();
														
							_Matrix						eqMatrix (rowDim,colDim,false,true);
							_Parameter					fe		  = 1./rowDim;
							
							for (;rowDim>=0;rowDim--)
								for (long c=0; c<colDim; c++)
									eqMatrix.Store(rowDim,c,fe);
																					
							*modelDef = *modelDef & '\n' & ident & '=' & _String((_String*)eqMatrix.toStr()) & ';';
						}
						break;
					case 3: // estimate
						{
							ident2 = ident & "_emp";
							
							if (!infer)
								FindUnusedObjectName (errMsg,ident2,variableNames,true);
							
							ident3 = *(_String*)dataSetFilterNamesList (partIDVal) & "_EstFreq";
							
							if (!infer)	
								FindUnusedObjectName (errMsg,ident3,variableNames,true);
							
							*modelDef = modelDef->Replace ("frequencyVariable",ident3,true);							
							*modelDef = *modelDef & "\nHarvestFrequencies("&ident2&','&
										*(_String*)dataSetFilterNamesList (partIDVal)&",1,1,1);";	
										
							*modelDef = *modelDef & '\n' & ident &" = 0;\ndummy = " & efvFunction & "(\"" & 
							ident & "\"," & ident2 & ");";
						}
						break;
					case 4: // model spec
						{
							*modelDef = modelDef->Replace (dataPartitionIDString,*(_String*)dataSetFilterNamesList (partIDVal),true);
							*modelDef = modelDef->Replace (dataSetIDString,*(_String*)dataSetNamesList (dataSetID),true);
							*modelDef = modelDef->Replace (modelEFVVector,ident,true);				
							*modelDef = *modelDef & "dummy=" & userModelEFV & "(0);\n";			
						}
						break;
				}
				// if codon, then create genetic code vector
				if (thisPartition->GetUnitLength()==3)
				{
					_String * genT = GetMatrixFromCode ((partData.lData[k]&HY_DATAPANEL_CODEMASK)>>16);
					*modelDef = *modelDef & '\n'  &modelGenCode & '=' & *genT &';' &  modelMatrixDimension & "=0;\n";
					DeleteObject (genT);				
				}
				
				ident2 = *(_String*)dataSetFilterNamesList (partIDVal) & '_'& *((_String*)(*theModel)(0));
				if (!infer)
					FindUnusedObjectName (ident, ident2, variableNames,true);
					
				ident3 = ident2 & "_model";
				if (!infer)
					FindUnusedObjectName (ident, ident3, modelNames);
					
				if (freqs<4)
					*modelDef = *modelDef & '\n' & ident2 &" = 0;\n"&multiplyByFrequencies&'=' & modelFunction & "(\"" & ident2 & "\"," & ident & ");";
				else
					*modelDef = *modelDef & '\n' & ident2 &" = 0;\ndummy=" & modelFunction & "(\"" & ident2 & "\"," & ident & ");";
					
				if (thisPartition->GetUnitLength()>1)
					*modelDef = *modelDef & '\n' & ident & '=' & buildCodonFrequencies & '(' & ident & ");";
				*modelDef = *modelDef & "\nModel "& ident3 & "=("& ident2 &','& ident & ','& multiplyByFrequencies&");";
				
				/*FILE * flfl = doFileOpen ("debug.dump","w");
				fprintf (flfl,"%s\n", modelDef->getStr());
				fclose (flfl);*/
								
				_Variable *thisTree;
				
				_String	  *treeString,
						  treeName;
						  
				if (treeRefID>=0)
				{
					thisTree = FetchVar(LocateVarByName (*(_String*)cachedTreeNames (k)));
					
					treeString = (_String*) thisTree->toStr();
					treeName   = *thisTree->GetName();
					
					if (!(infer && inferCache.lLength))
						if ((treesList.Find (&treeName)>=0)||(otherTreeList.Find (&treeName)>=0))
						{	
							_String  os (128L, true);
							os << "\nTree topology";
							os << treeName;
							os << " was cloned for partition ";
							os << ((_String*)dataSetFilterNamesList (dataPartitions.lData[k]));
							os.Finalize ();
							StringToConsole (os);
							
							FindUnusedObjectName (treeName,treeName,variableNames,true);
							pl->SetCellData (&treeName,k,DF_TREE_COLUMN,pl->cellTypes.lData[k*pl->horizontalSpaces.lLength+DF_TREE_COLUMN],true);
							pl->_MarkCellForUpdate (k*pl->horizontalSpaces.lLength+DF_TREE_COLUMN);
						}
					
				}
				else
				{
					thisTree = LocateVar (treeRef);
					
					treeString = (_String*) thisTree->toStr(),
					treeName   = *(_String*)dataSetNamesList (dataSetID) & '_' & *thisTree->GetName();
					
					if (inferCount < inferCache.lLength)
						treeName = *(_String*)inferCache (inferCount);
					else
					{
						if ((treesList.Find (&treeName)>=0)||(otherTreeList.Find (&treeName)>=0))
							FindUnusedObjectName (treeName,treeName,variableNames,true);
						
						inferCache && & treeName;
					}
					inferCount ++;
				}
				
				treesList && &treeName;
				*modelDef = *modelDef & "\nTree " & treeName & '=' & *treeString;
				DeleteObject (treeString);
				
				g = batchLanguageFunctionNames.lLength;

				_ExecutionList		ex;
				ex.BuildList (*modelDef);
				ex.ExecuteAndClean(g);
				if (terminateExecution)
				{
					terminateExecution = false;
					DeleteObject (modelDef);
					return;
				}
				
				if (treeRefID>=0)
					#ifndef USE_AVL_NAMES
						treeVarReferences.lData[k]=variableReindex.lData[LocateVarByName(treeName)];
					#else
						treeVarReferences.lData[k]=variableNames.GetXtra(LocateVarByName(treeName));					
					#endif
			}
		}
		// finally write the likelihood function
		if (!(theModel&&modelDef))
		{
			errMsg = _String ("I can't build the likelihood function because partition '")&
					 *(_String*)dataSetFilterNamesList(partIDVal)& _String("'s model definition file couldn't be found/read. ")&
					 ReportThisError;
			ProblemReport (errMsg,(Ptr)this);
			return;
		}
		thisPartition->SetDimensions();
		thisPartition->SetupConversion();
		DeleteObject (modelDef);
	}
	
	if ( lName )
		errMsg = *lName;
	else
		errMsg = *(_String*)dataSetNamesList (dataSetID)&"_LF";
		
	if (!infer)
		FindUnusedObjectName (errMsg, errMsg, likeFuncNamesList);	
		
	_String	  lfSetup ((unsigned long)32,true);
	lfSetup << "\nLikelihoodFunction ";
	lfSetup << &errMsg;
	lfSetup << '=';
	lfSetup << '(';
	
	inferCount = 0;
	
	for (l=0; l<goodPartitions.lLength; l++)
	{
		k = goodPartitions.lData[l];
		if (l)
			lfSetup<<',';
			
		long 		treeRefID = treeVarReferences.lData[k];
		
		if (treeRefID >= 0)
			treeRefID = dataPartitions.lData[k];
		else
			treeRefID = inferCacheDF.lData[inferCount++];

		lfSetup << (_String*)dataSetFilterNamesList(treeRefID);
		lfSetup << ',';
		lfSetup << (_String*)treesList(l);
		k *= pl->horizontalSpaces.lLength;
		k ++;
		if (!(pl->cellTypes.lData[k]&HY_TABLE_BOLD))
		{
			pl->cellTypes.lData[k]|=HY_TABLE_BOLD;
			pl->_MarkCellForUpdate (k);
		}
	}
	lfSetup << ')';
	lfSetup << ';';
	lfSetup.Finalize();
	_ExecutionList		ex;
	ex.BuildList (lfSetup);
	ex.Execute();
	if (terminateExecution)
	{
		if (!subset)
			terminateExecution = false;
		return;
	}
	
	_PaintThermRect();
	
	lfID = likeFuncNamesList.Find (&errMsg);
	
	if (!infer)
	{
		_LikelihoodFunction *lf = (_LikelihoodFunction*)likeFuncList (lfID);
		
		char	buffer [1024];
	
		sprintf (buffer,"\nCreated likelihood function '%s' with\n %d\tpartitions,\n %d\tshared parameters,\n %d\tlocal parameters,\n %d\tconstrained parameters.\n",
				  errMsg.getStr(), lf->CountObjects (0),  lf->CountObjects (1),  lf->CountObjects (2),  lf->CountObjects (3)); 
				  
		BufferToConsole (buffer);
	
		lf->ComputePruningEfficiency (k,g);
		sprintf (buffer,"\nPruning efficiency %d vs %d (%g %% savings)\n", g,k, 100.-g*100./k);
		BufferToConsole (buffer);
	}
	else
	{
		if (topConstr.sLength)
		{
			_String 	   dtc 	  (topConstr);
			_ExecutionList constr (dtc);
			constr.Execute();
			if (terminateExecution)
			{
				errMsg = "Problems with topology constraints; tree inference failed.";
				ProblemReport (errMsg);
			}
		}
	}
	
	long verbLebel = VerbosityLevel ();
	ApplyPreferences();
	setParameter (VerbosityLevelString,verbLebel);
	SetStatusBarValue (-1,1,0);
	SetStatusLine ("Idle");
	
	_HYButtonBar*   	bb = (_HYButtonBar*)GetObject (2);
	bb->EnableButton (8,true);
	
	_UpdateLFMenu();
	ReportAnalysisAsFinished("Likelihood function build finished");
	postLFSpawnEvent (GetID(),lfID);
}

//__________________________________________________________

long	_HYDataPanel::SpawnLikelihoodFunction (_DataSet* ds, _String* dsName, _List& partitionCache, _SimpleList& sequenceCache, _SimpleList* remapper, _SimpleList* seqMap)
{
	if (lfID>=0)
	{
		_String			errMsg;
		_DataSetFilter* thisPartition;
		long			k,g;
		int				model,options,freqs,classes;
		
		bool			restore = (partitionCache.lLength>0);

		
		_SimpleList		goodPartitions;
		GenerateGoodPartitions (goodPartitions);
		
		long			globalOffset = 0;
		
		for (long j=0; j<goodPartitions.lLength; j++)
		{
			k = goodPartitions.lData[j];
			thisPartition = (_DataSetFilter*)dataSetFilterList (dataPartitions.lData[k]);
			_SimpleList  vSpec,
						 hSpec,
						 excl (thisPartition->theExclusions);
						 
			if (restore)
			{
				hSpec.Duplicate (&sequenceCache);
				vSpec.Duplicate (partitionCache(j));
			}
			else
			{
				partitionCache && & thisPartition->theOriginalOrder;
				if (sequenceCache.lLength==0)
					sequenceCache.Duplicate (&thisPartition->theNodeMap);
					
				if (seqMap)
					hSpec.Duplicate (seqMap);
				else
				{	
					for (g=0; g<thisPartition->theNodeMap.lLength; g++)
						hSpec << g;	
				}

				if (remapper)
				{
					/* this will break if the null and the alternative
					   span different ranges of sites
					*/
					for (g = 0; g < thisPartition->theOriginalOrder.lLength; g++)
						vSpec << remapper->lData[thisPartition->theOriginalOrder[g]];				
				}
				else
					for (g=0; g<thisPartition->theOriginalOrder.lLength; g++)
						vSpec << globalOffset+g;
					
							
				globalOffset += thisPartition->theOriginalOrder.lLength;
			}
			
			
			thisPartition->SetFilter (ds,thisPartition->GetUnitLength(),hSpec,vSpec,false);
						
			if (excl.lLength)
			{
				_String *string = GetExclusionsFromExcList (&excl);
				thisPartition->SetExclusions (string, false);
				DeleteObject (string);	
			}	
			
			LongToModelData (modelReferences.lData[k],model,options,freqs,classes);
			
			if (freqs<2)
			{
				_List* 	   theModel = FindModelTemplate (model-1,thisPartition->GetDimension());
				_String*   modelDef = nil;
				
				if (theModel)
				{
					FILE * thisModel = doFileOpen (((_String*)(*theModel)(2))->getStr(),"rb");
					if (thisModel)
					{
						modelDef = new _String (thisModel);
						fclose (thisModel);
						
						_SimpleList modelID;
						_TheTree*   thisTree = (_TheTree*)LocateVar (treeVarReferences.lData[k]);
						thisTree->CompileListOfModels (modelID);
						
						g = modelFrequenciesIndices.lData[modelID.lData[0]];
						if (g<0)
							g = -g-1;
						
						_String ident = *LocateVar(g)->GetName();
						
						switch (freqs)
						{
							case 0: // partition
								*modelDef = *modelDef & "\nHarvestFrequencies("&ident&','&
											*(_String*)dataSetFilterNamesList (dataPartitions.lData[k])&",1,1,1);";
								break;
							case 1: // dataset
								*modelDef = *modelDef & "\nHarvestFrequencies("&ident&','&
											*dsName&",1,1,1);";
								break;		
						}
						// if codon, then create genetic code vector
						if (thisPartition->GetUnitLength()==3)
						{
							_String * genT = GetMatrixFromCode ((partData.lData[k]&HY_DATAPANEL_CODEMASK)>>16);
							*modelDef = *modelDef & '\n'  &modelGenCode & '=' & *genT &';';
							DeleteObject (genT);
						
						}
						if (thisPartition->GetUnitLength()==3)
							*modelDef = *modelDef & '\n' & ident & '=' & buildCodonFrequencies & '(' & ident & ");";
							
						g = batchLanguageFunctionNames.lLength;

						_ExecutionList		ex;
						ex.BuildList (*modelDef);
						ex.ExecuteAndClean(g);
						if (terminateExecution)
						{
							terminateExecution = false;
							DeleteObject (modelDef);
							return -1;
						}
					}
				}
				DeleteObject (modelDef);
			}
			thisPartition->SetDimensions();
			thisPartition->SetupConversion();
		}
		
		_LikelihoodFunction * thisLF = (_LikelihoodFunction*)likeFuncList (lfID);
		thisLF->Rebuild();
		long verbLebel = VerbosityLevel ();
		ApplyPreferences();
		setParameter (VerbosityLevelString,verbLebel);
		SetStatusBarValue (-1,1,0);
		SetStatusLine ("Idle");
		return	lfID;  
	}
	return -1;
}

//__________________________________________________________

void	_HYDataPanel::RefreshCategoryVars (void)
{
	if (lfID>=0)
	{
		_SimpleList * cVars = &((_LikelihoodFunction*)likeFuncList (lfID))->GetCategoryVars();
		for (long k=0; k<cVars->lLength; k++)
			((_CategoryVariable*)LocateVar(cVars->lData[k]))->Refresh(true);
	}
}


//__________________________________________________________

long	_HYDataPanel::SpawnLikelihoodFunctionNP (_List& cachedPartitions, bool permute)
{
	if (lfID>=0)
	{
		_String			errMsg;
		_DataSetFilter* thisPartition;
		_DataSet*		ds = (_DataSet*)dataSetList (dataSetID);
		long			k,g;
		int				model,options,freqs,classes;
		bool			restore = (cachedPartitions.lLength > 0);

		
		_SimpleList		goodPartitions;
		GenerateGoodPartitions (goodPartitions);
		
		for (long j=0; j<goodPartitions.lLength; j++)
		{
			k = goodPartitions.lData[j];
			thisPartition = (_DataSetFilter*)dataSetFilterList (dataPartitions.lData[k]);
			_SimpleList  vSpec 		,
						 hSpec 		(thisPartition->theNodeMap),
						 excl		(thisPartition->theExclusions);
						 
			if (restore)
				vSpec.Duplicate (cachedPartitions(j));
			else
			{
				vSpec.Duplicate (&thisPartition->theOriginalOrder);
				if (permute)
					vSpec.Permute (thisPartition->GetUnitLength());
				else
					vSpec.PermuteWithReplacement (thisPartition->GetUnitLength());
				cachedPartitions && &thisPartition->theOriginalOrder;
			}
			
			thisPartition->SetFilter (ds,thisPartition->GetUnitLength(),hSpec,vSpec,false);
						
			if (excl.lLength)
			{
				_String *string = GetExclusionsFromExcList (&excl);
				thisPartition->SetExclusions (string, false);
				DeleteObject (string);	
			}	
			
			LongToModelData (modelReferences.lData[k],model,options,freqs,classes);
			
			if (freqs==0)
			{
				_List* 	   theModel = FindModelTemplate (model-1,thisPartition->GetDimension());
				_String*   modelDef = nil;
				
				if (theModel)
				{
					FILE * thisModel = doFileOpen (((_String*)(*theModel)(2))->getStr(),"rb");
					if (thisModel)
					{
						modelDef = new _String (thisModel);
						fclose (thisModel);
						
						_SimpleList modelID;
						_TheTree*   thisTree = (_TheTree*)LocateVar (treeVarReferences.lData[k]);
						thisTree->CompileListOfModels (modelID);
	
						g = modelFrequenciesIndices.lData[modelID.lData[0]];
						if (g<0)
							g = -g-1;

						_String ident = *LocateVar(g)->GetName();
						
						*modelDef = *modelDef & "\nHarvestFrequencies("&ident&','&
									*(_String*)dataSetFilterNamesList (dataPartitions.lData[k])&",1,1,1);";

						// if codon, then create genetic code vector
						if (thisPartition->GetUnitLength()==3)
						{
							_String * genT = GetMatrixFromCode ((partData.lData[k]&HY_DATAPANEL_CODEMASK)>>16);
							*modelDef = *modelDef & '\n'  &modelGenCode & '=' & *genT &';';
							DeleteObject (genT);
						
						}
						if (thisPartition->GetUnitLength()==3)
							*modelDef = *modelDef & '\n' & ident & '=' & buildCodonFrequencies & '(' & ident & ");";
							
						g = batchLanguageFunctionNames.lLength;

						_ExecutionList		ex;
						ex.BuildList (*modelDef);
						ex.Execute();
						if (terminateExecution)
						{
							terminateExecution = false;
							DeleteObject (modelDef);
							return -1;
						}
						
						while (g<batchLanguageFunctionNames.lLength)
						{
							batchLanguageFunctionNames.Delete (g);
							batchLanguageFunctionParameters.Delete (g);
							batchLanguageFunctions.Delete(g);
							batchLanguageFunctionParameterLists.Delete(g);
							batchLanguageFunctionClassification.Delete(g);
						}

					}
				}
				DeleteObject (modelDef);
			}
			thisPartition->SetDimensions();
			thisPartition->SetupConversion();
		}
		
		_LikelihoodFunction * thisLF = (_LikelihoodFunction*)likeFuncList (lfID);
		thisLF->Rebuild();
		long verbLebel = VerbosityLevel ();
		ApplyPreferences();
		setParameter (VerbosityLevelString,verbLebel);
		SetStatusBarValue (-1,1,0);
		SetStatusLine ("Idle");
		return	lfID;  
	}
	return -1;
}
//__________________________________________________________
void	_HYDataPanel::ConstructDataTypeOptions (_List& options)
{
	if (dataType&HY_DATAPANEL_NUCDATA)
	{
		options && & nucDataType;
		options && & dinucDataType;
		options && & codonDataType;
	}
	else
		if (dataType&HY_DATAPANEL_PROTDATA)
			options && & proteinDataType;
		else
			if (dataType&HY_DATAPANEL_BINARYDATA)
			{
				options && & binaryDataType;
				options && & disequonDataType;
			}
			else
				options && & unknownDataType;
	
}	
//__________________________________________________________
void	_HYDataPanel::BuildDataPartitions (void)
{
	if (dataSetID>=0)
	{
		_DataSet *theDS = (_DataSet*)dataSetList (dataSetID);
		_HYSequencePane* sp = (_HYSequencePane*)GetObject (0);
				
		for (long k=0; k<dataSetFilterList.lLength; k++)
		{
			_DataSetFilter * thisDF = (_DataSetFilter*) dataSetFilterList (k);
			if (thisDF && thisDF->GetData()==theDS)
			{
				if (thisDF->theNodeMap.Equal (sp->speciesIndex))
					AddPartition(thisDF);
			}
		}
		/* add missing species to the 'omitted list'*/
		if (sp->speciesIndex.lLength<theDS->NoOfSpecies())
		{
			_SimpleList sorted;
			sorted.Duplicate (&sp->speciesIndex);
			sorted.Sort();
			long   shift = 0;
			for (long k=0; k<sorted.lLength; k++,shift++)
			{
				while (sorted.lData[k]>shift)
				{
					omittedSeqs<<shift;
					shift++;
				}
			}
		}
		GenerateStatusLine();
	}
}

//__________________________________________________________
bool	_HYDataPanel::ConfirmClose (void)
{
	if ((lfID>=0)&&(lfID == lockedLFID))
	{
		ProblemReport (lfCantKillWarning, (Ptr)this);
		return false;
	}
	
	if ((tainted)||((filePath.sLength==0)&&(savePath.sLength==0)))
	{
		_String warnMessage ("Would you like to save partitioning info before '");
		warnMessage = warnMessage & GetTitle() & "' is closed?";
		char    res = YesNoCancelPrompt(warnMessage);
		if (res!=2)
		{
			if (res==1)
			// generate the save
			{
				if (!SaveDataPanel (savePath.sLength))
					return false;
			}
		}
		else
			return false;
	}
	
	// check to see if standard analyses made use of this data set
	
	for (long k2 = likeFuncList.lLength-1; k2 >= 0; k2--)
	{
		if (k2 != lfID && likeFuncList(k2))
		{
			if (((_LikelihoodFunction*)likeFuncList(k2))-> DependOnDS (dataSetID)>=0)
			{
				postLFKillEvent (GetID(),k2);
				KillLFRecord (k2);
			}	
		}
	}
	
	if (lfID>=0)
	{
		postLFKillEvent (GetID(),lfID);
		KillLFRecord (lfID);
		lfID = -1;
	}
	

	for (long k=0; k<treeVarReferences.lLength; k++)
		if (treeVarReferences.lData[k]>=0)
		{
			_TheTree * thisTree = (_TheTree*)LocateVar(treeVarReferences.lData[k]);
			if (thisTree)
			{
				postTreeKillEvent (GetID(),treeVarReferences.lData[k]);
				handleGUI();
				DeleteVariable (*thisTree->GetName());
			}
		}

	for (long k=0; k<dataPartitions.lLength; k++)
		KillDataFilterRecord (dataPartitions.lData[k]);
		
	KillDataSetRecord (dataSetID);

	return true;
}	

//__________________________________________________________
_DataSet*	_HYDataPanel::GenerateOrderedDataSet (void)
{
	_DataSet* newDS = new _DataSet();
	checkPointer (newDS);
	
	if (lfID>=0)
	{
		_SimpleList    *lfFilters = &((_LikelihoodFunction*)likeFuncList(lfID))->GetTheFilters ();
		_DataSetFilter *dsf = (_DataSetFilter*)dataSetFilterList (lfFilters->lData[0]);
		_DataSet	   *cDS = dsf->GetData();
		
		long 		   i,
					   k,
					   sc;

		for (i=0; i<dsf->NumberSpecies(); i++)
			newDS->GetNames() << dsf->GetData()->GetNames() (((_SimpleList*)dsf->GetMap())->lData[i]);
			
		sc = newDS->GetNames().lLength;
			
		for (i=0; i<lfFilters->lLength; i++)
		{
			dsf = (_DataSetFilter*)dataSetFilterList (lfFilters->lData[i]);
			for (k=0; k<dsf->theOriginalOrder.lLength; k++)
			{
				_Site * thisS = cDS->GetSite(dsf->theOriginalOrder.lData[k]);
				newDS->AddSite(thisS->sData[0]);
			}
			
			for (long j=1; j<sc; j++)
			{
				_Site * thisS = cDS->GetSite(dsf->theOriginalOrder.lData[k]);
				for (k=0; k<dsf->theOriginalOrder.lLength; k++)
					newDS->Write2Site (k, thisS->sData[j]);
			}
		}

		newDS->SetTranslationTable (dsf->GetData());
	}
	
	newDS->Finalize();
	newDS->SetNoSpecies(newDS->GetNames().lLength);
	return newDS;	

}

//__________________________________________________________
bool	_HYDataPanel::SaveDataPanel (bool saveAs, _String* saveFile, _String* dsPath, bool saveAllStates, _DataSet* dsXtra, bool straightOrder)
{
	FILE* dsout = nil;

	long	includeData = 0;
	
	if ((filePath.sLength==0)||dsXtra)
	{
		if (dsPath||savePath.sLength||ProceedPrompt (saveDSPrompt,(Ptr)this))
		{
			if (dsPath)
			{
				if (dsXtra)
				{
					*dsPath = SavePartition (-1,dsPath, dsXtra);
					if (dsPath->sLength==0)
						return false;
				}
			}
			else
			{
				if (savePath.sLength)
					includeData = 1;
				else
				{
					filePath = SavePartition (-1);
					if (filePath.sLength==0)
						return false;
				}
			}
		}
		else
			return false;
	}
	
	 
	if (((savePath.sLength == 0)||saveAs)&&(!saveFile))
	// prompt for destination
	{
		_String ffPrompt ("Save partitioning info to:"),
				foPrompt ("Format:"),
				ffDefName(*((_String*)dataSetNamesList(dataSetID))),
				ffOption1("Do not include sequence data"),
				ffOption2("Include sequence data, NEXUS"),
				ffOption3("Export data, partitions and trees to NEXUS");
				
		_List	ffOptions;
				
		ffOptions && & ffOption1;
		ffOptions && & ffOption2;
		ffOptions && & ffOption3;
		
		
		includeData = SaveFileWithPopUp (savePath,ffPrompt,ffDefName,foPrompt,ffOptions);
		
		if (includeData < 0)
			return false;
			
		terminateExecution = false;
		if (savePath.sLength == 0) 
			return false;
	}
	
	_DataSet* thisDS = (_DataSet*)dataSetList (dataSetID);
	
	_String	  relPath;
	
	if (saveFile)
	{
		if (dsPath)
		{
			relPath = saveFile->PathSubtraction(*dsPath,1);
			if (relPath.sLength == 0) 
				relPath = *dsPath;	
		}
		else
		{
	 		relPath = savePath.PathSubtraction(filePath,1);
	 		if (relPath.sLength == 0)
	 			relPath = filePath;
	 	}
	}
	else
	{
	 	relPath = savePath.PathSubtraction(filePath,1);
		if (relPath.sLength == 0) 
			relPath = filePath;
	}
	
	_HYSequencePane* sp = (_HYSequencePane*)GetObject (0);
	
	dsout = doFileOpen ((saveFile?saveFile->sData:savePath.sData),"w");
	
	if (!dsout)
	{
		_String errMsg = _String("Trying to open '")&(saveFile?*saveFile:savePath)&"' for writing failed. Oops!";
		ProblemReport (errMsg,(Ptr)this);
		return false;
	}
	
	if (includeData == 2)
	{
		stashParameter (dataFilePrintFormat,4,true);
		_SimpleList		dummyH, dummyV;
		_DataSetFilter  dummyDF;
		
		stashParameter (dataFileTree,0,true);
		dummyDF.SetFilter ((_DataSet*)dataSetList(dataSetID),1,dummyH,dummyV);
		dummyDF.toFileStr (dsout);
		stashParameter (dataFileTree,0,false);
		
		if (dataPartitions.lLength)
		{
			fprintf (dsout, "\n\nBEGIN ASSUMPTIONS;");
			for (long k=0; k<dataPartitions.lLength; k++)
			{
				_DataSetFilter* theDF = (_DataSetFilter*)dataSetFilterList (dataPartitions.lData[k]);
				fprintf (dsout, "\n\tCHARSET %s = ", ((_String*)dataSetFilterNamesList (dataPartitions.lData[k]))->sData);
				
				_SimpleList offsetList;
				for (long k2 = 0; k2 < theDF->theOriginalOrder.lLength; k2=k2+1)
					offsetList << (theDF->theOriginalOrder.lData[k2]+1);
			
				_String partString ((_String*)offsetList.ListToPartitionString());
				fprintf (dsout, "%s;", partString.sData);
			}
			fprintf (dsout, "\nEND;\n\n");
			
			fprintf (dsout, "\n\nBEGIN TREES;\n\n");
			
			_SimpleList     alreadyDone;
			_AVLList		alreadyDoneAVL (&alreadyDone);
			
			for (long k=0; k<treeVarReferences.lLength; k++)
			{
				long k2 = treeVarReferences.lData[k];
				if (k2>=0 && alreadyDoneAVL.Find ((BaseRef)k2) < 0)
				{
					_TheTree* myTree = (_TheTree*)LocateVar (k2);
					myTree->StepWiseT(true);
					_String  *tStr = new _String  ((unsigned long)1024,true);
					myTree->SubTreeString (*tStr, false, -1, nil);
					tStr->Finalize();
					fprintf (dsout, "\n\tTREE %s = %s;", myTree->GetName()->sData, tStr->sData);
					DeleteObject (tStr);
					alreadyDoneAVL.Insert ((BaseRef)k2);
				}
			}
			fprintf (dsout, "\nEND;\n\n");
		}
	}
	else
	{
		if (includeData)
		{
			stashParameter (dataFilePrintFormat,4,true);
			_SimpleList		dummyH, dummyV;
			_DataSetFilter  dummyDF;
			
			dummyDF.SetFilter ((_DataSet*)dataSetList(dataSetID),1,dummyH,dummyV);
			dummyDF.toFileStr (dsout);
			
			fprintf (dsout, "\n\nBEGIN HYPHY;\n\n");
			
			stashParameter (dataFilePrintFormat,0,false);
			
			relPath = empty;
		}

		fprintf (dsout, "%s=\"%s\";\nDataSet %s = ReadDataFile (%s);\n",
						 dataPanelSourcePath.getStr(),
						 relPath.getStr(), 
						 ((_String*)dataSetNamesList(dataSetID))->getStr(),
						 (includeData?useNexusFileData.getStr():dataPanelSourcePath.getStr()));

		_String  *vertPart;
		
		if (sp->speciesIndex.lLength == thisDS->NoOfSpecies())
			vertPart = new _String(empty);
		else
			vertPart = (_String*)sp->speciesIndex.ListToPartitionString();
		
		_String			partInfo (16,true);
		_List	    	treeList;
		_SimpleList*	lfDFs = lfID>=0?(&((_LikelihoodFunction*)likeFuncList(lfID))->GetTheFilters()):nil;
		long			glOffset = 0;
		
		for (long k=0; k<dataPartitions.lLength; k++)
		{
			if (lfDFs && (!saveAllStates))
				if (lfDFs->Find (dataPartitions.lData[k])<0)
					continue;
			
			_DataSetFilter* theDF = (_DataSetFilter*)dataSetFilterList (dataPartitions.lData[k]);
			_String		  * horPart;
			
			if (straightOrder)
			{
				horPart = new _String(empty);
				checkPointer (horPart);
				*horPart = _String (glOffset) & '-' & _String((long)(glOffset+theDF->theOriginalOrder.lLength-1));
			}
			else
			 	horPart = (_String*)theDF->theOriginalOrder.ListToPartitionString();
			 	
			_String			exclusions (empty);
			if ((theDF->GetUnitLength()==3)&&theDF->theExclusions.lLength)
			{
				DFExclusionsToString	  (theDF,exclusions);
				exclusions = _String(",\"")&exclusions&'"';
			}

			fprintf (dsout, "DataSetFilter %s = CreateFilter (%s,%d,\"%s\",\"%s\",\"%s\");\n", 
					 ((_String*)dataSetFilterNamesList(dataPartitions.lData[k]))->getStr(),							
					 ((_String*)dataSetNamesList(dataSetID))->getStr(),
					 theDF->GetUnitLength(),
					 horPart->getStr(),
					 vertPart->getStr(),
					 exclusions.getStr());
					 
			DeleteObject (horPart);
			//partInfo << _String (modelReferences.lData[k]).getStr();
			int model, options, freqs, classes;
			LongToModelData (modelReferences.lData[k],model,options,freqs,classes);
			_List*   mList = FindModelTemplate (model-1,theDF->GetDimension());
			if (mList)
			{
				partInfo << (_String*)(*mList)(0);
				partInfo << ',';
				partInfo << _String (modelReferences.lData[k]).getStr();
			}
			else
			{
				partInfo << _String ((long)-1);
				partInfo << ',';
				partInfo << _String ((long)0);
			}
			partInfo << ',';
			partInfo << _String (partData.lData[k]).getStr();
			partInfo << ',';
			partInfo << _String (partitionColors.lData[k]).getStr();
			partInfo << ',';

			if (treeVarReferences.lData[k]>=0)
			{
				_String * tName = LocateVar (treeVarReferences.lData[k])->GetName();
				partInfo << tName;
				if (treeList.Find (tName) < 0)
					treeList << tName;
			}
			else
				partInfo << "No_tree";
				
			if (k<dataPartitions.lLength-1)
				partInfo << ';';
				
			glOffset += theDF->theOriginalOrder.lLength;
		}
		
		partInfo.Finalize();	
		for (long m=0; m<treeList.lLength; m++)
		{
			_TheTree* tVar = (_TheTree*)FetchVar(LocateVarByName (*(_String*)treeList(m)));
			_String* tStr = (_String*)tVar->toStr();
			fprintf (dsout, "Tree %s=%s\n", ((_String*)treeList(m))->getStr(), tStr->getStr());
			DeleteObject (tStr);
			/*tStr = tVar->TreeUserParams();
			if (tStr->sLength)
				fprintf (dsout, "%s\n", tStr->getStr());
			DeleteObject (tStr);*/
		}
		_String		panelState (16,true);
		panelState << _String((long)addedLines).getStr();
		panelState << ',';
		panelState << _String(referenceSequence).getStr();
		panelState << ',';
		panelState << _String(translatedSequence).getStr();
		panelState << ',';
		panelState << _String(sp->startColumn).getStr();
		panelState << ',';
		panelState << _String((long)sp->nameDisplayFlags).getStr();
		panelState << ',';
		panelState << _String((long)sp->showDots).getStr();
		panelState << ',';
		panelState << _String((long)sp->blockWidth).getStr();
		panelState.Finalize();
		//panelState << ',';
		
		if (savedLFNames.lLength && saveAllStates)
		{
			fprintf (dsout, "\n%s = {", savedLFMatrix.getStr());
			
			for (long idx = 0; idx<savedLFNames.lLength; idx++)
				fprintf (dsout, "\n{\"%s\",\n\"%s\"}", ((_String*)savedLFNames(idx))->getStr(),	
													   ((_String*)savedLFStates(idx))->getStr());
			fprintf (dsout, "\n};\n\n");
		}
		
		fprintf (dsout, "OpenDataPanel(%s,\"%s\",\"%s\",\"%s\"", 
						((_String*)dataSetNamesList(dataSetID))->getStr(),
						vertPart->getStr(),panelState.getStr(),partInfo.getStr());
		if (lfID>=0)
		{
			fprintf (dsout,",%s);\n",((_String*)likeFuncNamesList(lfID))->getStr());
			stashParameter (likefuncOutput,4.0,true);
			((_LikelihoodFunction *)likeFuncList (lfID))->toFileStr(dsout);
			stashParameter (likefuncOutput,0.0,false);
		}
		else
			fprintf (dsout,");\n");
		if (includeData)
		{
			fprintf (dsout, "\n\nEND;\n");
		}
		DeleteObject (vertPart);
	}
	fclose (dsout);
	tainted = false;
	return  true;
}
//__________________________________________________________
void	_HYDataPanel::RestorePartInfo (_String* pInfo)
{
	_List parts;
	pInfo->StripQuotes();
	_ElementaryCommand::ExtractConditions (*pInfo,0,parts,';');
	for (long k=0; (k<parts.lLength)&&(k<dataPartitions.lLength); k++)
	{
		_List	 thisPart;
		_String* thisInfo = (_String*)parts(k),
				 errMsg;
		_ElementaryCommand::ExtractConditions (*thisInfo,0,thisPart,',');
		if (thisPart.lLength == 5)
		{
			#ifndef USE_AVL_NAMES
				treeVarReferences.lData[k]= variableNames.BinaryFind((_String*)thisPart(4));
				if (treeVarReferences.lData[k]>=0)
					treeVarReferences.lData[k] = variableReindex.lData[treeVarReferences.lData[k]];
				else
					treeVarReferences.lData[k]=-1;
			#else
				treeVarReferences.lData[k]=  LocateVarByName (*(_String*)thisPart(4));
				if (treeVarReferences.lData[k]>=0)
					treeVarReferences.lData[k] = variableNames.GetXtra(treeVarReferences.lData[k]);
				else
					treeVarReferences.lData[k]=-1;			
			#endif
				
			_DataSetFilter* theDF = (_DataSetFilter*)dataSetFilterList (dataPartitions.lData[k]);
			long thisModel = FindModelTemplate (((_String*)thisPart(0)),theDF->GetDimension());
			if (thisModel>=0)
			{
				partData.lData[k]=((_String*)thisPart(2))->toNum();
				partitionColors.lData[k]=((_String*)thisPart(3))->toNum();
				int model, options, freqs, classes;
				LongToModelData (((_String*)thisPart(1))->toNum(),model,options,freqs,classes);
				
				_List* modelInfo = FindModelTemplate ((_String*)thisPart(0));
				long   modelOption = ((_SimpleList*)(*modelInfo)(1))->lData[0];
				
				bool   allowedStates [5] = 
					  {(!(modelOption&HY_DATAPANEL_MODEL_MODELS))||(!(modelOption&(HY_DATAPANEL_MODEL_GLOBAL|HY_DATAPANEL_MODEL_GLOBALG))),
					   modelOption&HY_DATAPANEL_MODEL_GLOBAL,
					   modelOption&HY_DATAPANEL_MODEL_GLOBALG,
					   false,
					   false};
					   
				if ((options>2)||(!allowedStates[options]))
				{
					errMsg = _String ("Invalid parameter options supplied for model ") & *(_String*)thisPart(0);
					ProblemReport (errMsg);
					continue;
				}
				allowedStates[0] = !(modelOption&HY_DATAPANEL_MODEL_MODELS);
				allowedStates[1] = !(modelOption&HY_DATAPANEL_MODEL_MODELS);
				allowedStates[2] = !(modelOption&HY_DATAPANEL_MODEL_MODELS);
				allowedStates[3] = modelOption&HY_DATAPANEL_MODEL_EFVEST;
				allowedStates[4] = modelOption&HY_DATAPANEL_MODEL_MODELS;
				
				if ((freqs>4)||(!allowedStates[freqs]))
				{
					errMsg = _String ("Invalid frequency options supplied for model ") & *(_String*)thisPart(0);
					ProblemReport (errMsg);
					continue;
				}
					
				if (classes>HY_DATAPANEL_MAX_CLASSES)
				{
					errMsg = "Too many rate classes requested.";
					ProblemReport (errMsg);
					continue;
				}
					
				modelReferences.lData[k] = ModelDataToLong (thisModel+1,options,freqs,classes);
				
				_List			modelList;
				GenerateModelList (modelList,k);
				
				RefreshPartRow (modelList,k,true);
					
			}
			else
			{
				_List dummy;
				RefreshPartRow (dummy,k,true,false);
			}
		}
	}
	BuildThermometer();
	BuildMarksPane();
}

//__________________________________________________________
void	_HYDataPanel::RestorePanelSettings (_String* pInfo)
{
	_List parts;
	pInfo->StripQuotes();
	_ElementaryCommand::ExtractConditions (*pInfo,0,parts,',');
	if (parts.lLength==7)
	{
		_String* thisArg = (_String*)parts(0);
		_HYSequencePane*  sp = (_HYSequencePane*)GetObject (0);
		_DataSet*		  ds = (_DataSet*)dataSetList (dataSetID);
		long	 savedAL = thisArg->toNum();
		
		if (savedAL&HY_DATAPANEL_CONSENSUS)
			AdjustStatusLine (0,false,0);
		if (savedAL&HY_DATAPANEL_TRANSLATION)
		{
			thisArg = (_String*)parts(2);
			translatedSequence = thisArg->toNum();
			if ((translatedSequence>=0)&&(translatedSequence<ds->NoOfSpecies()))
			{
				AdjustStatusLine (2,false,translatedSequence+1);
			}
			else
				translatedSequence = -1;
			
		}
		if (savedAL&HY_DATAPANEL_REFERENCE)
		{
			thisArg = (_String*)parts(1);
			referenceSequence = thisArg->toNum();
			if ((referenceSequence>=0)&&(referenceSequence<ds->NoOfSpecies()))
			{
				AdjustStatusLine (3,false,referenceSequence+1);
			}
			else
				referenceSequence = -1;
			
		}
				
		/*thisArg = (_String*)parts(2);
		translatedSequence = thisArg->toNum();	*/
		
		//thisArg = (_String*)parts(3);
		//sp->startColumn = thisArg->toNum();	
		//if ((sp->startColumn<0)||(sp->startColumn>=siteAssignments.lLength))
		//	sp->startColumn = 0;
		thisArg = (_String*)parts(4);
		sp->SetNameDisplayMode (thisArg->toNum(),false);
		thisArg = (_String*)parts(5);
		sp->showDots = thisArg->toNum();
		thisArg = (_String*)parts(6);
		sp->blockWidth = thisArg->toNum();
		if ((sp->blockWidth!=9)&&(sp->blockWidth!=10))
			sp->blockWidth = 10;
		sp->BuildPane();
	}
}


//__________________________________________________________
void	_HYDataPanel::BuildThermometer (_HYRect * printRect)
{
	_HYStretchCanvas* therm = (_HYStretchCanvas*)GetObject (1);
	_HYSequencePane*  sp    = (_HYSequencePane*)GetObject (0);
	
	_HYRect canvasDim = printRect?*printRect:therm->GetCanvasSize(),
			backupTherm;
	
	canvasDim.width = 1;
	
	if (!printRect)
	{
		therm->StartDraw();
		therm->SetColor (labelBColor);
		therm->FillRect (canvasDim);
	}
	else
	{
		backupTherm = thermRect;
	}
	
	canvasDim.bottom--;
	therm->SetColor (black);
	
	if (!printRect)
		therm->DrawRect(canvasDim);
	
	thermRect = canvasDim;
	
	if (!printRect)
	{
		thermRect.left		+=	sp->headerWidth;
		thermRect.right		-=	HY_DATAPANEL_THERM_HSPACE;
		thermRect.top		+=	HY_DATAPANEL_THERM_VSPACE;
		thermRect.bottom	-=	HY_DATAPANEL_THERM_VSPACE+2;	
		therm->SetColor (thermFill);
	}
	else
		therm->SetColor (_hyWhiteColor);
	
	therm->FillRect (thermRect);
	therm->SetColor (black);
	therm->DrawRect(thermRect);
	
	long    pixelWidth 	 = thermRect.right-thermRect.left,
			k,
			lastColor = -1,
			thisAss,
			lastOverlap=0;
			
	double	pixelsPerSite = pixelWidth/(double)siteAssignments.lLength, 
			currentFPixel = 0.0;
			
	
	_HYRect siteRect = {thermRect.top+1,thermRect.left+1,thermRect.bottom-1,thermRect.left+1,1};
	
	//if (pixelsPerSite>=1.0)
	{
		for (k=0;k<siteAssignments.lLength;k++,currentFPixel+=pixelsPerSite)
		{
			thisAss = siteAssignments.lData[k];
			if ((thisAss!=lastColor)||((thisAss==lastColor)&&(lastColor==-2)&&(lastOverlap!=overlaps.lData[k])))
			{
				siteRect.right = thermRect.left+1+currentFPixel;
				if (lastColor>=0)
				{
					long thisColor = partitionColors.lData[lastColor];
					therm->SetColor (LongToHYColor(thisColor));
					therm->FillRect (siteRect);
				}
				else
				if (lastColor == -2)
				{
					//therm->SetColor(black);
					//therm->FillRect (siteRect);
					long p4 = (lastOverlap&0xFF000000)>>24,
				     	 p3 = (lastOverlap&0x00FF0000)>>16,
				     	 p2 = (lastOverlap&0x0000FF00)>>8,
				     	 p1 = (lastOverlap&0x000000FF),
				     	 smallH;
				     	 
				    _HYRect smallRect = siteRect;
				    
				    if(p4>0) // four overlaps
				    {
				    	smallH = (smallRect.bottom-smallRect.top)/4;
				    	smallRect.bottom = smallRect.top+smallH;
				    	therm->SetColor(LongToHYColor(partitionColors.lData[p4]));
				    	therm->FillRect(smallRect);
				    	therm->SetColor(LongToHYColor(partitionColors.lData[p3]));
				    	smallRect.top = smallRect.bottom;
				    	smallRect.bottom = smallRect.top+smallH;
				    	therm->FillRect(smallRect);
				    }
				    else
					    if(p3>0)
					    {
					    	smallH = (smallRect.bottom-smallRect.top)/3;
					    	smallRect.bottom = smallRect.top+smallH;
					    	therm->SetColor(LongToHYColor(partitionColors.lData[p3]));
					    	therm->FillRect(smallRect);
					    }
					    else
					    {
					    	smallH = (smallRect.bottom-smallRect.top)/2;
					    	smallRect.bottom=smallRect.top;
					    }
			    	therm->SetColor(LongToHYColor(partitionColors.lData[p2]));
			    	smallRect.top = smallRect.bottom;
			    	smallRect.bottom = smallRect.top+smallH;
			    	therm->FillRect(smallRect);
			    	therm->SetColor(LongToHYColor(partitionColors.lData[p1]));
			    	smallRect.top = smallRect.bottom;
			    	smallRect.bottom = siteRect.bottom;
			    	therm->FillRect(smallRect);
				}
				lastColor = thisAss;
				siteRect.left = siteRect.right;
				lastOverlap = overlaps.lData[k];
			}
		}
		siteRect.right = thermRect.right-1;
		if (lastColor>=0)
		{
			long thisColor = partitionColors.lData[lastColor];
			therm->SetColor (LongToHYColor(thisColor));
			therm->FillRect (siteRect);
		}
		else
		{
			if (lastColor == -2)
			{
				long p4 = (lastOverlap&0xFF000000)>>24,
			     	 p3 = (lastOverlap&0x00FF0000)>>16,
			     	 p2 = (lastOverlap&0x0000FF00)>>8,
			     	 p1 = (lastOverlap&0x000000FF),
			     	 smallH;
			     	 
			    _HYRect smallRect = siteRect;
			    
			    if(p4>0) // four overlaps
			    {
			    	smallH = (smallRect.bottom-smallRect.top)/4;
			    	smallRect.bottom = smallRect.top+smallH;
			    	therm->SetColor(LongToHYColor(partitionColors.lData[p4]));
			    	therm->FillRect(smallRect);
			    	therm->SetColor(LongToHYColor(partitionColors.lData[p3]));
			    	smallRect.top = smallRect.bottom;
			    	smallRect.bottom = smallRect.top+smallH;
			    	therm->FillRect(smallRect);
			    }
			    else
			    if(p3>0)
			    {
			    	smallH = (smallRect.bottom-smallRect.top)/3;
			    	smallRect.bottom = smallRect.top+smallH;
			    	therm->SetColor(LongToHYColor(partitionColors.lData[p3]));
			    	therm->FillRect(smallRect);
			    }
			    else
			    {
			    	smallH = (smallRect.bottom-smallRect.top)/2;
			    	smallRect.bottom=smallRect.top;
			    }
		    	therm->SetColor(LongToHYColor(partitionColors.lData[p2]));
		    	smallRect.top = smallRect.bottom;
		    	smallRect.bottom = smallRect.top+smallH;
		    	therm->FillRect(smallRect);
		    	therm->SetColor(LongToHYColor(partitionColors.lData[p1]));
		    	smallRect.top = smallRect.bottom;
		    	smallRect.bottom = siteRect.bottom;
		    	therm->FillRect(smallRect);			
		    }
		}
	}

	if (printRect)
	{
		_SimpleList * p_coords,
					* p_count;
					
		if (dataPartitions.lLength)
		{
			checkPointer(p_coords = new _SimpleList (dataPartitions.lLength,0,0));
			checkPointer(p_count  = new _SimpleList (dataPartitions.lLength,0,0));
			for (k=0;k<siteAssignments.lLength;k++)
			{
				thisAss = siteAssignments.lData[k];
				if (thisAss==-2) // overlap
					break;
				if (thisAss>=0)
				{
					p_coords->lData[thisAss] += k;
					p_count->lData[thisAss]  ++;		
				}
			}
			if (k==siteAssignments.lLength)
			{
				for (k=0; k<dataPartitions.lLength; k++)
				{
					currentFPixel = (pixelsPerSite*p_coords->lData[k])/p_count->lData[k]+0.5*pixelsPerSite;
					_HYFont labelFont = {"Times", (thermRect.bottom-thermRect.top)*2/3,HY_FONT_PLAIN};
					_String partName = ((_String*)dataSetFilterNamesList (dataPartitions.lData[k]))->Replace ("_"," ",true);
					thisAss = GetVisibleStringWidth (partName,labelFont)/2+1;
					
					_HYColor partColor = LongToHYColor(partitionColors.lData[lastColor]);
					if (MAX(partColor.R,MAX(partColor.B,partColor.G))<=128)
						therm->SetColor (_hyWhiteColor);
					else
						therm->SetColor (black);
					
					therm->SetFont (labelFont);
					therm->DisplayText (partName, thermRect.top+labelFont.size*7/6, thermRect.left+currentFPixel-thisAss,true);
				}
			}
			DeleteObject (p_coords);
			DeleteObject (p_count);
		}
	}

	if (!printRect)
	{
		therm->EndDraw();
		therm->_MarkForUpdate();
	}
	else
	{
		thermRect = backupTherm;
	}
}

//__________________________________________________________
void	_HYDataPanel::BuildMarksPane (void)
{
	_HYStretchCanvas* marks = (_HYStretchCanvas*)GetObject (3);
	_HYSequencePane*  sp 	= (_HYSequencePane*) GetObject (0);
	
	_HYRect canvasDim = marks->GetCanvasSize(), 
			sRect = {canvasDim.bottom-2,0,canvasDim.bottom-2,canvasDim.right,1};
	
	marks->StartDraw();
	marks->EraseAll();		

	marks->SetColor(black);
	marks->DrawLine (sRect);
	
	long  k,
		  a=-1,
		  b,
		  c=-1,
		  d;
		  		  
	sRect.top     = canvasDim.top;
	sRect.bottom -= 2;
	sRect.left = sp->headerWidth;
	for (k=sp->startColumn;k<sp->endColumn-1;k++)
	{
		b = siteAssignments.lData[k];
		d = overlaps.lData[k];
		if ((a!=b)||((a==-2)&&(c!=d)))
		{
			if (a>=0)
			{
				marks->SetColor (LongToHYColor(partitionColors.lData[a]));
				sRect.right		=	sRect.left;
				sRect.top 		=   0;
				sRect.bottom	=	canvasDim.bottom-3;
				marks->DrawLine (sRect);
				sRect.right=(sRect.left--);
				sRect.top++;
				marks->DrawLine (sRect);
				sRect.right=(sRect.left--);
				sRect.top++;
				sRect.bottom--;
				marks->DrawLine (sRect);
				sRect.left+=2;
			}
			if (b>=0)
			{
				marks->SetColor (LongToHYColor(partitionColors.lData[b]));
				sRect.right=(sRect.left+=2);
				sRect.top = 0;
				sRect.bottom=canvasDim.bottom-3;
				marks->DrawLine (sRect);
				sRect.right=(sRect.left++);
				sRect.top++;
				marks->DrawLine (sRect);
				sRect.right=(sRect.left++);
				sRect.top++;
				sRect.bottom--;
				marks->DrawLine (sRect);
				sRect.left-=4;					
			}
		}
		sRect.left+=sp->charWidth;
		if (k&&(k%sp->blockWidth==0)) sRect.left+=2;
		a=b;
		c=d;
	}
	if (k==siteAssignments.lLength-1)
	{
		b = siteAssignments.lData[k];
		sRect.left+=sp->charWidth;
		if (k&&(k%sp->blockWidth==0)) sRect.left+=2;
		if (b>=0)
		{
			marks->SetColor (LongToHYColor(partitionColors.lData[b]));
			sRect.right=sRect.left;
			sRect.top = 0;
			sRect.bottom=canvasDim.bottom-3;
			marks->DrawLine (sRect);
			sRect.right=(sRect.left--);
			sRect.top++;
			marks->DrawLine (sRect);
			sRect.right=(sRect.left--);
			sRect.top++;
			sRect.bottom--;
			marks->DrawLine (sRect);
			sRect.left+=2;
		}
	}	
	marks->EndDraw();
	marks->_MarkForUpdate();	
}	

//__________________________________________________________
_HYRect  _HYDataPanel::ComputeNavRect (void)
{
	_HYRect 		   res;
	_HYSequencePane*   seqPane 		= (_HYSequencePane*)	GetObject (0);
	
	res.width = 1;
	
	_Parameter	nS = thermRect.right-thermRect.left+1;

	_Parameter visProp = (seqPane->endColumn-seqPane->startColumn)/(_Parameter)seqPane->columnStrings.lLength;
	res.left  = seqPane->startColumn*nS/(_Parameter)seqPane->columnStrings.lLength;
	res.right = res.left+visProp*nS;
	
	if (res.right<res.left+3) 
		res.right = res.left+3;
	
	visProp    = (seqPane->endRow-seqPane->startRow)/(_Parameter)seqPane->RowCount();
	res.top    = 0;
	res.bottom = HY_DATAPANEL_THERMWIDTH;
	return res;
}
//__________________________________________________________
void  _HYDataPanel::SetNavRectCenter (long h, long)
{
	_HYSequencePane * bigCanvas = (_HYSequencePane*)GetObject (0);
			 
	long		nS = thermRect.right-thermRect.left+1,
				t;
				
	t 	= navRect.right-navRect.left;
	h  -= t/2;
	
	if (h<0)
		h=0;

	if (h+t>nS)
		h=nS-t;

	h  = (h*(_Parameter)bigCanvas->columnStrings.lLength)/nS;
	if (h<0) 
		h = 0;
		
	if (h!=bigCanvas->startColumn)
	{
		bigCanvas->ProcessEvent(generateScrollEvent(h-bigCanvas->startColumn,0));
		ProcessEvent (generateScrollEvent(0,0));
	}
}

//__________________________________________________________
long  _HYDataPanel::FindUnusedColor (void)
{
	long tryColor;
	for (long d=1; d<4; d++)
	{
		for (long k=0;k<HY_DATAPANEL_DEF_COLORS;k++)
		{
			_HYColor thisColor = dataColors[k];
			switch (d)
			{
				case 2:
					thisColor.R/=1.5;
					thisColor.G/=1.5;
					thisColor.B/=1.5;
					break;
				case 3:
					thisColor.R = (thisColor.R+0x0000ffff)/2;
					thisColor.G = (thisColor.G+0x0000ffff)/2;
					thisColor.B = (thisColor.B+0x0000ffff)/2;
					break;
			}
			tryColor = HYColorToLong (thisColor);
			if (partitionColors.Find (tryColor)<0)
				return tryColor;
		}
	}
	_HYColor bailOut = {0,0,0};
	_String bailOutString ("Select a color for the partition");
	return HYColorToLong (SelectAColor(bailOut,bailOutString));
}

//__________________________________________________________
void  _HYDataPanel::MarkSites (_SimpleList& siteList,long k)
{
	long p;
	for (long j=0;j<siteList.lLength;j++)
	{
		long index = siteList.lData[j];

		p = siteAssignments.lData[index];
		if (p==-1)
		{
			siteAssignments.lData[index]=k;
		}
		else
		if (p>=0)
		{
			if (k<p)
				overlaps.lData[index] = (p<<8)+k;
			else
				overlaps.lData[index] = (k<<8)+p;
			siteAssignments.lData[index]=-2;
		}
		else
		{
			p = overlaps.lData[index];
			long p3 = (p&0x00FF0000)>>16, 
			     p2 = (p&0x0000FF00)>>8, 
			     p1 = p&0x000000FF,
			     p4;
			if (p3&&(k>p3))
				p4 = k;
			else
			if (p2&&(k>p2))
			{
				p4 = p3;
				p3 = k;
			}
			else
			if (k>p1)
			{
				p4=p3;
				p3=p2;
				p2=k;
			}
			else
			{
				p4=p3;
				p3=p2;
				p2=p1;
				p1=k;
			}
			overlaps.lData[index]=(p4<<24)+(p3<<16)+(p2<<8)+p1;		
		}
	}
}
//__________________________________________________________
void  _HYDataPanel::NavBarDblClick (long k)
{
	_HYTable* 		  dl    = (_HYTable*)GetObject (5);
	double	pixelsPerSite = (thermRect.right-thermRect.left)/(double)siteAssignments.lLength;
	long	index,t;
	k/=pixelsPerSite;
	_SimpleList newSelection;
	index = siteAssignments.lData[k];
	if (index>=0)
		newSelection<<index;
	else
		if (index==-2)
		{
			index = overlaps.lData[k];
			t = (index&0x000000FF);
			newSelection<<t;
			index=index>>8;
			t=index%256;
			while(t&&index)
			{
				newSelection<<t;
				index=index>>8;
				t=index%256;
			}
		}
	if (newSelection.lLength)
	{
		_SimpleList oldSelection;
		dl->GetRowSelection (oldSelection);
		
		if (!oldSelection.Equal(newSelection))
		{
			dl->SetRowSelection(newSelection);
			UpdatePartitionOperations();
		}
		if (newSelection.lLength==1)
		{
			EditPartitionProperties (newSelection.lData[0]);
		}
	}
}
//__________________________________________________________

void  _HYDataPanel::GenerateTreeList (_List& trees)
{
	long		 k, 
				 mT;
	
	if (dataWrapper)
		mT = dataWrapper->NumberSpecies();
	else
		mT = ((_DataSet*)dataSetList(dataSetID))->NoOfSpecies();
		
	trees && & none;
	trees && & makeNewTree;
	trees && & readTreeFile;
	trees && & inferTreeStr;
	
	for (k=0; k<variablePtrs.lLength;k++)
	{
		_Variable* thisVar = LocateVar(k);
		if (thisVar&&(thisVar->ObjectClass()==TREE))
		{
			_TheTree* thisTree = (_TheTree*)thisVar;
			_PMathObj tipCount = thisTree->TipCount();
			if (tipCount->Value()==mT || (mT==2 && tipCount->Value()==1))
			{
				if (trees.lLength==4)
					trees << & menuSeparator;
				trees << thisTree->GetName();
			}
			DeleteObject (tipCount);
		}
	}
	
}

//__________________________________________________________
void  _HYDataPanel::GenerateModelList (_List& modelList, long partID)
{

	long  mT;
		
	_DataSetFilter * theDF = (_DataSetFilter*)dataSetFilterList (dataPartitions.lData[partID]);
	
	mT = theDF->GetDimension();
	
	modelList && & none;
	
	for (long k=0; k<modelTemplates.lLength; k++)
	{
		_List * thisModel = (_List*)modelTemplates (k);
		_SimpleList* modelOptions = (_SimpleList*)(*thisModel)(1);
		if ((modelOptions->lData[1]==mT)||((modelOptions->lData[1]==64)&&(theDF->GetUnitLength()==3)))
			modelList <<  (_String*)(*thisModel)(0);
	}
}

//__________________________________________________________
void  _HYDataPanel::GenerateModelPOptionList (_List& optionsList, long partID)
{
	_DataSetFilter* thisPartition = (_DataSetFilter*)dataSetFilterList (dataPartitions.lData[partID]);
	_List* 	   		theModel = FindModelTemplate ((HY_DATAPANEL_MODELID&modelReferences.lData[partID])-1,thisPartition->GetDimension());

	if (theModel)
	{
		_SimpleList* opt = (_SimpleList*)(*theModel)(1);
		
		if (!(opt->lData[0]&HY_DATAPANEL_MODEL_MODELS))
		{
			optionsList && & parameterOption[0];
			if (opt->lData[0]&HY_DATAPANEL_MODEL_GLOBAL)
				optionsList && & parameterOption[1];
				
			if (opt->lData[0]&HY_DATAPANEL_MODEL_GLOBALG)
				optionsList && & parameterOption[2];
		}
		else
		{
			partID = 0;
			if (opt->lData[0]&HY_DATAPANEL_MODEL_GLOBALG)
				partID = 2;
			else
				if (opt->lData[0]&HY_DATAPANEL_MODEL_GLOBAL)
					partID = 1;
			optionsList && & parameterOption[partID];
		}
	}
}

//__________________________________________________________
void  _HYDataPanel::GenerateModelFOptionList (_List& optionsList, long partID)
{
	_DataSetFilter* thisPartition = (_DataSetFilter*)dataSetFilterList (dataPartitions.lData[partID]);
	_List* 	   theModel = FindModelTemplate ((modelReferences.lData[partID]&HY_DATAPANEL_MODELID)-1,thisPartition->GetDimension());

	if (theModel)
	{
		_SimpleList* opt = (_SimpleList*)(*theModel)(1);
		
		if (!(opt->lData[0]&HY_DATAPANEL_MODEL_MODELS))
			for (partID=0;partID<4; partID++)
				optionsList && & freqOption[partID];
		else
			optionsList && & freqOption[4];
	}
}

//__________________________________________________________
void  _HYDataPanel::UnmarkSites (_SimpleList& sL,long ind)
{
	_SimpleList		siteList;
	
	siteList.Duplicate (&sL);
	siteList.Sort();
	long p;
	for (long j=0;j<siteList.lLength;j++)
	{
		long index = siteList.lData[j];
		p = siteAssignments.lData[index];
		if (p==ind)
		{
			siteAssignments.lData[index]=-1;
		}
		else
		if (p==-2)
		{
			p = overlaps.lData[index];
			long p3 = (p&0x00FF0000)>>16, 
			     p2 = (p&0x0000FF00)>>8, 
			     p1 = p&0x000000FF,
			     p4 = (p&0xFF000000)>>24;
			if (ind==p1)
			{
				p1 = p2;
				p2 = p3;
				p3 = p4;
				p4 = 0;
			}
			else
			if (ind==p2)
			{
				p2=p3;
				p3=p4;
				p4=0;
			}
			else
			if (p3&&(ind==p3))
			{
				p3=p4;
				p4=0;
			}
			else
			if (p4&&(ind==p4))
			{
				p4=0;
			}
			
			p=(p4<<24)+(p3<<16)+(p2<<8)+p1;
			if (p<256)
			{
				siteAssignments.lData[index]=p;
				overlaps.lData[index]=0;
			}			
			else
				overlaps.lData[index]=p;
		}
	}
}

//__________________________________________________________
void  _HYDataPanel::CorrectSites (long ind)
{
	long p;
	for (long j=0;j<siteAssignments.lLength;j++)
	{
		p = siteAssignments.lData[j];
		if (p>ind)
		{
			siteAssignments.lData[j]--;
		}
		else
		if (p==-2)
		{
			p = overlaps.lData[j];
			long p3 = (p&0x00FF0000)>>16, 
			     p2 = (p&0x0000FF00)>>8, 
			     p1 = p&0x000000FF,
			     p4 = (p&0xFF000000)>>24;
			     
			if (p4>ind) p4--;
			if (p3>ind) p3--;
			if (p2>ind) p2--;
			if (p1>ind) p1--;
			
			p=(p4<<24)+(p3<<16)+(p2<<8)+p1;
			overlaps.lData[j]=p;
		}
	}
}
//__________________________________________________________
bool  _HYDataPanel::IsSelectionNonEmpty (void)
{
	return ((_HYSequencePane*)GetObject(0))->selection.lLength>0;
}

//__________________________________________________________

void  _HYDataPanel::UpdateConsensusSequence (_String& newDataString, bool applyChanges)
{
	if (dataWrapper)
		newDataString = dataWrapper->GenerateConsensusString();
	else
	{
		_DataSetFilter temp;
		_SimpleList	   empty1, empty2;
		temp.SetFilter ((_DataSet*)dataSetList (dataSetID), 1, empty1, empty2, false);
		newDataString = temp.GenerateConsensusString()&"    ";
	}
	
	if (applyChanges)
	{
		if (addedLines&HY_DATAPANEL_CONSENSUS)
		{
			_HYSequencePane* sp2 = (_HYSequencePane*)components(4);
			for (long k=0; k<newDataString.sLength; k++)
				((_String*)statusData(k))->sData[0] = newDataString[k];
			sp2->BuildPane();
			sp2->_MarkForUpdate();
		}
	}
}

//__________________________________________________________

void  _HYDataPanel::UpdateTranslationString (_String& newDataString, long index, bool applyChanges, long genCodeRef)
{	
	// go filter by filter; take each codon filter, apply properties and translate
	
	_Parameter     *freqVector = (_Parameter*)MemAllocate(sizeof(_Parameter)*64);
	long	k;
	
	_String spaces ((unsigned long)(siteAssignments.lLength+2), false);
	
	for (k=0; k<siteAssignments.lLength+4; k++)
		spaces.sData[k] = ' ';
		
	newDataString = spaces;
	
	spaces = "   ";
	
	if (genCodeRef == -1)
		for (k=0; k<dataPartitions.lLength; k++)
		{
			_DataSetFilter* df = (_DataSetFilter*)dataSetFilterList(dataPartitions.lData[k]);
			if (df->GetUnitLength() == 3) // a codon filter; proceed
			{
				_SimpleList		backup;
				backup.Duplicate (&df->theExclusions);
				df->theExclusions.Clear();
				bool rev;
				char offset;
				long genCode;
				LongToPartData (partData.lData[k],offset,rev,genCode);
				_SimpleList* gencode = ((_SimpleList*)((*(_List*)geneticCodes(genCode))(2)));
				if (rev) // reverse direction
					for (long m=df->theOriginalOrder.lLength-offset-1; m>=2; m-=3)
					{
						spaces.sData[0] = (*df->GetData())(df->theOriginalOrder.lData[m],index,1);
						spaces.sData[1] = (*df->GetData())(df->theOriginalOrder.lData[m-1],index,1);
						spaces.sData[2] = (*df->GetData())(df->theOriginalOrder.lData[m-2],index,1);
						CodeTo3AA (spaces,df->Translate2Frequencies(spaces,freqVector,false),gencode,freqVector);
						newDataString.sData[df->theOriginalOrder.lData[m-2]]=spaces.sData[0];
						newDataString.sData[df->theOriginalOrder.lData[m-1]]=spaces.sData[1];
						newDataString.sData[df->theOriginalOrder.lData[m]]=spaces.sData[2];
					}
				else
					for (long m=offset; m<df->theOriginalOrder.lLength-2; m+=3)
					{
						spaces.sData[0] = (*df->GetData())(df->theOriginalOrder.lData[m],index,1);
						spaces.sData[1] = (*df->GetData())(df->theOriginalOrder.lData[m+1],index,1);
						spaces.sData[2] = (*df->GetData())(df->theOriginalOrder.lData[m+2],index,1);
						CodeTo3AA (spaces,df->Translate2Frequencies(spaces,freqVector,false),gencode,freqVector);
						newDataString.sData[df->theOriginalOrder.lData[m]]=spaces.sData[0];
						newDataString.sData[df->theOriginalOrder.lData[m+1]]=spaces.sData[1];
						newDataString.sData[df->theOriginalOrder.lData[m+2]]=spaces.sData[2];
					}
				df->theExclusions.Duplicate (&backup);
			}	
		}
	else
	{
		_DataSet * ds = (_DataSet*)dataSetList (dataSetID);
		
		_DataSetFilter dummy;
		_SimpleList	   e1,
					   e2;
					   
		e1<<0;
		e1<<1;
		e1<<2;
					   
		dummy.SetFilter (ds,3, e1, e2);
		_SimpleList* gencode = ((_SimpleList*)((*(_List*)geneticCodes(genCodeRef))(2)));
		
		for (k=0; k<siteAssignments.lLength-2; k+=3)
		{
			spaces.sData[0] = (*ds)(k,index,1);
			spaces.sData[1] = (*ds)(k+1,index,1);
			spaces.sData[2] = (*ds)(k+2,index,1);
			CodeTo3AA (spaces,dummy.Translate2Frequencies(spaces,freqVector,false),gencode,freqVector);
			newDataString.sData[k]=spaces.sData[0];
			newDataString.sData[k+1]=spaces.sData[1];
			newDataString.sData[k+2]=spaces.sData[2];
		}
		
		for (; k<siteAssignments.lLength; k++)
			newDataString.sData[k] = ' ';
	}
	
	free (freqVector);
	
	if (applyChanges)
	{
		if (addedLines&HY_DATAPANEL_TRANSLATION)
		{
			long idx = ((_String*)statusData(0))->sLength-1;
			if (addedLines&HY_DATAPANEL_REFERENCE)
				idx--;
			_HYSequencePane* sp2 = (_HYSequencePane*)components(4);
			for (long k=0; k<newDataString.sLength; k++)
				((_String*)statusData(k))->sData[idx] = newDataString[k];
			sp2->BuildPane();
			sp2->_MarkForUpdate();
		}
	}
}

//__________________________________________________________

bool  _HYDataPanel::GetTranslationString (_String& newDataString, long index, char resolve, long filterID, _List * cachedResolutions)
{	
	// go filter by filter; take each codon filter, apply properties and translate
	
	_Parameter     *freqVector = (_Parameter*)MemAllocate(sizeof(_Parameter)*64);
	long	k,
			index2;
			
	_String newData (16, true), 
			spaces  (3, false);
	
	for (k=filterID>=0?filterID:0; k<dataPartitions.lLength; k++)
	{
		_DataSetFilter* df = (_DataSetFilter*)dataSetFilterList(dataPartitions.lData[k]);
		_String			gaps (3,false);
		gaps.sData[0] = df->GetData()->GetTT()->GetGapChar();
		gaps.sData[1] = gaps.sData[0];
		gaps.sData[2] = gaps.sData[0];
		
		if (df->GetUnitLength() == 3) // a codon filter; proceed
		{
			if 	(resolve)
				index2 = df->theNodeMap.Find (index);
			
			_SimpleList		backup;
			bool rev;
			char offset;
			long genCode;
			LongToPartData (partData.lData[k],offset,rev,genCode);
			
			_SimpleList* 	gencode 		= ((_SimpleList*)((*(_List*)geneticCodes(genCode))(2)));

			_List			*thisFilterMap = nil;
			
			if (resolve > 0)
			{
				if (cachedResolutions && cachedResolutions->lLength <= k)
				{
					_List * thisFilterResolution = new _List;
					for (long m=offset; m<df->theOriginalOrder.lLength-2; m+=3)
						thisFilterResolution->AppendNewInstance(df->CountAndResolve(df->duplicateMap[df->theOriginalOrder.lData[m]/3],freqVector,resolve-1));
					
					cachedResolutions->AppendNewInstance (thisFilterResolution);					
				}
				thisFilterMap = (_List*)(*cachedResolutions)(k);
			}
			else
			{
				backup.Duplicate (&df->theExclusions);
				df->theExclusions.Clear();
			}	
			
			if (rev) // reverse direction
				//for (long m=2+(df->theOriginalOrder.lLength-offset)%3; m<df->theOriginalOrder.lLength-offset; m+=3)
				for (long m=df->theOriginalOrder.lLength-offset-1; m>=2; m-=3)
				{
					spaces.sData[0] = (*df->GetData())(df->theOriginalOrder.lData[m],index,1);
					spaces.sData[1] = (*df->GetData())(df->theOriginalOrder.lData[m-1],index,1);
					spaces.sData[2] = (*df->GetData())(df->theOriginalOrder.lData[m-2],index,1);
					if (resolve && thisFilterMap)
						newData << CodeToAA (df->CorrectCode(((_SimpleList*)thisFilterMap->lData[m/3])->lData[index2]),gencode);
					else
						if (spaces.Equal(&gaps))	
							newData << '-';
						else
							newData << CodeToAA (df->Translate2Frequencies(spaces,freqVector,false),gencode,freqVector);
							
				}
			else
				for (long m=offset; m<df->theOriginalOrder.lLength-2; m+=3)
				{
					spaces.sData[0] = (*df->GetData())(df->theOriginalOrder.lData[m],index,1);
					spaces.sData[1] = (*df->GetData())(df->theOriginalOrder.lData[m+1],index,1);
					spaces.sData[2] = (*df->GetData())(df->theOriginalOrder.lData[m+2],index,1);
					if (resolve && thisFilterMap)
						newData << CodeToAA (df->CorrectCode(((_SimpleList*)thisFilterMap->lData[m/3])->lData[index2]),gencode);
					else
						if (spaces.Equal(&gaps))	
							newData << '-';
						else
							newData << CodeToAA (df->Translate2Frequencies(spaces,freqVector,false),gencode,freqVector);
				}
			if (resolve == 0)
				df->theExclusions.Duplicate (&backup);
			if (filterID>=0)
			{
				free (freqVector);
				newData.Finalize();
				newDataString = newData;
				return rev;
			}
		}
	}
	
	free (freqVector);
	newData.Finalize();
	newDataString = newData;
	
	return true;
}

//__________________________________________________________
void  _HYDataPanel::CreatePartition (_SimpleList& siteList, char unitLength, bool jump2New, _String* prefix)
{
	if (dataSetID>=0)
	{
		bool 			good = true;
		if ((partitionColors.lLength==255)&&(siteAssignments.Find(-2)>=0))
			good = false;
		else
			for (long k=0;k<siteList.lLength;k++)
				if (overlaps.lData[siteList.lData[k]]>0x00FFFFFF)
				{
					good = false;
					break;
				}
		if (!good)
		{
			_String errMsg = _String("Couldn't create new partition. ")&overlapWarning;
			ProblemReport (errMsg);
			return;
		}
		_DataSetFilter* newDF = new _DataSetFilter();
		_SimpleList		emptyList;
		if (dataWrapper)
			emptyList.Duplicate (&dataWrapper->theNodeMap);
		newDF->SetFilter((_DataSet*)dataSetList(dataSetID),unitLength,emptyList,siteList,false);
		if (terminateExecution)
		{
			terminateExecution = false;
			return;
		}
		_String filterName;
		if (!prefix)
		 	filterName = (*dataSetName)& "_part";
		else
			filterName = *prefix;
			
		dataPartitions<<AddFilterToList (filterName,newDF);
		treeVarReferences<<-1;
		modelReferences<<0;
		partData<<0;
		partitionColors<<FindUnusedColor();
		long nr = AddPartitionRow (&filterName, newDF->GetUnitLength());
		MarkSites (siteList,partitionColors.lLength-1);
		if (jump2New)
		{
			_HYTable* tl = (_HYTable*)GetObject (5);
			_SimpleList sl (nr);
			tl->SetRowSelection(sl);
			tl->ScrollToRow (nr);
			UpdatePartitionOperations();
		}
		tainted = true;
		BuildThermometer();
		BuildMarksPane();
	}
}

//__________________________________________________________
long  _HYDataPanel::AddPartitionRow (_String* partName, long partUnit)
{
	_HYTable * pT = (_HYTable*)GetObject (5);
	long	 k = pT->verticalSpaces.lLength;
	
	bool	 hasPadding = (((_String*)pT->GetCellData (0,k-1))->sLength==0);
			
	if (hasPadding)
	{
		k--;
		long s;
		if ((s=pT->GetRowSpacing(k))<20)
			pT->DeleteRow (k);
		else
			pT->SetRowSpacing (k,-20,false);
	}
	pT->AddRow 		(k,20,HY_TABLE_STATIC_TEXT);
	pT->SetCellData (partName,k,DF_ID_COLUMN,HY_TABLE_STATIC_TEXT,true);
	pT->SetCellData (&none,k,DF_TREE_COLUMN,HY_TABLE_STATIC_TEXT|HY_TABLE_PULLDOWN,true);
	pT->SetCellData (&none,k,DF_MODEL_COLUMN,HY_TABLE_STATIC_TEXT|HY_TABLE_PULLDOWN,true);
	_String	* dataTypeS = &unknownDataType;
	
	if (dataType&HY_DATAPANEL_NUCDATA)
		switch (partUnit)
		{
			case 1:
				dataTypeS = &nucDataType;
				break;
			case 2:
				dataTypeS = &dinucDataType;
				break;
			case 3:
				dataTypeS = &codonDataType;
				break;
		}
	else
		if (dataType&HY_DATAPANEL_PROTDATA)
			dataTypeS = &proteinDataType;
		else
			if (dataType&HY_DATAPANEL_BINARYDATA)
			{
				if (partUnit == 1)
					dataTypeS = &binaryDataType;
				else
					dataTypeS = &disequonDataType;
			}	
			else	
				dataTypeS = &unknownDataType;
	
	pT->SetCellData (dataTypeS,k,DF_TYPE_COLUMN,HY_TABLE_STATIC_TEXT|HY_TABLE_PULLDOWN,true);
	
	_SimpleList	colorRect;
	
	colorRect << partitionColors.lData[partitionColors.lLength-1];
	colorRect << 13;
	colorRect << 13;
	colorRect << HY_TABLE_COLOR_CIRCLE;
	
	pT->SetCellData (&colorRect,k,DF_COLOR_COLUMN,HY_TABLE_ICON,true);

	pT->_MarkContentsForUpdate();
	pT->SetVisibleSize (pT->rel);
	return k;
}


//__________________________________________________________
void  _HYDataPanel::DeletePartitionRow (long index)
{
	_HYTable * pT = (_HYTable*)GetObject (5);
	
	pT->DeleteRow(index);

	long 		h2 = pT->verticalSpaces.lData[pT->verticalSpaces.lLength-1];
	bool		hasPadding = (((_String*)pT->GetCellData (0,pT->verticalSpaces.lLength-1))->sLength == 0);
	
	if (h2<100)
	{
		if (hasPadding)
			pT->SetRowSpacing (pT->verticalSpaces.lLength-1,100-h2,false);
		else
			pT->AddRow (-1,100-h2,HY_TABLE_STATIC_TEXT|HY_TABLE_CANTSELECT);
	}
	pT->SetVisibleSize (pT->rel);
	pT->ScrollToRow(index>0?index-1:0);
	pT->_MarkForUpdate();
}

//__________________________________________________________
void  _HYDataPanel::AddPartition (_DataSetFilter* theDF)
{
	if (dataSetID>=0)
	{
		bool 			good = true;
		if ((partitionColors.lLength==255)&&(siteAssignments.Find(-2)>=0))
			good = false;
		else
			for (long k=0;k<theDF->theOriginalOrder.lLength;k++)
				if (overlaps.lData[theDF->theOriginalOrder.lData[k]]>0x00FFFFFF)
				{
					good = false;
					break;
				}
		if (!good)
		{
			_String errMsg = _String("Couldn't add new partition. ")&overlapWarning;
			ProblemReport (errMsg);
			return;
		}
		long			filterID = dataSetFilterList._SimpleList::Find ((long)theDF);
		if (filterID>=0)
		{
			_String*	filterName = (_String*)dataSetFilterNamesList (filterID);
			dataPartitions<<filterID;
			treeVarReferences<<-1;
			modelReferences<<0;
			/*filterID = 0;
			if (theDF->GetUnitLength () == 3) // codon
			{
				_String excl;
				DFExclusionsToString	  (theDF,excl);
				filterID = FindGeneticCodeByExcl (&excl);
				if (filterID<0)
					filterID = 0;
			}
			partData<<PartDataToLong (0,0,filterID);*/
			partData << 0;
			partitionColors<<FindUnusedColor();
			AddPartitionRow (filterName, theDF->GetUnitLength());
			MarkSites (theDF->theOriginalOrder,partitionColors.lLength-1);
			BuildThermometer();
			BuildMarksPane();
		}
	}
}

//__________________________________________________________

_String  _HYDataPanel::SavePartition (long index, _String* dsFile, _DataSet* dsXtra)
{
	/* stuff to add:
	
		window-based width?
		save attached tree
		
	*/
	if (dfFilterFormats.lLength==0)
	{
		_String format ("# sequential");
		dfFilterFormats && &format;
		format = "# interleaved";
		dfFilterFormats && &format;
		format = "PHYLIP Sequential";
		dfFilterFormats && &format;
		format = "PHYLIP Interleaved";
		dfFilterFormats && &format;
		format = "NEXUS Sequential[Labels]";
		dfFilterFormats && &format;
		format = "NEXUS Interleaved[Labels]";
		dfFilterFormats && &format;
		format = "NEXUS Sequential[No Labels]";
		dfFilterFormats && &format;
		format = "NEXUS Interleaved[No Labels]";
		dfFilterFormats && &format;
		format = "Comma Separated Character Data";
		dfFilterFormats && &format;
		format = "FASTA sequential";
		dfFilterFormats && &format;
	}
	if (dataSetID>=0)
	{
		bool		cleanUp = false;
		if (index>=(long)dataPartitions.lLength) return empty;
		long dfID, formatChoice;
		_DataSetFilter *theDF;
		_HYSequencePane* sp = (_HYSequencePane*)GetObject(0);
		
		if (index>=0)
		{
			dfID = dataPartitions.lData[index];
			theDF = (_DataSetFilter*)dataSetFilterList(dfID);
			_HYButtonBar*   	bb = (_HYButtonBar*)GetObject(2);
			bb->_UnpushButton();
			if (treeVarReferences.lData[index]>=0)
			{
				stashParameter (noInternalLabels,1.0,true);
				_String* treeString = (_String*)((_TheTree*)LocateVar(treeVarReferences.lData[index]))->toStr();
				_FString treeFS (*treeString);
				DeleteObject (treeString);
				setParameter (dataFileTreeString, &treeFS);
				setParameter (dataFileTree, 1.0);
				stashParameter (noInternalLabels,1.0,false);
			}
		}
		else
		{
			if (dsXtra)
			{
				theDF = new _DataSetFilter;
				checkPointer (theDF);
				_SimpleList e1;
				theDF->SetFilter (dsXtra,1,sp->speciesIndex,e1,false);
				cleanUp = true;
			}
			
			else
			{
				if (dataWrapper)
					theDF = dataWrapper;
				else
				{
					theDF = new _DataSetFilter;
					_DataSet* ds = (_DataSet*)dataSetList(dataSetID);
					checkPointer (theDF);
					_SimpleList e1;
					theDF->SetFilter (ds,1,sp->speciesIndex,e1,false);
					cleanUp = true;
				}
			}
		}
		_String		fileName,prompt ("Save partition as:"),
					defFileName (index>=0?(*(_String*)dataSetFilterNamesList(dfID)):
										  (*(_String*)dataSetNamesList(dataSetID))),
					listLabel("File Format:");
					
		if ((index>=0)||(!dsFile))
			formatChoice = SaveFileWithPopUp (fileName,prompt,defFileName,listLabel,dfFilterFormats);
		else
		{
			_Parameter     dsFormat;
			checkParameter (dataFilePrintFormat, dsFormat,0);
			formatChoice = dsFormat;
		}
			
		if (formatChoice>=0)
		{
			FILE*  outFile = doFileOpen ((dsFile?dsFile->sData:fileName.sData),"w");
			if (!outFile)
			{
				fileName = _String("Trying to open ")&(dsFile?*dsFile:fileName)&" for writing failed. Oops!";
				ProblemReport (fileName,(Ptr)this);
				return empty;
			}
			dfID = LocateVarByName (dataFilePrintFormat);
			if (dfID<0)
			{
				_Variable dummyVar (dataFilePrintFormat);
				dfID = LocateVarByName (dataFilePrintFormat);
			}
			_Variable* settingsVar = FetchVar(dfID);
			_Constant  dummy (formatChoice);
			_Parameter saveDefFormat = settingsVar->Value();
			
			settingsVar->SetValue (&dummy);
			_SimpleList	savedOrder;
			savedOrder.Duplicate(&theDF->theNodeMap);
			theDF->SetMap(sp->speciesIndex);
			theDF->toFileStr(outFile);
			dummy.SetValue (saveDefFormat);
			theDF->SetMap(savedOrder);
			settingsVar->SetValue (&dummy);
			
			fclose (outFile);
			if (cleanUp)
				DeleteObject (theDF);
			return (dsFile?*dsFile:fileName);
		}
		if (cleanUp)
			DeleteObject (theDF);
	}
	return empty;
}

//__________________________________________________________
void  _HYDataPanel::KillPartition (long index)
{
	if (dataSetID>=0)
	{
		if ((index<0)||(index>=dataPartitions.lLength)) return;
		if (!PurgeLFFilter(index))	
			return; 
		long dfID = dataPartitions.lData[index];
		_DataSetFilter *theDF = (_DataSetFilter*)dataSetFilterList(dfID);
		UnmarkSites (theDF->theOriginalOrder,index);
		dataPartitions.Delete(index);
		treeVarReferences.Delete(index);
		modelReferences.Delete(index);		
		partitionColors.Delete(index);
		partData.Delete (index);
		DeletePartitionRow (index);
		
		if (dataPartitions.lLength)
		{
			if (index==dataPartitions.lLength)
				index--;
			_SimpleList newSel;
			newSel<<index;
			_HYTable* pT = (_HYTable*)GetObject (5);
			pT->SetRowSelection (newSel);
		}
		CorrectSites (index);
		KillDataFilterRecord (dfID);
		UpdatePartitionOperations();
		BuildThermometer();
		BuildMarksPane();
		tainted = true;
	}
}

//__________________________________________________________
void  _HYDataPanel::CombPartition (long index)
{
	if (dataSetID>=0)
	{
		_SimpleList	res;
		bool		okcancel = false;
		
		_HYCombDialog* sd = new _HYCombDialog (&res,&okcancel, (Ptr)this);
		sd->Activate();
		while (windowObjectRefs.Find ((long)sd)>=0)
			handleGUI();
		
		if (okcancel)
		//if (ChooseCombingPattern(res))
		{
			if (res.lLength>=2)
			{
				_SimpleList		combed, steps;
				long		    k,m,n;
				
				for (k=0;k<res.lLength;k++)
					if (res.lData[k])
					{
						steps<<k;
					}
				if (index>=0)
				{
					if (!PurgeLFFilter (index))
						return;
					
					long dfID = dataPartitions.lData[index];
					_DataSetFilter *theDF = (_DataSetFilter*)dataSetFilterList(dfID);
					theDF->theOriginalOrder.Sort();
					k = 0;	
					while (k<theDF->theOriginalOrder.lLength)
					{
						m = 0;
						n = k+steps.lData[m];
						while (n<theDF->theOriginalOrder.lLength)
						{
							combed<<theDF->theOriginalOrder.lData[n];
							m++;
							if (m==steps.lLength)
								break;
							n = k+steps.lData[m];
						}
						
						k+=res.lLength;
					}
					UnmarkSites (theDF->theOriginalOrder,index);
					MarkSites (combed,index);
					steps.Clear();
					if (dataWrapper)
						steps.Duplicate (&dataWrapper->theNodeMap);
					theDF->SetFilter((_DataSet*)dataSetList(dataSetID),1,steps,combed,false);
				}
				else
				{
					_HYSequencePane* sp = (_HYSequencePane*)GetObject(0);
					if (sp->selection.lLength>=3)
					{
						k = 0;	
						while (k<sp->selection.lLength)
						{
							m = 0;
							n = k+steps.lData[m];
							while (n<sp->selection.lLength)
							{
								combed<<sp->selection.lData[n];
								m++;
								if (m==steps.lLength)
									break;
								n = k+steps.lData[m];
							}
							k+=res.lLength;
						}
						CreatePartition(combed);
					}
					
				
				}
				BuildThermometer();
				BuildMarksPane();
			}			
		}
		_HYButtonBar*   	bb = (_HYButtonBar*)GetObject (2);
		bb->_UnpushButton();
	}
}

//__________________________________________________________
void  _HYDataPanel::SplitPartition (long index, long split)
{
	if (dataSetID>=0)
	{
		if ((index<0)||(index>=dataPartitions.lLength)) 
			return;
		
		if ((partitionColors.lLength==255)&&(siteAssignments.Find(-2)>=0))
		{
			_String errMsg = _String("Couldn't perform partition splitting. ")&overlapWarning;
			ProblemReport (errMsg,(Ptr)this);
			return;
		}
		
		long 			 k,
						 m,
						 geneticCode = 0;
						 
		bool		 	 direction = true;
		char		 	 readFrame = 0;
		
		LongToPartData (partData.lData[index],readFrame,direction,geneticCode);
		if (!PurgeLFFilter(index))
			return;
			
		_DataSetFilter *theDF = (_DataSetFilter*)dataSetFilterList(dataPartitions(index));
		_SimpleList	  smallerPartition, secondPartition, emptyList;
		for (k=0; k<theDF->theOriginalOrder.lLength;k++)
		{
			m = theDF->theOriginalOrder.lData[k];
			if (m<split)
				smallerPartition<<m;
			else
				secondPartition<<m;
		}
		if (secondPartition.lLength%theDF->GetUnitLength())
		{
			_String warnMsg (
			"Your splitting position isn't valid, because the resulting partitions will NOT contain a whole number\
of evolutionary units (di-nucs or codons). Move your splitting position or change the datatype to nucleotide.");
			ProblemReport (warnMsg,(Ptr)this);
			_HYButtonBar*   	bb = (_HYButtonBar*)GetObject (2);
			bb->_UnpushButton();
			return;	
		}
		UnmarkSites (theDF->theOriginalOrder,index);
		MarkSites (smallerPartition,index);
		if (dataWrapper)
			emptyList.Duplicate (&dataWrapper->theNodeMap);
		theDF->SetFilter((_DataSet*)dataSetList(dataSetID),theDF->GetUnitLength(),emptyList,smallerPartition,false);
		CreatePartition(secondPartition, theDF->GetUnitLength());
		
		modelReferences.lData[modelReferences.lLength-1] = modelReferences.lData[index];
		partData.lData[partData.lLength-1] = partData.lData[index];

		if (theDF->GetUnitLength()==3)
		{		
			_DataSetFilter *DF2 = (_DataSetFilter*)dataSetFilterList (dataPartitions.lData[partData.lLength-1]);
			SetCodonExclusions (DF2,partData.lLength-1,geneticCode);
			SetCodonExclusions (theDF,index,geneticCode);
			theDF->SetDimensions();			
			DF2->SetDimensions();			
		}

		_List mOptions;
		GenerateModelList (mOptions,partData.lLength-1);
		RefreshPartRow(mOptions,partData.lLength-1,true);
	}
}

//__________________________________________________________
void  _HYDataPanel::JoinPartitions (long i, long j)
{
	if (dataSetID>=0)
	{
		long k,m;
		_DataSetFilter *theDF = (_DataSetFilter*)dataSetFilterList(dataPartitions(i)),
					   *DF2 = (_DataSetFilter*)dataSetFilterList(dataPartitions(j));
		bool		   updateTransl = false;
		if (partData.lData[i]!=partData.lData[j])
		{
			_String warnMsg ("The codon paritions that you are about to join have different properties (reading frame and/or directions and/or genetic codes). The resulting partition will inherit the settings of ");
			warnMsg = warnMsg & *(_String*)dataSetFilterNamesList (i) &'.';
			if (!ProceedPrompt (warnMsg,(Ptr)this))
				return;
			_HYButtonBar*   	bb = (_HYButtonBar*)GetObject (2);
			bb->_UnpushButton();	
			updateTransl = (HY_DATAPANEL_TRANSLATION&addedLines);
		}
		
		if (!(PurgeLFFilter(i)&&PurgeLFFilter (j)))
			return;
				
		long			 geneticCode = 0;
		bool		 	 direction = true;
		char		 	 readFrame = 0;
		
		LongToPartData (partData.lData[i],readFrame,direction,geneticCode);
		
		_SimpleList	  jointPartition, emptyList, colorSites;
		jointPartition.Duplicate (&theDF->theOriginalOrder);
		for (k=0; k<DF2->theOriginalOrder.lLength;k++)
		{
			m = DF2->theOriginalOrder.lData[k];
			if (siteAssignments.lData[m]>=0)
			{
				jointPartition<<m;
				siteAssignments.lData[m]=i;
				overlaps.lData[m]=0;
			}
			else
			{
				if (theDF->theOriginalOrder.Find(m)<0)
				{
					jointPartition<<m;
					colorSites<<m;
				}
			}
		}
		UnmarkSites (DF2->theOriginalOrder, j);
		MarkSites (colorSites,i);
		if (dataWrapper)
			emptyList.Duplicate (&dataWrapper->theNodeMap);
		theDF->SetFilter((_DataSet*)dataSetList(dataSetID),theDF->GetUnitLength(),emptyList,jointPartition,false);
		treeVarReferences.Delete(j);
		modelReferences.Delete(j);		
		partitionColors.Delete(j);
		partData.Delete (j);
		DeletePartitionRow (j);
		KillDataFilterRecord (dataPartitions(j));
		CorrectSites (j);
		dataPartitions.Delete(j);
		UpdatePartitionOperations();
		
		if (theDF->GetUnitLength()==3)
		{		
			SetCodonExclusions (theDF,i,geneticCode);
			theDF->SetDimensions();			
		}
		
		if (updateTransl)
		{
			_String dummy;
			UpdateTranslationString (dummy, translatedSequence,true);
		}
		BuildThermometer();
		BuildMarksPane();
	}	
}

//__________________________________________________________
void  _HYDataPanel::InterleavePartitions (long i, long j, bool dir1, bool dir2)
{
	if (dataSetID>=0)
	{
		long k,m;

		if (!(PurgeLFFilter(i)&&PurgeLFFilter(j)))
			return;

		_DataSetFilter *theDF = (_DataSetFilter*)dataSetFilterList(dataPartitions(i)),
					   *DF2 = (_DataSetFilter*)dataSetFilterList(dataPartitions(j));
				
		_SimpleList	  jointPartition, emptyList, colorSites;
		if (dir1&&dir2)
		// both reverse
		{
			m = DF2->theOriginalOrder.lLength-1;
			for (k=theDF->theOriginalOrder.lLength-1; (k>=0)&&(m>=0) ; k--, m--)
			{			
				jointPartition << theDF->theOriginalOrder.lData[k];
				jointPartition << DF2->theOriginalOrder.lData[m];
			}
			while (m>0)
				jointPartition << DF2->theOriginalOrder.lData[m--];
			while (k>0)
				jointPartition << theDF->theOriginalOrder.lData[k--];
				
		}
		else
		if (dir1)
		// first reverse
		{
			m = 0;
			for (k=theDF->theOriginalOrder.lLength-1; (k>=0)&&(m<DF2->theOriginalOrder.lLength) ; k--, m++)
			{			
				jointPartition << theDF->theOriginalOrder.lData[k];
				jointPartition << DF2->theOriginalOrder.lData[m];
			}
			while (m<DF2->theOriginalOrder.lLength)
				jointPartition << DF2->theOriginalOrder.lData[m++];
			while (k>0)
				jointPartition << theDF->theOriginalOrder.lData[k--];
				
		}		
		else
		if (dir2)
		{
			m = DF2->theOriginalOrder.lLength-1;
			for (k=0; (k<theDF->theOriginalOrder.lLength)&&(m>=0) ; k++, m--)
			{			
				jointPartition << theDF->theOriginalOrder.lData[k];
				jointPartition << DF2->theOriginalOrder.lData[m];
			}
			while (m>0)
				jointPartition << DF2->theOriginalOrder.lData[m--];
			while (k<theDF->theOriginalOrder.lLength)
				jointPartition << theDF->theOriginalOrder.lData[k++];
				
		}		
		else
		// both normal
		{
			m = 0;
			for (k=0; (k<theDF->theOriginalOrder.lLength)&&(m<DF2->theOriginalOrder.lLength) ; k++, m++)
			{			
				jointPartition << theDF->theOriginalOrder.lData[k];
				jointPartition << DF2->theOriginalOrder.lData[m];
			}
			while (m<DF2->theOriginalOrder.lLength)
				jointPartition << DF2->theOriginalOrder.lData[m++];
			while (k<theDF->theOriginalOrder.lLength)
				jointPartition << theDF->theOriginalOrder.lData[k++];				
		}		
		
		UnmarkSites (DF2->theOriginalOrder, j);
		MarkSites (DF2->theOriginalOrder,i);
		if (dataWrapper)
			emptyList.Duplicate (&dataWrapper->theNodeMap);
		theDF->SetFilter((_DataSet*)dataSetList(dataSetID),theDF->GetUnitLength(),emptyList,jointPartition,false);
		treeVarReferences.Delete(j);
		modelReferences.Delete(j);		
		partitionColors.Delete(j);
		partData.Delete (j);
		DeletePartitionRow (j);
		KillDataFilterRecord (dataPartitions(j));
		CorrectSites (j);
		dataPartitions.Delete(j);
		UpdatePartitionOperations();
		BuildThermometer();
		BuildMarksPane();
	}	
}

//__________________________________________________________
void  _HYDataPanel::SubtractPartitions (long i, long j)
{
	if (dataSetID>=0)
	{
		long k,m;
		_DataSetFilter *theDF = (_DataSetFilter*)dataSetFilterList(dataPartitions(i)),
					   *DF2   = (_DataSetFilter*)dataSetFilterList(dataPartitions(j));
					   
					   
		_SimpleList	  	splitPartition, 
						emptyList, 
						uncolorSites;
						
		for (k=0; k<theDF->theOriginalOrder.lLength;k++)
		{
			m = theDF->theOriginalOrder.lData[k];
			if (siteAssignments.lData[m]<0)
			{
				if (DF2->theOriginalOrder.Find(m)<0)
					splitPartition << m;
				else
					uncolorSites<<m;
			}
			else
				splitPartition << m;
		}
		UnmarkSites (uncolorSites, i);
		
		if (!PurgeLFFilter(i))
			return;
		
		if (splitPartition.lLength)
		{
			long			 geneticCode = 0;
			bool		 	 direction = true;
			char		 	 readFrame = 0;
			LongToPartData (partData.lData[i],readFrame,direction,geneticCode);
			
			if (dataWrapper)
				emptyList.Duplicate (&dataWrapper->theNodeMap);
				
			theDF->SetFilter((_DataSet*)dataSetList(dataSetID),theDF->GetUnitLength(),emptyList,splitPartition,false);
			if (theDF->GetUnitLength()==3)
			{		
				SetCodonExclusions (theDF,i,geneticCode);
				theDF->SetDimensions();			
			}
		}
		else
		{
			treeVarReferences.Delete(i);
			modelReferences.Delete(i);		
			partitionColors.Delete(i);
			partData.Delete (i);
			DeletePartitionRow (i);
			KillDataFilterRecord (dataPartitions(i));
			CorrectSites (i);
			dataPartitions.Delete(i);
		}
		UpdatePartitionOperations();
		BuildThermometer();
		BuildMarksPane();
	}	
}
//__________________________________________________________
void _HYDataPanel::UpdatePartitionOperations (void)
{
	_HYButtonBar* bb = (_HYButtonBar*)GetObject (2);
	_HYTable	* pl = (_HYTable*)GetObject (5);
	_SimpleList	  selp,
				  toggles;
				  
	pl->GetRowSelection (selp);
		
	bb->EnableButton(0,false);
	bb->EnableButton(1,false);
	bb->EnableButton(2,false);
	bb->EnableButton(4,false);
	bb->EnableButton(6,false);
	bb->EnableButton(7,false);

	bb->EnableButton(3,selp.lLength);

	if (selp.lLength == 1)
	{
		toggles<<1;
		bb->EnableButton(5,true);
		bb->EnableButton(0,CanSplit());
		bb->EnableButton(4,CanComb(selp.lData[0]));
		bb->EnableButton(7,true);
		
	}
	else
	{
		toggles<<0;
		bb->EnableButton(5,false);
		if (selp.lLength >= 2)
		{
			bool canJoinAll = true;
			for (long pi=1; pi<selp.lLength;pi++)
				if (!CanJoin (selp.lData[pi-1],selp.lData[pi]))
				{
					canJoinAll = false;
					break;
				}
				
			bb->EnableButton (1,canJoinAll);

		}
		if (selp.lLength == 2)
		{
			if (CanSubtract (selp.lData[0],selp.lData[1]))
				bb->EnableButton (2,true);
			if (CanInterleave (selp.lData[0],selp.lData[1]))
				bb->EnableButton (6,true);
		}
		else
			if (selp.lLength == 0)
				UpdateSelDepPartitionOperations();
	}
	_UpdatePartitionOperations(&toggles);
}

//__________________________________________________________
void _HYDataPanel::UpdateSelDepPartitionOperations (void)
{
	_HYTable*   		dl = (_HYTable*)		GetObject (5);
	_HYButtonBar*   	bb = (_HYButtonBar*)	GetObject (2);
	_HYSequencePane* 	sp = (_HYSequencePane*)	GetObject (0);
	
	_SimpleList 		selp;	
	dl->GetRowSelection (selp);
	
	if (selp.lLength!=1)
		bb->EnableButton(4,(sp->selection.lLength>=3));
		
	if (selp.lLength == 1)
		bb->EnableButton(0,CanSplit());
	else
		bb->EnableButton(0,false);
		
	statusBar.Trim (0,statusBar.Find (selectionStatPrefix)+selectionStatPrefix.sLength-1);
	_String newSB (statusBar);
	
	if (sp->selection.lLength==0)
	{
		if (sp->vselection.lLength)
		{
			bb->EnableButton(1,CanJoinSpecies(sp->vselection));
			bb->EnableButton (3,true);
		}
		else	
			bb->EnableButton (3,false);
		newSB = newSB &"empty";
	}
	else
	{
		long hi = sp->selection.lData[sp->selection.lLength-1]+1,
			 lo = sp->selection.lData[0]+1;
		
		if (hi-lo+1>sp->selection.lLength)
			newSB = newSB & " part of ";
			
		newSB = newSB & _String (lo) & '-' & _String (hi);
	}
	if (!statusBar.Equal(&newSB))
	{
		SetStatusBar (newSB);
		if (forceUpdateForScrolling)
		{
			_PaintStatusBar();
			_PaintLFStatus ();
		}
	}
}

//__________________________________________________________
void _HYDataPanel::JoinSpeciesDisplay (void)
{
	_HYSequencePane* 	sp = (_HYSequencePane*)GetObject(0);
	_HYButtonBar*   	bb = (_HYButtonBar*)GetObject (2);
	_SimpleList  newSpIndex;
	long		 k,m=1,p=sp->vselection.lData[0];
	
	for (k=0; k<=sp->vselection.lData[0];k++)
		newSpIndex<<sp->speciesIndex.lData[k];
	for (k=1; k<sp->vselection.lLength;k++)
		newSpIndex<<sp->speciesIndex.lData[sp->vselection.lData[k]];
	for (k=sp->vselection.lData[0]+1;k<sp->speciesIndex.lLength;k++)
	{
		if (k==sp->vselection.lData[m])
		{
			m++;
			if (m==sp->vselection.lLength)
				break;
		}
		else
			newSpIndex<<sp->speciesIndex.lData[k];
	}
	for (k++;k<sp->speciesIndex.lLength;k++)
	{
		newSpIndex<<sp->speciesIndex.lData[k];
	}
	sp->speciesIndex.Clear();
	sp->speciesIndex.Duplicate (&newSpIndex);
	m=sp->vselection.lLength;
	sp->vselection.Clear();
	for (k=p; k<p+m; k++)
	{
		sp->vselection<<k;
	}
	sp->BuildPane();
	sp->_MarkForUpdate();
	bb->EnableButton (1,false);
}

//__________________________________________________________
void _HYDataPanel::OmitSelectedSpecies (void)
{
	_HYSequencePane* 	sp = (_HYSequencePane*)GetObject(0);
	_HYButtonBar*   	bb = (_HYButtonBar*)GetObject (2);
	_HYTable	*   	pl = (_HYTable*)GetObject (5);
	
	long				k;
	
	if (sp->vselection.lLength==0) return;
	
	if (sp->speciesIndex.lLength-sp->vselection.lLength<2)
	{
		ProblemReport(cantOmitWarning,(Ptr)this);
		bb->_UnpushButton ();
		return;
	}
	
	if (!warnOmit)
	{
		bool prcd = ProceedPromptWithCheck (seqOmitWarning,donotWarnAgain,warnOmit, (Ptr)this);
		bb->_UnpushButton ();
		if (!prcd) return;
	}
	
	if (!PurgeLF())
		return;
		
	// remove tree references from the table
	bool	update = false;
	for (k=0; k<treeVarReferences.lLength; k++)
		if (treeVarReferences.lData[k]>=0)
		{
			treeVarReferences.lData[k]=-1;			
			long idx = k*pl->horizontalSpaces.lLength+DF_TREE_COLUMN;
			pl->SetCellData (&none, k,DF_TREE_COLUMN, pl->cellTypes.lData[idx],true);
			update = true;
		}
	
	if (update)
		pl->_MarkColumnForUpdate (2);

	_SimpleList			omit;
	for (k=0; k<sp->vselection.lLength; k++)
	{
		//omittedSeqs.InsertElement((BaseRef)(sp->speciesIndex.lData[sp->vselection.lData[k]]),-1,false,false);
		omittedSeqs<<sp->speciesIndex.lData[sp->vselection.lData[k]];
		omit<<sp->speciesIndex.lData[sp->vselection.lData[k]];
	}
	_OmitSelectedSpecies(omit);
	for (k=sp->vselection.lLength-1;k>=0;k--)
		sp->speciesIndex.Delete (sp->vselection.lData[k]);
		
	UpdateDataWrapper();
	sp->vselection.Clear();
	sp->BuildPane();
	sp->_MarkForUpdate();
	bb->EnableButton (3,false);
	_String newCS;
	UpdateConsensusSequence (newCS,true);
	SetWindowRectangle (top,left,bottom,right);
}

//__________________________________________________________
void _HYDataPanel::RestoreOmittedSequence (long index)
{
	_HYSequencePane* 	sp = (_HYSequencePane*)GetObject(0);
	
	if (!PurgeLF())
		return;
	
	_RestoreOmittedSequence(index);
	
	if (index>=0)
	{
		sp->speciesIndex << omittedSeqs.lData[index];
		omittedSeqs.Delete(index);
	}
	else
	{
		for (long k=0; k<omittedSeqs.lLength; k++)
			sp->speciesIndex << omittedSeqs.lData[k];
		omittedSeqs.Clear();
	}

	UpdateDataWrapper();
	sp->vselection.Clear();
	sp->BuildPane();
	sp->_MarkForUpdate();
	_String newCS;
	UpdateConsensusSequence (newCS,true);
	SetWindowRectangle (top,left,bottom,right);
}
//__________________________________________________________
void _HYDataPanel::UpdateDataWrapper()
{
	_DataSet 			*ds = (_DataSet*)dataSetList (dataSetID);
	_HYSequencePane* 	sp  = (_HYSequencePane*)GetObject(0);
	if (!dataWrapper)
		dataWrapper = new _DataSetFilter();

	if (!omittedSeqs.lLength)
	{
		DeleteObject (dataWrapper);
		dataWrapper = nil;
	}
	else
	{
		_SimpleList dummy;
		dataWrapper->SetFilter (ds,1,sp->speciesIndex,dummy,false);
	}
	
	for (long k=0; k<dataPartitions.lLength;k++)
	{
		_SimpleList vp;
		_DataSetFilter* dsf = (_DataSetFilter*)dataSetFilterList (dataPartitions.lData[k]);
		// check tree remaps here
		vp.Duplicate(&dsf->theOriginalOrder);
		dsf->SetFilter (ds,dsf->GetUnitLength(),sp->speciesIndex,vp,false);
		if (dsf->GetUnitLength() == 3)
		{
			long			 geneticCode = 0;
			bool		 	 direction = true;
			char		 	 readFrame = 0;
			LongToPartData 	(partData.lData[k],readFrame,direction,geneticCode);
			SetCodonExclusions (dsf,k,geneticCode);
			dsf->SetDimensions();	
		}
	}
	if (addedLines&HY_DATAPANEL_CONSENSUS)
	{
		_String newCS;
		UpdateConsensusSequence (newCS,true);
	}
	GenerateStatusLine();
}

//__________________________________________________________

bool _HYDataPanel::CanSplit (void)
{
	_HYTable*   		dl = (_HYTable*)GetObject (5);
	_HYSequencePane* 	sp = (_HYSequencePane*)GetObject(0);
		
	if (sp->selection.lLength!=1) 
		return false;

	_DataSetFilter*	  thisFilter = (_DataSetFilter*)dataSetFilterList (dataPartitions(dl->GetFirstRowSelection()));
	
	long k = sp->selection.lData[0];
	k = thisFilter->theOriginalOrder.Find(k);
	
	if (k>0)
		return true;
		
	return false;
}

//__________________________________________________________

bool _HYDataPanel::CanJoin (long i, long j)
{
	_DataSetFilter* fi = (_DataSetFilter*)dataSetFilterList (dataPartitions(i)),
				  * fj = (_DataSetFilter*)dataSetFilterList (dataPartitions(j));
	return fi->GetUnitLength()==fj->GetUnitLength();
}

//__________________________________________________________

bool _HYDataPanel::CanInterleave (long i, long j)
{
	_DataSetFilter* fi = (_DataSetFilter*)dataSetFilterList (dataPartitions(i)),
				  * fj = (_DataSetFilter*)dataSetFilterList (dataPartitions(j));
	if ((fi->GetUnitLength()==1)&&(fj->GetUnitLength()==1))
	{
		return (!CanSubtract(i,j));
	}
	return false;
}

//__________________________________________________________

bool _HYDataPanel::CanSubtract (long i, long j)
{
	_DataSetFilter* fi = (_DataSetFilter*)dataSetFilterList (dataPartitions(i)),
				  * fj = (_DataSetFilter*)dataSetFilterList (dataPartitions(j));
	if (fi->GetUnitLength()==fj->GetUnitLength())
	{
		for (long k=0; k<fi->theOriginalOrder.lLength; k++)
		{
			long m = fi->theOriginalOrder.lData[k];
			if (siteAssignments.lData[m]<0)
			{
				long p  = overlaps.lData[m],
					 p1 = p&0xFF,
					 p2 = (p&0xFF00)>>8,
					 p3 = (p&0xFF0000)>>16,
					 p4 = (p&0xFF000000)>>24;
				if (p3==0)
				{
					if (((i==p1)||(i==p2))&&((j==p1)||(j==p2)))
					{
						return true;
					}
					continue;
				}
				if (p4==0)
				{
					if (((i==p1)||(i==p2)||(i==p3))&&((j==p1)||(j==p2)||(j==p3)))
					{
						return true;
					}
					continue;
				}
				if (((i==p1)||(i==p2)||(i==p3)||(i==p4))&&((j==p1)||(j==p2)||(j==p3)||(j==p4)))
				{
					return true;
				}
			}
		}
	}
	return false;
}

//__________________________________________________________

bool _HYDataPanel::CanJoinSpecies (_SimpleList& vsel)
{
	if (vsel.lLength > 1)
	{
		for (long k=1; k<vsel.lLength;k++)
			if (vsel.lData[k]-vsel.lData[k-1]>1) return true;
	}
	return false;
}

//__________________________________________________________

bool _HYDataPanel::CanComb (long i)
{
	_DataSetFilter* fi = (_DataSetFilter*)dataSetFilterList (dataPartitions(i));
	return fi->GetUnitLength()==1;
}

//__________________________________________________________
void _HYDataPanel::SelectPartition (void)
{
	_HYTable*   	 dl = (_HYTable*)GetObject (5);
	_HYSequencePane* sp = (_HYSequencePane*)GetObject(0);
	
	long filterIndex = dataPartitions(dl->GetFirstRowSelection());
	
	_DataSetFilter*    dsf = (_DataSetFilter*)dataSetFilterList (filterIndex);
	
	sp->SelectRange (dsf->theOriginalOrder);
	sp->ScrollToHSelection();
	UpdateSelDepPartitionOperations();
	_UpdateSelectionChoices (true);
}

//__________________________________________________________
void _HYDataPanel::InvertSelection(void)
{
	_HYSequencePane* sp = (_HYSequencePane*)GetObject(0);
	
	_SimpleList		newSel,
			       *selList = nil;
			      
	long			shift = 0, 
					k = 0,
					upTo;
					
	if (sp->vselection.lLength)
	{
		selList = &sp->vselection;
		upTo    = sp->speciesIndex.lLength;
	}
	else
	{
		selList = &sp->selection;
		upTo	= sp->columnStrings.lLength;
		if (selList->lLength == 0 &&  ((_HYTable*)GetObject (5))->GetFirstRowSelection() >=0)
		{
			SelectPartition ();
			InvertSelection ();
			return;
		}
		
	}
	
	for (; (k<upTo)&&(shift<selList->lLength); k++)
	{
		if (k==selList->lData[shift])
			shift++;
		else
			newSel<<k;
	}
	
	for (;k<upTo;k++)
		newSel<<k;
		
	sp->SelectRange (newSel,sp->vselection.lLength);
	UpdateSelDepPartitionOperations();
	_UpdateSelectionChoices (true);
}

//__________________________________________________________
void _HYDataPanel::PartitionPropsMenu (void)
{
	_HYTable*   dl = (_HYTable*)GetObject (5);
	
	EditPartitionProperties (dl->GetFirstRowSelection());
}
//__________________________________________________________

void _HYDataPanel::ProcessContextualPopUpMain (long l, long t)
{
	_HYSequencePane*  sp = (_HYSequencePane*)GetObject (0);
	_List			  menuItems;
	
	_String			  buffer;
	
	if (sp->selection.lLength)
	{
		if (dataType&HY_DATAPANEL_NUCDATA)
			menuItems&& & contextString1;
		else
			if (dataType&HY_DATAPANEL_PROTDATA)
				menuItems&& & contextString2;
			else
				if (dataType&HY_DATAPANEL_BINARYDATA)
					menuItems&& & contextString1_5;
		
		menuItems&& & contextString13;
		if (dataType&HY_DATAPANEL_PROTDATA)
			menuItems&& & contextString17;
		
		menuItems&& & contextString15;
		
		if (sp->selection.lLength==1)
		{
			if (menuItems.lLength)
			{
				buffer = menuSeparator;
				menuItems&& & buffer;
			}
			menuItems && & contextString4;
			menuItems && & contextString14;
		}
		
		buffer = HandlePullDown (menuItems,l,t,0);
		
		if (buffer.Equal (&contextString4))
		{
			// select all sites like this one
			_SimpleList		matches;
			if (dataWrapper)
				dataWrapper->FindAllSitesLikeThisOne(sp->selection.lData[0],matches);
			else
			{
				_DataSet* ds = (_DataSet*)dataSetList (dataSetID);
				ds->FindAllSitesLikeThisOne(sp->selection.lData[0],matches);
			}
			_String fres (32L, true);
			if (matches.lLength>1)
			{
				sp->selection.Clear();
				sp->selection.Duplicate (&matches);
				sp->BuildPane();
				sp->_MarkForUpdate();
				UpdateSelDepPartitionOperations();
				fres << "Found ";
				fres << (long)matches.lLength;
				fres << " sites like the one selected\n";
			}
			else
				BufferToConsole ("Selected site is unique\n");
		}	
		else
			if (buffer.Equal (&contextString1) || buffer.Equal (&contextString2) || buffer.Equal (&contextString1_5) )
			{
				char	  bufferS [512];
				_DataSet* ds = (_DataSet*)dataSetList (dataSetID);
				_Matrix*  freqs = ds->HarvestFrequencies(1,1,false,sp->speciesIndex,sp->selection);
				long	  totCount = sp->speciesIndex.lLength*sp->selection.lLength;
				if (buffer.Equal (&contextString1))
				{
					sprintf (bufferS,"\nNucleotide counts and frequencies (%d sites)\nA:\t%g\t%g\nC:\t%g\t%g\nG:\t%g\t%g\nT:\t%g\t%g\n",sp->selection.lLength,
							(totCount*(*freqs)[0]),(*freqs)[0],(totCount*(*freqs)[1]),(*freqs)[1],(totCount*(*freqs)[2]),(*freqs)[2],(totCount*(*freqs)[3]),(*freqs)[3]);
					BufferToConsole (bufferS);
				}
				else
					if (buffer.Equal (&contextString2))
					{
						sprintf (bufferS,"\nAminoacid counts and frequencies (%d sites)",sp->selection.lLength);
						BufferToConsole (bufferS);
						for (long k=0; k<aminoAcidOneCharCodes.sLength; k++)
						{
							sprintf (bufferS,"\n%c:\t%g\t%g",aminoAcidOneCharCodes.sData[k],(totCount*(*freqs)[k]),(*freqs)[k]);
							BufferToConsole (bufferS);
						}
						NLToConsole ();
					}
					else
					{
						sprintf (bufferS,"\nBinary counts and frequencies (%d sites)\n0:\t%g\t%g\n1:\t%g\t%g\n",sp->selection.lLength,
								(totCount*(*freqs)[0]),(*freqs)[0],(totCount*(*freqs)[1]),(*freqs)[1]);
						BufferToConsole (bufferS);
					}
				DeleteObject(freqs);
			}
			else
			{
				if (buffer.Equal (&contextString13) || buffer.Equal (&contextString17))
				{
					_String *toCopy;
					if (buffer.Equal (&contextString13))
						toCopy = (_String*)sp->selection.ListToPartitionString();
					else
					{
						_SimpleList dup;
						dup.RequestSpace (sp->selection.lLength*3);
						for (long k=0; k<sp->selection.lLength;k++)
						{
							dup << 3*sp->selection.lData[k];
							dup << 3*sp->selection.lData[k]+1;
							dup << 3*sp->selection.lData[k]+2;
						}
						toCopy = (_String*)dup.ListToPartitionString();
					}

					PlaceStringInClipboard (*toCopy, (Ptr)this);
					DeleteObject (toCopy);
				}
				else
				{
					if (buffer.Equal (&contextString14))
					{
						_DataSet* 		ds    = (_DataSet*)dataSetList (dataSetID);
						_String 		aSite (*(_String*) (ds->_List::operator ()) (ds->GetTheMap().lData[sp->selection.lData[0]]));
						
						_SimpleList	    remap;
						DeleteObject	(aSite.Sort(&remap));
						sp->SetSequenceOrder (remap);
						
					}
					else
						if (buffer.Equal (&contextString15))
						{
							_String	 nhc ("Please choose a new highlight color");
							_HYColor newC = SelectAColor (sp->highlightColor,nhc);
							if (!(newC==sp->highlightColor))
								sp->SetHighliteColor (newC);
						}
				}
			}
	}
	else
	{
		_SimpleList		treeIDs;
		for (long l2=0; l2<treeVarReferences.lLength; l2++)
		{
			long t = treeVarReferences.lData[l2];
			if (t>=0)
			{
				buffer = _String("Order as in '")&*LocateVar(t)->GetName()&"'";
				menuItems && & buffer;
				treeIDs << l2;
			}
		}

		if (sp->vselection.lLength==1)
		{
			if (addedLines&HY_DATAPANEL_CONSENSUS)
			{
				if (menuItems.lLength)
					menuItems && & menuSeparator;
				menuItems && & contextString5;
				menuItems && & contextString10;
			}
			
			if ((addedLines&HY_DATAPANEL_REFERENCE)&&(genCodeID<0))
			{
				if (menuItems.lLength)
					menuItems && & menuSeparator;
				menuItems && & contextString8;
			}
		}
		
		if (sp->vselection.lLength)
		{
			menuItems && & copySeqsToClip;
			menuItems && & contextString11;
			menuItems && & contextString12;
			menuItems && & contextString16;
		}

		if (menuItems.lLength)
		{
			buffer = HandlePullDown (menuItems,l,t,0);
			if (buffer.sLength)
				if (buffer.Equal(&contextString5)||buffer.Equal(&contextString8))
				{
					_SimpleList noMatch;
					t = sp->speciesIndex.lData[sp->vselection[0]];
					long        sl = 0;
					if (buffer.Equal (&contextString8))
						sl = ((_String*)statusData(l))->sLength-1;
					for (l=0; l<sp->columnStrings.lLength;l++)
					{
						if (((_String*)sp->columnStrings(l))->sData[t]!=((_String*)statusData(l))->sData[sl])
							noMatch << l;
					}
					if (noMatch.lLength==0)
						BufferToConsole ("A perfect match!\n");
					else
					{
						sp->vselection.Clear();
						sp->selection.Clear();
						sp->selection.Duplicate (&noMatch);
						sp->BuildPane();
						sp->_MarkForUpdate();
						UpdateSelDepPartitionOperations();
						char buffer [256];
						sprintf (buffer,"%d (%g %%) mismatches found\n",noMatch.lLength, 
																noMatch.lLength*100./(_Parameter)sp->columnStrings.lLength);
						BufferToConsole (buffer);
					}
				}
				else
				{
					if (buffer.Equal(&contextString10))
					{
						_String consString (sp->columnStrings.lLength,true);
						for (long kl =0; kl<sp->columnStrings.lLength;kl++)
							consString << ((_String*)statusData(kl))->sData[0];
						consString.Finalize();
						PlaceStringInClipboard (consString, GetOSWindowData());
					}
					else
					{
						if (buffer.Equal (&copySeqsToClip))
						{
							_String theSeqs (sp->vselection.lLength*sp->columnStrings.lLength,true);
							
							for (long sc = 0; sc < sp->vselection.lLength; sc ++)
							{
								long seqIndex = sp->speciesIndex.lData[sp->vselection.lData[sc]];
								
								theSeqs << '>';
								theSeqs << *(_String*)sp->rowHeaders (seqIndex);
								theSeqs << '\n';
								for (long kl =0; kl<sp->columnStrings.lLength;kl++)
									theSeqs << ((_String*)sp->columnStrings(kl))->sData[seqIndex];
								
								theSeqs << '\n';
							}
							theSeqs.Finalize();
							
							PlaceStringInClipboard (theSeqs, GetOSWindowData());						
						}
						else
						{
							if (buffer.Equal (&contextString11))
							{
								_DataSet* newDS 	= new _DataSet (),
										* thisDS	= (_DataSet*)dataSetList (dataSetID);
										
								checkPointer (newDS);
								newDS->SetTranslationTable (thisDS->GetTT());
								for (long sc = 0; sc < sp->vselection.lLength; sc ++)
								{
									long seqIndex = sp->speciesIndex.lData[sp->vselection.lData[sc]];
									
									newDS->AddName (*(_String*)sp->rowHeaders (seqIndex)); 
									if (sc == 0)
									{
										for (long kl =0; kl<sp->columnStrings.lLength;kl++)
											newDS->AddSite (((_String*)sp->columnStrings(kl))->sData[seqIndex]);
									}
									else
										for (long kl =0; kl<sp->columnStrings.lLength;kl++)
											newDS->Write2Site (kl,((_String*)sp->columnStrings(kl))->sData[seqIndex]);
								}
								newDS->Finalize();
								newDS->SetNoSpecies (sp->vselection.lLength);
								_String newName (GetTitle()); 
								newName.Trim(newName.FirstSpaceIndex (0,-1,-1)+1,-1);
								AddDataSetToList (newName,newDS);
								_HYDataPanel* myTWindow = new _HYDataPanel(newName,newName);
								myTWindow->BringToFront();
							}
							else
							{
								if (buffer.Equal (&contextString12))
								{
									_String			  res,
													  prompt ("Minimum ambiguity %-age to be selected:");
													  
									if (EnterStringDialog (res,prompt,(Ptr)this))
									{
										_Parameter 	  gateValue = res.toNum();
										if (gateValue < 0.0 || gateValue > 100.0)
										{
											res = _String ("Invalid percentage: ") & res;
											ProblemReport (res,(Ptr)this);
										}
										else
										{
											_SimpleList			newVerticalSelection;
											_DataSet* thisDS	= (_DataSet*)dataSetList (dataSetID);
											_TranslationTable* myTT = thisDS->GetTT();
											
											long				charCount = myTT->LengthOfAlphabet(),
																*charSpool = new long [charCount],
																dsLength  = thisDS->NoOfColumns();
																
											checkPointer 		(charSpool);

											for (long sc = 0; sc < sp->vselection.lLength; sc ++)
											{
												long seqIndex   = sp->speciesIndex.lData[sp->vselection.lData[sc]],
													 ambigCount = 0;
													 
												for (long cc = 0; cc < dsLength; cc ++)
												{
													myTT->TokenCode ((*thisDS)(cc,seqIndex,1),charSpool);
													char countMe = 0;
													
													for (long cc2 = 0; cc2 < charCount && countMe < 2; cc2++)
														countMe += charSpool[cc2];
														
													if (countMe > 1)
														ambigCount ++;
												}
												
												if ((ambigCount*100.0)/dsLength >= gateValue)
													newVerticalSelection << sp->vselection.lData[sc];
											}
											
											sp->SelectRange (newVerticalSelection,true);
											delete charSpool;
											
										}
									}
								}
								else
								if (buffer.Equal (&contextString16))
								{
									long choice = SelectOpenWindowObjects (HY_WINDOW_KIND_TREE, "tree panel", treeSelectorPopulator, (Ptr)this);
									if (choice>=0)
									{
										_HYTreePanel * tp 		 = ((_HYTreePanel*)windowObjectRefs(choice));
										_List 		   selectedSequences;
										_String 	   treeNameString (*tp->LocateMyTreeVariable()->GetName()&'.');
										
										for (long sc = 0; sc < sp->vselection.lLength; sc ++)
										{
											_String seqName = treeNameString & *(_String*)sp->rowHeaders (sp->speciesIndex.lData[sp->vselection.lData[sc]]);
											selectedSequences && & seqName;
										}
										selectedSequences.Sort();
										tp->SelectRangeAndScroll(selectedSequences, false);
										
									}
								}
								else
								{
									t = treeIDs.lData[menuItems.Find(&buffer)];
									_DataSetFilter* fi = (_DataSetFilter*)dataSetFilterList (dataPartitions(t));
									if (((_TheTree*)LocateVar(treeVarReferences.lData[t]))->MatchLeavesToDF (treeIDs, fi, true))
										sp->SetSequenceOrder (treeIDs);
								}
							}
						}
					}
				}
		}
	}
	_UpdateSelectionChoices(IsSelectionNonEmpty());
}

//__________________________________________________________

void _HYDataPanel::ProcessContextualPopUpAux (long l, long t)
{
	_HYSequencePane*  sp = (_HYSequencePane*)GetObject (4);
	_List			  menuItems;
	
	_String			  buffer;

	if (addedLines&HY_DATAPANEL_RATECLASS)
	{
		long		 dataCount [37],
					 i,
					 accessIndex = 0,
					 m = 0;
					 
		char		 c;
					 
		if (addedLines&HY_DATAPANEL_CONSENSUS)
			accessIndex ++;
		
		for	(i=0; i<37; i++)
			dataCount [i] = 0;
		
		menuItems && & contextString6;
		menuItems && & menuSeparator;
		
		for (i=0; i<siteAssignments.lLength; i++)
		{
			c = ((_String*)sp->columnStrings(i))->sData[accessIndex];
			if (c=='*')
				dataCount [36]++;
			else
				if ((c>='0')&&(c<='9'))
					dataCount [c-'0']++;
				else
					dataCount [c-'A'+10]++;
		}
		

		for (i=0; i<10; i++)
		{
			if (dataCount[i]>0)
			{
				buffer = contextString7 & i;
				menuItems && & buffer;
				m = i;
			}
		}

		for (i=0; i<26; i++)
		{
			if (dataCount[i+10]>0)
			{
				buffer = contextString7 & (char)('A'+i);
				menuItems && & buffer;
				m = i+10;
			}
		}
		
		if (dataCount[36])
		{
			m = 36;
			buffer = contextString7 & "*";
			menuItems && & buffer;
		}
		
		menuItems && & menuSeparator;
		buffer = "Auto-partition by rate class";
		menuItems && & buffer;
		
		i = m+1;
		
		buffer = HandlePullDown (menuItems,l,t,0);
		
		m = menuItems.Find (&buffer);
		
		if (m==0)
		{
			_List		rowTitles;
			_String		entry ("Count");
			rowTitles   && & entry;
			entry = "Proportion";
			rowTitles   && & entry;
			
			entry		= _String("Frequency Data for ")&GetTitle();
			_String     tryMe (entry);
			
			m = 2;
			while (FindWindowByName (tryMe)>=0)
			{
				tryMe = entry & '[' & m &']';
				m++;
			}
			
			_Matrix		dataPoints (i+1,2,false,true);
			for (m=0; m<i; m++)
			{
				dataPoints.Store (m,0,dataCount[m]);
				dataPoints.Store (m,1,dataCount[m]/(_Parameter)siteAssignments.lLength);
			}
			
			_HYChartWindow* cw = new _HYChartWindow (tryMe,rowTitles, dataPoints, nil);
			checkPointer (cw);
			cw->SetChartType ("Bar Chart","Index","Count");
			
			cw->BringToFront();
		}
		else
		if (m>1)
		{
			if (m!=menuItems.lLength-1)
			{
				_SimpleList matches;
				c = ((_String*)menuItems (m))->sData[((_String*)menuItems (m))->sLength-1];
				for (i=0; i<siteAssignments.lLength; i++)
				{
					if (((_String*)sp->columnStrings(i))->sData[accessIndex] == c)
						matches<<i;
				}	
				if (matches.lLength)
				{
					sp = (_HYSequencePane*)components(0);
					sp->selection.Clear();
					sp->selection.Duplicate (&matches);
					sp->BuildPane();
					sp->_MarkForUpdate();
					UpdateSelDepPartitionOperations();
					char    buffer[128];
					sprintf (buffer,"Found %d sites in rate class %c\n",matches.lLength,c);
					BufferToConsole (buffer);
				}
			}
			else
			{
				_List  allMatches;
				for (m=0; m<i; m++)
					allMatches << new _SimpleList;
					
				for (i=0; i<siteAssignments.lLength; i++)
				{
					c = ((_String*)sp->columnStrings(i))->sData[accessIndex];
					_SimpleList * ss;
					if ((c<='9')&&(c>='0'))
						ss = (_SimpleList*)allMatches (c-'0');
					else
						if ((c>='A')&&(c<='Z'))
							ss = (_SimpleList*)allMatches (10+c-'A');
						else
							ss = (_SimpleList*)allMatches (36);
							
					(*ss) << i;
				}	
				
				for (m=0; m<allMatches.lLength; m++)
				{
					_String 	pName ("Class_");
					if (m<10)
						pName = pName & m;
					else
						if (m<36)
							pName = pName & (char)(m-10+'A');
						else
							pName = pName & "high";
							
					_SimpleList *ss1 = (_SimpleList*)allMatches (m);
					if (ss1->lLength)
						CreatePartition (*ss1,1,false,&pName);
				}
			}
		}
	}
	_UpdateSelectionChoices(IsSelectionNonEmpty());
}
//__________________________________________________________
_String  _HYDataPanel::TreeTopologyChange (long filterID, _String* menuChoice)
{
	long		newTreeID, 
				oldTreeID;
				
	_String		res (*menuChoice);
					
	oldTreeID = treeVarReferences.lData[filterID];
	
	if (menuChoice->Equal(&makeNewTree))
	{
		res = NewTreeWindow (dataPartitions.lData[filterID]);
		if (res.sLength)
			#ifndef USE_AVL_NAMES
			newTreeID = variableReindex.lData[LocateVarByName(res)];
			#else
			newTreeID = variableNames.GetXtra(LocateVarByName(res));			
			#endif
		else
			newTreeID = -2; // no change
	}
	else
		if (menuChoice->Equal(&none))
		{
			newTreeID = -1;
			res = none;
		}
		else
			if (menuChoice->Equal (&readTreeFile))
			{
				newTreeID = -2;
				_List	before,
						after;
						
				GenerateTreeList (before); 
				if (OpenTreeFile())
				{
					GenerateTreeList (after);
					if (after.lLength>before.lLength)
					{
						before.Sort();
						for (long k = 0; k<after.lLength; k++)
						{
							_String * thisTree = (_String*)after(k);
							if (thisTree->Equal (&menuSeparator))
								continue;
								
							if (before.BinaryFind (thisTree)<0)
							{
								#ifndef USE_AVL_NAMES
								newTreeID = variableReindex.lData[LocateVarByName(*thisTree)];
								#else
								newTreeID = variableNames.GetXtra(LocateVarByName(*thisTree));			
								#endif
								res = *thisTree;
								break;
							}	
						}
					}
				}
			}
			else
				if (menuChoice->Equal (&inferTreeStr))
					newTreeID = -3;
				else
					#ifndef USE_AVL_NAMES
					newTreeID = variableReindex.lData[LocateVarByName(*menuChoice)];
					#else
					newTreeID = variableNames.GetXtra(LocateVarByName(*menuChoice));			
					#endif
	
	// handle treeID change
	
	if (oldTreeID!=newTreeID)
	{
		if (newTreeID==-2)
		{
			newTreeID = oldTreeID;
			if (oldTreeID>=0)
				res = *LocateVar (oldTreeID)->GetName();	
			else
				if (oldTreeID == -3)
					res = inferTreeStr;
				else
					res = none;
		}
		else
		{
			if (!PurgeLFFilter(filterID)) 
				return *(_String*)((_HYTable*)GetObject(5))->GetCellData (DF_TREE_COLUMN,filterID);
			treeVarReferences[filterID]=newTreeID;
		}
	}
	return res;
}

//__________________________________________________________
void  _HYDataPanel::ModelChange (_String* menuChoice, long filterID, long newModelID)
{
	if (!PurgeLFFilter(filterID)) 
		return;
	
	long		 oldModelID;
				 				 
	oldModelID = modelReferences.lData[filterID];
	
	if (menuChoice->Equal(&none))
	{
		newModelID = 0;
		modelReferences.lData[filterID] = 0;
	}
	else
	{
		_List* 		  theList = FindModelTemplate (menuChoice);
		_SimpleList*  opt = (_SimpleList*)(*theList)(1);
		
		if (!(opt->lData[0]&HY_DATAPANEL_MODEL_MODELS))
		{
			modelReferences.lData[filterID] = ModelDataToLong (newModelID,0,0,0);
		}
		else
		{
			oldModelID = 0;
			if (opt->lData[0]&HY_DATAPANEL_MODEL_GLOBALG)
				oldModelID = 2;
			else
				if (opt->lData[0]&HY_DATAPANEL_MODEL_GLOBAL)
					oldModelID = 1;
					
			if (oldModelID<2)
				modelReferences.lData[filterID] = ModelDataToLong (newModelID,oldModelID,4,0);
			else
				modelReferences.lData[filterID] = ModelDataToLong (newModelID,2,4,4);
		}
	}
}

//__________________________________________________________
void  _HYDataPanel::ModelOptionChange (_String* oType, long filterID)
{
	if (!PurgeLFFilter(filterID)) 
			return;
	
	long		 newModelID; 
	
	for			 (newModelID=0; newModelID<2; newModelID++)
		if (parameterOption[newModelID].Equal(oType))
			break;

	int		model, 
			options, 
			freqs, 
			rates;

	LongToModelData (modelReferences.lData[filterID],model,options,freqs,rates);
	
	if (newModelID==2)
		rates = 4;
	else
		rates = 0;
	
	if (options != newModelID)
	{
		modelReferences.lData[filterID]=ModelDataToLong (model,newModelID,freqs,rates);
		//UpdatePartitionOperations();
	}
}

//__________________________________________________________
void  _HYDataPanel::ModelFreqChange (_String* oType, long filterID)
{
	if (!PurgeLFFilter(filterID)) 
			return;
	
	long		 newModelID; 
	
	for			 (newModelID=0; newModelID<4; newModelID++)
		if (freqOption[newModelID].Equal(oType))
			break;

	int		model, 
			options, 
			freqs, 
			rates;

	LongToModelData (modelReferences.lData[filterID],model,options,freqs,rates);
	
	if (freqs != newModelID)
	{
		modelReferences.lData[filterID]=ModelDataToLong (model,options,newModelID,rates);
		//UpdatePartitionOperations();
	}
}

//__________________________________________________________
void  _HYDataPanel::ModelRateClassChange (long filterID, long newModelID)
{
	int			 model, 
				 options, 
				 freqs, 
				 rates;

	LongToModelData (modelReferences.lData[filterID],model,options,freqs,rates);

	if (rates != newModelID)
	{
		modelReferences.lData[filterID]=ModelDataToLong (model,options,freqs,newModelID);
		if (lfID>=0)
		{
			if (lfID != lockedLFID)
			{
				_LikelihoodFunction* thisLF = (_LikelihoodFunction*)likeFuncList(lfID);
				if (thisLF->GetTheFilters().Find (dataPartitions.lData[filterID])<0)
					return;
					
				_CategoryVariable*   cv = thisLF->FindCategoryVar (filterID);
				if (!cv)
				{
					_String	warnMsg ("Internal error in ModelRateClassChange. Please rebuild the likelihood function.");
					PurgeLF();
					ProblemReport	(warnMsg,(Ptr)this);
				}
				else
				{
					cv->ChangeNumberOfIntervals (newModelID);
					thisLF->SetIthIndependent (0,thisLF->GetIthIndependent(0));
				}
			}
			else
				ProblemReport (lfCantKillWarning, (Ptr)this);
		}
	}
}
//__________________________________________________________
bool  _HYDataPanel::PurgeLF(bool all)
{
	if (lfID>=0)
	{
		if ((cantDeleteLF)||(lockedLFID == lfID))
		{
			ProblemReport (lfCantKillWarning, (Ptr)this);
			return false;
		}
		if ((!warnKillLF)&&all)
		{
			bool prcd = ProceedPromptWithCheck (lfKillWarning,donotWarnAgain,warnKillLF, (Ptr)this);
			if (!prcd) return false;
		}
		
		postLFKillEvent (GetID(),lfID);
		KillLFRecord(lfID);
		if (all)
		{
			lfID = -1;
			if (addedLines&HY_DATAPANEL_RATECLASS)
				AdjustStatusLine (1);
			_HYButtonBar*   	bb = (_HYButtonBar*)GetObject (2);
			bb->EnableButton (8,false);
			_UpdateLFMenu();
		}
		
		_HYTable* pl = (_HYTable*)GetObject (5);
		
		_SimpleList changedCells;
		
		for (long k=DF_ID_COLUMN; k<pl->cellTypes.lLength; k+=pl->horizontalSpaces.lLength)
		{
			pl->cellTypes.lData[k] &= (0xffffffff-HY_TABLE_BOLD);
			changedCells << k;
		}
		
		pl->_MarkCellsForUpdate (changedCells);
		
			
		savedLFNames.Clear();
		savedLFStates.Clear();
	}
	return true;
}

//__________________________________________________________
bool  _HYDataPanel::PurgeLFFilter (long fID)
{
	/*_HYTable* pl = (_HYTable*)GetObject (5);
	
	if (pl->cellTypes.lData[fID*pl->horizontalSpaces.lLength]&HY_TABLE_BOLD)
		return PurgeLF ();*/
		
	if (lfID>=0)
	{
		if (((_LikelihoodFunction*)likeFuncList (lfID))->GetTheFilters().Find (dataPartitions.lData[fID]) >=0)
			return PurgeLF();
		else
			if (cantDeleteLF && (treeVarReferences.lData[lfID] == -3))
				return PurgeLF();
	}
	
	return true;
}

//__________________________________________________________
void  _HYDataPanel::LongToPartData (long data, char& offset, bool& rev, long& code)
{
	offset = data&HY_DATAPANEL_OFFSETMASK;
	rev = data&HY_DATAPANEL_REVMASK;
	code = (data&HY_DATAPANEL_CODEMASK)>>16;
}

//__________________________________________________________
long _HYDataPanel::PartDataToLong (char offset, bool rev, long code)
{
	return (offset+((long)(rev)<<2)+(code<<16));
}

//__________________________________________________________
void  _HYDataPanel::LongToModelData (long data, int& model, int& options, int& freqs, int& classes)
{
	model = data&HY_DATAPANEL_MODELID;
	options = (data&HY_DATAPANEL_OPTIONS)>>16;
	freqs = (data&HY_DATAPANEL_FREQS)>>20;
	classes = (data&HY_DATAPANEL_RATES)>>24;
}

//__________________________________________________________
long _HYDataPanel::ModelDataToLong (int model, int options, int freqs, int classes)
{
	return (model+(options<<16)+(freqs<<20)+(classes<<24));
}


//__________________________________________________________
bool  _HYDataPanel::DataTypeChange (long filterID, long newType)
{
	long		 oldType;
	
	_DataSetFilter* 
				 theDF = (_DataSetFilter*)dataSetFilterList (dataPartitions.lData[filterID]);
	
	oldType = theDF->GetUnitLength();
	newType ++;
	bool	needToChop = false;
	long i;

	_SimpleList exH, exV, chopped;
	
	if ((newType>1)&&(newType!=oldType))
	{
		if (!EditPartitionProperties(filterID,newType))
			return false;
			
		bool		  direction = true;
		char		  readFrame = 0;
		
		LongToPartData (partData.lData[filterID],readFrame,direction,i);
		
		if  (direction)
		// reverse the order of sites in the partition
		{
			exH.Duplicate (&theDF->theOriginalOrder);
			exH.Flip();
			exV.Duplicate (&theDF->theNodeMap);
		}
		if (readFrame)
		// handle the reading frame
		{
			if (!exH.lLength)
			{
				exH.Duplicate (&theDF->theOriginalOrder);
				exV.Duplicate (&theDF->theNodeMap);
			}
			
			for (; readFrame; readFrame--)
				exH.Delete(0);
				
		}
		
		if (newType == 2)
			UnmarkSites(theDF->theOriginalOrder,filterID);
		
		if (exH.lLength)
			theDF->SetFilter(theDF->GetData(),newType,exV,exH,false);
			
		if (newType == 2)
		{
			MarkSites(theDF->theOriginalOrder,filterID);
			BuildThermometer();
			BuildMarksPane  ();	
		}
	}
	if (oldType!=newType)
	// need to change unit length
	{
		if (!PurgeLFFilter (filterID))
			return false;
		
		if (theDF->theOriginalOrder.lLength%newType)
		{
			_String warnMessage ("The number of sites in the filter ");
			warnMessage = warnMessage & *(_String*)dataSetFilterNamesList (dataPartitions.lData[filterID])
						  & " (" & _String ((long)theDF->theOriginalOrder.lLength) &") is not divisible by "
						  &_String (newType)&". I must chop off the last "&_String ((long)theDF->theOriginalOrder.lLength%newType)
						  & " sites to fix this problem.";
			if (!ProceedPrompt (warnMessage,(Ptr)this))
				return false;
			else
				needToChop = true;
		}
		exH.Duplicate (&theDF->theOriginalOrder);
		if (needToChop)
		{
			for (long k=(exH.lLength/newType)*newType; k<exH.lLength; k++)
				chopped<<exH.lData[k];
		}
		exV.Duplicate (&theDF->theNodeMap);

		theDF->SetFilter(theDF->GetData(),newType,exV,exH,false);
		if (newType == 3)
		// codon
		{
			SetCodonExclusions (theDF,filterID,i);
			theDF->SetDimensions();	
			partData.lData[filterID] = PartDataToLong (0,0,i);
		}
		else
		{
			partData.lData[filterID] = PartDataToLong (0,0,0);
			if (needToChop)
			{
				UnmarkSites(chopped,filterID);
				BuildMarksPane();
				BuildThermometer();
			}
		}
			
		modelReferences.lData[filterID] = 0;
		UpdatePartitionOperations();
		if (addedLines&HY_DATAPANEL_TRANSLATION)
		{
			_String dummy;
			if (newType==3) // codon data
				UpdateTranslationString(dummy,translatedSequence,true);
			else
				if (oldType==3)
				{
					for (filterID = 0; filterID<dataPartitions.lLength; filterID++)
					{
				 		theDF = (_DataSetFilter*)dataSetFilterList (dataPartitions.lData[filterID]);	
				 		if (theDF->GetUnitLength()==3)
				 			break;
					}
					if (filterID<dataPartitions.lLength)
						UpdateTranslationString(dummy,translatedSequence,true);	
					else
					// no codon partitions left; kill translation display
						AdjustStatusLine (2,true);
				}
		}
	}
	return true;
}
	
//__________________________________________________________
void  _HYDataPanel::ActivateInfoLines (bool onOff)
{
	if (onOff)
	{	
		SetTableDimensions (6,2);
		SetCell (0,0,GetObject (1));
		SetCell (0,1,GetObject (1));
		SetCell (1,0,GetObject (3));
		SetCell (1,1,GetObject (3));
		SetCell (3,0,GetObject (4));
		SetCell (3,1,GetObject (4));
		SetCell (2,0,GetObject (0));
		SetCell (2,1,GetObject (0));
		SetCell (4,0,GetObject (2));
		SetCell (4,1,GetObject (6));
		SetCell (5,0,GetObject (2));
		SetCell (5,1,GetObject (5));
		
	}
	else
	{
	//AddObject (sp);      // 0
	//AddObject (sc);      // 1
	//AddObject (b1);      // 2
	//AddObject (marks);   // 3
	//AddObject (sp2);     // 4
	//AddObject (partList);// 5
	//AddObject (partHead);// 6

		SetTableDimensions (5,2);
		SetCell (0,0,GetObject (1));
		SetCell (0,1,GetObject (1));
		SetCell (1,0,GetObject (3));
		SetCell (1,1,GetObject (3));
		SetCell (2,0,GetObject (0));
		SetCell (2,1,GetObject (0));
		SetCell (3,0,GetObject (2));
		SetCell (3,1,GetObject (6));
		SetCell (4,0,GetObject (2));
		SetCell (4,1,GetObject (5));
	}
}

//__________________________________________________________
bool  _HYDataPanel::AdjustStatusLine (long onOff, bool force, long preselected)
{
	_HYSequencePane*  sp2 = (_HYSequencePane*)components(4),
				   *  sp = (_HYSequencePane*)components(0);
					
	/* add consensus line */
	char 	oldAdded = addedLines, oldCount = sp2->speciesIndex.lLength;
	bool	ret = false;
	_String newDataString;
	if (preselected == -1)
		tainted = true;
	
	switch (onOff)
	{
		case 0:
		
			if (addedLines&HY_DATAPANEL_CONSENSUS)
			{
				addedLines &= 0xfe;
			}
			else
			{
				UpdateConsensusSequence (newDataString);
				addedLines |= HY_DATAPANEL_CONSENSUS;
				ret = true;
			}
			break;
			
		case 1:
			
			if (addedLines&HY_DATAPANEL_RATECLASS)
			{
				addedLines-=HY_DATAPANEL_RATECLASS;
			}
			else
			{
				if (lfID>=0)
				{
					if (lfID != lockedLFID)
					{
						_LikelihoodFunction *lf = (_LikelihoodFunction*)likeFuncList (lfID);
						_Matrix*			rates = lf->ConstructCategoryMatrix (true);
						_String 			rateS (32,true);
						long				k = 0,
						                    kk;
						
						for (long j=0; j<siteAssignments.lLength; j++)
						{
							if ((kk=siteAssignments.lData[j])>=0)
							{
								kk = ((_DataSetFilter*)dataSetFilterList (kk))->GetUnitLength();
								long v = (*rates)(0,k);
								char c;
								if (v>9)
									if (v<36)
										c = (char)('A'+v-10);
									else
										c = '*';
								else
									c = '0'+v;
								if (kk==1)
									rateS << c;
								else
								{
									for (long p = 0; p<kk; p++)	
										rateS << c;
									j+=kk-1;
								}
								k++;
							}
							else
								rateS << ' ';
						}
						rateS<<' ';
						rateS<<' ';
						rateS<<' ';
						rateS<<' ';
						rateS.Finalize();
						newDataString.Duplicate (&rateS);
						addedLines+=HY_DATAPANEL_RATECLASS;
						DeleteObject (rates);
						ret = true;
					}
					else
					{
						ProblemReport (lfCantKillWarning, (Ptr)this);
						return false;
					}
				}
			}
			break;
		
		case 2:
		{
			// check if can translate
			
			if (force)
			{
				addedLines &= 0xfb;
				break;
			}
			long k;
			_DataSetFilter* df;
			for (k=0; k<dataPartitions.lLength;k++)
			{
				df = (_DataSetFilter*)dataSetFilterList(dataPartitions.lData[k]);
				if (df->GetUnitLength()==3) break;
			}
			if (!(dataPartitions.lLength&&(k<dataPartitions.lLength)))
			{
				_String warnMsg ("You must define at least one codon partition before obtaining aminoacid translations. To translate an entire sequence, set the reference sequence to it, and choose the genetic code to translate with.");
				ProblemReport (warnMsg,(Ptr)this);
				return false;
			}
			// can translate, choose which sequence
			long sel;
			if (preselected<0)
			{
				_List choices, thisPair;
				_SimpleList	dummyChoices, dummySel, validChoices;
				dummyChoices<<0;
				dummyChoices<<1;
				_String desc ("All");
				thisPair && & desc;
				desc = "Translate all sequences and display the result in a new window";
				thisPair && & desc;
				choices && & thisPair;
				validChoices<<0;
				if (addedLines&HY_DATAPANEL_TRANSLATION)
				{
					thisPair.Clear();
					desc = "Turn off aminoacid translation display";
					validChoices<<1;
					thisPair && &none;
					thisPair && & desc;
					choices && & thisPair;
				}
				for (long k=0; k< sp->speciesIndex.lLength; k++)
				{
					thisPair.Clear();
					thisPair << sp->rowHeaders (sp->speciesIndex.lData[k]);
					desc = _String("Translate ")& *(_String*)thisPair(0);
					thisPair && & desc;
					validChoices << choices.lLength;
					choices && & thisPair;
				}
				sel = HandleListSelection (choices, dummyChoices, validChoices, "Choose sequence for translation", dummySel,1);
				if (sel<0) 
					return addedLines&HY_DATAPANEL_TRANSLATION;
			}
			else
				sel = preselected;
			
			
			if (addedLines&HY_DATAPANEL_TRANSLATION)
			{
				if (sel==1)
				{
					addedLines &= 0xfb;
					break;
				}
			}
			if (sel)
			{
				if (addedLines&HY_DATAPANEL_TRANSLATION)
					sel--;
				translatedSequence = sp->speciesIndex.lData[sel-1];
				UpdateTranslationString (newDataString,translatedSequence,false);
				ret = true;
				addedLines |= HY_DATAPANEL_TRANSLATION;
				if (!(oldAdded&HY_DATAPANEL_CONSENSUS))
					onOff--;
				if (!(oldAdded&HY_DATAPANEL_RATECLASS))
					onOff--;
			}
			else
			{
				
				_List			choices;
				
				_String mI 		("Map to missing data");
				choices && & mI;
				mI = "Any positions with ambiguities will be translated to an fully ambiguous amino-acid";
				choices && & mI;
				mI = "Most likely resolution";
				choices && & mI;
				mI = "Any ambiguity will be resolved to the most frequent one given other observed states in its alignment column; if an ambiguity can not be resolved to any of the observed states (e.g. a Y in a column of As and Gs), a random resolution will be provides.";
				choices && & mI;
				mI = "Random resolution";
				choices && & mI;
				mI = "Any ambiguity will be resolved randomly based on probabilities derived from observed states in its alignment column.";
				choices && & mI;
									
				char resolve = HandleListSelection (choices, "Ambiguities",(Ptr)this);	
				
				if (resolve	< 0) 
					return 0;
				
				_DataSet*		newDS = new _DataSet();
				_List* 			cache = new _List;
				GetTranslationString (newDataString, sp->speciesIndex.lData[0],resolve,-1,cache);
				for (k=0; k<newDataString.sLength; k++)
					newDS->AddSite (newDataString.sData[k]);
				newDS->AddName(*(_String*)sp->rowHeaders(sp->speciesIndex.lData[0]));
				for (k=1; k<sp->speciesIndex.lLength; k++)
				{
					GetTranslationString (newDataString, sp->speciesIndex.lData[k],resolve,-1,cache);
					newDS->AddName(*(_String*)sp->rowHeaders(sp->speciesIndex.lData[k]));
					for (long j=0; j<newDataString.sLength; j++)
							newDS->Write2Site (j,newDataString.sData[j]);
				}
				newDS->Finalize();
				DeleteObject (cache);
				
				_TranslationTable aaTable (defaultTranslationTable);
				aaTable.baseLength = 20;
				newDS->SetTranslationTable (&aaTable);
				newDS->SetNoSpecies (sp->speciesIndex.lLength);
				newDataString = (*(_String*)dataSetNamesList(dataSetID))&("_AA");
				AddDataSetToList (newDataString, newDS);
				_HYDataPanel* myTWindow = new _HYDataPanel(newDataString,newDataString);
				myTWindow->BringToFront();
				return addedLines&HY_DATAPANEL_TRANSLATION;
			}
			break;
		}

		case 3:
		{
			long sel,
				 sel2 = -1,
				 k;
				 
			if (preselected<0)
			{
				/*_List choices, thisPair;
				_SimpleList	dummyChoices, dummySel, validChoices;
				dummyChoices<<0;
				dummyChoices<<1;
				_String desc ("Turn off reference sequence display");
				if (addedLines&HY_DATAPANEL_REFERENCE)
				{
					validChoices<<0;
					thisPair && &none;
					thisPair && & desc;
					choices && & thisPair;
				}
				for (long k=0; k< sp->speciesIndex.lLength; k++)
				{
					thisPair.Clear();
					thisPair << sp->rowHeaders (sp->speciesIndex.lData[k]);
					desc = _String("Set reference sequence to ")& *(_String*)thisPair(0);
					thisPair && & desc;
					validChoices << choices.lLength;
					choices && & thisPair;
				}*/
				
				_List	theList,
						seqNames,
						gCodes;
					
				_String*iv;
						
				if (addedLines&HY_DATAPANEL_REFERENCE)
					iv = (_String*)sp->rowHeaders (referenceSequence);
				else
					iv = (_String*)sp->rowHeaders (0);
					

				if (addedLines&HY_DATAPANEL_REFERENCE)
					seqNames && & none;
				
				for (k=0; k< sp->speciesIndex.lLength; k++)
					seqNames << sp->rowHeaders (sp->speciesIndex.lData[k]);
				
				
				AddItemToPreferences (1|8,-1,"Sequence","Choose the reference sequence.","",nil,theList,false);
				AddItemToPreferences (0,PREFITEM_POPUP,"Sequence ID","Choose the reference sequence.",*iv,&seqNames,theList,false);
				
				if (dataType == HY_DATAPANEL_NUCDATA)
				{
					gCodes && & none;
					iv = &none;
					
					if (genCodeID>=0)
						iv = (_String*)(*((_List*)geneticCodes(genCodeID)))(0);
					
					for (k=0; k<geneticCodes.lLength; k++)
						gCodes << (_String*)(*((_List*)geneticCodes(k)))(0);
					
					AddItemToPreferences (1|8,-1,"Translation","Aminoacid translation options.","",nil,theList,false);
					AddItemToPreferences (0,PREFITEM_POPUP,"Genetic Code","Choose genetic code to translate with. Select \"None\" to display the nucleotide sequence.",
										   *iv,&gCodes,theList,false);		
				}
				
				if (HandlePreferences (theList, "Reference Sequence Setup", false))
				{
					sel  = seqNames.Find((*((_List*)theList.lData[4]))(1));
					if (dataType == HY_DATAPANEL_NUCDATA)
						sel2  = gCodes.Find((*((_List*)theList.lData[4]))(3))-1;
				}
				else
					return addedLines&HY_DATAPANEL_REFERENCE;				
			}
			else
				sel = preselected;
			if (addedLines&HY_DATAPANEL_REFERENCE)
			{
				if (sel==0) // turn off ref sequences
				{
					ret = false;
					addedLines &= 0xf7;
					genCodeID = -1;
				}
				else
				{
					sel = sp->speciesIndex.lData[sel-1];
					if ((referenceSequence!=sel)||(sel2!=genCodeID))
					{
						long m = ((_String*)statusData(0))->sLength-1;
						if (sel2<0)
						{
							for (long k=0; k<statusData.lLength; k++)
								((_String*)statusData(k))->sData[m] = ((_String*)sp->columnStrings(k))->sData[sel];
						}
						else
						{
							_String dummy;
							UpdateTranslationString (dummy, sel, false, sel2);
							for (long k=0; k<statusData.lLength; k++)
								((_String*)statusData(k))->sData[m] = dummy.sData[k];
						}
						referenceSequence = sel;
					}
					ret = true;
				}
			}
			else
			{
				addedLines|=HY_DATAPANEL_REFERENCE;
				referenceSequence = sp->speciesIndex.lData[sel];
				if (sel2<0)
				{
					_String newD ((unsigned long)sp->columnStrings.lLength,true);
					for (long k=0; k<sp->columnStrings.lLength; k++)
						newD << ((_String*)sp->columnStrings(k))->sData[referenceSequence];
					newD<<' ';
					newD<<' ';
					newD.Finalize();
					newDataString = newD;
				}
				else
					UpdateTranslationString (newDataString, referenceSequence, false, sel2);
				ret = true;
			}
						
			if (ret)
				genCodeID = sel2;
				
			if (!(oldAdded&HY_DATAPANEL_CONSENSUS))
				onOff--;
			if (!(oldAdded&HY_DATAPANEL_RATECLASS))
				onOff--;
			if (!(oldAdded&HY_DATAPANEL_TRANSLATION))
				onOff--;			
			break;
		}
	}
			
	if (ret)
	{
		if (oldAdded)
		{
			if (oldAdded!=addedLines)
				for (long k=0; k<newDataString.sLength; k++)
					((_String*)statusData(k))->Insert (newDataString.sData[k],onOff);
			else
				for (long k=0; k<newDataString.sLength; k++)
					((_String*)statusData(k))->sData[onOff] = newDataString.sData[k];				
		}				
		else
		{
			_String cs ('0');
			for (long k=0; k<newDataString.sLength; k++)
			{
				cs.sData[0] = newDataString.sData[k];
				statusData && & cs;
				sp2->InsertColumn ((_String*)statusData(k),-1);
			}
		}		
	}
	else
	{
		if (addedLines)
		{
			for (long k=0; k<statusData.lLength; k++)
				((_String*)statusData(k))->Delete(onOff,onOff);
		}
		else
			statusData.Clear();			
	}
	if (addedLines)
	{
		AdjustInfoNames();
		sp2->SetHeaders (&statusLines,false);
		sp2->BuildPane();
		sp2->_MarkForUpdate();
	}
	if (addedLines&&(!oldAdded))
	{
		sp2->blockWidth = sp->blockWidth;
		sp2->SetNameDisplayMode (sp->nameDisplayFlags,false);
		if (sp->startColumn!=0)
			sp2->startColumn = sp->startColumn;
		ActivateInfoLines (true);
	}
	else
		if (oldAdded&&(!addedLines))
		{
			sp2->columnStrings.Clear();
			sp2->rowHeaders.Clear();
			sp2->speciesIndex.Clear();
			ActivateInfoLines (false);
		}
		
	if (oldCount!=sp2->speciesIndex.lLength)
	{
		SetWindowRectangle (top,left,bottom,right);
	}
	return ret;
}
//__________________________________________________________

void  _HYDataPanel::AdjustInfoNames (void)
{
	_HYSequencePane*  sp = (_HYSequencePane*)components(0);
				   
	long			maxSpL = 0,t,k;
	
	for (long k=0; k<sp->rowHeaders.lLength;k++)
	{
		t = ((_String*)sp->rowHeaders(k))->sLength;
		if (t>maxSpL)
			maxSpL = t;
	}
	
	statusLines.Clear();
	if (addedLines&HY_DATAPANEL_CONSENSUS)
		statusLines&& &consensusInfoString;
	
	if (addedLines&HY_DATAPANEL_RATECLASS)
		statusLines&& &rateClassInfoString;
	
	if (addedLines&HY_DATAPANEL_TRANSLATION)
		statusLines&& sp->rowHeaders (translatedSequence);
	
	if (addedLines&HY_DATAPANEL_REFERENCE)
		statusLines&& sp->rowHeaders (referenceSequence);
	
	for (k=0;k<statusLines.lLength;k++)
	{
		_String *thisString = (_String*)statusLines(k);
		if (thisString->sLength>maxSpL)
			thisString->Trim (0,maxSpL-1);
		else
			if (k==0 && thisString->sLength<maxSpL)
			{
				_String paddedString ((unsigned long)maxSpL,true);
				paddedString<< thisString;
				for (long l = thisString->sLength; l<maxSpL; l++)
					paddedString<<' ';
				paddedString.Finalize();
				*thisString = paddedString;
			}
	}

}

//__________________________________________________________

void  _HYDataPanel::CodeTo3AA (_String& rec, long code, _SimpleList* genCode, _Parameter* freqVector)
{
	if (code<0)
	{
		if (freqVector)
		{
			long shift	   = 0,
				 codeB 	   = -1;
				 
			for (long k=0; k<64; k++)
				if (genCode->lData[k] == 10)
					shift++;
				else
					if (freqVector[k-shift] > 0.0)
						if (code < 0)
							code = genCode->lData[k];
						else
						{
							long code2 = genCode->lData[k];
							if (code2 != code)
							{
								if ((code == 13 && code2 == 15)
								    || (code == 15 && code2 == 13)
								    || (code == 12 && code2 == 16)
								    || (code == 16 && code2 == 12))
								{
									codeB = MIN(code,code2);
								   continue;
								}
							}
							rec = "???";
							return;
						}
			if (codeB >= 0)
				if (codeB == 13)
					code = 21;
				else
					code = 22;
						
		}
		else
		{
			rec = "???";
			return;
			
		}
	}
	else
		code = genCode->lData[code];
		
	switch (code)
	{
		case 0:
			rec = "Phe";
			break;
		case 1:
			rec = "Leu";
			break;
		case 2:
			rec = "Ile";
			break;
		case 3:
			rec = "Met";
			break;
		case 4:
			rec = "Val";
			break;
		case 5:
			rec = "Ser";
			break;
		case 6:
			rec = "Pro";
			break;
		case 7:
			rec = "Thr";
			break;
		case 8:
			rec = "Ala";
			break;
		case 9:
			rec = "Tyr";
			break;
		case 10:
			rec = "XXX";
			break;
		case 11:
			rec = "His";
			break;
		case 12:
			rec = "Gln";
			break;
		case 13:
			rec = "Asn";
			break;
		case 14:
			rec = "Lys";
			break;
		case 15:
			rec = "Asp";
			break;
		case 16:
			rec = "Glu";
			break;
		case 17:
			rec = "Cys";
			break;
		case 18:
			rec = "Trp";
			break;
		case 19:
			rec = "Arg";
			break;
		case 20:
			rec = "Gly";
			break;
		case 21:
			rec = "Asx";
			break;
		case 22:
			rec = "Glx";
			break;
	}
}

//__________________________________________________________

void  _HYDataPanel::IndexToCodon (long index, char* rec)
{
	char mapper[] = "ACGT";

	rec[0] = mapper[index/16];
	rec[1] = mapper[(index%16)/4];
	rec[2] = mapper[index%4];
}

//__________________________________________________________

_String*  _HYDataPanel::GetExclusionsFromCode (long index)
{
	_SimpleList* gencode = ((_SimpleList*)((*(_List*)geneticCodes(index))(2)));
	return 		 GetExclusionsFromList (gencode);
}

//__________________________________________________________

_String*  _HYDataPanel::GetExclusionsFromList (_SimpleList* gencode)
{
	_String* res = new _String (16,true);
	checkPointer (res);
	char     transl[4];
	transl[3] = 0;
	
	long index = 0;
	for (long k=0; k<64; k++)
		if (gencode->lData[k]==10)
		{
			IndexToCodon (k,transl);
			if (index)
				(*res)<<',';
			(*res)<< transl;
			index++;
		}
	res->Finalize();
	return  res;
}

//__________________________________________________________

_String*  _HYDataPanel::GetExclusionsFromExcList (_SimpleList* excl)
{
	_String* res = new _String (16,true);
	checkPointer (res);
	char     transl[4];
	transl[3] = 0;
	
	for (long k=0; k<excl->lLength; k++)
	{
		IndexToCodon (excl->lData[k],transl);
		if (k)
			(*res)<<',';
		(*res)<< transl;
	}
	res->Finalize();
	return  res;
}

//__________________________________________________________

void _HYDataPanel::SetCodonExclusions (_DataSetFilter* theDF, long filterID, long index)
{
	_String* exclusions = GetExclusionsFromCode (index);
	UnmarkSites (theDF->theOriginalOrder,filterID);
	theDF->SetExclusions(exclusions);
	MarkSites   (theDF->theOriginalOrder,filterID);
	BuildMarksPane();
	BuildThermometer();
	DeleteObject (exclusions);
}

//__________________________________________________________

_String*  _HYDataPanel::GetMatrixFromCode (long index)
{
	_String* res = new _String (16,true);
	checkPointer (res);
	
	_SimpleList* gencode = ((_SimpleList*)((*(_List*)geneticCodes(index))(2)));
	(*res)<<'{';
	(*res)<<'{';
	for (long k=0; k<64; k++)
	{
		_String c (gencode->lData[k]);
		if (k)
			(*res)<<',';
		(*res)<< &c;
	}
	(*res)<<'}';
	(*res)<<'}';
	res->Finalize();
	return  res;
}

//__________________________________________________________

void  _HYDataPanel::ShowConstantSites (bool deletions, bool relaxed, bool sequences)
{
	_HYSequencePane *sp = (_HYSequencePane*)GetObject (0);
	_HYTable		*pl = (_HYTable*)GetObject (5);
	
	_DataSetFilter  *df = (_DataSetFilter*)dataSetFilterList(dataPartitions.lData[pl->GetFirstRowSelection()]);
	_SimpleList		constantSites;
	
	if (deletions)
	{
		if (sequences)
		{
			_AVLList	  seqswithdels (&constantSites);
			
			for (long i = 0; i<df->theFrequencies.lLength; i++)
				df->HasDeletions(i,&seqswithdels);
				
			seqswithdels.ReorderList ();
			
			if (constantSites.lLength)
			{
				_SimpleList	 remap,
							 displayOrder;
							 
				_AVLListX	 reordered (&remap);
				
				for (long k = 0; k<sp->speciesIndex.lLength; k++)
					reordered.Insert ((BaseRef)sp->speciesIndex.lData[k],k);
				
				for (long k2 = 0; k2 < constantSites.lLength; k2++)
				{
					long f = reordered.Find (BaseRef(constantSites.lData[k2]));
					if (f>=0)
						displayOrder << reordered.GetXtra (f);
				}
				
				constantSites.Clear();
				constantSites.Duplicate (&displayOrder);
				constantSites.Sort();
			}
		}
		else
		{
			for (long i = 0; i<df->theFrequencies.lLength; i++)
			{
				if (df->HasDeletions(i))
					df->FindAllSitesLikeThisOne(df->theMap[i*df->GetUnitLength()],constantSites);
			}
			
			if (relaxed) // constant delete sites
			{
				_SimpleList sites2, sites1 (constantSites);
				for (long i = 0; i<df->theFrequencies.lLength; i++)
				{
					if (df->IsConstant(i,false))
						df->FindAllSitesLikeThisOne(df->theMap[i*df->GetUnitLength()],sites2);
				}
				sites1.Sort();
				sites2.Sort();
				constantSites.Intersect (sites1, sites2);
			}
		}
	}
	else
		for (long i = 0; i<df->theFrequencies.lLength; i++)
		{
			if (df->IsConstant(i,relaxed))
				df->FindAllSitesLikeThisOne(df->theMap[i*df->GetUnitLength()],constantSites);
		}
	
	_String  	outWord;
	_Parameter	fnd;
	long		fndc;
	if (sequences)
	{
		sp->selection.Clear();
		sp->vselection.Clear();
		sp->vselection.Duplicate(&constantSites);
		outWord = "sequences";
		fnd  = ((_Parameter)sp->vselection.lLength)/sp->speciesIndex.lLength*100.;
		fndc = sp->vselection.lLength;
	}	
	else
	{
		constantSites.Sort();
		sp->selection.Clear();
		sp->selection.Duplicate(&constantSites);
		outWord = "sites";
		fnd = 100.*((_Parameter)sp->selection.lLength)/df->GetFullLengthSpecies();
		fndc = sp->selection.lLength/df->GetUnitLength();
	}
	UpdateSelDepPartitionOperations ();
	_UpdateSelectionChoices (sp->selection.lLength);
	sp->BuildPane();
	sp->_MarkForUpdate();
	char 	buffer[128];
	sprintf (buffer,"%d(%4.4g%%) %s found\n",fndc,fnd,outWord.sData);
	BufferToConsole (buffer);
}

//__________________________________________________________

void  _HYDataPanel::ShowDuplicateSequences (bool relaxed)
{
	_HYSequencePane *sp = (_HYSequencePane*)GetObject (0);
	_HYTable		*pl = (_HYTable*)GetObject (5);
	
	_DataSetFilter  *df = (_DataSetFilter*)dataSetFilterList(dataPartitions.lData[pl->GetFirstRowSelection()]);
	
	_SimpleList 	duplicateSequences,
					otherInstance,
					stillToCheck;
					
	long			vd;
	
	for (vd = 0; vd < df->NumberSpecies(); vd++)
		stillToCheck << vd;
					
	long			sitesToCheck = df->GetFullLengthSpecies ()/df->GetUnitLength();
					
	vd	= df->GetDimension(true);
	
	_Parameter		*translatedVector = new _Parameter [vd],
					*translatedVector2= new _Parameter [vd];
					
	checkPointer    (translatedVector);
	checkPointer    (translatedVector2);
	
	_String			state1 (df->GetUnitLength(),false),
					state2 (df->GetUnitLength(),false);
	
	for (long k=0; k<stillToCheck.countitems()-1; k++)
	{
		for (long l=k+1; l<stillToCheck.countitems(); l++)
		{
			bool checkState = true;
			for (long m=0; m<sitesToCheck; m++)
			{
				df->RetrieveState (m,stillToCheck.lData[l], state1);
				df->RetrieveState (m,stillToCheck.lData[k], state2);
				long idx1 = df->Translate2Frequencies (state1, translatedVector,  true),
					 idx2 = df->Translate2Frequencies (state2, translatedVector2, true);
					 
				if (idx1>=0)
				{
					if (idx1==idx2)
						continue;
					else
					{
						checkState = false;
						break;
					}
				}
				else
				{
					if (relaxed)
					{
						bool first  = true,
							 second = true;
							 
						for (long t = 0; (first||second)&&(t<vd); t++)
						{
							if (translatedVector[t]>0.0)
								second &= (translatedVector2[t]>0.0);
							if (translatedVector2[t]>0.0)
								first  &= (translatedVector[t]>0.0);
						}
						
						if (!(first||second))
						{
							checkState = false;
							break;
						}
					}
					else
					{
						for (long t = 0; t<vd; t++)
							if (translatedVector[t]!=translatedVector2[t])
							{
								checkState = false;
								break;
							}
					}
				}
			}
			
			if (checkState)
			{
				duplicateSequences << stillToCheck.lData[l];
				otherInstance << stillToCheck.lData[k];
				stillToCheck.Delete (l);
				l--;
			}
		}	
	}
	
	delete (translatedVector);
	delete (translatedVector2);
	
	
	SortLists (&duplicateSequences, &otherInstance);
	
	sp->selection.Clear();
	sp->vselection.Clear();
	for (long k=0; k<duplicateSequences.lLength; k++)
		sp->vselection << sp->speciesIndex.Find(df->theNodeMap.lData[duplicateSequences.lData[k]]);
		
	sp->vselection.Sort();
	
	UpdateSelDepPartitionOperations ();
	_UpdateSelectionChoices (0);
	
	sp->BuildPane();
	sp->_MarkForUpdate();
	
	char	buffer[128];
	sprintf (buffer,"%d(%4.4g%%) duplicate sequences found\n", duplicateSequences.lLength ,((_Parameter)duplicateSequences.lLength)/sp->speciesIndex.lLength*100.);
	BufferToConsole (buffer);
	for (long idx = 0; idx < duplicateSequences.lLength; idx++)
	{
		BufferToConsole ("\t");
		StringToConsole (*df->GetSequenceName (duplicateSequences.lData[idx]));
		BufferToConsole (" matched ");
		StringToConsole (*df->GetSequenceName (otherInstance.lData[idx]));
		NLToConsole		();
	}
}


//__________________________________________________________

char  _HYDataPanel::CodeToAA (long code, _SimpleList* genCode, _Parameter* freqVector)
{
	if (code<0)
	{
		if (freqVector)
		{
			long shift	   = 0,
				 codeB 	   = -1;
				 
			for (long k=0; k<64; k++)
				if (genCode->lData[k] == 10)
					shift++;
				else
					if (freqVector[k-shift] > 0.0)
						if (code < 0)
							code = genCode->lData[k];
						else
						{
							long code2 = genCode->lData[k];
							if (code2 != code)
							{
								if ((code == 13 && code2 == 15)
								    || (code == 15 && code2 == 13)
								    || (code == 12 && code2 == 16)
								    || (code == 16 && code2 == 12))
								{
									codeB = MIN(code,code2);
								   continue;
								}
							}
							return '?';
						}
			if (codeB >= 0)
				if (codeB == 13)
					code = 21;
				else
					code = 22;
						
		}
		else
			return '?';
	}
	else
		code = genCode->lData[code];
		
	switch (code)
	{
		case 0:
			return 'F';
			break;
		case 1:
			return 'L';
			break;
		case 2:
			return 'I';
			break;
		case 3:
			return 'M';
			break;
		case 4:
			return 'V';
			break;
		case 5:
			return 'S';
			break;
		case 6:
			return 'P';
			break;
		case 7:
			return 'T';
			break;
		case 8:
			return 'A';
			break;
		case 9:
			return 'Y';
			break;
		case 10:
			return 'X';
			break;
		case 11:
			return 'H';
			break;
		case 12:
			return 'Q';
			break;
		case 13:
			return 'N';
			break;
		case 14:
			return 'K';
			break;
		case 15:
			return 'D';
			break;
		case 16:
			return 'E';
			break;
		case 17:
			return 'C';
			break;
		case 18:
			return 'W';
			break;
		case 19:
			return 'R';
			break;
		case 20:
			return 'G';
			break;
		case 21:
			return 'B';
			break;
		case 22:
			return 'Z';
			break;
	
	}
	
	return 'X';
}

//__________________________________________________________

void		_HYDataPanel::DFExclusionsToString	  (_DataSetFilter* theDF, _String& res)
{
	_String			exclusions (64,true);
	if ((theDF->GetUnitLength()==3)&&theDF->theExclusions.lLength)
	{
		exclusions <<  theDF->ConvertCodeToLetters (theDF->theExclusions.lData[0],3);
		for (long kk=1; kk<theDF->theExclusions.lLength; kk++)
		{
			exclusions << ',';
			exclusions <<  theDF->ConvertCodeToLetters (theDF->theExclusions.lData[kk],3);
		}
	}
	exclusions.Finalize();
	res = exclusions;
}


//__________________________________________________________

void  _HYDataPanel::DisplayParameterTable (void)
{
	if (lfID>=0)
	{
		_String		windowName;
		windowName = _String ("Likelihood parameters for ") & *dataSetName;
		long	k = FindWindowByName (windowName);
		if (k>=0)
		{
			#ifdef __MAC__
				_HYPlatformWindow* thisWindow = (_HYPlatformWindow*)windowObjects(k);
				thisWindow->_Activate();
			#else
				_HYWindow* thisWindow = (_HYWindow*)windowObjectRefs(k);
				thisWindow->BringToFront();
			#endif						
		}
		else
		{
			_HYParameterTable* newPT = new _HYParameterTable (windowName,lfID);
			newPT->_Zoom (true);
			newPT->BringToFront();
		}
	}
}

//__________________________________________________________

void  _HYDataPanel::OpenGeneralBSWindow (void)
{
	if (lfID>=0)
	{
		#ifndef __HYPHY_GTK__
			_String		windowName;
			windowName = _String ("Bootstrap setup ") & *dataSetName;
			long	k = FindWindowByName (windowName);
			if (k>=0)
			{
				_HYGeneralBootstrapWindow* thisWindow = (_HYGeneralBootstrapWindow*)windowObjectRefs(k);
				thisWindow->_Activate();
				thisWindow->SetNullLF (lfID);					
			}
			else
			{
				_HYGeneralBootstrapWindow* newPT = new _HYGeneralBootstrapWindow (windowName,lfID);
				newPT->BringToFront();
			}
		#endif
	}
}

//__________________________________________________________

_String*  _HYDataPanel::LFSnapshot (void)
{
	if (lfID>=0)
	{
		stashParameter (likefuncOutput,4,true);
		_String * result = (_String*)((_LikelihoodFunction*)likeFuncList(lfID))->toStr();
		stashParameter (likefuncOutput,4,false);
		tainted  = true;
		return   result;
	}
	return nil;
}


//__________________________________________________________

bool 	_HYDataPanel::LFRestore (long index)
{
	if (lfID>=0)
		if (lfID!=lockedLFID)
		{
			if (index<savedLFStates.lLength)
			{
				_String dupList (*(_String*)savedLFStates(index));
				_ExecutionList exl (dupList);
				exl.Execute();
				((_LikelihoodFunction*)likeFuncList(lfID))->RescanAllVariables();
				RefreshCategoryVars();
				return true;
			}
		}
		else
			ProblemReport (lfCantKillWarning);
	return false;
}

//__________________________________________________________

long  	_HYDataPanel::GetHypothesis (bool alt)
{
	_String suffix;
	if (alt)
		suffix = alterSuffix;
	else
		suffix = nullSuffix;
		
	for (long k=0; k<savedLFNames.lLength; k++)
		if (((_String*)savedLFNames(k))->endswith (suffix))
			return k;
			
	return -1;
}

//__________________________________________________________

void  	_HYDataPanel::SetHypothesis (_String* s, bool alt)
{
	_String suffix,
			nSuffix;
	if (alt)
	{
		suffix   = alterSuffix;
		nSuffix  = nullSuffix;
	}
	else
	{
		suffix   = nullSuffix;
		nSuffix  = alterSuffix;
	}
		
	long k = GetHypothesis (alt),
		 f = savedLFNames.Find (s);
		 
	if (k>=0)
		((_String*)savedLFNames(k))->Trim(0,((_String*)savedLFNames(k))->sLength-suffix.sLength-1);
		
	if (f>=0)
	{
		if (s->endswith (nSuffix))
			s->Trim (0, s->sLength-nSuffix.sLength-1);
		*((_String*)savedLFNames(f)) = *s & suffix;
	}
	
	tainted = true;
}	

//__________________________________________________________

long  	_HYDataPanel::FindLFState (_String s)
{
	long f = savedLFNames.Find (&s);
	if (f<0)
	{
		_String s2;
		s2 = s &  alterSuffix;
		f = savedLFNames.Find (&s2);
		if (f<0)
		{
			s2 = s & nullSuffix;
			f = savedLFNames.Find (&s2);
		
		}
	}
	return f;
}	

//__________________________________________________________

_String  	_HYDataPanel::GetLFStateName (long k)
{
	_String res;
	
	if ((k>=0)&&(k<savedLFNames.lLength))
	{
		res = *(_String*)savedLFNames (k);
		if (res.endswith (nullSuffix))
			return res.Cut (0, res.sLength-nullSuffix.sLength-1);
		
		if (res.endswith (alterSuffix))
			return res.Cut (0, res.sLength-alterSuffix.sLength-1);
	}
	
	return res;
}	

//__________________________________________________________

_String*  	_HYDataPanel::GetLFStateString (long k)
{
	if ((k>=0)&&(k<savedLFNames.lLength))
	{
		return (_String*)savedLFStates (k);
	}
	
	return nil;
}

//__________________________________________________________

void 	_HYDataPanel::FindFunction (void)
{
	_String tPrompt ("Find data:"),
			cPrompt ("Use regular expressions");
			
	_List	searchInOptions;
	
	_String opt ("Search sites");
	searchInOptions && & opt;
	opt = "Search sequences, show sequences";
	searchInOptions && & opt;
	opt = "Search sequences, show sites";
	searchInOptions && & opt;
	opt = "Search sequence names";
	searchInOptions && & opt;
	
	_HYSequencePane*sp = (_HYSequencePane*)GetObject(0);
	
	_DataSetFilter *selectedFilter = nil;
	
	_HYTable	   *pl = (_HYTable*)GetObject (5);
	_SimpleList	   filterSel;
	
	pl->GetRowSelection (filterSel);

	if (sp->vselection.lLength==1)
	{
		opt = "Search in selected sequence";
		searchInOptions && & opt;	
	}
	
	if (filterSel.lLength == 1)
	{
		selectedFilter = (_DataSetFilter*)dataSetFilterList (dataPartitions.lData[filterSel.lData[0]]);
		if (selectedFilter->GetUnitLength () == 3)
		{
			searchInOptions && & menuSeparator;	
			opt = "Search for amino-acid motifs";
			searchInOptions && & opt;	
		}
			
	}
	
	long	seqSearch;
			
	if (EnterStringDialogWithPulldown (dataPanelSearchTerm,tPrompt,cPrompt,seqSearch, searchInOptions, nil, regExpSearch,findPanelSelection,(Ptr)this)&&dataPanelSearchTerm.sLength)
	{
		_SimpleList 	matches;
		_DataSet*		myData = (_DataSet*)dataSetList (dataSetID);
		_SimpleList		eligibleSequences;
		
		findPanelSelection = seqSearch;
		
		if (seqSearch >= 4)
			regExpSearch = true;
		
		Ptr				regExpie = nil;
		
		if (regExpSearch)
		{
			int 	   errCode;
			regExpie = PrepRegExp(&dataPanelSearchTerm, errCode, false);
			if (errCode)
			{
				_String errMsg = GetRegExpError(errCode);
				ProblemReport (errMsg, (Ptr)this);
				return;
			}
		}
		
		if (seqSearch==3)
		{
			_List * seqNames = &myData->GetNames();
			
			for (long idx = 0; idx<sp->speciesIndex.lLength;idx++)
			{
				if (regExpie)
				{
					_SimpleList found;
					((_String*)(*seqNames)(sp->speciesIndex.lData[idx]))->RegExpMatch (regExpie, found);
					if (found.lLength)
						matches << idx;
				}
				else
					if (((_String*)(*seqNames)(sp->speciesIndex.lData[idx]))->FindAnyCase (dataPanelSearchTerm)>=0)
						matches << idx;
			}				
		}
		else
		{

			if (seqSearch >= 5) // a.a. motif
			{
				if (sp->vselection.lLength)
					eligibleSequences.Duplicate (&sp->vselection);
				else
					for (long k=0; k<sp->speciesIndex.lLength; k++)
						eligibleSequences << k;

				_String    stran (selectedFilter->theOriginalOrder.lLength/3,false);
				
				_AVLList   avl (&matches);
								
				for (long idx = 0; idx < eligibleSequences.lLength; idx++)
				{
					long		 			seqIdx =  sp->speciesIndex.lData[idx];
					bool					doRev  = GetTranslationString (stran,sp->speciesIndex.lData[eligibleSequences.lData[idx]],0,filterSel.lData[0]);
					
					_SimpleList 			found;
					stran.RegExpMatchAll (regExpie, found);	
					if (doRev)
						for (long idx3 = 0; idx3 < found.lLength; idx3+=2)
							for (long idx4 = found.lData[idx3]; idx4<=found.lData[idx3+1]; idx4++)
							{
								seqIdx = selectedFilter->theOriginalOrder.lLength-idx4*3;
								
								avl.Insert ((BaseRef)selectedFilter->theOriginalOrder.lData[seqIdx-1]);
								avl.Insert ((BaseRef)selectedFilter->theOriginalOrder.lData[seqIdx-2]);
								avl.Insert ((BaseRef)selectedFilter->theOriginalOrder.lData[seqIdx-3]);
							}
					else
						for (long idx3 = 0; idx3 < found.lLength; idx3+=2)
							for (long idx4 = found.lData[idx3]; idx4<=found.lData[idx3+1]; idx4++)
							{
								avl.Insert ((BaseRef)selectedFilter->theOriginalOrder.lData[idx4*3]);
								avl.Insert ((BaseRef)selectedFilter->theOriginalOrder.lData[idx4*3+1]);
								avl.Insert ((BaseRef)selectedFilter->theOriginalOrder.lData[idx4*3+2]);
							}
					
				}
				
				avl.ReorderList();
				seqSearch = 0;
			}
			else
			{
							
				_SimpleList	*theMap = &myData->GetTheMap(),
						matchedSites ((unsigned long)myData->NoOfUniqueColumns(),0,0);
							   						   
				//pl->GetRowSelection (filterSel);
					
				if (filterSel.lLength)
				{
					for (long d = 0; d<filterSel.lLength; d++)
					{
						_DataSetFilter * aFilter = (_DataSetFilter*)dataSetFilterList (dataPartitions.lData[filterSel.lData[d]]);
						for (long dd = 0; dd<aFilter->theMap.lLength; dd++)
							matchedSites.lData [theMap->lData[aFilter->theMap.lData[dd]]] = 1;
					}
				
					for (long k=0; k<matchedSites.lLength;matchedSites.lData[k++]=0)
						if (matchedSites.lData[k])
							eligibleSequences << k;
				}
				else
					for (long k=0; k<myData->NoOfUniqueColumns();k++)
						eligibleSequences << k;		
				
				if (seqSearch)
				{
					_AVLList*	  avm = nil;		
					if (seqSearch == 2)
					{
						avm = new _AVLList (&matches);
						checkPointer (avm);
					}
					for (long idx = (seqSearch==4)?sp->vselection.lData[0]:0; idx<sp->speciesIndex.lLength;idx++)
					{
						_String		 seqStringC (eligibleSequences.lLength, false),
									 seqString  (theMap->lLength, false);
						
						//long		 seqIdx =  dataWrapper?dataWrapper->theNodeMap[sp->speciesIndex.lData[idx]]:sp->speciesIndex.lData[idx];
						long		 seqIdx =  sp->speciesIndex.lData[idx];
						
						for (long idx2 = 0; idx2 < eligibleSequences.lLength; idx2++)
							seqStringC.sData[idx2] = ((_Site*)(myData->_List::operator ())(eligibleSequences.lData[idx2]))->sData[seqIdx];
							
						for (long idx2 = 0; idx2 < theMap->lLength; idx2++)
							seqString.sData[idx2] = seqStringC.sData[theMap->lData[idx2]];
						
						if (seqSearch == 2)
						{
							if (regExpie)
							{
								_SimpleList found;
								seqString.RegExpMatchAll (regExpie, found);	
								for (long idx3 = 0; idx3 < found.lLength; idx3+=2)
									for (long idx4 = found.lData[idx3]; idx4<=found.lData[idx3+1]; idx4++)
										avm->Insert ((BaseRef)idx4);

							}
							else
							{
								long f = seqString.FindAnyCase (dataPanelSearchTerm);
								while (f>=0)
								{
									for (long k2 = 0; k2 < dataPanelSearchTerm.sLength; k2++,f++)
										avm->Insert ((BaseRef)f);
										
									f = seqString.FindAnyCase (dataPanelSearchTerm,f,-1);
								}				
							}
						}
						
						else
						{
							if (regExpie)
							{
								_SimpleList found;
								if (seqSearch < 4)
								{
									seqString.RegExpMatch (regExpie, found);	
									if (found.lLength)
										matches << idx;
								}
								else
								{
									seqString.RegExpMatchAll (regExpie, found);	
									for (long idx3 = 0; idx3 < found.lLength; idx3+=2)
										for (long idx4 = found.lData[idx3]; idx4<=found.lData[idx3+1]; idx4++)
											matches << idx4;

									break;
								}
							}
							else
								if (seqString.FindAnyCase (dataPanelSearchTerm)>=0)
									matches << idx;
						}
					}		
					if (seqSearch == 2)
					{
						avm->ReorderList();
						DeleteObject (avm);
					}
					if (seqSearch!=1)
						seqSearch = 0;		
				}
				else
				{
					bool		   permuteString = false;
					
					for (long k=0; k<sp->speciesIndex.lLength;k++)
						if (sp->speciesIndex.lData[k]!=k)
						{
							permuteString = true;
							break;
						}
						
					if (permuteString)
					{
						for (long idx = 0; idx<eligibleSequences.lLength;idx++)
						{
							_String * col = ((_Site*)(myData->_List::operator ())(eligibleSequences.lData[idx])),
									pCol (10L, true);
									
							for (long idx2 = 0; idx2 < sp->speciesIndex.lLength; idx2++)
								pCol << col->sData[sp->speciesIndex.lData[idx2]];
								
							pCol.Finalize();
							if (regExpie)
							{
								_SimpleList	 found;
								pCol.RegExpMatch (regExpie, found);
								matchedSites.lData[eligibleSequences.lData[idx]] = (found.lLength>0);
							}
							else
								matchedSites.lData[eligibleSequences.lData[idx]] = (pCol.FindAnyCase (dataPanelSearchTerm)>=0);
								
						}
					}
					else
					{
						for (long idx = 0; idx<eligibleSequences.lLength;idx++)
						{
							if (regExpie)
							{
								_SimpleList	 found;
								((_Site*)(myData->_List::operator ())(eligibleSequences.lData[idx]))->RegExpMatch (regExpie, found);
								matchedSites.lData[eligibleSequences.lData[idx]] = (found.lLength>0);
							}
							else
								matchedSites.lData[eligibleSequences.lData[idx]] =
									 (((_Site*)(myData->_List::operator ())(eligibleSequences.lData[idx]))->FindAnyCase (dataPanelSearchTerm)>=0);
						}
					}
				}
				if (seqSearch == 0)
					for (long idx = 0; idx < theMap->lLength; idx++)
						if (matchedSites.lData[theMap->lData[idx]])
							matches << idx;
			}		
		}
		if (regExpie)
			 FlushRegExp (regExpie);
			 
		sp->SelectRange (matches,seqSearch);
		sp->ScrollToHSelection();
		UpdateSelDepPartitionOperations();
		_UpdateSelectionChoices (sp->selection.lLength);

	}
}



//__________________________________________________________

long	 _HYDataPanel::FindGeneticCodeByExcl 	(_String* excl)
{
	for (long k=0; k<geneticCodes.lLength; k++)
	{
		_String* 	  myExcl = GetExclusionsFromCode (k);
		if (excl->Equal (myExcl))
		{
			DeleteObject (myExcl);
			return k;
		}
		DeleteObject (myExcl);
	}
	return -1;
}

//__________________________________________________________

void	 _HYDataPanel::HandleFontChange 	(void)
{
	_HYSequencePane* sp = (_HYSequencePane*)components (0);		
			
	_HYFontDialog * fD = new _HYFontDialog (sp->GetFont(),this);
	fD->Activate();
}


//__________________________________________________________

void	 _HYDataPanel::SetFont 	(_HYFont& newFont)
{
	_HYSequencePane* sp = (_HYSequencePane*)components (0),
				   * sp2= (_HYSequencePane*)components (4);
		
	sp->SetFont (newFont);
	sp2->SetFont (newFont);
	SetWindowRectangle (top,left,bottom,right);
}

//__________________________________________________________

void ReadGeneticCodes (void)
{
	_List receptacle, 
		  addedNames;
		  
	_SimpleList addedIndex,
				codeIndex;
	
	_String line1,line2;
	line1 = "Universal";
	_List  uncode;
	_SimpleList untable;
	uncode && & line1;
	line2 = "Built-in Universal Genetic Code";
	uncode && & line2;
	untable<<14;
	untable<<13; 
	untable<<14; 
	untable<<13;
	untable<<7;  
	untable<<7;  
	untable<<7;  
	untable<<7; 
	untable<<19; 
	untable<<5; 
	untable<<19;  
	untable<<5; 
	untable<<2;  
	untable<<2; 	
	untable<<3;  
	untable<<2; 
	untable<<12; 
	untable<<11; 
	untable<<12; 
	untable<<11;
	untable<<6;  
	untable<<6;  
	untable<<6;  
	untable<<6; 
	untable<<19; 
	untable<<19; 
	untable<<19; 
	untable<<19;
	untable<<1;  
	untable<<1;  
	untable<<1;  
	untable<<1; 
	untable<<16; 
	untable<<15; 
	untable<<16; 
	untable<<15;
	untable<<8;  
	untable<<8;  
	untable<<8;  
	untable<<8; 
	untable<<20; 
	untable<<20; 
	untable<<20; 
	untable<<20;
	untable<<4;  
	untable<<4;  
	untable<<4;  
	untable<<4; 
	untable<<10; 
	untable<<9; 
	untable<<10; 
	untable<<9; 
	untable<<5;  
	untable<<5;  
	untable<<5;  
	untable<<5; 
	untable<<10; 
	untable<<17; 
	untable<<18; 
	untable<<17;
	untable<<1;  
	untable<<0;  
	untable<<1;  
	untable<<0; 		
	uncode && & untable;
	geneticCodes && & uncode;
	addedNames && & line1;
	addedIndex << 0;
	codeIndex << 0;
		
	_String pathToGeneticCodes;
	pathToGeneticCodes = baseDirectory&"GeneticCodes";
	
	ScanDirectoryForFileNames (pathToGeneticCodes,receptacle,true);
	
	for (long k=0; k<receptacle.lLength;k++)
	{
		FILE * thisFile = doFileOpen (((_String*)receptacle(k))->sData,"rb");
		if (thisFile)
		{	
			_String buffer (thisFile),
					line3,
					line4;
					
			long	curCodeIndex;
					
			fclose (thisFile);
			line4 = _ElementaryCommand::FindNextCommand(buffer);
			line1 = _ElementaryCommand::FindNextCommand(buffer);
			line2 = _ElementaryCommand::FindNextCommand(buffer);
			line3 = _ElementaryCommand::FindNextCommand(buffer);
			line1.Trim (0,line1.sLength-2);
			line2.Trim (0,line2.sLength-2);
			line3.Trim (0,line3.sLength-2);
			line4.Trim (0,line4.sLength-2);
			
			curCodeIndex = line4.toNum();
			
			if (line1.sLength&&line2.sLength&&line3.sLength&&(curCodeIndex>=0))
			{
				if (line1.IsValidIdentifier()&&(addedNames.Find(&line1)==-1)&&(codeIndex.Find(curCodeIndex)==-1))
				{
					_List partition;
					_SimpleList translations;
					_ElementaryCommand::ExtractConditions (line3,0,partition,',');
					for (long l=0; l<partition.lLength; l++)
					{
						_String* thisLine = (_String*)partition(l);
						long aa = thisLine->FirstNonSpaceIndex();
						if (aa>=0)
						{
							thisLine->Trim(aa,-1);
							if (thisLine->sLength)
							{
								aa = thisLine->toNum();
								if ((aa>=0)&&(aa<=20))
								{
									translations << aa;
								}
							}
						}
					}
					if (translations.lLength == 64)
					{
						_List goodCode;
						goodCode && & line1;
						goodCode && & line2;
						goodCode && & translations;
						addedIndex << geneticCodes.lLength;
						codeIndex << curCodeIndex;
						geneticCodes && & goodCode;
						addedNames && & line1;
					}
				}
			}
		}
	}
	if (geneticCodes.lLength>1)
	{
		SortLists (&codeIndex, &addedIndex);
		
		long	  k;
		
		for (k=1; k<codeIndex.lLength; k++)
		{
			if (codeIndex.lData[k] - codeIndex.lData[k-1] != 1)
				break;
		}
		
		if ((k<codeIndex.lLength)||(codeIndex.lData[1]!=1))
		{
			_String warnMsg ("There is a genetic code (or two) missing from the '");
			warnMsg = warnMsg & pathToGeneticCodes & "'. You may experience problems with saved GUI analyses which include codon data with non Universal codes.";
			StringToConsole (warnMsg);	
		}
		
		_List   sortedCodes;
		for (k=0; k<geneticCodes.lLength; k++)
			sortedCodes << geneticCodes (codeIndex.lData[k]);
			
		geneticCodes.Clear();
		geneticCodes.Duplicate (&sortedCodes);
		
		_String status = _String("Loaded ") & (long)geneticCodes.lLength-1 & " genetic code tables from " & pathToGeneticCodes.getStr() & '\n';		
		StringToConsole (status);
	}
	else
	{
		_String warnMsg ("I couldn't find any valid genetic code tables in '");
		warnMsg = warnMsg & pathToGeneticCodes & "'. You can still use the Universal code which is built-in.";
		StringToConsole (warnMsg);	
	}
}

//__________________________________________________________

void ReadModelTemplates (void)
{
	_List receptacle, addedNames;
	
	_String line1;
	_String pathToModelTemplates;
	_SimpleList	modelParams;

	pathToModelTemplates = baseDirectory&"SubstitutionModels";
	
	ScanDirectoryForFileNames (pathToModelTemplates,receptacle,true);
	
	for (long k=0; k<receptacle.lLength;k++)
	{
		FILE * thisFile = doFileOpen (((_String*)receptacle(k))->sData,"rb");
		if (thisFile)
		{	
			_String buffer (thisFile);
			fclose (thisFile);
			long   g = batchLanguageFunctionNames.lLength;
			_ExecutionList   thisList;
			terminateExecution = false;
			thisList.BuildList (buffer);
			thisList.Execute();

			if (terminateExecution==false)
			{
			/* check for variables and functions */
				long 	 popFunc = batchLanguageFunctionNames.Find (&modelFunction),
						 efvFunc = batchLanguageFunctionNames.Find (&efvFunction),
						 codonFunc = batchLanguageFunctionNames.Find (&buildCodonFrequencies),
						 tLong;
				if (popFunc>=0)
				{
					long var1 = LocateVarByName (modelName);
					if (var1>=0)
					{
						long var2 = LocateVarByName (modelOptions);
						if (var2>=0)
						{
							long var3 = LocateVarByName (modelDimension);
							if (var3>=0)
							{
								_Variable* v1 = FetchVar (var1),
										 * v2 = FetchVar (var2),
										 * v3 = FetchVar (var3);
								if ((v1->ObjectClass()==STRING)&&(v2->ObjectClass()==NUMBER)&&(v3->ObjectClass()==NUMBER))
								{
									_List thisList;
									
									modelParams.Clear();
									tLong = v2->Value();
									if (efvFunc==-1)
									{
										if (tLong&HY_DATAPANEL_MODEL_EFVEST)
											tLong -= HY_DATAPANEL_MODEL_EFVEST;
									}								
									modelParams << tLong;
									modelParams << (tLong=(long)v3->Value());
									if (((tLong==64)&&(codonFunc>=0))||(tLong!=64))
									{
									
										thisList && ((_FString*)v1->GetValue())->theString;
										thisList && & modelParams;
										thisList << receptacle(k);
										if (addedNames.Find (thisList(0))<0)
										{
											modelTemplates && & thisList;
											addedNames << thisList (0);
										}
									}
								}
							}
						}
					}
				}
			}
			while (g<batchLanguageFunctionNames.lLength)
			{
				batchLanguageFunctionNames.Delete (g);
				batchLanguageFunctionParameters.Delete (g);
				batchLanguageFunctions.Delete(g);
				batchLanguageFunctionClassification.Delete(g);
				batchLanguageFunctionParameterLists.Delete(g);
			}
		}
	}
	if (modelTemplates.lLength>1)
	{
		_String status = _String("Loaded ") & (long)modelTemplates.lLength & " model templates from " & pathToModelTemplates.getStr() & '\n';		
		StringToConsole (status);
	}
	else
	{
		_String warnMsg ("I couldn't find any valid model templates in '");
		warnMsg = warnMsg & pathToModelTemplates & "'. Please check your installation of HYPHY for missing files.";
		StringToConsole (warnMsg);
	}
}

//____________________________________________________________________________________________

void  NewGeneticCodeTable (long starting)
{
	static _String 		line1, line2;
	static _SimpleList	translationTable;

	if ((starting>=0)&&(starting<geneticCodes.lLength))
	{
		_List * startingTable = (_List*)geneticCodes(starting);
		line1 = *(_String*)(*startingTable)(0);
		line2 = *(_String*)(*startingTable)(1);
		translationTable.Duplicate ((*startingTable)(2));
	}
	else
		if (translationTable.lLength==0)
		{
			if (geneticCodes.lLength==0)
				ReadGeneticCodes();
			NewGeneticCodeTable (0);
			return;
		}

	_List	theList,protein;
	_String as, nucs("ACGT");
	as = "Phenylalanine";
	protein && & as;
	as = "Leucine";
	protein && & as;
	as = "Isoleucine";
	protein && & as;
	as = "Methionine";
	protein && & as;
	as = "Valine";
	protein && & as;
	as = "Serine";
	protein && & as;
	as = "Proline";
	protein && & as;
	as = "Threonine";
	protein && & as;
	as = "Alanine";
	protein && & as;
	as = "Tyrosine";
	protein && & as;
	as = "Stop Codon";
	protein && & as;
	as = "Histidine";
	protein && & as;
	as = "Glutamine";
	protein && & as;
	as = "Asparagine";
	protein && & as;
	as = "Lysine";
	protein && & as;
	as = "Aspartic Acid";
	protein && & as;
	as = "Glutamic Acid";
	protein && & as;
	as = "Cysteine";
	protein && & as;
	as = "Tryptophan";
	protein && & as;
	as = "Arginine";
	protein && & as;
	as = "Glycine";
	protein && & as;
	
	as = "AAA";
	AddItemToPreferences (1|8,-1,"Name and Description","Identifier and description for the genetic code.","",nil,theList,false);
	AddItemToPreferences (0,PREFITEM_TEXTBOX,"Code Identifier","This must be a valid unique HYPHY identifier (begins with a letter or an underscore, contains letters, numbers or underscores).",line1,nil, theList,false);
	AddItemToPreferences (0,PREFITEM_TEXTBOX,"Description","A brief description of the code table.",line2,nil,theList,false);
	AddItemToPreferences (1|8,-1,"Codon translations.","Define aminoacids encoded by each codon.","",nil,theList,false);

	for (long f = 0; f<4; f++)
	{
		as.sData[0] = nucs[f];
		for (long s = 0; s<4; s++)
		{
			as.sData[1] = nucs[s];
			for (long t = 0; t<4; t++)
			{
				as.sData[2] = nucs[t];
				AddItemToPreferences (0,PREFITEM_POPUP,as,_String("Target aminoacid for ")& as,
										*(_String*)protein(translationTable.lData[f*16+s*4+t]),&protein,theList,false);
			}
		}
	}
	
	if (HandlePreferences (theList, "New Translation Code", false))
	{
		_List *codeSettings = (_List*)theList.lData[4];
					  
		line1 = *(_String*)(*codeSettings)(1);
		long k;
		for (k=0; k<geneticCodes.lLength;k++)
		{
			if (line1.Equal((_String*)(*((_List*)geneticCodes(k)))(0)))
			{
				as = line1 & " is already in use. Please select another name for the new genetic code table.";
				if (ProceedPrompt(as))
				{
					NewGeneticCodeTable (-1);
				}
				return;
			}
		}
		
		line2 = *(_String*)(*codeSettings)(2);
		for (k = 4; k< 68 ; k++)
			translationTable.lData[k-4] = protein.Find((*codeSettings)(k));
			
		// add  the code to the table
		theList.Clear();
		theList&& &line1;
		theList&& &line2;
		theList&& &translationTable;	
		// save the code to a file
		
		#ifdef __MAC__
			as = baseDirectory & "GeneticCodes:";
		#endif
		#ifdef __WINDOZE__
			as = baseDirectory & "GeneticCodes\\";
		#endif
		
		as = as & line1;
		
		_String fName = as & ".cod";
		
		FILE* testMe;
		
		k = 2;
		
		while ((testMe=doFileOpen(fName.sData,"r")))
		{
			fclose (testMe);
			fName = as & _String (k++)& ".cod";
		}
		
		testMe = doFileOpen (fName.sData,"w");
		if (testMe)
		{
			fprintf(testMe,"%d;\n%s;\n%s;\n",geneticCodes.lLength,line1.sData,line2.sData);
			for (k=0; k<63; k++)
			{
				fprintf (testMe,"%d /*%c%c%c*/,\n",translationTable.lData[k],nucs.sData[k/16],nucs.sData[(k%16)/4],nucs.sData[k%4]);
			}
			fprintf (testMe,"%d /*TTT*/;\n",translationTable.lData[63]);
			fclose (testMe);
			geneticCodes && & theList;
		}
		else
		{
			as = "Sorry, but a file error foiled my attempts to save the new code table to disk.";
			ProblemReport (as);
		}
			
	}	
}

//____________________________________________________________________________________________

_List*	 FindModelTemplate (_String* modelName)
{
	for (long k=0; k<modelTemplates.lLength; k++)
	{
		_List* thisList = (_List*)modelTemplates(k);
		if (modelName->Equal((_String*)(*thisList)(0)))
			return thisList;
	}
	return nil;
}
//____________________________________________________________________________________________

_List*	 FindModelTemplate (long mID, long mSize)
{
	long m = 0;
	for (long k=0; k<modelTemplates.lLength; k++)
	{
		_List* thisList = (_List*)modelTemplates(k);
		if ((mSize==((_SimpleList*)(*thisList)(1))->lData[1])||
			(((((_SimpleList*)(*thisList)(1))->lData[1])==64)&&(mSize>45)))
		{
			if (mID==m)
				return thisList;
			m++;
		}
	}
	return nil;
}

//____________________________________________________________________________________________

long	 FindModelTemplate (_String* name, long mSize)
{
	long m = 0;
	for (long k=0; k<modelTemplates.lLength; k++)
	{
		_List* thisList = (_List*)modelTemplates(k);
		if ((mSize==((_SimpleList*)(*thisList)(1))->lData[1])||
			(((((_SimpleList*)(*thisList)(1))->lData[1])==64)&&(mSize>45)))
		{
			if (name->Equal((_String*)(*thisList)(0)))
				return m;
			m++;
		}
	}
	return -1;
}

//____________________________________________________________________________________________

long	 DimensionOfGenCode (long k)
{
	if ((k>=0)&&(k<geneticCodes.lLength))
	{
		_SimpleList * gC = (_SimpleList*)(*(_List*)geneticCodes(k))(2);
		long	dc =0;
		for (k=0; k<gC->lLength; k++)
			if (gC->lData[k]==10)
				dc++;
		return gC->lLength-dc;
	}
	return 0;
}

//____________________________________________________________________________________________

void	 SetModelMenus (int options, _HYPullDown * m1, _HYPullDown * m2)
{
	m1->EnableItem (1,options&HY_DATAPANEL_MODEL_GLOBAL);
	m1->EnableItem (2,options&HY_DATAPANEL_MODEL_GLOBALG);
	m2->EnableItem (3,options&HY_DATAPANEL_MODEL_EFVEST);
	if (options&HY_DATAPANEL_MODEL_MODELS)
	{
		m2->EnableItem (0,false);
		m2->EnableItem (1,false);
		m2->EnableItem (2,false);
		m2->EnableItem (3,false);
		m2->EnableItem (4,true);
		m1->EnableItem (0,!((options&HY_DATAPANEL_MODEL_GLOBAL)||(options&HY_DATAPANEL_MODEL_GLOBALG)));
	}
	else
	{
		m1->EnableItem (0,true);
		m2->EnableItem (0,true);
		m2->EnableItem (1,true);
		m2->EnableItem (2,true);
		m2->EnableItem (3,true);
		m2->EnableItem (4,false);
	}
}

//____________________________________________________________________________________________

bool	 RequestDataSetReplace (long k)
{
	for (long index = 0; index < windowObjectRefs.lLength; index++)
	{
		if (((_HYWindow*)windowObjectRefs(index))->WindowKind()==HY_WINDOW_KIND_DATAPANEL)
		{
			_HYDataPanel* thisPanel = (_HYDataPanel*)windowObjectRefs(index);
			if (thisPanel->GetDSID()==k)
			{
				_String errMsg ("Data Set '");
				errMsg = errMsg & *(_String*)dataSetNamesList (k) & "' can't be overwritten because it is in use by the window '" & 
						 thisPanel->GetTitle() & "'. Close this window, and try again.";
				ProblemReport (errMsg);
				return false;
			}
		}
	}
	return true;
}

//____________________________________________________________________________________________

bool	 RequestLFDeleteOrAlter (long k)
{
	for (long index = 0; index < windowObjectRefs.lLength; index++)
	{
		if (((_HYWindow*)windowObjectRefs(index))->WindowKind()==HY_WINDOW_KIND_DATAPANEL)
		{
			_HYDataPanel* thisPanel = (_HYDataPanel*)windowObjectRefs(index);
			if (thisPanel->GetLFID()==k)
			{
				_String errMsg ("The operation can't proceed because the likelihood function about to be affected is in use by the window '");
				errMsg = errMsg & thisPanel->GetTitle() & "'.";
				ProblemReport (errMsg);
				return false;
			}
		}
	}
	return true;
}

//__________________________________________________________

void ReadDataPanelProcessors (void)
{
	_String 	pathToModelTemplates;
	_List		receptacle;

	pathToModelTemplates	 = baseDirectory&"DatapanelAddIns";
	ScanDirectoryForFileNames (pathToModelTemplates,receptacle,false);
	
	for (long k=0; k<receptacle.lLength;k++)
	{
		FILE * thisFile = doFileOpen (((_String*)receptacle(k))->sData,"rb");
		if (thisFile)
		{	
			_String buffer (thisFile);
			fclose (thisFile);
			if (buffer.sLength)
				dataPanelProcessors << (_String*)receptacle(k);
		}
	}
}



//EOF