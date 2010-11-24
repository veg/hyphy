/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2009
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon    		   (apoon@cfenet.ubc.ca)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#include 	  "likefunc.h"
#include	  "scfg.h"

#if defined __AFYP_REWRITE_BGM__
#include	  "bayesgraph.h"
#else
#include 	  "bgm.h"
#endif

#if defined		  __UNIX__ && !defined __HEADLESS__ && !defined __HYPHY_NO_SQLITE__
	#include "SQLite/sqlite3.h"
#else
	#if defined  __HEADLESS__ && !defined __HYPHY_NO_SQLITE__
		#include 	  "sqlite3.h"	
	#else
		#ifndef __HYPHY_NO_SQLITE__
			#include 	  "sqlite3.h"
		#endif
		#if !defined		  __UNIX__  && !defined __HEADLESS__
			#include 	  "HYUtils.h"
		#endif
	#endif
#endif

#ifdef 	  __HYPHYDMALLOC__
	#include "dmalloc.h"
#endif



//____________________________________________________________________________________	
// global variables

_String		sqlOpen 				("SQL_OPEN"),
			sqlClose				("SQL_CLOSE"),
			sqlRowData				("SQL_ROW_DATA"),
			sqlColNames				("SQL_COLUMN_NAMES"),
			seqAlignMap				("SEQ_ALIGN_CHARACTER_MAP"),
			seqAlignScore			("SEQ_ALIGN_SCORE_MATRIX"),
			seqAlignGapChar			("SEQ_ALIGN_GAP_CHARACTER"),
			seqAlignGapOpen			("SEQ_ALIGN_GAP_OPEN"),
			seqAlignGapExtend		("SEQ_ALIGN_GAP_EXTEND"),
			seqAlignGapOpen2		("SEQ_ALIGN_GAP_OPEN2"),
			seqAlignGapExtend2		("SEQ_ALIGN_GAP_EXTEND2"),
			seqAlignGapLocal		("SEQ_ALIGN_NO_TP"),
			seqAlignGapAffine		("SEQ_ALIGN_AFFINE"),
			seqAlignGapLinearSpace	("SEQ_ALIGN_LINEAR_SPACE"),
			completeFlag 			("COMPLETE"),
			conditionalWeights		("WEIGHTS"),
			siteProbabilities		("SITE_LOG_LIKELIHOODS"),
			lastSetOfConstraints	("LAST_SET_OF_CONSTRAINTS"),
			deferConstrainAssignment("DEFER_CONSTRAINT_APPLICATION"),
			assertionBehavior		("ASSERTION_BEHAVIOR"),
			_hyStatusConditionProbsMatrix				
									("Constructing Conditional Probabilities Matrix"),

			isDynamicGraph			("BGM_DYNAMIC");
			
			
extern		_String					blDoSQL,
									blAlignSequences,
									blGetNeutralNull,
									blHBLProfile,
									blDeleteObject,
									timeStamp,
									versionString,
									lastModelParameterList,
									blGetString,
									blRequireVersion,
									blAssert;

_SimpleList sqlDatabases;

_List		scfgList,
			scfgNamesList,
			bgmList,
			bgmNamesList,
			_HY_GetStringGlobalTypesAux;

_AVLListX	 _HY_GetStringGlobalTypes (&_HY_GetStringGlobalTypesAux);



//____________________________________________________________________________________	

int  	 	 _HYSQLCallBack						(void* data,int callCount);
//int  	 	 _HYSQLBusyCallBack					(void* exList,int cc,char** rd,char** cn);
_Parameter	 AlignStrings						(_String*,_String*,_SimpleList&,_Matrix*,char,_Parameter,_Parameter,_Parameter,_Parameter,bool,bool,bool,_List&);
_Parameter	 CostOnly							(_String*,_String*, long, long, long, long, bool, bool, _SimpleList&, _Matrix*, _Parameter, _Parameter, _Parameter, _Parameter, bool, bool,_Matrix&, _Matrix*, _Matrix*, char = 0, char* = nil);
_Parameter   LinearSpaceAlign					(_String*,_String*, _SimpleList&, _Matrix*, _Parameter, _Parameter, _Parameter, _Parameter, bool, bool, _SimpleList&,_Parameter,long,long,long,long,_Matrix**, char, char*);
void 		 BacktrackAlign						(_SimpleList&, long&, long&, _Parameter, _Parameter, _Parameter);
void 		 MismatchScore						(_String*, _String*, long, long, _SimpleList&, _Matrix*, _Parameter&);
_String		 ProcessStringArgument				(_String*);
bool		 RecurseDownTheTree					(_SimpleList&, _List&, _List&, _List&, _SimpleList&);

//____________________________________________________________________________________	

void		_HBL_Init_Const_Arrays  (void)
{
		// init GetString lookups
		_HY_GetStringGlobalTypes.Insert(new _String("LikelihoodFunction"), 0);
		_HY_GetStringGlobalTypes.Insert(new _String("DataSet"), 1);
		_HY_GetStringGlobalTypes.Insert(new _String("DataSetFilter"), 2);
		_HY_GetStringGlobalTypes.Insert(new _String("UserFunction"), 3);
		_HY_GetStringGlobalTypes.Insert(new _String("Tree"), 4);
		_HY_GetStringGlobalTypes.Insert(new _String("SCFG"), 5);
		_HY_GetStringGlobalTypes.Insert(new _String("Variable"), 6);
}

//____________________________________________________________________________________	
void		 InsertVarIDsInList		(_AssociativeList* theList , _String theKey, _SimpleList& varIDs)
{
	_FString arrayKey (theKey, false);
	_Matrix *mxEntry = nil;
	
	if (varIDs.lLength)
	{
		_List     varNames;
		for (long i=0; i < varIDs.lLength; i++)
		{
			_Variable* v = LocateVar (varIDs.lData[i]);
			if (v)
				varNames << v->GetName();
		}
		mxEntry = new _Matrix (varNames);
	}
	else
		mxEntry = new _Matrix;
		
	checkPointer (mxEntry);
	theList->MStore (&arrayKey,mxEntry,false);
}

//____________________________________________________________________________________	
void		 InsertStringListIntoAVL	(_AssociativeList* theList , _String theKey, _SimpleList& stringsToPick, _List& theStrings)
{
	_FString arrayKey (theKey, false);
	_Matrix *mxEntry = nil;
	
	if (stringsToPick.lLength)
	{
		_List     theNames;
		for (long i=0; i < stringsToPick.lLength; i++)
		{
			_String * v = (_String*)theStrings (stringsToPick.lData[i]);
			if (v)
				theNames << v;
		}
		mxEntry = new _Matrix (theNames);
	}
	else
		mxEntry = new _Matrix;
		
	checkPointer (mxEntry);
	theList->MStore (&arrayKey,mxEntry,false);
}


//____________________________________________________________________________________	

_Matrix *   CheckMatrixArg 			(_String* mxName, bool onlyStrings)
{
	_Variable * mVar = FetchVar (LocateVarByName (*mxName));
	if (mVar && mVar->ObjectClass() == MATRIX)
	{
		_Matrix * mx = (_Matrix*)mVar->GetValue();
		if (onlyStrings && (!mx->IsAStringMatrix()))
			return nil;
		return mx;
	}
	return nil;
}

//____________________________________________________________________________________	

_AssociativeList *   CheckAssociativeListArg (_String* mxName)
{
	_Variable * mVar = FetchVar (LocateVarByName (*mxName));
	if (mVar && mVar->ObjectClass() == ASSOCIATIVE_LIST)
		return (_AssociativeList*)mVar->GetValue();
	return nil;
}

//____________________________________________________________________________________	

void	  _ElementaryCommand::ExecuteCase17 (_ExecutionList& chain)
{
	chain.currentCommand++;
	_String errMsg;
	if (parameters.lLength == 2)
	{
		_FString		* outLF = new _FString (new _String (8192L,1));	
		checkPointer	(outLF);
		_String			objectID (chain.AddNameSpaceToID(*(_String*)parameters(1)));
		_LikelihoodFunction * lf = FindLikeFuncByName (objectID);
		if (!lf)
		{
			long modelSpec = FindModelName (objectID);
			
			if (modelSpec<0)
			{
				long modelSpec = FindDataSetFilterName (objectID);
				if (modelSpec < 0)
				{
					WarnError (objectID & " does not refer to an existing likelihood function/model specification");
					outLF->theString->Finalize();
					DeleteObject (outLF);
					return ;					
				}
				else
				{
					outLF->theString->Finalize(); 
					DeleteObject (outLF->theString);
					checkPointer (outLF->theString = new _String ((_String*)((_DataSetFilter*)dataSetFilterList(modelSpec))->toStr()));
				}
			}
			else
			{
				SerializeModel (*outLF->theString,modelSpec,nil,true);
				outLF->theString->Finalize();
			}
		}
		else
		{
			lf->SerializeLF (*outLF->theString);
			outLF->theString->Finalize();
		}
		objectID = chain.AddNameSpaceToID(*(_String*)parameters(0));
		CheckReceptacleAndStore (&objectID, "Export", true, outLF, false);
	}
	else
	{
		_Matrix* m[2];
		for (long k=1;k<3;k++)
			if ((m[k-1] = (_Matrix*)FetchObjectFromVariableByType ((_String*)parameters(k),MATRIX)) == nil)
			{
				errMsg =  _String("Identifier ")&*(_String*)parameters(k)&" does not refer to a valid matrix variable";
				acknError (errMsg);
				return;
			}
			
		_String fName (*(_String*)parameters(0));
		fName.ProcessFileName();
		if (terminateExecution) return;
		FILE* 	theDump = doFileOpen (fName.getStr(),"w");
		if (!theDump)
		{
			WarnError (((_String)("File ")& fName &_String(" couldn't be open for writing.")));
			return;
		}
		m[1]->ExportMatrixExp(m[0],theDump);
		fclose (theDump);
	}
}

//____________________________________________________________________________________	

bool	_ElementaryCommand::ConstructDoSQL (_String&source, _ExecutionList&target)
// syntax: DoSQL (dbID,action string|file name,<callback ID>)
{
	_List pieces;
	ExtractConditions (source,blDoSQL.sLength,pieces,',');
	if (pieces.lLength!=3)
	{
		WarnError (_String ("Expected syntax:")& blDoSQL &"(dbID|" & sqlOpen & '|' & sqlClose & ",transaction string|file name,callback ID for an SQL transaction|where to store DB numeric ID)");
		return false;
	}
										
	_ElementaryCommand * dsql = new _ElementaryCommand (53);
	dsql->addAndClean(target,&pieces,0);
	return true;
}

//____________________________________________________________________________________	

bool	_ElementaryCommand::ConstructProfileStatement (_String&source, _ExecutionList&target)
// syntax: #profile START|PAUSE|RESUME|indetifier to dump in
{
	
	_List pieces;
	ExtractConditions (source,blHBLProfile.sLength+1,pieces,';');
	if (pieces.lLength!=2)
	{
		WarnError (_String ("Expected syntax:")& blHBLProfile &" START|PAUSE|RESUME|where to store)");
		return false;
	}
										
	_ElementaryCommand *sp = new _ElementaryCommand (58);
	sp->addAndClean(target,&pieces,0);
	return true;
}

//____________________________________________________________________________________	

bool	_ElementaryCommand::ConstructDeleteObject (_String&source, _ExecutionList&target)
// syntax: DeleteObject (object);
{
	_List pieces;
	ExtractConditions (source,blDeleteObject.sLength,pieces,';');
	if (pieces.lLength<1)
	{
		WarnError (_String ("Expected syntax:")& blDeleteObject &"comma separated list of objects to delete)");
		return false;
	}
	_ElementaryCommand *dobj = new _ElementaryCommand (59);
	dobj->addAndClean(target,&pieces,0);	
	return true;
}

//____________________________________________________________________________________	

bool	_ElementaryCommand::ConstructRequireVersion (_String&source, _ExecutionList&target)
// syntax: RequireVersion (stringObject);
{
	
	_List pieces;
	ExtractConditions (source,blRequireVersion.sLength,pieces,';');
	
	if (pieces.lLength != 1)
	{
		WarnError (_String ("Expected syntax:")& blRequireVersion &"a string object with HyPhy version)");
		return false;
	}
										
	_ElementaryCommand *rv = new _ElementaryCommand (60);
	rv->addAndClean(target,&pieces,0);	
	
	return true;
}



//____________________________________________________________________________________	

int  _HYSQLCallBack (void* exL,int cc, char** rd, char** cn)
{
	_ExecutionList * exList = (_ExecutionList *)exL;
	
	if (!terminateExecution)
		if (exList && cc && exList->lLength)
		{	
			_List	  rowData,
					  columnNames;
					  
			for (long cnt = 0; cnt < cc; cnt++)
			{
				if (rd[cnt])
					rowData.AppendNewInstance (new _String (rd[cnt]));
				else
					rowData.AppendNewInstance (new _String);

				if (cn[cnt])
					columnNames.AppendNewInstance (new _String (cn[cnt]));
				else
					columnNames.AppendNewInstance (new _String);
			}
					  
					  
			_Matrix * rowDataM     = new _Matrix (rowData),
					* columnNamesM = new _Matrix (columnNames);
					
			if (!(rowDataM && columnNamesM))
				checkPointer (nil);
			
			_Variable* rdv = CheckReceptacle (&sqlRowData, blDoSQL,false),
					 * cnv = CheckReceptacle (&sqlColNames, blDoSQL,false);
					 
			rdv->SetValue (rowDataM,false);
			cnv->SetValue (columnNamesM,false);
			
			exList->Execute();
			
		}
	return 0; 
}

//____________________________________________________________________________________	
/*
int  _HYSQLBusyCallBack (void* data, int callCount)
{
	if (callCount > 100)
		return 0;
	
#ifdef __WINDOZE__
	Sleep  (1 + genrand_real2()*100.);
#else	
	usleep (100 + genrand_real2()*100000.);
#endif
	return 1;
}*/


//____________________________________________________________________________________	

void	  _ElementaryCommand::ExecuteDataFilterCases (_ExecutionList& chain)
{
	chain.currentCommand++;
	
	_String dataObjectID = chain.AddNameSpaceToID(*(_String*)parameters(1));
	long dsID 			= (parameters.lLength>2)?FindDataSetName (dataObjectID):-1;
	bool isFilter		= false;
	
	if (dsID == -1)
	{
		dsID = (parameters.lLength>2)?FindDataSetFilterName (dataObjectID):-1;
		if (dsID == -1)
		{
			_AssociativeList * numericFilter = (_AssociativeList*)FetchObjectFromVariableByType(&dataObjectID, ASSOCIATIVE_LIST);
			if (numericFilter)
			{
				_String errCode;
				
				long	categoryCount = 1;
				
				if (parameters.lLength > 2) 
				// have multiple categories
					categoryCount = (long) ProcessNumericArgument((_String*)parameters(2),nil);
				
				_String namesKey ("FILTER_NAMES"),
						dataKey	 ("FILTER_ARRAYS"),
						freqKey	 ("FILTER_FREQS");
					
				_Matrix* sequenceNames = (_Matrix*)numericFilter->GetByKey (namesKey,MATRIX);
				_List seqNames;
				if (sequenceNames)
					sequenceNames->FillInList (seqNames);
				if (!sequenceNames || seqNames.lLength == 0)
					errCode = _String("Expected a non-empty string matrix as the ") & namesKey & " argument in call to CreateFilter";
				else
				{
					_AssociativeList * dataList = (_AssociativeList*)numericFilter->GetByKey (dataKey,ASSOCIATIVE_LIST);
					_Matrix			 * freqList = (_Matrix*)numericFilter->GetByKey (freqKey,MATRIX);
									 
					if (dataList && freqList)
					{
						_List 		goodSeqs;
						long 		sitePatterns 	= freqList->GetVDim(),
							  		categDim 		= -1;
							  		
						if (freqList->GetHDim() != 1 || sitePatterns < 1 || freqList->MatrixType() != 1 )
							errCode = _String("Expected a non-empty numeric ROW matrix as the ") & freqKey & " argument in call to CreateFilter";
						else
						{
							for (long k=0; k<seqNames.lLength; k=k+1)
							{
								_Matrix * dataMx = (_Matrix*)dataList->GetByKey (k,MATRIX);
								if (dataMx && dataMx->MatrixType() == 1 )
								{
									 if (categDim < 0)
									 {
									 	categDim = dataMx->GetVDim();
									 	if (categDim < 1)
									 		break;
									 }
									 else
									 	if (dataMx->GetVDim() != categDim)
									 		break;
									 if (dataMx->GetHDim () != sitePatterns*categoryCount)
									 	break;
									 
									 goodSeqs << dataMx;
									 continue;
								}
								break;
							}
							
							if (goodSeqs.lLength == seqNames.lLength)
							{
								_DataSet * dummyDS = new _DataSet;
								dummyDS->SetNoSpecies (seqNames.lLength);
								dummyDS->GetNames().Duplicate (&seqNames);
								dummyDS->GetTheMap().Populate (sitePatterns,0,1);
								errCode = (*(_String*)parameters(0)) & "_internal_ds";
								dsID = FindDataSetName (errCode);
								if (dsID < 0)
								{
									dataSetList			<<dummyDS;
									DeleteObject		(dummyDS);
									dataSetNamesList&&  &errCode;
								}
								else
									dataSetList.Replace (dsID,dummyDS,false);
								
								errCode = (*(_String*)parameters(0));
								_DataSetFilterNumeric * dsn = new _DataSetFilterNumeric (freqList,goodSeqs,dummyDS,categoryCount);
								checkPointer (dsn);
								dsID 	= FindDataSetFilterName (errCode);
									
								if (dsID < 0)
								{
									dataSetFilterList<<	 dsn;
									DeleteObject 		(dsn);
									dataSetFilterNamesList&& & errCode;
								}
								else
									dataSetFilterList.Replace (dsID,dsn,false);
								return; 
								
							}
							else
								errCode = _String ("Site frequency patterns/numeric vectors did not pass dimension checks in call to CreateFilter");
						}
					}
								
				}
				if (errCode)
				{
					WarnError(errCode);
					return;
				}
				
			}
			_String errMsg = (((_String)("DataSet(Filter)/Associative Array ")&dataObjectID&_String(" has not been properly initialized")));
			WarnError (errMsg);
			return;
		}
		isFilter = true;
	}
	
	// build the formula from the 2nd parameter (unit size)
	
	char	  			unit = ProcessNumericArgument((_String*)parameters(2),chain.nameSpacePrefix); 
					// here's our unit
	
	_String  			dataFilterID (chain.AddNameSpaceToID(*(_String*)parameters(0))),
						hSpecs, 
						vSpecs;
	
	long				status 	= FindDataSetFilterName (dataFilterID);
	
	_DataSetFilter		*thedf;
	
	if (status!=-1)
		thedf = (_DataSetFilter*)dataSetFilterList (status);
	else
	{
		thedf				= new _DataSetFilter();
		checkPointer		(thedf);
		AddFilterToList		(dataFilterID,thedf,false);
	}
	
	if (parameters.lLength>3)
		vSpecs = *(_String*)parameters(3);
	if (parameters.lLength>4)
		hSpecs = *(_String*)parameters(4);
	else
		hSpecs = empty;

	_DataSet 			*dataset;
	
	_SimpleList 		hL, 
						vL;
						
	hL.RequestSpace (1024);
	vL.RequestSpace (1024);
	
	if (!isFilter)
	{
		dataset = (_DataSet*)dataSetList(dsID);
		dataset -> ProcessPartition (hSpecs,hL,false);
		if (code!=6 && vSpecs.sLength==0)
				vSpecs = _String("0-")&_String(dataset->NoOfColumns()-1);
		dataset->ProcessPartition (vSpecs,vL,true);
	}
	else
	{
		_DataSetFilter * dataset1 = (_DataSetFilter*)dataSetFilterList(dsID);
		dataset1->GetData()->ProcessPartition (hSpecs,hL,false, &dataset1->theNodeMap, &dataset1->theOriginalOrder);
		
		if (code!=6 && vSpecs.sLength==0)
			vSpecs = _String("0-")&_String(dataset1->GetFullLengthSpecies()-1);
		
		dataset1->GetData()->ProcessPartition (vSpecs,vL,true,  &dataset1->theOriginalOrder, &dataset1->theNodeMap);
		dataset = (_DataSet*)dataset1;
	}
	
	if (code!=6)
	{
		if (vL.lLength%unit)
		{
		   vSpecs = (_String)"Unit size of "& unit & " doesn't divide the length of specified partition in call to ";
		   if (code==27) // Permute
		   		vSpecs = vSpecs & "Permute";
		   else
		   		vSpecs = vSpecs & "Bootstrap";

		   vSpecs = vSpecs & ". The partition has been trimmed at the end.";
		   ReportWarning (vSpecs);
		   for (status = vL.lLength%unit; status>0; status--)
		   		vL.Delete (vL.lLength-1);
		}
		if (code == 27)
			vL.Permute (unit);
		else
			vL.PermuteWithReplacement(unit);
	
	}
	
	thedf->SetFilter (dataset, unit, hL, vL, isFilter);
	
	if (parameters.lLength>5)
	{
		hSpecs = GetStringFromFormula((_String*)parameters(5),chain.nameSpacePrefix);
		thedf->SetExclusions(&hSpecs);
	}
	else
		if ((code!=6)&&isFilter)
		{
			_DataSetFilter * df1 = (_DataSetFilter*)dataSetFilterList(dsID);
			if (df1->theExclusions.lLength)
			{
				thedf->theExclusions.Duplicate (&df1->theExclusions);
				thedf->SetDimensions();
			}
		}
	
	thedf->SetDimensions();
	thedf->SetupConversion();
	
	SetDataFilterParameters (dataFilterID, thedf, true);							
}
//____________________________________________________________________________________	

void	  _ElementaryCommand::ExecuteCase21 (_ExecutionList& chain)
{
	chain.currentCommand++;

	SetStatusLine (_hyStatusConditionProbsMatrix);
	_String	errMsg,
			objectName 	  = 	chain.AddNameSpaceToID(*(_String*)parameters(1)),
			resultID      =     chain.AddNameSpaceToID(*(_String*)parameters(0));

	long objectID 		  = 	FindLikeFuncName (objectName, true);
	_PMathObj ob  		  = 	nil;
	
	if (objectID >=0) // likelihood function
	{
		_Matrix * partitionList			= nil;
		if (parameters.lLength>3)
		{
			_String  secondArg = *(_String*)parameters(3);
			partitionList = (_Matrix*)ProcessAnArgumentByType (&secondArg, chain.nameSpacePrefix, MATRIX);
		}
		_SimpleList						partsToDo;
		_LikelihoodFunction*			lf = (_LikelihoodFunction*)likeFuncList(objectID);
		if (lf->ProcessPartitionList(partsToDo, partitionList, _hyStatusConditionProbsMatrix))
		{
			char runMode = _hyphyLFConstructCategoryMatrixConditionals;
			if (parameters.lLength > 2)
			{
				if (((_String*)parameters(2))->Equal(&completeFlag))
					runMode = _hyphyLFConstructCategoryMatrixConditionals;
				else
					if (((_String*)parameters(2))->Equal(&conditionalWeights))
						runMode = _hyphyLFConstructCategoryMatrixWeights;
					else 
						if (((_String*)parameters(2))->Equal(&siteProbabilities))
							runMode = _hyphyLFConstructCategoryMatrixSiteProbabilities;
						else
							runMode = _hyphyLFConstructCategoryMatrixClasses;
			}
			ob = lf->ConstructCategoryMatrix(partsToDo,runMode,true, &resultID);
		}
	}
	else
	{
		_TheTree * testTree = (_TheTree*) FetchObjectFromVariableByType (&objectName, TREE);
		if (testTree)
		{
			long 	pid = 0;
			objectID = testTree->IsLinkedToALF (pid);
			if (objectID >= 0)
			{
				_LikelihoodFunction * anLF 		= (_LikelihoodFunction*) likeFuncList (objectID);
				_DataSetFilter		* dsf  		= (_DataSetFilter*) dataSetFilterList (anLF->GetTheFilters()(pid));
				anLF->PrepareToCompute();
				anLF->Compute		  ();
				objectID						= dsf->NumberDistinctSites();
				
				_Matrix				*condMx 	= new _Matrix	(2*objectID*(testTree->GetLeafCount() 
																+ testTree->GetINodeCount()) * testTree->categoryCount, 
																 testTree->GetCodeBase(), 
																 false, true);
				_List				leafNames,
									inodeNames;
				
				testTree->DepthWiseT(true);
				
				while (testTree->currentNode)
				{
					_String		  * bs = new _String; 	
					testTree->GetNodeName  (testTree->currentNode, *bs);
					if (testTree->IsCurrentNodeATip())
						leafNames << bs;
					else
						inodeNames << bs;
					DeleteObject (bs);
					testTree->DepthWiseT(false);
				}		
				
				leafNames << inodeNames;
				
				_Matrix	 *nodeNames = new _Matrix (leafNames);
				
				for (long siteC = 0; siteC < objectID; siteC ++)				
					testTree->RecoverNodeSupportStates (dsf,siteC,siteC-1,*condMx);
					
				anLF->DoneComputing   ();
				_AssociativeList *retMe = new _AssociativeList;
				retMe->MStore ("Nodes",nodeNames,false);
				retMe->MStore ("Values",condMx,false);
				ob = retMe;
			}
		}
	}
	
	if (ob)
		CheckReceptacleAndStore (&resultID, blConstructCM, true, ob, false);
	else
		WarnError (objectName & " must be either a likelihood function or a tree variable tied to a likelihood function.");
		
}

//____________________________________________________________________________________	
//	GetString
void	  _ElementaryCommand::ExecuteCase33 (_ExecutionList& chain)
{
	chain.currentCommand++;
	// first check to see if matrix parameters here are valid
	
	_String *currentArgument = (_String*)parameters(0), 
			 errMsg, 
			 result;
			 
	long 	f,
			sID,
			sID2 = -1;

	_Variable * theReceptacle = CheckReceptacle (&AppendContainerName(*currentArgument,chain.nameSpacePrefix),blGetString,true);
	
	if (!theReceptacle)
		return;

	sID = ProcessNumericArgument ((_String*)parameters(2),chain.nameSpacePrefix);
	if (parameters.lLength>3)
		sID2 = ProcessNumericArgument ((_String*)parameters(3),chain.nameSpacePrefix);
		
	currentArgument = (_String*)parameters(1);
	
	_String temp = ProcessStringArgument (currentArgument);
	
	if (temp.sLength)
		currentArgument = &temp;
	
	f = _HY_GetStringGlobalTypes.Find(currentArgument);
	if (f >=0 ) f = _HY_GetStringGlobalTypes.GetXtra (f);
	
	switch (f)
	{
		
		case 0: // LikelihoodFunction
			if (sID<likeFuncNamesList.lLength)
				result = *(_String*)likeFuncNamesList(sID);
			break;
		case 1: // DataSet
			if (sID<dataSetNamesList.lLength)
				result = *(_String*)dataSetNamesList(sID);
			break;
		case 2: // DataSetFilter
			if (sID<dataSetFilterNamesList.lLength)
				result = *(_String*)dataSetFilterNamesList(sID);
			break;
		case 3: // UserFunction
			if (sID<batchLanguageFunctions.lLength)
			{
				_AssociativeList * resAVL = (_AssociativeList *)checkPointer(new _AssociativeList);
				resAVL->MStore ("ID", new _FString (*(_String*)batchLanguageFunctionNames(sID)), false); 
				resAVL->MStore ("Arguments", new _Matrix(*(_List*)batchLanguageFunctionParameterLists(sID)), false); 
				theReceptacle->SetValue (resAVL,false);
				return;
			}
		case 4: // Tree
			{
				long tc = 0;
				_SimpleList nts;
				long		rt,
				vi = variableNames.Traverser (nts, rt, variableNames.GetRoot());
				
				for (; vi >= 0; vi = variableNames.Traverser (nts, rt))
					if (FetchVar(variableNames.GetXtra (vi))->ObjectClass () == TREE)
					{
						if (tc==sID)
						{
							result = *(_String*)variableNames.Retrieve(vi);
							break;
						}
						else	
							tc++;
					}
			
				break;
			}
		case 5: // SCFG
			if (sID<scfgNamesList.lLength)
				result = *(_String*)scfgNamesList(sID);
			break;
		default: // everything else...
		{
			// decide what kind of object current argument represents
			_String nmspaced = AppendContainerName(*currentArgument,chain.nameSpacePrefix);
			f = FindDataSetName (nmspaced);
			if (f>=0) // DataSet variable
			{
				_DataSet* dataSetObject = (_DataSet*)dataSetList(f);
				if (sID>=0 && sID<dataSetObject->NoOfSpecies())
					result = *(_String*)(dataSetObject->GetNames())(sID);
				else
				{
					theReceptacle->SetValue (new _Matrix (dataSetObject->GetNames()), false);
					return;
				}
			}
			else
			{
				f = FindDataSetFilterName (nmspaced);
				if (f>=0) // DateSetFilter variable
				{
					_DataSetFilter* dataSetFilterObject = (_DataSetFilter*)dataSetFilterList(f);
					if (sID >=0 && sID<dataSetFilterObject->NumberSpecies())
						result = *(_String*)(dataSetFilterObject->GetData()->GetNames())(dataSetFilterObject->theNodeMap(sID));
					else
					{
						_List filterSeqNames;
						for (f=0; f<dataSetFilterObject->NumberSpecies(); f++)
							filterSeqNames << (dataSetFilterObject->GetData()->GetNames())(dataSetFilterObject->theNodeMap(f));
						theReceptacle->SetValue (new _Matrix (filterSeqNames), false);
						return;
					}
					
				}
				else // not a data set filter
				{
					f = likeFuncNamesList.Find (&nmspaced);
					long g = -1;
						 
					if (f<0) // not a regular likelihood function
					{
	#if defined __AFYP_REWRITE_BGM__
						f = bgmNamesList.Find (&nmspaced);
						if (f >= 0)	// it's a BGM, export node score cache
						{
							_BayesianGraphicalModel	* this_bgm		= (_BayesianGraphicalModel *) bgmList (f);
							
							switch (sID)
							{
								case 0:		// return associative list containing node score cache
								{
									_AssociativeList		* export_alist	= new _AssociativeList;
									
									if (this_bgm -> ExportCache (export_alist))
									{
										theReceptacle -> SetValue (export_alist, false);
									}
									else
									{
										errMsg = _String("Failed to export node score cache.");
										ReportWarning (errMsg);
									}
									
									return;
								}
								case 1:		// return associative list with network structure and parameters
								{
									this_bgm -> SerializeBGM (result);
									_FString * resObj = new _FString (result);
									checkPointer (resObj);
									theReceptacle->SetValue (resObj,false);
									return;
								}	
								default:
								{
									errMsg = _String ("Unrecognized index ") & sID & "in call to GetString on BGM object";
									WarnError (errMsg);
									return;
								}
							}
							
						}
						else
						{
							g = FindSCFGName (nmspaced);
							f = g;
						}
	#else
						g = FindSCFGName (nmspaced);
						f = g;
	#endif
					}
					if (f>=0) // LikelihoodFunction
					{
						_LikelihoodFunction *lf = (_LikelihoodFunction*)(g>=0?scfgList(f):likeFuncList(f));
						if (sID>=0)
						{
							if (sID<lf->GetIndependentVars().lLength)
								result = *(LocateVar(lf->GetIndependentVars().lData[sID])->GetName());
							else
							{
								if (sID<lf->GetIndependentVars().lLength+lf->GetDependentVars().lLength)
									result = *(LocateVar(lf->GetDependentVars().lData[sID-lf->GetIndependentVars().lLength])->GetName());
								else
								{
									errMsg = _String(sID) & " is too high an index in call to GetString for likelihood function ";
									errMsg = errMsg&*(_String*)likeFuncNamesList(f);
									ReportWarning(errMsg);
								}
							}
						}
						else
						// return the AVL with parameters
						// AVL will have 8 entries 
							// "Categories"
							// "Global Independent"
							// "Global Constrained"
							// "Local Independent"
							// "Local Constrained"
							// "Trees"
							// "Models"
							// "Base frequencies"
							// "Datafilters"
							// "Compute Template"
						{
							_AssociativeList * resList = new _AssociativeList;
							
							_SimpleList		    *vl,
												v1,
												v2,
												v3,
												v4,
												v5,
												v6;
							_List				modelList;
												
							InsertVarIDsInList (resList, "Categories", lf->GetCategoryVars ());
							
							vl = &lf->GetIndependentVars ();
							for (long k=0; k<vl->lLength; k++)
								if (LocateVar (vl->lData[k])->IsGlobal())
									v1 << vl->lData[k];
								else
									v3 << vl->lData[k];
							
							vl = &lf->GetDependentVars ();
							for (long m=0; m<vl->lLength; m++)
								if (LocateVar (vl->lData[m])->IsGlobal())
									v2 << vl->lData[m];
								else
									v4 << vl->lData[m];
								
								
							vl = &lf->GetTheTrees ();
							for (long n=0; n<vl->lLength; n++)
							{
								v5 << vl->lData[n];
								_SimpleList partModels;
								((_TheTree*)FetchVar (vl->lData[n]))->CompileListOfModels(partModels);
								if (partModels.lLength == 1)
									modelList << modelNames (partModels.lData[0]);
								else
									modelList.AppendNewInstance(new _String ("__MULTIPLE__"));
							}
							

							vl = &lf->GetTheFilters ();
							for (long p=0; p<vl->lLength; p++)
								v6 << vl->lData[p];

							InsertVarIDsInList (resList, "Global Independent", v1);
							InsertVarIDsInList (resList, "Global Constrained", v2);
							InsertVarIDsInList (resList, "Local Independent", v3);
							InsertVarIDsInList (resList, "Local Constrained", v4);
							InsertVarIDsInList (resList, "Trees", v5);
							InsertVarIDsInList (resList, "Base frequencies", lf->GetBaseFreqs());
							InsertStringListIntoAVL (resList, "Datafilters", v6, dataSetFilterNamesList);
							{
								_SimpleList indexer (modelList.lLength,0,1);
								InsertStringListIntoAVL		(resList, "Models", indexer, modelList);
							}
							
							_FString		aKey,
											ct;
							_Formula		*computeT = lf->HasComputingTemplate();
							*ct.theString = computeT?(_String*)computeT->toStr():new _String;
							*aKey.theString = "Compute Template";
							resList->MStore (&aKey,&ct,true);
							if (g>=0)
								((Scfg*)lf)->AddSCFGInfo (resList);
							
							theReceptacle->SetValue (resList,false);
							return;
						}
					}
					else
					{
						f = currentArgument->Equal(&lastModelParameterList)?lastMatrixDeclared:modelNames.Find(&nmspaced);

						if (f>=0)
							// a model
						{
							if (sID>=0)
							{
								_Variable*		theMx = LocateVar (modelMatrixIndices.lData[f]);
								
								if (sID2 < 0) // get the sID's parameter name
								{
									_SimpleList		modelP;
									_AVLList		modelPA (&modelP);
									theMx->ScanForVariables(modelPA,false);
									modelPA.ReorderList();
									if (sID<modelP.lLength)
										result = *LocateVar(modelP.lData[sID])->GetName();
									else
									{
										errMsg = _String(sID) & " is too high a parameter index in call to GetString for model ";
										errMsg = errMsg&*(_String*)modelNames(f);
										ReportWarning(errMsg);
									}		
								}
								else // get the formula for cell (sID, sID2)
								{
									_Formula * cellFla = ((_Matrix*)theMx->GetValue())->GetFormula (sID,sID2);
									if (cellFla)
										result.CopyDynamicString ((_String*)cellFla->toStr(),true);
								}

							}
							else
							{
								_Variable	* tV, * tV2;
								bool		 mByF;
								RetrieveModelComponents (f, tV, tV2, mByF);
								
								if (tV)
								{
									if (sID == -1) // branch length expression
										result = ((_Matrix*)tV->GetValue())->BranchLengthExpression((_Matrix*)tV2->GetValue(),mByF);
									else 
									/* 
										returns an AVL with keys
										"RATE_MATRIX" - the ID of the rate matrix
										"EQ_FREQS"	  - the ID of eq. freq. vector
										"MULT_BY_FREQ" - a 0/1 flag to determine which format the matrix is in.
									*/
									{
										_AssociativeList * resList = new _AssociativeList;
										resList->MStore ("RATE_MATRIX",new _FString(*tV->GetName()),false);
										resList->MStore ("EQ_FREQS",new _FString(*tV2->GetName()),true);
										resList->MStore ("MULT_BY_FREQ",new _Constant (mByF),false);
										theReceptacle->SetValue (resList,false);
										return;
									}
								}
								else
								{
									theReceptacle->SetValue (new _FString(), false);
									return;
								}
							}
						}
						else
						{
							if (currentArgument->Equal(&versionString))
							{
								if (sID > 1.5)
								#ifdef __HEADLESS__
									result = _String ("Library version ") & __KERNEL__VERSION__;
								#else
									#ifdef __MAC__
										result = _String("Macintosh ") & __KERNEL__VERSION__;
									#else
										#ifdef __WINDOZE__
											result = _String("Windows ") & __KERNEL__VERSION__;
										#else
											result = _String("Source ") & __KERNEL__VERSION__;
										#endif
									#endif
								#endif
										else
											if (sID > 0.5)
												result = GetVersionString();
											else
												result = __KERNEL__VERSION__;
							}
							else
								if (currentArgument->Equal(&timeStamp))
									result = GetTimeStamp(sID < 0.5);
								else
								{
									f = FindBFFunctionName (*currentArgument);
									if (f >= 0)
									{
										_AssociativeList * resAVL = (_AssociativeList *)checkPointer(new _AssociativeList);
										resAVL->MStore ("ID", new _FString (*(_String*)batchLanguageFunctionNames(f)), false); 
										resAVL->MStore ("Arguments", new _Matrix(*(_List*)batchLanguageFunctionParameterLists(f)), false); 
										resAVL->MStore("Body", new _FString (((_ExecutionList*)batchLanguageFunctions(f))->sourceText,false),false);
										theReceptacle->SetValue (resAVL,false);
										return;										
									}
									else
									{
										_Variable* theVar = FetchVar(LocateVarByName (*currentArgument));
										if (theVar)
										{
											_String* theStr = nil;
											if (theVar->IsIndependent())
												theStr = (_String*)theVar->toStr();
											else
											{
												if (sID < -1.)
													// list of variables
												{
													_SimpleList vL; 
													_AVLList	vAVL (&vL);
													theVar->ScanForVariables (vAVL, sID > -2.5);
													vAVL.ReorderList();
													_AssociativeList   * resL = (_AssociativeList *) checkPointer (new _AssociativeList);
													InsertVarIDsInList (resL, "DEPENDANCIES", vL);
													theReceptacle->SetValue (resL,false);
													return;
												}
																									
												else	// formula string
													theStr = (_String*)theVar->GetFormulaString ();
											}
												
											result.CopyDynamicString(theStr, true);
										}
										else
										{
											errMsg = *currentArgument & " is not a validly defined object in call to GetString";
											ReportWarning (errMsg);
										}
									}
								}
						}				
					}
				}
			}
		}
	}
	_FString * resObj = new _FString (result);
	checkPointer (resObj);
	theReceptacle->SetValue (resObj,false);
}

//____________________________________________________________________________________	

void	  _ElementaryCommand::ExecuteCase53 (_ExecutionList& chain)
{
	chain.currentCommand++;

	#ifdef __HYPHY_NO_SQLITE__
		_String errStr ("SQLite commands can not be used in a HyPhy version built with the __HYPHY_NO_SQLITE__ flag");
		WarnError (errStr);
	#else
	
	_String arg1  (*(_String*)parameters(0)); 
						
	char  * errMsg = nil;
	_String errStr;
			
	if (arg1.Equal (&sqlOpen))
	{
		_Variable * dbVar = CheckReceptacle ((_String*)parameters(2), blDoSQL);
		
		if (dbVar)
		{
			_String arg2 (*(_String*)parameters(1));
			arg2.ProcessFileName(true,true,(Ptr)chain.nameSpacePrefix);
			int errCode  = SQLITE_OK;
			sqlite3 *aDB = nil;
			#ifdef __HYPHYXCODE__
				errCode = sqlite3_open (DoMacToPOSIX(arg2).getStr(),&aDB);
			#else
				errCode = sqlite3_open (arg2.sData,&aDB);
			#endif			
			if (errCode == SQLITE_OK)
				errCode = sqlite3_exec(aDB, "SELECT COUNT(*) FROM sqlite_master", _HYSQLCallBack, nil, nil);
			if (errCode != SQLITE_OK)
			{
				WarnError (sqlite3_errmsg(aDB));
				sqlite3_close(aDB);
				return;
			}
			else
			{
				long f = sqlDatabases.Find (0);
				if (f<0)
				{
					f = sqlDatabases.lLength;
					sqlDatabases << (long)aDB;
				}
				else
					sqlDatabases.lData[f] = (long)aDB;
				
				sqlite3_busy_timeout (aDB, 5000);
				
				dbVar->SetValue (new _Constant (f), false);
			}
		}
	}
	else
	{
		bool doClose = 	arg1.Equal (&sqlClose);

		long dbIdx = ProcessNumericArgument (doClose?(_String*)parameters(2):&arg1,chain.nameSpacePrefix);
		
		if (dbIdx<0 || dbIdx >= sqlDatabases.lLength || sqlDatabases.lData[dbIdx] == 0)
		{
			errStr = _String(dbIdx) & " is an invalid database index";
		}
		else
		{
			if (doClose)
			{
				sqlite3_close ((sqlite3*)sqlDatabases.lData[dbIdx]);
				sqlDatabases.lData[dbIdx] = 0;
			}
			else
			{
				_String arg3 (ProcessLiteralArgument((_String*)parameters(2),chain.nameSpacePrefix));
				
				_ExecutionList sqlProcessor (arg3,chain.nameSpacePrefix?(chain.nameSpacePrefix->GetName()):nil);
				if (!terminateExecution)
				{
					_String arg2 (ProcessLiteralArgument ((_String*)parameters(1),chain.nameSpacePrefix));
					
					if (sqlite3_exec((sqlite3*)sqlDatabases.lData[dbIdx], arg2.sData, _HYSQLCallBack, (Ptr)&sqlProcessor, &errMsg) != SQLITE_OK)
					{
						WarnError (sqlite3_errmsg((sqlite3*)sqlDatabases.lData[dbIdx]));
						return;							
					}
				}
			}
		}
			
	}
		
	if (errStr.sLength)
	{
		errStr = errStr & " in call to DoSQL";
		WarnError (errStr);
	}
	
	#endif
}

//____________________________________________________________________________________	

void	  _ElementaryCommand::ExecuteCase54 (_ExecutionList& chain)
{
	chain.currentCommand++;
	
	SetStatusLine (_String("Constructing Topology ")&*(_String*)parameters(0));
	
	_String  *treeSpec = ((_String*)parameters(1));
	treeSpec->ProcessParameter();
	_TreeTopology * tr = nil;
	
	if (treeSpec->sLength)
	{
		if (treeSpec->sData[0]!='(')
		{
			_Variable* testTree = FetchVar(LocateVarByName (AppendContainerName(*treeSpec,chain.nameSpacePrefix)));
			if (testTree && testTree->ObjectClass () == TREE)
				tr = new _TreeTopology ((_TheTree*)testTree);		
			else
			{
				_String   flaData (*treeSpec);
				_Formula  nameForm (flaData,chain.nameSpacePrefix);
				_PMathObj formRes = nameForm.Compute();
				if (formRes&&formRes->ObjectClass () == STRING)
						tr = new _TreeTopology (AppendContainerName(*(_String*)parameters(0),chain.nameSpacePrefix),
												*((_FString*)formRes)->theString,
												false);
			}
		}
		else
			tr = new _TreeTopology (AppendContainerName(*(_String*)parameters(0),chain.nameSpacePrefix),
									*(_String*)parameters(1),false);
	}
		
	if (!tr)
		WarnError ("Illegal right hand side in call to Topology id = ...; it must be a string, a Newick tree spec or a topology");
}

//_________________________________________________________________________

bool	ExpressionCalculator (void)
{
	_String data (StringFromConsole(false));
	
	#ifndef __UNIX__
		if (terminateExecution)
			return false;
		BufferToConsole (">");
		StringToConsole (data);
		BufferToConsole ("\n");
	#endif
	
	if (data.sLength == 4)
	{
		_String checkForExit (data);
		checkForExit.LoCase();
		if (checkForExit == _String ("exit"))
			return false;
	}
	
	_Formula  lhs,
			  rhs;
	
	long	   refV,
			   retCode = Parse(&lhs, data, refV, nil,nil);
	
	if (!terminateExecution)
	{
		if (retCode == HY_FORMULA_EXPRESSION)
		{
			_PMathObj formRes = lhs.Compute();
			if (!formRes)
			{
				BufferToConsole ("NULL\n");
			}
			else
			{
				_String * objValue = (_String*)formRes->toStr();
				StringToConsole (*objValue);
				BufferToConsole ("\n");
				DeleteObject    (objValue);
			}
		}
		else
		{
			BufferToConsole ("NO RETURN VALUE\n");
		}
	}
	terminateExecution = false;
	return true;
}

//_________________________________________________________________________

bool	PushFilePath (_String& pName)
{
	char c = GetPlatformDirectoryChar();
	
	long	f = pName.FindBackwards(_String(c),0,-1);
	if (f>=0)
	{
		_String newP = pName.Cut(0,f);
		pathNames && & newP;
		pName.Trim (f+1,-1);
		return true;
	}
	else
		if (pathNames.lLength)
			pathNames && pathNames(pathNames.lLength-1);
		else
			pathNames && & empty;
		
	return false;
}

//_________________________________________________________________________
void   PopFilePath (void)
{
	pathNames.Delete (pathNames.lLength-1);
}

//_________________________________________________________________________
void   ExecuteBLString (_String& BLCommand, _VariableContainer* theP)
{
	_ExecutionList ex;
	if (theP)
		ex.SetNameSpace(*theP->GetName());
	ex.BuildList   (BLCommand);
	terminateExecution = false;
	ex.Execute	    ();
	terminateExecution = false;
}

//____________________________________________________________________________________	

bool	_ElementaryCommand::ConstructAlignSequences (_String&source, _ExecutionList&target)
// syntax: AlignSequences (result, input string matrix,  options matrix)
{
	_List pieces;
	ExtractConditions (source,blAlignSequences.sLength,pieces,',');
	if (pieces.lLength!=3)
	{
		WarnError ("Expected syntax: AlignSequences(result, input string matrix, options list);");
		return false;
	}
										
	_ElementaryCommand * as = new _ElementaryCommand (55);
	as->addAndClean(target,&pieces,0);
	return true;
}

//____________________________________________________________________________________	

bool	_ElementaryCommand::ConstructGetNeutralNull (_String&source, _ExecutionList&target)
// syntax: GetNeutralNull (result, likelihood function, syn sub count matrix, non-syn sub count matrix, iterations per root state)
{
	_List pieces;
	ExtractConditions (source,blGetNeutralNull.sLength,pieces,',');
	if (pieces.lLength!=5)
	{
		WarnError ("Expected syntax: GetNeutralNull (result, likelihood function, syn sub count matrix, non-syn sub count matrix, iterations per root state);");
		return false;
	}
										
	_ElementaryCommand * gnn = new _ElementaryCommand (57);
	gnn->addAndClean(target,&pieces,0);
	return true;
}

//____________________________________________________________________________________	

void	  _ElementaryCommand::ExecuteCase55 (_ExecutionList& chain)
{
	chain.currentCommand++;

	_String errStr;
			
	_Variable * storeResultIn = CheckReceptacle (&AppendContainerName(*(_String*)parameters(0),chain.nameSpacePrefix), blAlignSequences, true);
	
	if (storeResultIn)
	{
		_Matrix * inStrings = CheckMatrixArg (&AppendContainerName(*(_String*)parameters(1),chain.nameSpacePrefix),true);
		if (inStrings && (inStrings->GetHDim()==1||inStrings->GetVDim()==1))
		{
			_AssociativeList * mappingTable = CheckAssociativeListArg (&AppendContainerName(*(_String*)parameters(2),chain.nameSpacePrefix));
			if (mappingTable)
			{
				// check for required parameters
				_FString * charVector = (_FString*)mappingTable->GetByKey (seqAlignMap, STRING);
				
				long	   charCount = 0;
				
				_SimpleList ccount (256,-1,0);
				
				if (charVector)
				{
					for (long cc = 0; cc < charVector->theString->sLength; cc++)
						if (ccount.lData[charVector->theString->sData[cc]]>=0)
						{
							charCount = 0;
							break;
						}
						else
						{
							ccount.lData[charVector->theString->sData[cc]] = cc;
							charCount ++;
						}
				}
				if (charVector && charCount)
				{
					// now check that all characters
					_Matrix * scoreMatrix = (_Matrix*)mappingTable->GetByKey (seqAlignScore, MATRIX);
					if (scoreMatrix && scoreMatrix->GetHDim () == charCount && scoreMatrix->GetVDim () == charCount)
					{
						scoreMatrix = (_Matrix*)scoreMatrix->ComputeNumeric();
						scoreMatrix->CheckIfSparseEnough(true);
						
						char		gapCharacter = '-';
						_FString	*gapC = (_FString*)mappingTable->GetByKey (seqAlignGapChar, STRING);
						
						_String		settingReport (128L,true);
						
						settingReport << "Running sequence alignment with the following options:";
						
						if (gapC && gapC->theString->sLength == 1)
							gapCharacter = gapC->theString->sData[0];
							
						settingReport << "\n\tGap character:";
						settingReport << gapCharacter;

						_Parameter	gapOpen   = 15.,
									gapOpen2  = 15.,
									gapExtend = 1.,
									gapExtend2= 1.;
									
						bool		doLocal 	= false,
									doAffine	= false,
									doLinear	= true;
									
									
						_PMathObj 	c = mappingTable->GetByKey (seqAlignGapOpen, NUMBER);
						if (c)
							gapOpen = c->Compute()->Value();
						
						settingReport << "\n\tGap open cost:";
						settingReport << _String (gapOpen);
						
						gapOpen2 = gapOpen;
						c = mappingTable->GetByKey (seqAlignGapOpen2, NUMBER);
						if (c)
							gapOpen2 = c->Compute()->Value();
							
						settingReport << "\n\tGap open cost 2:";
						settingReport << _String (gapOpen2);

						c = mappingTable->GetByKey (seqAlignGapExtend, NUMBER);
						if (c)
							gapExtend = c->Compute()->Value();
							
						settingReport << "\n\tGap extend cost:";
						settingReport << _String (gapExtend);

						gapExtend2 = gapExtend;
						c = mappingTable->GetByKey (seqAlignGapExtend2, NUMBER);
						if (c)
							gapExtend2 = c->Compute()->Value();
							
						settingReport << "\n\tGap extend cost:";
						settingReport << _String (gapExtend2);

						c = mappingTable->GetByKey (seqAlignGapLocal, NUMBER);
						if (c)
							doLocal = c->Compute()->Value() > 0.5;

						settingReport << "\n\tIgnore terminal gaps: ";
						if (doLocal)
							settingReport << "Yes";
						else
							settingReport << "No";

						c = mappingTable->GetByKey (seqAlignGapAffine, NUMBER);
						if (c)
							doAffine = c->Compute()->Value() > 0.5;

						settingReport << "\n\tAffine gap costs: ";
						if (doAffine)
							settingReport << "Yes";
						else
							settingReport << "No";

						c = mappingTable->GetByKey (seqAlignGapLinearSpace, NUMBER);
						if (c)
							doLinear = c->Compute()->Value() > 0.5;
							
						settingReport << "\n\tUse linear space routines: ";
						if (doLinear)
							settingReport << "Yes";
						else
							settingReport << "No";
							
						settingReport.Finalize();
						ReportWarning (settingReport);

						long stringCount = inStrings->GetHDim() * inStrings->GetVDim();
						
						_AssociativeList *alignedStrings = new _AssociativeList;
						checkPointer (alignedStrings);
						
						
						for (long s1 = 0; s1 < stringCount; s1++)
						{
							_String*  str1 = ((_FString*)inStrings->GetFormula(0,s1)->Compute())->theString;
							if (!str1)
							{
								errStr = _String("The ") & (s1+1) & "-th argument is not a string";
								break;
							}
							for (long s2 = s1+1; s2 < stringCount; s2++)
							{
								_String		  *string2 = ((_FString*)inStrings->GetFormula(0,s2)->Compute())->theString;
								if (!string2)
								{
									errStr = _String("The ") & (s2+1) & "-th argument is not a string";
									break;
								}
								_AssociativeList * pairwiseComp = new _AssociativeList;
								checkPointer (pairwiseComp);
								
								_Parameter 	  score = 0.0;
								
								if (doLinear == false)
								{
									_List 		  store;
									score = AlignStrings (str1,string2,ccount,scoreMatrix,gapCharacter,
																	gapOpen,gapExtend,gapOpen2,gapExtend2,doLocal,doAffine,doLinear,store);
									store.bumpNInst();
									pairwiseComp->MStore ("1", new _FString((_String*)store(0)), false);
									pairwiseComp->MStore ("2", new _FString((_String*)store(1)), false);
									
									if (store.lLength == 0)
									{
										errStr = "Unspecified error in AlignStrings";
										DeleteObject (pairwiseComp);
										s1 = stringCount;
										break;
									}
								}
								else
								{
									_Matrix		  scoreM	    (string2->sLength+1,1,false,true),
												  scoreM2	    (string2->sLength+1,1,false,true),
												  gap1Matrix    (string2->sLength+1,1,false,true),
												  gap2Matrix    (string2->sLength+1,1,false,true),
												  gap1Matrix2   (string2->sLength+1,1,false,true),
												  gap2Matrix2   (string2->sLength+1,1,false,true),
												  *buffers[6];
									
									char		  *alignmentRoute = new char[2*(string2->sLength+1)];
									
									alignmentRoute[0] = alignmentRoute[string2->sLength+1] = 0;
									buffers[0] = &scoreM; buffers[1] = &gap1Matrix; buffers[2] = &gap2Matrix;
									buffers[3] = &scoreM2; buffers[4] = &gap1Matrix2; buffers[5] = &gap2Matrix2;
									_SimpleList ops (str1->sLength+2,-2,0);
									ops.lData[str1->sLength+1] = string2->sLength;
									ops.lData[0]			   = -1;
									
									score = LinearSpaceAlign(str1,string2,ccount,scoreMatrix,
																		 gapOpen,gapExtend,gapOpen2,gapExtend2,
																		 doLocal,doAffine,ops,score,0,
																		 str1->sLength,0,string2->sLength,buffers,0,alignmentRoute);
									
									delete[]	alignmentRoute;

									_String		*result1 = new _String (str1->sLength+1, true),
												*result2 = new _String (string2->sLength+1, true);
									
									long		last_column		= ops.lData[ops.lLength-1];
									
									for (long position = str1->sLength-1; position>=0; position--)
									{
										long current_column		= ops.lData[position+1];
										
										if (current_column<0)
										{
											if (current_column == -2 /*|| (current_column == -3 && last_column == string2->sLength)*/)
												current_column = last_column;
											else
												if (current_column == -3)
												{
													// find the next matched char or a -1
													long	p	= position,
															s2p;
													while ((ops.lData[p+1]) < -1) p--;
													
													s2p = ops.lData[p+1];
													//if (last_column == string2->sLength)
													//	last_column = string2->sLength-1;
													
													//if (s2p < 0)
													//	s2p = 0;
													
													for (long j = last_column-1; j>s2p;)
													{
														(*result1) << gapCharacter;
														(*result2) << string2->sData[j--];
													}
													
													last_column     = s2p+1;
													
													for (;position>p;position--)
													{
														(*result2) << gapCharacter;
														(*result1) << str1->sData[position];
													}
													position ++;
													continue;
												}
												else
												{
													for (last_column--; last_column >=0; last_column--)
													{
														(*result1) << gapCharacter;
														(*result2) << string2->sData[last_column];								
													}
													while (position>=0)
													{
														(*result1) << str1->sData[position--];
														(*result2) << gapCharacter;										
													}
													break;
												}
										}
										
										if (current_column == last_column) // insert in sequence 2
										{
											(*result1) << str1->sData[position];
											(*result2) << gapCharacter;
										}
										else
										{
											last_column--;
										
											for (; last_column > current_column; last_column--) // insert in column 1
											{
												(*result2) << string2->sData[last_column];
												(*result1) << gapCharacter;												
											}
											(*result1) << str1->sData[position];
											(*result2) << string2->sData[current_column];
										}
										//printf ("%s\n%s\n", result1->sData, result2->sData);
									}
									
									for (last_column--; last_column >=0; last_column--)
									{
										(*result1) << gapCharacter;
										(*result2) << string2->sData[last_column];								
									}
									
									result1->Finalize(); result1->Flip ();
									result2->Finalize(); result2->Flip ();
									pairwiseComp->MStore ("1", new _FString(result1), false);
									pairwiseComp->MStore ("2", new _FString(result2), false);
								}
								/*
								long gap1c = 0,
									 gap2c = 0;
								 
								 _Parameter scoreCheck = 0.;
								 
								 for (long sp = 0; sp<result1->sLength; sp++)
								 {
									 char cs1 = result1->sData[sp],
										  cs2 = result2->sData[sp];
									 
									 if (cs1 == gapCharacter)
									 {
										 if (gap1c && doAffine)
											 scoreCheck -= gapExtend;
										 else
											 scoreCheck -= gapOpen;
										 gap2c = 0;
										 gap1c++;
									 }
									 else
									 if (cs2 == gapCharacter)
									 {
										 if (gap2c && doAffine)
											 scoreCheck -= gapExtend2;
										 else
											 scoreCheck -= gapOpen2;
										 gap1c = 0;
										 gap2c++;
									 }
									 else
									 {
										 gap1c = 0;
										 gap2c = 0;
										 long code1 = ccount.lData[cs1],
											  code2 = ccount.lData[cs2];
									 
										 if (code1 >=0 && code2 >=0 )
											 scoreCheck += (*scoreMatrix)(code1,code2);
									 }
								 }
								 if (doLocal)
								 {
									for (long k = 0; result1->sData[k] == gapCharacter; k++)
										if (doAffine)
											scoreCheck += k?gapExtend:gapOpen;
										else
											scoreCheck += gapOpen;
									 for (long k = 0; result2->sData[k] == gapCharacter; k++)
										 if (doAffine)
											 scoreCheck += k?gapExtend2:gapOpen2;
										 else
											 scoreCheck += gapOpen2;
									 for (long k = result1->sLength-1; result1->sData[k] == gapCharacter; k--)
										 if (doAffine)
											 scoreCheck += k==result1->sLength-1?gapOpen:gapExtend;
										 else
											 scoreCheck += gapOpen;
									 for (long k = result2->sLength-1; result2->sData[k] == gapCharacter; k--)
										 if (doAffine)
											 scoreCheck += k==result2->sLength-1?gapOpen2:gapExtend2;
										 else
											 scoreCheck += gapOpen2;
								 }*/
								 						
								
								
								pairwiseComp->MStore ("0", new _Constant (score), false);
								/*pairwiseComp->MStore ("3", new _Constant (score2), false);
								pairwiseComp->MStore ("4", new _FString(result1), false);
								pairwiseComp->MStore ("5", new _FString(result2), false);
								pairwiseComp->MStore ("6", new _FString((_String*)ops.toStr()), false);
								pairwiseComp->MStore ("7", new _Constant (scoreCheck), false);*/
								alignedStrings->MStore (_String(s1), pairwiseComp, false);
							}
						}
						
						storeResultIn->SetValue (alignedStrings, false);

					}
					else
						errStr = seqAlignScore & " is a required option, which must be a square matrix with dimension matching the size of " & seqAlignMap;
				}
				else
					errStr = seqAlignMap & " is a required option, which must be a non-empty string without repeating characters ";
			}
			else
				errStr = *(_String*)parameters(2) & " was expected to be an associative array of alignment options";
		}
		else
			errStr = *(_String*)parameters(1) & " was expected to be a vector of strings";
	}
		
	if (errStr.sLength)
	{
		errStr = errStr & " in call to " & blAlignSequences;
		WarnError (errStr);
	}
}

//____________________________________________________________________________________	

bool	RecurseDownTheTree (_SimpleList& theNodes, _List& theNames, _List&theConstraints, _List& theParts, _SimpleList& partIndex)
{
	_SimpleList localNodes;
	
	node<long>* firstNode = (node<long>*)theNodes(0), *otherNode;
	bool		doThisOne = (firstNode->get_parent()!=nil), good = true;
	long 		index, ind, i;
	
	/*if (!doThisOne)
	{
		BufferToConsole (_String((_String*)theParts.toStr()));
		NLToConsole();
		BufferToConsole (_String((_String*)partIndex.toStr()));
		NLToConsole();
	}*/
	
	// there are a few cases to consider
	for (ind = 1; ind<=firstNode->get_num_nodes(); ind++)// have children nodes
	{
		localNodes<< (long)firstNode->go_down(ind);
		for (index = 1; index<theNodes.lLength; index++)
		{
			otherNode = (node<long>*)theNodes(index);
			otherNode = otherNode->go_down(ind);
			if (!otherNode)
			{
				good = false;
				break;
			}
			localNodes<<(long)otherNode;
		}
		if (!good) break;
		good = RecurseDownTheTree (localNodes, theNames, theConstraints, theParts, partIndex);
		if (!good) break;
		localNodes.Clear();
	}
	
	// do this constraint now	
	
	if (doThisOne&&good) // not a root - so we apply the constraint
	{
		_CalcNode*  firstCNode = (_CalcNode*)LocateVar (firstNode->get_data());
		_SimpleList goodVars;
		_List		otherGoodVars;
		_Variable* firstVar;
		
		ind = 0;
		
		while ((firstVar=firstCNode->GetIthIndependent(ind)))
		{
			for (index = 0; index<partIndex.lLength; index++)
			{
				if (partIndex.lData[index]==0)
				{
					if (!firstVar->GetName()->EqualWithWildChar((_String*)theParts.lData[index],'?'))
						break;
				}
			}
			if (index==partIndex.lLength)
				goodVars<<ind;
			ind++;
		}
		
		for (i = 1; i<theNodes.lLength; i++)
		{
			otherNode = (node<long>*)theNodes(i);
			firstCNode = (_CalcNode*)LocateVar (otherNode->get_data());
			_SimpleList dummy;
			otherGoodVars && & dummy;
			long          theseInd = firstCNode->CountAll();
			_SimpleList   avVars;
			for (index = 0; index < theseInd; index++)
				avVars << index;
			for (index = 0; index<goodVars.countitems(); index++)
			{
				long j=0,k=0;
				bool found1 = false;
				for (k = 0; k<partIndex.lLength; k++)
					if (partIndex.lData[k]==i) break;
				
				for (; j<avVars.lLength; j++)
				{
					firstVar = firstCNode->GetIthParameter(avVars.lData[j]);						
					if (firstVar->GetName()->EqualWithWildChar((_String*)theParts.lData[k],'?'))
					{
						(*(_SimpleList*)(otherGoodVars(i-1))) << avVars.lData[j];
						avVars.Delete (j);
						found1 = true;
						break;
					}
				}
				if (!found1)
				{
					goodVars.Delete (index);
					for (long ff = 0; ff < i-1; ff++)
						((_SimpleList*)(otherGoodVars(i-1)))->Delete (index);
					index--;
				}
			}
		} 
		
		// now the constraints can be built
			
		for (index = 0; index < goodVars.lLength; index++)
		{
			_String newConstraint;
			for (ind = 0; ind < partIndex.lLength; ind++)
			{
				if (partIndex.lData[ind]<0)
					newConstraint = newConstraint & *(_String*)theParts(ind);
				else
				{
					otherNode = (node<long>*)theNodes(partIndex.lData[ind]);
					_CalcNode*	CNode = (_CalcNode*)LocateVar (otherNode->get_data());

					if (ind>0)
						newConstraint = newConstraint & 
									*(CNode->GetIthParameter((*(_SimpleList*)
									(otherGoodVars(partIndex.lData[ind]-1))).lData[index])->GetName());
					else
						newConstraint = newConstraint & 
									*(CNode->GetIthIndependent(goodVars.lData[index])->GetName());
				}
			
			}
			theConstraints&& &newConstraint;
		}
	}
	
	if (!good)
	{
		_String errMsg (*(LocateVar(firstNode->get_data())->GetName())& " is incompatible with "&
					   (*LocateVar(((node<long>*)theNodes(index-1))->get_data())->GetName()) & " in call to ReplicateConstraint");
		WarnError (errMsg);
		return false;
	}
	 
	return true; 

}

//____________________________________________________________________________________	

void	  _ElementaryCommand::ExecuteCase26 (_ExecutionList& chain)
{
	chain.currentCommand++;
	// we have to build a list of _CalcNodes to deal with
	// all of the trees/nodes in ReplicateConstraint must be of the same topology
	// the constraint will be processed by trying all of the subnodes of the given node
	// and within each - trying all of the variables to see if the constraint is matched
	// exactly the same operation will be repeated on each of the other parameters
	
	_String * 		replicateSource, 
					thisS, 
					prStr = GetStringFromFormula((_String*)parameters(0),chain.nameSpacePrefix);
	
	replicateSource = &prStr;
	
	_List			parts, 
					theConstraints;
				 
	_SimpleList  	thisIndex, 
					thisArgs;
	
	long			ind1	=	replicateSource->Find("this"), 
					ind2, 
					ind3, 
					ind4;
			
	if (ind1<0)
	{
		WarnError (*(_String*)parameters(0)&" has no 'this' references in call to ReplicateConstraint!");
		return ;
	}
	
	_SimpleList thisHits (parameters.lLength-1,0,0);
	
	while (ind1>=0) // references to 'this' still exist
	{			
		ind2 = ind1+4; // look forward to the number of 'this'
		while ('0'<=replicateSource->sData[ind2] && replicateSource->sData[ind2]<='9')
			ind2++;

		ind3  = replicateSource->Cut(ind1+4,ind2-1).toNum();
		ind2  = replicateSource->FindEndOfIdent (ind1,-1,'?');
		// now ind1-ind2 contains a reference with this...
		_String	newS  (*replicateSource,ind1,ind2);
		thisS = _String("this")&_String(ind3);
		if ((ind4 = ((_String*)parameters(ind3))->Find('.'))>=0) // branch argument
			newS = newS.Replace (thisS,((_String*)parameters(ind3))->Cut(0,ind4-1), true);
		else // non-branch argument
			newS = newS.Replace (thisS,*((_String*)parameters(ind3)), true);
		parts&& &newS;
		ind3--;
		thisIndex<<ind3; // sequence of references to this
		
		if (ind3<0 || ind3 >= thisHits.lLength)
		{
			WarnError (_String("Invalid reference to ") & thisS & " in the constraint specification");
			return ;			
		}
		thisHits.lData[ind3] = 1;
		
		if (ind2>=replicateSource->sLength-1) 
			break;
		ind1 = replicateSource->Find("this",ind2+1,-1);
		if (ind1==-1)
			newS = replicateSource->Cut(ind2+1,-1);
		else
			newS = replicateSource->Cut(ind2+1,ind1-1);
		parts&& &newS;
		thisIndex<<-1;
	}
	// now that the string is conveniently partritioned into blocks 
	// we will check the arguments and store references
	
	for (ind1 = 1; ind1<parameters.lLength; ind1++)
	{ 
		if (thisHits.lData[ind1-1] == 0)
		{
			WarnError (_String("Unused ") & ind1 & "-th reference variable: " & *(_String*)parameters(ind1));
			return ;			
		}
								  
		ind2 = LocateVarByName (*(_String*)parameters(ind1));
		if (ind2<0)
		{
			_String newS = *(_String*)parameters(ind1) & " is undefined in call to ReplicateConstraint.";
			acknError (newS);
			return	;
		}
		
		_Variable* thisNode = FetchVar (ind2);
		if (thisNode->ObjectClass()==TREE_NODE)
			thisArgs<< (long)((_CalcNode*)thisNode)->LocateMeInTree();
		else
			if (thisNode->ObjectClass()==TREE)
				thisArgs<< (long)&((_TheTree*)thisNode)->GetRoot();
			else
			{
				WarnError (*(_String*)parameters(ind1) & " is neither a tree nor a tree node in call to ReplicateConstraint.");
				return ;
			}
	}
	
	// now with this list ready we can recurse down the tree and produce the contsraints
	if (RecurseDownTheTree(thisArgs, parameters, theConstraints, parts, thisIndex))
	{			
		if (theConstraints.lLength)
		{
			ReportWarning  (_String("\nReplicateConstraint generated the following contsraints:"));
			_Parameter      doDeferSet;
			checkParameter (deferConstrainAssignment,doDeferSet,0.0);
			bool			applyNow = CheckEqual(doDeferSet,0.0);
			_String			*constraintAccumulator = (_String*)checkPointer(new _String(128L,true));
			
			if (applyNow)
			{
				deferSetFormula = new _SimpleList;
				checkPointer (deferSetFormula);
			}
			
			for (ind1 = 0; ind1 < theConstraints.lLength; ind1++)
			{
				replicateSource = (_String*)(theConstraints(ind1)->toStr());
				if (applyNow)
				{
					_Formula rhs, lhs;
					long	 varRef;
					ind2 = Parse (&rhs,*replicateSource,varRef,chain.nameSpacePrefix,&lhs);
					ExecuteFormula(&rhs,&lhs,ind2,varRef);
				}
				(*constraintAccumulator) << replicateSource;
				(*constraintAccumulator) << ';';
				(*constraintAccumulator) << '\n';
				//ReportWarning (*replicateSource);
				DeleteObject (replicateSource);
			}
			constraintAccumulator->Finalize();
			ReportWarning (*constraintAccumulator);
			CheckReceptacleAndStore (&lastSetOfConstraints,"ReplicateConstraint",false,new _FString(constraintAccumulator),false);
			if (applyNow)
				FinishDeferredSF();				
		}
	}
}

//____________________________________________________________________________________	

void	  _ElementaryCommand::ExecuteCase57 (_ExecutionList& chain)
{
	chain.currentCommand++;

	_String errStr;
			
	_Variable * storeResultIn = CheckReceptacle (&AppendContainerName(*(_String*)parameters(0),chain.nameSpacePrefix), blGetNeutralNull, true),
			  *	sv		  	  = FetchVar(LocateVarByName (AppendContainerName(*(_String*)parameters(2),chain.nameSpacePrefix))),
			  *	nsv	  	 	  = FetchVar(LocateVarByName (AppendContainerName(*(_String*)parameters(3),chain.nameSpacePrefix)));
			  
	_Parameter itCountV		  = ProcessNumericArgument ((_String*)parameters(4),chain.nameSpacePrefix);
	
	_String   * lfName		  = (_String*)parameters(1);
	
	long		f = FindLikeFuncName(AppendContainerName(*lfName,chain.nameSpacePrefix));
	
	if (f>=0)
	{
		if (sv && sv->ObjectClass () == MATRIX)
		{
			if (nsv && nsv->ObjectClass () == MATRIX)
			{
				_Matrix * sMatrix  = (_Matrix*)((_Matrix*)sv->Compute())->ComputeNumeric();
				_Matrix * nsMatrix = (_Matrix*)((_Matrix*)nsv->Compute())->ComputeNumeric();
				
				sMatrix->CheckIfSparseEnough (true);
				nsMatrix->CheckIfSparseEnough (true);
				
				if (   sMatrix->GetHDim()  == sMatrix->GetVDim() && 
					   nsMatrix->GetHDim() == nsMatrix->GetVDim() && 
					   sMatrix->GetHDim()  ==  nsMatrix->GetVDim() )
				{
					_LikelihoodFunction * theLF = (_LikelihoodFunction*)likeFuncList (f);
					
					if (((_DataSetFilter*)dataSetFilterList (theLF->GetTheFilters() (0)))->GetDimension (true) == sMatrix->GetHDim())
					{
						long itCount = itCountV;
						if (itCount>0)
						{
							_AssociativeList * res = theLF->SimulateCodonNeutral ((_Matrix*)sMatrix, (_Matrix*)nsMatrix, itCount);
						    storeResultIn->SetValue (res,false);
						}
						else
							errStr = "Invalid iterations per character state";
					}
					else
						errStr = "Incompatible data and cost matrices";
					
				}
				else
					errStr = "Incompatible syn and non-syn cost matrix dimensions";
			}		
			else
				errStr = "Invalid non-syn cost matrix argument";
		}
		else
			errStr = "Invalid syn cost matrix argument";
		
	}
	else
		errStr = _String("Likelihood function ") & *lfName & " has not been defined";
			   
	if (errStr.sLength)
	{
		errStr = errStr & " in call to " & blGetNeutralNull;
		WarnError (errStr);
	}
}

//____________________________________________________________________________________	

void	  _ElementaryCommand::ExecuteCase58 (_ExecutionList& chain)
{
	chain.currentCommand++;

	_String 	errStr;
	_String * 	profileCode   = (_String*)parameters(0);
			  
	if (*profileCode == _String ("START"))
	{
		if (chain.profileCounter)
			DeleteObject (chain.profileCounter);
		checkPointer(chain.profileCounter = new _Matrix (chain.lLength, 2, false, true));
		chain.doProfile = 1;
	}
	else
		if (*profileCode == _String ("PAUSE"))
			chain.doProfile = 2;
		else
			if (*profileCode == _String ("RESUME"))
				chain.doProfile = 1;
			else
			{
				_Variable * outVar = CheckReceptacle (&AppendContainerName(*profileCode,chain.nameSpacePrefix), blHBLProfile, true);
				if (outVar)
				{
					if (chain.profileCounter)
					{
						_AssociativeList * profileDump = new _AssociativeList;
						checkPointer 	 (profileDump);
						
						_SimpleList		 instructions;
						_List 			 descriptions;
						
						for (long k=1; k<2*chain.lLength; k+=2)
						{
							if (chain.profileCounter->theData[k] > 0.0)
							{
								instructions << k/2;
								_String * desc = (_String*)((_ElementaryCommand*)chain(k/2))->toStr();
								descriptions << desc;
								DeleteObject (desc);
							}
						}
						
						_Matrix 		* execProfile = new _Matrix (instructions.lLength,2,false,true),
										* instCounter = new _Matrix (instructions),
										* descList	  = new _Matrix (descriptions);

						checkPointer 	(execProfile);
						checkPointer 	(instCounter);
						checkPointer 	(descList);
										
						long k2 = 0;				
						for (long m=1; m<2*chain.lLength; m+=2)
						{
							if (chain.profileCounter->theData[m] > 0.0)
							{
								execProfile->theData[k2++] = chain.profileCounter->theData[m];
								execProfile->theData[k2++] = chain.profileCounter->theData[m-1];
							}
						}
						
						_FString  aKey;
						*aKey.theString = "INSTRUCTION INDEX";
						profileDump->MStore (&aKey, instCounter, false);
						*aKey.theString = "INSTRUCTION";
						profileDump->MStore (&aKey, descList, false);
						*aKey.theString = "STATS";
						profileDump->MStore (&aKey, execProfile, false);
						outVar->SetValue (profileDump,false);
						chain.doProfile = 0;
						DeleteObject (chain.profileCounter);
						chain.profileCounter = nil;
					}
					else
						errStr = "Profiler dump invoked before #profile START; "; 
				}
			}

}

//____________________________________________________________________________________	

void	  _ElementaryCommand::ExecuteCase59 (_ExecutionList& chain)
{
	chain.currentCommand++;

	for (long objCount = 0; objCount < parameters.lLength; objCount++)
	{
		long	  f;
		if 	((f = FindLikeFuncName(AppendContainerName(*(_String*)parameters(objCount),chain.nameSpacePrefix))) >= 0)
			KillLFRecord (f,true);
		
		
	}		   
}

//____________________________________________________________________________________	

void	  _ElementaryCommand::ExecuteCase60 (_ExecutionList& chain)
{
	chain.currentCommand++;
	_String	theVersion = ProcessLiteralArgument ((_String*)parameters (0),chain.nameSpacePrefix);
				
	if (__KERNEL__VERSION__.toNum() < theVersion.toNum())
		WarnError (_String ("Current batch file requires at least version :")& theVersion &" of HyPhy. Please download an updated version from http://www.hyphy.org and try again.");
}

//____________________________________________________________________________________	

void	  _ElementaryCommand::ExecuteCase65 (_ExecutionList& chain)
{
	chain.currentCommand++;
	
	_String assertion (*(_String*)parameters(0));

	_Formula rhs, lhs;
	long	 varRef;
	if (Parse (&rhs,assertion,varRef,chain.nameSpacePrefix,&lhs) == HY_FORMULA_EXPRESSION)
	{
		_PMathObj assertionResult = rhs.Compute();
		if (assertionResult && assertionResult->ObjectClass () == NUMBER)
		{
			if (CheckEqual(assertionResult->Value(),0.0))
			{
				_Parameter whatToDo;
				checkParameter (assertionBehavior, whatToDo, 0.0);
                
                _String errMsg;
                
                if (parameters.lLength == 1)
                {
                    errMsg = _String("Assertion '") & *(_String*)parameters(0) & "' failed.";
                }
				else
                {
                    errMsg = ProcessLiteralArgument((_String*)parameters(1),chain.nameSpacePrefix);
				}
                
				if (CheckEqual (whatToDo, 1.))
				{
					StringToConsole (errMsg);
					NLToConsole();
					chain.GoToLastInstruction ();
				}
				else
					WarnError (errMsg);	
			}
			return;
		}
	}
	WarnError ("Assertion statement could not be computed or was not numeric.");
}

//____________________________________________________________________________________	

void	  _ElementaryCommand::ExecuteCase61 (_ExecutionList& chain)
{
	chain.currentCommand++;
	
	_PMathObj			avl1 	= FetchObjectFromVariableByType (&AppendContainerName(*(_String*)parameters(1),chain.nameSpacePrefix),ASSOCIATIVE_LIST), 
						avl2 	= FetchObjectFromVariableByType (&AppendContainerName(*(_String*)parameters(2),chain.nameSpacePrefix),ASSOCIATIVE_LIST), 
						start	= parameters.lLength>3?FetchObjectFromVariableByType (&AppendContainerName(*(_String*)parameters(3),chain.nameSpacePrefix),NUMBER):nil;

	if (! (avl1 && avl2))
		WarnError (_String ("Both arguments (") & *(_String*)parameters(1) & " and " & *(_String*)parameters(2) & " in a call to SCFG = ... must be evaluate to associative arrays");
	else
	{
		Scfg    * scfg		= new Scfg ((_AssociativeList*)avl1,(_AssociativeList*)avl2,start?start->Value():0);
		_String scfgName    = AppendContainerName(*(_String*)parameters(0),chain.nameSpacePrefix);
		long 	f			= FindSCFGName (scfgName);

		if (f==-1)
		{
			for (f=0; f<scfgNamesList.lLength; f++)
				if (((_String*)scfgNamesList(f))->sLength==0)
					break;
					
			if (f==scfgNamesList.lLength)
			{
				scfgList << scfg;
				scfgNamesList&&(&scfgName);
				DeleteObject (scfg);
			}
			else
			{
				scfgNamesList.Replace(f,&scfgName,true);
				scfgList.lData[f] = (long)scfg;
			}
		}
		else
		{
			scfgNamesList.Replace(f,&scfgName,true);
			scfgList.Replace(f,scfg,false);
		}
	}	
}

//____________________________________________________________________________________	

void	  _ElementaryCommand::ExecuteCase63 (_ExecutionList& chain)
{
	chain.currentCommand++;
	
	/*_PMathObj			avl1 	= FetchObjectFromVariableByType ((_String*)parameters(1),ASSOCIATIVE_LIST), 
						avl2 	= FetchObjectFromVariableByType ((_String*)parameters(2),ASSOCIATIVE_LIST), 
						start	= parameters.lLength>3?FetchObjectFromVariableByType ((_String*)parameters(3),NUMBER):nil;

	if (! (avl1 && avl2))
	{
		_String errMsg = _String ("Both arguments (") & *(_String*)parameters(1) & " and " & *(_String*)parameters(2) & " in a call to SCFG = ... must be evaluate to associative arrays";
		WarnError (errMsg);
	}	
	else
	{
		Scfg    * scfg   = new Scfg ((_AssociativeList*)avl1,(_AssociativeList*)avl2,start?start->Value():0);
		_String * str    = (_String*)parameters(0);
		long 	f 		 = FindSCFGName (*str);

		if (f==-1)
		{
			for (f=0; f<scfgNamesList.lLength; f++)
				if (((_String*)scfgNamesList(f))->sLength==0)
					break;
					
			if (f==scfgNamesList.lLength)
			{
				scfgList << scfg;
				scfgNamesList&&(str);
				DeleteObject (scfg);
			}
			else
			{
				scfgNamesList.Replace(f,str,true);
				scfgList.lData[f] = (long)scfg;
			}
		}
		else
		{
			scfgNamesList.Replace(f,str,true);
			scfgList.Replace(f,scfg,false);
		}
	}	*/
}


//____________________________________________________________________________________	

void	_ElementaryCommand::ExecuteCase64 (_ExecutionList& chain)
{
	chain.currentCommand++;
#if defined __AFYP_REWRITE_BGM__
	_PMathObj	avl1	= FetchObjectFromVariableByType (&AppendContainerName(*(_String*)parameters(1),chain.nameSpacePrefix), ASSOCIATIVE_LIST);
	
	if (! (avl1))
	{
		WarnError (_String ("Argument (") & *(_String*)parameters(1) & " in call to BGM = ... must evaluate to associative array");
	}
	else
	{
		_BayesianGraphicalModel	* bgm	= new _BayesianGraphicalModel ((_AssociativeList *) avl1);
		
#else
	_PMathObj	avl1	= FetchObjectFromVariableByType (&AppendContainerName(*(_String*)parameters(1),chain.nameSpacePrefix), ASSOCIATIVE_LIST),
				avl2	= FetchObjectFromVariableByType (&AppendContainerName(*(_String*)parameters(2),chain.nameSpacePrefix), ASSOCIATIVE_LIST);
	
	if (! (avl1 && avl2))
	{
		WarnError (_String ("Both arguments (") & *(_String*)parameters(1) & " and " & *(_String*)parameters(2) & 
				   " in a call to BGM = ... must evaluate to associative arrays");
	}
	else
	{
		// is this a dynamic Bayesian network?
		_Parameter		dynamicArg;
		checkParameter (isDynamicGraph, dynamicArg, 0.);
		bool is_dynamic_graph = (dynamicArg > 0) ? TRUE : FALSE;
		
		Bgm		* bgm;	// pointer to base class
		
		long	num_nodes	= ((_Matrix *) avl1)->GetSize() + ((_Matrix *) avl2)->GetSize();
		
		if (num_nodes < 2)
		{
			WarnError (_String("Cannot construct a Bgm object on less than 2 nodes, received ") & num_nodes);
			return;
		}
		else
		{
			if (is_dynamic_graph)	bgm = new _DynamicBgm ((_AssociativeList*)avl1, (_AssociativeList*)avl2);
			else					bgm = new Bgm ((_AssociativeList*)avl1, (_AssociativeList*)avl2);
		}
#endif
		_String bgmName	    = AppendContainerName (*(_String *) parameters(0), chain.nameSpacePrefix);
		long	bgmIndex	= FindBgmName (bgmName);
		
		if (bgmIndex == -1)	// not found
		{
			for (bgmIndex = 0; bgmIndex < bgmNamesList.lLength; bgmIndex++)
			{
				// locate empty strings in name list
				if (((_String *)bgmNamesList(bgmIndex))->sLength == 0)
					break;
			}
			
			if (bgmIndex == bgmNamesList.lLength)
			{
				// reached end of list without finding empty string, append new string
				bgmList << bgm;
				bgmNamesList && (&bgmName);
				DeleteObject (bgm);
			}
			else
			{
				// replace empty string in list
				bgmNamesList.Replace (bgmIndex, &bgmName, true);
				bgmList.Replace (bgmIndex, bgm, false);
			}
		}
		else // 20070626: SLKP edit to deal with already existing BGMs
		{
			bgmNamesList.Replace(bgmIndex,&bgmName,true);
			bgmList.Replace(bgmIndex,bgm,false);
		}
	}
}


//____________________________________________________________________________________	
bool	_ElementaryCommand::ConstructSCFG (_String&source, _ExecutionList&target)
// syntax: SCFG ident = (Rules1, Rules2 <,start>)
{
	
	long	mark1 = source.FirstSpaceIndex(0,-1,1), 
			mark2 = source.Find ('=', mark1, -1);
	
	_String	scfgID (source, mark1+1,mark2-1);

	if (mark1==-1 || mark2==-1 || mark1+1>mark2-1 || !scfgID.IsValidIdentifier(true))
	{
		WarnError ("SCFG declaration missing a valid identifier");
		return false;
	}
			
	_List pieces;

	mark1 = source.Find ('(',mark2,-1);
	if (mark1 >= 0)
		ExtractConditions (source,mark1+1,pieces,',');
	
	if (pieces.lLength != 2 && pieces.lLength != 3)
	{
		WarnError ("Expected: SCFG ident = (Rules1, Rules2 <,start>)");
		return false;
	}
	
	_ElementaryCommand * scfg = new _ElementaryCommand (61);
	
	scfg->parameters	&&(&scfgID);
	scfg->addAndClean(target,&pieces,0);
	return true;
}

//____________________________________________________________________________________	
bool	_ElementaryCommand::ConstructNN (_String&source, _ExecutionList&target)
// syntax: NeuralNet ident = (InMatrix,OutMatrix,HiddenNodes)
{
	
	long	mark1 = source.FirstSpaceIndex(0,-1,1), 
			mark2 = source.Find ('=', mark1, -1);
	
	_String	nnID (source, mark1+1,mark2-1);

	if (mark1==-1 || mark2==-1 || mark1+1>mark2-1 || !nnID.IsValidIdentifier(true))
	{
		WarnError ("NeutalNet declaration missing a valid identifier");
		return false;
	}
		
	
	_List pieces;

	mark1 = source.Find ('(',mark2,-1);
	if (mark1 >= 0)
		ExtractConditions (source,mark1+1,pieces,',');
	
	if (pieces.lLength != 3)
	{
		WarnError ("NeuralNet ident = (InMatrix,OutMatrix,HiddenNodes)");
		return false;
	}
	
	_ElementaryCommand * nn = new _ElementaryCommand (63);
	nn->parameters	&& (&nnID);
	nn->addAndClean(target,&pieces,0);
	return true;
}
	
//____________________________________________________________________________________	
bool	_ElementaryCommand::ConstructBGM (_String&source, _ExecutionList&target)
// syntax: BGM ident = (<discrete nodes>, <continuous nodes>)
{
	// locate ident in HBL string
	long	mark1 = source.FirstSpaceIndex(0,-1,1), 
	mark2 = source.Find ('=', mark1, -1);
	
	// assign ident to _String variable
	_String	bgmID (source, mark1+1,mark2-1);
	
	if (mark1==-1 || mark2==-1 || mark1+1>mark2-1 || !bgmID.IsValidIdentifier(true))
	{
		WarnError ("BGM declaration missing a valid identifier");
		return false;
	}
	
	
	// extract arguments from remainder of HBL string
	_List pieces;
	
	mark1 = source.Find ('(',mark2,-1);
	if (mark1 >= 0)
		ExtractConditions (source,mark1+1,pieces,',');
	
#if defined __AFYP_REWRITE_BGM__
	if (pieces.lLength != 1)
	{
		WarnError ("Expected: BGM ident = (<nodes>)");
		return false;
	}
#else
	if (pieces.lLength < 2)
	{
		WarnError ("Expected: BGM ident = (<discrete nodes>, <continuous nodes>)");
		return false;
	}
#endif
	
	_ElementaryCommand * bgm = new _ElementaryCommand (64);
	bgm->parameters	&& (&bgmID);
	bgm->addAndClean(target,&pieces,0);
	return true;
}
	
	
//____________________________________________________________________________________	
bool	_ElementaryCommand::ConstructAssert (_String&source, _ExecutionList&target)
// syntax: assert (statement,message on failure)
{
	
	// extract arguments from remainder of HBL string
	_List pieces;
	
	ExtractConditions (source,blAssert.sLength,pieces,',');
	
	if (pieces.lLength != 2 && pieces.lLength != 1)
	{
		WarnError ("Expected: assert (statement,<message on failure>");
		return false;
	}

	_ElementaryCommand * bgm = new _ElementaryCommand (65);
	bgm->addAndClean(target,&pieces,0);
	return true;
}


//____________________________________________________________________________________	

_Parameter	 AlignStrings 	(_String* s1,_String* s2,_SimpleList& cmap,_Matrix* ccost,char gap,_Parameter gopen,_Parameter gextend,_Parameter gopen2,_Parameter gextend2,bool doLocal,bool doAffine,bool,_List& store)
{
	_String *res1 = new _String (s1->sLength+1, true),
			*res2 = new _String (s2->sLength+1, true);
			
	checkPointer (res1);
	checkPointer (res2);
	
	_Parameter	 score = 0.;
	
	if (s1->sLength)
	{
		if (s2->sLength)
		{
			_Matrix 			scoreMatrix (s1->sLength+1,s2->sLength+1, false, true),
							   *gapScore1 = nil,
							   *gapScore2 = nil;
			_SimpleList			editOps 	(MAX(s1->sLength,s2->sLength));
			long				colCount = 	s2->sLength+1;

			if (doAffine)
			{
				gapScore1 = new _Matrix (s1->sLength+1,s2->sLength+1, false, true),
				gapScore2 = new _Matrix (s1->sLength+1,s2->sLength+1, false, true);
				checkPointer (gapScore1);
				checkPointer (gapScore2);
			}
			
			if (!doLocal)
			// initialize gap costs in first column and first row
			// they are 0 for local alignments
			{
				if (doAffine)
				{
					_Parameter cost = -gopen;
					for (long k=1; k < colCount; k++, cost-=gextend)
					{
						scoreMatrix.theData[k] = cost;
						gapScore1->theData[k]  = cost;
						gapScore2->theData[k]  = cost;
					}
						
					cost = -gopen2;
					
					gapScore1->theData[0] = -gopen;
					gapScore2->theData[0] = -gopen2;
					
					for (long m=colCount; m < (s1->sLength+1)*colCount; m+=colCount, cost-=gextend2)
					{
						scoreMatrix.theData[m] = cost;
						gapScore1->theData [m] = cost;
						gapScore2->theData [m] = cost;
					}
				}				
				else
				{
					_Parameter cost = -gopen;
					for (long m=1; m < colCount; m++, cost-=gopen)
						scoreMatrix.theData[m] = cost;
						
					cost = -gopen2;
					for (long k=colCount; k < (s1->sLength+1)*colCount; k+=colCount, cost-=gopen2)
						scoreMatrix.theData[k] = cost;
				}
			}
			else
			{
				if (doAffine)
				{
					for (long m=1; m < colCount; m++)
						gapScore2->theData[m] = -gopen2;						
					for (long m=colCount; m < (s1->sLength+1)*colCount; m+=colCount)
						gapScore1->theData[m] = -gopen;						
				}	
			}
			
			if (doAffine)
			{
				//long upto1 = doLocal?s1->sLength-1:s1->sLength,
					// upto2 = doLocal?s2->sLength-1:s2->sLength;
					 
				long mapL = ccost->GetVDim();
				for (long r=1; r<=s1->sLength; r++)
				{
					long	  c1 = cmap.lData[s1->sData[r-1]];
					for (long c=1; c<=s2->sLength; c++)
					{
						long	   mIndex 	= r*colCount+c,
								   mIndex2	= mIndex-colCount;
						
						_Parameter gscore1  = MAX(scoreMatrix.theData[mIndex2]-gopen2,gapScore2->theData[mIndex2]-((r>1)?gextend2:gopen2)), 	// gap in 2nd 
								   gscore2  = MAX(scoreMatrix.theData[mIndex-1]-gopen,gapScore1->theData[mIndex-1]-((c>1)?gextend:gopen)),    // gap in 1st
								   gscore3  = scoreMatrix.theData[mIndex2-1];     
								   
								   
						if (c1>=0)
						{
							long	   c2 = cmap.lData[s2->sData[c-1]];
							
							if (c2>=0)
								gscore3 += ccost->theData[c1*mapL+c2];
						}
						
						scoreMatrix.theData[mIndex] = MAX(gscore1,MAX(gscore2,gscore3));
						gapScore2->theData[mIndex]  = gscore1;
						gapScore1->theData[mIndex]  = gscore2;
					}
				}
			}
			else
			// populate the cost matrix row by row
			{
				long upto1 = s1->sLength,
					 upto2 = s2->sLength;
					 
				for (long r=1; r<=upto1; r++)
				{
					long	  c1 = cmap.lData[s1->sData[r-1]];
					for (long c=1; c<=upto2; c++)
					{
						_Parameter score1 = scoreMatrix.theData[(r-1)*colCount+c] - gopen2, // gap in 2nd 
								   score2 = scoreMatrix.theData[r*colCount+c-1]   - gopen,  // gap in 1st
								   score3 = scoreMatrix.theData[(r-1)*colCount+c-1];     
								   
						if (c1>=0)
						{
							long	   c2 = cmap.lData[s2->sData[c-1]];
							
							if (c2>=0)
								score3 += (*ccost)(c1,c2);
						}
						
						scoreMatrix.theData[r*colCount+c] = MAX(score1,MAX(score2,score3));
					}
				}
			}
				
			long p1 = s1->sLength,
				 p2 = s2->sLength;

			if (doLocal)
			{		
				score = scoreMatrix.theData[(s1->sLength+1)*colCount-1];
				
				for (long mc2=colCount-1; mc2<s1->sLength*colCount; mc2+=colCount)
					if (scoreMatrix.theData[mc2]>score)
					{
						score = scoreMatrix.theData[mc2];
						p1 = mc2/colCount;
					}

				for (long mc=s1->sLength*colCount; mc<(s1->sLength+1)*colCount-1; mc++)
					if (scoreMatrix.theData[mc]>score)
					{
						score = scoreMatrix.theData[mc];
						p1 = s1->sLength;
						p2 = mc-s1->sLength*colCount;
					}
					
				for (long mp = p1;mp<s1->sLength; mp++)
					editOps << -1;
				for (long mp2 = p2;mp2<s2->sLength; mp2++)
					editOps << 1;
			}
			else
				score = scoreMatrix.theData[(s1->sLength+1)*colCount-1];
			
			// backtrack now
			
			if (doAffine)
			{
				while (p1 && p2)
				{
					_Parameter bscore1 = gapScore2->theData[p1*colCount+p2], 
							   bscore2 = gapScore1->theData[p1*colCount+p2],
							   bscore3 = scoreMatrix.theData[(p1-1)*colCount+p2-1];
							   
					MismatchScore (s1,s2,p1,p2,cmap,ccost,bscore3);
					if (bscore1>=bscore2 && bscore1>=bscore3)
					{
						long  ci = p1*colCount+p2;
						p1--;
						editOps << -1;
						while (p1 && (scoreMatrix.theData[ci-colCount]-gopen2 <= gapScore2->theData[ci-colCount]-gextend2))
						{
							p1--;
							editOps << -1;
							ci -= colCount;
						}
					}
					else
					{
						if (bscore2>=bscore1 && bscore2>=bscore3)
						{
							long  ci = p1*colCount+p2;
							p2--;
							editOps << 1;
							while (p2 && (scoreMatrix.theData[ci-1]-gopen <= gapScore1->theData[ci-1]-gextend))
							{
								p2--;
								editOps << 1;
								ci --;
							}
						}
						else
						{
							p1--;
							p2--;
							editOps << 0;
						}
					}					 
				}
			}
			else
			{
				while (p1 && p2)
				{
					_Parameter bscore1 = scoreMatrix.theData[(p1-1)*colCount+p2]-gopen2,
							   bscore2 = scoreMatrix.theData[p1*colCount+p2-1]-gopen,
							   bscore3 = scoreMatrix.theData[(p1-1)*colCount+p2-1];
							   							   
					MismatchScore (s1,s2,p1,p2,cmap,ccost,bscore3);
					BacktrackAlign (editOps, p1,p2,bscore1,bscore2,bscore3);
				}
			}
			
			while (p1>0)
			{
				p1--;
				editOps << -1;
			}

			while (p2>0)
			{
				p2--;
				editOps << 1;
			}
			
			for (long eo = editOps.lLength-1; eo>=0; eo--)
				switch (editOps.lData[eo])
				{
					case 0:
						(*res1) << s1->sData[p1++];
						(*res2) << s2->sData[p2++];
						break;
					case 1:
						(*res1) << gap;
						(*res2) << s2->sData[p2++];
						break;
					case -1:
						(*res2) << gap;
						(*res1) << s1->sData[p1++];
						break;
				}
				
			
			/*_String		alignDebug ("alignScoreMatrix");
			_Variable * ad = CheckReceptacle (&alignDebug, empty, false);
			ad->SetValue (&scoreMatrix, true);
			if (doAffine)
			{
				_String		alignDebug ("alignScoreMatrixG1");
				_Variable * ad = CheckReceptacle (&alignDebug, empty, false);
				ad->SetValue (gapScore1, true);
				alignDebug  = ("alignScoreMatrixG2");
				ad = CheckReceptacle (&alignDebug, empty, false);
				ad->SetValue (gapScore2, true);
			}*/
			
			DeleteObject (gapScore1);
			DeleteObject (gapScore2);
			
		}
		else
		{
			(*res1) << *s1;
			for (long s1i = 0; s1i < s1->sLength; s1i++)
				(*res2) << gap;
				
			if (!doLocal)
				if (doAffine)
					return -gopen2-(s1->sLength-1)*gextend2;
				else
					return s1->sLength*gopen2;
		}
	}
	else
		if (s2->sLength)
		{
			(*res2) << *s2;
			for (long s2i = 0; s2i < s2->sLength; s2i++)
				(*res1) << gap;
				
			if (!doLocal)
				if (doAffine)
					return -gopen-(s2->sLength-1)*gextend;
				else
					return s2->sLength*gopen;
		}
		
	res1->Finalize();
	res2->Finalize();
	
	
	/* verify the score */
			
	/*long gap1c = 0,
		 gap2c = 0;
		 
	_Parameter scoreCheck = 0.;
		 
	for (long sp = 0; sp<res1->sLength; sp++)
	{
		char cs1 = res1->sData[sp],
			 cs2 = res2->sData[sp];
			 
		if (cs1 == gap)
		{
			if (gap1c)
				scoreCheck -= gextend;
			else
				scoreCheck -= gopen;
			gap2c = 0;
			gap1c++;
		}
		else
			if (cs2 == gap)
			{
				if (gap2c)
					scoreCheck -= gextend2;
				else
					scoreCheck -= gopen2;
				gap1c = 0;
				gap2c++;
			}
			else
			{
				gap1c = 0;
				gap2c = 0;
				long code1 = cmap.lData[cs1],
					 code2 = cmap.lData[cs2];
				
				if (code1 >=0 && code2 >=0 )
					scoreCheck += (*ccost)(code1,code2);
			}
			
	}
	
	char checkScore [256];
	sprintf (checkScore, "\nScore check: %g\n", scoreCheck);
	BufferToConsole (checkScore);*/

	store.AppendNewInstance(res1);
	store.AppendNewInstance(res2);
	
	return score;
}

//____________________________________________________________________________________	

_Parameter		LinearSpaceAlign (_String *s1,					// first string
							  _String *s2,						// second string
							  _SimpleList& cmap,				// char -> position in scoring matrix mapper
							  _Matrix*    ccost,				// NxN matrix of edit distances on characters
							  _Parameter gopen,					// the cost of opening a gap in sequence 1
							  _Parameter gextend,				// the cost of extending a gap in sequence 1 (ignored unless doAffine == true)
							  _Parameter gopen2,				// the cost of opening a gap in sequence 2
							  _Parameter gextend2,				// the cost of opening a gap in sequence 2	 (ignored unless doAffine == true)
							  bool doLocal,						// ignore prefix and suffix gaps
							  bool doAffine,					// use affine gap penalties
							  _SimpleList& ops,					// edit operations for the optimal alignment
							  _Parameter   scoreCheck,			// check the score of the alignment
							  long		   from1,
							  long		   to1,
							  long		   from2,
							  long		   to2,
							  _Matrix	   **buffer,				// matrix storage,
							  char		   parentGapLink,
							  char		   *ha
							  )
{
	if (to2 == from2 || to1 == from1)
		return 0;
			 
	long					midpoint = (from1 + to1)/2,
							span	 = to2-from2,
							span1	 = to1-from1;
	
	if						(span1 > 1)
	{
		CostOnly				(s1,s2,from1,from2,midpoint,to2,false,false,cmap,ccost,gopen,gextend,gopen2,gextend2,doLocal,doAffine,*(buffer[0]), buffer[1], buffer[2], parentGapLink>=2, ha);
		CostOnly				(s1,s2,midpoint,from2,to1,to2,true,true,  cmap,ccost,gopen,gextend,gopen2,gextend2,doLocal,doAffine,*(buffer[3]), buffer[4], buffer[5],   2*(parentGapLink%2), ha+s2->sLength+1);
	}
	else
		CostOnly				(s1,s2,from1,from2,to1,to2,false,false,cmap,ccost,gopen,gextend,gopen2,gextend2,doLocal,doAffine,*(buffer[0]), buffer[1], buffer[2], (parentGapLink>=2), ha);
	
	_Parameter maxScore = -1e100;	
	long	   maxIndex = 0;
	bool	   gapLink	= false;
	char	   alignmentKind    = 0;

	_Parameter	  gapOffsetScore   = gopen2-gextend2;	
	if (!doAffine)
	{
		if (span1 > 1)
		{
			for (long k = 0; k <= span; k++)
			{
				_Parameter currentScore = buffer[0]->theData[k] + buffer[3]->theData[span-k];
				if (currentScore > maxScore)
				{
					maxScore = currentScore; 
					maxIndex = k;
				}
			}
		}
		else // handle the case of a single row span correctly
		{
			for (long k = 0; k <= span; k++)
			{
				_Parameter currentScore		= buffer[0]->theData[k];
				
				if (!doLocal || to1 != s1->sLength)
					currentScore -= gopen*(span-k);

				if (currentScore > maxScore)
				{
					maxScore		= currentScore; 
					alignmentKind	= ha[k];
					maxIndex = k;
				}
			}
		}
	}
	else
	{
		if (span1 > 1)
		{
			// two cases here: no-gap link 
			// or gap-to-gap link
			
			for (long k = 0; k <= span; k++)
			{
				_Parameter currentScoreNoGap	= buffer[0]->theData[k] + buffer[3]->theData[span-k],
						   currentScoreWithGap2	= buffer[2]->theData[k] + buffer[5]->theData[span-k] + gapOffsetScore;
				
		
				if (doAffine && ((from1 == 0 || from2==0) && k == 0 || (to1 == s1->sLength || to2 == s2->sLength) && k == span))
					currentScoreWithGap2 -= gapOffsetScore;
					
				if (currentScoreNoGap > maxScore)
				{
					maxScore = currentScoreNoGap; maxIndex = k;
					gapLink	 = false;
				}
				if (currentScoreWithGap2 > maxScore)
				{
					maxScore = currentScoreWithGap2; maxIndex = k;
					gapLink  = true;
				}
				/*printf ("[%d %d %d %d] (%d) %d %g %g: %g %g / %g %d\n", from1, to1, from2, to2, parentGapLink,  k, 
						buffer[0]->theData[k],  buffer[3]->theData[span-k], buffer[2]->theData[k],  buffer[5]->theData[span-k],
						maxScore, maxIndex);*/
				
			}

		}
		else // handle the case of a single row span correctly
		{
			if (parentGapLink == 1)
			{
				maxIndex	  = span;
				maxScore	  = buffer[2]->theData[span];
				alignmentKind = 1;
			}
			else
			{
				for (long k = 0; k <= span; k++)
				{
					_Parameter currentScoreNoGap	= buffer[0]->theData[k],
							currentScoreWithGap2	= buffer[2]->theData[k];
					
					if (!doLocal || to1 != s1->sLength) // indel in sequence 1
						if (span-k)
						{
							currentScoreNoGap		-= gopen;
							currentScoreWithGap2	-= gopen;
							if (span-k>1)
							{
								currentScoreNoGap    -= gextend*(span-k-1);
								currentScoreWithGap2 -= gextend*(span-k-1);
							}
						}
								
					/*printf ("[%d %d %d %d] %d %g %g: %g %g / %g %d\n", from1, to1, from2, to2, k, 
							buffer[0]->theData[k],  buffer[2]->theData[k], currentScoreNoGap, currentScoreWithGap2,
							maxScore, maxIndex);*/

					if (currentScoreNoGap > maxScore)
					{
						maxScore = currentScoreNoGap; maxIndex = k;
						alignmentKind	= ha[k];
					}
					if (currentScoreWithGap2 > maxScore)
					{
						maxScore = currentScoreWithGap2; maxIndex = k;
						alignmentKind	= 0;
					}
				}
			}
		}
	}
		
	if (span1 == 1)
	{
		if (alignmentKind == 2)
			ops.lData[from1+1] = from2+maxIndex-1;
		else
			if (alignmentKind == 0 && maxIndex == 0/*&& to2 == s2->sLength && to1 == s1->sLength*/)
				ops.lData[from1+1]  = -3;
	}
	else
	{
	
		_Parameter check1 = buffer[0]->theData[maxIndex],
				check2 = buffer[3]->theData[span-maxIndex];
	
		if (span1>1)
		{
			if (maxIndex > 0)
			{
				char gapCode = gapLink;
				if (parentGapLink >= 2)
					gapCode += 2;
				LinearSpaceAlign (s1,s2,cmap,ccost,gopen,gextend,gopen2,gextend2,doLocal,doAffine,ops,check1, from1, midpoint, from2, from2 + maxIndex, buffer, gapCode, ha); 
			}
			else
				if (from2 == 0)
					for (long k = from1; k < midpoint; k++) 
						ops.lData[k+1] = -3;
			
			if (maxIndex < span)
			{
				char gapCode = 2*gapLink;
				if (parentGapLink % 2 == 1)
					gapCode ++;
				LinearSpaceAlign (s1,s2,cmap,ccost,gopen,gextend,gopen2,gextend2,doLocal,doAffine,ops,check2, midpoint, to1, from2 + maxIndex, to2, buffer, gapCode, ha); 
			}
		}
	}
	return maxScore;
}
							  
//____________________________________________________________________________________	

#define _ALIGNMENT_NOLOCAL		0x00
#define _ALIGNMENT_LOCAL_START	0x01
#define _ALIGNMENT_LOCAL_END	0x02

//____________________________________________________________________________________	

_Parameter	 CostOnly 	(_String* s1,				// first string
						 _String* s2,				// second string
						 long     from1,			// start here in string1
						 long	  from2,			// start here in string2
						 long	  to1,				// up to here in string1 // not inclusive
						 long	  to2,				// up to here in string2 // not inclusive
						 bool	  rev1,				// reverse string1
						 bool	  rev2,				// reverse string2
						 _SimpleList& cmap,			// char -> position in scoring matrix mapper
						 _Matrix* ccost,			// NxN matrix of edit distances on characters
						 _Parameter gopen,			// the cost of opening a gap in sequence 1
						 _Parameter gextend,		// the cost of extending a gap in sequence 1 (ignored unless doAffine == true)
						 _Parameter gopen2,			// the cost of opening a gap in sequence 2
						 _Parameter gextend2,		// the cost of opening a gap in sequence 2	 (ignored unless doAffine == true)
						 bool doLocal,				// ignore prefix and suffix gaps
						 bool doAffine,				// use affine gap penalties
						 _Matrix& scoreMatrix,		// where to write the last row of the scoring matrix
						 _Matrix * gapScore1,		// where to write the last row of open gap in 1st sequence matrix (ignored unless doAffine == true)
						 _Matrix * gapScore2,		// same but for open gap in 2nd sequence matrix
						 char	   secondGap,
						 char	* howAchieved)
{	
	_Parameter	 score    = 0.;
	
	long		 s1Length = to1-from1,
				 s2Length = to2-from2;
	
	
	bool		 doLocal1S = false,
				 doLocal1E = false,
				 doLocal2S = false,
				 doLocal2E = false;
	
	if (doLocal)
	{
		if (rev1)
		{
			doLocal1S = (to1==s1->sLength);
			doLocal1E = from1 == 0;
			//doLocal1 = (to1==s1->sLength)*_ALIGNMENT_LOCAL_START + (from1 == 0)*_ALIGNMENT_LOCAL_END; 	
		}
		else
		{
			doLocal1E = (to1==s1->sLength);
			doLocal1S = from1 == 0;
			//doLocal1 = (from1==0)*_ALIGNMENT_LOCAL_START + (to1==s1->sLength)*_ALIGNMENT_LOCAL_END; 
		}
		if (rev2)
		{
			doLocal2E = from2 == 0;
			doLocal2S = (to2==s2->sLength);
			//doLocal2 = (to2==s2->sLength)*_ALIGNMENT_LOCAL_START + (from2 == 0)*_ALIGNMENT_LOCAL_END; 			
		}
		else
		{
			doLocal2S = from2 == 0;
			doLocal2E = (to2==s2->sLength);
			//doLocal2 = (from2==0)*_ALIGNMENT_LOCAL_START + (to2==s2->sLength)*_ALIGNMENT_LOCAL_END; 
		}
	}
	
	if (s1Length) 
	// first string not empty
	{
		if (s2Length) 
		// second string not empty
		{
			_Parameter			aux2;	
			long				colCount = s2Length+1;
			
			scoreMatrix.theData[0] = 0.;
			if (doAffine)
				gapScore1->theData[0] = gapScore2->theData[0] = 0.;
			
			
			if (doLocal1S == 0)
			{
				_Parameter cost = -gopen;
				if (doAffine)
				{
					for (long k=1; k < colCount; k++, cost-=gextend)
					{
						scoreMatrix.theData[k]  = cost;
						gapScore1->theData [k]  = cost;
						gapScore2->theData [k]  = cost;
					}
				}
				else
					for (long m=1; m < colCount; m++, cost-=gopen)
						scoreMatrix.theData[m] = cost;
			}
			else
			{
				for (long k=1; k < colCount; k++)
					scoreMatrix.theData[k] = 0.;
				
				if (doAffine)
				{
					for (long k=1; k < colCount; k++)
					{
						gapScore1->theData[k] = 0;
						gapScore2->theData[k] = -(secondGap==1?gextend2:gopen2); 
						// prefix gaps in the second sequence
					}
					gapScore1->theData[0] = -gopen;
				}				
			}
			
						
			long mapL = ccost->GetVDim(); // how many valid characters
			
			if (doAffine)
			{
				aux2 = 0.;
				
				if (doLocal2S == 0)
					gapScore1->theData[0] = gapScore2->theData[0] = -(secondGap==1?gextend2:gopen2); 
				
				from2 --;
				from1 --;
				for (long r=1; r<=s1Length; r++) // iterate by rows
				{
					long	  c1 = cmap.lData[s1->sData[rev1?(to1-r):(from1+r)]];
					
					if (doLocal2S)
						aux2		= 0.;
					else
					{
						if (r>1)
							aux2	     = -((r-2)*gextend2 + (secondGap==1?gextend2:gopen2)); 
						scoreMatrix.theData[0] = -((secondGap==1?gextend2:gopen2) + (r-1)*gextend2);
					}
					
					for (long c=1; c<=s2Length; c++) // iterate by columns
					{
						_Parameter gscore1  ,			// gap in 2nd 
								   gscore2  ,			// gap in 1st
								   gscore3  = aux2,		// no gap
								   t;
						
						// if secondGap == 2, then we MUST _start_ with a gap in the 2nd sequence
						
						
						if (doLocal1E && r == s1Length)
						{
							//gscore2 = MAX(scoreMatrix.theData[c-1],gapScore1->theData[c-1]);
							gscore2 = scoreMatrix.theData[c-1];
							if (gapScore1->theData[c-1] > gscore2)
								gscore2 = gapScore1->theData[c-1];
						}
						else
						{
							gscore2 = scoreMatrix.theData[c-1]-gopen;
							t       = gapScore1->theData[c-1]-((c>1)?gextend:gopen);
							if (t > gscore2)
								gscore2 = t;
						}
						
						if (doLocal2E && c == s2Length)
						{
							//gscore1 = MAX(scoreMatrix.theData[c],gapScore2->theData[c]);	
							gscore1 = scoreMatrix.theData[c];
							if (gscore1 < gapScore2->theData[c])
								gscore1 = gapScore2->theData[c];
						}
						else
						{
							//gscore1 = MAX(scoreMatrix.theData[c]-gopen2,gapScore2->theData[c]-((r>1)?gextend2:gopen2));
							gscore1 = scoreMatrix.theData[c]-gopen2;
							t		= gapScore2->theData[c]-((r>1)?gextend2:gopen2);
							if (t > gscore1)
								gscore1 = t;
						}
						// either open a new gap from a character; or continue an existing one
						// if this is the second row, then we start a gap in the second sequence -|
						
						if (c1>=0)
						{
							long	   c2 = cmap.lData[s2->sData[rev2?(to2-c):(from2+c)]];
							
							if (c2>=0)
								gscore3 += ccost->theData[c1*mapL+c2];
						}
						
						aux2					= scoreMatrix.theData[c];
						char					  option = 0;
						t						= gscore2;
						
						
						if (r > 1 || secondGap == 0)
						{
							if (gscore1 > gscore2)
							{
								t = gscore1;
								option				   = 1;
							}
							if (gscore3 > t)
							{
								t					   = gscore3;
								option				   = 2;
							}
						}
						scoreMatrix.theData[c] = t;
						if (howAchieved)
							howAchieved[c] = option;
						
						//if (rev2 && secondGap==2 && c == s2Length)
						//	gscore1 = MAX(scoreMatrix.theData[c]-gextend2,gapScore2->theData[c]-((r>1)?gextend2:gopen2));							

						gapScore2->theData [c]  = gscore1;
						gapScore1->theData [c]  = gscore2;
						
					}
					
					if (doLocal2S && r < s1Length)
					{
						gapScore1->theData[0]-=gextend2;gapScore2->theData[0]-=gextend2;
					}
				}
			}
			else
				// populate the cost matrix row by row
			{
				aux2 = 0.;
				for (long r=1; r<=s1Length; r++)
				{
					if (doLocal2S)
						aux2		= 0.;
					else
					{
						scoreMatrix.theData[0] = -(gopen2 * r);
						if (r>1)
							aux2	     = -((r-1)*gopen2); 
					}
					
					//printf ("%d: %g\t", r, scoreMatrix.theData[0]);
					long	  c1 = cmap.lData[s1->sData[rev1?(to1-r):(from1+r-1)]];
					
					for (long c=1; c<=s2Length; c++)
					{
						_Parameter score1 = scoreMatrix.theData[c], // gap in 2nd 
								   score2 = scoreMatrix.theData[c-1],  // gap in 1st
								   score3 = aux2;     
						
						if (c < s2Length || doLocal2E == 0)
							score1 -= gopen2;
						if (r < s1Length || doLocal1E == 0)
							score2 -= gopen;
						
						if (c1>=0)
						{
							long	   c2 = cmap.lData[s2->sData[rev2?(to2-c):(from2+c-1)]];
							
							if (c2>=0)
								score3 += ccost->theData[c1*mapL+c2];
						}
						
						aux2					= scoreMatrix.theData[c];
						char					option = 0;
						scoreMatrix.theData[c]  = score1;
						if (score2 > score1)
						{
							scoreMatrix.theData[c] = score2;
							option				   = 1;
						}
						if (score3 > scoreMatrix.theData[c])
						{
							scoreMatrix.theData[c] = score3;
							option				   = 2;
						}
						if (howAchieved)
							howAchieved[c] = option;
					}
					//printf ("\n");
				}
				
			}
			score = scoreMatrix.theData[s2Length];
		}
		else // 2nd string empty
		{
			if ((doLocal2S || doLocal2E) == false)
			{
				if (doAffine)
					score = gopen2+gextend2*s1Length;
				else
					score = gopen2 * s1Length;
			}
		}
	}
	else // first string empty
		if (s2Length) // second string not empty
		{
			if ((doLocal1S || doLocal1E) == false)
			{
				score = -gopen;
				
				scoreMatrix.theData[0] = 0.0;
				if (doAffine)
				{	
					gapScore1->theData[0] = gapScore2->theData[0] = 0.0;
					for (long k = 1; k <= s2Length; k++, score-=gextend)
						scoreMatrix.theData[k] = gapScore1->theData[k] = gapScore2->theData[k] = score;
					
					score += gextend;
				}
				else
				{
					for (long k = 1; k <= s2Length; k++, score-=gopen)
						scoreMatrix.theData[k] = score;
					score += gopen;
				}
			}
			else
			{
				for (long k = 0; k <= s2Length; k++)
					scoreMatrix.theData[k] = 0.;
				if (doAffine)
					for (long k = 0; k <= s2Length; k++)
					{
						gapScore1->theData[k] = 0.;
						gapScore2->theData[k] = 0.;
					}
			}

		}
	
	return score;
}


//____________________________________________________________________________________	

inline	void BacktrackAlign			(_SimpleList& editOps , long& p1, long& p2, _Parameter score1, _Parameter score2, _Parameter score3)
{
	if ((score1>=score2)&&(score1>=score3))
	{
		p1--;
		editOps << -1;
	}
	else
	{
		if ((score2>=score1)&&(score2>=score3))
		{
			p2--;
			editOps << 1;
		}
		else
		{
			p1--;
			p2--;
			editOps << 0;
		}
	}
}

//____________________________________________________________________________________	

inline	void MismatchScore			(_String* s1, _String*s2 , long p1, long p2, _SimpleList& cmap, _Matrix* ccost, _Parameter& score)
{
	long	  c1 = cmap.lData[s1->sData[p1-1]];
	if (c1>=0)
	{
		long	   c2 = cmap.lData[s2->sData[p2-1]];
		
		if (c2>=0)
			score += (*ccost)(c1,c2);
	}
}

//____________________________________________________________________________________	

void	RetrieveModelComponents (long mid, _Matrix*& mm, _Matrix*& fv, bool & mbf)
{
	if (mid >=0 && mid < modelTypeList.lLength)
	{
		if (modelTypeList.lData[mid] == 0)
			mm = (_Matrix*)FetchObjectFromVariableByTypeIndex(modelMatrixIndices.lData[mid],MATRIX);
		else
			mm = nil;
		
		long fvi = modelFrequenciesIndices.lData[mid];
		fv = (_Matrix*)FetchObjectFromVariableByTypeIndex(fvi>=0?fvi:(-fvi-1),MATRIX);
		mbf = (fvi>=0);
	}
	else
	{
		mm = fv = nil;
		mbf = false;
	}
}
	
//____________________________________________________________________________________	

void	RetrieveModelComponents (long mid, _Variable*& mm, _Variable*& fv, bool & mbf)
{
	if (modelTypeList.lData[mid] == 0)
		mm = LocateVar(modelMatrixIndices.lData[mid]);
	else
		mm = nil;
	
	long fvi = modelFrequenciesIndices.lData[mid];
	fv = LocateVar (fvi>=0?fvi:(-fvi-1));
	mbf = (fvi>=0);
}

//____________________________________________________________________________________	

bool	IsModelReversible (long mid)
{
	_Matrix *m = nil,
			*f = nil;
	bool	mbf;
	RetrieveModelComponents (mid, m, f, mbf);
	if (m&&f)
		return m->IsReversible(mbf?nil:f);
	return false;
}

	
//____________________________________________________________________________________	

void	ScanModelForVariables		 (long modelID, _AVLList& theReceptacle, bool inclG, long modelID2, bool inclCat)
{
	if (modelID != HY_NO_MODEL)
	{
		if (modelTypeList.lData[modelID] == 0)
			// standard rate matrix
			((_Matrix*) (LocateVar(modelMatrixIndices.lData[modelID])->GetValue()))->ScanForVariables2(theReceptacle,inclG,modelID2,inclCat);
		else
			// formula based
			((_Formula*)modelMatrixIndices.lData[modelID])->ScanFForVariables(theReceptacle, inclG, false, inclCat);
	}
}
	
 