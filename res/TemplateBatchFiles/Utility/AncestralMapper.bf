RequireVersion ("2.0020100914");

_ancestralRecoveryCache       = {};

/*******************************************
	call _buildAncestralCache function with the likelihood function ID 
	and a 0-based partition index to produce an
	internal structure storing internal states at 
	nodes;
	
	returns an integer index to reference
	the opaque structure for subsequent operations

*******************************************/

LoadFunctionLibrary			("TreeTools.ibf");

_bacCacheInstanceCounter = 0;

/*******************************************
	wrapper functions
*******************************************/

function _buildAncestralCache (_lfID, _lfComponentID)
{
	return _buildAncestralCacheInternal (_lfID, _lfComponentID, 0);
}

function _buildAncestralCacheSample (_lfID, _lfComponentID)
{
	return _buildAncestralCacheInternal (_lfID, _lfComponentID, 1);
}

/*******************************************
	internal function to do the work
*******************************************/

function _buildAncestralCacheInternal (_lfID, _lfComponentID, doSample)
{

/* 1; grab the information AVL from the likelihood function */
	ExecuteCommands ("GetString(_bac_lfInfo,"+_lfID+",-1);");
	if (Columns (_bac_lfInfo["Trees"]) <= _lfComponentID)
	{
		return -1;
	}
	
	_bac_treeID     = (_bac_lfInfo["Trees"])		[_lfComponentID];
	_bac_filterID 	= (_bac_lfInfo["Datafilters"])	[_lfComponentID];

/* 2; construct a temporary likelihood function with 
   the tree and the filter; and an AVL representation of the tree */
   
    if (doSample)
    {
		_ancestorCall = "SampleAncestors";
    }  
	else
	{
		_ancestorCall = "ReconstructAncestors";	
	}
	   
	ExecuteCommands ("_bac_tree_avl = " + _bac_treeID + "^0;"+
					 "DataSet _bac_ancDS = " + _ancestorCall + " ("+_lfID+",{{"+_lfComponentID+"}});"+
					 "GetString (_bacSequenceNames,"+ _bac_filterID +",-1);");
					 
/* 3; obtain ID->string mapping for the datafilter;
	  also deduce how many chars/state there are */

	ExecuteCommands 	  			("GetDataInfo (_bacCharHandles, `_bac_filterID`, \"CHARACTERS\")");
	ExecuteCommands 	  			("GetDataInfo (_bacCharProperties, `_bac_filterID`, \"PARAMETERS\")");
	ExecuteCommands 	  			("GetDataInfo (_bacFilterPatternMap, `_bac_filterID`)");
	_bacFilterDimension 		=   Columns (_bacCharHandles);
	_bacCharsPerState 			=   _bacCharProperties["ATOM_SIZE"];
	DataSetFilter	_bacAF		=   CreateFilter (_bac_ancDS,_bacCharsPerState,"","",_bacCharProperties["EXCLUSIONS"]);
	GetString						(_bacAncestralNames,_bacAF,-1);		
	GetDataInfo						(_bacAncestralPatternMap,_bacAF);
	
/* now start building a matrix of mapped states;
   
   each row in the matrix corresponds to a tree NODE in the same order as 
   they are listed in the tree AVL 
   
   each column correponds to a site
   
   the entries are as follows:
	   0.._bacFilterDimension-1 an unambigous character
	   -1 a deletion/complete ambiguity
	   less that -1 : -(entry)-2 maps to an AVL storing the state
					  and the complete character resolution (as an array)
 */
	 
	 _bacBranchCount		 	  = Abs(_bac_tree_avl)-1;
	 _bacMatrixOfResolutions 	  = {_bacBranchCount,_bacAF.sites};
	 _bacHandledResolutions  	  = {};
	 _bacHandledResolutionsChars  = {};
	 _bacHandledResolutionsAmbig  = {};
	 
 /* map where sequences are in the filter vs where they are in the 
	tree structure */
	    
	 _bacMapTreeNodeToDF = {};
	 _bacTreeAVLOrder 	 = {};
	 
	 for (_bacCounter = 1; _bacCounter < _bacBranchCount+1; _bacCounter = _bacCounter + 1)
	 {
	 	_bacTreeAVLOrder[(_bac_tree_avl[_bacCounter])["Name"]&&1] = _bacCounter;
	 }
	 
	 _bacFilterSequenceCount = Columns(_bacSequenceNames);
	 for (_bacCounter = 0; _bacCounter < _bacFilterSequenceCount; _bacCounter = _bacCounter + 1)
	 {
	 	_bacMapTreeNodeToDF[_bacTreeAVLOrder[_bacSequenceNames[_bacCounter]&&1]-1] = _bacCounter;
	 }
	 for (_bacCounter = 0; _bacCounter < Columns(_bacAncestralNames); _bacCounter = _bacCounter + 1)
	 {
	 	_bacMapTreeNodeToDF[_bacTreeAVLOrder[_bacAncestralNames[_bacCounter]&&1]-1] = _bacCounter+_bacFilterSequenceCount;
	 }
	 
 /* make some auxiliary variables */
 
 	_bacUnitRow     	= {1,Columns(_bacCharHandles)}["1"];
 	_bacSequenceRow     = {1,Columns(_bacCharHandles)}["_MATRIX_ELEMENT_COLUMN_"];

 /* loop over branches (rows) */
	 for (_bacBranchCounter = 1; _bacBranchCounter <= _bacBranchCount; _bacBranchCounter = _bacBranchCounter + 1)
	 {
	 	_bacRowIndex = _bacMapTreeNodeToDF[_bacBranchCounter-1];
		if (_bacRowIndex < _bacFilterSequenceCount) 
		/* an extant sequence */
		{
			ExecuteCommands ("GetDataInfo (_bacSequenceString,`_bac_filterID`,_bacRowIndex)");
		}
		else 
		/* an ancestral sequence */
		{
			GetDataInfo (_bacSequenceString,_bacAF,_bacRowIndex-_bacFilterSequenceCount);
		}

 /* loop over sites (columns) */
	 	for (_bacSiteCounter = 0; _bacSiteCounter < _bacAF.sites; _bacSiteCounter = _bacSiteCounter + 1)
	 	{
	 		_bacCurrentState = _bacSequenceString[_bacSiteCounter*_bacCharsPerState][(1+_bacSiteCounter)*_bacCharsPerState-1];
	 		_bacCurrentIndex = _bacHandledResolutions[_bacCurrentState];
	 		if (_bacCurrentIndex > 0)
	 		{
	 			_bacMatrixOfResolutions[_bacBranchCounter-1][_bacSiteCounter] = _bacCurrentIndex-1;
	 		}
	 		else
	 		{
	 			if (_bacCurrentIndex < 0)
	 			{
	 				_bacMatrixOfResolutions[_bacBranchCounter-1][_bacSiteCounter] = _bacCurrentIndex;
	 			}
	 			else
	 			/* handle the resolution */
	 			{
					if (_bacRowIndex < _bacFilterSequenceCount) 
					{
						ExecuteCommands ("GetDataInfo (_bacCharState,`_bac_filterID`,_bacRowIndex,_bacFilterPatternMap[_bacSiteCounter])");
					}
					else 
					{
						GetDataInfo (_bacCharState,_bacAF,_bacRowIndex-_bacFilterSequenceCount,_bacAncestralPatternMap[_bacSiteCounter]);
					}
					_bacResolutionCount = (_bacUnitRow*_bacCharState)[0];
	 				if (_bacResolutionCount == 1) /* fully resolved */
	 				{
	 					_bacHandledResolutions[_bacCurrentState] 					  = (_bacSequenceRow*_bacCharState)[0]+1;
		 				_bacMatrixOfResolutions[_bacBranchCounter-1][_bacSiteCounter] = _bacHandledResolutions[_bacCurrentState]-1;
	 				}
	 				else
	 				{
	 					if (_bacResolutionCount == Columns(_bacCharHandles)) /* gap/full ambig */
	 					{
	 						_bacHandledResolutions[_bacCurrentState] 							= -1;
	 						_bacMatrixOfResolutions[_bacBranchCounter-1][_bacSiteCounter]		= -1;
	 					}
	 					else
	 					{
	 						_bacHandledResolutions[_bacCurrentState] 							= -2-Abs(_bacHandledResolutionsAmbig);
	 						_bacHandledResolutionsAmbig[Abs(_bacHandledResolutionsAmbig)] 		= _bacCharState;
	 						_bacMatrixOfResolutions[_bacBranchCounter-1][_bacSiteCounter]		= _bacHandledResolutions[_bacCurrentState];
	 					}
	 				}
	 			}
	 		}
		}
	 }
	 
	 _bacAncestralCache = {};
	 _bacAncestralCache ["CHARS"]     	  = _bacCharHandles;
	 _bacAncestralCache ["MATRIX"]    	  = _bacMatrixOfResolutions;
	 _bacAncestralCache ["TREE_AVL"]  	  = _bac_tree_avl;
	 _bacAncestralCache ["AMBIGS"]    	  = _bacHandledResolutionsAmbig;
	 _bacAncestralCache ["AMBIGS_REV"]    = _bacReverseAmbigMap;
	
	 /* reverse map integer codes into ambig maps */

	 _bacCurrentState 		= Rows (_bacHandledResolutionsAmbig);
	 _bacResolutionCount	= Abs  (_bacHandledResolutionsAmbig);
	 _bacReverseAmbigMap	= {};
	 
	 for (_bacCounter = 0; _bacCounter < _bacResolutionCount; _bacCounter = _bacCounter + 1)
	 {
	 	_bacReverseAmbigMap[(_bacAncestralCache ["AMBIGS"])[_bacCurrentState[_bacCounter]]] = _bacCurrentState[_bacCounter];
	 }
	 
	 
	 _ancestralRecoveryCache[_bacCacheInstanceCounter] = _bacAncestralCache;
	 _bacCacheInstanceCounter = _bacCacheInstanceCounter + 1;
	 return _bacCacheInstanceCounter-1;
}

/*******************************************
	call this function the index returned by 
	_buildAncestralCache when the ancestral
	reconstruction is no longer needed

*******************************************/

	function _destroyAncestralCache (_ancID)
	{
		_ancestralRecoveryCache - _ancID;
		return 0;
	}

/*******************************************
	count subsitutions at a given site;
	returns an AVL with 
	
	["CHARS"]  - the 1xN matrix which maps index->character string
	["COUNTS"] - the NxN (N = Columns(returnValue["CHARS"]) matrix 
				 with counts of substitutions at the site
	
	pass this function the ID of ancestral cache
	and the index of the site to recover

*******************************************/

	function _substitutionsBySite (_ancID, _siteID)
	{
		if (Abs (_ancestralRecoveryCache[_ancID]))
		{
			if (_siteID >= 0 && _siteID < Columns ((_ancestralRecoveryCache[_ancID])["MATRIX"]))
			{
				_bacSiteC   = {};
				_thisColumn 		= ((_ancestralRecoveryCache[_ancID])["MATRIX"])[-1][_siteID];
				_bacSiteC ["CHARS"] = (_ancestralRecoveryCache[_ancID])["CHARS"];
				_bacSiteDim 		= Columns (_bacSiteC ["CHARS"]);
				_bacCounter 		= Rows (_thisColumn)-1;
				_bacSiteMx			= {_bacSiteDim,_bacSiteDim};
				
				for (_bacTreeIterator = 0; _bacTreeIterator < _bacCounter; _bacTreeIterator = _bacTreeIterator + 1)
				{
					_bacParentID = (((_ancestralRecoveryCache[_ancID])["TREE_AVL"])[_bacTreeIterator+1])["Parent"]-1;
					_myState	 = _thisColumn[_bacTreeIterator];
					_pState		 = _thisColumn[_bacParentID];
					_expandSubstitutionMap (_pState,_myState,_ancID,"_bacSiteMx");
				}
				
				_bacSiteC ["COUNTS"] = _bacSiteMx;
				return _bacSiteC;
			}
		}
		return {};
	}
	

/*******************************************
	count subsitutions at a given site;
	along a subset of branches. 
	
	["CHARS"]  - the 1xN matrix which maps index->character string
	["COUNTS"] - the NxN (N = Columns(returnValue["CHARS"]) matrix 
				 with counts of substitutions at the site
	
	pass this function the ID of ancestral cache
	and the index of the site to recover

*******************************************/

	function _substitutionsBySiteSubset (_ancID, _siteID, _branchSubset)
	{
		if (Abs (_ancestralRecoveryCache[_ancID]))
		{
			if (_siteID >= 0 && _siteID < Columns ((_ancestralRecoveryCache[_ancID])["MATRIX"]))
			{
				_bacSiteC   = {};
				_thisColumn 		= ((_ancestralRecoveryCache[_ancID])["MATRIX"])[-1][_siteID];
				_bacSiteC ["CHARS"] = (_ancestralRecoveryCache[_ancID])["CHARS"];
				_bacSiteDim 		= Columns (_bacSiteC ["CHARS"]);
				_bacCounter 		= Rows (_thisColumn)-1;
				_bacSiteMx			= {_bacSiteDim,_bacSiteDim};
				
				for (_bacTreeIterator = 0; _bacTreeIterator < _bacCounter; _bacTreeIterator = _bacTreeIterator + 1)
				{
					if (_branchSubset[(((_ancestralRecoveryCache[_ancID])["TREE_AVL"])[_bacTreeIterator+1])["Name"] && 1])
					{
						_bacParentID = (((_ancestralRecoveryCache[_ancID])["TREE_AVL"])[_bacTreeIterator+1])["Parent"]-1;
						_myState	 = _thisColumn[_bacTreeIterator];
						_pState		 = _thisColumn[_bacParentID];
						_expandSubstitutionMap (_pState,_myState,_ancID,"_bacSiteMx");
					}
				}
				
				_bacSiteC ["COUNTS"] = _bacSiteMx;
				return _bacSiteC;
			}
		}
		return {};
	}	

/*******************************************
	returns the reconstructed root state for a given site;
	returns an AVL with 
	
	["CHAR"]  - the root character
	["INDEX"] - the index of the root character (consisent with the data filter representation),
			  - this value will be negative if the root state is ambiguous (e.g. a gap)
	
	pass this function the ID of ancestral cache
	and the index of the site to recover

*******************************************/

	function _rootState (_ancID, _siteID)
	{
		if (Abs (_ancestralRecoveryCache[_ancID]))
		{
			if (_siteID >= 0 && _siteID < Columns ((_ancestralRecoveryCache[_ancID])["MATRIX"]))
			{
				_bacRootState   = {};
				_bacRootIndex   = (((_ancestralRecoveryCache[_ancID]) ["TREE_AVL"])[0])["Root"]-1;
				_rootStateIndex = ((_ancestralRecoveryCache[_ancID])["MATRIX"])[_bacRootIndex][_siteID];
				_bacRootState["INDEX"] = _rootStateIndex;
				if (_rootStateIndex>=0)
				{
					_bacRootState["CHAR"] = ((_ancestralRecoveryCache[_ancID]) ["CHARS"])[_rootStateIndex];
				}
				else
				{
					_bacRootState["CHAR"] = "-";	
				}
				return _bacRootState;
			}
		}
		return {};
	}	
	
/*******************************************
	map all site substitutions or character states to a tree; 
	'_scaled' != 0 will use branch lengths to draw the tree
	returns PostScript code for the image

*******************************************/
	function _mapSubstitutionsBySite (_ancID, _siteID, _scaled)
	{
		return _mapSubstitutionsBySiteAux (_ancID, _siteID, _scaled, 0);
	}

	/********************************************/

	function _mapCharactersBySite (_ancID, _siteID, _scaled)
	{
		return _mapSubstitutionsBySiteAux (_ancID, _siteID, _scaled, 1);
	}

	/********************************************/

	function _mapSNSBySite (_ancID, _siteID, _scaled)
	{
		return _mapSubstitutionsBySiteAux (_ancID, _siteID, _scaled, 2);
	}

	/********************************************/
	
	function _mapSubstitutionsBySiteAux (_ancID, _siteID, _scaled, mode)
	{
		if (Abs (_ancestralRecoveryCache[_ancID]))
		{
			if (_siteID >= 0 && _siteID < Columns ((_ancestralRecoveryCache[_ancID])["MATRIX"]))
			{
				TREE_OUTPUT_OPTIONS = {};
				_thisColumn 		= ((_ancestralRecoveryCache[_ancID])["MATRIX"])[-1][_siteID];
				_bacSiteC 		    = (_ancestralRecoveryCache[_ancID])["CHARS"];
				_bacSiteDim 		= Columns (_bacSiteC);
				_bacCounter 		= Rows (_thisColumn)-1;
				
				for (_bacTreeIterator = 0; _bacTreeIterator < _bacCounter; _bacTreeIterator = _bacTreeIterator + 1)
				{
					_bacParentID   = (((_ancestralRecoveryCache[_ancID])["TREE_AVL"])[_bacTreeIterator+1])["Parent"]-1;
					_myState	   = _thisColumn[_bacTreeIterator];
					_pState		   = _thisColumn[_bacParentID];
					_bacStateLabel = "";
					
					
					_bac_bn = (((_ancestralRecoveryCache[_ancID])["TREE_AVL"])[_bacTreeIterator+1])["Name"];
					TREE_OUTPUT_OPTIONS [_bac_bn] = {};
                    if (mode == 2)
                    {
                        (TREE_OUTPUT_OPTIONS[_bac_bn]) ["TREE_OUTPUT_BRANCH_THICKNESS"] = 1;
                    }
					if (_myState >= 0)
					{
						if (mode == 0)
						{
							 if (_myState != _pState && _pState >= 0)
							 {
								_bacStateLabel = _bacSiteC[_pState] + "->" + _bacSiteC[_myState];
							 }
						}
						else
						{
							_bacStateLabel = _bacSiteC[_myState];
							
						}
					}
					
					
					if (Abs(_bacStateLabel))
					{
						if (mode == 0)
						{
							(TREE_OUTPUT_OPTIONS[_bac_bn]) ["TREE_OUTPUT_OVER_BRANCH"] = "gsave 0.7 0.7 scale ("+_bacStateLabel+") drawletter grestore\n";
							(TREE_OUTPUT_OPTIONS[_bac_bn]) ["TREE_OUTPUT_BRANCH_LABEL"]  = "__FONT_SIZE__ 2 idiv\n__FONT_SIZE__ 3 idiv\nneg\nrmoveto\n (__NODE_NAME__) show";
							(TREE_OUTPUT_OPTIONS[_bac_bn]) ["TREE_OUTPUT_BRANCH_COLOR"] = {{1,0,0}};
						}
						else
						{
							if (mode == 1)
							{
								(TREE_OUTPUT_OPTIONS[_bac_bn]) ["TREE_OUTPUT_BRANCH_TLABEL"] = _bacStateLabel;	
							}
							else
							{
								if (mode == 2)
								{
									haveS 	= _OBSERVED_S_		[_pState][_myState];
									haveNS	= _OBSERVED_NS_	[_pState][_myState];
									if (haveNS >= 1)
									{
										if (haveS >= 1)
										{
											(TREE_OUTPUT_OPTIONS[_bac_bn]) ["TREE_OUTPUT_BRANCH_COLOR"] = {{1,0.7,0}};
										}
										else
										{
											(TREE_OUTPUT_OPTIONS[_bac_bn]) ["TREE_OUTPUT_BRANCH_COLOR"] = {{1,0,0}};										
										}
									}
									else
									{
										if (haveS >= 1)
										{
											(TREE_OUTPUT_OPTIONS[_bac_bn]) ["TREE_OUTPUT_BRANCH_COLOR"] = {{0,0,1}};																				
										}
									}
									
									(TREE_OUTPUT_OPTIONS[_bac_bn]) ["TREE_OUTPUT_BRANCH_THICKNESS"] = 1 + haveS + haveNS;
								}
							}
						}
					}
				}
				
                if (mode == 2)
                {
                    (TREE_OUTPUT_OPTIONS[(((_ancestralRecoveryCache[_ancID])["TREE_AVL"])[_bacTreeIterator+1])["Name"]]) = {"TREE_OUTPUT_BRANCH_THICKNESS" : 1};
                }

				_bacTreeString 	  = PostOrderAVL2StringDL ((_ancestralRecoveryCache[_ancID])["TREE_AVL"], _scaled);
				Tree _bacTempTree = _bacTreeString;
				_bac_bn = "";
				_bac_bn * 128;
				_bac_bn * "/drawletter {5 5 rmoveto 1 copy false charpath pathbbox 2 index 3 sub sub exch 3 index 3 sub sub exch  0.85 setgray 4 copy rectfill 0 setgray  3 index 3 index currentlinewidth 0.5 setlinewidth 7 3 roll rectstroke setlinewidth exch 1.5 add exch 1.5 add moveto show} def\n";
				
				if (_scaled)
				{
					_bac_bn * (PSTreeString (_bacTempTree,"STRING_SUPPLIED_LENGTHS",{{-1,-1}}));				
				}
				else
				{
					_bac_bn * (PSTreeString (_bacTempTree,"",{{-1,-1}}));
				}
				_bac_bn * 0;
				
				return _bac_bn;
			}
		}
		return "";
	}

/*******************************************
	map all site substitutions (or character states) to a tree; 
	'_scaled' != 0 will use branch lengths to draw the tree
	returns Newick code for the annotated tree

*******************************************/
	function _mapSubstitutionsBySiteNewick (_ancID, _siteID, _scaled)
	{
		return _mapSubstitutionsBySiteNewickAux (_ancID, _siteID, _scaled, 0);
	}

	function _mapCharactersBySiteNewick (_ancID, _siteID, _scaled)
	{
		return _mapSubstitutionsBySiteNewickAux (_ancID, _siteID, _scaled, 1);
	}

	function _mapSubstitutionsBySiteNewickAux (_ancID, _siteID, _scaled, mode)
	{
		if (Abs (_ancestralRecoveryCache[_ancID]))
		{
			if (_siteID >= 0 && _siteID < Columns ((_ancestralRecoveryCache[_ancID])["MATRIX"]))
			{
				_thisColumn 		= ((_ancestralRecoveryCache[_ancID])["MATRIX"])[-1][_siteID];
				_bacSiteC 		    = (_ancestralRecoveryCache[_ancID])["CHARS"];
				_bacSiteDim 		= Columns (_bacSiteC);
				_bacCounter 		= Rows (_thisColumn)-1;
				
				for (_bacTreeIterator = 0; _bacTreeIterator < _bacCounter; _bacTreeIterator = _bacTreeIterator + 1)
				{
					_bacParentID   = (((_ancestralRecoveryCache[_ancID])["TREE_AVL"])[_bacTreeIterator+1])["Parent"]-1;
					_myState	   = _thisColumn[_bacTreeIterator];
					_pState		   = _thisColumn[_bacParentID];
					_bacStateLabel = "";
					
					
					if (_myState >= 0)
					{
						if (mode == 0)
						{
							 if (_myState != _pState && _pState >= 0)
							 {
								_bacStateLabel = _bacSiteC[_pState] + "->" + _bacSiteC[_myState];
							 }
						}
						else
						{
							_bacStateLabel = _bacSiteC[_myState];
						}
					}
					
					(((_ancestralRecoveryCache[_ancID])["TREE_AVL"])[_bacTreeIterator+1])["SubLabel"] = _bacStateLabel;
				}
				
				return 	 PostOrderAVL2StringAnnotate ((_ancestralRecoveryCache[_ancID])["TREE_AVL"], _scaled, "SubLabel");
			}
		}
		return "";
	}

	
/*******************************************
	
	_filterDimensions returns a {{sites,branches}} matrix
	
*******************************************/

	function _filterDimensions (_ancID)
	{
		_result = {1,2};
		if (Abs (_ancestralRecoveryCache[_ancID]))
		{
			_sites    = Columns ((_ancestralRecoveryCache[_ancID])["MATRIX"]);
			_branches = Rows	((_ancestralRecoveryCache[_ancID])["MATRIX"])-1;
			_result [0] = _sites;
			_result [1] = _branches;
		}
		return _result;
	}

/*******************************************
	
	_countSubstitutionsByBranchSite returns a 
	binary vector (one entry per branch) which
	contains '1' if the branch has substitutions
	defined by '_filter' (a CxC matrix)

*******************************************/

	function _countSubstitutionsByBranchSite (_ancID, _siteID, _filter)
	{
		_result = {0,1};
		if (Abs (_ancestralRecoveryCache[_ancID])) {
			if (_siteID >= 0 && _siteID < Columns ((_ancestralRecoveryCache[_ancID])["MATRIX"])) {
				TREE_OUTPUT_OPTIONS = {};
				_thisColumn 		= ((_ancestralRecoveryCache[_ancID])["MATRIX"])[-1][_siteID];
				_bacSiteC 		    = (_ancestralRecoveryCache[_ancID])["CHARS"];
				_bacSiteDim 		= Columns (_bacSiteC);
				_bacCounter 		= Rows (_thisColumn)-1;
				_result				= {_bacCounter,1};

				for (_bacTreeIterator = 0; _bacTreeIterator < _bacCounter; _bacTreeIterator += 1) {
					_bacParentID   = (((_ancestralRecoveryCache[_ancID])["TREE_AVL"])[_bacTreeIterator+1])["Parent"]-1;
					_myState	   = _thisColumn[_bacTreeIterator];
					_pState		   = _thisColumn[_bacParentID];
					_bacSiteMx			= {_bacSiteDim,_bacSiteDim};
					_expandSubstitutionMap (_pState,_myState,_ancID,"_bacSiteMx");
					_result[_bacTreeIterator] = (+(_bacSiteMx$_filter))>=1;
				}
				
			}
		}
		return _result;
	}



/*******************************************
	
	_tabulateSubstitutionsAtSiteByBranch returns
	a dictionary with keys = branch names and
	values = counts of synonymous and non-synonymous
	substitutions.

*******************************************/

	function _tabulateSubstitutionsAtSiteByBranch (_ancID, _siteID)
	{
		_result = {};
		if (Abs (_ancestralRecoveryCache[_ancID]))
		{
			if (_siteID >= 0 && _siteID < Columns ((_ancestralRecoveryCache[_ancID])["MATRIX"]))
			{
				_thisColumn 		= ((_ancestralRecoveryCache[_ancID])["MATRIX"])[-1][_siteID];
				_bacSiteC 		    = (_ancestralRecoveryCache[_ancID])["CHARS"];
				_bacSiteDim 		= Columns (_bacSiteC);
				_bacCounter 		= Rows (_thisColumn)-1;


				for (_bacTreeIterator = 0; _bacTreeIterator < _bacCounter; _bacTreeIterator = _bacTreeIterator + 1)
				{
					_bacParentID   = (((_ancestralRecoveryCache[_ancID])["TREE_AVL"])[_bacTreeIterator+1])["Parent"]-1;
					_myState	   = _thisColumn[_bacTreeIterator];
					_pState		   = _thisColumn[_bacParentID];
					
					if (_pState >= 0 && _myState >=0) {
					    haveS 	= _OBSERVED_S_	[_pState][_myState];
					    haveNS	= _OBSERVED_NS_	[_pState][_myState];
					} else {
					    haveS = 0;
					    haveNS = 0;
					}
					_result [(((_ancestralRecoveryCache[_ancID])["TREE_AVL"])[_bacTreeIterator+1])["Name"]] = {{haveS__,haveNS__}};
                }
				
			}
		}
		return _result;
	}
	
	
/*******************************************/
	
	function _expandSubstitutionMap (_state1, _state2, _ancID, _resultMatrix&)
	{
		if (_state1 >= 0 && _state2 >= 0)
		/* both unique states */
		{
			_resultMatrix[_state1][_state2] = _resultMatrix[_state1][_state2]+1;
		}
		else
		{
			if (_state1 != (-1) && _state2 != (-1))
			/* at least one is NOT a gap state: otherwise nothing to map */
			{
				_thisDim = Rows (_resultMatrix);
				if (_state1 >= 0)
				{
					_vec1 = {_thisDim,1};
					_vec1 [_state1] = 1;
				}
				else
				{
					_vec1 = ((_ancestralRecoveryCache[_ancID])["AMBIGS"])[-_state1-2];
				}
				if (_state2 >= 0)
				{
					_vec2 = {_thisDim,1};
					_vec2 [_state2] = 1;
				}
				else
				{
					_vec2 = ((_ancestralRecoveryCache[_ancID])["AMBIGS"])[-_state2-2];
				}
				_vec1 = _vec1 * Transpose(_vec2);
				_vec2 = {_thisDim,1}["1"];
				_vec2 = (Transpose(_vec2)*_vec1*_vec2)[0];
				_resultMatrix = _resultMatrix + _vec1 * (1/_vec2); /* 1/K(=all possible pairwise resolutions) to each possible resolution */
			}
			
		}
		return 0;
	}
	
/*******************************************

*******************************************/
	
	function _convertSubstitutionToCharacters (_state1, _state2, _ancID, _resultMatrix&)
	{
		if (_state1 >= 0 && _state2 >= 0)
		/* both unique states */
		{
			_resultMatrix[_state1][_state2] = _resultMatrix[_state1][_state2]+1;
		}
		else
		{
			if (_state1 != (-1) && _state2 != (-1))
			/* at least one is NOT a gap state: otherwise nothing to map */
			{
				_thisDim = Rows (_resultMatrix);
				if (_state1 >= 0)
				{
					_vec1 = {_thisDim,1};
					_vec1 [_state1] = 1;
				}
				else
				{
					_vec1 = ((_ancestralRecoveryCache[_ancID])["AMBIGS"])[-_state1-2];
				}
				if (_state2 >= 0)
				{
					_vec2 = {_thisDim,1};
					_vec2 [_state2] = 1;
				}
				else
				{
					_vec2 = ((_ancestralRecoveryCache[_ancID])["AMBIGS"])[-_state2-2];
				}
				_vec1 = _vec1 * Transpose(_vec2);
				_vec2 = {_thisDim,1}["1"];
				_vec2 = (Transpose(_vec2)*_vec1*_vec2)[0];
				_resultMatrix = _resultMatrix + _vec1 * (1/_vec2); /* 1/K(=all possible pairwise resolutions*/ to each possible resolution */
			}
			
		}
		return 0;
	}	
