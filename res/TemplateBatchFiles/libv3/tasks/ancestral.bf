RequireVersion("2.3.11");

LoadFunctionLibrary("../UtilityFunctions.bf");

/** @module ancestral */

ancestral._ancestralRecoveryCache = {};

/*******************************************

Example Structure

{
 "DIMENSIONS":{
   "SITES":1,
   "SEQUENCES":10,
   "BRANCHES":17,
   "CHARS":61
  },
 "CHARS":{
  {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"}
  },
 "MATRIX":{
  {-1}
  {-1}
  {-1}
  {-1}
  {-1}
  {-1}
  {0}
  {0}
  {0}
  {0}
  {32}
  {0}
  {0}
  {0}
  {-1}
  {-1}
  {0}
  },
 "TREE_AVL":{
   "1":{
     "Name":"PIG",
     "Length":-1,
     "Depth":4,
     "Parent":3
    },
   "2":{
     "Name":"COW",
     "Length":-1,
     "Depth":4,
     "Parent":3
    },
    .
    .
    .
   "17":{
     "Name":"Node0",
     "Length":-1,
     "Depth":0,
     "Children":{
       "0":14,
       "1":15,
       "2":16
      }
    },
   "0":{
     "Name":"meme.site_tree_bsrel",
     "Root":17
    }
  },
 "AMBIGS":{
  },
 "MAPPING":{
   "-1":"---",
   "0":"AAA",
   "32":"GAA"
  }
}


*******************************************/


ancestral._bacCacheInstanceCounter = 0;

/**
 * Builds ancestral states
 * @name ancestral.build
 * @param {Number} _lfID - the likelihood function ID
 * @param {Number} _lfComponentID -
 * @param {options} options -
 * @returns an integer index to reference
 * the opaque structure for subsequent operations
 * @example
 */
function ancestral.build (_lfID, _lfComponentID, options) {
    return ancestral._buildAncestralCacheInternal(_lfID, _lfComponentID, options["sample"], options["marginal"]);
}


/**
 * internal function to do the work of ancestral.build
 * @name ancestral._buildAncestralCacheInternal
 * @private
 * @param {Number} _lfID - the likelihood function ID
 * @param {Number} _lfComponentID -
 * @param {options} options -
 * @returns an integer index to reference
 * the opaque structure for subsequent operations
 */
lfunction ancestral._buildAncestralCacheInternal(_lfID, _lfComponentID, doSample, doMarginal) {


    /* 1; grab the information AVL from the likelihood function */

    GetString(_bac_lfInfo, ^_lfID, -1);
    if (Columns(_bac_lfInfo["Trees"]) <= _lfComponentID) {
        return None;
    }

    _bac_treeID = (_bac_lfInfo["Trees"])[_lfComponentID];
    _bac_filterID = (_bac_lfInfo["Datafilters"])[_lfComponentID];

    /* 2; construct a temporary likelihood function with
       the tree and the filter; and an AVL representation of the tree */

    if (doSample) {
        DataSet _bac_ancDS = SampleAncestors( ^ _lfID, {
            {
                _lfComponentID
            }
        });
    } else {
        if (doMarginal) {
            DataSet _bac_ancDS = ReconstructAncestors( ^ _lfID, {
                {
                    _lfComponentID
                }
            }, MARGINAL);
        } else {
            DataSet _bac_ancDS = ReconstructAncestors( ^ _lfID, {
                {
                    _lfComponentID
                }
            });
        }
    }

    _bac_tree_avl = ( ^ _bac_treeID) ^ 0;
    
    GetString(_bacSequenceNames, ^ _bac_filterID, -1);
    

    /* 3; obtain ID->string mapping for the datafilter;
    	  also deduce how many chars/state there are */

    GetDataInfo(_bacCharHandles, ^ _bac_filterID, "CHARACTERS");
    GetDataInfo(_bacCharProperties, ^ _bac_filterID, "PARAMETERS");
    GetDataInfo(_bacFilterPatternMap, ^ _bac_filterID);
    _bacFilterDimension = Columns(_bacCharHandles);
    _bacCharsPerState = _bacCharProperties["ATOM_SIZE"];
    DataSetFilter _bacAF = CreateFilter(_bac_ancDS, _bacCharsPerState, "", "", _bacCharProperties["EXCLUSIONS"]);
    
    GetString(_bacAncestralNames, _bacAF, -1);
    GetDataInfo(_bacAncestralPatternMap, _bacAF);


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

    _bacBranchCount = Abs(_bac_tree_avl) - 1;
    _bacMatrixOfResolutions = {
        _bacBranchCount,
        _bacAF.sites
    };
    
    _bacHandledResolutions = {};
    _bacHandledResolutionsAmbig = {};

    reverse_mapping = {};
    
    for (i,_bacCounter; in; _bacCharHandles) {
        _bacHandledResolutions[_bacCounter] = i + 1;
        reverse_mapping [i] = _bacCounter;
    }
    
    
    

    /* map where sequences are in the filter vs where they are in the
	tree structure */

    _bacMapTreeNodeToDF = {};
    _bacTreeAVLOrder = {};

    for (_bacCounter = 1; _bacCounter <= _bacBranchCount; _bacCounter += 1) {
        _bacNodeName = (_bac_tree_avl[_bacCounter])["Name"];
        
         if (Abs ((_bac_tree_avl[_bacCounter])["Children"]) == 0) {
            _bacNodeName = _bacNodeName && 1;
         }
        _bacTreeAVLOrder[_bacNodeName] = _bacCounter;
    }


    _bacFilterSequenceCount = Columns(_bacSequenceNames);
    for (_bacCounter = 0; _bacCounter < _bacFilterSequenceCount; _bacCounter += 1) {
        _bacMapTreeNodeToDF[_bacTreeAVLOrder[_bacSequenceNames[_bacCounter] && 1] - 1] = _bacCounter;
    }

    for (_bacCounter = 0; _bacCounter < Columns(_bacAncestralNames); _bacCounter += 1) {
        _bacMapTreeNodeToDF[_bacTreeAVLOrder[_bacAncestralNames[_bacCounter]] - 1] = _bacCounter + _bacFilterSequenceCount;
    }

    /* make some auxiliary variables */

    _bacUnitRow = {
        1,
        Columns(_bacCharHandles)
    }["1"];
    _bacSequenceRow = {
        1,
        Columns(_bacCharHandles)
    }["_MATRIX_ELEMENT_COLUMN_"];


    /* loop over branches (rows) */
    for (_bacBranchCounter = 1; _bacBranchCounter <= _bacBranchCount; _bacBranchCounter += 1) {
        _bacRowIndex = _bacMapTreeNodeToDF[_bacBranchCounter - 1];
        if (_bacRowIndex < _bacFilterSequenceCount) {
            /* an extant sequence */
            GetDataInfo(_bacSequenceString, ^ _bac_filterID, _bacRowIndex);
            GetString(_bacSequenceName, ^ _bac_filterID, _bacRowIndex);
        } else {
            /* an ancestral sequence */
            GetDataInfo(_bacSequenceString, _bacAF, _bacRowIndex - _bacFilterSequenceCount);
            GetString(_bacSequenceName, _bacAF, _bacRowIndex - _bacFilterSequenceCount);
        }

        /* loop over sites (columns) */
        for (_bacSiteCounter = 0; _bacSiteCounter < _bacAF.sites; _bacSiteCounter += 1) {
            _bacCurrentState = _bacSequenceString[_bacSiteCounter * _bacCharsPerState][(1 + _bacSiteCounter) * _bacCharsPerState - 1];
            _bacCurrentIndex = _bacHandledResolutions[_bacCurrentState];
            if (_bacCurrentIndex > 0) {
                _bacMatrixOfResolutions[_bacBranchCounter - 1][_bacSiteCounter] = _bacCurrentIndex - 1;
            } else {
                if (_bacCurrentIndex < 0) {
                    _bacMatrixOfResolutions[_bacBranchCounter - 1][_bacSiteCounter] = _bacCurrentIndex;
                } else {
                    /* handle the resolution */
                    if (_bacRowIndex < _bacFilterSequenceCount) {
                        GetDataInfo(_bacCharState, ^ _bac_filterID, _bacRowIndex, _bacFilterPatternMap[_bacSiteCounter]);
                    } else {
                        GetDataInfo(_bacCharState, _bacAF, _bacRowIndex - _bacFilterSequenceCount, _bacAncestralPatternMap[_bacSiteCounter]);
                    }
                    _bacResolutionCount = +_bacCharState;
                    if (_bacResolutionCount == Columns(_bacCharHandles)) {
                            /* gap/full ambig */
                        if ( (reverse_mapping/(-1)) == FALSE) {
                            _bacHandledResolutions[_bacCurrentState] = -1;
                            reverse_mapping[-1] = _bacCurrentState;
                        }
                        _bacMatrixOfResolutions[_bacBranchCounter - 1][_bacSiteCounter] = -1;
                    } else {
                        _bacHandledResolutions[_bacCurrentState] = -2 - Abs(_bacHandledResolutionsAmbig);
                        _bacHandledResolutionsAmbig + _bacCharState;
                        _bacMatrixOfResolutions[_bacBranchCounter - 1][_bacSiteCounter] = _bacHandledResolutions[_bacCurrentState];
                        reverse_mapping [_bacHandledResolutions[_bacCurrentState]] = _bacCurrentState;
                    }
                 }
            }
        }
    }

    

    return {
        "DIMENSIONS": {
            "SITES": _bacAF.sites,
            "SEQUENCES": _bacFilterSequenceCount,
            "BRANCHES": _bacBranchCount,
            "CHARS": _bacFilterDimension
        },
        "CHARS": _bacCharHandles,
        "MATRIX": _bacMatrixOfResolutions,
        "TREE_AVL": _bac_tree_avl,
        "AMBIGS": _bacHandledResolutionsAmbig,
        "MAPPING": reverse_mapping
            //"AMBIGS_REV": utility.SwapKeysAndValues (_bacHandledResolutionsAmbig)
    };
}

/*******************************************/
/**
 * @name ancestral.ComputeSubstitutionCounts
   branch filter helper function
*/

lfunction ancestral._branch_filter_helper (ancestral_data, branch_filter, selected_branches, selected_branch_names) {
    coordinates = {{k-1, parent-1}};

    for (k = 1; k < Abs(ancestral_data["TREE_AVL"]); k+=1) {
        parent = ((ancestral_data["TREE_AVL"])[k])["Parent"];

        if (parent) {
            node_name = ((ancestral_data["TREE_AVL"])[k])["Name"];
            if (None != branch_filter) {
                if (Type (branch_filter) == "AssociativeList") {
                    if (branch_filter / node_name == FALSE) {
                        continue;
                    }
                } else {
                    if (Call (branch_filter, ((ancestral_data["TREE_AVL"])[k])) == FALSE) {
                        continue;
                    }
                }
            }
            selected_branches + Eval (coordinates);
            selected_branch_names + node_name;
        }
    }

    return Abs (selected_branches);
}

/*******************************************/


lfunction ancestral._select_internal (node) {
    return Abs (node["Children"]);
}

/*******************************************/
/**
 * @name ancestral.Sequences
 * @param {Dictionary} ancestral_data - the dictionary returned by ancestral.build
 
 * @returns
        {
         "Node Name" :  {String} inferred ancestral sequence,
        }

 */

lfunction ancestral.Sequences (ancestral_data) {
    selected_branches       = {};
    selected_branch_names   = {};
    
    branches =  ancestral._branch_filter_helper (ancestral_data, "ancestral._select_internal", selected_branches, selected_branch_names) + 1;    
    sites  = (ancestral_data["DIMENSIONS"])["SITES"];
    result = {};
    selected_branches + {{(Abs(ancestral_data["TREE_AVL"])-2),0}};
    selected_branch_names + "root";

    for (b = 0; b < branches; b += 1) {
        self   = (selected_branches[b])[0];    
        seq_string = ""; seq_string * sites;
        
        for (s = 0; s < sites; s += 1) {
            own_state    = (ancestral_data["MATRIX"])[self][s];
            seq_string * (ancestral_data["CHARS"])[own_state];
            
        }   
        seq_string * 0;
        result[selected_branch_names[b]] = seq_string;
    }
    
    return result;
}



/*******************************************/
/**
 * @name ancestral.ComputeSubstitutionCounts
 * @param {Dictionary} ancestral_data - the dictionary returned by ancestral.build
 * @param {Dictionary/Function/None} branch_filter - now to determine the subset of branches to count on
          None -- all branches
          Dictionary -- all branches that appear as keys in this dict
          Function -- all branches on which the function (called with branch name) returns 1

 * @param {Function/None} substitution_filter - how to decide if the substitution should count
          None -- different characters (except anything vs a gap) yields a 1
          Function -- callback (state1, state2, ancestral_data) will return the value

 * @param {Function/None} site_filter - how to decide which sites will be kept
          None     -- at least one substitution
          Function -- callback (substitution vector) will T/F

 * @returns
        {
         "Branches" :  {Matrix Nx1} names of selected branches,
         "Sites"    :  {Matrix Sx1} indices of sites passing filter,
         "Counts"   :  {Matrix NxS} of substitution counts,
        }

        N = number of selected branches
        S = number of sites passing filter

 */




lfunction ancestral.ComputeSubstitutionCounts (ancestral_data, branch_filter, substitution_filter, site_filter) {
    selected_branches       = {};
    selected_branch_names   = {};

    branches =  ancestral._branch_filter_helper (ancestral_data, branch_filter, selected_branches, selected_branch_names);

    sites  = (ancestral_data["DIMENSIONS"])["SITES"];
    counts = {branches,sites};
    retained_sites = {};


    for (b = 0; b < branches; b += 1) {
        self   = (selected_branches[b])[0];
        parent = (selected_branches[b])[1];
        for (s = 0; s < sites; s += 1) {
            own_state    = (ancestral_data["MATRIX"])[self][s];
            parent_state = (ancestral_data["MATRIX"])[parent][s];
            if (None == substitution_filter) {
                counts[b][s] = (own_state != parent_state) && (own_state != -1) && (parent_state != -1);
            } else {
                counts[b][s] = Call (substitution_filter, own_state, parent_state, ancestral_data);
            }
        }
    }

    for (s = 0; s < sites; s += 1) {
        site_counts = counts [-1][s];
        if (None == site_filter) {
            if (+site_counts == 0) {
                continue;
            }
        } else {
            if (Call(site_filter, site_counts) == FALSE) {
                continue;
            }
        }
        retained_sites + s;
    }

    retained_site_count = Abs (retained_sites);
    retained_counts = {branches, retained_site_count};

    for (s = 0; s < retained_site_count; s+=1) {
        full_index = retained_sites [s];
        for (b = 0; b < branches; b += 1) {
            retained_counts[b][s] = counts[b][full_index];
        }
    }

    counts = None;

    return  {
             "Branches"  : selected_branch_names,
             "Sites"     : retained_sites,
             "Counts"    : retained_counts
            };

}

/*******************************************
 **
 * @name ancestral.ComputeSubstitutionBySite
 * tabulates all substitutions occurring at a specified site, possibly restricted to a
 * subset of branches

 * @param {Dictionary} ancestral_data - the dictionary returned by ancestral.build
 * @param {Number} site - the 0-based index of a site

 * @param {Dictionary/Function/None} branch_filter - now to determine the subset of branches to count on
          None -- all branches
          Dictionary -- all branches that appear as keys in this dict
          Function -- all branches on which the function (called with branch name) returns 1

 * @returns
        {
         terms.trees.branches :  {Matrix Nx1} names of selected branches,
         terms.substitutions   :  {Dictionary} of substitution counts;
                              will look like "FROM" -> "TO" -> List of branches where this type of substitution occured
         terms.substitutions:{
               "N":{
                 "H":{
                   "0":"COW"
                  }
                },
               "K":{
                 "N":{
                   "0":"Node5"
                  },
                 "T":{
                   "0":"HORSE"
                  }
                },
               "R":{
                 "K":{
                   "0":"Node3"
                  },
                 "P":{
                   "0":"MOUSE"
                  }
                }
              }
        }

        N = number of selected branches
 */

/*******************************************/


lfunction ancestral.ComputeSubstitutionBySite (ancestral_data, site, branch_filter) {
    selected_branches       = {};
    selected_branch_names   = {};

    branches = ancestral._branch_filter_helper (ancestral_data, branch_filter, selected_branches, selected_branch_names);

    result   = {};


    for (b = 0; b < branches; b += 1) {
        self   = (selected_branches[b])[0];
        parent = (selected_branches[b])[1];
        own_state    = (ancestral_data["MATRIX"])[self][site];
        parent_state = (ancestral_data["MATRIX"])[parent][site];

        if  ((own_state != parent_state) && (own_state != -1) && (parent_state != -1)) {
            own_state = (ancestral_data["CHARS"])[own_state];
            parent_state = (ancestral_data["CHARS"])[parent_state];
            utility.EnsureKey (result, parent_state);
            utility.EnsureKey (result[parent_state], own_state);
            ((result[parent_state])[own_state]) + selected_branch_names[b];
        }
    }


    return  {
             ^"terms.trees.branches"  : selected_branch_names,
             ^"terms.substitutions"   : result
            };

}

/*******************************************
 **
 * @name ancestral.ComputeSiteComposition
 * computes the counts of individual characters (including ancestral states)
 * at a specific site, possibly filtered to a subset of branches
 * @param {Dictionary} ancestral_data - the dictionary returned by ancestral.build
 * @param {Number} site - the 0-based index of a site

 * @param {Dictionary/Function/None} branch_filter - now to determine the subset of branches to count on
          None -- all branches
          Dictionary -- all branches that appear as keys in this dict
          Function -- all branches on which the function (called with branch name) returns 1

 * @returns
        {
         terms.trees.branches :  {Matrix Nx1} names of selected branches,
         terms.data.composition   :  {Dictionary} of character counts;
                              will look like

        }

        N = number of selected branches
 */

/*******************************************/


lfunction ancestral.ComputeSiteComposition (ancestral_data, site, branch_filter) {
    selected_branches       = {};
    selected_branch_names   = {};


    branches =  ancestral._branch_filter_helper (ancestral_data, branch_filter, selected_branches, selected_branch_names);

    result   = {};


    for (b = 0; b < branches; b += 1) {
        self   = (selected_branches[b])[0];
        own_state    = (ancestral_data["MATRIX"])[self][site];

        if (own_state >= 0) {
            own_state    = (ancestral_data["CHARS"])[own_state];
            result [own_state] += 1;
        } else {
            if (own_state < (-1)) {
                resolution = (ancestral_data["AMBIGS"])[-own_state - 2];
                resolution = (resolution["(_MATRIX_ELEMENT_ROW_+1)*_MATRIX_ELEMENT_VALUE_"])[resolution];
                count = 1/utility.Array1D (resolution);
                utility.ForEach (resolution, "_value_", "`&result`[(`&ancestral_data`['CHARS'])[_value_-1]] += `&count`;");
            }
        }

    }


    return  {
             ^"terms.trees.branches"  : selected_branch_names,
             ^"terms.data.composition"   : result
            };

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

function _substitutionsBySiteSubset(_ancID, _siteID, _branchSubset) {
    if (Abs_ancID) {
        if (_siteID >= 0 && _siteID < Columns(_ancID["MATRIX"])) {
            _bacSiteC = {};
            _thisColumn = (_ancID["MATRIX"])[-1][_siteID];
            _bacSiteC["CHARS"] = _ancID["CHARS"];
            _bacSiteDim = Columns(_bacSiteC["CHARS"]);
            _bacCounter = Rows(_thisColumn) - 1;
            _bacSiteMx = {
                _bacSiteDim,
                _bacSiteDim
            };

            for (_bacTreeIterator = 0; _bacTreeIterator < _bacCounter; _bacTreeIterator = _bacTreeIterator + 1) {
                if (_branchSubset[((_ancID["TREE_AVL"])[_bacTreeIterator + 1])["Name"] && 1]) {
                    _bacParentID = ((_ancID["TREE_AVL"])[_bacTreeIterator + 1])["Parent"] - 1;
                    _myState = _thisColumn[_bacTreeIterator];
                    _pState = _thisColumn[_bacParentID];
                    _expandSubstitutionMap(_pState, _myState, _ancID, "_bacSiteMx");
                }
            }

            _bacSiteC["COUNTS"] = _bacSiteMx;
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

function _rootState(_ancID, _siteID) {
    if (Abs_ancID) {
        if (_siteID >= 0 && _siteID < Columns(_ancID["MATRIX"])) {
            _bacRootState = {};
            _bacRootIndex = ((_ancID["TREE_AVL"])[0])["Root"] - 1;
            _rootStateIndex = (_ancID["MATRIX"])[_bacRootIndex][_siteID];
            _bacRootState["INDEX"] = _rootStateIndex;
            if (_rootStateIndex >= 0) {
                _bacRootState["CHAR"] = (_ancID["CHARS"])[_rootStateIndex];
            } else {
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
function _mapSubstitutionsBySite(_ancID, _siteID, _scaled) {
    return _mapSubstitutionsBySiteAux(_ancID, _siteID, _scaled, 0);
}

/********************************************/

function _mapCharactersBySite(_ancID, _siteID, _scaled) {
    return _mapSubstitutionsBySiteAux(_ancID, _siteID, _scaled, 1);
}

/********************************************/

function _mapSNSBySite(_ancID, _siteID, _scaled) {
    return _mapSubstitutionsBySiteAux(_ancID, _siteID, _scaled, 2);
}

/********************************************/

function _mapSubstitutionsBySiteAux(_ancID, _siteID, _scaled, mode) {
    if (Abs_ancID) {
        if (_siteID >= 0 && _siteID < Columns(_ancID["MATRIX"])) {
            TREE_OUTPUT_OPTIONS = {};
            _thisColumn = (_ancID["MATRIX"])[-1][_siteID];
            _bacSiteC = _ancID["CHARS"];
            _bacSiteDim = Columns(_bacSiteC);
            _bacCounter = Rows(_thisColumn) - 1;

            for (_bacTreeIterator = 0; _bacTreeIterator < _bacCounter; _bacTreeIterator = _bacTreeIterator + 1) {
                _bacParentID = ((_ancID["TREE_AVL"])[_bacTreeIterator + 1])["Parent"] - 1;
                _myState = _thisColumn[_bacTreeIterator];
                _pState = _thisColumn[_bacParentID];
                _bacStateLabel = "";


                _bac_bn = ((_ancID["TREE_AVL"])[_bacTreeIterator + 1])["Name"];
                TREE_OUTPUT_OPTIONS[_bac_bn] = {};
                if (mode == 2) {
                    (TREE_OUTPUT_OPTIONS[_bac_bn])["TREE_OUTPUT_BRANCH_THICKNESS"] = 1;
                }
                if (_myState >= 0) {
                    if (mode == 0) {
                        if (_myState != _pState && _pState >= 0) {
                            _bacStateLabel = _bacSiteC[_pState] + "->" + _bacSiteC[_myState];
                        }
                    } else {
                        _bacStateLabel = _bacSiteC[_myState];

                    }
                }


                if (Abs(_bacStateLabel)) {
                    if (mode == 0) {
                        (TREE_OUTPUT_OPTIONS[_bac_bn])["TREE_OUTPUT_OVER_BRANCH"] = "gsave 0.7 0.7 scale (" + _bacStateLabel + ") drawletter grestore\n";
                        (TREE_OUTPUT_OPTIONS[_bac_bn])["TREE_OUTPUT_BRANCH_LABEL"] = "__FONT_SIZE__ 2 idiv\n__FONT_SIZE__ 3 idiv\nneg\nrmoveto\n (__NODE_NAME__) show";
                        (TREE_OUTPUT_OPTIONS[_bac_bn])["TREE_OUTPUT_BRANCH_COLOR"] = {
                            {
                                1,
                                0,
                                0
                            }
                        };
                    } else {
                        if (mode == 1) {
                            (TREE_OUTPUT_OPTIONS[_bac_bn])["TREE_OUTPUT_BRANCH_TLABEL"] = _bacStateLabel;
                        } else {
                            if (mode == 2) {
                                haveS = _OBSERVED_S_[_pState][_myState];
                                haveNS = _OBSERVED_NS_[_pState][_myState];
                                if (haveNS >= 1) {
                                    if (haveS >= 1) {
                                        (TREE_OUTPUT_OPTIONS[_bac_bn])["TREE_OUTPUT_BRANCH_COLOR"] = {
                                            {
                                                1,
                                                0.7,
                                                0
                                            }
                                        };
                                    } else {
                                        (TREE_OUTPUT_OPTIONS[_bac_bn])["TREE_OUTPUT_BRANCH_COLOR"] = {
                                            {
                                                1,
                                                0,
                                                0
                                            }
                                        };
                                    }
                                } else {
                                    if (haveS >= 1) {
                                        (TREE_OUTPUT_OPTIONS[_bac_bn])["TREE_OUTPUT_BRANCH_COLOR"] = {
                                            {
                                                0,
                                                0,
                                                1
                                            }
                                        };
                                    }
                                }

                                (TREE_OUTPUT_OPTIONS[_bac_bn])["TREE_OUTPUT_BRANCH_THICKNESS"] = 1 + haveS + haveNS;
                            }
                        }
                    }
                }
            }

            if (mode == 2) {
                (TREE_OUTPUT_OPTIONS[((_ancID["TREE_AVL"])[_bacTreeIterator + 1])["Name"]]) = {
                    "TREE_OUTPUT_BRANCH_THICKNESS": 1
                };
            }

            _bacTreeString = PostOrderAVL2StringDL(_ancID["TREE_AVL"], _scaled);
            Tree _bacTempTree = _bacTreeString;
            _bac_bn = "";
            _bac_bn * 128;
            _bac_bn * "/drawletter {5 5 rmoveto 1 copy false charpath pathbbox 2 index 3 sub sub exch 3 index 3 sub sub exch  0.85 setgray 4 copy rectfill 0 setgray  3 index 3 index currentlinewidth 0.5 setlinewidth 7 3 roll rectstroke setlinewidth exch 1.5 add exch 1.5 add moveto show} def\n";

            if (_scaled) {
                _bac_bn * (PSTreeString(_bacTempTree, "STRING_SUPPLIED_LENGTHS", {
                    {
                        -1, -1
                    }
                }));
            } else {
                _bac_bn * (PSTreeString(_bacTempTree, "", {
                    {
                        -1, -1
                    }
                }));
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
function _mapSubstitutionsBySiteNewick(_ancID, _siteID, _scaled) {
    return _mapSubstitutionsBySiteNewickAux(_ancID, _siteID, _scaled, 0);
}

function _mapCharactersBySiteNewick(_ancID, _siteID, _scaled) {
    return _mapSubstitutionsBySiteNewickAux(_ancID, _siteID, _scaled, 1);
}

function _mapSubstitutionsBySiteNewickAux(_ancID, _siteID, _scaled, mode) {
    if (Abs(_ancID)) {
        if (_siteID >= 0 && _siteID < Columns(_ancID["MATRIX"])) {
            _thisColumn = (_ancID["MATRIX"])[-1][_siteID];
            _bacSiteC = _ancID["CHARS"];
            _bacSiteDim = Columns(_bacSiteC);
            _bacCounter = Rows(_thisColumn) - 1;

            for (_bacTreeIterator = 0; _bacTreeIterator < _bacCounter; _bacTreeIterator = _bacTreeIterator + 1) {
                _bacParentID = ((_ancID["TREE_AVL"])[_bacTreeIterator + 1])["Parent"] - 1;
                _myState = _thisColumn[_bacTreeIterator];
                _pState = _thisColumn[_bacParentID];
                _bacStateLabel = "";


                if (_myState >= 0) {
                    if (mode == 0) {
                        if (_myState != _pState && _pState >= 0) {
                            _bacStateLabel = _bacSiteC[_pState] + "->" + _bacSiteC[_myState];
                        }
                    } else {
                        _bacStateLabel = _bacSiteC[_myState];
                    }
                }

                ((_ancID["TREE_AVL"])[_bacTreeIterator + 1])["SubLabel"] = _bacStateLabel;
            }

            return PostOrderAVL2StringAnnotate(_ancID["TREE_AVL"], _scaled, "SubLabel");
        }
    }
    return "";
}


/*******************************************

	_filterDimensions returns a {{sites,branches}} matrix

*******************************************/

function _filterDimensions(_ancID) {
    _result = {
        1,
        2
    };
    if (Abs_ancID) {
        _sites = Columns(_ancID["MATRIX"]);
        _branches = Rows(_ancID["MATRIX"]) - 1;
        _result[0] = _sites;
        _result[1] = _branches;
    }
    return _result;
}

/*******************************************

	_countSubstitutionsByBranchSite returns a
	binary vector (one entry per branch) which
	contains '1' if the branch has substitutions
	defined by '_filter' (a CxC matrix)

*******************************************/

function _countSubstitutionsByBranchSite(_ancID, _siteID, _filter) {
    _result = {
        0,
        1
    };
    if (Abs_ancID) {
        if (_siteID >= 0 && _siteID < Columns(_ancID["MATRIX"])) {
            TREE_OUTPUT_OPTIONS = {};
            _thisColumn = (_ancID["MATRIX"])[-1][_siteID];
            _bacSiteC = _ancID["CHARS"];
            _bacSiteDim = Columns(_bacSiteC);
            _bacCounter = Rows(_thisColumn) - 1;
            _result = {
                _bacCounter,
                1
            };

            for (_bacTreeIterator = 0; _bacTreeIterator < _bacCounter; _bacTreeIterator += 1) {
                _bacParentID = ((_ancID["TREE_AVL"])[_bacTreeIterator + 1])["Parent"] - 1;
                _myState = _thisColumn[_bacTreeIterator];
                _pState = _thisColumn[_bacParentID];
                _bacSiteMx = {
                    _bacSiteDim,
                    _bacSiteDim
                };
                _expandSubstitutionMap(_pState, _myState, _ancID, "_bacSiteMx");
                _result[_bacTreeIterator] = (+(_bacSiteMx$_filter)) >= 1;
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

function _tabulateSubstitutionsAtSiteByBranch(_ancID, _siteID) {
    _result = {};
    if (Abs_ancID) {
        if (_siteID >= 0 && _siteID < Columns(_ancID["MATRIX"])) {
            _thisColumn = (_ancID["MATRIX"])[-1][_siteID];
            _bacSiteC = _ancID["CHARS"];
            _bacSiteDim = Columns(_bacSiteC);
            _bacCounter = Rows(_thisColumn) - 1;


            for (_bacTreeIterator = 0; _bacTreeIterator < _bacCounter; _bacTreeIterator = _bacTreeIterator + 1) {
                _bacParentID = ((_ancID["TREE_AVL"])[_bacTreeIterator + 1])["Parent"] - 1;
                _myState = _thisColumn[_bacTreeIterator];
                _pState = _thisColumn[_bacParentID];

                if (_pState >= 0 && _myState >= 0) {
                    haveS = _OBSERVED_S_[_pState][_myState];
                    haveNS = _OBSERVED_NS_[_pState][_myState];
                } else {
                    haveS = 0;
                    haveNS = 0;
                }
                _result[((_ancID["TREE_AVL"])[_bacTreeIterator + 1])["Name"]] = {
                    {
                        haveS__,
                        haveNS__
                    }
                };
            }

        }
    }
    return _result;
}


/*******************************************/

lfunction ancestral._expandSubstitutionMap(_state1, _state2, ambig_map, _resultMatrix) {
    if (_state1 >= 0 && _state2 >= 0) {
        /* both unique states */
        _resultMatrix[_state1][_state2] = _resultMatrix[_state1][_state2] + 1;
    } else {
        if (_state1 != (-1) && _state2 != (-1)) {
            /* neither one is a gap state: otherwise nothing to map */
            _thisDim = Rows(_resultMatrix);
            if (_state1 >= 0) {
                _vec1 = {
                    _thisDim,
                    1
                };
                _vec1[_state1] = 1;
            } else {
                _vec1 = ambig_map[-_state1 - 2];
            }
            if (_state2 >= 0) {
                _vec2 = {
                    _thisDim,
                    1
                };
                _vec2[_state2] = 1;
            } else {
                _vec2 = ambig_map[-_state2 - 2];
            }
            _vec1 = _vec1 * Transpose(_vec2);
            _vec2 = +_vec1;
            _resultMatrix += _vec1 * (1 / _vec2); /* 1/K(=all possible pairwise resolutions) to each possible resolution */
        }

    }
}

/*******************************************

*******************************************/

function _convertSubstitutionToCharacters(_state1, _state2, _ancID, _resultMatrix & ) {
    if (_state1 >= 0 && _state2 >= 0)
    /* both unique states */
    {
        _resultMatrix[_state1][_state2] = _resultMatrix[_state1][_state2] + 1;
    } else {
        if (_state1 != (-1) && _state2 != (-1))
        /* at least one is NOT a gap state: otherwise nothing to map */
        {
            _thisDim = Rows(_resultMatrix);
            if (_state1 >= 0) {
                _vec1 = {
                    _thisDim,
                    1
                };
                _vec1[_state1] = 1;
            } else {
                _vec1 = (_ancID["AMBIGS"])[-_state1 - 2];
            }
            if (_state2 >= 0) {
                _vec2 = {
                    _thisDim,
                    1
                };
                _vec2[_state2] = 1;
            } else {
                _vec2 = (_ancID["AMBIGS"])[-_state2 - 2];
            }
            _vec1 = _vec1 * Transpose(_vec2);
            _vec2 = {
                _thisDim,
                1
            }["1"];
            _vec2 = (Transpose(_vec2) * _vec1 * _vec2)[0];
            _resultMatrix = _resultMatrix + _vec1 * (1 / _vec2); /* 1/K(=all possible pairwise resolutions*/
            to each possible resolution * /
        }

    }
    return 0;
}
