/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
   Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
   Art FY Poon    (apoon42@uwo.ca)
   Steven Weaver (sweaver@temple.edu)

Module Developers:
        Lance Hepler (nlhepler@gmail.com)
        Martin Smith (martin.audacis@gmail.com)

Significant contributions from:
  Spencer V Muse (muse@stat.ncsu.edu)
  Simon DW Frost (sdf22@cam.ac.uk)

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#include "dataset_filter.h"
#include "avllistxl_iterator.h"

//_________________________________________________________
// Data Set Filter/Numeric
//_________________________________________________________

_DataSetFilter::_DataSetFilter(void) {
  unitLength = 0;
  theData = NULL;
  accessCache = nil;
}
//_________________________________________________________
_DataSetFilter::_DataSetFilter(_DataSet *ds, char, _String &) {
  theData = ds;
  accessCache = nil;
}
//_________________________________________________________
_DataSetFilter::~_DataSetFilter(void) { DeleteObject(accessCache); }

//_______________________________________________________________________

void _DataSetFilter::CopyFilter (_DataSetFilter const *copyFrom) {
    memcpy ((char*)this, (char*)copyFrom, sizeof (_DataSetFilter));
    
    theFrequencies.Duplicate        (&copyFrom->theFrequencies);
    theNodeMap.Duplicate            (&copyFrom->theNodeMap);
    theMap.Duplicate                (&copyFrom->theMap);
    theOriginalOrder.Duplicate      (&copyFrom->theOriginalOrder);
    conversionCache.Duplicate       (&copyFrom->conversionCache);
    duplicateMap.Duplicate          (&copyFrom->duplicateMap);
    
    dimension               = copyFrom->dimension;
    undimension             = copyFrom->undimension;
    unitLength              = copyFrom->unitLength;
    accessCache             = nil;
    
}

//_______________________________________________________________________

BaseRef _DataSetFilter::makeDynamic (void) const {
    _DataSetFilter * r = new _DataSetFilter;
    r->CopyFilter   (this);
    
    return r;
}


//_______________________________________________________________________
void    _DataSetFilter::SetDimensions (void)
{
    dimension   = GetDimension(true);
    undimension = GetDimension(false);
}

//_______________________________________________________________________
unsigned long    _DataSetFilter::FindUniqueSequences  (_SimpleList& indices, _SimpleList& map, _SimpleList& counts, short mode) const {
    indices.Clear(); map.Clear(); counts.Clear();
    
    unsigned long             sites  = theMap.lLength,
    seqs   = theNodeMap.lLength,
    unit   = GetUnitLength();
    
    if (mode == 0) {
        _SimpleList hashSupport;
        _AVLListXL  sequenceHashes     (&hashSupport);
        
        for (unsigned long sequenceIndex = 0; sequenceIndex < seqs; sequenceIndex ++){
            _String * thisSequence = GetSequenceCharacters (sequenceIndex);
            
            long     sequenceHash   = thisSequence->Adler32(),
            f              = sequenceHashes.Find ((BaseRef)sequenceHash),
            rawSequenceIdx = theNodeMap.lData[sequenceIndex];
            
            DeleteObject (thisSequence);
            
            _SimpleList * sameScore = nil;
            if (f>=0) {
                sameScore = (_SimpleList*)sequenceHashes.GetXtra (f);
                for (long k = 0; k<sameScore->lLength; k++) {
                    bool fit = true;
                    f = sameScore->lData[k];
                    
                    long fRaw = theNodeMap.lData[indices.lData[f]];
                    
                    for (unsigned long site = 0; site < sites && fit; site++){
                        for (unsigned long unitIndex = 0; unitIndex < unit; unitIndex ++){
                            _Site * thisSite = theData->GetSite(theMap.lData[unit*site+unitIndex]);
                            if (thisSite->char_at(fRaw) != thisSite->char_at(rawSequenceIdx)) {
                                fit = false;
                                break;
                            }
                        }
                    }
                    
                    if (fit) {
                        map << f;
                        counts.lData[f] ++;
                        
                    } else {
                        f = -1;
                    }
                }
            }
            if (f==-1) { // fit failed or unique site
                if (!sameScore) {
                    sameScore = new _SimpleList;
                    sequenceHashes.Insert ((BaseRef)sequenceHash,(long)sameScore,false);
                }
                
                (*sameScore) << indices.lLength;
                map     << indices.lLength;
                indices << sequenceIndex;
                counts  << 1;
            }
        }
        
    }
    else{
        long             vd  = GetDimension(true);
        
        hyFloat      *translatedVector = new hyFloat [vd],
        *translatedVector2= new hyFloat [vd];
        
        _String         state1 (unit,false),
        state2 (unit,false);
        
        for (long sequenceIndex = 0; sequenceIndex < seqs; sequenceIndex++) {
            bool checkState = false;
            for (long idx=0; idx<indices.countitems(); idx++) {
                for (long m=0; m<sites; m++) {
                    RetrieveState (m,sequenceIndex, state1,false);
                    RetrieveState (m,indices.lData[idx], state2,false);
                    
                    checkState = true;
                    long idx1 = Translate2Frequencies (state1, translatedVector,  true),
                    idx2 = Translate2Frequencies (state2, translatedVector2, true);
                    
                    //printf ("(%ld, %ld) %ld = %ld %ld\n", sequenceIndex, indices.lData[idx], m, idx1, idx2);
                    
                    if (idx2 >=0 && idx1 >=0) {
                        if (idx1==idx2) {
                            continue;
                        } else {
                            checkState = false;
                            break;
                        }
                    } else {
                        
                        // check for equal ambigs
                        long k = 0;
                        for (; k < vd; k++){
                            if (translatedVector[k] != translatedVector2[k]){
                                break;
                            }
                        }
                        
                        if (k == vd)
                            continue;
                        
                        if (mode == 1){
                            
                            long count1 = 0,
                            count2 = 0;
                            
                            for (long t = 0; t<vd; t++) {
                                count1 += translatedVector[t]>0.0;
                                count2 += translatedVector2[t]>0.0;
                            }
                            
                            if (count1 < vd && count2 < vd) {
                                checkState = false;
                                break;
                            }
                            
                        } else {
                            bool first  = mode==2,
                            second = mode==2;
                            if (mode == 2){
                                for (long t = 0; (first||second)&&(t<vd); t++) {
                                    if (translatedVector[t]>0.0) {
                                        second &= (translatedVector2[t]>0.0);
                                    }
                                    if (translatedVector2[t]>0.0) {
                                        first  &= (translatedVector[t]>0.0);
                                    }
                                }
                                if (!(first||second)) {
                                    checkState = false;
                                    break;
                                }
                            } else {
                                for (long t = 0; t<vd; t++) {
                                    if (translatedVector[t]>0.0) {
                                        second |= (translatedVector2[t]>0.0);
                                    }
                                    if (translatedVector2[t]>0.0) {
                                        first  |= (translatedVector[t]>0.0);
                                    }
                                }
                                if (!(first&&second)) {
                                    checkState = false;
                                    break;
                                }
                            }
                        }
                    }
                }
                
                if (checkState) {
                    map << idx;
                    counts.lData[idx] ++;
                    break;
                }
            }
            
            if (!checkState){
                map     << indices.lLength;
                indices << sequenceIndex;
                counts  << 1;
            }
            
        }
        
        delete [] translatedVector;
        delete [] translatedVector2;
    }
    
    
    return indices.lLength;
}



//_______________________________________________________________________

void    _DataSetFilter::SetFilter (_DataSet const * ds, unsigned char unit, _SimpleList& horizontalList, _SimpleList& verticalList, bool isFilteredAlready) {
    // we must step thru the underlying dataset and recompute the frequencies
    // we will store the vertical map in theMap
    // and the horizontal map in theNodeMap
    // theFrequencies will be store the new frequencies
    // theOriginalOrder is the receptacle for the original site order in the data filter
    
    bool            copiedSelf = false; // tag if refiltering self
    
    _DataSetFilter* firstOne = nil;
    if (isFilteredAlready) {
        if ((hyPointer)this == (hyPointer)ds) {
            firstOne = (_DataSetFilter*)makeDynamic();
            copiedSelf = true;
        } else {
            firstOne = (_DataSetFilter*)ds;
        }
        ds       = firstOne->theData;
    }
    
    theMap.Clear();
    theNodeMap.Clear();
    theOriginalOrder.Clear();
    theFrequencies.Clear();
    theExclusions.Clear();
    conversionCache.Clear();
    duplicateMap.Clear();
    
    theData     = (_DataSet*)ds;
    unitLength  = unit;
    
    long        i,
    j;
    
    // security checks
    if (horizontalList.empty() || verticalList.lLength<unit) {
        ReportWarning (_String("Row and/or column partition is empty. All the data will be used by default."));
        if (horizontalList.empty()) {
            horizontalList.Populate (isFilteredAlready ? firstOne->theNodeMap.lLength : ds->NoOfSpecies(),0,1);
        }
        if (verticalList.lLength<unit) {
            verticalList.Clear();
            verticalList.Populate (isFilteredAlready ? firstOne->GetSiteCount() : ds->GetNoTypes(),0,1);
        }
    }
    
    if (!isFilteredAlready) {
        theNodeMap.Clear();
        theNodeMap.Duplicate (&horizontalList);
    } else {
        for (unsigned long k = 0UL; k<horizontalList.lLength; k++) {
            theNodeMap<<firstOne->theNodeMap.lData[horizontalList.lData[k]];
        }
        
        horizontalList.Clear();
        horizontalList.Duplicate(&verticalList);
        verticalList.Clear();
        verticalList.RequestSpace(firstOne->theOriginalOrder.lLength);
        
        for (i = 0; i<horizontalList.lLength; i++) {
            j = horizontalList.lData[i];
            if (j>=0 && j<firstOne->theOriginalOrder.lLength) {
                verticalList<<firstOne->theOriginalOrder.lData[j];
            } else {
                _String tooBig (j);
                if (j<0) {
                    ReportWarning  (tooBig &" is a negative site index and is ignored");
                } else {
                    ReportWarning  (tooBig &" exceeds the number of sites in the underlying data filter and is ignored");
                }
            }
        }
    }
    
    j = ds->NoOfSpecies();
    
    for (i=0; i<theNodeMap.lLength; i++) {
        if (theNodeMap.lData[i]>=j) {
            _String invalid(theNodeMap.lData[i]);
            ReportWarning ((invalid&" exceeds the number of species in the underlying dataset and is ignored"));
            theNodeMap.Delete(i);
            i--;
        }
    }
    
    j = ds->GetNoTypes();
    for (i=0; i<verticalList.lLength; i++)
        if (verticalList.lData[i]>=j) {
            _String invalid(verticalList.lData[i]);
            ReportWarning ((invalid&" exceeds the number of sites in the underlying dataset and is ignored"));
            verticalList.Delete(i);
            i--;
        }
    
    if (verticalList.lLength%unit!=0) {
        ReportWarning (_String("Number of sites in datasetfilter is not divisible by the unit - will truncate to the nearest integer"));
        while(verticalList.lLength%unit) {
            verticalList.Delete(verticalList.lLength-1);
        }
    }
    
    
    theOriginalOrder.Duplicate (&verticalList);
    
    
    // done with security checks
    
    _SimpleList indices;        // numeric indices intended to facilitate the reindexing
    _AVLListXL  siteIndices     (&indices);
    
    
    // sweep through the columns left to right
    
    duplicateMap.RequestSpace (verticalList.lLength/unit+1);
    
    _String      siteHolder   (unit*theNodeMap.lLength,false);
    
    //bool       startD = false;
    
    for (i=0; i<verticalList.lLength; i+=unit) {
        long colIndex = 0;
        
        for (j=0; j<unit; j++) // sweep within one block
            for (long k=0; k<theNodeMap.lLength; k++) // sweep down the columns
                                                      //colIndex+=
                                                      //(((_String**)ds->lData)[ds->theMap.lData[verticalList.lData[i+j]]])->sData[theNodeMap.lData[k]];
            {
                siteHolder[colIndex++] = (((_String**)ds->lData)[ds->theMap.lData[verticalList.lData[i+j]]])->char_at(theNodeMap.lData[k]);
            }
        
        colIndex = siteHolder.Adler32();
        
        long        f = siteIndices.Find ((BaseRef)colIndex);
        _SimpleList * sameScore = nil;
        
        if (f>=0) {
            sameScore = (_SimpleList*)siteIndices.GetXtra (f);
            for (long k = 0; k<sameScore->lLength; k++) {
                bool fit = true;
                f = sameScore->lData[k];
                for (long j=0; fit&&(j<unit); j++) { // sweep within one block
                    _Site * site1 = ds->GetSite(verticalList.lData[i+j]),
                    * site2 = ds->GetSite(theMap.lData[unit*f+j]);
                    
                    for (long k=0L; k<theNodeMap.lLength; k++) { // sweep down the columns
                        long mapped_index = theNodeMap.lData[k];
                        if (site1->char_at (mapped_index) != site2->char_at (mapped_index)) {
                            fit = false;
                            break;
                        }
                    }
                }
                
                if (fit) {
                    theFrequencies[f]++;
                    duplicateMap<<f;
                    f = 0;
                    break;
                } else {
                    f = -1;
                }
            }
        }
        if (f==-1) { // fit failed or unique site
            if (!sameScore) {
                sameScore = new _SimpleList;
                siteIndices.Insert ((BaseRef)colIndex,(long)sameScore,false);
            }
            
            (*sameScore) << theFrequencies.lLength;
            duplicateMap<<theFrequencies.lLength;
            theFrequencies<<1;
            for (j=0; j<unit; j++) {
                theMap<<verticalList.lData[i+j];
            }
        }
    }
    
    siteIndices.Clear();
    
    duplicateMap.TrimMemory();
    theOriginalOrder.TrimMemory();
    
    if (copiedSelf) {
        DeleteObject (firstOne);
    }
    
    SetDimensions();
    FilterDeletions();
    
}
//_______________________________________________________________________
long    _DataSetFilter::FindSpeciesName (_List& s, _SimpleList& r) const {
  
    // TODO SLKP 20171002: this is a generic list->list map; should be a method in _List
    r.Clear();
    
    _List           newNames;
    _AVLListX       matched (&newNames);
    
    for (unsigned long k=0UL; k<theNodeMap.lLength; k++) {
        long i = theNodeMap.lData[k];
        matched.Insert (new _String (((_String*)theData->theNames (i))->ChangeCase(kStringUpperCase)),i);
    }
    
    for (unsigned long m = 0UL; m < s.lLength; m++) {
        _String const key = ((_String*)s.GetItem (m))->ChangeCase(kStringUpperCase);
        long f = matched.Find (&key);
        if (f>=0L) {
            r << matched.GetXtra (f);
        } else {
            break;
        }
    }
    
    return r.lLength;
}

//_______________________________________________________________________

void    _DataSetFilter::FilterDeletions(_SimpleList *theExc) {
  
    // TODO SLKP 20171002: test this function because there were many semantic changes including the ability to handle both exclusion states and n-fold gaps
  
    bool  skip_nfolds = hy_env::EnvVariableTrue(hy_env::skip_omissions);
  
    if (skip_nfolds || theExc ) { // something to do
        _SimpleList patterns_to_be_removed;
        if (theExc) {
          hyFloat   *store_vec = new hyFloat [GetDimension(false)];
          patterns_to_be_removed = theFrequencies.FilterIndex (
              [&] (long value, unsigned long index) -> bool {
                long invalid_state = HasExclusions(index, theExc, store_vec);
                if (invalid_state != -1) {
                  ReportWarning ((*this)(index,invalid_state).Enquote () & " was encountered in sequence "& *GetSequenceName (invalid_state) & " at site pattern " & (long)index
                                 & ". All corresponding alignment columns will be removed from subsequent analyses.");
                  return true;
                }
                return false;
              }
          );
          delete [] store_vec;
        }
        if (skip_nfolds) {
          theFrequencies.Each (
                       [&] (long value, unsigned long index) -> void {
                         if( HasDeletions(index)) {
                           patterns_to_be_removed.BinaryInsert(index);
                         }
                       }
          );
        }
      
        
        if (patterns_to_be_removed.countitems() == GetPatternCount()) {
            ReportWarning("All the sites in the datafilter have deletions and removing them creates an emptyString filter");
        }
        
        _SimpleList data_sites_to_be_deleted, // e.g. nucleotides
                    filter_sites_to_be_deleted; // e.g. codons
      
      
        if (patterns_to_be_removed.countitems()) {
          /**
           
              Deleting all instances of selected patterns; also need to reindex 
              remaining duplicateMap to reference correct remaining indices
              For example, if the original site->pattern map was like this
           
              0,1,1,2,3,2,4 
           
              and we are deleting patterns 0 and 2, then the updated list should look like
           
              0,0,1,2
           
           */
          
          _SimpleList remapped_duplicates ((unsigned long)(duplicateMap.countitems() - filter_sites_to_be_deleted.countitems())),
                      running_indexer ((unsigned long)theFrequencies.countitems());
          
          long        skip_offset = 0L;
          
          
          duplicateMap.Each ( [&](long value, unsigned long index) -> void {
            long delete_this_entry = patterns_to_be_removed.BinaryFind(value);
            if (delete_this_entry != -1) {
              if (running_indexer.countitems() <= delete_this_entry) { // first time across this site pattern
                running_indexer << -1L;
                skip_offset ++;
              }
              data_sites_to_be_deleted.AppendRange (unitLength, index*unitLength, 1L);
            } else {
              if (running_indexer.countitems() <= delete_this_entry) { // first time across this site pattern
                running_indexer << value - skip_offset;
              }
              remapped_duplicates << running_indexer.get (value);
            }
          });
        }
      
        filter_sites_to_be_deleted.Clear();
        theOriginalOrder.DeleteList (data_sites_to_be_deleted);
        theFrequencies.DeleteList (patterns_to_be_removed);
      
        
        for (unsigned long i=0UL; i<patterns_to_be_removed.countitems(); i++) {
            long pattern_index = patterns_to_be_removed.get(i);
            
            for (unsigned long j=0UL; j<unitLength; j++) {
                theMap.lData[pattern_index*unitLength+j]=-1;
                filter_sites_to_be_deleted << pattern_index*unitLength+j;
            }
        }
        
        
        if (data_sites_to_be_deleted.countitems()) {
            /*allDeleted.Sort();*/
            
            _String     warnMsg ("The following sites are being omitted ");
          
            if (!theExc) {
                warnMsg = warnMsg & " (because of n-fold [unresolved] sites or indels)";
            }
            
            ReportWarning(warnMsg& _String ((_String*)data_sites_to_be_deleted.toStr()));
        }
        
        // this seems to be an old debugging code snippet
        // TODO SLKP 20171002 : review this
        _SimpleList saveMap (theMap);
        theMap.DeleteList (filter_sites_to_be_deleted);
        for (long k=0; k<theMap.lLength; k++) {
           if (theMap.get (k) < 0) {
              HandleApplicationError (_String ("Internal Error in ") & __PRETTY_FUNCTION__);
          }
        }
    }
}
//_______________________________________________________________________
_DataSetFilter*  _DataSetFilter::PairFilter (long index1, long index2, _DataSetFilter* result) {
    _SimpleList species;
    species<<theNodeMap(index1);
    species<<theNodeMap(index2);
    result->SetFilter (theData,unitLength,species,theMap);
    if (theExclusions.countitems()) {
        _String exclusions ((_String*)theExclusions.toStr());
        exclusions.StripQuotes('{','}');
        result->SetExclusions   (exclusions);
    }
    return result;
}

//_________________________________________________________

void    _DataSetFilter::MatchStartNEnd (_SimpleList& order, _SimpleList& positions, _SimpleList* parent) const {
    // SLKP 20171002: needs to be documented and reviewed
    // for example, the return representation will break if one has more than 2^16 species
  
    if (order.empty() == 0) {
        return;
    }
  
    if (hy_env::EnvVariableTrue(hy_env::use_traversal_heuristic)) {
        if (parent) {
            for (long i = 1; i < order.lLength; i++) {
                unsigned long
                j       = 0,
                n       = theNodeMap.lLength-1,
                p0  = parent->lData[i],
                p1  = order.lData[i];
                
                while (CompareTwoSites(p0,p1,j)) {
                    j++;
                }
                while (CompareTwoSites(p0,p1,n)) {
                    n--;
                }
                n = (n<<16) + j;
                positions << n;
            }
        }
        else {
          long p0 = order.get(0);
          for (long i = 1; i < order.lLength; i++) {
                unsigned long j = 0,
                n = theNodeMap.lLength-1,
                p1 = order.lData[i];
                
                while (CompareTwoSites(p0,p1,j)) {
                    j++;
                }
                while (CompareTwoSites(p0,p1,n)) {
                    n--;
                }
                n = (n<<16) + j;
                positions << n;
                p0 = p1;
          }
        }
    } else {
        for (long i = 1; i < order.lLength; i++) {
            unsigned long j = 0,
            n = theNodeMap.lLength-1;
            
            n = (n<<16) + j;
            positions << n;
        }
    }
  
}

//_______________________________________________________________________

void    _DataSetFilter::SetExclusions (_String const& exclusion_string, bool filter) {
  
    theExclusions.Clear();
    _String character_list = exclusion_string;
    character_list.StripQuotes();
    if (character_list.empty()) {
        return;
    }
    
    _AVLList     exclusions (&theExclusions);
  
    character_list.Tokenize(',').ForEach ([&] (BaseRefConst  exclusion_character) -> void {
      _String* kth_token = (_String*)exclusion_character;
      long character_index = MapStringToCharIndex(*kth_token);
      if (character_index < 0) {
        ReportWarning (_String("Exclusion request for '") & *kth_token &"' does not represent a unique state and will therefore be ignored.");
      } else {
        if (exclusions.InsertNumber(character_index) < 0) {
          ReportWarning (_String("Exclusion symbol for '") & *kth_token &"' is included more than once.");
        }
      }
    });

    exclusions.ReorderList();
    if (filter) {
        FilterDeletions (&theExclusions);
    }
    
 }

//_______________________________________________________________________

_String*    _DataSetFilter::GetExclusions (void) const {
    _StringBuffer * res = new _StringBuffer (16UL);
    
    if (theExclusions.empty() == false) {
        (*res) << ConvertCodeToLetters (theExclusions.get(0L), unitLength);
        for (long k=1; k<theExclusions.lLength; k++) {
            (*res) << ',' << ConvertCodeToLetters (theExclusions.get(k), unitLength);
        }
    }
    
    res->TrimSpace();
    return res;
}

//_______________________________________________________________________

unsigned long    _DataSetFilter::GetDimension (bool correct) const {
    unsigned long result = ComputePower (theData->theTT->LengthOfAlphabet(), unitLength);
    return correct ? result - theExclusions.lLength : result;
}

//_______________________________________________________________________

void    _DataSet::ProcessPartition (_String const & input2 , _SimpleList & target , bool isVertical, _SimpleList const* additionalFilter, _SimpleList const* otherDimension, _String const* scope) const {
    // TODO SLKP : 20170928 this needs serious cleanup and testing
    
    if (input2.empty()) {
        return;
    }
    // decide if the input is an enumeration or a formula
    long totalLength;
    
    if (additionalFilter) {
        totalLength = additionalFilter->lLength;
    } else {
        totalLength = isVertical?theMap.lLength:noOfSpecies;
    }
    
    _String input (input2);
    
    if (!input.IsALiteralArgument(true)) { // not a literal argument
    
        _Formula fmla, lhs;
        _FormulaParsingContext fpc;
        fpc.setScope (scope);
        
        long     outcome = Parse (&fmla, input, fpc,&lhs);
        
        if (outcome!=HY_FORMULA_EXPRESSION) {
            HandleApplicationError(input.Enquote() & _String(" is an invalid partition specification"));
            return;
        }
        _PMathObj   fV = fmla.Compute();
        if (fV && fV->ObjectClass()==STRING) {
            ProcessPartition (((_FString*)fV)->get_str().Enquote(), target, isVertical, additionalFilter, nil, scope);
        } else {
            _DataSet::MatchIndices (fmla, target, isVertical, totalLength, scope);
        }
    } else { // an explicit enumeration or a regular expression
        if (input (0) =='/' && input (-1) == '/') {
            // a regular expression
            input.Trim(1,input.length()-2);
            int   errCode;
            regex_t*   regex = _String::PrepRegExp (&input, errCode, true);
            if (errCode) {
                HandleApplicationError(_String::GetRegExpError(errCode));
                return;
            }
            // now set do the matching
            // using only the sites that are specced in the additionalFilter
            
            if (!isVertical) {
                
                for (long specCount = 0L; specCount < (additionalFilter ? additionalFilter->lLength : totalLength); specCount++) {
                    _String pattern ((unsigned long)theMap.countitems());
                    long    seqPos = additionalFilter ? additionalFilter->Element (specCount) : specCount;
                    
                    if (otherDimension)
                        for (long seqSlider = 0L; seqSlider < otherDimension->lLength; seqSlider ++) {
                            pattern.set_char(seqSlider, GetSite(otherDimension->Element(seqSlider))->get_char (seqPos));
                        }
                    else
                        for (long seqSlider = 0L; seqSlider < theMap.lLength; seqSlider ++) {
                            pattern.set_char(seqSlider, GetSite(seqSlider)->get_char (seqPos));
                        }
                    
                    if (pattern.RegExpMatch (regex, 0L).countitems()) {
                        target << specCount;
                    }
                }
            } else {
                bool         *eligibleMarks = new bool[lLength] {false};
                
                if (additionalFilter) {
                    for (long siteIndex = 0; siteIndex < additionalFilter->lLength; siteIndex ++) {
                        eligibleMarks[theMap.lData[additionalFilter->lData[siteIndex]]] = true;
                    }
                }
                else {
                    InitializeArray(eligibleMarks, lLength, true);
                }
                
                _SimpleList matches;
                _String     *tempString = nil;
                if (otherDimension) {
                    tempString = new _String (otherDimension->countitems());
                }
                
                for (long siteCounter = 0; siteCounter < lLength; siteCounter ++)
                    if (eligibleMarks[siteCounter]) {
                        matches.Clear();
                        if (otherDimension) {
                            _Site * aSite = ((_Site**)lData)[siteCounter];
                            for (long tc = 0; tc < otherDimension->lLength; tc++) {
                                tempString->set_char (tc, aSite->char_at(otherDimension->get(tc)));
                            }
                            matches = tempString->RegExpMatch (regex, 0L);
                        } else {
                            matches = ((_Site**)lData)[siteCounter]->RegExpMatch (regex, 0L);
                        }
                        if (matches.empty()) {
                            eligibleMarks[siteCounter] = false;
                        }
                    }
                
                DeleteObject (tempString);
                
                if (additionalFilter) {
                    for (long afi = 0; afi < additionalFilter->lLength; afi++)
                        if (eligibleMarks[theMap.lData[additionalFilter->lData[afi]]]) {
                            target << afi;
                        }
                } else {
                    for (long afi = 0; afi < theMap.lLength; afi++)
                        if (eligibleMarks[theMap.lData[afi]]) {
                            target << afi;
                        }
                }
                delete eligibleMarks;
            }
            _String::FlushRegExp (regex);
        } else {
            input = input.KillSpaces ();
            // now process the string
            long count = 0L,anchor;
            
            _SimpleList numbers,
            links;
            
            numbers.RequestSpace (1024);
            links.RequestSpace (1024);
            
            // first check if it is has a comb filter
            
            if ( input (0) =='<' && input (-1) =='>') {
                for (count=1; count<input.length()-1; count++) {
                    if (input.char_at(count) != '0') {
                        numbers<<count-1;
                    }
                }
                if (numbers.countitems()) {
                    long k = input.length()-2; // step size
                    anchor = 0;
                    if (totalLength == -1) {
                        totalLength = theMap.lLength;
                    }
                    while (anchor<totalLength-k) {
                        for (count = 0; count< numbers.lLength; count++) {
                            target<<anchor+numbers.lData[count];
                        }
                        anchor+=k;
                    }
                    if ( (k=totalLength-1-anchor) ) {
                        for (count = 0; count< numbers.lLength; count++) {
                            if (numbers.lData[count]>k) {
                                break;
                            }
                            target<<anchor+numbers.lData[count];
                        }
                    }
                    return;
                    
                }
            }
            
            while (count<input.length()) {
                anchor = count;
                
                for (; count<input.length() && isdigit(input.char_at (count)); count++) ;
                
                long    aNumber = (input.Cut (anchor,count-1)).to_long();
                
                if (aNumber < 0) {
                    _String warnMsg ("A negative number was found in partition specification: ");
                    ReportWarning (warnMsg & input.Cut (0,anchor-1) & '?' & input.Cut (anchor,-1));
                    target.Clear();
                    return;
                }
                numbers<< aNumber;
                
                char current_char = input.char_at (count);
                
                if (current_char == '<' || current_char == '>') {
                    ReportWarning (_String  ("A comb partition cannot be combined with other types. The entire partition is reset to first..last") & input.Cut (0,anchor-1) & '?' & input.Cut (anchor,-1));
                    target.Clear();
                    return;
                }
                
                if (current_char == '&') {
                    links << numbers.lLength;
                }
                
                // TODO SLKP 20171001 this needs to be checked for correctness
                if (current_char == ',' || count == input.length()) { // wrap it up dude
                    if (numbers.countitems() == 1) {
                        target<<numbers(0);
                    } else {
                        if (links.empty()) {
                            if (numbers[0]>numbers[1]) { // backward order
                                for (long k = numbers[0]; k>=numbers[1]; k--) {
                                    target<<k;
                                }
                            } else {
                                for (long k = numbers[0]; k<=numbers[1]; k++) {
                                    target<<k;
                                }
                            }
                        } else {
                            // linked locations
                            if (links.countitems() != (numbers.countitems()-2) / 2) {
                                ReportWarning ("A part of the partition specification has not been understood and is being skipped.");
                                target.Clear();
                                return;
                            } else {
                                _SimpleList signs;
                                signs<<(numbers(0)<numbers(1)?1:-1);
                                for (long k = 0; k<links.lLength; k+=2) {
                                    signs<<(numbers(links(k))<numbers(links(k+1))?1:-1);
                                }
                                
                                for (long k=numbers(0), l=0 ; signs(0)*k<=signs(0)*numbers(1); k+=signs(0), l++) {
                                    target<<numbers(0)+l*signs(0);
                                    for (long m=0; m<links.lLength; m++) {
                                        target<<numbers(links(m))+l*signs(m+1);
                                    }
                                }
                             }
                        }
                    }
                    numbers.Clear();
                    links.Clear();
                }
                count++;
            }
        }
    }
}
//_______________________________________________________________________

void     _DataSetFilter::SetMap  (_String const &s) {
    theNodeMap.Clear();
    if (s.nonempty()) {
        s.Tokenize(_String (",")).ForEach([&] (BaseRef piece) -> void {
            theNodeMap << ((_String*)piece)->to_long();
        });
    }
}

//_______________________________________________________________________
_String* _DataSetFilter::MakeSiteBuffer (void) const {
    return new _String ((unsigned long)unitLength, false);
}

//_______________________________________________________________________
void _DataSetFilter::retrieve_individual_site_from_raw_coordinates (_String& store, unsigned long site, unsigned long sequence) const {
  if (unitLength==1UL) {
    store.set_char (0, (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site]]])->char_at (sequence));
  } else {
    site*=unitLength;
    for (unsigned long k = 0UL; k<unitLength; k++) {
      store.set_char (k, ((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site++]]]->char_at(sequence));
    }
  }
}

//_______________________________________________________________________

_String&     _DataSetFilter::operator () (unsigned long site, unsigned long pos) {
    // TODO: 20171002 SLKP, this check needs to happen when unitLength changes (once!)
    if (!accessCache || accessCache->length() != unitLength) {
        if (accessCache) {
            DeleteObject (accessCache);
        }
        accessCache = new _String ((unsigned long)unitLength);
    }
  
    retrieve_individual_site_from_raw_coordinates (*accessCache, site, theNodeMap.get (pos));
    return *accessCache;
}

//_______________________________________________________________________

const _String     _DataSetFilter::RetrieveState (unsigned long site, unsigned long pos) const {
    _String state ((unsigned long)unitLength, false);
    RetrieveState (site, pos, state, false);
    return state;
}

//_______________________________________________________________________

void     _DataSetFilter::RetrieveState (unsigned long site, unsigned long pos, _String& reply, bool map) const {
    retrieve_individual_site_from_raw_coordinates (reply, map ? duplicateMap.get(site) : site, theNodeMap.get (pos));
}

//_______________________________________________________________________

_List *  _DataSetFilter::ComputePatternToSiteMap (void) const {
    _List * result = new _List ();
    
    for (unsigned long k = 0UL; k < theFrequencies.countitems(); k++) {
        (*result) < new _SimpleList;
    }
    for (unsigned long s = 0UL; s < duplicateMap.lLength; s++) {
        *((_SimpleList*)(result->GetItem(duplicateMap.get(s)))) << s;
    }
    return result;
}


//_______________________________________________________________________

char _DataSetFilter::direct_index_character (unsigned long site, unsigned long sequence) const {
  return (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site]]])->char_at(sequence);
}

  //_______________________________________________________________________

bool     _DataSetFilter::CompareTwoSites (unsigned long site1, unsigned long site2, unsigned long pos1) const {
  pos1 = theNodeMap.lData[pos1];
  
  
  switch (unitLength) {
    case 3: { // codon
      site1*=3;
      site2*=3;
      return direct_index_character (site1, pos1)     == direct_index_character (site2,     pos1) &&
      direct_index_character (site1 + 1, pos1) == direct_index_character (site2 + 1, pos1) &&
      direct_index_character (site1 + 2, pos1) == direct_index_character (site2 + 2, pos1);
    }
      
    case 1: { // nucs, protein, etc
      return direct_index_character (site1, pos1)     == direct_index_character (site2,     pos1);
    }
      
    case 2: { // di-nucs
      site1*=2;
      site2*=2;
      return direct_index_character (site1, pos1)     == direct_index_character (site2,     pos1) &&
      direct_index_character (site1 + 1, pos1) == direct_index_character (site2 + 1, pos1);
    }
      
      
    default: {
      site1*=unitLength;
      site2*=unitLength;
      
      
      for (unsigned long k = 0UL; k<unitLength; k++) {
        if (direct_index_character (site1 + k, pos1)     != direct_index_character (site2 + k,     pos1)) {
          return false;
        }
      }
      return true;
    }
      
  }
}

//_______________________________________________________________________
bool    _DataSetFilter::HasDeletions (unsigned long site, _AVLList* storage) const {
    // TODO 20171002: in this and other functions, replace local caches with
    // an instance of
  
    long        filter_dimension  = GetDimension(true),
                sequence_count    = theNodeMap.countitems()?theNodeMap.countitems():theData->NoOfSpecies();
  
    hyFloat* store    = new hyFloat [filter_dimension];
    _String buffer ((unsigned long) GetUnitLength());
  
    bool outcome = false;
    
    for (unsigned long k = 0UL; k<sequence_count; k++) {
      
        RetrieveState(site, k, buffer, false);
        Translate2Frequencies (buffer, store, false);
        
        bool has_ones = false,
             has_zeros = false;
        
        for (unsigned long j = 0UL; j<filter_dimension; j++) {
            if (store[j]==0.0) {
                has_zeros = true;
                if (has_ones) {
                  break;
                }
            } else if (store[j]==1.0) {
                has_ones = true;
                if (has_zeros) {
                  break;
                }
            }
        }
        if (!(has_ones && has_zeros )) {
            if (storage) {
                outcome = true;
                storage->InsertNumber(theNodeMap.get(k));
            } else {
                delete [] store;
                return true;
            }
        }
    }
    
    delete [] store;
    return outcome;
}

//_______________________________________________________________________

_Matrix*        _DataSetFilter::GetFilterCharacters (bool flip) const {
    unsigned long  unit_length = GetUnitLength (),
    seq_length  = flip?theFrequencies.lLength:(GetSiteCount () / unitLength),
    seq_count   = NumberSpecies();
    
    _List       result;
    
    _StringBuffer      char_buffer (unit_length);
    
    if (flip) {
        for (long k=0; k< seq_length; k++) {
            _StringBuffer *alignment_column = new _StringBuffer ((unsigned long)seq_count+1);
            for (long k2=0; k2< seq_count ; k2++) {
                RetrieveState(k,k2,char_buffer,false);
                (*alignment_column) << char_buffer;
            }
            alignment_column->TrimSpace();
            result  < alignment_column;
        }
    } else {
        for (long k=0; k < seq_count; k++) {
            result      < GetSequenceCharacters(k);
        }
    }
    
    return new _Matrix (result, false);
}

//_______________________________________________________________________

_String*        _DataSetFilter::GetSequenceCharacters (long seqID)  const{
    long            unitSizeL   = GetUnitLength();
    
    _StringBuffer * aSequence = new _StringBuffer (GetSiteCount());
  
  
    
    if (seqID >= 0 && seqID < theNodeMap.countitems()) {
        _StringBuffer      aState (unitSizeL);
        unsigned long        upTo = GetSiteCountInUnits();
        for (unsigned long k2=0UL; k2<upTo; k2++) {
            RetrieveState(k2,seqID,aState);
            (*aSequence) << aState;
        }
    }
    aSequence->TrimSpace ();
    return aSequence;
}


//_______________________________________________________________________
long    _DataSetFilter::HasExclusions (unsigned long site, _SimpleList* theExc, hyFloat*store ) const {
  
    if (theNodeMap.countitems()) {
        _String buffer ((unsigned long)GetUnitLength());
  
        for (unsigned long k = 0UL; k<theNodeMap.countitems(); k++) {
            RetrieveState           (site, k, buffer, false);
            Translate2Frequencies   (buffer, store, false);
            
          
            unsigned long   filter_dim = GetDimension(false);
            long            found_forbidden = -1;
          
            for (unsigned long character = 0UL ;  character < filter_dim;  character ++) {
                   if (store[character] > 0.0) {
                      if (theExc->Find(character) < 0) {
                          return -1; // found at least one non-excluded character (possibly partial)
                      } else {
                        found_forbidden = character;
                      }
                  }
            }
          
            return found_forbidden;
        }
    }
  
    return -1;
}

//_________________________________________________________
_Matrix* _DataSetFilter::ComputePairwiseDifferences (long i, long j, _hy_dataset_filter_ambiguity_resolution resolution_option) const {
    // TODO: 20171002, needs a more in-depth code review
    
    try {
        
        if (unitLength > 3) {
            throw _String("ComputePairwiseDifferences is not implemented for data filters with unit size > 3");
        }
        
        long    mxDim      = GetDimension (true);
        
        _Matrix     *res   = new _Matrix  (mxDim,mxDim,false,true);
        
        hyFloat  *sm1   = new hyFloat[mxDim],
        *sm2   = new hyFloat[mxDim];
        
        
        
        _String      state1 ((unsigned long)unitLength),
                     state2 ((unsigned long)unitLength);
        
        
        if (conversionCache.lLength == 0) {
            throw _String ("ComputePairwiseDifferences called on a filter with emptyString conversionCache");
        }
        
        long        *tcodes  = conversionCache.lData+89,
        *ccodes  = conversionCache.lData+1,
        ccount   = conversionCache.lData[0];
        
        for (unsigned long site_pattern = 0UL; site_pattern < theFrequencies.countitems(); site_pattern++) {
            long s1 = -1, s2 = -1;
          
            long seq_i = theNodeMap.get (i),
                 seq_j = theNodeMap.get (j);
          
            int c1, c2;
            
            c1 = direct_index_character (site_pattern, seq_i),
            c2 = direct_index_character (site_pattern, seq_j);
            
            if (unitLength == 1) {
                s1 = conversionCache.lData[(c1-40)*(undimension+1)+undimension],
                s2 = conversionCache.lData[(c2-40)*(undimension+1)+undimension];
            } else {
                int         c12 = direct_index_character (site_pattern + 1, seq_i),
                            c22 = direct_index_character (site_pattern + 1, seq_j);
                
                
                state1.set_char(0, c1);
                state1.set_char(1, c12);

                state2.set_char(0, c2);
                state2.set_char(1, c22);

              
                c1  = ccodes[c1-40];
                c12 = ccodes[c12-40];
                
                c2  = ccodes[c2-40];
                c22 = ccodes[c22-40];
                
                if (unitLength == 2) {
                    if (c1>=0 && c12>=0) {
                        s1 = tcodes[c1*ccount+c12];
                    }
                    
                    if (c2>=0 && c22>=0) {
                        s2 = tcodes[c2*ccount+c22];
                    }
                } else {
                    int         c13 = direct_index_character (site_pattern + 2, seq_i),
                    c23 = direct_index_character (site_pattern + 2, seq_j);
                    
                    //printf ("\n%c %c", c13, c23);
                    state1.set_char(2, c13);
                    state1.set_char(2, c23);
                  
                    c13 = ccodes[c13-40];
                    c23 = ccodes[c23-40];
                    
                    //printf (" %d %d %s %s\n", c13, c23, state1.sData, state2.sData);
                    
                    if (c1>=0 && c12>=0 && c13>=0) {
                        s1 = tcodes[ccount*(c1*ccount+c12)+c13];
                    }
                    
                    if (c2>=0 && c22>=0 && c23>=0) {
                        s2 = tcodes[ccount*(c2*ccount+c22)+c23];
                    }
                }
            }
          
            
            if (s1>=0 && s2>=0) { // one to one
                res->theData[s1*mxDim+s2] += site_frequency(site_pattern);
            } else {
                if (resolution_option != kAmbiguityHandlingSkip) {
                    _Matrix * freqsAtSite = nil;
                    
                    if (resolution_option != kAmbiguityHandlingResolve) {
                        _SimpleList   //seqList,
                        siteList;
                        
                        
                        for (long si = 0; si < unitLength; si++) {
                            siteList << theMap.lData[unitLength*site_pattern+si];
                        }
                        
                        _SimpleList copy_node_order (theNodeMap);
                        freqsAtSite     = theData->HarvestFrequencies (unitLength, unitLength, 0, copy_node_order, siteList);
                        if (theExclusions.lLength) {
                            long k = 0,
                            u = GetDimension (false);
                            
                            for (long i = 0; i<u; i++) {
                                if (i==theExclusions.lData[k] && k<theExclusions.lLength) {
                                    k++;
                                    continue;
                                }
                                freqsAtSite->theData[i-k] = freqsAtSite->theData[i];
                            }
                        }
                        //XferwCorrection (freqsAtSite->theData, freqsAtSite->theData, mxDim);
                    }
                    
                    if (s1>=0) {
                        // one to many
                        if (unitLength>1) {
                            Translate2Frequencies (state2,sm1,false);
                        } else {
                            Translate2Frequencies (c2,sm1,false);
                        }
                        
                        if (freqsAtSite) {
                            if (resolution_option == kAmbiguityHandlingAverageFrequencyAware) {
                                hyFloat totalW = 0.0;
                                
                                for  (long m=0; m<mxDim; m++)
                                    if (sm1[m]>0.0) {
                                        totalW += freqsAtSite->theData[m];
                                    }
                                
                                if (totalW>0.0) {
                                    s1 = s1*mxDim;
                                    
                                    for  (long m=0; m<mxDim; m++,s1++)
                                        if (sm1[m]>0.0) {
                                            res->theData[s1] += theFrequencies.lData[site_pattern]*freqsAtSite->theData[m]/totalW;
                                        }
                                }
                                
                            } else {
                                hyFloat maxW   = 0.0;
                                long       maxIdx = -1;
                                
                                for  (long m=0; m<mxDim; m++) {
                                    if (sm1[m]>0.0) {
                                        hyFloat myWeight = freqsAtSite->theData[m];
                                        if (myWeight > maxW) {
                                            maxW = myWeight;
                                            maxIdx = m;
                                        }
                                    }
                                }
                                
                                if (maxIdx>=0) {
                                    res->theData[s1*mxDim+maxIdx] += theFrequencies.lData[site_pattern];
                                }
                            }
                        } else {
                            /* adopt the following convention here:
                             - if ambig resolves to one s1 - count as a match
                             - otherwise - count all contributions equally
                             */
                            
                            if (sm1[s1] > 0.0) {
                                res->theData[s1*mxDim+s1] += theFrequencies.lData[site_pattern];
                            } else {
                                long ambCount = 0;
                                for  (long m=0; m<mxDim; m++) {
                                    if (sm1[m]>0.0) {
                                        ambCount ++;
                                    }
                                }
                                
                                s1 *= mxDim;
                                
                                hyFloat addFac = theFrequencies.lData[site_pattern]/(hyFloat)ambCount;
                                
                                for  (long m=0; m<mxDim; m++,s1++)
                                    if (sm1[m]>0.0) {
                                        res->theData[s1] += addFac;
                                    }
                            }
                        }
                    } else {
                        if (s2>=0)
                            // many to one
                        {
                            if (unitLength>1) {
                                Translate2Frequencies (state1,sm1,false);
                            } else {
                                Translate2Frequencies (c1,sm1,false);
                            }
                            
                            if (freqsAtSite) {
                                if (resolution_option == kAmbiguityHandlingAverageFrequencyAware) {
                                    hyFloat totalW = 0.0;
                                    
                                    for  (long m=0; m<mxDim; m++)
                                        if (sm1[m]>0.0) {
                                            totalW += freqsAtSite->theData[m];
                                        }
                                    
                                    if (totalW>0.0) {
                                        for  (long m=0; m<mxDim; m++,s2+=mxDim)
                                            if (sm1[m]>0.0) {
                                                res->theData[s2] += theFrequencies.lData[site_pattern]*freqsAtSite->theData[m]/totalW;
                                            }
                                    }
                                    
                                } else {
                                    hyFloat maxW   = 0.0;
                                    long       maxIdx = -1;
                                    
                                    for  (long m=0; m<mxDim; m++) {
                                        if (sm1[m]>0.0) {
                                            hyFloat myWeight = freqsAtSite->theData[m];
                                            if (myWeight > maxW) {
                                                maxW = myWeight;
                                                maxIdx = m;
                                            }
                                        }
                                    }
                                    
                                    if (maxIdx>=0) {
                                        res->theData[maxIdx*mxDim+s2] += theFrequencies.lData[site_pattern];
                                    }
                                }
                            } else {
                                if (sm1[s2] > 0.0) {
                                    res->theData[s2*mxDim+s2] += theFrequencies.lData[site_pattern];
                                } else {
                                    long ambCount = 0;
                                    for  (long m=0; m<mxDim; m++)
                                        if (sm1[m]>0.0) {
                                            ambCount ++;
                                        }
                                    
                                    hyFloat addFac = theFrequencies.lData[site_pattern]/(hyFloat)ambCount;
                                    {
                                        for  (long m=0; m<mxDim; m++,s2+=mxDim)
                                            if (sm1[m]>0.0) {
                                                res->theData[s2] += addFac;
                                            }
                                    }
                                }
                            }
                        } else
                            // many to many
                        {
                            if (unitLength>1) {
                                Translate2Frequencies (state1,sm1,false);
                                Translate2Frequencies (state2,sm2,false);
                            } else {
                                Translate2Frequencies (c1,sm1,false);
                                Translate2Frequencies (c2,sm2,false);
                            }
                            
                            if (freqsAtSite) {
                                if (resolution_option == kAmbiguityHandlingAverageFrequencyAware) {
                                    hyFloat totalW = 0.0;
                                    
                                    for  (long m=0; m<mxDim; m++)
                                        if (sm1[m]>0)
                                            for  (long m2=0; m2<mxDim; m2++)
                                                if (sm2[m2]>0) {
                                                    totalW += freqsAtSite->theData[m]*freqsAtSite->theData[m2];
                                                }
                                    
                                    if (totalW>0.0) {
                                        for  (long m=0; m<mxDim; m++)
                                            if (sm1[m]>0)
                                                for  (long m2=0; m2<mxDim; m2++)
                                                    if (sm2[m2]>0) {
                                                        res->theData[m*mxDim+m2] += theFrequencies.lData[site_pattern]*freqsAtSite->theData[m]*freqsAtSite->theData[m2]/totalW;
                                                    }
                                    }
                                    
                                } else {
                                    hyFloat maxW   = 0.0;
                                    long       maxIdx  = -1,
                                    maxIdx2 = -1;
                                    
                                    for  (long m=0; m<mxDim; m++)
                                        if (sm1[m]>0)
                                            for  (long m2=0; m2<mxDim; m2++)
                                                if (sm2[m2]>0) {
                                                    hyFloat myWeight = freqsAtSite->theData[m]*freqsAtSite->theData[m2];
                                                    if (myWeight > maxW) {
                                                        maxW = myWeight;
                                                        maxIdx  = m;
                                                        maxIdx2 = m2;
                                                    }
                                                }
                                    
                                    if (maxIdx>=0) {
                                        res->theData[maxIdx*mxDim+maxIdx2] += theFrequencies.lData[site_pattern];
                                    }
                                }
                            } else {
                                long ambCount  = 0,
                                ambCount2 = 0,
                                m         = 0;
                                
                                for  (; m<mxDim; m++) {
                                    if (sm1[m]>0.0) {
                                        if (sm2[m]>0.0) {
                                            break;
                                        } else {
                                            ambCount ++;
                                        }
                                    } else if (sm2[m]>0.0) {
                                        ambCount2 ++;
                                    }
                                }
                                
                                if (m==mxDim) {
                                    hyFloat addFac = theFrequencies.lData[site_pattern]/(hyFloat)(ambCount*ambCount2);
                                    
                                    for  (long m=0; m<mxDim; m++)
                                        if (sm1[m]>0)
                                            for  (long m2=0; m2<mxDim; m2++)
                                                if (sm2[m2]>0) {
                                                    res->theData[m*mxDim+m2] += addFac;
                                                }
                                }
                            }
                        }
                    }
                    DeleteObject (freqsAtSite);
                }
            }
        }
        
        delete[] sm1;
        delete[] sm2;
        
        return res;
    }
    catch (const _String error) {
        HandleApplicationError(error);
        return    new _Matrix (1,1,false,true);
        
    }
}

//_________________________________________________________

void _DataSetFilter::ComputePairwiseDifferences (_Matrix& target, long i, long j) const
// matrix of dimension nx4n containing pairwise distances as follows (n=number of species)
// first lower diag - count the same (AA,CC,GG,TT)
// first upper diag - count AC,CA
// 2nd   lower diag - count AG,GA
// 2nd   upper diag - count AT,TA
// 3rd   lower diag - count CG,GC
// 3rd   upper diag - count CT,TC
// 4th   lower diag - count GT,TG
{
    if ((target.GetHDim()!=1)||(target.GetVDim()!=7)) {
        _Matrix::CreateMatrix (&target,1,7,false,true,false);
    }
    
    if (!theData->theTT->IsStandardNucleotide()) {
        return;
    }
    long k,l;
    
    for (k=0; k<7; k++) {
        target[k] = 0;
    }
    
    k = theNodeMap.lData[i];
    l = theNodeMap.lData[j];
    
    if (l>k) {
        EXCHANGE (k,l);
    }
    
    for (unsigned long m=0; m < theMap.lLength; m++) {
        char const * thisSite = GetColumn (m);
        char a = thisSite[k],
        b = thisSite[l];
        
        long fc = theFrequencies.lData[m/unitLength];
        
        if (a>b) {
            EXCHANGE (a,b);
        }
        
        if (a==b) {
            target[0]+=fc;
        } else {
            if (a=='A') {
                switch (b) {
                    case 'C': {
                        target[1]+=fc;
                        break;
                    }
                    case 'G': {
                        target[2]+=fc;
                        break;
                    }
                    case 'T': {
                        target[3]+=fc;
                        break;
                    }
                }
            } else if (a=='C') {
                switch (b) {
                    case 'G': {
                        target[4]+=fc;
                        break;
                    }
                    case 'T': {
                        target[5]+=fc;
                        break;
                    }
                }
            } else if (a=='G') {
                if (b=='T') {
                    target[6]+=fc;
                }
            }
            
        }
    }
}

//_________________________________________________________

_Matrix * _DataSetFilter::HarvestFrequencies (char unit, char atom, bool posSpec, bool countGaps) const {
    _SimpleList copy_seqs (theNodeMap), copy_sites (theOriginalOrder);
    return theData->HarvestFrequencies (unit,atom, posSpec, copy_seqs, copy_sites, countGaps);
}


/*
 //_______________________________________________________________________
 long    _DataSetFilter::GetVectorCode(long site,long seq)
 {
 if (!symbolVector) return -1;
 long* fi = symbolVector->quickArrayAccess();
 return fi[*fi*site+seq+1];
 }
 //_______________________________________________________________________
 
 void    _DataSetFilter::ProduceSymbolVector(bool smear)
 {
 // compute the size of the vector cells
 hyFloat cellSize=log((hyFloat)theData->theTT->LengthOfAlphabet())*hyFloat(unitLength)/log(128.0);
 if (cellSize>2.0)
 {
 _String errMsg ("DataSetFilter has more than 32767 states, which is currently unsupported");
 FlagError(errMsg);
 }
 long intCellSize = cellSize>1.0?2:1;
 // now produce the conversion vector
 long sites = theMap.lLength, species= theNodeMap.lLength?theNodeMap.lLength:theData->NoOfSpecies();
 symbolVector = new _SimpleList ();
 checkPointer(symbolVector);
 //  (*symbolVector)<<intCellSize;
 (*symbolVector)<<species;
 // we will now speciate into byte and word size cases
 // the data will be stored column by column
 // if there is a unique code translation, we then store that code in the symbol vector for faster
 // processing during tree pruning business.
 // use a standard convert to frequencies function to check whether a character has a unique conversion
 if (intCellSize==1) // char based storage
 {
 union
 {
 long composite;
 char bytes[sizeof(long)];
 } converterb;
 char byteposition = 0, bytesPerLong = sizeof(long);
 for (long i=0;i<sites;i++)
 {
 for (long j=0; j<species; j++)
 {
 //              if (byteposition==bytesPerLong)
 //              {
 //                  byteposition = 0;
 //                  (*symbolVector)<<converterb.composite;
 //              }
 //              converterb.bytes[byteposition]=(char)Translate2Frequencies((*this)(i,j),nil,smear,false);
 //              byteposition++;
 (*symbolVector)<<Translate2Frequencies((*this)(i,j),nil,smear,false);
 }
 }
 if (byteposition)
 {
 for(long i=sizeof(long)-1;i>=byteposition;i--)
 {
 converterb.bytes[i]=0;
 }
 (*symbolVector)<<converterb.composite;
 }
 
 }
 else
 {
 union
 {
 long composite;
 short int words[sizeof(long)/2];
 } converterw;
 char wordposition = 0, wordsPerLong = sizeof(long);
 for (long i=0;i<sites;i++)
 {
 for (long j=0; j<species; j++)
 {
 if (wordposition==wordsPerLong)
 {
 wordposition = 0;
 (*symbolVector)<<converterw.composite;
 }
 converterw.words[wordposition]=(char)Translate2Frequencies((*this)(i,j),nil,smear,false);
 }
 }
 if (wordposition)
 {
 for(long i=sizeof(long)/2-1;i>=wordposition;i--)
 {
 converterw.words[i]=0;
 }
 (*symbolVector)<<converterw.composite;
 }
 }
 }*/

//_______________________________________________________________________

long    _DataSetFilter::CorrectCode (long code) const {
    if (theExclusions.lLength == 0) {
        return code;
    }
    return theExclusions.SkipCorrect (code);
}



//_______________________________________________________________________
long    _DataSetFilter::Translate2Frequencies (_String const& str, hyFloat* parvect, bool smear) const {
    // TODO : SLKP 20171002, this will break ugly if dimension > HYPHY_SITE_DEFAULT_BUFFER_SIZE
  
    long  store      [HYPHY_SITE_DEFAULT_BUFFER_SIZE],
                      resolution_count  = -1L;
    
    
    InitializeArray(parvect, dimension, 0.);
    
    if (unitLength == 1) {
        resolution_count = theData->theTT->TokenResolutions (str.char_at(0),store,smear);
    } else {
        resolution_count = theData->theTT->MultiTokenResolutions(str,store, smear);
    }
    
    long mapped_resolution_count = theExclusions.lLength ? theExclusions.CorrectForExclusions(store, resolution_count) : resolution_count;
    
    /* handle the cases when no unambiguous resolutions were available */
    for (long i = 0L; i < mapped_resolution_count; i++) {
        parvect[store[i]] = 1.;
    }
    
    if (mapped_resolution_count == 1L) {
        return store[0];
    }
    if (mapped_resolution_count == 0L && resolution_count == 0L && smear) {
        InitializeArray(parvect, dimension, 1.);
    }
    
    return -1L;
}


//_______________________________________________________________________

long    _DataSetFilter::MapStringToCharIndex (_String& str) const {
    
    long  store      [HYPHY_SITE_DEFAULT_BUFFER_SIZE],
    resolution_count  = -1L;
    
    
    if (unitLength == 1) {
        resolution_count = theData->theTT->TokenResolutions (str.char_at(0),store);
    } else {
        resolution_count = theData->theTT->MultiTokenResolutions(str,store);
    }
    
    long mapped_resolution_count = theExclusions.lLength ? theExclusions.CorrectForExclusions(store, resolution_count) : resolution_count;
    
    if (mapped_resolution_count == 1L) {
        return store[0];
    }
    
    return -1L;
}


//_______________________________________________________________________
long    _DataSetFilter::Translate2Frequencies (char s, hyFloat* parvect, bool smear) const {
    long  store      [HYPHY_SITE_DEFAULT_BUFFER_SIZE],
    resolution_count  = theData->theTT->TokenResolutions (s,store,smear);
    
    long mapped_resolution_count = theExclusions.lLength ? theExclusions.CorrectForExclusions(store, resolution_count) : resolution_count;
    
    if (mapped_resolution_count == 0L) {
        if (smear) {
            InitializeArray(parvect, dimension, 1.);
            return -1;
        }
    }
    
    InitializeArray (parvect, dimension, 0.);
    
    for (long i = 0L; i < mapped_resolution_count; i++) {
        parvect[store[i]] = 1.;
    }
    
    return resolution_count==1L?1L:-1L;
}

//_______________________________________________________________________
long    _DataSetFilter::LookupConversion (char s, hyFloat* parvect) const
{
    if (undimension==4) {
        long* cCache = conversionCache.lData+(s-40)*5;
        parvect[0] = cCache[0];
        parvect[1] = cCache[1];
        parvect[2] = cCache[2];
        parvect[3] = cCache[3];
        return cCache[4];
        
    } else {
        int idx = (s-40)*(undimension+1);
        for (long i=0; i<undimension; parvect[i++] = conversionCache.lData[idx++]) ;
        return conversionCache.lData[idx];
    }
}
//_______________________________________________________________________
bool   _DataSetFilter::ConfirmConversionCache() const {
    return conversionCache.lLength || unitLength > 3;
}

//_______________________________________________________________________
void    _DataSetFilter::SetupConversion (void) {
    if (conversionCache.countitems()) {
        return;
    }
    
    if ( unitLength==1 ) {
        char c = 40;
        hyFloat *temp    = new hyFloat [undimension+1UL];
        
        while(c<127) {
            //InitializeArray(temp, undimension + 1UL, 0.0);
            
            Translate2Frequencies(c, temp, true);
            
            long resolution_count = -1;
            
            
            for (unsigned long i=0UL; i<undimension; i++) {
                long character_code_resolution =  (long)temp[i];
                conversionCache << character_code_resolution;
                if (character_code_resolution) {
                    if (resolution_count == -1) {
                        resolution_count = i;
                    } else {
                        resolution_count = -2;
                    }
                }
            }
            
            conversionCache<<resolution_count;
            c++;
        }
        delete[] temp;
    } else {
        
        if (unitLength==2 || unitLength==3) {
            
            _String alphabet = theData->theTT->GetAlphabetString();
            unsigned long alphabet_dim = alphabet.length();
            
            
            long  ccache [88],
            uncorrected_dimension = GetDimension(false) ;
            
            conversionCache.RequestSpace (89+uncorrected_dimension);
            conversionCache << alphabet.length();
            
            for (unsigned long i=0UL; i<88; i++) {
                ccache[i] = -1;
            }
            for (unsigned long i=0UL; i<alphabet.length(); i++) {
                ccache [alphabet.char_at(i)-40] = i;
            }
            for (unsigned long i=0UL; i<88; i++) {
                conversionCache << ccache[i];
            }
            
            _String s ((unsigned long)unitLength);
            for (unsigned long char_index = 0; char_index < uncorrected_dimension; char_index++ ) {
                _SimpleList components = SplitIntoDigits (char_index, unitLength, alphabet_dim);
                for (unsigned long position = 0; position < unitLength; position ++) {
                    s.set_char(position, alphabet.char_at(components.get(position)));
                }
                conversionCache << MapStringToCharIndex(s);
            }
        }
    }
}

  //_________________________________________________________

_String const _DataSetFilter::GenerateConsensusString (_SimpleList* majority) const {
  
  if (unitLength > 3) {
    return kEmptyString;
  }
  
  _String     result ((unsigned long)theOriginalOrder.lLength),
  pattern_consensus  ((unsigned long)(unitLength*theFrequencies.lLength));
  
  long        char_states         = GetDimension(false),
  *translation_buffer = new long [char_states];
  
  hyFloat* count_buffer = new hyFloat [char_states];
  
  for (unsigned long site_pattern = 0UL; site_pattern<theFrequencies.lLength; site_pattern ++) {
    long    index_in_dataset = theMap.lData[site_pattern];
    
    InitializeArray (count_buffer, char_states, 0.);
    
    for (unsigned long sequence_index =0UL; sequence_index < theNodeMap.lLength; sequence_index ++) {
      long resolution_count = theData->theTT->TokenResolutions ((*theData)(index_in_dataset, theNodeMap.lData[sequence_index],1),translation_buffer, false);
      
      
      if (resolution_count>1L) {
        hyFloat equal_weight = 1./resolution_count;
        for (long resolution_index = 0L; resolution_index < resolution_count; resolution_index++) {
          count_buffer [translation_buffer[resolution_index]] += equal_weight;
        }
      } else {
        if (resolution_count == 1) {
          count_buffer [translation_buffer[0]] += 1.;
        }
      }
    }
    
      // find the residue with the highest frequency
    
    hyFloat       max_weight      = -1.;
    InitializeArray (translation_buffer, char_states, 0L);
    long             max_char_count  = 0L;
    
    for (unsigned long char_index = 0UL; char_index < char_states; char_index++) {
      if (StoreIfGreater(max_weight, count_buffer[char_index])) {
        max_char_count = 1;
        translation_buffer [0] = char_index;
      } else {
        if (count_buffer[char_index] == max_weight) {
          translation_buffer [max_char_count ++] = char_index;
        }
      }
    }
    
    if (max_char_count > 1L) {
      pattern_consensus.set_char(site_pattern, theData->theTT->AmbigToLetter(translation_buffer, max_char_count));
    } else {
      pattern_consensus.set_char(site_pattern, theData->theTT->ConvertCodeToLetters(translation_buffer[0],1) [0]);
    }
    if (majority) {
      (*majority) << max_weight;
    }
  }
  
  delete [] count_buffer;
  delete [] translation_buffer;
  
  for (unsigned long m=0UL; m<theOriginalOrder.lLength; m++) {
    result.set_char (m, pattern_consensus.char_at(duplicateMap.get(m)));
  }
  
  return result;
}


  //_________________________________________________________
void    _DataSetFilter::toFileStr (FILE*dest, unsigned long) {
    // write out the file with this dataset filter
  if (dest) {
      internalToStr (dest,nil);
  }
}

  //_________________________________________________________
void    _DataSetFilter::ConvertCodeToLettersBuffered (long code, unsigned char unit, _String& storage, _AVLListXL* lookup) const {
    // write out the file with this dataset filter
  long            lookupC     = lookup->FindLong (code);
  const char      *lookupV;
  if (lookupC>=0) {
    lookupV = ((_String*)lookup->GetXtra(lookupC))->get_str();
  } else {
    _String * newT = new _String (ConvertCodeToLetters (code,unit));
    lookup->Insert ((BaseRef)code, (long)newT, false);
    lookupV = newT->get_str();
  }
  
  if (unit == 1) {
    storage.set_char (0, lookupV[0]);
  } else {
    for (unsigned long k = 0UL; k < unit; k++) {
      storage.set_char (k,lookupV[k]);
    }
  }
}




  //_________________________________________________________

void    _DataSetFilter::internalToStr (FILE * file ,_StringBuffer * string_buffer) {
  
  
  auto trim_to_10 = [] (const _String& seq_name) -> _String const& {
    if (seq_name.length() >= 10) {
      return seq_name.Cut (0,9) & ' ';
    }
    return seq_name & _String (_String (" "), 11-seq_name.length ());
  };
  
  long outputFormat = hy_env::EnvVariableGetNumber(hy_env::data_file_print_format),
       printWidth   = hy_env::EnvVariableGetNumber(hy_env::data_file_default_width),
       gapWidth     = hy_env::EnvVariableGetNumber(hy_env::data_file_gap_width);
  
    // write out the file with this dataset filter
  
  unsigned long sequence_count = NumberSpecies(),
  site_count     = GetSiteCount();
  
  if (printWidth <= 0) {
    printWidth = hy_env::EnvVariableGetDefaultNumber(hy_env::data_file_default_width);
  }
  if (gapWidth <= 0) {
    gapWidth = hy_env::EnvVariableGetDefaultNumber(hy_env::data_file_gap_width);
  }
  
  StringFileWrapper write_here (file ? nil : string_buffer, file);
  
  if (outputFormat < 4 || outputFormat > 8) {
      // not NEXUS or serial
    if (!(theData->theTT->IsStandardNucleotide() || theData->theTT->IsStandardAA())) {
      _String * bSet = &theData->theTT->baseSet;
      
      write_here << "$BASESET:\""
      << *bSet
      << "\"\n";
      
      if (theData->theTT->tokensAdded.nonempty()) {
        for (long at = 0; at < theData->theTT->tokensAdded.length(); at++) {
          write_here << "$TOKEN:\""
          << theData->theTT->tokensAdded.char_at(at)
          << "\" = \""
          << theData->theTT->ExpandToken (theData->theTT->tokensAdded.char_at(at))
          << "\"\n";
        }
      }
    }
  }
  
  switch (outputFormat) {
    case 1: // hash-mark interleaved
    case 10: { // FASTA interleaved
      
      long sitesDone    = 0,
      upTo;
      
      char seqDelimiter = (outputFormat==1)?'#':'>';
      
      for (unsigned long i = 0UL; i<theNodeMap.lLength; i++) {
        write_here << seqDelimiter
        << GetSequenceName(i)
        << kStringFileWrapperNewLine;
      }
      
      while (sitesDone<theOriginalOrder.lLength) {
        
        write_here << kStringFileWrapperNewLine
        << kStringFileWrapperNewLine;
        
        
        upTo = sitesDone+printWidth;
        if (upTo>theOriginalOrder.lLength) {
          upTo = theOriginalOrder.lLength;
        }
        
        for (unsigned long i = 0UL; i<theNodeMap.lLength; i++) {
          for (unsigned long j = sitesDone; j<upTo; j++) {
            if ((j-sitesDone)%gapWidth==0) {
              write_here << ' ';
            }
            write_here << (*theData)(theOriginalOrder.lData[j],theNodeMap.lData[i],1);
          }
          
          write_here << kStringFileWrapperNewLine;
        }
        
        sitesDone = upTo;
      }
      break;
    }
      
    case 2:     // PHYLIP sequential
    case 11:    // PAML
    {
      
      write_here << _String((long)theNodeMap.lLength)
      << kStringFileWrapperTab
      << _String(theOriginalOrder.lLength)
      << kStringFileWrapperNewLine;
      
        // proceed to spool out the data
      for (unsigned long i = 0UL; i<theNodeMap.lLength; i++) {
        _String const * sequence_name = GetSequenceName(i);
        _String sequence_name_10;
        
        if (outputFormat == 2) { // PHYLIP
          sequence_name_10 = trim_to_10 (*sequence_name);
        } else {
          sequence_name_10 = *sequence_name & "  ";
        }
        
        write_here << sequence_name_10;
        
        for (unsigned long site_index = 0; site_index<theOriginalOrder.lLength; site_index++) {
          if ((site_index%printWidth==0)&&site_index) {
            write_here << "\n           ";
          }
          
          write_here << (*theData)(theOriginalOrder(site_index),theNodeMap(i),1);
          if (site_index%gapWidth==gapWidth-1) {
            write_here << ' ';
          }
        }
        write_here << kStringFileWrapperNewLine;
        
      }
      break;
    }
      
    case 3: { // phylip interleaved
              // print PHYLIP format header
              //fprintf (dest,"$FORMAT:\"PHYLIPI\"\n");
              // print number of species and sites
      
      write_here << _String((long)theNodeMap.lLength)
      << kStringFileWrapperTab
      << _String(theOriginalOrder.lLength)
      << kStringFileWrapperNewLine;
      
        // proceed to spool out the data
      for (unsigned long i = 0UL; i<theNodeMap.lLength; i++) {
        write_here << trim_to_10 (*GetSequenceName(i));
        
        for (unsigned long j = 0UL; j<theOriginalOrder.lLength; j++) {
          if (j==printWidth) {
            write_here << kStringFileWrapperNewLine;
          } else {
            if (j%gapWidth==0) {
              write_here << ' ';
            }
            write_here << (*theData)(theOriginalOrder.lData[j],theNodeMap.lData[i],1);
          }
        }
      }
      
      
      
      unsigned long completed = printWidth;
      
      while (completed<theOriginalOrder.lLength-1) {
        long upTo = completed+printWidth<theOriginalOrder.lLength?completed+printWidth:theOriginalOrder.lLength;
        for (unsigned long i = 0UL; i<theNodeMap.lLength; i++) {
          write_here << "\n           ";
          for (unsigned long j = completed; j<upTo; j++) {
            if ((j-completed)%gapWidth==0) {
              write_here <<  ' ';
            }
            write_here << (*theData)(theOriginalOrder.lData[j],theNodeMap.lData[i],1);
          }
        }
        completed+=printWidth;
        write_here << kStringFileWrapperNewLine;
      }
      
      
      break;
    }
      
        // various flavors of NEXUS
      
    case 4: // labels, sequential
    case 5: // labels, interleaved
    case 6: // no labels, sequential
    case 7: { // no labels, interleaved
              // write out the header
      
      write_here << "#NEXUS\n\nBEGIN TAXA;\n\tDIMENSIONS NTAX = "
      << _String ((long)sequence_count)
      << ";\n\tTAXLABELS\n\t\t";
      
      for (unsigned long i=0UL; i< sequence_count; i++) {
        write_here << GetSequenceName(i)->Enquote('\'') << ' ';
      }
      
      write_here << ";\nEND;\n\nBEGIN CHARACTERS;\n\tDIMENSIONS NCHAR = "
      << _String((long)theOriginalOrder.lLength)
      << ";\n\tFORMAT\n\t\t";
      
      if (theData->theTT->IsStandardNucleotide()) {
        write_here << "DATATYPE = DNA\n";
      } else {
        if (theData->theTT->IsStandardAA()) {
          write_here << "DATATYPE = PROTEIN\n";
        } else if (theData->theTT->IsStandardBinary()) {
          write_here << "DATATYPE = BINARY\n";
        } else {
          long alphabet_length = theData->theTT->baseSet.length();
          
          write_here << "\t\tSYMBOLS = \"";
          for (unsigned long bc = 0UL; bc < alphabet_length-1; bc++) {
            write_here << theData->theTT->baseSet.char_at (bc)
            << ' ';
          }
          write_here << theData->theTT->baseSet.char_at (alphabet_length-1)
          << "\"\n";
          
          if (theData->theTT->tokensAdded.nonempty())
            for (long at = 0; at < theData->theTT->tokensAdded.length(); at++) {
              write_here << "\nEQUATE =\""
              << theData->theTT->tokensAdded.char_at(at)
              << " = "
              << theData->theTT->ExpandToken(theData->theTT->tokensAdded.char_at(at))
              << "\"";
            }
        }
      }
      if (theData->theTT->GetGapChar()) {
        write_here << "\t\tGAP=" << theData->theTT->GetGapChar();
      }
      if (theData->theTT->GetSkipChar()) {
        write_here << "\n\t\tMISSING=" << theData->theTT->GetSkipChar();
      }
      if (outputFormat>5) {
        write_here << "\n\t\tNOLABELS";
      }
      if (outputFormat%2) {
        write_here << "\n\t\tINTERLEAVE";
      }
      
      write_here << "\n\t;\n\nMATRIX";
      
      
      
        //compute space alignment for different taxa names
        // two passes - one to locate the max length and 2nd to compute padding lengths
      
      unsigned long max_length = 0UL;
      
      for (unsigned long i=0UL; i<sequence_count; i++) {
        StoreIfGreater (max_length, GetSequenceName(i)->length());
      }
      
      _SimpleList taxaNamesPadding;
      
      for (unsigned long i=0UL; i<sequence_count; i++) {
        taxaNamesPadding <<  max_length - GetSequenceName(i)->length();
      }
      
      
      if (outputFormat%2==0) { // sequential
        for (unsigned long i=0UL; i< sequence_count; i++) {
          if (outputFormat == 4) { // labels
            write_here << "\n\t'"
            << GetSequenceName(i)
            << '\''
            << _String (" ", taxaNamesPadding (i));
            
            
          } else {
            write_here << kStringFileWrapperNewLine;
          }
          write_here << ' ';
          for (long site_index = 0UL; site_index < site_count; site_index++) {
            write_here << (*theData)(theOriginalOrder.lData[site_index],theNodeMap.lData[i],1);
          }
        }
      } else {
        long  sitesDone = 0, upTo;
        
        while (sitesDone< site_count) {
          upTo = sitesDone+printWidth;
          
          if (upTo>site_count) {
            upTo = site_count;
          }
          
          
          for (unsigned long i=0UL; i< sequence_count; i++) {
            if (outputFormat == 5) { // labels
              write_here << "\n\t'"
              << GetSequenceName(i)
              << '"'
              << _String (" ", taxaNamesPadding (i));
            } else {
              write_here << kStringFileWrapperNewLine;
            }
            
            write_here << ' ';
            for (long site_index = sitesDone; site_index < upTo; site_index++) {
              write_here << (*theData)(theOriginalOrder.lData[site_index],theNodeMap.lData[i],1);
            }
            
          }
          write_here << kStringFileWrapperNewLine << kStringFileWrapperNewLine;
          sitesDone = upTo;
        }
        
      }
      write_here << ";\nEND;";
      break;
    }
      
    case 8: {
      for (unsigned long i = 0UL; i< sequence_count; i++) {
        write_here << (*theData)(theOriginalOrder(0),theNodeMap(i),1);
        for (unsigned long j = 1UL; j<site_count; j++) {
          write_here << ',' << (*theData)(theOriginalOrder(j),theNodeMap(i),1);
        }
        write_here << kStringFileWrapperNewLine;
      }
      break;
    }
      
    default: { // hash-mark sequential
      char seqDelimiter = (outputFormat==9)?'>':'#';
      
      for (unsigned long i = 0UL; i< sequence_count; i++) {
        write_here << seqDelimiter << GetSequenceName(i);
        for (unsigned long j = 0UL; j<site_count; j++) {
          if (j % printWidth == 0) {
            write_here << kStringFileWrapperNewLine;
          }
          write_here << (*theData)(theOriginalOrder(j),theNodeMap(i),1);
        }
        write_here << kStringFileWrapperNewLine;
      }
      
        // finally see if we need to write out a tree
      
    }
  }
  
  if (outputFormat != 8) {
      if (hy_env::EnvVariableTrue(hy_env::data_file_tree)) {
          _PMathObj tree_var = hy_env::EnvVariableGet(hy_env::data_file_tree_string, HY_ANY_OBJECT);
          if (tree_var) {
            _String* treeString = (_String*)(tree_var->Compute())->toStr();
            switch (outputFormat) {
              case 0:
              case 1:
              case 9:
              case 10: {
                write_here << kStringFileWrapperNewLine
                << kStringFileWrapperNewLine
                << *treeString;
                break;
              }
              case 2:
              case 3: {
                write_here << "\n1\n" << *treeString;
                break;
              }
              default: {
                write_here << "\n\nBEGIN TREES;\n\tTREE tree = "
                << *treeString
                << ";\nEND;";
              }
            }
            DeleteObject (treeString);
          }
        }
  }
}



  //_________________________________________________________

_Matrix * _DataSet::HarvestFrequencies (unsigned char unit, unsigned char atom, bool posSpec, _SimpleList& hSegmentation, _SimpleList& vSegmentation, bool countGaps) const {
  
  if (hSegmentation,empty () || vSegmentation.countitems() < unit) { // revert to default (all data)
    if (hSegmentation.empty ()) {
      hSegmentation.Populate (NoOfSpecies(),0,1);
    }
    if (vSegmentation.countitems () <unit) {
      vSegmentation.Clear();
      vSegmentation.Populate (GetNoTypes(),0,1);
    }
  }
  
  if (unit%atom > 0) { // 20120814 SLKP: changed this behavior to throw errors
    HandleApplicationError (_String("Atom should divide unit in ") & _String (__PRETTY_FUNCTION__).Enquote() &" call");
    return new _Matrix (1,1);
  }
  
  _Matrix   *  out = new _Matrix (ComputePower (theTT->LengthOfAlphabet(), atom),
                                  posSpec?unit/atom:1,
                                  false,
                                  true);
  
  long     positions  =   unit/atom,
  static_store [HYPHY_SITE_DEFAULT_BUFFER_SIZE];
  
  _String unit_for_counting ((unsigned long)atom);
  
  for (unsigned long site_pattern = 0UL; site_pattern <vSegmentation.lLength;  site_pattern +=unit) { // loop over the set of segments
                                                                                                      // make sure the partition is kosher
    
    if (site_pattern + unit > vSegmentation.lLength) {
      break;
    }
    
    for (unsigned long primary_site = site_pattern; primary_site < site_pattern+unit; primary_site += atom) {
      
      long   index_in_pattern = (primary_site-site_pattern)/atom;
      
      for (unsigned long sequence_index = 0; sequence_index <hSegmentation.lLength; sequence_index ++) {
          // loop down each column
        
        unsigned long mapped_sequence_index = hSegmentation.lData[sequence_index];
          // build atomic probabilities
        
        for (unsigned long m = 0UL; m<atom; m++ ) {
          unit_for_counting.set_char (m, (*this)(vSegmentation.lData[primary_site+m],mapped_sequence_index,atom));
        }
        
        long resolution_count = theTT->MultiTokenResolutions(unit_for_counting, static_store, countGaps);
        
        if (resolution_count > 0UL) {
          
          hyFloat    normalized = 1./resolution_count;
          
          for (long resolution_index = 0UL; resolution_index < resolution_count; resolution_index ++) {
            out->theData[posSpec? static_store[resolution_index]*positions+index_in_pattern: static_store[resolution_index]] += normalized;
          }
        }
      }
    }
  }
  
    //scale the matrix now
  
  unsigned long row_count    = out->GetHDim(),
  column_count = out->GetVDim();
  
  for (unsigned long column =0UL; column < column_count; column++) { // normalize each _column_ to sum to 1.
    hyFloat sum = 0.0;
    
    for (unsigned long row = 0UL; row < row_count; row++) {
      sum += out->theData [row*column_count + column];
    }
    
    for (unsigned long row = 0UL; row < row_count; row++) {
      out->theData [row*column_count + column] /= sum;
    }
  }
  
  
  return out;
}

