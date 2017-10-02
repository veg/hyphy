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
    // we must step thru the underlying dataset and recompute the frequenices
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
        ReportWarning (_String("Row and/or column partition is emptyString. All the data will be used by default."));
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
    // MOD 12/16/03
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
extern _String skipOmissions;

void    _DataSetFilter::FilterDeletions(_SimpleList *theExc)
{
    hyFloat      skipo;
    checkParameter (skipOmissions,skipo,0.0);
    
    if (skipo>.5 || theExc ) { // delete omissions
                               //build up the list of "bad" sites
        _SimpleList sitesWithDeletions;
        if (!theExc) {
            for (long i=0; i<theFrequencies.lLength; i++)
                if (HasDeletions(i)) {
                    sitesWithDeletions<<i;
                }
        } else {
            hyFloat   *store_vec = (hyFloat*)checkPointer(new hyFloat [GetDimension(false)]);
            
            for (long i=0; i<theFrequencies.lLength; i++) {
                long pos = HasExclusions(i,theExc,store_vec);
                if  (pos != -1) {
                    sitesWithDeletions<<i;
                    _String warnMsg ((*this)(i,pos));
                    warnMsg = warnMsg & " was encountered in sequence "& *GetSequenceName (pos) & " at site pattern " & i
                    & ". All corresponding alignment columns will be removed from subsequent analyses.";
                    ReportWarning (warnMsg);
                }
            }
            
            delete [] store_vec;
        }
        
        if (sitesWithDeletions.lLength==theFrequencies.lLength) {
            _String errMsg ("All the sites in the datafilter have deletions and removing them creates an emptyString filter");
            ReportWarning(errMsg);
        }
        
        _SimpleList allDeleted,
        dupDeletes;
        
        for (long k=0; k < duplicateMap.lLength; k++)
            if (sitesWithDeletions.BinaryFind (duplicateMap.lData[k]) >= 0) {
                dupDeletes << k;
                for (long j = 0; j < unitLength; j++ ) {
                    allDeleted << k*unitLength + j;
                }
            }
        
        duplicateMap.DeleteList (dupDeletes);
        dupDeletes.Clear();
        theOriginalOrder.DeleteList (allDeleted);
        theFrequencies.DeleteList (sitesWithDeletions);
        
        
        for (long i=0; i<sitesWithDeletions.lLength; i++) {
            long sitePos = sitesWithDeletions.lData[i];
            
            for (long j=0; j<unitLength; j++) {
                theMap.lData[sitePos*unitLength+j]=-1;
                dupDeletes << sitePos*unitLength+j;
            }
        }
        
        
        if (allDeleted.lLength) {
            /*allDeleted.Sort();*/
            
            _String     warnMsg ("The following sites are being omitted:"),
            *s = (_String*)allDeleted.toStr();
            
            if (!theExc) {
                warnMsg = warnMsg & "(b/c of deletions/omissions)";
            }
            
            warnMsg = warnMsg&*s;
            DeleteObject(s);
            ReportWarning(warnMsg);
            
            _SimpleList shiftIdxBy (sitesWithDeletions.lLength+theFrequencies.lLength);
            
            long        shiftBy = sitesWithDeletions.lLength,
            marker  = sitesWithDeletions.lData[sitesWithDeletions.lLength-1],
            markerI = sitesWithDeletions.lLength-2;
            
            shiftIdxBy.lLength = sitesWithDeletions.lLength+theFrequencies.lLength;
            
            for (long i=shiftIdxBy.lLength-1; i>=0; i--) {
                if (i==marker) {
                    shiftBy--;
                    if (markerI>=0) {
                        marker = sitesWithDeletions.lData[markerI];
                        markerI --;
                    } else {
                        marker = -1;
                    }
                }
                shiftIdxBy.lData[i] = shiftBy;
            }
            {
                for (long i=0; i<duplicateMap.lLength; i++) {
                    duplicateMap.lData[i] -= shiftIdxBy.lData[duplicateMap.lData[i]];
                }
            }
        }
        
        // one final pass on theMap to clear it out
        /*for (long i=theMap.lLength-1;i>=0;i--)
         if (theMap(i)<0)
         theMap.Delete(i);*/
        _SimpleList saveMap (theMap);
        theMap.DeleteList (dupDeletes);
        {
            for (long k=0; k<theMap.lLength; k++)
                if (theMap.lData[k] < 0) {
                    saveMap.DeleteList (dupDeletes);
                    WarnError ("Internal Error in _DataSetFilter::FilterDeletions");
                }
        }
    }
    
}
//_______________________________________________________________________
_DataSetFilter*  _DataSetFilter::PairFilter (long index1, long index2, _DataSetFilter* result)
{
    _SimpleList species;
    species<<theNodeMap(index1);
    species<<theNodeMap(index2);
    result->SetFilter (theData,unitLength,species,theMap);
    if (theExclusions.lLength) {
        _String* s = (_String*)theExclusions.toStr();
        *s = s->Cut (1,s->Length()-2);
        result->SetExclusions   (s);
        DeleteObject(s);
    }
    return result;
}

//_________________________________________________________

void    _DataSetFilter::MatchStartNEnd (_SimpleList& order, _SimpleList& positions, _SimpleList* parent) const {
    if (order.lLength == 0) {
        return;
    }
    
    long p0 = order.lData[0];
    
    hyFloat uth;
    checkParameter (useTraversalHeuristic,uth,1.0);
    
    if (uth>.5) {
        if (parent)
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
        else
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
    } else
        for (long i = 1; i < order.lLength; i++) {
            unsigned long j = 0,
            n = theNodeMap.lLength-1;
            
            n = (n<<16) + j;
            positions << n;
        }
    
}

//_______________________________________________________________________

void    _DataSetFilter::SetExclusions (_String* theList, bool filter)
{
    
    theExclusions.Clear();
    theList->StripQuotes();
    
    if (theList->sLength == 0) {
        return;
    }
    
    _List        tokens (theList->Tokenize(','));
    _SimpleList  holder;
    _AVLList     exclusions (&holder);
    
    for (long k = 0; k < tokens.lLength; k++) {
        
        _String* kth_token = (_String*)tokens.GetItem(k);
        
        long posMarker = MapStringToCharIndex(*kth_token);
        
        if (posMarker < 0) {
            ReportWarning (_String("Exclusion request for '") & *kth_token &"' does not represent a unique state and will therefore be ignored.");
        } else {
            if (exclusions.Insert((BaseRef)posMarker) < 0) {
                ReportWarning (_String("Exclusion symbol for '") & *kth_token &"' is included more than once.");
            }
        }
    }
    
    exclusions.ReorderList();
    
    if (filter) {
        FilterDeletions (&holder);
    }
    
    theExclusions<<holder;
}

//_______________________________________________________________________

_String*    _DataSetFilter::GetExclusions (void) const {
    _StringBuffer * res = new _StringBuffer (16UL);
    
    if (theExclusions.empty() == false) {
        for (long k=0; k<theExclusions.lLength-1; k++) {
            (*res) << ConvertCodeToLetters (theExclusions.get(k), unitLength) << ',';
        }
        (*res) << ConvertCodeToLetters (theExclusions.get(theExclusions.lLength-1), unitLength);
    }
    
    res->TrimSpace();
    
    return res;
}

//_______________________________________________________________________

unsigned long    _DataSetFilter::GetDimension (bool correct) const {
    unsigned long result = ComputePower (theData->theTT->baseLength, unitLength);
    return correct ? result - theExclusions.lLength : result;
}

//_______________________________________________________________________

void    _DataSet::ProcessPartition (_String const & input2 , _SimpleList & target , bool isVertical, _SimpleList const* additionalFilter, _SimpleList const* otherDimension, _String const* scope) const {
    // TODO 20170928 this needs serious cleanup and testing
    
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

//_________________________________________________________

void    _DataSetFilter::FindAllSitesLikeThisOne (long index, _SimpleList& receptacle) const {
    long   oindex = theOriginalOrder.Find(index);
    
    if (oindex == kNotFound) {
        return;
    }
    
    if (theData->NoOfSpecies()==theNodeMap.countitems()) { // filter has all sequences
        long *matchMap = new long[unitLength];
        
        for (unsigned long m=0UL; m<unitLength; m++) {
            matchMap[m] = theData->theMap.lData[theOriginalOrder.lData[oindex+m]];
        }
        
        
        for (long k=0; k<theOriginalOrder.lLength; k+=unitLength) {
            unsigned long m = 0UL;
            for (; m<unitLength; m++) {
                if (theData->theMap.lData[theOriginalOrder.lData[k+m]]!=matchMap[m]) {
                    break;
                }
            }
            if (m==unitLength) {
                for (unsigned long m=0UL; m<unitLength; m++) {
                    receptacle<<theOriginalOrder.lData[k+m];
                }
            }
        }
        
        delete [] matchMap;
    } else {
        char ** matchMap = (char**)MemAllocate (sizeof (char*) * unitLength);
        checkPointer (matchMap);
        
        for (m=0; m<unitLength; m++) {
            matchMap[m] = ((_Site*)(((BaseRef*)theData->lData)[theData->theMap.lData[oindex+m]]))->sData;
        }
        for (long k=0; k<theOriginalOrder.lLength; k+=unitLength) {
            for (m=0; m<unitLength; m++) {
                char* checkStr = ((_Site*)(((BaseRef*)theData->lData)[theData->theMap.lData[k+m]]))->sData;
                long t;
                for (t = 0; t<theNodeMap.lLength; t++) {
                    if (checkStr[t]!=matchMap[m][t]) {
                        break;
                    }
                }
                if (t<theNodeMap.lLength) {
                    break;
                }
            }
            if (m==unitLength)
                for (m=0; m<unitLength; m++) {
                    receptacle<<theOriginalOrder.lData[k+m];
                }
        }
        delete matchMap;
    }
}

//_______________________________________________________________________
_String* _DataSetFilter::MakeSiteBuffer (void) const {
    return new _String ((unsigned long)unitLength, false);
}

//_______________________________________________________________________

_String&     _DataSetFilter::operator () (unsigned long site, unsigned long pos) {
    if (!accessCache || accessCache->sLength != unitLength) {
        if (accessCache) {
            DeleteObject (accessCache);
        }
        accessCache = new _String ((unsigned long)unitLength, false);
    }
    
    long vIndex = theNodeMap.lData[pos];
    if (unitLength==1) {
        accessCache->sData[0]=(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site]]])->sData[vIndex];
    } else {
        site*=unitLength;
        for (int k = 0; k<unitLength; k++) {
            accessCache->sData[k] = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site++]]])->sData[vIndex];
        }
    }
    return *accessCache;
}

//_______________________________________________________________________

const _String     _DataSetFilter::RetrieveState (unsigned long site, unsigned long pos) const {
    _String state ((unsigned long)unitLength, false);
    RetrieveState (site, pos, state, false);
    return state;
    
}

//_______________________________________________________________________

void     _DataSetFilter::RetrieveState (unsigned long site, unsigned long pos, _String& reply, bool map) const
{
    long vIndex = theNodeMap.lData[pos];
    if (map) {
        if (unitLength==1) {
            reply.sData[0]=(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[duplicateMap.lData[site]]]])->sData[vIndex];
        } else {
            site = unitLength*duplicateMap.lData[site];
            for (int k = 0; k<unitLength; k++) {
                reply.sData[k] = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site++]]])->sData[vIndex];
            }
        }
    } else {
        if (unitLength==1) {
            reply.sData[0]=(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site]]])->sData[vIndex];
        } else
            site*=unitLength;
        for (int k = 0; k<unitLength; k++) {
            reply.sData[k] = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site++]]])->sData[vIndex];
        }
    }
}

//_______________________________________________________________________

void _DataSetFilter::GrabSite (unsigned long site, unsigned long pos, _String& storage)
{
    
    long vIndex = theNodeMap.lData[pos];
    if (unitLength==1) {
        storage.sData[0]=(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site]]])->sData[vIndex];
    } else {
        site*=unitLength;
        for (int k = 0; k<unitLength; k++) {
            storage.sData[k] = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site++]]])->sData[vIndex];
        }
    }
}

//_______________________________________________________________________

void _DataSetFilter::GrabSite (unsigned long site, unsigned long pos, char * s)
{
    long vIndex = theNodeMap.lData[pos];
    if (unitLength==1) {
        s[0]=(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site]]])->sData[vIndex];
    } else {
        site*=unitLength;
        for (int k = 0; k<unitLength; k++) {
            s[k] = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site++]]])->sData[vIndex];
        }
    }
}

//_______________________________________________________________________

_SimpleList* _DataSetFilter::CountAndResolve (long pattern, hyFloat * storage, bool randomly)
// last cell in the list contains the count of distinct characters in the column
{
    _SimpleList* resList = new _SimpleList (theNodeMap.lLength+1,0,0),
    counts (dimension,0,0);
    
    checkPointer (resList);
    
    _List        ambStates;
    _String      aState  (unitLength, false);
    
    hyFloat*  freqStorage = storage;
    
    if (!freqStorage) {
        freqStorage = new hyFloat [undimension];
    }
    
    long    normalizingSum = 0,
    charCount      = 0;
    
    for (long k=0; k<theNodeMap.lLength; k++) {
        GrabSite (pattern, k, aState);
        long      characterRes = Translate2Frequencies (aState, freqStorage, true);
        if (characterRes>=0) {
            resList->lData[k] = characterRes;
            
            if (characterRes >= dimension) {
                WarnError (_String("Internal error in _DataSetFilter::CountAndResolve\n"));
            }
            
            if ((counts.lData[characterRes]++) == 0) {
                normalizingSum ++;
            }
            
            charCount ++;
        } else {
            _SimpleList * possibleResolutions = new _SimpleList;
            if (!possibleResolutions) {
                checkPointer (possibleResolutions);
            }
            
            (*possibleResolutions) << k;
            
            for (long m=0; m<dimension; m++)
                if (freqStorage[m]>0.) {
                    (*possibleResolutions) << m;
                }
            
            ambStates.AppendNewInstance (possibleResolutions);
        }
    }
    
    if (normalizingSum > 0) {
        if (ambStates.lLength) {
            _SimpleList  ambResolutions (dimension,0,0);
            for (long t=0; t<ambStates.lLength; t++) {
                _SimpleList * stateResolutions = (_SimpleList*)ambStates(t);
                
                if (!randomly) {
                    long          totalSum = 0,
                    idx = 0;
                    
                    for (long l=1; l<stateResolutions->lLength; l++) {
                        long tmp = counts.lData[stateResolutions->lData[l]];
                        if (tmp>totalSum) {
                            idx = l;
                            totalSum = tmp;
                        }
                    }
                    if (idx > 0)
                        // if no resolutions, resolve randomly
                    {
                        idx = stateResolutions->lData[idx];
                        resList->lData[stateResolutions->lData[0]] = idx;
                        ambResolutions.lData [idx] ++;
                        continue;
                    }
                    
                }
                
                long          totalSum = 0;
                for (long l=1; l<stateResolutions->lLength; l++) {
                    totalSum += counts.lData[stateResolutions->lData[l]];
                }
                
                if (totalSum > 0) {
                    long          randomN = genrand_real2() * totalSum - counts.lData[stateResolutions->lData[1]],
                    ind = 1;
                    
                    while (randomN > 0) {
                        randomN -= counts.lData[stateResolutions->lData[++ind]];
                    }
                    
                    totalSum = stateResolutions->lData[ind];
                } else {
                    long          randomN = genrand_real2() * charCount - counts.lData[0],
                    ind = 0;
                    
                    while (randomN > 0) {
                        randomN -= counts.lData[++ind];
                    }
                }
                resList->lData[stateResolutions->lData[0]] = totalSum;
                ambResolutions.lData [totalSum] ++;
            }
            
            for (long l=0; l<dimension; l++)
                if (ambResolutions.lData[l] && !counts.lData[l]) {
                    normalizingSum ++;
                }
        }
    }
    
    resList->lData[theNodeMap.lLength] = normalizingSum;
    
    if (freqStorage != storage) {
        delete freqStorage;
    }
    
    return       resList;
}

//_______________________________________________________________________

_Matrix* _DataSetFilter::PairwiseCompare (_SimpleList* s1, _SimpleList *s2, _List* labels)
// s1 and s2 are the lists produced by CountAndResolve
// if labels is not nil, then it will receive row and column labels in the contigency table
// the result matrix has rows labeled by states in s1, and columns - by states in s2
{
    long    * sort1 = new long[dimension],
    * sort2 = new long[dimension],
    c = s2->lData[s2->lLength-1];
    
    _Matrix * res   = new _Matrix (s1->lData[s1->lLength-1],c,false,true);
    
    if (sort1 && sort2 && res) {
        for (long k = 0; k<dimension; k++) {
            sort1[k] = -1;
            sort2[k] = -1;
        }
        
        long idx1 = 0,
        idx2 = 0;
        
        _SimpleList  *lbl1 = nil,
        *lbl2 = nil;
        
        if (labels) {
            lbl1 = new _SimpleList;
            lbl2 = new _SimpleList;
            
            checkPointer (lbl1);
            checkPointer (lbl2);
            
            (*labels) << lbl1;
            (*labels) << lbl2;
            
            DeleteObject (lbl1);
            DeleteObject (lbl2);
        }
        
        for (long k2 = 0; k2 < s1->lLength-1; k2++) {
            long c1 = s1->lData[k2],
            c2 = s2->lData[k2];
            
            if (sort1[c1] < 0) {
                sort1[c1] = idx1;
                if (lbl1) {
                    (*lbl1) << c1;
                }
                c1 = idx1++;
            } else {
                c1 = sort1[c1];
            }
            
            if (sort2[c2] < 0) {
                sort2[c2] = idx2;
                if (lbl2) {
                    (*lbl2) << c2;
                }
                c2 = idx2++;
            } else {
                c2 = sort2[c2];
            }
            
            /*if ((c1>=res->GetHDim())||(c2>=res->GetVDim()))
             {
             printf ("\nInternal Error\n");
             }*/
            
            res->theData[c1*c+c2] += 1.;
        }
        
        delete [] sort1;
        delete [] sort2;
    } else {
        checkPointer (nil);
    }
    
    return res;
}

//_______________________________________________________________________

_List *  _DataSetFilter::ComputePatternToSiteMap (void) const {
    _List * result = new _List ();
    
    for (unsigned long k = 0UL; k < theFrequencies.lLength; k++) {
        (*result) < new _SimpleList;
    }
    for (unsigned long s = 0UL; s < duplicateMap.lLength; s++) {
        *((_SimpleList**)result->lData)[duplicateMap.lData[s]] << s;
    }
    return result;
}

//_______________________________________________________________________

char     _DataSetFilter::GetChar (unsigned long site, unsigned long pos)
{
    //long vIndex = theNodeMap.lLength?theNodeMap.lData[pos]:pos;
    return (*theData)(theMap.lData[site],theNodeMap.lData[pos],1);
}
//_______________________________________________________________________

bool     _DataSetFilter::CompareTwoSites (unsigned long site1, unsigned long site2, unsigned long pos1) const {
    pos1 = theNodeMap.lData[pos1];
    
    if (unitLength == 3) { // codon
        site1*=3;
        site2*=3;
        return
        ((((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1]]])->sData[pos1]==
         ( ((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site2]]])->sData[pos1])
        &&((((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1+1]]])->sData[pos1]==
           (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site2+1]]])->sData[pos1])
        &&((((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1+2]]])->sData[pos1]==
           (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site2+2]]])->sData[pos1]);
    } else {
        site1*=unitLength;
        site2*=unitLength;
        unsigned long k = 0UL;
        
        /*if ((((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1]]])->sLength<=pos1)
         {
         printf ("(%d)%s\n(%d)%s\n",site1,(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1]]])->sData,
         site2,(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site2]]])->sData);
         FlagError ("Internal DataSetFilter bug\n");
         }*/
        
        for (; k<unitLength; k++) {
            if ((((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1]]])->sData[pos1]!=
                (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site2]]])->sData[pos1]) {
                break;
            }
            site1++;
            site2++;
        }
        if (k==unitLength) {
            return true;
        }
    }
    return false;
}

//_______________________________________________________________________

bool     _DataSetFilterNumeric::CompareTwoSites (unsigned long, unsigned long, unsigned long) const {
    return false;
}


//_______________________________________________________________________

bool     _DataSetFilter::CompareTwoSitesChar (unsigned long site1, unsigned long site2, unsigned long pos1) const {
    
    //  long *fI = theMap.quickArrayAccess();
    
    //  if (theNodeMap.lLength)
    //  {
    pos1 = theNodeMap(pos1);
    //  }
    //  else
    //  {
    //      vIndex1=pos1;
    //  }
    //  return ((*theData)(fI[site1],pos1, 1)==(*theData)(fI[site2],pos1, 1));
    return ((*theData)(theMap.lData[site1],pos1, 1)==(*theData)(theMap.lData[site2],pos1, 1));
}
//_______________________________________________________________________
long    _DataSetFilter::SiteFrequency (unsigned long site)
{
    return theFrequencies.lData[site];
}

//_______________________________________________________________________
bool    _DataSetFilter::HasDeletions (unsigned long site, _AVLList* storage)
{
    long        loopDim  = GetDimension();
    hyFloat* store    = new hyFloat [loopDim];
    
    long j,
    upTo = theNodeMap.lLength?theNodeMap.lLength:theData->NoOfSpecies();
    
    bool outcome = false;
    
    for (unsigned int k = 0; k<upTo; k++) {
        Translate2Frequencies ((*this)(site,k), store, false);
        
        bool oneF = false,
        zeroF = false;
        
        for (j=0; j<loopDim; j++) {
            if (store[j]==0.0) {
                zeroF = true;
            } else if (store[j]==1.0) {
                oneF = true;
            }
        }
        if (!(oneF&&zeroF)) {
            if (storage) {
                outcome = true;
                storage->Insert ((BaseRef)theNodeMap.lData[k]);
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
bool    _DataSetFilter::IsConstant (unsigned long site,bool relaxedDeletions)
{
    hyFloat *store = new hyFloat [GetDimension()],
    *store2 = new hyFloat [GetDimension()];
    
    unsigned long j,
    upTo = theNodeMap.lLength?theNodeMap.lLength:theData->NoOfSpecies(),
    loopDim = GetDimension();
    
    Translate2Frequencies ((*this)(site,0), store, false);
    
    if (relaxedDeletions) {
        for (unsigned long k = 1UL; k<upTo; k++) {
            Translate2Frequencies ((*this)(site,k), store2, false);
            for (j=0UL; j<loopDim; j++) {
                if (store2[j]==0.0) {
                    store[j]=0.0;
                }
            }
        }
        for (j=0UL; j<loopDim; j++)
            if (store[j]!=0.0) {
                break;
            }
        
        delete [] store;
        delete [] store2;
        return j!=loopDim;
        
    } else {
        for (unsigned long k = 1; k<upTo; k++) {
            Translate2Frequencies ((*this)(site,k), store2, false);
            for (j=0UL; j<loopDim; j++)
                
                if (store[j]!=store2[j]) {
                    delete [] store;
                    delete [] store2;
                    return false;
                }
        }
    }
    
    delete [] store;
    delete [] store2;
    return true;
}

//_______________________________________________________________________

_Matrix*        _DataSetFilter::GetFilterCharacters (bool flip) const {
    long        unit_length = GetUnitLength (),
    seq_length  = flip?theFrequencies.lLength:(GetSiteCount () / unitLength),
    seq_count   = NumberSpecies();
    
    _List       result;
    
    _String      char_buffer (unit_length ,false);
    
    if (flip) {
        for (long k=0; k< seq_length; k++) {
            _String *alignment_column = new _String (seq_count+1,true);
            for (long k2=0; k2< seq_count ; k2++) {
                RetrieveState(k,k2,char_buffer,false);
                (*alignment_column) << char_buffer;
            }
            alignment_column->Finalize();
            result  < alignment_column;
        }
    } else
        for (long k=0; k < seq_count; k++) {
            result      < GetSequenceCharacters(k);
        }
    
    return new _Matrix (result);
}

//_______________________________________________________________________

_String*        _DataSetFilter::GetSequenceCharacters (long seqID)  const{
    long            unitSizeL   = GetUnitLength();
    
    _String * aSequence = new _String (GetSiteCount(),true);
    
    if (seqID >= 0 && seqID < theNodeMap.lLength) {
        _String      aState (unitSizeL,false);
        unsigned long        upTo = GetSiteCountInUnits();
        for (unsigned long k2=0UL; k2<upTo; k2++) {
            RetrieveState(k2,seqID,aState);
            (*aSequence) << aState;
        }
    }
    aSequence->Finalize();
    return aSequence;
}

//_______________________________________________________________________

_String*        _DataSet::GetSequenceCharacters (long seqID)  const{
    
    unsigned long        upTo = NoOfColumns();
    _String * aSequence = new _String (upTo,true);
    
    if (seqID >= 0 && seqID < noOfSpecies) {
        for (unsigned long k2=0UL; k2<upTo; k2++) {
            (*aSequence) << GetSite (k2)->getChar (seqID);
        }
    }
    aSequence->Finalize();
    return aSequence;
}


//_______________________________________________________________________
long    _DataSetFilter::HasExclusions (unsigned long site, _SimpleList* theExc, hyFloat*store )
{
    long   filterDim = GetDimension(false);
    
    if (theNodeMap.lLength)
        for (unsigned long k = 0; k<theNodeMap.lLength; k++) {
            Translate2Frequencies   ((*this)(site,k), store, false);
            
            long                    j                       = 0,
            s                      = 0;
            
            for (j=0; j<filterDim; j++)
                if (store[j] > 0.0) {
                    s++;
                    if (theExc->Find(j) < 0) {
                        break;
                    }
                }
            
            if (j == filterDim && s) {
                return k;
            }
        }
    
    return -1;
}
//_______________________________________________________________________
void    _DataSetFilter::Freeze (long site)
{
    for (int k = 0; k<unitLength; k++) {
        _Site* tC = (_Site*)((*(_List*)theData)(theData->theMap(this->theMap(site*unitLength+k))));
        tC->SetRefNo(-1);
        tC->PrepareToUse();
    }
}

//_______________________________________________________________________
void    _DataSetFilter::UnFreeze (long site)
{
    for (int k = 0; k<unitLength; k++) {
        _Site* tC = (_Site*)((*(_List*)theData)(theData->theMap(this->theMap(site*unitLength+k))));
        tC->SetRefNo(0);
        //      tC->Archive();
    }
}

//_________________________________________________________
_Matrix* _DataSetFilter::ComputePairwiseDifferences (long i, long j, _hy_dataset_filter_ambiguity_resolution resolution_option) const {
    
    try {
        
        if (unitLength > 3) {
            throw _String("ComputePairwiseDifferences is not implemented for data filters with unit size > 3");
        }
        
        long    mxDim      = GetDimension (true);
        
        _Matrix     *res   = new _Matrix  (mxDim,mxDim,false,true);
        
        hyFloat  *sm1   = new hyFloat[mxDim],
        *sm2   = new hyFloat[mxDim];
        
        
        
        _String      state1 (unitLength,false),
        state2 (unitLength,false);
        
        
        if (conversionCache.lLength == 0) {
            throw _String ("ComputePairwiseDifferences called on a filter with emptyString conversionCache");
        }
        
        long        *tcodes  = conversionCache.lData+89,
        *ccodes  = conversionCache.lData+1,
        ccount   = conversionCache.lData[0];
        
        for (unsigned long site_pattern = 0UL; site_pattern < theFrequencies.lLength; site_pattern++) {
            long s1 = -1, s2 = -1;
            
            int c1, c2;
            
            c1 = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[unitLength*site_pattern]]])->sData[theNodeMap.lData[i]],
            c2 = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[unitLength*site_pattern]]])->sData[theNodeMap.lData[j]];
            
            if (unitLength == 1) {
                s1 = conversionCache.lData[(c1-40)*(undimension+1)+undimension],
                s2 = conversionCache.lData[(c2-40)*(undimension+1)+undimension];
            } else {
                int         c12 = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[unitLength*site_pattern+1]]])->sData[theNodeMap.lData[i]],
                c22 = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[unitLength*site_pattern+1]]])->sData[theNodeMap.lData[j]];
                
                
                state1.sData[0] = c1;
                state1.sData[1] = c12;
                
                state2.sData[0] = c2;
                state2.sData[1] = c22;
                
                c1  = ccodes[c1-40];
                c12 = ccodes[c12-40];
                
                c2  = ccodes[c2-40];
                c22 = ccodes[c22-40];
                
                if (unitLength == 2) {
                    if ((c1>=0)&&(c12>=0)) {
                        s1 = tcodes[c1*ccount+c12];
                    }
                    
                    if ((c2>=0)&&(c22>=0)) {
                        s2 = tcodes[c2*ccount+c22];
                    }
                } else {
                    int         c13 = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[unitLength*site_pattern+2]]])->sData[theNodeMap.lData[i]],
                    c23 = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[unitLength*site_pattern+2]]])->sData[theNodeMap.lData[j]];
                    
                    //printf ("\n%c %c", c13, c23);
                    
                    state1.sData[2] = c13;
                    state2.sData[2] = c23;
                    
                    c13 = ccodes[c13-40];
                    c23 = ccodes[c23-40];
                    
                    //printf (" %d %d %s %s\n", c13, c23, state1.sData, state2.sData);
                    
                    if ((c1>=0)&&(c12>=0)&&(c13>=0)) {
                        s1 = tcodes[ccount*(c1*ccount+c12)+c13];
                    }
                    
                    if ((c2>=0)&&(c22>=0)&&(c23>=0)) {
                        s2 = tcodes[ccount*(c2*ccount+c22)+c23];
                    }
                }
            }
            /*
             if (pcAmbiguitiesAverage.Equal (resFlag)) {
             res = dsf->ComputePairwiseDifferences (seq,site,1);
             } else if (pcAmbiguitiesResolve.Equal (resFlag)) {
             res = dsf->ComputePairwiseDifferences (seq,site,2);
             } else if (pcAmbiguitiesSkip.Equal (resFlag)) {
             res = dsf->ComputePairwiseDifferences (seq,site,3);
             } else {
             res = dsf->ComputePairwiseDifferences (seq,site,0);
             }*/
            
            
            if (s1>=0 && s2>=0) { // one to one
                res->theData[s1*mxDim+s2] += theFrequencies.lData[site_pattern];
            } else {
                if (resolution_option != kAmbiguityHandlingSkip) {
                    _Matrix * freqsAtSite = nil;
                    
                    if (resolution_option != kAmbiguityHandlingResolve) {
                        _SimpleList   //seqList,
                        siteList;
                        
                        
                        for (long si = 0; si < unitLength; si++) {
                            siteList << theMap.lData[unitLength*site_pattern+si];
                        }
                        
                        _SimpleList copy_node_oder (theNodeMap);
                        freqsAtSite     = theData->HarvestFrequencies (unitLength, unitLength, 0, copy_node_oder, siteList);
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
        WarnError (error);
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
        CreateMatrix (&target,1,7,false,true,false);
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
 // use a standard convert to frequencies function to check whether a character has a unique convertion
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
    long  store      [HYPHY_SITE_DEFAULT_BUFFER_SIZE],
    resolution_count  = -1L;
    
    
    InitializeArray(parvect, dimension, 0.);
    
    if (unitLength == 1) {
        resolution_count = theData->theTT->TokenResolutions (str.sData[0],store,smear);
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
        resolution_count = theData->theTT->TokenResolutions (str.sData[0],store);
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
void    _DataSetFilter::SetupConversion (void)
{
    if (conversionCache.lLength) {
        return;
    }
    
    if ( unitLength==1 ) { // do stuff
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
            unsigned long alphabet_dim = alphabet.sLength;
            
            
            long  ccache [88],
            uncorrected_dimension = GetDimension(false) ;
            
            conversionCache.RequestSpace (89+uncorrected_dimension);
            conversionCache << alphabet.sLength;
            
            for (unsigned long i=0UL; i<88; i++) {
                ccache[i] = -1;
            }
            for (unsigned long i=0UL; i<alphabet.sLength; i++) {
                ccache [alphabet.sData[i]-40] = i;
            }
            for (unsigned long i=0UL; i<88; i++) {
                conversionCache << ccache[i];
            }
            
            _String s (unitLength,false);
            for (unsigned long char_index = 0; char_index < uncorrected_dimension; char_index++ ) {
                _SimpleList components = SplitIntoDigits (char_index, unitLength, alphabet_dim);
                for (unsigned long position = 0; position < unitLength; position ++) {
                    s.sData[position] = alphabet.sData[components (position)];
                }
                conversionCache << MapStringToCharIndex(s);
            }
        }
    }
}

void printDSFilter(_DataSetFilter *d);
