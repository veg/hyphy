/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon    (apoon@cfenet.ubc.ca)
  Steven Weaver (sweaver@ucsd.edu)
  
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

#ifndef     __DATASETFILTER__
#define     __DATASETFILTER__ 

#include "dataset.h"
#include "site.h"


class _Site;

class _DataSetFilter:public BaseObj
{

public:

    _DataSetFilter               (void);
    _DataSetFilter              (_DataSet*, char, _String&);

    virtual                     ~_DataSetFilter (void);

    virtual  BaseRef            toStr           (void);  // convert to string
    virtual  void               toFileStr       (FILE*); // convert to string

    virtual  BaseRef            makeDynamic     (void);
    virtual  long               FreeUpMemory    (long);
    virtual  bool               IsNormalFilter  (void) {
        return true;
    }

    void                CopyFilter      (_DataSetFilter*);


//  void     SetFilter  (_DataSet*, char, _String&);

    void     SetFilter (_DataSet*, char, _SimpleList&, _SimpleList&, bool isFilteredAlready = false);

    void     SetExclusions (_String*, bool = true);

    _String* GetExclusions (void);

    void     SetMap  (_String&); // used to allow nonsequential maps to tree leaves

    void     SetMap  (_SimpleList& newMap) {
        theNodeMap.Clear();    // used to allow nonsequential maps to tree leaves
        theNodeMap.Duplicate(&newMap);
    }

    long     NumberDistinctSites (void) {
        return theFrequencies.lLength;
    }

    long     NumberSpecies (void) {
        return theNodeMap.lLength;
    }

    virtual  long
    GetFullLengthSpecies (void) {
        return theOriginalOrder.lLength;
    }

    virtual  long
    GetSiteCount        (void) {
        return duplicateMap.lLength;
    }

    long     GetFrequency  (long i) {
        return theFrequencies(i);
    }

    long     GetUnitLength  (void) {
        return unitLength;
    }

    virtual  long    GetDimension (bool correct = true);

    long     GetOriginalToShortMap (long i);

    void     ComputePairwiseDifferences (_Matrix&, long, long);
    _Matrix* ComputePairwiseDifferences (long, long, char = 0);

    BaseRef  GetMap (void) {
        return theNodeMap.lLength?&theNodeMap:NULL;
    }

    virtual  _String&   operator () (unsigned long site, unsigned long pos);
    // site indexes unique sites

    virtual  void   RetrieveState (unsigned long site, unsigned long pos, _String&, bool = true);
    // site indexes all sites, including duplicates

    void     FindAllSitesLikeThisOne (long, _SimpleList&);

    _String  GenerateConsensusString (_SimpleList* =nil);

    void     GrabSite (unsigned long,unsigned long,_String&);
    void     GrabSite (unsigned long,unsigned long,char*);

    virtual  char      GetChar(unsigned long site, unsigned long pos);
    long       SiteFrequency  (unsigned long site);
    bool       HasDeletions   (unsigned long site, _AVLList* = nil);
    long       HasExclusions  (unsigned long site, _SimpleList* theExc, _Parameter *buffer);
    bool       IsConstant     (unsigned long site,  bool relaxedDeletions = true);

    long     Translate2Frequencies (_String&, _Parameter*, bool);
    long     MapStringToCharIndex  (_String&);
    //long   Translate2Frequencies (char*,    _Parameter*, bool = true);
    long     Translate2Frequencies (char,     _Parameter*, bool);

    _Matrix* HarvestFrequencies (char unit, char atom, bool posSpec, bool = true);

    void     Freeze (long);

    void     UnFreeze (long);

    void     MatchStartNEnd (_SimpleList&, _SimpleList&, _SimpleList* = nil);

    _String* GetSequenceName(long idx) {
        return (_String*)(theData->GetNames ()(theNodeMap.lData[idx]));
    }

    _String* GetSequenceCharacters
    (long);

    _DataSet*                       GetData                     (void) {
        return theData;
    }
    void                            SetData                     (_DataSet* ds) {
        theData = ds;
    }
    _String                         ConvertCodeToLetters        (long code, char base) {
        return theData->theTT->ConvertCodeToLetters(code,base);
    }

    void                            ConvertCodeToLettersBuffered(long code, char base, char *, _AVLListXL* );
    // 20090212: SLKP
    // added this function to cache repeated character code -> string conversions
    // and to skip returning temp objects but simply writing to buffer

    /**
    * Find all unique sequences in the data filter. 
    *
    * \n Usage: FindDuplicateSequences(uniqueIndex, instanceCount, true);
    * @author SLKP
    * @param indices For each sequence - the list of indices corresponding to the unique strings
                     For example, if sequence 1 == sequence 3 and sequence 4 == sequence 5 this list 
                     will contain 0,1,3 
    * @param map The index of the unique string to which the current string is mapped
                     For example, if sequence 1 == sequence 3 and sequence 4 == sequence 5 this list 
                     will contain 0,1,0,2,2 
    * @param counts  The number of copies for each unique string found
                     For example, if sequence 1 == sequence 3 and sequence 4 == sequence 5 this list 
                     will contain 2,1,2
    * @param strict  Controls is the strings must match exactly (0), exactly + gap (1), via the superset rule (2) or via the partial match rule (3).
                     Nucleotide letters A and - (or ?) (IUPAC code for A or G) will count as mismatched for mode 0 and matched for mode 1.
                     Nucleotide letters A and R (IUPAC code for A or G) will count as mismatched for mode 0 and matched for
                     modes 2 and 3 because R is a superset of A. (note that R matches R in all modes, even though the letter is
                     an ambiguous nucleotide).
                     Nucleotide letters R (A or G) and M (A or C) will match under mode 3 (because they both encode A as an option), 
                     but not under modes 0-2. 
                     Match in mode (0) => match in mode (1) => match in mode (2) => match in mode (3).
                    
    * @return The number of unique sequences. 
    */
    unsigned long                   FindUniqueSequences      (_SimpleList&, _SimpleList&, _SimpleList&, short = 0);


    long                            CorrectCode                 (long code);
    virtual  bool                   CompareTwoSites             (unsigned long, unsigned long,unsigned long);
    bool                            CompareTwoSitesChar         (unsigned long, unsigned long,unsigned long);
    long                            FindSpeciesName             (_List&, _SimpleList&);
    _DataSetFilter*                 PairFilter                  (long, long, _DataSetFilter*);
    void                            SetDimensions               ();
    long                            LookupConversion            (char c, _Parameter* receptacle);
    void                            SetupConversion             (void);
    void                            FilterDeletions             (_SimpleList* theExc = nil);
    _Matrix*                        GetFilterCharacters         (bool = false);
    _SimpleList*                    CountAndResolve             (long, _Parameter* = nil, bool = false);
    _Matrix*                        PairwiseCompare             (_SimpleList*, _SimpleList*, _List* = nil);

    _List*                          ComputePatternToSiteMap     (void);
    // 20090206: SLKP
    // a utility function to return a _List of simplelists (one per unique site pattern) that provides an ordered list of
    //           the indices of all sites that have the same pattern in the original alignment

    void                            PatternToSiteMapper         (void*, void*, char, long);
    /*
        20090325: SLKP
        a function that takes per pattern values (source, argument 1)
        and maps them onto sites into target (argument 2)
        the third argument is 0 to treat the pointers as _Parameter*
        1 to treat them as long*
        2 and to treat them as _Parameter* and long*, respetively
        20090929: SLKP
        the fourth argument is used to speficy a padding-size,
            all values from the filter size up to that value are set to 1 (for mode 0) and 0 (for mode 1)
            this is needed to handle uneven data filters in SITE_LIKELIHOOD constructs
     */

    _SimpleList
    theFrequencies,
    theNodeMap,
    theMap,
    theOriginalOrder,
    theExclusions,
    duplicateMap;

    char*    GetColumn (long index) {
        return ((_Site*)(((BaseRef*)theData->lData)[theData->theMap.lData[theMap.lData[index]]]))->sData;
    }

    _SimpleList     conversionCache;

protected:

    char            unitLength;
    long            dimension;

private:

    void            internalToStr (FILE*,_String&);
    _String*        accessCache;

    long            undimension;

    _DataSet*       theData;
//      _SimpleList     conversionCache;

    void            XferwCorrection (_Matrix& , _Parameter*, long);
    void            XferwCorrection (_Parameter* , _Parameter*, long);
    void            XferwCorrection (long* , _Parameter*, long);
};

#endif
