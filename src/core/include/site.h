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

#ifndef _GENSITE_
#define _GENSITE_
//#pragma once
#include "sequence.h"
#include "parser.h"
#include "simplelist.h"
#include "list.h"
#include "avllist.h"
#include "avllistx.h"
#include "avllistxl.h"
#include "stdlib.h"

#define   NUCLEOTIDEDATA 0
#define   CODONDATA      1





//_________________________________________________________
class _TranslationTable:public BaseObj
{

public:

    _TranslationTable                       (void);
    _TranslationTable                       (char);
    _TranslationTable                       (_String&);
    /* 20100618: SLKP

            - new constructor (needed to handle ExecuteCase52 / Simulate properly)
              which takes an alphabet string and checks to see if it's a standard one
              DNA/RNA/Protein or Binary

     */
    _TranslationTable                       (_TranslationTable&);
    virtual ~_TranslationTable              (void) {
        if (checkTable) {
            free (checkTable);
        }
    }
    virtual BaseRef  makeDynamic            (void);

    long    TokenCode                       (char);
    char    CodeToLetter                    (long*);

    void    AddBaseSet                      (_String&);
    bool    TokenCode                       (char, long*, bool = true);
    void    SplitTokenCode                  (long, long*);

    void    AddTokenCode                    (char, _String&);
    void    PrepareForChecks                (void);
    bool    IsCharLegal                     (char);
    char    GetSkipChar                     (void);
    char    GetGapChar                      (void);
    _String ConvertCodeToLetters            (long, char);
    long    LengthOfAlphabet                (void);
    bool    IsStandardBinary                (void) {
        return baseLength==2 && baseSet.sLength==0;
    }
    bool    IsStandardNucleotide            (void) {
        return baseLength==4 && baseSet.sLength==0;
    }
    bool    IsStandardAA                    (void) {
        return baseLength==20&& baseSet.sLength==0;
    }
    _TranslationTable*
    MergeTables                     (_TranslationTable*);

    char                                    baseLength;
    // number of "fundamental" tokens
    //(4 for nucl, ACGT; 20 for amino acids)


    _String                                 tokensAdded,
                                            baseSet;

    _SimpleList                             translationsAdded;
    char*                                   checkTable;
    // if null - then assume default translation table;
};

//_________________________________________________________

// data set file state data struct
struct FileState {

    _TranslationTable* translationTable;
    long
    curSpecies,
    totalSpeciesRead,
    totalSitesRead,
    totalSpeciesExpected,
    totalSitesExpected,
    curSite,
    maxStringLength,
    pInSrc;
    bool
    acceptingCommands,
    allSpeciesDefined,
    interleaved,
    autoDetect,
    isSkippingInNEXUS;
    int
    fileType,
    baseLength;
    char
    repeat,
    skip;

    _String
    *theSource,
    *theNamespace;

    _SimpleList
    rawLinesFormat;
};
//_________________________________________________________

class _Site:public _CString   // compressible string
{

public:
    _Site (void);
    //does nothing
    _Site (_String&);
    // data constructor
    _Site (char);
    // data constructor
    _Site (long);
    // reference constructor

    virtual     ~_Site (void);
    //destructor

    void    Complete (void); // mark this site as complete and compress it

    virtual     BaseRef makeDynamic(void);
    virtual     void    Duplicate  (BaseRef);
    virtual     void    Clear  (void);

    void    PrepareToUse (void); // decompress the site preparing for intensive use
    void    Archive      (void); // archive the site for later use

    long    GetRefNo (void) {
        return refNo<0?-refNo-2:refNo-2;
    }

    bool    IsComplete (void) {
        return refNo<0;
    }

    void    SetRefNo (long r) {
        refNo = -r-2;
    }


private:

    long             refNo;      // if this site contains a reference to another one
    // if refNo is negative, then shows whether the definition of this datatype has been completed
};

//_________________________________________________________

// data set file state data struct
struct _DSHelper {

    _SimpleList characterPositions;
    _List       incompletePatternStorage;
    _AVLListX*  incompletePatterns;

    _DSHelper(void) {
        incompletePatterns = new _AVLListX (&incompletePatternStorage);
        checkPointer (incompletePatterns);
    }
    ~_DSHelper(void) {
        DeleteObject (incompletePatterns);
    }
};

//_________________________________________________________

class   _DataSet:public _List // a complete data set
{
public:

    _DataSet                (void);
    _DataSet                (long);
    _DataSet                (FILE*);
    // with estimated number of sites per file
    virtual             ~_DataSet               (void);

    virtual  BaseRef    makeDynamic             (void);

    void        AddSite                 (char);

    void        Write2Site              (long, char);
    void        CheckMapping            (long);

    void        Finalize                (void);
    // remove duplicate data types and compress

    long        GetNoTypes              (void);
    // return the number of unique sites

    long        GetCharDimension        (void);
    // return the size of the alphabet space

    long        GetFreqType             (long);
    // return the frequency of a site

    _Site*      GetSite                 (long index) {
        return ((_Site**)lData)[theMap.lData[index]];
    }

    long        ComputeSize             (void);
    // compute the size of this object in memory

    void        Clear                   (void);

    virtual  char       operator ()             (unsigned long, unsigned long, unsigned int);
    // retrieve element pos of site-th site

    virtual  BaseRef    toStr                   (void);
    // convert to string

    virtual  void       toFileStr               (FILE*dest);

    void        Compact                 (long);
    // release string overhead
    void        ConvertRepresentations  (void);

    _Matrix *           HarvestFrequencies      (char, char, bool, _SimpleList&, _SimpleList&, bool = true);
    // this function counts observed frequencies of elements in a data set
    // unit is the length of an info unit (nucl - 1, codons - 3)
    // atom is the "minimal" countable element (nucl - 1, codons - 1)
    // posSpec - if position of an atom within an item is to be accounted for
    // segmentation - partition of the underlying DataSet to look at
    // null for segmentation assumes the entire dataset

    void        MatchIndices            (_Formula&, _SimpleList& , bool , long );
    friend   void       printFileResults        (_DataSet* );
    char        InternalStorageMode     (void) {
        return useHorizontalRep;
    }

    long        NoOfSpecies             (void) {
        return noOfSpecies;
    }

    long        NoOfColumns             (void) {
        return theMap.lLength;
    }
    long        NoOfUniqueColumns       (void) {
        return lLength;
    }
    void        AddName                 (_String&);
    _List&     GetNames             (void) {
        return theNames;
    }
    _SimpleList&
    GetTheMap               (void) {
        return theMap;
    }
    void        FindAllSitesLikeThisOne (long, _SimpleList&);

    friend      class       _DataSetFilter;
    friend      _DataSet*    ReadDataSetFile        (FILE*,char,_String*,_String*, _String*,_TranslationTable*);
    friend      long         ProcessLine            (_String&s , FileState *fs, _DataSet& ds);

    static      _DataSet*    Concatenate            (_SimpleList);
    static      _DataSet*    Combine                (_SimpleList);

    static      _TranslationTable*
    CheckCompatibility(_SimpleList& ref, char concatOrCombine);


    void         ProcessPartition       (_String&, _SimpleList&,  bool, _SimpleList* = nil, _SimpleList* = nil);
    void         SetTranslationTable    (_DataSet *  newTT );
    void         SetTranslationTable    (_TranslationTable *  newTT );
    _TranslationTable*
    GetTT                   (void) {
        return theTT;
    }
    _Parameter   CheckAlphabetConsistency
    (void);

    void         SetNoSpecies           (long n) {
        noOfSpecies = n;
    }
    void         ResetIHelper           (void);
private:

     void        constructFreq           (long*, _Parameter *, char, long, long, int, int, int);

    _SimpleList theMap,
                theFrequencies;         // remapping vector, and the counter of frequencies

    unsigned int
    noOfSpecies;

    _TranslationTable*
    theTT;                  // translation Table, if any

    _List       theNames;               // Names of species
    FILE*       streamThrough;

    _DSHelper*  dsh;
    bool        useHorizontalRep;

};

//_________________________________________________________

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



//_________________________________________________________

class _DataSetFilterNumeric:public _DataSetFilter
{

public:

    _DataSetFilterNumeric                       (void) {
    }
    _DataSetFilterNumeric                       (_Matrix*, _List&,_DataSet*,long);

    virtual  bool                                       IsNormalFilter          (void) {
        return false;
    }


    virtual  BaseRef                                    makeDynamic             (void);
    virtual  long                                       GetDimension            (bool) {
        return dimension;
    }

    _Parameter*                 getProbabilityVector    (long,long,long = 0);
    virtual  bool               CompareTwoSites         (unsigned long, unsigned long,unsigned long);

    long                    shifter,
                            categoryShifter,
                            categoryCount;
    _Matrix                 probabilityVectors;
    // N x M dense matrix
    // N = spec count
    // M = unique sites * dimension

};


//_________________________________________________________

extern          _TranslationTable       defaultTranslationTable;

void            ReadNextLine            (FILE* fp, _String *s, FileState* fs, bool append = false, bool upCase = true);
_DataSet*       ReadDataSetFile         (FILE*, char = 0, _String* = nil, _String* = nil, _String* = nil,_TranslationTable* = &defaultTranslationTable);
void            fillDefaultCharTable    (void);
void            printFileResults        (_DataSet*);
void            printDSFilter           (_DataSetFilter* d);

bool            StoreADataSet           (_DataSet*, _String*);


extern _String  dataFileTree,
       dataFileTreeString,
       aminoAcidOneCharCodes,
       dnaOneCharCodes,
       rnaOneCharCodes,
       binaryOneCharCodes,
       nexusFileTreeMatrix,
       dataFilePartitionMatrix,
       defaultLargeFileCutoff,
       nexusBFBody;

extern _DataSet*lastNexusDataMatrix;

#endif
