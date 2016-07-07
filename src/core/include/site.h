/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
   Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
   Art FY Poon    (apoon@cfenet.ubc.ca)
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

#define   HY_TRANSLATION_TABLE_DNA      0x01
#define   HY_TRANSLATION_TABLE_RNA      0x02
#define   HY_TRANSLATION_TABLE_BINARY   0x04
#define   HY_TRANSLATION_TABLE_PROTEIN  0x08


//_________________________________________________________
class _TranslationTable:public BaseObj {
  
private:
    static _List _list_of_default_tables;

public:

    _TranslationTable                       (void);
    _TranslationTable                       (unsigned char);
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

    long    TokenCode                       (char) const;
    char    AmbigToLetter                    (long*, unsigned long ) const;

    void    AddBaseSet                      (_String&);
    void    SplitTokenCode                  (long, long*) const;

    long    TokenResolutions                (char token, long * buffer, bool gap_to_one = true) const;
    /**
      Given a character in this translation table, return the number of
      base alphabet characters that map to it, and populate buffer with
      their codes. For example `TokenResolutions ('T', buffer) will
      return 2 and set buffer[0] = 1, buffer[1] = 3, assuming the translation
      table has the standard IUPAC nucleotdie code
     
      @param token the character (unique or ambiguous) to translate
      @param buffer store the resolved characters (up to baseLength) here
      @param gap_to_one if `true`, map gaps (or equivalents) to fully ambiguous characters
     
      @return the number of resolutions
     
    */

    long    MultiTokenResolutions                (const _String& tokens, long * buffer, bool gap_to_one = true) const;
    /**
     Given a string of several in this translation table, return the unique resolutions of 
     characters to base, and populate buffer with
     their codes. For example `TokenResolutions ('ATR', buffer) will
     return 2 and set buffer[0] = 12 [0*16+3*4+0], buffer[1] = 14 [0*16+3*4+2], assuming the translation
     table has the standard IUPAC nucleotdie code. Passing NULL as buffer will result in 
     returning the code for the resolution (if UNIQUE), otherwise -1.
     
     @param tokens the characters (unique or ambiguous) to translate
     @param buffer store the resolved characters (up to baseLength) here [must be at least baseLength ^ length (token) long]
            can be set to NULL in which case the return behavior is modified
     @param gap_to_one if `true`, map gaps (or equivalents) to fully ambiguous characters
     
     @return the number of resolutions OR (if buffer == NULL) the code for the SINGLE resolution or -1 (multiple or invalid resolutions)
    */

    void    AddTokenCode                    (char, _String&);
    void    PrepareForChecks                (void);
    bool    IsCharLegal                     (char);
    char    GetSkipChar                     (void);
    char    GetGapChar                      (void) const;
    const _String ConvertCodeToLetters            (long, unsigned char) const;
    long    LengthOfAlphabet                (void) const;
    bool    IsStandardBinary                (void) const {
        return baseLength==2 && baseSet.sLength==0;
    }
    bool    IsStandardNucleotide            (void) const {
        return baseLength==4 && baseSet.sLength==0;
    }
    bool    IsStandardAA                    (void) const {
        return baseLength==20&& baseSet.sLength==0;
    }
  
    const _String&   ExpandToken            (char token) const;
  
    /** given a (possibly) ambiguous character 
        expand it to a string to equivalent characters
        e.g. ExpandToken ('R') -> "AG" for IUPAC nucleotide data
     
        @param token the character to expand
        @return a string of complete expansions
     
     */
  
    const _String& GetAlphabetString (void) const;

    /** return the alphabet string for this table, so that code 0 maps to result[0], etc
     *
     * @return the alphabet string (e.g. ACGT for standard DNA)
    */
  
    _TranslationTable*
    MergeTables                     (_TranslationTable*);
  
    static const _String&                   GetDefaultTable (long tableType);

    unsigned char                                    baseLength;
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

    long        GetNoTypes              (void) const;
    // return the number of unique sites

    long        GetCharDimension        (void);
    // return the size of the alphabet space

    long        GetFreqType             (long);
    // return the frequency of a site

    _Site*      GetSite                 (long index) const {
        return ((_Site**)lData)[theMap.lData[index]];
    }

    long        ComputeSize             (void);
    // compute the size of this object in memory

    void        Clear                   (bool = true);

    virtual  char       operator ()             (unsigned long, unsigned long, unsigned int) const;
    // retrieve element pos of site-th site

    virtual  BaseRef    toStr                   (unsigned long = 0UL);
    // convert to string

    virtual  void       toFileStr               (FILE*dest, unsigned long = 0UL);

    void        Compact                 (long);
    // release string overhead
    void        ConvertRepresentations  (void);

    _Matrix *           HarvestFrequencies      (unsigned char, unsigned char, bool, _SimpleList&, _SimpleList&, bool = true) const;
    // this function counts observed frequencies of elements in a data set
    // unit is the length of an info unit (nucl - 1, codons - 3)
    // atom is the "minimal" countable element (nucl - 1, codons - 1)
    // posSpec - if position of an atom within an item is to be accounted for
    // segmentation - partition of the underlying DataSet to look at
    // null for segmentation assumes the entire dataset

    void        MatchIndices            (_Formula&, _SimpleList& , bool , long, _String const* = nil) const;
    friend   void       printFileResults        (_DataSet* );
    char        InternalStorageMode     (void) const{
        return useHorizontalRep;
    }

    unsigned long        NoOfSpecies             (void) const{
        return noOfSpecies;
    }

    unsigned long        NoOfColumns             (void) const{
        return theMap.lLength;
    }
    unsigned long        NoOfUniqueColumns       (void) const{
        return lLength;
    }
    void        AddName                 (_String const&);
    void        InsertName              (_String const& name, long where);
  
    _String*   GetSequenceName      (unsigned long i) const {
      return (_String*)theNames.GetItem (i) ;
    }
  
    _List const&     GetNames             (void) const {
        return theNames;
    }
  
    void  ClearNames (void) {
        theNames.Clear();
    }
    
    bool   SetSequenceName (long index, _String * new_name) {
      if (index >= 0L && index < theNames.lLength) {
        theNames.Replace (index, new_name, false);
        return true;
      }
      return false;
    }
  
    void  SetNames (_List const& copy_from) {
      theNames.Clear(); theNames << copy_from;
    }
    
  
    _SimpleList&
    GetTheMap               (void) {
        return theMap;
    }
    void        FindAllSitesLikeThisOne (long, _SimpleList&);

    friend      class       _DataSetFilter;
    friend      _DataSet*    ReadDataSetFile        (FILE*,char,_String*,_String*, _String*,_TranslationTable*);
    friend      long         ProcessLine            (_String&s , FileState *fs, _DataSet& ds);

    static      _DataSet*    Concatenate            (const _SimpleList&);
    static      _DataSet*    Combine                (const _SimpleList&);

    static      _TranslationTable*
    CheckCompatibility(_SimpleList const& ref, char concatOrCombine);


    void         ProcessPartition       (_String const&, _SimpleList&,  bool, _SimpleList const* = nil, _SimpleList const* = nil, _String const* scope = nil) const;
  
    void         SetTranslationTable    (_DataSet *  newTT );
    void         SetTranslationTable    (_TranslationTable *  newTT );
    _TranslationTable*
    GetTT                   (void) {
        return theTT;
    }
    _Parameter   CheckAlphabetConsistency
    (void);

    void         SetNoSpecies           (unsigned long n) {
        noOfSpecies = n;
    }
    void         ResetIHelper           (void);
  
private:

 
    _SimpleList theMap,
                theFrequencies;         // remapping vector, and the counter of frequencies

    unsigned long
    noOfSpecies;

    _TranslationTable*
    theTT;                  // translation Table, if any

    _List       theNames;               // Names of species
    FILE*       streamThrough;

    _DSHelper*  dsh;
    bool        useHorizontalRep;

};

enum _hy_dataset_filter_ambiguity_resolution {
  kAmbiguityHandlingResolve,
  kAmbiguityHandlingResolveFrequencyAware,
  kAmbiguityHandlingAverageFrequencyAware,
  kAmbiguityHandlingSkip
};

//_________________________________________________________

class _DataSetFilter:public BaseObj {

public:

    _DataSetFilter               (void);
    _DataSetFilter              (_DataSet*, char, _String&);

    virtual                     ~_DataSetFilter (void);

    virtual  BaseRef            toStr           (unsigned long = 0UL);  // convert to string
    virtual  void               toFileStr       (FILE*, unsigned long = 0UL); // convert to string

    virtual  BaseRef            makeDynamic     (void);
    virtual  long               FreeUpMemory    (long);
    virtual  bool               IsNormalFilter  (void) const {
        return true;
    }

    void     CopyFilter         (_DataSetFilter*);
    void     SetFilter          (_DataSet*, char, _SimpleList&, _SimpleList&, bool isFilteredAlready = false);
    void     SetExclusions      (_String*, bool = true);

    _String* GetExclusions      (void) const;

    void     SetMap  (_String&); // used to allow nonsequential maps to tree leaves

    void     SetMap  (_SimpleList& newMap) {
        theNodeMap.Clear();    // used to allow nonsequential maps to tree leaves
        theNodeMap.Duplicate(&newMap);
    }

    unsigned long     NumberSpecies (void) const{
        return theNodeMap.lLength;
    }

    virtual  long GetSiteCount (void) const{
        return theOriginalOrder.lLength;
    }

    virtual  long GetSiteCountInUnits (void) const{
      return theOriginalOrder.lLength / unitLength;
    }

    virtual  long
    GetPatternCount        (void) const{
        return theFrequencies.lLength;
    }

    unsigned long     GetFrequency  (long i) const{
        return theFrequencies.Element(i);
    }

    long     GetUnitLength  (void) const{
        return unitLength;
    }

    virtual  unsigned long    GetDimension (bool correct = true) const;

    long     GetOriginalToShortMap (long i);

    void     ComputePairwiseDifferences (_Matrix&, long, long) const;
    _Matrix* ComputePairwiseDifferences (long, long, _hy_dataset_filter_ambiguity_resolution = kAmbiguityHandlingResolveFrequencyAware) const;

    BaseRefConst  GetMap (void) const {
        return theNodeMap.lLength?&theNodeMap:NULL;
    }

    BaseRefConst  GetDuplicateSiteMap (void) const {
      return duplicateMap.lLength?&duplicateMap:NULL;
    }

    virtual  _String&   operator () (unsigned long site, unsigned long pos);
    // site indexes unique sites

    const _String   RetrieveState (unsigned long site, unsigned long pos) const;
      // site indexes all sites, including duplicates

    virtual  void   RetrieveState (unsigned long site, unsigned long pos, _String&, bool = true) const;
    // site indexes all sites, including duplicates
  
    _TranslationTable const * GetTranslationTable (void) const {
        return theData->theTT;
    }
  
    _String* MakeSiteBuffer (void) const;
    /**
      * allocate a string of the right dimension to pass as an argument to RetrieveState
      * it is the responsibility of the caller to delete the returned string when done
      * @return the buffer string 
     */

    void     FindAllSitesLikeThisOne (long, _SimpleList&);

    _String const&  GenerateConsensusString (_SimpleList* =nil) const;

    void     GrabSite (unsigned long,unsigned long,_String&);
    void     GrabSite (unsigned long,unsigned long,char*);

    virtual  char      GetChar(unsigned long site, unsigned long pos);
    long       SiteFrequency  (unsigned long site);
    bool       HasDeletions   (unsigned long site, _AVLList* = nil);
    long       HasExclusions  (unsigned long site, _SimpleList* theExc, _Parameter *buffer);
    bool       IsConstant     (unsigned long site,  bool relaxedDeletions = true);

    long     Translate2Frequencies (_String const&, _Parameter*, bool) const;
    long     MapStringToCharIndex  (_String&) const;
    //long   Translate2Frequencies (char*,    _Parameter*, bool = true);
    long     Translate2Frequencies (char,     _Parameter*, bool) const;

    _Matrix* HarvestFrequencies (char unit, char atom, bool posSpec, bool = true) const;

    void     Freeze (long);

    void     UnFreeze (long);

    void     MatchStartNEnd (_SimpleList&, _SimpleList&, _SimpleList* = nil) const;

    _String* GetSequenceName(long idx) const {
        return theData->GetSequenceName(theNodeMap.lData[idx]);
    }

    _String* GetSequenceCharacters (long) const;

    _DataSet*                       GetData                     (void) const{
        return theData;
    }
    void                            SetData                     (_DataSet* ds) {
        theData = ds;
    }
    const _String                         ConvertCodeToLetters        (long code, char base) const {
        return theData->theTT->ConvertCodeToLetters(code,base);
    }

    void                            ConvertCodeToLettersBuffered(long code, char base, char *, _AVLListXL* ) const;
    // 20090212: SLKP
    // added this function to cache repeated character code -> string conversions
    // and to skip returning temp objects but simply writing to buffer

    /**
    * Find all unique sequences in the data filter. 
    *
    * \n Usage: FindUniqueSequences (uniqueIndex, instanceCount, true);
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
    * @param mode  Controls is the strings must match exactly (0), exactly + gap (1), via the superset rule (2) or via the partial match rule (3).
                     Nucleotide letters A and - (or ?) (IUPAC code for A or G) will count as mismatched for mode 0 and matched for mode 1.
                     Nucleotide letters A and R (IUPAC code for A or G) will count as mismatched for mode 0 and matched for
                     modes 2 and 3 because R is a superset of A. (note that R matches R in all modes, even though the letter is
                     an ambiguous nucleotide).
                     Nucleotide letters R (A or G) and M (A or C) will match under mode 3 (because they both encode A as an option), 
                     but not under modes 0-2. 
                     Match in mode (0) => match in mode (1) => match in mode (2) => match in mode (3).
                    
    * @return The number of unique sequences. 
    */
    unsigned long                   FindUniqueSequences      (_SimpleList& indices, _SimpleList& map, _SimpleList& counts, short mode = 0) const;


    long                            CorrectCode                 (long code) const;
    virtual  bool                   CompareTwoSites             (unsigned long, unsigned long,unsigned long) const;
    bool                            CompareTwoSitesChar         (unsigned long, unsigned long,unsigned long) const;
    long                            FindSpeciesName             (_List&, _SimpleList&) const;
    _DataSetFilter*                 PairFilter                  (long, long, _DataSetFilter*);
    void                            SetDimensions               ();
    long                            LookupConversion            (char c, _Parameter* receptacle) const;
    void                            SetupConversion             (void);
    bool                            ConfirmConversionCache      (void) const;
    void                            FilterDeletions             (_SimpleList* theExc = nil);
    _Matrix*                        GetFilterCharacters         (bool = false) const;
    _SimpleList*                    CountAndResolve             (long, _Parameter* = nil, bool = false);
    _Matrix*                        PairwiseCompare             (_SimpleList*, _SimpleList*, _List* = nil);

    _List*                          ComputePatternToSiteMap     (void) const;
  
    /**
     * a utility function to return a _List of simplelists (one per unique site pattern) that provides an ordered list of
     * the indices of all sites that have the same pattern in the original alignment
    */
  
    void                            PatternToSiteMapper         (void*, void*, char, long) const;
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

    char const*    GetColumn (long index) const {
        return ((_Site*)(((BaseRef*)theData->lData)[theData->theMap.lData[theMap.lData[index]]]))->sData;
    }

    _SimpleList     conversionCache;

protected:

    unsigned char   unitLength;
    long            dimension;

private:

    void            internalToStr (FILE*,_String&);
    _String*        accessCache;

    long            undimension;

    _DataSet*       theData;
//      _SimpleList     conversionCache;

};



//_________________________________________________________

class _DataSetFilterNumeric:public _DataSetFilter
{

public:

    _DataSetFilterNumeric                       (void) {
    }
    _DataSetFilterNumeric                       (_Matrix*, _List&,_DataSet*,long);

    virtual  bool                                       IsNormalFilter          (void)  const {
        return false;
    }


    virtual  BaseRef                                    makeDynamic             (void);
    virtual  unsigned long                              GetDimension            (bool) const{
        return dimension;
    }

    _Parameter*                 getProbabilityVector    (long,long,long = 0);
    virtual  bool               CompareTwoSites         (unsigned long, unsigned long,unsigned long) const;

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



extern _String  nexusBFBody;
extern _DataSet*lastNexusDataMatrix;

#endif
