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

#pragma once

#include "list.h"
#include "stdlib.h"
#include "dataset.h"
#include "string_file_wrapper.h"

//_________________________________________________________
enum _hy_dataset_filter_ambiguity_resolution {
  kAmbiguityHandlingResolve,
  kAmbiguityHandlingResolveFrequencyAware,
  kAmbiguityHandlingAverageFrequencyAware,
  kAmbiguityHandlingSkip
};

enum _hy_dataset_filter_unique_match {
  kUniqueMatchExact = 0L,
  kUniqueMatchExactOrGap = 1L,
  kUniqueMatchSuperset = 2L,
  kUniqueMatchPartialMatch = 3L
};

class _DataSetFilter : public BaseObj {

public:
  _DataSetFilter(void);
  _DataSetFilter(_DataSet *, char, _String &);

  virtual ~_DataSetFilter(void);

  virtual BaseRef toStr(unsigned long = 0UL);          // convert to string
  virtual void toFileStr(FILE *, unsigned long = 0UL); // convert to string

  virtual BaseRef makeDynamic(void) const;
  virtual void Duplicate(BaseRefConst);
  virtual bool IsNormalFilter(void) const { return true; }

  void CopyFilter(_DataSetFilter const *);
  void SetFilter(_DataSet const *, unsigned char, _SimpleList &, _SimpleList &,
                 bool isFilteredAlready = false);
  void SetExclusions(_String const&, bool = true);

  _String *GetExclusions(void) const;

  void SetMap(_String const&); // used to allow nonsequential maps to tree leaves

  void SetMap(_SimpleList &newMap) {
    theNodeMap.Clear(); // used to allow nonsequential maps to tree leaves
    theNodeMap.Duplicate(&newMap);
  }

  unsigned long NumberSpecies(void) const { return theNodeMap.lLength; }

  virtual long GetSiteCount(void) const { return theOriginalOrder.lLength; }

  virtual long GetSiteCountInUnits(void) const {
    return theOriginalOrder.lLength / unitLength;
  }

  virtual long GetPatternCount(void) const { return theFrequencies.lLength; }

  unsigned long GetFrequency(long i) const { return theFrequencies.Element(i); }

  unsigned long GetUnitLength(void) const { return unitLength; }

  virtual unsigned long GetDimension(bool correct = true) const;

  long GetOriginalToShortMap(long i);

  void ComputePairwiseDifferences(_Matrix &, long, long) const;
  _Matrix *
  ComputePairwiseDifferences(long, long,
                             _hy_dataset_filter_ambiguity_resolution =
                                 kAmbiguityHandlingResolveFrequencyAware) const;

  BaseRefConst GetMap(void) const {
    return theNodeMap.lLength ? &theNodeMap : NULL;
  }

  BaseRefConst GetDuplicateSiteMap(void) const {
    return duplicateMap.lLength ? &duplicateMap : NULL;
  }

  virtual _String &operator()(unsigned long site, unsigned long pos);
  // site indexes unique sites

  const _String RetrieveState(unsigned long site, unsigned long pos) const;
  // site indexes all sites, including duplicates

  virtual void RetrieveState(unsigned long site, unsigned long pos, _String &,
                             bool = true) const;
  // site indexes all sites, including duplicates

  _TranslationTable const *GetTranslationTable(void) const {
    return theData->theTT;
  }

  _String *MakeSiteBuffer(void) const;
  /**
   * allocate a string of the right dimension to pass as an argument to
   * RetrieveState it is the responsibility of the caller to delete the returned
   * string when done
   * @return the buffer string
   */

  _String const GenerateConsensusString(_SimpleList * = nil) const;


  long site_frequency (unsigned long site) const { return theFrequencies.get (site); }
  bool HasDeletions(unsigned long site, _AVLList * = nil) const;
  long HasExclusions(unsigned long site, _SimpleList *theExc, hyFloat *buffer) const;

  long Translate2Frequencies(_String const &, hyFloat *, bool, bool = true) const;
  long MapStringToCharIndex(_String &) const;
  // long   Translate2Frequencies (char*,    hyFloat*, bool = true);
  long Translate2Frequencies(char, hyFloat *, bool) const;

  _Matrix *HarvestFrequencies(char unit, char atom, bool posSpec,
                              bool = true) const;


  void MatchStartNEnd(_SimpleList &, _SimpleList &, _SimpleList * = nil) const;

  _String *GetSequenceName(long idx) const {
    return theData->GetSequenceName(theNodeMap.get(idx));
  }

  _String *GetSequenceCharacters(long) const;

  _DataSet *GetData(void) const { return theData; }
  void SetData(_DataSet *ds) { theData = ds; }
  const _String ConvertCodeToLetters(long code, char base) const {
    return theData->theTT->ConvertCodeToLetters(code, base);
  }

  void ConvertCodeToLettersBuffered(long code, unsigned char base, _String&,
                                    _AVLListXL *) const;
  // 20090212: SLKP
  // added this function to cache repeated character code -> string conversions
  // and to skip returning temp objects but simply writing to buffer

  /**
  * Find all unique sequences in the data filter.
  *
  * \n Usage: FindUniqueSequences (uniqueIndex, instanceCount, true);
  * @author SLKP
  * @param indices For each sequence - the list of indices corresponding to the
  unique strings For example, if sequence 1 == sequence 3 and sequence 4 ==
  sequence 5 this list will contain 0,1,3
  * @param map The index of the unique string to which the current string is
  mapped For example, if sequence 1 == sequence 3 and sequence 4 == sequence 5
  this list will contain 0,1,0,2,2
  * @param counts  The number of copies for each unique string found
                   For example, if sequence 1 == sequence 3 and sequence 4 ==
  sequence 5 this list will contain 2,1,2
  * @param mode  Controls is the strings must match exactly (0), exactly + gap
  (1), via the superset rule (2) or via the partial match rule (3). Nucleotide
  letters A and - (or ?) (IUPAC code for A or G) will count as mismatched for
  mode 0 and matched for mode 1. Nucleotide letters A and R (IUPAC code for A or
  G) will count as mismatched for mode 0 and matched for modes 2 and 3 because R
  is a superset of A. (note that R matches R in all modes, even though the
  letter is an ambiguous nucleotide). Nucleotide letters R (A or G) and M (A or
  C) will match under mode 3 (because they both encode A as an option), but not
  under modes 0-2. Match in mode (0) => match in mode (1) => match in mode (2)
  => match in mode (3).

  * @return The number of unique sequences.
  */
  unsigned long FindUniqueSequences(_SimpleList &indices, _SimpleList &map,
                                    _SimpleList &counts,  _hy_dataset_filter_unique_match mode = kUniqueMatchExact) const;

  long CorrectCode(long code) const;
  virtual bool CompareTwoSites(unsigned long, unsigned long,
                               unsigned long) const;

  template <typename FUNC> void CompareTwoSitesCallback(unsigned long site1, unsigned long site2,
                                                        FUNC&& callback) const {
                                                            
                                                    
        switch (unitLength) {
            case 3: {
                site1 = (site1 << 1) + site1;
                site2 = (site2 << 1) + site2;
                
                const char * site1c[3] = {GetColumn(site1),GetColumn(site1+1),GetColumn(site1+2)},
                           * site2c[3] = {GetColumn(site2),GetColumn(site2+1),GetColumn(site2+2)};
                
                theNodeMap.Each ([site1c,site2c,&callback] (long ci, unsigned long i) -> void {
                    if (site1c[0][ci] != site2c[0][ci] || site1c[1][ci] != site2c[1][ci] || site1c[2][ci] != site2c[2][ci]) {
                        callback (ci, i);
                    }
                });
             }
             break;
                
            case 2: {
                site1 = (site1 << 1);
                site2 = (site2 << 1);
                
                const char * site1c[2] = {GetColumn(site1),GetColumn(site1+1)},
                           * site2c[2] = {GetColumn(site2),GetColumn(site2+1)};
                
                theNodeMap.Each ([site1c,site2c,&callback] (long ci, unsigned long i) -> void {
                    if (site1c[0][ci] != site2c[0][ci] || site1c[1][ci] != site2c[1][ci]) {
                        callback (ci, i);
                    }
                });
             }
            break;

            case 1: {
                
                const char * site1c = GetColumn(site1),
                           * site2c = GetColumn(site2);
                
                theNodeMap.Each ([site1c,site2c,&callback] (long ci, unsigned long i) -> void {
                    if (site1c[ci] != site2c[ci]) {
                        callback (ci, i);
                    }
                });
                break;
            }
                
            default: {
                    site1 *= unitLength;
                    site2 *= unitLength;
                
                    const char **site1c = (const char **)alloca (sizeof (char*) * unitLength),
                               **site2c = (const char **)alloca (sizeof (char*) * unitLength);
                
                    for (int i = 0; i < unitLength; i++) {
                        site1c[i] = GetColumn (site1 + i);
                        site2c[i] = GetColumn (site2 + i);
                    }
                    
                    theNodeMap.Each ([site1c,site2c,this,&callback] (long ci, unsigned long ) -> void {
                            for (int i = 0; i < this->unitLength; i++) {
                                if (site1c[i][ci] != site2c[i][ci]) {
                                    callback (ci, i);
                                    break;
                                }
                            }
                    });
                    break;
            }
        }
  }
    
  long FindSpeciesName(_List &, _SimpleList &) const;
  _DataSetFilter *PairFilter(long, long, _DataSetFilter *);
  void SetDimensions();
  long LookupConversion(char c, hyFloat *receptacle) const;
  void SetupConversion(void);
  bool ConfirmConversionCache(void) const;
  void FilterDeletions(_SimpleList *theExc = nil);
  _Matrix *GetFilterCharacters(bool = false) const;

  _List *ComputePatternToSiteMap(void) const;

  /**
   * a utility function to return a _List of simplelists (one per unique site
   * pattern) that provides an ordered list of the indices of all sites that
   * have the same pattern in the original alignment
   */

  template<typename SOURCE_TYPE, typename TARGET_TYPE> void PatternToSiteMapper(SOURCE_TYPE const* source, TARGET_TYPE  * target,  long padup, TARGET_TYPE filler) const {
    
  
    duplicateMap.Each ([source, target] (long value, unsigned long site) -> void {
      target [site] = source[value];
    });
    
    if (padup > 0) {
      for (long site = duplicateMap.countitems(); site < padup; site++) {
        target[site] = filler;
      }
    }
  }
  
  /*
      20180920: SLKP changed to templates
      20090325: SLKP
   
      a function that takes per pattern values (source, argument 1)
      and maps them onto sites into target (argument 2)
      the third argument is 
      0 to treat the pointers as hyFloat*
      1 to treat them as long*
      2 and to treat them as hyFloat* and long*, respetively
      20090929: SLKP
      the third argument is used to speficy a padding-size,
          all values from the filter size up to that value are set to 1 (for
     mode 0) and 0 (for mode 1) this is needed to handle uneven data filters in
     SITE_LIKELIHOOD constructs
   */

  _SimpleList theFrequencies, theNodeMap, theMap, theOriginalOrder,
      theExclusions, duplicateMap;

  char const *GetColumn(long index) const {
    return (const char *)(*(_Site *)((
        (BaseRef *)
            theData->list_data)[theData->theMap.list_data[theMap.list_data[index]]]));
  }

  _SimpleList conversionCache;

protected:
  unsigned char unitLength;
  long dimension;

private:
  void internalToStr(FILE *, _StringBuffer *);
  
  
   inline void retrieve_individual_site_from_raw_coordinates (_String & store, unsigned long site, unsigned long sequence) const {
        if (unitLength==1UL) {
          store.set_char (0, (((_String**)theData->list_data)[theData->theMap.list_data[theMap.list_data[site]]])->char_at (sequence));
        } else {
          site*=unitLength;
          for (unsigned long k = 0UL; k<unitLength; k++) {
            store.set_char (k, ((_String**)theData->list_data)[theData->theMap.list_data[theMap.list_data[site++]]]->char_at(sequence));
          }
        }
   }
    inline char direct_index_character (unsigned long site, unsigned long sequence) const {
        return (((_String**)theData->list_data)[theData->theMap.list_data[theMap.list_data[site]]])->char_at(sequence);
    }
  
   _String *accessCache;

  long undimension;

  _DataSet *theData;
  //      _SimpleList     conversionCache;
};


