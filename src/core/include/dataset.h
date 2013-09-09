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

#ifndef __DATASET__
#define __DATASET__

#include "translationtable.h"
#include "list.h"
#include "avllistx.h"
#include "formula.h"
#include "matrix.h"

class _Site;

// data set file state data struct
struct FileState {

  _TranslationTable *translationTable;
  long curSpecies, totalSpeciesRead, totalSitesRead, totalSpeciesExpected,
      totalSitesExpected, curSite, maxStringLength, pInSrc;
  bool acceptingCommands, allSpeciesDefined, interleaved, autoDetect,
      isSkippingInNEXUS;
  int fileType, baseLength;
  char repeat, skip;

  _String *theSource, *theNamespace;

  _SimpleList rawLinesFormat;
};

// data set file state data struct
struct _DSHelper {

  _SimpleList characterPositions;
  _List incompletePatternStorage;
  _AVLListX *incompletePatterns;

  _DSHelper(void) {
    incompletePatterns = new _AVLListX(&incompletePatternStorage);
    checkPointer(incompletePatterns);
  }
  ~_DSHelper(void) { DeleteObject(incompletePatterns); }
};

class _DataSet : public virtual _List, public virtual _AssociativeList // a complete data set
                 {
public:

  _DataSet(void);
  _DataSet(long);
  _DataSet(FILE *);
  // with estimated number of sites per file
  virtual ~_DataSet(void);

  virtual BaseRef makeDynamic(void);
  
  virtual void Duplicate (BaseRef);

  void AddSite(char);

  void Write2Site(long, char);
  void CheckMapping(long);

  void Finalize(void);
  // remove duplicate data types and compress

  long GetNoTypes(void);
  // return the number of unique sites

  long GetCharDimension(void);
  // return the size of the alphabet space

  long GetFreqType(long);
  // return the frequency of a site

  _Site *GetSite(long index) { return ((_Site **)lData)[theMap.lData[index]]; }

  long ComputeSize(void);
  // compute the size of this object in memory

  void Clear(void);

  virtual char operator()(unsigned long, unsigned long, unsigned int);
  // retrieve element pos of site-th site

  virtual BaseRef toStr(void);
  // convert to string

  virtual void toFileStr(FILE *dest);

  void Compact(long);
  // release string overhead
  void ConvertRepresentations(void);

  _Matrix *HarvestFrequencies(char, char, bool, _SimpleList &, _SimpleList &,
                              bool = true);
  // this function counts observed frequencies of elements in a data set
  // unit is the length of an info unit (nucl - 1, codons - 3)
  // atom is the "minimal" countable element (nucl - 1, codons - 1)
  // posSpec - if position of an atom within an item is to be accounted for
  // segmentation - partition of the underlying DataSet to look at
  // null for segmentation assumes the entire dataset

  void MatchIndices(_Formula &, _SimpleList &, bool, long);
  friend void printFileResults(_DataSet *);
  char InternalStorageMode(void) { return useHorizontalRep; }

  long NoOfSpecies(void) { return noOfSpecies; }

  long NoOfColumns(void) { return theMap.lLength; }
  long NoOfUniqueColumns(void) { return lLength; }
  void AddName(_String &);
  _List &GetNames(void) { return theNames; }
  _SimpleList &GetTheMap(void) { return theMap; }
  void FindAllSitesLikeThisOne(long, _SimpleList &);

  friend class _DataSetFilter;
  friend _DataSet *ReadDataSetFile(FILE *, char, _String *, _String *,
                                   _String *, _TranslationTable *);
  friend long ProcessLine(_String &s, FileState *fs, _DataSet &ds);

  static _DataSet *Concatenate(_List&);
  static _DataSet *Combine(_List&);

  static _TranslationTable *CheckCompatibility(_List &objects,
                                               char concatOrCombine);

  void ProcessPartition(_String &, _SimpleList &, bool, _SimpleList * = nil,
                        _SimpleList * = nil);
  void SetTranslationTable(_DataSet *newTT);
  void SetTranslationTable(_TranslationTable *newTT);
  _TranslationTable *GetTT(void) { return theTT; }
  _Parameter CheckAlphabetConsistency(void);

  void SetNoSpecies(long n) { noOfSpecies = n; }
  void ResetIHelper(void);

private:

  void constructFreq(long *, _Parameter *, char, long, long, int, int, int);

  _SimpleList theMap,
      theFrequencies; // remapping vector, and the counter of frequencies

  unsigned int noOfSpecies;

  _TranslationTable *theTT; // translation Table, if any

  _List theNames; // Names of species
  FILE *streamThrough;

  _DSHelper *dsh;
  bool useHorizontalRep;

};

#endif
