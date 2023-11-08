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

#include "global_things.h"
#include "list.h"
#include "stdlib.h"
#include "translation_table.h"
#include "site.h"
#include "function_templates.h"

using namespace hy_global;

#define HYPHY_SITE_DEFAULT_BUFFER_SIZE 512
#define DATA_SET_SWITCH_THRESHOLD 100000

// data set file state data struct
struct _DSHelper {

    _SimpleList characterPositions;
    _List       incompletePatternStorage;
    _AVLListX*  incompletePatterns;

    _DSHelper(void) {
        incompletePatterns = new _AVLListX (&incompletePatternStorage);
    }
    ~_DSHelper(void) {
        DeleteObject (incompletePatterns);
    }
};


class _DataSet : public _List // a complete data set
{
public:
  _DataSet(void);
  _DataSet(long);
  _DataSet(hyFile *);
  // with estimated number of sites per file
  virtual ~_DataSet(void);

  virtual BaseRef makeDynamic(void) const;

  void AddSite(char);

  void Write2Site(long, char, char skip_char = '?');

  void Finalize(void);
  // remove duplicate data types and compress

  long GetNoTypes(void) const;
  // return the number of unique sites

  unsigned long GetCharDimension(void) const;
  // return the size of the alphabet space

  unsigned long GetFreqType(long) const;
  // return the frequency of a site

  _Site *GetSite(long index) const {
    return ((_Site **)list_data)[theMap.list_data[index]];
  }

  long ComputeSize(void);
  // compute the size of this object in memory

  void Clear(bool = true);

  virtual char operator()(unsigned long, unsigned long, unsigned int) const;
  // retrieve element pos of site-th site

  virtual BaseRef toStr(unsigned long = 0UL);
  // convert to string

  virtual void toFileStr(hyFile *dest, unsigned long = 0UL);

  void Compact(long);
  // release string overhead
  void ConvertRepresentations(void);

  _Matrix *HarvestFrequencies(unsigned char, unsigned char, bool, _SimpleList &,
                              _SimpleList &, bool = true) const;
  // this function counts observed frequencies of elements in a data set
  // unit is the length of an info unit (nucl - 1, codons - 3)
  // atom is the "minimal" countable element (nucl - 1, codons - 1)
  // posSpec - if position of an atom within an item is to be accounted for
  // segmentation - partition of the underlying DataSet to look at
  // null for segmentation assumes the entire dataset

  void MatchIndices(_Formula &, _SimpleList &, bool, long,
                    _String const * = nil) const;
  friend void printFileResults(_DataSet *);
  char InternalStorageMode(void) const { return useHorizontalRep; }
  unsigned long NoOfSpecies(void) const { return noOfSpecies; }
  unsigned long NoOfColumns(void) const { return theMap.lLength; }
  unsigned long NoOfUniqueColumns(void) const { return lLength; }
  void AddName(_String const &);
  void InsertName(_String const &name, long where);

  _String *GetSequenceName(unsigned long i) const {
    return (_String *)theNames.GetItem(i);
  }

  _List const &GetNames(void) const { return theNames; }

  void ClearNames(void) { theNames.Clear(); }

  _String *GetSequenceCharacters(long seqID) const;

  bool SetSequenceName(long index, _String *new_name) {
    if (index >= 0L && index < (long)theNames.lLength) {
      theNames.Replace(index, new_name, false);
      return true;
    }
    return false;
  }

  void SetNames(_List const &copy_from) {
    theNames.Clear();
    theNames << copy_from;
  }

  _SimpleList &GetTheMap(void) { return theMap; }

  _SimpleList const &DuplicateMap(void) const { return theMap; }

  friend class _DataSetFilter;
  friend _DataSet *ReadDataSetFile(hyFile *, char, _String *, _String *,
                                   _String *, _TranslationTable *,
                                   _ExecutionList *);
  friend long ProcessLine(_String &s, FileState *fs, _DataSet &ds);

  static _DataSet *Concatenate(const _SimpleList &);
  static _DataSet *Combine(const _SimpleList &);

  static _TranslationTable *CheckCompatibility(_SimpleList const &ref,
                                               char concatOrCombine);

  void ProcessPartition(_String const &, _SimpleList &, bool, int unit_length,
                        _SimpleList const * = nil, _SimpleList const * = nil,
                        _String const *scope = nil) const;

  void SetTranslationTable(_DataSet *newTT);
  void SetTranslationTable(_TranslationTable *newTT);
  _TranslationTable const *GetTT(void) const { return theTT; }
  hyFloat CheckAlphabetConsistency(void);

  void SetNoSpecies(unsigned long n) { noOfSpecies = n; }
  void ResetIHelper(void);

private:
  _SimpleList theMap,
      theFrequencies; // remapping vector, and the counter of frequencies

  unsigned long noOfSpecies;

  _TranslationTable *theTT; // translation Table, if any

  _List theNames; // Names of species
  hyFile *streamThrough;

  _DSHelper *dsh;
  bool useHorizontalRep;
};

void ReadNextLine(hyFile *fp, _StringBuffer *s, FileState *fs, bool append = false,
                  bool upCase = true);

_DataSet *ReadDataSetFile(hyFile *, char = 0, _String * = nil, _String * = nil,
                          _String * = nil,
                          _TranslationTable * = &hy_default_translation_table,
                          _ExecutionList *target = nil);


bool StoreADataSet(_DataSet *, _String *);
void    ReadNexusFile               (FileState& fState, hyFile*f, _DataSet& result);


extern _StringBuffer nexusBFBody;
extern _DataSet *lastNexusDataMatrix;

