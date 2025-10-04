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

#include <ctype.h>

#include "batchlan.h"
#include "dataset.h"
#include "global_object_lists.h"
#include "site.h"
#include "translation_table.h"

using namespace hyphy_global_objects;

//_________________________________________________________

/**
 * @brief Converts a character to uppercase.
 * @param in The character to convert.
 * @return The uppercase version of the character, or the original character if
 * it's not a lowercase letter.
 *
 * @example
 * @code
 *   char lower = 'a';
 *   char upper = _uppercase_char(lower);
 *   // upper is now 'A'
 * @endcode
 */
inline char _uppercase_char(char in) {
  if (in >= 'a' && in <= 'z') {
    return (char)(in - 'a' + 'A');
  }
  return in;
}

//_________________________________________________________

#define DATA_SET_SWITCH_THRESHOLD 100000

/**
 * @brief Default constructor for the _DataSet class.
 *
 * Initializes an empty dataset with a default translation table.
 */
_DataSet::_DataSet(void)
    : theTT(&hy_default_translation_table), streamThrough(nullptr),
      useHorizontalRep(false) {}

/**
 * @brief Constructs a _DataSet with a pre-allocated capacity.
 * @param l The initial capacity of the dataset.
 *
 * This constructor is used when the number of sites is roughly known
 * beforehand, allowing for pre-allocation of memory to improve performance.
 */
_DataSet::_DataSet(long l)
    : _List((unsigned long)l),
      theFrequencies(
          (unsigned long)l) // with estimated number of sites per file
{
  streamThrough = nullptr;
  theTT = &hy_default_translation_table;
  useHorizontalRep = false;
}

//_______________________________________________________________________

/**
 * @brief Constructs a _DataSet that streams data to a file.
 * @param f A pointer to a hyFile object to which the data will be written.
 *
 * This constructor is used for creating large datasets that are written
 * directly to a file to avoid storing them in memory.
 */
_DataSet::_DataSet(hyFile *f) {
  useHorizontalRep = false;
  theTT = &hy_default_translation_table;
  streamThrough = f;
  theMap << 0; // current sequence
  theMap << 0; // current site
  theMap << 0; // total sites
}

//_______________________________________________________________________

/**
 * @brief Destructor for the _DataSet class.
 *
 * Frees the memory allocated for the translation table if it's not the default
 * one.
 */
_DataSet::~_DataSet(void) {
  if (theTT != &hy_default_translation_table) {
    DeleteObject(theTT);
  }
}

//_______________________________________________________________________

/**
 * @brief Clears the contents of the dataset.
 *
 * Resets the dataset to an empty state, clearing all sequences, names,
 * frequencies, and the translation table.
 */
void _DataSet::Clear(bool) {
  _List::Clear();
  theMap.Clear();
  theFrequencies.Clear();
  theNames.Clear();
  if (theTT != &hy_default_translation_table) {
    DeleteObject(theTT);
    theTT = &hy_default_translation_table;
  }
  noOfSpecies = 0;
  useHorizontalRep = false;
}

//_______________________________________________________________________

/**
 * @brief Creates a dynamic copy of the dataset.
 * @return A new _DataSet object that is a copy of the current one.
 *
 * This function creates a new _DataSet instance and copies the essential
 * properties of the current dataset, such as the map, frequencies, names,
 * and translation table.
 * @note The caller is responsible for managing the memory of the returned
 * object.
 */
BaseRef _DataSet::makeDynamic(void) const {
  _DataSet *r = new _DataSet;
  r->theMap.Duplicate(&theMap);
  r->theFrequencies.Duplicate(&theFrequencies);
  if (theTT != &hy_default_translation_table) {
    r->theTT->AddAReference();
  }
  r->theNames.Duplicate(&theNames);
  r->streamThrough = streamThrough;

  r->useHorizontalRep = false;
  return r;
}

//_______________________________________________________________________

//_______________________________________________________________________

/**
 * @brief Converts the internal data storage to a horizontal representation.
 *
 * This function changes the data storage from a "vertical" representation (a
 * list of sites, where each site is a column of characters) to a "horizontal"
 * representation (a list of sequences, where each sequence is a string). This
 * can be more efficient for certain operations, especially when the number of
 * sites is large.
 */
void _DataSet::ConvertToHorizontalRepresentation(void) {
  if (useHorizontalRep == false) {
    _List horStrings;

    if (lLength == 0) {
      AppendNewInstance(new _StringBuffer(128UL));
    } else {
      _Site *aSite = (_Site *)list_data[0];

      // Create a new string buffer for each sequence
      for (unsigned long str = 0; str < aSite->length(); str++) {
        horStrings < new _StringBuffer(DATA_SET_SWITCH_THRESHOLD);
      }

      // Transpose the data from sites to sequences
      for (unsigned long s = 0; s < lLength; s++) {
        _Site *site_record = (_Site *)GetItem(s);

        if (site_record->length() > horStrings.lLength ||
            site_record->GetRefNo() != -1) {
          HandleApplicationError("Internal error: Inconsistent data in "
                                 "_DataSet::ConvertToHorizontalRepresentation.",
                                 true);
          return;
        }

        for (unsigned long s2 = 0; s2 < site_record->length(); s2++) {
          (*(_StringBuffer *)horStrings.list_data[s2])
              << site_record->get_char(s2);
        }
      }

      _List::Clear();
      theFrequencies.Clear();
      // Replace the old vertical data with the new horizontal data
      (*this) << horStrings;
    }
    useHorizontalRep = true;
  }
}

//_______________________________________________________________________

/**
 * @brief Adds a character to a new site in the dataset.
 * @param c The character to add.
 *
 * This function's behavior depends on the dataset's configuration:
 * - If the dataset is being streamed to a file (`streamThrough` is not null),
 * the character is written to the file.
 * - If the dataset is in memory and using the vertical representation, a new
 * site is created with the given character.
 * - If the dataset has been converted to the horizontal representation, the
 * character is appended to the first sequence.
 * - If the number of sites exceeds `DATA_SET_SWITCH_THRESHOLD`, the dataset is
 * converted to the horizontal representation.
 */
void _DataSet::AddSite(char c) {
  if (streamThrough) {
    if (theMap.get(0) == 0) {
      if (theMap.get(1) == 0) {
        if (theNames.lLength) {
          streamThrough->puts(((_String *)theNames(0))->get_str());
          streamThrough->fputc('\n');
        } else {
          streamThrough->puts(">Sequence 1");
        }
        AppendNewInstance(new _String(kEmptyString));
      }

      theMap.list_data[1]++;
      theMap.list_data[2]++;
      streamThrough->fputc(c);
    } else {
      HandleApplicationError(
          "Cannot add sites to a file-based dataset with multiple sequences.",
          true);
    }
  } else {
    if (useHorizontalRep == false) {
      if (lLength < DATA_SET_SWITCH_THRESHOLD) {
        theFrequencies << 1L;
        AppendNewInstance(new _Site(c));
        return;
      } else {
        ConvertToHorizontalRepresentation();
      }
    }

    (*GetSite(0)) << c;
  }
}
//_______________________________________________________________________

/**
 * @brief Writes a character to a specific position in the dataset.
 *
 * @param index The site index (column) to write to.
 * @param c The character to write.
 * @param skip_char The character to use for padding if a new sequence is
 * shorter than existing ones.
 *
 * This function writes a character `c` to the specified `index`. The exact
 * behavior depends on the internal state of the `_DataSet` object:
 *
 * - **File-based streaming (`streamThrough` is active):** The function handles
 * writing sequences one after another to the file. It ensures that sequences
 * have equal lengths.
 * - **Horizontal in-memory representation (`useHorizontalRep` is true):** The
 * function appends the character to the appropriate sequence string. It can
 * also handle creating new sequences.
 * - **Vertical in-memory representation (`useHorizontalRep` is false):** The
 * function appends the character to the `_Site` object at the given `index`.
 */
void _DataSet::Write2Site(long index, char c, char skip_char) {
  if (streamThrough) {
    if (index == 0) {
      if (theMap.get(2) == theMap.get(1)) {
        theMap.list_data[0]++;

        if ((long)theNames.lLength > theMap.list_data[0]) {
          streamThrough->puts("\n>");
          streamThrough->puts(
              ((_String *)theNames(theMap.list_data[0]))->get_str());
          streamThrough->fputc('\n');
        } else {
          streamThrough->puts("\n>");
          char buffer[64];
          snprintf(buffer, 64, "%ld", theMap.list_data[0] + 1);
          streamThrough->puts(buffer);
          streamThrough->fputc('\n');
        }

        theMap.list_data[1] = 0;
      } else {
        HandleApplicationError("Cannot write sequences of unequal lengths to a "
                               "file-based dataset.");
        return;
      }
    } else if (index != theMap.list_data[1]) {
      HandleApplicationError(
          "Sites must be written consecutively to a file-based dataset.");
      return;
    }

    theMap.list_data[1]++;
    streamThrough->fputc(c);
  } else {
    if (useHorizontalRep) {
      long currentWritten = ((_String *)list_data[0])->length();

      if (index >= currentWritten) {
        // When we enter here, this means that the current string is longer
        // that all other strings; they need to be padded
        for (unsigned long s = 0; s < lLength - 1; s++) {
          (*GetSite(s)) << skip_char;
        }

        (*GetSite(lLength - 1)) << c;
      } else {
        if (index == 0) {
          _StringBuffer *newString = new _StringBuffer(currentWritten);
          (*newString) << c;
          (*this) < newString;
        } else {
          // check to see if the last sequence has the right length

          _StringBuffer *last_record = (_StringBuffer *)GetItem(lLength - 1);

          if ((long)last_record->length() == index) {
            *last_record << c;
          } else {
            unsigned long s = 1;
            for (; s < lLength; s++) {
              _Site *site_s = GetSite(s);
              if ((long)site_s->length() == index) {
                (*site_s) << c;
                break;
              }
            }
            if (s == lLength) {
              HandleApplicationError(
                  "Internal error in 'Write2Site': no appropriate string to "
                  "write to (compact representation).");
              return;
            }
          }
        }
      }

    } else {
      if (index >= (long)lLength) {
        HandleApplicationError(
            _String(
                "Internal error in 'Write2Site': index is out of bounds. ") &
            index & '/' & (long)lLength);
        return;
      }
      (*(_Site *)GetItem(index)) << c;
    }
  }
}

//_______________________________________________________________________

unsigned long _DataSet::GetCharDimension(
    void) const { // return the size of the alphabet space
  return theTT->baseLength;
}

//_______________________________________________________________________

long _DataSet::GetNoTypes(void) const // return the number of unique columns
{
  return theMap.countitems();
}
//_______________________________________________________________________

unsigned long
_DataSet::GetFreqType(long index) const { // return the frequency of a site
  return theFrequencies(theMap(index));
}
//_______________________________________________________________________

void _DataSet::SetTranslationTable(_DataSet *newTT) {
  if (theTT && (theTT != &hy_default_translation_table)) {
    DeleteObject(theTT);
  }
  theTT = (_TranslationTable *)newTT->theTT->makeDynamic();
}

//_______________________________________________________________________

void _DataSet::SetTranslationTable(_TranslationTable *newTT) {
  if (theTT && (theTT != &hy_default_translation_table)) {
    DeleteObject(theTT);
  }
  theTT = (_TranslationTable *)newTT->makeDynamic();
}
//_______________________________________________________________________
/**
 * @brief Finalizes the dataset after all data has been added.
 *
 * This function is called after all sequences and sites have been added to the
 * dataset. Its main purpose is to process the raw data to identify unique site
 * patterns, calculate their frequencies, and compact the dataset to save
 * memory.
 *
 * The behavior depends on the internal state of the dataset:
 * - **File-based streaming:** The output file stream is closed.
 * - **Horizontal in-memory representation:** The function transposes the data
 * back to a vertical representation (list of sites), finds unique site patterns
 * using a checksum-based approach, and builds a map from original site indices
 * to unique pattern indices.
 * - **Vertical in-memory representation:** The function iterates through the
 * sites, finds duplicates using an AVL tree, and updates the site frequencies
 * and mapping accordingly.
 */
void _DataSet::Finalize(void) {
  if (streamThrough) {
    streamThrough->close();
    delete (streamThrough);
    theMap.Clear();
  } else {
    if (useHorizontalRep) {
      bool good = true;

      const unsigned long ref_length = ((_String *)list_data[0])->length();
      _Site **sequence_records = (_Site **)list_data;

      for (long s = 0; s < (long)lLength && good; s++) {
        good = good && ref_length == sequence_records[s]->length();
      }

      if (!good) {
        Clear();
        HandleApplicationError("Internal Error in _DataSet::Finalize. Unequal "
                               "sequence lengths in compact representation",
                               true);
        return;
      }

      _List uniquePats;

      _SimpleList _checkSumBins;
      _AVLListXL checkSumBins(&_checkSumBins);

      long siteCounter = ((_String *)list_data[0])->length();

      _String site_holder(lLength, nullptr);

      auto insert_new_pattern = [&](long index) -> void {
        _SimpleList *indices = (_SimpleList *)checkSumBins.GetXtra(index);
        _Site *tC = new _Site(lLength, -1);
        (*indices) << uniquePats.lLength;
        (*tC) << site_holder;
        uniquePats << tC;
        theMap << theFrequencies.lLength;
        theFrequencies << 1;
        DeleteObject(tC);
      };

      _SimpleList *sites_with_same_checksum = new _SimpleList;

      for (long i1 = 0L; i1 < siteCounter; i1++) {
        for (unsigned long i2 = 0L; i2 < lLength; i2++) {
          site_holder.set_char_no_check(i2, sequence_records[i2]->get_char(i1));
        }

        long patternChecksum = site_holder.Adler32();

        long checkSumExists =
            checkSumBins.Insert((BaseRef)patternChecksum,
                                (long)sites_with_same_checksum, false, false);
        if (checkSumExists >= 0) {
          insert_new_pattern(checkSumExists);
          sites_with_same_checksum = new _SimpleList;
        } else {

          checkSumExists = -checkSumExists - 1;

          _SimpleList *indices =
              (_SimpleList *)checkSumBins.GetXtra(checkSumExists);
          unsigned long i = 0;
          for (; i < indices->lLength; i++) {
            long pattern_index = indices->get(i);
            _Site *comp = (_Site *)uniquePats.GetItem(pattern_index);
            if (comp->Equal(site_holder)) {
              theMap << pattern_index;
              theFrequencies.list_data[pattern_index]++;
              break;
            }
          }
          if (i == indices->lLength) {
            insert_new_pattern(checkSumExists);
          }
        }
      }
      DeleteObject(sites_with_same_checksum);
      _List::Clear();
      _List::Duplicate(&uniquePats);
    } else {
      long j, k;

      _Site *tC;
      {
        _List dups;
        _AVLListX dupsAVL(&dups);

        for (unsigned long i1 = 0; i1 < lLength; i1++) {
          tC = (_Site *)list_data[i1];
          long ff = dupsAVL.Find(tC);
          if (ff < 0) {
            dupsAVL.Insert(tC, i1);
          } else {
            ff = dupsAVL.GetXtra(ff);
            tC->Clear();
            tC->SetRefNo(ff);
            theFrequencies.list_data[ff]++;
          }
        }
        dupsAVL.Clear(false);
      }

      _SimpleList refs(lLength), toDelete(lLength);
      j = 0;

      for (unsigned long i1 = 0; i1 < lLength; i1++) {
        tC = (_Site *)(*(_List *)this)(i1);
        k = tC->GetRefNo();
        if (k == -1) {
          refs << j++;
        } else {
          toDelete << i1;
          refs << -1;
        }
      }

      for (unsigned long i2 = 0; i2 < lLength; i2++) {
        tC = (_Site *)(*(_List *)this)(i2);
        k = tC->GetRefNo();
        if (k >= 0) {
          j = refs.list_data[k];
          if (j < 0) {
            HandleApplicationError(kErrorStringDatasetRefIndexError);
          } else {
            refs.list_data[i2] = j;
          }
        }
      }

      theMap.Clear();
      theMap.Duplicate(&refs);
      DeleteList(toDelete);
      theFrequencies.DeleteList(toDelete);

      for (unsigned long i3 = 0; i3 < lLength; i3++) {
        tC = (_Site *)GetItem(i3);
        tC->TrimSpace();
        tC->SetRefNo(0);
      }
    }
  }
}
//_______________________________________________________________________
/**
 * @brief Marks a site as a duplicate of another.
 * @param index The index of the site to mark as a duplicate.
 *
 * This function is part of the vertical-to-horizontal conversion process.
 * It handles the case where a site is found to be a duplicate of a previously
 * encountered site. It updates the reference of the current site to point to
 * the original one and increments the frequency count of the original site.
 */
void _DataSet::Compact(long index) {
  if (useHorizontalRep) {
    HandleApplicationError("Internal Error: _DataSet::Compact called on a "
                           "dataset already using a Compact representation",
                           true);
    return;
  }
  _Site *tC = (_Site *)GetItem(index);
  if (tC->GetRefNo() != -1) {
    // take care of double referencing
    _Site *tCC = tC;
    long lastRef, count = 0L;
    do {
      lastRef = tCC->GetRefNo();
      count++;
      tCC = (_Site *)GetItem(tCC->GetRefNo());
    } while (tCC->GetRefNo() != -1);

    if (count > 1L) {
      theFrequencies[lastRef]++;
    }

    tC->SetRefNo(lastRef);
  }
}

//_______________________________________________________________________
/**
 * @brief Provides access to a character in the dataset using function-call
 * syntax.
 * @param site The site index (column).
 * @param pos The sequence index (row).
 * @param (unused) An unused third parameter.
 * @return The character at the specified `(site, pos)` coordinate.
 *
 * @example
 * @code
 *   char character = myDataSet(5, 2); // Gets the character from site 5,
 * sequence 2.
 * @endcode
 */
inline char _DataSet::operator()(unsigned long site, unsigned long pos,
                                 unsigned int) const {
  return (((_String **)list_data)[theMap.list_data[site]])->get_char(pos);
}

//_________________________________________________________
/**
 * @brief Computes the approximate memory size of the dataset in bytes.
 * @return The estimated size of the dataset object and its associated data.
 */
long _DataSet::ComputeSize(void) {
  long res = sizeof(_DataSet);

  res += (theMap.lLength + lLength + theFrequencies.lLength) * sizeof(long);
  res += lLength * sizeof(_Site);

  for (unsigned long i = 0; i < lLength; i++) {
    res += ((_Site *)GetItem(i))->length();
  }
  return res;
}

//_________________________________________________________
/**
 * @brief Checks the consistency of the dataset against its translation table.
 * @return A float representing the proportion of characters in the dataset that
 * are valid according to the current translation table (excluding gaps). A
 * value of 1.0 indicates perfect consistency.
 *
 * This function iterates through all the data to verify that the characters
 * are defined in the current alphabet. It's useful for validating data
 * integrity.
 */
hyFloat _DataSet::CheckAlphabetConsistency(void) {
  long validCharacterCount = 0L, gaps = 0L, total = 0L;

  bool checks[256] = {false};

  char gapChar = theTT->GetGapChar();

  _String baseSymbols;

  if (theTT->baseSet.length()) {
    baseSymbols = theTT->baseSet;
  } else if (theTT->baseLength == 4) {
    baseSymbols = "ACGUT";
  } else if (theTT->baseLength == 20) {
    baseSymbols =
        _TranslationTable::GetDefaultTable(HY_TRANSLATION_TABLE_PROTEIN);
  } else {
    baseSymbols =
        _TranslationTable::GetDefaultTable(HY_TRANSLATION_TABLE_BINARY);
  }

  for (long i = 0; i < (long)baseSymbols.length(); i++) {
    checks[(unsigned char)baseSymbols.get_char(i)] = true;
  }

  for (unsigned long i = 0; i < lLength; i++) {
    _String *thisColumn = (_String *)GetItem(i);
    long w = theFrequencies.get(i);
    for (unsigned long j = 0; j < thisColumn->length(); j++)
      if (checks[(unsigned char)thisColumn->get_char(j)]) {
        validCharacterCount += w;
      } else if (gapChar == thisColumn->get_char(j)) {
        gaps += w;
      }

    total += w * thisColumn->length();
  }

  return (hyFloat)validCharacterCount / (total - gaps + 1.);
}

//___________________________________________________

/**
 * @brief Generates a string summary of the dataset.
 * @param (unused) An unused parameter.
 * @return A `_StringBuffer` object containing a summary of the dataset,
 * including the number of species and sites.
 */
BaseRef _DataSet::toStr(unsigned long) {
  _StringBuffer *s = new _StringBuffer(NoOfSpecies() * 30);

  (*s) << _String((long)NoOfSpecies()) << " species:";

  (*s) << (_String *)theNames.toStr();

  (*s) << ";\nTotal Sites:" << _String((long)GetNoTypes())
       << ";\nDistinct Sites:" << _String((long)theFrequencies.lLength);

  return s;
}

//___________________________________________________

/**
 * @brief Writes a string summary of the dataset to a file.
 * @param dest A pointer to the `hyFile` to write to.
 * @param padding A parameter for formatting the output of the sequence names.
 */
void _DataSet::toFileStr(hyFile *dest, unsigned long padding) {
  char buffer[512];
  snprintf(buffer, 512, "%ld species: ", NoOfSpecies());
  dest->puts(buffer);

  theNames.toFileStr(dest, padding);
  snprintf(buffer, 512, ";\nTotal Sites: %ld", GetNoTypes());
  dest->puts(buffer);

  snprintf(buffer, 512, ";\nDistinct Sites: %ld", theFrequencies.lLength);
  dest->puts(buffer);
}

//_________________________________________________________

/**
 * @brief Adds a sequence name to the dataset.
 * @param s The name of the sequence to add.
 *
 * The name is trimmed of trailing whitespace before being added.
 */
void _DataSet::AddName(_String const &s) {
  theNames.AppendNewInstance(
      new _String(s, 0, s.FirstNonSpaceIndex(0, -1, kStringDirectionBackward)));
}

//_________________________________________________________

/**
 * @brief Inserts a sequence name at a specific position.
 * @param name The name of the sequence to insert.
 * @param where The 0-based index at which to insert the name.
 */
void _DataSet::InsertName(_String const &name, long where) {
  theNames.InsertElement(new _String(name), where, false);
}

//_________________________________________________________

/**
 * @brief Filters indices by evaluating a formula.
 * @param f The formula to evaluate for each index.
 * @param receptacle A `_SimpleList` to store the indices for which the formula
 * is true.
 * @param isVert If true, iterates through site indices. If false, iterates
 * through species indices.
 * @param limit The upper bound of the indices to check.
 * @param scope The variable scope to use for evaluating the formula.
 *
 * This function iterates from 0 to `limit - 1`. In each iteration, it sets a
 * context variable
 * (`siteIndex` or `speciesIndex`) to the current index value, evaluates the
 * formula `f`, and if the result is non-zero, adds the index to the
 * `receptacle` list.
 */
void _DataSet::MatchIndices(_Formula &f, _SimpleList &receptacle, bool isVert,
                            long limit, _String const *scope) const {
  _String varName = isVert ? "siteIndex" : "speciesIndex";
  varName = AppendContainerName(varName, scope);
  _Variable *v = CheckReceptacle(&varName, kEmptyString, false);

  for (long i = 0L; i < limit; i++) {
    v->SetValue(new _Constant((hyFloat)i), false, i == 0, nullptr);
    HBLObjectRef res = f.Compute();
    if (res && !CheckEqual(res->Value(), 0.0)) {
      receptacle << i;
    }
  }
  v->SetValue(new _Constant(0.0), false, false, nullptr);
}

//_________________________________________________________

/**
 * @brief Checks if a list of datasets are compatible for a merge operation and
 * creates a merged translation table.
 * @param ref A `_SimpleList` containing the indices of the datasets to check.
 * @param concatOrCombine A flag indicating the type of merge: 1 for
 * concatenation (joining columns), 0 for combination (joining rows).
 * @return A new `_TranslationTable` object if all datasets are compatible,
 * otherwise `nullptr`.
 *
 * This function verifies two conditions:
 * 1. The translation tables of all datasets can be merged into a single,
 * consistent table.
 * 2. The dimensions of the datasets are compatible for the specified operation
 * (e.g., for concatenation, they must have the same number of species). If
 * compatible, it returns a new, merged translation table. Otherwise, it reports
 * an error and returns `nullptr`.
 */
_TranslationTable *_DataSet::CheckCompatibility(_SimpleList const &ref,
                                                char concatOrCombine) {
  _DataSet *currentSet = (_DataSet *)dataSetList(ref.Element(0));

  _TranslationTable *theEnd = new _TranslationTable(*(currentSet->theTT));

  long refNo =
      concatOrCombine ? currentSet->NoOfSpecies() : currentSet->NoOfColumns();
  char emptyStringChar = theEnd->GetSkipChar();

  for (unsigned long k = 1; k < ref.lLength; k++) {
    currentSet = (_DataSet *)dataSetList(ref.Element(k));

    _TranslationTable *tryMe = theEnd->MergeTables(currentSet->theTT);

    if (tryMe) {
      if (emptyStringChar) {
        DeleteObject(theEnd);
        theEnd = tryMe;
        continue;
      } else {
        if ((concatOrCombine && ((long)currentSet->NoOfSpecies() == refNo)) ||
            (!concatOrCombine && ((long)currentSet->NoOfColumns() == refNo))) {
          DeleteObject(theEnd);
          theEnd = tryMe;
          continue;
        }
      }
    }
    _String warningMessage("The data set ");
    warningMessage =
        warningMessage &
        ((_String *)dataSetNamesList(ref.Element(k)))->Enquote() &
        _String(" was found incompatible with one of the following data sets ");
    for (unsigned long i = 0; i < k; i++) {
      if (k) {
        warningMessage = warningMessage & ", ";
      }
      warningMessage = warningMessage &
                       ((_String *)dataSetNamesList(ref.Element(k)))->Enquote();
    }
    HandleApplicationError(warningMessage);
    DeleteObject(tryMe);
    DeleteObject(theEnd);
    return nullptr;
  }

  return theEnd;
}

//_________________________________________________________

/**
 * @brief Calculates frequencies of character patterns (e.g., codons, amino
 * acids) from the dataset.
 *
 * @param unit The number of characters that form a single unit (e.g., 3 for a
 * codon).
 * @param atom The number of characters that form an atomic pattern within the
 * unit (e.g., 1 for a nucleotide, 3 for a codon).
 * @param posSpec If true, frequencies are position-specific within the unit. If
 * false, frequencies are pooled across all positions.
 * @param hSegmentation A `_SimpleList` of sequence indices (rows) to include in
 * the calculation.
 * @param vSegmentation A `_SimpleList` of site indices (columns) to include in
 * the calculation.
 * @param countGaps If true, ambiguity codes that resolve to gaps are included
 * in the frequency counts.
 * @return A `_Matrix` where rows correspond to character patterns and columns
 * correspond to positions (if `posSpec` is true).
 *
 * This is a powerful and complex function for extracting detailed frequency
 * information from the dataset. It can be used to calculate things like codon
 * usage, amino acid frequencies, or nucleotide frequencies, optionally
 * constrained to specific subsets of sequences and sites.
 */
_Matrix *_DataSet::HarvestFrequencies(unsigned char unit, unsigned char atom,
                                      bool posSpec, _SimpleList &hSegmentation,
                                      _SimpleList &vSegmentation,
                                      bool countGaps) const {

  if (hSegmentation.empty() ||
      vSegmentation.countitems() < unit) { // revert to default (all data)
    if (hSegmentation.empty()) {
      hSegmentation.Populate(NoOfSpecies(), 0, 1);
    }
    if (vSegmentation.countitems() < unit) {
      vSegmentation.Clear();
      vSegmentation.Populate(GetNoTypes(), 0, 1);
    }
  }

  if (atom == 0 ||
      unit % atom > 0) { // 20120814 SLKP: changed this behavior to throw errors
    HandleApplicationError(
        _String("Atom must be non-zero and divide unit, had ") &
        _String((long)unit) & "/" & _String((long)atom));
    return new _Matrix(1, 1);
  }

  _Matrix *out = new _Matrix(ComputePower(theTT->LengthOfAlphabet(), atom),
                             posSpec ? unit / atom : 1, false, true);

  long positions = unit / atom, static_store[HYPHY_SITE_DEFAULT_BUFFER_SIZE];

  _String unit_for_counting((unsigned long)atom);

  for (unsigned long site_pattern = 0UL;
       site_pattern + unit <= vSegmentation.lLength;
       site_pattern += unit) { // loop over the set of segments

    for (unsigned long primary_site = site_pattern;
         primary_site < site_pattern + unit; primary_site += atom) {

      long index_in_pattern = (primary_site - site_pattern) / atom;

      for (unsigned long sequence_index = 0;
           sequence_index < hSegmentation.lLength; sequence_index++) {
        // loop down each column

        unsigned long mapped_sequence_index =
            hSegmentation.list_data[sequence_index];
        // build atomic probabilities

        for (unsigned long m = 0UL; m < atom; m++) {
          unit_for_counting.set_char(
              m, (*this)(vSegmentation.list_data[primary_site + m],
                         mapped_sequence_index, atom));
        }

        long resolution_count = theTT->MultiTokenResolutions(
            unit_for_counting, static_store, countGaps);

        if (resolution_count > 0L) {

          hyFloat normalized = 1. / resolution_count;

          for (long resolution_index = 0UL; resolution_index < resolution_count;
               resolution_index++) {
            out->theData[posSpec ? static_store[resolution_index] * positions +
                                       index_in_pattern
                                 : static_store[resolution_index]] +=
                normalized;
          }
        }
      }
    }
  }

  // scale the matrix now

  unsigned long row_count = out->GetHDim(), column_count = out->GetVDim();

  for (unsigned long column = 0UL; column < column_count;
       column++) { // normalize each _column_ to sum to 1.
    hyFloat sum = 0.0;

    for (unsigned long row = 0UL; row < row_count; row++) {
      sum += out->theData[row * column_count + column];
    }

    for (unsigned long row = 0UL; row < row_count; row++) {
      out->theData[row * column_count + column] /= sum;
    }
  }

  return out;
}

//_______________________________________________________________________

/**
 * @brief Parses a complex partition string to select a subset of sites or
 * sequences.
 *
 * @param input2 The partition string. This can be a simple list of indices
 * (e.g., "1,3,5-10"), a formula, a regular expression, or a callback function.
 * @param target A `_SimpleList` that will be populated with the indices of the
 * selected sites or sequences.
 * @param isVertical If true, the partition applies to sites (columns). If
 * false, it applies to sequences (rows).
 * @param unit_length The size of the unit for certain partitioning modes (e.g.,
 * 3 for codons).
 * @param additionalFilter An optional pre-filter to apply before this
 * partition.
 * @param otherDimension An optional filter for the other dimension (e.g., if
 * partitioning sites, this could filter sequences).
 * @param scope The variable scope for evaluating formulas.
 *
 * This is a highly flexible and powerful function for selecting data. It can
 * handle various complex selection criteria, making it a central part of the
 * data analysis capabilities.
 */
void _DataSet::ProcessPartition(_String const &input2, _SimpleList &target,
                                bool isVertical, int unit_length,
                                _SimpleList const *additionalFilter,
                                _SimpleList const *otherDimension,
                                _String const *scope) const {
  // TODO SLKP : 20170928 this needs serious cleanup and testing

  if (input2.empty()) {
    return;
  }
  // decide if the input is an enumeration or a formula
  long totalLength;

  if (additionalFilter) {
    totalLength = additionalFilter->countitems();
  } else {
    totalLength = isVertical ? theMap.countitems() : noOfSpecies;
  }

  _String input(input2);

  if (!input.IsALiteralArgument(true)) { // not a literal argument

    _Formula fmla, lhs;
    _FormulaParsingContext fpc;
    fpc.setScope(scope);

    long outcome = Parse(&fmla, input, fpc, &lhs);

    if (outcome != HY_FORMULA_EXPRESSION) {
      HandleApplicationError(input.Enquote() &
                             _String(" is an invalid partition specification"));
      return;
    }
    HBLObjectRef fV = fmla.Compute();
    if (fV && fV->ObjectClass() == STRING) {
      ProcessPartition(((_FString *)fV)->get_str().Enquote(), target,
                       isVertical, unit_length, additionalFilter, nullptr,
                       scope);
    } else {
      _DataSet::MatchIndices(fmla, target, isVertical, totalLength, scope);
    }
  } else { // an explicit enumeration or a regular expression

    // check to see if argument is a callback

    bool is_regexp = input(0) == '/' && input(-1) == '/';
    long is_hbl_function = -1L;

    if (!is_regexp) {
      is_hbl_function = hyphy_global_objects::FindBFFunctionName(input);
      if (is_hbl_function >= 0) {
        if (GetBFFunctionArgumentCount(is_hbl_function) != 2) {
          HandleApplicationError(
              input.Enquote() &
              _String(
                  " is not a valid callback function: must have two arguments "
                  "(name, sequence for sites; string, frequencies for sites)"));
          return;
        }
      }
    }

    if (is_regexp || is_hbl_function >= 0L) {
      // a regular expression or a callback
      std::regex *regex = nullptr;
      _Formula filter_formula;

      if (is_regexp) {
        input.Trim(1, input.length() - 2);
        int errCode = 0;
        std::regex *regex_ptr = _String::PrepRegExp(input, errCode, true);
        if (errCode) {
          HandleApplicationError(_String::GetRegExpError(errCode));
          return;
        }
        regex = regex_ptr;
      }
      // now set do the matching
      // using only the sites that are specced in the additionalFilter

      if (!isVertical) { // partitioning sequences
        _FString *string_object = nullptr, *string_name = nullptr;
        if (!is_regexp) {
          filter_formula.GetList() < new _Operation() < new _Operation() <
              new _Operation(kEmptyString, -is_hbl_function - 1L);

          string_object = new _FString;
          string_name = new _FString;
        }

        const long loop_limit =
            additionalFilter ? additionalFilter->lLength : totalLength;

        for (long specCount = 0L; specCount < loop_limit; specCount++) {
          _String pattern((unsigned long)theMap.countitems());
          long seqPos = additionalFilter ? additionalFilter->Element(specCount)
                                         : specCount;

          if (otherDimension) {
            for (unsigned long seqSlider = 0L;
                 seqSlider < otherDimension->lLength; seqSlider++) {
              pattern.set_char(seqSlider,
                               GetSite(otherDimension->Element(seqSlider))
                                   ->get_char(seqPos));
            }
          } else {
            for (unsigned long seqSlider = 0L; seqSlider < theMap.lLength;
                 seqSlider++) {
              pattern.set_char(seqSlider, GetSite(seqSlider)->get_char(seqPos));
            }
          }

          if (is_regexp) {
            if (pattern.RegExpMatch(regex, 0UL).countitems()) {
              target << specCount;
            }
          } else {
            string_object->SetStringContent(new _StringBuffer(pattern));
            string_name->SetStringContent(
                new _StringBuffer(*GetSequenceName(seqPos)));
            filter_formula.GetIthTerm(1)->SetNumber(string_object);
            filter_formula.GetIthTerm(0)->SetNumber(string_name);
            if (!CheckEqual(0., filter_formula.Compute()->Value())) {
              target << specCount;
            }
          }
        }

        if (!is_regexp) {
          filter_formula.GetIthTerm(0)->SetNumber(nullptr);
          filter_formula.GetIthTerm(1)->SetNumber(nullptr);
          DeleteObject(string_object);
          DeleteObject(string_name);
        }
      } else {

        auto map_site = [](const _Site *site, _String &buffer,
                           _SimpleList const *mapper) -> void {
          mapper->Each([&](long value, unsigned long index) -> void {
            buffer.set_char(index, site->char_at(value));
          });
        };

        bool *eligibleMarks = nullptr;

        if (is_regexp) {
          eligibleMarks = new bool[lLength];
          if (additionalFilter) {
            InitializeArray(eligibleMarks, lLength, false);
            for (unsigned long siteIndex = 0;
                 siteIndex < additionalFilter->lLength; siteIndex++) {
              eligibleMarks
                  [theMap.list_data[additionalFilter->list_data[siteIndex]]] =
                      true;
            }
          } else {
            InitializeArray(eligibleMarks, lLength, true);
          }

          _String *tempString = nullptr;
          _SimpleList matches;
          if (otherDimension) {
            tempString = new _String(otherDimension->countitems());
          }

          for (unsigned long siteCounter = 0; siteCounter < lLength;
               siteCounter++)
            if (eligibleMarks[siteCounter]) {
              matches.Clear();
              if (otherDimension) {
                map_site(GetSite(siteCounter), *tempString, otherDimension);
                matches = tempString->RegExpMatch(regex, 0UL);
              } else {
                matches = GetSite(siteCounter)->RegExpMatch(regex, 0UL);
              }
              if (matches.empty()) {
                eligibleMarks[siteCounter] = false;
              }
            }

          DeleteObject(tempString);
          if (additionalFilter) {
            for (unsigned long afi = 0; afi < additionalFilter->lLength;
                 afi++) {
              if (eligibleMarks
                      [theMap.list_data[additionalFilter->list_data[afi]]]) {
                target << afi;
              }
            }
          } else {
            theMap.Each([&target, eligibleMarks](long site_pattern,
                                                 unsigned long index) -> void {
              if (eligibleMarks[site_pattern]) {
                target << index;
              }
            });
          }
        } else {

          long freq_dimension = ComputePower(GetCharDimension(), unit_length);

          if (freq_dimension > 0xffff) {
            HandleApplicationError("The dimension of the character space is "
                                   "too high for callback filtering");
            return;
          }
          eligibleMarks = new bool[theMap.lLength];

          if (additionalFilter) {
            InitializeArray(eligibleMarks, theMap.lLength, false);
            for (unsigned long siteIndex = 0;
                 siteIndex < additionalFilter->lLength; siteIndex++) {
              eligibleMarks[additionalFilter->list_data[siteIndex]] = true;
            }
          } else {
            InitializeArray(eligibleMarks, theMap.lLength, true);
          }

          filter_formula.GetList() < new _Operation() < new _Operation() <
              new _Operation(kEmptyString, -is_hbl_function - 1L);

          _Matrix *strings = nullptr, *frequencies = nullptr;

          _List string_list, string_storage;

          if (otherDimension) {
            for (int i = 0; i < unit_length; i++) {
              string_storage < new _String(otherDimension->countitems());
            }
          }

          _SimpleList sites(unit_length, 0, 0), sequences;

          if (otherDimension) {
            sequences = *otherDimension;
          } else {
            sequences.Populate((unsigned long)NoOfSpecies(), 0, 1);
          }

          for (long siteCounter = 0L;
               siteCounter + unit_length <= (long)theMap.lLength;
               siteCounter += unit_length) {
            long unit_space = 0L;
            string_list.Clear();
            for (; unit_space < unit_length; unit_space++) {
              if (eligibleMarks[siteCounter + unit_space]) {
                sites[unit_space] = siteCounter + unit_space;
                if (otherDimension) {
                  map_site((_Site *)GetSite(siteCounter + unit_space),
                           *(_String *)string_storage.GetItem(unit_space),
                           otherDimension);
                } else {
                  string_list << GetSite(siteCounter + unit_space);
                }
              } else {
                break;
              }
            }

            if (unit_space == unit_length) {
              // eligibleMarks[siteCounter + unit_space]
              BatchDelete(strings, frequencies);
              if (otherDimension) {
                strings = new _Matrix(string_storage);
              } else {
                strings = new _Matrix(string_list);
              }

              filter_formula.GetIthTerm(0)->SetNumber(strings);

              //(unsigned char unit, unsigned char atom, bool posSpec,
              //_SimpleList& hSegmentation, _SimpleList& vSegmentation, bool
              // countGaps)
              frequencies = HarvestFrequencies((unsigned char)unit_space,
                                               (unsigned char)unit_space, false,
                                               sequences, sites, false);
              filter_formula.GetIthTerm(1)->SetNumber(frequencies);

              if (!CheckEqual(0., filter_formula.Compute()->Value())) {
                continue;
              }
            }

            for (unit_space = 0L; unit_space < unit_length; unit_space++) {
              eligibleMarks[siteCounter + unit_space] = false;
            }
          }

          theMap.Each(
              [&target, eligibleMarks](long, unsigned long index) -> void {
                if (eligibleMarks[index]) {
                  target << index;
                }
              });
          // strings && frequencies will be cleaned up by the destructor of
          // filter_formula
        }

        delete[] eligibleMarks;
      }
      if (regex) {
        _String::FlushRegExp(regex);
      }

    } else {
      input = input.KillSpaces();
      // now process the string
      long count = 0L, anchor;

      _SimpleList numbers, links;

      numbers.RequestSpace(1024);
      links.RequestSpace(1024);

      // first check if it is has a comb filter

      if (input(0) == '<' && input(-1) == '>') {
        for (count = 1; count < (long)input.length() - 1; count++) {
          if (input.char_at(count) != '0') {
            numbers << count - 1;
          }
        }
        if (numbers.countitems()) {
          long k = input.length() - 2; // step size
          anchor = 0;
          if (totalLength == -1) {
            totalLength = theMap.lLength;
          }
          while (anchor < totalLength - k) {
            for (count = 0; count < (long)numbers.lLength; count++) {
              target << anchor + numbers.list_data[count];
            }
            anchor += k;
          }
          if ((k = totalLength - 1 - anchor)) {
            for (count = 0; count < (long)numbers.lLength; count++) {
              if (numbers.list_data[count] > k) {
                break;
              }
              target << anchor + numbers.list_data[count];
            }
          }
          return;
        }
      }

      while (count < (long)input.length()) {
        anchor = count;

        for (; count < (long)input.length() && isdigit(input.char_at(count));
             count++)
          ;

        long aNumber = (input.Cut(anchor, count - 1)).to_long();

        if (aNumber < 0) {
          _String warnMsg(
              "A negative number was found in partition specification: ");
          HandleAlignmentValidationError(warnMsg & input.Cut(0, anchor - 1) &
                                         '?' & input.Cut(anchor, -1));
          target.Clear();
          return;
        }
        numbers << aNumber;

        char current_char = input.char_at(count);

        if (current_char == '<' || current_char == '>') {
          HandleAlignmentValidationError(
              _String("A comb partition cannot be combined with other types. "
                      "The entire partition is reset to first..last") &
              input.Cut(0, anchor - 1) & '?' & input.Cut(anchor, -1));
          target.Clear();
          return;
        }

        if (current_char == '&') {
          links << numbers.lLength;
        }

        // TODO SLKP 20171001 this needs to be checked for correctness
        if (current_char == ',' ||
            count == (long)input.length()) { // wrap it up dude
          if (numbers.countitems() == 1) {
            target << numbers(0);
          } else {
            if (links.empty()) {
              if (numbers[0] > numbers[1]) { // backward order
                for (long k = numbers[0]; k >= numbers[1]; k--) {
                  target << k;
                }
              } else {
                for (long k = numbers[0]; k <= numbers[1]; k++) {
                  target << k;
                }
              }
            } else {
              // linked locations
              if (links.countitems() != (numbers.countitems() - 2) / 2) {
                HandleAlignmentValidationError(
                    "A part of the partition specification has not "
                    "been understood and is being skipped.");
                target.Clear();
                return;
              } else {
                _SimpleList signs;
                signs << (numbers(0) < numbers(1) ? 1 : -1);
                for (unsigned long k = 0; k < links.lLength; k += 2) {
                  signs << (numbers(links(k)) < numbers(links(k + 1)) ? 1 : -1);
                }

                for (long k = numbers(0), l = 0;
                     signs(0) * k <= signs(0) * numbers(1);
                     k += signs(0), l++) {
                  target << numbers(0) + l * signs(0);
                  for (unsigned long m = 0; m < links.lLength; m++) {
                    target << numbers(links(m)) + l * signs(m + 1);
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

//_________________________________________________________

/**
 * @brief Concatenates multiple datasets into a single new dataset.
 * @param ref A `_SimpleList` of indices of the datasets to concatenate.
 * @return A new `_DataSet` object containing the concatenated data.
 *
 * This function joins multiple datasets side-by-side (column-wise).
 * If the datasets have a different number of sequences (rows), the resulting
 * dataset will have the same number of sequences as the largest dataset, and
 * the smaller datasets will be padded with the skip character. The translation
 * tables of all datasets are merged.
 *
 * @note The caller is responsible for managing the memory of the returned
 * `_DataSet` object.
 */
_DataSet *_DataSet::Concatenate(_SimpleList const &ref)

// concatenates (adds columns together) several datasets
// in case the number of species in the datasets are different the deficiencies
// will be padded by omission symbols in case translation tables are different,
// they will be merged, provided it can be done, otherwise the incompatible
// datasets will be ignored during this operation.

{
  _TranslationTable *jointTable = CheckCompatibility(ref, 1);
  if (!jointTable) {
    return new _DataSet;
  }

  _DataSet *bigDataSet = new _DataSet;

  bigDataSet->theTT = jointTable;

  // pass one - determine the max max number of species present and what dataset
  // are they coming from

  long maxSpecies = 0, indexOfLargestDataset = 0, siteIndex;

  _DataSet *currentSet;

  char emptyStringSlot = jointTable->GetSkipChar();

  for (unsigned long i = 0; i < ref.lLength; i++) {
    currentSet = (_DataSet *)dataSetList(ref(i));

    long specCount = currentSet->NoOfSpecies(),
         siteCount = currentSet->NoOfColumns();

    if (specCount > maxSpecies) {
      maxSpecies = specCount;
      indexOfLargestDataset = i;
    }
    for (long j = 0; j < siteCount; j++) {
      bigDataSet->AddSite((*currentSet)(j, 0, 1));
    }
  }

  for (long k = 1; k < maxSpecies; k++) {
    siteIndex = 0;
    for (unsigned long i = 0; i < ref.lLength; i++) {
      currentSet = (_DataSet *)dataSetList(ref.list_data[i]);

      long cns = currentSet->NoOfSpecies(), cnc = currentSet->NoOfColumns();

      if (cns <= k)
        for (long j = 0; j < cnc; j++, siteIndex++) {
          bigDataSet->Write2Site(siteIndex, emptyStringSlot);
        }
      else
        for (long j = 0; j < cnc; j++, siteIndex++) {
          bigDataSet->Write2Site(siteIndex, (*currentSet)(j, k, 1));
        }
    }
  }

  currentSet = (_DataSet *)dataSetList(ref(indexOfLargestDataset));
  for (long i = 0L; i < maxSpecies; i++) {
    bigDataSet->AddName(*currentSet->GetSequenceName(i));
  }

  bigDataSet->Finalize();
  bigDataSet->SetNoSpecies(maxSpecies);
  return bigDataSet;
}

//_________________________________________________________

/**
 * @brief Combines multiple datasets into a single new dataset.
 * @param ref A `_SimpleList` of indices of the datasets to combine.
 * @return A new `_DataSet` object containing the combined data.
 *
 * This function joins multiple datasets top-to-bottom (row-wise).
 * If the datasets have a different number of sites (columns), the resulting
 * dataset will have the same number of sites as the widest dataset, and the
 * narrower datasets will be padded with the skip character. The translation
 * tables of all datasets are merged.
 *
 * @note The caller is responsible for managing the memory of the returned
 * `_DataSet` object.
 */
_DataSet *_DataSet::Combine(_SimpleList const &ref) {

  // combines (adds rows together) several datasets
  // in case the number of species in the datasets are different the
  // deficiencies will be padded by omission symbols in case translation tables
  // are different, they will be merged, provided it can be done, otherwise the
  // incompatible datasets will be ignored during this operation.

  _TranslationTable *joint_table = CheckCompatibility(ref, 0);

  if (!joint_table) {
    return new _DataSet;
  }

  _DataSet *combined_data = new _DataSet;
  combined_data->theTT = joint_table;

  // pass one - determine the max max number of sites present and what dataset
  // are they coming from

  unsigned long max_sites = 0UL, total_species_count = 0UL;

  char emptyStringSlot = joint_table->GetSkipChar();

  for (unsigned long set_index = 0UL; set_index < ref.lLength; set_index++) {
    _DataSet const *current_data_set =
        (_DataSet const *)dataSetList(ref.Element(set_index));
    StoreIfGreater(max_sites, current_data_set->NoOfColumns());
    total_species_count += current_data_set->NoOfSpecies();
  }

  for (unsigned long set_index = 0UL; set_index < ref.lLength; set_index++) {
    _DataSet const *current_data_set =
        (_DataSet const *)dataSetList(ref.Element(set_index));
    unsigned long sites_in_this_set = current_data_set->NoOfColumns(),
                  sequences_in_this_set = current_data_set->NoOfSpecies();

    for (unsigned long seq_index = 0UL; seq_index < sequences_in_this_set;
         seq_index++) {
      combined_data->AddName(*current_data_set->GetSequenceName(seq_index));
      if (seq_index == 0UL && set_index == 0UL) {
        /* use AddSite write out the first sequence */
        unsigned long site_index = 0UL;
        for (site_index = 0UL; site_index < sites_in_this_set; site_index++) {
          combined_data->AddSite((*current_data_set)(site_index, 0UL, 1));
        }
        for (; site_index < max_sites; site_index++) {
          combined_data->AddSite(emptyStringSlot);
        }
      } else {
        /* use Write2Site to create subsequence sequences */
        unsigned long site_index = 0UL;
        for (site_index = 0UL; site_index < sites_in_this_set; site_index++) {
          combined_data->Write2Site(
              site_index, (*current_data_set)(site_index, seq_index, 1));
        }
        for (; site_index < max_sites; site_index++) {
          combined_data->Write2Site(site_index, emptyStringSlot);
        }
      }
    }
  }

  combined_data->Finalize();
  combined_data->SetNoSpecies(total_species_count);
  return combined_data;
}

//_______________________________________________________________________

_String *_DataSet::GetSequenceCharacters(long seqID) const {

  unsigned long upTo = NoOfColumns();
  _StringBuffer *aSequence = new _StringBuffer(upTo);

  if (seqID >= 0L && seqID < (long)noOfSpecies) {
    for (unsigned long k2 = 0UL; k2 < upTo; k2++) {
      (*aSequence) << GetSite(k2)->char_at(seqID);
    }
  }
  aSequence->TrimSpace();
  return aSequence;
}

//_________________________________________________________

bool StoreADataSet(_DataSet *ds, _String *setName) {
  if (!setName->IsValidIdentifier(fIDAllowCompound | fIDAllowFirstNumeric)) {
    HandleApplicationError(
        setName->Enquote() &
        " is not a valid identifier while constructing a DataSet");
    return false;
  }

  long type = HY_BL_DATASET, index;
  _DataSet *existing_ds = (_DataSet *)_HYRetrieveBLObjectByNameMutable(
      *setName, type, &index, false, false);

  if (!existing_ds) {
    dataSetNamesList << setName;
    dataSetList < ds;
  } else {

    bool isDifferent =
        existing_ds->NoOfSpecies() != ds->NoOfSpecies() ||
        existing_ds->NoOfColumns() != ds->NoOfColumns() ||
        existing_ds->NoOfUniqueColumns() != ds->NoOfUniqueColumns() ||
        (existing_ds->GetTT() != ds->GetTT() &&
         !(*existing_ds->GetTT() == *ds->GetTT()));

    for (AVLListXLIteratorKeyValue filter_key_value :
         ObjectIndexer(HY_BL_DATASET_FILTER)) {
      _DataSetFilter *filter = (_DataSetFilter *)filter_key_value.get_object();
      if (filter->GetData() == existing_ds) {
        if (isDifferent) {
          ReportWarning(
              _String("Overwriting dataset '") & *setName &
              "' caused DataSetFilter " &
              GetFilterName(filter_key_value.get_index())->Enquote('\'') &
              " to be deleted");
          DeleteDataFilter(filter_key_value.get_index());
        } else {
          filter->SetData(ds);
        }
      }
    }

    dataSetList.Replace(index, ds, false);
  }

  CheckReceptacleAndStore(*setName & ".mapping", kEmptyString, false,
                          new _MathObject, false);

  if (hy_env::EnvVariableTrue(hy_env::normalize_sequence_names)) {
    _List _id_mapping;
    _AVLListXL id_mapping(&_id_mapping);
    bool did_something = false;

    for (unsigned long i = 0UL; i < ds->NoOfSpecies(); i++) {
      _String *old_name = new _String(*ds->GetSequenceName(i));
      if (!old_name->IsValidIdentifier(fIDAllowFirstNumeric)) {
        *ds->GetSequenceName(i) =
            ds->GetSequenceName(i)->ConvertToAnIdent(fIDAllowFirstNumeric);
        did_something = true;
      }
      if (id_mapping.Find(ds->GetSequenceName(i)) >= 0) {
        _String new_name(*ds->GetSequenceName(i));
        long suffix = 1L;
        do {
          new_name = *ds->GetSequenceName(i) & "_" & suffix++;
        } while (id_mapping.Find(&new_name) >= 0);
        *ds->GetSequenceName(i) = new_name;
        did_something = true;
      }

      ds->GetSequenceName(i)->AddAReference();
      id_mapping.Insert(ds->GetSequenceName(i), (long)old_name, false, false);
    }

    if (did_something) {
      _AssociativeList *mapping = new _AssociativeList();

      for (AVLListXLIteratorKeyValue filter_key_value :
           AVLListXLIterator(&id_mapping)) {
        // printf ("%d => %s\n", filter_key_value.get_index(),
        // ((_String*)filter_key_value.get_object())->get_str());
        mapping->MStore(
            *(_String *)id_mapping.Retrieve(filter_key_value.get_index()),
            *(_String *)filter_key_value.get_object());
      }

      CheckReceptacleAndStore(*setName & ".mapping", kEmptyString, false,
                              mapping, false);
    }
  }

  CheckReceptacleAndStore(*setName & ".species", kEmptyString, false,
                          new _Constant(ds->NoOfSpecies()), false);
  CheckReceptacleAndStore(*setName & ".sites", kEmptyString, false,
                          new _Constant(ds->NoOfColumns()), false);
  CheckReceptacleAndStore(*setName & ".unique_sites", kEmptyString, false,
                          new _Constant(ds->NoOfUniqueColumns()), false);

  return true;
}

//_________________________________________________________
//_________________________________________________________
// reading the data set file in here

//_________________________________________________________
void checkTTStatus(FileState *fs) { // check whether the translation table needs
                                    // to be refreshed}
  if (fs->translationTable == &hy_default_translation_table) {
    fs->translationTable =
        (_TranslationTable *)hy_default_translation_table.makeDynamic();
  }
}
//_________________________________________________________
void processCommand(_String *s, FileState *fs) {

  static const _List _CommandList(new _String("BASESET"), new _String("FORMAT"),
                                  new _String("RAWLINE"), new _String("REPEAT"),
                                  new _String("TOKEN"));

  static const _String kBase20("BASE20"), kPHYLIPi("PHYLIPI"),
      kPHYLIPs("PHYLIPS"), kRAW("RAW");

  long f = -1, command_index;

  for (command_index = 0L; command_index < (long)_CommandList.countitems();
       ++command_index) {
    f = s->Find(*(_String *)_CommandList.GetItem(command_index));
    if (f != kNotFound) {
      break;
    }
  }

  try {
    if (f == -1) { // unrecognized command
      return;
    } else {
      // trim the string
      // s->Trim (f+((_String*)CommandList(i))->Length(),-1);

      f = s->Find(
          ":",
          f + 1L + ((_String *)_CommandList.GetItem(command_index))->length());

      if (f == kNotFound) { // poorly formed command
        throw(s->Enquote('[', ']') & " was not of the form COMMAND : DATA");
      }

      if (command_index >= 1 && command_index <= 3) {
        long start = f + 1;
        long end = s->ExtractEnclosedExpression(start, '"', '"', 0);

        if (end == kNotFound || end - start <= 2L) {
          throw(s->Enquote('[', ']') &
                " was not of the form COMMAND : \"DATA\" (missing quotes)");
        }
        s->Trim(start + 1L, end - 1L);
      } else {
        s->Trim(f + 1L, kStringEnd);
      }

      // 's' should now contain only the payload of the command

      switch (command_index) {

      case 4: { // new token
        checkTTStatus(fs);
        // attempt to extract a token. Looking for (e.g):   "c" = "AC"

        _SimpleList matches = s->RegExpMatch(
            "\\\"([a-z,A-Z])\\\"\\ *=\\ *\\\"([a-z,A-Z]+)\\\"", false, true);
        if (matches.countitems() == 6) {
          fs->translationTable->AddTokenCode(
              (char)matches.get(2), s->Cut(matches.get(3), matches.get(4)));
        } else {
          throw(s->Enquote('[', ']') &
                " was not of the form \"token\"=\"translation\"");
        }
      } break;

      case 0: { // new code set, e.g  "ACGU"
        checkTTStatus(fs);
        // erase previous char definitions
        fs->translationTable->translationsAdded.Clear();
        fs->translationTable->tokensAdded = "";
        if (*s != kBase20) {
          long start = 0;
          long end = s->ExtractEnclosedExpression(start, '"', '"', 0);
          if (end == kNotFound || end - start <= 3L) {
            throw(s->Enquote('[', ']') &
                  " was not of the form \"at least two letters\"");
          }
          // TODO : there is no check that baseset is actually valid (e.g. no
          // duplicate characters etc)
          fs->translationTable->AddBaseSet(s->Cut(start + 1, end - 1));
        } else {
          fs->translationTable->AddBaseSet(kEmptyString);
          fs->translationTable->baseLength = 20;
        }
      } break;

      case 1:                 // FORMAT
        if (*s == kPHYLIPi) { // PHYLIP Interleaved
          fs->fileType = 1;
          fs->interleaved = TRUE;
        } else if (*s == kPHYLIPs) { // PHYLIP sequential
          fs->fileType = 1;
          fs->interleaved = FALSE;
        }
        if (*s == kRAW) { // RAW Sequential Data (as in NEXUS)
          fs->fileType = 2;
          fs->interleaved = FALSE;
        }
        fs->autoDetect = false;
        break;

      case 3: // REPEAT CHAR
        fs->repeat = s->get_char(0);
        break;

      case 2: // RAWLINE template e.g 1,-1 skips one word at the beginning and
              // one word at the end
        _List chips(s, ',');
        chips.ForEach([&](BaseRef number, unsigned long) -> void {
          fs->rawLinesFormat << ((_String *)number)->to_long();
        });
        break;
      }
    }
  } catch (const _String &warning) {
    HandleAlignmentValidationError(warning);
  }
}
//_________________________________________________________

void FilterRawString(_String &s, FileState *fs, _DataSet &ds) {
  // TODO: SLKP 20180803 this needs to be tested or deprecated.
  s.CompressSpaces();
  _List words(s.Tokenize(" "));

  long current_start = 0L, current_end = (long)words.countitems() - 1L;

  fs->rawLinesFormat.Each([&](long word, unsigned long) -> void {
    if (word > 0L) {
      current_start += word;
    } else {
      if (word < 0L) {
        current_end += word;
      } else {
        if (current_start < current_end) {
          ds.AddName(*((_String *)words.GetItem(current_start)));
          current_start++;
        }
      }
    }
  });

  if (current_start >= current_end) {
    s = kEmptyString;
  } else {
    s = words.Join(" ", current_start, current_end);
  }
}
//_________________________________________________________________________________________________

//_________________________________________________________________________________________________

void ProcessTree(FileState *fState, hyFile *f, _StringBuffer &CurrentLine) {

  // TODO SLKP 20180921 this does extra work to read in the tree string multiple
  // times; the solution is to have a proper buffer wrapper, and to

  _StringBuffer *tree_string = new _StringBuffer(128L);
  long start_index = 0,
       end_index = CurrentLine.ExtractEnclosedExpression(
           start_index, '(', ')', fExtractRespectQuote | fExtractRespectEscape);

  while (start_index == kNotFound || end_index == kNotFound) {
    _StringBuffer next_line;
    ReadNextLine(f, &next_line, fState, false);
    CurrentLine = CurrentLine & next_line;
    start_index = 0L;
    end_index = CurrentLine.ExtractEnclosedExpression(
        start_index, '(', ')', fExtractRespectQuote | fExtractRespectEscape);
  }

  if (start_index == kNotFound || end_index == kNotFound) {
    HandleAlignmentValidationError(tree_string->Enquote() &
                                   " has mimatched '(' and ')'");
    DeleteObject(tree_string);
  } else {
    *tree_string << CurrentLine.Cut(start_index, end_index);
    tree_string->TrimSpace();
    CurrentLine.Trim(end_index + 1, kStringEnd);
    hy_env::EnvVariableSetNamespace(hy_env::data_file_tree,
                                    new HY_CONSTANT_TRUE, fState->theNamespace,
                                    false);
    hy_env::EnvVariableSetNamespace(hy_env::data_file_tree_string,
                                    new _FString(tree_string), nil, false);
  }
}

//_________________________________________________________

long ProcessLine(_String &s, FileState *fs, _DataSet &ds) {
  long sitesAttached = 0L;

  try {
    s.Each([&](char letter, unsigned long) -> void {
      letter = _uppercase_char(letter);
      if (fs->translationTable->IsCharLegal(letter)) { // go on
        if (fs->curSpecies == 0) {                     // add new column
          ds.AddSite(letter);
          sitesAttached++;
        } else { // append to exisiting column
                 // if (c == fs->skip) continue;
                 //  check to see if this species needs to be padded

          if (letter == fs->repeat) {
            if (fs->curSite + sitesAttached >=
                (long)ds.lLength) { // a dot not matched by a previously read
                                    // character; ignore
              throw sitesAttached;
            }

            letter =
                ((_Site *)(ds._List::operator()(fs->curSite + sitesAttached)))
                    ->get_char(0);
            if (letter == '\0') {
              letter = ((_Site *)(ds._List::operator()(
                            ((_Site *)(ds._List::operator()(fs->curSite +
                                                            sitesAttached)))
                                ->GetRefNo())))
                           ->get_char(0);
            }
          }

          if (fs->curSite + sitesAttached + 1 > fs->totalSitesRead) {
            // pad previous species to full length
            _Site *newS = new _Site(fs->skip);
            newS->AppendNCopies(fs->skip, fs->curSpecies - 1L);
            (*newS) << letter;

            ds.theFrequencies << 1;
            newS->SetRefNo(-1);

            ds < newS;
            fs->totalSitesRead++;
          } else {
            ds.Write2Site(fs->curSite + sitesAttached, letter, fs->skip);
          }

          sitesAttached++;
        }
      }
    });
  } catch (long e) {
    return e;
  }

  // make sure that this species has enough data in it, and if not - pad it with
  // '?'

  if (fs->curSite + sitesAttached < fs->totalSitesRead && fs->interleaved) {
    // pad this species to full length
    for (long j = fs->curSite + sitesAttached; j < fs->totalSitesRead; j++) {
      ds.Write2Site(j, fs->skip);
    }
  }
  if (!fs->curSpecies) {
    fs->totalSitesRead += sitesAttached;
  }
  return sitesAttached;
}

//_________________________________________________________
void PadLine(FileState &fState,
             _DataSet &result) { // make sure that there is enough data in this
                                 // line and if not - "pad" it with '?''s
  if (fState.curSite < fState.totalSitesRead) { // pad line if needed
    // printf ("\nPADLINE %d %d\n", fState.curSpecies,
    // fState.totalSitesRead-fState.curSite);
    for (long j = fState.curSite; j < fState.totalSitesRead; j++) {
      result.Write2Site(j, fState.skip);
    }
  }
}

//_________________________________________________________
void ISelector(FileState &fState, _StringBuffer &CurrentLine,
               _DataSet &result) {
  if (fState.interleaved) { // interleaved file
    if (fState.curSpecies &&
        (!((fState.curSpecies) %
           fState.totalSpeciesExpected))) { // read a chunk of all species
      if (fState.totalSitesRead && !result.InternalStorageMode()) {
        for (long i = fState.curSite; i < fState.totalSitesRead; i++) {
          result.Compact(i);
        }
      }
      fState.curSite = fState.totalSitesRead;
      fState.curSpecies = 0;
      ProcessLine(CurrentLine, &fState, result);
      fState.curSpecies = 1;
      if (!fState.curSite) {
        fState.totalSpeciesRead++;
      }
    } else {
      ProcessLine(CurrentLine, &fState, result);
      if (!fState.curSite) {
        fState.totalSpeciesRead++;
      }
      fState.curSpecies++;
    }
  } else {
    if (fState.curSpecies + 1 < fState.totalSpeciesExpected) {
      fState.curSpecies++;
    }
    if (fState.curSpecies == fState.totalSpeciesRead) {
      PadLine(fState, result);
      fState.curSite = 0;
    }
    if (fState.totalSpeciesRead < fState.totalSpeciesExpected) {
      fState.totalSpeciesRead++;
    }

    fState.curSite += ProcessLine(CurrentLine, &fState, result);
  }
}

//_________________________________________________________
bool SkipLine(_StringBuffer &theLine, FileState *fS) {

  if (theLine.char_at(0) == '/' && theLine.char_at(1) == '/') {
    return true;
  }

  char c = theLine.FirstNonSpace();

  if (c && (!(c == '$' && !fS->acceptingCommands))) {
    return false;
  }

  return true;
}

//_________________________________________________________
void ReadNextLine(hyFile *fp, _StringBuffer *s, FileState *fs, bool,
                  bool upCase) {

  fs->lineBuffer.Reset();
  fs->currentFileLine++;

  char lastc;

  if (fp) {
    lastc = (char)fp->read_char();
  } else {
    lastc = fs->pInSrc < (long)fs->theSource->length()
                ? fs->theSource->char_at(fs->pInSrc++)
                : 0;
  }

  if (fs->fileType != 3) { // not NEXUS - do not skip [..]
    if (fp)
      while (!fp->feof() && lastc != '\r' && lastc != '\n') {
        if (lastc) {
          fs->lineBuffer << lastc;
        }

        lastc = (char)fp->read_char();
      }
    else
      while (lastc && lastc != '\r' && lastc != '\n') {
        fs->lineBuffer << lastc;
        lastc = fs->theSource->char_at(fs->pInSrc++);
      }

  } else {
    if (upCase) {
      lastc = _uppercase_char(lastc);
    }

    while (((fp && !fp->feof()) ||
            (fs->theSource && (fs->pInSrc <= (long)fs->theSource->length()))) &&
           lastc != '\r' && lastc != '\n') {
      if (lastc == '[') {
        if (fs->isSkippingInNEXUS) {
          ReportWarning("Nested comments in NEXUS really shouldn't be used.");
        } else {
          fs->isSkippingInNEXUS = true;
        }
      }
      if (fs->isSkippingInNEXUS) {
        if (lastc == ']') {
          fs->isSkippingInNEXUS = false;
          fs->lineBuffer << ' ';
        }
      } else {
        fs->lineBuffer << lastc;
      }

      if (fp) {
        lastc = (char)fp->read_char();
        if (upCase) {
          lastc = _uppercase_char(lastc);
        }
      } else {
        if (upCase) {
          lastc = _uppercase_char(fs->theSource->char_at(fs->pInSrc++));
        } else {
          lastc = fs->theSource->char_at(fs->pInSrc++);
        }
      }
    }

    if (lastc == '\r' || lastc == '\n') {
      fs->lineBuffer << ' ';
    }
  }

  if ((fp && fp->feof()) ||
      (fs->theSource && fs->pInSrc >= (long)fs->theSource->length())) {
    if (fs->lineBuffer.empty()) {
      s->Clear();
      return;
    }
  }
  // tempBuffer.TrimSpace();
  *s = fs->lineBuffer;

  if (SkipLine(*s, fs)) {
    ReadNextLine(fp, s, fs, false, upCase);
  }

  if (s->nonempty() && s->char_at(s->length() - 1) == '\n') {
    s->Trim(0, (long)s->length() - 2L);
  }
}
//_________________________________________________________
void TrimPhylipLine(_String &CurrentLine, _DataSet &ds) {
  long fNS = CurrentLine.FirstNonSpaceIndex(),
       space2 = CurrentLine.FirstSpaceIndex(fNS + 1);

  if (space2 < 0 && CurrentLine.length() > 10) {
    _String sequence_name(CurrentLine, fNS, CurrentLine.length());
    ds.AddName(sequence_name);
    CurrentLine.Trim(CurrentLine.length(), -1);
  } else {
    // hack for PAML support
    if (space2 > fNS && isspace(CurrentLine.char_at(space2 + 1))) {
      _String sequence_name(CurrentLine, fNS, space2);
      CurrentLine.Trim(space2 + 2, -1); // chop out the name
      ds.AddName(sequence_name);
    } else {
      _String sequence_name(CurrentLine, fNS, fNS + 9);
      CurrentLine.Trim(fNS + 10, -1); // chop out the name
      ds.AddName(sequence_name);
    }
  }
}

//_________________________________________________________
_DataSet *ReadDataSetFile(hyFile *f, char execBF, _String *theS,
                          _String *bfName, _String *namespaceID,
                          _TranslationTable *dT, _ExecutionList *ex) {

  static const _String kNEXUS("#NEXUS"), kDefSeqNamePrefix("Species");

  bool doAlphaConsistencyCheck = true;
  _DataSet *result = new _DataSet;

  try {

    // if (f) flockfile (f);
    if (f)
      f->lock();

    hy_env::EnvVariableSet(hy_env::data_file_tree_string, new _Matrix, false);

    _String savedLine;

    /*_String         CurrentLine = hy_env::data_file_tree_string & "={{}};",
                    savedLine;

    _ExecutionList reset (CurrentLine);
    reset.Execute();*/

#ifdef __HYPHYMPI__
    if (hy_mpi_node_rank == 0L)
#endif
      terminate_execution = false;

    hy_env::EnvVariableSet(hy_env::data_file_tree, new HY_CONSTANT_FALSE,
                           false);

    // initialize the instance of a file state variable
    FileState fState;
    fState.translationTable = dT;
    fState.curSpecies = fState.totalSpeciesRead = fState.totalSitesRead =
        fState.totalSpeciesExpected = fState.totalSitesExpected =
            fState.curSite = fState.currentFileLine = fState.maxStringLength =
                0;
    fState.acceptingCommands = true;
    fState.allSpeciesDefined = false;
    fState.interleaved = false;
    fState.isSkippingInNEXUS = false;
    fState.autoDetect = true;
    fState.fileType = -1;
    fState.baseLength = 4;
    fState.repeat = '.';
    fState.skip = 0;
    fState.theSource = theS;
    fState.pInSrc = 0;
    fState.theNamespace = namespaceID;

    if (!f && !theS) {
      throw _String("ReadDataSetFile received null file AND string references. "
                    "At least one must be specified");
    }
    // done initializing

    if (f) {
      f->rewind();
    }

    _StringBuffer CurrentLine;

    // if (f==NULL) return (_DataSet*)result.makeDynamic();
    //  nothing to do

    // CurrentLine = kEmptyString;

    ReadNextLine(f, &CurrentLine, &fState);
    if (CurrentLine.empty()) {
      throw _String("Empty File Encountered By ReadDataSet.");
    } else {
      if (CurrentLine.BeginsWith(kNEXUS, false)) {
        ReadNexusFile(fState, f, (*result));
        doAlphaConsistencyCheck = false;
      } else {
        long i, j, k, filePosition = -1, saveSpecExpected = 0x7FFFFFFF;
        char c;
        while (CurrentLine.nonempty()) { // stuff to do
          // check if the line has a command in it

          c = CurrentLine.FirstNonSpace();
          while (1) {
            if (fState.acceptingCommands) {
              if (c == '$') { // command line
                processCommand(&CurrentLine, &fState);
                break;
              }
            }

            if (!fState.skip) {
              fState.skip = fState.translationTable->GetSkipChar();
            }
            fState.acceptingCommands = FALSE;

            if (fState.fileType ==
                -1) { // undecided file type - assume it is PHYLIP sequential
              if ((c == '#') || (c == '>')) { // hash-mark format
                fState.fileType = 0;
              } else { // assume this is a sequential PHYLIP file
                fState.fileType = 1;
                fState.interleaved = false;
              }
            }
            // decide what to do next
            // if format is PHYLIP and we do not know the expected dimensions,
            //   we must read those in first
            if (fState.fileType == 1) { // PHYLIP
              if ((filePosition < 0) && (fState.autoDetect)) {
                filePosition = (f ? f->tell()
#ifdef __WINDOZE__
                                        - 1
#endif
                                  : fState.pInSrc);
                savedLine = CurrentLine;
              }

              if (fState.totalSitesExpected == 0 ||
                  fState.totalSpeciesExpected ==
                      0) { // must read dimensions first
                i = CurrentLine.FirstNonSpaceIndex();
                j = CurrentLine.FirstSpaceIndex(i);
                if (j != kNotFound) {
                  k = CurrentLine.FirstNonSpaceIndex(j);
                  if (k != kNotFound) { // could have dimensions
                    saveSpecExpected = fState.totalSpeciesExpected =
                        CurrentLine.Cut(i, j - 1L).to_long();
                    fState.totalSitesExpected =
                        CurrentLine.Cut(k, kStringEnd).to_long();
                  }
                  if (CurrentLine.Find('I', k, kStringDirectionBackward) >=
                      0) { // interleaved
                    fState.interleaved = true;
                  }
                }
              } else
              // now for the data crunching part
              // detect a line, diagnose it and dispatch accordingly
              {
                if (fState.interleaved) {
                  if (fState.totalSpeciesRead < fState.totalSpeciesExpected) {
                    TrimPhylipLine(CurrentLine, (*result));
                  }
                  if (fState.curSite && fState.curSpecies >= saveSpecExpected &&
                      fState.totalSitesRead >= fState.totalSitesExpected) {
                    // reached the end of the data - see maybe there is a tree
                    ReadNextLine(f, &CurrentLine, &fState);
                    if (CurrentLine.nonempty()) {
                      if (CurrentLine.FirstNonSpace() ==
                          '(') { // could be a tree string
                        ProcessTree(&fState, f, CurrentLine);
                      }
                    }
                    break;
                  }

                } else {
                  if (fState.totalSitesRead > fState.totalSitesExpected)
                  // oops - autodetect incorrectly assumed that the file was
                  // sequential
                  {
                    fState.curSpecies = fState.totalSpeciesRead =
                        fState.totalSitesRead = fState.curSite =
                            fState.totalSpeciesExpected =
                                fState.totalSitesExpected =
                                    fState.maxStringLength = 0;
                    fState.allSpeciesDefined = false;
                    fState.interleaved = true;
                    fState.autoDetect = true;

                    if (f) {
                      f->seek(filePosition, SEEK_SET);
                    } else {
                      fState.pInSrc = filePosition;
                    }

                    CurrentLine = savedLine;
                    result->ForEach([](BaseRef site, unsigned long) -> void {
                      ((_Site *)site)->TrimSpace();
                    });

                    result->theNames.Clear();
                    result->theMap.Clear();
                    result->Clear();
                    result->theFrequencies.Clear();
                    continue;
                  }
                  if (fState.totalSpeciesRead == 0) {
                    fState.totalSpeciesExpected = 1;
                    if (!fState.curSite) {
                      TrimPhylipLine(CurrentLine, (*result));
                    }
                  }

                  else if (fState.curSite >= fState.totalSitesExpected) {
                    fState.totalSpeciesExpected++;
                    if (fState.totalSpeciesExpected > saveSpecExpected) {
                      // reached the end of the data - see maybe there is a tree
                      ReadNextLine(f, &CurrentLine, &fState);
                      if (CurrentLine.nonempty()) {
                        if (CurrentLine.FirstNonSpace() ==
                            '(') { // could be a tree string
                          ProcessTree(&fState, f, CurrentLine);
                        }
                      }
                      break;
                    }
                    TrimPhylipLine(CurrentLine, (*result));
                  }
                }

                ISelector(fState, CurrentLine, (*result));
              }
              break;
            }
            // that's all for PHYLIP

            // now handle raw data case
            if (fState.fileType == 2) { // raw data
              FilterRawString(CurrentLine, &fState, (*result));
              if (CurrentLine.nonempty()) {
                break;
              }
              if (ProcessLine(CurrentLine, &fState, (*result))) {
                fState.curSpecies++;
                fState.totalSpeciesRead++;
              }
              break;
            }

            // lastly, handle the auto-detect standard case

            // check to see if the string defines a name
            if (c == '#' || c == '>') {       // a name it is
              if (fState.allSpeciesDefined) { // can't define the species after
                                              // data
                break;
              } else {
                if ((!fState.totalSpeciesRead) &&
                    (fState.totalSpeciesExpected >= 1)) {
                  fState.interleaved = TRUE;
                } else {
                  fState.interleaved = FALSE;
                }
                fState.totalSpeciesExpected++;
                CurrentLine.Trim(CurrentLine.FirstNonSpaceIndex(1), kStringEnd);
                if (CurrentLine.char_at(0) == '#' ||
                    CurrentLine.char_at(0) == '>') {
                  CurrentLine =
                      kDefSeqNamePrefix & _String(fState.totalSpeciesExpected);
                }
                result->AddName(CurrentLine);
              }
              break;
            }
            // check to see if the string defines a tree
            if (c == '(') {
              ProcessTree(&fState, f, CurrentLine);
              ReadNextLine(f, &CurrentLine, &fState);
            }

            // check to see where to stick the incoming line

            if (fState.totalSpeciesExpected == 0) {
              // raw data fed before names defined - skip
              break;
            }
            if (fState.totalSpeciesExpected > 1 &&
                fState.totalSpeciesRead == 0) {
              fState.allSpeciesDefined = TRUE;
            }

            // repeat the structure of PHYLIP reader

            ISelector(fState, CurrentLine, (*result));

            break;
          }

          ReadNextLine(f, &CurrentLine, &fState);
        }
      }
    }

    if (fState.totalSitesRead && fState.interleaved &&
        !result->InternalStorageMode()) {
      for (long i = fState.curSite; i < fState.totalSitesRead; i++) {
        result->Compact(i);
      }
    }

    if ((!fState.interleaved) && (fState.fileType != 2)) {
      PadLine(fState, (*result));
    }

    // make sure interleaved duplications are handled correctly

    result->Finalize();
    result->noOfSpecies = fState.totalSpeciesRead;
    result->theTT = fState.translationTable;

    // check to see if result may be an amino-acid data
    if (doAlphaConsistencyCheck &&
        result->theTT == &hy_default_translation_table) {
      if (result->GetNoTypes() == 0)
      // emptyString data set
      // try binary data
      {
        _TranslationTable *trialTable =
            new _TranslationTable(hy_default_translation_table);
        trialTable->baseLength = 2;
        if (f)
          f->unlock();
        _DataSet *res2 =
            ReadDataSetFile(f, execBF, theS, bfName, namespaceID, trialTable);
        if (res2->GetNoTypes()) {
          DeleteObject(result);
          return res2;
        }
        DeleteObject(res2);
      } else
        // check it out
        if (result->CheckAlphabetConsistency() < 0.5)
        // less than 50% of the data in the alphabet is not in the basic
        // alphabet
        {
          _TranslationTable trialTable(hy_default_translation_table);
          trialTable.baseLength = 20;
          (*result).theTT = &trialTable;
          if ((*result).CheckAlphabetConsistency() < 0.5) {
            CurrentLine = "More than 50% of characters in the data are not in "
                          "the alphabet.";
            (*result).theTT = &hy_default_translation_table;
            ReportWarning(CurrentLine);
          } else {
            (*result).theTT = (_TranslationTable *)trialTable.makeDynamic();
          }
        }
    }
    if (nexusBFBody.nonempty()) {
      if (execBF == 1) {
        lastNexusDataMatrix = result;

        long bfl = GetBFFunctionCount();

        _ExecutionList *nexusBF = ex ? ex : new _ExecutionList;
        if (namespaceID) {
          nexusBF->SetNameSpace(*namespaceID);
        }

        bool do_profile = hy_env::EnvVariableTrue("_PROFILE_NEXUS_LOADS_");

        nexusBF->BuildList(nexusBFBody, nil, false, true);
        //_ExecutionList nexusBF (nexusBFBody,namespaceID);

        if (bfName) {
          nexusBF->sourceFile = *bfName;
        }

        if (do_profile) {
          nexusBF->StartProfile();
        }
        nexusBF->ExecuteAndClean(bfl);
        if (do_profile) {
          CheckReceptacleAndStore("_NEXUS_PROFILE_DATA_", kEmptyString, false,
                                  nexusBF->CollectProfile(), false);
        }

        if (nexusBF != ex) {
          DeleteObject(nexusBF);
        } else {
          ex->ClearExecutionList();
          ex->Clear();
        }
        nexusBFBody = kEmptyString;
      } else if (execBF == 0) {
        nexusBFBody = kEmptyString;
      }
    }
  } catch (const _String &err) {
    DeleteObject(result);
    if (f)
      f->unlock();
    HandleApplicationError(err);
    result = nil;
  }
  if (f)
    f->unlock();
  return result;
}
