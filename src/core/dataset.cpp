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

#include "dataset.h"
#include "translation_table.h"
#include "batchlan.h"
#include "site.h"
#include "global_object_lists.h"

using namespace hyphy_global_objects;


#define DATA_SET_SWITCH_THRESHOLD 100000

extern _TranslationTable defaultTranslationTable;

_DataSet::_DataSet(void) {
  theTT = &defaultTranslationTable;
  streamThrough = nil;
  dsh = nil;
  useHorizontalRep = false;
}

_DataSet::_DataSet(long l)
    : _List((unsigned long)l),
      theFrequencies(
          (unsigned long)l) // with estimated number of sites per file
{
  dsh = nil;
  streamThrough = nil;
  theTT = &defaultTranslationTable;
  useHorizontalRep = false;
}

//_______________________________________________________________________

_DataSet::_DataSet(FILE *f) {
  dsh = nil;
  useHorizontalRep = false;
  theTT = &defaultTranslationTable;
  streamThrough = f;
  theMap << 0; // current sequence
  theMap << 0; // current site
  theMap << 0; // total sites
}

//_______________________________________________________________________

_DataSet::~_DataSet(void) {
  if (theTT != &defaultTranslationTable) {
    DeleteObject(theTT);
  }
}

//_______________________________________________________________________

void _DataSet::Clear(bool) {
  _List::Clear();
  theMap.Clear();
  theFrequencies.Clear();
  theNames.Clear();
  if (theTT != &defaultTranslationTable) {
    DeleteObject(theTT);
    theTT = &defaultTranslationTable;
  }
  noOfSpecies = 0;
  if (dsh) {
    dsh->incompletePatterns->Clear(false);
    delete (dsh);
    dsh = nil;
  }
  useHorizontalRep = false;
}

//_______________________________________________________________________

BaseRef _DataSet::makeDynamic(void) const {
  _DataSet *r = new _DataSet;
  r->theMap.Duplicate(&theMap);
  r->theFrequencies.Duplicate(&theFrequencies);
  if (theTT != &defaultTranslationTable) {
    r->theTT->AddAReference();
  }
  r->theNames.Duplicate(&theNames);
  r->streamThrough = streamThrough;
  // 20170507: SLKP TODO why do we need an additional reference here?
  // nInstances++;
  r->dsh = nil;
  r->useHorizontalRep = false;
  return r;
}

//_______________________________________________________________________

void _DataSet::ResetIHelper(void) {
  if (dsh && dsh->characterPositions.lLength == 256)
    for (long k = 0; k < 256; k++) {
      dsh->characterPositions.lData[k] = -1;
    }
}

//_______________________________________________________________________

void _DataSet::ConvertRepresentations(void) {
  if (useHorizontalRep == false) {
    _List horStrings;

    if (lLength == 0) {
      AppendNewInstance(new _Site);
    } else {
      _Site *aSite = (_Site *)lData[0];

      for (long str = 0; str < aSite->length(); str++) {
        horStrings < new _StringBuffer (DATA_SET_SWITCH_THRESHOLD);
      }

      for (long s = 0; s < lLength; s++) {
        _Site *aSite = (_Site *)lData[s];
        if (aSite->length() > horStrings.lLength || aSite->GetRefNo() != -1) {
          HandleApplicationError("Irrecoverable internal error in "
                                 "_DataSet::ConvertRepresentations. Sorry "
                                 "about that.",
                                 true);
          return;
        }

        for (long s2 = 0L; s2 < aSite->length(); s2++) {
          (*(_StringBuffer *)horStrings.lData[s2]) << aSite->get_char(s2);
        }
      }

      _List::Clear();
      theFrequencies.Clear();
      {
        for (long s = 0; s < horStrings.lLength; s++) {
          (*this) << horStrings(s);
        }
      }
    }
    useHorizontalRep = true;
  }
}

//_______________________________________________________________________

void _DataSet::AddSite(char c) {
  if (streamThrough) {
    if (theMap.lData[0] == 0) {
      if (theMap.lData[1] == 0) {
        if (theNames.lLength) {
          fprintf(streamThrough, ">%s\n", ((_String *)theNames(0))->get_str());
        } else {
          fprintf(streamThrough, ">Sequence 1\n");
        }
        AppendNewInstance(new _String(kEmptyString));
      }

      theMap.lData[1]++;
      theMap.lData[2]++;
      fputc(c, streamThrough);
    } else {
      HandleApplicationError("Can't add more sites to a file based data set, "
                             "when more that one sequence has been written!",
                             true);
    }
  } else {
    if (useHorizontalRep == false) {
      if (lLength < DATA_SET_SWITCH_THRESHOLD) {
        _Site *nC = new _Site(c);
        theFrequencies << 1L;
        AppendNewInstance(nC);
        return;
      } else {
        ConvertRepresentations();
      }
    }

    (*((_StringBuffer *)lData[0])) << c;

    /*long  f;

    if (!dsh)
    {
        checkPointer (dsh = new _DSHelper);
        for (f=0; f<256; f++)
            dsh->characterPositions << -1;
    }

    f = dsh->characterPositions.lData[c];

    if (f!=-1)
    {
        _Site* nC = new _Site(f);
        checkPointer(nC);
        theFrequencies[f]++;
        theFrequencies<<0;
        (*this)<<nC;
        nC->nInstances --;
    }
    else
    {
        dsh->characterPositions.lData[c] = lLength;*/
    //}
  }
}
//_______________________________________________________________________

void _DataSet::Write2Site(long index, char c) {
  if (streamThrough) {
    if (index == 0) {
      if (theMap.lData[2] == theMap.lData[1]) {
        theMap.lData[0]++;

        if (theNames.lLength > theMap.lData[0]) {
          fprintf(streamThrough, "\n>%s\n",
                  ((_String *)theNames(theMap.lData[0]))->get_str());
        } else {
          fprintf(streamThrough, "\n>Sequence %ld\n", theMap.lData[0] + 1);
        }

        theMap.lData[1] = 0;
      } else {
        HandleApplicationError("Can't write sequences of unequal lengths to a "
                               "file based data set.");
        return;
      }
    } else if (index != theMap.lData[1]) {
      HandleApplicationError("Can't write sites which are not consecutive to a "
                             "file based data set.");
      return;
    }

    theMap.lData[1]++;
    fputc(c, streamThrough);
  } else {
    /*if (!dsh)
    {
        WarnError ("Internal Error in 'Write2Site' - called Write2Site before
    any AddSite calls"); return;
    }*/

    if (useHorizontalRep) {
      long currentWritten = ((_String *)lData[0])->length();

      if (index >= currentWritten) {
        HandleApplicationError("Internal Error in 'Write2Site' - index is too "
                               "high (using compact representation)");
        return;
      } else {
        if (index == 0) {
          _StringBuffer *newString = new _StringBuffer(currentWritten);
          (*newString) << c;
          (*this) < newString;
        } else {
          long s = 1;
          for (; s < lLength; s++) {
            _StringBuffer *aString = (_StringBuffer *)lData[s];
            if (aString->length() == index) {
              (*aString) << c;
              break;
            }
          }
          if (s == lLength) {
            HandleApplicationError("Internal Error in 'Write2Site' - no "
                                   "appropriate  string to write too (compact "
                                   "representation)");
            return;
          }
        }
      }
    } else {
      if (index >= lLength) {
        HandleApplicationError(
            "Internal Error in 'Write2Site' - index is too high");
        return;
      }
      _Site *s = (_Site *)lData[index];
      long rN = s->GetRefNo();
      if (rN == -1) { // independent site
        // dsh->incompletePatterns->Delete (s,false);
        (*s) << c;
        // dsh->incompletePatterns->Insert (s,index);
      } else {
        _Site *ss = (_Site *)lData[rN];
        long sL = ss->length() - 1;
        if (ss->get_char(sL) != c) { // appending distinct char
          s->Duplicate(ss);
          s->set_char(sL, c);
          theFrequencies.lData[rN]--;

          rN = dsh->incompletePatterns->Find(s);
          if (rN >= 0) {
            rN = dsh->incompletePatterns->GetXtra(rN);
            /*_Site* s2 = (_Site*)lData[rN];
            if (s2->GetRefNo() != -1 || !s->Equal(s2))
            {
                WarnError ("Mapping Error");
            }*/
            theFrequencies[rN]++;
            s->Clear();
            s->SetRefNo(rN);
          } else {
            theFrequencies[index]++;
            s->SetRefNo(-1);
            dsh->incompletePatterns->Insert(s, index);
          }
        }
      }
    }
  }
}

//_______________________________________________________________________

void _DataSet::CheckMapping(long index) {
  if (index >= lLength) {
    HandleApplicationError(
        "Internal Error in 'CheckMapping' - index is too high", true);
  }

  _Site *s = (_Site *)lData[index];

  for (long k = 0L; k < index; k++) {
    _Site *ss = (_Site *)lData[k];
    if (ss->GetRefNo() == -1) {
      if (s->Equal(ss)) {
        theFrequencies[index]--;
        theFrequencies[k]++;
        s->Clear();
        s->SetRefNo(k);
      }
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
  if (theTT && (theTT != &defaultTranslationTable)) {
    DeleteObject(theTT);
  }
  theTT = (_TranslationTable *)newTT->theTT->makeDynamic();
}

//_______________________________________________________________________

void _DataSet::SetTranslationTable(_TranslationTable *newTT) {
  if (theTT && (theTT != &defaultTranslationTable)) {
    DeleteObject(theTT);
  }
  theTT = (_TranslationTable *)newTT->makeDynamic();
}
//_______________________________________________________________________
void _DataSet::Finalize(void) {
  if (streamThrough) {
    fclose(streamThrough);
    streamThrough = nil;
    theMap.Clear();
  } else {
    if (useHorizontalRep) {
      bool good = true;
      for (long s = 0; s < lLength; s++) {
        good = good &&
               ((_String *)lData[0])->length() == ((_String *)lData[s])->length();
      }

      if (!good) {
        Clear();
        HandleApplicationError("Internal Error in _DataSet::Finalize. Unequal "
                               "sequence lengths in compact representation",
                               true);
        return;
      }

      _List dups;
      _List uniquePats;
      _AVLListX dupsAVL(&dups);

      long siteCounter = ((_String *)lData[0])->length();

      for (long i1 = 0L; i1 < siteCounter; i1++) {
        _Site *tC = new _Site();

        for (long i2 = 0L; i2 < lLength; i2++) {
          (*tC) << ((_String *)lData[i2])->get_char(i1);
        }

        long ff = dupsAVL.Find(tC);
        if (ff < 0) {
          uniquePats << tC;
          dupsAVL.Insert(tC, theFrequencies.lLength);
          theMap << theFrequencies.lLength;
          theFrequencies << 1;
        } else {
          ff = dupsAVL.GetXtra(ff);
          theMap << ff;
          theFrequencies.lData[ff]++;
        }

        DeleteObject(tC);
      }
      dupsAVL.Clear(false);
      _List::Clear();
      _List::Duplicate(&uniquePats);
    } else {
      long j, k;

      _Site *tC;
      {
        _List dups;
        _AVLListX dupsAVL(&dups);

        for (long i1 = 0; i1 < lLength; i1++) {
          tC = (_Site *)lData[i1];
          long ff = dupsAVL.Find(tC);
          if (ff < 0) {
            dupsAVL.Insert(tC, i1);
          } else {
            ff = dupsAVL.GetXtra(ff);
            tC->Clear();
            tC->SetRefNo(ff);
            theFrequencies.lData[ff]++;
          }
        }
        dupsAVL.Clear(false);
      }

      _SimpleList refs(lLength), toDelete(lLength);
      j = 0;

      for (long i1 = 0; i1 < lLength; i1++) {
        tC = (_Site *)(*(_List *)this)(i1);
        k = tC->GetRefNo();
        if (k == -1) {
          refs << j++;
        } else {
          toDelete << i1;
          refs << -1;
        }
      }

      for (long i2 = 0; i2 < lLength; i2++) {
        tC = (_Site *)(*(_List *)this)(i2);
        k = tC->GetRefNo();
        if (k >= 0) {
          j = refs.lData[k];
          if (j < 0) {
            HandleApplicationError(kErrorStringDatasetRefIndexError);
          } else {
            refs.lData[i2] = j;
          }
        }
      }

      theMap.Clear();
      theMap.Duplicate(&refs);
      DeleteList(toDelete);
      theFrequencies.DeleteList(toDelete);

      for (long i3 = 0; i3 < lLength; i3++) {
        tC = (_Site *)(*(_List *)this)(i3);
        tC->SetRefNo(0);
      }
      if (dsh) {
        dsh->incompletePatterns->Clear(false);
        delete (dsh);
        dsh = nil;
      }
    }
  }
}
//_______________________________________________________________________
void _DataSet::Compact(long index) {
  if (useHorizontalRep) {
    HandleApplicationError(
        "Internal Error: _DataSet::Compact called with compact represntation",
        true);
    return;
  }

  _Site *tC = (_Site *)(*(_List *)this)(index);
  if (tC->GetRefNo() != -1)
  // take care of double referencing
  {
    _Site *tCC = tC;
    long lastRef, count = 0;
    do {
      lastRef = tCC->GetRefNo();
      count++;
      tCC = (_Site *)(*(_List *)this)(tCC->GetRefNo());
    } while (tCC->GetRefNo() != -1);
    if (count > 1) {
      theFrequencies[lastRef]++;
    }

    tC->SetRefNo(lastRef);
  }
  /*if (tC->GetRefNo()==-1)
  {
   long f = dsh->incompletePatterns->Find(tC);
   if (f >= 0)
   {
          f = dsh->incompletePatterns->GetXtra (f);
          if (f<index)
          {
          theFrequencies[f]++;
          tC->Clear();
          tC->SetRefNo(f);
      }
      else
          tC->Finalize();
   }
  }*/
}

//_______________________________________________________________________
inline char _DataSet::operator()(unsigned long site, unsigned long pos,
                                 unsigned int) const {
  return (((_String **)lData)[theMap.lData[site]])->get_char(pos);
}

//_________________________________________________________
long _DataSet::ComputeSize(void) {
  long res = sizeof(_DataSet);

  res += (theMap.lLength + lLength + theFrequencies.lLength) * sizeof(long);
  res += lLength * sizeof(_Site);

  for (long i = 0; i < lLength; i++) {
    res += ((_Site *)(*(_List *)this)(i))->length();
  }

  return res;
}

//_________________________________________________________
hyFloat _DataSet::CheckAlphabetConsistency(void) {
  long charsIn = 0, gaps = 0, total = 0;

  bool checks[256];

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

  for (; charsIn < 256; charsIn++) {
    checks[charsIn] = false;
  }

  for (charsIn = 0; charsIn < baseSymbols.length(); charsIn++) {
    checks[(unsigned char)baseSymbols.get_char(charsIn)] = true;
  }

  charsIn = 0;

  for (long i = 0; i < lLength; i++) {
    _String *thisColumn = (_String *)lData[i];
    long w = theFrequencies.lData[i];
    for (long j = 0; j < thisColumn->length(); j++)
      if (checks[(unsigned char)thisColumn->get_char(j)]) {
        charsIn += w;
      } else if (gapChar == thisColumn->get_char(j)) {
        gaps += w;
      }

    total += w * thisColumn->length();
  }

  return (hyFloat)charsIn / (total - gaps + 1.);
}

//___________________________________________________

BaseRef _DataSet::toStr(unsigned long) {
  _StringBuffer *s = new _StringBuffer(NoOfSpecies() * 30), *str;

  (*s) << _String((long)NoOfSpecies()) << " species:";

  (*s) << (_String *)theNames.toStr();

  (*s) << ";\nTotal Sites:" << _String((long)GetNoTypes())
       << ";\nDistinct Sites:" << _String((long)theFrequencies.lLength);

  return s;
}

//___________________________________________________

void _DataSet::toFileStr(FILE *dest, unsigned long padding) {
  fprintf(dest, "%ld species: ", NoOfSpecies());
  theNames.toFileStr(dest, padding);

  fprintf(dest, ";\nTotal Sites: %ld", GetNoTypes());
  fprintf(dest, ";\nDistinct Sites: %ld", theFrequencies.lLength);

  /*  fprintf (dest,"\n");
      for (long j=0; j<noOfSpecies;j++)
      {
          fprintf (dest,"\n");
          for (long i=0; i<theMap.lLength; i++)
          {
              fprintf (dest,"%c",(*this)(i,j,1));
          }
      }*/
}

//_________________________________________________________

void _DataSet::AddName(_String const &s) {
  theNames.AppendNewInstance(
      new _String(s, 0, s.FirstNonSpaceIndex(0, -1, kStringDirectionForward)));
}

//_________________________________________________________

void _DataSet::InsertName(_String const &name, long where) {
  theNames.InsertElement(new _String(name), where, false);
}

//_________________________________________________________

void _DataSet::MatchIndices(_Formula &f, _SimpleList &receptacle, bool isVert,
                            long limit, _String const *scope) const {
  _String varName = isVert ? "siteIndex" : "speciesIndex";
  varName = AppendContainerName(varName, scope);
  _Variable *v = CheckReceptacle(&varName, kEmptyString, false);

  // fprintf (stderr, "\n_DataSet::MatchIndices %d %s [%s] %s\n", isVert, scope
  // ? scope->sData : "none", varName.sData, ((_String*)f.toStr())->sData);

  for (long i = 0L; i < limit; i++) {
    v->SetValue(new _Constant((hyFloat)i), nil);
    HBLObjectRef res = f.Compute();
    // fprintf (stderr, "%ld %g\n", i, res->Compute()->Value());
    if (res && !CheckEqual(res->Value(), 0.0)) {
      receptacle << i;
    }
  }
  v->SetValue(new _Constant(0.0), nil);
}

//_________________________________________________________

_TranslationTable *_DataSet::CheckCompatibility(_SimpleList const &ref,
                                                char concatOrCombine) {
  _DataSet *currentSet = (_DataSet *)dataSetList(ref.Element(0));

  _TranslationTable *theEnd = new _TranslationTable(*(currentSet->theTT));

  long refNo =
      concatOrCombine ? currentSet->NoOfSpecies() : currentSet->NoOfColumns();
  char emptyStringChar = theEnd->GetSkipChar();

  for (long k = 1; k < ref.lLength; k++) {
    currentSet = (_DataSet *)dataSetList(ref.Element(k));

    _TranslationTable *tryMe = theEnd->MergeTables(currentSet->theTT);

    if (tryMe) {
      if (emptyStringChar) {
        DeleteObject(theEnd);
        theEnd = tryMe;
        continue;
      } else {
        if ((concatOrCombine && (currentSet->NoOfSpecies() == refNo)) ||
            (!concatOrCombine && (currentSet->NoOfColumns() == refNo))) {
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
    for (long i = 0; i < k; i++) {
      if (k) {
        warningMessage = warningMessage & ", ";
      }
      warningMessage = warningMessage &
                       ((_String *)dataSetNamesList(ref.Element(k)))->Enquote();
    }
    HandleApplicationError(warningMessage);
    DeleteObject(tryMe);
    DeleteObject(theEnd);
    return nil;
  }

  return theEnd;
}

//_________________________________________________________

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

  long maxSpecies = 0, maxDataSet = 0, siteIndex;

  _DataSet *currentSet;

  char emptyStringSlot = jointTable->GetSkipChar();

  for (long i = 0; i < ref.lLength; i++) {
    currentSet = (_DataSet *)dataSetList(ref(i));

    long specCount = currentSet->NoOfSpecies(),
         siteCount = currentSet->NoOfColumns();

    if (specCount > maxSpecies) {
      maxSpecies = specCount;
      maxDataSet = i;
    }
    for (long j = 0; j < siteCount; j++) {
      bigDataSet->AddSite((*currentSet)(j, 0, 1));
    }
  }

  for (long k = 1; k < maxSpecies; k++) {
    siteIndex = 0;
    for (long i = 0; i < ref.lLength; i++) {
      currentSet = (_DataSet *)dataSetList(ref.lData[i]);

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

  currentSet = (_DataSet *)dataSetList(ref(maxDataSet));
  for (long i = 0L; i < maxSpecies; i++) {
    bigDataSet->AddName(*currentSet->GetSequenceName(i));
  }

  bigDataSet->Finalize();
  bigDataSet->SetNoSpecies(maxSpecies);
  return bigDataSet;
}

//_________________________________________________________

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

_String*        _DataSet::GetSequenceCharacters (long seqID)  const{
  
  unsigned long        upTo = NoOfColumns();
  _StringBuffer * aSequence = new _StringBuffer (upTo);
  
  if (seqID >= 0L && seqID < noOfSpecies) {
    for (unsigned long k2=0UL; k2<upTo; k2++) {
      (*aSequence) << GetSite (k2)->char_at(seqID);
    }
  }
  aSequence->TrimSpace ();
  return aSequence;
}

//_________________________________________________________

bool    StoreADataSet (_DataSet* ds, _String* setName) {
    if (!setName->IsValidIdentifier (true)) {
        HandleApplicationError (setName->Enquote() & " is not a valid identifier while constructing a DataSet");
        return false;
    }
    
    long type = HY_BL_DATASET, index;
    _DataSet * existing_ds = (_DataSet * )_HYRetrieveBLObjectByNameMutable (*setName, type, &index, false, false);
    
    
    if (! existing_ds) {
        dataSetNamesList << setName;
        dataSetList < ds;
    } else {
        
        
        bool isDifferent = existing_ds->NoOfSpecies () != ds->NoOfSpecies() ||
        existing_ds->NoOfColumns () != ds->NoOfColumns() ||
        existing_ds->NoOfUniqueColumns () != ds->NoOfUniqueColumns() ||
        existing_ds->GetTT () != ds->GetTT();
        
        
        
        for (AVLListXLIteratorKeyValue filter_key_value : ObjectIndexer (HY_BL_DATASET_FILTER)) {
            _DataSetFilter * filter = (_DataSetFilter*) filter_key_value.get_object();
            if (filter->GetData() == existing_ds) {
                if (isDifferent) {
                    ReportWarning (_String("Overwriting dataset '") & *setName & "' caused DataSetFilter " & GetFilterName(filter_key_value.get_index())->Enquote('\'') & " to be deleted");
                    DeleteDataFilter(filter_key_value.get_index());
                } else {
                    filter->SetData(ds);
                }
            }
        }
        
        dataSetList.Replace(index,ds,false);
    }
    
    CheckReceptacleAndStore (*setName&".mapping",kEmptyString,false, new _MathObject, false);

    if (hy_env::EnvVariableTrue(hy_env::normalize_sequence_names)) {
        _List _id_mapping;
        _AVLListXL id_mapping (&_id_mapping);
        bool       did_something = false;
        
        for (unsigned long i = 0UL; i < ds->NoOfSpecies(); i ++) {
            _String * old_name = new _String (*ds->GetSequenceName (i));
            if (! old_name->IsValidIdentifier(false) ) {
                ds->GetSequenceName (i)->ConvertToAnIdent(false);
                did_something = true;
            }
            if (id_mapping.Find (ds->GetSequenceName (i)) >= 0) {
                _String new_name (*ds->GetSequenceName (i));
                long suffix = 1L;
                do {
                    new_name = *ds->GetSequenceName (i) & "_" & suffix++;
                } while (id_mapping.Find (&new_name) >= 0);
                *ds->GetSequenceName (i) = new_name;
                did_something = true;
            }
            
            ds->GetSequenceName (i)->AddAReference();
            id_mapping.Insert (ds->GetSequenceName (i), (long)old_name, false, false);
        }
        
        if (did_something) {
            _AssociativeList * mapping = new _AssociativeList();
            
            _SimpleList history;
            long t,
            current_index = id_mapping.Traverser(history, t, id_mapping.GetRoot());
            
            while (current_index >= 0L) {
                mapping->MStore(*(_String*)_id_mapping.GetItem (current_index), *(_String*)id_mapping.GetXtra(current_index));
                current_index = id_mapping.Traverser(history, t);
            }
            
            CheckReceptacleAndStore (*setName&".mapping",kEmptyString,false, mapping, false);
        }
    }
    
    CheckReceptacleAndStore (*setName&".species",kEmptyString,false, new _Constant (ds->NoOfSpecies()), false);
    CheckReceptacleAndStore (*setName&".sites",kEmptyString,false, new _Constant (ds->NoOfColumns()), false);
    CheckReceptacleAndStore (*setName&".unique_sites",kEmptyString,false, new _Constant (ds->NoOfUniqueColumns()), false);
    
    return true;
}
