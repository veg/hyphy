/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2009
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon              (apoon@cfenet.ubc.ca)

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

#include "site.h"
#include "dataset.h"
#include "batchlan.h"
#include "ctype.h"

_DataSet::_DataSet(void) {
  theTT = &defaultTranslationTable;
  streamThrough = nil;
  dsh = nil;
  useHorizontalRep = false;
}

//______________________________________________________________________________
_DataSet::_DataSet(long l)
    : _List((unsigned long) l),
      theFrequencies((unsigned long)
                     l) // with estimated number of sites per file
      {
  dsh = nil;
  streamThrough = nil;
  theTT = &defaultTranslationTable;
  useHorizontalRep = false;
}

//______________________________________________________________________________
_DataSet::_DataSet(FILE *f) {
  dsh = nil;
  useHorizontalRep = false;
  theTT = &defaultTranslationTable;
  streamThrough = f;
  theMap << 0; // current sequence
  theMap << 0; // current site
  theMap << 0; // total sites
}

//______________________________________________________________________________
_DataSet::~_DataSet(void) {
  if (theTT != &defaultTranslationTable) {
    DeleteObject(theTT);
  }
}

//______________________________________________________________________________
void _DataSet::Clear(void) {

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

//______________________________________________________________________________
BaseRef _DataSet::makeDynamic(void) {
  _DataSet *r = new _DataSet;
  checkPointer(r);
  memcpy((char *)r, (char *)this, sizeof(_DataSet));
  r->nInstances = 1;
  r->theMap.Duplicate(&theMap);
  r->theFrequencies.Duplicate(&theFrequencies);
  if (theTT != &defaultTranslationTable) {
    r->theTT->nInstances++;
  }
  r->theNames.Duplicate(&theNames);
  r->streamThrough = streamThrough;
  nInstances++;
  r->dsh = nil;
  r->useHorizontalRep = false;
  return r;
}

//______________________________________________________________________________
void _DataSet::ResetIHelper(void) {
  if (dsh && dsh->characterPositions.lLength == 256)
    for (long k = 0; k < 256; k++) {
      dsh->characterPositions.lData[k] = -1;
    }
}

//______________________________________________________________________________
void _DataSet::ConvertRepresentations(void) {
  if (useHorizontalRep == false) {
    _List horStrings;

    if (lLength == 0) {
      AppendNewInstance(new _Site);
    } else {
      _Site *aSite = (_Site *)lData[0];

      for (long str = 0; str < aSite->sLength; str++) {
        _String *aString = new _String(DATA_SET_SWITCH_THRESHOLD, true);
        horStrings << aString;
        aString->nInstances--;
      }

      for (long s = 0; s < lLength; s++) {
        _Site *aSite = (_Site *)lData[s];
        if (aSite->sLength > horStrings.lLength || aSite->GetRefNo() != -1) {
          FlagError("Irrecoverable internal error in "
                    "_DataSet::ConvertRepresentations. Sorry about that.");
          return;
        }
        aSite->Finalize();
        for (long s2 = 0; s2 < aSite->sLength; s2++) {
          (*(_String *)horStrings.lData[s2]) << aSite->sData[s2];
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

//______________________________________________________________________________
void _DataSet::AddSite(char c) {
  if (streamThrough) {
    if (theMap.lData[0] == 0) {
      if (theMap.lData[1] == 0) {
        if (theNames.lLength) {
          fprintf(streamThrough, ">%s\n", ((_String *)theNames(0))->getStr());
        } else {
          fprintf(streamThrough, ">Sequence 1\n");
        }
        (*this) && &empty;
      }

      theMap.lData[1]++;
      theMap.lData[2]++;
      fputc(c, streamThrough);
    } else {
      WarnError("Can't add more sites to a file based data set, when more that "
                "one sequence has been written!");
    }
  } else {
    if (useHorizontalRep == false) {
      if (lLength < DATA_SET_SWITCH_THRESHOLD) {
        _Site *nC = new _Site(c);
        checkPointer(nC);
        theFrequencies << 1;
        (*this) << nC;
        nC->nInstances--;
        return;
      } else {
        ConvertRepresentations();
      }
    }

    (*((_String *)lData[0])) << c;

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

//______________________________________________________________________________
void _DataSet::Write2Site(long index, char c) {
  if (streamThrough) {
    if (index == 0) {
      if (theMap.lData[2] == theMap.lData[1]) {
        theMap.lData[0]++;

        if (theNames.lLength > theMap.lData[0]) {
          fprintf(streamThrough, "\n>%s\n",
                  ((_String *)theNames(theMap.lData[0]))->getStr());
        } else {
          fprintf(streamThrough, "\n>Sequence %ld\n", theMap.lData[0] + 1);
        }

        theMap.lData[1] = 0;
      } else {
        WarnError("Can't write sequences of unequal lengths to a file based "
                  "data set.");
        return;
      }
    } else if (index != theMap.lData[1]) {
      WarnError("Can't write sites which are not consecutive to a file based "
                "data set.");
      return;
    }

    theMap.lData[1]++;
    fputc(c, streamThrough);
  } else {
    /*if (!dsh)
    {
        WarnError ("Internal Error in 'Write2Site' - called Write2Site before
    any AddSite calls");
        return;
    }*/

    if (useHorizontalRep) {
      long currentWritten = ((_String *)lData[0])->sLength;

      if (index >= currentWritten) {
        WarnError("Internal Error in 'Write2Site' - index is too high (using "
                  "compact representation)");
        return;
      } else {
        if (index == 0) {
          _String *newString = new _String(currentWritten, true);
          (*newString) << c;
          (*this) << newString;
          newString->nInstances--;
        } else {
          long s = 1;
          for (; s < lLength; s++) {
            _String *aString = (_String *)lData[s];
            if (aString->sLength == index) {
              (*aString) << c;
              break;
            }
          }
          if (s == lLength) {
            WarnError("Internal Error in 'Write2Site' - no appropriate  string "
                      "to write too (compact representation)");
            return;
          }
        }
      }
    } else {
      if (index >= lLength) {
        WarnError("Internal Error in 'Write2Site' - index is too high");
        return;
      }
      _Site *s = (_Site *)lData[index];
      long rN = s->GetRefNo();
      if (rN == -1) { // independent site
                      //dsh->incompletePatterns->Delete (s,false);
        (*s) << c;
        //dsh->incompletePatterns->Insert (s,index);
      } else {
        _Site *ss = (_Site *)lData[rN];
        long sL = ss->sLength - 1;
        if (ss->sData[sL] != c) { // appending distinct char
          s->Duplicate(ss);
          s->sData[sL] = c;
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

//______________________________________________________________________________
void _DataSet::CheckMapping(long index) {
  if (index >= lLength) {
    FlagError("Internal Error in 'CheckMapping' - index is too high");
  }

  _Site *s = (_Site *)lData[index];

  for (long k = 0; k < index; k++) {
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

//______________________________________________________________________________
long _DataSet::GetCharDimension(void) {
  // return the size of the alphabet space
  return theTT->LengthOfAlphabet();
}

//______________________________________________________________________________
long _DataSet::GetNoTypes(void) // return the number of unique columns
    {
  return theMap.countitems();
}

//______________________________________________________________________________
long _DataSet::GetFreqType(long index) // return the frequency of a site
    {
  return theFrequencies(theMap(index));
}

//______________________________________________________________________________
void _DataSet::SetTranslationTable(_DataSet *newTT) {
  if (theTT && (theTT != &defaultTranslationTable)) {
    DeleteObject(theTT);
  }
  theTT = (_TranslationTable *)newTT->theTT->makeDynamic();
}

//______________________________________________________________________________
void _DataSet::SetTranslationTable(_TranslationTable *newTT) {
  if (theTT && (theTT != &defaultTranslationTable)) {
    DeleteObject(theTT);
  }
  theTT = (_TranslationTable *)newTT->makeDynamic();

}

//______________________________________________________________________________
void _DataSet::Finalize(void) {
  if (streamThrough) {
    fclose(streamThrough);
    streamThrough = nil;
    theMap.Clear();
  } else {
    if (useHorizontalRep) {
      bool good = true;
      for (long s = 0; s < lLength; s++) {
        ((_String *)lData[s])->Finalize();
        good = good &&
               ((_String *)lData[0])->sLength == ((_String *)lData[s])->sLength;
      }

      if (!good) {
        Clear();
        WarnError("Internal Error in _DataSet::Finalize. Unequal sequence "
                  "lengths in compact representation");
        return;
      }

      _List dups;
      _List uniquePats;
      _AVLListX dupsAVL(&dups);

      long siteCounter = ((_String *)lData[0])->sLength;

      for (long i1 = 0; i1 < siteCounter; i1++) {
        _Site *tC = new _Site();
        checkPointer(tC);

        for (long i2 = 0; i2 < lLength; i2++) {
          (*tC) << ((_String *)lData[i2])->sData[i1];
        }

        tC->Finalize();

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
            warnError(-171);
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
        tC->Finalize();
      }

      if (dsh) {
        dsh->incompletePatterns->Clear(false);
        delete (dsh);
        dsh = nil;
      }
    }
  }
}

//______________________________________________________________________________
void _DataSet::Compact(long index) {
  if (useHorizontalRep) {
    WarnError(
        "Internal Error: _DataSet::Compact called with compact represntation");
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

//______________________________________________________________________________
inline char _DataSet::operator()(unsigned long site, unsigned long pos,
                                 unsigned int) {
  return (((_String **)lData)[theMap.lData[site]])->sData[pos];
}

//______________________________________________________________________________
long _DataSet::ComputeSize(void) {
  long res = sizeof(_DataSet);

  res += (theMap.lLength + lLength + theFrequencies.lLength) * sizeof(long);
  res += lLength * sizeof(_Site);

  for (long i = 0; i < lLength; i++) {
    res += ((_Site *)(*(_List *)this)(i))->sLength;
  }

  return res;
}

//______________________________________________________________________________
_Parameter _DataSet::CheckAlphabetConsistency(void) {
  long charsIn = 0, gaps = 0, total = 0;

  char checks[256], gapChar = theTT->GetGapChar();

  _String baseSymbols(16L, true);
  baseSymbols << theTT->RetrieveCharacters();
  if (theTT->DetectType() == HY_TRANSLATION_TABLE_STANDARD_NUCLEOTIDE) {
    baseSymbols << 'U';
  }
  baseSymbols.Finalize();

  memset(checks, 0, 256);

  for (long c = 0L; c < baseSymbols.sLength; c++) {
    checks[baseSymbols.sData[c]] = 1;
  }

  charsIn = 0;

  for (unsigned long i = 0L; i < lLength; i++) {
    _String *thisColumn = (_String *)lData[i];
    long w = theFrequencies.lData[i];
    for (long j = 0; j < thisColumn->sLength; j++)
      if (checks[thisColumn->sData[j]]) {
        charsIn += w;
      } else if (gapChar == thisColumn->sData[j]) {
        gaps += w;
      }

    total += w * thisColumn->sLength;
  }

  return (_Parameter) charsIn / (total - gaps + 1.);

}

//______________________________________________________________________________
BaseRef _DataSet::toStr(void) {
  _String *s = new _String(NoOfSpecies() * 30, true), *str;

  checkPointer(s);
  (*s) << _String((long) NoOfSpecies());
  (*s) << " species:";

  str = (_String *)theNames.toStr();
  (*s) << *str;
  DeleteObject(str);

  (*s) << ";\nTotal Sites:";
  (*s) << _String((long) GetNoTypes());
  (*s) << ";\nDistinct Sites:";
  (*s) << _String((long) theFrequencies.lLength);

  s->Finalize();

  return s;
}

//______________________________________________________________________________
void _DataSet::toFileStr(FILE *dest) {
  fprintf(dest, "%ld species: ", NoOfSpecies());
  theNames.toFileStr(dest);

  fprintf(dest, ";\nTotal Sites: %ld", GetNoTypes());
  fprintf(dest, ";\nDistinct Sites: %ld", theFrequencies.lLength);
}

//______________________________________________________________________________
void _DataSet::AddName(_String &s) {
  theNames.AppendNewInstance(
      new _String(s, 0, s.FirstNonSpaceIndex(0, -1, -1)));
}

//______________________________________________________________________________
void _DataSet::MatchIndices(_Formula &f, _SimpleList &receptacle, bool isVert,
                            long limit) {
  _String varName = isVert ? "siteIndex" : "speciesIndex";
  _Variable *v = CheckReceptacle(&varName, empty, false);

  for (long i = 0; i < limit; i++) {
    v->SetValue(new _Constant(i), nil);
    _PMathObj res = f.Compute();
    if (res && !CheckEqual(res->Value(), 0.0)) {
      receptacle << i;
    }
  }
  v->SetValue(new _Constant(0.0), nil);
}

//______________________________________________________________________________
void _DataSet::FindAllSitesLikeThisOne(long index, _SimpleList &receptacle) {
  if (index >= 0 && index < theMap.lLength) {
    index = theMap.lData[index];
    for (long k = 0; k < theMap.lLength; k++)
      if (theMap.lData[k] == index) {
        receptacle << k;
      }
  }
}

//______________________________________________________________________________
_TranslationTable *_DataSet::CheckCompatibility(_SimpleList &ref,
                                                char concatOrCombine) {

  _DataSet *currentSet = (_DataSet *)dataSetList(ref(0));

  _TranslationTable *theEnd = new _TranslationTable(*(currentSet->theTT));
  checkPointer(theEnd);
  long refNo =
      concatOrCombine ? currentSet->NoOfSpecies() : currentSet->NoOfColumns();

  char emptyChar = theEnd->GetSkipChar();

  for (long k = 1; k < ref.lLength; k++) {

    currentSet = (_DataSet *)dataSetList(ref(k));
    _TranslationTable *tryMe = theEnd->MergeTables(currentSet->theTT);

    if (tryMe) {
      if (emptyChar) {
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
    _String warningMessage("The data set:");
    warningMessage =
        warningMessage & *((_String *)dataSetNamesList(ref(k))) &
        _String(" was found incompatible with one of the following data sets:");
    for (long i = 0; i < k; i++) {
      warningMessage = warningMessage & *((_String *)dataSetNamesList(ref(i))) &
                       _String(",");
    }
    warningMessage =
        warningMessage &
        _String(" and was dropped from the dataset merging operation");
    ReportWarning(warningMessage);
    ref.Delete(k);
    k--;
  }

  return theEnd;
}

//______________________________________________________________________________
// concatenates (adds columns together) several datasets
// in case the number of species in the datasets are different the
// deficiencies will be padded
// by omission symbols
// in case translation tables are different, they will be merged, provided
// it can be done,
// otherwise the incompatible datasets will be ignored during this
// operation.
_DataSet *_DataSet::Concatenate(_SimpleList ref) {

  _TranslationTable *jointTable;

  jointTable = CheckCompatibility(ref, 1);

  _DataSet *bigDataSet = new _DataSet;
  checkPointer(bigDataSet);

  bigDataSet->theTT = jointTable;

  // pass one - determine the max max number of species present and what dataset
  // are they coming from

  long maxSpecies = 0, maxDataSet = 0, siteIndex;

  _DataSet *currentSet;

  char emptySlot = jointTable->GetSkipChar();

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

      if (cns <= k) {
        for (long j = 0; j < cnc; j++, siteIndex++) {
          bigDataSet->Write2Site(siteIndex, emptySlot);
        }
      } else {
        for (long j = 0; j < cnc; j++, siteIndex++) {
          bigDataSet->Write2Site(siteIndex, (*currentSet)(j, k, 1));
        }
      }
    }
  }

  currentSet = (_DataSet *)dataSetList(ref(maxDataSet));
  {
    for (long i = 0; i < maxSpecies; i++) {
      bigDataSet->AddName(*((_String *)(currentSet->GetNames())(i)));
    }
  }

  bigDataSet->Finalize();
  bigDataSet->SetNoSpecies(maxSpecies);
  return bigDataSet;
}

//______________________________________________________________________________
_DataSet *_DataSet::Combine(_SimpleList ref) {

  // combines (adds rows together) several datasets
  // in case the number of species in the datasets are different the
  // deficiencies will be padded
  // by omission symbols
  // in case translation tables are different, they will be merged, provided
  // it can be done,
  // otherwise the incompatible datasets will be ignored during this
  // operation.

  _TranslationTable *jointTable;

  jointTable = CheckCompatibility(ref, 0);

  _DataSet *bigDataSet = new _DataSet;
  checkPointer(bigDataSet);
  bigDataSet->theTT = jointTable;

  // pass one - determine the max max number of sites present and what dataset
  // are they coming from

  long i, j, k, maxSites = 0, sitesAvail, nsc = 0;

  _DataSet *currentSet;

  char emptySlot = jointTable->GetSkipChar();

  for (i = 0; i < ref.lLength; i++) {
    currentSet = (_DataSet *)dataSetList(ref(i));
    if (currentSet->NoOfColumns() > maxSites) {
      maxSites = currentSet->NoOfColumns();
    }
    nsc += currentSet->NoOfSpecies();
  }

  for (k = 0; k < ref.lLength; k++) {
    currentSet = (_DataSet *)dataSetList(ref(k));
    sitesAvail = currentSet->NoOfColumns();

    long cns = currentSet->NoOfSpecies();
    for (i = 0; i < cns; i++) {
      bigDataSet->AddName(*((_String *)(currentSet->GetNames())(i)));

      if (!(k || i)) {
        for (j = 0; j < sitesAvail; j++) {
          bigDataSet->AddSite((*currentSet)(j, 0, 1));
        }

        for (; j < maxSites; j++) {
          bigDataSet->AddSite(emptySlot);
        }

      } else {
        for (j = 0; j < sitesAvail; j++) {
          bigDataSet->Write2Site(j, (*currentSet)(j, i, 1));
        }
        for (; j < maxSites; j++) {
          bigDataSet->Write2Site(j, emptySlot);
        }
      }
    }
  }

  bigDataSet->Finalize();
  bigDataSet->SetNoSpecies(nsc);
  return bigDataSet;
}

//______________________________________________________________________________
//  20110610: SLKP, some cleanup and refactoring
void _DataSet::ProcessPartition(_String &input2, _SimpleList &target,
                                bool isVertical, _SimpleList *additionalFilter,
                                _SimpleList *otherDimension) {
  if (!input2.sLength) {
    return;
  }

  // decide if the input is an enumeration or a formula
  long totalLength;

  if (additionalFilter) {
    totalLength = additionalFilter->lLength;
  } else {
    totalLength = isVertical ? theMap.lLength : noOfSpecies;
  }

  _String input(input2);

  if (!input.IsALiteralArgument(true)) { // not a literal argument
    _Formula fmla, lhs;

    _FormulaParsingContext fpc;
    long outcome = Parse(&fmla, input, fpc, &lhs);

    if (outcome != HY_FORMULA_EXPRESSION) {
      WarnError(input & _String(" is an invalid partition specification"));
      return;
    }
    _PMathObj fV = fmla.Compute();
    if (fV && fV->ObjectClass() == STRING) {
      _String newSpec(128L, true);
      newSpec << '"';
      newSpec << ((_FString *)fV)->theString;
      newSpec << '"';
      newSpec.Finalize();
      ProcessPartition(newSpec, target, isVertical, additionalFilter);
    } else {
      _DataSet::MatchIndices(fmla, target, isVertical, totalLength);
    }
  } else { // an explicit enumeration or a regular expression

    // a regular expression
    if (input.getChar(0) == '/' && input.getChar(input.sLength - 1) == '/') {
      input.Trim(1, input.sLength - 2);
      int errCode;
      Ptr regex = PrepRegExp(&input, errCode, true);
      if (errCode) {
        WarnError(GetRegExpError(errCode));
        return;
      }

      // now set do the matching
      // using only the sites that are specced in the additionalFilter

      if (!isVertical) {
        _SimpleList *eligibleSeqs;

        if (additionalFilter) {
          eligibleSeqs = additionalFilter;
        } else {
          eligibleSeqs = new _SimpleList(0, totalLength, 1);
        }

        _SimpleList matches;
        for (long specCount = 0; specCount < eligibleSeqs->lLength;
             specCount++) {
          _String pattern(theMap.lLength, false);
          long seqPos = eligibleSeqs->lData[specCount];

          if (otherDimension)
            for (long seqSlider = 0; seqSlider < otherDimension->lLength; seqSlider++) {
              pattern.sData[seqSlider] =
                  GetSite(otherDimension->lData[seqSlider])->sData[seqPos];
            }
          else
            for (long seqSlider = 0; seqSlider < theMap.lLength; seqSlider++) {
              pattern.sData[seqSlider] = GetSite(seqSlider)->sData[seqPos];
            }

          matches.Clear();
          pattern.RegExpMatch(regex, matches);
          if (matches.lLength) {
            target << specCount;
          }
        }

        if (eligibleSeqs != additionalFilter) {
          DeleteObject(eligibleSeqs);
        }
      } else {
        bool *eligibleMarks = new bool[lLength];
        checkPointer(eligibleMarks);

        for (long fillerID = 0; fillerID < lLength; fillerID++) {
          eligibleMarks[fillerID] = false;
        }

        if (additionalFilter) {
          for (long siteIndex = 0; siteIndex < additionalFilter->lLength; siteIndex++) {
            eligibleMarks[theMap.lData[additionalFilter->lData[siteIndex]]] =
                true;
          }
        } else {
          for (long siteIndex = 0; siteIndex < lLength; siteIndex++) {
            eligibleMarks[siteIndex] = true;
          }
        }
        _SimpleList matches;
        _String *tempString = nil;
        if (otherDimension) {
          tempString = new _String(otherDimension->lLength, false);
        }

        for (long siteCounter = 0; siteCounter < lLength; siteCounter++)
          if (eligibleMarks[siteCounter]) {
            matches.Clear();
            if (otherDimension) {
              _Site *aSite = ((_Site **)lData)[siteCounter];
              for (long tc = 0; tc < otherDimension->lLength; tc++) {
                tempString->sData[tc] = aSite->sData[otherDimension->lData[tc]];
              }
              tempString->RegExpMatch(regex, matches);
            } else {
              ((_Site **)lData)[siteCounter]->RegExpMatch(regex, matches);
            }
            if (matches.lLength == 0) {
              eligibleMarks[siteCounter] = false;
            }
          }

        DeleteObject(tempString);

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
      FlushRegExp(regex);
    } else {
      input.KillSpaces(input);

      // now process the string
      long count = 0, anchor, k;

      _SimpleList numbers, links;

      numbers.RequestSpace(1024);
      links.RequestSpace(1024);

      // first check if it is has a comb filter
      if ((input.sData[0] == '<') && (input.sData[input.sLength - 1] == '>')) {
        for (count = 1; count < input.sLength - 1; count++) {
          if (input.sData[count] != '0') {
            numbers << count - 1;
          }
        }
        if (numbers.lLength) {
          k = input.sLength - 2; // step size
          anchor = 0;
          if (totalLength == -1) {
            totalLength = theMap.lLength;
          }
          while (anchor < totalLength - k) {
            for (count = 0; count < numbers.lLength; count++) {
              target << anchor + numbers.lData[count];
            }
            anchor += k;
          }
          if ((k = totalLength - 1 - anchor)) {
            for (count = 0; count < numbers.lLength; count++) {
              if (numbers.lData[count] > k) {
                break;
              }
              target << anchor + numbers.lData[count];
            }
          }
          return;
        }
      }

      while (count < input.sLength) {
        anchor = count;
        for (;(count < input.sLength) && (isdigit(input[count])); count++)
          ;
        long aNumber = (input.Cut(anchor, count - 1)).toNum();
        if (aNumber < 0) {
          _String warnMsg(
              "A negative number was found in partition specification: ");
          ReportWarning(warnMsg & input.Cut(0, anchor - 1) & '?' &
                        input.Cut(anchor, -1));
          target.Clear();
          return;
        }
        numbers << aNumber;

        if ((input[count] == '<') || (input[count] == '>')) {
          _String warnMsg(
              "A comb partition cannot be combined with other types. The "
              "entire partition is reset to first..last");
          ReportWarning(warnMsg & input.Cut(0, anchor - 1) & '?' &
                        input.Cut(anchor, -1));
          target.Clear();
          return;
        }
        if (input[count] == '&') {
          links << numbers.lLength;
        }
        if ((input[count] == ',') ||
            (count == input.sLength)) { // wrap it up dude
          if (numbers.lLength == 1) {
            target << numbers(0);
          } else {
            if (links.lLength == 0) {
              if (numbers[0] > numbers[1]) { // backward order
                for (k = numbers[0]; k >= numbers[1]; k--) {
                  target << k;
                }
              } else {
                for (k = numbers[0]; k <= numbers[1]; k++) {
                  target << k;
                }
              }
            } else {
              // linked locations
              if (links.lLength != (numbers.lLength - 2) / 2) {
                _String errMsg("A part of the partition specification has not "
                               "been understood and has been skipped.");
                ReportWarning(errMsg);
                target.Clear();
                return;
              } else {
                _SimpleList signs;
                signs << (numbers(0) < numbers(1) ? 1 : -1);
                for (k = 0; k < links.lLength; k++) {
                  signs << (numbers(links(k)) < numbers(links(k) + 1) ? 1 : -1);
                }

                long l, m;
                for (k = numbers(0), l = 0; signs(0) *k <= signs(0) *numbers(1);
                     k += signs(0), l++) {
                  target << numbers(0) + l * signs(0);
                  for (m = 0; m < links.lLength; m++) {
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

//______________________________________________________________________________
_Matrix *_DataSet::HarvestFrequencies(char unit, char atom, bool posSpec,
                                      _SimpleList &hSegmentation,
                                      _SimpleList &vSegmentation,
                                      bool countGaps) {
                                      
                                      
  // revert to default (all data)
  if (hSegmentation.lLength == 0L || vSegmentation.lLength < unit) { 
    if (hSegmentation.lLength == 0) {
      hSegmentation.Populate(NoOfSpecies(), 0, 1);
    }
    if (vSegmentation.lLength < unit) {
      vSegmentation.Clear();
      vSegmentation.Populate(GetNoTypes(), 0, 1);
    }
  }

  if (unit % atom > 0) { // 20120814 SLKP: changed this behavior to throw errors
    WarnError("Atom should divide unit in HarvestFrequencies call");
    return new _Matrix(1, 1);
  }

  unsigned long vD, alph_dim = theTT->Dimension(),
                    hD = compute_power(alph_dim, atom);

  // create the output Matrix

  vD = posSpec ? unit / atom : 1;

  _Matrix *out = (_Matrix *)checkPointer(new _Matrix(hD, vD, false, true));

  long positions = unit / atom, *store = new long[atom * alph_dim];

  // loop over the set of segments
  // make sure the partition is kosher
  for (unsigned long i = 0; i < vSegmentation.lLength; i += unit) { 

    if (i + unit > vSegmentation.lLength) {
      break;
    }

    for (unsigned long jj = i; jj < i + unit; jj += atom) {
      long k = (jj - i) / atom;

      for (unsigned long ll = 0; ll < hSegmentation.lLength; ll++)
          // loop down each column
          {
        int l = hSegmentation.lData[ll];
        unsigned long count = 1L;
        // build atomic probabilities
        for (unsigned long m = 0; m < atom; m++) {
          theTT->TokenCode((*this)(vSegmentation.lData[jj + m], l, atom),
                           store + alph_dim * m, countGaps);
        }

        long index = 0, shifter = 1;
        for (int m = atom - 1; m >= 0; m--) {
          int smcount = 0;
          for (int n = 0; n < alph_dim; n++) {
            if (store[alph_dim * m + n]) {
              index += shifter * n;
              smcount++;
            }
          }
          shifter *= alph_dim;
          count *= smcount;
        }

        if (count > 1) {
          constructFreq(store, out->theData, posSpec ? positions : 1,
                        posSpec ? k : 0, count, atom - 1, 1, 0);
        } else {
          out->theData[posSpec ? index * positions + k : index] += count;
        }
      }
    }
  }

  delete[] store;
  //scale the matrix now

  hD = out->GetHDim();
  vD = out->GetVDim();

  // normalize each _column_ to sum to 1.
  for (unsigned long i = 0; i < vD; i++) { 
    _Parameter temp = 0.0;

    for (long r = hD - 1; r >= 0; r--) {
      temp += out->theData[r * vD + i];
    }

    for (long r = i; r < vD *hD; r += posSpec ? positions : 1) {
      out->theData[r] /= temp;
    }
  }

  return out;
}

//______________________________________________________________________________
void _DataSet::constructFreq(long *d, _Parameter *m, char positions,
                             long column, long counter, int level, int shifter,
                             int index) {

  const unsigned long alph_len = theTT->Dimension();

  for (unsigned i = 0; i < alph_len; i++) {
    if (d[level * alph_len + i]) {
      if (level) {
        constructFreq(d, m, positions, column, counter, level - 1,
                      shifter * alph_len, index + i * shifter);
      } else {
        m[(index + i * shifter) * positions + column] += 1.0 / counter;
      }
    }
  }
}
