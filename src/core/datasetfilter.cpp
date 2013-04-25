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

#include "datasetfilter.h"
#include "hy_globals.h"
#include "likefunc.h"
#include "site.h"

extern _String dataFileTree, dataFileTreeString, nexusFileTreeMatrix,
    dataFilePartitionMatrix, useTraversalHeuristic, defaultLargeFileCutoff,
    fileTreeString;

_DataSetFilter::_DataSetFilter(void) {
  unitLength = 0;
  theData = NULL;
  accessCache = nil;
}

//______________________________________________________________________________
_DataSetFilter::_DataSetFilter(_DataSet *ds, char, _String &) {
  theData = ds;
  accessCache = nil;
}

//______________________________________________________________________________
_DataSetFilter::~_DataSetFilter(void) {
  if (accessCache) {
    DeleteObject(accessCache);
  }
}

//______________________________________________________________________________

void _DataSetFilter::CopyFilter(_DataSetFilter *copyFrom) {

  memcpy((char *)this, (char *)copyFrom, sizeof(_DataSetFilter));

  theFrequencies.Duplicate(&copyFrom->theFrequencies);
  theNodeMap.Duplicate(&copyFrom->theNodeMap);
  theMap.Duplicate(&copyFrom->theMap);
  theOriginalOrder.Duplicate(&copyFrom->theOriginalOrder);
  conversionCache.Duplicate(&copyFrom->conversionCache);
  duplicateMap.Duplicate(&copyFrom->duplicateMap);

  nInstances = 1;
  dimension = copyFrom->dimension;
  undimension = copyFrom->undimension;
  unitLength = copyFrom->unitLength;
  accessCache = nil;

}

//______________________________________________________________________________

BaseRef _DataSetFilter::makeDynamic(void) {
  _DataSetFilter *r = new _DataSetFilter;
  checkPointer(r);
  r->CopyFilter(this);

  return r;
}

//______________________________________________________________________________
long _DataSetFilter::FreeUpMemory(long requestedBytes) {
  long res = 0;
  for (long i = 0;(i < theMap.lLength) && (res < requestedBytes); i++) {
    res += (theData->GetSite(theMap[i]))->FreeUpMemory(requestedBytes - res);
  }
  return res;
}

//______________________________________________________________________________
void _DataSetFilter::SetDimensions(void) {
  dimension = GetDimension(true);
  undimension = GetDimension(false);
}

//______________________________________________________________________________
unsigned long _DataSetFilter::FindUniqueSequences(_SimpleList &indices,
                                                  _SimpleList &map,
                                                  _SimpleList &counts,
                                                  short mode) {
  indices.Clear();
  map.Clear();
  counts.Clear();

  unsigned long sites = theMap.lLength, seqs = theNodeMap.lLength,
                unit = GetUnitLength();

  if (mode == 0) {
    _SimpleList hashSupport;
    _AVLListXL sequenceHashes(&hashSupport);

    for (unsigned long sequenceIndex = 0; sequenceIndex < seqs;
         sequenceIndex++) {
      _String *thisSequence = GetSequenceCharacters(sequenceIndex);

      long sequenceHash = thisSequence->Adler32(),
           f = sequenceHashes.Find((BaseRef) sequenceHash),
           rawSequenceIdx = theNodeMap.lData[sequenceIndex];

      DeleteObject(thisSequence);

      _SimpleList *sameScore = nil;
      if (f >= 0) {
        sameScore = (_SimpleList *)sequenceHashes.GetXtra(f);
        for (long k = 0; k < sameScore->lLength; k++) {
          bool fit = true;
          f = sameScore->lData[k];

          long fRaw = theNodeMap.lData[indices.lData[f]];

          for (unsigned long site = 0; site < sites && fit; site++) {
            for (unsigned long unitIndex = 0; unitIndex < unit; unitIndex++) {
              _Site *thisSite =
                  theData->GetSite(theMap.lData[unit * site + unitIndex]);
              if (thisSite->sData[fRaw] != thisSite->sData[rawSequenceIdx]) {
                fit = false;
                break;
              }
            }
          }

          if (fit) {
            map << f;
            counts.lData[f]++;

          } else {
            f = -1;
          }
        }
      }
      if (f == -1) { // fit failed or unique site
        if (!sameScore) {
          sameScore = (_SimpleList *)checkPointer(new _SimpleList);
          sequenceHashes.Insert((BaseRef) sequenceHash, (long) sameScore, false);
        }

        (*sameScore) << indices.lLength;
        map << indices.lLength;
        indices << sequenceIndex;
        counts << 1;
      }
    }

  } else {
    long vd = GetDimension(true);

    _Parameter *translatedVector =
                   (_Parameter *)checkPointer(new _Parameter[vd]),
               *translatedVector2 =
                   (_Parameter *)checkPointer(new _Parameter[vd]);

    _String state1(unit, false), state2(unit, false);

    for (long sequenceIndex = 0; sequenceIndex < seqs; sequenceIndex++) {
      bool checkState = false;
      for (long idx = 0; idx < indices.countitems(); idx++) {
        for (long m = 0; m < sites; m++) {
          RetrieveState(m, sequenceIndex, state1, false);
          RetrieveState(m, indices.lData[idx], state2, false);

          checkState = true;
          long idx1 = Translate2Frequencies(state1, translatedVector, true),
               idx2 = Translate2Frequencies(state2, translatedVector2, true);

          //printf ("(%ld, %ld) %ld = %ld %ld\n", sequenceIndex,
          //indices.lData[idx], m, idx1, idx2);

          if (idx2 >= 0 && idx1 >= 0) {
            if (idx1 == idx2) {
              continue;
            } else {
              checkState = false;
              break;
            }
          } else {

            // check for equal ambigs
            long k = 0;
            for (; k < vd; k++) {
              if (translatedVector[k] != translatedVector2[k]) {
                break;
              }
            }

            if (k == vd)
              continue;

            if (mode == 1) {

              long count1 = 0, count2 = 0;

              for (long t = 0; t < vd; t++) {
                count1 += translatedVector[t] > 0.0;
                count2 += translatedVector2[t] > 0.0;
              }

              if (count1 < vd && count2 < vd) {
                checkState = false;
                break;
              }

            } else {
              bool first = mode == 2, second = mode == 2;
              if (mode == 2) {
                for (long t = 0;(first || second) && (t < vd); t++) {
                  if (translatedVector[t] > 0.0) {
                    second &= (translatedVector2[t] > 0.0);
                  }
                  if (translatedVector2[t] > 0.0) {
                    first &= (translatedVector[t] > 0.0);
                  }
                }
                if (!(first || second)) {
                  checkState = false;
                  break;
                }
              } else {
                for (long t = 0; t < vd; t++) {
                  if (translatedVector[t] > 0.0) {
                    second |= (translatedVector2[t] > 0.0);
                  }
                  if (translatedVector2[t] > 0.0) {
                    first |= (translatedVector[t] > 0.0);
                  }
                }
                if (!(first && second)) {
                  checkState = false;
                  break;
                }
              }
            }
          }
        }

        if (checkState) {
          map << idx;
          counts.lData[idx]++;
          break;
        }
      }

      if (!checkState) {
        map << indices.lLength;
        indices << sequenceIndex;
        counts << 1;
      }

    }

    delete[] translatedVector;
    delete[] translatedVector2;
  }

  return indices.lLength;
}

//______________________________________________________________________________

void _DataSetFilter::SetFilter(_DataSet *ds, char unit,
                               _SimpleList &horizontalList,
                               _SimpleList &verticalList,
                               bool isFilteredAlready) {

  // we must step thru the underlying dataset and recompute the frequenices
  // we will store the vertical map in theMap
  // and the horizontal map in theNodeMap
  // theFrequencies will be store the new frequencies
  // theOriginalOrder is the receptacle for the original site order in the data
  // filter

  bool copiedSelf = false; // tag if refiltering self

  _DataSetFilter *firstOne = nil;

  if (isFilteredAlready) {
    if ((Ptr) this == (Ptr) ds) {
      firstOne = (_DataSetFilter *)makeDynamic();
      copiedSelf = true;
    } else {
      firstOne = (_DataSetFilter *)ds;
    }
    ds = firstOne->theData;
  }

  theMap.Clear();
  theNodeMap.Clear();
  theOriginalOrder.Clear();
  theFrequencies.Clear();
  theExclusions.Clear();
  conversionCache.Clear();
  duplicateMap.Clear();

  theData = ds;
  unitLength = unit;

  long i, j;

  // security checks
  if (!horizontalList.lLength || (verticalList.lLength < unit)) {
    ReportWarning(_String("Row and/or column partition is empty. All the data "
                          "will be used by default."));
    if (horizontalList.lLength == 0) {
      if (!isFilteredAlready) {
        j = ds->NoOfSpecies();
      } else {
        j = firstOne->theNodeMap.lLength;
      }
      horizontalList.Populate(j, 0, 1);
    }
    if (verticalList.lLength < unit) {
      verticalList.Clear();
      if (!isFilteredAlready) {
        j = ds->GetNoTypes();
      } else {
        j = firstOne->theOriginalOrder.lLength;
      }
      verticalList.Populate(j, 0, 1);
    }
  }

  if (!isFilteredAlready) {
    theNodeMap.Clear();
    theNodeMap.Duplicate(&horizontalList);
  } else {
    for (long k = 0; k < horizontalList.lLength; k++) {
      theNodeMap << firstOne->theNodeMap.lData[horizontalList.lData[k]];
    }
    horizontalList.Clear();
    horizontalList.Duplicate(&verticalList);
    verticalList.Clear();
    verticalList.RequestSpace(firstOne->theOriginalOrder.lLength);

    for (i = 0; i < horizontalList.lLength; i++) {
      j = horizontalList.lData[i];
      if (j >= 0 && j < firstOne->theOriginalOrder.lLength) {
        verticalList << firstOne->theOriginalOrder.lData[j];
      } else {
        _String tooBig(j);
        if (j < 0) {
          ReportWarning(tooBig & " is a negative site index and is ignored");
        } else {
          ReportWarning(tooBig & " exceeds the number of sites in the "
                                 "underlying data filter and is ignored");
        }
      }
    }
  }

  j = ds->NoOfSpecies();

  for (i = 0; i < theNodeMap.lLength; i++) {
    if (theNodeMap.lData[i] >= j) {
      _String invalid(theNodeMap.lData[i]);
      ReportWarning((invalid & " exceeds the number of species in the "
                               "underlying dataset and is ignored"));
      theNodeMap.Delete(i);
      i--;
    }
  }

  j = ds->GetNoTypes();
  for (i = 0; i < verticalList.lLength; i++)
    if (verticalList.lData[i] >= j) {
      _String invalid(verticalList.lData[i]);
      ReportWarning((invalid & " exceeds the number of sites in the underlying "
                               "dataset and is ignored"));
      verticalList.Delete(i);
      i--;
    }

  if (verticalList.lLength % unit != 0) {
    ReportWarning(
        _String("Number of sites in datasetfilter is not divisible by the unit "
                "- will truncate to the nearest integer"));
    while (verticalList.lLength % unit) {
      verticalList.Delete(verticalList.lLength - 1);
    }
  }

  theOriginalOrder.Duplicate(&verticalList);

  // done with security checks

  _SimpleList indices; // numeric indices intended to facilitate the reindexing
  _AVLListXL siteIndices(&indices);

  // sweep through the columns left to right
  duplicateMap.RequestSpace(verticalList.lLength / unit + 1);

  _String siteHolder(unit * theNodeMap.lLength, false);

  //bool startD = false;

  for (i = 0; i < verticalList.lLength; i += unit) {
    long colIndex = 0;

    for (j = 0; j < unit; j++)                      // sweep within one block
      // sweep down the columns
      for (long k = 0; k < theNodeMap.lLength; k++) {
        //colIndex+=
        //(((_String**)ds->lData)[ds->theMap.lData[verticalList.lData[i+j]]])->sData[theNodeMap.lData[k]];
        siteHolder[colIndex++] = ((
            (_String **)ds->lData)[ds->theMap.lData[verticalList.lData[i + j]]])->sData[theNodeMap.lData[k]];
      }

    colIndex = siteHolder.Adler32();

    long f = siteIndices.Find((BaseRef) colIndex);
    _SimpleList *sameScore = nil;

    if (f >= 0) {

      sameScore = (_SimpleList *)siteIndices.GetXtra(f);
      for (long k = 0; k < sameScore->lLength; k++) {
        bool fit = true;
        f = sameScore->lData[k];
        for (long j = 0; fit && (j < unit); j++) { // sweep within one block
          _Site *site1 = ds->GetSite(verticalList.lData[i + j]),
                *site2 = ds->GetSite(theMap.lData[unit * f + j]);

          // sweep down the columns
          for (long k = 0; k < theNodeMap.lLength; k++) 
            if (site1->sData[theNodeMap.lData[k]] !=
                site2->sData[theNodeMap.lData[k]]) {
              fit = false;
              break;
            }
        }

        if (fit) {
          theFrequencies[f]++;
          duplicateMap << f;
          f = 0;
          break;
        } else {
          f = -1;
        }
      }
    }

    if (f == -1) { // fit failed or unique site
      if (!sameScore) {
        sameScore = (_SimpleList *)checkPointer(new _SimpleList);
        siteIndices.Insert((BaseRef) colIndex, (long) sameScore, false);
      }

      (*sameScore) << theFrequencies.lLength;
      duplicateMap << theFrequencies.lLength;
      theFrequencies << 1;
      for (j = 0; j < unit; j++) {
        theMap << verticalList.lData[i + j];
      }
    }
  }

  siteIndices.Clear();

  duplicateMap.TrimMemory();
  theOriginalOrder.TrimMemory();

  if (copiedSelf) {
    DeleteObject(firstOne);
  }

  SetDimensions();
  FilterDeletions();

}

//______________________________________________________________________________
long _DataSetFilter::FindSpeciesName(_List &s, _SimpleList &r) {

  // MOD 12/16/03
  r.Clear();

  _List newNames;
  _AVLListX matched(&newNames);

  for (long k = 0; k < theNodeMap.lLength; k++) {
    long i = theNodeMap.lData[k];
    _String *uC = new _String(*(_String *)theData->theNames(i));
    uC->UpCase();
    matched.Insert(uC, i);
  }

  for (long m = 0; m < s.lLength; m++) {
    _String ts(*((_String *)s(m)));
    ts.UpCase();
    long f = matched.Find(&ts);
    if (f >= 0) {
      r << matched.GetXtra(f);
    } else {
      break;
    }
  }

  return r.lLength;
}

//______________________________________________________________________________
extern _String skipOmissions;

void _DataSetFilter::FilterDeletions(_SimpleList *theExc) {
  _Parameter skipo;
  checkParameter(skipOmissions, skipo, 0.0);

  if (skipo > .5 || theExc) { // delete omissions
                              //build up the list of "bad" sites
    _SimpleList sitesWithDeletions;
    if (!theExc) {
      for (long i = 0; i < theFrequencies.lLength; i++)
        if (HasDeletions(i)) {
          sitesWithDeletions << i;
        }
    } else {
      _Parameter *store_vec =
          (_Parameter *)checkPointer(new _Parameter[GetDimension(false)]);

      for (long i = 0; i < theFrequencies.lLength; i++) {
        long pos = HasExclusions(i, theExc, store_vec);
        if (pos != -1) {
          sitesWithDeletions << i;
          _String warnMsg((*this)(i, pos));
          warnMsg = warnMsg & " was encountered in sequence " &
                    *GetSequenceName(pos) & " at site pattern " & i &
                    ". All corresponding alignment columns will be removed "
                    "from subsequent analyses.";
          ReportWarning(warnMsg);
        }
      }

      delete[] store_vec;
    }

    if (sitesWithDeletions.lLength == theFrequencies.lLength) {
      _String errMsg("All the sites in the datafilter have deletions and "
                     "removing them creates an empty filter");
      ReportWarning(errMsg);
    }

    _SimpleList allDeleted, dupDeletes;

    for (long k = 0; k < duplicateMap.lLength; k++)
      if (sitesWithDeletions.BinaryFind(duplicateMap.lData[k]) >= 0) {
        dupDeletes << k;
        for (long j = 0; j < unitLength; j++) {
          allDeleted << k *unitLength + j;
        }
      }

    duplicateMap.DeleteList(dupDeletes);
    dupDeletes.Clear();
    theOriginalOrder.DeleteList(allDeleted);
    theFrequencies.DeleteList(sitesWithDeletions);

    for (long i = 0; i < sitesWithDeletions.lLength; i++) {
      long sitePos = sitesWithDeletions.lData[i];

      for (long j = 0; j < unitLength; j++) {
        theMap.lData[sitePos * unitLength + j] = -1;
        dupDeletes << sitePos *unitLength + j;
      }
    }

    if (allDeleted.lLength) {
      /*allDeleted.Sort();*/

      _String warnMsg("The following sites are being omitted:"),
          *s = (_String *)allDeleted.toStr();

      if (!theExc) {
        warnMsg = warnMsg & "(b/c of deletions/omissions)";
      }

      warnMsg = warnMsg & *s;
      DeleteObject(s);
      ReportWarning(warnMsg);

      _SimpleList shiftIdxBy(sitesWithDeletions.lLength +
                             theFrequencies.lLength);

      long shiftBy = sitesWithDeletions.lLength,
           marker = sitesWithDeletions.lData[sitesWithDeletions.lLength - 1],
           markerI = sitesWithDeletions.lLength - 2;

      shiftIdxBy.lLength = sitesWithDeletions.lLength + theFrequencies.lLength;

      for (long i = shiftIdxBy.lLength - 1; i >= 0; i--) {
        if (i == marker) {
          shiftBy--;
          if (markerI >= 0) {
            marker = sitesWithDeletions.lData[markerI];
            markerI--;
          } else {
            marker = -1;
          }
        }
        shiftIdxBy.lData[i] = shiftBy;
      }
      {
        for (long i = 0; i < duplicateMap.lLength; i++) {
          duplicateMap.lData[i] -= shiftIdxBy.lData[duplicateMap.lData[i]];
        }
      }
    }

    // one final pass on theMap to clear it out
    /*for (long i=theMap.lLength-1;i>=0;i--)
        if (theMap(i)<0)
            theMap.Delete(i);*/
    _SimpleList saveMap(theMap);
    theMap.DeleteList(dupDeletes);
    {
      for (long k = 0; k < theMap.lLength; k++)
        if (theMap.lData[k] < 0) {
          saveMap.DeleteList(dupDeletes);
          WarnError("Internal Error in _DataSetFilter::FilterDeletions");
        }
    }
  }

}

//______________________________________________________________________________
_DataSetFilter *_DataSetFilter::PairFilter(long index1, long index2,
                                           _DataSetFilter *result) {
  _SimpleList species;
  species << theNodeMap(index1);
  species << theNodeMap(index2);
  result->SetFilter(theData, unitLength, species, theMap);
  if (theExclusions.lLength) {
    _String *s = (_String *)theExclusions.toStr();
    *s = s->Cut(1, s->Length() - 2);
    result->SetExclusions(s);
    DeleteObject(s);
  }
  return result;
}

//______________________________________________________________________________

void _DataSetFilter::MatchStartNEnd(_SimpleList &order, _SimpleList &positions,
                                    _SimpleList *parent) {
  if (order.lLength == 0) {
    return;
  }

  long p0 = order.lData[0];

  _Parameter uth;
  checkParameter(useTraversalHeuristic, uth, 1.0);

  if (uth > .5) {
    if (parent)
      for (long i = 1; i < order.lLength; i++) {
        unsigned long j = 0, n = theNodeMap.lLength - 1, p0 = parent->lData[i],
                      p1 = order.lData[i];

        while (CompareTwoSites(p0, p1, j)) {
          j++;
        }
        while (CompareTwoSites(p0, p1, n)) {
          n--;
        }
        n = (n << 16) + j;
        positions << n;
      }
    else
      for (long i = 1; i < order.lLength; i++) {
        unsigned long j = 0, n = theNodeMap.lLength - 1, p1 = order.lData[i];

        while (CompareTwoSites(p0, p1, j)) {
          j++;
        }
        while (CompareTwoSites(p0, p1, n)) {
          n--;
        }
        n = (n << 16) + j;
        positions << n;
        p0 = p1;
      }
  } else
    for (long i = 1; i < order.lLength; i++) {
      unsigned long j = 0, n = theNodeMap.lLength - 1;
      n = (n << 16) + j;
      positions << n;
    }

}

//______________________________________________________________________________
void _DataSetFilter::SetExclusions(_String *theList, bool filter) {

  theExclusions.Clear();
  theList->StripQuotes();

  if (theList->sLength == 0) {
    return;
  }

  _List *tokens = theList->Tokenize(',');
  _SimpleList holder;
  _AVLList exclusions(&holder);

  for (long k = 0; k < tokens->lLength; k++) {
    long posMarker = MapStringToCharIndex(*(_String *)((*tokens)(k)));

    if (posMarker < 0) {
      ReportWarning(
          _String("Exclusion request for '") & *(_String *)((*tokens)(k)) &
          "' does not represent a unique state and will therefore be ignored.");
    } else {
      if (exclusions.Insert((BaseRef) posMarker) < 0) {
        ReportWarning(_String("Exclusion symbol for '") &
                      *(_String *)((*tokens)(k)) &
                      "' is included more than once.");
      }
    }
  }

  DeleteObject(tokens);
  exclusions.ReorderList();

  if (filter) {
    FilterDeletions(&holder);
  }

  theExclusions << holder;
}

//______________________________________________________________________________

_String *_DataSetFilter::GetExclusions(void) {
  _String *res = new _String(16L, true);
  checkPointer(res);

  if (theExclusions.lLength) {
    for (long k = 0; k < theExclusions.lLength - 1; k++) {
      (*res) << ConvertCodeToLetters(theExclusions.lData[k], unitLength);
      (*res) << ',';
    }

    (*res) << ConvertCodeToLetters(
                  theExclusions.lData[theExclusions.lLength - 1], unitLength);
  }

  res->Finalize();
  return res;
}

//______________________________________________________________________________

long _DataSetFilter::GetDimension(bool correct) {
  long result = theData->theTT->Dimension();
  for (long i = 1; i < unitLength; i++) {
    result *= theData->theTT->Dimension();
  }
  if (correct) {
    result -= theExclusions.lLength;
  }
  return result;
}

//______________________________________________________________________________

_String &_DataSetFilter::operator()(unsigned long site, unsigned long pos) {
  if (!accessCache || accessCache->sLength != unitLength) {
    if (accessCache) {
      DeleteObject(accessCache);
    }
    checkPointer(accessCache = new _String((unsigned long) unitLength, false));
  }

  long vIndex = theNodeMap.lData[pos];
  if (unitLength == 1) {
    accessCache->sData[0] = ((
        (_String **)theData->lData)[theData->theMap.lData[theMap.lData[site]]])
        ->sData[vIndex];
  } else {
    site *= unitLength;
    for (int k = 0; k < unitLength; k++) {
      accessCache->sData[k] =
          (((_String **)theData
                ->lData)[theData->theMap.lData[theMap.lData[site++]]])
              ->sData[vIndex];
    }
  }
  return *accessCache;
}

//______________________________________________________________________________

void _DataSetFilter::RetrieveState(unsigned long site, unsigned long pos,
                                   _String &reply, bool map) {
  long vIndex = theNodeMap.lData[pos];
  if (map) {
    if (unitLength == 1) {
      reply.sData[0] = (((_String **)theData->lData)[
          theData->theMap.lData[theMap.lData[duplicateMap.lData[site]]]])
          ->sData[vIndex];
    } else {
      site = unitLength * duplicateMap.lData[site];
      for (int k = 0; k < unitLength; k++) {
        reply.sData[k] =
            (((_String **)theData
                  ->lData)[theData->theMap.lData[theMap.lData[site++]]])
                ->sData[vIndex];
      }
    }
  } else {
    if (unitLength == 1) {
      reply.sData[0] =
          (((_String **)theData
                ->lData)[theData->theMap.lData[theMap.lData[site]]])
              ->sData[vIndex];
    } else
      for (int k = 0; k < unitLength; k++) {
        reply.sData[k] =
            (((_String **)theData
                  ->lData)[theData->theMap.lData[theMap.lData[site++]]])
                ->sData[vIndex];
      }
  }
}

//______________________________________________________________________________

void _DataSetFilter::GrabSite(unsigned long site, unsigned long pos,
                              _String &storage) {

  long vIndex = theNodeMap.lData[pos];
  if (unitLength == 1) {
    storage.sData[0] = ((
        (_String **)theData->lData)[theData->theMap.lData[theMap.lData[site]]])
        ->sData[vIndex];
  } else {
    site *= unitLength;
    for (int k = 0; k < unitLength; k++) {
      storage.sData[k] =
          (((_String **)theData
                ->lData)[theData->theMap.lData[theMap.lData[site++]]])
              ->sData[vIndex];
    }
  }
}

//______________________________________________________________________________

void _DataSetFilter::GrabSite(unsigned long site, unsigned long pos, char *s) {
  long vIndex = theNodeMap.lData[pos];
  if (unitLength == 1) {
    s[0] = (
        ((_String **)theData->lData)[theData->theMap.lData[theMap.lData[site]]])
        ->sData[vIndex];
  } else {
    site *= unitLength;
    for (int k = 0; k < unitLength; k++) {
      s[k] = (((_String **)theData
                   ->lData)[theData->theMap.lData[theMap.lData[site++]]])
          ->sData[vIndex];
    }
  }
}

//______________________________________________________________________________

_SimpleList *_DataSetFilter::CountAndResolve(long pattern, _Parameter *storage,
                                             bool randomly) {
  // last cell in the list contains the count of distinct characters in the
  // column
  _SimpleList *resList = new _SimpleList(theNodeMap.lLength + 1, 0, 0),
              counts(dimension, 0, 0);

  checkPointer(resList);

  _List ambStates;
  _String aState(unitLength, false);

  _Parameter *freqStorage = storage;

  if (!freqStorage) {
    freqStorage = new _Parameter[undimension];
  }

  long normalizingSum = 0, charCount = 0;

  for (long k = 0; k < theNodeMap.lLength; k++) {
    GrabSite(pattern, k, aState);
    long characterRes = Translate2Frequencies(aState, freqStorage, true);
    if (characterRes >= 0) {
      resList->lData[k] = characterRes;

      if (characterRes >= dimension) {
        WarnError(
            _String("Internal error in _DataSetFilter::CountAndResolve\n"));
      }

      if ((counts.lData[characterRes]++) == 0) {
        normalizingSum++;
      }
      charCount++;
    } else {
      _SimpleList *possibleResolutions = new _SimpleList;
      if (!possibleResolutions) {
        checkPointer(possibleResolutions);
      }

      (*possibleResolutions) << k;

      for (long m = 0; m < dimension; m++)
        if (freqStorage[m] > 0.) {
          (*possibleResolutions) << m;
        }

      ambStates.AppendNewInstance(possibleResolutions);

    }
  }

  if (normalizingSum > 0) {
    if (ambStates.lLength) {
      _SimpleList ambResolutions(dimension, 0, 0);
      for (long t = 0; t < ambStates.lLength; t++) {
        _SimpleList *stateResolutions = (_SimpleList *)ambStates(t);

        if (!randomly) {
          long totalSum = 0, idx = 0;

          for (long l = 1; l < stateResolutions->lLength; l++) {
            long tmp = counts.lData[stateResolutions->lData[l]];
            if (tmp > totalSum) {
              idx = l;
              totalSum = tmp;
            }
          }
          if (idx > 0)
              // if no resolutions, resolve randomly
              {
            idx = stateResolutions->lData[idx];
            resList->lData[stateResolutions->lData[0]] = idx;
            ambResolutions.lData[idx]++;
            continue;
          }

        }

        long totalSum = 0;
        for (long l = 1; l < stateResolutions->lLength; l++) {
          totalSum += counts.lData[stateResolutions->lData[l]];
        }

        if (totalSum > 0) {
          long randomN = genrand_real2() * totalSum -
                         counts.lData[stateResolutions->lData[1]],
               ind = 1;

          while (randomN > 0) {
            randomN -= counts.lData[stateResolutions->lData[++ind]];
          }

          totalSum = stateResolutions->lData[ind];
        } else {
          long randomN = genrand_real2() * charCount - counts.lData[0], ind = 0;

          while (randomN > 0) {
            randomN -= counts.lData[++ind];
          }
        }
        resList->lData[stateResolutions->lData[0]] = totalSum;
        ambResolutions.lData[totalSum]++;
      }

      for (long l = 0; l < dimension; l++)
        if (ambResolutions.lData[l] && !counts.lData[l]) {
          normalizingSum++;
        }
    }
  }

  resList->lData[theNodeMap.lLength] = normalizingSum;

  if (freqStorage != storage) {
    delete freqStorage;
  }

  return resList;
}

//______________________________________________________________________________

_Matrix *_DataSetFilter::PairwiseCompare(_SimpleList *s1, _SimpleList *s2,
                                         _List *labels) {

  // s1 and s2 are the lists produced by CountAndResolve
  // if labels is not nil, then it will receive row and column labels in the
  // contigency table
  // the result matrix has rows labeled by states in s1, and columns - by
  // states in s2

  long *sort1 = new long[dimension], *sort2 = new long[dimension],
       c = s2->lData[s2->lLength - 1];

  _Matrix *res = new _Matrix(s1->lData[s1->lLength - 1], c, false, true);

  if (sort1 && sort2 && res) {
    for (long k = 0; k < dimension; k++) {
      sort1[k] = -1;
      sort2[k] = -1;
    }

    long idx1 = 0, idx2 = 0;

    _SimpleList *lbl1 = nil, *lbl2 = nil;

    if (labels) {
      lbl1 = new _SimpleList;
      lbl2 = new _SimpleList;

      checkPointer(lbl1);
      checkPointer(lbl2);

      (*labels) << lbl1;
      (*labels) << lbl2;

      DeleteObject(lbl1);
      DeleteObject(lbl2);
    }

    for (long k2 = 0; k2 < s1->lLength - 1; k2++) {
      long c1 = s1->lData[k2], c2 = s2->lData[k2];

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

      res->theData[c1 * c + c2] += 1.;
    }

    delete sort1;
    delete sort2;
  } else {
    checkPointer(nil);
  }

  return res;
}

//______________________________________________________________________________
_List *_DataSetFilter::ComputePatternToSiteMap(void) {
  _List *result = new _List();
  for (long k = 0; k < theFrequencies.lLength; k++) {
    result->AppendNewInstance(new _SimpleList);
  }
  for (long s = 0; s < duplicateMap.lLength; s++) {
    *((_SimpleList **)result->lData)[duplicateMap.lData[s]] << s;
  }
  return result;
}

//______________________________________________________________________________
char _DataSetFilter::GetChar(unsigned long site, unsigned long pos) {
  //long vIndex = theNodeMap.lLength?theNodeMap.lData[pos]:pos;
  return (*theData)(theMap.lData[site], theNodeMap.lData[pos], 1);
}

//______________________________________________________________________________
bool _DataSetFilter::CompareTwoSites(unsigned long site1, unsigned long site2,
                                     unsigned long pos1) {
  pos1 = theNodeMap.lData[pos1];

  if (unitLength == 3) { // codon
    site1 *= 3;
    site2 *= 3;
    if (((((_String **)theData
               ->lData)[theData->theMap.lData[theMap.lData[site1]]])
             ->sData[pos1] == (((_String **)theData->lData)[
                                  theData->theMap.lData[theMap.lData[site2]]])
             ->sData[pos1]) &&
        ((((_String **)theData
               ->lData)[theData->theMap.lData[theMap.lData[site1 + 1]]])
             ->sData[pos1] ==
         (((_String **)theData
               ->lData)[theData->theMap.lData[theMap.lData[site2 + 1]]])
             ->sData[pos1]) &&
        ((((_String **)theData
               ->lData)[theData->theMap.lData[theMap.lData[site1 + 2]]])
             ->sData[pos1] ==
         (((_String **)theData->lData)[
             theData->theMap.lData[theMap.lData[site2 + 2]]])->sData[pos1])) {
      return true;
    }
  } else {
    site1 *= unitLength;
    site2 *= unitLength;
    long k;

    /*if
    ((((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1]]])->sLength<=pos1)
    {
        printf
    ("(%d)%s\n(%d)%s\n",site1,(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1]]])->sData,
                site2,(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site2]]])->sData);
        FlagError ("Internal DataSetFilter bug\n");
    }*/

    for (k = 0; k < unitLength; k++) {
      if ((((_String **)theData
                ->lData)[theData->theMap.lData[theMap.lData[site1]]])
              ->sData[pos1] !=
          (((_String **)theData->lData)[
              theData->theMap.lData[theMap.lData[site2]]])->sData[pos1]) {
        break;
      }
      site1++;
      site2++;
    }
    if (k == unitLength) {
      return true;
    }
  }
  return false;
}

//______________________________________________________________________________
bool _DataSetFilter::CompareTwoSitesChar(unsigned long site1,
                                         unsigned long site2,
                                         unsigned long pos1) {

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
  return ((*theData)(theMap.lData[site1], pos1, 1) ==
          (*theData)(theMap.lData[site2], pos1, 1));
}

//______________________________________________________________________________
long _DataSetFilter::SiteFrequency(unsigned long site) {
  return theFrequencies.lData[site];
}

//______________________________________________________________________________
bool _DataSetFilter::HasDeletions(unsigned long site, _AVLList *storage) {
  long loopDim = GetDimension();
  _Parameter *store = new _Parameter[loopDim];

  if (!store) {
    warnError(-108);
  }

  long j,
      upTo = theNodeMap.lLength ? theNodeMap.lLength : theData->NoOfSpecies();

  bool outcome = false;

  for (unsigned int k = 0; k < upTo; k++) {
    Translate2Frequencies((*this)(site, k), store, false);

    bool oneF = false, zeroF = false;

    for (j = 0; j < loopDim; j++) {
      if (store[j] == 0.0) {
        zeroF = true;
      } else if (store[j] == 1.0) {
        oneF = true;
      }
    }
    if (!(oneF && zeroF)) {
      if (storage) {
        outcome = true;
        storage->Insert((BaseRef) theNodeMap.lData[k]);
      } else {
        delete store;
        return true;
      }
    }
  }

  delete store;
  return outcome;
}

//______________________________________________________________________________
bool _DataSetFilter::IsConstant(unsigned long site, bool relaxedDeletions) {

  _Parameter *store = nil, *store2 = nil;

  store = new _Parameter[GetDimension()];
  store2 = new _Parameter[GetDimension()];

  if (!(store && store2)) {
    warnError(-108);
  }
  long j,
      upTo = theNodeMap.lLength ? theNodeMap.lLength : theData->NoOfSpecies(),
      loopDim = GetDimension();

  Translate2Frequencies((*this)(site, 0), store, false);

  if (relaxedDeletions) {
    for (unsigned int k = 1; k < upTo; k++) {
      Translate2Frequencies((*this)(site, k), store2, false);
      for (j = 0; j < loopDim; j++) {
        if (store2[j] == 0.0) {
          store[j] = 0.0;
        }
      }
    }
    for (j = 0; j < loopDim; j++)
      if (store[j] != 0.0) {
        delete store;
        delete store2;
        return true;
      }
    if (j == loopDim) {
      delete store;
      delete store2;
      return false;
    }
  } else {
    for (unsigned int k = 1; k < upTo; k++) {
      Translate2Frequencies((*this)(site, k), store2, false);
      for (j = 0; j < loopDim; j++)
        if (store[j] != store2[j]) {
          delete store;
          delete store2;
          return false;
        }
    }
  }

  return true;
}

//______________________________________________________________________________
_Matrix *_DataSetFilter::GetFilterCharacters(bool flip) {
  long unitLength = GetUnitLength(),
       seqLength = flip ? theFrequencies.lLength
                        : (GetFullLengthSpecies() / unitLength),
       f = NumberSpecies();

  _List result;

  _String aState((long) GetUnitLength(), false);

  if (flip) {
    for (long k = 0; k < seqLength; k++) {
      _String *aSite = new _String(128L, true);
      for (long k2 = 0; k2 < f; k2++) {
        RetrieveState(k, k2, aState, false);
        (*aSite) << aState;
      }
      aSite->Finalize();
      result << aSite;
      DeleteObject(aSite);
    }
  } else
    for (long k = 0; k < f; k++) {
      _String *fstr = GetSequenceCharacters(k);
      result << fstr;
      DeleteObject(fstr);
    }

  return new _Matrix(result);

}

//______________________________________________________________________________
_String *_DataSetFilter::GetSequenceCharacters(long seqID) {
  long unitSizeL = GetUnitLength();
  _String *aSequence = new _String(theOriginalOrder.lLength, true);

  if (seqID >= 0 && seqID < theNodeMap.lLength) {
    _String aState(unitSizeL, false);
    long upTo = theOriginalOrder.lLength / unitSizeL;
    for (long k2 = 0; k2 < upTo; k2++) {
      RetrieveState(k2, seqID, aState);
      (*aSequence) << aState;
    }
  }
  aSequence->Finalize();
  return aSequence;
}

//______________________________________________________________________________
long _DataSetFilter::HasExclusions(unsigned long site, _SimpleList *theExc,
                                   _Parameter *store) {
  long filterDim = GetDimension(false);

  if (theNodeMap.lLength)
    for (unsigned long k = 0; k < theNodeMap.lLength; k++) {
      Translate2Frequencies((*this)(site, k), store, false);

      long j = 0, s = 0;

      for (j = 0; j < filterDim; j++)
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

//______________________________________________________________________________
void _DataSetFilter::Freeze(long site) {
  for (int k = 0; k < unitLength; k++) {
    _Site *tC = (_Site *)((*(_List *)theData)(
        theData->theMap(this->theMap(site * unitLength + k))));
    tC->SetRefNo(-1);
    tC->PrepareToUse();
  }
}

//______________________________________________________________________________
void _DataSetFilter::UnFreeze(long site) {
  for (int k = 0; k < unitLength; k++) {
    _Site *tC = (_Site *)((*(_List *)theData)(
        theData->theMap(this->theMap(site * unitLength + k))));
    tC->SetRefNo(0);
    //      tC->Archive();
  }
}

//______________________________________________________________________________
_Matrix *_DataSetFilter::ComputePairwiseDifferences(long i, long j, char amb) {
  if (unitLength > 3) {
    WarnError("ComputePairwiseDifferences is not implemented for data filters "
              "with unit size > 3");
    return new _Matrix(1, 1, false, true);
  }

  long mxDim = GetDimension(true);

  _Matrix *res = new _Matrix(mxDim, mxDim, false, true);

  _Parameter *sm1 = new _Parameter[mxDim], *sm2 = new _Parameter[mxDim];

  checkPointer(res);
  checkPointer(sm1);
  checkPointer(sm2);

  _String state1(unitLength, false), state2(unitLength, false);

  if (!conversionCache.lLength) {
    SetupConversion();
  }

  long *tcodes = conversionCache.lData + 89,
       *ccodes = conversionCache.lData + 1, ccount = conversionCache.lData[0];

  for (long k = 0; k < theFrequencies.lLength; k++) {
    long s1 = -1, s2 = -1;

    int c1, c2;

    c1 = (((_String **)theData
               ->lData)[theData->theMap.lData[theMap.lData[unitLength * k]]])
        ->sData[theNodeMap.lData[i]],
    c2 = (((_String **)theData
               ->lData)[theData->theMap.lData[theMap.lData[unitLength * k]]])
        ->sData[theNodeMap.lData[j]];

    if (unitLength == 1) {
      s1 = conversionCache.lData[(c1 - 40) * (undimension + 1) + undimension],
      s2 = conversionCache.lData[(c2 - 40) * (undimension + 1) + undimension];
    } else {
      int c12 = (((_String **)theData->lData)[
              theData->theMap.lData[theMap.lData[unitLength * k + 1]]])
              ->sData[theNodeMap.lData[i]],
          c22 = (((_String **)theData->lData)[
              theData->theMap.lData[theMap.lData[unitLength * k + 1]]])
              ->sData[theNodeMap.lData[j]];

      state1.sData[0] = c1;
      state1.sData[1] = c12;

      state2.sData[0] = c2;
      state2.sData[1] = c22;

      c1 = ccodes[c1 - 40];
      c12 = ccodes[c12 - 40];

      c2 = ccodes[c2 - 40];
      c22 = ccodes[c22 - 40];

      if (unitLength == 2) {
        if ((c1 >= 0) && (c12 >= 0)) {
          s1 = tcodes[c1 * ccount + c12];
        }

        if ((c2 >= 0) && (c22 >= 0)) {
          s2 = tcodes[c2 * ccount + c22];
        }
      } else {
        int c13 = (((_String **)theData->lData)[
                theData->theMap.lData[theMap.lData[unitLength * k + 2]]])
                ->sData[theNodeMap.lData[i]],
            c23 = (((_String **)theData->lData)[
                theData->theMap.lData[theMap.lData[unitLength * k + 2]]])
                ->sData[theNodeMap.lData[j]];

        //printf ("\n%c %c", c13, c23);

        state1.sData[2] = c13;
        state2.sData[2] = c23;

        c13 = ccodes[c13 - 40];
        c23 = ccodes[c23 - 40];

        //printf (" %d %d %s %s\n", c13, c23, state1.sData, state2.sData);

        if ((c1 >= 0) && (c12 >= 0) && (c13 >= 0)) {
          s1 = tcodes[ccount * (c1 * ccount + c12) + c13];
        }

        if ((c2 >= 0) && (c22 >= 0) && (c23 >= 0)) {
          s2 = tcodes[ccount * (c2 * ccount + c22) + c23];
        }
      }
    }

    if (s1 >= 0 && s2 >= 0)
        // one to one
        {
      res->theData[s1 * mxDim + s2] += theFrequencies.lData[k];
    } else {
      if (amb < 3) {
        _Matrix *freqsAtSite = nil;
        if (amb) {
          _SimpleList //seqList,
              siteList;

          //seqList  << theNodeMap[i];
          //seqList  << theNodeMap[j];

          for (long si = 0; si < unitLength; si++) {
            siteList << theMap.lData[unitLength * k + si];
          }

          freqsAtSite = theData->HarvestFrequencies(unitLength, unitLength, 0,
                                                    theNodeMap, siteList);
          if (theExclusions.lLength) {
            long k = 0, u = GetDimension(false);

            for (long i = 0; i < u; i++) {
              if (i == theExclusions.lData[k] && k < theExclusions.lLength) {
                k++;
                continue;
              }
              freqsAtSite->theData[i - k] = freqsAtSite->theData[i];
            }
          }
          //XferwCorrection (freqsAtSite->theData, freqsAtSite->theData, mxDim);
        }

        if (s1 >= 0)
            // one to many
            {
          if (unitLength > 1) {
            Translate2Frequencies(state2, sm1, false);
          } else {
            Translate2Frequencies(c2, sm1, false);
          }

          if (freqsAtSite) {
            if (amb == 1) {
              _Parameter totalW = 0.0;

              for (long m = 0; m < mxDim; m++)
                if (sm1[m] > 0.0) {
                  totalW += freqsAtSite->theData[m];
                }

              if (totalW > 0.0) {
                s1 = s1 * mxDim;

                for (long m = 0; m < mxDim; m++, s1++)
                  if (sm1[m] > 0.0) {
                    res->theData[s1] += theFrequencies.lData[k] *
                                        freqsAtSite->theData[m] / totalW;
                  }
              }

            } else {
              _Parameter maxW = 0.0;
              long maxIdx = -1;

              for (long m = 0; m < mxDim; m++) {
                if (sm1[m] > 0.0) {
                  _Parameter myWeight = freqsAtSite->theData[m];
                  if (myWeight > maxW) {
                    maxW = myWeight;
                    maxIdx = m;
                  }
                }
              }

              if (maxIdx >= 0) {
                res->theData[s1 * mxDim + maxIdx] += theFrequencies.lData[k];
              }
            }
          } else {
            /* adopt the following convention here:
                - if ambig resolves to one s1 - count as a match
                - otherwise - count all contributions equally
            */

            if (sm1[s1] > 0.0) {
              res->theData[s1 * mxDim + s1] += theFrequencies.lData[k];
            } else {
              long ambCount = 0;
              {
                for (long m = 0; m < mxDim; m++)
                  if (sm1[m] > 0.0) {
                    ambCount++;
                  }
              }
              s1 *= mxDim;

              _Parameter addFac =
                  theFrequencies.lData[k] / (_Parameter) ambCount;

              for (long m = 0; m < mxDim; m++, s1++)
                if (sm1[m] > 0.0) {
                  res->theData[s1] += addFac;
                }
            }
          }
        } else {
          if (s2 >= 0)
              // many to one
              {
            if (unitLength > 1) {
              Translate2Frequencies(state1, sm1, false);
            } else {
              Translate2Frequencies(c1, sm1, false);
            }

            if (freqsAtSite) {
              if (amb == 1) {
                _Parameter totalW = 0.0;

                for (long m = 0; m < mxDim; m++)
                  if (sm1[m] > 0.0) {
                    totalW += freqsAtSite->theData[m];
                  }

                if (totalW > 0.0) {
                  for (long m = 0; m < mxDim; m++, s2 += mxDim)
                    if (sm1[m] > 0.0) {
                      res->theData[s2] += theFrequencies.lData[k] *
                                          freqsAtSite->theData[m] / totalW;
                    }
                }

              } else {
                _Parameter maxW = 0.0;
                long maxIdx = -1;

                for (long m = 0; m < mxDim; m++) {
                  if (sm1[m] > 0.0) {
                    _Parameter myWeight = freqsAtSite->theData[m];
                    if (myWeight > maxW) {
                      maxW = myWeight;
                      maxIdx = m;
                    }
                  }
                }

                if (maxIdx >= 0) {
                  res->theData[maxIdx * mxDim + s2] += theFrequencies.lData[k];
                }
              }
            } else {
              if (sm1[s2] > 0.0) {
                res->theData[s2 * mxDim + s2] += theFrequencies.lData[k];
              } else {
                long ambCount = 0;
                for (long m = 0; m < mxDim; m++)
                  if (sm1[m] > 0.0) {
                    ambCount++;
                  }

                _Parameter addFac =
                    theFrequencies.lData[k] / (_Parameter) ambCount;
                {
                  for (long m = 0; m < mxDim; m++, s2 += mxDim)
                    if (sm1[m] > 0.0) {
                      res->theData[s2] += addFac;
                    }
                }
              }
            }
          } else
              // many to many
              {
            if (unitLength > 1) {
              Translate2Frequencies(state1, sm1, false);
              Translate2Frequencies(state2, sm2, false);
            } else {
              Translate2Frequencies(c1, sm1, false);
              Translate2Frequencies(c2, sm2, false);
            }

            if (freqsAtSite) {
              if (amb == 1) {
                _Parameter totalW = 0.0;

                for (long m = 0; m < mxDim; m++)
                  if (sm1[m] > 0)
                    for (long m2 = 0; m2 < mxDim; m2++)
                      if (sm2[m2] > 0) {
                        totalW +=
                            freqsAtSite->theData[m] * freqsAtSite->theData[m2];
                      }

                if (totalW > 0.0) {
                  for (long m = 0; m < mxDim; m++)
                    if (sm1[m] > 0)
                      for (long m2 = 0; m2 < mxDim; m2++)
                        if (sm2[m2] > 0) {
                          res->theData[m * mxDim + m2] +=
                              theFrequencies.lData[k] *
                              freqsAtSite->theData[m] *
                              freqsAtSite->theData[m2] / totalW;
                        }
                }

              } else {
                _Parameter maxW = 0.0;
                long maxIdx = -1, maxIdx2 = -1;

                for (long m = 0; m < mxDim; m++)
                  if (sm1[m] > 0)
                    for (long m2 = 0; m2 < mxDim; m2++)
                      if (sm2[m2] > 0) {
                        _Parameter myWeight =
                            freqsAtSite->theData[m] * freqsAtSite->theData[m2];
                        if (myWeight > maxW) {
                          maxW = myWeight;
                          maxIdx = m;
                          maxIdx2 = m2;
                        }
                      }

                if (maxIdx >= 0) {
                  res->theData[maxIdx * mxDim + maxIdx2] +=
                      theFrequencies.lData[k];
                }
              }
            } else {
              long ambCount = 0, ambCount2 = 0, m = 0;

              for (; m < mxDim; m++) {
                if (sm1[m] > 0.0) {
                  if (sm2[m] > 0.0) {
                    break;
                  } else {
                    ambCount++;
                  }
                } else if (sm2[m] > 0.0) {
                  ambCount2++;
                }
              }

              if (m == mxDim) {
                _Parameter addFac = theFrequencies.lData[k] /
                                    (_Parameter)(ambCount * ambCount2);

                for (long m = 0; m < mxDim; m++)
                  if (sm1[m] > 0)
                    for (long m2 = 0; m2 < mxDim; m2++)
                      if (sm2[m2] > 0) {
                        res->theData[m * mxDim + m2] += addFac;
                      }
              }
            }
          }
        }
        DeleteObject(freqsAtSite);
      }
    }
  }

  delete[] sm1;
  delete[] sm2;

  return res;
}

//______________________________________________________________________________
void _DataSetFilter::ComputePairwiseDifferences(_Matrix &target, long i, long j) {

  // matrix of dimension nx4n containing pairwise distances as follows
  // (n=number of species)
  // first lower diag - count the same (AA,CC,GG,TT)
  // first upper diag - count AC,CA
  // 2nd   lower diag - count AG,GA
  // 2nd   upper diag - count AT,TA
  // 3rd   lower diag - count CG,GC
  // 3rd   upper diag - count CT,TC
  // 4th   lower diag - count GT,TG

  if ((target.GetHDim() != 1) || (target.GetVDim() != 7)) {
    CreateMatrix(&target, 1, 7, false, true, false);
  }

  if (theData->theTT->DetectType() !=
      HY_TRANSLATION_TABLE_STANDARD_NUCLEOTIDE) {
    return;
  }
  long k, l, m;

  for (k = 0; k < 7; k++) {
    target[k] = 0;
  }
  k = theNodeMap.lData[i];
  l = theNodeMap.lData[j];
  if (l > k) {
    m = l;
    l = k;
    k = m;
  }

  for (m = theMap.lLength - 1; m > -1; m--) {
    char *thisSite = GetColumn(m);
    char a = thisSite[k], b = thisSite[l], c;

    long fc = theFrequencies.lData[m / unitLength];

    if (a > b) {
      c = a;
      a = b;
      b = c;
    }
    if (a == b) {
      target[0] += fc;
    } else {
      if (a == 'A') {
        switch (b) {
        case 'C': {
          target[1] += fc;
          break;
        }
        case 'G': {
          target[2] += fc;
          break;
        }
        case 'T': {
          target[3] += fc;
          break;
        }
        }
      } else if (a == 'C') {
        switch (b) {
        case 'G': {
          target[4] += fc;
          break;
        }
        case 'T': {
          target[5] += fc;
          break;
        }
        }
      } else if (a == 'G') {
        if (b == 'T') {
          target[6] += fc;
        }
      }

    }
  }
}

//______________________________________________________________________________

_Matrix *_DataSetFilter::HarvestFrequencies(char unit, char atom, bool posSpec,
                                            bool countGaps) {

  return theData->HarvestFrequencies(unit, atom, posSpec, theNodeMap,
                                     theOriginalOrder, countGaps);
}

//______________________________________________________________________________
void _DataSetFilter::XferwCorrection(_Matrix &source, _Parameter *target,
                                     long _length) {
  long k = 0;
  _Parameter *mxdata = source.fastIndex();
  if (theExclusions.lLength == 0) {
    for (long i = 0; i < _length; i++) {
      target[i] = (mxdata[i] != 0.0);
    }
  } else {
    for (long i = 0; i < _length; i++) {
      if (k < theExclusions.lLength && i == theExclusions.lData[k]) {
        k++;
        continue;
      }
      target[i - k] = (mxdata[i] != 0.);
    }
  }
}

//______________________________________________________________________________
void _DataSetFilter::XferwCorrection(_Parameter *source, _Parameter *target,
                                     long _length) {
  long k = 0;
  if (theExclusions.lLength == 0) {
    for (long i = 0; i < _length; i++) {
      target[i] = (source[i] != 0.0);
    }
  } else {
    for (long i = 0; i < _length; i++) {
      if (i == theExclusions.lData[k] && k < theExclusions.lLength) {
        k++;
        continue;
      }
      target[i - k] = (source[i] != 0);
    }
  }
}

//______________________________________________________________________________
void _DataSetFilter::XferwCorrection(long *source, _Parameter *target,
                                     long _length) {
  long k = 0;
  if (theExclusions.lLength == 0) {
    for (long i = 0; i < _length; i++) {
      target[i] = source[i];
    }
  } else {
    for (long i = 0; i < _length; i++) {
      if (i == theExclusions[k]) {
        k++;
        continue;
      }
      target[i - k] = source[i];
    }
  }
}

/*
//______________________________________________________________________________
long    _DataSetFilter::GetVectorCode(long site,long seq)
{
    if (!symbolVector) return -1;
    long* fi = symbolVector->quickArrayAccess();
    return fi[*fi*site+seq+1];
}
//______________________________________________________________________________

void    _DataSetFilter::ProduceSymbolVector(bool smear)
{
    // compute the size of the vector cells
    _Parameter
cellSize=log((_Parameter)theData->theTT->LengthOfAlphabet())*_Parameter(unitLength)/log(128.0);
    if (cellSize>2.0)
    {
        _String errMsg ("DataSetFilter has more than 32767 states, which is
currently unsupported");
        FlagError(errMsg);
    }
    long intCellSize = cellSize>1.0?2:1;
    // now produce the conversion vector
    long sites = theMap.lLength, species=
theNodeMap.lLength?theNodeMap.lLength:theData->NoOfSpecies();
    symbolVector = new _SimpleList ();
    checkPointer(symbolVector);
//  (*symbolVector)<<intCellSize;
    (*symbolVector)<<species;
    // we will now speciate into byte and word size cases
    // the data will be stored column by column
    // if there is a unique code translation, we then store that code in the
symbol vector for faster
    // processing during tree pruning business.
    // use a standard convert to frequencies function to check whether a
character has a unique convertion
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
//
converterb.bytes[byteposition]=(char)Translate2Frequencies((*this)(i,j),nil,smear,false);
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

//______________________________________________________________________________
long _DataSetFilter::CorrectCode(long code) {
  if (theExclusions.lLength != 0) {
    for (long k = 0; k < theExclusions.lLength; k++)
      if (code >= theExclusions.lData[k]) {
        code++;
      }
  }
  return code;
}

//______________________________________________________________________________
long _DataSetFilter::Translate2Frequencies(_String &str, _Parameter *parvect,
                                           bool smear) {

  long count = 0, nonzeropos = 0, store[HYPHY_SITE_DEFAULT_BUFFER_SIZE];

  if (unitLength == 1) {
    theData->theTT->TokenCode(str.sData[0], store, smear);
    if (theExclusions.lLength == 0) {
      for (long i = 0; i < undimension; i++)
        if ((parvect[i] = store[i])) {
          nonzeropos = i;
          count++;
        }
    } else {
      long k = 0;
      for (long i = 0; i < undimension; i++) {
        if (i == theExclusions.lData[k] && k < theExclusions.lLength) {
          k++;
        } else if (store[i]) {
          nonzeropos = i;
          count++;
        }
        parvect[i - k] = store[i];
      }
    }
    if (count == 0) {
      if (smear)
        for (long i = 0; i < undimension; i++) {
          parvect[i] = 1.;
        }

      return -1;
    }

    if (count > 1) {
      return -1;
    }

    return nonzeropos;
  } else {
    //pull the frequencies out of the Translation table
    _Matrix out(undimension, 1, false, true);

    long m, n, index = 0, shifter = 1, *lp, *storeP,
               alphabet_dimension = theData->theTT->Dimension();

    _Parameter *fl;

    if (alphabet_dimension * unitLength >= HYPHY_SITE_DEFAULT_BUFFER_SIZE) {
      storeP = new long[alphabet_dimension * unitLength];
    } else {
      storeP = store;
    }

    count = 1;
    for (m = 0; m < unitLength; m++) {
      theData->theTT->TokenCode(str.sData[m], storeP + alphabet_dimension * m);
    }

    for (m = unitLength - 1; m >= 0; m--) {
      int smcount = 0;
      lp = storeP + alphabet_dimension * m;
      for (n = 0; n < alphabet_dimension; n++, lp++) {
        if (*lp) {
          index += shifter * n;
          smcount++;
        }
      }

      if ((smcount == 0) && smear) { // deletion -- replace with 1's
        lp = storeP + alphabet_dimension * m;
        for (n = 0; n < alphabet_dimension; n++, lp++) {
          *lp = 1;
        }
        smcount = alphabet_dimension;
      }

      shifter *= alphabet_dimension;
      count *= smcount;
    }

    if (count > 1) {
      theData->constructFreq(storeP, out.theData, 1, 0, count, unitLength - 1,
                             1, 0);
    } else if (count == 1) {
      out.theData[index] = count;
    }

    if (storeP != store) {
      delete[] storeP;
    }

    if (count == 1) {
      fl = out.theData;
      m = 0;
      if (theExclusions.lLength) {
        for (n = 0; n < undimension; n++, fl++) {
          if (m < theExclusions.lLength && n == theExclusions.lData[m]) {
            m++;
            continue;
          }

          if (*fl > 0.) {
            parvect[n - m] = 1.0;
            break;
          } else {
            parvect[n - m] = 0.0;
          }
        }
        if (n < undimension) {
          n -= m - 1;
          parvect += n;
          for (; n < dimension; n++, parvect++) {
            *parvect = 0.0;
          }
          return index - m;
        } else {
          if (smear)
            for (n = 0; n < dimension; n++, parvect++) {
              *parvect = 1.0;
            }
          return -1;
        }

      } else {
        if (count) {
          for (n = 0; n < undimension; n++, fl++) {
            parvect[n] = *fl > 0.0 ? 1.0 : 0.0;
          }
        }
        return index;
      }
    } else {
      XferwCorrection(out, parvect, undimension);

      if (smear) {
        for (count = 0; count < dimension; count++)
          if (parvect[count] > 0.0) {
            break;
          }
        if (count == dimension) {
          for (count = 0; count < dimension; count++) {
            parvect[count] = 1.0;
          }
        }
      }
      return -1;
    }
  }
  return 0;
}

//______________________________________________________________________________
long _DataSetFilter::MapStringToCharIndex(_String &str) {
  long count = 0, nonzeropos = 0, store[HYPHY_SITE_DEFAULT_BUFFER_SIZE];

  if (unitLength == 1) {
    theData->theTT->TokenCode(str.sData[0], store);

    if (theExclusions.lLength == 0) {
      for (long i = 0; i < undimension; i++)
        if (store[i]) {
          nonzeropos = i;
          count++;
        }
    } else {
      long k = 0;
      for (long i = 0; i < undimension; i++) {
        if (i == theExclusions.lData[k] && k < theExclusions.lLength) {
          k++;
        } else if (store[i]) {
          nonzeropos = i;
          count++;
        }
      }
    }

    if (count != 1) {
      return -1;
    }

    return nonzeropos;
  } else {
    long m, n, index = 0, shifter = 1, *lp, *storeP,
               alphabet_dimesion = theData->theTT->Dimension();

    count = 1;

    if (alphabet_dimesion * unitLength >= HYPHY_SITE_DEFAULT_BUFFER_SIZE) {
      storeP = new long[alphabet_dimesion * unitLength];
    } else {
      storeP = store;
    }

    for (m = 0; m < unitLength; m++) {
      theData->theTT->TokenCode(str.sData[m], storeP + alphabet_dimesion * m);
    }

    for (m = unitLength - 1; m >= 0; m--) {
      int smcount = 0;

      lp = storeP + alphabet_dimesion * m;

      for (n = 0; n < alphabet_dimesion; n++, lp++)
        if (*lp) {
          index += shifter * n;
          smcount++;
        }

      if (smcount == 0) {
        lp = storeP + alphabet_dimesion * m;

        for (n = 0; n < alphabet_dimesion; n++, lp++) {
          *lp = 1;
        }
        smcount = alphabet_dimesion;
      }
      shifter *= alphabet_dimesion;
      count *= smcount;
    }

    if (storeP != store) {
      delete[] storeP;
    }

    if (count == 1) {
      m = 0;
      if (theExclusions.lLength) {
        for (long exc = 0; exc < theExclusions.lLength; exc++)
          if (index == theExclusions.lData[exc]) {
            return -1;
          } else if (theExclusions.lData[exc] > index) {
            return index - exc;
          }

        return index - theExclusions.lLength;
      } else {
        return index;
      }
    } else {
      return -1;
    }
  }

  return 0;
}

//______________________________________________________________________________
long _DataSetFilter::Translate2Frequencies(char s, _Parameter *parvect,
                                           bool smear) {
  long count = 0, store[HYPHY_SITE_DEFAULT_BUFFER_SIZE];

  theData->theTT->TokenCode(s, store, smear);

  if (theExclusions.lLength) {
    long k = 0;
    for (long i = 0; i < undimension; i++) {
      if (i == theExclusions[k]) {
        k++;
      } else if (store[i]) {
        count++;
      }
    }
    if (count) {
      XferwCorrection(store, parvect, undimension);
    }
  } else {
    for (long i = 0; i < undimension; i++)
      if ((parvect[i] = store[i])) {
        count++;
      }
  }

  if (count == 0 && smear)
    for (long i = 0; i < undimension; i++) {
      parvect[i] = 1.0;
    }

  return count == 1 ? 1 : -1;
}

//______________________________________________________________________________
long _DataSetFilter::LookupConversion(char s, _Parameter *parvect) {
  if (undimension == 4) {
    //int idx = (s-40)*5;
    long *cCache = conversionCache.lData + (s - 40) * 5;
    /*parvect[0] = conversionCache.lData[idx++];
    parvect[1] = conversionCache.lData[idx++];
    parvect[2] = conversionCache.lData[idx++];
    parvect[3] = conversionCache.lData[idx++];
    return conversionCache.lData[idx];*/
    parvect[0] = *cCache;
    cCache++;
    parvect[1] = *cCache;
    cCache++;
    parvect[2] = *cCache;
    cCache++;
    parvect[3] = *cCache;
    return *(++cCache);

  } else {
    int idx = (s - 40) * (undimension + 1);
    for (long i = 0; i < undimension;
         parvect[i++] = conversionCache.lData[idx++])
      ;
    return conversionCache.lData[idx];
  }
}

//______________________________________________________________________________
void _DataSetFilter::SetupConversion(void) {
  if (conversionCache.lLength) {
    return;
  }

  if (unitLength == 1) { // do stuff
    char c = 40;

    long i, charCount, theCode;

    _Parameter *temp = new _Parameter[undimension + 1];

    while (c < 127) {
      for (i = 0; i < undimension; i++) {
        temp[i] = 0;
      }

      Translate2Frequencies(c, temp, true);

      charCount = -1;
      for (i = 0; i < undimension; i++) {
        theCode = (long) temp[i];
        conversionCache << theCode;
        if (theCode) {
          if (charCount == -1) {
            charCount = i;
          } else {
            charCount = -2;
          }
        }
      }
      conversionCache << charCount;
      c++;
    }
    delete[] temp;
  } else {
    if (unitLength == 2 || unitLength == 3) {
      _String alphabet(16, true);
      alphabet << theData->theTT->RetrieveCharacters();
      alphabet.Finalize();

      long ccache[88], i, k = GetDimension(false);

      conversionCache.RequestSpace(89 + k);
      conversionCache << alphabet.sLength;

      for (i = 0; i < 88; i++) {
        ccache[i] = -1;
      }
      for (i = 0; i < alphabet.sLength; i++) {
        ccache[alphabet.sData[i] - 40] = i;
      }
      for (i = 0; i < 88; i++) {
        conversionCache << ccache[i];
      }

      long *tcache = new long[k];
      checkPointer(tcache);
      //_Parameter *vcache = new _Parameter [k];
      //checkPointer (vcache);

      if (unitLength == 3) {
        _String s(3, false);
        i = 0;
        for (long a = 0; a < alphabet.sLength;
             a++, i += alphabet.sLength * alphabet.sLength) {
          s.sData[0] = alphabet.sData[a];
          for (long b = 0; b < alphabet.sLength; b++) {
            s.sData[1] = alphabet.sData[b];
            for (long c = 0; c < alphabet.sLength; c++) {
              s.sData[2] = alphabet.sData[c];
              //tcache [i+b*alphabet.sLength+c] = Translate2Frequencies
              //(s,vcache,true);
              tcache[i + b * alphabet.sLength + c] = MapStringToCharIndex(s);
            }
          }
        }
      } else {
        _String s(2, false);
        i = 0;
        for (long a = 0; a < alphabet.sLength; a++, i += alphabet.sLength) {
          s.sData[0] = alphabet.sData[a];
          for (long b = 0; b < alphabet.sLength; b++) {
            s.sData[1] = alphabet.sData[b];
            //tcache [i+b] = Translate2Frequencies (s,vcache,true);
            tcache[i + b] = MapStringToCharIndex(s);
          }
        }
      }
      for (i = 0; i < k; i++) {
        conversionCache << tcache[i];
      }

      //delete vcache;
      delete[] tcache;
    }
  }
}

//______________________________________________________________________________
BaseRef _DataSetFilter::toStr(void) {
  //return new _String("DataSetFilters only print to files");
  _String *res = new _String(4096L, true);
  checkPointer(res);
  internalToStr(nil, *res);
  res->Finalize();
  return res;
}

//______________________________________________________________________________
void _DataSetFilter::PatternToSiteMapper(void *source, void *target, char mode,
                                         long padup) {
  for (long site = 0; site < duplicateMap.lLength; site++)
    if (mode == 0) {
      ((_Parameter *)target)[site] =
          ((_Parameter *)source)[duplicateMap.lData[site]];
    } else if (mode == 1) {
      ((long *)target)[site] = ((long *)source)[duplicateMap.lData[site]];
    } else if (mode == 2) {
      ((long *)target)[site] = ((_Parameter *)source)[duplicateMap.lData[site]];
    }

  for (long site = duplicateMap.lLength; site < padup; site++)
    if (mode == 0) {
      ((_Parameter *)target)[site] = 1.;
    } else if (mode == 1) {
      ((long *)target)[site] = 0;
    }
}

//______________________________________________________________________________
long _DataSetFilter::GetOriginalToShortMap(long index) {
  long pos1 = theData->theMap.lData[theOriginalOrder.lData[index]], pos2;
  pos2 = theMap.Find(pos1);
  if (pos2 == -1) {
    for (long i = theData->theMap.lLength - 1; i >= 0; i--) {
      if (theData->theMap.lData[i] == pos1) {
        pos2 = theMap.Find(i);
        if (pos2 != -1) {
          break;
        }
      }
    }
  }
  return pos2;
}

//______________________________________________________________________________
_String _DataSetFilter::GenerateConsensusString(_SimpleList *majority) {
  if (unitLength > 3) {
    return empty;
  }

  _String result((unsigned long) theOriginalOrder.lLength),
      tRes((unsigned long)(unitLength * theFrequencies.lLength));

  long charStates = GetDimension(false),
       *translationBuffer = (long *)MemAllocate(sizeof(long) * charStates);

  _Parameter *countBuffer =
                 (_Parameter *)MemAllocate(sizeof(_Parameter) * charStates),
             nf;

  SetupConversion();

  for (long k = 0; k < theFrequencies.lLength; k++) {
    long m, t = theMap.lData[k], f;

    for (m = 0; m < charStates; m++) {
      countBuffer[m] = 0.0;
    }

    for (long p = 0; p < theNodeMap.lLength; p++) {
      theData->theTT
          ->TokenCode((*theData)(t, theNodeMap.lData[p], 1), translationBuffer);
      f = 0;

      for (m = 0; m < charStates; m++)
        if (translationBuffer[m]) {
          f++;
        }

      if (f > 1) {
        nf = 1. / f;
        for (m = 0; m < charStates; m++)
          if (translationBuffer[m]) {
            countBuffer[m] += nf;
          }
      } else {
        if (f == 1) {
          m = 0;
          while (!translationBuffer[m++])
            ;
          countBuffer[m - 1] += 1.;
        }
      }
    }

    nf = -1;
    f = 1;

    for (m = 0; m < charStates; m++) {
      if (countBuffer[m] > nf) {
        nf = countBuffer[m];
        t = m;
        f = 1;
      } else if (countBuffer[m] == nf) {
        f++;
      }
    }

    if (f > 1) {
      for (m = 0; m < charStates; m++) {
        if (countBuffer[m] == nf) {
          translationBuffer[m] = 1;
        } else {
          translationBuffer[m] = 0;
        }
      }
      tRes.sData[k] = theData->theTT->CodeToLetter(translationBuffer);
      if (majority) {
        //(*majority) << -1;
        (*majority) << nf;
      }
    } else {
      _String conv = theData->theTT->ConvertCodeToLetters(t, 1);
      tRes.sData[k] = conv.sData[0];
      if (majority) {
        (*majority) << nf;
      }
    }
  }

  free(countBuffer);
  free(translationBuffer);

  for (long m = 0; m < theOriginalOrder.lLength; m++) {
    result.sData[m] = tRes.sData[duplicateMap.lData[m]];
  }

  return result;
}

//______________________________________________________________________________
void _DataSetFilter::toFileStr(FILE *dest) {
  // write out the file with this dataset filter
  if (!dest) {
    return;
  }

  _String dummy;
  internalToStr(dest, dummy);
}

//______________________________________________________________________________
void _DataSetFilter::ConvertCodeToLettersBuffered(long code, char unit,
                                                  char *storage,
                                                  _AVLListXL *lookup) {
  // write out the file with this dataset filter
  long lookupC = lookup->Find((BaseRef) code);
  char *lookupV;
  if (lookupC >= 0) {
    lookupV = ((_String *)lookup->GetXtra(lookupC))->sData;
  } else {
    _String *newT = new _String(ConvertCodeToLetters(code, unit));
    lookup->Insert((BaseRef) code, (long) newT, false);
    lookupV = newT->sData;
  }

  for (long k = 0; k < unit; k++) {
    storage[k] = lookupV[k];
  }
}

//______________________________________________________________________________
void _DataSetFilter::internalToStr(FILE *dest, _String &rec) {
  // write out the file with this dataset filter
  checkParameter(dataFilePrintFormat, dFPrintFormat, 6.0);
  checkParameter(dataFileDefaultWidth, dFDefaultWidth, 50.0);
  _Parameter gW;

  long outputFormat = dFPrintFormat, printWidth = dFDefaultWidth, gapWidth;

  checkParameter(dataFileGapWidth, gW, 10.0);
  if (!printWidth) {
    printWidth = 50;
  }

  gapWidth = gW;
  if (gapWidth <= 0) {
    gapWidth = printWidth;
  }

  long i, j;

  if (outputFormat < 4 || outputFormat > 8) {
    if (!theData->theTT->CheckType(HY_TRANSLATION_TABLE_STANDARD_PROTEIN |
                                   HY_TRANSLATION_TABLE_STANDARD_NUCLEOTIDE)) {
      const _String *bSet = theData->theTT->RetrieveCharacters(),
                    *tokens_added = &theData->theTT->RetrieveAddedTokens();
      if (dest) {
        fprintf(dest, "$BASESET:\"%s\"\n", bSet->sData);

        for (long at = 0; at < tokens_added->sLength; at++) {
          char this_token = tokens_added->getChar(at);
          fprintf(dest, "$TOKEN:\"%c\" = \"", this_token);
          long buf[256];
          theData->theTT
              ->SplitTokenCode(theData->theTT->TokenCode(this_token), buf);
          for (long tc = 0; tc < bSet->sLength; tc++)
            if (buf[tc]) {
              fprintf(dest, "%c", bSet->sData[tc]);
            }
          fprintf(dest, "\"\n");
        }
      } else {
        rec << "$BASESET:\"";
        rec << *bSet;
        rec << "\"\n";
        for (long at = 0; at < tokens_added->sLength; at++) {
          char this_token = tokens_added->getChar(at);
          rec << "$TOKEN:\"";
          rec << this_token;
          rec << "\" = \"";
          long buf[256];
          theData->theTT
              ->SplitTokenCode(theData->theTT->TokenCode(this_token), buf);
          for (long tc = 0; tc < bSet->sLength; tc++)
            if (buf[tc]) {
              rec << bSet->sData[tc];
            }
          rec << "\"\n";
        }
      }
    }
  }

  switch (outputFormat) {
  case 1:    // hash-mark interleaved
  case 10: { // FASTA interleaved
    long sitesDone = 0, upTo;
    char seqDelimiter = (outputFormat == 1) ? '#' : '>';

    for (i = 0; i < theNodeMap.lLength; i++) {
      _String *curName = (_String *)theData->GetNames()(theNodeMap.lData[i]);
      if (dest) {
        fprintf(dest, "%c%s\n", seqDelimiter, curName->sData);
      } else {
        rec << *curName;
        rec << '\n';
      }
    }
    while (sitesDone < theOriginalOrder.lLength) {
      if (dest) {
        fprintf(dest, "\n\n");
      } else {
        rec << '\n';
        rec << '\n';
      }

      upTo = sitesDone + printWidth;
      if (upTo > theOriginalOrder.lLength) {
        upTo = theOriginalOrder.lLength;
      }

      if (dest)
        for (i = 0; i < theNodeMap.lLength; i++) {
          for (j = sitesDone; j < upTo; j++) {
            if ((j - sitesDone) % gapWidth == 0) {
              fprintf(dest, " ");
            }
            fprintf(dest, "%c", (*theData)(theOriginalOrder.lData[j],
                                           theNodeMap.lData[i], 1));
          }

          fprintf(dest, "\n");
        }
      else
        for (i = 0; i < theNodeMap.lLength; i++) {
          for (j = sitesDone; j < upTo; j++) {
            if ((j - sitesDone) % gapWidth == 0) {
              rec << ' ';
            }
            rec << (*theData)(theOriginalOrder.lData[j], theNodeMap.lData[i],
                              1);
          }

          rec << '\n';
        }

      sitesDone = upTo;
    }
    break;
  }

  case 2:  // phylip sequential
  case 11: // PAML
           {

    if (dest) {
      fprintf(dest, "%ld\t%ld\n", theNodeMap.lLength, theOriginalOrder.lLength);
    } else {
      rec << _String((long) theNodeMap.lLength);
      rec << '\t';
      rec << _String((long) theNodeMap.lLength, theOriginalOrder.lLength);
      rec << '\n';
    }
    // proceed to spool out the data
    for (i = 0; i < theNodeMap.lLength; i++) {
      _String *curName = (_String *)theData->GetNames()(theNodeMap(i)),
              choppedTo10Chars;
      if (outputFormat == 2) {
        if (curName->Length() >= 10) {
          choppedTo10Chars = curName->Cut(0, 9) & ' ';
        } else {
          choppedTo10Chars = *curName;
          while (choppedTo10Chars.Length() < 11) {
            choppedTo10Chars = choppedTo10Chars & ' ';
          }
        }
      } else {
        choppedTo10Chars = *curName & "  ";
      }

      if (dest) {
        fprintf(dest, "%s", choppedTo10Chars.sData);

        for (j = 0; j < theOriginalOrder.lLength; j++) {
          if ((j % printWidth == 0) && j) {
            fprintf(dest, "\n           ");
          }
          fprintf(dest, "%c",
                  (*theData)(theOriginalOrder(j), theNodeMap(i), 1));
          if (j % gapWidth == gapWidth - 1) {
            fprintf(dest, " ");
          }
        }
        fprintf(dest, "\n");
      } else {
        rec << choppedTo10Chars;

        for (j = 0; j < theOriginalOrder.lLength; j++) {
          if ((j % printWidth == 0) && j) {
            rec << "\n           ";
          }

          rec << (*theData)(theOriginalOrder(j), theNodeMap(i), 1);
          if (j % gapWidth == gapWidth - 1) {
            rec << ' ';
          }
        }
        rec << '\n';
      }

    }
    break;
  }
  case 3: { // phylip interleaved
            // print PHYLIP format header
            //fprintf (dest,"$FORMAT:\"PHYLIPI\"\n");
            // print number of species and sites
    if (dest) {
      fprintf(dest, "%ld\t%ld\n", theNodeMap.lLength, theOriginalOrder.lLength);
    } else {
      rec << _String((long) theNodeMap.lLength);
      rec << '\t';
      rec << _String((long) theNodeMap.lLength, theOriginalOrder.lLength);
      rec << '\n';
    }
    // proceed to spool out the data
    for (i = 0; i < theNodeMap.lLength; i++) {
      _String *curName = (_String *)theData->GetNames()(theNodeMap(i)),
              choppedTo10Chars;
      if (curName->Length() >= 10) {
        choppedTo10Chars = curName->Cut(0, 9) & ' ';
      } else {
        choppedTo10Chars = *curName;
        while (choppedTo10Chars.Length() < 11) {
          choppedTo10Chars = choppedTo10Chars & ' ';
        }
      }

      if (dest) {
        fprintf(dest, "%s", choppedTo10Chars.sData);

        for (j = 0; j < theOriginalOrder.lLength; j++) {
          if (j == printWidth) {
            fprintf(dest, "\n");
            break;
          }
          if (j % gapWidth == 0) {
            fprintf(dest, " ");
          }
          fprintf(dest, "%c", (*theData)(theOriginalOrder.lData[j],
                                         theNodeMap.lData[i], 1));
        }
      } else {
        rec << choppedTo10Chars;

        for (j = 0; j < theOriginalOrder.lLength; j++) {
          if (j == printWidth) {
            rec << '\n';
            break;
          }
          if (j % gapWidth == 0) {
            rec << ' ';
          }
          rec << (*theData)(theOriginalOrder.lData[j], theNodeMap.lData[i], 1);
        }
      }

    }

    long completed = printWidth;

    if (dest) {
      while (completed < theOriginalOrder.lLength - 1) {
        long upTo = completed + printWidth < theOriginalOrder.lLength
                        ? completed + printWidth
                        : theOriginalOrder.lLength;
        for (i = 0; i < theNodeMap.lLength; i++) {
          fprintf(dest, "\n           ");
          for (j = completed; j < upTo; j++) {
            if ((j - completed) % gapWidth == 0) {
              fprintf(dest, " ");
            }
            fprintf(dest, "%c", (*theData)(theOriginalOrder.lData[j],
                                           theNodeMap.lData[i], 1));
          }
        }
        completed += printWidth;
        fprintf(dest, "\n");
      }
    } else {
      while (completed < theOriginalOrder.lLength - 1) {
        long upTo = completed + printWidth < theOriginalOrder.lLength
                        ? completed + printWidth
                        : theOriginalOrder.lLength;
        for (i = 0; i < theNodeMap.lLength; i++) {
          rec << "\n           ";
          for (j = completed; j < upTo; j++) {
            if ((j - completed) % gapWidth == 0) {
              rec << ' ';
            }
            rec << (*theData)(theOriginalOrder.lData[j], theNodeMap.lData[i],
                              1);
          }
        }
        completed += printWidth;
        rec << '\n';
      }
    }

    break;
  }

  // various flavors of NEXUS

  case 4:   // labels, sequential
  case 5:   // labels, interleaved
  case 6:   // no labels, sequential
  case 7: { // no labels, interleaved
            // write out the header
    j = theNodeMap.lLength;
    if (dest) {
      fprintf(dest, "#NEXUS\n\n[\nGenerated by %s on %s\n]\n\nBEGIN "
                    "TAXA;\n\tDIMENSIONS NTAX = %ld;\n\tTAXLABELS\n\t\t",
              GetVersionString().getStr(), GetTimeStamp().getStr(), j);
      for (i = 0; i < j; i++) {
        fprintf(
            dest, "'%s' ",
            ((_String *)theData->GetNames().lData[theNodeMap.lData[i]])->sData);
      }

      fprintf(dest, ";\nEND;\n\nBEGIN CHARACTERS;\n\tDIMENSIONS NCHAR = "
                    "%ld;\n\tFORMAT\n\t\t",
              theOriginalOrder.lLength);
      if (theData->theTT->CheckType(HY_TRANSLATION_TABLE_STANDARD_NUCLEOTIDE)) {
        fprintf(dest, "DATATYPE = DNA\n");
      } else {
        if (theData->theTT->CheckType(HY_TRANSLATION_TABLE_STANDARD_PROTEIN)) {
          fprintf(dest, "DATATYPE = PROTEIN\n");
        } else if (theData->theTT->CheckType(
                       HY_TRANSLATION_TABLE_STANDARD_BINARY)) {
          fprintf(dest, "DATATYPE = BINARY\n");
        } else {
          const _String *bSet = theData->theTT->RetrieveCharacters(),
                        *tokens_added = &theData->theTT->RetrieveAddedTokens();
          fprintf(dest, "\n\tSYMBOLS = \"");
          for (long bc = 0; bc < bSet->sLength - 1; bc++) {
            fprintf(dest, "%c ", bSet->sData[bc]);
          }
          fprintf(dest, "%c\"\n", bSet->sData[bSet->sLength - 1]);
          for (long at = 0; at < tokens_added->sLength; at++) {
            char this_token = tokens_added->getChar(at);
            fprintf(dest, "\n\tEQUATE=\"%c = ", this_token);
            long buf[256];
            theData->theTT
                ->SplitTokenCode(theData->theTT->TokenCode(this_token), buf);
            for (long tc = 0; tc < bSet->sLength; tc++)
              if (buf[tc]) {
                fprintf(dest, "%c", bSet->sData[tc]);
              }
            fprintf(dest, "\"");
          }

        }
      }
      if (theData->theTT->GetGapChar()) {
        fprintf(dest, "\n\t\tGAP=%c", theData->theTT->GetGapChar());
      }
      if (theData->theTT->GetSkipChar()) {
        fprintf(dest, "\n\t\tMISSING=%c", theData->theTT->GetSkipChar());
      }

      if (outputFormat > 5) {
        fprintf(dest, "\n\t\tNOLABELS");
      }
      if (outputFormat % 2) {
        fprintf(dest, "\n\t\tINTERLEAVE");
      }

      fprintf(dest, "\n\t;\n\nMATRIX");
    } else {
      rec << "#NEXUS\n\nBEGIN TAXA;\n\tDIMENSIONS NTAX = ";
      rec << _String((long) j);
      rec << ";\n\tTAXLABELS\n\t\t";

      for (i = 0; i < j; i++) {
        rec << "'";
        rec << ((_String *)theData->GetNames().lData[theNodeMap.lData[i]]);
        rec << "' ";
      }
      rec << ";\nEND;\n\nBEGIN CHARACTERS;\n\tDIMENSIONS NCHAR = ";
      rec << _String((long) theOriginalOrder.lLength);
      rec << ";\n\tFORMAT\n\t\t";

      if (theData->theTT->CheckType(HY_TRANSLATION_TABLE_STANDARD_NUCLEOTIDE)) {
        rec << "DATATYPE = DNA\n";
      } else {
        if (theData->theTT->CheckType(HY_TRANSLATION_TABLE_STANDARD_PROTEIN)) {
          rec << "DATATYPE = PROTEIN\n";
        } else if (theData->theTT->CheckType(
                       HY_TRANSLATION_TABLE_STANDARD_BINARY)) {
          rec << "DATATYPE = BINARY\n";
        } else {
          rec << "\t\tSYMBOLS = \"";
          const _String *bSet = theData->theTT->RetrieveCharacters(),
                        *added_tokens = &theData->theTT->RetrieveAddedTokens();
          for (long bc = 0; bc < bSet->sLength - 1; bc++) {
            rec << bSet->sData[bc];
            rec << ' ';
          }
          rec << bSet->sData[bSet->sLength - 1];
          rec << "\"\n";
          for (long at = 0; at < added_tokens->sLength; at++) {
            char this_token = added_tokens->getChar(at);
            rec << "\nEQUATE =\"";
            rec << this_token;
            rec << " = ";
            long buf[256];
            theData->theTT
                ->SplitTokenCode(theData->theTT->TokenCode(this_token), buf);
            for (long tc = 0; tc < bSet->sLength; tc++)
              if (buf[tc]) {
                rec << bSet->sData[tc];
              }

            rec << "\"";
          }
        }
      }
      if (theData->theTT->GetGapChar()) {
        rec << "\t\tGAP=";
        rec << theData->theTT->GetGapChar();
      }
      if (theData->theTT->GetSkipChar()) {
        rec << "\n\t\tMISSING=";
        rec << theData->theTT->GetSkipChar();
      }
      if (outputFormat > 5) {
        rec << "\n\t\tNOLABELS";
      }
      if (outputFormat % 2) {
        rec << "\n\t\tINTERLEAVE";
      }

      rec << "\n\t;\n\nMATRIX";

    }

    //compute space alignment for different taxa names
    // two passes - one to locate the max length and 2nd to compute padding
    // lengths

    j = 0;
    for (i = 0; i < theNodeMap.lLength; i++) {
      if (((_String *)theData->GetNames()(theNodeMap(i)))->sLength > j) {
        j = ((_String *)theData->GetNames()(theNodeMap(i)))->sLength;
      }
    }

    _SimpleList taxaNamesPadding;

    for (i = 0; i < theNodeMap.lLength; i++) {
      taxaNamesPadding
          << (j - ((_String *)theData->GetNames()(theNodeMap(i)))->sLength);
    }

    if (outputFormat % 2 == 0) { // sequential
      if (dest)
        for (i = 0; i < theNodeMap.lLength; i++) {
          if (outputFormat == 4) { // labels
            fprintf(dest, "\n\t'%s'",
                    ((_String *)theData->GetNames()(theNodeMap(i)))->sData);
            for (j = 0; j <= taxaNamesPadding.lData[i]; j++) {
              fprintf(dest, " ");
            }
          } else {
            fprintf(dest, "\n");
          }
          fprintf(dest, " ");
          for (j = 0; j < theOriginalOrder.lLength; j++) {
            fprintf(dest, "%c", (*theData)(theOriginalOrder.lData[j],
                                           theNodeMap.lData[i], 1));
          }
        }
      else
        for (i = 0; i < theNodeMap.lLength; i++) {
          if (outputFormat == 4) { // labels
            rec << "\n\t'";
            rec << (*(_String *)theData->GetNames()(theNodeMap(i)));
            rec << "'";

            for (j = 0; j <= taxaNamesPadding.lData[i]; j++) {
              rec << ' ';
            }
          } else {
            rec << '\n';
          }
          rec << ' ';
          for (j = 0; j < theOriginalOrder.lLength; j++) {
            rec << (*theData)(theOriginalOrder.lData[j], theNodeMap.lData[i],
                              1);
          }
        }
    } else {
      long sitesDone = 0, upTo;

      while (sitesDone < theOriginalOrder.lLength) {
        upTo = sitesDone + printWidth;
        if (upTo > theOriginalOrder.lLength) {
          upTo = theOriginalOrder.lLength;
        }

        if (dest) {
          for (i = 0; i < theNodeMap.lLength; i++) {
            if (outputFormat == 5) { // labels
              fprintf(dest, "\n\t'%s'",
                      ((_String *)theData->GetNames()(theNodeMap(i)))->sData);
              for (j = 0; j <= taxaNamesPadding.lData[i]; j++) {
                fprintf(dest, " ");
              }
            } else {
              fprintf(dest, "\n");
            }
            fprintf(dest, " ");
            for (j = sitesDone; j < upTo; j++) {
              fprintf(dest, "%c", (*theData)(theOriginalOrder.lData[j],
                                             theNodeMap.lData[i], 1));
            }
          }
          fprintf(dest, "\n\n");
        } else {
          for (i = 0; i < theNodeMap.lLength; i++) {
            if (outputFormat == 5) { // labels
              rec << "\n\t'";
              rec << (*(_String *)theData->GetNames()(theNodeMap(i)));
              rec << "'";
              for (j = 0; j <= taxaNamesPadding.lData[i]; j++) {
                rec << ' ';
              }
            } else {
              rec << '\n';
            }

            rec << ' ';
            for (j = sitesDone; j < upTo; j++) {
              rec << (*theData)(theOriginalOrder.lData[j], theNodeMap.lData[i],
                                1);
            }
          }
          rec << '\n';
          rec << '\n';
        }

        sitesDone = upTo;
      }

    }
    if (dest) {
      fprintf(dest, ";\nEND;");
    } else {
      rec << ";\nEND;";
    }
    break;
  }

  case 8: {
    for (i = 0; i < theNodeMap.lLength; i++) {
      if (dest) {
        fprintf(dest, "%c", (*theData)(theOriginalOrder(0), theNodeMap(i), 1));
        for (j = 1; j < theOriginalOrder.lLength; j++) {
          fprintf(dest, ",%c",
                  (*theData)(theOriginalOrder(j), theNodeMap(i), 1));
        }
        fprintf(dest, "\n");
      } else {
        rec << (*theData)(theOriginalOrder(0), theNodeMap(i), 1);
        for (j = 1; j < theOriginalOrder.lLength; j++) {
          rec << ',';
          rec << (*theData)(theOriginalOrder(j), theNodeMap(i), 1);
        }
        rec << '\n';
      }
    }
    break;
  }

  default: { // hash-mark sequential
    char seqDelimiter = (outputFormat == 9) ? '>' : '#';

    for (i = 0; i < theNodeMap.lLength; i++) {
      _String *curName = (_String *)theData->GetNames()(theNodeMap(i));

      if (dest) {
        fprintf(dest, "%c%s", seqDelimiter, curName->sData);
        for (j = 0; j < theOriginalOrder.lLength; j++) {
          if (j % printWidth == 0) {
            fprintf(dest, "\n");
          }
          fprintf(dest, "%c",
                  (*theData)(theOriginalOrder(j), theNodeMap(i), 1));
        }
        fprintf(dest, "\n");
      } else {
        rec << seqDelimiter;
        rec << *curName;
        for (j = 0; j < theOriginalOrder.lLength; j++) {
          if (j % printWidth == 0) {
            rec << '\n';
          }
          rec << (*theData)(theOriginalOrder(j), theNodeMap(i), 1);
        }
        rec << '\n';
      }
    }
  }

    // finally see if we need to write out a tree

  }

  if (outputFormat != 8) {
    _Parameter treeDefined;
    checkParameter(dataFileTree, treeDefined, 0.0);
    if (treeDefined) {
      _Variable *treeVar = FetchVar(LocateVarByName(dataFileTreeString));
      if (treeVar) {
        _String *treeString = (_String *)(treeVar->Compute())->toStr();
        switch (outputFormat) {
        case 0:
        case 1:
        case 9:
        case 10: {
          if (dest) {
            fprintf(dest, "\n\n%s;", treeString->sData);
          } else {
            rec << '\n';
            rec << '\n';
            rec << *treeString;
          }
          break;
        }
        case 2:
        case 3: {
          if (dest) {
            fprintf(dest, "\n1\n%s;", treeString->sData);
          } else {
            rec << "\n1\n";
            rec << *treeString;
          }
          break;
        }
        default: {
          if (dest) {
            fprintf(dest, "\n\nBEGIN TREES;\n\tTREE tree = %s;\nEND;",
                    treeString->sData);
          } else {
            rec << "\n\nBEGIN TREES;\n\tTREE tree = ";
            rec << *treeString;
            rec << ";\nEND;";
          }
        }
        }
        DeleteObject(treeString);
      }
    }
  }

}

void _DataSetFilter::SetMap(_String &s) {
  theNodeMap.Clear();
  if (s.Length()) {
    long f, g = 0;
    //_String sc(",");
    f = s.Find(',');
    while (f != -1) {
      theNodeMap << (long) s.Cut(g, f - 1).toNum();
      g = f + 1;
      f = s.Find(',', f + 1, -1);
    }
    theNodeMap << (long) s.Cut(g, -1).toNum();
  }
}
