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

#include "hy_globals.h"
#include "helperfunctions.h"

#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include "likefunc.h"
#include "datasetfilter.h"
#include "executionlist.h"
#include "math.h"

#if !defined __UNIX__ || defined __HEADLESS__
#include "preferences.h"
#endif

//______________________________________________________________________________

_String dataFileTree("IS_TREE_PRESENT_IN_DATA"),
    dataFileTreeString("DATAFILE_TREE"),
    nexusFileTreeMatrix("NEXUS_FILE_TREE_MATRIX"),
    dataFilePartitionMatrix("DATA_FILE_PARTITION_MATRIX"),
    useTraversalHeuristic("USE_TRAVERSAL_HEURISTIC"),
    defaultLargeFileCutoff("USE_MEMORY_SAVING_DATA_STRUCTURES"), fileTreeString;

//______________________________________________________________________________
/* function declarations */
void checkTTStatus(FileState *fs);
void processCommand(_String *s, FileState *fs);
void FilterRawString(_String &s, FileState *fs, _DataSet &ds);
long ProcessLine(_String &s, FileState *fs, _DataSet &ds);
void PadLine(FileState &fState, _DataSet &result);
void ISelector(FileState &fState, _String &CurrentLine, _DataSet &result);
bool SkipLine(_String &theLine, FileState *fS);
void TrimPhylipLine(_String &CurrentLine, _DataSet &ds);
void ProcessTree(FileState *, FILE *, _String &);
void ReadNexusFile(FileState &fState, FILE *f, _DataSet &result);

//______________________________________________________________________________
_Site::_Site(void) : _CString(16, true) { refNo = -1; }

//______________________________________________________________________________
_Site::_Site(_String &s) : _CString(s.sLength, true) {
  refNo = -1;
  (*this) << &s;
}

//______________________________________________________________________________
_Site::_Site(char s) : _CString(16, true) {
  refNo = -1;
  (*this) << s;
}

//______________________________________________________________________________
_Site::_Site(long s)
    //:_CString (1, true)
    {
  SetRefNo(s);
}

//______________________________________________________________________________
_Site::~_Site(void) {}

//______________________________________________________________________________
void _Site::Complete(void) {
  if (refNo == -1) {
    _String::Finalize();
  }

  refNo = refNo < 0 ? -refNo : refNo;
}

//______________________________________________________________________________
BaseRef _Site::makeDynamic(void) {
  _Site *r = new _Site;
  checkPointer(r);

  memcpy((char *)r, (char *)this, sizeof(_Site));
  r->nInstances = 1;
  nInstances++;
  return r;
}

//______________________________________________________________________________
void _Site::Duplicate(BaseRef ref) {
  _Site *s = (_Site *)ref;
  sLength = s->sLength;
  if (sData) {
    free(sData);
  }
  sData = s->sData;
  allocatedSpace = s->allocatedSpace;
  //nInstances = ref->nInstances;
  if (sData) {
    /*long theLength = sLength/storageIncrement;
    if (!sLength||sLength%storageIncrement) theLength++;
    theLength*=storageIncrement;
    checkPointer (sData = (char*)MemAllocate (theLength));
    memcpy (sData, s->sData, sLength);*/
    if (allocatedSpace) {
      checkPointer(sData = (char *)MemAllocate(allocatedSpace * sizeof(char)));
    } else {
      checkPointer(sData = (char *)MemAllocate(sLength * sizeof(char)));
    }
    memcpy(sData, s->sData, sLength);
  }
  refNo = -1;
}

//______________________________________________________________________________
void _Site::Clear(void) {
  if (sData) {
    free(sData);
    sData = NULL;

    //nInstances = 0;
  }
  allocatedSpace = 0;
  sLength = 0;
}

//______________________________________________________________________________
void _Site::PrepareToUse(void) {
  if (IsCompressed()) {
    _String *s = Decompress();
    DuplicateErasing(s);
    ;
    DeleteObject(s);
    SetDecompressed();
  }
}

//______________________________________________________________________________
void _Site::Archive(void) {
  if ((!IsCompressed()) && (GetRefNo() >= 0)) {
    BestCompress(NUCLEOTIDEALPHABET);
  }
}

//______________________________________________________________________________
// reading the data set file in here
// check whether the translation table needs to be refreshed
void checkTTStatus(FileState *fs) {
  if (fs->translationTable == &defaultTranslationTable) {
    fs->translationTable =
        (_TranslationTable *)defaultTranslationTable.makeDynamic();
  }
}

//______________________________________________________________________________
void processCommand(_String *s, FileState *fs) {
  // loop thru understood values of commands
  static _List CommandList;
  if (CommandList.lLength == 0)
      // first time in, should init commands
      {
    _String command("BASESET");
    CommandList &&&command;
    command = "FORMAT";
    CommandList &&&command;
    command = "RAWLINE";
    CommandList &&&command;
    command = "REPEAT";
    CommandList &&&command;
    command = "TOKEN";
    CommandList &&&command;
  }

  long f = -1;
  long i, k = 0, l = 0, m;
  for (i = 0;(i < CommandList.lLength); i++) {
    f = s->Find(*(_String *)CommandList(i));
    if (f != -1) {
      break;
    }
  }

  if (f == -1) { // unrecognized command
    return;
  } else {
    // trim the string
    s->Trim(f + ((_String *)CommandList(i))->Length(), -1);
    f = s->Find(":");
    if (f == -1) { // poorly formed command
      return;
    } else {
      s->Trim(f + 1, -1);
    }

    if ((i >= 1) && (i <= 3)) {
      k = s->Find('\"');
      if (k == -1) {
        return;
      }
      l = s->Find('\"', k + 1, -1);
      if (l == -1) {
        return;
      }
      if (l <= k) {
        return;
      }
      s->Trim(k + 1, l - 1);
    }

    switch (i) {
      char c;
    case 4: // new token
      checkTTStatus(fs);
      // attempt to extract a token. Looking for (e.g):   "c" = "AC"
      k = s->Find('"');
      if (k == -1) {
        return;
      }
      if ((*s)[k + 2] != '"') {
        return;
      }
      l = s->Find('"', k + 3, -1);
      m = s->Find('"', l + 1, -1);
      if ((l == -1) || (m == -1)) {
        return;
      }

      c = (*s)[k + 1];
      s->Trim(l + 1, m - 1);
      fs->translationTable->AddTokenCode(c, *s);
      break;

    case 0: // new code set, e.g  "ACGU"
      checkTTStatus(fs);
      // erase previous char definitions
      fs->translationTable->Clear();
      if (*s != _String("BASE20")) {
        fs->translationTable->AddBaseSet(*s);
      } else {
        fs->translationTable
            ->SetStandardType(HY_TRANSLATION_TABLE_STANDARD_PROTEIN);
      }
      break;

    case 1:                           //FORMAT
      if (*s == _String("PHYLIPI")) { // PHYLIP Interleaved
        fs->fileType = 1;
        fs->interleaved = TRUE;
      } else if (*s == _String("PHYLIPS")) { // PHYLIP sequential
        fs->fileType = 1;
        fs->interleaved = FALSE;
      }
      if (*s == _String("RAW")) { // RAW Sequential Data (as in NEXUS)
        fs->fileType = 2;
        fs->interleaved = FALSE;
      }
      fs->autoDetect = false;
      break;

    case 3: // REPEAT CHAR
      fs->repeat = s->getChar(0);
      break;

    case 2:
        // RAWLINE template e.g 1,-1 skips one word at the beginning and one
        // word at the end
      _List chips(s, ',');
      for (int i = 0; i < chips.lLength; i++) {
        fs->rawLinesFormat << (long)(((_String *)chips(i))->toNum());
      }

    }
  }
}

//______________________________________________________________________________
void FilterRawString(_String &s, FileState *fs, _DataSet &ds) {
  int i;
  for (i = 0; i < fs->rawLinesFormat.lLength; i++) {
    long f = fs->rawLinesFormat(i), p = 0, l = 0;
    if (f > 0) {
      for (int j = 0;(j < f) && (p >= 0) && (l >= 0); j++) {
        p = s.FirstNonSpaceIndex(l, -1, 1);
        l = s.FirstSpaceIndex(p, -1, 1);
      }
      if (l < 0) {
        break;
      }
      p = s.FirstNonSpaceIndex(l, -1, 1);
      s.Trim(p, -1);
    } else {
      if (f != 0) {
        p = 0;
        l = 0;
        for (int j = 0;(j > f) && (p >= 0) && (l >= 0); j--) {
          p = s.FirstNonSpaceIndex(p, -1, -1);
          l = s.FirstSpaceIndex(0, p, -1);
        }
        if (l < 0) {
          break;
        }
        p = s.FirstNonSpaceIndex(0, l, -1);
        s.Trim(0, p);
      } else {
        // Name
        p = s.FirstNonSpaceIndex();
        l = s.FirstSpaceIndex(p + 1, -1, 1);
        if ((p < 0) || (l < 0)) {
          break;
        }
        _String Name = s.Cut(p, l - 1);
        ds.AddName(Name);
        s.Trim(s.FirstNonSpaceIndex(l, -1, 1), -1);
      }
    }

  }
  if (i != fs->rawLinesFormat.lLength) {
    s = "";
  }
}

//______________________________________________________________________________
void ProcessTree(FileState *fState, FILE *f, _String &CurrentLine) {
  long j = 0, i = 0; // parenthesis balance
  char c;
  _String treeString((unsigned long) 10, true);
  do {
    for (i = 0; i < CurrentLine.sLength; i++) {
      c = CurrentLine.sData[i];
      if (!isspace(c)) {
        treeString << c;
        if (c == ')') {
          j--;
          if (!j) {
            break;
          }
        } else if (c == '(') {
          j++;
        }
      }
    }
    ReadNextLine(f, &CurrentLine, fState, false);
  } while (j && CurrentLine.sLength);

  if (j) {
    _String errMsg(
        "Tree string found in data file had unbalanced parentheses: ");
    if (j > 0) {
      errMsg = errMsg & j & " too few closing parentheses.";
    } else {
      errMsg = errMsg & (-j) & " too many closing parentheses.";
    }
    ReportWarning(errMsg);
  } else {
    treeString.Finalize();
    setParameter(dataFileTree, 1.0, fState->theNamespace);
    setParameter(dataFileTreeString, new _FString(treeString), false);
  }

}

//______________________________________________________________________________
long ProcessLine(_String &s, FileState *fs, _DataSet &ds) {
  long sitesAttached = 0, sL = s.Length();

  for (long l = 0; l < sL; l++) {
    // see if it is a legal char
    char c = toupper(s.sData[l]);
    if (fs->translationTable->IsCharLegal(c)) { 
      if (fs->curSpecies == 0) {                
        // add new column
        ds.AddSite(c);
        sitesAttached++;
      } else { 
        //append to exisiting column
         //if (c == fs->skip) continue;
         // check to see if this species needs to be padded
        if (c == fs->repeat) {
          if (fs->curSite + sitesAttached >= ds.lLength) { 
            // a dot not matched by a previously read character;
            // ignore
            return sitesAttached;
          }

          c = ((_Site *)(ds._List::operator()(fs->curSite + sitesAttached)))
              ->getChar(0);
          if (c == 0)
            c = ((_Site *)(ds._List::operator()(((_Site *)(ds._List::operator()(
                fs->curSite + sitesAttached)))->GetRefNo())))->getChar(0);
        }

        if (fs->curSite + sitesAttached + 1 > fs->totalSitesRead) {
          // pad previous species to full length
          _Site *newS = new _Site(fs->skip);
          checkPointer(newS);
          for (long j = 1; j < fs->curSpecies; j++) {
            (*newS) << fs->skip;
          }

          (*newS) << c;

          /*long rN = ds.dsh->incompletePatterns->Find (newS);

                    if (rN>=0)
                    {
                        rN =  ds.dsh->incompletePatterns->GetXtra
          (rN);
                        ds.theFrequencies[rN]++;
                        newS->Clear();
                        newS->SetRefNo(rN);
                        ds.theFrequencies << 0;
                    }
                    else
                    {*/
          ds.theFrequencies << 1;
          newS->SetRefNo(-1);
          //}

          ds << newS;
          newS->nInstances--;

          fs->totalSitesRead++;
        } else {
          ds.Write2Site(fs->curSite + sitesAttached, c);
        }

        sitesAttached++;
      }
    }
  }

  // make sure that this species has enough data in it, and if not - pad it with
  // '?'
  if ((fs->curSite + sitesAttached < fs->totalSitesRead) && (fs->interleaved)) {
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

//______________________________________________________________________________
// make sure that there is enough data in this line
// and if not - "pad" it with '?''s
void PadLine(FileState &fState, _DataSet &result) {
  if (fState.curSite < fState.totalSitesRead) // pad line if needed
    for (long j = fState.curSite; j < fState.totalSitesRead; j++) {
      result.Write2Site(j, fState.skip);
    }
}

//______________________________________________________________________________
void ISelector(FileState &fState, _String &CurrentLine, _DataSet &result) {
  if (fState.interleaved) {                 
    // interleaved file
    if (fState.curSpecies && (!((fState.curSpecies) % fState.totalSpeciesExpected))) { 
      // read a chunk of all species
      if (fState.totalSitesRead && !result.InternalStorageMode()) {
        for (long i = fState.curSite; i < fState.totalSitesRead; i++) {
          result.Compact(i);
        }

        result.ResetIHelper();

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

//______________________________________________________________________________
bool SkipLine(_String &theLine, FileState *fS) {

  if (theLine.sData[0] == '/' && theLine.sData[1] == '/') {
    return true;
  }

  char c = theLine.FirstNonSpace();

  if (c && (!((c == '$') && (!fS->acceptingCommands)))) {
    return false;
  }

  return true;
}

//______________________________________________________________________________
#define READ_NEXT_LINE_BUFFER_SIZE 1024 * 1024

//______________________________________________________________________________
void ReadNextLine(FILE *fp, _String *s, FileState *fs, bool, bool upCase) {
  _String tempBuffer(1024L, true);

  char lastc;

  if (fp) {
    lastc = fgetc(fp);
  } else {
    lastc =
        fs->pInSrc < fs->theSource->sLength ? fs->theSource->sData[fs->pInSrc++]
                                            : 0;
  }

  if (fs->fileType != 3) { // not NEXUS - do not skip [..]
    if (fp)
      while (!feof(fp) && lastc != 10 && lastc != 13) {
        if (lastc) {
          tempBuffer << lastc;
        }
        lastc = fgetc(fp);
      }
    else
      while (lastc && lastc != 10 && lastc != 13) {
        tempBuffer << lastc;
        lastc = fs->theSource->sData[fs->pInSrc++];
      }
  } else {
    if (upCase) {
      lastc = toupper(lastc);
    }

    while (((fp && (!feof(fp))) ||
            (fs->theSource && (fs->pInSrc <= fs->theSource->sLength))) &&
           lastc != 10 && lastc != 13) {

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
          tempBuffer << ' ';
        }
      } else {
        tempBuffer << lastc;
      }

      if (fp) {
        if (upCase) {
          lastc = toupper(fgetc(fp));
        } else {
          lastc = fgetc(fp);
        }
      } else {
        if (upCase) {
          lastc = toupper(fs->theSource->sData[fs->pInSrc++]);
        } else {
          lastc = fs->theSource->sData[fs->pInSrc++];
        }
      }

    }

    if (lastc == 10 || lastc == 13) {
      tempBuffer << ' ';
    }
  }

  tempBuffer.Finalize();

  if ((fp && feof(fp)) ||
      (fs->theSource && fs->pInSrc >= fs->theSource->sLength))
    if (tempBuffer.sLength == 0) {
      *s = "";
      return;
    }

  if (s->nInstances > 1) {
    *s = tempBuffer;
  } else {
    Ptr saveData = s->sData;
    s->sData = tempBuffer.sData;
    tempBuffer.sData = saveData;

    s->sLength = tempBuffer.sLength;
  }

  if (SkipLine(*s, fs)) {
    ReadNextLine(fp, s, fs, false, upCase);
  }

  if (s->sLength && s->sData[s->sLength - 1] == '\n') {
    s->Trim(0, s->sLength - 2);
  }
}

//______________________________________________________________________________
void TrimPhylipLine(_String &CurrentLine, _DataSet &ds) {

  int fNS = CurrentLine.FirstNonSpaceIndex(),
      space2 = CurrentLine.FirstSpaceIndex(fNS + 1);

  // hack for PAML support
  if (space2 > fNS && isspace(CurrentLine.getChar(space2 + 1))) {
    _String sequence_name(CurrentLine, fNS, space2);
    CurrentLine.Trim(space2 + 2, -1); // chop out the name
    ds.AddName(sequence_name);
  } else {
    _String sequence_name(CurrentLine, fNS, fNS + 9);
    CurrentLine.Trim(fNS + 10, -1); // chop out the name
    ds.AddName(sequence_name);
  }
}

//______________________________________________________________________________
_DataSet *ReadDataSetFile(FILE *f, char execBF, _String *theS, _String *bfName,
                          _String *namespaceID, _TranslationTable *dT) {

  bool doAlphaConsistencyCheck = true;
  _String::storageIncrement = 16;
  _DataSet *result = new _DataSet;
  fileTreeString = "";

  _String CurrentLine = dataFilePartitionMatrix & "={{}};", savedLine;

  if (1) {
    _ExecutionList reset(CurrentLine);
    reset.Execute();
#ifdef __HYPHYMPI__
    if (_hy_mpi_node_rank == 0)
#endif
      terminateExecution = false;
  }

  // initialize the instance of a file state variable
  setParameter(dataFileTree, 0.0);
  FileState fState;
  fState.translationTable = dT;
  fState.curSpecies = fState.totalSpeciesRead = fState.totalSitesRead =
      fState.totalSpeciesExpected = fState.totalSitesExpected = fState.curSite =
          fState.maxStringLength = 0;
  fState.acceptingCommands = true;
  fState.allSpeciesDefined = false;
  fState.interleaved = false;
  fState.isSkippingInNEXUS = false;
  fState.autoDetect = true;
  fState.fileType = -1;
  fState.baseLength = 4;
  fState.repeat = '.', fState.skip = 0;
  fState.theSource = theS;
  fState.pInSrc = 0;
  fState.theNamespace = namespaceID;

  if (!(f || theS)) {
    CurrentLine = "ReadDataSetFile received null file AND string references. "
                  "At least one must be specified";
    warnError(CurrentLine);
  }

  // done initializing
  long fileLength = 0, lastDone = 10, cDone;

#ifdef __HYPHYMPI__
  if (_hy_mpi_node_rank == 0) {
#endif
    if (f) {
      fseek(f, 0, SEEK_END);
      fileLength = ftell(f);
      rewind(f);
    } else {
      fileLength = theS->sLength;
    }

#ifdef __HYPHYMPI__
  }
#endif

  //if (f==NULL) return (_DataSet*)result.makeDynamic();
  // nothing to do
  CurrentLine = empty;

  ReadNextLine(f, &CurrentLine, &fState);
  if (!CurrentLine.sLength) {
    CurrentLine = "Empty File Encountered By ReadDataSet.";
    WarnError(CurrentLine);
    return result;
  } else {
    if (CurrentLine.beginswith("#NEXUS", false)) {
      ReadNexusFile(fState, f, (*result));
      doAlphaConsistencyCheck = false;
    } else {
      long i, j, k, filePosition = -1, saveSpecExpected;
      char c;
      while (CurrentLine.sLength) { // stuff to do
                                    // check if the line has a command in it
#ifdef __HYPHYMPI__
        if (_hy_mpi_node_rank == 0) {
#endif
          if (f) {
            cDone = ftell(f) * 100. / fileLength;
          } else {
            cDone = fState.pInSrc * 100. / fileLength;
          }

          if (cDone > lastDone) {
            SetStatusBarValue(lastDone, 1, 0);
#ifdef __MAC__
            handleGUI(true);
#endif
            lastDone += 10;
          }
#ifdef __HYPHYMPI__
        }
#endif

        c = CurrentLine.FirstNonSpace();
        while (1) {
          if (fState.acceptingCommands) {
            if (c == '$') { 
              // command line
              processCommand(&CurrentLine, &fState);
              break;
            }
          }

          if (!fState.skip) {
            fState.skip = fState.translationTable->GetSkipChar();
          }
          fState.acceptingCommands = FALSE;

          if (fState.fileType == -1) { 
            // undecided file type - assume it is PHYLIP sequential
            if ((c == '#') || (c == '>')) { 
              // hash-mark format
              fState.fileType = 0;
            } else { 
              // assume this is a sequential PHYLIP file
              fState.fileType = 1;
              fState.interleaved = false;
            }

          }

          // decide what to do next
          // if format is PHYLIP and we do not know the expected dimensions,
          //   we must read those in first
          if (fState.fileType == 1) { // PHYLIP
            if ((filePosition < 0) && (fState.autoDetect)) {
              filePosition = (f ? ftell(f)
#ifdef __WINDOZE__
                                      - 1
#endif
                                : fState.pInSrc);
              savedLine = CurrentLine;
            }

            if ((fState.totalSitesExpected == 0) ||
                (fState.totalSpeciesExpected == 0)) { 
              // must read dimensions first
              i = CurrentLine.FirstNonSpaceIndex();
              j = CurrentLine.FirstSpaceIndex(i, -1, 1);
              if (j >= 0) {
                k = CurrentLine.FirstNonSpaceIndex(j, -1, 1);
                if (k >= 0) { // could have dimensions
                  saveSpecExpected = fState.totalSpeciesExpected =
                      CurrentLine.Cut(i, j - 1).toNum();
                  fState.totalSitesExpected = CurrentLine.Cut(k, -1).toNum();
                }
                if (CurrentLine.Find('I', k, -1) >= 0) { // interleaved
                  fState.interleaved = true;
                }
              }
            } else {
              // now for the data crunching part
              // detect a line, diagnose it and dispatch accordingly
              if (fState.interleaved) {
                if (fState.totalSpeciesRead < fState.totalSpeciesExpected) {
                  TrimPhylipLine(CurrentLine, (*result));
                }
                if ((fState.curSite) &&
                    (fState.curSpecies >= saveSpecExpected) &&
                    (fState.totalSitesRead >= fState.totalSitesExpected)) {
                  // reached the end of the data - see maybe there is a tree
                  ReadNextLine(f, &CurrentLine, &fState);
                  if (CurrentLine.sLength) {
                    if (CurrentLine.FirstNonSpace() == '(') { 
                      // could be a tree string
                      ProcessTree(&fState, f, CurrentLine);
                    }
                  }
                  break;
                }

              } else {
                if (fState.totalSitesRead > fState.totalSitesExpected) {
                  // oops - autodetect incorrectly assumed that the file was
                  // sequential
                  fState.curSpecies = fState.totalSpeciesRead =
                      fState.totalSitesRead = fState.curSite = fState
                          .totalSpeciesExpected = fState.totalSitesExpected =
                              fState.maxStringLength = 0;

                  fState.allSpeciesDefined = false;
                  fState.interleaved = true;
                  fState.autoDetect = true;

                  if (f) {
                    fseek(f, filePosition, SEEK_SET);
                  } else {
                    fState.pInSrc = filePosition;
                  }

                  CurrentLine = savedLine;
                  for (long idx = 0; idx < (*result).lLength; idx++) {
                    ((_Site *)(*result).lData[idx])->Finalize();
                  }

                  (*result).theNames.Clear();
                  (*result).theMap.Clear();
                  (*result).Clear();
                  (*result).theFrequencies.Clear();
                  if ((*result).dsh) {
                    (*result).dsh->incompletePatterns->Clear(false);
                    delete ((*result).dsh);
                    (*result).dsh = nil;
                  }
                  continue;
                }

                if (fState.totalSpeciesRead == 0) {
                  fState.totalSpeciesExpected = 1;
                  if (!fState.curSite) {
                    TrimPhylipLine(CurrentLine, (*result));
                  }
                } else if (fState.curSite >= fState.totalSitesExpected) {
                  fState.totalSpeciesExpected++;
                  if (fState.totalSpeciesExpected > saveSpecExpected) {
                    // reached the end of the data - see maybe there is a tree
                    ReadNextLine(f, &CurrentLine, &fState);
                    if (CurrentLine.sLength) {
                      if (CurrentLine.FirstNonSpace() == '(') { 
                        // could be a tree string
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
            if (CurrentLine.sLength) {
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
          if ((c == '#') || (c == '>')) { // a name it is
            if (fState.allSpeciesDefined) { 
              // can't define the species after data
              break;
            } else {
              if ((!fState.totalSpeciesRead) &&
                  (fState.totalSpeciesExpected >= 1)) {
                fState.interleaved = TRUE;
              } else {
                fState.interleaved = FALSE;
              }
              fState.totalSpeciesExpected++;
              CurrentLine.Trim(CurrentLine.FirstNonSpaceIndex(1, -1, 1), -1);
              if ((CurrentLine.sData[0] == '#') ||
                  (CurrentLine.sData[0] == '>')) {
                CurrentLine =
                    _String("Species") & _String(fState.totalSpeciesExpected);
              }
              (*result).AddName(CurrentLine);
            }
            break;
          }

          // check to see if the string defines a tree
          if (c == '(') {
            ProcessTree(&fState, f, CurrentLine);
            ReadNextLine(f, &CurrentLine, &fState);
          }

          // check to see where to stick the incoming line
          if (!fState.totalSpeciesExpected) {
            // raw data fed before names defined - skip
            break;
          }

          if ((fState.totalSpeciesExpected > 1) && (!fState.totalSpeciesRead)) {
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

  if (fState.totalSitesRead && fState.interleaved && !result->InternalStorageMode()) {
    for (int i = fState.curSite; i < fState.totalSitesRead; i++) {
      (*result).Compact(i);
    }
    (*result).ResetIHelper();
  }

  if ((!fState.interleaved) && (fState.fileType != 2)) {
    PadLine(fState, (*result));
  }

#ifdef __HYPHYMPI__
  if (_hy_mpi_node_rank == 0) {
#endif
    SetStatusBarValue(-1, 1, 0);
#ifdef __MAC__
    handleGUI(true);
#endif
#ifdef __HYPHYMPI__
  }
#endif

  // make sure interleaved duplications are handled correctly
  (*result).Finalize();
  _String::storageIncrement = 64;
  (*result).noOfSpecies = fState.totalSpeciesRead;
  (*result).theTT = fState.translationTable;

  // check to see if result may be an amino-acid data
  if (doAlphaConsistencyCheck && result->theTT == &defaultTranslationTable) {
    if (result->GetNoTypes() == 0) {
      // empty data set
      // try binary data
      _TranslationTable *trialTable =
          new _TranslationTable(defaultTranslationTable);
      trialTable->SetStandardType(HY_TRANSLATION_TABLE_STANDARD_BINARY);
      _DataSet *res2 =
          ReadDataSetFile(f, execBF, theS, bfName, namespaceID, trialTable);
      if (res2->GetNoTypes()) {
        DeleteObject(result);
        return res2;
      }
      DeleteObject(trialTable);
    } else
        // check it out
        if (result->CheckAlphabetConsistency() < 0.5)
        // less than 50% of the data in the alphabet is not in the basic
        // alphabet
        {
      _TranslationTable *trialTable =
          new _TranslationTable(defaultTranslationTable);
      trialTable->SetStandardType(HY_TRANSLATION_TABLE_STANDARD_PROTEIN);
      result->theTT = trialTable;
      if ((*result).CheckAlphabetConsistency() < 0.5) {
        CurrentLine =
            "More than 50% of characters in the data are not in the alphabet.";
        (*result).theTT = &defaultTranslationTable;
        DeleteObject(trialTable);
        ReportWarning(CurrentLine);
      }
    }

  }
  if (nexusBFBody.sLength) {
    if (execBF == 1) {
      lastNexusDataMatrix = result;

      long bfl = batchLanguageFunctions.lLength;

      _ExecutionList nexusBF(nexusBFBody, namespaceID);
      if (bfName) {
        nexusBF.sourceFile = *bfName;
      }

#ifndef __UNIX__
#ifdef __HYPHYMPI__
      if (_hy_mpi_node_rank == 0)
#endif
        ApplyPreferences();
#endif

      nexusBF.ExecuteAndClean(bfl);
      nexusBFBody = empty;
    } else if (execBF == 0) {
      nexusBFBody = empty;
    }
  }

  return result;
}

//______________________________________________________________________________
bool StoreADataSet(_DataSet *ds, _String *setName) {
  if (!setName->IsValidIdentifier(true)) {
    WarnError(*setName &
              " is not a valid identifier while constructing a DataSet");
    return false;
  }

  long pos;
  
  _DataSet *existingDS = _HY2DATASET (_HYRetrieveBLObjectByNameFixedType (*setName,
          HY_BL_DATASET, &pos, false, false));
  
      FindDataSetName(*setName);

  if (!existingDS) {
    dataSetNamesList << setName;
    dataSetList.AppendNewInstance(ds);
  } else {
#if !defined __UNIX__ && !defined __HEADLESS__
    if (!RequestDataSetReplace(pos)) {
      terminateExecution = true;
      DeleteObject(ds);
      return false;
    }
#endif

    bool isDifferent =
        existingDS->NoOfSpecies() != ds->NoOfSpecies() ||
        existingDS->NoOfColumns() != ds->NoOfColumns() ||
        existingDS->NoOfUniqueColumns() != ds->NoOfUniqueColumns() ||
        existingDS->GetTT() != ds->GetTT();

    for (long dfIdx = 0; dfIdx < dataSetFilterNamesList.lLength; dfIdx++)
      if (((_String *)dataSetFilterNamesList(dfIdx))->sLength) {
        _DataSetFilter *aDF = _HY2DATASETFILTER (dataSetFilterList(dfIdx));
        if (aDF->GetData() == existingDS) {
          if (isDifferent) {
            ReportWarning(_String("Overwriting dataset '") & *setName &
                          "' caused DataSetFilter '" &
                          *((_String *)dataSetFilterNamesList(dfIdx)) &
                          "' to be deleted");
            KillDataFilterRecord(dfIdx, false);
          } else {
            aDF->SetData(ds);
          }
        }
      }
    dataSetList.Replace(pos, ds, false);
  }

  CheckReceptacleAndStore(*setName & ".species", empty, false, new _Constant(ds->NoOfSpecies()), false);
  CheckReceptacleAndStore(*setName & ".sites", empty, false, new _Constant(ds->NoOfColumns()), false);
  CheckReceptacleAndStore(*setName & ".unique_sites", empty, false, new _Constant(ds->NoOfUniqueColumns()), false);
  return true;
}
