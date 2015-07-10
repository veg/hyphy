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

#include "site.h"
#include "string.h"
#include "ctype.h"
#include "stdlib.h"
#include "list.h"
#include "batchlan.h"

#include "math.h"
#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif
//_________________________________________________________

_DataSet* lastNexusDataMatrix = nil;
_String   nexusBFBody;


void    checkTTStatus               (FileState* fs);
void    processCommand              (_String*s, FileState*fs);
void    FilterRawString             (_String& s, FileState* fs, _DataSet & ds);
long    ProcessLine                 (_String&s , FileState *fs, _DataSet& ds);
void    PadLine                     (FileState& fState, _DataSet& result);
void    ISelector                   (FileState& fState, _String& CurrentLine, _DataSet& result);
bool    SkipLine                    (_String& theLine, FileState* fS);
void    TrimPhylipLine              (_String& CurrentLine, _DataSet& ds);
void    ReadNexusFile               (FileState& fState, FILE*f, _DataSet& result);
bool    ProcessNexusData            (FileState&, long,  FILE*, _String&, _DataSet&);
void    ProcessNexusHYPHY           (FileState&, long,  FILE*, _String&, _DataSet&);
void    ProcessNexusAssumptions     (FileState&, long,  FILE*, _String&, _DataSet&);
void    ProcessNexusTaxa            (FileState&,long, FILE*, _String&, _DataSet&);
void    ProcessNexusTrees           (FileState&, long, FILE*, _String&, _DataSet&);
bool    FindNextNexusToken          (FileState& fState, FILE* f, _String& CurrentLine, long pos);
bool    SkipUntilNexusBlockEnd      (FileState& fState, FILE* f, _String& CurrentLine, long pos);
bool    ReadNextNexusStatement      (FileState&, FILE* , _String&, long, _String&, bool, bool = true, bool = true, bool = false, bool = false, bool = false);
long    ReadNextNexusEquate         (FileState&, FILE* , _String&, long, _String&, bool = false, bool = true);
void    NexusParseEqualStatement    (_String&);


//_________________________________________________________

bool    FindNextNexusToken (FileState& fState, FILE* f, _String& CurrentLine, long pos)
{
    pos = CurrentLine.FirstNonSpaceIndex (pos,-1,1);
    if (pos==-1) {
        ReadNextLine(f,&CurrentLine,&fState,false);
        pos = CurrentLine.FirstNonSpaceIndex (0,-1,1);
        if (pos==-1) {
            return false;
        }
    }
    CurrentLine.Trim (pos, -1);
    return true;
}


//_________________________________________________________

bool    SkipUntilNexusBlockEnd (FileState& fState, FILE* file, _String& CurrentLine, long pos)
{
    static _String endMark ("END");
    pos = CurrentLine.Find (endMark,pos+1,-1);
    while (pos<0) {
        ReadNextLine(file,&CurrentLine,&fState,false);
        if (!CurrentLine.sLength) {
            return false;
        }
        pos = CurrentLine.Find (endMark,0,-1);
        if (pos>=0) {
            pos = CurrentLine.Find (';',pos+endMark.sLength,-1);
            if (pos>=0) {
                CurrentLine.Trim (pos+endMark.sLength, -1);
                if (!CurrentLine.sLength) {
                    ReadNextLine(file,&CurrentLine,&fState,false);
                }
            } else {
                _String errMsg ("Found END w/o a trailing semicolon. Assuming end of block and skipping the rest of the line.");
                ReportWarning (errMsg);
                ReadNextLine(file,&CurrentLine,&fState,false);
            }
            return true;
        }
    }
    return false;
}
//_________________________________________________________
void    NexusParseEqualStatement (_String& source)
{
    long f = source.Find('=');
    if (f>=0) {
        f = source.FirstNonSpaceIndex (f+1,-1);
        if (f>=0) {
            source.Trim (f,-1);
            return;
        }
    }
    source = "";

}
//_________________________________________________________

bool ReadNextNexusStatement (FileState& fState, FILE* f, _String& CurrentLine, long pos, _String& blank, bool stopOnSpace, bool stopOnComma, bool stopOnQuote, bool NLonly, bool preserveSpaces, bool preserveQuotes)
{
    bool done          = false,
         insideLiteral = false,
         startedReading = false;

    long newPos = pos>0?pos+1:pos;
    char c;

    while (1) {
        while (newPos<CurrentLine.sLength) {
            c = CurrentLine.sData[newPos];
            if (isspace(c)) {
                if (stopOnSpace && startedReading && (!insideLiteral) && (!NLonly || (NLonly && (c==10 || c==13)))) {
                    done = true;
                    break;
                } else {
                    if (insideLiteral||preserveSpaces) {
                        blank<<' ';
                    }
                }
            } else {
                if (c==';' && ! insideLiteral) { // terminate always
                    done = true;
                    newPos++;
                    break;
                } else if (stopOnComma && c==',' && ! insideLiteral) {
                    done = true;
                    newPos++;
                    break;
                } else if (! preserveQuotes && (c=='\'' || c=='"') ) {
                    if (c=='\'') {
                        if (newPos+1<CurrentLine.sLength)
                            // check for a double quote
                        {
                            c = CurrentLine.sData[newPos+1];
                            if (c=='\'') {
                                newPos += 2;
                                blank<<c;
                                continue;
                            }
                            //else
                            //  if (!startedReading || insideLiteral)
                            //      newPos--;
                        }
                    }
                    if (startedReading &&stopOnQuote) {
                        done = true;
                        newPos++;
                        break;
                    } else {
                        insideLiteral = !insideLiteral;
                    }
                } else {
                    startedReading = true;
                    blank<<c;
                }
            }
            newPos++;
        }
        if (!done) {
            if (NLonly&&(blank.FirstNonSpaceIndex(0,-1,1)>=0)) {
                break;
            }
            ReadNextLine(f,&CurrentLine,&fState,false);
            newPos = 0;
            if (!CurrentLine.sLength) {
                c=';';
                break;
            }
        } else {
            break;
        }

    }
    blank.Finalize();
    blank.UpCase();
    if (newPos<CurrentLine.sLength) {
        CurrentLine.Trim (newPos,-1);
    } else {
        CurrentLine = "";
    }
    return c==';';
}

//_________________________________________________________

long    ReadNextNexusEquate (FileState& fState, FILE* f, _String& CurrentLine, long pos2, _String& blank, bool resetP, bool demandSemicolon)
{
    long pos = blank.Find ('=',pos2,-1), res;
    if (pos>=0) {
        if (pos<blank.sLength-1) {
            blank.Trim (pos+1,-1);
            return 1;
        } else {
            _String blank2 ((unsigned long)10, true);
            res = ReadNextNexusStatement (fState, f, CurrentLine, resetP?0:pos, blank2, true, true, false,false,false);
            if (blank2.sLength) {
                blank = blank2;
                return res?2:1;
            }
        }
        return 0;
    } else {
        _String blank2 ((unsigned long)10, true);
        res = ReadNextNexusStatement (fState, f, CurrentLine, pos2, blank2, true, true, false,false,false)?2:1;
        if (res!=2 && demandSemicolon) {
            if((res=ReadNextNexusEquate (fState, f, CurrentLine, 0, blank2))) {
                blank = blank2;
                return res;
            }
        } else if((res = ReadNextNexusEquate (fState, f, CurrentLine, 0, blank2, resetP, false))) {
            blank = blank2;
            return res;
        } else {
            return 0;
        }
    }
    return 0;
}

//_________________________________________________________
void    ProcessNexusTaxa (FileState& fState, long pos, FILE*f, _String& CurrentLine, _DataSet& result)
{
    _String key1 = "DIMENSIONS", key2 = "NTAX", key3 = "TAXLABELS", keyEnd = "END";

    bool    done = false;

    long    speciesExpected = -1, offset;

    while (!done) {
        if (!FindNextNexusToken (fState, f, CurrentLine, pos)) {
            break;
        }
        // now that we've got stuff to work with see what it is

        if (CurrentLine.beginswith (keyEnd, false)) {
            pos = -1;
            break;
        }

        if (CurrentLine.beginswith (key1, false)) {
            if (result.GetNames().lLength) { // check the number of dimensions
                // some data already present
                key1 = "Only one taxa definition per NEXUS file is recognized, the others will be ignored.";
                ReportWarning (key1);
                SkipUntilNexusBlockEnd (fState, f,CurrentLine, pos);
                break;
            } else {
                _String blank ((unsigned long)10, true);
                ReadNextNexusStatement (fState, f, CurrentLine, key1.sLength, blank, false,true, true,false,false);
                // this will actually return '= number'
                NexusParseEqualStatement (blank);
                speciesExpected = blank.toNum();
            }
        } else if (CurrentLine.beginswith (key3, false)) {
            if (speciesExpected == -1) {
                key1 = "TAXLABELS must be preceded by a valid NTAX statement. Skipping the entire TAXA block.";
                ReportWarning (key1);
                SkipUntilNexusBlockEnd (fState, f,CurrentLine, pos);
                break;
            } else {
                offset = key3.sLength;
                do {
                    _String blank ((unsigned long)10, true);
                    if (ReadNextNexusStatement (fState, f, CurrentLine,offset, blank, true,true,true,false,false)) {
                        if (blank.sLength) {
                            result.GetNames()&& & blank;
                        }
                        break;
                    } else {
                        if (blank.sLength) {
                            result.GetNames()&& & blank;
                        }
                    }
                    offset = 0;

                } while (1);
                if (result.GetNames().lLength!=speciesExpected) {
                    key1 = "TAXALABELS provided ";
                    key1 = key1& _String ((long)result.GetNames().lLength) &" species, whereas the NTAX statement promised:";
                    key1 = key1& _String (speciesExpected) & ". HYPHY will use TAXALABELS count.";
                    ReportWarning (key1);
                }
                done = true;
            }
        } else {
            long offSet;

            _String errMsg = CurrentLine.Cut (0, offSet = CurrentLine.FirstSpaceIndex(1,-1)) & " is not used by HYPHY";
            ReportWarning (errMsg);
            while (!done) {
                _String blank ((unsigned long)10, true);
                done = ReadNextNexusStatement (fState, f, CurrentLine, offSet, blank, false, false,true,false,false);
            }
            done = false;
        }

        if (!done) {
            pos = 0;
        }

        //if (!done)
        //{
        //  ReadNextLine(f,&CurrentLine,&fState,false);
        //}
    }

    SkipUntilNexusBlockEnd (fState, f,CurrentLine, pos);
}

//_________________________________________________________

void    ProcessNexusAssumptions (FileState& fState, long pos, FILE*f, _String& CurrentLine, _DataSet&)
{
    _String key1 = "CHARSET", keyEnd = "END",
            errMsg;

    bool    done = false;

    _List   charSetIDs,
            charSetSpec;

    while (!done) {
        if (!FindNextNexusToken (fState, f, CurrentLine, pos)) {
            break;
        }
        // now that we've got stuff to work with see what it is

        if (CurrentLine.beginswith (keyEnd, false)) {
            pos = -1;
            break;
        }

        if (CurrentLine.beginswith (key1, false)) { // actual tree strings & idents
            _String blank ((unsigned long)10, true);
            if (!ReadNextNexusStatement (fState, f, CurrentLine, key1.sLength, blank, false, false, false,false,true)) {
                errMsg = _String("CHARSET construct not followed by ';'.");
                ReportWarning (errMsg);
                done = true;
                break;
            } else {
                pos = blank.Find ('=',1,-1);
                if (pos==-1) {
                    errMsg = blank&": is not of the form Charset ID = specification of the partition.";
                    ReportWarning (errMsg);
                } else {
                    long pos2 = blank.FirstNonSpaceIndex (0,pos-1,-1);
                    if (pos2>=0) {
                        long j = blank.FirstNonSpaceIndex (0,pos2-1,1);
                        if (j>=0) {
                            if (blank.sData[j]=='*') {
                                j = blank.FirstNonSpaceIndex (j+1,pos2-1,1);
                            }

                            if (j>=0) {
                                _String CharSetID (blank.Cut (j,pos2));

                                if (!CharSetID.IsValidIdentifier()) {
                                    _String dummy;
                                    errMsg = CharSetID & " is not a valid data filter identifier in HYPHY. Replacing with: ";
                                    CharSetID.ConvertToAnIdent();
                                    FindUnusedObjectName (dummy,CharSetID, dataSetFilterNamesList, true);
                                    FindUnusedObjectName (dummy,CharSetID,charSetIDs, false);
                                    errMsg = errMsg &  CharSetID;
                                    ReportWarning (errMsg);
                                } else {
                                    _String dummy;
                                    FindUnusedObjectName (dummy,CharSetID, dataSetFilterNamesList, true);
                                    FindUnusedObjectName (dummy,CharSetID,charSetIDs, false);
                                }

                                //  now get the rest of the tree string
                                pos2 = blank.FirstNonSpaceIndex(pos+1,-1);
                                pos  = blank.FirstNonSpaceIndex(pos2,-1,-1);
                                blank.Trim (pos2,pos);
                                blank.CompressSpaces ();
                                blank = blank & " ";

                                _String hpSpec (blank.sLength+1, true),
                                        numberOne,
                                        numberTwo,
                                        numberThree;

                                bool    spoolInto2nd = false,
                                        spoolInto3rd = false,
                                        hitASpace    = false,
                                        okFlag         = true,
                                        firstFlag  = true;

                                for (long k=0; k<blank.sLength; k++) {
                                    char ch = blank.sData[k];

                                    if ((ch>='0' && ch<='9') || ch=='.') {
                                        if (spoolInto2nd) {
                                            numberTwo = numberTwo & ch;
                                        } else if (spoolInto3rd) {
                                            numberThree = numberThree & ch;
                                        } else {
                                            numberOne = numberOne & ch;
                                        }
                                        hitASpace = false;
                                    }

                                    if (ch==' ') {
                                        if (numberTwo.sLength == 1 && numberTwo.sData[0] == '.') {
                                            numberTwo = (long)fState.totalSitesRead;
                                        }

                                        if (spoolInto3rd) {
                                            spoolInto3rd = false;
                                            // handle 'every' n-th
                                            //if (!firstFlag)
                                            //hpSpec << ',';


                                            long    from = numberOne.toNum()-1,
                                                    upto = numberTwo.toNum()-1,
                                                    step = numberThree.toNum();

                                            if ((upto>=from)&&(step>0)) {
                                                if (!firstFlag) {
                                                    hpSpec << ',';
                                                }
                                                hpSpec << _String(from);
                                                for (long kk = from+step; kk<=upto; kk+=step) {
                                                    hpSpec << ',';
                                                    hpSpec << (_String)(kk);
                                                }

                                                numberOne   = empty;
                                                numberTwo   = empty;
                                                numberThree = empty;
                                            } else {
                                                errMsg = _String("Invalid from-to\\step specification: ") & blank.Cut (0,k) & " <=? " & blank.Cut (k+1,-1);
                                                ReportWarning (errMsg);
                                                okFlag = false;
                                                break;
                                            }

                                            firstFlag = false;
                                        } else {
                                            if (spoolInto2nd) {
                                                spoolInto2nd = false;
                                                if (!firstFlag) {
                                                    hpSpec << ',';
                                                }

                                                numberOne = numberOne.toNum()-1;
                                                hpSpec << numberOne;
                                                numberOne = ch;
                                                hpSpec << '-';
                                                numberTwo = numberTwo.toNum()-1;
                                                hpSpec << numberTwo;
                                                numberTwo = empty;
                                                firstFlag = false;

                                            } else {
                                                if (numberOne.sLength && numberOne.toNum() > 0) {
                                                    numberOne = numberOne.toNum()-1;
                                                    if (!firstFlag) {
                                                        hpSpec << ',';
                                                    }
                                                    hpSpec << numberOne;
                                                }
                                                numberOne = empty;
                                                firstFlag = false;
                                            }
                                        }
                                        //hitASpace = true;
                                        hitASpace = false;
                                    } else if (ch=='-') {
                                        if (spoolInto2nd||spoolInto3rd) {
                                            errMsg = _String("Misplaced '-' in CHARSET specification: ") & blank.Cut (0,k) & " <=? " & blank.Cut (k+1,-1);
                                            ReportWarning (errMsg);
                                            okFlag = false;
                                            break;
                                        }
                                        spoolInto2nd = true;
                                    } else if (ch=='\\') {
                                        if ((!spoolInto2nd)||spoolInto3rd) {
                                            errMsg = _String("Misplaced '\\' in CHARSET specification: ") & blank;
                                            ReportWarning (errMsg);
                                            okFlag = false;
                                            break;
                                        }
                                        spoolInto2nd = false;
                                        spoolInto3rd = true;
                                    }
                                }

                                hpSpec.Finalize();

                                if (okFlag) {
                                    charSetIDs && & CharSetID;
                                    charSetSpec && & hpSpec;
                                }
                            }
                        }
                        if (j<0) {
                            errMsg = _String("Could not find a charset identifier in:")&blank;
                            ReportWarning (errMsg);
                        }
                    } else {
                        errMsg = blank&" is not of the form CharSetID = char set string";
                        ReportWarning (errMsg);
                    }
                }
            }
        } else {
            long offSet;

            _String errMsg = CurrentLine.Cut (0, offSet = CurrentLine.FirstSpaceIndex(1,-1)) & " is not used by HYPHY";
            ReportWarning (errMsg);
            while (!done) {
                _String blank ((unsigned long)10, true);
                done = ReadNextNexusStatement (fState, f, CurrentLine, offSet, blank, false, false,true,false,false);
            }
            done = false;
        }

        if (!done) {
            //ReadNextLine(f,&CurrentLine,&fState,false);
            pos = 0;
        }
    }

    if (charSetIDs.lLength) {
        _String defineCharsets (128L, true);

        defineCharsets << dataFilePartitionMatrix;
        defineCharsets << "={2,";
        defineCharsets << _String ((long)charSetIDs.lLength);
        defineCharsets << "};\n";

        for (long id = 0; id < charSetIDs.lLength; id++) {
            defineCharsets << dataFilePartitionMatrix;
            defineCharsets << "[0][";
            defineCharsets << _String (id);
            defineCharsets << "]:=\"";
            defineCharsets << (_String*)charSetIDs(id);
            defineCharsets << "\";\n";
            defineCharsets << dataFilePartitionMatrix;
            defineCharsets << "[1][";
            defineCharsets << _String (id);
            defineCharsets << "]:=\"";
            defineCharsets << (_String*)charSetSpec(id);
            defineCharsets << "\";\n";
        }
        defineCharsets.Finalize();
        _ExecutionList defMx (defineCharsets);
        defMx.Execute();
        terminateExecution = false;
    }

    SkipUntilNexusBlockEnd (fState, f,CurrentLine, pos);
}

//_________________________________________________________

void    ProcessNexusTrees (FileState& fState, long pos, FILE*f, _String& CurrentLine, _DataSet& result)
{
    _String key1 = "TRANSLATE", key2 = "TREE", errMsg, keyEnd = "END";

    bool    done = false, readResult, good;
    _List   translationsFrom, translationsTo;
    _List   treeIdents, treeStrings;
    long    treeSelected = 0, insertPos = 0;

    while (!done) {
        if (!FindNextNexusToken (fState, f, CurrentLine, pos)) {
            break;
        }
        // now that we've got stuff to work with see what it is

        if (CurrentLine.beginswith (keyEnd, false)) {
            pos = -1;
            break;
        }

        if (CurrentLine.beginswith (key1, false)) {
            // set up translations between nodes and data labels
            long offset = key1.sLength;
            do {
                _String blank ((unsigned long)10, true);
                readResult = ReadNextNexusStatement (fState, f, CurrentLine, offset, blank, true, true,true,false,false);
                if (blank.sLength) {
                    if (translationsTo.lLength<translationsFrom.lLength) {
                        good = (result.GetNames().Find(&blank)>=0);
                        if (good) {
                            translationsTo.InsertElement (&blank, insertPos);
                        } else {
                            errMsg = blank & " is not a valid taxon name for TRANSLATE";
                            ReportWarning (errMsg);
                            translationsFrom.Delete (insertPos);
                        }

                    } else {
                        if (!readResult) {
                            //translationsFrom && & blank;
                            insertPos = translationsFrom.BinaryInsert (&blank);
                        }
                    }
                }
                if (readResult) {
                    break;
                }
                if  ((f&&feof(f))||(fState.theSource&&(fState.theSource->sLength<=fState.pInSrc))) {
                    break;
                }
                offset = 0;

            } while (1);
        } else if (CurrentLine.beginswith (key2, false)) { // actual tree strings & idents
            _String blank ((unsigned long)10, true);
            if (!ReadNextNexusStatement (fState, f, CurrentLine, key2.sLength, blank, false, false, false,false,false, true)) {
                errMsg = _String("TREE construct not followed by ';'.");
                ReportWarning (errMsg);
                done = true;
                break;
            } else {
                // here goes the tree string in the form: treeID = treeString
                // pull the ID out first - check if it is a valid one
                // next crudely parse the tree string, extracting species names and
                pos = blank.Find ('=',1,-1);
                if (pos==-1) {
                    errMsg = blank&": is not of the form Tree Name = Tree String";
                    ReportWarning (errMsg);
                } else {
                    long pos2 = blank.FirstNonSpaceIndex (0,pos-1,-1);
                    if (pos2>=0) {
                        long j = blank.FirstNonSpaceIndex (0,pos2-1,1);
                        if (j>=0) {
                            if (blank.sData[j]=='*') {
                                j = blank.FirstNonSpaceIndex (j+1,pos2-1,1);
                                treeSelected = treeIdents.lLength;
                            }
                            if (j>=0) {
                                _String TreeID (blank.Cut (j,pos2));

                                if (!TreeID.IsValidIdentifier()) {
                                    _String dummy;
                                    errMsg = TreeID & " is not a valid variable (tree) identifier in HYPHY. Replacing with: ";
                                    TreeID = "Tree";
                                    FindUnusedObjectName (dummy, TreeID, variableNames, true);
                                    FindUnusedObjectName (dummy, TreeID, treeIdents, true);
                                    errMsg = errMsg & TreeID;
                                    ReportWarning (errMsg);
                                } else {
                                    _String dummy;
                                    FindUnusedObjectName (dummy, TreeID, variableNames, true);
                                    FindUnusedObjectName (dummy, TreeID, treeIdents, true);
                                }

                                treeIdents && & TreeID;
                                //  now get the rest of the tree string
                                pos2 = blank.FirstNonSpaceIndex(pos+1,-1);
                                blank.Trim (pos2,-1);
                                treeStrings && & blank;
                            }
                        }
                        if (j<0) {
                            errMsg = _String("Could not find a tree identifier in:")&blank;
                            ReportWarning (errMsg);
                        }
                    } else {
                        errMsg = blank&" is not of the form TreeID = TreeString";
                        ReportWarning (errMsg);
                    }
                }

            }
        } else {

            long offSet = 0;

            _String errMsg = CurrentLine.Cut (0,CurrentLine.FirstSpaceIndex(1,-1)) & " is not used by HYPHY in TREES block";
            ReportWarning (errMsg);
            while (!done) {
                _String blank ((unsigned long)10, true);
                done = ReadNextNexusStatement (fState, f, CurrentLine, offSet, blank, false, false,true,false,false);
            }
            done = false;
        }

        if (!done) {
            //ReadNextLine(f,&CurrentLine,&fState,false);
            pos = 0;
        }
    }

    // now we shall check the string and match up node names with those present in the file

    for (long id = 0; id<treeStrings.lLength; id++) {
        key1 = *((_String*)treeStrings.lData[id]);
        long    treeLevel = 0, lastNode, i;
        _String revisedTreeString ((unsigned long)10, true);
        char    c;
        for (i=0; i<key1.sLength; i++) {
            switch (key1.sData[i]) {
            case '(': { // creating a new internal node one level down
                treeLevel++;
                revisedTreeString<<'(';
                break;
            }

            case ',':
            case ')': { // creating a new node on the same level and finishes updating the list of parameters
                if (key1.sData[i]==')') { // also create a new node on the same level
                    treeLevel--;
                }
                revisedTreeString<<key1.sData[i];
                break;
            }

            case ':' : { // tree branch definition
                lastNode = i+1;
                revisedTreeString<<':';
                c = key1.sData[lastNode];
                while (((c<='9')&&(c>='0'))||(c=='.')||(c=='-')||(c=='e')||(c=='E')) {
                    if (lastNode<key1.sLength) {
                        lastNode++;
                        revisedTreeString<<c;
                        c = key1.sData[lastNode];
                    } else {
                        break;
                    }
                }
                i = lastNode-1;
                break;
            }

            default: { // node name
                lastNode = i;
                c = key1.sData[lastNode];
                if (isspace (c)) {
                    break;
                }
                if (!(isalnum(c)||(c=='_'))) {
                    errMsg = ((_String("Node names should begin with a letter, a number, or an underscore")&key1.Cut(0,i)&"?"&key1.Cut(i+1,-1)));
                    ReportWarning (errMsg);
                    i = key1.sLength+2;
                    break;
                }
                while (isalnum(c)||(c=='_')) {
                    if (lastNode<key1.Length()) {
                        lastNode++;
                        c = key1.sData[lastNode];
                    } else {
                        break;
                    }
                }
                key2 = key1.Cut(i,lastNode-1);
                i = lastNode-1;
                lastNode = translationsFrom.BinaryFind (&key2);
                if (lastNode>=0) {
                    revisedTreeString<< (_String*)translationsTo.lData[lastNode];
                } else {
                    revisedTreeString<< &key2;
                }
                break;
            }
            }
        }
        revisedTreeString.Finalize();
        if (treeLevel) {
            errMsg = _String("Unbalanced '(,)' in the tree string:") & revisedTreeString;
            ReportWarning (errMsg);
        } else if (i==key1.sLength) {
            ((_String*)treeStrings.lData[id])->DuplicateErasing(&revisedTreeString);
        }
    }

    if (treeSelected < treeStrings.lLength) {
        setParameter (dataFileTree,1.0,fState.theNamespace);
        setParameter (dataFileTreeString, new _FString((*(_String*)treeStrings.lData[treeSelected])),fState.theNamespace);
    }

    if (treeStrings.lLength) {
        _String initTreeMatrix (1024L, true);

        initTreeMatrix << nexusFileTreeMatrix;
        initTreeMatrix << "={";
        initTreeMatrix << _String ((long)treeStrings.lLength);
        initTreeMatrix << ",2};\n";


        for (long id = 0; id < treeStrings.lLength; id++) {
            initTreeMatrix << nexusFileTreeMatrix;
            initTreeMatrix << '[';
            initTreeMatrix << _String (id);
            initTreeMatrix << "][0]=\"";
            initTreeMatrix << (_String*)treeIdents(id);
            initTreeMatrix << "\";\n";
            initTreeMatrix << nexusFileTreeMatrix;
            initTreeMatrix << '[';
            initTreeMatrix << _String (id);
            initTreeMatrix << "][1]=\"";
            initTreeMatrix << (_String*)treeStrings(id);
            initTreeMatrix << "\";\n";
        }
        initTreeMatrix.Finalize();

        _ExecutionList el (initTreeMatrix);
        el.Execute();
        terminateExecution = false;
    }
    SkipUntilNexusBlockEnd (fState, f,CurrentLine, pos);
}

//_________________________________________________________

void    ProcessNexusHYPHY (FileState& fState, long pos, FILE*file, _String& CurrentLine, _DataSet&)
{
    _String endMark ("END;"),
            bfBody  (128L,true);

    long      p2 = pos;
    pos = CurrentLine.FindAnyCase (endMark,pos+1,-1);

    fState.fileType = 0;

    if (pos>=0) {
        bfBody << CurrentLine.Cut (p2,pos-1);
        CurrentLine.Trim(pos+endMark.sLength,-1);
    } else {
        bfBody << CurrentLine.Cut (p2,-1);
        while (pos<0) {
            ReadNextLine(file,&CurrentLine,&fState,false,false);
            if (!CurrentLine.sLength) {
                break;
            }

            pos = CurrentLine.FindAnyCase (endMark,0,-1);
            if (pos>=0) {
                if (pos>0) {
                    bfBody << CurrentLine.Cut (0,pos-1);
                }

                CurrentLine.Trim (pos+endMark.sLength, -1);
                if (!CurrentLine.sLength) {
                    ReadNextLine(file,&CurrentLine,&fState,false,false);
                }

                break;
            } else {
                bfBody << CurrentLine;
            }

        }
    }
    bfBody.Finalize();

    nexusBFBody = bfBody;

    fState.fileType = 3;

    CurrentLine.UpCase();

}

//_________________________________________________________

bool    ProcessNexusData (FileState& fState, long pos, FILE*f, _String& CurrentLine, _DataSet& result)
{
    _String key1 ("DIMENSIONS"), key11 ("NTAX"), key12 ("NCHAR"),
            key2 ("FORMAT"),key21 ("DATATYPE"), key22 ("MISSING"), key23 ("GAP"), key24 ("SYMBOLS"),
            key25 ("EQUATE"), key26 ("MATCHCHAR"), key27 ("NOLABELS"), key28 ("INTERLEAVE"), key3 ("MATRIX"), keyEnd ("END"),
            errMsg, newAlph;

    bool    done = false,
            labels = true;

    char    charState = 0;

    _List   translations;
    char    missing = '?', gap = '-' , repeat = '.', charSwitcher;

    long    offSet, count, spExp = result.GetNames().lLength, sitesExp = 0;

    while (!done) {
        if (!FindNextNexusToken (fState, f, CurrentLine, pos)) {
            break;
        }

        if (CurrentLine.beginswith (keyEnd, false)) {
            pos = -1;
            break;
        }

        if (CurrentLine.beginswith (key1, false)) {
            offSet = key1.sLength;
            while (!done) {
                _String blank ((unsigned long)10, true);
                done = ReadNextNexusStatement (fState, f, CurrentLine, offSet, blank, true, true,true,false,false);

                if (blank.beginswith(key11, false)) {
                    if (result.GetNames().lLength) {
                        errMsg = "NTAX will override the definition of taxa names from the TAXA block";
                        ReportWarning (errMsg);
                    }
                    if (!(count=ReadNextNexusEquate (fState,f,CurrentLine, 0 ,blank))) {
                        errMsg = "NTAX is not followed by '= number-of-taxa'";
                        ReportWarning (errMsg);
                        done = true;
                    } else {
                        done = done||(count>1);
                        spExp = blank.toNum();
                        if(spExp<=0) {
                            errMsg = "NTAX must be a positive number";
                            ReportWarning (errMsg);
                            done = true;
                            spExp = result.GetNames().lLength?result.GetNames().lLength:1;
                        }
                    }
                } else if (blank.beginswith(key12, false)) {
                    if (!(count=ReadNextNexusEquate (fState,f,CurrentLine, 0 ,blank))) {
                        errMsg = "NCHAR is not followed by '= number-of-charaters'";
                        ReportWarning (errMsg);
                        done = true;
                    } else {
                        done = done||(count>1);
                        sitesExp = blank.toNum();
                    }
                }
                offSet = 0;
            }
            done = false;
        } else if (CurrentLine.beginswith (key2, false)) { // format instruction
            offSet = key2.sLength;
            while (!done) {
                charSwitcher = 0;
                _String blank ((unsigned long)10, true);
                done = ReadNextNexusStatement (fState, f, CurrentLine, offSet, blank, true, true,true,false,false);
                offSet = 0;
                blank.Trim (blank.FirstNonSpaceIndex(),-1);
                if (blank.beginswith(key21)) { // datatype
                    if (!(count=ReadNextNexusEquate (fState,f,CurrentLine, 0 ,blank))) {
                        errMsg = "DATATYPE is not followed by '= DNA|RNA|NUCLEOTIDE|PROTEIN|BINARY'";
                        ReportWarning (errMsg);
                        done = true;
                    } else {
                        done = done||(count>1);
                        if ((blank==_String("DNA"))||(blank==_String("RNA"))||(blank==_String("NUCLEOTIDE"))) {
                            if (newAlph.sLength) {
                                errMsg = _String("DNA|RNA|NUCLEOTIDE datatype directive will over-ride the custom symbols definition: ") & newAlph;
                                newAlph = empty;
                                ReportWarning (errMsg);
                            }
                            if (done) {
                                done = false;
                                break;
                            }
                            continue;
                        } else if (blank==_String("PROTEIN") || blank == _String ("BINARY")) {
                            charState = 1+(blank==_String("BINARY"));
                            if (newAlph.sLength) {
                                errMsg = _String("PROTEIN|BINARY datatype directive will override the custom symbols definition: ") & newAlph;
                                newAlph = empty;
                                ReportWarning (errMsg);
                            }
                            if (done) {
                                done = false;
                                break;
                            }
                            continue;
                        } else {
                            errMsg = blank&" is not a recognized data type (DNA|RNA|NUCLEOTIDE|PROTEIN|BINARY are allowed).";
                            ReportWarning (errMsg);
                            done = false;
                        }
                    }
                } else if (blank.beginswith(key22, false)) { // MISSING
                    charSwitcher = 1;
                } else if (blank.beginswith(key23, false)) { // GAP
                    charSwitcher = 2;
                } else if (blank.beginswith(key26, false)) { // MATCHCHAR
                    charSwitcher = 3;
                } else if (blank.beginswith(key27, false)) { // NOLABELS
                    labels = false;
                } else if (blank.beginswith(key28, false)) { // INTERLEAVE
                    fState.interleaved = true;
                } else if (blank.beginswith(key24, false)) { // SYMBOLS
                    count=ReadNextNexusEquate (fState,f,CurrentLine, 0 ,blank, true,false);
                    if (blank.sLength == 0) {
                        errMsg = blank& _String("is not of the form SYMBOLS = \"sym1 sym2 ...\". The entire block is ignored.");
                        ReportWarning (errMsg);
                        done = true;
                        break;
                    }
                    _String tempNewAlpha ((unsigned long)10, true);
                    //bool  mult = false;
                    for (long pos1 = 0; pos1<blank.sLength; pos1++) {
                        charSwitcher = blank.sData[pos1];
                        if (!isspace(charSwitcher))
                            //{
                            //if (!mult)
                            //{
                        {
                            tempNewAlpha<<charSwitcher;
                        }
                        //mult = true;
                        //}
                        //else
                        //{
                        //done = true;
                        //errMsg = "Multicharacter symbols are not supported by HYPHY kernel. Skipping the entire data block";
                        //break;
                        //}
                        //}
                        //else
                        //mult = false;
                    }
                    if (done) {
                        break;
                    }
                    tempNewAlpha.Finalize();
                    newAlph = tempNewAlpha;
                    charSwitcher = 0;
                    done = done||(count>1);
                } else if (blank.beginswith(key25, false)) { // EQUATE
                    blank.Trim(key25.sLength,-1);
                    if (!(count=ReadNextNexusEquate (fState,f,CurrentLine, 0,blank,true,false))) {
                        errMsg = errMsg&" is not followed by '=char'";
                        ReportWarning (errMsg);
                        done = true;
                    }
                    done = done||(count>1);
                    // blank now contains a full list of the form token=(token)
                    _String symbol, meaning;
                    bool    symbolDefined = false, meaningDefined = false;
                    for (count=0; count<blank.sLength; count++) {
                        charSwitcher = blank.sData[count];
                        if (isspace(charSwitcher)) {
                            /*if (meaningDefined)
                            {
                                translations&& &symbol;
                                translations&& &meaning;
                                symbol = "";
                                meaning = "";
                                symbolDefined = false;
                                meaningDefined = false;
                            }   */
                            continue;
                        } else if (charSwitcher == '=') {
                            if (symbolDefined&&!meaningDefined) {
                                meaningDefined = true;
                            }
                        } else
                            //if (((charSwitcher>='A')&&(charSwitcher<='Z'))||((charSwitcher>='0')&&(charSwitcher<='9')))
                        {
                            if (!symbolDefined) {
                                symbolDefined = true;
                                symbol = charSwitcher;
                                continue;
                            }
                            if (!meaningDefined) {
                                errMsg = "EQUATE can only be used to define single-character tokens. Ignoring the EQUATE command.";
                                translations.Clear();
                                break;
                            }
                            meaning = meaning & charSwitcher;
                        }
                    }
                    if (symbol.sLength&&meaning.sLength) {
                        translations&& &symbol;
                        translations&& &meaning;
                    }
                    charSwitcher = 0;
                    blank = empty;
                }

                offSet = 0;
                if (charSwitcher) {
                    switch (charSwitcher) {
                    case 1:
                        errMsg = "MISSING";
                        break;
                    case 2:
                        errMsg = "GAP";
                        break;
                    case 3:
                        errMsg = "MATCHCHAR";
                        break;
                    }
                    if (!(count=ReadNextNexusEquate (fState,f,CurrentLine, 0 ,blank, true))) {
                        errMsg = errMsg&" is not followed by '=char'";
                        ReportWarning (errMsg);
                        done = true;
                    } else {
                        done = done||(count>1);
                        if (blank.sLength!=1) {
                            errMsg = blank&" is not a valid "&errMsg&" character.";
                            ReportWarning (errMsg);
                        }
                    }
                    switch (charSwitcher) {
                    case 1:
                        missing = blank.getChar(0);
                        if (gap == missing) {
                            gap = 0;
                        }
                        if (repeat == missing) {
                            repeat = 0;
                        }

                        break;
                    case 2:
                        gap = blank.getChar(0);
                        if (missing == gap) {
                            missing = 0;
                        }
                        if (repeat == gap) {
                            repeat= 0;
                        }

                        break;
                    case 3:
                        repeat = blank.getChar(0);
                        if (missing == repeat) {
                            missing = 0;
                        }
                        if (repeat == gap) {
                            gap = 0;
                        }

                        break;
                    }
                }

                if (done) {
                    done = false;
                    break;
                }
                done = false;
            }
        } else if (CurrentLine.beginswith (key3, false)) { // matrix instruction
            // if needed, set up a new symbol set
            offSet = key3.sLength;
            if (newAlph.sLength>1) { // a valid new alphabet set
                checkTTStatus (&fState);
                fState.translationTable->AddBaseSet (newAlph);
            } else {
                if (charState) {
                    checkTTStatus (&fState);
                    if (charState==1) {
                        newAlph = aminoAcidOneCharCodes;
                        fState.translationTable->baseLength = 20;
                    } else {
                        newAlph = binaryOneCharCodes;
                        fState.translationTable->baseLength = 2;
                    }
                } else {
                    newAlph = dnaOneCharCodes;
                }
            }
            // set up translations
            if (translations.lLength) {
                checkTTStatus (&fState);
            }

            for (long k = 0; k<translations.lLength; k+=2) {
                char c = ((_String*)translations(k))->sData[0];
                fState.translationTable->AddTokenCode (c,*((_String*)translations(k+1)));
            }

            if (fState.translationTable->GetSkipChar()!=missing) {
                checkTTStatus (&fState);
                fState.translationTable->AddTokenCode (missing,newAlph);
            }

            if (fState.translationTable->GetGapChar()!=gap) {
                checkTTStatus (&fState);
                newAlph = "";
                fState.translationTable->AddTokenCode (gap,newAlph);
            }

            if (repeat == missing) {
                repeat = 0;
            }

            fState.repeat               = repeat;
            fState.skip                 = missing;

            //fState.totalSitesExpected   = sitesExp;

            // now proceed to read the data

            long loopIterations = 0;
            if (labels == true) {
                result.GetNames().Clear();
            }


            done = false;
            while (1) {
                _String  blank  ((unsigned long)10, true),
                         blank2 ((unsigned long)10, true),
                         *source;

                done = ReadNextNexusStatement (fState, f, CurrentLine, offSet?offSet+1:0, blank, true, true,true,!labels,false);
                offSet = 0;
                // in each line that should produce first the name of the taxon
                // and then the data string for the taxon

                if (labels) {
                    if (result.GetNames().lLength<spExp) {
                        if ((spExp>0)&&(blank.sLength==0)) {
                            errMsg = _String("Could not find NTAX taxon names in the matrix. Read: ")&_String((long)result.GetNames().lLength) & " sequences.";
                            ReportWarning (errMsg);
                            done = true;
                            blank2.Finalize(); // dmalloc fix 06162005
                            break;
                        }

                        if (!(sitesExp&&fState.curSite&&(fState.curSite<sitesExp)&&(!fState.interleaved))) {
                            result.AddName(blank);
                            fState.totalSpeciesExpected++;
                        }
                    } else {
                        if (done) {
                            blank2.Finalize();  // dmalloc fix 06162005
                            break;
                        }
                    }

                    if (!(sitesExp&&fState.curSite&&(fState.curSite<sitesExp)&&(!fState.interleaved))) {
                        done = ReadNextNexusStatement (fState, f, CurrentLine, offSet, blank2, true, true,true,true,false);
                        source = &blank2;
                    } else {
                        blank2.Finalize();  // dmalloc fix 06162005
                        source = &blank;
                    }
                } else {
                    blank2.Finalize();  // dmalloc fix 06162005

                    if (loopIterations<spExp) {
                        if (!(sitesExp&&fState.curSite&&(fState.curSite<sitesExp)&&(!fState.interleaved))) {
                            fState.totalSpeciesExpected++;
                        } else {
                            loopIterations --;
                        }
                    }
                    source = &blank;
                }

                if (source->sLength==0) {
                    errMsg = _String("Could not find NTAX data strings in the matrix. Read: ")&_String((long)result.GetNames().lLength) & " sequences.";
                    ReportWarning (errMsg);
                    done = true;
                    break;
                }
                loopIterations++;
                ISelector (fState, *source, result);

                if (done)
                    if (loopIterations>=fState.totalSpeciesExpected) {
                        break;    // finished reading
                    }

                if  ((f&&feof(f))||(fState.theSource&&(fState.theSource->sLength<=fState.pInSrc))) {
                    break;
                }
            }


            if (result.GetNames().lLength!=spExp) {
                errMsg = _String ("Expected ")&spExp&" taxa, but found "&(long)result.GetNames().lLength;
                ReportWarning(errMsg);
            }
            if (result.lLength!=sitesExp && result.InternalStorageMode() == 0) {
                errMsg = _String ("Expected ")&sitesExp&" sites, but found "&(long)result.lLength;
                ReportWarning(errMsg);
            }
            if (loopIterations%spExp) {
                errMsg = _String ("There is an inconsistency between NTAX and the number of data strings in the matrix");
                ReportWarning(errMsg);
            }
            done = true;
        } else {
            errMsg = CurrentLine.Cut (0, CurrentLine.FirstSpaceIndex(1,-1)) & " is not used by HYPHY";
            ReportWarning (errMsg);
            while (!done) {
                _String blank ((unsigned long)10, true);
                done = ReadNextNexusStatement (fState, f, CurrentLine, offSet, blank, true, false,true,false,false);
            }
            done = false;
        }
        if (!done) {
            if (CurrentLine.sLength == 0) {
                ReadNextLine(f,&CurrentLine,&fState,false);
            }
            pos = 0;
            if (CurrentLine.sLength==0) {
                done = true;
            }
        }
    }

    SkipUntilNexusBlockEnd (fState, f,CurrentLine, pos);
    return true;
}

//_________________________________________________________

void    ReadNexusFile (FileState& fState, FILE*file, _DataSet& result)
{
    bool   dataRead = false, lookForEnd = false;
    long   f,g;

    fState.fileType = 3; // NEXUS
    _String CurrentLine, beginMark ("BEGIN"), endMark ("END"), blockName, data ("DATA"), chars ("CHARACTERS"),
            taxa ("TAXA"), trees ("TREES"), assumptions ("ASSUMPTIONS"), hyphy ("HYPHY"), sets ("SETS");

    ReadNextLine(file,&CurrentLine,&fState,false);
    while (CurrentLine.sLength) {
        f = 0;
        while ((f = CurrentLine.FindAnyCase(beginMark,f,-1 ))>=0) {
            f = CurrentLine.FirstNonSpaceIndex (f+beginMark.sLength,-1,1);
            if (f!=-1) { // process
                g = CurrentLine.Find (';', f, -1);
                if (g!=-1) {
                    blockName = CurrentLine.Cut (f,g-1);
                    // dispatch to block readers
                    if (blockName.iEqual(&data)) {
                        blockName = blockName &" block is now deprecated in NEXUS and should not be used.";
                        ReportWarning (blockName);

                        if (!dataRead) {
                            dataRead = ProcessNexusData (fState, g+1, file, CurrentLine, result);
                        }
                        //SkipUntilNexusBlockEnd (fState,file,CurrentLine,f);

                        else {
                            blockName = "Only one data set per NEXUS file is read by ReadDataSet - the 1st valid one.";
                            ReportWarning (blockName);
                        }
                    } else if (blockName.iEqual(&taxa)) {
                        if (!dataRead) {
                            ProcessNexusTaxa (fState, g+1, file, CurrentLine, result);
                        } else {
                            blockName = "The TAXA block was encountered after CHARACTER had been read and will be ignored.";
                            ReportWarning (blockName);
                        }
                    } else if (blockName.iEqual(&trees)) {
                        ProcessNexusTrees (fState, g+1, file, CurrentLine, result);
                    } else if (blockName.iEqual(&chars)) {
                        if (!dataRead) {
                            dataRead = ProcessNexusData (fState, g+1, file, CurrentLine, result);
                        } else {
                            blockName = "Only one data set per NEXUS file is read by ReadDataSet - the 1st valid one.";
                            ReportWarning (blockName);
                        }
                    } else if (blockName.iEqual(&assumptions)||blockName.iEqual(&sets)) {
                        ProcessNexusAssumptions (fState, g+1, file, CurrentLine, result);
                    } else if (blockName.iEqual(&hyphy)) {
                        ProcessNexusHYPHY (fState, g+1, file, CurrentLine, result);
                    } else {
                        blockName = _String("NEXUS blocks ")&blockName&(" are not used by HYPHY.");
                        ReportWarning (blockName);
                        lookForEnd = true;
                        break;
                        // now look for the end of this block
                    }

                } else {
                    break;
                }
            } else {
                blockName = _String ("NEXUS BEGIN must be followed by the name of the block. Skipping until next BEGIN statement.");
                ReportWarning (blockName);
                break;
            }
        }

        if (lookForEnd) {
            lookForEnd = false;
            SkipUntilNexusBlockEnd (fState,file,CurrentLine,f);
        } else {
            ReadNextLine(file,&CurrentLine,&fState,false);
        }

    }

}
