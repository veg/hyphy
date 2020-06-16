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

#include <math.h>
#include <string.h>
#include <ctype.h>

#include "hy_string_buffer.h"
#include "dataset.h"
#include "list.h"
#include "batchlan.h"
#include "hbl_env.h"
#include "global_object_lists.h"
#include "global_things.h"


using namespace hy_global;
using namespace hyphy_global_objects;


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
bool    ProcessNexusData            (FileState&, long,  FILE*, _String&, _DataSet&);
void    ProcessNexusHYPHY           (FileState&, long,  FILE*, _String&, _DataSet&);
void    ProcessNexusAssumptions     (FileState&, long,  FILE*, _String&, _DataSet&);
void    ProcessNexusTaxa            (FileState&,long, FILE*, _String&, _DataSet&);
void    ProcessNexusTrees           (FileState&, long, FILE*, _String&, _DataSet&);
bool    FindNextNexusToken          (FileState& fState, FILE* f, _String& CurrentLine, long pos);
bool    SkipUntilNexusBlockEnd      (FileState& fState, FILE* f, _String& CurrentLine, long pos);
bool    ReadNextNexusStatement      (FileState&, FILE* , _String&, long, _StringBuffer&, bool, bool = true, bool = true, bool = false, bool = false, bool = false);
long    ReadNextNexusEquate         (FileState&, FILE* , _String&, long, _String&, bool = false, bool = true);
void    NexusParseEqualStatement    (_String&);

static auto  error_conext = [] (_String const& buffer, long position) -> const _String {return (buffer.Cut (0,position) & " <=? " & buffer.Cut (position+1,kStringEnd)).Enquote();};


//_________________________________________________________

bool    FindNextNexusToken (FileState& fState, FILE* f, _String& CurrentLine, long pos) {
    pos = CurrentLine.FirstNonSpaceIndex (pos,-1,kStringDirectionForward);
    if (pos==kNotFound) {
        ReadNextLine(f,&CurrentLine,&fState,false);
        pos = CurrentLine.FirstNonSpaceIndex (0,-1,kStringDirectionForward);
        if (pos==kNotFound) {
            return false;
        }
    }
    CurrentLine.Trim (pos, kStringEnd);
    return true;
}


//_________________________________________________________

bool    SkipUntilNexusBlockEnd (FileState& fState, FILE* file, _String& CurrentLine, long pos) {
    static const _String endMark ("END");
    pos = CurrentLine.Find (endMark,pos+1,kStringEnd);
    while (pos == kNotFound) {
        ReadNextLine(file,&CurrentLine,&fState,false);
        if (CurrentLine.empty()) {
            return false;
        }
        pos = CurrentLine.Find (endMark,0,kStringEnd);
        if (pos != kNotFound) {
            pos = CurrentLine.Find (';',pos+endMark.length(),kStringEnd);
            if (pos != kNotFound) {
                CurrentLine.Trim (pos+endMark.length(), kStringEnd);
                if (CurrentLine.empty()) {
                    ReadNextLine(file,&CurrentLine,&fState,false);
                }
            } else {
                ReportWarning ("Found END w/o a trailing semicolon. Assuming end of block and skipping the rest of the line.");
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
    if (f != kNotFound) {
        f = source.FirstNonSpaceIndex (f+1,kStringEnd);
        if (f != kNotFound) {
            source.Trim (f,kStringEnd);
            return;
        }
    }
    source.Clear();

}
//_________________________________________________________

bool ReadNextNexusStatement (FileState& fState, FILE* f, _String& CurrentLine, long pos, _StringBuffer & blank, bool stopOnSpace, bool stopOnComma, bool stopOnQuote, bool NLonly, bool preserveSpaces, bool preserveQuotes) {
    bool done          = false,
         insideLiteral = false,
         startedReading = false;

    long newPos = pos>0?pos+1L:pos;
    char c = '\0';

    while (1) {
        while (newPos<CurrentLine.length()) {
            c = CurrentLine.char_at (newPos);
            if (isspace(c)) {
                if (stopOnSpace && startedReading && (!insideLiteral) && (!NLonly || (NLonly && (c=='\r' || c=='\n')))) {
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
                        if (newPos+1<CurrentLine.length())
                            // check for a double quote
                        {
                            c = CurrentLine.char_at (newPos+1);
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
            if (NLonly&&(blank.FirstNonSpaceIndex(0,kStringEnd,kStringDirectionForward)>=0)) {
                break;
            }
            ReadNextLine(f,&CurrentLine,&fState,false);
            newPos = 0;
            if (CurrentLine.empty()) {
                c=';';
                break;
            }
        } else {
            break;
        }

    }
    // TODO 20170821: SLKP, this needs to be case sensitive
    blank.ChangeCaseInPlace(kStringUpperCase);
    if (newPos<CurrentLine.length()) {
        CurrentLine.Trim (newPos,kStringEnd);
    } else {
        CurrentLine.Clear();
    }
    return c==';';
}

//_________________________________________________________

long    ReadNextNexusEquate (FileState& fState, FILE* f, _String& CurrentLine, long pos2, _String& blank, bool resetP, bool demandSemicolon) {
    long pos = blank.Find ('=',pos2,-1), res;
    if (pos>=0) {
        if (pos<blank.length()-1) {
            blank.Trim (pos+1,-1);
            return 1;
        } else {
            _StringBuffer buffer (128UL);
            res = ReadNextNexusStatement (fState, f, CurrentLine, resetP?0:pos, buffer, true, true, false,false,false);
            if (!buffer.empty()) {
                blank = buffer;
                return res?2:1;
            }
        }
        return 0;
    } else {
        _StringBuffer buffer (128UL);
        res = ReadNextNexusStatement (fState, f, CurrentLine, pos2, buffer, true, true, false,false,false)?2:1;
        if (res!=2 && demandSemicolon) {
            if((res=ReadNextNexusEquate (fState, f, CurrentLine, 0, buffer))) {
                blank = buffer;
                return res;
            }
        } else if((res = ReadNextNexusEquate (fState, f, CurrentLine, 0, buffer, resetP, false))) {
            blank = buffer;
            return res;
        } else {
            return 0;
        }
    }
    return 0;
}

//_________________________________________________________
void    ProcessNexusTaxa (FileState& fState, long pos, FILE*f, _String& CurrentLine, _DataSet& result) {
    static const _String key1 = "DIMENSIONS", key2 = "NTAX", key3 = "TAXLABELS", keyEnd = "END";

    bool    done = false;

    long    speciesExpected = -1, offset;

    while (!done) {
        if (!FindNextNexusToken (fState, f, CurrentLine, pos)) {
            break;
        }
        // now that we've got stuff to work with see what it is

        if (CurrentLine.BeginsWith (keyEnd, false)) {
            pos = -1;
            break;
        }

        if (CurrentLine.BeginsWith (key1, false)) {
            if (result.GetNames().lLength) { // check the number of dimensions
                // some data already present
                ReportWarning ("Only one taxa definition per NEXUS file is recognized, the others will be ignored.");
                SkipUntilNexusBlockEnd (fState, f,CurrentLine, pos);
                break;
            } else {
                _StringBuffer buffer (128UL);
                ReadNextNexusStatement (fState, f, CurrentLine, key1.length(), buffer, false,true, true,false,false);
                // this will actually return '= number'
                NexusParseEqualStatement (buffer);
                speciesExpected = buffer.to_long();
            }
        } else if (CurrentLine.BeginsWith (key3, false)) {
            if (speciesExpected == -1) {
                ReportWarning ("TAXLABELS must be preceded by a valid NTAX statement. Skipping the entire TAXA block.");
                SkipUntilNexusBlockEnd (fState, f,CurrentLine, pos);
                break;
            } else {
                offset = key3.length();
                do {
                    _StringBuffer buffer (128UL);
                    if (ReadNextNexusStatement (fState, f, CurrentLine,offset, buffer, true,true,true,false,false)) {
                        if (buffer.nonempty()) {
                          result.AddName(buffer);
                        }
                        break;
                    } else {
                        if (buffer.nonempty()) {
                          result.AddName(buffer);;
                        }
                    }
                    offset = 0;

                } while (1);
                if (result.GetNames().lLength!=speciesExpected) {
                    ReportWarning ( _String ("TAXALABELS provided ") &
                                    _String ((long)result.GetNames().lLength) &" species, whereas the NTAX statement promised:" &
                                    _String (speciesExpected) & ". HYPHY will use TAXALABELS count.");
                }
                done = true;
            }
        } else {
            long offSet = 0;

            ReportWarning (CurrentLine.Cut (0, CurrentLine.FirstSpaceIndex(1,kStringEnd)) & " is not used by HYPHY");
            while (!done) {
                _StringBuffer buffer (128UL);
                done = ReadNextNexusStatement (fState, f, CurrentLine, offSet, buffer, false, false,true,false,false);
            }
            done = false;
        }

        if (!done) {
            pos = 0;
        }
    }

    SkipUntilNexusBlockEnd (fState, f,CurrentLine, pos);
}

//_________________________________________________________

void    ProcessNexusAssumptions (FileState& fState, long pos, FILE*f, _String& CurrentLine, _DataSet&) {
    static const _String key1 = "CHARSET", keyEnd = "END";
 
    bool    done = false;

    _List   charSetIDs,
            charSetSpec;

    while (!done) {
        if (!FindNextNexusToken (fState, f, CurrentLine, pos)) {
            break;
        }
        // now that we've got stuff to work with see what it is

        if (CurrentLine.BeginsWith (keyEnd, false)) {
            pos = -1;
            break;
        }

        if (CurrentLine.BeginsWith (key1, false)) { // actual tree strings & idents
            _StringBuffer buffer (128UL);
            if (!ReadNextNexusStatement (fState, f, CurrentLine, key1.length(), buffer, false, false, false,false,true)) {
                ReportWarning ("CHARSET construct not followed by ';'.");
                break;
            } else {
                pos = buffer.Find ('=',1,kStringEnd);
                if (pos==-1) {
                    ReportWarning (buffer.Enquote() & " is not of the form Charset ID = specification of the partition.");
                } else {
                    long pos2 = buffer.FirstNonSpaceIndex (0,pos-1,kStringDirectionBackward);
                    if (pos2 != kNotFound) {
                        long j = buffer.FirstNonSpaceIndex (0,pos2-1,kStringDirectionForward);
                        if (j != kNotFound) {
                            if (buffer.char_at (j) == '*') {
                                j = buffer.FirstNonSpaceIndex (j+1,pos2-1,kStringDirectionForward);
                            }

                            if (j != kNotFound) {
                                _String nexus_name (buffer,j,pos2),
                                        charset_id (nexus_name);

                                if (!nexus_name.IsValidIdentifier(fIDAllowCompound)) {
                                    charset_id = nexus_name.ConvertToAnIdent();
                                }
                              
                                charset_id = charSetIDs.GenerateUniqueNameForList(GenerateUniqueObjectIDByType (nexus_name, HY_BL_DATASET_FILTER), false);
                              
                                if (charset_id != nexus_name) {
                                  ReportWarning(nexus_name.Enquote('\'') & " has been renamed to " & charset_id.Enquote('\'') & " to avoid naming conflicts and/or comply with HyPhy ID requirements");
                                }
                                

                                //  now get the rest of the tree string
                                pos2 = buffer.FirstNonSpaceIndex(pos+1,kStringEnd);
                                pos  = buffer.FirstNonSpaceIndex(pos2,kStringEnd,kStringDirectionBackward);
                                buffer.Trim (pos2,pos);
                                buffer  = buffer.CompressSpaces () & " ";

                                _StringBuffer hpSpec (buffer.length()+1UL);
                              
                                _String numberOne,
                                        numberTwo,
                                        numberThree;

                                bool    spoolInto2nd = false,
                                        spoolInto3rd = false,
                                        okFlag         = true,
                                        firstFlag  = true;

                                for (long k=0; k<buffer.length(); k++) {
                                    char ch = buffer.char_at(k);

                                    if ((ch>='0' && ch<='9') || ch=='.') {
                                        if (spoolInto2nd) {
                                            numberTwo = numberTwo & ch;
                                        } else if (spoolInto3rd) {
                                            numberThree = numberThree & ch;
                                        } else {
                                            numberOne = numberOne & ch;
                                        }
                                    }

                                    if (ch==' ') {
                                        if (numberTwo.length() == 1 && numberTwo.char_at (0) == '.') {
                                            numberTwo = (long)fState.totalSitesRead;
                                        }

                                        if (spoolInto3rd) {
                                            spoolInto3rd = false;
                                            // handle 'every' n-th
 

                                            long    from = numberOne.to_long()-1,
                                                    upto = numberTwo.to_long()-1,
                                                    step = numberThree.to_long();

                                            if ((upto>=from)&&(step>0)) {
                                                if (!firstFlag) {
                                                    hpSpec << ',';
                                                }
                                                hpSpec << _String(from);
                                                for (long kk = from+step; kk<=upto; kk+=step) {
                                                    hpSpec << ',' << (_String)(kk);
                                                }

                                                numberOne.Clear();
                                                numberTwo.Clear();
                                                numberThree.Clear();
                                              
                                            } else {
                                                 ReportWarning (_String("Invalid from-to\\step specification: ") & error_conext (buffer, k));
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

                                                numberOne = numberOne.to_long ()-1;
                                                hpSpec << numberOne;
                                                numberOne = ch;
                                                hpSpec << '-';
                                                numberTwo = numberTwo.to_long()-1;
                                                hpSpec << numberTwo;
                                                numberTwo.Clear();
                                                firstFlag = false;

                                            } else {
                                              long n1;
                                                if (numberOne.nonempty() && (n1 = numberOne.to_long() > 0)) {
                                                    numberOne = n1-1;
                                                    if (!firstFlag) {
                                                        hpSpec << ',';
                                                    }
                                                    hpSpec << numberOne;
                                                }
                                                numberOne.Clear();
                                                firstFlag = false;
                                            }
                                        }
                                        //hitASpace = true;
                                      
                                    } else if (ch=='-') {
                                        if (spoolInto2nd||spoolInto3rd) {
                                            ReportWarning (_String("Misplaced '-' in CHARSET specification: ") & error_conext (buffer, k));
                                            okFlag = false;
                                            break;
                                        }
                                        spoolInto2nd = true;
                                    } else if (ch=='\\') {
                                        if ((!spoolInto2nd)||spoolInto3rd) {
                                            ReportWarning (_String("Misplaced '\\' in CHARSET specification: ") & buffer.Enquote());
                                            okFlag = false;
                                            break;
                                        }
                                        spoolInto2nd = false;
                                        spoolInto3rd = true;
                                    }
                                }


                                if (okFlag) {
                                    charSetIDs  && & charset_id;
                                    charSetSpec && & hpSpec;
                                }
                            }
                        }
                        if (j<0) {
                            ReportWarning (_String("Could not find a charset identifier in: ")& buffer.Enquote());
                        }
                    } else {
                        ReportWarning (buffer.Enquote() &" is not of the form CharSetID = char set string");
                    }
                }
            }
        } else {
            long offSet = 0L;

            ReportWarning (CurrentLine.Cut (0, CurrentLine.FirstSpaceIndex(1,-1)) & " is not used by HYPHY");
            while (!done) {
                _StringBuffer buffer (128UL);
                done = ReadNextNexusStatement (fState, f, CurrentLine, offSet, buffer, false, false,true,false,false);
            }
            done = false;
        }

        if (!done) {
            //ReadNextLine(f,&CurrentLine,&fState,false);
            pos = 0;
        }
    }

    if (charSetIDs.lLength) {
        _StringBuffer defineCharsets (256UL);

        defineCharsets << hy_env::data_file_partition_matrix << "={2," << _String ((long)charSetIDs.lLength) << "};\n";

        for (long id = 0; id < charSetIDs.lLength; id++) {
            defineCharsets << hy_env::data_file_partition_matrix
             << "[0]["
             << _String (id)
             << "]:=\""
             << (_String*)charSetIDs(id)
             << "\";\n"
             << hy_env::data_file_partition_matrix
             << "[1]["
             << _String (id)
             << "]:=\""
             << (_String*)charSetSpec(id)
             << "\";\n";
        }
         _ExecutionList defMx (defineCharsets);
        defMx.Execute();
        terminate_execution = false;
    }

    SkipUntilNexusBlockEnd (fState, f,CurrentLine, pos);
}

//_________________________________________________________

void    ProcessNexusTrees (FileState& fState, long pos, FILE*f, _String& CurrentLine, _DataSet& result) {
    static _String const key1 = "TRANSLATE", key2 = "TREE", errMsg, keyEnd = "END";

    bool    done = false, readResult, good;
    _List   translationsFrom, translationsTo;
    _List   treeIdents, treeStrings;
    long    treeSelected = 0, insertPos = 0;

    while (!done) {
        
        if (!FindNextNexusToken (fState, f, CurrentLine, pos)) {
            break;
        }
        // now that we've got stuff to work with see what it is

        if (CurrentLine.BeginsWith (keyEnd, false)) {
            pos = -1;
            break;
        }

        if (CurrentLine.BeginsWith (key1, false)) {
            // set up translations between nodes and data labels
            long offset = key1.length();
            do {
                _StringBuffer buffer (128UL);
                readResult = ReadNextNexusStatement (fState, f, CurrentLine, offset, buffer, true, true,true,false,false);
                if (buffer.nonempty()) {
                    if (translationsTo.lLength<translationsFrom.lLength) {
                        good = (result.GetNames().FindObject(&buffer)>=0);
                        if (good) {
                            translationsTo.InsertElement (&buffer, insertPos);
                        } else {
                            ReportWarning (buffer.Enquote() & " is not a valid taxon name for TRANSLATE" );
                            translationsFrom.Delete (insertPos);
                        }

                    } else {
                        if (!readResult) {
                            insertPos = translationsFrom.BinaryInsert (&buffer);
                        }
                    }
                }
                if (readResult) {
                    break;
                }
                if  ((f&&feof(f))||(fState.theSource&&(fState.theSource->length()<=fState.pInSrc))) {
                    break;
                }
                offset = 0;

            } while (1);
        } else if (CurrentLine.BeginsWith (key2, false)) { // actual tree strings & idents
            _StringBuffer buffer (128UL);
            if (!ReadNextNexusStatement (fState, f, CurrentLine, key2.length(), buffer, false, false, false,false,false, true)) {
                ReportWarning ("TREE construct not followed by ';'.");
                break;
            } else {
                // here goes the tree string in the form: treeID = treeString
                // pull the ID out first - check if it is a valid one
                // next crudely parse the tree string, extracting species names and
                pos = buffer.Find ('=',1,kStringEnd);
                if (pos==kNotFound) {
                    ReportWarning (buffer.Enquote () &" is not of the form TreeID = TreeString");
                } else {
                    long pos2 = buffer.FirstNonSpaceIndex (0,pos-1,kStringDirectionBackward);
                    if (pos2 != kNotFound) {
                        long j = buffer.FirstNonSpaceIndex (0,pos2-1,kStringDirectionForward);
                        if (j != kNotFound ) {
                            if (buffer.char_at (j) == '*') {
                                j = buffer.FirstNonSpaceIndex (j+1,pos2-1,kStringDirectionForward);
                                treeSelected = treeIdents.lLength;
                            }
                            if (j != kNotFound) {
                                _String nexus_tree_id (buffer,j,pos2),
                                        tree_id (nexus_tree_id);
                              
                                if (!nexus_tree_id.IsValidIdentifier(fIDAllowCompound)) {
                                  tree_id = nexus_tree_id.ConvertToAnIdent();
                                }
                                
                                tree_id = treeIdents.GenerateUniqueNameForList(GenerateUniqueObjectIDByType (nexus_tree_id, HY_BL_TREE) ,false);
                              
                                if (tree_id != nexus_tree_id) {
                                  ReportWarning(nexus_tree_id.Enquote('\'') & " has been renamed to " & tree_id.Enquote('\'') & " to avoid naming conflicts and/or comply with HyPhy ID requirements");
                                }

                              
                                treeIdents && & tree_id;
                                //  now get the rest of the tree string
                                pos2 = buffer.FirstNonSpaceIndex(pos2,pos+1, kStringDirectionBackward);
                                buffer.Trim (pos2,kStringEnd);
                                treeStrings && & buffer;
                            }
                        }
                        if (j == kNotFound) {
                             ReportWarning (_String("Could not find a tree identifier in:") & buffer.Enquote());
                        }
                    } else {
                        ReportWarning (buffer.Enquote () &" is not of the form TreeID = TreeString");
                    }
                }

            }
        } else {

           long offSet = 0L;

           ReportWarning (CurrentLine.Cut (0, CurrentLine.FirstSpaceIndex(1,kStringEnd)) & " is not used by HYPHY in TREES block");
           while (!done) {
                _StringBuffer buffer (128UL);
                done = ReadNextNexusStatement (fState, f, CurrentLine, offSet, buffer, false, false,true,false,false);
           }
           done = false;
        }

        if (!done) {
            //ReadNextLine(f,&CurrentLine,&fState,false);
            pos = 0;
        }
    }

    // now we shall check the string and match up node names with those present in the file

    for (long id = 0L; id<treeStrings.lLength; id++) {
        _String const * file_tree_string = (_String const *) treeStrings (id);
        long    treeLevel = 0L,
                lastNode,
                i = 0L;
      
        _StringBuffer revisedTreeString (128L);
      
      // TODO SLKP 20170621: looks like this is a generic Newick parser; why duplicate it here?
        for (i=0; i<file_tree_string->length(); ++i) {
            char    cc = file_tree_string->char_at (i);

            switch (cc) {
                case '(': { // creating a new internal node one level down
                  treeLevel++;
                  revisedTreeString<<'(';
                  break;
                }
                  
                case ',':
                case ')': { // creating a new node on the same level and finishes updating the list of parameters
                  if (cc==')') { // also create a new node on the same level
                    treeLevel--;
                  }
                  revisedTreeString<<cc;
                  break;
                }
                  
                case ':' : { // tree branch definition
                  lastNode = i+1;
                  revisedTreeString<<':';
                  char c = file_tree_string->char_at (lastNode);
                  
                  while (isdigit (c) || c=='.' || c=='-' || c=='e' || c=='E') {
                    if (lastNode<file_tree_string->length()) {
                      lastNode++;
                      revisedTreeString<<c;
                      c = file_tree_string->char_at (lastNode);
                    } else {
                      break;
                    }
                  }
                  i = lastNode-1;
                  break;
                }
                  
                default: { // node name
                  lastNode = i;
                  char c = file_tree_string->char_at (lastNode);
                  if (isspace (c)) {
                    break;
                  }
                  if (!(isalnum(c)||(c=='_'))) {
                    ReportWarning (_String("Node names should begin with a letter, a number, or an underscore: ") & error_conext (*file_tree_string, i) );
                    i = file_tree_string->length() +2;
                    break;
                  }
                  while (isalnum(c)||(c=='_')) {
                    if (lastNode<file_tree_string->length()) {
                      lastNode++;
                      c = file_tree_string->char_at (lastNode);
                    } else {
                      break;
                    }
                  }
                  _String node_name (*file_tree_string, i, lastNode-1);
                  i = lastNode-1;
                  lastNode = translationsFrom.BinaryFindObject (&node_name);
                  if (lastNode != kNotFound) {
                    revisedTreeString<< (_String*)translationsTo.list_data[lastNode];
                  } else {
                    revisedTreeString<< node_name;
                  }
                  break;
                }
          }
        }
        if (treeLevel) {
            ReportWarning (_String("Unbalanced '(,)' in the tree string:") & revisedTreeString.Enquote());
        } else if (i==file_tree_string->length()) {
            *((_String*)treeStrings.list_data[id]) = revisedTreeString;
        }
    }

    if (treeSelected < treeStrings.lLength) {
        hy_env :: EnvVariableSetNamespace(hy_env::data_file_tree, new HY_CONSTANT_TRUE,fState.theNamespace, false);
        hy_env :: EnvVariableSetNamespace(hy_env::data_file_tree_string, new _FString(*(_String*)treeStrings.list_data[treeSelected], false),fState.theNamespace, false);
     }

    if (treeStrings.lLength) {
        _StringBuffer initTreeMatrix (1024UL);

        initTreeMatrix   << hy_env::nexus_file_tree_matrix
                         << "={"
                         << _String ((long)treeStrings.lLength)
                         << ",2};\n";


        for (long id = 0; id < treeStrings.lLength; id++) {
            initTreeMatrix   << hy_env::nexus_file_tree_matrix
                             << '['
                             << _String (id)
                             << "][0]=\""
                             << (_String*)treeIdents(id)
                             << "\";\n"
                             << hy_env::nexus_file_tree_matrix
                             << '['
                             << _String (id)
                             << "][1]=\""
                             << (_String*)treeStrings(id)
                             << "\";\n";
        }

        _ExecutionList el (initTreeMatrix);
        el.Execute();
        terminate_execution = false;
    }
    SkipUntilNexusBlockEnd (fState, f,CurrentLine, pos);
}

//_________________________________________________________

void    ProcessNexusHYPHY (FileState& fState, long pos, FILE*file, _String& CurrentLine, _DataSet&) {
    static _String const endMark ("END;");
    _StringBuffer bfBody  (128UL);

    long      p2 = pos;
    pos = CurrentLine.FindAnyCase (endMark,pos+1,kStringEnd);

    fState.fileType = 0;

    if (pos != kNotFound) {
        bfBody << CurrentLine.Cut (p2,pos-1);
        CurrentLine.Trim(pos+endMark.length(),-1);
    } else {
        bfBody << CurrentLine.Cut (p2,-1);
        while (pos == kNotFound) {
            ReadNextLine(file,&CurrentLine,&fState,false,false);
            if (CurrentLine.empty()) {
                break;
            }

            pos = CurrentLine.FindAnyCase (endMark,0,kStringEnd);
            if (pos != kNotFound) {
                if (pos > 0) {
                    bfBody << CurrentLine.Cut (0,pos-1);
                }

                CurrentLine.Trim (pos+endMark.length(), -1);
                if (CurrentLine.empty()) {
                    ReadNextLine(file,&CurrentLine,&fState,false,false);
                }

                break;
            } else {
                bfBody << CurrentLine;
            }

        }
    }
     nexusBFBody = bfBody;

    fState.fileType = 3;

    CurrentLine = CurrentLine.ChangeCase(kStringUpperCase);

}

//_________________________________________________________

bool    ProcessNexusData (FileState& fState, long pos, FILE*f, _String& CurrentLine, _DataSet& result) {
    static const _String key1 ("DIMENSIONS"), key11 ("NTAX"), key12 ("NCHAR"),
            key2 ("FORMAT"),key21 ("DATATYPE"), key22 ("MISSING"), key23 ("GAP"), key24 ("SYMBOLS"),
            key25 ("EQUATE"), key26 ("MATCHCHAR"), key27 ("NOLABELS"), key28 ("INTERLEAVE"), key3 ("MATRIX"), keyEnd ("END");
  
    _String newAlph;

    bool    done = false,
            labels = true;

    char    charState = 0;

    _List   translations;
    char    missing = '?', gap = '-' , repeat = '.', charSwitcher;

    long    offSet = 0L, count, spExp = result.GetNames().lLength, sitesExp = 0;

    while (!done) {
        if (!FindNextNexusToken (fState, f, CurrentLine, pos)) {
            break;
        }

        if (CurrentLine.BeginsWith (keyEnd, false)) {
            pos = -1;
            break;
        }

        if (CurrentLine.BeginsWith (key1, false)) {
          // DIMENSIONS
            offSet = key1.length ();
            while (!done) {
                _StringBuffer buffer (128UL);
                done = ReadNextNexusStatement (fState, f, CurrentLine, offSet, buffer, true, true,true,false,false);

                if (buffer.BeginsWith(key11, false)) {
                    if (result.GetNames().lLength) {
                        ReportWarning ("NTAX will override the definition of taxa names from the TAXA block");
                    }
                    if (!(count=ReadNextNexusEquate (fState,f,CurrentLine, 0 ,buffer))) {
                        ReportWarning ("NTAX is not followed by '= number-of-taxa'");
                        done = true;
                    } else {
                        done = done||(count>1);
                        spExp = buffer.to_long();
                        if(spExp<=0L) {
                            ReportWarning ("NTAX must be a positive number");
                            done = true;
                            spExp = result.GetNames().lLength?result.GetNames().lLength:1;
                        }
                    }
                } else if (buffer.BeginsWith(key12, false)) {
                    if (!(count=ReadNextNexusEquate (fState,f,CurrentLine, 0 ,buffer))) {
                        ReportWarning ("NCHAR is not followed by '= number-of-charaters'");
                        done = true;
                    } else {
                        done = done||(count>1);
                        sitesExp = buffer.to_long();
                    }
                }
                offSet = 0L;
            }
            done = false;
        } else if (CurrentLine.BeginsWith (key2, false)) {
            // FORMAT
            offSet = key2.length();
            while (!done) {
                charSwitcher = 0;
                _StringBuffer buffer (128UL);
                done = ReadNextNexusStatement (fState, f, CurrentLine, offSet, buffer, true, true,true,false,false);
                offSet = 0L;
                buffer.Trim (buffer.FirstNonSpaceIndex(),kStringEnd);
                if (buffer.BeginsWith(key21)) { // datatype
                    if (!(count=ReadNextNexusEquate (fState,f,CurrentLine, 0 ,buffer))) {
                        ReportWarning ("DATATYPE is not followed by '= DNA|RNA|NUCLEOTIDE|PROTEIN|BINARY'");
                        done = true;
                    } else {
                        done = done||(count>1);
                        if ( buffer == _String("DNA") || buffer == _String("RNA") || buffer ==  _String("NUCLEOTIDE" )) {
                            if (newAlph.nonempty()) {
                               ReportWarning (_String("DNA|RNA|NUCLEOTIDE datatype directive will over-ride the custom symbols definition: ") & newAlph.Enquote());
                              newAlph.Clear();
                            }
                            if (done) {
                                done = false;
                                break;
                            }
                            continue;
                        } else if (buffer==_String("PROTEIN") || buffer == _String ("BINARY")) {
                            charState = 1+(buffer==_String("BINARY"));
                            if (newAlph.nonempty()) {
                                 newAlph = kEmptyString;
                                 ReportWarning (_String("PROTEIN|BINARY datatype directive will override the custom symbols definition: ") & newAlph.Enquote());
                                newAlph.Clear();
                            }
                            if (done) {
                                done = false;
                                break;
                            }
                            continue;
                        } else {
                            ReportWarning (buffer.Enquote() &" is not a recognized data type (DNA|RNA|NUCLEOTIDE|PROTEIN|BINARY are allowed).");
                            done = false;
                        }
                    }
                } else if (buffer.BeginsWith (key22, false)) { // MISSING
                    charSwitcher = 1;
                } else if (buffer.BeginsWith (key23, false)) { // GAP
                    charSwitcher = 2;
                } else if (buffer.BeginsWith (key26, false)) { // MATCHCHAR
                    charSwitcher = 3;
                } else if (buffer.BeginsWith (key27, false)) { // NOLABELS
                    labels = false;
                } else if (buffer.BeginsWith (key28, false)) { // INTERLEAVE
                    fState.interleaved = true;
                } else if (buffer.BeginsWith(key24, false)) { // SYMBOLS
                    count=ReadNextNexusEquate (fState,f,CurrentLine, 0 ,buffer, true,false);
                    if (buffer.empty()) {
                        ReportWarning (buffer.Enquote() & _String("is not of the form SYMBOLS = \"sym1 sym2 ...\". The entire block is ignored."));
                        done = true;
                        break;
                    }
                    _StringBuffer tempNewAlpha (128UL);
                    for (long pos1 = 0; pos1<buffer.length (); pos1++) {
                        charSwitcher = buffer.char_at (pos1);
                        if (!isspace(charSwitcher)) {
                            tempNewAlpha<<charSwitcher;
                        }
 
                    }
                    if (done) {
                        break;
                    }
                    newAlph = tempNewAlpha;
                    charSwitcher = 0;
                    done = done||(count>1);
                } else if (buffer.BeginsWith(key25, false)) { // EQUATE
                    buffer.Trim(key25.length(),kStringEnd);
                    if (!(count=ReadNextNexusEquate (fState,f,CurrentLine, 0,buffer,true,false))) {
                        ReportWarning (buffer.Enquote ()&" is not followed by '=char'");
                        done = true;
                    }
                    done = done||(count>1);
                    // blank now contains a full list of the form token=(token)
                    _String symbol, meaning;
                    bool    symbolDefined = false, meaningDefined = false;
                    for (count=0; count<buffer.length(); count++) {
                        charSwitcher = buffer.char_at (count);
                        if (isspace(charSwitcher)) {
                            continue;
                        } else if (charSwitcher == '=') {
                            if (symbolDefined&&!meaningDefined) {
                                meaningDefined = true;
                            }
                        } else
                            if (!symbolDefined) {
                                symbolDefined = true;
                                symbol = charSwitcher;
                                continue;
                            }
                            if (!meaningDefined) {
                                ReportWarning("EQUATE can only be used to define single-character tokens. Ignoring the EQUATE command.");
                                translations.Clear();
                                break;
                            }
                            meaning = meaning & charSwitcher;
                     }
                    if (symbol.length () && meaning.length () ) {
                        translations < new _String (symbol);
                        translations < new _String (meaning);
                    }
                    charSwitcher = 0;
                  buffer.Clear();
                }

                offSet = 0;
              
                _String built_in;
              
                if (charSwitcher) {
                    switch (charSwitcher) {
                      case 1:
                          built_in = "MISSING";
                          break;
                      case 2:
                          built_in = "GAP";
                          break;
                      case 3:
                          built_in = "MATCHCHAR";
                          break;
                    }
                    if (!(count=ReadNextNexusEquate (fState,f,CurrentLine, 0 ,buffer, true))) {
                        ReportWarning (buffer.Enquote() & " is not followed by '=char'");
                        done = true;
                    } else {
                        done = done||(count>1);
                        if (buffer.length () !=1) {
                            ReportWarning (buffer.Enquote() &" is not a valid " & built_in &" character.");
                        }
                    }
                    switch (charSwitcher) {
                      case 1:
                          missing = buffer.char_at (0);
                          if (gap == missing) {
                              gap = 0;
                          }
                          if (repeat == missing) {
                              repeat = 0;
                          }

                          break;
                      case 2:
                          gap = buffer.char_at (0);
                          if (missing == gap) {
                              missing = 0;
                          }
                          if (repeat == gap) {
                              repeat= 0;
                          }

                          break;
                      case 3:
                          repeat = buffer.char_at(0);
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
        } else if (CurrentLine.BeginsWith (key3, false)) { // matrix instruction
            // if needed, set up a new symbol set
            offSet = key3.length();
            if (newAlph.length()>1) { // a valid new alphabet set
                checkTTStatus (&fState);
                fState.translationTable->AddBaseSet (newAlph);
            } else {
                if (charState) {
                    checkTTStatus (&fState);
                    if (charState==1) {
                        newAlph = _TranslationTable::GetDefaultTable(HY_TRANSLATION_TABLE_PROTEIN);
                        fState.translationTable->baseLength = 20;
                    } else {
                        newAlph = _TranslationTable::GetDefaultTable(HY_TRANSLATION_TABLE_BINARY);
                        fState.translationTable->baseLength = 2;
                    }
                } else {
                    newAlph = _TranslationTable::GetDefaultTable(HY_TRANSLATION_TABLE_DNA);
                }
            }
            // set up translations
            if (translations.lLength) {
                checkTTStatus (&fState);
            }

            for (long k = 0; k<translations.lLength; k+=2) {
                char c = ((_String*)translations(k))->char_at (0);
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
                result.ClearNames();
            }


            while (1) {
              
              _StringBuffer buffer   (128L),
                            buffer_2 (128L),
                            * source;
              
                done = ReadNextNexusStatement (fState, f, CurrentLine, offSet?offSet+1:0, buffer, true, true,true,!labels,false);
                offSet = 0;
                // in each line that should produce first the name of the taxon
                // and then the data string for the taxon

                if (labels) {
                    if (result.GetNames().lLength<spExp) {
                        if (spExp>0 && buffer.empty ()) {
                            ReportWarning (_String("Could not find NTAX taxon names in the matrix. Read: ")&_String((long)result.GetNames().lLength) & " sequences.");
                            break;
                        }

                        if (!(sitesExp&&fState.curSite&&(fState.curSite<sitesExp)&&(!fState.interleaved))) {
                            result.AddName(buffer);
                            fState.totalSpeciesExpected++;
                        }
                    } else {
                        if (done) {
                            break;
                        }
                    }

                    if (!(sitesExp&&fState.curSite&&(fState.curSite<sitesExp)&&(!fState.interleaved))) {
                        done = ReadNextNexusStatement (fState, f, CurrentLine, offSet, buffer_2, true, true,true,true,false);
                        source = &buffer_2;
                    } else {
                        source = &buffer;
                    }
                } else {
                    if (loopIterations<spExp) {
                        if (!(sitesExp&&fState.curSite&&(fState.curSite<sitesExp)&&(!fState.interleaved))) {
                            fState.totalSpeciesExpected++;
                        } else {
                            loopIterations --;
                        }
                    }
                    source = &buffer;
                }

                if (source->empty()) {
                    ReportWarning (_String("Could not find NTAX data strings in the matrix. Read: ")&_String((long)result.GetNames().lLength) & " sequences.");
                    break;
                }
                loopIterations++;
                ISelector (fState, *source, result);

                if (done)
                    if (loopIterations>=fState.totalSpeciesExpected) {
                        break;    // finished reading
                    }

                if  ((f&&feof(f))||(fState.theSource&&(fState.theSource->length()<=fState.pInSrc))) {
                    break;
                }
            }


            if (result.GetNames().lLength!=spExp) {
                ReportWarning(_String ("Expected ")&spExp&" taxa, but found "&(long)result.GetNames().lLength);
            }
            if (result.lLength!=sitesExp && result.InternalStorageMode() == 0) {
                ReportWarning(_String ("Expected ")&sitesExp&" sites, but found "&(long)result.lLength);
            }
            if (spExp && loopIterations%spExp) {
                ReportWarning(_String ("There is an inconsistency between NTAX and the number of data strings in the matrix"));
            }
            done = true;
        } else {
            ReportWarning (CurrentLine.Cut (0, CurrentLine.FirstSpaceIndex(1,kStringEnd)) & " is not used by HYPHY");
            while (!done) {
                _StringBuffer buffer (128L);
                done = ReadNextNexusStatement (fState, f, CurrentLine, offSet, buffer, true, false,true,false,false);
            }
            done = false;
        }
        if (!done) {
            if (CurrentLine.empty ()) {
                ReadNextLine(f,&CurrentLine,&fState,false);
            }
            pos = 0;
            if (CurrentLine.empty () ) {
                done = true;
            }
        }
    }

    SkipUntilNexusBlockEnd (fState, f,CurrentLine, pos);
    return true;
}

//_________________________________________________________

void    ReadNexusFile (FileState& fState, FILE*file, _DataSet& result) {
    bool   dataRead = false, lookForEnd = false;
    long   f,g, file_line = fState.currentFileLine;

    fState.fileType = 3; // NEXUS
    static const _String beginMark ("BEGIN"), endMark ("END"), data ("DATA"), chars ("CHARACTERS"),
            taxa ("TAXA"), trees ("TREES"), assumptions ("ASSUMPTIONS"), hyphy ("HYPHY"), sets ("SETS");
  
    _String CurrentLine, blockName;

    ReadNextLine(file,&CurrentLine,&fState,false);
    while (CurrentLine.nonempty()) {
        f = 0;
        /** TODO SLKP 20180921 : if any of the commands loads a new CurrentLine, the marker 'f' needs to be reset
            but currently we have no way of knowing whether or not a new line was loaded.
            For the time-being fixing by adding a line # tracker for fState
         */
        while ((f = CurrentLine.FindAnyCase(beginMark,file_line == fState.currentFileLine ? f : 0L,kStringEnd ))>=0) {
           file_line = fState.currentFileLine;
          
            f = CurrentLine.FirstNonSpaceIndex (f+beginMark.length(),kStringEnd,kStringDirectionForward);
            if (f!=-1) { // process
                g = CurrentLine.Find (';', f, -1);
                if (g!=kNotFound) {
                    blockName = CurrentLine.Cut (f,g-1);
                    // dispatch to block readers
                    if (blockName.EqualIgnoringCase(data)) {
                        ReportWarning (blockName.Enquote() & " block is now deprecated in NEXUS and should not be used.");

                        if (!dataRead) {
                            dataRead = ProcessNexusData (fState, g+1, file, CurrentLine, result);
                        }
                        //SkipUntilNexusBlockEnd (fState,file,CurrentLine,f);

                        else {
                            ReportWarning ("Only one data set per NEXUS file is read by ReadDataSet - the 1st valid one.");
                        }
                    } else if (blockName.EqualIgnoringCase(taxa)) {
                        if (!dataRead) {
                            ProcessNexusTaxa (fState, g+1, file, CurrentLine, result);
                        } else {
                            ReportWarning ("The TAXA block was encountered after CHARACTER had been read and will be ignored.");
                        }
                    } else if (blockName.EqualIgnoringCase(trees)) {
                        ProcessNexusTrees (fState, g+1, file, CurrentLine, result);
                    } else if (blockName.EqualIgnoringCase(chars)) {
                        if (!dataRead) {
                            dataRead = ProcessNexusData (fState, g+1, file, CurrentLine, result);
                        } else {
                            ReportWarning ("Only one data set per NEXUS file is read by ReadDataSet - the 1st valid one.");
                        }
                    } else if (blockName.EqualIgnoringCase(assumptions)||blockName.EqualIgnoringCase(sets)) {
                        ProcessNexusAssumptions (fState, g+1, file, CurrentLine, result);
                    } else if (blockName.EqualIgnoringCase(hyphy)) {
                        ProcessNexusHYPHY (fState, g+1, file, CurrentLine, result);
                    } else {
                        ReportWarning (_String("NEXUS blocks ")&blockName.Enquote()&(" are not used by HYPHY."));
                        lookForEnd = true;
                        break;
                        // now look for the end of this block
                    }

                } else {
                    break;
                }
            } else {
                ReportWarning (_String ("NEXUS BEGIN must be followed by the name of the block. Skipping until next BEGIN statement."));
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
