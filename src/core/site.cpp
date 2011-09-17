/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2009
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon              (apoon@cfenet.ubc.ca)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#define  HYPHY_SITE_DEFAULT_BUFFER_SIZE 256

#include "likefunc.h"

#define   DATA_SET_SWITCH_THRESHOLD     100000

#include "math.h"

#ifdef __MAC__
extern bool handleGUI(bool);
#endif

#if !defined  __UNIX__ && !defined __HEADLESS__
#include "HYDataPanel.h"
#endif

#if !defined  __UNIX__ || defined __HEADLESS__
#include "preferences.h"
#endif

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

#ifdef    __HYPHYMPI__
extern int _hy_mpi_node_rank;
#endif


_TranslationTable      defaultTranslationTable;

//_________________________________________________________

extern _Parameter dFPrintFormat,
       dFDefaultWidth;

//_________________________________________________________

_String           dataFileTree              ("IS_TREE_PRESENT_IN_DATA"),
                  dataFileTreeString       ("DATAFILE_TREE"),
                  aminoAcidOneCharCodes    ("ACDEFGHIKLMNPQRSTVWY"),
                  dnaOneCharCodes          ("ACGT"),
                  rnaOneCharCodes          ("ACGU"),
                  binaryOneCharCodes       ("01"),
                  nexusFileTreeMatrix      ("NEXUS_FILE_TREE_MATRIX"),
                  dataFilePartitionMatrix   ("DATA_FILE_PARTITION_MATRIX"),
                  useTraversalHeuristic        ("USE_TRAVERSAL_HEURISTIC"),
                  defaultLargeFileCutoff   ("USE_MEMORY_SAVING_DATA_STRUCTURES"),
                  fileTreeString;

//_________________________________________________________

/* function declarations */

void    checkTTStatus (FileState* fs);
void    processCommand (_String*s, FileState*fs);
void    FilterRawString (_String& s, FileState* fs, _DataSet & ds);
long    ProcessLine (_String&s , FileState *fs, _DataSet& ds);
void    PadLine (FileState& fState, _DataSet& result);
void    ISelector (FileState& fState, _String& CurrentLine, _DataSet& result);
bool    SkipLine (_String& theLine, FileState* fS);
void    TrimPhylipLine (_String& CurrentLine, _DataSet& ds);
void    ProcessTree    (FileState*, FILE*, _String&);
void    ReadNexusFile (FileState& fState, FILE*f, _DataSet& result);

//_________________________________________________________
_TranslationTable::_TranslationTable (void)
{
    baseLength = 4;
    checkTable = NULL;
}

//_________________________________________________________
_TranslationTable::_TranslationTable (char baseL)
{
    baseLength = (baseL==20)?20:4;
    checkTable = NULL;
}

//_________________________________________________________
_TranslationTable::_TranslationTable (_TranslationTable& t)
{
    tokensAdded = t.tokensAdded;
    baseLength = t.baseLength;
    baseSet = t.baseSet;
    translationsAdded.Duplicate (&t.translationsAdded);
    checkTable = NULL;
}

//_________________________________________________________
_TranslationTable::_TranslationTable (_String& alphabet)
{
    baseLength = alphabet.sLength;
    checkTable = NULL;
    if (!(alphabet.Equal (&dnaOneCharCodes) || alphabet.Equal (&rnaOneCharCodes) ||
            alphabet.Equal (&binaryOneCharCodes) || alphabet.Equal (&aminoAcidOneCharCodes))) {
        AddBaseSet (alphabet);
    }
}

//_________________________________________________________
BaseRef     _TranslationTable::makeDynamic (void)
{
    _TranslationTable * r = new _TranslationTable;
    checkPointer(r);

    memcpy ((char*)r, (char*)this, sizeof (_TranslationTable));
    r->nInstances = 1;
    r->tokensAdded.Duplicate (&tokensAdded);
    r->baseSet.Duplicate (&baseSet);
    r->translationsAdded.Duplicate (&translationsAdded);
    r->checkTable = NULL;
    return r;
}

//_________________________________________________________
long    _TranslationTable::TokenCode (char token)
{
    // standard translations
    long * receptacle       = new long[baseLength];
    if (!receptacle) {
        checkPointer    (receptacle);
    }
    TokenCode               (token,receptacle);

    long                    theCode         = 0,
                            shifter       = 1;

    for (int i = 0; i<baseLength; i++, shifter <<= 1) {
        theCode +=  shifter*receptacle[i];
    }

    delete receptacle;
    return theCode;
}

//_________________________________________________________
char    _TranslationTable::CodeToLetter (long* split)
// assumes a non-unique translation of split
// for unique - use ConvertCodeToLetters
{
    long     shifter = 1,
             trsl    = 0;

    for (long k=0; k<(baseSet.sLength?baseSet.sLength:baseLength); k++,shifter<<=1) {
        trsl+=shifter*split[k];
    }

    if (baseSet.sLength == 0)
        // one of the standard alphabers
    {
        if (baseLength==4)
            // nucleotides
        {
            switch (trsl) {
            case 3:
                return 'M';
            case 5:
                return 'S';
            case 6:
                return 'R';
            case 7:
                return 'V';
            case 9:
                return 'W';
            case 10:
                return 'Y';
            case 11:
                return 'H';
            case 12:
                return 'K';
            case 14:
                return 'B';
            }
        } else if (baseLength==20)
            // amino acids
        {
            switch (trsl) {
            case 2052:
                return 'B';
            case 8200:
                return 'Z';
            }
        }
    } else if (tokensAdded.sLength) {
        shifter = translationsAdded.Find(trsl);
        // linear search for (binary) translations
        if (shifter>=0) {
            return tokensAdded.sData[shifter];
        }
    }
    return '?';
}

//_________________________________________________________
void    _TranslationTable::SplitTokenCode (long code, long* receptacle)
{
    long shifter = 1;
    for (int i=0; i<baseLength; i++) {
        receptacle[i] = ((code&shifter)!=0);
        shifter*=2;
    }
}

//_________________________________________________________
long    _TranslationTable::LengthOfAlphabet (void)
{
    return baseSet.sLength?baseSet.sLength:baseLength;
}

//_________________________________________________________

bool    _TranslationTable::TokenCode (char token, long* receptacle, bool gapToOnes)
{

    long f = tokensAdded.sLength?tokensAdded.Find (token):-1;
    // check for custom translations
    // OPTIMIZE FLAG linear search:
    // SLKP 20071002 should really be a 256 char lookup table

    if (f != -1) {
        SplitTokenCode(translationsAdded.lData[f], receptacle);
        return true;
    }

    if (baseSet.sLength)
        // custom base alphabet
    {
        for (long k=0; k<baseLength; k++) {
            receptacle[k] = 0;
        }

        f = baseSet.Find(token);
        // OPTIMIZE FLAG linear search:
        // SLKP 20071002 should really be a 256 char lookup table

        if (f!=-1) {
            receptacle[f] = 1;
        }

        return true;
    }

    if (baseLength==4)
        // standard nucleotide
    {
        receptacle[0]=0;
        receptacle[1]=0;
        receptacle[2]=0;
        receptacle[3]=0;

        switch (token) {
        case 'A':
            receptacle[0]=1;
            break;

        case 'C':
            receptacle[1]=1;
            break;

        case 'G':
            receptacle[2]=1;
            break;

        case 'T':
        case 'U':
            receptacle[3]=1;
            break;

        case 'Y':
            receptacle[3]=1;
            receptacle[1]=1;
            break;

        case 'R':
            receptacle[0]=1;
            receptacle[2]=1;
            break;

        case 'W':
            receptacle[3]=1;
            receptacle[0]=1;
            break;

        case 'S':
            receptacle[1]=1;
            receptacle[2]=1;
            break;

        case 'K':
            receptacle[3]=1;
            receptacle[2]=1;
            break;

        case 'M':
            receptacle[1]=1;
            receptacle[0]=1;
            break;

        case 'B':
            receptacle[1]=1;
            receptacle[2]=1;
            receptacle[3]=1;
            break;

        case 'D':
            receptacle[0]=1;
            receptacle[2]=1;
            receptacle[3]=1;
            break;

        case 'H':
            receptacle[1]=1;
            receptacle[0]=1;
            receptacle[3]=1;
            break;

        case 'V':
            receptacle[1]=1;
            receptacle[2]=1;
            receptacle[0]=1;
            break;

        case 'X':
        case 'N':
        case '?':
        case '.':
        case '*':
            receptacle[1]=1;
            receptacle[2]=1;
            receptacle[3]=1;
            receptacle[0]=1;
            break;

        case '-':
            if (gapToOnes) {
                receptacle[1]=1;
                receptacle[2]=1;
                receptacle[3]=1;
                receptacle[0]=1;
                break;
            }
        }
    } else {
        if (baseLength==20) {
            for (long k=0; k<baseLength; k++) {
                receptacle[k] = 0;
            }

            switch (token) {
            case 'A':
                receptacle[0]=1;
                break;

            case 'B':
                receptacle[2]=1;
                receptacle[11]=1;
                break;

            case 'C':
                receptacle[1]=1;
                break;

            case 'D':
                receptacle[2]=1;
                break;

            case 'E':
                receptacle[3]=1;
                break;

            case 'F':
                receptacle[4]=1;
                break;

            case 'G':
                receptacle[5]=1;
                break;

            case 'H':
                receptacle[6]=1;
                break;

            case 'I':
                receptacle[7]=1;
                break;

            case 'K':
                receptacle[8]=1;
                break;

            case 'L':
                receptacle[9]=1;
                break;

            case 'M':
                receptacle[10]=1;
                break;

            case 'N':
                receptacle[11]=1;
                break;

            case 'P':
                receptacle[12]=1;
                break;

            case 'Q':
                receptacle[13]=1;
                break;

            case 'R':
                receptacle[14]=1;
                break;

            case 'S':
                receptacle[15]=1;
                break;

            case 'T':
                receptacle[16]=1;
                break;

            case 'V':
                receptacle[17]=1;
                break;

            case 'W':
                receptacle[18]=1;
                break;

            case 'Y':
                receptacle[19]=1;
                break;

            case 'Z':
                receptacle[3]=1;
                receptacle[13]=1;
                break;

            case 'X':
            case '?':
            case '.':
            case '*': {
                for (int j = 0; j<20; j++) {
                    receptacle[j] = 1;
                }
            }
            break;
            case '-': {
                if (gapToOnes)
                    for (int j = 0; j<20; j++) {
                        receptacle[j] = 1;
                    }
            }
            break;
            }
        } else
            // binary
        {
            receptacle[0] = 0;
            receptacle[1] = 0;
            switch (token) {
            case '0':
                receptacle[0]=1;
                break;

            case '1':
                receptacle[1]=1;
                break;

            case 'X':
            case '?':
            case '.':
            case '*': {
                receptacle[0] = 1;
                receptacle[1] = 1;
            }
            break;
            case '-': {
                if (gapToOnes) {
                    receptacle[0] = 1;
                    receptacle[1] = 1;
                }
            }
            break;
            }

        }
    }
    return false;

}
//_________________________________________________________
void    _TranslationTable::PrepareForChecks (void)
{
    if (checkTable == NULL) {
        checkTable = MemAllocate (256);
    }

    for (long i2=0; i2<256; i2++) {
        checkTable[i2]=0;
    }

    _String checkSymbols;
//  if (baseLength == 4)
//      checkSymbols = _String("ACGTUYRWSKMBDHVXN?O-.")&tokensAdded;
    if (baseSet.sLength) {
        checkSymbols = baseSet&tokensAdded;
    } else if (baseLength == 2) {
        checkSymbols = _String("01*?-.")&tokensAdded;
    } else {
        checkSymbols = _String("ABCDEFGHIJKLMNOPQRSTUVWXYZ*?-.")&tokensAdded;
    }

    for (long i=0; i<checkSymbols.sLength; i++) {
        checkTable[checkSymbols(i)] = 1;
    }
}

//_________________________________________________________
bool    _TranslationTable::IsCharLegal (char c)
{
    if (!checkTable) {
        PrepareForChecks();
    }
    return checkTable[c];
}
//___________________________________________

void    _TranslationTable::AddTokenCode (char token, _String& code)
{
    long    f,
            newCode = 0;

    bool    killBS = false;

    if (baseSet.sLength==0)
        // fill in baseSet for standard alphabets
    {
        if (baseLength == 4) {
            baseSet = dnaOneCharCodes;
        } else if (baseLength == 20) {
            baseSet = aminoAcidOneCharCodes;
        } else {
            baseSet = binaryOneCharCodes;
        }
        killBS = true;
    }


    if (baseSet.sLength) {
        long shifter = 1;
        for (int j = 0; j<baseSet.sLength; j++, shifter*=2)
            if (code.Find (baseSet.sData[j])>=0) {
                newCode += shifter;
            }
    }

    f = baseSet.Find (token);
    if (killBS) {
        baseSet = empty;
    }
    if (f>=0) {
        return;
    }
    // see if the character being added is a base
    // character; those cannot be redefined

    f = tokensAdded.Find (token,0,-1);
    // new definition or redefinition?

    if (f==-1) { // new
        tokensAdded             = tokensAdded&token;
        translationsAdded       << 0;
        f                       = tokensAdded.sLength-1;
    }

    translationsAdded.lData[f] = newCode;
}

//_________________________________________________________

void    _TranslationTable::AddBaseSet (_String& code)
{
    baseSet         = code;
    baseSet.StripQuotes();
    baseLength      = baseSet.sLength;
    if (baseLength > HY_WIDTH_OF_LONG)
        // longer than the bit size of 'long'
        // can't handle those
    {
        _String err = _String ("Alphabets with more than ")
                      & HY_WIDTH_OF_LONG &
                      " characters are not supported";
        WarnError (err);
    }

}

//_________________________________________________________

char    _TranslationTable::GetSkipChar (void)
{
    if ( baseSet.sLength==0 && translationsAdded.lLength==0 ) {
        return '?';    // this is the default
    }

    // see if there is a symbol
    // which maps to all '1'

    long    all     = 0,
            ul       = baseSet.sLength?baseSet.sLength:baseLength,
            shifter = 1;

    for  (long f=0; f<ul; f++, shifter <<= 1) {
        all |= shifter;
    }

    if  ((all = translationsAdded.Find(all))==-1) {
        return '?';
    } else {
        return tokensAdded[all];
    }
}

//_________________________________________________________

char    _TranslationTable::GetGapChar (void)
{
    if ( baseSet.sLength==0 && translationsAdded.lLength==0 ) {
        return '-';    // default gap character
    }

    long f = translationsAdded.Find(0L);

    if  (f==-1) {
        return 0;
    } else {
        return tokensAdded[f];
    }
}

//_________________________________________________________
_String _TranslationTable::ConvertCodeToLetters (long code, char base)
{

    _String res (base,false);
    if (code >= 0) {
        // OPTIMIZE FLAG; repeated memory allocation/deallocation
        if (baseSet.sLength)
            for (long k=1; k<=base; k++, code/=baseLength) {
                res.sData[base-k]=baseSet.sData[code%baseLength];
            }
        else if (baseLength==4) {
            for (long k=1; k<=base; k++, code/=baseLength) {
                switch (code%baseLength) {
                case 0:
                    res[base-k]='A';
                    break;
                case 1:
                    res[base-k]='C';
                    break;
                case 2:
                    res[base-k]='G';
                    break;
                case 3:
                    res[base-k]='T';
                    break;
                }
            }
        } else if (baseLength == 20) {
            for (long k=1; k<=base; k++, code/=baseLength) {
                char out = code%baseLength;
                if (out==0) {
                    res[base-k] = 'A';
                } else if (out<=7) {
                    res[base-k] = 'B'+out;
                } else if (out<=11) {
                    res[base-k] = 'C'+out;
                } else if (out<=16) {
                    res[base-k] = 'D'+out;
                } else if (out<=18) {
                    res[base-k] = 'E'+out;
                } else {
                    res[base-k]='Y';
                }
            }
        } else if (baseLength == 2)
            for (long k=1; k<=base; k++, code/=baseLength) {
                switch (code%baseLength) {
                case 0:
                    res[base-k]='0';
                    break;
                case 1:
                    res[base-k]='1';
                    break;
                }
            }
    } else {
        char c = GetGapChar();
        for (long k=0; k<base; k++) {
            res.sData[k] = c;
        }
    }
    return res;


}

//_________________________________________________________

_TranslationTable*  _TranslationTable::MergeTables (_TranslationTable* table2)
// merge the translation tables if they are compatible, return the result,
// otherwise return nil
{
    if (baseSet.sLength==table2->baseSet.sLength) {
        if (baseSet.sLength==0) { // standard alphabet
            if (baseLength!=table2->baseLength) {
                return nil;
            }
        } else if(!(baseSet.Equal (&table2->baseSet))) {
            return nil;
        }

        _TranslationTable* result = new _TranslationTable (*this);
        checkPointer     (result);
        if (table2->tokensAdded.sLength) {
            for (long i=0; i<table2->tokensAdded.sLength; i++) {
                long f = tokensAdded.Find (table2->tokensAdded[i]);
                if (f==-1) {
                    result->tokensAdded       && table2->tokensAdded[i];
                    // SLKP 20071002 added the next line;
                    // was not adding the translation for the new token
                    result->translationsAdded << table2->translationsAdded[i];
                } else if (translationsAdded.lData[f] != table2->translationsAdded.lData[i]) {
                    DeleteObject (result);
                    return nil;
                }
            }
            return result;
        } else {
            return result;
        }
    }
    return nil;
}

//_________________________________________________________

_Site::_Site (void):_CString(16,true)
{
    refNo = -1;
}

//_________________________________________________________
_Site::_Site (_String& s):_CString (s.sLength, true)
{
    refNo = -1;
    (*this)<<&s;
}

//_________________________________________________________
_Site::_Site (char s):_CString (16, true)
{
    refNo = -1;
    (*this)<<s;
}

//_________________________________________________________
_Site::_Site (long s)
//:_CString (1, true)
{
    SetRefNo (s);
}

//_________________________________________________________
_Site::~_Site (void)
{}

//_________________________________________________________
void    _Site::Complete (void)
{
    if (refNo==-1) {
        _String::Finalize();
    }

    refNo = refNo<0?-refNo:refNo;
}
//_________________________________________________________
BaseRef _Site::makeDynamic(void)
{
    _Site * r = new _Site;
    checkPointer(r);

    memcpy ((char*)r, (char*)this, sizeof (_Site));
    r->nInstances = 1;
    nInstances++;
    return r;
}

//_______________________________________________________________________
void    _Site::Duplicate (BaseRef ref)
{
    _Site * s = (_Site*)ref;
    sLength = s->sLength;
    if (sData) {
        free(sData);
    }
    sData      = s->sData;
    allocatedSpace = s->allocatedSpace;
    //nInstances = ref->nInstances;
    if (sData) {
        /*long theLength = sLength/storageIncrement;
        if (!sLength||sLength%storageIncrement) theLength++;
        theLength*=storageIncrement;
        checkPointer (sData = (char*)MemAllocate (theLength));
        memcpy (sData, s->sData, sLength);*/
        if (allocatedSpace) {
            checkPointer (sData = (char*)MemAllocate (allocatedSpace*sizeof(char)));
        } else {
            checkPointer (sData = (char*)MemAllocate (sLength*sizeof(char)));
        }
        memcpy (sData, s->sData, sLength);
    }
    refNo = -1;
}

//_______________________________________________________________________
void    _Site::Clear (void)
{
    if (sData) {
        free(sData);
        sData = NULL;

        //nInstances = 0;
    }
    allocatedSpace = 0;
    sLength = 0;
}
//_______________________________________________________________________
void    _Site::PrepareToUse (void)
{
    if (IsCompressed()) {
        _String * s = Decompress();
        DuplicateErasing(s);;
        DeleteObject (s);
        SetDecompressed();
    }
}
//_______________________________________________________________________
void    _Site::Archive (void)
{
    if ((!IsCompressed())&&(GetRefNo()>=0)) {
        BestCompress (NUCLEOTIDEALPHABET);
    }
}


//_______________________________________________________________________
//  _DataSet Functions
//_______________________________________________________________________

_DataSet::_DataSet (void)
{
    theTT                = &defaultTranslationTable;
    streamThrough        = nil;
    dsh                  = nil;
    useHorizontalRep     = false;
}

//_______________________________________________________________________

_DataSet::_DataSet (long l):_List((unsigned long)l),theFrequencies((unsigned long)l) // with estimated number of sites per file
{
    dsh                 = nil;
    streamThrough       = nil;
    theTT               = &defaultTranslationTable;
    useHorizontalRep     = false;
}

//_______________________________________________________________________

_DataSet::_DataSet (FILE *f)
{
    dsh                  = nil;
    useHorizontalRep     = false;
    theTT                = &defaultTranslationTable;
    streamThrough        = f;
    theMap << 0; // current sequence
    theMap << 0; // current site
    theMap << 0; // total sites
}

//_______________________________________________________________________

_DataSet::~_DataSet (void)
{
    if (theTT != &defaultTranslationTable) {
        DeleteObject (theTT);
    }
}

//_______________________________________________________________________

void _DataSet::Clear (void)
{
    _List::Clear();
    theMap.Clear();
    theFrequencies.Clear();
    theNames.Clear();
    if (theTT != &defaultTranslationTable) {
        DeleteObject (theTT);
        theTT = &defaultTranslationTable;
    }
    noOfSpecies = 0;
    if (dsh) {
        dsh->incompletePatterns->Clear(false);
        delete (dsh);
        dsh = nil;
    }
    useHorizontalRep     = false;
}

//_______________________________________________________________________

BaseRef _DataSet::makeDynamic (void)
{
    _DataSet * r = new _DataSet;
    checkPointer(r);
    memcpy ((char*)r, (char*)this, sizeof (_DataSet));
    r->nInstances = 1;
    r->theMap.Duplicate (&theMap);
    r->theFrequencies.Duplicate (&theFrequencies);
    if (theTT!=&defaultTranslationTable) {
        r->theTT->nInstances++;
    }
    r->theNames.Duplicate (&theNames);
    r->streamThrough = streamThrough;
    nInstances++;
    r->dsh = nil;
    r->useHorizontalRep      = false;
    return r;
}

//_______________________________________________________________________

void     _DataSet::ResetIHelper (void)
{
    if (dsh && dsh->characterPositions.lLength == 256)
        for (long k=0; k<256; k++) {
            dsh->characterPositions.lData[k] = -1;
        }
}

//_______________________________________________________________________

void     _DataSet::ConvertRepresentations (void)
{
    if (useHorizontalRep == false) {
        _List horStrings;

        if (lLength == 0) {
            AppendNewInstance (new _Site);
        } else {
            _Site * aSite = (_Site*)lData[0];

            for (long str = 0; str < aSite->sLength; str++) {
                _String * aString = new _String (DATA_SET_SWITCH_THRESHOLD,true);
                horStrings << aString;
                aString->nInstances --;
            }

            for  (long s = 0; s < lLength; s++) {
                _Site * aSite = (_Site*)lData[s];
                if (aSite->sLength>horStrings.lLength || aSite->GetRefNo() != -1) {
                    FlagError ("Irrecoverable internal error in _DataSet::ConvertRepresentations. Sorry about that.");
                    return;
                }
                aSite->Finalize();
                for (long s2 = 0; s2 < aSite->sLength; s2++) {
                    (*(_String*)horStrings.lData[s2]) << aSite->sData[s2];
                }
            }

            _List::Clear();
            theFrequencies.Clear();
            {
                for  (long s = 0; s < horStrings.lLength; s++) {
                    (*this) << horStrings(s);
                }
            }
        }
        useHorizontalRep = true;
    }
}

//_______________________________________________________________________

void     _DataSet::AddSite (char c)
{
    if (streamThrough) {
        if (theMap.lData[0] == 0) {
            if (theMap.lData[1] == 0) {
                if (theNames.lLength) {
                    fprintf (streamThrough,">%s\n",((_String*)theNames(0))->getStr());
                } else {
                    fprintf (streamThrough,">Sequence 1\n");
                }
                (*this) && & empty;
            }

            theMap.lData[1]++;
            theMap.lData[2]++;
            fputc (c,streamThrough);
        } else {
            WarnError ("Can't add more sites to a file based data set, when more that one sequence has been written!");
        }
    } else {
        if (useHorizontalRep == false) {
            if (lLength < DATA_SET_SWITCH_THRESHOLD) {
                _Site* nC = new _Site(c);
                checkPointer(nC);
                theFrequencies<<1;
                (*this)<<nC;
                nC->nInstances --;
                return;
            } else {
                ConvertRepresentations ();
            }
        }

        (*((_String*)lData[0])) << c;

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

void     _DataSet::Write2Site (long index, char c)
{
    if (streamThrough) {
        if (index == 0) {
            if (theMap.lData[2] == theMap.lData[1]) {
                theMap.lData[0] ++;

                if (theNames.lLength > theMap.lData[0]) {
                    fprintf (streamThrough,"\n>%s\n",((_String*)theNames(theMap.lData[0]))->getStr());
                } else {
                    fprintf (streamThrough,"\n>Sequence %ld\n",theMap.lData[0]+1);
                }

                theMap.lData[1] = 0;
            } else {
                WarnError ("Can't write sequences of unequal lengths to a file based data set.");
                return;
            }
        } else if (index != theMap.lData[1]) {
            WarnError ("Can't write sites which are not consecutive to a file based data set.");
            return;
        }

        theMap.lData[1] ++;
        fputc (c,streamThrough);
    } else {
        /*if (!dsh)
        {
            WarnError ("Internal Error in 'Write2Site' - called Write2Site before any AddSite calls");
            return;
        }*/

        if (useHorizontalRep) {
            long     currentWritten = ((_String*)lData[0])->sLength;

            if (index>=currentWritten) {
                WarnError ("Internal Error in 'Write2Site' - index is too high (using compact representation)");
                return;
            } else {
                if (index == 0) {
                    _String * newString = new _String (currentWritten,true);
                    (*newString) << c;
                    (*this) << newString;
                    newString->nInstances --;
                } else {
                    long s = 1;
                    for (; s<lLength; s++) {
                        _String *aString = (_String*)lData[s];
                        if (aString->sLength == index) {
                            (*aString) << c;
                            break;
                        }
                    }
                    if (s == lLength) {
                        WarnError ("Internal Error in 'Write2Site' - no appropriate  string to write too (compact representation)");
                        return;
                    }
                }
            }
        } else {
            if (index>=lLength) {
                WarnError ("Internal Error in 'Write2Site' - index is too high");
                return;
            }
            _Site* s = (_Site*)lData[index];
            long rN = s->GetRefNo();
            if (rN==-1) { // independent site
                //dsh->incompletePatterns->Delete (s,false);
                (*s)<<c;
                //dsh->incompletePatterns->Insert (s,index);
            } else {
                _Site *ss = (_Site*)lData[rN];
                long sL = ss->sLength-1;
                if (ss->sData[sL] != c) { // appending distinct char
                    s->Duplicate (ss);
                    s->sData[sL]=c;
                    theFrequencies.lData[rN]--;

                    rN = dsh->incompletePatterns->Find (s);
                    if (rN>=0) {
                        rN =  dsh->incompletePatterns->GetXtra (rN);
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
                        dsh->incompletePatterns->Insert (s,index);
                    }
                }
            }
        }
    }
}



//_______________________________________________________________________

void     _DataSet::CheckMapping (long index)
{
    if (index>=lLength) {
        FlagError ("Internal Error in 'CheckMapping' - index is too high");
    }

    _Site* s = (_Site*)lData[index];

    for (long k = 0; k < index; k ++) {
        _Site* ss = (_Site*)lData[k];
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

long     _DataSet::GetCharDimension         (void)  // return the size of the alphabet space
{
    return theTT->baseLength;
}

//_______________________________________________________________________

long     _DataSet::GetNoTypes (void)  // return the number of unique columns
{
    return theMap.countitems ();
}
//_______________________________________________________________________

long     _DataSet::GetFreqType (long index)  // return the frequency of a site
{
    return theFrequencies(theMap(index));
}
//_______________________________________________________________________

void    _DataSet:: SetTranslationTable (_DataSet *  newTT )
{
    if (theTT&&(theTT!= &defaultTranslationTable)) {
        DeleteObject (theTT);
    }
    theTT = (_TranslationTable*)newTT->theTT->makeDynamic();

}

//_______________________________________________________________________

void    _DataSet:: SetTranslationTable (_TranslationTable *  newTT )
{
    if (theTT&&(theTT!= &defaultTranslationTable)) {
        DeleteObject (theTT);
    }
    theTT = (_TranslationTable*)newTT->makeDynamic();

}
//_______________________________________________________________________
void    _DataSet::Finalize (void)
{
    if (streamThrough) {
        fclose (streamThrough);
        streamThrough = nil;
        theMap.Clear();
    } else {
        if (useHorizontalRep) {
            bool  good = true;
            for (long s = 0; s<lLength; s++) {
                ((_String*)lData[s])->Finalize();
                good = good && ((_String*)lData[0])->sLength == ((_String*)lData[s])->sLength;
            }

            if (!good) {
                Clear();
                WarnError ("Internal Error in _DataSet::Finalize. Unequal sequence lengths in compact representation");
                return;
            }

            _List               dups;
            _List               uniquePats;
            _AVLListX           dupsAVL (&dups);

            long  siteCounter = ((_String*)lData[0])->sLength;

            for (long i1 = 0; i1<siteCounter; i1++) {
                _Site * tC = new _Site ();
                checkPointer (tC);

                for (long i2 = 0; i2 < lLength; i2++) {
                    (*tC) << ((_String*)lData[i2])->sData[i1];
                }

                tC->Finalize();

                long ff = dupsAVL.Find(tC);
                if (ff<0) {
                    uniquePats << tC;
                    dupsAVL.Insert(tC,theFrequencies.lLength);
                    theMap         << theFrequencies.lLength;
                    theFrequencies << 1;
                } else {
                    ff =  dupsAVL.GetXtra(ff);
                    theMap << ff;
                    theFrequencies.lData[ff] ++;
                }

                DeleteObject (tC);
            }
            dupsAVL.Clear(false);
            _List::Clear();
            _List::Duplicate(&uniquePats);
        } else {
            long j,
                 k;

            _Site         *tC;
            {
                _List          dups;
                _AVLListX      dupsAVL (&dups);

                for (long i1 = 0; i1<lLength; i1++) {
                    tC = (_Site*)lData[i1];
                    long ff = dupsAVL.Find(tC);
                    if (ff<0) {
                        dupsAVL.Insert(tC,i1);
                    } else {
                        ff =  dupsAVL.GetXtra(ff);
                        tC->Clear();
                        tC->SetRefNo(ff);
                        theFrequencies.lData[ff] ++;
                    }
                }
                dupsAVL.Clear(false);
            }

            _SimpleList  refs(lLength),
                         toDelete(lLength);
            j = 0;

            for (long i1 = 0; i1<lLength; i1++) {
                tC = (_Site*)(*(_List*)this)(i1);
                k = tC->GetRefNo();
                if (k==-1) {
                    refs << j++;
                } else {
                    toDelete << i1;
                    refs << -1;
                }
            }

            for (long i2=0; i2<lLength; i2++) {
                tC = (_Site*)(*(_List*)this)(i2);
                k = tC->GetRefNo();
                if (k>=0) {
                    j = refs.lData[k];
                    if (j<0) {
                        warnError (-171);
                    } else {
                        refs.lData[i2]=j;
                    }
                }
            }

            theMap.Clear();
            theMap.Duplicate (&refs);
            DeleteList (toDelete);
            theFrequencies.DeleteList (toDelete);

            for (long i3 = 0; i3<lLength; i3++) {
                tC = (_Site*)(*(_List*)this)(i3);
                tC->SetRefNo (0);
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
//_______________________________________________________________________
void    _DataSet::Compact (long index)
{
    if (useHorizontalRep) {
        WarnError ("Internal Error: _DataSet::Compact called with compact represntation");
        return;
    }

    _Site* tC = (_Site*)(*(_List*)this)(index);
    if (tC->GetRefNo()!=-1)
        // take care of double referencing
    {
        _Site*tCC=tC;
        long lastRef, count = 0;
        do {
            lastRef = tCC->GetRefNo();
            count++;
            tCC = (_Site*)(*(_List*)this)(tCC->GetRefNo());
        } while (tCC->GetRefNo()!=-1);
        if (count>1) {
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
inline char  _DataSet::operator () (unsigned long site, unsigned long pos, unsigned int)
{
    return (((_String**)lData)[theMap.lData[site]])->sData[pos];
}


//_________________________________________________________
long _DataSet::ComputeSize(void)
{
    long res = sizeof (_DataSet);

    res+=(theMap.lLength+lLength+theFrequencies.lLength)*sizeof(long);
    res+=lLength*sizeof(_Site);

    for (long i=0; i<lLength; i++) {
        res+= ((_Site*)(*(_List*)this)(i))->sLength;
    }

    return res;
}

//_________________________________________________________
_Parameter _DataSet::CheckAlphabetConsistency(void)
{
    long        charsIn = 0,
                gaps    = 0,
                total   = 0;

    char        checks    [256],
                gapChar = theTT->GetGapChar();

    _String     baseSymbols;

    if (theTT->baseSet.sLength) {
        baseSymbols = theTT->baseSet;
    } else if (theTT->baseLength == 4) {
        baseSymbols = "ACGUT";
    } else if (theTT->baseLength == 20) {
        baseSymbols = "ACDEFGHIKLMNOPQRSTVWY";
    } else {
        baseSymbols = binaryOneCharCodes;
    }



    for (; charsIn<256; charsIn++) {
        checks[charsIn] = 0;
    }

    for (charsIn=0; charsIn<baseSymbols.sLength; charsIn++) {
        checks[baseSymbols.sData[charsIn]] = 1;
    }

    charsIn = 0;

    for (long i = 0; i<lLength; i++) {
        _String* thisColumn = (_String*)lData[i];
        long     w = theFrequencies.lData[i];
        for (long j = 0; j<thisColumn->sLength; j++)
            if (checks[thisColumn->sData[j]]) {
                charsIn+=w;
            } else if (gapChar == thisColumn->sData[j]) {
                gaps += w;
            }

        total += w*thisColumn->sLength;
    }

    return (_Parameter)charsIn/(total-gaps+1.);

}

//___________________________________________________

BaseRef _DataSet::toStr (void)
{
    _String * s = new _String(NoOfSpecies()*30, true),
    *str;

    checkPointer (s);
    (*s) << _String ((long)NoOfSpecies());
    (*s) << " species:";

    str = (_String*)theNames.toStr();
    (*s) << *str;
    DeleteObject(str);

    (*s) << ";\nTotal Sites:";
    (*s) << _String((long)GetNoTypes());
    (*s) << ";\nDistinct Sites:";
    (*s) << _String((long)theFrequencies.lLength);

    s->Finalize();

    return s;
}

//___________________________________________________

void    _DataSet::toFileStr (FILE* dest)
{
    fprintf (dest, "%ld species: ",NoOfSpecies());
    theNames.toFileStr(dest);

    fprintf (dest, ";\nTotal Sites: %ld",GetNoTypes()) ;
    fprintf (dest, ";\nDistinct Sites: %ld",theFrequencies.lLength);

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
void    _DataSet::constructFreq (long* d, _Matrix*m, char positions, long column, long counter, int level,
                                 int shifter, int index)
{
    for(int i=0; i<theTT->baseLength; i++) {
        if (d[level*theTT->baseLength+i]) {
            if (level) {
                constructFreq (d,m,positions,column,counter, level-1,shifter*theTT->baseLength,index + i*shifter);
            } else {
                m->theData[(index + i*shifter)*positions+column]+=1.0/counter;
            }
        }
    }
}

//_________________________________________________________
void    _DataSet::constructFreq (long* d, _Parameter *m, char positions, long column, long counter, int level,
                                 int shifter, int index)
{
    for(int i=0; i<theTT->baseLength; i++)
        if (d[level*theTT->baseLength+i]) {
            if (level) {
                constructFreq (d,m,positions,column,counter, level-1,shifter*theTT->baseLength,index + i*shifter);
            } else {
                m[(index + i*shifter)*positions+column]+=1.0/counter;
            }
        }
}

//_________________________________________________________

_Matrix * _DataSet::HarvestFrequencies (char unit, char atom, bool posSpec, _SimpleList& hSegmentation, _SimpleList& vSegmentation, bool countGaps)
{
    long    vD,
            hD = 1,
            i,
            j;

    if (hSegmentation.lLength==0||vSegmentation.lLength<unit) { // default segmentation
        if (hSegmentation.lLength==0) {
            j = NoOfSpecies();
            hSegmentation.RequestSpace (j);
            for (long i = 0; i<j; i++) {
                hSegmentation<<i;
            }
        }
        if (vSegmentation.lLength<unit) {
            vSegmentation.Clear();
            j = GetNoTypes();
            vSegmentation.RequestSpace (j);
            for (long i = 0; i<j; i++) {
                vSegmentation<<i;
            }
        }
    }

    if (unit%atom!=0) {
        ReportWarning (_String ("Atom should divide Unit in HarvestFrequencies call. Bailing out by setting atom = 1"));
        atom = 1; // default bailout
    }

    // create the output Matrix


    for (i=0; i<atom; i++) {
        hD*=theTT->baseLength;
    }

    vD = posSpec?unit/atom:1;

    _Matrix     out (hD, vD, false, true);

    long     positions  =   unit/atom,
             *store        = new long[atom*theTT->baseLength];

    checkPointer(store);

    for (i = 0; i<vSegmentation.lLength; i+=unit) { // loop over the set of segments
        // make sure the partition is kosher

        if (i+unit>vSegmentation.lLength) {
            break;
        }

        for (long jj=i; jj<i+unit; jj+=atom) {
            long   k = (jj-i)/atom,
                   count,
                   m,
                   ll;

            for (ll = 0; ll<hSegmentation.lLength; ll++)
                // loop down each column
            {
                int l = hSegmentation.lData[ll];
                count = 1;
                // build atomic probabilities
                for (m = 0; m<atom; m++ ) {
                    theTT->TokenCode ((*this)(vSegmentation.lData[jj+m],l,atom), store+theTT->baseLength*m,countGaps);
                }

                long index = 0, shifter = 1;
                for (int m = atom-1; m>=0; m--) {
                    int smcount = 0;
                    for (int n = 0; n<theTT->baseLength; n++) {
                        if (store[theTT->baseLength*m+n]) {
                            index += shifter*n;
                            smcount++;
                        }
                    }
                    shifter*=theTT->baseLength;
                    count *=smcount;
                }

                if (count>1) {
                    constructFreq (store, &out, posSpec?positions:1, posSpec?k:0, count, atom-1 , 1, 0);
                } else {
                    out.theData[posSpec?index*positions+k:index] += count;
                }
            }
        }
    }

    delete[] store;
    //scale the matrix now

    hD = out.GetHDim();
    vD = out.GetVDim();
    for (i=0; i<vD; i++) {
        _Parameter temp = 0.0;

        for (long r=hD-1; r>=0; r--) {
            temp+=out.theData[r*vD+i];
        }

        {
            for (long r=i; r<vD*hD; r+=posSpec?positions:1) {
                out.theData[r]/=temp;
            }
        }
    }


    return (_Matrix*)out.makeDynamic();
}

//_________________________________________________________

void    _DataSet::AddName (_String& s)
{
    s.Trim(0,s.FirstNonSpaceIndex (0,-1,-1));
    theNames&&(&s);
}



//_________________________________________________________

void    _DataSet::MatchIndices (_Formula&f, _SimpleList& receptacle, bool isVert, long limit)
{
    _String     varName  = isVert ? "siteIndex" : "speciesIndex";
    _Variable   *v       = CheckReceptacle (&varName, empty, false);

    for (long i=0; i<limit; i++) {
        v->SetValue (new _Constant(i), nil);
        _PMathObj res = f.Compute();
        if (res && !CheckEqual(res->Value(),0.0)) {
            receptacle<<i;
        }
    }
    v->SetValue (new _Constant(0.0), nil);
}

//_________________________________________________________

void    _DataSet::FindAllSitesLikeThisOne (long index, _SimpleList& receptacle)
{
    if (index>=0 && index<theMap.lLength) {
        index = theMap.lData[index];
        for (long k=0; k<theMap.lLength; k++)
            if (theMap.lData[k]==index) {
                receptacle << k;
            }
    }
}

//_________________________________________________________

_TranslationTable*      _DataSet::CheckCompatibility (_SimpleList& ref, char concatOrCombine)
{
    _DataSet* currentSet = (_DataSet*)dataSetList(ref(0));

    _TranslationTable* theEnd = new _TranslationTable (*(currentSet->theTT));
    checkPointer(theEnd);
    long    refNo     =  concatOrCombine?currentSet->NoOfSpecies():currentSet->NoOfColumns();
    char    emptyChar = theEnd->GetSkipChar();

    for (long k=1; k<ref.lLength; k++) {
        currentSet = (_DataSet*)dataSetList(ref(k));

        _TranslationTable* tryMe = theEnd->MergeTables (currentSet->theTT);

        if (tryMe) {
            if (emptyChar) {
                DeleteObject (theEnd);
                theEnd = tryMe;
                continue;
            } else {
                if ((concatOrCombine&&(currentSet->NoOfSpecies()==refNo))||(!concatOrCombine&&(currentSet->NoOfColumns()==refNo))) {
                    DeleteObject (theEnd);
                    theEnd = tryMe;
                    continue;
                }
            }
        }
        _String warningMessage ("The data set:");
        warningMessage = warningMessage & *((_String*)dataSetNamesList(ref(k))) & _String (" was found incompatible with one of the following data sets:");
        for (long i=0; i<k; i++) {
            warningMessage = warningMessage & *((_String*)dataSetNamesList(ref(i))) & _String (",");
        }
        warningMessage = warningMessage &  _String (" and was dropped from the dataset merging operation");
        ReportWarning (warningMessage);
        ref.Delete (k);
        k--;
    }

    return theEnd;
}

//_________________________________________________________

_DataSet*   _DataSet::Concatenate (_SimpleList ref)

// concatenates (adds columns together) several datasets
// in case the number of species in the datasets are different the deficiencies will be padded
// by omission symbols
// in case translation tables are different, they will be merged, provided it can be done,
// otherwise the incompatible datasets will be ignored during this operation.

{
    _TranslationTable  * jointTable;

    jointTable = CheckCompatibility (ref,1);

    _DataSet           * bigDataSet = new _DataSet;
    checkPointer(bigDataSet);

    bigDataSet->theTT = jointTable;

    // pass one - determine the max max number of species present and what dataset are they coming from

    long      maxSpecies=0,
              maxDataSet=0,
              siteIndex;

    _DataSet *currentSet;

    char     emptySlot = jointTable->GetSkipChar();

    for (long i=0; i<ref.lLength; i++) {
        currentSet = (_DataSet*)dataSetList(ref(i));

        long       specCount = currentSet->NoOfSpecies(),
                   siteCount = currentSet->NoOfColumns();


        if (specCount>maxSpecies) {
            maxSpecies = specCount;
            maxDataSet = i;
        }
        for (long j=0; j<siteCount; j++) {
            bigDataSet->AddSite((*currentSet)(j,0,1));
        }
    }

    for (long k=1; k<maxSpecies; k++) {
        siteIndex = 0;
        for (long i=0; i<ref.lLength; i++) {
            currentSet = (_DataSet*)dataSetList(ref.lData[i]);

            long       cns = currentSet->NoOfSpecies(),
                       cnc = currentSet->NoOfColumns();

            if (cns<=k)
                for (long j=0; j< cnc; j++, siteIndex++) {
                    bigDataSet->Write2Site(siteIndex,emptySlot);
                }
            else
                for (long j=0; j< cnc; j++, siteIndex++) {
                    bigDataSet->Write2Site(siteIndex,(*currentSet)(j,k,1));
                }
        }
    }

    currentSet = (_DataSet*)dataSetList(ref(maxDataSet));
    {
        for (long i=0; i<maxSpecies; i++) {
            bigDataSet->AddName (*((_String*)(currentSet->GetNames())(i)));
        }
    }

    bigDataSet->Finalize();
    bigDataSet->SetNoSpecies (maxSpecies);
    return bigDataSet;
}

//_________________________________________________________

_DataSet*   _DataSet::Combine (_SimpleList ref)

// combines (adds rows together) several datasets
// in case the number of species in the datasets are different the deficiencies will be padded
// by omission symbols
// in case translation tables are different, they will be merged, provided it can be done,
// otherwise the incompatible datasets will be ignored during this operation.

{
    _TranslationTable  * jointTable;

    jointTable = CheckCompatibility (ref,0);

    _DataSet           * bigDataSet = new _DataSet;
    checkPointer(bigDataSet);
    bigDataSet->theTT = jointTable;

    // pass one - determine the max max number of sites present and what dataset are they coming from

    long     i,
             j,
             k,
             maxSites=0,
             sitesAvail,
             nsc = 0;

    _DataSet *currentSet;

    char     emptySlot = jointTable->GetSkipChar();

    for (i=0; i<ref.lLength; i++) {
        currentSet = (_DataSet*)dataSetList(ref(i));
        if (currentSet->NoOfColumns()>maxSites) {
            maxSites = currentSet->NoOfColumns();
        }
        nsc += currentSet->NoOfSpecies();
    }


    for (k=0; k<ref.lLength; k++) {
        currentSet = (_DataSet*)dataSetList(ref(k));
        sitesAvail = currentSet->NoOfColumns();

        long     cns = currentSet->NoOfSpecies();
        for (i=0; i<cns; i++) {
            bigDataSet->AddName (*((_String*)(currentSet->GetNames())(i)));
            if (!(k||i)) {
                for (j=0; j<sitesAvail; j++) {
                    bigDataSet->AddSite ((*currentSet)(j,0,1));
                }
                for (; j<maxSites; j++) {
                    bigDataSet->AddSite (emptySlot);
                }
            } else {
                for (j=0; j<sitesAvail; j++) {
                    bigDataSet->Write2Site (j,(*currentSet)(j,i,1));
                }
                for (; j<maxSites; j++) {
                    bigDataSet->Write2Site (j,emptySlot);
                }
            }
        }
    }

    bigDataSet->Finalize();
    bigDataSet->SetNoSpecies(nsc);
    return bigDataSet;
}

//_________________________________________________________
// Data Set Filter/Numeric
//_________________________________________________________

_DataSetFilter::_DataSetFilter (void)
{
    unitLength = 0;
    theData = NULL;
    accessCache = nil;
}
//_________________________________________________________
_DataSetFilter::_DataSetFilter (_DataSet* ds, char, _String&)
{
    theData     = ds;
    accessCache = nil;
}
//_________________________________________________________
_DataSetFilter::~_DataSetFilter (void)
{
    if (accessCache) {
        DeleteObject (accessCache);
    }
}

//_______________________________________________________________________

_DataSetFilterNumeric::_DataSetFilterNumeric (_Matrix* freqs, _List& values, _DataSet *ds, long cc)
{
    unitLength      = 1;
    categoryCount   = cc;

    SetData         (ds);

    _SimpleList     baseFreqs;

    freqs->ConvertToSimpleList (baseFreqs);
    dimension =                ((_Matrix*)values(0))->GetVDim();

    theNodeMap.Populate         (ds->GetNames().lLength,0,1);
    theOriginalOrder.Populate   (((_Matrix*)values(0))->GetHDim()/categoryCount,0,1);

    //theMap.Populate           (theFrequencies.lLength,0,1);
    //theOriginalOrder.Populate     (theFrequencies.lLength,0,1);
    //duplicateMap.Populate     (theFrequencies.lLength,0,1);


    /*CreateMatrix (&probabilityVectors, theNodeMap.lLength, shifter,false,true, false);

    _Parameter   *storeHere = probabilityVectors.theData;
    for (long spec = 0; spec < theNodeMap.lLength; spec++)
    {
        _Matrix * specMatrix = (_Matrix*)values(spec);
        for (long site = 0; site < theFrequencies.lLength; site++)
            for (long state = 0; state < dimension; state++,storeHere++)
                //probabilityVectors.theData [shifter*spec + site*dimension+state]
                *storeHere= specMatrix->theData[site*dimension+state];
    }*/

    _List        siteScores;
    _AVLListXL   siteIndices(&siteScores);

    duplicateMap.RequestSpace (baseFreqs.lLength+1);

    //bool       startD = false;

    char buffer[255];

    for (long site =0; site <baseFreqs.lLength; site++) {
        _Parameter      testV = 0.0;

        for (long k=0; k<theNodeMap.lLength; k++) // sweep down the columns
            for (long state = 0; state < dimension; state++) {
                testV += ((_Matrix*)(((_Matrix**)values.lData)[k]))->theData[site*dimension+state];
            }

        sprintf     (buffer, "%20.18g", testV);
        _String     testS (buffer);
        long        f = siteIndices.Find (&testS);

        _SimpleList * sameScore = nil;

        if (f>=0) {
            sameScore = (_SimpleList*)siteIndices.GetXtra (f);
            for (long k = 0; k<sameScore->lLength; k++) {
                bool fit = true;
                f        = sameScore->lData[k];


                for (long spec=0; spec<theNodeMap.lLength && fit; spec++) { // sweep down the columns
                    _Matrix * specMatrix =(_Matrix*)(((_Matrix**)values.lData)[spec]);
                    for (long state = 0; state < dimension; state++)
                        if (specMatrix->theData[site*dimension+state]!=specMatrix->theData[theMap.lData[f]*dimension+state]) {
                            fit = false;
                            break;
                        }
                }

                if (fit) {
                    theFrequencies[f]+=baseFreqs[site];
                    duplicateMap<<f;
                    f = 0;
                    break;
                } else {
                    f = -1;
                }
            }
        }
        if (f==-1) {
            if (!sameScore) {
                checkPointer        (sameScore = new _SimpleList);
                if (siteIndices.Insert  (testS.makeDynamic(),(long)sameScore,false) < 0) {
                    _String eh ("WTF?");
                    StringToConsole(eh);
                }
            }

            (*sameScore) << theFrequencies.lLength;

            duplicateMap<<theFrequencies.lLength;
            theFrequencies<<baseFreqs[site];
            theMap<<site;
        }
    }

    siteIndices.Clear(true);
    shifter         = theFrequencies.lLength*dimension;
    categoryShifter = shifter*theNodeMap.lLength;

    CreateMatrix (&probabilityVectors, theNodeMap.lLength, shifter*categoryCount,false,true, false);
    _Parameter   *storeHere    = probabilityVectors.theData;

    long      refShifter = 0;
    for (long cc = 0; cc < categoryCount; cc++, refShifter += theOriginalOrder.lLength * dimension) {
        for (long spec = 0; spec < theNodeMap.lLength; spec++) {
            _Matrix * specMatrix = (_Matrix*)values(spec);
            for (long site = 0; site < theFrequencies.lLength; site++)
                for (long state = 0; state < dimension; state++,storeHere++) {
                    *storeHere = specMatrix->theData[refShifter + theMap.lData[site]*dimension+state];
                }
        }
    }
}


//_______________________________________________________________________

void _DataSetFilter::CopyFilter (_DataSetFilter *copyFrom)
{
    memcpy ((char*)this, (char*)copyFrom, sizeof (_DataSetFilter));

    theFrequencies.Duplicate        (&copyFrom->theFrequencies);
    theNodeMap.Duplicate            (&copyFrom->theNodeMap);
    theMap.Duplicate                (&copyFrom->theMap);
    theOriginalOrder.Duplicate      (&copyFrom->theOriginalOrder);
    conversionCache.Duplicate       (&copyFrom->conversionCache);
    duplicateMap.Duplicate          (&copyFrom->duplicateMap);

    nInstances              = 1;
    dimension               = copyFrom->dimension;
    undimension             = copyFrom->undimension;
    unitLength              = copyFrom->unitLength;
    accessCache             = nil;

}


//_______________________________________________________________________

BaseRef _DataSetFilter::makeDynamic (void)
{
    _DataSetFilter * r = new _DataSetFilter;
    checkPointer    (r);
    r->CopyFilter   (this);

    return r;
}


//_______________________________________________________________________

BaseRef _DataSetFilterNumeric::makeDynamic (void)
{
    _DataSetFilterNumeric * r = new _DataSetFilterNumeric();
    checkPointer            (r);
    r->CopyFilter           (this);
    r->probabilityVectors.Duplicate(&probabilityVectors);
    return r;
}

//_______________________________________________________________________

_Parameter * _DataSetFilterNumeric::getProbabilityVector (long spec, long site, long categoryID)
{
    return probabilityVectors.theData + categoryID * categoryShifter + spec * shifter + site * dimension;
}

//_______________________________________________________________________
long    _DataSetFilter::FreeUpMemory (long requestedBytes)
{
    long res = 0;
    for (long i=0; (i<theMap.lLength)&&(res<requestedBytes); i++) {
        res+=(theData->GetSite(theMap[i]))->FreeUpMemory(requestedBytes-res);
    }
    return res;
}


//_______________________________________________________________________
void    _DataSetFilter::SetDimensions (void)
{
    dimension   = GetDimension(true);
    undimension = GetDimension(false);
}

//_______________________________________________________________________
unsigned long    _DataSetFilter::FindUniqueSequences  (_SimpleList& indices, _SimpleList& map, _SimpleList& counts, short mode)
{
    indices.Clear(); map.Clear(); counts.Clear();
    
    unsigned long             sites  = theMap.lLength,
    seqs   = theNodeMap.lLength,
    unit   = GetUnitLength();
    
    if (mode == 0)
        {
        _SimpleList hashSupport;        
        _AVLListXL  sequenceHashes     (&hashSupport);
        
        for (unsigned long sequenceIndex = 0; sequenceIndex < seqs; sequenceIndex ++){
            _String * thisSequence = GetSequenceCharacters (sequenceIndex);
            
            long     sequenceHash   = thisSequence->Adler32(),
            f              = sequenceHashes.Find ((BaseRef)sequenceHash),
            rawSequenceIdx = theNodeMap.lData[sequenceIndex];
            
            DeleteObject (thisSequence);
            
            _SimpleList * sameScore = nil;
            if (f>=0) {
                sameScore = (_SimpleList*)sequenceHashes.GetXtra (f);
                for (long k = 0; k<sameScore->lLength; k++) {
                    bool fit = true;
                    f = sameScore->lData[k];
                    
                    long fRaw = theNodeMap.lData[indices.lData[f]];
                    
                    for (unsigned long site = 0; site < sites && fit; site++){
                        for (unsigned long unitIndex = 0; unitIndex < unit; unitIndex ++){
                            _Site * thisSite = theData->GetSite(theMap.lData[unit*site+unitIndex]);
                            if (thisSite->sData[fRaw] != thisSite->sData[rawSequenceIdx]){
                                fit = false;
                                break;
                            }
                        }
                    }
                    
                    if (fit) {
                        map << f;
                        counts.lData[f] ++;
                        
                    } else {
                        f = -1;
                    }
                }
            }
            if (f==-1) { // fit failed or unique site
                if (!sameScore) {
                    sameScore = (_SimpleList*)checkPointer(new _SimpleList);
                    sequenceHashes.Insert ((BaseRef)sequenceHash,(long)sameScore,false);
                }
                
                (*sameScore) << indices.lLength;
                map     << indices.lLength;
                indices << sequenceIndex;
                counts  << 1;
            }
        }
        
        }
    else{
        long             vd  = GetDimension(true);
        
        _Parameter      *translatedVector = (_Parameter*)checkPointer(new _Parameter [vd]),
        *translatedVector2= (_Parameter*)checkPointer(new _Parameter [vd]);
        
        _String         state1 (unit,false),
        state2 (unit,false);
        
        for (long sequenceIndex = 0; sequenceIndex < seqs; sequenceIndex++) {
            bool checkState = false;
            for (long idx=0; idx<indices.countitems(); idx++) {
                for (long m=0; m<sites; m++) {
                    RetrieveState (m,sequenceIndex, state1);
                    RetrieveState (m,indices.lData[idx], state2);
                    
                    checkState = true;
                    long idx1 = Translate2Frequencies (state1, translatedVector,  true),
                    idx2 = Translate2Frequencies (state2, translatedVector2, true);
                    
                    printf ("(%ld, %ld) %ld = %ld %ld\n", sequenceIndex, indices.lData[idx], m, idx1, idx2); 
                    
                    if (idx2 >=0 && idx1 >=0) {
                        if (idx1==idx2) {
                            continue;
                        } else {
                            checkState = false;
                            break;
                        }
                    } else {
                        
                        // check for equal ambigs 
                        long k = 0;
                        for (; k < vd; k++){
                            if (translatedVector[k] != translatedVector2[k]){
                                checkState = false;
                                break;
                            }
                        
                        if (checkState){
                            
                                long count1 = 0,
                                     count2 = 0;
                                
                                for (long t = 0; t<vd; t++) {
                                    count1 += translatedVector[t]>0.0;
                                    count2 += translatedVector2[t]>0.0;
                                }
                                
                                if (count1 < vd && count2 < vd) {
                                    checkState = false;
                                    break;
                                }
                                
                            } else {
                                bool first  = mode==2,
                                second = mode==2;
                                if (mode == 2){
                                    for (long t = 0; (first||second)&&(t<vd); t++) {
                                        if (translatedVector[t]>0.0) {
                                            second &= (translatedVector2[t]>0.0);
                                        }
                                        if (translatedVector2[t]>0.0) {
                                            first  &= (translatedVector[t]>0.0);
                                        }
                                    }
                                    if (!(first||second)) {
                                        checkState = false;
                                        break;
                                    }
                                } else {
                                    for (long t = 0; (first||second)&&(t<vd); t++) {
                                        if (translatedVector[t]>0.0) {
                                            second |= (translatedVector2[t]>0.0);
                                        }
                                        if (translatedVector2[t]>0.0) {
                                            first  |= (translatedVector[t]>0.0);
                                        }
                                    }                       
                                    if (!(first&&second)) {
                                        checkState = false;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                
                if (checkState) {
                    map << idx;
                    counts.lData[idx] ++;
                    break;
                }
            }
            
            if (!checkState){
                map     << indices.lLength;
                indices << sequenceIndex;
                counts  << 1;
            }
            
        }
        
        delete [] translatedVector;
        delete [] translatedVector2;   
    }
    
    
    return indices.lLength;
}

//_______________________________________________________________________

void    _DataSetFilter::SetFilter (_DataSet* ds, char unit, _SimpleList& horizontalList, _SimpleList& verticalList, bool isFilteredAlready)
{
    // we must step thru the underlying dataset and recompute the frequenices
    // we will store the vertical map in theMap
    // and the horizontal map in theNodeMap
    // theFrequencies will be store the new frequencies
    // theOriginalOrder is the receptacle for the original site order in the data filter

    bool            copiedSelf = false; // tag if refiltering self

    _DataSetFilter* firstOne = nil;
    if (isFilteredAlready) {
        if ((Ptr)this == (Ptr)ds) {
            firstOne = (_DataSetFilter*)makeDynamic();
            copiedSelf = true;
        } else {
            firstOne = (_DataSetFilter*)ds;
        }
        ds       = firstOne->theData;
    }

    theMap.Clear();
    theNodeMap.Clear();
    theOriginalOrder.Clear();
    theFrequencies.Clear();
    theExclusions.Clear();
    conversionCache.Clear();
    duplicateMap.Clear();

    theData     = ds;
    unitLength  = unit;

    long        i,
                j;

    // security checks
    if (!horizontalList.lLength||(verticalList.lLength<unit)) {
        ReportWarning (_String("Row and/or column partition is empty. All the data will be used by default."));
        if (horizontalList.lLength==0) {
            if (!isFilteredAlready) {
                j = ds->NoOfSpecies();
            } else {
                j = firstOne->theNodeMap.lLength;
            }

            horizontalList.Populate (j,0,1);
        }
        if (verticalList.lLength<unit) {
            verticalList.Clear();
            if (!isFilteredAlready) {
                j = ds->GetNoTypes();
            } else {
                j = firstOne->theOriginalOrder.lLength;
            }

            verticalList.Populate (j,0,1);
        }
    }

    if (!isFilteredAlready) {
        theNodeMap.Clear();
        theNodeMap.Duplicate (&horizontalList);
    } else {
        for (long k = 0; k<horizontalList.lLength; k++) {
            theNodeMap<<firstOne->theNodeMap.lData[horizontalList.lData[k]];
        }

        horizontalList.Clear();
        horizontalList.Duplicate(&verticalList);
        verticalList.Clear();
        verticalList.RequestSpace(firstOne->theOriginalOrder.lLength);

        for (i = 0; i<horizontalList.lLength; i++) {
            j = horizontalList.lData[i];
            if (j>=0 && j<firstOne->theOriginalOrder.lLength) {
                verticalList<<firstOne->theOriginalOrder.lData[j];
            } else {
                _String tooBig (j);
                if (j<0) {
                    ReportWarning  (tooBig &" is a negative site index and is ignored");
                } else {
                    ReportWarning  (tooBig &" exceeds the number of sites in the underlying data filter and is ignored");
                }
            }
        }
    }

    j = ds->NoOfSpecies();

    for (i=0; i<theNodeMap.lLength; i++) {
        if (theNodeMap.lData[i]>=j) {
            _String invalid(theNodeMap.lData[i]);
            ReportWarning ((invalid&" exceeds the number of species in the underlying dataset and is ignored"));
            theNodeMap.Delete(i);
            i--;
        }
    }

    j = ds->GetNoTypes();
    for (i=0; i<verticalList.lLength; i++)
        if (verticalList.lData[i]>=j) {
            _String invalid(verticalList.lData[i]);
            ReportWarning ((invalid&" exceeds the number of sites in the underlying dataset and is ignored"));
            verticalList.Delete(i);
            i--;
        }

    if (verticalList.lLength%unit!=0) {
        ReportWarning (_String("Number of sites in datasetfilter is not divisible by the unit - will truncate to the nearest integer"));
        while(verticalList.lLength%unit) {
            verticalList.Delete(verticalList.lLength-1);
        }
    }


    theOriginalOrder.Duplicate (&verticalList);


    // done with security checks

    _SimpleList indices;        // numeric indices intended to facilitate the reindexing
    _AVLListXL  siteIndices     (&indices);


    // sweep through the columns left to right

    duplicateMap.RequestSpace (verticalList.lLength/unit+1);

    _String      siteHolder   (unit*theNodeMap.lLength,false);

    //bool       startD = false;

    for (i=0; i<verticalList.lLength; i+=unit) {
        long colIndex = 0;

        for (j=0; j<unit; j++) // sweep within one block
            for (long k=0; k<theNodeMap.lLength; k++) // sweep down the columns
                //colIndex+=
                //(((_String**)ds->lData)[ds->theMap.lData[verticalList.lData[i+j]]])->sData[theNodeMap.lData[k]];
            {
                siteHolder[colIndex++] = (((_String**)ds->lData)[ds->theMap.lData[verticalList.lData[i+j]]])->sData[theNodeMap.lData[k]];
            }

        colIndex = siteHolder.Adler32();

        long        f = siteIndices.Find ((BaseRef)colIndex);
        _SimpleList * sameScore = nil;
        
        if (f>=0) {
            sameScore = (_SimpleList*)siteIndices.GetXtra (f);
            for (long k = 0; k<sameScore->lLength; k++) {
                bool fit = true;
                f = sameScore->lData[k];
                for (long j=0; fit&&(j<unit); j++) { // sweep within one block
                    _Site * site1 = ds->GetSite(verticalList.lData[i+j]),
                          * site2 = ds->GetSite(theMap.lData[unit*f+j]); 
                            
                    for (long k=0; k<theNodeMap.lLength; k++) // sweep down the columns
                        if (site1->sData[theNodeMap.lData[k]]!=site2->sData[theNodeMap.lData[k]]) {
                            fit = false;
                            break;
                        }
                }

                if (fit) {
                    theFrequencies[f]++;
                    duplicateMap<<f;
                    f = 0;
                    break;
                } else {
                    f = -1;
                }
            }
        }
        if (f==-1) { // fit failed or unique site
            if (!sameScore) {
                sameScore = (_SimpleList*)checkPointer(new _SimpleList);
                siteIndices.Insert ((BaseRef)colIndex,(long)sameScore,false);
             }

            (*sameScore) << theFrequencies.lLength;
            duplicateMap<<theFrequencies.lLength;
            theFrequencies<<1;
            for (j=0; j<unit; j++) {
                theMap<<verticalList.lData[i+j];
            }
        }
    }

    siteIndices.Clear();

    duplicateMap.TrimMemory();
    theOriginalOrder.TrimMemory();

    if (copiedSelf) {
        DeleteObject (firstOne);
    }

    SetDimensions();
    FilterDeletions();

}
//_______________________________________________________________________
long    _DataSetFilter::FindSpeciesName (_List& s, _SimpleList& r)
{
    // MOD 12/16/03
    r.Clear();

    _List           newNames;
    _AVLListX       matched (&newNames);

    for (long k=0; k<theNodeMap.lLength; k++) {
        long i = theNodeMap.lData[k];
        _String * uC = new _String (*(_String*)theData->theNames (i));
        uC->UpCase();
        matched.Insert (uC,i);
    }

    for (long m = 0; m < s.lLength; m++) {
        _String ts (*((_String*)s(m)));
        ts.UpCase();
        long f = matched.Find (&ts);
        if (f>=0) {
            r << matched.GetXtra (f);
        } else {
            break;
        }
    }

    return r.lLength;
}

//_______________________________________________________________________
extern _String skipOmissions;

void    _DataSetFilter::FilterDeletions(_SimpleList *theExc)
{
    _Parameter      skipo;
    checkParameter (skipOmissions,skipo,0.0);

    if (skipo>.5 || theExc ) { // delete omissions
        //build up the list of "bad" sites
        _SimpleList sitesWithDeletions;
        if (!theExc) {
            for (long i=0; i<theFrequencies.lLength; i++)
                if (HasDeletions(i)) {
                    sitesWithDeletions<<i;
                }
        } else {
            _Parameter   *store_vec = (_Parameter*)checkPointer(new _Parameter [GetDimension(false)]);

            for (long i=0; i<theFrequencies.lLength; i++) {
                long pos = HasExclusions(i,theExc,store_vec);
                if  (pos != -1) {
                    sitesWithDeletions<<i;
                    _String warnMsg ((*this)(i,pos));
                    warnMsg = warnMsg & " was encountered in sequence "& *GetSequenceName (pos) & " at site pattern " & i
                              & ". All corresponding alignment columns will be removed from subsequent analyses.";
                    ReportWarning (warnMsg);
                }
            }

            delete [] store_vec;
        }

        if (sitesWithDeletions.lLength==theFrequencies.lLength) {
            _String errMsg ("All the sites in the datafilter have deletions and removing them creates an empty filter");
            ReportWarning(errMsg);
        }

        _SimpleList allDeleted,
                    dupDeletes;

        for (long k=0; k < duplicateMap.lLength; k++)
            if (sitesWithDeletions.BinaryFind (duplicateMap.lData[k]) >= 0) {
                dupDeletes << k;
                for (long j = 0; j < unitLength; j++ ) {
                    allDeleted << k*unitLength + j;
                }
            }

        duplicateMap.DeleteList (dupDeletes);
        dupDeletes.Clear();
        theOriginalOrder.DeleteList (allDeleted);
        theFrequencies.DeleteList (sitesWithDeletions);


        for (long i=0; i<sitesWithDeletions.lLength; i++) {
            long sitePos = sitesWithDeletions.lData[i];

            for (long j=0; j<unitLength; j++) {
                theMap.lData[sitePos*unitLength+j]=-1;
                dupDeletes << sitePos*unitLength+j;
            }
        }


        if (allDeleted.lLength) {
            /*allDeleted.Sort();*/

            _String     warnMsg ("The following sites are being omitted:"),
                        *s = (_String*)allDeleted.toStr();

            if (!theExc) {
                warnMsg = warnMsg & "(b/c of deletions/omissions)";
            }

            warnMsg = warnMsg&*s;
            DeleteObject(s);
            ReportWarning(warnMsg);

            _SimpleList shiftIdxBy (sitesWithDeletions.lLength+theFrequencies.lLength);

            long        shiftBy = sitesWithDeletions.lLength,
                        marker  = sitesWithDeletions.lData[sitesWithDeletions.lLength-1],
                        markerI = sitesWithDeletions.lLength-2;

            shiftIdxBy.lLength = sitesWithDeletions.lLength+theFrequencies.lLength;

            for (long i=shiftIdxBy.lLength-1; i>=0; i--) {
                if (i==marker) {
                    shiftBy--;
                    if (markerI>=0) {
                        marker = sitesWithDeletions.lData[markerI];
                        markerI --;
                    } else {
                        marker = -1;
                    }
                }
                shiftIdxBy.lData[i] = shiftBy;
            }
            {
                for (long i=0; i<duplicateMap.lLength; i++) {
                    duplicateMap.lData[i] -= shiftIdxBy.lData[duplicateMap.lData[i]];
                }
            }
        }

        // one final pass on theMap to clear it out
        /*for (long i=theMap.lLength-1;i>=0;i--)
            if (theMap(i)<0)
                theMap.Delete(i);*/
        _SimpleList saveMap (theMap);
        theMap.DeleteList (dupDeletes);
        {
            for (long k=0; k<theMap.lLength; k++)
                if (theMap.lData[k] < 0) {
                    saveMap.DeleteList (dupDeletes);
                    WarnError ("Internal Error in _DataSetFilter::FilterDeletions");
                }
        }
    }

}
//_______________________________________________________________________
_DataSetFilter*  _DataSetFilter::PairFilter (long index1, long index2, _DataSetFilter* result)
{
    _SimpleList species;
    species<<theNodeMap(index1);
    species<<theNodeMap(index2);
    result->SetFilter (theData,unitLength,species,theMap);
    if (theExclusions.lLength) {
        _String* s = (_String*)theExclusions.toStr();
        *s = s->Cut (1,s->Length()-2);
        result->SetExclusions   (s);
        DeleteObject(s);
    }
    return result;
}

//_________________________________________________________

void    _DataSetFilter::MatchStartNEnd (_SimpleList& order, _SimpleList& positions, _SimpleList* parent)
{
    if (order.lLength == 0) {
        return;
    }

    long p0 = order.lData[0];

    _Parameter uth;
    checkParameter (useTraversalHeuristic,uth,1.0);

    if (uth>.5) {
        if (parent)
            for (long i = 1; i < order.lLength; i++) {
                unsigned long
                j       = 0,
                n       = theNodeMap.lLength-1,
                p0  = parent->lData[i],
                p1  = order.lData[i];

                while (CompareTwoSites(p0,p1,j)) {
                    j++;
                }
                while (CompareTwoSites(p0,p1,n)) {
                    n--;
                }
                n = (n<<16) + j;
                positions << n;
            }
        else
            for (long i = 1; i < order.lLength; i++) {
                unsigned long j = 0,
                              n = theNodeMap.lLength-1,
                              p1 = order.lData[i];

                while (CompareTwoSites(p0,p1,j)) {
                    j++;
                }
                while (CompareTwoSites(p0,p1,n)) {
                    n--;
                }
                n = (n<<16) + j;
                positions << n;
                p0 = p1;
            }
    } else
        for (long i = 1; i < order.lLength; i++) {
            unsigned long j = 0,
                          n = theNodeMap.lLength-1;

            n = (n<<16) + j;
            positions << n;
        }

}

//_______________________________________________________________________

void    _DataSetFilter::SetExclusions (_String* theList, bool filter)
{

    theExclusions.Clear();
    theList->StripQuotes();

    if (theList->sLength == 0) {
        return;
    }

    _List        *tokens = theList->Tokenize(',');
    _SimpleList  holder;
    _AVLList     exclusions (&holder);

    for (long k = 0; k < tokens->lLength; k++) {
        long posMarker = MapStringToCharIndex(*(_String*)((*tokens)(k)));

        if (posMarker < 0) {
            ReportWarning (_String("Exclusion request for '") & *(_String*)((*tokens)(k)) &"' does not represent a unique state and will therefore be ignored.");
        } else {
            if (exclusions.Insert((BaseRef)posMarker) < 0) {
                ReportWarning (_String("Exclusion symbol for '") & *(_String*)((*tokens)(k)) &"' is included more than once.");
            }
        }
    }

    DeleteObject (tokens);
    exclusions.ReorderList();

    if (filter) {
        FilterDeletions (&holder);
    }

    theExclusions<<holder;
}

//_______________________________________________________________________

_String*    _DataSetFilter::GetExclusions (void)
{
    _String * res = new _String (16L, true);
    checkPointer (res);

    if (theExclusions.lLength) {
        for (long k=0; k<theExclusions.lLength-1; k++) {
            (*res) << ConvertCodeToLetters (theExclusions.lData[k], unitLength);
            (*res) << ',';
        }

        (*res) << ConvertCodeToLetters (theExclusions.lData[theExclusions.lLength-1], unitLength);
    }

    res->Finalize();

    return res;
}

//_______________________________________________________________________

long    _DataSetFilter::GetDimension (bool correct)
{
    long result = theData->theTT->baseLength;
    for (long i=1; i<unitLength; i++) {
        result *= theData->theTT->baseLength;
    }
    if (correct) {
        result-=theExclusions.lLength;
    }
    return result;
}

//____________________________________________________________________________________
//  20110610: SLKP, some cleanup and refactoring

void    _DataSet::ProcessPartition (_String& input2 , _SimpleList& target , bool isVertical, _SimpleList* additionalFilter, _SimpleList* otherDimension)
{
    if (!input2.sLength) {
        return;
    }
    // decide if the input is an enumeration or a formula
    long totalLength;

    if (additionalFilter) {
        totalLength = additionalFilter->lLength;
    } else {
        totalLength = isVertical?theMap.lLength:noOfSpecies;
    }

    _String input (input2);

    if (!input.IsALiteralArgument(true)) { // not a literal argument
        _Formula fmla, lhs;

        long     varRef = 0,
                 outcome = Parse (&fmla, input, varRef, nil,&lhs);

        if (outcome!=HY_FORMULA_EXPRESSION) {
            WarnError (input & _String(" is an invalid partition specification"));
            return;
        }
        _PMathObj   fV = fmla.Compute();
        if (fV && fV->ObjectClass()==STRING) {
            _String newSpec (128L, true);
            newSpec << '"';
            newSpec << ((_FString*)fV)->theString;
            newSpec << '"';
            newSpec.Finalize();
            ProcessPartition (newSpec, target, isVertical, additionalFilter);
        } else {
            _DataSet::MatchIndices (fmla, target, isVertical, totalLength);
        }
    } else { // an explicit enumeration or a regular expression
        if (input.getChar(0)=='/' && input.getChar(input.sLength-1)=='/')
            // a regular expression
        {
            input.Trim(1,input.sLength-2);
            int   errCode;
            Ptr   regex = PrepRegExp (&input, errCode, true);
            if (errCode) {
                WarnError(GetRegExpError(errCode));
                return;
            }
            // now set do the matching
            // using only the sites that are specced in the additionalFilter

            if (!isVertical) {
                _SimpleList*        eligibleSeqs;

                if (additionalFilter) {
                    eligibleSeqs = additionalFilter;
                } else {
                    eligibleSeqs = new _SimpleList (0, totalLength, 1);
                }

                _SimpleList matches;
                for (long specCount = 0; specCount < eligibleSeqs->lLength; specCount++) {
                    _String pattern (theMap.lLength, false);
                    long    seqPos = eligibleSeqs->lData[specCount];

                    if (otherDimension)
                        for (long seqSlider = 0; seqSlider < otherDimension->lLength; seqSlider ++) {
                            pattern.sData[seqSlider] =  GetSite(otherDimension->lData[seqSlider])->sData[seqPos];
                        }
                    else
                        for (long seqSlider = 0; seqSlider < theMap.lLength; seqSlider ++) {
                            pattern.sData[seqSlider] =  GetSite(seqSlider)->sData[seqPos];
                        }

                    matches.Clear();
                    pattern.RegExpMatch (regex, matches);
                    if (matches.lLength) {
                        target << specCount;
                    }
                }

                if (eligibleSeqs != additionalFilter) {
                    DeleteObject (eligibleSeqs);
                }
            } else {
                bool         *eligibleMarks = new bool[lLength];
                checkPointer (eligibleMarks);

                for (long fillerID = 0; fillerID < lLength; fillerID++) {
                    eligibleMarks [fillerID] = false;
                }

                if (additionalFilter)
                    for (long siteIndex = 0; siteIndex < additionalFilter->lLength; siteIndex ++) {
                        eligibleMarks[theMap.lData[additionalFilter->lData[siteIndex]]] = true;
                    }
                else
                    for (long siteIndex = 0; siteIndex < lLength; siteIndex ++) {
                        eligibleMarks[siteIndex] = true;
                    }

                _SimpleList matches;
                _String     *tempString = nil;
                if (otherDimension) {
                    tempString = new _String (otherDimension->lLength,false);
                }

                for (long siteCounter = 0; siteCounter < lLength; siteCounter ++)
                    if (eligibleMarks[siteCounter]) {
                        matches.Clear();
                        if (otherDimension) {
                            _Site * aSite = ((_Site**)lData)[siteCounter];
                            for (long tc = 0; tc < otherDimension->lLength; tc++) {
                                tempString->sData[tc] = aSite->sData[otherDimension->lData[tc]];
                            }
                            tempString->RegExpMatch (regex, matches);
                        } else {
                            ((_Site**)lData)[siteCounter]->RegExpMatch (regex, matches);
                        }
                        if (matches.lLength == 0) {
                            eligibleMarks[siteCounter] = false;
                        }
                    }

                DeleteObject (tempString);

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
            FlushRegExp (regex);
        } else {
            input.KillSpaces (input);
            // now process the string
            long count = 0,anchor,k;

            _SimpleList numbers,
                        links;

            numbers.RequestSpace (1024);
            links.RequestSpace (1024);

            // first check if it is has a comb filter
            if ((input.sData[0]=='<')&&(input.sData[input.sLength-1]=='>')) {
                for (count=1; count<input.sLength-1; count++) {
                    if (input.sData[count]!='0') {
                        numbers<<count-1;
                    }
                }
                if (numbers.lLength) {
                    k = input.sLength-2; // step size
                    anchor = 0;
                    if (totalLength == -1) {
                        totalLength = theMap.lLength;
                    }
                    while (anchor<totalLength-k) {
                        for (count = 0; count< numbers.lLength; count++) {
                            target<<anchor+numbers.lData[count];
                        }
                        anchor+=k;
                    }
                    if ( (k=totalLength-1-anchor) ) {
                        for (count = 0; count< numbers.lLength; count++) {
                            if (numbers.lData[count]>k) {
                                break;
                            }
                            target<<anchor+numbers.lData[count];
                        }
                    }
                    return;

                }
            }

            while (count<input.sLength) {
                anchor = count;
                for (; (count<input.sLength)&&(isdigit(input[count])); count++) ;
                long    aNumber = (input.Cut (anchor,count-1)).toNum();
                if (aNumber < 0) {
                    _String warnMsg ("A negative number was found in partition specification: ");
                    ReportWarning (warnMsg & input.Cut (0,anchor-1) & '?' & input.Cut (anchor,-1));
                    target.Clear();
                    return;
                }
                numbers<< aNumber;

                if ((input[count]=='<')||(input[count]=='>')) {
                    _String warnMsg ("A comb partition cannot be combined with other types. The entire partition is reset to first..last");
                    ReportWarning (warnMsg & input.Cut (0,anchor-1) & '?' & input.Cut (anchor,-1));
                    target.Clear();
                    return;
                }
                if (input[count]=='&') {
                    links << numbers.lLength;
                }
                if ((input[count]==',')||(count==input.sLength)) { // wrap it up dude
                    if (numbers.lLength==1) {
                        target<<numbers(0);
                    } else {
                        if (links.lLength==0) {
                            if (numbers[0]>numbers[1]) { // backward order
                                for (k = numbers[0]; k>=numbers[1]; k--) {
                                    target<<k;
                                }
                            } else {
                                for (k = numbers[0]; k<=numbers[1]; k++) {
                                    target<<k;
                                }
                            }
                        } else {
                            // linked locations
                            if (links.lLength!=(numbers.lLength-2)/2) {
                                _String errMsg ("A part of the partition specification has not been understood and has been skipped.");
                                ReportWarning (errMsg);
                                target.Clear();
                                return;
                            } else {
                                _SimpleList signs;
                                signs<<(numbers(0)<numbers(1)?1:-1);
                                for (k = 0; k<links.lLength; k++) {
                                    signs<<(numbers(links(k))<numbers(links(k)+1)?1:-1);
                                }

                                long l,m;
                                for (k=numbers(0), l=0 ; signs(0)*k<=signs(0)*numbers(1); k+=signs(0), l++) {
                                    target<<numbers(0)+l*signs(0);
                                    for (m=0; m<links.lLength; m++) {
                                        target<<numbers(links(m))+l*signs(m+1);
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
//_______________________________________________________________________

void     _DataSetFilter::SetMap  (_String&s)
{
    theNodeMap.Clear();
    if (s.Length()) {
        long f,g=0;
        //_String sc(",");
        f = s.Find(',');
        while (f!=-1) {
            theNodeMap<<(long)s.Cut(g,f-1).toNum();
            g = f+1;
            f = s.Find (',',f+1,-1);
        }
        theNodeMap<<(long)s.Cut(g,-1).toNum();
    }
}

//_________________________________________________________

void    _DataSetFilter::FindAllSitesLikeThisOne (long index, _SimpleList& receptacle)
{
    long   oindex = theOriginalOrder.Find(index),m;

    if (oindex<0) {
        return;
    }

    if (theData->NoOfSpecies()==theNodeMap.lLength) {
        long *matchMap = new long[unitLength];

        checkPointer (matchMap);
        for (m=0; m<unitLength; m++)
            //matchMap[m] = theData->theMap.lData[theOriginalOrder.lData[oindex+1]];
        {
            matchMap[m] = theData->theMap.lData[theOriginalOrder.lData[oindex+m]];
        }


        for (long k=0; k<theOriginalOrder.lLength; k+=unitLength) {
            for (m=0; m<unitLength; m++) {
                if (theData->theMap.lData[theOriginalOrder.lData[k+m]]!=matchMap[m]) {
                    break;
                }
            }
            if (m==unitLength)
                for (m=0; m<unitLength; m++) {
                    receptacle<<theOriginalOrder.lData[k+m];
                }
        }

        delete matchMap;
    } else {
        char ** matchMap = (char**)MemAllocate (sizeof (char*) * unitLength);
        checkPointer (matchMap);

        for (m=0; m<unitLength; m++) {
            matchMap[m] = ((_Site*)(((BaseRef*)theData->lData)[theData->theMap.lData[oindex+m]]))->sData;
        }
        for (long k=0; k<theOriginalOrder.lLength; k+=unitLength) {
            for (m=0; m<unitLength; m++) {
                char* checkStr = ((_Site*)(((BaseRef*)theData->lData)[theData->theMap.lData[k+m]]))->sData;
                long t;
                for (t = 0; t<theNodeMap.lLength; t++) {
                    if (checkStr[t]!=matchMap[m][t]) {
                        break;
                    }
                }
                if (t<theNodeMap.lLength) {
                    break;
                }
            }
            if (m==unitLength)
                for (m=0; m<unitLength; m++) {
                    receptacle<<theOriginalOrder.lData[k+m];
                }
        }
        delete matchMap;
    }
}


//_______________________________________________________________________

_String&     _DataSetFilter::operator () (unsigned long site, unsigned long pos)
{
    if (!accessCache || accessCache->sLength != unitLength) {
        if (accessCache) {
            DeleteObject (accessCache);
        }
        checkPointer(accessCache = new _String ((unsigned long)unitLength, false));
    }

    long vIndex = theNodeMap.lData[pos];
    if (unitLength==1) {
        accessCache->sData[0]=(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site]]])->sData[vIndex];
    } else {
        site*=unitLength;
        for (int k = 0; k<unitLength; k++) {
            accessCache->sData[k] = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site++]]])->sData[vIndex];
        }
    }
    return *accessCache;
}

//_______________________________________________________________________

void     _DataSetFilter::RetrieveState (unsigned long site, unsigned long pos, _String& reply, bool map)
{
    long vIndex = theNodeMap.lData[pos];
    if (map) {
        if (unitLength==1) {
            reply.sData[0]=(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[duplicateMap.lData[site]]]])->sData[vIndex];
        } else {
            site = unitLength*duplicateMap.lData[site];
            for (int k = 0; k<unitLength; k++) {
                reply.sData[k] = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site++]]])->sData[vIndex];
            }
        }
    } else {
        if (unitLength==1) {
            reply.sData[0]=(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site]]])->sData[vIndex];
        } else
            for (int k = 0; k<unitLength; k++) {
                reply.sData[k] = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site++]]])->sData[vIndex];
            }
    }
}

//_______________________________________________________________________

void _DataSetFilter::GrabSite (unsigned long site, unsigned long pos, _String& storage)
{

    long vIndex = theNodeMap.lData[pos];
    if (unitLength==1) {
        storage.sData[0]=(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site]]])->sData[vIndex];
    } else {
        site*=unitLength;
        for (int k = 0; k<unitLength; k++) {
            storage.sData[k] = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site++]]])->sData[vIndex];
        }
    }
}

//_______________________________________________________________________

void _DataSetFilter::GrabSite (unsigned long site, unsigned long pos, char * s)
{
    long vIndex = theNodeMap.lData[pos];
    if (unitLength==1) {
        s[0]=(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site]]])->sData[vIndex];
    } else {
        site*=unitLength;
        for (int k = 0; k<unitLength; k++) {
            s[k] = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site++]]])->sData[vIndex];
        }
    }
}

//_______________________________________________________________________

_SimpleList* _DataSetFilter::CountAndResolve (long pattern, _Parameter * storage, bool randomly)
// last cell in the list contains the count of distinct characters in the column
{
    _SimpleList* resList = new _SimpleList (theNodeMap.lLength+1,0,0),
    counts (dimension,0,0);

    checkPointer (resList);

    _List        ambStates;
    _String      aState  (unitLength, false);

    _Parameter*  freqStorage = storage;

    if (!freqStorage) {
        freqStorage = new _Parameter [undimension];
    }

    long    normalizingSum = 0,
            charCount      = 0;

    for (long k=0; k<theNodeMap.lLength; k++) {
        GrabSite (pattern, k, aState);
        long      characterRes = Translate2Frequencies (aState, freqStorage, true);
        if (characterRes>=0) {
            resList->lData[k] = characterRes;

            if (characterRes >= dimension) {
                WarnError (_String("Internal error in _DataSetFilter::CountAndResolve\n"));
            }

            if ((counts.lData[characterRes]++) == 0) {
                normalizingSum ++;
            }

            charCount ++;
        } else {
            _SimpleList * possibleResolutions = new _SimpleList;
            if (!possibleResolutions) {
                checkPointer (possibleResolutions);
            }

            (*possibleResolutions) << k;

            for (long m=0; m<dimension; m++)
                if (freqStorage[m]>0.) {
                    (*possibleResolutions) << m;
                }

            ambStates.AppendNewInstance (possibleResolutions);
        }
    }

    if (normalizingSum > 0) {
        if (ambStates.lLength) {
            _SimpleList  ambResolutions (dimension,0,0);
            for (long t=0; t<ambStates.lLength; t++) {
                _SimpleList * stateResolutions = (_SimpleList*)ambStates(t);

                if (!randomly) {
                    long          totalSum = 0,
                                  idx = 0;

                    for (long l=1; l<stateResolutions->lLength; l++) {
                        long tmp = counts.lData[stateResolutions->lData[l]];
                        if (tmp>totalSum) {
                            idx = l;
                            totalSum = tmp;
                        }
                    }
                    if (idx > 0)
                        // if no resolutions, resolve randomly
                    {
                        idx = stateResolutions->lData[idx];
                        resList->lData[stateResolutions->lData[0]] = idx;
                        ambResolutions.lData [idx] ++;
                        continue;
                    }

                }

                long          totalSum = 0;
                for (long l=1; l<stateResolutions->lLength; l++) {
                    totalSum += counts.lData[stateResolutions->lData[l]];
                }

                if (totalSum > 0) {
                    long          randomN = genrand_real2() * totalSum - counts.lData[stateResolutions->lData[1]],
                                  ind = 1;

                    while (randomN > 0) {
                        randomN -= counts.lData[stateResolutions->lData[++ind]];
                    }

                    totalSum = stateResolutions->lData[ind];
                } else {
                    long          randomN = genrand_real2() * charCount - counts.lData[0],
                                  ind = 0;

                    while (randomN > 0) {
                        randomN -= counts.lData[++ind];
                    }
                }
                resList->lData[stateResolutions->lData[0]] = totalSum;
                ambResolutions.lData [totalSum] ++;
            }

            for (long l=0; l<dimension; l++)
                if (ambResolutions.lData[l] && !counts.lData[l]) {
                    normalizingSum ++;
                }
        }
    }

    resList->lData[theNodeMap.lLength] = normalizingSum;

    if (freqStorage != storage) {
        delete freqStorage;
    }

    return       resList;
}

//_______________________________________________________________________

_Matrix* _DataSetFilter::PairwiseCompare (_SimpleList* s1, _SimpleList *s2, _List* labels)
// s1 and s2 are the lists produced by CountAndResolve
// if labels is not nil, then it will receive row and column labels in the contigency table
// the result matrix has rows labeled by states in s1, and columns - by states in s2
{
    long    * sort1 = new long[dimension],
    * sort2 = new long[dimension],
    c = s2->lData[s2->lLength-1];

    _Matrix * res   = new _Matrix (s1->lData[s1->lLength-1],c,false,true);

    if (sort1 && sort2 && res) {
        for (long k = 0; k<dimension; k++) {
            sort1[k] = -1;
            sort2[k] = -1;
        }

        long idx1 = 0,
             idx2 = 0;

        _SimpleList  *lbl1 = nil,
                      *lbl2 = nil;

        if (labels) {
            lbl1 = new _SimpleList;
            lbl2 = new _SimpleList;

            checkPointer (lbl1);
            checkPointer (lbl2);

            (*labels) << lbl1;
            (*labels) << lbl2;

            DeleteObject (lbl1);
            DeleteObject (lbl2);
        }

        for (long k2 = 0; k2 < s1->lLength-1; k2++) {
            long c1 = s1->lData[k2],
                 c2 = s2->lData[k2];

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

            res->theData[c1*c+c2] += 1.;
        }

        delete sort1;
        delete sort2;
    } else {
        checkPointer (nil);
    }

    return res;
}

//_______________________________________________________________________

_List *  _DataSetFilter::ComputePatternToSiteMap (void)
{
    _List * result = new _List ();
    for (long k = 0; k < theFrequencies.lLength; k++) {
        result->AppendNewInstance (new _SimpleList);
    }
    for (long s = 0; s < duplicateMap.lLength; s++) {
        *((_SimpleList**)result->lData)[duplicateMap.lData[s]] << s;
    }
    return result;
}

//_______________________________________________________________________

char     _DataSetFilter::GetChar (unsigned long site, unsigned long pos)
{
    //long vIndex = theNodeMap.lLength?theNodeMap.lData[pos]:pos;
    return (*theData)(theMap.lData[site],theNodeMap.lData[pos],1);
}
//_______________________________________________________________________

bool     _DataSetFilter::CompareTwoSites (unsigned long site1, unsigned long site2, unsigned long pos1)
{
    pos1 = theNodeMap.lData[pos1];

    if (unitLength == 3) { // codon
        site1*=3;
        site2*=3;
        if (
            ((((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1]]])->sData[pos1]==
             ( ((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site2]]])->sData[pos1])
            &&((((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1+1]]])->sData[pos1]==
               (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site2+1]]])->sData[pos1])
            &&((((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1+2]]])->sData[pos1]==
               (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site2+2]]])->sData[pos1])) {
            return true;
        }
    } else {
        site1*=unitLength;
        site2*=unitLength;
        long k;

        /*if ((((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1]]])->sLength<=pos1)
        {
            printf ("(%d)%s\n(%d)%s\n",site1,(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1]]])->sData,
                    site2,(((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site2]]])->sData);
            FlagError ("Internal DataSetFilter bug\n");
        }*/

        for (k = 0; k<unitLength; k++) {
            if ((((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site1]]])->sData[pos1]!=
                    (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[site2]]])->sData[pos1]) {
                break;
            }
            site1++;
            site2++;
        }
        if (k==unitLength) {
            return true;
        }
    }
    return false;
}

//_______________________________________________________________________

bool     _DataSetFilterNumeric::CompareTwoSites (unsigned long, unsigned long, unsigned long)
{
    return false;
}


//_______________________________________________________________________

bool     _DataSetFilter::CompareTwoSitesChar (unsigned long site1, unsigned long site2, unsigned long pos1)
{

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
    return ((*theData)(theMap.lData[site1],pos1, 1)==(*theData)(theMap.lData[site2],pos1, 1));
}
//_______________________________________________________________________
long    _DataSetFilter::SiteFrequency (unsigned long site)
{
    return theFrequencies.lData[site];
}

//_______________________________________________________________________
bool    _DataSetFilter::HasDeletions (unsigned long site, _AVLList* storage)
{
    long        loopDim  = GetDimension();
    _Parameter* store    = new _Parameter [loopDim];

    if (!store) {
        warnError( -108);
    }

    long j,
         upTo = theNodeMap.lLength?theNodeMap.lLength:theData->NoOfSpecies();

    bool outcome = false;

    for (unsigned int k = 0; k<upTo; k++) {
        Translate2Frequencies ((*this)(site,k), store, false);

        bool oneF = false,
             zeroF = false;

        for (j=0; j<loopDim; j++) {
            if (store[j]==0.0) {
                zeroF = true;
            } else if (store[j]==1.0) {
                oneF = true;
            }
        }
        if (!(oneF&&zeroF)) {
            if (storage) {
                outcome = true;
                storage->Insert ((BaseRef)theNodeMap.lData[k]);
            } else {
                delete store;
                return true;
            }
        }
    }

    delete store;
    return outcome;
}


//_______________________________________________________________________
bool    _DataSetFilter::IsConstant (unsigned long site,bool relaxedDeletions)
{
    _Parameter* store = nil, *store2 = nil;

    store = new _Parameter [GetDimension()];
    store2 = new _Parameter [GetDimension()];

    if (!(store&&store2)) {
        warnError(-108);
    }
    long j,
         upTo = theNodeMap.lLength?theNodeMap.lLength:theData->NoOfSpecies(),
         loopDim = GetDimension();

    Translate2Frequencies ((*this)(site,0), store, false);

    if (relaxedDeletions) {
        for (unsigned int k = 1; k<upTo; k++) {
            Translate2Frequencies ((*this)(site,k), store2, false);
            for (j=0; j<loopDim; j++) {
                if (store2[j]==0.0) {
                    store[j]=0.0;
                }
            }
        }
        for (j=0; j<loopDim; j++)
            if (store[j]!=0.0) {
                delete store;
                delete store2;
                return true;
            }
        if (j==loopDim) {
            delete store;
            delete store2;
            return false;
        }
    } else {
        for (unsigned int k = 1; k<upTo; k++) {
            Translate2Frequencies ((*this)(site,k), store2, false);
            for (j=0; j<loopDim; j++)
                if (store[j]!=store2[j]) {
                    delete store;
                    delete store2;
                    return false;
                }
        }
    }

    return true;
}

//_______________________________________________________________________

_Matrix*        _DataSetFilter::GetFilterCharacters (bool flip)
{
    long        unitLength = GetUnitLength (),
                seqLength  = flip?theFrequencies.lLength:(GetFullLengthSpecies () / unitLength),
                f          = NumberSpecies();

    _List       result;

    _String      aState ((long)GetUnitLength(),false);

    if (flip) {
        for (long k=0; k<seqLength; k++) {
            _String *aSite = new _String (128L,true);
            for (long k2=0; k2<f; k2++) {
                RetrieveState(k,k2,aState,false);
                (*aSite) << aState;
            }
            aSite->Finalize();
            result      << aSite;
            DeleteObject (aSite);
        }
    } else
        for (long k=0; k<f; k++) {
            _String *fstr = GetSequenceCharacters(k);
            result      << fstr;
            DeleteObject (fstr);
        }

    return new _Matrix (result);
}

//_______________________________________________________________________

_String*        _DataSetFilter::GetSequenceCharacters (long seqID)
{
    long            unitSizeL   = GetUnitLength();
    _String * aSequence = new _String (theOriginalOrder.lLength,true);

    if (seqID >= 0 && seqID < theNodeMap.lLength) {
        _String      aState (unitSizeL,false);
        long        upTo = theOriginalOrder.lLength/unitSizeL;
        for (long k2=0; k2<upTo; k2++) {
            RetrieveState(k2,seqID,aState);
            (*aSequence) << aState;
        }
    }
    aSequence->Finalize();
    return aSequence;
}

//_______________________________________________________________________
long    _DataSetFilter::HasExclusions (unsigned long site, _SimpleList* theExc, _Parameter*store )
{
    long   filterDim = GetDimension(false);

    if (theNodeMap.lLength)
        for (unsigned long k = 0; k<theNodeMap.lLength; k++) {
            Translate2Frequencies   ((*this)(site,k), store, false);

            long                    j                       = 0,
                                    s                      = 0;

            for (j=0; j<filterDim; j++)
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
//_______________________________________________________________________
void    _DataSetFilter::Freeze (long site)
{
    for (int k = 0; k<unitLength; k++) {
        _Site* tC = (_Site*)((*(_List*)theData)(theData->theMap(this->theMap(site*unitLength+k))));
        tC->SetRefNo(-1);
        tC->PrepareToUse();
    }
}

//_______________________________________________________________________
void    _DataSetFilter::UnFreeze (long site)
{
    for (int k = 0; k<unitLength; k++) {
        _Site* tC = (_Site*)((*(_List*)theData)(theData->theMap(this->theMap(site*unitLength+k))));
        tC->SetRefNo(0);
//      tC->Archive();
    }
}

//_________________________________________________________
_Matrix* _DataSetFilter::ComputePairwiseDifferences (long i, long j, char amb)
{
    if (unitLength > 3) {
        WarnError ("ComputePairwiseDifferences is not implemented for data filters with unit size > 3");
        return    new _Matrix (1,1,false,true);
    }

    long    mxDim       = GetDimension (true);

    _Matrix*res         = new _Matrix  (mxDim,mxDim,false,true);

    _Parameter
    *sm1   = new _Parameter[mxDim],
    *sm2   = new _Parameter[mxDim];


    checkPointer (res);
    checkPointer (sm1);
    checkPointer (sm2);

    _String      state1 (unitLength,false),
                 state2 (unitLength,false);


    if (!conversionCache.lLength) {
        SetupConversion();
    }

    long        *tcodes  = conversionCache.lData+89,
                 *ccodes  = conversionCache.lData+1,
                  ccount   = conversionCache.lData[0];

    for (long k=0; k<theFrequencies.lLength; k++) {
        long s1 = -1,
             s2 = -1;

        int c1, c2;

        c1 = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[unitLength*k]]])->sData[theNodeMap.lData[i]],
        c2 = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[unitLength*k]]])->sData[theNodeMap.lData[j]];

        if (unitLength == 1) {
            s1 = conversionCache.lData[(c1-40)*(undimension+1)+undimension],
            s2 = conversionCache.lData[(c2-40)*(undimension+1)+undimension];
        } else {
            int         c12 = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[unitLength*k+1]]])->sData[theNodeMap.lData[i]],
                        c22 = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[unitLength*k+1]]])->sData[theNodeMap.lData[j]];


            state1.sData[0] = c1;
            state1.sData[1] = c12;

            state2.sData[0] = c2;
            state2.sData[1] = c22;

            c1  = ccodes[c1-40];
            c12 = ccodes[c12-40];

            c2  = ccodes[c2-40];
            c22 = ccodes[c22-40];

            if (unitLength == 2) {
                if ((c1>=0)&&(c12>=0)) {
                    s1 = tcodes[c1*ccount+c12];
                }

                if ((c2>=0)&&(c22>=0)) {
                    s2 = tcodes[c2*ccount+c22];
                }
            } else {
                int         c13 = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[unitLength*k+2]]])->sData[theNodeMap.lData[i]],
                            c23 = (((_String**)theData->lData)[theData->theMap.lData[theMap.lData[unitLength*k+2]]])->sData[theNodeMap.lData[j]];

                //printf ("\n%c %c", c13, c23);

                state1.sData[2] = c13;
                state2.sData[2] = c23;

                c13 = ccodes[c13-40];
                c23 = ccodes[c23-40];

                //printf (" %d %d %s %s\n", c13, c23, state1.sData, state2.sData);

                if ((c1>=0)&&(c12>=0)&&(c13>=0)) {
                    s1 = tcodes[ccount*(c1*ccount+c12)+c13];
                }

                if ((c2>=0)&&(c22>=0)&&(c23>=0)) {
                    s2 = tcodes[ccount*(c2*ccount+c22)+c23];
                }
            }
        }

        if (s1>=0 && s2>=0)
            // one to one
        {
            res->theData[s1*mxDim+s2] += theFrequencies.lData[k];
        } else {
            if (amb<3) {
                _Matrix * freqsAtSite = nil;
                if (amb) {
                    _SimpleList   //seqList,
                    siteList;

                    //seqList  << theNodeMap[i];
                    //seqList  << theNodeMap[j];

                    for (long si = 0; si < unitLength; si++) {
                        siteList << theMap.lData[unitLength*k+si];
                    }

                    freqsAtSite     = theData->HarvestFrequencies (unitLength, unitLength, 0, theNodeMap, siteList);
                    if (theExclusions.lLength) {
                        long k = 0,
                             u = GetDimension (false);

                        for (long i = 0; i<u; i++) {
                            if (i==theExclusions.lData[k] && k<theExclusions.lLength) {
                                k++;
                                continue;
                            }
                            freqsAtSite->theData[i-k] = freqsAtSite->theData[i];
                        }
                    }
                    //XferwCorrection (freqsAtSite->theData, freqsAtSite->theData, mxDim);
                }

                if (s1>=0)
                    // one to many
                {
                    if (unitLength>1) {
                        Translate2Frequencies (state2,sm1,false);
                    } else {
                        Translate2Frequencies (c2,sm1,false);
                    }

                    if (freqsAtSite) {
                        if (amb == 1) {
                            _Parameter totalW = 0.0;

                            for  (long m=0; m<mxDim; m++)
                                if (sm1[m]>0.0) {
                                    totalW += freqsAtSite->theData[m];
                                }

                            if (totalW>0.0) {
                                s1 = s1*mxDim;

                                for  (long m=0; m<mxDim; m++,s1++)
                                    if (sm1[m]>0.0) {
                                        res->theData[s1] += theFrequencies.lData[k]*freqsAtSite->theData[m]/totalW;
                                    }
                            }

                        } else {
                            _Parameter maxW   = 0.0;
                            long       maxIdx = -1;

                            for  (long m=0; m<mxDim; m++) {
                                if (sm1[m]>0.0) {
                                    _Parameter myWeight = freqsAtSite->theData[m];
                                    if (myWeight > maxW) {
                                        maxW = myWeight;
                                        maxIdx = m;
                                    }
                                }
                            }

                            if (maxIdx>=0) {
                                res->theData[s1*mxDim+maxIdx] += theFrequencies.lData[k];
                            }
                        }
                    } else {
                        /* adopt the following convention here:
                            - if ambig resolves to one s1 - count as a match
                            - otherwise - count all contributions equally
                        */

                        if (sm1[s1] > 0.0) {
                            res->theData[s1*mxDim+s1] += theFrequencies.lData[k];
                        } else {
                            long ambCount = 0;
                            {
                                for  (long m=0; m<mxDim; m++)
                                    if (sm1[m]>0.0) {
                                        ambCount ++;
                                    }
                            }
                            s1 *= mxDim;

                            _Parameter addFac = theFrequencies.lData[k]/(_Parameter)ambCount;

                            for  (long m=0; m<mxDim; m++,s1++)
                                if (sm1[m]>0.0) {
                                    res->theData[s1] += addFac;
                                }
                        }
                    }
                } else {
                    if (s2>=0)
                        // many to one
                    {
                        if (unitLength>1) {
                            Translate2Frequencies (state1,sm1,false);
                        } else {
                            Translate2Frequencies (c1,sm1,false);
                        }

                        if (freqsAtSite) {
                            if (amb == 1) {
                                _Parameter totalW = 0.0;

                                for  (long m=0; m<mxDim; m++)
                                    if (sm1[m]>0.0) {
                                        totalW += freqsAtSite->theData[m];
                                    }

                                if (totalW>0.0) {
                                    for  (long m=0; m<mxDim; m++,s2+=mxDim)
                                        if (sm1[m]>0.0) {
                                            res->theData[s2] += theFrequencies.lData[k]*freqsAtSite->theData[m]/totalW;
                                        }
                                }

                            } else {
                                _Parameter maxW   = 0.0;
                                long       maxIdx = -1;

                                for  (long m=0; m<mxDim; m++) {
                                    if (sm1[m]>0.0) {
                                        _Parameter myWeight = freqsAtSite->theData[m];
                                        if (myWeight > maxW) {
                                            maxW = myWeight;
                                            maxIdx = m;
                                        }
                                    }
                                }

                                if (maxIdx>=0) {
                                    res->theData[maxIdx*mxDim+s2] += theFrequencies.lData[k];
                                }
                            }
                        } else {
                            if (sm1[s2] > 0.0) {
                                res->theData[s2*mxDim+s2] += theFrequencies.lData[k];
                            } else {
                                long ambCount = 0;
                                for  (long m=0; m<mxDim; m++)
                                    if (sm1[m]>0.0) {
                                        ambCount ++;
                                    }

                                _Parameter addFac = theFrequencies.lData[k]/(_Parameter)ambCount;
                                {
                                    for  (long m=0; m<mxDim; m++,s2+=mxDim)
                                        if (sm1[m]>0.0) {
                                            res->theData[s2] += addFac;
                                        }
                                }
                            }
                        }
                    } else
                        // many to many
                    {
                        if (unitLength>1) {
                            Translate2Frequencies (state1,sm1,false);
                            Translate2Frequencies (state2,sm2,false);
                        } else {
                            Translate2Frequencies (c1,sm1,false);
                            Translate2Frequencies (c2,sm2,false);
                        }

                        if (freqsAtSite) {
                            if (amb == 1) {
                                _Parameter totalW = 0.0;

                                for  (long m=0; m<mxDim; m++)
                                    if (sm1[m]>0)
                                        for  (long m2=0; m2<mxDim; m2++)
                                            if (sm2[m2]>0) {
                                                totalW += freqsAtSite->theData[m]*freqsAtSite->theData[m2];
                                            }

                                if (totalW>0.0) {
                                    for  (long m=0; m<mxDim; m++)
                                        if (sm1[m]>0)
                                            for  (long m2=0; m2<mxDim; m2++)
                                                if (sm2[m2]>0) {
                                                    res->theData[m*mxDim+m2] += theFrequencies.lData[k]*freqsAtSite->theData[m]*freqsAtSite->theData[m2]/totalW;
                                                }
                                }

                            } else {
                                _Parameter maxW   = 0.0;
                                long       maxIdx  = -1,
                                           maxIdx2 = -1;

                                for  (long m=0; m<mxDim; m++)
                                    if (sm1[m]>0)
                                        for  (long m2=0; m2<mxDim; m2++)
                                            if (sm2[m2]>0) {
                                                _Parameter myWeight = freqsAtSite->theData[m]*freqsAtSite->theData[m2];
                                                if (myWeight > maxW) {
                                                    maxW = myWeight;
                                                    maxIdx  = m;
                                                    maxIdx2 = m2;
                                                }
                                            }

                                if (maxIdx>=0) {
                                    res->theData[maxIdx*mxDim+maxIdx2] += theFrequencies.lData[k];
                                }
                            }
                        } else {
                            long ambCount  = 0,
                                 ambCount2 = 0,
                                 m         = 0;

                            for  (; m<mxDim; m++) {
                                if (sm1[m]>0.0) {
                                    if (sm2[m]>0.0) {
                                        break;
                                    } else {
                                        ambCount ++;
                                    }
                                } else if (sm2[m]>0.0) {
                                    ambCount2 ++;
                                }
                            }

                            if (m==mxDim) {
                                _Parameter addFac = theFrequencies.lData[k]/(_Parameter)(ambCount*ambCount2);

                                for  (long m=0; m<mxDim; m++)
                                    if (sm1[m]>0)
                                        for  (long m2=0; m2<mxDim; m2++)
                                            if (sm2[m2]>0) {
                                                res->theData[m*mxDim+m2] += addFac;
                                            }
                            }
                        }
                    }
                }
                DeleteObject (freqsAtSite);
            }
        }
    }

    delete[] sm1;
    delete[] sm2;

    return res;
}

//_________________________________________________________

void _DataSetFilter::ComputePairwiseDifferences (_Matrix& target, long i, long j)
// matrix of dimension nx4n containing pairwise distances as follows (n=number of species)
// first lower diag - count the same (AA,CC,GG,TT)
// first upper diag - count AC,CA
// 2nd   lower diag - count AG,GA
// 2nd   upper diag - count AT,TA
// 3rd   lower diag - count CG,GC
// 3rd   upper diag - count CT,TC
// 4th   lower diag - count GT,TG
{
    if ((target.GetHDim()!=1)||(target.GetVDim()!=7)) {
        CreateMatrix (&target,1,7,false,true,false);
    }

    if (!theData->theTT->IsStandardNucleotide()) {
        return;
    }
    long k,l,m;

    for (k=0; k<7; k++) {
        target[k] = 0;
    }
    k = theNodeMap.lData[i];
    l = theNodeMap.lData[j];
    if (l>k) {
        m=l;
        l=k;
        k=m;
    }

    for (m=theMap.lLength-1; m>-1; m--) {
        char * thisSite = GetColumn (m);
        char a = thisSite[k],
             b = thisSite[l],
             c;

        long fc = theFrequencies.lData[m/unitLength];

        if (a>b) {
            c=a;
            a=b;
            b=c;
        }
        if (a==b) {
            target[0]+=fc;
        } else {
            if (a=='A') {
                switch (b) {
                case 'C': {
                    target[1]+=fc;
                    break;
                }
                case 'G': {
                    target[2]+=fc;
                    break;
                }
                case 'T': {
                    target[3]+=fc;
                    break;
                }
                }
            } else if (a=='C') {
                switch (b) {
                case 'G': {
                    target[4]+=fc;
                    break;
                }
                case 'T': {
                    target[5]+=fc;
                    break;
                }
                }
            } else if (a=='G') {
                if (b=='T') {
                    target[6]+=fc;
                }
            }

        }
    }
}

//_________________________________________________________

_Matrix * _DataSetFilter::HarvestFrequencies (char unit, char atom, bool posSpec, bool countGaps)
{
    return theData->HarvestFrequencies (unit,atom, posSpec, theNodeMap, theOriginalOrder, countGaps);
}
//_______________________________________________________________________
void    _DataSetFilter::XferwCorrection (_Matrix& source, _Parameter* target, long _length)
{
    long k=0;
    _Parameter* mxdata = source.fastIndex();
    if (theExclusions.lLength==0) {
        for (long i = 0; i<_length; i++) {
            target[i] = (mxdata[i]!=0.0);
        }
    } else {
        for (long i = 0; i<_length; i++) {
            if (k<theExclusions.lLength && i==theExclusions.lData[k]) {
                k++;
                continue;
            }
            target[i-k] = (mxdata[i] != 0.);
        }
    }
}
//_______________________________________________________________________
void    _DataSetFilter::XferwCorrection (_Parameter* source, _Parameter* target, long _length)
{
    long k=0;
    if (theExclusions.lLength==0) {
        for (long i = 0; i<_length; i++) {
            target[i] = (source[i]!=0.0);
        }
    } else {
        for (long i = 0; i<_length; i++) {
            if (i==theExclusions.lData[k] && k<theExclusions.lLength) {
                k++;
                continue;
            }
            target[i-k] = (source[i]!=0);
        }
    }
}
//_______________________________________________________________________

void    _DataSetFilter::XferwCorrection (long* source, _Parameter* target, long _length)
{
    long k=0;
    if (theExclusions.lLength==0) {
        for (long i = 0; i<_length; i++) {
            target[i] = source[i];
        }
    } else {
        for (long i = 0; i<_length; i++) {
            if (i==theExclusions[k]) {
                k++;
                continue;
            }
            target[i-k] = source[i];
        }
    }
}

/*
//_______________________________________________________________________
long    _DataSetFilter::GetVectorCode(long site,long seq)
{
    if (!symbolVector) return -1;
    long* fi = symbolVector->quickArrayAccess();
    return fi[*fi*site+seq+1];
}
//_______________________________________________________________________

void    _DataSetFilter::ProduceSymbolVector(bool smear)
{
    // compute the size of the vector cells
    _Parameter cellSize=log((_Parameter)theData->theTT->LengthOfAlphabet())*_Parameter(unitLength)/log(128.0);
    if (cellSize>2.0)
    {
        _String errMsg ("DataSetFilter has more than 32767 states, which is currently unsupported");
        FlagError(errMsg);
    }
    long intCellSize = cellSize>1.0?2:1;
    // now produce the conversion vector
    long sites = theMap.lLength, species= theNodeMap.lLength?theNodeMap.lLength:theData->NoOfSpecies();
    symbolVector = new _SimpleList ();
    checkPointer(symbolVector);
//  (*symbolVector)<<intCellSize;
    (*symbolVector)<<species;
    // we will now speciate into byte and word size cases
    // the data will be stored column by column
    // if there is a unique code translation, we then store that code in the symbol vector for faster
    // processing during tree pruning business.
    // use a standard convert to frequencies function to check whether a character has a unique convertion
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
//              converterb.bytes[byteposition]=(char)Translate2Frequencies((*this)(i,j),nil,smear,false);
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

//_______________________________________________________________________

long    _DataSetFilter::CorrectCode (long code)
{
    if (theExclusions.lLength!=0) {
        for (long k=0; k<theExclusions.lLength; k++)
            if (code>=theExclusions.lData[k]) {
                code++;
            }
    }
    return code;
}

//_______________________________________________________________________
long    _DataSetFilter::Translate2Frequencies (_String& str, _Parameter* parvect, bool smear)
{
    long count      = 0,
         nonzeropos = 0,
         store      [HYPHY_SITE_DEFAULT_BUFFER_SIZE];


    if (unitLength == 1) {
        theData->theTT->TokenCode (str.sData[0],store,smear);
        if (theExclusions.lLength==0) {
            for (long i = 0; i<undimension; i++)
                if ( (parvect[i]=store [i]) ) {
                    nonzeropos = i;
                    count++;
                }
        } else {
            long k=0;
            for (long i = 0; i<undimension; i++) {
                if (i==theExclusions.lData[k] && k<theExclusions.lLength) {
                    k++;
                } else if ( store [i] ) {
                    nonzeropos = i;
                    count++;
                }
                parvect[i-k]=store[i];
            }
        }
        if (count == 0) {
            if (smear)
                for (long i = 0; i<undimension; i++) {
                    parvect[i] = 1.;
                }

            return -1;
        }

        if (count>1) {
            return -1;
        }

        return nonzeropos;
    } else {
        //pull the frequencies out of the Translation table
        _Matrix     out (undimension,1,false,true);

        long        m,
                    n,
                    index = 0,
                    shifter = 1,
                    *lp,
                    *storeP;

        _Parameter* fl;

        if (theData->theTT->baseLength * unitLength >= HYPHY_SITE_DEFAULT_BUFFER_SIZE) {
            storeP = new long [theData->theTT->baseLength * unitLength];
        } else {
            storeP = store;
        }

        count = 1;
        for (m = 0; m<unitLength; m++ ) {
            theData->theTT->TokenCode (str.sData[m], storeP+theData->theTT->baseLength*m);
        }

        for (m = unitLength-1; m>=0; m--) {
            int smcount = 0;
            lp = storeP+theData->theTT->baseLength*m;
            for (n = 0; n<theData->theTT->baseLength; n++,lp++) {
                if (*lp) {
                    index += shifter*n;
                    smcount++;
                }
            }
            if ((smcount==0)&&smear) { // deletion -- replace with 1's
                lp = storeP+theData->theTT->baseLength*m;
                for (n = 0; n<theData->theTT->baseLength; n++,lp++) {
                    *lp=1;
                }
                smcount = theData->theTT->baseLength;
            }

            shifter*=theData->theTT->baseLength;
            count *=smcount;
        }


        if (count>1) {
            theData->constructFreq (storeP, &out, 1, 0, count, unitLength-1 , 1, 0);
        } else if (count == 1) {
            out.theData[index] = count;
        }

        if (storeP != store) {
            delete [] storeP;
        }

        if (count==1) {
            fl = out.theData;
            m = 0;
            if (theExclusions.lLength) {
                for (n = 0; n<undimension; n++, fl++) {
                    if (m <theExclusions.lLength && n==theExclusions.lData[m]) {
                        m++;
                        continue;
                    }


                    if (*fl>0.) {
                        parvect[n-m] = 1.0;
                        break;
                    } else {
                        parvect[n-m] = 0.0;
                    }
                }
                if (n<undimension) {
                    n-=m-1;
                    parvect+=n;
                    for (; n<dimension; n++,parvect++) {
                        *parvect = 0.0;
                    }
                    return index-m;
                } else {
                    if (smear)
                        for (n = 0; n<dimension; n++, parvect++) {
                            *parvect = 1.0;
                        }
                    return -1;
                }

            } else {
                if (count) {
                    for (n = 0; n<undimension; n++,fl++) {
                        parvect[n] = *fl>0.0?1.0:0.0;
                    }
                }
                return index;
            }
        } else {
            XferwCorrection(out, parvect,undimension);

            if (smear) {
                for (count = 0; count<dimension; count++)
                    if (parvect[count]>0.0) {
                        break;
                    }
                if (count==dimension) {
                    for (count = 0; count<dimension; count++) {
                        parvect[count]=1.0;
                    }
                }
            }
            return -1;
        }
    }
    return 0;
}


//_______________________________________________________________________

long    _DataSetFilter::MapStringToCharIndex (_String& str)
{
    long    count       = 0,
            nonzeropos  = 0,
            store       [HYPHY_SITE_DEFAULT_BUFFER_SIZE];

    if (unitLength == 1) {
        theData->theTT->TokenCode (str.sData[0],store);

        if (theExclusions.lLength==0) {
            for (long i = 0; i<undimension; i++)
                if (store [i]) {
                    nonzeropos = i;
                    count++;
                }
        } else {
            long k=0;
            for (long i = 0; i<undimension; i++) {
                if (i==theExclusions.lData[k] && k<theExclusions.lLength) {
                    k++;
                } else if ( store [i] ) {
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
        long    m,
                n,
                index = 0,
                shifter = 1,
                *lp,
                *storeP;

        count = 1;

        if (theData->theTT->baseLength * unitLength >= HYPHY_SITE_DEFAULT_BUFFER_SIZE) {
            storeP = new long [theData->theTT->baseLength * unitLength];
        } else {
            storeP = store;
        }

        for (m = 0; m<unitLength; m++ ) {
            theData->theTT->TokenCode (str.sData[m], storeP+theData->theTT->baseLength*m);
        }

        for (m = unitLength-1; m>=0; m--) {
            int smcount = 0;

            lp = storeP+theData->theTT->baseLength*m;

            for (n = 0; n<theData->theTT->baseLength; n++,lp++)
                if (*lp) {
                    index += shifter*n;
                    smcount++;
                }

            if (smcount==0) {
                lp = storeP+theData->theTT->baseLength*m;

                for (n = 0; n<theData->theTT->baseLength; n++,lp++) {
                    *lp=1;
                }

                smcount = theData->theTT->baseLength;
            }

            shifter *= theData->theTT->baseLength;

            count *=smcount;
        }


        if (storeP != store) {
            delete [] storeP;
        }

        if (count==1) {
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


//_______________________________________________________________________
long    _DataSetFilter::Translate2Frequencies (char s, _Parameter* parvect, bool smear)
{
    long     count =0,
             store [HYPHY_SITE_DEFAULT_BUFFER_SIZE];

    theData->theTT->TokenCode (s,store, smear);

    if (theExclusions.lLength) {
        long k = 0;
        for (long i = 0; i<undimension; i++) {
            if (i==theExclusions[k]) {
                k++;
            } else if ( store [i] ) {
                count++;
            }
        }
        if (count) {
            XferwCorrection (store, parvect,undimension);
        }
    } else {
        for (long i = 0; i<undimension; i++)
            if ( (parvect[i]=store [i]) ) {
                count++;
            }
    }

    if (count == 0 && smear)
        for (long i = 0; i<undimension; i++) {
            parvect[i] = 1.0;
        }

    return count==1?1:-1;
}

//_______________________________________________________________________
long    _DataSetFilter::LookupConversion (char s, _Parameter* parvect)
{
    if (undimension==4) {
        //int idx = (s-40)*5;
        long* cCache = conversionCache.lData+(s-40)*5;
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
        int idx = (s-40)*(undimension+1);
        for (long i=0; i<undimension; parvect[i++] = conversionCache.lData[idx++]) ;
        return conversionCache.lData[idx];
    }
}

//_______________________________________________________________________
void    _DataSetFilter::SetupConversion (void)
{
    if (conversionCache.lLength) {
        return;
    }

    if ( unitLength==1 ) { // do stuff
        char c = 40;

        long i,
             charCount,
             theCode;

        _Parameter *temp    = new _Parameter [undimension+1];

        while(c<127) {
            for (i=0; i<undimension; i++) {
                temp[i]=0;
            }

            Translate2Frequencies(c,temp,true);

            charCount = -1;
            for (i=0; i<undimension; i++) {
                theCode =  (long)temp[i];
                conversionCache <<theCode;
                if (theCode) {
                    if (charCount==-1) {
                        charCount = i;
                    } else {
                        charCount = -2;
                    }
                }
            }
            conversionCache<<charCount;
            c++;
        }
        delete[] temp;
    } else {
        if (unitLength==2 || unitLength==3) {
            _String alphabet (16,true);
            if (theData->theTT->baseSet.sLength == 0) {
                if (theData->theTT->baseLength == 4) {
                    alphabet << 'A';
                    alphabet << 'C';
                    alphabet << 'G';
                    alphabet << 'T';
                } else {
                    if (theData->theTT->baseLength == 20) {
                        alphabet << 'A';
                        alphabet << 'C';
                        alphabet << 'D';
                        alphabet << 'E';
                        alphabet << 'F';
                        alphabet << 'G';
                        alphabet << 'H';
                        alphabet << 'I';
                        alphabet << 'K';
                        alphabet << 'L';
                        alphabet << 'M';
                        alphabet << 'N';
                        alphabet << 'P';
                        alphabet << 'Q';
                        alphabet << 'R';
                        alphabet << 'S';
                        alphabet << 'T';
                        alphabet << 'V';
                        alphabet << 'W';
                        alphabet << 'Y';
                    } else {
                        alphabet << '0';
                        alphabet << '1';
                    }
                }
            } else {
                alphabet << &theData->theTT->baseSet;
            }


            alphabet.Finalize();


            long  ccache [88],
                  i,
                  k = GetDimension(false) ;

            conversionCache.RequestSpace (89+k);
            conversionCache << alphabet.sLength;

            for (i=0; i<88; i++) {
                ccache[i] = -1;
            }
            for (i=0; i<alphabet.sLength; i++) {
                ccache [alphabet.sData[i]-40] = i;
            }
            for (i=0; i<88; i++) {
                conversionCache << ccache[i];
            }

            long       *tcache = new long [k];
            checkPointer (tcache);
            //_Parameter *vcache = new _Parameter [k];
            //checkPointer (vcache);

            if (unitLength == 3) {
                _String s (3,false);
                i = 0;
                for (long a = 0; a<alphabet.sLength; a++,i+=alphabet.sLength*alphabet.sLength) {
                    s.sData[0] = alphabet.sData[a];
                    for (long b = 0; b<alphabet.sLength; b++) {
                        s.sData[1] = alphabet.sData[b];
                        for (long c = 0; c<alphabet.sLength; c++) {
                            s.sData[2] = alphabet.sData[c];
                            //tcache [i+b*alphabet.sLength+c] = Translate2Frequencies (s,vcache,true);
                            tcache [i+b*alphabet.sLength+c] = MapStringToCharIndex (s);
                        }
                    }
                }
            } else {
                _String s (2,false);
                i = 0;
                for (long a = 0; a<alphabet.sLength; a++,i+=alphabet.sLength) {
                    s.sData[0] = alphabet.sData[a];
                    for (long b = 0; b<alphabet.sLength; b++) {
                        s.sData[1] = alphabet.sData[b];
                        //tcache [i+b] = Translate2Frequencies (s,vcache,true);
                        tcache [i+b] = MapStringToCharIndex (s);
                    }
                }
            }
            for (i=0; i<k; i++) {
                conversionCache << tcache[i];
            }

            //delete vcache;
            delete [] tcache;
        }
    }
}



//_________________________________________________________
//_________________________________________________________
// reading the data set file in here



//_________________________________________________________
void    checkTTStatus (FileState* fs) // check whether the translation table needs to be refreshed
{
    if (fs->translationTable == &defaultTranslationTable) {
        fs->translationTable =  (_TranslationTable*)defaultTranslationTable.makeDynamic();
    }
}
//_________________________________________________________
void    processCommand (_String*s, FileState*fs)
{
    // loop thru understood values of commands
    static _List CommandList;
    if (CommandList.lLength == 0)
        // first time in, should init commands
    {
        _String command ("BASESET");
        CommandList&& & command;
        command="FORMAT";
        CommandList&& & command;
        command="RAWLINE";
        CommandList&& & command;
        command="REPEAT";
        CommandList&& & command;
        command="TOKEN";
        CommandList&& & command;
    }

    long f = -1;
    long i,k = 0,l = 0,m;
    for (i=0; (i<CommandList.lLength); i++) {
        f = s->Find (*(_String*)CommandList(i));
        if (f!=-1) {
            break;
        }
    }

    if (f==-1) { // unrecognized command
        return;
    } else {
        // trim the string
        s->Trim (f+((_String*)CommandList(i))->Length(),-1);
        f = s->Find (":");
        if (f==-1) { // poorly formed command
            return;
        } else {
            s->Trim (f+1,-1);
        }

        if ((i>=1)&&(i<=3)) {
            k = s->Find ('\"');
            if (k==-1) {
                return;
            }
            l = s->Find ('\"', k+1,-1);
            if (l==-1) {
                return;
            }
            if (l<=k) {
                return;
            }
            s->Trim (k+1,l-1);
        }

        switch (i) {
            char c;
        case 4: // new token
            checkTTStatus (fs);
            // attempt to extract a token. Looking for (e.g):   "c" = "AC"
            k = s->Find ('"');
            if (k==-1) {
                return;
            }
            if ((*s)[k+2]!='"') {
                return;
            }
            l = s->Find ('"',k+3,-1);
            m = s->Find ('"',l+1,-1);
            if ((l==-1)||(m==-1)) {
                return;
            }

            c = (*s)[k+1];
            s->Trim (l+1,m-1);
            fs->translationTable->AddTokenCode (c,*s);
            break;


        case 0:// new code set, e.g  "ACGU"
            checkTTStatus(fs);
            // erase previous char definitions
            fs->translationTable->translationsAdded.Clear();
            fs->translationTable->tokensAdded = "";
            if (*s!=_String("BASE20")) {
                fs->translationTable->AddBaseSet (*s);
            } else {
                fs->translationTable->AddBaseSet (empty);
                fs->translationTable->baseLength = 20;
            }
            break;

        case 1: //FORMAT
            if (*s==_String("PHYLIPI")) { // PHYLIP Interleaved
                fs->fileType = 1;
                fs->interleaved = TRUE;
            } else if (*s==_String("PHYLIPS")) { // PHYLIP sequential
                fs->fileType = 1;
                fs->interleaved = FALSE;
            }
            if (*s==_String("RAW")) { // RAW Sequential Data (as in NEXUS)
                fs->fileType = 2;
                fs->interleaved = FALSE;
            }
            fs->autoDetect = false;
            break;

        case 3: // REPEAT CHAR
            fs->repeat = s->getChar(0);
            break;

        case 2: // RAWLINE template e.g 1,-1 skips one word at the beginning and one word at the end
            _List chips (s,',');
            for (int i = 0; i<chips.lLength; i++) {
                fs->rawLinesFormat<<(long)(((_String*)chips(i))->toNum());
            }

        }
    }
}
//_________________________________________________________

void    FilterRawString (_String& s, FileState* fs, _DataSet & ds)
{
    int i;
    for (i = 0; i<fs->rawLinesFormat.lLength; i++) {
        long f = fs->rawLinesFormat (i),p=0,l=0;
        if (f>0) {
            for (int j = 0; (j<f)&&(p>=0)&&(l>=0); j++) {
                p = s.FirstNonSpaceIndex(l,-1,1);
                l = s.FirstSpaceIndex(p,-1,1);
            }
            if (l<0) {
                break;
            }
            p = s.FirstNonSpaceIndex(l,-1,1);
            s.Trim(p,-1);
        } else {
            if (f!=0) {
                p = 0;
                l = 0;
                for (int j = 0; (j>f)&&(p>=0)&&(l>=0); j--) {
                    p = s.FirstNonSpaceIndex(p,-1,-1);
                    l = s.FirstSpaceIndex(0,p,-1);
                }
                if (l<0) {
                    break;
                }
                p = s.FirstNonSpaceIndex(0,l,-1);
                s.Trim(0,p);
            } else {
                // Name
                p = s.FirstNonSpaceIndex();
                l = s.FirstSpaceIndex(p+1,-1,1);
                if ((p<0)||(l<0)) {
                    break;
                }
                _String Name = s.Cut (p,l-1);
                ds.AddName (Name);
                s.Trim (s.FirstNonSpaceIndex(l,-1,1),-1);
            }
        }

    }
    if (i!=fs->rawLinesFormat.lLength) {
        s = "";
    }
}
//_________________________________________________________
void    ProcessTree (FileState *fState, FILE* f, _String& CurrentLine)
{
    long j = 0, i=0; // parenthesis balance
    char c;
    _String treeString ((unsigned long)10, true);
    do {
        for (i=0; i<CurrentLine.sLength; i++) {
            c = CurrentLine.sData[i];
            if (!isspace (c)) {
                treeString<<c;
                if (c==')') {
                    j--;
                    if (!j) {
                        break;
                    }
                } else if (c=='(') {
                    j++;
                }
            }

        }
        ReadNextLine (f,&CurrentLine,fState, false);
    } while (j&&CurrentLine.sLength);

    if (j) {
        _String errMsg ("Tree string found in data file had unbalanced parentheses: ");
        if (j>0) {
            errMsg = errMsg & j & " too few closing parentheses.";
        } else {
            errMsg = errMsg & (-j) & " too many closing parentheses.";
        }
        ReportWarning (errMsg);
    } else {
        treeString.Finalize();
        setParameter (dataFileTree,1.0,fState->theNamespace);
        setParameter (dataFileTreeString, new _FString (treeString), false);
    }

}

//_________________________________________________________

long    ProcessLine (_String&s , FileState *fs, _DataSet& ds)
{
    long sitesAttached = 0,sL=s.Length();

    for (long l = 0; l<sL; l++) {
        // see if it is a legal char
        char c = toupper (s.sData[l]);
        if (fs->translationTable->IsCharLegal(c)) { // go on
            if (fs->curSpecies==0) { // add new column
                ds.AddSite (c);
                sitesAttached++;
            } else { //append to exisiting column
                //if (c == fs->skip) continue;
                // check to see if this species needs to be padded

                if (c == fs->repeat) {
                    if (fs->curSite+sitesAttached >= ds.lLength) { // a dot not matched by a previously read character; ignore
                        return sitesAttached;
                    }

                    c = ((_Site*)(ds._List::operator () (fs->curSite+sitesAttached)))->getChar(0);
                    if (c==0)
                        c = ((_Site*)(ds._List::operator ()
                                      (((_Site*)(ds._List::operator () (fs->curSite+sitesAttached)))->GetRefNo())))->getChar(0);
                }

                if (fs->curSite+sitesAttached+1>fs->totalSitesRead) {
                    // pad previous species to full length
                    _Site * newS = new _Site (fs->skip);
                    checkPointer (newS);

                    for (long j = 1; j<fs->curSpecies; j++) {
                        (*newS) << fs->skip;
                    }

                    (*newS) << c;

                    /*long rN = ds.dsh->incompletePatterns->Find (newS);

                    if (rN>=0)
                    {
                        rN =  ds.dsh->incompletePatterns->GetXtra (rN);
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
                    newS->nInstances --;

                    fs->totalSitesRead++;
                } else {
                    ds.Write2Site (fs->curSite+sitesAttached, c);
                }

                sitesAttached++;
            }
        }
    }
    // make sure that this species has enough data in it, and if not - pad it with '?'

    if ((fs->curSite+sitesAttached<fs->totalSitesRead)&&(fs->interleaved)) {
        // pad this species to full length
        for (long j = fs->curSite+sitesAttached; j<fs->totalSitesRead; j++) {
            ds.Write2Site (j, fs->skip);
        }
    }
    if (!fs->curSpecies) {
        fs->totalSitesRead+=sitesAttached;
    }
    return sitesAttached;
}

//_________________________________________________________
void    PadLine (FileState& fState, _DataSet& result) // make sure that there is enough data in this line
// and if not - "pad" it with '?''s
{
    if (fState.curSite<fState.totalSitesRead) // pad line if needed
        for (long j = fState.curSite; j<fState.totalSitesRead; j++) {
            result.Write2Site (j, fState.skip);
        }
}

//_________________________________________________________
void    ISelector (FileState& fState, _String& CurrentLine, _DataSet& result)
{
    if (fState.interleaved) { // interleaved file
        if (fState.curSpecies&&(!((fState.curSpecies)%fState.totalSpeciesExpected))) { // read a chunk of all species
            if (fState.totalSitesRead && !result.InternalStorageMode()) {
                for (long i = fState.curSite; i<fState.totalSitesRead; i++) {
                    result.Compact(i);
                }

                result.ResetIHelper();

            }
            fState.curSite = fState.totalSitesRead;
            fState.curSpecies = 0;
            ProcessLine (CurrentLine, &fState, result);
            fState.curSpecies = 1;
            if (!fState.curSite) {
                fState.totalSpeciesRead++;
            }
        } else {
            ProcessLine (CurrentLine, &fState, result);
            if (!fState.curSite) {
                fState.totalSpeciesRead++;
            }
            fState.curSpecies++;
        }
    } else {
        if (fState.curSpecies+1<fState.totalSpeciesExpected) {
            fState.curSpecies++;
        }
        if (fState.curSpecies == fState.totalSpeciesRead) {
            PadLine (fState, result);
            fState.curSite = 0;
        }
        if (fState.totalSpeciesRead<fState.totalSpeciesExpected) {
            fState.totalSpeciesRead++;
        }

        fState.curSite+=ProcessLine (CurrentLine, &fState, result);

    }
}

//_________________________________________________________
bool SkipLine (_String& theLine, FileState* fS)
{
    if ( theLine.sData[0]=='/' && theLine.sData[1]=='/' ) {
        return true;
    }

    char c = theLine.FirstNonSpace();

    if (c&& (!((c=='$')&&(!fS->acceptingCommands))) ) {
        return false;
    }

    return true;
}

//_________________________________________________________

#define     READ_NEXT_LINE_BUFFER_SIZE      1024*1024


//_________________________________________________________
void ReadNextLine (FILE* fp, _String *s, FileState* fs, bool, bool upCase)
{
    _String  tempBuffer (1024L, true);

    char lastc;

    if (fp) {
        lastc = fgetc(fp);
    } else {
        lastc = fs->pInSrc<fs->theSource->sLength?fs->theSource->sData[fs->pInSrc++]:0;
    }

    if (fs->fileType != 3) { // not NEXUS - do not skip [..]
        if (fp)
            while ( !feof(fp) && lastc!=10 && lastc!=13 ) {
                if (lastc) {
                    tempBuffer << lastc;
                }

                lastc = fgetc(fp);
            }
        else
            while (lastc && lastc!=10 && lastc!=13 ) {
                tempBuffer << lastc;
                lastc = fs->theSource->sData[fs->pInSrc++];
            }

    } else {
        if (upCase) {
            lastc = toupper(lastc);
        }

        while (((fp&&(!feof(fp)))||(fs->theSource&&(fs->pInSrc<=fs->theSource->sLength))) && lastc!=10 && lastc!=13) {
            if (lastc=='[') {
                if (fs->isSkippingInNEXUS) {
                    ReportWarning ("Nested comments in NEXUS really shouldn't be used.");
                } else {
                    fs->isSkippingInNEXUS = true;
                }
            }
            if (fs->isSkippingInNEXUS) {
                if (lastc==']') {
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

        if ( lastc==10 || lastc==13 ) {
            tempBuffer << ' ';
        }
    }

    tempBuffer.Finalize();

    if ( fp&&feof(fp) || fs->theSource&& fs->pInSrc>=fs->theSource->sLength )
        if (tempBuffer.sLength == 0) {
            *s = "";
            return;
        }

    if (s->nInstances > 1) {
        *s = tempBuffer;
    } else {
        Ptr         saveData = s->sData;
        s->sData    = tempBuffer.sData;
        tempBuffer.sData = saveData;

        s->sLength  = tempBuffer.sLength;
    }

    if (SkipLine (*s, fs)) {
        ReadNextLine(fp,s,fs,false,upCase);
    }

    if (s->sLength && s->sData[s->sLength-1]== '\n') {
        s->Trim (0,s->sLength-2);
    }
}
//_________________________________________________________
void    TrimPhylipLine (_String& CurrentLine, _DataSet& ds)
{
    int  fNS = CurrentLine.FirstNonSpaceIndex();
    _String     Name (CurrentLine.Cut (fNS, fNS+9));
    CurrentLine.Trim(fNS+10,-1); // chop out the name
    ds.AddName(Name);
}


//_________________________________________________________
_DataSet* ReadDataSetFile (FILE*f, char execBF, _String* theS, _String* bfName, _String* namespaceID, _TranslationTable* dT)
{

    bool     doAlphaConsistencyCheck = true;
    _String::storageIncrement = 16;
    _DataSet* result = new _DataSet;
    fileTreeString = "";


    _String         CurrentLine = dataFilePartitionMatrix & "={{}};",
                    savedLine;

    if (1) {
        _ExecutionList reset (CurrentLine);
        reset.Execute();
#ifdef __HYPHYMPI__
        if (_hy_mpi_node_rank == 0)
#endif
            terminateExecution = false;
    }

    // initialize the instance of a file state variable
    setParameter(dataFileTree, 0.0);
    FileState   fState;
    fState.translationTable =  dT;
    fState.curSpecies =
        fState.totalSpeciesRead =
            fState.totalSitesRead =
                fState.totalSpeciesExpected =
                    fState.totalSitesExpected =
                        fState.curSite =
                            fState.maxStringLength   = 0;
    fState.acceptingCommands = true;
    fState.allSpeciesDefined = false;
    fState.interleaved       = false;
    fState.isSkippingInNEXUS = false;
    fState.autoDetect        = true;
    fState.fileType          = -1;
    fState.baseLength        = 4;
    fState.repeat            = '.',
           fState.skip            = 0;
    fState.theSource         = theS;
    fState.pInSrc            = 0;
    fState.theNamespace      = namespaceID;

    if (!(f||theS)) {
        CurrentLine = "ReadDataSetFile received null file AND string references. At least one must be specified";
        warnError     (CurrentLine);
    }
    // done initializing

    long     fileLength = 0,
             lastDone = 10,
             cDone;

#ifdef __HYPHYMPI__
    if (_hy_mpi_node_rank == 0) {
#endif
        if       (f) {
            fseek    (f,0,SEEK_END);
            fileLength = ftell(f);
            rewind  (f);
        } else {
            fileLength = theS->sLength;
        }

#ifdef __HYPHYMPI__
    }
#endif



    //if (f==NULL) return (_DataSet*)result.makeDynamic();
    // nothing to do

    CurrentLine = empty;

    ReadNextLine (f,&CurrentLine,&fState);
    if (!CurrentLine.sLength) {
        CurrentLine = "Empty File Encountered By ReadDataSet.";
        WarnError (CurrentLine);
        return result;
    } else {
        if (CurrentLine.beginswith ("#NEXUS",false)) {
            ReadNexusFile (fState,f,(*result));
            doAlphaConsistencyCheck = false;
        } else {
            long i,j,k, filePosition = -1, saveSpecExpected;
            char c;
            while (CurrentLine.sLength) { // stuff to do
                // check if the line has a command in it
#ifdef __HYPHYMPI__
                if (_hy_mpi_node_rank == 0) {
#endif
                    if (f) {
                        cDone = ftell (f)*100./fileLength;
                    } else {
                        cDone = fState.pInSrc*100./fileLength;
                    }

                    if (cDone>lastDone) {
                        SetStatusBarValue (lastDone,1,0);
#ifdef __MAC__
                        handleGUI(true);
#endif
                        lastDone+=10;
                    }
#ifdef __HYPHYMPI__
                }
#endif

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

                    if (fState.fileType==-1) { // undecided file type - assume it is PHYLIP sequential
                        if ((c == '#')||(c=='>')) { // hash-mark format
                            fState.fileType = 0;
                        } else { // assume this is a sequential PHYLIP file
                            fState.fileType = 1;
                            fState.interleaved = false;
                        }

                    }
                    // decide what to do next
                    // if format is PHYLIP and we do not know the expected dimensions,
                    //   we must read those in first
                    if (fState.fileType==1) { // PHYLIP
                        if ((filePosition<0)&&(fState.autoDetect)) {
                            filePosition = (f?
                                            ftell (f)
#ifdef __WINDOZE__
                                            -1
#endif
                                            :fState.pInSrc);
                            savedLine = CurrentLine;
                        }

                        if ((fState.totalSitesExpected==0)||(fState.totalSpeciesExpected==0)) { // must read dimensions first
                            i = CurrentLine.FirstNonSpaceIndex();
                            j = CurrentLine.FirstSpaceIndex(i,-1,1);
                            if (j>=0) {
                                k = CurrentLine.FirstNonSpaceIndex(j,-1,1);
                                if (k>=0) { // could have dimensions
                                    saveSpecExpected = fState.totalSpeciesExpected=CurrentLine.Cut(i,j-1).toNum();
                                    fState.totalSitesExpected=CurrentLine.Cut(k,-1).toNum();
                                }
                                if (CurrentLine.Find ('I', k, -1)>=0) { // interleaved
                                    fState.interleaved = true;
                                }
                            }
                        } else
                            // now for the data crunching part
                            // detect a line, diagnose it and dispatch accordingly
                        {
                            if (fState.interleaved) {
                                if (fState.totalSpeciesRead<fState.totalSpeciesExpected) {
                                    TrimPhylipLine (CurrentLine, (*result));
                                }
                                if ((fState.curSite)&&(fState.curSpecies >= saveSpecExpected)&&
                                        (fState.totalSitesRead >= fState.totalSitesExpected)) {
                                    // reached the end of the data - see maybe there is a tree
                                    ReadNextLine (f,&CurrentLine,&fState);
                                    if (CurrentLine.sLength) {
                                        if (CurrentLine.FirstNonSpace()=='(') { // could be a tree string
                                            ProcessTree (&fState,f, CurrentLine);
                                        }
                                    }
                                    break;
                                }

                            } else {
                                if (fState.totalSitesRead > fState.totalSitesExpected)
                                    // oops - autodetect incorrectly assumed that the file was sequential
                                {
                                    fState.curSpecies =
                                        fState.totalSpeciesRead =
                                            fState.totalSitesRead =
                                                fState.curSite =
                                                    fState.totalSpeciesExpected =
                                                        fState.totalSitesExpected =
                                                            fState.maxStringLength = 0;
                                    fState.allSpeciesDefined = false;
                                    fState.interleaved = true;
                                    fState.autoDetect = true;

                                    if(f) {
                                        fseek (f, filePosition, SEEK_SET);
                                    } else {
                                        fState.pInSrc = filePosition;
                                    }

                                    CurrentLine = savedLine;
                                    for (long idx = 0; idx < (*result).lLength; idx++) {
                                        ((_Site*)(*result).lData[idx])->Finalize();
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
                                if (fState.totalSpeciesRead==0) {
                                    fState.totalSpeciesExpected = 1;
                                    if (!fState.curSite) {
                                        TrimPhylipLine (CurrentLine, (*result));
                                    }
                                }

                                else if (fState.curSite>=fState.totalSitesExpected) {
                                    fState.totalSpeciesExpected++;
                                    if (fState.totalSpeciesExpected>saveSpecExpected) {
                                        // reached the end of the data - see maybe there is a tree
                                        ReadNextLine (f,&CurrentLine,&fState);
                                        if (CurrentLine.sLength) {
                                            if (CurrentLine.FirstNonSpace()=='(') { // could be a tree string
                                                ProcessTree (&fState,f, CurrentLine);
                                            }
                                        }
                                        break;
                                    }
                                    TrimPhylipLine (CurrentLine, (*result));
                                }
                            }

                            ISelector (fState, CurrentLine, (*result));
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
                        if (ProcessLine (CurrentLine, &fState, (*result))) {
                            fState.curSpecies++;
                            fState.totalSpeciesRead++;
                        }
                        break;
                    }

                    // lastly, handle the auto-detect standard case

                    // check to see if the string defines a name
                    if ((c=='#')||(c=='>')) { // a name it is
                        if (fState.allSpeciesDefined) { // can't define the species after data
                            break;
                        } else {
                            if ((!fState.totalSpeciesRead)&&(fState.totalSpeciesExpected>=1)) {
                                fState.interleaved = TRUE;
                            } else {
                                fState.interleaved = FALSE;
                            }
                            fState.totalSpeciesExpected++;
                            CurrentLine.Trim(CurrentLine.FirstNonSpaceIndex(1,-1,1),-1);
                            if ((CurrentLine.sData[0]=='#')||(CurrentLine.sData[0]=='>')) {
                                CurrentLine = _String("Species")&_String(fState.totalSpeciesExpected);
                            }
                            (*result).AddName (CurrentLine);
                        }
                        break;
                    }
                    // check to see if the string defines a tree
                    if (c=='(') {
                        ProcessTree (&fState,f, CurrentLine);
                        ReadNextLine (f,&CurrentLine,&fState);
                    }

                    // check to see where to stick the incoming line

                    if (!fState.totalSpeciesExpected)
                        // raw data fed before names defined - skip
                    {
                        break;
                    }
                    if((fState.totalSpeciesExpected>1)&&(!fState.totalSpeciesRead)) {
                        fState.allSpeciesDefined = TRUE;
                    }

                    // repeat the structure of PHYLIP reader

                    ISelector (fState, CurrentLine, (*result));

                    break;
                }

                ReadNextLine (f,&CurrentLine,&fState);

            }
        }
    }



    if (fState.totalSitesRead && fState.interleaved && !result->InternalStorageMode()) {
        for (int i = fState.curSite; i<fState.totalSitesRead; i++) {
            (*result).Compact(i);
        }
        (*result).ResetIHelper();
    }

    if ((!fState.interleaved)&&(fState.fileType!=2)) {
        PadLine (fState, (*result));
    }


#ifdef __HYPHYMPI__
    if (_hy_mpi_node_rank == 0) {
#endif
        SetStatusBarValue (-1,1,0);
#ifdef __MAC__
        handleGUI(true);
#endif
#ifdef __HYPHYMPI__
    }
#endif

    // make sure interleaved duplications are handled correctly

    (*result).Finalize();
    _String::storageIncrement   = 64;
    (*result).noOfSpecies       = fState.totalSpeciesRead;
    (*result).theTT             = fState.translationTable;

    // check to see if result may be an amino-acid data
    if (doAlphaConsistencyCheck && result->theTT == &defaultTranslationTable) {
        if (result->GetNoTypes() == 0)
            // empty data set
            // try binary data
        {
            _TranslationTable *trialTable = new _TranslationTable (defaultTranslationTable);
            trialTable->baseLength = 2;
            _DataSet * res2 = ReadDataSetFile (f, execBF, theS, bfName, namespaceID, trialTable);
            if (res2->GetNoTypes()) {
                DeleteObject (result);
                return res2;
            }
            DeleteObject (trialTable);
        } else
            // check it out
            if (result->CheckAlphabetConsistency()<0.5)
                // less than 50% of the data in the alphabet is not in the basic alphabet
            {
                _TranslationTable trialTable (defaultTranslationTable);
                trialTable.baseLength = 20;
                (*result).theTT = &trialTable;
                if ((*result).CheckAlphabetConsistency()<0.5) {
                    CurrentLine = "More than 50% of characters in the data are not in the alphabet.";
                    (*result).theTT =  &defaultTranslationTable;
                    ReportWarning (CurrentLine);
                } else {
                    (*result).theTT = (_TranslationTable*)trialTable.makeDynamic();
                }

            }

    }
    if (nexusBFBody.sLength) {
        if (execBF == 1) {
            lastNexusDataMatrix = result;

            long            bfl = batchLanguageFunctions.lLength;

            _ExecutionList nexusBF (nexusBFBody,namespaceID);
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

            //DeleteObject (lastNexusDataMatrix);
            lastNexusDataMatrix = nil;
            nexusBFBody         = empty;
        } else if (execBF == 0) {
            nexusBFBody         = empty;
        }
    }

    //return (_DataSet*)result.makeDynamic();
    return result;
}

//_________________________________________________________

BaseRef _DataSetFilter::toStr (void)
{
    //return new _String("DataSetFilters only print to files");
    _String * res = new _String (4096L, true);
    checkPointer (res);
    internalToStr (nil,*res);
    res->Finalize();
    return res;
}

//_________________________________________________________

void    _DataSetFilter::PatternToSiteMapper (void* source, void* target, char mode, long padup)
{
    for (long site = 0; site < duplicateMap.lLength; site++)
        if (mode == 0) {
            ((_Parameter*)target)[site] = ((_Parameter*)source)[duplicateMap.lData[site]];
        } else if (mode == 1) {
            ((long*)target)[site] = ((long*)source)[duplicateMap.lData[site]];
        } else if (mode == 2) {
            ((long*)target)[site] = ((_Parameter*)source)[duplicateMap.lData[site]];
        }


    for (long site = duplicateMap.lLength; site < padup; site++)
        if (mode == 0) {
            ((_Parameter*)target)[site] = 1.;
        } else if (mode == 1) {
            ((long*)target)[site] = 0;
        }
}


//_________________________________________________________

long    _DataSetFilter::GetOriginalToShortMap(long index)
{
    long pos1=theData->theMap.lData[theOriginalOrder.lData[index]],pos2;
    pos2 = theMap.Find(pos1);
    if (pos2==-1) {
        for (long i=theData->theMap.lLength-1; i>=0; i--) {
            if (theData->theMap.lData[i]==pos1) {
                pos2 = theMap.Find(i);
                if (pos2!=-1) {
                    break;
                }
            }
        }
    }
    return pos2;
}

//_________________________________________________________

_String _DataSetFilter::GenerateConsensusString (_SimpleList* majority)
{
    if (unitLength > 3) {
        return empty;
    }

    _String     result ((unsigned long)theOriginalOrder.lLength),
                tRes   ((unsigned long)(unitLength*theFrequencies.lLength));

    long        charStates         = GetDimension(false),
                *translationBuffer = (long*)MemAllocate(sizeof(long)*charStates);

    _Parameter* countBuffer = (_Parameter*)MemAllocate(sizeof(_Parameter)*charStates),
                nf;

    SetupConversion ();

    for (long k=0; k<theFrequencies.lLength; k++) {
        long    m,
                t = theMap.lData[k],
                f;

        for (m=0; m<charStates; m++) {
            countBuffer[m] = 0.0;
        }

        for (long p=0; p<theNodeMap.lLength; p++) {
            theData->theTT->TokenCode ((*theData)(t, theNodeMap.lData[p],1),translationBuffer);
            f = 0;

            for (m=0; m<charStates; m++)
                if (translationBuffer[m]) {
                    f++;
                }

            if (f>1) {
                nf = 1./f;
                for (m=0; m<charStates; m++)
                    if (translationBuffer[m]) {
                        countBuffer[m]+=nf;
                    }
            } else {
                if (f==1) {
                    m=0;
                    while (!translationBuffer[m++]) ;
                    countBuffer[m-1]+=1.;
                }
            }
        }

        nf = -1;
        f  =  1;

        for (m=0; m<charStates; m++) {
            if (countBuffer[m]>nf) {
                nf = countBuffer[m];
                t = m;
                f = 1;
            } else if (countBuffer[m]==nf) {
                f++;
            }
        }

        if (f>1) {
            for (m=0; m<charStates; m++) {
                if (countBuffer[m]==nf) {
                    translationBuffer[m] = 1;
                } else {
                    translationBuffer[m] = 0;
                }
            }
            tRes.sData[k]=theData->theTT->CodeToLetter(translationBuffer);
            if (majority)
                //(*majority) << -1;
            {
                (*majority) << nf;
            }
        } else {
            _String conv = theData->theTT->ConvertCodeToLetters(t,1);
            tRes.sData[k] = conv.sData[0];
            if (majority) {
                (*majority) << nf;
            }
        }
    }

    free (countBuffer);
    free (translationBuffer);

    for (long m=0; m<theOriginalOrder.lLength; m++) {
        result.sData[m] = tRes.sData[duplicateMap.lData[m]];
    }

    return result;
}


//_________________________________________________________
void    _DataSetFilter::toFileStr (FILE*dest)
{
// write out the file with this dataset filter
    if (!dest) {
        return;
    }

    _String       dummy;
    internalToStr (dest,dummy);
}

//_________________________________________________________
void    _DataSetFilter::ConvertCodeToLettersBuffered (long code, char unit, char* storage, _AVLListXL* lookup)
{
    // write out the file with this dataset filter
    long      lookupC     = lookup->Find ((BaseRef)code);
    char      *lookupV;
    if (lookupC>=0) {
        lookupV = ((_String*)lookup->GetXtra(lookupC))->sData;
    } else {
        _String * newT = new _String (ConvertCodeToLetters (code,unit));
        lookup->Insert ((BaseRef)code, (long)newT, false);
        lookupV = newT->sData;
    }

    for (long k = 0; k < unit; k++) {
        storage[k] = lookupV[k];
    }
}


//_________________________________________________________

void    _DataSetFilter::internalToStr (FILE*dest,_String& rec)
{
// write out the file with this dataset filter
    checkParameter (dataFilePrintFormat,dFPrintFormat,6.0);
    checkParameter (dataFileDefaultWidth,dFDefaultWidth,50.0);
    _Parameter  gW;

    long outputFormat = dFPrintFormat,
         printWidth   = dFDefaultWidth,
         gapWidth;

    checkParameter (dataFileGapWidth,gW,10.0);
    if(!printWidth) {
        printWidth = 50;
    }

    gapWidth = gW;
    if (gapWidth<=0) {
        gapWidth = printWidth;
    }

    long i,
         j;

    if (outputFormat < 4 || outputFormat > 8)
        if (!(theData->theTT->IsStandardNucleotide() || theData->theTT->IsStandardAA())) {
            _String * bSet = &theData->theTT->baseSet;
            if (dest) {
                fprintf (dest,"$BASESET:\"%s\"\n",bSet->sData);
                if (theData->theTT->tokensAdded.sLength) {
                    for (long at = 0; at < theData->theTT->tokensAdded.sLength; at++) {
                        fprintf (dest, "$TOKEN:\"%c\" = \"", theData->theTT->tokensAdded.sData[at]);
                        long    buf [256];
                        theData->theTT->SplitTokenCode(theData->theTT->TokenCode(theData->theTT->tokensAdded.sData[at]), buf);
                        for (long tc = 0; tc < bSet->sLength; tc++)
                            if (buf[tc]) {
                                fprintf (dest, "%c", bSet->sData[tc]);
                            }
                        fprintf (dest, "\"\n");
                    }
                }
            } else {
                rec << "$BASESET:\"";
                rec << *bSet;
                rec << "\"\n";
                if (theData->theTT->tokensAdded.sLength) {
                    for (long at = 0; at < theData->theTT->tokensAdded.sLength; at++) {
                        rec << "$TOKEN:\"";
                        rec << theData->theTT->tokensAdded.sData[at];
                        rec << "\" = \"";
                        long    buf [256];
                        theData->theTT->SplitTokenCode(theData->theTT->TokenCode(theData->theTT->tokensAdded.sData[at]), buf);
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
    case 1: // hash-mark interleaved
    case 10: { // FASTA interleaved
        long sitesDone    = 0,
             upTo;
        char seqDelimiter = (outputFormat==1)?'#':'>';

        for ( i = 0; i<theNodeMap.lLength; i++) {
            _String * curName = (_String *)theData->GetNames() (theNodeMap.lData[i]);
            if (dest) {
                fprintf (dest, "%c%s\n", seqDelimiter,curName->sData);
            } else {
                rec << *curName;
                rec << '\n';
            }
        }
        while (sitesDone<theOriginalOrder.lLength) {
            if (dest) {
                fprintf (dest,"\n\n");
            } else {
                rec << '\n';
                rec << '\n';
            }

            upTo = sitesDone+printWidth;
            if (upTo>theOriginalOrder.lLength) {
                upTo = theOriginalOrder.lLength;
            }

            if (dest)
                for ( i = 0; i<theNodeMap.lLength; i++) {
                    for ( j = sitesDone; j<upTo; j++) {
                        if ((j-sitesDone)%gapWidth==0) {
                            fprintf (dest, " ");
                        }
                        fprintf (dest, "%c",(*theData)(theOriginalOrder.lData[j],theNodeMap.lData[i],1));
                    }

                    fprintf (dest, "\n");
                }
            else
                for ( i = 0; i<theNodeMap.lLength; i++) {
                    for ( j = sitesDone; j<upTo; j++) {
                        if ((j-sitesDone)%gapWidth==0) {
                            rec << ' ';
                        }
                        rec << (*theData)(theOriginalOrder.lData[j],theNodeMap.lData[i],1);
                    }

                    rec << '\n';
                }

            sitesDone = upTo;
        }
        break;
    }

    case 2: { // phylip sequential
        // print PHYLIP format header
        //fprintf (dest,"$FORMAT:\"PHYLIPS\"\n");
        // print number of species and sites
        if (dest) {
            fprintf (dest,"%ld\t%ld\n",theNodeMap.lLength,theOriginalOrder.lLength);
        } else {
            rec << _String((long)theNodeMap.lLength);
            rec << '\t';
            rec << _String((long)theNodeMap.lLength,theOriginalOrder.lLength);
            rec << '\n';
        }
        // proceed to spool out the data
        for ( i = 0; i<theNodeMap.lLength; i++) {
            _String * curName = (_String *)theData->GetNames() (theNodeMap(i)), choppedTo10Chars;
            if (curName->Length()>=10) {
                choppedTo10Chars = curName->Cut(0,9)&' ';
            } else {
                choppedTo10Chars = *curName;
                while (choppedTo10Chars.Length()<11) {
                    choppedTo10Chars=choppedTo10Chars&' ';
                }
            }

            if (dest) {
                fprintf (dest, "%s",choppedTo10Chars.sData);

                for ( j = 0; j<theOriginalOrder.lLength; j++) {
                    if ((j%printWidth==0)&&j) {
                        fprintf(dest,"\n           ");
                    }
                    fprintf (dest, "%c",(*theData)(theOriginalOrder(j),theNodeMap(i),1));
                    if (j%gapWidth==gapWidth-1) {
                        fprintf (dest, " ");
                    }
                }
                fprintf (dest, "\n");
            } else {
                rec << choppedTo10Chars;

                for ( j = 0; j<theOriginalOrder.lLength; j++) {
                    if ((j%printWidth==0)&&j) {
                        rec << "\n           ";
                    }

                    rec << (*theData)(theOriginalOrder(j),theNodeMap(i),1);
                    if (j%gapWidth==gapWidth-1) {
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
            fprintf (dest,"%ld\t%ld\n",theNodeMap.lLength,theOriginalOrder.lLength);
        } else {
            rec << _String((long)theNodeMap.lLength);
            rec << '\t';
            rec << _String((long)theNodeMap.lLength,theOriginalOrder.lLength);
            rec << '\n';
        }
        // proceed to spool out the data
        for ( i = 0; i<theNodeMap.lLength; i++) {
            _String * curName = (_String *)theData->GetNames() (theNodeMap(i)), choppedTo10Chars;
            if (curName->Length()>=10) {
                choppedTo10Chars = curName->Cut(0,9)&' ';
            } else {
                choppedTo10Chars = *curName;
                while (choppedTo10Chars.Length()<11) {
                    choppedTo10Chars=choppedTo10Chars&' ';
                }
            }

            if (dest) {
                fprintf (dest, "%s",choppedTo10Chars.sData);

                for ( j = 0; j<theOriginalOrder.lLength; j++) {
                    if (j==printWidth) {
                        fprintf(dest,"\n");
                        break;
                    }
                    if (j%gapWidth==0) {
                        fprintf (dest, " ");
                    }
                    fprintf (dest, "%c",(*theData)(theOriginalOrder.lData[j],theNodeMap.lData[i],1));
                }
            } else {
                rec << choppedTo10Chars;

                for ( j = 0; j<theOriginalOrder.lLength; j++) {
                    if (j==printWidth) {
                        rec << '\n';
                        break;
                    }
                    if (j%gapWidth==0) {
                        rec << ' ';
                    }
                    rec << (*theData)(theOriginalOrder.lData[j],theNodeMap.lData[i],1);
                }
            }

        }

        long completed = printWidth;

        if (dest) {
            while (completed<theOriginalOrder.lLength-1) {
                long upTo = completed+printWidth<theOriginalOrder.lLength?completed+printWidth:theOriginalOrder.lLength;
                for ( i = 0; i<theNodeMap.lLength; i++) {
                    fprintf(dest,"\n           ");
                    for ( j = completed; j<upTo; j++) {
                        if ((j-completed)%gapWidth==0) {
                            fprintf (dest, " ");
                        }
                        fprintf (dest, "%c",(*theData)(theOriginalOrder.lData[j],theNodeMap.lData[i],1));
                    }
                }
                completed+=printWidth;
                fprintf (dest,"\n");
            }
        } else {
            while (completed<theOriginalOrder.lLength-1) {
                long upTo = completed+printWidth<theOriginalOrder.lLength?completed+printWidth:theOriginalOrder.lLength;
                for ( i = 0; i<theNodeMap.lLength; i++) {
                    rec << "\n           ";
                    for ( j = completed; j<upTo; j++) {
                        if ((j-completed)%gapWidth==0) {
                            rec <<  ' ';
                        }
                        rec << (*theData)(theOriginalOrder.lData[j],theNodeMap.lData[i],1);
                    }
                }
                completed+=printWidth;
                rec << '\n';
            }
        }

        break;
    }

    // various flavors of NEXUS

    case 4: // labels, sequential
    case 5: // labels, interleaved
    case 6: // no labels, sequential
    case 7: { // no labels, interleaved
        // write out the header
        j = theNodeMap.lLength;
        if (dest) {
            fprintf (dest, "#NEXUS\n\n[\nGenerated by %s on %s\n]\n\nBEGIN TAXA;\n\tDIMENSIONS NTAX = %ld;\n\tTAXLABELS\n\t\t", GetVersionString().getStr() , GetTimeStamp().getStr(), j);
            for (i=0; i<j; i++) {
                fprintf (dest, "'%s' ", ((_String*)theData->GetNames().lData[theNodeMap.lData[i]])->sData);
            }

            fprintf (dest,";\nEND;\n\nBEGIN CHARACTERS;\n\tDIMENSIONS NCHAR = %ld;\n\tFORMAT\n\t\t",theOriginalOrder.lLength);
            if (theData->theTT->IsStandardNucleotide()) {
                fprintf (dest,"DATATYPE = DNA\n");
            } else {
                if (theData->theTT->IsStandardAA()) {
                    fprintf (dest,"DATATYPE = PROTEIN\n");
                } else if (theData->theTT->IsStandardBinary()) {
                    fprintf (dest,"DATATYPE = BINARY\n");
                } else {
                    _String * bSet = &theData->theTT->baseSet;
                    fprintf (dest, "\n\tSYMBOLS = \"");
                    for (long bc = 0; bc < bSet->sLength-1; bc++) {
                        fprintf (dest, "%c ", bSet->sData[bc]);
                    }
                    fprintf (dest, "%c\"\n", bSet->sData[bSet->sLength-1]);
                    if (theData->theTT->tokensAdded.sLength) {
                        for (long at = 0; at < theData->theTT->tokensAdded.sLength; at++) {
                            fprintf (dest,"\n\tEQUATE=\"%c = ", theData->theTT->tokensAdded.sData[at]);
                            long    buf [256];
                            theData->theTT->SplitTokenCode(theData->theTT->TokenCode(theData->theTT->tokensAdded.sData[at]), buf);
                            for (long tc = 0; tc < bSet->sLength; tc++)
                                if (buf[tc]) {
                                    fprintf (dest, "%c", bSet->sData[tc]);
                                }
                            fprintf (dest, "\"");
                        }
                    }
                }
            }
            if (theData->theTT->GetGapChar()) {
                fprintf (dest,"\n\t\tGAP=%c",theData->theTT->GetGapChar());
            }
            if (theData->theTT->GetSkipChar()) {
                fprintf (dest,"\n\t\tMISSING=%c",theData->theTT->GetSkipChar());
            }


            if (outputFormat>5) {
                fprintf (dest,"\n\t\tNOLABELS");
            }
            if (outputFormat%2) {
                fprintf (dest,"\n\t\tINTERLEAVE");
            }

            fprintf (dest,"\n\t;\n\nMATRIX");
        } else {
            rec << "#NEXUS\n\nBEGIN TAXA;\n\tDIMENSIONS NTAX = ";
            rec << _String ((long)j);
            rec << ";\n\tTAXLABELS\n\t\t";

            for (i=0; i<j; i++) {
                rec << "'";
                rec << ((_String*)theData->GetNames().lData[theNodeMap.lData[i]]);
                rec <<  "' ";
            }
            rec << ";\nEND;\n\nBEGIN CHARACTERS;\n\tDIMENSIONS NCHAR = ";
            rec << _String((long)theOriginalOrder.lLength);
            rec << ";\n\tFORMAT\n\t\t";

            if (theData->theTT->IsStandardNucleotide()) {
                rec << "DATATYPE = DNA\n";
            } else {
                if (theData->theTT->IsStandardAA()) {
                    rec << "DATATYPE = PROTEIN\n";
                } else if (theData->theTT->IsStandardBinary()) {
                    rec << "DATATYPE = BINARY\n";
                } else {
                    rec << "\t\tSYMBOLS = \"";
                    _String * bSet = &theData->theTT->baseSet;
                    for (long bc = 0; bc < bSet->sLength-1; bc++) {
                        rec << bSet->sData[bc];
                        rec << ' ';
                    }
                    rec << bSet->sData[bSet->sLength-1];
                    rec << "\"\n";
                    if (theData->theTT->tokensAdded.sLength)
                        for (long at = 0; at < theData->theTT->tokensAdded.sLength; at++) {
                            rec << "\nEQUATE =\"";
                            rec << theData->theTT->tokensAdded.sData[at];
                            rec << " = ";
                            long    buf [256];
                            theData->theTT->SplitTokenCode(theData->theTT->TokenCode(theData->theTT->tokensAdded.sData[at]), buf);
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
            if (outputFormat>5) {
                rec << "\n\t\tNOLABELS";
            }
            if (outputFormat%2) {
                rec << "\n\t\tINTERLEAVE";
            }

            rec << "\n\t;\n\nMATRIX";

        }

        //compute space alignment for different taxa names
        // two passes - one to locate the max length and 2nd to compute padding lengths

        j = 0;
        for (i=0; i<theNodeMap.lLength; i++) {
            if (((_String *)theData->GetNames() (theNodeMap(i)))->sLength>j) {
                j = ((_String *)theData->GetNames() (theNodeMap(i)))->sLength;
            }
        }

        _SimpleList taxaNamesPadding;

        for (i=0; i<theNodeMap.lLength; i++) {
            taxaNamesPadding<<(j - ((_String *)theData->GetNames() (theNodeMap(i)))->sLength);
        }

        if (outputFormat%2==0) { // sequential
            if (dest)
                for (i=0; i<theNodeMap.lLength; i++) {
                    if (outputFormat == 4) { // labels
                        fprintf (dest,"\n\t'%s'",((_String *)theData->GetNames() (theNodeMap(i)))->sData);
                        for (j=0; j<=taxaNamesPadding.lData[i]; j++) {
                            fprintf (dest," ");
                        }
                    } else {
                        fprintf (dest, "\n");
                    }
                    fprintf (dest," ");
                    for ( j = 0; j<theOriginalOrder.lLength; j++) {
                        fprintf (dest, "%c",(*theData)(theOriginalOrder.lData[j],theNodeMap.lData[i],1));
                    }
                }
            else
                for (i=0; i<theNodeMap.lLength; i++) {
                    if (outputFormat == 4) { // labels
                        rec << "\n\t'";
                        rec << (*(_String *)theData->GetNames() (theNodeMap(i)));
                        rec << "'";

                        for (j=0; j<=taxaNamesPadding.lData[i]; j++) {
                            rec << ' ';
                        }
                    } else {
                        rec << '\n';
                    }
                    rec << ' ';
                    for ( j = 0; j<theOriginalOrder.lLength; j++) {
                        rec << (*theData)(theOriginalOrder.lData[j],theNodeMap.lData[i],1);
                    }
                }
        } else {
            long  sitesDone = 0, upTo;

            while (sitesDone<theOriginalOrder.lLength) {
                upTo = sitesDone+printWidth;
                if (upTo>theOriginalOrder.lLength) {
                    upTo = theOriginalOrder.lLength;
                }


                if (dest) {
                    for (i = 0; i<theNodeMap.lLength; i++) {
                        if (outputFormat == 5) { // labels
                            fprintf (dest,"\n\t'%s'",((_String *)theData->GetNames() (theNodeMap(i)))->sData);
                            for (j=0; j<=taxaNamesPadding.lData[i]; j++) {
                                fprintf (dest," ");
                            }
                        } else {
                            fprintf (dest,"\n");
                        }
                        fprintf (dest," ");
                        for ( j = sitesDone; j<upTo; j++) {
                            fprintf (dest, "%c",(*theData)(theOriginalOrder.lData[j],theNodeMap.lData[i],1));
                        }
                    }
                    fprintf (dest,"\n\n");
                } else {
                    for (i = 0; i<theNodeMap.lLength; i++) {
                        if (outputFormat == 5) { // labels
                            rec << "\n\t'";
                            rec << (*(_String *)theData->GetNames() (theNodeMap(i)));
                            rec << "'";
                            for (j=0; j<=taxaNamesPadding.lData[i]; j++) {
                                rec << ' ';
                            }
                        } else {
                            rec << '\n';
                        }

                        rec << ' ';
                        for ( j = sitesDone; j<upTo; j++) {
                            rec << (*theData)(theOriginalOrder.lData[j],theNodeMap.lData[i],1);
                        }
                    }
                    rec << '\n';
                    rec << '\n';
                }

                sitesDone = upTo;
            }

        }
        if (dest) {
            fprintf (dest,";\nEND;");
        } else {
            rec << ";\nEND;";
        }
        break;
    }

    case 8: {
        for ( i = 0; i<theNodeMap.lLength; i++) {
            if (dest) {
                fprintf (dest, "%c",(*theData)(theOriginalOrder(0),theNodeMap(i),1));
                for ( j = 1; j<theOriginalOrder.lLength; j++) {
                    fprintf (dest, ",%c",(*theData)(theOriginalOrder(j),theNodeMap(i),1));
                }
                fprintf (dest, "\n");
            } else {
                rec << (*theData)(theOriginalOrder(0),theNodeMap(i),1);
                for ( j = 1; j<theOriginalOrder.lLength; j++) {
                    rec << ',';
                    rec << (*theData)(theOriginalOrder(j),theNodeMap(i),1);
                }
                rec << '\n';
            }
        }
        break;
    }

    default: { // hash-mark sequential
        char seqDelimiter = (outputFormat==9)?'>':'#';

        for ( i = 0; i<theNodeMap.lLength; i++) {
            _String * curName = (_String *)theData->GetNames() (theNodeMap(i));

            if (dest) {
                fprintf (dest, "%c%s", seqDelimiter ,curName->sData);
                for ( j = 0; j<theOriginalOrder.lLength; j++) {
                    if (j%printWidth==0) {
                        fprintf(dest,"\n");
                    }
                    fprintf (dest, "%c",(*theData)(theOriginalOrder(j),theNodeMap(i),1));
                }
                fprintf (dest, "\n");
            } else {
                rec << seqDelimiter;
                rec << *curName;
                for ( j = 0; j<theOriginalOrder.lLength; j++) {
                    if (j%printWidth==0) {
                        rec << '\n';
                    }
                    rec << (*theData)(theOriginalOrder(j),theNodeMap(i),1);
                }
                rec << '\n';
            }
        }
    }

    // finally see if we need to write out a tree

    }

    if (outputFormat != 8) {
        _Parameter  treeDefined;
        checkParameter (dataFileTree, treeDefined,0.0);
        if (treeDefined) {
            _Variable *treeVar = FetchVar(LocateVarByName (dataFileTreeString));
            if (treeVar) {
                _String* treeString = (_String*)(treeVar->Compute())->toStr();
                switch (outputFormat) {
                case 0:
                case 1:
                case 9:
                case 10: {
                    if (dest) {
                        fprintf (dest,"\n\n%s;",treeString->sData);
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
                        fprintf (dest,"\n1\n%s;",treeString->sData);
                    } else {
                        rec << "\n1\n";
                        rec << *treeString;
                    }
                    break;
                }
                default: {
                    if (dest) {
                        fprintf (dest,"\n\nBEGIN TREES;\n\tTREE tree = %s;\nEND;",treeString->sData);
                    } else {
                        rec << "\n\nBEGIN TREES;\n\tTREE tree = ";
                        rec << *treeString;
                        rec << ";\nEND;";
                    }
                }

                }
                DeleteObject (treeString);
            }
        }
    }

}

//_________________________________________________________

bool    StoreADataSet (_DataSet* ds, _String* setName)
{
    if (!setName->IsValidIdentifier (true)) {
        WarnError (*setName & " is not a valid identifier while constructing a DataSet");
        return false;
    }

    long pos = FindDataSetName (*setName);

    if (pos==-1) {
        dataSetNamesList << setName;
        dataSetList<<ds;
        DeleteObject (ds);
    } else {
#if !defined __UNIX__ && ! defined __HEADLESS__
        if (!RequestDataSetReplace (pos)) {
            terminateExecution = true;
            DeleteObject (ds);
            return false;
        }
#endif

        _DataSet* existingDS = (_DataSet*)dataSetList (pos);

        bool isDifferent = existingDS->NoOfSpecies () != ds->NoOfSpecies() ||
                           existingDS->NoOfColumns () != ds->NoOfColumns() ||
                           existingDS->NoOfUniqueColumns () != ds->NoOfUniqueColumns() ||
                           existingDS->GetTT () != ds->GetTT();


        for (long dfIdx = 0; dfIdx < dataSetFilterNamesList.lLength; dfIdx++)
            if (((_String*)dataSetFilterNamesList(dfIdx))->sLength) {
                _DataSetFilter * aDF = (_DataSetFilter*)dataSetFilterList(dfIdx);
                if (aDF->GetData() == existingDS) {
                    if (isDifferent) {
                        ReportWarning (_String("Overwriting dataset '") & *setName & "' caused DataSetFilter '" & *((_String*)dataSetFilterNamesList(dfIdx)) & "' to be deleted");
                        KillDataFilterRecord(dfIdx, false);
                    } else {
                        aDF->SetData(ds);
                    }
                }
            }
        dataSetList.Replace(pos,ds,false);
    }


    CheckReceptacleAndStore (*setName&".species",empty,false, new _Constant (ds->NoOfSpecies()), false);
    CheckReceptacleAndStore (*setName&".sites",empty,false, new _Constant (ds->NoOfColumns()), false);
    CheckReceptacleAndStore (*setName&".unique_sites",empty,false, new _Constant (ds->NoOfUniqueColumns()), false);

    return true;
}
