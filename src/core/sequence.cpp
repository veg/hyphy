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

#include "sequence.h"
#include "errorfns.h"
#include "stdio.h"
#include "string.h"
#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

_String   NuclAlphabet          = "ACGT-?",
          CodonAlphabet       = "ABCDEFGHIJKLMNOPQRSTUVWXYZ*?-.",
          FullAlphabet,
          CompleteNuclAlphabet  = "AGCTUYRWSKMBDHVXN?0-.";

long     countCompress = 0,
         countDecompress = 0;

void     initFullAlphabet  (void);
void     WriteBitsToString (_String&, long&, char);

//_________________________________________________________

void    initFullAlphabet (void)
{
    _String fA ((unsigned long)256);
    for (long i=0; i<256; i++) {
        fA[i]=i;
    }
    FullAlphabet = fA;
}

//_________________________________________________________

_CString::_CString (void)
{
    allocatedSpace = 0;
    if (FullAlphabet.sLength==0) {
        initFullAlphabet();
    }
    compressionType = NOCOMPRESSION;
}

//_________________________________________________________

_CString::_CString (_String&s): _String(s)
{
    allocatedSpace = 0;
    if (FullAlphabet.sLength==0) {
        initFullAlphabet();
    }
    compressionType = NOCOMPRESSION;
}
//_________________________________________________________
_CString::_CString (char* s): _String(s)
{
    allocatedSpace = 0;
    if (FullAlphabet.sLength==0) {
        initFullAlphabet();
    }
    compressionType = NOCOMPRESSION;
}
//_________________________________________________________
_CString::_CString (char s): _String(s)
{
    allocatedSpace = 0;
    if (FullAlphabet.sLength==0) {
        initFullAlphabet();
    }
    compressionType = NOCOMPRESSION;
}

//_________________________________________________________
_CString::_CString (unsigned long sL, bool flag)
{
    if (flag) {
        sLength = 0;
        if (sL<storageIncrement) {
            sL = storageIncrement;
        }
        sData = (char*)MemAllocate (sL*sizeof (char));
        allocatedSpace = sL;
        if (!sData) {
            warnError( -108);
        }
    } else {
        allocatedSpace = 0;
        sLength = sL;
        sData = (char*)MemAllocate (sL+1);
        if (sData) {
            memset (sData,0,sL+1);
        } else {
            sLength = 0;
            isError(0);
        }
    }
    compressionType = NOCOMPRESSION;
}

//_________________________________________________________
_CString::~_CString(void)
{
}
//_________________________________________________________
long    _CString::FreeUpMemory(long)
{
    if (!IsCompressed()) {
        _Parameter comprratio = BestCompress (NUCLEOTIDEALPHABET);
        if (comprratio == 1) {
            comprratio = BestCompress (CODONALPHABET);
        }
        return sLength*(1/comprratio-1);
    }
    return 0;

}

//_______________________________________________________________________
// append operator
void _CString::Finalize (void)
{
    sData = MemReallocate (sData,sLength+1);

    if (!sData) {
        warnError(-108);
    }

    sData[sLength]=0;
    allocatedSpace = 0;
}

//_______________________________________________________________________
// append operator
void _CString::operator << (char c)
{
    if (allocatedSpace <= sLength) {
        unsigned long incBy = ((storageIncrement*8 > sLength)? storageIncrement: (sLength/8+1));

        allocatedSpace+=incBy;

        sData = (char*)MemReallocate((char*)sData, allocatedSpace*sizeof(char));

        if (!sData) {
            checkPointer (sData);
        }
    }

    sData[sLength++]=c;
}

//_______________________________________________________________________
// append operator
void _CString::operator << (_String* s)
{
    if ( s && s->sLength) {
        if (allocatedSpace < sLength + s->sLength) {
            unsigned long incBy = sLength + s->sLength - nInstances;

            if (incBy < storageIncrement) {
                incBy = storageIncrement;
            }

            if (incBy < sLength/8) {
                incBy = sLength/8;
            }

            allocatedSpace+=incBy;

            sData = (char*)MemReallocate((char*)sData, allocatedSpace*sizeof(char));

            if (!sData) {
                checkPointer (sData);
            }
        }

        memcpy(sData+sLength,s->sData,s->sLength);
        sLength+=s->sLength;
    }
}

//_________________________________________________________

BaseRef   _CString::makeDynamic (void)
{
    _CString* res = (_CString*)new _CString;
    checkPointer(res);

    _String::Duplicate (res);

    res->compressionType = compressionType;

    return  res;

}

//_________________________________________________________

void      _CString::Duplicate (BaseRef res)
{
    _String::Duplicate (res);

    ((_CString*)res)->compressionType = compressionType;
}

//_________________________________________________________

unsigned char powersOf2[9]= {0,2,6,14,30,62,126,254,0};
unsigned char realPowersOf2[8]= {1,2,4,8,16,32,64,128};

// auxiliary string bit write function

//_________________________________________________________

void    WriteBitsToString (_String&s, long& bitAt, char lengthToWrite)
{
    long leftOver = 8-bitAt%8, curPos = bitAt/8;
    if (leftOver >= lengthToWrite) { // will fit in current byte
        unsigned char value = (unsigned char)s[curPos];
        value += powersOf2[leftOver-1]-powersOf2[leftOver-lengthToWrite];
        s[curPos]=value;
    } else {
        unsigned char value = (unsigned char)s[curPos];
        value += powersOf2[leftOver-1]+1;
        s[curPos]=value;
        char fullBytes = (lengthToWrite-leftOver-1)/8;
        while (fullBytes) {
            s[++curPos]=255;
            fullBytes--;
        }
        s[++curPos]=254-powersOf2[8-(lengthToWrite-leftOver)%8];
    }
    bitAt+=lengthToWrite;
}

//_________________________________________________________

/*
void    PrintStringBitWise (char* str, long theL)
{
    unsigned int theByte;
    char    byteout[9];
    byteout[8] = 0;
    long    i = 0;
    for (;i<theL;i++)
    {
        if(i) printf(", ");
        theByte = (unsigned int) (str[i]);
        for (int k=0; k<8; k++)
        {
            byteout[7-k]='0'+theByte%2;
            theByte/=2;
        }
        printf ("%s", byteout);
    }
}
*/



//_________________________________________________________
_String* _CString::SelectAlpha (unsigned char alpha)
{
    alpha&=0xf0;
    switch (alpha) {
    case NUCLEOTIDEALPHABET:
        return &NuclAlphabet;
    case CODONALPHABET:
        return &CodonAlphabet;
    case FULLNUCLALPHABET:
        return &CompleteNuclAlphabet;
    }
    return &FullAlphabet;
}


//_________________________________________________________

_Parameter      _CString::FrequencyCompress(unsigned char theAlpha,bool doit)
{

    _String* theAlphabet = SelectAlpha (theAlpha);
    if ((*theAlphabet).sLength>31) {
        return 1;    // can't do much - the alphabet is too large
    }
    char codeLength[256];
    long freqs [256],j,t;
    long maxOccurences[256], locationsOfMaxSymbols[256] ; //simply ensures that we
    // won't have symbols out of the alphabet


    //analyze the frequency distribution of alphabetic symbols

    for (j=0; j<256; freqs[j]=0,codeLength[j]=0,maxOccurences[j]=0, j++ ) {}


    for (j=0; j<sLength; j++) {
        freqs[sData[j]]++;
    }

    t = 0;
    for (j=0; j<theAlphabet->sLength; j++) {
        freqs[NuclAlphabet[j]]*=-1;
    }

    //make sure that the alphabet is "large" enough for the nucleotide case
    // NEW 03/29/98
    for (j=0; j<256; j++)
        if (freqs[j]>0) {
            t = 1;
            break;
        } else {
            freqs[j]*=-1;
        }
    if (t) {
        if (theAlphabet == &NuclAlphabet) {
            return FrequencyCompress (FULLNUCLALPHABET, doit);
        } else {
            return 1;
        }
    }



    // now build the prefix code for the alphabet
    // fisrt find four most frequently occurring symbols

    for (j=0; j<(*theAlphabet).sLength; j++) {
        for (long k = 0; k<(*theAlphabet).sLength; k++)
            if (freqs[(*theAlphabet)[j]]>=maxOccurences[k]) {
                for (long l=(*theAlphabet).sLength-1; l>=k+1; l--) {
                    maxOccurences[l]=maxOccurences[l-1];
                    locationsOfMaxSymbols[l]=locationsOfMaxSymbols[l-1];
                }
                maxOccurences[k]=freqs[(*theAlphabet)[j]];
                locationsOfMaxSymbols[k]=(*theAlphabet)[j];
                break;
            }
    }


    // compute efficiency
    //j will store the predicted bit length of the compressed string

    j = (*theAlphabet).sLength*5; // translation table size
    j=8*((j%8)?(j/8+1):j/8);

    // we are also ready to build the code table

    for (long k = 0; k<(*theAlphabet).sLength; k++) {
        long l;
        for (l=0; l<(*theAlphabet).sLength; l++)
            if ((*theAlphabet)[k]==locationsOfMaxSymbols[l]) {
                j+=(l+1)*freqs[(*theAlphabet)[k]];
                codeLength [locationsOfMaxSymbols[l]] = l+1;
                break;
            }
    }

//  if (j>Length()*8) return 1;
// no compression could be performed
    if (!doit) {
        return j/8.0/sLength;
    }

    _String result ((unsigned long)(j%8?j/8+1:j/8)); // allocate output string


// let's roll!!
    long csize = 0; //will indicate the current bit position in the target string
    t = 0; // current position in the string
    //first we must write out the encoding table as 5 bits of length per each

    for (j=0; j<(*theAlphabet).sLength; j++, csize+=5, t = csize/8) {
        long leftover = 8-csize%8;
        if (leftover>=5) {
            unsigned char value = result[t];
            switch (leftover) {
            case 5:
                value+=codeLength[(*theAlphabet)[j]];
                break;
            case 6:
                value+=codeLength[(*theAlphabet)[j]]*2;
                break;
            case 7:
                value+=codeLength[(*theAlphabet)[j]]*4;
                break;
            default:
                value+=codeLength[(*theAlphabet)[j]]*8;
            }
            result[t]=value;
        } else {
            result[t]+=codeLength[(*theAlphabet)[j]]/realPowersOf2[5-leftover];
            result[++t]=(codeLength[(*theAlphabet)[j]]%realPowersOf2[5-leftover])*realPowersOf2[3+leftover];
        }
    }



    //  result[++t]=0;
    // mark the end of tabular encoding
    t++;
    // now encode the actual sequence
    t*=8;
    //t+=8;

    for (j=0; j<sLength; j++) {
        WriteBitsToString (result,t,codeLength[(unsigned char)sData[j]]);
    }

    // pad the rest of the last byte in the string by ones

    if (t%8) {
        unsigned char value = result [t/8];
        value += powersOf2[7-t%8]+1;
        result[t/8]=value;
        t++;
    }

    // yahoo! we are done - store compression flag and replace the string with compressed string
    _Parameter factor = result.sLength/(_Parameter)sLength;
    if (factor<1) { // compression took place
        DuplicateErasing(&result);
        SetFlag( FREQCOMPRESSION);
        SetFlag (theAlpha);
    }
    return factor;

}

//_________________________________________________________

_String*    _CString::DecompressFrequency(void)
{

    _String* theAlphabet = SelectAlpha (compressionType);

    if (!IsFlag(FREQCOMPRESSION)) {
        return nil;    // wrong compression type nothing to do
    }
    unsigned char *codeMaps = new unsigned char [(*theAlphabet).sLength];

    if (!codeMaps) {
        warnError( -108);    // no memory
    }

    unsigned int i,j,k,l,t; // temporary vars


    // read in the alphabet encoding
    i=0;
    t=0;

    for (j=0; j<(*theAlphabet).sLength; j++, i+=5, t = i/8) {
        long leftover = 8-i%8;
        unsigned char value = sData[t];
        if (leftover>=5) {
            switch (leftover) {
            case 5:
                value%=32;
                break;
            case 6:
                value=(value%64)/2;
                break;
            case 7:
                value=(value%128)/4;
                break;
            default:
                value=value/8;
            }

        } else {
            value=((unsigned char) sData[t])%realPowersOf2[leftover]*realPowersOf2[5-leftover];
            value+=((unsigned char) sData[t+1])/realPowersOf2[3+leftover];
        }
        codeMaps[value-1]=j;
    }

    if(i%8) {
        t++;    //possible rounding error
    }

    // now read in the data
    // first we must guess the correct size of the string, but to be safe we'll just set it to
    // the maximum possible value (i.e. if all the bytes compressed to 1 bit)

    _String     result (10,true);

    // k will count the actual length
    // j,l will be used for relative positions of 0's

    for (k=0,j=t*8;; ) { // go by bits
        l=j; // look for the next zero bit in the stream
        while (1) {
            for (i=8-l%8; i>0; i--) {
                if (!(sData[t]&realPowersOf2[i-1])) {
                    break;
                }
            }
            if (i) {
                l+=(8-i)-l%8;
                break;
            }
            if (t<sLength-1)    {
                t++;
            } else {
                l=0;
                break;
            }
            l=t*8;
        }
        if (!l) {
            break;
        }
        l++;
        _String addOn (theAlphabet->getChar(codeMaps[l-j-1]));
        result<<&addOn;
        if ((t=l/8)>=sLength) {
            break;
        }
        j=l;
        k++;
    }
    result.Finalize();
    delete [] codeMaps;
    return (_String*)(_String (result.getStr())).makeDynamic();

}

//_________________________________________________________
inline unsigned int     ToLZWCode (long l)
{
    return (l>127)?l|0x8000:l;
}
#include "stdio.h"

//_________________________________________________________

_Parameter      _CString::LZWCompress (unsigned char theAlpha)
{
    _List theTable;
    _SimpleList theCodes;

    _String* theAlphabet = SelectAlpha (theAlpha);

    _String output (*this), curString(""), testString;
    long k = 0, pos, codeMax = 0, pos1;
    char checkTable [256];

    for (; k<256; checkTable[k++]=0) {}
    for (k=0; k<theAlphabet->sLength; checkTable[(*theAlphabet)[k++]]=1) {}

    // init the table

    for (long j = 0; j<(*theAlphabet).sLength; j++) {
        _String a((*theAlphabet)[j]);
        theCodes.InsertElement ((BaseRef)ToLZWCode(j),theTable.BinaryInsert(&a),false,false);
    }

    codeMax = (*theAlphabet).sLength;

    for (long p=0; p<sLength; p++) {
        if (!checkTable[sData[p]]) {
            return 1;    // symbol not in alphabet - can't compress
        }
        testString = curString&sData[p];
        pos = theTable.BinaryFind(&testString);
        if (pos<0) {
            pos=theTable.BinaryInsert(&testString);
            theCodes.InsertElement ((BaseRef)ToLZWCode(codeMax++),pos,false,false);
            curString = sData[p];
            long acode = theCodes(pos1);
            pos1 = theTable.BinaryFind(&curString);
            if (acode>127) {
                output[k+1]=acode%256;
                output[k]=acode/256;
                k+=2;
            } else {
                output[k++]=acode;
            }

        } else {
            pos1=pos;
            curString = testString;
        }
    }

    long acode = theCodes(pos1);
    if (acode>127) {
        output[k+1]=acode%256;
        output[k]=acode/256;
        k+=2;
    } else {
        output[k++]=acode;
    }

    output[k]=0;

    // debugging info
    /*  _String* bs = (_String*)theTable.toStr(), *bss = (_String*)theCodes.toStr();
        printf ("\nString Table:%s\nCode Table:%s\n",bs->getStr(), bss->getStr());
        DeleteObject(bs);
        DeleteObject(bss);
        for (long j=0;j<k; j++)
        {
            unsigned int value = output[j];
            if (value>127)
            {
                value&0x7f;
                value*=256;
                value+=output[++j];
            }
            printf ("%u,",value);
        }*/
    // end debugging
    output.SetLength(k+1);
    _Parameter factor = k/_Parameter(sLength);
    if (factor<1) {
        DuplicateErasing(&output);
        SetFlag( LZWCOMPRESSION);
        SetFlag (theAlpha);
    }
    return factor;
}
//_________________________________________________________
inline unsigned int     GetNextCode (_String& s, long& p)
{
    if (s.sData[p]<0) {
        unsigned int val = s.sData[p++]&0x7f;
        val *= 256;
        val += (unsigned char)(s.sData[p++]);
        return val;
    }
    return s[p++];
}

//_________________________________________________________

_String*        _CString::DecompressLZW (void)
{
    _String* theAlphabet = SelectAlpha (this->compressionType);


    if (!sLength||(!IsFlag( LZWCOMPRESSION))) {
        return nil;    //nothing to do!
    }
    _List theTable;

    _String output(storageIncrement,true), testString;
    long codeMax = 0, oldCode;

    // init the table

    for (long j = 0; j<(*theAlphabet).sLength; j++) {
        _String a((*theAlphabet)[j]);
        theTable&&(&a);
    }

    long p = 0;

    oldCode = GetNextCode(*this,p);
//  output = output&*  (_String*)theTable(oldCode);
    output <<  (_String*)theTable(oldCode);

    for (; p<sLength-1;) {
        codeMax=GetNextCode(*this,p);
        if (theTable.countitems()-1>=codeMax) {
            output << (_String*)theTable(codeMax);
//          output = output& *(_String*)theTable(codeMax);
            _String addOn((*(_String*)theTable(oldCode)));
            addOn =addOn&*((_String*)theTable(codeMax))[0];
            theTable&& &addOn;
            oldCode = codeMax;
        } else {
            testString =  *(_String*)theTable(oldCode);
            testString = testString&testString.getChar(0);
            theTable&&(&testString);
            output << &testString;
//          output = output & testString;
            oldCode = codeMax;
        }
    }
    output.Finalize();

    // debugging info
//  _String* bs = (_String*)theTable.toStr();
//  printf ("\nString Table:%s\n",bs->getStr());
//  DeleteObject(bs);

    // end debugging

    return (_String*)output.makeDynamic();
}

//_________________________________________________________

_Parameter _CString::BestCompress(unsigned char theAlpha, long triggerSize)
{
    countCompress++;
    _Parameter freqcomp = FrequencyCompress(theAlpha, false), lzwcomp = 1;
    _CString test(*this);
    if ((triggerSize>=sLength)||(triggerSize==-1)) {
        lzwcomp  = test.LZWCompress (theAlpha);
    }
    if ((freqcomp<1)||(lzwcomp<1)) { // stuff to do
        if (freqcomp<lzwcomp) {
            FrequencyCompress(theAlpha, true);
            return freqcomp;
        } else {
            DuplicateErasing (&test);
            compressionType = test.compressionType;
            return lzwcomp;
        }
    }
    SetDecompressed ();
    return 1;
}

//_________________________________________________________

_String* _CString::Decompress(void)
{
    countDecompress++;
    if (IsFlag(LZWCOMPRESSION)) {
        return DecompressLZW();
    }
    if (IsFlag(FREQCOMPRESSION)) {
        return DecompressFrequency();
    }
    if (compressionType==NOCOMPRESSION) {
        return (_String*)toStr();
    }
    return nil;
}
