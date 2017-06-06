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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>

#ifdef   __UNIX__
#if !defined __MINGW32__
#include <sys/utsname.h>
#endif
#include <unistd.h>
#endif

#include "hy_strings.h"
#include "gnuregex.h"
#include "batchlan.h"
#include "global_things.h"
#include "mersenne_twister.h"
#include "function_templates.h"

using namespace hy_global;


#define MOD_ADLER 65521

_String   compileDate = __DATE__,
          __HYPHY__VERSION__ = _String ("2.4") & compileDate.Cut (7,10) & compileDate.Cut (0,2).Replace("Jan", "01", true).
                                                                                                  Replace("Feb", "02", true).
                                                                                                  Replace("Mar", "03", true).
                                                                                                  Replace("Apr", "04", true).
                                                                                                  Replace("May", "05", true).
                                                                                                  Replace("Jun", "06", true).
                                                                                                  Replace("Jul", "07", true).
                                                                                                  Replace("Aug", "08", true).
                                                                                                  Replace("Sep", "09", true).
                                                                                                  Replace("Oct", "10", true).
                                                                                                  Replace("Nov", "11", true).
                                                                                                  Replace("Dec", "12", true)
                                                                                                  & compileDate.Cut (4,5).Replace (" ", "0", true) & "alpha";

 
_String     emptyAssociativeList ("{}"),
            hyphyCiteString ("\nPlease cite S.L. Kosakovsky Pond, S. D. W. Frost and S.V. Muse. (2005) HyPhy: hypothesis testing using phylogenies. Bioinformatics 21: 676-679 if you use HyPhy in a publication\nIf you are a new HyPhy user, the tutorial located at http://www.hyphy.org/docs/HyphyDocs.pdf may be a good starting point.\n");

// TODO SLKP 20170517 Possibly deprecate these global strings


const struct _hyValidIDCharsType {
    bool valid_chars[256];
    _hyValidIDCharsType(void) {
        
        InitializeArray (valid_chars, 256, false);
        
        for (unsigned char c='a'; c<='z'; c++) {
                valid_chars[c] = true;
        }
        for (unsigned char c='A'; c<='Z'; c++) {
            valid_chars[c] = true;
        }
        for (unsigned char c='0'; c<='9'; c++) {
                valid_chars[c] = true;
        }
        valid_chars[(unsigned char)'_'] = true;
    }
}
_hyValidIDChars;


/*
==============================================================
Constructors/Destructors/Copiers
==============================================================
*/

_String::_String (void) {
  Initialize();
}

//=============================================================

void _String::Initialize (bool) {
    s_length = 0UL;
    s_data = nil;
}

//=============================================================

_String::_String(long const number) {
    char s[64];
    s_length = snprintf(s, sizeof(s), "%ld", number);
    AllocateAndCopyString (s, s_length);
}

//=============================================================

_String::_String(const unsigned long sL) {
    s_length = sL;
    s_data = (char*) MemAllocate(sL + 1L, true);
}

//=============================================================

_String::_String(const hy_float val, const char *format) {
    char s_val[128];
    s_length = snprintf(s_val, 128, format ? format : PRINTF_FORMAT_STRING, val);
    AllocateAndCopyString (s_val, s_length);
}

//=============================================================

_String::_String(const _String &s) {
    Initialize ();
    Duplicate(& s);
}

//=============================================================

_String::_String(_String *s) {
    if (s->CanFreeMe ()) {
        s_data       = s->s_data;
        s_length     = s->s_length;
        s->s_data    = nil;
        DeleteObject (s);
    } else {
        AllocateAndCopyString (s_data, s_length);
        s->RemoveAReference();
    }
}


//=============================================================
_String::_String(const _String &source, long start, long end) {
    if (source.s_length) {
        
        long requested_range = source.NormalizeRange(start, end);
        
        if (requested_range > 0L) {
            AllocateAndCopyString (source.s_data + start, requested_range);
            return;
            
        }
    }
    
    s_length = 0UL;
    s_data = (char *)MemAllocate(1UL);
    s_data[0] = '\0';
    
}

//=============================================================

_String::_String(const char *c_string) {
    AllocateAndCopyString (c_string, strlen(c_string));
}

//=============================================================
_String::_String(const wchar_t *wc_string) {
    unsigned long allocated = wcslen (wc_string);
    s_length = 0UL;
    s_data = (char *)MemAllocate(allocated + 1UL);
    for (unsigned long cid = 0UL; cid < allocated; cid ++) {
        int this_char = wctob (wc_string[cid]);
        if (this_char != WEOF) {
            s_data[s_length++] = (char) this_char;
        }
    }
    if (s_length != allocated) {
        s_data = (char *)MemReallocate((char *)s_data, (s_length+1) * sizeof(char));
    }
    s_data[s_length] = '\0';
}

//=============================================================
_String::_String(const char c) {
    s_length = 1UL;
    s_data = (char *)MemAllocate(2UL);
    s_data[0] = c;
    s_data[1] = '\0';
}

//=============================================================
_String::_String (const _String& str, unsigned long copies) {
    s_length = copies * str.s_length;
    s_data = (char*)MemAllocate (s_length+1UL);
    if (s_length > 0UL) {
        for (unsigned long i = 0UL; i < copies; i++) {
            memcpy (s_data + i * str.s_length, s_data, str.s_length);
        }
    }
    s_data[s_length]='\0';
}

//=============================================================
_String::_String(FILE * file) {
    Initialize ();
    if (file) {
        fseek(file, 0, SEEK_END);
        s_length = (unsigned long) ftell(file);
        s_data = (char *)MemAllocate(s_length + 1UL);
        rewind(file);
        fread(s_data, 1, s_length, file);
        s_data[s_length] = '\0';
    }
}

//=============================================================
_String::~_String(void)
{
    if (CanFreeMe()) {
        if (s_data) {
            free (s_data);
            s_data = nil;
        }
        s_length = 0UL;
    } else {
        RemoveAReference();
    }
}

//=============================================================
BaseRef _String::makeDynamic (void) const {
    _String * r = new _String;
    r->Duplicate(this);
    return r;
}

//=============================================================
void    _String::Duplicate (BaseRefConst ref) {
    if (s_data) {
        free (s_data);
    }
    
    _String const * s = (_String const*)ref;
    
    s_length = s->s_length;
    s_data   = s->s_data;
    
    if (s_data) {
        AllocateAndCopyString (s->s_data, s_length);
        s_data = (char*)MemAllocate (s_length+1L);
        memcpy (s_data, s->s_data, s_length+1);
    }
}

//=============================================================
void _String::operator = (_String const& s) {
    Duplicate (&s);
}



/*
 ==============================================================
 Private helpers
 ==============================================================
 */

long  _String::NormalizeRange(long & from, long & to) const {
    
    if (s_length == 0UL) {
        return 0L;
    }
    
    if (from < 0L) {
        from = 0L;
    }
    
    if (to < 0L || to >= s_length ) {
        to = s_length - 1UL;
    }
    
    return to - from + 1L;
    
}

//=============================================================

void _String::AllocateAndCopyString (const char * source_string, unsigned long length) {
    s_length = length;
    s_data = (char*) MemAllocate (length+1UL);
    if (s_length) {
        memcpy (s_data, source_string, length);
    }
    s_data [length] = '\0';
}


/*
==============================================================
Getters and setters
==============================================================
*/

char& _String::operator [] (long index) {
    if (index < s_length && index >= 0L) {
        return s_data[index];
    }
    HandleApplicationError (_String ("Internal error at ") & __PRETTY_FUNCTION__ & ": an invalid index requested");
}

//=============================================================


char _String::operator () (long index) {
    return get_char(index);
}

//=============================================================

const char _String::get_char(long index) const {
    if (index >= 0L && index < s_length) {
        return s_data[index];
    }
    return default_return;
}

//=============================================================

void _String::set_char (unsigned long index, char const data) {
    if (index < s_length) {
        s_data[index] = data;
    }
}

/*
 ==============================================================
 Comparisons
 ==============================================================
 */

_HY_COMPARISON_TYPE _String::Compare(_String const& rhs) const {
    unsigned long up_to = MIN (s_length, rhs.s_length);
    
    for (unsigned long i = 0UL; i < up_to; i++) {
        if (s_data[i] < rhs.s_data[i]) {
            return kCompareLess;
        }
        if (s_data[i] > rhs.s_data[i]) {
            return kCompareGreater;
        }
    }

    if (s_length == rhs.s_length) {
        return kCompareEqual;
    }
    
    return s_length < rhs.s_length ? kCompareLess : kCompareGreater;
}

_HY_COMPARISON_TYPE _String::CompareIgnoringCase(_String const& rhs) const {
    unsigned long up_to = MIN (s_length, rhs.s_length);
    
    for (unsigned long i = 0UL; i < up_to; i++) {
       
        char llhs = tolower (s_data[i]), lrhs = tolower (rhs.s_data[i]);
        
        
        if (llhs < lrhs) {
            return kCompareLess;
        }
        if (llhs > lrhs) {
            return kCompareGreater;
        }
    }
    
    if (s_length == rhs.s_length) {
        return kCompareEqual;
    }
    
    return s_length < rhs.s_length ? kCompareLess : kCompareGreater;
}

bool _String::operator==(const _String& s) const { return Compare  (s) == kCompareEqual; }
bool _String::operator>(const _String & s) const { return Compare  (s) == kCompareGreater; }
bool _String::operator<=(const _String & s) const { return Compare (s) != kCompareGreater; }
bool _String::operator>=(const _String & s) const { return Compare (s) != kCompareLess; }
bool _String::operator!=(const _String & s) const { return Compare (s) != kCompareEqual; }
bool _String::operator<(const _String & s) const { return Compare  (s) == kCompareLess; }

bool _String::Equal(const _String& s) const { return Compare  (s) == kCompareEqual; };
bool _String::EqualIgnoringCase(const _String& s) const { return CompareIgnoringCase  (s) == kCompareEqual; };
bool _String::Equal(const char c) const {
    return s_length == 1UL && s_data[0] == c;
}


bool _String::EqualWithWildChar(const _String& pattern, const char wildchar, unsigned long start_this, unsigned long start_pattern) const {
    // wildcards only matter in the second string
    
    if (pattern.s_length > start_pattern && wildchar != '\0') {
        unsigned long   match_this_char = start_pattern;
        // the position we are currently trying to match in the pattern
        
        bool            is_wildcard = pattern.s_data[match_this_char] == wildchar,
        scanning_pattern = is_wildcard;
        
        unsigned long i = start_this;
        // the position we are currently trying to match in *this
        
        while (i <= s_length) {
            
            if (scanning_pattern) { // skip consecutive wildcards in "pattern"
                scanning_pattern = pattern.s_data[++match_this_char] == wildchar;
            } else {
                if (s_data[i] == pattern.s_data[match_this_char]) {
                    
                    if (is_wildcard) {
                        // could either match the next char or consume it into the wildcard
                        if (EqualWithWildChar (pattern, wildchar, i, match_this_char)) {
                            // matching worked
                            return true;
                        } else {
                            i++;
                            continue;
                        }
                    } else {
                        // try character match
                        // note that the terminal '0' characters will always match, so
                        // this is where we terminate
                        i++;
                        match_this_char++;
                        if (i > s_length || match_this_char > pattern.s_length) {
                            break;
                        }
                        is_wildcard =  pattern.s_data[match_this_char] == wildchar;
                        scanning_pattern = is_wildcard;
                    }
                } else { // match wildcard
                    if (!is_wildcard) {
                        return false;
                    }
                    scanning_pattern = false;
                    i++;
                }
            }
        }
        
        return match_this_char > pattern.s_length;
    } else {
        return s_length == start_this;
    }
    
    return false;
}


/*
 ==============================================================
 Content-modification and extraction methods
 ==============================================================
 */




//Append operator
const _String _String::operator & (const _String& rhs) const {
    unsigned long combined_length = s_length + rhs.s_length;
    
    if (combined_length == 0UL) {
        return kEmptyString;
    }
    
    _String res(combined_length);
    
    if (s_length) {
        memcpy(res.s_data, s_data, s_length);
    }
    
    if (rhs.s_length) {
        memcpy(res.s_data + s_length, rhs.s_data, rhs.s_length);
    }
    
    res.s_data[res.s_length] = '\0';
    return res;
}

const _String _String::Chop(long start, long end) const{
    
    long resulting_length = NormalizeRange(start,end);
    
    if (resulting_length > 0L) {
        _String res((unsigned long)(s_length - resulting_length));
        if (start > 0L) {
            memcpy(res.s_data, s_data, start);
        }
        if (end + 1L < s_length) {
            memcpy(res.s_data + start, s_data + end + 1L, s_length - end - 1L);
        }
        
        return res;
    }
    
    return *this;
    
}

const _String _String::Cut(long start, long end) const {
    return _String (*this, start, end);
}

void _String::Delete(long start, long end) {
    long resulting_length = NormalizeRange(start,end);
    
    if (resulting_length > 0L) {
        if (end < (long)s_length - 1UL) {
            memmove(s_data + start, s_data + end + 1L, s_length - end - 1L);
        }
        s_length -= resulting_length;
        s_data = (char*)MemReallocate(s_data, sizeof(char) * (s_length + 1UL));
        s_data[s_length] = '\0';
    }
}

void _String::Flip(void) {
    for (unsigned long i = 0UL; i < (s_length >> 1); i++) {
        char c;
        SWAP  (s_data[i], s_data[s_length - 1 - i], c);
    }
}


const _String _String::Reverse(void) const {
    
    _String result (*this);
    for (unsigned long s = 0UL, e = s_length - 1L;  s < s_length; s++, e--) {
        result.s_data[s] = s_data[e];
    }
    return result;
}

void _String::Insert(char c, long where) {
    if (where < 0L || where >= s_length) {
        where = s_length;
    }
    
    s_data = (char*)MemReallocate(s_data, sizeof(char) * (s_length + 2UL));
    
    if (where < s_length) {
        memmove(s_data + where + 1UL, s_data + where, s_length - where);
    }
    
    s_data[where] = c;
    s_data[++s_length] = '\0';
}

void _String::Trim(long start, long end) {
    
    long resulting_length = NormalizeRange(start, end);
    
    if (resulting_length > 0L) {
        if (start > 0L) {
            memmove(s_data, s_data + start, resulting_length);
        }
        s_length = resulting_length;
        s_data = (char*)MemReallocate(s_data, resulting_length + 1UL);
        s_data[resulting_length] = '\0';
    } else {
        s_length = 0UL;
        s_data = (char*)MemReallocate(s_data, 1UL);
        s_data[0] = '\0';
    }
}

/*
==============================================================
Methods
==============================================================
*/

//Return good ole char*
_String::operator const char* (void) const {
    return sData;
}

// Compute Adler-32 CRC for a string
// Implementation shamelessly lifted from http://en.wikipedia.org/wiki/Adler-32
long _String::Adler32(void)
{
    unsigned char  *data = (unsigned char*)sData;

    unsigned long len = sLength,
                  a   = 1,
                  b      = 0;

    while (len) {
        unsigned tlen = len > 5550 ? 5550 : len;
        len -= tlen;
        do {
            a += *data++;
            b += a;
        } while (--tlen);
        a = (a & 0xffff) + (a >> 16) * (65536-MOD_ADLER);
        b = (b & 0xffff) + (b >> 16) * (65536-MOD_ADLER);
    }

    if (a >= MOD_ADLER) {
        a -= MOD_ADLER;
    }

    b = (b & 0xffff) + (b >> 16) * (65536-MOD_ADLER);

    if (b >= MOD_ADLER) {
        b -= MOD_ADLER;
    }

    return b << 16 | a;
}


bool _String::ContainsSubstring(_String& s)
// -1, indicates that search term has not been found
{
    if (!sLength) {
        return false;
    }
    if (sLength<s.sLength) {
        return false;
    }
    char *sP = sData, *ssP = s.sData;
    for (long i=0; i<sLength-s.sLength; i++,sP++) {
        long j = 0;
        for (; (sP[j]==ssP[j])&&(j<s.sLength); j++) ;
        if (j==s.sLength) {
            return true;
        }
    }
    return false;
}

//Cut string from, to (-1 for any means from beginning/to end)
const _String _String::Cut(long from, long to) const
{
    if (sLength) {
      if (from == -1) {
          from = 0;
      }
      if (to == -1 || to >= sLength) {
          to = sLength-1UL;
      }
    
      if (to>=from) {
        _String res (to-from+1UL, false);
      
        for (unsigned long index = from; index <= to; index++) {
          res.sData[index-from] = sData[index];
        }
      
        return res;
      }
    }
    return kEmptyString;
}

//Delete range char operator
void _String::Delete (long from, long to)
{
    if (from<0) {
        from = 0;
    }

    if (to<0) {
        to = sLength-1;
    }

    if (to<sLength-1) {
        memmove (sData+from, sData+to+1, sLength-to-1);
    }
    sLength -= to-from+1;
    sData = (char*)MemReallocate (sData,sizeof(char)*(sLength+1));
    sData[sLength]=0;
}




//Replace string 1 with string 2, all occurences true/false
void _String::FormatTimeString(long time_diff)
{
    long secs = time_diff,
         mins = secs/60,
         hrs  = mins/60;

    mins = mins%60;
    secs = secs%60;
    if (hrs<10) {
        (*this) = _String('0')&hrs;
    } else {
        (*this) = _String(hrs);
    }
    (*this) = (*this) &':';
    if (mins<10) {
        (*this) = (*this)&_String('0')&mins;
    } else {
        (*this) = (*this)&_String(mins);
    }
    (*this) = (*this) &':';
    if (secs<10) {
        (*this) = (*this)&_String('0')&secs;
    } else {
        (*this) = (*this)&_String(secs);
    }
}

//Finalize buffer string
void _String::Finalize (void)
{
    if (!(sData = (char*)MemReallocate (sData, sLength+1))) {
      return;
    }

    sData[sLength]  = 0;
    buffer_allocation      = 0UL;

}

long _String::FindEndOfIdent(long start, long end, char wild)
{
    if(sLength==0) {
        return -1;
    }

    if (start == -1) {
        start = ((long)sLength)-1;
    }
    if (end == -1) {
        end = ((long)sLength)-1;
    }

    long i = start;

    for (; i<=end; i++)
        if (!(isalnum(sData[i])||sData[i]=='.'||sData[i]==wild||sData[i]=='_')) {
            break;
        }

    if (i>start+2 && sData[i-1] == '_' && sData[i-2] == '_') {
        return i-3;
    }

    return i-1;
}

// find first occurence of the string between from and to
long _String::Find(const _String s, long from, long to) const
// -1, indicates that search term has not been found
{
    if (!sLength) {
        return -1;
    }
    if (from == -1) {
        from = 0;
    }
    if (to == -1) {
        to = ((long)sLength)-1;
    }
    if (to<from) {
        return -1;
    }
    if (to-from+1<s.sLength) {
        return -1;
    }
    char *sP = sData+from, *ssP = s.sData;
    for (long i=from; i<=to-s.sLength+1; i++,sP++) {
        long j;
        for (j = 0; (sP[j]==ssP[j])&&(j<s.sLength); j++) ;
        if (j==s.sLength) {
            return i;
        }
    }
    return -1;
}

long _String::FindKMP(_String s, long from, long to)
// -1, indicates that search term has not been found
{
    //Reproduced from http://en.wikipedia.org/wiki/Knuth–Morris–Pratt_algorithm

    if (!sLength) {
        return -1;
    }
    if (from == -1) {
        from = 0;
    }
    if (to == -1) {
        to = ((long)sLength)-1;
    }
    if (to<from) {
        return -1;
    }
    if (to-from+1<s.sLength) {
        return -1;
    }

    char *sP = sData+from; //Start of Haystack substring
    char *ssP = s.sData;  //Start of Needle substring
    int m = 0; //beginning of the current match in haystack
    int i = 0; //the position of the current character in needle

    while(m+i < (to-m+i+1)) {
        if(ssP[i] == sP[m+i]) {
            if (i == (s.sLength-1)) {
                return m;
            }
            ++i;
        }

        else {
            m = m + i - this->kmpTable[i];
            if(this->kmpTable[i] > -1) {
                i = this->kmpTable[i];
            } else {
                i = 0;
            }
        }
    }

    return -1;
}

// Construct a KMP table
void _String::buildKmpTable(_String s)
{
    //Reproduced from http://en.wikipedia.org/wiki/Knuth–Morris–Pratt_algorithm
    int pos = 2;
    int cnd = 0;

    this->kmpTable = new int[sizeof(int) * sLength];

    this->kmpTable[0] = -1;
    this->kmpTable[1] =  0;

    char *ssP = s.sData;  //Start of Needle substring

    while (pos < s.sLength) {
        if(ssP[pos-1] == ssP[cnd]) {
            ++cnd;
            this->kmpTable[pos] = cnd;
            ++pos;
        }

        else if(cnd > 0) {
            cnd = this->kmpTable[cnd];
        }

        else {
            this->kmpTable[pos] = 0;
            ++pos;
        }
    }
}

// find first occurence of the string between from and to
// case insensitive
long _String::FindAnyCase (_String s, long from, long to)
// -1, indicates that search term has not been found
{
    if (!sLength) {
        return -1;
    }
    if (from == -1) {
        from = 0;
    }
    if (to == -1) {
        to = ((long)sLength)-1;
    }
    if (to<from) {
        return -1;
    }
    if (to-from+1<s.sLength) {
        return -1;
    }

    s.UpCase();
    char *sP = sData+from, *ssP = s.sData;
    for (long i=from; i<=to-s.sLength+1; i++,sP++) {
        long j;
        for (j = 0; (toupper(sP[j])==ssP[j])&&(j<s.sLength); j++) ;
        if (j==s.sLength) {
            return i;
        }
    }
    return -1;
}

long _String::ExtractEnclosedExpression (long& from, char open, char close, bool respectQuote, bool respectEscape)
{
    long   currentPosition = from,
           currentLevel    = 0;

    char       isQuote = 0;
    bool       doEscape = false;

    while (currentPosition < sLength) {
        char thisChar = sData[currentPosition];

        if (!doEscape) {
            if (thisChar == '"' && respectQuote && !doEscape) {
                isQuote = !isQuote;
            } else if (thisChar == open && !isQuote) {
                // handle the case when close and open are the same
                if (currentLevel == 1 && open == close && from < currentPosition) {
                    return currentPosition;
                }
                currentLevel++;
                if (currentLevel == 1) {
                    from = currentPosition;
                }
            } else if (thisChar == close && !isQuote) {
                currentLevel--;
                if (currentLevel == 0 && from < currentPosition) {
                    return currentPosition;
                }
                if (currentLevel < 0) {
                    return -1;
                }
            } else if (thisChar == '\\' && respectEscape && isQuote && !doEscape) {
                doEscape = true;
            }
        }

        else {
            doEscape = false;
        }

        currentPosition++;
    }

    return -1;
}


//Find first occurence of the string between from and to
long _String::Find(char s, long from, long to) const
{
    if (!sLength) {
        return -1;
    }
    if (from == -1) {
        from = 0;
    }
    if (to == -1) {
        to = ((long)sLength)-1;
    }
    if (to<from) {
        return -1;
    }
    //if (to-from<0) return -1;

    for (long i=from; i<=to; i++)
        if (sData[i]==s) {
            return i;
        }

    return -1;
}

//Find first occurence of the string between from and to
long _String::FindBackwards(_String const & s, long from, long to) const {
// -1, indicates that search term has not been found
    if (!sLength) {
        return -1;
    }
    if (from == -1) {
        from = 0;
    }
    if (to == -1) {
        to = ((long)sLength)-1;
    }
    if (to<from) {
        return -1;
    }
    if (to-from+1<s.sLength) {
        return -1;
    }
    char *sP = sData, *ssP = (s.sData);
    for (long i=to-s.sLength+1; i>=(long)from; i--) {
        long j;
        for (j = 0; (sP[i+j]==ssP[j])&&(j<s.sLength); j++) ;
        if (j==s.sLength) {
            return i;
        }
    }
    return -1;
}

//Find first occurence of the string
long _String::FindBinary(char s)
// -1, indicates that search term has not been found
{
    for (long i=0; i < sLength; ++i) {
        if (sData[i] == s) {
            return i;
        }
    }
    return -1;
}

long _String::FindTerminator (long from, _String const& terminators) const {
    long   currentPosition  = from,
           currentCurly     = 0L,
           currentSquare    = 0L,
           currentParen     = 0L;

    bool   isQuote = false,
           doEscape = false;

    while (currentPosition < sLength) {
        char thisChar = sData[currentPosition];
        if (!doEscape) {
            if (thisChar == '"' && !doEscape) {
                isQuote = !isQuote;
            } else {
                if (!isQuote) {
                    if (thisChar == '{') {
                        currentCurly ++;
                    } else if (thisChar == '[') {
                        currentSquare ++;
                    } else if (thisChar == '(') {
                        currentParen ++;
                    }
                    if (currentCurly > 0 && thisChar == '}') {
                        currentCurly --;
                    } else if (currentSquare > 0 && thisChar == ']') {
                        currentSquare --;
                    } else if (currentParen > 0 && thisChar == ')') {
                        currentParen --;
                    } else if (currentParen == 0 && currentSquare == 0 && currentCurly == 0)
                        for (long s = 0; s < terminators.sLength; s++)
                            if (thisChar == terminators.sData[s]) {
                                return currentPosition;
                            }
                } else {
                    if (thisChar == '\\' && isQuote && !doEscape) {
                        doEscape = true;
                    }
                }
            }
        } else {
            doEscape = false;
        }

        currentPosition++;
    }

    return -1;
}



//s[0]...s[sLength-1] => s[sLength-1]...s[0]
const _String _String::Flip(void) {
  
  for (long s = 0L, e = (long)sLength-1L;  s < e; s++, e--) {
    char c = sData[s];
    sData[s] = sData[e];
    sData[e] = c;
  }
  return *this;
}



// Return good ole char*
const char * _String::getStr (void) const
{
    return sData;
}



//Insert char operator
void _String::Insert (char c, long pos)
{
    if (pos<0) {
        pos = sLength;
    }

    sData = (char*)MemReallocate (sData,sizeof(char)*(sLength+2));

    if (pos<sLength) {
        memmove(sData+pos+1,sData+pos, sLength-pos);
    }
    sData[pos] = c;
    sLength++;
    sData[sLength] = 0;
}

long _String::LempelZivProductionHistory (_SimpleList* rec)
{
    if (rec) {
        rec->Clear();
    }

    if (sLength == 0) {
        return 0;
    }

    if (rec) {
        (*rec) << 0;
    }

    long   cp = 1,
           pH = 1;

    while (cp<sLength) {
        long maxExtension = 0;

        for (long ip = 0; ip < cp; ip++) {
            long sp = ip,
                 mp = cp;

            while ((mp<sLength) && (sData[mp] == sData[sp])) {
                mp++;
                sp++;
            }

            if (mp==sLength) {
                maxExtension = sLength-cp;
                break;
            } else {
                if ((mp = mp-cp+1)>maxExtension) {
                    maxExtension = mp;
                }
            }
        }

        cp = cp+maxExtension;
        if (rec) {
            (*rec) << cp-1;
        } else {
            pH ++;
        }
    }

    if (rec) {
        return rec->lLength;
    }

    return pH;
}



//Replace string 1 with string 2, all occurences true/false
const _String _String::Replace(const _String s, const _String d, bool flag) const
{
    if (!sLength) {
        return kEmptyString;
    }
    if (sLength<s.sLength) {
        return *this;
    }
    if (s.sLength == 0) {
        return (*this);
    }

    if (flag) { // replace all
        // max possible number of replaces
        unsigned long t = sLength, cp=0UL;

        // allocate space for positions of substring s in this
      
        long *finds = new long [t],
              curSlot = 0L;


        // find all substrings s in this
        finds[0]=Find(s);
        if (finds[0]!=-1) {
            curSlot++;
            while ((finds[curSlot]=Find(s,finds[curSlot-1]+s.sLength,-1))!=-1) {
                curSlot++;
            }
        }

        // calculate the length of resulting string

        _String Res(sLength-(s.sLength-d.sLength)*curSlot, false);

        if (!curSlot) { // not found
            delete [] finds;
            return *this;
        }

        char *rP = (Res.sData), *dsP =(d.sData), *sP=(sData);

        if (finds[0]) {
            memcpy(rP,sP,finds[0]);    //head of the string;
        }
        cp+=finds[0];

        for (t=0; t<curSlot-1; t++) { // do the replacing
            if (d.sLength) {
                memcpy(rP+cp,dsP,d.sLength);
            }
            cp+=d.sLength;
            if (finds[t+1]-finds[t]-s.sLength) {
                memcpy(rP+cp,sP+finds[t]+s.sLength,finds[t+1]-finds[t]-s.sLength);
            }
            cp+=finds[t+1]-finds[t]-s.sLength;
        }
        if (d.sLength) {
            memcpy(rP+cp,dsP,d.sLength);
        }
        cp+=d.sLength;
        if(sLength-finds[curSlot-1]-s.sLength) {
            memcpy(rP+cp,sP+finds[curSlot-1]+s.sLength,sLength-finds[curSlot-1]-s.sLength);
        }
        //tail
        delete [] finds;
        return Res;
    }

    //first occurrence replace
    long t = Find(s),cp=0;
    if (t==-1) {
        return *this;
    }
    // substring not found

    _String Res(sLength-(s.sLength-d.sLength));

    char *rP = Res.sData, *dsP =d.sData, *sP=sData;
    if(t) {
        memcpy(rP,sP,t);    //head of the string;
    }
    cp+=t;
    if (d.sLength) {
        memcpy(rP+cp,dsP,d.sLength);
    }
    cp+=d.sLength;
    if (sLength-t-s.sLength) {
        memcpy(rP+cp,sP+t+s.sLength,sLength-t-s.sLength);
    }
    //tail
    return Res;

}


_String* _String::Sort (_SimpleList* index)
{
    if (index) {
        index->Clear();
    }

    if (sLength) {
        _SimpleList charList (sLength);
        if (index) {
            for (unsigned long i=0; i<sLength; i++) {
                charList << sData[i];
                (*index) << i;
            }
            SortLists (&charList, index);
        } else {
            for (unsigned long i=0; i<sLength; i++) {
                charList << sData[i];
            }

            charList.Sort();
        }
        _String * sorted = new _String (sLength);
        for (unsigned long i=0; i<sLength; i++) {
            sorted->sData[i] = charList.lData[i];
        }

        return sorted;
    }

    return new _String;
}

void    _String::StripQuotes (void) {
    if (sLength&&(sData[sLength-1]=='"')&&(sData[0]=='"')) {
        Trim(1,sLength-2);
    }    
}

const _String _String::Enquote (char quote_char) const {
  _String quoted (2UL + sLength, true);
  quoted << quote_char  << *this << quote_char;
  quoted.Finalize();
  return quoted;
}

const _List _String::Tokenize (_String const& s) const {
    _List pieces;
  
     if (s.sLength > 0L) {
        long cp=0,cpp;
        while ((cpp = Find(s,cp,-1))!=-1) {
            if (cpp>cp) {
                pieces.AppendNewInstance (new _String (*this,cp,cpp-1));
            } else {
                pieces.AppendNewInstance(new _String());
            }

            cp=cpp+s.sLength;
        }

        pieces.AppendNewInstance (new _String (*this,cp,-1));
    }
  
    return pieces;
}

hy_float _String::toNum (void) const {
    if (sLength == 0UL) {
        return 0.;
    }
    char * endP;
    return strtod(sData,&endP);
}

long _String::toLong (void) const {
  if (sLength == 0) {
    return 0;
  }
  char * endP;
  return strtol(sData,&endP,10);
}

//Return good ole char*
BaseRef _String::toStr (unsigned long) {
    AddAReference();
    return this;
}

hy_float  _String::ProcessTreeBranchLength (void)
{
    hy_float res = -1.; 

    if (sLength) {
        if (sData[0]==':') {
            res = Cut(1,-1).toNum();
        } else {
            res = toNum();
        }    


        if (res < 1e-10) {
            res = 1e-10;
        }    
    }    

    return res; 
}


bool    _String::IsALiteralArgument (bool stripQuotes) {
    if (sLength >= 2) { 
        long from = 0, 
             to = ExtractEnclosedExpression (from,'"','"',false,true);

        if (from == 0 && to == sLength - 1) { 
            if (stripQuotes){
                Trim (1, sLength-2);
            }    
            return true;
        }    
    }    
    return false;
}

//TODO: This is a global function.
bool    hyIDValidator (_String* s)
{
    return s->IsValidIdentifier(false);
}

/*
==============================================================
Space Methods
==============================================================
*/

//Replace all space runs with a single space
void _String::CompressSpaces (void)
{
    _String temp (sLength+1,true);
    bool    skipping = false;

    for (long k=0; k<sLength; k++)
        if (isspace (sData[k])) {
            if (!skipping) {
                skipping = true;
                temp << ' ';
            }
        } else {
            temp << sData[k];
            skipping = false;
        }
    temp.Finalize();
    *this = temp;
}

//Locate the first non-space charachter of the string
char _String::FirstNonSpace(long start, long end, char direction)
{
    long r = FirstNonSpaceIndex(start,end,direction);
    return r==-1?0:sData[r];
}

//Locate the first non-space charachter of the string
long _String::FirstNonSpaceIndex(long start, long end, char direction) const {
    if (start == -1) {
        start = ((long)sLength)-1;
    }
    if (end == -1) {
        end = ((long)sLength)-1;
    }
    if (direction<0) {
        //long t = start;
        start = end;
        end = start;
    }
    if (sLength&&(start<sLength)&&(!isspace (sData[start]))) {
        return start;    // first char is non-space
    }
 
    for (int i = start; i<=end; i+=direction) {
      if (!isspace (sData[i])) {
        return i;
      }
    }

    return -1;
}

//Locate the first non-space charachter of the string
long _String::FirstSpaceIndex(long start, long end, char direction) const {
    if (start == -1) {
        start = ((long)sLength)-1;
    }
    if (end == -1) {
        end = ((long)sLength)-1;
    }
    if (direction<0) {
        //long t = start;
        start = end;
        end = start;
    }
    if (sLength&&(isspace (sData[start]))) {
        return start;    // first char is non-space
    }
  
    for (int i = start; i<=end; i+=direction) {
        if (isspace (sData[i])) {
          return i;
        }
    }
    return -1;
}


long _String::FirstNonSpaceFollowingSpace(long start, long end, char direction) const {
  long first_space = FirstSpaceIndex (start, end, direction);
  if (first_space >= 0) {
    first_space = FirstNonSpaceIndex(first_space, end, direction);
  }
  return first_space;
}



//Remove all spaces
void _String::KillSpaces (_String& result)
{
    _String temp (sLength+1,true);
    for (long k=0; k<sLength; k++)
        if (!isspace (sData[k])) {
            temp << sData[k];
        }
    temp.Finalize();
    result = temp;
}

//Cut string from, to (-1 for any means from beginning/to end)
void _String::Trim(long from, long to, bool softTrim)
{
    if (!sLength) {
        return;
    }
    if (from < 0) {
        from = 0;
    } else if (from>=sLength) {
        from = ((long)sLength)-1;
    }
    if (to < 0) {
        to = ((long)sLength)-1;
    } else if (to>=sLength) {
        to = ((long)sLength)-1;
    }

    if (softTrim) {
        sData += from;
        sLength = to-from+1;
    } else if (to-from+1>0) {
        if (from) {
            memmove (sData,sData+from,  to-from+1);
        }

        sLength = to-from+1;
        sData = (char*)MemReallocate (sData, to-from+2);
        sData[to-from+1]=0;
    } else {
        sLength = 0;
        sData = (char*)MemReallocate (sData, 1);
        sData [0] = 0;
    }
}

/*
==============================================================
Lexicographic Comparison Methods
==============================================================
*/

bool _String::contains (_String s)
{
    return Find(s)!=-1;
}

bool _String::contains (char c)
{
    return Find(c)!=-1;
}



/*
==============================================================
Begins and Ends With Methods
==============================================================
*/

//Begins with string
bool _String::beginswith (_String const s, bool caseSensitive) const
{
    if (sLength<s.sLength) {
        return FALSE;
    }
    if (caseSensitive) {
        for (unsigned long i=0UL; i<s.sLength; i++)
            if (s.sData[i]!=sData[i]) {
                return FALSE;
            }
    } else
        for (long i=0; i<s.sLength; i++)
            if (toupper(s.sData[i])!=toupper(sData[i])) {
                return FALSE;
            }


    return TRUE;
}

//Begins with string
bool _String::startswith (_String const& s) const
{
    if (sLength<s.sLength) {
        return FALSE;
    }

    char *sP  = sData,
          *ssP = s.sData;

    for (; *ssP; sP++,ssP++)
        if (*sP!=*ssP) {
            return false;
        }

    return true;
}

//Begins with string
bool _String::startswith_noident (_String const& s) const
{
  
  if (startswith (s)) {
    if (sLength > s.sLength) {
      char next_char = get_char (s.sLength);
      //printf ("Next char %c (%d)\n", next_char, next_char);
      if (isalnum(next_char) || next_char == '.' || next_char == '_' || next_char == '&') {
        return false;
      }
    }
    return true;
  }
  return false;
}


//Ends with string
bool _String::endswith (_String s, bool caseSensitive)
{
    if (sLength<s.sLength) {
        return FALSE;
    }
    char *sP = sData+sLength-s.sLength,
          *ssP = (s.sData),
           *ssP2 = s.sData+s.sLength;

    if (caseSensitive) {
        for (; ssP!=ssP2; ssP++,sP++)
            if (*sP-*ssP) {
                return FALSE;
            }
    } else
        for (; ssP!=ssP2; ssP++,sP++)
            if (toupper(*sP)!=toupper(*ssP)) {
                return FALSE;
            }

    return TRUE;
}

/*
==============================================================
Case Methods
==============================================================
*/

void    _String::UpCase (void)
{
    for (unsigned long i = 0; i<sLength; i++) {
        sData[i] = toupper (sData[i]);
    }
}

void    _String::LoCase (void)
{
    for (unsigned long i = 0; i<sLength; i++) {
        sData[i] = tolower (sData[i]);
    }
}

void    _String::ProcessParameter(void)
{
    if (Equal(&getDString)) {
        *this = ReturnDialogInput();
    }
}

//==============================================================
//Filename and Platform Methods
//==============================================================

bool    _String::ProcessFileName (bool isWrite, bool acceptStringVars, hy_pointer theP, bool assume_platform_specific, _ExecutionList * caller)
{
    _String errMsg;
    
    try {
        if (Equal(&getFString) || Equal (&tempFString)) { // prompt user for file
            if (Equal (&tempFString)) {
                #if not defined __MINGW32__ && not defined __WINDOZE__
                    #ifdef __MAC__
                        char tmpFileName[] = "HYPHY-XXXXXX";
                    #else
                        char tmpFileName[] = "/tmp/HYPHY-XXXXXX";
                    #endif
                    
                    int fileDescriptor = mkstemp(tmpFileName);
                    if (fileDescriptor == -1){
                        throw ("Failed to create a temporary file name");
                    }
                    *this = tmpFileName;
                    CheckReceptacleAndStore(&hy_env::last_file_path,kEmptyString,false, new _FString (*this, false), false);
                    close (fileDescriptor);
                    return true;
                #else
                    throw (tempFString & " is not implemented for this platform");
                #endif
            } else {
                if (!isWrite) {
                    *this = ReturnFileDialogInput();
                } else {
                    *this = WriteFileDialogInput ();
                }
            }
            ProcessFileName(false,false,theP,
            #if defined __MAC__ || defined __WINDOZE__
                true
            #else
                false
            #endif
            ,caller);
            
            CheckReceptacleAndStore(&hy_env::last_file_path,kEmptyString,false, new _FString (*this, false), false);
            return true;
        }

        if (acceptStringVars) {
            *this = ProcessLiteralArgument (this,(_VariableContainer*)theP, caller);
            if (caller && caller->IsErrorState()) {
                return false;
            }
 
        } else {
            StripQuotes();
        }

        if (!sLength) {
           CheckReceptacleAndStore(&hy_env::last_file_path,kEmptyString,false, new _FString (*this, false), false);
           return true;
        }
    }
    
    catch (_String errmsg) {
        if (caller) {
            caller->ReportAnExecutionError(errMsg);
        } else {
            HandleApplicationError (errMsg);
        }
        return false;
    }


#if (defined __UNIX__ || defined __HYPHY_GTK__) && !defined __MINGW32__
//UNIX LINES HERE
    if (Find('\\')!=-1) { // DOS (ASSUME RELATIVE) PATH
        *this = Replace ("\\","/",true);
    } else if (Find(':')!=-1) { // Mac (Assume Relative) PATH
        *this = Replace ("::",":../", true);
        if (get_char(0)==':') {
            Trim(1,-1);
        }
        *this = Replace (':','/',true);
    }

    if (get_char(0) != '/') { // relative path
        if (pathNames.lLength) {
            _String*    lastPath = (_String*)pathNames(pathNames.lLength-1);
            long        f = lastPath->sLength-2,
                        k = 0;

            // check the last stored absolute path and reprocess this relative path into an absolute.
            while (beginswith("../")) {
                if ( (f = lastPath->FindBackwards('/',0,f)-1) ==-1) {
                    return true;
                }
                Trim(3,-1);
                k++;
            }
            if (k==0) {
                *this = *lastPath& (*this);
            } else {
                *this = lastPath->Cut(0,f+1)& (*this);
            }
        } 
    }
#endif

#if defined __WINDOZE__ || defined __MINGW32__ // WIN/DOS code'
  
    if (Find('/')!=-1) { // UNIX PATH
        if (get_char(0)=='/') {
            Trim(1,-1);
        }
        *this = Replace ("/","\\",true);
    } else {
        if (Find('\\')==-1) {
            // check to see if this is a relative path
            *this = Replace ("::",":..\\", true);
            if ((sData[0]==':')) {
                Trim(1,-1);
            }
            *this = Replace (':','\\',true);
        }
    }

  if (Find(':') < 0 && Find("\\\\",0,1)==-1) { // relative path

        if (pathNames.lLength) {
            _String* lastPath = (_String*)pathNames(pathNames.lLength-1);
            long f = lastPath->sLength-2, k = 0;
            // check the last stored absolute path and reprocess this relative path into an absolute.
            while (beginswith("..\\")) {
                f = lastPath->FindBackwards('\\',0,f)-1;
                if (f==-1) {
                    return false;
                }
                Trim(3,-1);
                k++;
            }
            if (k==0) {
                if (lastPath->sData[lastPath->sLength-1]!='\\') {
                    *this = *lastPath&'\\'& (*this);
                } else {
                    *this = *lastPath& (*this);
                }
            } else {
                *this = lastPath->Cut(0,f+1)& (*this);
            }
        } 

    }



   _String escapedString (sLength, true);
    for (long stringIndex = 0; stringIndex < sLength; stringIndex ++) {
        char currentChar = get_char (stringIndex);
        //char b[256];
        //snprintf (b, sizeof(b),"%c %d\n", currentChar, currentChar);
        //BufferToConsole (b);
        switch (currentChar) {
        case '\t':
            escapedString << '\\';
            escapedString << 't';
            break;
        case '\n':
            escapedString << '\\';
            escapedString << 'n';
            break;
        default:
            escapedString << currentChar;
        }
    }
    escapedString.Finalize();
    (*this) = escapedString;

#endif

#ifdef __MAC__
    if (!assume_platform_specific && Find('/')!=-1) { // UNIX PATH
        bool rootPath = false;
        if (sData[0]=='/') {
            rootPath = true;
            *this = volumeName&Cut(1,-1);
        }

        if (beginswith("..")) {
            *this = _String('/')&Cut(2,-1);
        }

        *this = Replace ("/",":",true);
        *this = Replace ("..","",true);

        if (sData[0]!=':' && !rootPath) {
            *this = _String(':')&*this;
        }
    } else {
        if (!assume_platform_specific && Find('\\')!=-1) { // DOS PATH (ASSUME PARTIAL)
            if (beginswith("..")) {
                *this = _String('\\')&Cut(2,-1);
            }
            *this = Replace ("\\",":",true);
            *this = Replace ("..","",true);
            if (Find(':')!=-1) {
                *this = _String(':')&*this;
            }
        } else { // MAC PATH
            if (Find(':')!=-1) {
                if (sData[0]!=':') {
                    if (!beginswith(volumeName)) {
                        if (pathNames.lLength) {
                            _String* lastPath = (_String*)pathNames(pathNames.lLength-1);
                            if (!beginswith (lastPath->Cut (0,lastPath->Find(':')))) {
                                *this = _String(':')&*this;
                            }
                        } else {
                            *this = _String(':')&*this;
                        }
                    }
                }
            } else {
                *this = _String(':')&*this;
            }
        }
    }

    if (sData[0]==':') { // relative path
        long f = -1, k = 0;
        if (pathNames.lLength) {
            _String* lastPath = (_String*)pathNames(pathNames.lLength-1);
            // check the last stored absolute path and reprocess this relative path into an absolute.
            while (sData[k]==':') {
                f = lastPath->FindBackwards(':',0,f)-1;
                if (f==-1) {
                    return true;
                }
                k++;
            }
            *this = lastPath->Cut(0,f+1)& Cut(k,-1);
        } else {
            *this = kEmptyString;
        }
    }
#endif
    CheckReceptacleAndStore(&hy_env::last_file_path,kEmptyString,false, new _FString (*this, false), false);
    return true;
}

//Compose two UNIX paths (abs+rel)
_String const _String::PathComposition (_String const relPath) const {
    if (relPath.sData[0]!='/') { // relative path
        long f = -1L, k = 0;
        f = sLength-2;
        _String result = *this,
        relative_path = relPath;

        while (relative_path.beginswith("../")) {

            //Cut Trim relPath
            f = FindBackwards('/',0,f)-1;

            relative_path = relative_path.Chop(0,2);
            result.Trim(0,f+1);

            if (f==-1) {
                return kEmptyString;
            }
            k++;

        }

        return result&relative_path;
    }

    else {
        return relPath;
    }
    return kEmptyString;
}

//Mac only so far
_String const _String::PathSubtraction (_String const p2, char) const {
    _String result;
    char separator = get_platform_directory_char();

    //if (pStyle == 0)
    //    separator = ':';
    long k;
    for (k=0; (k<sLength)&&(k<p2.sLength)&&(sData[k]==p2.sData[k]); k++) ;
    if (k>0) {
        while (sData[k]!=separator) {
            k--;
        }
        if (k>0) {
            long m=k+1, levels = 0;
            for (; m<sLength; m++)
                if (sData[m]==separator) {
                    levels++;
                }
            if (levels) {
                result = separator;
                while (levels) {
                    result.Insert (separator,-1);
                    levels--;
                }
            }
            result = result & p2.Cut(k+1,-1);
            return result;
        }
    }
    return kEmptyString;
}

//TODO: These are global methods. Should they even be here?


_String GetVersionString (void)
{
    _String theMessage = _String("HYPHY ")&__HYPHY__VERSION__;
#ifdef __MP__
    theMessage = theMessage & "(MP)";
#endif
#ifdef __HYPHYMPI__
    theMessage = theMessage & "(MPI)";
#endif
    theMessage = theMessage & " for ";
#ifdef __MAC__
    theMessage = theMessage & "MacOS";
#ifdef TARGET_API_MAC_CARBON
    theMessage = theMessage & "(Carbon)";
#endif
#endif
#ifdef __WINDOZE__
    theMessage = theMessage & "Windows (Win32)";
#endif
#ifdef __UNIX__
#if !defined __HEADLESS_WIN32__ && ! defined __MINGW32__
    struct      utsname      name;
    uname       (&name);
    theMessage = theMessage & name.sysname & " on " & name.machine;
#endif
#if defined __MINGW32__
    theMessage = theMessage & "MinGW ";// " & __MINGW32_VERSION;
#endif
#endif
    return theMessage;
}

_String GetTimeStamp (bool doGMT)
{
    time_t cTime;
    time (&cTime);

    if (doGMT) {
        tm* gmt = gmtime (&cTime);
        return _String ((long)1900+gmt->tm_year) & '/' & _String (1+(long)gmt->tm_mon) & '/'
               & _String ((long)gmt->tm_mday) & ' ' & _String ((long)gmt->tm_hour) & ':' & _String ((long)gmt->tm_min);
    }

    tm*     localTime = localtime (&cTime);

    return  asctime (localTime);

}

/*
==============================================================
Identifier Methods
==============================================================
*/

bool    _String::IsValidIdentifier (bool strict) const
{
    if (sLength == 0) { 
        return false;
    }    

    if (strict) {
        if (!(isalpha(sData[0]) || sData[0]=='_' )) { 
            return false;
        }    
    } else if (!(isalnum(sData[0]) || sData[0]=='_' )) { 
        return false;
    }    


    for(unsigned long p = 1; p<sLength; p++) {
        char c = sData[p];
        if (!(isalnum(c)|| c=='_' || (strict && c == '.'))) {
            return false;
        }    
    }    

    // check to see if it's not a keyword / function name etc

    return hyReservedWords.FindObject (this) == -1;
}

bool    _String::IsValidRefIdentifier (void) const
{
    if (sLength<2) {
        return false;
    }
    if (sData[sLength-1]=='&') {
        return Cut(0,sLength-2).IsValidIdentifier();
    }
    return false;
}

//Convert a string to a valid ident
void _String::ConvertToAnIdent (bool strict)
{
    _String * result = new _String ((unsigned long)sLength+1,true);

    if (sLength) {
        if (strict) {
            if (((sData[0]>='a')&&(sData[0]<='z'))||((sData[0]>='A')&&(sData[0]<='Z'))||(sData[0]=='_')) {
                (*result)<<sData[0];
            } else {
                (*result)<<'_';
            }
        } else {
            if (((sData[0]>='a')&&(sData[0]<='z'))||((sData[0]>='A')&&(sData[0]<='Z'))||(sData[0]=='_')||((sData[0]>='0')&&(sData[0]<='9'))) {
                (*result)<<sData[0];
            } else {
                (*result)<<'_';
            }
        }

        long l = 0;
        for (long k=1; k<sLength; k++) {
            unsigned char c = sData[k];
            if (_hyValidIDChars.valid_chars[c]) {
                (*result)<<c;
                l++;
            } else if (result->sData[l] != '_') {
                (*result)<<'_';
                l++;
            }
        }
    }
    result->Finalize();

    CopyDynamicString (result, true);
}

_String _String::ShortenVarID (_String& containerID)
{
    long matched=-1,
         upTo = sLength<containerID.sLength?sLength:containerID.sLength,
         k;

    for (k=0; k<upTo; k++) {
        if (sData[k]!=containerID.sData[k]) {
            break;
        } else if (sData[k] == '.') {
            matched = k;
        }
    }

    if ((upTo==containerID.sLength)&&(upTo<sLength)&&(k==upTo)&&(sData[upTo]=='.')) {
        matched = upTo;
    }

    return Cut (matched+1,-1);
}

/*
==============================================================
Regular Expression Methods
==============================================================
*/

_String GetRegExpError(int error)
{
    char buffer [512];
    buffer[regerror (error, nil, buffer, 511)] = 0;
    return _String("Regular Expression error:")&buffer;
}

void FlushRegExp(hy_pointer regExpP)
{
    regex_t*        regEx = (regex_t*)regExpP;
    regfree        (regEx);
    delete          regEx;
}

hy_pointer PrepRegExp(_String* source, int& errCode, bool caseSensitive)
{
    regex_t  * res = new regex_t;

    errCode = regcomp (res, source->sData, REG_EXTENDED|(caseSensitive?0:REG_ICASE));

    if (errCode) {
        FlushRegExp ((hy_pointer)res);
        return nil;
    }
    return (hy_pointer)res;
}

void _String::RegExpMatch(hy_pointer pattern, _SimpleList& matchedPairs)
{
    if (sLength) {
        regex_t*        regEx = (regex_t*)pattern;

        regmatch_t*     matches = new regmatch_t [regEx->re_nsub+1];
        int             errNo = regexec (regEx, sData,regEx->re_nsub+1, matches, 0);
        if (errNo == 0) {
            for (long k=0; k<=regEx->re_nsub; k++) {
                matchedPairs << matches[k].rm_so;
                matchedPairs << matches[k].rm_eo-1;
            }
        }
        delete      []  matches;
    }
}

void _String::RegExpMatchAll(hy_pointer pattern, _SimpleList& matchedPairs)
{
    if (sLength) {
        regex_t*        regEx = (regex_t*)pattern;

        regmatch_t*     matches = new regmatch_t [regEx->re_nsub+1];
        int             errNo =  regexec (regEx, sData,regEx->re_nsub+1, matches, 0);
        while (errNo == 0) {
            long         offset = matchedPairs.lLength?matchedPairs.lData[matchedPairs.lLength-1]+1:0;

            matchedPairs << matches[0].rm_so+offset;
            matchedPairs << matches[0].rm_eo-1+offset;

            offset += matches[0].rm_eo;
            if (offset < sLength) {
                errNo =  regexec (regEx, sData+offset ,regEx->re_nsub+1, matches, 0);
            } else {
                break;
            }
        }
        delete  []      matches;
    }
}

void _String::RegExpMatchOnce(_String* pattern, _SimpleList& matchedPairs, bool caseSensitive, bool handleErrors)
{
    if (sLength) {
        int errNo = 0;
        hy_pointer regex = PrepRegExp (pattern, errNo, caseSensitive);
        if (regex) {
            RegExpMatch (regex, matchedPairs);
            FlushRegExp (regex);
        } else if (handleErrors) {
            HandleApplicationError (GetRegExpError (errNo));
        }
    }
}

_String _String::Random(const unsigned long length, const _String * alphabet)
{
    _String random (length + 1, true);
    unsigned long alphabet_length = alphabet?alphabet->sLength:127;
    if (length > 0 && alphabet_length > 0) {
        for (unsigned long c = 0; c < length; c++) {
            unsigned long idx = genrand_int32 () % alphabet_length;
            if (alphabet) {
                random << alphabet->sData[idx];
            } else {
                random << (char)(1+idx);
            }
        }
    }
    
    random.Finalize();
    return random;
}

void    _String::AppendNCopies   (_String const& value, unsigned long copies) {
  for (unsigned long i = 0UL; i < copies; i++) {
    (*this) << value;
  }
}


unsigned char _String::ProcessVariableReferenceCases (_String& referenced_object, _String const * context) const {
    char first_char    = get_char(0);
    bool is_func_ref  = get_char(sLength-1) == '&';
         
    if (first_char == '*' || first_char == '^') {
        if (is_func_ref) {
            referenced_object = kEmptyString;
            return kStringInvalidReference;
        }
        bool is_global_ref = first_char == '^';
        _String   choppedVarID (*this, 1, -1);
      
        if (choppedVarID.IsValidIdentifier()) {
          if (context) {
              choppedVarID = *context & '.' & choppedVarID;
          }
          _FString * dereferenced_value = (_FString*)FetchObjectFromVariableByType(&choppedVarID, STRING);
          if (dereferenced_value && dereferenced_value->theString->ProcessVariableReferenceCases (referenced_object) == kStringDirectReference) {
              if (!is_global_ref && context) {
                  referenced_object = *context & '.' & referenced_object;
              }
              return is_global_ref?kStringGlobalDeference:kStringLocalDeference;
          }
        } else {
          
          _String try_as_expression;
          if (context) {
            _VariableContainer ctxt (*context);
            try_as_expression = ProcessLiteralArgument (&choppedVarID, &ctxt);
          } else {
            try_as_expression = ProcessLiteralArgument (&choppedVarID, nil);
          }
           if (try_as_expression.ProcessVariableReferenceCases (referenced_object) == kStringDirectReference) {
             if (!is_global_ref && context) {
              referenced_object = *context & '.' & try_as_expression;
            }
             
            return is_global_ref?kStringGlobalDeference:kStringLocalDeference;
          }
        }
    }
    
    if (is_func_ref) {
        referenced_object = Cut (0, sLength-2);
        if (referenced_object.IsValidIdentifier()) {
            referenced_object = (context? (*context & '.' & referenced_object): (referenced_object)) & '&';
            return kStringDirectReference;
        }    
    }
    else {
        if (IsValidIdentifier()) {
          if (context) {
            _String cdot = *context & '.';
            referenced_object = startswith(cdot) ? *this : (cdot & *this);
          } else {
            referenced_object = *this;
          }
          return kStringDirectReference;
        }
    }
    
    referenced_object = kEmptyString;
    return kStringInvalidReference;
}


///// MOVE TO STRING BUFFER

// append operator
_String& _String::operator << (const _String* s) {
    if ( s && s->sLength) {
        if (buffer_allocation < sLength + s->sLength) {
            unsigned long incBy = sLength + s->sLength - buffer_allocation;
            
            if (incBy < storageIncrement) {
                incBy = storageIncrement;
            }
            
            if (incBy < sLength >> 3) {
                incBy = sLength >> 3;
            }
            
            buffer_allocation+=incBy;
            
            if (! (sData = (char*)MemReallocate((char*)sData, buffer_allocation*sizeof(char)))) {
                return *this;
            }
            
        }
        
        for (unsigned long k = 0UL; k < s->sLength; k++) {
            sData[sLength+k] = s->sData[k];
        }
        
        //memcpy(sData+sLength,s->sData,s->sLength);
        sLength+=s->sLength;
    }
    
    return *this;
}

// append operator
_String& _String::operator << (const _String& s) {
    return (*this) << &s;
}

//Append operator
_String& _String::operator << (const char* str) {
    _String conv (str);
    return (*this)<<&conv;
}

//Append operator
_String& _String::operator << (const char c) {
    if (buffer_allocation <= sLength) {
        buffer_allocation  += ((storageIncrement*8 > sLength)? storageIncrement: (sLength/8+1));
        sData = (char*)MemReallocate((char*)sData, buffer_allocation*sizeof(char));
    }
    
    sData[sLength++]=c;
    return *this;
}

//Append operator
void _String::EscapeAndAppend (const char c, char mode)
{
    if (mode == 2) {
        (*this) << c;
        switch (c) {
            case '\'':
                (*this) << c;
        }
        return;
    } else {
        if (mode == 1) {
            switch (c) {
                case '(':
                case ')':
                case '%':
                    (*this) << '\\';
                    (*this) << c;
                    return;
            }
        } else {
            if (mode == 4) {
                switch (c) {
                    case '"':
                        (*this) << "&quot;";
                        break;
                    case '\'':
                        (*this) << "&apos;";
                        break;
                    case '<':
                        (*this) << "&lt;";
                        break;
                    case '>':
                        (*this) << "&gt;";
                        break;
                    case '&':
                        (*this) << "&amp;";
                        break;
                    default:
                        (*this) << c;
                }
                return;
            } else {
                if (mode == 5) { // regexp
                    switch (c) {
                        case '[':
                        case '^':
                        case '$':
                        case '.':
                        case '|':
                        case '?':
                        case '*':
                        case '+':
                        case '(':
                        case ')':
                            (*this) << '\\';
                            (*this) << c;
                            break;
                        case '\\':
                            (*this) << "\\\\";
                            break;
                        default:
                            (*this) << c;
                    }
                    return;
                    
                }
            }
        }
    }
    switch (c) {
        case '\n':
            (*this) << '\\';
            (*this) << 'n';
            break;
        case '\t':
            (*this) << '\\';
            (*this) << 't';
            break;
        case '"':
            (*this) << '\\';
            (*this) << '"';
            break;
        case '\\':
            (*this) << '\\';
            (*this) << '\\';
            break;
        default:
            (*this) << c;
    }
}

//Append operator
void _String::EscapeAndAppend (const _String & s, char mode)
{
    for (long i=0; i<s.sLength;  i++) {
        EscapeAndAppend (s.sData[i], mode);
    }
}

//Append and delete operator
void _String::AppendNewInstance (_String* s) {
    (*this) << s;
    DeleteObject (s);
}

void _String::AppendAnAssignmentToBuffer(_String* id, _String *value, unsigned long flags)
{
    
    if (flags & kAppendAnAssignmentToBufferGlobal) {
        (*this) << "global ";
    }
    
    (*this) << id;
    
    if (flags & kAppendAnAssignmentToBufferAssignment) {
        (*this) << ':';
    }
    
    (*this) << '=';
    
    if (flags & kAppendAnAssignmentToBufferQuote) {
        (*this) << '"';
    }
    
    (*this) << value;
    
    if (flags & kAppendAnAssignmentToBufferQuote) {
        (*this) << '"';
    }
    
    (*this) << ";\n";
    
    if (flags & kAppendAnAssignmentToBufferFree) {
        DeleteObject (value);
    }
    
}


void _String::AppendVariableValueAVL (_String* id, _SimpleList& varNumbers)
{
    for (long k=0; k<varNumbers.lLength; k++) {
        _Variable *tiv = LocateVar(varNumbers.lData[k]);
        if (tiv) {
            (*this) << id;
            (*this) << "[\"";
            (*this) << tiv->GetName();
            (*this) << "\"]=";
            _PMathObj varValue = tiv->Compute();
            switch (varValue->ObjectClass()) {
                case NUMBER:
                    (*this) << _String (varValue->Value());
                    break;
                case STRING:
                    (*this) << '"';
                    EscapeAndAppend (*((_FString*)varValue)->theString);
                    (*this) << '"';
                    break;
                default:
                    AppendNewInstance ((_String*)(varValue->toStr()));
                    break;
                    
            }
            (*this) << ";\n";
        }
    }
}


