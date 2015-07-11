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

#include "baseobj.h"
#include "errorfns.h"
#include "hy_strings.h"
#include "site.h"
#include "batchlan.h"
#include "category.h"
#include "likefunc.h"

#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef   __HYPHYMPI__
#include "likefunc.h"
extern int _hy_mpi_node_rank;
#endif

#if defined   __UNIX__ || defined __HYPHY_GTK__
#include <sys/time.h>
#include <unistd.h>
#endif

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

#ifdef   __HYPHYXCODE__
#include "HYUtils.h"
#endif

#ifdef __WINDOZE__
    #include <Windows.h>
#endif



bool        terminateExecution  = false;

#include    "batchlan.h"

//____________________________________________________________________________________

FILE*     globalErrorFile   = nil,
          *   globalMessageFile = nil;

extern    bool          isInFunction;

extern    _String   scanfLastFilePath;

extern    _SimpleList   freeSlots;


//____________________________________________________________________________________

long            globalRandSeed;

//____________________________________________________________________________________

_String         errorFileName   ("errors.log"),
                messageFileName ("messages.log");


//____________________________________________________________________________________

BaseObj::BaseObj()
{
    nInstances=1;
}


//____________________________________________________________________________________
BaseRef   BaseObj::toStr (void)
{
    return new _String ("null");
}

//____________________________________________________________________________________
BaseRef   BaseObj::toErrStr (void)
{
    return toStr();
}

//____________________________________________________________________________________
void     BaseObj::toFileStr (FILE* dest)
{
    _String* s = (_String*)toStr();
    fwrite(s->sData,1,s->Length(),dest);
    DeleteObject (s);
}

//____________________________________________________________________________________
BaseObj*  BaseObj::makeDynamic(void)
{
    warnError(-112);
    return nil;
}

//____________________________________________________________________________________

FILE *      doFileOpen (const char * fileName, const char * mode, bool warn)
{
    FILE    *daFile = nil;

    if (fileName) {
        _String fnstr (fileName);
#ifdef __HYPHYXCODE__
        daFile = fopen (DoMacToPOSIX(fnstr).getStr(),mode);
#else
        daFile = fopen (fileName,mode);
#endif
        if (!daFile && warn) {
            WarnError (_String("Could not open file '") & *fileName & "' with mode '" & *mode & "'.");
        }
    }
    return daFile;
}


//____________________________________________________________________________________
bool    GlobalStartup (void)
{
    SetupOperationLists     ();
    unsigned long seed_init = 0L;
    time_t                  k;
    time                    (&k);
    seed_init               = k;
    
#ifdef __HYPHYMPI__
    _hyApplicationGlobals.Insert(new _String (mpiNodeID));
    _hyApplicationGlobals.Insert(new _String (mpiNodeCount));
    _hyApplicationGlobals.Insert(new _String (mpiLastSentMsg));
    seed_init += _hy_mpi_node_rank;
#endif

    
#ifndef __WINDOZE__
    seed_init               += getpid();
#else
    seed_init               += long(GetCurrentProcess());
#endif
    init_genrand            (seed_init);
    globalRandSeed          = seed_init;
    setParameter            (randomSeed,globalRandSeed);
    long                    p   = 1;

    _hyApplicationGlobals.Insert(new _String (dataFileTree));
    _hyApplicationGlobals.Insert(new _String (dataFileTreeString));
    _hyApplicationGlobals.Insert(new _String (siteWiseMatrix));
    _hyApplicationGlobals.Insert(new _String (blockWiseMatrix));
    _hyApplicationGlobals.Insert(new _String (selectionStrings));
    _hyApplicationGlobals.Insert(new _String (randomSeed));
    _hyApplicationGlobals.Insert(new _String (statusBarUpdateString));
    _hyApplicationGlobals.Insert(new _String (statusBarProgressValue));
    _hyApplicationGlobals.Insert(new _String (hyphyBaseDirectory));
    _hyApplicationGlobals.Insert(new _String (hyphyLibDirectory));
    _hyApplicationGlobals.Insert(new _String (platformDirectorySeparator));
    _hyApplicationGlobals.Insert(new _String (pathToCurrentBF));

    _String             dd (GetPlatformDirectoryChar());

    standardLibraryPaths.AppendNewInstance      (new _String(libDirectory & "TemplateBatchFiles" & dd));
    standardLibraryPaths.AppendNewInstance      (new _String(libDirectory & "TemplateBatchFiles" & dd & "TemplateModels" & dd ));
    standardLibraryPaths.AppendNewInstance      (new _String(libDirectory & "TemplateBatchFiles" & dd & "Utility" & dd));
    standardLibraryPaths.AppendNewInstance      (new _String(libDirectory & "UserAddIns" & dd));
    standardLibraryPaths.AppendNewInstance      (new _String(libDirectory & "TemplateBatchFiles" & dd & "Distances" & dd));

    standardLibraryExtensions.AppendNewInstance (new _String (""));
    standardLibraryExtensions.AppendNewInstance (new _String (".bf"));
    standardLibraryExtensions.AppendNewInstance (new _String (".ibf"));
    standardLibraryExtensions.AppendNewInstance (new _String (".def"));
    standardLibraryExtensions.AppendNewInstance (new _String (".mdl"));

    _HBL_Init_Const_Arrays  ();


#if not defined (__HYPHY_MPI_MESSAGE_LOGGING__) && defined (__HYPHYMPI__)
    if (_hy_mpi_node_rank == 0) {
#endif



#ifndef __HEADLESS__ // do not create log files for _HEADLESS_
#ifndef __HYPHYMPI__
    _String fileName(errorFileName);
#if defined __HYPHYXCODE__ || defined __WINDOZE__
    fileName = baseDirectory & fileName;
#endif
#else
    _String fileName = errorFileName & ".mpinode" & (long)_hy_mpi_node_rank;
#endif

    globalErrorFile = doFileOpen (fileName.sData,"w+");
    while (globalErrorFile == nil && p<10) {
        fileName = errorFileName&'.'&_String(p);
#if defined __HYPHYXCODE__ || defined __WINDOZE__
        fileName = baseDirectory & fileName;
#endif
        globalErrorFile = doFileOpen (fileName.sData,"w+");
        p++;
    }
    errorFileName = fileName;

    p=1;
#ifndef __HYPHYMPI__
    fileName = messageFileName;
#if defined __HYPHYXCODE__ || defined __WINDOZE__
    fileName = baseDirectory & fileName;
#endif
#else
    fileName = messageFileName & ".mpinode" & (long)_hy_mpi_node_rank;
#endif

    globalMessageFile = doFileOpen (fileName.getStr(),"w+");

    while (globalMessageFile == nil && p<10) {
        fileName = messageFileName&'.'&_String(p);
#if defined __HYPHYXCODE__ || defined __WINDOZE__
        fileName = baseDirectory & fileName;
#endif
        globalMessageFile = doFileOpen (fileName.sData,"w+");
        p++;
    }
    messageFileName = fileName;
#endif
#if not defined (__HYPHY_MPI_MESSAGE_LOGGING__) && defined (__HYPHYMPI__)
    }
#endif

    return globalErrorFile && globalMessageFile;
}


//____________________________________________________________________________________
bool    GlobalShutdown (void)
{
    bool res = true;

#if defined (__HYPHY_GTK__) || defined (__HEADLESS__)
#elif defined __UNIX__
    if (needExtraNL)
        printf ("\n");
#endif

#ifdef  __HYPHYMPI__
    int     flag, size;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    PurgeAll (true);
    // force manual clear to help debuggin 'exit' crashes
    ReportWarning ("PurgeAll was successful");
    if (_hy_mpi_node_rank == 0) {
        for (long count = 1; count < size; count++) {
            ReportWarning (_String ("Sending shutdown command to node ") & count & '.');
            MPISendString(empty,count);
        }
    }
#endif

#ifdef  __HYPHYMPI__
    // MPI_Barrier (MPI_COMM_WORLD);
    ReportWarning ("Calling MPI_Finalize");
#ifdef __USE_ABORT_HACK__
    MPI_Abort(MPI_COMM_WORLD,0);
#else
    MPI_Finalized(&flag);
    if (!flag)
        MPI_Finalize();
#endif
    ReportWarning ("Returned from MPI_Finalize");
#endif


    if (globalErrorFile) {
        fflush (globalErrorFile);
        fseek(globalErrorFile,0,SEEK_END);
        unsigned long fileSize = ftell(globalErrorFile);
        if (fileSize) {
#if defined (__MAC__) || defined (__WINDOZE__) || defined (__HYPHYMPI__) || defined (__HYPHY_GTK__)

#else
            fprintf (stderr, "\nCheck %s for details on execution errors.\n",errorFileName.getStr());
#endif
            res = false;
            fclose (globalErrorFile);

        } else {
            fclose (globalErrorFile);
#ifdef __HYPHYXCODE__
            remove (DoMacToPOSIX(errorFileName).getStr());
#else
            remove (errorFileName.getStr());
#endif
        }
    }
    if (globalMessageFile) {
        if (ftell(globalMessageFile)) {
#if defined (__MAC__) || defined (__WINDOZE__) || defined (__HYPHYMPI__) || defined (__HYPHY_GTK__)
#else
            fprintf (stderr, "\nCheck %s details of this run.\n",messageFileName.getStr());
#endif
            fclose (globalMessageFile);
        } else {
            fclose (globalMessageFile);
#ifdef __HYPHYXCODE__
            remove (DoMacToPOSIX(messageFileName).getStr());
#else
            remove (messageFileName.getStr());
#endif
        }
    }
    
    _SimpleList  hist;
    long         ls,
                 cn = _HY_HBLCommandHelper.Traverser (hist,ls,_HY_HBLCommandHelper.GetRoot());

    while (cn >= 0) {
        delete ((_HBLCommandExtras*)_HY_HBLCommandHelper.GetXtra(cn));
        cn = _HY_HBLCommandHelper.Traverser (hist,ls);
    }
    _HY_HBLCommandHelper.Clear();
    _HY_ValidHBLExpressions.Clear();

    return res;
}

//____________________________________________________________________________________

void    PurgeAll (bool all)
{
    batchLanguageFunctions.Clear();
    batchLanguageFunctionNames.Clear();
    batchLanguageFunctionParameterLists.Clear();
    batchLanguageFunctionParameters.Clear();
    batchLanguageFunctionClassification.Clear();
    executionStack.Clear();
    loadedLibraryPaths.Clear(true);
    _HY_HBL_Namespaces.Clear();
    if (all) {
        likeFuncList.Clear();
        likeFuncNamesList.Clear();
        dataSetFilterList.Clear();
        dataSetFilterNamesList.Clear();
        dataSetList.Clear();
        dataSetNamesList.Clear();
        compiledFormulaeParameters.Clear();
        modelNames.Clear();
        KillExplicitModelFormulae ();
        modelMatrixIndices.Clear();
        modelFrequenciesIndices.Clear();
        modelTypeList.Clear();
        listOfCompiledFormulae.Clear();
        variablePtrs.Clear();
        freeSlots.Clear();
        lastMatrixDeclared = -1;
        {
            variableNames.Clear(true);
        }
        _x_ = nil;
        _n_ = nil;
        pathNames.Clear();
    }
    scanfLastFilePath = empty;
    setParameter (randomSeed,globalRandSeed);
    isInFunction        = false;
    isDefiningATree     = 0;
#ifdef __HYPHYMPI__
    int            size;

    MPI_Comm_size   (MPI_COMM_WORLD, &size);
    setParameter    (mpiNodeID, (_Parameter)_hy_mpi_node_rank);
    setParameter    (mpiNodeCount, (_Parameter)size);
#endif
}

//____________________________________________________________________________________
void    DeleteObject (BaseRef theObject)
{
    if (theObject) {
        if (theObject->nInstances<=1) {
            delete (theObject);
        } else {
            theObject->nInstances--;
        }
    }
}


#ifndef __HYPHYDMALLOC__
//____________________________________________________________________________________
char* MemAllocate (long chunk)
{
    char* result = (char*)malloc (chunk);
    if (!result) {
        _String   errMsg ("Failed to allocate ");
        errMsg = errMsg & chunk & " bytes.";
        FlagError (errMsg);
    }
    return result;
}


//____________________________________________________________________________________
char* MemReallocate (Ptr oldP, long chunk)
{
    char* result  = (char*) realloc (oldP, chunk);

    if (!result) {
        warnError(-108);
    }

    return result;
}

#endif


//____________________________________________________________________________________


#ifdef __HEADLESS__
void        yieldCPUTime(void)
{
    if (globalInterfaceInstance) {
        terminateExecution = !(*globalInterfaceInstance->GetCallbackHandler())
                             (_THyPhyGetStringStatus(),_THyPhyGetLongStatus(),_THyPhyGetDoubleStatus());
    }
}

#endif 

//____________________________________________________________________________________

// time differencing function

double      TimerDifferenceFunction (bool doRetrieve)
{
    double timeDiff = 0.0;
#ifdef __MAC__
    static UnsignedWide microsecsIn;
    UnsignedWide microsecsOut;

    if (doRetrieve) {
        Microseconds      (&microsecsOut);
        timeDiff = 0.000001*((microsecsOut.hi-microsecsIn.hi)*(_Parameter)0xffffffff
                             + (microsecsOut.lo-microsecsIn.lo));
    } else {
        Microseconds      (&microsecsIn);
    }

#endif

#ifdef  __WINDOZE__
    static          char            canRunTimer    = 0;
    static          _Parameter      winTimerScaler = 0.0;
    static          LARGE_INTEGER   tIn;

    LARGE_INTEGER   tOut;

    if (canRunTimer == 0) {
        if (QueryPerformanceFrequency(&tIn) && tIn.QuadPart) {
            winTimerScaler = 1./tIn.QuadPart;
            canRunTimer    = 1;
        } else {
            canRunTimer = -1;
        }

    }

    if (canRunTimer == 1)
        if (doRetrieve) {
            QueryPerformanceCounter (&tOut);
            timeDiff   = (tOut.QuadPart-tIn.QuadPart) * winTimerScaler;

        } else {
            QueryPerformanceCounter (&tIn);
        }
#endif

#if !defined __WINDOZE__ && !defined __MAC__
    /*static    clock_t         clockIn;
            clock_t         clockOut;

    if (doRetrieve)
    {
        clockOut = clock();
        timeDiff   = (clockOut-clockIn) * 1.0 / CLOCKS_PER_SEC;
    }
    else
        clockIn  = clock();*/
    static timeval clockIn, clockOut;
    if (doRetrieve) {
        gettimeofday (&clockOut,nil);
        timeDiff = (clockOut.tv_sec-clockIn.tv_sec) + (clockOut.tv_usec-clockIn.tv_usec)*0.000001;
    } else {
        gettimeofday (&clockIn,nil);
    }

#endif


    return timeDiff;
}


//____________________________________________________________________________________
//       RANDOM NUMBER GENERATOR CODE
//____________________________________________________________________________________

/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/

#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* prototypes */

void        init_by_array       (unsigned long[], unsigned long);
long        genrand_int31       (void);
double      genrand_real1       (void);
double      genrand_real3       (void);
double      genrand_res53       (void);


/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
            (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void init_by_array(unsigned long init_key[], unsigned long key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1;
    j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
                + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        j++;
        if (i>=N) {
            mt[0] = mt[N-1];
            i=1;
        }
        if (j>=key_length) {
            j=0;
        }
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
                - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) {
            mt[0] = mt[N-1];
            i=1;
        }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]= {0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1) { /* if init_genrand() has not been called, */
            init_genrand(5489UL);    /* a default initial seed is used */
        }

        for (kk=0; kk<N-M; kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (; kk<N-1; kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0);
    /* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return ((double)genrand_int32())*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void)
{
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6;
    return(a*67108864.0+b)*(1.0/9007199254740992.0);
}
/* These real versions are due to Isaku Wada, 2002/01/09 added */




//EOF
