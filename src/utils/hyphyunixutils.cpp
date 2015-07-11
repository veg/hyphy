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

#include <stdio.h>
#include "likefunc.h"



#ifndef __HYPHY_NO_CURL__
//#define __HYPHYCURL__
#endif

#ifdef  __HYPHYCURL__
#include <curl/curl.h>
#endif

#ifdef          __HYPHYMPI__
_String         preserveSlaveNodeState ("PRESERVE_SLAVE_NODE_STATE"),
                MPI_NEXUS_FILE_RETURN  ("MPI_NEXUS_FILE_RETURN");

int             _hy_mpi_node_rank;

void            mpiNormalLoop    (int, int, _String &);
void            mpiOptimizerLoop (int, int);

void            mpiBgmLoop (int, int);
_SimpleList     mpiNodesThatCantSwitch;
#endif


//_________________________________________________________________________

size_t url2File   (void *ptr, size_t size, size_t nmemb, void *stream)
{
    return fwrite (ptr, size, nmemb, (FILE*)stream);
}

//_________________________________________________________________________

size_t url2String (void *ptr, size_t size, size_t nmemb, void *stream)
{
    _String * s = (_String*)stream;
    char    * p = (char*)ptr;

    for (unsigned long k=0; k<size*nmemb; k++) {
        (*s) << p[k];
    }

    return size*nmemb;
}

//_________________________________________________________________________

bool    Get_a_URL (_String& urls, _String* fileName)
{
#ifdef __HYPHYCURL__
    CURL *curl;
    CURLcode res ;
    curl = curl_easy_init ();
    FILE   * f = nil;
    _String* s = nil;
    char cErr [CURL_ERROR_SIZE+1];
    if(curl) {
        if (fileName) {
            f = fopen (fileName->sData,"wb");
            if (!f) {
                urls = _String ("Failed to open ") & *fileName & " for writing";
                return false;
            }
        } else {
            s = new _String (8192, true);
            checkPointer (s);
        }

        curl_easy_setopt (curl, CURLOPT_URL, urls.sData );
        curl_easy_setopt (curl, CURLOPT_ERRORBUFFER, cErr);

        //Do not check peer certificate, since we only ever get urls
        curl_easy_setopt (curl, CURLOPT_SSL_VERIFYPEER, 0);
        curl_easy_setopt (curl, CURLOPT_SSL_VERIFYHOST, 0);
        
        if (f) {
            curl_easy_setopt (curl, CURLOPT_FILE, (void*)f);
        } else {
            curl_easy_setopt (curl, CURLOPT_FILE, (void*)s);
        }

        _String ver (GetVersionString());
        curl_easy_setopt (curl, CURLOPT_USERAGENT, ver.sData);
        //curl_easy_setopt (curl, CURLOPT_VERBOSE, 1);
        curl_easy_setopt (curl, CURLOPT_WRITEFUNCTION, (void*)(f?url2File:url2String));
        _Parameter vbl = 0.0;
        checkParameter (VerbosityLevelString,vbl,0.0);
        if (vbl<0.5) {
            curl_easy_setopt (curl,CURLOPT_NOPROGRESS,1);
        }
        res = curl_easy_perform (curl);
        curl_easy_cleanup (curl);

        if (f) {
            fclose (f);
        } else {
            s->Finalize();
            urls = *s;
            DeleteObject (s);
        }
        if (!res) {
            return true;
        }
    } else {
        urls = "Failed to initialize CURL object";
        return false;
    }
    urls = _String ("CURL error:") & (long)res & "." & cErr;
    return false;
#else
    urls = "This feature requires libcurl";
    return false;
#endif
}

#ifndef __HYPHY_GTK__

//____________________________________________________________________________________

_String*    StringFromConsole   (bool)
{
    _String * returnme = new _String (32L, true);
#if not defined __HEADLESS__ && not defined _MINGW32_MEGA_
    int       readAChar;
    while    ((readAChar = getc(stdin)) != '\n') {
        if (readAChar == EOF) {
            CheckReceptacleAndStore (&hasEndBeenReached,empty,false,new _Constant (1.), false);
            break;
        }
        *returnme << readAChar;
    }
#else
    WarnError ("Unhandled standard input interaction in StringFromConsole for headless HyPhy");
    return NULL;
#endif
    returnme->Finalize ();
    return returnme;
}

//__________________________________________________________________________________

void    StringToConsole (_String & s,  _SimpleList *)
{
    BufferToConsole ((const char*)s.sData);
}


//__________________________________________________________________________________

void    BufferToConsole (const char* s, _SimpleList *)
{
#ifdef __HYPHYMPI__
    if (_hy_mpi_node_rank == 0)
#endif
#ifdef __HEADLESS__
        if (globalInterfaceInstance) {
            _String st (s);
            globalInterfaceInstance->PushOutString(&st);
        }
#else
        printf ("%s",s);
#endif
}

//__________________________________________________________________________________

void    NLToConsole (void)
{
    BufferToConsole ("\n");
}

#endif

#ifdef __HYPHYMPI__

//__________________________________________________________________________________
void mpiNormalLoop    (int rank, int size, _String & baseDir)
{
    long         senderID = 0;

    ReportWarning ("[MPI] Entered mpiNormalLoop");

    _String* theMessage     = MPIRecvString (-1,senderID),  // listen for messages from any node
             * resStr        = nil;

    while (theMessage->sLength) {
        setParameter    (mpiNodeID,    (_Parameter)rank);
        setParameter    (mpiNodeCount, (_Parameter)size);

        //ReportWarning (*theMessage);
        DeleteObject (resStr);
        resStr       = nil;
        if (theMessage->startswith (mpiLoopSwitchToOptimize) ) {
            hyphyMPIOptimizerMode   = theMessage->Cut(mpiLoopSwitchToOptimize.sLength,-1).toNum();

            ReportWarning           (_String("[MPI] Switched to mpiOptimizer loop with mode ") & hyphyMPIOptimizerMode);
            MPISendString           (mpiLoopSwitchToOptimize,senderID);
            mpiOptimizerLoop        (rank,size);
            ReportWarning           ("[MPI] Returned from mpiOptimizer loop");
            hyphyMPIOptimizerMode   = _hyphyLFMPIModeNone;
            pathNames               && & baseDir;
        } else if ( theMessage->Equal (&mpiLoopSwitchToBGM) ) {
            ReportWarning       ("[MPI] Received signal to switch to mpiBgmLoop");
            MPISendString       (mpiLoopSwitchToBGM, senderID); // feedback to source to confirm receipt of message
            mpiBgmLoop          (rank, size);
            ReportWarning       ("[MPI] Returned from mpiBgmLoop");
        } else {
            if (theMessage->beginswith ("#NEXUS")) {
                _String             msgCopy (*theMessage);
                ReportWarning       ("[MPI] Received a function to optimize");
                ReadDataSetFile     (nil,true,theMessage);
                ReportWarning       ("[MPI] Done with the optimization");
                _Variable*          lfName = FetchVar(LocateVarByName(MPI_NEXUS_FILE_RETURN));

                if (lfName) {
                    resStr = (_String*)(lfName->Compute()->toStr());
                } else {
                    _FString        *lfID = (_FString*)FetchObjectFromVariableByType (&lf2SendBack, STRING);

                    if (!lfID) {
                        FlagError (_String("[MPI] Malformed MPI likelihood function optimization request - did not specify the LF name to return in variable ") & lf2SendBack & ".\n\n\n" );
                        break;
                    }

                    long f = likeFuncNamesList.Find (lfID->theString);

                    if (f<0) {
                        FlagError ("[MPI] Malformed MPI likelihood function optimization request - LF name to return did not refer to a well-defined likelihood function.\n\n\n");
                        break;
                    }
                    _Parameter      pv;
                    checkParameter (shortMPIReturn, pv ,0);
                    resStr       = (_String*)checkPointer(new _String (1024L,true));
                    ((_LikelihoodFunction*)likeFuncList (f))->SerializeLF(*resStr,pv>0.5?_hyphyLFSerializeModeShortMPI:_hyphyLFSerializeModeLongMPI);
                    resStr->Finalize();
                }
            } else {
                _ExecutionList exL (*theMessage);
                _PMathObj res = exL.Execute();
                resStr = res?(_String*)res->toStr():new _String ("0");
            }

            checkPointer (resStr);
            MPISendString(*resStr,senderID);

            _Parameter      keepState = 0.0;
            checkParameter  (preserveSlaveNodeState, keepState, 0.0);

            if (keepState < 0.5) {
                PurgeAll (true);
                pathNames && & baseDir;
            }
        }
        DeleteObject (theMessage);
        theMessage = MPIRecvString (-1,senderID);
    }
    /*MPISendString(empty,senderID);*/
    DeleteObject (resStr);
    DeleteObject (theMessage);
}

//__________________________________________________________________________________
void mpiOptimizerLoop (int rank, int size)
{
    long         senderID = 0;

    ReportWarning (_String ("[MPI] Node:") & (long)rank & " is ready for MPIParallelOptimizer tasks");

    if (hyphyMPIOptimizerMode == _hyphyLFMPIModePartitions) {
        ReportWarning ("[MPI] MPI Partitions mode");
    }

    //printf ("Node %d waiting for a string\n", rank);
    _String* theMessage = MPIRecvString (-1,senderID);
    while (theMessage->sLength) {
        if (theMessage->beginswith ("#NEXUS")) {
            //ReportWarning (*theMessage);
            ReadDataSetFile (nil,true,theMessage);
            if (likeFuncNamesList.lLength!=1) {
                FlagError ("[MPI] Malformed MPI likelihood function paraller optimizer startup command. Exactly ONE valid LF must be defined.n\n\n");
                break;
            }

            // send back the list of independent variables

            _LikelihoodFunction * theLF = (_LikelihoodFunction*)likeFuncList (0);
            if (hyphyMPIOptimizerMode == _hyphyLFMPIModeREL && theLF->CountObjects (4)) {
                FlagError (_String("[MPI] Likelihood functions spawned off to slave MPI nodes can't have category variables.n\n\n"));
                break;
            }

            _SimpleList* ivl = & theLF->GetIndependentVars();
            _String      variableSpec (128L, true);

            (variableSpec) << LocateVar(ivl->lData[0])->GetName();

            for (long kk = 1; kk < ivl->lLength; kk++) {
                (variableSpec) << ';';
                (variableSpec) << LocateVar(ivl->lData[kk])->GetName();
            }
            variableSpec.Finalize();
            ReportWarning         (_String("[MPI] Sending back the following variable list\n") & variableSpec);
            MPISendString         (variableSpec,senderID);
            theLF->PrepareToCompute();
            theLF->MPI_LF_Compute (senderID, !(hyphyMPIOptimizerMode == _hyphyLFMPIModeREL ||
                                               hyphyMPIOptimizerMode  == _hyphyLFMPIModeSiteTemplate));
            theLF->DoneComputing();
            PurgeAll (true);
        }
        DeleteObject (theMessage);
        theMessage = MPIRecvString (-1,senderID);
    }
    DeleteObject (theMessage);
}


//__________________________________________________________________________________
void mpiBgmLoop (int rank, int size)
{
    long        senderID    = 0;
    _String *   resStr      = nil;

    ReportWarning (_String ("MPI Node:") & (long)rank & " is ready for MPIBgmCacheNodeScores tasks");

    // receive serialized Bgm
    _String* theMessage = MPIRecvString (-1, senderID);

    while (theMessage->sLength) {
        _ExecutionList  exL (*theMessage);
        _PMathObj       res = exL.Execute();    // should send this process into CacheNodeScores()

        resStr = res ? (_String*)res->toStr() : new _String ("0");
        ReportWarning (_String ("MPI Node: ") & (long)rank & " executed HBL with result:\n" & resStr);

        if (bgmNamesList.lLength < 1) {
            _String errMsg ("Malformed HBL. No valid BGM has been defined.\n");
            FlagError (errMsg);
            break;
        }
    }

    DeleteObject (theMessage);
}
#endif

