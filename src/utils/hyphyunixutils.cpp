/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (spond@temple.edu)
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

#include <stdio.h>


#ifdef  __HYPHYCURL__
    #include <curl/curl.h>
#endif

#include "global_object_lists.h"
#include "likefunc.h"
#include "hy_string_buffer.h"
#include "hbl_env.h"

using    namespace hy_global;


#ifdef          __HYPHYMPI__
_String         preserveSlaveNodeState ("PRESERVE_SLAVE_NODE_STATE"),
                MPI_NEXUS_FILE_RETURN  ("MPI_NEXUS_FILE_RETURN");


void            mpiNormalLoop    (int, int, _String &);
void            mpiOptimizerLoop (int, int);

void            mpiBgmLoop (int, int);
#endif

extern _List batchLanguageFunctionNames;


//_________________________________________________________________________

size_t url2File   (void *ptr, size_t size, size_t nmemb, void *stream) {
    return fwrite (ptr, size, nmemb, (FILE*)stream);
}

//_________________________________________________________________________

size_t url2String (void *ptr, size_t size, size_t nmemb, void *stream) {
    _StringBuffer * s = (_StringBuffer*)stream;
    char    * p = (char*)ptr;

    for (unsigned long k=0UL; k<size*nmemb; k++) {
        (*s) << p[k];
    }

    return size*nmemb;
}

//_________________________________________________________________________

bool    Get_a_URL (_String& urls, _String*
#ifdef __HYPHYCURL__
                   fileName
#endif
                   ) {
#ifdef __HYPHYCURL__
    CURL *curl;
    CURLcode res ;
    curl = curl_easy_init ();
    FILE   * f = nil;
    _StringBuffer * s = nil;
    char cErr [CURL_ERROR_SIZE+1];
    if(curl) {
        if (fileName) {
            f = fopen (fileName->get_str(),"wb");
            if (!f) {
                urls = _String ("Failed to open ") & *fileName & " for writing";
                return false;
            }
        } else {
            s = new _StringBuffer (8192UL);
        }

        curl_easy_setopt (curl, CURLOPT_URL, urls.get_str() );
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
        curl_easy_setopt (curl, CURLOPT_USERAGENT, ver.get_str());
        //curl_easy_setopt (curl, CURLOPT_VERBOSE, 1);
        curl_easy_setopt (curl, CURLOPT_WRITEFUNCTION, (void*)(f?url2File:url2String));
      
        /*if (hy_env::EnvVariableTrue (VerbosityLevelString)) {
            curl_easy_setopt (curl,CURLOPT_NOPROGRESS,1);
        }*/
        
        res = curl_easy_perform (curl);
        curl_easy_cleanup (curl);

        if (f) {
            fclose (f);
        } else {
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

_String *    StringFromConsole   () {
    fflush(stdout);
    _StringBuffer * returnme = new _StringBuffer (32UL);
#if not defined __HEADLESS__ && not defined _MINGW32_MEGA_
    int       readAChar;
    while    ((readAChar = getc(stdin)) != '\n') {
        if (readAChar == EOF) {
            hy_env::EnvVariableSet(hy_env::end_of_file, new HY_CONSTANT_TRUE, false);
            break;
        }
        *returnme << (char)readAChar;
    }
#else
    HandleApplicationError("Unhandled standard input interaction in StringFromConsole for headless HyPhy");
    return NULL;
#endif
    return returnme;
}

#ifdef _USE_EMSCRIPTEN_
#include "emscripten.h"

EM_JS(void, _jsSendStatusUpdateStdout, (const char *status_update), {
    // Send a message back to main thread from WebWorker
    postMessage({
        type: "biowasm",
        value: {
        text: Module.AsciiToString(status_update),
        type : "print"
        }
    });
})
#endif
//__________________________________________________________________________________

void    StringToConsole (_String const & s,  void * extra) {
    BufferToConsole ((const char*)s, extra);
}


//__________________________________________________________________________________

void    BufferToConsole (const char* s, void * extra) {
#ifdef __HYPHYMPI__
    if (hy_mpi_node_rank == 0)
#endif
#ifdef __HEADLESS__
        if (globalInterfaceInstance) {
            _String st (s);
            globalInterfaceInstance->PushOutString(&st);
        }
#else
  if (extra) {
    fprintf ((FILE*)extra, "%s",s);
  } else {
#ifdef _USE_EMSCRIPTEN_
    if (!extra) {
        _jsSendStatusUpdateStdout (s);
    }
#endif
    printf ("%s",s);
    fflush(stdout);
  }
#endif
}

//__________________________________________________________________________________

void    NLToConsole (void * extra) {
    BufferToConsole ("\n", extra);
}

//__________________________________________________________________________________

void    ObjectToConsole (BaseRef obj,  void * extra) {
    StringToConsole (_String (((_String*)obj->toStr())), extra );
}

#endif

#ifdef __HYPHYMPI__

//__________________________________________________________________________________
void mpiNormalLoop    (int  rank, int size, _String & baseDir)
{
    long         senderID = 0;

    ReportWarning ("[MPI] Entered mpiNormalLoop");

    _String * theMessage     = MPIRecvString (-1,senderID);   // listen for messages from any node
    _StringBuffer * resStr        = nil;
    
    //int loop_count = 0;

    while (theMessage->nonempty()) {
        hy_env :: EnvVariableSet (hy_env ::mpi_node_id, new _Constant (rank), false);
        hy_env :: EnvVariableSet (hy_env ::mpi_node_count, new _Constant (size), false);
 
        DeleteObject (resStr);
        resStr       = nil;
        if (theMessage->BeginsWith (mpiLoopSwitchToOptimize) ) {
            hyphyMPIOptimizerMode   = theMessage->Cut(mpiLoopSwitchToOptimize.length(),kStringEnd).to_long();

            ReportWarning           (_String("[MPI] Switched to mpiOptimizer loop with mode ") & hyphyMPIOptimizerMode);
            MPISendString           (mpiLoopSwitchToOptimize,senderID);
            mpiOptimizerLoop        (rank,size);
            ReportWarning           ("[MPI] Returned from mpiOptimizer loop");
            hyphyMPIOptimizerMode   = _hyphyLFMPIModeNone;
            PushFilePath(baseDir, false, false);
        } else if ( *theMessage == mpiLoopSwitchToBGM) {
            ReportWarning       ("[MPI] Received signal to switch to mpiBgmLoop");
            MPISendString       (mpiLoopSwitchToBGM, senderID); // feedback to source to confirm receipt of message
            mpiBgmLoop          (rank, size);
            ReportWarning       ("[MPI] Returned from mpiBgmLoop");
        } else {
            if (theMessage->BeginsWith ("#NEXUS")) {
                _String             msgCopy (*theMessage);
                //printf ("\n\nMPI Node %d; message \n %s\n", rank, theMessage->get_str());
                ReportWarning       ("[MPI] Received a likelihood function");
                //ReportWarning       (msgCopy);
                ReadDataSetFile     (nil,true,theMessage);
                ReportWarning       ("[MPI] Read/optimized the likelihood function");
                _Variable*          lfName = FetchVar(LocateVarByName(MPI_NEXUS_FILE_RETURN));

                if (lfName) {
                    resStr = new _StringBuffer ((_String*)lfName->Compute()->toStr());
                 } else {
                    _FString        *lfID = (_FString*)FetchObjectFromVariableByType (&lf2SendBack, STRING);

                    if (!lfID) {
                        HandleApplicationError (_String("[MPI] Malformed MPI likelihood function optimization request - did not specify the LF name to return in variable ") & lf2SendBack & ".\n\n\n" );
                        break;
                    }

                    long type = HY_BL_LIKELIHOOD_FUNCTION, index;
                    
                    _LikelihoodFunction *lf = (_LikelihoodFunction *)hyphy_global_objects::_HYRetrieveBLObjectByName    (lfID->get_str(), type, &index, false, false);

                    if (lf == nil) {
                        HandleApplicationError (_String("[MPI] Malformed MPI likelihood function optimization request - '") & lfID->get_str() &"' did not refer to a well-defined likelihood function.\n\n\n");
                        break;
                    }
                  
                    resStr       = new _StringBuffer (1024UL);
                    lf->SerializeLF(*resStr,hy_env::EnvVariableTrue (hy_env::short_mpi_return) ? _hyphyLFSerializeModeShortMPI:_hyphyLFSerializeModeLongMPI);
                }
            } else {
                //ReportWarning(_String ("[MPI] Received commands\n") & *theMessage & "\n");
                //printf ("\n\nMPI Node %d; message \n %s\n", rank, theMessage->get_str());
                //printf ("\n\nMPI Node %d; likefuncs \n %s\n", rank, ((_String*)likeFuncNamesList.toStr())->get_str());
                _StringBuffer code (*theMessage);
		_ExecutionList exL (code);
                
                 
                //ReportWarning (_String ((_String*)batchLanguageFunctionNames.toStr()));
                HBLObjectRef res = exL.Execute();
                //printf ("\n==>%d\n", hy_env::EnvVariableTrue (preserveSlaveNodeState));
                if (res) {
                  resStr = new _StringBuffer ((_String*)res->toStr());
                } else {
                  resStr = new _StringBuffer ("0");
                }
            }

            MPISendString(*resStr,senderID);
            
            /*loop_count++;
            printf ("\n\nMPI Node %d; variable array length %d\n", rank, variableNames.countitems());
            if (loop_count == 20) {
                printf ("%s\n", ((_String*)variableNames.toStr())->get_str());
                abort();
            }*/
            
            
            if (hy_env::EnvVariableTrue (preserveSlaveNodeState) == false) {
                printf ("\n\nMPI Node %d; PURGING\n", rank);
                PurgeAll (true);
                InitializeGlobals ();
                PushFilePath(baseDir, false, false);
                ReportWarning("Reset node state");
            } else {
                ReportWarning("Preserved node state");
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
    while (theMessage->nonempty()) {
        if (theMessage->BeginsWith ("#NEXUS")) {
            //ReportWarning (*theMessage);
            ReadDataSetFile (nil,true,theMessage);
            if (likeFuncNamesList.lLength!=1) {
                HandleApplicationError ("[MPI] Malformed MPI likelihood function paraller optimizer startup command. Exactly ONE valid LF must be defined.n\n\n");
                break;
            }

            // send back the list of independent variables

            _LikelihoodFunction * theLF = (_LikelihoodFunction*)likeFuncList (0);
            if (hyphyMPIOptimizerMode == _hyphyLFMPIModeREL && theLF->CountObjects (kLFCountCategoryVariables)) {
                HandleApplicationError (_String("[MPI] Likelihood functions spawned off to slave MPI nodes can't have category variables.n\n\n"));
                break;
            }

            _SimpleList const * ivl = & theLF->GetIndependentVars();
  
          
            _StringBuffer      variableSpec (128UL);


            (variableSpec) << LocateVar(ivl->list_data[0])->GetName();
            
            for (long kk = 1; kk < ivl->lLength; kk++) {
              (variableSpec) << ';';
              (variableSpec) << LocateVar(ivl->list_data[kk])->GetName();
            }
          
            ReportWarning         (_String("[MPI] Sending back the following variable list\n") & variableSpec);
            MPISendString         (variableSpec,senderID);
            theLF->PrepareToCompute();
            theLF->MPI_LF_Compute (senderID, !(hyphyMPIOptimizerMode == _hyphyLFMPIModeREL ||
                                               hyphyMPIOptimizerMode  == _hyphyLFMPIModeSiteTemplate));
            theLF->DoneComputing();
            PurgeAll (true);
            InitializeGlobals ();
            ReportWarning("Reset node state at the end of MPI optimizaer loop");
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

    while (theMessage->nonempty()) {

	_StringBuffer code (*theMessage);
        _ExecutionList  exL (code);
        HBLObjectRef       res = exL.Execute();    // should send this process into CacheNodeScores()

        resStr = res ? (_String*)res->toStr() : new _String ("0");
        ReportWarning (_String ("MPI Node: ") & (long)rank & " executed HBL with result:\n" & resStr);

        if (bgmNamesList.lLength < 1) {
            HandleApplicationError ("Malformed HBL. No valid BGM has been defined.\n");
            break;
        }
    }

    DeleteObject (theMessage);
}
#endif
