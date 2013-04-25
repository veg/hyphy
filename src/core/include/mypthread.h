/*  HyPhy - Hypothesis Testing Using Phylogenies.  Copyright (C) 1997-now
 * Core Developers: Sergei L Kosakovsky Pond (spond@ucsd.edu) Art FY Poon
 * (apoon@cfenet.ubc.ca) Steven Weaver (sweaver@ucsd.edu)  Module
 * Developers: Lance Hepler (nlhepler@gmail.com) Martin Smith
 * (martin.audacis@gmail.com)  Significant contributions from: Spencer V Muse
 * (muse@stat.ncsu.edu) Simon DW Frost (sdf22@cam.ac.uk)  Permission is
 * hereby granted, free of charge, to any person obtaining a copy of this
 * software and associated documentation files (the "Software"), to deal in the
 * Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:  The above
 * copyright notice and this permission notice shall be included in all copies
 * or substantial portions of the Software.  THE SOFTWARE IS PROVIDED "AS IS",
 * WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
 * TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  */ //#define
                                                                //_BSD_CLOCK_T_DEFINED_ //#define
                                                                                                          //_BSD_SIZE_T_DEFINED_ //#define
                                                                                                                                               //_BSD_SSIZE_T_DEFINED_ //#define
                                                                                                                                                                                         //_BSD_TIME_T_DEFINED_ //#undef
                                                                                                                                                                                                                                  //__dest_os //#define       __dest_os
                                                                                                                                                                                                                                                                //__mac_os_x #include
    "time.h" #include<
        pthread
            .h> //#undef            __dest_os //#define       __dest_os __mac_os #include
    "CFBundle.h" #include "CFString.h" #include
    "Folders.h" /*#define            _MAC#ifndef __POSIX_LIB__struct
                _pthread_handler_rec{    void           (*routine)(void *);
                void           *arg;                 /    struct
                _pthread_handler_rec *next;};extern "C"{#define
                __PTHREAD_SIZE__           596#define __PTHREAD_ATTR_SIZE__
                36#define __PTHREAD_MUTEXATTR_SIZE__ 8#define
                __PTHREAD_MUTEX_SIZE__     40#define __PTHREAD_CONDATTR_SIZE__
                4#define __PTHREAD_COND_SIZE__      24#define
                __PTHREAD_ONCE_SIZE__      4typedef struct _opaque_pthread_t
                { long sig; struct _pthread_handler_rec  *cleanup_stack; char
                opaque[__PTHREAD_SIZE__];} *pthread_t;typedef struct
                _opaque_pthread_attr_t { long sig; char
                opaque[__PTHREAD_ATTR_SIZE__]; } pthread_attr_t;typedef struct
                _opaque_pthread_mutexattr_t { long sig; char
                opaque[__PTHREAD_MUTEXATTR_SIZE__]; }
                pthread_mutexattr_t;typedef struct _opaque_pthread_mutex_t {
                long sig; char opaque[__PTHREAD_MUTEX_SIZE__]; }
                pthread_mutex_t;typedef struct _opaque_pthread_condattr_t {
                long sig; char opaque[__PTHREAD_CONDATTR_SIZE__]; }
                pthread_condattr_t;typedef struct _opaque_pthread_cond_t {
                long sig;  char opaque[__PTHREAD_COND_SIZE__]; }
                pthread_cond_t;typedef struct { long sig; char
                opaque[__PTHREAD_ONCE_SIZE__]; }
                pthread_once_t;#endiftypedef unsigned long
                pthread_key_t;}*/ /*extern "C"{typedef int
                                   (*pthread_join_ptr) (pthread_t thread, void
                                   **value_ptr);typedef int
                                   (*pthread_create_ptr)(pthread_t *thread,
                                   const pthread_attr_t *attr,
                                   void *(*start_routine)(void *),
                                   void *arg);typedef int
                                   (*pthread_mutex_lock_ptr)   (pthread_mutex_t
                                   *mutex);typedef int
                                   (*pthread_mutex_unlock_ptr) (pthread_mutex_t
                                   *mutex);extern  pthread_join_ptr
                                   pthread_join_p;extern  pthread_create_ptr
                                   pthread_create_p;extern
                                   pthread_mutex_lock_ptr
                                   pthread_mutex_lock_p;extern
                                   pthread_mutex_unlock_ptr
                                   pthread_mutex_unlock_p;}*/ typedef void *(
                                       *PThreadHook)(void *data);
void *MachOFunctionPointerForCFMFunctionPointer(void *cfmfp);
void *ThreadReleafFunctionAA(void *arg);
void *ThreadReleafFunctionCodon(void *arg);
void *ThreadReleafFunctionNuc(void *arg);
void *ThreadReleafFunctionMNuc(void *arg);
void *ThreadReleafFunctionMAA(void *arg);
void *ThreadReleafFunctionMCodon(void *arg);
void *MatrixUpdateFunction(void *arg);
extern PThreadHook ThreadReleafFunctionAAHook, ThreadReleafFunctionCodonHook,
    ThreadReleafFunctionNucHook, ThreadReleafFunctionMNucHook,
    ThreadReleafFunctionMAAHook, ThreadReleafFunctionMCodonHook,
    MatrixUpdateFunctionHook, StateCounterMPHook;
OSStatus LoadFrameworkBundle(CFStringRef framework, CFBundleRef *bundlePtr);
void ImportBundle(void);