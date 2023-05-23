/*

HyPhy - Hypothesis Testing Using Phylogenies.

This file implements data and interface classes
used for linking HyPhy with MEGA

Written by SL Kosakovsky Pond; June 2007
Dedicated to Comet (http://www.hyphy.org/comet.jpg)
?/?/1999-06/05/2007

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


#include "global_object_lists.h"
#include "global_things.h"
#include "THyPhy.h"
#include "batchlan.h"
#include "string.h"

using namespace hy_global;
using namespace hyphy_global_objects;

/* declare some utility functions */


/* global variables */

_String _tHYPHYAskFor           ("_THyPhyAskFor"),
        _tHYPHYNotHandled       ("_THyPhy_NOT_HANDLED_"),
        _tHYPHYCurrentStatus;

long    _tHYPHYDone             = 0;
double  _tHYPHYValue            = 0.0;

_THyPhy * globalInterfaceInstance = nil;

//_________________________________________________________
// default callback hanlder
bool _tHyPhyDefaultHandler (const char*,int,double) {
    return true;
}

//_________________________________________________________
// default callback hanlder
//void SetGlobalInterfaceInstance (const _THyPhy* hi)
//{
//  globalInterfaceInstance = hi;
//}


/* Basic return types supported by the interface object */

//_________________________________________________________
_THyPhyString* _THyPhyReturnObject::castToString(void) {
    if (myType() == THYPHY_TYPE_STRING) {
        return (_THyPhyString*)this;
    }
    return NULL;
}

//_________________________________________________________
_THyPhyNumber* _THyPhyReturnObject::castToNumber(void)
{
    if (myType() == THYPHY_TYPE_NUMBER) {
        return (_THyPhyNumber*)this;
    }
    return NULL;
}

//_________________________________________________________
_THyPhyMatrix* _THyPhyReturnObject::castToMatrix(void)
{
    if (myType() == THYPHY_TYPE_MATRIX) {
        return (_THyPhyMatrix*)this;
    }
    return NULL;
}


//_________________________________________________________
_THyPhyString::_THyPhyString(const char* characters, long length) {
    if (characters) {
        if (length == 0) {
            while (characters[length++]) ;
            length --;
        }
        sData = (char*)MemAllocate (length+1);
        memcpy       (sData,characters,length+1);
    } else {
        sData   = nil;
    }
    sLength = length;
}

//_________________________________________________________
_THyPhyString::~_THyPhyString(void) {
    if (sData) {
        free (sData);
    }
};

//_________________________________________________________
_THyPhyNumber::_THyPhyNumber(double v) {
    nValue = v;
}


//_________________________________________________________
_THyPhyMatrix::_THyPhyMatrix(void) {
    mData = nil;
    mRows = 0;
    mCols = 0;
}

//_________________________________________________________
_THyPhyMatrix::_THyPhyMatrix(const long r, const long c, const double* d) {
    mData = (double*)MemAllocate (r*c*sizeof(double));
    mRows = r;
    mCols = c;
    for (long i = 0; i < r*c; i++) {
        mData[i] = d[i];
    }

}

//_________________________________________________________
_THyPhyMatrix::~_THyPhyMatrix(void) {
    if (mData) {
        free (mData);
    }
};

//_________________________________________________________
double _THyPhyMatrix::MatrixCell(long r, long c) {
    return mData[mCols*r+c];
}

//_________________________________________________________
// Begin _THyPhy definitions
//_________________________________________________________

_THyPhy::_THyPhy(_ProgressCancelHandler* mHandler, const char* baseDirPath, long cpuCount) {
    InitTHyPhy (mHandler, baseDirPath, cpuCount);
}

//_________________________________________________________


_THyPhy::_THyPhy(const char* baseDirPath, long cpuCount) {
    InitTHyPhy (_tHyPhyDefaultHandler, baseDirPath, cpuCount);
}

//_________________________________________________________

_THyPhy::~_THyPhy           (void) {
    if (currentResultHolder) {
        delete (_THyPhyReturnObject*)currentResultHolder;
    }
    if (baseDirectoryInstance) {
        delete (baseDirectoryInstance);
    }
    ClearAll();
    DeleteObject        ((_StringBuffer*)errors);
    DeleteObject        ((_StringBuffer*)warnings);
    DeleteObject        ((_StringBuffer*)textout);
    if (globalInterfaceInstance == this) {
        globalInterfaceInstance = nil;
    }

    //PurgeAll(true);
    GlobalShutdown();
}

//_________________________________________________________

void _THyPhy::InitTHyPhy (_ProgressCancelHandler* mHandler, const char* baseDirPath, long cpuCount)
{
    char dirSlash = get_platform_directory_char ();
    SetCallbackHandler (mHandler);
    currentResultHolder = new _THyPhyString;
    askFID = -1;
    system_CPU_count = MAX (cpuCount, 1);
    if (baseDirPath) {
        // set base directory
        hy_base_directory = baseDirPath;
        if (hy_base_directory.get_char(hy_base_directory.length()-1) != dirSlash) {
            hy_base_directory = hy_base_directory & dirSlash;
        }
        baseDirectoryInstance = new _THyPhyString (hy_base_directory.get_str());
        hy_base_directory = baseDirectoryInstance->sData;
        pathNames && &hy_base_directory;
    }

#ifdef _HYPHY_LIBDIRECTORY_    
    hy_lib_directory = _HYPHY_LIBDIRECTORY_;
    if (hy_lib_directory.get_char(hy_lib_directory.length()-1) != dirSlash) {
        hy_lib_directory = hy_lib_directory & dirSlash;
    }
#else
    if (hy_base_directory)
        hy_lib_directory = hy_base_directory;
#endif

    pathNames && &hy_lib_directory;
    GlobalStartup();
    errors   = nil;
    warnings = nil;
    textout  = nil;
    globalInterfaceInstance = this;

}

//_________________________________________________________

_THyPhyString * _THyPhy::ExecuteBF (const char * buffer, bool doPurge) {
        
    if (doPurge) {
        PurgeAll            (true);    // cleanup results of previous analysis
    }

    _String             dd (get_platform_directory_char());

    _FString            bp  (hy_base_directory, false),
                        lp  (hy_lib_directory, false),
                        ds  (dd),
                        cfp (pathNames.lLength?*(_String*)pathNames(pathNames.countitems()-1):kEmptyString),
                        * stashed = (_FString*)FetchObjectFromVariableByType (&hy_env::path_to_current_bf, STRING);

    setParameter        (hy_env::directory_separator_char, &ds);
    setParameter        (hy_env::base_directory, &bp);
    setParameter        (hy_env::lib_directory, &lp);

    if (stashed) {
        stashed = (_FString*)stashed->makeDynamic();
    }
    setParameter        (hy_env::path_to_current_bf,&cfp);

    _StringBuffer             commandString (buffer);

    if (commandString.BeginsWith ("#NEXUS"),false) {
        lastNexusDataMatrix = ReadDataSetFile (nil, 2, &commandString);
        commandString = nexusBFBody;
    }

    _ExecutionList      compiledCode  (commandString);

    BatchDelete        ((_StringBuffer*)errors, (_StringBuffer*)warnings, (_StringBuffer*)textout);
    
    errors              = new _StringBuffer (128);
    warnings            = new _StringBuffer (128L);
    textout             = new _StringBuffer (128L);

    compiledCode.Execute ();
    askFID = FindBFFunctionName (_tHYPHYAskFor, NULL);
    
    //printf ("\n_THyPhy::ExecuteBF %d\n", askFID);
    
    HBLObjectRef bfReturn  = compiledCode.GetResult ();

    ((_StringBuffer*)errors)->TrimSpace();
    ((_StringBuffer*)warnings)->TrimSpace();
    ((_StringBuffer*)textout)->TrimSpace();

    if (currentResultHolder->sData) {
        free (currentResultHolder->sData);
        currentResultHolder->sData = nil;
    }
    if (bfReturn) {
        _String * serializedReturn   = (_String*) bfReturn->toStr();
        currentResultHolder->sData   = (char*)serializedReturn->get_str();
        //serializedReturn->sData      = nil;
        currentResultHolder->sLength = serializedReturn->length();
        serializedReturn->Initialize();
    }
    return currentResultHolder;
}

//_________________________________________________________

void _THyPhy::ClearAll (void) {
    PurgeAll(true);
}

//_________________________________________________________

void* _THyPhy::AskFor (const char* resultID) {
    if (resultID && askFID >= 0) {
        _StringBuffer theCommand (128L);
        theCommand << "return "
         <<  _tHYPHYAskFor
         << "(\"";
        theCommand.SanitizeAndAppend (resultID) << "\");";
        //printf ("\n_THyPhy::AskFor %s\n", theCommand.get_str());
        _ExecutionList      compiledCode  (theCommand);
        compiledCode.Execute();
        HBLObjectRef retResult = compiledCode.GetResult ();
        if (retResult && retResult->ObjectClass() == STRING) {
            _FString * checkHandled = (_FString*)retResult;
            if (checkHandled->get_str() == _tHYPHYNotHandled) {
                return nil;
            }
        }
        return retResult->makeDynamic();
    }
    return nil;
}

//_________________________________________________________

void _THyPhy::DumpResult (void* aResult) {
    if (aResult) {
        DeleteObject ((BaseRef)aResult);
    }
}

//_________________________________________________________

void _THyPhy::SetCallbackHandler (_ProgressCancelHandler* newHandler) {
    theHandler = newHandler;
}

//_________________________________________________________

_ProgressCancelHandler* _THyPhy::GetCallbackHandler (void) {
    return theHandler;
}

//_________________________________________________________

void        _THyPhy::PushWarning (const void * o) {
    if (warnings) {
        *((_StringBuffer*)warnings) << *(const _String*)o;
    }
}

//_________________________________________________________

void        _THyPhy::PushError (const void * o) {
    if (errors) {
        *((_StringBuffer*)errors) << *(const _String*)o;
    }
}

//_________________________________________________________

void        _THyPhy::PushOutString (const void * o) {
    if (textout) {
        *((_StringBuffer*)textout) << *(const _String*)o;
    }
}

//_________________________________________________________

_THyPhyString   * _THyPhy::ConvertHyPhyString (void * o) {
    return new _THyPhyString (((_String*)o)->get_str(),((_String*)o)->length());
}

//_________________________________________________________

bool _THyPhy::CanCast (const void* theObject, const int requestedType) {
    if (theObject) {
        switch (((HBLObjectRef)theObject)->ObjectClass()) {
        case NUMBER:
            return requestedType!=THYPHY_TYPE_JSON;;
            // can cast a number to everything
        case STRING:
            return requestedType==THYPHY_TYPE_NUMBER || requestedType==THYPHY_TYPE_STRING;
            // can cast anything to a string
                
        case MATRIX:
            return requestedType!=THYPHY_TYPE_NUMBER;
            // can not cast matrix to number

        case TREE:
        case TOPOLOGY:
            return requestedType==THYPHY_TYPE_STRING;

        case ASSOCIATIVE_LIST:
            return requestedType==THYPHY_TYPE_JSON;

        }
    }
    return false;
}

//_________________________________________________________

_THyPhyReturnObject* _THyPhy::CastResult (const void* theObject, const int requestedType) {
    static const _String kUseJSONForMatrix ("USE_JSON_FOR_MATRIX");
    
    _THyPhyReturnObject * convertedObject = nil;
    if (CanCast(theObject,requestedType)) {
        int hyphyObjClass = ((HBLObjectRef)theObject)->ObjectClass();
        switch (requestedType) {
            case THYPHY_TYPE_NUMBER: {
                if (hyphyObjClass == NUMBER) {
                    return new _THyPhyNumber (((HBLObjectRef)theObject)->Compute()->Value());
                }
                if (hyphyObjClass == STRING) {
                    _String sV ((_String*)((_FString*)theObject)->toStr());
                    return new _THyPhyNumber (sV.to_float());
                }
            }
            case THYPHY_TYPE_STRING: {
                _String sV ((_String*)((HBLObjectRef)theObject)->toStr());
                return new _THyPhyString (sV.get_str(),sV.length());
            }
            case THYPHY_TYPE_MATRIX: {
                if (hyphyObjClass == NUMBER) {
                    double evaluate = ((HBLObjectRef)theObject)->Compute()->Value();
                    return new _THyPhyMatrix (1,1,&evaluate);
                }

                if (hyphyObjClass == MATRIX) {
                    _Matrix * evalutedNumeric =  (_Matrix*)((_Matrix*)(((HBLObjectRef)theObject)->Compute()))
                                                 ->ComputeNumeric();

                    return new _THyPhyMatrix (evalutedNumeric->GetHDim(),evalutedNumeric->GetVDim(),evalutedNumeric->theData);
                }
            }
            case THYPHY_TYPE_JSON: {
                HBLObjectRef stash = hy_env :: EnvVariableGet (kUseJSONForMatrix, HY_BL_ANY);
                if (stash) {
                    stash->AddAReference();
                } else {
                    stash = new _MathObject;
                }
                hy_env :: EnvVariableSet (kUseJSONForMatrix, new _Constant (1.),false);
                if (hyphyObjClass == ASSOCIATIVE_LIST) {
                    _String sV ((_String*)((_AssociativeList*)(((HBLObjectRef)theObject)))->toStr());
                    hy_env :: EnvVariableSet (kUseJSONForMatrix, stash , false);
                    return new _THyPhyString (sV.get_str(),sV.length());
                }

                if (hyphyObjClass == MATRIX) {
                    _String sV ((_String*)((_Matrix*)(((HBLObjectRef)theObject)))->toStr());
                    hy_env :: EnvVariableSet (kUseJSONForMatrix, stash , false);
                    return new _THyPhyString (sV.get_str(),sV.length());
                }
                hy_env :: EnvVariableSet (kUseJSONForMatrix, stash , false);

            }
        }
    }

    return convertedObject;
}

//_________________________________________________________________________

void    SetStatusLine (_String arg) {
    _tHYPHYCurrentStatus = arg;
}

//_________________________________________________________________________

void    SetStatusBarValue (long l, hyFloat max, hyFloat rate) {
    _tHYPHYDone = l;
    _tHYPHYCurrentStatus   = _String ("LF Optimization. Value=") & _String (max) &", "&_String (rate) & " evals/sec.";
    _tHYPHYValue = max;
}

//_________________________________________________________________________

long    _THyPhyGetLongStatus        (void) {
    return _tHYPHYDone;
}

//_________________________________________________________________________

double  _THyPhyGetDoubleStatus      (void) {
    return _tHYPHYValue;
}

//_________________________________________________________________________

const char*   _THyPhyGetStringStatus      (void) {
    return _tHYPHYCurrentStatus.get_str();
}
