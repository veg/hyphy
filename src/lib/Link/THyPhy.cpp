/*

HyPhy - Hypothesis Testing Using Phylogenies.

This file implements data and interface classes
used for linking HyPhy with MEGA

Written by SL Kosakovsky Pond; June 2007
Dedicated to Comet (http://www.hyphy.org/comet.jpg)
?/?/1999-06/05/2007

Copyright (C) 1997-2007
Primary Development:
  Sergei L Kosakovsky Pond (sergeilkp@mac.com)
Significant contributions from:
  Spencer V Muse (muse@stat.ncsu.edu)
  Simon DW Frost (sdfrost@ucsd.edu)
  Art FY Poon    (apoon@biomail.ucsd.edu)

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


#include "THyPhy.h"
#include "batchlan.h"
#include "string.h"


/* declare some utility functions */

void    ReadPreferences         (void);
void    ApplyPreferences        (void);

/* global variables */

_String _tHYPHYAskFor           ("_THyPhyAskFor"),
        _tHYPHYNotHandled       ("_THyPhy_NOT_HANDLED_"),
        _tHYPHYCurrentStatus;

long    _tHYPHYDone             = 0;
double  _tHYPHYValue            = 0.0;

extern long systemCPUCount;

_THyPhy * globalInterfaceInstance = nil;

//_________________________________________________________
// default callback hanlder
bool _tHyPhyDefaultHandler (char*,int,double)
{
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
_THyPhyString* _THyPhyReturnObject::castToString(void)
{
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
_THyPhyString::_THyPhyString(const char* characters, long length)
{
    if (characters) {
        if (length == 0) {
            while (characters[length++]) ;
            length --;
        }
        checkPointer (sData = (char*)MemAllocate (length+1));
        memcpy       (sData,characters,length+1);

    } else {
        sData   = nil;
    }
    sLength = length;
}

//_________________________________________________________
_THyPhyString::~_THyPhyString(void)
{
    if (sData) {
        free (sData);
    }
};

//_________________________________________________________
_THyPhyNumber::_THyPhyNumber(double v)
{
    nValue = v;
}


//_________________________________________________________
_THyPhyMatrix::_THyPhyMatrix(void)
{
    mData = nil;
    mRows = 0;
    mCols = 0;
}

//_________________________________________________________
_THyPhyMatrix::_THyPhyMatrix(const long r, const long c, const double* d)
{
    checkPointer (mData = (double*)MemAllocate (r*c*sizeof(double)));
    mRows = r;
    mCols = c;
    for (long i = 0; i < r*c; i++,d++) {
        mData[i] = *d;
    }

}

//_________________________________________________________
_THyPhyMatrix::~_THyPhyMatrix(void)
{
    if (mData) {
        free (mData);
    }
};

//_________________________________________________________
double _THyPhyMatrix::MatrixCell(long r, long c)
{
    return mData[mCols*r+c];
}

//_________________________________________________________
// Begin _THyPhy definitions
//_________________________________________________________

_THyPhy::_THyPhy(_ProgressCancelHandler* mHandler, const char* baseDirPath, long cpuCount)
{
    InitTHyPhy (mHandler, baseDirPath, cpuCount);
}

//_________________________________________________________


_THyPhy::_THyPhy(const char* baseDirPath, long cpuCount)
{
    InitTHyPhy (_tHyPhyDefaultHandler, baseDirPath, cpuCount);
}

//_________________________________________________________

_THyPhy::~_THyPhy           (void)
{
    if (currentResultHolder) {
        delete (_THyPhyReturnObject*)currentResultHolder;
    }
    if (baseDirectoryInstance) {
        delete (baseDirectoryInstance);
    }
    ClearAll();
    DeleteObject        ((_String*)errors);
    DeleteObject        ((_String*)warnings);
    DeleteObject        ((_String*)textout);
    if (globalInterfaceInstance == this) {
        globalInterfaceInstance = nil;
    }

    PurgeAll(true);
    GlobalShutdown();
}

//_________________________________________________________

void _THyPhy::InitTHyPhy (_ProgressCancelHandler* mHandler, const char* baseDirPath, long cpuCount)
{
    char dirSlash = GetPlatformDirectoryChar ();
    systemCPUCount = cpuCount;
    SetCallbackHandler (mHandler);
    checkPointer (currentResultHolder = new _THyPhyString);
    askFID = -1;
    if (baseDirPath)
        // set base directory
    {
        baseDirectory = baseDirPath;
        if (baseDirectory.getChar(baseDirectory.sLength-1) != dirSlash) {
            baseDirectory = baseDirectory & dirSlash;
        }
        baseDirectoryInstance = new _THyPhyString (baseDirectory.sData);
        baseDirectory = baseDirectoryInstance->sData;
        pathNames && &baseDirectory;
        ReadPreferences ();
    }

#ifdef _HYPHY_LIBDIRECTORY_    
    libDirectory = _HYPHY_LIBDIRECTORY_;
    if (libDirectory.getChar(libDirectory.sLength-1) != dirSlash) {
        libDirectory = libDirectory & dirSlash;
    }
    pathNames && &libDirectory;
#else
    if (baseDirectory)
        libDirectory = baseDirectory;
#endif

    GlobalStartup();
    errors   = nil;
    warnings = nil;
    textout  = nil;
    globalInterfaceInstance = this;

}

//_________________________________________________________

_THyPhyString * _THyPhy::ExecuteBF (const char * buffer, bool doPurge)
{
    if (doPurge) {
        PurgeAll            (true);    // cleanup results of previous analysis
    }
    _String             commandString (buffer);

    if (commandString.beginswith ("#NEXUS"),false) {
        lastNexusDataMatrix = ReadDataSetFile (nil, 2, &commandString);
        commandString = nexusBFBody;
    }

    _ExecutionList      compiledCode  (commandString);

    ApplyPreferences    ();

    DeleteObject        ((_String*)errors);
    DeleteObject        ((_String*)warnings);
    DeleteObject        ((_String*)textout);

    errors              = new _String (128L,true);
    warnings            = new _String (128L,true);
    textout             = new _String (128L,true);

    askFID              = compiledCode.ExecuteAndClean (0x7ffffff,&_tHYPHYAskFor);
    _PMathObj bfReturn  = compiledCode.GetResult ();

    ((_String*)errors)->Finalize();
    ((_String*)warnings)->Finalize();
    ((_String*)textout)->Finalize();

    if (currentResultHolder->sData) {
        free (currentResultHolder->sData);
        currentResultHolder->sData = nil;
    }
    if (bfReturn) {
        _String * serializedReturn   = (_String*) bfReturn->toStr();
        currentResultHolder->sData   = serializedReturn->sData;
        serializedReturn->sData      = nil;
        currentResultHolder->sLength = serializedReturn->sLength;
    }
    return currentResultHolder;
}

//_________________________________________________________

void _THyPhy::ClearAll (void)
{
    PurgeAll(true);
}

//_________________________________________________________

void* _THyPhy::AskFor (const char* resultID)
{
    if (resultID && askFID >= 0) {
        _String theCommand (128L,true);
        theCommand << "return ";
        theCommand <<  _tHYPHYAskFor;
        theCommand << "(\"";
        theCommand.EscapeAndAppend (resultID);
        theCommand << "\");";
        theCommand.Finalize();
        _ExecutionList      compiledCode  (theCommand);
        compiledCode.ExecuteAndClean (0x7ffffff);
        _PMathObj retResult = compiledCode.GetResult ();
        if (retResult && retResult->ObjectClass() == STRING) {
            _FString * checkHandled = (_FString*)retResult;
            if (checkHandled->theString->Equal (&_tHYPHYNotHandled)) {
                return nil;
            }
        }
        return retResult->makeDynamic();
    }
    return nil;
}

//_________________________________________________________

void _THyPhy::DumpResult (void* aResult)
{
    if (aResult) {
        DeleteObject ((BaseRef)aResult);
    }
}

//_________________________________________________________

void _THyPhy::SetCallbackHandler (_ProgressCancelHandler* newHandler)
{
    theHandler = newHandler;
}

//_________________________________________________________

_ProgressCancelHandler* _THyPhy::GetCallbackHandler (void)
{
    return theHandler;
}

//_________________________________________________________

void        _THyPhy::PushWarning (void * o)
{
    if (warnings) {
        *((_String*)warnings) << *(_String*)o;
    }
}

//_________________________________________________________

void        _THyPhy::PushError (void * o)
{
    if (errors) {
        *((_String*)errors) << *(_String*)o;
    }
}

//_________________________________________________________

void        _THyPhy::PushOutString (void * o)
{
    if (textout) {
        *((_String*)textout) << *(_String*)o;
    }
}

//_________________________________________________________

_THyPhyString   * _THyPhy::ConvertHyPhyString (void * o)
{
    return new _THyPhyString (((_String*)o)->sData,((_String*)o)->sLength);
}

//_________________________________________________________

bool _THyPhy::CanCast (const void* theObject, const int requestedType)
{
    if (theObject) {
        _PMathObj mo = (_PMathObj) theObject;
        switch (((_PMathObj)theObject)->ObjectClass()) {
        case NUMBER:
            return true;
            // can cast a number to everything
        case STRING:
            return requestedType!=THYPHY_TYPE_MATRIX;
            // can cast anything to a string
        case MATRIX:
            return requestedType!=THYPHY_TYPE_NUMBER;
            // can not cast matrix to number

        case TREE:
        case TOPOLOGY:
            return requestedType==THYPHY_TYPE_STRING;

        }
    }
    return false;
}

//_________________________________________________________

_THyPhyReturnObject* _THyPhy::CastResult (const void* theObject, const int requestedType)
{
    _THyPhyReturnObject * convertedObject = nil;
    if (CanCast(theObject,requestedType)) {
        int hyphyObjClass = ((_PMathObj)theObject)->ObjectClass();
        switch (hyphyObjClass) {
        case NUMBER: {
            if (hyphyObjClass == NUMBER) {
                return new _THyPhyNumber (((_PMathObj)theObject)->Compute()->Value());
            }
            if (hyphyObjClass == STRING) {
                _String sV ((_String*)((_FString*)theObject)->toStr());
                return new _THyPhyNumber (sV.toNum());
            }
        }
        case STRING: {
            _String sV ((_String*)((_PMathObj)theObject)->toStr());
            return new _THyPhyString (sV.sData,sV.sLength);
        }
        case MATRIX: {
            if (hyphyObjClass == NUMBER) {
                double evaluate = ((_PMathObj)theObject)->Compute()->Value();
                return new _THyPhyMatrix (1,1,&evaluate);
            }

            if (hyphyObjClass == MATRIX) {
                _Matrix * evalutedNumeric =  (_Matrix*)((_Matrix*)(((_PMathObj)theObject)->Compute()))
                                             ->ComputeNumeric();

                return new _THyPhyMatrix (evalutedNumeric->GetHDim(),evalutedNumeric->GetVDim(),evalutedNumeric->theData);
            }
        }
        }
    }

    return convertedObject;
}

//_________________________________________________________________________

void    SetStatusLine (_String arg)
{
    _tHYPHYCurrentStatus = arg;
    yieldCPUTime();
}

//_________________________________________________________________________

void    SetStatusBarValue (long l,_Parameter max, _Parameter rate)
{
    _tHYPHYDone = l;
    _tHYPHYCurrentStatus   = _String ("LF Optimization. Value=") & _String (max) &", "&_String (rate) & " evals/sec.";
    _tHYPHYValue = max;
    yieldCPUTime();
}

//_________________________________________________________________________

long    _THyPhyGetLongStatus        (void)
{
    return _tHYPHYDone;
}

//_________________________________________________________________________

double  _THyPhyGetDoubleStatus      (void)
{
    return _tHYPHYValue;
}

//_________________________________________________________________________

char*   _THyPhyGetStringStatus      (void)
{
    return _tHYPHYCurrentStatus.getStr();
}
