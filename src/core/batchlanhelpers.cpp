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

#include      "batchlan.h"
#include      "defines.h"


_List   templateModelList;

_String    _HYStandardDirectory (const unsigned long which_one) 
{
    _String dirSpacer (GetPlatformDirectoryChar());

    switch (which_one) {

        case HY_HBL_DIRECTORY_TEMPLATE_MODELS:
            return libDirectory & "TemplateBatchFiles" & dirSpacer & "TemplateModels" & dirSpacer;
    }

    return empty;
}

void    ReadModelList(void) 
{
    if (templateModelList.lLength > 0) return; 
    
    _String     modelListFile (_HYStandardDirectory (HY_HBL_DIRECTORY_TEMPLATE_MODELS) & "models.lst");
    
    FILE* modelList = doFileOpen (modelListFile.getStr(),"rb");
    if (!modelList) {
        return;
    }
    _String theData (modelList);
    fclose (modelList);
    if (theData.sLength) {
        _ElementaryCommand::ExtractConditions(theData,0,templateModelList);
        for (unsigned long i = 0; i<templateModelList.countitems(); i++) {
            _String* thisString = (_String*)templateModelList(i);
            _List   thisModel;
            _ElementaryCommand::ExtractConditions(*thisString,thisString->FirstNonSpaceIndex(),thisModel,',');
            if (thisModel.lLength!=5) {
                templateModelList.Delete(i);
                i--;
                continue;
            }
            for (long j = 0; j<5; j++) {
                ((_String*)thisModel(j))->StripQuotes();
            }
            ((_String*)thisModel(0))->UpCase();
            templateModelList.Replace(i,&thisModel,true);
        }
    }
}

bool ExpressionCalculator (_String data)
{
#ifndef __UNIX__
    if (terminateExecution) {
        return false;
    }
    BufferToConsole (">");
    StringToConsole (data);
    BufferToConsole ("\n");
#endif

    if (data.sLength == 4) {
        _String checkForExit (data);
        checkForExit.LoCase();
        if (checkForExit == _String ("exit")) {
            return false;
        }
    }

    _Formula  lhs,
              rhs;

    long       refV,
               retCode = Parse(&lhs, data, refV, nil,nil);

    if (!terminateExecution) {
        if (retCode == HY_FORMULA_EXPRESSION) {
            _PMathObj formRes = lhs.Compute();
            if (!formRes) {
                BufferToConsole ("NULL\n");
            } else {
                _String * objValue = (_String*)formRes->toStr();
                StringToConsole (*objValue);
                //BufferToConsole ("\n");
                DeleteObject    (objValue);
            }
        } else {
            BufferToConsole ("NO RETURN VALUE");
        }
    }
    terminateExecution = false;
    return true;
}

bool    ExpressionCalculator (void)
{
    _String data (StringFromConsole(false));

#ifndef __UNIX__
    if (terminateExecution) {
        return false;
    }
    BufferToConsole (">");
    StringToConsole (data);
    BufferToConsole ("\n");
#endif

    if (data.sLength == 4) {
        _String checkForExit (data);
        checkForExit.LoCase();
        if (checkForExit == _String ("exit")) {
            return false;
        }
    }

    _Formula  lhs,
              rhs;

    long       refV,
               retCode = Parse(&lhs, data, refV, nil,nil);

    if (!terminateExecution) {
        if (retCode == HY_FORMULA_EXPRESSION) {
            _PMathObj formRes = lhs.Compute();
            if (!formRes) {
                BufferToConsole ("NULL\n");
            } else {
                _String * objValue = (_String*)formRes->toStr();
                StringToConsole (*objValue);
                //BufferToConsole ("\n");
                DeleteObject    (objValue);
            }
        } else {
            BufferToConsole ("NO RETURN VALUE");
        }
    }
    terminateExecution = false;
    return true;
}

bool    PushFilePath (_String& pName, bool trim)
{
    char c = GetPlatformDirectoryChar();

    long    f = pName.FindBackwards(_String(c),0,-1);
    if (f>=0) {
        _String newP = pName.Cut(0,f);
        pathNames && & newP;
        if (trim)
            pName.Trim (f+1,-1);
        return true;
    } else if (pathNames.lLength) {
        pathNames && pathNames(pathNames.lLength-1);
    } else {
        pathNames && & empty;
    }

    return false;
}

void   PopFilePath (void)
{
    pathNames.Delete (pathNames.lLength-1);
}

void   ExecuteBLString (_String& BLCommand, _VariableContainer* theP)
{
    _ExecutionList ex;
    if (theP) {
        ex.SetNameSpace(*theP->GetName());
    }
    ex.BuildList   (BLCommand);
    terminateExecution = false;
    ex.Execute      ();
    terminateExecution = false;
}
