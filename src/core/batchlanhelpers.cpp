/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon@cfenet.ubc.ca)
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

#include      "batchlan.h"
#include      "defines.h"
#include      "global_object_lists.h"
#include      "function_templates.h"

#ifdef __HYPHYQT__
#include "hyphy_qt_helpers.h"
#endif

#if defined __MAC__ || defined __WINDOZE__ || defined __HYPHY_GTK__
#include "HYUtils.h"
#endif

_Trie   _HY_HBL_Namespaces;
_List   templateModelList;

extern  _List batchLanguageFunctionNames;
extern  _String markdownOutput;

//____________________________________________________________________________________

_String    _HYGenerateANameSpace () {
    _String nmsp,
            capLetters ("ABCDEFGHIJKLMNOPQRSTUVWXYZ_abcdefghijklmnopqrstuvwxyz");
    do {
        nmsp = _String::Random (8, &capLetters);
      
    } while (_HY_HBL_Namespaces.FindKey (nmsp) != HY_TRIE_NOTFOUND);
    
    _HY_HBL_Namespaces.Insert (nmsp, 0);
    return nmsp;
}

//____________________________________________________________________________________

_String    _HYStandardDirectory (const unsigned long which_one) {
    _String dirSpacer (GetPlatformDirectoryChar());

    switch (which_one) {

        case HY_HBL_DIRECTORY_TEMPLATE_MODELS:
            return libDirectory & "TemplateBatchFiles" & dirSpacer & "TemplateModels" & dirSpacer;
    }

    return emptyString;
}

//____________________________________________________________________________________


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

//____________________________________________________________________________________


bool ExpressionCalculator (_String data)
{
    //Checking for exit
    #ifndef __HYPHYQT__
        if (data.sLength == 4) {
            _String checkForExit (data);
            checkForExit.LoCase();
            if (checkForExit == _String ("exit")) {
                return false;
            }
        }
    #endif

    _Formula   lhs,
               rhs;
              
    _String    errMsg;
    _FormulaParsingContext fpc (&errMsg, nil);
    
    long       retCode = Parse(&lhs, data, fpc, nil);

    if (retCode != HY_FORMULA_FAILED) {
        if (retCode == HY_FORMULA_EXPRESSION) {
            _PMathObj formRes = lhs.Compute(0,nil,nil,&errMsg);
            if (errMsg.sLength) {
                WarnError(errMsg);
            } else {
                _String * objValue = (_String*)formRes->toStr();
                StringToConsole(*objValue);
                DeleteObject(objValue);
            }
        } else {
            BufferToConsole ("NO RETURN VALUE");
        }
    } else {
        WarnError(errMsg);
    }
    return true;
}

//____________________________________________________________________________________


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

    _FormulaParsingContext fpc;
    long retCode = Parse(&lhs, data, fpc, nil);

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
    NLToConsole();
    terminateExecution = false;
    return true;
}

//____________________________________________________________________________________


bool    PushFilePath (_String& pName, bool trim, bool process) {
  
    //fprintf (stderr, "\nPushing %s\n", pName.sData);
  
    long f;
  
    if (process) {
      _String dir_sep (GetPlatformDirectoryChar());
      pName.ProcessFileName();
      f = pName.FindBackwards(dir_sep,0,-1);
    } else {
      f = pName.Length();
    }
  
  
  
    if (f>=0) {
         pathNames < new _String (pName, 0, f);
        if (trim)
            pName.Trim (f+1,-1);
        return true;
    } else if (pathNames.lLength) {
        pathNames && pathNames(pathNames.lLength-1);
    } else {
        pathNames && & emptyString;
    }

    return false;
}

//____________________________________________________________________________________


const _String   PopFilePath (void) {
    if (pathNames.empty()) {
      return emptyString;
    }
    _String top = *(_String*)pathNames.GetElement (-1L);
    pathNames.Delete (pathNames.lLength-1);
    //fprintf (stderr, "\nPopping %s\n", top.sData);
    return top;
}

//____________________________________________________________________________________


const _String   GetPathStack (const _String spacer)  {
  return _String ((_String*) pathNames.Join (&spacer));
}

//____________________________________________________________________________________


const _String *  PeekFilePath (void) {
  if (pathNames.empty()) {
    return nil;
  }
  return (_String*)pathNames.GetElement (-1L);
}


//____________________________________________________________________________________


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

//____________________________________________________________________________________

_String ReturnDialogInput(bool dispPath)
{
    long do_markdown;
    checkParameter (markdownOutput, do_markdown, 0L);
  
    NLToConsole ();
  
    if (do_markdown) {
      BufferToConsole("\n>");
    }

    StringToConsole (dialogPrompt);
  
    if (dispPath) {
      BufferToConsole (" (`");
      if (PeekFilePath()) {
            StringToConsole(*PeekFilePath());
        } else {
            StringToConsole (baseDirectory);
        }
      BufferToConsole ("`)");
    }
    BufferToConsole (" ");
    return StringFromConsole();
}


//____________________________________________________________________________________

_String ReturnFileDialogInput(void)
{
    if (currentExecutionList && currentExecutionList->stdinRedirect) {
        _String outS (currentExecutionList->FetchFromStdinRedirect());
        if (outS.sLength) {
            return outS;
        }
    }
    
    _String resolvedFilePath;
    
#ifdef __HEADLESS__
    WarnError ("Unhandled standard input call in headless HYPHY. Only redirected standard input (via ExecuteAFile) is allowed");
    return empty;
#endif

#ifdef __MAC__
    resolvedFilePath =  MacSimpleFileOpen();
#endif

#ifdef __WINDOZE__
    resolvedFilePath =  ReturnFileDialogSelectionWin(false);
#endif

#ifdef __HYPHY_GTK__
    if (PopUpFileDialog (dialogPrompt)) {
        resolvedFilePath = *argFileName;
    }
#endif 

#ifdef __HYPHYQT__  
    resolvedFilePath = _hyQTFileDialog (dialogPrompt,empty, false);
#endif
    
#if defined __UNIX__ && ! defined __HYPHYQT__ && ! defined __HYPHY_GTK__
    resolvedFilePath = ReturnDialogInput(true);
#endif

    
    if (resolvedFilePath.sLength == 0) {
        terminateExecution = true;
    }
    
    return resolvedFilePath;
}

//____________________________________________________________________________________

_String ProcessStringArgument (_String* data) {
    if (data->sLength>2) {
        if (data->sData[data->sLength-1]=='_' && data->sData[data->sLength-2]=='_') {
            _String varName (*data,0,data->sLength-3);
            _FString* theVar = (_FString*)FetchObjectFromVariableByType(&varName,STRING);
            if (theVar) {
                return *theVar->theString;
            }
        }
    }
    return emptyString;
}

//____________________________________________________________________________________

_String WriteFileDialogInput(void) {
    if (currentExecutionList && currentExecutionList->stdinRedirect) {
        _String outS (currentExecutionList->FetchFromStdinRedirect());
        if (outS.sLength) {
            return outS;
        }
    }
    
    defFileNameValue = ProcessLiteralArgument (&defFileString,nil);
    _String resolvedFilePath;
    
#ifdef __HEADLESS__
    WarnError ("Unhandled standard input call in headless HYPHY. Only redirected standard input (via ExecuteAFile) is allowed");
    return empty;
#else
  #ifdef __MAC__
      resolvedFilePath =  MacSimpleFileSave();
  #endif

  #ifdef __WINDOZE__

      resolvedFilePath = ReturnFileDialogSelectionWin(true);
  #endif

  #ifdef __HYPHY_GTK__
      if (PopUpFileDialog (dialogPrompt)) {
          resolvedFilePath = *argFileName;
      } 
  #endif
    #ifdef __HYPHYQT__  
        resolvedFilePath = _hyQTFileDialog (dialogPrompt,defFileNameValue, true);
    #endif
            
    #if defined __UNIX__ && ! defined __HYPHYQT__ && ! defined __HYPHY_GTK__
        resolvedFilePath = ReturnDialogInput(true);
    #endif
#endif
    
    if (resolvedFilePath.sLength == 0) {
        terminateExecution = true;
    }
    defFileNameValue = emptyString;
    return resolvedFilePath;

}



