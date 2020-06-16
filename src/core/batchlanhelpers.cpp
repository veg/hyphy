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

#include      "batchlan.h"
#include      "defines.h"
#include      "global_object_lists.h"
#include      "global_things.h"

using namespace hy_global;

_Trie   _HY_HBL_Namespaces;
_List   templateModelList;

extern  _List batchLanguageFunctionNames;

//____________________________________________________________________________________

_String const   _HYGenerateANameSpace () {
     
    _String nmsp,
            capLetters ("ABCDEFGHIJKLMNOPQRSTUVWXYZ_abcdefghijklmnopqrstuvwxyz");
    do {
        nmsp = _String::Random (8, &capLetters);

    } while (_HY_HBL_Namespaces.FindKey (nmsp) != kNotFound);

    _HY_HBL_Namespaces.Insert (nmsp, 0L);
    
    return nmsp;
}

//____________________________________________________________________________________

void  _HYClearANameSpace (const _String& nm) {
    _HY_HBL_Namespaces.Delete(nm);
 }

//____________________________________________________________________________________


void    ReadModelList(void)
{
    if (templateModelList.empty() == false) return;

    _String     modelListFile (GetStandardDirectory (HY_HBL_DIRECTORY_TEMPLATE_MODELS) & "models.lst"),
                theData (doFileOpen (modelListFile.get_str(),"rb"));
  
    if (theData.nonempty()) {
        _ElementaryCommand::ExtractConditions(theData,0,templateModelList);
        for ( long i = 0L; i<templateModelList.countitems(); i++) {
            _String* thisString = (_String*)templateModelList(i);
            _List   thisModel;
            _ElementaryCommand::ExtractConditions(*thisString,thisString->FirstNonSpaceIndex(),thisModel,',');
            if (thisModel.countitems() != 5UL) {
                templateModelList.Delete(i);
                i--;
                continue;
            }
          
            for (long j = 0L; j<5UL; j++) {
                ((_String*)thisModel.GetItem (j))->StripQuotes();
            }
            *(_String*)thisModel(0) = ((_String*)thisModel(0))->ChangeCase(kStringUpperCase);
            templateModelList.Replace(i,&thisModel,true);
        }
    }
}

//____________________________________________________________________________________


bool    ExpressionCalculator (void) {
  
  const static _String kExit ("exit");
  
  _String data (StringFromConsole());
  
    //Checking for exit
  if (data.CompareIgnoringCase(kExit) == kCompareEqual) {
    return false;
  }
  
  _Formula  lhs,
  rhs;
  
  _FormulaParsingContext fpc;
  long retCode = Parse(&lhs, data, fpc, nil);
    
  //printf ("%s\n", _String ((_String*)lhs.GetList().toStr(kFormulaStringConversionNormal)).get_str());
  
  if (!terminate_execution) {
    if (retCode == HY_FORMULA_EXPRESSION) {
      HBLObjectRef formRes = lhs.Compute();
      if (!formRes) {
        BufferToConsole ("NULL\n");
      } else {
        StringToConsole (_String ((_String*)formRes->toStr()));
      }
    } else {
      BufferToConsole ("NO RETURN VALUE");
    }
  }
  NLToConsole();
  terminate_execution = false;
  return true;
}

//____________________________________________________________________________________


bool    PushFilePath (_String& pName, bool trim, bool process) {
  
  long f;
  
  if (process) {
    ProcessFileName(pName);
    f = pName.FindBackwards(get_platform_directory_char());
  } else {
    f = pName.length();
  }
  
  if (f>=0L) {
    pathNames < new _String (pName, 0, f);
    if (trim)
      pName.Trim (f+1L,kStringEnd);
    return true;
  } else if (pathNames.lLength) {
    pathNames && pathNames(pathNames.lLength-1);
  } else {
    pathNames.AppendNewInstance (new _String(kEmptyString));
  }
  
  return false;
}

//____________________________________________________________________________________


const _String   PopFilePath (void) {
    if (pathNames.empty()) {
      return kEmptyString;
    }
    _String top = *(_String*)pathNames.GetElement (-1L);
    pathNames.Delete (pathNames.lLength-1);
    //fprintf (stderr, "\nPopping %s\n", top.sData);
    return top;
}

//____________________________________________________________________________________


const _String   GetPathStack (const _String& spacer)  {
  return _String ((_String*) pathNames.Join (spacer));
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
    terminate_execution = false;
    ex.Execute      ();
    terminate_execution = false;
}

//____________________________________________________________________________________

_String const ReturnDialogInput(bool dispPath, _String const * rel_path) {
    bool do_markdown     = hy_env :: EnvVariableTrue(hy_env :: produce_markdown_output);
    NLToConsole ();

    if (do_markdown) {
      BufferToConsole("\n>");
    }
    StringToConsole (dialogPrompt);
    if (dispPath) {
      BufferToConsole (" (`");
        if (rel_path) {
            StringToConsole (*rel_path);
        } else {
              if (PeekFilePath()) {
                    StringToConsole(*PeekFilePath());
                } else {
                    StringToConsole (hy_base_directory);
                }
        }
        BufferToConsole ("`)");
    }
    BufferToConsole (" ");
    return StringFromConsole();
}


//____________________________________________________________________________________

_String const ReturnFileDialogInput(_String const * rel_path) {
    if (currentExecutionList && (currentExecutionList->has_stdin_redirect() || currentExecutionList->has_keyword_arguments())) {
        try {
            _String dialog_string (currentExecutionList->FetchFromStdinRedirect());
            if (dialog_string.nonempty()) {
                hy_env::EnvVariableSet (hy_env::last_raw_file_prompt, new _FString (dialog_string), false);
                return dialog_string;
            }
        } catch (_String const& e) {
            if (e != kNoKWMatch) {
                HandleApplicationError (e);
                hy_env::EnvVariableSet (hy_env::last_raw_file_prompt, new _FString(), false);
                return kEmptyString;
            }
        }
    }

    _String file_path = ReturnDialogInput(true, rel_path);
    hy_env::EnvVariableSet (hy_env::last_raw_file_prompt, new _FString (file_path), false);
    terminate_execution = file_path.empty();
    return file_path;
}


//____________________________________________________________________________________

_String const WriteFileDialogInput(_String const * rel_path) {
  return ReturnFileDialogInput(rel_path);
}



