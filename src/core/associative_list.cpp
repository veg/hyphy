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

/*#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <limits.h>

#include "batchlan.h"
#include "polynoml.h"
#include "likefunc.h"

#include "avllist.h"
#include "avllistx.h"
#include "avllistxl.h"

#include "function_templates.h"
#include "mersenne_twister.h"
#include "global_things.h"
#include "string_file_wrapper.h"*/

#include "global_things.h"
#include "global_object_lists.h"
#include "associative_list.h"
#include "batchlan.h"
#include "avllistxl_iterator.h"


using namespace hy_global;
using namespace hyphy_global_objects;

#define         kAssociativeListDefaultReturn (new _Constant (0.0));




//_____________________________________________________________________________________________
// AssociativeList
//_____________________________________________________________________________________________

_AssociativeList::_AssociativeList (void):avl(&theData) {
}

//_____________________________________________________________________________________________

BaseRef _AssociativeList::makeDynamic (void) const {
    _AssociativeList * newAL = new _AssociativeList ();
    newAL->Duplicate (this);
    return newAL;
}

//_____________________________________________________________________________________________

bool _AssociativeList::ParseStringRepresentation (_String& serialized_form, _FormulaParsingContext& fpc ) {
    _List               splitKeys;
    bool                doErrors = fpc.errMsg() == nil,
                        compute_keys_values = fpc.buildComplexObjects();
   _VariableContainer const* theP = fpc.formulaScope();
  
    _ElementaryCommand::ExtractConditions (serialized_form, 0, splitKeys, ',' , false);
  
    try {
        for (unsigned long k = 0UL; k < splitKeys.countitems(); k ++) {
            _List key_value_pair;
            _ElementaryCommand::ExtractConditions (*(_String*)splitKeys(k), 0, key_value_pair, ':' , false);
            if (key_value_pair.countitems() == 2UL) {
                
                _String  key        (compute_keys_values ? ProcessLiteralArgument((_String*)key_value_pair(0),theP) : *(_String*)key_value_pair(0)),
                         errMsg;
              
                if (key.empty()) {
                  key = *(_String*)key_value_pair(0);
                }
                
                _Formula value      (*(_String*)key_value_pair(1),theP, doErrors?nil :&errMsg);
                _PMathObj   valueC  = compute_keys_values ? value.Compute() : new _MathObject;
              
                if (valueC) {
                    MStore (key, valueC, compute_keys_values);
                } else {
                    throw (((_String*)key_value_pair(1))->Enquote() & " could not be evaluated");

                }
            } else {
                throw (((_String*)splitKeys(k))->Enquote() & " does not appear to specify a valid key:value pair");
            }
        }
    } catch (const _String err) {
        if (doErrors) {
            HandleApplicationError(err);
        }
        return false;
    }
    return true;
}

//_____________________________________________________________________________________________

BaseRef _AssociativeList::toStr (unsigned long padding) {
    return Serialize  (padding);
}

//_____________________________________________________________________________________________

void _AssociativeList::Duplicate (BaseRefConst br) {
    if (!SingleReference ()) {
        HandleApplicationError(_String (__PRETTY_FUNCTION__).Enquote() & " called from an object with multiple references");
    }
    _AssociativeList const * copyMe = (_AssociativeList const*)br;
    theData.Duplicate (&copyMe->theData);
    avl.leftChild.Duplicate (&copyMe->avl.leftChild);
    avl.rightChild.Duplicate (&copyMe->avl.rightChild);
    avl.balanceFactor.Duplicate (&copyMe->avl.balanceFactor);
    avl.emptySlots.Duplicate (&copyMe->avl.emptySlots);
    avl.xtraD.Duplicate (&copyMe->avl.xtraD);
    avl.root = copyMe->avl.root;
}

//_____________________________________________________________________________________________

_PMathObj _AssociativeList::MCoord (_PMathObj p) {
    return new _FString ((_String*)p->toStr());
}

//_____________________________________________________________________________________________
_PMathObj _AssociativeList::MAccess (_PMathObj p) {
    long        f;

    if (p->ObjectClass() == STRING) {
        f = avl.Find (&((_FString*)p)->get_str());
    } else {
        _String s ((_String*)p->toStr());
        f = avl.Find (&s);
    }
    if (f>=0L) {
        _PMathObj res = (_PMathObj)avl.GetXtra (f);
        res->AddAReference();
        return res;
    } else {
        return kAssociativeListDefaultReturn;
    }
}

//_____________________________________________________________________________________________
_PMathObj _AssociativeList::MIterator (_PMathObj p, _PMathObj p2) {
    
    const _String     kAVLIteratorOrder          = "INDEXORDER",
                      kAVLIteratorOrderValue     = "VALUEINDEXORDER";

    long done = 0;

    _List reference_manager;
    
    try {

        if (p->ObjectClass() == STRING && p2->ObjectClass() == STRING) {

            long avlRoot = avl.GetRoot();

            if (avlRoot >= 0) {
                
                _String * callback_id  = (_String*)p->toStr(),
                        * filter_id    = (_String*)p2->toStr();
                
                reference_manager < callback_id < filter_id;

                long    callback  = FindBFFunctionName (*callback_id),
                        filter    = FindBFFunctionName (*filter_id);

                if (callback < 0L || GetBFFunctionArgumentCount(callback) != 2L) {
                    throw ("The first argument in an iterator call for Associative Arrays must be a valid identifier of a function taking two arguments (key, value)");
                } else {
                    if (filter >= 0L && GetBFFunctionArgumentCount (filter) != 1L) {
                        throw ("The second argument in an iterator call for Associative Arrays must be either empty or a valid identifier of a function taking a single argument");
                    }

                    _Formula      testFormula,
                                  actionFormula;

                    actionFormula.GetList() < new _Operation()
                                            < new _Operation()
                                            < new _Operation(kEmptyString,-callback-1L);

                    if (filter >= 0L) {
                        testFormula.GetList() < new _Operation() < new _Operation(kEmptyString,-filter-1L);
                    }

                    _FString fKey;
                    
                    for (AVLListXLIteratorKeyValue filter_key_value : AVLListXLIterator (&avl)) {
                        _String * current_key = (_String *)avl.Retrieve (filter_key_value.get_index());
                        if (current_key) {
                            fKey.SetStringContent (new _StringBuffer (*current_key));
                            if (filter >= 0L) {
                                testFormula.GetIthTerm(0)->SetNumber(&fKey);
                                if (CheckEqual(testFormula.Compute()->Value(),0.0)) {
                                    continue;
                                }
                            }
                            actionFormula.GetIthTerm(0)->SetNumber(&fKey);
                            actionFormula.GetIthTerm(1)->SetNumber((_PMathObj)filter_key_value.get_object());
                            actionFormula.Compute();
                            done ++;
                        }
                    }
                    
                    actionFormula.GetIthTerm(0)->SetNumber(nil);
                    actionFormula.GetIthTerm(1)->SetNumber(nil);
                    if (filter >= 0) {
                        testFormula.GetIthTerm(0)->SetNumber(nil);
                    }
                    
                }
             }
        } else if (p->ObjectClass () == STRING && p2->ObjectClass () == NUMBER) {
            _String * mode  = (_String*)p->toStr();
            reference_manager < mode;
            _PMathObj result = nil;

            if (*mode == kAVLIteratorOrder|| *mode == kAVLIteratorOrderValue) {
                long index = avl.GetByIndex(p2->Compute()->Value());
              
                if (index >= 0L) {
                  result = (*mode == kAVLIteratorOrder )? (new _FString(*(_String*)avl.Retrieve(index),false)): ((_PMathObj)avl.GetXtra (index)->makeDynamic());
                } else {
                  throw ("Index out of bounds in call to AVL iterator (by index)");
                }
            }
          
            if (result) {
              return result;
            }
        } else {
            throw ("Both arguments must be Strings (or a String Literal and a number) in an iterator call for Associative Arrays");
        }
    } catch (const _String err) {
        HandleApplicationError (err);
    }
    return new _Constant (done);
}

//_____________________________________________________________________________________________
_PMathObj _AssociativeList::GetByKey (_String const& key, long objType) const {
    _PMathObj res = GetByKey(key);
    if (res && ((res->ObjectClass() & objType) > 0L)) {
      return res;
    }
    return nil;
}

//_____________________________________________________________________________________________
_PMathObj _AssociativeList::GetByKey (_String const& key) const {
    return (_PMathObj)avl.GetDataByKey(&key);
}

//_____________________________________________________________________________________________
_PMathObj _AssociativeList::GetByKey (long nKey, long objType) const {
    return GetByKey (_String(nKey), objType);
}

  //_____________________________________________________________________________________________
void _AssociativeList::Clear (void) {
  avl.Clear(true);
}


//_____________________________________________________________________________________________
void _AssociativeList::DeleteByKey (_PMathObj p) {
    if (p->ObjectClass() == STRING) {
        avl.Delete (&((_FString*)p)->get_str(),true);
    } else {
        if (p->ObjectClass() == ASSOCIATIVE_LIST) {
            _List * keys2remove = ((_AssociativeList*)p)->GetKeys();
            for (long ki = 0; ki < keys2remove->lLength; ki++) {
                avl.Delete ((_String*)(*keys2remove)(ki),true);
            }
            DeleteObject (keys2remove);
        } else {
            _String * s = (_String*)p->toStr();
            avl.Delete (s,true);
            DeleteObject (s);
        }
    }
}


//_____________________________________________________________________________________________
void _AssociativeList::MStore (_PMathObj p, _PMathObj inObject, bool repl, long opCode) {
    if (!p) {
        return;
    }

    _FString * index = (_FString*)p;
    long       f     = avl.Find (&index->get_str());

    if (f>=0) { // already exists - replace
        if (opCode == HY_OP_CODE_ADD) {
            _List arguments;
            arguments << inObject;
          
            _PMathObj newObject = ((_PMathObj)avl.GetXtra(f))->ExecuteSingleOp(HY_OP_CODE_ADD,&arguments);
            if (repl == false) {
                DeleteObject (inObject);
            } else {
                repl = false;
            }
            inObject = newObject;
        }
        avl.xtraD.Replace (f, inObject, repl);
    } else { // insert new
        if (repl) {
            BaseRef br = inObject->makeDynamic();
            avl.Insert (index->get_str().makeDynamic(),(long)br,false);
            //br->nInstances--;
        } else {
            avl.Insert (index->get_str().makeDynamic(),(long)inObject,false);
        }
    }
}

//_____________________________________________________________________________________________
void _AssociativeList::MStore (const _String& obj, _PMathObj inObject, bool repl) {
    _FString f (obj);
    MStore (&f,inObject, repl);
}

//_____________________________________________________________________________________________

_AssociativeList &  _AssociativeList:: operator <<     (_associative_list_key_value pair) {
  if (pair.key) {
    MStore (pair.key, pair.payload, true);
  } else {
    _String next_index ((long)Length());
    MStore (next_index, pair.payload, true);
  }
  return *this;
}

//_____________________________________________________________________________________________


_AssociativeList &  _AssociativeList:: operator <     (_associative_list_key_value pair) {
  if (pair.key) {
    MStore (pair.key, pair.payload, false);
  } else {
    _String next_index ((long)Length());
    MStore (next_index, pair.payload, false);
  }
  return *this;
}


//_____________________________________________________________________________________________
void _AssociativeList::MStore (const _String& obj, const _String& info) {
    _FString *  inf = new _FString (info);
    MStore (obj, inf, false);
}


//_____________________________________________________________________________________________
_StringBuffer * _AssociativeList::Serialize (unsigned long padding) const {
  
    _StringBuffer  * out_string = new _StringBuffer (1024L);
    _String          padder (" ", padding);
  
  
  
    (*out_string) << "{";
    bool        doComma = false;
    _List * meKeys = GetKeys();
    for (long k = 0L; k < meKeys->countitems(); k++) {
        _String   *thisKey  = (_String*)(*meKeys)(k);
        if (thisKey) {
            if (doComma) {
                (*out_string) << ',';
             }
          
            (*out_string) << '\n' << padder << ' ';
          
            (*out_string) << '"';
            out_string->SanitizeAndAppend (*thisKey);
            (*out_string) << '"';

            _PMathObj anObject = GetByKey (*thisKey);

            (*out_string) << ':';
            if (anObject->ObjectClass() == STRING) {
                (*out_string) << '"';
                out_string->SanitizeAndAppend(_String ((_String*)anObject->toStr(padding+2)));
                (*out_string) << '"';
            } else {
                out_string->AppendNewInstance((_String*)anObject->toStr(padding+2));
            }
            doComma = true;
        }
    }
    
    DeleteObject (meKeys);
    
    (*out_string) << '\n' << padder << "}";
    out_string->TrimSpace ();
    
    return out_string;
}


//__________________________________________________________________________________
_PMathObj _AssociativeList::Compute (void) {
    return this;
}

//__________________________________________________________________________________
_List* _AssociativeList::GetKeys (void) const {
    
    _List * keys = new _List;
    
    for (AVLListXLIteratorKeyValue key_value : AVLListXLIterator (&avl)) {
        (*keys) << avl.Retrieve(key_value.get_index());
    }

    return keys;
    
}

//_____________________________________________________________________________________________
void        _AssociativeList::FillInList (_List& fill_me) {
  
    unsigned long ll = fill_me.countitems();
    // checkpoint the length of the list
    try {
      // try filling in in numerical order first 0, 1, ...
      // failing that (if the keys are not numeric)
      // fill in by alphanumeric keys
      unsigned long  my_length = avl.countitems();
      
      for (long index = 0L; index < my_length; index++) {
        _String key (index);
        if (_PMathObj value = GetByKey (key)) {
          fill_me.AppendNewInstance(value->toStr());
        } else {
          throw (1);
        }
      }
    } catch (int e) {
        while (fill_me.countitems () > ll) {
          fill_me.Delete(fill_me.countitems ()-1);
        }

        for (AVLListXLIteratorKeyValue key_value : AVLListXLIterator (&avl)) {
            _String* aKey = (_String*)avl.Retrieve(key_value.get_index());
            if (aKey) {
                fill_me.AppendNewInstance(key_value.get_object()->toStr());
            }
        }
    }
}


//_____________________________________________________________________________________________
void        _AssociativeList::Merge (_PMathObj p) {
    //SW20111207: I don't think we should ever have to worry about avl traversing 
    //here as long as the other methods are implemented properly


    if(p==this){
        return;
    }

    if (p && p->ObjectClass() == ASSOCIATIVE_LIST) {
     _AssociativeList *rhs = (_AssociativeList*) p;
      if (rhs->avl.countitems()) {
          for (AVLListXLIteratorKeyValue key_value : AVLListXLIterator (&rhs->avl)) {
              MStore(*(_String*)rhs->avl.Retrieve (key_value.get_index()),(_PMathObj)key_value.get_object(),true);
          }
      }
    }
    else {
        HandleApplicationError ("Associative list merge operation requires an associative list argument.");
    }
}

  //_____________________________________________________________________________________________
_PMathObj        _AssociativeList::ExtremeValue (bool do_minimum) const {
  _String const * best_key = nil;
  hyFloat best_value = do_minimum ? INFINITY : -INFINITY;
  
  if (avl.countitems()) {
      for (AVLListXLIteratorKeyValue key_value : AVLListXLIterator (&avl)) {
          _PMathObj value = (_PMathObj)key_value.get_object();
          switch (value->ObjectClass()){
              case NUMBER:
                  hyFloat number = ((_Constant*)value)->Value();
                  if (do_minimum) {
                      if (number < best_value) {
                          best_value = number;
                          best_key   = (_String const*)avl.Retrieve (key_value.get_index());
                      }
                  } else {
                      if (number > best_value) {
                          best_value = number;
                          best_key   = (_String const*)avl.Retrieve (key_value.get_index());
                      }
                  }
                  break;
          }
      }
  }

  
  _AssociativeList * result = new _AssociativeList;
  (*result) < _associative_list_key_value {"key", best_key ? new _FString (*best_key, false) : new _MathObject}
            < _associative_list_key_value {"value", new _Constant (best_value)};
  return result;
  
}

//_____________________________________________________________________________________________
_PMathObj        _AssociativeList::Sum (void) {
    hyFloat sum = 0.;
        
    for (AVLListXLIteratorKeyValue key_value : AVLListXLIterator (&avl)) {
        _PMathObj value = (_PMathObj)key_value.get_object();
        switch (value->ObjectClass()){
            case NUMBER:
                sum += ((_Constant*)value)->Value();
                break;
            case STRING:
                sum += ((_FString*)value)->get_str().to_float();
                break;
            case MATRIX: {
                _Constant * sumOfValue = (_Constant*) ((_Matrix*)value->Compute())->Sum();
                sum += sumOfValue->Value();
                DeleteObject (sumOfValue);
                break;
            }
            case ASSOCIATIVE_LIST: {
                _Constant * sumOfValue = (_Constant*) ((_AssociativeList*)value->Compute())->Sum();
                sum += sumOfValue->Value();
                DeleteObject (sumOfValue);
                break;
            }
        }
    }
    
    return new _Constant (sum);
}

//__________________________________________________________________________________


_PMathObj _AssociativeList::ExecuteSingleOp (long opCode, _List* arguments, _hyExecutionContext* context)  {
  
  switch (opCode) {
    case HY_OP_CODE_ABS:
      return new _Constant (Length());
      
    case HY_OP_CODE_EVAL:
      return (_PMathObj) makeDynamic();
      
    case HY_OP_CODE_COLUMNS: {
      // Columns -- get all unique values (as strings)
      _List    unique_values_aux;
      _AVLList unique_values (&unique_values_aux);
      
      for (unsigned long k=0UL; k<avl.dataList->lLength; k++) {
        BaseRef anItem = ((BaseRef*)avl.dataList->lData)[k];
        if (anItem) {
          _String* string_value = (_String*) avl.GetXtra(k)->toStr();
          if (unique_values.Insert (string_value, 0L, false) < 0) {
            DeleteObject(string_value);
          }
        }
      }
      unique_values.ReorderList();
      return new _Matrix (*(_List*)unique_values.dataList);
    }
      
    case HY_OP_CODE_ROWS: {
      // Rows - get keys
      if (avl.emptySlots.lLength) {
        _List  dataListCompact;
        for (long k=0; k<avl.dataList->lLength; k++) {
          BaseRef anItem = ((BaseRef*)avl.dataList->lData)[k];
          if (anItem) {
            dataListCompact << anItem;
          }
        }
        return new _Matrix (dataListCompact);
      }
      return new _Matrix (*(_List*)avl.dataList);
    }
      
    case HY_OP_CODE_TYPE: // Type
      return Type();
      
    case HY_OP_CODE_MAX: // Max
      return ExtremeValue (false);
      
    case HY_OP_CODE_MIN: // Max
      return ExtremeValue (true);
     

      
  }
  
  _MathObject * arg0 = _extract_argument (arguments, 0UL, false);
  
  switch (opCode) { // next check operations without arguments or with one argument
    case HY_OP_CODE_ADD: // +
      if (arg0) {
        MStore (_String((long)avl.countitems()), arg0, true);
        return new _Constant (avl.countitems());
      }
      return Sum ();
  }
  

  if (arg0) {
    switch (opCode) { // operations that require exactly one argument
      case HY_OP_CODE_MCOORD: // MCoord
        return MCoord (arg0);
        
      case HY_OP_CODE_MUL: // merge
        Merge (arg0);
        return new _Constant (avl.countitems());
        
      case HY_OP_CODE_SUB:
        DeleteByKey (arg0);
        return new _Constant (avl.countitems());
        
      case HY_OP_CODE_DIV:
        
        if (arg0->ObjectClass () == STRING) {
          if (avl.Find (&((_FString*)arg0)->get_str()) >= 0L) {
            return new _Constant (1.0);
          }
        } else {
          _String serialized ((_String*)arg0->toStr());
          if (avl.Find (&serialized) >= 0L) {
            return new _Constant (1.0);
          }
        }
        return new _Constant (0.0);
        
    }
    _MathObject * arg1 = _extract_argument (arguments, 1UL, false);
    
    switch (opCode) { //  check operations with 1 or 2 arguments
      case HY_OP_CODE_MACCESS: // MAccess
        if (arg1) {
          return MIterator (arg0,arg1);
        } else {
          return MAccess   (arg0);
        }
    }
  }
  
  
  switch (opCode) {
    case HY_OP_CODE_TYPE:
    case HY_OP_CODE_ADD:
    case HY_OP_CODE_MCOORD:
    case HY_OP_CODE_MUL:
    case HY_OP_CODE_SUB:
    case HY_OP_CODE_MACCESS:
    case HY_OP_CODE_DIV:
      WarnWrongNumberOfArguments (this, opCode,context, arguments);
      break;
    default:
      WarnNotDefined (this, opCode,context);
  }
  
  return new _MathObject;
  
}


