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

#include "associativelist.h"
#include "batchlan.h"

extern _String MATRIX_AGREEMENT, ANAL_COMP_FLAG, ANAL_MATRIX_TOLERANCE,
    PROFILE_MEAN_VAR_MULT, CACHE_FORMULA_DEPENDANCY, BRANCH_LENGTH_STENCIL,
    AVL_ITERATOR_ORDER, AVL_ITERATOR_ORDER_VALUE;

_AssociativeList::_AssociativeList(void) : avl(&theData) {}


//______________________________________________________________________________
BaseRef _AssociativeList::makeDynamic(void) {
  _AssociativeList *newAL = new _AssociativeList();
  newAL->Duplicate(this);
  return newAL;
}

//______________________________________________________________________________
bool _AssociativeList::ParseStringRepresentation(_String &serializedForm,
                                                 bool doErrors,
                                                 _VariableContainer *theP) {
  _List splitKeys;
  _ElementaryCommand::ExtractConditions(serializedForm, 0, splitKeys, ',', false);

  for (unsigned long k = 0; k < splitKeys.lLength; k++) {

    _List aPair;
    _ElementaryCommand::ExtractConditions(*(_String *)splitKeys(k), 0, aPair,':', false);

    if (aPair.lLength == 2) {

      _String key(ProcessLiteralArgument((_String *)aPair(0), theP)), errMsg;
      _Formula value(*(_String *)aPair(1), theP, doErrors ? nil : &errMsg);
      _PMathObj valueC = value.Compute();

      if (valueC) {
        MStore(key, valueC, true);
      } else {
        if (doErrors) {
          WarnError(*(_String *)aPair(1) & " could not be evaluated");
        }
        return false;
      }
    } else {
      if (doErrors) {
        WarnError(*(_String *)splitKeys(k) &
                  " does not appear to specify a valid key:value pair");
      }
      return false;
    }
  }
  return true;
}

//______________________________________________________________________________
_AssociativeList::_AssociativeList (_PMathObj definition): avl(&theData) {
    _SimpleList* defs = (_SimpleList*) definition; 
    
    for (unsigned long i = 0L; i < defs->lLength; i+=2) {
    
        _Formula * key   = (_Formula*) defs->GetElement(i),
                 * obj = (_Formula*) defs->GetElement(i+1);
        
        _PMathObj  key_value = key->Compute(), 
                   obj_value = obj->Compute();
                   
         if (key_value && obj_value) {
            _String string_key ((_String*)key_value->toStr());
            MStore(string_key, obj_value, true);
              
         } else {
           break;
         }
    }
}

//______________________________________________________________________________
BaseRef _AssociativeList::toStr(void) {
  _String defName("_hyphyAssociativeArray");
  return Serialize(defName);
}

//______________________________________________________________________________
void _AssociativeList::Duplicate(BaseRef br) {
  nInstances = 1;
  _AssociativeList *copyMe = (_AssociativeList *)br;
  theData.Duplicate(&copyMe->theData);
  avl.leftChild.Duplicate(&copyMe->avl.leftChild);
  avl.rightChild.Duplicate(&copyMe->avl.rightChild);
  avl.balanceFactor.Duplicate(&copyMe->avl.balanceFactor);
  avl.emptySlots.Duplicate(&copyMe->avl.emptySlots);
  avl.xtraD.Duplicate(&copyMe->avl.xtraD);
  avl.root = copyMe->avl.root;
}

//______________________________________________________________________________
_PMathObj _AssociativeList::MCoord(_PMathObj p) {
  return new _FString((_String *)p->toStr());
}

//______________________________________________________________________________
_PMathObj _AssociativeList::MAccess(_PMathObj p) {

  long f;

  if (p->ObjectClass() == STRING) {
    f = avl.Find(((_FString *)p)->theString);
  } else {
    _String s((_String *)p->toStr());
    f = avl.Find(&s);
  }
  if (f >= 0) {
    _PMathObj res = (_PMathObj) avl.GetXtra(f);
    res->nInstances++;
    return res;
  } else {
    return new _Constant(0.0);
  }

}

//______________________________________________________________________________
_PMathObj _AssociativeList::MIterator(_PMathObj p, _PMathObj p2) {
  long done = 0;

  if (p->ObjectClass() == STRING && p2->ObjectClass() == STRING) {

    long avlRoot = avl.GetRoot();

    if (avlRoot >= 0) {
      _String *s = (_String *)p->toStr(), *s2 = (_String *)p2->toStr();

      long fID = FindBFFunctionName(*s), fID2 = FindBFFunctionName(*s2);

      if (fID < 0 || batchLanguageFunctionParameters.lData[fID] != 2) {
        WarnError("The first argument in an iterator call for Associative "
                  "Arrays must be a valid identifier of a function taking two "
                  "arguments (key, value)");
      } else {
        if (fID2 >= 0 && batchLanguageFunctionParameters.lData[fID2] != 1) {
          WarnError("The second argument in an iterator call for Associative "
                    "Arrays must be either empty or a valid identifier of a "
                    "function taking a single argument");
        }

        _Formula testFormula, actionFormula;

        actionFormula.GetList().AppendNewInstance(new _Operation());
        actionFormula.GetList().AppendNewInstance(new _Operation());
        actionFormula.GetList()
            .AppendNewInstance(new _Operation( 
              _HY_OPERATION_FUNCTION_CALL, fID, 2, NULL));

        if (fID2 >= 0) {
          testFormula.GetList().AppendNewInstance(new _Operation());
          testFormula.GetList()
              .AppendNewInstance(new _Operation( 
              _HY_OPERATION_FUNCTION_CALL, fID2, 1, NULL));
        }

        _SimpleList hist;
        long ls, cn = avl.Traverser(hist, ls, avlRoot);

        _FString *fKey = new _FString;
        while (cn >= 0) {
          _String *aKey = ((_String **)avl.dataList->lData)[cn];
          if (aKey) {
            DeleteObject(fKey->theString);
            fKey->theString = (_String *)aKey->toStr();
            if (fID2 >= 0) {
              ((_Operation **)testFormula.GetList().lData)[0]->SetPayload(fKey);
              if (CheckEqual(testFormula.Compute()->Value(), 0.0)) {
                cn = avl.Traverser(hist, ls);
                continue;
              }
            }
            ((_Operation **)actionFormula.GetList().lData)[0]->SetPayload(fKey);
            ((_Operation **)actionFormula.GetList().lData)[1]
                ->SetPayload((_PMathObj) avl.GetXtra(cn));
            actionFormula.Compute();
            done++;
          }
          cn = avl.Traverser(hist, ls);
        }

        DeleteObject(fKey);

        ((_Operation **)actionFormula.GetList().lData)[0]->SetPayload(nil);
        ((_Operation **)actionFormula.GetList().lData)[1]->SetPayload(nil);
        if (fID2 >= 0) {
          ((_Operation **)testFormula.GetList().lData)[0]->SetPayload(nil);
        }

      }
      DeleteObject(s);
      DeleteObject(s2);
    }
  } else if (p->ObjectClass() == STRING && p2->ObjectClass() == NUMBER) {
    _String *s = (_String *)p->toStr();

    if (s->Equal(&AVL_ITERATOR_ORDER) || s->Equal(&AVL_ITERATOR_ORDER_VALUE)) {
      long index = avl.GetByIndex(p2->Compute()->Value());
      if (index >= 0) {
        return s->Equal(&AVL_ITERATOR_ORDER)
                   ? (new _FString(*((_String **)avl.dataList->lData)[index],false))
                   : ((_PMathObj) avl.GetXtra(index)->makeDynamic());
      } else {
        WarnError("Index out of bounds in call to AVL iterator (by index)");
      }
    }

    DeleteObject(s);
  } else {
    WarnError("Both arguments must be Strings (or a String Literal and a "
              "number) in an iterator call for Associative Arrays");
  }
  return new _Constant(done);
}

//______________________________________________________________________________
_PMathObj _AssociativeList::GetByKey(_String &key, long objType) {

  long f = avl.Find(&key);

  if (f >= 0) {
    _PMathObj res = (_PMathObj) avl.GetXtra(f);
    if (res->ObjectClass() == objType) {
      return res;
    }
  }

  return nil;
}

//______________________________________________________________________________
_PMathObj _AssociativeList::GetByKey(_String &key) {
  long f = avl.Find(&key);

  if (f >= 0) {
    return (_PMathObj) avl.GetXtra(f);
  }

  return nil;
}

//______________________________________________________________________________
_PMathObj _AssociativeList::GetByKey(long nKey, long objType) {
  _String key(nKey);
  return GetByKey(key, objType);
}

//______________________________________________________________________________
void _AssociativeList::DeleteByKey(_PMathObj p) {
  if (p->ObjectClass() == STRING) {
    avl.Delete(((_FString *)p)->theString, true);
  } else {
    if (p->ObjectClass() == ASSOCIATIVE_LIST) {
      _List *keys2remove = ((_AssociativeList *)p)->GetKeys();
      for (long ki = 0; ki < keys2remove->lLength; ki++) {
        avl.Delete((_String *)(*keys2remove)(ki), true);
      }
    } else {
      _String *s = (_String *)p->toStr();
      avl.Delete(s, true);
      DeleteObject(s);
    }
  }
}

//______________________________________________________________________________
void _AssociativeList::MStore(_PMathObj p, _PMathObj inObject, bool repl,
                              long opCode) {
  if (!p) {
    return;
  }

  _FString *index = (_FString *)p;
  long f = avl.Find(index->theString);

  if (f >= 0) { // already exists - replace
    if (opCode == HY_OP_CODE_ADD) {
      _PMathObj newObject =
          ((_PMathObj) avl.GetXtra(f))->Execute(HY_OP_CODE_ADD, inObject);
      if (repl == false) {
        DeleteObject(inObject);
      } else {
        repl = false;
      }
      inObject = newObject;
    }
    avl.xtraD.Replace(f, inObject, repl);
  } else { // insert new
    if (repl) {
      BaseRef br = inObject->makeDynamic();
      avl.Insert(index->theString->makeDynamic(), (long) br, false);
      //br->nInstances--;
    } else {
      avl.Insert(index->theString->makeDynamic(), (long) inObject, false);
    }
  }
}

//______________________________________________________________________________
void _AssociativeList::MStore(_String obj, _PMathObj inObject, bool repl) {
  _FString f(obj);
  MStore(&f, inObject, repl);
}

//______________________________________________________________________________
void _AssociativeList::MStore(_String obj, _String info) {
  _FString inf(info);
  MStore(obj, &inf, true);
}

//______________________________________________________________________________
_String *_AssociativeList::Serialize(_String &avlName) {
  _String *outString = new _String(1024L, true);
  checkPointer(outString);

  (*outString) << "{";
  bool doComma = false;
  _List *meKeys = GetKeys();
  for (long k = 0; k < meKeys->lLength; k = k + 1) {
    _String *thisKey = (_String *)(*meKeys)(k);
    if (thisKey) {
      if (doComma) {
        (*outString) << ',';
        (*outString) << '\n';
      }

      (*outString) << '"';
      outString->EscapeAndAppend(*thisKey, false);
      (*outString) << '"';

      _PMathObj anObject = GetByKey(*thisKey);

      (*outString) << ':';
      if (anObject->ObjectClass() == STRING) {
        (*outString) << '"';
        outString->EscapeAndAppend(_String((_String *)anObject->toStr()), 0);
        (*outString) << '"';
      } else {
        (*outString) << _String((_String *)anObject->toStr());
      }
      doComma = true;
    }
  }
  (*outString) << "}";
  outString->Finalize();
  return outString;
}

//______________________________________________________________________________
_PMathObj _AssociativeList::Compute(void) { return this; }

//______________________________________________________________________________
_List *_AssociativeList::GetKeys(void) { return (_List *)avl.dataList; }

//______________________________________________________________________________
void _AssociativeList::FillInList(_List &fillMe) {
  _SimpleList hist;
  long ls, cn = avl.Traverser(hist, ls, avl.GetRoot());

  while (cn >= 0) {
    _String *aKey = ((_String **)avl.dataList->lData)[cn];
    if (aKey) {
      fillMe.AppendNewInstance(avl.GetXtra(cn)->toStr());
    }
    cn = avl.Traverser(hist, ls);
  }
}

//______________________________________________________________________________
void _AssociativeList::Merge(_PMathObj p) {

  if (p == this) {
    return;
  }

  if (p && p->ObjectClass() == ASSOCIATIVE_LIST) {

    _AssociativeList *rhs = (_AssociativeList *)p;

    _SimpleList hist;
    long ls, cn = rhs->avl.Traverser(hist, ls, rhs->avl.GetRoot());

    /*   SLKP20120111: we need to skip over "blanks" (e.g. resulting from
         previous delete operations)
         here; using the traversal of the second list is the easiest way to go.
         */

    while (cn >= 0) {
      MStore(*(_String *)(*(_List *)rhs->avl.dataList)(cn),
             (_PMathObj) rhs->avl.GetXtra(cn), true);
      cn = rhs->avl.Traverser(hist, ls);
    }
  } else {
    WarnError("Associative list merge operation requires an associative list "
              "argument.");
  }
}

//______________________________________________________________________________
_PMathObj _AssociativeList::Sum(void) {
  _Parameter sum = 0.;

  _SimpleList hist;
  long ls, cn = avl.Traverser(hist, ls, avl.GetRoot());

  while (cn >= 0) {
    _PMathObj value = (_PMathObj) avl.GetXtra(cn);
    switch (value->ObjectClass()) {
    case NUMBER:
      sum += ((_Constant *)value)->Value();
      break;
    case STRING:
      sum += ((_FString *)value)->theString->toNum();
      break;
    case MATRIX: {
      _Constant *sumOfValue = (_Constant *)((_Matrix *)value->Compute())->Sum();
      sum += sumOfValue->Value();
      DeleteObject(sumOfValue);
      break;
    }
    case ASSOCIATIVE_LIST: {
      _Constant *sumOfValue =
          (_Constant *)((_AssociativeList *)value->Compute())->Sum();
      sum += sumOfValue->Value();
      DeleteObject(sumOfValue);
      break;
    }
    }
    cn = avl.Traverser(hist, ls);
  }

  return new _Constant(sum);
}

//______________________________________________________________________________
// execute this operation with the second arg if necessary
_PMathObj _AssociativeList::Execute(long opCode, _PMathObj p, _PMathObj p2,
                                    _hyExecutionContext *context) {
  switch (opCode) {
    case HY_OP_CODE_ADD: // +
      if (p) {
        MStore(_String((long) avl.countitems()), p, true);
        return new _Constant(avl.countitems());
      } else {
        return Sum();
      }
  
    case HY_OP_CODE_MUL: // merge
      Merge(p);
      return new _Constant(avl.countitems());
  
    case HY_OP_CODE_SUB:
    case HY_OP_CODE_ABS:
      if (opCode == HY_OP_CODE_SUB) {
        DeleteByKey(p);
      }
      return new _Constant(avl.countitems());
      break;
  
    case HY_OP_CODE_MACCESS: // MAccess
      if (p2) {
        return MIterator(p, p2);
      } else {
        return MAccess(p);
      }
      break;
  
    case HY_OP_CODE_MCOORD: // MCoord
      return MCoord(p);
      break;
  
    case HY_OP_CODE_ROWS: // Rows - get keys
      if (avl.emptySlots.lLength) {
        _List dataListCompact;
        for (long k = 0; k < avl.dataList->lLength; k++) {
          BaseRef anItem = ((BaseRef *)avl.dataList->lData)[k];
          if (anItem) {
            dataListCompact << anItem;
          }
        }
        return new _Matrix(dataListCompact);
      } else {
        return new _Matrix(*(_List *)avl.dataList);
      }
      break;
  
    case HY_OP_CODE_TYPE: // Type
      return Type();
      break;
  }

  WarnNotDefined(this, opCode, context);
  return nil;

}
