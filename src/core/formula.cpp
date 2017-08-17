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

#include <math.h>
#include <float.h>
#include "defines.h"
#include "formula.h"


#include "parser.h"
#include "batchlan.h"
#include "function_templates.h"

//Constants
extern _Parameter twoOverSqrtPi;

extern  _Variable*  _x_, *_n_;

extern _Parameter machineEps;
extern _Parameter tolerance;

extern _String intPrecFact;
extern _String intMaxIter;

#define         maxRombergSteps  8L
#define         integrationPrecisionFactor  1.e-5

_Formula::_Formula (void)
{
    theTree     = nil;
    resultCache = nil;
    call_count = 0UL;
    recursion_calls = nil;
}
//__________________________________________________________________________________

_Formula::_Formula (_PMathObj p, bool isAVar)
{
    theTree     = nil;
    resultCache = nil;
    recursion_calls = nil;
    call_count = 0UL;
  
    if (!isAVar) {
        theFormula.AppendNewInstance (new _Operation (p));
    } else {
        _Variable* v = (_Variable*)p;
        theFormula.AppendNewInstance (new _Operation (true,*v->GetName(),v->IsGlobal(), nil));
    }
}
//__________________________________________________________________________________
void _Formula::Initialize (void) {}

//__________________________________________________________________________________
void _Formula::Duplicate  (BaseRef f)
{
    _Formula * f_cast = (_Formula*) f;

    theFormula.Duplicate       (& f_cast->theFormula);
    theStack.theStack.Duplicate(& f_cast->theStack.theStack);
    call_count = f_cast->call_count;
    if (f_cast->recursion_calls) {
      recursion_calls = (_PMathObj)f_cast->recursion_calls->makeDynamic();
    } else {
      recursion_calls = nil;
    }
    if (f_cast->theTree) {
        theTree = f_cast->theTree->duplicate_tree();
    } else {
        theTree = nil;
    }

    if (f_cast->resultCache) {
        resultCache = (_List*)f_cast->resultCache->makeDynamic();
    } else {
        resultCache = nil;
    }
}

//__________________________________________________________________________________
void _Formula::DuplicateReference  (const _Formula* f)
{
    for (unsigned long i=0; i<f->theFormula.lLength; i++) {
        _Operation *theO = ((_Operation**)f->theFormula.lData)[i];
        if (theO->GetAVariable()==-2) {
            theFormula.AppendNewInstance(new _Operation ((_PMathObj)LocateVar (-theO->GetNoTerms()-1)->Compute()->makeDynamic()));
        } else {
            theFormula&& theO;
        }
    }
}

//__________________________________________________________________________________
BaseRef _Formula::makeDynamic (void)
{
    _Formula * res = new _Formula;
    res->Duplicate((BaseRef)this);
    return (BaseRef)res;
}
//__________________________________________________________________________________

_Formula::~_Formula (void)
{
  
    Clear();
}

//__________________________________________________________________________________
void    _Formula::Clear (void)
{
    if (theTree) {
        theTree->delete_tree();
        delete theTree;
    }
    theTree = nil;
    if (resultCache) {
        DeleteObject (resultCache);
    }

    theFormula.Clear();
    if (recursion_calls) {
      delete (recursion_calls);
    }
  
//  theStack.Clear();
}

//__________________________________________________________________________________
BaseRef _Formula::toStr (_List* matchedNames, bool dropTree)
{
    ConvertToTree(false);

    _String * result = new _String(16UL,true);

    long          savepd = printDigits;
    printDigits          = 0;

    if (theTree) { // there is something to do
        internalToStr (*result, theTree, -1, matchedNames);
    } else {
        if (theFormula.lLength) {
            (*result) << "RPN:";
            internalToStr (*result,nil,0,nil,(_Operation*)(theFormula(0)));
            for (unsigned long k=1; k<theFormula.lLength; k++) {
                (*result)<<'|';
                internalToStr (*result,nil,0,nil,(_Operation*)(theFormula(k)));
            }
        }
    }

    printDigits = savepd;
    result->Finalize ();
    if (theTree && dropTree) {
        theTree->delete_tree();
        delete theTree;
        theTree = nil;
    }
    return result;
}
//__________________________________________________________________________________
node<long>* _Formula::DuplicateFormula (node<long>* src, _Formula& tgt)
{
    node<long>* resNode = new node<long>;
    checkPointer (resNode);

    tgt.theFormula && (_Operation*) theFormula (src->in_object);

    resNode->in_object = tgt.theFormula.lLength-1;

    for (long k=1; k<=src->get_num_nodes(); k++) {
        resNode->add_node (*DuplicateFormula (src->go_down (k), tgt));
    }

    return     resNode;
}

//__________________________________________________________________________________
_Formula* _Formula::Differentiate (_String varName, bool bail, bool convert_from_tree)
{
    long          varID = LocateVarByName (varName);

    if (varID<0) {
        return new _Formula (new _Constant (0.0));
    }

    varID = variableNames.GetXtra (varID);

    _Formula*     res = new _Formula ();
 
     ConvertToTree    ();

    _SimpleList  varRefs,
                 dydx;


    _AVLList al (&varRefs);
    ScanFForVariables (al, true, true, true);
    al.ReorderList ();
  
    for (unsigned long k=0UL; k < varRefs.lLength; k++) {
        _Variable* thisVar = LocateVar (varRefs.lData[k]);
        _Formula * dYdX;
      
        if (thisVar->IsIndependent()) {
            dYdX = new _Formula ((thisVar->GetName()->Equal (&varName))?new _Constant (1.0):new _Constant (0.0));
            dYdX->ConvertToTree();
            dydx << (long)dYdX;
        } else {
            dYdX = thisVar->varFormula->Differentiate (varName, bail, false);
          
            if (dYdX->theFormula.lLength == 0) {
                delete (dYdX);
                return res;
            }
            dydx << (long)dYdX;
        }
      }

    SortLists             (&varRefs, &dydx);
    node<long>*           dTree = nil;

    if (!(dTree = InternalDifferentiate (theTree, varID, varRefs, dydx, *res))) {
        for (unsigned long k=0UL; k<dydx.lLength; k++) {
            delete ((_Formula*)dydx.lData[k]);
        }

        if (bail) {
            WarnError    (_String ("Differentiation of ") & _String((_String*)toStr()) & " failed.");
            res->Clear();
            return       res;
        } else {
            delete res;
            return nil;
        }
    }

    for (unsigned long k=0UL; k<dydx.lLength; k++) {
        delete ((_Formula*)dydx.lData[k]);
    }

    res->theFormula.AppendNewInstance (new _Operation(new _Constant (0.0))) ;
    res->theTree         = dTree;
  
    // consistency check
  
    /*_SimpleList preOrderNodes;
    node_iterator<long> iterator (dTree, _HY_TREE_TRAVERSAL_PREORDER);
    
    while (node <long> * ni = iterator.Next()) {
      printf ("Differentiate %x\n", ni);
      if (preOrderNodes.Find ((long)ni) >= 0) {
        WarnError ("Tree construction error!");
        return;
      }
      preOrderNodes << (long)ni;
    }
  
    StringToConsole (_String("\nNode count:") & (long) preOrderNodes.countitems());
    NLToConsole();*/
  
    res->InternalSimplify (dTree);
  
    if (convert_from_tree)
      res->ConvertFromTree  ();
    return res;

}


//__________________________________________________________________________________
bool _Formula::InternalSimplify (node<long>* startNode)
// returns true if the subexpression at
// and below startnode is constant
{
    long        numChildren = startNode->get_num_nodes();

    _Operation* op = GetIthTerm(startNode->get_data());
  
    if  (numChildren == 0) {
      return !op->IsAVariable();
    }

  
    bool        isConstant  = true,
                firstConst  = true,
                secondConst = (numChildren>1L);
  
    long        prune_this_child = -1;

    _Parameter  theVal      = 0.0;
    _PMathObj   newVal      = nil;


    //printf ("InternalSimplify %x\n", startNode);
  
    for  (unsigned long k=1UL; k<=numChildren; k++) {
        if (k==1UL) {
            firstConst = InternalSimplify (startNode->go_down(k));
        } else if (k==2UL) {
            secondConst = InternalSimplify (startNode->go_down(k));
        } else {
          if (!InternalSimplify (startNode->go_down(k))) {
            isConstant = false;
          }
        }
    }

    isConstant = isConstant && firstConst && (numChildren==1 || secondConst);

    if (op->opCode > HY_OP_CODE_NONE) {
        if (isConstant) { // this executes the subxpression starting at the current node
            _Stack scrap;
            for  (unsigned long k=1UL; k<=numChildren; k++) {
                ((_Operation*)theFormula (startNode->go_down(k)->get_data()))->Execute (scrap);
            }
            op->Execute (scrap);
            newVal = (_PMathObj)scrap.Pop();//->makeDynamic();
        } else {
            if (firstConst||secondConst) {
              
                _PMathObj constant_value = ((_Operation*)theFormula (startNode->go_down(firstConst?1:2)->get_data()))->GetANumber();
              
              
                if (constant_value->ObjectClass() == NUMBER) {
                  theVal  = constant_value->Value();

                  switch (op->opCode) {
                      case HY_OP_CODE_MUL: { // *
                          if (CheckEqual (theVal,0.0)) { // *0 => 0
                                                         //printf ("*0\n");
                              newVal = new _Constant (0.0);
                              break;
                          }
                          if (CheckEqual (theVal,1.0)) { // x*1 => x
                                                         //printf ("*1\n");
                              prune_this_child = firstConst?1:2;
                              break;
                          }
                      }
                      break;

                      case HY_OP_CODE_ADD: { // +
                          if (CheckEqual (theVal,0.0)) { // x+0 => x
                                                         //printf ("+0\n");
                              prune_this_child = firstConst?1:2;
                          }
                          break;
                      }

                      case HY_OP_CODE_SUB: { // x-0 => x
                                             // 0-x => -x
                        
                          if (CheckEqual (theVal,0.0)) {
                              //printf ("-0\n");
                             prune_this_child = firstConst? -2 : 2;
                          }
                          break;
                      }

                      case HY_OP_CODE_DIV: { // /
                          if (firstConst&&CheckEqual (theVal,0.0)) { // 0/x => 0
                              newVal = new _Constant (0.0);
                              //printf ("0/\n");

                              break;
                          }
                          if (secondConst&&CheckEqual (theVal,1.0)) { // x/1 => x
                                                                      //printf ("/1\n");
                              prune_this_child = 2;
                              break;
                          }
                      }
                      break;

                      case HY_OP_CODE_POWER: { // ^
                          if (firstConst&&CheckEqual (theVal,1.0)) { // 1^x => 1
                                                                     //printf ("1^\n");
                              newVal = new _Constant (1.0);
                              break;
                          }
                          if (secondConst&&CheckEqual (theVal,1.0)) { // x^1 => ?
                                                                      //printf ("^1\n");
                              prune_this_child = 1;
                              break;
                          }
                      }
                      break;
                  }
                }
            }
        }
    }

    if (newVal) {
        for  (int k=1; k <= numChildren; k++) {
            startNode->go_down(k)->delete_tree(true);
        }
        startNode->kill_all_nodes();
        startNode->in_object = theFormula.lLength;
        theFormula < (new _Operation(newVal));
    } else {

      if (prune_this_child !=- 1L) {
          if (prune_this_child > 0L) {


              startNode->go_down(prune_this_child)->delete_tree(true);
              startNode->kill_node (prune_this_child);
              node <long>*    replaceWith = startNode->go_down(1);

              startNode->kill_all_nodes();

              for (unsigned long k=1; k<=replaceWith->get_num_nodes(); k++) {
                  startNode->add_node (*replaceWith->go_down(k));
              }
              startNode->in_object = replaceWith->in_object;
              delete (replaceWith);

          } else { // 0-? => -?
              startNode->go_down(1)->delete_tree(true);
              startNode->kill_node(1);
              op->SetTerms(1);
             //startNode->kill_node(1);
          }
      }
    }
    return isConstant;
}


//__________________________________________________________________________________
void _Formula::internalToStr (_String& result, node<long>* currentNode, char opLevel, _List* matchNames, _Operation* thisNodeOperation)
{
    if (!thisNodeOperation) {
        thisNodeOperation = (_Operation*)theFormula (currentNode->get_data());
    }

    // decide what to do about this operation

    if (thisNodeOperation->IsAVariable(false))
        // this operation is just a variable - add ident to string and return
    {
        if (subNumericValues) {
            if (subNumericValues == 2) {
                _Variable* theV = LocateVar(thisNodeOperation->GetAVariable());
                if  (_x_&&(theV->GetAVariable()==_x_->GetAVariable())) {
                    result << _x_->GetName();
                    return;
                }
            }

            _Variable *thisVariable = LocateVar(thisNodeOperation->GetAVariable());
            _PMathObj subThisValue = thisVariable->Compute();

            if (subThisValue->ObjectClass () == NUMBER) {
                if (subNumericValues == 3) {
                    result << LocateVar(thisNodeOperation->GetAVariable())->GetName();
                    result << '[';
                    result.AppendNewInstance(new _String (subThisValue->Value()));
                    result << ':';
                    result.AppendNewInstance(new _String (thisVariable->GetLowerBound()));
                    result << '-';
                    result.AppendNewInstance(new _String (thisVariable->GetUpperBound()));
                    result << ']';

                } else {
                    result.AppendNewInstance(new _String (subThisValue->Value()));
                }
            } else if (subThisValue->ObjectClass () == STRING) {
                result.AppendNewInstance((_String*)subThisValue->toStr());
            } else {
                result << LocateVar(thisNodeOperation->GetAVariable())->GetName();
            }
        } else {
            long variableIDX = thisNodeOperation->GetAVariable();
            if (variableIDX>=0) {
                _String * vName = LocateVar(variableIDX)->GetName();
                if (matchNames) {
                    _List * p1 = (_List*)(*matchNames)(0),
                            * p2 = (_List*)(*matchNames)(1);

                    long  f = p1->FindObject (vName);

                    if (f<0) {
                        result<<vName;
                    } else {
                        result<<(_String*)(*p2)(f);
                    }
                } else {
                    result<<vName;
                }
            }
        }
        return;
    }

    long nOps = thisNodeOperation->GetNoTerms();
    if (nOps>0)
        // a built-in operation or a function call
    {
        // check if it's a built-in binary operation
        _String  opString (thisNodeOperation->GetCode());
        long f = BinOps.Find ((opString.sLength>1)?(opString.sData[0]*256+opString.sData[1]):opString.sData[0]);
        if (f!=-1)
            // indeed - a binary operation is what we have. check if need to wrap the return in parentheses
        {
            if (!currentNode || currentNode->get_num_nodes()==2 ) {
                char tOpLevel  = opPrecedence(f),
                     tOpLevel2 = tOpLevel;

                if (associativeOps.Find(f)<0) {
                    tOpLevel2 ++;
                }
                if (opLevel>=0) { // need to worry about op's precedence
                    {
                        bool parens = opPrecedence(f)<opLevel;

                        if (parens&&currentNode) { // put parentheses around the return of this expression
                            result<<'(';
                            internalToStr (result, currentNode->go_down(1),tOpLevel,matchNames);
                            result<<&thisNodeOperation->GetCode();
                            internalToStr (result, currentNode->go_down(2),tOpLevel2,matchNames);
                            result<<')';
                            return;
                        }
                    }
                }
                if (currentNode) {
                    internalToStr (result, currentNode->go_down(1),tOpLevel,matchNames);
                }
                result<<&thisNodeOperation->GetCode();
                if (currentNode) {
                    internalToStr (result, currentNode->go_down(2),tOpLevel2,matchNames);
                }
                return;
            } else { // mixed binary-unary operation
                result<<&thisNodeOperation->GetCode();
                if (currentNode) {
                    result<<'(';
                    internalToStr (result, currentNode->go_down(1),opPrecedence(f),matchNames);
                    result<<')';
                }
                return;
            }
        } else {
            _String mac ("MAccess");
            if (!thisNodeOperation->GetCode().Equal(&mac)) {
                result<<&thisNodeOperation->GetCode();
                if (currentNode) {
                    result<<'(';
                    for (long k=1; k<=nOps; k++) {
                        if (k>1) {
                            result<<',';
                        }
                        internalToStr (result, currentNode->go_down(k),-1,matchNames);
                    }
                    result<<')';
                }
            } else { // matrix element access - treat specially
                if (currentNode) {
                    internalToStr (result, currentNode->go_down(1),-1,matchNames);
                    for (long k=2; k<=nOps; k++) {
                        result<<'[';
                        internalToStr (result, currentNode->go_down(k),-1,matchNames);
                        result<<']';
                    }
                }
            }
        }
        return;
    }
    if (nOps<0)
        // a user-defined function
    {
        long func_id = thisNodeOperation->UserFunctionID();
        result<< & GetBFFunctionNameByIndex(func_id);
        if (currentNode) {
            result<<'(';
            long argument_count = GetBFFunctionArgumentCount(func_id);
            for (long k=1L; k<=argument_count; k++) {
                if(k>1L) {
                    result<<',';
                }
                internalToStr (result, currentNode->go_down(k),-1,matchNames);
            }
            result<<')';
        }
        return;
    }
    _PMathObj opValue = thisNodeOperation->GetANumber();
    _String* conv = (_String*)opValue->toStr();
    if (opValue->ObjectClass()==STRING) {
        result<<'"';
        result<<conv;
        result<<'"';
    } else {
        if (opValue->ObjectClass() == NUMBER && opValue->Value() < 0.0) {
            result<<'(';
            result<<conv;
            result<<')';
        } else {
            result<<conv;
        }
    }
    DeleteObject(conv);
}
//__________________________________________________________________________________
bool     _Formula::IsEmpty(void) const {
// is there anything in the formula
    return bool(!theFormula.lLength);
}

//__________________________________________________________________________________
_Parameter   _Formula::Newton(_Formula& derivative, _Variable* unknown, _Parameter targetValue, _Parameter left, _Parameter right)
// find a root of the formulaic expression, using Newton's method, given the derivative and a bracketed root.
// will alternate between bisections and Newton iterations based on what is fatser
{
  // check that there is indeed a sign change on the interval
  _Parameter    func_left, func_right, // keeps track of function values in the current interval, [left,right]
  root_guess, func_root_guess,
  lastCorrection = 100.,
  newCorrection;
  
  
  unknown->SetValue(left);
  func_left = Compute()->Value()-targetValue;
  if (func_left==0.0) {
    return left;
  }
  unknown->SetValue(right);
  func_right = Compute()->Value()-targetValue;
  if (func_right==0.0) {
    return right;
  }
  
  if (func_left*func_right>0.0) { // bracket fail
    ReportWarning (_String((_String*)toStr())&"="&_String(targetValue)&" has no (or multiple) roots in ["&_String(left)&",Inf)");
    return    left;
  }
  // else all is good we can start the machine
  bool useNewton = false;
  
  root_guess = (right+left) * .5;
  
  for (unsigned long iterCount  = 0L; iterCount < 200UL; iterCount++) {
    if ( fabs(right-left)/MAX(left,right) > machineEps*10. ) { // stuff to do
      
      unknown->SetValue(root_guess);
      func_root_guess = Compute()->Value()-targetValue;
      if (func_root_guess == 0.) {
        return root_guess;
      }
      // get the correction term from the derivative
      _Parameter df_dx = derivative.Compute()->Value(),
      adjusted_root_guess;
      
      useNewton = true;
      if (df_dx==0.0) {
        useNewton = false;
      } else {
        newCorrection = -func_root_guess/df_dx;
        
        if (fabs(newCorrection/func_root_guess)<machineEps*2. || fabs(newCorrection)<machineEps*2.) { // correction too small - the root has been found
            return root_guess;
        }
        
        if (fabs(newCorrection/lastCorrection)>4.) { // Newton correction too large - revert to bisection
          useNewton = false;
        }
        
        adjusted_root_guess = root_guess +newCorrection;
        if (adjusted_root_guess<=left || adjusted_root_guess >=right) {
          useNewton = false;
        } else {
          lastCorrection = newCorrection;
        }
      }
      
      if (useNewton) {
        root_guess = adjusted_root_guess;
      } else {
        if (func_root_guess==0.0) {
          return root_guess;
        }
        if (func_root_guess*func_left > 0.0) { // move to the right half
          func_left   = func_root_guess;
          left  = root_guess;
        } else {
          right = root_guess;
        }
        root_guess = (right+left) * .5;
      }
      
    } else {
      break;
    }
  }
  return root_guess;
}


//__________________________________________________________________________________
_Parameter   _Formula::Brent(_Variable* unknown, _Parameter a, _Parameter b, _Parameter tol, _List* storeEvals, _Parameter rhs)
// find a root of the formulaic expression, using Brent's method and a bracketed root.
{
    // check that there is indeed a sign change on the interval

    _Parameter  fa = 0.0,fb = 0.0,fc,d = b-a,e = b-a ,min1,min2,xm,p,q,r,s,tol1,
                c = b;

    min1 = unknown->GetLowerBound();
    min2 = unknown->GetUpperBound();

    long        it = 0;

    if (a>b) { // autobracket to the left
        unknown->SetValue(b);
        fb = Compute()->Value();
        if (storeEvals) {
            storeEvals->AppendNewInstance(new _Constant (b));
            storeEvals->AppendNewInstance(new _Constant (fb));
        }

        if (b<0.00001 && b>-0.00001) {
            a = b-0.0001;
        } else {
            a = b-fabs(b)*0.1;
        }

        if (a<min1) {
            a = min1+0.5*(b-min1);
        }

        unknown->SetValue(a);
        fa = Compute()->Value()-rhs;
        if (storeEvals) {
          storeEvals->AppendNewInstance(new _Constant (a));
          storeEvals->AppendNewInstance(new _Constant (fa));
        }

        for (long k=0; k<50; k++) {
            if (fb*fa<0.0) {
                break;
            }

            d  = (b-a)*GOLDEN_RATIO;
            b  = a;
            a -= d;

            if (a<min1) {
                if (b>min1) {
                    a = min1;
                } else {
                    break;
                }
            }

            fb = fa;

            unknown->SetValue(a);
            fa = Compute()->Value()-rhs;
            if (storeEvals) {
              storeEvals->AppendNewInstance(new _Constant (a));
              storeEvals->AppendNewInstance(new _Constant (fa));
            }
        }
    } else if (CheckEqual (a,b)) { // autobracket to the right
        unknown->SetValue(a);
        fa = Compute()->Value()-rhs;

        if (storeEvals) {
          storeEvals->AppendNewInstance(new _Constant (a));
          storeEvals->AppendNewInstance(new _Constant (fa));
        }
        a = b;

        if ((b<0.00001)&&(b>-0.00001)) {
            b = b+0.0001;
        } else {
            b = b+fabs(b)*0.1;
        }

        if (b>min2) {
            b = a+0.5*(min2-a);
        }

        unknown->SetValue(b);
        fb = Compute()->Value()-rhs;
      
        if (storeEvals) {
          storeEvals->AppendNewInstance(new _Constant (b));
          storeEvals->AppendNewInstance(new _Constant (fb));
        }

        for (long k=0; k<50; k++) {
            if (fb*fa<0.0) {
                break;
            }

            d  = (b-a)*GOLDEN_RATIO;
            a  = b;
            b += d;

            if (b>min2) {
                if (a<min2) {
                    b = min2;
                } else {
                    break;
                }
            }

            fa = fb;

            unknown->SetValue(b);
            fb = Compute()->Value()-rhs;
            if (storeEvals) {
              storeEvals->AppendNewInstance(new _Constant (b));
              storeEvals->AppendNewInstance(new _Constant (fb));
            }
        }
    }


    if (fa == 0.0) {
        unknown->SetValue(a);
        fa = Compute()->Value()-rhs;
        if (storeEvals) {
          storeEvals->AppendNewInstance(new _Constant (a));
          storeEvals->AppendNewInstance(new _Constant (fa));
        }
        if (fa == 0.0) {
            return a;
        }
    }

    if (fb == 0.0) {
        unknown->SetValue(b);
        fb = Compute()->Value()-rhs;
        if (storeEvals) {
          storeEvals->AppendNewInstance(new _Constant (b));
          storeEvals->AppendNewInstance(new _Constant (fb));
        }
        if (fb == 0.0) {
            return b;
        }
    }

    if (fa*fb<0.0) {
        fc = fb;

        for (it = 0; it < MAX_BRENT_ITERATES; it++) {
            if (fb*fc>0.0) {
                fc = fa;
                c  = a;
                e = d = b-a;
            }

            if (fabs (fc) < fabs (fb)) {
                a     = b;
                b     = c;
                c     = a;
                fa    = fb;
                fb    = fc;
                fc    = fa;
            }

            tol1 = 2.*fabs(b)*machineEps+.5*tol;

            xm = .5*(c-b);

            if (fabs(xm)<=tol1 || fb == 0.0) {
                return b;
            }

            if (fabs(e)>=tol1 && fabs (fa) > fabs (fb)) {
                s = fb/fa;
                if (a==c) {
                    p = 2.*xm*s;
                    q = 1.-s;
                } else {
                    q = fa/fc;
                    r = fb/fc;
                    p = s*(2.*xm*q*(q-r)-(b-a)*(r-1.0));
                    q = (q-1.)*(r-1.)*(s-1.);
                }
                if (p>0.0) {
                    q = -q;
                }

                if (p<0.0) {
                    p = -p;
                }

                min1 = 3.*xm*q-fabs(tol1*q);
                min2 = fabs (e*q);
                if (2.*p < (min1<min2?min1:min2)) {
                    e = d;
                    d = p/q;
                } else {
                    d = xm;
                    e = d;
                }
            } else {
                d = xm;
                e = d;
            }
            a = b;
            fa = fb;
            if (fabs(d)>tol1) {
                b+=d;
            } else {
                b+=fabs(tol1)*(1.*2.*(xm>=.0));
            }

            unknown->SetValue(b);
            fb = Compute()->Value()-rhs;
            if (storeEvals) {
              storeEvals->AppendNewInstance(new _Constant (b));
              storeEvals->AppendNewInstance(new _Constant (fb));
            }
        }
    }

    subNumericValues = 2;
  
    /*for (long i = 0; i < theFormula.lLength; i++) {
      _Operation *op_i = GetIthTerm(i);
      printf ("%ld: %s\n", i+1, _String((_String*)op_i->toStr()).sData);
    }*/
  
    _String msg ((_String*)toStr());
    subNumericValues = 0;
    msg = msg & "=" & rhs;
    if (it < MAX_BRENT_ITERATES) {
        msg =   msg & " has no (or multiple) roots in ["&_String(a)&","&_String(b)&"]";
    } else {
        msg =   msg & " failed to converge to sufficient precision in " & MAX_BRENT_ITERATES &" iterates.";
    }

    ReportWarning (msg);
    return    b;
}
//__________________________________________________________________________________
_Parameter   _Formula::Newton(_Formula& derivative, _Parameter targetValue, _Parameter left, _Parameter max_right, _Variable* unknown)
// given a monotone function and a left bracket bound, found the right bracket bound and solve
{
    // check that there is indeed a sign change on the interval
    _Constant   dummy;
    dummy.SetValue (left);
    unknown->SetValue(&dummy);
    _Parameter  t1 = Compute()->Value(), right = left, t2, step = 1.0;

    if (max_right-left < step * 100) {
        step = (max_right-left)/100;
    }
    if  (step==0.0) {
        return left;
    }
    do {
        right += step;
        if (right>max_right) { // function doesn't seem to have a root
            subNumericValues = 2;
            _String *s = (_String*)toStr();
            subNumericValues = 0;
            _String msg = *s&"="&_String(targetValue)&" has no (or multiple) roots in ["&_String(left)&","&_String(right)&")";
            ReportWarning (msg);
            DeleteObject (s);
            return    left;
        }
        dummy.SetValue (right);
        unknown->SetValue(&dummy);
        t2 = Compute()->Value();
        step*=2;
        if (right+step>max_right)
            if (right<max_right) {
                step = max_right-right;
            }
    } while ((targetValue-t1)*(targetValue-t2)>0);
    return Newton (derivative, unknown, targetValue, left, right);

}

//__________________________________________________________________________________
_Parameter   _Formula::Newton(_Variable* unknown, _Parameter targetValue, _Parameter x_min, _Parameter left, _Parameter right)
{
    // check that there is indeed a sign change on the interval
    _Constant   dummy;
    _Parameter  t1,t2,t3,t4,t5,lastCorrection = 100, newCorrection;
    _String     msg;
    t1 =Integral(unknown, x_min, left)-targetValue;
    if (t1==0.0) {
        return left;
    }
    t2 = t1+Integral(unknown, left, right);
    if (t2==0) {
        return right;
    }
    if (t1*t2>0.0) {
        subNumericValues = 2;
        _String *s = (_String*)toStr();
        subNumericValues = 0;
        _String msg = *s&"="&_String(targetValue)&" has no (or multiple) roots in ["&_String(left)&","
                      &_String (right)&"]";
        ReportWarning (msg);
        DeleteObject (s);
        return    left;
    }
    // else all is good we can start the machine
    bool useNewton = false;
  
    t3 = (right + left) * 0.5;
  
    while (right-left>1e-6) { // stuff to do
        if (!useNewton) {
            t3 = (right+left) * 0.5;
        }
        dummy.SetValue(t3);
        unknown->SetValue(&dummy);
        t4 = Integral(unknown, x_min, t3)-targetValue;
        // get the correction term from the derivative
        dummy.SetValue(t3);
        unknown->SetValue(&dummy);
        t5 = Compute()->Value();
        useNewton = true;
        if (t5==0.0) {
            useNewton = false;
        } else {
            newCorrection = -t4/t5;
            if (fabs(newCorrection)<1e-5) { // correction too small - the root has been found
                return t3;
            }
            if (fabs(newCorrection/lastCorrection)>4) { // Newton correction too large - revert to bisection
                useNewton = false;
            }
            t5 = t3+newCorrection;
            if ((t5<=left)||(t5>=right)) {
                useNewton = false;
            } else {
                lastCorrection = newCorrection;
            }
        }
        if (useNewton) {
            t3 = t5;
        } else {
            t4 = Integral(unknown, x_min, t3)-targetValue;
            if (t4==0.0) {
                return t3;
            }
            if (t4*t1 >0) {
                t1 = t4;
                left = t3;
            } else {
                right = t3;
            }
        }

    }
    return t3;
}

//__________________________________________________________________________________
_Parameter   _Formula::Newton( _Variable* unknown, _Parameter targetValue,_Parameter x_min, _Parameter left)
// given a monotone function and a left bracket bound, found the right bracket bound and solve
{
    // check that there is indeed a sign change on the interval
    _Parameter  t1 = Integral(unknown, x_min, left), right = left, t2, step = 1.0;
    do {
        right += step;
        t2 = Integral(unknown, right-step, right);
        step*=2;
        if (right>=1e10) { // function doesn't seem to have a root
            subNumericValues = 2;
            _String *s = (_String*)toStr();
            subNumericValues = 0;
            _String msg = *s&"="&_String(targetValue)&" has no (or multiple) roots in ["&_String(left)&",Inf)";
            WarnError (msg);
            DeleteObject (s);
            return    0.0;
        }
    } while ((targetValue-t1)*(targetValue-t2-t1)>=0);
    return Newton ( unknown, targetValue,x_min, left, right);

}

_Parameter   _Formula::Integral(_Variable* dx, _Parameter left, _Parameter right, bool infinite)
// uses Romberg's intergation method
{
    if (infinite) { // tweak "right" here
        _Parameter value = 1.0, step = 1.0, right1=-1;
        right = left;
        while (value>machineEps) {
            right+=step;
            _Constant dummy (right);
            dx->SetValue(&dummy);
            value = fabs(Compute()->Value());
            if ((value<0.001)&&(right1<0)) {
                right1 = right;
            }
            step *= 2;
            if (step>100000) { // integrand decreasing too slowly
                _String msg, *s = (_String*)toStr();
                msg = *s & " decreases too slowly to be integrated over an infinite interval";
                DeleteObject(s);
                WarnError (msg);
                return 0.0;
            }
        }
        if (right1<right-machineEps) {
            return Integral(dx,left,right1,false)+Integral(dx,right1,right,false);
        } else {
            return Integral(dx,left,right1,false);
        }
    }
  
    _Parameter       localPrecisionFactor;
    long             localIntegrationLoops;
  
    checkParameter (intPrecFact,localPrecisionFactor,integrationPrecisionFactor);
    checkParameter (intMaxIter, localIntegrationLoops,maxRombergSteps);

    _Parameter ss,
               dss,
               *s,
               *h;

    s = new _Parameter [(long)maxRombergSteps];
    h = new _Parameter [(long)(maxRombergSteps+1)];

    h[0]=1.0;

    long         interpolateSteps = 5,
                 stackDepth = 0;

    _SimpleList  fvidx,
                 changingVars,
                 idxMap;

    _Parameter   * ic = new _Parameter[interpolateSteps],
    * id = new _Parameter[interpolateSteps];

    _SimpleFormulaDatum
    * stack = nil,
      * vvals = nil;

    checkPointer (ic);
    checkPointer (id);

    if (AmISimple (stackDepth,fvidx)) {
        stack = new _SimpleFormulaDatum [stackDepth];
        vvals = new _SimpleFormulaDatum [fvidx.lLength];
        ConvertToSimple (fvidx);
        stackDepth = dx->GetAVariable();
        for (long vi = 0; vi < fvidx.lLength; vi++) {
            _Variable * checkvar = LocateVar (fvidx.lData[vi]);
            if (checkvar->CheckFForDependence (stackDepth,true)) {
                changingVars << fvidx.lData[vi];
                idxMap << vi;
            }
            vvals[vi].value = checkvar->Compute()->Value();
        }
        changingVars.InsertElement ((BaseRef)stackDepth,0,false,false);
        idxMap.InsertElement ((BaseRef)fvidx.Find(stackDepth),0,false,false);
    } else {
        stackDepth = -1;
    }

    for (long j=0; j<(long)maxRombergSteps; j++) {
        if (stackDepth>=0) {
            s[j] = TrapezoidLevelKSimple(*this, dx, left, right, j+1, stack, vvals,changingVars,idxMap);
        } else {
            s[j] = TrapezoidLevelK(*this, dx, left, right, j+1);
        }
        if (j>=4) {
            ss = InterpolateValue(&h[j-4],&s[j-4],interpolateSteps,ic,id,0.0, dss);
            if (fabs(dss)<= integrationPrecisionFactor * fabs(ss)) {
                delete s;
                delete h;
                delete ic;
                delete id;
                if (stackDepth>=0) {
                    ConvertFromSimple(fvidx);
                    delete (stack);
                    delete (vvals);
                }
                return ss;
            }
        }
        h[j+1] =  h[j]/9.0;
    }

    if (stackDepth>=0) {
        ConvertFromSimple(fvidx);
        delete [] stack;
        delete [] vvals;
    }
    _String *str = (_String*)toStr(),
             msg = _String("Integral of ")&*str & " over ["&_String(left)&","&_String(right)&"] converges slowly, loss of precision may occur. Change either INTEGRATION_PRECISION_FACTOR or INTEGRATION_MAX_ITERATES";

    DeleteObject (str);
    ReportWarning (msg);

    delete [] s;
    delete [] h;
    delete [] ic;
    delete [] id;
    return ss;
}

//__________________________________________________________________________________
_Parameter   _Formula::MeanIntegral(_Variable* dx, _Parameter left, _Parameter right, bool infinite)
{
    _Formula newF;
    _String tim ("*");
    _Operation times (tim,2);
    _Operation term (true, *(dx->GetName()));
    newF.Duplicate((BaseRef)this);
    newF.theFormula&& (&term);
    newF.theFormula&& (& times);
    return newF.Integral (dx,left,right, infinite);
}

//__________________________________________________________________________________
long     _Formula::NumberOperations(void) const
// number of operations in the formula
{
    return theFormula.lLength;
}

//__________________________________________________________________________________

long      _Formula::ExtractMatrixExpArguments (_List* storage) {
    long count = 0;
    
    if (resultCache && resultCache->lLength) {
        long cacheID     = 0;     
        bool cacheUpdated = false; 
            // whether or not a cached result was used 
            

        for (unsigned long i=0; i<theFormula.lLength; i++) {
            _Operation* thisOp ((_Operation*)(((BaseRef**)theFormula.lData)[i]));
            if (i < theFormula.lLength-1) {
                _Operation* nextOp  ((_Operation*)(((BaseRef**)theFormula.lData)[i+1]));

                if (! cacheUpdated && nextOp->CanResultsBeCached(thisOp)) {
                    _Stack temp;
                    thisOp->Execute (temp);

                    _Matrix *currentArg  = (_Matrix*)temp.Pop(true),
                            *cachedArg   = (_Matrix*)((_PMathObj)(*resultCache)(cacheID)),
                            *diff        = nil;

                    if (cachedArg->ObjectClass() == MATRIX) {
                        diff =  (_Matrix*)cachedArg->SubObj(currentArg);
                    }

                    if (diff && diff->MaxElement() <= 1e-12) {
                        cacheID += 2;
                        i ++;
                    } else {
                        cacheUpdated = true;
                        cacheID++;
                        if (nextOp->CanResultsBeCached(thisOp, true)) {
                            storage->AppendNewInstance(currentArg);
                            count ++;
                        }
                    }
                    DeleteObject (diff);
                    continue;
                }
                if (cacheUpdated) {
                    cacheID++;
                    cacheUpdated = false;
                }            
            }
        }
    }
    
    return count;
}

//__________________________________________________________________________________
_Variable * _Formula::Dereference (bool ignore_context, _hyExecutionContext* theContext) {
    _Variable * result = nil;
    _PMathObj computedValue = Compute (0, theContext->GetContext(), nil, theContext->GetErrorBuffer());
    if (computedValue && computedValue->ObjectClass() == STRING) {
        result =  (_Variable*)((_FString*)computedValue)->Dereference(ignore_context, theContext, true);
    }
    
    if (!result) {
        theContext->ReportError( (_String ("Failed to dereference '") & _String ((_String*)toStr()) & "' in the " & (ignore_context ? "global" : "local") & " context"));
    }
    
    return result;
}

  //unsigned long ticker = 0UL;

//__________________________________________________________________________________
_PMathObj _Formula::Compute (long startAt, _VariableContainer const * nameSpace, _List* additionalCacheArguments, _String* errMsg, long valid_type)
// compute the value of the formula
{
    _Stack * scrap_here;
    if (theFormula.lLength == 0) {
        theStack.theStack.Clear();
        theStack.Push (new _MathObject, false);
        scrap_here = &theStack;
    } else {
        bool wellDone = true;
      
      
        if (call_count++) {
          scrap_here = new _Stack;
        } else {
          scrap_here = &theStack;
          if (startAt == 0) {
              theStack.Reset();
          }
        }
      
      
      
        /*ticker++;
        if (ticker >= 1462440) {
          printf ("\n_Formula::Compute (%x, %d)  %ld terms, stack depth %ld\n", this, ticker, theFormula.lLength, theStack.theStack.lLength);
        }*/
      
        if (startAt == 0 && resultCache && resultCache->lLength) {
            long cacheID     = 0;     
                // where in the cache are we currently looking
            bool cacheUpdated = false; 
                // whether or not a cached result was used 

            for (unsigned long i=0; i<theFormula.lLength; i++) {
                _Operation* thisOp ((_Operation*)(((BaseRef**)theFormula.lData)[i]));
                if (i < theFormula.lLength-1) {
                    _Operation* nextOp  ((_Operation*)(((BaseRef**)theFormula.lData)[i+1]));

                    if (! cacheUpdated && nextOp->CanResultsBeCached(thisOp)) {
                        if (!thisOp->Execute(*scrap_here,nameSpace, errMsg)) {
                            wellDone = false;
                            break;
                        }

                        _Matrix *currentArg  = (_Matrix*)scrap_here->Pop(false),
                                *cachedArg   = (_Matrix*)((_PMathObj)(*resultCache)(cacheID)),
                                *diff        = nil;

                        if (cachedArg->ObjectClass() == MATRIX) {
                            diff =  (_Matrix*)cachedArg->SubObj(currentArg);
                        }
                        
                        bool    no_difference = diff && diff->MaxElement() <= 1e-12;

                        if (no_difference || (additionalCacheArguments && additionalCacheArguments->lLength && nextOp->CanResultsBeCached(thisOp,true))) {
                            DeleteObject  (scrap_here->Pop  ());
                            if (no_difference) {
                                scrap_here->Push ((_PMathObj)(*resultCache)(cacheID+1));
                            } else {
 
                                scrap_here->Push ((_PMathObj)additionalCacheArguments->GetItem (0));
                                resultCache->Replace(cacheID,scrap_here->Pop(false),true);
                                resultCache->Replace(cacheID+1,(_PMathObj)additionalCacheArguments->GetItem (0),false);
                                additionalCacheArguments->Delete (0, false);                              
                                //printf ("_Formula::Compute additional arguments %ld\n", additionalCacheArguments->lLength);
                           }
                            cacheID += 2;
                            i ++;
                            //printf ("Used formula cache %s\n", _String((_String*)nextOp->toStr()).sData);
                        } else {
                            cacheUpdated = true;
                            resultCache->Replace(cacheID++,scrap_here->Pop(false),true);
                            //printf ("Updated formula cache %s\n", _String((_String*)nextOp->toStr()).sData);
                       }
                       DeleteObject (diff);
                       continue;
                    }
                }
                if (!thisOp->Execute(*scrap_here,nameSpace, errMsg)) { // does this always get executed?
                    wellDone = false;
                    break;
                }
                if (cacheUpdated) {
                    resultCache->Replace(cacheID++,scrap_here->Pop(false),true);
                    cacheUpdated = false;
                }
            }
        } else {
          
            for (unsigned long i=startAt; i<theFormula.lLength; i++) {
                  _Operation * this_step =((_Operation*)(((BaseRef**)theFormula.lData)[i]));
              
              
                  /*
                   if (ticker >= 1462435) {
                    printf ("[FORMULA %x] Step %ld, Stack %ld, Op %s (top stack %s)\n", this, i, theStack.theStack.lLength, _String ((_String*)this_step->toStr()).sData, theStack.theStack.lLength ? _String((_String*)theStack.Pop(false)->toStr()).sData : "empty" );
                  } 
                   */
              
                  if (! this_step->Execute(*scrap_here, nameSpace, errMsg)) {
                      wellDone = false;
                      break;
                  }
              
                  /*
                  if (ticker >= 1462435) {
                    printf ("[DONE %x] Step %ld, Stack %ld, top stack %s\n", this, i, theStack.theStack.lLength, theStack.theStack.lLength ? _String((_String*)theStack.Pop(false)->toStr()).sData : "empty" );
                  }
                   */
            }
        }
        if (scrap_here->StackDepth() != 1L || !wellDone) {
            _String errorText = _String ("'") & _String((_String*)toStr()) & _String("' evaluated with errors.");
          
            if (scrap_here->StackDepth() > 1 && wellDone) {
                errorText = errorText & " Unconsumed values on the stack";
                for (long stack_id = scrap_here->StackDepth()-1; stack_id >= 0; stack_id --) {
                  errorText = errorText & "\n[" & (stack_id+1) & "]------------------\n" & (_String*) scrap_here->Pop(false)->toStr();
                }
                errorText & "\n------------------\n";
            }
          
            if (errMsg) {
                *errMsg = *errMsg & errorText;
            }
            else {
                WarnError (errorText);
            }
            scrap_here->Reset();
            scrap_here->Push (new _Constant (0.0), false);
         }
    }

  
    _PMathObj return_value = scrap_here->Pop(false);
  
    if (theFormula.lLength) {
         DeleteObject (recursion_calls);
         if (--call_count) {
          recursion_calls = return_value;
          return_value->AddAReference();
          delete scrap_here;
          
        } else {
          recursion_calls = nil;
        }
    }
    return valid_type == HY_ANY_OBJECT ? return_value : ((return_value->ObjectClass() & valid_type) ? return_value : nil);
}

//__________________________________________________________________________________
bool _Formula::CheckSimpleTerm (_PMathObj thisObj)
{
    if (thisObj) {
        long oc = thisObj->ObjectClass();
        if (oc !=NUMBER) {
            if (oc ==MATRIX) {
                _Matrix * mv = (_Matrix*)thisObj->Compute();
                if (mv->IsIndependent () && !mv->SparseDataStructure()) {
                    return true;
                }
            }
        } else {
            return true;
        }
    }
    return false;
}

//__________________________________________________________________________________
_Formula _Formula::operator+ (const _Formula& operand2) {
    _Formula joint;
    return joint.PatchFormulasTogether (*this,operand2,HY_OP_CODE_ADD);
}

//__________________________________________________________________________________
_Formula _Formula::operator- (const _Formula& operand2) {
    _Formula joint;
    return joint.PatchFormulasTogether (*this,operand2,HY_OP_CODE_SUB);
}

//__________________________________________________________________________________
_Formula _Formula::operator* (const _Formula& operand2) {
    _Formula joint;
    return joint.PatchFormulasTogether (*this,operand2,HY_OP_CODE_MUL);
}

//__________________________________________________________________________________
_Formula _Formula::operator/ (const _Formula& operand2) {
    _Formula joint;
    return joint.PatchFormulasTogether (*this,operand2,HY_OP_CODE_DIV);
}

//__________________________________________________________________________________
_Formula _Formula::operator^ (const _Formula& operand2) {
    _Formula joint;
    return joint.PatchFormulasTogether (*this,operand2,HY_OP_CODE_POWER);
}

//__________________________________________________________________________________
_Formula& _Formula::PatchFormulasTogether (_Formula& target, const _Formula& operand2, const char op_code) {
    target.Clear();
    target.DuplicateReference(this);
    target.DuplicateReference(&operand2);
    target.theFormula.AppendNewInstance(new _Operation (op_code, 2));
    return target;
}


//__________________________________________________________________________________
void _Formula::ConvertMatrixArgumentsToSimpleOrComplexForm (bool makeComplex)
{
  if (makeComplex) {
    if (resultCache) {
      DeleteObject (resultCache);
      resultCache = nil;
    }
  } else {
    if (!resultCache) {
      resultCache = new _List();
      for (int i=1; i<theFormula.lLength; i++) {
        _Operation* thisOp = ((_Operation*)(((BaseRef**)theFormula.lData)[i]));
        if (thisOp->CanResultsBeCached(((_Operation*)(((BaseRef**)theFormula.lData)[i-1])))) {
          resultCache->AppendNewInstance(new _MathObject());
          resultCache->AppendNewInstance(new _MathObject());
        }
      }
    }
  }
  
  for (int i=0; i<theFormula.lLength; i++) {
    _Operation* thisOp = ((_Operation*)(((BaseRef**)theFormula.lData)[i]));
    
    _Matrix   * thisMatrix = nil;
    
    if (thisOp->theNumber) {
      if (thisOp->theNumber->ObjectClass() == MATRIX) {
        thisMatrix = (_Matrix*)thisOp->theNumber;
      }
    } else {
      if (thisOp->theData>-1) {
        _Variable* thisVar = LocateVar (thisOp->theData);
        if (thisVar->ObjectClass() == MATRIX) {
          thisMatrix = (_Matrix*)thisVar->GetValue();
        }
      }
    }
    
    if (thisMatrix) {
      if (makeComplex) {
        thisMatrix->MakeMeGeneral();
      } else {
        thisMatrix->MakeMeSimple();
      }
    }
  }
}

//__________________________________________________________________________________
long _Formula::StackDepth (long from, long to) const {
  _SimpleList::NormalizeCoordinates(from, to, NumberOperations());
  long result = 0L;
  
  for ( long i = from; i <= to; i++) {
     result += GetIthTerm (i)->StackDepth ();
  }
  
  return result;
  
}


//__________________________________________________________________________________
bool _Formula::AmISimple (long& stackDepth, _SimpleList& variableIndex)
{
    if (!theFormula.lLength) {
        return true;
    }

    long locDepth = 0L;

    for (int i=0; i<theFormula.lLength; i++) {
        _Operation* thisOp = ((_Operation*)(((BaseRef**)theFormula.lData)[i]));
        locDepth++;
        if ( thisOp->theData<-2 || thisOp->numberOfTerms<0) {
            if (thisOp->theData < -2 && i == 0) {
              if (variableIndex.Find (thisOp->GetAVariable())==-1) {
                variableIndex<<thisOp->GetAVariable();
              }
              continue;
            }
            return false;
        }

        if (thisOp->theNumber) {
            if (thisOp->theNumber->ObjectClass() != NUMBER) {
                return false;
            }
        } else {
            if (thisOp->theData>-1) {
                _Variable* thisVar = LocateVar (thisOp->theData);
                if (thisVar->ObjectClass()!=NUMBER) {
                    _PMathObj cv = thisVar->GetValue();
                    if (!CheckSimpleTerm (cv)) {
                        return false;
                    }
                }
                if (variableIndex.Find (thisOp->theData)==-1) {
                    variableIndex<<thisOp->theData;
                }
            } else {
                if (simpleOperationCodes.Find(thisOp->opCode)==-1) {
                    return false;
                } else if ((thisOp->opCode == HY_OP_CODE_MACCESS || thisOp->opCode == HY_OP_CODE_MCOORD || thisOp->opCode == HY_OP_CODE_MUL) && thisOp->numberOfTerms != 2) {
                    return false;
                }

                locDepth-=thisOp->numberOfTerms;
            }
        }
        if (locDepth>stackDepth) {
            stackDepth = locDepth;
        } else if (locDepth==0L) {
            _String errStr = _String("Invalid formula (no return value) passed to _Formula::AmISimple") & _String ((_String*)toStr());
            WarnError (errStr);
            return false;
        }
    }
    return true;
}

//__________________________________________________________________________________
bool _Formula::IsArrayAccess (void){
    if (theFormula.lLength) {
        return ((_Operation*)(theFormula(theFormula.lLength-1)))->GetCode().Equal((_String*)BuiltInFunctions(HY_OP_CODE_MACCESS));
    }
    return false;
}

//__________________________________________________________________________________
bool _Formula::ConvertToSimple (_SimpleList& variableIndex)
{
    bool has_volatiles = false;
    if (theFormula.lLength)
        for (unsigned long i=0; i<theFormula.lLength; i++) {
            _Operation* thisOp = ((_Operation*)(((BaseRef**)theFormula.lData)[i]));
            if (thisOp->theNumber) {
                continue;
            } else if (thisOp->theData >= 0) {
                thisOp->theData = variableIndex.Find (thisOp->theData);
            } else if (thisOp->opCode == HY_OP_CODE_SUB && thisOp->numberOfTerms == 1) {
                thisOp->opCode = (long)MinusNumber;
            } else {
                if (thisOp->opCode == HY_OP_CODE_MACCESS) {
                    thisOp->numberOfTerms = -2;
                } else {
                  if (thisOp->opCode == HY_OP_CODE_MCOORD) {
                    thisOp->numberOfTerms = -3;
                  }
                }
                if (thisOp->opCode == HY_OP_CODE_RANDOM || thisOp->opCode == HY_OP_CODE_TIME)
                    has_volatiles = true;
                thisOp->opCode = simpleOperationFunctions(simpleOperationCodes.Find(thisOp->opCode));
            }

        }
    return has_volatiles;
}

//__________________________________________________________________________________
void _Formula::ConvertFromSimple (_SimpleList& variableIndex)
{
    if (!theFormula.lLength) {
        return;
    }

    for (int i=0; i<theFormula.lLength; i++) {
        _Operation* thisOp = ((_Operation*)(((BaseRef**)theFormula.lData)[i]));
        if (thisOp->theNumber) {
            continue;
        } else {
            if (thisOp->theData>-1) {
                thisOp->theData = variableIndex[thisOp->theData];
            } else if (thisOp->opCode == (long)MinusNumber) {
                thisOp->opCode = HY_OP_CODE_SUB;
            } else {
                if (thisOp->opCode == (long)FastMxAccess) {
                    thisOp->numberOfTerms = 2;
                } else {
                  if (thisOp->opCode == (long)FastMxWrite) {
                    thisOp->numberOfTerms = 2;
                  }
                }
                thisOp->opCode = simpleOperationCodes(simpleOperationFunctions.Find(thisOp->opCode));
            }
        }
    }
}

//__________________________________________________________________________________
_Parameter _Formula::ComputeSimple (_SimpleFormulaDatum* stack, _SimpleFormulaDatum* varValues)
{
    if (!theFormula.lLength) {
        return 0.0;
    }

    long stackTop = 0;

    for (int i=0; i<theFormula.lLength; i++) {
        _Operation* thisOp = ((_Operation*)(((BaseRef*)theFormula.lData)[i]));
        if (thisOp->theNumber) {
            stack[stackTop++].value = thisOp->theNumber->Value();
            continue;
        } else {
            if (thisOp->theData>-1) {
                stack[stackTop++] = varValues[thisOp->theData];
            } else {
                stackTop--;
                if (thisOp->numberOfTerms==2) {
                    _Parameter  (*theFunc) (_Parameter, _Parameter);
                    theFunc = (_Parameter(*)(_Parameter,_Parameter))thisOp->opCode;
                    if (stackTop<1L) {
                        WarnError ("Internal error in _Formula::ComputeSimple - stack underflow.)");
                        return 0.0;
                    }
                    stack[stackTop-1].value = (*theFunc)(stack[stackTop-1].value,stack[stackTop].value);
                } else {
                  switch (thisOp->numberOfTerms) {
                    case -2 : {
                        _Parameter  (*theFunc) (Ptr,_Parameter);
                        theFunc = (_Parameter(*)(Ptr,_Parameter))thisOp->opCode;
                        if (stackTop<1L) {
                            WarnError ("Internal error in _Formula::ComputeSimple - stack underflow.)");
                            return 0.0;
                        }
                        stack[stackTop-1].value = (*theFunc)(stack[stackTop-1].reference,stack[stackTop].value);
                        break;
                      }
                    case -3 : {
                      void  (*theFunc) (Ptr,_Parameter,_Parameter);
                      theFunc = (void(*)(Ptr,_Parameter,_Parameter))thisOp->opCode;
                      if (stackTop != 2 || i != theFormula.lLength - 1) {
                        WarnError ("Internal error in _Formula::ComputeSimple - stack underflow or MCoord command is not the last one.)");
                        
                        return 0.0;
                      }
                      //stackTop = 0;
                      // value, reference, index
                      (*theFunc)(stack[1].reference,stack[2].value, stack[0].value);
                      break;
                    }
                    default: {
                        _Parameter  (*theFunc) (_Parameter);
                        theFunc = (_Parameter(*)(_Parameter))thisOp->opCode;
                        stack[stackTop].value = (*theFunc)(stack[stackTop].value);
                        ++stackTop;
                    }
                  }
                
                }
            }
        }
    }

    return stack[0].value;
}

//__________________________________________________________________________________
bool _Formula::EqualFormula (_Formula* f)
{
    if (theFormula.lLength == f->theFormula.lLength) {
        for (int i=0; i<theFormula.lLength; i++) {
            if (!((_Operation*)(((BaseRef**)theFormula.lData)[i]))->EqualOp (((_Operation*)(((BaseRef**)f->theFormula.lData)[i])))) {
                return false;
            }
        }
        return true;
    }
    return false;
}

//__________________________________________________________________________________
_PMathObj _Formula::ConstructPolynomial (void) // compute the value of the formula
{
    theStack.Reset();
    bool wellDone = true;
    _String errMsg; 

    for (long i=0; i<theFormula.lLength; i++)
        if (!GetIthTerm(i)->ExecutePolynomial(theStack, nil, &errMsg)) {
            wellDone = false;
            break;
        }

    if (theStack.theStack.lLength !=1 || !wellDone) {
        return    nil;
    }

    return theStack.Pop(false);
}
//__________________________________________________________________________________
bool _Formula::HasChanged (bool ingoreCats) {
    _Operation *thisOp;
    for (unsigned long  i = 0UL; i<theFormula.lLength; i++) {
        long dataID;
        thisOp = (_Operation*)((BaseRef**)theFormula.lData)[i];
        if (thisOp->IsAVariable()) {
            dataID = thisOp->GetAVariable();
            if (dataID>=0) {
                if (((_Variable*)(((BaseRef*)(variablePtrs.lData))[dataID]))->HasChanged(ingoreCats)) {
                    return true;
                }
            } else if (thisOp->theNumber->HasChanged()) {
                return true;
            }
        } else if (thisOp->opCode == HY_OP_CODE_BRANCHLENGTH||thisOp->opCode == HY_OP_CODE_RANDOM||thisOp->opCode == HY_OP_CODE_TIME)
            // time, random or branch length
        {
            return true;
        } else if (thisOp->numberOfTerms<0L) {
            dataID = -thisOp->numberOfTerms-2;
            if (IsBFFunctionIndexValid (dataID)) {
                if (GetBFFunctionType (dataID) == BL_FUNCTION_SKIP_UPDATE) {
                    continue;
                }
            }
            return true;
        }
    }
    return false;
}


//__________________________________________________________________________________

void _Formula::ScanFormulaForHBLFunctions (_AVLListX& collection , bool recursive) {
  
  
  auto handle_function_id = [&collection, recursive] (const long hbl_id) -> void {
    if (IsBFFunctionIndexValid(hbl_id)) {
      _String function_name = GetBFFunctionNameByIndex(hbl_id);
      
      if (collection.Find(&function_name) < 0) {
        collection.Insert (new _String (function_name), HY_BL_HBL_FUNCTION, false, false);
        if (recursive) {
          GetBFFunctionBody(hbl_id).BuildListOfDependancies(collection, true);
        }
      }
    }
  };
  
  ConvertToTree();
  
  if (theTree) {
    
    InternalSimplify(theTree);
    node_iterator<long> ni (theTree, _HY_TREE_TRAVERSAL_PREORDER);
    
    while (node<long>* iterator = ni.Next()) {
      _Operation *this_op = GetIthTerm(iterator->get_data());
      
      long hbl_id = -1L;
      
      if (this_op -> IsHBLFunctionCall()) {
        hbl_id = this_op -> GetHBLFunctionID();
      } else {
        if (this_op->opCode == HY_OP_CODE_CALL) {
          node <long>* function_to_call = iterator->go_down (1);
          _Operation * function_to_call_value = GetIthTerm (function_to_call->get_data());
          
          if (function_to_call->get_num_nodes() == 0 && function_to_call_value -> IsConstantOfType(STRING)) {
              hbl_id = FindBFFunctionName (*((_FString*)function_to_call_value->theNumber->Compute())->theString);
          } else {
              ReportWarning ("Cannot export Call function arguments which are run-time dependent");
          }
        } else {
          if (this_op->opCode == HY_OP_CODE_MACCESS) { // handle AVL iterators
            if (this_op->GetNoTerms() == 3) { // [][]
              if (iterator->go_down(2)->get_num_nodes() == 0 && iterator->go_down(3)->get_num_nodes() == 0) {
                _Operation* bracket_1 = GetIthTerm (iterator->go_down(2)->get_data());
                _Operation* bracket_2 = GetIthTerm (iterator->go_down(3)->get_data());
                if (bracket_1->IsConstantOfType(STRING) && bracket_2->IsConstantOfType(STRING)) {
                  handle_function_id (FindBFFunctionName (*((_FString*)bracket_1->theNumber->Compute())->theString));
                  handle_function_id (FindBFFunctionName (*((_FString*)bracket_2->theNumber->Compute())->theString));
                  continue;
                }
              }
              ReportWarning ("Potentially missed dependence on a function in [][]; arguments are run-time dependent");
            }
          }
        }
      }
      
      handle_function_id (hbl_id);
      
      
      
    }
  }
}

//__________________________________________________________________________________
bool _Formula::HasChangedSimple (_SimpleList& variableIndex)
{
    _Operation *thisOp;
    for (unsigned long i = 0; i<theFormula.lLength; i++) {
        thisOp = (_Operation*)((BaseRef**)theFormula.lData)[i];
        if (thisOp->theNumber) {
            continue;
        } else if (thisOp->theData >= 0) {
            if (((_Variable*)(((BaseRef*)(variablePtrs.lData))[variableIndex.lData[thisOp->theData]]))->HasChanged(false)) {
                return true;
            }
        } else {
            if (thisOp->opCode == (long) RandomNumber) {
                return true;
            }
        }
    }
    return false;
}

//__________________________________________________________________________________
void _Formula::ScanFForVariables (_AVLList&l, bool includeGlobals, bool includeAll, bool includeCategs, bool skipMatrixAssignments, _AVLListX* tagger, long weight)
{
    for (unsigned long i = 0; i<theFormula.lLength; i++) {
        _Operation* theObj = ((_Operation**)theFormula.lData)[i];
        if (theObj->IsAVariable()) {
            if (!includeGlobals)
                // This change was part of a commit that introduced an optimizer bug (suspected location:
                // src/core/batchlan2.cpp:2220). This change is suspicious as well (removed and undocumented condition).
                //if ((((_Variable*)LocateVar(theObj->GetAVariable()))->IsGlobal())||
                 //       (((_Variable*)LocateVar(theObj->GetAVariable()))->ObjectClass()!=NUMBER)) {
                if (((_Variable*)LocateVar(theObj->GetAVariable()))->IsGlobal()) {
                    continue;
                }

            long f = theObj->GetAVariable();

            if (f>=0) {
                _Variable * v = LocateVar(f);

                if (v->IsCategory()&&includeCategs) {
                    v->ScanForVariables (l,includeGlobals,tagger, weight);
                }

                if(includeAll || v->ObjectClass()==NUMBER) {
                    l.Insert ((BaseRef)f);
                    if (tagger) {
                        tagger -> UpdateValue((BaseRef)f, weight, 0);
                    }
                }

                if (skipMatrixAssignments) {
                    if (v->ObjectClass()!=MATRIX || !theObj->AssignmentVariable()) {
                        v->ScanForVariables(l,includeGlobals,tagger, weight);
                    }
                } else if (!v->IsIndependent()) {
                    v->ScanForVariables(l,includeGlobals,tagger);
                }
            } else if (theObj->theNumber)
                if (theObj->theNumber->ObjectClass()==MATRIX) {
                    ((_Matrix*)theObj->theNumber)->ScanForVariables(l,includeGlobals,tagger, weight);
                }
        }
    }
}

//__________________________________________________________________________________
void _Formula::ScanFForType (_SimpleList &l, int type)
{
    for (long i = 0; i<theFormula.lLength; i++) {
        _Operation* theObj = ((_Operation**)theFormula.lData)[i];
        if (theObj->IsAVariable()) {
            long f = theObj->GetAVariable();

            if (f>=0) {
                _Variable * v = LocateVar(f);

                if(v->ObjectClass()==type) {
                    l << f;
                }

            }
        }
    }
}

//__________________________________________________________________________________
bool _Formula::CheckFForDependence (long varID, bool checkAll)
{
    for (int i = 0; i<theFormula.lLength; i++) {
        _Operation* theObj = (_Operation*)theFormula(i);
        long f;
        if (theObj->IsAVariable()) {
            f = theObj->GetAVariable();
            if (f>=0) {
                if (f == varID) {
                    return true;
                }
                if (checkAll) {
                    _Variable * v = LocateVar(f);
                    if (!v->IsIndependent())
                        if (v->CheckFForDependence(varID)) {
                            return true;
                        }
                }
            }
        }
    }
    return false;
}

//__________________________________________________________________________________
void  _Formula::LocalizeFormula (_Formula& ref, _String& parentName, _SimpleList& iv, _SimpleList& iiv, _SimpleList& dv, _SimpleList& idv)
{
    for (int i = 0; i<ref.theFormula.lLength; i++) {
        if (((_Operation*)ref.theFormula(i))->IsAVariable()) {
            long     vIndex = ((_Operation*)ref.theFormula(i))->GetAVariable();
            _Variable* theV = LocateVar (vIndex);
            if (theV->IsGlobal()) {
                theFormula&& ref.theFormula(i);
                continue;
            }
            if (theV->IsContainer()) {
                continue;
            }
            _String  fullName = parentName&"."&*(theV->GetName());
            long lvIndex = LocateVarByName(fullName);
            if (lvIndex==-1)
                // local variable doesn't yet exist - create
            {
                _Variable dummy (fullName);
                lvIndex = LocateVarByName(fullName);
                if (theV->IsIndependent()) {
                    iv<<lvIndex;
                    iiv<<vIndex;
                } else {
                    dv<<lvIndex;
                    idv<<vIndex;
                }
            }
            _Operation newVar (true, fullName);
            theFormula && &newVar;
        } else {
            theFormula&& ref.theFormula(i);
        }
    }
}

//__________________________________________________________________________________
bool _Formula::DependsOnVariable (long idx)
{
    for (long f = 0; f<theFormula.lLength; f++) {
        _Operation * anOp = ((_Operation**)theFormula.lData)[f];
        if (anOp->IsAVariable() && anOp->GetAVariable() == idx) {
            return true;
        }
    }

    return false;
}

//__________________________________________________________________________________
_Operation* _Formula::GetIthTerm (long idx) const {
    if (idx >= 0 && idx < theFormula.lLength) {
        return (_Operation*)theFormula.GetItem(idx);
    }

    return nil;
}

//__________________________________________________________________________________
bool _Formula::IsAConstant (void)
{
    for (unsigned long i = 0; i<theFormula.lLength; i++)
        if (((_Operation*)((BaseRef**)theFormula.lData)[i])->IsAVariable()) {
            return false;
        }

    return true;
}

//__________________________________________________________________________________
bool _Formula::IsConstant (void)
{
    for (unsigned long i = 0; i<theFormula.lLength; i++)
        if (((_Operation*)((BaseRef**)theFormula.lData)[i])->IsConstant() == false) {
            return false;
        }

    return true;
}

//__________________________________________________________________________________
void _Formula::PushTerm (BaseRef object) {
  _List * test_list;
  if ((test_list = dynamic_cast<_List*>(object))) {
    theFormula << *test_list;
  } else {
    theFormula << object;
  }
}

//__________________________________________________________________________________
void _Formula::SimplifyConstants (void){
  ConvertToTree    ();
  InternalSimplify (theTree);
  ConvertFromTree ();
}

//__________________________________________________________________________________
_PMathObj _Formula::GetTheMatrix (void)
{
    if (theFormula.lLength==1) {
        _Operation* firstOp = (_Operation*)theFormula(0);
        _PMathObj   ret = firstOp->GetANumber();
        if (ret&&(ret->ObjectClass()==MATRIX)) {
            return ret;
        } else {
            if (firstOp->theData!=-1) {
                _Variable* firstVar = LocateVar (firstOp->GetAVariable());
                ret = firstVar->GetValue();
                if (ret&&(ret->ObjectClass()==MATRIX)) {
                    return ret;
                }
            }
        }
    }
    return nil;
}

//__________________________________________________________________________________
long _Formula::ObjectClass (void)
{
    if (theStack.theStack.lLength) {
        return ((_PMathObj)theStack.theStack.lData[0])->ObjectClass();
    }

    _PMathObj res =   Compute();

    if (res) {
        return res->ObjectClass();
    }

    return HY_UNDEFINED;
}

//__________________________________________________________________________________
_Formula::_Formula (_String const &s, _VariableContainer const* theParent, _String* reportErrors) {
    theTree     = nil;
    resultCache = nil;
    recursion_calls = nil;
    call_count = 0UL;

    _FormulaParsingContext fpc (reportErrors, theParent);
  
    _String formula_copy (s);
  
    if (Parse (this, formula_copy, fpc, nil) != HY_FORMULA_EXPRESSION) {
        Clear();
    }
}
//__________________________________________________________________________________
void    _Formula::ConvertToTree (bool err_msg) {
    if (!theTree && theFormula.lLength) { // work to do
        _SimpleList nodeStack;
        _Operation* currentOp;
      
        for (unsigned long i=0UL; i<theFormula.lLength; i++) {
          
            currentOp = GetIthTerm(i);
          
            if (currentOp->theNumber || currentOp->theData >= 0L || currentOp->theData < -2L) { // a data bit
                node<long>* leafNode = new node<long>;
                leafNode->init(i);
                nodeStack<<(long)leafNode;
            } else { // an operation
                long nTerms = currentOp->GetNoTerms();
                if (nTerms<0L) {
                    nTerms = GetBFFunctionArgumentCount(currentOp->opCode);
                }

                if (nTerms>nodeStack.lLength) {
                    if (err_msg) {
                        WarnError (_String ("Insufficient number of arguments for a call to ") & _String ((_String*)currentOp->toStr()));
                    }
                    theTree = nil;
                    return;
                }

                node<long>* operationNode = new node<long>;
                operationNode->init(i);
                for (long j=0; j<nTerms; j++) {
                    operationNode->prepend_node(*((node<long>*)nodeStack.Pop()));
                }
                nodeStack<<(long)operationNode;
            }
        }
        if (nodeStack.lLength!=1) {
            if (err_msg) {
                WarnError ((_String)"The expression '" & _String ((_String*)toStr()) & "' has " & (long)nodeStack.lLength & " terms left on the stack after evaluation");
            }
            theTree = nil;
        } else {
            theTree = (node<long>*)nodeStack(0);
        }
      
      
    }
}
//__________________________________________________________________________________
void    _Formula::ConvertFromTree (void)
{
    if (theTree) { // work to do
        _SimpleList termOrder;
        node_iterator<long> ni (theTree, _HY_TREE_TRAVERSAL_POSTORDER);
      
      while (node<long>* iterator = ni.Next()) {
            termOrder<<iterator->get_data();
        }
      
        if (termOrder.lLength!=theFormula.lLength) { // something has changed
            _List newFormula;
            for (long i=0; i<termOrder.lLength; i++) {
                newFormula<<theFormula(termOrder(i));
            }
            theFormula.Clear();
            theFormula.Duplicate(&newFormula);
            //ConvertToTree();
        }
        theTree->delete_tree();
        delete (theTree);
        theTree = nil;
    }
}

//__________________________________________________________________________________
node<long>* _Formula::InternalDifferentiate (node<long>* currentSubExpression, long varID, _SimpleList& varRefs, _SimpleList& dydx, _Formula& tgt)
{
    _Operation * op = (_Operation*)theFormula (currentSubExpression->in_object);

    if (op->theData!=-1) {
        long k     = varRefs.BinaryFind (op->GetAVariable());
        if (k<0) {
            return nil;
        }

        _Formula* dYdX = (_Formula*)dydx(k);
        return dYdX->DuplicateFormula (dYdX->theTree, tgt);
    }

    if (op->theNumber) {
        _Formula src (new _Constant (0.0));
        src.ConvertToTree ();
        return   src.DuplicateFormula (src.theTree, tgt);
    }

    node<long>* newNode = (node<long>*)checkPointer (new node<long>);


    switch (op->opCode) {
    case HY_OP_CODE_MUL: {
      
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt),
                  * b2 = InternalDifferentiate (currentSubExpression->go_down(2), varID, varRefs, dydx, tgt);

        if (!b1 || !b2) {
            newNode->delete_tree(true);
            if (b1) {
                b1->delete_tree (true);
            }
            if (b2) {
                b2->delete_tree (true);
            }
            return nil;
        }

        _String           opC  ('*'),
                          opC2 ('+');

        _Operation*       newOp  = new _Operation (opC2,2),
        *         newOp2 = new _Operation (opC ,2),
        *         newOp3 = new _Operation (opC ,2);

      
      
      
        node<long>*       newNode2 = (node<long>*)checkPointer(new node<long>);
        node<long>*       newNode3 = (node<long>*)checkPointer(new node<long>);

        newNode2->add_node (*b1);
        newNode2->add_node (*DuplicateFormula (currentSubExpression->go_down(2),tgt));

        newNode3->add_node (*b2);
        newNode3->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));

        newNode->add_node  (*newNode2);
        newNode->add_node  (*newNode3);

        newNode3->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance (newOp3);
        newNode2->in_object = tgt.theFormula.lLength;
        tgt.theFormula. AppendNewInstance (newOp2);
        newNode->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance (newOp);

        return          newNode;
    }
    break;

    case HY_OP_CODE_ADD: // +
    case HY_OP_CODE_SUB: { // -
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt),
                    * b2 = nil;

        if (!b1) {
            newNode->delete_tree (true);
            return nil;
        }

        long      isUnary = (currentSubExpression->get_num_nodes()==1);

        if (!isUnary) {
            b2 = InternalDifferentiate (currentSubExpression->go_down(2), varID, varRefs, dydx, tgt);
            if (!b2) {
                b1->delete_tree      (true);
                newNode->delete_tree (true);
                return nil;
            }
        }


        _Operation*   newOp = new _Operation ();
        checkPointer  (newOp);
        newOp->Duplicate (op);
        newNode->add_node (*b1);
        if (!isUnary) {
            newNode->add_node (*b2);
        }
        newNode->in_object = tgt.theFormula.lLength;

        tgt.theFormula.AppendNewInstance(newOp);
        return          newNode;
    }
    break;

    case HY_OP_CODE_DIV: { // /
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt),
                    * b2 = InternalDifferentiate (currentSubExpression->go_down(2), varID, varRefs, dydx, tgt);

        if (!b1 || !b2) {
            newNode->delete_tree(true);
            if (b1) {
                b1->delete_tree (true);
            }
            if (b2) {
                b2->delete_tree (true);
            }
            return nil;
        }

        _String           opC  ('*'),
                          opC2 ('-'),
                          opC3 ('/'),
                          opC4 ('^');

        _Operation*       newOp  = new _Operation (opC3 ,2),
        *         newOp2 = new _Operation (opC4 ,2),
        *         newOp3 = new _Operation (opC2 ,2),
        *         newOp4 = new _Operation (opC  ,2),
        *         newOp5 = new _Operation (opC  ,2),
        *         newOp6 = new _Operation (new _Constant (2.0));

 
        node<long>*       newNode2 = new node<long>;
        node<long>*       newNode3 = new node<long>;
        node<long>*       newNode4 = new node<long>;
        node<long>*       newNode5 = new node<long>;
        node<long>*       newNode6 = new node<long>;

        newNode6->add_node (*b1);
        newNode6->add_node (*DuplicateFormula (currentSubExpression->go_down(2),tgt));

        newNode5->add_node (*b2);
        newNode5->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));

        newNode4->add_node  (*newNode6);
        newNode4->add_node  (*newNode5);

        newNode2->add_node (*DuplicateFormula (currentSubExpression->go_down(2),tgt));
        newNode2->add_node (*newNode3);

        newNode->add_node  (*newNode4);
        newNode->add_node  (*newNode2);

        newNode6->in_object = tgt.theFormula.lLength;
      
      
        tgt.theFormula.AppendNewInstance(newOp5);
        newNode5->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp4);
        newNode4->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp3);
        newNode3->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp6);
        newNode2->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp2);
        newNode->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp);


        return          newNode;
    }
    break;

    case HY_OP_CODE_ARCTAN: { // Arctan
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);

        if (!b1) {
            newNode->delete_tree (true);
            return nil;
        }

        _String           opC  ('/'),
                          opC2 ('+'),
                          opC3 ('^');

        _Operation*       newOp  = new _Operation (opC ,2),
        *         newOp2 = new _Operation (opC2 ,2),
        *         newOp3 = new _Operation (new _Constant (1.0)),
        *         newOp4 = new _Operation (opC3 ,2),
        *         newOp5 = new _Operation (new _Constant (2.0));

        checkPointer      (newOp);
        checkPointer      (newOp2);
        checkPointer      (newOp3);
        checkPointer      (newOp4);

        node<long>*       newNode2 = new node<long>;
        node<long>*       newNode3 = new node<long>;
        node<long>*       newNode4 = new node<long>;
        node<long>*       newNode5 = new node<long>;

        checkPointer      (newNode2);
        checkPointer      (newNode3);
        checkPointer      (newNode4);
        checkPointer      (newNode5);

        newNode4->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));
        newNode4->add_node (*newNode5);

        newNode2->add_node (*newNode3);
        newNode2->add_node (*newNode4);

        newNode->add_node (*b1);
        newNode->add_node (*newNode2);

        newNode5->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp5);
        newNode4->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp4);
        newNode3->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp3);
        newNode2->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp2);
        newNode->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp);

        return          newNode;
    }
    break;

    case HY_OP_CODE_COS: { // Cos
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);

        if (!b1) {
            newNode->delete_tree (true);
            return nil;
        }

        _String           opC  ('*'),
                          opC2 ('-'),
                          opC3 ("Sin");

        _Operation*       newOp  = new _Operation (opC  ,2),
        *         newOp2 = new _Operation (opC2 ,1),
        *         newOp3 = new _Operation (opC3 ,1);

        checkPointer      (newOp);
        checkPointer      (newOp2);
        checkPointer      (newOp3);

        node<long>*       newNode2 = new node<long>;
        node<long>*       newNode3 = new node<long>;

        checkPointer      (newNode2);
        checkPointer      (newNode3);

        newNode3->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));

        newNode2->add_node (*newNode3);

        newNode->add_node  (*newNode2);
        newNode->add_node  (*b1);

        newNode3->in_object = tgt.theFormula.lLength;
        tgt.theFormula << newOp3;
        newNode2->in_object = tgt.theFormula.lLength;
        tgt.theFormula << newOp2;
        newNode->in_object = tgt.theFormula.lLength;
        tgt.theFormula << newOp;

        DeleteObject    (newOp);
        DeleteObject    (newOp2);
        DeleteObject    (newOp3);

        return          newNode;
    }
    break;

    case HY_OP_CODE_ERF: { // Erf
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);

        if (!b1) {
            newNode->delete_tree (true);
            return nil;
        }

        _String           opC  ('*'),
                          opC2 ('/'),
                          opC3 ("Exp"),
                          opC4 ("-"),
                          opC5 ("^");

        _Operation*       newOp  = new _Operation (opC  ,2),
        *         newOp2 = new _Operation (opC2 ,2),
        *         newOp3 = new _Operation (new _Constant (twoOverSqrtPi)),
        *         newOp4 = new _Operation (opC3 ,1),
        *         newOp5 = new _Operation (opC4 ,1),
        *         newOp6 = new _Operation (opC5 ,2),
        *         newOp7 = new _Operation (new _Constant (2.0));

        checkPointer      (newOp);
        checkPointer      (newOp2);
        checkPointer      (newOp3);
        checkPointer      (newOp4);
        checkPointer      (newOp5);
        checkPointer      (newOp6);
        checkPointer      (newOp7);

        node<long>*       newNode2 = new node<long>;
        node<long>*       newNode3 = new node<long>;
        node<long>*       newNode4 = new node<long>;
        node<long>*       newNode5 = new node<long>;
        node<long>*       newNode6 = new node<long>;
        node<long>*       newNode7 = new node<long>;

        checkPointer      (newNode2);
        checkPointer      (newNode3);
        checkPointer      (newNode4);
        checkPointer      (newNode5);
        checkPointer      (newNode6);
        checkPointer      (newNode7);

        newNode6->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));
        newNode6->add_node (*newNode7);

        newNode5->add_node (*newNode6);
        newNode4->add_node (*newNode5);
        newNode2->add_node (*newNode4);
        newNode2->add_node (*newNode3);

        newNode->add_node  (*b1);
        newNode->add_node  (*newNode2);

        newNode7->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp7);
        newNode6->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp6);
        newNode5->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp5);
        newNode4->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp4);
        newNode3->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp3);
        newNode2->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp2);
        newNode->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp);

        return          newNode;
    }
    break;

    case HY_OP_CODE_EXP: // HY_OP_CODE_EXP
    case HY_OP_CODE_SIN: { // HY_OP_CODE_SIN
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);

        if (!b1) {
            newNode->delete_tree (true);
            return nil;
        }

        _String           opC  ('*'),
                          opC2;

        if (op->opCode==HY_OP_CODE_SIN) {
            opC2 = *(_String*)BuiltInFunctions(HY_OP_CODE_COS);
        } else {
            opC2 = *(_String*)BuiltInFunctions(HY_OP_CODE_EXP);
        }

        _Operation*       newOp  = new _Operation (opC  ,2),
        *         newOp2 = new _Operation (opC2 ,1);

        checkPointer      (newOp);
        checkPointer      (newOp2);

        node<long>*       newNode2 = new node<long>;

        checkPointer      (newNode2);

        newNode2->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));

        newNode->add_node  (*newNode2);
        newNode->add_node  (*b1);

        newNode2->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp2);
        newNode->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp);

        return           newNode;
    }
    break;

    case HY_OP_CODE_LOG: { // Log
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);

        if (!b1) {
            newNode->delete_tree (true);
            return nil;
        }
        _String           opC  ('/');

        _Operation*       newOp  = new _Operation (opC  ,2);

        checkPointer      (newOp);

        newNode->add_node  (*b1);
        newNode->add_node  (*DuplicateFormula (currentSubExpression->go_down(1),tgt));

        newNode->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp);

        return           newNode;
    }
    break;

    case HY_OP_CODE_SQRT: { // Sqrt
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);

        if (!b1) {
            newNode->delete_tree (true);
            return nil;
        }

        _String           opC  ('/'),
                          opC2 ('*'),
                          opC3 (*(_String*)BuiltInFunctions(HY_OP_CODE_SQRT));

        _Operation*       newOp  = new _Operation (opC  ,2),
        *         newOp2 = new _Operation (opC2 ,2),
        *         newOp3 = new _Operation (opC3 ,1),
        *         newOp4 = new _Operation (new _Constant (2.0));

        checkPointer      (newOp);
        checkPointer      (newOp2);
        checkPointer      (newOp3);
        checkPointer      (newOp4);

        node<long>*       newNode2 = new node<long>;
        node<long>*       newNode3 = new node<long>;
        node<long>*       newNode4 = new node<long>;

        checkPointer      (newNode2);
        checkPointer      (newNode3);
        checkPointer      (newNode4);

        newNode3->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));

        newNode2->add_node (*newNode4);
        newNode2->add_node (*newNode3);

        newNode->add_node  (*b1);
        newNode->add_node  (*newNode2);

        newNode4->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp4);
        newNode3->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp3);
        newNode2->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp2);
        newNode->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp);

        return          newNode;
    }
    break;

    case HY_OP_CODE_TAN: { // Tan
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);

        if (!b1) {
            newNode->delete_tree (true);
            return nil;
        }

        _String           opC  ('/'),
                          opC2 ('^'),
                          opC3 (*(_String*)BuiltInFunctions(HY_OP_CODE_COS));

        _Operation*       newOp  = new _Operation (opC ,2),
        *         newOp2 = new _Operation (opC2,2),
        *         newOp3 = new _Operation (new _Constant (2.0)),
        *         newOp4 = new _Operation (opC3,1);

        checkPointer      (newOp);
        checkPointer      (newOp2);
        checkPointer      (newOp3);

        node<long>*       newNode2 = new node<long>;
        node<long>*       newNode3 = new node<long>;
        node<long>*       newNode4 = new node<long>;

        checkPointer      (newNode2);
        checkPointer      (newNode3);
        checkPointer      (newNode4);

        newNode4->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));

        newNode2->add_node (*newNode4);
        newNode2->add_node (*newNode3);

        newNode->add_node (*b1);
        newNode->add_node (*newNode2);

        newNode4->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp4);
        newNode3->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp3);
        newNode2->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp2);
        newNode->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp);

        return          newNode;
    }
    break;

    case HY_OP_CODE_POWER: // ^
        // f[x]^g[x] (g'[x] Log[f[x]] + f'[x]g[x]/f[x])
    {
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);
        if (!b1) {
            newNode->delete_tree(true);
            return nil;
        }
        node<long>* b2 = InternalDifferentiate (currentSubExpression->go_down(2), varID, varRefs, dydx, tgt);
        if (!b2) {
            newNode->delete_tree(true);
            return nil;
        }
        
        long opCodes[7] = {HY_OP_CODE_MUL, HY_OP_CODE_POWER, HY_OP_CODE_ADD, HY_OP_CODE_DIV, HY_OP_CODE_MUL, HY_OP_CODE_MUL, HY_OP_CODE_LOG};
        long opArgs [7] = {2,2,2,2,2,2,1};

        node<long> * newNodes [7];
        newNodes[0] = newNode;
        for (long k = 1; k < 7; k++) {
            newNodes[k] = new node <long>;
        }

        newNodes[6]->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));

        newNodes[5]->add_node (*b2);
        newNodes[5]->add_node (*newNodes[6]);

        newNodes[4]->add_node (*b1);
        newNodes[4]->add_node (*DuplicateFormula (currentSubExpression->go_down(2),tgt));

        newNodes[3]->add_node (*newNodes[4]);
        newNodes[3]->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));

        newNodes[2]->add_node (*newNodes[5]);
        newNodes[2]->add_node (*newNodes[3]);

        newNodes[1]->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));
        newNodes[1]->add_node (*DuplicateFormula (currentSubExpression->go_down(2),tgt));

        newNode->add_node  (*newNodes[1]);
        newNode->add_node  (*newNodes[2]);
        
        for (long k = 6; k >=0 ; k--) {
            newNodes[k]->in_object = tgt.theFormula.lLength;
            tgt.theFormula .AppendNewInstance (new _Operation (opCodes[k], opArgs[k]));
        } 

        return          newNode;
    }
    }

    delete (newNode);
    return nil;
}



//__________________________________________________________________________________
_FormulaParsingContext::_FormulaParsingContext (_String* err, _VariableContainer const* scope) {
    assignment_ref_id   = -1;
    assignment_ref_type = HY_STRING_DIRECT_REFERENCE;
    is_volatile = false;
    in_assignment = false;
    build_complex_objects = true;
    err_msg = err;
    formula_scope = scope;
}

//__________________________________________________________________________________
_String const _FormulaParsingContext::contextualizeRef (_String& ref) {
    if (formula_scope) {
        return *formula_scope->GetName () & '.' & ref;
    }
    return ref;
}

//__________________________________________________________________________________
void _FormulaParsingContext::setScope(const _String *scope) {
  if (scope && scope->sLength) {
    _VariableContainer vc (*scope);
    formula_scope = (_VariableContainer*)FetchVar(vc.GetAVariable());
  } else {
    formula_scope = nil;
  }
}
