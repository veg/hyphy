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

#include    "scfg.h"
#include "function_templates.h"
#include "global_things.h"
#include "vector.h"

using namespace hy_global;

const _String     kSCFG_TERM_T       ("T"),
                  kSCFG_TERM_P       ("P"),
                  kSCFG_TERM_L       ("L"),
                  kSCFG_TERM_NT_1    ("1"),
                  kSCFG_TERM_NT_2    ("2");


#include    <math.h>


// _String      covariancePrecision ("COVARIANCE_PRECISION");

#define     _HYSCFG_CHARACTER_MASK_  0x000000ff // used to extract the character part from the 'long' data field of the parse tree
#define     _HYSCFG_TERMINAL_MASK_   0xffffff00 // used to extract the terminal index part from the 'long' data field of the parse tree
#define     _HYSCFG_NT_LHS_          0x00000001 // a flag used to indicate that a given non-terminal appears on the LHS of a production rule
#define     _HYSCFG_NT_DTERM_        0x00000002 // a flag used to indicate that a given non-terminal is on the LHS of a <NT>->Terminal production rule
#define     _HYSCFG_NT_START_        0x00000004 // a flag used to indicate that a given non-terminal can be reached from a start symbol after some productions
#define     _HYSCFG_NT_TERM_         0x00000008 // a flag used to indicate that a given non-terminal can be resolved to at least one terminal after some productions
#define     _HYSCFG_UPPER_BOUND_     2147483648.// maximum indexable; constraining the maximum length of input strings 

#define     SCFG_OPTIMIZATION_THRESHOLD     0.0001  // minimum difference between log-likelihood to terminate
#define     MAX_ITERATIONS_OPTIMIZE         1000        // maximum number of times to iterate Optimize() main loop

bitMasks bitMaskArray;







/*          ================                                                        */
/* ========| SCFG FUNCTIONS |====================================================== */
/*          ================                                                        */

#ifdef                  _USE_HYPHY_HOOKS_

Scfg::Scfg  (_AssociativeList* T_Rules,  _AssociativeList* NT_Rules, long ss) {
  
     _SimpleList     parsedFormulas;
    // stash pointers to processed production expressions formulas here
    
    startSymbol     = ss;
    insideCalls     = 0;
    outsideCalls    = 0;
    
    // initialize the pointers
    
    
    long          termRules     = T_Rules->countitems(),
    nonTermRules  = NT_Rules->countitems();
    
    try {
    
      if (termRules == 0L) {
        throw _String("A SCFG can not be constructed from an empty set of <Nonterminal>-><Terminal> production rules");
      }
      
      _List         ruleProbabilities,            // an auxiliary list used to store strings describing production probabilities
      alreadySeenX;                 // a list of production rules that have already been added
                                    // it is used to keep track of duplicate production rules
      
      _SimpleList   foundNT,                      // an array to keep track of all 'declared' non-terminals
      ntFlags;                      // an array of flags of non-terminal properties
      
      _AVLList      tempTerminals (&terminals),       // an auxiliary wrapper around the set of terminal symbols used to check for duplication
      alreadySeen   (&alreadySeenX);    // and a wrapper for the alreadySeen list
      
      _AVLListX     tempNT        (&foundNT);         // and a wrapper for the set of non-terminals
      
      
      
      for (long tc = 0L; tc < termRules; tc++) {
          // check Terminal rules for consistency
          // traverse the T_Rules AVL, which is assumed to be indexed by 0..termRules-1
        _AssociativeList * aRule = (_AssociativeList*)T_Rules->GetByKey (tc, ASSOCIATIVE_LIST);
        
        if (aRule == nil) {
          throw _String("Each production rule must be specified as an associative array, but rule ") & tc & " was not.";
        }
        _FString    * literal       = (_FString*)aRule->GetByKey (kSCFG_TERM_T,STRING),
        * expression   = (_FString*)aRule->GetByKey (kSCFG_TERM_P, STRING);
        
        _Constant   * lhs           = (_Constant*)aRule->GetByKey (kSCFG_TERM_L, NUMBER);
        
        if (! (literal && lhs && !literal->empty())) {
          throw _String("Each terminal rule must have a non-empty target terminal and a left-hand side non-terminal. Rule ") & tc & " did not comply.";
        }
        
        
        long index = tempTerminals.Insert (literal->get_str().makeDynamic(), -1, false, true);
          // insert to the list of terminals if not seen before
        if  (index>=0) {// new terminal; added to list
                       // now we process the terminal into the parse tree
          if (terminal_symbols.FindKey(literal->get_str(), nil, true) != kNotFound) {
             throw  _String("Terminal symbols must form a prefix code, but ") & *literal->get_str().Enquote() & " contains another terminal symbol as a prefix.";
          }
          
          terminal_symbols.Insert(literal->get_str(), index);
        } else {
          index = -index-1;    // if this terminal has already been added, use the index
        }
        
        long          nt_index  = (long)(lhs->Compute()->Value()),
        avl_index = tempNT.Insert ((BaseRef)nt_index); // store the integer index of the LHS if needed
        
        if (avl_index<0) { // nt_index already exists; correct to positive range
          avl_index = -avl_index - 1;
        }
        
          // update status flags for this non-terminal
        tempNT.SetXtra (avl_index, tempNT.GetXtra (avl_index)|_HYSCFG_NT_LHS_|_HYSCFG_NT_DTERM_|_HYSCFG_NT_TERM_);
        
        
          // first ensure the rule is not a duplicate
        _String         ruleString = _String (nt_index) & ",[" & index & ']';
          // create a new record for the rule of the form "nt index,[t index]"
        long            seenMe     = alreadySeen.Insert (ruleString.makeDynamic(), -1, false, true);
        if (seenMe < 0L) { // already seen
          throw _String ("Duplicate production rule:" ) & GetRuleString (-seenMe-1);
        }
        
        rules < & ((*new _SimpleList ()) << nt_index << index);
        ProcessAFormula (expression, ruleProbabilities, parsedFormulas);

      } // done checking terminal rules
      
      for (long tc = 0L; tc < nonTermRules; tc++) {
          // check Non-terminal rules for consistency
          // traverse the NT_Rules AVL, which is assumed to be indexed by 0..nonTermRules-1
        _AssociativeList * aRule = (_AssociativeList*)NT_Rules->GetByKey (tc, ASSOCIATIVE_LIST);
        if (!aRule) {
          throw _String("Each rule must be specified as an associative array, but rule ") & tc & " was not.";
        }
        _FString    * expression    = (_FString*)aRule->GetByKey     (kSCFG_TERM_P, STRING);
        
        _Constant   * lhs           = (_Constant*)aRule->GetByKey    (kSCFG_TERM_L, NUMBER),
        * rhs1          = (_Constant*)aRule->GetByKey    (kSCFG_TERM_NT_1, NUMBER),
        * rhs2          = (_Constant*)aRule->GetByKey    (kSCFG_TERM_NT_2, NUMBER);
        
        if (! (lhs && rhs1 && rhs2)) {
          throw _String("Each non-terminal rule must have two right-hand side non-terminals and a left-hand side non-terminal. Rule ") & tc & " did not comply.";
        }
        
        ProcessAFormula (expression, ruleProbabilities, parsedFormulas);
        
        _SimpleList goodNTRule;
        long          nt_index = (long)(lhs->Compute()->Value());
        long avl_index = tempNT.Insert ((BaseRef)nt_index); // store the integer index of the LHS if needed
        if (avl_index<0) {
          avl_index = -avl_index - 1;
        }
        tempNT.SetXtra (avl_index, tempNT.GetXtra (avl_index)|_HYSCFG_NT_LHS_); // update status flags for this non-terminal
        
        _StringBuffer     ruleString;
        ruleString << _String(nt_index) << ',';
        goodNTRule << nt_index;
        
        nt_index    = (long)(rhs1->Compute()->Value());
        tempNT.Insert ((BaseRef)nt_index);
        goodNTRule << nt_index;
        ruleString << _String(nt_index) << ',';
        
        nt_index    = (long)(rhs2->Compute()->Value());
        tempNT.Insert ((BaseRef)nt_index);
        goodNTRule << nt_index;
        ruleString <<_String (nt_index);
        
        long            seenMe     = alreadySeen.Insert (new _String(ruleString), -1, false, true);
        if (seenMe < 0) { // already seen
          throw _String ("Duplicate production rule: ") & tc & " : " & ruleString.Enquote();
        }
        rules < new _SimpleList (goodNTRule); // append the new rule to the list of existing rules
      }
      
       // all rules were good; next steps:
       // (a). Validate the grammar
       // (1).   Each non-terminal must appear in at least one LHS
       // (2).   Each non-terminal must be resolvable to a terminal at some point
       // (3).   Each non-terminal must be reachable from the start symbol
        
        ntToTerminalMap.Populate (rules.lLength*terminals.lLength,-1,0);
        
        bool         continueLoops = true;
      
          // prepopulate the list of rules stratified by the LHS non-terminal by empty lists
        for (long ntC = 0L; ntC < foundNT.lLength; ntC ++) {
          byNT2 < new _SimpleList;
          byNT3 < new _SimpleList;
          byRightNT1 < new _SimpleList;
          byRightNT2 < new _SimpleList;
        }
        
        for (long ruleIdx = 0L; ruleIdx < rules.lLength; ruleIdx ++) {
          _SimpleList *aList = (_SimpleList*)rules(ruleIdx);      // retrieve two- or three-integer list
          
          if (aList->countitems() == 3UL) { // NT->NT NT
            *((_SimpleList*)byNT3 (aList->get(0))) << ruleIdx;
            /* addition by AFYP, 2006-06-20 */
            *((_SimpleList*)byRightNT1 (aList->get(1))) << ruleIdx;
            *((_SimpleList*)byRightNT2 (aList->get(2))) << ruleIdx;
            /* ---------- end add --------- */
          } else {
            *((_SimpleList*)byNT2 (aList->get(0))) << ruleIdx;
            ntToTerminalMap [indexNT_T (aList->list_data[0],aList->get(1))] = ruleIdx;
          }
          
        }
        
        while (continueLoops) { // set status flags for all NT symbols based on production rules
          continueLoops = false;
          for (long ruleIdx = 0L; ruleIdx < rules.lLength; ruleIdx ++) {
            _SimpleList *aList = (_SimpleList*)rules(ruleIdx);
            if (aList->countitems() == 3UL) { // NT->NT NT
              continueLoops = continueLoops || CheckANT (aList->get(0),aList->get(1), aList->get(2), tempNT, startSymbol);
            }
          }
        }
        
        
          // now iterate over the list of declared NT and verify that they all comply to the three conditions above
        for (long ntC = 0L; ntC < foundNT.lLength; ntC ++) {
          long ntFlag = tempNT.GetXtra (ntC);
          if ((ntFlag & _HYSCFG_NT_LHS_) == 0L) {
            throw _String ("Non-terminal symbol ") & foundNT.list_data[ntC] & " does not appear on the left-hand side of any production rules.";
          }
          if ((ntFlag & _HYSCFG_NT_START_) == 0L) {
            throw _String ("Non-terminal symbol ") & foundNT.list_data[ntC] & " can not be derived from the start symbol.";
          }
          if ((ntFlag & _HYSCFG_NT_TERM_) == 0L) {
            throw _String ("Non-terminal symbol ") & foundNT.list_data[ntC] & " can not be used to derive any terminal symbols.";
          }
        }
      
          _Matrix::CreateMatrix (&probabilities, parsedFormulas.lLength, 1, false, true, false);
          probabilities.Convert2Formulas ();
          parsedFormulas.Each ([&] (long fl_ref, unsigned long index) ->  void {
            probabilities.StoreFormula (index, 0, *(_Formula*)fl_ref);
          });
      
        auto get_rule_idx = [&] (long i, long j) -> long {
          return ((_SimpleList*)rules.GetItem (i))->get(j);
        };
      
          
          long countT  = terminals.countitems(),
               countNT = byNT2.countitems();
          
          firstArray.Populate (countNT*countT,0,0);
          lastArray.Populate (countNT*countT,0,0);
          precursorArray.Populate (countNT*countT,0,0);
          followArray.Populate (countNT*countT,0,0);
          
          for (long i = 0L; i < countNT; i++) {
            _SimpleList * myRules = (_SimpleList*)byNT2.GetItem (i);
            for (long i2 = 0L; i2 < myRules->lLength; i2++) {    // for all i->m productions
              long flatIndex = indexNT_T (i, get_rule_idx (myRules->get (i2), 1));
              firstArray  [flatIndex]  = 1L;
              lastArray   [flatIndex]  = 1L;
            }
          }
          
          continueLoops = true;
          while (continueLoops) {
            continueLoops = false;
            for (long i = 0L; i < countNT; i++) {
              _SimpleList * myRules = (_SimpleList*)byNT3.GetItem(i);
              for (long i2 = 0L; i2 < myRules->lLength; i2++) {
                
                long rhs1 = get_rule_idx (myRules->get (i2), 1),
                     rhs2 = get_rule_idx (myRules->get (i2), 2);
                
                for (long i3 = 0L; i3 < countT; i3++) {
                  long anIndex = indexNT_T (i,i3),
                  rhs1v = firstArray.get(indexNT_T (rhs1,i3)),
                  lhs1v = firstArray.get(anIndex),
                  rhs2v = lastArray.get (indexNT_T (rhs2,i3)),
                  lhs2v = lastArray.get (anIndex);
                  
                  if (lhs1v == 0 && rhs1v) {
                    continueLoops = true;
                    firstArray[anIndex] = 1;
                  }
                  if (lhs2v == 0 && rhs2v) {
                    continueLoops = true;
                    lastArray[anIndex] = 1;
                  }
                }
              }
            }
          }
      
          for (long i = 0L; i < countNT; i++) { // initialize Precursor and Follow
            _SimpleList * myRules = (_SimpleList*)byNT3.GetItem(i);
            for (long i2 = 0L; i2 < myRules->lLength; i2++) {
              long rhs1 = get_rule_idx (myRules->get (i2), 1),
                   rhs2 = get_rule_idx (myRules->get (i2), 2);
              
              for (long i3 = 0L; i3 < countT; i3++) {
                long idx1 = indexNT_T (rhs1,i3),
                     idx2 = indexNT_T (rhs2,i3);
                
                if (lastArray.get (idx1)) {
                  precursorArray[idx2] = 1;
                }
                if (firstArray.get(idx2)) {
                  followArray[idx1]    = 1;
                }
              }
            }
          }
      
          
          continueLoops = true;
          while (continueLoops) { // populate Precursor and Follow
            continueLoops = false;
            for (long i = 0L; i < countNT; i++) {
              _SimpleList * myRules = (_SimpleList*)byNT3.GetItem(i);
              for (long i2 = 0L; i2 < myRules->lLength; i2++) {
                long rhs1 = get_rule_idx (myRules->get (i2), 1),
                     rhs2 = get_rule_idx (myRules->get (i2), 2);
                
                for (long i3 = 0L; i3 < countT; i3++) {
                  long anIndex = indexNT_T (i,i3),
                  idx1    = indexNT_T (rhs1,i3),
                  idx2    = indexNT_T (rhs2,i3),
                  rhs1v   = precursorArray.get(idx1),
                  rhs2v   = followArray.get   (idx2),
                  lhs1v   = precursorArray.get(anIndex),
                  lhs2v   = followArray.get   (anIndex);
                  
                  if (lhs1v && rhs1v == 0) {
                    continueLoops = true;
                    precursorArray[idx1] = 1;
                  }
                  if (lhs2v && rhs2v == 0) {
                    continueLoops = true;
                    followArray[idx2] = 1;
                  }
                }
              }
            }
          }
      
          parsedFormulas.ClearFormulasInList();
    } catch (const _String & err) {
      HandleApplicationError (err);
    }
  
    ScanAllVariables   ();
    // RandomSampleVerify (100);
    /* temporarily removed this for node censoring (degenerate grammar) -- AFYP, August 30, 2006 */
}


/*--------------------------------------------------------------------------------------------------------------------------------*/

long        Scfg::indexNT_T (long nt, long t) const {
    return nt*terminals.lLength+t;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/
 void    Scfg::ProcessAFormula  (_FString* expression, _List & ruleProbabilities, _SimpleList& parsedFormulas) {
    _Formula *aFormula ;
    if (expression) { // probabilistic rule
        aFormula = new _Formula;

        
        _String  anExpression (expression->get_str());

        _Formula lhs;
        _FormulaParsingContext fpc;

        if (Parse   (aFormula, anExpression, fpc, &lhs) != HY_FORMULA_EXPRESSION) { // not a valid expression
            throw _String ("Invalid probability expression ") & expression->get_str().Enquote();
        } else {
            ruleProbabilities < new _String(expression->get_str());
        }
    } else { // determininstic rule (prob = 1.0)
        aFormula = new _Formula (new HY_CONSTANT_TRUE, false); // constant 1.0
        ruleProbabilities < new _String (kSCFG_TERM_NT_1);
    }

    parsedFormulas << (long)aFormula;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

bool    Scfg::CheckANT  (long lhs, long rhs1, long rhs2, _AVLListX& tempNT, long startSymbol) const {

  
    long avl_index1 = tempNT.Find ((BaseRef)lhs),
         avl_index2 = tempNT.Find ((BaseRef)rhs1),
         avl_index3 = tempNT.Find ((BaseRef)rhs2),
         old_flag1  = tempNT.GetXtra (avl_index1),
         new_flag1  = old_flag1,
         old_flag2  = tempNT.GetXtra (avl_index2),
         new_flag2  = old_flag2,
         old_flag3  = tempNT.GetXtra (avl_index3),
         new_flag3  = old_flag3;

    if  (lhs == startSymbol || (old_flag1 & _HYSCFG_NT_START_)) {
        new_flag1 |= _HYSCFG_NT_START_;
        new_flag2 |= _HYSCFG_NT_START_;
        new_flag3 |= _HYSCFG_NT_START_;
    }

    if  ((old_flag2&_HYSCFG_NT_TERM_) || (old_flag3&_HYSCFG_NT_TERM_)) {
        new_flag1 |= _HYSCFG_NT_TERM_;
    }

    tempNT.SetXtra (avl_index1, new_flag1);
    if (rhs1 != lhs) {
        tempNT.SetXtra (avl_index2, new_flag2);
    }
    if (rhs2 != lhs && rhs2 != rhs1) {
        tempNT.SetXtra (avl_index3, new_flag3);
    }

    return (old_flag1!=new_flag1 || old_flag2!=new_flag2 || old_flag3!=new_flag3);
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void    Scfg::ScanAllVariables  (void) {
    indexInd.Clear();
    indexDep.Clear();
    indexCat.Clear();

    _SimpleList allVariables;
    _AVLList    scannerList(&allVariables);

    for (long formCount = 0; formCount < probabilities.GetHDim(); formCount++) {
        probabilities.GetFormula (formCount,0)->ScanFForVariables(scannerList,true,false,true,true);
    }

    scannerList.ReorderList(); // sort all scanned variables
    
    allVariables.Each ([this] (long var_index, unsigned long) -> void {
        _Variable * parameter = LocateVar (var_index);
        if (parameter->IsCategory()) {
            this->indexCat << var_index;
        } else if (parameter->IsIndependent()) {
            this->indexInd << var_index;
        } else {
            this->indexDep << var_index;
        }
    });
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void    Scfg::VerifyValues  (void) {
    
    _Matrix * sampled_values = (_Matrix*)probabilities.Compute();
    
    sampled_values->ForEachCellNumeric([this] (
        hyFloat value, unsigned long index, unsigned long, unsigned long) -> void {
            if (value < 0.0 || value > 1.0) {
                throw _String ("Probability value for rule ") & _String (GetRuleString (index)) & " is not within [0,1] " & value;
            }
    });
  

  
    byNT2.ForEach([this, sampled_values] (BaseRef object, unsigned long index) -> void {
      hyFloat checksum = 0.;
      auto compute_sum = [&checksum,sampled_values] (long value, unsigned long idx) -> void {
        checksum += sampled_values->directIndex (value);
      };
      ((_SimpleList*)object)->Each (compute_sum);
      ((_SimpleList*)byNT3.GetItem(index))->Each (compute_sum);
      
      if (!CheckEqual (checksum, 1.0)) { // check within reasonable system precision
        throw _String (_String ("Probability values for non-terminal ") & (index+1) & " do not appear to add up to one: " & checksum);
      }
    });
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void    Scfg::RandomSampleVerify  (long samples) {
  try {
    if (samples>0L) {
        long        paramCount = GetIndependentVars().countitems();
        hyFloat  step_factor = 1./samples;

        if (paramCount > 0L) { // some adjustable parameters
            _Matrix parameterBounds (paramCount,2,false,true); // used to store the lower parameter bound
            // the sampling step
            // and the current value
            _Matrix stash;
            GetAllIndependent (stash);
          
            for (long var = 0L; var < paramCount; var++) {
                hyFloat bounds [2] = {GetIthIndependentBound (var, true), GetIthIndependentBound (var)};
                parameterBounds.Store (var, 0, bounds[0]);
                parameterBounds.Store (var, 1, (bounds[1]-bounds[0]) * step_factor);
            }

            _SimpleList zeroThruNm1 (samples-1,0,1);

            for (long it = 0L; it < samples; ++it) {
                zeroThruNm1.Permute (1);
                for (long var = 0L; var < paramCount; var++) {
                    SetIthIndependent (var, parameterBounds(var,0) + parameterBounds(var,1) * zeroThruNm1.list_data[var]);
                }
              
                VerifyValues ();
            }
          
            SetAllIndependent(&stash);
        } else { // check fixed values
            VerifyValues ();
        }
    }
  } catch (const _String & err) {
    HandleApplicationError (err);
  }
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void    Scfg::SetStringCorpus  (_Matrix* stringMatrix) {
  try{
    corpusChar.Clear();
    corpusInt.Clear();
    DumpComputeStructures ();
  
    for (long stringRow = 0L; stringRow < stringMatrix->GetHDim(); stringRow++) {
        for (long stringColumn = 0L; stringColumn < stringMatrix->GetVDim(); stringColumn++) {
            _FString    * aString   = (_FString *)stringMatrix->GetFormula (stringRow,stringColumn)->Compute();
            _SimpleList * tokenized = new _SimpleList;
            TokenizeString (aString->get_str(), *tokenized);
            corpusChar < new _StringBuffer (aString->get_str());
            corpusInt  < tokenized;
        }

    }
    InitComputeStructures();
  } catch (const _String & err) {
    HandleApplicationError (err);
  }
}

/*--------------------------------------------------------------------------------------------------------------------------------*/
void    Scfg::TokenizeString    (_StringBuffer const& inString, _SimpleList& outTokens) const {
    // TODO: SLKP 20180820 Changed the algorithm for tokenization using _Trie calls; confirm that this still works
  
    if (inString.empty()) {
       throw _String ("Empty strings are not allowed as SCFG input.");
    }

    // check to see if the string is too long

    if ((inString.length()+1.0)*inString.length()*0.5*byNT3.countitems() > _HYSCFG_UPPER_BOUND_) {
        throw _String ("The input string is too long.");
    }
  
    unsigned long current_index = 0L;
  
    while (current_index < inString.length()) {
      unsigned long last_index = current_index;
      long next_token_location = terminal_symbols.FindKey (inString, NULL, true, &current_index);
      if (next_token_location < 0L) {
        throw _String ("Invalid terminal symbol in the input string starting with ") & inString.Chop (last_index, current_index).Enquote ();
      }
      outTokens << terminal_symbols.GetValue (next_token_location);
    }
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

Scfg::~Scfg  (void) {
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

_StringBuffer *   Scfg::GetRuleString (long rule_index) const {
    if (rule_index < 0 || rule_index >= rules.countitems()) {
        return new _StringBuffer;
    }

    _StringBuffer * result = new _StringBuffer (64UL);
    _SimpleList * aRule = (_SimpleList*) rules.GetItem (rule_index);
 
    (*result) << "{" << _String(aRule->get(0)) << "}->";
    if (aRule->countitems() == 2) { // NT->T
        (*result) << "\"" << *(_String*)terminals.GetItem (aRule->get (1)) << "\" : ";
    } else {
        (*result) << "{" << _String(aRule->get(1)) << "}{" << _String (aRule->get (2)) << "} : ";
    }

    result->AppendNewInstance((_String*)probabilities.GetFormula (rule_index,0)->toStr(kFormulaStringConversionNormal));
    result->TrimSpace ();
    return result;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

BaseRef     Scfg::toStr  (unsigned long)
{
    _StringBuffer * result = new _StringBuffer (128UL); // allocate a buffer with the initial length of 128 characters
  
    for (long i = 0UL; i < rules.countitems(); i++) {
        result->AppendNewInstance (new _String (GetRuleString(i))) << '\n';
    }

    result->TrimSpace ();
    return result;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void        Scfg::DumpComputeStructures (void) {
    insideProbs.Clear       ();
    insideProbsT.Clear      ();
    outsideProbs.Clear      ();
    outsideProbsT.Clear     ();
    storedInsideP.Clear     ();
    storedOutsideP.Clear    ();
    computeFlagsI.Clear     ();
    computeFlagsO.Clear     ();
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

hyFloat      Scfg::Compute (void) {

    bool first = computeFlagsI.nonempty();

    hyFloat  res = 0.0, ip;
             //useJP = hy_env::EnvVariableGetNumber(kSCFGUseJeffreysPrior, 0.0);

    // update the probability matrix
    probabilities.Compute ();


    for (long stringID = 0L; stringID < corpusChar.countitems(); stringID++) {
        // assuming that production probabilities have changed since last Compute(),
        // retain 0 and 1's but clear buffer of stored 0<Pr<1 values.
        _Matrix * cachedProbs = (_Matrix*)(storedInsideP(stringID));
      
        for (long cid = 0L; cid < cachedProbs->GetHDim(); cid++) { // reset 0,1 values to -1
            cachedProbs->Store (cid,0,-1.);
        }

        /*
        snprintf (buf, sizeof(buf), "\nComputing inside prob for string: %s\n", (const char *) *(_String *)corpusChar(stringID));
        BufferToConsole (buf);
        snprintf (buf, sizeof(buf), "\tstart symbol = %d\n\tfirst = %d\n", startSymbol, first);
        BufferToConsole (buf);
        */

        hyFloat  temp = ComputeInsideProb (0,((_String*)corpusChar(stringID))->length()-1,stringID,startSymbol, first);
 
        if (temp == 0.) {
            ReportWarning (_String("Underflow detected for string ") & stringID & ". Spiking optimizer to avoid this region of parameter space.");
            return (-INFINITY);
        }

        ip = log (temp);
        // ip = log( ComputeInsideProb (0,((_String*)corpusChar(stringID))->sLength-1,stringID,startSymbol, first) );
        // note log-transformation

        res += ip;
        if (first) { // reset the big strings
            computeFlagsI.Populate (computeFlagsI.lLength,0,0);
            // computeFlagsO.Populate (computeFlagsO.lLength,0,0);
        }
    }

    if (first) {
        computeFlagsI.Clear();
        // computeFlagsO.Clear();
    }

    insideCalls = 0;    // reset counter

    return res;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void        Scfg::InitComputeStructures (void) {
    unsigned long maxStringLength = 0UL;
    for (unsigned long stringCount = 0UL; stringCount < corpusInt.countitems(); stringCount++) {
        StoreIfGreater(maxStringLength, ((_SimpleList*)corpusInt.GetItem (stringCount))->countitems());

        insideProbsT < new _SimpleList;
        outsideProbsT < new _SimpleList;

        insideProbs  < new _AVLListX ((_SimpleList*)insideProbsT(stringCount));
        outsideProbs < new _AVLListX ((_SimpleList*)outsideProbsT(stringCount));
      
        storedInsideP < new _Vector;
        storedInsideP < new _Vector;

    }
    maxStringLength = (maxStringLength * (maxStringLength+1) * byNT2.countitems() / 64)+1;
    // SLKP 20180820: TODO is 64 here to represent a hardcoded sizeof (long)?
    computeFlagsI.Populate (maxStringLength,0,0);
    computeFlagsO.Populate (maxStringLength,0,0);
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

inline  long    Scfg::scfgIndexIntoAnArray            (long start,long end,long nt,long stringLength) {
    //return (2*stringLength-start-1)*start/2 + end + nt*(stringLength+1)*stringLength/2;

    return end + ((((stringLength << 1) - start -1L) * start + nt * (stringLength + 1L) * stringLength) >> 1);

    /* Each non terminal is given a triangular array for all substrings of a string, indexed by start, end (0-based, of course),
       and subject to end >= start

              end
              0 1 2 ... stringLength

       start =0 1 2 ... stringLength
              x 1 2 ... stringLength
              x x 2 ... stringLength
                    ...
              x x x ... stringLength

        To find the offest of the (s,e) substring from the 0-th element of this array we compute (L = string length)

        (e-s) + sum_{k=0}{s-1} (L-k) = e-s + s/2 (L + L - s + 1) = e - s + (2L - s + 1) s /2 =

     */



    //return (2*stringLength-start+1)*start/2 + end + nt*(stringLength+1)*stringLength/2;       // AFYP June 21, 2006
}


/*--------------------------------------------------------------------------------------------------------------------------------*/

void    Scfg::setIndexBit             (long start,long end,long nt,long stringLength,_SimpleList& theArray) {
    //long theIndex = (2*stringLength-start-1)*start/2 + end + nt*(stringLength+1)*stringLength/2;
    long array_index = scfgIndexIntoAnArray (start, end, nt, stringLength);
    //theArray.list_data[theIndex/32] |= bitMaskArray.masks[theIndex%32];
    theArray[array_index >> 5] |= bitMaskArray.masks [array_index - (array_index << 5 >> 5)];
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

bool    Scfg::getIndexBit             (long start,long end,long nt,long stringLength,_SimpleList& theArray) {
    //long theIndex = (2*stringLength-start-1)*start/2 + end + nt*(stringLength+1)*stringLength/2;
    long array_index = scfgIndexIntoAnArray (start, end, nt, stringLength);

    /*char str255 [255];
    snprintf (str255, sizeof(str255),"Fetch %d %d %d %d from %x to give %d\n", start, end, nt, stringLength, theArray.list_data[theIndex/32], (theArray.list_data[theIndex/32] & (bitMaskArray.masks[theIndex%32])) > 0);
    BufferToConsole (str255);*/
    //return (theArray.list_data[theIndex/32] & bitMaskArray.masks[theIndex%32]) > 0;
  
    return (theArray.get (array_index >> 5) & bitMaskArray.masks [array_index - (array_index << 5 >> 5)]) > 0L;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/

/* ---- SCFG:  Inside Probability -------------------------- */

hyFloat   Scfg::ComputeInsideProb(long from, long to, long stringIndex, long ntIndex, bool firstPass) {

    // static long insideCalls = 0;
    insideCalls ++;


    /* quickly handle extreme cases */
    if (to > from) { // more than a single terminal in the substring
        if (((_SimpleList*)byNT3.GetItem (ntIndex))->empty()) {
            // ntIndex non-terminal can only generate single terminals
            return 0.;
        }
    } else {    // to == from [single terminal]
        if (((_SimpleList*)byNT2.GetItem (ntIndex))->empty()) {
            // ntIndex non-terminal can only generate other non-terminals
            return 0.;
        }
    }

    /* look up the triple in the cache and decide if it's already been done */

    _AVLListX *             theAVL = (_AVLListX*)insideProbs.GetItem (stringIndex);
    // see if this is already done
  
    long    stringL         = ((_SimpleList*)corpusInt.GetItem (stringIndex))->countitems(),
            tripletIndex    = scfgIndexIntoAnArray(from,to,ntIndex,stringL),
            avlIndex        = theAVL->FindLong (tripletIndex),
            matrixIndex     = -1;

    if (avlIndex < 0) { // the triplet is not in the search tree
        if (firstPass) { // has this been done?
            if (getIndexBit (from, to, ntIndex, stringL, computeFlagsI)) { // already done
                return 0.0;
            }
        } else {
            return 0.0;
        }
    } else { // is in the tree
        matrixIndex = theAVL->GetXtra (avlIndex);
        if (matrixIndex < 0) { // identically 1
            return 1.0;
        } else { // have we already computed this?
            hyFloat current_value = ((_Vector*)storedInsideP.GetItem (stringIndex))->directIndex(matrixIndex);
            if (current_value >= 0.0) {
                return current_value;
            }
            // else was reset to -1, re-compute
        }
    }

    /* have to compute stuff */

    hyFloat      insideProbValue = 0.0;

    if (to == from) { // single terminal substring; direct lookup
        //long ruleIndex = ntToTerminalMap.list_data[indexNT_T(ntIndex,((_SimpleList**)corpusInt.list_data)[stringIndex]->list_data[to])];
        long rule_index = ntToTerminalMap.get (indexNT_T(ntIndex,((_SimpleList*)corpusInt.GetItem (stringIndex))->get(to)));
        if (rule_index >= 0L) {
            insideProbValue = LookUpRuleProbability (rule_index);
        }
    } else {
        if (firstPass) {
            // try the heuristics
            _SimpleList* stringList = (_SimpleList*)corpusInt.GetItem (stringIndex);
            // if any of these conditions are TRUE, then both inside and outside Pr(s,t,i) = 0
            if (firstArray.get(indexNT_T(ntIndex, stringList->get(from))) ||
                    lastArray.get(indexNT_T(ntIndex, stringList->get(to)))
                   || (from == 0L || precursorArray.get(indexNT_T(ntIndex, stringList->get(from-1L))))
                   || (to == stringL-1L || followArray.get(indexNT_T(ntIndex, stringList->get(to+1L))))) {
                setIndexBit (from, to, ntIndex, stringL, computeFlagsI);
                return 0.0;
            }
        }
    }

    ((_SimpleList*)byNT3.GetItem (ntIndex))->Each (
         [this,from, to, &insideProbValue,stringIndex,firstPass] (long rule_index, unsigned long array_index) -> void {
           hyFloat    rule_prob = LookUpRuleProbability(rule_index);
           if (rule_prob > 0.0) {
             _SimpleList * currentRule = (_SimpleList *)rules.GetItem(rule_index);
             long          nt1          = currentRule->get(1),
                           nt2          = currentRule->get(2),
                           halfway      = from+ ((to-from) >> 1)+1;
             
             for (long bp = from+1L; bp <= halfway; bp++){
               hyFloat t = ComputeInsideProb (from,bp-1,stringIndex,nt1,firstPass);
               if (t>0.0) {
                 insideProbValue += t*ComputeInsideProb (bp,to,stringIndex,nt2,firstPass)*rule_prob;
               }
             }
             
             for (long bp = halfway+1L; bp <= to; bp++) { // now loop over all breakpoints
               hyFloat t = ComputeInsideProb (bp,to,stringIndex,nt2,firstPass);
               if (t > 0.0) {
                 
                 insideProbValue += t * ComputeInsideProb (from,bp-1,stringIndex,nt1,firstPass)*rule_prob;
               }
             }
           }
         }
    );

    /* decide what to do with the probability:

    0 - do nothing
    1 - store -1 into the AVL
    2 - store into the buffer matrix and the index into the AVL

    */


    if (insideProbValue > 0.0) {
        /*
        char str255[255];
        snprintf (str255, sizeof(str255), "%d->%d:%d \t p = %g\n", ntIndex, from, to, insideProbValue);
        BufferToConsole (str255);
        */

        if (avlIndex < 0L) {
            long    mxID = -1L;
            if (insideProbValue < 1.0) {
                mxID = ((_Vector*)storedInsideP.GetItem (stringIndex))->Store(insideProbValue);
            }
            theAVL->Insert ((BaseRef)tripletIndex, mxID);
        } else { // update values
            ((_Vector*)storedInsideP.GetItem (stringIndex))->_Matrix::Store (matrixIndex,0,insideProbValue);
        }
    }



    if (firstPass) {
        setIndexBit (from, to, ntIndex, stringL, computeFlagsI);

        /* ******** DEBUGGING ******** */
        /*
        char buf [255];
        snprintf (buf, sizeof(buf), "inside %d\t%d\t%d\t%lf\n", from, to, ntIndex, insideProbValue);
        BufferToConsole(buf);
         */
    }


    return  insideProbValue;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/
void        Scfg::AddSCFGInfo (_AssociativeList* theList) {
  static _String const      // various constants used in AddSCFGInfo
                kAddSCFGInfoStats           ("STATISTICS"),
                kAddSCFGInfoProductions     ("PRODUCTIONS"),
                kAddSCFGInfoTerminals       ("TERMINALS"),
                kAddSCFGInfoProbabilities   ("PROBABILITIES");

  
    _SimpleList indexer (corpusChar.lLength, 0, 1);
    InsertStringListIntoAVL (theList, hy_env::kSCFGCorpus, indexer, corpusChar);
  
    _List       ruleStrings;
    for (long i=0; i<rules.countitems(); i++) {
        ruleStrings.AppendNewInstance(GetRuleString (i));
    }

    indexer.Populate (rules.countitems(), 0, 1);
    InsertStringListIntoAVL (theList, kAddSCFGInfoProductions, indexer, ruleStrings);
    indexer.Populate (terminals.countitems (), 0, 1);
    InsertStringListIntoAVL (theList, kAddSCFGInfoTerminals, indexer, terminals);

    _Matrix * stats = new _Matrix (corpusChar.lLength,3,false,true);

    for (long k=0L; k<corpusChar.countitems(); k++) {
        long       strL    = ((_String*)corpusChar.GetItem (k))->length(),
                   pNot0   = ((_AVLListX*)insideProbs.GetItem(k))->countitems(),
                   p01     = ((_Vector*)storedInsideP.GetItem(k))->get_used();

        hyFloat totalPR = strL*(strL+1.)*0.5*byNT2.countitems();
        stats->Store (k,0,totalPR); // total number of tuples (i,j,v)
        stats->Store (k,1,pNot0-p01);   // number of tuples with probability = 1
        stats->Store (k,2,p01);     // number of tuples with probability in interval (0,1)
    }
    theList->MStore (kAddSCFGInfoStats, stats, false);
    stats = (_Matrix*)probabilities.Compute();
    theList->MStore (kAddSCFGInfoProbabilities, stats, true);
}


/*--------------------------------------------------------------------------------------------------------------------------------*/

hyFloat   Scfg::ComputeOutsideProb(long from, long to, long stringIndex, long ntIndex, bool firstOutside, bool firstInside) {
    outsideCalls++;

  
    /* get length of this string */
    const long            stringL = ((_SimpleList*)corpusInt.GetItem (stringIndex))->countitems();

    // quickly handle extreme cases
    if (from == 0L && to == stringL-1L) { // dealing with entire string
        return ntIndex == 0L ? 1.0 : 0.0; // unity if the current symbol is the start symbol
    }

    if (to > from) {
        if (((_SimpleList*)byNT3.GetItem (ntIndex))->empty()) {
           return 0.;  // NT can only emit terminal; inside is 0, so outside is useless
        }
    } else {    // from==to
        if (((_SimpleList*)byNT2.GetItem (ntIndex))->empty()) {
          return 0.;  // NT can only generate NT NT
        }
    }


    if (!firstInside && firstOutside) { // re-use information from inside probs.
        if (ComputeInsideProb(from,to,stringIndex,ntIndex,firstInside) == 0.) {
            setIndexBit (from,to,ntIndex,stringL,computeFlagsO);
            return 0.;
        }
    }

    /*
    char buf [256];
    snprintf (buf, sizeof(buf), "%d\t%d\t%d\n", from, to, ntIndex);
    BufferToConsole(buf);
    */

    /* look up (s,t,i) in cache and decide if it has already been computed */
    _AVLListX *     theAVL          = (_AVLListX*)outsideProbs.GetItem (stringIndex);

    long            tripletIndex    = scfgIndexIntoAnArray(from,to,ntIndex,stringL),    // convert (s,t,i)
                    avlIndex     = theAVL->FindLong (tripletIndex),              // retrieve key
                    matrixIndex        = -1;

    if (avlIndex < 0L) { // the triplet is not in the search tree
        if (firstOutside) { // has this been computed?
            if (getIndexBit (from, to, ntIndex, stringL, computeFlagsO)) {  // if TRUE, then already done and equals 0
                return 0.;
            }
            // else FALSE, and need to compute
        } else {            // otherwise, absence from AVL indicates Pr = 0
            return 0.;
        }
    } else {    // is in the tree
        matrixIndex = theAVL->GetXtra (avlIndex);   // retrieve index into stored outside probs (0<p<1)
        if (matrixIndex < 0) {  // identically 1
            return 1.0;
        } else { // have we already computed this?
            hyFloat currentValue = ((_Vector*)storedOutsideP.GetItem (stringIndex))->directIndex(matrixIndex);
            if (currentValue >= 0.0) {
                return currentValue;
            }
        }
    }

    /* have to compute stuff; note, nothing below this line should get called after the first pass */
    if (firstOutside) {
        // try the heuristics
 
        _SimpleList* stringList = (_SimpleList*)corpusInt.GetItem (stringIndex);
        // if any of these conditions are TRUE, then both inside and outside Pr(s,t,i) = 0
        if (firstArray.get(indexNT_T(ntIndex, stringList->get(from))) ||
            lastArray.get(indexNT_T(ntIndex, stringList->get(to)))
            || (from == 0L || precursorArray.get(indexNT_T(ntIndex, stringList->get(from-1L))))
            || (to == stringL-1L || followArray.get(indexNT_T(ntIndex, stringList->get(to+1L))))) {
          setIndexBit (from, to, ntIndex, stringL, computeFlagsO);
          return 0.0;
        }

    }

    hyFloat      outsideProbValue = 0.0;
  
  

    /* loop over all j-> k  {ntIndex} productions; */
    ((_SimpleList*)byRightNT2.GetItem (ntIndex))->Each (
       [this, from, to, firstInside,firstOutside, stringIndex, &outsideProbValue] (long rule_index, unsigned long) -> void {
         hyFloat  rule_prob = LookUpRuleProbability(rule_index);
         if (rule_prob > 0.0) {
           _SimpleList *   current_rule = (_SimpleList *)rules.GetItem(rule_index);  // rule = list of two integers
           long            j           = current_rule->get(0),
                           k           = current_rule->get(1);
           
           for (long leftBisect = 0; leftBisect < from; leftBisect++) {
             hyFloat  t = ComputeInsideProb(leftBisect,from-1,stringIndex,k,firstInside);
             if (t > 0.0)
               outsideProbValue += t * ComputeOutsideProb(leftBisect,to,stringIndex,j,firstOutside,firstInside) * rule_prob;
           }
         }
       }
    );
  
  /* loop over all j-> {ntIndex} k productions; */
    ((_SimpleList*)byRightNT1.GetItem (ntIndex))->Each (
           [this, from, to, firstInside,firstOutside, stringIndex, stringL, &outsideProbValue] (long rule_index, unsigned long) -> void {
             hyFloat  rule_prob = LookUpRuleProbability(rule_index);
             if (rule_prob > 0.0) {
               _SimpleList *   current_rule = (_SimpleList *)rules.GetItem(rule_index);  // rule = list of two integers
               long            j           = current_rule->get(0),
                               k           = current_rule->get(2);
               
                for (long rightBisect = to+1; rightBisect < stringL; rightBisect++) {
                 hyFloat  t = ComputeInsideProb(to+1,rightBisect,stringIndex,k,firstInside);
                 if (t > 0.0)
                   outsideProbValue += t * ComputeOutsideProb(from,rightBisect,stringIndex,j,firstOutside,firstInside) * rule_prob;
               }
             }

           }
    );
                                                   


    /* decide what to do with the probability:

            0 - do nothing
            1 - store -1 into the AVL
            2 - store into the buffer matrix and the index into the AVL

    */

    if (outsideProbValue > 0.0) {
        /* DEBUGGING */
        //char str255[255];
        //snprintf (str255, sizeof(str255), "%d->%d:%d \t p = %g\n", ntIndex, from, to, outsideProbValue);
        //BufferToConsole (str255);
        /* --------- */

        if (avlIndex < 0) { // triplet (s,t,i) not in search tree
            long    mxID = -1;
            if (outsideProbValue < 1.0) {
                mxID = ((_Vector*)(storedOutsideP.GetItem(stringIndex)))->Store (outsideProbValue);
            }
            theAVL->Insert((BaseRef)tripletIndex, mxID);    // if Pr = 1.0, stores -1
        } else {    // update values
            ((_Vector*)(storedOutsideP(stringIndex)))->_Matrix::Store (matrixIndex,0,outsideProbValue);
        }
    }

    if (firstOutside) { // flag (s,t,i) as computed
        setIndexBit (from,to,ntIndex,stringL,computeFlagsO);
    }

    return outsideProbValue;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

_Matrix*     Scfg::Optimize (_AssociativeList const* options)  /* created by AFYP, 2006-06-20 */ {
  // SLKP 20180820 : TODO check to see if this is functional
  
  const static _String kSCFGOptimizationMethod  ("SCFG_OPTIMIZATION_METHOD");
  
  if (hy_env::EnvVariableGetNumber(kSCFGOptimizationMethod, 0.0)) {
    // fall back to default optimization; do not use EM
    return _LikelihoodFunction::Optimize(options);
  }
 
    //char    buf [256];
    // debugging buffer




    // populated with zeros for longest string by InitComputeStructures();
    // inside flags cleared by calling Scfg::Compute()
  
    bool    firstInside     = computeFlagsI.nonempty(),
            firstOutside    = computeFlagsO.nonempty();


    // update probability matrix -- converts formulas to numerical values
    probabilities.Compute ();


    // create container for new production rule probability estimates
    //  numerator in column 1, denominator in column 2
  
    long        nRules = probabilities.GetSize();
    _Matrix     nextProbs(nRules, 2, FALSE, TRUE);



    // scan production rule formulas for linking constraints, done once for entire optimization
    _List           links;
  
  // TODO SLKP 20180820: this is quite inefficient; should be optimized unless the whole
  // function is being deprecated
    for (long ruleCount = 0L; ruleCount < nRules; ruleCount++) {
        _Formula *  thisRule = probabilities.GetFormula (ruleCount, 0);

        if ( !(thisRule->IsAConstant()) ) {
            _SimpleList * linked_rules = new _SimpleList;
            *linked_rules << ruleCount;

            for (long nextrule = ruleCount+1L; nextrule < nRules; nextrule++) {
                _Formula *  thatRule = probabilities.GetFormula (nextrule, 0);
                if ( thisRule->EqualFormula (thatRule) ) {
                    *linked_rules << nextrule;
                    //snprintf (buf, sizeof(buf), "linked rule %ld to %ld\n", ruleCount+1, nextrule+1);
                    //BufferToConsole (buf);
                }
            }

            if (linked_rules->countitems() > 1) {
                links < linked_rules;    // append duplicate
            } else {
              delete linked_rules;
            }
        }
    }


    /* calculate current corpus log-likelihood */
    hyFloat  newLk = Compute(),
                oldLk;
  
    long        rep = 0;


    /* ========================= */
    /*  start optimization loop  */
    /* ========================= */

    do {
        oldLk = newLk;
      
         for (long ruleCount = 0; ruleCount < nRules; ruleCount++) {
            nextProbs.Store(ruleCount, 0, 0.);   // zero all entries in container matrix
            nextProbs.Store(ruleCount, 1, 0.);
        }



        // loop over strings in corpus
        for (long stringID = 0; stringID < corpusChar.countitems(); stringID++) {
            long    stringL     = ((_String*)corpusInt.GetItem(stringID))->length(),
                    countNT      = byNT2.countitems();

            hyFloat  denom       = 0.,
                     stringProb = ComputeInsideProb(0, stringL-1, stringID, startSymbol, firstInside);
            // by calculating this, firstInside should be reset to FALSE

            _SimpleList*    tokenized_string = (_SimpleList*)corpusInt.GetItem(stringID);

            /* After computing string Pr., absence in AVL indicates inside Pr = 0   */
            /* so treat firstInside as FALSE below...                               */

            if (stringProb > 0.) {
                for (long ntIndex = 0L; ntIndex < countNT; ntIndex++) {  // loop over all non-terminal symbols

                    /* calculate denominator (number of times i-th NT was used in string */
                    denom = 0.;
                  bool  has_NTNT = ((_SimpleList*)byNT3.GetItem (ntIndex))->nonempty();
                  
                    for (long from = 0; from < stringL; from++) {
                        if (!has_NTNT) { // only i->m
                            denom += ComputeInsideProb(from,from,stringID,ntIndex,FALSE) *
                                     ComputeOutsideProb(from,from,stringID,ntIndex,firstOutside,FALSE) /
                                     stringProb;
                        } else {    // NT has i->jk productions
                            //for (long to = from; to < stringL; to++)
                            for (long to = from + 1L; to < stringL; to++) {  // count only i->jk

                                denom += ComputeInsideProb(from,to,stringID,ntIndex,FALSE) *
                                         ComputeOutsideProb(from,to,stringID,ntIndex,firstOutside,FALSE) /
                                         stringProb;


                                /*  ^^^ It became necesssary to use this calcalation rather than
                                    the short-cut below, because it choked on node-censored grammars. */

                                /*
                                long        tripletIndex    = scfgIndexIntoAnArray(from,to,ntIndex,stringL),
                                avlIndex        = theAVL->Find ((BaseRef)tripletIndex),
                                matrixIndex     = -1;

                                if (avlIndex > -1)  // triplet has non-zero inside Pr
                                {
                                    matrixIndex = theAVL->GetXtra (avlIndex);
                                    if (matrixIndex > -1)   // skip productions with Pr identically 1, don't train.
                                                            //  This also applies to node censoring -- afyp, Aug 30, 2006
                                    {
                                        hyFloat  ip  = ((_GrowingVector**)storedInsideP.list_data)[stringID]->theData[matrixIndex],
                                        op  = ComputeOutsideProb(from,to,stringID,ntIndex,firstOutside,FALSE);
                                        if (op > 0)
                                        {
                                            denom += ip*op/stringProb;
                                            // denom += exp( log(ip) + log(op) - log(stringProb) );
                                        }
                                    }
                                }
                                */
                            }
                        }
                    } // end loop over [s]

 
                    if (denom > 0.) {
                        // loop over i->m productions of current NT
                      
                      ((_SimpleList*)byNT2.GetItem (ntIndex))->Each([this, ntIndex, stringID, firstOutside, stringProb, nRules, denom, &nextProbs] (long rule_index, unsigned long) -> void {
                          hyFloat         numer       = 0.;
                         _SimpleList *    currentRule =  (_SimpleList * )rules.GetItem (rule_index);
                          long            termIndex   = currentRule->get(1);
                          hyFloat         rule_prob    = LookUpRuleProbability (rule_index);
                        
                          if (rule_prob == 1.0 || rule_prob == 0.0) {
                            return;
                          }
                        
                        ((_SimpleList*) corpusInt.GetItem (stringID))->Each ([this, termIndex, stringID, ntIndex, stringProb, firstOutside, &numer] (long token_code, unsigned long from) -> void {
                          if (token_code == termIndex) {
                            hyFloat  tryProb = ComputeInsideProb(from,from,stringID,ntIndex,FALSE);
                            if (tryProb > 0.0) {
                              numer += tryProb *
                              ComputeOutsideProb(from,from,stringID,ntIndex,firstOutside,FALSE) /
                              stringProb;
                            }
                          }
                        });
                        nextProbs[rule_index] += numer;
                        nextProbs[nRules + rule_index] += denom;
                     });
                      
 

                      

                        // loop over all i->jk productions of current NT
                    ((_SimpleList*)byNT3.GetItem (ntIndex))->Each ([this, nRules, stringL, ntIndex, tokenized_string, stringID, firstOutside, stringProb, denom, &nextProbs] (long rule_index, unsigned long) {
                        hyFloat      rule_prob    = LookUpRuleProbability (rule_index),
                                    numer        = 0.;
                      
                      if (rule_prob == 0. || rule_prob == 1.) {
                        return;
                      }
                      
                      _SimpleList *   currentRule = (_SimpleList*)rules.GetItem (rule_index);
                      long            jIndex      = currentRule->get(1),
                                      kIndex      = currentRule->get(2);

                      for (long from = 0L; from < stringL-1; from++) { /* calculate equation (20) */
                        // use heuristics to skip zero productions
                        
                        if ( ! (firstArray.get(indexNT_T(ntIndex,tokenized_string->get(from)))
                                && (from == 0 || precursorArray.get(indexNT_T(ntIndex,tokenized_string->get(from-1))))) ) {
                          continue;
                        }
                        
                        
                        for (long to = from; to < stringL; to++) {  /* BUG? <- start at s=t or s=t+1 ??? */
                          // more short-cuts from heuristics
                          
                          if (! (lastArray.get (indexNT_T(ntIndex, tokenized_string->get(to)))
                                 && (to == stringL-1 || followArray.get(indexNT_T(ntIndex,tokenized_string->get(to+1)))))) {
                            continue;
                          }
                          
                
                          hyFloat  op  = ComputeOutsideProb(from,to,stringID,ntIndex,firstOutside,FALSE);
                          
                          if (op == 0.) {
                            continue;    // yet another short-cut
                          }
                          
                          
                          for (long bisect = from; bisect < to; bisect++) {
                            hyFloat ip = ComputeInsideProb(from,bisect,stringID,jIndex,FALSE);
                            if (ip > 0.) {
                              numer += rule_prob * ip *
                              ComputeInsideProb(bisect+1,to,stringID,kIndex,FALSE) * op /
                              stringProb;
                            }
                          }
                        }
                      }
                      
                      /*
                       snprintf (buf, sizeof(buf), "storing %f/%f into ruleIndex %d\n", numer, denom, ruleIndex);
                       BufferToConsole(buf);
                       */
                      
                      nextProbs.theData[rule_index] += numer;
                      nextProbs.theData[nRules + rule_index] += denom;
                    });

                        /*
                         snprintf (buf, sizeof(buf), "nt %d has %d ijk rules\n", ntIndex, ijkRules->lLength);
                         BufferToConsole(buf);
                         */

                    } // endif denom > 0


                } // end loop over NT


            }

            // repopulate with zeros for next string
            if (firstInside) {
                computeFlagsI.Populate(computeFlagsI.countitems(),0,0);
            }
            if (firstOutside) {
                computeFlagsO.Populate(computeFlagsO.countitems(),0,0);
            }



            // clear stored inside and outside probabilities -- need to re-estimate and re-fill AVLs
          
          _Matrix * reset[2] = {(_Matrix*)(storedInsideP.GetItem(stringID)), (_Matrix*)(storedOutsideP.GetItem(stringID))};
          for (_Matrix * mp : reset) {
              for (long cid = 0L; cid < mp->GetHDim(); cid++) {
                mp->Store(cid,0,-1.);
              }
          }
        } // end loop over strings



        // search through rules and sum linked probabilities

        for (long linkIndex = 0; linkIndex < links.lLength; linkIndex++) {
            _SimpleList *   thisLink    = (_SimpleList*)links.list_data[linkIndex];
            hyFloat      linkNumer   = 0.,
                            linkDenom = 0.;
            /*
            snprintf (buf, sizeof(buf), "Link %d contains ", linkIndex);
            BufferToConsole (buf);
             */
            {
                for (long lcount = 0; lcount < thisLink->lLength; lcount++) {
                    /*
                    snprintf (buf, sizeof(buf), "%d(%3.3f/%3.3f) ", thisLink->list_data[lcount], nextProbs.theData[thisLink->list_data[lcount]],
                            nextProbs.theData[nRules+thisLink->list_data[lcount]]);
                    BufferToConsole (buf);
                     */
                    linkNumer += nextProbs.theData[thisLink->list_data[lcount]];
                    linkDenom += nextProbs.theData[nRules+thisLink->list_data[lcount]];
                }
            }

            for (long lcount = 0; lcount < thisLink->lLength; lcount++) {
                // update linked probabilities
                nextProbs.theData[thisLink->list_data[lcount]] = linkNumer / thisLink->lLength;
                nextProbs.theData[nRules+thisLink->list_data[lcount]] = linkDenom / thisLink->lLength;
                /*
                 snprintf (buf, sizeof(buf), "\t all set to %3.3f/%3.3f\n", linkNumer, linkDenom);
                 BufferToConsole(buf);
                 */
            }
        }


        // update rule probabilities with new estimates
        {
            for (long ruleCount = 0; ruleCount < probabilities.GetHDim(); ruleCount++) {
                _Formula *      thisFormula = probabilities.GetFormula (ruleCount,0);
                if (thisFormula->IsAConstant()) {
                    continue;
                }

                _SimpleList     itsVariables;
                _AVLList        scannerList(&itsVariables);

                thisFormula->ScanFForVariables(scannerList,TRUE,FALSE,TRUE,TRUE);
                // for the time being, assume that the first variable in the formula is the one to adjust
                _Variable *     theParameter    = LocateVar (itsVariables.list_data[0]);
                _Constant       nextValue (nextProbs.theData[ruleCount] / nextProbs.theData[nRules+ruleCount]);


                if (nextProbs.theData[nRules+ruleCount] > 0.0) {    // skip productions with denominators = 0
                    // Temporary fix -- AFYP, Sept 27, 2006
                    /*
                    snprintf (buf, sizeof(buf), "set rule %d from %f to %f\n", ruleCount, LookUpRuleProbability(ruleCount),
                         nextProbs.theData[ruleCount] / nextProbs.theData[nRules+ruleCount]);
                    BufferToConsole(buf);
                     */
                    theParameter->SetValue (nextValue.Compute(),true,true,NULL);
                }

            }
        }

        if (firstInside) {
            firstInside = false;    // reset flags
        }
        if (firstOutside) {
            firstOutside = false;
        }

        probabilities.Compute ();   // convert updated formulas to numerical values

        // calculated log-likelihood of corpus using updated production Pr
        // and compare to previous value.
        newLk = Compute();
        ++rep;

        /*
        snprintf (buf, sizeof(buf), "Iteration %d, Log L = %f\n", rep, newLk);
        BufferToConsole (buf);
        */

        if (rep > MAX_ITERATIONS_OPTIMIZE) {
            break;
        }

        /*if (newLk < oldLk) {
            snprintf (buf, sizeof(buf), "\nERROR: log L = %f is lower than previous value %f at step %ld \n", newLk, oldLk, rep);
            BufferToConsole (buf);
            // break;
        }*/

    } while (fabs(newLk - oldLk) > SCFG_OPTIMIZATION_THRESHOLD);    // should be non-absolute - afyp 11/30/2006

    //snprintf (buf, sizeof(buf), "Used %ld iterations in expectation maximization.\n", rep);
    //BufferToConsole (buf);



    // SLKP:  TBI this actually needs to be populated still
    /* return new _Matrix (1,1,false,true); */
    return NULL;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

_String* Scfg::SpawnRandomString(long ntIndex, _SimpleList* storageString) {
  try {
    if (ntIndex < 0L) {
        storageString = new _SimpleList;
        SpawnRandomString (startSymbol, storageString);
        _StringBuffer      *backString = new _StringBuffer (storageString->countitems());
        storageString->Each ([&backString, this] (long terminal_index, unsigned long) -> void {
          (*backString) << *(_String*)terminals.GetItem(terminal_index);
        });
      
        delete storageString;
        backString->TrimSpace();
        return backString;
    }

    hyFloat      randomValue = genrand_real2 (),
                 sum         = 0.;

    
    long            ruleIndex   = 0;
    
    _SimpleList*    aList       = (_SimpleList*)byNT2(ntIndex),
               *    aRule;
    
    // loop through terminal rules (X->x) and sum probabilities
    for (; ruleIndex<aList->lLength && sum < randomValue ; ruleIndex++) {
      sum += LookUpRuleProbability (aList->list_data[ruleIndex]);
    }
    
    if (sum >= randomValue) {
      aRule = (_SimpleList*)rules(aList->get(ruleIndex-1));
      (*storageString) << aRule->list_data[1];
      return nil;
    }
    
    ruleIndex = 0;
    aList       = (_SimpleList*)byNT3(ntIndex);
    for (; ruleIndex<aList->lLength && sum < randomValue ; ruleIndex++) {
      sum += LookUpRuleProbability (aList->list_data[ruleIndex]);
    }
    
    if (sum >= randomValue) {
      aRule = (_SimpleList*)rules(aList->get(ruleIndex-1));
      SpawnRandomString (aRule->get(1),storageString);
      SpawnRandomString (aRule->get(2),storageString);
    } else {
      throw _String ("SCFG::SpawnRandomString() randomValue ") & randomValue & " exceeded sum " & sum & ": nt=" & ntIndex & " stor=" & *(_String *) storageString->toStr();
    }
    
   
  } catch (const _String & err) {
    HandleApplicationError(err);
    return new _String;
  }
  
  return nil;

}


/*--------------------------------------------------------------------------------------------------------------------------------*/
/*  Implementation of CYK algorithm by AFYP     2006-07-12                                                                        */
/*--------------------------------------------------------------------------------------------------------------------------------*/
_StringBuffer *   Scfg::BestParseTree(void) {
    long            countNT         = byNT2.countitems();
    bool            firstPass       = computeFlagsI.nonempty();

    _StringBuffer *       parseTreeString = new _StringBuffer();
    //char          buf [4096];

    for (long stringIndex = 0L; stringIndex < corpusInt.countitems(); stringIndex++) {
        long    stringL     = ((_SimpleList*)corpusInt.GetItem(stringIndex))->countitems();

        // initialize AVL tree for storing non-zero subtrees (i,j,v) to refer to _SimpleList
        _SimpleList     triplets;
        _AVLListX *     theAVL = new _AVLListX (&triplets);

        _SimpleList     argMaxYZK;      // stores (y,z,k) keyed by argmax(i,j,v) in AVL
        _Vector         *theMatrix = new _Vector;     // stores likelihood for subtree (i,j,v)

        for (long from = 0; from < stringL; from++) {   // initialization
            for (long ntIndex = 0; ntIndex < countNT; ntIndex++) {
                long        tripletIndex    = scfgIndexIntoAnArray (from,from,ntIndex,stringL),
                            mxID         = -1;
                hyFloat  lk = ComputeInsideProb(from,from,stringIndex,ntIndex,firstPass);

                if (lk > 0.) {
                    mxID = theMatrix->Store (lk);
                    theAVL->Insert ((BaseRef)tripletIndex, mxID);
                    argMaxYZK.AppendRange (3,0,0);
                    for (long count = 0; count < 3; count++) {
                        argMaxYZK << 0;    // append (0,0,0)
                    }
                }
            }
        }

        for (long from = 0; from < stringL-1; from++) { // iterate over all substrings and non-terminals
            for (long to = from+1L; to < stringL; to++) {
                for (long ntIndex = 0L; ntIndex < countNT; ntIndex++) {
                    hyFloat      maxLk = 0.;
                  
                    long            maxLeft    = -1L,
                                    maxRight   = -1L,
                                    maxBisect  = -1L;
                  
                    _SimpleList *   itsRules = (_SimpleList*)byNT3.GetItem (ntIndex);

                    for (unsigned long ruleIdx = 0UL; ruleIdx < itsRules->countitems(); ruleIdx++) {    // iterate over all productions
                        long            currentRuleIndex    = itsRules->get(ruleIdx);
                        _SimpleList *   currentRule         = (_SimpleList*)rules.GetItem (currentRuleIndex);
                        hyFloat         ruleProb            = LookUpRuleProbability(currentRuleIndex);
                        long            leftNT              = currentRule->get(1),
                                        rightNT             = currentRule->get(2);

                        if (ruleProb > 0.) {
                            for (long bisect = from; bisect < to; bisect++) {       // iterate over all bisects of substring
                                hyFloat  tryProb = ComputeInsideProb(from,bisect,stringIndex,leftNT,firstPass);
                                if (tryProb > 0.) {
                                    hyFloat  lk  = ruleProb * tryProb *
                                                      ComputeInsideProb(bisect+1,to,stringIndex,rightNT,firstPass);
                                    if (lk > maxLk) {
                                        maxLk = lk;
                                        maxLeft = leftNT;
                                        maxRight = rightNT;
                                        maxBisect = bisect;
                                    }
                                }
                            }
                        }
                    }

                    long    tripletIndex    = scfgIndexIntoAnArray (from,to,ntIndex,stringL),
                            mxID         = -1,
                            insertFlag;

                    if (maxLk > 0.) {
                        mxID = theMatrix->Store (maxLk);    // store most likely production and bisect for triplet

                        // snprintf (buf, sizeof(buf), "stored triplet into matrix ID %d\n", mxID);
                        // BufferToConsole(buf);

                        insertFlag = theAVL->Insert ((BaseRef)tripletIndex, mxID);
                        if (insertFlag > -1) {  // store (y,z,k)
                            argMaxYZK << maxLeft;
                            argMaxYZK << maxRight;
                            argMaxYZK << maxBisect;
                            /*
                             snprintf (buf, sizeof(buf), "(%d,%d,%d) stored (%d,%d,%d) with L=%f\n", from,to,ntIndex, maxLeft,maxRight,maxBisect, maxLk);
                             BufferToConsole(buf);
                             */
                        }
                    }
                }
            }
        }

        // apply stack to reconstructing best parse
        CykTraceback (0,stringL-1,0,stringIndex,theAVL,&argMaxYZK,theMatrix,*parseTreeString);

        (*parseTreeString) = (*parseTreeString) & "\n"; // separate tree strings
    }


    // BufferToConsole((char *)(const char *) (*parseTreeString));
    parseTreeString -> TrimSpace();
    return (parseTreeString);
}


/*------------------------------------------------------------------------------------------------------------------------------*/
/*  Converts stack of best-parse productions into a tree string.- AFYP                                                          */
void    Scfg::CykTraceback (long i, long j, long v, long stringIndex, _AVLListX const * theAVL, _SimpleList const * theYZKs, _Vector const * theMatrix, _StringBuffer& parseTreeString) const {
    //  for sub-string (i,j) and non-terminal (v), string of length (stringL)
    _SimpleList * corpusString  = ((_SimpleList*)corpusInt.GetItem (stringIndex));

    long    stringL         = corpusString->countitems(),
            tripletIndex    = scfgIndexIntoAnArray(i,j,v,stringL),
            avlIndex        = theAVL -> Find ((BaseRef)tripletIndex);


    if (avlIndex < 0) { // value not in tree
        ReportWarning (_String("ERROR: Unknown triplet encountered in CYK traceback: (") & i & "," & j & "," & v & ")");
    } else {
        long        matrixIndex = theAVL->GetXtra (avlIndex) * 3,
                    y           = theYZKs->get(matrixIndex),
                    z           = theYZKs->get(matrixIndex + 1),
                    k           = theYZKs->get(matrixIndex + 2);

        if (y==0 && z==0 && k==0) { // node terminates
            parseTreeString << '(' << v << "terminal " << corpusString->get (i) << ')';
            //(*parseTreeString) = (*parseTreeString) & "(" & v & " " & corpusString->sData[i] & ")";
            /*
            snprintf (buf, sizeof(buf), "(%d %s)", v, corpusString->sData[i]);
            (*parseString) << (const char *) buf;
             */
        } else {
            parseTreeString << "(" << v << " ";
            /*
            snprintf (buf, sizeof(buf), "(%d ", v);
            (*parseString) << (const char *) buf;
            */

            CykTraceback(i,k,y,stringIndex,theAVL,theYZKs,theMatrix,parseTreeString);

            //(*parseTreeString) = (*parseTreeString) & " ";

            CykTraceback(k+1,j,z,stringIndex,theAVL,theYZKs,theMatrix,parseTreeString);

            parseTreeString << ")";

            /*
            snprintf (buf, sizeof(buf), ")");
            (*parseString) << (const char *) buf;
             */
        }
    }
}

#endif

/*--------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------*/



