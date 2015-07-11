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

#include    "scfg.h"

#ifdef      _USE_HYPHY_HOOKS_

#include    "math.h"

_String     _HYSCFG_TERM_KEY_T  ("T"),
            _HYSCFG_KEY_P       ("P"),
            _HYSCFG_KEY_L       ("L"),
            _HYSCFG_NT_KEY_1    ("1"),
            _HYSCFG_NT_KEY_2    ("2");

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


_String     // various constants used in AddSCFGInfo
_addSCFGInfoStats           ("STATISTICS"),
                            _addSCFGInfoProductions     ("PRODUCTIONS"),
                            _addSCFGInfoTerminals       ("TERMINALS"),
                            _addSCFGInfoProbabilities   ("PROBABILITIES"),

                            useJeffreysPrior            ("USE_JEFFREYS_PRIOR"),     // added Nov. 28, 2006 by afyp
                            scfgOptimizationMethod      ("SCFG_OPTIMIZATION_METHOD");

extern _String  scfgCorpus;


#ifdef      __MACPROFILE__
#include "profiler.h"
#endif

#else

#define     MAX_WARNINGS        3
#define     MAX_LINE_LEN        10000

#endif

/*          ================                                                        */
/* ========| SCFG FUNCTIONS |====================================================== */
/*          ================                                                        */

#ifdef                  _USE_HYPHY_HOOKS_

Scfg::Scfg  (_AssociativeList* T_Rules,  _AssociativeList* NT_Rules, long ss)
{
    _String         errorMessage;
    _SimpleList     parsedFormulas;
    // stash pointers to processed formulas here

    startSymbol     = ss;
    insideCalls     = 0;
    outsideCalls    = 0;

    // initialize the pointers
    parseTree     = nil;

    long          termRules     = T_Rules->avl.countitems(),
                  nonTermRules  = NT_Rules->avl.countitems();

    if (termRules == 0) {
        errorMessage = "A SCFG can not be constructed from an empty set of <Nonterminal>-><Terminal> production rules";
    } else {
        _List         ruleProbabilities,            // an auxiliary list used to store strings describing production probabilities
                      alreadySeenX;                 // a list of production rules that have already been added
        // it is used to keep track of duplicate production rules

        _SimpleList   foundNT,                      // an array to keep track of all 'declared' non-terminals
                      ntFlags;                      // an array of flags of non-terminal properties

        _AVLList      tempTerminals (&terminals),       // an auxiliary wrapper around the set of terminal symbols used to check for duplication
                      alreadySeen   (&alreadySeenX);    // and a wrapper for the alreadySeen list

        _AVLListX     tempNT        (&foundNT);         // and a wrapper for the set of non-terminals

        // build a prefix parse tree as we go along
        // begin by allocating memory for the root data structure
        checkPointer (parseTree = new node<long>* [256]);

        for (long it = 0; it < 256; it++) {
            parseTree [it] = (node<long>*)nil;
        }

        for (long tc = 0; tc < termRules; tc++)
            // check Terminal rules for consistency
            // traverse the T_Rules AVL, which is assumed to be indexed by 0..termRules-1
        {
            _AssociativeList * aRule = (_AssociativeList*)T_Rules->GetByKey (tc, ASSOCIATIVE_LIST);
            if (aRule) {
                _FString    * literal       = (_FString*)aRule->GetByKey (_HYSCFG_TERM_KEY_T,STRING),
                              * expression  = (_FString*)aRule->GetByKey (_HYSCFG_KEY_P, STRING);

                _Constant   * lhs           = (_Constant*)aRule->GetByKey (_HYSCFG_KEY_L, NUMBER);

                if (literal && lhs && literal->theString->sLength) {
                    long index = tempTerminals.Insert (literal->theString);
                    // insert to the list of terminals if not seen before
                    if  (index>=0) // new terminal; added to list
                        // now we process the terminal into the parse tree
                    {
                        literal->theString->nInstances++; // increase reference counter for the string object
                        // add    the literal to the parse tree
                        // handle the first character separately

                        char        currentCharacter = literal->theString->getChar(0);
                        node<long>* currentTreeNode  = parseTree[currentCharacter];

                        bool        addedRootStub    = false;

                        if (currentTreeNode == nil) {
                            checkPointer (currentTreeNode = new node<long>);
                            currentTreeNode->init(0);
                            parseTree[currentCharacter] = currentTreeNode;
                            addedRootStub = true;
                        }

                        long charP = 1;

                        for  (; charP < literal->theString->sLength; charP++) {
                            currentCharacter        = literal->theString->getChar(charP);

                            long   availableNodes   = currentTreeNode->get_num_nodes (),
                                   nodeCounter      = (availableNodes>0);

                            // SLKP: 20100630 need to check for a one-character existing prefix

                            if (availableNodes)
                                for (; nodeCounter<=availableNodes; nodeCounter++) {
                                    node <long> * tryANode = currentTreeNode->go_down (nodeCounter);
                                    if ((tryANode->get_data () & _HYSCFG_CHARACTER_MASK_) == currentCharacter)
                                        // can spell the 'literal' string down this branch
                                    {
                                        if (tryANode->get_num_nodes () == 0) // ERROR: not a prefix code; the terminal currently being added
                                            // contains another terminal as a prefix
                                        {
                                            availableNodes = 0;
                                        } else {
                                            currentTreeNode = tryANode;
                                        }
                                        break;
                                    }
                                }

                            if (availableNodes || (addedRootStub && charP == 1)) // no error set
                                // SLKP: 20100630
                                // looks like there was a bug here, when a partial prefix would
                                // not simply be reused, but rather re-added to the parent node
                            {
                                if (availableNodes == nodeCounter)
                                    // and no matching child node has been found
                                {
                                    // insert the new child
                                    node<long> *    addANode = (node<long>*) checkPointer(new node<long>);
                                    addANode->init ((long)currentCharacter);
                                    currentTreeNode->add_node (*addANode);
                                    currentTreeNode = addANode;
                                }
                            } else {
                                errorMessage = _String("Terminal symbols must form a prefix code, but '") & *literal->theString & "' contains another terminal symbol as a prefix.";
                                break;
                            }
                        }

                        if (errorMessage.sLength == 0) { // no error set; check to see if this terminal is not a prefix of something else
                            if (currentTreeNode->get_num_nodes () != 0) {
                                errorMessage = _String("Terminal symbols must form a prefix code. ") & *literal->theString & " is a prefix of another terminal.";
                            } else {
                                currentTreeNode->init (currentTreeNode->get_data() | (index << 8));
                            }
                        }
                    } else {
                        index = -index-1;    // if this terminal has already been added, use the index
                    }

                    if (errorMessage.sLength == 0)
                        // a valid production rule
                    {
                        long          nt_index  = (long)(lhs->Compute()->Value()),
                                      avl_index = tempNT.Insert ((BaseRef)nt_index); // store the integer index of the LHS if needed

                        if (avl_index<0) { // nt_index already exists; correct to positive range
                            avl_index = -avl_index - 1;
                        }

                        // update status flags for this non-terminal
                        tempNT.SetXtra (avl_index, tempNT.GetXtra (avl_index)|_HYSCFG_NT_LHS_|_HYSCFG_NT_DTERM_|_HYSCFG_NT_TERM_);


                        // first ensure the rule is not a duplicate
                        _String         ruleString = _String (nt_index) & ",[" & index & ']';
                        long            seenMe     = alreadySeen.Insert (ruleString.makeDynamic());
                        if (seenMe < 0) { // already seen
                            errorMessage = _String ("Duplicate production rule:" ) & GetRuleString (-seenMe-1);
                        } else {
                            // create a new record for the rule of the form [nt index] -> [t index]
                            _SimpleList *goodTRule = (_SimpleList*) checkPointer (new _SimpleList ((long)nt_index));
                            (*goodTRule) << index;
                            rules.AppendNewInstance (goodTRule); // append the new rule to the list of existing rules

                            // process the formula
                            ProcessAFormula (expression, ruleProbabilities, parsedFormulas, errorMessage);
                        }

                        if (errorMessage.sLength == 0) {
                            continue;    // good rule! go on to check the next one
                        }
                    }
                } else {
                    errorMessage = _String("Each terminal rule must have a non-empty target terminal and a left-hand side non-terminal. Rule ") & tc & " did not comply.";
                }
            } else {
                errorMessage = _String("Each rule must be specified as an associative array, but rule ") & tc & " was not.";
            }
            break;
        } // done checking terminal rules

        if (errorMessage.sLength == 0) { // all terminal rules were good; now we can check the non-terminal rules
            for (long tc = 0; tc < nonTermRules; tc++)
                // check Non-terminal rules for consistency
                // traverse the NT_Rules AVL, which is assumed to be indexed by 0..nonTermRules-1
            {
                _AssociativeList * aRule = (_AssociativeList*)NT_Rules->GetByKey (tc, ASSOCIATIVE_LIST);
                if (aRule) {
                    _FString    * expression    = (_FString*)aRule->GetByKey     (_HYSCFG_KEY_P, STRING);

                    _Constant   * lhs           = (_Constant*)aRule->GetByKey    (_HYSCFG_KEY_L, NUMBER),
                                  * rhs1            = (_Constant*)aRule->GetByKey    (_HYSCFG_NT_KEY_1, NUMBER),
                                    * rhs2           = (_Constant*)aRule->GetByKey    (_HYSCFG_NT_KEY_2, NUMBER);

                    if (lhs && rhs1 && rhs2) {
                        ProcessAFormula (expression, ruleProbabilities, parsedFormulas, errorMessage);
                        if (errorMessage.sLength == 0) {
                            _SimpleList goodNTRule;
                            long          nt_index = (long)(lhs->Compute()->Value());
                            long avl_index = tempNT.Insert ((BaseRef)nt_index); // store the integer index of the LHS if needed
                            if (avl_index<0) {
                                avl_index = -avl_index - 1;
                            }
                            tempNT.SetXtra (avl_index, tempNT.GetXtra (avl_index)|_HYSCFG_NT_LHS_); // update status flags for this non-terminal
                            _String     ruleString = _String (nt_index) & ',';
                            goodNTRule << nt_index;
                            nt_index    = (long)(rhs1->Compute()->Value());
                            tempNT.Insert ((BaseRef)nt_index);
                            goodNTRule << nt_index;
                            ruleString = ruleString & _String (nt_index) & ',';
                            nt_index    = (long)(rhs2->Compute()->Value());
                            tempNT.Insert ((BaseRef)nt_index);
                            goodNTRule << nt_index;
                            ruleString = ruleString & _String (nt_index);
                            long            seenMe     = alreadySeen.Insert (ruleString.makeDynamic());
                            if (seenMe < 0) { // already seen
                                //errorMessage = _String ("Duplicate production rule:" ) & GetRuleString (-seenMe-1);
                                errorMessage = _String ("Duplicate production rule: ") & tc & " : " & ruleString;
                            } else {
                                rules     && & goodNTRule; // append the new rule to the list of existing rules
                                continue;
                            }
                        }
                    } else {
                        errorMessage = _String("Each non-terminal rule must have two right-hand side non-terminals and a left-hand side non-terminal. Rule ") & tc & " did not comply.";
                    }
                } else {
                    errorMessage = _String("Each rule must be specified as an associative array, but rule ") & tc & " was not.";
                }
                break;
            }
        }

        if (errorMessage.sLength == 0) { // all rules were good; next steps:
            // (a). Validate the grammar
            // (1).   Each non-terminal must appear in at least one LHS
            // (2).   Each non-terminal must be resolvable to a terminal at some point
            // (3).   Each non-terminal must be reachable from the start symbol

            ntToTerminalMap.Populate (rules.lLength*terminals.lLength,-1,0);

            bool         continueLoops = true;

            // prepopulate the list of rules stratified by the LHS non-terminal by empty lists
            for (long ntC = 0; ntC < foundNT.lLength; ntC ++) {
                _SimpleList emptyList;
                byNT3 && & emptyList;
                byNT2 && & emptyList;

                byRightNT1 && & emptyList;  // addition by AFYP, 2006-06-20
                byRightNT2 && & emptyList;
            }

            for (long ruleIdx = 0; ruleIdx < rules.lLength; ruleIdx ++) {
                _SimpleList *aList = (_SimpleList*)rules(ruleIdx);      // retrieve two- or three-integer list
                if (aList->lLength == 3) { // NT->NT NT
                    *((_SimpleList*)byNT3 (aList->lData[0])) << ruleIdx;
                    /* addition by AFYP, 2006-06-20 */
                    *((_SimpleList*)byRightNT1 (aList->lData[1])) << ruleIdx;
                    *((_SimpleList*)byRightNT2 (aList->lData[2])) << ruleIdx;
                    /* ---------- end add --------- */
                } else {
                    *((_SimpleList*)byNT2 (aList->lData[0])) << ruleIdx;
                    ntToTerminalMap.lData [indexNT_T (aList->lData[0],aList->lData[1])] = ruleIdx;
                }

            }

            while (continueLoops) { // set status flags for all NT symbols based on production rules
                continueLoops = false;
                for (long ruleIdx = 0; ruleIdx < rules.lLength; ruleIdx ++) {
                    _SimpleList *aList = (_SimpleList*)rules(ruleIdx);
                    if (aList->lLength == 3) { // NT->NT NT
                        continueLoops = continueLoops || CheckANT (aList->lData[0],aList->lData[1],aList->lData[2], tempNT, startSymbol);
                    }
                }
            }


            // now iterate over the list of declared NT and verify that they all comply to the three conditions above
            {
                for (long ntC = 0; ntC < foundNT.lLength; ntC ++) {
                    long ntFlag = tempNT.GetXtra (ntC);
                    if ((ntFlag & _HYSCFG_NT_LHS_) == 0) {
                        errorMessage = _String ("Non-terminal symbol ") & foundNT.lData[ntC] & " does not appear on the left-hand side of any production rules.";
                        break;
                    }
                    if ((ntFlag & _HYSCFG_NT_START_) == 0) {
                        errorMessage = _String ("Non-terminal symbol ") & foundNT.lData[ntC] & " can not be derived from the start symbol.";
                        break;
                    }
                    if ((ntFlag & _HYSCFG_NT_TERM_) == 0) {
                        errorMessage = _String ("Non-terminal symbol ") & foundNT.lData[ntC] & " can not be used to derive any terminal symbols.";
                        break;
                    }
                }
            }

            if (errorMessage.sLength == 0) // all non-terminals checked out
                // populate the matrix of formulas
            {
                CreateMatrix (&probabilities, parsedFormulas.lLength, 1, false, true, false);
                probabilities.Convert2Formulas ();
                for (long formC = 0; formC < parsedFormulas.lLength; formC++) {
                    probabilities.StoreFormula (formC,0,* ((_Formula**)parsedFormulas.lData)[formC]);
                }

                long countT  = terminals.lLength,
                     countNT = byNT2.lLength;

                // populate firstArray
                firstArray.Populate (countNT*countT,0,0);
                lastArray.Populate (countNT*countT,0,0);
                precursorArray.Populate (countNT*countT,0,0);
                followArray.Populate (countNT*countT,0,0);

                for (long i = 0; i < countNT; i++) {
                    _SimpleList * myRules = ((_SimpleList**)byNT2.lData)[i];
                    for (long i2 = 0; i2 < myRules->lLength; i2++) {    // for all i->m productions
                        long flatIndex = indexNT_T (i,((_List**)rules.lData)[myRules->lData[i2]]->lData[1]);
                        firstArray.lData[flatIndex] = 1;
                        lastArray.lData [flatIndex] = 1;
                    }
                }

                continueLoops = true;
                while (continueLoops) {
                    continueLoops = false;
                    for (long i = 0; i < countNT; i++) {
                        _SimpleList * myRules = ((_SimpleList**)byNT3.lData)[i];
                        for (long i2 = 0; i2 < myRules->lLength; i2++) {
                            long rhs1 = ((_List**)rules.lData)[myRules->lData[i2]]->lData[1],
                                 rhs2 = ((_List**)rules.lData)[myRules->lData[i2]]->lData[2];

                            for (long i3 = 0; i3 < countT; i3++) {
                                long anIndex = indexNT_T (i,i3),
                                     rhs1v = firstArray.lData[indexNT_T (rhs1,i3)],
                                     lhs1v = firstArray.lData[anIndex],
                                     rhs2v = lastArray.lData [indexNT_T (rhs2,i3)],
                                     lhs2v = lastArray.lData [anIndex];

                                if (lhs1v == 0 && rhs1v) {
                                    continueLoops = true;
                                    firstArray.lData[anIndex] = 1;
                                }
                                if (lhs2v == 0 && rhs2v) {
                                    continueLoops = true;
                                    lastArray.lData[anIndex] = 1;
                                }
                            }
                        }
                    }
                }
                {
                    for (long i = 0; i < countNT; i++) { // initialize Precursor and Follow
                        _SimpleList * myRules = ((_SimpleList**)byNT3.lData)[i];
                        for (long i2 = 0; i2 < myRules->lLength; i2++) {
                            long rhs1 = ((_List**)rules.lData)[myRules->lData[i2]]->lData[1],
                                 rhs2 = ((_List**)rules.lData)[myRules->lData[i2]]->lData[2];

                            for (long i3 = 0; i3 < countT; i3++) {
                                long idx1 = indexNT_T (rhs1,i3),
                                     idx2 = indexNT_T (rhs2,i3);

                                if (lastArray.lData [idx1]) {
                                    precursorArray.lData[idx2] = 1;
                                }
                                if (firstArray.lData[idx2]) {
                                    followArray.lData[idx1]    = 1;
                                }
                            }
                        }
                    }
                }

                continueLoops = true;
                while (continueLoops) { // populate Precursor and Follow
                    continueLoops = false;
                    for (long i = 0; i < countNT; i++) {
                        _SimpleList * myRules = ((_SimpleList**)byNT3.lData)[i];
                        for (long i2 = 0; i2 < myRules->lLength; i2++) {
                            long rhs1 = ((_List**)rules.lData)[myRules->lData[i2]]->lData[1],
                                 rhs2 = ((_List**)rules.lData)[myRules->lData[i2]]->lData[2];

                            for (long i3 = 0; i3 < countT; i3++) {
                                long anIndex = indexNT_T (i,i3),
                                     idx1    = indexNT_T (rhs1,i3),
                                     idx2    = indexNT_T (rhs2,i3),
                                     rhs1v   = precursorArray.lData[idx1],
                                     rhs2v   = followArray.lData   [idx2],
                                     lhs1v   = precursorArray.lData[anIndex],
                                     lhs2v   = followArray.lData   [anIndex];

                                if (lhs1v && rhs1v == 0) {
                                    continueLoops = true;
                                    precursorArray.lData[idx1] = 1;
                                }
                                if (lhs2v && rhs2v == 0) {
                                    continueLoops = true;
                                    followArray.lData[idx2] = 1;
                                }
                            }
                        }
                    }
                }
                /*
                for (long i = 0; i < countNT; i++)
                {
                    for (long i2 = 0; i2 < countT; i2++)
                    {
                        char buf [255];
                        snprintf (buf, sizeof(buf), "%d=>%d %s %s %s %s\n", i, i2, firstArray.lData[indexNT_T (i,i2)]?"Yes":"No ",
                                                                lastArray.lData[indexNT_T (i,i2)]?"Yes":"No ",
                                                                precursorArray.lData[indexNT_T (i,i2)]?"Yes":"No ",
                                                                followArray.lData[indexNT_T (i,i2)]?"Yes":"No ");
                        BufferToConsole (buf);
                    }
                }
                */
            }
        }
    }

    for (long clearFormulas = 0; clearFormulas < parsedFormulas.lLength; clearFormulas++)
        // clean up memory from parsed probability formulas
    {
        delete ((_Formula**)parsedFormulas.lData)[clearFormulas];
    }

    if (errorMessage.sLength) {
        WarnError      (errorMessage);
        ClearParseTree ();
    } else {
        ScanAllVariables   ();
        // RandomSampleVerify (100);
        /* temporarily removed this for node censoring (degenerate grammar) -- AFYP, August 30, 2006 */
    }

}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void        Scfg::ClearParseTree    (void)
{
    if (parseTree) {
        for (long pt = 0; pt < 256; pt++) {
            node<long>* aNode = parseTree [pt];
            if (aNode) {
                aNode->delete_tree ();
                delete (aNode);
            }
        }
        delete [] parseTree;
        parseTree = nil;
    }
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

long        Scfg::indexNT_T (long nt, long t)
{
    return nt*terminals.lLength+t;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void    Scfg::ProcessAFormula  (_FString* expression, _List & ruleProbabilities, _SimpleList& parsedFormulas, _String& errorMessage)
{
    _Formula *aFormula ;
    if (expression) { // probabilistic rule
        checkPointer (aFormula = new _Formula);

        _String  anExpression = *expression->theString;

        _Formula lhs;
        _FormulaParsingContext fpc;

        if (Parse   (aFormula, anExpression, fpc, &lhs) != HY_FORMULA_EXPRESSION) { // not a valid expression
            errorMessage = _String ("Invalid probability expression: ") & expression->theString;
        } else {
            ruleProbabilities << expression->theString;
        }
    } else { // determininstic rule (prob = 1.0)
        checkPointer (aFormula = new _Formula (new _Constant (1.0), false)); // constant 1.0
        ruleProbabilities && & _HYSCFG_NT_KEY_1;
    }

    if (errorMessage.sLength == 0) {
        parsedFormulas << (long)aFormula;
    }
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

bool    Scfg::CheckANT  (long lhs, long rhs1, long rhs2, _AVLListX& tempNT, long startSymbol)
{
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

void    Scfg::ScanAllVariables  (void)
{
    GetIndependentVars().Clear();
    GetDependentVars().Clear();
    GetCategoryVars().Clear();

    _SimpleList allVariables;
    _AVLList    scannerList(&allVariables);

    for (long formCount = 0; formCount < probabilities.GetHDim(); formCount++) {
        probabilities.GetFormula (formCount,0)->ScanFForVariables(scannerList,true,false,true,true);
    }

    scannerList.ReorderList(); // sort all scanned variables

    for (long varID = 0; varID < allVariables.lLength; varID++) {
        _Variable * aParameter = LocateVar (allVariables.lData[varID]);
        if (aParameter->IsCategory()) {
            GetCategoryVars() << allVariables.lData[varID];
        } else if (aParameter->IsIndependent()) {
            GetIndependentVars() << allVariables.lData[varID];
        } else {
            GetDependentVars() << allVariables.lData[varID];
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

_String*    Scfg::VerifyValues  (void)
{
    _Matrix * probValues = (_Matrix*)probabilities.Compute(); // initialize all the probability values
    for (long k=0; k<rules.lLength; k++) { // check that all probabilities are in [0,1]
        _Parameter  aValue = (*probValues)(k,0);
        /*
        _SimpleList *   r = (_SimpleList*)rules(k);
        char buf [256];
        if (r->lLength==2)
        {
            snprintf (buf, sizeof(buf), "rule %d [%d->%d] Pr %f\n", k,r->lData[0],r->lData[1],aValue);
        } else {
            snprintf (buf, sizeof(buf), "rule %d [%d->%d,%d] Pr %f\n", k,r->lData[0],r->lData[1],r->lData[2],aValue);
        }
        BufferToConsole(buf);
        */
        if (aValue < 0.0 || aValue > 1.0) {
            return (_String*)(_String ("Probability value for rule ") & _String (GetRuleString (k)) & " is not within [0,1]: " & aValue).makeDynamic();
        }
    }

    // now check that for each non-terminal, the sum of all probabilities is 1
    {
        for (long k=0; k<byNT2.lLength; k++) {
            _Parameter    p_sum = 0.0;
            _SimpleList * l2 = (_SimpleList*)byNT2(k),
                          * l3 = (_SimpleList*)byNT3(k);

            for (long r = 0; r < l2->lLength ; r++) {
                p_sum += (*probValues)(l2->lData[r],0);
            }

            {
                for (long r = 0; r < l3->lLength ; r++) {
                    p_sum += (*probValues)(l3->lData[r],0);
                }
            }

            if (!CheckEqual (p_sum, 1.0)) { // check within reasonable system precision
                return (_String*)(_String ("Probability values for non-terminal ") & (k+1) & " do not appear to add up to one: " & p_sum).makeDynamic();
            }
        }
    }

    return nil;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void    Scfg::RandomSampleVerify  (long samples)
{
    if (samples>0) {
        _String *   errMsg     = nil;
        long        paramCount = GetIndependentVars().lLength;
        _Parameter  stepFactor = 1./samples;

        if (paramCount > 0) { // some adjustable parameters
            _Matrix parameterBounds (paramCount,3,false,true); // used to store the lower parameter bound
            // the sampling step
            // and the current value
            for (long var = 0; var < paramCount; var++) {
                _Variable * aVar = LocateVar (GetIndependentVars()(var));
                parameterBounds.Store (var,0,aVar->GetLowerBound());
                parameterBounds.Store (var,1,(aVar->GetUpperBound()-parameterBounds(var,0))*stepFactor);
                parameterBounds.Store (var,2,aVar->Compute()->Value());
            }

            _SimpleList zeroThruNm1 (samples-1,0,1);

            for (long it = 0; it < samples; it = it+1) {
                zeroThruNm1.Permute (1);
                for (long var = 0; var < paramCount; var++) {
                    SetIthIndependent (var, parameterBounds(var,0) + parameterBounds(var,1) * zeroThruNm1.lData[var]);
                }

                if ((errMsg=VerifyValues ())) {
                    char buf [256];
                    snprintf (buf, sizeof(buf), "Breaking from RandomSampleVerify() on iteration %ld of %ld", it, samples);
                    BufferToConsole(buf);

                    break;
                }
            }
            {
                for (long var = 0; var < paramCount; var++) { // restore initial parameter values
                    SetIthIndependent (var,parameterBounds(var,2));
                }
            }
        } else { // check fixed values
            errMsg = VerifyValues ();
        }

        if (errMsg) {
            WarnError (_String(errMsg));
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void    Scfg::SetStringCorpus  (_String* varID)
{


    _PMathObj    theMatrix  = FetchObjectFromVariableByType (varID, MATRIX),
                 theString  = nil;

    if (!theMatrix) {   // if the variable is not a matrix, then try treating it as a string
        theString = FetchObjectFromVariableByType (varID, STRING);
    }

    if (theMatrix && ((_Matrix*)theMatrix)->IsAStringMatrix ()) {
        SetStringCorpus ((_Matrix*)theMatrix);
        return;
    } else if (theString) {
        _List t;
        t << ((_FString*)theString)->theString;
#ifdef _NEVER_DEFINED_
        char buf [255];     // DEBUG
        snprintf (buf, sizeof(buf), "\nSetStringCorpus() string = %s\n", (const char *) *((_FString*)theString)->theString);
        BufferToConsole (buf);
#endif
        _Matrix wrapper (t);
        SetStringCorpus (&wrapper);
        return;
    }

    _String    errMsg = *varID & " must refer either to a matrix of strings or to a single string when setting the corpus for a SCFG.";
    WarnError (errMsg);
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void    Scfg::SetStringCorpus  (_Matrix* stringMatrix)
{
    corpusChar.Clear();
    corpusInt.Clear();
    DumpComputeStructures ();

    for (long stringRow = 0; stringRow < stringMatrix->GetHDim(); stringRow++) {
        for (long stringColumn = 0; stringColumn < stringMatrix->GetVDim(); stringColumn++) {
            _FString    * aString   = (_FString *)stringMatrix->GetFormula (stringRow,stringColumn)->Compute();
            _SimpleList * tokenized = new _SimpleList;
            checkPointer (tokenized);
            _String      * errMsg   = TokenizeString (*aString->theString, *tokenized);
            if (errMsg) {
                WarnError (_String(errMsg));
                return;
            }
            corpusChar << aString->theString;
            corpusInt  << tokenized;

            /*for (long nt=0; nt < byNT2.lLength; nt++)
                for (long start = 0; start < aString->theString->sLength; start++)
                    for (long end = start; end < aString->theString->sLength; end++)
                    {
                        char buf [255];
                        snprintf (buf, sizeof(buf), "%2d %2d %2d => %4d\n", start, end, nt, scfgIndexIntoAnArray (start,end,nt,aString->theString->sLength));
                        BufferToConsole (buf);
                    }*/
            DeleteObject (tokenized);
        }

    }
    /*
    char buf [255];
    for (long c = 0; c < corpusChar.lLength; c++)
    {
        snprintf (buf, sizeof(buf), "string %d in corpusChar = %s\n", c, (const char *) *((_String *) corpusChar.lData[c]));
        BufferToConsole (buf);
    }
    */
    InitComputeStructures();
}

/*--------------------------------------------------------------------------------------------------------------------------------*/
_String*    Scfg::TokenizeString    (_String& inString, _SimpleList& outTokens)
{
    if (inString.sLength == 0) {
        return new _String ("Empty strings are not allowed as SCFG input.");
    }

    // check to see if the string is too long

    if ((inString.sLength+1.0)*inString.sLength*0.5*byNT3.lLength > _HYSCFG_UPPER_BOUND_) {
        return new _String ("The input string is too long.");
    }


    node<long>* currentTreeNode = nil; // used to keep track of where we are in the tree
    long        stringIndex     = 0;

    for (; stringIndex < inString.sLength; stringIndex++) {
        char currentChar  = inString.getChar (stringIndex);
        if (currentTreeNode == nil) { // root of the tree
            if   (!(currentTreeNode = parseTree[currentChar])) {
                break;
            }
        } else {
            long child = 1,
                 upTo  = currentTreeNode->get_num_nodes();

            for (; child <= upTo; child++)
                if ((currentTreeNode->go_down(child)->get_data() & _HYSCFG_CHARACTER_MASK_) == currentChar) {
                    currentTreeNode = currentTreeNode->go_down(child);
                    break;
                }
            if (child == upTo) {
                break;
            }
        }

        if (currentTreeNode->get_num_nodes() == 0) { // at the leaf; emit a token
            outTokens << ((currentTreeNode->get_data()&_HYSCFG_TERMINAL_MASK_) >> 8);
            currentTreeNode = nil; // reset to the root
        }
    }

    // check error conditions
    if (currentTreeNode) { // stuck in the middle of the tree
        return new _String ("Premature string end: incomplete terminal");
    }

    if (stringIndex < inString.sLength)
        return new _String (_String("Invalid terminal symbol in the input string between '") & inString.Cut (stringIndex-10, stringIndex-1) & "' and '"
                            & inString.Cut (stringIndex, stringIndex+10) & "'.");

    return nil;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

Scfg::~Scfg  (void)
{
    ClearParseTree ();
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

_String *   Scfg::GetRuleString (long ruleIdx)
{
    if (ruleIdx < 0 || ruleIdx >= rules.lLength) {
        return new _String;
    }

    _String * result = new _String (64L, true);
    _SimpleList * aRule = (_SimpleList*) rules (ruleIdx);
    _String     * temp  = (_String*)probabilities.GetFormula (ruleIdx,0)->toStr();

    (*result) << "{";
    (*result) << _String(aRule->lData[0]);
    (*result) << "}->";
    if (aRule->lLength == 2) { // NT->T
        (*result) << "\"";
        (*result) << *(_String*)terminals(aRule->lData[1]);
        (*result) << "\" : ";
    } else {
        (*result) << "{";
        (*result) << _String(aRule->lData[1]);
        (*result) << "}{";
        (*result) << _String(aRule->lData[2]);
        (*result) << "} : ";
    }

    (*result) << temp;
    DeleteObject (temp);
    result->Finalize ();
    return result;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

BaseRef     Scfg::toStr  (void)
{
    _String * result = new _String (128L, true); // allocate a buffer with the initial length of 128 characters

    for (long i = 0; i < rules.lLength; i++) {
        (*result) << new _String (GetRuleString(i));
        (*result) << "\n";
    }

    result->Finalize ();
    return result;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void        Scfg::DumpComputeStructures (void)
{
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

_Parameter      Scfg::Compute (void)
{
#ifdef      __MACPROFILE__
    ProfilerInit(collectDetailed,bestTimeBase,1000,500);
#endif

    bool first = computeFlagsI.lLength;

    _Parameter  res = 0.0,
                ip,
                useJP;

    checkParameter (useJeffreysPrior, useJP, 0.0);

    // update the probability matrix
    probabilities.Compute ();


    for (long stringID = 0; stringID < corpusChar.lLength; stringID++) {
        // assuming that production probabilities have changed since last Compute(),
        // retain 0 and 1's but clear buffer of stored 0<Pr<1 values.
        _Matrix * cachedProbs = (_Matrix*)(storedInsideP(stringID));
        for (long cid = 0; cid < cachedProbs->GetHDim(); cid++) { // reset 0,1 values to -1
            cachedProbs->Store (cid,0,-1.);
        }

        /*
        snprintf (buf, sizeof(buf), "\nComputing inside prob for string: %s\n", (const char *) *(_String *)corpusChar(stringID));
        BufferToConsole (buf);
        snprintf (buf, sizeof(buf), "\tstart symbol = %d\n\tfirst = %d\n", startSymbol, first);
        BufferToConsole (buf);
        */

        _Parameter  temp = ComputeInsideProb (0,((_String*)corpusChar(stringID))->sLength-1,stringID,startSymbol, first);
        if (temp == 0.) {
            ReportWarning (_String("Underflow detected for string ") & stringID & ". Spiking optimizer to avoid this region of parameter space.");
            return (-A_LARGE_NUMBER);
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

    /*char buf [256];
    snprintf (buf, sizeof(buf), "Compute() total inside calls = %d\n", insideCalls);
    BufferToConsole(buf);*/

    insideCalls = 0;    // reset counter

#ifdef      __MACPROFILE__
    ProfilerDump("\pSCFG Profile");
    ProfilerTerm();
#endif


#ifdef __NEVER_DEFINED__
    if (useJP >= 1.0) { // using Jeffrey's prior (the square root of the determinant of Fisher's information matrix)
        setParameter(useJeffreysPrior, 0.0);        // Compute() is called by CovarianceMatrix(), prevent infinite loop


        _Matrix *   fisherinfo  = (_Matrix *) _LikelihoodFunction::CovarianceMatrix (nil, FALSE);
        // returns the Hessian matrix
        // i.e. square matrix of second partial derivatives of log-likelihood, obtained from Compute().
        // This is the empirical Fisher information matrix.
        // Matrix must be non-negative definite, i.e. all eigenvalues are positive.
        // The determinant of the matrix is the product of eigenvalues.


        // need to remove rows and columns from Fisher's information matrix
        _Matrix     fi (fisherinfo->GetHDim()-1, fisherinfo->GetVDim()-1, false, true);

        for (long i = 0; i < fisherinfo->GetHDim() - 1; i++) {  // display matrix
            for (long j = 0; j < fisherinfo->GetVDim() - 1; j++) {
                fi.Store(i, j, (*fisherinfo) (i,j));
            }
        }


        setParameter (useJeffreysPrior, 1.0);       // reset flag

        _AssociativeList *  eigen       = (_AssociativeList *) fi.Eigensystem();



        _Matrix *           eigenvalues = (_Matrix *)eigen->GetByKey(0, MATRIX);
        _Parameter          logdet      = 0.0,      // avoid overflow from large determinants
                            checkzero   = 0.0;


        for (long i = 0; i < eigenvalues->GetHDim(); i++) { // display eigenvalues
            for (long j = 0; j < eigenvalues->GetVDim(); j++) {
                _Parameter      temp = (*eigenvalues) (i,j);

                checkzero += temp;

                if (temp < 0.0) {       // can't allow any negative eigenvalues!
                    snprintf (buf, sizeof(buf), "WARNING: negative eigenvalue.\n");
                    BufferToConsole (buf);

                    return res - 100.0;     // return a large negative number (worse log-likelihood)
                }
                logdet += log(temp);                // product of eigenvalues is the determinant

            }
        }

        /*
        //if (checkzero == 0.0) {       // all eigenvalues are zero, QR probably barfed
            for (long i = 0; i < fi.GetHDim(); i++)
            {
                for (long j = 0; j < fi.GetVDim(); j++)
                {
                    snprintf (buf, sizeof(buf), "%lf\t", fi (i, j) );
                    BufferToConsole (buf);
                }
                snprintf (buf, sizeof(buf), "\n");
                BufferToConsole (buf);
            }
            snprintf (buf, sizeof(buf), "\n");
            BufferToConsole (buf);

        //  return res - 100.0;
        //}


        snprintf (buf, sizeof(buf), "log determinant = %lf\n", logdet);
        BufferToConsole (buf);
        */

        _Parameter  jp  = 0.5 * logdet;     // this is Jeffrey's (1946) invariant prior applied to log-likelihood

        return res + jp;    // else


        // calculate the determinant
        _Matrix *   ludFisher   = (_Matrix *) fisherinfo->LUDecompose();        // returns n x (n+1) matrix, where (n+1)th column contains vector
        // of row interchanges, and the remaining nxn is the LUD.

        for (long i = 0; i < ludFisher->GetHDim(); i++) {   // display matrix
            for (long j = 0; j < ludFisher->GetVDim(); j++) {
                _Parameter  temp = (*ludFisher) (i,j);
                snprintf (buf, sizeof(buf), "%5.3lf\t", temp);
                BufferToConsole (buf);
            }
            snprintf (buf, sizeof(buf), "\n");
            BufferToConsole (buf);
        }

        snprintf (buf, sizeof(buf), "\n");
        BufferToConsole (buf);

        _Parameter  trace       = 1.0;

        for (long diag = 0; diag < ludFisher->GetHDim(); diag++) {          // diagonal of LU contains the eigenvalues
            _Parameter  temp = (*ludFisher) (diag,diag);
            snprintf (buf, sizeof(buf), "%5.3lf\t", temp);
            BufferToConsole (buf);

            trace *= temp;      // product of eigenvalues is the determinant
        }

        snprintf (buf, sizeof(buf), "<- trace\n");
        BufferToConsole (buf);

        /*
        snprintf (buf, sizeof(buf), "penalize %f\n", 0.5 * log(fabs(trace)) );
        BufferToConsole (buf);
        */

        fisherinfo->Clear();
        ludFisher->Clear();

        return res + 0.5 * log(trace);      // penalized log-likelihood

    } else {
        return res;
    }
#endif

    return res;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void        Scfg::InitComputeStructures (void)
{
    long maxStringLength = 0;
    for (long stringCount = 0; stringCount < corpusChar.lLength; stringCount++) {
        _SimpleList             emptyList;
        _GrowingVector  *aMatrix;
        _AVLListX               *searchTree;

        long                    maxDimension = ((_String*)corpusChar(stringCount))->sLength;
        maxStringLength         = MAX (maxStringLength, maxDimension);

        // unused variable? -AFYP 2006-07-07
        maxDimension = (maxDimension*(maxDimension+1)/2*byNT2.lLength/32+1)*32;

        insideProbsT  && & emptyList;
        outsideProbsT && & emptyList;

        checkPointer (searchTree = new _AVLListX ((_SimpleList*)insideProbsT(stringCount)));
        insideProbs << searchTree;
        DeleteObject (searchTree);

        checkPointer (searchTree = new _AVLListX ((_SimpleList*)outsideProbsT(stringCount)));
        outsideProbs << searchTree;
        DeleteObject (searchTree);

        checkPointer (aMatrix = new _GrowingVector);
        storedInsideP << aMatrix;
        DeleteObject (aMatrix);

        checkPointer (aMatrix = new _GrowingVector);
        storedOutsideP << aMatrix;
        DeleteObject (aMatrix);

    }
    maxStringLength = (maxStringLength * (maxStringLength+1) * byNT2.lLength / 64)+1;
    computeFlagsI.Populate (maxStringLength,0,0);
    computeFlagsO.Populate (maxStringLength,0,0);
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

inline  long    scfgIndexIntoAnArray            (long start,long end,long nt,long stringLength)
{
    return (2*stringLength-start-1)*start/2 + end + nt*(stringLength+1)*stringLength/2;

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

        (e-s) + sum_{k=0}{s-1} (L-k) = e-s + s/2 (L + L - s + 1) = e - s + (2L - s + 1) s /2 = e + (2L - s - 1)

     */



    //return (2*stringLength-start+1)*start/2 + end + nt*(stringLength+1)*stringLength/2;       // AFYP June 21, 2006
}


/*--------------------------------------------------------------------------------------------------------------------------------*/

void    setIndexBit             (long start,long end,long nt,long stringLength,_SimpleList& theArray)
{
    long theIndex = (2*stringLength-start-1)*start/2 + end + nt*(stringLength+1)*stringLength/2;
    theArray.lData[theIndex/32] |= bitMaskArray.masks[theIndex%32];
    /*char str255 [255];
    snprintf (str255, sizeof(str255),"Store %d %d %d %d into %x\n", start, end, nt, stringLength, theIndex, theArray.lData[theIndex/32]);
    BufferToConsole (str255);*/
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

bool    getIndexBit             (long start,long end,long nt,long stringLength,_SimpleList& theArray)
{
    long theIndex = (2*stringLength-start-1)*start/2 + end + nt*(stringLength+1)*stringLength/2;
    /*char str255 [255];
    snprintf (str255, sizeof(str255),"Fetch %d %d %d %d from %x to give %d\n", start, end, nt, stringLength, theArray.lData[theIndex/32], (theArray.lData[theIndex/32] & (bitMaskArray.masks[theIndex%32])) > 0);
    BufferToConsole (str255);*/
    return (theArray.lData[theIndex/32] & bitMaskArray.masks[theIndex%32]) > 0;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/

/* ---- SCFG:  Inside Probability -------------------------- */

_Parameter   Scfg::ComputeInsideProb(long from, long to, long stringIndex, long ntIndex, bool firstPass)
{

    // static long insideCalls = 0;
    insideCalls ++;


    /* quickly handle extreme cases */
    if (to > from) { // more than a single terminal in the substring
        if (((_SimpleList**)byNT3.lData)[ntIndex]->lLength == 0)
            // ntIndex non-terminal can only generate terminals
        {
            //if (firstPass) // mark this triplet as computed
            //setIndexBit (from, to, ntIndex, ((_String*)corpusChar.lData[stringIndex])->sLength, computeFlagsI);
            return 0.;
        }
    } else {    // to == from
        if (((_SimpleList**)byNT2.lData)[ntIndex]->lLength == 0)
            // ntIndex non-terminal can only generate other non-terminals
        {
            //if (firstPass) // mark this triplet as computed
            //setIndexBit (from, to, ntIndex, ((_String*)corpusChar.lData[stringIndex])->sLength, computeFlagsI);

            return 0.;
        }
    }

    /* look up the triple in the cache and decide if it's already been done */

    _AVLListX *             theAVL = (_AVLListX*)insideProbs(stringIndex);
    // see if this is already done
    long    stringL         = ((_String**)corpusChar.lData)[stringIndex]->sLength,
            tripletIndex   = scfgIndexIntoAnArray(from,to,ntIndex,stringL),
            avlIndex        = theAVL->FindLong (tripletIndex),
            matrixIndex        = -1;

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
            _Parameter currentValue = ((_GrowingVector**)storedInsideP.lData)[stringIndex]->theData[matrixIndex];
            if (currentValue >= 0.0) {
                return currentValue;
            }
            // else was reset to -1, re-compute
        }
    }

    /* have to compute stuff */

    _Parameter      insideProbValue = 0.0;

    if (to == from) { // single terminal substring; direct lookup
        long ruleIndex = ntToTerminalMap.lData[indexNT_T(ntIndex,((_SimpleList**)corpusInt.lData)[stringIndex]->lData[to])];
        if (ruleIndex >= 0) {
            insideProbValue = LookUpRuleProbability (ruleIndex);
        }
    } else {
        if (firstPass)
            // try the heuristics
        {
            _SimpleList* stringList = ((_SimpleList**)corpusInt.lData)[stringIndex];
            // if any of these conditions are TRUE, then both inside and outside Pr(s,t,i) = 0
            if (! (firstArray.lData[indexNT_T(ntIndex, stringList->lData[from])] &&
                    lastArray.lData[indexNT_T(ntIndex, stringList->lData[to])]
                    && (from == 0 || precursorArray.lData[indexNT_T(ntIndex, stringList->lData[from-1])])
                    && (to == stringL-1 || followArray.lData[indexNT_T(ntIndex, stringList->lData[to+1])]))) {
                setIndexBit (from, to, ntIndex, stringL, computeFlagsI);
                return 0.0;
            }
        }

        _SimpleList     * myNTNTRules = ((_SimpleList**)byNT3.lData)[ntIndex];

        for (long ruleIdx = 0; ruleIdx < myNTNTRules->lLength; ruleIdx++) { // loop over all NT-> NT NT rules
            long          currentRuleIndex = myNTNTRules->lData[ruleIdx];
            _Parameter    ruleProb = LookUpRuleProbability(currentRuleIndex);

            if (ruleProb > 0.0) {
                _SimpleList * currentRule = ((_SimpleList **)rules.lData)[currentRuleIndex];
                long          nt1         = currentRule->lData[1],
                              nt2        = currentRule->lData[2],
                              halfway      = from+(to-from)/2+1;

                {
                    for (long bp = from+1; bp <= halfway; bp++) { // now loop over all breakpoints
                        _Parameter t = ComputeInsideProb (from,bp-1,stringIndex,nt1,firstPass);
                        if (t>0.0) {
                            insideProbValue += t*
                                               ComputeInsideProb (bp,to,stringIndex,nt2,firstPass)*
                                               ruleProb;
                            /*
                            insideProbValue += exp( log(t) +
                                                    log(ComputeInsideProb (bp,to,stringIndex,nt2,firstPass)) +
                                                    log(ruleProb) );
                             */
                            // attempts to prevent likelihood underflow -- AFYP Aug 4, 2006
                        }
                    }
                }


                for (long bp = halfway+1; bp <= to; bp++) { // now loop over all breakpoints
                    _Parameter t = ComputeInsideProb (bp,to,stringIndex,nt2,firstPass);
                    if (t > 0.0) {

                        insideProbValue += t *
                                           ComputeInsideProb (from,bp-1,stringIndex,nt1,firstPass)*
                                           ruleProb;
                        /*
                        insideProbValue += exp( log(t) +
                                                log(ComputeInsideProb (from,bp-1,stringIndex,nt1,firstPass)) +
                                                log(ruleProb));
                         */
                    }
                }
            }
        }
    }

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

        if (avlIndex < 0) {
            long    mxID = -1;
            if (insideProbValue < 1.0) {
                mxID = ((_GrowingVector*)(storedInsideP(stringIndex)))->Store (insideProbValue);
            }
            theAVL->Insert ((BaseRef)tripletIndex, mxID);
        } else { // update values
            ((_GrowingVector*)(storedInsideP(stringIndex)))->_Matrix::Store (matrixIndex,0,insideProbValue);
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
void        Scfg::AddSCFGInfo (_AssociativeList* theList)
{
    _SimpleList indexer (corpusChar.lLength, 0, 1);
    InsertStringListIntoAVL (theList, scfgCorpus, indexer, corpusChar);
    _List       ruleStrings;
    for (long i=0; i<rules.lLength; i++) {
        ruleStrings.AppendNewInstance(GetRuleString (i));
    }

    indexer.Populate (rules.lLength, 0, 1);
    InsertStringListIntoAVL (theList, _addSCFGInfoProductions, indexer, ruleStrings);
    indexer.Populate (terminals.lLength, 0, 1);
    InsertStringListIntoAVL (theList, _addSCFGInfoTerminals, indexer, terminals);

    _Matrix * stats = new _Matrix (corpusChar.lLength,3,false,true);
    checkPointer (stats);

    for (long k=0; k<corpusChar.lLength; k++) {
        long       strL    = ((_String*)corpusChar(k))->sLength,
                   pNot0   = ((_AVLListX*)insideProbs(k))->dataList->lLength,
                   p01     = ((_GrowingVector*)storedInsideP(k))->GetUsed();

        _Parameter totalPR = strL*(strL+1.)*0.5*byNT2.lLength;
        stats->Store (k,0,totalPR); // total number of tuples (i,j,v)
        stats->Store (k,1,pNot0-p01);   // number of tuples with probability = 1
        stats->Store (k,2,p01);     // number of tuples with probability in interval (0,1)
    }
    theList->MStore (_addSCFGInfoStats, stats, false);
    stats = (_Matrix*)probabilities.Compute();
    theList->MStore (_addSCFGInfoProbabilities, stats, true);
}


/*--------------------------------------------------------------------------------------------------------------------------------*/

_Parameter   Scfg::ComputeOutsideProb(long from, long to, long stringIndex, long ntIndex, bool firstOutside, bool firstInside)
{

    /* _String nyi ("This function is under development - AFYP");
    WarnError (nyi); */

    // static long      outsideCalls = 0;

    /* get length of this string */
    long            stringL = ((_String**)corpusChar.lData)[stringIndex]->sLength;

    outsideCalls++;

    /*
    if (firstOutside)
    {
        char buf [256];
        snprintf (buf, sizeof(buf), "outside %d\t%d\t%d\n", from,to,ntIndex);
        BufferToConsole(buf);
    }
    */

    // quickly handle extreme cases
    if (from == 0 && to == stringL-1) { // dealing with entire string
        if (ntIndex == 0) { // for start symbol
            //if(firstPass) // mark this triplet as computed
            //setIndexBit(from, to, ntIndex, ((_String*)corpusChar.lData[stringIndex])->sLength, computeFlagsO);
            return 1.0;
        } else {            // not the start symbol
            //if(firstPass) // mark this triplet as computed
            //setIndexBit(from, to, ntIndex, ((_String*)corpusChar.lData[stringIndex])->sLength, computeFlagsO);
            return 0.;
        }
    }

    if (to > from) {
        if (( (_SimpleList**)byNT3.lData)[ntIndex]->lLength == 0) {
            return 0.;  // NT can only emit terminal; inside is 0, so outside is useless
        }
    } else {    // from==to
        if (( (_SimpleList**)byNT2.lData)[ntIndex]->lLength == 0) {
            return 0.;  // NT can only generate NT NT
        }
    }


    if (!firstInside && firstOutside) { // re-use information from inside probs.
        if (ComputeInsideProb(from,to,stringIndex,ntIndex,firstInside) == 0) {
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
    _AVLListX *     theAVL          = (_AVLListX*)outsideProbs(stringIndex);

    long            tripletIndex    = scfgIndexIntoAnArray(from,to,ntIndex,stringL),    // convert (s,t,i)
                    avlIndex     = theAVL->FindLong (tripletIndex),              // retrieve key
                    matrixIndex        = -1;

    if (avlIndex < 0) { // the triplet is not in the search tree
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
            _Parameter currentValue = ((_GrowingVector**)storedOutsideP.lData)[stringIndex]->theData[matrixIndex];
            if (currentValue >= 0.0) {
                return currentValue;
            }
        }
    }

    /* have to compute stuff; note, nothing below this line should get called after the first pass */
    if (firstOutside)
        // try the heuristics
    {
        _SimpleList* stringList = ((_SimpleList**)corpusInt.lData)[stringIndex];
        // if any of these conditions are TRUE, then both inside and outside Pr(s,t,i) = 0
        if (! (firstArray.lData[indexNT_T(ntIndex, stringList->lData[from])] &&
                lastArray.lData[indexNT_T(ntIndex, stringList->lData[to])]
                && (from == 0 || precursorArray.lData[indexNT_T(ntIndex, stringList->lData[from-1])])
                && (to == stringL-1 || followArray.lData[indexNT_T(ntIndex, stringList->lData[to+1])]))) {
            setIndexBit (from, to, ntIndex, stringL, computeFlagsO);
            return 0.0;
        }
    }

    _Parameter      outsideProbValue = 0.0;

    /* loop over all j->ki productions */
    _SimpleList *   myNTNTRules = ((_SimpleList**)byRightNT2.lData)[ntIndex];
    for (long ruleIdx = 0; ruleIdx < myNTNTRules->lLength; ruleIdx++) {
        long        currentRuleIndex = myNTNTRules->lData[ruleIdx];
        _Parameter  ruleProb = LookUpRuleProbability(currentRuleIndex);     // Pr(j->ki)

        if (ruleProb > 0.0) {
            _SimpleList *   currentRule = ((_SimpleList **)rules.lData)[currentRuleIndex];  // rule = list of three integers
            long            j           = currentRule->lData[0],
                            k           = currentRule->lData[1];

            for (long leftBisect = 0; leftBisect < from; leftBisect++) {
                _Parameter  t = ComputeInsideProb(leftBisect,from-1,stringIndex,k,firstInside);
                if (t > 0.0)
                    outsideProbValue += t *
                                        ComputeOutsideProb(leftBisect,to,stringIndex,j,firstOutside,firstInside) *
                                        ruleProb;
            }
        }
    }

    /* loop over all j->ik productions */
    myNTNTRules = ((_SimpleList**)byRightNT1.lData)[ntIndex];
    {
        for (long ruleIdx = 0; ruleIdx < myNTNTRules->lLength; ruleIdx++) {
            long        currentRuleIndex = myNTNTRules->lData[ruleIdx];
            _Parameter  ruleProb = LookUpRuleProbability(currentRuleIndex);     // Pr(j->ik)

            if (ruleProb > 0.0) {
                _SimpleList *   currentRule = ((_SimpleList **)rules.lData)[currentRuleIndex];
                long            j           = currentRule->lData[0],
                                k           = currentRule->lData[2];

                for (long rightBisect = to+1; rightBisect < stringL; rightBisect++) {
                    _Parameter  t = ComputeInsideProb(to+1,rightBisect,stringIndex,k,firstInside);
                    if (t > 0.0)
                        outsideProbValue += t *
                                            ComputeOutsideProb(from,rightBisect,stringIndex,j,firstOutside,firstInside) *
                                            ruleProb;
                }
            }
        }
    }

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
                mxID = ((_GrowingVector*)(storedOutsideP(stringIndex)))->Store (outsideProbValue);
            }
            theAVL->Insert((BaseRef)tripletIndex, mxID);    // if Pr = 1.0, stores -1
        } else {    // update values
            ((_GrowingVector*)(storedOutsideP(stringIndex)))->_Matrix::Store (matrixIndex,0,outsideProbValue);
        }
    }

    if (firstOutside) { // flag (s,t,i) as computed
        setIndexBit (from,to,ntIndex,stringL,computeFlagsO);
    }

    return outsideProbValue;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

_Matrix*     Scfg::Optimize (void)  /* created by AFYP, 2006-06-20 */
{
    _Parameter  useJP,
                sOM;

    checkParameter (useJeffreysPrior, useJP, 0.0);
    checkParameter (scfgOptimizationMethod, sOM, 0.0);

    if (useJP > 0.0 || sOM > 0.0) {
        return _LikelihoodFunction::Optimize();    // Jeffrey's prior violates EM
    }

    char    buf [256];      // debugging buffer


#ifdef      __MACPROFILE__
    ProfilerInit(collectDetailed,bestTimeBase,1000,500);
#endif


    // populated with zeros for longest string by InitComputeStructures();
    // inside flags cleared by calling Scfg::Compute()
    bool    firstInside     = computeFlagsI.lLength,
            firstOutside = computeFlagsO.lLength;


    // update probability matrix -- converts formulas to numerical values
    probabilities.Compute ();


    // create container for new production rule probability estimates
    //  numerator in column 1, denominator in column 2
    long        nRules = (unsigned long) probabilities.GetHDim();
    _Matrix     nextProbs(nRules, 2, FALSE, TRUE);



    // scan production rule formulas for linking constraints, done once for entire optimization
    _List           links;
    _SimpleList *   thisLink = new _SimpleList();

    for (long ruleCount = 0; ruleCount < nRules; ruleCount++) {
        _Formula *  thisRule = probabilities.GetFormula (ruleCount, 0);

        if ( !(thisRule->IsAConstant()) ) {
            *thisLink << ruleCount;

            for (long nextrule = ruleCount+1; nextrule < nRules; nextrule++) {
                _Formula *  thatRule = probabilities.GetFormula (nextrule, 0);
                if ( thisRule->EqualFormula (thatRule) ) {
                    *thisLink << nextrule;
                    snprintf (buf, sizeof(buf), "linked rule %ld to %ld\n", ruleCount+1, nextrule+1);
                    BufferToConsole (buf);
                }
            }

            if (thisLink->lLength > 1) {
                links && thisLink;    // append duplicate
            }
            thisLink->Clear();
        }
    }


    /* calculate current corpus log-likelihood */
    _Parameter  newLk = Compute(),
                oldLk = newLk;
    long        rep = 0;


    /* ========================= */
    /*  start optimization loop  */
    /* ========================= */

    do {
        oldLk = newLk;

        for (long ruleCount = 0; ruleCount < nRules; ruleCount++) {
            nextProbs.Store(ruleCount, 0, 0);   // zero all entries in container matrix
            nextProbs.Store(ruleCount, 1, 0);
        }



        // loop over strings in corpus
        for (long stringID = 0; stringID < corpusChar.lLength; stringID++) {
            long    stringL     = ((_String**)corpusChar.lData)[stringID]->sLength,
                    countNT      = byNT2.lLength;

            _Parameter  denom       = 0,
                        stringProb = ComputeInsideProb(0, stringL-1, stringID, startSymbol, firstInside);
            // by calculating this, firstInside should be reset to FALSE

            _SimpleList*    stringList = ((_SimpleList**)corpusInt.lData)[stringID];


            /* After computing string Pr., absence in AVL indicates inside Pr = 0   */
            /* so treat firstInside as FALSE below...                               */

            if (stringProb > 0.) {
                for (long ntIndex = 0; ntIndex < countNT; ntIndex++) {  // loop over all non-terminal symbols

                    /* calculate denominator (number of times i-th NT was used in string */
                    denom = 0;

                    for (long from = 0; from < stringL; from++) {
                        if (( (_SimpleList**)byNT3.lData)[ntIndex]->lLength == 0) { // only i->m
                            denom += ComputeInsideProb(from,from,stringID,ntIndex,FALSE) *
                                     ComputeOutsideProb(from,from,stringID,ntIndex,firstOutside,FALSE) /
                                     stringProb;
                        } else {    // NT has i->jk productions
                            //for (long to = from; to < stringL; to++)
                            for (long to = from + 1; to < stringL; to++) {  // count only i->jk

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
                                        _Parameter  ip  = ((_GrowingVector**)storedInsideP.lData)[stringID]->theData[matrixIndex],
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

                    /*
                    snprintf (buf, sizeof(buf), "string = %d\tnt = %d\tdenom = %f\n", stringID, ntIndex, denom);
                    BufferToConsole (buf);
                     */

                    if (denom > 0.) {
                        // loop over i->m productions of current NT
                        _SimpleList *   imRules = ((_SimpleList**)byNT2.lData)[ntIndex];

                        for (long rC = 0; rC < imRules->lLength; rC++) {
                            _Parameter      numer       = 0.;
                            long            ruleIndex   = imRules->lData[rC];
                            _SimpleList *   currentRule = ((_SimpleList**)rules.lData)[ruleIndex];
                            long            termIndex   = currentRule->lData[1];
                            _Parameter      ruleProb    = LookUpRuleProbability (ruleIndex);

                            if (ruleProb == 1.0 || ruleProb == 0.0) {
                                continue;
                            }

                            char *  thisString = ((_String**)corpusInt.lData)[stringID]->getStr();

                            for (long from = 0; from < stringL; from++) {
                                long    thisSymbol = thisString[from];
                                if (thisSymbol == termIndex) {
                                    _Parameter  tryProb = ComputeInsideProb(from,from,stringID,ntIndex,FALSE);
                                    if (tryProb > 0.0) {
                                        numer += tryProb *
                                                 ComputeOutsideProb(from,from,stringID,ntIndex,firstOutside,FALSE) /
                                                 stringProb;
                                    }
                                }
                            }
                            // update i->m production probability

                            /*
                             snprintf (buf, sizeof(buf), "storing %f/%f into ruleIndex %d\n", numer, denom, ruleIndex);
                             BufferToConsole(buf);
                             */

                            nextProbs.theData[ruleIndex] += numer;
                            nextProbs.theData[nRules + ruleIndex] += denom;
                        }



                        // loop over all i->jk productions of current NT
                        _SimpleList *   ijkRules = ((_SimpleList**)byNT3.lData)[ntIndex];

                        /*
                         snprintf (buf, sizeof(buf), "nt %d has %d ijk rules\n", ntIndex, ijkRules->lLength);
                         BufferToConsole(buf);
                         */
                        {
                            for (long rC = 0; rC < ijkRules->lLength; rC++) {
                                long            ruleIndex   = ijkRules->lData[rC];
                                _Parameter      ruleProb    = LookUpRuleProbability (ruleIndex),
                                                numer        = 0.;

                                _SimpleList *   currentRule = ((_SimpleList**)rules.lData)[ruleIndex];
                                long            jIndex      = currentRule->lData[1],
                                                kIndex      = currentRule->lData[2];

                                if (ruleProb == 0.0 || ruleProb == 1.0) {
                                    continue;
                                }

                                /*
                                snprintf (buf, sizeof(buf), "rule %d (%d->%d,%d)\tPr = %f\n", ruleIndex, ntIndex, jIndex, kIndex, ruleProb);
                                BufferToConsole (buf);
                                 */


                                for (long from = 0; from < stringL-1; from++) { /* calculate equation (20) */
                                    // use heuristics to skip zero productions

                                    if ( ! (firstArray.lData[indexNT_T(ntIndex,stringList->lData[from])]
                                            && (from == 0 || precursorArray.lData[indexNT_T(ntIndex,stringList->lData[from-1])])) ) {
                                        continue;
                                    }


                                    for (long to = from; to < stringL; to++) {  /* BUG? <- start at s=t or s=t+1 ??? */
                                        // more short-cuts from heuristics

                                        if (! (lastArray.lData[indexNT_T(ntIndex, stringList->lData[to])]
                                                && (to == stringL-1 || followArray.lData[indexNT_T(ntIndex, stringList->lData[to+1])])) ) {
                                            continue;
                                        }


                                        _Parameter  op  = ComputeOutsideProb(from,to,stringID,ntIndex,firstOutside,FALSE);

                                        if (op == 0) {
                                            continue;    // yet another short-cut
                                        }


                                        for (long bisect = from; bisect < to; bisect++) {
                                            _Parameter ip = ComputeInsideProb(from,bisect,stringID,jIndex,FALSE);
                                            if (ip > 0.) {
                                                numer += ruleProb * ip *
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

                                nextProbs.theData[ruleIndex] += numer;
                                nextProbs.theData[nRules + ruleIndex] += denom;
                            }
                        }

                    } // endif denom > 0


                } // end loop over NT


            } else {    // stringProb <= 0
                snprintf (buf, sizeof(buf), "WARNING: String probability <= 0 (%f)\n", stringProb);
                BufferToConsole (buf);
            }

            // repopulate with zeros for next string
            if (firstInside) {
                computeFlagsI.Populate(computeFlagsI.lLength,0,0);
            }
            if (firstOutside) {
                computeFlagsO.Populate(computeFlagsO.lLength,0,0);
            }



            // clear stored inside and outside probabilities -- need to re-estimate and re-fill AVLs
            _Matrix * cachedProbs   = (_Matrix*)(storedInsideP(stringID));
            {
                for (long cid = 0; cid < cachedProbs->GetHDim(); cid++) {
                    cachedProbs->Store(cid,0,-1.);
                }
            }
            cachedProbs = (_Matrix*)(storedOutsideP(stringID));
            {
                for (long cid = 0; cid < cachedProbs->GetHDim(); cid++) {
                    cachedProbs->Store(cid,0,-1.);
                }
            }
        } // end loop over strings



        // search through rules and sum linked probabilities

        for (long linkIndex = 0; linkIndex < links.lLength; linkIndex++) {
            _SimpleList *   thisLink    = (_SimpleList*)links.lData[linkIndex];
            _Parameter      linkNumer   = 0.,
                            linkDenom = 0.;
            /*
            snprintf (buf, sizeof(buf), "Link %d contains ", linkIndex);
            BufferToConsole (buf);
             */
            {
                for (long lcount = 0; lcount < thisLink->lLength; lcount++) {
                    /*
                    snprintf (buf, sizeof(buf), "%d(%3.3f/%3.3f) ", thisLink->lData[lcount], nextProbs.theData[thisLink->lData[lcount]],
                            nextProbs.theData[nRules+thisLink->lData[lcount]]);
                    BufferToConsole (buf);
                     */
                    linkNumer += nextProbs.theData[thisLink->lData[lcount]];
                    linkDenom += nextProbs.theData[nRules+thisLink->lData[lcount]];
                }
            }

            for (long lcount = 0; lcount < thisLink->lLength; lcount++) {
                // update linked probabilities
                nextProbs.theData[thisLink->lData[lcount]] = linkNumer / thisLink->lLength;
                nextProbs.theData[nRules+thisLink->lData[lcount]] = linkDenom / thisLink->lLength;
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
                _Variable *     theParameter    = LocateVar (itsVariables.lData[0]);
                _Constant       nextValue (nextProbs.theData[ruleCount] / nextProbs.theData[nRules+ruleCount]);


                if (nextProbs.theData[nRules+ruleCount] > 0.0) {    // skip productions with denominators = 0
                    // Temporary fix -- AFYP, Sept 27, 2006
                    /*
                    snprintf (buf, sizeof(buf), "set rule %d from %f to %f\n", ruleCount, LookUpRuleProbability(ruleCount),
                         nextProbs.theData[ruleCount] / nextProbs.theData[nRules+ruleCount]);
                    BufferToConsole(buf);
                     */
                    theParameter->SetValue (nextValue.Compute());
                }

            }
        }

        if (firstInside) {
            firstInside = FALSE;    // reset flags
        }
        if (firstOutside) {
            firstOutside = FALSE;
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

        if (newLk < oldLk) {
            snprintf (buf, sizeof(buf), "\nERROR: log L = %f is lower than previous value %f at step %ld \n", newLk, oldLk, rep);
            BufferToConsole (buf);
            // break;
        }

    } while (fabs(newLk - oldLk) > SCFG_OPTIMIZATION_THRESHOLD);    // should be non-absolute - afyp 11/30/2006

    snprintf (buf, sizeof(buf), "Used %ld iterations in expectation maximization.\n", rep);
    BufferToConsole (buf);

    delete  thisLink;   // release allocated memory

#ifdef      __MACPROFILE__
    ProfilerDump("\pSCFGProfile");
    ProfilerTerm();
#endif

    // SLKP:  TBI this actually needs to be populated still
    /* return new _Matrix (1,1,false,true); */
    return NULL;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

_String* Scfg::SpawnRandomString(long ntIndex, _SimpleList* storageString)
{
    if (ntIndex < 0) {
        checkPointer (storageString = new _SimpleList);
        SpawnRandomString (startSymbol, storageString);
        _String      *backString = new _String (storageString->lLength, true);
        checkPointer (backString);
        for (long k=0; k<storageString->lLength; k++) {
            (*backString) << (_String*)terminals (storageString->lData[k]);
        }
        backString->Finalize();
        DeleteObject (storageString);
        return backString;
    }

    _Parameter      randomValue = genrand_real2 (),
                    sum         = 0.;

    long            ruleIndex   = 0;
    _SimpleList*    aList       = (_SimpleList*)byNT2(ntIndex),
                    *  aRule;

    // loop through terminal rules (X->x) and sum probabilities
    for (; ruleIndex<aList->lLength && sum < randomValue ; ruleIndex++) {
        sum += LookUpRuleProbability (aList->lData[ruleIndex]);
    }

    if (sum >= randomValue) {
        aRule = (_SimpleList*)rules(aList->lData[ruleIndex-1]);
        (*storageString) << aRule->lData[1];
        return nil;
    }

    ruleIndex = 0;
    aList       = (_SimpleList*)byNT3(ntIndex);
    for (; ruleIndex<aList->lLength && sum < randomValue ; ruleIndex++) {
        sum += LookUpRuleProbability (aList->lData[ruleIndex]);
    }

    if (sum >= randomValue) {
        aRule = (_SimpleList*)rules(aList->lData[ruleIndex-1]);
        SpawnRandomString (aRule->lData[1],storageString);
        SpawnRandomString (aRule->lData[2],storageString);
    } else {
        _String oops ("SCFG::SpawnRandomString() randomValue ");
        oops = oops & randomValue & " exceeded sum " & sum;
        oops = oops & ": nt=" & ntIndex & " stor=" & (_String *) storageString->toStr();
        WarnError (oops);
    }

    return nil;

}



/*--------------------------------------------------------------------------------------------------------------------------------*/
/*  Implementation of CYK algorithm by AFYP     2006-07-12                                                                        */
/*--------------------------------------------------------------------------------------------------------------------------------*/
_String *   Scfg::BestParseTree(void)
{
    long            countNT         = byNT2.lLength;
    bool            firstPass       = computeFlagsI.lLength;

    _String *       parseTreeString = new _String();
    //char          buf [4096];


    for (long stringIndex = 0; stringIndex < corpusChar.lLength; stringIndex++) {
        long    stringL     = ((_String**)corpusChar.lData)[stringIndex]->sLength;

        // initialize AVL tree for storing non-zero subtrees (i,j,v) to refer to _SimpleList
        _SimpleList     triplets;
        _AVLListX *     theAVL;

        checkPointer (theAVL = new _AVLListX (&triplets));

        _SimpleList             argMaxYZK;      // stores (y,z,k) keyed by argmax(i,j,v) in AVL
        _GrowingVector  *theMatrix;     // stores likelihood for subtree (i,j,v)

        checkPointer (theMatrix = new _GrowingVector);

        {
            for (long from = 0; from < stringL; from++) {   // initialization
                for (long ntIndex = 0; ntIndex < countNT; ntIndex++) {
                    long        tripletIndex    = scfgIndexIntoAnArray (from,from,ntIndex,stringL),
                                mxID         = -1;
                    _Parameter  lk = ComputeInsideProb(from,from,stringIndex,ntIndex,firstPass);

                    if (lk > 0.) {
                        mxID = theMatrix->Store (lk);
                        theAVL->Insert ((BaseRef)tripletIndex, mxID);
                        for (long count = 0; count < 3; count++) {
                            argMaxYZK << 0;    // append (0,0,0)
                        }
                    }
                }
            }
        }

        for (long from = 0; from < stringL-1; from++) { // iterate over all substrings and non-terminals
            for (long to = from+1; to < stringL; to++) {
                for (long ntIndex = 0; ntIndex < countNT; ntIndex++) {
                    _Parameter      maxLk = 0;
                    long            maxLeft,
                                    maxRight,
                                    maxBisect;
                    _SimpleList *   itsRules = ((_SimpleList **) byNT3.lData)[ntIndex];

                    for (long ruleIdx = 0; ruleIdx < itsRules->lLength; ruleIdx++) {    // iterate over all productions
                        long            currentRuleIndex    = itsRules->lData[ruleIdx];
                        _SimpleList *   currentRule         = ((_SimpleList**)rules.lData)[currentRuleIndex];
                        _Parameter      ruleProb            = LookUpRuleProbability(currentRuleIndex);
                        long            leftNT              = currentRule->lData[1],
                                        rightNT               = currentRule->lData[2];

                        if (ruleProb > 0.) {
                            for (long bisect = from; bisect < to; bisect++) {       // iterate over all bisects of substring
                                _Parameter  tryProb = ComputeInsideProb(from,bisect,stringIndex,leftNT,firstPass);
                                if (tryProb > 0.) {
                                    _Parameter  lk  = ruleProb * tryProb *
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

                    if (maxLk > 0) {
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
        CykTraceback (0,stringL-1,0,stringIndex,theAVL,&argMaxYZK,theMatrix,parseTreeString);

        (*parseTreeString) = (*parseTreeString) & "\n"; // separate tree strings
    }


    // BufferToConsole((char *)(const char *) (*parseTreeString));
    parseTreeString -> Finalize();
    return (parseTreeString);
}


/*------------------------------------------------------------------------------------------------------------------------------*/
/*  Converts stack of best-parse productions into a tree string.- AFYP                                                          */
void    Scfg::CykTraceback (long i, long j, long v, long stringIndex, _AVLListX * theAVL, _SimpleList * theYZKs, _GrowingVector * theMatrix, _String * parseTreeString)
{
    //  for sub-string (i,j) and non-terminal (v), string of length (stringL)

    long    stringL         = ((_String**)corpusChar.lData)[stringIndex]->sLength,
            tripletIndex  = scfgIndexIntoAnArray(i,j,v,stringL),
            avlIndex       = theAVL -> Find ((BaseRef)tripletIndex);

    _String * corpusString  = ((_String**)corpusChar.lData)[stringIndex];

    if (avlIndex < 0) { // value not in tree
        ReportWarning (_String("ERROR: Unknown triplet encountered in CYK traceback: (") & i & "," & j & "," & v & ")");
    } else {
        long        matrixIndex = theAVL->GetXtra (avlIndex),
                    y           = theYZKs->lData[matrixIndex * 3],
                    z           = theYZKs->lData[matrixIndex * 3 + 1],
                    k           = theYZKs->lData[matrixIndex * 3 + 2];
#ifdef __NEVER_DEFINED__
        if (y == 0 && z == 0 && k == 0) {   // node terminates
            snprintf (buf, sizeof(buf), "%d:%d", i, v);
            // BufferToConsole(buf);
            (*parseString) << (const char *) buf;
        } else {                        // node spawns two children

            snprintf (buf, sizeof(buf), "(");
            // BufferToConsole(buf);
            (*parseString) << (const char *) buf;

            CykTraceback(i,k,y,stringIndex,theAVL,theYZKs,theMatrix,parseString);

            snprintf (buf, sizeof(buf), ",");
            // BufferToConsole(buf);
            (*parseString) << (const char *) buf;

            CykTraceback(k+1,j,z,stringIndex,theAVL,theYZKs,theMatrix,parseString);

            snprintf (nodename, sizeof(nodename), "Node%d", tripletIndex);
            snprintf (buf, sizeof(buf), ")%s:%d", nodename, v);
            // BufferToConsole(buf);
            (*parseString) << (const char *) buf;

        }
#else
        if (y==0 && z==0 && k==0) { // node terminates
            (*parseTreeString) = (*parseTreeString) & "(" & v & " " & corpusString->sData[i] & ")";
            /*
            snprintf (buf, sizeof(buf), "(%d %s)", v, corpusString->sData[i]);
            (*parseString) << (const char *) buf;
             */
        } else {
            (*parseTreeString) = (*parseTreeString) & "(" & v & " ";
            /*
            snprintf (buf, sizeof(buf), "(%d ", v);
            (*parseString) << (const char *) buf;
            */

            CykTraceback(i,k,y,stringIndex,theAVL,theYZKs,theMatrix,parseTreeString);

            //(*parseTreeString) = (*parseTreeString) & " ";

            CykTraceback(k+1,j,z,stringIndex,theAVL,theYZKs,theMatrix,parseTreeString);

            (*parseTreeString) = (*parseTreeString) & ")";

            /*
            snprintf (buf, sizeof(buf), ")");
            (*parseString) << (const char *) buf;
             */
        }
#endif
    }
}

#endif

/*--------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------*/



