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

#ifndef  _HY_INCLUDE_SCFG_

#define  _HY_INCLUDE_SCFG_
#define  _USE_HYPHY_HOOKS_
// if defined, the code below will use HyPhy backend classes instead of std namespace things
// everything inside _USE_HYPHY_HOOKS_ ifdefs is written by SLKP and needs to be checked :)

#ifndef  _USE_HYPHY_HOOKS_
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "randomc.h"
#include "strings.h"

#define     MAX_NAME_LENGTH         50
#define     TRUE                    1
#define     FALSE                   0

using namespace std;
#else
#include "list.h"
#include "simplelist.h"
#include "avllist.h"
#include "avllistx.h"
#include "classes.h"
#include "likefunc.h"
#endif

#ifdef __NEVER_DEFINED__

/* ==================== */
/*      STRUCTURES      */
/* ==================== */

struct Prob {       // inner or outer probability stored with sub-string indices

    int     i;      // non-terminal index
    int     s;      // left index in string
    int     t;      // right index in string

    double  p;      // non-zero values only
};


struct Production {     // identifies production rules
    int     i, j, k, m;
};


struct ProductionLink {     // for linking production rule probabilities
    int             size;
    Production *    pvec;
};


struct Triplet {    // for trace-backs in CYK algorithm
    int i1, i2, i3;
};



/* ================================================================ */
/*  CLASS: ParseNode                                                */
/* ---------------------------------------------------------------- */
/*  Root node allocates memory for 95 printable ASCII characters.   */
/*  Character not occurring in grammar corresponds to NULL.         */
/*  Non-root nodes, push-back vector as required.                   */
/*                                                                  */
/*  itsIndex = terminal symbol integer (0-index)                    */
/*  itsChar = terminal symbol                                       */
/* ================================================================ */

class ParseNode
{
private:
    ParseNode *             itsParent;
    vector <ParseNode *>    itsChildren;
    char                    itsChar;
    int                     itsIndex;

public:
    ParseNode() {   // for root node
        itsParent = NULL;
        itsChar = '\0';
        itsIndex = -1;
        for (int i = 0; i < 95; i++) {
            itsChildren.push_back(NULL);
        }
    }

    ParseNode(ParseNode * p, char c) {
        itsParent = p;
        itsChar = c;
        itsIndex = -1;
    }

    // default destructor
    ~ParseNode()    { }

    // member accessor functions
    void SetParent(ParseNode * p)   {
        itsParent = p;
    }
    ParseNode * GetParent()         {
        return itsParent;
    }

    ParseNode * AddChild(char c) {
        ParseNode * pn = new ParseNode(this, c);
        itsChildren.push_back(pn);
        return pn;
    }

    ParseNode * SetChild(int i, char c) {
        ParseNode * pn = new ParseNode(this, c);
        itsChildren.at(i) = pn;
        return pn;
    }

    ParseNode * GetChild(int i) {
        return itsChildren.at(i);
    }

    void SetIndex(int i)        {
        itsIndex = i;
    }
    int GetIndex()              {
        return itsIndex;
    }

    void SetChar(char c)        {
        itsChar = c;
    }
    char GetChar()              {
        return itsChar;
    }

    int HowManyChildren()       {
        return itsChildren.size();
    }

    ParseNode * FindChild(char c) {
        for (int i = 0; i < itsChildren.size(); i++)
            if (c == itsChildren.at(i)->GetChar()) {
                return itsChildren.at(i);
            }
        return NULL;
    }
};


#endif

// SLKP 20070525: Moved _GrowingVector to matrix.h
// and renamed it _GrowingVector


/* ================================================================ */
/*  CLASS: Scfg (Stochastic context-free grammar)                   */
/* ---------------------------------------------------------------- */
/*  Required functions:                                             */
/*      - generate string from grammar                              */
/*      - parse string (membership problem)                         */
/*      - calculate inside/outside probabilities of a given string  */
/* ================================================================ */

class Scfg
#ifdef _USE_HYPHY_HOOKS_
    : public _LikelihoodFunction
#endif
{

private:

#ifdef _NEVER_DEFINED_
    /* GRAMMAR SYMBOLS */

    int                         numNT;      // number of non-terminal symbols
    int                         numT;       //  "   "   terminal        "

    int                         firstPreTerm;   // index of first non-terminal symbol
    // in string of form i->m

    ParseNode *                 itsParseTree;

    /* PRODUCTION PROBABILITIES */

    /* use symbol indices for production probability matrices   */
    /* ijk initialized N x N x N, where N = numNT               */
    /* im initialized N x M, where M = numT                     */
    double *                    itsIJKProbs;    // i->jk production probabilities
    double *                    itsIMProbs;     // i->m         "       "


    /* STRING PROBABILITY STORAGE */

    // stores whether the inside and outside probabilities for
    // the production of sub-string (s,t) from the i-th non-terminal symbol
    // has been evaluated;
    // and if so whether it takes the value 0, 1, or something inbetween.

    /*      character key                                       */
    /* -------------------------------------------------------- */
    /*                              inside                      */
    /*  o               NA      0       1       0 < p < 1       */
    /*  u   ----------------------------------------------      */
    /*  t   NA          @       A       B       C               */
    /*  s   0           D       E       F       G               */
    /*  i   1           H       I       J       K               */
    /*  d   0 < p < 1   L       M       N       O               */
    /*  e   ----------------------------------------------      */



    char **                     itsStringProbs;

    // stores inside and outside probabilities between 0 and 1
    vector <Prob>               itsInsideProbs;
    vector <Prob>               itsOutsideProbs;

    vector <ProductionLink>     itsLinked;

    /* STRING AND CORPUS */

    String *                    itsString;
    vector <String *>           itsCorpus;

    int                         lengthString;

    TRandomMersenne *           pMersenne;

    ifstream                    infile;
    ofstream                    outfile;

#endif

#ifdef                  _USE_HYPHY_HOOKS_

#endif

public:
    // general specification of grammar by string
#ifndef                     _USE_HYPHY_HOOKS_
    Scfg(char * grammar, int nt, int t, TRandomMersenne * pMers);
#else
    Scfg(_AssociativeList*, _AssociativeList *, long = 0);
    /* SLKP:

        the first  associative list defines: i->m  type terminal rules
        the second one                     : i->jk non-terminal  rules
        the third one                      : the index of the start NT symbol

        the start symbol is assumed to be the non-terminal indexed '0'

        Both arrays are assumed to have integer keys (most naturally: from 0 to the number of rules - 1)
        Each of the lists has associative arrays as entries.

        (1). Terminal rules are defined as:
             Term_rule ["L"] = integer - left hand side non-terminal;
             Term_rule ["T"] = "non-terminal literal";
             Term_rule ["P"] = either empty (0) which is taken to mean that the rule is deterministic
                             = or a string with the ID of the parameter controlling the production probability

        (2). Non-terminal rules are defined as:
             NT_rule ["L"]   = integer - left hand side non-terminal;
             NT_rule ["1"]   = integer - first  right hand side non-terminal symbol, referred to by its index in the NT associative array;
             NT_rule ["2"]   = integer - second right hand side non-terminal symbol, referred to by its index in the NT associative array;
             NT_rule ["P"]   = either empty (0) which is taken to mean that the rule is deterministic
                             = or a string with the ID of the parameter controlling the production probability


        The constructor does the following:
        (a). Check each rule for validity
        (b). Check that the terminals form a proper prefix code and build a parse tree for the input; build a map of 'index' <-> 'literal'
        (c). Verify several condition on the rules involving non-terminal symbols to check for possible grammar specification issues
        (d). Build a list of production rules, and associated probabilities
        (e). Compile a list of variables involved as probabilities in the production rules


    */

    virtual         ~Scfg(void); // virtual destructor because of inheritance

    BaseRef             toStr                   (void); // a string representation of the SCFG

    virtual void        ScanAllVariables        (void);
    virtual _String*    GetRuleString           (long); // return a string representation of a derivation rule
    virtual _String*    VerifyValues            (void);
    // for a given set of parameter values, check that all probabilities are
    // (a) in [0,1]
    // (b) for a fixed NT 'i' the sum of all probabilities for the rules
    //     involving the NT is 1
    // returns nil if all is good

    virtual void        RandomSampleVerify      (long);
    // draw a number of values using LHC sampling on parameter bounds
    // and check that they all define valid sets of probabilities

    virtual void        SetStringCorpus         (_String*);
    virtual void        SetStringCorpus         (_Matrix*);
    // the strings to set SCFG string corpus
    // the first one does it from a HyPhy variable identifier, which is assumed
    // to reference either a string or a matrix of strings (it calls the second)
    // the second one does the work; lexing the input strings and converting

    virtual _Parameter          Compute                 (void);
    // compute the derivation probability of the current corpus
    virtual _Matrix*            Optimize                ();
    // train the grammar using current corpus

    virtual void                AddSCFGInfo             (_AssociativeList*);
    // store various diagnostics about the grammar into the list

    virtual _String*            SpawnRandomString       (long = -1, _SimpleList* = nil);
    // generate a random string using current production probabiltities
    // the argument is the index of a non-terminal (or -1 for the start)
    // from which to generate a substring (called recursively)
    // the string argument is the storage, created by the first call

    virtual _String *           BestParseTree           (void);

    virtual void                CykTraceback            (long,long,long,long,_AVLListX *,_SimpleList *,_GrowingVector *,_String *);

#endif

#ifdef _NEVER_DEFINED_

    // destructor
    ~Scfg()     { }


    // member access functions
    int     GetNumberNTs()                      {
        return numNT;
    }
    int     GetNumberTerms()                    {
        return numT;
    }

    void    SetFirstPT(int fpt)     {
        firstPreTerm = fpt;
    }
    int     GetFirstPreTerm()       {
        return firstPreTerm;
    }


    // allocate memory for production probabilities
    void    InitIJKProbs();
    void    InitIMProbs();


    /* -------- PRODUCTION PROBABILITY FUNCTIONS ------------- */
    void    SetProb(int i, int m, double prob) {
        itsIMProbs[i*numT + m] = prob;
    }

    void    SetProb(int i, int j, int k, double prob) {
        itsIJKProbs[i * numNT * numNT + j * numNT + k] = prob;
    }

    void    CreateProduction(int i, int m, double prob) {
        itsIMProbs[i * numT + m] = prob;
    }

    void    CreateProduction(int i, int j, int k, double prob) {
        itsIJKProbs[i * numNT * numNT + j * numNT + k] = prob;
    }


    double  GetProb(int i, int m)   {
        return itsIMProbs[i * numT + m];
    }

    double  GetProb(int i, int j, int k) {
        return (itsIJKProbs[i*numNT*numNT + j*numNT + k]);
    }


    void    LinkProductions(int nlink, int* linkArray);

    /* list the production rules in this grammar */
    void    ShowProductions(bool detFlag);



    /* -------- TERMINAL PARSE TREE FUNCTIONS ---------- */

    void        SetParseTree(char * terminals);
    String *    ParseString(char * s);


    /* -------- MEMBER STRING FUNCTIONS ------------- */

    /* for non-member access to a string */
    String *        AccessString()      {
        return itsString;
    }

    /* start a new string */
    void    InitiateString() {
        itsString = new String(0, FALSE, numNT, numT);
    }
    void    InitiateString(int nt) {
        itsString = new String(nt, FALSE, numNT, numT);
    }


    /* converts non-terminal to terminal or other non-
        terminals symbols using production rules */
    bool    EmitSymbol(String * theString, int s);

    /* applies ResolveElement to all non-terminals in string */
    void    GenerateRandomString(bool verbose)  {
        GenerateRandomString(verbose, itsString);
    }
    void    GenerateRandomString(bool verbose, String * thisString);

    void    ReadStringsFromFile(char * filename);



    /* --------- CORPUS FUNCTIONS ------------------ */

    /* funcitons for managing corpus (group of strings) */
    void        AddStringToCorpus() {
        itsCorpus.push_back(itsString);
    }

    void        AddStringToCorpus(String * string) {
        itsCorpus.push_back(string);
    }

    String *    GetStringFromCorpus(int i)      {
        return itsCorpus.at(i);
    }

    void    ClearCorpus()       {
        itsCorpus.clear();
    }

    void    PrintCorpus();

    int     GetCorpusSize()     {
        return itsCorpus.size();
    }

    void    WriteCorpusToFile(char * filename);



    /* -------- RANDOM NUMBER GENERATOR FUNCTIONS ---------- */

    /* uniform random number generator */
    float   Urn()           {
        return pMersenne->Random();
    }

    /* for non-member access to the random-number generator */
    TRandomMersenne *   AccessMersenne()    {
        return pMersenne;
    }



    /* --------- PRODUCTION PROBABILITY FUNCTIONS ---------- */

    /* calculate inside probability of sub-string */
    double  InsideProb(String * thisString, int i, int s, int t);

    /* calculate outside probability of sub-string */
    double  OutsideProb(String * thisString, int i, int s, int t);

    /* initial calculation of inside and outside probabilities */
    void    CalculateProbs(String * stringPtr);

    // initialize storage array
    void    InitStringProbs(String * thisString);

    // release allocated memory and clear stored probabilities
    void    ReleaseEvals();

    char    GetStringProb(String * thisString, int i, int s, int t);

    /* returns evaluated and stored probabilities not equal to 0 or 1 */
    double  LookUpProb(bool isInner, int i, int s, int t);
    void    ShowProbs(bool inside);

    /* print out all stored probabilities not equal to 0 or 1 */
    void    PrintNonZeroProbs(bool inside);
    void    PrintStringProbs(String * thisString);
    void    PrintStringProbs(int s)     {
        PrintStringProbs(itsCorpus.at(s));
    }



    /* --------- PARSING FUNCTIONS ------------------------ */

    /* get likelihood of a string */
    double  GetLikString(String * stringPtr);

    double  GetLikCorpus();

    /* train production rule probabilities given a corpus */
    double  TrainCorpus();

    void    BestParseCorpus();

    void    CykFillMatrices(String *, double *, Triplet *);
    void    CykTraceback(String *, Triplet *);
#endif

#ifdef      _USE_HYPHY_HOOKS_               // SLKP June 2006

    /* SLKP SCFG data members */

    _List               terminals,          // stores the list of terminal strings. Allows the map from an
                        // internal terminal index to the actual literal
                        // as terminals[index]

                        rules,              // a list containing processed rules - lists as well, with
                        // 2 integer entries for NT->T rules
                        // 3 integer entries for NT->NT NT rules

                        byNT3,              // a list containing a _SimpleList for each non-terminal
                        // Each _SimpleList is a collection of NT->NT NT rules (indexed into 'rules')
                        // which have the corresponding NT on the LHS

                        byNT2,              // a list containing a _SimpleList for each non-terminal
                        // Each _SimpleList is a collection of NT->Terminal rules (indexed into 'rules')
                        // which have the corresponding NT on the LHS

                        /* modification by AFYP on June 20, 2006 */
                        byRightNT1,         // a list containing a _SimpleList for each non-terminal NT
                        // Each _SimpleList is a collection of NT->NT NT rules (indexed into 'rules')
                        // which have the corresponding NT on the RHS, appearing first

                        byRightNT2,         // a list containing a _SimpleList for each non-terminal NT
                        // Each _SimpleList is a collection of NT->NT NT rules (indexed into 'rules')
                        // which have the corresponding NT on the RHS, appearing second

                        /* -----  end mod ----- */

                        corpusChar,         // a list of strings representing current string corpus
                        corpusInt,          // the same set of strings, but converted into terminal tokens
                        // each represented by a _SimpleList

                        insideProbs,        // A list of AVL tree for each string in the corpus.
                        // Each AVL tree indexed using s (start of the substring), t (end of the substring) and i - index of the nonterminal
                        // and stores the inside probability of generating this substring (s,t) from non-terminal i
                        //  -- (s,t,i) is mapped onto an integer using scfgIndexIntoAnArray (see below)
                        //  -- if (s,t,i) is not in the tree, then the inside probability is 0
                        //  -- if (s,t,i) is in the tree, the value associated with the key is either
                        //      -- an index into the appropriate entry of storedInsideP (>=0); these are the only values which will be
                        //      recomputed when rule probabilities change
                        //      -- if < 0, then the corresponding inside probability is identically 1 for all parameter values

                        insideProbsT,       // auxiliary list used by insideProbs
                        outsideProbsT,      // auxiliary list used by outsideProbs

                        outsideProbs,       // similar to InsideProbs, but used for storing outside probabilities

                        storedInsideP,      // storage containers for the 0<p<1 inside probabilities
                        storedOutsideP;     // storage containers for the 0<p<1 outside probabilities

    _SimpleList         ntToTerminalMap,    // an array which maps each pair (nt, terminal) to
                        // index of the rule nt->terminal if it exists
                        // -1, if no such rule exists

                        computeFlagsI,      // a _SimpleList treated as bit strings, and used to flag
                        // whether a particular (s,t,i) has been computed for inside probabilities
                        // if this list is not empty, then no computations on the current corpus have
                        // been done yet; after the first computation this (and the next) list are
                        // flushed
                        // It is allocated with enought space to accomodate the longest
                        // string

                        computeFlagsO,      // same, but for outside probabilities

                        firstArray,         // an NTxT array (stored as a flat vector), where (i,j)-th element is set to one iff
                        // i-th non terminal can be used to generate a substing beginning with the j-th terminal
                        // otherwise the element is 0 (i.e. {i} => j ....)

                        lastArray,          // an NTxT array (stored as a flat vector), where (i,j)-th element is set to one iff
                        // i-th non terminal can be used to generate a substing ending with the j-th terminal
                        // otherwise the element is 0 (i.e. {i} => .... j)

                        precursorArray,     // an NTxT array (stored as a flat vector), where (i,j)-th element is set to one iff
                        // i-th non terminal can be used to generate an incompletely resolved substing beginning
                        // with the j-th terminal, followed by the i-th non-terminal (i.e. {start} => ... j {i} ... )

                        followArray;        // an NTxT array (stored as a flat vector), where (i,j)-th element is set to one iff
    // i-th non terminal can be used to generate an incompletely resolved substing ending
    // with the i-th non-terminal, followed by the j-th terminal (i.e. {start} => ... {i} j ...)

    _Matrix             probabilities;      // This is a matrix array of formulae which describe production rule probabilities
    // For deterministic rules, the formulae are set to constant 1
    // the indexing for the rules is the same as in the 'rules' list


    node<long>**        parseTree;          // maintains a parse tree which maps character input (ASCII) to the set of
    // terminal characters expressed as integer indices.
    // Each node is associated with the branch which terminates in it
    // The 3 highest order bytes of the <long> data filed are only used for leaves,
    // and represent the index of the terminal symbol, whilst the lowest order byte
    // associates a character with the branch
    // The 'root' of the parseTree is represented by the array of 256
    // node<long>* for faster indexing

    long                startSymbol,
                        insideCalls,
                        outsideCalls;

protected:

    void        ClearParseTree  (void);
    /* SLKP: this function is used to clear the data structures holding the input parse tree */
    void        ProcessAFormula (_FString*, _List&, _SimpleList&, _String&);
    /* SLKP: utility function to process a probability expression */
    bool        CheckANT        (long,long,long, _AVLListX&, long);
    /* SLKP: utility function which checks conditions on rules involving non-terminals
             returning true if status flags for some of the non-terminals were modified*/
    _String*    TokenizeString  (_String&, _SimpleList&);
    /* SLKP: convert a string (1st argument) into a series of terminal tokens (stored into the 2nd argument)
             Returns nil if all is good, or an error string if something went wrong
    */

    void        DumpComputeStructures (void);
    // clear all computational structures associated with a corpus

    void        InitComputeStructures (void);
    // initialize compute structures for a new corpus
    // by populating appropriate data structures with empties

    _Parameter  ComputeInsideProb     (long, long, long, long, bool);
    // compute the inside probability for substring from s (arg1) to t (arg2) - both zero based -
    // in corpus string j (arg3) derived from non-terminal i (arg4). The bool flag shows whether or
    // not this is the first call into a given corpus and that computeFlagsI should be consulted
    // during computation

    _Parameter  ComputeOutsideProb    (long, long, long, long, bool, bool);
    // compute the outside probability for substring from s (arg1) to t (arg2) - both zero based -
    // in corpus string j (arg3) derived from non-terminal i (arg4). The FIRST bool flag shows whether or
    // not this is the first call for outside probabilities into a given corpus and that computeFlagsO
    // should be consulted during computation.  The SECOND bool indicates first call for inside
    // probabilities and that computeFlagsI should be consulted.

    long        indexNT_T             (long, long);
    // index (nt, term) pairs into ntToTerminalMap

    _Parameter  LookUpRuleProbability (long index) {
        return ((_Matrix*)probabilities.RetrieveNumeric())->theData[index];
    }

#endif
};

#ifdef      _USE_HYPHY_HOOKS_               // SLKP June 2006

inline  long    scfgIndexIntoAnArray            (long,long,long,long);
/* this function indexes the triple start_substring, end_substring,non-terminal index given the string length (last argument)
   into a linear array
*/

_SimpleList *   arrayIntoScfgIndex              (long,long);
/* this function indexes from a linear array to the triple start_substring, end_substring, non-terminal index
    given the string length -AFYP
*/

inline  bool    getIndexBit                             (long,long,long,long,_SimpleList&);
inline  void    setIndexBit                             (long,long,long,long,_SimpleList&);
/* this function indexes the triple start_substring, end_substring,non-terminal index given the string length (4th argument)
   into a bit array stored in the last argument
*/

#endif

#define    _HY_BITMASK_WIDTH_ (8*sizeof (unsigned long))

struct      bitMasks {
    unsigned long masks[_HY_BITMASK_WIDTH_];
    bitMasks (void) {
        unsigned long aBit = 1;
        for (long k=0; k<_HY_BITMASK_WIDTH_; k++) {
            masks[k] = aBit;
            aBit = aBit << 1;
        }
    }
};

extern bitMasks bitMaskArray;

#endif
