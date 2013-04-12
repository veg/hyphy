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

#include      "alignment.h"
#include      "trie.h"
#include      "likefunc.h"
#include      "scfg.h"
#include      <ctype.h>

#include      "bayesgraph.h"
#include      "hy_globals.h"
#include      "executionlist.h"

#ifndef __HYPHY_NO_SQLITE__
#include "sqlite3.h"
#endif

#if !defined __HEADLESS__ && !defined __UNIX__
#include      "HYUtils.h"
#endif

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif



//____________________________________________________________________________________
// global variables

_String     sqlOpen                 ("SQL_OPEN"),
            sqlClose                ("SQL_CLOSE"),
            sqlRowData              ("SQL_ROW_DATA"),
            sqlColNames             ("SQL_COLUMN_NAMES"),
            seqAlignMap             ("SEQ_ALIGN_CHARACTER_MAP"),
            seqAlignScore           ("SEQ_ALIGN_SCORE_MATRIX"),
            seqAlignScoreCodon2     ("SEQ_ALIGN_SCORE_MATRIX_PARTIAL_CODON_2"),
            seqAlignScoreCodon1     ("SEQ_ALIGN_SCORE_MATRIX_PARTIAL_CODON_1"),
            seqAlignGapChar         ("SEQ_ALIGN_GAP_CHARACTER"),
            seqAlignGapOpen         ("SEQ_ALIGN_GAP_OPEN"),
            seqAlignGapExtend       ("SEQ_ALIGN_GAP_EXTEND"),
            seqAlignGapOpen2        ("SEQ_ALIGN_GAP_OPEN2"),
            seqAlignGapExtend2      ("SEQ_ALIGN_GAP_EXTEND2"),
            seqAlignFrameShift      ("SEQ_ALIGN_FRAMESHIFT"),
            seqAlignGapLocal        ("SEQ_ALIGN_NO_TP"),
            seqAlignGapAffine       ("SEQ_ALIGN_AFFINE"),
            seqAlignCodonAlign      ("SEQ_ALIGN_CODON_ALIGN"),
            seqAlignGapLinearSpace  ("SEQ_ALIGN_LINEAR_SPACE"),
            seqAlignGapCodon3x1     ("SEQ_ALIGN_PARTIAL_3x1_SCORES"),
            seqAlignGapCodon3x2     ("SEQ_ALIGN_PARTIAL_3x2_SCORES"),
            seqAlignGapCodon3x4     ("SEQ_ALIGN_PARTIAL_3x4_SCORES"),
            seqAlignGapCodon3x5     ("SEQ_ALIGN_PARTIAL_3x5_SCORES"),
            seqAlignDoLocal         ("SEQ_ALIGN_LOCAL_ALIGNMENT"),
            completeFlag            ("COMPLETE"),
            conditionalWeights      ("WEIGHTS"),
            siteProbabilities       ("SITE_LOG_LIKELIHOODS"),
            lastSetOfConstraints    ("LAST_SET_OF_CONSTRAINTS"),
            deferConstrainAssignment("DEFER_CONSTRAINT_APPLICATION"),
            assertionBehavior       ("ASSERTION_BEHAVIOR"),
            _hyStatusConditionProbsMatrix
            ("Constructing Conditional Probabilities Matrix"),

            isDynamicGraph          ("BGM_DYNAMIC");


extern      _String                 blDoSQL,
            blAlignSequences,
            blGetNeutralNull,
            blHBLProfile,
            blDeleteObject,
            timeStamp,
            versionString,
            lastModelParameterList,
            blGetString,
            blRequireVersion,
            blAssert;

_SimpleList sqlDatabases,
            _HY_HBLCommandHelperAux;
            
_List        scfgList,
             scfgNamesList,
             bgmList,
             bgmNamesList,
             _HY_GetStringGlobalTypesAux;
            
_Trie        _HY_ValidHBLExpressions;

_AVLListX    _HY_GetStringGlobalTypes (&_HY_GetStringGlobalTypesAux),
             _HY_HBLCommandHelper     (&_HY_HBLCommandHelperAux);
             

//____________________________________________________________________________________

//int          _HYSQLCallBack                     (void* data,int callCount);
//int        _HYSQLBusyCallBack                 (void* exList,int cc,char** rd,char** cn);
// this is a ReplicateConstraint helper

//____________________________________________________________________________________

_HBLCommandExtras* _hyInitCommandExtras (const long cut, const long conditions, _String commandInvocation, const char sep, const bool doTrim, const bool isAssignment, const bool needsVerb, _SimpleList * conditionList) {
    
    struct _HBLCommandExtras * commandInfo           = new _HBLCommandExtras();
    commandInfo->cut_string                          = cut;
    if (conditions < 0 && conditionList) 
        commandInfo->extract_conditions              << *conditionList;
    else
        commandInfo->extract_conditions              << conditions;
    commandInfo->extract_condition_separator         = sep;
    commandInfo->do_trim                             = doTrim;
    commandInfo->is_assignment                       = isAssignment;
    commandInfo->needs_verb                          = needsVerb;
    commandInfo->command_invocation                  && & commandInvocation;

    return                                             commandInfo;
    
}

//____________________________________________________________________________________

void        _HBL_Init_Const_Arrays  (void)
{
    // init GetString lookups
    _HY_GetStringGlobalTypes.Insert(new _String("LikelihoodFunction"), HY_BL_LIKELIHOOD_FUNCTION);
    _HY_GetStringGlobalTypes.Insert(new _String("DataSet"), HY_BL_DATASET);
    _HY_GetStringGlobalTypes.Insert(new _String("DataSetFilter"), HY_BL_DATASET_FILTER);
    _HY_GetStringGlobalTypes.Insert(new _String("UserFunction"), HY_BL_HBL_FUNCTION);
    _HY_GetStringGlobalTypes.Insert(new _String("Tree"), HY_BL_TREE);
    _HY_GetStringGlobalTypes.Insert(new _String("SCFG"), HY_BL_SCFG);
    _HY_GetStringGlobalTypes.Insert(new _String("Variable"), HY_BL_VARIABLE);
	_HY_GetStringGlobalTypes.Insert(new _String("BayesianGraphicalModel"), HY_BL_BGM);


    _HY_ValidHBLExpressions.Insert ("function ",                            HY_HBL_COMMAND_FUNCTION);
    _HY_ValidHBLExpressions.Insert ("ffunction ",                           HY_HBL_COMMAND_FFUNCTION);
    _HY_ValidHBLExpressions.Insert ("return ",                              HY_HBL_COMMAND_RETURNSPACE);
    _HY_ValidHBLExpressions.Insert ("return(",                              HY_HBL_COMMAND_RETURNPAREN);
    _HY_ValidHBLExpressions.Insert ("if(",                                  HY_HBL_COMMAND_IF);
    _HY_ValidHBLExpressions.Insert ("else",                                 HY_HBL_COMMAND_ELSE);
    _HY_ValidHBLExpressions.Insert ("do{",                                  HY_HBL_COMMAND_DO);
    _HY_ValidHBLExpressions.Insert ("break;",                               HY_HBL_COMMAND_BREAK);
    _HY_ValidHBLExpressions.Insert ("continue;",                            HY_HBL_COMMAND_CONTINUE);
    _HY_ValidHBLExpressions.Insert ("#include",                             HY_HBL_COMMAND_INCLUDE);
    _HY_ValidHBLExpressions.Insert ("DataSet ",                             HY_HBL_COMMAND_DATA_SET);
    _HY_ValidHBLExpressions.Insert ("DataSetFilter ",                       HY_HBL_COMMAND_DATA_SET_FILTER);
    _HY_ValidHBLExpressions.Insert ("ConstructCategoryMatrix(",				HY_HBL_COMMAND_CONSTRUCT_CATEGORY_MATRIX);
    _HY_ValidHBLExpressions.Insert ("Tree ",                                HY_HBL_COMMAND_TREE);
    _HY_ValidHBLExpressions.Insert ("LikelihoodFunction ",					HY_HBL_COMMAND_LIKELIHOOD_FUNCTION);
    _HY_ValidHBLExpressions.Insert ("LikelihoodFunction3 ",					HY_HBL_COMMAND_LIKELIHOOD_FUNCTION_3);
    _HY_ValidHBLExpressions.Insert ("MolecularClock(",                      HY_HBL_COMMAND_MOLECULAR_CLOCK);
    _HY_ValidHBLExpressions.Insert ("fscanf(",                              HY_HBL_COMMAND_FSCANF);
    _HY_ValidHBLExpressions.Insert ("sscanf(",                              HY_HBL_COMMAND_SSCANF);
    _HY_ValidHBLExpressions.Insert ("ReplicateConstraint(",					HY_HBL_COMMAND_REPLICATE_CONSTRAINT);
    //_HY_ValidHBLExpressions.Insert ("Import(",                              HY_HBL_COMMAND_IMPORT);
    _HY_ValidHBLExpressions.Insert ("category ",                            HY_HBL_COMMAND_CATEGORY);
    _HY_ValidHBLExpressions.Insert ("Model ",                               HY_HBL_COMMAND_MODEL);
    _HY_ValidHBLExpressions.Insert ("ChoiceList(",                          HY_HBL_COMMAND_SET_CHOICE_LIST);
    _HY_ValidHBLExpressions.Insert ("OpenDataPanel(",                       HY_HBL_COMMAND_OPEN_DATA_PANEL);
    _HY_ValidHBLExpressions.Insert ("GetInformation(",                      HY_HBL_COMMAND_GET_INFORMATION);
    _HY_ValidHBLExpressions.Insert ("ExecuteCommands(",                     HY_HBL_COMMAND_EXECUTE_COMMANDS);
    _HY_ValidHBLExpressions.Insert ("ExecuteAFile(",                        HY_HBL_COMMAND_EXECUTE_A_FILE);
    _HY_ValidHBLExpressions.Insert ("LoadFunctionLibrary(",					HY_HBL_COMMAND_LOAD_FUNCTION_LIBRARY);
    _HY_ValidHBLExpressions.Insert ("OpenWindow(",                          HY_HBL_COMMAND_OPEN_WINDOW);
    _HY_ValidHBLExpressions.Insert ("SpawnLikelihoodFunction(",				HY_HBL_COMMAND_SPAWN_LIKELIHOOD_FUNCTION);
    _HY_ValidHBLExpressions.Insert ("FindRoot(",                            HY_HBL_COMMAND_FIND_ROOT);
    _HY_ValidHBLExpressions.Insert ("MPIReceive(",                          HY_HBL_COMMAND_MPI_RECEIVE);
    _HY_ValidHBLExpressions.Insert ("MPISend(",                             HY_HBL_COMMAND_MPI_SEND);
    _HY_ValidHBLExpressions.Insert ("GetDataInfo(",                         HY_HBL_COMMAND_GET_DATA_INFO);
    _HY_ValidHBLExpressions.Insert ("StateCounter(",                        HY_HBL_COMMAND_STATE_COUNTER);
    _HY_ValidHBLExpressions.Insert ("Integrate(",                           HY_HBL_COMMAND_INTEGRATE);
    _HY_ValidHBLExpressions.Insert ("DoSQL(",                               HY_HBL_COMMAND_DO_SQL);
    _HY_ValidHBLExpressions.Insert ("Topology ",                            HY_HBL_COMMAND_TOPOLOGY);
    _HY_ValidHBLExpressions.Insert ("AlignSequences(",                      HY_HBL_COMMAND_ALIGN_SEQUENCES);
    _HY_ValidHBLExpressions.Insert ("GetNeutralNull(",                      HY_HBL_COMMAND_GET_NEUTRAL_NULL);
    _HY_ValidHBLExpressions.Insert ("#profile",                             HY_HBL_COMMAND_PROFILE);
    _HY_ValidHBLExpressions.Insert ("SCFG ",                                HY_HBL_COMMAND_SCFG);
    _HY_ValidHBLExpressions.Insert ("NeuralNet ",                           HY_HBL_COMMAND_NEURAL_NET);
    _HY_ValidHBLExpressions.Insert ("BGM ",                                 HY_HBL_COMMAND_BGM);
    _HY_ValidHBLExpressions.Insert ("SimulateDataSet",                      HY_HBL_COMMAND_SIMULATE_DATA_SET);
/*
const long cut, const long conditions, const char sep, const bool doTrim, const bool isAssignment, const bool needsVerb, length options
*/

    _SimpleList lengthOptions;
    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_FOR, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("for(", HY_HBL_COMMAND_FOR,false),3, "for (<initialization>;<condition>;<increment>) {loop body}"));

    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_WHILE,
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("while(", HY_HBL_COMMAND_WHILE,false),1, "while (<condition>) {loop body}"));
                                                            

    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_SET_DIALOG_PROMPT, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("SetDialogPrompt(", HY_HBL_COMMAND_SET_DIALOG_PROMPT,false),
                                    1, 
                                    "SetDialogPrompt(<prompt string>);"));
    
    lengthOptions.Clear();lengthOptions.Populate (3,5,1);
    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_HARVEST_FREQUENCIES, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("HarvestFrequencies(", HY_HBL_COMMAND_HARVEST_FREQUENCIES,false),
                                    -1, 
                                    "HarvestFrequencies(<receptacle>, <DataSet or DataSetFilter>, <atom INTEGER>, <unit INTEGER <= atom>, <position aware 0 or 1>, [optional site partion], [optional sequence partition] (only for DataSetArguments)",
                                        ',',
                                        true,
                                        false,
                                        false,
                                        &lengthOptions));
    
    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_OPTIMIZE, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("Optimize(", HY_HBL_COMMAND_OPTIMIZE,false),
                                                                2, 
                                                                "Optimize (<receptacle>, <likelihood function/scfg/bgm>)",','));

    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_LFCOMPUTE, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("LFCompute(", HY_HBL_COMMAND_LFCOMPUTE,false),
                                                                2, 
                                                                "LFCompute (<likelihood function/scfg/bgm>,<LF_START_COMPUTE|LF_DONE_COMPUTE|receptacle>)",','));


    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_COVARIANCE_MATRIX, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("CovarianceMatrix(", HY_HBL_COMMAND_COVARIANCE_MATRIX,false),
                                                                2, 
                                                                "CovarianceMatrix (<receptacle>, <likelihood function/scfg/bgm>)",','));

     
    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_SELECT_TEMPLATE_MODEL, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("SelectTemplateModel(", HY_HBL_COMMAND_SELECT_TEMPLATE_MODEL,false),
                                    1, 
                                    "SelectTemplateModel(<DataSetFilter>);"));


    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_USE_MODEL, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("UseModel(", HY_HBL_COMMAND_USE_MODEL,false),
                                                                1, 
                                                                "UseModel (<model ID>)",','));

    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_SET_PARAMETER, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("SetParameter(", HY_HBL_COMMAND_SET_PARAMETER,false),
                                                                3, 
                                                                "SetParameter(<object>, <parameter index>, <value>)",','));


    lengthOptions.Clear();lengthOptions.Populate (2,1,1);
    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_ASSERT, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("assert(", HY_HBL_COMMAND_ASSERT,false),
                                    -1, 
                                    "assert (<statement>,[optional message on failure]>",
                                        ',',
                                        true,
                                        false,
                                        false,
                                        &lengthOptions));


    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_REQUIRE_VERSION, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("RequireVersion(", HY_HBL_COMMAND_REQUIRE_VERSION,false),
                                                                1, 
                                                                "RequireVersion (<version string>)",','));

    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_DELETE_OBJECT, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("DeleteObject(", HY_HBL_COMMAND_DELETE_OBJECT,false),
                                                                -1, 
                                                                "DeleteObject(<object 1> [optional ,<object 2>, <object 3>, ..., <object N>])",','));
    
    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_CLEAR_CONSTRAINTS, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("ClearConstraints(", HY_HBL_COMMAND_CLEAR_CONSTRAINTS,false),
                                                                -1, 
                                                                "ClearConstraints(<object 1> [optional ,<object 2>, <object 3>, ..., <object N>])",','));

    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_MOLECULAR_CLOCK, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("MolecularClock(", HY_HBL_COMMAND_MOLECULAR_CLOCK,false),
                                                                -2, 
                                                                "MolecularClock(tree or tree node, local variable 1 [optional ,<local variable 2>, ..., <local variable N>])",','));


    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_FPRINTF,
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("fprintf(", HY_HBL_COMMAND_FPRINTF,false),
                                                                -2,
                                                                "fprintf(stdout|MESSAGE_LOG|TEMP_FILE_NAME|PROMPT_FOR_FILE|file path, object 1 [optional ,<object 2>, ..., <object N>])",','));

    
    lengthOptions.Clear();lengthOptions.Populate (1,2,1); // 2
    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_EXPORT, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("Export(", HY_HBL_COMMAND_EXPORT,false),
                                                                -1, 
                                                                "Export (<string variable ID>, <object ID>)",
                                                                ',',
                                                                true,
                                                                false,
                                                                false,
                                                                &lengthOptions));


    lengthOptions.Clear();lengthOptions.Populate (2,2,1);
    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_GET_URL, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("GetURL(", HY_HBL_COMMAND_GET_URL,false),
                                    -1, 
                                    "GetURL (<receptacle>,<URL>[, SAVE_TO_FILE])",
                                        ',',
                                        true,
                                        false,
                                        false,
                                        &lengthOptions));

    lengthOptions.Clear();lengthOptions.Populate (2,3,1); // 3 or 4
    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_GET_STRING, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("GetString(", HY_HBL_COMMAND_GET_STRING,false),
                                    -1, 
                                    "GetString(<receptacle>,<object>,<index>,[optional <second index>])",
                                        ',',
                                        true,
                                        false,
                                        false,
                                        &lengthOptions));
                                        

    lengthOptions.Clear();lengthOptions.Populate (2,3,1); // 3 or 4
    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_DIFFERENTIATE, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("Differentiate(", HY_HBL_COMMAND_DIFFERENTIATE,false),
                                    -1, 
                                    "Differentiate(<receptacle>, <the expression to differentiate>, <variable to differentiate>[, number of times, default = 1])",
                                    ',',
                                    true,
                                    false,
                                    false,
                                    &lengthOptions));
                                        
// matrix global arrays


    _HY_MatrixRandomValidPDFs.Insert ("Dirichlet", _HY_MATRIX_RANDOM_DIRICHLET);
    _HY_MatrixRandomValidPDFs.Insert ("Gaussian", _HY_MATRIX_RANDOM_GAUSSIAN);
    _HY_MatrixRandomValidPDFs.Insert ("Wishart", _HY_MATRIX_RANDOM_WISHART);
    _HY_MatrixRandomValidPDFs.Insert ("InverseWishart", _HY_MATRIX_RANDOM_INVERSE_WISHART);
    _HY_MatrixRandomValidPDFs.Insert ("Multinomial", _HY_MATRIX_RANDOM_MULTINOMIAL);


}

//____________________________________________________________________________________
void         InsertVarIDsInList     (_AssociativeList* theList , _String theKey, _SimpleList& varIDs)
{
    _FString arrayKey (theKey, false);
    _Matrix *mxEntry = nil;

    if (varIDs.lLength) {
        _List     varNames;
        for (unsigned long i=0; i < varIDs.lLength; i++) {
            _Variable* v = LocateVar (varIDs.lData[i]);
            if (v) {
                varNames << v->GetName();
            }
        }
        mxEntry = new _Matrix (varNames);
    } else {
        mxEntry = new _Matrix;
    }

    checkPointer (mxEntry);
    theList->MStore (&arrayKey,mxEntry,false);
}

//____________________________________________________________________________________
void         InsertStringListIntoAVL    (_AssociativeList* theList , _String theKey, _SimpleList& stringsToPick, _List& theStrings)
{
    _FString arrayKey (theKey, false);
    _Matrix *mxEntry = nil;

    if (stringsToPick.lLength) {
        _List     theNames;
        for (unsigned long i=0; i < stringsToPick.lLength; i++) {
            _String * v = (_String*)theStrings (stringsToPick.lData[i]);
            if (v) {
                theNames << v;
            }
        }
        mxEntry = new _Matrix (theNames);
    } else {
        mxEntry = new _Matrix;
    }

    checkPointer (mxEntry);
    theList->MStore (&arrayKey,mxEntry,false);
}


//____________________________________________________________________________________

_Matrix *   CheckMatrixArg          (_String* mxName, bool onlyStrings)
{
    _Variable * mVar = FetchVar (LocateVarByName (*mxName));
    if (mVar && mVar->ObjectClass() == MATRIX) {
        _Matrix * mx = (_Matrix*)mVar->GetValue();
        if (onlyStrings && (!mx->IsAStringMatrix())) {
            return nil;
        }
        return mx;
    }
    return nil;
}

//____________________________________________________________________________________

_AssociativeList *   CheckAssociativeListArg (_String* mxName)
{
    _Variable * mVar = FetchVar (LocateVarByName (*mxName));
    if (mVar && mVar->ObjectClass() == ASSOCIATIVE_LIST) {
        return (_AssociativeList*)mVar->GetValue();
    }
    return nil;
}


//____________________________________________________________________________________
/*
int  _HYSQLBusyCallBack (void* data, int callCount)
{
    if (callCount > 100)
        return 0;

#ifdef __WINDOZE__
    Sleep  (1 + genrand_real2()*100.);
#else
    usleep (100 + genrand_real2()*100000.);
#endif
    return 1;
}*/

//____________________________________________________________________________________

bool    RecurseDownTheTree (_SimpleList& theNodes, _List& theNames, _List&theConstraints, _List& theParts, _SimpleList& partIndex)
{
    _SimpleList localNodes;

    node<long>* firstNode = (node<long>*)theNodes(0), *otherNode;
    bool        doThisOne = (firstNode->get_parent()!=nil), good = true;
    long        index, ind, i;

    /*if (!doThisOne)
    {
        BufferToConsole (_String((_String*)theParts.toStr()));
        NLToConsole();
        BufferToConsole (_String((_String*)partIndex.toStr()));
        NLToConsole();
    }*/

    // there are a few cases to consider
    for (ind = 1; ind<=firstNode->get_num_nodes(); ind++) { // have children nodes
        localNodes<< (long)firstNode->go_down(ind);
        for (index = 1; index<theNodes.lLength; index++) {
            otherNode = (node<long>*)theNodes(index);
            otherNode = otherNode->go_down(ind);
            if (!otherNode) {
                good = false;
                break;
            }
            localNodes<<(long)otherNode;
        }
        if (!good) {
            break;
        }
        good = RecurseDownTheTree (localNodes, theNames, theConstraints, theParts, partIndex);
        if (!good) {
            break;
        }
        localNodes.Clear();
    }

    // do this constraint now

    if (doThisOne&&good) { // not a root - so we apply the constraint
        _CalcNode*  firstCNode = (_CalcNode*)LocateVar (firstNode->get_data());
        _SimpleList goodVars;
        _List       otherGoodVars;
        _Variable* firstVar;

        ind = 0;

        while ((firstVar=firstCNode->GetIthIndependent(ind))) {
            for (index = 0; index<partIndex.lLength; index++) {
                if (partIndex.lData[index]==0) {
                    if (!firstVar->GetName()->EqualWithWildChar((_String*)theParts.lData[index],'?')) {
                        break;
                    }
                }
            }
            if (index==partIndex.lLength) {
                goodVars<<ind;
            }
            ind++;
        }

        for (i = 1; i<theNodes.lLength; i++) {
            otherNode = (node<long>*)theNodes(i);
            firstCNode = (_CalcNode*)LocateVar (otherNode->get_data());
            _SimpleList dummy;
            otherGoodVars && & dummy;
            long          theseInd = firstCNode->CountAll();
            _SimpleList   avVars;
            for (index = 0; index < theseInd; index++) {
                avVars << index;
            }
            for (index = 0; index<goodVars.countitems(); index++) {
                long j=0,k=0;
                bool found1 = false;
                for (k = 0; k<partIndex.lLength; k++)
                    if (partIndex.lData[k]==i) {
                        break;
                    }

                for (; j<avVars.lLength; j++) {
                    firstVar = firstCNode->GetIthParameter(avVars.lData[j]);
                    if (firstVar->GetName()->EqualWithWildChar((_String*)theParts.lData[k],'?')) {
                        (*(_SimpleList*)(otherGoodVars(i-1))) << avVars.lData[j];
                        avVars.Delete (j);
                        found1 = true;
                        break;
                    }
                }
                if (!found1) {
                    goodVars.Delete (index);
                    for (long ff = 0; ff < i-1; ff++) {
                        ((_SimpleList*)(otherGoodVars(i-1)))->Delete (index);
                    }
                    index--;
                }
            }
        }

        // now the constraints can be built

        for (index = 0; index < goodVars.lLength; index++) {
            _String newConstraint;
            for (ind = 0; ind < partIndex.lLength; ind++) {
                if (partIndex.lData[ind]<0) {
                    newConstraint = newConstraint & *(_String*)theParts(ind);
                } else {
                    otherNode = (node<long>*)theNodes(partIndex.lData[ind]);
                    _CalcNode*  CNode = (_CalcNode*)LocateVar (otherNode->get_data());

                    if (ind>0)
                        newConstraint = newConstraint &
                                        *(CNode->GetIthParameter((*(_SimpleList*)
                                                (otherGoodVars(partIndex.lData[ind]-1))).lData[index])->GetName());
                    else
                        newConstraint = newConstraint &
                                        *(CNode->GetIthIndependent(goodVars.lData[index])->GetName());
                }

            }
            theConstraints&& &newConstraint;
        }
    }

    if (!good) {
        _String errMsg (*(LocateVar(firstNode->get_data())->GetName())& " is incompatible with "&
                        (*LocateVar(((node<long>*)theNodes(index-1))->get_data())->GetName()) & " in call to ReplicateConstraint");
        WarnError (errMsg);
        return false;
    }

    return true;

}


//____________________________________________________________________________________

void    RetrieveModelComponents (long mid, _Matrix*& mm, _Matrix*& fv, bool & mbf)
{
    if (mid >=0 && mid < modelTypeList.lLength) {
        if (modelTypeList.lData[mid] == 0) {
            mm = (_Matrix*)FetchObjectFromVariableByTypeIndex(modelMatrixIndices.lData[mid],MATRIX);
        } else {
            mm = nil;
        }

        long fvi = modelFrequenciesIndices.lData[mid];
        fv = (_Matrix*)FetchObjectFromVariableByTypeIndex(fvi>=0?fvi:(-fvi-1),MATRIX);
        mbf = (fvi>=0);
    } else {
        mm = fv = nil;
        mbf = false;
    }
}

//____________________________________________________________________________________

void    RetrieveModelComponents (long mid, _Variable*& mm, _Variable*& fv, bool & mbf)
{
    if (mid >= 0 && modelTypeList.lData[mid] == 0) {
        mm = LocateVar(modelMatrixIndices.lData[mid]);
    } else {
        mm = nil;
    }

    long fvi = modelFrequenciesIndices.lData[mid];
    fv = LocateVar (fvi>=0?fvi:(-fvi-1));
    mbf = (fvi>=0);
}

//____________________________________________________________________________________

bool    IsModelReversible (long mid)
{
    _Matrix *m = nil,
             *f = nil;
    bool    mbf;
    RetrieveModelComponents (mid, m, f, mbf);
    if (m&&f) {
        return m->IsReversible(mbf?nil:f);
    }
    return false;
}

//____________________________________________________________________________________

bool    IsModelOfExplicitForm (long modelID) {
    if (modelID != HY_NO_MODEL) {
        return modelTypeList.lData[modelID] != 0;
    }
    return false;
}

//____________________________________________________________________________________

void    ScanModelForVariables        (long modelID, _AVLList& theReceptacle, bool inclG, long modelID2, bool inclCat)
{
    if (modelID != HY_NO_MODEL) {
        // standard rate matrix
        if (modelTypeList.lData[modelID] == 0) {
            ((_Matrix*) (LocateVar(modelMatrixIndices.lData[modelID])->GetValue()))->ScanForVariables2(theReceptacle,inclG,modelID2,inclCat);
        } else {
        // formula based
            // inclG was replaced with false in a previous commit. This caused problems in the optimizer and in
            // likelihood reporting (it was consistently worse than optimizer results)
            ((_Formula*)modelMatrixIndices.lData[modelID])->ScanFForVariables(theReceptacle, inclG, false, inclCat);
        }
    }
}

//____________________________________________________________________________________

_String _HYHBLTypeToText (long type) {
    _String result (128L,true);
    if (type & HY_BL_DATASET) {
        result << "DataSet|";
    }
    
    if (type & HY_BL_DATASET_FILTER) {
        result << "DataSetFilter|";
    }
    
    if (type & HY_BL_LIKELIHOOD_FUNCTION) {
        result << "LikelihoodFunction|";
    }
    
    if (type & HY_BL_SCFG) {
        result << "SCFG|";
    }
    
    if (type & HY_BL_BGM) {
        result << "BGM|";
    }
    
    if (type & HY_BL_MODEL) {
        result << "Model|";
    }
    
    if (type & HY_BL_HBL_FUNCTION) {
        result << "function|";
    }
    
    result.Finalize();
    result.Trim (0,result.sLength-2);
    return result;
}


//____________________________________________________________________________________

BaseRef _HYRetrieveBLObjectByName    (_String& name, long& type, long *index, bool errMsg, bool tryLiteralLookup)
{
    long loc = -1;
    if (type & HY_BL_DATASET) {
        loc = FindDataSetName (name);
        if (loc >= 0) {
            type = HY_BL_DATASET;
            if (index) {
                *index = loc;
            }
            return dataSetList (loc);
        }
    }

    if (type & HY_BL_DATASET_FILTER) {
        loc = FindDataSetFilterName (name);
        if (loc >= 0) {
            type = HY_BL_DATASET_FILTER;
            if (index) {
                *index = loc;
            }
            return dataSetFilterList (loc);
        }
    }

    if (type & HY_BL_LIKELIHOOD_FUNCTION) {
        loc = FindLikeFuncName (name);
        if (loc >= 0) {
            type = HY_BL_LIKELIHOOD_FUNCTION;
            if (index) {
                *index = loc;
            }
            return likeFuncList (loc);
        }
    }

    if (type & HY_BL_SCFG) {
        loc = FindSCFGName (name);
        if (loc >= 0) {
            type = HY_BL_SCFG;
            if (index) {
                *index = loc;
            }
            return scfgList (loc);
        }
    }

    if (type & HY_BL_BGM) {
        loc = FindBgmName (name);
        if (loc >= 0) {
            type = HY_BL_BGM;
            if (index) {
                *index = loc;
            }
            return bgmList (loc);
        }
    }

    if (type & HY_BL_MODEL) {
        loc = FindModelName(name);
        if (loc < 0 && (name.Equal (&lastModelParameterList) || name.Equal (&useLastModel))) {
            loc = lastMatrixDeclared;
        }
        if (loc >= 0) {
            type = HY_BL_MODEL;
            if (index) {
                *index = loc;
            }
            if (IsModelOfExplicitForm(loc)) {
                return (BaseRef)modelMatrixIndices.lData[loc];
            }
            return LocateVar (modelMatrixIndices.lData[loc]);
        }
    }

    if (type & HY_BL_HBL_FUNCTION) {
        loc = FindBFFunctionName(name);
        if (loc >= 0) {
            type = HY_BL_HBL_FUNCTION;
            if (index) {
                *index = loc;
            }
            return batchLanguageFunctions (loc);
        }
    }
    
    if (tryLiteralLookup) {
        _String nameIDRef = ProcessLiteralArgument(&name, nil);
        return _HYRetrieveBLObjectByName (nameIDRef, type, index, errMsg, false);
    }

    if (errMsg) {
        WarnError (_String ("'") & name & "' does not refer to an existing object of type " & _HYHBLTypeToText (type));
    }
    type = HY_BL_NOT_DEFINED;
    return nil;
}


