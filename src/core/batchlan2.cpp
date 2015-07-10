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

int          _HYSQLCallBack                     (void* data,int callCount);
//int        _HYSQLBusyCallBack                 (void* exList,int cc,char** rd,char** cn);
bool         RecurseDownTheTree                 (_SimpleList&, _List&, _List&, _List&, _SimpleList&);
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

bool      _ElementaryCommand::ExtractValidateAddHBLCommand (_String& current_stream,const long command_code,  _List* pieces, _HBLCommandExtras* command_spec, _ExecutionList& command_list)

{
    if (command_spec->is_assignment) {
        // TBA
    } else {
        // by default push all of the 'pieces' arguments to the "argument" list
        _ElementaryCommand *cv = new _ElementaryCommand (command_code);
        cv->addAndClean  (command_list, pieces, 0);        
    }
    
    return true;
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

bool    _ElementaryCommand::ConstructDoSQL (_String&source, _ExecutionList&target)
// syntax: DoSQL (dbID,action string|file name,<callback ID>)
{
    _List pieces;
    ExtractConditions (source,blDoSQL.sLength,pieces,',');
    if (pieces.lLength!=3) {
        WarnError (_String ("Expected syntax:")& blDoSQL &"(dbID|" & sqlOpen & '|' & sqlClose & ",transaction string|file name,callback ID for an SQL transaction|where to store DB numeric ID)");
        return false;
    }

    _ElementaryCommand * dsql = new _ElementaryCommand (53);
    dsql->addAndClean(target,&pieces,0);
    return true;
}

//____________________________________________________________________________________

bool    _ElementaryCommand::ConstructProfileStatement (_String&source, _ExecutionList&target)
// syntax: #profile START|PAUSE|RESUME|indetifier to dump in
{

    _List pieces;
    ExtractConditions (source,blHBLProfile.sLength+1,pieces,';');
    if (pieces.lLength!=2) {
        WarnError (_String ("Expected syntax:")& blHBLProfile &" START|PAUSE|RESUME|where to store)");
        return false;
    }

    _ElementaryCommand *sp = new _ElementaryCommand (58);
    sp->addAndClean(target,&pieces,0);
    return true;
}


//____________________________________________________________________________________

int  _HYSQLCallBack (void* exL,int cc, char** rd, char** cn)
{
    _ExecutionList * exList = (_ExecutionList *)exL;

    if (!terminateExecution)
        if (exList && cc && exList->lLength) {
            _List     rowData,
                      columnNames;

            for (long cnt = 0; cnt < cc; cnt++) {
                if (rd[cnt]) {
                    rowData.AppendNewInstance (new _String (rd[cnt]));
                } else {
                    rowData.AppendNewInstance (new _String);
                }

                if (cn[cnt]) {
                    columnNames.AppendNewInstance (new _String (cn[cnt]));
                } else {
                    columnNames.AppendNewInstance (new _String);
                }
            }


            _Matrix * rowDataM     = new _Matrix (rowData),
            * columnNamesM = new _Matrix (columnNames);

            if (!(rowDataM && columnNamesM)) {
                checkPointer (nil);
            }

            _Variable* rdv = CheckReceptacle (&sqlRowData, blDoSQL,false),
                       * cnv = CheckReceptacle (&sqlColNames, blDoSQL,false);

            rdv->SetValue (rowDataM,false);
            cnv->SetValue (columnNamesM,false);

            exList->Execute();

        }
    return 0;
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

void      _ElementaryCommand::ExecuteDataFilterCases (_ExecutionList& chain)
{
    chain.currentCommand++;

    _String dataObjectID = chain.AddNameSpaceToID(*(_String*)parameters(1));
    
    long dsID           = (parameters.lLength>2)?FindDataSetName (dataObjectID):-1;
    bool isFilter       = false;

    if (dsID == -1) {
        dsID = (parameters.lLength>2)?FindDataSetFilterName (dataObjectID):-1;
        if (dsID == -1) {
            _AssociativeList * numericFilter = (_AssociativeList*)FetchObjectFromVariableByType(&dataObjectID, ASSOCIATIVE_LIST);
            if (numericFilter) {
                _String errCode;

                long    categoryCount = 1;

                if (parameters.lLength > 2)
                    // have multiple categories
                {
                    categoryCount = (long) ProcessNumericArgument((_String*)parameters(2),nil);
                }

                _String namesKey ("FILTER_NAMES"),
                        dataKey  ("FILTER_ARRAYS"),
                        freqKey  ("FILTER_FREQS");

                _Matrix* sequenceNames = (_Matrix*)numericFilter->GetByKey (namesKey,MATRIX);
                _List seqNames;
                if (sequenceNames) {
                    sequenceNames->FillInList (seqNames);
                }
                if (!sequenceNames || seqNames.lLength == 0) {
                    errCode = _String("Expected a non-empty string matrix as the ") & namesKey & " argument in call to CreateFilter";
                } else {
                    _AssociativeList * dataList = (_AssociativeList*)numericFilter->GetByKey (dataKey,ASSOCIATIVE_LIST);
                    _Matrix          * freqList = (_Matrix*)numericFilter->GetByKey (freqKey,MATRIX);

                    if (dataList && freqList) {
                        _List       goodSeqs;
                        long        sitePatterns    = freqList->GetVDim(),
                                    categDim      = -1;

                        if (freqList->GetHDim() != 1 || sitePatterns < 1 || freqList->MatrixType() != 1 ) {
                            errCode = _String("Expected a non-empty numeric ROW matrix as the ") & freqKey & " argument in call to CreateFilter";
                        } else {
                            for (long k=0; k<seqNames.lLength; k=k+1) {
                                _Matrix * dataMx = (_Matrix*)dataList->GetByKey (k,MATRIX);
                                if (dataMx && dataMx->MatrixType() == 1 ) {
                                    if (categDim < 0) {
                                        categDim = dataMx->GetVDim();
                                        if (categDim < 1) {
                                            break;
                                        }
                                    } else if (dataMx->GetVDim() != categDim) {
                                        break;
                                    }
                                    if (dataMx->GetHDim () != sitePatterns*categoryCount) {
                                        break;
                                    }

                                    goodSeqs << dataMx;
                                    continue;
                                }
                                break;
                            }

                            if (goodSeqs.lLength == seqNames.lLength) {
                                _DataSet * dummyDS = new _DataSet;
                                dummyDS->SetNoSpecies (seqNames.lLength);
                                dummyDS->GetNames().Duplicate (&seqNames);
                                dummyDS->GetTheMap().Populate (sitePatterns,0,1);
                                errCode = (*(_String*)parameters(0)) & "_internal_ds";
                                dsID = FindDataSetName (errCode);
                                if (dsID < 0) {
                                    dataSetList         <<dummyDS;
                                    DeleteObject        (dummyDS);
                                    dataSetNamesList&&  &errCode;
                                } else {
                                    dataSetList.Replace (dsID,dummyDS,false);
                                }

                                errCode = (*(_String*)parameters(0));
                                _DataSetFilterNumeric * dsn = new _DataSetFilterNumeric (freqList,goodSeqs,dummyDS,categoryCount);
                                checkPointer (dsn);
                                dsID    = FindDataSetFilterName (errCode);

                                if (dsID < 0) {
                                    dataSetFilterList<<  dsn;
                                    DeleteObject        (dsn);
                                    dataSetFilterNamesList&& & errCode;
                                } else {
                                    dataSetFilterList.Replace (dsID,dsn,false);
                                }
                                return;

                            } else {
                                errCode = _String ("Site frequency patterns/numeric vectors did not pass dimension checks in call to CreateFilter");
                            }
                        }
                    }

                }
                if (errCode) {
                    WarnError(errCode);
                    return;
                }

            }
            _String errMsg = (((_String)("DataSet(Filter)/Associative Array ")&dataObjectID&_String(" has not been properly initialized")));
            WarnError (errMsg);
            return;
        }
        isFilter = true;
    }

    // build the formula from the 2nd parameter (unit size)

    char                unit = ProcessNumericArgument((_String*)parameters(2),chain.nameSpacePrefix);
    // here's our unit

    _String             dataFilterID (chain.AddNameSpaceToID(*(_String*)parameters(0))),
                        hSpecs,
                        vSpecs;

    long                status  = FindDataSetFilterName (dataFilterID);

    _DataSetFilter      *thedf;

    if (status!=-1) {
        thedf = (_DataSetFilter*)dataSetFilterList (status);
    } else {
        thedf               = new _DataSetFilter();
        checkPointer        (thedf);
        AddFilterToList     (dataFilterID,thedf,false);
    }

    if (parameters.lLength>3) {
        vSpecs = *(_String*)parameters(3);
    }
    if (parameters.lLength>4) {
        hSpecs = *(_String*)parameters(4);
    } else {
        hSpecs = empty;
    }

    _DataSet            *dataset;

    _SimpleList         hL,
                        vL;

    hL.RequestSpace (1024);
    vL.RequestSpace (1024);

    if (!isFilter) {
        dataset = (_DataSet*)dataSetList(dsID);
        dataset -> ProcessPartition (hSpecs,hL,false);
        if (code!=6 && vSpecs.sLength==0) {
            vSpecs = _String("0-")&_String(dataset->NoOfColumns()-1);
        }
        dataset->ProcessPartition (vSpecs,vL,true);
    } else {
        _DataSetFilter * dataset1 = (_DataSetFilter*)dataSetFilterList(dsID);
        dataset1->GetData()->ProcessPartition (hSpecs,hL,false, &dataset1->theNodeMap, &dataset1->theOriginalOrder);

        if (code!=6 && vSpecs.sLength==0) {
            vSpecs = _String("0-")&_String(dataset1->GetFullLengthSpecies()-1);
        }

        dataset1->GetData()->ProcessPartition (vSpecs,vL,true,  &dataset1->theOriginalOrder, &dataset1->theNodeMap);
        dataset = (_DataSet*)dataset1;
    }

    if (code!=6) {
        if (vL.lLength%unit) {
            vSpecs = (_String)"Unit size of "& unit & " doesn't divide the length of specified partition in call to ";
            if (code==27) { // Permute
                vSpecs = vSpecs & "Permute";
            } else {
                vSpecs = vSpecs & "Bootstrap";
            }

            vSpecs = vSpecs & ". The partition has been trimmed at the end.";
            ReportWarning (vSpecs);
            for (status = vL.lLength%unit; status>0; status--) {
                vL.Delete (vL.lLength-1);
            }
        }
        if (code == 27) {
            vL.Permute (unit);
        } else {
            vL.PermuteWithReplacement(unit);
        }

    }

    thedf->SetFilter (dataset, unit, hL, vL, isFilter);

    if (parameters.lLength>5) {
        hSpecs = GetStringFromFormula((_String*)parameters(5),chain.nameSpacePrefix);
        thedf->SetExclusions(&hSpecs);
    } else if ((code!=6)&&isFilter) {
        _DataSetFilter * df1 = (_DataSetFilter*)dataSetFilterList(dsID);
        if (df1->theExclusions.lLength) {
            thedf->theExclusions.Duplicate (&df1->theExclusions);
            thedf->SetDimensions();
        }
    }

    thedf->SetDimensions();
    thedf->SetupConversion();

    SetDataFilterParameters (dataFilterID, thedf, true);
}
//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase21 (_ExecutionList& chain)
{
    chain.currentCommand++;

    SetStatusLine (_hyStatusConditionProbsMatrix);
    _String errMsg,
            objectName    =     chain.AddNameSpaceToID(*(_String*)parameters(1)),
            resultID      =     chain.AddNameSpaceToID(*(_String*)parameters(0));

    long objectID         =     FindLikeFuncName (objectName, true);
    _PMathObj ob          =     nil;

    if (objectID >=0) { // likelihood function
        _Matrix * partitionList         = nil;
        if (parameters.lLength>3) {
            _String  secondArg = *(_String*)parameters(3);
            partitionList = (_Matrix*)ProcessAnArgumentByType (&secondArg, chain.nameSpacePrefix, MATRIX);
        }
        _SimpleList                     partsToDo;
        _LikelihoodFunction*            lf = (_LikelihoodFunction*)likeFuncList(objectID);
        if (lf->ProcessPartitionList(partsToDo, partitionList, _hyStatusConditionProbsMatrix)) {
            char runMode = _hyphyLFConstructCategoryMatrixConditionals;
            if (parameters.lLength > 2) {
                if (((_String*)parameters(2))->Equal(&completeFlag)) {
                    runMode = _hyphyLFConstructCategoryMatrixConditionals;
                } else if (((_String*)parameters(2))->Equal(&conditionalWeights)) {
                    runMode = _hyphyLFConstructCategoryMatrixWeights;
                } else if (((_String*)parameters(2))->Equal(&siteProbabilities)) {
                    runMode = _hyphyLFConstructCategoryMatrixSiteProbabilities;
                } else {
                    runMode = _hyphyLFConstructCategoryMatrixClasses;
                }
            }
            ob = lf->ConstructCategoryMatrix(partsToDo,runMode,true, &resultID);
        }
    } else {
        _TheTree * testTree = (_TheTree*) FetchObjectFromVariableByType (&objectName, TREE);
        if (testTree) {
            long    pid = 0;
            objectID = testTree->IsLinkedToALF (pid);
            if (objectID >= 0) {
                _LikelihoodFunction * anLF      = (_LikelihoodFunction*) likeFuncList (objectID);
                _DataSetFilter      * dsf       = (_DataSetFilter*) dataSetFilterList (anLF->GetTheFilters()(pid));
                anLF->PrepareToCompute();
                anLF->Compute         ();
                objectID                        = dsf->NumberDistinctSites();

                _Matrix             *condMx     = new _Matrix   (2*objectID*(testTree->GetLeafCount()
                        + testTree->GetINodeCount()) * testTree->categoryCount,
                        testTree->GetCodeBase(),
                        false, true);
                _List               leafNames,
                                    inodeNames;

                testTree->DepthWiseT(true);

                while (testTree->currentNode) {
                    _String       * bs = new _String;
                    testTree->GetNodeName  (testTree->currentNode, *bs);
                    if (testTree->IsCurrentNodeATip()) {
                        leafNames << bs;
                    } else {
                        inodeNames << bs;
                    }
                    DeleteObject (bs);
                    testTree->DepthWiseT(false);
                }

                leafNames << inodeNames;

                _Matrix  *nodeNames = new _Matrix (leafNames);

                for (long siteC = 0; siteC < objectID; siteC ++) {
                    testTree->RecoverNodeSupportStates (dsf,siteC,siteC-1,*condMx);
                }

                anLF->DoneComputing   ();
                _AssociativeList *retMe = new _AssociativeList;
                retMe->MStore ("Nodes",nodeNames,false);
                retMe->MStore ("Values",condMx,false);
                ob = retMe;
            }
        }
    }

    if (ob) {
        CheckReceptacleAndStore (&resultID, blConstructCM, true, ob, false);
    } else {
        WarnError (objectName & " must be either a likelihood function or a tree variable tied to a likelihood function.");
    }

}


//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase53 (_ExecutionList& chain)
{
    chain.currentCommand++;

#ifdef __HYPHY_NO_SQLITE__
    _String errStr ("SQLite commands can not be used in a HyPhy version built with the __HYPHY_NO_SQLITE__ flag");
    WarnError (errStr);
#else

    _String arg1  (*(_String*)parameters(0));

    char  * errMsg = nil;
    _String errStr;

    if (arg1.Equal (&sqlOpen)) {
        _Variable * dbVar = CheckReceptacle ((_String*)parameters(2), blDoSQL);

        if (dbVar) {
            _String arg2 (*(_String*)parameters(1));
            arg2.ProcessFileName(true,true,(Ptr)chain.nameSpacePrefix);
            int errCode  = SQLITE_OK;
            sqlite3 *aDB = nil;
#ifdef __HYPHYXCODE__
            errCode = sqlite3_open (DoMacToPOSIX(arg2).getStr(),&aDB);
#else
            errCode = sqlite3_open (arg2.sData,&aDB);
#endif
            if (errCode == SQLITE_OK) {
                errCode = sqlite3_exec(aDB, "SELECT COUNT(*) FROM sqlite_master", _HYSQLCallBack, nil, nil);
            }
            if (errCode != SQLITE_OK) {
                WarnError (sqlite3_errmsg(aDB));
                sqlite3_close(aDB);
                return;
            } else {
                long f = sqlDatabases.Find (0);
                if (f<0) {
                    f = sqlDatabases.lLength;
                    sqlDatabases << (long)aDB;
                } else {
                    sqlDatabases.lData[f] = (long)aDB;
                }

                sqlite3_busy_timeout (aDB, 5000);

                dbVar->SetValue (new _Constant (f), false);
            }
        }
    } else {
        bool doClose =  arg1.Equal (&sqlClose);

        long dbIdx = ProcessNumericArgument (doClose?(_String*)parameters(2):&arg1,chain.nameSpacePrefix);

        if (dbIdx<0 || dbIdx >= sqlDatabases.lLength || sqlDatabases.lData[dbIdx] == 0) {
            errStr = _String(dbIdx) & " is an invalid database index";
        } else {
            if (doClose) {
                sqlite3_close ((sqlite3*)sqlDatabases.lData[dbIdx]);
                sqlDatabases.lData[dbIdx] = 0;
            } else {
                _String arg3 (ProcessLiteralArgument((_String*)parameters(2),chain.nameSpacePrefix));

                _ExecutionList sqlProcessor (arg3,chain.nameSpacePrefix?(chain.nameSpacePrefix->GetName()):nil);
                if (!terminateExecution) {
                    _String arg2 (ProcessLiteralArgument ((_String*)parameters(1),chain.nameSpacePrefix));

                    if (sqlite3_exec((sqlite3*)sqlDatabases.lData[dbIdx], arg2.sData, _HYSQLCallBack, (Ptr)&sqlProcessor, &errMsg) != SQLITE_OK) {
                        WarnError (sqlite3_errmsg((sqlite3*)sqlDatabases.lData[dbIdx]));
                        return;
                    }
                }
            }
        }

    }

    if (errStr.sLength) {
        errStr = errStr & " in call to DoSQL";
        WarnError (errStr);
    }

#endif
}

//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase54 (_ExecutionList& chain)
{
    chain.currentCommand++;

    SetStatusLine (_String("Constructing Topology ")&*(_String*)parameters(0));

    _String  *treeSpec = ((_String*)parameters(1));
    treeSpec->ProcessParameter();
    _TreeTopology * tr = nil;

    if (treeSpec->sLength) {
        if (treeSpec->sData[0]!='(') {
            _Variable* testTree = FetchVar(LocateVarByName (AppendContainerName(*treeSpec,chain.nameSpacePrefix)));
            if (testTree && testTree->ObjectClass () == TREE) {
                tr = new _TreeTopology ((_TheTree*)testTree);
            } else {
                _String   flaData (*treeSpec);
                _Formula  nameForm (flaData,chain.nameSpacePrefix);
                _PMathObj formRes = nameForm.Compute();
                if (formRes&&formRes->ObjectClass () == STRING)
                    tr = new _TreeTopology (AppendContainerName(*(_String*)parameters(0),chain.nameSpacePrefix),
                                            *((_FString*)formRes)->theString,
                                            false);
            }
        } else
            tr = new _TreeTopology (AppendContainerName(*(_String*)parameters(0),chain.nameSpacePrefix),
                                    *(_String*)parameters(1),false);
    }

    if (!tr) {
        WarnError ("Illegal right hand side in call to Topology id = ...; it must be a string, a Newick tree spec or a topology");
    }
}


//____________________________________________________________________________________

bool    _ElementaryCommand::ConstructAlignSequences (_String&source, _ExecutionList&target)
// syntax: AlignSequences (result, input string matrix,  options matrix)
{
    _List pieces;
    ExtractConditions (source,blAlignSequences.sLength,pieces,',');
    if (pieces.lLength!=3) {
        WarnError ("Expected syntax: AlignSequences(result, input string matrix, options list);");
        return false;
    }

    _ElementaryCommand * as = new _ElementaryCommand (55);
    as->addAndClean(target,&pieces,0);
    return true;
}

//____________________________________________________________________________________

bool    _ElementaryCommand::ConstructGetNeutralNull (_String&source, _ExecutionList&target)
// syntax: GetNeutralNull (result, likelihood function, syn sub count matrix, non-syn sub count matrix, iterations per root state)
{
    _List pieces;
    ExtractConditions (source,blGetNeutralNull.sLength,pieces,',');
    if (pieces.lLength!=5) {
        WarnError ("Expected syntax: GetNeutralNull (result, likelihood function, syn sub count matrix, non-syn sub count matrix, iterations per root state);");
        return false;
    }

    _ElementaryCommand * gnn = new _ElementaryCommand (57);
    gnn->addAndClean(target,&pieces,0);
    return true;
}

//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase55 (_ExecutionList& chain)
{
    chain.currentCommand++;

    _String errStr;

    _Variable * storeResultIn = CheckReceptacle (&AppendContainerName(*(_String*)parameters(0),chain.nameSpacePrefix), blAlignSequences, true);

    if (storeResultIn) {
        _Matrix * inStrings = CheckMatrixArg (&AppendContainerName(*(_String*)parameters(1),chain.nameSpacePrefix),true);
        if (inStrings && (inStrings->GetHDim()==1||inStrings->GetVDim()==1)) {
            _AssociativeList * mappingTable = CheckAssociativeListArg (&AppendContainerName(*(_String*)parameters(2),chain.nameSpacePrefix));
            if (mappingTable) {
                // check for required parameters

                _FString * charVector = (_FString*)mappingTable->GetByKey (seqAlignMap, STRING);

                long       charCount = 0;

                _SimpleList ccount (256,-1,0);

                if (charVector) {
                    for (long cc = 0; cc < charVector->theString->sLength; cc++)
                        if (ccount.lData[(unsigned char)charVector->theString->sData[cc]]>=0) {
                            charCount = 0; // this is an error condition for
                            // duplicate characters in the string
                            break;
                        } else {
                            ccount.lData[(unsigned char)charVector->theString->sData[cc]] = cc;
                            charCount ++;
                        }
                }

                if (charVector && charCount) {
                    // now check that all characters
                    bool        doLocal      = false,
                                doAffine     = false,
                                doLinear     = true,
                                doCodon      = false,
                                doFullLocal  = false;


                    long        codonCount = charCount*charCount*charCount;

                    _PMathObj   c = mappingTable->GetByKey (seqAlignCodonAlign, NUMBER);
                    if (c) {
                        doCodon = c->Compute()->Value() > 0.5;
                    }


                    _Matrix * scoreMatrix = (_Matrix*)mappingTable->GetByKey (seqAlignScore, MATRIX);
                    if (scoreMatrix && scoreMatrix->GetHDim () == (doCodon?codonCount+1:charCount) && scoreMatrix->GetVDim () == scoreMatrix->GetHDim ()) {
                        scoreMatrix = (_Matrix*)scoreMatrix->ComputeNumeric();
                        scoreMatrix->CheckIfSparseEnough(true);

                        char        gapCharacter = '-';
                        _FString    *gapC = (_FString*)mappingTable->GetByKey (seqAlignGapChar, STRING);

                        _Matrix     *codon3x5 = nil,
                                    *codon3x4 = nil,
                                    *codon3x2 = nil,
                                    *codon3x1 = nil;

                        if (doCodon) {
                            codon3x5 = (_Matrix*)mappingTable->GetByKey (seqAlignGapCodon3x5, MATRIX);
                            codon3x4 = (_Matrix*)mappingTable->GetByKey (seqAlignGapCodon3x4, MATRIX);
                            codon3x2 = (_Matrix*)mappingTable->GetByKey (seqAlignGapCodon3x2, MATRIX);
                            codon3x1 = (_Matrix*)mappingTable->GetByKey (seqAlignGapCodon3x1, MATRIX);
                            if ( codon3x5 && codon3x4 && codon3x2 && codon3x1
                              && codon3x5->GetHDim() == codonCount+1
                              && codon3x4->GetHDim() == codonCount+1
                              && codon3x2->GetHDim() == codonCount+1
                              && codon3x1->GetHDim() == codonCount+1
                              && codon3x5->GetVDim() == charCount*charCount*charCount*10
                              && codon3x4->GetVDim() == charCount*charCount*charCount*4
                              && codon3x2->GetVDim() == charCount*charCount*3
                              && codon3x1->GetVDim() == charCount*3) {
                                codon3x5 = (_Matrix*)codon3x5->ComputeNumeric();
                                codon3x5 -> CheckIfSparseEnough(true);
                                codon3x4 = (_Matrix*)codon3x4->ComputeNumeric();
                                codon3x4 -> CheckIfSparseEnough(true);
                                codon3x2 = (_Matrix*)codon3x2->ComputeNumeric();
                                codon3x2 -> CheckIfSparseEnough(true);
                                codon3x1 = (_Matrix*)codon3x1->ComputeNumeric();
                                codon3x1-> CheckIfSparseEnough(true);
                            } else {
                                errStr = ( seqAlignGapCodon3x5 & ", "
                                         & seqAlignGapCodon3x4 & ", "
                                         & seqAlignGapCodon3x2 & ", or "
                                         & seqAlignGapCodon3x1 & " matrices are missing or have incorrect dimensions" );
                            }

                        }

                        if (errStr.sLength == 0) {
                            _String     settingReport (128L,true);

                            settingReport << "Running sequence alignment with the following options:";

                            if (gapC && gapC->theString->sLength == 1) {
                                gapCharacter = gapC->theString->sData[0];
                            }

                            settingReport << "\n\tGap character:";
                            settingReport << gapCharacter;

                            _Parameter  gapOpen       = 15.,
                                        gapOpen2      = 15.,
                                        gapExtend     = 1.,
                                        gapExtend2    = 1.,
                                        gapFrameshift = 50.;




                            c = mappingTable->GetByKey (seqAlignGapOpen, NUMBER);
                            if (c) {
                                gapOpen = c->Compute()->Value();
                            }

                            settingReport << "\n\tGap open cost:";
                            settingReport << _String (gapOpen);

                            gapOpen2 = gapOpen;
                            c = mappingTable->GetByKey (seqAlignGapOpen2, NUMBER);
                            if (c) {
                                gapOpen2 = c->Compute()->Value();
                            }

                            settingReport << "\n\tGap open cost 2:";
                            settingReport << _String (gapOpen2);

                            c = mappingTable->GetByKey (seqAlignGapExtend, NUMBER);
                            if (c) {
                                gapExtend = c->Compute()->Value();
                            }

                            settingReport << "\n\tGap extend cost:";
                            settingReport << _String (gapExtend);

                            gapExtend2 = gapExtend;
                            c = mappingTable->GetByKey (seqAlignGapExtend2, NUMBER);
                            if (c) {
                                gapExtend2 = c->Compute()->Value();
                            }

                            settingReport << "\n\tGap extend cost 2:";
                            settingReport << _String (gapExtend2);

                            c = mappingTable->GetByKey (seqAlignFrameShift, NUMBER);
                            if (c) {
                                gapFrameshift = c->Compute()->Value();
                            }

                            settingReport << "\n\tCodon frameshift cost:";
                            settingReport << _String (gapFrameshift);


                            c = mappingTable->GetByKey (seqAlignGapLocal, NUMBER);
                            if (c) {
                                doLocal = c->Compute()->Value() > 0.5;
                            }

                            settingReport << "\n\tIgnore terminal gaps: ";
                            settingReport << (doLocal?"Yes":"No");

                            settingReport << "\n\tUse codon alignment with frameshift routines: ";
                            if (doCodon) {
                                for (long i = 0; i < 256; i++)
                                    if (ccount.lData[i] < 0) {
                                        ccount.lData[i] = -codonCount - 1;
                                    }
                                settingReport << "Yes";
                                doLinear = false;

                                settingReport << "\n\t Linear space routines  are not implemented";
                            } else {
                                settingReport << "No";
                            }
                            
                            c = mappingTable->GetByKey (seqAlignDoLocal, NUMBER);
                            if (c) {
                                doFullLocal = c->Compute()->Value ()>0.5;
                            }
                            settingReport << "\n\tLocal alignment: ";
                            settingReport << (doFullLocal?"Yes":"No");
                            if (!doCodon && doFullLocal) {
                                 settingReport << "\n\t Local alignment is currently available for the codon aligner only.";                           
                            }

                            c = mappingTable->GetByKey (seqAlignGapAffine, NUMBER);
                            if (c) {
                                doAffine = c->Compute()->Value() > 0.5;
                            }
                            settingReport << "\n\tAffine gap costs: ";
                            settingReport << (doAffine?"Yes":"No");

                            c = mappingTable->GetByKey (seqAlignGapLinearSpace, NUMBER);
                            if (c) {
                                doLinear = c->Compute()->Value() > 0.5;
                            }

                            settingReport << "\n\tUse linear space routines: ";
                            settingReport << (doLinear?"Yes":"No");

                            settingReport.Finalize();
                            ReportWarning (settingReport);

                            long stringCount = inStrings->GetHDim() * inStrings->GetVDim();

                            _AssociativeList *alignedStrings = new _AssociativeList;
                            checkPointer (alignedStrings);


                            for (long s1 = 0; s1 < stringCount; s1++) {
                                _String*  str1 = ((_FString*)inStrings->GetFormula(0,s1)->Compute())->theString;
                                if (!str1) {
                                    errStr = _String("The ") & (s1+1) & "-th argument is not a string";
                                    break;
                                }
                                for (long s2 = s1+1; s2 < stringCount; s2++) {
                                    _String       *string2 = ((_FString*)inStrings->GetFormula(0,s2)->Compute())->theString;
                                    if (!string2) {
                                        errStr = _String("The ") & (s2+1) & "-th argument is not a string";
                                        break;
                                    }
                                    _AssociativeList * pairwiseComp = new _AssociativeList;
                                    checkPointer (pairwiseComp);

                                    _Parameter    score = 0.0;

                                    if (doLinear == false) {
                                        char * str1r = NULL,
                                             * str2r = NULL;
                                        _List         store;
                                        score = AlignStrings (str1->sData,string2->sData,str1r,str2r,ccount.lData,scoreMatrix->theData,scoreMatrix->GetVDim(),
                                                              gapCharacter,gapOpen,gapExtend,gapOpen2,gapExtend2,gapFrameshift,doLocal,doAffine,doCodon,
                                                              charCount, codon3x5->theData, codon3x4->theData, codon3x2->theData, codon3x1->theData, doFullLocal);

                                        if ( str1r && str2r ) {
                                            _String * r_res = ( _String * ) checkPointer( new _String( str1r ) ),
                                                    * q_res = ( _String * ) checkPointer( new _String( str2r ) );
                                            delete [] str1r;
                                            delete [] str2r;
                                            r_res->Finalize();
                                            q_res->Finalize();
                                            store.AppendNewInstance( r_res );
                                            store.AppendNewInstance( q_res );
                                        } else
                                            WarnError( "Internal Error in AlignStrings" );

                                        store.bumpNInst();

                                        if (store.lLength == 0) {
                                            errStr = "Unspecified error in AlignStrings";
                                            DeleteObject (pairwiseComp);
                                            s1 = stringCount;
                                            break;
                                        } else {
                                            pairwiseComp->MStore ("1", new _FString((_String*)store(0)), false);
                                            pairwiseComp->MStore ("2", new _FString((_String*)store(1)), false);
                                        }
                                    } else {
                                        _Matrix       scoreM        (string2->sLength+1,1,false,true),
                                                      scoreM2       (string2->sLength+1,1,false,true),
                                                      gap1Matrix    (string2->sLength+1,1,false,true),
                                                      gap2Matrix    (string2->sLength+1,1,false,true),
                                                      gap1Matrix2   (string2->sLength+1,1,false,true),
                                                      gap2Matrix2   (string2->sLength+1,1,false,true),
                                                      *buffers[6];

                                        char          *alignmentRoute = new char[2*(string2->sLength+1)];

                                        alignmentRoute[0] = alignmentRoute[string2->sLength+1] = 0;
                                        buffers[0] = &scoreM;
                                        buffers[1] = &gap1Matrix;
                                        buffers[2] = &gap2Matrix;
                                        buffers[3] = &scoreM2;
                                        buffers[4] = &gap1Matrix2;
                                        buffers[5] = &gap2Matrix2;
                                        _SimpleList ops (str1->sLength+2,-2,0);
                                        ops.lData[str1->sLength+1] = string2->sLength;
                                        ops.lData[0]               = -1;

                                        score = LinearSpaceAlign(str1,string2,ccount,scoreMatrix,
                                                                 gapOpen,gapExtend,gapOpen2,gapExtend2,
                                                                 doLocal,doAffine,ops,score,0,
                                                                 str1->sLength,0,string2->sLength,buffers,0,alignmentRoute);

                                        delete[]    alignmentRoute;

                                        _String     *result1 = new _String (str1->sLength+1, true),
                                        *result2 = new _String (string2->sLength+1, true);

                                        long        last_column     = ops.lData[ops.lLength-1];

                                        for (long position = str1->sLength-1; position>=0; position--) {
                                            long current_column     = ops.lData[position+1];

                                            if (current_column<0) {
                                                if (current_column == -2 /*|| (current_column == -3 && last_column == string2->sLength)*/) {
                                                    current_column = last_column;
                                                } else if (current_column == -3) {
                                                    // find the next matched char or a -1
                                                    long    p   = position,
                                                            s2p;
                                                    while ((ops.lData[p+1]) < -1) {
                                                        p--;
                                                    }

                                                    s2p = ops.lData[p+1];
                                                    //if (last_column == string2->sLength)
                                                    //  last_column = string2->sLength-1;

                                                    //if (s2p < 0)
                                                    //  s2p = 0;

                                                    for (long j = last_column-1; j>s2p;) {
                                                        (*result1) << gapCharacter;
                                                        (*result2) << string2->sData[j--];
                                                    }

                                                    last_column     = s2p+1;

                                                    for (; position>p; position--) {
                                                        (*result2) << gapCharacter;
                                                        (*result1) << str1->sData[position];
                                                    }
                                                    position ++;
                                                    continue;
                                                } else {
                                                    for (last_column--; last_column >=0; last_column--) {
                                                        (*result1) << gapCharacter;
                                                        (*result2) << string2->sData[last_column];
                                                    }
                                                    while (position>=0) {
                                                        (*result1) << str1->sData[position--];
                                                        (*result2) << gapCharacter;
                                                    }
                                                    break;
                                                }
                                            }

                                            if (current_column == last_column) { // insert in sequence 2
                                                (*result1) << str1->sData[position];
                                                (*result2) << gapCharacter;
                                            } else {
                                                last_column--;

                                                for (; last_column > current_column; last_column--) { // insert in column 1
                                                    (*result2) << string2->sData[last_column];
                                                    (*result1) << gapCharacter;
                                                }
                                                (*result1) << str1->sData[position];
                                                (*result2) << string2->sData[current_column];
                                            }
                                            //printf ("%s\n%s\n", result1->sData, result2->sData);
                                        }

                                        for (last_column--; last_column >=0; last_column--) {
                                            (*result1) << gapCharacter;
                                            (*result2) << string2->sData[last_column];
                                        }

                                        result1->Finalize();
                                        result1->Flip ();
                                        result2->Finalize();
                                        result2->Flip ();
                                        pairwiseComp->MStore ("1", new _FString(result1), false);
                                        pairwiseComp->MStore ("2", new _FString(result2), false);
                                    }
                                    /*
                                    long gap1c = 0,
                                         gap2c = 0;

                                     _Parameter scoreCheck = 0.;

                                     for (long sp = 0; sp<result1->sLength; sp++)
                                     {
                                         char cs1 = result1->sData[sp],
                                              cs2 = result2->sData[sp];

                                         if (cs1 == gapCharacter)
                                         {
                                             if (gap1c && doAffine)
                                                 scoreCheck -= gapExtend;
                                             else
                                                 scoreCheck -= gapOpen;
                                             gap2c = 0;
                                             gap1c++;
                                         }
                                         else
                                         if (cs2 == gapCharacter)
                                         {
                                             if (gap2c && doAffine)
                                                 scoreCheck -= gapExtend2;
                                             else
                                                 scoreCheck -= gapOpen2;
                                             gap1c = 0;
                                             gap2c++;
                                         }
                                         else
                                         {
                                             gap1c = 0;
                                             gap2c = 0;
                                             long code1 = ccount.lData[cs1],
                                                  code2 = ccount.lData[cs2];

                                             if (code1 >=0 && code2 >=0 )
                                                 scoreCheck += (*scoreMatrix)(code1,code2);
                                         }
                                     }
                                     if (doLocal)
                                     {
                                        for (long k = 0; result1->sData[k] == gapCharacter; k++)
                                            if (doAffine)
                                                scoreCheck += k?gapExtend:gapOpen;
                                            else
                                                scoreCheck += gapOpen;
                                         for (long k = 0; result2->sData[k] == gapCharacter; k++)
                                             if (doAffine)
                                                 scoreCheck += k?gapExtend2:gapOpen2;
                                             else
                                                 scoreCheck += gapOpen2;
                                         for (long k = result1->sLength-1; result1->sData[k] == gapCharacter; k--)
                                             if (doAffine)
                                                 scoreCheck += k==result1->sLength-1?gapOpen:gapExtend;
                                             else
                                                 scoreCheck += gapOpen;
                                         for (long k = result2->sLength-1; result2->sData[k] == gapCharacter; k--)
                                             if (doAffine)
                                                 scoreCheck += k==result2->sLength-1?gapOpen2:gapExtend2;
                                             else
                                                 scoreCheck += gapOpen2;
                                     }*/



                                    pairwiseComp->MStore ("0", new _Constant (score), false);
                                    /*pairwiseComp->MStore ("3", new _Constant (score2), false);
                                    pairwiseComp->MStore ("4", new _FString(result1), false);
                                    pairwiseComp->MStore ("5", new _FString(result2), false);
                                    pairwiseComp->MStore ("6", new _FString((_String*)ops.toStr()), false);
                                    pairwiseComp->MStore ("7", new _Constant (scoreCheck), false);*/
                                    alignedStrings->MStore (_String(s1), pairwiseComp, false);
                                }
                            }

                            storeResultIn->SetValue (alignedStrings, false);
                        }

                    } else {
                        errStr = seqAlignScore & " is a required option, which must be a square matrix with dimension matching the size of " & seqAlignMap;
                    }
                } else {
                    errStr = seqAlignMap & " is a required option, which must be a non-empty string without repeating characters ";
                }
            } else {
                errStr = *(_String*)parameters(2) & " was expected to be an associative array of alignment options";
            }
        } else {
            errStr = *(_String*)parameters(1) & " was expected to be a vector of strings";
        }
    }

    if (errStr.sLength) {
        errStr = errStr & " in call to " & blAlignSequences;
        WarnError (errStr);
    }
}

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

void      _ElementaryCommand::ExecuteCase26 (_ExecutionList& chain)
{
    chain.currentCommand++;
    // we have to build a list of _CalcNodes to deal with
    // all of the trees/nodes in ReplicateConstraint must be of the same topology
    // the constraint will be processed by trying all of the subnodes of the given node
    // and within each - trying all of the variables to see if the constraint is matched
    // exactly the same operation will be repeated on each of the other parameters

    _String *       replicateSource,
            thisS,
            prStr = GetStringFromFormula((_String*)parameters(0),chain.nameSpacePrefix);

    replicateSource = &prStr;

    _List           parts,
                    theConstraints;

    _SimpleList     thisIndex,
                    thisArgs;

    long            ind1    =   replicateSource->Find("this"),
                    ind2,
                    ind3,
                    ind4;

    if (ind1<0) {
        WarnError (*(_String*)parameters(0)&" has no 'this' references in call to ReplicateConstraint!");
        return ;
    }

    _SimpleList thisHits (parameters.lLength-1,0,0);

    while (ind1>=0) { // references to 'this' still exist
        ind2 = ind1+4; // look forward to the number of 'this'
        while ('0'<=replicateSource->sData[ind2] && replicateSource->sData[ind2]<='9') {
            ind2++;
        }

        ind3  = replicateSource->Cut(ind1+4,ind2-1).toNum();
        ind2  = replicateSource->FindEndOfIdent (ind1,-1,'?');
        // now ind1-ind2 contains a reference with this...
        _String newS  (*replicateSource,ind1,ind2);
        thisS = _String("this")&_String(ind3);
        if ((ind4 = ((_String*)parameters(ind3))->Find('.'))>=0) { // branch argument
            newS = newS.Replace (thisS,((_String*)parameters(ind3))->Cut(0,ind4-1), true);
        } else { // non-branch argument
            newS = newS.Replace (thisS,*((_String*)parameters(ind3)), true);
        }
        parts&& &newS;
        ind3--;
        thisIndex<<ind3; // sequence of references to this

        if (ind3<0 || ind3 >= thisHits.lLength) {
            WarnError (_String("Invalid reference to ") & thisS & " in the constraint specification");
            return ;
        }
        thisHits.lData[ind3] = 1;

        if (ind2>=replicateSource->sLength-1) {
            break;
        }
        ind1 = replicateSource->Find("this",ind2+1,-1);
        if (ind1==-1) {
            newS = replicateSource->Cut(ind2+1,-1);
        } else {
            newS = replicateSource->Cut(ind2+1,ind1-1);
        }
        parts&& &newS;
        thisIndex<<-1;
    }
    // now that the string is conveniently partritioned into blocks
    // we will check the arguments and store references

    for (ind1 = 1; ind1<parameters.lLength; ind1++) {
        if (thisHits.lData[ind1-1] == 0) {
            WarnError (_String("Unused ") & ind1 & "-th reference variable: " & *(_String*)parameters(ind1));
            return ;
        }

        ind2 = LocateVarByName (*(_String*)parameters(ind1));
        if (ind2<0) {
            _String newS = *(_String*)parameters(ind1) & " is undefined in call to ReplicateConstraint.";
            acknError (newS);
            return  ;
        }

        _Variable* thisNode = FetchVar (ind2);
        if (thisNode->ObjectClass()==TREE_NODE) {
            thisArgs<< (long)((_CalcNode*)thisNode)->LocateMeInTree();
        } else if (thisNode->ObjectClass()==TREE) {
            thisArgs<< (long)&((_TheTree*)thisNode)->GetRoot();
        } else {
            WarnError (*(_String*)parameters(ind1) & " is neither a tree nor a tree node in call to ReplicateConstraint.");
            return ;
        }
    }

    // now with this list ready we can recurse down the tree and produce the contsraints
    if (RecurseDownTheTree(thisArgs, parameters, theConstraints, parts, thisIndex)) {
        if (theConstraints.lLength) {
            ReportWarning  (_String("\nReplicateConstraint generated the following contsraints:"));
            _Parameter      doDeferSet;
            checkParameter (deferConstrainAssignment,doDeferSet,0.0);
            bool            applyNow = CheckEqual(doDeferSet,0.0);
            _String         *constraintAccumulator = (_String*)checkPointer(new _String(128L,true));

            if (applyNow) {
                deferSetFormula = new _SimpleList;
                checkPointer (deferSetFormula);
            }

            for (ind1 = 0; ind1 < theConstraints.lLength; ind1++) {
                replicateSource = (_String*)(theConstraints(ind1)->toStr());
                if (applyNow) {
                    _Formula rhs, lhs;
                    _FormulaParsingContext fpc (nil, chain.nameSpacePrefix);
                    ind2 = Parse (&rhs,*replicateSource,fpc,&lhs);
                    ExecuteFormula(&rhs,&lhs,ind2,fpc.assignmentRefID(),chain.nameSpacePrefix,fpc.assignmentRefType());
                }
                (*constraintAccumulator) << replicateSource;
                (*constraintAccumulator) << ';';
                (*constraintAccumulator) << '\n';
                //ReportWarning (*replicateSource);
                DeleteObject (replicateSource);
            }
            constraintAccumulator->Finalize();
            ReportWarning (*constraintAccumulator);
            CheckReceptacleAndStore (&lastSetOfConstraints,"ReplicateConstraint",false,new _FString(constraintAccumulator),false);
            if (applyNow) {
                FinishDeferredSF();
            }
        }
    }
}

//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase57 (_ExecutionList& chain)
{
    chain.currentCommand++;

    _String errStr;

    _Variable * storeResultIn = CheckReceptacle (&AppendContainerName(*(_String*)parameters(0),chain.nameSpacePrefix), blGetNeutralNull, true),
                *   sv            = FetchVar(LocateVarByName (AppendContainerName(*(_String*)parameters(2),chain.nameSpacePrefix))),
                    *  nsv           = FetchVar(LocateVarByName (AppendContainerName(*(_String*)parameters(3),chain.nameSpacePrefix)));

    _Parameter itCountV       = ProcessNumericArgument ((_String*)parameters(4),chain.nameSpacePrefix);

    _String   * lfName        = (_String*)parameters(1);

    long        f = FindLikeFuncName(AppendContainerName(*lfName,chain.nameSpacePrefix));

    if (f>=0) {
        if (sv && sv->ObjectClass () == MATRIX) {
            if (nsv && nsv->ObjectClass () == MATRIX) {
                _Matrix * sMatrix  = (_Matrix*)((_Matrix*)sv->Compute())->ComputeNumeric();
                _Matrix * nsMatrix = (_Matrix*)((_Matrix*)nsv->Compute())->ComputeNumeric();

                sMatrix->CheckIfSparseEnough (true);
                nsMatrix->CheckIfSparseEnough (true);

                if (   sMatrix->GetHDim()  == sMatrix->GetVDim() &&
                        nsMatrix->GetHDim() == nsMatrix->GetVDim() &&
                        sMatrix->GetHDim()  ==  nsMatrix->GetVDim() ) {
                    _LikelihoodFunction * theLF = (_LikelihoodFunction*)likeFuncList (f);

                    if (((_DataSetFilter*)dataSetFilterList (theLF->GetTheFilters() (0)))->GetDimension (true) == sMatrix->GetHDim()) {
                        long itCount = itCountV;
                        if (itCount>0) {
                            _AssociativeList * res = theLF->SimulateCodonNeutral ((_Matrix*)sMatrix, (_Matrix*)nsMatrix, itCount);
                            storeResultIn->SetValue (res,false);
                        } else {
                            errStr = "Invalid iterations per character state";
                        }
                    } else {
                        errStr = "Incompatible data and cost matrices";
                    }

                } else {
                    errStr = "Incompatible syn and non-syn cost matrix dimensions";
                }
            } else {
                errStr = "Invalid non-syn cost matrix argument";
            }
        } else {
            errStr = "Invalid syn cost matrix argument";
        }

    } else {
        errStr = _String("Likelihood function ") & *lfName & " has not been defined";
    }

    if (errStr.sLength) {
        errStr = errStr & " in call to " & blGetNeutralNull;
        WarnError (errStr);
    }
}

//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase58 (_ExecutionList& chain)
{
    chain.currentCommand++;

    _String     errStr;
    _String *   profileCode   = (_String*)parameters(0);

    if (*profileCode == _String ("START")) {
        if (chain.profileCounter) {
            DeleteObject (chain.profileCounter);
        }
        checkPointer(chain.profileCounter = new _Matrix (chain.lLength, 2, false, true));
        chain.doProfile = 1;
    } else if (*profileCode == _String ("PAUSE")) {
        chain.doProfile = 2;
    } else if (*profileCode == _String ("RESUME")) {
        chain.doProfile = 1;
    } else {
        _Variable * outVar = CheckReceptacle (&AppendContainerName(*profileCode,chain.nameSpacePrefix), blHBLProfile, true);
        if (outVar) {
            if (chain.profileCounter) {
                _AssociativeList * profileDump = new _AssociativeList;
                checkPointer     (profileDump);

                _SimpleList      instructions;
                _List            descriptions;

                for (long k=1; k<2*chain.lLength; k+=2) {
                    if (chain.profileCounter->theData[k] > 0.0) {
                        instructions << k/2;
                        _String * desc = (_String*)((_ElementaryCommand*)chain(k/2))->toStr();
                        descriptions << desc;
                        DeleteObject (desc);
                    }
                }

                _Matrix         * execProfile = new _Matrix (instructions.lLength,2,false,true),
                * instCounter = new _Matrix (instructions),
                * descList    = new _Matrix (descriptions);

                checkPointer    (execProfile);
                checkPointer    (instCounter);
                checkPointer    (descList);

                long k2 = 0;
                for (long m=1; m<2*chain.lLength; m+=2) {
                    if (chain.profileCounter->theData[m] > 0.0) {
                        execProfile->theData[k2++] = chain.profileCounter->theData[m];
                        execProfile->theData[k2++] = chain.profileCounter->theData[m-1];
                    }
                }

                _FString  aKey;
                *aKey.theString = "INSTRUCTION INDEX";
                profileDump->MStore (&aKey, instCounter, false);
                *aKey.theString = "INSTRUCTION";
                profileDump->MStore (&aKey, descList, false);
                *aKey.theString = "STATS";
                profileDump->MStore (&aKey, execProfile, false);
                outVar->SetValue (profileDump,false);
                chain.doProfile = 0;
                DeleteObject (chain.profileCounter);
                chain.profileCounter = nil;
            } else {
                errStr = "Profiler dump invoked before #profile START; ";
            }
        }
    }

}


//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase61 (_ExecutionList& chain)
{
    chain.currentCommand++;

    _PMathObj           avl1    = FetchObjectFromVariableByType (&AppendContainerName(*(_String*)parameters(1),chain.nameSpacePrefix),ASSOCIATIVE_LIST),
                        avl2  = FetchObjectFromVariableByType (&AppendContainerName(*(_String*)parameters(2),chain.nameSpacePrefix),ASSOCIATIVE_LIST),
                        start   = parameters.lLength>3?FetchObjectFromVariableByType (&AppendContainerName(*(_String*)parameters(3),chain.nameSpacePrefix),NUMBER):nil;

    if (! (avl1 && avl2)) {
        WarnError (_String ("Both arguments (") & *(_String*)parameters(1) & " and " & *(_String*)parameters(2) & " in a call to SCFG = ... must be evaluate to associative arrays");
    } else {
        Scfg    * scfg      = new Scfg ((_AssociativeList*)avl1,(_AssociativeList*)avl2,start?start->Value():0);
        _String scfgName    = AppendContainerName(*(_String*)parameters(0),chain.nameSpacePrefix);
        long    f           = FindSCFGName (scfgName);

        if (f==-1) {
            for (f=0; f<scfgNamesList.lLength; f++)
                if (((_String*)scfgNamesList(f))->sLength==0) {
                    break;
                }

            if (f==scfgNamesList.lLength) {
                scfgList << scfg;
                scfgNamesList&&(&scfgName);
                DeleteObject (scfg);
            } else {
                scfgNamesList.Replace(f,&scfgName,true);
                scfgList.lData[f] = (long)scfg;
            }
        } else {
            scfgNamesList.Replace(f,&scfgName,true);
            scfgList.Replace(f,scfg,false);
        }
    }
}

//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase63 (_ExecutionList& chain)
{
    chain.currentCommand++;

    /*_PMathObj         avl1    = FetchObjectFromVariableByType ((_String*)parameters(1),ASSOCIATIVE_LIST),
                        avl2    = FetchObjectFromVariableByType ((_String*)parameters(2),ASSOCIATIVE_LIST),
                        start   = parameters.lLength>3?FetchObjectFromVariableByType ((_String*)parameters(3),NUMBER):nil;

    if (! (avl1 && avl2))
    {
        _String errMsg = _String ("Both arguments (") & *(_String*)parameters(1) & " and " & *(_String*)parameters(2) & " in a call to SCFG = ... must be evaluate to associative arrays";
        WarnError (errMsg);
    }
    else
    {
        Scfg    * scfg   = new Scfg ((_AssociativeList*)avl1,(_AssociativeList*)avl2,start?start->Value():0);
        _String * str    = (_String*)parameters(0);
        long    f        = FindSCFGName (*str);

        if (f==-1)
        {
            for (f=0; f<scfgNamesList.lLength; f++)
                if (((_String*)scfgNamesList(f))->sLength==0)
                    break;

            if (f==scfgNamesList.lLength)
            {
                scfgList << scfg;
                scfgNamesList&&(str);
                DeleteObject (scfg);
            }
            else
            {
                scfgNamesList.Replace(f,str,true);
                scfgList.lData[f] = (long)scfg;
            }
        }
        else
        {
            scfgNamesList.Replace(f,str,true);
            scfgList.Replace(f,scfg,false);
        }
    }   */
}


//____________________________________________________________________________________

void    _ElementaryCommand::ExecuteCase64 (_ExecutionList& chain)
{
	ReportWarning (_String("ExecuteCase64()"));
    chain.currentCommand++;

    _PMathObj   avl1    = FetchObjectFromVariableByType (&AppendContainerName(*(_String*)parameters(1),chain.nameSpacePrefix), ASSOCIATIVE_LIST);

    if (! (avl1)) {
        WarnError (_String ("Argument (") & *(_String*)parameters(1) & " in call to BGM = ... must evaluate to associative array");
    } else {
        _BayesianGraphicalModel * bgm   = new _BayesianGraphicalModel ((_AssociativeList *) avl1);

        _String bgmName     = AppendContainerName (*(_String *) parameters(0), chain.nameSpacePrefix);
        long    bgmIndex    = FindBgmName (bgmName);
		
        if (bgmIndex == -1) {   // not found
            for (bgmIndex = 0; bgmIndex < bgmNamesList.lLength; bgmIndex++) {
                // locate empty strings in name list
                if (((_String *)bgmNamesList(bgmIndex))->sLength == 0) {
                    break;
                }
            }

            if (bgmIndex == bgmNamesList.lLength) {
                // reached end of list without finding empty string, append new string
                bgmList.AppendNewInstance(bgm);
                bgmNamesList && (&bgmName);
            } else {
                // replace empty string in list
                bgmNamesList.Replace (bgmIndex, &bgmName, true);
                bgmList.Replace (bgmIndex, bgm, false);
            }
        } else { // 20070626: SLKP edit to deal with already existing BGMs
            bgmNamesList.Replace(bgmIndex,&bgmName,true);
            bgmList.Replace(bgmIndex,bgm,false);
        }
		
		ReportWarning(_String("Created BGM ") & bgmName & " at index " & bgmIndex);
    }
}


//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructSCFG (_String&source, _ExecutionList&target)
// syntax: SCFG ident = (Rules1, Rules2 <,start>)
{

    long    mark1 = source.FirstSpaceIndex(0,-1,1),
            mark2 = source.Find ('=', mark1, -1);

    _String scfgID (source, mark1+1,mark2-1);

    if (mark1==-1 || mark2==-1 || mark1+1>mark2-1 || !scfgID.IsValidIdentifier(true)) {
        WarnError ("SCFG declaration missing a valid identifier");
        return false;
    }

    _List pieces;

    mark1 = source.Find ('(',mark2,-1);
    if (mark1 >= 0) {
        ExtractConditions (source,mark1+1,pieces,',');
    }

    if (pieces.lLength != 2 && pieces.lLength != 3) {
        WarnError ("Expected: SCFG ident = (Rules1, Rules2 <,start>)");
        return false;
    }

    _ElementaryCommand * scfg = new _ElementaryCommand (61);

    scfg->parameters    &&(&scfgID);
    scfg->addAndClean(target,&pieces,0);
    return true;
}

//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructBGM (_String&source, _ExecutionList&target)
// syntax: BGM ident = (<nodes>)
{
	ReportWarning(_String("ConstructBGM()"));
    // locate ident in HBL string
    long    mark1 = source.FirstSpaceIndex(0,-1,1),
            mark2 = source.Find ('=', mark1, -1);

    // assign ident to _String variable
    _String bgmID (source, mark1+1,mark2-1);

    if (mark1==-1 || mark2==-1 || mark1+1>mark2-1 || !bgmID.IsValidIdentifier(true)) {
        WarnError ("BGM declaration missing a valid identifier");
        return false;
    }

    // extract arguments from remainder of HBL string
    _List pieces;

    mark1 = source.Find ('(',mark2,-1);
    if (mark1 >= 0) {
        ExtractConditions (source,mark1+1,pieces,',');
    }

    if (pieces.lLength != 1) {
        WarnError ("Expected: BGM ident = (<nodes>)");
        return false;
    }

    _ElementaryCommand * bgm = new _ElementaryCommand (64);
    bgm->parameters && (&bgmID);
    bgm->addAndClean(target,&pieces,0);
	
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


