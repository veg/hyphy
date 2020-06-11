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

#include      <ctype.h>

#include      "alignment.h"
#include      "hy_string_buffer.h"
#include      "global_things.h"
#include      "trie.h"
#include      "likefunc.h"
#include      "scfg.h"
#include      "global_object_lists.h"

#include      "bayesgraph.h"



using namespace hyphy_global_objects;
using namespace hy_global;

//____________________________________________________________________________________
// global variables

_String     isDynamicGraph          ("BGM_DYNAMIC"),
            treeNodeNameMapping     ("TREE_NODE_NAME_MAPPING");


extern      _String                 blDoSQL,
            blAlignSequences,
            blGetNeutralNull,
            blHBLProfile,
            blDeleteObject,
            blGetString,
            blRequireVersion,
            blAssert;

_SimpleList _HY_HBLCommandHelperAux;
            
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

//____________________________________________________________________________________

_HBLCommandExtras* _hyInitCommandExtras (const long cut, const long conditions, _String const& commandInvocation, const char sep, const bool doTrim, const bool isAssignment, const bool needsVerb, _SimpleList * conditionList) {
    
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
    commandInfo->command_invocation                  < new _String (commandInvocation);

    return                                             commandInfo;
    
}

//____________________________________________________________________________________

bool      _ElementaryCommand::ExtractValidateAddHBLCommand (_String& current_stream,const long command_code,  _List* pieces, _HBLCommandExtras* command_spec, _ExecutionList& command_list)

{
    if (command_spec->is_assignment) {
        // TBA
    } else {
        // by default push all of the 'pieces' arguments to the "argument" list
        (new _ElementaryCommand (command_code))->addAndClean  (command_list, pieces, 0);
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
    _HY_ValidHBLExpressions.Insert ("ChoiceList(",                          HY_HBL_COMMAND_CHOICE_LIST);
    _HY_ValidHBLExpressions.Insert ("GetInformation(",                      HY_HBL_COMMAND_GET_INFORMATION);
    _HY_ValidHBLExpressions.Insert ("ExecuteCommands(",                     HY_HBL_COMMAND_EXECUTE_COMMANDS);
    _HY_ValidHBLExpressions.Insert ("ExecuteAFile(",                        HY_HBL_COMMAND_EXECUTE_A_FILE);
    _HY_ValidHBLExpressions.Insert ("LoadFunctionLibrary(",					HY_HBL_COMMAND_LOAD_FUNCTION_LIBRARY);
    _HY_ValidHBLExpressions.Insert ("FindRoot(",                            HY_HBL_COMMAND_FIND_ROOT);
    _HY_ValidHBLExpressions.Insert ("MPIReceive(",                          HY_HBL_COMMAND_MPI_RECEIVE);
    _HY_ValidHBLExpressions.Insert ("MPISend(",                             HY_HBL_COMMAND_MPI_SEND);
    _HY_ValidHBLExpressions.Insert ("GetDataInfo(",                         HY_HBL_COMMAND_GET_DATA_INFO);
    _HY_ValidHBLExpressions.Insert ("StateCounter(",                        HY_HBL_COMMAND_STATE_COUNTER);
    _HY_ValidHBLExpressions.Insert ("Integrate(",                           HY_HBL_COMMAND_INTEGRATE);
    _HY_ValidHBLExpressions.Insert ("DoSQL(",                               HY_HBL_COMMAND_DO_SQL);
    _HY_ValidHBLExpressions.Insert ("Topology ",                            HY_HBL_COMMAND_TOPOLOGY);
    _HY_ValidHBLExpressions.Insert ("AlignSequences(",                      HY_HBL_COMMAND_ALIGN_SEQUENCES);
    _HY_ValidHBLExpressions.Insert ("#profile",                             HY_HBL_COMMAND_PROFILE);
    _HY_ValidHBLExpressions.Insert ("SCFG ",                                HY_HBL_COMMAND_SCFG);
    _HY_ValidHBLExpressions.Insert ("NeuralNet ",                           HY_HBL_COMMAND_NEURAL_NET);
    _HY_ValidHBLExpressions.Insert ("BGM ",                                 HY_HBL_COMMAND_BGM);
    _HY_ValidHBLExpressions.Insert ("SimulateDataSet",                      HY_HBL_COMMAND_SIMULATE_DATA_SET);
    _HY_ValidHBLExpressions.Insert ("KeywordArgument",                      HY_HBL_COMMAND_KEYWORD_ARGUMENT);
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

    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_REPLICATE_CONSTRAINT,
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("ReplicateConstraint(", HY_HBL_COMMAND_REPLICATE_CONSTRAINT,false),
                                                                -2,
                                                                "ReplicateConstraint(<constraint pattern in terms of 'this1', 'this2',...>, <an argument to replace 'this*', for each 'thisN' in the pattern);",','));


    lengthOptions.Clear();lengthOptions.Populate (3,1,1);
    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_EXECUTE_COMMANDS,
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("ExecuteCommands(", HY_HBL_COMMAND_EXECUTE_COMMANDS,false),
                                                                -1,
                                                                "ExecuteCommands(<source code>, [optional <'compiled' | (input redirect , [optional <namespace>]) ])",
                                                                ',',
                                                                true,
                                                                false,
                                                                false,
                                                                &lengthOptions));

    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_EXECUTE_A_FILE,
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("ExecuteAFile(", HY_HBL_COMMAND_EXECUTE_A_FILE,false),
                                                                -1,
                                                                "ExecuteAFile(<file path>, [optional <'compiled' | (input redirect , [optional <namespace>]) ])",
                                                                ',',
                                                                true,
                                                                false,
                                                                false,
                                                                &lengthOptions));

    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_LOAD_FUNCTION_LIBRARY,
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("LoadFunctionLibrary(", HY_HBL_COMMAND_EXECUTE_A_FILE,false),
                                                                -1,
                                                                "LoadFunctionLibrary(<file path | library name>, [optional <'compiled' | (input redirect , [optional <namespace>]) ])",
                                                                ',',
                                                                true,
                                                                false,
                                                                false,
                                                                &lengthOptions));

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
    

    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_GET_INFORMATION,
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("GetInformation(", HY_HBL_COMMAND_GET_INFORMATION,false),
                                                                2,
                                                                "GetInformation(<receptacle>, <DataSet or DataSetFilter or LikelihoodFunction or Model or Variable or Regexp or String>",
                                                                ','));

    lengthOptions.Clear();lengthOptions.Populate (4,2,1); // 2, 3, 4, 5
    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_GET_DATA_INFO,
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("GetDataInfo(", HY_HBL_COMMAND_GET_DATA_INFO,false),
                                                                -1,
                                                                "GetDataInfo(<receptacle>, <DataSet or DataSetFilter>, [optional <sequence ref, site ref | sequence 1 , sequence 2, DISTANCES>])",
                                                                ',',
                                                                true,
                                                                false,
                                                                false,
                                                                &lengthOptions));

    lengthOptions.Clear();lengthOptions.Populate (2,2,1); // 2, 3
    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_MPI_SEND,
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("MPISend(", HY_HBL_COMMAND_MPI_SEND,false),
                                                                -1,
                                                                "MPISend(<node id>, <string | likelihood function ID | filename [in conjuction with argument 3]>, [if specified, treat the second argument as a script path, and use the dict supplied here as input options to the script])",
                                                                ',',
                                                                true,
                                                                false,
                                                                false,
                                                                &lengthOptions));

  _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_CONSTRUCT_CATEGORY_MATRIX,
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("ConstructCategoryMatrix(", HY_HBL_COMMAND_CONSTRUCT_CATEGORY_MATRIX,false),
                                                                -1,
                                                                "ConstructCategoryMatrix(<receptacle>, <Likelihood Function|Tree>, [optional <COMPLETE|SHORT|WEIGHTS|CLASSES (default = COMPLETE)> , matrix argument with partitions to include (defaut = all)>])",
                                                                ',',
                                                                true,
                                                                false,
                                                                false,
                                                                &lengthOptions));

    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_OPTIMIZE,
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("Optimize(", HY_HBL_COMMAND_OPTIMIZE,false),
                                                                -1,
                                                                "Optimize (<receptacle>, <likelihood function/scfg/bgm>, [optional dictionary of arguments]",',',
                                                                true,
                                                                false,
                                                                false,
                                                                &lengthOptions));

    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_FIND_ROOT,
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("FindRoot(", HY_HBL_COMMAND_FIND_ROOT,false),
                                                                5,
                                                                "FindRoot (<receptacle>, <expression>, <variable to solve for>,<left bound>,<right bound>)",','));


    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_INTEGRATE,
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("Integrate(", HY_HBL_COMMAND_INTEGRATE,false),
                                                                5,
                                                                "Integrate (<receptacle>, <expression>, <variable to integrate over for>,<left bound>,<right bound>)",','));


    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_MPI_RECEIVE,
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("MPIReceive(", HY_HBL_COMMAND_MPI_RECEIVE,false),
                                                                3,
                                                                "MPIReceive (<from node; or -1 to receive from any>, <message storage>, <sender index storage>)",','));

    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_LFCOMPUTE,
                                      (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("LFCompute(", HY_HBL_COMMAND_LFCOMPUTE,false),
                                                                  2, 
                                                                  "LFCompute (<likelihood function/scfg/bgm>,<LF_START_COMPUTE|LF_DONE_COMPUTE|receptacle>)",','));


    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_COVARIANCE_MATRIX, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("CovarianceMatrix(", HY_HBL_COMMAND_COVARIANCE_MATRIX,false),
                                                                2, 
                                                                "CovarianceMatrix (<receptacle>, <likelihood function/scfg/bgm>)",','));

    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_DO_SQL,
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("DoSQL(", HY_HBL_COMMAND_DO_SQL,false),
                                                                3,
                                                                "DoSQL (<dbID | SQL_OPEN | SQL_CLOSE>, <transaction string | file name>, <ID here | result here>)",','));

  
    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_SELECT_TEMPLATE_MODEL,
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("SelectTemplateModel(", HY_HBL_COMMAND_SELECT_TEMPLATE_MODEL,false),
                                    1, 
                                    "SelectTemplateModel(<DataSetFilter>);"));


    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_USE_MODEL, 
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("UseModel(", HY_HBL_COMMAND_USE_MODEL,false),
                                                                1, 
                                                                "UseModel (<model ID>)",','));


    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_ALIGN_SEQUENCES,
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("AlignSequences(", HY_HBL_COMMAND_ALIGN_SEQUENCES,false),
                                                                3,
                                                                "AlignSequences (result, sequences, options)",','));

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

    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_FSCANF,
                                  (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("fscanf(", HY_HBL_COMMAND_FSCANF,false),
                                                              -3,
                                                              _String ("fscanf(file path (string),<optional 'REWIND'>,'type 1 (") & _String ((_String*)_ElementaryCommand::fscanf_allowed_formats.Join('|')) & ")[optional , <type 2>, ... <type N>]', var1 [optional , <var 2>, ... <var N>])",','));

    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_SSCANF,
                                  (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("sscanf(", HY_HBL_COMMAND_SSCANF,false),
                                                              -3,
                                                              _String ("sscanf(string,<optional 'REWIND'>,'type 1 (") & _String ((_String*)_ElementaryCommand::fscanf_allowed_formats.Join('|')) & ")[optional , <type 2>, ... <type N>]', var1 [optional , <var 2>, ... <var N>])",','));
  

    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_CHOICE_LIST,
                                  (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("ChoiceList(", HY_HBL_COMMAND_CHOICE_LIST,false),
                                                              -5,
                                                              _String ("ChoiceList(store_here {ID}, title {String}, how many choices (0 for any number >= 1) {Integer}, [NO_SKIP or SKIP_NONE | list of indices not to show as options], [source object | comma separated list of 'key', 'description' pairs]"),
                                                              ','));

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

    lengthOptions.Clear();lengthOptions.Populate (3,2,1);
    _HY_HBLCommandHelper.Insert    ((BaseRef)HY_HBL_COMMAND_KEYWORD_ARGUMENT,
                                    (long)_hyInitCommandExtras (_HY_ValidHBLExpressions.Insert ("KeywordArgument(", HY_HBL_COMMAND_KEYWORD_ARGUMENT,false),
                                                                -1,
                                                                "KeywordArgument (keyword, description, [default value, [dialog reference]])",
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


    _HY_MatrixRandomValidPDFs.Insert ("Dirichlet", _HY_MATRIX_RANDOM_DIRICHLET,
                                      "Gaussian", _HY_MATRIX_RANDOM_GAUSSIAN,
                                      "Wishart", _HY_MATRIX_RANDOM_WISHART,
                                      "InverseWishart", _HY_MATRIX_RANDOM_INVERSE_WISHART,
                                      "Multinomial", _HY_MATRIX_RANDOM_MULTINOMIAL);
  
  
  _List keywords;
  keywords << &blReturn << &blDataSet << &blDataSetFilter << &blTree << &blTopology << &blLF << &blLF3 << &blSCFG;
  
  
  for (long key = 0; key < keywords.lLength; key++) {
    _String* key_string = (_String*)keywords.GetItem (key);
    _HY_HBL_KeywordsPreserveSpaces.Insert (key_string->Reverse().Cut (1,kStringEnd), key_string->length()-1);
  }
  
}

//____________________________________________________________________________________
void         InsertVarIDsInList     (_AssociativeList* theList , _String const& theKey, _SimpleList const& varIDs) {
    _FString arrayKey (theKey, false);
    _Matrix *mxEntry = nil;

    if (varIDs.lLength) {
        _List     varNames;
        for (unsigned long i=0; i < varIDs.lLength; i++) {
            _Variable* v = LocateVar (varIDs.list_data[i]);
            if (v) {
                varNames << v->GetName();
            }
        }
        mxEntry = new _Matrix (varNames);
    } else {
        mxEntry = new _Matrix;
    }

    theList->MStore (&arrayKey,mxEntry,false);
}

//____________________________________________________________________________________
void         InsertStringListIntoAVL    (_AssociativeList* theList , _String const& theKey, _SimpleList const& stringsToPick, _List const& theStrings) {
    _FString arrayKey (theKey, false);
    _Matrix *mxEntry = nil;

    if (stringsToPick.lLength) {
        _List     theNames;
        for (unsigned long i=0; i < stringsToPick.lLength; i++) {
            _String * v = (_String*)theStrings.GetItem(stringsToPick.list_data[i]);
            if (v) {
                theNames << v;
            }
        }
        mxEntry = new _Matrix (theNames);
    } else {
        mxEntry = new _Matrix;
    }

    theList->MStore (&arrayKey,mxEntry,false);
}

//____________________________________________________________________________________

bool    _ElementaryCommand::ConstructProfileStatement (_String&source, _ExecutionList&target)
// syntax: #profile START|PAUSE|RESUME|indetifier to dump in
{

    _List pieces;
    ExtractConditions (source,blHBLProfile.length()+1,pieces,';');
    if (pieces.lLength!=2) {
        HandleApplicationError (_String ("Expected syntax:")& blHBLProfile &" START|PAUSE|RESUME|where to store)");
        return false;
    }

    _ElementaryCommand *sp = new _ElementaryCommand (58);
    sp->addAndClean(target,&pieces,0);
    return true;
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

void      _ElementaryCommand::ExecuteDataFilterCases (_ExecutionList& chain) {
  
    using namespace hyphy_global_objects;
  
    chain.currentCommand++;

    _String dataObjectID = chain.AddNameSpaceToID(*(_String*)parameters(1));
  
    /*printf ("%s raw : '%s' processed : '%s'\n", __PRETTY_FUNCTION__,  (const char*)*(_String*)parameters(1), (const char*) dataObjectID);
    */
  
    long dsID           = (parameters.lLength>2)?FindDataSetName (dataObjectID):-1;
    bool isFilter       = false;

    if (dsID == -1) {
        dsID = (parameters.lLength>2)? FindDataFilter (dataObjectID):-1;
        if (dsID == -1) {
            _AssociativeList * numericFilter = (_AssociativeList*)FetchObjectFromVariableByType(&dataObjectID, ASSOCIATIVE_LIST);
            if (numericFilter) {
                _String errCode;

                long    categoryCount = 1;

                if (parameters.lLength > 2) {
                    // have multiple categories
                    categoryCount = (long) ProcessNumericArgument((_String*)parameters(2),nil);
                }

                _String const namesKey ("FILTER_NAMES"),
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
                                _String             dataFilterID (chain.AddNameSpaceToID(*(_String*)parameters(0))),
                                                    dataSetName  = dataFilterID & "_internal_ds";
 
                                _DataSet * dummyDS = new _DataSet;
                              
                                dummyDS->SetNoSpecies (seqNames.lLength);
                                dummyDS->SetNames     (seqNames);
                                dummyDS->GetTheMap().Populate (sitePatterns,0,1);
 
                                AddDataSetToList  (dataSetName, dummyDS);

                                _DataSetFilterNumeric * dsn = new _DataSetFilterNumeric (freqList,goodSeqs,dummyDS,categoryCount);
                              
                                if (StoreDataFilter (dataFilterID, dsn) >= 0) {
                                  return;
                                }
                                
                                errCode = _String ("Failed to store _DataSetFilterNumeric '") & dataFilterID & "'";
                              
                            } else {
                                errCode = _String ("Site frequency patterns/numeric vectors did not pass dimension checks in call to CreateFilter");
                            }
                        }
                    }

                }
                if (errCode.nonempty()) {
                    HandleApplicationError(errCode);
                    return;
                }

            }
            HandleApplicationError (((_String)("DataSet(Filter)/Associative Array ")&dataObjectID&_String(" has not been properly initialized")));
            return;
        }
        isFilter = true;
    }

    // build the formula from the 2nd parameter (unit size)

    unsigned char                unit = ProcessNumericArgument((_String*)parameters(2),chain.nameSpacePrefix);
    // here's our unit

    _String             dataFilterID (chain.AddNameSpaceToID(*(_String*)parameters(0))),
                        site_partition,
                        sequence_partition;



    if (parameters.countitems()>3) {
        sequence_partition = *(_String*)parameters(3);
    }
    
    if (parameters.countitems()>4) {
        site_partition = *(_String*)parameters(4);
    }

    _DataSet            *dataset;

    _SimpleList         site_list,
                        sequence_list;

    site_list.RequestSpace (1024UL);
    sequence_list.RequestSpace (1024UL);
    
    auto ensure_site_partition = [] (_String& site_part, const long code, const long max_site) -> void {
        if (code != 6 && site_part.empty () ) {
            site_part = _String ("\"0-") & _String (max_site - 1L) & ("\"");
        }
    };

    if (!isFilter) {
        dataset = (_DataSet*)dataSetList(dsID);
        dataset -> ProcessPartition (site_partition,site_list,false, unit, nil, nil, chain.GetNameSpace());
        ensure_site_partition (sequence_partition, code, dataset->NoOfColumns());
        dataset->ProcessPartition (sequence_partition,sequence_list,true, unit, nil, nil, chain.GetNameSpace());

    } else {
        const _DataSetFilter * dataset1 = GetDataFilter (dsID);
        dataset1->GetData()->ProcessPartition (site_partition,site_list,false, unit, &dataset1->theNodeMap, &dataset1->theOriginalOrder, chain.GetNameSpace());
        ensure_site_partition (sequence_partition, code, dataset1->GetSiteCount());
        dataset1->GetData()->ProcessPartition (sequence_partition,sequence_list ,true,  unit, &dataset1->theOriginalOrder, &dataset1->theNodeMap, chain.GetNameSpace());
        dataset = (_DataSet*)dataset1;
    }

    if (code!=6) {
        if (sequence_list.countitems() % unit) {
            ReportWarning ((_String)"Unit size of "& _String((long)unit) & " doesn't divide the length of specified partition. The partition has been trimmed at the end.");
            for (long chop = sequence_list.countitems(); chop>0; chop--) {
                sequence_list.Pop();
            }
        }
        if (code == 27) {
            sequence_list.Permute (unit);
        } else {
            sequence_list.PermuteWithReplacement(unit);
        }

    }

    _DataSetFilter * thedf = new _DataSetFilter;
    thedf->SetFilter (dataset, unit, site_list, sequence_list, isFilter);

    if (parameters.lLength>5) {
        thedf->SetExclusions (GetStringFromFormula((_String*)parameters(5),chain.nameSpacePrefix));
    } else if ( code!=6 && isFilter ) {
        const _DataSetFilter * df1 = GetDataFilter (dsID);
        if (df1->theExclusions.nonempty()) {
            thedf->theExclusions << df1->theExclusions;
            thedf->SetDimensions();
        }
    }

    thedf->SetDimensions();
    thedf->SetupConversion();
  
    StoreDataFilter(dataFilterID, thedf, true);
}



//____________________________________________________________________________________

void      _ElementaryCommand::ExecuteCase54 (_ExecutionList& chain) {
    chain.currentCommand++;

    SetStatusLine (_String("Constructing Topology ")&*(_String*)parameters(0));

    _String  *treeSpec = ((_String*)parameters(1));
    _TreeTopology * tr = nil;
  
    _AssociativeList* mapping = (_AssociativeList*)FetchObjectFromVariableByType(&treeNodeNameMapping, ASSOCIATIVE_LIST);

    if (treeSpec->nonempty()) {
        if (treeSpec->char_at (0)!='(') {
            _Variable* testTree = FetchVar(LocateVarByName (AppendContainerName(*treeSpec,chain.nameSpacePrefix)));
            if (testTree && testTree->ObjectClass () == TREE) {
                tr = new _TreeTopology ((_TheTree*)testTree);
            } else {
                _String   flaData (*treeSpec);
                _Formula  nameForm (flaData,chain.nameSpacePrefix);
                HBLObjectRef formRes = nameForm.Compute();
                if (formRes&&formRes->ObjectClass () == STRING)
                    tr = new _TreeTopology (AppendContainerName(*(_String*)parameters(0),chain.nameSpacePrefix),
                                            ((_FString*)formRes)->get_str(),
                                            false, mapping);
            }
        } else
            tr = new _TreeTopology (AppendContainerName(*(_String*)parameters(0),chain.nameSpacePrefix),
                                    *(_String*)parameters(1),false);
    }

    if (!tr) {
        HandleApplicationError ("Illegal right hand side in call to Topology id = ...; it must be a string, a Newick tree spec or a topology");
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
        chain.profileCounter = new _Matrix (chain.lLength, 2, false, true);
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

                _SimpleList      instructions;
                _List            descriptions;

                for (unsigned long k=1UL; k<2*chain.lLength; k+=2UL) {
                    if (chain.profileCounter->theData[k] > 0.0) {
                        instructions << k/2;
                        descriptions.AppendNewInstance((_String*)((_ElementaryCommand*)chain(k/2))->toStr());
                    }
                }

                _Matrix         * execProfile = new _Matrix (instructions.lLength,2,false,true),
                * instCounter = new _Matrix (instructions),
                * descList    = new _Matrix (descriptions);

                unsigned long k2 = 0UL;
                for (unsigned long m=1UL; m<2*chain.lLength; m+=2UL) {
                    if (chain.profileCounter->theData[m] > 0.0) {
                        execProfile->theData[k2++] = chain.profileCounter->theData[m];
                        execProfile->theData[k2++] = chain.profileCounter->theData[m-1];
                    }
                }

                profileDump->MStore ("INSTRUCTION INDEX", instCounter, false);
                profileDump->MStore ("INSTRUCTION", descList, false);
                profileDump->MStore ("STATS", execProfile, false);
                outVar->SetValue (profileDump,false,true,NULL);
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

    HBLObjectRef           avl1    = FetchObjectFromVariableByType (&AppendContainerName(*(_String*)parameters(1),chain.nameSpacePrefix),ASSOCIATIVE_LIST),
                        avl2  = FetchObjectFromVariableByType (&AppendContainerName(*(_String*)parameters(2),chain.nameSpacePrefix),ASSOCIATIVE_LIST),
                        start   = parameters.lLength>3?FetchObjectFromVariableByType (&AppendContainerName(*(_String*)parameters(3),chain.nameSpacePrefix),NUMBER):nil;

    if (! (avl1 && avl2)) {
        HandleApplicationError (_String ("Both arguments (") & *(_String*)parameters(1) & " and " & *(_String*)parameters(2) & " in a call to SCFG = ... must be evaluate to associative arrays");
    } else {
        Scfg    * scfg      = new Scfg ((_AssociativeList*)avl1,(_AssociativeList*)avl2,start?start->Value():0);
        _String scfgName    = AppendContainerName(*(_String*)parameters(0),chain.nameSpacePrefix);
        long    f           = FindSCFGName (scfgName);

        if (f==-1) {
            for (f=0; f<scfgNamesList.lLength; f++)
                if (((_String*)scfgNamesList(f))->empty()) {
                    break;
                }

            if (f==scfgNamesList.lLength) {
                scfgList << scfg;
                scfgNamesList&&(&scfgName);
                DeleteObject (scfg);
            } else {
                scfgNamesList.Replace(f,&scfgName,true);
                scfgList.list_data[f] = (long)scfg;
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
                scfgList.list_data[f] = (long)scfg;
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

    HBLObjectRef   avl1    = FetchObjectFromVariableByType (&AppendContainerName(*(_String*)parameters(1),chain.nameSpacePrefix), ASSOCIATIVE_LIST);

    if (! (avl1)) {
        HandleApplicationError (_String ("Argument (") & *(_String*)parameters(1) & " in call to BGM = ... must evaluate to associative array");
    } else {
        _BayesianGraphicalModel * bgm   = new _BayesianGraphicalModel ((_AssociativeList *) avl1);

        _String bgmName     = AppendContainerName (*(_String *) parameters(0), chain.nameSpacePrefix);
        long    bgmIndex    = FindBgmName (bgmName);
		
        if (bgmIndex == -1) {   // not found
            for (bgmIndex = 0; bgmIndex < bgmNamesList.lLength; bgmIndex++) {
                // locate empty strings in name list
                if (((_String *)bgmNamesList(bgmIndex))->empty() == 0) {
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


void   _ElementaryCommand::appendCompiledFormulae(_Formula* f, _Formula *f2) {
  if (f || f2) {
    _SimpleList*        varList = new _SimpleList;
    _AVLList            varListA (varList);
    if (f)
      f->ScanFForVariables (varListA, true, true, true, true);
    if (f2)
      f2->ScanFForVariables(varListA, true, true);
    varListA.ReorderList();
    listOfCompiledFormulae<<(long)this;
    compiledFormulaeParameters.AppendNewInstance(varList);
  }
}

//____________________________________________________________________________________
bool    _ElementaryCommand::DecompileFormulae (void) {
  switch (code) {
    case 0:
      if (simpleParameters.nonempty()) {
        _Formula* f = (_Formula*)simpleParameters.list_data[1],
                *f2 = (_Formula*)simpleParameters.list_data[2] ;
        if (f) {
          delete f;
        }
        if (f2) {
          delete f2;
        }
        simpleParameters.Clear();
        
        return true;
      }
      break;
    case 4: {
      if (parameters.lLength && simpleParameters.lLength == 3) {
        _Formula* f = (_Formula*)simpleParameters.list_data[2];
        if (f) {
          delete f;
        }
        simpleParameters.Delete (2);
        return true;
      }
      break;
    }
    case 14: {
      if (parameters.lLength && simpleParameters.lLength == 2) {
        _Formula* f = (_Formula*)simpleParameters.list_data[1];
        if (f) {
          delete f;
        }
        simpleParameters.Delete (1);
        return true;
      }
      break;
    }
      

  }
  return false;
}

//____________________________________________________________________________________
bool    _ElementaryCommand::ConstructSCFG (_String&source, _ExecutionList&target) {
// syntax: SCFG ident = (Rules1, Rules2 <,start>)

    long    mark1 = source.FirstNonSpaceFollowingSpace (),
            mark2 = mark1 > 0 ? source.FindTerminator (mark1 + 1, "=") : 0;

 
    if ( mark1==-1 || mark2==-1 || mark1+1 > mark2  ) {
        HandleApplicationError ("SCFG declaration missing a valid identifier");
        return false;
    }

    _String scfgID (source,mark1,mark2-1);
 
    _List pieces;
    mark2 ++;
    mark1 = source.ExtractEnclosedExpression(mark2, '(', ')', fExtractRespectQuote | fExtractRespectEscape);

    ExtractConditions (source,mark2+1,pieces,',');
 
    if ( mark1==-1 || mark2==-1 || mark1<mark2 || (pieces.lLength != 2 && pieces.lLength != 3)) {
        HandleApplicationError ("Expected: SCFG ident = (Rules1, Rules2 <,start>)");
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
    //ReportWarning(_String("ConstructBGM()"));
    // locate ident in HBL string
    long    mark1 = source.FirstNonSpaceFollowingSpace (),
            mark2 = mark1 > 0 ? source.FindTerminator (mark1 + 1, "=") : 0;

    // assign ident to _String variable

    if ( mark1==-1 || mark2==-1 || mark1+1 > mark2  ) {
        HandleApplicationError ("BGM declaration missing a valid identifier");
        return false;
    }

    _String bgmID (source,mark1,mark2-1);

    _List pieces;
    mark2 ++;
    mark1 = source.ExtractEnclosedExpression(mark2, '(', ')', fExtractRespectQuote | fExtractRespectEscape);
    
    ExtractConditions (source,mark2+1,pieces,',');

    if ( mark1==-1 || mark2==-1 || mark1<mark2 || (pieces.lLength != 1)) {
        HandleApplicationError ("Expected: BGM ident = (<nodes>)");
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
        if (modelTypeList.list_data[mid] == 0) {
            mm = (_Matrix*)FetchObjectFromVariableByTypeIndex(modelMatrixIndices.list_data[mid],MATRIX);
        } else {
            mm = nil;
        }

        long fvi = modelFrequenciesIndices.list_data[mid];
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
    if (mid >= 0 && modelTypeList.list_data[mid] == 0) {
        mm = LocateVar(modelMatrixIndices.list_data[mid]);
    } else {
        mm = nil;
    }

    long fvi = modelFrequenciesIndices.list_data[mid];
    fv = LocateVar (fvi>=0?fvi:(-fvi-1));
    mbf = (fvi>=0);
}



//____________________________________________________________________________________

void    ScanModelForVariables        (long modelID, _AVLList& theReceptacle, bool inclG, long modelID2, bool inclCat)
{
    if (modelID != HY_NO_MODEL) {
        // standard rate matrix
        if (modelTypeList.list_data[modelID] == 0) {
            ((_Matrix*) (LocateVar(modelMatrixIndices.list_data[modelID])->GetValue()))->ScanForVariables2(theReceptacle,inclG,modelID2,inclCat);
        } else {
        // formula based
            // inclG was replaced with false in a previous commit. This caused problems in the optimizer and in
            // likelihood reporting (it was consistently worse than optimizer results)
            ((_Formula*)modelMatrixIndices.list_data[modelID])->ScanFForVariables(theReceptacle, inclG, false, inclCat);
        }
    }
}

//____________________________________________________________________________________

_String const _HYHBLTypeToText (long type) {
    _StringBuffer result (128L);
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
    
    result.Trim (0,result.length()-2L);
    return result;
}




//____________________________________________________________________________________


void _ElementaryCommand::ScanStringExpressionForHBLFunctions (_String* expression, _ExecutionList const& chain, bool recursive, _AVLListX& collection ) {
  
  _Formula f, f2;
  
  _String err_msg;
  _FormulaParsingContext fpc (&err_msg, chain.nameSpacePrefix);
  fpc.buildComplexObjects() = false;
  
  long     parseCode = Parse(&f,*expression,fpc,&f2);
  
  if (parseCode != HY_FORMULA_FAILED ) {
    f.ScanFormulaForHBLFunctions (collection, recursive);
    f2.ScanFormulaForHBLFunctions(collection, recursive);
  }

  
}

//____________________________________________________________________________________

void      _ElementaryCommand::BuildListOfDependancies    (_AVLListX & collection, bool recursive, _ExecutionList const & chain) {
  
  switch (code) {
      
    case 0:
    case 4:
    case 14:
    {
      if (parameters.lLength) {
        ScanStringExpressionForHBLFunctions((_String*)parameters (0), chain, recursive, collection);
      }
      break;
    }
      
    
    
  }
}


