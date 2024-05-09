/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon    (apoon42@uwo.ca)
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

#ifndef __BATCHLANGUAGE__
#define __BATCHLANGUAGE__


#include "global_things.h"
#include "parser.h"
#include "site.h"
#include "dataset.h"
#include "trie.h"
#include "associative_list.h"
#include "hy_string_buffer.h"

#include <stdio.h>


#define  HY_BL_ERROR_HANDLING_DEFAULT 0
#define  HY_BL_ERROR_HANDLING_SOFT    1

//____________________________________________________________________________________
struct    _CELInternals {
    _SimpleFormulaDatum     * values,
                            * stack;
    
    bool                    *is_compiled;

    _SimpleList             varList,
                            storeResults;

};

//____________________________________________________________________________________
struct    _HBLCommandExtras {
    long                cut_string;
    char                extract_condition_separator;
    _SimpleList         extract_conditions;
    _List               command_invocation;
    
    bool                do_trim,
                        is_assignment,
                        needs_verb;
};

class _ElementaryCommand;

//____________________________________________________________________________________
class   _ExecutionList: public _List // a sequence of commands to be executed
{
    
public:
    _ExecutionList (); // doesn't do much
    _ExecutionList (_StringBuffer&, _String* = nil, bool = false, bool* = nil);
    void Init (_String* = nil);

    virtual     ~_ExecutionList (void);

    virtual     BaseRef     makeDynamic (void) const;
    virtual     BaseRef     toStr (unsigned long = 0UL);
    virtual     void        Duplicate       (BaseRefConst);
    bool        BuildList                   (_StringBuffer&, _SimpleList* = nil, bool = false, bool = false);
    
    void        SetKWArgs   (_AssociativeList*);

    HBLObjectRef   Execute                     (_ExecutionList* parent = nil, bool ignore_parent_kwargs = false);
        // if parent is specified, copy stdin redirects from it
        // run this execution list
    HBLObjectRef   GetResult                   (void) {
        return result;
    }
    void        ExecuteSimple               (_ExecutionList * parent = nil);             // run a simple compiled list
    bool        TryToMakeSimple             (bool partial_ok = false);             // see if a list can be made into a compiled version

    void        ExecuteAndClean             (long);

    void        ResetFormulae               (void);             // decompile formulas (for reference functions)
    void        ResetNameSpace              (void);
    void        SetNameSpace                (_String const &);
    _String const
                GetFileName                 (void) const;
    _String*    GetNameSpace                (void);
    _String const    AddNameSpaceToID            (_String const&, _String const * = nil);
    _String     TrimNameSpaceFromID         (_String&);
  
    bool        has_stdin_redirect         (void) const {return stdinRedirect != nil;}
    bool        has_keyword_arguments      (void) const {return (kwargs && kwargs -> countitems()) || (kwarg_tags && kwarg_tags->countitems());}
  
    _String*    FetchFromStdinRedirect     (_String const * dialog_tag = nil, bool handle_multi_choice = false, bool do_echo = false);
    
    _ElementaryCommand* FetchLastCommand (void) {
        if (currentCommand - 1L < (long)lLength && currentCommand > 0L) {
            return (_ElementaryCommand*)(*this)(currentCommand - 1L);
        }
        return nil;
    }
    _ElementaryCommand* GetIthCommand       (long i) const {
        return (_ElementaryCommand*) GetItem (i);
    }

    void        GoToLastInstruction         (void) {
        currentCommand = MAX(currentCommand,(long)lLength-1L);
    }
    
    _StringBuffer const GenerateHelpMessage         (_AVLList * scanned_functions = nil) const;
    
    bool        IsErrorState    (void)     {
            return errorState;
    }

    void ClearExecutionList (void);


    void              ReportAnExecutionError (_String errMsg, bool doCommand = true, bool appendToExisting = false);
    /**
     * Handle an error message according to the reporting policy of this execution list (defined by errorHandlingMode)
     * @param errMsg -- the current command text stream
     * @param doCommand -- add standard text about the current command
     * @param appendToExisting -- append text to existing error
     
     */
  
    void        BuildListOfDependancies   (_AVLListX & collection, bool recursive = true, bool help_mode = false);
  
    /**
     
     Scan the body of this function/code for dependancies on various objects
     (currently only supports HBL functions), and store them in `collection`.
     
     If recursive is true, then new HBL functions will be scanned for dependancies
     as well.
     
    */

    /** Advance program counter */
    void      advance (void) {currentCommand ++;}
    
    bool      is_compiled (long idx = -1) const {if (cli) {if (idx < 0L) return true; else return cli->is_compiled[idx];} return false;}
    
    void      CopyCLIToVariables (void);
    void      StartProfile (void);
    _AssociativeList*   CollectProfile (void);
  
    // data fields
    // _____________________________________________________________

    long                            currentCommand,
                                    currentKwarg;
    
    char                            doProfile;
    int                             errorHandlingMode; // how does this execution list handle errors
    bool                            errorState;

    HBLObjectRef                    result;

    _VariableContainer*             nameSpacePrefix;
    _AssociativeList*               kwargs;

    _AVLListXL                      *stdinRedirect;
    _List                           *stdinRedirectAux,
                                    *kwarg_tags;
    
    /** SLKP 20190223
        kwarg_tags, if set, stores the ordered list of tagged inputs, which are created by
        invoking the
            KeywordArgument ("keyword", "description", "default");
        procedure. They are inherited down the call chain, like stdinRedirect and stdinRedirectAux
     */

    _String                         sourceFile,
                                    sourceText,
                                    enclosingNamespace;
    
    _SimpleList                     callPoints,
                                    lastif;

    _Matrix                         *profileCounter;

    _CELInternals                   *cli;
  
    protected:
        void       BuildExecuteCommandInstruction (_List * pieces, long code);
        void       BuildFscanf                    (_List * pieces, long code);
        void       BuildChoiceList                (_List * pieces, long code);


};

//____________________________________________________________________________________
// an elementary command



class   _ElementaryCommand: public _String // string contains the literal for this command
{
public:

    _ElementaryCommand (void); //dummy default constructor
    _ElementaryCommand (long); // with operation code
    _ElementaryCommand (_String& command); // process this string (and maybe an entire scope)
    // starting at a given position
    virtual                  ~_ElementaryCommand (void);

    virtual   BaseRef        makeDynamic (void) const;
    virtual   void           Duplicate (BaseRefConst);
    virtual   BaseRef        toStr (unsigned long = 0UL);

    bool      Execute        (_ExecutionList&); // perform this command in a given list
    void      ExecuteCase0   (_ExecutionList&);
    void      ExecuteCase4   (_ExecutionList&);
    void      ExecuteCase5   (_ExecutionList&);
    void      ExecuteDataFilterCases (_ExecutionList&);
    void      ExecuteCase11  (_ExecutionList&);
    void      ExecuteCase12  (_ExecutionList&);
    void      ExecuteCase31  (_ExecutionList&); // model construction
    void      ExecuteCase38  (_ExecutionList&, bool); // Reconstruct Ancestors
    void      ExecuteCase47  (_ExecutionList&); // ConstructStateCounter
    void      ExecuteCase52  (_ExecutionList&); // Simulate
    void      ExecuteCase53  (_ExecutionList&); // DoSQL
    void      ExecuteCase54  (_ExecutionList&); // Topology
    void      ExecuteCase58  (_ExecutionList&); // Profile Code
    void      ExecuteCase61  (_ExecutionList&); // SCFG
    void      ExecuteCase63  (_ExecutionList&); // NN; currently not functional
    void      ExecuteCase64  (_ExecutionList&); // BGM
    
    bool      HandleReplicateConstraint             (_ExecutionList&);
    bool      HandleAlignSequences                  (_ExecutionList&);
    bool      HandleConstructCategoryMatrix         (_ExecutionList&);
    bool      HandleGetDataInfo                     (_ExecutionList&);
    bool      HandleGetInformation                  (_ExecutionList&);
    bool      HandleFprintf                         (_ExecutionList&);
    bool      HandleHarvestFrequencies              (_ExecutionList&);
    bool      HandleOptimizeCovarianceMatrix        (_ExecutionList&, bool);
    bool      HandleComputeLFFunction               (_ExecutionList&);
    bool      HandleSelectTemplateModel             (_ExecutionList&);
    bool      HandleUseModel                        (_ExecutionList&);
    bool      HandleSetParameter                    (_ExecutionList&);
    bool      HandleAssert                          (_ExecutionList&);
    bool      HandleRequireVersion                  (_ExecutionList&);
    bool      HandleDeleteObject                    (_ExecutionList&);
    bool      HandleClearConstraints                (_ExecutionList&);
    bool      HandleMolecularClock                  (_ExecutionList&);
    bool      HandleGetURL                          (_ExecutionList&);
    bool      HandleGetString                       (_ExecutionList&);
    bool      HandleKeywordArgument                 (_ExecutionList&);
    bool      HandleExport                          (_ExecutionList&);
    bool      HandleDifferentiate                   (_ExecutionList&);
    bool      HandleFindRootOrIntegrate             (_ExecutionList&, bool do_integrate = false);
    bool      HandleConvertBranchLength             (_ExecutionList&);
    bool      HandleMPISend                         (_ExecutionList&);
    bool      HandleMPIReceive                      (_ExecutionList&);
    bool      HandleExecuteCommandsCases            (_ExecutionList&, bool do_load_from_file = false, bool do_load_library = false);
    bool      HandleDoSQL                           (_ExecutionList&);
    bool      HandleFscanf                          (_ExecutionList&, bool is_sscanf = false);
    bool      HandleChoiceList                      (_ExecutionList&);
    bool      HandleInitializeIterator              (_ExecutionList&);
    bool      HandleAdvanceIterator                 (_ExecutionList&);
    
    long      get_code                              (void) const { return code; };
    unsigned  long parameter_count                  (void) const { return parameters.countitems();}
    
    static    void      FindNextCommand       (_StringBuffer&, _StringBuffer&);
    // finds & stores the next command from _String into _StringBuffer
    // chops the input to remove the newly found line

    static  long      ExtractConditions     (_StringBuffer const& , long , _List&, char delimeter = ';', bool includeEmptyConditions = true);
    // used to extract the loop, if-then conditions

    static  bool      ExtractValidateAddHBLCommand (_StringBuffer& current_stream, const long command_code, _List* pieces, _HBLCommandExtras* command_spec, _ExecutionList& command_list);
    /**
     * Take a command from the current command stream, extract it, make an _ElementaryCommand and add it to the execution list
     * @param current_stream -- the current command text stream
     * @param command_code   -- the numerical code (from HY_HBL_COMMAND_*)
     * @param pieces         -- the list of parameters extracted from the () part of the command
     * @param command_spec   -- command specification structure
     * @param command_list   -- the command list object to append the command to
     * @return success/failure. 
     */
   

    static  bool      BuildFor              (_StringBuffer&, _ExecutionList&, _List*);
    // builds the for loop starting from
    // the beginning of input
    // this will process the loop header
    // and the entire scope afterwards

    static  bool      BuildIfThenElse       (_StringBuffer&, _ExecutionList&, _SimpleList*);
    // builds the if-then-else construct starting from
    // the beginning of input
    // this will process the loop header
    // and the entire scope afterwards

    static  bool      BuildWhile            (_StringBuffer&, _ExecutionList&, _List*);
    // builds the while(..) construct starting from
    // the beginning of input
    // this will process the loop header
    // and the entire scope afterwards

    static  bool      BuildDoWhile          (_StringBuffer&, _ExecutionList&);
    // builds the do {} while(..); construct starting from
    // the beginning of input
    // this will process the loop header
    // and the entire scope afterwards

    static  bool      ProcessInclude        (_StringBuffer&, _ExecutionList&);
    // processes the include command


    static  bool      ConstructDataSet      (_StringBuffer&, _ExecutionList&);
    // construct a dataset from the string


    static  bool      ConstructDataSetFilter(_StringBuffer&, _ExecutionList&);
    // construct a dataset filter from the string

    static  bool      ConstructTree         (_StringBuffer&, _ExecutionList&);
    // construct a tree


 
    static  bool      ConstructLF           (_StringBuffer&, _ExecutionList&);
    // construct a likelihood function


    static  bool      ConstructFunction     (_StringBuffer&, _ExecutionList&);
    // construct a fprintf command

    static  bool      ConstructReturn       (_StringBuffer&, _ExecutionList&);
    // construct a fprintf command

    static  bool      ConstructCategory     (_StringBuffer&, _ExecutionList&);
    // construct a category variable

    static  bool      ConstructChoiceList   (_StringBuffer&, _ExecutionList&);
    // construct a category variable

    static  bool      ConstructModel        (_StringBuffer&, _ExecutionList&);

   
    static  bool      ConstructProfileStatement
    (_StringBuffer&, _ExecutionList&);

 
    static  bool      ConstructSCFG         (_StringBuffer&, _ExecutionList&);

    static  bool      ConstructBGM          (_StringBuffer&, _ExecutionList&);

  
 
    static  bool      MakeGeneralizedLoop      (_String*, _String*, _String* , bool , _StringBuffer&, _ExecutionList&);
  
    bool              DecompileFormulae        (void);
  
    void              BuildListOfDependancies  (_AVLListX & collection, bool recursive, _ExecutionList const& chain, bool help_mode = false);
    
    
    
    /**
     
     Check this command for
     (currently only supports HBL functions), and store them in `collection`.
     
     If recursive is true, then new HBL functions will be scanned for dependancies
     as well.
     
     */
  
     static      const _List        fscanf_allowed_formats;


protected:
  
    static    void ScanStringExpressionForHBLFunctions (_String*, _ExecutionList const&, bool, _AVLListX& , bool help_mode = false);

    _String  *   GetIthParameter       (unsigned long i, bool range_check = true) const {
        BaseRef p = parameters.GetItemRangeCheck(i);
        if (!p && range_check) {
            hy_global::HandleApplicationError("Internal error in ElemenaryCommand::GetIthParameter", true);
            return new _String;
        }
        return (_String *)p;
    }


    bool      MakeJumpCommand       (_String*,  long, long, _ExecutionList&);
    // internal command used
    // to build a jump command
    // with two branches
    // and a condition
  
    void       addAndClean            (_ExecutionList&, _List* = nil, long = 0);
    void       appendCompiledFormulae (_Formula *, _Formula* = nil);


    friend  class     _ExecutionList;
    friend  void      DeleteVariable     (long, bool);
    friend  void      UpdateChangingFlas (long);
    friend  void      UpdateChangingFlas (_SimpleList const&);


private:
    _Variable*      _ValidateStorageVariable (_ExecutionList& program, unsigned long argument_index = 0UL) const;
    
    /**
        Extract the identifier from an expression like
        Type <id> = (...);
     
        @param ([in] _String) the source text to scan
        @param ([out] long) a positive integer to the first character following the '='; kNotFound if missing
        @param ([in] bool) if true, validate the identifier
        @param ([in] bool) if true, throw _String exceptions on errors
        @param ([in] long) start searching at this position in the string
     
        @return the identifier, or empty string if failed
        @version 0.1 SLKP 20190211
    */
    static const _String   ExtractStatementAssignment (_String const& source, long& end_at, const bool validate = true, const bool exceptions = true, const long offset = 0L);

    /**
     Process declaration of the form
     Type <id> = procedure (...);
     
     @param ([in] _String) the source text to scan
     @param ([out] procedure) the string for the procedure name (e.g. CreateFilter)
     @param ([out] pieces) comma separated arguments from the parentheses
     
     @return the identifier, or empty string if failed
     
     Throws _String exceptions
     
     @version 0.1 SLKP 20190212
     */
    static const _String   ProcessProcedureCall (_String const& source, _String& procedure, _List& pieces);

protected:  // data members

    _List       parameters;        // a list of parameters
    _SimpleList simpleParameters;  // a list of numeric parameters
    int         code;              // code describing this command
  

};

//____________________________________________________________________________________

_ElementaryCommand               * makeNewCommand       (long);

//____________________________________________________________________________________

#ifdef __HYPHYMPI__
#include <mpi.h>


#define  HYPHY_MPI_SIZE_TAG     111
#define  HYPHY_MPI_STRING_TAG   112
#define  HYPHY_MPI_DONE_TAG     113
#define  HYPHY_MPI_VARS_TAG     114
#define  HYPHY_MPI_DATA_TAG     115

#define  HYPHY_MPI_DIE_TAG      666


void     ReportMPIError         (int, bool);
void     MPISendString          (_String const&,long,bool=false);
_String* MPIRecvString          (long,long&);

#endif
//____________________________________________________________________________________

  //  EXTERN GLOBALS TODO SLKP 20180920 : these need to be reviewed and removed


extern  _List

dataSetList,
dataSetNamesList,
likeFuncList,
templateModelList,
scfgNamesList,
scfgList,
batchLanguageFunctions,

bgmNamesList,       // modified by afyp
bgmList,

likeFuncNamesList,
modelNames,
executionStack,
compiledFormulaeParameters,
standardLibraryPaths,
standardLibraryExtensions;


extern  _SimpleList
modelMatrixIndices,
modelTypeList,
// SLKP: 20100313 this list stores 0 for  normal (rate-matrix based models),
//       vs expression based matrices, for which the dimension is stored.
modelFrequenciesIndices,
listOfCompiledFormulae;


extern  _String

useLastFString,
mpiMLELFValue,
lf2SendBack,
defFileNameValue,
defFileString,
blConstructCM,
blFprintfRedirect               ,
blFprintfDevNull                ,
globalPolynomialCap             ,
enforceGlobalPolynomialCap      ,
dropPolynomialTerms             ,
maxPolyTermsPerVariable         ,
maxPolyExpIterates              ,
polyExpPrecision                ,
systemVariableDump              ,
selfDump                        ,
printDigitsSpec                 ,
explicitFormMExp                ,
multByFrequencies               ,
getDString                      ,
useLastFString                  ,
getFString                      ,
defFileString                   ,
VerbosityLevelString            ,
clearFile                       ,
keepFileOpen                    ,
closeFile                       ,
useLastDefinedMatrix            ,
MessageLogging                  ,
selectionStrings                ,
stdoutDestination               ,
messageLogDestination           ,
dataPanelSourcePath             ,
windowTypeTree                  ,
windowTypeClose                 ,
windowTypeTable                 ,
windowTypeDistribTable          ,
windowTypeDatabase              ,
screenWidthVar                  ,
screenHeightVar                 ,
useNexusFileData                ,
mpiMLELFValue                   ,
lf2SendBack                     ,
lfStartCompute                  ,
lfDoneCompute                   ,
getURLFileFlag                  ,
versionString                   ,
timeStamp                       ,
listLoadedLibraries             ,
simulationFilter                ,
prefixDS                        ,
prefixDF                        ,
prefixLF                        ,
replaceTreeStructure            ,
hyphyBaseDirectory              ,
platformDirectorySeparator      ,
covarianceParameterList         ,
matrixEvalCount                 ,
scfgCorpus                      ,
pathToCurrentBF                 ,
errorReportFormatExpression     ,
errorReportFormatExpressionStr  ,
errorReportFormatExpressionStack,
errorReportFormatExpressionStdin,
deferConstrainAssignment        ,
kBGMData                        ,
bgmConstraintMx                 ,
bgmParameters                   ,
bgmWeights                      ,
assertionBehavior               ,
dialogPrompt                    ,
_hyLastExecutionError           ,
_hyExecutionErrorMode           ,
blReturn                        ,
blDataSet                       ,
blDataSetFilter                 ,
blLF                            ,
blLF3                           ,
blTree                          ,
blTopology                      ,
blSCFG                          ;

extern  _ExecutionList              *currentExecutionList;

extern  _AVLList                    loadedLibraryPaths;
extern  _AVLListX                   _HY_HBLCommandHelper,
                                    _HY_GetStringGlobalTypes;


extern  _Trie                       _HY_ValidHBLExpressions,
                                    _HY_HBL_Namespaces,
                                    _HY_HBL_KeywordsPreserveSpaces;
 



const _String&  GetBFFunctionNameByIndex        (long);
long      GetBFFunctionArgumentCount  (long);
_List&    GetBFFunctionArgumentList   (long);
_SimpleList&    GetBFFunctionArgumentTypes   (long);
hyBLFunctionType
         GetBFFunctionType            (long);
_ExecutionList&
          GetBFFunctionBody           (long);

_String const
          ExportBFFunction            (long, bool = true);


void      ClearBFFunctionLists        (long = -1L);
bool      IsBFFunctionIndexValid      (long);
long      GetBFFunctionCount          (void);


void    ScanModelForVariables        (long modelID, _AVLList& theReceptacle, bool inclG, long modelID2, bool inclCat);
/* 20100316 SLKP:
    factored out a function call to scan a particular model
    for variables to permit the use of explicit (formula-based) model definitions
 */

void    ReadBatchFile                      (_String&, _ExecutionList&);
_String const ReturnDialogInput            (bool dispPath = false, _String const * rel_path = nil);
_String const ReturnFileDialogInput        (_String const* rel_path = nil);
_String const WriteFileDialogInput         (_String const* rel_path = nil);


hyFloat
_ProcessNumericArgumentWithExceptions (_String&,_VariableContainer const*);

hyFloat
ProcessNumericArgument                (_String*,_VariableContainer const*, _ExecutionList* = nil);
const _String ProcessLiteralArgument (_String const*,_VariableContainer const*, _ExecutionList* = nil);
_AssociativeList*
ProcessDictionaryArgument (_String* data, _VariableContainer* theP, _ExecutionList* = nil);

const _String GetStringFromFormula         (_String const*,_VariableContainer*);

void    SerializeModel               (_StringBuffer &,long,_AVLList* = nil, bool = false, _AssociativeList* options = nil);
bool    Get_a_URL                    (_String&,_String* = nil);

long    AddDataSetToList             (_String&,_DataSet*);
void    KillLFRecord                 (long, bool = true);
void    KillDataSetRecord            (long);
void    KillModelRecord              (long);
void    KillExplicitModelFormulae    (void);
bool    PushFilePath                 (_String&, bool = true, bool process = true);
_String const    PopFilePath         (void);
_String const *  PeekFilePath        (void);
_String const    GetPathStack        (const _String& spacer = ",");

void    RetrieveModelComponents      (long, _Matrix*&,     _Matrix*&, bool &);
void    RetrieveModelComponents      (long, _Variable*&, _Variable*&, bool &);
long    RetrieveModelFreq            (long);

void    ReadModelList                (void);
_String ProcessStringArgument        (_String* data);

const _String _hblCommandAccessor               (_ExecutionList*, long);
_String const _HYGenerateANameSpace             (void);
void          _HYClearANameSpace                (_String const&);

HBLObjectRef
ProcessAnArgumentByType      (_String const*, _VariableContainer const*, long, _ExecutionList* = nil);

void    _HBL_Init_Const_Arrays       (void);




/**
 An accessor function which attempts to retrieve a reference to a HyPhy Batch Language Object
 by name. A list of acceptable object classes can be specified in the type parameter. Note that
 types will be searched in the following order:
 
 HY_BL_DATASET,HY_BL_DATASET_FILTER,HY_BL_LIKELIHOOD_FUNCTION,HY_BL_SCFG,HY_BL_BGM,HY_BL_MODEL,HY_BL_HBL_FUNCTION
 
 i.e. if there is a dataset named 'foo' and a likelihood function named 'foo', then the dataset
 will be returned.
 
 
 
 @param   name provides a string with the name of the object to be retrieved.
 @param   type [in] which types of objects will be searched.
 [out] which type of object was retrieved (HY_BL_NOT_DEFINED if not found)
 @param   index (if not nil) will receive the index of the found object in the corresponding array
 @param   errMsg if set to True, will cause the function to report an error if no object of corresponding type could be found
 @param   tryLiteralLookup if set to True, will cause the function to, upon a failed lookup, to also try interpreting name as a string variable ID
 @return  pointer to the retrieved object or nil if not found
 @author  SLKP
 @version 20120324
 */


_String const _HYHBLTypeToText                (long type);

_HBLCommandExtras* _hyInitCommandExtras (const long = 0, const long = 0, const _String& = hy_global::kEmptyString, const char = ';', const bool = true, const bool = false, const bool = false, _SimpleList* = nil);


extern  bool                        numericalParameterSuccessFlag;
extern  hyFloat                  messageLogFlag;


extern enum       _hy_nested_check {
  _HY_NO_FUNCTION,
  _HY_FUNCTION,
  _HY_NAMESPACE } isInFunction;


#endif
