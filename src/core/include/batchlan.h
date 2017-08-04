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


#include "parser.h"
#include "site.h"
#include "trie.h"
#include "global_things.h"

#include <stdio.h>


#define  HY_BL_ERROR_HANDLING_DEFAULT 0
#define  HY_BL_ERROR_HANDLING_SOFT    1

//____________________________________________________________________________________
struct    _CELInternals {
    _SimpleFormulaDatum     * values,
                            * stack;

    _SimpleList       varList,
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
    _ExecutionList (_String&, _String* = nil, bool = false, bool* = nil);
    void Init (_String* = nil);

    virtual
    ~_ExecutionList (void);

    virtual
    BaseRef     makeDynamic (void) const;

    virtual
    BaseRef     toStr (unsigned long = 0UL);

    virtual
    void        Duplicate                   (BaseRefConst);
    bool        BuildList                   (_String&, _SimpleList* = nil, bool = false, bool = false);

    _PMathObj   Execute                     (_ExecutionList* parent = nil);
        // if parent is specified, copy stdin redirects from it
        // run this execution list
    _PMathObj   GetResult                   (void) {
        return result;
    }
    void        ExecuteSimple               (void);             // run a simple compiled list
    bool        TryToMakeSimple             (void);             // see if a list can be made into a compiled version

    long        ExecuteAndClean             (long,_String* = nil);

    void        ResetFormulae               (void);             // decompile formulas (for reference functions)
    void        ResetNameSpace              (void);
    void        SetNameSpace                (_String const &);
    _String const
                GetFileName                 (void) const;
    _String*    GetNameSpace                (void);
    _String const    AddNameSpaceToID            (_String const&, _String const * = nil);
    _String     TrimNameSpaceFromID         (_String&);
    _String*    FetchFromStdinRedirect      (void);
    _ElementaryCommand* FetchLastCommand (void) {
        if (currentCommand - 1 < lLength && currentCommand > 0) {
            return (_ElementaryCommand*)(*this)(currentCommand - 1);
        }
        return nil;
    }
    _ElementaryCommand* GetIthCommand       (long i) const {
        return (_ElementaryCommand*) GetItem (i);
    }

    void        GoToLastInstruction         (void) {
        currentCommand = MAX(currentCommand,lLength-1);
    }
    
    bool        IsErrorState    (void)     {
            return errorState;
    }

    void              ReportAnExecutionError (_String errMsg, bool doCommand = true, bool appendToExisting = false);
    /**
     * Handle an error message according to the reporting policy of this execution list (defined by errorHandlingMode)
     * @param errMsg -- the current command text stream
     * @param doCommand -- add standard text about the current command
     * @param appendToExisting -- append text to existing error
     
     */
  
    void        BuildListOfDependancies   (_AVLListX & collection, bool recursive = true);
  
    /**
     
     Scan the body of this function/code for dependancies on various objects
     (currently only supports HBL functions), and store them in `collection`.
     
     If recursive is true, then new HBL functions will be scanned for dependancies
     as well.
     
    */

    /** Advance program counter */
    void      advance (void) {currentCommand ++;}
  
    // data fields
    // _____________________________________________________________

    long                            currentCommand;
    char                            doProfile;
    int                             errorHandlingMode; // how does this execution list handle errors
    bool                            errorState;

    _PMathObj                       result;

    _VariableContainer*             nameSpacePrefix;

    _AVLListXL                      *stdinRedirect;

    _List                           *stdinRedirectAux;

    _String                         sourceFile,
                                    sourceText,
                                    enclosingNamespace;

    _SimpleList                     callPoints,
                                    lastif;

    _Matrix                         *profileCounter;

    _CELInternals                   *cli;

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
    void      ExecuteCase25  (_ExecutionList&, bool = false); // fscanf
    void      ExecuteCase31  (_ExecutionList&); // model construction
    void      ExecuteCase32  (_ExecutionList&); // list selection handler
    void      ExecuteCase34  (_ExecutionList&); // CovarianceMatrix
    void      ExecuteCase38  (_ExecutionList&, bool); // Reconstruct Ancestors
    void      ExecuteCase39  (_ExecutionList&); // Execute Commands
    void      ExecuteCase47  (_ExecutionList&); // ConstructStateCounter
    void      ExecuteCase52  (_ExecutionList&); // Simulate
    void      ExecuteCase53  (_ExecutionList&); // DoSQL
    void      ExecuteCase54  (_ExecutionList&); // Topology
    void      ExecuteCase57  (_ExecutionList&); // GetNeutralNull
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
    bool      HandleExport                          (_ExecutionList&);
    bool      HandleDifferentiate                   (_ExecutionList&);
    bool      HandleFindRootOrIntegrate             (_ExecutionList&, bool do_integrate = false);
    bool      HandleMPISend                         (_ExecutionList&);
    bool      HandleMPIReceive                      (_ExecutionList&);
  
    long      get_code                              (void) const { return code; };
    unsigned  long parameter_count                  (void) const { return parameters.countitems();}
    
    static  const _String   FindNextCommand       (_String&);
    // finds & returns the next command block in input
    // chops the input to remove the newly found line

    static  long      ExtractConditions     (_String& , long , _List&, char delimeter = ';', bool includeEmptyConditions = true);
    // used to extract the loop, if-then conditions

    static  bool      ExtractValidateAddHBLCommand (_String& current_stream, const long command_code, _List* pieces, _HBLCommandExtras* command_spec, _ExecutionList& command_list);
    /**
     * Take a command from the current command stream, extract it, make an _ElementaryCommand and add it to the execution list
     * @param current_stream -- the current command text stream
     * @param command_code   -- the numerical code (from HY_HBL_COMMAND_*)
     * @param pieces         -- the list of parameters extracted from the () part of the command
     * @param command_spec   -- command specification structure
     * @param command_list   -- the command list object to append the command to
     * @return success/failure. 
     */
   

    static  bool      BuildFor              (_String&, _ExecutionList&, _List*);
    // builds the for loop starting from
    // the beginning of input
    // this will process the loop header
    // and the entire scope afterwards

    static  bool      BuildIfThenElse       (_String&, _ExecutionList&, _SimpleList*);
    // builds the if-then-else construct starting from
    // the beginning of input
    // this will process the loop header
    // and the entire scope afterwards

    static  bool      BuildWhile            (_String&, _ExecutionList&, _List*);
    // builds the while(..) construct starting from
    // the beginning of input
    // this will process the loop header
    // and the entire scope afterwards

    static  bool      BuildDoWhile          (_String&, _ExecutionList&);
    // builds the do {} while(..); construct starting from
    // the beginning of input
    // this will process the loop header
    // and the entire scope afterwards

    static  bool      ProcessInclude        (_String&, _ExecutionList&);
    // processes the include command


    static  bool      ConstructDataSet      (_String&, _ExecutionList&);
    // construct a dataset from the string

    static  bool      ConstructExport       (_String&, _ExecutionList&);
    // construct a matrix export command

    static  bool      ConstructGetString    (_String&, _ExecutionList&);
    // construct a matrix import command

    static  bool      ConstructDataSetFilter(_String&, _ExecutionList&);
    // construct a dataset filter from the string

    static  bool      ConstructTree         (_String&, _ExecutionList&);
    // construct a tree

    static  bool      ConstructFscanf       (_String&, _ExecutionList&);
    // construct a fscanf command

    static  bool      ConstructExecuteCommands
    (_String&, _ExecutionList&);
    // construct a fscanf command

    static  bool      ConstructLF           (_String&, _ExecutionList&);
    // construct a likelihood function


    static  bool      ConstructFunction     (_String&, _ExecutionList&);
    // construct a fprintf command

    static  bool      ConstructReturn       (_String&, _ExecutionList&);
    // construct a fprintf command

    static  bool      ConstructSetParameter (_String&, _ExecutionList&);
    // construct a set parameter clause

    static  bool      ConstructCategory     (_String&, _ExecutionList&);
    // construct a category variable

    static  bool      ConstructChoiceList   (_String&, _ExecutionList&);
    // construct a category variable

    static  bool      ConstructModel        (_String&, _ExecutionList&);

    static  bool      ConstructMPIReceive   (_String&, _ExecutionList&);

    static  bool      ConstructStateCounter (_String&, _ExecutionList&);

    static  bool      ConstructDoSQL        (_String&, _ExecutionList&);

    static  bool      ConstructGetNeutralNull
    (_String&, _ExecutionList&);

    static  bool      ConstructProfileStatement
    (_String&, _ExecutionList&);

    static  bool      ConstructDeleteObject
    (_String&, _ExecutionList&);

    static  bool      ConstructSCFG         (_String&, _ExecutionList&);

    static  bool      ConstructBGM          (_String&, _ExecutionList&);

    static  bool      ConstructAssert       (_String&, _ExecutionList&);

    static  bool      SelectTemplateModel   (_String&, _ExecutionList&);

    static  bool      MakeGeneralizedLoop      (_String*, _String*, _String* , bool , _String&, _ExecutionList&);
  
    bool              DecompileFormulae        (void);
  
    void              BuildListOfDependancies  (_AVLListX & collection, bool recursive, _ExecutionList& chain);
    
    /**
     
     Check this command for
     (currently only supports HBL functions), and store them in `collection`.
     
     If recursive is true, then new HBL functions will be scanned for dependancies
     as well.
     
     */


protected:
  
    static    void ScanStringExpressionForHBLFunctions (_String*, _ExecutionList&, bool, _AVLListX& );

    _String  *   GetIthParameter       (unsigned long i, bool range_check = true) const {
        BaseRef p = parameters.GetItemRangeCheck(i);
        if (!p && range_check) {
            hy_global::HandleApplicationError("Internal error in ElemenaryCommand::GetIthParameter", true);
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
    friend  void      UpdateChangingFlas (_SimpleList&);
    
private:
    _Variable* _ValidateStorageVariable (_ExecutionList& program, unsigned long argument_index = 0UL) const;

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

extern  _List

dataSetList,
dataSetNamesList,
likeFuncList,
templateModelList,
scfgNamesList,
scfgList,

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
hasEndBeenReached               ,
clearFile                       ,
keepFileOpen                    ,
closeFile                       ,
useLastDefinedMatrix            ,
MessageLogging                  ,
selectionStrings                ,
useNoModel                      ,
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
lastModelUsed                   ,
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

extern  long                        matrixExpCount;
 

long      FindDataSetName                 (_String const&);
long      FindSCFGName                    (_String const&);
long      FindBFFunctionName              (_String const&, _VariableContainer const* = nil);
long    FindBgmName                       (_String const&);
// added by afyp, March 18, 2007

long    FindLikeFuncName                  (_String const&, bool = false);
long    FindModelName                     (_String const&);


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

void    ReadBatchFile                (_String&, _ExecutionList&);
_String ReturnDialogInput            (bool dispPath = false);
_String ReturnFileDialogInput        (void);
_String*ProcessCommandArgument       (_String*);
_String WriteFileDialogInput         (void);


hyFloat
_ProcessNumericArgumentWithExceptions (_String&,_VariableContainer const*);

hyFloat
ProcessNumericArgument                (_String*,_VariableContainer const*, _ExecutionList* = nil);
const _String ProcessLiteralArgument (_String const*,_VariableContainer const*, _ExecutionList* = nil);
_AssociativeList*
ProcessDictionaryArgument (_String* data, _VariableContainer* theP, _ExecutionList* = nil);

const _String GetStringFromFormula         (_String const*,_VariableContainer*);
void    ExecuteBLString              (_String&,_VariableContainer*);

void    SerializeModel               (_String&,long,_AVLList* = nil, bool = false);
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
bool    IsModelReversible            (long);
bool    IsModelOfExplicitForm        (long);

void    ReadModelList                (void);
_String ProcessStringArgument        (_String* data);

const _String _hblCommandAccessor          (_ExecutionList*, long);
_String _HYGenerateANameSpace             (void);

_PMathObj
ProcessAnArgumentByType      (_String const*, _VariableContainer const*, long, _ExecutionList* = nil);

void    _HBL_Init_Const_Arrays       (void);




BaseRefConst _HYRetrieveBLObjectByName              (_String const& name, long& type, long* index = nil, bool errMsg = false, bool tryLiteralLookup = false);
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

BaseRef _HYRetrieveBLObjectByNameMutable       (_String const& name, long& type, long* index = nil, bool errMsg = false, bool tryLiteralLookup = false);

_String const _HYHBLTypeToText                (long type);

_HBLCommandExtras* _hyInitCommandExtras (const long = 0, const long = 0, const _String& = hy_global::kEmptyString, const char = ';', const bool = true, const bool = false, const bool = false, _SimpleList* = nil);


extern  bool                        numericalParameterSuccessFlag;
extern  hyFloat                  messageLogFlag;

extern enum       _hy_nested_check {
  _HY_NO_FUNCTION,
  _HY_FUNCTION,
  _HY_NAMESPACE } isInFunction;


#endif
