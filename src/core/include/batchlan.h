/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2011
Primary Development:
  Sergei L Kosakovsky Pond (sergeilkp@mac.com)
Significant contributions from:
  Spencer V Muse (muse@stat.ncsu.edu)
  Simon DW Frost (sdfrost@ucsd.edu)
  Art FY Poon    (apoon@biomail.ucsd.edu)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#ifndef __BATCHLANGUAGE__
#define __BATCHLANGUAGE__


#include "parser.h"
#include "site.h"
#include "stdio.h"

#define   BL_FUNCTION_ALWAYS_UPDATE     0
#define   BL_FUNCTION_NORMAL_UPDATE     1


//!  Batch Language 'Object' type codes
/*!
     20110608 SLKP introduced.

     Bit flag style type tags and masks
     This is primarily to be used for object retrieval
     using the _HYRetrieveBLObjectByName function
*/

#define   HY_BL_NOT_DEFINED             0
#define   HY_BL_DATASET                 1
#define   HY_BL_DATASET_FILTER          2
#define   HY_BL_LIKELIHOOD_FUNCTION     4
#define   HY_BL_SCFG                    8
#define   HY_BL_BGM                     16
#define   HY_BL_MODEL                   32
#define   HY_BL_HBL_FUNCTION            64

#define   HY_BL_ANY                     1023


//____________________________________________________________________________________
struct    _CELInternals {
    _SimpleFormulaDatum     * values,
                            * stack;

    _SimpleList       varList,
                      storeResults;

};

//____________________________________________________________________________________
class   _ExecutionList: public _List // a sequence of commands to be executed
{
public:
    _ExecutionList (); // doesn't do much
    _ExecutionList (_String&, _String* = nil, bool = false);

    virtual
    ~_ExecutionList (void);

    virtual
    BaseRef     makeDynamic (void);

    virtual
    BaseRef     toStr (void);

    virtual
    void        Duplicate                   (BaseRef);
    bool        BuildList                   (_String&, _SimpleList* = nil, bool = false);

    _PMathObj   Execute                     (void);             // run this execution list
    _PMathObj   GetResult                   (void) {
        return result;
    }
    void        ExecuteSimple               (void);             // run a simple compiled list
    bool        TryToMakeSimple             (void);             // see if a list can be made into a compiled version

    long        ExecuteAndClean             (long,_String* = nil);

    void        ResetFormulae               (void);             // decompile formulas (for reference functions)
    void        ResetNameSpace              (void);
    void        SetNameSpace                (_String);
    _String*    GetNameSpace                (void);
    _String     AddNameSpaceToID            (_String&);
    _String     TrimNameSpaceFromID         (_String&);
    _String*    FetchFromStdinRedirect      (void);

    void        GoToLastInstruction         (void) {
        currentCommand = MAX(currentCommand,lLength-1);
    }



    // data fields
    // _____________________________________________________________

    long                            currentCommand;
    char                            doProfile;

    _PMathObj                       result;

    _VariableContainer*             nameSpacePrefix;

    _AVLListXL                      *stdinRedirect;

    _List                           *stdinRedirectAux;

    _String                         sourceFile,
                                    sourceText;

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

    virtual   BaseRef        makeDynamic (void);
    virtual   void           Duplicate (BaseRef);
    virtual   BaseRef        toStr (void);

    bool      Execute        (_ExecutionList&); // perform this command in a given list
    void      ExecuteCase0   (_ExecutionList&);
    void      ExecuteCase4   (_ExecutionList&);
    void      ExecuteCase5   (_ExecutionList&);
    void      ExecuteDataFilterCases
    (_ExecutionList&);
    void      ExecuteCase8   (_ExecutionList&);
    void      ExecuteCase11  (_ExecutionList&);
    void      ExecuteCase12  (_ExecutionList&);
    void      ExecuteCase17  (_ExecutionList&);
    void      ExecuteCase21  (_ExecutionList&);
    void      ExecuteCase24  (_ExecutionList&);
    void      ExecuteCase25  (_ExecutionList&, bool = false); // fscanf
    void      ExecuteCase26  (_ExecutionList&); // ReplicateConstraint
    void      ExecuteCase31  (_ExecutionList&); // model construction
    void      ExecuteCase32  (_ExecutionList&); // list selection handler
    void      ExecuteCase33  (_ExecutionList&); // index string selector
    void      ExecuteCase34  (_ExecutionList&); // CovarianceMatrix
    void      ExecuteCase35  (_ExecutionList&); // SetParameter
    void      ExecuteCase36  (_ExecutionList&); // OpenDataPanel
    void      ExecuteCase37  (_ExecutionList&); // GetInformation
    void      ExecuteCase38  (_ExecutionList&, bool); // Reconstruct Ancestors
    void      ExecuteCase39  (_ExecutionList&); // Execute Commands
    void      ExecuteCase40  (_ExecutionList&); // Open Window
    void      ExecuteCase41  (_ExecutionList&); // Spawn LF
    void      ExecuteCase42  (_ExecutionList&); // Differentiate
    void      ExecuteCase43  (_ExecutionList&); // FindRoot
    void      ExecuteCase44  (_ExecutionList&); // MPISend
    void      ExecuteCase45  (_ExecutionList&); // MPIReceive
    void      ExecuteCase46  (_ExecutionList&); // GetDataInfo
    void      ExecuteCase47  (_ExecutionList&); // ConstructStateCounter
    void      ExecuteCase49  (_ExecutionList&); // LFCompute
    void      ExecuteCase51  (_ExecutionList&); // GetURL
    void      ExecuteCase52  (_ExecutionList&); // Simulate
    void      ExecuteCase53  (_ExecutionList&); // DoSQL
    void      ExecuteCase54  (_ExecutionList&); // Topology
    void      ExecuteCase55  (_ExecutionList&); // AlignSequences
    void      ExecuteCase57  (_ExecutionList&); // GetNeutralNull
    void      ExecuteCase58  (_ExecutionList&); // Profile Code
    void      ExecuteCase59  (_ExecutionList&); // DeleteObject
    void      ExecuteCase60  (_ExecutionList&); // RequireVersion
    void      ExecuteCase61  (_ExecutionList&); // SCFG
    void      ExecuteCase63  (_ExecutionList&); // NN; currently not functional
    void      ExecuteCase64  (_ExecutionList&); // BGM
    void      ExecuteCase65  (_ExecutionList&); // assert

    static  _String   FindNextCommand       (_String&, bool = false);
    // finds & returns the next command block in input
    // chops the input to remove the newly found line

    static  long      ExtractConditions     (_String& , long , _List&, char delimeter = ';', bool includeEmptyConditions = true);
    // used to extract the loop, if-then conditions

    static  bool      BuildFor              (_String&, _ExecutionList&);
    // builds the for loop starting from
    // the beginning of input
    // this will process the loop header
    // and the entire scope afterwards

    static  bool      BuildIfThenElse       (_String&, _ExecutionList&, _SimpleList*);
    // builds the if-then-else construct starting from
    // the beginning of input
    // this will process the loop header
    // and the entire scope afterwards

    static  bool      BuildWhile            (_String&, _ExecutionList&);
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

    static  bool      ConstructImport       (_String&, _ExecutionList&);
    // construct a matrix import command

    static  bool      ConstructGetString    (_String&, _ExecutionList&);
    // construct a matrix import command

    static  bool      ConstructDataSetFilter(_String&, _ExecutionList&);
    // construct a dataset filter from the string

    static  bool      ConstructTree         (_String&, _ExecutionList&);
    // construct a tree

    static  bool      ConstructFprintf      (_String&, _ExecutionList&);
    // construct a fprintf command

    static  bool      ConstructFscanf       (_String&, _ExecutionList&);
    // construct a fscanf command

    static  bool      ConstructExecuteCommands
    (_String&, _ExecutionList&);
    // construct a fscanf command

    static  bool      ConstructReplicateConstraint
    (_String&, _ExecutionList&);
    // construct a replicate constraint command

    static  bool      ConstructOptimize     (_String&, _ExecutionList&);


    static  bool      ConstructLF           (_String&, _ExecutionList&);
    // construct a likelihood function


    static  bool      ConstructHarvestFreq  (_String&, _ExecutionList&);
    // construct a fprintf command

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

    static  bool      ConstructClearConstraints (_String&, _ExecutionList&);
    // construct a clear constraints command

    static  bool      ConstructMolecularClock (_String&, _ExecutionList&);
    // construct a molecular clock constraint

    static  bool      ConstructCategoryMatrix (_String&, _ExecutionList&);
    // construct a category matrix for the optimized like func

    static  bool      ConstructOpenDataPanel (_String&, _ExecutionList&);
    // open data panel with given settings

    static  bool      ConstructUseMatrix    (_String&, _ExecutionList&);

    static  bool      ConstructOpenWindow   (_String&, _ExecutionList&);

    static  bool      ConstructSpawnLF      (_String&, _ExecutionList&);

    static  bool      ConstructDifferentiate(_String&, _ExecutionList&);

    static  bool      ConstructFindRoot     (_String&, _ExecutionList&);

    static  bool      ConstructGetInformation
    (_String&, _ExecutionList&);

    static  bool      ConstructModel        (_String&, _ExecutionList&);

    static  bool      ConstructMPISend      (_String&, _ExecutionList&);

    static  bool      ConstructMPIReceive   (_String&, _ExecutionList&);

    static  bool      ConstructGetDataInfo  (_String&, _ExecutionList&);

    static  bool      ConstructStateCounter (_String&, _ExecutionList&);

    static  bool      SetDialogPrompt       (_String&, _ExecutionList&);

    static  bool      ConstructGetURL       (_String&, _ExecutionList&);

    static  bool      ConstructDoSQL        (_String&, _ExecutionList&);

    static  bool      ConstructAlignSequences
    (_String&, _ExecutionList&);

    static  bool      ConstructGetNeutralNull
    (_String&, _ExecutionList&);

    static  bool      ConstructProfileStatement
    (_String&, _ExecutionList&);

    static  bool      ConstructDeleteObject
    (_String&, _ExecutionList&);

    static  bool      ConstructRequireVersion
    (_String&, _ExecutionList&);

    static  bool      ConstructSCFG         (_String&, _ExecutionList&);

    static  bool      ConstructNN           (_String&, _ExecutionList&);

    static  bool      ConstructBGM          (_String&, _ExecutionList&);

    static  bool      ConstructAssert       (_String&, _ExecutionList&);

    static  bool      SelectTemplateModel   (_String&, _ExecutionList&);

    static  bool      MakeGeneralizedLoop   (_String*, _String*, _String* , bool , _String&, _ExecutionList&);

protected:

    bool      MakeJumpCommand       (_String*,  long, long, _ExecutionList&);
    // internal command used
    // to build a jump command
    // with two branches
    // and a condition

    void       addAndClean          (_ExecutionList&, _List* = nil, long = 0);


    friend  class     _ExecutionList;
    friend  void      DeleteVariable     (long, bool);
    friend  void      UpdateChangingFlas (long);
    friend  void      UpdateChangingFlas (_SimpleList&);

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

extern   _String                mpiNodeID,
         mpiNodeCount;

#define  HYPHY_MPI_SIZE_TAG     111
#define  HYPHY_MPI_STRING_TAG   112
#define  HYPHY_MPI_DONE_TAG     113
#define  HYPHY_MPI_VARS_TAG     114
#define  HYPHY_MPI_DATA_TAG     115

#define  HYPHY_MPI_DIE_TAG      666


void     ReportMPIError         (int, bool);
void     MPISendString          (_String&,long,bool=false);
_String* MPIRecvString          (long,long&);

#endif
//____________________________________________________________________________________

extern  _List

batchLanguageFunctions,
batchLanguageFunctionNames,
batchLanguageFunctionParameterLists,
dataSetList,
dataSetNamesList,
likeFuncList,
dataSetFilterList,
dataSetFilterNamesList,
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

batchLanguageFunctionParameters,
batchLanguageFunctionClassification,
modelMatrixIndices,
modelTypeList,
// SLKP: 20100313 this list stores 0 for  normal (rate-matrix based models),
//       vs expression based matrices, for which the dimension is stored.
modelFrequenciesIndices,
listOfCompiledFormulae;


extern  _String

getDString,
getFString,
baseDirectory,
libDirectory,
useLastFString,
mpiMLELFValue,
lf2SendBack,
hyphyBaseDirectory,
platformDirectorySeparator,
defFileNameValue,
defFileString,
blConstructCM,
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
useLastModel                    ,
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
lastModelParameterList          ,
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
pcAmbiguitiesResolve            ,
pcAmbiguitiesAverage            ,
pcAmbiguitiesSkip               ,
lfStartCompute                  ,
lfDoneCompute                   ,
getURLFileFlag                  ,
versionString                   ,
timeStamp                       ,
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
bgmData                         ,
bgmWeights                      ,
pathToCurrentBF                 ,
statusBarUpdateString           ,
statusBarProgressValue          ,
errorReportFormatExpression     ,
errorReportFormatExpressionStr  ,
errorReportFormatExpressionStack,
errorReportFormatExpressionStdin,

#ifdef      __HYPHYMPI__
mpiNodeID                       ,
mpiNodeCount                    ,
mpiLastSentMsg                  ,
#endif
hfCountGap                      ;

extern  _ExecutionList              *currentExecutionList;

extern  _AVLList
loadedLibraryPaths;


long    FindDataSetName              (_String&);
long    FindDataSetFilterName        (_String&);
long    FindSCFGName                 (_String&);
long    FindBFFunctionName           (_String&, _VariableContainer* = nil);

long    FindBgmName                  (_String &);
// added by afyp, March 18, 2007

long    FindLikeFuncName             (_String&, bool = false);
long    FindModelName                (_String&);
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
_Parameter
ProcessNumericArgument       (_String*,_VariableContainer*);
_String ProcessLiteralArgument       (_String*,_VariableContainer*);
_AssociativeList*
ProcessDictionaryArgument (_String* data, _VariableContainer* theP);

_String GetStringFromFormula         (_String*,_VariableContainer*);
void    ExecuteBLString              (_String&,_VariableContainer*);

void    SerializeModel               (_String&,long,_AVLList* = nil, bool = false);
bool    Get_a_URL                    (_String&,_String* = nil);

long    AddFilterToList              (_String&,_DataSetFilter*,bool = false);
long    AddDataSetToList             (_String&,_DataSet*);
void    SetDataFilterParameters      (_String&, _DataSetFilter*, bool);
void    KillDataFilterRecord         (long, bool = false);
void    KillLFRecord                 (long, bool = true);
void    KillDataSetRecord            (long);
void    KillModelRecord              (long);
void    KillExplicitModelFormulae    (void);
bool    PushFilePath                 (_String&);
void    PopFilePath                  (void);
_Matrix*CheckMatrixArg               (_String*, bool);
_AssociativeList *
CheckAssociativeListArg      (_String*);
void    RetrieveModelComponents      (long, _Matrix*&,     _Matrix*&, bool &);
void    RetrieveModelComponents      (long, _Variable*&, _Variable*&, bool &);
bool    IsModelReversible            (long);

_PMathObj
ProcessAnArgumentByType      (_String*, _VariableContainer*, long);

void    _HBL_Init_Const_Arrays       (void);
void    ReturnCurrentCallStack       (_List&, _List&);

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
    @return  pointer to the retrieved object or nil if not found
    @author  SLKP
    @version 20110608
*/

BaseRef _HYRetrieveBLObjectByName    (_String& name, long& type, long* index = nil);


extern  bool                        numericalParameterSuccessFlag;
extern  _Parameter                  messageLogFlag;


#endif
